MODULE p4zfechem
   !!======================================================================
   !!                         ***  MODULE p4zfechem  ***
   !! TOP :   PISCES Compute iron chemistry and scavenging
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, A. Tagliabue, C. Ethe) Original code
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p4z_fechem       : Compute remineralization/scavenging of iron
   !!   p4z_fechem_init  : Initialisation of parameters for remineralisation
   !!   p4z_fechem_alloc : Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zche          ! chemical model
   USE p4zsbc          ! Boundary conditions from sediments
   USE prtctl_trc      ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_fechem        ! called in p4zbio.F90
   PUBLIC   p4z_fechem_init   ! called in trcsms_pisces.F90

   LOGICAL          ::   ln_ligvar    !: boolean for variable ligand concentration following Tagliabue and voelker
   REAL(wp), PUBLIC ::   xlam1        !: scavenging rate of Iron 
   REAL(wp), PUBLIC ::   xlamdust     !: scavenging rate of Iron by dust 
   REAL(wp), PUBLIC ::   ligand       !: ligand concentration in the ocean 
   REAL(wp), PUBLIC ::   kfep         !: rate constant for nanoparticle formation
   LOGICAL          ::  ln_culigvar    !: CAM boolean for variable cu ligandconcentration following Tagliabue and voelker
   REAL(wp), PUBLIC :: culig
   REAL(wp), PUBLIC :: cukdp
   REAL(wp), PUBLIC :: cukdg
   REAL(wp), PUBLIC ::  znlig        !:
   REAL(wp), PUBLIC ::  znkdp         !: Kd for Zn scav (1/(mmol particle /m3)
   REAL(wp), PUBLIC ::  znkdg         !: Kd for Zn scav (1/(mmol particle /m3)
!!!!!! Co model
   REAL(wp), PUBLIC :: coscav
   REAL(wp), PUBLIC :: colig
   REAL(wp), PUBLIC :: coscav_ox
   REAL(wp), PUBLIC :: codiss_ox
   REAL(wp), PUBLIC :: kbacts
   REAL(wp), PUBLIC :: kpars
   REAL(wp), PUBLIC :: fnan
   REAL(wp), PUBLIC :: ko2s
!!!!!! Mn model
   REAL(wp), PUBLIC :: mnscav
   REAL(wp), PUBLIC :: mnlig
   REAL(wp), PUBLIC :: mnscav_ox
   REAL(wp), PUBLIC :: mndiss_ox
   REAL(wp), PUBLIC :: kbacts_mn
   REAL(wp), PUBLIC :: kpars_mn
   REAL(wp), PUBLIC :: ko2s_mn
   REAL(wp), PUBLIC :: rdiss

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zfechem.F90 10416 2018-12-19 11:45:43Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_fechem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_fechem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of iron
      !!
      !! ** Method  :   A simple chemistry model of iron from Aumont and Bopp (2006)
      !!                based on one ligand and one inorganic form
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk, jic, jn
      REAL(wp) ::   zdep, zlam1a, zlam1b, zlamfac
      REAL(wp) ::   zkeq, zfeequi, zfesatur, zfecoll, fe3sol
      REAL(wp) ::   zdenom1, zscave, zaggdfea, zaggdfeb, zcoag
      REAL(wp) ::   ztrc, zdust
      REAL(wp) ::   zdenom2
      REAL(wp) ::   zzFeL1, zzFeL2, zzFe2, zzFeP, zzFe3, zzstrn2
      REAL(wp) ::   zrum, zcodel, zargu, zlight
      REAL(wp) ::   zkox, zkph1, zkph2, zph, zionic, ztligand
      REAL(wp) ::   za, zb, zc, zkappa1, zkappa2, za0, za1, za2
      REAL(wp) ::   zxs, zfunc, zp, zq, zd, zr, zphi, zfff, zp3, zq2
      REAL(wp) ::   ztfe, zoxy, zhplus, zxlam
      REAL(wp) ::   zaggliga, zaggligb
      REAL(wp) ::   dissol, zligco
      REAL(wp) :: zrfact2
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zTL1, zFe3, ztotlig, precip, zFeL1
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zcoll3d, zscav3d, zlcoll3d
      REAL(wp) :: tcu, ztcu, fcu, zCuscav, zcukeq, zcufesatur, zscavecu  !Cam copper model
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zCufree, zcul, cuflux, zcutotlig ! Cam Cu model
      REAL(wp) :: ztrc1, ztrc2
      REAL(wp) :: tzn, ztzn, fzn, zZnscav, zzkeq, zzfesatur, znkd1
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zZnfree, zznl, znflux ! Cam Zn model
      REAL(wp) :: zlcob, zcoscav, znanof, kpar, ko2, kbact, lam, minlam, o2a
      REAL(wp) :: zcodesor, tgfun, zcopcp
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zlcot, cscav, cdiss
      INTEGER  :: ikt, ikt_mn
!!! Cam Mn model: I took all the CO variables and copied them to Mn
      REAL(wp) :: zlmnb, zmnscav, znanof_mn, kpar_mn, ko2_mn, kbact_mn, lam_mn, o2a_mn
      REAL(wp) :: zmndesor, tgfun_mn, zmnpcp
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zlmnt, mscav, mdiss
      REAL(wp) :: rdiss, rdiss2
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_fechem')
      !
      zFe3 (:,:,:) = 0.
      zFeL1(:,:,:) = 0.
      zTL1 (:,:,:) = 0.
IF (ln_copper) THEN !Cam
zCufree(:,:,:) = 0.
zcul(:,:,:) = 0.
cuflux(:,:,:) = 0.
zcutotlig(:,:,:) = 0.
ENDIF

IF (ln_cobalt) THEN 
zlcot(:,:,:) =0.
cscav(:,:,:) =0.
cdiss(:,:,:) =0.
zcofree(:,:,:) = 0.
ENDIF

IF (ln_manganese) THEN
zlmnt(:,:,:) =0.
mscav(:,:,:) =0.
mdiss(:,:,:) =0.
zmnfree(:,:,:) = 0.
ENDIF

IF (ln_zinc) THEN 
zZnfree(:,:,:) = 0.
zznl(:,:,:) =0.
znflux(:,:,:) =0.
ENDIF
      ! Total ligand concentration : Ligands can be chosen to be constant or variable
      ! Parameterization from Tagliabue and Voelker (2011)
      ! -------------------------------------------------
      IF( ln_ligvar ) THEN
         ztotlig(:,:,:) =  0.09 * trb(:,:,:,jpdoc) * 1E6 + ligand * 1E9
         ztotlig(:,:,:) =  MIN( ztotlig(:,:,:), 10. )
      ELSE
        IF( ln_ligand ) THEN  ;   ztotlig(:,:,:) = trb(:,:,:,jplgw) * 1E9
        ELSE                  ;   ztotlig(:,:,:) = ligand * 1E9
        ENDIF
      ENDIF

      ! ------------------------------------------------------------
      !  from Aumont and Bopp (2006)
      ! This model is based on one ligand and Fe' 
      ! Chemistry is supposed to be fast enough to be at equilibrium
      ! ------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zTL1(ji,jj,jk)  = ztotlig(ji,jj,jk)
               zkeq            = fekeq(ji,jj,jk)
               zfesatur        = zTL1(ji,jj,jk) * 1E-9
               ztfe            = trb(ji,jj,jk,jpfer) 
               ! Fe' is the root of a 2nd order polynom
               zFe3 (ji,jj,jk) = ( -( 1. + zfesatur * zkeq - zkeq * ztfe )               &
                  &              + SQRT( ( 1. + zfesatur * zkeq - zkeq * ztfe )**2       &
                  &              + 4. * ztfe * zkeq) ) / ( 2. * zkeq )
               zFe3 (ji,jj,jk) = zFe3(ji,jj,jk) * 1E9
               zFeL1(ji,jj,jk) = MAX( 0., trb(ji,jj,jk,jpfer) * 1E9 - zFe3(ji,jj,jk) )
           END DO
         END DO
      END DO
         !

      zdust = 0.         ! if no dust available
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Scavenging rate of iron. This scavenging rate depends on the load of particles of sea water. 
               ! This parameterization assumes a simple second order kinetics (k[Particles][Fe]).
               ! Scavenging onto dust is also included as evidenced from the DUNE experiments.
               ! --------------------------------------------------------------------------------------
               zhplus  = max( rtrn, hi(ji,jj,jk) )
               fe3sol  = fesol(ji,jj,jk,1) * ( zhplus**3 + fesol(ji,jj,jk,2) * zhplus**2  &
               &         + fesol(ji,jj,jk,3) * zhplus + fesol(ji,jj,jk,4)     &
               &         + fesol(ji,jj,jk,5) / zhplus )
               !
               zfeequi = zFe3(ji,jj,jk) * 1E-9
               zhplus  = max( rtrn, hi(ji,jj,jk) )
               fe3sol  = fesol(ji,jj,jk,1) * ( zhplus**3 + fesol(ji,jj,jk,2) * zhplus**2  &
                  &         + fesol(ji,jj,jk,3) * zhplus + fesol(ji,jj,jk,4)     &
                  &         + fesol(ji,jj,jk,5) / zhplus )
               zfecoll = 0.5 * zFeL1(ji,jj,jk) * 1E-9
               ! precipitation of Fe3+, creation of nanoparticles
               precip(ji,jj,jk) = MAX( 0., ( zFe3(ji,jj,jk) * 1E-9 - fe3sol ) ) * kfep * xstep
               !
               ztrc   = ( trb(ji,jj,jk,jppoc) + trb(ji,jj,jk,jpgoc) + trb(ji,jj,jk,jpcal) + trb(ji,jj,jk,jpgsi) ) * 1.e6 
               IF( ln_dust )  zdust  = dust(ji,jj) / ( wdust / rday ) * tmask(ji,jj,jk) &
               &  * EXP( -gdept_n(ji,jj,jk) / 540. )
               IF (ln_ligand) THEN
                  zxlam  = xlam1 * MAX( 1.E-3, EXP(-2 * etot(ji,jj,jk) / 10. ) * (1. - EXP(-2 * trb(ji,jj,jk,jpoxy) / 100.E-6 ) ))
               ELSE
                  zxlam  = xlam1 * 1.0
               ENDIF
               zlam1b = 3.e-5 + xlamdust * zdust + zxlam * ztrc
               zscave = zfeequi * zlam1b * xstep

               ! Compute the different ratios for scavenging of iron
               ! to later allocate scavenged iron to the different organic pools
               ! ---------------------------------------------------------
               zdenom1 = zxlam * trb(ji,jj,jk,jppoc) / zlam1b
               zdenom2 = zxlam * trb(ji,jj,jk,jpgoc) / zlam1b

               !  Increased scavenging for very high iron concentrations found near the coasts 
               !  due to increased lithogenic particles and let say it is unknown processes (precipitation, ...)
               !  -----------------------------------------------------------
               zlamfac = MAX( 0.e0, ( gphit(ji,jj) + 55.) / 30. )
               zlamfac = MIN( 1.  , zlamfac )
               zdep    = MIN( 1., 1000. / gdept_n(ji,jj,jk) )
               zcoag   = 1E-4 * ( 1. - zlamfac ) * zdep * xstep * trb(ji,jj,jk,jpfer)

               !  Compute the coagulation of colloidal iron. This parameterization 
               !  could be thought as an equivalent of colloidal pumping.
               !  It requires certainly some more work as it is very poorly constrained.
               !  ----------------------------------------------------------------
               zlam1a   = ( 0.369  * 0.3 * trb(ji,jj,jk,jpdoc) + 102.4  * trb(ji,jj,jk,jppoc) ) * xdiss(ji,jj,jk)    &
                   &      + ( 114.   * 0.3 * trb(ji,jj,jk,jpdoc) )
               zaggdfea = zlam1a * xstep * zfecoll
               !
               zlam1b   = 3.53E3 * trb(ji,jj,jk,jpgoc) * xdiss(ji,jj,jk)
               zaggdfeb = zlam1b * xstep * zfecoll
               !
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zscave - zaggdfea - zaggdfeb &
               &                     - zcoag - precip(ji,jj,jk)
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * zdenom1 + zaggdfea
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zscave * zdenom2 + zaggdfeb
               zscav3d(ji,jj,jk)   = zscave
               zcoll3d(ji,jj,jk)   = zaggdfea + zaggdfeb
               !
            END DO
         END DO
      END DO
      !
      !  Define the bioavailable fraction of iron
      !  ----------------------------------------
      biron(:,:,:) = trb(:,:,:,jpfer) 
      !
      IF( ln_ligand ) THEN
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zlam1a   = ( 0.369  * 0.3 * trb(ji,jj,jk,jpdoc) + 102.4  * trb(ji,jj,jk,jppoc) ) * xdiss(ji,jj,jk)    &
                      &    + ( 114.   * 0.3 * trb(ji,jj,jk,jpdoc) )
                  !
                  zlam1b   = 3.53E3 *   trb(ji,jj,jk,jpgoc) * xdiss(ji,jj,jk)
                  zligco   = 0.5 * trn(ji,jj,jk,jplgw)
                  zaggliga = zlam1a * xstep * zligco
                  zaggligb = zlam1b * xstep * zligco
                  tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) - zaggliga - zaggligb
                  zlcoll3d(ji,jj,jk)  = zaggliga + zaggligb
               END DO
            END DO
         END DO
         !
         plig(:,:,:) =  MAX( 0., ( ( zFeL1(:,:,:) * 1E-9 ) / ( trb(:,:,:,jpfer) +rtrn ) ) )
         !
      ENDIF

      IF( ln_cobalt ) THEN
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                znanof = min( phyrel(ji,jj,jk) , fnan )
                kpar = MIN( 0.9 ,  etot(ji,jj,jk)**2 / (  etot(ji,jj,jk)**2 + kpars**2 ) )
                zlcot(ji,jj,jk) = min( 150.E-9 , trb(ji,jj,jk,jpdco) * znanof ) ! set maximum
                zlcot(ji,jj,jk) = max( (colig*(1-kpar))*1000. , zlcot(ji,jj,jk) ) ! set background Co ligand concentration
                zcofree(ji,jj,jk) = max( 0., ( trn(ji,jj,jk,jpdco) - zlcot(ji,jj,jk) ) )
            END DO
         END DO
      END DO
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                kpar = MIN( 1. ,  etot(ji,jj,jk)**2 / (  etot(ji,jj,jk)**2 + kpars**2))
                o2a = MAX( 0. , trb(ji,jj,jk,jpoxy) - coscav_ox )
                o2a = MAX( 10.E-6 , o2a )
                ko2 = MAX(0.25, o2a**2 / (o2a**2 + ko2s**2))
                kbact = bacbio(ji,jj,jk)**2 / ( bacbio(ji,jj,jk)**2 + kbacts**2)
                tgfun = 0.5 * EXP( 0.1012 * tsn(ji,jj,jk,jp_tem) ) ! Q10 = 2.75 Lee and Fisher
                lam = MAX( 0., (1 - kpar) * ko2 * kbact )
                minlam = 0.025 * 10 * ( trb(ji,jj,jk,jpdco)**4 / (trb(ji,jj,jk,jpdco)**4 + 200.E-9**4 ) )
                minlam = MAX( MIN(minlam, 0.025) , 0.001 )
                IF( gdept_n(ji,jj,jk) > hmld(ji,jj) ) lam = max( lam , minlam )
                zlam1b = coscav * tgfun * lam ! cobalt scavenging rate
                zcoscav = zlam1b * xstep * zcofree(ji,jj,jk) ! loss of dco via scavenging
! augment scavenging of cobalt in shallow systems where lots of labile cobalt is
! present
               ikt = mbkt(ji,jj)
               zdep    = MIN( 1., 1000. / gdept_n(ji,jj,ikt+1) )
               zcopcp  = (coscav*10) * MAX( 0., zcofree(ji,jj,jk) - 100.E-9 ) * zdep
!! Dissolution of scavenged Co, define new ko2:
                o2a = MAX( 0. , trb(ji,jj,jk,jpoxy) - codiss_ox )
                o2a = MAX( 10.E-6 , o2a )
                ko2 = MAX( 0.5 , o2a**2 / (o2a**2 + ko2s**2) )
!                zlam1b = (20*coscav) * max( (1-ko2), kpar)
                zlam1b = (rdiss*coscav) * max( (1-ko2), kpar)
                zcodesor = zlam1b * xstep * trb(ji,jj,jk,jpsco)
                zcodesor = min( zcodesor, trb(ji,jj,jk,jpsco) )
                cscav(ji,jj,jk) = zcoscav
                cdiss(ji,jj,jk) = zcodesor
                tra(ji,jj,jk,jpdco) = tra(ji,jj,jk,jpdco)  - zcoscav + zcodesor - zcopcp
                tra(ji,jj,jk,jpsco) = tra(ji,jj,jk,jpsco)  + zcoscav - zcodesor
            END DO
         END DO
      END DO
      ENDIF ! ln_cobalt logical

      IF( ln_manganese ) THEN
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                zlmnt(ji,jj,jk) = mnlig ! set background Mn ligand concentration
                zmnfree(ji,jj,jk) = max( 0., ( trn(ji,jj,jk,jpdmn) - zlmnt(ji,jj,jk) ) )
            END DO
         END DO
      END DO
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                kpar_mn = MIN( 1. ,  etot(ji,jj,jk)**2 / (  etot(ji,jj,jk)**2 + kpars_mn**2))
                o2a_mn = MAX( 0. , trb(ji,jj,jk,jpoxy) - mnscav_ox )
                o2a_mn = MAX( 10.E-6 , o2a_mn )
!                ko2_mn = MAX(0.5, o2a_mn**2 / (o2a_mn**2 + ko2s_mn**2) )
                ko2_mn = MAX(0.25, o2a_mn**2 / (o2a_mn**2 + ko2s_mn**2) )
                kbact_mn = bacbio(ji,jj,jk)**2 / ( bacbio(ji,jj,jk)**2 + kbacts_mn**2)
                tgfun_mn = 0.5 * EXP( 0.0718 * tsn(ji,jj,jk,jp_tem) ) ! Q10 = 1.55 Lee and Fisher
                lam_mn = MAX( 0., (1 - kpar_mn) * ko2_mn * kbact_mn )
                lam = 0.1 * ( trb(ji,jj,jk,jpdmn)**4 / (trb(ji,jj,jk,jpdmn)**4 + 5.E-9**4 ) )
                minlam = MAX( MIN(minlam, 0.1) , 0.01 )
                minlam = 0.1
                IF( gdept_n(ji,jj,jk) > hmld(ji,jj) )  lam_mn = max( lam_mn , minlam )
                zlam1b = mnscav * tgfun_mn * lam_mn ! Mn scavenging rate
                zmnscav = zlam1b * xstep * zmnfree(ji,jj,jk) ! loss of dmn via scavenging
! Dissolution of scavenged Mn, define new ko2:
                o2a_mn = MAX( 0. , trb(ji,jj,jk,jpoxy) - mndiss_ox )
                o2a_mn = MAX( 10.E-6 , o2a_mn )
                ko2_mn = MAX( 0.5 , o2a_mn**2 / (o2a_mn**2 + ko2s_mn**2) )
                zlam1b = ( max( xdiss(ji,jj,jk), 0.5) * rdiss ) * mnscav * max( (1-ko2_mn), kpar_mn)
                zmndesor = zlam1b * xstep * trb(ji,jj,jk,jpsmn)
                zmndesor = min( zmndesor, trb(ji,jj,jk,jpsmn) )
               ikt = mbkt(ji,jj)
               zdep    = MIN( 1., 500. / gdept_n(ji,jj,ikt+1) )
               !zmnpcp  = (mnscav*10) * MAX( 0., zmnfree(ji,jj,jk) - 10.E-9 ) * zdep
                lam = trb(ji,jj,jk,jpdmn)**4 / (trb(ji,jj,jk,jpdmn)**4 + 5.E-9**4)
                zmnpcp  = (mnscav*10) * zmnfree(ji,jj,jk) * lam * zdep * xstep
                mscav(ji,jj,jk) = zmnscav
                mdiss(ji,jj,jk) = zmndesor
                tra(ji,jj,jk,jpdmn) = tra(ji,jj,jk,jpdmn)  - zmnscav + zmndesor - zmnpcp
                tra(ji,jj,jk,jpsmn) = tra(ji,jj,jk,jpsmn)  + zmnscav - zmndesor
            END DO
         END DO
      END DO
     ENDIF ! ln_manganese logical


      IF( ln_zinc ) THEN
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
!------------------------------------------------------------------------------
! Zn speciation
!------------------------------------------------------------------------------
! compute free zinc:
                zznl(ji,jj,jk) = znlig
                zzkeq           = 10**10.5   ! Lohan
                zzfesatur       = znlig  ! moles/L
                ztzn           = MAX( 0., trb(ji,jj,jk,jpdzn) )
                ! free Zn is the root of a 2nd order polynomial
                zZnfree(ji,jj,jk)  = ( -( 1. + zzfesatur * zzkeq - zzkeq * ztzn ) &
                     &             + SQRT( ( 1. + zzfesatur * zzkeq - zzkeq * ztzn)**2       &
                     &               + 4. * ztzn * zzkeq) ) / ( 2. * zzkeq )
                zZnfree(ji,jj,jk)  = MAX( 0., zZnfree(ji,jj,jk) ) ! Zn' in moles/L
                zznl(ji,jj,jk) = trb(ji,jj,jk,jpdzn) - zZnfree(ji,jj,jk)
                bzinc(ji,jj,jk) = zZnfree(ji,jj,jk) / 2.1 ! To derive Zn2+ from Zn' using 
                                                          ! alpha from Ellwood
                                                          ! and Van den Berg,
                                                          ! 2000
            END DO
         END DO
      END DO
               IF(znkdp > 0) THEN ! if a non zero value for the partition coefficient is used
               IF(znkdg > 0) THEN
     DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

!------------------------------------------------------------------------------
! Zn scavenging
!------------------------------------------------------------------------------
! use partition coefficient approach to calculate dissolved and adsorbed Zn:
               !IF( ln_kdexp ) THEN ! Explicit calculation of fluxes
               ztrc1 = ( trb(ji,jj,jk,jppoc) * 1e6 ) * znkdp
               ztrc2 = ( trb(ji,jj,jk,jpgoc) * 1e6 ) * znkdg
               fzn = ( ztrc1 + ztrc2 ) / ( 1. +  ztrc1 + ztrc2 )
               fzn = MIN( 1., MAX( 0., fzn ))
               tzn = trb(ji,jj,jk,jpszp) + trb(ji,jj,jk,jpszg) + rtrn
               zdenom1 = ( ztrc1 / ( ztrc1 + 1. ) * zZnfree(ji,jj,jk) ) - trb(ji,jj,jk,jpszp)
               zdenom2 = ( ztrc2 / ( ztrc2 + 1. ) * zZnfree(ji,jj,jk) ) - trb(ji,jj,jk,jpszg)
               znflux(ji,jj,jk) = zdenom1 + zdenom2
               tra(ji,jj,jk,jpdzn) = tra(ji,jj,jk,jpdzn) - znflux(ji,jj,jk)
               tra(ji,jj,jk,jpszp) = tra(ji,jj,jk,jpszp) + zdenom1
               tra(ji,jj,jk,jpszg) = tra(ji,jj,jk,jpszg) + zdenom2
            END DO
         END DO
      END DO
               ELSE
     DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztrc1 = ( trb(ji,jj,jk,jppoc) * 1e6 ) * znkdp
               fzn = ztrc1 / ( 1. +  ztrc1 )
               fzn = MIN( 1., MAX( 0., fzn ))
               tzn = trb(ji,jj,jk,jpszp) + rtrn
               znflux(ji,jj,jk) = ( ztrc1 / ( ztrc1 + 1. ) * zZnfree(ji,jj,jk) ) - trb(ji,jj,jk,jpszp)
               tra(ji,jj,jk,jpdzn) = tra(ji,jj,jk,jpdzn) - znflux(ji,jj,jk)
               tra(ji,jj,jk,jpszp) = tra(ji,jj,jk,jpszp) + znflux(ji,jj,jk)
            END DO
         END DO
      END DO
                ENDIF ! end of if znkdg > 0 logical
                ENDIF ! if kdznp > 0 logical
      ENDIF ! ln_zinc logical


      IF( ln_copper ) THEN !Cam

        IF( ln_culigvar ) THEN  !Cam key for variable ligands
                zcutotlig(:,:,:) =  0.09 * trb(:,:,:,jpdoc) * 1E6 + culig * 1E9
                zcutotlig(:,:,:) =  MIN( zcutotlig(:,:,:), 10. )
                zcutotlig(:,:,:) = zcutotlig(:,:,:) * 1E-9
        ELSE
                zcutotlig(:,:,:) = culig
        ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
!------------------------------------------------------------------------------
! Cu speciation
!------------------------------------------------------------------------------
! compute free copper:
                zcul(ji,jj,jk) = zcutotlig(ji,jj,jk)
                zcukeq           = 10**13.5   ! Arthur Gourain 
                zcufesatur       = zcul(ji,jj,jk)  ! moles/L
                ztcu           = MAX( 0., trb(ji,jj,jk,jpdcu) )
                ! free Zn is the root of a 2nd order polynomial
                zCufree(ji,jj,jk)  = ( -( 1. + zcufesatur * zcukeq - zcukeq * ztcu ) &
                     &             + SQRT( ( 1. + zcufesatur * zcukeq - zcukeq * ztcu)**2       &
                     &               + 4. * ztcu * zcukeq) ) / ( 2. * zcukeq )
                zCufree(ji,jj,jk)  = MAX( 0., zCufree(ji,jj,jk) ) ! Zn' in moles/L
                zcul(ji,jj,jk) = trb(ji,jj,jk,jpdcu) - zCufree(ji,jj,jk)
!                bcopper(ji,jj,jk) = zCufree(ji,jj,jk) 
                bcopper(ji,jj,jk) = trb(ji,jj,jk,jpdcu)      ! Cu2+ = 0.1% * Cu'

            END DO
         END DO
      END DO
               IF(cukdp > 0) THEN ! if a non zero value for the partition coefficient is used
               IF(cukdg > 0) THEN
     DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

!------------------------------------------------------------------------------
! Cu scavenging
!------------------------------------------------------------------------------
! use partition coefficient approach to calculate dissolved and adsorbed Zn:
               !IF( ln_kdexp ) THEN ! Explicit calculation of fluxes

               ztrc1 = ( trb(ji,jj,jk,jppoc) * 1e6 ) * cukdp
               ztrc2 = ( trb(ji,jj,jk,jpgoc) * 1e6 ) * cukdg

!                IF (ln_sulfprecip )  THEN
!                zoxy2 = trb(ji,jj,jk,jpoxy)
!                sulf_o2 = zoxy2 /( zoxy2 + k_o2sulf)
!                sulf_o2 = MAX(0., 1.-sulf_o2) ! modulate cukd by o2
!                ztrc1 = ztrc1 * sulf_o2 + ztrc1
!                ztrc2 = ztrc2 * sulf_o2 + ztrc2
!                ztrc1 = max(ztrc1, min_cukd1)
!                ztrc2 = max(ztrc2, min_cukd2)
!                ENDIF

               fcu = ( ztrc1 + ztrc2 ) / ( 1. +  ztrc1 + ztrc2 )
               fcu = MIN( 1., MAX( 0., fcu ))
               tcu = trb(ji,jj,jk,jpscup) + trb(ji,jj,jk,jpscug) + rtrn
               zdenom1 = ( ztrc1 / ( ztrc1 + 1. ) * zCufree(ji,jj,jk) ) - trb(ji,jj,jk,jpscup)
               zdenom2 = ( ztrc2 / ( ztrc2 + 1. ) * zCufree(ji,jj,jk) ) - trb(ji,jj,jk,jpscug)
               cuflux(ji,jj,jk) = zdenom1 + zdenom2
!               IF (ln_irrevscav) THEN !Cam irreersible scavenging
!               zlam1b = 3.e-5 + xlamdust * zdust(ji,jj,jk) + xlam * ztrc
!               zscavecu = zCufree(ji,jj,jk)  * zlam1b * xstep
!               tra(ji,jj,jk,jpdcu) = tra(ji,jj,jk,jpdcu) - zscavecu
!               tra(ji,jj,jk,jpcup) = tra(ji,jj,jk,jpcup) + zscavecu * zdenom1
!               tra(ji,jj,jk,jpcug) = tra(ji,jj,jk,jpcug) + zscavecu * zdenom2

!              ELSE

               tra(ji,jj,jk,jpdcu) = tra(ji,jj,jk,jpdcu) - cuflux(ji,jj,jk)
               tra(ji,jj,jk,jpscup) = tra(ji,jj,jk,jpscup) + zdenom1
               tra(ji,jj,jk,jpscug) = tra(ji,jj,jk,jpscug) + zdenom2
!               ENDIF
            END DO
         END DO
      END DO

               ELSE
     DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztrc = ( trb(ji,jj,jk,jppoc) * 1e6 ) * cukdp
               fcu = ztrc1 / ( 1. +  ztrc1 )
               fcu = MIN( 1., MAX( 0., fcu ))
               tcu = trb(ji,jj,jk,jpscup) + rtrn
               cuflux(ji,jj,jk) = ( ztrc1 / ( ztrc1 + 1. ) * zCufree(ji,jj,jk) ) - trb(ji,jj,jk,jpscup)
               tra(ji,jj,jk,jpdcu) = tra(ji,jj,jk,jpdcu) - cuflux(ji,jj,jk)
               tra(ji,jj,jk,jpscup) = tra(ji,jj,jk,jpscup) + cuflux(ji,jj,jk)
            END DO
         END DO
      END DO
               ENDIF
               ENDIF
        ENDIF !Cam ln_copper
 






      !  Output of some diagnostics variables
      !     ---------------------------------
      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
            zrfact2 = 1.e3 * rfact2r  ! conversion from mol/L/timestep into mol/m3/s
            IF( iom_use("Fe3")    )  CALL iom_put("Fe3"    , zFe3   (:,:,:)       * tmask(:,:,:) )   ! Fe3+
            IF( iom_use("FeL1")   )  CALL iom_put("FeL1"   , zFeL1  (:,:,:)       * tmask(:,:,:) )   ! FeL1
            IF( iom_use("TL1")    )  CALL iom_put("TL1"    , zTL1   (:,:,:)       * tmask(:,:,:) )   ! TL1
            IF( iom_use("Totlig") )  CALL iom_put("Totlig" , ztotlig(:,:,:)       * tmask(:,:,:) )   ! TL
            IF( iom_use("Biron")  )  CALL iom_put("Biron"  , biron  (:,:,:)  * 1e9 * tmask(:,:,:) )   ! biron
            IF( iom_use("FESCAV") )  CALL iom_put("FESCAV" , zscav3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
            IF( iom_use("FECOLL") )  CALL iom_put("FECOLL" , zcoll3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
            IF( iom_use("LGWCOLL"))  CALL iom_put("LGWCOLL", zlcoll3d(:,:,:) * 1e9 * tmask(:,:,:) * zrfact2 )
IF (ln_copper) THEN 
         IF( iom_use("CUPR"))  CALL iom_put("CUPR"  , zCufree(:,:,:) * 1e9 * tmask(:,:,:) ) ! Free copper
         IF( iom_use("CuL"))   CALL iom_put("CuL"   , zcul  (:,:,:) * 1e9 * tmask(:,:,:) ) ! Zn-L
         IF( iom_use("CUSCAV")) CALL iom_put("CUSCAV" , cuflux(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 ) ! Cu scav flux
         IF( iom_use("Totculig")) CALL iom_put("Totculig" , zcutotlig(:,:,:)  * 1e9 * tmask(:,:,:)  ) ! Cu Ligands
ENDIF
IF (ln_zinc) THEN 
         IF( iom_use("ZNPR"))  CALL iom_put("ZNPR"  , zZnfree(:,:,:) * 1e9 * tmask(:,:,:) ) ! Free zinc
         IF( iom_use("ZNL"))   CALL iom_put("ZNL"   , zznl  (:,:,:) * 1e9 * tmask(:,:,:) ) ! Zn-L
         IF( iom_use("ZNSCAV")) CALL iom_put("ZNSCAV" , znflux(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 ) ! Zn scav flux
ENDIF
IF (ln_cobalt) THEN
         IF( iom_use("Copr")) CALL iom_put("Copr"   , zcofree(:,:,:) * tmask(:,:,:) * 1.e9 ) !FreeCobgalt in pmol/L
         IF( iom_use("CoLig")) CALL iom_put("CoLig"  , zlcot(:,:,:) * tmask(:,:,:) * 1.e9) ! Co ligands in pmol/L
         IF( iom_use("CoSCAV")) CALL iom_put("CoSCAV" , cscav(:,:,:) * tmask(:,:,:) * zrfact2 ) !mmol m-2 s-1
         IF( iom_use("CoDISS")) CALL iom_put("CoDISS" , cdiss(:,:,:) * tmask(:,:,:) * zrfact2 ) ! mmol m-2 s-1
ENDIF
IF (ln_manganese) THEN
         IF( iom_use("Mnpr")) CALL iom_put("Mnpr"   , zmnfree(:,:,:) * tmask(:,:,:) * 1.e9 ) !FreeCobgalt in pmol/L
         IF( iom_use("MnLig")) CALL iom_put("MnLig"  , zlmnt(:,:,:) * tmask(:,:,:) * 1.e9) ! Co ligands in pmol/L
         IF( iom_use("MnSCAV")) CALL iom_put("MnSCAV" , mscav(:,:,:) * tmask(:,:,:) * zrfact2 ) !mmol m-2 s-1
         IF( iom_use("MnDISS")) CALL iom_put("MnDISS" , mdiss(:,:,:) * tmask(:,:,:) * zrfact2 ) ! mmol m-2 s-1
ENDIF

         ENDIF
      ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('fechem')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_fechem')
      !
   END SUBROUTINE p4z_fechem


   SUBROUTINE p4z_fechem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_fechem_init  ***
      !!
      !! ** Purpose :   Initialization of iron chemistry parameters
      !!
      !! ** Method  :   Read the nampisfer namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisfer
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer 
      !!
      NAMELIST/nampisfer/ ln_ligvar, ln_culigvar, xlam1, xlamdust, ligand, kfep &
&                          , culig, cukdp, cukdg, znlig, znkdp, znkdg , coscav, colig, fnan, ko2s, kbacts, kpars, coscav_ox, &
      &                    codiss_ox,  mnscav, mnlig, ko2s_mn, kbacts_mn, kpars_mn, mnscav_ox, &
      &                    mndiss_ox, rdiss 
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_rem_init : Initialization of iron chemistry parameters'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )            ! Namelist nampisfer in reference namelist : Pisces iron chemistry
      READ  ( numnatp_ref, nampisfer, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisfer in reference namelist', lwp )
      REWIND( numnatp_cfg )            ! Namelist nampisfer in configuration namelist : Pisces iron chemistry
      READ  ( numnatp_cfg, nampisfer, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisfer in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, nampisfer )

      IF(lwp) THEN                     ! control print
         WRITE(numout,*) '   Namelist : nampisfer'
         WRITE(numout,*) '      variable concentration of ligand          ln_ligvar    =', ln_ligvar
!         WRITE(numout,*) '    irreversible Cu scavenging                  ln_irrevscav =', ln_irrevscav
         WRITE(numout,*) '    variable concentration of cu ligand         ln_culigvar  =', ln_culigvar
         WRITE(numout,*) '      scavenging rate of Iron                   xlam1        =', xlam1
         WRITE(numout,*) '      scavenging rate of Iron by dust           xlamdust     =', xlamdust
         WRITE(numout,*) '      ligand concentration in the ocean         ligand       =', ligand
         WRITE(numout,*) '      rate constant for nanoparticle formation  kfep         =', kfep
         WRITE(numout,*) '  Cu ligand concentration                       culig        =', culig
         WRITE(numout,*) ' Kd for Cu reversible scavenging                cukdp        =', cukdp
         WRITE(numout,*) ' Kd for Cu reversible scavenging                cukdg        =', cukdg
         WRITE(numout,*) '  Zn ligand concentration                       znlig =', znlig
         WRITE(numout,*) ' Kd for Zn reversible scavenging                znkdg =', znkdp
         WRITE(numout,*) ' Kd for Zn reversible scavenging                znkdp =', znkdg
         WRITE(numout,*) '    scavenging rate of Cobalt                   coscav    =', coscav
         WRITE(numout,*) '    Background Cobalt ligand concentration      colig     =', colig
         WRITE(numout,*) '    Maximum nano fraction for Cobalt ligands    fnan      =', fnan
         WRITE(numout,*) '  Half saturation constant for o2 - Co scav     ko2s      =', ko2s
         WRITE(numout,*) '  Half saturation constant for bact - Co scav   kbacts    =', kbacts
         WRITE(numout,*) '  Half saturation constant for par - Co scav    kpars     =', kpars
         WRITE(numout,*) '  O2 threshold for Co scavenging                coscav_ox =', coscav_ox
         WRITE(numout,*) '  O2 threshold for Co dissolution               codiss_ox =', codiss_ox
         WRITE(numout,*) '    scavenging rate of Manganes                 mnscav    =', mnscav
         WRITE(numout,*) '    Background Manganese ligand concentration   mnlig     =', mnlig
         WRITE(numout,*) '  Half saturation constant for o2 - Mn scav     ko2s_mn   =', ko2s_mn
         WRITE(numout,*) '  Half saturation constant for bact - Mn scav   kbacts_mn =', kbacts_mn
         WRITE(numout,*) '  Half saturation constant for par - Mn scav    kpars_mn  =', kpars_mn
         WRITE(numout,*) '  O2 threshold for Mn scavenging                mnscav_ox =', mnscav_ox
         WRITE(numout,*) '  O2 threshold for mn dissolution               mndiss_ox =', mndiss_ox
         WRITE(numout,*) ' relative rate of Co/Mn oxide dissolution           rdiss =', rdiss
      ENDIF
      ! 
   END SUBROUTINE p4z_fechem_init
   
   !!======================================================================
END MODULE p4zfechem
