MODULE p4zlim
   !!======================================================================
   !!                         ***  MODULE p4zlim  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota 
   !!----------------------------------------------------------------------
   !!   p4z_lim        :   Compute the nutrients limitation terms 
   !!   p4z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------
   USE oce_trc         ! Shared ocean-passive tracers variables
   USE trc             ! Tracers defined
   USE sms_pisces      ! PISCES variables
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC p4z_lim    
   PUBLIC p4z_lim_init    
   PUBLIC p4z_lim_alloc

   !! * Shared module variables
   REAL(wp), PUBLIC ::  concnno3    !:  NO3, PO4 half saturation   
   REAL(wp), PUBLIC ::  concdno3    !:  Phosphate half saturation for diatoms  
   REAL(wp), PUBLIC ::  concnnh4    !:  NH4 half saturation for phyto  
   REAL(wp), PUBLIC ::  concdnh4    !:  NH4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concnfer    !:  Iron half saturation for nanophyto 
   REAL(wp), PUBLIC ::  concdfer    !:  Iron half saturation for diatoms  
   REAL(wp), PUBLIC ::  concbno3    !:  NO3 half saturation  for bacteria 
   REAL(wp), PUBLIC ::  concbnh4    !:  NH4 half saturation for bacteria
   REAL(wp), PUBLIC ::  xsizedia    !:  Minimum size criteria for diatoms
   REAL(wp), PUBLIC ::  xsizephy    !:  Minimum size criteria for nanophyto
   REAL(wp), PUBLIC ::  xsizern     !:  Size ratio for nanophytoplankton
   REAL(wp), PUBLIC ::  xsizerd     !:  Size ratio for diatoms
   REAL(wp), PUBLIC ::  xksi1       !:  half saturation constant for Si uptake 
   REAL(wp), PUBLIC ::  xksi2       !:  half saturation constant for Si/C 
   REAL(wp), PUBLIC ::  xkdoc       !:  2nd half-sat. of DOC remineralization  
   REAL(wp), PUBLIC ::  concbfe     !:  Fe half saturation for bacteria 
   REAL(wp), PUBLIC ::  oxymin      !:  half saturation constant for anoxia
   REAL(wp), PUBLIC ::  qnfelim     !:  optimal Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qdfelim     !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  caco3r      !:  mean rainratio 
   REAL(wp), PUBLIC ::  kscun      !: Cam half saturation for Cu
   REAL(wp), PUBLIC ::  kscud
   REAL(wp), PUBLIC ::  ksznn
   REAL(wp), PUBLIC ::  ksznd
   REAL(wp), PUBLIC ::  kscon , k1con
   REAL(wp), PUBLIC ::  kscod, k1cod
   REAL(wp), PUBLIC ::  ksmnn 
   REAL(wp), PUBLIC ::  ksmnd
   LOGICAL , PUBLIC :: ln_colim
   LOGICAL , PUBLIC :: ln_colim_int
   REAL(wp), PUBLIC ::  qmnmin
   REAL(wp), PUBLIC ::  mnchln
   REAL(wp), PUBLIC ::  mnchld
   REAL(wp), PUBLIC ::  qcoznmin
   REAL(wp), PUBLIC ::  qcozncan
   REAL(wp), PUBLIC ::  qcozncad
   REAL(wp), PUBLIC ::  kco2
   REAL(wp), PUBLIC ::  qcoznapn
   REAL(wp), PUBLIC ::  qcoznapd
   REAL(wp), PUBLIC ::   kpo4

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanono3   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatno3   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanonh4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatnh4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanopo4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatpo4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimphy    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdia    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnfe    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdfe    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimsi     !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimbac    !: ??
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimbacl   !: ??
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concdfe    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concnfe    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concdcu    !: Cam Copper model
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concncu    !: 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concdzn    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concnzn    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concdco    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concnco    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concdmn    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concnmn    !: ???
   !! colimitation for diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   zmnued
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   zcoznued
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   coznlimd
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   mnlimd
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   colimd
   !!
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   zmnuen
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   zcoznuen
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   coznlimn
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   mnlimn
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   colimn

   ! Coefficient for iron limitation
   REAL(wp) ::  xcoef1   = 0.0016  / 55.85  
   REAL(wp) ::  xcoef2   = 1.21E-5 * 14. / 55.85 / 7.625 * 0.5 * 1.5
   REAL(wp) ::  xcoef3   = 1.15E-4 * 14. / 55.85 / 7.625 * 0.5 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zlim.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_lim( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!              for the various phytoplankton species
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt, knt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zlim1, zlim2, zlim3, zlim4, zno3, zferlim
      REAL(wp) ::   zconcd, zconcd2, zconcn, zconcn2
      REAL(wp) ::   z1_trbdia, z1_trbphy, ztem1, ztem2, zetot1, zetot2
      REAL(wp) ::   zdenom, zratio, zironmin
      REAL(wp) ::   zconc1d, zconc1dnh4, zconc0n, zconc0nnh4   
      REAL(wp) ::   pco2, po4, zchlc, zqmn, zqcozn, zmu, mnued, coznued
      REAL(wp) ::   mnuen, coznuen
     REAL(wp), DIMENSION(jpi,jpj,jpk) :: qfepsn, qfeno3n, qferesn, xlimnnit
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: qfepsd, qfeno3d, qferesd, xlimdnit
      !!---------------------------------------------------------------------
      !
      qfepsn(:,:,:) = 0._wp  ;      qfeno3n(:,:,:) = 0._wp
      qferesn(:,:,:) = 0._wp ;      qfepsd(:,:,:) = 0._wp
      qfeno3d(:,:,:) = 0._wp ;      qferesd(:,:,:) = 0._wp
      xnanonh4(:,:,:) = 0._wp ;     xdiatnh4(:,:,:) = 0._wp
      xnanono3(:,:,:) = 0._wp ;     xdiatno3(:,:,:) = 0._wp
      xlimnnit(:,:,:) = 0._wp ;     xlimdnit(:,:,:) = 0._wp
      !
      IF( ln_timing )   CALL timing_start('p4z_lim')
      !
      IF (ln_cobalt) THEN
               k1con = kscon * 1000
               k1cod = kscod * 1000
      ENDIF
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               
               ! Tuning of the iron concentration to a minimum level that is set to the detection limit
               !-------------------------------------
               zno3    = trb(ji,jj,jk,jpno3) / 40.e-6
               zferlim = MAX( 3e-11 * zno3 * zno3, 5e-12 )
               zferlim = MIN( zferlim, 7e-11 )
               trb(ji,jj,jk,jpfer) = MAX( trb(ji,jj,jk,jpfer), zferlim )

               ! Computation of a variable Ks for iron on diatoms taking into account
               ! that increasing biomass is made of generally bigger cells
               !------------------------------------------------
               zconcd   = MAX( 0.e0 , trb(ji,jj,jk,jpdia) - xsizedia )
               zconcd2  = trb(ji,jj,jk,jpdia) - zconcd
               zconcn   = MAX( 0.e0 , trb(ji,jj,jk,jpphy) - xsizephy )
               zconcn2  = trb(ji,jj,jk,jpphy) - zconcn
               z1_trbphy   = 1. / ( trb(ji,jj,jk,jpphy) + rtrn )
               z1_trbdia   = 1. / ( trb(ji,jj,jk,jpdia) + rtrn )

               concdfe(ji,jj,jk) = MAX( concdfer, ( zconcd2 * concdfer + concdfer * xsizerd * zconcd ) * z1_trbdia )
               zconc1d           = MAX( concdno3, ( zconcd2 * concdno3 + concdno3 * xsizerd * zconcd ) * z1_trbdia )
               zconc1dnh4        = MAX( concdnh4, ( zconcd2 * concdnh4 + concdnh4 * xsizerd * zconcd ) * z1_trbdia )

               concnfe(ji,jj,jk) = MAX( concnfer, ( zconcn2 * concnfer + concnfer * xsizern * zconcn ) * z1_trbphy )
               zconc0n           = MAX( concnno3, ( zconcn2 * concnno3 + concnno3 * xsizern * zconcn ) * z1_trbphy )
               zconc0nnh4        = MAX( concnnh4, ( zconcn2 * concnnh4 + concnnh4 * xsizern * zconcn ) * z1_trbphy )
               IF (ln_copper) THEN
               concdcu(ji,jj,jk) = MAX( kscud,    ( zconcd2 * kscud    + kscud * xsizerd * zconcd ) * z1_trbdia )
               concncu(ji,jj,jk) = MAX( kscun,    ( zconcn2 * kscun    + kscun * xsizern * zconcn ) * z1_trbphy )
               ENDIF
               IF (ln_zinc) THEN
               concdzn(ji,jj,jk) = MAX( ksznd,    ( zconcd2 * ksznd    + ksznd * xsizerd * zconcd ) * z1_trbdia )
               concnzn(ji,jj,jk) = MAX( ksznn,    ( zconcn2 * ksznn    + ksznn * xsizern * zconcn ) * z1_trbphy )
               ENDIF
               IF (ln_cobalt) THEN
               concdco(ji,jj,jk) = MAX( k1cod,    ( zconcd2 * k1cod    + k1cod * xsizerd * zconcd ) * z1_trbdia )
               concnco(ji,jj,jk) = MAX( k1con,    ( zconcn2 * k1con    + k1con * xsizern * zconcn ) * z1_trbphy )
               ENDIF
               IF (ln_manganese) THEN
               concdmn(ji,jj,jk) = MAX( ksmnd,    ( zconcd2 * ksmnd    + ksmnd * xsizerd * zconcd ) * z1_trbdia )
               concnmn(ji,jj,jk) = MAX( ksmnn,    ( zconcn2 * ksmnn    + ksmnn * xsizern * zconcn ) * z1_trbphy )
               ENDIF

               ! Michaelis-Menten Limitation term for nutrients Small bacteria
               ! -------------------------------------------------------------
               zdenom = 1. /  ( concbno3 * concbnh4 + concbnh4 * trb(ji,jj,jk,jpno3) + concbno3 * trb(ji,jj,jk,jpnh4) )
               xnanono3(ji,jj,jk) = trb(ji,jj,jk,jpno3) * concbnh4 * zdenom
               xnanonh4(ji,jj,jk) = trb(ji,jj,jk,jpnh4) * concbno3 * zdenom
               !
               zlim1    = xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk)
               zlim2    = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + concbnh4 )
               zlim3    = trb(ji,jj,jk,jpfer) / ( concbfe + trb(ji,jj,jk,jpfer) )
               zlim4    = trb(ji,jj,jk,jpdoc) / ( xkdoc   + trb(ji,jj,jk,jpdoc) )
               xlimbacl(ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
               xlimbac (ji,jj,jk) = MIN( zlim1, zlim2, zlim3 ) * zlim4

               ! Michaelis-Menten Limitation term for nutrients Small flagellates
               ! -----------------------------------------------
               zdenom = 1. /  ( zconc0n * zconc0nnh4 + zconc0nnh4 * trb(ji,jj,jk,jpno3) + zconc0n * trb(ji,jj,jk,jpnh4) )
               xnanono3(ji,jj,jk) = trb(ji,jj,jk,jpno3) * zconc0nnh4 * zdenom
               xnanonh4(ji,jj,jk) = trb(ji,jj,jk,jpnh4) * zconc0n    * zdenom
               !
               zlim1    = xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk)
               zlim2    = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + zconc0nnh4 )
               xlimnnit(ji,jj,jk) = zlim1
               zratio   = trb(ji,jj,jk,jpnfe) * z1_trbphy 
               zironmin = xcoef1 * trb(ji,jj,jk,jpnch) * z1_trbphy + xcoef2 * zlim1 + xcoef3 * xnanono3(ji,jj,jk)
               zlim3    = MAX( 0.,( zratio - zironmin ) / qnfelim )
               xnanopo4(ji,jj,jk) = zlim2
               xlimnfe (ji,jj,jk) = MIN( 1., zlim3 )
               qfepsn  (ji,jj,jk) = xcoef1 * trb(ji,jj,jk,jpnch) * z1_trbphy
               qfeno3n (ji,jj,jk) = xcoef3 * xnanono3(ji,jj,jk)
               qferesn (ji,jj,jk) = xcoef2 * zlim1
               xlimphy (ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
               IF (ln_manganese .AND. ln_zinc .AND. ln_cobalt) THEN
               IF ( ln_colim ) THEN
               po4     = trb(ji,jj,jk,jppo4) * po4r + rtrn
!               pco2    = pco2s(ji,jj) + rtrn 
               pco2    = co2aq(ji,jj,jk) 
               zchlc   = trb(ji,jj,jk,jpnch) * z1_trbphy
               zmu     =  mu_nm(ji,jj,jk) + rtrn
               zqmn    = trb(ji,jj,jk,jpmnn) * z1_trbphy
               zqcozn   = ( trb(ji,jj,jk,jpznn) + 1.E-3 * trb(ji,jj,jk,jpcon) ) * z1_trbphy
               zmnuen(ji,jj,jk)   = (1/rday) / ( qmnmin + ( mnchln * zchlc ) )
               zcoznuen(ji,jj,jk) = (1/rday) / ( qcoznmin + ( qcozncan * ( kco2**2 /  &
               &                    ( kco2**2 + pco2**2 ) ) )   &
               &                  + ( qcoznapn * ( kpo4**3 / ( kpo4**3 + po4**3 ) ) ) )
               !! UE is in Q per second, mulitply by quota to get per second,
               !then divide by mu_max to get unitless
               mnlimn(ji,jj,jk)   = MIN( 1., ( zmnuen(ji,jj,jk)  * zqmn ) / zmu )
               coznlimn(ji,jj,jk) = MIN( 1., ( zcoznuen(ji,jj,jk) * zqcozn ) / zmu )
               colimn(ji,jj,jk)   = mnlimn(ji,jj,jk) * coznlimn(ji,jj,jk)
               IF ( ln_colim_int ) THEN
               xlimphy (ji,jj,jk) = MIN( xlimphy(ji,jj,jk), mnlimn(ji,jj,jk), coznlimn(ji,jj,jk) )
               ELSE
               mnlimn(ji,jj,jk)   = 1.0
               coznlimn(ji,jj,jk) = 1.0
               colimn(ji,jj,jk)   = 1.0
               ENDIF
               ENDIF
               ENDIF
               !
               !   Michaelis-Menten Limitation term for nutrients Diatoms
               !   ----------------------------------------------
               zdenom   = 1. / ( zconc1d * zconc1dnh4 + zconc1dnh4 * trb(ji,jj,jk,jpno3) + zconc1d * trb(ji,jj,jk,jpnh4) )
               xdiatno3(ji,jj,jk) = trb(ji,jj,jk,jpno3) * zconc1dnh4 * zdenom
               xdiatnh4(ji,jj,jk) = trb(ji,jj,jk,jpnh4) * zconc1d    * zdenom
               !
               zlim1    = xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk)
               xlimdnit(ji,jj,jk) = zlim1
               zlim2    = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + zconc1dnh4  )
               zlim3    = trb(ji,jj,jk,jpsil) / ( trb(ji,jj,jk,jpsil) + xksi(ji,jj) )
               zratio   = trb(ji,jj,jk,jpdfe) * z1_trbdia
               zironmin = xcoef1 * trb(ji,jj,jk,jpdch) * z1_trbdia + xcoef2 * zlim1 + xcoef3 * xdiatno3(ji,jj,jk)
               zlim4    = MAX( 0., ( zratio - zironmin ) / qdfelim )
               xdiatpo4(ji,jj,jk) = zlim2
               xlimdfe (ji,jj,jk) = MIN( 1., zlim4 )
               qfepsd  (ji,jj,jk) = xcoef1 * trb(ji,jj,jk,jpdch) * z1_trbdia
               qfeno3d (ji,jj,jk) = xcoef3 * xdiatno3(ji,jj,jk)
               qferesd (ji,jj,jk) = xcoef2 * zlim1
               xlimdia (ji,jj,jk) = MIN( zlim1, zlim2, zlim3, zlim4 )
               xlimsi  (ji,jj,jk) = MIN( zlim1, zlim2, zlim4 )
               IF (ln_manganese .AND. ln_zinc .AND. ln_cobalt) THEN
               IF ( ln_colim ) THEN
               po4     = trb(ji,jj,jk,jppo4) * po4r + rtrn
!               pco2    = pco2s(ji,jj) + rtrn
               pco2    = co2aq(ji,jj,jk)
               zchlc   = trb(ji,jj,jk,jpdch) * z1_trbdia
               zmu     =  mu_dm(ji,jj,jk) + rtrn
               zqmn    = trb(ji,jj,jk,jpmnd) * z1_trbdia
               zqcozn   = ( trb(ji,jj,jk,jpznd) + 1.E-3 * trb(ji,jj,jk,jpcod) ) * z1_trbdia
               zmnued(ji,jj,jk)   = (1/rday) / ( qmnmin + ( mnchld * zchlc ) )
               zcoznued(ji,jj,jk) = (1/rday) / ( qcoznmin + ( qcozncad * ( kco2**2 /  &  
               &                    ( kco2**2 + pco2**2 ) ) )   &
               &                  + ( qcoznapd * ( kpo4**3 / ( kpo4**3 + po4**3 ) ) ) ) 
               !! UE is in Q per second, mulitply by quota to get per second,
               !then divide by mu_max to get unitless
               mnlimd(ji,jj,jk)   = MIN( 1., ( zmnued(ji,jj,jk)  * zqmn ) / zmu )
               coznlimd(ji,jj,jk) = MIN( 1., ( zcoznued(ji,jj,jk) * zqcozn ) / zmu )
               colimd(ji,jj,jk)   = mnlimd(ji,jj,jk) * coznlimd(ji,jj,jk)
               IF ( ln_colim_int ) THEN
               xlimdia (ji,jj,jk) = MIN( xlimdia(ji,jj,jk), mnlimd(ji,jj,jk), coznlimd(ji,jj,jk) )
               ELSE
               mnlimd(ji,jj,jk)   = 1.0
               coznlimd(ji,jj,jk) = 1.0
               colimd(ji,jj,jk)   = 1.0
               ENDIF
               ENDIF
               ENDIF
           END DO
         END DO
      END DO

      ! --------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zlim1 =  ( trb(ji,jj,jk,jpno3) * concnnh4 + trb(ji,jj,jk,jpnh4) * concnno3 )    &
                  &   / ( concnno3 * concnnh4 + concnnh4 * trb(ji,jj,jk,jpno3) + concnno3 * trb(ji,jj,jk,jpnh4) ) 
               zlim2  = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + concnnh4 )
               zlim3  = trb(ji,jj,jk,jpfer) / ( trb(ji,jj,jk,jpfer) +  5.E-11   )
               ztem1  = MAX( 0., tsn(ji,jj,jk,jp_tem) )
               ztem2  = tsn(ji,jj,jk,jp_tem) - 10.
               zetot1 = MAX( 0., etot_ndcy(ji,jj,jk) - 1.) / ( 4. + etot_ndcy(ji,jj,jk) ) 
               zetot2 = 30. / ( 30. + etot_ndcy(ji,jj,jk) ) 

               xfracal(ji,jj,jk) = caco3r * MIN( zlim1, zlim2, zlim3 )                  &
                  &                       * ztem1 / ( 0.1 + ztem1 )                     &
                  &                       * MAX( 1., trb(ji,jj,jk,jpphy) * 1.e6 / 2. )  &
                  &                       * zetot1 * zetot2               &
                  &                       * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )         &
                  &                       * MIN( 1., 50. / ( hmld(ji,jj) + rtrn ) )
               xfracal(ji,jj,jk) = MIN( 0.8 , xfracal(ji,jj,jk) )
               xfracal(ji,jj,jk) = MAX( 0.02, xfracal(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! denitrification factor computed from O2 levels
               nitrfac(ji,jj,jk) = MAX(  0.e0, 0.4 * ( 6.e-6  - trb(ji,jj,jk,jpoxy) )    &
                  &                                / ( oxymin + trb(ji,jj,jk,jpoxy) )  )
               nitrfac(ji,jj,jk) = MIN( 1., nitrfac(ji,jj,jk) )
               !
               ! denitrification factor computed from NO3 levels
               nitrfac2(ji,jj,jk) = MAX( 0.e0,       ( 1.E-6 - trb(ji,jj,jk,jpno3) )  &
                  &                                / ( 1.E-6 + trb(ji,jj,jk,jpno3) ) )
               nitrfac2(ji,jj,jk) = MIN( 1., nitrfac2(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN        ! save output diagnostics
        IF( iom_use( "xfracal" ) )   CALL iom_put( "xfracal", xfracal(:,:,:) * tmask(:,:,:) )  ! euphotic layer deptht
        IF( iom_use( "LNnut"   ) )   CALL iom_put( "LNnut"  , xlimphy(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LDnut"   ) )   CALL iom_put( "LDnut"  , xlimdia(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LNFe"    ) )   CALL iom_put( "LNFe"   , xlimnfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LDFe"    ) )   CALL iom_put( "LDFe"   , xlimdfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "MNUED"   ) )   CALL iom_put( "MNUED"  , zmnued(:,:,:) * tmask(:,:,:)  )  ! Iron limitation term
        IF( iom_use( "COZNUED" ) )   CALL iom_put( "COZNUED", zcoznued(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "MNLIMD"  ) )   CALL iom_put( "MNLIMD" , mnlimd(:,:,:) * tmask(:,:,:)  )  !
        IF( iom_use( "COZNLIMD") )   CALL iom_put( "COZNLIMD", coznlimd(:,:,:) * tmask(:,:,:) )  !
        IF( iom_use( "COLIMD"  ) )   CALL iom_put( "COLIMD" , colimd(:,:,:) * tmask(:,:,:)  )  !
        IF( iom_use( "MNUED"   ) )   CALL iom_put( "MNUEN"  , zmnuen(:,:,:) * tmask(:,:,:)  )  ! Iron limitation term
        IF( iom_use( "COZNUED" ) )   CALL iom_put( "COZNUEN", zcoznuen(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "MNLIMD"  ) )   CALL iom_put( "MNLIMN" , mnlimn(:,:,:) * tmask(:,:,:)  )  !
        IF( iom_use( "COZNLIMD") )   CALL iom_put( "COZNLIMN", coznlimn(:,:,:) * tmask(:,:,:) )  !
        IF( iom_use( "COLIMD"  ) )   CALL iom_put( "COLIMN" , colimn(:,:,:) *tmask(:,:,:)  )  !
        IF( iom_use( "PCO2"  ) )   CALL iom_put( "PCO2" , pco2s(:,:) *tmask(:,:,1)  )  !
        IF( iom_use( "CO2AQ"  ) )   CALL iom_put( "CO2AQ" , co2aq(:,:,:) *tmask(:,:,:)  )  !
        IF( iom_use( "QFePSN"  ) )   CALL iom_put( "QFePSN" , qfepsn(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "QFeNO3N" ) )   CALL iom_put( "QFeNO3N", qfeno3n(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "QFeRESN" ) )   CALL iom_put( "QFeRESN", qferesn(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "QFePSD"  ) )   CALL iom_put( "QFePSD" , qfepsd(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "QFeNO3D" ) )   CALL iom_put( "QFeNO3D", qfeno3d(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "QFeRESD" ) )   CALL iom_put( "QFeRESD", qferesd(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LNN"    ) )   CALL iom_put( "LNN"   , xlimnnit(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LDN"    ) )   CALL iom_put( "LDN"   , xlimdnit(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LNP"    ) )   CALL iom_put( "LNP"   , xnanopo4(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LDP"    ) )   CALL iom_put( "LDP"   , xdiatpo4(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_lim')
      !
   END SUBROUTINE p4z_lim


   SUBROUTINE p4z_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of nutrient limitation parameters
      !!
      !! ** Method  :   Read the nampislim namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampislim
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zlim/ concnno3, concdno3, concnnh4, concdnh4, concnfer, concdfer, concbfe,   &
         &                concbno3, concbnh4, xsizedia, xsizephy, xsizern, xsizerd,          & 
         &                xksi1, xksi2, xkdoc, qnfelim, qdfelim, caco3r, oxymin,kscun, kscud,  &
         &                ksznn, ksznd, kscon, kscod, ksmnn, ksmnd, & 
         &                ln_colim, ln_colim_int, qmnmin, mnchln, mnchld, qcoznmin, qcozncan, qcozncad, &
         &                kco2, qcoznapn, qcoznapd, kpo4
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_lim_init : initialization of nutrient limitations'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampislim in reference namelist : Pisces nutrient limitation parameters
      READ  ( numnatp_ref, namp4zlim, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zlim in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampislim in configuration namelist : Pisces nutrient limitation parameters 
      READ  ( numnatp_cfg, namp4zlim, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zlim in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp4zlim )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zlim'
         WRITE(numout,*) '      mean rainratio                           caco3r    = ', caco3r
         WRITE(numout,*) '      NO3 half saturation of nanophyto         concnno3  = ', concnno3
         WRITE(numout,*) '      NO3 half saturation of diatoms           concdno3  = ', concdno3
         WRITE(numout,*) '      NH4 half saturation for phyto            concnnh4  = ', concnnh4
         WRITE(numout,*) '      NH4 half saturation for diatoms          concdnh4  = ', concdnh4
         WRITE(numout,*) '      half saturation constant for Si uptake   xksi1     = ', xksi1
         WRITE(numout,*) '      half saturation constant for Si/C        xksi2     = ', xksi2
         WRITE(numout,*) '      half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(numout,*) '      Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(numout,*) '      Iron half saturation for diatoms         concdfer  = ', concdfer
         WRITE(numout,*) '      size ratio for nanophytoplankton         xsizern   = ', xsizern
         WRITE(numout,*) '      size ratio for diatoms                   xsizerd   = ', xsizerd
         WRITE(numout,*) '      NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(numout,*) '      NH4 half saturation for bacteria         concbnh4  = ', concbnh4
         WRITE(numout,*) '      Minimum size criteria for diatoms        xsizedia  = ', xsizedia
         WRITE(numout,*) '      Minimum size criteria for nanophyto      xsizephy  = ', xsizephy
         WRITE(numout,*) '      Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(numout,*) '      halk saturation constant for anoxia       oxymin   =' , oxymin
         WRITE(numout,*) '      optimal Fe quota for nano.               qnfelim   = ', qnfelim
         WRITE(numout,*) '      Optimal Fe quota for diatoms             qdfelim   = ', qdfelim
         WRITE(numout,*) '      Cu half saturation of nanophyt             kscun   = ', kscun
         WRITE(numout,*) '      Cu half saturation of diatoms              kscud   = ', kscud
         WRITE(numout,*) '      Zn half saturation of nanophyt             ksznn   = ', ksznn
         WRITE(numout,*) '      Zn half saturation of diatoms              ksznd   = ', ksznd
         WRITE(numout,*) '    Co half saturation of nanophyt               kscon   = ', kscon
         WRITE(numout,*) '    Co half saturation of diatoms                kscod   = ', kscod
         WRITE(numout,*) '    Mn half saturation of nanophyt               ksmnn   = ', ksmnn
         WRITE(numout,*) '    Mn half saturation of diatoms                ksmnd   = ', ksmnd
         WRITE(numout,*) ' qmnmin=', qmnmin
         WRITE(numout,*) ' mnchln=', mnchln
         WRITE(numout,*) ' mnchld=', mnchld
         WRITE(numout,*) ' qcoznmin=',qcoznmin
         WRITE(numout,*) ' qcozncan=',qcozncan
         WRITE(numout,*) ' qcozncad=',qcozncad
         WRITE(numout,*) ' kco2=',kco2
         WRITE(numout,*) ' qcoznapn=',qcoznapn
         WRITE(numout,*) ' qcoznapd=',qcoznapd
         WRITE(numout,*) ' kpo4=',kpo4
      ENDIF
      !
      nitrfac (:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_lim_init


   INTEGER FUNCTION p4z_lim_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_lim_alloc  ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_stop
      !!----------------------------------------------------------------------

      !*  Biological arrays for phytoplankton growth
      ALLOCATE( xnanono3(jpi,jpj,jpk), xdiatno3(jpi,jpj,jpk),       &
         &      xnanonh4(jpi,jpj,jpk), xdiatnh4(jpi,jpj,jpk),       &
         &      xnanopo4(jpi,jpj,jpk), xdiatpo4(jpi,jpj,jpk),       &
         &      xlimphy (jpi,jpj,jpk), xlimdia (jpi,jpj,jpk),       &
         &      xlimnfe (jpi,jpj,jpk), xlimdfe (jpi,jpj,jpk),       &
         &      xlimbac (jpi,jpj,jpk), xlimbacl(jpi,jpj,jpk),       &
         &      concnfe (jpi,jpj,jpk), concdfe (jpi,jpj,jpk),       &
         &      xlimsi  (jpi,jpj,jpk), STAT=p4z_lim_alloc )
      IF (ln_copper) THEN !Cam
      ALLOCATE( concdcu(jpi,jpj,jpk), concncu(jpi,jpj,jpk), STAT=p4z_lim_alloc )
      ENDIF
      IF (ln_zinc) THEN !Cam
      ALLOCATE( concdzn(jpi,jpj,jpk), concnzn(jpi,jpj,jpk), STAT=p4z_lim_alloc )
      ENDIF
      IF (ln_cobalt) THEN
      ALLOCATE( concdco(jpi,jpj,jpk), concnco(jpi,jpj,jpk),  STAT=p4z_lim_alloc )
      ENDIF
      IF (ln_manganese) THEN
      ALLOCATE( concdmn(jpi,jpj,jpk), concnmn(jpi,jpj,jpk),         &
         &      zmnued(jpi,jpj,jpk),  zcoznued(jpi,jpj,jpk),        &
         &      zmnuen(jpi,jpj,jpk),  zcoznuen(jpi,jpj,jpk),        &
         &      coznlimn(jpi,jpj,jpk), mnlimn(jpi,jpj,jpk), colimn(jpi,jpj,jpk),  &
         &      coznlimd(jpi,jpj,jpk), mnlimd(jpi,jpj,jpk), colimd(jpi,jpj,jpk), STAT=p4z_lim_alloc )
      ENDIF

      !
      IF( p4z_lim_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_lim_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_lim_alloc

   !!======================================================================
END MODULE p4zlim
