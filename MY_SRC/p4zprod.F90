MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!----------------------------------------------------------------------
   !!   p4z_prod       : Compute the growth Rate of the two phytoplanktons groups
   !!   p4z_prod_init  : Initialization of the parameters for growth
   !!   p4z_prod_alloc : Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zlim          ! Co-limitations of differents nutrients
   USE prtctl_trc      ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_prod_alloc

   REAL(wp), PUBLIC ::   pislopen     !:
   REAL(wp), PUBLIC ::   pisloped     !:
   REAL(wp), PUBLIC ::   xadap        !:
   REAL(wp), PUBLIC ::   excretn      !:
   REAL(wp), PUBLIC ::   excretd      !:
   REAL(wp), PUBLIC ::   bresp        !:
   REAL(wp), PUBLIC ::   chlcnm       !:
   REAL(wp), PUBLIC ::   chlcdm       !:
   REAL(wp), PUBLIC ::   chlcmin      !:
   REAL(wp), PUBLIC ::   fecnm        !:
   REAL(wp), PUBLIC ::   fecdm        !:
   REAL(wp), PUBLIC ::   grosip       !:
   REAL(wp), PUBLIC ::  cupmn           !: Cam Copper model
   REAL(wp), PUBLIC ::  cupmd
   REAL(wp), PUBLIC ::  znpmn
   REAL(wp), PUBLIC ::  znpmd
   LOGICAL , PUBLIC ::  ln_fezn
   LOGICAL , PUBLIC ::  ln_cozn
   REAL(wp), PUBLIC ::  kcozn
   REAL(wp), PUBLIC ::  copmn
   REAL(wp), PUBLIC ::  copmd
   REAL(wp), PUBLIC ::  mnpmn
   REAL(wp), PUBLIC ::  mnpmd
   LOGICAL , PUBLIC ::  ln_comnzn
   LOGICAL , PUBLIC ::  ln_comnzn_simple
   REAL(wp), PUBLIC ::  mntmaxd
   REAL(wp), PUBLIC ::  kmnd
   REAL(wp), PUBLIC ::  kznd
   REAL(wp), PUBLIC ::  kcoznd
   REAL(wp), PUBLIC ::  cozntmaxd
   REAL(wp), PUBLIC ::  kb12d
   REAL(wp), PUBLIC ::  qmnmaxd
   REAL(wp), PUBLIC ::  qznmaxd
   LOGICAL , PUBLIC ::  ln_mnzn_int
   LOGICAL , PUBLIC ::  ln_cozn_int
   LOGICAL , PUBLIC ::  ln_mnqzn_int
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotad   !: proxy of N quota in diatomee
   
   REAL(wp) ::   r1_rday    ! 1 / rday
   REAL(wp) ::   texcretn   ! 1 - excretn 
   REAL(wp) ::   texcretd   ! 1 - excretd        

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zprod.F90 10401 2018-12-17 15:57:34Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zsilfac, znanotot, zdiattot, zconctemp, zconctemp2
      REAL(wp) ::   zratio, zmax, zsilim, ztn, zadap, zlim, zsilfac2, zsiborn
      REAL(wp) ::   zprod, zproreg, zproreg2, zprochln, zprochld
      REAL(wp) ::   zmaxday, zdocprod, zpislopen, zpisloped
      REAL(wp) ::   zmxltst, zmxlday
      REAL(wp) ::   zrum, zcodel, zargu, zval, zfeup, chlcnm_n, chlcdm_n
      REAL(wp) ::   zfact
      CHARACTER (len=25) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), DIMENSION(jpi,jpj    ) :: zstrn, zmixnano, zmixdiat
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprmaxn,zprmaxd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpislopeadn, zpislopeadd, zysopt  
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprdia, zprbio, zprdch, zprnch   
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprorcan, zprorcad, zprofed, zprofen
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpronewn, zpronewd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmxl_fac, zmxl_chl
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpligprod1, zpligprod2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprocun, zprocud
      REAL(wp) :: dcu
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zproznn, zproznd, zcoznt
      REAL(wp) :: znf1, dzn, zn_fun, dco, dmn
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprocon, zprocod
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpromnn, zpromnd, zmnt
      REAL(wp) :: zqmn, zqzn, zqcozn, zmax2, vmzn, vmmn, vmco, vm
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: mnznd, mnznn, znmnd, znmnn, coznd, mncod
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmax_mnznn, zmax_mnznd
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_prod')
      !
      !  Allocate temporary workspace
      !
      zprorcan(:,:,:) = 0._wp ; zprorcad(:,:,:) = 0._wp ; zprofed (:,:,:) = 0._wp
      zprofen (:,:,:) = 0._wp ; zysopt  (:,:,:) = 0._wp
      zpronewn(:,:,:) = 0._wp ; zpronewd(:,:,:) = 0._wp ; zprdia  (:,:,:) = 0._wp
      zprbio  (:,:,:) = 0._wp ; zprdch  (:,:,:) = 0._wp ; zprnch  (:,:,:) = 0._wp 
      zmxl_fac(:,:,:) = 0._wp ; zmxl_chl(:,:,:) = 0._wp 
      IF( ln_ligand ) THEN
           zpligprod1(:,:,:)  = 0._wp     ;      zpligprod2(:,:,:) = 0._wp
      ENDIF
     IF( ln_copper ) THEN
         !ALLOCATE(zprocun(jpi,jpj,jpk), zprocud(jpi,jpj,jpk))
         zprocun(:,:,:) = 0._wp ; zprocud(:,:,:) = 0._wp
     ENDIF
     IF( ln_cobalt ) THEN
      zprocon(:,:,:) = 0._wp ; zprocod(:,:,:) = 0._wp
     ENDIF
     IF( ln_manganese ) THEN
      zpromnn(:,:,:) = 0._wp ; zpromnd(:,:,:) = 0._wp
     ENDIF

    IF (ln_zinc) THEN
         zproznn(:,:,:) = 0._wp ; zproznd(:,:,:) = 0._wp ; mnznd(:,:,:) = 0._wp ; mnznn(:,:,:) = 0._wp
         znmnd(:,:,:) = 0._wp ; znmnn(:,:,:) = 0._wp ; coznd(:,:,:) = 0._wp ; mncod(:,:,:) = 0._wp
         zmax_mnznn(:,:,:) = 1._wp ; zmax_mnznd(:,:,:) = 1._wp
    ENDIF
      zcoznt(:,:,:) = 0. ; zmnt(:,:,:) = 0.
      ! Computation of the optimal production
      zprmaxn(:,:,:) = 0.8_wp * r1_rday * tgfunc(:,:,:)
      zprmaxd(:,:,:) = zprmaxn(:,:,:)

      ! compute the day length depending on latitude and the day
      zrum = REAL( nday_year - 80, wp ) / REAL( nyear_len(1), wp )
      zcodel = ASIN(  SIN( zrum * rpi * 2._wp ) * SIN( rad * 23.5_wp )  )

      ! day length in hours
      zstrn(:,:) = 0.
      DO jj = 1, jpj
         DO ji = 1, jpi
            zargu = TAN( zcodel ) * TAN( gphit(ji,jj) * rad )
            zargu = MAX( -1., MIN(  1., zargu ) )
            zstrn(ji,jj) = MAX( 0.0, 24. - 2. * ACOS( zargu ) / rad / 15. )
         END DO
      END DO

      ! Impact of the day duration and light intermittency on phytoplankton growth
      DO jk = 1, jpkm1
         DO jj = 1 ,jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  zval = MAX( 1., zstrn(ji,jj) )
                  IF( gdept_n(ji,jj,jk) <= hmld(ji,jj) ) THEN
                     zval = zval * MIN(1., heup_01(ji,jj) / ( hmld(ji,jj) + rtrn ))
                  ENDIF
                  zmxl_chl(ji,jj,jk) = zval / 24.
                  zmxl_fac(ji,jj,jk) = 1.5 * zval / ( 12. + zval )
               ENDIF
            END DO
         END DO
      END DO

      zprbio(:,:,:) = zprmaxn(:,:,:) * zmxl_fac(:,:,:)
      zprdia(:,:,:) = zprmaxd(:,:,:) * zmxl_fac(:,:,:)

      ! Maximum light intensity
      WHERE( zstrn(:,:) < 1.e0 ) zstrn(:,:) = 24.

      ! Computation of the P-I slope for nanos and diatoms
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  ztn         = MAX( 0., tsn(ji,jj,jk,jp_tem) - 15. )
                  zadap       = xadap * ztn / ( 2.+ ztn )
                  zconctemp   = MAX( 0.e0 , trb(ji,jj,jk,jpdia) - xsizedia )
                  zconctemp2  = trb(ji,jj,jk,jpdia) - zconctemp
                  !
                  zpislopeadn(ji,jj,jk) = pislopen * ( 1.+ zadap  * EXP( -0.25 * enano(ji,jj,jk) ) )  &
                  &                   * trb(ji,jj,jk,jpnch) /( trb(ji,jj,jk,jpphy) * 12. + rtrn)
                  !
                  zpislopeadd(ji,jj,jk) = (pislopen * zconctemp2 + pisloped * zconctemp) / ( trb(ji,jj,jk,jpdia) + rtrn )   &
                  &                   * trb(ji,jj,jk,jpdch) /( trb(ji,jj,jk,jpdia) * 12. + rtrn)
               ENDIF
            END DO
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                   ! Computation of production function for Carbon
                   !  ---------------------------------------------
                   zpislopen = zpislopeadn(ji,jj,jk) / ( ( r1_rday + bresp * r1_rday ) &
                   &            * zmxl_fac(ji,jj,jk) * rday + rtrn)
                   zpisloped = zpislopeadd(ji,jj,jk) / ( ( r1_rday + bresp * r1_rday ) &
                   &            * zmxl_fac(ji,jj,jk) * rday + rtrn)
                   zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1.- EXP( -zpislopen * enano(ji,jj,jk) )  )
                   zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediat(ji,jj,jk) )  )
                   !  Computation of production function for Chlorophyll
                   !--------------------------------------------------
                   zpislopen = zpislopeadn(ji,jj,jk) / ( zprmaxn(ji,jj,jk) * zmxl_chl(ji,jj,jk) * rday + rtrn )
                   zpisloped = zpislopeadd(ji,jj,jk) / ( zprmaxd(ji,jj,jk) * zmxl_chl(ji,jj,jk) * rday + rtrn )
                   zprnch(ji,jj,jk) = zprmaxn(ji,jj,jk) * ( 1.- EXP( -zpislopen * enanom(ji,jj,jk) ) )
                   zprdch(ji,jj,jk) = zprmaxd(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediatm(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      !  Computation of a proxy of the N/C ratio
      !  ---------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                zval = MIN( xnanopo4(ji,jj,jk), ( xnanonh4(ji,jj,jk) + xnanono3(ji,jj,jk) ) )   &
                &      * zprmaxn(ji,jj,jk) / ( zprbio(ji,jj,jk) + rtrn )
                quotan(ji,jj,jk) = MIN( 1., 0.2 + 0.8 * zval )
                zval = MIN( xdiatpo4(ji,jj,jk), ( xdiatnh4(ji,jj,jk) + xdiatno3(ji,jj,jk) ) )   &
                &      * zprmaxd(ji,jj,jk) / ( zprdia(ji,jj,jk) + rtrn )
                quotad(ji,jj,jk) = MIN( 1., 0.2 + 0.8 * zval )
            END DO
         END DO
      END DO


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

                IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                   !    Si/C of diatoms
                   !    ------------------------
                   !    Si/C increases with iron stress and silicate availability
                   !    Si/C is arbitrariliy increased for very high Si concentrations
                   !    to mimic the very high ratios observed in the Southern Ocean (silpot2)
                  zlim  = trb(ji,jj,jk,jpsil) / ( trb(ji,jj,jk,jpsil) + xksi1 )
                  zsilim = MIN( zprdia(ji,jj,jk) / ( zprmaxd(ji,jj,jk) + rtrn ), xlimsi(ji,jj,jk) )
                  zsilfac = 4.4 * EXP( -4.23 * zsilim ) * MAX( 0.e0, MIN( 1., 2.2 * ( zlim - 0.5 ) )  ) + 1.e0
                  zsiborn = trb(ji,jj,jk,jpsil) * trb(ji,jj,jk,jpsil) * trb(ji,jj,jk,jpsil)
                  IF (gphit(ji,jj) < -30 ) THEN
                    zsilfac2 = 1. + 2. * zsiborn / ( zsiborn + xksi2**3 )
                  ELSE
                    zsilfac2 = 1. +      zsiborn / ( zsiborn + xksi2**3 )
                  ENDIF
                  zysopt(ji,jj,jk) = grosip * zlim * zsilfac * zsilfac2
              ENDIF
            END DO
         END DO
      END DO

      !  Mixed-layer effect on production 
      !  Sea-ice effect on production

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
            END DO
         END DO
      END DO

      ! Computation of the various production terms 
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto. (C)
                  zprorcan(ji,jj,jk) = zprbio(ji,jj,jk)  * xlimphy(ji,jj,jk) * trb(ji,jj,jk,jpphy) * rfact2
                  zpronewn(ji,jj,jk)  = zprorcan(ji,jj,jk)* xnanono3(ji,jj,jk) / ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) + rtrn )
                  !
!                  zprofen(ji,jj,jk) = fecnm * zprmaxn(ji,jj,jk) * ( 1.0 - fr_i(ji,jj) )  &
                     vm = qminfen(ji,jj,jk) + ( (fecnm) - qminfen(ji,jj,jk) ) * ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) )
                     vm = vm + ( MIN( 40.E-6, fecnm-vm) * ( biron(ji,jj,jk)**2 / ( biron(ji,jj,jk)**2 + 2.E-9**2 ) ) )
!                  zratio = trb(ji,jj,jk,jpnfe) / ( trb(ji,jj,jk,jpphy) * fecnm + rtrn )
                  zratio = trb(ji,jj,jk,jpnfe) / ( trb(ji,jj,jk,jpphy) * vm + rtrn )
                  zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )

                  zprofen(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * ( 1.0 - fr_i(ji,jj) )  &
                  &             * ( 4. - 4.5 * xlimnfe(ji,jj,jk) / ( xlimnfe(ji,jj,jk) + 0.5 ) )    &
                  &             * biron(ji,jj,jk) / ( biron(ji,jj,jk) + concnfe(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
                  !  production terms for diatoms (C)
                  zprorcad(ji,jj,jk) = zprdia(ji,jj,jk) * xlimdia(ji,jj,jk) * trb(ji,jj,jk,jpdia) * rfact2
                  zpronewd(ji,jj,jk) = zprorcad(ji,jj,jk) * xdiatno3(ji,jj,jk) / ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) + rtrn )
                  !
                  zratio = trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) * fecdm + rtrn )
                  zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) ) 
!                  zprofed(ji,jj,jk) = fecdm * zprmaxd(ji,jj,jk) * ( 1.0 - fr_i(ji,jj) )  &
                     vm = qminfed(ji,jj,jk) + ( (fecdm) - qminfed(ji,jj,jk) ) * ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) )
                     vm = vm + ( MIN( 40.E-6, fecdm-vm) * ( biron(ji,jj,jk)**2 / ( biron(ji,jj,jk)**2 + 2.E-9**2 ) ) )
!                  zratio = trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) * fecdm + rtrn )
                  zratio = trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) * vm + rtrn )
                  zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )

                  zprofed(ji,jj,jk) = ( vm )* zprmaxd(ji,jj,jk) * ( 1.0 - fr_i(ji,jj) )  &
                  &             * ( 4. - 4.5 * xlimdfe(ji,jj,jk) / ( xlimdfe(ji,jj,jk) + 0.5 ) )    &
                  &             * biron(ji,jj,jk) / ( biron(ji,jj,jk) + concdfe(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
                  mu_nm(ji,jj,jk) = zprbio(ji,jj,jk)
                  mu_n(ji,jj,jk) = zprbio(ji,jj,jk) * xlimphy(ji,jj,jk)
                  mu_dm(ji,jj,jk) = zprdia(ji,jj,jk)
                  mu_d(ji,jj,jk) = zprdia(ji,jj,jk) * xlimdia(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the chlorophyll production terms
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto. ( chlorophyll )
                  znanotot = enanom(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod    = rday * zprorcan(ji,jj,jk) * zprnch(ji,jj,jk) * xlimphy(ji,jj,jk)
                  zprochln = chlcmin * 12. * zprorcan (ji,jj,jk)
                  chlcnm_n   = MIN ( chlcnm, ( chlcnm / (1. - 1.14 / 43.4 *tsn(ji,jj,jk,jp_tem))) * (1. - 1.14 / 43.4 * 20.))
                  zprochln = zprochln + (chlcnm_n-chlcmin) * 12. * zprod / &
                                        & (  zpislopeadn(ji,jj,jk) * znanotot +rtrn)
                  !  production terms for diatoms ( chlorophyll )
                  zdiattot = ediatm(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod    = rday * zprorcad(ji,jj,jk) * zprdch(ji,jj,jk) * xlimdia(ji,jj,jk)
                  zprochld = chlcmin * 12. * zprorcad(ji,jj,jk)
                  chlcdm_n   = MIN ( chlcdm, ( chlcdm / (1. - 1.14 / 43.4 * tsn(ji,jj,jk,jp_tem))) * (1. - 1.14 / 43.4 * 20.))
                  zprochld = zprochld + (chlcdm_n-chlcmin) * 12. * zprod / &
                                        & ( zpislopeadd(ji,jj,jk) * zdiattot +rtrn )
                  !   Update the arrays TRA which contain the Chla sources and sinks
                  tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln * texcretn
                  tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) + zprochld * texcretd
               ENDIF
            END DO
         END DO
      END DO
      IF ( ln_copper ) THEN !Cam
      ! Computation of the copper terms
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                     dcu = max( bcopper(ji,jj,jk) , 0. )
! Nanos:
                     zratio = cupmn * po4r
                     zratio = trb(ji,jj,jk, jpcun) / ( trb(ji,jj,jk,jpphy) * zratio + rtrn )
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )
                     vm = ( cupmn * po4r ) * MAX( 0.2, (xnanono3(ji,jj,jk) +xnanonh4(ji,jj,jk) ) )
                     zprocun(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dcu / ( dcu + concncu(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
! Diatoms:
                     zratio = cupmd * po4r
                     zratio = trb(ji,jj,jk, jpcud) / ( trb(ji,jj,jk,jpdia) * zratio + rtrn )
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )
                     vm = ( cupmn * po4r ) * MAX( 0.2 , (xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) ) )
                     zprocud(ji,jj,jk) = ( vm ) *zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk))  &
                  &             * dcu / ( dcu + concdcu(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
               ENDIF
            END DO
         END DO
      END DO
      ENDIF ! ln_copper no interaction with fe

      IF ( ln_comnzn ) THEN ! interactive scheme

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
! set bioavailable nutrient pools:

!                     dmn = max( zmnfree(ji,jj,jk), 0. )
                     dmn = max( trb(ji,jj,jk,jpdmn), 0. )
                     dzn = max( bzinc(ji,jj,jk) , 0. )
!                     dzn = trb(ji,jj,jk,jpdzn)
                     dco = max( zcofree(ji,jj,jk) * 1.E-3,  0. )
!! diatoms
! quotas
                 zqzn = trb(ji,jj,jk,jpznd) / ( trb(ji,jj,jk,jpdia) + rtrn )
                 zqmn = trb(ji,jj,jk,jpmnd) / ( trb(ji,jj,jk,jpdia) + rtrn )
                 zqcozn = ( trb(ji,jj,jk,jpznd) + 1.E-3 * trb(ji,jj,jk,jpcod) ) / ( trb(ji,jj,jk,jpdia) + rtrn )
! number of Mn transporters:
!                zmnt(ji,jj,jk) = (1 - ( zqmn**2 / ( zqmn**2 + qmnmaxd**2 ) ) ) & 
!                &              * (1 - ( zqzn**2 / ( zqzn**2 + qznmaxd**2 ) ) ) &
                zmnt(ji,jj,jk) = (1 - ( zqzn**2 / ( zqzn**2 + qznmaxd**2 ) ) ) &
                &              * (1 - ( zqmn**4 / ( zqmn**4 + ( 2 * mu_d(ji,jj,jk) &
                &              / zmnued(ji,jj,jk) )**4 + rtrn ) ) )
                zmnt(ji,jj,jk) = MAX(0., zmnt(ji,jj,jk) )
! number of high affinity Co and Zn transporters:
                zcoznt(ji,jj,jk) = ( 1 - ( zqcozn**4 / ( zqcozn**4 + ( 2 * mu_d(ji,jj,jk) &
                &                  / zcoznued(ji,jj,jk) )**4 + rtrn ) ) )
                zcoznt(ji,jj,jk) = MAX(0., zcoznt(ji,jj,jk) )
!!! uptake terms
! Mn uptake rate:
                     zratio = trb(ji,jj,jk,jpmnd) / ( trb(ji,jj,jk,jpdia) + rtrn )
                     zratio = zratio / ( mnpmd * po4r) ! 
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )

                 zpromnd(ji,jj,jk) =  ( mntmaxd * zprmaxd(ji,jj,jk) * zmnt(ji,jj,jk) ) &
                 &  * ( (kmnd*dmn) / ( 1 + (kmnd*dmn)  + (kznd*dzn) ) ) & 
                 &  * (1.0 - fr_i(ji,jj)) * trb(ji,jj,jk,jpdia) * rfact2
! Zn uptake
                     zratio = trb(ji,jj,jk, jpznd) / ( trb(ji,jj,jk,jpdia) + rtrn )
                     zratio = zratio / ( znpmd * po4r) !
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )

                 zproznd(ji,jj,jk) = ( (mntmaxd * zprmaxd(ji,jj,jk) * zmnt(ji,jj,jk))   &
                 & * ( (kznd*dzn) / (1 + (kmnd*dmn)  + (kznd*dzn) ) ) &
                 & + ( ( cozntmaxd * zprmaxd(ji,jj,jk) * zcoznt(ji,jj,jk) ) &
                 & * ( (kcozn * dzn) / ( 1 + ( kcoznd * ( dco + dzn ) ) ) ) ) )  & 
                 &   * (1.0 - fr_i(ji,jj)) * trb(ji,jj,jk,jpdia) * rfact2
! Co uptake
                     zratio = trb(ji,jj,jk, jpcod) / ( trb(ji,jj,jk,jpdia) * 1000. + rtrn )
                     zratio = zratio / ( copmd * po4r) ! molCo / mol C : mol Co / mol C
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )! no units
                 zprocod(ji,jj,jk) = ( ( cozntmaxd * zprmaxd(ji,jj,jk) * zcoznt(ji,jj,jk) ) &
                 & * ( (kcoznd * dco) / (1 + ( kcoznd * ( dco + dzn ) ) ) )  &
!                & + ( zprmaxd(ji,jj,jk) * kb12d * dco ) ) & 
                & ) &
                 & * (1.0 - fr_i(ji,jj)) * trb(ji,jj,jk,jpdia) * rfact2
!!!********************************************************************
! nanos, copy from stnd:
!                     dmn = max(zmnfree(ji,jj,jk), 0. )
                     dmn = max( trb(ji,jj,jk,jpdmn), 0. )
                     zratio = mnpmn * po4r
                     zratio = trb(ji,jj,jk, jpmnn) / ( trb(ji,jj,jk,jpphy) * zratio + rtrn )
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )
                     zpromnn(ji,jj,jk) = ( mnpmn * po4r ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dmn / ( dmn + concnmn(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
                     dzn = max( bzinc(ji,jj,jk) , 0. )
                     zratio = znpmn * po4r
                     zratio = trb(ji,jj,jk, jpznn) / ( trb(ji,jj,jk,jpphy) * zratio + rtrn )
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )
      IF( ln_fezn ) THEN
                     zproznn(ji,jj,jk) = ( znpmn * po4r ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dzn / ( dzn + concnzn(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * xlimnfe(ji,jj,jk) / ( xlimnfe(ji,jj,jk) + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ELSE
                     zproznn(ji,jj,jk) = ( znpmn * po4r ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dzn / ( dzn + concnzn(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ENDIF
                     dco = max( trb(ji,jj,jk, jpdco)-1.E-9, 0. )
                     zn_fun = dzn / ( dzn + kcozn ) !
                     zratio = trb(ji,jj,jk, jpcon) / ( trb(ji,jj,jk,jpphy) * 1000. + rtrn )
                     zratio = zratio / ( copmn * po4r) ! molCo / mol C : mol Co / mol C
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )! no units
                     zprocon(ji,jj,jk) = ( copmn * po4r) * zprmaxn(ji,jj,jk)  &
                  &             * dco / ( dco + concnco(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
               ENDIF
            END DO
         END DO
      END DO


      ELSE ! no ln_comnzn, allows 'stnd' uptake scheme or a simple coznmn interaction scheme

      IF ( ln_zinc ) THEN
      ! Computation of the zinc terms
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                     dzn = max( bzinc(ji,jj,jk) , 0. )
! Nanos:
!                     vm = ( znpmn * po4r ) * MAX( 0.2, (xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) ) )
                     vm = qminczn(ji,jj,jk) + ( (znpmn * po4r) - qminczn(ji,jj,jk) ) * ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) )
                     zratio = znpmn * po4r
                     zratio = vm
                     zratio = trb(ji,jj,jk, jpznn) / ( trb(ji,jj,jk,jpphy) * zratio + rtrn )
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )
      IF ( ln_comnzn_simple) THEN
                 zqcozn = ( trb(ji,jj,jk,jpznn) + 1.E-3 * trb(ji,jj,jk,jpcon) ) / ( trb(ji,jj,jk,jpphy) + rtrn )
                 zmax2 = MIN( 1., ( zcoznuen(ji,jj,jk)  * zqcozn ) / (mu_nm(ji,jj,jk) + rtrn) )
      IF (ln_mnqzn_int) THEN
                 zmax_mnznn(ji,jj,jk) = MAX( 0.1, ( 1. - zratio ) / ABS( 1.05 - zratio ) )
      ENDIF
      IF( ln_fezn ) THEN
                     zproznn(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dzn / ( dzn + concnzn(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * xlimnfe(ji,jj,jk) / ( xlimnfe(ji,jj,jk) + 0.5 ) )    &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ELSE
                     zproznn(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dzn / ( dzn + concnzn(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ENDIF ! ln_fezn
      ELSE ! ln_comnzn_simple
      IF( ln_fezn ) THEN
                     zproznn(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dzn / ( dzn + concnzn(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * xlimnfe(ji,jj,jk) / ( xlimnfe(ji,jj,jk) + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ELSE
                     zproznn(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dzn / ( dzn + concnzn(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ENDIF ! ln_fezn
      ENDIF ! ln_comnzn_simple   
! Diatoms:
!                     vm = ( znpmd * po4r ) * MAX ( 0.2, (xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) ) )
                     vm = qminczd(ji,jj,jk) + ( (znpmd * po4r) - qminczd(ji,jj,jk) ) * ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) )
                     zratio = znpmd * po4r
                     zratio = vm
                     zratio = trb(ji,jj,jk, jpznd) / ( trb(ji,jj,jk,jpdia) * zratio + rtrn )
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )
      IF ( ln_comnzn_simple) THEN
                 zqcozn = ( trb(ji,jj,jk,jpznd) + 1.E-3 * trb(ji,jj,jk,jpcod) ) / ( trb(ji,jj,jk,jpdia) + rtrn ) 
                 zmax2 = MIN( 1., ( zcoznued(ji,jj,jk)  * zqcozn ) / (mu_dm(ji,jj,jk) + rtrn) )
      IF (ln_mnqzn_int) THEN
                zmax_mnznd(ji,jj,jk) = MAX( 0.1, ( 1. - zratio ) / ABS( 1.05 - zratio ) )
      ENDIF
      IF( ln_fezn ) THEN
                     zproznd(ji,jj,jk) = ( vm ) * zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk))  &
                  &             * dzn / ( dzn + concdzn(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * xlimdfe(ji,jj,jk) / ( xlimdfe(ji,jj,jk) + 0.5 ) )    &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
      ELSE
                     zproznd(ji,jj,jk) = ( vm ) * zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk))  &
                  &             * dzn / ( dzn + concdzn(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
      ENDIF ! ln_fezn
      ELSE  ! ln_comnzn_simple
      IF( ln_fezn ) THEN
                     zproznd(ji,jj,jk) = ( vm ) * zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk))  &
                  &             * dzn / ( dzn + concdzn(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * xlimdfe(ji,jj,jk) / ( xlimdfe(ji,jj,jk) + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
      ELSE
                     zproznd(ji,jj,jk) = ( vm ) * zprmaxd(ji,jj,jk) *(1.0 - fr_i(ji,jk))  &
                  &             * dzn / ( dzn + concdzn(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
      ENDIF ! ln_fezn
      ENDIF ! ln_comnzn_simple
               ENDIF
            END DO
         END DO
      END DO
      ENDIF ! ln_zinc

      IF ( ln_manganese ) THEN !Cam
      ! Computation of the Mn terms
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
!                     dmn = max(zmnfree(ji,jj,jk), 0. )
                     dmn = max( trb(ji,jj,jk,jpdmn), 0. )
                     dzn = max( bzinc(ji,jj,jk) , 0. )
! Nanos:
!                     vm = ( mnpmn * po4r ) * MAX( 0.2, (xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) ) )
                     vm = qminmnn(ji,jj,jk) + ( (mnpmn * po4r) - qminmnn(ji,jj,jk) ) * ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) )                     
                     zratio = mnpmn * po4r
                     zratio = vm
                     zratio = trb(ji,jj,jk, jpmnn) / ( trb(ji,jj,jk,jpphy) * zratio + rtrn )
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )
      IF ( ln_comnzn_simple) THEN
                 zqmn = trb(ji,jj,jk,jpmnn) / ( trb(ji,jj,jk,jpphy) + rtrn )
                 zmax2 = MIN( 1., mnlimn(ji,jj,jk) )
                     kmnd = 1. / concnmn(ji,jj,jk)                     
      IF (ln_mnzn_int) THEN
                     zpromnn(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * (dmn*kmnd) / (1. + (dmn*kmnd) + (kznd*dzn) ) &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * zmax_mnznn(ji,jj,jk) * trb(ji,jj,jk,jpphy) * rfact2
      ELSE
                     zpromnn(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dmn / ( dmn + concnmn(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ENDIF ! ln_mnzn_int
      ELSE
                     zpromnn(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jj))  &
                  &             * dmn / ( dmn + concnmn(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ENDIF ! ln_comnzn_simple
! Diatoms:
!                     vm = ( mnpmd * po4r ) * MAX( 0.2, (xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) ) )
                     vm = qminmnd(ji,jj,jk) + ( (mnpmd * po4r) - qminmnd(ji,jj,jk) ) * ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) )
                     zratio = mnpmd * po4r
                     zratio = vm
                     zratio = trb(ji,jj,jk, jpmnd) / ( trb(ji,jj,jk,jpdia) * zratio + rtrn )
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )
      IF ( ln_comnzn_simple) THEN
                 zqmn = trb(ji,jj,jk,jpmnd) / ( trb(ji,jj,jk,jpdia) + rtrn )
!                 zmax2 = MIN( 1., ( zmnued(ji,jj,jk)  * zqmn ) / (mu_dm(ji,jj,jk) + rtrn) )
                 zmax2 = MIN( 1., mnlimd(ji,jj,jk) )
                     kmnd = 1. / concdmn(ji,jj,jk)
      IF (ln_mnzn_int) THEN
                     zpromnd(ji,jj,jk) = ( vm ) *zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk))  &
                  &             * (dmn*kmnd) / (1. + (dmn*kmnd) + (kznd*dzn) ) &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * zmax_mnznd(ji,jj,jk) * trb(ji,jj,jk,jpdia) * rfact2
      ELSE
                     zpromnd(ji,jj,jk) = ( vm ) *zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk))  &
                  &             * dmn / ( dmn + concdmn(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
      ENDIF ! ln_mnzn_int
      ELSE
                     zpromnd(ji,jj,jk) = ( vm ) *zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk))  &
                  &             * dmn / ( dmn + concdmn(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
      ENDIF ! ln_comnzn_simple
               ENDIF
            END DO
         END DO
      END DO
      ENDIF ! ln_manganese

      IF ( ln_cobalt ) THEN
      ! Computation of the cobalt terms
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
               IF ( ln_zinc ) THEN
                IF ( ln_cozn ) THEN
                 dzn = max( trb(ji,jj,jk,jpdzn) , 0. )
                ELSE
                 dzn = 0.065 * ( trb(ji,jj,jk,jpsil) * 1.E6 ) + 0.183 ! nmol / L,Lohan unpub
                 dzn = max( dzn * 1.E-9 , 0. )  ! mol/L
                ENDIF ! ln_cozn
               ELSE
                dzn = 0.065 * ( trb(ji,jj,jk,jpsil) * 1.E6 ) + 0.183 ! nmol / L,Lohan unpub
                dzn = max( dzn * 1.E-9 , 0. )  ! mol/L
               ENDIF ! ln_zinc
! Cobalt precision is *1000 but here is is computed in 'normal' mol/l and
! adjusted by F(1000) below when adjusting SMS
! nanos:
                     dco = max( trb(ji,jj,jk, jpdco)-1.E-9, 0. )
                     zn_fun = dzn / ( dzn + kcozn ) !
                     zratio = trb(ji,jj,jk, jpcon) / ( trb(ji,jj,jk,jpphy) * 1000. + rtrn )
                     zratio = zratio / ( copmn * po4r) ! molCo / mol C : mol Co / mol C
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )! no units
                     vm = ( copmn * po4r ) * MAX( 0.2, (xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) ) )
      IF ( ln_comnzn_simple) THEN
                 zqcozn = ( trb(ji,jj,jk,jpznn) + 1.E-3 * trb(ji,jj,jk,jpcon) ) / ( trb(ji,jj,jk,jpphy) + rtrn )
                 zmax2 = MIN( 1., coznlimn(ji,jj,jk) )
                     kmnd = 1. / concnco(ji,jj,jk)
                     zprocon(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jk)) &
                  &             * dco / ( dco + concnco(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ELSE
                     zprocon(ji,jj,jk) = ( vm ) * zprmaxn(ji,jj,jk) * (1.0 - fr_i(ji,jk)) &
                  &             * dco / ( dco + concnco(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
      ENDIF
! Diatoms:
                     dco = max( zcofree(ji,jj,jk)-1.E-9, 0. ) ! diatoms only acquire free dCo
                     zn_fun = dzn / ( dzn + kcozn ) !
                     zratio = trb(ji,jj,jk, jpcod) / ( trb(ji,jj,jk,jpdia) * 1000. + rtrn )
                     zratio = zratio / ( copmd * po4r)
                     zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) )
                     vm = ( copmd * po4r ) * MAX( 0.2, (xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) ) )
      IF ( ln_comnzn_simple) THEN
                     dzn = max( bzinc(ji,jj,jk) , 0. )
                 zqcozn = ( trb(ji,jj,jk,jpznd) + 1.E-3 * trb(ji,jj,jk,jpcod) ) / ( trb(ji,jj,jk,jpdia) + rtrn ) 
                 zmax2 = MIN( 1., coznlimn(ji,jj,jk) )
                     kmnd = 1. / concdco(ji,jj,jk)
      IF (ln_cozn_int) THEN
                     zprocod(ji,jj,jk) = ( vm ) * zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk)) &
                  &             * ( ( kcoznd * ( 1.E-3 * dco ) ) / & 
                  &             ( 1. + kcoznd * ( ( 1.E-3 * dco ) + dzn ) ) ) &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
      ELSE
                     zprocod(ji,jj,jk) = ( vm ) * zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk)) &
                  &             * dco / ( dco + concdco(ji,jj,jk) )  &
                  &             * ( 4. - 4.5 * zmax2 / ( zmax2 + 0.5 ) )    &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
      ENDIF
      ELSE
                     zprocod(ji,jj,jk) = ( vm ) * zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jk)) &
                  &             * dco / ( dco + concdco(ji,jj,jk) )  &
                  &             * MAX( 0.10 , ( 3.*(1-zn_fun) ) ) & ! impact of Zn on Co uptake for diatoms
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
      ENDIF 
              ENDIF
            END DO
         END DO
      END DO
      ENDIF ! ln_cobalt

      ENDIF ! for ln_comnzn logical

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                 zproreg  = zprorcan(ji,jj,jk) - zpronewn(ji,jj,jk)
                 zproreg2 = zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk)
                 zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)
                 tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk)
                 tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronewn(ji,jj,jk) - zpronewd(ji,jj,jk)
                 tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproreg - zproreg2
                 tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorcan(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) + zprofen(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) + zprofed(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zprorcad(ji,jj,jk) * zysopt(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zdocprod
                 tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + o2ut * ( zproreg + zproreg2) &
                 &                   + ( o2ut + o2nit ) * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) )
                 !
                 zfeup = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk)
                 tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zfeup
                 tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) - texcretd * zprorcad(ji,jj,jk) * zysopt(ji,jj,jk)
                 tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk)
                 tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) ) &
                 &                                         - rno3 * ( zproreg + zproreg2 )
              ENDIF
           END DO
        END DO
     END DO
     !
     IF( ln_ligand ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji =1 ,jpi
                 IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                    zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)
                    zfeup    = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk)
                    tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + zdocprod * ldocp    &
                    &       - zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) + 2.E3 * (1.0 - plig(ji,jj,jk) ) ) 
                    zpligprod1(ji,jj,jk) = zdocprod * ldocp
                    zpligprod2(ji,jj,jk) = zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) + 2.E3 * (1.0 - plig(ji,jj,jk) ) )
                 ENDIF
              END DO
           END DO
        END DO
     ENDIF
     IF( ln_copper ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji =1 ,jpi
                 IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
              tra(ji,jj,jk,jpdcu) = tra(ji,jj,jk,jpdcu) - ( ( zprocun(ji,jj,jk) * texcretn ) &
                 &                                      + ( zprocud(ji,jj,jk) * texcretd ) )
              tra(ji,jj,jk,jpcun) = tra(ji,jj,jk,jpcun) + zprocun(ji,jj,jk) *texcretn
              tra(ji,jj,jk,jpcud) = tra(ji,jj,jk,jpcud) + zprocud(ji,jj,jk) * texcretd
                 ENDIF
              END DO
           END DO
        END DO
     ENDIF

     IF( ln_manganese ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji =1 ,jpi
                 IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
              tra(ji,jj,jk,jpdmn) = tra(ji,jj,jk,jpdmn) - ( ( zpromnn(ji,jj,jk) * texcretn ) &
                 &                                      + ( zpromnd(ji,jj,jk) * texcretd ) )
              tra(ji,jj,jk,jpmnn) = tra(ji,jj,jk,jpmnn) + zpromnn(ji,jj,jk) *texcretn
              tra(ji,jj,jk,jpmnd) = tra(ji,jj,jk,jpmnd) + zpromnd(ji,jj,jk) * texcretd
                 ENDIF
              END DO
           END DO
        END DO
     ENDIF

     IF( ln_cobalt ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji =1 ,jpi
                 IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
              tra(ji,jj,jk,jpdco) = tra(ji,jj,jk,jpdco) - ( ( zprocon(ji,jj,jk) * texcretn ) &
                 &                                      + ( zprocod(ji,jj,jk) * texcretd ) ) * 1000. ! 1000x for precision
              tra(ji,jj,jk,jpcon) = tra(ji,jj,jk,jpcon) + zprocon(ji,jj,jk) *texcretn * 1000. ! 1000x for precision
              tra(ji,jj,jk,jpcod) = tra(ji,jj,jk,jpcod) + zprocod(ji,jj,jk) * texcretd * 1000. ! 1000x for precision
                 ENDIF
              END DO
           END DO
        END DO
     ENDIF

     IF( ln_zinc ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji =1 ,jpi
                 IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                 IF( ln_znf ) THEN
                 znf1 = 0.35 * (xlimdia(ji,jj,jk) + rtrn ) / ( xlimdfe(ji,jj,jk)+ rtrn )
                 znf1 = max( 0.03, znf1 )
                 ELSE
                 znf1=znf
                 ENDIF
                 tra(ji,jj,jk,jpdzn) = tra(ji,jj,jk,jpdzn) - ( ( zproznn(ji,jj,jk) * texcretn ) &
                 &                                      + ( zproznd(ji,jj,jk) * texcretd ) )
                 tra(ji,jj,jk,jpznn) = tra(ji,jj,jk,jpznn) + zproznn(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpznd) = tra(ji,jj,jk,jpznd) + zproznd(ji,jj,jk) * (1 - znf1) * texcretd
                 tra(ji,jj,jk,jpzfd) = tra(ji,jj,jk,jpzfd) + zproznd(ji,jj,jk) * znf1 * texcretd
                 ENDIF
              END DO
           END DO
        END DO
     ENDIF

    ! Total primary production per year
    IF( iom_use( "tintpp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
         & tpp = glob_sum( 'p4zprod', ( zprorcan(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) )

    IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          ALLOCATE( zw2d(jpi,jpj), zw3d(jpi,jpj,jpk) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "PPPHYN" ) .OR. iom_use( "PPPHYD" ) )  THEN
              zw3d(:,:,:) = zprorcan(:,:,:) * zfact * tmask(:,:,:)  ! primary production by nanophyto
              CALL iom_put( "PPPHYN"  , zw3d )
              !
              zw3d(:,:,:) = zprorcad(:,:,:) * zfact * tmask(:,:,:)  ! primary production by diatomes
              CALL iom_put( "PPPHYD"  , zw3d )
          ENDIF
          IF( iom_use( "PPNEWN" ) .OR. iom_use( "PPNEWD" ) )  THEN
              zw3d(:,:,:) = zpronewn(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by nanophyto
              CALL iom_put( "PPNEWN"  , zw3d )
              !
              zw3d(:,:,:) = zpronewd(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by diatomes
              CALL iom_put( "PPNEWD"  , zw3d )
          ENDIF
          IF( iom_use( "PBSi" ) )  THEN
              zw3d(:,:,:) = zprorcad(:,:,:) * zfact * tmask(:,:,:) * zysopt(:,:,:) ! biogenic silica production
              CALL iom_put( "PBSi"  , zw3d )
          ENDIF
          IF( iom_use( "PFeN" ) .OR. iom_use( "PFeD" ) )  THEN
              zw3d(:,:,:) = zprofen(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by nanophyto
              CALL iom_put( "PFeN"  , zw3d )
              !
              zw3d(:,:,:) = zprofed(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by  diatomes
              CALL iom_put( "PFeD"  , zw3d )
          ENDIF
          IF( iom_use( "LPRODP" ) )  THEN
              zw3d(:,:,:) = zpligprod1(:,:,:) * 1e9 * zfact * tmask(:,:,:)
              CALL iom_put( "LPRODP"  , zw3d )
          ENDIF
          IF( iom_use( "LDETP" ) )  THEN
              zw3d(:,:,:) = zpligprod2(:,:,:) * 1e9 * zfact * tmask(:,:,:)
              CALL iom_put( "LDETP"  , zw3d )
          ENDIF
          IF( iom_use( "Mumax" ) )  THEN
              zw3d(:,:,:) = zprmaxn(:,:,:) * tmask(:,:,:)   ! Maximum growth rate
              CALL iom_put( "Mumax"  , zw3d )
          ENDIF
          IF( iom_use( "MuN" ) .OR. iom_use( "MuD" ) )  THEN
              zw3d(:,:,:) = zprbio(:,:,:) * xlimphy(:,:,:) * tmask(:,:,:)  ! Realized growth rate for nanophyto
              CALL iom_put( "MuN"  , zw3d )
              !
              zw3d(:,:,:) =  zprdia(:,:,:) * xlimdia(:,:,:) * tmask(:,:,:)  ! Realized growth rate for diatoms
              CALL iom_put( "MuD"  , zw3d )
          ENDIF
          IF( iom_use( "LNlight" ) .OR. iom_use( "LDlight" ) )  THEN
              zw3d(:,:,:) = zprbio (:,:,:) / (zprmaxn(:,:,:) + rtrn) * tmask(:,:,:) ! light limitation term
              CALL iom_put( "LNlight"  , zw3d )
              !
              zw3d(:,:,:) = zprdia (:,:,:) / (zprmaxd(:,:,:) + rtrn) * tmask(:,:,:)  ! light limitation term
              CALL iom_put( "LDlight"  , zw3d )
          ENDIF
          IF( iom_use( "TPP" ) )  THEN
              zw3d(:,:,:) = ( zprorcan(:,:,:) + zprorcad(:,:,:) ) * zfact * tmask(:,:,:)  ! total primary production
              CALL iom_put( "TPP"  , zw3d )
          ENDIF
          IF( iom_use( "PCUN" ) )  THEN
              zw3d(:,:,:) = zprocun(:,:,:) * zfact * tmask(:,:,:)
              CALL iom_put( "PCUN"  , zw3d )
          ENDIF
          IF( iom_use( "PCUD" ) )  THEN
              zw3d(:,:,:) = zprocud(:,:,:) * zfact * tmask(:,:,:)
              CALL iom_put( "PCUD"  , zw3d )
          ENDIF
          IF( iom_use( "PZNN" ) )  THEN
!              zw3d(:,:,:) = vmzn(:,:,:) * zfact * tmask(:,:,:)
              zw3d(:,:,:) = zproznn(:,:,:) * zfact * tmask(:,:,:)
              CALL iom_put( "PZNN"  , zw3d )
          ENDIF
          IF( iom_use( "PZND" ) )  THEN
              zw3d(:,:,:) = zproznd(:,:,:) * zfact * tmask(:,:,:)
              CALL iom_put( "PZND"  , zw3d )
          ENDIF
          IF( iom_use( "MNZNN" ) )  THEN
              zw3d(:,:,:) = mnznn(:,:,:) * tmask(:,:,:)
              CALL iom_put( "MNZNN"  , zw3d )
          ENDIF
          IF( iom_use( "COZND" ) )  THEN
              zw3d(:,:,:) = coznd(:,:,:) * tmask(:,:,:)
              CALL iom_put( "COZND"  , zw3d )
          ENDIF
          IF( iom_use( "MNCOD" ) )  THEN
              zw3d(:,:,:) = mncod(:,:,:) * tmask(:,:,:)
              CALL iom_put( "MNCOD"  , zw3d )
          ENDIF
          IF( iom_use( "MNZND" ) )  THEN
              zw3d(:,:,:) = mnznd(:,:,:) * tmask(:,:,:)
              CALL iom_put( "MNZND"  , zw3d )
          ENDIF
          IF( iom_use( "ZNMNN" ) )  THEN
              zw3d(:,:,:) = znmnn(:,:,:) * tmask(:,:,:)
              CALL iom_put( "ZNMNN"  , zw3d )
          ENDIF
          IF( iom_use( "ZNMND" ) )  THEN
              zw3d(:,:,:) = znmnd(:,:,:) * tmask(:,:,:)
              CALL iom_put( "ZNMND"  , zw3d )
          ENDIF
          IF( iom_use( "PCON" ) )  THEN
              zw3d(:,:,:) = zprocon(:,:,:) * zfact * tmask(:,:,:)
              CALL iom_put( "PCON"  , zw3d )
          ENDIF
          IF( iom_use( "PCOD" ) )  THEN
              zw3d(:,:,:) = zprocod(:,:,:) * zfact * tmask(:,:,:)
              CALL iom_put( "PCOD"  , zw3d )
          ENDIF
          IF( iom_use( "PMNN" ) )  THEN
              zw3d(:,:,:) = zpromnn(:,:,:) * zfact * tmask(:,:,:)
              CALL iom_put( "PMNN"  , zw3d )
          ENDIF
          IF( iom_use( "PMND" ) )  THEN
              zw3d(:,:,:) = zpromnd(:,:,:) * zfact * tmask(:,:,:)
              CALL iom_put( "PMND"  , zw3d )
          ENDIF
          IF( iom_use( "COZNT" ) )  THEN
              zw3d(:,:,:) =zcoznt(:,:,:) * tmask(:,:,:)
              CALL iom_put( "COZNT"  , zw3d )
          ENDIF
          IF( iom_use( "MNT" ) )  THEN
              zw3d(:,:,:) =zmnt(:,:,:) * tmask(:,:,:)
              CALL iom_put( "MNT"  , zw3d )
          ENDIF
          IF( iom_use( "TPNEW" ) )  THEN
              zw3d(:,:,:) = ( zpronewn(:,:,:) + zpronewd(:,:,:) ) * zfact * tmask(:,:,:)  ! total new production
              CALL iom_put( "TPNEW"  , zw3d )
          ENDIF
          IF( iom_use( "TPBFE" ) )  THEN
              zw3d(:,:,:) = ( zprofen(:,:,:) + zprofed(:,:,:) ) * zfact * tmask(:,:,:)  ! total biogenic iron production
              CALL iom_put( "TPBFE"  , zw3d )
          ENDIF
          IF( iom_use( "INTPPPHYN" ) .OR. iom_use( "INTPPPHYD" ) ) THEN  
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
               zw2d(:,:) = zw2d(:,:) + zprorcan(:,:,jk) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk)  ! vert. integrated  primary produc. by nano
             ENDDO
             CALL iom_put( "INTPPPHYN" , zw2d )
             !
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + zprorcad(:,:,jk) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk) ! vert. integrated  primary produc. by diatom
             ENDDO
             CALL iom_put( "INTPPPHYD" , zw2d )
          ENDIF
          IF( iom_use( "INTPP" ) ) THEN   
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + ( zprorcan(:,:,jk) + zprorcad(:,:,jk) ) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk) ! vert. integrated pp
             ENDDO
             CALL iom_put( "INTPP" , zw2d )
          ENDIF
          IF( iom_use( "INTPNEW" ) ) THEN    
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + ( zpronewn(:,:,jk) + zpronewd(:,:,jk) ) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk)  ! vert. integrated new prod
             ENDDO
             CALL iom_put( "INTPNEW" , zw2d )
          ENDIF
          IF( iom_use( "INTPBFE" ) ) THEN           !   total biogenic iron production  ( vertically integrated )
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + ( zprofen(:,:,jk) + zprofed(:,:,jk) ) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk) ! vert integr. bfe prod
             ENDDO
            CALL iom_put( "INTPBFE" , zw2d )
          ENDIF
          IF( iom_use( "INTPBSI" ) ) THEN           !   total biogenic silica production  ( vertically integrated )
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + zprorcad(:,:,jk) * zysopt(:,:,jk) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk)  ! vert integr. bsi prod
             ENDDO
             CALL iom_put( "INTPBSI" , zw2d )
          ENDIF
          IF( iom_use( "tintpp" ) )  CALL iom_put( "tintpp" , tpp * zfact )  !  global total integrated primary production molC/s
          !
          DEALLOCATE( zw2d, zw3d )
       ENDIF
     ENDIF

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
     ENDIF
      !
      IF( ln_timing )  CALL timing_stop('p4z_prod')
      !
   END SUBROUTINE p4z_prod


   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zprod/ pislopen, pisloped, xadap, bresp, excretn, excretd,  &
         &                 chlcnm, chlcdm, chlcmin, fecnm, fecdm, grosip, cupmn, cupmd, &
         &                 znpmn, znpmd, ln_fezn,  copmn, copmd, ln_cozn, kcozn, mnpmn, mnpmd, &
         &                 ln_comnzn, ln_comnzn_simple, mntmaxd, kmnd, kznd, kcoznd, cozntmaxd, &
         &                 kb12d, qmnmaxd, qznmaxd, ln_mnzn_int, ln_mnqzn_int, ln_cozn_int
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_prod_init : phytoplankton growth'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampisprod in reference namelist : Pisces phytoplankton production
      READ  ( numnatp_ref, namp4zprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zprod in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisprod in configuration namelist : Pisces phytoplankton production
      READ  ( numnatp_cfg, namp4zprod, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zprod in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp4zprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zprod'
         WRITE(numout,*) '      mean Si/C ratio                           grosip       =', grosip
         WRITE(numout,*) '      P-I slope                                 pislopen     =', pislopen
         WRITE(numout,*) '      Acclimation factor to low light           xadap        =', xadap
         WRITE(numout,*) '      excretion ratio of nanophytoplankton      excretn      =', excretn
         WRITE(numout,*) '      excretion ratio of diatoms                excretd      =', excretd
         WRITE(numout,*) '      basal respiration in phytoplankton        bresp        =', bresp
         WRITE(numout,*) '      Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
         WRITE(numout,*) '      P-I slope  for diatoms                    pisloped     =', pisloped
         WRITE(numout,*) '      Minimum Chl/C in nanophytoplankton        chlcnm       =', chlcnm
         WRITE(numout,*) '      Minimum Chl/C in diatoms                  chlcdm       =', chlcdm
         WRITE(numout,*) '      Maximum Fe/C in nanophytoplankton         fecnm        =', fecnm
         WRITE(numout,*) '      Minimum Fe/C in diatoms                   fecdm        =', fecdm
         WRITE(numout,*) '    Maximum Cu/P in phytoplankton               cupmn        =', cupmn
         WRITE(numout,*) '    Maximum Cu/P in diatoms                     cupmd        =', cupmd
         WRITE(numout,*) '    Maximum Zn/P in diatoms                     znpmd        =', znpmd
         WRITE(numout,*) '    Enable param. of Fe impacting Zn/P ratio  (T/F) ln_fezn  =', ln_fezn
         WRITE(numout,*) '    Maximum Co/P in phytoplankton                copmn       =', copmn
         WRITE(numout,*) '    Maximum Co/P in diatoms                      copmd       =', copmd
         WRITE(numout,*) '    Dynamic Co-Zn interaction?               (T/F) ln_cozn   =', ln_cozn
         WRITE(numout,*) '    Shapefunction for Co-Zn dependancy, moles Zn, kcozn      =', kcozn
         WRITE(numout,*) '    Maximum Mn/P in phytoplankton                mnpmn       =', mnpmn
         WRITE(numout,*) '    Maximum Mn/P in diatoms                      mnpmd       =', mnpmd
         WRITE(numout,*) '    use full interactive  scheme or not ln_comnzn      =',ln_comnzn
         WRITE(numout,*) '    use simple interactive uptake scheme or not ln_comnzn =',ln_comnzn_simple
         WRITE(numout,*) '    maximum transporters                       mntmaxd =',mntmaxd
         WRITE(numout,*) '    half saturation                              kmnd =',kmnd
         WRITE(numout,*) '    half saturation                              kznd =',kznd
         WRITE(numout,*) '    kcoznd=',kcoznd
         WRITE(numout,*) '    cozntmaxd=',cozntmaxd
         WRITE(numout,*) '    kb12d=',kb12d
         WRITE(numout,*) '    qmnmaxd=',qmnmaxd
         WRITE(numout,*) '    qznmaxd=',qznmaxd
         WRITE(numout,*) ' include transporter interaction between Zn and Mn uptake, ln_mnzn_int =',ln_mnzn_int
         WRITE(numout,*) ' include transporter interaction between Zn and Co uptake, ln_cozn_int =',ln_cozn_int
         WRITE(numout,*) ' include Zn hyperaccumulationn interaction between Zn and Mn uptake, ln_mnqzn_int =',ln_mnqzn_int
      ENDIF
      !
      r1_rday   = 1._wp / rday 
      texcretn  = 1._wp - excretn
      texcretd  = 1._wp - excretd
      tpp       = 0._wp
      !
   END SUBROUTINE p4z_prod_init


   INTEGER FUNCTION p4z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( quotan(jpi,jpj,jpk), quotad(jpi,jpj,jpk), STAT = p4z_prod_alloc )
      !
      IF( p4z_prod_alloc /= 0 ) CALL ctl_warn('p4z_prod_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_prod_alloc

   !!======================================================================
END MODULE p4zprod
