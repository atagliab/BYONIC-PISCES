MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE trcsink         !  General routine to compute sedimentation
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_sink_alloc

   REAL(wp), PUBLIC,ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking, sinking2  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
   REAL(wp), PUBLIC, SAVE,ALLOCATABLE, DIMENSION(:,:,:) ::   sinkingn, sinking2n  !: POC sinking fluxes 
   REAL(wp), PUBLIC,ALLOCATABLE,  SAVE, DIMENSION(:,:,:) ::   sinkingp, sinking2p  !: POC sinking fluxes 
   REAL(wp), PUBLIC, SAVE,ALLOCATABLE, DIMENSION(:,:,:) ::   sinkcal, sinksil   !: CaCO3 and BSi sinking fluxes
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE,DIMENSION(:,:,:) ::   sinkfer            !: Small BFe sinking fluxes
   REAL(wp), PUBLIC,  SAVE,ALLOCATABLE, DIMENSION(:,:,:) ::   sinkfer2           !: Big iron sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE,SAVE, DIMENSION(:,:,:) ::   sinkcu, sinkcu2, sinkscup, sinkscug !: Cam copper model
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkzn, sinkzn2, sinkszp, sinkszg, sinkznf
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkco, sinkco2, sinksco
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkmn, sinkmn2, sinksmn
   INTEGER  :: ik100

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsink.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   'standard sinking parameterisation'                  ???
   !!----------------------------------------------------------------------

   SUBROUTINE p4z_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk
      CHARACTER (len=25) :: charout
      REAL(wp) :: zmax, zfact
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_sink')

      ! Initialization of some global variables
      ! ---------------------------------------
      prodpoc(:,:,:) = 0.
      conspoc(:,:,:) = 0.
      prodgoc(:,:,:) = 0.
      consgoc(:,:,:) = 0.

      !
      !    Sinking speeds of detritus is increased with depth as shown
      !    by data and from the coagulation theory
      !    -----------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / wsbio2scale
               wsbio4(ji,jj,jk) = wsbio2 + MAX(0., ( wsbio2max - wsbio2 )) * zfact
            IF ( ln_cobalt ) wssco(ji,jj,jk) = wsbio + MAX(0.,  10. - wsbio ) * zfact
            IF ( ln_manganese ) wssmn(ji,jj,jk) = wsbio + MAX(0.,  10. - wsbio ) * zfact
            END DO
         END DO
      END DO

      ! limit the values of the sinking speeds to avoid numerical instabilities  
      wsbio3(:,:,:) = wsbio

      !
      !  Initializa to zero all the sinking arrays 
      !   -----------------------------------------
      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      sinkfer (:,:,:) = 0.e0
      sinksil (:,:,:) = 0.e0
      sinkfer2(:,:,:) = 0.e0
 
      IF ( ln_copper ) THEN
      sinkcu(:,:,:) = 0.e0
      sinkcu2(:,:,:) = 0.e0
      sinkscup(:,:,:) = 0.e0
      sinkscug(:,:,:) = 0.e0
      ENDIF
      IF ( ln_zinc ) THEN
      sinkzn(:,:,:) = 0.e0
      sinkzn2(:,:,:) = 0.e0
      sinkszp(:,:,:) = 0.e0
      sinkszg(:,:,:) = 0.e0
      sinkznf(:,:,:) = 0.e0
      ENDIF
      IF ( ln_cobalt ) THEN
      sinkco(:,:,:) = 0.e0
      sinkco2(:,:,:) = 0.e0
      sinksco(:,:,:) = 0.e0
      ENDIF
      IF ( ln_manganese ) THEN
      sinkmn(:,:,:) = 0.e0
      sinkmn2(:,:,:) = 0.e0
      sinksmn(:,:,:) = 0.e0
      ENDIF

      !   Compute the sedimentation term using p4zsink2 for all the sinking particles
      !   -----------------------------------------------------
      CALL trc_sink( kt, wsbio3, sinking , jppoc, rfact2 )
      CALL trc_sink( kt, wsbio3, sinkfer , jpsfe, rfact2 )
      CALL trc_sink( kt, wsbio4, sinking2, jpgoc, rfact2 )
      CALL trc_sink( kt, wsbio4, sinkfer2, jpbfe, rfact2 )
      CALL trc_sink( kt, wsbio4, sinksil , jpgsi, rfact2 )
      CALL trc_sink( kt, wsbio4, sinkcal , jpcal, rfact2 )
         IF ( ln_copper ) THEN

        CALL trc_sink(kt, wsbio3, sinkcu  , jpcup, rfact2 )
        CALL trc_sink( kt, wsbio3, sinkscup , jpscup, rfact2 )

        CALL trc_sink( kt, wsbio4, sinkcu2, jpcug, rfact2 )
        CALL trc_sink( kt, wsbio4, sinkscug, jpscug, rfact2 )
           ENDIF ! ln_copper

         IF ( ln_zinc ) THEN

        CALL trc_sink(kt, wsbio3, sinkzn  , jpznp, rfact2 )
        CALL trc_sink( kt, wsbio3, sinkszp , jpszp, rfact2 )

        CALL trc_sink( kt, wsbio4, sinkzn2, jpzng, rfact2 )
        CALL trc_sink( kt, wsbio4, sinkszg, jpszg, rfact2 )
        CALL trc_sink( kt, wsbio4, sinkznf, jpznf, rfact2 )

           ENDIF ! ln_zinc

           IF (ln_cobalt) THEN 
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / wsbio2scale
!               wssco(ji,jj,jk) = wsbio + MAX(0.,  10. - wsbio ) * zfact
               wssco(ji,jj,jk) = (wsbio/2) + MAX(0.,  10. - wsbio ) * zfact
!               wssco(ji,jj,jk) = (wsbio/2) + MAX(0.,  5. - wsbio ) * zfact
!                  IF( tmask(ji,jj,jk) == 1 ) THEN
!!                    zwsmax = 0.5 * e3t_n(ji,jj,jk) / xstep
!                    wssco(ji,jj,jk) = MIN( wssco(ji,jj,jk), zwsmax * REAL( iiter1, wp ) )
!                  ENDIF
               END DO
            END DO
         END DO
        CALL trc_sink(kt, wsbio3, sinkco  , jpcop, rfact2 )
        CALL trc_sink( kt, wssco, sinksco , jpsco, rfact2 )

        CALL trc_sink( kt, wsbio4, sinkco2, jpcog, rfact2 )
          ENDIF 

           IF (ln_manganese) THEN
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / wsbio2scale
               wssmn(ji,jj,jk) = (wsbio/2) + MAX(0.,  10. - wsbio ) * zfact
               END DO
            END DO
         END DO
        CALL trc_sink(kt, wsbio3, sinkmn  , jpmnp, rfact2 )
        CALL trc_sink( kt, wssmn, sinksmn , jpsmn, rfact2 )
        CALL trc_sink( kt, wsbio4, sinkmn2, jpmng, rfact2 )
          ENDIF
 
      IF( ln_p5z ) THEN
         sinkingn (:,:,:) = 0.e0
         sinking2n(:,:,:) = 0.e0
         sinkingp (:,:,:) = 0.e0
         sinking2p(:,:,:) = 0.e0

         !   Compute the sedimentation term using p4zsink2 for all the sinking particles
         !   -----------------------------------------------------
         CALL trc_sink( kt, wsbio3, sinkingn , jppon, rfact2 )
         CALL trc_sink( kt, wsbio3, sinkingp , jppop, rfact2 )
         CALL trc_sink( kt, wsbio4, sinking2n, jpgon, rfact2 )
         CALL trc_sink( kt, wsbio4, sinking2p, jpgop, rfact2 )
      ENDIF

     ! Total carbon export per year
     IF( iom_use( "tcexp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
        &   t_oce_co2_exp = glob_sum( 'p4zsink', ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * e1e2t(:,:) * tmask(:,:,1) )
     !
     IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          ALLOCATE( zw2d(jpi,jpj), zw3d(jpi,jpj,jpk) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "EPC100" ) )  THEN
              zw2d(:,:) = ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of carbon at 100m
              CALL iom_put( "EPC100"  , zw2d )
          ENDIF
          IF( iom_use( "EPFE100" ) )  THEN
              zw2d(:,:) = ( sinkfer(:,:,ik100) + sinkfer2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of iron at 100m
              CALL iom_put( "EPFE100"  , zw2d )
          ENDIF
          IF( iom_use( "EPCAL100" ) )  THEN
              zw2d(:,:) = sinkcal(:,:,ik100) * zfact * tmask(:,:,1) ! Export of calcite at 100m
              CALL iom_put( "EPCAL100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSI100" ) )  THEN
              zw2d(:,:) =  sinksil(:,:,ik100) * zfact * tmask(:,:,1) ! Export of bigenic silica at 100m
              CALL iom_put( "EPSI100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPC" ) )  THEN
              zw3d(:,:,:) = ( sinking(:,:,:) + sinking2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of carbon in the water column
              CALL iom_put( "EXPC"  , zw3d )
          ENDIF
          IF( iom_use( "EXPFE" ) )  THEN
              zw3d(:,:,:) = ( sinkfer(:,:,:) + sinkfer2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of iron 
              CALL iom_put( "EXPFE"  , zw3d )
          ENDIF
          IF( iom_use( "EXPCAL" ) )  THEN
              zw3d(:,:,:) = sinkcal(:,:,:) * zfact * tmask(:,:,:) ! Export of calcite 
              CALL iom_put( "EXPCAL"  , zw3d )
          ENDIF
          IF( iom_use( "EXPSI" ) )  THEN
              zw3d(:,:,:) = sinksil(:,:,:) * zfact * tmask(:,:,:) ! Export of bigenic silica
              CALL iom_put( "EXPSI"  , zw3d )
          ENDIF
          IF( iom_use( "tcexp" ) )  CALL iom_put( "tcexp" , t_oce_co2_exp * zfact )   ! molC/s
          ! 
IF (ln_copper) THEN
          IF( iom_use( "EPCU100" ) )  THEN
              zw2d(:,:) = ( sinkcu(:,:,ik100) + sinkcu2(:,:,ik100) + sinkscup(:,:,ik100) &
              &      + sinkscug(:,:,ik100) ) * 1.E6 * zfact * tmask(:,:,1) ! Export of zinc
              CALL iom_put( "EPCU100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSCU100" ) )  THEN
              zw2d(:,:) = (sinkscup(:,:,ik100) + sinkscug(:,:,ik100) ) * 1.E6 * zfact * tmask(:,:,1)
              CALL iom_put( "EPSCU100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPCU" ) )  THEN
              zw3d(:,:,:) = ( sinkcu(:,:,:) + sinkcu2(:,:,:) + sinkscup(:,:,:) &
              &      + sinkscug(:,:,:) ) * 1.E6 * zfact * tmask(:,:,:) ! Export of zinc
              CALL iom_put( "EXPCU"  , zw3d )
          ENDIF
ENDIF
IF (ln_zinc) THEN
          IF( iom_use( "EPZNF100" ) )  THEN
              zw2d(:,:) = sinkznf(:,:,ik100) * 1.E6 * zfact * tmask(:,:,1)
              CALL iom_put( "EPZNF100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSZN100" ) )  THEN
              zw2d(:,:) = (sinkszp(:,:,ik100) + sinkszg(:,:,ik100) ) * 1.E6 * zfact * tmask(:,:,1)
              CALL iom_put( "EPSZN100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPZN" ) )  THEN
              zw3d(:,:,:) = ( sinkzn(:,:,:) + sinkzn2(:,:,:) + sinkszp(:,:,:) &
              &      + sinkznf(:,:,:) + sinkszg(:,:,:) ) * 1.E6 * zfact * tmask(:,:,:) ! Export of zinc
              CALL iom_put( "EXPZN"  , zw3d )
          ENDIF
ENDIF
IF (ln_cobalt) THEN 
          IF( iom_use( "EPCO100" ) )  THEN
              zw2d(:,:) = ( sinkco(:,:,ik100) + sinkco2(:,:,ik100) + sinksco(:,:,ik100) &
              &      ) * 1.E6 * zfact * tmask(:,:,1) ! Export of cobalt
              CALL iom_put( "EPCO100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSCO100" ) )  THEN
              zw2d(:,:) = sinksco(:,:,ik100) * 1.E6 * zfact * tmask(:,:,1)
              CALL iom_put( "EPSCO100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPCO" ) )  THEN
              zw3d(:,:,:) = ( sinkco(:,:,:) + sinkco2(:,:,:) + sinksco(:,:,:) &
              &      ) * 1.E6 * zfact * tmask(:,:,:) ! Export of cobalt
              CALL iom_put( "EXPCO"  , zw3d )
          ENDIF
ENDIF
IF (ln_manganese) THEN
          IF( iom_use( "EPMN100" ) )  THEN
              zw2d(:,:) = ( sinkmn(:,:,ik100) + sinkmn2(:,:,ik100) + sinksmn(:,:,ik100) &
              &      ) * 1.E6 * zfact * tmask(:,:,1) ! Export of manganese
              CALL iom_put( "EPMN100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSMN100" ) )  THEN
              zw2d(:,:) = sinksmn(:,:,ik100) * 1.E6 * zfact * tmask(:,:,1)
              CALL iom_put( "EPSMN100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPMN" ) )  THEN
              zw3d(:,:,:) = ( sinkmn(:,:,:) + sinkmn2(:,:,:) + sinksmn(:,:,:) &
              &      ) * 1.E6 * zfact * tmask(:,:,:) ! Export of manganese
              CALL iom_put( "EXPMN"  , zw3d )
          ENDIF
ENDIF

          DEALLOCATE( zw2d, zw3d )
        ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink


   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!----------------------------------------------------------------------
      INTEGER :: jk
      !!----------------------------------------------------------------------
      !
      ik100 = 10        !  last level where depth less than 100 m
      DO jk = jpkm1, 1, -1
         IF( gdept_1d(jk) > 100. )  ik100 = jk - 1
      END DO
      IF (lwp) WRITE(numout,*)
      IF (lwp) WRITE(numout,*) ' Level corresponding to 100m depth ',  ik100 + 1
      IF (lwp) WRITE(numout,*)
      !
      t_oce_co2_exp = 0._wp
      !
   END SUBROUTINE p4z_sink_init

   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(3)
      !!----------------------------------------------------------------------
      !
      ierr(:) = 0
      !
      ALLOCATE( sinking(jpi,jpj,jpk) , sinking2(jpi,jpj,jpk)                    ,     &                
         &      sinkcal(jpi,jpj,jpk) , sinksil (jpi,jpj,jpk)                    ,     &                
         &      sinkfer2(jpi,jpj,jpk)                                           ,     &                
         &      sinkfer(jpi,jpj,jpk)                                            , STAT=ierr(1) )                
         !
      IF( ln_p5z    ) ALLOCATE( sinkingn(jpi,jpj,jpk), sinking2n(jpi,jpj,jpk)   ,     &
         &                      sinkingp(jpi,jpj,jpk), sinking2p(jpi,jpj,jpk)   , STAT=ierr(2) )
      
      p4z_sink_alloc = MAXVAL( ierr )
      IF( p4z_sink_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_sink_alloc : failed to allocate arrays.' )
      !
     IF ( ln_copper ) ALLOCATE( sinkcu(jpi,jpj,jpk), sinkcu2(jpi,jpj,jpk)   , &
         &                     sinkscup(jpi,jpj,jpk), sinkscug(jpi,jpj,jpk) , STAT=ierr(3) )
      !
      IF ( ln_zinc ) ALLOCATE( sinkzn(jpi,jpj,jpk), sinkzn2(jpi,jpj,jpk)   , &
         &                     sinkszp(jpi,jpj,jpk), sinkszg(jpi,jpj,jpk), sinkznf(jpi,jpj,jpk) , STAT=ierr(4) )
      IF ( ln_cobalt ) ALLOCATE( sinkco(jpi,jpj,jpk), sinkco2(jpi,jpj,jpk) ,     &
         &                       sinksco(jpi,jpj,jpk) , STAT=ierr(5) )
      IF ( ln_manganese ) ALLOCATE( sinkmn(jpi,jpj,jpk), sinkmn2(jpi,jpj,jpk) , &
         &                       sinksmn(jpi,jpj,jpk) , STAT=ierr(5) )

   END FUNCTION p4z_sink_alloc
   
   !!======================================================================
END MODULE p4zsink
