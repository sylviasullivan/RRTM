MODULE mo_lrtm_rtrnmr

! >> sylvia_20200318
! Convert this guy back into a module also
USE mo_lrtm_setup,   ONLY: ntbl, bpade, exp_tbl, tfn_tbl, tau_tbl
! >> sylvia_20200402, Adding back in mo_lrtm_par.f90
USE mo_lrtm_par,     ONLY : nbndlw, delwave, ngs
! << sylvia_20200402

! >> sylvia_20200318
! Set up the wp type that is often used below. nbndlw pulled from
! mo_lrtm_par.f90
INTEGER, PARAMETER   :: pd =  12
INTEGER, PARAMETER   :: rd = 307
INTEGER, PARAMETER   :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
INTEGER, PARAMETER   :: wp = dp
REAL (wp), PARAMETER ::  pi = 3.14159265358979323846264338327950288_wp
REAL(wp), PARAMETER  :: fluxfac = 2.0e+04_wp * pi
REAL(wp), PARAMETER  :: tblint = REAL(ntbl, wp)
! << sylvia_20200318

CONTAINS
  SUBROUTINE lrtm_rtrnmr(&
                & nlayers, istart, iend, iout, semiss, ncbands, &
                & cldfrac, taucloud, planklay, planklev, plankbnd, &
                & pwvcm, fracs, taut, &
                & totuflux, totdflux, fnet, &
                & totuclfl, totdclfl, fnetc, &
                & idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )
        !-----------------------------------------------------------------------------
        !
        !  Original version:   E. J. Mlawer, et al. RRTM_V3.0
        !  Revision for GCMs:  Michael J. Iacono; October, 2002
        !  Revision for F90:  Michael J. Iacono; June, 2006
        !  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, November 2009
        !  Version for In-Order architectures: Thomas Jahns, DKRZ 2016
        !
        !  This program calculates the upward fluxes, downward fluxes, and
        !  heating rates for an arbitrary clear or cloudy atmosphere.  The input
        !  to this program is the atmospheric profile, all Planck function
        !  information, and the cloud fraction by layer.  A variable diffusivity
        !  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9
        !  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of
        !  the column water vapor, and other bands use a value of 1.66.  The Gaussian
        !  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that
        !  use of the emissivity angle for the flux integration can cause errors of
        !  1 to 4 W/m2 within cloudy layers.
        !  Clouds are treated with a maximum-random cloud overlap method.
        !  This subroutine also provides the optional capability to calculate
        !  the derivative of upward flux respect to surface temperature using
        !  the pre-tabulated derivative of the Planck function with respect to
        !  temperature integrated over each spectral band.
        !
        !  Setting the preprocessor switch LRTM_FULL_VECTORIZATION selects
        !  an equivalent (bit-reproducible results with NAG Fortran)
        !  algorithm that vectorizes much more thoroughly and gives
        !  parallelization-invariant results on platforms with FMA where
        !  the divergent branches of the default implementation lead to
        !  different contractions per branch. But it's only faster on
        !  architectures where branches and loop overhead are relatively
        !  expensive compared to computation and memory accesses. The
        !  alternative, incompletely-vectorized version was found to run
        !  faster on Intel Sandy Bridge and Ivy Bridge Xeons and about
        !  equally fast on later AVX2 CPUs (Haswell/Broadwell).

        !***************************************************************************

        ! ------- Declarations -------

        ! ----- Input -----
        INTEGER, INTENT(in) :: nlayers         ! total number of layers
        INTEGER, INTENT(in) :: istart          ! beginning band of calculation
        INTEGER, INTENT(in) :: iend            ! ending band of calculation
        INTEGER, INTENT(in) :: iout            ! output option flag

        ! Atmosphere
        REAL(wp), INTENT(in) :: pwvcm           ! precipitable water vapor (cm)
        !    Dimensions: (0:nlayers)
        REAL(wp), INTENT(in) :: semiss(:)        ! lw surface emissivity
        !    Dimensions: (nbndlw)
        REAL(wp), INTENT(in) :: planklay(:,:)      !
        !    Dimensions: (nlayers,nbndlw)
        REAL(wp), INTENT(in) :: planklev(0:,:)     !
        !    Dimensions: (0:nlayers,nbndlw)
        REAL(wp), INTENT(in) :: plankbnd(:)        !
        !    Dimensions: (nbndlw)
        REAL(wp), INTENT(in) :: fracs(:,:)         !
        !    Dimensions: (nlayers,ngptw)
        REAL(wp), INTENT(in) :: taut(:,:)          ! gaseous + aerosol optical depths
        !    Dimensions: (nlayers,ngptlw)

        ! Clouds
        INTEGER, INTENT(in) :: ncbands          ! number of cloud spectral bands
        REAL(wp), INTENT(in) :: cldfrac(:)       ! layer cloud fraction
        !    Dimensions: (nlayers)
        REAL(wp), INTENT(in) :: taucloud(:,:)      ! layer cloud optical depth
        !    Dimensions: (nlayers,nbndlw)
        INTEGER, INTENT(in) :: idrv            ! flag for calculation of dF/dt from
        ! Planck derivative [0=off, 1=on]
        REAL(wp), INTENT(in) :: dplankbnd_dt(:)    ! derivative of Planck function wrt temp
        !    Dimensions: (nbndlw)

        ! ----- Output -----
        REAL(wp), INTENT(out) :: totuflux(0:)      ! upward longwave flux (w/m2)
        !    Dimensions: (0:nlayers)
        REAL(wp), INTENT(out) :: totdflux(0:)      ! downward longwave flux (w/m2)
        !    Dimensions: (0:nlayers)
        REAL(wp), INTENT(out) :: fnet(0:)          ! net longwave flux (w/m2)
        !    Dimensions: (0:nlayers)
        REAL(wp), INTENT(out) :: totuclfl(0:)      ! clear sky upward longwave flux (w/m2)
        !    Dimensions: (0:nlayers)
        REAL(wp), INTENT(out) :: totdclfl(0:)      ! clear sky downward longwave flux (w/m2)
        !    Dimensions: (0:nlayers)
        REAL(wp), INTENT(out) :: fnetc(0:)         ! clear sky net longwave flux (w/m2)
        !    Dimensions: (0:nlayers)
        REAL(wp), INTENT(out) :: dtotuflux_dt(0:)  ! change in upward longwave flux (w/m2/k)
        ! with respect to surface temperature
        !    Dimensions: (0:nlayers)
        REAL(wp), INTENT(out) :: dtotuclfl_dt(0:)  ! change in upward longwave flux (w/m2/k)
        ! with respect to surface temperature
        !    Dimensions: (0:nlayers)
        
        ! >> sylvia_20200310
        ! Commenting out the lines below.
        !#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
        !CONTIGUOUS :: semiss, planklay, planklev, plankbnd, fracs, taut, &
        !        cldfrac, taucloud, dplankbnd_dt, totuflux, totdflux, fnet, totuclfl, &
        !        fnetc, dtotuflux_dt, dtotuclfl_dt
        !#endif
        ! << sylvia_20200310

        
        ! Vectorized version implemented by Guenther Zaengl, DWD
        ! See above for variable descriptions
        ! ----- Local -----
        ! Declarations for radiative transfer
        REAL(wp) :: atot(nlayers)
        REAL(wp) :: atrans(nlayers)
        REAL(wp) :: bbugas(nlayers)
        REAL(wp) :: bbutot(nlayers)
        REAL(wp) :: clrurad(0:nlayers)
        REAL(wp) :: clrdrad(0:nlayers)
        REAL(wp) :: uflux
        REAL(wp) :: dflux
        REAL(wp) :: urad(0:nlayers)
        REAL(wp) :: drad(0:nlayers)
        REAL(wp) :: uclfl
        REAL(wp) :: dclfl

        REAL(wp) :: secdiff(nbndlw)          ! secant of diffusivity angle
        REAL(wp) :: dplankup(nlayers), dplankdn(nlayers)
        REAL(wp), PARAMETER :: rec_6 = 0.166667_wp, &
                ! diffusivity angle adjustment coefficients
        a0(nbndlw) = (/ 1.66_wp,  1.55_wp,  1.58_wp,  1.66_wp, &
                & 1.54_wp, 1.454_wp,  1.89_wp,  1.33_wp, &
                & 1.668_wp,  1.66_wp,  1.66_wp,  1.66_wp, &
                & 1.66_wp,  1.66_wp,  1.66_wp,  1.66_wp /),&
                a1(nbndlw) = (/ 0.00_wp,  0.25_wp,  0.22_wp,  0.00_wp, &
                & 0.13_wp, 0.446_wp, -0.10_wp,  0.40_wp, &
                & -0.006_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
                & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /),&
                a2(nbndlw) = (/ 0.00_wp, -12.0_wp, -11.7_wp,  0.00_wp, &
                & -0.72_wp,-0.243_wp,  0.19_wp,-0.062_wp, &
                & 0.414_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
                & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /), &
                ! This secant and weight corresponds to the standard diffusivity
        ! angle.  This initial value is redefined below for some bands.
        wtdiff = 0.5_wp
        REAL(wp) :: radld, radclrd, plfrac
        REAL(wp) :: odepth, odtot, gassrc, ttot
        REAL(wp) :: bbd
        REAL(wp) :: rad0, reflect, radlu, radclru

        REAL(wp) :: duflux_dt
        REAL(wp) :: duclfl_dt
        REAL(wp) :: d_urad_dt(0:nlayers)
        REAL(wp) :: d_clrurad_dt(0:nlayers)
        REAL(wp) :: d_rad0_dt, d_radlu_dt, d_radclru_dt

        LOGICAL :: lcldlyr(0:nlayers+1)             ! flag for cloud in layer
        ! >> sylvia_20200318, Remove ib below.
        INTEGER :: ibnd, iband, lay, lev        ! loop indices
        INTEGER :: igc                          ! g-point interval counter
        LOGICAL :: iclddn                       ! flag for cloud in down path
        INTEGER :: ittot, itgas                 ! lookup table indices
        INTEGER, PARAMETER :: ipat(16,0:2) = RESHAPE((/&
             ! These arrays indicate the spectral 'region' (used in the
             ! calculation of ice cloud optical depths) corresponding
             ! to each spectral band.  See cldprop.f for more details.
             & 1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
             & 1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/), (/16, 3/))

        INTEGER :: ibv

        ! Declarations for cloud overlap adjustment
        REAL(wp) :: faccld1(nlayers+1),faccld2(nlayers+1)
        REAL(wp) :: facclr1(nlayers+1),facclr2(nlayers+1)
        REAL(wp) :: faccmb1(nlayers+1),faccmb2(nlayers+1)
        ! >> sylvia_20200318
        ! Taking out the 0: from the arrays below.
        REAL(wp) :: faccld1d(0:nlayers),faccld2d(0:nlayers)
        REAL(wp) :: facclr1d(0:nlayers),facclr2d(0:nlayers)
        REAL(wp) :: faccmb1d(0:nlayers),faccmb2d(0:nlayers)
        !--------------------------------------------------------------------------
        ! Maximum/Random cloud overlap variables
        ! for upward radiative transfer
        !  facclr2  fraction of clear radiance from previous layer that needs to
        !           be switched to cloudy stream
        !  facclr1  fraction of the radiance that had been switched in the previous
        !           layer from cloudy to clear that needs to be switched back to
        !           cloudy in the current layer
        !  faccld2  fraction of cloudy radiance from previous layer that needs to
        !           be switched to clear stream
        !  faccld1  fraction of the radiance that had been switched in the previous
        !           layer from clear to cloudy that needs to be switched back to
        !           clear in the current layer
        ! for downward radiative transfer
        !  facclr2d fraction of clear radiance from previous layer that needs to
        !           be switched to cloudy stream
        !  facclr1d fraction of the radiance that had been switched in the previous
        !           layer from cloudy to clear that needs to be switched back to
        !           cloudy in the current layer
        !  faccld2d fraction of cloudy radiance from previous layer that needs to
        !           be switched to clear stream
        !  faccld1d fraction of the radiance that had been switched in the previous
        !           layer from clear to cloudy that needs to be switched back to
        !           clear in the current layer
        !--------------------------------------------------------------------------

        REAL(wp) :: cldsrc, oldclr, oldcld, radmod
        ! >> sylvia_20200310
        ! Below I am assuming that LRTM_FULL_VECTORIZATION is not defined.
        !#ifndef LRTM_FULL_VECTORIZATION
        ! >> sylvia_20200318
        ! Uncommenting the two lines below again....
        REAL(wp) :: odepth_rec, odtot_rec, tblind, tfactot, tfacgas, bbdtot, &
                transc, tausfac
        !! >> sylvia_20200306
        !! Changing from DIMENSION:: var(kproma) to REAL(wp) here.
        !REAL(wp) :: gassrc, cldsrc, bbdtot, oldcld, oldclr, odepth, odtot, radmod
        !! << sylvia_20200306
        !INTEGER :: itr, icld1, npoints1, npoints2, npoints3, ilist1, ilist2, ilist3
        !! >> sylvia_20200306
        !! Removing the kproma dimension from the ilist* variables.
        !#else
        ! >> sylvia_20200318
        ! Remove unused odepth_rec_or_tfacgas and odtot_rec_or_tfac_tot below.
        !REAL(wp) :: odepth_rec_or_tfacgas, odtot_rec_or_tfactot, &
        !REAL(wp) :: clrradd_temp, cldradd_temp, >> sylvia_20200318
        !#endif
        ! << sylvia_20200318

        ! >> sylvia_20200306
        ! Removing the DIMENSION(kproma) below.
        REAL(wp) :: clrradd, cldradd, clrradu, cldradu, rad
        ! << sylvia_20200306

        ! >> sylvia_20200310
        ! Again assuming that LRTM_FULL_VECTORIZATION is not defined
        !#ifdef LRTM_FULL_VECTORIZATION
        !LOGICAL :: branch_od1, branch_od2
        !#endif
        ! << sylvia_20200310

        ! ------- Definitions -------
        ! input
        !    nlayers                      ! number of model layers
        !    ngptlw                       ! total number of g-point subintervals
        !    nbndlw                       ! number of longwave spectral bands
        !    ncbands                      ! number of spectral bands for clouds
        !    secdiff                      ! diffusivity angle
        !    wtdiff                       ! weight for radiance to flux conversion
        !    pavel                        ! layer pressures (mb)
        !    tavel                        ! layer temperatures (k)
        !    tz                           ! level (interface) temperatures(mb)
        !    tbound                       ! surface temperature (k)
        !    cldfrac                      ! layer cloud fraction
        !    taucloud                     ! layer cloud optical depth
        !    lcldlyr                      ! flag for cloudy layers
        !    iclddn                       ! flag for cloud in column at any layer
        !    semiss                       ! surface emissivities for each band
        !    reflect                      ! surface reflectance
        !    bpade                        ! 1/(pade constant)
        !    exp_tbl                      ! exponential look-up table for transmittance
        !    tfn_tbl                      ! tau transition function look-up table

        ! local
        !    atrans                       ! gaseous absorptivity
        !    atot                         ! combined gaseous and cloud absorptivity
        !    odclr                        ! clear sky (gaseous) optical depth
        !    odcld                        ! cloud optical depth
        !    odtot                        ! optical depth of gas and cloud
        !    tfacgas                      ! gas-only pade factor, used for planck fn
        !    tfactot                      ! gas and cloud pade factor, used for planck fn
        !    bbdgas                       ! gas-only planck function for downward rt
        !    bbugas                       ! gas-only planck function for upward rt
        !    bbdtot                       ! gas and cloud planck function for downward rt
        !    bbutot                       ! gas and cloud planck function for upward calc.
        !    gassrc                       ! source radiance due to gas only
        !    radlu                        ! spectrally summed upward radiance
        !    radclru                      ! spectrally summed clear sky upward radiance
        !    urad                         ! upward radiance by layer
        !    clrurad                      ! clear sky upward radiance by layer
        !    radld                        ! spectrally summed downward radiance
        !    radclrd                      ! spectrally summed clear sky downward radiance
        !    drad                         ! downward radiance by layer
        !    clrdrad                      ! clear sky downward radiance by layer
        !    d_radlu_dt                   ! spectrally summed upward radiance
        !    d_radclru_dt                 ! spectrally summed clear sky upward radiance
        !    d_urad_dt                    ! upward radiance by layer
        !    d_clrurad_dt                 ! clear sky upward radiance by layer

        ! output
        !    totuflux                     ! upward longwave flux (w/m2)
        !    totdflux                     ! downward longwave flux (w/m2)
        !    fnet                         ! net longwave flux (w/m2)
        !    totuclfl                     ! clear sky upward longwave flux (w/m2)
        !    totdclfl                     ! clear sky downward longwave flux (w/m2)
        !    fnetc                        ! clear sky net longwave flux (w/m2)
        !    dtotuflux_dt                 ! change in upward longwave flux (w/m2/k)
        !                                 ! with respect to surface temperature
        !    dtotuclfl_dt                 ! change in clear sky upward longwave flux (w/m2/k)
        !

        ! Local variables for cloud / no cloud index lists
        !    icld_ind(1:n_cloudpoints(lev),lev) stores indices of cloudy points
        ! >> sylvia_20200306
        ! Commenting out the icld_ind below.
        !INTEGER :: icld_ind(nlayers)
        INTEGER :: iclear, n_cloudpoints(nlayers) !, icld
        ! << sylvia_20200306


        ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
        ! and 1.80) as a function of total column water vapor.  The function
        ! has been defined to minimize flux and cooling rate errors in these bands
        ! over a wide range of precipitable water values.
        ! >> sylvia_20200325
        ! Set the ngs values from mo_lrtm_setup.
        ngs(:) = (/10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,140/)
        ! << sylvia_20200325
        DO ibnd = 1,nbndlw
            IF (ibnd.EQ.1 .OR. ibnd.EQ.4 .OR. ibnd.GE.10) THEN
                secdiff(ibnd) = 1.66_wp
            ELSE
                secdiff(ibnd) = a0(ibnd) + a1(ibnd)*EXP(a2(ibnd)*pwvcm)
                IF (secdiff(ibnd) .GT. 1.80_wp) secdiff(ibnd) = 1.80_wp
                IF (secdiff(ibnd) .LT. 1.50_wp) secdiff(ibnd) = 1.50_wp
            ENDIF
        ENDDO

        !CDIR BEGIN COLLAPSE
        ! >> sylvia_20200306 For all of the following (:,:) changed to (:).
        urad(:)     = 0.0_wp
        drad(:)     = 0.0_wp
        totuflux(:) = 0.0_wp
        totdflux(:) = 0.0_wp
        clrurad(:)  = 0.0_wp
        clrdrad(:)  = 0.0_wp
        totuclfl(:) = 0.0_wp
        totdclfl(:) = 0.0_wp
        ! << sylvia_20200306
        !CDIR END

        IF (idrv .EQ. 1) THEN
            !CDIR BEGIN COLLAPSE
            ! >> sylvia_20200306 For all of the following (:,:) changed to (:).
            d_urad_dt(:)    = 0.0_wp
            d_clrurad_dt(:) = 0.0_wp
            dtotuflux_dt(:) = 0.0_wp
            dtotuclfl_dt(:) = 0.0_wp
            !CDIR END
        ENDIF

        lcldlyr(0) = .FALSE.
        DO lay = 1, nlayers

           icld   = 0
           iclear = 0

           ! >> sylvia_20200306
           ! Remove the icld_ind statements below. Only use this to set icld, iclear, and n_cloudpoints.
           IF (cldfrac(lay) .GE. 1.e-6_wp) THEN
               lcldlyr(lay) = .TRUE.
               icld = icld + 1
               !icld_ind(icld,lay) = jl
           ELSE
                lcldlyr(lay) = .FALSE.
                !icld_ind(kproma - iclear,lay) = jl
                iclear = iclear + 1
           ENDIF
           ! << sylvia_20200306
           n_cloudpoints(lay) = icld

        ENDDO
        lcldlyr(nlayers+1) = .FALSE.

        ! Maximum/Random cloud overlap parameter

        ! >> sylvia_20200306
        ! Removing the sixth input here (icld_ind)
        ! >> sylvia_20200318
        ! Removing the cloudpoints input.
        CALL cloud_overlap( nlayers, 1, cldfrac, &
                lcldlyr, 1, nlayers, 1, &
                faccld1, faccld2, facclr1, facclr2, faccmb1, faccmb2)

        ! istcld(i,j) = .true. if j == 1 or .not. lcldlyr(i, j - 1)
        ! istcldd(i,j) = .true. if j == nlayers or .not. lcldlyr(i, j + 1)

        CALL cloud_overlap( nlayers, 0, cldfrac, &
                lcldlyr, nlayers, 1, -1, &
                faccld1d, faccld2d, facclr1d, facclr2d, faccmb1d, faccmb2d)
        ! << sylvia_20200306

        igc = 1
        ! Loop over frequency bands.
        DO iband = istart, iend
          IF (ncbands .EQ. 1) THEN
              ibv = ipat(iband,0)
          ELSEIF (ncbands .EQ. 16) THEN
              ibv = ipat(iband,2)
          ELSEIF (ncbands .EQ.  5) THEN
              ibv = ipat(iband,1)
          ENDIF

          DO lev = 1, nlayers
             !DIR$ SIMD
             dplankup(lev) = planklev(lev,iband) - planklay(lev,iband)
             dplankdn(lev) = planklev(lev-1,iband) - planklay(lev,iband)
          ENDDO

          ! Reinitialize g-point counter for each band if output for each band is requested.
          IF (iout.GT.0.AND.iband.GE.2) igc = ngs(iband-1)+1

          ! Loop over g-channels.
          1000    CONTINUE

          ! Radiative transfer starts here.
          radld   = 0._wp
          radclrd = 0._wp
          iclddn  = .FALSE.
          ! >> sylvia_20200310
          !#ifdef LRTM_FULL_VECTORIZATION
          !clrradd = 0._wp
          !cldradd = 0._wp
          !#endif
          ! << sylvia_20200310

          ! Downward radiative transfer loop.
          DO lev = nlayers, 1, -1
            ! >> sylvia_20200310 Removing the ~100 lines below with the assumption that 
            ! LRTM_FULL_VECTORIZATION is not defined.
            ! >> sylvia_20200306
            ! Here we can simplify the IF statement by removing the option with both cloudy 
            ! and clear pixels. In the initial IF n_cloudpoints(lev) ==  to IF n_cloudpoint(lev) == 1.
            IF (n_cloudpoints(lev) == 1) THEN ! all points are cloudy

                plfrac = fracs(lev,igc)
                odepth = MAX(0.0_wp, secdiff(iband) * taut(lev,igc))

                iclddn = .TRUE.
                odtot = odepth + secdiff(ibv) * taucloud(lev,ibv)
                IF (odtot .LT. 0.06_wp) THEN
                    atrans(lev) = odepth - 0.5_wp*odepth*odepth
                    odepth_rec = rec_6*odepth
                    gassrc = plfrac*(planklay(lev,iband) &
                            + dplankdn(lev)*odepth_rec)*atrans(lev)

                    atot(lev) = odtot - 0.5_wp*odtot*odtot
                    odtot_rec = rec_6*odtot
                    bbdtot =  plfrac * (planklay(lev,iband)+dplankdn(lev)*odtot_rec)
                    bbd = plfrac*(planklay(lev,iband)+dplankdn(lev)*odepth_rec)

                    bbugas(lev) =  plfrac * (planklay(lev,iband)+dplankup(lev)*odepth_rec)
                    bbutot(lev) =  plfrac * (planklay(lev,iband)+dplankup(lev)*odtot_rec)
                ELSEIF (odepth .LE. 0.06_wp) THEN
                    atrans(lev) = odepth - 0.5_wp*odepth*odepth
                    odepth_rec = rec_6*odepth
                    gassrc = plfrac*(planklay(lev,iband) &
                            + dplankdn(lev)*odepth_rec)*atrans(lev)

                    tblind = odtot/(bpade+odtot)
                    ittot = INT(tblint*tblind + 0.5_wp)
                    tfactot = tfn_tbl(ittot)
                    bbdtot = plfrac * (planklay(lev,iband) + tfactot*dplankdn(lev))
                    bbd = plfrac*(planklay(lev,iband)+dplankdn(lev)*odepth_rec)
                    atot(lev) = 1._wp - exp_tbl(ittot)

                    bbugas(lev) = plfrac * (planklay(lev,iband) + dplankup(lev)*odepth_rec)
                    bbutot(lev) = plfrac * (planklay(lev,iband) + tfactot * dplankup(lev))
                ELSE
                    tblind = odepth/(bpade+odepth)
                    itgas = INT(tblint*tblind+0.5_wp)
                    odepth = tau_tbl(itgas)
                    atrans(lev) = 1._wp - exp_tbl(itgas)
                    tfacgas = tfn_tbl(itgas)
                    gassrc = plfrac * (planklay(lev,iband) &
                            + tfacgas*dplankdn(lev)) * atrans(lev)

                    tblind = odtot/(bpade+odtot)
                    ittot = INT(tblint*tblind + 0.5_wp)
                    tfactot = tfn_tbl(ittot)
                    bbdtot = plfrac * (planklay(lev,iband) + tfactot*dplankdn(lev))
                    bbd = plfrac*(planklay(lev,iband)+tfacgas*dplankdn(lev))
                    atot(lev) = 1._wp - exp_tbl(ittot)

                    bbugas(lev) = plfrac * (planklay(lev,iband) + tfacgas * dplankup(lev))
                    bbutot(lev) = plfrac * (planklay(lev,iband) + tfactot * dplankup(lev))
                ENDIF

                IF (.NOT. lcldlyr(lev+1)) THEN
                    cldradd = cldfrac(lev) * radld
                    clrradd = radld - cldradd
                    oldcld = cldradd
                    oldclr = clrradd
                    rad = 0._wp
                ENDIF
                ttot = 1._wp - atot(lev)
                cldsrc = bbdtot * atot(lev)
                cldradd = cldradd * ttot + cldfrac(lev) * cldsrc
                clrradd = clrradd * (1._wp-atrans(lev)) + &
                    & (1._wp-cldfrac(lev))*gassrc
                radld = cldradd + clrradd
                drad(lev-1) = drad(lev-1) + radld

                radmod = rad * &
                    & (facclr1d(lev-1) * (1._wp-atrans(lev)) + &
                    & faccld1d(lev-1) *  ttot) - &
                    & faccmb1d(lev-1) * gassrc + &
                    & faccmb2d(lev-1) * cldsrc

                oldcld = cldradd - radmod
                oldclr = clrradd + radmod
                rad = -radmod + facclr2d(lev-1)*oldclr -&
                    &  faccld2d(lev-1)*oldcld
                cldradd = cldradd + rad
                clrradd = clrradd - rad

            ELSE 
                ! >> sylvia_20200306
                ! This now corresponds to the ELSE IF (n_cloudpoints(lev) == 0) THEN ! all points are clear

                plfrac = fracs(lev,igc)
                odepth = MAX(0.0_wp, secdiff(iband) * taut(lev,igc))

                IF (odepth .LE. 0.06_wp) THEN
                    atrans(lev) = odepth-0.5_wp*odepth*odepth
                    odepth_rec = rec_6*odepth
                    bbd = plfrac*(planklay(lev,iband)+dplankdn(lev)*odepth_rec)
                    bbugas(lev) = plfrac*(planklay(lev,iband)+dplankup(lev)*odepth_rec)
                ELSE
                    tblind = odepth/(bpade+odepth)
                    itr = INT(tblint*tblind+0.5_wp)
                    transc = exp_tbl(itr)
                    atrans(lev) = 1._wp-transc
                    tausfac = tfn_tbl(itr)
                    bbd = plfrac*(planklay(lev,iband)+tausfac*dplankdn(lev))
                    bbugas(lev) = plfrac * (planklay(lev,iband) + tausfac * dplankup(lev))
                ENDIF
                radld = radld + (bbd-radld)*atrans(lev)
                drad(lev-1) = drad(lev-1) + radld

            ENDIF

            !  Set clear sky stream to total sky stream as long as layers
            !  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
            !  and clear sky stream must be computed separately from that point.
            ! >> sylvia_20200306
            ! Remove the loop over nproma below.
            IF (iclddn) THEN
                radclrd = radclrd + (bbd-radclrd) * atrans(lev)
                clrdrad(lev-1) = clrdrad(lev-1) + radclrd
            ELSE
                radclrd = radld
                clrdrad(lev-1) = drad(lev-1)
            ENDIF
            ! << sylvia_20200306 and sylvia_20200310 for #endif below
            !#endif
          ENDDO

          ! >> sylvia_20200306
          ! Remove the loop over nproma below.
          ! Spectral emissivity & reflectance
          !  Include the contribution of spectrally varying longwave emissivity
          !  and reflection from the surface to the upward radiative transfer.
          !  Note: Spectral and Lambertian reflection are identical for the
          !  diffusivity angle flux integration used here.
          !  Note: The emissivity is applied to plankbnd and dplankbnd_dt when
          !  they are defined in subroutine setcoef.

          rad0 = fracs(1,igc) * plankbnd(iband)

          !  Add in reflection of surface downward radiance.
          reflect = 1._wp - semiss(iband)
          radlu = rad0 + reflect * radld
          radclru = rad0 + reflect * radclrd

          ! Upward radiative transfer loop.

          urad(0) = urad(0) + radlu
          clrurad(0) = clrurad(0) + radclru

          IF (idrv .EQ. 1) THEN
            !DIR$ SIMD
            d_rad0_dt = fracs(1,igc) * dplankbnd_dt(iband)
            d_radlu_dt = d_rad0_dt
            d_urad_dt(0) = d_urad_dt(0) + d_radlu_dt
            d_radclru_dt = d_rad0_dt
            d_clrurad_dt(0) = d_clrurad_dt(0) + d_radclru_dt
          ENDIF

          DO lev = 1, nlayers
            ! >> sylvia_20200310
            ! Removing the ~30 lines below with the assumption that
            ! LRTM_FULL_VECTORIZATION is not defined.
            ! >> sylvia_20200306
            ! As for the downward radiative transfer section, here we adjust the IF
            ! statement to be n_cloudpoints(lev) == 1 or 0. There is no frational cloudiness.
            IF (n_cloudpoints(lev) == 1) THEN ! all points are cloudy

                gassrc = bbugas(lev) * atrans(lev)
                IF (.NOT. lcldlyr(lev-1)) THEN
                    cldradu = cldfrac(lev) * radlu
                    clrradu = radlu - cldradu
                    oldcld = cldradu
                    oldclr = clrradu
                    rad = 0._wp
                ENDIF
                ttot = 1._wp - atot(lev)
                cldsrc = bbutot(lev) * atot(lev)
                cldradu = cldradu * ttot + cldfrac(lev) * cldsrc
                clrradu = clrradu * (1.0_wp-atrans(lev))+(1._wp-cldfrac(lev))*gassrc
                ! Total sky radiance
                radlu = cldradu + clrradu
                urad(lev) = urad(lev) + radlu
                radmod = rad * &
                    & (facclr1(lev+1)*(1.0_wp-atrans(lev))+ &
                    & faccld1(lev+1) *  ttot) - &
                    & faccmb1(lev+1) * gassrc + &
                    & faccmb2(lev+1) * cldsrc
                oldcld = cldradu - radmod
                oldclr = clrradu + radmod
                rad = -radmod + facclr2(lev+1)*oldclr - faccld2(lev+1)*oldcld
                cldradu = cldradu + rad
                clrradu = clrradu - rad
                IF (idrv .EQ. 1) THEN
                    d_radlu_dt = d_radlu_dt * cldfrac(lev) * (1.0_wp - atot(lev)) + &
                        & d_radlu_dt * (1.0_wp - cldfrac(lev)) * (1.0_wp - atrans(lev))
                    d_urad_dt(lev) = d_urad_dt(lev) + d_radlu_dt
                ENDIF

            ELSE 
                ! >> sylvia_20200306
                ! Would correspond to IF (n_cloudpoints(lev) == 0) THEN ! all points are clear

                radlu = radlu + (bbugas(lev)-radlu)*atrans(lev)
                urad(lev) = urad(lev) + radlu
                IF (idrv .EQ. 1) THEN
                    d_radlu_dt = d_radlu_dt * (1.0_wp - atrans(lev))
                    d_urad_dt(lev) = d_urad_dt(lev) + d_radlu_dt
                ENDIF

            ENDIF

            ! >> sylvia_20200310
            !#endif
            !  Set clear sky stream to total sky stream as long as all layers
            !  are clear (iclddn=true).  Streams must be calculated separately at
            !  all layers when a cloud is present (iclddn=false), because surface
            !  reflectance is different for each stream.
            !DIR$ SIMD
            radclru = MERGE(radclru + (bbugas(lev)-radclru)*atrans(lev), &
                radlu, iclddn)
            clrurad(lev) = MERGE(clrurad(lev) + radclru, urad(lev), &
                iclddn)

            IF (idrv .EQ. 1) THEN
                !DIR$ SIMD
                d_radlu_dt = d_radlu_dt * MERGE(cldfrac(lev) * (1.0_wp - atot(lev)) + &
                    & (1.0_wp - cldfrac(lev)) * (1.0_wp - atrans(lev)), &
                    & (1.0_wp - atrans(lev)), lcldlyr(lev))
                d_urad_dt(lev) = d_urad_dt(lev) + d_radlu_dt
                d_radclru_dt = MERGE(d_radclru_dt * (1.0_wp - atrans(lev)), &
                    d_radlu_dt, iclddn)
                d_clrurad_dt(lev) = MERGE(d_clrurad_dt(lev) + d_radclru_dt, &
                    d_urad_dt(lev), iclddn)
            ENDIF
          ENDDO

          ! Increment g-point counter
          igc = igc + 1
          ! Return to continue radiative transfer for all g-channels in present band
          IF (igc .LE. ngs(iband)) GO TO 1000

          ! Process longwave output from band.
          ! Calculate upward, downward, and net flux.
          DO lev = nlayers, 0, -1
            !DIR$ SIMD
            uflux = urad(lev)*wtdiff
            dflux = drad(lev)*wtdiff
            urad(lev) = 0.0_wp
            drad(lev) = 0.0_wp
            totuflux(lev) = totuflux(lev) + uflux * delwave(iband)
            totdflux(lev) = totdflux(lev) + dflux * delwave(iband)
            uclfl = clrurad(lev)*wtdiff
            dclfl = clrdrad(lev)*wtdiff
            clrurad(lev) = 0.0_wp
            clrdrad(lev) = 0.0_wp
            totuclfl(lev) = totuclfl(lev) + uclfl * delwave(iband)
            totdclfl(lev) = totdclfl(lev) + dclfl * delwave(iband)
          ENDDO

          ! Calculate total change in upward flux wrt surface temperature
          IF (idrv .EQ. 1) THEN
            DO lev = nlayers, 0, -1
              !DIR$ SIMD
              duflux_dt = d_urad_dt(lev) * wtdiff
              d_urad_dt(lev) = 0.0_wp
              dtotuflux_dt(lev) = dtotuflux_dt(lev) + duflux_dt * delwave(iband) * fluxfac

              duclfl_dt = d_clrurad_dt(lev) * wtdiff
              d_clrurad_dt(lev) = 0.0_wp
              dtotuclfl_dt(lev) = dtotuclfl_dt(lev) + duclfl_dt * delwave(iband) * fluxfac
            ENDDO
          ENDIF

        ! End spectral band loop
        ENDDO

        ! Calculate fluxes at surface (lev==0) and model levels
        !PREVENT_INCONSISTENT_IFORT_FMA
        DO lev = 0, nlayers
          !DIR$ SIMD
          totuflux(lev) = totuflux(lev) * fluxfac
          totdflux(lev) = totdflux(lev) * fluxfac
          fnet(lev) = totuflux(lev) - totdflux(lev)
          totuclfl(lev) = totuclfl(lev) * fluxfac
          totdclfl(lev) = totdclfl(lev) * fluxfac
          fnetc(lev) = totuclfl(lev) - totdclfl(lev)
        ENDDO

        ! end vectorized version
  END SUBROUTINE lrtm_rtrnmr

  SUBROUTINE cloud_overlap(nlayers, ofs, cldfrac, &
       lcldlyr, start_lev, end_lev, lev_incr, &
       faccld1, faccld2, facclr1, facclr2, faccmb1, faccmb2)
    INTEGER, INTENT(in) :: nlayers         ! total number of layers
    INTEGER, INTENT(in) :: ofs             ! layer offset to use for
                                           ! references to facc*
    ! >> sylvia_20200318
    ! Comment the below out.
    !INTEGER, INTENT(in) :: n_cloudpoints(nlayers)
    ! >> sylvia_20200306
    ! Removing this input
    !INTEGER, INTENT(in) :: icld_ind(nlayers)
    ! << sylvia_20200306
    INTEGER, INTENT(in) :: start_lev, end_lev, lev_incr
    LOGICAL, INTENT(in) :: lcldlyr(0:nlayers+1)

    REAL(wp), INTENT(inout) :: faccld1(nlayers+1), &
         faccld2(nlayers+1), &
         facclr1(nlayers+1), &
         facclr2(nlayers+1), &
         faccmb1(nlayers+1), &
         faccmb2(nlayers+1)
    REAL(wp), INTENT(in) :: cldfrac(:)       ! layer cloud fraction
    ! >> sylvia_20200318
    ! Assuming that HAVE_FC_ATTRIBUTE_CONTINGUOUS is not defined.
!#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
!    CONTIGUOUS :: cldfrac
!#endif

    REAL(wp) :: fmax, fmin, rat1, rat2
    INTEGER :: lev, blev, clev, olev ! icld, iclr

    ! >> sylvia_20200306
    ! Removed the first dimension here.
    faccld1(start_lev-lev_incr+ofs) = 0.0_wp
    faccld2(start_lev-lev_incr+ofs) = 0.0_wp
    facclr1(start_lev-lev_incr+ofs) = 0.0_wp
    facclr2(start_lev-lev_incr+ofs) = 0.0_wp
    faccmb1(start_lev-lev_incr+ofs) = 0.0_wp
    faccmb2(start_lev-lev_incr+ofs) = 0.0_wp
    ! << sylvia_20200306

    DO lev = start_lev, end_lev - lev_incr, lev_incr
      olev = lev + ofs
      blev = lev + 1 - ofs
      clev = lev + lev_incr
      ! >> sylvia_20200306
      ! Removing the IF statement here (determining whether all points are cloudy or not).
      ! The interior of the two statements was identical, just differently iterating.
      IF (cldfrac(clev) .GE. cldfrac(lev)) THEN
        faccld1(olev) = 0._wp
        faccld2(olev) = 0._wp
        IF (.NOT. lcldlyr(lev-lev_incr)) THEN
           facclr1(olev) = 0._wp
           facclr2(olev) = 0._wp
           IF (cldfrac(lev) .LT. 1._wp) facclr2(olev) = &
              & (cldfrac(clev)-cldfrac(lev))/(1._wp-cldfrac(lev))
           facclr2(blev) = 0._wp
           faccld2(blev) = 0._wp
        ELSE
           fmax = MAX(cldfrac(lev),cldfrac(lev-lev_incr))
           IF (cldfrac(clev) .GT. fmax) THEN
             facclr1(olev) = rat2
             facclr2(olev) = (cldfrac(clev)-fmax)/(1._wp-fmax)
           ELSEIF (cldfrac(clev) .LT. fmax) THEN
             facclr1(olev) = (cldfrac(clev)-cldfrac(lev))/ &
               & (cldfrac(lev-lev_incr)-cldfrac(lev))
             facclr2(olev) = 0._wp
           ELSE
             facclr1(olev) = rat2
             facclr2(olev) = 0._wp
           ENDIF
        ENDIF
        rat1 = MERGE(1._wp, 0._wp, &
            facclr1(olev) > 0._wp .OR. facclr2(olev) > 0._wp)
        rat2 = 0._wp
      ELSE
        facclr1(olev) = 0._wp
        facclr2(olev) = 0._wp
        IF (.NOT. lcldlyr(lev-lev_incr)) THEN
           faccld1(olev) = 0._wp
           faccld2(olev) = (cldfrac(lev)-cldfrac(clev))/cldfrac(lev)

           facclr2(blev) = 0._wp
           faccld2(blev) = 0._wp
        ELSE
           fmin = MIN(cldfrac(lev),cldfrac(lev-lev_incr))
           IF (cldfrac(clev) .LE. fmin) THEN
             faccld1(olev) = rat1
             faccld2(olev) = (fmin-cldfrac(clev))/fmin
           ELSE
             faccld1(olev) = (cldfrac(lev)-cldfrac(clev))/(cldfrac(lev)-fmin)
             faccld2(olev) = 0._wp
           ENDIF
        ENDIF
        rat2 = MERGE(1._wp, 0._wp, &
             faccld1(olev) > 0._wp .OR. faccld2(olev) > 0._wp)
        rat1 = 0._wp
      ENDIF
      IF (lev == start_lev) THEN
         faccmb1(olev) = 0._wp
         faccmb2(olev) = faccld1(olev) * facclr2(blev)
      ELSE
         faccmb1(olev) = facclr1(olev) * faccld2(blev) * cldfrac(lev-lev_incr)
         faccmb2(olev) = faccld1(olev) * facclr2(blev) * (1._wp - cldfrac(lev-lev_incr))
      ENDIF
    ENDDO
    
    olev = end_lev + ofs
    faccld1(olev) = 0.0!_wp
    faccld2(olev) = 0.0!_wp
    facclr1(olev) = 0.0!_wp
    facclr2(olev) = 0.0!_wp
    faccmb1(olev) = 0.0!_wp
    faccmb2(olev) = 0.0!_wp

  END SUBROUTINE cloud_overlap

END MODULE mo_lrtm_rtrnmr
