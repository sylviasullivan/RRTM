!>
!! Calculates GEMS ozone climatology.
!! Taken and adapted from ECMWF's IFS (37r2).
!!
!! @par Revision History
!! Initial Release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-10-18)
!!
! >> sylvia_20200305, Changing the input from mtime_datetime to day, month, hour, minute.
MODULE o3_gems

! >> sylvia_20200318
! Make this back into a module
USE o3_gems_data,   ONLY: rghg7
IMPLICIT NONE
PRIVATE
PUBLIC :: calc_o3_gems
  
CONTAINS
! << sylvia_20200318

  SUBROUTINE calc_o3_gems(month, day, hour, minute, latt, klev, pp_hl, dpres_mc, xm_o3)

    ! >> sylvia_20200309
    ! Set up the wp type that is often used below.
    ! >> sylvia_20200313
    INTEGER, PARAMETER :: pd =  12
    INTEGER, PARAMETER :: rd = 307
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
    INTEGER, PARAMETER :: wp = dp
    REAL(wp), PARAMETER :: pi = 3.14159265358979323846264338327950288_wp
    REAL(wp), PARAMETER ::  rad2deg = 180.0_wp/pi
    REAL(wp), PARAMETER :: deg2rad = pi/180.0_wp
    ! << sylvia_20200309
    
    INTEGER, INTENT(in)  ::        &
      &  klev,                     & !< number of levels
      &  month,                    & !< month associated with this profile, sylvia_20200309
      &  day,                      & !< day associated with this profile, sylvia_20200305
      &  hour,                     & !< hour associated with this profile
      &  minute                      !< minute associated with this profile
    ! >> sylvia_20200318
    REAL(wp), INTENT(in) ::        &
      & latt                         ! >> sylvia_20200318, Input latitude.
    REAL(wp)  ::                   & ! >> sylvia_20200310, Ideally I would want INTENT(in) here...
      &  dpres_mc(klev),           & !< layer thickness with respect to pressure
      &  pp_hl(klev+1)               !< half level pressure in Pa

    ! >> sylvia_20200305
    ! Changing the TYPE(t_external_data), INTENT(inout) :: ext_data to just REAL, INTENT(out) :: xm_o3
    REAL(wp), INTENT(out) ::       &
      &  xm_o3(klev)
    ! << sylvia_20200305

    CHARACTER(len=*), PARAMETER :: routine = 'calc_o3_gems'
    INTEGER , PARAMETER :: ilat = 64, nlev_gems = 91
    

    
    
    ! Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).
    REAL(wp), PARAMETER     :: zrefp(nlev_gems) = (/ &
      & 0.10000000E+01_wp, 0.29900000E+01_wp, 0.56834998E+01_wp, 0.10147500E+02_wp, &
      & 0.17160500E+02_wp, 0.27682501E+02_wp, 0.42848999E+02_wp, 0.63956501E+02_wp, &
      & 0.92440987E+02_wp, 0.12985049E+03_wp, 0.17781149E+03_wp, 0.23799651E+03_wp, &
      & 0.31209000E+03_wp, 0.40175449E+03_wp, 0.50860199E+03_wp, 0.63416602E+03_wp, &
      & 0.77987903E+03_wp, 0.94705548E+03_wp, 0.11368750E+04_wp, 0.13503740E+04_wp, &
      & 0.15884360E+04_wp, 0.18517920E+04_wp, 0.21410149E+04_wp, 0.24565259E+04_wp, &
      & 0.27986001E+04_wp, 0.31673640E+04_wp, 0.35628110E+04_wp, 0.39848059E+04_wp, &
      & 0.44330962E+04_wp, 0.49073169E+04_wp, 0.54070068E+04_wp, 0.59314971E+04_wp, &
      & 0.64797832E+04_wp, 0.70505981E+04_wp, 0.76428970E+04_wp, 0.82572246E+04_wp, &
      & 0.88957646E+04_wp, 0.95614326E+04_wp, 0.10257570E+05_wp, 0.10988080E+05_wp, &
      & 0.11757620E+05_wp, 0.12571580E+05_wp, 0.13435160E+05_wp, 0.14352070E+05_wp, &
      & 0.15325200E+05_wp, 0.16357520E+05_wp, 0.17452211E+05_wp, 0.18612539E+05_wp, &
      & 0.19841900E+05_wp, 0.21144000E+05_wp, 0.22522650E+05_wp, 0.23981760E+05_wp, &
      & 0.25525529E+05_wp, 0.27158301E+05_wp, 0.28884641E+05_wp, 0.30709301E+05_wp, &
      & 0.32637240E+05_wp, 0.34673641E+05_wp, 0.36823859E+05_wp, 0.39093570E+05_wp, &
      & 0.41487949E+05_wp, 0.44010520E+05_wp, 0.46651148E+05_wp, 0.49386160E+05_wp, &
      & 0.52190051E+05_wp, 0.55035488E+05_wp, 0.57894629E+05_wp, 0.60745699E+05_wp, &
      & 0.63574230E+05_wp, 0.66367703E+05_wp, 0.69114523E+05_wp, 0.71801031E+05_wp, &
      & 0.74412922E+05_wp, 0.76938641E+05_wp, 0.79368109E+05_wp, 0.81689688E+05_wp, &
      & 0.83891563E+05_wp, 0.85965023E+05_wp, 0.87903594E+05_wp, 0.89700203E+05_wp, &
      & 0.91347422E+05_wp, 0.92841422E+05_wp, 0.94181656E+05_wp, 0.95367977E+05_wp, &
      & 0.96399797E+05_wp, 0.97280203E+05_wp, 0.98034609E+05_wp, 0.98677672E+05_wp, &
      & 0.99198242E+05_wp, 0.99594977E+05_wp, 0.99881500E+05_wp /)

    ! Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).
    REAL(wp), PARAMETER     :: zmday(12) = (/ &
      &  31._wp,    59.25_wp,  90.25_wp, 120.25_wp, 151.25_wp, 181.25_wp, &
      & 212.25_wp, 243.25_wp, 273.25_wp, 304.25_wp, 334.25_wp, 365.25_wp /)

    ! Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).
    REAL(wp), PARAMETER     :: zytime(12) = (/ &
      &  22320._wp,  64980._wp, 107640._wp, 151560._wp, 195480._wp, 239400._wp, &
      & 283320._wp, 327960._wp, 371880._wp, 415800._wp, 459720._wp, 503640._wp /)    

    ! local fields
    ! >> sylvia_20200309
    ! Below syntax like idx(0:klev) or zlat(0:ilat+1) caused errors.
    ! Rewritten as simply idx(klev) or zlat(ilat+1)
    ! >> sylvia_20200313
    ! Removing unused dtdz,deltaz,o3_clim variables below. Replacing l_found by lfound.
    INTEGER  :: idx0(0:klev)
    REAL(wp) :: zlat(0:ilat+1)
    REAL(wp) :: zozn(0:ilat+1,1:nlev_gems)
    REAL(wp) :: zpresh(0:nlev_gems)
    REAL(wp) :: rclpr(0:nlev_gems)
    REAL(wp) :: zo3(1:nlev_gems)
    REAL(wp) :: zviozo(0:klev)
    REAL(wp) :: zozovi(0:nlev_gems) 
    ! >> sylvia_20200324, changing klev to nlev_gems for zozovi

    ! local scalars
    ! >> sylvia_20200310, Replaced LOGICAL :: lfound_all below with LOGICAL :: lfound
    ! >> sylvia_20200313, Removing jk1,zadd_o3 below.
    INTEGER  :: jk,jkk,jl !loop indices
    INTEGER  :: idy,im,imn,im1,im2,jk_start
    REAL(wp) :: ztimi,zxtime,zjl,zlatint,zint
    LOGICAL  :: lfound
    
    ! >> sylvia_20200305
    ! Copying some parameters from rrtm_interface
    REAL(wp), PARAMETER   :: amd       = 28.970_wp        !> [g/mol] dry air
    REAL(wp), PARAMETER   :: amo3      = 47.9982_wp       !! [g/mol] O3
    ! << sylvia_20200305
    
    ! >> sylvia_20200305
    ! Changing mtime_datetime%date%day and mtime_datetime%date%month to day, month.
    IDY = day - 1      ! NDD(KINDAT)-1
    IMN = month        ! NMM(KINDAT)
    ! << sylvia_20200305
    
    IF (IMN == 1) THEN
      ! >> sylvia_20200305
      ! Changing mtime_datetime%date%hour and mtime_datetime%date%minute to hour, minute.
      ZXTIME=REAL(IDY*1440 + hour*60 + minute, wp)      
      ! << sylvia_20200305
    ELSEIF (IMN == 2) THEN
      IF(IDY == 28) IDY=IDY-1
      ! A DAY IN FEB. IS 28.25*24*60/28=1452.8571min LONG.
      ZXTIME=44640._wp+REAL(IDY,wp)*1452.8571_wp + REAL(hour*60 + minute, wp)
    ELSE
      ZXTIME=(ZMDAY(IMN-1)+REAL(IDY,KIND(ZXTIME)))*1440._wp + REAL(hour*60 + minute, wp)
    ENDIF
    ! 525960=MINUTES IN A SIDERAL YEAR (365.25d)
    ZXTIME=MOD(ZXTIME,525960._wp)

    IM1=0
    IM2=0
    IF (ZXTIME <= ZYTIME(1)) THEN
      IM1=12
      IM2=1
      ZTIMI=(ZYTIME(1)-ZXTIME)/44640._wp
    ELSEIF(ZXTIME > ZYTIME(12)) THEN
      IM1=12
      IM2=1
      ZTIMI=(548280._wp-ZXTIME)/44640._wp
      ! 548280.=(365.25d + 15.5d)*24*60
    ELSE
      DO IM=1,11
        IF (ZXTIME > ZYTIME(IM) .AND. ZXTIME <= ZYTIME(IM+1)) THEN
          IM1=IM
          IM2=IM+1
          ZTIMI=(ZXTIME-ZYTIME(IM2))/(ZYTIME(IM1)-ZYTIME(IM2))
        ENDIF
      ENDDO
      IF ( IM1 == 0 .OR. IM2 == 0 ) THEN
        ! >> sylvia_20200318
        ! Remove the CALL finish(TRIM(routine)) and replace with print
        WRITE(*,*) 'calc_o3_gems: Problem with time interpolation in suecozc!'
        ! << sylvia_20200318
      ENDIF
    ENDIF

    ! Pressure levels of climatological ozone fields. From su_ghgclim.F90 of ECMWF's IFS (37r2).

    ZPRESH(0)=0.0_wp
    RCLPR(0) =0.0_wp
    DO JK=1,NLEV_GEMS-1
      ZPRESH(JK)=(ZREFP(JK)+ZREFP(JK+1))*0.5_wp
      RCLPR(JK) =ZPRESH(JK)
    ENDDO
    ZPRESH(NLEV_GEMS)=110000._wp
    RCLPR(nlev_gems) =ZPRESH(nlev_gems)
    ! Preparations for latitude interpolations

    zlatint=180._wp/REAL(ilat,wp)

    DO jl=0,ilat+1
      zlat(jl) = (-90._wp+0.5_wp*zlatint+(jl-1)*zlatint)*deg2rad
    ENDDO

    ! volume mixing ratio to ozone pressure thickness

    ! >> sylvia_20200305
    ! Extracting only case 7 below as irad_o3 == 7.
    DO jk=1,nlev_gems
      DO jl=1,ilat
        zozn(JL,JK) = amo3/amd * (RGHG7(JL,JK,IM2)&
          & +ZTIMI*(RGHG7(JL,JK,IM1)-RGHG7(JL,JK,IM2)))
        zozn(JL,JK) = zozn(JL,JK) * (ZPRESH(JK)-ZPRESH(JK-1))
      ENDDO
    ENDDO

    DO jk=1,nlev_gems
      zozn(0,JK)      = zozn(1,jk)
      zozn(ilat+1,jk) = zozn(ilat,jk)
    ENDDO

    ! >> sylvia_20200305
    ! Removing i_startblk and i_endblk here, as well as the associated DO loop.
    ! Latitude interpolation
    
    ! Removing i_startidx, i_endidx here. Assume one fixed latitude for the input profile.
    DO jkk = 1, nlev_gems
      zjl = 1._wp + (rad2deg*latt+90._wp-.5_wp*zlatint)/zlatint
      !just select nearest value
      !zo3(jc,jkk,jb) = zozn(NINT(zjl),jkk)

      !linear interpolation
      zo3(jkk) = zozn(INT(zjl),jkk) &
        & + (zozn(INT(zjl) + 1,jkk) - zozn(INT(zjl),jkk)) &
        &  /  (zlat(INT(zjl)) - zlat(INT(zjl)+1)) &
        &  * (zlat(INT(zjl)) - latt)

    ENDDO !jkk

    ! ACCUMULATE FROM TOP TO BOTTOM THE LATITUDE INTERPOLATED FIELDS
    ! From radghg.F90 of ECMWF's IFS.
    zozovi(0) = 0._wp
    DO jkk = 1, nlev_gems
      zozovi(jkk) = zozovi(jkk-1) + zo3(jkk)
    ENDDO

    ! REDISTRIBUTE THE VERTIC. INTEGR. CLIM. O3 ON THE MODEL GRID
    ! Adapted from radghg.F90 of ECMWF's IFS.
    jk = 0
    zviozo(0) = 0._wp
    jk_start = 0

    DO jk = 0, klev
      lfound = .FALSE.
      ! >> sylvia_20200310, Removing the lfound_all = .FALSE. initialization.
      DO jkk = jk_start,nlev_gems-1
        IF( pp_hl(jk+1) >= RCLPR(jkk)  &
          & .AND. pp_hl(jk+1) < RCLPR(jkk+1)) THEN
          ZINT = (pp_hl(jk+1) - RCLPR(jkk))/(RCLPR(jkk+1) - RCLPR(jkk)) 
          ZVIOZO(JK) = ZOZOVI(jkk) + ZINT * (ZOZOVI(jkk+1) - ZOZOVI(jkk))
          lfound = .TRUE.
          idx0(jk) = jkk
        ELSEIF ( pp_hl(jk+1) > RCLPR(nlev_gems) ) THEN
          lfound = .TRUE.
          idx0(jk) = nlev_gems
        ENDIF

        ! >> sylvia_20200310
        ! Taking the ALL operator off of ALL(l_found) below and removing lfound_all = .TRUE.
        IF (lfound) THEN
          EXIT
        ENDIF
        ! << sylvia_20200310
      ENDDO !jkk
      ! >> sylvia_20200310
      ! Changed lfound_all to lfound below. Removed MINVAL(idx0(jk)) below as there's no array.
      IF (lfound) THEN
        jk_start = MIN(idx0(jk),nlev_gems-1)
      ENDIF
      ! << sylvia_20200310
    ENDDO !jk

    ! COMPUTE THE MASS MIXING RATIO
    DO jk = 1, klev
        xm_o3(jk) = (ZVIOZO(jk) - ZVIOZO(jk-1)) / dpres_mc(jk)
    ! Assuming that ltuning_ozone is 0 and removing the associated IF here.
    ENDDO

    ! >> sylvia_20200305
    ! Also assuming that icpl_o3_tp parameter equals 0 by default. It was not been specified in the nwp_phy_nml.
    ! ==> No ozone adaptation around the extratropical tropopause
    ! << sylvia_20200305
 
  END SUBROUTINE calc_o3_gems

END MODULE o3_gems
