!>
!! @brief optical depth calculations for RRTM Shortwave Radiation
!!
!! @par Description
!!     Compute the optical depth by interpolating in ln(pressure),
!!     temperature, and appropriate species.  Below LAYTROP, the water
!!     vapor self-continuum is interpolated (in temperature) separately.
!!
!! @author Eli J. Mlawer, Atmospheric & Environmental Research.
!!
!! $ID: n/a$
!!
!! @par Revisions
!!  M.Hamrud      01-Oct-2003 CY28 Cleaning
!!  JJMorcrette 2003-02-24 adapted to ECMWF environment
!!  D.Salmond  31-Oct-2007 Vector version
!!  B. Stevens (MPI) 04-10-2009 Rewritten as module for ECHAM
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!
! >> sylvia_20200313
! I don't know what this means.
! << sylvia_20200313

MODULE mo_srtm_taumol

  ! >> sylvia_20200311
  ! Pulling the following constants from mo_srtm_config
  ! >> sylvia_20200318
  ! Putting the USE statement back in
  USE mo_srtm_config,  ONLY : &
     &  jpg , ng16, ng17, ng18, ng19, ng20, ng21, ng22, &
     &  ng23, ng24, ng25, ng26, ng27, ng28, ng29, nspa, nspb
  ! << sylvia_20200311

  IMPLICIT NONE

  ! >> sylvia_20200311
  ! Replacing USE mo_kind,  ONLY : wp,i4 with the following.
  ! >> sylvia_20200311
  ! Set up the wp type that is often used below.
  ! Pull the nbndlw and ngptlw constants from mo_lrtm_par.f90
  INTEGER, PARAMETER :: pd = 12
  INTEGER, PARAMETER :: rdd = 307
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rdd) !< double precision
  INTEGER, PARAMETER :: wp = dp
  INTEGER, PARAMETER :: pi4 =  9
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)   !< at least 4 byte integer
  ! << sylvia_20200311
  
 PRIVATE
 PUBLIC :: srtm_taumol16, srtm_taumol17, srtm_taumol18, srtm_taumol19, &
       &    srtm_taumol20, srtm_taumol21, srtm_taumol22, srtm_taumol23, &
       &    srtm_taumol24, srtm_taumol25, srtm_taumol26, srtm_taumol27, &
       &    srtm_taumol28, srtm_taumol29


CONTAINS

  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 16:  2600-3250 cm-1 (low - H2O,CH4; high - CH4)
  !!
  ! >> sylvia_20200311
  ! Deleting icount, kbdim as inputs to the subroutine and in their declarations below
  ! for all subroutines.
  ! << sylvia_20200311
  SUBROUTINE srtm_taumol16 &
       & ( klev,&
       & p_fac00   , p_fac01   , p_fac10   , p_fac11,&
       & k_jp      , k_jt      , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colch4  , p_colmol,&
       & k_laytrop , p_selffac , p_selffrac, k_indself , p_forfac  , &
       & p_forfrac , k_indfor  , p_sfluxzen, p_taug    , p_taur &
       & )

    ! >> sylvia_20200313
    ! Adding the constants from MODULE yoesrta16 in mo_srtm_kgs.f90
    ! >> sylvia_20200319, Replace with the USE statement.
    USE yoesrta16, ONLY : absa, absb, forrefc, selfrefc &
      & , sfluxrefc, rayl, layreffr, strrat1
    !INTEGER, PARAMETER :: ng16 = 16
    !REAL(wp) :: absa(585,ng16), absb(235,ng16), forrefc(3,ng16)
    !REAL(wp) :: selfrefc(10,ng16), sfluxrefc(ng16)
    !REAL(wp) :: rayl, strrat1
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colch4(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    ! >> sylvia_20200313
    !  Remove the iplon declaration.
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max


    ! >> sylvia_20200313
    ! These lines with min and max do not make sense anymore as k_laytrop is a scalar.
    laytrop_min = k_laytrop        ! MINVAL
    laytrop_max = k_laytrop        ! MAXVAL
    ! << sylvia_20200313

    i_nlayers = klev
    i_laysolfr = i_nlayers  ! >> sylvia_20200313, Remove (:) from i_laysolfr.

    ! >> sylvia_20200313
    ! Remove the loop over icount.
    DO i_lay = 1, laytrop_min
         z_speccomb = p_colh2o(i_lay) + strrat1*p_colch4(i_lay)
         z_specparm = p_colh2o(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(16) + js
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(16) + js
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG16
         DO ig = 1, ng16
           p_taug(i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                & ) +                                                       &
                & p_colh2o(i_lay) *                                   &
                & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO
    ! << sylvia_20200313

    ! >> sylvia_20200313
    ! Remove the loop over icount.
    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            z_speccomb = p_colh2o(i_lay) + strrat1*p_colch4(i_lay)
            z_specparm = p_colh2o(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(16) + js
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(16) + js
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG16
            DO ig = 1, ng16
              p_taug(i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                   & ) +                                                       &
                   & p_colh2o(i_lay) *                                   &
                   & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(i_lay-1) < layreffr &
                 &  .AND. k_jp(i_lay) >= layreffr)  i_laysolfr = i_lay
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(16) + 1
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(16)+ 1
            z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG16
            DO ig = 1, ng16
              p_taug(i_lay,ig) = p_colch4(i_lay) * &
                   & (p_fac00(i_lay) * absb(ind0  ,ig) + &
                   & p_fac10(i_lay) * absb(ind0+1,ig)  + &
                   & p_fac01(i_lay) * absb(ind1  ,ig)  + &
                   & p_fac11(i_lay) * absb(ind1+1,ig))
              IF (i_lay == i_laysolfr) THEN
                  p_sfluxzen(ig) = sfluxrefc(ig)
              ENDIF
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO
    ! << sylvia_20200313

!     write(0,*) "laytrop_max+1,i_nlayers=",laytrop_max+1,i_nlayers
    ! >> sylvia_20200313
    ! Remove the loop over icount.
    DO i_lay = laytrop_max+1, i_nlayers
         IF (k_jp(i_lay-1) < layreffr &
              &  .AND. k_jp(i_lay) >= layreffr)  i_laysolfr = i_lay
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(16) + 1
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(16)+ 1
         z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG16
         DO ig = 1, ng16
            p_taug(i_lay,ig) = p_colch4(i_lay) * &
                & (p_fac00(i_lay) * absb(ind0  ,ig) + &
                & p_fac10(i_lay) * absb(ind0+1,ig)  + &
                & p_fac01(i_lay) * absb(ind1  ,ig)  + &
                & p_fac11(i_lay) * absb(ind1+1,ig))
           IF (i_lay == i_laysolfr) THEN 
               p_sfluxzen(ig) = sfluxrefc(ig)
               
           ENDIF
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

!     write(0,*) "srtm_taumol16 ends!"
  END SUBROUTINE srtm_taumol16
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)
  !
  SUBROUTINE srtm_taumol17 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colco2 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, &
       & p_forfrac , k_indfor , p_sfluxzen, p_taug   , p_taur &
       & )
     
      
    ! >> sylvia_20200313
    ! Adding the constants from MODULE yoesrta17 in mo_srtm_kgs.f90
    ! >> sylvia_20200318, Replacing with the USE statement again.
    USE yoesrta17, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, rayl , layreffr, strrat
    !INTEGER, PARAMETER :: ng17 = 16
    !REAL(wp) :: absa(585,ng17), absb(1175,ng17)
    !REAL(wp) :: forrefc(4,ng17), selfrefc(10,ng17)
    !REAL(wp) :: sfluxrefc(ng17,5), rayl, strrat
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colco2(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    ! >> sylvia_20200313
    ! Removing iplon variable.
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    ! >> sylvia_20200313
    ! Remove min and max operators.
    laytrop_min = k_laytrop  ! MINVAL
    laytrop_max = k_laytrop  ! MAXVAL
    ! << sylvia_20200313
    
    i_nlayers = klev
    i_laysolfr = i_nlayers

    ! >> sylvia_20200313
    ! Remove the loop over icount.
    DO i_lay = 1, laytrop_min
         z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
         z_specparm = p_colh2o(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(17) + js
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(17) + js
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG17
         DO ig = 1, ng17
           p_taug(i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                & )  +                                                      &
                & p_colh2o(i_lay) *                                   &
                & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))

           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
            z_specparm = p_colh2o(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(17) + js
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(17) + js
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG17
            DO ig = 1, ng17
              p_taug(i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                   & )  +                                                      &
                   & p_colh2o(i_lay) *                                   &
                   & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))

              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(i_lay-1) < layreffr &
                 & .AND. k_jp(i_lay) >= layreffr) i_laysolfr = i_lay
            z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
            z_specparm = p_colh2o(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 4._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(17)+ js
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(17)+js
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG17
            DO ig = 1, ng17
              p_taug(i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absb(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absb(ind0+5,ig) * p_fac10(i_lay) +  &
                   &                 absb(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absb(ind1+5,ig) * p_fac11(i_lay))+  &
                   & z_fs        * ( absb(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absb(ind0+6,ig) * p_fac10(i_lay) +  &
                   &                 absb(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absb(ind1+6,ig) * p_fac11(i_lay) )  &
                   & ) +                                                       &
                   & p_colh2o(i_lay) *                                   &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig)))
              IF (i_lay == i_laysolfr) p_sfluxzen(ig) = sfluxrefc(ig,js) &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
         IF (k_jp(i_lay-1) < layreffr &
              & .AND. k_jp(i_lay) >= layreffr) i_laysolfr = i_lay
         z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
         z_specparm = p_colh2o(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 4._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(17)+ js
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(17)+js
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG17
         DO ig = 1, ng17
           p_taug(i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absb(ind0,ig) * p_fac00(i_lay) +    &
                &                 absb(ind0+5,ig) * p_fac10(i_lay) +  &
                &                 absb(ind1,ig) * p_fac01(i_lay) +    &
                &                 absb(ind1+5,ig) * p_fac11(i_lay))+  &
                & z_fs        * ( absb(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absb(ind0+6,ig) * p_fac10(i_lay) +  &
                &                 absb(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absb(ind1+6,ig) * p_fac11(i_lay) )  &
                & ) +                                                       &
                & p_colh2o(i_lay) *                                   &
                & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig)))
           IF (i_lay == i_laysolfr) p_sfluxzen(ig) = sfluxrefc(ig,js) &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol17
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)
  !
  SUBROUTINE srtm_taumol18 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colch4 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing with the USE statement again.
    USE yoesrta18, ONLY : absa, absb, forrefc, selfrefc &
        & , sfluxrefc, rayl , layreffr, strrat
    ! Adding the constants from MODULE yoesrta18 in mo_srtm_kgs.f90
    !INTEGER, PARAMETER :: ng18 = 16
    !REAL(wp) :: absa(585,ng18), absb(235,ng18)
    !REAL(wp) :: forrefc(3,ng18), selfrefc(10,ng18)
    !REAL(wp) :: sfluxrefc(ng18,9), rayl, strrat
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colch4(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    ! >> sylvia_20200313
    ! Last set of comments, below is the same. No iplon.
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    ! >> sylvia_20200313
    ! Last set of comments, below is the same. Removing min and max.
    laytrop_min = k_laytrop       ! MINVAL
    laytrop_max = k_laytrop       ! MAXVAL
    ! << sylvia_20200313

    i_nlayers = klev
    i_laysolfr = k_laytrop

    ! >> sylvia_20200313
    ! Last set of comments, below is the same. Removing the iplon loop.
    DO i_lay = 1, laytrop_min
         IF (k_jp(i_lay) < layreffr                            &
              &    .AND. k_jp(i_lay+1) >= layreffr)            &
              &    i_laysolfr = MIN(i_lay+1,k_laytrop)
         z_speccomb = p_colh2o(i_lay) + strrat*p_colch4(i_lay)
         z_specparm = p_colh2o(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(18) + js
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(18) + js
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG18
         DO ig = 1, ng18
           p_taug(i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                & ) +                                                       &
                & p_colh2o(i_lay) *                                   &
                & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           IF (i_lay == i_laysolfr)                           &
                &    p_sfluxzen(ig) = sfluxrefc(ig,js)         &
                &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO
    ! << sylvia_20200313

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            IF (k_jp(i_lay) < layreffr                            &
                 &    .AND. k_jp(i_lay+1) >= layreffr)            &
                 &    i_laysolfr = MIN(i_lay+1,k_laytrop)
            z_speccomb = p_colh2o(i_lay) + strrat*p_colch4(i_lay)
            z_specparm = p_colh2o(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(18) + js
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(18) + js
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG18
            DO ig = 1, ng18
              p_taug(i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                   & ) +                                                       &
                   & p_colh2o(i_lay) *                                   &
                   & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              IF (i_lay == i_laysolfr)                           &
                   &    p_sfluxzen(ig) = sfluxrefc(ig,js)         &
                   &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(18) + 1
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(18)+ 1
            z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG18
            DO ig = 1, ng18
              p_taug(i_lay,ig) = p_colch4(i_lay) * &
                   & (p_fac00(i_lay) * absb(ind0,ig)  +  &
                   & p_fac10(i_lay) * absb(ind0+1,ig) +  &
                   & p_fac01(i_lay) * absb(ind1,ig)   +  &
                   & p_fac11(i_lay) * absb(ind1+1,ig))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(18) + 1
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(18)+ 1
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG18
         DO ig = 1, ng18
           p_taug(i_lay,ig) = p_colch4(i_lay) * &
                & (p_fac00(i_lay) * absb(ind0,ig)  +  &
                & p_fac10(i_lay) * absb(ind0+1,ig) +  &
                & p_fac01(i_lay) * absb(ind1,ig)   +  &
                & p_fac11(i_lay) * absb(ind1+1,ig))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol18
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)
  !
  SUBROUTINE srtm_taumol19 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colco2 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! Adding the constants from MODULE yoesrta19 in mo_srtm_kgs.f90
    ! >> sylvia_20200319, Replacing with the USE statement again.
    USE yoesrta19, ONLY : absa, absb, forrefc, selfrefc &
       & , sfluxrefc, rayl, layreffr, strrat
    !INTEGER, PARAMETER :: ng19 = 16
    !REAL(wp) :: absa(585,ng19), absb(235,ng19)
    !REAL(wp) :: forrefc(3,ng19), selfrefc(10,ng19)
    !REAL(wp) :: sfluxrefc(ng19,9), rayl, strrat
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colco2(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr
    REAL(wp) ::  z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop    ! MINVAL
    laytrop_max = k_laytrop    ! MAXVAL

    i_nlayers = klev
    i_laysolfr = k_laytrop

    DO i_lay = 1, laytrop_min
         IF (k_jp(i_lay) < layreffr                         &
              & .AND. k_jp(i_lay+1) >= layreffr)            &
              & i_laysolfr = MIN(i_lay+1,k_laytrop)
         z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
         z_specparm = p_colh2o(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(19) + js
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(19) + js
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG19
         DO ig = 1 , ng19
           p_taug(i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                & ) +                                                       &
                & p_colh2o(i_lay) *                                   &
                & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           IF (i_lay == i_laysolfr)                           &
                &    p_sfluxzen(ig) = sfluxrefc(ig,js)         &
                &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            IF (k_jp(i_lay) < layreffr                         &
                 & .AND. k_jp(i_lay+1) >= layreffr)            &
                 & i_laysolfr = MIN(i_lay+1,k_laytrop)
            z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
            z_specparm = p_colh2o(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(19) + js
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(19) + js
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG19
            DO ig = 1 , ng19
              p_taug(i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                   & ) +                                                       &
                   & p_colh2o(i_lay) *                                   &
                   & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              IF (i_lay == i_laysolfr)                           &
                   &    p_sfluxzen(ig) = sfluxrefc(ig,js)         &
                   &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(19) + 1
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(19)+ 1
            z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG19
            DO ig = 1 , ng19
              p_taug(i_lay,ig) = p_colco2(i_lay) * &
                   & (p_fac00(i_lay) * absb(ind0,ig) +   &
                   & p_fac10(i_lay) * absb(ind0+1,ig) +  &
                   & p_fac01(i_lay) * absb(ind1,ig) +    &
                   & p_fac11(i_lay) * absb(ind1+1,ig))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(19) + 1
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(19)+ 1
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG19
         DO ig = 1 , ng19
           p_taug(i_lay,ig) = p_colco2(i_lay) * &
                & (p_fac00(i_lay) * absb(ind0,ig) +   &
                & p_fac10(i_lay) * absb(ind0+1,ig) +  &
                & p_fac01(i_lay) * absb(ind1,ig) +    &
                & p_fac11(i_lay) * absb(ind1+1,ig))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol19
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)
  !
  SUBROUTINE srtm_taumol20 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , &
       & p_colh2o  , p_colch4 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing with the USE statement again.
    USE yoesrta20, ONLY : absa, absb, forrefc, selfrefc &
       & , sfluxrefc, absch4c, rayl, layreffr
    !INTEGER, PARAMETER :: ng20 = 16
    !REAL(wp) :: absa(65,ng20), absb(235,ng20)
    !REAL(wp) :: absch4c(ng20)
    !REAL(wp) :: forrefc(4,ng20), selfrefc(10,ng20)
    !REAL(wp) :: sfluxrefc(ng20), rayl
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colch4(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, inds, indf, i_lay, i_nlayers, i_laysolfr
    REAL(wp) :: z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop         ! MINVAL
    laytrop_max = k_laytrop         ! MAXVAL

    i_nlayers = klev
    i_laysolfr = k_laytrop

    DO i_lay = 1, laytrop_min
         IF (k_jp(i_lay) < layreffr                           &
              &    .AND. k_jp(i_lay+1) >= layreffr)           &
              &    i_laysolfr = MIN(i_lay+1,k_laytrop)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(20) + 1
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(20) + 1
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG20
         DO ig = 1 , ng20
           p_taug(i_lay,ig) = p_colh2o(i_lay) *      &
                & ((p_fac00(i_lay) * absa(ind0,ig) +       &
                & p_fac10(i_lay) * absa(ind0+1,ig) +       &
                & p_fac01(i_lay) * absa(ind1,ig) +         &
                & p_fac11(i_lay) * absa(ind1+1,ig)) +      &
                & p_selffac(i_lay) * (selfrefc(inds,ig) +  &
                & p_selffrac(i_lay) *                      &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +   &
                & p_forfac(i_lay) * (forrefc(indf,ig) +    &
                & p_forfrac(i_lay) *                       &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))      &
                & + p_colch4(i_lay) * absch4c(ig)
           p_taur(i_lay,ig) = z_tauray
           IF(i_lay == i_laysolfr) p_sfluxzen(ig)=sfluxrefc(ig)
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            IF (k_jp(i_lay) < layreffr                           &
                 &    .AND. k_jp(i_lay+1) >= layreffr)           &
                 &    i_laysolfr = MIN(i_lay+1,k_laytrop)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(20) + 1
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(20) + 1
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG20
            DO ig = 1 , ng20
              p_taug(i_lay,ig) = p_colh2o(i_lay) *      &
                   & ((p_fac00(i_lay) * absa(ind0,ig) +       &
                   & p_fac10(i_lay) * absa(ind0+1,ig) +       &
                   & p_fac01(i_lay) * absa(ind1,ig) +         &
                   & p_fac11(i_lay) * absa(ind1+1,ig)) +      &
                   & p_selffac(i_lay) * (selfrefc(inds,ig) +  &
                   & p_selffrac(i_lay) *                      &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +   &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +    &
                   & p_forfrac(i_lay) *                       &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))      &
                   & + p_colch4(i_lay) * absch4c(ig)
              p_taur(i_lay,ig) = z_tauray
              IF(i_lay == i_laysolfr) p_sfluxzen(ig)=sfluxrefc(ig)
            ENDDO
          ELSE
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(20) + 1
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(20) +1
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG20
            DO ig = 1 , ng20
              p_taug(i_lay,ig) = p_colh2o(i_lay) *   &
                   & (p_fac00(i_lay) * absb(ind0,ig) +     &
                   & p_fac10(i_lay) * absb(ind0+1,ig) +    &
                   & p_fac01(i_lay) * absb(ind1,ig) +      &
                   & p_fac11(i_lay) * absb(ind1+1,ig) +    &
                   & p_forfac(i_lay) * (forrefc(indf,ig) + &
                   & p_forfrac(i_lay) *                    &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig)))) + &
                   & p_colch4(i_lay) * absch4c(ig)
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(20) + 1
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(20) +1
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG20
         DO ig = 1 , ng20
           p_taug(i_lay,ig) = p_colh2o(i_lay) *   &
                & (p_fac00(i_lay) * absb(ind0,ig) +     &
                & p_fac10(i_lay) * absb(ind0+1,ig) +    &
                & p_fac01(i_lay) * absb(ind1,ig) +      &
                & p_fac11(i_lay) * absb(ind1+1,ig) +    &
                & p_forfac(i_lay) * (forrefc(indf,ig) + &
                & p_forfrac(i_lay) *                    &
                & (forrefc(indf+1,ig) - forrefc(indf,ig)))) + &
                & p_colch4(i_lay) * absch4c(ig)
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol20
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
  !
  SUBROUTINE srtm_taumol21 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colco2 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing with USE statement again.
    USE yoesrta21, ONLY : absa, absb, forrefc, selfrefc &
      & , sfluxrefc, rayl, layreffr, strrat
    !INTEGER, PARAMETER :: ng21 = 16
    !REAL(wp) :: absa(585,ng21), absb(1175,ng21)
    !REAL(wp) :: forrefc(4,ng21), selfrefc(10,ng21)
    !REAL(wp) :: sfluxrefc(ng21,9), rayl, strrat
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colco2(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm,  z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop
    laytrop_max = k_laytrop

    i_nlayers = klev
    i_laysolfr = k_laytrop

    DO i_lay = 1, laytrop_min
         IF (k_jp(i_lay) < layreffr                           &
              &    .AND. k_jp(i_lay+1) >= layreffr)           &
              &    i_laysolfr = MIN(i_lay+1,k_laytrop)
         z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
         z_specparm = p_colh2o(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(21) + js
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(21) + js
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG21
         DO ig = 1 , ng21
           p_taug(i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                & ) +                                                       &
                & p_colh2o(i_lay) *                                   &
                & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           IF (i_lay == i_laysolfr)                        &
                &    p_sfluxzen(ig) = sfluxrefc(ig,js)      &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            IF (k_jp(i_lay) < layreffr                           &
                 &    .AND. k_jp(i_lay+1) >= layreffr)           &
                 &    i_laysolfr = MIN(i_lay+1,k_laytrop)
            z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
            z_specparm = p_colh2o(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(21) + js
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(21) + js
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG21
            DO ig = 1 , ng21
              p_taug(i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                   & ) +                                                       &
                   & p_colh2o(i_lay) *                                   &
                   & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              IF (i_lay == i_laysolfr)                        &
                   &    p_sfluxzen(ig) = sfluxrefc(ig,js)      &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
            z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
            z_specparm = p_colh2o(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 4._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(21) +js
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(21)+js
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG21
            DO ig = 1 , ng21
              p_taug(i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp - z_fs) * ( absb(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absb(ind0+5,ig) * p_fac10(i_lay) +  &
                   &                 absb(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absb(ind1+5,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absb(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absb(ind0+6,ig) * p_fac10(i_lay) +  &
                   &                 absb(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absb(ind1+6,ig) * p_fac11(i_lay) )  &
                   & ) +                                                       &
                   & p_colh2o(i_lay) *                                   &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig)))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
         z_speccomb = p_colh2o(i_lay) + strrat*p_colco2(i_lay)
         z_specparm = p_colh2o(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 4._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(21) +js
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(21)+js
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG21
         DO ig = 1 , ng21
           p_taug(i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp - z_fs) * ( absb(ind0,ig) * p_fac00(i_lay) +    &
                &                 absb(ind0+5,ig) * p_fac10(i_lay) +  &
                &                 absb(ind1,ig) * p_fac01(i_lay) +    &
                &                 absb(ind1+5,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absb(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absb(ind0+6,ig) * p_fac10(i_lay) +  &
                &                 absb(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absb(ind1+6,ig) * p_fac11(i_lay) )  &
                & ) +                                                       &
                & p_colh2o(i_lay) *                                   &
                & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig)))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol21
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 29:  7700-8050 cm-1 (low - H2O,O2; high - O2)
  !
  SUBROUTINE srtm_taumol22 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colmol , p_colo2,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing with USE statement again.
    USE yoesrta22, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, rayl, layreffr, strrat
    !INTEGER, PARAMETER :: ng22 = 16
    !REAL(wp) :: absa(585,ng22), absb(235,ng22)
    !REAL(wp) :: forrefc(3,ng22), selfrefc(10,ng22)
    !REAL(wp) :: sfluxrefc(ng22,9), rayl, strrat
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_colo2(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray, z_o2adj , &
         &      z_o2cont
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop
    laytrop_max = k_laytrop

    i_nlayers = klev

    !     The following factor is the ratio of total O2 band intensity (lines
    !     and Mate continuum) to O2 band intensity (line only).  It is needed
    !     to adjust the optical depths since the k's include only lines.
    z_o2adj = 1.6_wp
    i_laysolfr = k_laytrop

    DO i_lay = 1, laytrop_min
         IF (k_jp(i_lay) < layreffr                           &
              &    .AND. k_jp(i_lay+1) >= layreffr)           &
              &    i_laysolfr = MIN(i_lay+1,k_laytrop)
         z_o2cont = 4.35e-4_wp*p_colo2(i_lay)/(350.0_wp*2.0_wp)
         z_speccomb = p_colh2o(i_lay) +          &
              &     z_o2adj*strrat*p_colo2(i_lay)
         z_specparm = p_colh2o(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(22) + js
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(22) + js
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG22
         DO ig = 1 , ng22
           p_taug(i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                & ) +                                                       &
                & p_colh2o(i_lay) *                                   &
                & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))                 &
                & + z_o2cont
           IF (i_lay == i_laysolfr)                        &
                & p_sfluxzen(ig) = sfluxrefc(ig,js)         &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            IF (k_jp(i_lay) < layreffr                           &
                 &    .AND. k_jp(i_lay+1) >= layreffr)           &
                 &    i_laysolfr = MIN(i_lay+1,k_laytrop)
            z_o2cont = 4.35e-4_wp*p_colo2(i_lay)/(350.0_wp*2.0_wp)
            z_speccomb = p_colh2o(i_lay) +          &
                 &     z_o2adj*strrat*p_colo2(i_lay)
            z_specparm = p_colh2o(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(22) + js
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(22) + js
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG22
            DO ig = 1 , ng22
              p_taug(i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                   & ) +                                                       &
                   & p_colh2o(i_lay) *                                   &
                   & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))                 &
                   & + z_o2cont
              IF (i_lay == i_laysolfr)                        &
                   & p_sfluxzen(ig) = sfluxrefc(ig,js)         &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
            z_o2cont = 4.35e-4_wp*p_colo2(i_lay)/(350.0_wp*2.0_wp)
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(22) + 1
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(22)+ 1
            z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG22
            DO ig = 1 , ng22
              p_taug(i_lay,ig) = p_colo2(i_lay) * z_o2adj * &
                   & (p_fac00(i_lay) * absb(ind0,ig) +            &
                   & p_fac10(i_lay) * absb(ind0+1,ig) +           &
                   & p_fac01(i_lay) * absb(ind1,ig) +             &
                   & p_fac11(i_lay) * absb(ind1+1,ig)) +          &
                   & z_o2cont
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
         z_o2cont = 4.35e-4_wp*p_colo2(i_lay)/(350.0_wp*2.0_wp)
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(22) + 1
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(22)+ 1
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG22
         DO ig = 1 , ng22
           p_taug(i_lay,ig) = p_colo2(i_lay) * z_o2adj * &
                & (p_fac00(i_lay) * absb(ind0,ig) +            &
                & p_fac10(i_lay) * absb(ind0+1,ig) +           &
                & p_fac01(i_lay) * absb(ind1,ig) +             &
                & p_fac11(i_lay) * absb(ind1+1,ig)) +          &
                & z_o2cont
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol22
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
  !
  SUBROUTINE srtm_taumol23 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , &
       & p_colh2o  , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing with USE statement again.
    USE yoesrta23, ONLY : absa, forrefc, selfrefc &
         & , sfluxrefc, raylc, layreffr, givfac
    !INTEGER, PARAMETER :: ng23 = 16
    !REAL(wp) :: absa(65,ng23)
    !REAL(wp) :: forrefc(3,ng23), selfrefc(10,ng23)
    !REAL(wp) :: sfluxrefc(ng23), raylc(ng23), givfac
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, inds, indf, i_lay, i_nlayers, i_laysolfr
    REAL(wp) ::  z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop
    laytrop_max = k_laytrop

    i_nlayers = klev
    i_laysolfr = k_laytrop

    DO i_lay = 1, laytrop_min
         IF (k_jp(i_lay) < layreffr                            &
              &    .AND. k_jp(i_lay+1) >= layreffr)            &
              &    i_laysolfr = MIN(i_lay+1,k_laytrop)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(23) + 1
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(23) + 1
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)

!CDIR EXPAND=NG23
         DO ig = 1 , ng23
           z_tauray = p_colmol(i_lay) * raylc(ig)
           p_taug(i_lay,ig) = p_colh2o(i_lay) *         &
                & (givfac * (p_fac00(i_lay) * absa(ind0,ig) + &
                & p_fac10(i_lay) * absa(ind0+1,ig) +          &
                & p_fac01(i_lay) * absa(ind1,ig) +            &
                & p_fac11(i_lay) * absa(ind1+1,ig)) +         &
                & p_selffac(i_lay) * (selfrefc(inds,ig) +     &
                & p_selffrac(i_lay) *                         &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +      &
                & p_forfac(i_lay) * (forrefc(indf,ig) +       &
                & p_forfrac(i_lay) *                          &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           IF (i_lay == i_laysolfr) &
                p_sfluxzen(ig) = sfluxrefc(ig)
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            IF (k_jp(i_lay) < layreffr                            &
                 &    .AND. k_jp(i_lay+1) >= layreffr)            &
                 &    i_laysolfr = MIN(i_lay+1,k_laytrop)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(23) + 1
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(23) + 1
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)

!CDIR EXPAND=NG23
            DO ig = 1 , ng23
              z_tauray = p_colmol(i_lay) * raylc(ig)
              p_taug(i_lay,ig) = p_colh2o(i_lay) *         &
                   & (givfac * (p_fac00(i_lay) * absa(ind0,ig) + &
                   & p_fac10(i_lay) * absa(ind0+1,ig) +          &
                   & p_fac01(i_lay) * absa(ind1,ig) +            &
                   & p_fac11(i_lay) * absa(ind1+1,ig)) +         &
                   & p_selffac(i_lay) * (selfrefc(inds,ig) +     &
                   & p_selffrac(i_lay) *                         &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +      &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +       &
                   & p_forfrac(i_lay) *                          &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              IF (i_lay == i_laysolfr) &
                   p_sfluxzen(ig) = sfluxrefc(ig)
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
!CDIR EXPAND=NG23
            DO ig = 1 , ng23
              p_taug(i_lay,ig) = 0.0_wp
              p_taur(i_lay,ig) = p_colmol(i_lay) * raylc(ig)
            ENDDO
          ENDIF
    ENDDO

    DO ig = 1 , ng23
      DO i_lay = laytrop_max+1, i_nlayers
          p_taug(i_lay,ig) = 0.0_wp
          p_taur(i_lay,ig) = p_colmol(i_lay) * raylc(ig)
      ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol23
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 24:  12850-16000 cm-1 (low - H2O,O2; high - O2)
  !
  SUBROUTINE srtm_taumol24 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colmol , p_colo2   , p_colo3,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing the USE statement again.
    USE yoesrta24, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, abso3ac, abso3bc, raylac, raylbc, layreffr, strrat
    !INTEGER, PARAMETER :: ng24 = 16
    !REAL(wp) :: absa(585,ng24), absb(235,ng24)
    !REAL(wp) :: forrefc(3,ng24), selfrefc(10,ng24)
    !REAL(wp) :: sfluxrefc(ng24,9), strrat
    !REAL(wp) :: raylac(ng24,9), raylbc(ng24)
    !REAL(wp) :: abso3ac(ng24), abso3bc(ng24)
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_colo2(klev)
    REAL(wp)   ,INTENT(in)    :: p_colo3(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr
    REAL(wp) ::  z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop
    laytrop_max = k_laytrop

    i_nlayers = klev
    i_laysolfr = k_laytrop

    DO i_lay = 1, laytrop_min
         IF (k_jp(i_lay) < layreffr                            &
              &    .AND. k_jp(i_lay+1) >= layreffr)            &
              &    i_laysolfr = MIN(i_lay+1,k_laytrop)
         z_speccomb = p_colh2o(i_lay) + strrat*p_colo2(i_lay)
         z_specparm = p_colh2o(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(24) + js
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(24) + js
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)

!CDIR EXPAND=NG24
         DO ig = 1 , ng24
           z_tauray = p_colmol(i_lay) * (raylac(ig,js) + &
                & z_fs * (raylac(ig,js+1) - raylac(ig,js)))
           p_taug(i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                & ) +                                                       &
                & p_colo3(i_lay) * abso3ac(ig) +                      &
                & p_colh2o(i_lay) *                                   &
                & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           IF (i_lay == i_laysolfr)                        &
                & p_sfluxzen(ig) = sfluxrefc(ig,js)         &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(i_lay,ig) =  z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            IF (k_jp(i_lay) < layreffr                            &
                 &    .AND. k_jp(i_lay+1) >= layreffr)            &
                 &    i_laysolfr = MIN(i_lay+1,k_laytrop)
            z_speccomb = p_colh2o(i_lay) + strrat*p_colo2(i_lay)
            z_specparm = p_colh2o(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(24) + js
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(24) + js
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)

!CDIR EXPAND=NG24
            DO ig = 1 , ng24
              z_tauray = p_colmol(i_lay) * (raylac(ig,js) + &
                   & z_fs * (raylac(ig,js+1) - raylac(ig,js)))
              p_taug(i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                   & ) +                                                       &
                   & p_colo3(i_lay) * abso3ac(ig) +                      &
                   & p_colh2o(i_lay) *                                   &
                   & (p_selffac(i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              IF (i_lay == i_laysolfr)                        &
                   & p_sfluxzen(ig) = sfluxrefc(ig,js)         &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(i_lay,ig) =  z_tauray
            ENDDO
          ELSE
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(24) + 1
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(24)+ 1

!CDIR EXPAND=NG24
            DO ig = 1 , ng24
              z_tauray = p_colmol(i_lay) * raylbc(ig)
              p_taug(i_lay,ig) = p_colo2(i_lay) *  &
                   & (p_fac00(i_lay) * absb(ind0,ig) +   &
                   & p_fac10(i_lay) * absb(ind0+1,ig) +  &
                   & p_fac01(i_lay) * absb(ind1,ig) +    &
                   & p_fac11(i_lay) * absb(ind1+1,ig)) + &
                   & p_colo3(i_lay) * abso3bc(ig)
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(24) + 1
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(24)+ 1

!CDIR EXPAND=NG24
         DO ig = 1 , ng24
           z_tauray = p_colmol(i_lay) * raylbc(ig)
           p_taug(i_lay,ig) = p_colo2(i_lay) *  &
                & (p_fac00(i_lay) * absb(ind0,ig) +   &
                & p_fac10(i_lay) * absb(ind0+1,ig) +  &
                & p_fac01(i_lay) * absb(ind1,ig) +    &
                & p_fac11(i_lay) * absb(ind1+1,ig)) + &
                & p_colo3(i_lay) * abso3bc(ig)
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol24
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 25:  16000-22650 cm-1 (low - H2O; high - nothing)
  !
  SUBROUTINE srtm_taumol25 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , &
       & p_colh2o  , p_colmol , p_colo3,&
       & k_laytrop,&
       & p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing the USE statement again.
    USE yoesrta25, ONLY : absa, sfluxrefc, abso3ac, abso3bc, raylc, layreffr    
    !INTEGER, PARAMETER :: ng25 = 16
    !REAL(wp) :: absa(65,ng21)
    !REAL(wp) :: sfluxrefc(ng25), raylc(ng25)
    !REAL(wp) :: abso3ac(ng25), abso3bc(ng25)
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_colo3(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, i_lay, i_nlayers, i_laysolfr
    REAL(wp) ::  z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop
    laytrop_max = k_laytrop

    i_nlayers = klev
    i_laysolfr = k_laytrop

    DO i_lay = 1, laytrop_min
         IF (k_jp(i_lay) < layreffr .AND.   &
              &    k_jp(i_lay+1) >= layreffr) &
              &    i_laysolfr = MIN(i_lay+1,k_laytrop)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(25) + 1
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(25) + 1
!CDIR EXPAND=NG25
         DO ig = 1 , ng25
           z_tauray = p_colmol(i_lay) * raylc(ig)
           p_taug(i_lay,ig) = p_colh2o(i_lay) * &
                & (p_fac00(i_lay) * absa(ind0,ig)   + &
                & p_fac10(i_lay) * absa(ind0+1,ig)  + &
                & p_fac01(i_lay) * absa(ind1,ig)    + &
                & p_fac11(i_lay) * absa(ind1+1,ig)) + &
                & p_colo3(i_lay) * abso3ac(ig)
           IF(i_lay == i_laysolfr) p_sfluxzen(ig)=sfluxrefc(ig)
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            IF (k_jp(i_lay) < layreffr .AND.   &
                 &    k_jp(i_lay+1) >= layreffr) &
                 &    i_laysolfr = MIN(i_lay+1,k_laytrop)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(25) + 1
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(25) + 1
!CDIR EXPAND=NG25
            DO ig = 1 , ng25
              z_tauray = p_colmol(i_lay) * raylc(ig)
              p_taug(i_lay,ig) = p_colh2o(i_lay) * &
                   & (p_fac00(i_lay) * absa(ind0,ig)   + &
                   & p_fac10(i_lay) * absa(ind0+1,ig)  + &
                   & p_fac01(i_lay) * absa(ind1,ig)    + &
                   & p_fac11(i_lay) * absa(ind1+1,ig)) + &
                   & p_colo3(i_lay) * abso3ac(ig)
              IF(i_lay == i_laysolfr) p_sfluxzen(ig)=sfluxrefc(ig)
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
!CDIR EXPAND=NG25
            DO ig = 1 , ng25
              z_tauray = p_colmol(i_lay) * raylc(ig)
              p_taug(i_lay,ig) = p_colo3(i_lay) * abso3bc(ig)
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO

    DO ig = 1 , ng25
      DO i_lay = laytrop_max+1, i_nlayers
          z_tauray = p_colmol(i_lay) * raylc(ig)
          p_taug(i_lay,ig) = p_colo3(i_lay) * abso3bc(ig)
          p_taur(i_lay,ig) = z_tauray
      ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol25
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)
  !
  SUBROUTINE srtm_taumol26 &
       & ( klev,&
       & p_colmol,&
       & k_laytrop , &
       & p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing the USE statement.
    USE yoesrta26, ONLY : sfluxrefc, raylc
    !INTEGER, PARAMETER :: ng26 = 16
    !REAL(wp) :: sfluxrefc(ng26), raylc(ng26)
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_laytrop
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, i_lay, i_nlayers, i_laysolfr
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop
    laytrop_max = k_laytrop

    i_nlayers = klev
    i_laysolfr = k_laytrop

    DO i_lay = 1, laytrop_min
!CDIR EXPAND=NG26
         DO ig = 1 , ng26
           IF(i_lay == i_laysolfr) p_sfluxzen(ig)=sfluxrefc(ig)
           p_taug(i_lay,ig) = 0.0_wp
           p_taur(i_lay,ig) = p_colmol(i_lay) * raylc(ig)
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
!CDIR EXPAND=NG26
            DO ig = 1 , ng26
              IF(i_lay == i_laysolfr) p_sfluxzen(ig)=sfluxrefc(ig)
              p_taug(i_lay,ig) = 0.0_wp
              p_taur(i_lay,ig) = p_colmol(i_lay) * raylc(ig)
            ENDDO
          ELSE
!CDIR EXPAND=NG26
            DO ig = 1 , ng26
              p_taug(i_lay,ig) = 0.0_wp
              p_taur(i_lay,ig) = p_colmol(i_lay) * raylc(ig)
            ENDDO
          ENDIF
    ENDDO

    DO ig = 1 , ng26
       DO i_lay = laytrop_max+1, i_nlayers
           p_taug(i_lay,ig) = 0.0_wp
           p_taur(i_lay,ig) = p_colmol(i_lay) * raylc(ig)
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol26
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 27:  29000-38000 cm-1 (low - O3; high - O3)
  !
  SUBROUTINE srtm_taumol27 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_colmol  , p_colo3,&
       & k_laytrop, p_sfluxzen, p_taug    , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing the USE statement again.
    USE yoesrta27, ONLY : absa, absb, sfluxrefc, raylc, layreffr, scalekur    
    !INTEGER, PARAMETER :: ng27 = 16
    !REAL(wp) :: absa(65,ng27), absb(235,ng27)
    !REAL(wp) :: sfluxrefc(ng27), raylc(ng27)
    !REAL(wp) :: scalekur
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_colo3(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, i_lay, i_nlayers, i_laysolfr
    REAL(wp) :: z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop
    laytrop_max = k_laytrop

    i_nlayers = klev
    i_laysolfr = i_nlayers

    DO i_lay = 1, laytrop_min
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(27) + 1
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(27) + 1
!CDIR EXPAND=NG27
         DO ig = 1 , ng27
           z_tauray = p_colmol(i_lay) * raylc(ig)
           p_taug(i_lay,ig) = p_colo3(i_lay) * &
                & (p_fac00(i_lay) * absa(ind0,ig)  + &
                & p_fac10(i_lay) * absa(ind0+1,ig) + &
                & p_fac01(i_lay) * absa(ind1,ig) +   &
                & p_fac11(i_lay) * absa(ind1+1,ig))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(27) + 1
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(27) + 1
!CDIR EXPAND=NG27
            DO ig = 1 , ng27
              z_tauray = p_colmol(i_lay) * raylc(ig)
              p_taug(i_lay,ig) = p_colo3(i_lay) * &
                   & (p_fac00(i_lay) * absa(ind0,ig)  + &
                   & p_fac10(i_lay) * absa(ind0+1,ig) + &
                   & p_fac01(i_lay) * absa(ind1,ig) +   &
                   & p_fac11(i_lay) * absa(ind1+1,ig))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(i_lay-1) < layreffr &
                 &    .AND. k_jp(i_lay) >= layreffr) i_laysolfr = i_lay
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(27) + 1
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(27)+ 1
!CDIR EXPAND=NG27
            DO ig = 1 , ng27
              z_tauray = p_colmol(i_lay) * raylc(ig)
              p_taug(i_lay,ig) = p_colo3(i_lay) * &
                   & (p_fac00(i_lay) * absb(ind0,ig)  + &
                   & p_fac10(i_lay) * absb(ind0+1,ig) + &
                   & p_fac01(i_lay) * absb(ind1,ig) +   &
                   & p_fac11(i_lay) * absb(ind1+1,ig))
              IF (i_lay == i_laysolfr) &
                   & p_sfluxzen(ig) = scalekur * sfluxrefc(ig)
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
         IF (k_jp(i_lay-1) < layreffr &
              &    .AND. k_jp(i_lay) >= layreffr) i_laysolfr = i_lay
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(27) + 1
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(27)+ 1
!CDIR EXPAND=NG27
         DO ig = 1 , ng27
           z_tauray = p_colmol(i_lay) * raylc(ig)
           p_taug(i_lay,ig) = p_colo3(i_lay) * &
                & (p_fac00(i_lay) * absb(ind0,ig)  + &
                & p_fac10(i_lay) * absb(ind0+1,ig) + &
                & p_fac01(i_lay) * absb(ind1,ig) +   &
                & p_fac11(i_lay) * absb(ind1+1,ig))
           IF (i_lay == i_laysolfr) &
                & p_sfluxzen(ig) = scalekur * sfluxrefc(ig)
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol27
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 28:  38000-50000 cm-1 (low - O3,O2; high - O3,O2)
  !
  SUBROUTINE srtm_taumol28 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colmol  , p_colo2  , p_colo3,&
       & k_laytrop,&
       & p_sfluxzen, p_taug   , p_taur &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing the USE statement below.
    USE yoesrta28, ONLY : absa, absb, sfluxrefc, rayl, layreffr, strrat    
    !INTEGER, PARAMETER :: ng28 = 16
    !REAL(wp) :: absa(65,ng28), absb(1175,ng28)
    !REAL(wp) :: sfluxrefc(ng28,5), rayl, strrat
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_colo2(klev)
    REAL(wp)   ,INTENT(in)    :: p_colo3(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, js, i_lay, i_nlayers, i_laysolfr
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop
    laytrop_max = k_laytrop

    i_nlayers = klev
    i_laysolfr = i_nlayers

    DO i_lay = 1, laytrop_min
         z_speccomb = p_colo3(i_lay) + strrat*p_colo2(i_lay)
         z_specparm = p_colo3(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(28) + js
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(28) + js
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG28
         DO ig = 1 , ng28
           p_taug(i_lay,ig) = z_speccomb * &
                & (&
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                & )
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            z_speccomb = p_colo3(i_lay) + strrat*p_colo2(i_lay)
            z_specparm = p_colo3(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(28) + js
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(28) + js
            z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG28
            DO ig = 1 , ng28
              p_taug(i_lay,ig) = z_speccomb * &
                   & (&
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(i_lay) ) &
                   & )
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(i_lay-1) < layreffr &
                 &    .AND. k_jp(i_lay) >= layreffr) i_laysolfr = i_lay
            z_speccomb = p_colo3(i_lay) + strrat*p_colo2(i_lay)
            z_specparm = p_colo3(i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus,z_specparm)
            z_specmult = 4._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(28)+ js
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(28)+js
            z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG28
            DO ig = 1 , ng28
              p_taug(i_lay,ig) = z_speccomb * &
                   & (&
                   & (1._wp- z_fs) * ( absb(ind0,ig) * p_fac00(i_lay) +    &
                   &                 absb(ind0+5,ig) * p_fac10(i_lay) +  &
                   &                 absb(ind1,ig) * p_fac01(i_lay) +    &
                   &                 absb(ind1+5,ig) * p_fac11(i_lay) )+ &
                   & z_fs        * ( absb(ind0+1,ig) * p_fac00(i_lay) +  &
                   &                 absb(ind0+6,ig) * p_fac10(i_lay) +  &
                   &                 absb(ind1+1,ig) * p_fac01(i_lay) +  &
                   &                 absb(ind1+6,ig) * p_fac11(i_lay) )  &
                   & )
              IF (i_lay == i_laysolfr) p_sfluxzen(ig) = sfluxrefc(ig,js) &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO


    DO i_lay = laytrop_max+1, i_nlayers
         IF (k_jp(i_lay-1) < layreffr &
              &    .AND. k_jp(i_lay) >= layreffr) i_laysolfr = i_lay
         z_speccomb = p_colo3(i_lay) + strrat*p_colo2(i_lay)
         z_specparm = p_colo3(i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus,z_specparm)
         z_specmult = 4._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(28)+ js
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(28)+js
         z_tauray = p_colmol(i_lay) * rayl

!CDIR EXPAND=NG28
         DO ig = 1 , ng28
           p_taug(i_lay,ig) = z_speccomb * &
                & (&
                & (1._wp- z_fs) * ( absb(ind0,ig) * p_fac00(i_lay) +    &
                &                 absb(ind0+5,ig) * p_fac10(i_lay) +  &
                &                 absb(ind1,ig) * p_fac01(i_lay) +    &
                &                 absb(ind1+5,ig) * p_fac11(i_lay) )+ &
                & z_fs        * ( absb(ind0+1,ig) * p_fac00(i_lay) +  &
                &                 absb(ind0+6,ig) * p_fac10(i_lay) +  &
                &                 absb(ind1+1,ig) * p_fac01(i_lay) +  &
                &                 absb(ind1+6,ig) * p_fac11(i_lay) )  &
                & )
           IF (i_lay == i_laysolfr) p_sfluxzen(ig) = sfluxrefc(ig,js) &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol28
  !-----------------------------------------------------------------------------
  !>
  !! @brief BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)
  !
  SUBROUTINE srtm_taumol29 &
       & ( klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     ,&
       & p_colh2o  , p_colco2 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, &
       & k_indfor, p_sfluxzen, p_taug   , p_taur  &
       & )

    ! >> sylvia_20200313
    ! >> sylvia_20200319, Replacing the USE statement below.
    USE yoesrta29, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, absh2oc, absco2c, rayl, layreffr
    !INTEGER, PARAMETER :: ng29 = 16
    !REAL(wp) :: absa(65,ng29), absb(235,ng29)
    !REAL(wp) :: forrefc(4,ng29), selfrefc(10,ng29)
    !REAL(wp) :: sfluxrefc(ng29), absh2oc(ng29), absco2c(ng29)
    !REAL(wp) :: rayl
    !INTEGER  :: layreffr
    ! << sylvia_20200313

    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(klev)
    INTEGER,INTENT(in)    :: k_jt(klev)
    INTEGER,INTENT(in)    :: k_jt1(klev)
    INTEGER,INTENT(in)    :: k_laytrop
    INTEGER,INTENT(in)    :: k_indself(klev)
    INTEGER,INTENT(in)    :: k_indfor(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(klev)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(klev)
    REAL(wp)   ,INTENT(in)    :: p_colco2(klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(klev,jpg)

    INTEGER :: ig, ind0, ind1, inds, indf, i_lay, i_nlayers
    INTEGER :: i_laysolfr
    REAL(wp) ::  z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = k_laytrop
    laytrop_max = k_laytrop

    i_nlayers = klev
    i_laysolfr = i_nlayers

    DO i_lay = 1, laytrop_min
         ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(29) + 1
         ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(29) + 1
         inds = k_indself(i_lay)
         indf = k_indfor(i_lay)
         z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG29
         DO ig = 1, ng29
           p_taug(i_lay,ig) = p_colh2o(i_lay) *     &
                & ((p_fac00(i_lay) * absa(ind0,ig) +      &
                & p_fac10(i_lay) * absa(ind0+1,ig) +      &
                & p_fac01(i_lay) * absa(ind1,ig) +        &
                & p_fac11(i_lay) * absa(ind1+1,ig)) +     &
                & p_selffac(i_lay) * (selfrefc(inds,ig) + &
                & p_selffrac(i_lay) *                     &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +  &
                & p_forfac(i_lay) * (forrefc(indf,ig) +   &
                & p_forfrac(i_lay) *                      &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))     &
                & + p_colco2(i_lay) * absco2c(ig)
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
          IF (i_lay <= k_laytrop) THEN
            ind0 = ((k_jp(i_lay)-1)*5+(k_jt(i_lay)-1))*nspa(29) + 1
            ind1 = (k_jp(i_lay)*5+(k_jt1(i_lay)-1))*nspa(29) + 1
            inds = k_indself(i_lay)
            indf = k_indfor(i_lay)
            z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG29
            DO ig = 1, ng29
              p_taug(i_lay,ig) = p_colh2o(i_lay) *     &
                   & ((p_fac00(i_lay) * absa(ind0,ig) +      &
                   & p_fac10(i_lay) * absa(ind0+1,ig) +      &
                   & p_fac01(i_lay) * absa(ind1,ig) +        &
                   & p_fac11(i_lay) * absa(ind1+1,ig)) +     &
                   & p_selffac(i_lay) * (selfrefc(inds,ig) + &
                   & p_selffrac(i_lay) *                     &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +  &
                   & p_forfac(i_lay) * (forrefc(indf,ig) +   &
                   & p_forfrac(i_lay) *                      &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))     &
                   & + p_colco2(i_lay) * absco2c(ig)
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(i_lay-1) < layreffr                                &
                 &   .AND. k_jp(i_lay) >= layreffr)  i_laysolfr = i_lay
            ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(29) + 1
            ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(29)+ 1
            z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG29
            DO ig = 1 , ng29
              p_taug(i_lay,ig) = p_colco2(i_lay) * &
                   & (p_fac00(i_lay) * absb(ind0,ig) +   &
                   & p_fac10(i_lay) * absb(ind0+1,ig) +  &
                   & p_fac01(i_lay) * absb(ind1,ig) +    &
                   & p_fac11(i_lay) * absb(ind1+1,ig))   &
                   & + p_colh2o(i_lay) * absh2oc(ig)
              IF (i_lay == i_laysolfr) p_sfluxzen(ig) = sfluxrefc(ig)
              p_taur(i_lay,ig) = z_tauray
            ENDDO
          ENDIF
    ENDDO



    DO i_lay = laytrop_max+1, i_nlayers
         IF (k_jp(i_lay-1) < layreffr                                &
              &   .AND. k_jp(i_lay) >= layreffr)  i_laysolfr = i_lay
         ind0 = ((k_jp(i_lay)-13)*5+(k_jt(i_lay)-1))*nspb(29) + 1
         ind1 = ((k_jp(i_lay)-12)*5+(k_jt1(i_lay)-1))*nspb(29)+ 1
         z_tauray = p_colmol(i_lay) * rayl
!CDIR EXPAND=NG29
         DO ig = 1 , ng29
           p_taug(i_lay,ig) = p_colco2(i_lay) * &
                & (p_fac00(i_lay) * absb(ind0,ig) +   &
                & p_fac10(i_lay) * absb(ind0+1,ig) +  &
                & p_fac01(i_lay) * absb(ind1,ig) +    &
                & p_fac11(i_lay) * absb(ind1+1,ig))   &
                & + p_colh2o(i_lay) * absh2oc(ig)
           IF (i_lay == i_laysolfr) p_sfluxzen(ig) = sfluxrefc(ig)
           p_taur(i_lay,ig) = z_tauray
         ENDDO
    ENDDO

    !-----------------------------------------------------------------------
  END SUBROUTINE srtm_taumol29

END MODULE mo_srtm_taumol
