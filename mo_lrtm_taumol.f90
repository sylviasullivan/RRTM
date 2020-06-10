! >> sylvia_20200318
!#include "mod1.inc" Copying and pasting the contents of this here.
! << sylvia_20200318
module mo_lrtm_taumol

  implicit none
  private
  PUBLIC  :: chi_mls, lrtm_taumol

  ! >> sylvia_20200310
  ! Set up the wp type that is often used below.
  INTEGER, PARAMETER :: pd =  12
  INTEGER, PARAMETER :: rd = 307
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
  INTEGER, PARAMETER :: wp = dp
  ! << sylvia_20200310

  ! >> sylvia_20200313
  ! Replacing actual variables with these values.
  !use mo_lrtm_par, only : nspa, nspb
  ! nspa = For the lower atmosphere, the number of reference atmospheres
  ! that are stored for each spectral band per pressure level and temperature.
  ! Each of these atmospheres has different relative amounts of the key 
  ! species for the band.
  INTEGER, PARAMETER :: jpband = 16
  INTEGER :: nspa(jpband)
  ! nspb = Same as nspa for the upper atmosphere.
  INTEGER :: nspb(jpband)
  ! << sylvia_20200313

  real(wp), parameter :: oneminus = 1.0_wp - 1.0e-06_wp
  real(wp) :: chi_mls(7,59)
  
CONTAINS

  ! >> sylvia_20200318
  ! Define the MOD1 function instead of an include file.
  function MOD1(x) 
    REAL(wp), INTENT(IN) :: x
    REAL(wp)             :: MOD1

    ! Below AINT truncates its argument to a whole number.
    ! So this function is giving the departure from the next whole number.
    MOD1 = (x) - AINT((x))
  end function MOD1
  ! << sylvia_20200318

  subroutine lrtm_taumol(nlayers, pavel, wx, coldry, &
               & laytrop, jp, jt, jt1, &
               & colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
               & colbrd, fac00, fac01, fac10, fac11, &
               & rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
               & rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
               & rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
               & selffac, selffrac, indself, forfac, forfrac, indfor, &
               & minorfrac, scaleminor, scaleminorn2, indminor, &
               & fracs, taug )
        !----------------------------------------------------------------------------

        ! *******************************************************************************
        ! *                                                                             *
        ! *                  Optical depths developed for the                           *
        ! *                                                                             *
        ! *                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
        ! *                                                                             *
        ! *                                                                             *
        ! *            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
        ! *                        131 HARTWELL AVENUE                                  *
        ! *                        LEXINGTON, MA 02421                                  *
        ! *                                                                             *
        ! *                                                                             *
        ! *                           ELI J. MLAWER                                     *
        ! *                         JENNIFER DELAMERE                                   *
        ! *                         STEVEN J. TAUBMAN                                   *
        ! *                         SHEPARD A. CLOUGH                                   *
        ! *                                                                             *
        ! *                                                                             *
        ! *                                                                             *
        ! *                                                                             *
        ! *                       email:  mlawer@aer.com                                *
        ! *                       email:  jdelamer@aer.com                              *
        ! *                                                                             *
        ! *        The authors wish to acknowledge the contributions of the             *
        ! *        following people:  Karen Cady-Pereira, Patrick D. Brown,             *
        ! *        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom.    *
        ! *                                                                             *
        ! *******************************************************************************
        ! *                                                                             *
        ! *  Revision for g-point reduction: Michael J. Iacono, AER, Inc.               *
        ! *                                                                             *
        ! *******************************************************************************
        ! *     TAUMOL                                                                  *
        ! *                                                                             *
        ! *     This file contains the subroutines TAUGBn (where n goes from            *
        ! *     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
        ! *     per g-value and layer for band n.                                       *
        ! *                                                                             *
        ! *  Output:  optical depths (unitless)                                         *
        ! *           fractions needed to compute Planck functions at every layer       *
        ! *               and g-value                                                   *
        ! *                                                                             *
        ! *     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
        ! *     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
        ! *                                                                             *
        ! *  Input                                                                      *
        ! *                                                                             *
        ! *     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
        ! *     COMMON /PRECISE/  ONEMINUS                                              *
        ! *     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
        ! *     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
        ! *     COMMON /PROFDATA/ LAYTROP,                                              *
        ! *    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
        ! *    &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),             *
        ! *    &                  COLO2(MXLAY)
        ! *     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
        ! *    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
        ! *     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
        ! *     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
        ! *                                                                             *
        ! *     Description:                                                            *
        ! *     NG(IBAND) - number of g-values in band IBAND                            *
        ! *     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
        ! *                   atmospheres that are stored for band IBAND per            *
        ! *                   pressure level and temperature.  Each of these            *
        ! *                   atmospheres has different relative amounts of the         *
        ! *                   key species for the band (i.e. different binary           *
        ! *                   species parameters).                                      *
        ! *     NSPB(IBAND) - same for upper atmosphere                                 *
        ! *     ONEMINUS - since problems are caused in some cases by interpolation     *
        ! *                parameters equal to or greater than 1, for these cases       *
        ! *                these parameters are set to this value, slightly < 1.        *
        ! *     PAVEL - layer pressures (mb)                                            *
        ! *     TAVEL - layer temperatures (degrees K)                                  *
        ! *     PZ - level pressures (mb)                                               *
        ! *     TZ - level temperatures (degrees K)                                     *
        ! *     LAYTROP - layer at which switch is made from one combination of         *
        ! *               key species to another                                        *
        ! *     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
        ! *               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
        ! *               respectively (molecules/cm**2)                                *
        ! *     FACij(LAY) - for layer LAY, these are factors that are needed to        *
        ! *                  compute the interpolation factors that multiply the        *
        ! *                  appropriate reference k-values.  A value of 0 (1) for      *
        ! *                  i,j indicates that the corresponding factor multiplies     *
        ! *                  reference k-value for the lower (higher) of the two        *
        ! *                  appropriate temperatures, and altitudes, respectively.     *
        ! *     JP - the index of the lower (in altitude) of the two appropriate        *
        ! *          reference pressure levels needed for interpolation                 *
        ! *     JT, JT1 - the indices of the lower of the two appropriate reference     *
        ! *               temperatures needed for interpolation (for pressure           *
        ! *               levels JP and JP+1, respectively)                             *
        ! *     SELFFAC - scale factor needed for water vapor self-continuum, equals    *
        ! *               (water vapor density)/(atmospheric density at 296K and        *
        ! *               1013 mb)                                                      *
        ! *     SELFFRAC - factor needed for temperature interpolation of reference     *
        ! *                water vapor self-continuum data                              *
        ! *     INDSELF - index of the lower of the two appropriate reference           *
        ! *               temperatures needed for the self-continuum interpolation      *
        ! *     FORFAC  - scale factor needed for water vapor foreign-continuum.        *
        ! *     FORFRAC - factor needed for temperature interpolation of reference      *
        ! *                water vapor foreign-continuum data                           *
        ! *     INDFOR  - index of the lower of the two appropriate reference           *
        ! *               temperatures needed for the foreign-continuum interpolation   *
        ! *                                                                             *
        ! *  Data input                                                                 *
        ! *     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),*
        ! *                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'                        *
        ! *        (note:  n is the band number,'MGAS' is the species name of the minor *
        ! *         gas)                                                                *
        ! *                                                                             *
        ! *     Description:                                                            *
        ! *     KA - k-values for low reference atmospheres (key-species only)          *
        ! *          (units: cm**2/molecule)                                            *
        ! *     KB - k-values for high reference atmospheres (key-species only)         *
        ! *          (units: cm**2/molecule)                                            *
        ! *     KA_M'MGAS' - k-values for low reference atmosphere minor species        *
        ! *          (units: cm**2/molecule)                                            *
        ! *     KB_M'MGAS' - k-values for high reference atmosphere minor species       *
        ! *          (units: cm**2/molecule)                                            *
        ! *     SELFREF - k-values for water vapor self-continuum for reference         *
        ! *               atmospheres (used below LAYTROP)                              *
        ! *               (units: cm**2/molecule)                                       *
        ! *     FORREF  - k-values for water vapor foreign-continuum for reference      *
        ! *               atmospheres (used below/above LAYTROP)                        *
        ! *               (units: cm**2/molecule)                                       *
        ! *                                                                             *
        ! *     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
        ! *     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
        ! *                                                                             *
        !*******************************************************************************

        ! ------- Declarations -------

        ! ----- Input -----
        integer, intent(in) :: nlayers         ! total number of layers
        real(wp), intent(in) :: pavel(:)     ! layer pressures (mb)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: wx(:,:)      ! cross-section amounts (mol/cm2)
        !    Dimensions: (maxxsec,nlayers)
        real(wp), intent(in) :: coldry(:)    ! column amount (dry air)
        !    Dimensions: (nlayers)

        integer, intent(in) :: laytrop        ! tropopause layer index
        integer, intent(in) :: jp(:)           !
        !    Dimensions: (nlayers)
        integer, intent(in) :: jt(:)           !
        !    Dimensions: (nlayers)
        integer, intent(in) :: jt1(:)          !
        !    Dimensions: (nlayers)

        real(wp), intent(in) :: colh2o(:)          ! column amount (h2o)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: colco2(:)          ! column amount (co2)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: colo3(:)           ! column amount (o3)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: coln2o(:)          ! column amount (n2o)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: colco(:)           ! column amount (co)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: colch4(:)          ! column amount (ch4)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: colo2(:)           ! column amount (o2)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: colbrd(:)          ! column amount (broadening gases)
        !    Dimensions: (nlayers)

        integer, intent(in) :: indself(:)
        !    Dimensions: (nlayers)
        integer, intent(in) :: indfor(:)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: selffac(:)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: selffrac(:)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: forfac(:)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: forfrac(:)
        !    Dimensions: (nlayers)

        integer, intent(in) :: indminor(:)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: minorfrac(:)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: scaleminor(:)
        !    Dimensions: (nlayers)
        real(wp), intent(in) :: scaleminorn2(:)
        !    Dimensions: (nlayers)

        real(wp), intent(in) :: &                  !
                fac00(:), fac01(:), &             !    Dimensions: (nlayers)
                fac10(:), fac11(:)
        real(wp), intent(in) :: &                  !
                rat_h2oco2(:),rat_h2oco2_1(:), &
                rat_h2oo3(:),rat_h2oo3_1(:), & !    Dimensions: (nlayers)
                rat_h2on2o(:),rat_h2on2o_1(:), &
                rat_h2och4(:),rat_h2och4_1(:), &
                rat_n2oco2(:),rat_n2oco2_1(:), &
                rat_o3co2(:),rat_o3co2_1(:)

        ! ----- Output -----
        ! >> sylvia_20200306
        ! Remove the first nproma dimension.
        real(wp), intent(out) :: fracs(:,:)        ! planck fractions
        !    Dimensions: (nlayers,ngptlw)
        real(wp), intent(out) :: taug(:,:)         ! gaseous optical depth
        !    Dimensions: (nlayers,ngptlw)
        ! << sylvia_20200306

        ! >> sylvia_20200318, Remove jc declaration.
        integer :: icl, ich, lay

        !     local integer arrays
        INTEGER :: laytrop_min, laytrop_max
        integer :: ixc(nlayers)
        ! >> sylvia_20200306
        ! Removing these variables: ixlow(nlayers), ixhigh(nlayers)

        ! >> sylvia_20200310
        ! These lines with MINVAL MAXVAL do not make sense anymore as laytrop is a scalar.
        laytrop_min = laytrop    ! MINVAL
        laytrop_max = laytrop    ! MAXVAL
        ! << sylvia_20200310
        ixc    = 0

        ! create index lists for mixed layers
        ! >> sylvia_20200306
        ! Removing the do jc = 1, kproma loop here. 
        ! Because of the above this loop should actually never run...
        do lay = laytrop_min+1, laytrop_max
        icl = 0
        ich = 0
        if ( lay <= laytrop ) then
                icl = icl + 1
                ! >> sylvia_20200306 Removing this statemnent:
                ! ixlow(icl,lay) = jc
                ! >> sylvia_20200310 Removing the (jc) index on laytrop.
        else
                ich = ich + 1
                ! >> sylvia_20200306 Removing this statemnent:
                ! ixhigh(ich,lay) = jc
        endif
        ixc(lay) = icl
        enddo
        ! << sylvia_20200306


        ! Calculate gaseous optical depth and planck fractions for each spectral band.

        call taugb1
        call taugb2
        call taugb3
        call taugb4
        call taugb5
        call taugb6
        call taugb7
        call taugb8
        call taugb9
        call taugb10
        call taugb11
        call taugb12
        call taugb13
        call taugb14
        call taugb15
        call taugb16

contains

        !----------------------------------------------------------------------------
        subroutine taugb1
                !----------------------------------------------------------------------------

                ! ------- Modifications -------
                !  Written by Eli J. Mlawer, Atmospheric & Environmental Research.
                !  Revised by Michael J. Iacono, Atmospheric & Environmental Research.
                !
                !     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
                !                          (high key - h2o; high minor - n2)
                !
                !     note: previous versions of rrtm band 1:
                !           10-250 cm-1 (low - h2o; high - h2o)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200313, commenting the three lines below out, not used.
                !REAL(wp) :: kao(5,13,no1)
                !REAL(wp) :: kbo(5,13:59,no1)
                !REAL(wp) :: kao_mn2(19,no1) , kbo_mn2(19,no1)
                !REAL(wp) :: selfrefo(10,no1), forrefo(4,no1)
                
                ! >> sylvia_20200320
                ! Putting the USE statements back in.
                use rrlw_kg01,   only : fracrefa, fracrefb, absa, absb, &
                         ka_mn2, kb_mn2, selfref, forref
                !REAL(wp) :: fracrefa(ng1)  , fracrefb(ng1)
                !REAL(wp) :: absa(65,ng1)   , absb(235,ng1)
                !REAL(wp) :: ka_mn2(19,ng1) , kb_mn2(19,ng1)
                
                ! >> sylvia_20200310, Replace use mo_lrtm_par, only : ng1
                ! with INTEGER, PARAMETER :: ng1  = 10 from the parameter
                ! module.
                ! Pulling only necessary data declarations from MODULE
                ! rrlw_kg01.
                ! fracrefa, fracrefb, absa, absb, ka_mn2, kb_mn2, selfref,
                ! forref
                INTEGER, PARAMETER :: ng1  = 10
                INTEGER, PARAMETER :: no1  = 16
                ! << sylvia_20200310

                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318
                ! Removing ixc0 and ixp declarations.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                real(wp) :: pp, corradj, scalen2, tauself, taufor, taun2


                ! Minor gas mapping levels:
                !     lower - n2, p = 142.5490 mbar, t = 215.70 k
                !     upper - n2, p = 142.5490 mbar, t = 215.70 k

                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature.  Below laytrop, the water vapor self-continuum and
                ! foreign continuum is interpolated (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min
                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(1) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(1) + 1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)
                pp = pavel(lay)
                corradj =  1.0_wp
                if (pp .lt. 250._wp) then
                    corradj = 1._wp - 0.15_wp * (250._wp-pp) / 154.4_wp
                endif

                scalen2 = colbrd(lay) * scaleminorn2(lay)
                !CDIR EXPAND=NG1
                ! >> sylvia_20200306
                ! Note these iteractions with ig are over "reduced g intervals per spectral band".
                do ig = 1, ng1
                   tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                            (selfref(inds+1,ig) - selfref(inds,ig)))
                   taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                            (forref(indf+1,ig) -  forref(indf,ig)))
                   taun2 = scalen2*(ka_mn2(indm,ig) + &
                           minorfrac(lay) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
                   taug(lay,ig) = corradj * (colh2o(lay) * &
                            (fac00(lay) * absa(ind0,ig) + &
                            fac10(lay) * absa(ind0+1,ig) + &
                            fac01(lay) * absa(ind1,ig) + &
                            fac11(lay) * absa(ind1+1,ig)) &
                            + tauself + taufor + taun2)
                    fracs(lay,ig) = fracrefa(ig)
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers
                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(1) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(1) + 1
                indf = indfor(lay)
                indm = indminor(lay)
                pp = pavel(lay)
                corradj =  1._wp - 0.15_wp * (pp / 95.6_wp)

                scalen2 = colbrd(lay) * scaleminorn2(lay)
                !CDIR EXPAND=NG1
                do ig = 1, ng1
                taufor = forfac(lay) * (forref(indf,ig) + &
                        forfrac(lay) * (forref(indf+1,ig) - forref(indf,ig)))
                taun2 = scalen2*(kb_mn2(indm,ig) + &
                        minorfrac(lay) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
                taug(lay,ig) = corradj * (colh2o(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + taufor + taun2)
                fracs(lay,ig) = fracrefb(ig)
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max
                !CDIR NODEP,VOVERTAKE,VOB
                ! >> sylvia_20200306
                ! Removing ixc0 = ixc(lay); do ixp = 1, ixc0; jl = ixlow(ixp,lay)
                ! Code that would have set the first index / horizontal dimension of the following.
                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(1) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(1) + 1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)
                pp = pavel(lay)
                corradj =  1.0_wp
                if (pp .lt. 250._wp) then
                        corradj = 1._wp - 0.15_wp * (250._wp-pp) / 154.4_wp
                endif
                ! >> sylvia_20200306

                scalen2 = colbrd(lay) * scaleminorn2(lay)
                !CDIR EXPAND=NG1

                do ig = 1, ng1
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) -  forref(indf,ig)))
                taun2 = scalen2*(ka_mn2(indm,ig) + &
                        minorfrac(lay) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
                taug(lay,ig) = corradj * (colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) + &
                        fac11(lay) * absa(ind1+1,ig)) &
                        + tauself + taufor + taun2)
                fracs(lay,ig) = fracrefa(ig)
                enddo

                ! Upper atmosphere part
                ! >> sylvia_20200306
                ! Remvoing ixc0 = kproma - ixc0;  do ixp = 1, ixc0; jl = ixhigh(ixp,lay)
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(1) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(1) + 1
                indf = indfor(lay)
                indm = indminor(lay)
                pp = pavel(lay)
                corradj =  1._wp - 0.15_wp * (pp / 95.6_wp)

                scalen2 = colbrd(lay) * scaleminorn2(lay)
                !CDIR EXPAND=NG1
                do ig = 1, ng1
                taufor = forfac(lay) * (forref(indf,ig) + &
                        forfrac(lay) * (forref(indf+1,ig) - forref(indf,ig)))
                taun2 = scalen2*(kb_mn2(indm,ig) + &
                        minorfrac(lay) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
                taug(lay,ig) = corradj * (colh2o(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + taufor + taun2)
                fracs(lay,ig) = fracrefb(ig)
                enddo

                enddo

        end subroutine taugb1

        !----------------------------------------------------------------------------
        subroutine taugb2
                !----------------------------------------------------------------------------
                !
                !     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
                !
                !     note: previous version of rrtm band 2:
                !           250 - 500 cm-1 (low - h2o; high - h2o)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Putting the USE statements back in.
                use rrlw_kg02,   only : fracrefa, fracrefb, absa, absb, &
                        selfref, forref
                !REAL(wp) :: fracrefa(ng2)  , fracrefb(ng2)
                !REAL(wp) :: absa(65,ng2), absb(235,ng2)
                !REAL(wp) :: selfref(10,ng2), forref(4,ng2)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                ! Pulling ng2, ngs1 from mo_lrtm_par.f90
                ! fracrefa, fracrefb, absa, absb, selfref, forrek from
                ! mo_lrtm_kgs.f90
                INTEGER, PARAMETER :: ng2   = 12
                INTEGER, PARAMETER :: ngs1  = 10
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf
                real(wp) :: pp, corradj, tauself, taufor


                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature.  Below laytrop, the water vapor self-continuum and
                ! foreign continuum is interpolated (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(2) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(2) + 1
                inds = indself(lay)
                indf = indfor(lay)
                pp = pavel(lay)
                corradj = 1._wp - .05_wp * (pp - 100._wp) / 900._wp
                !CDIR EXPAND=NG2
                do ig = 1, ng2
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs1+ig) = corradj * (colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) + &
                        fac11(lay) * absa(ind1+1,ig)) &
                        + tauself + taufor)
                fracs(lay,ngs1+ig) = fracrefa(ig)
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(2) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(2) + 1
                indf = indfor(lay)
                !CDIR EXPAND=NG2
                do ig = 1, ng2
                taufor =  forfac(lay) * (forref(indf,ig) + &
                        forfrac(lay) * (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs1+ig) = colh2o(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + taufor
                fracs(lay,ngs1+ig) = fracrefb(ig)
                enddo

                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max
                ! >> sylvia_20200306
                ! Removing ixc0 = ixc(lay); do ixp = 1, ixc0j; l = ixlow(ixp,lay)
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(2) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(2) + 1
                inds = indself(lay)
                indf = indfor(lay)
                pp = pavel(lay)
                corradj = 1._wp - .05_wp * (pp - 100._wp) / 900._wp
                !CDIR EXPAND=NG2
                do ig = 1, ng2
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs1+ig) = corradj * (colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) + &
                        fac11(lay) * absa(ind1+1,ig)) &
                        + tauself + taufor)
                fracs(lay,ngs1+ig) = fracrefa(ig)
                enddo

                ! Upper atmosphere part
                ! >> sylvia_20200306
                ! Removing ixc0 = kproma - ixc0; do ixp = 1, ixc0; jl = ixhigh(ixp,lay)
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(2) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(2) + 1
                indf = indfor(lay)
                !CDIR EXPAND=NG2
                do ig = 1, ng2
                taufor =  forfac(lay) * (forref(indf,ig) + &
                        forfrac(lay) * (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs1+ig) = colh2o(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + taufor
                fracs(lay,ngs1+ig) = fracrefb(ig)
                enddo

                enddo

        end subroutine taugb2

        !----------------------------------------------------------------------------
        subroutine taugb3
                !----------------------------------------------------------------------------
                !
                !     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
                !                           (high key - h2o,co2; high minor - n2o)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Putting the USE statements back in.
                use rrlw_kg03,   only : fracrefa, fracrefb, absa, absb, &
                      ka_mn2o, kb_mn2o, selfref, forref
                !REAL(wp) :: fracrefa(ng3,9) ,fracrefb(ng3,5)
                !REAL(wp) :: absa(585,ng3),  absb(1175,ng3)
                !REAL(wp) :: ka_mn2o(9,19,ng3), kb_mn2o(5,19,ng3)
                !REAL(wp) :: selfref(10,ng3)
                !REAL(wp) :: forref(4,ng3)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                ! ng3 and ngs2 from mo_lrtm_par module. Others from mo_lrtm_kgs.
                INTEGER, PARAMETER :: ng3  = 16
                INTEGER, PARAMETER :: ngs2  = 22
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                integer :: js, js1, jmn2o, jpl
                real(wp) :: speccomb, specparm, specmult, fs
                real(wp) :: speccomb1, specparm1, specmult1, fs1
                real(wp) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, &
                        fmn2o, chi_n2o, ratn2o, adjfac, adjcoln2o
                real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
                real(wp) :: p, p4, fk0, fk1, fk2
                real(wp) :: fac000, fac100, fac200
                real(wp) :: fac010, fac110, fac210
                real(wp) :: fac001, fac101, fac201
                real(wp) :: fac011, fac111, fac211
                real(wp) :: tauself, taufor, n2om1, n2om2, absn2o
                real(wp) :: refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b
                real(wp) :: tau_major(ng3), tau_major1(ng3)


                ! Minor gas mapping levels:
                !     lower - n2o, p = 706.272 mbar, t = 278.94 k
                !     upper - n2o, p = 95.58 mbar, t = 215.7 k

                !  P = 212.725 mb
                refrat_planck_a = chi_mls(1,9)/chi_mls(2,9)

                !  P = 95.58 mb
                refrat_planck_b = chi_mls(1,13)/chi_mls(2,13)

                !  P = 706.270mb
                refrat_m_a = chi_mls(1,3)/chi_mls(2,3)

                !  P = 95.58 mb
                refrat_m_b = chi_mls(1,13)/chi_mls(2,13)

                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature, and appropriate species.  Below laytrop, the water vapor
                ! self-continuum and foreign continuum is interpolated (in temperature)
                ! separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mn2o = colh2o(lay) + refrat_m_a*colco2(lay)
                specparm_mn2o = MIN(colh2o(lay)/speccomb_mn2o,oneminus)
                specmult_mn2o = 8._wp*specparm_mn2o
                jmn2o = 1 + int(specmult_mn2o)
                fmn2o = MOD1(specmult_mn2o)
                !  In atmospheres where the amount of N2O is too great to be considered
                !  a minor species, adjust the column amount of N2O by an empirical factor
                !  to obtain the proper contribution.
                chi_n2o = coln2o(lay)/coldry(lay)
                ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(lay)+1)
                if (ratn2o .gt. 1.5_wp) then
                        adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
                        adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcoln2o = coln2o(lay)
                endif

                speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(3) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(3) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG3
                        tau_major(1:ng3) = speccomb *    &
                                (fac000 * absa(ind0,1:ng3)    + &
                                fac100 * absa(ind0+1,1:ng3)  + &
                                fac200 * absa(ind0+2,1:ng3)  + &
                                fac010 * absa(ind0+9,1:ng3)  + &
                                fac110 * absa(ind0+10,1:ng3) + &
                                fac210 * absa(ind0+11,1:ng3))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG3
                        tau_major(1:ng3) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng3) + &
                                fac100 * absa(ind0,1:ng3)   + &
                                fac000 * absa(ind0+1,1:ng3) + &
                                fac210 * absa(ind0+8,1:ng3) + &
                                fac110 * absa(ind0+9,1:ng3) + &
                                fac010 * absa(ind0+10,1:ng3))
                else
                        !CDIR EXPAND=NG3
                        tau_major(1:ng3) = speccomb *   &
                                (fac000 * absa(ind0,1:ng3)   + &
                                fac100 * absa(ind0+1,1:ng3) + &
                                fac010 * absa(ind0+9,1:ng3) + &
                                fac110 * absa(ind0+10,1:ng3))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG3
                        tau_major1(1:ng3) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng3)    + &
                                fac101 * absa(ind1+1,1:ng3)  + &
                                fac201 * absa(ind1+2,1:ng3)  + &
                                fac011 * absa(ind1+9,1:ng3)  + &
                                fac111 * absa(ind1+10,1:ng3) + &
                                fac211 * absa(ind1+11,1:ng3))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG3
                        tau_major1(1:ng3) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng3) + &
                                fac101 * absa(ind1,1:ng3)   + &
                                fac001 * absa(ind1+1,1:ng3) + &
                                fac211 * absa(ind1+8,1:ng3) + &
                                fac111 * absa(ind1+9,1:ng3) + &
                                fac011 * absa(ind1+10,1:ng3))
                else
                        !CDIR EXPAND=NG3
                        tau_major1(1:ng3) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng3)   + &
                                fac101 * absa(ind1+1,1:ng3) + &
                                fac011 * absa(ind1+9,1:ng3) + &
                                fac111 * absa(ind1+10,1:ng3))
                endif

                !CDIR EXPAND=NG3
                do ig = 1, ng3
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                        (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
                n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                        (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
                absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)

                taug(lay,ngs2+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + adjcoln2o*absn2o
                fracs(lay,ngs2+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 4._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 4._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                fac000 = (1._wp - fs) * fac00(lay)
                fac010 = (1._wp - fs) * fac10(lay)
                fac100 = fs * fac00(lay)
                fac110 = fs * fac10(lay)
                fac001 = (1._wp - fs1) * fac01(lay)
                fac011 = (1._wp - fs1) * fac11(lay)
                fac101 = fs1 * fac01(lay)
                fac111 = fs1 * fac11(lay)

                speccomb_mn2o = colh2o(lay) + refrat_m_b*colco2(lay)
                specparm_mn2o = MIN(colh2o(lay)/speccomb_mn2o,oneminus)
                specmult_mn2o = 4._wp*specparm_mn2o
                jmn2o = 1 + int(specmult_mn2o)
                fmn2o = MOD1(specmult_mn2o)
                !  In atmospheres where the amount of N2O is too great to be considered
                !  a minor species, adjust the column amount of N2O by an empirical factor
                !  to obtain the proper contribution.
                chi_n2o = coln2o(lay)/coldry(lay)
                ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(lay)+1)
                if (ratn2o .gt. 1.5_wp) then
                        adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
                        adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcoln2o = coln2o(lay)
                endif

                speccomb_planck = colh2o(lay)+refrat_planck_b*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 4._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(3) + js
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(3) + js1
                indf = indfor(lay)
                indm = indminor(lay)
                !CDIR EXPAND=NG3
                do ig = 1, ng3
                taufor = forfac(lay) * (forref(indf,ig) + &
                        forfrac(lay) * (forref(indf+1,ig) - forref(indf,ig)))
                n2om1 = kb_mn2o(jmn2o,indm,ig) + fmn2o * &
                        (kb_mn2o(jmn2o+1,indm,ig)-kb_mn2o(jmn2o,indm,ig))
                n2om2 = kb_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                        (kb_mn2o(jmn2o+1,indm+1,ig)-kb_mn2o(jmn2o,indm+1,ig))
                absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)
                taug(lay,ngs2+ig) = speccomb * &
                        (fac000 * absb(ind0,ig) + &
                        fac100 * absb(ind0+1,ig) + &
                        fac010 * absb(ind0+5,ig) + &
                        fac110 * absb(ind0+6,ig)) &
                        + speccomb1 * &
                        (fac001 * absb(ind1,ig) +  &
                        fac101 * absb(ind1+1,ig) + &
                        fac011 * absb(ind1+5,ig) + &
                        fac111 * absb(ind1+6,ig))  &
                        + taufor &
                        + adjcoln2o*absn2o
                fracs(lay,ngs2+ig) = fracrefb(ig,jpl) + fpl * &
                        (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max

                ! >> sylvia_20200306
                ! Removing ixc0 = ixc(lay); do ixp = 1, ixc0; jl = ixlow(ixp,lay)

                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mn2o = colh2o(lay) + refrat_m_a*colco2(lay)
                specparm_mn2o = MIN(colh2o(lay)/speccomb_mn2o,oneminus)
                specmult_mn2o = 8._wp*specparm_mn2o
                jmn2o = 1 + int(specmult_mn2o)
                fmn2o = MOD1(specmult_mn2o)
                !  In atmospheres where the amount of N2O is too great to be considered
                !  a minor species, adjust the column amount of N2O by an empirical factor
                !  to obtain the proper contribution.
                chi_n2o = coln2o(lay)/coldry(lay)
                ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(lay)+1)
                if (ratn2o .gt. 1.5_wp) then
                        adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
                        adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcoln2o = coln2o(lay)
                endif

                speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(3) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(3) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG3
                        tau_major(1:ng3) = speccomb *    &
                                (fac000 * absa(ind0,1:ng3)    + &
                                fac100 * absa(ind0+1,1:ng3)  + &
                                fac200 * absa(ind0+2,1:ng3)  + &
                                fac010 * absa(ind0+9,1:ng3)  + &
                                fac110 * absa(ind0+10,1:ng3) + &
                                fac210 * absa(ind0+11,1:ng3))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG3
                        tau_major(1:ng3) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng3) + &
                                fac100 * absa(ind0,1:ng3)   + &
                                fac000 * absa(ind0+1,1:ng3) + &
                                fac210 * absa(ind0+8,1:ng3) + &
                                fac110 * absa(ind0+9,1:ng3) + &
                                fac010 * absa(ind0+10,1:ng3))
                else
                        !CDIR EXPAND=NG3
                        tau_major(1:ng3) = speccomb *   &
                                (fac000 * absa(ind0,1:ng3)   + &
                                fac100 * absa(ind0+1,1:ng3) + &
                                fac010 * absa(ind0+9,1:ng3) + &
                                fac110 * absa(ind0+10,1:ng3))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG3
                        tau_major1(1:ng3) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng3)    + &
                                fac101 * absa(ind1+1,1:ng3)  + &
                                fac201 * absa(ind1+2,1:ng3)  + &
                                fac011 * absa(ind1+9,1:ng3)  + &
                                fac111 * absa(ind1+10,1:ng3) + &
                                fac211 * absa(ind1+11,1:ng3))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG3
                        tau_major1(1:ng3) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng3) + &
                                fac101 * absa(ind1,1:ng3)   + &
                                fac001 * absa(ind1+1,1:ng3) + &
                                fac211 * absa(ind1+8,1:ng3) + &
                                fac111 * absa(ind1+9,1:ng3) + &
                                fac011 * absa(ind1+10,1:ng3))
                else
                        !CDIR EXPAND=NG3
                        tau_major1(1:ng3) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng3)   + &
                                fac101 * absa(ind1+1,1:ng3) + &
                                fac011 * absa(ind1+9,1:ng3) + &
                                fac111 * absa(ind1+10,1:ng3))
                endif

                !CDIR EXPAND=NG3
                do ig = 1, ng3
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                        (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
                n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                        (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
                absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)

                taug(lay,ngs2+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + adjcoln2o*absn2o
                fracs(lay,ngs2+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                ! until here

                ! Upper atmosphere part
                ! >> sylvia_20200306
                ! Removing ixc0 = kproma - ixc0; do ixp = 1, ixc0; jl = ixhigh(ixp,lay)
                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 4._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 4._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                fac000 = (1._wp - fs) * fac00(lay)
                fac010 = (1._wp - fs) * fac10(lay)
                fac100 = fs * fac00(lay)
                fac110 = fs * fac10(lay)
                fac001 = (1._wp - fs1) * fac01(lay)
                fac011 = (1._wp - fs1) * fac11(lay)
                fac101 = fs1 * fac01(lay)
                fac111 = fs1 * fac11(lay)

                speccomb_mn2o = colh2o(lay) + refrat_m_b*colco2(lay)
                specparm_mn2o = MIN(colh2o(lay)/speccomb_mn2o,oneminus)
                specmult_mn2o = 4._wp*specparm_mn2o
                jmn2o = 1 + int(specmult_mn2o)
                fmn2o = MOD1(specmult_mn2o)
                !  In atmospheres where the amount of N2O is too great to be considered
                !  a minor species, adjust the column amount of N2O by an empirical factor
                !  to obtain the proper contribution.
                chi_n2o = coln2o(lay)/coldry(lay)
                ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(lay)+1)
                if (ratn2o .gt. 1.5_wp) then
                        adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
                        adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcoln2o = coln2o(lay)
                endif

                speccomb_planck = colh2o(lay)+refrat_planck_b*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 4._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(3) + js
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(3) + js1
                indf = indfor(lay)
                indm = indminor(lay)
                !CDIR EXPAND=NG3
                do ig = 1, ng3
                taufor = forfac(lay) * (forref(indf,ig) + &
                        forfrac(lay) * (forref(indf+1,ig) - forref(indf,ig)))
                n2om1 = kb_mn2o(jmn2o,indm,ig) + fmn2o * &
                        (kb_mn2o(jmn2o+1,indm,ig)-kb_mn2o(jmn2o,indm,ig))
                n2om2 = kb_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                        (kb_mn2o(jmn2o+1,indm+1,ig)-kb_mn2o(jmn2o,indm+1,ig))
                absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)
                taug(lay,ngs2+ig) = speccomb * &
                        (fac000 * absb(ind0,ig) + &
                        fac100 * absb(ind0+1,ig) + &
                        fac010 * absb(ind0+5,ig) + &
                        fac110 * absb(ind0+6,ig)) &
                        + speccomb1 * &
                        (fac001 * absb(ind1,ig) +  &
                        fac101 * absb(ind1+1,ig) + &
                        fac011 * absb(ind1+5,ig) + &
                        fac111 * absb(ind1+6,ig))  &
                        + taufor &
                        + adjcoln2o*absn2o
                fracs(lay,ngs2+ig) = fracrefb(ig,jpl) + fpl * &
                        (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
                enddo
                enddo

        end subroutine taugb3

        !----------------------------------------------------------------------------
        subroutine taugb4
                !----------------------------------------------------------------------------
                !
                !     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Putting the USE statements back in.
                use rrlw_kg04,   only : fracrefa, fracrefb, absa, absb, &
                        selfref, forref
                !REAL(wp) :: fracrefa(ng4,9), fracrefb(ng4,5)
                !REAL(wp) :: absa(585,ng4), absb(1175,ng4)
                !REAL(wp) :: selfref(10,ng4), forref(4,ng4)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                ! Adding in variables from mo_lrtm_par and mo_lrtm_kgs
                ! >> sylvia_20200318
                ! I missed kb_mn2o here.
                INTEGER, PARAMETER :: ng4  = 14
                INTEGER, PARAMETER :: ngs3  = 16
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf
                integer :: js, js1, jpl
                real(wp) :: speccomb, specparm, specmult, fs
                real(wp) :: speccomb1, specparm1, specmult1, fs1
                real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
                real(wp) :: p, p4, fk0, fk1, fk2
                real(wp) :: fac000, fac100, fac200
                real(wp) :: fac010, fac110, fac210
                real(wp) :: fac001, fac101, fac201
                real(wp) :: fac011, fac111, fac211
                real(wp) :: tauself, taufor
                real(wp) :: refrat_planck_a, refrat_planck_b
                real(wp) :: tau_major(ng4), tau_major1(ng4)


                ! P =   142.5940 mb
                refrat_planck_a = chi_mls(1,11)/chi_mls(2,11)

                ! P = 95.58350 mb
                refrat_planck_b = chi_mls(3,13)/chi_mls(2,13)

                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature, and appropriate species.  Below laytrop, the water
                ! vapor self-continuum and foreign continuum is interpolated (in temperature)
                ! separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(4) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(4) + js1
                inds = indself(lay)
                indf = indfor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG4
                        tau_major(1:ng4) = speccomb *    &
                                (fac000 * absa(ind0,1:ng4)    + &
                                fac100 * absa(ind0+1,1:ng4)  + &
                                fac200 * absa(ind0+2,1:ng4)  + &
                                fac010 * absa(ind0+9,1:ng4)  + &
                                fac110 * absa(ind0+10,1:ng4) + &
                                fac210 * absa(ind0+11,1:ng4))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG4
                        tau_major(1:ng4) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng4) + &
                                fac100 * absa(ind0,1:ng4)   + &
                                fac000 * absa(ind0+1,1:ng4) + &
                                fac210 * absa(ind0+8,1:ng4) + &
                                fac110 * absa(ind0+9,1:ng4) + &
                                fac010 * absa(ind0+10,1:ng4))
                else
                        !CDIR EXPAND=NG4
                        tau_major(1:ng4) = speccomb *   &
                                (fac000 * absa(ind0,1:ng4)   + &
                                fac100 * absa(ind0+1,1:ng4) + &
                                fac010 * absa(ind0+9,1:ng4) + &
                                fac110 * absa(ind0+10,1:ng4))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG4
                        tau_major1(1:ng4) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng4)    + &
                                fac101 * absa(ind1+1,1:ng4)  + &
                                fac201 * absa(ind1+2,1:ng4)  + &
                                fac011 * absa(ind1+9,1:ng4)  + &
                                fac111 * absa(ind1+10,1:ng4) + &
                                fac211 * absa(ind1+11,1:ng4))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG4
                        tau_major1(1:ng4) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng4) + &
                                fac101 * absa(ind1,1:ng4)   + &
                                fac001 * absa(ind1+1,1:ng4) + &
                                fac211 * absa(ind1+8,1:ng4) + &
                                fac111 * absa(ind1+9,1:ng4) + &
                                fac011 * absa(ind1+10,1:ng4))
                else
                        !CDIR EXPAND=NG4
                        tau_major1(1:ng4) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng4)   + &
                                fac101 * absa(ind1+1,1:ng4) + &
                                fac011 * absa(ind1+9,1:ng4) + &
                                fac111 * absa(ind1+10,1:ng4))
                endif

                !CDIR EXPAND=NG4
                do ig = 1, ng4
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))

                taug(lay,ngs3+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor
                fracs(lay,ngs3+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers

                speccomb = colo3(lay) + rat_o3co2(lay)*colco2(lay)
                specparm = MIN(colo3(lay)/speccomb,oneminus)
                specmult = 4._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colo3(lay) + rat_o3co2_1(lay)*colco2(lay)
                specparm1 = MIN(colo3(lay)/speccomb1,oneminus)
                specmult1 = 4._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                fac000 = (1._wp - fs) * fac00(lay)
                fac010 = (1._wp - fs) * fac10(lay)
                fac100 = fs * fac00(lay)
                fac110 = fs * fac10(lay)
                fac001 = (1._wp - fs1) * fac01(lay)
                fac011 = (1._wp - fs1) * fac11(lay)
                fac101 = fs1 * fac01(lay)
                fac111 = fs1 * fac11(lay)

                speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay)
                specparm_planck = MIN(colo3(lay)/speccomb_planck,oneminus)
                specmult_planck = 4._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(4) + js
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(4) + js1
                !CDIR EXPAND=NG4
                do ig = 1, ng4
                taug(lay,ngs3+ig) =  speccomb * &
                        (fac000 * absb(ind0,ig) + &
                        fac100 * absb(ind0+1,ig) + &
                        fac010 * absb(ind0+5,ig) + &
                        fac110 * absb(ind0+6,ig)) &
                        + speccomb1 * &
                        (fac001 * absb(ind1,ig) +  &
                        fac101 * absb(ind1+1,ig) + &
                        fac011 * absb(ind1+5,ig) + &
                        fac111 * absb(ind1+6,ig))
                fracs(lay,ngs3+ig) = fracrefb(ig,jpl) + fpl * &
                        (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
                enddo
                enddo

                ! Empirical modification to code to improve stratospheric cooling rates
                ! for co2.  Revised to apply weighting for g-point reduction in this band.
                ! >> sylvia_20200306
                ! Removing inner loop up to nproma.
                do lay = laytrop_max+1, nlayers
                taug(lay,ngs3+8)=taug(lay,ngs3+8)*0.92_wp
                taug(lay,ngs3+9)=taug(lay,ngs3+9)*0.88_wp
                taug(lay,ngs3+10)=taug(lay,ngs3+10)*1.07_wp
                taug(lay,ngs3+11)=taug(lay,ngs3+11)*1.1_wp
                taug(lay,ngs3+12)=taug(lay,ngs3+12)*0.99_wp
                taug(lay,ngs3+13)=taug(lay,ngs3+13)*0.88_wp
                taug(lay,ngs3+14)=taug(lay,ngs3+14)*0.943_wp
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                DO lay = laytrop_min+1, laytrop_max

                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(4) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(4) + js1
                inds = indself(lay)
                indf = indfor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG4
                        tau_major(1:ng4) = speccomb *    &
                                (fac000 * absa(ind0,1:ng4)    + &
                                fac100 * absa(ind0+1,1:ng4)  + &
                                fac200 * absa(ind0+2,1:ng4)  + &
                                fac010 * absa(ind0+9,1:ng4)  + &
                                fac110 * absa(ind0+10,1:ng4) + &
                                fac210 * absa(ind0+11,1:ng4))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG4
                        tau_major(1:ng4) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng4) + &
                                fac100 * absa(ind0,1:ng4)   + &
                                fac000 * absa(ind0+1,1:ng4) + &
                                fac210 * absa(ind0+8,1:ng4) + &
                                fac110 * absa(ind0+9,1:ng4) + &
                                fac010 * absa(ind0+10,1:ng4))
                else
                        !CDIR EXPAND=NG4
                        tau_major(1:ng4) = speccomb *   &
                                (fac000 * absa(ind0,1:ng4)   + &
                                fac100 * absa(ind0+1,1:ng4) + &
                                fac010 * absa(ind0+9,1:ng4) + &
                                fac110 * absa(ind0+10,1:ng4))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG4
                        tau_major1(1:ng4) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng4)    + &
                                fac101 * absa(ind1+1,1:ng4)  + &
                                fac201 * absa(ind1+2,1:ng4)  + &
                                fac011 * absa(ind1+9,1:ng4)  + &
                                fac111 * absa(ind1+10,1:ng4) + &
                                fac211 * absa(ind1+11,1:ng4))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG4
                        tau_major1(1:ng4) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng4) + &
                                fac101 * absa(ind1,1:ng4)   + &
                                fac001 * absa(ind1+1,1:ng4) + &
                                fac211 * absa(ind1+8,1:ng4) + &
                                fac111 * absa(ind1+9,1:ng4) + &
                                fac011 * absa(ind1+10,1:ng4))
                else
                        !CDIR EXPAND=NG4
                        tau_major1(1:ng4) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng4)   + &
                                fac101 * absa(ind1+1,1:ng4) + &
                                fac011 * absa(ind1+9,1:ng4) + &
                                fac111 * absa(ind1+10,1:ng4))
                endif

                !CDIR EXPAND=NG4
                do ig = 1, ng4
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))

                taug(lay,ngs3+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor
                fracs(lay,ngs3+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo

                ! Upper atmosphere part
                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colo3(lay) + rat_o3co2(lay)*colco2(lay)
                specparm = MIN(colo3(lay)/speccomb,oneminus)
                specmult = 4._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colo3(lay) + rat_o3co2_1(lay)*colco2(lay)
                specparm1 = MIN(colo3(lay)/speccomb1,oneminus)
                specmult1 = 4._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                fac000 = (1._wp - fs) * fac00(lay)
                fac010 = (1._wp - fs) * fac10(lay)
                fac100 = fs * fac00(lay)
                fac110 = fs * fac10(lay)
                fac001 = (1._wp - fs1) * fac01(lay)
                fac011 = (1._wp - fs1) * fac11(lay)
                fac101 = fs1 * fac01(lay)
                fac111 = fs1 * fac11(lay)

                speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay)
                specparm_planck = MIN(colo3(lay)/speccomb_planck,oneminus)
                specmult_planck = 4._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(4) + js
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(4) + js1
                !CDIR EXPAND=NG4
                do ig = 1, ng4
                taug(lay,ngs3+ig) =  speccomb * &
                        (fac000 * absb(ind0,ig) + &
                        fac100 * absb(ind0+1,ig) + &
                        fac010 * absb(ind0+5,ig) + &
                        fac110 * absb(ind0+6,ig)) &
                        + speccomb1 * &
                        (fac001 * absb(ind1,ig) +  &
                        fac101 * absb(ind1+1,ig) + &
                        fac011 * absb(ind1+5,ig) + &
                        fac111 * absb(ind1+6,ig))
                fracs(lay,ngs3+ig) = fracrefb(ig,jpl) + fpl * &
                        (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
                enddo
                enddo

                ! Empirical modification to code to improve stratospheric cooling rates
                ! for co2.  Revised to apply weighting for g-point reduction in this band.
                !CDIR NODEP,VOVERTAKE,VOB
                ! >> sylvia_20200306
                ! Removing the section below as ixp is only used to set jl which is no longer
                ! a relevant index.
                !do ixp = 1, ixc0
                !  jl = ixhigh(ixp,lay)
                !  taug(lay,ngs3+8)=taug(lay,ngs3+8)*0.92_wp
                !  taug(lay,ngs3+9)=taug(lay,ngs3+9)*0.88_wp
                !  taug(lay,ngs3+10)=taug(lay,ngs3+10)*1.07_wp
                !  taug(lay,ngs3+11)=taug(lay,ngs3+11)*1.1_wp
                !  taug(lay,ngs3+12)=taug(lay,ngs3+12)*0.99_wp
                !  taug(lay,ngs3+13)=taug(lay,ngs3+13)*0.88_wp
                !  taug(lay,ngs3+14)=taug(lay,ngs3+14)*0.943_wp
                !enddo

        end subroutine taugb4



        !----------------------------------------------------------------------------
        subroutine taugb5
                !----------------------------------------------------------------------------
                !
                !     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
                !                           (high key - o3,co2)
                !----------------------------------------------------------------------------

                ! ------- Modules -------
                ! >> sylvia_20200320
                ! Adding back in USE statements.
                use rrlw_kg05,   only : fracrefa, fracrefb, absa, absb, &
                         ka_mo3, selfref, forref, ccl4
                !REAL(wp) :: fracrefa(ng5,9) ,fracrefb(ng5,5)
                !REAL(wp) :: absa(585,ng5)
                !REAL(wp) :: absb(1175,ng5)
                !REAL(wp) :: ka_mo3(9,19,ng5)
                !REAL(wp) :: selfref(10,ng5)
                !REAL(wp) :: forref(4,ng5)
                !REAL(wp) :: ccl4(ng5)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng5  = 16
                INTEGER, PARAMETER :: ngs4 = 52
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                integer :: js, js1, jmo3, jpl
                real(wp) :: speccomb, specparm, specmult, fs
                real(wp) :: speccomb1, specparm1, specmult1, fs1
                real(wp) :: speccomb_mo3, specparm_mo3, specmult_mo3, fmo3
                real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
                real(wp) :: p, p4, fk0, fk1, fk2
                real(wp) :: fac000, fac100, fac200
                real(wp) :: fac010, fac110, fac210
                real(wp) :: fac001, fac101, fac201
                real(wp) :: fac011, fac111, fac211
                real(wp) :: tauself, taufor, o3m1, o3m2, abso3
                real(wp) :: refrat_planck_a, refrat_planck_b, refrat_m_a
                real(wp) :: tau_major(ng5), tau_major1(ng5)


                ! Minor gas mapping level :
                !     lower - o3, p = 317.34 mbar, t = 240.77 k
                !     lower - ccl4

                ! Calculate reference ratio to be used in calculation of Planck
                ! fraction in lower/upper atmosphere.

                ! P = 473.420 mb
                refrat_planck_a = chi_mls(1,5)/chi_mls(2,5)

                ! P = 0.2369 mb
                refrat_planck_b = chi_mls(3,43)/chi_mls(2,43)

                ! P = 317.3480
                refrat_m_a = chi_mls(1,7)/chi_mls(2,7)

                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature, and appropriate species.  Below laytrop, the
                ! water vapor self-continuum and foreign continuum is
                ! interpolated (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mo3 = colh2o(lay) + refrat_m_a*colco2(lay)
                specparm_mo3 = MIN(colh2o(lay)/speccomb_mo3,oneminus)
                specmult_mo3 = 8._wp*specparm_mo3
                jmo3 = 1 + int(specmult_mo3)
                fmo3 = MOD1(specmult_mo3)

                speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(5) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(5) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG5
                        tau_major(1:ng5) = speccomb *    &
                                (fac000 * absa(ind0,1:ng5)    + &
                                fac100 * absa(ind0+1,1:ng5)  + &
                                fac200 * absa(ind0+2,1:ng5)  + &
                                fac010 * absa(ind0+9,1:ng5)  + &
                                fac110 * absa(ind0+10,1:ng5) + &
                                fac210 * absa(ind0+11,1:ng5))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG5
                        tau_major(1:ng5) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng5) + &
                                fac100 * absa(ind0,1:ng5)   + &
                                fac000 * absa(ind0+1,1:ng5) + &
                                fac210 * absa(ind0+8,1:ng5) + &
                                fac110 * absa(ind0+9,1:ng5) + &
                                fac010 * absa(ind0+10,1:ng5))
                else
                        !CDIR EXPAND=NG5
                        tau_major(1:ng5) = speccomb *   &
                                (fac000 * absa(ind0,1:ng5)   + &
                                fac100 * absa(ind0+1,1:ng5) + &
                                fac010 * absa(ind0+9,1:ng5) + &
                                fac110 * absa(ind0+10,1:ng5))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG5
                        tau_major1(1:ng5) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng5)    + &
                                fac101 * absa(ind1+1,1:ng5)  + &
                                fac201 * absa(ind1+2,1:ng5)  + &
                                fac011 * absa(ind1+9,1:ng5)  + &
                                fac111 * absa(ind1+10,1:ng5) + &
                                fac211 * absa(ind1+11,1:ng5))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG5
                        tau_major1(1:ng5) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng5) + &
                                fac101 * absa(ind1,1:ng5)   + &
                                fac001 * absa(ind1+1,1:ng5) + &
                                fac211 * absa(ind1+8,1:ng5) + &
                                fac111 * absa(ind1+9,1:ng5) + &
                                fac011 * absa(ind1+10,1:ng5))
                else
                        !CDIR EXPAND=NG5
                        tau_major1(1:ng5) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng5)   + &
                                fac101 * absa(ind1+1,1:ng5) + &
                                fac011 * absa(ind1+9,1:ng5) + &
                                fac111 * absa(ind1+10,1:ng5))
                endif

                !CDIR EXPAND=NG5
                do ig = 1, ng5
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                o3m1 = ka_mo3(jmo3,indm,ig) + fmo3 * &
                        (ka_mo3(jmo3+1,indm,ig)-ka_mo3(jmo3,indm,ig))
                o3m2 = ka_mo3(jmo3,indm+1,ig) + fmo3 * &
                        (ka_mo3(jmo3+1,indm+1,ig)-ka_mo3(jmo3,indm+1,ig))
                abso3 = o3m1 + minorfrac(lay)*(o3m2-o3m1)

                taug(lay,ngs4+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + abso3*colo3(lay) &
                        + wx(1,lay) * ccl4(ig)
                fracs(lay,ngs4+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers

                speccomb = colo3(lay) + rat_o3co2(lay)*colco2(lay)
                specparm = MIN(colo3(lay)/speccomb,oneminus)
                specmult = 4._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colo3(lay) + rat_o3co2_1(lay)*colco2(lay)
                specparm1 = MIN(colo3(lay)/speccomb1,oneminus)
                specmult1 = 4._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                fac000 = (1._wp - fs) * fac00(lay)
                fac010 = (1._wp - fs) * fac10(lay)
                fac100 = fs * fac00(lay)
                fac110 = fs * fac10(lay)
                fac001 = (1._wp - fs1) * fac01(lay)
                fac011 = (1._wp - fs1) * fac11(lay)
                fac101 = fs1 * fac01(lay)
                fac111 = fs1 * fac11(lay)

                speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay)
                specparm_planck = MIN(colo3(lay)/speccomb_planck,oneminus)
                specmult_planck = 4._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(5) + js
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(5) + js1
                !CDIR EXPAND=NG5
                do ig = 1, ng5
                taug(lay,ngs4+ig) = speccomb * &
                        (fac000 * absb(ind0,ig) + &
                        fac100 * absb(ind0+1,ig) + &
                        fac010 * absb(ind0+5,ig) + &
                        fac110 * absb(ind0+6,ig)) &
                        + speccomb1 * &
                        (fac001 * absb(ind1,ig) + &
                        fac101 * absb(ind1+1,ig) + &
                        fac011 * absb(ind1+5,ig) + &
                        fac111 * absb(ind1+6,ig))  &
                        + wx(1,lay) * ccl4(ig)
                fracs(lay,ngs4+ig) = fracrefb(ig,jpl) + fpl * &
                        (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max

                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mo3 = colh2o(lay) + refrat_m_a*colco2(lay)
                specparm_mo3 = MIN(colh2o(lay)/speccomb_mo3,oneminus)
                specmult_mo3 = 8._wp*specparm_mo3
                jmo3 = 1 + int(specmult_mo3)
                fmo3 = MOD1(specmult_mo3)

                speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(5) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(5) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG5
                        tau_major(1:ng5) = speccomb *    &
                                (fac000 * absa(ind0,1:ng5)    + &
                                fac100 * absa(ind0+1,1:ng5)  + &
                                fac200 * absa(ind0+2,1:ng5)  + &
                                fac010 * absa(ind0+9,1:ng5)  + &
                                fac110 * absa(ind0+10,1:ng5) + &
                                fac210 * absa(ind0+11,1:ng5))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG5
                        tau_major(1:ng5) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng5) + &
                                fac100 * absa(ind0,1:ng5)   + &
                                fac000 * absa(ind0+1,1:ng5) + &
                                fac210 * absa(ind0+8,1:ng5) + &
                                fac110 * absa(ind0+9,1:ng5) + &
                                fac010 * absa(ind0+10,1:ng5))
                else
                        !CDIR EXPAND=NG5
                        tau_major(1:ng5) = speccomb *   &
                                (fac000 * absa(ind0,1:ng5)   + &
                                fac100 * absa(ind0+1,1:ng5) + &
                                fac010 * absa(ind0+9,1:ng5) + &
                                fac110 * absa(ind0+10,1:ng5))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG5
                        tau_major1(1:ng5) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng5)    + &
                                fac101 * absa(ind1+1,1:ng5)  + &
                                fac201 * absa(ind1+2,1:ng5)  + &
                                fac011 * absa(ind1+9,1:ng5)  + &
                                fac111 * absa(ind1+10,1:ng5) + &
                                fac211 * absa(ind1+11,1:ng5))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG5
                        tau_major1(1:ng5) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng5) + &
                                fac101 * absa(ind1,1:ng5)   + &
                                fac001 * absa(ind1+1,1:ng5) + &
                                fac211 * absa(ind1+8,1:ng5) + &
                                fac111 * absa(ind1+9,1:ng5) + &
                                fac011 * absa(ind1+10,1:ng5))
                else
                        !CDIR EXPAND=NG5
                        tau_major1(1:ng5) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng5)   + &
                                fac101 * absa(ind1+1,1:ng5) + &
                                fac011 * absa(ind1+9,1:ng5) + &
                                fac111 * absa(ind1+10,1:ng5))
                endif

                !CDIR EXPAND=NG5
                do ig = 1, ng5
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                o3m1 = ka_mo3(jmo3,indm,ig) + fmo3 * &
                        (ka_mo3(jmo3+1,indm,ig)-ka_mo3(jmo3,indm,ig))
                o3m2 = ka_mo3(jmo3,indm+1,ig) + fmo3 * &
                        (ka_mo3(jmo3+1,indm+1,ig)-ka_mo3(jmo3,indm+1,ig))
                abso3 = o3m1 + minorfrac(lay)*(o3m2-o3m1)

                taug(lay,ngs4+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + abso3*colo3(lay) &
                        + wx(1,lay) * ccl4(ig)
                fracs(lay,ngs4+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo

                ! Upper atmosphere part
                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colo3(lay) + rat_o3co2(lay)*colco2(lay)
                specparm = MIN(colo3(lay)/speccomb,oneminus)
                specmult = 4._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colo3(lay) + rat_o3co2_1(lay)*colco2(lay)
                specparm1 = MIN(colo3(lay)/speccomb1,oneminus)
                specmult1 = 4._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                fac000 = (1._wp - fs) * fac00(lay)
                fac010 = (1._wp - fs) * fac10(lay)
                fac100 = fs * fac00(lay)
                fac110 = fs * fac10(lay)
                fac001 = (1._wp - fs1) * fac01(lay)
                fac011 = (1._wp - fs1) * fac11(lay)
                fac101 = fs1 * fac01(lay)
                fac111 = fs1 * fac11(lay)

                speccomb_planck = colo3(lay)+refrat_planck_b*colco2(lay)
                specparm_planck = MIN(colo3(lay)/speccomb_planck,oneminus)
                specmult_planck = 4._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(5) + js
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(5) + js1
                !CDIR EXPAND=NG5
                do ig = 1, ng5
                taug(lay,ngs4+ig) = speccomb * &
                        (fac000 * absb(ind0,ig) + &
                        fac100 * absb(ind0+1,ig) + &
                        fac010 * absb(ind0+5,ig) + &
                        fac110 * absb(ind0+6,ig)) &
                        + speccomb1 * &
                        (fac001 * absb(ind1,ig) + &
                        fac101 * absb(ind1+1,ig) + &
                        fac011 * absb(ind1+5,ig) + &
                        fac111 * absb(ind1+6,ig))  &
                        + wx(1,lay) * ccl4(ig)
                fracs(lay,ngs4+ig) = fracrefb(ig,jpl) + fpl * &
                        (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
                enddo
                enddo

        end subroutine taugb5

        !----------------------------------------------------------------------------
        subroutine taugb6
                !----------------------------------------------------------------------------
                !
                !     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
                !                           (high key - nothing; high minor - cfc11, cfc12)
                !----------------------------------------------------------------------------

                ! ------- Modules -------
                ! >> sylvia_20200320
                ! Add the USE statements in again.
                use rrlw_kg06,   only : fracrefa, absa, ka_mco2, &
                        selfref, forref, cfc11adj, cfc12
                !REAL(wp) , DIMENSION(ng6) :: fracrefa
                !REAL(wp) :: absa(65,ng6)
                !REAL(wp) :: ka_mco2(19,ng6)
                !REAL(wp) :: selfref(10,ng6)
                !REAL(wp) :: forref(4,ng6)
                !REAL(wp) , DIMENSION(ng6) :: cfc11adj
                !REAL(wp) , DIMENSION(ng6) :: cfc12
                ! << sylvia_20200310

                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng6  = 8
                INTEGER, PARAMETER :: ngs5  = 68
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                real(wp) :: chi_co2, ratco2, adjfac, adjcolco2
                real(wp) :: tauself, taufor, absco2


                ! Minor gas mapping level:
                !     lower - co2, p = 706.2720 mb, t = 294.2 k
                !     upper - cfc11, cfc12

                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature. The water vapor self-continuum and foreign continuum
                ! is interpolated (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                ! In atmospheres where the amount of CO2 is too great to be considered
                ! a minor species, adjust the column amount of CO2 by an empirical factor
                ! to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.77_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(6) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(6) + 1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)
                !CDIR EXPAND=NG6
                do ig = 1, ng6
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                absco2 =  (ka_mco2(indm,ig) + minorfrac(lay) * &
                        (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
                taug(lay,ngs5+ig) = colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) +  &
                        fac11(lay) * absa(ind1+1,ig))  &
                        + tauself + taufor &
                        + adjcolco2 * absco2 &
                        + wx(2,lay) * cfc11adj(ig) &
                        + wx(3,lay) * cfc12(ig)
                fracs(lay,ngs5+ig) = fracrefa(ig)
                enddo
                enddo

                ! Upper atmosphere loop
                ! Nothing important goes on above laytrop in this band.
                do ig = 1, ng6
                do lay = laytrop_max+1, nlayers
                taug(lay,ngs5+ig) = 0.0_wp &
                        + wx(2,lay) * cfc11adj(ig) &
                        + wx(3,lay) * cfc12(ig)
                fracs(lay,ngs5+ig) = fracrefa(ig)
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max
                !CDIR NODEP,VOVERTAKE,VOB

                ! In atmospheres where the amount of CO2 is too great to be considered
                ! a minor species, adjust the column amount of CO2 by an empirical factor
                ! to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.77_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(6) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(6) + 1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)
                !CDIR EXPAND=NG6
                do ig = 1, ng6
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                absco2 =  (ka_mco2(indm,ig) + minorfrac(lay) * &
                        (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
                taug(lay,ngs5+ig) = colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) +  &
                        fac11(lay) * absa(ind1+1,ig))  &
                        + tauself + taufor &
                        + adjcolco2 * absco2 &
                        + wx(2,lay) * cfc11adj(ig) &
                        + wx(3,lay) * cfc12(ig)
                fracs(lay,ngs5+ig) = fracrefa(ig)
                enddo

                ! Upper atmosphere part
                ! Nothing important goes on above laytrop in this band.

                do ig = 1, ng6
                !CDIR NODEP,VOVERTAKE,VOB
                taug(lay,ngs5+ig) = 0.0_wp &
                        + wx(2,lay) * cfc11adj(ig) &
                        + wx(3,lay) * cfc12(ig)
                fracs(lay,ngs5+ig) = fracrefa(ig)
                enddo

                enddo

        end subroutine taugb6

        !----------------------------------------------------------------------------
        subroutine taugb7
                !----------------------------------------------------------------------------
                !
                !     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
                !                            (high key - o3; high minor - co2)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Add in the USE statements again.
                use rrlw_kg07,   only : fracrefa, fracrefb, absa, absb, &
                         ka_mco2, kb_mco2, selfref, forref
                !REAL(wp) , DIMENSION(ng7) :: fracrefb
                !REAL(wp) :: fracrefa(ng7,9)
                !REAL(wp) :: absa(585,ng7)
                !REAL(wp) :: absb(235,ng7)
                !REAL(wp) :: ka_mco2(9,19,ng7)
                !REAL(wp) :: kb_mco2(19,ng7)
                !REAL(wp) :: selfref(10,ng7)
                !REAL(wp) :: forref(4,ng7)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng7  = 12
                INTEGER, PARAMETER :: ngs6  = 76
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                integer :: js, js1, jmco2, jpl
                real(wp) :: speccomb, specparm, specmult, fs
                real(wp) :: speccomb1, specparm1, specmult1, fs1
                real(wp) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
                real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
                real(wp) :: p, p4, fk0, fk1, fk2
                real(wp) :: fac000, fac100, fac200
                real(wp) :: fac010, fac110, fac210
                real(wp) :: fac001, fac101, fac201
                real(wp) :: fac011, fac111, fac211
                real(wp) :: tauself, taufor, co2m1, co2m2, absco2
                real(wp) :: chi_co2, ratco2, adjfac, adjcolco2
                real(wp) :: refrat_planck_a, refrat_m_a
                real(wp) :: tau_major(ng7), tau_major1(ng7)


                ! Minor gas mapping level :
                !     lower - co2, p = 706.2620 mbar, t= 278.94 k
                !     upper - co2, p = 12.9350 mbar, t = 234.01 k

                ! Calculate reference ratio to be used in calculation of Planck
                ! fraction in lower atmosphere.

                ! P = 706.2620 mb
                refrat_planck_a = chi_mls(1,3)/chi_mls(3,3)

                ! P = 706.2720 mb
                refrat_m_a = chi_mls(1,3)/chi_mls(3,3)

                ! Compute the optical depth by interpolating in ln(pressure),
                ! temperature, and appropriate species.  Below laytrop, the water
                ! vapor self-continuum and foreign continuum is interpolated
                ! (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                speccomb = colh2o(lay) + rat_h2oo3(lay)*colo3(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oo3_1(lay)*colo3(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mco2 = colh2o(lay) + refrat_m_a*colo3(lay)
                specparm_mco2 = MIN(colh2o(lay)/speccomb_mco2,oneminus)
                specmult_mco2 = 8._wp*specparm_mco2

                jmco2 = 1 + int(specmult_mco2)
                fmco2 = MOD1(specmult_mco2)

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 3.0_wp+(ratco2-3.0_wp)**0.79_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                speccomb_planck = colh2o(lay)+refrat_planck_a*colo3(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(7) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(7) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG7
                        tau_major(1:ng7) = speccomb *    &
                                (fac000 * absa(ind0,1:ng7)    + &
                                fac100 * absa(ind0+1,1:ng7)  + &
                                fac200 * absa(ind0+2,1:ng7)  + &
                                fac010 * absa(ind0+9,1:ng7)  + &
                                fac110 * absa(ind0+10,1:ng7) + &
                                fac210 * absa(ind0+11,1:ng7))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG7
                        tau_major(1:ng7) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng7) + &
                                fac100 * absa(ind0,1:ng7)   + &
                                fac000 * absa(ind0+1,1:ng7) + &
                                fac210 * absa(ind0+8,1:ng7) + &
                                fac110 * absa(ind0+9,1:ng7) + &
                                fac010 * absa(ind0+10,1:ng7))
                else
                        !CDIR EXPAND=NG7
                        tau_major(1:ng7) = speccomb *   &
                                (fac000 * absa(ind0,1:ng7)   + &
                                fac100 * absa(ind0+1,1:ng7) + &
                                fac010 * absa(ind0+9,1:ng7) + &
                                fac110 * absa(ind0+10,1:ng7))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG7
                        tau_major1(1:ng7) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng7)    + &
                                fac101 * absa(ind1+1,1:ng7)  + &
                                fac201 * absa(ind1+2,1:ng7)  + &
                                fac011 * absa(ind1+9,1:ng7)  + &
                                fac111 * absa(ind1+10,1:ng7) + &
                                fac211 * absa(ind1+11,1:ng7))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG7
                        tau_major1(1:ng7) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng7) + &
                                fac101 * absa(ind1,1:ng7)   + &
                                fac001 * absa(ind1+1,1:ng7) + &
                                fac211 * absa(ind1+8,1:ng7) + &
                                fac111 * absa(ind1+9,1:ng7) + &
                                fac011 * absa(ind1+10,1:ng7))
                else
                        !CDIR EXPAND=NG7
                        tau_major1(1:ng7) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng7)   + &
                                fac101 * absa(ind1+1,1:ng7) + &
                                fac011 * absa(ind1+9,1:ng7) + &
                                fac111 * absa(ind1+10,1:ng7))
                endif

                !CDIR EXPAND=NG7
                do ig = 1, ng7
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                        (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
                co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                        (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
                absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1)

                taug(lay,ngs6+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + adjcolco2*absco2
                fracs(lay,ngs6+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.79_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(7) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(7) + 1
                indm = indminor(lay)
                !CDIR EXPAND=NG7
                do ig = 1, ng7
                absco2 = kb_mco2(indm,ig) + minorfrac(lay) * &
                        (kb_mco2(indm+1,ig) - kb_mco2(indm,ig))
                taug(lay,ngs6+ig) = colo3(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + adjcolco2 * absco2
                fracs(lay,ngs6+ig) = fracrefb(ig)
                enddo
                enddo

                ! Empirical modification to code to improve stratospheric cooling rates
                ! for o3.  Revised to apply weighting for g-point reduction in this band.
                do lay = laytrop_max+1, nlayers
                taug(lay,ngs6+6)=taug(lay,ngs6+6)*0.92_wp
                taug(lay,ngs6+7)=taug(lay,ngs6+7)*0.88_wp
                taug(lay,ngs6+8)=taug(lay,ngs6+8)*1.07_wp
                taug(lay,ngs6+9)=taug(lay,ngs6+9)*1.1_wp
                taug(lay,ngs6+10)=taug(lay,ngs6+10)*0.99_wp
                taug(lay,ngs6+11)=taug(lay,ngs6+11)*0.855_wp
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max

                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colh2o(lay) + rat_h2oo3(lay)*colo3(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oo3_1(lay)*colo3(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mco2 = colh2o(lay) + refrat_m_a*colo3(lay)
                specparm_mco2 = MIN(colh2o(lay)/speccomb_mco2,oneminus)
                specmult_mco2 = 8._wp*specparm_mco2

                jmco2 = 1 + int(specmult_mco2)
                fmco2 = MOD1(specmult_mco2)

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 3.0_wp+(ratco2-3.0_wp)**0.79_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                speccomb_planck = colh2o(lay)+refrat_planck_a*colo3(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(7) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(7) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG7
                        tau_major(1:ng7) = speccomb *    &
                                (fac000 * absa(ind0,1:ng7)    + &
                                fac100 * absa(ind0+1,1:ng7)  + &
                                fac200 * absa(ind0+2,1:ng7)  + &
                                fac010 * absa(ind0+9,1:ng7)  + &
                                fac110 * absa(ind0+10,1:ng7) + &
                                fac210 * absa(ind0+11,1:ng7))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG7
                        tau_major(1:ng7) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng7) + &
                                fac100 * absa(ind0,1:ng7)   + &
                                fac000 * absa(ind0+1,1:ng7) + &
                                fac210 * absa(ind0+8,1:ng7) + &
                                fac110 * absa(ind0+9,1:ng7) + &
                                fac010 * absa(ind0+10,1:ng7))
                else
                        !CDIR EXPAND=NG7
                        tau_major(1:ng7) = speccomb *   &
                                (fac000 * absa(ind0,1:ng7)   + &
                                fac100 * absa(ind0+1,1:ng7) + &
                                fac010 * absa(ind0+9,1:ng7) + &
                                fac110 * absa(ind0+10,1:ng7))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG7
                        tau_major1(1:ng7) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng7)    + &
                                fac101 * absa(ind1+1,1:ng7)  + &
                                fac201 * absa(ind1+2,1:ng7)  + &
                                fac011 * absa(ind1+9,1:ng7)  + &
                                fac111 * absa(ind1+10,1:ng7) + &
                                fac211 * absa(ind1+11,1:ng7))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG7
                        tau_major1(1:ng7) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng7) + &
                                fac101 * absa(ind1,1:ng7)   + &
                                fac001 * absa(ind1+1,1:ng7) + &
                                fac211 * absa(ind1+8,1:ng7) + &
                                fac111 * absa(ind1+9,1:ng7) + &
                                fac011 * absa(ind1+10,1:ng7))
                else
                        !CDIR EXPAND=NG7
                        tau_major1(1:ng7) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng7)   + &
                                fac101 * absa(ind1+1,1:ng7) + &
                                fac011 * absa(ind1+9,1:ng7) + &
                                fac111 * absa(ind1+10,1:ng7))
                endif

                !CDIR EXPAND=NG7
                do ig = 1, ng7
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                        (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
                co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                        (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
                absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1)

                taug(lay,ngs6+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + adjcolco2*absco2
                fracs(lay,ngs6+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo

                ! Upper atmosphere part
                !CDIR NODEP,VOVERTAKE,VOB

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.79_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(7) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(7) + 1
                indm = indminor(lay)
                !CDIR EXPAND=NG7
                do ig = 1, ng7
                absco2 = kb_mco2(indm,ig) + minorfrac(lay) * &
                        (kb_mco2(indm+1,ig) - kb_mco2(indm,ig))
                taug(lay,ngs6+ig) = colo3(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + adjcolco2 * absco2
                fracs(lay,ngs6+ig) = fracrefb(ig)
                enddo

                ! Empirical modification to code to improve stratospheric cooling rates
                ! for o3.  Revised to apply weighting for g-point reduction in this band.

                !CDIR NODEP,VOVERTAKE,VOB
                taug(lay,ngs6+6)=taug(lay,ngs6+6)*0.92_wp
                taug(lay,ngs6+7)=taug(lay,ngs6+7)*0.88_wp
                taug(lay,ngs6+8)=taug(lay,ngs6+8)*1.07_wp
                taug(lay,ngs6+9)=taug(lay,ngs6+9)*1.1_wp
                taug(lay,ngs6+10)=taug(lay,ngs6+10)*0.99_wp
                taug(lay,ngs6+11)=taug(lay,ngs6+11)*0.855_wp

                enddo

        end subroutine taugb7

        !----------------------------------------------------------------------------
        subroutine taugb8
                !----------------------------------------------------------------------------
                !
                !     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
                !                             (high key - o3; high minor - co2, n2o)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                use rrlw_kg08,   only : fracrefa, fracrefb, absa, absb, &
                      ka_mco2, ka_mn2o, ka_mo3, kb_mco2, kb_mn2o, &
                      selfref, forref, cfc12, cfc22adj
                !REAL(wp) , DIMENSION(ng8) :: fracrefa
                !REAL(wp) , DIMENSION(ng8) :: fracrefb
                !REAL(wp) , DIMENSION(ng8) :: cfc12
                !REAL(wp) , DIMENSION(ng8) :: cfc22adj

                !REAL(wp) :: absa(65,ng8)
                !REAL(wp) :: absb(235,ng8)
                !REAL(wp) :: ka_mco2(19,ng8)
                !REAL(wp) :: ka_mn2o(19,ng8)
                !REAL(wp) :: ka_mo3(19,ng8)
                !REAL(wp) :: kb_mco2(19,ng8)
                !REAL(wp) :: kb_mn2o(19,ng8)
                !REAL(wp) :: selfref(10,ng8)
                !REAL(wp) :: forref(4,ng8)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng8  = 8
                INTEGER, PARAMETER :: ngs7  = 88
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                real(wp) :: tauself, taufor, absco2, abso3, absn2o
                real(wp) :: chi_co2, ratco2, adjfac, adjcolco2


                ! Minor gas mapping level:
                !     lower - co2, p = 1053.63 mb, t = 294.2 k
                !     lower - o3,  p = 317.348 mb, t = 240.77 k
                !     lower - n2o, p = 706.2720 mb, t= 278.94 k
                !     lower - cfc12,cfc11
                !     upper - co2, p = 35.1632 mb, t = 223.28 k
                !     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature, and appropriate species.  Below laytrop, the water vapor
                ! self-continuum and foreign continuum is interpolated (in temperature)
                ! separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.65_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(8) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(8) + 1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)
                !CDIR EXPAND=NG8
                do ig = 1, ng8
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                absco2 =  (ka_mco2(indm,ig) + minorfrac(lay) * &
                        (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
                abso3 =  (ka_mo3(indm,ig) + minorfrac(lay) * &
                        (ka_mo3(indm+1,ig) - ka_mo3(indm,ig)))
                absn2o =  (ka_mn2o(indm,ig) + minorfrac(lay) * &
                        (ka_mn2o(indm+1,ig) - ka_mn2o(indm,ig)))
                taug(lay,ngs7+ig) = colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) +  &
                        fac11(lay) * absa(ind1+1,ig)) &
                        + tauself + taufor &
                        + adjcolco2*absco2 &
                        + colo3(lay) * abso3 &
                        + coln2o(lay) * absn2o &
                        + wx(3,lay) * cfc12(ig) &
                        + wx(4,lay) * cfc22adj(ig)
                fracs(lay,ngs7+ig) = fracrefa(ig)
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/coldry(lay)
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.65_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1) * coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(8) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(8) + 1
                indm = indminor(lay)
                !CDIR EXPAND=NG8
                do ig = 1, ng8
                absco2 =  (kb_mco2(indm,ig) + minorfrac(lay) * &
                        (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
                absn2o =  (kb_mn2o(indm,ig) + minorfrac(lay) * &
                        (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
                taug(lay,ngs7+ig) = colo3(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + adjcolco2*absco2 &
                        + coln2o(lay)*absn2o &
                        + wx(3,lay) * cfc12(ig) &
                        + wx(4,lay) * cfc22adj(ig)
                fracs(lay,ngs7+ig) = fracrefb(ig)
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max
                !CDIR NODEP,VOVERTAKE,VOB

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.65_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(8) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(8) + 1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)
                !CDIR EXPAND=NG8
                do ig = 1, ng8
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                absco2 =  (ka_mco2(indm,ig) + minorfrac(lay) * &
                        (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
                abso3 =  (ka_mo3(indm,ig) + minorfrac(lay) * &
                        (ka_mo3(indm+1,ig) - ka_mo3(indm,ig)))
                absn2o =  (ka_mn2o(indm,ig) + minorfrac(lay) * &
                        (ka_mn2o(indm+1,ig) - ka_mn2o(indm,ig)))
                taug(lay,ngs7+ig) = colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) +  &
                        fac11(lay) * absa(ind1+1,ig)) &
                        + tauself + taufor &
                        + adjcolco2*absco2 &
                        + colo3(lay) * abso3 &
                        + coln2o(lay) * absn2o &
                        + wx(3,lay) * cfc12(ig) &
                        + wx(4,lay) * cfc22adj(ig)
                fracs(lay,ngs7+ig) = fracrefa(ig)
                enddo

                ! Upper atmosphere loop
                !CDIR NODEP,VOVERTAKE,VOB

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/coldry(lay)
                ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(lay)+1)
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.65_wp
                        adjcolco2 = adjfac*chi_mls(2,jp(lay)+1) * coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(8) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(8) + 1
                indm = indminor(lay)
                !CDIR EXPAND=NG8
                do ig = 1, ng8
                absco2 =  (kb_mco2(indm,ig) + minorfrac(lay) * &
                        (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
                absn2o =  (kb_mn2o(indm,ig) + minorfrac(lay) * &
                        (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
                taug(lay,ngs7+ig) = colo3(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + adjcolco2*absco2 &
                        + coln2o(lay)*absn2o &
                        + wx(3,lay) * cfc12(ig) &
                        + wx(4,lay) * cfc22adj(ig)
                fracs(lay,ngs7+ig) = fracrefb(ig)
                enddo
                enddo

        end subroutine taugb8

        !----------------------------------------------------------------------------
        subroutine taugb9
                !----------------------------------------------------------------------------
                !
                !     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
                !                             (high key - ch4; high minor - n2o)
                !----------------------------------------------------------------------------

                ! ------- Modules -------
                ! >> sylvia_20200320
                ! Adding back in USE statements.
                use rrlw_kg09,   only : fracrefa, fracrefb, absa, absb, &
                        ka_mn2o, kb_mn2o, selfref, forref
                !REAL(wp) , DIMENSION(ng9) :: fracrefb
                !REAL(wp) :: fracrefa(ng9,9)
                !REAL(wp) :: absa(585,ng9)
                !REAL(wp) :: absb(235,ng9)
                !REAL(wp) :: ka_mn2o(9,19,ng9)
                !REAL(wp) :: kb_mn2o(19,ng9)
                !REAL(wp) :: selfref(10,ng9)
                !REAL(wp) :: forref(4,ng9)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng9  = 12
                INTEGER, PARAMETER :: ngs8  = 96
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                integer :: js, js1, jmn2o, jpl
                real(wp) :: speccomb, specparm, specmult, fs
                real(wp) :: speccomb1, specparm1, specmult1, fs1
                real(wp) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o
                real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
                real(wp) :: p, p4, fk0, fk1, fk2
                real(wp) :: fac000, fac100, fac200
                real(wp) :: fac010, fac110, fac210
                real(wp) :: fac001, fac101, fac201
                real(wp) :: fac011, fac111, fac211
                real(wp) :: tauself, taufor, n2om1, n2om2, absn2o
                real(wp) :: chi_n2o, ratn2o, adjfac, adjcoln2o
                real(wp) :: refrat_planck_a, refrat_m_a
                real(wp) :: tau_major(ng9), tau_major1(ng9)


                ! Minor gas mapping level :
                !     lower - n2o, p = 706.272 mbar, t = 278.94 k
                !     upper - n2o, p = 95.58 mbar, t = 215.7 k

                ! Calculate reference ratio to be used in calculation of Planck
                ! fraction in lower/upper atmosphere.

                ! P = 212 mb
                refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)

                ! P = 706.272 mb
                refrat_m_a = chi_mls(1,3)/chi_mls(6,3)

                ! Compute the optical depth by interpolating in ln(pressure),
                ! temperature, and appropriate species.  Below laytrop, the water
                ! vapor self-continuum and foreign continuum is interpolated
                ! (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                speccomb = colh2o(lay) + rat_h2och4(lay)*colch4(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2och4_1(lay)*colch4(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mn2o = colh2o(lay) + refrat_m_a*colch4(lay)
                specparm_mn2o = MIN(colh2o(lay)/speccomb_mn2o,oneminus)
                specmult_mn2o = 8._wp*specparm_mn2o
                jmn2o = 1 + int(specmult_mn2o)
                fmn2o = MOD1(specmult_mn2o)

                !  In atmospheres where the amount of N2O is too great to be considered
                !  a minor species, adjust the column amount of N2O by an empirical factor
                !  to obtain the proper contribution.
                chi_n2o = coln2o(lay)/(coldry(lay))
                ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(lay)+1)
                if (ratn2o .gt. 1.5_wp) then
                        adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
                        adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcoln2o = coln2o(lay)
                endif

                speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(9) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(9) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG9
                        tau_major(1:ng9) = speccomb *    &
                                (fac000 * absa(ind0,1:ng9)    + &
                                fac100 * absa(ind0+1,1:ng9)  + &
                                fac200 * absa(ind0+2,1:ng9)  + &
                                fac010 * absa(ind0+9,1:ng9)  + &
                                fac110 * absa(ind0+10,1:ng9) + &
                                fac210 * absa(ind0+11,1:ng9))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG9
                        tau_major(1:ng9) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng9) + &
                                fac100 * absa(ind0,1:ng9)   + &
                                fac000 * absa(ind0+1,1:ng9) + &
                                fac210 * absa(ind0+8,1:ng9) + &
                                fac110 * absa(ind0+9,1:ng9) + &
                                fac010 * absa(ind0+10,1:ng9))
                else
                        !CDIR EXPAND=NG9
                        tau_major(1:ng9) = speccomb *   &
                                (fac000 * absa(ind0,1:ng9)   + &
                                fac100 * absa(ind0+1,1:ng9) + &
                                fac010 * absa(ind0+9,1:ng9) + &
                                fac110 * absa(ind0+10,1:ng9))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG9
                        tau_major1(1:ng9) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng9)    + &
                                fac101 * absa(ind1+1,1:ng9)  + &
                                fac201 * absa(ind1+2,1:ng9)  + &
                                fac011 * absa(ind1+9,1:ng9)  + &
                                fac111 * absa(ind1+10,1:ng9) + &
                                fac211 * absa(ind1+11,1:ng9))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG9
                        tau_major1(1:ng9) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng9) + &
                                fac101 * absa(ind1,1:ng9)   + &
                                fac001 * absa(ind1+1,1:ng9) + &
                                fac211 * absa(ind1+8,1:ng9) + &
                                fac111 * absa(ind1+9,1:ng9) + &
                                fac011 * absa(ind1+10,1:ng9))
                else
                        !CDIR EXPAND=NG9
                        tau_major1(1:ng9) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng9)   + &
                                fac101 * absa(ind1+1,1:ng9) + &
                                fac011 * absa(ind1+9,1:ng9) + &
                                fac111 * absa(ind1+10,1:ng9))
                endif

                !CDIR EXPAND=NG9
                do ig = 1, ng9
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                        (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
                n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                        (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
                absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)

                taug(lay,ngs8+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + adjcoln2o*absn2o
                fracs(lay,ngs8+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers

                !  In atmospheres where the amount of N2O is too great to be considered
                !  a minor species, adjust the column amount of N2O by an empirical factor
                !  to obtain the proper contribution.
                chi_n2o = coln2o(lay)/(coldry(lay))
                ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(lay)+1)
                if (ratn2o .gt. 1.5_wp) then
                        adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
                        adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcoln2o = coln2o(lay)
                endif

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(9) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(9) + 1
                indm = indminor(lay)
                !CDIR EXPAND=NG9
                do ig = 1, ng9
                absn2o = kb_mn2o(indm,ig) + minorfrac(lay) * &
                        (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
                taug(lay,ngs8+ig) = colch4(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) +  &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + adjcoln2o*absn2o
                fracs(lay,ngs8+ig) = fracrefb(ig)
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max
                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colh2o(lay) + rat_h2och4(lay)*colch4(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2och4_1(lay)*colch4(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mn2o = colh2o(lay) + refrat_m_a*colch4(lay)
                specparm_mn2o = MIN(colh2o(lay)/speccomb_mn2o,oneminus)
                specmult_mn2o = 8._wp*specparm_mn2o
                jmn2o = 1 + int(specmult_mn2o)
                fmn2o = MOD1(specmult_mn2o)

                !  In atmospheres where the amount of N2O is too great to be considered
                !  a minor species, adjust the column amount of N2O by an empirical factor
                !  to obtain the proper contribution.
                chi_n2o = coln2o(lay)/(coldry(lay))
                ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(lay)+1)
                if (ratn2o .gt. 1.5_wp) then
                        adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
                        adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcoln2o = coln2o(lay)
                endif

                speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl= 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(9) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(9) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG9
                        tau_major(1:ng9) = speccomb *    &
                                (fac000 * absa(ind0,1:ng9)    + &
                                fac100 * absa(ind0+1,1:ng9)  + &
                                fac200 * absa(ind0+2,1:ng9)  + &
                                fac010 * absa(ind0+9,1:ng9)  + &
                                fac110 * absa(ind0+10,1:ng9) + &
                                fac210 * absa(ind0+11,1:ng9))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG9
                        tau_major(1:ng9) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng9) + &
                                fac100 * absa(ind0,1:ng9)   + &
                                fac000 * absa(ind0+1,1:ng9) + &
                                fac210 * absa(ind0+8,1:ng9) + &
                                fac110 * absa(ind0+9,1:ng9) + &
                                fac010 * absa(ind0+10,1:ng9))
                else
                        !CDIR EXPAND=NG9
                        tau_major(1:ng9) = speccomb *   &
                                (fac000 * absa(ind0,1:ng9)   + &
                                fac100 * absa(ind0+1,1:ng9) + &
                                fac010 * absa(ind0+9,1:ng9) + &
                                fac110 * absa(ind0+10,1:ng9))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG9
                        tau_major1(1:ng9) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng9)    + &
                                fac101 * absa(ind1+1,1:ng9)  + &
                                fac201 * absa(ind1+2,1:ng9)  + &
                                fac011 * absa(ind1+9,1:ng9)  + &
                                fac111 * absa(ind1+10,1:ng9) + &
                                fac211 * absa(ind1+11,1:ng9))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG9
                        tau_major1(1:ng9) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng9) + &
                                fac101 * absa(ind1,1:ng9)   + &
                                fac001 * absa(ind1+1,1:ng9) + &
                                fac211 * absa(ind1+8,1:ng9) + &
                                fac111 * absa(ind1+9,1:ng9) + &
                                fac011 * absa(ind1+10,1:ng9))
                else
                        !CDIR EXPAND=NG9
                        tau_major1(1:ng9) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng9)   + &
                                fac101 * absa(ind1+1,1:ng9) + &
                                fac011 * absa(ind1+9,1:ng9) + &
                                fac111 * absa(ind1+10,1:ng9))
                endif

                !CDIR EXPAND=NG9
                do ig = 1, ng9
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                        (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
                n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                        (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
                absn2o = n2om1 + minorfrac(lay) * (n2om2 - n2om1)

                taug(lay,ngs8+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + adjcoln2o*absn2o
                fracs(lay,ngs8+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo

                ! Upper atmosphere part
                !CDIR NODEP,VOVERTAKE,VOB

                !  In atmospheres where the amount of N2O is too great to be considered
                !  a minor species, adjust the column amount of N2O by an empirical factor
                !  to obtain the proper contribution.
                chi_n2o = coln2o(lay)/(coldry(lay))
                ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(lay)+1)
                if (ratn2o .gt. 1.5_wp) then
                        adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
                        adjcoln2o = adjfac*chi_mls(4,jp(lay)+1)*coldry(lay)*1.e-20_wp
                else
                        adjcoln2o = coln2o(lay)
                endif

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(9) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(9) + 1
                indm = indminor(lay)
                !CDIR EXPAND=NG9
                do ig = 1, ng9
                absn2o = kb_mn2o(indm,ig) + minorfrac(lay) * &
                        (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
                taug(lay,ngs8+ig) = colch4(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) +  &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + adjcoln2o*absn2o
                fracs(lay,ngs8+ig) = fracrefb(ig)
                enddo
                enddo

        end subroutine taugb9

        !----------------------------------------------------------------------------
        subroutine taugb10
                !----------------------------------------------------------------------------
                !
                !     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! Adding in USE statements
                use rrlw_kg10,   only : fracrefa, fracrefb, absa, absb, &
                          selfref, forref
                !REAL(wp) , DIMENSION(ng10) :: fracrefa
                !REAL(wp) , DIMENSION(ng10) :: fracrefb

                !REAL(wp) :: absa(65,ng10)
                !REAL(wp) :: absb(235,ng10)
                !REAL(wp) :: selfref(10,ng10)
                !REAL(wp) :: forref(4,ng10)
                ! << sylvia_20200310

                ! ------- Declarations -------
                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng10 = 6
                INTEGER, PARAMETER :: ngs9 = 108

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf
                real(wp) :: tauself, taufor


                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature.  Below laytrop, the water vapor self-continuum and
                ! foreign continuum is interpolated (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(10) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(10) + 1
                inds = indself(lay)
                indf = indfor(lay)
                !CDIR EXPAND=NG10
                do ig = 1, ng10
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs9+ig) = colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) + &
                        fac11(lay) * absa(ind1+1,ig))  &
                        + tauself + taufor
                fracs(lay,ngs9+ig) = fracrefa(ig)
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(10) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(10) + 1
                indf = indfor(lay)
                !CDIR EXPAND=NG10
                do ig = 1, ng10
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs9+ig) = colh2o(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) +  &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + taufor
                fracs(lay,ngs9+ig) = fracrefb(ig)
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(10) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(10) + 1
                inds = indself(lay)
                indf = indfor(lay)
                !CDIR EXPAND=NG10
                do ig = 1, ng10
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs9+ig) = colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) + &
                        fac11(lay) * absa(ind1+1,ig))  &
                        + tauself + taufor
                fracs(lay,ngs9+ig) = fracrefa(ig)
                enddo

                ! Upper atmosphere part
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(10) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(10) + 1
                indf = indfor(lay)
                !CDIR EXPAND=NG10
                do ig = 1, ng10
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs9+ig) = colh2o(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) +  &
                        fac11(lay) * absb(ind1+1,ig)) &
                        + taufor
                fracs(lay,ngs9+ig) = fracrefb(ig)
                enddo
                enddo

        end subroutine taugb10

        !----------------------------------------------------------------------------
        subroutine taugb11
                !----------------------------------------------------------------------------
                !
                !     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
                !                              (high key - h2o; high minor - o2)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Add back in USE statements.
                use rrlw_kg11,   only : fracrefa, fracrefb, absa, absb, &
                         ka_mo2, kb_mo2, selfref, forref
                !REAL(wp) , DIMENSION(ng11) :: fracrefa
                !REAL(wp) , DIMENSION(ng11) :: fracrefb

                !REAL(wp) :: absa(65,ng11)
                !REAL(wp) :: absb(235,ng11)
                !REAL(wp) :: ka_mo2(19,ng11)
                !REAL(wp) :: kb_mo2(19,ng11)
                !REAL(wp) :: selfref(10,ng11)
                !REAL(wp) :: forref(4,ng11)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng11 = 8
                INTEGER, PARAMETER :: ngs10 = 114
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                real(wp) :: scaleo2, tauself, taufor, tauo2


                ! Minor gas mapping level :
                !     lower - o2, p = 706.2720 mbar, t = 278.94 k
                !     upper - o2, p = 4.758820 mbarm t = 250.85 k

                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature.  Below laytrop, the water vapor self-continuum and
                ! foreign continuum is interpolated (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(11) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(11) + 1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)
                scaleo2 = colo2(lay)*scaleminor(lay)
                !CDIR EXPAND=NG11
                do ig = 1, ng11
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                tauo2 =  scaleo2 * (ka_mo2(indm,ig) + minorfrac(lay) * &
                        (ka_mo2(indm+1,ig) - ka_mo2(indm,ig)))
                taug(lay,ngs10+ig) = colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) + &
                        fac11(lay) * absa(ind1+1,ig)) &
                        + tauself + taufor &
                        + tauo2
                fracs(lay,ngs10+ig) = fracrefa(ig)
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(11) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(11) + 1
                indf = indfor(lay)
                indm = indminor(lay)
                scaleo2 = colo2(lay)*scaleminor(lay)
                !CDIR EXPAND=NG11
                do ig = 1, ng11
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                tauo2 =  scaleo2 * (kb_mo2(indm,ig) + minorfrac(lay) * &
                        (kb_mo2(indm+1,ig) - kb_mo2(indm,ig)))
                taug(lay,ngs10+ig) = colh2o(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig))  &
                        + taufor &
                        + tauo2
                fracs(lay,ngs10+ig) = fracrefb(ig)
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(11) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(11) + 1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)
                scaleo2 = colo2(lay)*scaleminor(lay)
                !CDIR EXPAND=NG11
                do ig = 1, ng11
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                tauo2 =  scaleo2 * (ka_mo2(indm,ig) + minorfrac(lay) * &
                        (ka_mo2(indm+1,ig) - ka_mo2(indm,ig)))
                taug(lay,ngs10+ig) = colh2o(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) + &
                        fac11(lay) * absa(ind1+1,ig)) &
                        + tauself + taufor &
                        + tauo2
                fracs(lay,ngs10+ig) = fracrefa(ig)
                enddo

                ! Upper atmosphere part
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(11) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(11) + 1
                indf = indfor(lay)
                indm = indminor(lay)
                scaleo2 = colo2(lay)*scaleminor(lay)
                !CDIR EXPAND=NG11
                do ig = 1, ng11
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                tauo2 =  scaleo2 * (kb_mo2(indm,ig) + minorfrac(lay) * &
                        (kb_mo2(indm+1,ig) - kb_mo2(indm,ig)))
                taug(lay,ngs10+ig) = colh2o(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig))  &
                        + taufor &
                        + tauo2
                fracs(lay,ngs10+ig) = fracrefb(ig)
                enddo
                enddo

        end subroutine taugb11

        !----------------------------------------------------------------------------
        subroutine taugb12
                !----------------------------------------------------------------------------
                !
                !     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Add in USE statements again.
                use rrlw_kg12,   only : fracrefa, absa, &
                       selfref, forref
                !REAL(wp) :: fracrefa(ng12,9)
                !REAL(wp) :: absa(585,ng12)
                !REAL(wp) :: selfref(10,ng12)
                !REAL(wp) :: forref(4,ng12)
                ! << sylvia_20200310

                ! ------- Declarations -------
                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng12 = 8
                INTEGER, PARAMETER :: ngs11 = 122
                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf
                integer :: js, js1, jpl
                real(wp) :: speccomb, specparm, specmult, fs
                real(wp) :: speccomb1, specparm1, specmult1, fs1
                real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
                real(wp) :: p, p4, fk0, fk1, fk2
                real(wp) :: fac000, fac100, fac200
                real(wp) :: fac010, fac110, fac210
                real(wp) :: fac001, fac101, fac201
                real(wp) :: fac011, fac111, fac211
                real(wp) :: tauself, taufor
                real(wp) :: refrat_planck_a
                real(wp) :: tau_major(ng12), tau_major1(ng12)


                ! Calculate reference ratio to be used in calculation of Planck
                ! fraction in lower/upper atmosphere.

                ! P =   174.164 mb
                refrat_planck_a = chi_mls(1,10)/chi_mls(2,10)

                ! Compute the optical depth by interpolating in ln(pressure),
                ! temperature, and appropriate species.  Below laytrop, the water
                ! vapor self-continuum adn foreign continuum is interpolated
                ! (in temperature) separately.

                ! Lower atmosphere loop
                DO lay = 1, laytrop_min

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(12) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(12) + js1
                inds = indself(lay)
                indf = indfor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG12
                        tau_major(1:ng12) = speccomb *    &
                                (fac000 * absa(ind0,1:ng12)    + &
                                fac100 * absa(ind0+1,1:ng12)  + &
                                fac200 * absa(ind0+2,1:ng12)  + &
                                fac010 * absa(ind0+9,1:ng12)  + &
                                fac110 * absa(ind0+10,1:ng12) + &
                                fac210 * absa(ind0+11,1:ng12))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG12
                        tau_major(1:ng12) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng12) + &
                                fac100 * absa(ind0,1:ng12)   + &
                                fac000 * absa(ind0+1,1:ng12) + &
                                fac210 * absa(ind0+8,1:ng12) + &
                                fac110 * absa(ind0+9,1:ng12) + &
                                fac010 * absa(ind0+10,1:ng12))
                else
                        !CDIR EXPAND=NG12
                        tau_major(1:ng12) = speccomb *   &
                                (fac000 * absa(ind0,1:ng12)   + &
                                fac100 * absa(ind0+1,1:ng12) + &
                                fac010 * absa(ind0+9,1:ng12) + &
                                fac110 * absa(ind0+10,1:ng12))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG12
                        tau_major1(1:ng12) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng12)    + &
                                fac101 * absa(ind1+1,1:ng12)  + &
                                fac201 * absa(ind1+2,1:ng12)  + &
                                fac011 * absa(ind1+9,1:ng12)  + &
                                fac111 * absa(ind1+10,1:ng12) + &
                                fac211 * absa(ind1+11,1:ng12))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG12
                        tau_major1(1:ng12) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng12) + &
                                fac101 * absa(ind1,1:ng12)   + &
                                fac001 * absa(ind1+1,1:ng12) + &
                                fac211 * absa(ind1+8,1:ng12) + &
                                fac111 * absa(ind1+9,1:ng12) + &
                                fac011 * absa(ind1+10,1:ng12))
                else
                        !CDIR EXPAND=NG12
                        tau_major1(1:ng12) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng12)   + &
                                fac101 * absa(ind1+1,1:ng12) + &
                                fac011 * absa(ind1+9,1:ng12) + &
                                fac111 * absa(ind1+10,1:ng12))
                endif

                !CDIR EXPAND=NG12
                do ig = 1, ng12
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))

                taug(lay,ngs11+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor
                fracs(lay,ngs11+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                ENDDO

                ! Upper atmosphere loop
                do ig = 1, ng12
                do lay = laytrop_max+1, nlayers
                taug(lay,ngs11+ig) = 0.0_wp
                fracs(lay,ngs11+ig) = 0.0_wp
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max
                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colh2o(lay) + rat_h2oco2(lay)*colco2(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_planck = colh2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(12) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(12) + js1
                inds = indself(lay)
                indf = indfor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG12
                        tau_major(1:ng12) = speccomb *    &
                                (fac000 * absa(ind0,1:ng12)    + &
                                fac100 * absa(ind0+1,1:ng12)  + &
                                fac200 * absa(ind0+2,1:ng12)  + &
                                fac010 * absa(ind0+9,1:ng12)  + &
                                fac110 * absa(ind0+10,1:ng12) + &
                                fac210 * absa(ind0+11,1:ng12))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG12
                        tau_major(1:ng12) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng12) + &
                                fac100 * absa(ind0,1:ng12)   + &
                                fac000 * absa(ind0+1,1:ng12) + &
                                fac210 * absa(ind0+8,1:ng12) + &
                                fac110 * absa(ind0+9,1:ng12) + &
                                fac010 * absa(ind0+10,1:ng12))
                else
                        !CDIR EXPAND=NG12
                        tau_major(1:ng12) = speccomb *   &
                                (fac000 * absa(ind0,1:ng12)   + &
                                fac100 * absa(ind0+1,1:ng12) + &
                                fac010 * absa(ind0+9,1:ng12) + &
                                fac110 * absa(ind0+10,1:ng12))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG12
                        tau_major1(1:ng12) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng12)    + &
                                fac101 * absa(ind1+1,1:ng12)  + &
                                fac201 * absa(ind1+2,1:ng12)  + &
                                fac011 * absa(ind1+9,1:ng12)  + &
                                fac111 * absa(ind1+10,1:ng12) + &
                                fac211 * absa(ind1+11,1:ng12))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG12
                        tau_major1(1:ng12) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng12) + &
                                fac101 * absa(ind1,1:ng12)   + &
                                fac001 * absa(ind1+1,1:ng12) + &
                                fac211 * absa(ind1+8,1:ng12) + &
                                fac111 * absa(ind1+9,1:ng12) + &
                                fac011 * absa(ind1+10,1:ng12))
                else
                        !CDIR EXPAND=NG12
                        tau_major1(1:ng12) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng12)   + &
                                fac101 * absa(ind1+1,1:ng12) + &
                                fac011 * absa(ind1+9,1:ng12) + &
                                fac111 * absa(ind1+10,1:ng12))
                endif

                !CDIR EXPAND=NG12
                do ig = 1, ng12
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))

                taug(lay,ngs11+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor
                fracs(lay,ngs11+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo

                ! Upper atmosphere part

                do ig = 1, ng12
                !CDIR NODEP,VOVERTAKE,VOB

                taug(lay,ngs11+ig) = 0.0_wp
                fracs(lay,ngs11+ig) = 0.0_wp
                enddo
                enddo

        end subroutine taugb12

        !----------------------------------------------------------------------------
        subroutine taugb13
                !----------------------------------------------------------------------------
                !
                !     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Add back in USE statement.
                use rrlw_kg13,   only : fracrefa, fracrefb, absa, &
                       ka_mco2, ka_mco, kb_mo3, selfref, forref
                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng13 = 4
                INTEGER, PARAMETER :: ngs12 = 130
                !REAL(wp) , DIMENSION(ng13) :: fracrefb

                !REAL(wp) :: fracrefa(ng13,9)
                !REAL(wp) :: absa(585,ng13)
                !REAL(wp) :: ka_mco2(9,19,ng13)
                !REAL(wp) :: ka_mco(9,19,ng13)
                !REAL(wp) :: kb_mo3(19,ng13)
                !REAL(wp) :: selfref(10,ng13)
                !REAL(wp) :: forref(4,ng13)
                ! << sylvia_20200310

                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                integer :: js, js1, jmco2, jmco, jpl
                real(wp) :: speccomb, specparm, specmult, fs
                real(wp) :: speccomb1, specparm1, specmult1, fs1
                real(wp) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
                real(wp) :: speccomb_mco, specparm_mco, specmult_mco, fmco
                real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
                real(wp) :: p, p4, fk0, fk1, fk2
                real(wp) :: fac000, fac100, fac200
                real(wp) :: fac010, fac110, fac210
                real(wp) :: fac001, fac101, fac201
                real(wp) :: fac011, fac111, fac211
                real(wp) :: tauself, taufor, co2m1, co2m2, absco2
                real(wp) :: com1, com2, absco, abso3
                real(wp) :: chi_co2, ratco2, adjfac, adjcolco2
                real(wp) :: refrat_planck_a, refrat_m_a, refrat_m_a3
                real(wp) :: tau_major(ng13), tau_major1(ng13)


                ! Minor gas mapping levels :
                !     lower - co2, p = 1053.63 mb, t = 294.2 k
                !     lower - co, p = 706 mb, t = 278.94 k
                !     upper - o3, p = 95.5835 mb, t = 215.7 k

                ! Calculate reference ratio to be used in calculation of Planck
                ! fraction in lower/upper atmosphere.

                ! P = 473.420 mb (Level 5)
                refrat_planck_a = chi_mls(1,5)/chi_mls(4,5)

                ! P = 1053. (Level 1)
                refrat_m_a = chi_mls(1,1)/chi_mls(4,1)

                ! P = 706. (Level 3)
                refrat_m_a3 = chi_mls(1,3)/chi_mls(4,3)

                ! Compute the optical depth by interpolating in ln(pressure),
                ! temperature, and appropriate species.  Below laytrop, the water
                ! vapor self-continuum and foreign continuum is interpolated
                ! (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                speccomb = colh2o(lay) + rat_h2on2o(lay)*coln2o(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2on2o_1(lay)*coln2o(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mco2 = colh2o(lay) + refrat_m_a*coln2o(lay)
                specparm_mco2 = MIN(colh2o(lay)/speccomb_mco2,oneminus)
                specmult_mco2 = 8._wp*specparm_mco2
                jmco2 = 1 + int(specmult_mco2)
                fmco2 = MOD1(specmult_mco2)

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/3.55e-4_wp
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.68_wp
                        adjcolco2 = adjfac*3.55e-4_wp*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                speccomb_mco = colh2o(lay) + refrat_m_a3*coln2o(lay)
                specparm_mco = MIN(colh2o(lay)/speccomb_mco,oneminus)
                specmult_mco = 8._wp*specparm_mco
                jmco = 1 + int(specmult_mco)
                fmco = MOD1(specmult_mco)

                speccomb_planck = colh2o(lay)+refrat_planck_a*coln2o(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(13) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(13) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG13
                        tau_major(1:ng13) = speccomb *    &
                                (fac000 * absa(ind0,1:ng13)    + &
                                fac100 * absa(ind0+1,1:ng13)  + &
                                fac200 * absa(ind0+2,1:ng13)  + &
                                fac010 * absa(ind0+9,1:ng13)  + &
                                fac110 * absa(ind0+10,1:ng13) + &
                                fac210 * absa(ind0+11,1:ng13))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG13
                        tau_major(1:ng13) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng13) + &
                                fac100 * absa(ind0,1:ng13)   + &
                                fac000 * absa(ind0+1,1:ng13) + &
                                fac210 * absa(ind0+8,1:ng13) + &
                                fac110 * absa(ind0+9,1:ng13) + &
                                fac010 * absa(ind0+10,1:ng13))
                else
                        !CDIR EXPAND=NG13
                        tau_major(1:ng13) = speccomb *   &
                                (fac000 * absa(ind0,1:ng13)   + &
                                fac100 * absa(ind0+1,1:ng13) + &
                                fac010 * absa(ind0+9,1:ng13) + &
                                fac110 * absa(ind0+10,1:ng13))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG13
                        tau_major1(1:ng13) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng13)    + &
                                fac101 * absa(ind1+1,1:ng13)  + &
                                fac201 * absa(ind1+2,1:ng13)  + &
                                fac011 * absa(ind1+9,1:ng13)  + &
                                fac111 * absa(ind1+10,1:ng13) + &
                                fac211 * absa(ind1+11,1:ng13))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG13
                        tau_major1(1:ng13) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng13) + &
                                fac101 * absa(ind1,1:ng13)   + &
                                fac001 * absa(ind1+1,1:ng13) + &
                                fac211 * absa(ind1+8,1:ng13) + &
                                fac111 * absa(ind1+9,1:ng13) + &
                                fac011 * absa(ind1+10,1:ng13))
                else
                        !CDIR EXPAND=NG13
                        tau_major1(1:ng13) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng13)   + &
                                fac101 * absa(ind1+1,1:ng13) + &
                                fac011 * absa(ind1+9,1:ng13) + &
                                fac111 * absa(ind1+10,1:ng13))
                endif

                !CDIR EXPAND=NG13
                do ig = 1, ng13
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                        (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
                co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                        (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
                absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1)
                com1 = ka_mco(jmco,indm,ig) + fmco * &
                        (ka_mco(jmco+1,indm,ig) - ka_mco(jmco,indm,ig))
                com2 = ka_mco(jmco,indm+1,ig) + fmco * &
                        (ka_mco(jmco+1,indm+1,ig) - ka_mco(jmco,indm+1,ig))
                absco = com1 + minorfrac(lay) * (com2 - com1)

                taug(lay,ngs12+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + adjcolco2*absco2 &
                        + colco(lay)*absco
                fracs(lay,ngs12+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers
                indm = indminor(lay)
                !CDIR EXPAND=NG13
                do ig = 1, ng13
                abso3 = kb_mo3(indm,ig) + minorfrac(lay) * &
                        (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
                taug(lay,ngs12+ig) = colo3(lay)*abso3
                fracs(lay,ngs12+ig) =  fracrefb(ig)
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max

                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = colh2o(lay) + rat_h2on2o(lay)*coln2o(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2on2o_1(lay)*coln2o(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mco2 = colh2o(lay) + refrat_m_a*coln2o(lay)
                specparm_mco2 = MIN(colh2o(lay)/speccomb_mco2,oneminus)
                specmult_mco2 = 8._wp*specparm_mco2
                jmco2 = 1 + int(specmult_mco2)
                fmco2 = MOD1(specmult_mco2)

                !  In atmospheres where the amount of CO2 is too great to be considered
                !  a minor species, adjust the column amount of CO2 by an empirical factor
                !  to obtain the proper contribution.
                chi_co2 = colco2(lay)/(coldry(lay))
                ratco2 = 1.e20_wp*chi_co2/3.55e-4_wp
                if (ratco2 .gt. 3.0_wp) then
                        adjfac = 2.0_wp+(ratco2-2.0_wp)**0.68_wp
                        adjcolco2 = adjfac*3.55e-4_wp*coldry(lay)*1.e-20_wp
                else
                        adjcolco2 = colco2(lay)
                endif

                speccomb_mco = colh2o(lay) + refrat_m_a3*coln2o(lay)
                specparm_mco = MIN(colh2o(lay)/speccomb_mco,oneminus)
                specmult_mco = 8._wp*specparm_mco
                jmco = 1 + int(specmult_mco)
                fmco = MOD1(specmult_mco)

                speccomb_planck = colh2o(lay)+refrat_planck_a*coln2o(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(13) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(13) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG13
                        tau_major(1:ng13) = speccomb *    &
                                (fac000 * absa(ind0,1:ng13)    + &
                                fac100 * absa(ind0+1,1:ng13)  + &
                                fac200 * absa(ind0+2,1:ng13)  + &
                                fac010 * absa(ind0+9,1:ng13)  + &
                                fac110 * absa(ind0+10,1:ng13) + &
                                fac210 * absa(ind0+11,1:ng13))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG13
                        tau_major(1:ng13) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng13) + &
                                fac100 * absa(ind0,1:ng13)   + &
                                fac000 * absa(ind0+1,1:ng13) + &
                                fac210 * absa(ind0+8,1:ng13) + &
                                fac110 * absa(ind0+9,1:ng13) + &
                                fac010 * absa(ind0+10,1:ng13))
                else
                        !CDIR EXPAND=NG13
                        tau_major(1:ng13) = speccomb *   &
                                (fac000 * absa(ind0,1:ng13)   + &
                                fac100 * absa(ind0+1,1:ng13) + &
                                fac010 * absa(ind0+9,1:ng13) + &
                                fac110 * absa(ind0+10,1:ng13))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG13
                        tau_major1(1:ng13) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng13)    + &
                                fac101 * absa(ind1+1,1:ng13)  + &
                                fac201 * absa(ind1+2,1:ng13)  + &
                                fac011 * absa(ind1+9,1:ng13)  + &
                                fac111 * absa(ind1+10,1:ng13) + &
                                fac211 * absa(ind1+11,1:ng13))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG13
                        tau_major1(1:ng13) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng13) + &
                                fac101 * absa(ind1,1:ng13)   + &
                                fac001 * absa(ind1+1,1:ng13) + &
                                fac211 * absa(ind1+8,1:ng13) + &
                                fac111 * absa(ind1+9,1:ng13) + &
                                fac011 * absa(ind1+10,1:ng13))
                else
                        !CDIR EXPAND=NG13
                        tau_major1(1:ng13) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng13)   + &
                                fac101 * absa(ind1+1,1:ng13) + &
                                fac011 * absa(ind1+9,1:ng13) + &
                                fac111 * absa(ind1+10,1:ng13))
                endif

                !CDIR EXPAND=NG13
                do ig = 1, ng13
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor = forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                        (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
                co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                        (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
                absco2 = co2m1 + minorfrac(lay) * (co2m2 - co2m1)
                com1 = ka_mco(jmco,indm,ig) + fmco * &
                        (ka_mco(jmco+1,indm,ig) - ka_mco(jmco,indm,ig))
                com2 = ka_mco(jmco,indm+1,ig) + fmco * &
                        (ka_mco(jmco+1,indm+1,ig) - ka_mco(jmco,indm+1,ig))
                absco = com1 + minorfrac(lay) * (com2 - com1)

                taug(lay,ngs12+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + adjcolco2*absco2 &
                        + colco(lay)*absco
                fracs(lay,ngs12+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo

                ! Upper atmosphere part
                !CDIR NODEP,VOVERTAKE,VOB

                indm = indminor(lay)
                !CDIR EXPAND=NG13
                do ig = 1, ng13
                abso3 = kb_mo3(indm,ig) + minorfrac(lay) * &
                        (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
                taug(lay,ngs12+ig) = colo3(lay)*abso3
                fracs(lay,ngs12+ig) =  fracrefb(ig)
                enddo
                enddo

        end subroutine taugb13

        !----------------------------------------------------------------------------
        subroutine taugb14
                !----------------------------------------------------------------------------
                !
                !     band 14:  2250-2380 cm-1 (low - co2; high - co2)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Add back in USE statements.
                use rrlw_kg14,   only : fracrefa, fracrefb, absa, absb, &
                     selfref, forref
                !REAL(wp) , DIMENSION(ng14) :: fracrefa
                !REAL(wp) , DIMENSION(ng14) :: fracrefb

                !REAL(wp) :: absa(65,ng14)
                !REAL(wp) :: absb(235,ng14)
                !REAL(wp) :: selfref(10,ng14)
                !REAL(wp) :: forref(4,ng14)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng14 = 2
                INTEGER, PARAMETER :: ngs13 = 134
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf
                real(wp) :: tauself, taufor


                ! Compute the optical depth by interpolating in ln(pressure) and
                ! temperature.  Below laytrop, the water vapor self-continuum
                ! and foreign continuum is interpolated (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(14) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(14) + 1
                inds = indself(lay)
                indf = indfor(lay)
                !CDIR EXPAND=NG14
                do ig = 1, ng14
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs13+ig) = colco2(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) + &
                        fac11(lay) * absa(ind1+1,ig)) &
                        + tauself + taufor
                fracs(lay,ngs13+ig) = fracrefa(ig)
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers
                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(14) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(14) + 1
                !CDIR EXPAND=NG14
                do ig = 1, ng14
                taug(lay,ngs13+ig) = colco2(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig))
                fracs(lay,ngs13+ig) = fracrefb(ig)
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(14) + 1
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(14) + 1
                inds = indself(lay)
                indf = indfor(lay)
                !CDIR EXPAND=NG14
                do ig = 1, ng14
                tauself = selffac(lay) * (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                taug(lay,ngs13+ig) = colco2(lay) * &
                        (fac00(lay) * absa(ind0,ig) + &
                        fac10(lay) * absa(ind0+1,ig) + &
                        fac01(lay) * absa(ind1,ig) + &
                        fac11(lay) * absa(ind1+1,ig)) &
                        + tauself + taufor
                fracs(lay,ngs13+ig) = fracrefa(ig)
                enddo

                ! Upper atmosphere part
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(14) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(14) + 1
                !CDIR EXPAND=NG14
                do ig = 1, ng14
                taug(lay,ngs13+ig) = colco2(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig))
                fracs(lay,ngs13+ig) = fracrefb(ig)
                enddo
                enddo

        end subroutine taugb14

        !----------------------------------------------------------------------------
        subroutine taugb15
                !----------------------------------------------------------------------------
                !
                !     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
                !                              (high - nothing)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Add in USE statements again.
                use rrlw_kg15,   only : fracrefa, absa, &
                       ka_mn2, selfref, forref
                !REAL(wp) :: fracrefa(ng15,9)
                !REAL(wp) :: absa(585,ng15)
                !REAL(wp) :: ka_mn2(9,19,ng15)
                !REAL(wp) :: selfref(10,ng15)
                !REAL(wp) :: forref(4,ng15)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng15 = 2
                INTEGER, PARAMETER :: ngs14 = 136

                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf, indm
                integer :: js, js1, jmn2, jpl
                real(wp) :: speccomb, specparm, specmult, fs
                real(wp) :: speccomb1, specparm1, specmult1, fs1
                real(wp) :: speccomb_mn2, specparm_mn2, specmult_mn2, fmn2
                real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
                real(wp) :: p, p4, fk0, fk1, fk2
                real(wp) :: fac000, fac100, fac200
                real(wp) :: fac010, fac110, fac210
                real(wp) :: fac001, fac101, fac201
                real(wp) :: fac011, fac111, fac211
                real(wp) :: scalen2, tauself, taufor, n2m1, n2m2, taun2
                real(wp) :: refrat_planck_a, refrat_m_a
                real(wp) :: tau_major(ng15), tau_major1(ng15)


                ! Minor gas mapping level :
                !     Lower - Nitrogen Continuum, P = 1053., T = 294.

                ! Calculate reference ratio to be used in calculation of Planck
                ! fraction in lower atmosphere.
                ! P = 1053. mb (Level 1)
                refrat_planck_a = chi_mls(4,1)/chi_mls(2,1)

                ! P = 1053.
                refrat_m_a = chi_mls(4,1)/chi_mls(2,1)

                ! Compute the optical depth by interpolating in ln(pressure),
                ! temperature, and appropriate species.  Below laytrop, the water
                ! vapor self-continuum and foreign continuum is interpolated
                ! (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                speccomb = coln2o(lay) + rat_n2oco2(lay)*colco2(lay)
                specparm = MIN(coln2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = coln2o(lay) + rat_n2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(coln2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mn2 = coln2o(lay) + refrat_m_a*colco2(lay)
                specparm_mn2 = MIN(coln2o(lay)/speccomb_mn2,oneminus)
                specmult_mn2 = 8._wp*specparm_mn2
                jmn2 = 1 + int(specmult_mn2)
                fmn2 = MOD1(specmult_mn2)

                speccomb_planck = coln2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(coln2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(15) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(15) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                scalen2 = colbrd(lay)*scaleminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG15
                        tau_major(1:ng15) = speccomb *    &
                                (fac000 * absa(ind0,1:ng15)    + &
                                fac100 * absa(ind0+1,1:ng15)  + &
                                fac200 * absa(ind0+2,1:ng15)  + &
                                fac010 * absa(ind0+9,1:ng15)  + &
                                fac110 * absa(ind0+10,1:ng15) + &
                                fac210 * absa(ind0+11,1:ng15))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG15
                        tau_major(1:ng15) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng15) + &
                                fac100 * absa(ind0,1:ng15)   + &
                                fac000 * absa(ind0+1,1:ng15) + &
                                fac210 * absa(ind0+8,1:ng15) + &
                                fac110 * absa(ind0+9,1:ng15) + &
                                fac010 * absa(ind0+10,1:ng15))
                else
                        !CDIR EXPAND=NG15
                        tau_major(1:ng15) = speccomb *   &
                                (fac000 * absa(ind0,1:ng15)   + &
                                fac100 * absa(ind0+1,1:ng15) + &
                                fac010 * absa(ind0+9,1:ng15) + &
                                fac110 * absa(ind0+10,1:ng15))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG15
                        tau_major1(1:ng15) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng15)    + &
                                fac101 * absa(ind1+1,1:ng15)  + &
                                fac201 * absa(ind1+2,1:ng15)  + &
                                fac011 * absa(ind1+9,1:ng15)  + &
                                fac111 * absa(ind1+10,1:ng15) + &
                                fac211 * absa(ind1+11,1:ng15))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG15
                        tau_major1(1:ng15) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng15) + &
                                fac101 * absa(ind1,1:ng15)   + &
                                fac001 * absa(ind1+1,1:ng15) + &
                                fac211 * absa(ind1+8,1:ng15) + &
                                fac111 * absa(ind1+9,1:ng15) + &
                                fac011 * absa(ind1+10,1:ng15))
                else
                        !CDIR EXPAND=NG15
                        tau_major1(1:ng15) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng15)   + &
                                fac101 * absa(ind1+1,1:ng15) + &
                                fac011 * absa(ind1+9,1:ng15) + &
                                fac111 * absa(ind1+10,1:ng15))
                endif

                !CDIR EXPAND=NG15
                do ig = 1, ng15
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                n2m1 = ka_mn2(jmn2,indm,ig) + fmn2 * &
                        (ka_mn2(jmn2+1,indm,ig) - ka_mn2(jmn2,indm,ig))
                n2m2 = ka_mn2(jmn2,indm+1,ig) + fmn2 * &
                        (ka_mn2(jmn2+1,indm+1,ig) - ka_mn2(jmn2,indm+1,ig))
                taun2 = scalen2 * (n2m1 + minorfrac(lay) * (n2m2 - n2m1))

                taug(lay,ngs14+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + taun2
                fracs(lay,ngs14+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                enddo

                ! Upper atmosphere loop
                do ig = 1, ng15
                do lay = laytrop_max+1, nlayers
                taug(lay,ngs14+ig) = 0.0_wp
                fracs(lay,ngs14+ig) = 0.0_wp
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max

                !CDIR NODEP,VOVERTAKE,VOB

                speccomb = coln2o(lay) + rat_n2oco2(lay)*colco2(lay)
                specparm = MIN(coln2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = coln2o(lay) + rat_n2oco2_1(lay)*colco2(lay)
                specparm1 = MIN(coln2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_mn2 = coln2o(lay) + refrat_m_a*colco2(lay)
                specparm_mn2 = MIN(coln2o(lay)/speccomb_mn2,oneminus)
                specmult_mn2 = 8._wp*specparm_mn2
                jmn2 = 1 + int(specmult_mn2)
                fmn2 = MOD1(specmult_mn2)

                speccomb_planck = coln2o(lay)+refrat_planck_a*colco2(lay)
                specparm_planck = MIN(coln2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(15) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(15) + js1
                inds = indself(lay)
                indf = indfor(lay)
                indm = indminor(lay)

                scalen2 = colbrd(lay)*scaleminor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG15
                        tau_major(1:ng15) = speccomb *    &
                                (fac000 * absa(ind0,1:ng15)    + &
                                fac100 * absa(ind0+1,1:ng15)  + &
                                fac200 * absa(ind0+2,1:ng15)  + &
                                fac010 * absa(ind0+9,1:ng15)  + &
                                fac110 * absa(ind0+10,1:ng15) + &
                                fac210 * absa(ind0+11,1:ng15))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG15
                        tau_major(1:ng15) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng15) + &
                                fac100 * absa(ind0,1:ng15)   + &
                                fac000 * absa(ind0+1,1:ng15) + &
                                fac210 * absa(ind0+8,1:ng15) + &
                                fac110 * absa(ind0+9,1:ng15) + &
                                fac010 * absa(ind0+10,1:ng15))
                else
                        !CDIR EXPAND=NG15
                        tau_major(1:ng15) = speccomb *   &
                                (fac000 * absa(ind0,1:ng15)   + &
                                fac100 * absa(ind0+1,1:ng15) + &
                                fac010 * absa(ind0+9,1:ng15) + &
                                fac110 * absa(ind0+10,1:ng15))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG15
                        tau_major1(1:ng15) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng15)    + &
                                fac101 * absa(ind1+1,1:ng15)  + &
                                fac201 * absa(ind1+2,1:ng15)  + &
                                fac011 * absa(ind1+9,1:ng15)  + &
                                fac111 * absa(ind1+10,1:ng15) + &
                                fac211 * absa(ind1+11,1:ng15))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG15
                        tau_major1(1:ng15) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng15) + &
                                fac101 * absa(ind1,1:ng15)   + &
                                fac001 * absa(ind1+1,1:ng15) + &
                                fac211 * absa(ind1+8,1:ng15) + &
                                fac111 * absa(ind1+9,1:ng15) + &
                                fac011 * absa(ind1+10,1:ng15))
                else
                        !CDIR EXPAND=NG15
                        tau_major1(1:ng15) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng15)   + &
                                fac101 * absa(ind1+1,1:ng15) + &
                                fac011 * absa(ind1+9,1:ng15) + &
                                fac111 * absa(ind1+10,1:ng15))
                endif

                !CDIR EXPAND=NG15
                do ig = 1, ng15
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))
                n2m1 = ka_mn2(jmn2,indm,ig) + fmn2 * &
                        (ka_mn2(jmn2+1,indm,ig) - ka_mn2(jmn2,indm,ig))
                n2m2 = ka_mn2(jmn2,indm+1,ig) + fmn2 * &
                        (ka_mn2(jmn2+1,indm+1,ig) - ka_mn2(jmn2,indm+1,ig))
                taun2 = scalen2 * (n2m1 + minorfrac(lay) * (n2m2 - n2m1))

                taug(lay,ngs14+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor &
                        + taun2
                fracs(lay,ngs14+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo

                ! Upper atmosphere part
                do ig = 1, ng15
                !CDIR NODEP,VOVERTAKE,VOB
                taug(lay,ngs14+ig) = 0.0_wp
                fracs(lay,ngs14+ig) = 0.0_wp
                enddo

                enddo

        end subroutine taugb15

        !----------------------------------------------------------------------------
        subroutine taugb16
                !----------------------------------------------------------------------------
                !
                !     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
                !----------------------------------------------------------------------------

                ! ------- Modules -------

                ! >> sylvia_20200320
                ! Add in the USE statements again.
                use rrlw_kg16,   only : fracrefa, fracrefb, absa, absb, &
                       selfref, forref
                !REAL(wp) , DIMENSION(ng16) :: fracrefb

                !REAL(wp) :: fracrefa(ng16,9)
                !REAL(wp) :: absa(585,ng16)
                !REAL(wp) :: absb(235,ng16)
                !REAL(wp) :: selfref(10,ng16)
                !REAL(wp) :: forref(4,ng16)
                ! << sylvia_20200310

                ! >> sylvia_20200310
                INTEGER, PARAMETER :: ng16 = 2
                INTEGER, PARAMETER :: ngs15 = 138
                ! ------- Declarations -------

                ! Local
                ! >> sylvia_20200318, Remove ixc0 and ixp.
                integer :: lay, ig
                integer :: ind0, ind1, inds, indf
                integer :: js, js1, jpl
                real(wp) :: speccomb, specparm, specmult, fs
                real(wp) :: speccomb1, specparm1, specmult1, fs1
                real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
                real(wp) :: p, p4, fk0, fk1, fk2
                real(wp) :: fac000, fac100, fac200
                real(wp) :: fac010, fac110, fac210
                real(wp) :: fac001, fac101, fac201
                real(wp) :: fac011, fac111, fac211
                real(wp) :: tauself, taufor
                real(wp) :: refrat_planck_a
                real(wp) :: tau_major(ng16), tau_major1(ng16)


                ! Calculate reference ratio to be used in calculation of Planck
                ! fraction in lower atmosphere.

                ! P = 387. mb (Level 6)
                refrat_planck_a = chi_mls(1,6)/chi_mls(6,6)

                ! Compute the optical depth by interpolating in ln(pressure),
                ! temperature,and appropriate species.  Below laytrop, the water
                ! vapor self-continuum and foreign continuum is interpolated
                ! (in temperature) separately.

                ! Lower atmosphere loop
                do lay = 1, laytrop_min

                speccomb = colh2o(lay) + rat_h2och4(lay)*colch4(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2och4_1(lay)*colch4(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(16) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(16) + js1
                inds = indself(lay)
                indf = indfor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG16
                        tau_major(1:ng16) = speccomb *    &
                                (fac000 * absa(ind0,1:ng16)    + &
                                fac100 * absa(ind0+1,1:ng16)  + &
                                fac200 * absa(ind0+2,1:ng16)  + &
                                fac010 * absa(ind0+9,1:ng16)  + &
                                fac110 * absa(ind0+10,1:ng16) + &
                                fac210 * absa(ind0+11,1:ng16))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG16
                        tau_major(1:ng16) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng16) + &
                                fac100 * absa(ind0,1:ng16)   + &
                                fac000 * absa(ind0+1,1:ng16) + &
                                fac210 * absa(ind0+8,1:ng16) + &
                                fac110 * absa(ind0+9,1:ng16) + &
                                fac010 * absa(ind0+10,1:ng16))
                else
                        !CDIR EXPAND=NG16
                        tau_major(1:ng16) = speccomb *   &
                                (fac000 * absa(ind0,1:ng16)   + &
                                fac100 * absa(ind0+1,1:ng16) + &
                                fac010 * absa(ind0+9,1:ng16) + &
                                fac110 * absa(ind0+10,1:ng16))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG16
                        tau_major1(1:ng16) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng16)    + &
                                fac101 * absa(ind1+1,1:ng16)  + &
                                fac201 * absa(ind1+2,1:ng16)  + &
                                fac011 * absa(ind1+9,1:ng16)  + &
                                fac111 * absa(ind1+10,1:ng16) + &
                                fac211 * absa(ind1+11,1:ng16))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG16
                        tau_major1(1:ng16) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng16) + &
                                fac101 * absa(ind1,1:ng16)   + &
                                fac001 * absa(ind1+1,1:ng16) + &
                                fac211 * absa(ind1+8,1:ng16) + &
                                fac111 * absa(ind1+9,1:ng16) + &
                                fac011 * absa(ind1+10,1:ng16))
                else
                        !CDIR EXPAND=NG16
                        tau_major1(1:ng16) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng16)   + &
                                fac101 * absa(ind1+1,1:ng16) + &
                                fac011 * absa(ind1+9,1:ng16) + &
                                fac111 * absa(ind1+10,1:ng16))
                endif

                !CDIR EXPAND=NG16
                do ig = 1, ng16
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))

                taug(lay,ngs15+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor
                fracs(lay,ngs15+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo
                enddo

                ! Upper atmosphere loop
                do lay = laytrop_max+1, nlayers
                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(16) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(16) + 1
                !CDIR EXPAND=NG16
                do ig = 1, ng16
                taug(lay,ngs15+ig) = colch4(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig))
                fracs(lay,ngs15+ig) = fracrefb(ig)
                enddo
                enddo

                IF (laytrop_max == laytrop_min) RETURN
                ! Mixed loop
                ! Lower atmosphere part
                do lay = laytrop_min+1, laytrop_max

                !CDIR NODEP,VOVERTAKE,VOB
                speccomb = colh2o(lay) + rat_h2och4(lay)*colch4(lay)
                specparm = MIN(colh2o(lay)/speccomb,oneminus)
                specmult = 8._wp*(specparm)
                js = 1 + int(specmult)
                fs = MOD1(specmult)

                speccomb1 = colh2o(lay) + rat_h2och4_1(lay)*colch4(lay)
                specparm1 = MIN(colh2o(lay)/speccomb1,oneminus)
                specmult1 = 8._wp*(specparm1)
                js1 = 1 + int(specmult1)
                fs1 = MOD1(specmult1)

                speccomb_planck = colh2o(lay)+refrat_planck_a*colch4(lay)
                specparm_planck = MIN(colh2o(lay)/speccomb_planck,oneminus)
                specmult_planck = 8._wp*specparm_planck
                jpl = 1 + int(specmult_planck)
                fpl = MOD1(specmult_planck)

                ind0 = ((jp(lay)-1)*5+(jt(lay)-1))*nspa(16) + js
                ind1 = (jp(lay)*5+(jt1(lay)-1))*nspa(16) + js1
                inds = indself(lay)
                indf = indfor(lay)

                if (specparm .lt. 0.125_wp) then
                        p = fs - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else if (specparm .gt. 0.875_wp) then
                        p = -fs
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac000 = fk0*fac00(lay)
                        fac100 = fk1*fac00(lay)
                        fac200 = fk2*fac00(lay)
                        fac010 = fk0*fac10(lay)
                        fac110 = fk1*fac10(lay)
                        fac210 = fk2*fac10(lay)
                else
                        fac000 = (1._wp - fs) * fac00(lay)
                        fac010 = (1._wp - fs) * fac10(lay)
                        fac100 = fs * fac00(lay)
                        fac110 = fs * fac10(lay)
                        fac200 = 0._wp
                        fac210 = 0._wp
                endif

                if (specparm1 .lt. 0.125_wp) then
                        p = fs1 - 1._wp
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else if (specparm1 .gt. 0.875_wp) then
                        p = -fs1
                        p4 = p**4
                        fk0 = p4
                        fk1 = 1._wp - p - 2.0_wp*p4
                        fk2 = p + p4
                        fac001 = fk0*fac01(lay)
                        fac101 = fk1*fac01(lay)
                        fac201 = fk2*fac01(lay)
                        fac011 = fk0*fac11(lay)
                        fac111 = fk1*fac11(lay)
                        fac211 = fk2*fac11(lay)
                else
                        fac001 = (1._wp - fs1) * fac01(lay)
                        fac011 = (1._wp - fs1) * fac11(lay)
                        fac101 = fs1 * fac01(lay)
                        fac111 = fs1 * fac11(lay)
                        fac201 = 0._wp
                        fac211 = 0._wp
                endif

                if (specparm .lt. 0.125_wp) then
                        !CDIR EXPAND=NG16
                        tau_major(1:ng16) = speccomb *    &
                                (fac000 * absa(ind0,1:ng16)    + &
                                fac100 * absa(ind0+1,1:ng16)  + &
                                fac200 * absa(ind0+2,1:ng16)  + &
                                fac010 * absa(ind0+9,1:ng16)  + &
                                fac110 * absa(ind0+10,1:ng16) + &
                                fac210 * absa(ind0+11,1:ng16))
                else if (specparm .gt. 0.875_wp) then
                        !CDIR EXPAND=NG16
                        tau_major(1:ng16) = speccomb *   &
                                (fac200 * absa(ind0-1,1:ng16) + &
                                fac100 * absa(ind0,1:ng16)   + &
                                fac000 * absa(ind0+1,1:ng16) + &
                                fac210 * absa(ind0+8,1:ng16) + &
                                fac110 * absa(ind0+9,1:ng16) + &
                                fac010 * absa(ind0+10,1:ng16))
                else
                        !CDIR EXPAND=NG16
                        tau_major(1:ng16) = speccomb *   &
                                (fac000 * absa(ind0,1:ng16)   + &
                                fac100 * absa(ind0+1,1:ng16) + &
                                fac010 * absa(ind0+9,1:ng16) + &
                                fac110 * absa(ind0+10,1:ng16))
                endif

                if (specparm1 .lt. 0.125_wp) then
                        !CDIR EXPAND=NG16
                        tau_major1(1:ng16) = speccomb1 *  &
                                (fac001 * absa(ind1,1:ng16)    + &
                                fac101 * absa(ind1+1,1:ng16)  + &
                                fac201 * absa(ind1+2,1:ng16)  + &
                                fac011 * absa(ind1+9,1:ng16)  + &
                                fac111 * absa(ind1+10,1:ng16) + &
                                fac211 * absa(ind1+11,1:ng16))
                else if (specparm1 .gt. 0.875_wp) then
                        !CDIR EXPAND=NG16
                        tau_major1(1:ng16) = speccomb1 * &
                                (fac201 * absa(ind1-1,1:ng16) + &
                                fac101 * absa(ind1,1:ng16)   + &
                                fac001 * absa(ind1+1,1:ng16) + &
                                fac211 * absa(ind1+8,1:ng16) + &
                                fac111 * absa(ind1+9,1:ng16) + &
                                fac011 * absa(ind1+10,1:ng16))
                else
                        !CDIR EXPAND=NG16
                        tau_major1(1:ng16) = speccomb1 * &
                                (fac001 * absa(ind1,1:ng16)   + &
                                fac101 * absa(ind1+1,1:ng16) + &
                                fac011 * absa(ind1+9,1:ng16) + &
                                fac111 * absa(ind1+10,1:ng16))
                endif

                !CDIR EXPAND=NG16
                do ig = 1, ng16
                tauself = selffac(lay)* (selfref(inds,ig) + selffrac(lay) * &
                        (selfref(inds+1,ig) - selfref(inds,ig)))
                taufor =  forfac(lay) * (forref(indf,ig) + forfrac(lay) * &
                        (forref(indf+1,ig) - forref(indf,ig)))

                taug(lay,ngs15+ig) = tau_major(ig) + tau_major1(ig) &
                        + tauself + taufor
                fracs(lay,ngs15+ig) = fracrefa(ig,jpl) + fpl * &
                        (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
                enddo

                ! Upper atmosphere part
                !CDIR NODEP,VOVERTAKE,VOB

                ind0 = ((jp(lay)-13)*5+(jt(lay)-1))*nspb(16) + 1
                ind1 = ((jp(lay)-12)*5+(jt1(lay)-1))*nspb(16) + 1
                !CDIR EXPAND=NG16
                do ig = 1, ng16
                taug(lay,ngs15+ig) = colch4(lay) * &
                        (fac00(lay) * absb(ind0,ig) + &
                        fac10(lay) * absb(ind0+1,ig) + &
                        fac01(lay) * absb(ind1,ig) + &
                        fac11(lay) * absb(ind1+1,ig))
                fracs(lay,ngs15+ig) = fracrefb(ig)
                enddo
                enddo

        end subroutine taugb16

  end subroutine lrtm_taumol
  
end module mo_lrtm_taumol
