!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModFaceValue

  use BATL_lib, ONLY: &
       test_start, test_stop, iTest, jTest, kTest, iBlockTest, iVarTest, &
       iDimTest

  use ModSize, ONLY: nI, nJ, nK, nG, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       x_, y_, z_, nDim, jDim_, kDim_
  use ModVarIndexes
  use ModAdvance, ONLY: UseFDFaceFlux, UseLowOrder, &
       UseLowOrderRegion,IsLowOrderOnly_B, UseAdaptiveLowOrder
  use omp_lib

  implicit none

  private ! except

  public:: read_face_value_param
  public:: calc_face_value
  public:: correct_monotone_restrict
  public:: set_low_order_face
  public:: calc_cell_norm_velocity

  logical, public :: UseAccurateResChange = .false.
  logical, public :: UseTvdResChange      = .true.
  logical, public :: DoLimitMomentum      = .false.
  logical, public :: UseVolumeIntegral4   = .false.
  logical, public :: UseFaceIntegral4     = .false.
  logical, public :: UseLimiter4          = .false.
  integer, public :: nGUsed               = nG

  real,             public :: BetaLimiter = 1.0
  character(len=6), public :: TypeLimiter = 'minmod'
  character(len=6), public :: TypeLimiter5= 'mp'

  logical, public :: UseAccurateExtremum = .true.
  integer:: nLowOrder = 2
  real:: VelCrit, pCritLow, pCritHigh

  ! Logical switch for 5th order scheme: use cweno5 or mp5 scheme.
  logical, public ::  UseCweno = .false.

  logical, public:: UsePerVarLimiter = .false. ! Variable for CWENO5
  integer, public:: iVarSmooth_V(nVar), iVarSmoothIndex_I(nVar)

  ! Region parameters for low order scheme
  character(len=200), public:: StringLowOrderRegion = 'none'
  integer, allocatable, public:: iRegionLowOrder_I(:)

  ! Local variables -----------------

  ! Parameters for the limiter applied near resolution changes
  ! Scheme has no effect if BetaLimierReschange is larger than BetaLimiter
  real    :: BetaLimiterResChange  = 2.0
  integer :: nFaceLimiterResChange = 2

  ! Parameters for limiting the logarithm of variables
  logical :: UseLogLimiter    = .false., UseLogLimiter_V(nVar) = .false.
  logical :: UseLogRhoLimiter = .false.
  logical :: UseLogPLimiter   = .false.

  ! Parameters for limiting the total pressure (p + p_e + p_wave)
  logical :: UsePtotalLtd     = .false.
  logical :: UsePtotalLimiter = .false.

  ! Parameters for limiting the variable divided by density
  logical :: UseScalarToRhoRatioLtd = .false.
  integer :: nVarLimitRatio
  integer, allocatable, save:: iVarLimitRatio_I(:)

  ! Colella's flattening scheme
  logical :: UseFlattening = .true.
  logical :: UseDuFlat     = .false.
  real    :: FlatDelta     = 0.33
  real    :: FlatRatioMin  = 0.75
  real    :: FlatRatioMax  = 0.85

  ! Maximum length of the stencil in 1D
  integer, parameter:: MaxIJK = max(nI,nJ,nK)

  ! index ranges for optimized slope limiter calculations
  integer, parameter:: Lo2_=nVar+1, Hi2_=nVar+nVar, Lo3_=Hi2_+1, Hi3_=nVar+Hi2_

  ! local constants
  real, parameter:: cThird = 1./3., cTwoThird = 2./3., cSixth=1./6.
  real, parameter:: c7over12 = 7.0/12.0, c1over12 = 1.0/12.0

  ! primitive variables
  real, allocatable, save:: Primitive_VG(:,:,:,:)

  ! Variables for "body" blocks with masked cells
  logical:: UseTrueCell
  logical:: IsTrueCell_I(1-nG:MaxIJK+nG)

  ! Low order switch for 1D stencil
  logical:: UseLowOrder_I(1:MaxIJK+1)

  ! variables used for TVD limiters
  real:: dVarLimR_VI(1:nVar,0:MaxIJK+1) ! limited slope for right state
  real:: dVarLimL_VI(1:nVar,0:MaxIJK+1) ! limited slope for left state
  real:: Primitive_VI(1:nVar,1-nG:MaxIJK+nG)

  ! variables for the PPM4 limiter
  integer:: iMin, iMax, jMin, jMax, kMin, kMax
  real:: Cell_I(1-nG:MaxIJK+nG)
  real:: Cell2_I(1-nG:MaxIJK+nG)
  real:: Face_I(0:MaxIJK+2)
  real, allocatable:: FaceL_I(:), FaceR_I(:)
  real:: Prim_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

  real:: LowOrderCrit_I(1:MaxIJK+1)
  !$omp threadprivate( LowOrderCrit_I )

  ! The weight of the four low order polynomials of cweno5
  real, allocatable:: WeightL_II(:,:), WeightR_II(:,:)

  ! OpenMP declaration
!!$omp threadprivate( iVarSmooth_V, iVarSmoothIndex_I ) !questionable?
!!$omp threadprivate( iRegionLowOrder_I ) !this is questionable? 
!!$omp threadprivate( iVarLimitRatio_I )
  !$omp threadprivate( Primitive_VG )
  !$omp threadprivate( UseTrueCell, IsTrueCell_I )
  !$omp threadprivate( UseLowOrder_I )
  !$omp threadprivate( dVarLimR_VI, dVarLimL_VI, Primitive_VI )
  !$omp threadprivate( iMin, iMax, jMin, jMax, kMin, kMax )
  !$omp threadprivate( Cell_I, Cell2_I, Face_I, FaceL_I, FaceR_I, Prim_VG )
  !$omp threadprivate( WeightL_II, WeightR_II)

contains
  !============================================================================
  subroutine read_face_value_param(NameCommand)

    use ModReadParam,  ONLY: read_var, lStringLine
    use ModUtilities,  ONLY: split_string
    use ModVarIndexes, ONLY: NameVar_V

    character(len=*), intent(in) :: NameCommand

    integer :: i, iVar
    character(len=10) :: NameVar_I(nVar)
    character(len=lStringLine) :: NameVarLimitRatio

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_face_value_param'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    select case(NameCommand)
    case('#RESOLUTIONCHANGE')
       call read_var('UseAccurateResChange',  UseAccurateResChange)
       call read_var('UseTvdResChange',       UseTvdResChange)
       call read_var('BetaLimiterResChange',  BetaLimiterResChange)
       call read_var('nFaceLimiterResChange', nFaceLimiterResChange)

       ! Make sure only one of UseAccurateResChange and UseTvdResChange is true
       if(UseAccurateResChange) UseTvdResChange=.false.

    case('#RESCHANGE')
       call read_var('UseAccurateResChange',UseAccurateResChange)
       if(UseAccurateResChange) UseTvdResChange=.false.

    case('#TVDRESCHANGE')
       call read_var('UseTvdResChange', UseTvdResChange)
       if(UseTvdResChange) UseAccurateResChange = .false.

    case("#LIMITER")
       call read_var('UseLogRhoLimiter', UseLogRhoLimiter)
       call read_var('UseLogPLimiter',   UseLogPLimiter)
       call read_var('UseScalarPerRhoLimiter', UseScalarToRhoRatioLtd)
       if(.not. UseScalarToRhoRatioLtd)RETURN
       call read_var('NameVarLimitRatio', NameVarLimitRatio)
       call split_string(NameVarLimitRatio, nVar, NameVar_I, nVarLimitRatio, &
            UseArraySyntaxIn=.true.)
       if(allocated(iVarLimitRatio_I)) deallocate(iVarLimitRatio_I)
       allocate(iVarLimitRatio_I(nVarLimitRatio))
       do i = 1, nVarLimitRatio
          do iVar = 1, nVar
             if(NameVar_V(iVar) == NameVar_I(i))then
                iVarLimitRatio_I(i) = iVar
                EXIT
             end if
          end do
          if(iVar > nVar) call stop_mpi(NameSub// &
               ' could not find NameVarLimitRatio='//NameVar_I(i))
       end do

    case("#LIMITPTOTAL")
       call read_var('UsePtotalLtd', UsePtotalLtd)

    case("#FLATTENING")
       call read_var('UseFlattening', UseFlattening)
       if(UseFlattening)then
          call read_var('UseDuFlat',     UseDuFlat)
          call read_var('FlatDelta',     FlatDelta)
          call read_var('FlatRatioMin',  FlatRatioMin)
          call read_var('FlatRatioMax',  FlatRatioMax)
       end if
    case("#LOWORDERREGION")
       call read_var('StringLowOrderRegion', StringLowOrderRegion)
       UseLowOrderRegion = .true.
    case("#ADAPTIVELOWORDER")
       call read_var('UseAdaptiveLowOrder', UseAdaptiveLowOrder)
       if(UseAdaptiveLowOrder) then
          call read_var('nLowOrder', nLowOrder)
          call read_var('pCritLow',  pCritLow)
          call read_var('pCritHigh', pCritHigh)
          call read_var('VelCrit',   VelCrit)
       endif
       if(.not.(nLowOrder==1 .or. nLowOrder==2)) &
            call stop_mpi(NameSub//' nLowOrder should be 1 or 2!!')
       if(pCritLow < (pCritHigh-1e-15)) &
            call stop_mpi(NameSub//' pCritLow should be >= pCritHigh')
    case default
       call stop_mpi(NameSub//' invalid command='//trim(NameCommand))
    end select

    call test_stop(NameSub, DoTest)
  end subroutine read_face_value_param
  !============================================================================
  subroutine tvd_reschange_body(&
       Coarse2_V         ,& ! State in the coarser ghostcell, 2nd layer
       Coarse1_V         ,& ! State in the coarser ghostcell, 1st layer
       Fine1_VII         ,& ! States in 4 fine physical cells,1st layer
       Fine2_VII         ,& ! States in 4 fine physical cells,2nd layer
       CoarseToFineF_VII ,& ! Values at subfaces, in the coarse ghostcell
       FineToCoarseF_VII ,& ! Facevalues in phys. cell looking at coarser cell
       FineF_VII         ,& ! Facevalues in phys. cell looking at phys. cell
       IsTrueCoarse2     ,& ! True if coarser ghostcell of 2nd layer is true
       IsTrueCoarse1     ,& ! True if coarser ghostcell of 1st layer is true
       IsTrueFine1       ,& ! True if all physical cells of 1st layer are true
       IsTrueFine2_II)      ! True for true physical cell of the 2nd layer
    !_____________!_____________!_______!_______!_
    !             !         CToF! FToC FF!
    ! C2_V        ! C1_V       _!_F1_V__!__F2_V_!_
    !             !         CToF! FToC FF!
    !_____________!_____________!_F1_V__!__F2_V_!_
    !             !             !       !       !
    real,dimension(nVar),intent(in):: Coarse2_V, Coarse1_V
    real,dimension(nVar,2,2),intent(in):: Fine1_VII, Fine2_VII
    real,dimension(nVar,2,2),intent(out):: &
         CoarseToFineF_VII, FineToCoarseF_VII, FineF_VII
    logical,intent(in):: IsTrueCoarse2, IsTrueCoarse1, IsTrueFine1
    logical,dimension(2,2),intent(in):: IsTrueFine2_II
    integer::iVar,i2,j2
    real,dimension(nVar):: AveragedFine1_V,GradNormal_V,SignGradNormal_V
    real,dimension(nVar):: GradNormalLtd_V  ! Ltd stands for "Limited"

    real:: Beta

    character(len=*), parameter:: NameSub = 'tvd_reschange_body'
    !--------------------------------------------------------------------------

    ! Calculate averaged Fine1_VII
    do iVar=1,nVar
       AveragedFine1_V(iVar) = 0.25*sum(Fine1_VII(iVar,:,:))
    end do
    GradNormal_V= AveragedFine1_V-Coarse1_V

    ! Save gradients squared
    SignGradNormal_V=sign(1.0,GradNormal_V)

    Beta = min(BetaLimiterResChange, BetaLimiter)

    if(IsTrueCoarse2.and.IsTrueCoarse1.and.IsTrueFine1)then
       ! Limit gradient in the first coarser cell
       GradNormalLtd_V= SignGradNormal_V*&
            max(0.0,&
            min(Beta*cTwoThird*abs(GradNormal_V),&
            Beta*0.5*SignGradNormal_V*(Coarse1_V-Coarse2_V),&
            cThird*abs(GradNormal_V)+&
            0.25*SignGradNormal_V*(Coarse1_V-Coarse2_V)))

       do j2=1,2;do i2=1,2
          ! Limit transverse gradients, if they are larger than the normal one
          ! The unlimited transverse gradients are Fine1_VII-AveragedFine1_V
          CoarseToFineF_VII(:,i2,j2)=Coarse1_V+GradNormalLtd_V+&
               sign(min(abs(GradNormalLtd_V), &
               abs(Fine1_VII(:,i2,j2)-AveragedFine1_V)),&
               Fine1_VII(:,i2,j2)-AveragedFine1_V)
       end do;end do
    else
       do j2=1,2;do i2=1,2
          ! First order scheme
          CoarseToFineF_VII(:,i2,j2)=Coarse1_V
       end do;end do
    end if
    if(.not.(IsTrueCoarse1.and.IsTrueFine1))then
       do j2=1,2;do i2=1,2
          ! First order scheme
          FineToCoarseF_VII(:,i2,j2)=Fine1_VII(:,i2,j2)
          FineF_VII(:,i2,j2)=Fine1_VII(:,i2,j2)
       end do;end do
    else
       do j2=1,2;do i2=1,2
          if(IsTrueFine2_II(i2,j2))then
             ! Limit gradient in the first layer of finer cells
             GradNormalLtd_V = SignGradNormal_V*&
                  max(0.0,&
                  min(Beta*cThird*abs(GradNormal_V),&
                  SignGradNormal_V*(Fine1_VII(:,i2,j2)-Coarse1_V),&
                  Beta*0.5*SignGradNormal_V*&
                  (Fine2_VII(:,i2,j2)-Fine1_VII(:,i2,j2)),&
                  cSixth*abs(GradNormal_V) + 0.25*SignGradNormal_V*&
                  (Fine2_VII(:,i2,j2)-Fine1_VII(:,i2,j2))))
          else
             ! First order scheme
             GradNormalLtd_V=0.0
          end if
          FineToCoarseF_VII(:,i2,j2)=Fine1_VII(:,i2,j2)-GradNormalLtd_V
          FineF_VII(:,i2,j2)=Fine1_VII(:,i2,j2)+GradNormalLtd_V
       end do;end do
    end if

  end subroutine tvd_reschange_body
  !============================================================================
  subroutine tvd_reschange(&
       Coarse2_V         ,& ! State in the coarser ghostcell, 2nd layer
       Coarse1_V         ,& ! State in the coarser ghostcell, 1st layer
       Fine1_VII         ,& ! States in 4 fine physical cells,1st layer
       Fine2_VII         ,& ! States in 4 fine physical cells,2nd layer
       CoarseToFineF_VII ,& ! Values at subfaces, in the coarse ghostcell
       FineToCoarseF_VII ,& ! Values in phys. cell looking at coarser cell
       FineF_VII)           ! Facevalues in the physical cell,
    !                         looking at another physical cell

    !_____________!_____________!________!_______!_
    !             !         CToF! FToC FF!
    ! C2_V        ! C1_V       _!__F1_V__!__F2_V_!_
    !             !         CToF! FToC FF!
    !_____________!_____________!__F1_V__!__F2_V_!_
    !             !             !        !       !

    real,dimension(nVar),intent(in):: Coarse2_V,Coarse1_V
    real,dimension(nVar,2,2),intent(in):: Fine1_VII,Fine2_VII
    real,dimension(nVar,2,2),intent(inout)::&
         CoarseToFineF_VII ,FineToCoarseF_VII , FineF_VII
    integer::iVar,i2,j2
    real,dimension(nVar):: AveragedFine1_V
    real,dimension(nVar):: GradNormal_V,SignGradNormal_V
    real,dimension(nVar):: GradNormalLtd_V  ! Ltd stands for "Limited"

    real :: Beta
    character(len=*), parameter:: NameSub = 'tvd_reschange'
    !--------------------------------------------------------------------------
    ! Calculate averaged Fine1_VII
    do iVar=1,nVar
       AveragedFine1_V(iVar) = 0.25*sum(Fine1_VII(iVar,:,:))
    end do
    GradNormal_V= AveragedFine1_V-Coarse1_V

    ! Save gradients squared
    SignGradNormal_V=sign(1.0,GradNormal_V)

    Beta = min(BetaLimiterResChange, BetaLimiter)

    ! Limit gradient in the first coarser cell
    GradNormalLtd_V= SignGradNormal_V*max(0.0, min( &
         Beta*cTwoThird*abs(GradNormal_V),&
         Beta*0.5*SignGradNormal_V*(Coarse1_V - Coarse2_V), &
         cThird*abs(GradNormal_V) + 0.25*&
         SignGradNormal_V*(Coarse1_V - Coarse2_V)))

    do j2=1,2;do i2=1,2
       ! Limit transverse gradients, if they are larger than the normal one
       ! Before limiting the transverse gradients are equal to
       ! Fine1_VII-AveragedFine1V
       CoarseToFineF_VII(:,i2,j2) = Coarse1_V + GradNormalLtd_V + &
            sign(min(abs(GradNormalLtd_V),&
            abs(Fine1_VII(:,i2,j2) - AveragedFine1_V)),&
            Fine1_VII(:,i2,j2) - AveragedFine1_V)
    end do;end do

    do j2=1,2;do i2=1,2
       ! Limit gradient in the first layer of finer cells
       GradNormalLtd_V = SignGradNormal_V*max(0.0,min( &
            SignGradNormal_V*(Fine1_VII(:,i2,j2)-Coarse1_V),&
            Beta*cThird*abs(GradNormal_V),&
            Beta*0.5*SignGradNormal_V &
            *(Fine2_VII(:,i2,j2)-Fine1_VII(:,i2,j2)), &
            cSixth*abs(GradNormal_V)+0.25*SignGradNormal_V&
            *(Fine2_VII(:,i2,j2)-Fine1_VII(:,i2,j2)) &
            ))
       FineToCoarseF_VII(:,i2,j2) = Fine1_VII(:,i2,j2)-GradNormalLtd_V
       FineF_VII(:,i2,j2)         = Fine1_VII(:,i2,j2)+GradNormalLtd_V

    end do;end do

  end subroutine tvd_reschange
  !============================================================================
  subroutine accurate_reschange3d(&
       Coarse2_V         ,& ! State in the coarser ghostcell,  2nd layer
       Coarse1_VII       ,& ! State in the coarser ghostcells, 1st layer
       Fine1_VII         ,& ! States in 4 fine physical cells, 1st layer
       Fine2_VII         ,& ! States in 4 fine physical cells, 2nd layer
       CoarseToFineF_VII ,& ! Values at subfaces, in the coarse ghostcell
       FineToCoarseF_VII ,& ! Values in phys. cell looking at coarser cell
       FineF_VII)           ! Values in the phys. cell,
    !                         looking at another physical cell

    !             ! C1_V        !       !       !
    !_____________!_____________!_______!_______!_
    !             !         CToF! FToC FF!       !
    ! C2_V        ! C1_V       _!_F1_V__!__F2_V_!_
    !             !         CToF! FToC FF!       !
    !_____________!_____________!_F1_V__!__F2_V_!_
    !             !             !       !       !
    !             ! C1_V        !       !       !

    real, intent(in) :: Coarse2_V(nVar)
    real, intent(in) :: Coarse1_VII(nVar,-1:4,-1:4)
    real, intent(in) :: Fine1_VII(nVar,2,2)
    real, intent(in) :: Fine2_VII(nVar,2,2)

    real, intent(inout), dimension(nVar, 2, 2)::&
         CoarseToFineF_VII ,FineToCoarseF_VII , FineF_VII

    integer::iVar,i2,j2
    real, dimension(nVar):: AveragedFine1_V, Slope1_V, Slope2_V
    real, dimension(nVar):: GradNormal_V, SignGradNormal_V
    real, dimension(nVar):: GradNormalLtd_V, FaceMiddle_V, Transverse_V
    real, dimension(nVar, 2, 2):: Coarse1Max_VII, Coarse1Min_VII

    real :: AverageOrig_V(nVar), AverageOrig
    real :: Coarse, Middle, FaceAverage, FaceTmp_II(2,2), Alpha, Alpha1
    real :: Beta
    real :: Denominator

    character(len=*), parameter:: NameSub = 'accurate_reschange3d'
    !--------------------------------------------------------------------------

    ! Calculate averaged Fine1_VII
    do iVar=1,nVar
       AveragedFine1_V(iVar) = 0.25*sum(Fine1_VII(iVar,:,:))
    end do
    GradNormal_V = AveragedFine1_V - Coarse1_VII(:,1,1)

    ! Save sign of the gradient
    SignGradNormal_V=sign(1.0,GradNormal_V)

    ! Limit gradient in the first coarser cell
    Slope1_V = cTwoThird*abs(GradNormal_V)
    Slope2_V = 0.5*SignGradNormal_V*(Coarse1_VII(:,1,1) - Coarse2_V)

    Beta = min(BetaLimiterResChange, BetaLimiter)

    GradNormalLtd_V= SignGradNormal_V*max(0.0,min( &
         Beta*Slope1_V, Beta*Slope2_V, 0.5*(Slope1_V+Slope2_V)))

    ! Add limited normal gradient to obtain the middle value for the fine face
    FaceMiddle_V = Coarse1_VII(:,1,1) + GradNormalLtd_V

    do j2=1,2; do i2=1,2
       ! Calculate transverse gradient between coarse cells
       do iVar = 1, nVar
          ! TransverseSlope = ( (Cside1 - Ccenter) + (Cside2 - Ccenter) ) / 4
          Transverse_V(iVar) = 0.0625* &
               ( sum(Coarse1_VII(iVar,4*i2-5:4*i2-4,1:2)) &
               + sum(Coarse1_VII(iVar,1:2,4*j2-5:4*j2-4)) &
               ) - 0.5*Coarse1_VII(iVar,1,1)
       end do

       ! Bound the face value by Coarse1, Coarse1+Transverse and Fine1
       Coarse1Max_VII(:,i2,j2) = Coarse1_VII(:,1,1) + max(0.0, Transverse_V)
       Coarse1Min_VII(:,i2,j2) = Coarse1_VII(:,1,1) + min(0.0, Transverse_V)

       ! Add transverse gradient and limit it
       CoarseToFineF_VII(:,i2,j2) = &
            max( min(Coarse1Min_VII(:,i2,j2), Fine1_VII(:,i2,j2)), &
            min( max(Coarse1Max_VII(:,i2,j2), Fine1_VII(:,i2,j2)), &
            FaceMiddle_V + Transverse_V) )

    end do; end do

    ! The average face value
    AverageOrig_V = 0.25* &
         ( CoarseToFineF_VII(:,1,1) + CoarseToFineF_VII(:,1,2) &
         + CoarseToFineF_VII(:,2,1) + CoarseToFineF_VII(:,2,2) )

    ! For each variable fix the face values if necessary
    do iVar = 1, nVar

       AverageOrig = AverageOrig_V(iVar)
       Coarse      = Coarse1_VII(iVar, 1, 1)
       Middle      = FaceMiddle_V(iVar)

       ! Check if the |L-M| <= |M-C| condition is satisfied
       if(abs(AverageOrig - Middle) <=  abs(Coarse - Middle) ) CYCLE

       ! Calculate the fixed average value L = Lorig + sgn(Lorig-M)*|M-C|
       FaceAverage = Middle + &
            sign( abs(Middle - Coarse), AverageOrig - Middle )

       ! Correct face values either upward or downward
       if(AverageOrig < FaceAverage)then
          FaceTmp_II = &
               max(Coarse1Max_VII(iVar,:,:), CoarseToFineF_VII(iVar,:,:))
       else
          FaceTmp_II = &
               min(Coarse1Min_VII(iVar,:,:), CoarseToFineF_VII(iVar,:,:))
       end if

       ! Calculate interpolation coefficient needed to satisfy the condition
       ! Avoid zero denominator.
       Denominator = 0.25*sum(FaceTmp_II) - AverageOrig
       if(abs(Denominator) < 1e-30) then
          CoarseToFineF_VII(iVar,:,:) = FaceTmp_II
       else
          Alpha = (FaceAverage - AverageOrig) / Denominator
          Alpha1 = 1.0 - Alpha

          ! Interpolate
          CoarseToFineF_VII(iVar,:,:) = &
               Alpha*FaceTmp_II + Alpha1*CoarseToFineF_VII(iVar,:,:)
       endif
    end do

    ! The face is half the distance in the fine cell
    Slope1_V = 0.5*Slope1_V
    do j2=1,2; do i2=1,2
       ! Limit gradient in the first layer of finer cells
       Slope2_V = 0.5*SignGradNormal_V*(Fine2_VII(:,i2,j2)-Fine1_VII(:,i2,j2))

       ! The first limiting ensures that the FineToCoarse face value
       ! remains between the Fine1 and Coarse values
       GradNormalLtd_V = SignGradNormal_V*max(0.0,min( &
            SignGradNormal_V*(Fine1_VII(:,i2,j2) - Coarse1_VII(:,1,1)), &
            Beta*Slope1_V, Beta*Slope2_V, 0.5*(Slope1_V + Slope2_V)))

       FineToCoarseF_VII(:,i2,j2) = Fine1_VII(:,i2,j2) - GradNormalLtd_V
       FineF_VII(:,i2,j2) = Fine1_VII(:,i2,j2) + GradNormalLtd_V

    end do; end do

  end subroutine accurate_reschange3d
  !============================================================================
  subroutine accurate_reschange2d(&
       Coarse2_V         ,& ! State in the coarser ghostcell,  2nd layer
       Coarse1_VI        ,& ! State in the coarser ghostcells, 1st layer
       Fine1_VI          ,& ! States in 2 fine physical cells, 1st layer
       Fine2_VI          ,& ! States in 2 fine physical cells, 2nd layer
       CoarseToFineF_VI  ,& ! Values at subfaces, in the coarse ghostcell
       FineToCoarseF_VI  ,& ! Values in phys. cell looking at coarser cell
       FineF_VI)            ! Facevalues in the physical cell,
    !                         looking at another physical cell

    !             ! C1_V        !       !       !
    !_____________!_____________!_______!_______!_
    !             !         CToF! FToC FF!       !
    ! C2_V        ! C1_V       _!_F1_V__!__F2_V_!_
    !             !         CToF! FToC FF!       !
    !_____________!_____________!_F1_V__!__F2_V_!_
    !             !             !       !       !
    !             ! C1_V        !       !       !

    real, intent(in) :: Coarse2_V(nVar)
    real, intent(in) :: Coarse1_VI(nVar,-1:4)
    real, intent(in) :: Fine1_VI(nVar,2)
    real, intent(in) :: Fine2_VI(nVar,2)

    real, intent(inout), dimension(nVar, 2)::&
         CoarseToFineF_VI ,FineToCoarseF_VI , FineF_VI

    integer:: iVar, i2
    real, dimension(nVar):: AveragedFine1_V, Slope1_V, Slope2_V
    real, dimension(nVar):: GradNormal_V, SignGradNormal_V
    real, dimension(nVar):: GradNormalLtd_V, FaceMiddle_V, Transverse_V
    real, dimension(nVar,2):: Coarse1Max_VI, Coarse1Min_VI

    real :: AverageOrig_V(nVar), AverageOrig
    real :: Coarse, Middle, FaceAverage, FaceTmp_I(2), Alpha, Alpha1
    real :: Beta
    real :: Denominator

    character(len=*), parameter:: NameSub = 'accurate_reschange2d'
    !--------------------------------------------------------------------------

    ! Calculate averaged Fine1_VI
    do iVar=1,nVar
       AveragedFine1_V(iVar) = 0.5*sum(Fine1_VI(iVar,:))
    end do
    GradNormal_V = AveragedFine1_V - Coarse1_VI(:,1)

    ! Save sign of the gradient
    SignGradNormal_V=sign(1.0,GradNormal_V)

    ! Limit gradient in the first coarser cell
    Slope1_V = cTwoThird*abs(GradNormal_V)
    Slope2_V = 0.5*SignGradNormal_V*(Coarse1_VI(:,1) - Coarse2_V)

    Beta = min(BetaLimiterResChange, BetaLimiter)

    GradNormalLtd_V= SignGradNormal_V*max(0.0,min( &
         Beta*Slope1_V, Beta*Slope2_V, 0.5*(Slope1_V+Slope2_V)))

    ! Add limited normal gradient to obtain the middle value for the fine face
    FaceMiddle_V = Coarse1_VI(:,1) + GradNormalLtd_V

    do i2 = 1, 2
       ! Calculate transverse gradient between coarse cells
       do iVar = 1, nVar
          ! TransverseSlope = (Cside1 - Ccenter) / 4
          Transverse_V(iVar) = &
               0.125*sum(Coarse1_VI(iVar,4*i2-5:4*i2-4)) &
               - 0.25*Coarse1_VI(iVar,1)
       end do

       ! Bound the face value by Coarse1, Coarse1+Transverse and Fine1
       Coarse1Max_VI(:,i2) = Coarse1_VI(:,1) + max(0.0, Transverse_V)
       Coarse1Min_VI(:,i2) = Coarse1_VI(:,1) + min(0.0, Transverse_V)

       ! Add transverse gradient and limit it
       CoarseToFineF_VI(:,i2) = &
            max( min(Coarse1Min_VI(:,i2), Fine1_VI(:,i2)), &
            min( max(Coarse1Max_VI(:,i2), Fine1_VI(:,i2)), &
            FaceMiddle_V + Transverse_V) )

    end do

    ! The average face value
    AverageOrig_V = 0.5*( CoarseToFineF_VI(:,1) + CoarseToFineF_VI(:,2) )

    ! For each variable fix the face values if necessary
    do iVar = 1, nVar

       AverageOrig = AverageOrig_V(iVar)
       Coarse      = Coarse1_VI(iVar,1)
       Middle      = FaceMiddle_V(iVar)

       ! Check if the |L-M| <= |M-C| condition is satisfied
       if(abs(AverageOrig - Middle) <=  abs(Coarse - Middle) ) CYCLE

       ! Calculate the fixed average value L = Lorig + sgn(Lorig-M)*|M-C|
       FaceAverage = Middle + &
            sign( abs(Middle - Coarse), AverageOrig - Middle )

       ! Correct face values either upward or downward
       if(AverageOrig < FaceAverage)then
          FaceTmp_I = max(Coarse1Max_VI(iVar,:), CoarseToFineF_VI(iVar,:))
       else
          FaceTmp_I = min(Coarse1Min_VI(iVar,:), CoarseToFineF_VI(iVar,:))
       end if

       ! Calculate interpolation coefficient needed to satisfy the condition
       ! Avoid zero denominator.
       Denominator = 0.5*sum(FaceTmp_I) - AverageOrig
       if(abs(Denominator) < 1e-30) then
          CoarseToFineF_VI(iVar,:) = FaceTmp_I
       else
          Alpha = (FaceAverage - AverageOrig) / Denominator
          Alpha1 = 1.0 - Alpha
          ! Interpolate
          CoarseToFineF_VI(iVar,:) = &
               Alpha*FaceTmp_I + Alpha1*CoarseToFineF_VI(iVar,:)
       endif

    end do

    ! The face is half the distance in the fine cell
    Slope1_V = 0.5*Slope1_V
    do i2 = 1, 2
       ! Limit gradient in the first layer of finer cells
       Slope2_V = 0.5*SignGradNormal_V*(Fine2_VI(:,i2) - Fine1_VI(:,i2))

       ! The first limiting ensures that the FineToCoarse face value
       ! remains between the Fine1 and Coarse values
       GradNormalLtd_V = SignGradNormal_V*max(0.0,min( &
            SignGradNormal_V*(Fine1_VI(:,i2) - Coarse1_VI(:,1)), &
            Beta*Slope1_V, Beta*Slope2_V, 0.5*(Slope1_V + Slope2_V)))

       FineToCoarseF_VI(:,i2) = Fine1_VI(:,i2) - GradNormalLtd_V
       FineF_VI(:,i2) = Fine1_VI(:,i2) + GradNormalLtd_V

    end do

  end subroutine accurate_reschange2d
  !============================================================================
  subroutine accurate_reschange1d(&
       Coarse2_V         ,& ! State in the coarser ghostcell,  2nd layer
       Coarse1_V         ,& ! State in the coarser ghostcells, 1st layer
       Fine1_V           ,& ! States in 2 fine physical cells, 1st layer
       Fine2_V           ,& ! States in 2 fine physical cells, 2nd layer
       CoarseToFineF_V   ,& ! Values at face, in the coarse ghostcell
       FineToCoarseF_V   ,& ! Values in the phys. cell looking at coarser cell
       FineF_V)             ! Facevalues in the physical cell,
    !                         looking at another physical cell

    !_____________!_____________!_______!_______!_
    !             !         CToF! FToC FF!       !
    ! C2_V        ! C1_V        ! F1_V  !  F2_V !
    !_____________!_____________!_______!_______!_

    real, intent(in) :: Coarse2_V(nVar)
    real, intent(in) :: Coarse1_V(nVar)
    real, intent(in) :: Fine1_V(nVar)
    real, intent(in) :: Fine2_V(nVar)

    real, intent(inout), dimension(nVar)::&
         CoarseToFineF_V ,FineToCoarseF_V , FineF_V

    real, dimension(nVar):: Slope1_V, Slope2_V
    real, dimension(nVar):: GradNormal_V, SignGradNormal_V, GradNormalLtd_V

    real :: Beta

    character(len=*), parameter:: NameSub = 'accurate_reschange1d'
    !--------------------------------------------------------------------------

    ! Calculate averaged Fine1_VI
    GradNormal_V = Fine1_V - Coarse1_V

    ! Save sign of the gradient
    SignGradNormal_V=sign(1.0,GradNormal_V)

    ! Limit gradient in the first coarser cell
    Slope1_V = cTwoThird*abs(GradNormal_V)
    Slope2_V = 0.5*SignGradNormal_V*(Coarse1_V - Coarse2_V)

    Beta = min(BetaLimiterResChange, BetaLimiter)

    GradNormalLtd_V = SignGradNormal_V*max(0.0,min( &
         Beta*Slope1_V, Beta*Slope2_V, 0.5*(Slope1_V+Slope2_V)))

    ! Add limited normal gradient to obtain the middle value for the fine face
    CoarseToFineF_V = Coarse1_V + GradNormalLtd_V

    ! The face is half the distance in the fine cell
    Slope1_V = 0.5*Slope1_V
    ! Limit gradient in the first layer of finer cells
    Slope2_V = 0.5*SignGradNormal_V*(Fine2_V - Fine1_V)

    ! The first limiting ensures that the FineToCoarse face value
    ! remains between the Fine1 and Coarse values
    GradNormalLtd_V = SignGradNormal_V*max(0.0,min( &
         SignGradNormal_V*(Fine1_V - Coarse1_V), &
         Beta*Slope1_V, Beta*Slope2_V, 0.5*(Slope1_V + Slope2_V)))

    FineToCoarseF_V = Fine1_V - GradNormalLtd_V
    FineF_V         = Fine1_V + GradNormalLtd_V

  end subroutine accurate_reschange1d
  !============================================================================
  subroutine calc_face_value(DoResChangeOnly, iBlock)

    use ModMultiFluid, ONLY: nIonFluid, iRho, iUx, iUz, iUx_I, iUz_I

    ! The subroutine calculates right and left face values (primitive
    ! variables) LeftState_VX .. RightState_VZ for block iBlock from
    ! the cell centered State_VGB.
    !
    ! If DoResChangeOnly is true, only facevalues next to a coarser
    ! neighbor block are calculated.

    use ModMain,     ONLY: nOrder, nOrderProlong, UseB0, &
         UseConstrainB, nIFace, nJFace, nKFace, &
         iMinFace, iMaxFace, jMinFace, jMaxFace, kMinFace, kMaxFace, &
         iMinFace2, iMaxFace2, jMinFace2, jMaxFace2, kMinFace2, kMaxFace2, &
         UseHighResChange

    use ModGeometry, ONLY : true_cell, body_BLK
    use ModPhysics, ONLY: GammaWave, C2light, InvClight2, Gamma_I
    use ModB0
    use ModAdvance, ONLY: State_VGB, Energy_GBI, &
         DoInterpolateFlux, FluxLeft_VGD, FluxRight_VGD, &
         Flux_VX, Flux_VY, Flux_VZ, &
         uDotArea_XI, uDotArea_YI, uDotArea_ZI, &
         UseElectronPressure, UseWavePressure, UseAnisoPressure, UseAnisoPe, &
         LeftState_VX,      &  ! Face Left  X
         RightState_VX,     &  ! Face Right X
         LeftState_VY,      &  ! Face Left  Y
         RightState_VY,     &  ! Face Right Y
         LeftState_VZ,      &  ! Face Left  Z
         RightState_VZ,     &  ! Face Right Z
         LowOrderCrit_XB, LowOrderCrit_YB, LowOrderCrit_ZB


    use ModParallel, ONLY : &
         neiLEV,neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth

    use ModEnergy, ONLY: calc_pressure

    use ModViscosity, ONLY: UseArtificialVisco

    use BATL_lib, ONLY: CellFace_DB

    logical, intent(in):: DoResChangeOnly
    integer, intent(in):: iBlock

    integer:: i, j, k, iSide, iFluid, iFlux
    real:: RhoInv

    real:: RhoC2Inv, BxFull, ByFull, BzFull, B2Full, uBC2Inv, Ga2Boris
    real:: B0_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    ! Number of cells needed to get the face values
    integer:: nStencil

    ! Variables related to 4th order finite volume scheme
    real, parameter:: c24th = 1.0/24.0
    real:: Laplace_V(nVar), Laplace

    real:: State_V(nVar), Energy
    integer:: iVarSmoothLast, iVarSmooth

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_face_value'
    !--------------------------------------------------------------------------
    if(.not. DoResChangeOnly)then
       call test_start(NameSub, DoTest, iBlock)
    else
       DoTest=.false.
    end if

    if(DoTest)then
       write(*,*) NameSub,' starting with DoResChangeOnly=', DoResChangeOnly
       if(iDimTest==0 .or. iDimTest==1)then
          write(*,*)'TestVar(iTest-nG:iTest+nG)=', &
               State_VGB(iVarTest,iTest-nG:iTest+nG,jTest,kTest,iBlockTest)
          if(.not.all(true_cell(iTest-nG:iTest+nG,jTest,kTest,iBlockTest))) &
               write(*,*)'true_cell(iTest-nG:iTest+nG)=',&
               true_cell(iTest-nG:iTest+nG,jTest,kTest,iBlockTest)
       end if
       if(nDim > 1 .and. (iDimTest==0 .or. iDimTest==2))then
          write(*,*)'TestVar(jTest-nG:jTest+nG)=', &
               State_VGB(iVarTest,iTest,jTest-nG:jTest+nG,kTest,iBlockTest)
          if(.not.all(true_cell(iTest,jTest-nG:jTest+nG,kTest,iBlockTest))) &
               write(*,*)'true_cell(jTest-nG:jTest+nG)=',&
               true_cell(iTest,jTest-nG:jTest+nG,kTest,iBlockTest)
       end if
       if(nDim > 2 .and. (iDimTest==0 .or. iDimTest==3)) then
          write(*,*)'TestVar(kTest-nG:kTest+nG)=', &
               State_VGB(iVarTest,iTest,jTest,kTest-nG:kTest+nG,iBlockTest)
          if(.not.all(true_cell(iTest,jTest,kTest-nG:kTest+nG,iBlockTest))) &
               write(*,*)'true_cell(kTest-nG:kTest+nG)=', &
               true_cell(iTest,jTest,kTest-nG:kTest+nG,iBlockTest)
       end if
    end if

    if(.not.allocated(Primitive_VG))&
         allocate(Primitive_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

    if(.not. allocated(FaceL_I)) then
       allocate(FaceL_I(1:MaxIJK+2))
       allocate(FaceR_I(0:MaxIJK+1))
       allocate(WeightL_II(-2:2,0:MaxIJK+1))
       allocate(WeightR_II(-2:2,0:MaxIJK+1))
    endif

    UseTrueCell = body_BLK(iBlock)

    UseLogLimiter   = nOrder > 1 .and. (UseLogRhoLimiter .or. UseLogPLimiter)
    UseLogLimiter_V = .false.
    if(UseLogLimiter)then
       if(UseLogRhoLimiter)then
          do iFluid = 1, nFluid
             UseLogLimiter_V(iRho_I(iFluid)) = .true.
          end do
       end if
       if(UseLogPLimiter)then
          do iFluid = 1, nFluid
             UseLogLimiter_V(iP_I(iFluid))   = .true.
          end do
          if(UseAnisoPressure)then
             do iFluid = IonFirst_, IonLast_
                UseLogLimiter_V(iPparIon_I(iFluid)) = .true.
             end do
          end if
          if(UseElectronPressure) UseLogLimiter_V(Pe_) = .true.
          if(UseAnisoPe)          UseLogLimiter_V(Pepar_) = .true.
       end if
    end if

    UsePtotalLimiter = nOrder > 1 .and. nIonFluid == 1 .and. UsePtotalLtd

    if(.not.DoResChangeOnly & ! In order not to call it twice
         .and. nOrder > 1   & ! Is not needed for nOrder=1
         .and. (UseAccurateResChange .or. UseTvdResChange)) &
         call correct_monotone_restrict(iBlock)

    ! first, calculate the CELL values for the variables to be limited
    ! for non-boris corrections they are: density, velocity, pressure
    ! for boris correction momentum is used instead of the velocity

    ! Number of cells away from the cell center
    if(nOrder == 5)then
       nStencil = 3
    elseif(nOrder == 4)then
       nStencil = nG
    else
       nStencil = nOrder
    end if

    if(DoLimitMomentum)then
       if(UseB0)then
          B0_DG=B0_DGB(:,:,:,:,iBlock)
       else
          B0_DG=0.00
       end if
       do k = kMinFace, kMaxFace
          do j = jMinFace, jMaxFace
             do i=1-nStencil,nI+nStencil
                call calc_primitives_boris    ! needed for x-faces
             end do
          end do
       end do
       if(nJ > 1)then
          do k = kMinFace, kMaxFace; do i = iMinFace, iMaxFace
             do j=1-nStencil,jMinFace-1
                call calc_primitives_boris    ! additional calculations for
             end do                           ! y -faces
             do j=jMaxFace+1,nJ+nStencil
                call calc_primitives_boris    ! additional calculations for
             end do	                      ! y-faces
          end do; end do
       end if
       if(nK > 1)then
          do j=jMinFace,jMaxFace; do i=iMinFace,iMaxFace
             do k=1-nStencil,kMinFace-1
                call calc_primitives_boris    ! additional calculations for
             end do                           ! z-faces
             do k=kMaxFace+1,nK+nStencil
                call calc_primitives_boris    ! additional calculations for
             end do	                      ! z-faces
          end do; end do
       end if
    else
       if(UseAccurateResChange .or. nOrder==4)then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             call calc_primitives_MHD         ! all cells
          end do; end do; end do
          if(nOrder == 4 .and. UseVolumeIntegral4)then
             ! Calculate 4th order accurate cell averaged primitive variables

             ! First get 4th order accurate cell centered conservative vars
             iMin = MinI + 1; iMax = MaxI - 1
             jMin = MinJ + jDim_; jMax = MaxJ - jDim_
             kMin = MinK + kDim_; kMax = MaxK - kDim_

             ! Store primitive and conservative values based on cell averages
             ! These are used to do the Laplace operators for corrections
             Prim_VG = Primitive_VG

             ! Convert to pointwise conservative variable (eq. 12)
             do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax

                ! Store cell averaged value
                State_V = State_VGB(:,i,j,k,iBlock)

                ! Calculate 4th order accurate cell center value (eq 12)
                Laplace_V = -2*nDim*State_V + &
                     State_VGB(:,i-1,j,k,iBlock) + State_VGB(:,i+1,j,k,iBlock)
                if(nJ > 1) Laplace_V = Laplace_V + &
                     State_VGB(:,i,j-1,k,iBlock) + State_VGB(:,i,j+1,k,iBlock)
                if(nK > 1) Laplace_V = Laplace_V + &
                     State_VGB(:,i,j,k-1,iBlock) + State_VGB(:,i,j,k+1,iBlock)
                State_VGB(:,i,j,k,iBlock) = State_V - c24th*Laplace_V

                do iFluid = 1, nFluid
                   ! Store cell averaged energy
                   Energy = Energy_GBI(i,j,k,iBlock,iFluid)

                   ! Calculate 4th order accurate cell center energy
                   Laplace = -2*nDim*Energy &
                        + Energy_GBI(i-1,j,k,iBlock,iFluid) &
                        + Energy_GBI(i+1,j,k,iBlock,iFluid)
                   if(nJ > 1) Laplace = Laplace &
                        + Energy_GBI(i,j-1,k,iBlock,iFluid) &
                        + Energy_GBI(i,j+1,k,iBlock,iFluid)
                   if(nK > 1) Laplace = Laplace &
                        + Energy_GBI(i,j,k-1,iBlock,iFluid) &
                        + Energy_GBI(i,j,k+1,iBlock,iFluid)
                   Energy_GBI(i,j,k,iBlock,iFluid) = Energy - c24th*Laplace
                   ! check positivity !!!

                   ! Get 4th order accurate cell center pressure
                   call calc_pressure(i,i,j,j,k,k,iBlock,iFluid,iFluid)

                   ! Restore cell averaged energy
                   Energy_GBI(i,j,k,iBlock,iFluid) = Energy
                end do

!!! check positivity ???
                ! Convert to pointwise primitive variables
                call calc_primitives_MHD

                ! Convert to cell averaged primitive variables (eq. 16)
                Laplace_V = Prim_VG(:,i-1,j,k) + Prim_VG(:,i+1,j,k) &
		     - 2*nDim*Prim_VG(:,i,j,k)
                if(nJ > 1) Laplace_V = Laplace_V &
                     + Prim_VG(:,i,j-1,k) + Prim_VG(:,i,j+1,k)
                if(nK > 1) Laplace_V = Laplace_V &
                     + Prim_VG(:,i,j,k-1) + Prim_VG(:,i,j,k+1)

                Primitive_VG(:,i,j,k) = Primitive_VG(:,i,j,k) + c24th*Laplace_V

                ! Restore cell averaged state
                State_VGB(:,i,j,k,iBlock) = State_V
             end do; end do; end do
          end if
       else
          do k=kMinFace,kMaxFace
             do j=jMinFace,jMaxFace
                do i=1-nStencil,nI+nStencil
                   call calc_primitives_MHD   ! needed for x-faces
                end do
             end do
          end do
          if(nJ > 1)then
             do k=kMinFace,kMaxFace; do i=iMinFace,iMaxFace
                do j=1-nStencil,jMinFace-1
                   call calc_primitives_MHD   ! additional calculations for
                end do                        ! y -faces
                do j=jMaxFace+1,nJ+nStencil
                   call calc_primitives_MHD   ! additional calculations for
                end do	                      ! y-faces
             end do; end do
          end if
          if(nK > 1)then
             do j=jMinFace,jMaxFace; do i=iMinFace,iMaxFace
                do k=1-nStencil,kMinFace-1
                   call calc_primitives_MHD   ! additional calculations for
                end do                        ! z-faces
                do k=kMaxFace+1,nK+nStencil
                   call calc_primitives_MHD   ! additional calculations for
                end do	                      ! z-faces
             end do; end do
          end if
       end if
    end if

    if(UseArtificialVisco) call calc_face_div_u(iBlock)



    !\
    ! Now the first or second order face values are calcuted
    !/
    select case(nOrder)
    case(1)
       ! First order reconstruction
       if (.not.DoResChangeOnly) then
          call get_faceX_first(&
               1,nIFace,jMinFace,jMaxFace,kMinFace,kMaxFace)
          if(nJ > 1) call get_faceY_first(&
               iMinFace,iMaxFace,1,nJFace,kMinFace,kMaxFace)
          if(nK > 1) call get_faceZ_first(&
               iMinFace,iMaxFace,jMinFace,jMaxFace,1,nKFace)
       else
          if(neiLeast(iBlock)==+1)&
               call get_faceX_first(1,1,1,nJ,1,nK)
          if(neiLwest(iBlock)==+1)&
               call get_faceX_first(nIFace,nIFace,1,nJ,1,nK)
          if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
               call get_faceY_first(1,nI,1,1,1,nK)
          if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
               call get_faceY_first(1,nI,nJFace,nJFace,1,nK)
          if(nK > 1 .and. neiLbot(iBlock)==+1) &
               call get_faceZ_first(1,nI,1,nJ,1,1)
          if(nK > 1 .and. neiLtop(iBlock)==+1) &
               call get_faceZ_first(1,nI,1,nJ,nKFace,nKFace)
       end if
    case default

       if (.not.DoResChangeOnly)then
          ! Calculate all face values with high order scheme
          if(nOrder==2 .or. IsLowOrderOnly_B(iBlock))then
             ! Second order scheme
             call get_faceX_second(&
                  1,nIFace,jMinFace,jMaxFace,kMinFace,kMaxFace)
             if(nJ > 1) call get_faceY_second(&
                  iMinFace,iMaxFace,1,nJFace,kMinFace,kMaxFace)
             if(nK > 1) call get_faceZ_second(&
                  iMinFace,iMaxFace,jMinFace,jMaxFace,1,nKFace)
          else
             ! High order scheme
             call get_facex_high(&
                  1,nIFace,jMinFace2,jMaxFace2,kMinFace2,kMaxFace2)
             if(nJ > 1) call get_facey_high(&
                  iMinFace2,iMaxFace2,1,nJFace,kMinFace2,kMaxFace2)
             if(nK > 1) call get_facez_high(&
                  iMinFace2,iMaxFace2,jMinFace2,jMaxFace2,1,nKFace)
          end if
       end if

       ! Now take care of faces at resolution changes
       if(nOrderProlong==1 .and..not.UseConstrainB &
            .and. .not.UseHighResChange)then

          ! If nOrderProlong is 1 then use TVD reschange or accurate reschange
          ! scheme and overwrite face values at resolution changes

          if(nJ == 1 .and. (UseAccurateResChange .or. UseTvdResChange))then
             do iSide = 1, 2
                if(neiLev(iSide,iBlock) == 1)call get_face_accurate1d(iSide)
             end do
          elseif (UseAccurateResChange)then
             if(nK == 1)then
                do iSide = 1, 4
                   if(neilev(iSide,iBlock) == 1)call get_face_accurate2d(iSide)
                end do
             else
                do iSide = 1, 6
                   if(neilev(iSide,iBlock) == 1)call get_face_accurate3d(iSide)
                end do
             end if
          else if(UseTvdResChange)then
             do iSide=1,6
                if(neilev(iSide,iBlock) == +1)call get_face_tvd(iSide)
             end do
          else
             ! First order facevalues at resolution change
             if(neiLeast(iBlock)==+1)&
                  call get_faceX_first(1,1,1,nJ,1,nK)
             if(neiLwest(iBlock)==+1)&
                  call get_faceX_first(nIFace,nIFace,1,nJ,1,nK)
             if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
                  call get_faceY_first(1,nI,1,1,1,nK)
             if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
                  call get_faceY_first(1,nI,nJFace,nJFace,1,nK)
             if(nK > 1 .and. neiLbot(iBlock)==+1) &
                  call get_faceZ_first(1,nI,1,nJ,1,1)
             if(nK > 1 .and. neiLtop(iBlock)==+1) &
                  call get_faceZ_first(1,nI,1,nJ,nKFace,nKFace)
          end if

       else if(DoResChangeOnly) then
          if(nOrder==2 .or. IsLowOrderOnly_B(iBlock))then
             ! Second order face values at resolution changes
             if(neiLeast(iBlock)==+1)&
                  call get_faceX_second(1,1,1,nJ,1,nK)
             if(neiLwest(iBlock)==+1)&
                  call get_faceX_second(nIFace,nIFace,1,nJ,1,nK)
             if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
                  call get_faceY_second(1,nI,1,1,1,nK)
             if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
                  call get_faceY_second(1,nI,nJFace,nJFace,1,nK)
             if(nK > 1 .and. neiLbot(iBlock)==+1) &
                  call get_faceZ_second(1,nI,1,nJ,1,1)
             if(nK > 1 .and. neiLtop(iBlock)==+1) &
                  call get_faceZ_second(1,nI,1,nJ,nKFace,nKFace)
          else
             ! High order face values at resolution changes
             if(neiLeast(iBlock)==+1)&
                  call get_faceX_high(1,1,1,nJ,1,nK)
             if(neiLwest(iBlock)==+1)&
                  call get_faceX_high(nIFace,nIFace,1,nJ,1,nK)
             if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
                  call get_faceY_high(1,nI,1,1,1,nK)
             if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
                  call get_faceY_high(1,nI,nJFace,nJFace,1,nK)
             if(nK > 1 .and. neiLbot(iBlock)==+1) &
                  call get_faceZ_high(1,nI,1,nJ,1,1)
             if(nK > 1 .and. neiLtop(iBlock)==+1) &
                  call get_faceZ_high(1,nI,1,nJ,nKFace,nKFace)
          end if
       endif

       if(UseLogLimiter .and. .not.DoLimitMomentum)then
          if(DoResChangeOnly)then
             if(neiLeast(iBlock)==+1) &
                  call logfaceX_to_faceX(1,1,1,nJ,1,nK)
             if(neiLwest(iBlock)==+1) &
                  call logfaceX_to_faceX(nIFace,nIFace,1,nJ,1,nK)
             if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
                  call logfaceY_to_faceY(1,nI,1,1,1,nK)
             if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
                  call logfaceY_to_faceY(1,nI,nJFace,nJFace,1,nK)
             if(nK > 1 .and. neiLbot(iBlock)==+1) &
                  call logfaceZ_to_faceZ(1,nI,1,nJ,1,1)
             if(nK > 1 .and. neiLtop(iBlock)==+1) &
                  call logfaceZ_to_faceZ(1,nI,1,nJ,nKFace,nKFace)
          else
             call logfaceX_to_faceX(1,nIFace,jMinFace,jMaxFace, &
                  kMinFace,kMaxFace)
             if(nJ > 1) call logfaceY_to_faceY(iMinFace,iMaxFace,1,nJFace, &
                  kMinFace,kMaxFace)
             if(nK > 1) call logfaceZ_to_faceZ(iMinFace,iMaxFace, &
                  jMinFace,jMaxFace,1,nKFace)
          end if
       end if

       if(UseScalarToRhoRatioLtd)then
          if(DoResChangeOnly)then
             if(neiLeast(iBlock)==+1) &
                  call ratio_to_scalar_faceX(1,1,1,nJ,1,nK)
             if(neiLwest(iBlock)==+1) &
                  call ratio_to_scalar_faceX(nIFace,nIFace,1,nJ,1,nK)
             if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
                  call ratio_to_scalar_faceY(1,nI,1,1,1,nK)
             if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
                  call ratio_to_scalar_faceY(1,nI,nJFace,nJFace,1,nK)
             if(nK > 1 .and. neiLbot(iBlock)==+1) &
                  call ratio_to_scalar_faceZ(1,nI,1,nJ,1,1)
             if(nK > 1 .and. neiLtop(iBlock)==+1) &
                  call ratio_to_scalar_faceZ(1,nI,1,nJ,nKFace,nKFace)
          else
             call ratio_to_scalar_faceX(1,nIFace,jMinFace,jMaxFace, &
                  kMinFace,kMaxFace)
             if(nJ > 1) call ratio_to_scalar_faceY(iMinFace,iMaxFace,1,nJFace,&
                  kMinFace,kMaxFace)
             if(nK > 1) call ratio_to_scalar_faceZ(iMinFace,iMaxFace, &
                  jMinFace,jMaxFace,1,nKFace)
          end if
       end if

       if(UsePtotalLimiter)then
          if(DoResChangeOnly)then
             if(neiLeast(iBlock)==+1) &
                  call ptotal_to_p_faceX(1,1,1,nJ,1,nK)
             if(neiLwest(iBlock)==+1) &
                  call ptotal_to_p_faceX(nIFace,nIFace,1,nJ,1,nK)
             if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
                  call ptotal_to_p_faceY(1,nI,1,1,1,nK)
             if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
                  call ptotal_to_p_faceY(1,nI,nJFace,nJFace,1,nK)
             if(nK > 1 .and. neiLbot(iBlock)==+1) &
                  call ptotal_to_p_faceZ(1,nI,1,nJ,1,1)
             if(nK > 1 .and. neiLtop(iBlock)==+1) &
                  call ptotal_to_p_faceZ(1,nI,1,nJ,nKFace,nKFace)
          else
             call ptotal_to_p_faceX(1,nIFace,jMinFace,jMaxFace, &
                  kMinFace,kMaxFace)
             if(nJ > 1) call ptotal_to_p_faceY(iMinFace,iMaxFace,1,nJFace, &
                  kMinFace,kMaxFace)
             if(nK > 1) call ptotal_to_p_faceZ(iMinFace,iMaxFace, &
                  jMinFace,jMaxFace,1,nKFace)
          end if
       end if

       if(nOrder==4 .and. UseFlattening .and. .not.DoResChangeOnly)then
          if(UseVolumeIntegral4)then
             call flatten(Prim_VG)
          else
             call flatten(Primitive_VG)
          end if
       end if

    end select  ! nOrder

    if(DoTest)then
       write(*,*) NameSub,' finishing with DoResChangeOnly=', DoResChangeOnly
       if(iDimTest==0 .or. iDimTest==1) &
            write(*,*)'Left,Right(i-1/2),Left,Right(i+1/2)=', &
            LeftState_VX(iVarTest,iTest,jTest,kTest), &
            RightState_VX(iVarTest,iTest,jTest,kTest), &
            LeftState_VX(iVarTest,iTest+1,jTest,kTest), &
            RightState_VX(iVarTest,iTest+1,jTest,kTest)

       if(nDim > 1 .and. (iDimTest==0 .or. iDimTest==2)) &
            write(*,*)'Left,Right(j-1/2),Left,Right(j+1/2)=', &
            LeftState_VY(iVarTest,iTest,jTest,kTest), &
            RightState_VY(iVarTest,iTest,jTest,kTest), &
            LeftState_VY(iVarTest,iTest,jTest+1,kTest), &
            RightState_VY(iVarTest,iTest,jTest+1,kTest)

       if(nDim > 2 .and. (iDimTest==0 .or. iDimTest==3)) &
            write(*,*)'Left,Right(k-1/2),Left,Right(k+1/2)=', &
            LeftState_VZ(iVarTest,iTest,jTest,kTest), &
            RightState_VZ(iVarTest,iTest,jTest,kTest), &
            LeftState_VZ(iVarTest,iTest,jTest,kTest+1), &
            RightState_VZ(iVarTest,iTest,jTest,kTest+1)

    end if

    call test_stop(NameSub, DoTest, iBlock)
  contains
    !==========================================================================
    subroutine limit_var(lMin, lMax, iVar, DoCalcWeightIn)

      ! Switch between various possibilities for high order variable limiter

      integer, intent(in):: lMin, lMax, iVar
      logical, optional, intent(in):: DoCalcWeightIn
      logical:: DoCalcWeight
      !------------------------------------------------------------------------

      DoCalcWeight = .false.
      if(present(DoCalcWeightIn)) DoCalcWeight = DoCalcWeightIn

      if(nOrder == 4)then
         call limiter_ppm4(lMin, lMax, iVar)
      elseif(UseCweno) then
         if (UsePerVarLimiter .or. DoCalcWeight) &
              call calc_cweno_weight(lMin, lMax)
         call limiter_cweno5(lMin, lMax, Cell_I, Cell_I, iVar)
      else
         call limiter_mp(lMin, lMax, Cell_I, Cell_I, iVar)
      end if

    end subroutine limit_var
    !==========================================================================
    subroutine limit_flux(lMin, lMax)

      ! Switch between various possibilities for 5th order flux limiter

      integer, intent(in):: lMin, lMax
      !------------------------------------------------------------------------

      if(UseCweno) then
         call limiter_cweno5(lMin, lMax, Cell_I, Cell2_I)
      else
         call limiter_mp(lMin, lMax, Cell_I, Cell2_I)
      end if

    end subroutine limit_flux
    !==========================================================================
    subroutine logfacex_to_facex(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer:: iVar

      !------------------------------------------------------------------------
      do iVar = 1, nVar
         if(.not.UseLogLimiter_V(iVar))CYCLE

         LeftState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(LeftState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         RightState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(RightState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
      end do

    end subroutine logfacex_to_facex
    !==========================================================================
    subroutine logfacey_to_facey(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer:: iVar

      !------------------------------------------------------------------------
      do iVar = 1, nVar
         if(.not.UseLogLimiter_V(iVar))CYCLE

         LeftState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(LeftState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         RightState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(RightState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
      end do

    end subroutine logfacey_to_facey
    !==========================================================================
    subroutine logfacez_to_facez(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer:: iVar

      !------------------------------------------------------------------------
      do iVar=1,nVar
         if(.not.UseLogLimiter_V(iVar))CYCLE

         LeftState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(LeftState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         RightState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(RightState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
      end do

    end subroutine logfacez_to_facez
    !==========================================================================
    subroutine ratio_to_scalar_facex(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer :: i, j, k

      !------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         LeftState_VX(iVarLimitRatio_I,i,j,k) = &
              LeftState_VX(iVarLimitRatio_I,i,j,k)*LeftState_VX(Rho_,i,j,k)
         RightState_VX(iVarLimitRatio_I,i,j,k) = &
              RightState_VX(iVarLimitRatio_I,i,j,k)*RightState_VX(Rho_,i,j,k)
      end do; end do; end do

    end subroutine ratio_to_scalar_facex
    !==========================================================================
    subroutine ratio_to_scalar_facey(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer :: i, j, k

      !------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         LeftState_VY(iVarLimitRatio_I,i,j,k) = &
              LeftState_VY(iVarLimitRatio_I,i,j,k)*LeftState_VY(Rho_,i,j,k)
         RightState_VY(iVarLimitRatio_I,i,j,k) = &
              RightState_VY(iVarLimitRatio_I,i,j,k)*RightState_VY(Rho_,i,j,k)
      end do; end do; end do

    end subroutine ratio_to_scalar_facey
    !==========================================================================
    subroutine ratio_to_scalar_facez(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer :: i, j, k

      !------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         LeftState_VZ(iVarLimitRatio_I,i,j,k) = &
              LeftState_VZ(iVarLimitRatio_I,i,j,k)*LeftState_VZ(Rho_,i,j,k)
         RightState_VZ(iVarLimitRatio_I,i,j,k) = &
              RightState_VZ(iVarLimitRatio_I,i,j,k)*RightState_VZ(Rho_,i,j,k)
      end do; end do; end do

    end subroutine ratio_to_scalar_facez
    !==========================================================================
    subroutine ptotal_to_p_facex(iMin, iMax, jMin, jMax, kMin, kMax)

      integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax
      integer :: i, j, k

      !------------------------------------------------------------------------
      if(UseElectronPressure)then
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            LeftState_VX(p_,i,j,k) = LeftState_VX(p_,i,j,k) &
                 - LeftState_VX(Pe_,i,j,k)
            RightState_VX(p_,i,j,k) = RightState_VX(p_,i,j,k) &
                 - RightState_VX(Pe_,i,j,k)
         end do; end do; end do
      end if
      if(UseWavePressure)then
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            LeftState_VX(p_,i,j,k) = LeftState_VX(p_,i,j,k) &
                 - (GammaWave-1)*sum(LeftState_VX(WaveFirst_:WaveLast_,i,j,k))
            RightState_VX(p_,i,j,k) = RightState_VX(p_,i,j,k) &
                 - (GammaWave-1)*sum(RightState_VX(WaveFirst_:WaveLast_,i,j,k))
         end do; end do; end do
      end if

    end subroutine ptotal_to_p_facex
    !==========================================================================
    subroutine ptotal_to_p_facey(iMin, iMax, jMin, jMax, kMin, kMax)

      integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax
      integer :: i, j, k

      !------------------------------------------------------------------------
      if(UseElectronPressure)then
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            LeftState_VY(p_,i,j,k) = LeftState_VY(p_,i,j,k) &
                 - LeftState_VY(Pe_,i,j,k)
            RightState_VY(p_,i,j,k) = RightState_VY(p_,i,j,k) &
                 - RightState_VY(Pe_,i,j,k)
         end do; end do; end do
      end if
      if(UseWavePressure)then
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            LeftState_VY(p_,i,j,k) = LeftState_VY(p_,i,j,k) &
                 - (GammaWave-1)*sum(LeftState_VY(WaveFirst_:WaveLast_,i,j,k))
            RightState_VY(p_,i,j,k) = RightState_VY(p_,i,j,k) &
                 - (GammaWave-1)*sum(RightState_VY(WaveFirst_:WaveLast_,i,j,k))
         end do; end do; end do
      end if

    end subroutine ptotal_to_p_facey
    !==========================================================================
    subroutine ptotal_to_p_facez(iMin, iMax, jMin, jMax, kMin, kMax)

      integer, intent(in) :: iMin, iMax, jMin, jMax, kMin, kMax
      integer :: i, j, k

      !------------------------------------------------------------------------
      if(UseElectronPressure)then
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            LeftState_VZ(p_,i,j,k) = LeftState_VZ(p_,i,j,k) &
                 - LeftState_VZ(Pe_,i,j,k)
            RightState_VZ(p_,i,j,k) = RightState_VZ(p_,i,j,k) &
                 - RightState_VZ(Pe_,i,j,k)
         end do; end do; end do
      end if
      if(UseWavePressure)then
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            LeftState_VZ(p_,i,j,k) = LeftState_VZ(p_,i,j,k) &
                 - (GammaWave-1)*sum(LeftState_VZ(WaveFirst_:WaveLast_,i,j,k))
            RightState_VZ(p_,i,j,k) = RightState_VZ(p_,i,j,k) &
                 - (GammaWave-1)*sum(RightState_VZ(WaveFirst_:WaveLast_,i,j,k))
         end do; end do; end do
      end if

    end subroutine ptotal_to_p_facez
    !==========================================================================

    subroutine calc_primitives_boris
      ! momentum is limited

      ! rhoU_Boris = rhoU - ((U x B) x B)/c^2
      !            = rhoU + (U B^2 - B U.B)/c^2
      !            = rhoU*(1+BB/(rho*c2)) - B UdotB/c^2
      !------------------------------------------------------------------------
      Primitive_VG(:,i,j,k) = &
           State_VGB(1:nVar,i,j,k,iBlock)
      BxFull = B0_DG(x_,i,j,k) + Primitive_VG(Bx_,i,j,k)
      ByFull = B0_DG(y_,i,j,k) + Primitive_VG(By_,i,j,k)
      BzFull = B0_DG(z_,i,j,k) + Primitive_VG(Bz_,i,j,k)
      B2Full = BxFull**2 + ByFull**2 + BzFull**2
      RhoC2Inv  = 1/(Primitive_VG(rho_,i,j,k)*C2light)
      uBC2Inv= (Primitive_VG(rhoUx_,i,j,k)*BxFull + &
           Primitive_VG(rhoUy_,i,j,k)*ByFull + &
           Primitive_VG(rhoUz_,i,j,k)*BzFull)*RhoC2Inv
      Ga2Boris= 1 + B2Full*RhoC2Inv

      Primitive_VG(Ux_,i,j,k)= Primitive_VG(rhoUx_,i,j,k)*&
           Ga2Boris - BxFull*uBC2Inv
      Primitive_VG(Uy_,i,j,k)= Primitive_VG(rhoUy_,i,j,k)*&
           Ga2Boris - ByFull*uBC2Inv
      Primitive_VG(Uz_,i,j,k)= Primitive_VG(rhoUz_,i,j,k)*&
           Ga2Boris - BzFull*uBC2Inv

      ! Transform p to Ptotal
      if(UsePtotalLimiter)then
         if(UseElectronPressure)then
            Primitive_VG(p_,i,j,k) = Primitive_VG(p_,i,j,k) &
                 + Primitive_VG(Pe_,i,j,k)
         end if
         if(UseWavePressure)then
            Primitive_VG(p_,i,j,k) = Primitive_VG(p_,i,j,k) &
                 + (GammaWave-1)*sum(Primitive_VG(WaveFirst_:WaveLast_,i,j,k))
         end if
      end if

      if(UseScalarToRhoRatioLtd) Primitive_VG(iVarLimitRatio_I,i,j,k) = &
           Primitive_VG(iVarLimitRatio_I,i,j,k)/Primitive_VG(Rho_,i,j,k)

    end subroutine calc_primitives_boris
    !==========================================================================
    subroutine get_facex_high(iMin,iMax,jMin,jMax,kMin,kMax)
      use ModAdvance, ONLY: FaceDivU_IX

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      real, allocatable, save:: State_VX(:,:,:,:)
      integer:: iVar, iSort
      logical:: IsSmoothIndictor

      !------------------------------------------------------------------------

      if(TypeLimiter == 'no')then
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            LeftState_VX(:,i,j,k) = &
                 c7over12*(Primitive_VG(:,i-1,j,k) + Primitive_VG(:,i,j,k)) - &
                 c1over12*(Primitive_VG(:,i-2,j,k) + Primitive_VG(:,i+1,j,k))
         end do; end do; end do
      else

         do k = kMin, kMax; do j = jMin, jMax
            if(UseLowOrder) then
               Primitive_VI(:,iMin-2:iMax+1)=Primitive_VG(:,iMin-2:iMax+1,j,k)
               if(nLowOrder==2) then
                  ! IsTrueCell needed by limiter_body and ppm4 limiter
                  IsTrueCell_I(iMin-nG:iMax-1+nG) = &
                       true_cell(iMin-nG:iMax-1+nG,j,k,iBlock)
                  ! Get 2nd order limited slopes
                  if(UseTrueCell)then
                     call limiter_body(iMin, iMax, BetaLimiter)
                  else
                     call limiter(iMin, iMax, BetaLimiter)
                  end if
               else
                  dVarLimL_VI = 0; dVarLimR_VI = 0
               endif
               ! Store 1st/2nd order accurate face values
               do i = iMin, iMax
                  ! if(UseLowOrder_I(i))then...
                  LeftState_VX(:,i,j,k) =Primitive_VI(:,i-1)+dVarLimL_VI(:,i-1)
                  RightState_VX(:,i,j,k)=Primitive_VI(:,i)  -dVarLimR_VI(:,i)                 
               end do
            endif

            if(.not.DoInterpolateFlux)then
               if(UseLowOrder) then
                  LowOrderCrit_I(iMin:iMax) = &
                       LowOrderCrit_XB(iMin:iMax,j,k,iBlock) 
               endif

               iVarSmoothLast = 0
               do iSort = 1, nVar
                  if(UseCweno) then
                     iVar = iVarSmoothIndex_I(iSort)
                     iVarSmooth = iVarSmooth_V(iVar)
                     if(iVarSmooth /= iVarSmoothLast) then
                        IsSmoothIndictor = .true.
                        iVarSmoothLast = iVarSmooth
                     else
                        IsSmoothIndictor = .false.
                     endif
                  else
                     iVar = iSort
                  endif

                  ! Copy points along i direction into 1D array
                  Cell_I(iMin-nG:iMax-1+nG) = &
                       Primitive_VG(iVar,iMin-nG:iMax-1+nG,j,k)

                  if(UseLowOrder)then
                     ! Use 2nd order face values where high order is skipped
                     ! where(UseLowOrder_I(iMin:iMax)...
                     FaceL_I(iMin:iMax) = LeftState_VX(iVar,iMin:iMax,j,k)
                     FaceR_I(iMin:iMax) = RightState_VX(iVar,iMin:iMax,j,k)
                  end if


                  call limit_var(iMin, iMax, iVar, &
                       DoCalcWeightIn = IsSmoothIndictor)

                  ! Copy back the results into the 3D arrays
                  LeftState_VX(iVar,iMin:iMax,j,k)  = FaceL_I(iMin:iMax)
                  RightState_VX(iVar,iMin:iMax,j,k) = FaceR_I(iMin:iMax)
               end do
            else
               if(UseCweno)then
                  ! Use Rho as the smooth indicator
                  Cell_I(iMin-nG:iMax-1+nG) = &
                       Primitive_VG(Rho_,iMin-nG:iMax-1+nG,j,k)
                  call calc_cweno_weight(iMin, iMax)
               end if

               ! Get face value for Ux
               ! Copy points along i direction into 1D array
               Cell_I(iMin-nG:iMax-1+nG) = &
                    Primitive_VG(Ux_,iMin-nG:iMax-1+nG,j,k)
               call limit_var(iMin, iMax, Ux_)
               uDotArea_XI(iMin:iMax,j,k,1) = CellFace_DB(1,iBlock) &
                    *0.5*(FaceL_I(iMin:iMax) + FaceR_I(iMin:iMax))

               ! Interpolate cell centered split fluxes to the face
               do iFlux = 1, nVar + nFluid
                  ! Copy left fluxes along i direction into 1D array
                  Cell_I(iMin-nG:iMax-1+nG) = &
                       FluxLeft_VGD(iFlux,iMin-nG:iMax-1+nG,j,k,x_)

                  Cell2_I(iMin-nG:iMax-1+nG) = &
                       FluxRight_VGD(iFlux,iMin-nG:iMax-1+nG,j,k,x_)

                  call limit_flux(iMin, iMax)

                  Flux_VX(iFlux,iMin:iMax,j,k) = &
                       FaceL_I(iMin:iMax) + FaceR_I(iMin:iMax)
               end do
            endif
         end do; end do

      end if

      if(nOrder == 4 .and. UseFaceIntegral4)then
         if(.not.allocated(State_VX)) allocate( &
              State_VX(nVar,nI+1,jMinFace2:jMaxFace2,kMinFace2:kMaxFace2))

         ! Convert from face averaged to face centered variables (eq 18)
         State_VX = LeftState_VX
         do k = kMinFace,kMaxFace; do j = jMinFace,jMaxFace; do i = iMin,iMax
            Laplace_V = -2*(nDim-1)*State_VX(:,i,j,k) &
                 + State_VX(:,i,j-1,k) + State_VX(:,i,j+1,k)
            if(nK>1) Laplace_V = Laplace_V &
                 + State_VX(:,i,j,k-1) + State_VX(:,i,j,k+1)
            LeftState_VX(:,i,j,k) = State_VX(:,i,j,k) - c24th*Laplace_V
            ! Keep positivity
            do iVar = 1, nVar
               if(DefaultState_V(iVar)>0 .and. LeftState_VX(iVar,i,j,k)<0) &
                    LeftState_VX(iVar,i,j,k) = State_VX(iVar,i,j,k)
            end do
         end do; end do; end do
         if(TypeLimiter /= 'no')then
            State_VX = RightState_VX
            do k=kMinFace,kMaxFace; do j=jMinFace,jMaxFace; do i=iMin,iMax
               Laplace_V = -2*(nDim-1)*State_VX(:,i,j,k) &
                    + State_VX(:,i,j-1,k) + State_VX(:,i,j+1,k)
               if(nK>1) Laplace_V = Laplace_V &
                    + State_VX(:,i,j,k-1) + State_VX(:,i,j,k+1)
               RightState_VX(:,i,j,k) = State_VX(:,i,j,k) - c24th*Laplace_V

               ! Keep positivity
               do iVar = 1, nVar
                  if(DefaultState_V(iVar)>0.and.RightState_VX(iVar,i,j,k)<0) &
                       RightState_VX(iVar,i,j,k) = State_VX(iVar,i,j,k)
               end do
            end do; end do; end do
         end if
      end if

      if(TypeLimiter == 'no') RightState_VX = LeftState_VX

      if(DoLimitMomentum)call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceX(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facex_high
    !==========================================================================
    subroutine get_facey_high(iMin,iMax,jMin,jMax,kMin,kMax)
      use ModAdvance, ONLY: FaceDivU_IY
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax

      real, allocatable, save:: State_VY(:,:,:,:)
      integer:: iVar, iSort
      logical:: IsSmoothIndictor
      !------------------------------------------------------------------------

      if(TypeLimiter == 'no')then
         do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
            LeftState_VY(:,i,j,k) = &
                 c7over12*(Primitive_VG(:,i,j-1,k) + Primitive_VG(:,i,j,k)) - &
                 c1over12*(Primitive_VG(:,i,j-2,k) + Primitive_VG(:,i,j+1,k))
         end do; end do; end do
      else
         do k = kMin, kMax; do i = iMin, iMax
            if(UseLowOrder) then
               Primitive_VI(:,jMin-2:jMax+1)=Primitive_VG(:,i,jMin-2:jMax+1,k)

               if(nLowOrder==2) then
                  IsTrueCell_I(jMin-nG:jMax-1+nG) = &
                       true_cell(i,jMin-nG:jMax-1+nG,k,iBlock)
                  if(UseTrueCell)then
                     call limiter_body(jMin, jMax, BetaLimiter)
                  else
                     call limiter(jMin, jMax, BetaLimiter)
                  end if
               else
                  dVarLimL_VI = 0; dVarLimR_VI = 0
               endif
               do j = jMin, jMax
                  LeftState_VY(:,i,j,k) =Primitive_VI(:,j-1)+dVarLimL_VI(:,j-1)
                  RightState_VY(:,i,j,k)=Primitive_VI(:,j)  -dVarLimR_VI(:,j)
               end do
            endif

            if(.not.DoInterpolateFlux)then
               if(UseLowOrder) then
                  LowOrderCrit_I(jMin:jMax) = &
                       LowOrderCrit_YB(i,jMin:jMax,k,iBlock)
               endif

               iVarSmoothLast = 0
               do iSort = 1, nVar
                  if(UseCweno) then
                     ! The variables use the same smooth indicator are
                     ! calculated one by one. And the smooth indicator
                     ! itself is calculated first.
                     iVar = iVarSmoothIndex_I(iSort)
                     iVarSmooth = iVarSmooth_V(iVar)
                     if(iVarSmooth /= iVarSmoothLast) then
                        IsSmoothIndictor = .true.
                        iVarSmoothLast = iVarSmooth
                     else
                        IsSmoothIndictor = .false.
                     endif
                  else
                     iVar = iSort
                  endif

                  ! Copy points along j direction into 1D array
                  Cell_I(jMin-nG:jMax-1+nG) = &
                       Primitive_VG(iVar,i,jMin-nG:jMax-1+nG,k)

                  if(UseLowOrder)then
                     ! Use 2nd order face values where high order is skipped
                     FaceL_I(jMin:jMax) = LeftState_VY(iVar,i,jMin:jMax,k)
                     FaceR_I(jMin:jMax) = RightState_VY(iVar,i,jMin:jMax,k)
                  end if

                  call limit_var(jMin, jMax, iVar, &
                       DoCalcWeightIn = IsSmoothIndictor)

                  ! Copy back the results into the 3D arrays
                  LeftState_VY(iVar,i,jMin:jMax,k)  = FaceL_I(jMin:jMax)
                  RightState_VY(iVar,i,jMin:jMax,k) = FaceR_I(jMin:jMax)
               end do
            else
               if(UseCweno)then
                  ! Use Rho as the smooth indicator
                  Cell_I(jMin-nG:jMax-1+nG) = &
                       Primitive_VG(Rho_,i,jMin-nG:jMax-1+nG,k)
                  call calc_cweno_weight(jMin, jMax)
               end if

               ! Get face value for Uy
               ! Copy points along i direction into 1D array
               Cell_I(jMin-nG:jMax-1+nG) = &
                    Primitive_VG(Uy_,i,jMin-nG:jMax-1+nG,k)
               call limit_var(jMin, jMax, Uy_)
               uDotArea_YI(i,jMin:jMax,k,1) = CellFace_DB(2,iBlock) &
                    *0.5*(FaceL_I(jMin:jMax) + FaceR_I(jMin:jMax))

               ! Interpolate cell centered split fluxes to the face
               do iFlux = 1, nVar + nFluid
                  ! Copy left fluxes along i direction into 1D array
                  Cell_I(jMin-nG:jMax-1+nG) = &
                       FluxLeft_VGD(iFlux,i,jMin-nG:jMax-1+nG,k,y_)

                  Cell2_I(jMin-nG:jMax-1+nG) = &
                       FluxRight_VGD(iFlux,i,jMin-nG:jMax-1+nG,k,y_)

                  call limit_flux(jMin, jMax)

                  Flux_VY(iFlux,i,jMin:jMax,k) = &
                       FaceL_I(jMin:jMax) + FaceR_I(jMin:jMax)
               end do
            endif

         end do; end do
      end if

      if(nOrder == 4 .and. UseFaceIntegral4)then
         if(.not.allocated(State_VY)) allocate( &
              State_VY(nVar,iMinFace2:iMaxFace2,nJ+1,kMinFace2:kMaxFace2))

         ! Convert from face averaged to face centered variables (eq 18)
         State_VY = LeftState_VY
         do k = kMinFace,kMaxFace; do j = jMin,jMax; do i = iMinFace,iMaxFace;
            Laplace_V = -2*(nDim-1)*State_VY(:,i,j,k) &
                 + State_VY(:,i-1,j,k) + State_VY(:,i+1,j,k)
            if(nK>1) Laplace_V = Laplace_V &
                 + State_VY(:,i,j,k-1) + State_VY(:,i,j,k+1)
            LeftState_VY(:,i,j,k) = State_VY(:,i,j,k) - c24th*Laplace_V
            ! Keep positivity
            do iVar = 1, nVar
               if(DefaultState_V(iVar)>0.and.LeftState_VY(iVar,i,j,k)<0) &
                    LeftState_VY(iVar,i,j,k) = State_VY(iVar,i,j,k)
            end do
         end do; end do; end do
         if(TypeLimiter /= 'no')then
            State_VY = RightState_VY
            do k=kMinFace,kMaxFace; do j=jMin,jMax; do i=iMinFace,iMaxFace;
               Laplace_V = -2*(nDim-1)*State_VY(:,i,j,k) &
                    + State_VY(:,i-1,j,k) + State_VY(:,i+1,j,k)
               if(nK>1) Laplace_V = Laplace_V &
                    + State_VY(:,i,j,k-1) + State_VY(:,i,j,k+1)
               RightState_VY(:,i,j,k) = State_VY(:,i,j,k) - c24th*Laplace_V
               ! Keep positivity
               do iVar = 1, nVar
                  if(DefaultState_V(iVar)>0.and.RightState_VY(iVar,i,j,k)<0) &
                       RightState_VY(iVar,i,j,k) = State_VY(iVar,i,j,k)
               end do
            end do; end do; end do
         end if
      end if

      if(TypeLimiter == 'no') RightState_VY = LeftState_VY
      if(DoLimitMomentum)call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceY(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facey_high
    !==========================================================================
    subroutine get_facez_high(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax

      real, allocatable, save:: State_VZ(:,:,:,:)
      integer:: iVar, iSort
      logical:: IsSmoothIndictor
      !------------------------------------------------------------------------

      if(TypeLimiter == 'no')then
         do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
            LeftState_VZ(:,i,j,k) = &
                 c7over12*(Primitive_VG(:,i,j,k-1) + Primitive_VG(:,i,j,k)) - &
                 c1over12*(Primitive_VG(:,i,j,k-2) + Primitive_VG(:,i,j,k+1))

            RightState_VZ(:,i,j,k)=LeftState_VZ(:,i,j,k)
         end do; end do; end do
      else
         do j = jMin, jMax; do i = iMin, iMax

            if(UseLowOrder)then
               Primitive_VI(:,kMin-2:kMax+1)=Primitive_VG(:,i,j,kMin-2:kMax+1)

               if(nLowOrder==2) then
                  IsTrueCell_I(kMin-nG:kMax-1+nG) = &
                       true_cell(i,j,kMin-nG:kMax-1+nG,iBlock)
                  if(UseTrueCell)then
                     call limiter_body(kMin, kMax, BetaLimiter)
                  else
                     call limiter(kMin, kMax, BetaLimiter)
                  end if
               else
                  dVarLimL_VI = 0; dVarLimR_VI = 0
               endif
               do k = kMin, kMax
                  LeftState_VZ(:,i,j,k) =Primitive_VI(:,k-1)+dVarLimL_VI(:,k-1)
                  RightState_VZ(:,i,j,k)=Primitive_VI(:,k)  -dVarLimR_VI(:,k)                  
               end do
            end if

            if(.not.DoInterpolateFlux)then
               if(UseLowOrder) then
                  LowOrderCrit_I(kMin:kMax) = &
                       LowOrderCrit_ZB(i,j,kMin:kMax,iBlock) 

               endif

               iVarSmoothLast = 0
               do iSort = 1, nVar
                  if(UseCweno) then
                     iVar = iVarSmoothIndex_I(iSort)
                     iVarSmooth = iVarSmooth_V(iVar)
                     if(iVarSmooth /= iVarSmoothLast) then
                        IsSmoothIndictor = .true.
                        iVarSmoothLast = iVarSmooth
                     else
                        IsSmoothIndictor = .false.
                     endif
                  else
                     iVar = iSort
                  endif

                  ! Copy points along k direction into 1D array
                  Cell_I(kMin-nG:kMax-1+nG) = &
                       Primitive_VG(iVar,i,j,kMin-nG:kMax-1+nG)

                  if(UseLowOrder)then
                     ! Use 2nd order face values where high order is skipped
                     FaceL_I(kMin:kMax) = LeftState_VZ(iVar,i,j,kMin:kMax)
                     FaceR_I(kMin:kMax) = RightState_VZ(iVar,i,j,kMin:kMax)
                  end if

                  call limit_var(kMin, kMax, iVar, &
                       DoCalcWeightIn = IsSmoothIndictor)

                  ! Copy back the results into the 3D arrays
                  LeftState_VZ(iVar,i,j,kMin:kMax)  = FaceL_I(kMin:kMax)
                  RightState_VZ(iVar,i,j,kMin:kMax) = FaceR_I(kMin:kMax)
               end do
            else
               if (UseCweno) then
                  Cell_I(kMin-nG:kMax-1+nG) = &
                       Primitive_VG(Rho_,i,j,kMin-nG:kMax-1+nG)
                  call calc_cweno_weight(kMin, kMax)
               end if

               ! Copy points along i direction into 1D array
               Cell_I(kMin-nG:kMax-1+nG) = &
                    Primitive_VG(Uz_,i,j,kMin-nG:kMax-1+nG)

               call limit_var(kMin, kMax, Uz_)
               uDotArea_ZI(i,j,kMin:kMax,1) = CellFace_DB(3,iBlock) &
                    *0.5*(FaceL_I(kMin:kMax) + FaceR_I(kMin:kMax))

               ! Interpolate cell centered split fluxes to the face
               do iFlux = 1, nVar + nFluid
                  ! Copy left fluxes along i direction into 1D array
                  Cell_I(kMin-nG:kMax-1+nG) = &
                       FluxLeft_VGD(iFlux,i,j,kMin-nG:kMax-1+nG,z_)

                  Cell2_I(kMin-nG:kMax-1+nG) = &
                       FluxRight_VGD(iFlux,i,j,kMin-nG:kMax-1+nG,z_)

                  call limit_flux(kMin, kMax)

                  Flux_VZ(iFlux,i,j,kMin:kMax) = &
                       FaceL_I(kMin:kMax) + FaceR_I(kMin:kMax)
               end do
            endif
         end do; end do
      end if

      if(nOrder == 4 .and. UseFaceIntegral4)then
         if(.not.allocated(State_VZ)) allocate( &
              State_VZ(nVar,iMinFace2:iMaxFace2,jMinFace2:jMaxFace2,nK+1))

         ! Convert from face averaged to face centered variables (eq 18)
         State_VZ = LeftState_VZ
         do k = kMin,kMax; do j = jMinFace,jMaxFace; do i = iMinFace,iMaxFace
            Laplace_V = -4*State_VZ(:,i,j,k) &
                 + State_VZ(:,i-1,j,k) + State_VZ(:,i+1,j,k) &
                 + State_VZ(:,i,j-1,k) + State_VZ(:,i,j+1,k)
            LeftState_VZ(:,i,j,k) = State_VZ(:,i,j,k) - c24th*Laplace_V
            ! Keep positivity
            do iVar = 1, nVar
               if(DefaultState_V(iVar)>0 .and. LeftState_VZ(iVar,i,j,k)<0) &
                    LeftState_VZ(iVar,i,j,k) = State_VZ(iVar,i,j,k)
            end do
         end do; end do; end do
         if(TypeLimiter /= 'no')then
            State_VZ = RightState_VZ
            do k=kMin,kMax; do j=jMinFace,jMaxFace; do i=iMinFace,iMaxFace
               Laplace_V = -4*State_VZ(:,i,j,k) &
                    + State_VZ(:,i-1,j,k) + State_VZ(:,i+1,j,k) &
                    + State_VZ(:,i,j-1,k) + State_VZ(:,i,j+1,k)
               RightState_VZ(:,i,j,k) = State_VZ(:,i,j,k) - c24th*Laplace_V
               ! Keep positivity
               do iVar = 1, nVar
                  if(DefaultState_V(iVar)>0.and.RightState_VZ(iVar,i,j,k)<0) &
                       RightState_VZ(iVar,i,j,k) = State_VZ(iVar,i,j,k)
               end do
            end do; end do; end do
         end if
      end if

      if(DoLimitMomentum)call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceZ(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facez_high
    !==========================================================================

    subroutine get_facex_first(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         LeftState_VX(:,i,j,k)=Primitive_VG(:,i-1,j,k)
         RightState_VX(:,i,j,k)=Primitive_VG(:,i,j,k)
      end do; end do; end do

      if(DoLimitMomentum)call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceX(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facex_first
    !==========================================================================
    subroutine get_facey_first(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         LeftState_VY(:,i,j,k)=Primitive_VG(:,i,j-1,k)
         RightState_VY(:,i,j,k)=Primitive_VG(:,i,j,k)
      end do; end do; end do

      if(DoLimitMomentum) call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceY(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facey_first
    !==========================================================================
    subroutine get_facez_first(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         LeftState_VZ(:,i,j,k)=Primitive_VG(:,i,j,k-1)
         RightState_VZ(:,i,j,k)=Primitive_VG(:,i,j,k)
      end do; end do; end do

      if(DoLimitMomentum)call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceZ(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facez_first
    !==========================================================================
    subroutine calc_primitives_MHD
      integer:: iVar
      !------------------------------------------------------------------------
      Primitive_VG(:,i,j,k) = State_VGB(1:nVar,i,j,k,iBlock)
      RhoInv = 1/Primitive_VG(Rho_,i,j,k)
      Primitive_VG(Ux_:Uz_,i,j,k)=RhoInv*Primitive_VG(RhoUx_:RhoUz_,i,j,k)
      do iFluid = 2, nFluid
         iRho = iRho_I(iFluid); iUx = iUx_I(iFluid); iUz = iUz_I(iFluid)
         RhoInv = 1/Primitive_VG(iRho,i,j,k)
         Primitive_VG(iUx:iUz,i,j,k)=RhoInv*Primitive_VG(iUx:iUz,i,j,k)
      end do

      if(UsePtotalLimiter)then
         ! Transform p to Ptotal
         if(UseElectronPressure)then
            Primitive_VG(p_,i,j,k) = Primitive_VG(p_,i,j,k) &
                 + Primitive_VG(Pe_,i,j,k)
         end if
         if(UseWavePressure)then
            Primitive_VG(p_,i,j,k) = Primitive_VG(p_,i,j,k) &
                 + (GammaWave-1)*sum(Primitive_VG(WaveFirst_:WaveLast_,i,j,k))
         end if
      end if

      if(UseScalarToRhoRatioLtd) Primitive_VG(iVarLimitRatio_I,i,j,k) = &
           RhoInv*Primitive_VG(iVarLimitRatio_I,i,j,k)

      if(UseLogLimiter)then
         do iVar = 1, nVar
            if(UseLogLimiter_V(iVar)) &
                 Primitive_VG(iVar,i,j,k) = log(Primitive_VG(iVar,i,j,k))
         end do
      end if

    end subroutine calc_primitives_MHD
    !==========================================================================
    subroutine BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      ! Convert face centered Boris momenta to MHD velocities

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      ! U_Boris=rhoU_Boris/rho
      ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

         ! Left face values
         RhoInv = 1/LeftState_VX(rho_,i,j,k)
         if(UseB0)then
            BxFull = B0_DX(x_,i,j,k) + LeftState_VX(Bx_,i,j,k)
            ByFull = B0_DX(y_,i,j,k) + LeftState_VX(By_,i,j,k)
            BzFull = B0_DX(z_,i,j,k) + LeftState_VX(Bz_,i,j,k)
         else
            BxFull = LeftState_VX(Bx_,i,j,k)
            ByFull = LeftState_VX(By_,i,j,k)
            BzFull = LeftState_VX(Bz_,i,j,k)
         end if
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  = InvClight2*RhoInv
         LeftState_VX(Ux_,i,j,k)=LeftState_VX(Ux_,i,j,k)*RhoInv
         LeftState_VX(Uy_,i,j,k)=LeftState_VX(Uy_,i,j,k)*RhoInv
         LeftState_VX(Uz_,i,j,k)=LeftState_VX(Uz_,i,j,k)*RhoInv
         uBC2Inv= (LeftState_VX(Ux_,i,j,k)*BxFull + &
              LeftState_VX(Uy_,i,j,k)*ByFull + &
              LeftState_VX(Uz_,i,j,k)*BzFull)*RhoC2Inv

         ! gammaA^2 = 1/[1+BB/(rho c^2)]
         Ga2Boris= 1/(1 + B2Full*RhoC2Inv)

         LeftState_VX(Ux_,i,j,k) = &
              Ga2Boris * (LeftState_VX(Ux_,i,j,k)+uBC2Inv*BxFull)
         LeftState_VX(Uy_,i,j,k) = &
              Ga2Boris * (LeftState_VX(Uy_,i,j,k)+uBC2Inv*ByFull)
         LeftState_VX(Uz_,i,j,k) = &
              Ga2Boris * (LeftState_VX(Uz_,i,j,k)+uBC2Inv*BzFull)

         ! Right face values
         RhoInv = 1/RightState_VX(rho_,i,j,k)
         if(UseB0)then
            BxFull = B0_DX(x_,i,j,k) + RightState_VX(Bx_,i,j,k)
            ByFull = B0_DX(y_,i,j,k) + RightState_VX(By_,i,j,k)
            BzFull = B0_DX(z_,i,j,k) + RightState_VX(Bz_,i,j,k)
         else
            BxFull = RightState_VX(Bx_,i,j,k)
            ByFull = RightState_VX(By_,i,j,k)
            BzFull = RightState_VX(Bz_,i,j,k)
         end if
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  =InvClight2*RhoInv
         RightState_VX(Ux_,i,j,k)=RightState_VX(Ux_,i,j,k)*RhoInv
         RightState_VX(Uy_,i,j,k)=RightState_VX(Uy_,i,j,k)*RhoInv
         RightState_VX(Uz_,i,j,k)=RightState_VX(Uz_,i,j,k)*RhoInv
         uBC2Inv= (RightState_VX(Ux_,i,j,k)*BxFull + &
              RightState_VX(Uy_,i,j,k)*ByFull + &
              RightState_VX(Uz_,i,j,k)*BzFull)*RhoC2Inv

         ! gammaA^2 = 1/[1+BB/(rho c^2)]
         Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

         RightState_VX(Ux_,i,j,k) = &
              Ga2Boris * (RightState_VX(Ux_,i,j,k)+uBC2Inv*BxFull)
         RightState_VX(Uy_,i,j,k) = &
              Ga2Boris * (RightState_VX(Uy_,i,j,k)+uBC2Inv*ByFull)
         RightState_VX(Uz_,i,j,k) = &
              Ga2Boris * (RightState_VX(Uz_,i,j,k)+uBC2Inv*BzFull)

      end do; end do; end do
    end subroutine BorisFaceXtoMHD
    !==========================================================================
    subroutine BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)
      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax

      ! U_Boris=rhoU_Boris/rho
      ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

      !------------------------------------------------------------------------
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

         ! Left face values
         RhoInv = 1/LeftState_VY(rho_,i,j,k)
         if(UseB0)then
            BxFull = B0_DY(x_,i,j,k) + LeftState_VY(Bx_,i,j,k)
            ByFull = B0_DY(y_,i,j,k) + LeftState_VY(By_,i,j,k)
            BzFull = B0_DY(z_,i,j,k) + LeftState_VY(Bz_,i,j,k)
         else
            BxFull = LeftState_VY(Bx_,i,j,k)
            ByFull = LeftState_VY(By_,i,j,k)
            BzFull = LeftState_VY(Bz_,i,j,k)
         end if
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  =InvClight2*RhoInv
         LeftState_VY(Ux_,i,j,k)=LeftState_VY(Ux_,i,j,k)*RhoInv
         LeftState_VY(Uy_,i,j,k)=LeftState_VY(Uy_,i,j,k)*RhoInv
         LeftState_VY(Uz_,i,j,k)=LeftState_VY(Uz_,i,j,k)*RhoInv
         uBC2Inv= (LeftState_VY(Ux_,i,j,k)*BxFull + &
              LeftState_VY(Uy_,i,j,k)*ByFull + &
              LeftState_VY(Uz_,i,j,k)*BzFull)*RhoC2Inv

         ! gammaA^2 = 1/[1+BB/(rho c^2)]
         Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

         LeftState_VY(Ux_,i,j,k) = &
              Ga2Boris * (LeftState_VY(Ux_,i,j,k)+uBC2Inv*BxFull)
         LeftState_VY(Uy_,i,j,k) = &
              Ga2Boris * (LeftState_VY(Uy_,i,j,k)+uBC2Inv*ByFull)
         LeftState_VY(Uz_,i,j,k) = &
              Ga2Boris * (LeftState_VY(Uz_,i,j,k)+uBC2Inv*BzFull)

         ! Right face values
         RhoInv = 1/RightState_VY(rho_,i,j,k)
         if(UseB0)then
            BxFull = B0_DY(x_,i,j,k) + RightState_VY(Bx_,i,j,k)
            ByFull = B0_DY(y_,i,j,k) + RightState_VY(By_,i,j,k)
            BzFull = B0_DY(z_,i,j,k) + RightState_VY(Bz_,i,j,k)
         else
            BxFull = RightState_VY(Bx_,i,j,k)
            ByFull = RightState_VY(By_,i,j,k)
            BzFull = RightState_VY(Bz_,i,j,k)
         end if
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv=InvClight2*RhoInv
         RightState_VY(Ux_,i,j,k)=RightState_VY(Ux_,i,j,k)*RhoInv
         RightState_VY(Uy_,i,j,k)=RightState_VY(Uy_,i,j,k)*RhoInv
         RightState_VY(Uz_,i,j,k)=RightState_VY(Uz_,i,j,k)*RhoInv
         uBC2Inv= (RightState_VY(Ux_,i,j,k)*BxFull + &
              RightState_VY(Uy_,i,j,k)*ByFull + &
              RightState_VY(Uz_,i,j,k)*BzFull)*RhoC2Inv

         ! gammaA^2 = 1/[1+BB/(rho c^2)]
         Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

         RightState_VY(Ux_,i,j,k) = &
              Ga2Boris * (RightState_VY(Ux_,i,j,k)+uBC2Inv*BxFull)
         RightState_VY(Uy_,i,j,k) = &
              Ga2Boris * (RightState_VY(Uy_,i,j,k)+uBC2Inv*ByFull)
         RightState_VY(Uz_,i,j,k) = &
              Ga2Boris * (RightState_VY(Uz_,i,j,k)+uBC2Inv*BzFull)
      end do; end do; end do

    end subroutine BorisFaceYtoMHD
    !==========================================================================
    subroutine BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)
      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax

      ! Convert face centered Boris momenta/rho to MHD velocities
      !------------------------------------------------------------------------
      ! U_Boris=rhoU_Boris/rho
      ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

         ! Left face values
         RhoInv = 1/LeftState_VZ(rho_,i,j,k)
         if(UseB0)then
            BxFull = B0_DZ(x_,i,j,k) + LeftState_VZ(Bx_,i,j,k)
            ByFull = B0_DZ(y_,i,j,k) + LeftState_VZ(By_,i,j,k)
            BzFull = B0_DZ(z_,i,j,k) + LeftState_VZ(Bz_,i,j,k)
         else
            BxFull = LeftState_VZ(Bx_,i,j,k)
            ByFull = LeftState_VZ(By_,i,j,k)
            BzFull = LeftState_VZ(Bz_,i,j,k)
         end if
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  =InvClight2*RhoInv
         LeftState_VZ(Ux_,i,j,k)=LeftState_VZ(Ux_,i,j,k)*RhoInv
         LeftState_VZ(Uy_,i,j,k)=LeftState_VZ(Uy_,i,j,k)*RhoInv
         LeftState_VZ(Uz_,i,j,k)=LeftState_VZ(Uz_,i,j,k)*RhoInv
         uBC2Inv= (LeftState_VZ(Ux_,i,j,k)*BxFull + &
              LeftState_VZ(Uy_,i,j,k)*ByFull + &
              LeftState_VZ(Uz_,i,j,k)*BzFull)*RhoC2Inv

         ! gammaA^2 = 1/[1+BB/(rho c^2)]
         Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

         LeftState_VZ(Ux_,i,j,k) = &
              Ga2Boris * (LeftState_VZ(Ux_,i,j,k)+uBC2Inv*BxFull)
         LeftState_VZ(Uy_,i,j,k) = &
              Ga2Boris * (LeftState_VZ(Uy_,i,j,k)+uBC2Inv*ByFull)
         LeftState_VZ(Uz_,i,j,k) = &
              Ga2Boris * (LeftState_VZ(Uz_,i,j,k)+uBC2Inv*BzFull)

         ! Right face values
         RhoInv = 1/RightState_VZ(rho_,i,j,k)
         if(UseB0)then
            BxFull = B0_DZ(x_,i,j,k) + RightState_VZ(Bx_,i,j,k)
            ByFull = B0_DZ(y_,i,j,k) + RightState_VZ(By_,i,j,k)
            BzFull = B0_DZ(z_,i,j,k) + RightState_VZ(Bz_,i,j,k)
         else
            BxFull = RightState_VZ(Bx_,i,j,k)
            ByFull = RightState_VZ(By_,i,j,k)
            BzFull = RightState_VZ(Bz_,i,j,k)
         end if
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  =InvClight2*RhoInv
         RightState_VZ(Ux_,i,j,k)=RightState_VZ(Ux_,i,j,k)*RhoInv
         RightState_VZ(Uy_,i,j,k)=RightState_VZ(Uy_,i,j,k)*RhoInv
         RightState_VZ(Uz_,i,j,k)=RightState_VZ(Uz_,i,j,k)*RhoInv
         uBC2Inv= (RightState_VZ(Ux_,i,j,k)*BxFull + &
              RightState_VZ(Uy_,i,j,k)*ByFull + &
              RightState_VZ(Uz_,i,j,k)*BzFull)*RhoC2Inv

         ! gammaA^2 = 1/[1+BB/(rho c^2)]
         Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

         RightState_VZ(Ux_,i,j,k) = &
              Ga2Boris * (RightState_VZ(Ux_,i,j,k)+uBC2Inv*BxFull)
         RightState_VZ(Uy_,i,j,k) = &
              Ga2Boris * (RightState_VZ(Uy_,i,j,k)+uBC2Inv*ByFull)
         RightState_VZ(Uz_,i,j,k) = &
              Ga2Boris * (RightState_VZ(Uz_,i,j,k)+uBC2Inv*BzFull)
      end do; end do; end do

    end subroutine BorisFaceZtoMHD
    !==========================================================================
    subroutine get_face_accurate3d(iSideIn)
      integer, intent(in):: iSideIn

      !------------------------------------------------------------------------
      select case(iSideIn)
      case(1)
         do k=1,nK,2; do j=1,nJ,2
            if(  all(true_cell(-1:2,j:j+1,k:k+1,iBlock)) .and. &
                 all(true_cell(0,j-2:j+3,k-2:k+3,iBlock)) ) then
               call accurate_reschange3d(&
                    Coarse2_V    =    Primitive_VG(:,-1,j,k)           ,&
                    Coarse1_VII  =    Primitive_VG(:, 0,j-2:j+3,k-2:k+3),&
                    Fine1_VII    =    Primitive_VG(:, 1,j:j+1,k:k+1)   ,&
                    Fine2_VII    =    Primitive_VG(:, 2,j:j+1,k:k+1)   ,&
                    CoarseToFineF_VII= LeftState_VX(:,1,j:j+1,k:k+1)   ,&
                    FineToCoarseF_VII=RightState_VX(:,1,j:j+1,k:k+1)   ,&
                    FineF_VII        = LeftState_VX(:,2,j:j+1,k:k+1))
            else
               call tvd_reschange_body(                                 &
                    Coarse2_V    =    Primitive_VG(:,-1,j,k)           ,&
                    Coarse1_V    =    Primitive_VG(:, 0,j,k)           ,&
                    Fine1_VII    =    Primitive_VG(:, 1,j:j+1,k:k+1)   ,&
                    Fine2_VII    =    Primitive_VG(:, 2,j:j+1,k:k+1)   ,&
                    CoarseToFineF_VII= LeftState_VX(:,1,j:j+1,k:k+1)   ,&
                    FineToCoarseF_VII=RightState_VX(:,1,j:j+1,k:k+1)   ,&
                    FineF_VII        = LeftState_VX(:,2,j:j+1,k:k+1)   ,&
                    IsTrueCoarse2    = true_cell(-1,j,k,iBlock)        ,&
                    IsTrueCoarse1    = true_cell( 0,j,k,iBlock)        ,&
                    IsTrueFine1  =all(true_cell( 1,j:j+1,k:k+1,iBlock)),&
                    IsTrueFine2_II   = true_cell( 2,j:j+1,k:k+1,iBlock))
            end if
         end do; end do
      case(2)
         do k=1,nK,2; do j=1,nJ,2
            if(  all(true_cell(nI-1:nI+2,j:j+1,k:k+1,iBlock)).and. &
                 all(true_cell(nI+1,j-2:j+3,k-2:k+3,iBlock)) ) then
               call accurate_reschange3d(&
                    Coarse2_V    =    Primitive_VG(:, nI+2,j,k)         ,&
                    Coarse1_VII  =    Primitive_VG(:, nI+1,j-2:j+3,k-2:k+3) ,&
                    Fine1_VII    =    Primitive_VG(:, nI,j:j+1,k:k+1)  ,&
                    Fine2_VII    =    Primitive_VG(:, nI-1,j:j+1,k:k+1),&
                    CoarseToFineF_VII=RightState_VX(:,nI+1,j:j+1,k:k+1),&
                    FineToCoarseF_VII=LeftState_VX(:,nI+1,j:j+1,k:k+1) ,&
                    FineF_VII       =RightState_VX(:,nI,j:j+1,k:k+1))
            else
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,nI+2,j,k)         ,&
                    Coarse1_V    =    Primitive_VG(:, nI+1,j,k)        ,&
                    Fine1_VII    =    Primitive_VG(:, nI,j:j+1,k:k+1)  ,&
                    Fine2_VII    =    Primitive_VG(:, nI-1,j:j+1,k:k+1),&
                    CoarseToFineF_VII=RightState_VX(:,nI+1,j:j+1,k:k+1),&
                    FineToCoarseF_VII=LeftState_VX(:,nI+1,j:j+1,k:k+1) ,&
                    FineF_VII        =RightState_VX(:,nI,j:j+1,k:k+1)  ,&
                    IsTrueCoarse2    = true_cell(nI+2,j,k,iBlock)      ,&
                    IsTrueCoarse1    = true_cell(nI+1,j,k,iBlock)      ,&
                    IsTrueFine1  =all(true_cell(nI,j:j+1,k:k+1,iBlock)),&
                    IsTrueFine2_II      =true_cell(nI-1,j:j+1,k:k+1,iBlock))
            end if
         end do; end do
      case(3)
         do k=1,nK,2; do i=1,nI,2
            if(  all(true_cell(i:i+1,-1:2,k:k+1,iBlock)) .and. &
                 all(true_cell(i-2:i+3,0,k-2:k+3,iBlock)) ) then
               call accurate_reschange3d(&
                    Coarse2_V    =    Primitive_VG(:,i,-1,k)           ,&
                    Coarse1_VII  =    Primitive_VG(:,i-2:i+3,0,k-2:k+3),&
                    Fine1_VII    =    Primitive_VG(:,i:i+1, 1,k:k+1)   ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1, 2,k:k+1)   ,&
                    CoarseToFineF_VII= LeftState_VY(:,i:i+1,1,k:k+1)   ,&
                    FineToCoarseF_VII=RightState_VY(:,i:i+1,1,k:k+1)   ,&
                    FineF_VII       = LeftState_VY(:,i:i+1,2,k:k+1))
            else
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,i,-1,k)           ,&
                    Coarse1_V    =    Primitive_VG(:,i, 0,k)           ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1, 1,k:k+1)   ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1, 2,k:k+1)   ,&
                    CoarseToFineF_VII= LeftState_VY(:,i:i+1,1,k:k+1)   ,&
                    FineToCoarseF_VII=RightState_VY(:,i:i+1,1,k:k+1)   ,&
                    FineF_VII        = LeftState_VY(:,i:i+1,2,k:k+1)   ,&
                    IsTrueCoarse2    = true_cell(i,-1,k,iBlock)        ,&
                    IsTrueCoarse1    = true_cell(i, 0,k,iBlock)        ,&
                    IsTrueFine1  =all(true_cell(i:i+1, 1,k:k+1,iBlock)),&
                    IsTrueFine2_II      =true_cell(i:i+1, 2,k:k+1,iBlock))
            end if
         end do; end do
      case(4)
         do k=1,nK,2; do i=1,nI,2
            if(  all(true_cell(i:i+1,nJ-1:nJ+2,k:k+1,iBlock)) .and. &
                 all(true_cell(i-2:i+3,nJ+1,k-2:k+3,iBlock)) ) then
               call accurate_reschange3d(&
                    Coarse2_V    =    Primitive_VG(:,i,nJ+2,k)         ,&
                    Coarse1_VII  =    Primitive_VG(:,i-2:i+3,nJ+1,k-2:k+3),&
                    Fine1_VII    =    Primitive_VG(:,i:i+1, nJ,k:k+1)  ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1, nJ-1,k:k+1),&
                    CoarseToFineF_VII=RightState_VY(:,i:i+1,nJ+1,k:k+1),&
                    FineToCoarseF_VII=LeftState_VY(:,i:i+1,nJ+1,k:k+1) ,&
                    FineF_VII        =RightState_VY(:,i:i+1,nJ,k:k+1))
            else
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,i,nJ+2,k)         ,&
                    Coarse1_V    =    Primitive_VG(:,i,nJ+1,k)         ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1, nJ,k:k+1)  ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1, nJ-1,k:k+1),&
                    CoarseToFineF_VII=RightState_VY(:,i:i+1,nJ+1,k:k+1),&
                    FineToCoarseF_VII=LeftState_VY(:,i:i+1,nJ+1,k:k+1) ,&
                    FineF_VII        =RightState_VY(:,i:i+1,nJ,k:k+1)  ,&
                    IsTrueCoarse2    = true_cell(i,nJ+2,k,iBlock)      ,&
                    IsTrueCoarse1    = true_cell(i,nJ+1,k,iBlock)      ,&
                    IsTrueFine1  =all(true_cell(i:i+1,nJ,k:k+1,iBlock)),&
                    IsTrueFine2_II      =true_cell(i:i+1,nJ-1,k:k+1,iBlock))
            end if
         end do; end do
      case(5)
         do j=1,nJ,2; do i=1,nI,2
            if(  all(true_cell(i:i+1,j:j+1,-1:2,iBlock)) .and. &
                 all(true_cell(i-2:i+3,j-2:j+3,0,iBlock)) ) then
               call accurate_reschange3d(&
                    Coarse2_V    =    Primitive_VG(:,i,j,-1)           ,&
                    Coarse1_VII  =    Primitive_VG(:,i-2:i+3,j-2:j+3,0),&
                    Fine1_VII    =    Primitive_VG(:,i:i+1,j:j+1, 1)   ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1,j:j+1, 2)   ,&
                    CoarseToFineF_VII= LeftState_VZ(:,i:i+1,j:j+1,1)   ,&
                    FineToCoarseF_VII=RightState_VZ(:,i:i+1,j:j+1,1)   ,&
                    FineF_VII        = LeftState_VZ(:,i:i+1,j:j+1,2))
            else
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,i,j,-1)           ,&
                    Coarse1_V    =    Primitive_VG(:,i,j, 0)           ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1,j:j+1, 1)   ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1,j:j+1, 2)   ,&
                    CoarseToFineF_VII= LeftState_VZ(:,i:i+1,j:j+1,1)   ,&
                    FineToCoarseF_VII=RightState_VZ(:,i:i+1,j:j+1,1)   ,&
                    FineF_VII        = LeftState_VZ(:,i:i+1,j:j+1,2)   ,&
                    IsTrueCoarse2    = true_cell(i,j,-1,iBlock)        ,&
                    IsTrueCoarse1    = true_cell(i,j, 0,iBlock)        ,&
                    IsTrueFine1 =all(true_cell(i:i+1,j:j+1, 1,iBlock)),&
                    IsTrueFine2_II      =true_cell(i:i+1,j:j+1, 2,iBlock))
            end if
         end do; end do
      case(6)
         do j=1,nJ,2; do i=1,nI,2
            if(  all(true_cell(i:i+1,j:j+1,nK-1:nK+2,iBlock)) .and. &
                 all(true_cell(i-2:i+3,j-2:j+3,nK+1,iBlock)) ) then
               call accurate_reschange3d(&
                    Coarse2_V    =    Primitive_VG(:,i,j,nK+2)         ,&
                    Coarse1_VII  =    Primitive_VG(:,i-2:i+3,j-2:j+3,nK+1),&
                    Fine1_VII    =    Primitive_VG(:,i:i+1,j:j+1, nK)  ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1,j:j+1, nK-1),&
                    CoarseToFineF_VII=RightState_VZ(:,i:i+1,j:j+1,nK+1),&
                    FineToCoarseF_VII=LeftState_VZ(:,i:i+1,j:j+1,nK+1) ,&
                    FineF_VII        =RightState_VZ(:,i:i+1,j:j+1,nK))
            else
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,i,j,nK+2)         ,&
                    Coarse1_V    =    Primitive_VG(:,i,j,nK+1)         ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1,j:j+1, nK)  ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1,j:j+1, nK-1),&
                    CoarseToFineF_VII=RightState_VZ(:,i:i+1,j:j+1,nK+1),&
                    FineToCoarseF_VII=LeftState_VZ(:,i:i+1,j:j+1,nK+1) ,&
                    FineF_VII        =RightState_VZ(:,i:i+1,j:j+1,nK)  ,&
                    IsTrueCoarse2    =true_cell(i,j,nK+2,iBlock)       ,&
                    IsTrueCoarse1    =true_cell(i,j,nK+1,iBlock)       ,&
                    IsTrueFine1  =all(true_cell(i:i+1,j:j+1,nK,iBlock)),&
                    IsTrueFine2_II  =true_cell(i:i+1,j:j+1,nK-1,iBlock))
            end if
         end do; end do
      end select
    end subroutine get_face_accurate3d
    !==========================================================================
    subroutine get_face_accurate1d(iSideIn)
      integer, intent(in):: iSideIn

      !------------------------------------------------------------------------
      select case(iSideIn)
      case(1)
         call accurate_reschange1d(&
              Coarse2_V       = Primitive_VG(:,-1,1,1)     ,&
              Coarse1_V       = Primitive_VG(:, 0,1,1)     ,&
              Fine1_V         = Primitive_VG(:, 1,1,1)     ,&
              Fine2_V         = Primitive_VG(:, 2,1,1)     ,&
              CoarseToFineF_V = LeftState_VX(:, 1,1,1)     ,&
              FineToCoarseF_V =RightState_VX(:, 1,1,1)     ,&
              FineF_V         = LeftState_VX(:, 2,1,1))
      case(2)
         call accurate_reschange1d(&
              Coarse2_V       = Primitive_VG(:,nI+2,1,1)   ,&
              Coarse1_V       = Primitive_VG(:,nI+1,1,1)   ,&
              Fine1_V         = Primitive_VG(:,nI  ,1,1)   ,&
              Fine2_V         = Primitive_VG(:,nI-1,1,1)   ,&
              CoarseToFineF_V =RightState_VX(:,nI+1,1,1)   ,&
              FineToCoarseF_V = LeftState_VX(:,nI+1,1,1)   ,&
              FineF_V         =RightState_VX(:,nI  ,1,1))
      end select

    end subroutine get_face_accurate1d
    !==========================================================================
    subroutine get_face_accurate2d(iSideIn)
      integer, intent(in):: iSideIn

      !------------------------------------------------------------------------
      select case(iSideIn)
      case(1)
         do j=1,nJ,2
            call accurate_reschange2d(&
                 Coarse2_V       = Primitive_VG(:,-1,j,1)         ,&
                 Coarse1_VI      = Primitive_VG(:, 0,j-2:j+3,1)   ,&
                 Fine1_VI        = Primitive_VG(:, 1,j:j+1,1)     ,&
                 Fine2_VI        = Primitive_VG(:, 2,j:j+1,1)     ,&
                 CoarseToFineF_VI= LeftState_VX(:, 1,j:j+1,1)     ,&
                 FineToCoarseF_VI=RightState_VX(:, 1,j:j+1,1)     ,&
                 FineF_VI        = LeftState_VX(:, 2,j:j+1,1))
         end do
      case(2)
         do j=1,nJ,2
            call accurate_reschange2d(&
                 Coarse2_V       = Primitive_VG(:,nI+2,j,1)       ,&
                 Coarse1_VI      = Primitive_VG(:,nI+1,j-2:j+3,1) ,&
                 Fine1_VI        = Primitive_VG(:,nI  ,j:j+1,1)   ,&
                 Fine2_VI        = Primitive_VG(:,nI-1,j:j+1,1)   ,&
                 CoarseToFineF_VI=RightState_VX(:,nI+1,j:j+1,1)   ,&
                 FineToCoarseF_VI= LeftState_VX(:,nI+1,j:j+1,1)   ,&
                 FineF_VI        =RightState_VX(:,nI  ,j:j+1,1))
         end do
      case(3)
         do i=1,nI,2
            call accurate_reschange2d(&
                 Coarse2_V       = Primitive_VG(:,i,-1,1)         ,&
                 Coarse1_VI      = Primitive_VG(:,i-2:i+3,0,1)    ,&
                 Fine1_VI        = Primitive_VG(:,i:i+1,1,1)      ,&
                 Fine2_VI        = Primitive_VG(:,i:i+1,2,1)      ,&
                 CoarseToFineF_VI= LeftState_VY(:,i:i+1,1,1)      ,&
                 FineToCoarseF_VI=RightState_VY(:,i:i+1,1,1)      ,&
                 FineF_VI        = LeftState_VY(:,i:i+1,2,1))
         end do
      case(4)
         do i=1,nI,2
            call accurate_reschange2d(&
                 Coarse2_V       = Primitive_VG(:,i,nJ+2,1)       ,&
                 Coarse1_VI      = Primitive_VG(:,i-2:i+3,nJ+1,1) ,&
                 Fine1_VI        = Primitive_VG(:,i:i+1,nJ,1)     ,&
                 Fine2_VI        = Primitive_VG(:,i:i+1,nJ-1,1)   ,&
                 CoarseToFineF_VI=RightState_VY(:,i:i+1,nJ+1,1)   ,&
                 FineToCoarseF_VI= LeftState_VY(:,i:i+1,nJ+1,1)   ,&
                 FineF_VI        =RightState_VY(:,i:i+1,nJ,1))
         end do
      end select
    end subroutine get_face_accurate2d
    !==========================================================================
    subroutine get_face_tvd(iSideIn)
      integer,intent(in)::iSideIn

      !------------------------------------------------------------------------
      select case(iSideIn)
      case(1)
         do k=1,nK,2; do j=1,nJ,2
            if(.not.all(true_cell(-1:2,j:j+1,k:k+1,iBlock)))then
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,-1,j,k)           ,&
                    Coarse1_V    =    Primitive_VG(:, 0,j,k)           ,&
                    Fine1_VII    =    Primitive_VG(:, 1,j:j+1,k:k+1)   ,&
                    Fine2_VII    =    Primitive_VG(:, 2,j:j+1,k:k+1)   ,&
                    CoarseToFineF_VII= LeftState_VX(:,1,j:j+1,k:k+1)   ,&
                    FineToCoarseF_VII=RightState_VX(:,1,j:j+1,k:k+1)   ,&
                    FineF_VII        = LeftState_VX(:,2,j:j+1,k:k+1)   ,&
                    IsTrueCoarse2    = true_cell(-1,j,k,iBlock)        ,&
                    IsTrueCoarse1    = true_cell( 0,j,k,iBlock)        ,&
                    IsTrueFine1  =all(true_cell( 1,j:j+1,k:k+1,iBlock)),&
                    IsTrueFine2_II      =true_cell( 2,j:j+1,k:k+1,iBlock))
            else
               call tvd_reschange(&
                    Coarse2_V    =    Primitive_VG(:,-1,j,k)           ,&
                    Coarse1_V    =    Primitive_VG(:, 0,j,k)           ,&
                    Fine1_VII    =    Primitive_VG(:, 1,j:j+1,k:k+1)   ,&
                    Fine2_VII    =    Primitive_VG(:, 2,j:j+1,k:k+1)   ,&
                    CoarseToFineF_VII= LeftState_VX(:,1,j:j+1,k:k+1)   ,&
                    FineToCoarseF_VII=RightState_VX(:,1,j:j+1,k:k+1)   ,&
                    FineF_VII        = LeftState_VX(:,2,j:j+1,k:k+1))
            end if
         end do; end do
      case(2)
         do k=1,nK,2; do j=1,nJ,2
            if(.not.all(true_cell(nI-1:nI+2,j:j+1,k:k+1,iBlock)))then
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,nI+2,j,k)         ,&
                    Coarse1_V    =    Primitive_VG(:, nI+1,j,k)        ,&
                    Fine1_VII    =    Primitive_VG(:, nI,j:j+1,k:k+1)  ,&
                    Fine2_VII    =    Primitive_VG(:, nI-1,j:j+1,k:k+1),&
                    CoarseToFineF_VII=RightState_VX(:,nI+1,j:j+1,k:k+1),&
                    FineToCoarseF_VII=LeftState_VX(:,nI+1,j:j+1,k:k+1) ,&
                    FineF_VII        =RightState_VX(:,nI,j:j+1,k:k+1)  ,&
                    IsTrueCoarse2    = true_cell(nI+2,j,k,iBlock)      ,&
                    IsTrueCoarse1    = true_cell(nI+1,j,k,iBlock)      ,&
                    IsTrueFine1  =all(true_cell(nI,j:j+1,k:k+1,iBlock)),&
                    IsTrueFine2_II      =true_cell(nI-1,j:j+1,k:k+1,iBlock))
            else
               call tvd_reschange(&
                    Coarse2_V    =    Primitive_VG(:,nI+2,j,k)         ,&
                    Coarse1_V    =    Primitive_VG(:, nI+1,j,k)        ,&
                    Fine1_VII    =    Primitive_VG(:, nI,j:j+1,k:k+1)  ,&
                    Fine2_VII    =    Primitive_VG(:, nI-1,j:j+1,k:k+1),&
                    CoarseToFineF_VII=RightState_VX(:,nI+1,j:j+1,k:k+1),&
                    FineToCoarseF_VII=LeftState_VX(:,nI+1,j:j+1,k:k+1) ,&
                    FineF_VII       =RightState_VX(:,nI,j:j+1,k:k+1))
            end if
         end do; end do
      case(3)
         do k=1,nK,2; do i=1,nI,2
            if(.not.all(true_cell(i:i+1,-1:2,k:k+1,iBlock)))then
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,i,-1,k)           ,&
                    Coarse1_V    =    Primitive_VG(:,i, 0,k)           ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1, 1,k:k+1)   ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1, 2,k:k+1)   ,&
                    CoarseToFineF_VII= LeftState_VY(:,i:i+1,1,k:k+1)   ,&
                    FineToCoarseF_VII=RightState_VY(:,i:i+1,1,k:k+1)   ,&
                    FineF_VII        = LeftState_VY(:,i:i+1,2,k:k+1)   ,&
                    IsTrueCoarse2    = true_cell(i,-1,k,iBlock)        ,&
                    IsTrueCoarse1    = true_cell(i, 0,k,iBlock)        ,&
                    IsTrueFine1  =all(true_cell(i:i+1, 1,k:k+1,iBlock)),&
                    IsTrueFine2_II      =true_cell(i:i+1, 2,k:k+1,iBlock))
            else
               call tvd_reschange(&
                    Coarse2_V    =    Primitive_VG(:,i,-1,k)           ,&
                    Coarse1_V    =    Primitive_VG(:,i, 0,k)           ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1, 1,k:k+1)   ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1, 2,k:k+1)   ,&
                    CoarseToFineF_VII= LeftState_VY(:,i:i+1,1,k:k+1)   ,&
                    FineToCoarseF_VII=RightState_VY(:,i:i+1,1,k:k+1)   ,&
                    FineF_VII       = LeftState_VY(:,i:i+1,2,k:k+1))
            end if
         end do; end do
      case(4)
         do k=1,nK,2; do i=1,nI,2
            if(.not.all(true_cell(i:i+1,nJ-1:nJ+2,k:k+1,iBlock)))then
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,i,nJ+2,k)         ,&
                    Coarse1_V    =    Primitive_VG(:,i,nJ+1,k)         ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1, nJ,k:k+1)  ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1, nJ-1,k:k+1),&
                    CoarseToFineF_VII=RightState_VY(:,i:i+1,nJ+1,k:k+1),&
                    FineToCoarseF_VII=LeftState_VY(:,i:i+1,nJ+1,k:k+1) ,&
                    FineF_VII        =RightState_VY(:,i:i+1,nJ,k:k+1)  ,&
                    IsTrueCoarse2    = true_cell(i,nJ+2,k,iBlock)      ,&
                    IsTrueCoarse1    = true_cell(i,nJ+1,k,iBlock)      ,&
                    IsTrueFine1  =all(true_cell(i:i+1,nJ,k:k+1,iBlock)),&
                    IsTrueFine2_II      =true_cell(i:i+1,nJ-1,k:k+1,iBlock))
            else
               call tvd_reschange(&
                    Coarse2_V    =    Primitive_VG(:,i,nJ+2,k)         ,&
                    Coarse1_V    =    Primitive_VG(:,i,nJ+1,k)         ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1, nJ,k:k+1)  ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1, nJ-1,k:k+1),&
                    CoarseToFineF_VII=RightState_VY(:,i:i+1,nJ+1,k:k+1),&
                    FineToCoarseF_VII=LeftState_VY(:,i:i+1,nJ+1,k:k+1) ,&
                    FineF_VII        =RightState_VY(:,i:i+1,nJ,k:k+1))
            end if
         end do; end do
      case(5)
         do j=1,nJ,2; do i=1,nI,2
            if(.not.all(true_cell(i:i+1,j:j+1,-1:2,iBlock)))then
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,i,j,-1)           ,&
                    Coarse1_V    =    Primitive_VG(:,i,j, 0)           ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1,j:j+1, 1)   ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1,j:j+1, 2)   ,&
                    CoarseToFineF_VII= LeftState_VZ(:,i:i+1,j:j+1,1)   ,&
                    FineToCoarseF_VII=RightState_VZ(:,i:i+1,j:j+1,1)   ,&
                    FineF_VII        = LeftState_VZ(:,i:i+1,j:j+1,2)   ,&
                    IsTrueCoarse2    = true_cell(i,j,-1,iBlock)        ,&
                    IsTrueCoarse1    = true_cell(i,j, 0,iBlock)        ,&
                    IsTrueFine1 =all(true_cell(i:i+1,j:j+1, 1,iBlock)),&
                    IsTrueFine2_II      =true_cell(i:i+1,j:j+1, 2,iBlock))
            else
               call tvd_reschange(&
                    Coarse2_V    =    Primitive_VG(:,i,j,-1)           ,&
                    Coarse1_V    =    Primitive_VG(:,i,j, 0)           ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1,j:j+1, 1)   ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1,j:j+1, 2)   ,&
                    CoarseToFineF_VII= LeftState_VZ(:,i:i+1,j:j+1,1)   ,&
                    FineToCoarseF_VII=RightState_VZ(:,i:i+1,j:j+1,1)   ,&
                    FineF_VII        = LeftState_VZ(:,i:i+1,j:j+1,2))
            end if
         end do; end do
      case(6)
         do j=1,nJ,2; do i=1,nI,2
            if(.not.all(true_cell(i:i+1,j:j+1,nK-1:nK+2,iBlock)))then
               call tvd_reschange_body(&
                    Coarse2_V    =    Primitive_VG(:,i,j,nK+2)         ,&
                    Coarse1_V    =    Primitive_VG(:,i,j,nK+1)         ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1,j:j+1, nK)  ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1,j:j+1, nK-1),&
                    CoarseToFineF_VII=RightState_VZ(:,i:i+1,j:j+1,nK+1),&
                    FineToCoarseF_VII=LeftState_VZ(:,i:i+1,j:j+1,nK+1) ,&
                    FineF_VII        =RightState_VZ(:,i:i+1,j:j+1,nK)  ,&
                    IsTrueCoarse2    =true_cell(i,j,nK+2,iBlock)       ,&
                    IsTrueCoarse1    =true_cell(i,j,nK+1,iBlock)       ,&
                    IsTrueFine1  =all(true_cell(i:i+1,j:j+1,nK,iBlock)),&
                    IsTrueFine2_II  =true_cell(i:i+1,j:j+1,nK-1,iBlock))
            else
               call tvd_reschange(&
                    Coarse2_V    =    Primitive_VG(:,i,j,nK+2)         ,&
                    Coarse1_V    =    Primitive_VG(:,i,j,nK+1)         ,&
                    Fine1_VII    =    Primitive_VG(:,i:i+1,j:j+1, nK)  ,&
                    Fine2_VII    =    Primitive_VG(:,i:i+1,j:j+1, nK-1),&
                    CoarseToFineF_VII=RightState_VZ(:,i:i+1,j:j+1,nK+1),&
                    FineToCoarseF_VII=LeftState_VZ(:,i:i+1,j:j+1,nK+1) ,&
                    FineF_VII        =RightState_VZ(:,i:i+1,j:j+1,nK))
            end if
         end do; end do
      end select
    end subroutine get_face_tvd
    !==========================================================================
    subroutine get_facex_second(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      integer::i1, iMinSharp, iMaxSharp

      !------------------------------------------------------------------------
      iMinSharp = iMin
      iMaxSharp = iMax
      if(BetaLimiter > BetaLimiterResChange)then
         if(neiLeast(iBlock) /= 0) iMinSharp = &
              max(iMin, min(iMax + 1,      1 + nFaceLimiterResChange))
         if(neiLwest(iBlock) /= 0) iMaxSharp = &
              min(iMax, max(iMin - 1, nI + 1 - nFaceLimiterResChange))
      endif

      do k=kMin, kMax; do j=jMin, jMax;
         Primitive_VI(:,iMin-2:iMax+1)=Primitive_VG(:,iMin-2:iMax+1,j,k)
         if(UseTrueCell)then
            IsTrueCell_I(iMin-2:iMax+1) = true_cell(iMin-2:iMax+1,j,k,iBlock)
            if(iMinSharp <= iMaxSharp) &
                 call limiter_body(iMinSharp, iMaxSharp, BetaLimiter)
            if(iMin < iMinSharp) &
                 call limiter_body(iMin, iMinSharp-1, BetaLimiterResChange)
            if(iMax > iMaxSharp) &
                 call limiter_body(iMaxSharp+1, iMax, BetaLimiterResChange)
         else
            if(iMinSharp <= iMaxSharp) &
                 call limiter(iMinSharp, iMaxSharp, BetaLimiter)
            if(iMin < iMinSharp) &
                 call limiter(iMin, iMinSharp-1, BetaLimiterResChange)
            if(iMax > iMaxSharp) &
                 call limiter(iMaxSharp+1, iMax, BetaLimiterResChange)
         end if
         do i=iMin,iMax
            i1=i-1
            LeftState_VX(:,i,j,k)  = Primitive_VI(:,i1) + dVarLimL_VI(:,i1)
            RightState_VX(:,i,j,k) = Primitive_VI(:,i ) - dVarLimR_VI(:,i )
         end do
      end do; end do

      if(DoLimitMomentum) call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facex_second
    !==========================================================================
    subroutine get_facey_second(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      integer::j1, jMinSharp, jMaxSharp

      !------------------------------------------------------------------------
      jMinSharp = jMin
      jMaxSharp = jMax
      if(BetaLimiter > BetaLimiterResChange)then
         if(neiLsouth(iBlock) /= 0) jMinSharp = &
              max(jMin, min(jMax + 1,      1 + nFaceLimiterResChange))
         if(neiLnorth(iBlock) /= 0) jMaxSharp = &
              min(jMax, max(jMin - 1, nJ + 1 - nFaceLimiterResChange))
      endif

      do k=kMin, kMax; do i=iMin,iMax
         Primitive_VI(:,jMin-2:jMax+1)=Primitive_VG(:,i,jMin-2:jMax+1,k)
         if(UseTrueCell)then
            IsTrueCell_I(jMin-2:jMax+1) = true_cell(i,jMin-2:jMax+1,k,iBlock)
            if(jMinSharp <= jMaxSharp) &
                 call limiter_body(jMinSharp, jMaxSharp, BetaLimiter)
            if(jMin < jMinSharp) &
                 call limiter_body(jMin, jMinSharp-1, BetaLimiterResChange)
            if(jMax > jMaxSharp) &
                 call limiter_body(jMaxSharp+1, jMax, BetaLimiterResChange)
         else
            if(jMinSharp <= jMaxSharp) &
                 call limiter(jMinSharp, jMaxSharp, BetaLimiter)
            if(jMin < jMinSharp) &
                 call limiter(jMin, jMinSharp-1, BetaLimiterResChange)
            if(jMax > jMaxSharp) &
                 call limiter(jMaxSharp+1, jMax, BetaLimiterResChange)
         end if
         do j=jMin, jMax
            j1=j-1
            LeftState_VY(:,i,j,k)  = Primitive_VI(:,j1) + dVarLimL_VI(:,j1)
            RightState_VY(:,i,j,k) = Primitive_VI(:,j ) - dVarLimR_VI(:,j )
         end do
      end do; end do

      if(DoLimitMomentum) call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facey_second
    !==========================================================================
    subroutine get_facez_second(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer::k1, kMinSharp, kMaxSharp

      !------------------------------------------------------------------------
      kMinSharp = kMin
      kMaxSharp = kMax
      if(BetaLimiter > BetaLimiterResChange)then
         if(neiLbot(iBlock) /= 0) kMinSharp = &
              max(kMin, min(kMax + 1,      1 + nFaceLimiterResChange))
         if(neiLtop(iBlock) /= 0) kMaxSharp = &
              min(kMax, max(kMin - 1, nK + 1 - nFaceLimiterResChange))
      endif
      do j=jMin,jMax; do i=iMin,iMax;
         Primitive_VI(:,kMin-2:kMax+1)=Primitive_VG(:,i,j,kMin-2:kMax+1)
         if(UseTrueCell)then
            IsTrueCell_I(kMin-2:kMax+1) = true_cell(i,j,kMin-2:kMax+1,iBlock)
            if(kMinSharp <= kMaxSharp) &
                 call limiter_body(kMinSharp, kMaxSharp, BetaLimiter)
            if(kMin < kMinSharp) &
                 call limiter_body(kMin, kMinSharp-1, BetaLimiterResChange)
            if(kMax > kMaxSharp) &
                 call limiter_body(kMaxSharp+1, kMax, BetaLimiterResChange)
         else
            if(kMinSharp <= kMaxSharp) &
                 call limiter(kMinSharp, kMaxSharp, BetaLimiter)
            if(kMin < kMinSharp) &
                 call limiter(kMin, kMinSharp-1, BetaLimiterResChange)
            if(kMax > kMaxSharp) &
                 call limiter(kMaxSharp+1, kMax, BetaLimiterResChange)
         end if
         do k=kMin,kMax
            k1=k-1
            LeftState_VZ(:,i,j,k)  = Primitive_VI(:,k1) + dVarLimL_VI(:,k1)
            RightState_VZ(:,i,j,k) = Primitive_VI(:,k ) - dVarLimR_VI(:,k )
         end do
      end do; end do

      if(DoLimitMomentum) call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facez_second
    !==========================================================================

    subroutine flatten(Prim_VG)

      use ModMultiFluid, ONLY: iFluid, iRho, iUx, iUy, iUz, iP, select_fluid

      real, intent(in):: Prim_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

      real:: InvRatioRange

      real:: pL, pR, Dp, Ratio, Coef1, Coef

      real:: FlatCoef_I(-1:MaxIJK+2)
      real, allocatable:: FlatCoef_G(:,:,:)

      integer:: i, j, k
      !------------------------------------------------------------------------
      ! call timing_start('flatten')

      if(.not.allocated(FlatCoef_G)) &
           allocate(FlatCoef_G(0:nI+1,1-jDim_:nJ+jDim_,1-kDim_:nK+kDim_))

      InvRatioRange = 1.0/(FlatRatioMax - FlatRatioMin)

      FLUIDLOOP: do iFluid = 1, nFluid

         call select_fluid

         do k = 1-kDim_,nK+kDim_; do j = 1-jDim_,nJ+jDim_
            do i = -1, nI+2

               ! Coef = 1 preserves the high order face value
               FlatCoef_I(i) = 1.0

               if(UseDuFlat)then
                  ! Check if there is compression
                  ! Note: Balsara suggests to look at rarefactions too...)
                  if(Prim_VG(iUx,i-1,j,k) - Prim_VG(iUx,i+1,j,k) <= 0)CYCLE
               end if

               pL = Prim_VG(iP,i-1,j,k)
               pR = Prim_VG(iP,i+1,j,k)
               Dp = abs(pR - pL)

               ! Check the shock strength. Nothing to do for weak shock
               if(Dp < FlatDelta*min(pL, pR)) CYCLE

               ! Calculate the shock width parameter
               pL = Prim_VG(iP,i-2,j,k)
               pR = Prim_VG(iP,i+2,j,k)
               Ratio = Dp / max(1e-30, abs(pR - pL))

               if(Ratio > FlatRatioMax)then
                  FlatCoef_I(i) = 0
               elseif(Ratio > FlatRatioMin)then
                  FlatCoef_I(i) = InvRatioRange*(FlatRatioMax - Ratio)
               end if
            end do

            do i = 0, nI+1
               FlatCoef_G(i,j,k) = minval(FlatCoef_I(i-1:i+1))
            end do
         end do; end do

         if(nDim > 1)then
            do k = 1-kDim_,nK+kDim_; do i = 0, nI+1
               do j = -1, nJ+2
                  FlatCoef_I(j) = 1.0

                  if(UseDuFlat)then
                     if(Prim_VG(iUy,i,j-1,k) - Prim_VG(iUy,i,j+1,k) <= 0)CYCLE
                  end if

                  pL = Prim_VG(iP,i,j-1,k)
                  pR = Prim_VG(iP,i,j+1,k)
                  Dp = abs(pR - pL)

                  if(Dp < FlatDelta*min(pL, pR)) CYCLE

                  pL = Prim_VG(iP,i,j-2,k)
                  pR = Prim_VG(iP,i,j+2,k)
                  Ratio = Dp / max(1e-30, abs(pR - pL))

                  if(Ratio > FlatRatioMax)then
                     FlatCoef_I(j) = 0
                  elseif(Ratio > FlatRatioMin)then
                     FlatCoef_I(j) = InvRatioRange*(FlatRatioMax - Ratio)
                  end if
               end do

               do j = 0, nJ+1
                  FlatCoef_G(i,j,k) = &
                       min(FlatCoef_G(i,j,k), minval(FlatCoef_I(j-1:j+1)))
               end do
            end do; end do
         end if

         if(nDim > 2)then
            do j = 0, nJ+1; do i = 0, nI+1
               do k = -1, nK+2
                  FlatCoef_I(k) = 1.0
                  if(UseDuFlat)then
                     if(Prim_VG(iUz,i,j,k-1) - Prim_VG(iUz,i,j,k+1) <= 0)CYCLE
                  end if

                  pL = Prim_VG(iP,i,j,k-1)
                  pR = Prim_VG(iP,i,j,k+1)
                  Dp = abs(pR - pL)

                  if(Dp < FlatDelta*min(pL, pR)) CYCLE

                  pL = Prim_VG(iP,i,j,k-2)
                  pR = Prim_VG(iP,i,j,k+2)
                  Ratio = Dp / max(1e-30, abs(pR - pL))

                  if(Ratio > FlatRatioMax)then
                     FlatCoef_I(k) = 0
                  elseif(Ratio > FlatRatioMin)then
                     FlatCoef_I(k) = InvRatioRange*(FlatRatioMax - Ratio)
                  end if
               end do

               do k = 0, nK+1
                  FlatCoef_G(i,j,k) = &
                       min(FlatCoef_G(i,j,k), minval(FlatCoef_I(k-1:k+1)))
               end do
            end do; end do
         end if

         do k = kMinFace, kMaxFace; do j = jMinFace, jMaxFace; do i = 0, nI+1
            Coef = FlatCoef_G(i,j,k)

            ! Coef is the final flattening parameter in eq. 34a,b
            if(Coef > 1 - 1e-12) CYCLE

            Coef1 = 1.0 - Coef
            if(i<=nI) LeftState_VX(iRho:iP,i+1,j,k) = &
                 Coef*LeftState_VX(iRho:iP,i+1,j,k)   &
                 + Coef1*Prim_VG(iRho:iP,i,j,k)
            if(i> 0 ) RightState_VX(iRho:iP,i,j,k)  = &
                 Coef*RightState_VX(iRho:iP,i,j,k)    &
                 + Coef1*Prim_VG(iRho:iP,i,j,k)
         end do; end do; end do

         if(nDim == 1) CYCLE FLUIDLOOP

         do k = kMinFace, kMaxFace; do j = 0, nJ+1; do i = iMinFace,iMaxFace
            Coef = FlatCoef_G(i,j,k)
            if(Coef > 1 - 1e-12) CYCLE
            Coef1 = 1.0 - Coef
            if(j<=nJ) LeftState_VY(iRho:iP,i,j+1,k) = &
                 Coef*LeftState_VY(iRho:iP,i,j+1,k)   &
                 + Coef1*Prim_VG(iRho:iP,i,j,k)
            if(j> 0 ) RightState_VY(iRho:iP,i,j,k)  = &
                 Coef*RightState_VY(iRho:iP,i,j,k)    &
                 + Coef1*Prim_VG(iRho:iP,i,j,k)
         end do; end do; end do

         if(nDim == 2) CYCLE FLUIDLOOP

         do k = 0, nK+1; do j = jMinFace, jMaxFace; do i = iMinFace,iMaxFace
            Coef = FlatCoef_G(i,j,k)
            if(Coef > 1 - 1e-12) CYCLE
            Coef1 = 1.0 - Coef
            if(k<=nK) LeftState_VZ(iRho:iP,i,j,k+1) = &
                 Coef*LeftState_VZ(iRho:iP,i,j,k+1)   &
                 + Coef1*Prim_VG(iRho:iP,i,j,k)
            if(k> 0 ) RightState_VZ(iRho:iP,i,j,k)  = &
                 Coef*RightState_VZ(iRho:iP,i,j,k)    &
                 + Coef1*Prim_VG(iRho:iP,i,j,k)
         end do; end do; end do

      end do FLUIDLOOP

      ! call timing_stop('flatten')

    end subroutine flatten
    !==========================================================================

  end subroutine calc_face_value
  !============================================================================

  subroutine set_low_order_face

    use ModMain,  ONLY: MaxBlock, iMinFace, iMaxFace, jMinFace, jMaxFace, &
         kMinFace, kMaxFace, nIFace, nJFace, nKFace, nOrder, nBlock
    use ModGeometry, ONLY : true_cell
    use ModB0
    use ModPhysics, ONLY: Gamma_I
    use ModAdvance, ONLY: LowOrderCrit_XB, LowOrderCrit_YB, LowOrderCrit_ZB, &
         State_VGB
    use BATL_lib, ONLY: block_inside_regions, Unused_B, MaxDim

    ! Set which faces should use low (up to second) order scheme
    ! Set logicals for the current block

    integer:: iBlock
    logical:: UseLowOrderOnly

    integer:: i, j, k
    real:: LowOrderCrit_X(nIFace,nJ,nK)
    real:: LowOrderCrit_Y(nI,nJFace,nK)
    real:: LowOrderCrit_Z(nI,nJ,nKFace)
    integer:: cLowOrder = 1
    real:: cSmall = 1e-6

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'set_low_order_face'
    !------------------------------------------------------------------------

    call test_start(NameSub, DoTest)

    if(nOrder<=2) then
       IsLowOrderOnly_B = .true.
    else
       IsLowOrderOnly_B = .false.
       UseLowOrder = UseTrueCell .or. UseLowOrderRegion .or. UseAdaptiveLowOrder
       if(UseLowOrder) then
          do iBlock = 1, nBlock
             if (Unused_B(iBlock)) CYCLE
             UseLowOrderOnly = .false.

             if(.not. allocated(LowOrderCrit_XB)) then
                allocate(LowOrderCrit_XB( &
                     nI+1,jMinFace:jMaxFace,kMinFace:kMaxFace,MaxBlock))
                allocate(LowOrderCrit_YB( &
                     iMinFace:iMaxFace,nJ+1,kMinFace:kMaxFace,MaxBlock))
                allocate(LowOrderCrit_ZB( &
                     iMinFace:iMaxFace,jMinFace:jMaxFace,nK+1,MaxBlock))
             endif

             LowOrderCrit_XB(:,:,:,iBlock) = 0
             LowOrderCrit_YB(:,:,:,iBlock) = 0
             LowOrderCrit_ZB(:,:,:,iBlock) = 0         

             ! Get the LowOrderCrit based on physical criteria. The vale of 
             ! 'LowOrderCrit_*B' is the ratio of the low order face value. 
             if(UseAdaptiveLowOrder) call set_physics_based_low_order_face

             if(allocated(iRegionLowOrder_I)) then

                call block_inside_regions(iRegionLowOrder_I, iBlock, &
                     size(LowOrderCrit_X), 'xface',Value_I=LowOrderCrit_X)
                LowOrderCrit_XB(:,:,:,iBlock) = &
                     max(LowOrderCrit_X,LowOrderCrit_XB(:,:,:,iBlock))
                UseLowOrderOnly = all(LowOrderCrit_XB(:,:,:,iBlock) &
                     >= cLowOrder-cSmall)

                if(nDim > 1)then
                   call block_inside_regions(iRegionLowOrder_I, iBlock, &
                        size(LowOrderCrit_Y), 'yface', Value_I=LowOrderCrit_Y)
                   ! Use low order if the face is inside the low order region OR
                   ! satisfies the physical condition.
                   LowOrderCrit_YB(:,:,:,iBlock) = &
                        max(LowOrderCrit_Y,LowOrderCrit_YB(:,:,:,iBlock))
                   UseLowOrderOnly = UseLowOrderOnly .and. &
                        all(LowOrderCrit_YB(:,:,:,iBlock) >= cLowOrder-cSmall)
                end if

                if(nDim > 2)then
                   call block_inside_regions(iRegionLowOrder_I, iBlock, &
                        size(LowOrderCrit_Z), 'zface', Value_I=LowOrderCrit_Z)
                   LowOrderCrit_ZB(:,:,:,iBlock) = &
                        max(LowOrderCrit_Z,LowOrderCrit_ZB(:,:,:,iBlock))
                   UseLowOrderOnly = UseLowOrderOnly .and. &
                        all(LowOrderCrit_ZB(:,:,:,iBlock) >= cLowOrder-cSmall)
                end if

             else if(UseTrueCell)then
                ! The 5th order schemes need 3 cells on both sides of the face
                do k=1, nK; do j=1, nJ; do i = 1, nI+1
                   if(.not.all(true_cell(i-3:i+2,j,k,iBlock))) &
                        LowOrderCrit_XB(i,j,k,iBlock) = cLowOrder
                end do; end do; end do
                if(nDim > 1)then
                   do k=1, nK; do j=1, nJ+1; do i = 1, nI
                      if(.not.all(true_cell(i,j-3:j+2,j,iBlock))) &
                           LowOrderCrit_YB(i,j,k,iBlock) = cLowOrder
                   end do; end do; end do
                end if
                if(nDim > 2)then
                   do k=1, nK+1; do j=1, nJ; do i = 1, nI
                      if(.not.all(true_cell(i,j,k-3:k+2,iBlock))) &
                           LowOrderCrit_ZB(i,j,k,iBlock) = cLowOrder
                   end do; end do; end do
                end if
             end if

             if(UseLowOrderOnly) IsLowOrderOnly_B(iBlock) = .true.

          enddo ! iBlock
       endif ! UseLowOrder
    endif ! nOrder
    call test_stop(NameSub, DoTest)

  contains

    subroutine set_physics_based_low_order_face
      use ModMain, ONLY:MaxBlock
      use ModAdvance, ONLY:Vel_IDGB
      use ModMultiFluid, ONLY: nFluid
      use ModPhysics, ONLY: Gamma_I
      !--------------------------------------------------------------------

      integer :: iFace, jFace, kFace
      real:: State_VI(nVar,-3:2)
      logical:: DoTest
      character(len=*), parameter:: NameSub = 'set_physics_based_low_order_face'
      !--------------------------------------------------------------------------

      call test_start(NameSub, DoTest)

      ! Face along x-direction
      do kFace=kMinFace,kMaxFace; do jFace=jMinFace,jMaxFace; do iFace=1,nIFace
         State_VI = State_VGB(:,iFace-3:iFace+2,jFace,kFace,iBlock)
         if(UseB0) State_VI(Bx_:Bz_,:) = State_VI(Bx_:Bz_,:) + &
              B0_DGB(:,iFace-3:iFace+2,jFace,kFace,iBlock)

         !write(*,*)'iFace = ',iFace
         LowOrderCrit_XB(iFace,jFace,kFace,iBlock) = &
              low_order_face_criteria( State_VI, &
              Vel_IDGB(:,x_,iFace-3:iFace+2,jFace,kFace,iBlock))
      enddo; enddo; enddo

      if(nDim>1) then
         ! Face along y-direction
         do kFace=kMinFace,kMaxFace; do jFace=1,nJFace; do iFace=iMinFace,iMaxFace
            State_VI = State_VGB(:,iFace,jFace-3:jFace+2,kFace,iBlock)
            if(UseB0) State_VI(Bx_:Bz_,:) = State_VI(Bx_:Bz_,:) + &
                 B0_DGB(:,iFace,jFace-3:jFace+2,kFace,iBlock)
            LowOrderCrit_YB(iFace,jFace,kFace,iBlock) = &
                 low_order_face_criteria(State_VI, &
                 Vel_IDGB(:,y_,iFace,jFace-3:jFace+2,kFace,iBlock))
         enddo; enddo; enddo
      endif

      if(nDim>2) then
         ! Face along z-direction
         do kFace=kMinFace,kMaxFace; do jFace=1,nJFace; do iFace=iMinFace,iMaxFace
            State_VI = State_VGB(:,iFace,jFace,kFace-3:kFace+2,iBlock)
            if(UseB0) State_VI(Bx_:Bz_,:) = State_VI(Bx_:Bz_,:) + &
                 B0_DGB(:,iFace,jFace,kFace-3:kFace+2,iBlock)
            LowOrderCrit_ZB(iFace,jFace,kFace,iBlock) = &
                 low_order_face_criteria(State_VI, &
                 Vel_IDGB(:,z_,iFace,jFace,kFace-3:kFace+2,iBlock))
         enddo; enddo; enddo       
      endif

      call test_stop(NameSub, DoTest)

    end subroutine set_physics_based_low_order_face
    !============================================================================
    real function low_order_face_criteria(State_VI, Vel_II)      
      use ModMain, ONLY: UseB0, UseB
      real, intent(in):: State_VI(nVar,-3:2)
      real, intent(in):: Vel_II(nFluid,-3:2)

      integer:: iFluid, iRho, iRhoUx, iRhoUy, iRhoUz, iP
      real:: pTotal_I(-3:2), Sound_I(-3:2)
      real:: crit, pRatio, VelRatio, SoundMin
      !--------------------------------------------------------------------

      pTotal_I = 0
      VelRatio = 0
      do iFluid = 1, nFluid
         iRho   = iRho_I(iFluid)
         iRhoUx = iRhoUx_I(iFluid)
         iRhoUy = iRhoUy_I(iFluid)
         iRhoUz = iRhoUz_I(iFluid)
         iP     = iP_I(iFluid)

         pTotal_I = pTotal_I + State_VI(iP,:)

         Sound_I  = sqrt(State_VI(iP,:)/State_VI(iRho,:)*Gamma_I(iFluid))
         SoundMin = minval(Sound_I)
         VelRatio = max(VelRatio, &
              (maxval(Vel_II(iFluid,:))-minval(Vel_II(iFluid,:)))/SoundMin)
      enddo

      crit = 0
      if(VelRatio > VelCrit) then
         if(UseB) then
            pTotal_I = pTotal_I + 0.5*(State_VI(Bx_,:)**2 + &
                 State_VI(By_,:)**2 + State_VI(Bz_,:)**2)
         endif

         pRatio = maxval(pTotal_I)/minval(pTotal_I)
         ! pRatio >= pCritLow    -> low order
         ! pRatio < pCritHigh    -> high order
         ! otherwise             -> linear combination of low and high order
         ! This function returns the ratio of the low order face value
         if(pRatio >= pCritLow) then
            crit = 1
         elseif(pRatio < pCritHigh) then
            crit = 0
         else
            crit = (pRatio - pCritHigh)/(pCritLow - pCritHigh)
         endif
      endif

      low_order_face_criteria = crit 

    end function low_order_face_criteria
    !============================================================================

  end subroutine set_low_order_face
  !==========================================================================

  subroutine calc_cell_norm_velocity
    use ModMain, ONLY: MaxBlock, nBlock
    use BATL_lib, ONLY: Xyz_DGB, IsCartesian, Unused_B, MaxDim, &
         iDim_, jDim_, kDim_
    use ModAdvance, ONLY: Vel_IDGB, State_VGB
    use ModMultiFluid, ONLY: nFluid

    !--------------------------------------------------------------------
    integer:: iBlock

    integer :: iMin, iMax, jMin, jMax, kMin, kMax, i, j, k, iDim
    real :: ijkDir_DD(MaxDim, MaxDim)

    integer :: iFluid, iRho, iRhoUx, iRhoUy, iRhoUz

    logical :: DoTest
    character(len=*), parameter:: NameSub = 'calc_cell_norm_velocity'
    !--------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(.not. allocated(Vel_IDGB)) &
         allocate(Vel_IDGB(nFluid,MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    iMin = 1 - nG;       iMax = nI + nG
    jMin = 1 - nG*jDim_; jMax = nJ + nG*jDim_
    kMin = 1 - nG*kDim_; kMax = nK + nG*kDim_

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE
       do iFluid = 1, nFluid
          iRho = iRho_I(iFluid)
          iRhoUx = iRhoUx_I(iFluid)
          iRhoUy = iRhoUy_I(iFluid)
          iRhoUz = iRhoUz_I(iFluid)

          if(IsCartesian) then
             Vel_IDGB(iFluid,x_,:,:,:,iBlock) = &
                  State_VGB(iRhoUx,:,:,:,iBlock)/State_VGB(iRho,:,:,:,iBlock)
             Vel_IDGB(iFluid,y_,:,:,:,iBlock) = &
                  State_VGB(iRhoUy,:,:,:,iBlock)/State_VGB(iRho,:,:,:,iBlock)
             Vel_IDGB(iFluid,z_,:,:,:,iBlock) = &
                  State_VGB(iRhoUz,:,:,:,iBlock)/State_VGB(iRho,:,:,:,iBlock)
          else 
             ! Convert cell center velocity in xyz direction into grid 
             ! aligned direction.
             do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
                ijkDir_DD = 0

                ijkDir_DD(:,x_) = Xyz_DGB(:,min(i+iDim_,iMax),j,k,iBlock) - &
                     Xyz_DGB(:,max(i-iDim_,iMin),j,k,iBlock)
                ijkDir_DD(:,x_) = ijkDir_DD(:,x_)/sqrt(sum(ijkDir_DD(:,x_)**2))

                if(nDim>1) then
                   ijkDir_DD(:,y_) = Xyz_DGB(:,i,min(j+jDim_,jMax),k,iBlock) - &
                        Xyz_DGB(:,i,max(j-jDim_,jMin),k,iBlock)
                   ijkDir_DD(:,y_) = ijkDir_DD(:,y_)/sqrt(sum(ijkDir_DD(:,y_)**2))
                endif

                if(nDim>2) then
                   ijkDir_DD(:,z_) = Xyz_DGB(:,i,j,min(k+kDim_,kMax),iBlock) - &
                        Xyz_DGB(:,i,j,max(k-kDim_,kMin),iBlock)
                   ijkDir_DD(:,z_) = ijkDir_DD(:,z_)/sqrt(sum(ijkDir_DD(:,z_)**2))
                endif

                do iDim = 1, MaxDim
                   Vel_IDGB(iFluid,iDim,i,j,k,iBlock) =  &                       
                        sum(State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock)* &
                        ijkDir_DD(:,iDim))/State_VGB(iRho,i,j,k,iBlock)
                enddo

             enddo; enddo; enddo

          endif ! IsCartesian
       enddo ! iFluid
    enddo ! iBlock

    call test_stop(NameSub, DoTest)
  end subroutine calc_cell_norm_velocity
  !============================================================================

  subroutine calc_face_div_u(iBlock)
    ! This subroutine 'estimates' div(V)*dl at the cell faces, where
    ! 'dl' is the cell size. 
    ! The algorithm is implemented based on the paper of 
    ! P. McCorquodale and P. Colella (2010). See section 2.52 of this paper 
    ! for more details. 

    use ModAdvance, ONLY: State_VGB, FaceDivU_IX, FaceDivU_IY, FaceDivU_IZ, &
         Vel_IDGB
    use BATL_size,   ONLY: nDim, MaxDim, iDim_, jDim_, kDim_
    use BATL_lib, ONLY: Xyz_DGB, IsCartesian
    use ModMain,  ONLY: iMinFace, iMaxFace, jMinFace, jMaxFace, kMinFace, &
         kMaxFace, nIFace, nJFace, nKFace
    integer, intent(in)::iBlock

    integer :: iRho, iRhoUx, iRhoUy, iRhoUz
    integer :: iFluid, iFace, jFace, kFace
    integer :: iMin, iMax, jMin, jMax, kMin, kMax, i, j, k, iDim
    real :: ijkDir_DD(MaxDim, MaxDim)

    real:: Vel_DG(x_:z_,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    character(len=*), parameter:: NameSub = 'calc_face_div_u'
    !------------------------------------------------------------------------

    call timing_start(NameSub)

    iMin = 1 - nG;       iMax = nI + nG
    jMin = 1 - nG*jDim_; jMax = nJ + nG*jDim_
    kMin = 1 - nG*kDim_; kMax = nK + nG*kDim_

    do iFluid = 1, nFluid
       iRho = iRho_I(iFluid)
       iRhoUx = iRhoUx_I(iFluid)
       iRhoUy = iRhoUy_I(iFluid)
       iRhoUz = iRhoUz_I(iFluid)

       Vel_DG = Vel_IDGB(iFluid,:,:,:,:,iBlock)

       ! Assume dl ~ dx ~ dy ~ dz, then FaceDivU_IX/Y/Z = div(U)*dl
       ! See eq(35) of P. McCorquodale and P. Colella (2010)
       do kFace=kMinFace,kMaxFace; do jFace=jMinFace,jMaxFace; do iFace=1,nIFace
          FaceDivU_IX(iFluid,iFace,jFace,kFace) = &
               (Vel_DG(x_,iFace,jFace,kFace) - &
               Vel_DG(x_,iFace-1,jFace,kFace)) 

          if(nDim>1) FaceDivU_IX(iFluid,iFace,jFace,kFace) = &
               FaceDivU_IX(iFluid,iFace,jFace,kFace)+ &
               (Vel_DG(y_,iFace,jFace+1,kFace) - &
               Vel_DG(y_,iFace,jFace-1,kFace) + &
               Vel_DG(y_,iFace-1,jFace+1,kFace) - &
               Vel_DG(y_,iFace-1,jFace-1,kFace))/4


          if(nDim>2) FaceDivU_IX(iFluid,iFace,jFace,kFace) = &
               FaceDivU_IX(iFluid,iFace,jFace,kFace)+ &
               (Vel_DG(z_,iFace,jFace,kFace+1) - &
               Vel_DG(z_,iFace,jFace,kFace-1) + &
               Vel_DG(z_,iFace-1,jFace,kFace+1) - &
               Vel_DG(z_,iFace-1,jFace,kFace-1))/4

       enddo; enddo; enddo

       If(nDim>1) then
          do kFace=kMinFace,kMaxFace; do jFace=1,nJFace; do iFace=iMinFace,iMaxFace
             FaceDivU_IY(iFluid,iFace,jFace,kFace) = &
                  (Vel_DG(y_,iFace,jFace,kFace) - &
                  Vel_DG(y_,iFace,jFace-1,kFace)) + &
                  (Vel_DG(x_,iFace+1,jFace,kFace) - &
                  Vel_DG(x_,iFace-1,jFace,kFace) + &
                  Vel_DG(x_,iFace+1,jFace-1,kFace) - &
                  Vel_DG(x_,iFace-1,jFace-1,kFace))/4

             if(nDim>2) FaceDivU_IY(iFluid,iFace,jFace,kFace) = &
                  FaceDivU_IY(iFluid,iFace,jFace,kFace) + &
                  (Vel_DG(z_,iFace,jFace,kFace+1) - &
                  Vel_DG(z_,iFace,jFace,kFace-1) + &
                  Vel_DG(z_,iFace,jFace-1,kFace+1) - &
                  Vel_DG(z_,iFace,jFace-1,kFace-1))/4
          enddo; enddo; enddo
       endif


       if(nDim>2) then
          do kFace=1,nKFace; do jFace=jMinFace,jMaxFace; do iFace=iMinFace,iMaxFace
             FaceDivU_IZ(iFluid,iFace,jFace,kFace) = &
                  (Vel_DG(z_,iFace,jFace,kFace) - &
                  Vel_DG(z_,iFace,jFace,kFace-1)) + &
                  (Vel_DG(x_,iFace+1,jFace,kFace) - &
                  Vel_DG(x_,iFace-1,jFace,kFace) + &
                  Vel_DG(x_,iFace+1,jFace,kFace-1) - &
                  Vel_DG(x_,iFace-1,jFace,kFace-1))/4 + &
                  (Vel_DG(y_,iFace,jFace+1,kFace) - &
                  Vel_DG(y_,iFace,jFace-1,kFace) + &
                  Vel_DG(y_,iFace,jFace+1,kFace-1) - &
                  Vel_DG(y_,iFace,jFace-1,kFace-1))/4
          enddo; enddo; enddo
       endif

    enddo ! iFluid

    call timing_stop(NameSub)
  end subroutine calc_face_div_u
  !============================================================================

  subroutine limiter_mp(lMin, lMax, Cell_I, Cell2_I, iVar)

    integer, intent(in):: lMin, lMax  ! face index range, e.g. 1...nI+1
    real,    intent(in):: Cell_I(1-nG:MaxIJK+nG)
    real,    intent(in):: Cell2_I(1-nG:MaxIJK+nG)

    integer, optional, intent(in):: iVar        ! variable index

    ! Apply 5th order MP limiter to calculate face values.
    ! Use Cell_I to calculate FaceL_I and Cell2_I for FaceR_I.
    ! If iVar is present, apply the positivity check

    ! Coefficient for 5th order accurate interpolation
    real, parameter:: &
         c1 = 2/60., c2 = -13/60., c3 = 47/60., c4 = 27/60., c5 = -3/60.

    ! If the cell center value is used, the interpolatio for the face value
    ! is implemented with the coefficients below
    real, parameter:: &
         d1 = 3./128, d2 = -20./128, d3 = 90./128, d4=60./128, d5=-5./128

    real, parameter:: cFourThird = 4./3., c6 = 0.6

    ! Cell centered values at l, l+1, l+2, l-1, l-2
    real:: Cell, Cellp, Cellpp, Cellm, Cellmm

    ! Second derivatives
    real:: D2_I(-3:MaxIJK+4)

    ! Limited second derivatives at l+1/2 and l-1/2
    real:: D2p, D2m

    ! Various face values
    real:: FaceOrig, FaceMp, UpperLimit, Average, Median, LargeCurve
    real:: FaceMin, FaceMax

    integer:: l

    real:: sMin, sMax, SoundMin, diff1, diff2_I(3)

    real:: beta = 0.3, r0, RatioVel, FaceLowOrder

    character(len=*), parameter:: NameSub = 'limiter_mp'
    !--------------------------------------------------------------------------

    ! Second derivative based on cell values (3 ghost cells are needed)
    do l = lMin-2, lMax+1
       D2_I(l) = Cell_I(l+1) - 2*Cell_I(l) + Cell_I(l-1)
    end do

    ! Limit left face first. Loop index l is for cell center, and face l+1/2
    do l = lMin-1, lMax-1

       Cellmm = Cell_I(l-2)
       Cellm  = Cell_I(l-1)
       Cell   = Cell_I(l)
       Cellp  = Cell_I(l+1)
       Cellpp = Cell_I(l+2)

       if(UseLowOrder) then
          if(LowOrderCrit_I(l+1) >=1) then
             CYCLE
          else
             FaceLowOrder = FaceL_I(l+1)
          endif

       endif

       if(UseFDFaceFlux) then
          ! Get the 5th order face value.
          FaceOrig = d1*Cellmm + d2*Cellm + d3*Cell + d4*Cellp + d5*Cellpp
       else
          ! 5th order interpolation
          FaceOrig = c1*Cellmm + c2*Cellm + c3*Cell + c4*Cellp + c5*Cellpp
       end if

       ! This is a quick check if there is a need to do any limiting
       FaceMp = Cell + minmod(Cellp - Cell, 4*(Cell - Cellm))

       if( (FaceOrig - Cell)*(FaceOrig - FaceMp) <= 1e-12)then
          FaceL_I(l+1) = FaceOrig
       else

          UpperLimit = Cell + 4*(Cell - Cellm)
          Average    = 0.5*(Cell + Cellp)

          diff1 = max(abs(Cell2_I(l)-Cell2_I(l-1)),1e-14)
          diff2_I(1) = abs(0.5*(Cell2_I(l-2) + Cell2_I(l) - 2*Cell2_I(l-1)))
          diff2_I(2) = abs(0.5*(Cell2_I(l-1) + Cell2_I(l+1) - 2*Cell2_I(l)))
          diff2_I(3) = abs(0.5*(Cell2_I(l) + Cell2_I(l+2) - 2*Cell2_I(l+1)))

          if(UseAccurateExtremum .or. maxval(diff2_I)/diff1 <0.5) then
             D2p = minmod4(4*D2_I(l) - D2_I(l+1),4*D2_I(l+1) - D2_I(l), &
                  D2_I(l), D2_I(l+1))
             D2m = minmod4(4*D2_I(l) - D2_I(l-1),4*D2_I(l-1) - D2_I(l), &
                  D2_I(l), D2_I(l-1))

             Median     = Average - 0.5*D2p
             LargeCurve = Cell + 0.5*(Cell - Cellm) + cFourThird*D2m
          else 
             Median     = Cell 
             LargeCurve = UpperLimit 
          endif

          ! Note: FaceMin <= Cell, FaceMax >= Cell, so FaceMin <= FaceMax
          FaceMin = max(min(Cell, Cellp, Median), &
               min(Cell, UpperLimit, LargeCurve))

          FaceMax = min(max(Cell, Cellp, Median), &
               max(Cell, UpperLimit, LargeCurve))

          ! FaceL = median(FaceOrig, FaceMin, FaceMax)
          FaceL_I(l+1) = min(FaceMax, max(FaceMin, FaceOrig))

       end if


       if(UseLowOrder) then
          FaceL_I(l+1) = (1-LowOrderCrit_I(l+1))*FaceL_I(l+1) + &
               LowOrderCrit_I(l+1)*FaceLowOrder                    
       endif

       ! If the face value is a very small fraction of the cell
       ! then switch to first order scheme. This can occur at
       ! a shock or a smooth minimum very close to zero.
       if(present(iVar))then
          if(DefaultState_V(iVar) > 0.0 .and. &
               FaceL_I(l+1) < c6*Cell) FaceL_I(l+1) = Cell
       end if

    end do

    do l = lMin-2, lMax+1
       D2_I(l) = Cell2_I(l+1) - 2*Cell2_I(l) + Cell2_I(l-1)
    end do

    ! Limit right face. Loop index l is for cell center, and face l-1/2
    do l = lMin, lMax

       Cellmm = Cell2_I(l-2)
       Cellm  = Cell2_I(l-1)
       Cell   = Cell2_I(l)
       Cellp  = Cell2_I(l+1)
       Cellpp = Cell2_I(l+2)


       if(UseLowOrder) then
          if(LowOrderCrit_I(l) >=1) then
             CYCLE
          else
             FaceLowOrder = FaceR_I(l)
          endif
       endif

       if(UseFDFaceFlux) then
          ! FaceOrig
          FaceOrig = d1*Cellpp + d2*Cellp + d3*Cell + d4* Cellm + d5*Cellmm
       else
          ! 5th order interpolation
          FaceOrig = c1*Cellpp + c2*Cellp + c3*Cell + c4*Cellm + c5*Cellmm
       end if

       ! This is a quick check if there is a need to do any limiting
       FaceMp = Cell + minmod(Cellm - Cell, 4*(Cell - Cellp))

       if( (FaceOrig - Cell)*(FaceOrig - FaceMp) <= 1e-12)then
          FaceR_I(l) = FaceOrig
       else

          UpperLimit = Cell + 4*(Cell - Cellp)
          Average    = 0.5*(Cell + Cellm)

          diff1 = max(abs(Cell2_I(l)-Cell2_I(l-1)),1e-14)
          diff2_I(1) = abs(0.5*(Cell2_I(l-2) + Cell2_I(l) - 2*Cell2_I(l-1)))
          diff2_I(2) = abs(0.5*(Cell2_I(l-1) + Cell2_I(l+1) - 2*Cell2_I(l)))
          diff2_I(3) = abs(0.5*(Cell2_I(l) + Cell2_I(l+2) - 2*Cell2_I(l+1)))

          if(UseAccurateExtremum .or. maxval(diff2_I)/diff1 <0.5 ) then
             D2p = minmod4(4*D2_I(l) - D2_I(l+1),4*D2_I(l+1) - D2_I(l), &
                  D2_I(l), D2_I(l+1))
             D2m = minmod4(4*D2_I(l) - D2_I(l-1),4*D2_I(l-1) - D2_I(l), &
                  D2_I(l), D2_I(l-1))

             Median     = Average - 0.5*D2m
             LargeCurve = Cell + 0.5*(Cell - Cellp) + cFourThird*D2p
          else
             Median     = Cell
             LargeCurve = UpperLimit
          endif

          ! Note: FaceMin <= Cell, FaceMax >= Cell, so FaceMin <= FaceMax
          FaceMin = max(min(Cell, Cellm, Median), &
               min(Cell, UpperLimit, LargeCurve))

          FaceMax = min(max(Cell, Cellm, Median), &
               max(Cell, UpperLimit, LargeCurve))

          ! FaceR = median(FaceOrig, FaceMin, FaceMax)
          FaceR_I(l) = min(FaceMax, max(FaceMin, FaceOrig))

       end if

       if(UseLowOrder) then
          FaceR_I(l) = (1-LowOrderCrit_I(l))*FaceR_I(l) + &
               LowOrderCrit_I(l)*FaceLowOrder
       endif

       ! Check fraction limit
       if(present(iVar))then
          if(DefaultState_V(iVar) > 0.0 .and. &
               FaceR_I(l) < c6*Cell) FaceR_I(l) = Cell
       end if

    end do

  contains
    !==========================================================================
    real function minmod(a, b)
      real, intent(in):: a, b
      !------------------------------------------------------------------------
      minmod = (sign(0.5,a) + sign(0.5,b))*min(abs(a), abs(b))
    end function minmod
    !==========================================================================

    real function minmod4(a, b, c, d)
      real, intent(in):: a, b, c, d
      real:: SignSum

      !------------------------------------------------------------------------
      SignSum = sign(0.25,a) + sign(0.25,b) + sign(0.25,c) + sign(0.25,d)
      if(abs(SignSum) < 0.9)then
         minmod4 = 0.0
      else
         minmod4 = SignSum*min(abs(a), abs(b), abs(c), abs(d))
      end if
    end function minmod4
    !==========================================================================

  end subroutine limiter_mp
  !============================================================================
  subroutine calc_cweno_weight(lMin, lMax)

    integer, intent(in):: lMin, lMax

    ! Constants in eq(24)
    real, parameter:: c7over120 = 7./120, c1over6 = 1./6, c1over3 = 1./3
    real, parameter:: c5over6 = 5./6, c21over40 = 21./40, c11over6 = 11./6
    real, parameter:: c73over120 = 73./120, c7over6 = 7./6, c1over60 = 1./60
    real, parameter:: c13over3 = 13./3

    ! Constants for ECHO scheme(UseFDFaceFlux)
    real, parameter:: c3over8 = 3./8, c1over8 = 1./8, c6over8 = 3./4
    real, parameter:: c10over8 = 10./8, c15over8 = 15./8
    real, parameter:: c3over64 = 3./64, c1over16 = 1./16
    real, parameter:: c15over32 = 15./32, c9over16 = 9./16

    ! Generic cell index
    integer:: l

    ! Scalars for Cell_I(l-2:l+2)
    real:: Cellmm, Cellm, Cell, Cellp, Cellpp

    ! Polynomial coefficients for the linear and quadratic terms for the
    ! three second order polynomals using the l-2:l, l-1:l, l:l+2 stencils
    real:: a1_I(3), a2_I(3)

    ! linear coefficients of four low order polynamials. eq (13)
    real, parameter:: LinearCoeff_I(4) = (/0.125, 0.25, 0.125, 0.5/)

    ! Smoothness indicators. eq (19)
    real:: ISLocal_I(4)
    ! eq (24)
    real:: AlphaIS_I(4), Weight_I(4), w1, w2, w3, w4
    real:: ISmin, ISmax, Epsilon

    character(len=*), parameter:: NameSub = 'calc_cweno_weight'
    !--------------------------------------------------------------------------
    do l = lMin - 1, lMax

       Cellmm = Cell_I(l-2)
       Cellm  = Cell_I(l-1)
       Cell   = Cell_I(l)
       Cellp  = Cell_I(l+1)
       Cellpp = Cell_I(l+2)

       ! eq(8)
       a1_I(1) = 1.5*Cell - 2*Cellm + 0.5*Cellmm
       a2_I(1) = 0.5*Cell -   Cellm + 0.5*Cellmm

       ! eq(9)
       a1_I(2) = 0.5*Cellp - 0.5*Cellm
       a2_I(2) = 0.5*Cellp - Cell + 0.5*Cellm

       ! eq(10)
       a1_I(3) = -1.5*Cell + 2*Cellp - 0.5*Cellpp
       a2_I(3) =  0.5*Cell -   Cellp + 0.5*Cellpp

       ! Calculate indicator of smoothness
       ! eq(20)
       ISLocal_I(1:3) = a1_I**2 + c13over3*a2_I**2

       ! Based on G. Capedeville's code sent to us on Sept 2 2013
       ISmin = minval(ISLocal_I(1:3))
       ISmax = maxval(ISLocal_I(1:3))
       ISLocal_I(4) = IsMax

       ! This expression is from G. Capdeville's code
       Epsilon = sqrt(((ISmin + 1e-12)/(ISmax + 1e-12))**3)

       ! the following is just as fast with NAG 5.1
       ! Epsilon = ((ISmin + 1e-12)/(ISmax + 1e-12))**1.5

       ! This alternative expression comes from:
       ! A.A.I. pEER, et al., Appl. Math. Lett. 22 (2009) 1730-1733
       ! Epsilon = 1e-6*min(1.0, ((ISmin+1e-28)/(ISmax-ISmin+1e-30))) + 1e-99

       ! eq(24)
       AlphaIS_I = LinearCoeff_I/(Epsilon + ISLocal_I)**2
       ! eq(23)
       Weight_I = AlphaIS_I/sum(AlphaIS_I)
       w1 = Weight_I(1)
       w2 = Weight_I(2)
       w3 = Weight_I(3)
       w4 = Weight_I(4)

       if(UseFDFaceFlux) then
          WeightL_II(-2,l) = c3over8*w1 - c3over64*w4
          WeightL_II(-1,l) = -(c10over8*w1 + c1over8*w2 - c1over16*w4)
          WeightL_II( 0,l) = c15over8*w1 + c6over8*w2 + c3over8*w3 + &
               c15over32*w4
          WeightL_II(+1,l) = c3over8*w2 + c6over8*w3 + c9over16*w4
          WeightL_II(+2,l) = -(c1over8*w3 + c3over64*w4)

          WeightR_II(-2,l) = -(c1over8*w1 + c3over64*w4)
          WeightR_II(-1,l) = c6over8*w1 + c3over8*w2 + c9over16*w4
          WeightR_II( 0,l) = c3over8*w1 + c6over8*w2 + c15over8*w3 + &
               c15over32*w4
          WeightR_II(+1,l) = -(c1over8*w2 + c10over8*w3 - c1over16*w4)
          WeightR_II(+2,l) = c3over8*w3 - c3over64*w4

       else
          ! Calculate interpolation weights used in eq (34) to obtain
          ! left and right face values
          WeightL_II(-2,l) = c1over3*w1 - c1over60*w4
          WeightL_II(-1,l) = -(c1over6*w2 + c7over6*w1 + c7over120*w4)
          WeightL_II( 0,l) = c5over6*w2 + c1over3*w3 + c11over6*w1 + &
               c73over120*w4
          WeightL_II(+1,l) = c1over3*w2 + c5over6*w3 + c21over40*w4
          WeightL_II(+2,l) = -(c1over6*w3 + c7over120*w4)

          WeightR_II(+2,l) = c1over3*w3 -c1over60*w4
          WeightR_II(+1,l) = -(c1over6*w2 + c7over6*w3 + c7over120*w4)
          WeightR_II( 0,l) = c5over6*w2 + c1over3*w1 + c11over6*w3 + &
               c73over120*w4
          WeightR_II(-1,l) = c1over3*w2 + c5over6*w1 + c21over40*w4
          WeightR_II(-2,l) = -(c1over6*w1 + c7over120*w4)
       endif
    end do

  end subroutine calc_cweno_weight
  !============================================================================
  subroutine limiter_cweno5(lMin, lMax, Cell_I, Cell2_I, iVar)

    ! G. Capdeville, A central WENO scheme for solving hyperbolic conservation
    ! laws on non-uniform meshes, J. Comput. Phys. 227 (2008) 2977-3014

    integer, intent(in):: lMin, lMax
    real,    intent(in):: Cell_I(1-nG:MaxIJK+nG)
    real,    intent(in):: Cell2_I(1-nG:MaxIJK+nG)
    integer, optional, intent(in):: iVar

    integer:: l

    character(len=*), parameter:: NameSub = 'limiter_cweno5'
    !--------------------------------------------------------------------------
    do l = lMin-1, lMax-1
       if(UseLowOrder)then
          if(UseLowOrder_I(l+1)) CYCLE
       end if

       ! eq (34)
       FaceL_I(l+1) = sum(WeightL_II(-2:2,l)*Cell_I(l-2:l+2))
    end do

    do l = lMin, lMax
       if(UseLowOrder)then
          if(UseLowOrder_I(l)) CYCLE
       end if

       ! eq (34)
       FaceR_I(l)   = sum(WeightR_II(-2:2,l)*Cell2_I(l-2:l+2))
    end do

    if(present(iVar))then
       if (DefaultState_V(iVar) > 0.0) then
          do l = lMin-1, lMax-1
             ! Make sure positive variables remain positive
             if(FaceL_I(l+1) < 0)&
                  FaceL_I(l+1) = 0.5*(Cell_I(l+1) + Cell_I(l))
          end do

          do l = lMin, lMax
             if(FaceR_I(l) < 0) &
                  FaceR_I(l) = 0.5*(Cell_I(l-1) + Cell_I(l))
          end do
       end if
    endif

  end subroutine limiter_cweno5
  !============================================================================
  subroutine limiter_ppm4(lMin, lMax, iVar)

    integer, intent(in):: lMin, lMax  ! face index range, e.g. 1...nI+1
    integer, intent(in):: iVar        ! variable to check for positivity

    ! Apply 4th order PPM limiter as described by
    !
    ! "A high-order finite volume method for hyperbolic conservation laws
    ! on locally-refined grids",
    ! P. McCorquodale and P. Colella, 2010, LBNL document
    !
    ! Input: cell centered primitive variables Cell_I(lMin-nG:lMax-1+nG)
    ! Output: limited 4th order accurate left  StateL_I(lMin:lMax)
    !                          and right face  StateR_I(lMin:lMax)
    !
    ! The code may also set StateR_I(lMin-1) and StateL_I(lMax+1)
    ! (this could be excluded).
    !
    ! For now the implementation is done per variable
    ! Optimization can be done once it works

    ! Various constants
    real, parameter:: c0 = 1e-12, c2 = 1.25, c3 = 0.1, c6=6.0

    ! Second derivative based on cell center values
    real:: D2c_I(2-nG:MaxIJK+nG-1)

    ! Third derivative is needed between cells where D2c is known
    real:: D3Face_I(3-nG:MaxIJK+nG-1), D3max, D3min

    real:: Dfm, Dfp, D2f, D2lim, D2Ratio

    integer:: l

    ! Fourth order interpolation scheme
    ! Fill in lMin-1 and lMax+1, because the limiter needs these face values

    character(len=*), parameter:: NameSub = 'limiter_ppm4'
    !--------------------------------------------------------------------------
    do l = lMin - 1, lMax + 1
       Face_I(l) = c7over12*(Cell_I(l-1) + Cell_I(l)) &
            -      c1over12*(Cell_I(l-2) + Cell_I(l+1))
    end do

    ! Second derivative based on cell values
    do l = lMin+1-nG, lMax-2+nG
       D2c_I(l) = Cell_I(l+1) - 2*Cell_I(l) + Cell_I(l-1)
    end do

    if(UseLimiter4 .and. .not.UseTrueCell)then
       ! Third derivative at face based on cell values
       ! Don't use this in "body" blocks, so the stencil is smaller
       do l = lMin-2, lMax+2
          D3Face_I(l) = D2c_I(l) - D2c_I(l-1)
       end do
    end if

    ! Loop through cells and modify FaceL and FaceR values if needed
    ! Start  at lMin-1 so that FaceL_I(lMin) gets set.
    ! Finish at lMax   so that FaceR_I(lMax) gets set.
    do l = lMin - 1, lMax

       if(UseTrueCell)then
          ! The PPM limiter seems to use these cells only
          if(.not.all(IsTrueCell_I(l-2:l+2))) CYCLE
       end if

       ! Set unlimited values as default
       FaceR_I(l)   = Face_I(l)
       FaceL_I(l+1) = Face_I(l+1)

       ! Definitions at bottom of page 6
       Dfm = Cell_I(l) - Face_I(l)
       Dfp = Face_I(l+1) - Cell_I(l)

       ! Check for local extremum on two different stencils
       ! Eqs. (24) and (25)
       if(Dfm*Dfp < 0 .or. &
            (Cell_I(l+2) - Cell_I(l)  )* &
            (Cell_I(l)   - Cell_I(l-2)) < 0)then

          ! Second derivative based on face value (22.3)
          D2f = c6*(Face_I(l+1) + Face_I(l) - 2*Cell_I(l))

          ! Check if second derivative is almost zero (26.1)
          if(abs(D2f) <= c0*maxval(abs(Cell_I(l-2:l+2))))then
             D2Ratio = 0.0
          else
             ! Get limited value for second derivative (26)
             ! First assume that all second derivatives are positive
             D2lim = min(c2*minval(D2c_I(l-1:l+1)), D2f)
             ! If this is negative, try the all negative case, else set 0
             if(D2lim < 0.0) &
                  D2lim = min(0.0, max(c2*maxval(D2c_I(l-1:l+1)), D2f))

             D2Ratio = D2lim / D2f

             ! If D2Ratio is close to 1, no need to limit
             if(D2Ratio >= 1 - c0) CYCLE

             if(UseLimiter4 .and. .not.UseTrueCell)then
                ! Check 3rd derivative condition (28)
                D3min = minval(D3Face_I(l-1:l+2))
                D3max = maxval(D3Face_I(l-1:l+2))
                if(c3*max(abs(D3max),abs(D3min)) > D3max - D3min) CYCLE
             end if
          end if
          if(Dfm*Dfp <= 0)then
             ! Eqs. (29) and (30)
             FaceR_I(l)   = Cell_I(l) - D2Ratio*Dfm
             FaceL_I(l+1) = Cell_I(l) + D2Ratio*Dfp
          elseif(abs(Dfm) >= 2*abs(Dfp))then
             ! Eq. (31)
             FaceR_I(l)   = Cell_I(l) - 2*(1-D2Ratio)*Dfp - D2Ratio*Dfm
          elseif(abs(Dfp) >= 2*abs(Dfm))then
             ! Eq. (32)
             FaceL_I(l+1) = Cell_I(l) + 2*(1-D2Ratio)*Dfm + D2Ratio*Dfp
          end if
       elseif(abs(Dfm) >= 2*abs(Dfp))then
          ! Eq. (33)
          FaceR_I(l)  = Cell_I(l) - 2*Dfp
       elseif(abs(Dfp) >= 2*abs(Dfm))then
          ! Eq. (34)
          FaceL_I(l+1) = Cell_I(l) + 2*Dfm
       end if

       ! Make sure positive variables remain positive
       if(DefaultState_V(iVar) > 0.0)then
          if(FaceR_I(l) < 0) &
               FaceR_I(l) = 0.5*(Cell_I(l-1) + Cell_I(l))
          if(FaceL_I(l+1) < 0)&
               FaceL_I(l+1) = 0.5*(Cell_I(l+1) + Cell_I(l))
       end if

    end do

  end subroutine limiter_ppm4
  !============================================================================
  ! TVD limiters:
  ! mimod limiter:
  !                slim = minmod(s1,s2)
  !
  ! generalized MC (monotonized central) limiter:
  !                slim = minmod(beta*s1, beta*s2, (s1+s2)/2)
  !
  ! Koren limiter (mc3):
  !                slimL = minmod(beta*s1, beta*s2, (s1+2*s2)/3)
  !                slimR = minmod(beta*s1, beta*s2, (2*s1+s2)/3)
  !
  ! beta-limiter   slim = maxmod(minmod(beta*s1,s2), minmod(beta*s2,s1))
  !
  !                (see C.Hirsch, Numerical computation of
  !                 internal and external flows, Volume 2, page 544-545.)
  !
  ! where s1 and s2 are the unlimited slopes, the minmod function
  ! is zero if the arguments have different signs, otherwise it
  ! select the argument with the smallest absolute value, while
  ! the maxmod function selects the argument with the largest absolute value.
  !
  ! For beta=1.0 the MC and the beta limiters coincide with minmod.
  ! For beta=2.0 the beta limiter becomes the superbee limiter.
  !
  ! Note: the subroutines limiter() and limiter_body() calculate the
  !       HALF of the limited difference, so it can be applied simply as
  !
  !       left_face  = central_value - limited_slope_left
  !       right_face = central_value + limited_slope_right
  subroutine limiter_body(lMin, lMax, Beta)

    integer, intent(in):: lMin, lMax
    real,    intent(in):: Beta

    real,dimension(Hi3_):: dVar1_I, dVar2_I !
    real,dimension(nVar):: dVar1_V, dVar2_V ! unlimited left and right slopes
    integer::l

    character(len=*), parameter:: NameSub = 'limiter_body'
    !--------------------------------------------------------------------------
    select case(TypeLimiter)
    case('beta')
       dVar1_I(1:nVar)=Primitive_VI(:,lMax+1)-Primitive_VI(:,lMax)
       dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
       dVar1_I(Lo3_:Hi3_)=Beta*dVar1_I(Lo2_:Hi2_)
       dVar1_I(1:nVar)=sign(0.25,dVar1_I(1:nVar))
       do l=lMax,lMin-1,-1
          dVar2_I=dVar1_I
          dVar1_I(1:nVar)=Primitive_VI(:,l)-Primitive_VI(:,l-1)
          dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
          dVar1_I(Lo3_:Hi3_)=Beta*dVar1_I(Lo2_:Hi2_)
          dVar1_I(1:nVar)=sign(0.25,dVar1_I(1:nVar))
          if(all(IsTrueCell_I(l-1:l+1)))then
             dVar2_I(1:nVar)=dVar2_I(1:nVar)+dVar1_I(1:nVar)
             dVar2_I(Lo2_:Hi2_)=min(dVar2_I(Lo2_:Hi2_),dVar1_I(Lo3_:Hi3_))
             dVar2_I(Lo3_:Hi3_)=min(dVar1_I(Lo2_:Hi2_),dVar2_I(Lo3_:Hi3_))
             dVar2_I(Lo2_:Hi2_)=max(dVar2_I(Lo2_:Hi2_),dVar2_I(Lo3_:Hi3_))
             dVarLimR_VI(:,l) = dVar2_I(1:nVar)*dVar2_I(Lo2_:Hi2_)
          else
             dVarLimR_VI(:,l) = 0.0
          end if
          dVarLimL_VI(:,l) = dVarLimR_VI(:,l)
       end do
    case('minmod')
       dVar1_I(1:nVar)=Primitive_VI(:,lMax+1)-Primitive_VI(:,lMax)
       dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
       dVar1_I(1:nVar)=sign(0.25,dVar1_I(1:nVar))
       do l=lMax,lMin-1,-1
          dVar2_I(1:Hi2_)=dVar1_I(1:Hi2_)
          dVar1_I(1:nVar)=Primitive_VI(:,l)-Primitive_VI(:,l-1)
          dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
          dVar1_I(1:nVar)=sign(0.25,dVar1_I(1:nVar))
          if(all(IsTrueCell_I(l-1:l+1)))then
             dVar2_I(1:nVar)=dVar2_I(1:nVar)+dVar1_I(1:nVar)
             dVar2_I(Lo2_:Hi2_)=min(dVar2_I(Lo2_:Hi2_),dVar1_I(Lo2_:Hi2_))
             dVarLimR_VI(:,l) = dVar2_I(1:nVar)*dVar2_I(Lo2_:Hi2_)
          else
             dVarLimR_VI(:,l) = 0.0
          end if
          dVarLimL_VI(:,l) = dVarLimR_VI(:,l)
       end do
    case('mc')
       dVar1_V = Primitive_VI(:,lMax+1) - Primitive_VI(:,lMax)
       do l=lMax,lMin-1,-1
          dVar2_V = dVar1_V
          dVar1_V = Primitive_VI(:,l) - Primitive_VI(:,l-1)
          if(all(IsTrueCell_I(l-1:l+1)))then
             dVarLimR_VI(:,l) = &
                  (sign(0.25,dVar1_V)+sign(0.25,dVar2_V))*&
                  min(Beta*abs(dVar1_V), Beta*abs(dVar2_V), &
                  0.5*abs(dVar1_V+dVar2_V))
          else
             dVarLimR_VI(:,l) = 0.0
          end if
          dVarLimL_VI(:,l) = dVarLimR_VI(:,l)
       end do
    case('mc3')
       dVar1_V = Primitive_VI(:,lMax+1) - Primitive_VI(:,lMax)
       do l=lMax,lMin-1,-1
          dVar2_V = dVar1_V
          dVar1_V = Primitive_VI(:,l) - Primitive_VI(:,l-1)
          if(all(IsTrueCell_I(l-1:l+1)))then
             dVarLimR_VI(:,l) = &
                  (sign(0.25,dVar1_V)+sign(0.25,dVar2_V))*&
                  min(Beta*abs(dVar1_V), Beta*abs(dVar2_V), &
                  cThird*abs(2*dVar1_V+dVar2_V))
             dVarLimL_VI(:,l) = &
                  (sign(0.25,dVar1_V)+sign(0.25,dVar2_V))*&
                  min(Beta*abs(dVar1_V), Beta*abs(dVar2_V), &
                  cThird*abs(dVar1_V+2*dVar2_V))
          else
             dVarLimR_VI(:,l) = 0.0
             dVarLimL_VI(:,l) = 0.0
          end if
       end do
    case default
       call stop_mpi('limiter_body: unknown TypeLimiter='//TypeLimiter)
    end select

  end subroutine limiter_body
  !============================================================================

  subroutine limiter(lMin, lMax, Beta)

    integer, intent(in):: lMin, lMax
    real,    intent(in):: Beta

    real,dimension(Hi3_):: dVar1_I, dVar2_I
    real,dimension(nVar):: dVar1_V, dVar2_V
    integer::l

    character(len=*), parameter:: NameSub = 'limiter'
    !--------------------------------------------------------------------------
    select case(TypeLimiter)
    case('beta')
       dVar1_I(1:nVar)=Primitive_VI(:,lMax+1)-Primitive_VI(:,lMax)
       dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
       dVar1_I(Lo3_:Hi3_)=Beta*dVar1_I(Lo2_:Hi2_)
       dVar1_I(1:nVar)=sign(0.25,dVar1_I(1:nVar))
       do l=lMax,lMin-1,-1
          dVar2_I=dVar1_I
          dVar1_I(1:nVar)=Primitive_VI(:,l)-Primitive_VI(:,l-1)
          dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
          dVar1_I(Lo3_:Hi3_)=Beta*dVar1_I(Lo2_:Hi2_)
          dVar1_I(1:nVar)=sign(0.25,dVar1_I(1:nVar))
          dVar2_I(1:nVar)=dVar2_I(1:nVar)+dVar1_I(1:nVar)
          dVar2_I(Lo2_:Hi2_)=min(dVar2_I(Lo2_:Hi2_),dVar1_I(Lo3_:Hi3_))
          dVar2_I(Lo3_:Hi3_)=min(dVar1_I(Lo2_:Hi2_),dVar2_I(Lo3_:Hi3_))
          dVar2_I(Lo2_:Hi2_)=max(dVar2_I(Lo2_:Hi2_),dVar2_I(Lo3_:Hi3_))
          dVarLimR_VI(:,l) = dVar2_I(1:nVar)*dVar2_I(Lo2_:Hi2_)

          dVarLimL_VI(:,l) = dVarLimR_VI(:,l)
       end do
    case('minmod')
       dVar1_I(1:nVar)=Primitive_VI(:,lMax+1)-Primitive_VI(:,lMax)
       dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
       dVar1_I(1:nVar)=sign(0.25,dVar1_I(1:nVar))
       do l=lMax,lMin-1,-1
          dVar2_I(1:Hi2_)=dVar1_I(1:Hi2_)
          dVar1_I(1:nVar)=Primitive_VI(:,l)-Primitive_VI(:,l-1)
          dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
          dVar1_I(1:nVar)=sign(0.25,dVar1_I(1:nVar))
          dVar2_I(1:nVar)=dVar2_I(1:nVar)+dVar1_I(1:nVar)
          dVar2_I(Lo2_:Hi2_)=min(dVar2_I(Lo2_:Hi2_),dVar1_I(Lo2_:Hi2_))
          dVarLimR_VI(:,l) = dVar2_I(1:nVar)*dVar2_I(Lo2_:Hi2_)

          dVarLimL_VI(:,l) = dVarLimR_VI(:,l)
       end do
    case('mc')
       ! Calculate rightmost unlimited slope
       dVar1_V = Primitive_VI(:,lMax+1) - Primitive_VI(:,lMax)
       do l=lMax,lMin-1,-1
          ! Propagate old left slope to become the right slope
          dVar2_V = dVar1_V
          ! Calculate left slope
          dVar1_V = Primitive_VI(:,l) - Primitive_VI(:,l-1)
          ! Calculate the limited slope
          dVarLimR_VI(:,l) = (sign(0.25,dVar1_V)+sign(0.25,dVar2_V))* &
               min(Beta*abs(dVar1_V), Beta*abs(dVar2_V), &
               0.5*abs(dVar1_V+dVar2_V))

          dVarLimL_VI(:,l) = dVarLimR_VI(:,l)
       end do
    case('mc3')
       ! Calculate rightmost unlimited slope
       dVar1_V = Primitive_VI(:,lMax+1) - Primitive_VI(:,lMax)
       do l=lMax,lMin-1,-1
          ! Propagate old left slope to become the right slope
          dVar2_V = dVar1_V
          ! Calculate left slope
          dVar1_V = Primitive_VI(:,l) - Primitive_VI(:,l-1)
          ! Calculate the limited slopes
          dVarLimR_VI(:,l) = (sign(0.25,dVar1_V)+sign(0.25,dVar2_V))* &
               min(Beta*abs(dVar1_V), Beta*abs(dVar2_V), &
               cThird*abs(2*dVar1_V+dVar2_V))
          dVarLimL_VI(:,l) = (sign(0.25,dVar1_V)+sign(0.25,dVar2_V))* &
               min(Beta*abs(dVar1_V), Beta*abs(dVar2_V), &
               cThird*abs(dVar1_V+2*dVar2_V))
       end do
    case default
       call stop_mpi('limiter: unknown TypeLimiter='//TypeLimiter)
    end select

  end subroutine limiter
  !============================================================================

  subroutine correct_monotone_restrict(iBlock)

    ! Correct the result of the first order monotone restriction by modifying
    ! the coarse ghost cell values such that the slope remains the same
    ! while the cell center is moved from the fine cell distance to the
    ! coarse cell distance.
    !
    ! This operation should be performed at most once per exchange messages.
    ! To avoid multiple modifications in schemes which switch off some blocks,
    ! the correction is not done if any of the finer block neighbors are unused

    use ModSize
    use ModVarIndexes, ONLY: DefaultState_V, nVar, &
         iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I
    use ModAdvance,    ONLY: State_VGB
    use ModParallel,   ONLY: neiLEV, &
         neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth, &
         neiBtop, neiBbot, neiBeast, neiBwest, neiBnorth, neiBsouth, &
         neiPtop, neiPbot, neiPeast, neiPwest, neiPnorth, neiPsouth
    use BATL_lib,  ONLY: Unused_BP

    ! For debugging

    integer, intent(in) :: iBlock
    integer             :: i, j, k

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'correct_monotone_restrict'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(all(neiLEV(:,iBlock) /= -1))RETURN

    if(DoTest)write(*,*)NameSub, ' state before: ',&
         State_VGB(iVarTest, nI:nI+1, jTest, kTest, iBlock)

    if(.not.DoLimitMomentum)then
       ! Convert momenta to velocities (that will be limited)
       do k = k0_, nKp1_; do j = j0_, nJp1_; do i = 0, nI+1
          State_VGB(iRhoUx_I,i,j,k,iBlock)=State_VGB(iRhoUx_I,i,j,k,iBlock) &
               / State_VGB(iRho_I,i,j,k,iBlock)
          State_VGB(iRhoUy_I,i,j,k,iBlock)=State_VGB(iRhoUy_I,i,j,k,iBlock) &
               / State_VGB(iRho_I,i,j,k,iBlock)
          State_VGB(iRhoUz_I,i,j,k,iBlock)=State_VGB(iRhoUz_I,i,j,k,iBlock) &
               / State_VGB(iRho_I,i,j,k,iBlock)
       end do; end do; end do
    end if

    if(neiLeast(iBlock) == -1)then
       if(        .not.Unused_BP(neiBeast(1,iBlock),neiPeast(1,iBlock)) &
            .and. .not.Unused_BP(neiBeast(2,iBlock),neiPeast(2,iBlock)) &
            .and. .not.Unused_BP(neiBeast(3,iBlock),neiPeast(3,iBlock)) &
            .and. .not.Unused_BP(neiBeast(4,iBlock),neiPeast(4,iBlock)))then
          do k=1,nK;do j=1,nJ
             State_VGB(1:nVar,0,j,k,iBlock) = &
                  State_VGB(1:nVar,0,j,k,iBlock) + cThird*(&
                  State_VGB(1:nVar,0,j,k,iBlock) - &
                  State_VGB(1:nVar,1,j,k,iBlock))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,0,j,k,iBlock) = &
                  max(State_VGB(1:nVar,0,j,k,iBlock), 1e-30)
          end do; end do
       end if
    end if
    if(neiLwest(iBlock) == -1)then
       if(        .not.Unused_BP(neiBwest(1,iBlock),neiPwest(1,iBlock)) &
            .and. .not.Unused_BP(neiBwest(2,iBlock),neiPwest(2,iBlock)) &
            .and. .not.Unused_BP(neiBwest(3,iBlock),neiPwest(3,iBlock)) &
            .and. .not.Unused_BP(neiBwest(4,iBlock),neiPwest(4,iBlock)))then
          do k=1,nK;do j=1,nJ
             State_VGB(1:nVar,nI+1,j,k,iBlock) = &
                  State_VGB(1:nVar,nI+1,j,k,iBlock) + cThird*( &
                  State_VGB(1:nVar,nI+1,j,k,iBlock) - &
                  State_VGB(1:nVar,nI,  j,k,iBlock))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,nI+1,j,k,iBlock) = &
                  max(State_VGB(1:nVar,nI+1,j,k,iBlock), 1e-30)
          end do; end do
       end if
    end if
    if(neiLsouth(iBlock) == -1 .and. nJ > 1)then
       if(        .not.Unused_BP(neiBsouth(1,iBlock),neiPsouth(1,iBlock)) &
            .and. .not.Unused_BP(neiBsouth(2,iBlock),neiPsouth(2,iBlock)) &
            .and. .not.Unused_BP(neiBsouth(3,iBlock),neiPsouth(3,iBlock)) &
            .and. .not.Unused_BP(neiBsouth(4,iBlock),neiPsouth(4,iBlock)))then
          do k=1,nK;do i=1,nI
             State_VGB(1:nVar,i,0,k,iBlock) = &
                  State_VGB(1:nVar,i,0,k,iBlock) + cThird*( &
                  State_VGB(1:nVar,i,0,k,iBlock) - &
                  State_VGB(1:nVar,i,1,k,iBlock))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,i,0,k,iBlock) = &
                  max(State_VGB(1:nVar,i,0,k,iBlock), 1e-30)
          end do; end do
       end if
    end if
    if(neiLnorth(iBlock) == -1 .and. nJ > 1)then
       if(     .not.Unused_BP(neiBnorth(1,iBlock),neiPnorth(1,iBlock)) &
            .and. .not.Unused_BP(neiBnorth(2,iBlock),neiPnorth(2,iBlock)) &
            .and. .not.Unused_BP(neiBnorth(3,iBlock),neiPnorth(3,iBlock)) &
            .and. .not.Unused_BP(neiBnorth(4,iBlock),neiPnorth(4,iBlock)))then
          do k=1,nK;do i=1,nI
             State_VGB(1:nVar,i,nJ+1,k,iBlock) =&
                  State_VGB(1:nVar,i,nJ+1,k,iBlock) + cThird*( &
                  State_VGB(1:nVar,i,nJ+1,k,iBlock) - &
                  State_VGB(1:nVar,i,nJ,k,iBlock))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,i,nJ+1,k,iBlock) = &
                  max(State_VGB(1:nVar,i,nJ+1,k,iBlock), 1e-30)
          end do; end do
       end if
    end if
    if(neiLbot(iBlock) == -1 .and. nK > 1)then
       if(        .not.Unused_BP(neiBbot(1,iBlock),neiPbot(1,iBlock)) &
            .and. .not.Unused_BP(neiBbot(2,iBlock),neiPbot(2,iBlock)) &
            .and. .not.Unused_BP(neiBbot(3,iBlock),neiPbot(3,iBlock)) &
            .and. .not.Unused_BP(neiBbot(4,iBlock),neiPbot(4,iBlock)))then
          do j=1,nJ;do i=1,nI
             State_VGB(1:nVar,i,j,0,iBlock) = &
                  State_VGB(1:nVar,i,j,0,iBlock) + cThird*( &
                  State_VGB(1:nVar,i,j,0,iBlock) - &
                  State_VGB(1:nVar,i,j,1,iBlock))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,i,j,0,iBlock) = &
                  max(State_VGB(1:nVar,i,j,0,iBlock), 1e-30)
          end do; end do
       end if
    end if
    if(neiLtop(iBlock) == -1 .and. nK > 1)then
       if(        .not.Unused_BP(neiBtop(1,iBlock),neiPtop(1,iBlock)) &
            .and. .not.Unused_BP(neiBtop(2,iBlock),neiPtop(2,iBlock)) &
            .and. .not.Unused_BP(neiBtop(3,iBlock),neiPtop(3,iBlock)) &
            .and. .not.Unused_BP(neiBtop(4,iBlock),neiPtop(4,iBlock)))then
          do j=1,nJ;do i=1,nI
             State_VGB(1:nVar,i,j,nK+1,iBlock) = &
                  State_VGB(1:nVar,i,j,nK+1,iBlock) + cThird*(&
                  State_VGB(1:nVar,i,j,nK+1,iBlock) - &
                  State_VGB(1:nVar,i,j,nK,iBlock))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,i,j,nK+1,iBlock) = &
                  max(State_VGB(1:nVar,i,j,nK+1,iBlock), 1e-30)
          end do; end do
       end if
    end if

    if(.not.DoLimitMomentum)then
       ! Convert velocities back to momenta
       do k = k0_, nKp1_; do j = j0_, nJp1_; do i = 0, nI+1
          State_VGB(iRhoUx_I,i,j,k,iBlock)=State_VGB(iRhoUx_I,i,j,k,iBlock) &
               * State_VGB(iRho_I,i,j,k,iBlock)
          State_VGB(iRhoUy_I,i,j,k,iBlock)=State_VGB(iRhoUy_I,i,j,k,iBlock) &
               * State_VGB(iRho_I,i,j,k,iBlock)
          State_VGB(iRhoUz_I,i,j,k,iBlock)=State_VGB(iRhoUz_I,i,j,k,iBlock) &
               * State_VGB(iRho_I,i,j,k,iBlock)
       end do; end do; end do
    end if

    if(DoTest)write(*,*)NameSub, ' state after: ',&
         State_VGB(iVarTest, nI:nI+1, jTest, kTest, iBlock)

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine correct_monotone_restrict
  !============================================================================


end module ModFaceValue
!==============================================================================

