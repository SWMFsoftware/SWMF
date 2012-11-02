!This code is a copyright protected software (c) 2002- University of Michigan
module ModFaceValue

  use ModSize, ONLY: nI, nJ, nK, nG, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       x_, y_, z_, nDim, jDim_, kDim_
  use ModVarIndexes

  use ModMain, ONLY: iTest

  implicit none

  private ! except

  logical, public :: UseAccurateResChange = .false.
  logical, public :: UseTvdResChange      = .true.
  logical, public :: DoLimitMomentum      = .false.
  logical, public :: UseVolumeIntegral4   = .false.
  logical, public :: UseFaceIntegral4     = .false.
  logical, public :: UseLimiter4          = .false.
  integer, public :: nGUsed               = nG

  real,             public :: BetaLimiter = 1.0
  character(len=6), public :: TypeLimiter = 'minmod'

  public :: read_face_value_param, calc_face_value, correct_monotone_restrict

  ! Local variables:
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
  logical:: IsTrueCell_I(1-nG:MaxIJK+nG)
  logical:: UseTrueCell

  ! variables used for TVD limiters
  real:: dVarLimR_VI(1:nVar,0:MaxIJK+1) ! limited slope for right state
  real:: dVarLimL_VI(1:nVar,0:MaxIJK+1) ! limited slope for left state
  real:: Primitive_VI(1:nVar,1-nG:MaxIJK+nG)

  ! variables for the PPM4 limiter
  integer:: iMin, iMax, jMin, jMax, kMin, kMax
  real:: Cell_I(1-nG:MaxIJK+nG)
  real:: Face_I(0:MaxIJK+2)
  real:: FaceL_I(1:MaxIJK+2)
  real:: FaceR_I(0:MaxIJK+1)
  real:: Prim_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

  integer :: iVar

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
    character(len=*), parameter :: NameSub = 'read_face_value_param'
    !--------------------------------------------------------------------------
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
          if(iVar > nVar) call CON_stop(NameSub// &
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
    case default
       call CON_stop(NameSub//' invalid command='//trim(NameCommand))
    end select

  end subroutine read_face_value_param
  !===========================================================================
  subroutine tvd_reschange_body(& 
       Coarse2_V         ,& !State in the coarser ghostcell, 2nd layer
       Coarse1_V         ,& !State in the coarser ghostcell, 1st layer
       Fine1_VII         ,& !States in 4 fine physical cells,1st layer
       Fine2_VII         ,& !States in 4 fine physical cells,2nd layer
       CoarseToFineF_VII ,& !Values at subfaces, in the coarse ghostcell
       FineToCoarseF_VII ,& !Facevalues in phys. cell looking at coarser cell 
       FineF_VII         ,& !Facevalues in phys. cell looking at phys. cell
       IsTrueCoarse2     ,& !True if coarser ghostcell of 2nd layer is true
       IsTrueCoarse1     ,& !True if coarser ghostcell of 1st layer is true
       IsTrueFine1       ,& !True if all physical cells of 1st layer are true
       IsTrueFine2_II)      !True for true physical cell of the 2nd layer
    !_____________!_____________!_______!_______!_
    !             !         CToF!FToC FF!
    ! C2_V        ! C1_V       _!_F1_V__!__F2_V_!_
    !             !         CToF!FToC FF!      
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
    real,dimension(nVar):: GradNormalLtd_V  !Ltd stands for "Limited"

    real:: Beta
    !-------------------------------------------------------------------------
    !Calculate averaged Fine1_VII 
    do iVar=1,nVar
       AveragedFine1_V(iVar) = 0.25*sum(Fine1_VII(iVar,:,:))
    end do
    GradNormal_V= AveragedFine1_V-Coarse1_V

    ! Save gradients squared
    SignGradNormal_V=sign(1.0,GradNormal_V)

    Beta = min(BetaLimiterResChange, BetaLimiter)

    if(IsTrueCoarse2.and.IsTrueCoarse1.and.IsTrueFine1)then
       !Limit gradient in the first coarser cell
       GradNormalLtd_V= SignGradNormal_V*&
            max(0.0,&
            min(Beta*cTwoThird*abs(GradNormal_V),&
            Beta*0.5*SignGradNormal_V*(Coarse1_V-Coarse2_V),&
            cThird*abs(GradNormal_V)+&
            0.25*SignGradNormal_V*(Coarse1_V-Coarse2_V)))

       do j2=1,2;do i2=1,2
          !Limit transverse gradients, if they are larger than the normal one
          !The unlimited transverse gradients are Fine1_VII-AveragedFine1_V
          CoarseToFineF_VII(:,i2,j2)=Coarse1_V+GradNormalLtd_V+&
               sign(min(abs(GradNormalLtd_V), &
               abs(Fine1_VII(:,i2,j2)-AveragedFine1_V)),&
               Fine1_VII(:,i2,j2)-AveragedFine1_V)
       end do;end do
    else
       do j2=1,2;do i2=1,2
          !First order scheme
          CoarseToFineF_VII(:,i2,j2)=Coarse1_V
       end do;end do
    end if
    if(.not.(IsTrueCoarse1.and.IsTrueFine1))then
       do j2=1,2;do i2=1,2
          !First order scheme
          FineToCoarseF_VII(:,i2,j2)=Fine1_VII(:,i2,j2)
          FineF_VII(:,i2,j2)=Fine1_VII(:,i2,j2)
       end do;end do
    else
       do j2=1,2;do i2=1,2
          if(IsTrueFine2_II(i2,j2))then
             !Limit gradient in the first layer of finer cells
             GradNormalLtd_V = SignGradNormal_V*&
                  max(0.0,&
                  min(Beta*cThird*abs(GradNormal_V),&
                  SignGradNormal_V*(Fine1_VII(:,i2,j2)-Coarse1_V),&
                  Beta*0.5*SignGradNormal_V*&
                  (Fine2_VII(:,i2,j2)-Fine1_VII(:,i2,j2)),&
                  cSixth*abs(GradNormal_V) + 0.25*SignGradNormal_V*&
                  (Fine2_VII(:,i2,j2)-Fine1_VII(:,i2,j2))))
          else
             !First order scheme
             GradNormalLtd_V=0.0
          end if
          FineToCoarseF_VII(:,i2,j2)=Fine1_VII(:,i2,j2)-GradNormalLtd_V
          FineF_VII(:,i2,j2)=Fine1_VII(:,i2,j2)+GradNormalLtd_V
       end do;end do
    end if
  end subroutine tvd_reschange_body
  !===========================================================================
  subroutine tvd_reschange(&
       Coarse2_V         ,& !State in the coarser ghostcell, 2nd layer
       Coarse1_V         ,& !State in the coarser ghostcell, 1st layer
       Fine1_VII         ,& !States in 4 fine physical cells,1st layer
       Fine2_VII         ,& !States in 4 fine physical cells,2nd layer
       CoarseToFineF_VII ,& !Values at subfaces, in the coarse ghostcell
       FineToCoarseF_VII ,& !Facevalues in the physical cell,
                                !                        looking at the coarser cell 
       FineF_VII)           !Facevalues in the physical cell,
    !   looking at another physical cell

    !_____________!_____________!_______!_______!_
    !             !         CToF!FToC FF!
    ! C2_V        ! C1_V       _!_F1_V__!__F2_V_!_
    !             !         CToF!FToC FF!      
    !_____________!_____________!_F1_V__!__F2_V_!_
    !             !             !       !       !

    real,dimension(nVar),intent(in):: Coarse2_V,Coarse1_V
    real,dimension(nVar,2,2),intent(in):: Fine1_VII,Fine2_VII
    real,dimension(nVar,2,2),intent(inout)::&
         CoarseToFineF_VII ,FineToCoarseF_VII , FineF_VII
    integer::iVar,i2,j2
    real,dimension(nVar):: AveragedFine1_V
    real,dimension(nVar):: GradNormal_V,SignGradNormal_V
    real,dimension(nVar):: GradNormalLtd_V  !Ltd stands for "Limited"

    real :: Beta
    !-------------------------------------------------------------------------

    !Calculate averaged Fine1_VII 
    do iVar=1,nVar
       AveragedFine1_V(iVar) = 0.25*sum(Fine1_VII(iVar,:,:))
    end do
    GradNormal_V= AveragedFine1_V-Coarse1_V

    !Save gradients squared
    SignGradNormal_V=sign(1.0,GradNormal_V)

    Beta = min(BetaLimiterResChange, BetaLimiter)

    !Limit gradient in the first coarser cell
    GradNormalLtd_V= SignGradNormal_V*max(0.0, min( &
         Beta*cTwoThird*abs(GradNormal_V),&
         Beta*0.5*SignGradNormal_V*(Coarse1_V - Coarse2_V), &
         cThird*abs(GradNormal_V) + 0.25*&
         SignGradNormal_V*(Coarse1_V - Coarse2_V)))

    do j2=1,2;do i2=1,2
       !Limit transverse gradients, if they are larger than the normal one
       !Before limiting the transverse gradients are equal to
       !Fine1_VII-AveragedFine1V
       CoarseToFineF_VII(:,i2,j2) = Coarse1_V + GradNormalLtd_V + &
            sign(min(abs(GradNormalLtd_V),&
            abs(Fine1_VII(:,i2,j2) - AveragedFine1_V)),&
            Fine1_VII(:,i2,j2) - AveragedFine1_V)
    end do;end do

    do j2=1,2;do i2=1,2
       !Limit gradient in the first layer of finer cells
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
  !===========================================================================
  subroutine accurate_reschange3d(&
       Coarse2_V         ,& ! State in the coarser ghostcell,  2nd layer
       Coarse1_VII       ,& ! State in the coarser ghostcells, 1st layer
       Fine1_VII         ,& ! States in 4 fine physical cells, 1st layer
       Fine2_VII         ,& ! States in 4 fine physical cells, 2nd layer
       CoarseToFineF_VII ,& ! Values at subfaces, in the coarse ghostcell
       FineToCoarseF_VII ,& ! Facevalues in the physical cell,
                                !   looking at the coarser cell 
       FineF_VII)           ! Facevalues in the physical cell,
    !                         looking at another physical cell


    !             ! C1_V        !       !       !
    !_____________!_____________!_______!_______!_
    !             !         CToF!FToC FF!       !
    ! C2_V        ! C1_V       _!_F1_V__!__F2_V_!_
    !             !         CToF!FToC FF!       !
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
    !-------------------------------------------------------------------------

    !Calculate averaged Fine1_VII 
    do iVar=1,nVar
       AveragedFine1_V(iVar) = 0.25*sum(Fine1_VII(iVar,:,:))
    end do
    GradNormal_V = AveragedFine1_V - Coarse1_VII(:,1,1)

    !Save sign of the gradient
    SignGradNormal_V=sign(1.0,GradNormal_V)

    !Limit gradient in the first coarser cell
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
       Alpha = (FaceAverage - AverageOrig) / &
            ( 0.25*sum( FaceTmp_II ) - AverageOrig)
       Alpha1 = 1.0 - Alpha

       ! Interpolate
       CoarseToFineF_VII(iVar,:,:) = &
            Alpha*FaceTmp_II + Alpha1*CoarseToFineF_VII(iVar,:,:)

    end do

    ! The face is half the distance in the fine cell
    Slope1_V = 0.5*Slope1_V
    do j2=1,2; do i2=1,2
       !Limit gradient in the first layer of finer cells
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
  !===========================================================================
  subroutine accurate_reschange2d(&
       Coarse2_V         ,& ! State in the coarser ghostcell,  2nd layer
       Coarse1_VI        ,& ! State in the coarser ghostcells, 1st layer
       Fine1_VI          ,& ! States in 2 fine physical cells, 1st layer
       Fine2_VI          ,& ! States in 2 fine physical cells, 2nd layer
       CoarseToFineF_VI  ,& ! Values at subfaces, in the coarse ghostcell
       FineToCoarseF_VI  ,& ! Facevalues in the physical cell,
                                !   looking at the coarser cell 
       FineF_VI)            ! Facevalues in the physical cell,
    !                         looking at another physical cell


    !             ! C1_V        !       !       !
    !_____________!_____________!_______!_______!_
    !             !         CToF!FToC FF!       !
    ! C2_V        ! C1_V       _!_F1_V__!__F2_V_!_
    !             !         CToF!FToC FF!       !
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
    !-------------------------------------------------------------------------

    !Calculate averaged Fine1_VI
    do iVar=1,nVar
       AveragedFine1_V(iVar) = 0.5*sum(Fine1_VI(iVar,:))
    end do
    GradNormal_V = AveragedFine1_V - Coarse1_VI(:,1)

    !Save sign of the gradient
    SignGradNormal_V=sign(1.0,GradNormal_V)

    !Limit gradient in the first coarser cell
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
       Alpha = (FaceAverage - AverageOrig) / (0.5*sum(FaceTmp_I) - AverageOrig)
       Alpha1 = 1.0 - Alpha

       ! Interpolate
       CoarseToFineF_VI(iVar,:) = &
            Alpha*FaceTmp_I + Alpha1*CoarseToFineF_VI(iVar,:)

    end do

    ! The face is half the distance in the fine cell
    Slope1_V = 0.5*Slope1_V
    do i2 = 1, 2
       !Limit gradient in the first layer of finer cells
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
  !===========================================================================
  subroutine accurate_reschange1d(&
       Coarse2_V         ,& ! State in the coarser ghostcell,  2nd layer
       Coarse1_V         ,& ! State in the coarser ghostcells, 1st layer
       Fine1_V           ,& ! States in 2 fine physical cells, 1st layer
       Fine2_V           ,& ! States in 2 fine physical cells, 2nd layer
       CoarseToFineF_V   ,& ! Values at face, in the coarse ghostcell
       FineToCoarseF_V   ,& ! Facevalues in the physical cell,
                                !   looking at the coarser cell 
       FineF_V)             ! Facevalues in the physical cell,
    !                         looking at another physical cell


    !_____________!_____________!_______!_______!_
    !             !         CToF!FToC FF!       !
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
    !-------------------------------------------------------------------------

    !Calculate averaged Fine1_VI
    GradNormal_V = Fine1_V - Coarse1_V

    !Save sign of the gradient
    SignGradNormal_V=sign(1.0,GradNormal_V)

    !Limit gradient in the first coarser cell
    Slope1_V = cTwoThird*abs(GradNormal_V)
    Slope2_V = 0.5*SignGradNormal_V*(Coarse1_V - Coarse2_V)

    Beta = min(BetaLimiterResChange, BetaLimiter)

    GradNormalLtd_V = SignGradNormal_V*max(0.0,min( &
         Beta*Slope1_V, Beta*Slope2_V, 0.5*(Slope1_V+Slope2_V)))

    ! Add limited normal gradient to obtain the middle value for the fine face
    CoarseToFineF_V = Coarse1_V + GradNormalLtd_V

    ! The face is half the distance in the fine cell
    Slope1_V = 0.5*Slope1_V
    !Limit gradient in the first layer of finer cells
    Slope2_V = 0.5*SignGradNormal_V*(Fine2_V - Fine1_V)

    ! The first limiting ensures that the FineToCoarse face value
    ! remains between the Fine1 and Coarse values
    GradNormalLtd_V = SignGradNormal_V*max(0.0,min( &
         SignGradNormal_V*(Fine1_V - Coarse1_V), &
         Beta*Slope1_V, Beta*Slope2_V, 0.5*(Slope1_V + Slope2_V)))

    FineToCoarseF_V = Fine1_V - GradNormalLtd_V
    FineF_V         = Fine1_V + GradNormalLtd_V

  end subroutine accurate_reschange1d
  !===========================================================================
  subroutine calc_face_value(DoResChangeOnly, iBlock)

    use ModMultiFluid, ONLY: nIonFluid, iRho, iUx, iUz, iUx_I, iUz_I

    ! The subroutine calculates right and left face values.
    ! If DoResChangeOnly is true, only facevalues next to a coarser 
    ! neighbor block are calculated.
    ! To improve the code stability, ONLY those values are calculated, which
    ! are actually used by the calc_facefluxes.

    use ModMain,     ONLY: nOrder, nOrderProlong, BlkTest, UseB0, &
         UseConstrainB, nIFace, nJFace, nKFace, &
         iMinFace, iMaxFace, jMinFace, jMaxFace, kMinFace, kMaxFace, &
         iMinFace2, iMaxFace2, jMinFace2, jMaxFace2, kMinFace2, kMaxFace2

    use ModGeometry, ONLY : true_cell, body_BLK
    use ModPhysics, ONLY: GammaWave, c2LIGHT, inv_c2LIGHT
    use ModB0
    use ModAdvance, ONLY: State_VGB, Energy_GBI, &
         UseElectronPressure, UseWavePressure, &
         LeftState_VX,      &  ! Face Left  X
         RightState_VX,     &  ! Face Right X
         LeftState_VY,      &  ! Face Left  Y
         RightState_VY,     &  ! Face Right Y
         LeftState_VZ,      &  ! Face Left  Z
         RightState_VZ         ! Face Right Z

    use ModParallel, ONLY : &
         neiLEV,neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth

    use ModEnergy, ONLY: calc_pressure

    logical, intent(in):: DoResChangeOnly
    integer, intent(in):: iBlock

    integer:: i, j, k, iSide, iFluid
    real:: RhoInv

    real:: RhoC2Inv, BxFull, ByFull, BzFull, B2Full, uBC2Inv, Ga2Boris
    real:: B0_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    ! Number of cells needed to get the face values
    integer:: nStencil

    ! Variables related to 4th order finite volume scheme
    real, parameter:: c24th = 1.0/24.0
    real:: Laplace_V(nVar), Laplace

    real:: State_V(nVar), Energy

    logical::DoTest,DoTestMe
    character(len=*), parameter :: NameSub = 'calc_face_value'
    !-------------------------------------------------------------------------
    if(iBlock==BLKtest .and. .not. DoResChangeOnly )then
       call set_oktest('calc_face_value', DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    if(.not.allocated(Primitive_VG))&
         allocate(Primitive_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

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
       end if
    end if

    UsePtotalLimiter = nOrder > 1 .and. nIonFluid == 1 .and. UsePtotalLtd

    if(.not.DoResChangeOnly & !In order not to call it twice
         .and. nOrder > 1   & !Is not needed for nOrder=1
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
                call calc_primitives_boris    !needed for x-faces
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
          if(UseVolumeIntegral4)then
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
                   call calc_primitives_MHD   !needed for x-faces
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
          if(nOrder==2)then
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
       if(nOrderProlong==1 .and..not.UseConstrainB)then

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
          if(nOrder==2)then
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
             if(nJ > 1) call ratio_to_scalar_faceY(iMinFace,iMaxFace,1,nJFace, &
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

  contains

    !==========================================================================
    subroutine logfaceX_to_faceX(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      do iVar=1,nVar
         if(.not.UseLogLimiter_V(iVar))CYCLE

         LeftState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(LeftState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         RightState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(RightState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
      end do

    end subroutine logfaceX_to_faceX
    !==========================================================================
    subroutine logfaceY_to_faceY(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      do iVar=1,nVar
         if(.not.UseLogLimiter_V(iVar))CYCLE

         LeftState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(LeftState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         RightState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(RightState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
      end do

    end subroutine logfaceY_to_faceY
    !==========================================================================
    subroutine logfaceZ_to_faceZ(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      do iVar=1,nVar
         if(.not.UseLogLimiter_V(iVar))CYCLE

         LeftState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(LeftState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         RightState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              exp(RightState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
      end do

    end subroutine logfaceZ_to_faceZ
    !==========================================================================
    subroutine ratio_to_scalar_faceX(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer :: i, j, k
      !------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         LeftState_VX(iVarLimitRatio_I,i,j,k) = &
              LeftState_VX(iVarLimitRatio_I,i,j,k)*LeftState_VX(Rho_,i,j,k)
         RightState_VX(iVarLimitRatio_I,i,j,k) = &
              RightState_VX(iVarLimitRatio_I,i,j,k)*RightState_VX(Rho_,i,j,k)
      end do; end do; end do

    end subroutine ratio_to_scalar_faceX
    !==========================================================================
    subroutine ratio_to_scalar_faceY(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer :: i, j, k
      !------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         LeftState_VY(iVarLimitRatio_I,i,j,k) = &
              LeftState_VY(iVarLimitRatio_I,i,j,k)*LeftState_VY(Rho_,i,j,k)
         RightState_VY(iVarLimitRatio_I,i,j,k) = &
              RightState_VY(iVarLimitRatio_I,i,j,k)*RightState_VY(Rho_,i,j,k)
      end do; end do; end do

    end subroutine ratio_to_scalar_faceY
    !==========================================================================
    subroutine ratio_to_scalar_faceZ(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer :: i, j, k
      !------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         LeftState_VZ(iVarLimitRatio_I,i,j,k) = &
              LeftState_VZ(iVarLimitRatio_I,i,j,k)*LeftState_VZ(Rho_,i,j,k)
         RightState_VZ(iVarLimitRatio_I,i,j,k) = &
              RightState_VZ(iVarLimitRatio_I,i,j,k)*RightState_VZ(Rho_,i,j,k)
      end do; end do; end do

    end subroutine ratio_to_scalar_faceZ
    !==========================================================================
    subroutine ptotal_to_p_faceX(iMin, iMax, jMin, jMax, kMin, kMax)

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

    end subroutine ptotal_to_p_faceX
    !==========================================================================
    subroutine ptotal_to_p_faceY(iMin, iMax, jMin, jMax, kMin, kMax)

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

    end subroutine ptotal_to_p_faceY
    !==========================================================================
    subroutine ptotal_to_p_faceZ(iMin, iMax, jMin, jMax, kMin, kMax)

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

    end subroutine ptotal_to_p_faceZ

    !==========================================================================
    subroutine calc_primitives_boris                                       
      !momentum is limited

      ! rhoU_Boris = rhoU - ((U x B) x B)/c^2
      !            = rhoU + (U B^2 - B U.B)/c^2
      !            = rhoU*(1+BB/(rho*c2)) - B UdotB/c^2
      Primitive_VG(:,i,j,k) = &
           State_VGB(1:nVar,i,j,k,iBlock)
      BxFull = B0_DG(x_,i,j,k) + Primitive_VG(Bx_,i,j,k)
      ByFull = B0_DG(y_,i,j,k) + Primitive_VG(By_,i,j,k)
      BzFull = B0_DG(z_,i,j,k) + Primitive_VG(Bz_,i,j,k)
      B2Full = BxFull**2 + ByFull**2 + BzFull**2
      RhoC2Inv  = 1/(Primitive_VG(rho_,i,j,k)*c2LIGHT)
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
    !=========================================================================
    subroutine get_facex_high(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax

      real, allocatable, save:: State_VX(:,:,:,:)
      !-----------------------------------------------------------------------
      if(TypeLimiter == 'no')then
         do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            LeftState_VX(:,i,j,k) = &
                 c7over12*(Primitive_VG(:,i-1,j,k) + Primitive_VG(:,i,j,k)) - &
                 c1over12*(Primitive_VG(:,i-2,j,k) + Primitive_VG(:,i+1,j,k))
         end do; end do; end do
      else
         do k = kMin, kMax; do j = jMin, jMax
            if(UseTrueCell)then
               IsTrueCell_I(iMin-nG:iMax-1+nG) = &
                    true_cell(iMin-nG:iMax-1+nG,j,k,iBlock)
               Primitive_VI(:,iMin-2:iMax+1)=Primitive_VG(:,iMin-2:iMax+1,j,k)
               call limiter_body(iMin, iMax, BetaLimiter)
               do i = iMin, iMax
                  LeftState_VX(:,i,j,k) =Primitive_VI(:,i-1)+dVarLimL_VI(:,i-1)
                  RightState_VX(:,i,j,k)=Primitive_VI(:,i)  -dVarLimR_VI(:,i)
               end do
            end if

            do iVar = 1, nVar
               ! Copy points along i direction into 1D array
               Cell_I(iMin-nG:iMax-1+nG) = &
                    Primitive_VG(iVar,iMin-nG:iMax-1+nG,j,k)

               if(UseTrueCell)then
                  ! Use the second order face values where high order is skipped
                  FaceL_I(iMin:iMax) = LeftState_VX(iVar,iMin:iMax,j,k)
                  FaceR_I(iMin:iMax) = RightState_VX(iVar,iMin:iMax,j,k)
               end if
               if(nOrder == 4)then
                  call limiter_ppm4(iMin, iMax)
               else
                  call limiter_mp(iMin, iMax, iVar)
               end if
               ! Copy back the results into the 3D arrays
               LeftState_VX(iVar,iMin:iMax,j,k)  = FaceL_I(iMin:iMax)
               RightState_VX(iVar,iMin:iMax,j,k) = FaceR_I(iMin:iMax)
            end do
         end do; end do

      end if

      if(UseFaceIntegral4 .and. nDim>1)then
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
         end do; end do; end do
         if(TypeLimiter /= 'no')then
            State_VX = RightState_VX
            do k=kMinFace,kMaxFace; do j=jMinFace,jMaxFace; do i=iMin,iMax
               Laplace_V = -2*(nDim-1)*State_VX(:,i,j,k) &
                    + State_VX(:,i,j-1,k) + State_VX(:,i,j+1,k)
               if(nK>1) Laplace_V = Laplace_V &
                    + State_VX(:,i,j,k-1) + State_VX(:,i,j,k+1)
               RightState_VX(:,i,j,k) = State_VX(:,i,j,k) - c24th*Laplace_V
            end do; end do; end do
         end if
      end if

      if(TypeLimiter == 'no') RightState_VX = LeftState_VX

      if(DoLimitMomentum)call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceX(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facex_high
    !========================================================================
    subroutine get_facey_high(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax

      real, allocatable, save:: State_VY(:,:,:,:)
      !-----------------------------------------------------------------------
      if(TypeLimiter == 'no')then
         do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
            LeftState_VY(:,i,j,k) = &
                 c7over12*(Primitive_VG(:,i,j-1,k) + Primitive_VG(:,i,j,k)) - &
                 c1over12*(Primitive_VG(:,i,j-2,k) + Primitive_VG(:,i,j+1,k))
         end do; end do; end do
      else
         do k = kMin, kMax; do i = iMin, iMax
            if(UseTrueCell)then
               IsTrueCell_I(jMin-nG:jMax-1+nG) = &
                    true_cell(i,jMin-nG:jMax-1+nG,k,iBlock)
               Primitive_VI(:,jMin-2:jMax+1)=Primitive_VG(:,i,jMin-2:jMax+1,k)
               call limiter_body(jMin, jMax, BetaLimiter)
               do j = jMin, jMax
                  LeftState_VY(:,i,j,k) =Primitive_VI(:,j-1)+dVarLimL_VI(:,j-1)
                  RightState_VY(:,i,j,k)=Primitive_VI(:,j)  -dVarLimR_VI(:,j)
               end do
            end if

            do iVar = 1, nVar
               ! Copy points along j direction into 1D array
               Cell_I(jMin-nG:jMax-1+nG) = &
                    Primitive_VG(iVar,i,jMin-nG:jMax-1+nG,k)

               if(UseTrueCell)then
                  ! Use the second order face values where high order is skipped
                  FaceL_I(jMin:jMax) = LeftState_VY(iVar,i,jMin:jMax,k)
                  FaceR_I(jMin:jMax) = RightState_VY(iVar,i,jMin:jMax,k)
               end if

               if(nOrder == 4)then
                  call limiter_ppm4(jMin, jMax)
               else
                  call limiter_mp(jMin, jMax, iVar)
               end if
               ! Copy back the results into the 3D arrays
               LeftState_VY(iVar,i,jMin:jMax,k)  = FaceL_I(jMin:jMax)
               RightState_VY(iVar,i,jMin:jMax,k) = FaceR_I(jMin:jMax)
            end do
         end do; end do
      end if

      if(UseFaceIntegral4)then
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
         end do; end do; end do
         if(TypeLimiter /= 'no')then
            State_VY = RightState_VY
            do k=kMinFace,kMaxFace; do j=jMin,jMax; do i=iMinFace,iMaxFace; 
               Laplace_V = -2*(nDim-1)*State_VY(:,i,j,k) &
                    + State_VY(:,i-1,j,k) + State_VY(:,i+1,j,k)
               if(nK>1) Laplace_V = Laplace_V &
                    + State_VY(:,i,j,k-1) + State_VY(:,i,j,k+1)
               RightState_VY(:,i,j,k) = State_VY(:,i,j,k) - c24th*Laplace_V
            end do; end do; end do
         end if
      end if

      if(TypeLimiter == 'no') RightState_VY = LeftState_VY
      if(DoLimitMomentum)call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceY(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facey_high
    !========================================================================
    subroutine get_facez_high(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax

      real, allocatable, save:: State_VZ(:,:,:,:)
      !-----------------------------------------------------------------------
      if(TypeLimiter == 'no')then
         do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
            LeftState_VZ(:,i,j,k) = &
                 c7over12*(Primitive_VG(:,i,j,k-1) + Primitive_VG(:,i,j,k)) - &
                 c1over12*(Primitive_VG(:,i,j,k-2) + Primitive_VG(:,i,j,k+1))

            RightState_VZ(:,i,j,k)=LeftState_VZ(:,i,j,k)
         end do; end do; end do
      else
         do j = jMin, jMax; do i = iMin, iMax

            if(UseTrueCell)then
               IsTrueCell_I(kMin-nG:kMax-1+nG) = &
                    true_cell(i,j,kMin-nG:kMax-1+nG,iBlock)
               Primitive_VI(:,kMin-2:kMax+1)=Primitive_VG(:,i,j,kMin-2:kMax+1)
               call limiter_body(kMin, kMax, BetaLimiter)
               do k = kMin, kMax
                  LeftState_VZ(:,i,j,k) =Primitive_VI(:,k-1)+dVarLimL_VI(:,k-1)
                  RightState_VZ(:,i,j,k)=Primitive_VI(:,k)  -dVarLimR_VI(:,k)
               end do
            end if

            do iVar = 1, nVar
               ! Copy points along k direction into 1D array
               Cell_I(kMin-nG:kMax-1+nG) = &
                    Primitive_VG(iVar,i,j,kMin-nG:kMax-1+nG)
               
               if(UseTrueCell)then
                  ! Use the second order face values where high order is skipped
                  FaceL_I(kMin:kMax) = LeftState_VZ(iVar,i,j,kMin:kMax)
                  FaceR_I(kMin:kMax) = RightState_VZ(iVar,i,j,kMin:kMax)
               end if

               if(nOrder == 4)then
                  call limiter_ppm4(kMin, kMax)
               else
                  call limiter_mp(kMin, kMax, iVar)
               end if
               ! Copy back the results into the 3D arrays
               LeftState_VZ(iVar,i,j,kMin:kMax)  = FaceL_I(kMin:kMax)
               RightState_VZ(iVar,i,j,kMin:kMax) = FaceR_I(kMin:kMax)
            end do
         end do; end do
      end if

      if(UseFaceIntegral4)then
         if(.not.allocated(State_VZ)) allocate( &
              State_VZ(nVar,iMinFace2:iMaxFace2,jMinFace2:jMaxFace2,nK+1))

         ! Convert from face averaged to face centered variables (eq 18)
         State_VZ = LeftState_VZ
         do k = kMin,kMax; do j = jMinFace,jMaxFace; do i = iMinFace,iMaxFace
            Laplace_V = -4*State_VZ(:,i,j,k) &
                 + State_VZ(:,i-1,j,k) + State_VZ(:,i+1,j,k) &
                 + State_VZ(:,i,j-1,k) + State_VZ(:,i,j+1,k)
            LeftState_VZ(:,i,j,k) = State_VZ(:,i,j,k) - c24th*Laplace_V
         end do; end do; end do
         if(TypeLimiter /= 'no')then
            State_VZ = RightState_VZ
            do k=kMin,kMax; do j=jMinFace,jMaxFace; do i=iMinFace,iMaxFace
               Laplace_V = -4*State_VZ(:,i,j,k) &
                    + State_VZ(:,i-1,j,k) + State_VZ(:,i+1,j,k) &
                    + State_VZ(:,i,j-1,k) + State_VZ(:,i,j+1,k)
               RightState_VZ(:,i,j,k) = State_VZ(:,i,j,k) - c24th*Laplace_V
            end do; end do; end do
         end if
      end if

      if(DoLimitMomentum)call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceZ(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_facez_high

    !========================================================================
    subroutine get_faceX_first(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      !-----------------------------------------------------------------------
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         LeftState_VX(:,i,j,k)=Primitive_VG(:,i-1,j,k)
         RightState_VX(:,i,j,k)=Primitive_VG(:,i,j,k)              
      end do; end do; end do

      if(DoLimitMomentum)call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceX(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_faceX_first
    !========================================================================
    subroutine get_faceY_first(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      !----------------------------------------------------------------------
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         LeftState_VY(:,i,j,k)=Primitive_VG(:,i,j-1,k)
         RightState_VY(:,i,j,k)=Primitive_VG(:,i,j,k)              
      end do; end do; end do

      if(DoLimitMomentum) call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceY(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_faceY_first
    !========================================================================
    subroutine get_faceZ_first(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      !----------------------------------------------------------------------
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         LeftState_VZ(:,i,j,k)=Primitive_VG(:,i,j,k-1)
         RightState_VZ(:,i,j,k)=Primitive_VG(:,i,j,k)
      end do; end do; end do

      if(DoLimitMomentum)call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceZ(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_faceZ_first
    !==========================================================================
    subroutine calc_primitives_MHD

      Primitive_VG(:,i,j,k) = State_VGB(1:nVar,i,j,k,iBlock)
      RhoInv = 1/Primitive_VG(Rho_,i,j,k)
      Primitive_VG(Ux_:Uz_,i,j,k)=RhoInv*Primitive_VG(RhoUx_:RhoUz_,i,j,k)
      do iFluid = 2, nFluid
         iRho = iRho_I(iFluid); iUx = iUx_I(iFluid); iUz = iUz_I(iFluid)
         RhoInv = 1/Primitive_VG(iRho,i,j,k)
         Primitive_VG(iUx:iUz,i,j,k)=RhoInv*Primitive_VG(iUx:iUz,i,j,k)
      end do

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
           RhoInv*Primitive_VG(iVarLimitRatio_I,i,j,k)

      if(UseLogLimiter)then
         do iVar=1,nVar
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
         BxFull = B0_DX(x_,i,j,k) + LeftState_VX(Bx_,i,j,k)
         ByFull = B0_DX(y_,i,j,k) + LeftState_VX(By_,i,j,k)
         BzFull = B0_DX(z_,i,j,k) + LeftState_VX(Bz_,i,j,k)
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  = inv_c2LIGHT*RhoInv
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
         BxFull = B0_DX(x_,i,j,k) + RightState_VX(Bx_,i,j,k)
         ByFull = B0_DX(y_,i,j,k) + RightState_VX(By_,i,j,k)
         BzFull = B0_DX(z_,i,j,k) + RightState_VX(Bz_,i,j,k)
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  =inv_c2LIGHT*RhoInv
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
    !=========================================================================
    subroutine BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)
      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax

      ! U_Boris=rhoU_Boris/rho
      ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

         ! Left face values
         RhoInv = 1/LeftState_VY(rho_,i,j,k)
         BxFull = B0_DY(x_,i,j,k) + LeftState_VY(Bx_,i,j,k)
         ByFull = B0_DY(y_,i,j,k) + LeftState_VY(By_,i,j,k)
         BzFull = B0_DY(z_,i,j,k) + LeftState_VY(Bz_,i,j,k)
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  =inv_c2LIGHT*RhoInv
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
         BxFull = B0_DY(x_,i,j,k) + RightState_VY(Bx_,i,j,k)
         ByFull = B0_DY(y_,i,j,k) + RightState_VY(By_,i,j,k)
         BzFull = B0_DY(z_,i,j,k) + RightState_VY(Bz_,i,j,k)
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv=inv_c2LIGHT*RhoInv
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
         BxFull = B0_DZ(x_,i,j,k) + LeftState_VZ(Bx_,i,j,k)
         ByFull = B0_DZ(y_,i,j,k) + LeftState_VZ(By_,i,j,k)
         BzFull = B0_DZ(z_,i,j,k) + LeftState_VZ(Bz_,i,j,k)
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  =inv_c2LIGHT*RhoInv
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
         BxFull = B0_DZ(x_,i,j,k) + RightState_VZ(Bx_,i,j,k)
         ByFull = B0_DZ(y_,i,j,k) + RightState_VZ(By_,i,j,k)
         BzFull = B0_DZ(z_,i,j,k) + RightState_VZ(Bz_,i,j,k)
         B2Full = BxFull**2 + ByFull**2 + BzFull**2
         RhoC2Inv  =inv_c2LIGHT*RhoInv
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
    !=======================================================================
    subroutine get_face_accurate1d(iSideIn)
      integer, intent(in):: iSideIn

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
    !=======================================================================
    subroutine get_face_accurate2d(iSideIn)
      integer, intent(in):: iSideIn

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
    !=======================================================================
    subroutine get_face_tvd(iSideIn)
      integer,intent(in)::iSideIn

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
    !=========================================================================
    subroutine get_faceX_second(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      integer::i1, iMinSharp, iMaxSharp
      !----------------------------------------------------------------------
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

    end subroutine get_faceX_second
    !==========================================================================
    subroutine get_faceY_second(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      integer::j1, jMinSharp, jMaxSharp
      !----------------------------------------------------------------------
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

    end subroutine get_faceY_second
    !==========================================================================
    subroutine get_faceZ_second(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer::k1, kMinSharp, kMaxSharp
      !----------------------------------------------------------------------
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

    end subroutine get_faceZ_second

    !========================================================================
    subroutine flatten(Prim_VG)

      use ModMultiFluid, ONLY: iFluid, iRho, iUx, iUy, iUz, iP, select_fluid

      real, intent(in):: Prim_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

      real:: InvRatioRange

      real:: pL, pR, Dp, Ratio, Coef1, Coef

      real:: FlatCoef_I(-1:MaxIJK+2)
      real, allocatable:: FlatCoef_G(:,:,:)

      integer:: i, j, k
      !----------------------------------------------------------------------
      !call timing_start('flatten')

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

      !call timing_stop('flatten')

    end subroutine flatten

  end subroutine calc_face_value

  !===========================================================================
  subroutine limiter_mp(lMin, lMax, iVar)

    integer, intent(in):: lMin, lMax  ! face index range, e.g. 1...nI+1

    integer, intent(in):: iVar        ! variable index

    ! Apply 5th order MP limiter

    ! Coefficient for 5th order accurate interpolation
    real, parameter:: &
         c1 = 2/60., c2 = -13/60., c3 = 47/60., c4 = 27/60., c5 = -3/60.

    real, parameter:: cFourThird = 4./3., c6 = 0.6

    ! Cell centered values at l, l+1, l+2, l-1, l-2
    real:: Cell, Cellp, Cellpp, Cellm, Cellmm

    ! Second derivatives
    real:: D2_I(-1:MaxIJK+2)

    ! Limited second derivatives at l+1/2 and l-1/2
    real:: D2p, D2m

    ! Various face values
    real:: FaceOrig, FaceMp, UpperLimit, Average, Median, LargeCurve
    real:: FaceMin, FaceMax

    integer:: l
    !------------------------------------------------------------------------

    ! Second derivative based on cell values (3 ghost cells are needed)
    do l = lMin-2, lMax+1
       D2_I(l) = Cell_I(l+1) - 2*Cell_I(l) + Cell_I(l-1) 
    end do

    ! Limit left face first. Loop index l is for cell center, and face l+1/2
    do l = lMin-1, lMax-1

       if(UseTrueCell)then
          ! The left face uses l-2:l+2, while the right face need l-1:l+3
          if(.not.all(IsTrueCell_I(l-2:l+3))) CYCLE
       end if

       Cellmm = Cell_I(l-2)
       Cellm  = Cell_I(l-1)
       Cell   = Cell_I(l)
       Cellp  = Cell_I(l+1)
       Cellpp = Cell_I(l+2)

       ! 5th order interpolation
       FaceOrig = c1*Cellmm + c2*Cellm + c3*Cell + c4*Cellp + c5*Cellpp

       ! This is a quick check if there is a need to do any limiting
       FaceMp = Cell + minmod(Cellp - Cell, 4*(Cell - Cellm))

       if( (FaceOrig - Cell)*(FaceOrig - FaceMp) <= 1e-12)then
          FaceL_I(l+1) = FaceOrig

       else

          D2p = minmod4(4*D2_I(l) - D2_I(l+1),4*D2_I(l+1) - D2_I(l), &
               D2_I(l), D2_I(l+1))

          D2m = minmod4(4*D2_I(l) - D2_I(l-1),4*D2_I(l-1) - D2_I(l), &
               D2_I(l), D2_I(l-1))
          
          UpperLimit = Cell + 4*(Cell - Cellm)
          Average    = 0.5*(Cell + Cellp)
          Median     = Average - 0.5*D2p
          LargeCurve = Cell + 0.5*(Cell - Cellm) + cFourThird*D2m

          ! Note: FaceMin <= Cell, FaceMax >= Cell, so FaceMin <= FaceMax
          FaceMin = max(min(Cell, Cellp, Median), &
               min(Cell, UpperLimit, LargeCurve))

          FaceMax = min(max(Cell, Cellp, Median), &
               max(Cell, UpperLimit, LargeCurve))

          ! FaceL = median(FaceOrig, FaceMin, FaceMax)
          FaceL_I(l+1) = min(FaceMax, max(FaceMin, FaceOrig))

       end if

       ! If the face value is a very small fraction of the cell
       ! then switch to first order scheme. This can occur at
       ! a shock or a smooth minimum very close to zero.
       if(DefaultState_V(iVar) > 0.0 .and. &
            FaceL_I(l+1) < c6*Cell) FaceL_I(l+1) = Cell

       !if(iVar == 1 .and. FaceL_I(l+1) < 0.0)then
       !   write(*,*)'!!! Negative FaceL for l=', l
       !   write(*,*)'Cell_I(l-2:l+2)=', Cell_I(l-2:l+2)
       !   write(*,*)'FaceL,FaceOrig,FaceMp=',FaceL_I(l+1), FaceOrig, FaceMp
       !   if( (FaceOrig - Cell)*(FaceOrig - FaceMp) > 1e-12)then
       !      write(*,*)'D2_I(l-1:l+1)=', D2_I(l-1:l+1)
       !      write(*,*)'D2p, D2m     =',D2p, D2m
       !      write(*,*)'UpperLimit   =', UpperLimit 
       !      write(*,*)'Average      =', Average
       !      write(*,*)'Median       =', Median
       !      write(*,*)'LargeCurve   =', LargeCurve
       !      write(*,*)'FaceMin      =', FaceMin
       !      write(*,*)'FaceMax      =', FaceMax
       !   end if
       !end if

    end do

    ! Limit right face. Loop index l is for cell center, and face l-1/2
    do l = lMin, lMax

       if(UseTrueCell)then
          ! The right face uses l-2:l+2, while the left face needs l-3:l+1
          if(.not.all(IsTrueCell_I(l-3:l+2))) CYCLE
       end if

       Cellmm = Cell_I(l-2)
       Cellm  = Cell_I(l-1)
       Cell   = Cell_I(l)
       Cellp  = Cell_I(l+1)
       Cellpp = Cell_I(l+2)

       ! 5th order interpolation
       FaceOrig = c1*Cellpp + c2*Cellp + c3*Cell + c4*Cellm + c5*Cellmm

       ! This is a quick check if there is a need to do any limiting
       FaceMp = Cell + minmod(Cellm - Cell, 4*(Cell - Cellp))

       if( (FaceOrig - Cell)*(FaceOrig - FaceMp) <= 1e-12)then
          FaceR_I(l) = FaceOrig
       else

          D2p = minmod4(4*D2_I(l) - D2_I(l+1),4*D2_I(l+1) - D2_I(l), &
               D2_I(l), D2_I(l+1))

          D2m = minmod4(4*D2_I(l) - D2_I(l-1),4*D2_I(l-1) - D2_I(l), &
               D2_I(l), D2_I(l-1))
          
          UpperLimit = Cell + 4*(Cell - Cellp)
          Average    = 0.5*(Cell + Cellm)
          Median     = Average - 0.5*D2m
          LargeCurve = Cell + 0.5*(Cell - Cellp) + cFourThird*D2p

          ! Note: FaceMin <= Cell, FaceMax >= Cell, so FaceMin <= FaceMax
          FaceMin = max(min(Cell, Cellm, Median), &
               min(Cell, UpperLimit, LargeCurve))

          FaceMax = min(max(Cell, Cellm, Median), &
               max(Cell, UpperLimit, LargeCurve))

          ! FaceR = median(FaceOrig, FaceMin, FaceMax)
          FaceR_I(l) = min(FaceMax, max(FaceMin, FaceOrig))

       end if

       ! Check fraction limit
       if(DefaultState_V(iVar) > 0.0 .and. &
            FaceR_I(l) < c6*Cell) FaceR_I(l) = Cell

       !if(iVar == 1 .and. FaceR_I(l) < 0.0)then
       !   write(*,*)'!!! Negative FaceR for l=', l
       !   write(*,*)'Cell_I(l-2:l+2)=', Cell_I(l-2:l+2)
       !   write(*,*)'FaceR,FaceOrig,FaceMp=',FaceR_I(l), FaceOrig, FaceMp
       !end if

    end do

  contains
    !========================================================================
    real function minmod(a, b)
      real, intent(in):: a, b
      !---------------------------------------------------------------------
      minmod = (sign(0.5,a) + sign(0.5,b))*min(abs(a), abs(b))
    end function minmod

    !========================================================================
    real function minmod4(a, b, c, d)
      real, intent(in):: a, b, c, d
      real:: SignSum
      !---------------------------------------------------------------------
      SignSum = sign(0.25,a) + sign(0.25,b) + sign(0.25,c) + sign(0.25,d)
      if(abs(SignSum) < 0.9)then
         minmod4 = 0.0
      else
         minmod4 = SignSum*min(abs(a), abs(b), abs(c), abs(d))
      end if
    end function minmod4
    !========================================================================

  end subroutine limiter_mp

  !===========================================================================
  subroutine limiter_ppm4(lMin, lMax)

    integer, intent(in):: lMin, lMax  ! face index range, e.g. 1...nI+1

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

    character(len=*), parameter:: NameSub = 'limiter_ppm4'
    !-------------------------------------------------------------------------
    ! Fourth order interpolation scheme
    ! Fill in lMin-1 and lMax+1, because the limiter needs these face values
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
  !===========================================================================
  subroutine limiter_body(lMin, lMax, Beta)

    integer, intent(in):: lMin, lMax
    real,    intent(in):: Beta

    real,dimension(Hi3_):: dVar1_I, dVar2_I ! 
    real,dimension(nVar):: dVar1_V, dVar2_V ! unlimited left and right slopes
    integer::l
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

  !===========================================================================

  subroutine limiter(lMin, lMax, Beta)

    integer, intent(in):: lMin, lMax
    real,    intent(in):: Beta

    real,dimension(Hi3_):: dVar1_I, dVar2_I
    real,dimension(nVar):: dVar1_V, dVar2_V
    integer::l
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

  !=======================================================================

  subroutine correct_monotone_restrict(iBLK)

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
    use ModProcMH, ONLY: iProc
    use ModMain, ONLY: VarTest, ProcTest, BlkTest, jTest, kTest

    integer, intent(in) :: iBLK
    integer             :: i, j, k
    logical :: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'correct_monotone_restrict'
    !--------------------------------------------------------------------------
    if(all(neiLEV(:,iBLK) /= -1))RETURN

    if(iProc == ProcTest .and. iBlk == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    if(DoTestMe)write(*,*)NameSub, ' state before: ',&
         State_VGB(VarTest, nI:nI+1, jTest, kTest, iBlk)

    if(.not.DoLimitMomentum)then
       ! Convert momenta to velocities (that will be limited)
       do k = k0_, nKp1_; do j = j0_, nJp1_; do i = 0, nI+1
          State_VGB(iRhoUx_I,i,j,k,iBLK)=State_VGB(iRhoUx_I,i,j,k,iBLK) &
               / State_VGB(iRho_I,i,j,k,iBLK)
          State_VGB(iRhoUy_I,i,j,k,iBLK)=State_VGB(iRhoUy_I,i,j,k,iBLK) &
               / State_VGB(iRho_I,i,j,k,iBLK)
          State_VGB(iRhoUz_I,i,j,k,iBLK)=State_VGB(iRhoUz_I,i,j,k,iBLK) &
               / State_VGB(iRho_I,i,j,k,iBLK)
       end do; end do; end do
    end if

    if(neiLeast(iBLK) == -1)then
       if(        .not.Unused_BP(neiBeast(1,iBLK),neiPeast(1,iBLK)) &
            .and. .not.Unused_BP(neiBeast(2,iBLK),neiPeast(2,iBLK)) &
            .and. .not.Unused_BP(neiBeast(3,iBLK),neiPeast(3,iBLK)) &
            .and. .not.Unused_BP(neiBeast(4,iBLK),neiPeast(4,iBLK)))then
          do k=1,nK;do j=1,nJ
             State_VGB(1:nVar,0,j,k,iBLK) = &
                  State_VGB(1:nVar,0,j,k,iBLK) + cThird*(&
                  State_VGB(1:nVar,0,j,k,iBLK) - &
                  State_VGB(1:nVar,1,j,k,iBLK))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,0,j,k,iBLK) = &
                  max(State_VGB(1:nVar,0,j,k,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLwest(iBLK) == -1)then
       if(        .not.Unused_BP(neiBwest(1,iBLK),neiPwest(1,iBLK)) &
            .and. .not.Unused_BP(neiBwest(2,iBLK),neiPwest(2,iBLK)) &
            .and. .not.Unused_BP(neiBwest(3,iBLK),neiPwest(3,iBLK)) &
            .and. .not.Unused_BP(neiBwest(4,iBLK),neiPwest(4,iBLK)))then
          do k=1,nK;do j=1,nJ
             State_VGB(1:nVar,nI+1,j,k,iBLK) = &
                  State_VGB(1:nVar,nI+1,j,k,iBLK) + cThird*( &
                  State_VGB(1:nVar,nI+1,j,k,iBLK) - &
                  State_VGB(1:nVar,nI,  j,k,iBLK))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,nI+1,j,k,iBLK) = &
                  max(State_VGB(1:nVar,nI+1,j,k,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLsouth(iBLK) == -1 .and. nJ > 1)then
       if(        .not.Unused_BP(neiBsouth(1,iBLK),neiPsouth(1,iBLK)) & 
            .and. .not.Unused_BP(neiBsouth(2,iBLK),neiPsouth(2,iBLK)) & 
            .and. .not.Unused_BP(neiBsouth(3,iBLK),neiPsouth(3,iBLK)) & 
            .and. .not.Unused_BP(neiBsouth(4,iBLK),neiPsouth(4,iBLK)))then
          do k=1,nK;do i=1,nI
             State_VGB(1:nVar,i,0,k,iBLK) = &
                  State_VGB(1:nVar,i,0,k,iBLK) + cThird*( &
                  State_VGB(1:nVar,i,0,k,iBLK) - &
                  State_VGB(1:nVar,i,1,k,iBLK))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,i,0,k,iBLK) = &
                  max(State_VGB(1:nVar,i,0,k,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLnorth(iBLK) == -1 .and. nJ > 1)then
       if(     .not.Unused_BP(neiBnorth(1,iBLK),neiPnorth(1,iBLK)) &
            .and. .not.Unused_BP(neiBnorth(2,iBLK),neiPnorth(2,iBLK)) &
            .and. .not.Unused_BP(neiBnorth(3,iBLK),neiPnorth(3,iBLK)) &
            .and. .not.Unused_BP(neiBnorth(4,iBLK),neiPnorth(4,iBLK)))then
          do k=1,nK;do i=1,nI
             State_VGB(1:nVar,i,nJ+1,k,iBLK) =&
                  State_VGB(1:nVar,i,nJ+1,k,iBLK) + cThird*( &
                  State_VGB(1:nVar,i,nJ+1,k,iBLK) - &
                  State_VGB(1:nVar,i,nJ,k,iBLK))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,i,nJ+1,k,iBLK) = &
                  max(State_VGB(1:nVar,i,nJ+1,k,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLbot(iBLK) == -1 .and. nK > 1)then
       if(        .not.Unused_BP(neiBbot(1,iBLK),neiPbot(1,iBLK)) &
            .and. .not.Unused_BP(neiBbot(2,iBLK),neiPbot(2,iBLK)) &
            .and. .not.Unused_BP(neiBbot(3,iBLK),neiPbot(3,iBLK)) &
            .and. .not.Unused_BP(neiBbot(4,iBLK),neiPbot(4,iBLK)))then
          do j=1,nJ;do i=1,nI
             State_VGB(1:nVar,i,j,0,iBLK) = &
                  State_VGB(1:nVar,i,j,0,iBLK) + cThird*( &
                  State_VGB(1:nVar,i,j,0,iBLK) - &
                  State_VGB(1:nVar,i,j,1,iBLK))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,i,j,0,iBLK) = &
                  max(State_VGB(1:nVar,i,j,0,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLtop(iBLK) == -1 .and. nK > 1)then
       if(        .not.Unused_BP(neiBtop(1,iBLK),neiPtop(1,iBLK)) &
            .and. .not.Unused_BP(neiBtop(2,iBLK),neiPtop(2,iBLK)) &
            .and. .not.Unused_BP(neiBtop(3,iBLK),neiPtop(3,iBLK)) &
            .and. .not.Unused_BP(neiBtop(4,iBLK),neiPtop(4,iBLK)))then
          do j=1,nJ;do i=1,nI
             State_VGB(1:nVar,i,j,nK+1,iBLK) = &
                  State_VGB(1:nVar,i,j,nK+1,iBLK) + cThird*(&
                  State_VGB(1:nVar,i,j,nK+1,iBLK) - &
                  State_VGB(1:nVar,i,j,nK,iBLK))
             where(DefaultState_V(1:nVar) > 0.0) &
                  State_VGB(1:nVar,i,j,nK+1,iBLK) = &
                  max(State_VGB(1:nVar,i,j,nK+1,iBLK), 1e-30)
          end do; end do
       end if
    end if

    if(.not.DoLimitMomentum)then
       ! Convert velocities back to momenta
       do k = k0_, nKp1_; do j = j0_, nJp1_; do i = 0, nI+1
          State_VGB(iRhoUx_I,i,j,k,iBLK)=State_VGB(iRhoUx_I,i,j,k,iBLK) &
               * State_VGB(iRho_I,i,j,k,iBLK)
          State_VGB(iRhoUy_I,i,j,k,iBLK)=State_VGB(iRhoUy_I,i,j,k,iBLK) &
               * State_VGB(iRho_I,i,j,k,iBLK)
          State_VGB(iRhoUz_I,i,j,k,iBLK)=State_VGB(iRhoUz_I,i,j,k,iBLK) &
               * State_VGB(iRho_I,i,j,k,iBLK)
       end do; end do; end do
    end if

    if(DoTestMe)write(*,*)NameSub, ' state after: ',&
         State_VGB(VarTest, nI:nI+1, jTest, kTest, iBlk)

  end subroutine correct_monotone_restrict

end module ModFaceValue

