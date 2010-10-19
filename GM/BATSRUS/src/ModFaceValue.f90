!^CFG COPYRIGHT UM
module ModFaceValue

  use ModSize, ONLY: nI, nJ, nK
  use ModVarIndexes

  implicit none

  private ! except

  logical, public :: UseAccurateResChange = .false.
  logical, public :: UseTvdResChange      = .true.
  logical, public :: DoLimitMomentum      = .false.

  real,             public :: BetaLimiter = 1.0
  character(len=6), public :: TypeLimiter = 'minmod'

  public :: read_face_value_param, calc_face_value, correct_monotone_restrict, &
       calc_face_value_adjoint

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
  logical :: NameV
  integer :: nVarLimitRatio
  integer, allocatable :: iVarLimitRatio_I(:)

  ! Maximum length of the stencil in 1D
  integer,parameter:: MaxIJK = max(nI,nJ,nK)

  ! index ranges for optimized slope limiter calculations
  integer, parameter:: Lo2_=nVar+1, Hi2_=nVar+nVar, Lo3_=Hi2_+1, Hi3_=nVar+Hi2_

  ! local constants
  real, parameter :: cThird = 1./3., cTwoThird = 2./3., cSixth=1./6.

  ! variables used for TVD limiters
  real   :: dVarLimR_VI(1:nVar,0:MaxIJK+1) ! limited slope for right state
  real   :: dVarLimL_VI(1:nVar,0:MaxIJK+1) ! limited slope for left state
  real   :: Primitive_VI(1:nVar,-1:MaxIJK+2)
  logical:: IsTrueCell_I(-1:MaxIJK+2)

  !primitive variables (here/global for adjoint purposes)
  real, dimension(nVar,-1:nI+2,-1:nJ+2,-1:nK+2):: Primitive_VG
  real, dimension(nVar,-1:nI+2,-1:nJ+2,-1:nK+2):: AdjPrimitive_VG !ADJOINT SPECIFIC

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

    ! The subroutine calculates right and left face values.
    ! If DoResChangeOnly is true, only facevalues next to a coarser 
    ! neighbor block are calculated.
    ! To improve the code stability, ONLY those values are calculated, which
    ! are actually used by the calc_facefluxes.

    use ModMultiFluid, ONLY: nIonFluid
    use ModMain,     ONLY: nOrder, prolong_order, BlkTest, UseB0, &
         UseConstrainB, nIFace, nJFace, nKFace, &
         jMinFaceX, jMaxFaceX, kMinFaceX, kMaxFaceX,&
         iMinFaceY, iMaxFaceY, kMinFaceY, kMaxFaceY,&
         iMinFaceZ, iMaxFaceZ, jMinFaceZ, jMaxFaceZ

    use ModGeometry, ONLY : true_cell,body_BLK
    use ModNumConst
    use ModPhysics, ONLY: GammaWave
    use ModPhysics, ONLY: c2LIGHT,inv_c2LIGHT  !^CFG IF BORISCORR BEGIN
    use ModB0                                  !^CFG END BORISCORR
    use ModAdvance, ONLY: State_VGB, UseElectronPressure, UseWavePressure, &
         LeftState_VX,      &  ! Face Left  X
         RightState_VX,     &  ! Face Right X
         LeftState_VY,      &  ! Face Left  Y
         RightState_VY,     &  ! Face Right Y
         LeftState_VZ,      &  ! Face Left  Z
         RightState_VZ         ! Face Right Z

    use ModParallel, ONLY : &
         neiLEV,neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth

    implicit none

    logical, intent(in):: DoResChangeOnly
    integer, intent(in):: iBlock

    integer:: i,j,k,iSide,iFluid
    real:: RhoInv

    real:: RhoC2Inv, BxFull, ByFull, BzFull, B2Full,& !^CFG IF BORISCORR
         uBC2Inv,Ga2Boris                           !^CFG IF BORISCORR
    real:: B0_DG(3,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn)

    logical::DoTest,DoTestMe
    character(len=*), parameter :: NameSub = 'calc_face_value'
    !-------------------------------------------------------------------------

    if(iBlock==BLKtest .and. .not. DoResChangeOnly )then
       call set_oktest('calc_face_value', DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    UseLogLimiter   = nOrder == 2 .and. (UseLogRhoLimiter .or. UseLogPLimiter)
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

    UsePtotalLimiter = nOrder == 2 .and. nIonFluid == 1 .and. UsePtotalLtd

    if(.not.DoResChangeOnly & !In order not to call it twice
         .and.nOrder==2     & !Is not needed for nOrder=1
         .and. (UseAccurateResChange .or. UseTvdResChange)) &
         call correct_monotone_restrict(iBlock)

    ! first, calculate the CELL values for the variables to be limited
    ! for non-boris corrections they are: density, velocity, pressure
    ! for boris correction momentum is used instead of the velocity
    !
    if(DoLimitMomentum)then                     !^CFG IF BORISCORR BEGIN
       if(UseB0)then
          B0_DG=B0_DGB(:,:,:,:,iBlock)
       else
          B0_DG=0.00
       end if
       do k=kMinFaceX,kMaxFaceX
          do j=jMinFaceX,jMaxFaceX
             do i=1-nOrder,nI+nOrder
                call calc_primitives_boris    !needed for x-faces
             end do
          end do
       end do
       if(nJ > 1)then
          do k=kMinFaceY,kMaxFaceY; do i=iMinFaceY,iMaxFaceY
             do j=1-nOrder,jMinFaceX-1
                call calc_primitives_boris    ! additional calculations for 
             end do                           ! y -faces
             do j=jMaxFaceX+1,nJ+nOrder
                call calc_primitives_boris    ! additional calculations for 
             end do	                      ! y-faces
          end do; end do
       end if
       if(nK > 1)then
          do j=jMinFaceZ,jMaxFaceZ; do i=iMinFaceZ,iMaxFaceZ
             do k=1-nOrder,kMinFaceX-1
                call calc_primitives_boris    ! additional calculations for 
             end do                           ! z-faces
             do k=kMaxFaceX+1,nK+nOrder
                call calc_primitives_boris    ! additional calculations for 
             end do	                      ! z-faces
          end do; end do
       end if
    else                                         !^CFG END BORISCORR
       if(UseAccurateResChange)then
          do k=-1,nK+2; do j=-1, nJ+2; do i=-1, nI+2
             call calc_primitives_MHD         !needed for x-faces
          end do; end do; end do
       else
          do k=kMinFaceX,kMaxFaceX
             do j=jMinFaceX,jMaxFaceX
                do i=1-nOrder,nI+nOrder
                   call calc_primitives_MHD   !needed for x-faces
                end do
             end do
          end do
          if(nJ > 1)then
             do k=kMinFaceY,kMaxFaceY; do i=iMinFaceY,iMaxFaceY
                do j=1-nOrder,jMinFaceX-1
                   call calc_primitives_MHD   ! additional calculations for 
                end do                        ! y -faces
                do j=jMaxFaceX+1,nJ+nOrder
                   call calc_primitives_MHD   ! additional calculations for 
                end do	                      ! y-faces
             end do; end do
          end if
          if(nK > 1)then
             do j=jMinFaceZ,jMaxFaceZ; do i=iMinFaceZ,iMaxFaceZ
                do k=1-nOrder,kMinFaceX-1
                   call calc_primitives_MHD   ! additional calculations for 
                end do                        ! z-faces
                do k=kMaxFaceX+1,nK+nOrder
                   call calc_primitives_MHD   ! additional calculations for 
                end do	                      ! z-faces
             end do; end do
          end if
       end if
    end if                                       !^CFG IF BORISCORR

    !\
    ! Now the first or second order face values are calcuted
    !/
    select case(nOrder)
    case(1)
       !\
       ! First order reconstruction
       !/
       if (.not.DoResChangeOnly) then
          call get_faceX_first(&
               1,nIFace,jMinFaceX,jMaxFaceX,kMinFaceX,kMaxFaceX)
          if(nJ > 1) call get_faceY_first(&
               iMinFaceY,iMaxFaceY,1,nJFace,kMinFaceY,kMaxFaceY)
          if(nK > 1) call get_faceZ_first(&
               iMinFaceZ,iMaxFaceZ,jMinFaceZ,jMaxFaceZ,1,nKFace)
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
    case(2)
       ! For second order scheme (nOrder==2)   
       ! use second order limited reconstruction.
       ! When prolong_order==1 we can use first order reconstruction 
       ! at resolution changes.
       ! However, constrained transport requires facevalues !^CFG IF CONSTRAINB
       ! to be independent of the resolution changes.       !^CFG IF CONSTRAINB

       if (.not.DoResChangeOnly)then
          call get_faceX_second(&
               1,nIFace,jMinFaceX,jMaxFaceX,kMinFaceX,kMaxFaceX)
          if(nJ > 1) call get_faceY_second(&
               iMinFaceY,iMaxFaceY,1,nJFace,kMinFaceY,kMaxFaceY)
          if(nK > 1) call get_faceZ_second(&
               iMinFaceZ,iMaxFaceZ,jMinFaceZ,jMaxFaceZ,1,nKFace)
       end if

       if(prolong_order==1 &
            .and..not.UseConstrainB & !^CFG IF CONSTRAINB
            )then

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
             do iSide=east_,top_
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
             call logfaceX_to_faceX(1,nIFace,jMinFaceX,jMaxFaceX, &
                  kMinFaceX,kMaxFaceX)
             if(nJ > 1) call logfaceY_to_faceY(iMinFaceY,iMaxFaceY,1,nJFace, &
                  kMinFaceY,kMaxFaceY)
             if(nK > 1) call logfaceZ_to_faceZ(iMinFaceZ,iMaxFaceZ, &
                  jMinFaceZ,jMaxFaceZ,1,nKFace)
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
             call ratio_to_scalar_faceX(1,nIFace,jMinFaceX,jMaxFaceX, &
                  kMinFaceX,kMaxFaceX)
             if(nJ > 1) call ratio_to_scalar_faceY(iMinFaceY,iMaxFaceY,1,nJFace, &
                  kMinFaceY,kMaxFaceY)
             if(nK > 1) call ratio_to_scalar_faceZ(iMinFaceZ,iMaxFaceZ, &
                  jMinFaceZ,jMaxFaceZ,1,nKFace)
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
             call ptotal_to_p_faceX(1,nIFace,jMinFaceX,jMaxFaceX, &
                  kMinFaceX,kMaxFaceX)
             if(nJ > 1) call ptotal_to_p_faceY(iMinFaceY,iMaxFaceY,1,nJFace, &
                  kMinFaceY,kMaxFaceY)
             if(nK > 1) call ptotal_to_p_faceZ(iMinFaceZ,iMaxFaceZ, &
                  jMinFaceZ,jMaxFaceZ,1,nKFace)
          end if
       end if

    end select  !end second order

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
    !^CFG IF BORISCORR BEGIN
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
      RhoC2Inv  =cOne/(Primitive_VG(rho_,i,j,k)*c2LIGHT)
      uBC2Inv= (Primitive_VG(rhoUx_,i,j,k)*BxFull + &
           Primitive_VG(rhoUy_,i,j,k)*ByFull + &
           Primitive_VG(rhoUz_,i,j,k)*BzFull)*RhoC2Inv
      Ga2Boris=cOne+B2Full*RhoC2Inv

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
    !^CFG END BORISCORR
    subroutine get_faceX_first(iMin,iMax,jMin,jMax,kMin,kMax)

      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      !-----------------------------------------------------------------------
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         LeftState_VX(:,i,j,k)=Primitive_VG(:,i-1,j,k)
         RightState_VX(:,i,j,k)=Primitive_VG(:,i,j,k)              
      end do; end do; end do

      !^CFG IF BORISCORR BEGIN
      if(DoLimitMomentum)call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 
      !^CFG END BORISCORR

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

      !^CFG IF BORISCORR BEGIN
      if(DoLimitMomentum) call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)
      !^CFG END BORISCORR

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

      !^CFG IF BORISCORR BEGIN
      if(DoLimitMomentum)call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 
      !^CFG END BORISCORR

      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceZ(&
           iMin,iMax,jMin,jMax,kMin,kMax)

    end subroutine get_faceZ_first
    !==========================================================================
    subroutine calc_primitives_MHD
      use ModMultiFluid
      Primitive_VG(:,i,j,k) = State_VGB(1:nVar,i,j,k,iBlock)
      ! TEMPORARY
      !if (iBlock == 1) then
      !   write(*,*)'(i,j,k): state = ',i,j,k,State_VGB(1,i,j,k,iBlock)
      !end if
      RhoInv=cOne/Primitive_VG(Rho_,i,j,k)
      Primitive_VG(Ux_:Uz_,i,j,k)=RhoInv*Primitive_VG(RhoUx_:RhoUz_,i,j,k)
      do iFluid = 2, nFluid
         call select_fluid
         RhoInv=cOne/Primitive_VG(iRho,i,j,k)
         Primitive_VG(iUx:iUz,i,j,k)=RhoInv*Primitive_VG(iRhoUx:iRhoUz,i,j,k)
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
    !^CFG IF BORISCORR BEGIN
    subroutine BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

      ! Convert face centered Boris momenta to MHD velocities

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      ! U_Boris=rhoU_Boris/rho
      ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

         ! Left face values
         RhoInv=cOne/LeftState_VX(rho_,i,j,k)
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
         Ga2Boris=cOne/(cOne+ B2Full*RhoC2Inv)

         LeftState_VX(Ux_,i,j,k) = &
              Ga2Boris * (LeftState_VX(Ux_,i,j,k)+uBC2Inv*BxFull)
         LeftState_VX(Uy_,i,j,k) = &
              Ga2Boris * (LeftState_VX(Uy_,i,j,k)+uBC2Inv*ByFull)
         LeftState_VX(Uz_,i,j,k) = &
              Ga2Boris * (LeftState_VX(Uz_,i,j,k)+uBC2Inv*BzFull)

         ! Right face values
         RhoInv=cOne/RightState_VX(rho_,i,j,k)
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
         Ga2Boris=cOne/(cOne+B2Full*RhoC2Inv)

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
         RhoInv=cOne/LeftState_VY(rho_,i,j,k)
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
         Ga2Boris=cOne/(cOne+ B2Full*RhoC2Inv)

         LeftState_VY(Ux_,i,j,k) = &
              Ga2Boris * (LeftState_VY(Ux_,i,j,k)+uBC2Inv*BxFull)
         LeftState_VY(Uy_,i,j,k) = &
              Ga2Boris * (LeftState_VY(Uy_,i,j,k)+uBC2Inv*ByFull)
         LeftState_VY(Uz_,i,j,k) = &
              Ga2Boris * (LeftState_VY(Uz_,i,j,k)+uBC2Inv*BzFull)

         ! Right face values
         RhoInv=cOne/RightState_VY(rho_,i,j,k) 
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
         Ga2Boris=cOne/(cOne + B2Full*RhoC2Inv)

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
         RhoInv=cOne/LeftState_VZ(rho_,i,j,k)
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
         Ga2Boris=cOne/(cOne + B2Full*RhoC2Inv)

         LeftState_VZ(Ux_,i,j,k) = &
              Ga2Boris * (LeftState_VZ(Ux_,i,j,k)+uBC2Inv*BxFull)
         LeftState_VZ(Uy_,i,j,k) = &
              Ga2Boris * (LeftState_VZ(Uy_,i,j,k)+uBC2Inv*ByFull)
         LeftState_VZ(Uz_,i,j,k) = &
              Ga2Boris * (LeftState_VZ(Uz_,i,j,k)+uBC2Inv*BzFull)

         ! Right face values
         RhoInv=cOne/RightState_VZ(rho_,i,j,k)
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
         Ga2Boris=cOne/(cOne + B2Full*RhoC2Inv)

         RightState_VZ(Ux_,i,j,k) = &
              Ga2Boris * (RightState_VZ(Ux_,i,j,k)+uBC2Inv*BxFull)
         RightState_VZ(Uy_,i,j,k) = &
              Ga2Boris * (RightState_VZ(Uy_,i,j,k)+uBC2Inv*ByFull)
         RightState_VZ(Uz_,i,j,k) = &
              Ga2Boris * (RightState_VZ(Uz_,i,j,k)+uBC2Inv*BzFull)
      end do; end do; end do

    end subroutine BorisFaceZtoMHD
    !^CFG END BORISCORR
    !==========================================================================
    subroutine get_face_accurate3d(iSideIn)
      integer, intent(in):: iSideIn

      select case(iSideIn)
      case(east_)
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
      case(west_)
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
      case(south_)
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
      case(north_)
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
      case(bot_)
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
      case(top_)
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
      case(south_)
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
      case(north_)
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
      case(east_)
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
      case(west_)
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
      case(south_)
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
      case(north_)
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
      case(bot_)
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
      case(top_)
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
         if(body_BLK(iBlock))then
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
      !^CFG IF BORISCORR BEGIN
      if(DoLimitMomentum) call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 
      !^CFG END BORISCORR
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
         if(body_BLK(iBlock))then
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
      !^CFG IF BORISCORR BEGIN
      if(DoLimitMomentum) call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 
      !^CFG END BORISCORR
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
         if(body_BLK(iBlock))then
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
      !^CFG IF BORISCORR BEGIN
      if(DoLimitMomentum) call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)
      !^CFG END BORISCORR
    end subroutine get_faceZ_second
  end subroutine calc_face_value


  !ADJOINT SPECIFIC BEGIN
  !===========================================================================
  subroutine calc_face_value_adjoint(DoResChangeOnly, iBlock)

    ! Adjoint version

    use ModMain
    use ModVarIndexes
    use ModGeometry, ONLY : true_cell,body_BLK
    use ModNumConst
    use ModPhysics, ONLY: c2LIGHT,inv_c2LIGHT  !^CFG IF BORISCORR BEGIN
    use ModB0                                  !^CFG END BORISCORR
    use ModAdvance, ONLY: State_VGB,&
         LeftState_VX,      &  ! Face Left  X
         RightState_VX,     &  ! Face Right X
         LeftState_VY,      &  ! Face Left  Y
         RightState_VY,     &  ! Face Right Y
         LeftState_VZ,      &  ! Face Left  Z
         RightState_VZ         ! Face Right Z

    use ModParallel, ONLY : &
         neiLEV,neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth

    use ModAdjoint

    implicit none

    logical, intent(in):: DoResChangeOnly
    integer, intent(in):: iBlock


    integer:: i,j,k,iSide,iFluid
    real:: RhoInv
    real:: AdjRhoInv

    real:: RhoC2Inv, BxFull, ByFull, BzFull, B2Full,& !^CFG IF BORISCORR
         uBC2Inv,Ga2Boris                           !^CFG IF BORISCORR
    real:: B0_DG(3,1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn)
    !primitive variables
    real, dimension(nVar,-1:nI+2,-1:nJ+2,-1:nK+2):: Primitive_VG
    !primitive adjoint variables
    real, dimension(nVar,-1:nI+2,-1:nJ+2,-1:nK+2):: AdjPrimitive_VG

    logical::DoTest,DoTestMe
    character(len=*), parameter :: NameSub = 'calc_face_value_adjoint'
    !-------------------------------------------------------------------------

    UseLogLimiter   = nOrder == 2 .and. (UseLogRhoLimiter .or. UseLogPLimiter)
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

    if(.not.DoResChangeOnly & !In order not to call it twice
         .and.nOrder==2     & !Is not needed for nOrder=1
         .and. (UseAccurateResChange .or. UseTvdResChange)) &
         call stop_mpi(NameSub // ' Not yet supported')


    ! State Primitive variables are already filled in from forward call
    ! LeftState, RightState are also valid, as calculated by the end
    !  of calc_face_value

    ! Initialize adjoint variables
    AdjPrimitive_VG = 0.0

    !\
    ! Now the first or second order face values are calculated
    ! Adjoint: this is first
    !/
    select case(nOrder)
    case(1)
       !\
       ! First order reconstruction
       !/
       if (.not.DoResChangeOnly) then
          call get_faceX_first_adjoint(&
               1,nIFace,jMinFaceX,jMaxFaceX,kMinFaceX,kMaxFaceX)
          if(nJ > 1) call get_faceY_first_adjoint(&
               iMinFaceY,iMaxFaceY,1,nJFace,kMinFaceY,kMaxFaceY)
          if(nK > 1) call get_faceZ_first_adjoint(&
               iMinFaceZ,iMaxFaceZ,jMinFaceZ,jMaxFaceZ,1,nKFace)
       else
          if(neiLeast(iBlock)==+1)&
               call get_faceX_first_adjoint(1,1,1,nJ,1,nK)
          if(neiLwest(iBlock)==+1)&
               call get_faceX_first_adjoint(nIFace,nIFace,1,nJ,1,nK)
          if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
               call get_faceY_first_adjoint(1,nI,1,1,1,nK)
          if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
               call get_faceY_first_adjoint(1,nI,nJFace,nJFace,1,nK)
          if(nK > 1 .and. neiLbot(iBlock)==+1) &
               call get_faceZ_first_adjoint(1,nI,1,nJ,1,1)
          if(nK > 1 .and. neiLtop(iBlock)==+1) &
               call get_faceZ_first_adjoint(1,nI,1,nJ,nKFace,nKFace)
       end if
    case(2)
       ! For second order scheme (nOrder==2)   
       ! use second order limited reconstruction.
       ! When prolong_order==1 we can use first order reconstruction 
       ! at resolution changes.
       ! However, constrained transport requires facevalues !^CFG IF CONSTRAINB
       ! to be independent of the resolution changes.       !^CFG IF CONSTRAINB

       
       if(UseScalarToRhoRatioLtd)then
          if(DoResChangeOnly)then
             if(neiLeast(iBlock)==+1) &
                  call ratio_to_scalar_faceX_adjoint(1,1,1,nJ,1,nK)
             if(neiLwest(iBlock)==+1) &
                  call ratio_to_scalar_faceX_adjoint(nIFace,nIFace,1,nJ,1,nK)
             if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
                  call ratio_to_scalar_faceY_adjoint(1,nI,1,1,1,nK)
             if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
                  call ratio_to_scalar_faceY_adjoint(1,nI,nJFace,nJFace,1,nK)
             if(nK > 1 .and. neiLbot(iBlock)==+1) &
                  call ratio_to_scalar_faceZ_adjoint(1,nI,1,nJ,1,1)
             if(nK > 1 .and. neiLtop(iBlock)==+1) &
                  call ratio_to_scalar_faceZ_adjoint(1,nI,1,nJ,nKFace,nKFace)
          else
             call ratio_to_scalar_faceX_adjoint(1,nIFace,jMinFaceX,jMaxFaceX, &
                  kMinFaceX,kMaxFaceX)
             if(nJ > 1) call ratio_to_scalar_faceY_adjoint(iMinFaceY,iMaxFaceY,1,nJFace, &
                  kMinFaceY,kMaxFaceY)
             if(nK > 1) call ratio_to_scalar_faceZ_adjoint(iMinFaceZ,iMaxFaceZ, &
                  jMinFaceZ,jMaxFaceZ,1,nKFace)
          end if
       end if


       if(UseLogLimiter .and. .not.DoLimitMomentum)then
          if(DoResChangeOnly)then
             call stop_mpi(NameSub // ' Not yet supported')
          else
             call logfaceX_to_faceX_adjoint(1,nIFace,jMinFaceX,jMaxFaceX, &
                  kMinFaceX,kMaxFaceX)
             if(nJ > 1) call logfaceY_to_faceY_adjoint(iMinFaceY,iMaxFaceY,1,nJFace, &
                  kMinFaceY,kMaxFaceY)
             if(nK > 1) call logfaceZ_to_faceZ_adjoint(iMinFaceZ,iMaxFaceZ, &
                  jMinFaceZ,jMaxFaceZ,1,nKFace)
          end if
       end if

       if(prolong_order==1 &
            .and..not.UseConstrainB & !^CFG IF CONSTRAINB
            )then
          call stop_mpi(NameSub // ' Not yet supported')
       else if(DoResChangeOnly) then
          ! Second order face values at resolution changes
          if(neiLeast(iBlock)==+1)&
               call get_faceX_second_adjoint(1,1,1,nJ,1,nK)
          if(neiLwest(iBlock)==+1)&
               call get_faceX_second_adjoint(nIFace,nIFace,1,nJ,1,nK)
          if(nJ > 1 .and. neiLsouth(iBlock)==+1) &
               call get_faceY_second_adjoint(1,nI,1,1,1,nK)
          if(nJ > 1 .and. neiLnorth(iBlock)==+1) &
               call get_faceY_second_adjoint(1,nI,nJFace,nJFace,1,nK)
          if(nK > 1 .and. neiLbot(iBlock)==+1) &
               call get_faceZ_second_adjoint(1,nI,1,nJ,1,1)
          if(nK > 1 .and. neiLtop(iBlock)==+1) &
               call get_faceZ_second_adjoint(1,nI,1,nJ,nKFace,nKFace)
       endif

       if (.not.DoResChangeOnly)then
          call get_faceX_second_adjoint(&
               1,nIFace,jMinFaceX,jMaxFaceX,kMinFaceX,kMaxFaceX)
          if(nJ > 1) call get_faceY_second_adjoint(&
               iMinFaceY,iMaxFaceY,1,nJFace,kMinFaceY,kMaxFaceY)
          if(nK > 1) call get_faceZ_second_adjoint(&
               iMinFaceZ,iMaxFaceZ,jMinFaceZ,jMaxFaceZ,1,nKFace)
       end if
    end select  !end second order


    ! first, calculate the CELL values for the variables to be limited
    ! for non-boris corrections they are: density, velocity, pressure
    ! for boris correction momentum is used instead of the velocity
    ! Adjoint: this is second
    if(DoLimitMomentum)then                     !^CFG IF BORISCORR BEGIN
       call stop_mpi(NameSub // ' Not yet supported')
    else                                         !^CFG END BORISCORR
       if(UseAccurateResChange)then
          call stop_mpi(NameSub // ' Not yet supported')
       else
          do k=kMinFaceX,kMaxFaceX
             do j=jMinFaceX,jMaxFaceX
                do i=1-nOrder,nI+nOrder
                   call calc_primitives_MHD_adjoint   !needed for x-faces
                end do
             end do
          end do
          if(nJ > 1)then
             do k=kMinFaceY,kMaxFaceY; do i=iMinFaceY,iMaxFaceY
                do j=1-nOrder,jMinFaceX-1
                   call calc_primitives_MHD_adjoint   ! additional calculations for 
                end do                        ! y -faces
                do j=jMaxFaceX+1,nJ+nOrder
                   call calc_primitives_MHD_adjoint   ! additional calculations for 
                end do	                      ! y-faces
             end do; end do
          end if
          if(nK > 1)then
             do j=jMinFaceZ,jMaxFaceZ; do i=iMinFaceZ,iMaxFaceZ
                do k=1-nOrder,kMinFaceX-1
                   call calc_primitives_MHD_adjoint   ! additional calculations for 
                end do                        ! z-faces
                do k=kMaxFaceX+1,nK+nOrder
                   call calc_primitives_MHD_adjoint   ! additional calculations for 
                end do	                      ! z-faces
             end do; end do
          end if
       end if
    end if                                       !^CFG IF BORISCORR

  contains

    ! NOTE: in the adjoint versions of the below routines,
    ! Left/RightState_V* is assumed to be post-forward routine.  For
    ! example, in the log limiter case, the exponential has already
    ! been taken.  These routines perform the reverse of the forward step
    ! on the state to prepare the state for use in get_face*_*_adjoint

    !==========================================================================
    subroutine logfaceX_to_faceX_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      do iVar=1,nVar
         if(.not.UseLogLimiter_V(iVar))CYCLE
         ! Reverse operation on state
         LeftState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              log(LeftState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         RightState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              log(RightState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         ! adjoint operation
         AdjLeftState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              AdjLeftState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax)*&
              (1.0+exp(LeftState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax)))
         AdjRightState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              AdjRightState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax)*&
              (1.0+exp(RightState_VX(iVar,iMin:iMax,jMin:jMax,kMin:kMax)))
      end do

    end subroutine logfaceX_to_faceX_adjoint
    !==========================================================================
    subroutine logfaceY_to_faceY_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      do iVar=1,nVar
         if(.not.UseLogLimiter_V(iVar))CYCLE
         ! Reverse operation on state
         LeftState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              log(LeftState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         RightState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              log(RightState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         ! adjoint operation
         AdjLeftState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              AdjLeftState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax)*&
              (1.0+exp(LeftState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax)))
         AdjRightState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              AdjRightState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax)*&
              (1.0+exp(RightState_VY(iVar,iMin:iMax,jMin:jMax,kMin:kMax)))
      end do

    end subroutine logfaceY_to_faceY_adjoint
    !==========================================================================
    subroutine logfaceZ_to_faceZ_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      !------------------------------------------------------------------------
      do iVar=1,nVar
         if(.not.UseLogLimiter_V(iVar))CYCLE
         ! Reverse operation on state
         LeftState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              log(LeftState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         RightState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              log(RightState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax))
         ! adjoint operation
         AdjLeftState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              AdjLeftState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax)*&
              (1.0+exp(LeftState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax)))
         AdjRightState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax) = &
              AdjRightState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax)*&
              (1.0+exp(RightState_VZ(iVar,iMin:iMax,jMin:jMax,kMin:kMax)))
      end do

    end subroutine logfaceZ_to_faceZ_adjoint
    !==========================================================================
    subroutine ratio_to_scalar_faceX_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer :: i, j, k
      !------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         ! Reverse operation on state
         LeftState_VX(iVarLimitRatio_I,i,j,k) = &
              LeftState_VX(iVarLimitRatio_I,i,j,k)/LeftState_VX(Rho_,i,j,k)
         RightState_VX(iVarLimitRatio_I,i,j,k) = &
              RightState_VX(iVarLimitRatio_I,i,j,k)/RightState_VX(Rho_,i,j,k)
         ! adjoint operation
         AdjLeftState_VX(Rho_,i,j,k) = AdjLeftState_VX(Rho_,i,j,k) + &
              dot_product(LeftState_VX(iVarLimitRatio_I,i,j,k),&
              AdjLeftState_VX(iVarLimitRatio_I,i,j,k))
         AdjLeftState_VX(iVarLimitRatio_I,i,j,k) = &
              AdjLeftState_VX(iVarLimitRatio_I,i,j,k)*(1.0+Leftstate_VX(Rho_,i,j,k))
         AdjRightState_VX(Rho_,i,j,k) = AdjRightState_VX(Rho_,i,j,k) + &
              dot_product(RightState_VX(iVarLimitRatio_I,i,j,k),&
              AdjRightState_VX(iVarLimitRatio_I,i,j,k))
         AdjRightState_VX(iVarLimitRatio_I,i,j,k) = &
              AdjRightState_VX(iVarLimitRatio_I,i,j,k)*(1.0+Rightstate_VX(Rho_,i,j,k))
      end do; end do; end do

    end subroutine ratio_to_scalar_faceX_adjoint
    !==========================================================================
    subroutine ratio_to_scalar_faceY_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer :: i, j, k
      !------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         ! Reverse operation on state
         LeftState_VY(iVarLimitRatio_I,i,j,k) = &
              LeftState_VY(iVarLimitRatio_I,i,j,k)/LeftState_VY(Rho_,i,j,k)
         RightState_VY(iVarLimitRatio_I,i,j,k) = &
              RightState_VY(iVarLimitRatio_I,i,j,k)/RightState_VY(Rho_,i,j,k)
         ! adjoint operation
         AdjLeftState_VY(Rho_,i,j,k) = AdjLeftState_VY(Rho_,i,j,k) + &
              dot_product(LeftState_VY(iVarLimitRatio_I,i,j,k),&
              AdjLeftState_VY(iVarLimitRatio_I,i,j,k))
         AdjLeftState_VY(iVarLimitRatio_I,i,j,k) = &
              AdjLeftState_VY(iVarLimitRatio_I,i,j,k)*(1.0+Leftstate_VY(Rho_,i,j,k))
         AdjRightState_VY(Rho_,i,j,k) = AdjRightState_VY(Rho_,i,j,k) + &
              dot_product(RightState_VY(iVarLimitRatio_I,i,j,k),&
              AdjRightState_VY(iVarLimitRatio_I,i,j,k))
         AdjRightState_VY(iVarLimitRatio_I,i,j,k) = &
              AdjRightState_VY(iVarLimitRatio_I,i,j,k)*(1.0+Rightstate_VY(Rho_,i,j,k))
      end do; end do; end do

    end subroutine ratio_to_scalar_faceY_adjoint
    !==========================================================================
    subroutine ratio_to_scalar_faceZ_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)

      integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer :: i, j, k
      !------------------------------------------------------------------------
      do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
         ! Reverse operation on state
         LeftState_VZ(iVarLimitRatio_I,i,j,k) = &
              LeftState_VZ(iVarLimitRatio_I,i,j,k)/LeftState_VZ(Rho_,i,j,k)
         RightState_VZ(iVarLimitRatio_I,i,j,k) = &
              RightState_VZ(iVarLimitRatio_I,i,j,k)/RightState_VZ(Rho_,i,j,k)
         ! adjoint operation
         AdjLeftState_VZ(Rho_,i,j,k) = AdjLeftState_VZ(Rho_,i,j,k) + &
              dot_product(LeftState_VZ(iVarLimitRatio_I,i,j,k),&
              AdjLeftState_VZ(iVarLimitRatio_I,i,j,k))
         AdjLeftState_VZ(iVarLimitRatio_I,i,j,k) = &
              AdjLeftState_VZ(iVarLimitRatio_I,i,j,k)*(1.0+Leftstate_VZ(Rho_,i,j,k))
         AdjRightState_VZ(Rho_,i,j,k) = AdjRightState_VZ(Rho_,i,j,k) + &
              dot_product(RightState_VZ(iVarLimitRatio_I,i,j,k),&
              AdjRightState_VZ(iVarLimitRatio_I,i,j,k))
         AdjRightState_VZ(iVarLimitRatio_I,i,j,k) = &
              AdjRightState_VZ(iVarLimitRatio_I,i,j,k)*(1.0+Rightstate_VZ(Rho_,i,j,k))
      end do; end do; end do

    end subroutine ratio_to_scalar_faceZ_adjoint

    !==========================================================================
    subroutine calc_primitives_MHD
      use ModMultiFluid
      Primitive_VG(:,i,j,k) = State_VGB(1:nVar,i,j,k,iBlock)
      ! TEMPORARY
      !if (iBlock == 1) then
      !   write(*,*)'(i,j,k): state = ',i,j,k,State_VGB(1,i,j,k,iBlock)
      !end if
      RhoInv=cOne/Primitive_VG(Rho_,i,j,k)
      Primitive_VG(Ux_:Uz_,i,j,k)=RhoInv*Primitive_VG(RhoUx_:RhoUz_,i,j,k)
      do iFluid = 2, nFluid
         call select_fluid
         RhoInv=cOne/Primitive_VG(iRho,i,j,k)
         Primitive_VG(iUx:iUz,i,j,k)=RhoInv*Primitive_VG(iRhoUx:iRhoUz,i,j,k)
      end do
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
    subroutine calc_primitives_MHD_adjoint

      use ModMultiFluid

      real, dimension(nVar,-1:nI+2,-1:nJ+2,-1:nK+2):: Primitive_VG_A

      Primitive_VG(:,i,j,k) = State_VGB(1:nVar,i,j,k,iBlock)
      RhoInv=cOne/Primitive_VG(Rho_,i,j,k)
      Primitive_VG(Ux_:Uz_,i,j,k)=RhoInv*Primitive_VG(RhoUx_:RhoUz_,i,j,k)

      Primitive_VG_A = Primitive_VG   ! store in A

      do iFluid = 2, nFluid
         call select_fluid
         RhoInv=cOne/Primitive_VG(iRho,i,j,k)
         Primitive_VG(iUx:iUz,i,j,k)=RhoInv*Primitive_VG(iRhoUx:iRhoUz,i,j,k)
      end do
      if(UseScalarToRhoRatioLtd) Primitive_VG(iVarLimitRatio_I,i,j,k) = &
           RhoInv*Primitive_VG(iVarLimitRatio_I,i,j,k)

      if(UseLogLimiter)then
         do iVar=1,nVar
            if(UseLogLimiter_V(iVar)) &
                 Primitive_VG(iVar,i,j,k) = log(Primitive_VG(iVar,i,j,k))
         end do
      end if

      ! Adjoint section
      AdjRhoInv = 0.0

      if(UseLogLimiter) call stop_mpi(NameSub // ' Not yet supported')

      if(UseScalarToRhoRatioLtd)then
         Primitive_VG = Primitive_VG_A  ! restore from A
         AdjRhoInv = AdjRhoInv + dot_product(Primitive_VG(iVarLimitRatio_I,i,j,k),&
              AdjPrimitive_VG(iVarLimitRatio_I,i,j,k))
         AdjPrimitive_VG(iVarLimitRatio_I,i,j,k) = &
              AdjPrimitive_VG(iVarLimitRatio_I,i,j,k) + &
              RhoInv*AdjPrimitive_VG(iVarLimitRatio_I,i,j,k)
      end if

      ! mutifluid not yet supported
      if (iFluid > 1)call stop_mpi(NameSub // ' Not yet supported')

      ! reset Primitive_VG
      Primitive_VG(:,i,j,k) = State_VGB(1:nVar,i,j,k,iBlock)
      RhoInv=cOne/Primitive_VG(Rho_,i,j,k)
      !Primitive_VG(Ux_:Uz_,i,j,k)=RhoInv*Primitive_VG(RhoUx_:RhoUz_,i,j,k)  
      AdjRhoInv = AdjRhoInv + dot_product(Primitive_VG(RhoUx_:RhoUz_,i,j,k),&
           AdjPrimitive_VG(RhoUx_:RhoUz_,i,j,k))
      AdjPrimitive_VG(RhoUx_:RhoUz_,i,j,k) = AdjPrimitive_VG(RhoUx_:RhoUz_,i,j,k) &
           + RhoInv*AdjPrimitive_VG(Ux_:Uz_,i,j,k)


      AdjPrimitive_VG(Rho_,i,j,k) =  AdjPrimitive_VG(Rho_,i,j,k) &
           - RhoInv**2 * AdjRhoInv

      ! Update AdjointVGB
      Adjoint_VGB(1:nVar,i,j,k,iBlock) = Adjoint_VGB(1:nVar,i,j,k,iBlock) &
           + AdjPrimitive_VG(:,i,j,k)

    end subroutine calc_primitives_MHD_adjoint
    !=========================================================================
    subroutine get_faceX_first_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceX_adjoint(&
           iMin,iMax,jMin,jMax,kMin,kMax)
      !^CFG IF BORISCORR BEGIN
      if(DoLimitMomentum)call stop_mpi(NameSub // ' Not yet supported')
      !^CFG END BORISCORR
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         AdjPrimitive_VG(:,i-1,j,k) = AdjPrimitive_VG(:,i-1,j,k) + AdjLeftState_VX(:,i,j,k)
         AdjPrimitive_VG(:,i  ,j,k) = AdjPrimitive_VG(:,i  ,j,k) + AdjRightState_VX(:,i,j,k)
      end do; end do; end do
    end subroutine get_faceX_first_adjoint
    !========================================================================
    subroutine get_faceY_first_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceY_adjoint(&
           iMin,iMax,jMin,jMax,kMin,kMax)
      !^CFG IF BORISCORR BEGIN
      if(DoLimitMomentum) call stop_mpi(NameSub // ' Not yet supported')
      !^CFG END BORISCORR
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         AdjPrimitive_VG(:,i,j-1,k) = AdjPrimitive_VG(:,i,j-1,k) + AdjLeftState_VY(:,i,j,k)
         AdjPrimitive_VG(:,i,j  ,k) = AdjPrimitive_VG(:,i,j  ,k) + AdjRightState_VY(:,i,j,k)
      end do; end do; end do
    end subroutine get_faceY_first_adjoint
    !========================================================================
    subroutine get_faceZ_first_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      if(UseScalarToRhoRatioLtd)call ratio_to_scalar_faceZ_adjoint(&
           iMin,iMax,jMin,jMax,kMin,kMax)
      !^CFG IF BORISCORR BEGIN
      if(DoLimitMomentum)call stop_mpi(NameSub // ' Not yet supported')
      !^CFG END BORISCORR
      do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
         AdjPrimitive_VG(:,i,j,k-1) = AdjPrimitive_VG(:,i,j,k-1) + AdjLeftState_VZ(:,i,j,k)
         AdjPrimitive_VG(:,i,j,k  ) = AdjPrimitive_VG(:,i,j,k  ) + AdjRightState_Vz(:,i,j,k)
      end do; end do; end do
    end subroutine get_faceZ_first_adjoint
    !=========================================================================
    subroutine get_faceX_second_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      integer::i1, iMinSharp, iMaxSharp
      !----------------------------------------------------------------------
      call stop_mpi(NameSub // ' Not yet supported')
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
         if(body_BLK(iBlock))then
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
      !^CFG IF BORISCORR BEGIN
      !if(DoLimitMomentum) call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 
      !^CFG END BORISCORR
    end subroutine get_faceX_second_adjoint
    !==========================================================================
    subroutine get_faceY_second_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
      integer::j1, jMinSharp, jMaxSharp
      !----------------------------------------------------------------------
      call stop_mpi(NameSub // ' Not yet supported')
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
         if(body_BLK(iBlock))then
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
      !^CFG IF BORISCORR BEGIN
      !if(DoLimitMomentum) call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) 
      !^CFG END BORISCORR
    end subroutine get_faceY_second_adjoint
    !==========================================================================
    subroutine get_faceZ_second_adjoint(iMin,iMax,jMin,jMax,kMin,kMax)
      integer,intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
      integer::k1, kMinSharp, kMaxSharp
      !----------------------------------------------------------------------
      call stop_mpi(NameSub // ' Not yet supported')
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
         if(body_BLK(iBlock))then
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
      !^CFG IF BORISCORR BEGIN
      !if(DoLimitMomentum) call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)
      !^CFG END BORISCORR
    end subroutine get_faceZ_second_adjoint
  end subroutine calc_face_value_adjoint
  !ADJOINT SPECIFIC END



  !============================================================================
  ! Limiters:
  ! mimod limiter:
  !                slim = minmod(s1,s2)
  !
  ! generalized MC (monotonized central) limiter:
  !                slim = minmod(beta*s1, beta*s2, (s1+s2)/2)
  !
  ! beta-limiter   slim = maxmod(minmod(beta*s1,s2), minmod(beta*s2,s1))
  !
  !                (see C.Hirsch, Numerical Computation of
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
  !       left_face  = central_value - limited_slope
  !       right_face = central_value + limited_slope
  !===========================================================================
  subroutine limiter_body(lMin,lMax,Beta)

    integer, intent(in):: lMin,lMax
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
       ! Calculate right most unlimited slope
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
       ! Calculate right most unlimited slope
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
    use ModVarIndexes, ONLY: DefaultState_V, nVar, nFluid, &
         iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I
    use ModAdvance,    ONLY: State_VGB
    use ModAMR,        ONLY: unusedBlock_BP
    use ModParallel,   ONLY: neiLEV, &
         neiLtop, neiLbot, neiLeast, neiLwest, neiLnorth, neiLsouth, &
         neiBtop, neiBbot, neiBeast, neiBwest, neiBnorth, neiBsouth, &
         neiPtop, neiPbot, neiPeast, neiPwest, neiPnorth, neiPsouth
    use ModNumConst, ONLY: cTiny

    ! For debugging
    use ModProcMH, ONLY: iProc
    use ModMain, ONLY: VarTest, ProcTest, BlkTest, iTest, jTest, kTest

    implicit none

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

    if(nFluid > 1 .and. .not.DoLimitMomentum)then
       ! Convert momenta to velocities (that will be limited)
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          State_VGB(iRhoUx_I,i,j,k,iBLK)=State_VGB(iRhoUx_I,i,j,k,iBLK) &
               / State_VGB(iRho_I,i,j,k,iBLK)
          State_VGB(iRhoUy_I,i,j,k,iBLK)=State_VGB(iRhoUy_I,i,j,k,iBLK) &
               / State_VGB(iRho_I,i,j,k,iBLK)
          State_VGB(iRhoUz_I,i,j,k,iBLK)=State_VGB(iRhoUz_I,i,j,k,iBLK) &
               / State_VGB(iRho_I,i,j,k,iBLK)
       end do; end do; end do
    end if

    if(neiLnorth(iBLK) == -1)then
       if(     .not.unusedBlock_BP(neiBnorth(1,iBLK),neiPnorth(1,iBLK)) &
            .and. .not.unusedBlock_BP(neiBnorth(2,iBLK),neiPnorth(2,iBLK)) &
            .and. .not.unusedBlock_BP(neiBnorth(3,iBLK),neiPnorth(3,iBLK)) &
            .and. .not.unusedBlock_BP(neiBnorth(4,iBLK),neiPnorth(4,iBLK)))then
          do k=1,nK;do i=1,nI
             State_VGB(1:nVar,i,nJ+1,k,iBLK) =&
                  State_VGB(1:nVar,i,nJ+1,k,iBLK) + cThird*( &
                  State_VGB(1:nVar,i,nJ+1,k,iBLK) - &
                  State_VGB(1:nVar,i,nJ,k,iBLK))
             where(DefaultState_V(1:nVar) > cTiny) &
                  State_VGB(1:nVar,i,nJ+1,k,iBLK) = &
                  max(State_VGB(1:nVar,i,nJ+1,k,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLsouth(iBLK) == -1)then
       if(        .not.unusedBlock_BP(neiBsouth(1,iBLK),neiPsouth(1,iBLK)) & 
            .and. .not.unusedBlock_BP(neiBsouth(2,iBLK),neiPsouth(2,iBLK)) & 
            .and. .not.unusedBlock_BP(neiBsouth(3,iBLK),neiPsouth(3,iBLK)) & 
            .and. .not.unusedBlock_BP(neiBsouth(4,iBLK),neiPsouth(4,iBLK)))then
          do k=1,nK;do i=1,nI
             State_VGB(1:nVar,i,0,k,iBLK) = &
                  State_VGB(1:nVar,i,0,k,iBLK) + cThird*( &
                  State_VGB(1:nVar,i,0,k,iBLK) - &
                  State_VGB(1:nVar,i,1,k,iBLK))
             where(DefaultState_V(1:nVar) > cTiny) &
                  State_VGB(1:nVar,i,0,k,iBLK) = &
                  max(State_VGB(1:nVar,i,0,k,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLeast(iBLK) == -1)then
       if(     .not.unusedBlock_BP(neiBeast(1,iBLK),neiPeast(1,iBLK)) &
            .and. .not.unusedBlock_BP(neiBeast(2,iBLK),neiPeast(2,iBLK)) &
            .and. .not.unusedBlock_BP(neiBeast(3,iBLK),neiPeast(3,iBLK)) &
            .and. .not.unusedBlock_BP(neiBeast(4,iBLK),neiPeast(4,iBLK)))then
          do k=1,nK;do j=1,nJ
             State_VGB(1:nVar,0,j,k,iBLK) = &
                  State_VGB(1:nVar,0,j,k,iBLK) + cThird*(&
                  State_VGB(1:nVar,0,j,k,iBLK) - &
                  State_VGB(1:nVar,1,j,k,iBLK))
             where(DefaultState_V(1:nVar) > cTiny) &
                  State_VGB(1:nVar,0,j,k,iBLK) = &
                  max(State_VGB(1:nVar,0,j,k,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLwest(iBLK) == -1)then
       if(        .not.unusedBlock_BP(neiBwest(1,iBLK),neiPwest(1,iBLK)) &
            .and. .not.unusedBlock_BP(neiBwest(2,iBLK),neiPwest(2,iBLK)) &
            .and. .not.unusedBlock_BP(neiBwest(3,iBLK),neiPwest(3,iBLK)) &
            .and. .not.unusedBlock_BP(neiBwest(4,iBLK),neiPwest(4,iBLK)))then
          do k=1,nK;do j=1,nJ
             State_VGB(1:nVar,nI+1,j,k,iBLK) = &
                  State_VGB(1:nVar,nI+1,j,k,iBLK) + cThird*( &
                  State_VGB(1:nVar,nI+1,j,k,iBLK) - &
                  State_VGB(1:nVar,nI,  j,k,iBLK))
             where(DefaultState_V(1:nVar)>cTiny) &
                  State_VGB(1:nVar,nI+1,j,k,iBLK) = &
                  max(State_VGB(1:nVar,nI+1,j,k,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLtop(iBLK) == -1)then
       if(        .not.unusedBlock_BP(neiBtop(1,iBLK),neiPtop(1,iBLK)) &
            .and. .not.unusedBlock_BP(neiBtop(2,iBLK),neiPtop(2,iBLK)) &
            .and. .not.unusedBlock_BP(neiBtop(3,iBLK),neiPtop(3,iBLK)) &
            .and. .not.unusedBlock_BP(neiBtop(4,iBLK),neiPtop(4,iBLK)))then
          do j=1,nJ;do i=1,nI
             State_VGB(1:nVar,i,j,nK+1,iBLK) = &
                  State_VGB(1:nVar,i,j,nK+1,iBLK) + cThird*(&
                  State_VGB(1:nVar,i,j,nK+1,iBLK) - &
                  State_VGB(1:nVar,i,j,nK,iBLK))
             where(DefaultState_V(1:nVar) > cTiny) &
                  State_VGB(1:nVar,i,j,nK+1,iBLK) = &
                  max(State_VGB(1:nVar,i,j,nK+1,iBLK), 1e-30)
          end do; end do
       end if
    end if
    if(neiLbot(iBLK) == -1)then
       if(        .not.unusedBlock_BP(neiBbot(1,iBLK),neiPbot(1,iBLK)) &
            .and. .not.unusedBlock_BP(neiBbot(2,iBLK),neiPbot(2,iBLK)) &
            .and. .not.unusedBlock_BP(neiBbot(3,iBLK),neiPbot(3,iBLK)) &
            .and. .not.unusedBlock_BP(neiBbot(4,iBLK),neiPbot(4,iBLK)))then
          do j=1,nJ;do i=1,nI
             State_VGB(1:nVar,i,j,0,iBLK) = &
                  State_VGB(1:nVar,i,j,0,iBLK) + cThird*( &
                  State_VGB(1:nVar,i,j,0,iBLK) - &
                  State_VGB(1:nVar,i,j,1,iBLK))
             where(DefaultState_V(1:nVar) > cTiny) &
                  State_VGB(1:nVar,i,j,0,iBLK) = &
                  max(State_VGB(1:nVar,i,j,0,iBLK), 1e-30)
          end do; end do
       end if
    end if

    if(nFluid > 1 .and. .not.DoLimitMomentum)then
       ! Convert velocities back to momenta
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
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

