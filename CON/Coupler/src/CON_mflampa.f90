!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_mflampa
  ! allocate the grid used in this model
  use ModUtilities,      ONLY: check_allocate
  use ModConst,          ONLY: cBoltzmann
  implicit none
  ! Number of variables in the state vector and the identifications
  integer, parameter :: Lower_=0, Upper_=1, nDim = 3, nMHData = 13,   &
       LagrID_     = 0, & ! Lagrangian id           ^saved/   ^set to 0
       X_          = 1, & !                         |read in  |in copy_
       Y_          = 2, & ! Cartesian coordinates   |restart  |old_stat
       Z_          = 3, & !                         v/        |saved to
       Rho_        = 4, & ! Background plasma density         |mhd1
       T_          = 5, & ! Background temperature            |
       Ux_         = 6, & !                                   |may be
       Uy_         = 7, & ! Background plasma bulk velocity   |read from
       Uz_         = 8, & !                                   |mhd1
       Bx_         = 9, & !                                   |or
       By_         =10, & ! Background magnetic field         |received
       Bz_         =11, & !                                   |from
       Wave1_      =12, & !                                   |coupler
       Wave2_      =13    ! Alfven wave turbulence            v
  !
  ! State vector is a pointer, which is joined to a target array
  !
  real, allocatable, target:: MHData_VIB(:,:,:)
  !
  ! nParticle_B is a pointer, which is joined to a target array
  ! For stand alone version the target array is allocated here
  !
  integer, allocatable, target:: nParticle_B(:)
  !
  integer, allocatable :: iOffset_B(:)
  ! Grid integer parameters:
  integer :: nParticleMax, nBlock
  ! Misc:
  integer:: iError
  ! coupling parameters:
  ! domain boundaries
  real, public :: rInterfaceMin, rInterfaceMax
  ! buffer boundaries located near lower (Lo) or upper (Up) boudanry of domain
  real, public :: rBufferLo, rBufferUp
  !Coefficient to transform energy
  real         :: EnergySi2Io = 1.0/cBoltzmann
contains
  !============================================================================
  subroutine set_state_pointer(rPointer_VIB, nPointer_B, &
       nBlockIn, nParticleIn, EnergySi2IoIn)

    real,    intent(inout), pointer :: rPointer_VIB(:,:,:)
    integer, intent(inout), pointer :: nPointer_B(:)
    integer, intent(in)             :: nBlockIn, nParticleIn
    real,    optional,   intent(in) :: EnergySi2IoIn
    integer :: iParticle
    character(len=*), parameter:: NameSub = 'set_state_pointer'
    !--------------------------------------------------------------------------
    !
    ! Store nBlock and nParticleMax
    nBlock    = nBlockIn
    nParticleMax = nParticleIn
    allocate(MHData_VIB(LagrID_:nMHData, 1:nParticleMax, 1:nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'MHData_VIB')
    rPointer_VIB => MHData_VIB
    !
    MHData_VIB = 0.0
    !
    ! reset lagrangian ids
    !
    do iParticle = 1, nParticleMax
       MHData_VIB(LagrID_, iParticle, 1:nBlock) = real(iParticle)
    end do
    !
    allocate(nParticle_B(1:nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'nParticle_B')
    nPointer_B => nParticle_B
    !
    nParticle_B = 0
    if(present(EnergySi2IoIn))EnergySi2Io = EnergySi2IoIn
    allocate(iOffset_B(nBlock)); iOffset_B = 0
  end subroutine set_state_pointer
  !============================================================================
  subroutine get_bounds(iModelIn, rMinIn, rMaxIn, &
       rBufferLoIn, rBufferUpIn)
    integer,        intent(in) :: iModelIn
    real,           intent(in) :: rMinIn, rMaxIn
    real, optional, intent(in) :: rBufferLoIn, rBufferUpIn
    ! set domain boundaries
    rInterfaceMin = rMinIn; rInterfaceMax = rMaxIn
    ! set buffer boundaries
    if(present(rBufferLoIn))then
       rBufferLo = rBufferLoIn
    else
       rBufferLo = rMinIn
    end if
    if(present(rBufferUpIn))then
       rBufferUp = rBufferUpIn
    else
       rBufferUp = rMaxIn
    end if
  end subroutine get_bounds
  !============================================================================
  subroutine SP_put_from_mh(nPartial,iPutStart,Put,W,DoAdd,Buff_I,nVar)
    use CON_coupler, ONLY: &
         iVar_V, DoCoupleVar_V, &
         Density_, RhoCouple_, Pressure_, PCouple_, &
         Momentum_, RhoUxCouple_, RhoUzCouple_, &
         BField_, BxCouple_, BzCouple_, &
         Wave_, WaveFirstCouple_, WaveLastCouple_
    use CON_router, ONLY: IndexPtrType, WeightPtrType
    use ModConst,   ONLY: cProtonMass 
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer:: iRho, iP, iMx, iMz, iBx, iBz, iWave1, iWave2
    integer:: i, iBlock
    integer:: iPartial
    real:: Weight
    real:: R, Aux

    character(len=100):: StringError
    character(len=*), parameter:: NameSub = 'SP_put_from_mh'
    !--------------------------------------------------------------------------
    ! check consistency of DoCoupleVar_V
    if(.not. DoCoupleVar_V(Density_) .and. &
         (DoCoupleVar_V(Pressure_) .or. DoCoupleVar_V(Momentum_)))&
         call CON_Stop(NameSub//': pressure or momentum is coupled,'&
         //' but density is not')
    ! indices of variables in the buffer
    iRho  = iVar_V(RhoCouple_)
    iP    = iVar_V(PCouple_)
    iMx   = iVar_V(RhoUxCouple_)
    iMz   = iVar_V(RhoUzCouple_)
    iBx   = iVar_V(BxCouple_)
    iBz   = iVar_V(BzCouple_)
    iWave1= iVar_V(WaveFirstCouple_)
    iWave2= iVar_V(WaveLastCouple_)
    Aux = 0
    if(DoAdd)Aux = 1.0
    do iPartial = 0, nPartial-1
       ! cell and block indices
       i      = Put%iCB_II(1, iPutStart + iPartial)
       iBlock = Put%iCB_II(4, iPutStart + iPartial)
       ! interpolation weight
       Weight = W%Weight_I(   iPutStart + iPartial)
       if(is_in_buffer_xyz(Lower_,MHData_VIB(X_:Z_,i,iBlock)))then
          R = sqrt(sum(MHData_VIB(X_:Z_,i,iBlock)**2))
          Aux = 1.0
          Weight = Weight * (0.50 + 0.50*tanh(2*(2*R - &
               RBufferLo - RInterfaceMin)/(RBufferLo - RInterfaceMin)))
       end if
       if(is_in_buffer_xyz(Upper_,MHData_VIB(X_:Z_,i,iBlock)))then
          R = sqrt(sum(MHData_VIB(X_:Z_,i,iBlock)**2))
          Weight = Weight * (0.50 - 0.50*tanh(2*(2*R - &
               RInterfaceMax - RBufferUp)/(RInterfaceMax - RBufferUp)))
       end if
       ! put the data
       if(DoCoupleVar_V(Density_))&
            MHData_VIB(Rho_,i,iBlock) = Aux*MHData_VIB(Rho_,i,iBlock) &
            + Buff_I(iRho)/cProtonMass*Weight
       if(DoCoupleVar_V(Pressure_))&
            MHData_VIB(T_,i,iBlock) = Aux*MHData_VIB(T_,i,iBlock) + &
            Buff_I(iP)*(cProtonMass/Buff_I(iRho))*EnergySi2Io*Weight
       if(DoCoupleVar_V(Momentum_))&
            MHData_VIB(Ux_:Uz_,i,iBlock) = Aux*MHData_VIB(Ux_:Uz_,i,iBlock) + &
            Buff_I(iMx:iMz) / Buff_I(iRho) * Weight
       if(DoCoupleVar_V(BField_))&
            MHData_VIB(Bx_:Bz_,i,iBlock) = Aux*MHData_VIB(Bx_:Bz_,i,iBlock) + &
            Buff_I(iBx:iBz) * Weight
       if(DoCoupleVar_V(Wave_))&
            MHData_VIB(Wave1_:Wave2_,i,iBlock) = &
            Aux*MHData_VIB(Wave1_:Wave2_,i,iBlock) + &
            Buff_I(iWave1:iWave2)*Weight
    end do
  end subroutine SP_put_from_mh
  !============================================================================
  function is_in_buffer_r(iBuffer, R) Result(IsInBuffer)
    integer,intent(in) :: iBuffer
    real,   intent(in) :: R
    logical:: IsInBuffer
    !--------------------------------------------------------------------------
    select case(iBuffer)
    case(Lower_)
       IsInBuffer = R >= rInterfaceMin .and. R < rBufferLo
    case(Upper_)
       IsInBuffer = R >= rBufferUp .and. R < rInterfaceMax
    case default
       call CON_stop("ERROR: incorrect call of SP_wrapper:is_in_buffer")
     end select
  end function is_in_buffer_r
  !============================================================================
  function is_in_buffer_xyz(iBuffer, Xyz_D) Result(IsInBuffer)
    integer,intent(in):: iBuffer
    real,   intent(in) :: Xyz_D(nDim)
    logical:: IsInBuffer
    real:: R
    !--------------------------------------------------------------------------
    R = sqrt(sum(Xyz_D**2))
    IsInBuffer = is_in_buffer_r(iBuffer, R)
  end function is_in_buffer_xyz
  !==========================================================================
end module CON_mflampa

