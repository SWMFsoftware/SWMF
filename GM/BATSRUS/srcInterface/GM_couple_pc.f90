!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE PC

module GM_couple_pc

  implicit none

  private ! except
  public:: GM_get_for_pc_dt
  public:: GM_get_for_pc_init
  public:: GM_get_for_pc
  public:: GM_put_from_pc

contains

  subroutine GM_get_for_pc_dt(DtSi)

    use ModMain,            ONLY: Dt, time_accurate, Time_Simulation
    use ModPhysics,         ONLY: No2Si_V, UnitT_
    use ModPIC,             ONLY: DnCouplePic
    use ModTimeStepControl, ONLY: set_global_timestep

    implicit none

    real, intent(out) ::  DtSi

    character(len=*), parameter :: NameSub='GM_get_for_pc_dt'
    !--------------------------------------------------------------------------

    ! use 1.0e12 >> dt  for limiting time 
    call set_global_timestep(1.0e12)

    DtSi = Dt*No2Si_V(UnitT_)*DnCouplePic

  end subroutine GM_get_for_pc_dt

  !============================================================================

  subroutine GM_get_for_pc_init(ParamInt_I, nParamReal, ParamReal_I)

    use ModMultiFluid, ONLY: nIonFluid, MassIon_I, ChargeIon_I
    use ModPhysics,    ONLY: No2Si_V, UnitX_, PePerPtotal
    use ModPIC,        ONLY: XyzMinPic_DI, XyzMaxPic_DI, nRegionPiC, &
         DxyzPic_DI, xUnitPicSi, uUnitPicSi, mUnitPicSi, UseFileCoupling
    use BATL_lib,      ONLY: x_, y_, z_, nDim
    use ModMain,       ONLY: lVerbose

    implicit none

    integer :: iFluid, iRegionPIC 
    integer :: i, n ! help indexes 

    integer, intent(inout) :: ParamInt_I(4)
    integer, intent(in)  :: nParamReal ! Size of ParamReal_I 
    real, optional, intent(inout) :: ParamReal_I(nParamReal)

    character(len=*), parameter :: NameSub='GM_get_for_pc_init'
    !--------------------------------------------------------------------------

    if(nParamReal < 1 ) then
       !! SEND FIRST
       ParamInt_I(1) = nIonFluid 
       ParamInt_I(2) = nRegionPic
       ParamInt_I(3) = nDim 
       !!! REMOVE
       if(UseFileCoupling) then 
          ParamInt_I(4) = 1
       else
          ParamInt_I(4) = 0
       end if

    else
       ! First part of ParamReal is the Nregion block, then Nspecis block
       n = 1
       do iRegionPic = 1, nRegionPic
          ParamReal_I(n) = XyzMinPic_DI(x_,iRegionPic); n = n+1
          ParamReal_I(n) = XyzMaxPic_DI(x_,iRegionPic); n = n+1
          ParamReal_I(n) = DxyzPic_DI(x_,iRegionPic)  ; n = n+1
          if(nDim > 1) then
             ParamReal_I(n) = XyzMinPic_DI(y_,iRegionPic); n = n+1
             ParamReal_I(n) = XyzMaxPic_DI(y_,iRegionPic); n = n+1
             ParamReal_I(n) = DxyzPic_DI(y_,iRegionPic); n = n+1
          else
             ! Single cell in Y direction with dy = dx
             ParamReal_I(n) = 0.0; n = n+1
             ParamReal_I(n:n+1) = DxyzPic_DI(x_,iRegionPic); n = n+2
          end if
          if(nDim > 2) then
             ParamReal_I(n) = XyzMinPic_DI(z_,iRegionPic); n = n+1
             ParamReal_I(n) = XyzMaxPic_DI(z_,iRegionPic); n = n+1
             ParamReal_I(n) = DxyzPic_DI(z_,iRegionPic); n = n+1
          else
             ! Single cell in Z direction with dz = max(dx, dy)
             ParamReal_I(n)     = 0.0; n = n+1
             ParamReal_I(n:n+1) = maxval(DxyzPic_DI(x_:y_,iRegionPic)); n = n+2
          endif
       end do

       ! convert to SI units
       do i = 1, n -  1
          ParamReal_I(i) = ParamReal_I(i) * No2Si_V(UnitX_)
       end do

       ! Q/Qi, M/Mi
       do iFluid = 0, nIonFluid-1
          ParamReal_I(n) = ChargeIon_I(iFluid+1); n = n+1
          ParamReal_I(n) = MassIon_I(iFluid+1)/MassIon_I(1); n = n+1
       end do

       ! Electron pressure ratio
       ParamReal_I(n) = PePerPtotal; n = n+1

       ! Units
       ParamReal_I(n) = xUnitPicSi; n = n+1
       ParamReal_I(n) = uUnitPicSi; n = n+1
       ParamReal_I(n) = mUnitPicSi
    end if

  end subroutine GM_get_for_pc_init
  !============================================================================

  subroutine GM_get_for_pc(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    ! Get magnetic field data from GM to PC

    use ModProcMH,  ONLY: iProc
    use ModPhysics, ONLY: Si2No_V, UnitX_, No2Si_V, iUnitCons_V
    use ModAdvance, ONLY: State_VGB, Bx_, Bz_, nVar, UseAnisoPressure
    use ModVarIndexes, ONLY: nVar, NamePrimitiveVar
    use ModB0,      ONLY: UseB0, get_b0
    use BATL_lib,   ONLY: nDim, MaxDim, MinIJK_D, MaxIJK_D, find_grid_block
    use ModInterpolate, ONLY: interpolate_vector

    implicit none

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    real:: Xyz_D(MaxDim), B0_D(MaxDim)
    real:: Dist_D(MaxDim), State_V(nVar)
    integer:: iCell_D(MaxDim)

    integer:: iPoint, iBlock, iProcFound

    character(len=*), parameter :: NameSub='GM_get_for_pc'
    !--------------------------------------------------------------------------

    Dist_D = -1.0
    Xyz_D  =  0.0

    do iPoint = 1, nPoint

       Xyz_D(1:nDim) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
       call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D, Dist_D, &
            UseGhostCell = .true.)

       if(iProcFound /= iProc)then
          write(*,*)NameSub,' ERROR: Xyz_D, iProcFound=', Xyz_D, iProcFound
          call stop_mpi(NameSub//' could not find position on this proc')
       end if

       State_V = interpolate_vector(State_VGB(:,:,:,:,iBlock), nVar, nDim, &
            MinIJK_D, MaxIJK_D, iCell_D=iCell_D, Dist_D=Dist_D)

       if(UseB0)then
          call get_b0(Xyz_D, B0_D)
          State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + B0_D
       end if

       Data_VI(1:nVar,iPoint) = State_V*No2Si_V(iUnitCons_V)

    end do

  end subroutine GM_get_for_pc
  !============================================================================
  subroutine GM_get_regions(NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

    ! This subrutine is called 3 times:
    !   First time we set the number of cell centered grid points 
    !       on this processor. 
    !   Second time we get the cordinates to the same time.
    !   Third time we set the state variables to the value given by PC.

    use CON_coupler, ONLY: i_proc, n_proc
    use BATL_lib,    ONLY: nDim, Xyz_DGB, CellSize_DB, nBlock, Unused_B, &
         nI, nJ, nK, MaxDim, find_grid_block
    use ModPIC,      ONLY: XyzMaxPic_DI, XyzPic0_DI, DxyzPic_DI, & 
         nRegionPic, nGhostPic
    use ModPhysics,  ONLY: No2Si_V, UnitX_
    use ModEnergy,   ONLY: calc_energy_cell
    !use ModProcMH,   ONLY: iProc

    implicit none

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real, intent(inout), allocatable, optional :: Pos_DI(:,:)  ! Positions

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    character(len=*), parameter :: NameSub='GM_get_regions'
    logical :: havenPoint

    integer :: i, j, k, iL, iR, jL, jR, kL, kR, iBlock, iPoint, iRegion
    integer :: nPicBlock

    real    :: XyzMin_D(nDim), XyzMax_D(nDim) ! min/max position of block
    real    :: XyzMinRegion_D(nDim), XyzMaxRegion_D(nDim) 
    real    :: Xyz_D(MaxDim)

    integer :: x_= 1, y_ = 1, z_ = 1 
    !--------------------------------------------------------------------------

    ! This function will be called 3 times :
    !
    ! 1) We count the number of points in all the different regions and return the value
    !
    ! 2) We received the position array where we will store the position of each point
    !
    ! 3) we revive the state variables associated with the position and a indexing array
    !    to get then in the same order as the original position array given by 2)

    ! to avoid index issues when nDim < 3 
    if(ndim > 1) y_ = 2
    if(ndim > 2) z_ = 3

    Xyz_D   = 0.0

    havenPoint = .true.
    ! first time we only find the number of points
    if(nPoint < 1 ) then
       havenPoint = .false.
       nPicBlock = 0
       nPoint = 0
    end if

    ! left and right index in the block 
    iL = 1; iR = 1; jL = 1; jR = 1; kL =1; kR = 1
    iPoint = 1
    do iRegion = 1, nRegionPic

       ! (nGhostPic +1) where +1 is from the IPIC3D ghost layor
       XyzMaxRegion_D = XyzMaxPic_DI(1:nDim,iRegion) - &
            (nGhostPic + 0.9)*DxyzPic_DI(:,iRegion)  

       XyzMinRegion_D = XyzPic0_DI(1:nDim,iRegion) + &
            (nGhostPic - 0.1)*DxyzPic_DI(:,iRegion)  


       do iBlock=1, nBlock
          if(Unused_B(iBlock)) CYCLE

          XyzMin_D = Xyz_DGB(1:nDim,1,1,1,iBlock)
          XyzMax_D = Xyz_DGB(1:nDim,nI,nJ,nK,iBlock)

          ! blocks that have no points inside the pic region iRegion
          if(any(XyzMax_D < XyzMinRegion_D).or. &
               any(XyzMin_D > XyzMaxRegion_D )) CYCLE

          ! blocks complitly inside the pic region iRegion
          if( all(XyzMax_D <= XyzMaxRegion_D .and. &
               XyzMin_D >= XyzMinRegion_D )) then

             if(present(Data_VI)) then
                call setStateVGB(1, nI, 1, nJ, 1, nK)
                call calc_energy_cell(iBlock)
             else if(havenPoint) then
                do i=1,nI; do j=1,nJ; do k=1,nK
                   ! PC domain start from origo
                   Pos_DI(1:nDim,iPoint) = &
                        Xyz_DGB(1:nDim,i,j,k,iBlock)*No2Si_V(UnitX_)
                   iPoint = iPoint + 1
                end do; end do; end do
             else
                nPicBlock = nPicBlock + 1
             end if
             CYCLE
          end if

          ! what we have left is blocks intersecting the PIC region iRegion
          iL = max(ceiling((XyzMinRegion_D(x_) - XyzMin_D(x_)) &
               /CellSize_DB(x_,iBlock)), 0) + 1
          iR = max(ceiling((XyzMaxRegion_D(x_) - XyzMin_D(x_)) &
               /CellSize_DB(x_,iBlock)), 0) 
          if(iR > nI ) iR = nI
          if(nDim > 1) then
             jL = max(ceiling((XyzMinRegion_D(y_) - XyzMin_D(y_)) &
                  /CellSize_DB(y_,iBlock)), 0) + 1
             jR = max(ceiling((XyzMaxRegion_D(y_) - XyzMin_D(y_)) &
                  /CellSize_DB(y_,iBlock)), 0) 
             if(jR > nJ ) jR = nJ
          end if
          if(nDim > 2) then
             kL = max(ceiling((XyzMinRegion_D(z_) - XyzMin_D(z_)) &
                  /CellSize_DB(z_,iBlock)), 0) + 1
             kR = max(ceiling((XyzMaxRegion_D(z_) - XyzMin_D(z_)) &
                  /CellSize_DB(z_,iBlock)), 0) 
             if(kR > nK ) kR = nK
          end if

          if(present(Data_VI)) then
             call setStateVGB(iL, iR, jL, jR, kL, kR)
             call calc_energy_cell(iBlock)
          else if(havenPoint) then
             do i=iL,iR; do j=jL,jR; do k=kL,kR
                ! PC domain start from origo
                Pos_DI(1:nDim,iPoint) = &
                     Xyz_DGB(1:nDim,i,j,k,iBlock)*No2Si_V(UnitX_)
                iPoint = iPoint + 1
             end do; end do; end do
          else
             nPoint = nPoint + (iR - iL + 1)*(jR - jL + 1)*(kR - kL + 1)
          end if

       end do ! block

    end do ! region

    if(.not. havenPoint) then
       nPoint = nPoint +  nPicBlock*nI*nJ*nK
    end if

  contains

    !==========================================================================
    subroutine setStateVGB(iLx,iRx,jLy,jRy,kLz,kRz)

      use ModAdvance,  ONLY: Bx_, Bz_
      use ModPhysics,  ONLY: Si2No_V, iUnitCons_V
      use ModAdvance,  ONLY: State_VGB
      use ModB0,       ONLY: B0_DGB, useB0

      integer, intent(in) :: iLx,iRx,jLy,jRy,kLz,kRz

      !------------------------------------------------------------------------
      do i=iLx,iRx; do j=jLy,jRy; do k=kLz,kRz         
         State_VGB(1:nVar,i,j,k,iBlock) = &
              Data_VI(1:nVar,iPoint_I(iPoint))*Si2No_V(iUnitCons_V)
         iPoint = iPoint + 1
      end do; end do; end do

      ! remove constatnt magnetic field
      if(UseB0)then
         do i=iLx,iRx; do j=jLy,jRy; do k=kLz,kRz
            State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                 State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
         end do; end do; end do
      end if

    end subroutine setStateVGB

  end  subroutine GM_get_regions


  !===========================================================================
  subroutine GM_put_from_pc( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    use CON_coupler, ONLY: i_proc, n_proc
    use BATL_lib,    ONLY: nDim
    !use ModProcMH,   ONLY: iProc

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter :: NameSub='GM_put_from_pc'
    !--------------------------------------------------------------------------

    if(.not. present(Data_VI))then
       nPoint=0;
       ! get nPoint
       call GM_get_regions(NameVar, nVar, nPoint, Pos_DI)

       if(allocated(Pos_DI)) deallocate(Pos_DI)
       allocate(Pos_DI(nDim,nPoint))

       ! get Pos_DI
       call GM_get_regions(NameVar, nVar, nPoint, Pos_DI)

       RETURN
    end if

    ! set State variables
    call GM_get_regions(NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

  end subroutine GM_put_from_pc

end module GM_couple_pc
