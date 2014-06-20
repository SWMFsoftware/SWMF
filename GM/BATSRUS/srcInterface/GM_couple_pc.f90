!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module GM_couple_pc

  implicit none

  private ! except
  public:: GM_get_for_pc_dt
  public:: GM_get_for_pc_init
  public:: GM_get_for_pc
  public:: GM_put_from_pc

contains

  subroutine GM_get_for_pc_dt(DtSi)

    use ModMain,      ONLY: Dt
    use ModPhysics,   ONLY: No2Si_V, UnitT_
    use ModPIC,       ONLY: DnCouplePic

    implicit none

    real, intent(out) ::  DtSi

    character(len=*), parameter :: NameSub='GM_get_for_pc_dt'
    !--------------------------------------------------------------------------

    DtSi = Dt*No2Si_V(UnitT_)*DnCouplePic

  end subroutine GM_get_for_pc_dt

  !==============================================================================

  subroutine GM_get_for_pc_init(ParamInt_I, n, ParamReal_I)

    use ModMultiFluid, ONLY: nIonFluid, MassIon_I, ChargeIon_I
    use ModPhysics,    ONLY: No2Si_V, UnitX_, ElectronTemperatureRatio
    use ModPIC,        ONLY: XyzMinPic_DI, XyzMaxPic_DI, nRegionPiC, &
         DxyzPic_DI, xUnitPicSi, uUnitPicSi, mUnitPicSi, UseFileCoupling
    use BATL_lib,      ONLY: x_, y_, z_, nDim
    use ModMain,       ONLY: lVerbose

    implicit none

    integer :: iFluid, iRegionPIC 
    integer :: i, idx ! help indexes 

    integer, intent(inout) :: ParamInt_I(4)
    integer, intent(in)  :: n 
    real, optional, intent(inout) :: ParamReal_I(n)

    character(len=*), parameter :: NameSub='GM_get_for_pc_init'
    !--------------------------------------------------------------------------

    if(n .lt. 1 ) then
       !! SEND FIRST
       ParamInt_I(1) = nIonFluid 
       ParamInt_I(2) = nRegionPic
       ParamInt_I(3) = nDim 
       if(UseFileCoupling) then 
          ParamInt_I(4) = 1
       else
          ParamInt_I(4) = 0
       end if

    else
       ParamReal_I(1) = XyzMinPic_DI(x_, 1)

       !! First part of ParamReal is the Nregion block, then Nspecis block
       idx = 1
       do iRegionPIC = 1, nRegionPic
          ParamReal_I(idx) = XyzMinPic_DI(x_, iRegionPIC)
          ParamReal_I(idx + 1) = XyzMaxPic_DI(x_, iRegionPIC) 
          ParamReal_I(idx + 2) = DxyzPic_DI(x_, iRegionPIC)
          if(nDim > 1) then
             ParamReal_I(idx + 3) = XyzMinPic_DI(y_, iRegionPIC)
             ParamReal_I(idx + 4) = XyzMaxPic_DI(y_, iRegionPIC)
             ParamReal_I(idx + 5) = DxyzPic_DI(y_, iRegionPIC)
          else         
             ParamReal_I(idx + 3) = 0.0
             ParamReal_I(idx + 4) = DxyzPic_DI(x_,iRegionPIC)
             ParamReal_I(idx + 5) = DxyzPic_DI(x_, iRegionPIC)
          end if
          if(nDim > 2) then
             ParamReal_I(idx + 6) = XyzMinPic_DI(z_, iRegionPIC)
             ParamReal_I(idx + 7) = XyzMaxPic_DI(z_, iRegionPIC)
             ParamReal_I(idx + 8) = DxyzPic_DI(z_, iRegionPIC)
          else         
             ParamReal_I(idx + 6) = 0.0
             ParamReal_I(idx + 7) = maxval(DxyzPic_DI(x_:y_, iRegionPIC))
             ParamReal_I(idx + 8) = maxval(DxyzPic_DI(x_:y_, iRegionPIC))
          endif
          idx = idx + 9
       end do

       ! convert to SI units
       do i = 1, idx -  1
          ParamReal_I(i) =  ParamReal_I(i) * No2Si_V(UnitX_)
       end do

       ! Q/Qi , M/Mi T/Ti
       do iFluid = 0, nIonFluid-1
          ParamReal_I(idx + iFluid*3)     = ChargeIon_I(iFluid+1)
          ParamReal_I(idx + iFluid*3 + 1) = MassIon_I(iFluid+1)
          ParamReal_I(idx + iFluid*3 + 2) = (1.0 - ElectronTemperatureRatio/nIonFluid)
       end do

       idx = idx + nIonFluid*3
       ParamReal_I(idx)   = xUnitPicSi
       ParamReal_I(idx+1) = uUnitPicSi
       ParamReal_I(idx+2) = mUnitPicSi
    end if

  end subroutine GM_get_for_pc_init
  !==============================================================================

  subroutine GM_get_for_pc(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    ! Get magnetic field data from GM to PC

    use ModProcMH,  ONLY: iProc
    use ModPhysics, ONLY: Si2No_V, UnitX_, No2Si_V, iUnitCons_V
    use ModAdvance, ONLY: State_VGB, Bx_, Bz_, nVar, UseAnisoPressure
    use ModVarIndexes, ONLY: nVar, NamePrimitiveVar
    use ModB0,      ONLY: UseB0, get_b0
    use BATL_lib,   ONLY: nDim, MaxDim, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, find_grid_block
    use ModInterpolate, ONLY: linear, bilinear, trilinear

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

       select case(nDim)
       case (1)
          State_V = linear(State_VGB(:,:,MinJ,MinK,iBlock), &
               nVar, MinI, MaxI, Xyz_D(1), iCell = iCell_D(1), Dist = Dist_D(1))
       case (2)
          State_V = bilinear(State_VGB(:,:,:,MinK,iBlock), &
               nVar, MinI, MaxI, MinJ, MaxJ, Xyz_D(1:2), iCell_D = iCell_D(1:2), Dist_D = Dist_D(1:2))
       case (3)
          State_V = trilinear(State_VGB(:,:,:,:,iBlock), &
               nVar, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, Xyz_D, iCell_D = iCell_D, Dist_D = Dist_D)
       end select

       if(UseB0)then
          call get_b0(Xyz_D, B0_D)
          State_V(Bx_:Bz_) = State_V(Bx_:Bz_) + B0_D
       end if

       Data_VI(1:nVar,iPoint) = State_V*No2Si_V(iUnitCons_V)

    end do

  end subroutine GM_get_for_pc
  !==============================================================================
  subroutine GM_get_regions(NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

    ! This subrutine is called 3 times, first time we get the number of cell centerd 
    ! grid points on this processor. Secound time we get the cordinates to the same time.
    ! Thierd time we set the state variables to the same value as given by PC componets

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
    real, intent(inout), allocatable, optional :: Pos_DI(:,:)               ! Position vectors

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    character(len=*), parameter :: NameSub='GM_get_regions'
    logical :: havenPoint

    integer :: i, j, k, iL, iR, jL, jR, kL, kR, iBlock, iPoint, iRegion
    integer :: nPicBlock

    real    :: XyzMin_D(nDim), XyzMax_D(nDim) ! min/max position of block
    real    :: XyzMinRegin_D(nDim), XyzMaxRegin_D(nDim) ! min/max positon of pic region
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
       XyzMaxRegin_D = XyzMaxPic_DI(1:nDim,iRegion) - &
            nGhostPic*DxyzPic_DI(:,iRegion)  

       XyzMinRegin_D = XyzPic0_DI(1:nDim,iRegion) + &
            nGhostPic*DxyzPic_DI(:,iRegion)  

       do iBlock=1, nBlock
          if(Unused_B(iBlock)) CYCLE

          XyzMin_D = Xyz_DGB(1:nDim,1,1,1,iBlock)
          XyzMax_D = Xyz_DGB(1:nDim,nI,nJ,nK,iBlock)

          ! blocks that have no points inside the pic region iRegion
          if(any(XyzMax_D < XyzMinRegin_D).or. &
               any(XyzMin_D > XyzMaxRegin_D )) CYCLE

          ! blocks complitly inside the pic region iRegion
          if( all(XyzMax_D <= XyzMaxRegin_D .and. &
               XyzMin_D >= XyzMinRegin_D )) then

             if(present(Data_VI)) then
                call setStateVGB(1, nI, 1, nJ, 1, nK)
                call calc_energy_cell(iBlock)
             else if(havenPoint) then
                do i=1,nI; do j=1,nJ; do k=1,nK
                   ! PC domain start from origo
                   Pos_DI(1:nDim,iPoint) = &
                        (Xyz_DGB(1:nDim,i,j,k,iBlock) - &
                        XyzPic0_DI(1:nDim,iRegion))*No2Si_V(UnitX_)
                   iPoint = iPoint + 1
                end do; end do; end do
             else
                nPicBlock = nPicBlock + 1
             end if
             CYCLE
          end if

          ! what we have left is blocks intercection the pic region iRegion
          iL = max(ceiling((XyzMinRegin_D(x_) - XyzMin_D(x_))/CellSize_DB(x_,iBlock)), 0) + 1
          iR = max(ceiling((XyzMaxRegin_D(x_) - XyzMin_D(x_))/CellSize_DB(x_,iBlock)), 0) + 1
          if(iR > nI ) iR = nI
          if(nDim > 1) then
             jL = max(ceiling((XyzMinRegin_D(y_) - XyzMin_D(y_))/CellSize_DB(y_,iBlock)), 0) + 1
             jR = max(ceiling((XyzMaxRegin_D(y_) - XyzMin_D(y_))/CellSize_DB(y_,iBlock)), 0) + 1
             if(jR > nJ ) jR = nJ
          end if
          if(nDim > 2) then
             kL = max(ceiling((XyzMinRegin_D(z_) - XyzMin_D(z_))/CellSize_DB(z_,iBlock)), 0) + 1
             kR = max(ceiling((XyzMaxRegin_D(z_) - XyzMin_D(z_))/CellSize_DB(z_,iBlock)), 0) + 1
             if(kR > nK ) kR = nK
          end if

          if(present(Data_VI)) then
             call setStateVGB(iL, iR, jL, jR, kL, kR)
             call calc_energy_cell(iBlock)
          else if(havenPoint) then
             do i=iL,iR; do j=jL,jR; do k=kL,kR
                ! PC domain start from origo
                Pos_DI(1:nDim,iPoint) = &
                     (Xyz_DGB(1:nDim,i,j,k,iBlock) - &
                     XyzPic0_DI(1:nDim,iRegion))*No2Si_V(UnitX_)
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

    !===========================================================================
    subroutine setStateVGB(iLx,iRx,jLy,jRy,kLz,kRz)

      use ModAdvance,  ONLY: Bx_, Bz_
      use ModPhysics,  ONLY: Si2No_V, iUnitCons_V
      use ModAdvance,  ONLY: State_VGB
      use ModB0,       ONLY: B0_DGB, useB0

      integer, intent(in) :: iLx,iRx,jLy,jRy,kLz,kRz

      !--------------------------------------------------------------------------
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

    !  logical,          intent(in)   :: UseData ! true when data is transferred
    ! false if positions are asked
    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:)               ! Position vectors

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
