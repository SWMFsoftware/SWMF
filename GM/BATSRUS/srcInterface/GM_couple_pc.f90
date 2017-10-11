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

    ! Calculate the global time step for PC
    
    use ModMain,            ONLY: Dt
    use ModPhysics,         ONLY: No2Si_V, UnitT_
    use ModTimeStepControl, ONLY: set_global_timestep

    real, intent(out) ::  DtSi
    !--------------------------------------------------------------------------
    ! use -1.0 so that no limit is applied on Dt
    call set_global_timestep(TimeSimulationLimit=-1.0)

    DtSi = Dt*No2Si_V(UnitT_)

  end subroutine GM_get_for_pc_dt

  !============================================================================

  subroutine GM_get_for_pc_init(nParamInt, nParamReal, iParam_I, Param_I)

    use ModVarIndexes, ONLY: MassSpecies_V, Pe_, Bx_, Ex_, &
         SpeciesFirst_, SpeciesLast_, nVar
    use ModMultiFluid, ONLY: nIonFluid, MassIon_I, ChargeIon_I, &
         iRhoIon_I, iRhoUxIon_I, iPparIon_I, iPIon_I
    use ModAdvance,    ONLY: UseMultiSpecies, nSpecies
    use ModPhysics,    ONLY: No2Si_V, UnitX_, PePerPtotal, rPlanetSi
    use ModPIC,        ONLY: XyzMinPic_DI, nRegionPiC, &
         DxyzPic_DI, xUnitPicSi, uUnitPicSi, mUnitPicSi, LenPic_DI, R_DDI
    use BATL_lib,      ONLY: x_, y_, z_, nDim

    integer, intent(inout) :: nParamInt, nParamReal

    integer, optional, intent(out):: iParam_I(nParamInt)
    real,    optional, intent(out):: Param_I(nParamReal)

    integer :: iFluid, iSpecies, iRegion 
    integer :: i, j, n
    !--------------------------------------------------------------------------
    if(.not.present(iParam_I))then
       ! nDim, nRegion
       ! nVar, nIonFluid, nSpecies
       ! iPe, iBx, iEx
       ! (iRho, iRhoUx, iPpar, iP) for each ion fluid
       nParamInt = 8 + nIonFluid*4

       ! Charge and mass per species/fluid
       ! XyzMin_D + LenPic_D + Dxzy_D + R_DD = 18 variables for each region
       ! PePerPtotal + xUnitPicSi + uUnitPicSi + mUnitPicSi + rPlanet
       ! + No2Si_V(UnitX_) = 6 variables
       nParamReal = max(nSpecies, nIonFluid)*2 + nRegionPic*18 + 6

       RETURN
    end if

    ! Set the integer parameters
    iParam_I(1:8) = (/nDim, nRegionPic, nVar, nIonFluid, nSpecies, Pe_, Bx_, Ex_/)
    n = 9
    iParam_I(n:n+nIonFluid-1) = iRhoIon_I;   n = n + nIonFluid
    iParam_I(n:n+nIonFluid-1) = iRhoUxIon_I; n = n + nIonFluid
    iParam_I(n:n+nIonFluid-1) = iPparIon_I;  n = n + nIonFluid
    iParam_I(n:n+nIonFluid-1) = iPIon_I

    ! First part of ParamReal defines the PIC region geometry
    ! The second part describes the species/fluids
    n = 1
    do iRegion = 1, nRegionPic
       Param_I(n) = XyzMinPic_DI(x_,iRegion)*No2Si_V(UnitX_); n = n+1
       Param_I(n) = LenPic_DI(x_,iRegion)*No2Si_V(UnitX_);    n = n+1
       Param_I(n) = DxyzPic_DI(x_,iRegion)*No2Si_V(UnitX_);   n = n+1
       if(nDim > 1) then
          Param_I(n) = XyzMinPic_DI(y_,iRegion)*No2Si_V(UnitX_); n = n+1
          Param_I(n) = LenPic_DI(y_,iRegion)*No2Si_V(UnitX_);    n = n+1
          Param_I(n) = DxyzPic_DI(y_,iRegion)*No2Si_V(UnitX_);   n = n+1
       else
          ! Single cell in Y direction with dy = dx
          Param_I(n) = 0.0; n = n+1
          Param_I(n:n+1) = DxyzPic_DI(x_,iRegion)*No2Si_V(UnitX_); n = n+2
       end if
       if(nDim > 2) then
          Param_I(n) = XyzMinPic_DI(z_,iRegion)*No2Si_V(UnitX_); n = n+1
          Param_I(n) = LenPic_DI(z_,iRegion)*No2Si_V(UnitX_); n = n+1
          Param_I(n) = DxyzPic_DI(z_,iRegion)*No2Si_V(UnitX_); n = n+1
       else
          ! Single cell in Z direction with dz = max(dx, dy)
          Param_I(n)     = 0.0; n = n+1
          Param_I(n:n+1) = maxval(DxyzPic_DI(x_:y_,iRegion))*No2Si_V(UnitX_); n = n+2
       endif

       ! The rotation matrix
       do i = 1, 3; do j = 1, 3
          if (i <= nDim .and. j<= nDim) then
             Param_I(n) = R_DDI(j,i,iRegion)
          else
             Param_I(n) = 0
          endif
          n = n+1
       enddo; enddo              
    end do

    ! Send charge and mass of each species/fluids
    if(UseMultiSpecies)then
       do iSpecies = SpeciesFirst_, SpeciesLast_
          Param_I(n) = 1.0; n = n+1                ! Charge is always 1
          Param_I(n) = MassSpecies_V(iSpecies); n = n+1
       end do
    else
       do iFluid = 1, nIonFluid
          Param_I(n) = ChargeIon_I(iFluid); n = n+1
          Param_I(n) = MassIon_I(iFluid); n = n+1
       end do
    end if

    ! Electron pressure ratio
    Param_I(n) = PePerPtotal; n = n+1

    ! Units
    Param_I(n) = xUnitPicSi; n = n+1
    Param_I(n) = uUnitPicSi; n = n+1
    Param_I(n) = mUnitPicSi; n = n+1
    Param_I(n) = rPlanetSi;  n = n+1
    Param_I(n) = No2Si_V(UnitX_)
    
  end subroutine GM_get_for_pc_init
  !============================================================================

  subroutine GM_get_for_pc(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    ! Interpolate Data_VI from GM at the list of positions Xyz_DI 
    ! required by PC

    use ModProcMH,  ONLY: iProc
    use ModPhysics, ONLY: Si2No_V, UnitX_, No2Si_V, iUnitCons_V
    use ModAdvance, ONLY: State_VGB, Bx_, Bz_, nVar
    use ModVarIndexes, ONLY: nVar
    use ModB0,      ONLY: UseB0, get_b0
    use BATL_lib,   ONLY: nDim, MaxDim, MinIJK_D, MaxIJK_D, find_grid_block
    use ModInterpolate, ONLY: interpolate_vector
    use ModIO, ONLY: iUnitOut

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

    integer, allocatable, save:: iBlockCell_DI(:,:)
    real,    allocatable, save:: Dist_DI(:,:)

    integer:: iPoint, iBlock, iProcFound

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='GM_get_for_pc'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! If nDim < MaxDim, make sure that all elements are initialized
    Dist_D = -1.0
    Xyz_D  =  0.0

    if(IsNew)then
       if(DoTest)write(iUnitOut,*) NameSub,': iProc, nPoint=', iProc, nPoint

       if(allocated(iBlockCell_DI)) deallocate(iBlockCell_DI, Dist_DI)
       allocate(iBlockCell_DI(0:nDim,nPoint), Dist_DI(nDim,nPoint))

       do iPoint = 1, nPoint

          Xyz_D(1:nDim) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)
          call find_grid_block(Xyz_D, iProcFound, iBlock, iCell_D, Dist_D, &
               UseGhostCell = .true.)

          if(iProcFound /= iProc)then
             write(*,*) NameSub,' ERROR: Xyz_D, iProcFound=', Xyz_D, iProcFound
             call stop_mpi(NameSub//' could not find position on this proc')
          end if

          ! Store block and cell indexes and distances for interpolation
          iBlockCell_DI(0,iPoint)      = iBlock
          iBlockCell_DI(1:nDim,iPoint) = iCell_D(1:nDim)
          Dist_DI(:,iPoint)            = Dist_D(1:nDim)

       end do
    end if

    do iPoint = 1, nPoint

       Xyz_D(1:nDim) = Xyz_DI(:,iPoint)*Si2No_V(UnitX_)

       ! Use stored block and cell indexes and distances
       iBlock          = iBlockCell_DI(0,iPoint)
       iCell_D(1:nDim) = iBlockCell_DI(1:nDim,iPoint)
       Dist_D(1:nDim)  = Dist_DI(:,iPoint)

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
    use BATL_lib,     ONLY: nDim, Xyz_DGB, nBlock, Unused_B, &
         nI, nJ, nK
    use ModPIC,       ONLY: DxyzPic_DI, LenPic_DI, & 
         nRegionPic, nGhostPic, mhd_to_pic_vec
    use ModPhysics,   ONLY: No2Si_V, UnitX_, Si2No_V, iUnitCons_V
    use ModMain,      ONLY: UseB0, UseHyperbolicDivB
    use ModB0,        ONLY: B0_DGB
    use ModAdvance,   ONLY: State_VGB, Bx_, Bz_, Hyp_, HypE_
    use ModMultiFluid,ONLY: nIonFluid
    use ModEnergy,    ONLY: calc_energy
    use ModVarIndexes,ONLY: DefaultState_V
    !use ModProcMH,   ONLY: iProc

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real, intent(inout), allocatable, optional :: Pos_DI(:,:)  ! Positions

    real,    intent(in), optional:: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)! Order of data

    logical :: DoCountOnly
    integer :: i, j, k, iBlock, iPoint, iRegion, iVar
    real    :: XyzMinRegion_D(nDim), XyzMaxRegion_D(nDim) 

    real    :: State_V(nVar)

    real    :: CoordPic_D(nDim)
    
    character(len=*), parameter :: NameSub='GM_get_regions'
    !--------------------------------------------------------------------------

    ! This function will be called 3 times :
    !
    ! 1) Count points inside all the PIC regions and return the value in nPoint
    !
    ! 2) Return the Xyz_DGB coordinates of grid cells (points) for PC in Pos_DI
    !
    ! 3) Recieve Data_VI from PC for each point and put them into State_VGB.
    !    The indexing array iPoint_I is needed to maintain the same order as
    !    the original position array Pos_DI was given in 2)

    ! to avoid index issues when nDim < 3 

    DoCountOnly = nPoint < 1
    if(DoCountOnly) nPoint = 0

    iPoint = 1
    do iRegion = 1, nRegionPic

       ! (nGhostPic +1) where +1 is from the IPIC3D ghost layer
       ! XyzMaxRegion_D and XyzMinRegion_D are in the rotated PIC coordinates.
       XyzMaxRegion_D = LenPic_DI(1:nDim,iRegion) - &
            (nGhostPic + 0.9)*DxyzPic_DI(:,iRegion)  

       XyzMinRegion_D = (nGhostPic + 0.9)*DxyzPic_DI(:,iRegion)  

       do iBlock=1, nBlock
          if(Unused_B(iBlock)) CYCLE
          
          do k = 1,nK; do j = 1,nJ; do i=1,nI

             ! Intersection with body is not handled yet ???

             ! Check if cell center is inside the PIC region
             call mhd_to_pic_vec(iRegion, Xyz_DGB(1:nDim,i,j,k,iBlock), &
                  CoordPic_D)             
             if(any(CoordPic_D > XyzMaxRegion_D)) CYCLE             
             if(any(CoordPic_D < XyzMinRegion_D)) CYCLE

             if(DoCountOnly)then
                nPoint = nPoint + 1
                CYCLE
             end if

             if(present(Data_VI))then
                ! Put Data_VI obtained from PC into State_VGB

                State_V = State_VGB(:,i,j,k,iBlock)
                do iVar = 1, nVar
                   ! Skip hyperbolic scalar variables
                   if(UseHyperbolicDivB .and. iVar==Hyp_)  CYCLE
                   if(HypE_ > 1         .and. iVar==Hype_) CYCLE
                   State_VGB(iVar,i,j,k,iBlock) = &
                        Data_VI(iVar,iPoint_I(iPoint))*Si2No_V(iUnitCons_V(iVar))
                end do
                if(UseB0) State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                     State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)

                do iVar = 1, nVar
                   ! Check for positivity
                   if(DefaultState_V(iVar) > 0 .and. State_VGB(iVar,i,j,k,iBlock) <= 0) then
                      ! Use original MHD state if PC state is not positive
                      State_VGB(:,i,j,k,iBlock) = State_V
                      EXIT
                   endif
                enddo                
                
                call calc_energy(i,i,j,j,k,k,iBlock,1,nIonFluid)
             else
                ! Provide position to PC
                Pos_DI(1:nDim,iPoint) = &
                     Xyz_DGB(1:nDim,i,j,k,iBlock)*No2Si_V(UnitX_)
             end if
             iPoint = iPoint + 1

          end do; end do; end do
       end do

    end do

  end subroutine GM_get_regions
  !===========================================================================
  subroutine GM_put_from_pc( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

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
