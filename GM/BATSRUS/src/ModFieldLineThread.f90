module ModFieldLineThread
  use ModMagnetogram,   ONLY: get_magnetogram_field
  !\
  !Input arguments: x,y,z, output - B0 vector in SI units
  !/
  use ModMain, ONLY: UseFieldLineThreads
  implicit none
  save
  logical, public, allocatable :: DoThreads_B(:), IsAllocatedThread_B(:)
  ! Named indexes for local use only

  !In the boundary blocks the physical cells near the boundary are connected 
  !to the photosphere with these threads
  type BoundaryThreads
     !\
     !Three coordinates for each point of each thread (enumerated by j and k) 
     !/
     real,pointer :: Xyz_DIII(:,:,:,:)

     !\ 
     ! The thread length, \int{ds}, from the photoshere to the given point 
     !/
     real,pointer :: Length_III(:,:,:)
     !\ 
     ! The integral, \int{B ds}, from the photoshere to the given point 
     !/
     real,pointer :: BLength_III(:,:,:)
     !\ 
     ! The integral, \int{1/sqrt(B) ds} 
     !/
     real,pointer :: Length2SqrtB_III(:,:,:)
     !\
     ! Magnetic field intensity and heliocentric distance
     !/
     real,pointer :: B_III(:,:,:),R_III(:,:,:)
     !\
     ! number of points
     !/
     integer,pointer :: nPoint_II(:,:)
     !\
     ! For a given thread, if the direction from the photoshere toward
     ! the corona low boundary is parallel (SignBr_II = +1) or 
     ! antiparallel (SignBr_II = -1) to the magnetic field.
     !/
     real, pointer :: SignBr_II(:,:)
     !\
     ! For a given thread, the derivative of a temperature gradient
     ! along the thread at the starting point over the temperature
     ! in the ghost cell. Used to derive the temperature in the ghost
     ! cell from the temperature gradient along the thread at 
     ! the starting point to be found by solving the MHD equations
     ! on the thread.
     !/
     real, pointer :: DGradTeOverGhostTe_DII(:,:,:)
     !\
     ! B0 field at the boundary
     !/
     real, pointer :: B0Face_DII(:,:,:)
  end type BoundaryThreads
  type(BoundaryThreads), pointer :: BoundaryThreads_B(:)

  integer :: nPointInThreadMax
  real    :: DsThreadMin
  !\
  ! Parameters of turbulence
  !/
  real :: HeatCondParSi, HeatCondPar
contains
  subroutine read_threads(iSession)
    use ModSize, ONLY: MaxBlock
    use ModReadParam, ONLY: read_var
    integer, intent(in):: iSession
    integer :: iBlock
    !-------------

    call read_var('UseFieldLineThreads', UseFieldLineThreads)
    if(UseFieldLineThreads)then
       if(iSession/=1)call CON_stop(&
            'UseFieldLineThreads can only be set ON during the first session')
       if(.not.allocated(DoThreads_B))then
          allocate(        DoThreads_B(MaxBlock))
          allocate(IsAllocatedThread_B(MaxBlock))
          DoThreads_B = .false.
          IsAllocatedThread_B = .false.
          allocate(BoundaryThreads_B(1:MaxBlock))
          do iBlock = 1, MaxBlock
             call nullify_thread_b(iBlock)
          end do
       end if
       call read_var('nPointInThreadMax', nPointInThreadMax)
       call read_var('DsThreadMin', DsThreadMin)
    else
       if(allocated(DoThreads_B))then
          deallocate(DoThreads_B)
          do iBlock = 1, MaxBlock
             if(IsAllocatedThread_B(iBlock))&
                  call deallocate_thread_b(iBlock)
          end do
          deallocate(IsAllocatedThread_B)
          deallocate(  BoundaryThreads_B)
       end if
    end if
  end subroutine read_threads
  !=========================
  subroutine nullify_thread_b(iBlock)
    integer, intent(in) :: iBlock
    nullify(BoundaryThreads_B(iBlock) % Xyz_DIII)
    nullify(BoundaryThreads_B(iBlock) % Length_III)
    nullify(BoundaryThreads_B(iBlock) % BLength_III)
    nullify(BoundaryThreads_B(iBlock) % Length2SqrtB_III)
    nullify(BoundaryThreads_B(iBlock) % B_III)
    nullify(BoundaryThreads_B(iBlock) % R_III)
    nullify(BoundaryThreads_B(iBlock) % nPoint_II)
    nullify(BoundaryThreads_B(iBlock) % SignBr_II)
    nullify(BoundaryThreads_B(iBlock) % DGradTeOverGhostTe_DII)
    nullify(BoundaryThreads_B(iBlock) % B0Face_DII)
  end subroutine nullify_thread_b
  !=========================
  subroutine deallocate_thread_b(iBlock)
    integer, intent(in) :: iBlock
    !------
    deallocate(BoundaryThreads_B(iBlock) % Xyz_DIII)
    deallocate(BoundaryThreads_B(iBlock) % Length_III)
    deallocate(BoundaryThreads_B(iBlock) % BLength_III)
    deallocate(BoundaryThreads_B(iBlock) % Length2SqrtB_III)
    deallocate(BoundaryThreads_B(iBlock) % B_III)
    deallocate(BoundaryThreads_B(iBlock) % R_III)
    deallocate(BoundaryThreads_B(iBlock) % nPoint_II)
    deallocate(BoundaryThreads_B(iBlock) % SignBr_II)
    deallocate(BoundaryThreads_B(iBlock) % DGradTeOverGhostTe_DII)
    deallocate(BoundaryThreads_B(iBlock) % B0Face_DII)
    IsAllocatedThread_B(iBlock) = .false.
    call nullify_thread_b(iBlock)
  end subroutine deallocate_thread_b
  !========================
  subroutine set_threads
    use ModMain,     ONLY: MaxBlock, Unused_B, body1_,&
         nDim, nJ, nK
    use ModGeometry, ONLY: Xyz_DGB, IsBoundaryBlock_IB
    integer:: iBlock
    !--------
    do iBlock = 1, MaxBlock
       if(Unused_B(iBlock))then
          DoThreads_B(iBlock) = .false.
          if(IsAllocatedThread_B(iBlock))&
               call deallocate_thread_b(iBlock)
          CYCLE
       end if
       if(.not.DoThreads_B(iBlock))CYCLE
       DoThreads_B(iBlock) = .false.
       !\
       ! Check if the block is at the inner boundary
       ! Otherwise CYCLE
       !/
       if(.not.IsBoundaryBlock_IB(body1_,iBlock))then
          if(IsAllocatedThread_B(iBlock))&
               call deallocate_thread_b(iBlock)
          CYCLE
       end if
       !\
       ! Allocate threads if needed
       !/
       if(.not.IsAllocatedThread_B(iBlock))then
          allocate(BoundaryThreads_B(iBlock) % Xyz_DIII(&
               1:nDim,-nPointInThreadMax:0,1:nJ,1:nK))
          allocate(BoundaryThreads_B(iBlock) % Length_III(&
               -nPointInThreadMax:0,1:nJ,1:nK))
          allocate(BoundaryThreads_B(iBlock) % BLength_III(&
               -nPointInThreadMax:0,1:nJ,1:nK))
          allocate(BoundaryThreads_B(iBlock) % Length2SqrtB_III(&
               -nPointInThreadMax:0,1:nJ,1:nK))
          allocate(BoundaryThreads_B(iBlock) % B_III(&
               -nPointInThreadMax:0,1:nJ,1:nK))
          allocate(BoundaryThreads_B(iBlock) % R_III(&
               -nPointInThreadMax:0,1:nJ,1:nK))
          allocate(BoundaryThreads_B(iBlock) % nPoint_II(&
               1:nJ,1:nK))
          allocate(BoundaryThreads_B(iBlock) % SignBr_II(&
               1:nJ,1:nK))
          allocate(BoundaryThreads_B(iBlock) % DGradTeOverGhostTe_DII(&
               1:nDim,1:nJ,1:nK))
          allocate(BoundaryThreads_B(iBlock) % B0Face_DII(&
               1:nDim,1:nJ,1:nK))
          IsAllocatedThread_B(iBlock) = .true.
       end if
       !\
       ! The threads are now set in a just created block, or
       ! on updating B_0 field
       !/
       call set_threads_b(iBlock)
    end do
  end subroutine set_threads
  !=====================
  subroutine set_threads_b(iBlock)
    use ModGeometry, ONLY: Xyz_DGB
    use ModPhysics,  ONLY: Si2No_V, No2Si_V, UnitB_
    use ModMain,     ONLY: nDim, nJ, nK, jTest, kTest, BlkTest
    use ModNumConst, ONLY: cToleranceOrig=>cTolerance
    integer, intent(in) :: iBlock
    !\
    ! Locals:
    !/
    !Loop variable: (j,k) enumerate the cells at which
    !the threads starts, iPoint starts from negative
    !values at the photospheric end and the maximal 
    !value of this index is 0 for the thread point
    !at the physical cell center. 
    integer :: j, k, iPoint

    !\
    ! rBody here is set to one keeping a capability to set
    ! the face-formulated boundary condition by modifying
    ! rBody
    !/
    real, parameter :: rBody = 1.0

    !Length interval, !Heliocentric distance
    real :: Ds, R
    !coordinates, field vector and modulus  
    real :: Xyz_D(nDim), B0_D(nDim), B0
    !Same stored for the starting point
    real :: XyzStart_D(nDim), B0Start_D(nDim),  B0Start
    ! 1 for ourward directed field, -1 otherwise
    real ::SignBr
    !Coordinates and magnetic field in the midpoint 
    !within the framework of the Runge-Kutta scheme
    real :: XyzAux_D(nDim), B0Aux_D(nDim)
    !Aux
    real :: ROld, Aux
    !\
    !The magnetic field amplitude is fixed if close to zero.
    !Here, the magnetic field in SI units is used, therefore,
    !SI field tolerance should be used
    !/
    real:: cTolerance != cToleranceOrig*No2Si_V(UnitB_)

    integer, parameter:: R_ = 1
    real :: Dxyz_D(3) 
    !logical :: DoTest=.false., DoTestMe=.false.
    !-------------
    !call set_oktest('set_threads_b', DoTest, DoTestMe)

    cTolerance = cToleranceOrig*No2Si_V(UnitB_)
    !\
    ! Initialize threads
    !/
    BoundaryThreads_B(iBlock) % Xyz_DIII = 0.0
    BoundaryThreads_B(iBlock) % Length_III = 0.0
    BoundaryThreads_B(iBlock) % BLength_III = 0.0
    BoundaryThreads_B(iBlock) % Length2SqrtB_III = 0.0
    BoundaryThreads_B(iBlock) % B_III = 0.0
    BoundaryThreads_B(iBlock) % R_III = 0.0
    BoundaryThreads_B(iBlock) % nPoint_II = 0
    BoundaryThreads_B(iBlock) % SignBr_II = 0.0
    BoundaryThreads_B(iBlock) % DGradTeOverGhostTe_DII = 0.0
    BoundaryThreads_B(iBlock) % B0Face_DII = 0.0
    !Loop over the thread starting points
    do k = 1, nK; do j = 1, nJ
       !\
       !First, take magnetic field in the ghost cell
       !/
       XyzStart_D = Xyz_DGB(:, 0, j, k, iBlock)
       !\
       ! Magnetic field in SI!
       !/
       call get_magnetogram_field(&
            XyzStart_D(1), XyzStart_D(2), XyzStart_D(3), B0Start_D)
       BoundaryThreads_B(iBlock) % B0Face_DII(1:nDim, j, k) = B0Start_D  
       !\
       ! Starting points for all threads are in the centers
       ! of  physical cells near the boundary
       !/
       XyzStart_D = Xyz_DGB(:, 1, j, k, iBlock)
       BoundaryThreads_B(iBlock) % Xyz_DIII(&
                  1:nDim, 0, j, k) = XyzStart_D
       !\
       ! Magnetic field in SI!
       !/
       call get_magnetogram_field(&
            XyzStart_D(1), XyzStart_D(2), XyzStart_D(3), B0Start_D)
       !\
       ! Calculate and save face field
       !/
       BoundaryThreads_B(iBlock) % B0Face_DII(1:nDim, j, k) = &
            (BoundaryThreads_B(iBlock) % B0Face_DII(1:nDim, j, k) + &
            B0Start_D) * 0.5 * Si2No_V(UnitB_)

       SignBr = sign(1.0, sum(XyzStart_D*B0Start_D) )
       BoundaryThreads_B(iBlock) % SignBr_II(j, k) = SignBr

       B0Start = sqrt( sum( B0Start_D**2 ) )
       BoundaryThreads_B(iBlock) % B_III(0, j, k) = &
            B0Start*Si2No_V(UnitB_)

       R = sqrt( sum( XyzStart_D**2 ) )
       BoundaryThreads_B(iBlock) % R_III(0, j, k) = R
            
       Ds = 0.50*DsThreadMin ! To enter the grid coarsening loop
       COARSEN: do
          !\
          ! Set initial Ds or increase Ds, if previous trial fails
          !/
          Ds = Ds * 2
          iPoint = 0
          Xyz_D = XyzStart_D
          B0 = B0Start
          B0_D = B0Start_D
          POINTS: do
             iPoint = iPoint + 1
             !\
             !If the number of gridpoints in the theads is too 
             !high, coarsen the grid
             !/
             if(iPoint > nPointInThreadMax)CYCLE COARSEN
             !\
             !For the previous point given are Xyz_D, B0_D, B0
             !R is only used near the photospheric end.
             !/ 
             !Two stage Runge-Kutta
             !1. Point at the half of length interval:
             XyzAux_D = Xyz_D - 0.50*Ds*SignBr*B0_D/max(B0, cTolerance)
             
             !2. Magnetic field in this point:
             call get_magnetogram_field(&
                  XyzAux_D(1), XyzAux_D(2), XyzAux_D(3), B0Aux_D)
             
             !3. New grid point:
             Xyz_D = Xyz_D - Ds*SignBr*B0Aux_D/max(&
                  sqrt(sum(B0Aux_D**2)), cTolerance)
             R = sqrt( sum( Xyz_D**2 ) )
             if(R <= rBody)EXIT COARSEN
             !\
             ! Store a point
             !/
             BoundaryThreads_B(iBlock) % Xyz_DIII(&
                  1:nDim, -iPoint, j, k) = Xyz_D
             BoundaryThreads_B(iBlock) % R_III(-iPoint, j, k) = R
             call get_magnetogram_field(&
                  Xyz_D(1), Xyz_D(2), Xyz_D(3), B0_D)
             B0 = sqrt( sum( B0_D**2 ) )
             BoundaryThreads_B(iBlock) % B_III(-iPoint, j, k) = &
                  B0*Si2No_V(UnitB_)
          end do POINTS
       end do COARSEN
       !Calculate more accurately the intersection point
       !with the photosphere surface
       ROld = BoundaryThreads_B(iBlock) % R_III(1-iPoint, j, k)
       Aux = (ROld - RBody) / (ROld -R)
       Xyz_D =(1 - Aux)* BoundaryThreads_B(iBlock) % Xyz_DIII(&
            :, 1-iPoint, j, k) +  Aux*Xyz_D 
       !\
       ! Store the last point
       !/
       BoundaryThreads_B(iBlock) % Xyz_DIII(&
            1:nDim, -iPoint, j, k) = Xyz_D
       BoundaryThreads_B(iBlock) % R_III(-iPoint, j, k) = RBody
       call get_magnetogram_field(&
            Xyz_D(1), Xyz_D(2), Xyz_D(3), B0_D)
       B0 = sqrt( sum( B0_D**2 ) )
       BoundaryThreads_B(iBlock) % B_III(-iPoint, j, k) = &
            B0*Si2No_V(UnitB_)
       !Store the number of points
       BoundaryThreads_B(iBlock) % nPoint_II(j,k) = iPoint
       !\
       !Store the lengths
       !/
       BoundaryThreads_B(iBlock) % Length_III(&
            1-iPoint, j, k) = Ds*Aux
       BoundaryThreads_B(iBlock) % BLength_III(&
            1-iPoint, j, k) = Ds*Aux*0.50*(&
            BoundaryThreads_B(iBlock) % B_III(-iPoint, j, k) +&
            BoundaryThreads_B(iBlock) % B_III(1-iPoint, j, k) )
       BoundaryThreads_B(iBlock) % Length2SqrtB_III(&
            1-iPoint, j, k) = Ds*Aux*(&
            0.5/sqrt( BoundaryThreads_B(iBlock) % B_III(-iPoint, j, k)) +&
            0.5/sqrt( BoundaryThreads_B(iBlock) % B_III(1-iPoint, j, k)) )
       iPoint = iPoint - 1
       do while(iPoint>0)
          BoundaryThreads_B(iBlock) % Length_III(&
               1-iPoint, j, k) = BoundaryThreads_B(iBlock) % Length_III(&
               -iPoint, j, k)+ Ds
          BoundaryThreads_B(iBlock) % BLength_III(&
               1-iPoint, j, k) = BoundaryThreads_B(iBlock) % BLength_III(&
               -iPoint, j, k) + Ds*0.50*(&
               BoundaryThreads_B(iBlock) % B_III(-iPoint, j, k) +&
               BoundaryThreads_B(iBlock) % B_III(1-iPoint, j, k) )
          BoundaryThreads_B(iBlock) % Length2SqrtB_III(&
               1-iPoint, j, k) = BoundaryThreads_B(iBlock) % Length2SqrtB_III(&
               -iPoint, j, k) + Ds*(&
               0.5/sqrt( BoundaryThreads_B(iBlock) % B_III(-iPoint, j, k)) +&
               0.5/sqrt( BoundaryThreads_B(iBlock) % B_III(1-iPoint, j, k)) )
          iPoint = iPoint - 1
       end do
       !\
       !Evaluate the deravitive of temperature gradient over Te in the
       !ghost cell
       !/
       Dxyz_D = Xyz_DGB(:,1,j,k,iBlock) - Xyz_DGB(:,0,j,k,iBlock)
       BoundaryThreads_B(iBlock) % DGradTeOverGhostTe_DII(:,j,k) = &
            - Dxyz_D(1:nDim)/sum(Dxyz_D(1:nDim)**2)
            
    end do; end do
    !if(DoTestMe.and.iBlock==BlkTest)then
    !   write(*,'(a,3es18.10)')'Thread starting at the point  ',&
    !        Xyz_DGB(:,1,jTest,kTest,iBlock)
    !   write(*,'(a,3es18.10)')'Derivative of the Grad Te over ghost Te  ',&
    !        BoundaryThreads_B(iBlock) % DGradTeOverGhostTe_D(:,jTest,kTest)
    !   write(*,'(a)')'x y z B R Length BLength Length2SqrtB'
    !   do iPoint = 0, -BoundaryThreads_B(iBlock) % nPoint_II(jTest,kTest),-1
    !      write(*,'(8es18.10)')&
    !           BoundaryThreads_B(iBlock) % Xyz_DIII(:,&
    !           iPoint, jTest,kTest),&
    !           BoundaryThreads_B(iBlock) % B_III(&
    !           iPoint, jTest,kTest),&
    !           BoundaryThreads_B(iBlock) % R_III(&
    !           iPoint, jTest,kTest),&
    !           BoundaryThreads_B(iBlock) % Length_III(&
    !           iPoint, jTest, kTest),&
    !           BoundaryThreads_B(iBlock) % BLength_III(&
    !           iPoint, jTest, kTest),&
    !           BoundaryThreads_B(iBlock) % Length2SqrtB_III(&
    !           iPoint, jTest, kTest)
    !   end do
    !end if
  end subroutine set_threads_b
  !=========================================================================
  subroutine check_tr_table(iComm,TypeFileIn)
    use ModConst,      ONLY: cBoltzmann, cElectronMass, cProtonMass, &
         cEps, cElectronCharge, cTwoPi
    use ModLookupTable, ONLY: Table_I, TableType, &
         i_lookup_table, init_lookup_table, make_lookup_table

    integer, optional, intent(in) :: iComm
    character(LEN=*),optional,intent(in)::TypeFileIn

    integer:: iTable
    real,parameter:: CoulombLog = 20.0
    character(len=5)::TypeFile

    character(len=*), parameter:: NameSub = 'check_TR_table'
    !------------------------------------------------------------------------
    if(present(TypeFileIn))then
       TypeFile = TypeFileIn
    else
       TypeFile = 'ascii'
    end if
    ! electron heat conduct coefficient for single charged ions
    ! = 9.2e-12 W/(m*K^(7/2))
    HeatCondParSi = 3.2*3.0*cTwoPi/CoulombLog &
         *sqrt(cTwoPi*cBoltzmann/cElectronMass)*cBoltzmann &
         *((cEps/cElectronCharge)*(cBoltzmann/cElectronCharge))**2

    iTable =  i_lookup_table('TR')
    if(iTable >= 0)return

    ! initialize the EOS table with the default parameters
    call init_lookup_table(                                      &
         NameTable = 'TR',                                       &
         NameCommand = 'save',                                   &
         NameVar = 'logTe logNe LPe UHeat',                      &
         nIndex_I = (/500,2/),                                   &
         IndexMin_I =(/1.0e4, 1.0e8/),                           &
         IndexMax_I =(/1.0e8, 1.0e18/),                          &
         NameFile = 'TR.dat',                                    &
         TypeFile = TypeFile,                                    &
         StringDescription = &
         'Model for transition region: [K] [1/m3] [N/m] [m/s]')
    
    iTable = i_lookup_table('TR')
         

    ! The table is now initialized. 

    !Fill in the table
    call make_lookup_table(iTable, calc_tr_table)
  end subroutine check_tr_table
  !===========================================================================
  subroutine calc_tr_table(iTableIn, Arg1, Arg2, Value_V)
    use ModLookupTable, ONLY: interpolate_lookup_table, i_lookup_table
    integer, intent(in):: iTableIn
    real, intent(in)   :: Arg1, Arg2
    real, intent(out)  :: Value_V(:)
    integer:: iTable, iTe
    real   :: DeltaLogTe, LambdaCgs_V(1)
    real,dimension(1:500):: TeSi_I, LambdaSi_I, LPe_I, UHeat_I

    logical, save:: IsFirstCall = .true.
    ! at the moment, radcool not a function of Ne, but need a dummy 2nd
    ! index, and might want to include Ne dependence in table later.
    ! Table variable should be normalized to radloss_cgs * 10E+22
    ! since we don't want to deal with such tiny numbers 
    real, parameter :: RadNorm = 1.0E+22
    real, parameter :: Cgs2SiEnergyDens = &
       1.0e-7&   !erg = 1e-7 J
       /1.0e-6    !cm3 = 1e-6 m3 
    real, parameter :: Radcool2Si = 1.0e-12 & ! (cm-3=>m-3)**2
         /RadNorm*Cgs2SiEnergyDens
    !--------------------------------------------------------------------------
 
      
    !-------------------------------------------------------------------------
    if(IsFirstCall)then
       IsFirstCall=.false.
       iTable = i_lookup_table('radcool')
       if(iTable<=0)&
            call CON_stop('Transition region table needs radcool table')
       DeltaLogTe = log(1.0e8/1.0e4)/499 !Log(MaxTe/MinTe)/(nPoint-1)

       TeSi_I(1) = 1.0e4; LPe_I(1) = 0.0; UHeat_I(1) = 0.0
       do iTe = 2,500 
          TeSi_I(iTe) = TeSi_I(iTe-1)*exp(DeltaLogTe) 
       end do

       do iTe = 1,500
          call interpolate_lookup_table(iTable,&
               TeSi_I(iTe), 1.0e2, LambdaCgs_V)
          LambdaSi_I(iTe) = LambdaCgs_V(1)*Radcool2Si
          if(iTe==1)CYCLE
          !\
          ! Integrate \sqrt{2\int{\kappa_0\Lambda Te**1.5 d(log T)}}/k_B
          !/
          UHeat_I(iTe) = UHeat_I(iTe-1) + &
               (LambdaSi_I(iTe-1)*TeSi_I(iTe-1)**1.50 + &
               LambdaSi_I(iTe)*TeSi_I(iTe)**1.50)*DeltaLogTe
       end do
       UHeat_I = sqrt(HeatCondParSi*UHeat_I)/cBoltzmann

       !\
       ! Temporal fix to avoid a singularity in the first point
       !/
       UHeat_I(1) = UHeat_I(2)
       do iTe = 2,500
          !\
          ! Integrate \int{\kappa_0\Lambda Te**3.5 d(log T)/UHeat}
          !/
          LPe_I(iTe) = LPe_I(iTe-1) + 0.5*DeltaLogTe* &
               ( TeSi_I(iTe-1)**3.50/UHeat_I(iTe-1) + &
               TeSi_I(iTe)**3.50/UHeat_I(iTe) )
       end do
       LPe_I = LPe_I*HeatCondParSi
       UHeat_I(1) = 0.0
    end if
    iTe = 1 + nint(log(Arg1/1.0e4)/DeltaLogTe)
    Value_V = (/ LPe_I(iTe), UHeat_I(iTe) /)
  end subroutine calc_tr_table
  !=========================================================================
  subroutine set_field_line_thread_bc(nGhost, iBlock, nVarState, State_VG, &
               iImplBlock)
    use BATL_lib, ONLY:  MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    integer, intent(in):: nGhost
    integer, intent(in):: iBlock
    integer, intent(in):: nVarState
    real, intent(inout):: State_VG(nVarState,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    ! Optional arguments when called by semi-implicit scheme
    integer, optional, intent(in):: iImplBlock
  end subroutine set_field_line_thread_bc
end module ModFieldLineThread
