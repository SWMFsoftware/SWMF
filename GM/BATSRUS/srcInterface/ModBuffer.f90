module ModBuffer
  use ModMain,     ONLY: nPhiBuff, nThetaBuff, RBuffMin, RBuffMax
  use ModNumConst, ONLY: cPi, cTwoPi
  use CON_global_vector,   ONLY: ubound_vector, point_state_v
  use CON_grid_descriptor, ONLY: GridDescriptorType, Nodes_, &
       bilinear_interpolation, DomainDecompositionType, is_proc0, &
       init_decomposition, get_root_decomposition, complete_grid, &
       set_standard_grid_descriptor, bcast_decomposition
  use CON_coupler,         ONLY: SC_, IH_
  implicit none
  save

  !This is a name of a global vector, which
  !is in the global vector storage
  character(LEN=10)::NameBuffer


  type(DomainDecompositionType),target::&
       LocalBufferDD
  type(GridDescriptorType)::LocalBufferGD

  integer:: SourceID_ = SC_, TargetID_ = IH_
  
  integer::nVarBuff=-1
  logical::DoInit =.true.
  real,dimension(:),allocatable::Buffer_V

contains

  subroutine set_buffer_name(NameIn,IDIn)
    character(LEN=*),intent(in)::NameIn
    integer,optional,intent(in)::IDIn
    !-------------

    NameBuffer=NameIn
    if(present(IDIn))SourceID_ = IDIn
  end subroutine set_buffer_name
  !====================================================
  subroutine set_spher_buffer_grid(DD,CompID_,IsLocal)
    type(DomainDecompositionType),&
         intent(out)::DD
    integer,intent(in)::CompID_
    logical,intent(in)::IsLocal
    !----------------
    call init_decomposition(&
         DD,CompID_,nDim=3,IsLocal=IsLocal)

    if(is_proc0(CompID_).or.IsLocal)call get_root_decomposition(&
         DD,&
         iRootMapDim_D=(/1,1,1/),&
         XyzMin_D=(/RBuffMin, 0.0, 0.0/),&
         XyzMax_D=(/RBuffMax,cTwoPi,cPi/), &
         nCells_D=(/1,nPhiBuff,nThetaBuff/),&
         PE_I=(/0/))
    if(.not.IsLocal)then
       call bcast_decomposition(DD)
    else
       call complete_grid(DD)
       call set_standard_grid_descriptor(DD,Standard_=Nodes_,&
            nGhostGridPoints=1,  &
            GridDescriptor=LocalBufferGD)
    end if
  end subroutine set_spher_buffer_grid
end module ModBuffer

!=============================================================!

subroutine get_from_spher_buffer_grid(XyzTarget_D,nVar,State_V)
  use ModBuffer
  use ModMain,       ONLY: nDim, R_, Phi_, Theta_, x_, y_, z_,&
                           TypeCoordSystem, Time_Simulation
  use ModAdvance,    ONLY: UseElectronPressure, UseAnisoPressure
  use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Ux_, Uz_, Bx_, Bz_, p_, &
                           WaveFirst_, WaveLast_, Pe_, Ppar_, nFluid, &
                           UseMultiSpecies
  use CON_coupler,   ONLY: Grid_C, DoCoupleVar_V, iVar_V, nVarCouple,&
                           Bfield_, ElectronPressure_, AnisoPressure_, Wave_,&
                           MultiFluid_, MultiSpecie_, &
                           RhoCouple_, RhoUxCouple_, RhoUzCouple_, PCouple_, &
                           BxCouple_, BzCouple_, PeCouple_, PparCouple_, &
                           WaveFirstCouple_, WaveLastCouple_
  use CON_axes,      ONLY: transform_matrix, transform_velocity
  use ModPhysics,    ONLY: No2Si_V,Si2No_V,UnitRho_,UnitU_,UnitB_,UnitP_,UnitX_
  use ModPhysics,    ONLY: UnitEnergyDens_
  use ModMultiFluid, ONLY: IsFullyCoupledFluid
  implicit none

  integer,intent(in)              :: nVar
  real,intent(in)                 :: XyzTarget_D(nDim)
  real,dimension(nVar),intent(out)::State_V
  real,dimension(nDim)            ::Sph_D

  character (len=3) :: TypeCoordSource
  real, save        :: SourceTarget_DD(nDim, nDim)
  real              :: TimeSimulationLast = -1.0
  real              :: XyzSource_D(nDim)

  !---------------------------------------------------------------------------
  if(DoInit)then
     DoInit = .false.
     allocate(Buffer_V(nVarCouple))
  end if

  TypeCoordSource = Grid_C(SourceID_) % TypeCoord

  if(TypeCoordSource /= TypeCoordSystem) then
     ! Convert target coordinates to the coordiante system of the model

     if(Time_Simulation > TimeSimulationLast)then
        SourceTarget_DD = &
             transform_matrix(Time_Simulation, TypeCoordSystem, TypeCoordSource)
        TimeSimulationLast = Time_Simulation
     end if
     XyzSource_D = matmul(SourceTarget_DD, XyzTarget_D)
  else
     XyzSource_D = XyzTarget_D
  end if

  ! Convert to left handed spherical coordinates !!!
  call xyz_to_spherical(XyzSource_D(x_),XyzSource_D(y_),XyzSource_D(z_),&
       Sph_D(R_),Sph_D(Phi_),Sph_D(Theta_))

  ! Get the target state from the spherical buffer grid
  Buffer_V=point_state_v(&
       NameBuffer,&
       nVarCouple,&
       nDim,    &
       Sph_D,   &
       LocalBufferGD,&
       bilinear_interpolation)

  State_V(Rho_) = Buffer_V(iVar_V(RhoCouple_))
  !Transform to primitive variables
  State_V(Ux_:Uz_) = &
       Buffer_V(iVar_V(RhoUxCouple_):iVar_V(RhoUzCouple_))&
       /Buffer_V(iVar_V(RhoCouple_))
  State_V(Bx_:Bz_) = Buffer_V(iVar_V(BxCouple_):iVar_V(BzCouple_))

  ! Transform vector variables from SC to IH
  if(TypeCoordSource /= TypeCoordSystem)then
     State_V(Ux_:Uz_) = transform_velocity(Time_Simulation,&
          State_V(Ux_:Uz_), XyzSource_D * No2Si_V(UnitX_), &
          TypeCoordSource, TypeCoordSystem)
     if(DoCoupleVar_V(Bfield_)) &
          State_V(Bx_:Bz_) = matmul( State_V(Bx_:Bz_), SourceTarget_DD)
  end if

  !Convert from SI units to normalized units
  State_V(rho_)      = State_V(rho_)   *Si2No_V(UnitRho_)
  State_V(Ux_:Uz_)   = State_V(Ux_:Uz_)*Si2No_V(UnitU_)
  if(DoCoupleVar_V(Bfield_)) &
       State_V(Bx_:Bz_)   = State_V(Bx_:Bz_)*Si2No_V(UnitB_)

  if(DoCoupleVar_V(Wave_))State_V(WaveFirst_:WaveLast_) = &
       Buffer_V(iVar_V(WaveFirstCouple_):iVar_V(WaveLastCouple_)) &
       *Si2No_V(UnitEnergyDens_)

  State_V(p_)  = Buffer_V(iVar_V(PCouple_))*Si2No_V(UnitP_)
  if(DoCoupleVar_V(ElectronPressure_))then
     State_V(Pe_) = Buffer_V(iVar_V(PeCouple_))*Si2No_V(UnitP_)
  else if(UseElectronPressure)then
     State_V(Pe_) = 0.5*State_V(p_)
     State_V(p_)  = State_V(Pe_)
  end if

  if(DoCoupleVar_V(AnisoPressure_))then
     State_V(Ppar_) = Buffer_V(iVar_V(PparCouple_))*Si2No_V(UnitP_)
  else if(UseAnisoPressure)then
     State_V(Ppar_) = Buffer_V(iVar_V(PCouple_))*Si2No_V(UnitP_)
  end if

  if( .not. DoCoupleVar_V(MultiFluid_)  .and. nFluid > 1 .or. &
      .not. DoCoupleVar_V(MultiSpecie_) .and. UseMultiSpecies)then
     ! Values for neutrals / ions should be prescribed in set_BCs.f90
     IsFullyCoupledFluid = .false.
  else
     IsFullyCoupledFluid = .true.
  end if
end subroutine get_from_spher_buffer_grid
          
  
