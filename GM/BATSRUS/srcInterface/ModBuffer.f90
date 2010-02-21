module ModBuffer
  use ModMain,     ONLY: nPhiBuff, nThetaBuff, RBuffMin, RBuffMax
  use ModNumConst, ONLY: cPi, cTwoPi
  use CON_global_vector,   ONLY: ubound_vector, point_state_v
  use CON_grid_descriptor, ONLY: GridDescriptorType, Nodes_, &
       bilinear_interpolation, DomainDecompositionType, is_proc0, &
       init_decomposition, get_root_decomposition, complete_grid, &
       set_standard_grid_descriptor, bcast_decomposition
  use CON_coupler,         ONLY: SC_
  implicit none
  save

  !This is a name of a global vector, which
  !is in the global vector storage
  character(LEN=10)::NameBuffer


  type(DomainDecompositionType),target::&
       LocalBufferDD
  type(GridDescriptorType)::LocalBufferGD

  integer:: SourceID_ = SC_

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
  use ModAdvance,    ONLY: Rho_, RhoUx_, RhoUz_, Ux_, Uz_, Bx_, Bz_, p_
  use ModVarIndexes, ONLY: WaveFirst_, WaveLast_, Pe_
  use CON_coupler,   ONLY: Grid_C
  use CON_axes,      ONLY: transform_matrix, transform_velocity
  use ModPhysics,    ONLY: No2Si_V,Si2No_V,UnitRho_,UnitU_,UnitB_,UnitP_,UnitX_
  use ModPhysics,    ONLY: UnitEnergyDens_
  implicit none
  integer,intent(in)::nVar
  real,intent(in)::XyzTarget_D(nDim)
  real,dimension(nVar),intent(out)::State_V
  real,dimension(nDim)::Sph_D

  character (len=3) :: TypeCoordSource
  real, save        :: SourceTarget_DD(nDim, nDim)
  real              :: TimeSimulationLast = -1.0
  real              :: XyzSource_D(nDim)

  ! logicals and buffer indices for waves and electron pressure
  logical :: UseWave, UseElectronPressure
  integer :: iBuffWaveFirst, iBuffWaveLast, iBuffPe

  !Misc
  integer:: UBound_I(2)
  !---------------------------------------------------------------------------

  !Set nVarBuff
  if(DoInit)then
     DoInit = .false.
     UBound_I=ubound_vector(NameBuffer)
     nVarBuff = UBound_I(1)
     allocate(Buffer_V(nVarBuff))
  end if

  UseWave = nVar>=10
  iBuffWaveFirst = 8
  iBuffWaveLast  = 8 + WaveLast_ - WaveFirst_

  UseElectronPressure = nVar/=2*(nVar/2)
  iBuffPe = nVarBuff - 1

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
       nVarBuff,&
       nDim,    &
       Sph_D,   &
       LocalBufferGD,&
       bilinear_interpolation)

  State_V(Rho_) = Buffer_V(Rho_)
  !Transform to primitive variables
  State_V(Ux_:Uz_) = Buffer_V(rhoUx_:rhoUz_)/Buffer_V(rho_)
  State_V(Bx_:Bz_) = Buffer_V(Bx_:Bz_)
  State_V(P_)      = Buffer_V(nVarBuff)

  ! Transform vector variables from SC to IH
  if(TypeCoordSource /= TypeCoordSystem)then
     State_V(Ux_:Uz_) = transform_velocity(Time_Simulation,&
          State_V(Ux_:Uz_), XyzSource_D * No2Si_V(UnitX_), &
          TypeCoordSource, TypeCoordSystem)
     State_V(Bx_:Bz_) = matmul( State_V(Bx_:Bz_), SourceTarget_DD)
  end if

  !Convert from SI units to normalized units
  State_V(rho_)      = State_V(rho_)   *Si2No_V(UnitRho_)
  State_V(Ux_:Uz_)   = State_V(Ux_:Uz_)*Si2No_V(UnitU_)
  State_V(Bx_:Bz_)   = State_V(Bx_:Bz_)*Si2No_V(UnitB_)
  State_V(P_)        = State_V(P_)     *Si2No_V(UnitP_)
  if(UseWave) State_V(WaveFirst_:WaveLast_) = &
       Buffer_V(iBuffWaveFirst:iBuffWaveLast)*Si2No_V(UnitEnergyDens_)
  if(UseElectronPressure) State_V(Pe_) = Buffer_V(iBuffPe)*Si2No_V(UnitP_)

end subroutine get_from_spher_buffer_grid
          
  
