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
  !==========================================================================
  subroutine set_buffer_name(NameIn,IDIn)
    character(LEN=*),intent(in)::NameIn
    integer,optional,intent(in)::IDIn
    !------------------------------------------------------------------------

    NameBuffer=NameIn
    if(present(IDIn))SourceID_ = IDIn
  end subroutine set_buffer_name
  !==========================================================================
  subroutine set_spher_buffer_grid(DD,CompID_,IsLocal)
    type(DomainDecompositionType),&
         intent(out)::DD
    integer,intent(in)::CompID_
    logical,intent(in)::IsLocal
    !-----------------------------------------------------------------------
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

!==============================================================================

subroutine get_from_spher_buffer_grid(XyzTarget_D,nVar,State_V)
  use ModBuffer
  use ModMain,       ONLY: &
       MaxDim, x_, y_, z_,&
       TypeCoordSystem, Time_Simulation, DoThinCurrentSheet
  use ModAdvance,    ONLY: &
       UseElectronPressure, UseAnisoPressure
  use ModVarIndexes, ONLY: &
       Rho_, Ux_, Uz_, Bx_, Bz_, p_, &
       WaveFirst_, WaveLast_, Pe_, Ppar_, nFluid, &
       UseMultiSpecies, SignB_
  use CON_coupler,   ONLY: &
       Grid_C, DoCoupleVar_V, iVar_V, nVarCouple,&
       Bfield_, ElectronPressure_, AnisoPressure_, Wave_,&
       MultiFluid_, MultiSpecie_, &
       RhoCouple_, RhoUxCouple_, RhoUzCouple_, PCouple_, &
       BxCouple_, BzCouple_, PeCouple_, PparCouple_, &
       WaveFirstCouple_, WaveLastCouple_, &
       UseGlobalMpiCoupler
  use CON_axes,      ONLY: transform_matrix, transform_velocity
  use ModPhysics,    ONLY: No2Si_V,Si2No_V,UnitRho_,UnitU_,UnitB_,UnitP_,UnitX_
  use ModPhysics,    ONLY: UnitEnergyDens_
  use ModMultiFluid, ONLY: IsFullyCoupledFluid
  implicit none

  integer,intent(in)              :: nVar
  real,intent(in)                 :: XyzTarget_D(MaxDim)
  real,dimension(nVar),intent(out)::State_V
  real,dimension(MaxDim)            ::Sph_D

  character (len=3) :: TypeCoordSource
  real, save        :: SourceTarget_DD(MaxDim, MaxDim)
  real              :: TimeSimulationLast = -1.0
  real              :: XyzSource_D(MaxDim)

  !---------------------------------------------------------------------------
  if(DoInit)then
     DoInit = .false.
     allocate(Buffer_V(nVarCouple))
  end if
  State_V = 0.0
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
       Sph_D(1),Sph_D(2),Sph_D(3))

  ! Get the target state from the spherical buffer grid
  if(UseGlobalMpiCoupler) then       
     call interpolate_from_global_buffer(Sph_D, nVarCouple, Buffer_V)
  else
     Buffer_V=point_state_v(&
          NameBuffer,&
          nVarCouple,&
          MaxDim,    &
          Sph_D,     &
          LocalBufferGD,&
          bilinear_interpolation)
  end if

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

  if(SignB_>1)then
     if(DoThinCurrentSheet)then
        ! In both IH and OH we have no B0, so we ignore that !         
        if(sum(State_V(Bx_:Bz_)*XyzTarget_D) < 0.0)then
           State_V(Bx_:Bz_) = -State_V(Bx_:Bz_)
           State_V(SignB_)=-1.0
        else
           State_V(SignB_)= 1.0
        end if
     else
        State_V(SignB_) = 0.0
     end if
  end if

end subroutine get_from_spher_buffer_grid

!============================================================================
subroutine interpolate_from_global_buffer(SphSource_D, nVar, Buffer_V)

  ! DESCRIPTION
  ! This subroutine is used to interpolate from  state variables defined on a
  ! spherical buffer grid into the input point SphSource_D.
  ! The buffer grid overlaps some part of the computational grid of a 
  ! source component that is coupled to this component.
  ! The buffer grid  has the same coordinate system as the source component
  ! (but may have a different grid resolution).
  ! It is assumed that the buffer grid was filled with the state vector from 
  ! the source component at some earlier stage.

  ! INPUT:
  ! SphSource_D is associated with a point in the target component, and it
  ! is assumed that is was already converted to the source coordinate system.
 
  ! nVar is the number of state variables used in coupling the two components.

  ! Implicit inputs to this subroutine are the buffer grid size, points 
  ! and the state vector at each point (USEd from BATSRUS).

  ! OUTPUT:
  ! Buffer_V is the state vector resulting from the interpolation.

  use ModInterpolate, ONLY: trilinear
  use ModMain,        ONLY: BufferState_VG, BufferMin_D,&
                            nRBuff, nPhiBuff, nThetaBuff, dSphBuff_D
 
  implicit none

  ! Input and output variables
  real,intent(in)    :: SphSource_D(3)
  integer,intent(in) :: nVar
  real,intent(out)   :: Buffer_V(nVar)

  real    :: NormSph_D(3)
  ! logical :: DoTest, DoTestMe

  character(len=*), parameter :: NameSub = 'interpolate_from_global_buffer'
  !-------------------------------------------------------------------------
  !  call CON_set_do_test(NameSub,DoTest, DoTestMe)

  ! Convert to normalized coordinates. 
  ! Radial is node centered, theta and phi are cell centered.
  NormSph_D = (SphSource_D - BufferMin_D)/dSphBuff_D + (/ 1.0, 0.5, 0.5 /)

  Buffer_V = trilinear(BufferState_VG, nVar, 1, nRBuff,0, nPhiBuff+1, &
       0, nThetaBuff+1, NormSph_D, DoExtrapolate=.true.)

end subroutine interpolate_from_global_buffer
          
  
