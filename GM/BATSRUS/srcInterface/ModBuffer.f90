module ModBuffer
  use ModMain,ONLY:&
       nPhiBuff,nThetaBuff,RBuffMin,RBuffMax
  use CON_global_vector
  use CON_grid_descriptor
  use CON_coupler,   ONLY: SC_
  implicit none
  save

  !This is a name of a global vector, which
  !is in the global vector storage
  character(LEN=10)::NameBuffer


  type(DomainDecompositionType),target::&
       LocalBufferDD
  type(GridDescriptorType)::LocalBufferGD

  integer:: TargetID_ = SC_
  logical::DoInit
contains
  subroutine set_buffer_name(NameIn,IDIn)
    character(LEN=*),intent(in)::NameIn
    integer,optional,intent(in)::IDIn
    !-------------

    NameBuffer=NameIn
    if(present(IDIn))TargetID_ = IDIn
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
         XyzMin_D=(/RBuffMin,cZero,cZero/),&
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

subroutine get_from_spher_buffer_grid(Xyz_D,nVar,State_V)
  use ModBuffer
  use ModMain,       ONLY: nDim, R_, Phi_, Theta_, x_, y_, z_,&
       TypeCoordSystem, Time_Simulation
  use ModAdvance,    ONLY: Rho_, RhoUx_, RhoUz_, Ux_, Uz_, Bx_, Bz_, p_
  use CON_coupler,   ONLY: Grid_C
  use CON_axes,      ONLY: transform_matrix, transform_velocity
  use ModPhysics,    ONLY: No2Si_V,Si2No_V,UnitRho_,UnitU_,UnitB_,UnitP_,UnitX_

  implicit none
  integer,intent(in)::nVar
  real,intent(in)::Xyz_D(nDim)
  real,dimension(nVar),intent(out)::State_V
  real,dimension(nDim)::Sph_D

  character (len=3) :: TypeCoordSc
  real, save        :: ScIh_DD(nDim, nDim)
  real              :: TimeSimulationLast = -1.0
  real              :: XyzSc_D(nDim)
  !---------------------------------------------------------------------------

  TypeCoordSc = Grid_C(TargetID_) % TypeCoord

  if(TypeCoordSc /= TypeCoordSystem) then
     ! Convert target coordinates to the coordiante system of the model

     if(Time_Simulation > TimeSimulationLast)then
        ScIh_DD = &
             transform_matrix(Time_Simulation, TypeCoordSystem, TypeCoordSc)
        TimeSimulationLast = Time_Simulation
     end if
     XyzSc_D = matmul(ScIh_DD, Xyz_D)
  else
     XyzSc_D = Xyz_D
  end if

  ! Convert to left handed spherical coordinates !!!
  call xyz_to_spherical(XyzSc_D(x_),XyzSc_D(y_),XyzSc_D(z_),&
       Sph_D(R_),Sph_D(Phi_),Sph_D(Theta_))

  ! Get the target state from the spherical buffer grid
  State_V=point_state_v(&
       NameBuffer,&
       nVar,    &
       nDim,    &
       Sph_D,   &
       LocalBufferGD,&
       bilinear_interpolation)

  !Transform to primitive variables
  State_V(Ux_:Uz_) = State_V(rhoUx_:rhoUz_)/State_V(rho_)

  ! Transform vector variables from SC to IH
  if(TypeCoordSc /= TypeCoordSystem)then
     State_V(Ux_:Uz_) = transform_velocity(Time_Simulation,&
          State_V(Ux_:Uz_), XyzSc_D * No2Si_V(UnitX_), &
          TypeCoordSc, TypeCoordSystem)
     State_V(Bx_:Bz_) = matmul( State_V(Bx_:Bz_), ScIh_DD)
  end if

  !Convert from SI units to normalized units
  State_V(rho_)      = State_V(rho_)   *Si2No_V(UnitRho_)
  State_V(Ux_:Uz_)   = State_V(Ux_:Uz_)*Si2No_V(UnitU_)
  State_V(Bx_:Bz_)   = State_V(Bx_:Bz_)*Si2No_V(UnitB_)
  State_V(P_)        = State_V(P_)     *Si2No_V(UnitP_)

end subroutine get_from_spher_buffer_grid
          
  
