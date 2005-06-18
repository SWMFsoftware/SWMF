module ModBuffer
  use ModMain,ONLY:&
       nPhiBuff,nThetaBuff,RBuffMin,RBuffMax
  use CON_global_vector
  use CON_grid_descriptor
  use ModNumConst,ONLY:cUnit_DD
  implicit none
  save
  character(LEN=10)::NameBuffer
  type(DomainDecompositionType),target::&
       LocalBufferDD
  type(GridDescriptorType)::LocalBufferGD
  logical::DoInit
  logical::UseRotatingBufferGrid=.false.
  real,dimension(3,3)::BuffToMh_DD
contains
  subroutine set_buffer_name(NameIn)
    character(LEN=*),intent(in)::NameIn
    NameBuffer=NameIn
  end subroutine set_buffer_name
  !====================================================
  subroutine set_spher_buffer_grid(DD,CompID_,IsLocal)
    type(DomainDecompositionType),&
         intent(out)::DD
    integer,intent(in)::CompID_
    logical,intent(in)::IsLocal
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
    BuffToMh_DD=cUnit_DD
  end subroutine set_spher_buffer_grid
end module ModBuffer
!=============================================================!
subroutine set_rotate_buffer_grid(Time)
  use ModBuffer
  use ModMain,       ONLY: nDim,TypeBc_I,&
       TypeCoordSystem, Time_Simulation,UseRotatingBc
  use CON_coupler,   ONLY: SC_, Grid_C
  use CON_axes,      ONLY: transform_matrix
  use ModNumConst,ONLY:cUnit_DD
  implicit none
  real,intent(in)::Time
  real:: TimeLast=-1.0
  if(any(TypeBc_I=='coronatoih').and.UseRotatingBc.and.Time/=TimeLast)then
    BuffToMh_DD = &
          transform_matrix(Time,Grid_C(SC_) % TypeCoord,TypeCoordSystem)
    TimeLast=Time
    UseRotatingBufferGrid=.true.
  else
      TimeLast=Time
      UseRotatingBufferGrid=.false.
      BuffToMh_DD=cUnit_DD
  end if
  

end subroutine set_rotate_buffer_grid

subroutine get_from_spher_buffer_grid(XyzMh_D,nVar,State_V)
  use ModBuffer
  use ModMain,       ONLY: nDim, R_, Phi_, Theta_, x_, y_, z_
  use ModAdvance,    ONLY: Rho_, RhoUx_, RhoUz_, Ux_, Uz_, Bx_, Bz_, p_
  use ModPhysics,    ONLY: UnitSI_rho,UnitSI_U,UnitSI_B,UnitSI_p,UnitSI_X

  implicit none
  integer,intent(in)::nVar
  real,intent(in)::XyzMh_D(nDim)
  real,dimension(nVar),intent(out)::State_V
  real,dimension(nDim)::Sph_D,XyzBuff_D

  !---------------------------------------------------------------------------
  if(UseRotatingBufferGrid)then
     XyzBuff_D=matmul(XyzMh_D,BuffToMh_DD)
  else
     XyzBuff_D=XyzMh_D
  end if
  ! Convert to left handed spherical coordinates !!!
  call xyz_to_spherical(XyzBuff_D(x_),XyzBuff_D(y_),XyzBuff_D(z_),&
       Sph_D(R_),Sph_D(Phi_),Sph_D(Theta_))

  ! Get the state from the spherical buffer grid
  State_V=point_state_v(&
       NameBuffer,&
       nVar,    &
       nDim,    &
       Sph_D,   &
       LocalBufferGD,&
       bilinear_interpolation)

  !Transform to primitive variables
  State_V(Ux_:Uz_) = State_V(rhoUx_:rhoUz_)/State_V(rho_)

  !Convert from SI units to normalized units
  State_V(rho_)      = State_V(rho_)   /UnitSI_rho
  State_V(Ux_:Uz_)   = State_V(Ux_:Uz_)/UnitSI_U
  State_V(Bx_:Bz_)   = State_V(Bx_:Bz_)/UnitSI_B
  State_V(P_)        = State_V(P_)     /UnitSI_P
  
  if(UseRotatingBufferGrid)then
    State_V(Ux_:Uz_)= matmul(BuffToMh_DD,State_V(Ux_:Uz_))
    State_V(Bx_:Bz_)= matmul(BuffToMh_DD,State_V(Bx_:Bz_))
 end if

end subroutine get_from_spher_buffer_grid
          
  
