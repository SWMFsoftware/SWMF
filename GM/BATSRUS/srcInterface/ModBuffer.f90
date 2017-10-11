!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModBuffer
  use ModMain,     ONLY: nPhiBuff, nThetaBuff, BufferMin_D, BufferMax_D
  use ModNumConst, ONLY: cPi, cTwoPi
  use CON_global_vector,   ONLY: point_state_v
  use CON_grid_descriptor, ONLY: GridDescriptorType, Nodes_, &
       bilinear_interpolation, DomainDecompositionType, is_proc0, &
       init_decomposition, get_root_decomposition, complete_grid, &
       set_standard_grid_descriptor, bcast_decomposition
  use CON_coupler,         ONLY: SC_, IH_, nVarIndexCouple, nCoupleVarGroup
  implicit none
  save

  !This is a name of a global vector, which
  !is in the global vector storage
  character(LEN=10)::NameBuffer

  type(DomainDecompositionType),target:: LocalBufferDD
  type(GridDescriptorType)::LocalBufferGD

  integer:: SourceID_ = SC_, TargetID_ = IH_

  integer::nVarBuff=-1
  logical::DoInit =.true.
  real,dimension(:),allocatable::Buffer_V

  integer, public:: nVarCouple
  integer, public:: iVar_V(nVarIndexCouple)
  logical, public:: DoCoupleVar_V(nCoupleVarGroup)  

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
         XyzMin_D= BufferMin_D,  &
         XyzMax_D= BufferMax_D, &
         nCells_D=(/1,nPhiBuff,nThetaBuff/),&
         PE_I=(/0/))
    if(.not.IsLocal)then
       call bcast_decomposition(DD)
    else
       call complete_grid(DD)
       call set_standard_grid_descriptor(DD, &
            nGhostGridPoints=1,  &
            iStandard=Nodes_,&
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
       UseMultiSpecies, SignB_, Ehot_
  use CON_coupler,   ONLY: &
       Grid_C, &
       Bfield_, ElectronPressure_, AnisoPressure_, Wave_,&
       MultiFluid_, MultiSpecie_, CollisionlessHeatFlux_, &
       RhoCouple_, RhoUxCouple_, RhoUzCouple_, PCouple_, &
       BxCouple_, BzCouple_, PeCouple_, PparCouple_, &
       WaveFirstCouple_, WaveLastCouple_, EhotCouple_

  use CON_axes,      ONLY: transform_matrix, transform_velocity
  use ModPhysics,    ONLY: No2Si_V,Si2No_V,UnitRho_,UnitU_,UnitB_,UnitP_,UnitX_
  use ModPhysics,    ONLY: UnitEnergyDens_
  use ModMultiFluid, ONLY: IsFullyCoupledFluid
  implicit none

  integer,intent(in) :: nVar
  real,   intent(in) :: XyzTarget_D(MaxDim)
  real,   intent(out):: State_V(nVar)

  real             :: Sph_D(MaxDim)

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
        SourceTarget_DD = transform_matrix(&
             Time_Simulation, TypeCoordSystem, TypeCoordSource)
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
  if(.true.)then
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

  if(DoCoupleVar_V(CollisionlessHeatFlux_))then
     State_V(Ehot_) = Buffer_V(iVar_V(EhotCouple_))*Si2No_V(UnitEnergyDens_)
  endif

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
!==============================
subroutine plot_buffer(iFile)
  use ModPlotFile, ONLY: save_plot_file
  use ModNumConst,   ONLY: cDegToRad
  use ModAdvance,    ONLY: &
       UseElectronPressure, UseAnisoPressure
  use ModVarIndexes, ONLY: &
       nVar,Rho_, Ux_, Uz_, Bx_, Bz_, p_, &
       WaveFirst_, WaveLast_, Pe_, Ppar_, nFluid, &
       UseMultiSpecies, SignB_, Ehot_
  use ModIO,            ONLY: NamePrimitiveVarOrig
  use ModTimeConvert,   ONLY: time_int_to_real, time_real_to_int
  use ModMain,          ONLY: StartTime, Time_Simulation, x_, y_, z_, n_step
  use ModMain, ONLY: nPhiBuff, nThetaBuff, BufferMin_D, BufferMax_D, BuffR_
  use ModPhysics,    ONLY: No2Si_V, UnitRho_, UnitU_, UnitB_, UnitP_, UnitX_,&
       UnitEnergyDens_
  use ModProcMH,     ONLY: iProc
  implicit none
  integer, intent(in)::iFile
  integer:: iTimePlot_I(7), iR, iPhi, iTheta
  real   :: R, Theta, Phi, CosTheta, SinTheta
  real, allocatable,dimension(:,:,:):: State_VII, Coord_DII
  character(LEN=30)::NameFile
  !-----------------------
  if(iProc/=0)RETURN !May be improved.
  allocate(State_VII(3 + nVar, 0:nPhiBuff, 0:nThetaBuff))
  allocate(Coord_DII(2, 0:nPhiBuff, 0:nThetaBuff))
  call time_real_to_int(StartTime + Time_Simulation, iTimePlot_I)
  do iR = 1,2
     R = BufferMin_D(BuffR_)*(2 - iR) + BufferMax_D(BuffR_)*(iR - 1)
     write(NameFile,'(a,i2.2,a,i4.4,a,5(i2.2,a))')'R=',nint(R),'Rs_',&
          iTimePlot_I(1),'_',iTimePlot_I(2),'_',iTimePlot_I(3),'_',&
          iTimePlot_I(4),'_',iTimePlot_I(5),'_',iTimePlot_I(6),'.out'
     do iTheta = 0, nThetaBuff
        Theta = -89.999 + (179.999/nThetaBuff)*iTheta
        Coord_DII(2,:,iTheta) = Theta
        CosTheta = cos(Theta*cDegToRad)
        SinTheta = sin(Theta*cDegToRad)
        do iPhi = 0, nPhiBuff
           Phi = (360.0/nPhiBuff)*iPhi
           Coord_DII( 1,iPhi,iTheta) = Phi
           State_VII(x_,iPhi,iTheta) = R*CosTheta*cos(Phi*cDegToRad)
           State_VII(y_,iPhi,iTheta) = R*CosTheta*sin(Phi*cDegToRad)
           State_VII(z_,iPhi,iTheta) = R*SinTheta
           call get_from_spher_buffer_grid(&
                State_VII(x_:z_,iPhi,iTheta), nVar,&
                State_VII(z_+1:z_+nVar,iPhi,iTheta))

        end do
     end do
     !Convert from normalized units to SI
     State_VII(  x_:z_,:,:) = State_VII(  x_:z_,:,:)*No2Si_V(UnitX_)
     State_VII(z_+rho_,:,:) = State_VII(z_+rho_,:,:)*No2Si_V(UnitRho_)
     State_VII(z_+Ux_:z_+Uz_,:,:)   = State_VII(z_+Ux_:z_+Uz_,:,:)*&
          No2Si_V(UnitU_)
     State_VII(z_+Bx_:z_+Bz_,:,:)   = State_VII(z_+Bx_:z_+Bz_,:,:)*&
          No2Si_V(UnitB_)

     if(WaveFirst_ > 1)State_VII(z_+WaveFirst_:z_+WaveLast_,:,:) = &
          State_VII(z_+WaveFirst_:z_+WaveLast_,:,:)* &
          No2Si_V(UnitEnergyDens_)

     State_VII(z_+p_,:,:)  = State_VII(z_+p_,:,:)*No2Si_V(UnitP_)
     if(UseElectronPressure)State_VII(z_+Pe_,:,:)  = &
          State_VII(z_+Pe_,:,:)*No2Si_V(UnitP_)

     if(UseAnisoPressure)State_VII(z_+Ppar_,:,:)  = &
          State_VII(z_+Ppar_,:,:)*No2Si_V(UnitP_)

     if(Ehot_>1)State_VII(z_+Ehot_,:,:) = &
          State_VII(z_+Ehot_,:,:)*No2Si_V(UnitEnergyDens_)

     call save_plot_file(NameFile,&
          StringHeaderIn=&
          'SC-IH interface, longitude and latitude are in deg, other in SI',&
          NameVarIn    = &
          'Long Lat x y z '//NamePrimitiveVarOrig//' R',&
          nDimIn=2,      &
          nStepIn=n_step, TimeIn=Time_Simulation,&
          ParamIn_I=(/R*No2Si_V(UnitX_)/), &
          CoordIn_DII=Coord_DII, &
          VarIn_VII=State_VII)
  end do
end subroutine plot_buffer
  
