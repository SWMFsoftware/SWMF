!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE SP

module CON_couple_mh_sp

  ! This coupler employs the following global arrays.
  ! SP_Xyz_DI - is the array of the Lagrangian points.
  ! A part of this array, as well as the mask 'SP_IsInIH' is only availbale
  ! at the PE set, at which either IH or SP run. 
  ! This coordinates are expressed in terms of
  ! the length units and with respect to the frame of reference defined in IH.
  ! Another part of this array, as well as the mask 'SP_IsInSC' 
  ! is only available at the PE set, at which either SC or SP run. 
  ! This coordinates are expressed in terms of
  ! the length units and with respect to the frame of reference defined in SC.
  ! SP_XyzSP - in the array of all Lagrangian points, in units and in the 
  ! frame of reference defined at SP, is available only at the PE set, 
  ! at which SP runs.

  use CON_coupler

  use IH_wrapper, ONLY: IH_synchronize_refinement, &        !^CMP IF IH
       IH_get_line, IH_get_for_sp, IH_get_a_line_point      !^CMP IF IH

  use SC_wrapper, ONLY: SC_synchronize_refinement, &        !^CMP IF SC
       SC_get_line, SC_get_for_sp, SC_get_a_line_point      !^CMP IF SC

  use CON_global_message_pass
  use CON_axes

  use SP_wrapper, ONLY: &
       SP_put_from_mh, SP_put_input_time, &
       SP_put_line, SP_request_line, SP_get_grid_descriptor_param, &
       SP_get_line_all

  implicit none

  private !Except
  public::couple_mh_sp_init
  public::couple_ih_sp              !^CMP IF IH
  public::couple_sc_sp              !^CMP IF SC

  type(GridDescriptorType),save::SP_GridDescriptor !Target

  type(GridDescriptorType),save::IH_GridDescriptor !Source  !^CMP IF IH
  type(RouterType),save,private::RouterIhSp                 !^CMP IF IH

  type(GridDescriptorType),save::SC_GridDescriptor !Source  !^CMP IF SC
  type(RouterType),save,private::RouterScSp                 !^CMP IF SC

  logical,save::DoInit=.true.
  real,allocatable,dimension(:,:)::XyzTemp_DI

  real,dimension(:,:),pointer ::Xyz_DI
  logical,dimension(:),pointer :: Is_I
  integer,parameter::nPointMax=5000
  integer::nPoint=0
  integer::iPoint
  integer::iError
  real::bDxyz_I(1:6)!The interpolated values of full B and DXyz
  real::DsResolution,XyzLine_D(3)
  real,save::rBoundIh=21.0                !^CMP IF IH
  real,save::rBoundSc=1.20                !^CMP IF SC
  logical::DoTest,DoTestMe
  character(LEN=*),parameter::NameSub='couple_mh_sp'
  real,dimension(3,3)::ScToIh_DD,ScToSp_DD,IhToSp_DD,SpToSc_DD
  real :: tNow

contains
  !==================================================================
  subroutine couple_mh_sp_init
    use CON_physics, ONLY: get_time
    use ModConst

    logical::DoneRestart
    integer::nLine

    integer:: iGridMin_D(3), iGridMax_D(3), ierror
    real:: Disp_D(3)
    real, pointer:: CoordMisc_DI(:,:)

    ! available directions of interface between SP and MH 
    ! (see subroutine exchange_lines below)
    integer, parameter:: &
         iInterfaceBegin = -1, iInterfaceOrigin = 0, iInterfaceEnd = 1

    character(len=*), parameter:: NameSub = 'couple_mh_sp_init'
    !----------------------------------------------------------------------
    if(.not.DoInit)return
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    DoInit=.false.
    !The initialization can be done only once

    call get_time(tSimulationOut=tNow)

    !\
    ! Set grid descriptors for components
    ! Initialize routers
    !/
    call SP_get_grid_descriptor_param(iGridMin_D, iGridMax_D, Disp_D)
    call set_grid_descriptor_id(SP_,&
         nDim = 3, &
         iGridPointMin_D = iGridMin_D, &
         iGridPointMax_D = iGridMax_D, &
         Displacement_D  = Disp_D, &
         GridDescriptor  = SP_GridDescriptor)

    if(use_comp(SC_))then  
       ! Set pair SC-SP
       call set_standard_grid_descriptor(SC_,GridDescriptor=&
            SC_GridDescriptor)
       call init_router(SC_GridDescriptor,SP_GridDescriptor,&
            RouterScSp)
    end if

    if(use_comp(IH_))then
       ! Set pair IH-SP
       call set_standard_grid_descriptor(IH_,GridDescriptor=&
            IH_GridDescriptor)
       call init_router(IH_GridDescriptor,SP_GridDescriptor,&
            RouterIhSp)
    end if


    !\
    ! Allocate storage for field line requests
    !/
    nLine = SP_GridDescriptor%DD%Ptr%nBlockAll
    call allocate_vector('SP_Xyz_DI',&
         SP_GridDescriptor%DD%Ptr%nDim, nLine)

    !\
    ! Extract and exchange initial data
    !/
    call exchange_lines(iInterfaceOrigin, SC_, SC_GridDescriptor)
    call exchange_lines(iInterfaceEnd,    IH_, IH_GridDescriptor)

    ! reserve memeory for SP grid
    call deallocate_vector('SP_Xyz_DI')
    call allocate_vector('SP_Xyz_DI', &
                  SP_GridDescriptor%DD%Ptr%nDim, &
                  product(SP_GridDescriptor%DD%Ptr%nCells_D)*nLine)
    call associate_with_global_vector(CoordMisc_DI, 'SP_Xyz_DI')
    if(is_proc(SP_))&
         call SP_get_line_all(CoordMisc_DI)
    nullify(CoordMisc_DI)
    call bcast_global_vector('SP_Xyz_DI',i_proc0(SP_),i_comm())

    if(use_comp(SC_))then
       if(RouterScSp%IsProc)then
          call allocate_mask('SP_IsInSC','SP_Xyz_DI')
          call set_mask('SP_IsInSC','SP_Xyz_DI',is_in_sc)
       end if
    end if

    if(use_comp(IH_))then
       if(RouterIhSp%IsProc)then
          call allocate_mask('SP_IsInIH','SP_Xyz_DI')
          call set_mask('SP_IsInIH','SP_Xyz_DI',is_in_ih)
       end if
    end if
    
  contains
    !================================================================    
    subroutine exchange_lines(iInterfaceType, iMH, MH_GridDescriptor)
      ! MH extracts and sends field lines requested by SP;
      !----------------------------------------------------------------
      ! direction interface between SP and MH component on field lines:
      ! -1 -> the beginning of lines
      !  0 -> origin points of lines as defined at SP
      !  1 -> the end of lines
      integer, intent(in):: iInterfaceType
      !----------------------------------------------------------------
      ! index of MH component
      integer, intent(in):: iMH
      !----------------------------------------------------------------
      ! grid descriptor of MH component
      type(GridDescriptorType), intent(in)::MH_GridDescriptor

      ! conversion matrix between SP and MH coordinates
      real:: Convert_DD(3,3)
      ! request coordinates, one per line
      real, pointer:: CoordMisc_DI(:,:)
      ! number of particles per line
      integer,allocatable:: nParticleAtLine_I(:)
      ! requested variables
      integer:: nVar
      character(len=100):: NameVar
      ! particle data
      real, allocatable:: Particle_II(:,:)
      ! field line index in particle data
      integer:: iFLIndex
      ! MPI
      integer:: MH_iProcFrom, SP_iProcTo
      integer:: iProcTo_I(1), iProcFrom_I(1)
      integer:: nParticleThisProc
      integer, allocatable:: nParticleRecv_I(:), nParticleSend_I(:)
      integer, allocatable:: iStatus_II(:,:), iRequestS_I(:), iRequestR_I(:)
      integer:: nRequestS, nRequestR
      integer:: iTag = 0
      real, allocatable:: BuffRecv_I(:), BuffSend_I(:)
      ! loop variables
      integer:: iLine, iBuff, iParticle
      !----------------------------------------------------------------
      ! allocate arrays for non-blocking communcations
      allocate(nParticleRecv_I(0:n_proc()-1))
      allocate(nParticleSend_I(0:n_proc()-1))
      allocate(iRequestS_I(n_proc()))
      allocate(iRequestR_I(n_proc()))
      allocate(iStatus_II(MPI_STATUS_SIZE, n_proc()))
      !\
      ! SP requests field lines from MH specifying:
      ! - direction of interface (begin/origin/end)variables
      ! - number & names of variables
      !/
      if(is_proc(SP_))then
         call associate_with_global_vector(CoordMisc_DI,'SP_Xyz_DI')
         call SP_request_line(NameVar, nVar, iInterfaceType, CoordMisc_DI)
         nullify(CoordMisc_DI)
      end if
      call bcast_global_vector('SP_Xyz_DI',i_proc0(SP_),i_comm())
      call MPI_bcast(NameVar, len(NameVar), MPI_CHARACTER, &
           i_proc0(SP_), i_comm(), iError)
      call MPI_bcast(nVar, 1, MPI_INTEGER, i_proc0(SP_), i_comm(), iError)
      !\
      ! determine which index corresponds to field line index
      !/
      iFLIndex = index(NameVar, 'fl')/3 + 1
      if(iFLIndex == 0) call CON_stop(NameSub//&
           ': set of variables requested is not sufficient for coupling')
      !\
      ! Extract field lines at MH and send to SP
      !/
      nRequestS = 0
      nRequestR = 0
      if(is_proc(iMH))then
         ! get request locations
         call associate_with_global_vector(CoordMisc_DI,'SP_Xyz_DI')
         Convert_DD = transform_matrix(tNow,&
              Grid_C(SP_)%TypeCoord, Grid_C(iMH)%TypeCoord)

         ! convert locations to MH coordinates
         do iLine = 1, nLine
            CoordMisc_DI(:,iLine) = matmul(Convert_DD, CoordMisc_DI(:,iLine))
         end do
         ! extract field line data
         allocate(nParticleAtLine_I(nLine))
         select case(iMH)
         case(SC_)
            call SC_get_line(nLine,CoordMisc_DI, iInterfaceType,nVar,NameVar, &
                 nParticleThisProc, Particle_II)
         case(IH_)
            call IH_get_line(nLine,CoordMisc_DI, iInterfaceType,nVar,NameVar, &
                 nParticleThisProc, Particle_II)
         case default
            call CON_stop(NameSub//': incorrect MH component')
         end select
         nullify(CoordMisc_DI)
         ! count number of particles per field line
         do iLine = 1, nLine
            nParticleAtLine_I(iLine)=&                 
                 count(nint(Particle_II(iFLIndex,:))==iLine)
         end do

         !\
         ! Send # of particles to SP: 
         ! on MH by index fl_ can find recepient on SP 
         !/
         do SP_iProcTo = 0, n_proc(SP_)-1
            ! translate proc at SP to global
            call MPI_Group_translate_ranks(&
                 i_group(SP_), 1, (/SP_iProcTo/), &
                 i_group(),            iProcTo_I, iError)
            ! count how many points will be sent
            nParticleSend_I(SP_iProcTo) = sum(nParticleAtLine_I, &
                 MASK = &
                 SP_GridDescriptor%DD%Ptr%iDecomposition_II(PE_,:) == iProcTo_I(1))
            ! send number of particles to be transfered
            nRequestS = nRequestS + 1
            call MPI_Isend(nParticleSend_I(SP_iProcTo), 1, MPI_INTEGER, &
                 iProcTo_I(1), iTag, i_comm(), iRequestS_I(nRequestS), iError)
         end do
      end if

      if(is_proc(SP_))then
         !\
         ! Recv # of particles from MH
         !/
         do MH_iProcFrom = 0, n_proc(iMH)-1
            ! translate proc at SC to global
            call MPI_Group_translate_ranks(&
                 i_group(iMH), 1, (/MH_iProcFrom/), &
                 i_group(),            iProcFrom_I, iError)
            ! recv # of particles to be received
            nRequestR = nRequestR + 1
            call MPI_Irecv(nParticleRecv_I(MH_iProcFrom), 1, MPI_INTEGER,&
                 iProcFrom_I(1), iTag, i_comm(), iRequestR_I(nRequestR),iError)
         end do
      end if
      ! finalize transfer
      call MPI_waitall(nRequestR, iRequestR_I, iStatus_II, iError)
      call MPI_waitall(nRequestS, iRequestS_I, iStatus_II, iError)
      !\
      ! send the actual data
      !/
      nRequestS = 0
      nRequestR = 0
      if(is_proc(iMH))then
         allocate(BuffSend_I(nVar * nParticleThisProc))
         iBuff = 1
         do SP_iProcTo = 0, n_proc(SP_)-1
            if(nParticleSend_I(SP_iProcTo) == 0) CYCLE
            ! translate proc at SP to global
            call MPI_Group_translate_ranks(&
                 i_group(SP_), 1, (/SP_iProcTo/), &
                 i_group(),            iProcTo_I, iError)
            ! prepare data
            do iParticle = 1, nParticleThisProc
               if(SP_GridDescriptor%DD%Ptr%iDecomposition_II(PE_,&
                    nint(Particle_II(iFLIndex,iParticle)))/=iProcTo_I(1)) CYCLE
               BuffSend_I(iBuff:iBuff + nVar - 1) = Particle_II(:,iParticle)
               iBuff = iBuff + nVar
            end do
            iBuff = iBuff - nVar*nParticleSend_I(SP_iProcTo)
            ! transfer data
            nRequestS = nRequestS + 1
            call MPI_Isend(BuffSend_I(iBuff),nVar*nParticleSend_I(SP_iProcTo),&
                 MPI_REAL, &
                 iProcTo_I(1), iTag, i_comm(), iRequestS_I(nRequestS), iError)
            iBuff = iBuff + nVar*nParticleSend_I(SP_iProcTo)
         end do
         deallocate(nParticleAtLine_I)
         deallocate(Particle_II)
      end if
      !\
      ! recv the actual data
      !/
      if(is_proc(SP_))then
         allocate(BuffRecv_I(nVar * sum(nParticleRecv_I)))
         iBuff = 1
         do MH_iProcFrom = 0, n_proc(iMH)-1
            if(nParticleRecv_I(MH_iProcFrom)==0)CYCLE
            ! translate proc at SC to global
            call MPI_Group_translate_ranks(&
                 i_group(iMH), 1, (/MH_iProcFrom/), &
                 i_group(),            iProcFrom_I, iError)
            ! recv data
            nRequestR = nRequestR + 1
            call MPI_Irecv(BuffRecv_I(iBuff),&
                 nVar*nParticleRecv_I(MH_iProcFrom), &
                 MPI_REAL,&
                 iProcFrom_I(1), iTag, i_comm(), iRequestR_I(nRequestR),iError)
            iBuff = iBuff + nVar*nParticleRecv_I(MH_iProcFrom)

         end do
      end if
      ! finalize transfer
      call MPI_waitall(nRequestR, iRequestR_I, iStatus_II, iError)
      call MPI_waitall(nRequestS, iRequestS_I, iStatus_II, iError)
      !\
      ! put data
      !/
      if(is_proc(SP_))then
         Convert_DD = transform_matrix(tNow,&
              Grid_C(iMH)%TypeCoord, Grid_C(SP_)%TypeCoord)
         call SP_put_line(NameVar, nVar, sum(nParticleRecv_I),&
              reshape(BuffRecv_I,(/nVar, sum(nParticleRecv_I)/)),&
              iInterfaceType, Convert_DD)
      end if

      ! deallocate arrays for non-blocking communications
      if(is_proc(SP_))&
           deallocate(BuffRecv_I)
      if(is_proc(iMH))&
           deallocate(BuffSend_I)
      deallocate(nParticleRecv_I)
      deallocate(nParticleSend_I)
      deallocate(iRequestS_I)
      deallocate(iRequestR_I)
      deallocate(iStatus_II)
    end subroutine exchange_lines

  end subroutine couple_mh_sp_init
  !==================================================================!

  subroutine transform_from_cartesian(iComp)
    use ModCoordTransform, ONLY: xyz_to_rlonlat
    integer,intent(in)::iComp
    real,pointer,dimension(:,:)::Coord_DI
    character(len=100) :: TypeGeometry
    real:: Coord_D(3), Rho
    real, parameter:: cTol = 1E-8
    integer::nU_I(2), iParticle
    real:: SpToMh_DD(3,3)
    !------------------------------------------
    nU_I = ubound_vector('SP_Xyz_DI')
    call associate_with_global_vector(Coord_DI, 'SP_Xyz_DI')
    ! convert from cartesian coordinates if necessary                        
    TypeGeometry = Grid_C(iComp)%TypeGeometry
    SpToMh_DD=transform_matrix(tNow,&
         Grid_C(SP_)%TypeCoord,&
         Grid_C(iComp)%TypeCoord)
    if( index(TypeGeometry, 'spherical_lnr') > 0 )then
       do iParticle = 1, nU_I(2)
          ! rotate cartesian
          Coord_D(:) = matmul(SpToMh_DD, Coord_DI(:,iParticle))
          call xyz_to_rlonlat(Coord_D, Coord_DI(1:3,iParticle))
          ! convert radius to log(radius)
          if(Coord_DI(1,iParticle) > 0)&
               Coord_DI(1,iParticle) = log(Coord_DI(1,iParticle))
       end do
    elseif(index(TypeGeometry, 'cartesian') > 0 )then
       do iParticle = 1, nU_I(2)
          ! rotate cartesian
          Coord_DI(:,iParticle) = matmul(SpToMh_DD, Coord_DI(:,iParticle))
       end do
    end if
    nullify(Coord_DI)
  end subroutine transform_from_cartesian

  !==================================================================!
  subroutine transform_to_cartesian(iComp)
    use ModCoordTransform, ONLY: rlonlat_to_xyz
    integer,intent(in)::iComp
    real,pointer,dimension(:,:)::Coord_DI
    character(len=100) :: TypeGeometry
    real:: Coord_D(3), Rho
    real, parameter:: cTol = 1E-8
    integer::nU_I(2), iParticle
    real:: MhToSp_DD(3,3)
    !------------------------------------------
    nU_I = ubound_vector('SP_Xyz_DI')
    call associate_with_global_vector(Coord_DI, 'SP_Xyz_DI')
    ! convert from cartesian coordinates if necessary                        
    TypeGeometry = Grid_C(iComp)%TypeGeometry
    MhToSp_DD=transform_matrix(tNow,&
         Grid_C(iComp)%TypeCoord,&
         Grid_C(SP_)%TypeCoord)
    if( index(TypeGeometry, 'spherical_lnr') > 0 )then
       do iParticle = 1, nU_I(2)
          ! convert log(radius) to radius
          if(Coord_DI(1,iParticle) > 0)&
               Coord_DI(1,iParticle) = exp(Coord_DI(1,iParticle))
          call rlonlat_to_xyz(Coord_DI(1:3,iParticle), Coord_D)
          ! rotate cartesian
          Coord_DI(:,iParticle) = Coord_D
          Coord_DI(:,iParticle) = matmul(MhToSp_DD, Coord_DI(:,iParticle))
       end do
    elseif(index(TypeGeometry, 'cartesian') > 0 )then
       do iParticle = 1, nU_I(2)
          ! rotate cartesian
          Coord_DI(:,iParticle) = matmul(MhToSp_DD, Coord_DI(:,iParticle))
       end do
    end if
    nullify(Coord_DI)
  end subroutine transform_to_cartesian

  !==================================================================!
  subroutine transform_to_sp_from(iComp)
    integer,intent(in)::iComp
    real,pointer,dimension(:,:)::SP_LocalXyz_DI
    logical,pointer,dimension(:)::Is_I
    integer::nU_I(2),i
    real,dimension(3,3)::MhToSp_DD
    character(LEN=2)::NameComp
    real::LengthRatio
    call get_comp_info(iComp,Name=NameComp)
    MhToSp_DD=transform_matrix(tNow,&
         Grid_C(iComp)%TypeCoord,&
         Grid_C(SP_)%TypeCoord)
    if(DoTest)write(*,*)'Transform SP coordinates from '//NameComp
    call associate_with_global_mask(Is_I,'SP_IsIn'//NameComp)
    call associate_with_global_vector(SP_LocalXyz_DI,'SP_Xyz_DI')
    nU_I=ubound(SP_LocalXyz_DI)
    if(DoTest)write(*,*)nU_I
    LengthRatio=Grid_C(iComp)%UnitX/Grid_C(SP_)%UnitX
    do i=1,nU_I(2)
       if(.not.Is_I(i))CYCLE
       SP_LocalXyz_DI(:,i)=matmul(MhToSp_DD,&
            point_state_v('SP_Xyz_DI',3,i))*LengthRatio
    end do
    nullify(SP_LocalXyz_DI,Is_I)
  end subroutine transform_to_sp_from

  !==================================================================
  !^CMP IF IH BEGIN
  subroutine couple_ih_sp(DataInputTime)     

    use CON_global_message_pass

    real,intent(in)::DataInputTime
    real,dimension(3)::Xyz_D
    !-------------------------------------------------------------------------
    if(.not.RouterIhSp%IsProc)return

    tNow=DataInputTime
    IhToSp_DD=transform_matrix(tNow,&
         Grid_C(IH_)%TypeCoord, Grid_C(SP_)%TypeCoord)
    ScToIh_DD=transform_matrix(tNow,&                   !^CMP IF SC
         Grid_C(SC_)%TypeCoord, Grid_C(IH_)%TypeCoord)  !^CMP IF SC


    call IH_synchronize_refinement(RouterIhSp%iProc0Source,&
         RouterIhSp%iCommUnion)
    call bcast_global_vector('SP_Xyz_DI',&
         RouterIhSp%iProc0Source,&
         RouterIhSp%iCommUnion,&
         'SP_IsInIH')
    call transform_from_cartesian(IH_)
    call set_router(& 
         GridDescriptorSource=IH_GridDescriptor,&
         GridDescriptorTarget=SP_GridDescriptor,&
         Router=RouterIhSp,&
         NameMappingVector='SP_Xyz_DI',&
         NameMask='SP_IsInIH',&
         interpolate=interpolation_amr_gc)

    call transform_to_cartesian(IH_)
    if(is_proc(SP_))then
!       call SP_put_input_time(DataInputTime)
!       call transform_to_sp_from(IH_)
    end if

    call global_message_pass(RouterIhSp,&
         nVar=11,&
         fill_buffer=IH_get_for_sp_and_transform,&
         apply_buffer=SP_put_from_mh)
    !^CMP IF SC BEGIN
    !This coupler is performed after SC-SP coupling, so that 
    !on SP the updated coordinates are available for those
    !points which passed from SC to IH

!    if(use_comp(SC_))then              
!       if(is_proc0(SP_))then
!          !Check the points which passed from SC to IH:
!          call associate_with_global_vector(Xyz_DI,'SP_Xyz_DI')
!          call associate_with_global_mask(Is_I,'SP_IsInIH')
!
!          do iPoint=1,nPoint
!             if(Is_I(iPoint))CYCLE
!             !Consider a point which before was   not in IH
!             Xyz_D=matmul(ScToIh_DD,Xyz_DI(:,iPoint))*&
!                  (Grid_C(SC_)%UnitX/Grid_C(IH_)%UnitX) 
!             if(is_in_ih(Xyz_D))& ! Now the point is in IH
!                  Xyz_DI(:,iPoint)=Xyz_D
!             !..that is why we convert it to IH coordinates
!          end do
!          nullify(Xyz_DI)
!          nullify(Is_I)
!       end if
!       call bcast_global_vector('SP_Xyz_DI',&
!            RouterIhSp%iProc0Target,&
!            RouterIhSp%iCommUnion)           !^CMP END SC
!    end if
!    call set_mask('SP_IsInIH','SP_Xyz_DI',is_in_ih)
  end subroutine couple_ih_sp
  !==================================================================
  logical function is_in_ih(Xyz_D)
    real,dimension(:),intent(in)::Xyz_D
    real:: R2
    R2 = dot_product(Xyz_D,Xyz_D)
    is_in_ih=R2>=rBoundIh**2.and.&
         all(Xyz_D < xyz_max_d(IH_)).and.all(Xyz_D >= xyz_min_d(IH_))
  end function is_in_ih
  !==================================================================!        
  subroutine IH_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    integer, parameter :: Rho_=1, Ux_=2, Uz_=4, Bx_=5, Bz_=7,&
         BuffX_    =9,BuffZ_=11
    !------------------------------------------------------------
    call IH_get_for_sp(&
         nPartial,iGetStart,Get,w,State_V,nVar)

    State_V(Ux_:Uz_)=&
         transform_velocity(tNow,&
         State_V(Ux_:Uz_),&
         State_V(BuffX_:BuffZ_),&
         Grid_C(IH_)%TypeCoord,Grid_C(SP_)%TypeCoord)

    State_V(Bx_:Bz_)=matmul(IhToSp_DD,State_V(Bx_:Bz_))

    ! transfrom coordinates
    State_V(9:11) =matmul(IhToSp_DD,State_V(9:11))

  end subroutine IH_get_for_sp_and_transform
  !^CMP END IH
  !=========================================================================
  !^CMP IF SC BEGIN
  subroutine couple_sc_sp(DataInputTime)
    use CON_global_message_pass

   real,intent(in)::DataInputTime
   !-------------------------------------------------------
    if(.not.RouterScSp%IsProc)return

    tNow=DataInputTime
    ScToSP_DD=transform_matrix(tNow,&
         Grid_C(SC_)%TypeCoord, Grid_C(SP_)%TypeCoord)

    call SC_synchronize_refinement(RouterScSp%iProc0Source,&
         RouterScSp%iCommUnion)
    call bcast_global_vector('SP_Xyz_DI',&
         RouterScSp%iProc0Source,&
         RouterScSp%iCommUnion,&
         'SP_IsInSC')
    call transform_from_cartesian(SC_)
    call set_router(& 
         GridDescriptorSource=SC_GridDescriptor,&
         GridDescriptorTarget=SP_GridDescriptor,&
         Router=RouterScSp,&
         NameMappingVector='SP_Xyz_DI',&
         NameMask='SP_IsInSC',&
         interpolate=interpolation_amr_gc)
    
    call transform_to_cartesian(SC_)
    if(is_proc(SP_))then
!       call SP_put_input_time(DataInputTime)  
!       call transform_to_sp_from(SC_)
    end if

    call global_message_pass(RouterScSp,&
         nVar=11,&
         fill_buffer=SC_get_for_sp_and_transform,&
         apply_buffer=SP_put_from_mh)
    call set_mask('SP_IsInSC','SP_Xyz_DI',is_in_sc)
  end subroutine couple_sc_sp
  !-------------------------------------------------------------------------
  logical function is_in_sc(Xyz_D)
    real,dimension(:),intent(in)::Xyz_D
    real::R2
    real,save:: RSCMin2 = -1
    real,save:: RSCMax2 = -1
    if(RSCMin2 == -1) RSCMin2 = dot_product(exp(xyz_min_d(SC_)),(/1,0,0/))**2
    if(RSCMax2 == -1) RSCMax2 = dot_product(exp(xyz_max_d(SC_)),(/1,0,0/))**2
    R2=dot_product(Xyz_D,Xyz_D)
    if(use_comp(IH_))then            !^CMP IF IH BEGIN
       is_in_sc=R2>=rBoundSc**2.and.R2<rBoundIh**2.and.&
            R2<RSCMax2.and.R2>=RSCMin2
    else                             !^CMP END IH
       is_in_sc=R2>=rBoundSc**2.and.&
            R2<RSCMax2.and.R2>=RSCMin2
    end if                           !^CMP IF IH
  end function is_in_sc
  !--------------------------------------------------------------------------
  subroutine SC_get_for_sp_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    integer, parameter :: Rho_=1, Ux_=2, Uz_=4, Bx_=5, Bz_=7,&
         BuffX_    =9,BuffZ_=11
    !------------------------------------------------------------
    call SC_get_for_sp(&
         nPartial,iGetStart,Get,w,State_V,nVar)

    State_V(Ux_:Uz_)=&
         transform_velocity(tNow,&
         State_V(Ux_:Uz_),&
         State_V(BuffX_:BuffZ_),&
         Grid_C(SC_)%TypeCoord,Grid_C(SP_)%TypeCoord)

    State_V(Bx_:Bz_)=matmul(ScToSp_DD,State_V(Bx_:Bz_))

    ! transform coordinates
    State_V(9:11) =matmul(ScToSp_DD,State_V(9:11))
  end subroutine SC_get_for_sp_and_transform
  !=========================================================================
  !^CMP END SC
end Module CON_couple_mh_sp
