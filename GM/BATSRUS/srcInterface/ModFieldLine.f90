!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModFieldLine
  ! the module contains subroutine for extracting magnetic field lines 
  ! for passing to other codes;
  ! field line intergration is performed with use of BATL library including
  ! - continuous AMR interpolation
  ! - particle methods
  !============================================================================
  
  use BATL_lib, ONLY: &
       iProc, &
       MaxDim, nDim, &
       interpolate_grid_amr_gc, check_interpolate_amr_gc, &
       Particle_I, CellSize_DB, &
       message_pass_particles, remove_undefined
  use ModAdvance, ONLY: State_VGB
  use ModVarIndexes, ONLY: Bx_, By_, Bz_
  use ModB0, ONLY: get_b0

  implicit none
  save


  ! sorts of particles used to generate a magnetic field line
  integer, parameter:: &
       SortEnd_     = 1, &
       SortReg_ = 2

  ! variable in the state vector of a particle
  integer, parameter:: &
       ! coordinates of a particle
       x_    = 1, y_    = 2, z_    = 3, & 
       ! auxilary position, e.g. middle step in Runge-Kutta 2 method
       AuxX_ = 3, AuxY_ = 4, AuxZ_ = 5, & 
       ! auxilary field, e.g. stepsize
       Aux_   = 6, & 
       ! field line this particle lays on
       fl_    = 7, &
       ! index of the particle along this field line
       Index_ = 8


  contains
    
    subroutine field_line_extract(nFieldLine, XyzInit_DI)
      ! extract nFieldLine magnetic field lines starting at XyzInit_DI;
      ! the whole field lines are extracted, i.e. they are traced forward
      ! and backward up until it reaches boundaries of the domain
      !------------------------------------------------------------------------
      integer, intent(in):: nFieldLine
      real,    intent(in):: XyzInit_DI(MaxDim, nFieldLine)

      ! indicators of the direction of integration
      integer, parameter:: iDirForward = 1, iDirBackward = -1

      integer :: iFieldLine ! loop variable
      integer :: iBlock
      integer :: iParticle ! loop variable
      logical :: IsCompleted_I(nFieldLine)

      ! current position
      real:: XyzCurrent_D(MaxDim)
      real:: Dir_D(MaxDim)

      ! parameters of end particles
      real,    pointer:: StateEnd_VI(:,:)
      integer, pointer:: iBlockEnd_I(:)
      ! parameters of regular particles
      real,    pointer:: StateReg_VI(:,:)
      integer, pointer:: iBlockReg_I(:)
      !------------------------------------------------------------------------
      ! set a pointers to parameters of end particles
      nullify(StateEnd_VI); nullify(iBlockEnd_I)
      StateEnd_VI => Particle_I(SortEnd_)%State_VI
      iBlockEnd_I => Particle_I(SortEnd_)%iBlock_I

      ! set a pointers to parameters of regular particles
      nullify(StateReg_VI); nullify(iBlockReg_I)
      StateReg_VI => Particle_I(SortReg_)%State_VI
      iBlockReg_I => Particle_I(SortReg_)%iBlock_I

      ! initialize field lines
      Particle_I(SortEnd_)%nParticle = 0
      do iFieldLine = 1, nFieldLine
         call init_field_line(XyzInit_DI(:, iFieldLine), iFieldLine)
      end do
      
      !\
      ! Trace field lines forward, i.e. in the direction of the magnetic field
      !/
      ! trace the field lines
      TRACEFORWARD: do

         ! copy last known coordinates to Auxilary 
         StateEnd_VI(AuxX_:AuxZ_,1:Particle_I(SortEnd_)%nParticle) = &
              StateEnd_VI(x_:z_, 1:Particle_I(SortEnd_)%nParticle)

         call copy_end_to_regular

         !\
         ! predictor step
         !/
         do iParticle = 1, Particle_I(SortEnd_)%nParticle
            ! get the direction of the magnetic field at original location
            call get_b_dir(StateEnd_VI(x_:z_, iParticle), Dir_D)
            ! find the step size 
            StateEnd_VI(Aux_, iParticle) = &
                 0.1 * sum(CellSize_DB(1:nDim,iBlockEnd_I(iParticle))**2 )**0.5
            ! get middle location
            StateEnd_VI(x_:z_, iParticle) = StateEnd_VI(x_:z_, iParticle) + &
                 0.5 * StateEnd_VI(Aux_, iParticle) * Dir_D
         end do

         !\
         ! Message pass: some particles may have moved to different procs
         !/
         call message_pass_particles
         ! also remove particle that went outside of the domain
         call remove_undefined(SortEnd_)
         
         !\
         ! Corrector step
         !/
         do iParticle = 1, Particle_I(SortEnd_)%nParticle
            ! get the direction of the magnetic field in the middle
            call get_b_dir(StateEnd_VI(x_:z_, iParticle), Dir_D)
            ! get final location
            StateEnd_VI(x_:z_,iParticle) = StateEnd_VI(AuxX_:AuxZ_,iParticle)+&
                 StateEnd_VI(Aux_, iParticle) * Dir_D
         end do

         !\
         ! Message pass: some particles may have moved to different procs
         !/
         call message_pass_particles
         ! also remove particle that went outside of the domain
         call remove_undefined(SortEnd_)

         ! increase particle index
         StateEnd_VI(fl_,1:Particle_I(SortEnd_)%nParticle) = &
              StateEnd_VI(fl_,1:Particle_I(SortEnd_)%nParticle) + 1
              
         ! check if all field lines have been completed
         if(Particle_I(SortEnd_)%nParticle == 0) EXIT TRACEFORWARD
      end do TRACEFORWARD
      !\
      ! Trace field line backward, i.e. in the direction of the magnetic field
      !/
      !      TRACEBACKWARD: do
      !      end do TRACEBACKWARD
      
    contains

      subroutine init_field_line(XyzStart_D, iFieldLine)
        real,    intent(in) :: XyzStart_D(MaxDim)
        integer, intent(in) :: iFieldLine

        real   :: Coord_D(MaxDim) ! generalized coordinates
        integer:: iProcOut, iBlockOut
        ! variables to call check_interpolate_amr_gc
        logical:: IsPossible, IsBoundary
        !----------------------------------------------------------------------
        ! find block and processor suitable for interpolation
        call check_interpolate_amr_gc(XyzStart_D, &
             1, IsPossible, & ! input block ID & output IsPossible don't matter
             iProcOut, iBlockOut, IsBoundary)
        !call find_block_to_interpolate_gc(Coord_D, iProcOut, iBlockOut)

        ! check whether point is outside of the domain
        if(IsBoundary) RETURN

        !\
        ! Assign particle to an appropriate processor
        !/
        if(iProc /= iProcOut) RETURN

        Particle_I(SortEnd_)%nParticle = &
             Particle_I(SortEnd_)%nParticle + 1         

        StateEnd_VI(x_:z_, Particle_I(SortEnd_)%nParticle) = XyzStart_D
        StateEnd_VI(fl_,   Particle_I(SortEnd_)%nParticle) = iFieldLine
        StateEnd_VI(Index_,Particle_I(SortEnd_)%nParticle) = 1
        iBlockEnd_I(Particle_I(SortEnd_)%nParticle) = iBlockOut
      end subroutine init_field_line
      !========================================================================
      subroutine copy_end_to_regular
        ! copies known end particles to regular particles
        !----------------------------------------------------------------------
        StateReg_VI(x_:z_,&
             Particle_I(SortReg_)%nParticle+1:&
             Particle_I(SortReg_)%nParticle+Particle_I(SortEnd_)%nParticle) =&
             StateEnd_VI(x_:z_, 1:Particle_I(SortEnd_)%nParticle)
        iBlockReg_I(&
             Particle_I(SortReg_)%nParticle+1:&
             Particle_I(SortReg_)%nParticle+Particle_I(SortEnd_)%nParticle) =&
             iBlockEnd_I(1:Particle_I(SortEnd_)%nParticle)
        Particle_I(SortReg_)%nParticle = &
             Particle_I(SortReg_)%nParticle + &
             Particle_I(SortEnd_)%nParticle
      end subroutine copy_end_to_regular
      !========================================================================
      subroutine get_b_dir(Xyz_D, Dir_D)
        ! returns the direction of magnetic field 
        ! as well as the block used for interpolation
        real,    intent(in) :: Xyz_D(MaxDim)
        real,    intent(out):: Dir_D(MaxDim)

        ! magnetic field
        real   :: B_D(MaxDim) = 0.0
        ! inteprolation data: number of cells, cell indices, weights
        integer:: nCell, iCell_II(0:nDim, 2**nDim)
        real   :: Weight_I(2**nDim)
        integer:: iCell ! loop variable
        integer:: i_D(MaxDim) = 1
        integer:: iBlock
        !----------------------------------------------------------------------
        ! get potential part of the magnetic field at the current location
        call get_b0(Xyz_D, B_D)
        ! get the remaining part of the maagnetic field
        call interpolate_grid_amr_gc(Xyz_D, nCell, iCell_II, Weight_I)
        ! interpolate magnetic field value
        do iCell = 1, nCell
           i_D(1:nDim) = iCell_II(1:nDim, iCell)
           iBlock      = iCell_II(0, iCell)
           B_D = B_D + &
                State_VGB(Bx_:Bz_,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell)
        end do
        ! normalize vector to unity
        Dir_D(1:nDim) = B_D(1:nDim) / sum(B_D(1:nDim)**2)**0.5
      end subroutine get_b_dir

    end subroutine field_line_extract
  
end module ModFieldLine
