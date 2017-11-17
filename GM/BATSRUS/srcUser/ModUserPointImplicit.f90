!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModUser

  use BATL_lib, ONLY: &
       test_start, test_stop
  ! This is an example user module demonstrating the use of the
  ! point implicit scheme for the user source terms.
  ! Please also see the documentation, and the ModPointImplicit.f90 file.

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_calc_sources,               &
       IMPLEMENTED3 => user_init_point_implicit

  include 'user_module.h' ! list of public methods

  !\
  ! Here you must define a user routine Version number and a
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'EXAMPLE FOR POINT IMPLICIT SOURCE, Toth'

  ! Local variables with default values
  real :: rFriction = 1.0, TauFriction = 1.0

contains
  !============================================================================
  subroutine user_read_inputs

    use ModReadParam
    character (len=100) :: NameCommand

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#FRICTION")
          call read_var('rFriction',rFriction)
          call read_var('TauFriction',TauFriction)
       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi('ERROR in ModUserPointImplicit: unknown command='//&
               NameCommand)
       end select
    end do

    call test_stop(NameSub, DoTest)
  end subroutine user_read_inputs
  !============================================================================

  subroutine user_calc_sources(iBlock)

    ! Evaluate the explicit or implicit or both source terms.
    ! If there is no explicit source term, the subroutine user_expl_source
    ! and the corresponding calls can be removed.

    use ModPointImplicit, ONLY:  UsePointImplicit, UsePointImplicit_B, &
         IsPointImplSource, iVarPointImpl_I, IsPointImplMatrixSet, DsDu_VVC
    use ModMain,    ONLY: nI, nJ, nK
    use ModAdvance, ONLY: State_VGB, Source_VC, &
         Rho_, RhoUx_, RhoUy_, RhoUz_, Energy_
    use ModGeometry, ONLY: r_BLK, rMin_BLK

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real    :: Coef
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_calc_sources'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    ! Only blocks within radius rFriction need to be point implicit
    ! This array is also used for load balancing
    ! if DoBalancePointImplicit is set to true in user_init_point_implicit.
    UsePointImplicit_B(iBlock) = rMin_BLK(iBlock) < rFriction

    ! Check if point implicit scheme is on (part implicit may switch it off)
    ! Also check if this particular block is point implicit or not
    if(.not.UsePointImplicit)then
       ! Add all source terms if we do not use the point implicit scheme
       call user_expl_source
       if(UsePointImplicit_B(iBlock)) call user_impl_source
    elseif(IsPointImplSource)then
       ! Add implicit sources only
       if(UsePointImplicit_B(iBlock)) call user_impl_source
    else
       ! Add explicit sources only
       call user_expl_source
    end if

    call test_stop(NameSub, DoTest, iBlock)
  contains
    !==========================================================================
    subroutine user_expl_source

      ! Add explicit source here
      ! The energy source is only needed in the explicit source term.

      ! In this example a simple friction term is added to the momentum and
      ! energy equations.
      !------------------------------------------------------------------------
      Coef = 1.0/TauFriction

      ! Here come the explicit source terms
      do k=1,nK; do j=1,nJ; do i=1,nI
         if(r_BLK(i,j,k,iBlock) > rFriction ) CYCLE

         ! Kinetic energy changes by F.v = (1/TauFriction)*rhoU^2/rho
         Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
              - Coef * sum(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)**2) &
              / State_VGB(Rho_,i,j,k,iBlock)

      end do; end do; end do

    end subroutine user_expl_source
    !==========================================================================
    subroutine user_impl_source

      ! This is a test and example for using point implicit source terms
      ! Apply friction relative to some medium at rest
      ! The friction force is proportional to the velocity and the density.

      ! Add implicit source here
      ! In this example a simple friction term is added to the momentum
      ! equation. Note that the energy is a dependent variable in the
      ! point implicit scheme, so there is no energy source here.
      ! The pressure source is zero.
      !------------------------------------------------------------------------
      Coef = 1.0/TauFriction

      do k=1,nK; do j=1,nJ; do i=1,nI
         if(r_BLK(i,j,k,iBlock) > rFriction ) CYCLE

         ! Friction force F = (1/TauFriction)*RhoU
         Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
              - Coef * State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)

      end do; end do; end do

      if(IsPointImplMatrixSet)then
         ! Set the non-zero dS/dU matrix elements here
         do k=1,nK; do j=1,nJ; do i=1,nI
            if( r_BLK(i,j,k,iBlock) > rFriction ) CYCLE
            DsDu_VVC(RhoUx_,RhoUx_,i,j,k) = - Coef
            DsDu_VVC(RhoUy_,RhoUy_,i,j,k) = - Coef
            DsDu_VVC(RhoUz_,RhoUz_,i,j,k) = - Coef
         end do; end do; end do
      end if

    end subroutine user_impl_source
    !==========================================================================
  end subroutine user_calc_sources
  !============================================================================

  subroutine user_init_point_implicit

    use ModVarIndexes, ONLY: RhoUx_, RhoUy_, RhoUz_
    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet, &
         DoBalancePointImplicit, IsDynamicPointImplicit
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_init_point_implicit'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Allocate and set iVarPointImpl_I
    ! In this example there are 3 implicit variables
    allocate(iVarPointImpl_I(3))

    ! In this example the implicit variables are the 3 momenta
    iVarPointImpl_I = (/RhoUx_, RhoUy_, RhoUz_/)

    ! Note that energy is not an independent variable for the
    ! point implicit scheme. The pressure is an independent variable,
    ! and in this example there is no implicit pressure source term.

    ! Tell the point implicit scheme if dS/dU will be set analytically
    ! If this is set to true the DsDu_VVC matrix has to be set below.
    IsPointImplMatrixSet = .true.

    ! Tell the load balancing scheme if point implicit blocks should be
    ! load balanced. If set to true, UsePointImplicit_B should be set in
    ! user_calc_sources. Default is false.
    DoBalancePointImplicit = .true.

    ! Tell the load balancing scheme if the assignment of UsePointImplicit_B
    ! keeps changing so load balancing is needed every time step.
    ! Default is false.
    IsDynamicPointImplicit = .true.

    call test_stop(NameSub, DoTest)
  end subroutine user_init_point_implicit
  !============================================================================

end module ModUser
!==============================================================================
