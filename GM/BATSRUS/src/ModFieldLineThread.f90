module ModFieldLineThread
  use ModMagnetogram,ONLY:interpolate_field
  !Input argument: R,Phi_,Theta_ coordinates
  !Theta is the colatitude
  !Output arguments: B_R, B_Phi, B_Theta,
  !Positive B_Theta is directed from the point Theta=0 (northern pole) 
  !toward Theta = Pi (southern pole) 

  use ModMagnetogram,ONLY:correct_angles
  !Returns the corrected values of the polar angles 
  !(if Theta<0, or Theta>cPi or Phi<0, or Phi>2Pi)

  use ModMagnetogram, ONLY: nR, nPhi, nTheta, &
       Dr, Dphi, dSinTheta
  use ModMain, ONLY: UseFieldLineThreads
  implicit none
  save
  logical, public, allocatable:: DoThreads_B(:)
  ! Named indexes for local use only
  integer, parameter:: r_=1, Phi_=2, Theta_=3

  integer::iIteration
contains
  subroutine read_threads
    use BATL_size, ONLY: MaxBlock
    use ModReadParam, ONLY: read_var
    !-------------
    call read_var('UseFieldLineTraeds', UseFieldLineThreads)
    if(UseFieldLineThreads)then
       if(.not.allocated(DoThreads_B))&
            allocate(DoThreads_B(MaxBlock))
       DoThreads_B = .false.
    else
       if(allocated(DoThreads_B))deallocate(DoThreads_B)
    end if
  end subroutine read_threads
  !=========================
  subroutine set_threads
  end subroutine set_threads
  ! This fucnction calculates the value of 
  ! F(i)= B(i)/|B|/(1,r*sin(colatitude),r)
  ! It is used in order to construct the line in 
  ! spherical coordinates, r_,Phi_,Theta_, where 
  ! Theta_ is the colatidude.
  !==========================================================================
  function f_d(RIn_D)
    ! This fucnction calculates the value of 
    ! F(i)= B(i)/|B|/(1,r*sin(colatitude),r)
    real :: f_d(3)
    real, intent(in) :: RIn_D(3)

    real, parameter :: cTol= 1.0e-10

    !Get the vector (B_r,B_phi,B_theta)
    call interpolate_field(RIn_D,f_d)

    !Divide by the metric coefficients, to obtain
    !the vector ||B|| d (r,phi,theta)/dS along the field line

    f_d=f_d/(/1.0,RIn_D(R_)*max(sin(RIn_D(Theta_)),cTol),&
         & RIn_D(R_)/)

    !Divide by some scale, to limit the displacement within the
    ! integration 
    !step
    f_d=f_d/sqrt(sum(f_d**2))
  end function f_d
  !==========================================================================
  subroutine advance_line_point(RInOut_D, Dir)
    real,intent(inout):: RInOut_D(3)
    real,intent(in)   :: Dir
    real:: dS
    !-----------------
    dS=0.25*min(dR,dPhi,dSinTheta,1.0)*2.0**(iIteration/ (20&
         &*max(nR,nPhi,nTheta)))
    !To avoid the line bouncing near null points
    RInOut_D=RInOut_D+Dir*dS*f_d( RInOut_D+Dir*dS*0.5&
         &*f_d(RInOut_D))
    call correct_angles(RInOut_D)
  end subroutine advance_line_point
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
