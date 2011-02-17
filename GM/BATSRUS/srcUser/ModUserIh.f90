!^CFG COPYRIGHT UM
!Revision history:
 
! 02 Oct. 2010 by R. Oran:
! Added user_read_inputs, user_set_ics to allow testing of different
! initial conditions.
! 10 Feb 2010 by R. Oran : added Parker initial conditions, user_face_bcs.
!                          Removed IC's options that are  now in ModUserWaves.f90 
!========================================================================
module ModUser
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_boundary_cells, &
       IMPLEMENTED2 => user_read_inputs,        &
       IMPLEMENTED3 => user_set_ics,            &
       IMPLEMENTED4 => user_face_bcs
 
  include 'user_module.h' !list of public methods

  real, parameter   :: VersionUserModule = 2.0
  character(len=*), parameter :: &
       NameUserModule = 'MHD HELIOSPHERE, Sokolov'

  ! For Parker spiral IC's:
  ! Input from PARAM.in file, solar wind values at ref. point in IO units:
  real    :: BrRefIo, NRefIo, TRefIo, UrSwIo, rRefIo
  logical :: UseParkerIcs = .false.

contains
  !============================================================================
  subroutine user_read_inputs

    use ModMain,      ONLY: lVerbose, TypeCoordSystem
    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_line, read_command, read_var
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    implicit none

    character(len=100) :: NameCommand
    character(len=*), parameter :: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input INNER HELIOSPHERE starts'
    endif
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case('#PARKERICS')
          if (TypeCoordSystem == 'HGR') then
             write(*,*) 'ERROR in PARAM.in: Parker ICs in HGR frame are not implemented'
             call CON_stop('Correct PARAM.in')
          else
             UseParkerIcs = .true.
             call read_var('rRef  [AU]',   rRefIo )
             call read_var('BrRef [nT]',   BrRefIo)
             call read_var('NRef  [n/cc]', NRefIo )
             call read_var('TRef  [K]',    TRefIo )
             call read_var('UrSw  [km/s]', UrSwIo )
          end if
                 
       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input INNER HELIOSPHERE ends'
          endif
          EXIT

       case default
          if(iProc == 0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs
  !========================================================================
  subroutine user_set_ics
    
    use ModMain,        ONLY: globalBLK, unusedBLK   
    use ModSize,        ONLY: nI, nJ, nK, gcn
    use ModGeometry,    ONLY: x_BLK, y_BLK, z_BLK, r_BLK
    use ModAdvance,     ONLY: State_VGB
    use ModPhysics,     ONLY: rBody
    use ModVarIndexes

    implicit none

    integer   :: i,j,k,iBlock
    real      :: x, y, z, r, State_V(nVar)
    character(len=*), parameter :: NameSub = 'user_set_ics'
    !----------------------------------------------------------------------
    iBlock = globalBLK
  
    if (unusedBLK(iBlock)) RETURN
  
    if(UseParkerIcs) then
       do k=1-gcn,nK+gcn ; do j=1-gcn,nJ+gcn ; do i=1-gcn,nI+gcn

          ! make variable names a little shorter
          x = x_BLK(i,j,k,iBlock)
          y = y_BLK(i,j,k,iBlock)
          z = z_BLK(i,j,k,iBlock)
          r = r_BLK(i,j,k,iBlock)

          if (r .le. rBody) CYCLE
       
          call get_parker_sln_cell(x,y,z,r,State_V)   
          State_VGB(:,i,j,k,iBlock) = State_V

       end do ; end do ; end do
    else
       write(*,*) 'You are trying to use user_set_ics in an unspecified manner.'
       call CON_stop('Correct PARAM.in')
    end if

  end subroutine user_set_ics
 !============================================================================
  subroutine get_parker_sln_cell(x, y, z, r, State_V)

    ! Set state to the steady state adiabatic Parker spiral solution.     
    ! R. Oran 

    ! OUTLINE:
    ! This solution is appropriate only for the domain lying outside the sonic point. 
    ! The solar wind is assumed to have reached its terminal speed, which is purely
    ! radial in the inertial frame HGI.
    
    ! The magnetic field is given by the Parker spiral.
    ! It is assumed to be purely radial at  r~10Rs.
   
    ! The density is derived from conservation of mass: rho*Ur*r^2 = Const.
    ! The pressure is given by the adiabatic EOS: p = p0(r0/r)^2gamma

    ! Depending on plasma parameters at rRef, the solution everywhere is found. 
  
    ! INPUTS (to be read from PARAM.in file):
    ! rRef - heliocentric radius of reference point
    ! bRef
    ! Rhoref
    ! pRef
    ! UrSW

    ! OUTPUT:
    ! The state variables for rho, P, U, and B are returned.
      
    use ModMain,           ONLY: TypeCoordSystem
    use ModVarIndexes,     ONLY: Rho_, RhoUx_, RhoUz_,Bx_,Bz_,p_, nVar
    use ModCoordTransform, ONLY: rot_xyz_sph
    use ModNumConst,       ONLY: cTwoPi
    use ModConst,          ONLY: RotationPeriodSun, cAU, cProtonMass
    use ModPhysics,        ONLY: rBody, Gamma0, &
                                 Si2No_V, Io2No_V, No2Si_V, &
                                 UnitX_, UnitN_, UnitB_, UnitU_, UnitRho_, &
                                 UnitP_, UnitTemperature_, UnitT_

    implicit none 

    real,intent(in)  :: x, y, z, r
    real,intent(out) :: State_V(nVar)
    real     :: OmegaSun
    real     :: SinTheta             ! polar angle in spherical coordinates
    real     :: b_D(3), u_D(3),r_D(3)! Cartesian vectors
    real     :: bSph_D(3), uSph_D(3) ! Spherical B and velocity vectors
    real     :: XyzSph_DD(3,3)       ! rotation matrix from spherical to xyz

    ! Solar wind values at reference point, in normalized units.
    real     :: Br0, Rho0, p0, r0, UrSw, N0, T0, MassFlux
    real     :: RotationSpeed

    character(len=*),parameter :: NameSub = 'get_parker_sln_cell'
    !--------------------------------------------------------------------
    ! Normalize and convert input parameters for rho, u, B, p
    r0   = rRefIo * cAu * Si2No_V(UnitX_)
    UrSw = UrSwIo*1000*Si2No_V(UnitU_)
    Br0  = BrRefIo * 1e-9 * Si2No_V(UnitB_)
    N0   = NRefIo * 1e6 * Si2No_V(UnitN_)
    Rho0 = cProtonMass * NRefIo * 1e6 * Si2No_V(UnitRho_)
    T0   = TRefIo *Si2No_V(UnitTemperature_  )
    p0   = Rho0 * T0

    ! Make sure the flow is supersonic at reference point
    if (UrSw .lt. sqrt(Gamma0*p0/Rho0) ) then
       write(*,*) 'UrSw is subsonic at reference point!'
       call CON_stop('Check #PARKERICS command!')
    end if

    MassFlux = Rho0 * UrSw * r0**2 ! conserved!     
    OmegaSun = cTwoPi/(RotationPeriodSun*Si2No_V(UnitT_))
    RotationSpeed = r*No2Si_V(UnitX_)*(cTwoPi/RotationPeriodSun)*&
         Si2No_V(UnitU_)
    !\ 
    ! Calculate parker solution for given position
    !/
    r_D = (/x,y,z/)
     
    ! sine of polar angle
    SinTheta = sqrt(x**2 + y**2)/r

    ! Transformation matrix from Spherical to Cartesian
    ! at this location
    XyzSph_DD = rot_xyz_sph(x,y,z) 
        
    !\
    ! Magnetic field
    !/         
    ! Magnetic field components in spherical coordinates as
    ! given by the Parker solution. Note the azimuthal component:
    ! Should depend on (r - rSourceSurface), since r >> rSourceSurface
    ! it can be neglected.
    bSph_D(1) = -sign(1.0,z)*Br0*(r0/r)**2              ! B_r
    bSph_D(2) = 0.0                                     ! B_theta
    bSph_D(3) = -bSph_D(1)*RotationSpeed*SinTheta/UrSw  ! B_phi

    ! Convert magnetic field vector to Cartesian coordinates
    b_D = matmul(XyzSph_DD,bSph_D)
    ! Init state variable for B
    State_V(Bx_:Bz_) = b_D

    !\
    ! Velocity
    !/
    ! constant Ur in inertial (HGI) frame
    uSph_D(1) = UrSw       ! Ur
    uSph_D(2) = 0.0        ! Utheta
    uSph_D(3) = 0.0        ! Uphi

    ! Convert velocity vector to Cartesian components
    u_D = matmul(XyzSph_DD, uSph_D)

    if (TypeCoordSystem == 'HGC') then
       ! Find the azimuthal component of the velocity
       ! In the co-rotating frame U || B, thus Uphi = Ur(Bphi/Br)
       u_D(1) = u_D(1) + OmegaSun*y
       u_D(2) = u_D(2) - OmegaSun*x
            
    end if
         
    ! Find density from conservation of mass 
    State_V(Rho_) = MassFlux / (sum(r_D*u_D) * r) 

    ! RhoU
    State_V(RhoUx_:RhoUz_) = u_D * State_V(Rho_)

    !\
    ! pressure
    !/
    State_V(p_) = p0 * (r0 / r)**(2*Gamma0)


    ! Make sure the resulting flow is supersonic
    if (UrSw  .lt. sqrt(Gamma0*State_V(p_)/State_V(Rho_)) ) then
       write(*,*) 'Ur is subsomic at: ',x, y, z
       call CON_stop('ERROR in get_parker_sln_cell')
    end if
      
  end subroutine get_parker_sln_cell
  !============================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModMain,        ONLY: x_, y_, z_
    use ModFaceBc,      ONLY: FaceCoords_D, VarsTrueFace_V
    use ModVarIndexes,  ONLY: nVar, Rho_, Ux_, Uz_

    real, intent(out) :: VarsGhostFace_V(nVar)

    real :: x, y, z, r, State_V(nVar), U_D(3)
    !--------------------------------------------------------------------------
    x = FaceCoords_D(x_)
    y = FaceCoords_D(y_)
    z = FaceCoords_D(z_)
    r = sqrt(x**2+y**2+z**2)
    call get_parker_sln_cell(x, y, z, r, State_V)
    VarsGhostFace_V = State_V

    !VarsGhostFace_V(Rho_) =  2.0*State_V(Rho_) - VarsTrueFace_V(Rho_)
  
  end subroutine user_face_bcs
  !========================================================================
  subroutine user_set_boundary_cells(iBLK)

    ! Set the boundary cell information IsBoundaryCell_GI(:,:,:,ExtraBc_) 
    ! for a sphere of radius rBody around the origin.
    ! Allow resolution change.

    use ModGeometry
    use ModBoundaryCells,ONLY:SaveBoundaryCells
    use ModPhysics,ONLY:rBody
    implicit none
    integer, intent(in):: iBLK

    IsBoundaryCell_GI(:,:,:,ExtraBc_) = R_BLK(:,:,:,iBLK)<rBody
    if(SaveBoundaryCells)return
    call stop_mpi('Set SaveBoundaryCells=.true. in PARAM.in file')
  end subroutine user_set_boundary_cells

end module ModUser

