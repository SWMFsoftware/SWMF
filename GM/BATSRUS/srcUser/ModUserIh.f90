!^CFG COPYRIGHT UM
!Revision history:
 
! 02 Oct. 2010 by R. Oran:
! Added user_read_inputs, user_set_ics to allow testing of different
! initial conditions. 
!========================================================================
module ModUser
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_boundary_cells, &
       IMPLEMENTED2 => user_read_inputs, &
       IMPLEMENTED3 => user_set_ics
 
  include 'user_module.h' !list of public methods


  real, parameter   :: VersionUserModule = 1.0
  character(len=*), parameter :: &
       NameUserModule = 'HELIOSPHERE, Sokolov'

  ! These variables are assigned values read from the PARAM.in file
  ! Used to set the initial conditions of test problems

  character(len=20) :: TypeIcs
  ! Uniform background  
  real              :: RhoInitialNo, pInitialNo

  ! For ICs containing uniform flow ('UniformU' , 'SphereAdvect')
  real              :: uInitialNo, FlowAngle ! Flow limited to XY plane

  ! Constant density sphere in uniform flow, initially at origin
  real              :: rSphere, RhoSphereNo         

 ! For Parker spiral IC's:
  ! Input from PARAM.in file, solar wind values at 1AU in IO units:
  real              :: BrOneAuIo, RhoOneAuIo, pOneAuIo, UrOneAuIo

  ! Solar wind values scaled to r = rSource, in Normalized units
  real              :: BrSourceNo, RhoSourceNo, pSourceNo
  real              :: rSourceNo, UrSourceNo

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
       case('#USERICS')
          call read_var('TypeIcs', TypeIcs)
          call read_var('RhoInitialNo', RhoInitialNo)
          call read_var('pInitialNo', pInitialNo)
          select case(TypeIcs)
          case('UniformU','SphereAdvect')
             call read_var('uInitialNo',uInitialNo)
             call read_var('FlowAngle', FlowAngle)  
             if (TypeIcs=='SphereAdvect') then
                call read_var('rSphere', rSphere)
                call read_var('RhoSphereNo',RhoSphereNo)
             else
                rSphere = 0.0
                RhoSphereNo = RhoInitialNo
             end if

          case('parker')
             if (TypeCoordSystem /= 'HGR') then
                write(*,*) 'ERROR in PARAM.in: Cannot use Parker solution in HGR'
                call CON_stop('Correct PARAM.in')
             else
                call read_var('rSourceNo',   rSourceNo )
                call read_var('BrOneAuSi',  BrOneAuIo)
                call read_var('RhoOneAuSi', RhoOneAuIo)
                call read_var('pOneAuSi',   pOneAuIo)
                call read_var('UrOneAuSi',  UrOneAuIo)
             end if
          end select
         
       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input INNERHELIOSPHERE ends'
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
    
    use ModAdvance,    ONLY: State_VGB
    use ModMain,       ONLY: globalBLK, unusedBLK, TypeCoordSystem
    use ModVarIndexes
    use ModPhysics,    ONLY: BodyRho_I, BodyP_I, No2Io_V, UnitRho_, UnitU_,UnitP_
    use ModGeometry,   ONLY: R_BLK, x_BLK,y_BLK
    use ModSize,       ONLY: nI, nJ, nK, gcn
    use ModNumConst,   ONLY: cDegToRad
    use ModConst,      ONLY: RotationPeriodSun
    implicit none

    integer                     :: i, j, k, iBlock
    real                        :: FlowAngleRad, Ux0, Uy0
    real                        :: Omega = RotationPeriodSun
    character(len=*), parameter :: NameSub = 'user_set_ics'
    !----------------------------------------------------------------------
    iBlock = globalBLK

     if (unusedBLK(iBlock)) RETURN
   
     select case(TypeIcs)
    case('UniformU','SphereAdvect')
       ! These cases describe an IC with uniform 1D flow of plasma with no
       ! density or pressure gradients and no magnetic field.
       ! The 'SphereAdvect' case also includes an embedded higher/lower
       ! density sphere, initially at the origin.  
       ! Calculate Ux and Uy in inertial frame.
       ! Flow angle is measured from the x axis
       FlowAngleRad = cDegToRad*FlowAngle
       Ux0 = uInitialNo*cos(FlowAngleRad)
       Uy0 = uInitialNo*sin(FlowAngleRad)
       !\
       ! Start filling in cells (including ghost cells)
       !/ 
       if (TypeIcs =='SphereAdvect') then
          do k= 1-gcn,nK+gcn ; do j= 1-gcn,nJ+gcn ; do i=1-gcn,nI+gcn          
             if (R_BLK(i,j,k,iBlock) .le. rSphere)then
                ! inside the sphere
                State_VGB(rho_,i,j,k,iBlock) = RhoSphereNo       
             else
                ! in background flow
                State_VGB(rho_,i,j,k,iBlock) = RhoInitialNo
             end if
          end do; end do ; end do
       else 
          State_VGB(rho_,:,:,:,iBlock) = RhoInitialNo
       end if
       ! velocity
       State_VGB(RhoUx_,:,:,:,iBlock) = Ux0*&
            State_VGB(rho_,:,:,:,iBlock)
       State_VGB(RhoUy_,:,:,:,iBlock) = Uy0*&
            State_VGB(rho_,:,:,:,iBlock)
       if(TypeCoordSystem=='HGR') then
          ! transform velocity to rotating frame
          ! rotation axis is parallel to Z axis.
          State_VGB(RhoUx_,:,:,:,iBlock) = State_VGB(RhoUx_,:,:,:,iBlock)&
               + Omega*y_BLK(:,:,:,iBlock)*State_VGB(Rho_,:,:,:,iBlock)
          State_VGB(RhoUy_,:,:,:,iBlock) = State_VGB(RhoUy_,:,:,:,iBlock)&
               - Omega*x_BLK(:,:,:,iBlock)*State_VGB(Rho_,:,:,:,iBlock)
       end if

       State_VGB(RhoUz_,:,:,:,iBlock) = 0.0
       State_VGB(Bx_:Bz_,:,:,:,iBlock) = 0.0
       State_VGB(p_,  :,:,:,iBlock) = pInitialNo*BodyP_I(1) 
     
    case('parker')
       call set_parker_spiral
          
    case default
       write(*,*) 'You are using user_set_ics with invalid TypeIcs =',TypeIcs
       call CON_stop('Correct PARAM.in')
    end select

   
  contains
    subroutine set_parker_spiral

      ! Set state to the steady state isothermal Parker spiral solution.
     
      ! IMPORTANT: This initial condition is valid for the co-rotating frame
      ! so it must be transformed when working in the HGR coordinate system.
      ! A CON_stop will be issued otherwise).

      ! OUTLINE:
      ! We choose a spherical source surface outside the sonic point, where
      ! the magnetic field is assumed to be purely radial.
      ! Depending on plasma parameters at that location, the Parker solution
      ! everywhere is derived.
      ! Note: the source surface here is not to be confused with the
      ! source surface defined for the PFSS model used for the solar corona.

      ! The magnetic and velocity fields are given by the Parker spiral solution
      ! in the co-rotating frame. The density is then derived from conservation
      ! of mass.

      ! INPUTS (to be read from PARAM.in file):
      ! rSource - heliocentric radius of source surface
      ! BrOneAu
      ! RhoOneAu
      ! pOneAu
      ! UrOneAu
      ! uSwFinal - final speed of solar wind, at infinity
     
      ! OUTPUT:
      ! The state variables for rho, P, U, and B are initialized with
      ! the appropriate value in each cell, including ghost cells.
      
      use ModMain,           ONLY: globalBLK
      use ModAdvance,        ONLY: State_VGB
      use ModVarIndexes
      use ModSize,           ONLY: nI, nJ, nK, gcn
      use ModGeometry,       ONLY: x_BLK, y_BLK, z_BLK, r_BLK, true_cell
      use ModCoordTransform, ONLY: rot_xyz_sph
      use ModPhysics,        ONLY: gBody, Si2No_V, Io2No_V, &
                                   UnitB_, UnitU_,UnitRho_,UnitP_,UnitX_
      use ModConst,          ONLY: RotationPeriodSun, rSun, cAU

      implicit none

      integer  :: i, j, k, iBlock
      real     :: x, y, z, r
      real     :: SinTheta             ! polar angle in spherical coordinates
      real     :: b_D(3), v_D(3)       ! Cartesian B and velocity vectors.
      real     :: bSph_D(3), vSph_D(3) ! Spherical B and velocity vectors
      real     :: XyzSph_DD(3,3)       ! rotation matrix from spherical to xyz

      ! Variables for Parker solution
      real     :: uSwFinalNo, BrSourceNo, RhoSourceNo, pSourceNo, UrSourceNo
      real     :: UsoundNo, rTransonicNo, rSunNo
      
      real     :: V0, V1, MassFlux, Rho, Const, dv
      integer  :: Iteration
      real     :: cAuNo

      character(len=*),parameter :: NameSub = 'set_parker_spiral'
      !--------------------------------------------------------------------
      iBlock = globalBLK

      cAuNo = cAU*Si2No_V(UnitX_)

      ! scale and normalize 
      UrSourceNo  = Io2No_V(UnitU_  ) * UrOneAuIo
      BrSourceNo  = Io2No_V(UnitB_  ) * BrOneAuIo*(cAuNo/rSourceNo)**2
      RhoSourceNo = Io2No_V(UnitRho_) * RhoOneAuIo*(cAuNo/rSourceNo)**2
      pSourceNo   = Io2No_V(UnitP_  ) * pOneAuIo
      ! get rSun in current normalized unit, needed for components that are not
      ! normalized by rSun.
      rSunNo      = Si2No_V(UnitX_  ) * rSun

      MassFlux = RhoSourceNo*UrSourceNo*rSourceNo**2 ! conserved!
      
      ! Speed of sound (isothermal, gamma=1)
      UsoundNo = sqrt(pSourceNo/RhoSourceNo)
    
      !Transonic point, where U=Usound
      rTransonicNo = 0.25*(-gBody**2)/(rSunNo*USoundNo**2)
      
      ! check that the source surface is outside the sonic point
      if (rTransonicNo <  rSourceNo) then
         write(*,*) 'Error in ',NameSub
         call CON_stop('rSource is inside the sonic point. Check PARAM.in')
      end if

      do k=1-gcn,nK+gcn ; do j=1-gcn,nJ+gcn ; do i=1-gcn,nI+gcn

         ! make variable names a little shorter
         x = x_BLK(i,j,k,iBlock)
         y = y_BLK(i,j,k,iBlock)
         z = z_BLK(i,j,k,iBlock)
         r = r_BLK(i,j,k,iBlock)

         ! calculate sine of polar angle
         SinTheta = sqrt(x**2 + y**2)/r
         ! calculate the transformation matrix from Spherical to Cartesian
         ! for each cell
         XyzSph_DD = rot_xyz_sph(x,y,z) 
        
         ! Magnetic field components in spherical coordinates as
         ! given by the Parker solution.
         bSph_D(1) = BrSourceNo*(rSourceNo/r)**2          ! B_r
         bSph_D(2) = 0.0                                  ! B_theta
         bSph_D(3) = -bSph_D(1)*RotationPeriodSun* &      ! B_phi          
              SinTheta/uSwFinalNo 
        
         ! Convert magnetic field vector to Cartesian coordinates
         b_D = matmul(XyzSph_DD,bSph_D)
         ! Init state variable for B
         State_VGB(Bx_:Bz_,i,j,k,iBlock) = b_D


         ! Newton's method to find Ur(r)
         ! Modified version of the code by van det Holst, found in ModUserScHeat.f90
         ! Note: V0, V1 below represent Ur/Usound 
         V0 = 1.0
         Const = 4.0*log(r/rTransonicNo)+3.0 ! part of Parker solution
                                            
         Iteration = 0            
         do
            Iteration = Iteration + 1
            V1 = V0 - (V0**2 -2*log(V0)-Const)*V0/(2.0*(V0**2-1))
            dV = abs(V1 - V0)
            if(dV < 1.0e-7 ) then
               ! The solution was found
               vSph_D(1) = V1*UsoundNo
               EXIT
            elseif(Iteration < 1000) then
               V0 = V1
               CYCLE
            else
               call CON_stop('Parker velocity solution exceeds 1000 iterations')
            end if
         end do

         ! Find the azimuthal component of the velocity
         ! In the co-rotating frame U || B, thus Uphi = Ur(Bphi/Br)
         vSph_D(3) = vSph_D(1)*bSph_D(3)/bSph_D(1)

         ! Convert velocity vector to Cartesian components
         v_D = matmul(XyzSph_DD, vSph_D)
        
         ! Calculate density from mass conservation
         ! rho*Ur*r**2 = constant
         Rho = MassFlux/(vSph_D(1)*r**2)
         State_VGB(Rho_,i,j,k,iBlock) = Rho
         ! Set momentum state variables
         State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = Rho*v_D

         ! pressure
         State_VGB(p_,i,j,k,iBlock) = Rho*pSourceNo/RhoSourceNo
      end do ; end do ; end do

    end subroutine set_parker_spiral

  end subroutine user_set_ics
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

