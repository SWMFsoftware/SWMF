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

  ! These values are to be read from the PARAM.in file and used to
  ! set the initial conditions in IH
  character(len=10) :: TypeIcs
  real              :: RhoInitialNo, pInitialNo, rMinShell, rMaxShell
  real              :: uSwFinalSi, BrSourceSi, RhoSourceSi, pSourceSi
  real              :: rSourceNo, UrSourceSi, pShellNo, RhoShellNo
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
          select case(TypeIcs)
          case('uniform','shell')
             call read_var('RhoInitialNo', RhoInitialNo)
             call read_var('pInitialNo', pInitialNo)
             if(TypeIcs=='shell') then
                call read_var('RhoShellNo', RhoShellNo)
                call read_var('pShellNo',   pShellNo)
                call read_var('rMinShell',  rMinShell)
                call read_var('rMaxShell',  rMaxShell)
             end if
          case('parker')
             if (TypeCoordSystem /= 'HGR') then
                write(*,*) 'ERROR in PARAM.in: Cannot use Parker solution in HGR'
                call CON_stop('Correct PARAM.in')
             else
                call read_var('rSourceNo',   rSourceNo )
                call read_var('BrSourceSi',  BrSourceSi)
                call read_var('RhoSourceSi', RhoSourceSi)
                call read_var('pSourceSi',   pSourceSi)
                call read_var('UrSourceSi',  UrSourceSi)
                call read_var('uSwFinalSi', uSwFinalSi )
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
    use ModMain,       ONLY: globalBLK
    use ModVarIndexes, ONLY: rho_, p_
    use ModPhysics,    ONLY: BodyRho_I, BodyP_I
    use ModGeometry,   ONLY: R_BLK
    use ModSize,       ONLY: nI, nJ, nK, gcn

    implicit none

    integer                     :: i, j, k, iBlock
    character(len=*), parameter :: NameSub = 'user_set_ics'
    !----------------------------------------------------------------------
    iBlock = globalBLK

    select case(TypeIcs)
    case('uniform')
       ! This case describes an IC of a plasma at rest with no
       ! density or pressure gradients and no magnetic field.
       ! The initial density and pressure are controled by the user
       ! through a combination of the normalized quantity RhoInitialNo
       ! and the #BODY command.

       State_VGB(:,:,:,:,iBlock)    = 0.0
       State_VGB(rho_,:,:,:,iBlock) = RhoInitialNo*BodyRho_I(1)       
       State_VGB(p_,  :,:,:,iBlock) = pInitialNo*BodyP_I(1)
    case('shell')
       ! this case describes a shell with density and/or pressure
       ! embedded in a uniform background at rest, without a magnetic field.
       ! The background density and pressure are derived from the #BODY command.
       
       State_VGB(:,:,:,:,iBlock)  = 0.0
       do k=-1-gcn,nK+gcn ; do j= 1-gcn,nJ+gcn ; do i=1-gcn,nI+gcn
          
          if ( R_BLK(i,j,k,iBlock) .ge. rMinShell .and. &
               R_BLK(i,j,k,iBlock) .le. rMaxShell ) then
             State_VGB(rho_,i,j,k,iBlock) = RhoShellNo*BodyRho_I(1)       
             State_VGB(p_,  i,j,k,iBlock) = pShellNo*BodyP_I(1)
          else
             State_VGB(rho_,i,j,k,iBlock) = RhoInitialNo*BodyRho_I(1)       
             State_VGB(p_,  i,j,k,iBlock) = pInitialNo*BodyP_I(1)
          end if  
       end do; end do ; end do

    case('parker')
       call set_parker_spiral
          
    case default
       write(*,*) NameSub, ' WARNING: No TypeIcs was selected. Using defaults' 
    end select

  contains
    subroutine set_parker_spiral

      ! Set state to the steady state isothermal Parker spiral solution.
     
      ! IMPORTANT: This initial condition is valid for the co-rotating frame
      ! so it must be used only when working in the HGR coordinate system.
      ! A CON_stop will be issued otherwise).

      ! OUTLINE:
      ! We choose a spherical source surface outside the sonic point, where
      ! the magnetic field is assumed to be purely radial.
      ! Depending on plasma parameters at that location, the Parker solution
      ! everywhere is derived.
      ! Note: the source surface here is not to be confused with the
      ! source surface used for the PFSS model used for the solar corona.
      ! The magnetic and velocity fields are given by the Parker spiral solution
      ! in the co-rotating frame. The density is then derived from conservation
      ! of mass.

      ! INPUTS (to be read from PARAM.in file):
      ! rSource - heliocentric radius of source surface
      ! BrSource
      ! Rho Source
      ! pSource
      ! UrSource
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
      use ModPhysics,        ONLY: gBody, Si2No_V, &
                                   UnitB_, UnitU_,UnitRho_,UnitP_,UnitX_
      use ModConst,          ONLY: RotationPeriodSun, rSun

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

      character(len=*),parameter :: NameSub = 'set_parker_spiral'
      !--------------------------------------------------------------------
      iBlock = globalBLK

      ! normalize 
      uSwFinalNo  = Si2No_V(UnitU_  ) * uSwFinalSi
      BrSourceNo  = Si2No_V(UnitB_  ) * BrSourceSi
      RhoSourceNo = Si2No_V(UnitRho_) * RhoSourceSi
      pSourceNo   = Si2No_V(UnitP_  ) * pSourceSi
      UrSourceNo  = Si2No_v(UnitU_  ) * UrSourceSi
      rSunNo      = Si2No_V(UnitX_  ) * rSun
      MassFlux = RhoSourceNo*UrSourceNo*rSourceNo**2
      
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
    ! for a sphere of radius 1 around the origin.
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

