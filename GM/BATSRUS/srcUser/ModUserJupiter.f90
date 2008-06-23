!========================================================================
Module ModUser
  use ModVarIndexes, ONLY: rho_, Ux_, Uy_, Uz_,p_,Bx_, By_, Bz_, Energy_, &
       rhoUx_,rhoUy_,rhoUz_
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_calc_sources

  include 'user_module.h' !list of public methods
 
  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Saturn Magnetosphere, Hansen, Nov, 2006'

  real :: MassLoadingRate
  logical :: UseMassLoading = .false., AccelerateMassLoading

  !\
  ! The following are needed in user_sources::
  !/
  real, dimension(1:nI,1:nJ,1:nK):: &
       Srho,SrhoUx,SrhoUy,SrhoUz,SBx,SBy,SBz,Sp,SE


contains

  !=====================================================================
  subroutine user_calc_sources
    use ModAdvance, ONLY: Source_VC
    use ModMain, ONLY: iTest, jTest, kTest, ProcTest, BlkTest, &
        GLOBALBLK
    use ModProcMH,   ONLY: iProc
    implicit none

    logical :: oktest,oktest_me
    !------------------------------------------------------------------------  
    if(iProc==PROCtest .and. globalBLK==BLKtest)then
       call set_oktest('user_calc_sources',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    ! Set the source arrays for this block to zero
    Srho   = 0.0
    SrhoUx = 0.0
    SrhoUy = 0.0
    SrhoUz = 0.0
    SBx    = 0.0
    SBy    = 0.0
    SBz    = 0.0
    SP     = 0.0
    SE     = 0.0

    call user_sources

    Source_VC(rho_   ,:,:,:) = Srho   + Source_VC(rho_   ,:,:,:)
    Source_VC(rhoUx_ ,:,:,:) = SrhoUx + Source_VC(rhoUx_ ,:,:,:)
    Source_VC(rhoUy_ ,:,:,:) = SrhoUy + Source_VC(rhoUy_ ,:,:,:)
    Source_VC(rhoUz_ ,:,:,:) = SrhoUz + Source_VC(rhoUz_ ,:,:,:)
    Source_VC(Bx_    ,:,:,:) = SBx    + Source_VC(Bx_    ,:,:,:)
    Source_VC(By_    ,:,:,:) = SBy    + Source_VC(By_    ,:,:,:)
    Source_VC(Bz_    ,:,:,:) = SBz    + Source_VC(Bz_    ,:,:,:)
    Source_VC(P_     ,:,:,:) = SP     + Source_VC(P_     ,:,:,:)
    Source_VC(Energy_,:,:,:) = SE     + Source_VC(Energy_,:,:,:)

  end subroutine user_calc_sources


  !=====================================================================
  subroutine user_sources
    use ModMain, ONLY: iTest, jTest, kTest, ProcTest, BlkTest, &
         GLOBALBLK
    use ModVarIndexes
    use ModAdvance, ONLY: State_VGB
    use ModGeometry, ONLY: x_Blk, y_Blk, z_Blk, r_Blk, rMin_Blk
    use ModConst
    use ModPlanetConst
    use ModPhysics
    use ModProcMH, ONLY: iProc

    implicit none

    ! Variables required by this user subroutine
    integer :: i,j,k
    logical :: oktest,oktest_me

    !\
    ! Variable meanings:
    !   Srho: Source terms for the continuity equation
    !   SE,SP: Source terms for the energy (conservative) and presure (primative)
    !          equations
    !   SrhoUx,SrhoUy,SrhoUz:  Source terms for the momentum equation
    !   SBx,SBy,SBz:  Souce terms for the magnetic field equations
    !/

    ! User declared local variables go here

    real :: Ux,Uy,Uz,Usq,Unx,Uny,Unz,Un,Unsq,Urelsq
    real :: Rxy,R1,R2,RO1,RO2,RO3
    real :: Hr1,Hr2,Hr3,Hz
    real :: P0,P1,L0,L
    real :: rhodot0,rhodot,rhodotnorm
    real :: accelerationfactor
    real :: alpha_rec,CXsource,Tfrac 
    real :: alphaCOROTATE,ScorotateUx,ScorotateUy,ScorotateUz,ScorotateE
    real :: uCOR,xHat,yHat,zHat,FdotB,Bmag,sinTheta
    real :: gradP_external,gradP_external_X,gradP_external_Y,gradP_external_Z

    !---------------------------------------------------------------------------
    if (.not. UseMassLoading) RETURN
    if (Rmin_BLK(globalBLK) > 15.0)return

    if(iProc==PROCtest .and. globalBLK==BLKtest)then
       call set_oktest('user_sources',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       ! Calculate only in the region where the contribution is not
       ! negligible.  

       ! Note that densities are given in 1/cm^3, 
       ! Velocities, distances and scale heights are normalized
       ! (unitless).

       if (((Z_BLK(i,j,k,globalBLK) < 6.0) .and.     &
            (Z_BLK(i,j,k,globalBLK) > -6.0)) .and.   &
            (R_Blk(i,j,k,globalBLK) < 15.0) .and.     &
            (R_BLK(i,j,k,globalBLK) > 3.0)) then


          Ux  = State_VGB(rhoUx_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
          Uy  = State_VGB(rhoUy_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
          Uz  = State_VGB(rhoUz_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
          Usq = sqrt(Ux**2 + Uy**2 + Uz**2)


          ! compute the cylindrical radius for later use
          Rxy = sqrt(X_BLK(i,j,k,globalBLK)**2 + &
               Y_BLK(i,j,k,globalBLK)**2)   

          ! The neutral velocity is taken to be the orbital velocity
          ! Note that Gbody/r**2 is a unitless acceleration.  So Gbody/r
          ! is a uniless u**2.  Gbody is negative so that gravity pulls 
          ! inward.  Take abs() here to get magnitude of the orbital 
          ! velocity.

          Unsq = abs(Gbody)/Rxy
          Un = sqrt(unsq)
          Unx = -un*Y_BLK(i,j,k,globalBLK)/Rxy
          Uny =  un*X_BLK(i,j,k,globalBLK)/Rxy
          Unz =  0.0

          Urelsq = (unx-Ux)**2 + (uny-Uy)**2 + (unz-Uz)**2

          ! now code the simple model taken from Schreier, 1998
          ! corresponds to a total production rate of 3.2e28 s^-1
          ! analytically and 3.0e28 s^-1 for a grid of .25 spacing 
          ! when rhodot0=0.5e-3.  Not that is not the peak rhodot	 
          ! but is really a scaling factor.  The real peak is 
          ! closer to 1.97e-2 s^-1cm^-3

          ! Variables preceeded by H are scale heights

          ! set up values for this mass loading model

          rhodot0 = 4.4e-4

          R1  = 5.71
          RO1 = 5.71
          R2  = 5.875
          RO2 = 5.685
          RO3 = 5.9455
          Hr1 = .2067
          Hr2 = .1912
          Hr3 = 0.053173*Rxy + 0.55858

          if (Rxy <= R1) then
             rhodot =  rhodot0*60.0*exp((Rxy-RO1)/Hr1)
          elseif (R1 < Rxy .and. Rxy <=R2) then
             rhodot =  rhodot0*60.0*exp(-(Rxy-RO2)/Hr2)
          else
             rhodot =  rhodot0*19.9*exp(-(Rxy-RO3)/Hr3)
          end if

          ! the scale height is 10 degree or Hz = tan(10*pi/180)*Rxy 
          ! Hz = .176327*Rxy
          ! the scale height is 5 degree or Hz = tan(5*pi/180)*Rxy 
          Hz = .087489*Rxy
          rhodot = rhodot*exp(-(Z_BLK(i,j,k,globalBLK)**2)/Hz**2) 

          ! now get everything normalized - note that rhodot is in 
          ! units particles/cm^3/s.  To get a unitless rhodot we simply 
          ! need to multiply by the mass of a proton times the 
          ! #amu/partical (22.0) and 1.0e6 (to get 1/m^3) and then 
          ! divide by (dimensions(rho)/dimensions(t)) which is 
          ! No2Si_V(UnitRho_)/No2Si_V(UnitT_) to get the normalized
          ! form.

          rhodotNorm =  No2Si_V(UnitT_)/No2Si_V(UnitRho_)

          ! testing only
          Alpha_rec = 0.0
          CXsource = 0.0
          Tfrac = 1.0

          rhodot = 22.0*cProtonMass*1.0e6*rhodot*rhodotNorm

          if (AccelerateMassLoading) then
             accelerationFactor = 10.0
             rhodot   = accelerationFactor*rhodot
          end if

          ! Before loading the source term arrays load source terms 
          ! associated with the corotational electric field acceleration
          ! of the newly created plasma ions. This is the extra term 
          ! currently under investigation and in dispute with referees.

          alphaCOROTATE = 1.0
          sinTheta = Rxy/R_BLK(i,j,k,globalBLK)
          uCOR = alphaCOROTATE*OMEGAbody*R_BLK(i,j,k,globalBLK)*sinTheta
          xHat = -Y_BLK(i,j,k,globalBLK)/Rxy
          yHat =  X_BLK(i,j,k,globalBLK)/Rxy
          ! ScorotateUx = rhodot*(uCOR-un)*xHat
          ! ScorotateUy = rhodot*(uCOR-un)*yHat
          ! ScorotateUz = 0.0
          ! ScorotateE  = rhodot*(uCOR-un)*uCOR

          ScorotateUx = 0.0
          ScorotateUy = 0.0
          ScorotateUz = 0.0
          ScorotateE  = 0.0

          ! now load the source terms into the right hand side

          Srho(i,j,k)   = Srho(i,j,k) +                                      & 
               (rhodot-Alpha_rec*State_VGB(rho_,i,j,k,globalBLK))  
          SrhoUx(i,j,k) = SrhoUx(i,j,k)                                      &
               + (rhodot + CXsource*State_VGB(rho_,i,j,k,globalBLK))*unx &
               - (CXsource + Alpha_rec)*                          &
               State_VGB(rhoUx_,i,j,k,globalBLK)                    &
               + ScorotateUx					        
          SrhoUy(i,j,k) = SrhoUy(i,j,k)                                      &
               + (rhodot + CXsource*State_VGB(rho_,i,j,k,globalBLK))*uny &
               - (CXsource + Alpha_rec)*                          &
               State_VGB(rhoUy_,i,j,k,globalBLK)                    &
               + ScorotateUy					        
          SrhoUz(i,j,k) = SrhoUz(i,j,k)                                      &  
               - (CXsource + Alpha_rec)*                          &
               State_VGB(rhoUz_,i,j,k,globalBLK)                    &
               + ScorotateUz
          SE(i,j,k)     = SE(i,j,k)                                               & 
               + 0.5*(rhodot + CXsource*State_VGB(rho_,i,j,k,globalBLK))*unsq &
               - (CXsource + Alpha_rec)*                               &
               (0.5*usq*State_VGB(rho_,i,j,k,globalBLK)                  &
               + 1.5*State_VGB(p_,i,j,k,globalBLK))                 &
               + 1.5*CXsource*State_VGB(p_,i,j,k,globalBLK)*Tfrac             &
               + ScorotateE  
          SP(i,j,k)     = SP(i,j,k)                                                 &
               + 0.5*(rhodot + CXsource*State_VGB(rho_,i,j,k,globalBLK))*urelsq &
               - 1.5*CXsource*State_VGB(p_,i,j,k,globalBLK)*(1.0-Tfrac)         &
               - 1.5*Alpha_rec*State_VGB(p_,i,j,k,globalBLK)                  

          ! Output to look at rates
          if (oktest_me .and. globalBLK==BLKtest  .and. &
               i==Itest .and. j==Jtest .and. k==Ktest ) then
             write(*,*) '----------Inner Torus Mass Loading Rates-------------------------'
             write(*,'(a,3(1X,E13.5))') 'X, Y, Z:', &
                  X_BLK(i,j,k,globalBLK), Y_BLK(i,j,k,globalBLK), Z_BLK(i,j,k,globalBLK)
             write(*,'(a,5(1X,i6))')    'i,j,k,globalBLK,iProc:', &
                  i,j,k,globalBLK,iProc
             write(*,'(a,3(1X,E13.5))') 'rhodot(unnormalized):', &
                  rhodot/rhodotNorm
             write(*,'(a,4(1X,E13.5))') 'rho(amu/cc),un,usq(km/s):', &
                  State_VGB(rho_,i,j,k,globalBLK)*&
                  No2Io_V(UnitRho_),un*No2Io_V(UnitU_),&
                  usq*No2Io_V(UnitU_)
             write(*,'(a,3(1X,E13.5))') 'rhodot (normalized):', &
                  rhodot
             write(*,'(a,4(1X,E13.5))') 'rho, un, usq (normalized):', &
                  State_VGB(rho_,i,j,k,globalBLK), un, usq
             write(*,*) '-----------------------------------------------------------------'
          end if


          ! now code the simple model of high energy particles
          ! as a pressure source taken from Mauk, 1996,1998.
          ! This is an extra grad P which is then found in both
          ! The energy equation and the momentum equation.  The
          ! term is much like gravity and will basically appear
          ! anywhere that gravity does in the MHD equations.
          ! Note that the pressure is constant along field lines so
          ! that the pressure gradient is perpendicular to field lines.

          P1 = 7.0482E-8
          L0 = 7.5419
          Hr1 = .79922

          L   = R_BLK(i,j,k,globalBLK)**3/Rxy**2;

          gradP_external   = -(2.0*P1*(L-L0)/Hr1**2)*exp(-(L-L0)**2/Hr1**2)
          gradP_external_X = gradP_external* &
               (3*X_BLK(i,j,k,globalBLK)*R_BLK(i,j,k,globalBLK)/Rxy**2 &
               - 2*X_BLK(i,j,k,globalBLK)*R_BLK(i,j,k,globalBLK)**3/Rxy**4)
          gradP_external_Y = gradP_external* &
               (3*Y_BLK(i,j,k,globalBLK)*R_BLK(i,j,k,globalBLK)/Rxy**2 &
               - 2*Y_BLK(i,j,k,globalBLK)*R_BLK(i,j,k,globalBLK)**3/Rxy**4)
          gradP_external_Z = gradP_external* &
               (3*Z_BLK(i,j,k,globalBLK)*R_BLK(i,j,k,globalBLK)/Rxy**2)

          ! now get everything normalized - note that gradP_external is in Pascal/Rj.
          ! To get a unitless gradP_external we simply need to divide 
          ! by No2Si_V(UnitP_)
          ! since the derivative was computed using the unitless Rj.

          gradP_external_X = gradP_external_X/No2Si_V(UnitP_)
          gradP_external_Y = gradP_external_Y/No2Si_V(UnitP_)
          gradP_external_Z = gradP_external_Z/No2Si_V(UnitP_)

          ! Now compute the source terms but do so in a way that ensures
          ! no field aligned forces.  In other words compute
          ! F' = F - F.b

          ! First compute the B unit vector b
          Bmag = sqrt(State_VGB(Bx_,i,j,k,globalBLK)**2 + &
               State_VGB(By_,i,j,k,globalBLK)**2 + &
               State_VGB(Bz_,i,j,k,globalBLK)**2 ) 
          xHat = State_VGB(Bx_,i,j,k,globalBLK)/Bmag
          yHat = State_VGB(By_,i,j,k,globalBLK)/Bmag
          zHat = State_VGB(Bz_,i,j,k,globalBLK)/Bmag

          ! now compute the F.b
          FdotB = gradP_external_x*xHat + &
               gradP_external_y*yHat + &
               gradP_external_z*zHat

          ! gradP_external_X =  gradP_external_X - FdotB*xHat 
          ! gradP_external_Y =  gradP_external_Y - FdotB*yHat 
          ! gradP_external_Z =  gradP_external_Z - FdotB*zHat 

          ! now load the source terms into the right hand side of the momentum and
          ! energy equations.  There are no terms in the other equations.

          SrhoUx(i,j,k) = SrhoUx(i,j,k) - gradP_external_X
          SrhoUy(i,j,k) = SrhoUy(i,j,k) - gradP_external_Y
          SrhoUz(i,j,k) = SrhoUz(i,j,k) - gradP_external_Z
          SE(i,j,k)     = SE(i,j,k) - gradP_external_X*Ux   &
               - gradP_external_Y*Uy   & 
               - gradP_external_Z*UZ 

       end if     ! location test

    end do; end do; end do     ! end the i,j,k loops

  end subroutine user_sources

  !=====================================================================
  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut
    implicit none
  
    integer:: i
    character (len=100) :: NameCommand
    !---------------------------------------------------------------------------
  
    if(iProc==0.and.lVerbose > 0)then
        call write_prefix; write(iUnitOut,*)'User read_input SATURN starts'
    endif
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case("#MASSLOADING")
          call read_var('UseMassLoading',UseMassLoading)
  	  call read_var('MassLoadingRate (#/s)', MassLoadingRate)
          call read_var('DoAccelerateMassLoading',AccelerateMassLoading) 

       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix; write(iUnitOut,*)'User read_input SATURN ends'
          endif
          EXIT
       case default
          if(iProc==0) then
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


end module ModUser

