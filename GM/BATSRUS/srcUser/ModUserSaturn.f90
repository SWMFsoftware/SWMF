!========================================================================
Module ModUser
  use ModVarIndexes, ONLY: rho_, Ux_, Uy_, Uz_,p_,Bx_, By_, Bz_, Energy_, &
       rhoUx_,rhoUy_,rhoUz_
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_calc_sources,               &
       IMPLEMENTED3 => user_set_boundary_cells,         &
       IMPLEMENTED4 => user_face_bcs

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Saturn Magnetosphere, Hansen, Nov, 2006'

  real :: MassLoadingRate
  logical :: UseMassLoading, AccelerateMassLoading

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
    use CON_axes

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
    real :: Rxy,R0,R1,Dtorus
    real :: Hr1,Hr2,Hneutral,Helectron
    real :: mn, N0,nN, nElec
    real :: rhodot0,rhodot,rhodotnorm,CXNorm
    real :: accelerationfactor, sourcefactor
    real :: alpha_rec,CXsource,Tfrac 
    real :: Telectron, EVelectron, LogEVelectron,eta,alphaE
    real :: LTerm,Part1,Part2
    real :: rSaturn, rTitan_Orbit,omegaTitan_orbit

    ! arrays for computing the tilt of the equatorial plane the therefore the tilt
    ! of the mass loading

    real, dimension(3) :: xGSE,xSMG,vGSE,vSMG
    real, dimension(3,3) :: SMG_GSE_mat

    !---------------------------------------------------------------------------
    if(iProc==PROCtest .and. globalBLK==BLKtest)then
       call set_oktest('user_sources',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    rSaturn = rPlanet_I(Saturn_)
    rTitan_Orbit  = rOrbitPlanet_I(Titan_)
    omegaTitan_orbit = 2.0*cPi/OrbitalPeriodPlanet_I(Titan_)


    if (UseMassLoading) then

       ! first rotate the coordinate system.  The mass loading	
       ! functions are centered at the rotational equator.  However,
       ! BATSRUS/GM solves in GSE so the z=0 plane is not the equatorial
       ! plane when doing the "real" saturn.  The fuctions here are coded
       ! assuming that the Z=0 plane is the same as the equatorial plane.
       ! This is the SMG frame (see CON_axes).  So, I need to transform
       ! x,y,z from GSE to SMG.  Do this using CON_planet.

       SMG_GSE_mat = transform_matrix(0.0, 'GSE','SMG')

       if (Rmin_BLK(globalBLK) > 45.0)return

       do k = 1, nK
          do j = 1, nJ
             do i = 1, nI


                ! Do the rotations.  Note that we want to compute the positions
                ! in the SMG coordinate system.  However, we need to compute the
                ! velocities in the GSE coordinate system where GM/BATSRUS works.
                ! this makes life a little complicated.

                xGSE(1) = x_BLK(i,j,k,globalBLK)
                xGSE(2) = y_BLK(i,j,k,globalBLK)
                xGSE(3) = z_BLK(i,j,k,globalBLK)
                xSMG = matmul(SMG_GSE_mat,xGSE)


                ! first calculate the mass loading source terms for the inner
                ! torus of H2O and products.  Althought the rates
                ! used fall of at infinity and could be calculated everywhere, we
                ! will only caluclate in the region where the contribution is not
                ! negligible.  

                ! Note that Velocities, distances and scale heights are normalized 
                ! (unitless).  Density is outlined below

                if (((xSMG(3) < 2.0) .and.    &
                     (xSMG(3) > -2.0)) .and.  &
                     (R_Blk(i,j,k,globalBLK) < 14.0) .and.    &
                     (R_BLK(i,j,k,globalBLK) > 3.0)) then

                   Ux  = State_VGB(rhoUx_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
                   Uy  = State_VGB(rhoUy_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
                   Uz  = State_VGB(rhoUz_,i,j,k,globalBLK)/State_VGB(rho_,i,j,k,globalBLK)
                   Usq = sqrt(Ux**2 + Uy**2 + Uz**2)


                   ! compute the cylindrical radius for later use
                   Rxy = sqrt(xSMG(1)**2 + xSMG(2)**2)   

                   ! The neutral velocity is taken to be the orbital velocity
                   ! Note that Gbody/r**2 is a unitless acceleration.  So Gbody/r
                   ! is a uniless u**2.  Gbody is negative so that gravity pulls 
                   ! inward.  Take abs() here to get magnitude of the orbital 
                   ! velocity.

                   !compute in the SMG frame and then rotate back to the GSE frame
                   ! where we really want the velocities.
                   unsq = abs(Gbody)/Rxy
                   un = sqrt(unsq)
                   unx = -un*xSMG(2)/Rxy
                   uny =  un*xSMG(1)/Rxy
                   unz =  0.0

                   vSMG(1) = unx
                   vSMG(2) = uny
                   vSMG(3) = unz
                   vGSE = matmul(vSMG,SMG_GSE_mat)

                   unx = vGSE(1)
                   uny = vGSE(2)
                   unz = vGSE(3)
                   urelsq = (unx-Ux)**2 + (uny-Uy)**2 + (unz-Uz)**2


                   ! now code the model taken from Richardson, 90,98.  See my
                   ! dissertation on p. 125-129.
                   !
                   ! In the dissertation rhodot0 = 7.2e-4 amu/cm^3/s.  This
                   ! rate corresponds to a total production rate of 8.6E26 s^-1
                   ! analytically.  Using an average mass of 16.6 amu per ion
                   ! the correponding mass production rate is 1.43e28 amu/s.
                   ! Note that rhodot0 is in amu/cm^3/s so it ALREADY TAKES INTO
                   ! ACCOUNT THE 16.6 AMU/ION.  It is already a MASS DENSITY 
                   ! not simply a number density!
                   !
                   ! Since Richardson gives a rate of about 1.3e27 s^-1 we multiply by
                   ! a factor of 2 to get a number close to his.  This means that
                   ! we typically use rhodot0 = 2.0* 7.2e-4 amu/cm^3/s.  Using
                   ! this rate and integrating on different size grids gives the
                   ! the following rates
                   !
                   ! dx		ndot		rhodot
                   !
                   ! 1/16		1.44E27		2.40E28
                   ! 1/8		1.16E27		1.92E28
                   ! 0.195	8.55E26		1.42E28	 
                   ! 1/4		6.45E26		1.07E28
                   ! 
                   ! So, to repeat.  Below, rhodot0 and rhodot are in units of 
                   ! amu/cm^3/s and already take into account the average mass
                   ! of an ion that is being added.  You do not need to multiply by
                   ! 16.6 to get the rhodot.		 
                   !
                   !
                   ! NOTE:  here we introduce a multiplication factor which we will
                   ! use to control the mass addition rate.  A factor of 1.0 corresponds
                   ! to a rate of rhodot0 = 7.2e-4 amu/cm^3/s.  A factor of 2.0 
                   ! obviously corresponds to twice this and therefore corresponds
                   ! to the rates in the table above!
                   !
                   !-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
                   !-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
                   !-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

                   !nominal rate of ~0.9E27 on a typical grid (3/16=0.195)
                   !this corresponds to the above numbers for a 3/16 grid)
                   ! sourcefactor = 2.0

                   !rate of ~0.3*1e28 on a typical grid (3/16=0.195) 
                   ! This number is taken from the most recent Richardson
                   ! paper which I just reviewed (Nov, 2004).  The paper
                   ! gives 1e28 as the neutral source rate of which .3 is
                   ! converted to ions.  The rest escapes the system as
                   ! neutrals.   The 3.333 comes from .3E28/.9E27
                   ! sourcefactor = 2.0*3.333333333

                   !rate of ~1e28 on a typical grid  (3/16=0.195) 
                   ! This number is the current neutral source rate being
                   ! used by Richardson as well as the Cassini people.  As
                   ! such this is the largest mass loading rate you might expect.
                   ! The 11.1 comes from 1e28/.9e27
                   ! sourcefactor = 2.0*11.1

                   !rate which can be adjusted by the user from the input file
                   ! The user inputs the MassLoadingRate in #ofparticles/cm^3.
                   ! We still assume that the average mass is 16.6 .  Now, so that
                   ! we don't have to change anything below, here we caluclate
                   ! the sourcefacter from the MassLoadingRate.  As pointed out
                   ! above a sourcefactor=2.0  corresponds to a rate of ~0.9E27.
                   ! So, to calcuate the sourcefactor we use
                   !    sourcefactor_A/2.0  =  MassLoadingRate/0.9E27
                   ! so
                   !    sourcefactor_A = 2.0*MassLoadingRate/0.9e27
                   ! or
                   !    sourcefactor_A = 2.22222E-27 * MassLoadingRate
                   !
                   ! Just to test.  If MassLoadingRate = 1.8E27 then
                   !
                   !    sourcefactor_A = 4.0
                   !
                   ! which is twice the value for a rate of 0.9E27! Good.
                   !
                   sourcefactor = 2.222222E-27 * MassLoadingRate





                   !-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
                   !-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
                   !-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


                   ! Variables preceeded by H are scale heights


                   ! If we are using accelerated mass loading then calculate 
                   ! the factor by which we should increase the rates

                   ! set up values for this mass loading model

                   Hneutral  = .45
                   Helectron = .6 + .2*(Rxy-3.0)

                   rhodot0 = sourcefactor*7.2E-4

                   R0 = 5.35
                   Hr1 = 1.2
                   Hr2 = 2.25 + 0.075*(Rxy-R0)

                   if (Rxy < R0) then
                      rhodot =  rhodot0*exp(-(Rxy-R0)**2/Hr1**2)
                   else
                      rhodot =  rhodot0*exp(-(Rxy-R0)**2/Hr2**2)
                   end if

                   rhodot = rhodot*exp(-(xSMG(3)**2)/Hneutral**2)&
                        *exp(-(xSMG(3)**2)/Helectron**2)


                   ! get the CX friction source terms
                   !
                   ! The charge exchange friction is also calculated from 
                   ! Richardson in my dissertation.  The charge exchange 
                   ! friction depends on the neutral density and the plasma
                   ! density.  We have coded eta so that the neutral density
                   ! dependence is built in.  So, if the neutral density is
                   ! doubled the eta should be doubled.  Since above we change	
                   ! the mass loading rate by a source factor which is general
                   ! because of a neutral rate, we will use the same factor here
                   ! to control the eta scaling.

                   R0 = 4.6
                   n0 = 0.4*1200
                   Hr1 = .8
                   Hr2 = 1.2 + 0.035*Rxy

                   if (Rxy < R0) then
                      Cxsource =  exp((Rxy-R0)/Hr1)
                   else
                      CXsource =  exp(-(Rxy-R0)/Hr2)
                   end if

                   CXsource = 0.6*((-0.0330 + 0.1636*Rxy)*1E-8)              &
                        *n0*CXsource*exp(-(xSMG(3)**2)/ &
                        Hneutral**2)

                   CXsource = sourcefactor*CXsource

                   ! get the electron recombination terms.  First get the 
                   ! electron density.  This is gotten by first converting from 
                   ! unitless rho to amr/cm^3 using No2Io_V(UnitRho_).  Then, by 
                   ! dividing by the average amu per particle. In that case we 
                   ! use 16.6 for O, OH, H2O.  Note that this term does not
                   ! depend on the neutrals in any way so we do not need to use the
                   ! sourcefactor to scale as we raise and lower the mass loading
                   ! rate.


                   nElec = State_VGB(rho_,i,j,k,globalBLK)*No2Io_V(UnitRho_)/16.6

                   ! We now need the electron temperature to calculate some 
                   ! of the ionization rates.  The electron temperature that 
                   ! comes out of this calculation have units of K (kelvin).
                   ! This is an electron temperature calculation that takes
                   ! Te = 1/15 Tp (Tp = Ti+Te). The EVelectron is the electron 
                   ! temperature converted to an energy in electron volts (eV).
                   ! This is used in a formula to get the recombination rate 
                   ! in s.

                   Tfrac = 1.0/15.0;

                   Telectron = (16.6*cProtonMass/cBoltzmann)*(State_VGB(p_,i,j,k,globalBLK)*     &
                        No2Si_V(UnitP_))/(State_VGB(rho_,i,j,k,globalBLK)*No2Si_V(UnitRho_))* &
                        Tfrac
                   EVelectron = cBoltzmann/cElectronCharge*Telectron
                   LogEVelectron = Log(EVelectron)

                   Alpha_rec = exp(-17.6168-1.3982*LogEVelectron+0.1439* &
                        LogEVelectron**2)*nElec
                   if (Telectron < 1.0) Alpha_rec = 2.234E-8*nElec
                   if (Telectron > 100) Alpha_rec = 7.552E-10*nElec

                   ! now get everything normalized - note that rhodot is in 
                   ! units amu/cm^3/s already (the mass is taken into account
                   ! above.  To get a unitless rhodot we simply need to multiply 
                   ! by the mass of a proton times the 1.0e6 (to get 1/m^3) and then 
                   ! divide by (dimensions(rho)/dimensions(t)) which is 
                   ! No2Si_V(UnitRho_)/No2Si_V(UnitT_) to get the normalized (dimensionless) 
                   ! form.

                   rhodotNorm =  No2Si_V(UnitT_)/No2Si_V(UnitRho_)
                   CXNorm = No2Si_V(UnitT_)

                   rhodot = cProtonMass*1.0e6*rhodot*rhodotNorm
                   CXsource = CXsource*CXNorm 
                   Alpha_rec = Alpha_rec*CXNorm

                   if (AccelerateMassLoading) then
                      accelerationFactor = 10.0
                      rhodot   = accelerationFactor*rhodot
                      CXsource = accelerationFactor*CXsource
                   end if

                   !                   ! testing only
                   !                   Alpha_rec = 0.0
                   !                   CXsource = 0.0

                   ! now load the source terms into the right hand side

                   Srho(i,j,k)   = Srho(i,j,k)                                        &  
                        + (rhodot-Alpha_rec*State_VGB(rho_,i,j,k,globalBLK))  
                   SrhoUx(i,j,k) = SrhoUx(i,j,k)                                      &
                        + (rhodot + CXsource*State_VGB(rho_,i,j,k,globalBLK))*unx &
                        - (CXsource + Alpha_rec)*                          &
                        State_VGB(rhoUx_,i,j,k,globalBLK)                    
                   SrhoUy(i,j,k) = SrhoUy(i,j,k)                                      &
                        + (rhodot + CXsource*State_VGB(rho_,i,j,k,globalBLK))*uny &
                        - (CXsource + Alpha_rec)*                          &
                        State_VGB(rhoUy_,i,j,k,globalBLK)                    
                   SrhoUz(i,j,k) = SrhoUz(i,j,k)                                      &  
                        - (CXsource + Alpha_rec)*                          &
                        State_VGB(rhoUz_,i,j,k,globalBLK)                    
                   SE(i,j,k)     = SE(i,j,k)                                               & 
                        + 0.5*(rhodot + CXsource*State_VGB(rho_,i,j,k,globalBLK))*unsq &
                        - (CXsource + Alpha_rec)*                               &
                        (0.5*usq*State_VGB(rho_,i,j,k,globalBLK)                  &
                        + 1.5*State_VGB(p_,i,j,k,globalBLK))                  &
                        + 1.5*CXsource*State_VGB(p_,i,j,k,globalBLK)*Tfrac             
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
                      write(*,'(a,3(1X,E13.5))') 'Telectron, EVelectron, LogEVelectron:', &
                           Telectron, EVelectron, LogEVelectron
                      write(*,'(a,3(1X,E13.5))') 'rhodot, CXsource, Alpha_rec (unnormalized):', &
                           rhodot/rhodotNorm, CXsource/CXNorm,Alpha_rec/CXNorm
                      write(*,'(a,4(1X,E13.5))') 'rho(amu/cc),nElec(1/cc),un,usq(km/s):', &
                           State_VGB(rho_,i,j,k,globalBLK)*No2Io_V(UnitRho_), nElec, un*No2Io_V(UnitU_), usq*No2Io_V(UnitU_)
                      write(*,'(a,3(1X,E13.5))') 'rhodot, CXsource, Alpha_rec (normalized):', &
                           rhodot, CXsource,Alpha_rec
                      write(*,'(a,4(1X,E13.5))') 'rho, nElec, un, usq (normalized):', &
                           State_VGB(rho_,i,j,k,globalBLK), nElec/No2Io_V(UnitRho_), un, usq
                      write(*,*) '-----------------------------------------------------------------'
                   end if

                end if ! inner torus - location test

                ! Now calculate the source terms for the neutral nitrogen torus
                ! centered around Titan.

                ! Calculate the radial distance from the center of the torus,
                ! if inside, compute the mass loading, if outside do nothing.


                Dtorus = sqrt( (sqrt(xSMG(1)**2+   &
                     xSMG(2)**2) - &
                     Rtitan_orbit/Rsaturn)**2 +          &
                     xSMG(3)**2)

                if (Dtorus < 10.0) then

                   usq = (State_VGB(rhoUx_,i,j,k,globalBLK)**2 +   &
                        State_VGB(rhoUy_,i,j,k,globalBLK)**2 +   &
                        State_VGB(rhoUz_,i,j,k,globalBLK)**2 ) / &
                        (State_VGB(rho_,i,j,k,globalBLK)**2)


                   ! compute the cylindrical radius for later use
                   Rxy = sqrt(xSMG(1)**2 + xSMG(2)**2)   

                   ! The neutral torus is taken to rotate ridgidly with Titan
                   unsq = (Rxy**2)*((OMEGAtitan_orbit*No2Si_V(UnitT_))**2)
                   un = sqrt(unsq)

                   rhodotNorm =  No2Si_V(UnitT_)/No2Si_V(UnitRho_)

                   ! See my dissertation to find the functional form of what
                   ! is being used below and how it was derived.  This is
                   ! calculated somewhat differently than the icy satellites.
                   ! Here we clearly calculate the nuetral density  and then
                   ! get rhodot (amu/cm^3/s) my using all the correct factors
                   ! including the average mass of the ion (14.0 nitrogen).
                   ! note that Dtorus**2/2.0 = Dtorus**2/(sqrt(2)**2) with a
                   ! resulting scale height of sqrt(2).  The resulting mass
                   ! loading rates for different grid sizes are:
                   !
                   ! dx		ndot		rhodot
                   !
                   ! 1/16		5.46E25		7.64E26
                   ! 1/8		5.21E25		7.29E26
                   ! 0.195	4.87E25		6.82E26	 
                   ! 1/4		4.61E25		6.45E26
                   ! 
                   !

                   mn = 14.0
                   nN = 10.*exp(-Dtorus**2/2.0)
                   rhodot = (1.0E+6*cProtonMass*(14.0*nN/(3.0E7)))*rhodotNorm
                   eta = 2.14287E-4

                   ! this is an electron temperature calculation that takes
                   ! Te = 1/2 Tp = 1/2 (Ti+Te).  

                   Tfrac = 1.0/2.0;

                   Telectron = (mn*cProtonMass/cBoltzmann)*(State_VGB(p_,i,j,k,globalBLK)*       &
                        No2Si_V(UnitP_))/(State_VGB(rho_,i,j,k,globalBLK)*No2Si_V(UnitRho_))* &
                        Tfrac
                   EVelectron = cBoltzmann/cElectronCharge*Telectron
                   LogEVelectron = Log(EVelectron)

                   ! Define alpha
                   if (EVelectron < 200.) then
                      alphaE = 7.E-7*sqrt(300./EVelectron)
                   else
                      alphaE = 2.342*7.E-7*   &
                           EVelectron**(0.2553-0.1633*log10(EVelectron))
                   end if

                   ! note that the term in the second bracket is the
                   ! electron density calculated by assuming the average
                   ! mass of the implanted ions. Futher note that the 0.1
                   ! is the number density in cm^-3.

                   Lterm = alphaE*(No2Si_V(UnitT_))*   &
                        (No2Io_V(UnitRho_)*State_VGB(rho_,i,j,k,globalBLK)/mn)

                   part1 = rhodot + rhodot*eta*State_VGB(rho_,i,j,k,globalBLK)
                   part2 = rhodot*eta + Lterm

                   ! now load the source terms into the right hand side

                   Srho(i,j,k)   = Srho(i,j,k)   + (rhodot - Lterm*State_VGB(rho_,i,j,k,globalBLK))

                   SrhoUx(i,j,k) = SrhoUx(i,j,k) + (part1*(-OMEGAtitan_orbit*No2Si_V(UnitT_)*  &
                        Y_BLK(i,j,k,globalBLK)) -  &
                        part2*State_VGB(rhoUx_,i,j,k,globalBLK))  
                   SrhoUy(i,j,k) = SrhoUy(i,j,k) + (part1*(OMEGAtitan_orbit*No2Si_V(UnitT_)*   &
                        X_BLK(i,j,k,globalBLK)) -   &
                        part2*State_VGB(rhoUy_,i,j,k,globalBLK))   
                   SrhoUz(i,j,k) = SrhoUz(i,j,k) + (part1*(0.0) - &
                        part2*State_VGB(rhoUz_,i,j,k,globalBLK))  
                   SE(i,j,k)     = SE(i,j,k)     + (part1*0.5*unsq -   &
                        part2*((0.5*State_VGB(rho_,i,j,k,globalBLK)*usq)+ &
                        ((3./2.)*State_VGB(p_,i,j,k,globalBLK))))  

                end if    ! end calculation of source terms inside the torus.


             end do     ! end the i loop
          end do     ! end the j loop
       end do     ! end the k loop

    end if     ! usemassloading test

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

  !============================================================================

  subroutine user_set_boundary_cells(iBlock)
    use ModGeometry,      ONLY: ExtraBc_, IsBoundaryCell_GI, x_Blk, x2
    use ModBoundaryCells, ONLY: SaveBoundaryCells

    implicit none

    integer, intent(in):: iBlock

    character (len=*), parameter :: Name='user_set_boundary_cells'
    !--------------------------------------------------------------------------
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = x_Blk(:,:,:,iBlock) > x2

    if(SaveBoundaryCells) return
    call stop_mpi('Set SaveBoundaryCells=.true. in PARAM.in file')

  end subroutine user_set_boundary_cells

  !============================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize, ONLY: x_
    use ModVarIndexes, ONLY: nVar, Bx_, Bz_
    use ModFaceBc, ONLY: TimeBc, iFace, jFace, kFace, FaceCoords_D
    use ModSolarwind, ONLY: get_solar_wind_point
    use ModB0, ONLY: B0_DX

    real, intent(out):: VarsGhostFace_V(nVar)
    !------------------------------------------------------------------------
    call get_solar_wind_point(TimeBc, FaceCoords_D(x_), VarsGhostFace_V)
    VarsGhostFace_V(Bx_:Bz_) = VarsGhostFace_V(Bx_:Bz_) - B0_DX(:,iFace, jFace, kFace)

  end subroutine user_face_bcs

end module ModUser

