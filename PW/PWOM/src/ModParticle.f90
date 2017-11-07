Module ModParticle
  Use ModRandomNumber, ONLY: random_real
  implicit none
  
  private
  
  !basic particle type
  type particle
     integer :: iSpecies
     integer :: iCell ! cell index for particle
     real    :: vpar, vperp !velocity in [cm/s]
     real    :: Alt ! position in configurational space [cm]
     logical :: IsOpen   ! defines if particle is in domain or index is avail.
     real    :: NumPerParticle !the weight of the particle (how many real 
                                !particles one macro particle represents).
  end type particle

  !type for holding pointer for particles of a specific species in a cell
  type particleCellSpecies
     type(particle),pointer :: Particle 
  end type particleCellSpecies

  !type for holding array of particles (mostly for bury and disinter)
  type particleHolder
     integer :: nParticleOnLine
     type(particle), allocatable :: SavedParticles_I(:) 
  end type particleHolder

  !array that holds the particles to push
  type(particle), target,allocatable :: Particles_I(:)
  
  !pointer to reference the particles of a particular species in a particular 
  !cell.
  type(particleCellSpecies), allocatable:: SortParticles_III(:,:,:)
  !number of particles of particluar type in cell for a given species
  integer,allocatable :: nSortedParticle_II(:,:)

  !total number of particles in the simulation
  integer :: nParticle
  
  !maximum number of particles before throwing an error (not set yet)
  integer :: MaxParticles
  
  !how many lines do we have
  integer,public :: nLine=1
  integer, public, allocatable :: iLineGlobal_I(:)
  integer :: iLineCurrent
  
  !Frequency of outputs
  real :: DtSaveProfile=300
  real :: DtSaveDF=60.0

  real :: DtSplitJoin=60.0
  
  !How many and which altitudes should the DF be saved
  !integer,parameter :: nSaveDfAlts=3
  !integer::iAltsDF_I(nSaveDfAlts)=(/26,76,151/)
  !real,allocatable :: SaveDfAlts_I(:)

  integer,parameter :: nSaveDfAlts=14
  integer::iAltsDF_I(nSaveDfAlts)=(/1,13,25,37,49,61,73,85,97,109,121,133,145,157/)

!  integer,parameter :: nSaveDfAlts=1
!  integer::iAltsDF_I(nSaveDfAlts)=(/1/)
  !hold the buried lines
  type(particleHolder),allocatable :: BuriedParticles_I(:) 
  
  ! grid variables
  integer :: nAlt    ! points on grid
  integer :: nCells  ! number of cells including ghost cells
  real,allocatable  :: Alt_G(:)  ! position of cell center with ghost [cm]
  real,allocatable  :: dAlt_G(:) ! width of cell with ghost [cm]
  real,allocatable  :: AltBot_F(:)  ! position of cell face [cm]
  real,allocatable  :: AltTop_F(:)  ! position of cell face [cm]
  real,allocatable  :: Volume_G(:)  ! Cell volume
  real :: AreaInlet ! area of the botom of the first computational cell

  !array to hold lower ghostcell
  real,allocatable :: DensityBC_I(:),VelocityBC_I(:),TemperatureBC_I(:)

  !variables for overlap region
  real,allocatable :: DensityOverlap_IC(:,:), VelocityOverlap_IC(:,:),&
       TemperatureOverlap_IC(:,:), FluidFrac_C(:)
  integer,public :: nOverlap=3
  logical,public :: UseOverlapRegion=.false.

  ! array relating cell index to particle index
  !integer :: iCell_II(:,:)

  ! Electric field variables [statV/cm]
  real,allocatable  :: Efield_G(:) 
  
  ! Time step for moving particles [s]
  real :: DtMove=1.0,DtMoveMax=1.0
  
  !time for colliding particles
  real :: DtCollide
  integer,parameter :: nCollide=1 !subcycle collisions
  
  !particle variables
  real,allocatable :: Mass_I(:)
  integer :: O_, H_, He_,H3_,H2_
  real,allocatable::ReducedMass_II(:,:)

  !variable for WPI
  real,allocatable :: Dperp_I(:),Dexp_I(:)
  logical,public :: UseWPI=.false.
  character(len=100),public :: TypeWPI='Barakat'
  real,public :: FracLeftHand=0.125
  
  real,public :: SpectralIndexAur=1.7
  real,public :: rWaveRefAur=7375e5
  real,public :: E2waveRefAur=1.2e-6 !V^2 m^-2 Hz^-1
  real,public :: fWaveRefAur=5.6    !Hz

  real,public :: SpectralIndexCap=1.7
  real,public :: rWaveRefCap=7375e5
  real,public :: E2waveRefCap=1.2e-6 !V^2 m^-2 Hz^-1
  real,public :: fWaveRefCap=5.6    !Hz

  
  !should the output associated with particles be verbose
  logical,public :: IsVerboseParticle = .false.

  integer :: nSpecies

  character(len=2),allocatable :: NameSpecies_I(:) 
  
  integer :: iSeed = 1
  
  ! use a fixed number of particles per cell
  integer, allocatable :: nParticlePerCell_I(:)



  real, parameter :: cBoltzmannCGS = 1.3807E-16

  real :: Time=0.0

  ! public routines
  public :: init_particle
  public :: put_to_particles
  public :: get_from_particles
  public :: bury_line
  public :: disinter_line
  public :: run_particles
  public :: write_restart_particle,read_restart_particle
  ! public unit tests
  public :: test_sample
  public :: test_split_join
  public :: test_pusher
  public :: test_coulomb_collision
  public :: test_wpi
  public :: test_combine_fluid_particle
contains

  !============================================================================
  subroutine init_particle(nAltIn,AltMin,AltMax,TypeGrid)
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I, rPlanet_I
    integer,intent(in) :: nAltIn
    real,   intent(in) :: AltMin, AltMax
    character(len=100),intent(in):: TypeGrid
    real, parameter :: cGramsPerAMU=1.66054e-24
    real :: rPlanetCM, alpha, Area, AreaBot,AreaTop, dAlt
    integer :: iAlt,iSpecies,jSpecies
    real, parameter :: cMtoCm=1e2

    real :: dFrac
    !--------------------------------------------------------------------------
    !set number of species and mass
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       nSpecies=3
       O_=1
       H_=2
       He_=3

       allocate(Mass_I(nSpecies))
       Mass_I(O_) = cGramsPerAMU*15.994
       Mass_I(H_) = cGramsPerAMU*1.00797
       Mass_I(He_) = cGramsPerAMU*4.0026
              
       !Mass_I(He_) = cGramsPerAMU*4.0
       allocate(NameSpecies_I(nSpecies))
       NameSpecies_I(O_)='O_'
       NameSpecies_I(H_)='H_'
       NameSpecies_I(He_)='He'

       !init the target particle per cell numbers
       allocate(nParticlePerCell_I(nSpecies))
       nParticlePerCell_I(O_)=20000
       nParticlePerCell_I(H_)=3000
       nParticlePerCell_I(He_)=1000

       
       !set WPI variables
       allocate(Dperp_I(nSpecies))
       allocate(Dexp_I(nSpecies))
    case('JUPITER')
       nSpecies=3
       H3_=1
       H_=2
       H2_=3

       allocate(Mass_I(nSpecies))
       Mass_I(H3_) = cGramsPerAMU*3.0237
       Mass_I(H_) = cGramsPerAMU*1.00797
       Mass_I(H2_) = cGramsPerAMU*2.0159
              
       !Mass_I(He_) = cGramsPerAMU*4.0
       allocate(NameSpecies_I(nSpecies))
       NameSpecies_I(H3_)='H3'
       NameSpecies_I(H_)='H_'
       NameSpecies_I(H2_)='H2'

       !init the target particle per cell numbers
       allocate(nParticlePerCell_I(nSpecies))
       nParticlePerCell_I(H3_)=5000
       nParticlePerCell_I(H_)=5000
       nParticlePerCell_I(H2_)=5000

       
       !set WPI variables
       allocate(Dperp_I(nSpecies))
       allocate(Dexp_I(nSpecies))
       

    case DEFAULT
       call con_stop('particles not for planet')
    end select
    !allocate arrays to hold ghost cell moments
    allocate(DensityBC_I(nSpecies),VelocityBC_I(nSpecies),&
         TemperatureBC_I(nSpecies))
    
    if(UseOverlapRegion) then
       !allocate arrays for overlap region
       allocate(DensityOverlap_IC(nSpecies,nOverlap),&
            VelocityOverlap_IC(nSpecies,nOverlap),&
            TemperatureOverlap_IC(nSpecies,nOverlap),&
            FluidFrac_C(nOverlap))
       
       ! set the fraction of the fluid in each overlap cell starting with mostly
       ! fluid at bottom to mostly particle at top
       dFrac=1.0/(real(nOverlap)+1.0)
       
       do iAlt=1,nOverlap 
          FluidFrac_C(iAlt)=1.0-dFrac*real(iAlt)
       enddo
    endif
    
    !save the reduced mass which is useful for collisions
    allocate(ReducedMass_II(nSpecies,nSpecies))
    do iSpecies=1,nSpecies
       do jSpecies=1,nSpecies
          ReducedMass_II(iSpecies,jSpecies)=Mass_I(iSpecies)*Mass_I(jSpecies)&
               /(Mass_I(iSpecies)+Mass_I(jSpecies))
       enddo
    enddo



    ! Initialize the grid, set size and allocate arrays. 
    ! note -1 and nAlt+2 ghost cells only needed for E field interpolation.
    ! 0 and nAlt ghost cells needed for sampling.
    nAlt = nAltIn
    nCells=nAlt+4
    
    allocate(Alt_G(-1:nAlt+2))
    allocate(dAlt_G(-1:nAlt+2))
    allocate(AltBot_F(-1:nAlt+2))
    allocate(AltTop_F(-1:nAlt+2))
    allocate(Volume_G(-1:nAlt+2))
    
    allocate(Efield_G(-1:nAlt+2))
    
    ! allocate array to tell number of particles of a given type in each cell
    write(*,*) nAlt+1
    allocate(nSortedParticle_II(nSpecies,0:nAlt+1))


    ! based on grid type set the grid
    select case(TypeGrid)
    case('Uniform')
       dAlt=(AltMax-AltMin)/(nAlt-1)
       dAlt_G=dAlt

       ! set area coef if A=alpha r^3 assuming crossection of 1 cm2 at base
       rPlanetCM=rPlanet_I(Planet_)*cMtoCm

       alpha= 1.0/(rPlanet_I(Planet_)*cMtoCm+AltMin)**3
       
       do iAlt=-1,nAlt+2
          ! set location of cell center and top and bottom faces
          Alt_G(iAlt)=AltMin-dAlt+dAlt*iAlt
          AltBot_F(iAlt)=Alt_G(iAlt)-0.5*dAlt
          AltTop_F(iAlt)=Alt_G(iAlt)+0.5*dAlt
          
          ! calculate the cell volume assuming calculated by assuming 
          !an area function in the form A=alpha r^3
          !and then assuming each cell is a truncated cone.
          AreaBot = alpha*(rPlanetCM+AltBot_F(iAlt))**3
          AreaTop = alpha*(rPlanetCM+AltTop_F(iAlt))**3
          
          Volume_G(iAlt) = 1.0/3.0 * dAlt *&
               ( AreaBot + AreaTop + (AreaBot*AreaTop)**0.5 )
       enddo

       ! set the area at the inlet of the simulation domain. Needed for 
       ! injecting particles when not using ghost cell filling
       AreaInlet = alpha*(rPlanetCM+AltBot_F(1))**3
       

    case DEFAULT
       call con_stop('Gridtype not recognized')
    end select
  end subroutine init_particle
  
  !============================================================================
  ! push guiding center of particles
  subroutine push_guiding_center
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I,mPlanet_I,rPlanet_I    
    use ModConst,ONLY:cGravitation
    use ModInterpolate, only: linear
    integer :: iParticle,iAlt
    real :: rCoord, Efield,acceleration,gravity,AltStart
    
    !interpolation variables
    real :: xAlt, Dx1, Dx2
    integer :: iMin,iMax

    real, parameter :: cMtoCm=1e2
    real, parameter :: cElecChargeCGS= 4.80320425e-10 !statcoulombs
    !--------------------------------------------------------------------------

    !$OMP PARALLEL DO PRIVATE(AltStart,rCoord,gravity,Efield,acceleration,iAlt,&
    !$OMP xAlt,iMin,iMax,Dx1,Dx2)
    do iParticle=1,nParticle
       AltStart=Particles_I(iParticle)%Alt
       !set the radial distance and the gravitational acceleration 
       ! at particle location. note need rCoord in cm but modplanetconst in SI
       rCoord=rPlanet_I(Planet_)*cMtoCm+AltStart
       gravity=cMtoCm**3*cGravitation*mPlanet_I(Planet_)/rCoord**2
       

       ! interpolate electric field to particle (note ModInterpolate cannot be 
       !used since it is not safe for multiple OpenMp threads. 
       !note current interpolation only works for uniform grid
       xAlt=(AltStart-Alt_G(-1))/dAlt_G(-1)-1.0
       iMin=floor(xAlt)
       iMax=ceiling(xAlt)
       !set interpolation weights
       Dx1=xAlt-iMin ; Dx2=1.0-Dx1
       !interpolate
       Efield=Dx2*Efield_G(iMin)+Dx1*Efield_G(iMax)
       
       !Efield= linear(Efield_G(:),-1,nAlt+1,Particles_I(iParticle)%Alt,Alt_G)
       !Efield=0.     
       !determine acceleration (mirror force-gravity-eField)
       acceleration=1.5*Particles_I(iParticle)%vperp**2/rCoord - gravity&
            +Efield*(cElecChargeCGS/Mass_I(Particles_I(iParticle)%iSpecies))
       
       
       !from initial velocity and acceleration update state
       Particles_I(iParticle)%vpar=Particles_I(iParticle)%vpar&
            +acceleration*DtMove
       Particles_I(iParticle)%Alt=AltStart+Particles_I(iParticle)%vpar*DtMove
       Particles_I(iParticle)%vperp=&
            Particles_I(iParticle)%vperp&
            *(AltStart/Particles_I(iParticle)%Alt)**1.5


       if (Particles_I(iParticle)%Alt<AltBot_F(1) .or. &
            Particles_I(iParticle)%Alt>AltTop_F(nAlt)) then
          Particles_I(iParticle)%IsOpen = .true.
       endif
       
       !Assign cell index to particle
       CELL_ASSIGN: do iAlt=0,nAlt+1
          if(Particles_I(iParticle)%Alt>AltBot_F(iAlt) &
               .and. Particles_I(iParticle)%Alt<AltTop_F(iAlt)) then
             Particles_I(iParticle)%iCell=iAlt
             exit CELL_ASSIGN
          endif
       end do CELL_ASSIGN
    end do
    !$OMP END PARALLEL DO
    
  end subroutine push_guiding_center

  !============================================================================
  ! push guiding center of particles using rk4 method
  subroutine push_guiding_center_rk4
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I,mPlanet_I,rPlanet_I    
    use ModConst,ONLY:cGravitation
    use ModInterpolate, only: linear
    integer :: iParticle,iAlt
    real :: rCoord, Efield,acceleration,gravity,AltStart
    
    !interpolation variables
    real :: xAlt, Dx1, Dx2
    integer :: iMin,iMax

    !variables for the rk4 method
    real :: dAlt1,dAlt2,dAlt3,dAlt4,dVpar1,dVpar2,dVpar3,dVpar4,Vperp,VperpStart
    real :: VparStart, Mass
    real, parameter ::OneOverSix=0.166666666666666666666666666666666666
    
    !conversion variables
    real, parameter :: cMtoCm=1e2
    real, parameter :: cElecChargeCGS= 4.80320425e-10 !statcoulombs
    
    !coef for calculation of gravity (to reduce repeated multiplies)
    real :: GravCoef
    !--------------------------------------------------------------------------

    GravCoef = cMtoCm**3*cGravitation*mPlanet_I(Planet_)

    !$OMP PARALLEL DO PRIVATE(AltStart,rCoord,gravity,Efield,acceleration,iAlt,&
    !$OMP xAlt,iMin,iMax,Dx1,Dx2,&
    !$OMP dAlt1,dAlt2,dAlt3,dAlt4,dVpar1,dVpar2,dVpar3,dVpar4,Vperp,VperpStart,&
    !$OMP VparStart, Mass)
    
    do iParticle=1,nParticle
       AltStart=Particles_I(iParticle)%Alt
       VperpStart=Particles_I(iParticle)%vperp
       VparStart=Particles_I(iParticle)%vpar
       Mass = Mass_I(Particles_I(iParticle)%iSpecies)
       
       !set the radial distance and the gravitational acceleration 
       ! at particle location. note need rCoord in cm but modplanetconst in SI
       rCoord=rPlanet_I(Planet_)*cMtoCm+AltStart
       gravity=GravCoef/rCoord**2
       

       ! interpolate electric field to particle (note ModInterpolate cannot be 
       !used since it is not safe for multiple OpenMp threads. 
       !note current interpolation only works for uniform grid
       xAlt=(AltStart-Alt_G(-1))/dAlt_G(-1)-1.0
       iMin=floor(xAlt)
       iMax=ceiling(xAlt)
       
       !set interpolation weights
       Dx1=xAlt-iMin ; Dx2=1.0-Dx1
       !interpolate
       Efield=Dx2*Efield_G(iMin)+Dx1*Efield_G(iMax)
       
       !Efield= linear(Efield_G(:),-1,nAlt+1,Particles_I(iParticle)%Alt,Alt_G)
       !Efield=0.     
       !determine acceleration (mirror force-gravity-eField)
       acceleration=1.5*VperpStart**2/rCoord - gravity&
            +Efield*(cElecChargeCGS/Mass)
       
       !save RK4 step 1
       dVpar1 = acceleration*DtMove
       dAlt1 = VparStart*DtMove

       
       !start to get RK4 step2
       Vperp = VperpStart&
            *(AltStart/(AltStart+0.5*dAlt1))**1.5
       !set the radial distance and the gravitational acceleration 
       ! at particle location. note need rCoord in cm but modplanetconst in SI
       rCoord=rPlanet_I(Planet_)*cMtoCm+(AltStart+0.5*dAlt1)
       gravity=GravCoef/rCoord**2
       

       ! interpolate electric field to particle (note ModInterpolate cannot be 
       !used since it is not safe for multiple OpenMp threads. 
       !note current interpolation only works for uniform grid
       xAlt=((AltStart+0.5*dAlt1)-Alt_G(-1))/dAlt_G(-1)-1.0
       iMin=floor(xAlt)
       iMax=ceiling(xAlt)
       !set interpolation weights
       Dx1=xAlt-iMin ; Dx2=1.0-Dx1
       !interpolate
       Efield=Dx2*Efield_G(iMin)+Dx1*Efield_G(iMax)
       
       !Efield= linear(Efield_G(:),-1,nAlt+1,Particles_I(iParticle)%Alt,Alt_G)
       !Efield=0.     
       !determine acceleration (mirror force-gravity-eField)
       acceleration=1.5*Vperp**2/rCoord - gravity&
            +Efield*(cElecChargeCGS/Mass)

       !save RK4 step 2
       dVpar2 = acceleration*DtMove
       dAlt2 = (VparStart+0.5*dVpar1)*DtMove
       
       !get RK4 step 3
       Vperp = VperpStart&
            *(AltStart/(AltStart+0.5*dAlt2))**1.5

       !set the radial distance and the gravitational acceleration 
       ! at particle location. note need rCoord in cm but modplanetconst in SI
       rCoord=rPlanet_I(Planet_)*cMtoCm+(AltStart+0.5*dAlt2)
       gravity=GravCoef/rCoord**2
       

       ! interpolate electric field to particle (note ModInterpolate cannot be 
       !used since it is not safe for multiple OpenMp threads. 
       !note current interpolation only works for uniform grid
       xAlt=((AltStart+0.5*dAlt2)-Alt_G(-1))/dAlt_G(-1)-1.0
       iMin=floor(xAlt)
       iMax=ceiling(xAlt)
       !set interpolation weights
       Dx1=xAlt-iMin ; Dx2=1.0-Dx1
       !interpolate
       Efield=Dx2*Efield_G(iMin)+Dx1*Efield_G(iMax)
       
       !Efield= linear(Efield_G(:),-1,nAlt+1,Particles_I(iParticle)%Alt,Alt_G)
       !Efield=0.     
       !determine acceleration (mirror force-gravity-eField)
       acceleration=1.5*Vperp**2/rCoord - gravity&
            +Efield*(cElecChargeCGS/Mass)


       !save RK4 step 3
       dVpar3 = acceleration*DtMove
       dAlt3 = (VparStart+0.5*dVpar2)*DtMove

       !get RK4 step 4
       Vperp = VperpStart&
            *(AltStart/(AltStart+dAlt3))**1.5

       !set the radial distance and the gravitational acceleration 
       ! at particle location. note need rCoord in cm but modplanetconst in SI
       rCoord=rPlanet_I(Planet_)*cMtoCm+(AltStart+dAlt3)
       gravity=GravCoef/rCoord**2
       

       ! interpolate electric field to particle (note ModInterpolate cannot be 
       !used since it is not safe for multiple OpenMp threads. 
       !note current interpolation only works for uniform grid
       xAlt=((AltStart+dAlt3)-Alt_G(-1))/dAlt_G(-1)-1.0
       iMin=floor(xAlt)
       iMax=ceiling(xAlt)
       !set interpolation weights
       Dx1=xAlt-iMin ; Dx2=1.0-Dx1
       !interpolate
       Efield=Dx2*Efield_G(iMin)+Dx1*Efield_G(iMax)
       
       !Efield= linear(Efield_G(:),-1,nAlt+1,Particles_I(iParticle)%Alt,Alt_G)
       !Efield=0.     
       !determine acceleration (mirror force-gravity-eField)
       acceleration=1.5*Vperp**2/rCoord - gravity&
            +Efield*(cElecChargeCGS/Mass)


       !save RK4 step 4
       dVpar4 = acceleration*DtMove
       dAlt4 = (VparStart+dVpar3)*DtMove


       ! from initial velocity and acceleration update state using 
       !rk4 method
       Particles_I(iParticle)%vpar=VparStart&
            +(dVpar1+2.0*dVpar2+2.0*dVpar3+dVpar4)*OneOverSix
       Particles_I(iParticle)%Alt=AltStart&
            +(dAlt1+2.0*dAlt2+2.0*dAlt3+dAlt4)*OneOverSix
       Particles_I(iParticle)%vperp=&
            VperpStart&
            *(AltStart/Particles_I(iParticle)%Alt)**1.5


       !check is particle leaves computational domain
       if (Particles_I(iParticle)%Alt<AltBot_F(1) .or. &
            Particles_I(iParticle)%Alt>AltTop_F(nAlt)) then
          Particles_I(iParticle)%IsOpen = .true.
       endif
       
       !Assign cell index to particle
       CELL_ASSIGN: do iAlt=0,nAlt+1
          if(Particles_I(iParticle)%Alt>AltBot_F(iAlt) &
               .and. Particles_I(iParticle)%Alt<AltTop_F(iAlt)) then
             Particles_I(iParticle)%iCell=iAlt
             exit CELL_ASSIGN
          endif
       end do CELL_ASSIGN
    end do
    !$OMP END PARALLEL DO
    
  end subroutine push_guiding_center_rk4



  !=============================================================================
  ! sample maxwellian in cell
  subroutine sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature)
    use ModNumConst, ONLY: cPi,cTwoPi
    ! index of cell to sample
    integer,intent(in):: iCell

    ! Species to Sample
    integer,intent(in):: iSpecies
    
    ! parameters of maxwellian for sampling
    real, intent(in) :: density ![cm-3]
    real, intent(in) :: uBulk ![cm/s]
    real, intent(in) :: Temperature ![k]
    
    real :: RandNum1,RandNum2
    real :: uTherm,uMax,uMin,uRange, uRand, uMostProb,RandNum, PitchAngle
    real :: uMagRel, uPar,uPerp
    real :: fmax, ftemp
    !variable to hold new particle info before inserting it to particles array
    type(particle),allocatable :: NewParticle_I(:)
    integer :: nNew,iParticle, iParticleOld
    logical :: DoTest=.true.
    
    !variables for assignment of particles
    type(particle),allocatable ::ParticlesOld_I(:)
    integer :: nParticleOld,nAvail,i,j
    real :: NumPerParticle
    integer,allocatable :: IndexAvail_I(:)
    !--------------------------------------------------------------------------
    nNew=nParticlePerCell_I(iSpecies)

    ! Find number of particles represented by a macro particle by !
    !taking Ntrue=density*volume and dividing by nParticlesPerCell
    NumPerParticle=Density * Volume_G(iCell)/real(nNew)
    !write(*,*) nNew

    ! allocate array to hold new particles
    allocate(NewParticle_I(nNew))
   
    !use accept-reject algorithm for distribution sampling
    iParticle=1
    sample_loop: do while (iParticle<=nNew)
       
       ! randomly choose velocity in range
       RandNum1=random_real(iSeed)
       RandNum2=random_real(iSeed)

       uPar=uBulk+sqrt(-2.0*(cBoltzmannCGS/Mass_I(iSpecies))*Temperature*log(RandNum1))*cos(cTwoPi*RandNum2)
       uPerp=sqrt(-2.0*(cBoltzmannCGS/Mass_I(iSpecies))*Temperature*log(RandNum1))*sin(cTwoPi*RandNum2)
          
       NewParticle_I(iParticle)%vpar =uPar
       NewParticle_I(iParticle)%vperp=abs(uPerp)
       
       !randomly place particle in cell
       RandNum=random_real(iSeed)
       NewParticle_I(iParticle)%Alt=RandNum*dAlt_G(iCell)+AltBot_F(iCell)
       
       !assign remaining particle properties
       NewParticle_I(iParticle)%iSpecies=iSpecies
       NewParticle_I(iParticle)%iCell=iCell
       NewParticle_I(iParticle)%IsOpen=.false.
       NewParticle_I(iParticle)%NumPerParticle=NumPerParticle
       

       !increment particle counter
       iParticle=iParticle+1
    end do sample_loop
!    write(*,*) 'nNew',nNew
!    if(DoTest) then 
!       call plot_distribution_cell(iSpecies,iCell,nNew,NewParticle_I)
!       return
!    endif

    ! Now newly sampled particles need to be put into main particle array
    !on the first call there is no Particles array allocated so allocate
    if (.not.allocated(Particles_I))then
       nParticle=nNew
       allocate(Particles_I(nParticle))
       Particles_I=NewParticle_I
    else
       ! when Particles_I already allocated (usual case) then calculate the
       !number of open spots (nAvail) then save Particles_I and allocate a new 
       !Particles_I array to save the good particles and the newly created ones
       allocate(IndexAvail_I(nParticle))
       where(Particles_I%IsOpen)
          IndexAvail_I=1
       elsewhere
          IndexAvail_I=0
       end where
       nAvail=sum(IndexAvail_I)
       deallocate(IndexAvail_I)
       
       !save old particle array information
       nParticleOld=nParticle
       allocate(ParticlesOld_I(nParticleOld))
       ParticlesOld_I=Particles_I
       deallocate(Particles_I)
       
       !allocate new Particles_I array
       nParticle=nParticleOld-nAvail+nNew
       allocate(Particles_I(nParticle))

       !now fill new Particle_I array with old particles that are not open and 
       !the newly created particles
       iParticle=1
       do iParticleOld=1,nParticleOld
          if(.not.ParticlesOld_I(iParticleOld)%IsOpen) then
             Particles_I(iParticle)=ParticlesOld_I(iParticleOld)
             iParticle=iParticle+1
          endif
       end do
       Particles_I(nParticleOld-nAvail+1:nParticle)=NewParticle_I

       deallocate(ParticlesOld_I)
    endif

    !deallocate to save memory
    deallocate(NewParticle_I)
  end subroutine sample_maxwellian_cell_boxmuller

  !=============================================================================
  real function maxwellian(mass,vel,Temp)
    use ModNumConst, ONLY: cTwoPi  
    real,intent(in) :: mass,vel,Temp
    
    !maxwellian = sqrt((mass/(cTwoPi*cBoltzmannCGS*Temp))**3)*2.0*cTwoPi*vel**2&
    !     *exp(-mass*vel**2/(2.0*cBoltzmannCGS*Temp))

    maxwellian = (mass/(cTwoPi*cBoltzmannCGS*Temp))**1.5&
         *exp(-mass*vel**2/(2.0*cBoltzmannCGS*Temp))
  end function maxwellian
  !=============================================================================
  ! plot a altitude profile of the integrated moments for a given species
  subroutine plot_profile(iSpecies)
    use ModNumConst, ONLY: cPi,cTwoPi
    use ModPlotFile,   ONLY: save_plot_file
    integer, intent(in) :: iSpecies
    
    integer :: iCell
    integer :: nParticleInCell
    
    real :: density,uBulkPar,uBulkPerp,Pressure,Temp, uTherm
    real :: dVel, vParMin, vParMax, vPerpMin, vPerpMax
    real :: Tpar,Tperp,Hpar,Hperp
    integer, parameter :: nVel = 100
    real :: vPar_C(nVel), vPerp_C(nVel)
    real, allocatable   :: Coord_I(:), PlotState_IV(:,:)
    !grid parameters
    integer, parameter :: nDim =1

    !plot variables
    character(len=100) :: NamePlot

    real,parameter :: cCmToKm=1.0e-5
    character(len=*),parameter :: NameHeader='moments'
    character(len=5) :: TypePlot='ascii'
    logical, save :: IsFirstCall=.true.
    logical,allocatable,save :: IsFirstCall_I(:)

    !indices for plots
    integer, parameter :: iDen_=1,iVel_=2,iPres_=3, iTemp_=4,iTpar_=5,iTperp_=6,&
         iHpar_=7,iHperp_=8
    integer, parameter :: nVar = 8

    !name for plot output
    character(len=150)::NameProfilePlotVar=&
         'Alt[km] n[cm-3] u[km/s] p T[k] Tpar[k] '&
         //'Tperp[k] Hpar Hperp g r'
    !---------------------------------------------------------------------------
    if (IsFirstCall) then
       allocate(IsFirstCall_I(nLine))
       IsFirstCall_I(:)=.true.
       IsFirstCall = .false.
    endif
    
    allocate(Coord_I(0:nAlt),PlotState_IV(0:nAlt,nVar))    
    do iCell=0,nAlt
       Coord_I(iCell)=Alt_G(iCell)*1e-5
       nParticleInCell=nSortedParticle_II(iSpecies,iCell)
       if(nParticleInCell==0) then
          PlotState_IV(iCell,iDen_)  = 0.0
          PlotState_IV(iCell,iVel_) = 0.0
          PlotState_IV(iCell,iPres_) = 0.0
          PlotState_IV(iCell,iTemp_) = 0.0
          PlotState_IV(iCell,iTpar_) = 0.0
          PlotState_IV(iCell,iTperp_) = 0.0
          PlotState_IV(iCell,iHpar_) = 0.0
          PlotState_IV(iCell,iHperp_) = 0.0
       else
          ! get moments in cell. note weighted calculation falls apart 
          ! next to ghost cell so revert to basic calcuation
          if (iCell<=1 .or. iCell==nAlt) then
             call calc_moments_cell(iSpecies,iCell,&
                  density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp,Hpar,Hperp)
          else
             call calc_moments_cell_weighted(iSpecies,iCell,&
                  density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp,Hpar,Hperp)
          endif
          PlotState_IV(iCell,iDen_)  = density
          PlotState_IV(iCell,iVel_) = uBulkPar*cCmToKm
          PlotState_IV(iCell,iPres_) = Pressure
          PlotState_IV(iCell,iTemp_) = Temp
          PlotState_IV(iCell,iTpar_) = Tpar
          PlotState_IV(iCell,iTperp_) = Tperp
          PlotState_IV(iCell,iHpar_) = Hpar
          PlotState_IV(iCell,iHperp_) = Hperp
       end if
    end do
    !Plot 
    !write(NamePlot,"(a)") 'Profile.out'
    write(NamePlot,"(a,i5.5,a)") &
         'PW/plots/Profile_',iLineGlobal_I(iLineCurrent),'.out'

    if(IsFirstCall_I(iLineCurrent)) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NameProfilePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
       IsFirstCall_I(iLineCurrent) = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NameProfilePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/))
    endif
    deallocate(Coord_I, PlotState_IV)
  end subroutine plot_profile
  
  !============================================================================
  subroutine plot_distribution_cell(iSpecies,iCell)
    use ModNumConst, ONLY: cPi,cTwoPi
    use ModPlotFile,   ONLY: save_plot_file
    integer,intent(in) :: iSpecies,iCell
    integer :: nParticleInCell
    
    real :: density,uBulkPar,uBulkPerp,Pressure,Temp, uTherm
    real :: dVel, vParMin, vParMax, vPerpMin, vPerpMax
    real :: Tpar,Tperp,Hpar,Hperp
    integer, parameter :: nVel = 100
    real :: vPar_C(nVel), vPerp_C(nVel)
    integer :: iVel, iParticle, iVpar,iVperp
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    !grid parameters
    integer, parameter :: nDim =2, Vperp_=1,Vpar_=2,nVar=1, PSD_=1
    !plot variables
    character(len=100) :: NamePlot
    character(len=100),parameter :: NamePlotVar='Vperp[cm/s] Vpar[cm/s] Particles g r'
    character(len=*),parameter :: NameHeader='distribution function'
    character(len=5) :: TypePlot='ascii'
    real :: TrueParticles

    integer,save :: iCounter=0
    !---------------------------------------------------------------------------
    TrueParticles=0
        
    nParticleInCell=nSortedParticle_II(iSpecies,iCell)
    if(nParticleInCell==0)return
    allocate(Coord_DII(nDim,nVel,nVel),PlotState_IIV(nVel,nVel,nVar))
    
    ! get moments in cell so we can calculate thermal velocity
    call calc_moments_cell(iSpecies,iCell,&
         density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp,Hpar,Hperp)
    
    write(*,*) 'density,uBulkPar,uBulkPerp,Pressure,Temp'&
         ,density,uBulkPar,uBulkPerp,Pressure,Temp
    ! calculate thermal velocity

    uTherm=sqrt(8.0*cBoltzmannCGS*Temp/Mass_I(iSpecies)/cPi)

    ! discretize velocity space, center around bulk velocity to 
    !5 uTherm in every direction with grid size .1 uTherm
    vParMin=uBulkPar-5.0*uTherm
    vParMax=uBulkPar+5.0*uTherm
    vPerpMin=uBulkPerp-5.0*uTherm
    vPerpMax=uBulkPerp+5.0*uTherm
    dVel=0.1*uTherm
    do iVel=1,nVel
       vPar_C(iVel)=vParMin+iVel*dVel
       vPerp_C(iVel)=vPerpMin+iVel*dVel
    enddo

    do iVpar=1,nVel
       do iVperp=1,nVel
          Coord_DII(Vperp_,iVperp,iVpar)=VPerp_C(iVperp)
          Coord_DII(Vpar_,iVperp,iVpar)=VPar_C(iVpar)
       enddo
    enddo
    
    PlotState_IIV=0.0
    !Sort particles into bins
    do iParticle=1,nParticleInCell
       iVpar =floor((SortParticles_III(iParticle,iSpecies,iCell)%Particle%vpar &
            -vParMin )/dVel)
       iVperp=floor((SortParticles_III(iParticle,iSpecies,iCell)%Particle%vperp&
            -vPerpMin)/dVel)
       
       if (iVpar>0 .and.iVpar<=nVel .and.iVperp>0 .and.iVperp<=nVel) then
          
          PlotState_IIV(iVperp,iVpar,PSD_)=&
               PlotState_IIV(iVperp,iVpar,PSD_)&
               +SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
          TrueParticles=TrueParticles&
               +SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
       endif
    enddo
    
    ! divide the particles in each velocity bin by the density and the 
    !velocity space cell volume. d3v = dvpar*dvperp*2pi*vperp assuming 
    !rings in velocity space around the vperp axis.
    if (nParticleInCell>0) then
       PlotState_IIV(:,:,PSD_) = PlotState_IIV(:,:,PSD_)&
            /TrueParticles!/(cTwoPi*Coord_DII(Vperp_,:,:)**2*dVel**2)
    else
       PlotState_IIV(:,:,PSD_) = 0.0
    endif

    !Plot 
    write(NamePlot,"(a,a,a,i5.5,a,i5.5,a)") &
         'PW/plots/DistFunc_',NameSpecies_I(iSpecies),'Alt',&
         floor(Alt_G(iCell)*1e-5),'km_iLine',iLineGlobal_I(iLineCurrent),'.out'

    if(iCounter<(nSaveDfAlts*nLine)) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       iCounter=iCounter+1
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
    endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_distribution_cell 
  !=============================================================================
  subroutine plot_distribution_cell_orig(iSpecies,iCell,nParticleInCell,CellParticle_I)
    use ModNumConst, ONLY: cPi,cTwoPi
    use ModPlotFile,   ONLY: save_plot_file
    integer,intent(in) :: iSpecies,iCell,nParticleInCell
    type(particle), intent(in) :: CellParticle_I(nParticleInCell)
    
    real :: density,uBulkPar,uBulkPerp,Pressure,Temp, uTherm
    real :: dVel, vParMin, vParMax, vPerpMin, vPerpMax
    integer, parameter :: nVel = 100
    real :: vPar_C(nVel), vPerp_C(nVel)
    integer :: iVel, iParticle, iVpar,iVperp
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    !grid parameters
    integer, parameter :: nDim =2, Vperp_=1,Vpar_=2,nVar=1, PSD_=1
    !plot variables
    character(len=100) :: NamePlot
    character(len=100),parameter :: NamePlotVar='Vperp[cm/s] Vpar[cm/s] Particles g r'
    real :: TrueParticles
    character(len=*),parameter :: NameHeader='distribution function'
    character(len=5) :: TypePlot='ascii'
    logical, save :: IsFirstCall=.true.
    !---------------------------------------------------------------------------
    TrueParticles=0
    
    allocate(Coord_DII(nDim,nVel,nVel),PlotState_IIV(nVel,nVel,nVar))

    ! get moments in cell so we can calculate thermal velocity
    call calc_moments_cell_orig(iSpecies,iCell,nParticleInCell,&
         CellParticle_I,density,uBulkPar,uBulkPerp,Pressure,Temp)
    
!    write(*,*) 'density,uBulkPar,uBulkPerp,Pressure,Temp'&
!         ,density,uBulkPar,uBulkPerp,Pressure,Temp
    ! calculate thermal velocity

    uTherm=sqrt(8.0*cBoltzmannCGS*Temp/Mass_I(iSpecies)/cPi)

    ! discretize velocity space, center around bulk velocity to 
    !5 uTherm in every direction with grid size .1 uTherm
    vParMin=uBulkPar-5.0*uTherm
    vParMax=uBulkPar+5.0*uTherm
    vPerpMin=uBulkPerp-5.0*uTherm
    vPerpMax=uBulkPerp+5.0*uTherm
    dVel=0.1*uTherm
    do iVel=1,nVel
       vPar_C(iVel)=vParMin+iVel*dVel
       vPerp_C(iVel)=vPerpMin+iVel*dVel
    enddo

    do iVpar=1,nVel
       do iVperp=1,nVel
          Coord_DII(Vperp_,iVperp,iVpar)=VPerp_C(iVperp)
          Coord_DII(Vpar_,iVperp,iVpar)=VPar_C(iVpar)
       enddo
    enddo
    
    PlotState_IIV=0.0
    !Sort particles into bins
    do iParticle=1,nParticleInCell
       iVpar =floor((CellParticle_I(iParticle)%vpar -vParMin )/dVel)
       iVperp=floor((CellParticle_I(iParticle)%vperp-vPerpMin)/dVel)
       
       if (iVpar>0 .and.iVpar<=nVel .and.iVperp>0 .and.iVperp<=nVel) then
          PlotState_IIV(iVperp,iVpar,PSD_)=&
               PlotState_IIV(iVperp,iVpar,PSD_)&
               +CellParticle_I(iParticle)%NumPerParticle

          TrueParticles=TrueParticles+CellParticle_I(iParticle)%NumPerParticle
       endif
    enddo
    
    ! divide the particles in each velocity bin by the density and the 
    !velocity space cell volume. d3v = dvpar*dvperp*2pi*vperp assuming 
    !rings in velocity space around the vperp axis.
    PlotState_IIV(:,:,PSD_) = PlotState_IIV(:,:,PSD_)&
         /TrueParticles!/(cTwoPi*Coord_DII(Vperp_,:,:)**2*dVel**2)

    !Plot 
    write(NamePlot,"(a,i4.4,a)") 'DistFuncCell_',iCell,'.out'
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=0,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
    endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_distribution_cell_orig

  !=============================================================================
  subroutine calc_moments_cell(iSpecies,iCell,&
       density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp,Hpar,Hperp)
    integer,intent(in) :: iSpecies,iCell
    integer :: nParticleInCell,iParticle
    real :: uParTmp,uPerpTmp
    real :: Ppar,Pperp
    real :: weight
    real, intent(out) :: density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp
    real, intent(out) :: Hpar,Hperp
    real :: TrueParticles
    !--------------------------------------------------------------------------
    
    TrueParticles=0
    nParticleInCell=nSortedParticle_II(iSpecies,iCell)
    
    if(nParticleInCell==0)then
       density=0.0
       Pressure=0.0
       uBulkPar=0.0
       uBulkPerp=0.0
       Temp=0.0
       Tpar=0.0
       Tperp=0.0
       Hpar=0.0
       Hperp=0.0
       return
    endif

    density=0.0
    !density is particles in cell over volume
    do iParticle=1,nParticleInCell
       density=density+&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
       TrueParticles=TrueParticles&
            +SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
    enddo
    density=density/Volume_G(iCell)

    !bulk velocity is average of velocity
    uBulkPar=0.0
    do iParticle=1,nParticleInCell
       !write(*,*) iParticle,nParticleInCell
       uParTmp =SortParticles_III(iParticle,iSpecies,iCell)%Particle%vpar&
            *SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
       uBulkPar  = uBulkPar + uParTmp        
    end do
    uBulkPar=uBulkPar/TrueParticles
!    uBulkPerp = sum(CellParticle_I(:)%vperp)/nParticleInCell
    uBulkPerp = 0.0
    
    !from eq6.12 Gombosi Gas Kinetic Theory book
    Pressure=0.0
    Ppar=0.0
    Pperp=0.0
    Hpar =0.0
    Hperp=0.0
    do iParticle=1,nParticleInCell
       uParTmp =SortParticles_III(iParticle,iSpecies,iCell)%Particle%vpar
       uPerpTmp=SortParticles_III(iParticle,iSpecies,iCell)%Particle%vperp
       weight=SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
       Pressure = Pressure+weight*(&
            (uParTmp-uBulkPar)**2&
            +2.0*(uPerpTmp-uBulkPerp)**2)
       Ppar = Ppar&
            +weight*((uParTmp-uBulkPar)**2)
       
       Pperp = Pperp&
            +weight*((uPerpTmp-uBulkPerp)**2)

       !from eq 37-38 in demars et al 1992 
       Hpar = Hpar+weight*(&
            (uParTmp-uBulkPar)*(uParTmp-uBulkPar)**2)
       Hperp= Hperp+weight*(&
            (uParTmp-uBulkPar)*(uPerpTmp)**2)
    enddo
    Pressure=Pressure*Mass_I(iSpecies)*density/(3.0*TrueParticles)
    Ppar=Ppar*Mass_I(iSpecies)*density/(TrueParticles)
    Pperp=Pperp*Mass_I(iSpecies)*density/(TrueParticles)
    Hpar=Hpar*Mass_I(iSpecies)*density/(TrueParticles)
    Hperp=0.5*Hperp*Mass_I(iSpecies)*density/(TrueParticles)
    
    !get temperature from ideal gas law P=nkT
    Temp = Pressure/(density*cBoltzmannCGS)
    Tpar = Ppar/(density*cBoltzmannCGS)
    Tperp = Pperp/(density*cBoltzmannCGS)


  end subroutine calc_moments_cell
  !============================================================================

!  !=============================================================================
!  ! same as calc_moments_cell, but if statistics are poor combine with 
!  ! neighbor cells until reasonable statistics are reached
!  subroutine calc_moments_cell_smooth(iSpecies,iCell,&
!       density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp)
!    integer,intent(in) :: iSpecies,iCell
!    integer :: nParticleInCell,iParticle,nNumPerParticle,iAlt
!    integer :: MinCell, MaxCell
!    real :: uParTmp,uPerpTmp
!    real :: Ppar,Pperp,Volume
!    real, intent(out) :: density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp
!    integer, parameter :: MinNumParticle=50.0
!    !--------------------------------------------------------------------------
!    nNumPerParticle=nNumPerParticle_I(iSpecies)    
!
!    nParticleInCell=nSortedParticle_II(iSpecies,iCell)
!    Volume=Volume_G(iCell)
!    ! if not enough particles in cell, then add in cells around iCell until 
!    !statistics are good.
!    MinCell=iCell
!    MaxCell=iCell
!    if(nParticleInCell<MinNumParticle)then
!       do while (nParticleInCell<MinNumParticle)
!          !set min and max cell
!          MinCell=max(MinCell-1,1)
!          MaxCell=min(MaxCell+1,nAlt)
!          !find numer of particles in combined cell
!          nParticleInCell=0
!          Volume=0.0
!          do iAlt=MinCell,MaxCell
!             nParticleInCell=nParticleInCell+nSortedParticle_II(iSpecies,iAlt)
!             Volume=Volume+Volume_G(iAlt)
!          end do
!       end do
!    endif
!
!    !density is particles in cell over volume
!    density=nNumPerParticle/Volume*nParticleInCell
!
!    !bulk velocity is average of velocity
!    uBulkPar=0.0
!    do iAlt=MinCell,MaxCell
!       do iParticle=1,nSortedParticle_II(iSpecies,iAlt)
!          !write(*,*) iParticle,nParticleInCell
!          uParTmp =SortParticles_III(iParticle,iSpecies,iAlt)%Particle%vpar
!          uBulkPar  = uBulkPar + uParTmp        
!       end do
!    enddo
!    uBulkPar=uBulkPar/nParticleInCell
!    !    uBulkPerp = sum(CellParticle_I(:)%vperp)/nParticleInCell
!    uBulkPerp = 0.0
!    
!    !from eq6.12 Gombosi Gas Kinetic Theory book
!    Pressure=0.0
!    Ppar=0.0
!    Pperp=0.0
!    do iAlt=MinCell,MaxCell
!       do iParticle=1,nSortedParticle_II(iSpecies,iAlt)
!          uParTmp =SortParticles_III(iParticle,iSpecies,iAlt)%Particle%vpar
!          uPerpTmp=SortParticles_III(iParticle,iSpecies,iAlt)%Particle%vperp
!          Pressure = Pressure+&
!               (uParTmp-uBulkPar)**2&
!               +2.0*(uPerpTmp-uBulkPerp)**2
!          Ppar = Ppar&
!               +(uParTmp-uBulkPar)**2
!          
!          Pperp = Pperp&
!               +(uPerpTmp-uBulkPerp)**2
!       enddo
!    enddo
!    Pressure=Pressure*Mass_I(iSpecies)*density/(3.0*nParticleInCell)
!    Ppar=Ppar*Mass_I(iSpecies)*density/(nParticleInCell)
!    Pperp=Pperp*Mass_I(iSpecies)*density/(nParticleInCell)
!    
!    !get temperature from ideal gas law P=nkT
!    Temp = Pressure/(density*cBoltzmannCGS)
!    Tpar = Ppar/(density*cBoltzmannCGS)
!    Tperp = Pperp/(density*cBoltzmannCGS)
!    
!  end subroutine calc_moments_cell_smooth
  !============================================================================
  subroutine calc_moments_cell_orig(iSpecies,iCell,nParticleInCell,CellParticle_I,&
       density,uBulkPar,uBulkPerp,Pressure,Temp)
    integer,intent(in) :: iSpecies,iCell, nParticleInCell
    type(particle), intent(in) :: CellParticle_I(nParticleInCell)
    real, intent(out) :: density,uBulkPar,uBulkPerp,Pressure,Temp
    real :: TrueParticles
    integer :: iParticle
    !--------------------------------------------------------------------------
    TrueParticles=0

    density=0.0
    !density is particles in cell over volume
    do iParticle=1,nParticleInCell
       density=density+&
            CellParticle_I(iParticle)%NumPerParticle
       TrueParticles=TrueParticles&
            +CellParticle_I(iParticle)%NumPerParticle
    enddo
    density=density/Volume_G(iCell)

    uBulkPar=0.0
    do iParticle=1,nParticleInCell
       uBulkPar = uBulkPar+CellParticle_I(iParticle)%vpar&
            *CellParticle_I(iParticle)%NumPerParticle
    enddo
    
    uBulkPar  = uBulkPar/TrueParticles
       
!    uBulkPerp = sum(CellParticle_I(:)%vperp)/nParticleInCell
    uBulkPerp = 0.0
    
    !from eq6.12 Gombosi Gas Kinetic Theory book 
    do iParticle=1,nParticleInCell
       Pressure = Pressure+CellParticle_I(iParticle)%NumPerParticle*(&
            (CellParticle_I(iParticle)%vpar-uBulkPar)**2&
            +2.0*(CellParticle_I(iParticle)%vperp-uBulkPerp)**2)
    enddo
    Pressure=Pressure*Mass_I(iSpecies)*density/(3.0*TrueParticles)

    !get temperature from ideal gas law P=nkT
    Temp = Pressure/(density*cBoltzmannCGS)

  end subroutine calc_moments_cell_orig

  !=============================================================================
  subroutine calc_moments_cell_weighted(iSpecies,iCell,&
       density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp,Hpar,Hperp)
    integer,intent(in) :: iSpecies,iCell
    integer :: nParticleInCell,iParticle,iAlt
    real :: uParTmp,uPerpTmp,weight,TotalWeight,Alt
    real :: Ppar,Pperp
    real, intent(out) :: density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp
    real, intent(out) :: Hpar,Hperp
    real :: TrueParticles
    !--------------------------------------------------------------------------
    TrueParticles=0
    TotalWeight=0.0
    density=0.0
    Pressure=0.0
    Ppar=0.0
    Pperp=0.0
    Hpar=0.0
    Hperp=0.0
    uBulkPar=0.0
    uBulkPerp=0.0
    Temp=0.0
    Tpar=0.0
    Tperp=0.0
    ALT_LOOP:do iAlt=iCell-1,iCell+1
       if (iAlt<0 .or. iAlt>nAlt+1) cycle ALT_LOOP
       nParticleInCell=nSortedParticle_II(iSpecies,iAlt)
       if (nParticleInCell==0) cycle ALT_LOOP
       PARTICLE_LOOP:do iParticle=1,nParticleInCell
          !get weight
          Alt=SortParticles_III(iParticle,iSpecies,iAlt)%Particle%Alt
          !check if no weight contributes
          if(Alt<Alt_G(iCell-1) .or. Alt>Alt_G(iCell+1)) cycle PARTICLE_LOOP
          weight=1.0-abs(Alt_G(iCell)-Alt)/dAlt_G(iCell) !assumes uniform grid!
          
!          write(*,*) iParticle
          !density is particles in cell over volume
!          if (iParticle==1005)then
!             write(*,*) ' '
!             write(*,*) iCell,iSpecies
!             write(*,*)SortParticles_III(iParticle,iSpecies,iCell)%Particle%Alt
!             write(*,*)SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
!          endif
             
          density=density+weight&
               *SortParticles_III(iParticle,iSpecies,iAlt)%Particle&
               %NumPerParticle/Volume_G(iCell)
          TrueParticles=TrueParticles+weight*&
            SortParticles_III(iParticle,iSpecies,iAlt)%Particle%NumPerParticle

          !bulk velocity is average of velocity
          uParTmp =SortParticles_III(iParticle,iSpecies,iAlt)%Particle%vpar
          uBulkPar  = uBulkPar + weight*uParTmp&
              *SortParticles_III(iParticle,iSpecies,iAlt)%Particle&
              %NumPerParticle
          
          TotalWeight=TotalWeight+weight
       end do PARTICLE_LOOP        
    end do ALT_LOOP
    
    if (TotalWeight<1e-10) then
       return
    endif
    uBulkPar=uBulkPar/TrueParticles
    uBulkPerp = 0.0

    !now loop again to get pressure  and temperature
    ALT_LOOP2:do iAlt=iCell-1,iCell+1
       if (iAlt<0 .or. iAlt>nAlt+1) cycle ALT_LOOP2
       nParticleInCell=nSortedParticle_II(iSpecies,iAlt)
       if (nParticleInCell==0) cycle ALT_LOOP2
       PARTICLE_LOOP2:do iParticle=1,nParticleInCell
          !get weight
          Alt=SortParticles_III(iParticle,iSpecies,iAlt)%Particle%Alt
          !check if no weight contributes
          if(Alt<Alt_G(iCell-1) .or. Alt>Alt_G(iCell+1)) cycle PARTICLE_LOOP2
          weight=1.0-abs(Alt_G(iCell)-Alt)/dAlt_G(iCell) !assumes uniform grid!

          !from eq6.12 Gombosi Gas Kinetic Theory book
          uParTmp =SortParticles_III(iParticle,iSpecies,iAlt)%Particle%vpar
          uPerpTmp=SortParticles_III(iParticle,iSpecies,iAlt)%Particle%vperp
          Pressure = Pressure+&
               ((uParTmp-uBulkPar)**2&
               +2.0*(uPerpTmp-uBulkPerp)**2)*weight*&
               SortParticles_III(iParticle,iSpecies,iAlt)%Particle%NumPerParticle
          Ppar = Ppar&
               +(uParTmp-uBulkPar)**2*weight*&
               SortParticles_III(iParticle,iSpecies,iAlt)%Particle%NumPerParticle
          Pperp = Pperp&
               +(uPerpTmp-uBulkPerp)**2*weight*&
               SortParticles_III(iParticle,iSpecies,iAlt)%Particle%NumPerParticle
          !from eq 37-38 in demars et al 1992
          Hpar = Hpar+weight*(&
               (uParTmp-uBulkPar)*(uParTmp-uBulkPar)**2)*&
               SortParticles_III(iParticle,iSpecies,iAlt)%Particle%NumPerParticle
          Hperp = Hperp+weight*(&
               (uParTmp-uBulkPar)*(uPerpTmp)**2)*&
               SortParticles_III(iParticle,iSpecies,iAlt)%Particle%NumPerParticle

       end do PARTICLE_LOOP2
    end do ALT_LOOP2

    Pressure=Pressure*Mass_I(iSpecies)*density/(3.0*TrueParticles)
    Ppar=Ppar*Mass_I(iSpecies)*density/(TrueParticles)
    Pperp=Pperp*Mass_I(iSpecies)*density/(TrueParticles)
    Hpar=Hpar*Mass_I(iSpecies)*density/(TrueParticles)
    Hperp=0.5*Hperp*Mass_I(iSpecies)*density/(TrueParticles)
    
    !get temperature from ideal gas law P=nkT
    Temp = Pressure/(density*cBoltzmannCGS)
    Tpar = Ppar/(density*cBoltzmannCGS)
    Tperp = Pperp/(density*cBoltzmannCGS)



  end subroutine calc_moments_cell_weighted
  
  !============================================================================
  ! sort particles so they can be referenced by species, cell, and particle 
  ! index in cell. Note that our sorting uses a derive type holding pointers 
  ! so all sorted data is referenced ultimately to the Particles_I array.
  ! so modification to the data in the sort should change the data in the 
  ! main particle array.
  subroutine sort_particles
    
    !integer and cell indices for which to assign to pointer
    integer :: iCell, iSpecies
    
    integer,allocatable :: Index_I(:),iParticleCell_II(:,:)
    integer :: iParticle,MaxSortedParticle
    integer :: nThreads,iThread
    
    integer, allocatable :: nOffset_TIG(:,:,:),nParticleCell_TIG(:,:,:)

!$    integer,external :: OMP_get_max_threads,OMP_get_thread_num, &
!$         OMP_get_num_threads

    !--------------------------------------------------------------------------




    !\
    ! the new openmp version
    !/
    
    allocate(iParticleCell_II(nSpecies,0:nAlt+1))
    iParticleCell_II=1

    ! go through particle list once to find how many particle in a given cell
    !are on a given thread. This is used to build an offset list

    !get the number of threads
    nThreads=1
    !$ nThreads = OMP_get_max_threads()
    
    ! allocate an offset array
    allocate(nOffset_TIG(0:nThreads-1,nSpecies,0:nAlt+1))
    nOffset_TIG=0
    
    ! allocate counter array of number of particles on a given thread 
    allocate(nParticleCell_TIG(0:nThreads-1,nSpecies,0:nAlt+1))
    nParticleCell_TIG=0

    !$OMP PARALLEL PRIVATE(iCell, iSpecies, iThread,iParticle)
    
    !$OMP DO SCHEDULE(STATIC)
    PARTICLE_LOOP1: do iParticle=1,nParticle
       iCell=Particles_I(iParticle)%iCell
       if(iCell<0 .or.iCell>nAlt+1) cycle PARTICLE_LOOP1
       iSpecies=Particles_I(iParticle)%iSpecies

       !get the thread number
       iThread=0
       !$ iThread=OMP_get_thread_num()


       
       ! count number of particles on given thread of given species in a given
       ! cell
       nParticleCell_TIG(iThread,iSpecies,iCell)=&
            nParticleCell_TIG(iThread,iSpecies,iCell)+1
       
    enddo PARTICLE_LOOP1
    !$OMP END DO 


    !First thread will always have no offset
    nOffset_TIG(0,:,:) = 0

    !find offset for other threads
    !$OMP DO 
    do iThread=1,nThreads-1
       do iSpecies=1,nSpecies
          do iCell=0,nAlt+1
             nOffset_TIG(iThread,iSpecies,iCell) = &
                  sum(nParticleCell_TIG(0:iThread-1,iSpecies,iCell))
          enddo
       enddo
    enddo
    !$OMP END DO 

    !find maximum number of particles in a given cell
    !$OMP DO 
    do iCell=0,nAlt+1
       do iSpecies=1,nSpecies
          nSortedParticle_II(iSpecies,iCell) = &
               sum(nParticleCell_TIG(:,iSpecies,iCell))
       enddo
    enddo
    !$OMP END DO

    !$OMP END PARALLEL

    MaxSortedParticle=maxval(nSortedParticle_II)
    
    !deallocate previous pointer and reallocate
    if (allocated(SortParticles_III)) then
       deallocate(SortParticles_III)
    endif
    allocate(SortParticles_III(maxval(nSortedParticle_II),nSpecies,0:nAlt+1))


    !$OMP PARALLEL PRIVATE(iCell, iSpecies, iThread,iParticleCell_II,iParticle)

    iParticleCell_II = 1
    !now bucket sort the particles    
    
    !$OMP DO SCHEDULE(STATIC)
    PARTICLE_LOOP2: do iParticle=1,nParticle
       iCell=Particles_I(iParticle)%iCell
       if(iCell<0 .or.iCell>nAlt+1) cycle PARTICLE_LOOP2
       iSpecies=Particles_I(iParticle)%iSpecies
!       if(iSpecies==0)write(*,*)iParticle,nParticle
       
       !get the thread number
       iThread=0
       !$ iThread=OMP_get_thread_num()


       !write(*,*) iThread,nOffset_TIG(iThread,iSpecies,iCell),iParticleCell_II(iSpecies,iCell)
!       if (nOffset_TIG(iThread,iSpecies,iCell)&
!            +iParticleCell_II(iSpecies,iCell) > maxval(nSortedParticle_II))&
!            write(*,*) 'ha', iThread,nOffset_TIG(iThread,iSpecies,iCell),&
!            iParticleCell_II(iSpecies,iCell),maxval(nSortedParticle_II)
       
       !sort the particles
       SortParticles_III(nOffset_TIG(iThread,iSpecies,iCell)&
            +iParticleCell_II(iSpecies,iCell),&
            iSpecies,iCell)%Particle& 
            =>Particles_I(iParticle)

       iParticleCell_II(iSpecies,iCell)=iParticleCell_II(iSpecies,iCell)+1
       
    enddo PARTICLE_LOOP2
    !$OMP END DO
    
    !$OMP END PARALLEL

    deallocate(iParticleCell_II,nOffset_TIG,nParticleCell_TIG)



  end subroutine sort_particles

  !=============================================================================
  ! Generate random permutation based on Knuth's algorithm
  subroutine get_permutation_index_array(nIndices, Index_I)
    integer, intent(in) :: nIndices
    integer, intent(out):: Index_I(nIndices)

    integer :: i, j, k
    integer :: temp
    real    :: RandNum
    !---------------------------------------------------------------------------
    do i=1,nIndices
       Index_I(i)=i
    enddo

    do j=nIndices,2,-1
       RandNum=random_real(iSeed)
       k = floor(j*RandNum) + 1

       ! exchange k and j index
       temp = Index_I(k)
       Index_I(k) = Index_I(j)
       Index_I(j) = temp

    end do

  end subroutine get_permutation_index_array

  !============================================================================
  ! apply coulomb collisions to cell using approach of Takizuka and Abe 1977 
  ! adapted to multiple weight particles by Nanbu and Yonemura 1998
  subroutine apply_coulomb_collision(iCell,iSpeciesIn,jSpeciesIn)
    use ModNumConst, ONLY: cPi,cTwoPi
    use ModConst, ONLY: cEps,cElectronCharge
    integer,intent(in) :: iCell,iSpeciesIn, jSpeciesIn

    real :: Theta,Phi
!    real :: RandNum,RandNum1,RandNum2,RandNum3,RandNum4,RandNum5
    real :: variance

    real,parameter::cCm3ToM3=1e6,cGtoKg=1e-3,cCmToM=1e-2

    real:: uBulkPar,uBulkPerp,Pressure,DensityTotal
    
    real::Velx1,Vely1,Velz1,vmag1,Velx2,Vely2,Velz2,vmag2
    real::vperp1,vpar1,vperp2,vpar2
    real::theta1,phi1,theta2,phi2
    real::DeltaRelVelx,DeltaRelVely,DeltaRelVelz
    real::RelVelx,RelVely,RelVelz,RelVelperp,RelVelMag
    integer :: iParticle,iSpecies1,iSpecies2,iSpecies,jSpecies
    integer ::  iCollider1,iCollider2
    integer :: nCollisions, iCounter1, iCounter2, iCollision
    real :: Density1, Density2,Temp1, Temp2, DensityMin
    real :: Density12,factor
    
    integer :: iParticleSpec,nParticleInCell
    real,allocatable::RandNum_I(:),RandNum2_I(:), &
         RandNum3_I(:),RandNum4_I(:), &
         RandNum5_I(:),RandNum6_I(:)
    integer,allocatable:: IndexPermuted_I(:)
    integer,allocatable:: iCollider1_I(:),iCollider2_I(:)
    real :: CoulLog
    integer :: MaxParticle, MinParticle, nParticleI, nParticleJ

    !variables to hold the weights and things to calculate them
    real :: N11,N12,N21,N22,N11w,N12w,N21w,N22w,P1,P2
    real :: weight_II(2,2), weight1, weight2
    real :: SimParticleRatio, DensityRatio
    real :: Tpar,Tperp,Hpar,Hperp
    logical :: DoTest=.false.

    real :: Delta,sintheta,OneMinusCosTheta
    
    !---------------------------------------------------------------------------
    
    !define the maximum number of particles of either species

    nParticleI=nSortedParticle_II(iSpeciesIn,iCell)
    nParticleJ=nSortedParticle_II(jSpeciesIn,iCell)
    
    MaxParticle=max(nParticleI,nParticleJ)
    MinParticle=min(nParticleI,nParticleJ)
    
    !define collision pairs if like collisions or unlike  
    If (iSpeciesIn==jSpeciesIn) then
       !case of like collisions
       iSpecies1=iSpeciesIn
       iSpecies2=iSpeciesIn
       
       ! Now we need to pair particles randomly. First create index array and 
       ! a permuted index array, These define the collision pairs.
      
       allocate(IndexPermuted_I(MaxParticle))
             
       !define the number of collisions
       if (mod(MaxParticle,2)==0) then 
          ! even number so collisions
          nCollisions=MaxParticle/2
       else
          nCollisions=MaxParticle/2+1
       end if
       
       !get the permuted index
       call get_permutation_index_array(MaxParticle, IndexPermuted_I)
    
       !allocate arrays to hold the colider indicies
       allocate(iCollider1_I(nCollisions),iCollider2_I(nCollisions))
       
       ! Now pair following approach of Takizuka and Abe 1977  where we pair 
       ! in order from the permuted index. Odd number we will just collide the 
       ! first particle again with the remaining particle
       iCounter1=-1
       iCounter2=0
       do iParticle=1,nCollisions
          !Update the counter index
          iCounter1=iCounter1+2
          iCounter2=iCounter2+2

          !define the particle index of the colliders
          iCollider1_I(iParticle)=IndexPermuted_I(iCounter1)

          if(iCounter2<=MaxParticle)then
             iCollider2_I(iParticle)=IndexPermuted_I(iCounter2)
          else
             iCollider2_I(iParticle)=IndexPermuted_I(1)
          endif
       end do
       
    else
       
       ! When particles are of different species. Now we follow approach of Nambu
       ! and randomly permute smaller array
              
       allocate(IndexPermuted_I(MaxParticle))

       !define the number of collisions
       nCollisions=MaxParticle

       !allocate arrays to hold the colider indicies
       allocate(iCollider1_I(nCollisions),iCollider2_I(nCollisions))
       
       !get the permuted index
       call get_permutation_index_array(MinParticle, IndexPermuted_I)

       ! Now pair following approach of Nanbu  where we pair a particle from 
       ! the species 1 with a randomly chosen particle from species 2. 
       ! The approach is to choose the population with a smaller number 
       ! from the from the permuted array. 
       ! Note that we wrap our selection when we reach the end of the permuted  
       ! array. 
       
       !start by selecting species one or two based on which is larger
       if (nParticleI<=nParticleJ) then
          iSpecies1 = jSpeciesIn
          iSpecies2 = iSpeciesIn
       else
          iSpecies1 = iSpeciesIn
          iSpecies2 = jSpeciesIn
       endif
       
       !now assign colider indices
       iCounter2=1
       do iParticle=1,nCollisions
          iCollider1_I(iParticle)=iParticle
          iCollider2_I(iParticle)=IndexPermuted_I(iCounter2)
          
          !Update the counter index
          iCounter2=iCounter2+1
          if (iCounter2>MinParticle) iCounter2=1
       enddo
          
    end If

    !cacluate the Coulomb Logarithm based on formula for mixed ion-ion 
    !collisions the the NRL plasma formulary (p. 34). To do this we first 
    !need momements. Also save the min density which is needed later for the 
    !variance of the scattering angle
    
    !set the coul log
    if (iSpeciesIn==jSpeciesIn) then
       call calc_moments_cell(iSpeciesIn,iCell,&
            Density1,uBulkPar,uBulkPerp,Pressure,Temp1,&
            Tpar,Tperp,Hpar,Hperp)
       Density2=Density1
       Temp2=Density1
       
       DensityMin = Density1
       
       CoulLog  =23.0&
            -log((Mass_I(iSpeciesIn)+Mass_I(jSpeciesIn))&
            /(Mass_I(iSpeciesIn)*Temp1&
            +Mass_I(jSpeciesIn)*Temp2)&
            *sqrt(cCm3ToM3*Density1/Temp1&
            +cCm3ToM3*Density2/Temp2))
       
    else
       call calc_moments_cell(iSpeciesIn,iCell,&
            Density1,uBulkPar,uBulkPerp,Pressure,Temp1,&
            Tpar,Tperp,Hpar,Hperp)
       call calc_moments_cell(jSpeciesIn,iCell,&
            Density2,uBulkPar,uBulkPerp,Pressure,Temp2,&
            Tpar,Tperp,Hpar,Hperp)
       
       DensityMin = min(Density1,Density2)
       
       CoulLog  =23.0&
            -log((Mass_I(iSpeciesIn)+Mass_I(jSpeciesIn))&
            /(Mass_I(iSpeciesIn)*Temp1&
            +Mass_I(jSpeciesIn)*Temp2)&
            *sqrt(cCm3ToM3*Density1/Temp1&
            +cCm3ToM3*Density2/Temp2))
    end if

    if(Temp1<1e-10.or.Temp2<1e-10 &
         .or. Density1<1e-10 .or. Density2<1e-10) then
       !In case of zero temperature or density return with no collisions
       return
    endif
 
    ! Adjust timestep for use in variance calculation. In Nanbu and Yonemura 1998 
    ! this is explained by the use of average time step per real particle.
    ! see eq 12 to get density12.
    Density12=0.0
    do iCollision=1,nCollisions
       iCollider1=iCollider1_I(iCollision)
       iCollider2=iCollider2_I(iCollision)
       weight1 = &
            SortParticles_III(iCollider1,iSpecies1,iCell)%Particle%NumPerParticle
       weight2 = &
            SortParticles_III(iCollider2,iSpecies2,iCell)%Particle%NumPerParticle
       if (weight1==0.0 .or. weight2==0.0) then
          write(*,*)weight1,iCell,iSpecies1,iCollider1
          write(*,*)weight2,iCell,iSpecies2,iCollider2
       endif
       Density12=Density12 + (weight1*weight2)/max(weight1,weight2)
    enddo
    if (iSpecies1==iSpecies2) Density12=2.0*Density12
    factor=(max(Density1,Density2)*Volume_G(iCell))/Density12
    !write(*,*) 'factor',factor
    
    !precompute random numbers for collisions (for optimizing openmp loop)
    allocate(RandNum_I(nCollisions),RandNum2_I(nCollisions), &
         RandNum3_I(nCollisions),RandNum4_I(nCollisions), &
         RandNum5_I(nCollisions),RandNum6_I(nCollisions))
    do iCollision=1,nCollisions
       RandNum_I (iCollision) =random_real(iSeed)
       RandNum2_I(iCollision)=random_real(iSeed)
       RandNum3_I(iCollision)=random_real(iSeed)
       RandNum4_I(iCollision)=random_real(iSeed)
       RandNum5_I(iCollision)=random_real(iSeed)
       RandNum6_I(iCollision)=random_real(iSeed)
    enddo


    !$OMP PARALLEL DO  PRIVATE(Velx2,Vely2,Velz2,Velx1,Vely1,&
    !$OMP Velz1,Delta,sintheta,OneMinusCosTheta,DeltaRelVelx,&
    !$OMP DeltaRelVely,DeltaRelVelz,Theta,variance,Phi,&
    !$OMP RelVelPerp,RelVelMag,RelVelx,RelVely,RelVelz,phi1,phi2,&
    !$OMP theta1,theta2,vmag1,vmag2,vperp1,vperp2,vpar1,vpar2,&
    !$OMP P1,P2, weight1, weight2,iCollider1,iCollider2)
    
    ! Now for each pairing calculate the change in velocity for each particle 
    !using the method of Takizuka and Abe [1977] as modified by 
    !Nambu  for particles not having equal weights.
    do iCollision=1,nCollisions
       !get index of colliding particles
       
       iCollider1=iCollider1_I(iCollision)
       iCollider2=iCollider2_I(iCollision)
    
            
       
       !Based on the weights set the probability factors in the momentum 
       !exchange equations.
       !No use approach described in Sentoku and Kemp, [2008] after eq 20 
       !whereby the lighter particle always collides but heavier particle only 
       !sometimes collides with prob w1/w2 where w2 is heavier particle
       weight1=&
            SortParticles_III(iCollider1,iSpecies1,iCell)%Particle%NumPerParticle
       weight2=&
            SortParticles_III(iCollider2,iSpecies2,iCell)%Particle%NumPerParticle

       !use a rejection scheme when weights are unequal to reject some collisions
       If(weight2>weight1) then
          P1=1.0
          if(RandNum6_I(iCollision)<=(weight1/weight2)) then
             P2=1.0
          else
             P2=0.0
          endif   
       else
          P2=1.0
          if(RandNum6_I(iCollision)<=(weight2/weight1)) then
             P1=1.0
          else
             P1=0.0
          endif
       endif
       
       !save vpar and vperp for each particle
       vpar1=&
            SortParticles_III(iCollider1,iSpecies1,iCell)%Particle%vpar
       vpar2=&
            SortParticles_III(iCollider2,iSpecies2,iCell)%Particle%vpar

       vperp1=&
            SortParticles_III(iCollider1,iSpecies1,iCell)%Particle%vperp
       vperp2=&
            SortParticles_III(iCollider2,iSpecies2,iCell)%Particle%vperp
       
       !calculate the velocity mag for each particle
       vmag1=sqrt(vpar1**2+vperp1**2)
       vmag2=sqrt(vpar2**2+vperp2**2)
       
       !calculate the pitchange (theta) for each particle
       theta1=acos(max(min(vpar1/vmag1,0.999999999999999),-0.999999999999999))
       theta2=acos(max(min(vpar2/vmag2,0.999999999999999),-0.999999999999999))

       !randomly choose phase (azimuthal angle) for each particle
       phi1=RandNum_I(iCollision)*cTwoPi
       phi2=RandNum2_I(iCollision)*cTwoPi
       
       !get the vx,vy, and vz component for each particle z aligned w/ B
       Velx1=vmag1*sin(theta1)*cos(phi1)
       Vely1=vmag1*sin(theta1)*sin(phi1)
       Velz1=vmag1*cos(theta1)
       
       
       Velx2=vmag2*sin(theta2)*cos(phi2)
       Vely2=vmag2*sin(theta2)*sin(phi2)
       Velz2=vmag2*cos(theta2)
       
       !write(*,*) 'at start2'
       !write(*,*) (Velx1**2+Vely1**2+Velz1**2)+(Velx2**2+Vely2**2+Velz2**2)

       !what do we know about the particles
       if(DoTest) then
          write(*,*)'iCollider1:',iCollider1
          write(*,*)vpar1,vperp1,Velz1,sqrt(Velx1**2+Vely1**2)
          
          write(*,*)'iCollider2:',iCollider2
          write(*,*)vpar2,vperp2,Velz2,sqrt(Velx2**2+Vely2**2)
       endif

       !get the relative velocity in the laboratory frame
       RelVelx=Velx1-Velx2
       RelVely=Vely1-Vely2
       RelVelz=Velz1-Velz2
       
       RelVelMag=sqrt(RelVelx**2.0+RelVely**2.0+RelVelz**2.0)
       RelVelPerp=sqrt(RelVelx**2.0+RelVely**2.0)
       
       !choose the post collision phase in the frame of relative velocity
       Phi=RandNum3_I(iCollision)*cTwoPi
       
       !choose the scattering angle in the post collision (Theta). 
       !This is done by 
       !setting a variation <d^2>=e1^2*e2^2*nL*lambda*dt&
       !/(8*pi*epsilon0^2*m12^2*u^3)
       !where nL is the lower density, u is relative velocity, 
       !m12 is the reduced mass (m1*m2)/(m1+m2), and lamda is the 
       !coulomb logarithm. Once d is chosen Theta is recovered 
       !by d=tan(Theta/2)
       
       ! set variance (note we only assume singly charged particles
       !also note we use SI for this one and only formula instead of 
       !cgs). Note Nambu suggests a change of the variance so as to  
       !correct for different particle weights and pairing 
       !statistics.
       variance=(cElectronCharge**4*cCm3ToM3*DensityMin*CoulLog)&
            /(8.*cPi*cEps**2*(cGtoKg*ReducedMass_II(iSpecies1,iSpecies2))**2&
            *(RelVelMag*cCmToM)**3)*DtCollide*factor

       !sample theta from the random distribution
       if (variance<1.0) then
          !Delta = sqrt(-2.0*variance*log(RandNum4_I(iCollision)))&
          !     *cos(cTwoPi*RandNum5_I(iCollision))
          Delta=sqrt(abs(variance*log(1.0-RandNum4_I(iCollision))))&
               *cos(cTwoPi*RandNum5_I(iCollision))
          !Theta=2.0*atan(sqrt(-2.0*variance*log(RandNum4_I(iCollision)))&
          !     *cos(cTwoPi*RandNum5_I(iCollision)))
          
          Theta=2.0*atan(Delta)
          ! eq 7a and 7b from TA77
          sintheta=max(min(2.0*Delta/(1.0+Delta**2.0),1.0),-1.0)
          OneMinusCosTheta=max(min(2.0*Delta**2.0/(1.0+Delta**2),2.0),0.0)
       else
          Theta=cPi*RandNum4_I(iCollision)
          sintheta=sin(Theta)
          OneMinusCosTheta=1.0-cos(Theta)
          Delta=-1.0
       endif

       ! need to check if it should be 2.0 or 1.o in front
       !Theta=1.0*atan(sqrt(-2.0*variance*log(RandNum4_I(iCollision)))&
       !     *cos(cTwoPi*RandNum5_I(iCollision)))
       
!       if(DoTest) then
!          write(*,*) 'CoulLog_II(iSpecies1,iSpecies2)',CoulLog_II(iSpecies1,iSpecies2)
!          write(*,*) 'cElectronCharge',cElectronCharge
!          write(*,*) 'cCm3ToM3*DensityTotal',cCm3ToM3*DensityTotal
!          write(*,*) 'cEps',cEps
!          write(*,*) 'cGtoKg*Mass_I(iSpecies1)',cGtoKg*Mass_I(iSpecies1)
!          write(*,*) 'cGtoKg*Mass_I(iSpecies2)',cGtoKg*Mass_I(iSpecies2)
!          write(*,*) 'cGtoKg*ReducedMass_II(iSpecies1,iSpecies2)',cGtoKg*ReducedMass_II(iSpecies1,iSpecies2)
!          write(*,*) 'DtMove',DtMove
!          write(*,*) 'RelVelMag*cCmToM',RelVelMag*cCmToM
          !write(*,*)'theta,variance',theta*180./3.14,variance
!          stop
!       endif
       
       !now that we know the post collision scatter and phase angle compute 
       !the velocity change for each component (eq. 4a-4d' of Takizuka and Abe)
       if(RelVelperp<1e-20)then
          !DeltaRelVelx=RelVelMag*sin(Theta)*cos(Phi)
          DeltaRelVelx=RelVelMag*sintheta*cos(Phi)
          !DeltaRelVely=RelVelMag*sin(Theta)*sin(Phi)
          DeltaRelVely=RelVelMag*sintheta*sin(Phi)
          !DeltaRelVelz=-RelVelMag*(1.0-cos(Theta))
          DeltaRelVelz=-RelVelMag*OneMinusCosTheta
       else
          !DeltaRelVelx=(RelVelx/RelVelperp)*RelVelz*sin(Theta)*cos(Phi)&
          !     -(RelVely/RelVelperp)*RelVelmag*sin(Theta)*sin(Phi)&
          !     -RelVelx*(1.0-cos(Theta))
          DeltaRelVelx=(RelVelx/RelVelperp)*RelVelz*sintheta*cos(Phi)&
               -(RelVely/RelVelperp)*RelVelmag*sintheta*sin(Phi)&
               -RelVelx*OneMinusCosTheta
          !DeltaRelVely=(RelVely/RelVelperp)*RelVelz*sin(Theta)*cos(Phi)&
          !     +(RelVelx/RelVelperp)*RelVelmag*sin(Theta)*sin(Phi)&
          !     -RelVely*(1.0-cos(Theta))
          DeltaRelVely=(RelVely/RelVelperp)*RelVelz*sintheta*cos(Phi)&
               +(RelVelx/RelVelperp)*RelVelmag*sintheta*sin(Phi)&
               -RelVely*OneMinusCosTheta
          !DeltaRelVelz=-RelVelperp*sin(Theta)*cos(Phi)&
          !     -RelVelz*(1.0-cos(Theta))
          DeltaRelVelz=-RelVelperp*sintheta*cos(Phi)&
               -RelVelz*OneMinusCosTheta
       endif

       
       !we can now update the post collision velocity of each particle
       Velx1=Velx1&
            +ReducedMass_II(iSpecies1,iSpecies2)/Mass_I(iSpecies1)&
            *P1*DeltaRelVelx
       Vely1=Vely1&
            +ReducedMass_II(iSpecies1,iSpecies2)/Mass_I(iSpecies1)&
            *P1*DeltaRelVely
       Velz1=Velz1&
            +ReducedMass_II(iSpecies1,iSpecies2)/Mass_I(iSpecies1)&
            *P1*DeltaRelVelz

       Velx2=Velx2&
            -ReducedMass_II(iSpecies1,iSpecies2)/Mass_I(iSpecies2)&
            *P2*DeltaRelVelx
       Vely2=Vely2&
            -ReducedMass_II(iSpecies1,iSpecies2)/Mass_I(iSpecies2)&
            *P2*DeltaRelVely
       Velz2=Velz2&
            -ReducedMass_II(iSpecies1,iSpecies2)/Mass_I(iSpecies2)&
            *P2*DeltaRelVelz

       

 
       !These updated velocities can now be converted back to vpar and vperp 
       !and stored back in their respective partiles
       SortParticles_III(iCollider1,iSpecies1,iCell)%Particle%vpar=&
            Velz1
       SortParticles_III(iCollider2,iSpecies2,iCell)%Particle%vpar=&
            Velz2

       SortParticles_III(iCollider1,iSpecies1,iCell)%Particle%vperp=&
            sqrt(Velx1**2+Vely1**2)
       SortParticles_III(iCollider2,iSpecies2,iCell)%Particle%vperp=&
            sqrt(Velx2**2+Vely2**2)


  !     write(*,*) 'at end',P1,P2
  !     write(*,*) (Velx1**2+Vely1**2+velz1**2)&
  !          +(Velx2**2+Vely2**2+Velz2**2)
       !stop
    enddo
    
    !$OMP END PARALLEL DO
    
    !deallocate
    deallocate(IndexPermuted_I)
    deallocate(RandNum_I,RandNum2_I, &
         RandNum3_I,RandNum4_I, &
         RandNum5_I,RandNum6_I)

    !deallocate
    deallocate(iCollider1_I,iCollider2_I)
  end subroutine apply_coulomb_collision
  
  !=============================================================================
  ! 
  subroutine apply_wave_particle_interaction(rWaveRef)
    use ModPlanetConst, ONLY: Planet_,rPlanet_I
    use ModNumConst, ONLY: cTwoPi
    real, intent(in):: rWaveRef
    real :: rPlanetCM, rCoord, Dperp,variance
    integer :: iSpecies, iParticle
    real :: vpar, vperp,vmag,dVx,dVy,theta,phi,Velx,Vely
    !random numbers for sampling
    real,allocatable :: RandNum1_I(:),RandNum2_I(:),RandNum3_I(:),&
         RandNum4_I(:),RandNum5_I(:)
    real, parameter :: cMtoCm=1e2
    !--------------------------------------------------------------------------

    !precompute random numbers for collisions (for optimizing openmp loop)
    allocate(RandNum1_I(nParticle), RandNum2_I(nParticle), &
         RandNum3_I(nParticle),RandNum4_I(nParticle),RandNum5_I(nParticle))

    do iParticle=1,nParticle
       RandNum1_I(iParticle)=random_real(iSeed)
       RandNum2_I(iParticle)=random_real(iSeed)
       RandNum3_I(iParticle)=random_real(iSeed)
       RandNum4_I(iParticle)=random_real(iSeed)
       RandNum5_I(iParticle)=random_real(iSeed)
    enddo    

    rPlanetCM=rPlanet_I(Planet_)*cMtoCm
    !$OMP PARALLEL PRIVATE(iParticle, iSpecies,vpar,vperp,vmag,rCoord,Dperp,&
    !$OMP variance,dVx,dVy,Velx,Vely,theta,phi)
    
    !$OMP DO  
    do iParticle=1,nParticle
       iSpecies=Particles_I(iParticle)%iSpecies
       vpar =Particles_I(iParticle)%vpar
       vperp=Particles_I(iParticle)%vperp
       vmag =sqrt(vpar**2+vperp**2)

       rCoord=rPlanetCM+Particles_I(iParticle)%Alt
       Dperp=Dperp_I(iSpecies)*(rCoord/rWaveRef)**Dexp_I(iSpecies)
       variance=2*Dperp*DtMove
       dVx=sqrt(-2.0*variance*log(RandNum1_I(iParticle)))&
            *cos(cTwoPi*RandNum2_I(iParticle))
       dVy=sqrt(-2.0*variance*log(RandNum3_I(iParticle)))&
            *cos(cTwoPi*RandNum4_I(iParticle))

       !calculate the pitchange (theta) for each particle
       theta=acos(vpar/vmag)
       !randomly choose phase (azimuthal angle) for particle
       phi=RandNum5_I(iParticle)*cTwoPi
       
       Velx=vmag*sin(theta)*cos(phi)
       Vely=vmag*sin(theta)*sin(phi)

       Particles_I(iParticle)%vperp=sqrt((Velx+dVx)**2+(Vely+dVy)**2)
    end do
    !$OMP END DO 
    
    !$OMP END PARALLEL
  end subroutine apply_wave_particle_interaction

  !============================================================================
  ! loop over all cells in the computational domain and check if number of 
  ! particles in the cell are equal to the target number of particles per cell 
  ! if outside the target then split or join particles as required to reach 
  ! target number of particles in a cell. Then resort when done
  subroutine split_join_particles(Tolerance)
    !integer, intent(in) :: iSpecies
    real   , intent(in) :: Tolerance
    integer :: iCell, nParticleInCell
    integer,parameter :: LowerLimit=10
    integer :: iSpecies
    !variables for splitting joining and holding new paritcles
    integer :: nJoin, nSplit, nNew, nNewOld
    type(particle),allocatable :: NewParticles_I(:), NewParticlesOld_I(:)
    type(particle),allocatable :: NewParticlesTmp_I(:)
    logical :: IsJoinSuccess
    integer :: nJoinRevise

    !variables for assignment of particles
    type(particle),allocatable ::ParticlesOld_I(:)
    integer :: nParticleOld,nAvail, iParticle,iParticleOld
    integer,allocatable :: IndexAvail_I(:)
    !--------------------------------------------------------------------------
    
    nNew=0
    nNewOld=0
    
    species_loop: do iSpecies=1,nSpecies
       !loop over cells
       cell_loop: do iCell=1,nAlt
          ! get the number of particles in the cell of a given species
          nParticleInCell=nSortedParticle_II(iSpecies,iCell)
          
          ! if particle in cell are less than some lower limit cycle
          if (nParticleInCell <= LowerLimit) cycle cell_loop
          
          !Compare to tolerance and if inside then cycle
          if(nParticleInCell>=nParticlePerCell_I(iSpecies)*(1.0-Tolerance) &
               .and. nParticleInCell<=nParticlePerCell_I(iSpecies)&
               *(1.0+Tolerance)) then
             cycle cell_loop
          endif
          
          !Save the number of new particles
          nNewOld=nNew
          
          ! on first time through loop nNewOld is zero so special case
          if (nNewOld>0) then
             allocate(NewParticlesOld_I(nNewOld))
             NewParticlesOld_I=NewParticles_I
          endif
          
          !Check if we need to split or join
          if (nParticleInCell>nParticlePerCell_I(iSpecies)) then
             !set number of joined particles to create never letting the number 
             !exceed 10% of the particles at a time
             nJoin=min(abs(nParticleInCell-nParticlePerCell_I(iSpecies)),&
                  floor(0.1*nParticleInCell))
             !nJoin=abs(nParticleInCell-nParticlePercell)
             
             !need loop to deal with case when not all particles can be joined
             IsJoinSuccess=.false.
             join: do while (.not.IsJoinSuccess)
                if(.not.allocated(NewParticlesTmp_I))&
                     allocate(NewParticlesTmp_I(nJoin))
                call join_particles_cell(iCell,iSpecies,nJoin,NewParticlesTmp_I,&
                     nJoinRevise)
                ! check if nJoin particles could be combined while keeping 
                ! reasonable descretization in velocity space
                if (nJoin == nJoinRevise) then
                   IsJoinSuccess=.true.
                else
                   deallocate(NewParticlesTmp_I)
                   nJoin=nJoinRevise
                endif
             end do join
             nNew=nNewOld+nJoin
          else
             !set number of particles to split never letting the number 
             !exceed 10% of the particles at a time
             nSplit=min(abs(nParticleInCell-nParticlePerCell_I(iSpecies)),&
                  floor(0.4*nParticleInCell))
             !nSplit=abs(nParticleInCell-nParticlePercell)
             nNew=nNewOld+2*nSplit
             if(.not.allocated(NewParticlesTmp_I))&
                  allocate(NewParticlesTmp_I(2*nSplit))
             call split_particles_cell(iCell,iSpecies,nSplit,NewParticlesTmp_I)
          endif
          
          ! append new particles to new particle list
          
          if (nNewOld>0) then
             deallocate(NewParticles_I)
             allocate(NewParticles_I(nNew))
             NewParticles_I(1:nNewOld)=NewParticlesOld_I
             NewParticles_I(nNewOld+1:nNew)=NewParticlesTmp_I
             deallocate(NewParticlesOld_I)
          else
             if(allocated(NewParticles_I)) deallocate(NewParticles_I)
             allocate(NewParticles_I(nNew))
             NewParticles_I=NewParticlesTmp_I
          endif
          
          
          !deallocate these arrays as they need to be allocated the next time
          !around the loop
          deallocate(NewParticlesTmp_I)
          
          ! now add new particles to the list, note that the old particles 
          ! involved in spliting or joining are already removed from the main 
          ! particle list
          
       end do cell_loop
    end do species_loop
    ! rebuild master particle list including new particles
    if (nNew >0) then

       !Calculate the
       !number of open spots (nAvail) then save Particles_I and allocate a new 
       !Particles_I array to save the good particles and the newly created ones
       allocate(IndexAvail_I(nParticle))
       where(Particles_I%IsOpen)
          IndexAvail_I=1
       elsewhere
          IndexAvail_I=0
       end where
       nAvail=sum(IndexAvail_I)
       deallocate(IndexAvail_I)
       
       !save old particle array information
       nParticleOld=nParticle
       allocate(ParticlesOld_I(nParticleOld))
       ParticlesOld_I=Particles_I
       deallocate(Particles_I)
       
       !allocate new Particles_I array
       nParticle=nParticleOld-nAvail+nNew
       allocate(Particles_I(nParticle))
       
       !now fill new Particle_I array with old particles that are not open and 
       !the newly created particles
       iParticle=1
       do iParticleOld=1,nParticleOld
          if(.not.ParticlesOld_I(iParticleOld)%IsOpen) then
             Particles_I(iParticle)=ParticlesOld_I(iParticleOld)
             iParticle=iParticle+1
          endif
       end do
       Particles_I(nParticleOld-nAvail+1:nParticle)=NewParticles_I
       
       deallocate(ParticlesOld_I)
       
    end if
    !deallocate to save memory
    if(allocated(NewParticles_I))deallocate(NewParticles_I)

    !clean up any particles split too far so weight is zero
    call clean_zero_weight_particles

  end subroutine split_join_particles
  !============================================================================
  ! Split nSplit particles of species iSpecies in cell iCell
  subroutine split_particles_cell(iCell,iSpecies,nSplit,NewParticle_I)
    use ModSort, ONLY: sort_quick
    integer,intent(in) :: iCell,iSpecies,nSplit
    type(particle),intent(out):: NewParticle_I(2*nSplit)
    integer :: nParticleInCell,iParticle,iParticleOld
    real,allocatable :: weights_I(:)
    integer,allocatable:: IndexSort_I(:)
    integer :: iSplit1, iSplit2

    !variables for assignment of particles
    type(particle),allocatable ::ParticlesOld_I(:)
    integer :: nParticleOld,nAvail,i,j
    integer,allocatable :: IndexAvail_I(:)

    !---------------------------------------------------------------------------
    nParticleInCell=nSortedParticle_II(iSpecies,iCell)
    allocate(weights_I(nParticleInCell))
    allocate(IndexSort_I(nParticleInCell))
    
    !get array of weights
    do iParticle=1,nParticleInCell
       weights_I(iParticle)=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
    enddo
    
    !get index array that sorts weights_I from lowest to biggest
    call sort_quick(nParticleInCell,weights_I,IndexSort_I)
    
    !now split the nSplit heaviest particles in the cell
    iSplit1=-1
    iSplit2=0
    do iParticle=nParticleInCell,nParticleInCell-nSplit+1,-1
       iSplit1=iSplit1+2
       iSplit2=iSplit2+2

       !copy old particle to new particle with half the weight and shift 
       !x+-deltaX/nParticleInCell in position
       !NewParticle_I(iSplit1)=&
       !     SortParticles_III(iParticle,iSpecies,iCell)%Particle
       !NewParticle_I(iSplit2)=&
       !     SortParticles_III(iParticle,iSpecies,iCell)%Particle
       
       NewParticle_I(iSplit1)%iSpecies=iSpecies
       NewParticle_I(iSplit2)%iSpecies=iSpecies

       NewParticle_I(iSplit1)%iCell=iCell
       NewParticle_I(iSplit2)%iCell=iCell            

       NewParticle_I(iSplit1)%NumPerParticle=&
            0.5*&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
       NewParticle_I(iSplit2)%NumPerParticle=&
            0.5*&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle

       NewParticle_I(iSplit1)%vpar=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%vpar
       NewParticle_I(iSplit2)%vpar=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%vpar

       NewParticle_I(iSplit1)%vperp=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%vperp
       NewParticle_I(iSplit2)%vperp=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%vperp


       NewParticle_I(iSplit1)%Alt=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%Alt&
            +dAlt_G(iCell)/real(nParticleInCell)
            !min(SortParticles_III(iParticle,iSpecies,iCell)%Particle%Alt&
            !+dAlt_G(iCell)/real(nParticleInCell)&
            !,AltTop_F(iCell))
       NewParticle_I(iSplit2)%Alt=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%Alt&
            -dAlt_G(iCell)/real(nParticleInCell)
            !max(SortParticles_III(iParticle,iSpecies,iCell)%Particle%Alt&
            !-dAlt_G(iCell)/real(nParticleInCell)&
            !,AltBot_F(iCell))

       !kludge set vperp such that first invarient is conserved
!       NewParticle_I(iSplit1)%vperp=&
!            SortParticles_III(iParticle,iSpecies,iCell)%Particle%vperp&
!            *(SortParticles_III(iParticle,iSpecies,iCell)%Particle%Alt&
!            /NewParticle_I(iSplit1)%Alt)**1.5
!
!       NewParticle_I(iSplit2)%vperp=&
!            SortParticles_III(iParticle,iSpecies,iCell)%Particle%vperp&
!            *(SortParticles_III(iParticle,iSpecies,iCell)%Particle%Alt&
!            /NewParticle_I(iSplit2)%Alt)**1.5
       
       
       !check if new particle is in computational domain
       if (NewParticle_I(iSplit1)%Alt<AltBot_F(1) .or. &
            NewParticle_I(iSplit1)%Alt>AltTop_F(nAlt)) then
          NewParticle_I(iSplit1)%IsOpen=.true.
       else
          NewParticle_I(iSplit1)%IsOpen=.false.
       endif

       if (NewParticle_I(iSplit2)%Alt<AltBot_F(1) .or. &
            NewParticle_I(iSplit2)%Alt>AltTop_F(nAlt)) then
          NewParticle_I(iSplit2)%IsOpen=.true.
       else
          NewParticle_I(iSplit2)%IsOpen=.false.
       endif

!          NewParticle_I(iSplit1)%IsOpen=.false.
!          NewParticle_I(iSplit2)%IsOpen=.false.


       !set particle to be split to open
       SortParticles_III(iParticle,iSpecies,iCell)%Particle%IsOpen=.true.
    enddo

!    !recreate particle array inserting new particles and removing open ones
!    !Calculate the
!    !number of open spots (nAvail) then save Particles_I and allocate a new 
!    !Particles_I array to save the good particles and the newly created ones
!    allocate(IndexAvail_I(nParticle))
!    where(Particles_I%IsOpen)
!       IndexAvail_I=1
!    elsewhere
!       IndexAvail_I=0
!    end where
!    nAvail=sum(IndexAvail_I)
!    deallocate(IndexAvail_I)
!    
!    !save old particle array information
!    nParticleOld=nParticle
!    allocate(ParticlesOld_I(nParticleOld))
!    ParticlesOld_I=Particles_I
!    deallocate(Particles_I)
!    
!    !allocate new Particles_I array
!    nParticle=nParticleOld-nAvail+2*nSplit
!    allocate(Particles_I(nParticle))
!    
!    !now fill new Particle_I array with old particles that are not open and 
!    !the newly created particles
!    iParticle=1
!    do iParticleOld=1,nParticleOld
!       if(.not.ParticlesOld_I(iParticleOld)%IsOpen) then
!          Particles_I(iParticle)=ParticlesOld_I(iParticleOld)
!          iParticle=iParticle+1
!       endif
!    end do
!    Particles_I(nParticleOld-nAvail+1:nParticle)=NewParticle_I
!    
!    deallocate(ParticlesOld_I)
!    
    
    !deallocate to save memory
    deallocate(weights_I, IndexSort_I)
  end subroutine split_particles_cell

  !============================================================================
  ! Join nJoin particles. Idea is to sort all particles in a cell into bins in 
  ! phase space. Also sort particles by weight. Then select smallest weight 
  ! particle and join it to another particle in the bin. Lapenta 2002 details the 
  ! idea behind this although the sorting into bins is not described there but 
  ! used here to ensure that all joined particles are close in phase space. 
  subroutine join_particles_cell(iCell,iSpecies,nJoin,NewParticle_I,nJoinRevise)
    use ModSort, ONLY: sort_quick
    use ModNumConst, ONLY: cPi,cTwoPi
    
    integer,intent(in) :: iSpecies,iCell,nJoin
    type(particle),intent(out):: NewParticle_I(nJoin)
    integer,intent(out):: nJoinRevise
    integer :: nParticleInCell,iJoin
    logical :: FoundParticle
    
    real :: density,uBulkPar,uBulkPerp,Pressure,Temp, uTherm
    real :: dVel, vParMin, vPerpMin
    real :: Tpar,Tperp,Hpar,Hperp
    integer, parameter :: nVel = 100
    real :: vPar_C(nVel), vPerp_C(nVel)
    integer :: iVel, iParticle, iVpar,iVperp
    integer, allocatable   :: IndexBinParticle_III(:,:,:)
    integer, allocatable   :: IndexVpar_I(:),IndexVperp_I(:)
    integer, allocatable   :: nParticleBin_II(:,:)
    integer, allocatable   :: nAvailBin_II(:,:)   
    logical, allocatable   :: IsAvail_I(:)
    integer, allocatable   :: jParticle_II(:,:)
    integer :: jParticle
    
    !grid parameters
    integer, parameter :: nDim =2, Vperp_=1,Vpar_=2,nVar=1, PSD_=1
    
    real :: TrueParticles
    real,allocatable :: weights_I(:)
    integer,allocatable:: IndexSort_I(:)

    integer,save :: iCounter=0
    
    !tmp variables for joining
    real :: weight1,weight2, Alt1, Alt2, Vpar1, Vpar2, Vperp1, Vperp2
    !variables for assignment of particles
    type(particle),allocatable ::ParticlesOld_I(:)
    integer :: nParticleOld,nAvail,i,j,iParticleOld
    integer,allocatable :: IndexAvail_I(:)
    logical :: IsEnough
    !---------------------------------------------------------------------------
    TrueParticles=0
    nJoinRevise=nJoin
    
    nParticleInCell=nSortedParticle_II(iSpecies,iCell)
    if(nParticleInCell==0)return
    allocate(IndexBinParticle_III(nVel,nVel,nParticleInCell))
    allocate(nParticleBin_II(nVel,nVel))
    allocate(jParticle_II(nVel,nVel))
    allocate(nAvailBin_II(nVel,nVel))
    allocate(IsAvail_I(nParticleInCell)) 
    allocate(IndexVpar_I(nParticleInCell),IndexVperp_I(nParticleInCell))
    IsAvail_I=.true.

    ! first get index array that sorts particles by weight from lowest to highest 
    ! so we can join lowest weight particles
    allocate(weights_I(nParticleInCell))
    allocate(IndexSort_I(nParticleInCell))
    !get array of weights
    do iParticle=1,nParticleInCell
       weights_I(iParticle)=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
    enddo
    
    !get index array that sorts weights_I from lowest to biggest
    call sort_quick(nParticleInCell,weights_I,IndexSort_I)

    !Now we discretize velocity space grid based on thermal velocity

    ! get moments in cell so we can calculate thermal velocity
    call calc_moments_cell(iSpecies,iCell,&
         density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp,Hpar,Hperp)
    
    !write(*,*) 'density,uBulkPar,uBulkPerp,Pressure,Temp'&
    !     ,density,uBulkPar,uBulkPerp,Pressure,Temp
    ! calculate thermal velocity based on max of tpar or tperp
    uTherm=max(sqrt(8.0*cBoltzmannCGS*Tpar/Mass_I(iSpecies)/cPi),&
         sqrt(8.0*cBoltzmannCGS*Tperp/Mass_I(iSpecies)/cPi))

    ! discretize velocity space, center around bulk velocity to 
    !5 uTherm in every direction with grid size .2 uTherm
    vParMin=uBulkPar-5.0*uTherm
    vPerpMin=0.0
    dVel=0.1*uTherm
    
    ! discretize velocity space but make sure enough particle are in enough 
    ! bins to join particles
    IsEnough=.false.
    discretize: do while (.not.IsEnough)
       do iVel=1,nVel
          vPar_C(iVel)=vParMin+iVel*dVel
          vPerp_C(iVel)=vPerpMin+iVel*dVel
       enddo
       
       
       !Sort particles into bins
       jParticle_II(:,:)=1
       nParticleBin_II(:,:)=0
       do iParticle=1,nParticleInCell
          iVpar =floor((SortParticles_III(iParticle,iSpecies,iCell)%Particle%vpar &
               -vParMin )/dVel)
          iVperp=floor((SortParticles_III(iParticle,iSpecies,iCell)%Particle%vperp&
               -vPerpMin)/dVel)
          
          !save the bin indices for each particle
          IndexVpar_I(iParticle)=iVpar
          IndexVperp_I(iParticle)=iVperp
          
          !cycle if either index is outside range
          if (iVpar<1 .or. iVpar>nVel .or.iVperp<1 .or. iVperp>nVel) then
             cycle
          endif
          
          !append the particle index for the bin
          jParticle=jParticle_II(iVpar,iVperp)
          IndexBinParticle_III(iVpar,iVperp,jParticle)=iParticle
          nParticleBin_II(iVpar,iVperp)=nParticleBin_II(iVpar,iVperp)+1
          jParticle_II(iVpar,iVperp)=jParticle_II(iVpar,iVperp)+1
          
       enddo
       
       !Initially all particle are available for joining
       nAvailBin_II = nParticleBin_II
       
       !check if there are enough particles in the bins for joining
       iCounter=0
       do iVpar=1,nVel
          do iVperp=1,nVel
             iCounter = iCounter+nAvailBin_II(iVpar,iVperp)/2
          enddo
       enddo
       
       !write(*,*) iCounter,nJoin,dVel/uTherm
       if (iCounter>= nJoin)then
          !the discretization is coarse enough! 
          IsEnough=.true.
       else
          !not coarse enough so coarsen by 2
          dVel=2.*dVel
          !check how coarse
          If (dVel>0.1*uTherm) then
             !write(*,*) 'PW ERROR: iCell,iSpecies',iCell,iSpecies
             !call con_stop('PW ERROR: To coarse in Join')
             !if too coarse then revise the number to join and try again
             nJoinRevise=iCounter
             return
          endif
       endif
    end do discretize
    ! now start joining particle starting from lowest weight. Proceedure is 
    !to select lowest weight particle, check how many are available for joining 
    !in its velocity space bin and then join it to the first available particle 
    !in the bin. Update the availablility list and move on to the next particle
    iCounter=0
    JOIN_LOOP: do iJoin=1,nJoin
       FoundParticle=.false.
       select_particle: do while (.not.FoundParticle)
          iCounter=iCounter+1
          if (iCounter>nParticleInCell) then
             !failure case! 
             write(*,*) 'not enough particles close enough to join'
             write(*,*) iJoin,nJoin,maxval(nAvailBin_II)
          endif
          iParticle=IndexSort_I(iCounter)
                    
          !check if particle is available for joining and another particle is 
          !available in the same bin
          if (IsAvail_I(iParticle)) then
             iVpar =IndexVpar_I(iParticle)
             iVperp=IndexVperp_I(iParticle)
             if (iVpar<1 .or. iVpar>nVel .or.iVperp<1 .or. iVperp>nVel) then
                cycle select_particle
             endif
             if (nAvailBin_II(iVpar,iVperp)>=2) FoundParticle=.true.
          endif
       end do select_particle
       !update availablility 
       IsAvail_I(iParticle)=.false.
    
       FoundParticle=.false.
       i=1
       select_particle2: do while (.not.FoundParticle)
          jParticle=IndexBinParticle_III(iVpar,iVperp,i)
          !check if available
          if (IsAvail_I(jParticle)) FoundParticle=.true.
          i=i+1
       end do select_particle2
       !update availability 
       IsAvail_I(jParticle)=.false.
       nAvailBin_II(iVpar,iVperp)=nAvailBin_II(iVpar,iVperp)-2
   
       !Now we can join iParticle and jParticle
       !first set these to open so they can be removed from the particle list
       SortParticles_III(iParticle,iSpecies,iCell)%Particle%IsOpen=.true.
       SortParticles_III(jParticle,iSpecies,iCell)%Particle%IsOpen=.true.
       
       !set some tmp variables so lines are not so long
       weight1=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
       weight2=&
            SortParticles_III(jParticle,iSpecies,iCell)%Particle%NumPerParticle
       Alt1=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%Alt
       Alt2=&
            SortParticles_III(jParticle,iSpecies,iCell)%Particle%Alt
       Vpar1=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%Vpar
       Vpar2=&
            SortParticles_III(jParticle,iSpecies,iCell)%Particle%Vpar
       Vperp1=&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%Vperp
       Vperp2=&
            SortParticles_III(jParticle,iSpecies,iCell)%Particle%Vperp
       
       !create the newly joined particle
       NewParticle_I(iJoin)%iSpecies=iSpecies
       NewParticle_I(iJoin)%iCell=iCell
       NewParticle_I(iJoin)%vpar=&
            (weight1*Vpar1+weight2*Vpar2)/(weight1+weight2)
       NewParticle_I(iJoin)%vperp=&
            (weight1*Vperp1+weight2*Vperp2)/(weight1+weight2)
       NewParticle_I(iJoin)%Alt=(weight1*Alt1+weight2*Alt2)/(weight1+weight2)
       NewParticle_I(iJoin)%IsOpen=.false.
       NewParticle_I(iJoin)%NumPerParticle=weight1+weight2
    end do JOIN_LOOP
       
    
!    !recreate particle array inserting new particles and removing open ones
!    !Calculate the
!    !number of open spots (nAvail) then save Particles_I and allocate a new 
!    !Particles_I array to save the good particles and the newly created ones
!    allocate(IndexAvail_I(nParticle))
!    where(Particles_I%IsOpen)
!       IndexAvail_I=1
!    elsewhere
!       IndexAvail_I=0
!    end where
!    nAvail=sum(IndexAvail_I)
!    deallocate(IndexAvail_I)
!    
!    !save old particle array information
!    nParticleOld=nParticle
!    allocate(ParticlesOld_I(nParticleOld))
!    ParticlesOld_I=Particles_I
!    deallocate(Particles_I)
!    
!    !allocate new Particles_I array
!    nParticle=nParticleOld-nAvail+nJoin
!    allocate(Particles_I(nParticle))
!    
!    !now fill new Particle_I array with old particles that are not open and 
!    !the newly created particles
!    iParticle=1
!    do iParticleOld=1,nParticleOld
!       if(.not.ParticlesOld_I(iParticleOld)%IsOpen) then
!          Particles_I(iParticle)=ParticlesOld_I(iParticleOld)
!          iParticle=iParticle+1
!       endif
!    end do
!    Particles_I(nParticleOld-nAvail+1:nParticle)=NewParticle_I
!    
!    deallocate(ParticlesOld_I)
!    deallocate(NewParticle_I)

    deallocate(IndexBinParticle_III)
    deallocate(nParticleBin_II)
    deallocate(jParticle_II)
    deallocate(nAvailBin_II)
    deallocate(IsAvail_I)
    deallocate(weights_I, IndexSort_I)
    deallocate(IndexVpar_I,IndexVperp_I)
  end subroutine join_particles_cell

  !============================================================================
  ! Bury particles on a line. Useful when considering multiple lines on a proc
  subroutine bury_line(iLine)
    integer, intent(in) :: iLine
    logical,save :: IsFirstCall=.true.
    !----------------------------------------------------------------------------

    !on first call allocate array to bury particles
    if (IsFirstCall) then
       allocate(BuriedParticles_I(nLine))
       IsFirstCall=.false.
    endif
    
    !deallocate currently buried line in the array
    if (allocated(BuriedParticles_I(iLine)%SavedParticles_I)) &
         deallocate(BuriedParticles_I(iLine)%SavedParticles_I)
    
    !reallocate with the current number of particles
    allocate(BuriedParticles_I(iLine)%SavedParticles_I(nParticle))

    !now save particles and number on line for later use
    BuriedParticles_I(iLine)%SavedParticles_I=Particles_I
    BuriedParticles_I(iLine)%nParticleOnLine=nParticle
    
    !deallocate the Particles_I array now that those particles are buried 
    deallocate(Particles_I)
  end subroutine bury_line

  !============================================================================
  ! Disinter particles on a line.Useful when considering multiple lines on a proc
  subroutine disinter_line(iLine)
    integer, intent(in) :: iLine
    !---------------------------------------------------------------------------

    !set total number of particles to buried value
    nParticle = BuriedParticles_I(iLine)%nParticleOnLine
    
    !deallocate current particle array and reallocate with number of particle on 
    !the line
    if (allocated(Particles_I)) &
         deallocate(Particles_I)

    allocate(Particles_I(nParticle))

    !now load particles and number on line for later use
    Particles_I=BuriedParticles_I(iLine)%SavedParticles_I
    
    ! save the current index for the line we are working on
    iLineCurrent = iLine
  end subroutine disinter_line

  !============================================================================
  ! routine for saving current particles for restart. note that if not in 
  ! standalone PWOM mode particles must be disintered first
  subroutine write_restart_particle(iLine)
    use ModIoUnit, ONLY: UnitTmp_
    integer, intent(in) :: iLine
    character(len=150) :: NameRestart
    integer :: iParticle
    !---------------------------------------------------------------------------
    !set restart name
    write(NameRestart,"(a,i4.4,a)") &
         'PW/restartOUT/restart_particles_iline',iLineGlobal_I(iLine),'.dat'
    
    open(UnitTmp_, FILE=NameRestart,status="replace", form="unformatted")    
   
    WRITE (UnitTmp_) nParticle
    do iParticle=1,nParticle
       WRITE (UnitTmp_) Particles_I(iParticle)%iSpecies,&
            Particles_I(iParticle)%iCell,Particles_I(iParticle)%vpar,&
            Particles_I(iParticle)%vperp,Particles_I(iParticle)%Alt,&
            Particles_I(iParticle)%IsOpen,Particles_I(iParticle)%NumPerParticle
    enddo
    close(UnitTmp_)
    
  end subroutine write_restart_particle

  !============================================================================
  ! routine for reading current particles for restart. note that if not in 
  ! standalone PWOM mode particles must be buried after
  subroutine read_restart_particle(iLine)
    use ModIoUnit, ONLY: UnitTmp_
    integer, intent(in) :: iLine
    character(len=150) :: NameRestart
    integer :: iParticle
    
    !---------------------------------------------------------------------------
    !set restart name
    write(NameRestart,"(a,i4.4,a)") &
         'PW/restartIN/restart_particles_iline',iLineGlobal_I(iLine),'.dat'
    
    open(UnitTmp_, FILE=NameRestart,status="OLD", form="unformatted")    
   
    read (UnitTmp_) nParticle
    !reallocate Particles_I array to restart size
    if (allocated(Particles_I)) deallocate(Particles_I)
    allocate(Particles_I(nParticle))
    do iParticle=1,nParticle
       read(UnitTmp_) Particles_I(iParticle)%iSpecies,&
            Particles_I(iParticle)%iCell,Particles_I(iParticle)%vpar,&
            Particles_I(iParticle)%vperp,Particles_I(iParticle)%Alt,&
            Particles_I(iParticle)%IsOpen,Particles_I(iParticle)%NumPerParticle
    enddo
    close(UnitTmp_)
    

    call clean_zero_weight_particles

    
  end subroutine read_restart_particle
  
  !============================================================================
  ! remove any particle split al the way to zero weight 
  subroutine clean_zero_weight_particles
    integer,allocatable::Index_I(:)
    integer :: nZero, nParticleOld, iParticle,iParticleOld
    type(particle),allocatable ::ParticlesOld_I(:)
    !--------------------------------------------------------------------------
    
    !get number of zero weight particles
    allocate(Index_I(nParticle))
    Index_I=0
    !$OMP PARALLEL DO  
    do iParticle=1,nParticle
       if(Particles_I(iParticle)%NumPerParticle==0.0) then
          Index_I(iParticle)=1
       else
          Index_I(iParticle)=0
       endif
    end do
    !$OMP END PARALLEL DO

    nZero=sum(Index_I)
    
    !if no zero particles then return
    if(nZero==0) return

    !save old particle array information
    nParticleOld=nParticle
    allocate(ParticlesOld_I(nParticleOld))
    ParticlesOld_I=Particles_I
    deallocate(Particles_I)
       
    !allocate new Particles_I array
    nParticle=nParticleOld-nZero
    allocate(Particles_I(nParticle))
    
    !fill in particle array with non-zero particles
    iParticle=1
    do iParticleOld=1,nParticleOld
       if(ParticlesOld_I(iParticleOld)%NumPerParticle/=0) then
          Particles_I(iParticle)=ParticlesOld_I(iParticleOld)
          iParticle=iParticle+1
       endif
    end do

    deallocate(Index_I, ParticlesOld_I)
  end subroutine clean_zero_weight_particles
  !============================================================================
  ! combine a fluid solution and a particle solution with some percent from the 
  ! fluid and some from the particle. Idea is to sample the maxwellian 
  ! representing the fluid for one part and reduce the weight of the particles
  ! for the other part.
  subroutine combine_fluid_particle_cell(iCell,iSpecies,Density,Velocity,&
       Temperature,FracFluid)
    integer, intent(in) :: iCell,iSpecies
    real, intent(in) :: Density,Velocity,Temperature,FracFluid
    integer :: nParticlePerCellSaved
    integer :: nParticleInCell, iParticle
    !--------------------------------------------------------------------------
    
    nParticleInCell=nSortedParticle_II(iSpecies,iCell)
    
    ! reduce weight of particles to 1-FracFluid. This basically says the 
    ! existing particle solution accounts for some part of the density and 
    ! the fluid accounts for the rest
    do iParticle=1,nParticleInCell
       SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle=&
            (1.0-FracFluid)*&
            SortParticles_III(iParticle,iSpecies,iCell)%Particle%NumPerParticle
    enddo
    
    !save the existing nParticlePerCell_I(iSpecies)
    nParticlePerCellSaved=nParticlePerCell_I(iSpecies)

    !reduce the target number of particle per cell so our new fluid sample 
    !does not overwhelm our ability to join particles later
    nParticlePerCell_I(iSpecies)=nParticlePerCell_I(iSpecies)*0.5
    
    !sample the fluid and add those particles of the main particle array
    call sample_maxwellian_cell_boxmuller(iCell,iSpecies,&
                  FracFluid*Density,Velocity,Temperature)

    !reset the nParticle per cell
    nParticlePerCell_I(iSpecies)=nParticlePerCellSaved
  end subroutine combine_fluid_particle_cell

  !============================================================================
  ! create a source for the first computational cell based on the boundary cell 
  ! determine the associated face flux and inject particles. This is an 
  ! alternative to the standard approach of just sampling ghost cell.
  subroutine inject_source_boundary(iSpecies)
    integer,intent(in) :: iSpecies
    real:: SoundSpeed,Flux, FluxL, FluxR, DenToInject,VelToInject, TempToInject
    real:: SoundSpeedL,SoundSpeedR
    real:: Density,VelocityPar,VelocityPerp,Temp,TparIon,TperpIon,Hpar,Hperp,Pressure
    real,parameter :: gamma=1.66666666666666666666
    !--------------------------------------------------------------------------
        
    !Determine the sound speed
    SoundSpeedL=sqrt(gamma*cBoltzmannCGS*TemperatureBC_I(iSpecies)&
         /(Mass_I(iSpecies)))
    
    ! The flux is the density times the sound speed + bulk speed. This gives 
    ! the bulk and diffusive flux contribution from the left side of interface
    FluxL = (VelocityBC_I(iSpecies))*DensityBC_I(iSpecies)
    
    !get flux on right side of interface
    call calc_moments_cell(iSpecies,1,&
         Density,VelocityPar,VelocityPerp,Pressure,Temp,&
         TparIon,TperpIon,Hpar,Hperp)

    !Determine the sound speed
    SoundSpeedR=sqrt(gamma*cBoltzmannCGS*Temp&
         /(Mass_I(iSpecies)))
    
    
   ! The flux is the density times the sound speed + bulk speed. This gives 
    ! the bulk and diffusive flux contribution from the left side of interface
    FluxR = (VelocityPar)*Density

!    Flux=0.5*(FluxL+FluxR)&
!         -0.5*max(VelocityBC_I(iSpecies)+SoundSpeedL,VelocityPar+SoundSpeedR)&
!         *(Density-DensityBC_I(iSpecies))

    Flux=0.5*(FluxL+FluxR)&
         -0.5*max(abs(VelocityBC_I(iSpecies)+SoundSpeedL),abs(VelocityPar+SoundSpeedR))&
         *(Density-DensityBC_I(iSpecies))


    ! get the number density to inject
    DenToInject=Flux * AreaInlet * DtMove / Volume_G(1) 

    ! get the velocity to inject
    VelToInject=VelocityBC_I(iSpecies)!0.5*(abs(FluxL)*VelocityBC_I(iSpecies)+abs(FluxR)*VelocityPar)/(abs(FluxL)+abs(FluxR))
    
    ! get the Temperature to inject
    TempToInject=TemperatureBC_I(iSpecies)

    !inject the particles. Only do this if flux is > 0. Note that interface 
    !flux out of the cell is 
    if (Flux>0) then
       call sample_maxwellian_cell_boxmuller(1,iSpecies,&
            DenToInject,VelToInject,TempToInject)
    endif
  end subroutine inject_source_boundary
  !============================================================================
  ! 
  subroutine get_from_particles(nAltIn,nSpeciesIn,AltIn_C,&
       DensityOut_IC,VelocityOut_IC,TemperatureOut_IC,HeatFluxOut_IC)
    use ModInterpolate, only: linear
    integer, intent(in) :: nAltIn,nSpeciesIn
    real,    intent(in) :: AltIn_C(nAltIn)
    
    !for each species, the state variables as a function of alt
    real,    intent(out) :: DensityOut_IC(nSpeciesIn,nAltIn)
    real,    intent(out) :: VelocityOut_IC(nSpeciesIn,nAltIn)
    real,    intent(out) :: TemperatureOut_IC(nSpeciesIn,nAltIn)
    real,    intent(out) :: HeatFluxOut_IC(nSpeciesIn,nAltIn)
    
    integer :: iSpecies, nParticleInCell, iAlt
    real :: density, Velocity,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,&
         Tperp,Temperature,Hpar,Hperp
    real, allocatable :: Density_IC(:,:),Velocity_IC(:,:),Temperature_IC(:,:)
    real, allocatable :: HeatFlux_IC(:,:)
    
    !min values for density and temperature
    !    real, parameter :: DensityMin=1e-4 !cm-3
    real, parameter :: DensityMin=1e-8 !cm-3
    !    real, parameter :: TemperatureMin=100.0 !k
        real, parameter :: TemperatureMin=1.0 !k
!    logical,parameter :: UseSmooth=.true.
    !---------------------------------------------------------------------------
    
    allocate(Density_IC(nSpeciesIn,nAlt),Velocity_IC(nSpeciesIn,nAlt),&
         Temperature_IC(nSpeciesIn,nAlt),HeatFlux_IC(nSpeciesIn,nAlt))
    
    do iAlt=1,nAlt
       do iSpecies=1,nSpeciesIn
 !         if (UseSmooth) then
 !            call calc_moments_cell_smooth(iSpecies,iAlt,&
 !                 density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp)
 !            Density_IC(iSpecies,iAlt)    =density
 !            Velocity_IC(iSpecies,iAlt)   =uBulkPar
 !            Temperature_IC(iSpecies,iAlt)=Temp
 !         else
             nParticleInCell=nSortedParticle_II(iSpecies,iAlt)
             if(nParticleInCell==0) then
                Density_IC(iSpecies,iAlt)    =0.0
                Velocity_IC(iSpecies,iAlt)   =0.0
                Temperature_IC(iSpecies,iAlt)=0.0
                HeatFlux_IC(iSpecies,iAlt)=0.0
             else
                ! get moments in each cell so we can later interpolate to fluid 
                !call calc_moments_cell(iSpecies,iAlt,&
                !     density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp)
                !note that weighed calculation fall appart next to ghost cell
                ! so revert to regualr moment calculation
                if (iAlt==1 .or. iAlt==nAlt) then
                   call calc_moments_cell(iSpecies,iAlt,&
                        density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp,Hpar,Hperp)
                else
                   call calc_moments_cell_weighted(iSpecies,iAlt,&
                        density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp,Hpar,Hperp)
                endif
                Density_IC(iSpecies,iAlt)    =density
                Velocity_IC(iSpecies,iAlt)   =uBulkPar
                Temperature_IC(iSpecies,iAlt)=Temp
                HeatFlux_IC(iSpecies,iAlt)=Hpar
             endif
!          end if
       enddo
    end do
   
    !Now interpolate results from particle grid onto fluid grid
    do iAlt=1,nAltIn
       do iSpecies=1,nSpeciesIn
          
          !check bounds 
          if(AltIn_C(iAlt)<Alt_G(1) .or. AltIn_C(iAlt)>Alt_G(nAlt)) then
             DensityOut_IC(iSpecies,iAlt) = DensityMin
             VelocityOut_IC(iSpecies,iAlt)    = 0.0
             TemperatureOut_IC(iSpecies,iAlt) = TemperatureMin
             HeatFlux_IC(iSpecies,iAlt) = 0.0
          else
             !interpolate
             DensityOut_IC(iSpecies,iAlt)  = &
                  max(linear(Density_IC(iSpecies,:),1,     &
                  nAlt,AltIn_C(iAlt),Alt_G(1:nAlt)),DensityMin)
             VelocityOut_IC(iSpecies,iAlt) = &
                  linear(Velocity_IC(iSpecies,:),1,     &
                  nAlt,AltIn_C(iAlt),Alt_G(1:nAlt))
             TemperatureOut_IC(iSpecies,iAlt) = &
                  max(linear(Temperature_IC(iSpecies,:),1,     &
                  nAlt,AltIn_C(iAlt),Alt_G(1:nAlt)),TemperatureMin)
             HeatFluxOut_IC(iSpecies,iAlt) = &
                  linear(HeatFlux_IC(iSpecies,:),1,     &
                  nAlt,AltIn_C(iAlt),Alt_G(1:nAlt))
          endif
          !if (iSpecies==1) write(*,*) iSpecies,AltIn_C(iAlt),&
          !     DensityOut_IC(iSpecies,iAlt),VelocityOut_IC(iSpecies,iAlt),&
          !     TemperatureOut_IC(iSpecies,iAlt)
          
       enddo
    enddo

    deallocate(Density_IC,Velocity_IC,&
         Temperature_IC,HeatFlux_IC)    
    
  end subroutine get_from_particles

  !============================================================================
  ! 
  subroutine put_to_particles(nAltIn,nSpeciesIn,AltIn_C,&
       DoInitAlt,DensityIn_IC,VelocityIn_IC,TemperatureIn_IC,EfieldIn_C)
    use ModInterpolate, only: linear
    integer, intent(in) :: nAltIn,nSpeciesIn
    real,    intent(in) :: AltIn_C(nAltIn)
    logical, intent(in) :: DoInitAlt
    !for each species, the state variables as a function of alt
    real,    intent(in) :: DensityIn_IC(nSpeciesIn,nAltIn)
    real,    intent(in) :: VelocityIn_IC(nSpeciesIn,nAltIn)
    real,    intent(in) :: TemperatureIn_IC(nSpeciesIn,nAltIn)
    real, optional,   intent(in) :: EfieldIn_C(nAltIn)

    !local variables
    real :: Density,Velocity,Temperature
    real :: uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp
    integer :: iAlt, iSpecies
    logical, save :: IsFirstCall=.true.
    !--------------------------------------------------------------------------
    
    !interpolate the incomming efield 
    If(Present(EfieldIn_C)) then
       do iAlt=-1,nAlt+1
          if (Alt_G(iAlt)<AltIn_C(1) .or. Alt_G(iAlt)>AltIn_C(nAltIn) ) then
             Efield_G(iAlt)=0.0
          else
             Efield_G(iAlt)= linear(EfieldIn_C(:),1,nAltIn,Alt_G(iAlt),AltIn_C)
          endif
       enddo
    endif

    !when DoInitAlt then sample particles at each altitude according to 
    !the initial PWOM condition
    if (DoInitAlt) then
       do iSpecies=1,nSpecies
          do iAlt=1,nAlt
             !interpolate statevariables to alt position and sample
             Density  = linear(DensityIn_IC(iSpecies,:),1,     &
                  nAltIn,Alt_G(iAlt),AltIn_C)
             Velocity = linear(VelocityIn_IC(iSpecies,:),1,    &
                  nAltIn,Alt_G(iAlt),AltIn_C)
             Temperature = linear(TemperatureIn_IC(iSpecies,:),1, &
                  nAltIn,Alt_G(iAlt),AltIn_C)
             
             !now Sample
             if (IsVerboseParticle) &
                  write(*,*) 'sample iSpecies,Alt', &
                  iSpecies,Alt_G(iAlt),Density,Velocity,Temperature
             call sample_maxwellian_cell_boxmuller(iAlt,iSpecies,&
                  Density,Velocity,Temperature)
          enddo
       enddo
    endif

    ! Get the density, velocity, and Temperature
    !write(*,*) 'TEST',Alt_G(0)
    !if (IsFirstCall)then
!    call sort_particles
    !endif
    do iSpecies=1,nSpecies
       DensityBC_I(iSpecies)     = linear(DensityIn_IC(iSpecies,:),1,     &
            nAltIn,Alt_G(0),AltIn_C)
       VelocityBC_I(iSpecies)    = linear(VelocityIn_IC(iSpecies,:),1,    &
            nAltIn,Alt_G(0),AltIn_C)
       TemperatureBC_I(iSpecies) = linear(TemperatureIn_IC(iSpecies,:),1, &
            nAltIn,Alt_G(0),AltIn_C)

       !       !kludge try characteristic based BCs
!       call calc_moments_cell(iSpecies,1,&
!            Density,uBulkPar,uBulkPerp,Pressure,Temp,Tpar,Tperp)
!       
!       if (VelocityBC_I(iSpecies)+uBulkPar>0.0) then
!          !inlet BC impose density and velocity, float pressure
!          TemperatureBC_I(iSpecies)=Temp*Density/DensityBC_I(iSpecies)
!       else
!          !outlet BC impose pressure, float density and velocity
!          DensityBC_I(iSpecies)=Density
!          VelocityBC_I(iSpecies)=uBulkPar
!       endif
       
       If(UseOverlapRegion) then
          !now fill fluid values in overlap cells
          do iAlt=1,nOverlap
             DensityOverlap_IC(iSpecies,iAlt)= &
                  linear(DensityIn_IC(iSpecies,:),1,nAltIn,Alt_G(iAlt),AltIn_C)
             VelocityOverlap_IC(iSpecies,iAlt)= &
                  linear(VelocityIn_IC(iSpecies,:),1,nAltIn,Alt_G(iAlt),AltIn_C)
             TemperatureOverlap_IC(iSpecies,iAlt)= &
                  linear(DensityIn_IC(iSpecies,:),1,nAltIn,Alt_G(iAlt),AltIn_C)
          enddo
       endif
       !
    enddo
!    call sort_particles
!    call plot_profile
!    write(*,*) 'TEST STOP'
  end subroutine put_to_particles
  !============================================================================
  ! advance the particle solution for some DtAdvance
  subroutine run_particles(DtAdvance,IsCuspOrAurora,SmLat)
    use ModConst,           ONLY: cElectronCharge
    use ModPlanetConst,     ONLY: Planet_, NamePlanet_I, rPlanet_I  
    real, intent(in) :: DtAdvance,SmLat
    logical, intent(in) :: IsCuspOrAurora

    integer :: iAlt, iSpecies, jSpecies, nTime,iTime,iAltPlot,iCollide
    integer, parameter :: iAltBC=0
    real :: TimeAdvance
    real :: WaveCoef, fci, SpectralIndex, E2waveRef,fWaveRef, rWaveRef
    character(len=100):: TypeGrid

    real,parameter :: cGtoKg = 1.0e-3, cMtoCm=1e2
    !---------------------------------------------------------------------------
    
    TimeAdvance=0.0
    
    if(IsVerboseParticle) write(*,*) 'nParticle=',nParticle
    TIMELOOP:do 
       !check stopping condition
       if (TimeAdvance >=DtAdvance) exit TIMELOOP
       
       !set timestep
       DtMove=min(DtAdvance-TimeAdvance,DtMoveMax)
       
       ! stop if within tolerance of stopping time
       if(DtMove <1.0e-6) exit TIMELOOP
       
       !set collision timestep
       DtCollide=DtMove/real(nCollide)

       !resample ghost cell
       do iSpecies=1,nSpecies
          call timing_start('sample_maxwellian_cell_boxmuller')
          call sample_maxwellian_cell_boxmuller(iAltBC,iSpecies,&
               DensityBC_I(iSpecies),VelocityBC_I(iSpecies),&
               TemperatureBC_I(iSpecies))    
          call timing_stop('sample_maxwellian_cell_boxmuller')
       end do

       !set overlap region
       if(UseOverlapRegion) then
          do iSpecies=1,nSpecies
             do iAlt=1,nOverlap
                !Sort the particles
                call timing_start('sort_particles')
                call sort_particles
                call timing_stop('sort_particles')
                
                call combine_fluid_particle_cell(iAlt,iSpecies,&
                     DensityOverlap_IC(iSpecies,iAlt),&
                     VelocityOverlap_IC(iSpecies,iAlt),&
                     TemperatureOverlap_IC(iSpecies,iAlt),FluidFrac_C(iAlt))
             enddo
          enddo
       endif

       !push the particles
       call timing_start('push_guiding_center')
       call push_guiding_center
       call timing_stop('push_guiding_center')

       !Sort the particles
       call timing_start('sort_particles')
       call sort_particles
       call timing_stop('sort_particles')
 
       
       !apply the collisions
       do iAlt=1,nAlt
          !check if we are above altitude where collisions can be neglected
          !if (Alt_G(iAlt)>7500.0e5) exit
          if (iAlt==nAlt-1) exit
          !write(*,*) 'iAlt',iAlt
          do iCollide=1,nCollide
             do iSpecies=1,nSpecies
                do jSpecies=iSpecies,nSpecies
                   call timing_start('apply_coulomb_collision')
                   call apply_coulomb_collision(iAlt,iSpecies,jSpecies)
                   call timing_stop('apply_coulomb_collision')
                enddo
             enddo
          enddo
       enddo
       
       !apply the WPI
       if (UseWPI) then
          !set appropriate WPI coefficients
          select case(TypeWPI)
          case('Barakat')
             !Uses the coeficients published by Barghouthi 1997 and used in
             !Barakat and Schunk 2001
             !from DE data survey. This is only valid for Earth
             if(IsCuspOrAurora) then
                Dperp_I(O_)=6.94e5   !auroral/cusp value
                Dexp_I(O_)=13.3
                
                Dperp_I(H_)=4.45e7   !auroral/cusp value
                Dexp_I(H_)=7.95

                Dperp_I(He_)=0.0  
                Dexp_I(He_)=7.95
             else
                Dperp_I(O_)=9.55e2  !polar cap value
                Dexp_I(O_)=13.3
                
                Dperp_I(H_)=5.77e3  !polar cap value
                Dexp_I(H_)=7.95

                Dperp_I(He_)=0.0  
                Dexp_I(He_)=7.95
             endif
             rWaveRef= rPlanet_I(Planet_)*cMtoCm
          case('General')
             ! uses approach described by Crew et al. [1990], as well as
             ! Retterer et al [1987],...

             if(IsCuspOrAurora) then
                SpectralIndex=SpectralIndexAur
                E2waveRef    = E2waveRefAur
                fWaveRef     = fWaveRefAur
                rWaveRef   = rWaveRefAur
             else
                SpectralIndex=SpectralIndexCap
                E2waveRef    = E2waveRefCap
                fWaveRef     = fWaveRefCap
                rWaveRef   = rWaveRefCap
             endif
             
             WaveCoef = E2waveRef*fWaveRef**SpectralIndex
             do iSpecies=1,nSpecies
                call get_fci(rWaveRef, iSpecies,SmLat, fci)
                Dperp_I(iSpecies) = &
                     (FracLeftHand*cElectronCharge**2)&
                     /(4.0*Mass_I(iSpecies)*cGtoKg) &
                     * WaveCoef*fci**(-SpectralIndex)*(cMtoCm**2)
                Dexp_I(iSpecies) = 3.0*SpectralIndex
             enddo
             
          end select
             
          call timing_start('apply_wave_particle_interaction')
          call apply_wave_particle_interaction(rWaveRef)
          call timing_stop('apply_wave_particle_interaction')
       endif

       
       ! split and join particles to be within 5% of the target particles/cell
       ! for speed up, only split and join every DtSplitJoin 
       if (floor((Time+1.0e-5)/DtSplitJoin) &
            /=floor((Time+1.0e-5-DtMove)/DtSplitJoin) )then 
          call timing_start('split_join_particles')
          call split_join_particles(0.05)
          call timing_stop('split_join_particles')
          
          call timing_start('sort_particles')
          call sort_particles
          call timing_stop('sort_particles')
       endif

       !advance the time
       Time=Time+DtMove
       TimeAdvance=TimeAdvance+DtMove
       
       !plot profile of moments
       if (floor((Time+1.0e-5)/DtSaveProfile) &
            /=floor((Time+1.0e-5-DtMove)/DtSaveProfile) )then 
          do iSpecies=1,nSpecies
             call plot_profile(iSpecies)
          enddo
       endif

       !plot DF                                                                 
       if (floor((Time+1.0e-5)/DtSaveDF) &
            /=floor((Time+1.0e-5-DtMove)/DtSaveDF) )then
          do iSpecies=1,nSpecies
             do iAltPlot=1,nSaveDfAlts
                call plot_distribution_cell(iSpecies,iAltsDF_I(iAltPlot))
             enddo
          enddo
       endif

    enddo TIMELOOP
    
    
  end subroutine run_particles

  !============================================================================
  ! calculate the ion cyclotron frequency in hertz
  subroutine get_fci(AltRef, iIon,SmLat, fci)
    use ModPlanetConst,     ONLY: Earth_,DipoleStrengthPlanet_I,rPlanet_I
    use ModNumConst,        ONLY: cDegToRad
    use ModConst,           ONLY: cElectronCharge
    real, intent(in):: AltRef !incomming reference alt [cm]
    real, intent(in):: SmLat !Lat in SM at foot of field line [degrees]
    integer, intent(in):: iIon !index of ion species
    real, intent(out)   :: fci

    real  :: B0ref
    
    real, parameter :: cCmToM=1.0e-2, cGtoKg=1.0e-3
    real    :: Lshell, rPlanet, dipmom, Lat
    real    :: rRef ! reference radius
    !--------------------------------------------------------------------------
    
    rPlanet = rPlanet_I(Earth_)                            ! planet's radius (m)
    dipmom  = abs(DipoleStrengthPlanet_I(Earth_)*rPlanet**3)  ! planet's dipole 
    
    !set the reference radius
    rRef = rPlanet+AltRef*cCmToM
    ! find corresponding l-shell
    Lshell = 1.0/(cos(SmLat*cDegToRad))**2.0
    
    ! find corresponding latitude for location on l-shell
    Lat = acos(sqrt(rRef*cCmToM/(Lshell*rPlanet)))
    
    ! get the magnetic field of the reference altitude
    B0ref = &
         dipmom*sqrt(1+3.0*(sin(Lat))**2.0)/(rRef)**3.0

    !calculate the gyro freq, fci
    fci = cElectronCharge/(Mass_I(iIon)*cGtoKg)*B0ref
    
  end subroutine get_fci

  !============================================================================
  ! unit test subroutine for sampling
  subroutine test_sample
    integer :: nAltIn, iCell, iSpecies
    real :: AltMin, AltMax, Density, uBulk, Temperature
    character(len=100):: TypeGrid
    real :: densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,&
                     TempTmp,TparTmp,TperpTmp,Hpar,Hperp
    !--------------------------------------------------------------------------
    !set global line info
    allocate(iLineGlobal_I(nLine))
    iLineGlobal_I(1)=1
    iLineCurrent=1
    nParticlePerCell_I(1)=5000

    nAltIn=2
    AltMin=1000.0e5
    AltMax=1020.0e5
    TypeGrid='Uniform'
    iCell = 1
    iSpecies=1
    Density=1e5
    uBulk=100000.0
    Temperature=1000.0
    write(*,*) 'init_particle'
    call init_particle(nAltIn,AltMin,AltMax,TypeGrid)
    
    write(*,*) 'sample_maxwellian_cell'
    
       call sample_maxwellian_cell_boxmuller(iCell-1,iSpecies,Density,uBulk,Temperature)
    call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature)

    call sample_maxwellian_cell_boxmuller(iCell+1,iSpecies,Density,uBulk,Temperature)

    call sample_maxwellian_cell_boxmuller(iCell,iSpecies+1,Density,uBulk,Temperature)

    call sample_maxwellian_cell_boxmuller(iCell+1,iSpecies+1,Density,uBulk,Temperature)
    
    write(*,*) 'test pointer extraction'
!    Particles_I%vperp=Particles_I%vperp*2.0
    call sort_particles
    !now plot from sorted
    call plot_distribution_cell(iSpecies,iCell)
    call plot_distribution_cell(iSpecies,iCell+1)
!    call plot_distribution_cell(iSpecies,iCell,nSortedParticle_II,&
!         SortParticles_III(iSpecies)%CellParticle)
    call calc_moments_cell_weighted(iSpecies,iCell,&
         densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,&
         TempTmp,TparTmp,TperpTmp,Hpar,Hperp)    
    write(*,*) 'moments',densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,&
         TempTmp,TparTmp,TperpTmp

    call calc_moments_cell(iSpecies,iCell,&
         densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,&
         TempTmp,TparTmp,TperpTmp,Hpar,Hperp)    
    write(*,*) 'moments',densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,&
         TempTmp,TparTmp,TperpTmp


    call plot_profile(iSpecies)
  end subroutine test_sample

  !============================================================================
  ! unit test subroutine for splitting and joining
  subroutine test_split_join
    integer :: nAltIn, iCell, iSpecies
    real :: AltMin, AltMax, Density, uBulk, Temperature
    integer :: nSplit, nJoin,iCount
    character(len=100):: TypeGrid
    real :: densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempTmp,&
         TparTmp,TperpTmp,Hpar,Hperp    
    !--------------------------------------------------------------------------
    !set global line info
    allocate(iLineGlobal_I(nLine))
    iLineGlobal_I(1)=1
    iLineCurrent=1
    
    nAltIn=2
    AltMin=1000.0e5
    AltMax=1020.0e5
    TypeGrid='Uniform'
    iCell = 1
    iSpecies=1
    Density=1e5
    uBulk=1.0e5
    Temperature=1000.0
    
    write(*,*) 'init_particle'
    call init_particle(nAltIn,AltMin,AltMax,TypeGrid)
    
    write(*,*) 'sample_maxwellian_cell'
    
    !set target particles per cell 
    nParticlePerCell_I(iSpecies)=10000

    call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature)
    call sort_particles
    !now plot from sorted
    call plot_distribution_cell(iSpecies,iCell)
    
    write(*,*) 'nParticle at start:', nParticle

    ! test repeated splitting (10 times)
    do iCount=1,50
       !Split 5% of particles
       nSplit=floor(nParticle*0.05)
       
       !update target particles per cell 
       nParticlePerCell_I(iSpecies)=nParticlePerCell_I(iSpecies)+nSplit

       !now plot from sorted
       call split_join_particles(0.03)
       call sort_particles
       
       call plot_distribution_cell(iSpecies,iCell)
       call calc_moments_cell(iSpecies,iCell,&
            densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempTmp,&
            TparTmp,TperpTmp,Hpar,Hperp) 
       write(*,*) 'After Split:'
       write(*,*) '   ','nParticle   ', nParticle
       write(*,*) '   ','densityTmp  ', densityTmp   
       write(*,*) '   ','uBulkParTmp ', uBulkParTmp  
       write(*,*) '   ','uBulkPerpTmp', uBulkPerpTmp 
       !write(*,*) '   ','PressureTmp ', PressureTmp  
       write(*,*) '   ','TempTmp     ', TempTmp      
       write(*,*) '   ','TparTmp     ', TparTmp      
       write(*,*) '   ','TperpTmp    ', TperpTmp     
    enddo
    
    ! test repeated joining (10 times)
    do iCount=1,50
       !Join 5% of particles
       nJoin=floor(nParticle*0.05)


       !update target particles per cell 
       nParticlePerCell_I(iSpecies)=nParticlePerCell_I(iSpecies)-nJoin
       
       call split_join_particles(0.03)

       call sort_particles
       !now plot from sorted
       call plot_distribution_cell(iSpecies,iCell)
       
       call calc_moments_cell(iSpecies,iCell,&
            densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempTmp,&
            TparTmp,TperpTmp,Hpar,Hperp) 
       write(*,*) 'After Join:'
       write(*,*) '   ','nParticle   ', nParticle
       write(*,*) '   ','densityTmp  ', densityTmp   
       write(*,*) '   ','uBulkParTmp ', uBulkParTmp  
       write(*,*) '   ','uBulkPerpTmp', uBulkPerpTmp 
       !write(*,*) '   ','PressureTmp ', PressureTmp  
       write(*,*) '   ','TempTmp     ', TempTmp      
       write(*,*) '   ','TparTmp     ', TparTmp      
       write(*,*) '   ','TperpTmp    ', TperpTmp     
    enddo
    

!    call plot_profile(iSpecies)
  end subroutine test_split_join


  !============================================================================
  ! unit test subroutine for sampling
  subroutine test_pusher
    integer :: nAltIn, iCell, iSpecies,nTime,iTime,iAltBC
    real :: AltMin, AltMax, Density, uBulk, Temperature,DtSavePlot
    character(len=100):: TypeGrid
    !--------------------------------------------------------------------------
    !set global line info
    allocate(iLineGlobal_I(nLine))
    iLineGlobal_I(1)=1
    iLineCurrent=1
    DtMove=1.0e0
    nTime=floor(1000./DtMove)
    
    !nTime=1000
!    nTime=10000
    DtSavePlot=10.
    nAltIn=10
    AltMin=1000.0e5
    AltMax=1200.0e5
    TypeGrid='Uniform'
    iAltBC = 0
    iSpecies=1
    Density=1e2
!    uBulk=1.0e5
    uBulk=0.0
    Temperature=3000.0
    nParticlePerCell_I(1)=5000
    write(*,*) 'init_particle'
    
    call timing_start('init_particle')
    call init_particle(nAltIn,AltMin,AltMax,TypeGrid)
    call timing_stop('init_particle')


    
    DensityBC_I(iSpecies)     = Density
    VelocityBC_I(iSpecies)    = uBulk
    TemperatureBC_I(iSpecies) = Temperature
    
    Efield_G(:)=0.0
    
    write(*,*) 'initialize the ghost cell'
    call timing_start('sample_maxwellian_cell_boxmuller')
    call sample_maxwellian_cell_boxmuller(iAltBC,iSpecies,Density,uBulk,Temperature)
!    call inject_source_boundary(iSpecies)
!    call combine_fluid_particle_cell(1,iSpecies,DensityBC_I(iSpecies),&
!         VelocityBC_I(iSpecies),TemperatureBC_I(iSpecies),0.5)

    call timing_stop('sample_maxwellian_cell_boxmuller')
    !push the guiding center 100 times and reinitialize ghost cell each time
    do iTime=1,nTime
       call timing_start('push_guiding_center')
       !call push_guiding_center
       call push_guiding_center_rk4
       call timing_stop('push_guiding_center')
       write(*,*) 'test1'
       !advance the time
       Time=Time+DtMove
       
       DtCollide=DtMove
       !resample ghost cell
       call timing_start('sample_maxwellian_cell_boxmuller')
       call sample_maxwellian_cell_boxmuller(iAltBC,iSpecies,Density,uBulk,Temperature)    
       !call inject_source_boundary(iSpecies)
!       write(*,*) 'test2'
!       call combine_fluid_particle_cell(1,iSpecies,DensityBC_I(iSpecies),&
!            VelocityBC_I(iSpecies),TemperatureBC_I(iSpecies),0.5)
       call timing_stop('sample_maxwellian_cell_boxmuller')

       !Sort the particles
       call timing_start('sort_particles')
       call sort_particles
       call timing_stop('sort_particles')
       
       call timing_start('split_join_particles')
       call split_join_particles(0.05)
       call timing_stop('split_join_particles')
       
       !Sort the particles
       call timing_start('sort_particles')
       call sort_particles
       call timing_stop('sort_particles')

       !kludge apply collisions here
       if (Time>50.) then
          do iCell=1,nAlt
             call apply_coulomb_collision(iCell,1,1)
             !call apply_coulomb_collision(iCell,2,2)
             !interspecies collisions
             !call apply_coulomb_collision(iCell,1,2)
          end do
       endif

       !plot profile of moments
       if (floor((Time+1.0e-5)/DtSavePlot) &
            /=floor((Time+1.0e-5-DtMove)/DtSavePlot) )then 
          call plot_profile(iSpecies)
       endif
       
       write(*,*) 'nParticle=',nParticle
    enddo
!    do iCell=1,nAlt
!       call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature*2*(iCell))    
!    enddo
    


    !now plot distribution function at each alt
    do iCell=0,nAltIn
       write(*,*)nSortedParticle_II(iSpecies,iCell)
       call plot_distribution_cell(iSpecies,iCell)
    enddo

  end subroutine test_pusher
  !============================================================================
  ! unit test for the application of coulomb collisions
  ! use the relaxation problem outlined in Miller and Combi as a test case
  ! this involves electrons and ions with reduced mass ratio so the 
  ! initialization needs to be somewhat overwritten.
  subroutine test_coulomb_collision
    use ModNumConst, ONLY: cPi,cTwoPi
    use ModConst, ONLY: cEps,cElectronCharge,cBoltzmann
    integer :: nAltIn, iCell,nTime,iTime,iAltBC
    integer :: iSpeciesIon,iSpeciesElec
    integer :: iSpecies,jSpecies
    real :: TempIon0,TempElec0,TempIon,TempElec
    real :: AltMin, AltMax, Density, uBulk, DtSavePlot
    character(len=100):: TypeGrid
    real :: densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempTmp
    real,parameter::cCm3ToM3=1e6,cGtoKg=1e-3,cCmToM=1e-2,ckToEv=.00008617328
    real,parameter ::cElectronChargeCGS=1.602176487e-20
    real :: nu0, nueq,dt,CoulLog
    real :: TparIon,TperpIon,TparElec,TperpElec,Hpar,Hperp
    !--------------------------------------------------------------------------
!    nTime=10
!    nTime=1000
    nTime=10000
    DtSavePlot=10
    nAltIn=10
    AltMin=1000.0e5
    AltMax=1200.0e5
    TypeGrid='Uniform'
    iAltBC = 0
    iSpeciesIon=1
    iSpeciesElec=2
    Density=5e5
    uBulk=0.0
!    TempIon0 =200.0/ckToEv!10000.0
    TempIon0 =10000.0
    TempElec0=0.5*TempIon0
    !TempElec0=TempIon0
    

    write(*,*) 'init_particle'
    Time=0.0
    
    call timing_start('init_particle')
    call init_particle(nAltIn,AltMin,AltMax,TypeGrid)
    call timing_stop('init_particle')

    !set global line info
    allocate(iLineGlobal_I(nLine))
    iLineGlobal_I(1)=1
    iLineCurrent=1

    !overwrite masses from init since we are using H+ and e- 
    Mass_I(iSpeciesIon)=1.66054e-24
    Mass_I(iSpeciesElec)=0.25*Mass_I(iSpeciesIon)

    !recalculate the reduced mass with these new masses
    do iSpecies=1,nSpecies
       do jSpecies=1,nSpecies
          ReducedMass_II(iSpecies,jSpecies)=Mass_I(iSpecies)*Mass_I(jSpecies)&
               /(Mass_I(iSpecies)+Mass_I(jSpecies))
       enddo
    enddo
    
    CoulLog=23.0&
         -log((Mass_I(iSpeciesIon)+Mass_I(iSpeciesElec))&
         /(Mass_I(iSpeciesIon)*TempIon0&
         +Mass_I(iSpeciesElec)*TempElec0)&
         *sqrt(cCm3ToM3*Density/TempIon0&
         +cCm3ToM3*Density/TempElec0))

    write(*,*) CoulLog,cElectronCharge,cCm3ToM3,Density,&
            cPi,cEps,cGtoKg,Mass_I(iSpeciesElec),TempElec0

    !set the relaxation frequencies
    nu0= cElectronCharge**4*cCm3ToM3*Density*CoulLog&
            /(8.*sqrt(2.0)*cPi*cEps**2*(cGtoKg*Mass_I(iSpeciesElec))**0.5&
            *(cBoltzmann*TempElec0)**1.5)

    nueq=(8.0/(3.0*sqrt(cPi)))*(Mass_I(iSpeciesElec)/Mass_I(iSpeciesIon))&
         *(1.0+(Mass_I(iSpeciesElec)/Mass_I(iSpeciesIon))&
         *(TempIon0/TempElec0))**(-1.5)*nu0

    DtCollide=.001/nu0
    !DtCollide=.01

    write(*,*) 'nu0,nueq,DtCollide',nu0,nueq,DtCollide
    
    write(*,*) 'initialize the ghost cell'
    call timing_start('sample_maxwellian_cell_boxmuller')
    nParticlePerCell_I(iSpecies)=15000    
    call sample_maxwellian_cell_boxmuller(iAltBC,iSpeciesIon,0.5*Density,uBulk,TempIon0)
    nParticlePerCell_I(iSpecies)=30000    
    call sample_maxwellian_cell_boxmuller(iAltBC,iSpeciesIon,0.5*Density,uBulk,TempIon0)

    nParticlePerCell_I(iSpecies)=15000    
    call sample_maxwellian_cell_boxmuller(iAltBC,iSpeciesElec,0.5*Density,uBulk,TempElec0)
    nParticlePerCell_I(iSpecies)=30000    
    call sample_maxwellian_cell_boxmuller(iAltBC,iSpeciesElec,0.5*Density,uBulk,TempElec0)
    call timing_stop('sample_maxwellian_cell_boxmuller')
    
    write(*,*)'nParticle=',nParticle
    !Particles_I%vperp=Particles_I%vperp*2.0
    
    !Sort the particles
    call timing_start('sort_particles')
    call sort_particles
    call timing_stop('sort_particles')

    !plot initial distirbution before collisions
    call plot_distribution_cell(iSpeciesIon,iAltBC)
    call plot_distribution_cell(iSpeciesElec,iAltBC)

    call calc_moments_cell(iSpeciesIon,iAltBC,&
         densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempIon,&
         TparIon,TperpIon,Hpar,Hperp)
    write(*,*) 'Initial paramerters Ion:'
    write(*,*) '   ',densityTmp,uBulkParTmp,PressureTmp,TempIon

    call calc_moments_cell(iSpeciesElec,iAltBC,&
         densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempElec,&
         TparElec,TperpElec,Hpar,Hperp)
    write(*,*) 'Initial paramerters Elec:'
    write(*,*) '   ',densityTmp,uBulkParTmp,PressureTmp,TempElec
    

    write(*,*) 'start: nu0*Time, (TIon-Te)/(TIon0-Te0), theory'
!    write(*,*) nu0*Time,(TempIon-TempElec)/(TempIon0-TempElec0),exp(-2.0*nueq*Time)

    !push the guiding center 100 times and reinitialize ghost cell each time
    do iTime=1,nTime
       call timing_start('apply_coulomb_collision')
       !self collisions
       call apply_coulomb_collision(iAltBC,1,1)
       call apply_coulomb_collision(iAltBC,2,2)
       !interspecies collisions
       call apply_coulomb_collision(iAltBC,1,2)
       call timing_stop('apply_coulomb_collision')
       
       !advance the time
       Time=Time+DtCollide
       
       !Sort the particles
       call timing_start('sort_particles')
       call sort_particles
       call timing_stop('sort_particles')

       call calc_moments_cell(iSpeciesIon,iAltBC,&
            densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempIon,&
            TparIon,TperpIon,Hpar,Hperp)
       call calc_moments_cell(iSpeciesElec,iAltBC,&
            densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempElec,&
            TparElec,TperpElec,Hpar,Hperp)

    write(*,*) nu0*Time,(TempIon-TempElec)/(TempIon0-TempElec0),exp(-2.0*nueq*Time)
!    write(*,*) TempIon,TempElec

!    if (floor((Time+1.0e-5)/(10.*DtCollide)) &
!         /=floor((Time+1.0e-5-DtCollide)/(10.*DtCollide)) )then 
!       call plot_distribution_cell(iSpeciesIon,iAltBC)
!       call plot_distribution_cell(iSpeciesElec,iAltBC)
!    endif

       !plot profile of moments
!       if (floor((Time+1.0e-5)/DtSavePlot) &
!            /=floor((Time+1.0e-5-DtCollide)/DtSavePlot) )then 
!          !call plot_distribution_cell(iSpeciesIon,iAltBC)
!          
!          call calc_moments_cell(iSpeciesIon,iAltBC,&
!               densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempTmp)
!          write(*,*) densityTmp,uBulkParTmp,TempTmp
!       endif

    enddo
!    do iCell=1,nAlt
!       call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature*2*(iCell))    
!    enddo
    


  end subroutine test_coulomb_collision

  !============================================================================
  ! unit test for the wpi. Basically apply the wpi repeatedly for a stationary 
  ! distribution and see perpendicular heating rate
  subroutine test_wpi
    use ModNumConst, ONLY: cPi,cTwoPi
    use ModConst, ONLY: cEps,cElectronCharge,cBoltzmann
    integer :: nAltIn, iCell,nTime,iTime,iAltBC
    integer :: iSpecies
    real :: AltMin, AltMax, Density, uBulk, DtSavePlot
    character(len=100):: TypeGrid
    real :: densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempTmp
    real,parameter::cCm3ToM3=1e6,cGtoKg=1e-3,cCmToM=1e-2,ckToEv=.00008617328
    real,parameter ::cElectronChargeCGS=1.602176487e-20
    real :: dt,uTherm
    real :: Temp,Tpar,Tperp,Hpar,Hperp
    !--------------------------------------------------------------------------
!    nTime=10
!    nTime=1000
    nTime=100
    DtSavePlot=10
    nAltIn=10
    AltMin=8000.0e5
    AltMax=8200.0e5
    TypeGrid='Uniform'
    iAltBC = 0
    iSpecies=1
    Density=5e5
    uBulk=0.0

    Temp =1000.0
        
    write(*,*) 'init_particle'
    Time=0.0
    
    call timing_start('init_particle')
    call init_particle(nAltIn,AltMin,AltMax,TypeGrid)
    call timing_stop('init_particle')

    !set global line info
    allocate(iLineGlobal_I(nLine))
    iLineGlobal_I(1)=1

    write(*,*) 'initialize the ghost cell'
    call timing_start('sample_maxwellian_cell_boxmuller')
    call sample_maxwellian_cell_boxmuller(iAltBC,iSpecies,Density,uBulk,Temp)
    call timing_stop('sample_maxwellian_cell_boxmuller')
    
    write(*,*)'nParticle=',nParticle
    !Particles_I%vperp=Particles_I%vperp*2.0
    
    !Sort the particles
    call timing_start('sort_particles')
    call sort_particles
    call timing_stop('sort_particles')

    !plot initial distirbution before collisions
    call plot_distribution_cell(iSpecies,iAltBC)

    call calc_moments_cell(iSpecies,iAltBC,&
         densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempTmp,&
         Tpar,Tperp,Hpar,Hperp)

    uTherm=sqrt(8.0*cBoltzmannCGS*TempTmp/Mass_I(iSpecies)/cPi)
    write(*,*) 0.01*uTherm**2/(9.55e2*((6375e5+AltMin)/6375.0e5)**13.3)

    !push the guiding center 100 times and reinitialize ghost cell each time
    do iTime=1,nTime
       call timing_start('apply_wpi')
       call apply_wave_particle_interaction(6375.0e5)
       call timing_stop('apply_wpi')
       
       !advance the time
       Time=Time+DtMove
       
       !Sort the particles
       call timing_start('sort_particles')
       call sort_particles
       call timing_stop('sort_particles')

       call calc_moments_cell(iSpecies,iAltBC,&
            densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,TempTmp,&
            Tpar,Tperp,Hpar,Hperp)
       write(*,*) densityTmp,uBulkParTmp,TempTmp,Tpar,Tperp

       if (floor((Time+1.0e-5)/DtSavePlot) &
            /=floor((Time+1.0e-5-DtMove)/DtSavePlot) )then 
          call plot_distribution_cell(iSpecies,iAltBC)
       endif
    enddo

  end subroutine test_wpi

  !============================================================================
  ! unit test subroutine for sampling
  subroutine test_combine_fluid_particle
    integer :: nAltIn, iCell, iSpecies,iCounter
    real :: AltMin, AltMax, Density, uBulk, Temperature
    real :: DensityFluid, uBulkFluid, TemperatureFluid
    character(len=100):: TypeGrid
    real :: densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,&
                     TempTmp,TparTmp,TperpTmp,Hpar,Hperp
    !--------------------------------------------------------------------------
    !set global line info
    allocate(iLineGlobal_I(nLine))
    iLineGlobal_I(1)=1
    iLineCurrent=1
    nParticlePerCell_I(1)=5000

    nAltIn=2
    AltMin=1000.0e5
    AltMax=1020.0e5
    TypeGrid='Uniform'
    iCell = 1
    iSpecies=1
    Density=1e5
    uBulk=1.0e5
    Temperature=1000.0


    DensityFluid=1e5
    uBulkFluid=1.0e5
    TemperatureFluid=1000.0

    write(*,*) 'init_particle'
    call init_particle(nAltIn,AltMin,AltMax,TypeGrid)

    !set the initial particle solution
    call sample_maxwellian_cell_boxmuller(iCell,iSpecies,Density,uBulk,Temperature)

    !sort the particles
    call sort_particles

    call plot_distribution_cell(iSpecies,iCell)    
    !stop


    do iCounter=1,10
       !combine with fluid
       call combine_fluid_particle_cell(iCell,iSpecies,DensityFluid,uBulkFluid,&
            TemperatureFluid,0.5)
       
       !sort the particles
       call sort_particles

       call split_join_particles(0.05)       

       call sort_particles

       !now plot from sorted
       call plot_distribution_cell(iSpecies,iCell)
       
       call calc_moments_cell(iSpecies,iCell,&
            densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,&
            TempTmp,TparTmp,TperpTmp,Hpar,Hperp)    
       !write(*,*) 'combined moments',densityTmp,uBulkParTmp,uBulkPerpTmp,PressureTmp,&
        !    TempTmp,TparTmp,TperpTmp

       write(*,*) 'nParticle',nParticle
    end do
       !    call plot_profile
  end subroutine test_combine_fluid_particle


end Module ModParticle
