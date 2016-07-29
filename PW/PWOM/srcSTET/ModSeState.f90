Module ModSeState
  implicit none
  
  private !except
  
  real, public :: Time
  !\
  ! FLUX VARIABLES; cm-2 s-1 eV-1 sr-1
  !/
  
  ! Flux in the plasmasphere from first to second (up) ionosphere 
  real, public,allocatable :: phiup(:,:,:,:),phidn(:,:,:,:)
  
  ! Flux in the ionosphere from first to second (up) ionosphere 
  real, public,allocatable :: iphiup(:,:,:,:),iphidn(:,:,:,:)
  
  ! Flux from previous time step
  real, public,allocatable :: lphiup(:,:,:,:),lphidn(:,:,:,:)
  real, public,allocatable :: liphiup(:,:,:,:),liphidn(:,:,:,:)
  
  !Omnidirectional flux for upward/downward-directed hemisphere
  real, public,allocatable :: specup(:,:,:),specdn(:,:,:)
  
  real, public :: delt,epsilon
  
  ! Ring current parameters for coulomb collisions (note that ring.dat would 
  ! need to be read in the future to use this feature.
  integer,parameter :: IRC=0 
  real :: NRC(IRC),TRC(IRC),MRC(IRC)          !,Iterm1,Iterm2

  real,allocatable :: SRC(:)
  integer,parameter :: NRing=6 ! ring current Maxwellian fits

  ! ilocal defines the number of points at the bottom of each ionosphere where 
  ! transport is explicitly turned off. default is shown
  integer :: ilocal = 12

  ! Arrays for electron production in ionosphere. These are only used to store 
  ! these values for output. They are already added into Qstar for the purposes 
  ! of calculation.
  real,allocatable :: Qestar_IICI(:,:,:,:) ! secondary production in ionosphere 
  real,allocatable :: Qpstar_ICI(:,:,:) ! Beam-induced production 
                                        ! (non-degredating beam)

  real,allocatable :: Qstar_ICI(:,:,:) ! Total e production spec in ionosphere 
  
  !Arrays for integrated output variables
  real,public,allocatable :: HeatingRate_IC(:,:)! volume heating rate [eV/cm3/s]
  real,public,allocatable :: NumberDens_IC(:,:) ! number density of SE [/cm3]
  real,public,allocatable :: NumberFlux_IC(:,:) ! number flux of SE [/cm2/s]
  real,public,allocatable :: SecondaryIonRate_IIC(:,:,:) 
  
  ! Total electron production rate (matches ionization rate) [/cm3/s]
  real,public,allocatable :: TotalIonizationRate_IC(:,:) 

  !Precipitation info
  logical, public :: UsePrecipitation = .false.
  real   , public :: PrecipEmin, PrecipEmax,PrecipEmean,PrecipEflux
  !Polar Rain info
  logical, public :: UsePolarRain = .false.
  real   , public :: PolarRainEmin=80.0,  PolarRainEmax=300.0, &
       PolarRainEmean=100.0, PolarRainEflux=1.0e-2
  !Ovation Precipitation info
  logical, public :: UseOvation
  real   , public :: OvationEmin,OvationEmax
  real   , public :: EMeanDiff,EFluxDiff,EMeanWave,EFluxWave,EMeanMono,EFluxMono

  !public methods
  public :: allocate_state_arrays
  public :: check_time,check_time_pot
  public :: initplas, initplas_pot,initiono,initiono_pot
  public :: update_se_state,update_se_state_pot
  public :: update_se_state_iono,update_se_state_iono_pot
  public :: calc_integrated_output
contains
  SUBROUTINE update_se_state(iLine,eThermalDensity_C,eThermalTemp_C,IsOpen)
    !*  This subroutine used to be the main program, until other operators
    !*  added (3/23/95). Now it is one of several operators called by the
    !*  driver program, 'stet1.f' (Super_Thermal_electron_Transport_Model).
    
    USE ModSeGrid, only: nEnergy, EnergyGrid_I, DeltaE_I, EnergyMin, EnergyMax,&
         EqAngleGrid_IG,nThetaAlt_II,dThetaEnd_II,mu_III,&
         FieldLineGrid_IC,nTop,Bfield_IC, BFieldEq_I, &
         nIono, nPlas, nPoint,nAngle
    
    use ModMath, only: midpnt_int
    
    use ModNumConst,    ONLY: cPi

    IMPLICIT NONE
!    INCLUDE 'numbers.h'
    integer, intent(in) :: iLine
    real   , intent(in) :: eThermalDensity_C(nPoint)
    real   , intent(in) :: eThermalTemp_C(nPoint)
    logical, intent(in) :: IsOpen
    REAL h,p,Fcheck, sigmaO,muO,velt,coef, lbeta(nPoint),alpha(nAngle), &
         sigma(nAngle),Flastj,kk, del1,del2, newphi,beta(nPoint),sigO1, &
         Fsum, delE,Flasti, Flastt,Qstar, Theta,s1,s2
    INTEGER i,j,k,iPlas,jj,Ist, &
         flag,count,warning,SPick,space,ii

    !set maximum iterations
    integer, parameter :: countmax=1200

    real :: cascade, lossum
    
    !Altitude range variables
    integer :: nAltMin, nAltMax
    
    real,parameter :: cEVtoCMperS = 5.88e7 ! convert energy to velocity

    !---------------------------------------------------------------------------
    
    !initialize lbeta to 0
    lbeta(:)=0.0
    
    !      initialize warning to 0
    warning =0
    alpha(1)=1.
    sigma(1)=0.
    kk=1.
    !  Start the energy loop
    !      PRINT *, 'MainPlas, t=',t
    ENERGY: DO j=nEnergy,1,-1
       delE=.5*(DeltaE_I(j+1)+DeltaE_I(j))
       IF (j.EQ.nEnergy) delE=DeltaE_I(j)
       velt=cEVtoCMperS*SQRT(EnergyGrid_I(j))*delt
       flag=1
       count=0
       !  Start the iteration loop
       ITERATION: DO WHILE ((flag.EQ.1).AND.(count.LT.countmax))
          flag=0
          count=count+1
          IF (count.GT.countmax-2) WRITE (10,*) 'Count: ',count
          !*  Upward Region #2: the plasmasphere
          do k=0,nThetaAlt_II(iLine,nIono)
             phiup(iLine,k,0,j)=iphiup(iLine,k,nIono,j)
          end do
          h=FieldLineGrid_IC(iLine,nIono+1)-FieldLineGrid_IC(iLine,nIono)
          
          !check if field line is open and set bounds for up and down loop
          if (IsOpen) then
             ! open line so detach hemispheres
             nAltMin = nIono+1
             nAltMax = nTop
          else
             ! closed line so keep hemispheres attached
             nAltMin = nIono+1
             nAltMax = nIono+nPlas
          endif
          
          FIELDLINE_UPWARD: DO i=nAltMin,nAltMax
!             write(*,*) 'Start FieldLine up'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_C(i),eThermalDensity_C(i),EnergyGrid_I(j))

             !KLUDGE kill pitchangle scattering
             !sigO1=0.0
             !END KLUDGE
             iPlas=i-nIono
             coef=BFieldEq_I(iLine)/Bfield_IC(iLine,i)
             CALL IonoVar3(Qstar, cascade, lossum)
             do k=1,nThetaAlt_II(iLine,i)-1
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,EqAngleGrid_IG(iLine,0:nAngle))
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=dThetaEnd_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=phiup(iLine,k,iPlas-1,j)
                Flastj=phiup(iLine,k,iPlas,j+1)
                Flastt=lphiup(iLine,k,iPlas,j)
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     mu_III(iLine,k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2,&
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,1,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,1,*9999)
             end do
             IF (i.EQ.nTop) THEN
                p=-2.*s1/dThetaEnd_II(iLine,i)
             ELSE
                p=-s1/SQRT(ABS(1./coef-1.))
             END IF
             Flastt=lphiup(iLine,nThetaAlt_II(iLine,i),iPlas,j)
             Fsum=phiup(iLine,nThetaAlt_II(iLine,i)-1,iPlas,j)+phidn(iLine,nThetaAlt_II(iLine,i)-1,iPlas,j)
             Fcheck=phiup(iLine,nThetaAlt_II(iLine,i),iPlas,j)
             Flastj=phiup(iLine,nThetaAlt_II(iLine,i),iPlas,j+1)
             newphi=(Flastt/velt-p*Fsum/(2.*dThetaEnd_II(iLine,i)) &
                  +s2*(lbeta(i)* Flastj/delE+Qstar+cascade))/ &
                  (1/velt-p/dThetaEnd_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,2,*9999)
             
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             phiup(iLine,nThetaAlt_II(iLine,i),iPlas,j)=newphi
             do k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*phiup(iLine,k,iPlas,j)+sigma(k)
                Fcheck=phiup(iLine,k-1,iPlas,j)
                !            CALL CheckWarn(newphi,warning,3,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                phiup(iLine,k-1,iPlas,j)=newphi
             end do
             CALL midpnt_int(specup(iLine,j,i),phiup(iLine,0,iPlas,j),&
                  mu_III(iLine,0,i),1, nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specup(iLine,j,i)=-.5*specup(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) time,j,i,iPlas,nThetaAlt_II(iLine,i)/4,&
                     flag,(phiup(iLine,k,iPlas,j),k=0,nThetaAlt_II(iLine,i))
                WRITE (10,*) time,j,i,iPlas,nThetaAlt_II(iLine,i),flag,&
                     (phiup(iLine,k,iPlas,j),k=nThetaAlt_II(iLine,i)/2+1, &
                     nThetaAlt_II(iLine,i),2)
             END IF
             h=FieldLineGrid_IC(iLine,i+1)-FieldLineGrid_IC(iLine,i)            
             
! End Upward Region #2 loop
          end DO FIELDLINE_UPWARD
          !*  Downward Region #2: the plasmasphere

          !check if field line is open and set bounds for up and down loop
          if (IsOpen) then
             ! open line so detach hemispheres
             nAltMin = nIono+1
             nAltMax = nTop
             !set step
             h = &
                  FieldLineGrid_IC(iLine,nTop) &
                  -FieldLineGrid_IC(iLine,nTop+1)
             !fill bc from iono
             do k=0,nThetaAlt_II(iLine,nTop)
                phidn(iLine,k,nTop+1,j) = 0.0
             end do
          else
             ! closed line so keep hemispheres attached
             nAltMin = nIono+1
             nAltMax = nIono+nPlas
             !set step
             h = &
                  FieldLineGrid_IC(iLine,nIono+nPlas) &
                  -FieldLineGrid_IC(iLine,nIono+nPlas+1)
             !fill bc from iono
             do k=0,nThetaAlt_II(iLine,nIono)
                phidn(iLine,k,nPlas+1,j)=iphidn(iLine,k,nIono+1,j)
             end do
          endif
          
          FIELDLINE_DOWN: DO i=nAltMax,nAltMin,-1
!             write(*,*) 'Start FieldLine Down'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_C(i),eThermalDensity_C(i),EnergyGrid_I(j))
             !KLUDGE kill pitchangle scattering
             !sigO1=0.0
             !END KLUDGE
             iPlas=i-nIono
             coef=BFieldEq_I(iLine)/Bfield_IC(iLine,i)
             CALL IonoVar3(Qstar, cascade, lossum)
             DO k=1,nThetaAlt_II(iLine,i)-1
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,EqAngleGrid_IG(iLine,0:nAngle))
                Theta=cPi-Theta
                del1=-del1
                del2=-del2
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=-dThetaEnd_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=phidn(iLine,k,iPlas+1,j)
                Flastj=phidn(iLine,k,iPlas,j+1)
                Flastt=lphidn(iLine,k,iPlas,j)
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     -mu_III(iLine,k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2, &
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,4,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,4,*9999)
             enddo
             IF (i.EQ.nTop) THEN
                p=2.*s1/dThetaEnd_II(iLine,i)
             ELSE
                p=s1/SQRT(ABS(1./coef-1.))
             END IF
             Flastt=phidn(iLine,nThetaAlt_II(iLine,i),iPlas,j)
             Fsum=phidn(iLine,nThetaAlt_II(iLine,i)-1,iPlas,j)+phiup(iLine,nThetaAlt_II(iLine,i)-1,iPlas,j)
             Fcheck=phidn(iLine,nThetaAlt_II(iLine,i),iPlas,j)
             Flastj=phidn(iLine,nThetaAlt_II(iLine,i),iPlas,j+1)
             newphi=(Flastt/velt+p*Fsum/(2.*dThetaEnd_II(iLine,i))+s2*(lbeta(i)* &
                  Flastj/delE+Qstar+cascade))&
                  /(1/velt+p/dThetaEnd_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,5,*9999)
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             phidn(iLine,nThetaAlt_II(iLine,i),iPlas,j)=newphi
             DO k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*phidn(iLine,k,iPlas,j)+sigma(k)
                Fcheck=phidn(iLine,k-1,iPlas,j)
                !            CALL CheckWarn(newphi,warning,6,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                phidn(iLine,k-1,iPlas,j)=newphi
             end DO
             CALL midpnt_int(specdn(iLine,j,i),phidn(iLine,0,iPlas,j),&
                  mu_III(iLine,0,i),1, nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specdn(iLine,j,i)=.5*specdn(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) time,j,i,iPlas,flag,&
                     (phidn(iLine,k,iPlas,j),k=0,nThetaAlt_II(iLine,i)/2,2)
                WRITE (10,*) time,j,i,iPlas,flag,&
                     (phidn(iLine,k,iPlas,j),k=nThetaAlt_II(iLine,i)/2+1, &
                     nThetaAlt_II(iLine,i),2)
             END IF
             h=FieldLineGrid_IC(iLine,i-1)-FieldLineGrid_IC(iLine,i)            ! End Downward Region #2 loop
          end DO FIELDLINE_DOWN
       END DO ITERATION                  ! End of iteration loop
       CALL CheckFlag(flag,1,warning,-1,*9999)
       !  Save the calculated values for the next time and energy steps
       do i=nIono+1,nPoint-nIono
          lbeta(i)=beta(i)
       end do
       !      PRINT *, 'Energy convergence:',j,EnergyGrid_I(j),count
    end DO ENERGY
    
9999 IF (warning.GT.0) THEN
       PRINT *, 'Negative Densities Occurred:',newphi,warning
       PRINT *, time,j,i,k,Bfield_IC(iLine,i),BFieldEq_I(iLine),&
            BFieldEq_I(iLine),FieldLineGrid_IC(iLine,i),h,Theta,muO,&
            mu_III(iLine,k,i),del1, dThetaEnd_II(iLine,i),beta(i),sigmaO,p,count,&
            Qstar
       PRINT *, phiup(iLine,k,iPlas,j),phidn(iLine,k,iPlas,j)
       STOP
    ELSE IF (warning.LT.0) THEN
       PRINT *, 'No convergence: ',count,Time,warning,j
       STOP
    END IF
    !  Format for the nonconvergent output
!!!27  FORMAT (100E10.3,4I4,1P,100E10.3)
    
    RETURN
  END SUBROUTINE update_se_state
  
  !============================================================================
  SUBROUTINE update_se_state_iono(iLine,IsIono1,eThermalDensity_C,&
       eThermalTemp_C,nNeutral,NeutralDens_IC,&
       SIGS,SIGI,SIGA,ePhotoProdSpec_IC)
    !*  Same as update_se_state but for ionosphere, includes sources
    
    USE ModSeGrid, only: nEnergy, EnergyGrid_I, DeltaE_I, EnergyMin, EnergyMax,&
         EqAngleGrid_IG,nThetaAlt_II,dThetaEnd_II,mu_III,&
         FieldLineGrid_IC,nTop,Bfield_IC, BFieldEq_I, &
         nIono, nPlas, nPoint,nAngle
    
    use ModMath, only: midpnt_int
    
    use ModNumConst,    ONLY: cPi

    IMPLICIT NONE
!    INCLUDE 'numbers.h'
    integer, intent(in) :: iLine
    logical, intent(in) :: IsIono1
    
    !thermal density and temperature
    real   , intent(in) :: eThermalDensity_C(nPoint)
    real   , intent(in) :: eThermalTemp_C(nPoint)

    !neutral atmosphere inputs
    integer, intent(in) :: nNeutral
    real   , intent(in) :: NeutralDens_IC(nNeutral,nIono)

    !incomming crossections
    real   , intent(in) :: SIGS(nNeutral,nEnergy)
    real   , intent(in) :: SIGI(nNeutral,nEnergy,nEnergy)
    real   , intent(in) :: SIGA(nNeutral,nEnergy,nEnergy)
    
    ! incomming photo electron production spectrum
    real   , intent(in) :: ePhotoProdSpec_IC(nEnergy,nIono)
    
    REAL h,p,Fcheck, sigmaO,muO,velt,coef, lbeta(nPoint),alpha(nAngle), &
         sigma(nAngle),Flastj,kk, del1,del2, newphi,beta(nPoint),sigO1, &
         Fsum, delE,Flasti, Flastt,Qstar, Theta,s1,s2
    INTEGER i,j,k,jj,Ist, &
         flag,count,warning,SPick,space,ii
    !inidicies that store full ionosphere grid(2*nIono points) and half 
    ! ionosphere grid (nIono points). These convert position along the field 
    ! line into these sub grids.
    integer :: iIono, iIonoHalf

    !set maximum iterations
    integer, parameter :: countmax=1200

    real :: cascade, lossum
    
    !Altitude range variables
    integer :: nAltMin, nAltMax
    
    ! array to hold temporarily the values being overwritten by BCs
    real :: fhold(0:nAngle)
    real,parameter :: cEVtoCMperS = 5.88e7 ! convert energy to velocity

    !---------------------------------------------------------------------------

    !initialize lbeta to 0
    lbeta(:)=0.0
    
    !      initialize warning to 0
    warning =0
    alpha(1)=1.
    sigma(1)=0.
    kk=1.
    !  Start the energy loop
    !      PRINT *, 'MainPlas, t=',t
    ENERGY: DO j=nEnergy,1,-1
       delE=.5*(DeltaE_I(j+1)+DeltaE_I(j))
       IF (j.EQ.nEnergy) delE=DeltaE_I(j)
       velt=cEVtoCMperS*SQRT(EnergyGrid_I(j))*delt
       flag=1
       count=0
       !  Start the iteration loop
       ITERATION: DO WHILE ((flag.EQ.1).AND.(count.LT.countmax))
          flag=0
          count=count+1
          IF (count.GT.countmax-2) WRITE (10,*) 'Count: ',count

          ! Upward Region

          ! set BC when in iono2. Put BC in iphiup at top of iono1 and hold 
          ! the current value to put it back after the calculation
          if (.not.IsIono1)then
             do k=0,nThetaAlt_II(iLine,nIono)
                Fhold(k)=iphiup(iLine,k,nIono,j)
             end do

             ! fill BC from the the plasmasphere
             do k=0,nThetaAlt_II(iLine,nIono)
                iphiup(iLine,k,nIono,j)=phiup(iLine,k,nPlas,j)
             end do

             ! Add precip info here
          endif
             
          !set bounds for ionosphere loop depending if you are in iono 1 or 2
          if (IsIono1) then
             ! set upward bounds for iono1
             nAltMin = 1
             nAltMax = nIono
             h=FieldLineGrid_IC(iLine,nAltMin+1)-FieldLineGrid_IC(iLine,nAltMin)
          else
             ! set upward bounds for iono2
             nAltMin = nPoint-nIono+1
             nAltMax = nPoint
             h=FieldLineGrid_IC(iLine,nAltMin)-FieldLineGrid_IC(iLine,nAltMin-1)
          endif


          FIELDLINE_UPWARD: DO i=nAltMin,nAltMax
!             write(*,*) 'Start FieldLine up'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_C(i),eThermalDensity_C(i),EnergyGrid_I(j))

             ! set iIono and iIonoHalf indices based on 1st or 2nd ionosphere
             ! remember, iIono is index spanning both ionospheres, 
             ! iIonoHalf is index spanning just one ionosphere
             if(IsIono1) then
                iIono=i
                iIonoHalf=i
             else
                iIono=i-nPlas
                iIonoHalf = nPoint-i+1
             end if
             
             !set the b-field ratio coeficient
             coef=BFieldEq_I(iLine)/Bfield_IC(iLine,i)
             
             ! get sigma0 and electron production (primary+secondary)
             CALL get_sigma0_and_eprod(j,nNeutral,sigO1,&
                  ePhotoProdSpec_IC(j,iIonoHalf), &
                  NeutralDens_IC(:,iIonoHalf),SIGS,SIGI,SIGA,Qstar, &
                  specup(iLine,:,i),specdn(iLine,:,i),&
                  Qestar_IICI(iLine,:,iIono,j),Qpstar_ICI(iLine,iIono,j))
             Qstar_ICI(iLine,iIono,j)=Qstar
                       
             do k=1,nThetaAlt_II(iLine,i)-1
                call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                     NeutralDens_IC(:,iIonoHalf),j,iIono,k,cascade,lossum,1)
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,&
                     EqAngleGrid_IG(iLine,0:nAngle))
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=dThetaEnd_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=iphiup(iLine, k,iIono-1,j)
                Flastj=iphiup(iLine, k,iIono,j+1)
                Flastt=liphiup(iLine,k,iIono,j)
                
                !When we are in the ilocal region kill the transport terms in 
                ! coeficients by setting kk to 0
                kk=1
                if (iIonoHalf <= iLocal) kk=0
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     mu_III(iLine,k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2,&
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,1,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,1,*9999)
                
             end do

             p=-s1/SQRT(ABS(1./coef-1.))

             Flastt=liphiup(iLine,nThetaAlt_II(iLine,i),iIono,j)
             Fsum=iphiup(iLine,nThetaAlt_II(iLine,i)-1,iIono,j)&
                  +iphidn(iLine,nThetaAlt_II(iLine,i)-1,iIono,j)
             Fcheck=iphiup(iLine,nThetaAlt_II(iLine,i),iIono,j)
             Flastj=iphiup(iLine,nThetaAlt_II(iLine,i),iIono,j+1)
             call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                  NeutralDens_IC(:,iIonoHalf),j,iIono,nThetaAlt_II(iLine,i),&
                  cascade,lossum,1)
             newphi=(Flastt/velt-p*Fsum/(2.*dThetaEnd_II(iLine,i)) &
                  +s2*(lbeta(i)* Flastj/delE+Qstar+cascade))/ &
                  (1/velt-p/dThetaEnd_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,2,*9999)
             
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             iphiup(iLine,nThetaAlt_II(iLine,i),iIono,j)=newphi
             
             do k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*iphiup(iLine,k,iIono,j)+sigma(k)
                Fcheck=iphiup(iLine,k-1,iIono,j)
                !            CALL CheckWarn(newphi,warning,3,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                
                
                iphiup(iLine,k-1,iIono,j)=newphi
             end do
             CALL midpnt_int(specup(iLine,j,i),iphiup(iLine,:,iIono,j),&
                  mu_III(iLine,:,i),1, nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specup(iLine,j,i)=-.5*specup(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) 'Upward region'
                WRITE (10,*) time,j,i,iIono,nThetaAlt_II(iLine,i),&
                     flag,(iphiup(iLine,k,iIono,j),k=0,nThetaAlt_II(iLine,i))
             END IF
             
             if (i < nPoint) then
                h=FieldLineGrid_IC(iLine,i+1)-FieldLineGrid_IC(iLine,i)
             else
                h=FieldLineGrid_IC(iLine,i)-FieldLineGrid_IC(iLine,i-1)
             endif
             ! End Upward Region of ionosphere
          end DO FIELDLINE_UPWARD
          
          !put back the iphiup at top of iono1 (which was storing plasmasph bc)
          if(.not.IsIono1)then
             do k=0,nThetaAlt_II(iLine,nIono)
                iphiup(iLine,k,nIono,j)=Fhold(k)
             end do
          endif
          
          !*  Downward Region of the ionosphere

          ! set BC when in iono1. Put BC in iphidn at top of iono2 and hold 
          ! the current value to put it back after the calculation
          if (IsIono1)then
             do k=0,nThetaAlt_II(iLine,nIono)
                Fhold(k)=iphidn(iLine,k,nIono+1,j)
             end do

             ! fill BC from the the plasmasphere
             do k=0,nThetaAlt_II(iLine,nIono)
                iphidn(iLine,k,nIono+1,j)=phidn(iLine,k,1,j)
             end do


          endif

          !set bounds for ionosphere loop depending if you are in iono 1 or 2
          if (IsIono1) then
             ! set upward bounds for iono1
             nAltMin = 1
             nAltMax = nIono
             h=FieldLineGrid_IC(iLine,nAltMax)-FieldLineGrid_IC(iLine,nAltMax+1)
          else
             ! set upward bounds for iono2
             nAltMin = nPoint-nIono+1
             nAltMax = nPoint
             h=FieldLineGrid_IC(iLine,nAltMax-1)-FieldLineGrid_IC(iLine,nAltMax)
          endif
          
          FIELDLINE_DOWN: DO i=nAltMax,nAltMin,-1
!             write(*,*) 'Start FieldLine Down'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_C(i),eThermalDensity_C(i),EnergyGrid_I(j))

             ! set iIono and iIonoHalf indices based on 1st or 2nd ionosphere
             ! remember, iIono is index spanning both ionospheres, 
             ! iIonoHalf is index spanning just one ionosphere
             if(IsIono1) then
                iIono=i
                iIonoHalf=i
             else
                iIono=i-nPlas
                iIonoHalf = nPoint-i+1
             end if

             !set the b-field ratio coeficient
             coef=BFieldEq_I(iLine)/Bfield_IC(iLine,i)

             ! get sigma0 and electron production (primary+secondary)
             CALL get_sigma0_and_eprod(j,nNeutral,sigO1,&
                  ePhotoProdSpec_IC(j,iIonoHalf), &
                  NeutralDens_IC(:,iIonoHalf),SIGS,SIGI,SIGA,Qstar, &
                  specup(iLine,:,i),specdn(iLine,:,i),&
                  Qestar_IICI(iLine,:,iIono,j),Qpstar_ICI(iLine,iIono,j))
             Qstar_ICI(iLine,iIono,j)=Qstar
             
             DO k=1,nThetaAlt_II(iLine,i)-1
                call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                     NeutralDens_IC(:,iIonoHalf),j,iIono,k,cascade,lossum,2)
                CALL ThetaVar1(nAngle,Theta,del1,del2,k, &
                     EqAngleGrid_IG(iLine,0:nAngle))
                Theta=cPi-Theta
                del1=-del1
                del2=-del2
                IF (k.EQ.nThetaAlt_II(iLine,i)-1) del2=-dThetaEnd_II(iLine,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=iphidn(iLine,k,iIono+1,j)
                Flastj=iphidn(iLine,k,iIono,j+1)
                Flastt=liphidn(iLine,k,iIono,j)
                
                !When we are in the ilocal region kill the transport terms in 
                ! coeficients by setting kk to 0
                kk=1
                if (iIonoHalf <= iLocal) kk=0
                CALL NumCalcVar(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     -mu_III(iLine,k,i),coef,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2, &
                     cascade,lossum)
                !            CALL CheckWarn(sigma(k+1),warning,4,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,4,*9999)
             enddo

             p=s1/SQRT(ABS(1./coef-1.))

             Flastt=liphidn(iLine,nThetaAlt_II(iLine,i),iIono,j)
             Fsum=iphidn(iLine,nThetaAlt_II(iLine,i)-1,iIono,j)&
                  +iphiup(iLine,nThetaAlt_II(iLine,i)-1,iIono,j)
             Fcheck=iphidn(iLine,nThetaAlt_II(iLine,i),iIono,j)
             Flastj=iphidn(iLine,nThetaAlt_II(iLine,i),iIono,j+1)
             call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                  NeutralDens_IC(:,iIonoHalf),j,iIono,nThetaAlt_II(iLine,i), &
                  cascade,lossum,2)
             newphi=(Flastt/velt+p*Fsum/(2.*dThetaEnd_II(iLine,i))+s2*(lbeta(i)* &
                  Flastj/delE+Qstar+cascade))&
                  /(1/velt+p/dThetaEnd_II(iLine,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,5,*9999)
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             iphidn(iLine,nThetaAlt_II(iLine,i),iIono,j)=newphi
             DO k=nThetaAlt_II(iLine,i),1,-1
                newphi=alpha(k)*iphidn(iLine,k,iIono,j)+sigma(k)
                Fcheck=iphidn(iLine,k-1,iIono,j)
                !            CALL CheckWarn(newphi,warning,6,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                iphidn(iLine,k-1,iIono,j)=newphi
             end DO
             CALL midpnt_int(specdn(iLine,j,i),iphidn(iLine,:,iIono,j),&
                  mu_III(iLine,:,i),1, nThetaAlt_II(iLine,i)+1,nAngle+1,1)
             specdn(iLine,j,i)=.5*specdn(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) 'Downward region'
                WRITE (10,*) time,j,i,iIono,flag,&
                     (iphidn(iLine,k,iIono,j),k=0,nThetaAlt_II(iLine,i))
             END IF
             if (i==1) then
                h=FieldLineGrid_IC(iLine,i)-FieldLineGrid_IC(iLine,i+1)        
             else
                h=FieldLineGrid_IC(iLine,i-1)-FieldLineGrid_IC(iLine,i)         
             endif
             ! End Downward Region of ionosphere
          end DO FIELDLINE_DOWN
          !put back the iphidn at top of iono2 (which was storing plasmasph bc)
          if (IsIono1)then
             do k=0,nThetaAlt_II(iLine,nIono)
                iphidn(iLine,k,nIono+1,j)=Fhold(k)
             end do
          endif
       END DO ITERATION                  ! End of iteration loop
       CALL CheckFlag(flag,1,warning,-1,*9999)
       !  Save the calculated values for the next time and energy steps
       if (IsIono1) then
          lbeta(1:nIono)=beta(1:nIono)
       else
          lbeta(nPoint-nIono+1:nPoint)=beta(nPoint-nIono+1:nPoint)
       endif
    end DO ENERGY
    
9999 IF (warning.GT.0) THEN
       PRINT *, 'Negative Densities Occurred:',newphi,warning
       PRINT *, time,j,i,k,Bfield_IC(iLine,i),BFieldEq_I(iLine),&
            BFieldEq_I(iLine),FieldLineGrid_IC(iLine,i),h,Theta,muO,&
            mu_III(iLine,k,i),del1, dThetaEnd_II(iLine,i),beta(i),sigmaO,p,count,&
            Qstar
       PRINT *, iphiup(iLine,k,iIono,j),iphidn(iLine,k,iIono,j)
       STOP
    ELSE IF (warning.LT.0) THEN
       PRINT *, 'No convergence: ',count,Time,warning,j
       STOP
    END IF
    !  Format for the nonconvergent output
!!!27  FORMAT (100E10.3,4I4,1P,100E10.3)
    
    RETURN
  END SUBROUTINE update_se_state_iono
  
  !=============================================================================
  ! same as the update_se_state except now we include the effect of a 
  ! potential. This is the plasmasphere solution
  SUBROUTINE update_se_state_pot(iLine,eThermalDensity_C,eThermalTemp_C,IsOpen)
    
    USE ModSeGrid, only: nEnergy, EnergyGrid_I, DeltaE_I, EnergyMin, EnergyMax,&
         EqAngleGrid_IG,nThetaAlt_IIC,dThetaEnd_III,mu_IIIC,&
         FieldLineGrid_IC,nTop,Bfield_IC, MaxAlt_IC,nMu0RefAlt_II,&
         BField0_II, DeltaPot_IC, DeltaPot0_II, KineticEnergy_IIC,&
         nIono, nPlas, nPoint,nAngle,nPwRegion
    
    use ModMath, only: midpnt_int
    
    use ModNumConst,    ONLY: cPi

    IMPLICIT NONE
!    INCLUDE 'numbers.h'
    integer, intent(in) :: iLine
    real   , intent(in) :: eThermalDensity_C(nPoint)
    real   , intent(in) :: eThermalTemp_C(nPoint)
    logical, intent(in) :: IsOpen
    REAL h,p,Fcheck, sigmaO,muO,velt,eta, lbeta(nPoint),alpha(nAngle), &
         sigma(nAngle),Flastj,kk, del1,del2, newphi,beta(nPoint),sigO1, &
         Fsum, delE,Flasti, Flastt,Qstar, Theta,s1,s2
    INTEGER i,j,k,iPlas,jj,Ist, &
         flag,count,warning,SPick,space,ii

    !set maximum iterations
    integer, parameter :: countmax=2500

    real :: cascade, lossum
    
    !Altitude range variables
    integer :: nAltMin, nAltMax
  
    ! ML variables
    real :: mui, muj ! mu at i-1/2 and i+1/2
    logical, parameter :: UseMLfix=.false.

    real,parameter :: cEVtoCMperS = 5.88e7 ! convert energy to velocity

    !---------------------------------------------------------------------------
    
    !initialize lbeta to 0
    lbeta(:)=0.0
    
    !      initialize warning to 0
    warning =0
    alpha(1)=1.
    sigma(1)=0.
    kk=1.
    !  Start the energy loop
    !      PRINT *, 'MainPlas, t=',t
    ENERGY: DO j=nEnergy,1,-1
       delE=.5*(DeltaE_I(j+1)+DeltaE_I(j))
       IF (j.EQ.nEnergy) delE=DeltaE_I(j)
       flag=1
       count=0
       !  Start the iteration loop
       ITERATION: DO WHILE ((flag.EQ.1).AND.(count.LT.countmax))
          flag=0
          count=count+1
          IF (count.GT.countmax-2) WRITE (10,*) 'Count: ',count
          !*  Upward Region #2: the plasmasphere
          
          do k=0,nThetaAlt_IIC(iLine,j,nIono)
             phiup(iLine,k,0,j)=iphiup(iLine,k,nIono,j)
          end do
                    
          h=FieldLineGrid_IC(iLine,nIono+1)-FieldLineGrid_IC(iLine,nIono)
          
          !check if field line is open and set bounds for up and down loop
          if (IsOpen) then
             ! open line so detach hemispheres
             nAltMin = nIono+1
             !nAltMax = MaxAlt_IC(iLine,j)
             nAltMax = min(MaxAlt_IC(iLine,j),nIono+nPwRegion)
          else
             ! closed line so keep hemispheres attached. Note that should 
             ! they get detached for energies less than the potential we will
             ! cycle over those altitude steps
             nAltMin = nIono+1
             nAltMax = nIono+nPlas
          endif
          
          FIELDLINE_UPWARD: DO i=nAltMin,nAltMax
             !write(*,*) 'start upward',i,j
             
             !when MaxAlt_IC is not equal to nTop, then the hemispheres are 
             ! detached for these energies between MaxAlt_IC(iLine,j)+1 and 
             ! nPoint-MaxAlt_IC. note that we need to lower the side on the 
             ! second hemisphere by one so we can take the BC from the phidn
             ! so we have reflection
             if ((MaxAlt_IC(iLine,j) /= nTop) .and. (i > MaxAlt_IC(iLine,j)) &
                  .and. (i <= nPoint-MaxAlt_IC(iLine,j))) then
                !set the reflection bc for hemisphere 2
                do k=0,nThetaAlt_IIC(iLine,j,nPoint-MaxAlt_IC(iLine,j))
                   phiup(iLine,k,nPoint-MaxAlt_IC(iLine,j)-nIono,j)=&
                        phidn(iLine,k,nPoint-MaxAlt_IC(iLine,j)-nIono,j)
                end do
                cycle FIELDLINE_UPWARD
             endif
             
             !velt is now altitude dependent
             velt=cEVtoCMperS*SQRT(KineticEnergy_IIC(iLine,j,i))*delt

!             write(*,*) 'Start FieldLine up'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_C(i),eThermalDensity_C(i),KineticEnergy_IIC(iLine,j,i))
             
             !KLUDGE kill pitchangle scattering
             !sigO1=0.0
             !END KLUDGE
             iPlas=i-nIono
             eta = (Bfield0_II(iLine,j)/Bfield_IC(iLine,i)) &
                  *(EnergyGrid_I(j)-DeltaPot_IC(iLine,i))   &
                  /(EnergyGrid_I(j)-DeltaPot0_II(iLine,j))

             CALL IonoVar3(Qstar, cascade, lossum)
             do k=1,nThetaAlt_IIC(iLine,j,i)-1
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,EqAngleGrid_IG(iLine,0:nAngle))
                IF (k.EQ.nThetaAlt_IIC(iLine,j,i)-1) del2=dThetaEnd_III(iLine,j,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=phiup(iLine,k,iPlas-1,j)
                
                !if (UseMLfix) then
                if(k>nThetaAlt_IIC(iLine,j,i-1)) Flasti=phidn(iLine,k,iPlas,j)
                !endif
                
                Flastj=phiup(iLine,k,iPlas,j+1)
                Flastt=lphiup(iLine,k,iPlas,j)
                !write(*,*) 'test1: i,j,k,ko',i,j,k,nThetaAlt_IIC(iLine,j,i)
                CALL NumCalcVar_pot(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     mu_IIIC(iLine,k,j,i),eta,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2,&
                     cascade,lossum,EnergyGrid_I(j),DeltaPot0_II(iLine,j), &
                     KineticEnergy_IIC(iLine,j,i),eThermalDensity_C(i))
                !            CALL CheckWarn(sigma(k+1),warning,1,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,1,*9999)
                
             end do

             ! Note that the reference altitude can occur in two places when 
             ! a potential is included
             IF (i==nMu0RefAlt_II(iLine,j) .or. &
                  i == nPoint-nMu0RefAlt_II(iLine,j)) THEN
                p=-2.*s1/dThetaEnd_III(iLine,j,i)
             ELSE
                p=pend_pot(kk,s1,eta,EnergyGrid_I(j),DeltaPot0_II(iLine,j), &
                     KineticEnergy_IIC(iLine,j,i),eThermalDensity_C(i))
             END IF
             
             Flastt=lphiup(iLine,nThetaAlt_IIC(iLine,j,i),iPlas,j)
             Fsum=phiup(iLine,nThetaAlt_IIC(iLine,j,i)-1,iPlas,j)+phidn(iLine,nThetaAlt_IIC(iLine,j,i)-1,iPlas,j)
             Fcheck=phiup(iLine,nThetaAlt_IIC(iLine,j,i),iPlas,j)
             Flastj=phiup(iLine,nThetaAlt_IIC(iLine,j,i),iPlas,j+1)
             
             if (UseMLfix) then
                Flasti=phiup(iLine,nThetaAlt_IIC(iLine,j,i),iPlas-1,j)
                mui=0.5*(mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i) &
                     +mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i-1))
                muj=0.5*(mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i) &
                     +mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i+1))
                if(mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i+1)==0) muj=mui
                if(mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i-1)==0) then
                   Flasti=phidn(iLine,nThetaAlt_IIC(iLine,j,i),iPlas+1,j)
                   mui=muj
                endif
                newphi=(Flastt/velt-p*Fsum/(2.*dThetaEnd_III(iLine,j,i)) &
                     +s2*(lbeta(i)* Flastj/delE+mui*Flasti/h+Qstar+cascade))/ &
                     (1/velt-p/dThetaEnd_III(iLine,j,i)+s2*(beta(i)/delE+muj/h+lossum))       
                !newphi=max(newphi,0.0)         
             else
             newphi=(Flastt/velt-p*Fsum/(2.*dThetaEnd_III(iLine,j,i)) &
                  +s2*(lbeta(i)* Flastj/delE+Qstar+cascade))/ &
                  (1/velt-p/dThetaEnd_III(iLine,j,i)+s2*(beta(i)/delE+lossum))
             newphi=max(newphi,0.0)             
          endif
             !            CALL CheckWarn(newphi,warning,2,*9999)
             

             CALL CheckConv(newphi,Fcheck,epsilon,flag)

             phiup(iLine,nThetaAlt_IIC(iLine,j,i),iPlas,j)=newphi
             do k=nThetaAlt_IIC(iLine,j,i),1,-1
                newphi=alpha(k)*phiup(iLine,k,iPlas,j)+sigma(k)

                Fcheck=phiup(iLine,k-1,iPlas,j)
                !            CALL CheckWarn(newphi,warning,3,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)

                phiup(iLine,k-1,iPlas,j)=newphi
             end do

             CALL midpnt_int(specup(iLine,j,i),phiup(iLine,0:nAngle,iPlas,j),&
                  mu_IIIC(iLine,0:nAngle,j,i),1, nThetaAlt_IIC(iLine,j,i)+1,nAngle+1,1)
             specup(iLine,j,i)=-.5*specup(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) time,j,i,iPlas,nThetaAlt_IIC(iLine,j,i)/4,&
                     flag,(phiup(iLine,k,iPlas,j),k=0,nThetaAlt_IIC(iLine,j,i))
                WRITE (10,*) time,j,i,iPlas,nThetaAlt_IIC(iLine,j,i),flag,&
                     (phiup(iLine,k,iPlas,j),k=nThetaAlt_IIC(iLine,j,i)/2+1, &
                     nThetaAlt_IIC(iLine,j,i),2)
             END IF
             h=FieldLineGrid_IC(iLine,i+1)-FieldLineGrid_IC(iLine,i)            
             
! End Upward Region #2 loop
          end DO FIELDLINE_UPWARD
          !*  Downward Region #2: the plasmasphere

          !check if field line is open and set bounds for up and down loop
          if (IsOpen) then
             ! open line so detach hemispheres
             nAltMin = nIono+1
             !nAltMax = MaxAlt_IC(iLine,j)
             nAltMax = min(MaxAlt_IC(iLine,j),nIono+nPwRegion)
             !set step
             h = &
                  FieldLineGrid_IC(iLine,nAltMax) &
                  -FieldLineGrid_IC(iLine,nAltMax+1)
             !fill bc from iono
          !   write(*,*) j,nTop, nThetaAlt_IIC(iLine,j,nTop)
             phidn(iLine,:,nTop+1,j) = 0.0
          else
             ! closed line so keep hemispheres attached. Note that should 
             ! they get detached for energies less than the potential we will
             ! cycle over those altitude steps
             nAltMin = nIono+1
             nAltMax = nIono+nPlas
             !set step
             h = &
                  FieldLineGrid_IC(iLine,nIono+nPlas) &
                  -FieldLineGrid_IC(iLine,nIono+nPlas+1)
             !fill bc from iono
             do k=0,nThetaAlt_IIC(iLine,j,nIono)
                phidn(iLine,k,nPlas+1,j)=iphidn(iLine,k,nIono+1,j)
             end do
          endif
          
          FIELDLINE_DOWN: DO i=nAltMax,nAltMin,-1
             !write(*,*) 'start downward',i,j,MaxAlt_IC(iLine,j),nPoint-MaxAlt_IC(iLine,j)

             !when MaxAlt_IC is not equal to nTop, then the hemispheres are 
             ! detached for these energies between MaxAlt_IC(iLine,j)+1 and 
             ! nPoint-MaxAlt_IC. note that we need to lower the side on the 
             ! first hemisphere by one so we can take the BC from the phiup
             ! so we have reflection
             if ((MaxAlt_IC(iLine,j) /= nTop) .and. (i >= MaxAlt_IC(iLine,j)) &
                  .and. (i < nPoint-MaxAlt_IC(iLine,j))) then
                !set the reflection bc for hemisphere 2
                do k=0,nThetaAlt_IIC(iLine,j,MaxAlt_IC(iLine,j))
                   phidn(iLine,k,MaxAlt_IC(iLine,j)-nIono,j)=&
                        phiup(iLine,k,MaxAlt_IC(iLine,j)-nIono,j)
                end do
                cycle FIELDLINE_DOWN
             endif
             
             iPlas=i-nIono
             
             !velt is now altitude dependent
             velt=cEVtoCMperS*SQRT(KineticEnergy_IIC(iLine,j,i))*delt
             
!             write(*,*) 'Start FieldLine Down'
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_C(i),eThermalDensity_C(i),KineticEnergy_IIC(iLine,j,i))
             !KLUDGE kill pitchangle scattering
             !sigO1=0.0
             !END KLUDGE

             eta = (Bfield0_II(iLine,j)/Bfield_IC(iLine,i)) &
                  *(EnergyGrid_I(j)-DeltaPot_IC(iLine,i))   &
                  /(EnergyGrid_I(j)-DeltaPot0_II(iLine,j))
             CALL IonoVar3(Qstar, cascade, lossum)
             DO k=1,nThetaAlt_IIC(iLine,j,i)-1
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,EqAngleGrid_IG(iLine,0:nAngle))
                Theta=cPi-Theta
                del1=-del1
                del2=-del2
                IF (k.EQ.nThetaAlt_IIC(iLine,j,i)-1) &
                     del2=-dThetaEnd_III(iLine,j,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=phidn(iLine,k,iPlas+1,j)
                
                !if(UseMLfix) then
                   if(k>nThetaAlt_IIC(iLine,j,i+1)) then
                      !write(*,*) 'HA!'
                      Flasti=phiup(iLine,k,iPlas,j)
                   endif
                !endif
                
                Flastj=phidn(iLine,k,iPlas,j+1)
                Flastt=lphidn(iLine,k,iPlas,j)
                !write(*,*) 'test2: i,j,k',i,j,k
                CALL NumCalcVar_pot(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     -mu_IIIC(iLine,k,j,i),eta,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2, &
                     cascade,lossum,EnergyGrid_I(j),DeltaPot0_II(iLine,j), &
                     KineticEnergy_IIC(iLine,j,i),eThermalDensity_C(i))
                !            CALL CheckWarn(sigma(k+1),warning,4,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,4,*9999)
             enddo
             IF (i==nMu0RefAlt_II(iLine,j) .or. &
                  i == nPoint-nMu0RefAlt_II(iLine,j)) THEN
                p=2.*s1/dThetaEnd_III(iLine,j,i)
             ELSE
                p=-pend_pot(kk,s1,eta,EnergyGrid_I(j),DeltaPot0_II(iLine,j), &
                     KineticEnergy_IIC(iLine,j,i),eThermalDensity_C(i))
             END IF
             Flastt=phidn(iLine,nThetaAlt_IIC(iLine,j,i),iPlas,j)
             Fsum=phidn(iLine,nThetaAlt_IIC(iLine,j,i)-1,iPlas,j)+phiup(iLine,nThetaAlt_IIC(iLine,j,i)-1,iPlas,j)
             Fcheck=phidn(iLine,nThetaAlt_IIC(iLine,j,i),iPlas,j)
             Flastj=phidn(iLine,nThetaAlt_IIC(iLine,j,i),iPlas,j+1)
!             newphi=(Flastt/velt+p*Fsum/(2.*dThetaEnd_III(iLine,j,i))+s2*(lbeta(i)* &
!                  Flastj/delE+Qstar+cascade))&
!                  /(1/velt+p/dThetaEnd_III(iLine,j,i)+s2*(beta(i)/delE+lossum))
             if(UseMLfix) then
                Flasti=phidn(iLine,nThetaAlt_IIC(iLine,j,i),iPlas+1,j)
                mui=0.5*(mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i) &
                     +mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i+1))
                muj=0.5*(mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i) &
                     +mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i-1))
                if(mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i-1)==0) muj=mui
                if(mu_IIIC(iLine,nThetaAlt_IIC(iLine,j,i),j,i+1)==0) then
                   Flasti=phiup(iLine,nThetaAlt_IIC(iLine,j,i),iPlas-1,j)
                   mui=muj
                endif
                newphi=(Flastt/velt+p*Fsum/(2.*dThetaEnd_III(iLine,j,i)) &
                     +s2*(lbeta(i)* Flastj/delE-mui*Flasti/h+Qstar+cascade))/ &
                     (1/velt+p/dThetaEnd_III(iLine,j,i)+s2*(beta(i)/delE-muj/h+lossum))         
                !newphi=max(newphi,0.0)         
             else
                newphi=(Flastt/velt+p*Fsum/(2.*dThetaEnd_III(iLine,j,i)) &
                     +s2*(lbeta(i)*Flastj/delE+Qstar+cascade))&
                     /(1/velt+p/dThetaEnd_III(iLine,j,i)+s2*(beta(i)/delE &
                     +lossum))
                newphi=max(newphi,0.0)         
             end if

             !            CALL CheckWarn(newphi,warning,5,*9999)
         
             CALL CheckConv(newphi,Fcheck,epsilon,flag)

             phidn(iLine,nThetaAlt_IIC(iLine,j,i),iPlas,j)=newphi
             
             DO k=nThetaAlt_IIC(iLine,j,i),1,-1
                newphi=alpha(k)*phidn(iLine,k,iPlas,j)+sigma(k)
                Fcheck=phidn(iLine,k-1,iPlas,j)
                !            CALL CheckWarn(newphi,warning,6,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)

                phidn(iLine,k-1,iPlas,j)=newphi
             end DO
             CALL midpnt_int(specdn(iLine,j,i),phidn(iLine,0:nAngle,iPlas,j),&
                  mu_IIIC(iLine,0:nAngle,j,i),1, nThetaAlt_IIC(iLine,j,i)+1,nAngle+1,1)
             specdn(iLine,j,i)=.5*specdn(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) time,j,i,iPlas,flag,&
                     (phidn(iLine,k,iPlas,j),k=0,nThetaAlt_IIC(iLine,j,i)/2,2)
                WRITE (10,*) time,j,i,iPlas,flag,&
                     (phidn(iLine,k,iPlas,j),k=nThetaAlt_IIC(iLine,j,i)/2+1, &
                     nThetaAlt_IIC(iLine,j,i),2)
             END IF
             h=FieldLineGrid_IC(iLine,i-1)-FieldLineGrid_IC(iLine,i)            ! End Downward Region #2 loop
          end DO FIELDLINE_DOWN
       END DO ITERATION                  ! End of iteration loop
                       
       CALL CheckFlag(flag,1,warning,-1,*9999)
       !  Save the calculated values for the next time and energy steps
       do i=nIono+1,nPoint-nIono
          lbeta(i)=beta(i)
       end do
       !      PRINT *, 'Energy convergence:',j,EnergyGrid_I(j),count
    end DO ENERGY
    
9999 IF (warning.GT.0) THEN
       PRINT *, 'Negative Densities Occurred:',newphi,warning
       PRINT *, time,j,i,k,Bfield_IC(iLine,i),Bfield0_II(iLine,j),&
            FieldLineGrid_IC(iLine,i),h,Theta,muO,&
            mu_IIIC(iLine,k,j,i),del1, dThetaEnd_III(iLine,j,i),beta(i),sigmaO,p,count,&
            Qstar
       PRINT *, phiup(iLine,k,iPlas,j),phidn(iLine,k,iPlas,j)
       STOP
    ELSE IF (warning.LT.0) THEN
       PRINT *, 'No convergence: ',count,Time,warning,j
       STOP
    END IF
    !  Format for the nonconvergent output
!!!27  FORMAT (100E10.3,4I4,1P,100E10.3)
    
    RETURN
  END SUBROUTINE update_se_state_pot
  
  !============================================================================

  !============================================================================
  SUBROUTINE update_se_state_iono_pot(iLine,IsIono1,eThermalDensity_C,&
       eThermalTemp_C,nNeutral,NeutralDens_IC,&
       SIGS,SIGI,SIGA,ePhotoProdSpec_IC)
    !*  Same as update_se_state but for ionosphere, includes sources
    
    USE ModSeGrid, only: nEnergy, EnergyGrid_I, DeltaE_I, EnergyMin, EnergyMax,&
         EqAngleGrid_IG,nThetaAlt_IIC,dThetaEnd_III,mu_IIIC,&
         FieldLineGrid_IC,nTop,Bfield_IC, BField0_II, DeltaPot_IC, &
         DeltaPot0_II, KineticEnergy_IIC,&
         nIono, nPlas, nPoint,nAngle
    
    use ModMath, only: midpnt_int
    
    use ModNumConst,    ONLY: cPi

    IMPLICIT NONE
!    INCLUDE 'numbers.h'
    integer, intent(in) :: iLine
    logical, intent(in) :: IsIono1
    
    !thermal density and temperature
    real   , intent(in) :: eThermalDensity_C(nPoint)
    real   , intent(in) :: eThermalTemp_C(nPoint)

    !neutral atmosphere inputs
    integer, intent(in) :: nNeutral
    real   , intent(in) :: NeutralDens_IC(nNeutral,nIono)

    !incomming crossections
    real   , intent(in) :: SIGS(nNeutral,nEnergy)
    real   , intent(in) :: SIGI(nNeutral,nEnergy,nEnergy)
    real   , intent(in) :: SIGA(nNeutral,nEnergy,nEnergy)
    
    ! incomming photo electron production spectrum
    real   , intent(in) :: ePhotoProdSpec_IC(nEnergy,nIono)
  
    
    REAL h,p,Fcheck, sigmaO,muO,velt,eta, lbeta(nPoint),alpha(nAngle), &
         sigma(nAngle),Flastj,kk, del1,del2, newphi,beta(nPoint),sigO1, &
         Fsum, delE,Flasti, Flastt,Qstar, Theta,s1,s2
    INTEGER i,j,k,jj,Ist, &
         flag,count,warning,SPick,space,ii
    !inidicies that store full ionosphere grid(2*nIono points) and half 
    ! ionosphere grid (nIono points). These convert position along the field 
    ! line into these sub grids.
    integer :: iIono, iIonoHalf

    !set maximum iterations
    integer, parameter :: countmax=1200

    real :: cascade, lossum
    
    !Altitude range variables
    integer :: nAltMin, nAltMax
    
    ! array to hold temporarily the values being overwritten by BCs
    real :: fhold(0:nAngle)
    real,parameter :: cEVtoCMperS = 5.88e7 ! convert energy to velocity

    ! precipitation 
    real :: PrecipCoef,PolarRainCoef
    !---------------------------------------------------------------------------

    !initialize lbeta to 0
    lbeta(:)=0.0
    
    !      initialize warning to 0
    warning =0
    alpha(1)=1.
    sigma(1)=0.
    kk=1.
    !  Start the energy loop
    !      PRINT *, 'MainPlas, t=',t
    ENERGY: DO j=nEnergy,1,-1
       delE=.5*(DeltaE_I(j+1)+DeltaE_I(j))
       IF (j.EQ.nEnergy) delE=DeltaE_I(j)
       
       flag=1
       count=0
       !  Start the iteration loop
       ITERATION: DO WHILE ((flag.EQ.1).AND.(count.LT.countmax))
          flag=0
          count=count+1
          IF (count.GT.countmax-2) WRITE (10,*) 'Count: ',count

          ! Upward Region

          ! set BC when in iono2. Put BC in iphiup at top of iono1 and hold 
          ! the current value to put it back after the calculation
          if (.not.IsIono1)then
             do k=0,nThetaAlt_IIC(iLine,j,nIono)
                Fhold(k)=iphiup(iLine,k,nIono,j)
             end do

             ! fill BC from the the plasmasphere
             do k=0,nThetaAlt_IIC(iLine,j,nIono)
                iphiup(iLine,k,nIono,j)=phiup(iLine,k,nPlas,j)
             end do

             ! Add precip info here
          endif
             
          !set bounds for ionosphere loop depending if you are in iono 1 or 2
          if (IsIono1) then
             ! set upward bounds for iono1
             nAltMin = 1
             nAltMax = nIono
             h=FieldLineGrid_IC(iLine,nAltMin+1)-FieldLineGrid_IC(iLine,nAltMin)
          else
             ! set upward bounds for iono2
             nAltMin = nPoint-nIono+1
             nAltMax = nPoint
             h=FieldLineGrid_IC(iLine,nAltMin)-FieldLineGrid_IC(iLine,nAltMin-1)
          endif


          FIELDLINE_UPWARD: DO i=nAltMin,nAltMax
             !write(*,*) 'Start FieldLine up',j,i
             !velt is now altitude dependent
             velt=cEVtoCMperS*SQRT(KineticEnergy_IIC(iLine,j,i))*delt

             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_C(i),eThermalDensity_C(i),KineticEnergy_IIC(iLine,j,i))

             ! set iIono and iIonoHalf indices based on 1st or 2nd ionosphere
             ! remember, iIono is index spanning both ionospheres, 
             ! iIonoHalf is index spanning just one ionosphere
             if(IsIono1) then
                iIono=i
                iIonoHalf=i
             else
                iIono=i-nPlas
                iIonoHalf = nPoint-i+1
             end if

             !set the b-field ratio coeficient
             eta = (Bfield0_II(iLine,j)/Bfield_IC(iLine,i)) &
                  *(EnergyGrid_I(j)-DeltaPot_IC(iLine,i))   &
                  /(EnergyGrid_I(j)-DeltaPot0_II(iLine,j))
             
             ! get sigma0 and electron production (primary+secondary) 
             ! note: we currently assume no E|| in ionosphere so KE=totalE
             ! therefore this subroutine doesn't have to change. It would 
             ! if we start allowing E|| in ionosphere portion
             CALL get_sigma0_and_eprod(j,nNeutral,sigO1,&
                  ePhotoProdSpec_IC(j,iIonoHalf), &
                  NeutralDens_IC(:,iIonoHalf),SIGS,SIGI,SIGA,Qstar, &
                  specup(iLine,:,i),specdn(iLine,:,i),&
                  Qestar_IICI(iLine,:,iIono,j),Qpstar_ICI(iLine,iIono,j))
             Qstar_ICI(iLine,iIono,j)=Qstar
             
             do k=1,nThetaAlt_IIC(iLine,j,i)-1
                call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                     NeutralDens_IC(:,iIonoHalf),j,iIono,k,cascade,lossum,1)
                CALL ThetaVar1(nAngle,Theta,del1,del2,k,&
                     EqAngleGrid_IG(iLine,0:nAngle))
                IF (k.EQ.nThetaAlt_IIC(iLine,j,i)-1) del2=dThetaEnd_III(iLine,j,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=iphiup(iLine, k,iIono-1,j)
                Flastj=iphiup(iLine, k,iIono,j+1)
                Flastt=liphiup(iLine,k,iIono,j)
                
                !When we are in the ilocal region kill the transport terms in 
                ! coeficients by setting kk to 0
                kk=1
                if (iIonoHalf <= iLocal) kk=0
                !write(*,*) 'test3: i,j,k',i,j,k
                CALL NumCalcVar_pot(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     mu_IIIC(iLine,k,j,i),eta,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2,&
                     cascade,lossum,EnergyGrid_I(j),DeltaPot0_II(iLine,j), &
                     KineticEnergy_IIC(iLine,j,i),eThermalDensity_C(i))
                !            CALL CheckWarn(sigma(k+1),warning,1,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,1,*9999)
                
             end do

             ! when nThetaAlt_IIC fills the entire reference PA, as when 
             ! the main potential drop occurs just above the ionosphere 
             ! reflecting all low energy particles we need to treate the pi/2
             ! point like we do in the plasmasphere for the reference altitude
             if (nThetaAlt_IIC(iLine,j,i) == nAngle) then
                p=-2.*s1/dThetaEnd_III(iLine,j,i)
             else
                p=pend_pot(kk,s1,eta,EnergyGrid_I(j),DeltaPot0_II(iLine,j), &
                     KineticEnergy_IIC(iLine,j,i),eThermalDensity_C(i))
             endif

             Flastt=liphiup(iLine,nThetaAlt_IIC(iLine,j,i),iIono,j)
             Fsum=iphiup(iLine,nThetaAlt_IIC(iLine,j,i)-1,iIono,j)&
                  +iphidn(iLine,nThetaAlt_IIC(iLine,j,i)-1,iIono,j)
             Fcheck=iphiup(iLine,nThetaAlt_IIC(iLine,j,i),iIono,j)
             Flastj=iphiup(iLine,nThetaAlt_IIC(iLine,j,i),iIono,j+1)
             
             call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                  NeutralDens_IC(:,iIonoHalf),j,iIono,nThetaAlt_IIC(iLine,j,i),&
                  cascade,lossum,1)
             newphi=(Flastt/velt-p*Fsum/(2.*dThetaEnd_III(iLine,j,i)) &
                  +s2*(lbeta(i)* Flastj/delE+Qstar+cascade))/ &
                  (1/velt-p/dThetaEnd_III(iLine,j,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,2,*9999)
             
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             iphiup(iLine,nThetaAlt_IIC(iLine,j,i),iIono,j)=newphi

             do k=nThetaAlt_IIC(iLine,j,i),1,-1
                newphi=alpha(k)*iphiup(iLine,k,iIono,j)+sigma(k)
                Fcheck=iphiup(iLine,k-1,iIono,j)
                !            CALL CheckWarn(newphi,warning,3,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                
                
                iphiup(iLine,k-1,iIono,j)=newphi
             end do
             CALL midpnt_int(specup(iLine,j,i),iphiup(iLine,:,iIono,j),&
                  mu_IIIC(iLine,:,j,i),1, nThetaAlt_IIC(iLine,j,i)+1,nAngle+1,1)
             specup(iLine,j,i)=-.5*specup(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) 'Upward region'
                WRITE (10,*) time,j,i,iIono,nThetaAlt_IIC(iLine,j,i),&
                     flag,(iphiup(iLine,k,iIono,j),k=0,nThetaAlt_IIC(iLine,j,i))
             END IF
             
             if (i < nPoint) then
                h=FieldLineGrid_IC(iLine,i+1)-FieldLineGrid_IC(iLine,i)
             else
                h=FieldLineGrid_IC(iLine,i)-FieldLineGrid_IC(iLine,i-1)
             endif
             ! End Upward Region of ionosphere
          end DO FIELDLINE_UPWARD
          
          !put back the iphiup at top of iono1 (which was storing plasmasph bc)
          if(.not.IsIono1)then
             do k=0,nThetaAlt_IIC(iLine,j,nIono)
                iphiup(iLine,k,nIono,j)=Fhold(k)
             end do
          endif
          
          !*  Downward Region of the ionosphere

          ! set BC when in iono1. Put BC in iphidn at top of iono2 and hold 
          ! the current value to put it back after the calculation
          if (IsIono1)then
             do k=0,nThetaAlt_IIC(iLine,j,nIono)
                Fhold(k)=iphidn(iLine,k,nIono+1,j)
             end do

             ! fill BC from the the plasmasphere
             do k=0,nThetaAlt_IIC(iLine,j,nIono)
                iphidn(iLine,k,nIono+1,j)=phidn(iLine,k,1,j)
             end do

             ! Add precip info here
             if(UsePrecipitation .and. EnergyGrid_I(j) > PrecipEmin &
                  .and. .not.UseOvation) then 
                do k=0,nThetaAlt_IIC(iLine,j,nIono)
                   !These current values are for soft electron precipitation 
                   ! from strangeway et al 2005. In future should come from 
                   ! PWOM
                   !PrecipCoef=get_precip_norm(400.0,100.0,1000.0,2.0)
                   PrecipCoef=get_precip_norm(PrecipEmean,PrecipEmin,&
                        PrecipEmax,PrecipEflux)
                   iphidn(iLine,k,nIono+1,j)= &
                        PrecipCoef*EnergyGrid_I(j)*exp(-EnergyGrid_I(j)&
                        /PrecipEmean)
                end do
             elseif(UseOvation .and. EnergyGrid_I(j) > OvationEmin &
                  .and. .not.UsePrecipitation) then 
                !Start with diffuse aurora
                PrecipCoef=get_precip_norm(EmeanDiff,OvationEmin,&
                     OvationEmax,EfluxDiff)
                do k=0,nThetaAlt_IIC(iLine,j,nIono)
                   iphidn(iLine,k,nIono+1,j)= &
                        PrecipCoef*EnergyGrid_I(j)*exp(-EnergyGrid_I(j)&
                        /EmeanDiff)
                enddo
                !add in wave aurora
                PrecipCoef=get_precip_norm(EmeanWave,OvationEmin,&
                     OvationEmax,EfluxWave)
                do k=0,nThetaAlt_IIC(iLine,j,nIono)
                   iphidn(iLine,k,nIono+1,j)= iphidn(iLine,k,nIono+1,j) + &
                        PrecipCoef*EnergyGrid_I(j)*exp(-EnergyGrid_I(j)&
                        /EmeanWave)
                enddo
                !add in Monoenergetic aurora
                PrecipCoef=get_precip_norm(EmeanMono,OvationEmin,&
                     OvationEmax,EfluxMono)
                do k=0,nThetaAlt_IIC(iLine,j,nIono)
                   iphidn(iLine,k,nIono+1,j)= iphidn(iLine,k,nIono+1,j) + &
                        PrecipCoef*EnergyGrid_I(j)*exp(-EnergyGrid_I(j)&
                        /EmeanMono)
                enddo
             end if
             
             ! add in polar rain
             if(UsePolarRain .and. EnergyGrid_I(j) > PolarRainEmin ) then 
                do k=0,nThetaAlt_IIC(iLine,j,nIono)
                   PolarRainCoef=get_precip_norm(PolarRainEmean,PolarRainEmin,&
                        PolarRainEmax,PolarRainEflux)
                   if (UsePrecipitation.and. EnergyGrid_I(j) > PrecipEmin) then
                      !add polar rain on top of the precipitation set
                      iphidn(iLine,k,nIono+1,j)= iphidn(iLine,k,nIono+1,j) + &
                           PolarRainCoef*EnergyGrid_I(j)*exp(-EnergyGrid_I(j)&
                           /PolarRainEmean)
                   else
                      !set polar rain values instead of adding
                      iphidn(iLine,k,nIono+1,j)= &
                           PolarRainCoef*EnergyGrid_I(j)*exp(-EnergyGrid_I(j)&
                           /PolarRainEmean)
                   endif
                end do
             end if
          endif

          !set bounds for ionosphere loop depending if you are in iono 1 or 2
          if (IsIono1) then
             ! set upward bounds for iono1
             nAltMin = 1
             nAltMax = nIono
             h=FieldLineGrid_IC(iLine,nAltMax)-FieldLineGrid_IC(iLine,nAltMax+1)
          else
             ! set upward bounds for iono2
             nAltMin = nPoint-nIono+1
             nAltMax = nPoint
             h=FieldLineGrid_IC(iLine,nAltMax-1)-FieldLineGrid_IC(iLine,nAltMax)
          endif
          
          FIELDLINE_DOWN: DO i=nAltMax,nAltMin,-1
             !             write(*,*) 'Start FieldLine Down'
                          !velt is now altitude dependent
             velt=cEVtoCMperS*SQRT(KineticEnergy_IIC(iLine,j,i))*delt
             CALL CoulVar(i,beta(i),sigO1,eThermalTemp_C(i),eThermalDensity_C(i),KineticEnergy_IIC(iLine,j,i))

             ! set iIono and iIonoHalf indices based on 1st or 2nd ionosphere
             ! remember, iIono is index spanning both ionospheres, 
             ! iIonoHalf is index spanning just one ionosphere
             if(IsIono1) then
                iIono=i
                iIonoHalf=i
             else
                iIono=i-nPlas
                iIonoHalf = nPoint-i+1
             end if

             !set the b-field ratio coeficient
             eta = (Bfield0_II(iLine,j)/Bfield_IC(iLine,i)) &
                  *(EnergyGrid_I(j)-DeltaPot_IC(iLine,i))   &
                  /(EnergyGrid_I(j)-DeltaPot0_II(iLine,j))
             ! get sigma0 and electron production (primary+secondary)
             CALL get_sigma0_and_eprod(j,nNeutral,sigO1,&
                  ePhotoProdSpec_IC(j,iIonoHalf), &
                  NeutralDens_IC(:,iIonoHalf),SIGS,SIGI,SIGA,Qstar, &
                  specup(iLine,:,i),specdn(iLine,:,i),&
                  Qestar_IICI(iLine,:,iIono,j),Qpstar_ICI(iLine,iIono,j))
             Qstar_ICI(iLine,iIono,j)=Qstar

             DO k=1,nThetaAlt_IIC(iLine,j,i)-1
                call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                     NeutralDens_IC(:,iIonoHalf),j,iIono,k,cascade,lossum,2)
                CALL ThetaVar1(nAngle,Theta,del1,del2,k, &
                     EqAngleGrid_IG(iLine,0:nAngle))
                Theta=cPi-Theta
                del1=-del1
                del2=-del2
                IF (k.EQ.nThetaAlt_IIC(iLine,j,i)-1) del2=-dThetaEnd_III(iLine,j,i)
                muO=cos(Theta)
                sigmaO=sigO1
                s1=sigmaO
                s2=1.
                Flasti=iphidn(iLine,k,iIono+1,j)
                Flastj=iphidn(iLine,k,iIono,j+1)
                Flastt=liphidn(iLine,k,iIono,j)
                
                !When we are in the ilocal region kill the transport terms in 
                ! coeficients by setting kk to 0
                kk=1
                if (iIonoHalf <= iLocal) kk=0
                !write(*,*) 'test4: i,j,k',i,j,k
                CALL NumCalcVar_pot(alpha(k+1),alpha(k),sigma(k+1),sigma(k), &
                     -mu_IIIC(iLine,k,j,i),eta,muO,del1,del2,beta(i),lbeta(i), &
                     Flasti,Flastj,Flastt,velt,h,delE,Qstar,kk,s1,s2, &
                     cascade,lossum,EnergyGrid_I(j),DeltaPot0_II(iLine,j), &
                     KineticEnergy_IIC(iLine,j,i),eThermalDensity_C(i))
                !            CALL CheckWarn(sigma(k+1),warning,4,*9999)
                !            CALL CheckWarn(alpha(k+1),warning,4,*9999)
             enddo
             ! when nThetaAlt_IIC fills the entire reference PA, as when 
             ! the main potential drop occurs just above the ionosphere 
             ! reflecting all low energy particles we need to treate the pi/2
             ! point like we do in the plasmasphere for the reference altitude
             if (nThetaAlt_IIC(iLine,j,i) == nAngle) then
                p=2.*s1/dThetaEnd_III(iLine,j,i)
             else
                p=-pend_pot(kk,s1,eta,EnergyGrid_I(j),DeltaPot0_II(iLine,j), &
                     KineticEnergy_IIC(iLine,j,i),eThermalDensity_C(i))
             endif
             Flastt=liphidn(iLine,nThetaAlt_IIC(iLine,j,i),iIono,j)
             Fsum=iphidn(iLine,nThetaAlt_IIC(iLine,j,i)-1,iIono,j)&
                  +iphiup(iLine,nThetaAlt_IIC(iLine,j,i)-1,iIono,j)
             Fcheck=iphidn(iLine,nThetaAlt_IIC(iLine,j,i),iIono,j)
             Flastj=iphidn(iLine,nThetaAlt_IIC(iLine,j,i),iIono,j+1)
             call get_cascade_and_lossum(iLine,nNeutral,SIGA,&
                  NeutralDens_IC(:,iIonoHalf),j,iIono,nThetaAlt_IIC(iLine,j,i), &
                  cascade,lossum,2)
             newphi=(Flastt/velt+p*Fsum/(2.*dThetaEnd_III(iLine,j,i))+s2*(lbeta(i)* &
                  Flastj/delE+Qstar+cascade))&
                  /(1/velt+p/dThetaEnd_III(iLine,j,i)+s2*(beta(i)/delE+lossum))
             !            CALL CheckWarn(newphi,warning,5,*9999)
             CALL CheckConv(newphi,Fcheck,epsilon,flag)
             iphidn(iLine,nThetaAlt_IIC(iLine,j,i),iIono,j)=newphi
             DO k=nThetaAlt_IIC(iLine,j,i),1,-1
                newphi=alpha(k)*iphidn(iLine,k,iIono,j)+sigma(k)
                Fcheck=iphidn(iLine,k-1,iIono,j)
                !            CALL CheckWarn(newphi,warning,6,*9999)
                CALL CheckConv(newphi,Fcheck,epsilon,flag)
                iphidn(iLine,k-1,iIono,j)=newphi
             end DO
             CALL midpnt_int(specdn(iLine,j,i),iphidn(iLine,:,iIono,j),&
                  mu_IIIC(iLine,:,j,i),1, nThetaAlt_IIC(iLine,j,i)+1,nAngle+1,1)
             specdn(iLine,j,i)=.5*specdn(iLine,j,i)
             !  Write to a file (if the solution isn't converging)
             IF (count.GT.countmax-2) THEN
                WRITE (10,*) 'Downward region'
                WRITE (10,*) time,j,i,iIono,flag,&
                     (iphidn(iLine,k,iIono,j),k=0,nThetaAlt_IIC(iLine,j,i))
             END IF
             if (i==1) then
                h=FieldLineGrid_IC(iLine,i)-FieldLineGrid_IC(iLine,i+1)        
             else
                h=FieldLineGrid_IC(iLine,i-1)-FieldLineGrid_IC(iLine,i)         
             endif
             ! End Downward Region of ionosphere
          end DO FIELDLINE_DOWN
          !put back the iphidn at top of iono2 (which was storing plasmasph bc)
          if (IsIono1)then
             do k=0,nThetaAlt_IIC(iLine,j,nIono)
                iphidn(iLine,k,nIono+1,j)=Fhold(k)
             end do
          endif
       END DO ITERATION                  ! End of iteration loop
       CALL CheckFlag(flag,1,warning,-1,*9999)
       !  Save the calculated values for the next time and energy steps
       if (IsIono1) then
          lbeta(1:nIono)=beta(1:nIono)
       else
          lbeta(nPoint-nIono+1:nPoint)=beta(nPoint-nIono+1:nPoint)
       endif

    end DO ENERGY
    
9999 IF (warning.GT.0) THEN
       PRINT *, 'Negative Densities Occurred:',newphi,warning
       PRINT *, time,j,i,k,Bfield_IC(iLine,i),&
            BField0_II(iLine,j),FieldLineGrid_IC(iLine,i),h,Theta,muO,&
            mu_IIIC(iLine,k,j,i),del1, dThetaEnd_III(iLine,j,i),beta(i),sigmaO,p,count,&
            Qstar
       PRINT *, iphiup(iLine,k,iIono,j),iphidn(iLine,k,iIono,j)
       STOP
    ELSE IF (warning.LT.0) THEN
       PRINT *, 'No convergence: ',count,Time,warning,j
       STOP
    END IF
    !  Format for the nonconvergent output
!!!27  FORMAT (100E10.3,4I4,1P,100E10.3)
    
    RETURN
  END SUBROUTINE update_se_state_iono_pot


  !=============================================================================
  !* ------------------------------------------------------------------ **
!  Subroutine NumCalcVar finds alpha(k) and sigma(k).
!*  VARIABLE DESCRIPTIONS
!*      alpha   Flux attenuator for this pitch angle
!*      alp      Flux attenuator for the previous pitch angle
!*      sigma   Flux offset for this angle: [flux]
!*      sig      Flux offset for the previous angle; [flux]
!*      del1    Pitch angle step between this and the previous step; rad
!*      del2    Pitch angle step between this and the next step; rad
!*      muO     Cosine of the equatorial pitch angle
!*      mu      Cosine of the actual pitch angle
!*      BB      Magnetic field strength for this step; G
!*      BFieldEq_I(iLine)      Magnetic field strength at the equator; G
!*      beta    Energy loss parameters; ????
!*      lbeta   Previous energy step's beta; ????
!*      Fj      Flux for previous E for same space, time, and angle step
!*      Fi      Flux for previous space for same angle, time, and E step
!*      Ft      Flux for previous time for same E, space, and angle step
!*      vt      Velocity of an electron at this energy * delt; cm
!*      h       Size of this step; cm
!*      dE      Energy step size; eV
!*      Qstar   Electron production rate(primary plus secondary;
!*                cm-3 s-1 rad-1 ????
!*      p
!*      q
!*      a,b,c,d
!*
  SUBROUTINE NumCalcVar(alpha,alp,sigma,sig,mu,coef,muO,del1, &
       del2,beta,lbeta,Fi,Fj,Ft,vt,h,dE,Qstar,kk,s1,s2,cascade,lossum)
    REAL alpha,alp,sigma,sig,mu,coef,muO,del1,del2,beta,lbeta,s1,s2, &
         Fi,Fj,Ft,vt,h,dE,a,b,c,d,p,q,Qstar,lossum,cascade,kk
    
    !write(*,*) 'muO, BB, BFieldEq_I(iLine)',muO, BB, BFieldEq_I(iLine)
    p=kk*s1*(muO**4+coef-1)/((muO**3)*(1-muO**2)**.5)
    q=kk*s1*ABS((coef-1+muO**2)/muO**2)
    a=(2*q-p*del2)/(del1*(del1+del2))
    b=(2*q+p*del1)/(del2*(del1+del2))
    c=(2*q+(del1-del2)*p)/(del1*del2)+s2*(kk*mu/h+beta/dE+lossum)+1/vt
    d=s2*(lbeta*Fj/dE+kk*mu*Fi/h+Qstar+cascade)+Ft/vt

    IF (b.LT.0.) b=0.
    alpha=b/(c-a*alp)
    sigma=(a*sig+d)/(c-a*alp)
    IF ((d.LT.0).OR.(d-d.NE.0)) THEN
       PRINT *,'NEG D ',p,q,a,b,c,d,sigma,sig,alpha,alp
       PRINT *,mu,h,dE,vt,Fi,Fj,Ft,Qstar,cascade
       call con_stop('')
    ELSE IF ((sigma.LT.0).OR.(sigma-sigma.NE.0)) THEN
       PRINT *,'NEG SIGMA ',p,q,a,b,c,d,sigma,sig,alpha,alp
       PRINT *,mu,h,dE,vt,Fi,Fj,Ft,Qstar
    ELSE IF ((alpha.LT.0).OR.(alpha-alpha.NE.0)) THEN
       PRINT *,'NEG ALPHA ',p,q,a,b,c,d,sigma,sig,alpha,alp
       PRINT *,mu,h,dE,vt,Fi,Fj,Ft,Qstar
    END IF
    !      PRINT 10, BB,BFieldEq_I(iLine)/BB,muO,del1,del2,p,q,a,b,alpha
10  FORMAT (10(1PG11.4,1X))
    RETURN
  END SUBROUTINE NumCalcVar

  !=============================================================================
  !* ------------------------------------------------------------------ **
!  Subroutine NumCalcVar_pot finds alpha(k) and sigma(k) in the case of 
!   field aligned potential.
!*  VARIABLE DESCRIPTIONS
!*      alpha   Flux attenuator for this pitch angle
!*      alp      Flux attenuator for the previous pitch angle
!*      sigma   Flux offset for this angle: [flux]
!*      sig      Flux offset for the previous angle; [flux]
!*      del1    Pitch angle step between this and the previous step; rad
!*      del2    Pitch angle step between this and the next step; rad
!*      muO     Cosine of the equatorial pitch angle
!*      mu      Cosine of the actual pitch angle
!*      BB      Magnetic field strength for this step; G
!*      BFieldEq_I(iLine)      Magnetic field strength at the equator; G
!*      beta    Energy loss parameters; ????
!*      lbeta   Previous energy step's beta; ????
!*      Fj      Flux for previous E for same space, time, and angle step
!*      Fi      Flux for previous space for same angle, time, and E step
!*      Ft      Flux for previous time for same E, space, and angle step
!*      vt      Velocity of an electron at this energy * delt; cm
!*      h       Size of this step; cm
!*      dE      Energy step size; eV
!*      Qstar   Electron production rate(primary plus secondary;
!*                cm-3 s-1 rad-1 ????
!*      p
!*      q
!*      a,b,c,d
!*
  SUBROUTINE NumCalcVar_pot(alpha,alp,sigma,sig,mu,eta,muO,del1, &
       del2,beta,lbeta,Fi,Fj,Ft,vt,h,dE,Qstar,kk,s1,s2,cascade,lossum,&
       TotalE,DeltaPot0, KE, ne)
    REAL alpha,alp,sigma,sig,mu,eta,muO,del1,del2,beta,lbeta,s1,s2, &
         Fi,Fj,Ft,vt,h,dE,a,b,c,d,p,q,Qstar,lossum,cascade,kk, pstar
    real, intent(in) ::TotalE,DeltaPot0, KE, ne
    real, parameter  :: Acoef = 2.6e-12 !ev^2cm^2
    !--------------------------------------------------------------------------
    !write(*,*) 'muO, BB, BFieldEq_I(iLine)',muO, BB, BFieldEq_I(iLine)
    p=kk*s1*(muO**4+eta-1.0)/((muO**3)*(1-muO**2)**.5)
    ! pstar was zero when there is no potential
    pstar = kk*Acoef*ne/2.0/(KE**2)*(sqrt(1.0-muO**2)/muO)&
         *((KE/(TotalE-DeltaPot0))-1.0)
    q=kk*s1*ABS((eta-1+muO**2)/muO**2)

    a=(2*q+(pstar-p)*del2)/(del1*(del1+del2))
    b=(2*q-(pstar-p)*del1)/(del2*(del1+del2))
    c=(2*q+(del1-del2)*(p-pstar))/(del1*del2)+s2*(kk*mu/h+beta/dE+lossum)+1/vt
    d=s2*(lbeta*Fj/dE+kk*mu*Fi/h+Qstar+cascade)+Ft/vt

    IF (b.LT.0.) b=0.
    alpha=b/(c-a*alp)
    sigma=(a*sig+d)/(c-a*alp)
    IF ((d.LT.0).OR.(d-d.NE.0)) THEN
       PRINT *,'NEG D: '
       write(*,*) 'd=',d
       write(*,*) 'd is made of:s2,lbeta,Fj,DE,kk,mu,Fi,h,Qstar,cascade,Ft',&
            s2,lbeta,Fj,DE,kk,mu,Fi,h,Qstar,cascade,Ft
       call con_stop('')
    ELSE IF ((sigma.LT.0).OR.(sigma-sigma.NE.0)) THEN
       PRINT *,'WARNING NEG SIGMA:'
       write(*,*) 'sigma=',sigma
       write(*,*) 'sigma is made of:a,sig,d,c,a,alp',&
            a,sig,d,c,a,alp
       write(*,*) 'a is made of:q,pstar,p,del1,del2',&
            q,pstar,p,del1,del2
       write(*,*) 'c is made of:q,pstar,p,del1,del2,kk,mu,h,beta,dE,lossum,vt',&
            q,pstar,p,del1,del2,kk,mu,h,beta,dE,lossum,vt
       write(*,*) 'd is made of:s2,lbeta,Fj,dE,kk,mu,Fi,h,Qstar,cascade,Ft,vt',&
            s2,lbeta,Fj,dE,kk,mu,Fi,h,Qstar,cascade,Ft,vt
       write(*,*) 'pstar is made of:Acoef,ne,KE,muO',&
            Acoef,ne,KE,muO
       !call con_stop('')
       write(*,*) 'setting sigma to 0 to proceed'
       sigma=0.0
    ELSE IF ((alpha.LT.0).OR.(alpha-alpha.NE.0)) THEN
       PRINT *,'NEG ALPHA ',p,q,a,b,c,d,sigma,sig,alpha,alp
       PRINT *,mu,h,dE,vt,Fi,Fj,Ft,Qstar
    END IF
    !      PRINT 10, BB,BFieldEq_I(iLine)/BB,muO,del1,del2,p,q,a,b,alpha
10  FORMAT (10(1PG11.4,1X))
    RETURN
  END SUBROUTINE NumCalcVar_pot
  !=============================================================================
  ! find p at local pi/2 when not at reference point (case for potential)
  real function pend_pot(kk,s1,eta,TotalE,DeltaPot0, KE, ne)
    real, intent(in) ::kk,s1,eta,TotalE,DeltaPot0, KE, ne
    real :: muO, p, pstar
    real, parameter  :: Acoef = 2.6e-12 !ev^2cm^2
    !--------------------------------------------------------------------------
    ! set mu0 = sqrt(1-eta) which is the case for local pi/2
    muO = sqrt(max(1.0-eta,0.0))
    
    !get p and pstar at this location
    p=kk*s1*(muO**4+eta-1.0)/((muO**3)*(1-muO**2)**.5)
    
    pstar = Acoef*ne/2.0/(KE**2)*(sqrt(1.0-muO**2)/muO)&
         *((KE/(TotalE-DeltaPot0))-1.0)

    ! now set pend as defined on page 150 of Mike Leimohn's thesis
    pend_pot = p-pstar
  end function pend_pot
  !=============================================================================
  !* ------------------------------------------------------------------ **
  !  Subroutine get_sigmaO_and_eprod (was IonoVar1) 
  !*  calculates beta, sigmaO, and Qstar for the ionosphere.
  !*  VARIABLE DESCRIPTIONS
  !*      EnIon   Energies for ionization states for each species; eV
  !*      SIGS    Elastic collision cross sections, species, energy; cm2
  !*      SIGI    Differential cross sections, state,species,sec,pri; cm2
  !*      ionsum  Sum of energy losses due to ionization; ???? [beta]
  !*      elassum Sum of elastic scattering cross sections; ???? [beta]
  !*      exsum   Sum of energy losses due to excitation; ???? [beta]
  !*      prodsum Sum of secondary electron production; cm-2 s-1 eV-1 sr-1
  !*      BINNUM  Function that finds the array position for a given E
  !*      jj,m,n  Energy, state, and species loop counters
  !*      PE      Photoelectron production rate in the ionos.; cm-3 s-1
  !*      pi      3.1415926536
  !*      j      Current energy increment
  !*      dens    Thermal electron density for this step; cm-3
  !*      Qstar   Electron production rate(primary plus secondary;
  !*                cm-3 s-1 rad-1 ????
  !*  FLUX VARIABLES; cm-2 s-1 eV-1 sr-1
  !*      phiup  Flux in the ionospheres (moving from first to second)
  !*      phidn  Flux in the ionospheres (moving from second to first)
  !*
  SUBROUTINE get_sigma0_and_eprod(iEnergyIn,nNeutral,&
       sigO1,PE,NeutralDens_I,SIGS,SIGI,SIGA,Qstar,&
       OmniDirFluxUp_I,OmniDirFluxDn_I,Qe,Qp)
    use ModSeGrid,      ONLY: FieldLineGrid_IC,nIono,nEnergy, nPoint, &
         DeltaE_I,EnergyGrid_I,BINNUM
    use ModMath,        ONLY: midpnt_int
    use ModNumConst,    ONLY: cPi

    IMPLICIT NONE
    integer, intent(in) :: iEnergyIn,nNeutral
    real   , intent(inout) :: sigO1
    real   , intent(in) :: PE,NeutralDens_I(nNeutral)
    !incomming crossections
    real   , intent(in) :: SIGS(nNeutral,nEnergy)
    real   , intent(in) :: SIGI(nNeutral,nEnergy,nEnergy)
    real   , intent(in) :: SIGA(nNeutral,nEnergy,nEnergy)
    !outgoing electron production rate
    real   , intent(out):: Qstar
    
    real   , intent(in) :: OmniDirFluxUp_I(nEnergy),OmniDirFluxDn_I(nEnergy)
    !outgoing electron production rate from secondary production
    real   , intent(out):: Qe(nNeutral)
    !incomming electron production rate from precip
    real   , intent(in) :: Qp

    real   :: NetFlux, SecProd(nNeutral,nEnergy)
    REAL   :: elassum
    
    ! Minimum ionization threshold
    real,parameter :: Eplus = 12.0 
    INTEGER m,n,LL,jj
    !---------------------------------------------------------------------------
        
    ! this part gives the sum elastic scattering cross sections
    elassum=0.
    Qe(:)=0.
    SecProd(:,:)=0.

    ! start with photoelectron and incoming e- production rate
    Qstar=PE/(4*cPi*DeltaE_I(iEnergyIn))+Qp

    do n=1,nNeutral
       elassum=elassum+NeutralDens_I(n)*SIGS(n,iEnergyIn)
    end do
    sigO1=sigO1+.5*elassum
    
    ! this calculates the electron total production
    DO  n=1,nNeutral
       IF (2*EnergyGrid_I(iEnergyIn)+Eplus &
            < EnergyGrid_I(nEnergy)+.5*DeltaE_I(nEnergy)) THEN
          LL=BINNUM(2*EnergyGrid_I(iEnergyIn)+Eplus)

          DO jj=LL,nEnergy
             NetFlux=OmniDirFluxUp_I(jj)-OmniDirFluxDn_I(jj)
             SecProd(n,jj)=SecProd(n,jj)+ &
                  NeutralDens_I(n)*SIGI(n,iEnergyIn,jj)*NetFlux
          enddo
       END IF
       ! energy-integrated secondary production per species into Qe
       CALL midpnt_int(Qe(n),SecProd(n,:),DeltaE_I,1,nEnergy,nEnergy,2)
       ! add in species secondary production to total e- production
       Qstar = Qstar+Qe(n)
    enddo
    
    RETURN
  end SUBROUTINE get_sigma0_and_eprod


  !=============================================================================
  !  Subroutine get_cascade_and_lossum (was IonoVar2) finds cascade and lossum.
  !     m       Specifies the use of iphiup or iphidn
  SUBROUTINE get_cascade_and_lossum(iLine,nNeutral,SIGA,NeutralDens_I, &
       iEnergyIn,iIono,iAngle,cascade,lossum,m)

    use ModSeGrid,      ONLY: FieldLineGrid_IC,nIono,nEnergy, nPoint, &
         DeltaE_I,EnergyGrid_I, EnergyMin

    IMPLICIT NONE    
    integer, intent(in) :: iLine, nNeutral
    real   , intent(in) :: SIGA(nNeutral,nEnergy,nEnergy)
    real   , intent(in) :: NeutralDens_I(nNeutral)
    integer, intent(in) :: iEnergyIn,iIono,iAngle
    real   , intent(out):: cascade,lossum
    integer, intent(in) :: m
    
    real :: flux
    integer :: jj,n,LL
    !---------------------------------------------------------------------------
    lossum=0.
    cascade=0.
    
    ! New method (with Swartz fix and cross section splitting)           
    do jj=1,iEnergyIn
       do n=1,nNeutral
          lossum=lossum+NeutralDens_I(n)*SIGA(n,jj,iEnergyIn)
       end do
    end do
    IF (iEnergyIn.LT.nEnergy) THEN
       ! what this is doing: For each energy step jj above the current energy 
       ! iEnergyIn, calculate and sum up the flux that has cascaded down into
       ! this energy bin
       do jj=iEnergyIn+1,nEnergy
          LL=jj-iEnergyIn
          ! IF (m.EQ.1) flux=AngIonJ(i,j,k,i,jj,iphiup(0,i1,jj))
          ! IF (m.EQ.2) flux=AngIonJ(i,j,k,i,jj,iphidn(0,i1,jj))
          IF (m.EQ.1) flux=iphiup(iLine,iAngle,iIono,jj)
          IF (m.EQ.2) flux=iphidn(iLine,iAngle,iIono,jj)
          do n=1,nNeutral
             cascade=cascade+NeutralDens_I(n)*SIGA(n,LL,jj)*flux
          end do
       end do
    END IF
    
    RETURN
  end SUBROUTINE get_cascade_and_lossum


  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Subroutine IonoVar3 sets collisional energy loss terms to zero,
  !*  the plasmasperic case.
  SUBROUTINE IonoVar3(Qstar,lossum,cascade)
    real, intent(out) :: Qstar,lossum,cascade
    Qstar=0.
    lossum=0.
    cascade=0.
  END SUBROUTINE IonoVar3

  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Finds beta and sigO, the coulomb collision coefficients.
  !NOTE: CHECK HOW SRC is defined
  SUBROUTINE CoulVar(i,beta,sigO1,TE,ZE,KE)
    use ModMath, only: G,erf
    REAL beta,sigO1,TE,ZE,KE,A,MC,X
    INTEGER n,i                              !,Iterm1,Iterm2
    
    DATA A/2.6E-12/, MC/1837./
    !write(*,*) 'KE, TE',KE, TE
    
    !kludge set the sigO1 and beta as in Mike's thesis
    !beta  = A*ZE/KE
    !sigO1 = 0.5*A*sigO1/(KE*KE)
    !return
    !endkludge
    
    X=SQRT(KE/TE)
    beta=ZE/TE*G(X)                  ! ion infl. --> zero
    sigO1=ZE*(1.+erf(X)-G(X))            ! ion infl. = ZE (thus the 1.
    do n=1,IRC                  ! Ring current influence
       X=SQRT(KE*MC*MRC(n)/TRC(n))
       beta=beta+SRC(i)*NRC(n)/TRC(n)*G(X)
       sigO1=sigO1+SRC(i)*NRC(n)*(erf(X)-G(X))
    end do
    beta=2.*A*beta
    sigO1=.25*A*sigO1/(KE*KE)
    RETURN
  END SUBROUTINE CoulVar

  !=============================================================================
  !* ------------------------------------------------------------------ **
  !  Subroutine ThetaVar1 finds Theta, del1, and del2.
  SUBROUTINE ThetaVar1(nAngle,Theta,del1,del2,k,grid)
    INTEGER, intent(in) :: nAngle
    REAL Theta,del1,del2,grid(0:nAngle)
    INTEGER k,kk
    !*
    kk=k
    Theta=grid(kk)
    del1=grid(kk)-grid(kk-1)
    del2=grid(kk+1)-grid(kk)
    RETURN
  END SUBROUTINE ThetaVar1
  !=============================================================================
  ! function to return normalization value for Maxwellian precipitation 
  real function  get_precip_norm(E0,E1,E2,eFlux)
    real, intent(in) :: E0, E1, E2, eFlux ! average, min, and max energy of precip
                                    ! eflux is integrated energy flux in 
                                    ! ergs/cm^2
    get_precip_norm = 1.98774e11*eFlux/&
         (exp(-E1/E0)*E0*(E1**2.0+2.0*E1*E0+2.0*E0**2.0)&
         -exp(-E2/E0)*E0*(E2**2.0+2.0*E2*E0+2.0*E0**2.0))
    return
  end function get_precip_norm
  
  !============================================================================
  ! function to return normalization value for Gaussian precipitation  
  real function  get_precip_norm_gaussian(E0,E1,E2,eFlux)
    real, intent(in) :: E0, E1, E2, eFlux ! average, min, and max energy of precip
    real sigmaG, x1, x2 ! width of the Gaussian, substitutions from integration
    ! eflux is integrated energy flux in                               
    ! ergs/cm^2                                                        
    sigmaG = E0/10.0
    x1 = (E1 - E0)/(sqrt(2.0)*sigmaG)
    x2 = (E2 - E0)/(sqrt(2.0)*sigmaG)
    
    
    get_precip_norm_gaussian = 1.98774e11*eFlux/&
         ((sqrt(2.0)*sigmaG/2.0)*(sqrt(2.0)*sigmaG*(exp(-x1**2.0) &
         - exp(-x2**2.0)) + E0*1.77245*(erf(x2) - erf(x1))))
    ! 1.77245 is sqrt(pi)                                                
    return
  end function get_precip_norm_gaussian


  !=============================================================================


  
  !* ------------------------------------------------------------------ **
  !  Subroutine CheckConv sees if a flux has converged or not.
  !*  VARIABLE DESCRIPTIONS
  !*      epsil   Convergence parameter in the iteration loop
  !*      flag    Indicates whether the fluxes have converged
  !*  FLUX VARIABLES; cm-2 s-1 eV-1 sr-1
  !*        flux      Current flux value to be checked
  !*      oldflux      Previous value to check against
  !*
  SUBROUTINE CheckConv(flux,oldflux,epsil,flag, DoReportError, ReturnError)
    REAL   ,  intent(in) :: flux,oldflux,epsil
    INTEGER, intent(out) :: flag
    LOGICAL,optional,intent(in) :: DoReportError
    real,optional,intent(out) :: ReturnError

    real :: error
    !--------------------------------------------------------------------------
    
!    IF (flux.GT.1E-20) THEN
    IF (flux.GT.1E-19) THEN
       error=ABS(flux-oldflux)/flux
       IF (error.GT.epsil) flag=1
    END IF
    if (present(DoReportError))then
       if(DoReportError) write(*,*) 'Error is:',error,'Epsilon is:',epsil
    endif
    
    if (present(ReturnError)) ReturnError=error
    RETURN
  END SUBROUTINE CheckConv

  !=============================================================================


  
  !* ------------------------------------------------------------------ **
  !  Subroutine CheckConv_full sees if a flux has converged or not.
  !  uses twice the convergence criterial for full flux as opposed to 
  !  criteria for pi/2
  !*  VARIABLE DESCRIPTIONS
  !*      epsil   Convergence parameter in the iteration loop
  !*      flag    Indicates whether the fluxes have converged
  !*  FLUX VARIABLES; cm-2 s-1 eV-1 sr-1
  !*        flux      Current flux value to be checked
  !*      oldflux      Previous value to check against
  !*
  SUBROUTINE CheckConv_full(flux,oldflux,epsil,flag, DoReportError, ReturnError)
    REAL   ,  intent(in) :: flux,oldflux,epsil
    INTEGER, intent(out) :: flag
    LOGICAL,optional,intent(in) :: DoReportError
    real,optional,intent(out) :: ReturnError

    real :: error
    !--------------------------------------------------------------------------
    
    error = 0.0
    IF (flux.GT.1E-19) THEN
       error = ABS(flux-oldflux)/flux
       IF (error > 2*epsil) flag=1
    END IF
    if (present(DoReportError) .and. DoReportError) then
       write(*,*) 'Error is:',error,'Epsilon is:',epsil
    endif
    
    if (present(ReturnError)) ReturnError=error

  END SUBROUTINE CheckConv_full
  !=============================================================================
  
  !* ------------------------------------------------------------------ **
  !  Subroutine CheckFlag sees if the loop  ended due to convergence or
  !  not.
  !*  VARIABLE DESCRIPTIONS
  !*      flag    Indicates whether the fluxes have converged
  !*      check      The number to check flag against
  !*      warn    Indicates (nonzero) negative densities have occurred
  !*      mark      Value warn is set to if check fails (for diagnostics)
  !*
  SUBROUTINE CheckFlag(flag,check,warn,mark,*)
    INTEGER flag,check,warn,mark
    IF (flag.EQ.check) THEN
       warn=mark
       PRINT *, 'Flag=',flag
       RETURN 1
    END IF
    RETURN
  END SUBROUTINE CheckFlag

  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Pass fluxes to "last time step fluxes" and set flux arrays to zero.
  SUBROUTINE InitPlas(iLine,DoSavePreviousAndReset)
    USE ModSeGrid, only: nEnergy, nPlas,nIono,nThetaAlt_II
    IMPLICIT NONE

    integer, intent(in) :: iLine
    logical, intent(in) :: DoSavePreviousAndReset
    real :: c
    INTEGER iPlas,j,iAngle,iAlt

    ! Set upper energy boundary condition on the flux. 
    ! Set fluxes at nEnergy+1 to fluxes at nEnergy times a constant <= 1 
    do iPlas=1,nPlas
       iAlt=iPlas+nIono
       do iAngle=0,nThetaAlt_II(iLine,iAlt)
          
          ! set bc for fluxes up
          c=0.
          IF (phiup(iLine,iAngle,iPlas,nEnergy-1).GT.0.) &
               c=min(1., .75*phiup(iLine,iAngle,iPlas,nEnergy)&
               /phiup(iLine,iAngle,iPlas,nEnergy-1))
         
          phiup(iLine,iAngle,iPlas,nEnergy+1)=&
               phiup(iLine,iAngle,iPlas,nEnergy)*c

          ! set bc for fluxes down
          c=0.
          IF (phidn(iLine,iAngle,iPlas,nEnergy-1).GT.0.) &
               c=min(1., .75*phidn(iLine,iAngle,iPlas,nEnergy)&
               /phidn(iLine,iAngle,iPlas,nEnergy-1))
         
          phidn(iLine,iAngle,iPlas,nEnergy+1)=&
               phidn(iLine,iAngle,iPlas,nEnergy)*c
       end do
    end do
    
    ! Only update nEnergy+1 values?
    if (.not.DoSavePreviousAndReset) RETURN

!*  Move current time fluxes to previous time fluxes
    lphiup(iLine,:,1:nPlas,1:nEnergy) = phiup(iLine,:,1:nPlas,1:nEnergy)
    lphidn(iLine,:,1:nPlas,1:nEnergy) = phidn(iLine,:,1:nPlas,1:nEnergy)

!*  Reset current time fluxes
    phiup(iLine,:,:,:)=0.0
    phidn(iLine,:,:,:)=0.0

    RETURN
  END SUBROUTINE InitPlas
  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Pass fluxes to "last time step fluxes" and set flux arrays to zero.
  SUBROUTINE InitPlas_pot(iLine,DoSavePreviousAndReset)
    USE ModSeGrid, only: nEnergy, nPlas,nIono,nThetaAlt_IIC
    IMPLICIT NONE

    integer, intent(in) :: iLine
    logical, intent(in) :: DoSavePreviousAndReset
    real :: c
    INTEGER iPlas,j,iAngle,iAlt

    ! Set upper energy boundary condition on the flux. 
    ! Set fluxes at nEnergy+1 to fluxes at nEnergy times a constant <= 1 
    do iPlas=1,nPlas
       iAlt=iPlas+nIono
       do iAngle=0,nThetaAlt_IIC(iLine,nEnergy,iAlt)
          
          ! set bc for fluxes up
          c=0.
          IF (phiup(iLine,iAngle,iPlas,nEnergy-1).GT.0.) &
               c=min(1., .75*phiup(iLine,iAngle,iPlas,nEnergy)&
               /phiup(iLine,iAngle,iPlas,nEnergy-1))
         
          phiup(iLine,iAngle,iPlas,nEnergy+1)=&
               phiup(iLine,iAngle,iPlas,nEnergy)*c

          ! set bc for fluxes down
          c=0.
          IF (phidn(iLine,iAngle,iPlas,nEnergy-1).GT.0.) &
               c=min(1., .75*phidn(iLine,iAngle,iPlas,nEnergy)&
               /phidn(iLine,iAngle,iPlas,nEnergy-1))
         
          phidn(iLine,iAngle,iPlas,nEnergy+1)=&
               phidn(iLine,iAngle,iPlas,nEnergy)*c
       end do
    end do
    
    ! Only update nEnergy+1 values?
    if (.not.DoSavePreviousAndReset) RETURN

!*  Move current time fluxes to previous time fluxes
    lphiup(iLine,:,1:nPlas,1:nEnergy) = phiup(iLine,:,1:nPlas,1:nEnergy)
    lphidn(iLine,:,1:nPlas,1:nEnergy) = phidn(iLine,:,1:nPlas,1:nEnergy)

!*  Reset current time fluxes
    phiup(iLine,:,:,:)=0.0
    phidn(iLine,:,:,:)=0.0

    RETURN
  END SUBROUTINE InitPlas_pot
  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Pass fluxes to "last time step fluxes" and set flux arrays to zero.
  SUBROUTINE InitIono(iLine,DoSavePreviousAndReset)
    USE ModSeGrid, only: nEnergy, nPlas,nIono,nThetaAlt_II
    IMPLICIT NONE

    integer, intent(in) :: iLine
    logical, intent(in) :: DoSavePreviousAndReset
    real :: c
    INTEGER iIono,j,iAngle,iAlt

    ! Set upper energy boundary condition on the flux. 
    ! Set fluxes at nEnergy+1 to fluxes at nEnergy times a constant <= 1 
    do iIono=1,nIono
       iAlt=iIono
       do iAngle=0,nThetaAlt_II(iLine,iAlt)
          
          ! set bc for fluxes up
          c=0.
          IF (iphiup(iLine,iAngle,iIono,nEnergy-1).GT.0.) c=min(1., &
               .75*iphiup(iLine,iAngle,iIono,nEnergy) &
               /iphiup(iLine,iAngle,iIono,nEnergy-1))
          iphiup(iLine,iAngle,iIono,nEnergy+1)=&
               iphiup(iLine,iAngle,iIono,nEnergy)*c

          ! set bc for fluxes down
          c=0.
          IF (iphidn(iLine,iAngle,iIono,nEnergy-1).GT.0.) c=min(1., &
               .75*iphidn(iLine,iAngle,iIono,nEnergy) &
               /iphidn(iLine,iAngle,iIono,nEnergy-1))
          iphidn(iLine,iAngle,iIono,nEnergy+1)=&
               iphidn(iLine,iAngle,iIono,nEnergy)*c

       end do
    end do
    
    ! Only update nEnergy+1 values?
    if (.not.DoSavePreviousAndReset) RETURN

!*  Move current time fluxes to previous time fluxes
    liphiup(iLine,:,1:2*nIono,1:nEnergy) = iphiup(iLine,:,1:2*nIono,1:nEnergy)
    liphidn(iLine,:,1:2*nIono,1:nEnergy) = iphidn(iLine,:,1:2*nIono,1:nEnergy)

!*  Reset current time fluxes
    iphiup(iLine,:,:,:)=0.0
    iphidn(iLine,:,:,:)=0.0

    RETURN
  END SUBROUTINE InitIono
  !=============================================================================
  !* ------------------------------------------------------------------ **
  !*  Pass fluxes to "last time step fluxes" and set flux arrays to zero.
  SUBROUTINE InitIono_pot(iLine,DoSavePreviousAndReset)
    USE ModSeGrid, only: nEnergy, nPlas,nIono,nThetaAlt_IIC
    IMPLICIT NONE

    integer, intent(in) :: iLine
    logical, intent(in) :: DoSavePreviousAndReset
    real :: c
    INTEGER iIono,j,iAngle,iAlt

    ! Set upper energy boundary condition on the flux. 
    ! Set fluxes at nEnergy+1 to fluxes at nEnergy times a constant <= 1 
    do iIono=1,nIono
       iAlt=iIono
       do iAngle=0,nThetaAlt_IIC(iLine,nEnergy,iAlt)
          ! set bc for fluxes up
          c=0.
          IF (iphiup(iLine,iAngle,iIono,nEnergy-1).GT.0.) c=min(1., &
               .75*iphiup(iLine,iAngle,iIono,nEnergy) &
               /iphiup(iLine,iAngle,iIono,nEnergy-1))
          iphiup(iLine,iAngle,iIono,nEnergy+1)=&
               iphiup(iLine,iAngle,iIono,nEnergy)*c

          ! set bc for fluxes down
          c=0.
          IF (iphidn(iLine,iAngle,iIono,nEnergy-1).GT.0.) c=min(1., &
               .75*iphidn(iLine,iAngle,iIono,nEnergy) &
               /iphidn(iLine,iAngle,iIono,nEnergy-1))
          iphidn(iLine,iAngle,iIono,nEnergy+1)=&
               iphidn(iLine,iAngle,iIono,nEnergy)*c

       end do
    end do
    
    ! Only update nEnergy+1 values?
    if (.not.DoSavePreviousAndReset) RETURN

!*  Move current time fluxes to previous time fluxes
    liphiup(iLine,:,1:2*nIono,1:nEnergy) = iphiup(iLine,:,1:2*nIono,1:nEnergy)
    liphidn(iLine,:,1:2*nIono,1:nEnergy) = iphidn(iLine,:,1:2*nIono,1:nEnergy)

!*  Reset current time fluxes
    iphiup(iLine,:,:,:)=0.0
    iphidn(iLine,:,:,:)=0.0

    RETURN
  END SUBROUTINE InitIono_pot

  !============================================================================
  !  This subroutine checks for "time" convergence for a particular line
  SUBROUTINE check_time_pot(iLine,flag)
    use ModSeGrid, only: nAngle, nPlas, nIono, nEnergy, nPoint, nLine,&
         nThetaAlt_IIC,IsVerbose

    integer,intent(in) :: iLine
    integer,intent(out):: flag
    integer :: iEnergy, iAngle, iAlt, iIono,iPlas
    integer ibad
    real    :: flux, oldflux
    
    real ::   error, ErrorMax, FluxAtMax, OldFluxAtMax
    integer:: iAltAtMax,iAngleAtMAx,iEnergyAtMax
    !--------------------------------------------------------------------------

    Ibad=0
    flag=0
    errormax=0.0
    ! loop over all energy, altitudes, and angles to evaluate convergence
    do iEnergy=nEnergy,1,-1
       do iAlt=1,nPoint
          do iAngle=0,nThetaAlt_IIC(iLine,nEnergy,iAlt)
             !choose what region to set the flux from
             IF (iAlt <= nIono) THEN
                !ionosphere 1
                flux=iphiup(iLine,iAngle,iAlt,iEnergy)
                oldflux=liphiup(iLine,iAngle,iAlt,iEnergy)
             elseif(iAlt >nIono .and. iAlt <=nPoint-nIono) then
                !plasmasphere
                iPlas=iAlt-nIono
                flux=phiup(iLine,iAngle,iPlas,iEnergy)
                oldflux=lphiup(iLine,iAngle,iPlas,iEnergy)
             else
                !ionosphere 2
                iIono=iAlt-nPlas
                flux=iphiup(iLine,iAngle,iIono,iEnergy)
                oldflux=liphiup(iLine,iAngle,iIono,iEnergy)
             endif
             
             ! check the convergence
!             call CheckConv(flux,oldflux,epsilon,flag,DoReportError=.false.,&
!                  ReturnError=error)
             call CheckConv_full(flux,oldflux,epsilon,flag,DoReportError=.false.,&
                  ReturnError=error)
             
             !find maximum error and associated values to report at end of check
             ErrorMax = max(error,errormax)
             ! If a new errormax update the associated iAltAtMax and fluxes
             if (ErrorMax==error) then
                iAltAtMax   = iAlt
                iAngleAtMax = iAngle
                iEnergyAtMax= iEnergy
                FluxAtMax   = flux
                OldFluxAtMax= oldflux
             endif

             if (flag ==1 .and. oldflux > 0.0) Ibad=Ibad+1
             
          end do
       end do
    end do
    !report results of check_time
    if(IsVerbose) then
       write(*,*) 'FINISHED check_time, Report for iLine:',iLine
       write(*,*) 'check_time: Ibad=',Ibad
       write(*,*) 'check_time: ErrorMax=',ErrorMax
       write(*,*) 'check_time: iAltAtMax=',iAltAtMax
       write(*,*) 'check_time: iAngleAtMax=',iAngleAtMax
       write(*,*) 'check_time: nThetaAtMax=',nThetaAlt_IIC(iLine,iEnergyAtMax,iAltAtMax)
       write(*,*) 'check_time: iEnergyAtMax=',iEnergyAtMax
       write(*,*) 'check_time: FluxAtMax=',FluxAtMax
       write(*,*) 'check_time: OldFluxMax=',OldFluxAtMax
    endif
    RETURN
  END SUBROUTINE check_time_pot
  !============================================================================
  !  This subroutine checks for "time" convergence for a particular line
  SUBROUTINE check_time(iLine,flag)
    use ModSeGrid, only: nAngle, nPlas, nIono, nEnergy, nPoint, nLine,&
         nThetaAlt_II,IsVerbose

    integer,intent(in) :: iLine
    integer,intent(out):: flag
    integer :: iEnergy, iAngle, iAlt, iIono,iPlas
    integer ibad
    real    :: flux, oldflux
    
    real ::   error, ErrorMax, FluxAtMax, OldFluxAtMax
    integer:: iAltAtMax,iAngleAtMAx,iEnergyAtMax
    !--------------------------------------------------------------------------

    Ibad=0
    flag=0
    errormax=0.0
    ! loop over all energy, altitudes, and angles to evaluate convergence
    do iEnergy=nEnergy,1,-1
       do iAlt=1,nPoint
          do iAngle=0,nThetaAlt_II(iLine,iAlt)
             !choose what region to set the flux from
             IF (iAlt <= nIono) THEN
                !ionosphere 1
                flux=iphiup(iLine,iAngle,iAlt,iEnergy)
                oldflux=liphiup(iLine,iAngle,iAlt,iEnergy)
             elseif(iAlt >nIono .and. iAlt <=nPoint-nIono) then
                !plasmasphere
                iPlas=iAlt-nIono
                flux=phiup(iLine,iAngle,iPlas,iEnergy)
                oldflux=lphiup(iLine,iAngle,iPlas,iEnergy)
             else
                !ionosphere 2
                iIono=iAlt-nPlas
                flux=iphiup(iLine,iAngle,iIono,iEnergy)
                oldflux=liphiup(iLine,iAngle,iIono,iEnergy)
             endif
             
             ! check the convergence
!             call CheckConv(flux,oldflux,epsilon,flag,DoReportError=.false.,&
!                  ReturnError=error)
             call CheckConv_full(flux,oldflux,epsilon,flag,DoReportError=.false.,&
                  ReturnError=error)
             
             !find maximum error and associated values to report at end of check
             ErrorMax = max(error,errormax)
             ! If a new errormax update the associated iAltAtMax and fluxes
             if (ErrorMax==error) then
                iAltAtMax   = iAlt
                iAngleAtMax = iAngle
                iEnergyAtMax= iEnergy
                FluxAtMax   = flux
                OldFluxAtMax= oldflux
             endif

             if (flag ==1 .and. oldflux > 0.0) Ibad=Ibad+1
             
          end do
       end do
    end do
    !report results of check_time
    if(IsVerbose) then
       write(*,*) 'FINISHED check_time, Report:'
       write(*,*) 'check_time: Ibad=',Ibad
       write(*,*) 'check_time: ErrorMax=',ErrorMax
       write(*,*) 'check_time: iAltAtMax=',iAltAtMax
       write(*,*) 'check_time: iAngleAtMax=',iAngleAtMax
       write(*,*) 'check_time: nThetaAtMax=',nThetaAlt_II(iLine,iAltAtMax)
       write(*,*) 'check_time: iEnergyAtMax=',iEnergyAtMax
       write(*,*) 'check_time: FluxAtMax=',FluxAtMax
       write(*,*) 'check_time: OldFluxMax=',OldFluxAtMax
    endif
    RETURN
  END SUBROUTINE check_time
  !============================================================================
  ! A subroutine to get integrated output quantities including number density,
  ! number flux, and volume heating rate
  subroutine calc_integrated_output(iLine,nNeutral,eThermalDensity_C)
    use ModSeGrid, only: nAngle, nPlas, nIono, nEnergy, nPoint,&
         KineticEnergy_IIC, DoIncludePotential,nThetaAlt_II,nThetaAlt_IIC, &
         mu_III,mu_IIIC, EnergyGrid_I, DeltaE_I,MaxAlt_IC,UsePwRegion,nPwRegion
    use ModMath, only: midpnt_int
    use ModNumConst,    ONLY: cPi
    implicit none
    
    integer, intent(in) :: iLine,nNeutral
    real   , intent(in) :: eThermalDensity_C(nPoint)
    real,allocatable :: spec(:)   !integrand of heating rate
    real,allocatable :: NumDensIntegrand_I(:) !integrand of number density
    real,allocatable :: IntegrandFluxUp_I(:) !integrand of flux
    real,allocatable :: IntegrandFluxDn_I(:) !integrand of flux
    real,allocatable :: FluxUp_IC(:,:) !fluxUp for each energy and alt
    real,allocatable :: FluxDn_IC(:,:) !fluxDn for each energy and alt
    real,allocatable :: NetFlux_IC(:,:) !net flux for each energy and alt
    real,allocatable :: delKE_I(:)!energy spacing of KE array at a given totalE
    real, parameter  :: Acoef = 2.6e-12 !ev^2cm^2

    !loop variables
    integer :: iPoint, iEnergy, iAngle, iPlas, iIono, iNeutral
    !---------------------------------------------------------------------------
    !allocate temporary arrays
    if (.not.allocated(spec)) allocate(spec(nEnergy))
    if (.not.allocated(NumDensIntegrand_I))allocate(NumDensIntegrand_I(nEnergy))
    if (.not.allocated(delKE_I)) allocate(delKE_I(nEnergy))

    !allocate and initialize integrand arrays
    if (.not.allocated(IntegrandFluxUp_I))allocate(IntegrandFluxUp_I(nAngle))
    if (.not.allocated(IntegrandFluxDn_I))allocate(IntegrandFluxDn_I(nAngle))
    if (.not.allocated(FluxUp_IC))        allocate(FluxUp_IC(nEnergy,nPoint))
    if (.not.allocated(FluxDn_IC))        allocate(FluxDn_IC(nEnergy,nPoint))
    if (.not.allocated(NetFlux_IC))       allocate(NetFlux_IC(nEnergy,nPoint))
    IntegrandFluxUp_I(:)=0.0
    IntegrandFluxDn_I(:)=0.0
    FluxUp_IC(:,:)=0.0
    FluxDn_IC(:,:)=0.0
    NetFlux_IC(:,:)=0.0
    NumberFlux_IC(iLine,:)=0.0

    HeatingRate_IC(iLine,:) = 0.0
    NumberDens_IC(iLine,:)=0.0
    
    ALONG_LINE: do iPoint=1,nPoint
       if (UsePwRegion .and. iPoint > nIono+nPwRegion) EXIT ALONG_LINE
       !\
       ! calculate the volume heating rate
       !/
       ENERGY: do iEnergy=1,nEnergy
          if (DoIncludePotential) then
             !when above maxalt for a given energy integrand is zero and cycle
             if (iPoint>MaxAlt_IC(iLine,iEnergy) &
                  .and. iPoint<=nPoint-MaxAlt_IC(iLine,iEnergy)) then
                Spec(iEnergy)=0.0
                cycle ENERGY
             endif

             if(KineticEnergy_IIC(iLine,iEnergy,iPoint) < 1e-30)then
                Spec(iEnergy) = 0.0
             else
                Spec(iEnergy)=&
                     (specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint))&
                     / KineticEnergy_IIC(iLine,iEnergy,iPoint)
             end if
             !set the delta kinetic energy array
             if (iEnergy<nEnergy) then
                delKE_I(iEnergy) = KineticEnergy_IIC(iLine,iEnergy+1,iPoint) &
                     - KineticEnergy_IIC(iLine,iEnergy,iPoint)
             else
                delKE_I(iEnergy) = delKE_I(iEnergy-1)
             endif
          else
             Spec(iEnergy)=&
                  (specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint))&
                  / EnergyGrid_I(iEnergy)
             delKE_I(iEnergy)=DeltaE_I(iEnergy)
          endif
       end do ENERGY
       CALL midpnt_int(HeatingRate_IC(iLine,iPoint),Spec,delKE_I,1,nEnergy,&
            nEnergy,2)
       HeatingRate_IC(iLine,iPoint)=&
            4.*cPi*Acoef*eThermalDensity_C(iPoint)*HeatingRate_IC(iLine,iPoint)
       !\
       ! calculate the total electron production rate for each altitude
       !/
       if (iPoint <= nIono .or. iPoint > nPoint-nIono)then 
          ! set iIono index
          if (iPoint <= nIono) then
             iIono = iPoint
          else
             iIono = iPoint-nIono-nPlas
          endif

          DO iNeutral=1,nNeutral
             CALL midpnt_int(SecondaryIonRate_IIC(iLine,iNeutral,iPoint),&
                  Qestar_IICI(iLine,iNeutral,iIono,:),delKE_I,1,&
                  nEnergy,nEnergy,2)
             SecondaryIonRate_IIC(iLine,iNeutral,iPoint) = &
                  4.0*cPi*SecondaryIonRate_IIC(iLine,iNeutral,iPoint)
          END DO
          
          CALL midpnt_int(TotalIonizationRate_IC(iLine,iPoint),&
               Qstar_ICI(iLine,iIono,:),delKE_I,1,nEnergy,nEnergy,2)
          TotalIonizationRate_IC(iLine,iPoint) = &
               4.0*cPi*TotalIonizationRate_IC(iLine,iPoint)
       else
          TotalIonizationRate_IC(iLine,iPoint)=0.0
       endif
       
       !\
       ! calculate number density
       !/
       ENERGY2: do iEnergy=1,nEnergy
          if (DoIncludePotential) then
             !when above maxalt for a given energy integrand is zero and cycle
             if (iPoint>MaxAlt_IC(iLine,iEnergy) &
                  .and. iPoint<=nPoint-MaxAlt_IC(iLine,iEnergy)) then
                NumDensIntegrand_I(iEnergy) = 0.0
                cycle ENERGY2
             endif
             NumDensIntegrand_I(iEnergy) = &
                  (specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint))&
                  /sqrt(KineticEnergy_IIC(iLine,iEnergy,iPoint))
          else
             NumDensIntegrand_I(iEnergy) = &
                  (specup(iLine,iEnergy,iPoint)-specdn(iLine,iEnergy,iPoint))&
                  /sqrt(EnergyGrid_I(iEnergy))
          endif
       enddo ENERGY2
       CALL midpnt_int(NumberDens_IC(iLine,iPoint),NumDensIntegrand_I,&
            delKE_I,1,nEnergy,nEnergy,2)
       NumberDens_IC(iLine,iPoint)=4.*cPi*1.7E-8*NumberDens_IC(iLine,iPoint)
       
       !\
       ! calculate the number flux
       !/
       ! start by finding the netfluxfor each energy
       ! for each altitude find the net flux for each energy, then integrate 
       ! over energy

       !set up the integrand
       ENERGY3: do iEnergy=1,nEnergy
          !Get the up flux 
          if (.not.DoIncludePotential) then
             ! when not including potential use mu_III and appropriate theta 
             ! range
             ! for each energy integrandflux=int(mu0*phi dmu) from 0 to mu0max
             do iAngle=1,nThetaAlt_II(iLine,iPoint)
                ! Choose iphi or phi based on position on field line
                if (iPoint <= nIono .or. iPoint > nPoint-nIono)then 
                   ! set iIono index
                   if (iPoint <= nIono) then
                      iIono = iPoint
                   else
                      iIono = iPoint-nIono-nPlas
                   endif
                   ! Integral is +
                   IntegrandFluxUp_I(iAngle)=-mu_III(iLine,iAngle,iPoint)&
                        * iphiup(iLine,iAngle,iIono,iEnergy)
                   ! Integral is -
                   IntegrandFluxDn_I(iAngle)=mu_III(iLine,iAngle,iPoint)&
                        * iphidn(iLine,iAngle,iIono,iEnergy)
                else
                   ! set iPlas index
                   iPlas = iPoint-nIono
                   ! Integral is +
                   IntegrandFluxUp_I(iAngle)=-mu_III(iLine,iAngle,iPoint)&
                        * phiup(iLine,iAngle,iPlas,iEnergy)
                   ! Integral is -
                   IntegrandFluxDn_I(iAngle)=mu_III(iLine,iAngle,iPoint)&
                        * phidn(iLine,iAngle,iPlas,iEnergy)
                endif
             end do
             ! integrate to get the up and down flux
             CALL midpnt_int(FluxUp_IC(iEnergy,iPoint),IntegrandFluxUp_I,&
                  mu_III(iLine,:,iPoint),1,nThetaAlt_II(iLine,iPoint),&
                  nAngle,1)
             CALL midpnt_int(FluxDn_IC(iEnergy,iPoint),IntegrandFluxDn_I,&
                  mu_III(iLine,:,iPoint),1,nThetaAlt_II(iLine,iPoint),&
                  nAngle,1)
          else
             ! when including potential use mu_IIIC and appropriate theta 
             ! range
             ! for each energy integrandflux=int(mu0*phi dmu) from 0 to mu0max
             
             !when above maxalt for a given energy integrand is zero and cycle
             if (iPoint>MaxAlt_IC(iLine,iEnergy) &
                  .and. iPoint<=nPoint-MaxAlt_IC(iLine,iEnergy)) then
                NetFlux_IC(iEnergy,iPoint) = 0.0
                cycle ENERGY3
             endif
             
             do iAngle=1,nThetaAlt_IIC(iLine,iEnergy,iPoint)
                ! Choose iphi or phi based on position on field line
                if (iPoint <= nIono .or. iPoint > nPoint-nIono)then 
                   ! set iIono index
                   if (iPoint <= nIono) then
                      iIono = iPoint
                   else
                      iIono = iPoint-nIono-nPlas
                   endif
                   ! Integral is +
                   IntegrandFluxUp_I(iAngle)=&
                        -mu_IIIC(iLine,iAngle,iEnergy,iPoint)&
                        * iphiup(iLine,iAngle,iIono,iEnergy)
                   ! Integral is -
                   IntegrandFluxDn_I(iAngle)=&
                        mu_IIIC(iLine,iAngle,iEnergy,iPoint)&
                        * iphidn(iLine,iAngle,iIono,iEnergy)
                else
                   ! set iPlas index
                   iPlas = iPoint-nIono
                   ! Integral is +
                   IntegrandFluxUp_I(iAngle)=&
                        -mu_IIIC(iLine,iAngle,iEnergy,iPoint)&
                        * phiup(iLine,iAngle,iPlas,iEnergy)
                   ! Integral is -
                   IntegrandFluxDn_I(iAngle)=&
                        mu_IIIC(iLine,iAngle,iEnergy,iPoint)&
                        * phidn(iLine,iAngle,iPlas,iEnergy)
                endif
             enddo
             ! integrate to get the up and down flux
             CALL midpnt_int(FluxUp_IC(iEnergy,iPoint),IntegrandFluxUp_I,&
                  mu_IIIC(iLine,:,iEnergy,iPoint),1,&
                  nThetaAlt_IIC(iLine,iEnergy,iPoint),nAngle,1)
             CALL midpnt_int(FluxDn_IC(iEnergy,iPoint),IntegrandFluxDn_I,&
                  mu_IIIC(iLine,:,iEnergy,iPoint),1,&
                  nThetaAlt_IIC(iLine,iEnergy,iPoint),nAngle,1)
          endif
          
          ! Get the net flux
          NetFlux_IC(iEnergy,iPoint)=&
               2.0*cPi*(FluxUp_IC(iEnergy,iPoint)+FluxDn_IC(iEnergy,iPoint))
       end do ENERGY3
       
       ! calculate the number flux
       CALL midpnt_int(NumberFlux_IC(iLine,iPoint),NetFlux_IC(:,iPoint),&
            delKE_I,1,nEnergy,nEnergy,2)
    end do ALONG_LINE
    
    
    ! deallocate to save memory
    deallocate(spec,delKE_I,NumDensIntegrand_I,IntegrandFluxUp_I,&
         IntegrandFluxDn_I,NetFlux_IC)
  end subroutine calc_integrated_output
  !============================================================================
  subroutine allocate_state_arrays
    use ModSeGrid, only: nAngle, nPlas, nIono, nEnergy, nPoint, nLine
    use ModSeBackground, only: nNeutralSpecies
    
    if(.not.allocated(phiup)) allocate(phiup(nLine,0:nAngle,0:nPlas+1,nEnergy+1))
    if(.not.allocated(phidn)) allocate(phidn(nLine,0:nAngle,0:nPlas+1,nEnergy+1))

    if(.not.allocated(iphiup)) &
         allocate(iphiup(nLine,0:nAngle,0:2*nIono+1,nEnergy+1))
    if(.not.allocated(iphidn)) &
         allocate(iphidn(nLine,0:nAngle,0:2*nIono+1,nEnergy+1))

    if(.not.allocated(specup)) allocate(specup(nLine,nEnergy,nPoint))
    if(.not.allocated(specdn)) allocate(specdn(nLine,nEnergy,nPoint))

    if(.not.allocated(lphiup)) allocate(lphiup(nLine,0:nAngle,nPlas,nEnergy))
    if(.not.allocated(lphidn)) allocate(lphidn(nLine,0:nAngle,nPlas,nEnergy))

    if(.not.allocated(liphiup))allocate(liphiup(nLine,0:nAngle,2*nIono,nEnergy))
    if(.not.allocated(liphidn))allocate(liphidn(nLine,0:nAngle,2*nIono,nEnergy))

    if(.not.allocated(SRC))    allocate(SRC(nPoint))

    if(.not.allocated(Qestar_IICI)) &
         allocate(Qestar_IICI(nLine,nNeutralSpecies,2*nIono,nEnergy))
    if(.not.allocated(Qpstar_ICI))allocate(Qpstar_ICI(nLine,2*nIono,nEnergy))
    if(.not.allocated(Qstar_ICI)) allocate(Qstar_ICI (nLine,2*nIono,nEnergy))
    ! Initiallize Qpstar to 0 as we do not include this will get updated 
    ! if a separate primary beam is included
    Qpstar_ICI(:,:,:)=0

    !allocate arrays to hold integrated output
    if(.not.allocated(HeatingRate_IC)) allocate(HeatingRate_IC(nLine,nPoint))
    if(.not.allocated(NumberDens_IC))  allocate(NumberDens_IC(nLine,nPoint))
    if(.not.allocated(NumberFlux_IC))  allocate(NumberFlux_IC(nLine,nPoint))
    if(.not.allocated(SecondaryIonRate_IIC))&
         allocate(SecondaryIonRate_IIC(nLine,nNeutralSpecies,nPoint))
    if(.not.allocated(TotalIonizationRate_IC))&
         allocate(TotalIonizationRate_IC(nLine,nPoint))

    ! Initializ to zero in case only a part of the line is used (with PWOM)
    specup                =0
    specdn                =0
    HeatingRate_IC        =0
    NumberDens_IC         =0
    NumberFlux_IC         =0
    SecondaryIonRate_IIC  =0
    TotalIonizationRate_IC=0
    
  end subroutine allocate_state_arrays

end Module ModSeState
