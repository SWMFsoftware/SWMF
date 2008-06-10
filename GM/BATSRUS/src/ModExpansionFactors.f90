!^CFG COPYRIGHT UM 
module ModExpansionFactors
  use ModMpi
  use ModIoUnit,   ONLY: io_unit_new
  use ModMagnetogram
  use ModConst

  !Dependecies to be removed
  use ModProcMH,   ONLY: iProc,nProc,iComm
  use ModIO, ONLY: iUnitOut, write_prefix
  implicit none
  save

  !Gravity potential, m^2/s^2
  real,parameter :: cSunGravitySI=cGravitation*mSun/Rsun
  real :: T0=3.5E+6 !in K

  !Gravity potential of a proton, in K
  real,parameter :: cSunGravityK =cSunGravitySI*cProtonMass/cBoltzmann

  character(len=20) :: TypeModel='WSA' ! Type of SW model
  !Distribution of the solar wind model parameters: 

  real,allocatable,dimension(:,:,:) :: FiskFactor_N
  ! The value of this factor at a given grid point
  ! is equal to the value of |B_R(r=r_Sun)|/|B(r=R_Sun)|}, 
  ! where B_R is taken at the "photospheric footpoint" of the 
  ! magnetic field line, passing through this point:
  !\               Grid point
  ! \R_Sun       iR,iPhi,iTheta   !R_surface
  !  -------------+---------------!
  ! /                             !
  !/  
  ! Field line predicted by the source surface model 
  ! (this is a real curved magnetic field, not a ray Theta=const
  ! The value of the Fisk factor at the considered grid point is 
  ! defined as
  ! 
  ! FiskFactor_N(iR,iPhi,iTheta)=|B_R(R=R_Sun)|/|B|, 
  !
  ! if the magnetic field line is open, 
  ! zero otherwise.

  real,allocatable,dimension(:,:,:) :: ExpansionFactorInv_N
  ! The expansion factor. !!!INVERTED!!!!
  ! The value of this factor at a given grid point
  ! is equal to the value of  
  ! B(R=R_{SourceSurface}/B(R=R_{Sun})*(R_SourceSurface/R_Sun)^2
  ! where the ratio of the magnetic field intensities is taken at the 
  ! two marginal point of the magnetic field line, passing through 
  ! the considered grid point:
  !
  !\               Grid point
  ! \R_Sun       iR,iPhi,iTheta   !R_SourceSurface
  !  -------------+---------------!
  ! /                             !
  !/ 
  ! Field line predicted by the source surface model 
  ! (this is a real curved magnetic field, not a ray Theta=const.
  ! The definition of the expansion factor is taken from:
  ! Wang & Shheley, ApJ, 1990, 355:726 and from 
  ! Arge & Pizzo, JGR, 2000, 105(A5):10,465    
  ! We use the referenced definition to define the expansion factor 
  ! for open field lines. If the field line is close, the expansion
  ! factor is set to be zero.
  real,allocatable,dimension(:,:,:) :: ThetaB_N

  ! Speed distribution extracted from Wang-Sheeley-Arge 
  ! model (Arge et al. 2006):
  real,allocatable,dimension(:,:,:) :: WSAspeed_N

  ! Speed distribution extracted from Fisk model
  real,allocatable,dimension(:,:,:) :: Fiskspeed_N
  real::UMin=265.0
  real::GammaSS=1.10
contains
  subroutine set_expansion_factors
    real :: dS,dSMax
    real,dimension(nDim) :: R_D !The vector r,phi,theta
    real,dimension(nDim) :: BTemp_D,BSun_D,BSS_D
    real,dimension(nDim) :: RSS_D,RSun_D,RPlusEnd_D,RMinusEnd_D
    integer :: iR,iPhi,iTheta
    integer :: iBcast, iStart, nSize, iError,iIteration
    real,allocatable,dimension(:,:) :: Phi_IJ,Theta_IJ
    real,parameter :: WSAPowerIndex=2.0/7.0
    real,parameter :: cLoopLifetime=8.0*cSecondPerHour
    real,parameter :: cFiskQ=1.64E+16 !m^2/s

    ! Allocte factors arrays
    if(allocated(ExpansionFactorInv_N))deallocate(ExpansionFactorInv_N)
    allocate(ExpansionFactorInv_N(-nRExt:nR,0:nPhi,0:nTheta))
    if(allocated(FiskFactor_N))deallocate(FiskFactor_N)
    allocate(FiskFactor_N(-nRExt:nR,0:nPhi,0:nTheta))
    if(allocated(ThetaB_N))deallocate(ThetaB_N)
    allocate(ThetaB_N(-nRExt:nR,0:nPhi,0:nTheta))

    !Initalize arrays:
    ExpansionFactorInv_N=cZero
    FiskFactor_N=cZero
    ThetaB_N=cZero

    if(allocated(Phi_IJ))deallocate(Phi_IJ)
    allocate(Phi_IJ(0:nPhi,0:nTheta))
    if(allocated(Theta_IJ))deallocate(Theta_IJ)
    allocate(Theta_IJ(0:nPhi,0:nTheta))

    Phi_IJ=cZero
    Theta_IJ=cZero
    do iTheta=0,nTheta
       do iPhi=0,nPhi
          Phi_IJ(iPhi,iTheta)=iPhi*dPhi
          Theta_IJ(iPhi,iTheta)=colatitude(iTheta)
       end do
    end do


    dSMax=cHalf*(Rs_PFSSM-Ro_PFSSM)

    !Loop by theta, each processor treats a separate part of the grid
    do iTheta=iProc*nThetaPerProc,(iProc+1)*nThetaPerProc-1 
       if(iTheta>nTheta)EXIT !Some processors may have less amount of
       ! work
       do iPhi=0,nPhi
          do iR=-nRExt,nR
             ! Define the location of the grid point
             call start_at_grid_point(iR,iPhi,iTheta)
             iIteration=0
             !Integrate in the positive direction of the magnetic
             ! field
             do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                !Integrate the solution per local dS
                iIteration=iIteration+1
                call advance_line_point(R_D,+1.0)
             end do

             ! Save the positive end of the field line
             RPlusEnd_D = R_D

             ! Define the location of the grid point
             call start_at_grid_point(iR,iPhi,iTheta)
             iIteration=0
             !Integrate in the negative direction of the magnetic
             ! field
             do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                !Integrate the solution per local dS
                iIteration=iIteration+1
                call advance_line_point(R_D,-1.0)             
             end do
             ! Save the negative end of the field line
             RMinusEnd_D = R_D


             ! Check if the field line end points are at the same 
             ! radius (closed or hanging field line). If the field is 
             ! closed, the inv_expansion factor and Fisk factor are 
             ! set to zero. 
             if(abs(RPlusEnd_D(R_)-RMinusEnd_D(R_)) <= dSMax)then
                ExpansionFactorInv_N(iR,iPhi,iTheta) = cZero
                FiskFactor_N(iR,iPhi,iTheta) = cZero
                ! 
             else
                ! Check which end of the field line is at the
                ! photosphere
                ! and which one is at the source surface. Then get
                ! the field
                ! components for both ends and calculate the value of 
                ! the factors.
                if(RPlusEnd_D(R_) > RMinusEnd_D(R_))then
                   RSS_D  = RPlusEnd_D
                   RSun_D = RMinusEnd_D
                else
                   RSS_D  = RMinusEnd_D
                   RSun_D = RPlusEnd_D
                end if
                call interpolate_field(RSS_D,BSS_D)
                call interpolate_field(RSun_D,BSun_D)
                ! Get factors for the grid point
                ExpansionFactorInv_N(iR,iPhi,iTheta)=&
                     & (sqrt(dot_product(BSS_D,BSS_D)) &
                     &/sqrt(dot_product(BSun_D,BSun_D)))* (Rs_PFSSM&
                     &/Ro_PFSSM)**2
                FiskFactor_N(iR,iPhi,iTheta) = abs(BSun_D(R_))/&
                     & sqrt(dot_product(BSun_D,BSun_D))
             end if
          end do
       end do
    end do

    if(nProc>1)then
       do iBcast=0,nProc-1
          iStart=iBcast*nThetaPerProc
          if(iStart>nTheta)EXIT
          nSize=min(nThetaPerProc,nTheta+1-iStart)*(nPhi+1)* (nR+1&
               &+nRExt)
          call MPI_bcast(ExpansionFactorInv_N(-nRExt,0,iStart) ,nSize&
               &,MPI_REAL,iBcast,iComm,iError)
          call MPI_bcast(FiskFactor_N(-nRExt,0,iStart) ,nSize&
               &,MPI_REAL,iBcast,iComm,iError)
       end do
    end if
    ! Calculate ThetaB_N
    do iTheta=iProc*nThetaPerProc,(iProc+1)*nThetaPerProc-1
       if(iTheta>nTheta)EXIT !Some processors may have less amount of
       ! work
       do iPhi=0,nPhi
          do iR=-nRExt,nR
             ! Define the location of the grid point
             call start_at_grid_point(iR,iPhi,iTheta)
             iIteration=0
             !Integrate in the positive direction of the magnetic
             ! field
             do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                !Integrate the solution per local dS
                iIteration=iIteration+1
                call advance_line_point(R_D,+1.0)
             end do
             ! Save the positive end of the field line
             RPlusEnd_D = R_D

             call start_at_grid_point(iR,iPhi,iTheta)
             iIteration=0
             !Integrate in the negative direction of the magnetic
             ! field
             do while (R_D(R_) <= Rs_PFSSM .and. R_D(R_) >= Ro_PFSSM)
                !Integrate the solution per local dS
                iIteration=iIteration+1
                call advance_line_point(R_D,-1.0)
             end do
             ! Save the negative end of the field line
             RMinusEnd_D = R_D


             ! Check if the field line end points are at the same 
             ! radius (closed or hanging field line). If the field is 
             ! closed, the inv_expansion factor and Fisk factor are 
             ! set to zero.                
             if(abs(RPlusEnd_D(R_)-RMinusEnd_D(R_)) <= dSMax)then
                ThetaB_N(iR,iPhi,iTheta) = cZero
             else
                ! Check which end of the field line is at the
                ! photosphere
                ! and which one is at the source surface. Then get
                ! the field
                ! components for both ends and calculate the value of 
                ! the factors.
                if(RPlusEnd_D(R_) > RMinusEnd_D(R_))then
                   RSS_D  = RPlusEnd_D
                   RSun_D = RMinusEnd_D
                else
                   RSS_D  = RMinusEnd_D
                   RSun_D = RPlusEnd_D
                end if
                ! Get ThetaB_N for the grid point
                ThetaB_N(iR,iPhi,iTheta)= theta_b(RSun_D(Phi_)&
                     &,RSun_D(Theta_))
             end if
          end do
       end do
    end do
    if(nProc>1)then
       do iBcast=0,nProc-1
          iStart=iBcast*nThetaPerProc
          if(iStart>nTheta)EXIT
          nSize=min(nThetaPerProc,nTheta+1-iStart)*(nPhi+1)* (nR+1&
               &+nRExt)
          call MPI_bcast(ThetaB_N(-nRExt,0,iStart) ,nSize,MPI_REAL&
               &,iBcast,iComm,iError)
       end do
    end if
    !Transform to Deg
    ThetaB_N=ThetaB_N*cRadToDeg

    ! Get WSA speed
    if(allocated(WSAspeed_N))deallocate(WSAspeed_N)
    allocate(WSAspeed_N(-nRExt:nR,0:nPhi,0:nTheta)) 
    WSAspeed_N=cZero

    ! Calculate WSA speed distribution using eq. 1 in Arge et al.
    ! 2004:
    !WSAspeed_N(:,:,:)=(265.0+&
    !     1.5*ExpansionFactorInv_N(:,:,:)**(1.0/3.0)/&
    !     ( cOne+ExpansionFactorInv_N(:,:,:) )**(1.0/3.0)* &
    !     (5.9-1.5*exp( 1.0-(ThetaB_N(:,:,:)/7.0)**(5.0/2.0) ) &
    !     )**(7.0/2.0) ) &    !km/s so far
    !     *cE3                 !To get the result in SI
    
    ! Calculate WSA speed distribution using eq. obtained
    ! by personal communication with N. Arge (2006)                                              
    WSAspeed_N(:,:,:)=(240.0+&
         (675.0*ExpansionFactorInv_N(:,:,:)**(1./4.5)/(&
         ExpansionFactorInv_N(:,:,:)+1.0)**(1./4.5))*&
         (1.0-0.8*exp(1.0-(ThetaB_N(:,:,:)/2.8)**1.25)/exp(1.0))**3.0)*cE3

    ! Calculate WSA speed distribution using eq. 2 in Arge et al.
    ! 2003:
    !WSAspeed_N(:,:,:)=(265.0+25.0* exp(log(ExpansionFactorInv_N(:,:&
    !     &,:)+cTiny)*WSAPowerIndex)* (5.0-1.1*exp(1.0-(ThetaB_N(:,:&
    !     &,:)/4.0)**2))**2)& !km/s so far
    !     *cE3         !To get the result in SI

    ! Get Fisk speed
    if(allocated(Fiskspeed_N))deallocate(Fiskspeed_N)
    allocate(Fiskspeed_N(-nRExt:nR,0:nPhi,0:nTheta)) 
    Fiskspeed_N=cZero

    ! Calculate Fisk final speed using the eq.:
    ! u_f=sqrt(2*(Q-G))
    ! where Q=4.55555e+11/T and G=g*MSun/RSun=1.9e+11 m^2/s^2
    ! The temperature in the loops, T, being  T=0.8/Fisk_factor
    ! in million degrees K.
    !The speed is limited to be greater than 265 km/s

    Fiskspeed_N(:,:,:)=100.0!FiskFactor_N(:,:,:)

    !Fiskspeed_N(:,:,:)=sqrt(max(2.0* (cFiskQ*& !m^2/s 
    !     max(FiskFactor_N(:,:,:)**2,cHalf**2)/cLoopTemp &
    !     &-cSunGravitySI),(265.0*cE3)**2))

    ! Finding the minimum value of the final speed
    select case(TypeModel) 
    case('WSA')
       UMin=minval(WSAspeed_N)
    case('Fisk')
       UMin=minval(Fiskspeed_N)
    end select

    ! Finding the maximum surface value of gamma (related to the minimum speed)
    gammaSS=( (cHalf*UMin**2+cSunGravitySI)/(T0*cBoltzmann/cProtonMass) ) &
         /( (cHalf*UMin**2+cSunGravitySI)/(T0*cBoltzmann/cProtonMass)-cOne )
  contains
    !----------------------------------------------------------------
    !----------
    subroutine advance_line_point(RInOut_D,Dir)
      real,intent(inout),dimension(nDim) :: RInOut_D
      real,intent(in) :: Dir
      dS=0.25*min(dR,dPhi,dSinTheta,cOne)*cTwo**(iIteration/ (20&
           &*max(nR,nPhi,nTheta)))
      !To avoid the line bouncing near null points
      RInOut_D=RInOut_D+Dir*dS*f_d( RInOut_D+Dir*dS*cHalf&
           &*f_d(RInOut_D))
      call correct_angles(RInOut_D)
    end subroutine advance_line_point
    !----------------------------------------------------------------
    !--
    subroutine start_at_grid_point(iR,iPhi,iTheta)
      integer,intent(in) :: iR,iPhi,iTheta
      R_D(R_)=Ro_PFSSM+real(iR)*dR
      R_D(Phi_)=real(iPhi)*dPhi
      R_D(Theta_)=colatitude(iTheta)
      R_D(R_)=min(max(R_D(R_),Ro_PFSSM+cQuarter*dR),Rs_PFSSM-cQuarter&
           &*dR)
    end subroutine start_at_grid_point
    ! This fucnction calculates the value of 
    ! F(i)= B(i)/|B|/(1,r*sin(colatitude),r)
    function f_d(RIn_D)
      real,dimension(nDim) :: f_d
      real,dimension(nDim),intent(in) :: RIn_D
      real,parameter :: cTol=cOne/(cE9*cE1)

      !Get the vector (B_r,B_phi,B_theta)
      call interpolate_field(RIn_D,f_d)

      !Divide by the metric coefficients, to obtain
      !the vector ||B|| d (r,phi,theta)/dS along the field line

      f_d=f_d/(/cOne,RIn_D(R_)*max(sin(RIn_D(Theta_)),cTol),&
           & RIn_D(R_)/)

      !Divide by some scale, to limit the displacement within the
      ! integration 
      !step
      f_d=f_d/sqrt(sum(f_d**2))
    end function f_d

    function theta_b(Phi,Theta)
      real,intent(in) :: Phi,Theta
      real :: theta_b

      theta_b=sqrt(minval((Phi-Phi_IJ(:,:))**2+ (Theta-Theta_IJ(:,:))&
           &**2, mask=ExpansionFactorInv_N(0,:,:)<0.001))

    end function theta_b

  end subroutine set_expansion_factors
end module ModExpansionFactors
!=================================set_empirical_model=============
subroutine set_empirical_model(TypeRead,BodyT0)
  use ModExpansionFactors
  implicit none

  character(LEN=*),intent(in) :: TypeRead
  real, intent(in) :: BodyT0
  !------------------------------------------------------------------
  TypeModel=trim(TypeRead)
     T0 = BodyT0
     call set_expansion_factors
     if(iProc==0)call write_expansion_tec

   end subroutine set_empirical_model
   !====================================================================
   subroutine get_bernoulli_integral(xInput,yInput,zInput,Output)
     use ModExpansionFactors
     implicit none
     real, intent(in)  :: xInput,yInput,zInput
     real, intent(out) :: Output
     real :: Rin_PFSSM,Theta_PFSSM,Phi_PFSSM

     integer :: Node_D(nDim)
     real :: Res_D(nDim)

     real :: Weight_III(0:1,0:1,0:1)
     real :: R_PFSSM

     !------------------------------------------------------------------
     !\
     ! Calculate cell-centered spherical coordinates::
     !/
     Rin_PFSSM   = sqrt(xInput**2+yInput**2+zInput**2)
     !\
     ! Avoid calculating inside a critical radius = 0.5*Rsun
     !/
     if (Rin_PFSSM <max(Ro_PFSSM-dR*nRExt,0.90*Ro_PFSSM)) then
        Output= cZero
        RETURN
     end if
     Theta_PFSSM = acos(zInput/Rin_PFSSM)
     Phi_PFSSM   = atan2(yInput,xInput)

     !\
     ! Set the source surface radius::
     ! The inner boundary in the simulations starts at a height
     ! H_PFSSM above that of the magnetic field measurements!
     !/

     R_PFSSM =min(Rin_PFSSM+H_PFSSM, Rs_PFSSM)


     !\
     ! Transform Phi_PFSSM from the component's frame to the
     ! magnetogram's frame.
     !/

     Phi_PFSSM = Phi_PFSSM - Phi_Shift*cDegToRad

     !\
     ! Take a residual for the bi-linear interpolation
     !/
     Res_D=(/R_PFSSM,Phi_PFSSM,Theta_PFSSM/)

     !Limit a value of R:
     Res_D(R_)=max(min(Res_D(R_),Rs_PFSSM-cTiny),Ro_PFSSM-nRExt*dR+cTiny)

     Res_D(R_)=Res_D(R_)-Ro_PFSSM

     call correct_angles(Res_D)
     Res_D(Theta_)=cos(Res_D(Theta_)) &!This is sin(latitude)
          -sin_latitude(0)     !the same for the iTheta=0 node
     ! of the grid
     Res_D=Res_D*dInv_D
     Node_D=floor(Res_D)
     if(Node_D(R_)==nR)Node_D(R_)=Node_D(R_)-1
     Res_D=Res_D-real(Node_D)
     if(Node_D(Phi_)==nPhi)Node_D(Phi_)=0

     if(Node_D(Theta_)>=nTheta)then
        Node_D(Theta_)=nTheta-1
        Res_D(Theta_)=cOne
     elseif(Node_D(Theta_)<=-1)then
        Node_D(Theta_)=0
        Res_D(Theta_)=cZero
     end if

     Weight_III(0,:,:)=cOne-Res_D(R_)
     Weight_III(1,:,:)=Res_D(R_)
     Weight_III(:,0,:)=Weight_III(:,0,:)*(cOne-Res_D(Phi_))
     Weight_III(:,1,:)=Weight_III(:,1,:)*Res_D(Phi_)
     Weight_III(:,:,0)=Weight_III(:,:,0)*(cOne-Res_D(Theta_))
     Weight_III(:,:,1)=Weight_III(:,:,1)*Res_D(Theta_)

     Output= sum(Weight_III*WSAspeed_N( Node_D(R_):Node_D(R_)+1,&
          & Node_D(Phi_):Node_D(Phi_)+1,&
          & Node_D(Theta_):Node_D(Theta_)+1))

   end subroutine get_bernoulli_integral

   !==========================================================================

   subroutine get_gamma_emp(xx,yy,zz,gammaOut)

     ! Subroutine get_gamma_emp
     ! Provides the distribution of the polytropic index, complying with
     ! the WSA or Fisk semi-empirical models

     use ModExpansionFactors
     use ModNumConst
     implicit none

     real, intent(in) :: xx,yy,zz
     real, intent(out)   :: gammaOut 
     real :: RR,Uf,BernoulliFactor
     real, parameter :: gammaIH=1.5
     real, parameter :: R1=2.50,R2=12.50
     integer,parameter::nPowerIndex=2
     !------------------------------------------------------------------
     !--
     !\
     ! Calculate cell-centered spherical coordinates::
     RR   = sqrt(xx**2+yy**2+zz**2)
     !\
     ! Avoid calculating inside a critical radius = 0.5*Rsun
     !/
     if (RR <max(Ro_PFSSM-dR*nRExt,0.90*Ro_PFSSM)) then 
        gammaOut= gammaSS
        RETURN
     end if

     ! Calculate gamma
     if(RR >= R2)then
        gammaOut=gammaIH
     else if(RR >= R1)then

        gammaOut=gammaSS+(RR-R1)*(gammaIH-gammaSS)/(R2-R1)
     else
        call get_bernoulli_integral(xx,yy,zz,Uf)
        BernoulliFactor=(cHalf*Uf**2+cSunGravitySI)/&
             (T0*cBoltzmann/cProtonMass/min(Uf/UMin,cTwo))&
             *(R1-RR)*&
             & (Ro_PFSSM/RR)**nPowerIndex/ (R1-Ro_PFSSM)+ GammaSS&
             &/(GammaSS-cOne)*(cOne- (R1-RR)*(Ro_PFSSM/RR)&
             &**nPowerIndex/ (R1-Ro_PFSSM))
        gammaOut = BernoulliFactor/(BernoulliFactor-cOne)
     end if

   end subroutine get_gamma_emp
   !====================================================================
   ! Subroutine write_expansion_tec generates a 2D tecplot output file,
   ! which displays the ModExpansionFactors parameters. All variables 
   ! are displayed in the MAGNETOGRAM frame of reference - 
   ! Mag_Long (Magnetogram longitude) and Mag_Lat (magnetogram latitude).
   !====================================================================
   subroutine write_expansion_tec
     use ModExpansionFactors
     implicit none

     real :: GammaR0, GammaRS
     integer :: iError,iPhi,iTheta,iUnit
     real :: xx,yy,zz,rLatitude

     iUnit=io_unit_new()
     call write_prefix;write(iUnitOut,*)'Writing PFSSM factors  output&
          & file'
     open(unit = iUnit, file = 'SC/IO2/PFSSM_Factors.dat', form =&
          & 'formatted', access = 'sequential', status = 'replace',&
          & iostat = iError )
     ! Tecplot file header  
     if ( iError /= 0 ) then
        call write_prefix;write(iUnitOut, '(a)' ) ' '
        call write_prefix;write(iUnitOut, '(a)' ) 'TECPLOT_WRITE_OPEN -&
             & Fatal error!'
        call write_prefix;write(iUnitOut, '(a)' ) '  Could not open the&
             & output file.'
        call stop_mpi('')
     end if

     write ( iUnit, '(a)' ) 'Title = "'     // trim ('PFSSM_Br') // '"'
     write ( iUnit, '(a)' ) 'Variables = ' // trim ( '"Mag_Long [Deg]",&
          & "Mag_Lat [Deg]",  "f_s","Fisk_factor","Theta_b","U_WSA [Km&
          &/s]","U_Fisk [Km/s]","Gamma0","GammaSS"')
     write ( iUnit, '(a)' ) ' '
     write ( iUnit, '(a,i6,a,i6,a)' ) 'Zone I = ', nPhi+1, ', J=',&
          & nTheta+1, ', F=point' 
     !\
     ! Writing parameters maps:
     ! List of parameters:
     ! 1) 1/fs(Rs), 2) Fisk_factor(Rs), 3) Theta_b(Rs), 4) Final WSA speed
     ! (at Rs),
     ! 5) Final Fisk speed (at Rs), 6) Gamma (R0), 7) Gamma (Rs)  
     !/
     do iTheta=0,nTheta
        rLatitude=r_latitude(iTheta)
        do iPhi=0,nPhi
           xx=Ro_PFSSM*cos(rLatitude)* cos(real(iPhi)*dPhi+Phi_Shift&
                &*cDegToRad)
           yy=Ro_PFSSM*cos(rLatitude)* sin(real(iPhi)*dPhi+Phi_Shift&
                &*cDegToRad)
           zz=Ro_PFSSM*sin(rLatitude) 
           call get_gamma_emp(xx,yy,zz,GammaR0)
           xx=Rs_PFSSM*cos(rLatitude)* cos(real(iPhi)*dPhi+Phi_Shift&
                &*cDegToRad)
           yy=Rs_PFSSM*cos(rLatitude)* sin(real(iPhi)*dPhi+Phi_Shift&
                &*cDegToRad)
           zz=Rs_PFSSM*sin(rLatitude)
           call get_gamma_emp(xx,yy,zz,GammaRS)
           write ( iUnit, '(6f10.3)' )real(iPhi)*dPhi/cDegToRad,&
                rLatitude/cDegToRad,&
                ExpansionFactorInv_N(nR,iPhi,iTheta),&
                FiskFactor_N(nR,iPhi,iTheta),&
                ThetaB_N(nR,iPhi,iTheta),&
                WSAspeed_N(nR,iPhi,iTheta)/cE3,&   !in Km/s
                Fiskspeed_N(nR,iPhi,iTheta)/cE3,&   !in Km/s
                GammaR0,&                       !At the solar surface
                GammaRS                        !At the source surface

        end do
     end do
     close(iUnit)

   end subroutine write_expansion_tec
   !==========================================================================
