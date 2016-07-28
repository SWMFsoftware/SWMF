Module ModSeGrid
  implicit none
  
  private !except

  !field line grid
  real, public, allocatable    :: FieldLineGrid_IC(:,:) 

  ! Grid spacings in each zone in the ionosphere
  real, public                 :: DrIono1,DrIono2,DrIono3,DrIono4

  ! Index Size of Each ionospheric zone
  integer, public, allocatable :: nIono1,nIono2,nIono3,nIono4

  ! Index Size of total ionosphere
!  integer, public, allocatable :: nIono
  
  ! Altitude of Ionospheric and Plasmaspheric base 
  real, public, allocatable :: BaseAltIono, BaseAltPlas

  
  ! Angle Grid
!  real, public, allocatable :: AngleGrid_ICI(:,:,:) 
  real, public, allocatable :: EqAngleGrid_IG(:,:) 
 
  integer, public,allocatable ::nTheta_II(:,:)! PA steps added for each zone on 
                                              !  each line
  real, public, allocatable   :: ThetaZone_II(:,:)! PA width for each zone 
                                                  !  on each line
  real, public,allocatable    :: dTheta_II(:,:)   ! PA spacing in each zone 
                                                  !  on each line
  integer, public,allocatable :: nThetaAlt_II(:,:)! number of PA points at 
                                                  !  each alt
  
  !B Grid
  real, public, allocatable :: Bfield_IC(:,:) ! Bfield in G at each alt step and
                                             ! for each line
  real, public, allocatable :: Lshell_I(:) !Lshell foreach line

  real, public :: nAngle=135 ! number of points in equatorial angle
  real, public :: nPoint=200 ! total number of points on grid
  real, public :: nIono =30  ! number of points in each ionosphere
  real, public :: nPlas =140  ! number of points in plasmasphere
  real, public :: nPlasHalf   ! number of points from BasePlas to equator
  real, public :: nTop   ! number of points to top (in open) 
                         !   or equatorial (in closed)

  real, public :: nLine=1! number of field lines
  real, public :: nZone=4! number of pitchangle zones
  real, public :: MaxTheta

  real, allocatable :: BFieldEq_I(:), BFieldIono_I(:)


  !E Grid
  character(len=10) :: TypeGridE
  integer,public :: nEnergy=80
  real, allocatable, public :: DeltaE_I(:),EnergyGrid_I(:)
  real,   public :: EnergyMin, EnergyMax, DeltaE
  
  public :: se_grid_test



  real :: rPlanetCM
contains
  !============================================================================
  subroutine init_se_grid
    use ModPlanetConst, ONLY: Earth_, rPlanet_I
    use ModNumConst,    ONLY: cPi
    real, parameter :: cMtoCM = 1.0e2
    integer :: iLine
    real    :: Biono, Beq, MLAT1, PhiBasePlas, QO, SphiO
    !--------------------------------------------------------------------------

    nIono=nIono1+nIono2+nIono3+nIono4
    BaseAltPlas = &
         DrIono1*nIono1+DrIono2*nIono2+DrIono3*nIono3+DrIono4*nIono4
    BaseAltIono = 9e6 !cm
    
    nPlas=nPoint-2*nIono
    nTop = nPoint/2
    nPlasHalf=nPlas/2
    
    do iLine = 1, nLine
       rPlanetCM=rPlanet_I(Earth_)*cMtoCM
       ! Set field line info
       !       write(*,*) 'rPlanetCM,BaseAltPlas,Lshell_I(iLine)',rPlanetCM,BaseAltPlas,Lshell_I(iLine)
       PhiBasePlas=&
            ACOS(SQRT((rPlanetCM+BaseAltPlas)/(Lshell_I(iLine)*rPlanetCM)))
       write(*,*) 'PhiBasePlas*180./cPi',PhiBasePlas*180./cPi
       SphiO=SIN(PhiBasePlas)
       QO=SQRT(1+3*SphiO**2)
       MLAT1=PhiBasePlas*180./cPi
       Beq   = 0.31/Lshell_I(iLine)**3
       Biono = 0.31*QO/(Lshell_I(iLine)**3*(1-SphiO**2)**3)
       write(*,*) 'calc bfield and s grid for iLine = ', iLine
       call calc_bfield_sgrid(iLine,Biono,PhiBasePlas)
       write(*,*) 'calc equatorial PA grid', iLine
       call calc_equatorial_pitchangle(iLine,Beq,Biono)
    enddo
  end subroutine init_se_grid

  !============================================================================
  subroutine calc_bfield_sgrid(iLine, Biono,PhiBasePlas)
    use ModNumConst,    ONLY: cDegToRad,cRadToDeg
    integer, intent(in) :: iLine
    real   , intent(in) :: Biono,PhiBasePlas
    integer :: iAlt
    real :: Phi, LastPhi
    real :: DeltaS, dPhi
    !length of iono regions
    real :: LengthIono1,LengthIono2,LengthIono3,LengthIono4 

!    if(.not.allocated(Bfield_IC)) allocate(Bfield_IC(nLine,nPoint))
    
    ! set dPhi
    dPhi=PhiBasePlas/nPlasHalf
!    write(*,*) 'nPlasHalf,dPhi*cRadToDeg',nPlasHalf,dPhi*cRadToDeg

    !set Length of iono regions
    LengthIono1 = DrIono1*(nIono1-1)
    LengthIono2 = DrIono2*(nIono2)
    LengthIono3 = DrIono3*(nIono3)
    LengthIono4 = DrIono4*(nIono4)

    
    do iAlt = 1, nPoint
       write(*,*) iAlt
       if (iAlt <= nIono1) then
          !set alt zone 1 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = BaseAltIono+(iAlt-1)*DrIono1
       elseif(iAlt <= nIono1+nIono2) then
          !set alt zone 2 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               BaseAltIono+LengthIono1+(iAlt-nIono1)*DrIono2
       elseif(iAlt <=nIono1+nIono2+nIono3) then
          !set alt zone 3 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               BaseAltIono+LengthIono1+LengthIono2 &
               +(iAlt-nIono1-nIono2)*DrIono3
       elseif(iAlt <=nIono1+nIono2+nIono3+nIono4) then
          !set alt zone 4 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
                  BaseAltIono+LengthIono1+LengthIono2+LengthIono3&
                  +(iAlt-nIono1-nIono2-nIono3)*DrIono4
          
       elseif(iAlt < nPoint-nIono) then
          ! set plasmaspheric region
          ! Get polar angle at fieldline point
          Phi = abs (PhiBasePlas-(iAlt-nIono)*dPhi)
          LastPhi = abs (PhiBasePlas-(iAlt-nIono-1)*dPhi)
          call get_b_deltaS_point(Phi,LastPhi,DeltaS,Lshell_I(iLine),Bfield_IC(iLine,iAlt))
          ! Get B and Delta s at each point
!          write(*,*) 'Bfield_IC(iLine,iAlt)',Bfield_IC(iLine,iAlt)
          FieldLineGrid_IC(iLine,iAlt)=FieldLineGrid_IC(iLine,iAlt-1)+DeltaS
          
       elseif(iAlt < nPoint-nIono1-nIono2-nIono3) then
          !set alt zone 4 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono-1)&
               +(iAlt-nPlas-nIono)*DrIono4
       elseif(iAlt < nPoint-nIono1-nIono2) then
          !set alt zone 3 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono1-nIono2-nIono3-1)&
               +(iAlt-nIono4-nPlas-nIono)*DrIono3
       elseif(iAlt < nPoint-nIono1)then
          !set alt zone 2 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono1-nIono2-1)&
               +(iAlt-nIono4-nIono3-nPlas-nIono)*DrIono3
       elseif(iAlt <= nPoint)then
          !set alt zone 1 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono1-nIono2-1)&
               +(iAlt-nIono4-nIono3-nIono2-nPlas-nIono)*DrIono3
       endif
       
!       write(*,*) 'iAlt,FieldLineGrid_IC(iLine,iAlt)',iAlt,FieldLineGrid_IC(iLine,iAlt)/1e5
    end do
    
    ! Set magnetic field in ionosphere to value at top of ionosphere
    Bfield_IC(iLine,1:nIono)= Biono
    Bfield_IC(iLine,nPoint-nIono:nPoint)= Biono
!    write(*,*) 'Biono',Biono
  end subroutine calc_bfield_sgrid
  !============================================================================
  
  subroutine calc_equatorial_pitchangle(iLine,Beq,Biono)
    use ModNumConst,    ONLY: cPi
    integer, intent(in) :: iLine
    real,    intent(in) :: Beq, Biono
    integer :: iAlt, iAngle
    real    :: LocalThetaMax

    real, allocatable :: tmp_array(:)
    !--------------------------------------------------------------------------
    
    
    ! Set Maximum number of points in equatorial PA
    MaxTheta = sum(nTheta_II(iLine,:))
    
    if (.not.allocated(tmp_array)) allocate(tmp_array(MaxTheta))
    
    ! Set the theta range for each zone and each line
    !  (currently assumes 4 zones)
    ThetaZone_II(iLine,1) = asin(sqrt(Beq/Biono)) !zone 1 is the loss cone
    ThetaZone_II(iLine,2) = ThetaZone_II(iLine,1)
    ThetaZone_II(iLine,3) = 1.2-2.0*ThetaZone_II(iLine,1)
    ThetaZone_II(iLine,4) = 0.5*cPi-1.2
    
    ! Get DeltaTheta for each zone
    dTheta_II(iLine,:) = ThetaZone_II(iLine,:)/nTheta_II(iLine,:)
    
    ! Fill in compuational region of equatorial PA grid
    do iAngle=1,MaxTheta
       if (iAngle <= nTheta_II(iLine,1)) then 
          EqAngleGrid_IG(iLine,iAngle) = dTheta_II(iLine,1) * iAngle
       elseif (iAngle <= sum(nTheta_II(iLine,1:2))) then
          EqAngleGrid_IG(iLine,iAngle) = &
               dTheta_II(iLine,2) * (iAngle-nTheta_II(iLine,1)) &
               + ThetaZone_II(iLine,1)
       elseif (iAngle <= sum(nTheta_II(iLine,1:3))) then
          EqAngleGrid_IG(iLine,iAngle) = &
               dTheta_II(iLine,3) * (iAngle-sum(nTheta_II(iLine,1:2))) &
               + sum(ThetaZone_II(iLine,1:2))
       else 
          EqAngleGrid_IG(iLine,iAngle) = &
               dTheta_II(iLine,4) * (iAngle-sum(nTheta_II(iLine,1:3))) &
               + sum(ThetaZone_II(iLine,1:3))
       end if
 
       write(*,*) 'iAngle,EqAngleGrid_IG(iLine,iAngle)',iAngle,EqAngleGrid_IG(iLine,iAngle)
   end do
    
    ! Fill in Ghost Cells of Equatorial Angle Grid
    EqAngleGrid_IG(iLine,0) = 0.0
    EqAngleGrid_IG(iLine,MaxTheta) = 0.5 * cPi
    
    ! Fill in nThetaAlt array (Number of points in theta for each alt)
    do iAlt = 1, nPoint
       ! find equatorial PA corresponding to local Pi/2  
       LocalThetaMax = asin(sqrt(min(1.0,Beq/Bfield_IC(iLine,iAlt))))
       
       ! From LocalThetaMax and EqAngleGrid find nThetaAlt
       Tmp_array = 0
       where(EqAngleGrid_IG(iLine,:) <= LocalThetaMax)
          Tmp_array = 1 
       end where
       nThetaAlt_II(iLine,iAlt) = sum(Tmp_array)
    enddo
    
  end subroutine calc_equatorial_pitchangle

  !============================================================================
  ! Subroutine SpaceVar calculates h and Bfield for a given spatial step.
  !  VARIABLE DESCRIPTIONS
  !      Q       Used in distance calculations, SQRT(1+3*(SIN(phi))**2)
  !      lQ      Same as Q, but for lphi
  !      phi     Azimuthal angle for the base of the plasmasphere; rad
  !      L       Equatorial distance in Earth radii of flux tube
  !      Re      Radius of the Earth; 8.378E8 cm
  !      lphi    Azimuthal angle for the previous step; radians
  !      h`       Size of this step; cm
  !      B       Magnetic field strength for this step; G
  !
  SUBROUTINE get_b_deltaS_point(phi,lphi,h,L,B)
    use ModPlanetConst, ONLY: Earth_
    REAL phi,lphi,h,L,B,Q,lQ
    !---------------------------------------------------------------------------
    
    Q=SQRT(1+3*(SIN(phi))**2)
    lQ=SQRT(1+3*(SIN(lphi))**2)
    h=ABS(.5*L*rPlanetCM*(lQ*SIN(lphi)-Q*SIN(phi)+1/SQRT(3.)*LOG((SQRT(3.) &
         *SIN(lphi)+lQ)/(SQRT(3.)*SIN(phi)+Q))))
    B=0.31*Q/(L**3*(COS(phi))**6)
    RETURN
  END SUBROUTINE get_b_deltaS_point
  !============================================================================
  subroutine calc_energy_grid
    implicit none
    
    private !except
    
  !field line grid
  real, public, allocatable    :: FieldLineGrid_IC(:,:) 

  ! Grid spacings in each zone in the ionosphere
  real, public                 :: DrIono1,DrIono2,DrIono3,DrIono4

  ! Index Size of Each ionospheric zone
  integer, public, allocatable :: nIono1,nIono2,nIono3,nIono4

  ! Index Size of total ionosphere
!  integer, public, allocatable :: nIono
  
  ! Altitude of Ionospheric and Plasmaspheric base 
  real, public, allocatable :: BaseAltIono, BaseAltPlas

  
  ! Angle Grid
!  real, public, allocatable :: AngleGrid_ICI(:,:,:) 
  real, public, allocatable :: EqAngleGrid_IG(:,:) 
 
  integer, public,allocatable ::nTheta_II(:,:)! PA steps added for each zone on 
                                              !  each line
  real, public, allocatable   :: ThetaZone_II(:,:)! PA width for each zone 
                                                  !  on each line
  real, public,allocatable    :: dTheta_II(:,:)   ! PA spacing in each zone 
                                                  !  on each line
  integer, public,allocatable :: nThetaAlt_II(:,:)! number of PA points at 
                                                  !  each alt
  
  !B Grid
  real, public, allocatable :: Bfield_IC(:,:) ! Bfield in G at each alt step and
                                             ! for each line
  real, public, allocatable :: Lshell_I(:) !Lshell foreach line

  real, public :: nAngle=135 ! number of points in equatorial angle
  real, public :: nPoint=200 ! total number of points on grid
  real, public :: nIono =30  ! number of points in each ionosphere
  real, public :: nPlas =140  ! number of points in plasmasphere
  real, public :: nPlasHalf   ! number of points from BasePlas to equator
  real, public :: nTop   ! number of points to top (in open) 
                         !   or equatorial (in closed)

  real, public :: nLine=1! number of field lines
  real, public :: nZone=4! number of pitchangle zones
  real, public :: MaxTheta

  real, allocatable :: BFieldEq_I(:), BFieldIono_I(:)


  !E Grid
  character(len=10) :: TypeGridE
  integer,public :: nEnergy=80
  real, allocatable, public :: DeltaE_I(:)
  real,   public :: EnergyMin, EnergyMax
  
  public :: se_grid_test



  real :: rPlanetCM
contains
  !============================================================================
  subroutine init_se_grid
    use ModPlanetConst, ONLY: Earth_, rPlanet_I
    use ModNumConst,    ONLY: cPi
    real, parameter :: cMtoCM = 1.0e2
    integer :: iLine
    real    :: Biono, Beq, MLAT1, PhiBasePlas, QO, SphiO
    !--------------------------------------------------------------------------

    nIono=nIono1+nIono2+nIono3+nIono4
    BaseAltPlas = &
         DrIono1*nIono1+DrIono2*nIono2+DrIono3*nIono3+DrIono4*nIono4
    BaseAltIono = 9e6 !cm
    
    nPlas=nPoint-2*nIono
    nTop = nPoint/2
    nPlasHalf=nPlas/2
    
    do iLine = 1, nLine
       rPlanetCM=rPlanet_I(Earth_)*cMtoCM
       ! Set field line info
!       write(*,*) 'rPlanetCM,BaseAltPlas,Lshell_I(iLine)',rPlanetCM,BaseAltPlas,Lshell_I(iLine)
       PhiBasePlas=&
            ACOS(SQRT((rPlanetCM+BaseAltPlas)/(Lshell_I(iLine)*rPlanetCM)))
       write(*,*) 'PhiBasePlas*180./cPi',PhiBasePlas*180./cPi
       SphiO=SIN(PhiBasePlas)
       QO=SQRT(1+3*SphiO**2)
       MLAT1=PhiBasePlas*180./cPi
       Beq   = 0.31/Lshell_I(iLine)**3
       Biono = 0.31*QO/(Lshell_I(iLine)**3*(1-SphiO**2)**3)
       write(*,*) 'calc bfield and s grid for iLine = ', iLine
       call calc_bfield_sgrid(iLine,Biono,PhiBasePlas)
       write(*,*) 'calc equatorial PA grid', iLine
       call calc_equatorial_pitchangle(iLine,Beq,Biono)
    enddo
  end subroutine init_se_grid

  !============================================================================
  subroutine calc_bfield_sgrid(iLine, Biono,PhiBasePlas)
    use ModNumConst,    ONLY: cDegToRad,cRadToDeg
    integer, intent(in) :: iLine
    real   , intent(in) :: Biono,PhiBasePlas
    integer :: iAlt
    real :: Phi, LastPhi
    real :: DeltaS, dPhi
    !length of iono regions
    real :: LengthIono1,LengthIono2,LengthIono3,LengthIono4 

!    if(.not.allocated(Bfield_IC)) allocate(Bfield_IC(nLine,nPoint))
    
    ! set dPhi
    dPhi=PhiBasePlas/nPlasHalf
!    write(*,*) 'nPlasHalf,dPhi*cRadToDeg',nPlasHalf,dPhi*cRadToDeg

    !set Length of iono regions
    LengthIono1 = DrIono1*(nIono1-1)
    LengthIono2 = DrIono2*(nIono2)
    LengthIono3 = DrIono3*(nIono3)
    LengthIono4 = DrIono4*(nIono4)

    
    do iAlt = 1, nPoint
       write(*,*) iAlt
       if (iAlt <= nIono1) then
          !set alt zone 1 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = BaseAltIono+(iAlt-1)*DrIono1
       elseif(iAlt <= nIono1+nIono2) then
          !set alt zone 2 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               BaseAltIono+LengthIono1+(iAlt-nIono1)*DrIono2
       elseif(iAlt <=nIono1+nIono2+nIono3) then
          !set alt zone 3 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               BaseAltIono+LengthIono1+LengthIono2 &
               +(iAlt-nIono1-nIono2)*DrIono3
       elseif(iAlt <=nIono1+nIono2+nIono3+nIono4) then
          !set alt zone 4 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
                  BaseAltIono+LengthIono1+LengthIono2+LengthIono3&
                  +(iAlt-nIono1-nIono2-nIono3)*DrIono4
          
       elseif(iAlt < nPoint-nIono) then
          ! set plasmaspheric region
          ! Get polar angle at fieldline point
          Phi = abs (PhiBasePlas-(iAlt-nIono)*dPhi)
          LastPhi = abs (PhiBasePlas-(iAlt-nIono-1)*dPhi)
          call get_b_deltaS_point(Phi,LastPhi,DeltaS,Lshell_I(iLine),Bfield_IC(iLine,iAlt))
          ! Get B and Delta s at each point
!          write(*,*) 'Bfield_IC(iLine,iAlt)',Bfield_IC(iLine,iAlt)
          FieldLineGrid_IC(iLine,iAlt)=FieldLineGrid_IC(iLine,iAlt-1)+DeltaS
          
       elseif(iAlt < nPoint-nIono1-nIono2-nIono3) then
          !set alt zone 4 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono-1)&
               +(iAlt-nPlas-nIono)*DrIono4
       elseif(iAlt < nPoint-nIono1-nIono2) then
          !set alt zone 3 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono1-nIono2-nIono3-1)&
               +(iAlt-nIono4-nPlas-nIono)*DrIono3
       elseif(iAlt < nPoint-nIono1)then
          !set alt zone 2 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono1-nIono2-1)&
               +(iAlt-nIono4-nIono3-nPlas-nIono)*DrIono3
       elseif(iAlt <= nPoint)then
          !set alt zone 1 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono1-nIono2-1)&
               +(iAlt-nIono4-nIono3-nIono2-nPlas-nIono)*DrIono3
       endif
       
!       write(*,*) 'iAlt,FieldLineGrid_IC(iLine,iAlt)',iAlt,FieldLineGrid_IC(iLine,iAlt)/1e5
    end do
    
    ! Set magnetic field in ionosphere to value at top of ionosphere
    Bfield_IC(iLine,1:nIono)= Biono
    Bfield_IC(iLine,nPoint-nIono:nPoint)= Biono
!    write(*,*) 'Biono',Biono
  end subroutine calc_bfield_sgrid
  !============================================================================
  
  subroutine calc_equatorial_pitchangle(iLine,Beq,Biono)
    use ModNumConst,    ONLY: cPi
    integer, intent(in) :: iLine
    real,    intent(in) :: Beq, Biono
    integer :: iAlt, iAngle
    real    :: LocalThetaMax

    real, allocatable :: tmp_array(:)
    !--------------------------------------------------------------------------
    
    
    ! Set Maximum number of points in equatorial PA
    MaxTheta = sum(nTheta_II(iLine,:))
    
    if (.not.allocated(tmp_array)) allocate(tmp_array(MaxTheta))
    
    ! Set the theta range for each zone and each line
    !  (currently assumes 4 zones)
    ThetaZone_II(iLine,1) = asin(sqrt(Beq/Biono)) !zone 1 is the loss cone
    ThetaZone_II(iLine,2) = ThetaZone_II(iLine,1)
    ThetaZone_II(iLine,3) = 1.2-2.0*ThetaZone_II(iLine,1)
    ThetaZone_II(iLine,4) = 0.5*cPi-1.2
    
    ! Get DeltaTheta for each zone
    dTheta_II(iLine,:) = ThetaZone_II(iLine,:)/nTheta_II(iLine,:)
    
    ! Fill in compuational region of equatorial PA grid
    do iAngle=1,MaxTheta
       if (iAngle <= nTheta_II(iLine,1)) then 
          EqAngleGrid_IG(iLine,iAngle) = dTheta_II(iLine,1) * iAngle
       elseif (iAngle <= sum(nTheta_II(iLine,1:2))) then
          EqAngleGrid_IG(iLine,iAngle) = &
               dTheta_II(iLine,2) * (iAngle-nTheta_II(iLine,1)) &
               + ThetaZone_II(iLine,1)
       elseif (iAngle <= sum(nTheta_II(iLine,1:3))) then
          EqAngleGrid_IG(iLine,iAngle) = &
               dTheta_II(iLine,3) * (iAngle-sum(nTheta_II(iLine,1:2))) &
               + sum(ThetaZone_II(iLine,1:2))
       else 
          EqAngleGrid_IG(iLine,iAngle) = &
               dTheta_II(iLine,4) * (iAngle-sum(nTheta_II(iLine,1:3))) &
               + sum(ThetaZone_II(iLine,1:3))
       end if
 
       write(*,*) 'iAngle,EqAngleGrid_IG(iLine,iAngle)',iAngle,EqAngleGrid_IG(iLine,iAngle)
   end do
    
    ! Fill in Ghost Cells of Equatorial Angle Grid
    EqAngleGrid_IG(iLine,0) = 0.0
    EqAngleGrid_IG(iLine,MaxTheta) = 0.5 * cPi
    
    ! Fill in nThetaAlt array (Number of points in theta for each alt)
    do iAlt = 1, nPoint
       ! find equatorial PA corresponding to local Pi/2  
       LocalThetaMax = asin(sqrt(min(1.0,Beq/Bfield_IC(iLine,iAlt))))
       
       ! From LocalThetaMax and EqAngleGrid find nThetaAlt
       Tmp_array = 0
       where(EqAngleGrid_IG(iLine,:) <= LocalThetaMax)
          Tmp_array = 1 
       end where
       nThetaAlt_II(iLine,iAlt) = sum(Tmp_array)
    enddo
    
  end subroutine calc_equatorial_pitchangle

  !============================================================================
  ! Subroutine SpaceVar calculates h and Bfield for a given spatial step.
  !  VARIABLE DESCRIPTIONS
  !      Q       Used in distance calculations, SQRT(1+3*(SIN(phi))**2)
  !      lQ      Same as Q, but for lphi
  !      phi     Azimuthal angle for the base of the plasmasphere; rad
  !      L       Equatorial distance in Earth radii of flux tube
  !      Re      Radius of the Earth; 8.378E8 cm
  !      lphi    Azimuthal angle for the previous step; radians
  !      h`       Size of this step; cm
  !      B       Magnetic field strength for this step; G
  !
  SUBROUTINE get_b_deltaS_point(phi,lphi,h,L,B)
    use ModPlanetConst, ONLY: Earth_
    REAL phi,lphi,h,L,B,Q,lQ
    !---------------------------------------------------------------------------
    
    Q=SQRT(1+3*(SIN(phi))**2)
    lQ=SQRT(1+3*(SIN(lphi))**2)
    h=ABS(.5*L*rPlanetCM*(lQ*SIN(lphi)-Q*SIN(phi)+1/SQRT(3.)*LOG((SQRT(3.) &
         *SIN(lphi)+lQ)/(SQRT(3.)*SIN(phi)+Q))))
    B=0.31*Q/(L**3*(COS(phi))**6)
    RETURN
  END SUBROUTINE get_b_deltaS_point
  !============================================================================
  subroutine calc_energy_grid

    select case(TypeGridE)
    case('ConstDE')
       EnergyMin=EnergyMax-nEnergy*DeltaE_I(1)
       do iEnergy = 1, nEnergy
          !Constant DeltaE grid
          Energy_I(iEnergy)=EnergyMin+DeltaE*(iEnergy-0.5)
          DeltaE_I(iEnergy)=DeltaE
       end do
    case('ConstLogDE')
       ! Constant log(DeltaE) grid
       EnergySide2=alog(EnergyMax)+0.5*DeltaE-DeltaE*nEnergy
       EnergyMin=exp(EnergySide2)
       do iEnergy =1, nEnergy
          EnergySide1=EnergySide2
          EnergySide2=EnergySide1+DeltaE
          DeltaE_I(iEnergy) = exp(EnergySide2)-exp(EnergySide1)
          Energy_I(iEnergy) = exp(EnergySide1)+0.5*DeltaE_I(iEnergy)
       end do
       Energymax=exp(EnergySide2)
    case('Staggered')
       !Staggered dE grid. Two regions with 2eV and 4eV spacing
       EnergyMin=1.0			! Input Emax is Emax
       Energy_I(1)=1.5
       DeltaE_I(1)=1.0
       do iEnergy=2,20	
          Energy_I(iEnergy)=REAL(2*iEnergy)-1.	
	  DeltaE_I(iEnergy)=2.0
       end do
       
       do iEnergy=21,21+nEnergy
          Energy_I(iEnergy)=40.+(0.5+iEnergy-21)*DeltaE
	  DeltaE_I(iEnergy)=DeltaE
       end do
       EnergyMax=Energy_I(nEnergy)+0.5*DeltaE_I(nEnergy)
       
    case('Geometric')
       EnergyMin=2.			! Given unless input file is changed
       Delta_I(1)=DeltaE
       Energy_I(1)=EnergyMin+0.5*(1)
       do iEnergy=2,nEnergy
          ! Input Emax is dE growth factor
	  DeltaE(iEnergy)=EnergyMax*DeltaE(iEnergy-1)
   Energy_I(iEnergy)=&
        Energy_I(iEnergy-1)+.5*(DeltaE_I(iEnergy-1)+DeltaE_I(iEnergy))
	end do
	Emax=ener(Jo)+.5*del(Jo)
      END IF

  end subroutine calc_energy_grid
  !==========================================================================
  subroutine allocate_grid_arrays
    
    ! Allocate PA related grid
    if(.not.allocated(nTheta_II))       allocate(nTheta_II(nLine,nZone))
    if(.not.allocated(ThetaZone_II))    allocate(ThetaZone_II(nLine,nZone))
    if(.not.allocated(dTheta_II))       allocate(dTheta_II(nLine,nZone))
    if(.not.allocated(nThetaAlt_II))    allocate(nThetaAlt_II(nLine,nPoint))
    if(.not.allocated(EqAngleGrid_IG))  allocate(EqAngleGrid_IG(nLine,nPoint))
    if(.not.allocated(Bfield_IC))       allocate(Bfield_IC(nLine,nPoint))
    if(.not.allocated(Lshell_I))        allocate(Lshell_I(nLine))
    if(.not.allocated(FieldLineGrid_IC))allocate(FieldLineGrid_IC(nLine,nPoint))
    if(.not.allocated(DeltaE_I))        allocate(DeltaE_I(nEnergy))
    if(.not.allocated(EnergyGrid_I))    allocate(EnergyGrid_I(nEnergy))
  end subroutine allocate_grid_arrays

  !============================================================================
  ! save grid plot for verification
  subroutine plot_grid
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    integer, parameter :: nDim =2, nVar=2, S_=1, PA_=2,B_=1
    character(len=100),parameter :: NamePlotVar='S PA B PA g r'
    character(len=12) :: NamePlot='GridPlot.out'
    character(len=*),parameter :: NameHeader='Grid output'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iLine,iPoint
    !--------------------------------------------------------------------------
    allocate(Coord_DII(nDim,nPoint,nAngle), PlotState_IIV(nPoint,nAngle,nVar))

    do iLine=1,nLine
       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
          do iPoint=1,nPoint
             do iAngle=1,nAngle
             Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)/6375.0e5
             Coord_DII(PA_,iPoint,iAngle)= EqAngleGrid_IG(iLine,iAngle)
             if (iAngle <= nThetaAlt_II(iLine,iPoint) )then
                PlotState_IIV(iPoint,iAngle,B_) = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = EqAngleGrid_IG(iLine,iAngle)
             else
                PlotState_IIV(iPoint,iAngle,:)=0.0
             endif
          enddo
       enddo
            
       !Plot grid for given line
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= 1,TimeIn=1.0,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
    end do
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_grid
  !============================================================================
  ! UNIT test for SE grid
  subroutine se_grid_test
    
    DrIono1 = 1e6
    nIono1  = 12
    DrIono2 = 2e6
    nIono2  = 10
    DrIono3 = 3e6
    nIono3  = 10
    DrIono4 = 5e6
    nIono4  = 2    

    ! Allocated the grid arrays and populate the bfield, sgrid, and PA grid  
    write(*,*) 'allocating arrays'
    call allocate_grid_arrays

    Lshell_I(1)=2.5
    nTheta_II(1,1)=5
    nTheta_II(1,2)=20
    nTheta_II(1,3)=90
    nTheta_II(1,4)=20

    nAngle = sum(nTheta_II(1,:))

    write(*,*) nLine,Lshell_I
    write(*,*) 'calling init_se_grid'
    call init_se_grid
    write(*,*) 'calling plot_grid'
    call plot_grid
  end subroutine se_grid_test
end Module ModSeGrid








  !============================================================================
  ! UNIT test for SE grid
  subroutine se_grid_test
    
    DrIono1 = 1e6
    nIono1  = 12
    DrIono2 = 2e6
    nIono2  = 10
    DrIono3 = 3e6
    nIono3  = 10
    DrIono4 = 5e6
    nIono4  = 2    

    ! Allocated the grid arrays and populate the bfield, sgrid, and PA grid  
    write(*,*) 'allocating arrays'
    call allocate_grid_arrays

    Lshell_I(1)=2.5
    nTheta_II(1,1)=5
    nTheta_II(1,2)=20
    nTheta_II(1,3)=90
    nTheta_II(1,4)=20

    nAngle = sum(nTheta_II(1,:))

    write(*,*) nLine,Lshell_I
    write(*,*) 'calling init_se_grid'
    call init_se_grid
    write(*,*) 'calling plot_grid'
    call plot_grid
  end subroutine se_grid_test
end Module ModSeGrid
