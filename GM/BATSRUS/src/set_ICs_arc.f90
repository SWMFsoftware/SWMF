!^CFG COPYRIGHT UM
!^CFG FILE NOT SIMPLE
!---------------------------------------------------------------------------
!Subroutine                   arcade_init
!---------------------------------------------------------------------------
!\
! Calculates magnetic field for sheared magnetic arcade that is 
! unstable to erupting resulting in a CME
! sheared intrusions of the analytical form phi = exp(-z)*(1.0+cos(x))
! see Low and Manchester, AstroPhysical Journal, Vol 528, p. 1026 for 
! forther details
!
! Written by Chip Manchester Oct 11, 2000
!/
!
!===========================================================================
subroutine set_ICs_arc
  use ModMain, ONLY : nI,nJ,nK,gcn,globalBLK
  use ModVarIndexes,ONLY:rho_,rhoUx_,rhoUy_,&
       rhoUz_,Bx_,By_,Bz_,P_
  use ModAdvance, ONLY : B0xCell_BLK,B0yCell_BLK,B0zCell_BLK, &
       B0xFace_x_BLK,B0yFace_x_BLK,B0zFace_x_BLK, &
       B0xFace_y_BLK,B0yFace_y_BLK,B0zFace_y_BLK, &
       B0xFace_z_BLK,B0yFace_z_BLK,B0zFace_z_BLK, &
       State_VGB
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK
  use ModPhysics
  implicit none

  integer :: i, j, k
  real ::  x, y, z
  real ::  phi, zs
  !
  !  PARAMETER LIST:  expArc,widthArc,phi0Arc,Gsun,muArc,BArcDim,ByArcDim,
  !                   TArc,RhoArc,UzPertArc
  !  Definition of Parameters used for the initial state
  !  expArc    = parameter of analytical expression that defines flux funtion
  !  widthArc   = sets width of arcade at -infinity = 2pi*widthArc*P_height
  !  phi0Arc     = flux cooridnate for the top-most field line of the arcade 
  !  Gsun     = gravitational acceleration at sun's surface = 2.734e4 cm/s**2
  !  muArc       = the average mass of the gas atoms in proton mass = 1.3
  !  BArcDim   = magnetic field strength parameter of the arcades
  !  ByArcDim  = strenght of the pervailing longitudinal B field
  !  TArc     = the isothermal temperature for the model
  !  RhoArc   = density of corona at z = 0
  !  UzPertArc   = magnitude of velocity perturbation
  !  Phtscl   = 1.0  scaled quantities for P_height in BATSRUS
  !  RHOscl   = 1.0                        RhoArc
  !  SSPscl   = 1.0                        SSPsun
  !  Vscl     = UzPertArc/SSPsun           UzPertArc

  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////////
  !=====================================================================
  !  if(iProc==0)then
  !     write(6,*) 'arcade_init called'
  !     write(6,*) 'gamma = ',g
  !     write(6,*) 'Gsun  = ',Gsun
  !     write(6,*) 'TArc  = ',TArc
  !     write(6,*) 'SSPsun = ',SSPsun
  !     write(6,*) 'RhoArc = ',RhoArc
  !     write(6,*) 'Scale height = ',P_height
  !     write(6,*) 'BArcDim = ',BArcDim
  !     write(6,*) 'ByArcDim = ',ByArcDim
  !  endif
  !
  do i=1-gcn,nI+gcn
     do j=1-gcn,nJ+gcn
        do k=1-gcn,nK+gcn
           !CALCULATE CELL CENTER SCALAR AND VECTOR FIELD VALUES FOR THE ARCADE
           x = x_BLK(i,j,k,globalBLK)
           y = y_BLK(i,j,k,globalBLK)
           z = z_BLK(i,j,k,globalBLK)
           phi = exp(-z/Phtscl)*( 1.0 + cos(x/(widthArc*Phtscl)) )
           !PLASMA DENSITY AND PRESSURE  
           !             Magnetized Region
           if( (phi .gt. phi0Arc) .and. (abs(x) .lt.cPi*widthArc*Phtscl) ) then
              State_VGB(rho_,i,j,k,globalBLK) = exp(-z/Phtscl)*( RHOscl  -  &
                   (phi/SSPscl**2)*((B0_scl/widthArc)*(phi -phi0Arc))**2 / &
                   phi**(2*expArc) )
              State_VGB(P_,i,j,k,globalBLK) = &
                   State_VGB(rho_,i,j,k,globalBLK) * &
                   SSPscl**2
              !             Non magnetic region
           else
              State_VGB(rho_,i,j,k,globalBLK) = exp(-z/Phtscl)*RHOscl
              State_VGB(P_,i,j,k,globalBLK) = exp(-z/Phtscl)*RHOscl*SSPscl**2
           endif
           !VELOCITY PERTURBATION
           State_VGB(rhoUx_,i,j,k,globalBLK) = 0.0
           State_VGB(rhoUy_,i,j,k,globalBLK) = 0.0
           zs = Phtscl*log( abs( (1.0 + cos(x/(widthArc*Phtscl))) /phi0Arc ) )
           if  (abs(x) .lt. 0.25*widthArc*cPi*Phtscl)  then 
              State_VGB(rhoUz_,i,j,k,globalBLK) = Vscl  *  &
                   State_VGB(rho_,i,j,k,globalBLK) *  & 
                   abs(cos((x+cTwoPi*Phtscl)/(2.0*Phtscl))) *  &
                   exp(-abs(4.0*(z-zs)/Phtscl))
           else 
              State_VGB(rhoUz_,i,j,k,globalBLK) = 0.0
           endif
           !X COMPONENT OF MAGNETIC FIELD
           !             Magnetized Region
           if( (phi .gt. phi0Arc) .and. (abs(x) .lt.cPi*widthArc*Phtscl)) then
              State_VGB(Bx_,i,j,k,globalBLK) = -B0_scl  *  & 
                   phi*(phi - phi0Arc)/(phi**expArc) 
              !             Non Magnetic Region            
           else
              State_VGB(Bx_,i,j,k,globalBLK) = 0.0
           endif
           !Y COMPONENT OF MAGNETIC FIELD
           !             Magnetized Region
           if( (phi .gt. phi0Arc) .and. (abs(x) .lt.cPi*widthArc*Phtscl) ) then
              State_VGB(By_,i,j,k,globalBLK) = sqrt( B0y_scl**2 +  & 
                   (1.0/widthArc**2 - 1.0)*(B0_scl*phi*(phi-phi0Arc)/phi**expArc)**2 ) 
              !             Non Magnetic Region
           else
              State_VGB(By_,i,j,k,globalBLK) = B0y_scl
           endif
           !Z COMPONENT OF MAGNETIC FIELD 
           !             Magnetized Region
           if( (phi .gt. phi0Arc) .and. (abs(x) .lt.cPi*widthArc*Phtscl) ) then
              State_VGB(Bz_,i,j,k,globalBLK) = (B0_scl/widthArc)  *  &
                   ((phi - phi0Arc)/phi**expArc)*exp(-z/Phtscl)*sin(x/(widthArc*Phtscl))
              !             Non Magnetic Region            
           else
              State_VGB(Bz_,i,j,k,globalBLK) = 0.0
           endif
           !
           !B0 THE "INTRINSIC" MAGNETIC FIELD cell centered
           B0xCell_BLK(i,j,k,globalBLK) = 0.0
           B0yCell_BLK(i,j,k,globalBLK) = 0.0
           B0zCell_BLK(i,j,k,globalBLK) = 0.0
           !PRESSURE TEST
           if(State_VGB(P_,i,j,k,globalBLK) .lt. 0.0) then
              write(*,*) 'ARCADE: negative pressure at', i,j,k,globalBLK
              stop
           end if
        end do
     end do
  end do

  do i= 0,nI+1
     do j=0,nJ+1
        do k=0,nK+1
           !B0 THE "INTRINSIC" MAGNETIC FIELD face centered
           B0xFace_x_BLK(i,j,k,globalBLK) = 0.0
           B0yFace_x_BLK(i,j,k,globalBLK) = 0.0
           B0zFace_x_BLK(i,j,k,globalBLK) = 0.0
           B0xFace_y_BLK(i,j,k,globalBLK) = 0.0
           B0yFace_y_BLK(i,j,k,globalBLK) = 0.0
           B0zFace_y_BLK(i,j,k,globalBLK) = 0.0
           B0xFace_z_BLK(i,j,k,globalBLK) = 0.0
           B0yFace_z_BLK(i,j,k,globalBLK) = 0.0
           B0zFace_z_BLK(i,j,k,globalBLK) = 0.0
        end do
     end do
  end do

end subroutine set_ICs_arc
