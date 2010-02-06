subroutine trace_dipole(Re,LatStart,nStep,nStepInside,&
     LineLength_I,Bfield_I,RadialDist_I,L)
  use ModCrcmPlanet,  ONLY: MagCoef => dipmom
  implicit none
  
  integer,intent(in) :: nStep, nStepInside
  real,   intent(in) :: Re,LatStart
  real,   intent(inout):: LineLength_I(nStep),Bfield_I(nStep),L,&
                        RadialDist_I(nStep)
  
!  real, parameter :: MagCoef = 7.84e15   ! T m^3
  real, parameter :: rBoundary = 2.5
  real    :: Lat, dLat, Length, LatMax, Length1,LatOuter,LatInner 
  integer :: iStep
  !----------------------------------------------------------------------------

  L  = 1.0 / cos(LatStart)**2.0
  
  ! If no data from MHD solution fill entire space with dipole solution
  if (nStep == 2*nStepInside) then
     dLat = (2.0*LatStart)/nStep
  else
     LatMax   = acos(sqrt(rBoundary/L))
     dLat     = (LatStart-LatMax)/nStepInside
  endif

  Lat = LatStart

  Length = 0.0

  ! Set values inside body to dipole
  do iStep = 1, nStepInside
     Length = calc_dipole_length(L,LatStart,Lat)
     LineLength_I(iStep) = Length
     
     ! Equation of Earth Dipole field strength:
     !B = mu * M /4/pi/r^3 *sqrt[1+3*sin(Lat)^2]
     Bfield_I(iStep) = MagCoef * sqrt(1+3.0*(sin(Lat))**2.0) &
          / (Re*L*cos(Lat)**2.0)**3.0  
     
     RadialDist_I(iStep) = L*cos(Lat)**2.0
     if(iStep < nStepInside) Lat = Lat - dLat
  enddo

  if (nStep > 2*nStepInside) then
     !Loop over part outside boundary
     LatInner = acos(sqrt(RadialDist_I(nStepInside)/L))
     Length1  = calc_dipole_length(L,LatStart,LatInner)
     do iStep = nStepInside+1, nStep-nStepInside
        LineLength_I(iStep) = LineLength_I(iStep) + Length1
     enddo
     LatOuter = -acos(sqrt(RadialDist_I(nStep - nStepInside)/L))
  else
     ! Continue with other half
     LatOuter = Lat
  endif
  
  !Loop over final part of boundary
  Lat      = -abs(Lat)
  do iStep = nStep-nStepInside+1, nStep
     Length = calc_dipole_length(L,LatOuter,Lat)
     LineLength_I(iStep) = Length + LineLength_I(nStep-nStepInside)
     
     ! Equation of Earth Dipole field strength:
     !B = mu * M /4/pi/r^3 *sqrt[1+3*sin(Lat)^2]
     Bfield_I(iStep) = MagCoef * sqrt(1+3.0*(sin(Lat))**2.0) &
          / (Re*L*cos(Lat)**2.0)**3.0  
     
     RadialDist_I(iStep) = L*cos(Lat)**2.0
     
     Lat = Lat - dLat
  enddo
  
contains
  
  !============================================================================
  real function calc_dipole_length(L,LatStart,LatEnd)
    
    implicit none
    real, intent(in) :: L,LatStart,LatEnd
    real :: x1,x2
    !--------------------------------------------------------------------------
    
    x1 = sin(LatStart)
    x2 = sin(LatEnd)
    
    calc_dipole_length =&
         (0.5*Sqrt(3.0*x2**2.0+1.0)*x2 + asinh(sqrt(3.0)*x2)/2.0/sqrt(3.0)) &
         - (0.5*Sqrt(3.0*x1**2.0+1.0)*x1 + asinh(sqrt(3.0)*x1)/2.0/sqrt(3.0)) 
    
    calc_dipole_length = abs(L*calc_dipole_length)
    
  end function calc_dipole_length
  
  !============================================================================
  real function asinh(x)
    
    implicit none
    real, intent(in) :: x
    
    !--------------------------------------------------------------------------
    
    asinh = log( x + sqrt(x**2.0 + 1.0) )
  end function asinh
end subroutine trace_dipole

!==============================================================================

subroutine trace_dipole_test
  implicit none
  real, parameter   :: Re = 6378000.0
  real              :: LatStart = 1.10714871779
  integer,parameter :: nStep = 20, nStepInside=10
!  integer,parameter :: nStep = 20, nStepInside=3
  real              :: LineLength_I(nStep)= 0.0, &
                       Bfield_I(nStep)    = 0.0, &
                       RadialDist_I(nStep)= 0.0
  real              :: L
  integer           :: i
  !----------------------------------------------------------------------------


!  LineLength_I(4:17) = (/1.6605227,2.2915752,2.9301178,3.5609026,4.1732516,&
!       4.7624826 ,5.3307447, 6.4423418,7.0106039,7.5998344, &
!       8.2121840, 8.8429689,9.4815121,10.1125641/)
!  
!  LineLength_I(4:17) = LineLength_I(4:17)-1.6605227
!  
!  RadialDist_I(4:17) = (/2.5519667,3.0996408,3.6180336,4.0818324,4.4683876,&
!       4.7588248,4.9389615,4.9389615,4.7588248,4.4683876,4.0818324,&
!       3.6180336,3.0996408,2.5519667/)
!
!  Bfield_I(4:17) = (/2.8568163E-06,1.4844210E-06,8.6292482E-07,5.5333624E-07,&
!       3.8897971E-07,2.9999231E-07,2.5536775E-07,2.5536775E-07,&
!       2.9999231E-07,3.8897971E-07,5.5333624E-07,8.6292482E-07,1.4844210E-06,&
!       2.8568163E-06/)
  
  call trace_dipole(Re,LatStart,nStep,nStepInside,&
     LineLength_I,Bfield_I,RadialDist_I,L)  
  
  write(*,*) 'Latitude:',LatStart
  write(*,*) 'i,LineLength_I, RadialDist_I,Bfield_I'
  do i=1,nStep
     write(*,*) i, LineLength_I(i), RadialDist_I(i), Bfield_I(i)
  enddo
end subroutine trace_dipole_test
