!                                                                      C
!**********************************************************************C
!                                                                      C
!  Neutral atmosphere using JGITM 1-D output
!                                                                      C
!**********************************************************************C
!

SUBROUTINE JupiterAtmos (nH2,nH,nH2O,nCH4,Temp)
  use ModParameters
  use ModCommonVariables, ONLY: GRAVTY,ALTD,ALTMIN,ALTMAX,NDIM
  use ModInterpolate, ONLY: linear
  use ModIoUnit, ONLY : UnitTmp_
  implicit none
  
  real, intent(out) :: nH2(0:MaxGrid+1),nH(0:MaxGrid+1),nH2O(0:MaxGrid+1)
  real, intent(out) :: nCH4(0:MaxGrid+1),Temp(0:MaxGrid+1)

  character (len=50), parameter :: DatafileName='PW/JGITM-1D-atmos.dat'
  integer,            parameter :: nSpecies = 4   ! 4 species in file
  integer,            parameter :: nFileAlt = 10000 ! file altitude grid
  real,               parameter :: rPlanet = 71492.0e3 ! m
  real,               parameter :: MPlanet = 1898.3e24 ! kg

  character (len=189) :: line1
  real :: Alt(0:NDIM+1) ! cm
  real :: AtmosArray(5+nSpecies,nFileAlt)
  real :: nHe(0:NDIM+1)
  integer :: iAlt,maxJGITMAlt
  real :: Scaleheight_H2,Scaleheight_H,Scaleheight_H2O,Scaleheight_CH4
  !------------------------------------------------------------------------
  
  open(UnitTmp_,FILE=DatafileName,STATUS='OLD')
  
  read(UnitTmp_,'(a)') line1
  read(UnitTmp_,*) AtmosArray
  
  close(UnitTmp_)

  Alt(0) = ALTMIN
  Alt(NDIM+1) = ALTMAX
  Alt(1:NDIM) = ALTD(1:nDim)

  print *,'Top of JGITM grid: ',AtmosArray(1,nFileAlt) ! km

  do iAlt=0,NDIM+1
     if (Alt(iAlt)/1.0e5.le.AtmosArray(1,nFileAlt)) then
        nH2(iAlt) = linear(AtmosArray(5,:), &        ! H2
             1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
        nHe(iAlt) = linear(AtmosArray(6,:), &        ! He
             1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
        nH(iAlt) = linear(AtmosArray(7,:), &        ! H
             1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
        nCH4(iAlt) = linear(AtmosArray(8,:), &        ! CH4
             1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
        Temp(iAlt) = linear(AtmosArray(2,:), &        ! Temp
             1,nFileAlt,Alt(iAlt)/1e5,AtmosArray(1,:))
        nH2O(iAlt) = 0.0
        maxJGITMAlt = iAlt
     else

        Temp(iAlt) = Temp(maxJGITMAlt)   ! constant temp
        
        ! Scale height in cm
        if (iAlt==nDim+1)then
           Scaleheight_H2 =1.380658e-19*Temp(iAlt)/(3.3452462e-27*GRAVTY(iAlt-1))
           Scaleheight_H  =1.380658e-19*Temp(iAlt)/(1.6726231e-27*GRAVTY(iAlt-1))
           Scaleheight_H2O=Scaleheight_H2
           Scaleheight_CH4=1.380658e-19*Temp(iAlt)/(2.65686432e-26*GRAVTY(iAlt-1))
           !        Scaleheight_H2O=1.380658e-19*Temp(iAlt)/(2.99e-26*GRAVTY(iAlt))
        else
           Scaleheight_H2 =abs(1.380658e-19*Temp(iAlt)/(3.3452462e-27*GRAVTY(iAlt)))
           Scaleheight_H  =abs(1.380658e-19*Temp(iAlt)/(1.6726231e-27*GRAVTY(iAlt)))
           Scaleheight_H2O=Scaleheight_H2
           Scaleheight_CH4=abs(1.380658e-19*Temp(iAlt)/(2.65686432e-26*GRAVTY(iAlt)))
           !        Scaleheight_H2O=1.380658e-19*Temp(iAlt)/(2.99e-26*GRAVTY(iAlt))
        endif
        nH2(iAlt)  = nH2(iAlt-1) *exp(-(Alt(iAlt)-Alt(iAlt-1))/Scaleheight_H2)
        nH(iAlt)   = nH(iAlt-1)  *exp(-(Alt(iAlt)-Alt(iAlt-1))/Scaleheight_H)
        nH2O(iAlt) = nH2O(iAlt-1)*exp(-(Alt(iAlt)-Alt(iAlt-1))/Scaleheight_H2O)
        nCH4(iAlt) = nCH4(iAlt-1)*exp(-(Alt(iAlt)-Alt(iAlt-1))/Scaleheight_CH4)
     endif
  end do

END SUBROUTINE JupiterAtmos
