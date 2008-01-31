subroutine RB_plot_potential
  use rbe_time,    ONLY: Time=>t
  use rbe_grid,    ONLY: nLat => ir, nLon => ip
  use rbe_cgrid,   ONLY: Lat_I => xlati, Lon_I => phi
  use rbe_convect, ONLY: Potential_II => potent
  use ModNumConst, ONLY: cTwoPi,cPi,cHalfPi
  use ModInterpolate, ONLY: bilinear
  use ModIoUnit, ONLY: UnitTmp_
  implicit none

  character(len=*), parameter :: NameSub='RB_plot_potential'

  !INPUT ARGUMENTS:
  integer :: iLat, iLon,iTime
  real,dimension(:,:),allocatable  :: x_II,y_II,z_II
  Character(len=100) :: NameElectrodynamics
  !----------------------------------------------------------------------------
  !get integer time for filename
  iTime = int(Time)
  write(NameElectrodynamics,"(a,i8.8,a)") &
       'RB/Electrodynamics_Time',iTime,'.dat'
  
  allocate(x_II(nLat,nLon),y_II(nLat,nLon),z_II(nLat,nLon))
  
  open(UnitTmp_,FILE=NameElectrodynamics)
  write(UnitTmp_,*) &
       'VARIABLES = "X", "Y", "Z", "V"'
  
  write(UnitTmp_,*) 'Zone I=', nLon, ', J=', nLat,', DATAPACKING=POINT'
  do iLon=1,nLon
     do iLat=1,nLat
        x_II(iLat,iLon)  =  &
             1.0*sin(cHalfPi-Lat_I(iLat))*cos(Lon_I(iLon))
        y_II(iLat,iLon)  =  &
             1.0*sin(cHalfPi-Lat_I(iLat))*sin(Lon_I(iLon))
        z_II(iLat,iLon)  =  &
             1.0*cos(cHalfPi-Lat_I(iLat))
        write(UnitTmp_,*) &
             x_II(iLat,iLon),y_II(iLat,iLon),z_II(iLat,iLon),Potential_II(iLat,iLon)
     enddo
  enddo
  deallocate(x_II,y_II,z_II)
  
end subroutine RB_plot_potential
