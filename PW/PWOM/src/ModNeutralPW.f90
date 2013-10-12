!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModNeutralPW

  integer :: nlon, nlat, nAlt
  Parameter (nlon=73)
  Parameter (nlat=72)
  Parameter (nAlt=40)
  real,dimension(:,:,:),allocatable :: &
       X, Y, Z, Longitude,                    &
       Latitude,Altitude, Rho,                &
       Temperature,DensityO,DensityO2,        &
       DensityN2,DensityNO,DensityN_4S,       &
       DensityN_2D,uNeutralEast,uNeutralNorth,&
       uNeutralUp,DensityOp,DensityNp,        &
       DensityO2p,DensityN2p,DensityNOp,      &
       Densitye,eTemp,iTemp,uIonEast,         &
       uIonNorth, uIonUp,Potential,           &
       uNeutralUpO,uNeutralUpO2,uNeutralUpN2, &
       uNeutralUpN,uExBeast,uExBnorth,        &
       uExBup,ExBmag,IonizationEUV,           &
       IonizationAurora,JouleHeating,HeatingRate

contains

  subroutine init_mod_neutral
    if(allocated(X)) return
    allocate( &
         X(nLon,nLat,nAlt), Y(nLon,nLat,nAlt), Z(nLon,nLat,nAlt), &
         Longitude(nLon,nLat,nAlt), Latitude(nLon,nLat,nAlt), &
         Altitude(nLon,nLat,nAlt), &
         Rho(nLon,nLat,nAlt), Temperature(nLon,nLat,nAlt), &
         DensityO(nLon,nLat,nAlt),DensityO2(nLon,nLat,nAlt), &
         DensityN2(nLon,nLat,nAlt),DensityNO(nLon,nLat,nAlt), &
         DensityN_4S(nLon,nLat,nAlt), DensityN_2D(nLon,nLat,nAlt), &
         uNeutralEast(nLon,nLat,nAlt),uNeutralNorth(nLon,nLat,nAlt), &
         uNeutralUp(nLon,nLat,nAlt), DensityOp(nLon,nLat,nAlt), &
         DensityNp(nLon,nLat,nAlt), DensityO2p(nLon,nLat,nAlt), &
         DensityN2p(nLon,nLat,nAlt), DensityNOp(nLon,nLat,nAlt),      &
         Densitye(nLon,nLat,nAlt), &
         eTemp(nLon,nLat,nAlt), iTemp(nLon,nLat,nAlt), &
         uIonEast(nLon,nLat,nAlt), uIonNorth(nLon,nLat,nAlt), &
         uIonUp(nLon,nLat,nAlt), Potential(nLon,nLat,nAlt), &
         uNeutralUpO(nLon,nLat,nAlt), uNeutralUpO2(nLon,nLat,nAlt), &
         uNeutralUpN2(nLon,nLat,nAlt), uNeutralUpN(nLon,nLat,nAlt), &
         uExBeast(nLon,nLat,nAlt), uExBnorth(nLon,nLat,nAlt),       &
         uExBup(nLon,nLat,nAlt), ExBmag(nLon,nLat,nAlt), &
         IonizationEUV(nLon,nLat,nAlt), IonizationAurora(nLon,nLat,nAlt), &
         JouleHeating(nLon,nLat,nAlt), HeatingRate(nLon,nLat,nAlt))

  end subroutine init_mod_neutral

end Module ModNeutralPW
