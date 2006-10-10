
Module ModNeutralPW
  integer :: nlon, nlat, nAlt
  Parameter (nlon=73)
  Parameter (nlat=72)
  Parameter (nAlt=40)
  real,dimension(1:nlon,1:nlat,1:nAlt) ::  X, Y, Z, Longitude,            &
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


  
end Module ModNeutralPW
