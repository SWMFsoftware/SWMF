!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
subroutine GetNeutralData

  ! Read in an input file produced by GITM and store the values in a module

  use ModNeutralPW

  Logical       :: IsFirstCall = .true.
  Character*100 :: NameGITM
  integer       :: iUnitGITM
  NameGITM='gitm.dat'
  iUnitGITM=30
  
  If (IsFirstCall) Then
     call init_mod_neutral
     open(iUnitGITM,file=NameGITM)
     do iLon=1,nlon 
        do iLat=1,nLat 
           do iAlt=1,nAlt
              read(unit=iUnitGITM,fmt='(39(1PE14.6))')                        &
                   X(iLon,iLat,iAlt), Y(iLon,iLat,iAlt), Z(iLon,iLat,iAlt),   &
                   Longitude(iLon,iLat,iAlt),Latitude(iLon,iLat,iAlt),        &
                   Altitude(iLon,iLat,iAlt), Rho(iLon,iLat,iAlt),             &
                   Temperature(iLon,iLat,iAlt),DensityO(iLon,iLat,iAlt),      &
                   DensityO2(iLon,iLat,iAlt),DensityN2(iLon,iLat,iAlt),       &
                   DensityNO(iLon,iLat,iAlt),DensityN_4S(iLon,iLat,iAlt),     &
                   DensityN_2D(iLon,iLat,iAlt),uNeutralEast(iLon,iLat,iAlt),  &
                   uNeutralNorth(iLon,iLat,iAlt),uNeutralUp(iLon,iLat,iAlt),  &
                   DensityOp(iLon,iLat,iAlt),DensityNp(iLon,iLat,iAlt),       &
                   DensityO2p(iLon,iLat,iAlt),DensityN2p(iLon,iLat,iAlt),     &
                   DensityNOp(iLon,iLat,iAlt),Densitye(iLon,iLat,iAlt),       &
                   eTemp(iLon,iLat,iAlt),iTemp(iLon,iLat,iAlt),               &
                   uIonEast(iLon,iLat,iAlt),uIonNorth(iLon,iLat,iAlt),        &
                   uIonUp(iLon,iLat,iAlt),Potential(iLon,iLat,iAlt),          &
                   uExBeast(iLon,iLat,iAlt),uExBnorth(iLon,iLat,iAlt),        &
                   uExBup(iLon,iLat,iAlt),ExBmag(iLon,iLat,iAlt),             &
                   IonizationEUV(iLon,iLat,iAlt),                             &
                   IonizationAurora(iLon,iLat,iAlt),                          &
                   JouleHeating(iLon,iLat,iAlt),HeatingRate(iLon,iLat,iAlt),  &
                   uNeutralUpO(iLon,iLat,iAlt),uNeutralUpO2(iLon,iLat,iAlt),  &
                   uNeutralUpN2(iLon,iLat,iAlt), uNeutralUpN(iLon,iLat,iAlt)
           enddo
        enddo
     enddo
     close(iUnitGITM)
     IsFirstCall=.false.
  endif

end subroutine GetNeutralData

!******************************************************************************
! Get_Neutrals searches the Neutral Atmosphere array from GITM and returns
! the Neutral densities and other parameters at the input point
!******************************************************************************

subroutine Get_Neutrals (glatin,glonin,gAltin,Temperature_Out,DensityO_Out,&
                     DensityO2_Out,DensityN2_Out,DensityNO_Out,            &
                     IonizationEUV_Out, IonizationAurora_Out)
  use ModNeutralPW
  
  integer  ::  iLon,iLat,iAlt, iGlat, iGlon,iGAlt

!These are output variables
  real,intent(out)::Temperature_Out,DensityO_Out,DensityO2_Out,        &
                  DensityN2_Out,DensityNO_Out,IonizationEUV_Out,       &
                  IonizationAurora_Out                        !  &
                  !,DensityN_4S_Out,                                    &
                  !DensityN_2D_Out,uNeutralEast_Out,uNeutralNorth_Out,  &
                  !uNeutralUp_Out,uIonEast_Out,                         &
                  !uIonNorth_Out, uIonUp_Out,Potential_Out,             &
                  !uNeutralUpO_Out,uNeutralUpO2_Out,uNeutralUpN2_Out,   &
                  !uNeutralUpN_Out


  real, intent(in) :: glatin,glonin,gAltin
  real             :: glat,glon,gAlt
  
  real    ::  ScaleHeightO, ScaleHeightO2, ScaleHeightN2, ScaleHeightNO 


  !convert input lat and lon to radians and alt to m
  glat=glatin*3.14159265358979/180.0
  glon=glonin*3.14159265358979/180.0
  gAlt=gAltin*1.0e3

  


! Find the associated latitude and longitude index and then extract an 
! altitude line-slice of parameters


  if (glat .le.Latitude(1,floor(nLat/2.0),1)) Then

     if (glat .le. Latitude(1,floor(nLat/4.0),1)) Then
        do iLat=1,floor(nLat/4.0)
           if (glat .le. Latitude(1,iLat,1)) Then
              iGlat=iLat
              exit
           endif
        enddo
     else

        do iLat=floor(nLat/4.0)+1,floor(nLat/2.0)
           if (glat .le. Latitude(1,iLat,1)) Then
              iGlat=iLat
              exit
           endif
        enddo
     endif
  else
     if (glat .le. Latitude(1,floor(3.0*nLat/4.0),1)) Then
        do iLat=floor(nLat/2.0),floor(3.0*nLat/4.0)
           if (glat .le. Latitude(1,iLat,1)) Then
              iGlat=iLat
              exit
           endif
        enddo
     else
        do iLat=floor(3.0*nLat/4.0)+1,nLat
           if (glat .le. Latitude(1,iLat,1)) Then
              iGlat=iLat
              exit
           endif
        enddo
     endif

  endif




  if (glon .le.Longitude(floor(nLon/2.0),1,1)) Then

     if (glon .le. Longitude(floor(nLon/4.0),1,1)) Then
        do iLon=1,floor(nLon/4.0)
           if (glon .le. Longitude(iLon,1,1)) Then
              iGlon=iLon
              exit
           endif
        enddo
     else

        do iLon=floor(nLon/4.0)+1,floor(nLon/2.0)
           if (glon .le. Longitude(iLon,1,1)) Then
              iGlon=iLon
              exit
           endif
        enddo
     endif
  else
     if (glon .le. Longitude(floor(3.0*nLon/4.0),1,1)) Then
        do iLon=floor(nLon/2.0),floor(3.0*nLon/4.0)
           if (glon .le. Longitude(iLon,1,1)) Then
              iGlon=iLon
              exit
           endif
        enddo
     else
        do iLon=floor(3.0*nLon/4.0)+1,nLon
           if (glon .le. Longitude(iLon,1,1)) Then
              iGlon=iLon
              exit
           endif
        enddo
     endif

  endif

  if (gAlt .le.Altitude(1,1,floor(nAlt/2.0))) Then

     if (gAlt .le. Altitude(1,1,floor(nAlt/4.0))) Then
        do iAlt=1,floor(nAlt/4.0)
           if (gAlt .le. Altitude(1,1,iAlt)) Then
              iGAlt=iAlt
              exit
           endif
        enddo
     else

        do iAlt=floor(nAlt/4.0)+1,floor(nAlt/2.0)
           if (gAlt .le. Altitude(1,1,iAlt)) Then
              iGAlt=iAlt
              exit
           endif
        enddo
     endif
  else if (gAlt .le. Altitude(1,1,nAlt)) Then
     if (gAlt .le. Altitude(1,1,floor(3.0*nAlt/4.0))) Then
        do iAlt=floor(nAlt/2.0),floor(3.0*nAlt/4.0)
           if (gAlt .le. Altitude(1,1,iAlt)) Then
              iGAlt=iAlt
              exit
           endif
        enddo
     else
        do iAlt=floor(3.0*nAlt/4.0)+1,nAlt
           if (gAlt .le. Altitude(1,1,iAlt)) Then
              iGAlt=iAlt
              exit
           endif
        enddo
     endif
  else
     iGAlt=nAlt+1
  endif


 

  
!  write(*,*) 'lat', iGlat, Latitude(1,iGlat,1)*180.0/3.14159265358979, &
!       Latitude(iGlon,iGlat,1)*180.0/3.14159265358979&
!       ,gLat*180.0/3.14159265358979
!  
!  write(*,*) 'lon', iGlon, Longitude(iGlon,1,1)*180.0/3.14159265358979, &
!       Longitude(iGlon,iGlat,1)*180.0/3.14159265358979,&
!       gLon*180.0/3.14159265358979
!  
!  write(*,*) 'Alt', iGAlt, Altitude(1,1,iGAlt), &
!       Altitude(iGlon,iGlat,iGAlt),&
!       gAlt

! convert output densities to cm^-3
  
  If (iGAlt .le. nAlt) then
     Temperature_Out=Temperature(iGlon,iGlat,iGAlt)  
     DensityO_Out   = DensityO  (iGlon,iGlat,iGAlt)*1.0e-6 
     DensityO2_Out  = DensityO2 (iGlon,iGlat,iGAlt)*1.0e-6
     DensityN2_Out  = DensityN2 (iGlon,iGlat,iGAlt)*1.0e-6
     DensityNO_Out  = DensityNO (iGlon,iGlat,iGAlt)*1.0e-6
     IonizationEUV_Out = IonizationEUV(iGlon,iGlat,iGAlt) 
     IonizationAurora_Out= IonizationAurora(iGlon,iGlat,iGAlt)
  else
     Temperature_Out=Temperature(iGlon,iGlat,nAlt)  
     
     ScaleHeightO   = 1.38e-23*Temperature_Out/(2.67e-26*9.8) 
     ScaleHeightO2  = 1.38e-23*Temperature_Out/(5.31e-26*9.8) 
     ScaleHeightN2  = 1.38e-23*Temperature_Out/(4.65e-26*9.8) 
     ScaleHeightNO  = 1.38e-23*Temperature_Out/(4.98e-26*9.8) 

     DensityO_Out   = DensityO(iGlon,iGlat,nAlt)*1.0e-6&
                      * exp((Altitude(iGlon,iGlat,nAlt)-GAlt)/ScaleHeightO)
     DensityO2_Out  = DensityO2(iGlon,iGlat,nAlt)*1.0e-6&
                      * exp((Altitude(iGlon,iGlat,nAlt)-GAlt)/ScaleHeightO2)
     DensityN2_Out  = DensityN2(iGlon,iGlat,nAlt)*1.0e-6&
                      * exp((Altitude(iGlon,iGlat,nAlt)-GAlt)/ScaleHeightN2)
     DensityNO_Out  = DensityNO(iGlon,iGlat,nAlt)*1.0e-6&
                      * exp((Altitude(iGlon,iGlat,nAlt)-GAlt)/ScaleHeightN0)
 
     IonizationEUV_Out = IonizationEUV(iGlon,iGlat,nAlt) &
          * exp((Altitude(iGlon,iGlat,nAlt)-GAlt)/ScaleHeightN0)

     IonizationAurora_Out= IonizationAurora(iGlon,iGlat,nAlt)&
          * exp((Altitude(iGlon,iGlat,nAlt)-GAlt)/ScaleHeightN0)


!     write(*,*) ScaleHeightO, DensityO_Out,DensityO(iGlon,iGlat,nAlt),&
!          GAlt,Altitude(iGlon,iGlat,nAlt)-GAlt,exp((Altitude(iGlon,iGlat,nAlt)-GAlt)/ScaleHeightO)
  endif


  


end subroutine Get_Neutrals
