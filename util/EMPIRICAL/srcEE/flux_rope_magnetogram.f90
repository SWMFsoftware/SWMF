!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
program magnetogram
  use EEE_ModMain
  use EEE_ModCommonVariables
  use ModConst, ONLY: cPi, cTwoPi, cDegToRad, cRadToDeg
  use ModPlotFile, ONLY:save_plot_file
  use ModReadParam, ONLY: read_file, read_init, &
         read_line, read_command
  use ModUtilities, ONLYL CON_stop

  implicit none

  integer,parameter:: nLong = 360, nLat = 180
  integer:: Long0, i,j
  real:: Longitude_I(nLong), SinLong_I(nLong), CosLong_I(nLong), LongitudeRad
  real:: Latitude_I(nLat), SinLat_I(nLat), CosLat_I(nLat), LatitudeRad
  real:: B_DC(3,nLong,nLat) = 0
  real::Xyz_D(3),SinLong,CosLong,SinLat,CosLat, Rho, P, U_D(3), B_D(3)
  integer:: iLong
  character(LEN=3)::TypeLatAxis !(sin or uni)
  character(LEN=100):: StringLine, NameCommand
  !----------------
  Io2Si_V = 1; Si2Io_V = 1; Io2No_V = 1 
  No2Io_V = 1; Si2No_V = 1; No2Si_V = 1
  read(*,*)Long0
  write(*,'(a,i3)')prefix,Long0
  read(*,*)TypeLatAxis
  write(*,'(a)')prefix//TypeLatAxis
  do i = 1, nLong
     Longitude_I(i) = i - 0.5
     LongitudeRad = Longitude_I(i)*cDegToRad
     SinLong_I(i) = sin(LongitudeRad)
     CosLong_I(i) = cos(LongitudeRad)
  end do
  select case(TypeLatAxis)
  case('uni')
     do j = 1, nLat
        Latitude_I(j) = j - 90.5
        LatitudeRad = Latitude_I(j)*cDegToRad
        SinLat_I(j) = sin(LatitudeRad)
        CosLat_I(j) = cos(LatitudeRad)
     end do
  case('sin')
     do j = 1, nLat
        SinLat_I(j) = (2.0*j - nLat - 1.0)/nLat
        CosLat_I(j) = sqrt(1.0 - SinLat_I(j)**2)
        LatitudeRad = asin(SinLat_I(j))
        Latitude_I(j) =LatitudeRad*cRadToDeg
     end do
  case default
     call CON_stop('TypeLatAxis='//TypeLatAxis//'... is not supported')
  end select
  call read_file('CME.in')
  call read_init
  READPARAM: do
     if(.not.read_line(StringLine) ) EXIT READPARAM
     if(.not.read_command(NameCommand)) CYCLE READPARAM
     select case(NameCommand)
     case('#END')
        EXIT READPARAM
     case('#CME')
        call EEE_set_parameters(NameCommand)
     case default
        call CON_stop('Unknown command:'//NameCommand)
     end select
  end do READPARAM
  do j = 1,nLat
     SinLat = SinLat_I(j)
     CosLat = CosLat_I(j)
     do i = 1,nLong
        iLong = i + Long0
        iLong = 1 + modulo(iLong -1,nLong)
        CosLong = CosLong_I(iLong)
        SinLong = SinLong_I(iLong)
        Xyz_D(x_) = CosLat*CosLong
        Xyz_D(y_) = CosLat*SinLong
        Xyz_D(z_) = SinLat
        call EEE_get_state_BC(Xyz_D,Rho,U_D,B_D,p,0.0,0,0)
        B_DC(1,i,j) = sum(Xyz_D*B_D)
        B_DC(2,i,j) =  B_D(y_)*CosLong - B_D(x_)*SinLong
        B_DC(3,i,j) = (B_D(x_)*CosLong + B_D(y_)*SinLong)*(-SinLat) &
                     + B_D(z_)*CosLat
     end do
  end do
  call save_plot_file(&
       NameFile = 'FRMagnetogram.out',&
       TypeFileIn = 'ascii',          &
       StringHeaderIn = 'Flux rope magnetic field at 1 R_s',&
       nDimIn = 2,                    &
       ParamIn_I = (/float(Long0)/),         &
       CoordIn_I = Longitude_I,       &
       Coord2In_I= Latitude_I,        &
       NameVarIn = 'Longitude Latitude Br  BPhi BTheta Long0',& 
       VarIn_VII = B_DC)
  stop
end program magnetogram
!============================================================================
