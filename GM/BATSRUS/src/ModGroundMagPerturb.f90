!^CFG COPYRIGHT UM                                                                               
!==============================================================================                   
module ModGroundMagPerturb

  use ModPlanetConst
  use ModPhysics
  use ModCoordTransform, ONLY:sph_to_xyz
  
  implicit none
  save

  public:: read_mag_input_file
  public:: open_magnetometer_output_file
  public:: close_magnetometer_output_file
  public:: write_magnetometers

  logical, public    :: save_magnetometer_data = .false.
    
  integer, parameter :: Max_MagnetometerNumber = 20
  integer            :: iUnitMag_I = -1
  integer            :: nMagnetometer = 0
  real, dimension(2,Max_MagnetometerNumber) :: Pos_Magnetometer, NewPos_Magnetometer
  character(len=50)  :: MagInputFile
  character(len=3)   :: Mag_name(Max_MagnetometerNumber), MagInCoord


contains

!==================================================
!\
! This subroutine is used to calculate the ground magnetic perturbations, 
! at given point in GM cells.
!/
subroutine ground_mag_perturb(Xyz_D, MagPerturb_D)
 
  use ModSize,           ONLY: nI, nJ, nK, nBLK, gcn
  use ModGeometry,       ONLY: x_BLK, y_BLK, z_BLK, R_BLK, &
                               x1, x2, y1, y2, z1, z2
  use ModAdvance,        ONLY: Tmp1_BLK
  use ModMain,           ONLY: x_, y_, z_, unusedBLK, nBlock, Time_Simulation, &
       r_,phi_,theta_
  use ModNumConst,       ONLY: cPi
  use Con_axes
  implicit none
 
  integer  :: i,j,k,iBLK
  real     :: Re = 6378000.00
  real     :: Xyz_D(3)
  real     :: x, y, z
  real     :: MagPerturb_D(3), Current_D(3)
  real, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn, nBLK) :: &
       temp_BLK_x, temp_BLK_y, temp_BLK_z
  real, external :: integrate_BLK

  !--------------------------------------------------------------------
  !\
  ! Calculate the magnetic perturbations in cartesian coordinates
  !/

  
  x = Xyz_D(1)/rplanet_I(Earth_)
  y = Xyz_D(2)/rPlanet_I(Earth_)
  z = Xyz_D(3)/rPlanet_I(Earth_)

  do iBLK=1, nBlock
     if (unusedBLK(iBLK))cycle
     do k=1, nK; do j=1, nJ; do i=1, nI
        if ( r_BLK(i,j,k,iBLK) < rCurrents .or. &
             x_BLK(i+1,j,k,iBLK) > x2 .or.      &
             x_BLK(i-1,j,k,iBLK) < x1 .or.      &
             y_BLK(i,j+1,k,iBLK) > y2 .or.      &
             y_BLK(i,j-1,k,iBLK) < y1 .or.      &
             z_BLK(i,j,k+1,iBLK) > z2 .or.      &
             z_BLK(i,j,k-1,iBLK) < z1 ) then
              temp_BLK_x(i,j,k,iBLK)=0.0
              temp_BLK_y(i,j,k,iBLK)=0.0
              temp_BLK_z(i,j,k,iBLK)=0.0
              CYCLE
           end if
           
           call get_current(i,j,k,iBLK,Current_D)


           temp_BLK_x(i,j,k,iBLK) = ( &
                - current_D(y_)* (z_BLK(i,j,k,iBLK)-z)  &
                + current_D(z_)* (y_BLK(i,j,k,iBLK)-y) ) &
                / (sqrt( &
                (z_BLK(i,j,k,iBLK)-z)**2 + &
                (y_BLK(i,j,k,iBLK)-y)**2 + &
                (x_BLK(i,j,k,iBLK)-x)**2 ))**3
           
           temp_BLK_y(i,j,k,iBLK) = ( &
                - current_D(z_)* (x_BLK(i,j,k,iBLK)-x)   &
                + current_D(x_)* (z_BLK(i,j,k,iBLK)-z) ) & 
                / (sqrt( &
                (z_BLK(i,j,k,iBLK)-z)**2 + &
                (y_BLK(i,j,k,iBLK)-y)**2 + &
                (x_BLK(i,j,k,iBLK)-x)**2 ))**3
           
           temp_BLK_z(i,j,k,iBLK) = ( &
                - current_D(x_)* (y_BLK(i,j,k,iBLK)-y)   &
                + current_D(y_)* (x_BLK(i,j,k,iBLK)-x) ) & 
                / (sqrt( &
                (z_BLK(i,j,k,iBLK)-z)**2 + &
                (y_BLK(i,j,k,iBLK)-y)**2 + &
                (x_BLK(i,j,k,iBLK)-x)**2 ))**3
        end do; end do; end do
     end do

     MagPerturb_D(x_) = Integrate_BLK(1, temp_BLK_x) / (4*cPi)
     MagPerturb_D(y_) = Integrate_BLK(1, temp_BLK_y) / (4*cPi)
     MagPerturb_D(z_) = Integrate_BLK(1, temp_BLK_z) / (4*cPi)

end subroutine ground_mag_perturb

!=====================================================================
subroutine read_mag_input_file

  use ModProcMH
  use ModMain
  use ModIO
  use CON_axes
  use ModMpi
  use ModKind
  use ModIoUnit  
 
  implicit none

  integer :: iError, nstat
  
  ! One line of input
  character (len=100) :: line
  character(len=3) :: iMag_name
  real             :: iMag_mlat, iMag_mlon, Xyz_D(3)
  real, dimension(Max_MagnetometerNumber)      :: &
       Mag_glat, Mag_glon, Mag_mlat, Mag_mlon
  integer          :: iMag
  character(len=*), parameter :: NameSub = 'read_magnetometer_input_files'
  logical          :: DoTest, DoTestMe
!---------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)
  
  ! Read file on the root processor
  if (iProc == 0) then
        
     filename = MagInputFile    
        
     if(lVerbose>0)then
        call write_prefix; write(iUnitOut,*) NameSub, &
             " reading: ",trim(filename)
     end if
               
     open(unit_tmp, file=filename, status="old", iostat = iError)
     if (iError /= 0) call stop_mpi(NameSub // &
          ' ERROR: unable to open file ' // trim(filename))

     nstat = 0
        
     ! Read the file: read #COOR TypeCoord, #START 
     READFILE: do
           
        read(unit_tmp,'(a)', iostat = iError ) line
           
        if (iError /= 0) EXIT READFILE
           
        if(index(line,'#COOR')>0) then
           read(unit_tmp,'(a)') line
           if (index(line, 'MAG')>0) then
              MagInCoord = 'MAG'
              write(*,*) 'Coordinates= ', MagInCoord
           else if(index(line, 'SMG')>0)then
              MagInCoord = 'SMG'
              write(*,*) 'Coordinates= ', MagInCoord
           else
              call write_prefix;
              write(*,*) NameSub, &
                   ' WARNING: Please Use', &
                   ' MAG(geomagnetic coordinates) for SMG Magnetometer InputFile!!'
              EXIT READFILE
           end if
        end if

        if(index(line,'#START')>0)then
           READPOINTS: do
              read(unit_tmp,*, iostat=iError) iMag_name, iMag_mlat, iMag_mlon
              if (iError /= 0) EXIT READFILE
      
              if (nstat >= Max_MagnetometerNumber) then
                 call write_prefix;
                 write(*,*) NameSub,' WARNING: magnetometers file: ',&
                      trim(filename),' contains too many stations! '
                 call write_prefix; write(*,*) NameSub, &
                      ': max number of stations =',Max_MagnetometerNumber
                 EXIT READFILE
              endif
                 
              !Add new points
              nstat = nstat + 1
                 
              !Store the locations and name of the stations
              Mag_mlat(nstat)    = iMag_mlat
              Mag_mlon(nstat)    = iMag_mlon
              Mag_name(nstat)    = iMag_name

           end do READPOINTS
           
        end if
        
     end do READFILE
        
     close(unit_tmp)
        
     if(DoTest)write(*,*) NameSub,': nstat=',nstat
        
     ! Number of magnetometers 
     nMagnetometer = nstat

     write(*,*) NameSub, ': Number of Magnetometers: ', nMagnetometer
 
     ! Save the positions (maglatitude, maglongitude)
     do iMag=1, nMagnetometer
        Pos_Magnetometer(1,iMag) = Mag_mlat(iMag)
        Pos_Magnetometer(2,iMag) = Mag_mlon(iMag)
     end do

  end if
     
  ! Tell the coordinates to the other processors
  call MPI_Bcast(MagInCoord, 3, MPI_CHARACTER, 0, iComm, iError)
  ! Tell the number of magnetometers to the other processors
  call MPI_Bcast(nMagnetometer, 1, MPI_INTEGER, 0, iComm, iError)
  ! Tell the magnetometer name to the other processors
  call MPI_Bcast(Mag_name, nMagnetometer*3, MPI_CHARACTER,0,iComm,iError)
  ! Tell the other processors the coordinates
  call MPI_Bcast(Pos_Magnetometer, 2*nMagnetometer, MPI_REAL, 0, &
       iComm, iError)
  !-----------------------------------------------------------------------
end subroutine read_mag_input_file

!=====================================================================
subroutine open_magnetometer_output_file

  use ModMain,   ONLY: n_step
  use ModIoUnit, ONLY: io_unit_new
  use ModIO,     ONLY: filename, NamePlotDir!, MagOutputFile

  implicit none

 
  logical :: oktest, oktest_me

! Open the output file 
  call set_oktest('open_magnetometer_output_files', oktest, oktest_me)
  
!  filename = MagOutputFile

  write(filename,'(a,i6.6,a)')trim(NamePlotDir)//&
       'mag_n',n_step,'.dat'

  if(oktest) then
     write(*,*) 'open_magnetometer_output_files: filename:', filename
  end if

  iUnitMag_I = io_unit_new()
  open(iUnitMag_I, file=filename)

  ! Write the header
  write(iUnitMag_I, '(i5,a)') nMagnetometer, ' magnetometers'
  write(iUnitMag_I, '(a)')  &
       'year mo dy hr mn sc msc station X Y Z dBn dBe dBd facdBn facdBe facdBn'

end subroutine open_magnetometer_output_file

!=====================================================================
subroutine write_magnetometers
  
  use ModProcMH
  use CON_axes
  use ModMain

  implicit none

  integer:: iMag
  real   :: MagIntoGsm_DD(3,3), GsmtoSmg_DD(3,3)
  character(len=3)  :: MagSystemCoord

  MagSystemCoord = TypeCoordSystem
  MagIntoGsm_DD = transform_matrix(Time_Simulation, MagInCoord, MagSystemCoord)
  
  ! Matrix between the two coordinates
  GsmtoSmg_DD = transform_matrix(Time_Simulation,MagSystemCoord,'SMG')
  
  ! Write out the Magnetometer Data at Each Output Time
  do iMag = 1, nMagnetometer
     call write_magnetometer(iMag)
  end do

  
contains

!=====================================================================
  subroutine write_magnetometer(iMag)
  
    use ModProcMH
    use ModMain
    use ModIO
    use ModMpi
    use ModUtilities
    use ModNumConst, ONLY: cTwoPi
    
    implicit none
    integer           :: n_=1, e_=2, d_=3
    integer           :: iMag, iTime_I(7) 
    !year,month,day,hour,minute,second,msecond
    real, dimension(3):: Xyz_D,Xyz_Magnetometer, MagPerturb_D, &
         MagVar_gm, MagVarSum, MagVar_fac, &
         MagPerturb_D_sph_gm,  MagPerturb_D_gm, tmp_d_sph, &
         MagPerturb_D_sph_fac, MagPerturb_D_fac, fac_tmp_d_sph
    integer           :: iError
    character(len=3)  :: station_name, MagOutCoord
    real              :: rPlanet_Mag, Phi_smg,Theta_smg, rPlanet_tmp

    !--------------------------------------------------------
    
    station_name= Mag_name(iMag)

    !\
    ! Transform the Radius position into cartesian coordinates. 
    ! Transform the magnetometer position from MagInCorrd to GSM
    !/

    ! (360,360) is for the station at the center of the planet
    if (Pos_Magnetometer(1,iMag) == 360 .and. &
         Pos_Magnetometer(2,iMag) == 360) then 
       rPlanet_Mag = 0.0
       call sph_to_xyz(rPlanet_Mag,                  &
            (90-Pos_Magnetometer(1,iMag))*cDegToRad, &
            Pos_Magnetometer(2,iMag)*cDegToRad,      &
            Xyz_D)
    else 
       rPlanet_Mag = rPlanet_I(Earth_) !not in normalized unit
       call  sph_to_xyz(rPlanet_Mag,                 &
            (90-Pos_Magnetometer(1,iMag))*cDegToRad, &
            Pos_Magnetometer(2,iMag)*cDegToRad,      &
            Xyz_D)
    end if

    !-------------------------------------------------------------------
    ! Calculate the perturbations in GSM coordinates from GM currents;
    ! Convert Magnetic perturbations into SMG coordinates
    !-------------------------------------------------------------------
    Xyz_Magnetometer = matmul(MagIntoGsm_DD, Xyz_D)
    call ground_mag_perturb(Xyz_Magnetometer, MagPerturb_D) 
    MagPerturb_D_gm = matmul(GsmtoSmg_DD, MagPerturb_D)

    !\
    ! Collect the variables from all the PEs
    !/
    MagVarSum = 0.0
    if(nProc>1)then 
       call MPI_reduce(MagPerturb_D_gm, MagVarSum, 3, MPI_REAL, MPI_SUM, 0, &
            iComm, iError)
       
       if(iProc==0)MagPerturb_D_gm = MagVarSum
    end if

    !-------------------------------------------------------------------
    ! Calculate the perturbations from Field Aligned currents in Gap Region;
    ! It is in SMG coordinate. This calculation only works from PE 0.
    !-------------------------------------------------------------------
    Xyz_D = matmul(GsmtoSmg_DD, Xyz_Magnetometer)

    call ground_mag_perturb_fac(Xyz_D, MagPerturb_D_fac)
       
    !\
    ! Collect the variables from all the PEs
    !/
    MagVarSum = 0.0
    if(nProc>1)then 
       call MPI_reduce(MagPerturb_D_fac, MagVarSum, 3, MPI_REAL, MPI_SUM, 0, &
            iComm, iError)
       
       if(iProc==0)MagPerturb_D_fac = MagVarSum
    end if
    
    !-------------------------------------------------------------------
    !\
    ! Transform to spherical coordinates (r, phi, theta)
    !/
    !-------------------------------------------------------------------
    if(iProc==0)then      
       if (Xyz_D(x_) == 0.0 .and. Xyz_D(y_) == 0.0 .and. Xyz_D(z_) == 0.0) then 
          ! the magnetometer at the center (0,0,0)
          ! the perturbations are actually in SMG cartesian coordinate
          MagPerturb_D_sph_gm  = MagPerturb_D_gm
          MagPerturb_D_sph_fac = MagPerturb_D_fac
       else
          
          call xyz_to_sph(Xyz_D(x_),Xyz_D(y_),Xyz_D(z_), rPlanet_tmp, Theta_smg, Phi_smg)
          
          ! the transform for the GM currents'perturbations
          tmp_D_sph(r_)   = Magperturb_D_gm(x_)*sin(Theta_smg)*cos(Phi_smg) &
               + Magperturb_D_gm(y_)*sin(Theta_smg)*sin(Phi_smg) &
               + Magperturb_D_gm(z_)*cos(Theta_smg) 
                    
          tmp_D_sph(phi_) = Magperturb_D_gm(x_)*cos(Theta_smg)*cos(Phi_smg) &
               + Magperturb_D_gm(y_)*cos(Theta_smg)*sin(Phi_smg) &
               - Magperturb_D_gm(z_)*sin(Theta_smg) 
          
          tmp_D_sph(theta_) = -Magperturb_D_gm(x_)*sin(Phi_smg) &
               + Magperturb_D_gm(y_)*cos(Phi_smg) 
          
          ! the transform for the FACs' perturbations
          fac_tmp_D_sph(r_)   = Magperturb_D_fac(x_)*sin(Theta_smg)*cos(Phi_smg) &
               + Magperturb_D_fac(y_)*sin(Theta_smg)*sin(Phi_smg) &
               + Magperturb_D_fac(z_)*cos(Theta_smg)

          fac_tmp_D_sph(phi_) = Magperturb_D_fac(x_)*cos(Theta_smg)*cos(Phi_smg) &
               + Magperturb_D_fac(y_)*cos(Theta_smg)*sin(Phi_smg) &
               - Magperturb_D_fac(z_)*sin(Theta_smg)
          
          fac_tmp_D_sph(theta_) = -Magperturb_D_fac(x_)*sin(Phi_smg) &
               + Magperturb_D_fac(y_)*cos(Phi_smg)


          !\
          ! Transform to spherical coordinates (north, east, down)
          !/
          MagPerturb_D_sph_gm(n_)  = -tmp_D_sph(phi_) 
          MagPerturb_D_sph_gm(e_)  =  tmp_D_sph(theta_) 
          MagPerturb_D_sph_gm(d_)  = -tmp_D_sph(r_) 

          MagPerturb_D_sph_fac(n_)  = -fac_tmp_D_sph(phi_) 
          MagPerturb_D_sph_fac(e_)  =  fac_tmp_D_sph(theta_) 
          MagPerturb_D_sph_fac(d_)  = -fac_tmp_D_sph(r_) 
         
       endif
       
       
    !   call flush_unit(6)
    !   call MPI_BARRIER(iComm,iError) ! ----------- BARRIER ------
       
       MagVar_gm = MagPerturb_D_sph_gm
       MagVar_fac = MagPerturb_D_sph_fac
       
       !\
       ! Write variables
       !/
       call get_date_time(iTime_I)
       
       !normalize the variable to I/O unit...
       MagVar_gm  = MagVar_gm  * No2Io_V(UnitB_)
       MagVar_fac = MagVar_fac * No2Io_V(UnitB_)
       ! Write date time
       write(iUnitMag_I,'(i5)',ADVANCE='NO') n_step
       write(iUnitMag_I,'(i5,5(1X,i2.2),1X,i3.3)',ADVANCE='NO') iTime_I
       write(iUnitMag_I,'(1X,a)', ADVANCE='NO')  station_name
       ! Write position of magnetometer in SGM Coords
       write(iUnitMag_I,'(3es13.5)',ADVANCE='NO') Xyz_D
       ! Get the Mag_perturb data and Write out
       write(iUnitMag_I, '(3es13.5)', ADVANCE='NO') MagVar_gm
       ! Write the FACs' perturbations
       write(iUnitMag_I, '(3es13.5)') MagVar_fac

       call flush_unit(iUnitMag_I)
       
    end if

  end subroutine write_magnetometer
  !==================================================================  
end subroutine write_magnetometers
!=====================================================================

subroutine close_magnetometer_output_file

  implicit none
  
  close(iUnitMag_I)
  
end subroutine close_magnetometer_output_file
  !============================================================================
  !^CFG END RAYTRACE

end module ModGroundMagPerturb
