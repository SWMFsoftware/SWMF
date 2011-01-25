!^CFG COPYRIGHT UM
!==============================================================================
module ModGroundMagPerturb

  use ModPlanetConst,    ONLY: rPlanet_I, Earth_
  use ModPhysics,        ONLY: rCurrents, No2Io_V, Si2No_V, UnitB_, UnitJ_
  use ModCoordTransform, ONLY: sph_to_xyz, rot_xyz_sph, cross_product,xyz_to_sph
  use ModConst,          ONLY: cHalfPi, cDegToRad

  implicit none
  save

  public:: read_mag_input_file
  public:: open_magnetometer_output_file
  public:: close_magnetometer_output_file
  public:: write_magnetometers

  logical, public    :: save_magnetometer_data = .false.

  integer, parameter :: Max_MagnetometerNumber = 100
  integer            :: iUnitMag = -1
  integer            :: nMagnetometer = 1
  real, dimension(2,Max_MagnetometerNumber) :: &
       PosMagnetometer_II
  character(len=50)  :: MagInputFile
  character(len=3)   :: MagName_I(Max_MagnetometerNumber), MagInCoord='MAG'

contains

  !===========================================================================
  subroutine ground_mag_perturb(Xyz_DI, MagPerturb_DI)

    ! This subroutine is used to calculate the ground magnetic perturbations, 
    ! at given point in GM cells.

    use ModSize,           ONLY: nI, nJ, nK, nBLK, gcn
    use ModGeometry,       ONLY: x_BLK, y_BLK, z_BLK, R_BLK, &
         x1, x2, y1, y2, z1, z2
    use ModAdvance,        ONLY: Tmp1_BLK
    use ModMain,           ONLY: x_, y_, z_, r_, phi_, theta_, &
         unusedBLK, nBlock, Time_Simulation
    use ModNumConst,       ONLY: cPi
    use ModCurrent,        ONLY: get_current

    implicit none

    real, intent(in), dimension(3,nMagnetometer) :: Xyz_DI
    real, intent(out),dimension(3,nMagnetometer) :: MagPerturb_DI
    integer  :: i,j,k,iBLK
    real     :: Re = 6378000.00, r3, x, y, z
    real, dimension(3):: Xyz_D,Xyz_BLK, Temp_D,Current_D
    integer  :: iMag
    real, external :: integrate_BLK
    real, allocatable, dimension(:,:,:,:) :: Temp_BLK_x,Temp_BLK_y,Temp_BLK_z

    !--------------------------------------------------------------------

    if(.not.allocated(Temp_BLK_x))&
         allocate(Temp_BLK_x(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn,nBLK), &
         Temp_BLK_y(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn,nBLK), &
         temp_BLK_z(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn,nBLK))

    Temp_BLK_x = 0.0
    Temp_BLK_y = 0.0
    Temp_BLK_z = 0.0
    !\
    ! Calculate the magnetic perturbations in cartesian coordinates
    !/

    do iMag = 1, nMagnetometer
       Xyz_D = Xyz_DI(:,iMag)/rPlanet_I(Earth_)

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
                Temp_BLK_x(i,j,k,iBLK)=0.0
                Temp_BLK_y(i,j,k,iBLK)=0.0
                Temp_BLK_z(i,j,k,iBLK)=0.0
                CYCLE
             end if

             call get_current(i,j,k,iBLK,Current_D)

             Xyz_BLK = (/x_BLK(i,j,k,iBLK),y_BLK(i,j,k,iBLK), &
                  z_BLK(i,j,k,iBLK)/)

             r3 = (sqrt(sum((Xyz_D-Xyz_BLK)**2)))**3

             Temp_D = cross_product(Current_D, Xyz_D-Xyz_BLK)/r3 

             Temp_BLK_x(i,j,k,iBLK) = Temp_D(1)
             Temp_BLK_y(i,j,k,iBLK) = Temp_D(2)
             Temp_BLK_z(i,j,k,iBLK) = Temp_D(3)
          end do; end do; end do
       end do

       MagPerturb_DI(x_,iMag) = Integrate_BLK(1,Temp_BLK_x)/(4*cPi) 
       MagPerturb_DI(y_,iMag) = Integrate_BLK(1,Temp_BLK_y)/(4*cPi)
       MagPerturb_DI(z_,iMag) = Integrate_BLK(1,Temp_BLK_z)/(4*cPi)

    end do

    deallocate(Temp_BLK_x,Temp_BLK_y,Temp_BLK_z)
  end subroutine ground_mag_perturb

  !=====================================================================
  subroutine ground_mag_perturb_fac(Xyz_DI, MagPerturb_DI)

    use ModProcMH,         ONLY: iProc, nProc, iComm
    use ModMain,           ONLY: Time_Simulation
    use CON_planet_field,  ONLY: get_planet_field, map_planet_field
    use ModNumConst,       ONLY: cPi, cTwoPi, cRadToDeg
    use ModConst,          ONLY: cMu
    use ModCurrent,        ONLY: calc_field_aligned_current
    use ModMpi

    implicit none

    real, intent(in)      :: Xyz_DI(3,nMagnetometer)
    real, intent(out)     :: MagPerturb_DI(3,nMagnetometer)


    real, parameter       :: rIonosphere = 1.01725 ! rEarth + iono_height
    integer, parameter    :: nTheta =181, nPhi =181,nCuts = 800

    integer               :: i, j, k, iHemisphere, iError
    integer               :: iTheta,iPhi,iLine,iMag
    real                  :: dR_Trace, Theta, Phi, Jr, r_tmp
    real                  :: dL, dS, dTheta, dPhi ,iLat, SinTheta
    real                  :: b, Fac, bRcurrents,JrRcurrents
    real, dimension(3) :: Xyz_D, b_D, bRcurrents_D, XyzRcurrents_D, XyzTmp_D, &
         j_D, temp_D
    real :: FacRcurrents_II(nTheta,nPhi)
    real :: bRcurrents_DII(3,nTheta,nPhi)
    !------------------------------------------------------------------

    MagPerturb_DI= 0.0

    dTheta = cPi    / (nTheta-1)
    dPhi   = cTwoPi / (nPhi-1)
    dR_Trace = (rCurrents - rIonosphere) / nCuts

    ! Get the current and B at the ionosphere boundary
    call calc_field_aligned_current(nTheta,nPhi,rCurrents, &
         FacRcurrents_II, bRcurrents_DII)

    if(nProc>1)then
       call MPI_Bcast(FacRcurrents_II, nTheta*nPhi,MPI_REAL,0,iComm,iError)
       call MPI_Bcast(bRcurrents_DII, 3*nTheta*nPhi,MPI_REAL,0,iComm,iError)
    end if

    ! only need those currents that are above certain threshold
    where(abs(FacRcurrents_II) * No2Io_V(UnitJ_) < 1.0E-4 &
         .or. abs(FacRcurrents_II) * No2Io_V(UnitJ_) >  1.0e3)&
         FacRcurrents_II = 0.0

    iLine=-1
    do iTheta = 1, nTheta
       Theta = (iTheta-1) * dTheta
       if(iTheta==1 .or. iTheta == nTheta)Theta = 1.0E-4
       SinTheta = sin(Theta)          

       do iPhi=1, nPhi
          Phi = (iPhi-1) * dPhi

          ! if the FAC is under certain threshold, do nothing
          if (FacRcurrents_II(iTheta,iPhi) ==0.0)CYCLE

          iLine = iLine +1
          ! do parallel computation among the processors
          if(mod(iLine,nProc) /= iProc)CYCLE

          call sph_to_xyz(rCurrents,Theta,Phi,XyzRcurrents_D)

          ! extract the field aligned current and B field
          JrRcurrents = FacRcurrents_II(iTheta,iPhi)

          bRcurrents_D= bRcurrents_DII(:,iTheta,iPhi)
          bRcurrents = sqrt(sum(bRcurrents_D**2))

          do k=1,nCuts

             r_tmp = rCurrents - dR_Trace * k
             ! get next position along the field line
             call map_planet_field(Time_Simulation,XyzRcurrents_D,'SMG NORM', &
                  r_tmp, XyzTmp_D,iHemisphere)

             ! get the B field at this position
             call get_planet_field(Time_Simulation, XyzTmp_D,'SMG NORM',B_D)
             B_D = B_D *Si2No_V(UnitB_)
             b = sqrt(sum(B_D**2))

             ! get the field alinged current at this position
             Fac = b/bRcurrents * JrRcurrents
             ! get the (x,y,z) components of the Jr
             j_D = Fac * B_D/b

             ! the length of the field line between two cuts
             iLat = abs(asin(XyzTmp_D(3)/sqrt(sum(XyzTmp_D**2))))
             dL = dR_Trace * sqrt(1+3*(sin(iLat))**2)/(2*sin(iLat))
             ! the cross section area by conversation of magnetic flux
             dS = bRcurrents/ b * rCurrents**2 * SinTheta * dTheta *dPhi

             do iMag=1,nMagnetometer
                Xyz_D = Xyz_DI(:,iMag)/rplanet_I(Earth_)

                if(Xyz_D(3) > 0 .and. Theta > cHalfPi &
                     .or. Xyz_D(3) < 0 .and. Theta < cHalfPi) CYCLE

                ! Do the Biot-Savart integral JxR/|R|^3 dV for all the magnetometers
                temp_D = cross_product(j_D, Xyz_D-XyzTmp_D) & 
                     * dL * dS /(sqrt(sum((XyzTmp_D-Xyz_D)**2)))**3

                MagPerturb_DI(:,iMag)=MagPerturb_DI(:,iMag)+temp_D/(4*cPi)

             end do

          end do
       end do
    end do

  end subroutine ground_mag_perturb_fac
  !================================================================
  subroutine read_mag_input_file

    use ModProcMH, ONLY: iProc, iComm
    use ModMain, ONLY: lVerbose
    use ModIO, ONLY: FileName, Unit_Tmp, iUnitOut, write_prefix
    use ModMpi

    implicit none

    integer :: iError, nStat

    ! One line of input
    character (len=100) :: Line
    character(len=3) :: iMagName
    real             :: iMagmLat, iMagmLon, Xyz_D(3)
    real, dimension(Max_MagnetometerNumber)      :: &
         MagmLat_I, MagmLon_I
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

       nStat = 0

       ! Read the file: read #COORD TypeCoord, #START 
       READFILE: do

          read(unit_tmp,'(a)', iostat = iError ) Line

          if (iError /= 0) EXIT READFILE

          if(index(Line,'#COORD')>0) then
             read(unit_tmp,'(a)') MagInCoord
             select case(MagInCoord)
             case('MAG','SMG')
                call write_prefix;
                write(iUnitOut,'(a)') 'Magnetometer Coordinates= '//MagInCoord
             case default
                call stop_mpi(NameSub//' invalid MagInCoord='//MagInCoord)
             end select
          endif

          if(index(Line,'#START')>0)then
             READPOINTS: do
                read(unit_tmp,*, iostat=iError) iMagName, iMagmLat, iMagmLon
                if (iError /= 0) EXIT READFILE

                if (nStat >= Max_MagnetometerNumber) then
                   call write_prefix;
                   write(*,*) NameSub,' WARNING: magnetometers file: ',&
                        trim(filename),' contains too many stations! '
                   call write_prefix; write(*,*) NameSub, &
                        ': max number of stations =',Max_MagnetometerNumber
                   EXIT READFILE
                endif

                !Add new points
                nStat = nStat + 1

                !Store the locations and name of the stations
                MagmLat_I(nStat)    = iMagmLat
                MagmLon_I(nStat)    = iMagmLon
                MagName_I(nStat)    = iMagName

             end do READPOINTS

          end if

       end do READFILE

       close(unit_tmp)

       if(DoTest)write(*,*) NameSub,': nstat=',nStat

       ! Number of magnetometers 
       nMagnetometer = nStat

       write(*,*) NameSub, ': Number of Magnetometers: ', nMagnetometer

       ! Save the positions (maglatitude, maglongitude)
       do iMag=1, nMagnetometer
          PosMagnetometer_II(1,iMag) = MagmLat_I(iMag)
          PosMagnetometer_II(2,iMag) = MagmLon_I(iMag)
       end do

    end if

    ! Tell the coordinates to the other processors
    call MPI_Bcast(MagInCoord, 3, MPI_CHARACTER, 0, iComm, iError)
    ! Tell the number of magnetometers to the other processors
    call MPI_Bcast(nMagnetometer, 1, MPI_INTEGER, 0, iComm, iError)
    ! Tell the magnetometer name to the other processors
    call MPI_Bcast(MagName_I, nMagnetometer*3, MPI_CHARACTER,0,iComm,iError)
    ! Tell the other processors the coordinates
    call MPI_Bcast(PosMagnetometer_II, 2*nMagnetometer, MPI_REAL, 0, &
         iComm, iError)
    !-----------------------------------------------------------------------
  end subroutine read_mag_input_file

  !=====================================================================
  subroutine open_magnetometer_output_file

    use ModMain,   ONLY: n_step
    use ModIoUnit, ONLY: io_unit_new
    use ModIO,     ONLY: filename, NamePlotDir

    implicit none

    integer :: iMag
    logical :: oktest, oktest_me

    ! Open the output file 
    call set_oktest('open_magnetometer_output_files', oktest, oktest_me)

    write(filename,'(a,i8.8,a)')trim(NamePlotDir)//&
         'GM_mag_n', n_step, '.mag'

    if(oktest) then
       write(*,*) 'open_magnetometer_output_files: filename:', filename
    end if

    iUnitMag= io_unit_new()
    open(iUnitMag, file=filename, status="replace")

    ! Write the header
    write(iUnitMag, '(i5,a)',ADVANCE="NO") nMagnetometer, ' magnetometers: '
    do iMag=1,nMagnetometer-1 
       write(iUnitMag, '(1X,a)', ADVANCE='NO') MagName_I(iMag)
    end do
    write(iUnitMag, '(1X,a)') MagName_I(nMagnetometer)
    write(iUnitMag, '(a)')  &
         'nstep year mo dy hr mn sc msc station X Y Z '// &
         'dBn dBe dBd facdBn facdBe facdBd'

  end subroutine open_magnetometer_output_file

  !=====================================================================
  subroutine write_magnetometers

    use ModProcMH,ONLY: iProc, nProc, iComm
    use CON_axes, ONLY: transform_matrix
    use ModMain,  ONLY: x_, y_, z_, r_, theta_, phi_, n_step, time_simulation,&
         TypeCoordSystem
    use ModUtilities, ONLY: flush_unit
    use ModNumConst, ONLY: cTwoPi
    use ModMpi

    implicit none

    integer           :: iMag, iTime_I(7), iError
    !year,month,day,hour,minute,second,msecond
    real, dimension(3):: Xyz_D, &
         TmpGmSph_D,TmpFacSph_D, &
         MagPerturbGmSph_D, MagPerturbFacSph_D,&
         MagVarFac_D, MagVarGm_D
    real, dimension(3,nMagnetometer):: MagPerturb_DI, & 
         MagPerturbGm_DI, MagPerturbFac_DI, &
         MagGsmXyz_DI, MagSmgXyz_DI,MagVarSum_DI
    real:: XyzSph_DD(3,3), MagtoGsm_DD(3,3), GsmtoSmg_DD(3,3)
    real:: rplanet_tmp, phi_smg,theta_smg
    !--------------------------------------------------------

    ! Matrix between two coordinates
    MagtoGsm_DD = transform_matrix(Time_Simulation, &
         MagInCoord, TypeCoordSystem)
    GsmtoSmg_DD = transform_matrix(Time_Simulation, TypeCoordSystem, 'SMG')

    !\
    ! Transform the Radius position into cartesian coordinates. 
    ! Transform the magnetometer position from MagInCorrd to GSM/SM
    !/

    do iMag=1,nMagnetometer
       ! (360,360) is for the station at the center of the planet
       if ( nint(PosMagnetometer_II(1,iMag)) == 360 .and. &
            nint(PosMagnetometer_II(2,iMag)) == 360) then 
          Xyz_D = 0.0
       else 
          call  sph_to_xyz(rPlanet_I(Earth_),             &
               (90-PosMagnetometer_II(1,iMag))*cDegToRad, &
               PosMagnetometer_II(2,iMag)*cDegToRad,      &
               Xyz_D)
          Xyz_D = matmul(MagtoGsm_DD, Xyz_D)
       end if

       MagGsmXyz_DI(:,iMag) = Xyz_D
       MagSmgXyz_DI(:,iMag) = matmul(GsmtoSmg_DD,Xyz_D)
    end do

    !-------------------------------------------------------------------
    ! Calculate the perturbations in GSM coordinates from GM currents;
    ! Convert Magnetic perturbations into SM coordinates
    !------------------------------------------------------------------

    call ground_mag_perturb(MagGsmXyz_DI, MagPerturb_DI) 

    do iMag=1,nMagnetometer 
       MagPerturbGm_DI(:,iMag) = matmul(GsmtoSmg_DD, MagPerturb_DI(:,iMag))
    end do
    !\
    ! Collect the variables from all the PEs
    !/
    MagVarSum_DI = 0.0
    if(nProc>1)then 
       call MPI_reduce(MagPerturbGm_DI, MagVarSum_DI, 3*nMagnetometer, &
            MPI_REAL, MPI_SUM, 0, iComm, iError)

       if(iProc==0)MagPerturbGm_DI = MagVarSum_DI
    end if

    !-------------------------------------------------------------------
    ! Calculate the perturbations from Field Aligned currents in Gap Region;
    ! It is in SM coordinate.
    !-------------------------------------------------------------------

    call ground_mag_perturb_fac(MagSmgXyz_DI, MagPerturbFac_DI)

    !\
    ! Collect the variables from all the PEs
    !/
    MagVarSum_DI = 0.0
    if(nProc>1)then 
       call MPI_reduce(MagPerturbFac_DI, MagVarSum_DI, 3*nMagnetometer, &
            MPI_REAL, MPI_SUM, 0, iComm, iError)

       if(iProc==0)MagPerturbFac_DI = MagVarSum_DI
    end if

    !-------------------------------------------------------------------
    ! Transform to spherical coordinates (r, phi, theta)
    !-------------------------------------------------------------------
    if(iProc==0)then      
       do iMag=1,nMagnetometer

          if (MagSmgXyz_DI(x_,iMag) == 0.0 .and. MagSmgXyz_DI(y_,iMag) == 0.0 &
               .and. MagSmgXyz_DI(z_,iMag) == 0.0) then 
             ! the magnetometer at the center (0,0,0)
             ! the perturbations are actually in SM cartesian coordinate
             MagPerturbGmSph_D  = MagPerturbGm_DI(:,iMag)
             MagPerturbFacSph_D = MagPerturbFac_DI(:,iMag)
          else

             ! Transform to spherical coordinates (r,theta,phi) components
             XyzSph_DD = rot_xyz_sph(MagSmgXyz_DI(:,iMag))

             TmpGmSph_D = matmul(MagPerturbGm_DI(:,iMag), XyzSph_DD)
             TmpFacSph_D= matmul(MagPerturbFac_DI(:,iMag), XyzSph_DD)

             ! Transform to spherical coordinates (north, east, down) components
             MagPerturbGmSph_D(1)  = -TmpGmSph_D(phi_) 
             MagPerturbGmSph_D(2)  =  TmpGmSph_D(theta_) 
             MagPerturbGmSph_D(3)  = -TmpGmSph_D(r_) 

             MagPerturbFacSph_D(1)  = -TmpFacSph_D(phi_) 
             MagPerturbFacSph_D(2)  =  TmpFacSph_D(theta_) 
             MagPerturbFacSph_D(3)  = -TmpFacSph_D(r_) 

          endif

          !\
          ! Write variables
          !/
          call get_date_time(iTime_I)

          !normalize the variable to I/O unit...
          MagVarGm_D  = MagPerturbGMSph_D * No2Io_V(UnitB_)
          MagVarFac_D = MagPerturbFacSph_D * No2Io_V(UnitB_)


          write(iUnitMag,'(i5)',ADVANCE='NO') n_step
          write(iUnitMag,'(i5,5(1X,i2.2),1X,i3.3)',ADVANCE='NO') iTime_I
          write(iUnitMag,'(1X,i2)', ADVANCE='NO')  iMag

          ! Write position of magnetometer in SGM Coords
          write(iUnitMag,'(3es13.5)',ADVANCE='NO') MagSmgXyz_DI(:,iMag)
          ! Get the Mag_perturb data and Write out
          write(iUnitMag, '(3es13.5)', ADVANCE='NO') MagVarGm_D
          ! Write the FACs' perturbations
          write(iUnitMag, '(3es13.5)') MagVarFac_D

          call flush_unit(iUnitMag)
       end do
    end if

  end subroutine write_magnetometers

  !=====================================================================

  subroutine close_magnetometer_output_file

    implicit none

    close(iUnitMag)

  end subroutine close_magnetometer_output_file
  !============================================================================

end module ModGroundMagPerturb
