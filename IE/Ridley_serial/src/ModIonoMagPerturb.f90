module ModIonoMagPerturb

  use ModCoordTransform, ONLY: sph_to_xyz, xyz_to_sph, cross_product, rot_xyz_sph
  use ModConst, ONLY: cDegToRad, cMu
  use ModProcIE
  use ModIonosphere
  use IE_ModIO
  use IE_ModMain, ONLY: Time_Simulation, time_array, nSolve
  use ModFiles, ONLY: MagInputFile
  implicit none

  save

  public:: read_mag_input_file
  public:: open_iono_magperturb_file
  public:: close_iono_magperturb_file
  public:: write_iono_magperturb_file

  logical, public    :: save_magnetometer_data = .false.,Initialized_Mag_File=.false.
  integer            :: nMagnetometer = 1
  integer            :: iUnitMag = -1
  integer, parameter :: Max_MagnetometerNumber = 100
  real, dimension(2,Max_MagnetometerNumber) :: &
       PosMagnetometer_II
  character(len=3)   :: MagName_I(Max_MagnetometerNumber), MagInCoord


contains
  !======================================================================
  subroutine iono_mag_perturb(Xyz0_DI, JhMagPerturb_DI, JpMagPerturb_DI)

    use CON_planet_field

    implicit none
    integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi

    real, intent(in), dimension(3,nMagnetometer):: Xyz0_DI
    real, intent(out), dimension(3,nMagnetometer) :: JhMagPerturb_DI, JpMagPerturb_DI
    real, dimension(nTheta*2, nPsi):: Phi, Theta, Psi, ETh, EPs, SigmaH, SigmaP, &
         sintheta, sinphi, costheta, cosphi, Br, Bn, b !,dTheta, dPsi
    real :: dTheta(nTheta*2),dPsi(nPsi)
    real, dimension(nTheta*2, nPsi, 3) :: Jh_DII, Jp_DII, Xyz_DII, eIono_DII
    real, dimension(3) :: bIono_D,  Xyz0_D, MagJh_D, MagJp_D, &
         XyzIono_D, tempJh_dB, tempJp_dB, TempMagJh_D,TempMagJp_D
    real :: sinclats, cosclats, sinlons, coslons, dv
    real :: XyzSph_DD(3,3)

    integer :: i, j, iMag
    !\
    ! calculate the magnetic perturbations at the location of (SMLat, SMLon)
    ! in the SM coordinates, by integrating over the Hall and Perdersen
    ! current systems.
    !/
    JhMagPerturb_DI = 0.0
    JpMagPerturb_DI = 0.0

    Phi(1:nTheta,:)   = Iono_North_PHi
    Theta(1:nTheta,:) = Iono_North_Theta
    Psi(1:nTheta,:)   = Iono_North_Psi
    SigmaH(1:nTheta,:)= Iono_North_SigmaH
    SigmaP(1:nTheta,:)= Iono_North_SigmaP
    dTheta(1:nTheta)  = dTheta_North
    dPsi = dPsi_North

    Phi(nTheta+1:nTheta*2,:)   = Iono_South_PHi
    Theta(nTheta+1:nTheta*2,:) = Iono_South_Theta
    Psi(nTheta+1:nTheta*2,:)   = Iono_South_Psi
    SigmaH(nTheta+1:nTheta*2,:)= Iono_South_SigmaH
    SigmaP(nTheta+1:nTheta*2,:)= Iono_South_SigmaP
    dTheta(nTheta+1:nTheta*2)  = dTheta_South
    dPsi = dPsi_South

    sintheta = sin(Theta)
    sinphi   = sin(Psi)
    costheta = cos(Theta)
    cosphi   = cos(Psi)

    ! dTheta at the poles is 1 degree; the rest are 2 degrees.
    ! dPsi is 4 degrees.
    dTheta(2:nTheta*2-1) = dTheta(2:nTheta*2-1)/2.0
    dPsi = dPsi/2.0

    do j = 1, nPsi
       if ( j<nPsi ) then

          do i = 1, nTheta*2-1
             ETh(i,j) = -(PHI(i+1,j)-PHI(i,j))/                     &
                  (dTheta(i)*Radius)
             EPs(i,j) = -(PHI(i,j+1)-PHI(i,j))/                     &
                  (dPsi(j)*Radius*sinTheta(i,j))
          end do
          ETh(nTheta*2,j)   = ETh(nTheta*2-1,j)
          EPs(nTheta*2,j)   = EPs(nTheta*2-1,j)

       else 
          do i = 1, nTheta*2 -1
             ETh(i,j) = -(PHI(i+1,j)-PHI(i,j))/                     &
                  (dTheta(i)*Radius)
             EPs(i,j) = -(PHI(i,1)-PHI(i,j))/                       &
                  (dPsi(j)*Radius*sinTheta(i,j))
          end do

          ETh(nTheta*2,j)   = ETh(nTheta*2-1,j)
          EPs(nTheta*2,j)   = EPs(nTheta*2-1,j)
       end if

    end do

    ! convert to xyz coords
    eIono_DII(:,:,1) =  ETh*costheta*cosphi - EPs*sinphi
    eIono_DII(:,:,2) =  ETh*costheta*sinphi + EPs*cosphi
    eIono_DII(:,:,3) = -ETh*sintheta

    do i = 1, nTheta*2
       do j = 1, nPsi
          call sph_to_xyz(Radius, Theta(i,j), Psi(i,j), XyzIono_D)
          Xyz_DII(i,j,:) = XyzIono_D
          ! get the magnetic field in SMG coords.
          call get_planet_field(Time_simulation, XyzIono_D, 'SMG', bIono_D)
          bIono_D = bIono_D/sqrt(sum(bIono_D**2))

          ! get the Hall and Perdersen currents in xyz coords
          Jh_DII(i,j,:) = cross_product(bIono_D, eIono_DII(i,j,:)) * SigmaH(i,j)
          Jp_DII(i,j,:) = eIono_DII(i,j,:) * SigmaP(i,j)
       end do
    end do

    do iMag=1, nMagnetometer

       Xyz0_D = Xyz0_DI(:,iMag)

       MagJh_D = 0.0
       MagJp_D = 0.0
       ! Biot-Savart integral to calculate the magnetic perturbations
       if (Xyz0_D(3) < 0) then           
          ! southern hemisphere
          if (iproc/=nproc-1)CYCLE
          do i = nTheta+1, nTheta*2
             do j = 1, nPsi

                tempJh_dB = cross_product(Jh_DII(i,j,:), Xyz0_D-Xyz_DII(i,j,:)) &
                     / (sqrt( sum( (Xyz_DII(i,j,:)-Xyz0_D)**2 )) )**3

                tempJp_dB = cross_product(Jp_DII(i,j,:), Xyz0_D-Xyz_DII(i,j,:)) &
                     / (sqrt( sum( (Xyz_DII(i,j,:)-Xyz0_D)**2 )) )**3

                dv = cMu/(4*cPi) * Radius**2 * dTheta(i) * dPsi(j) * sintheta(i,j)

                MagJh_D = MagJh_D + tempJh_dB * dv
                MagJp_D = MagJp_D + tempJp_dB * dv                

             end do
          end do

       else
          ! northern hemisphere
          if(iProc/=0)CYCLE
          do i = 1, nTheta
             do j = 1, nPsi

                tempJh_dB = cross_product(Jh_DII(i,j,:), Xyz0_D-Xyz_DII(i,j,:)) &
                     / (sqrt(sum((Xyz_DII(i,j,:)-Xyz0_D)**2)))**3

                tempJp_dB = cross_product(Jp_DII(i,j,:), Xyz0_D-Xyz_DII(i,j,:)) &
                     / (sqrt(sum((Xyz_DII(i,j,:)-Xyz0_D)**2)))**3

                dv = cMu/(4*cPi) * Radius**2 * dTheta(i) * dPsi(j) * sintheta(i,j)

                MagJh_D = MagJh_D + tempJh_dB * dv
                MagJp_D = MagJp_D + tempJp_dB * dv

             end do
          end do
       end if


       ! transform to spherical coords (r, theta, phi)
       XyzSph_DD = rot_xyz_sph(Xyz0_D)

       TempMagJh_D = matmul(MagJh_D,XyzSph_DD)
       TempMagJp_D = matmul(MagJp_D,XyzSph_DD)

       ! convert to (north, east, downward) coordinates in unit of nT
       JhMagPerturb_DI(1,iMag) = -TempMagJh_D(2) * 1.0e9 ! north
       JhMagPerturb_DI(2,iMag) =  TempMagJh_D(3) * 1.0e9 ! east
       JhMagPerturb_DI(3,iMag) = -TempMagJh_D(1) * 1.0e9 ! down

       JpMagPerturb_DI(1,iMag) = -TempMagJp_D(2) * 1.0e9 ! north
       JpMagPerturb_DI(2,iMag) =  TempMagJp_D(3) * 1.0e9 ! east
       JpMagPerturb_DI(3,iMag) = -TempMagJp_D(1) * 1.0e9 ! down

    end do


  end subroutine iono_mag_perturb


  !=====================================================================
  subroutine open_iono_magperturb_file

    use ModIoUnit, ONLY: io_unit_new
    implicit none

    integer :: iMag

    ! Open the output file 
    write(*,*) '=> writing magnetic perturbation output.'  
    
    write(NameFile,'(a,3i2.2,"_",3i2.2,a)')trim(NameIonoDir)//"IE_mag_it", &
         mod(time_array(1),100),time_array(2:6), ".dat"

    iUnitMag= io_unit_new()
    open(iUnitMag, file=NameFile)

    ! Write the header
    write(iUnitMag, '(i5,a)') nMagnetometer, ' magnetometers:'
    do iMag=1, nMagnetometer-1
       write(iUnitMag, '(1x,a)', ADVANCE='NO') MagName_I(iMag) 
    end do
    write(iUnitMag, '(1x,a)')MagName_I(nMagnetometer) 
    write(iUnitMag, '(a)')  &
         'nsolve year mo dy hr mn sc msc station X Y Z '// &
         'JhdBn JhdBe JhdBd JpBn JpBe JpBd'
  end subroutine open_iono_magperturb_file

  !======================================================================
  subroutine read_mag_input_file

    use ModMpi

    implicit none

    integer :: iError, nStat, iComm=-1 

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
    filename = MagInputFile

    call write_prefix; write(*,*) NameSub, &
         " reading: ", trim(filename)

    open(unit=iunit, file=filename, status="old", iostat = iError)
    if (iError /= 0) call stop_mpi(NameSub // &
         ' ERROR: unable to open file ' // trim(filename))

    nStat = 0
    ! Read the file: read #COORD TypeCoord, #START                                
    READFILE: do

       read(iunit,'(a)', iostat = iError ) Line

       if (iError /= 0) EXIT READFILE

       if(index(Line,'#COORD')>0) then
          read(iunit,'(a)') Line
          if (index(Line, 'MAG')>0) then
             MagInCoord = 'MAG'
             write(*,*) 'Coordinates= ', MagInCoord
          else if(index(Line, 'SMG')>0)then
             MagInCoord = 'SMG'
             write(*,*) 'Coordinates= ', MagInCoord
          else
             call write_prefix;
             write(*,*) NameSub, &
                  ' WARNING: Please Use', &
                  ' MAG(geomagnetic) coords or SMG coords for Magnetometer InputFile!!'
             EXIT READFILE
          end if
       end if

       if(index(Line,'#START')>0)then
          READPOINTS: do
             read(iunit,*, iostat=iError) iMagName, iMagmLat, iMagmLon
             if (iError /= 0) EXIT READFILE

             !Add new points                                                     
             nStat = nStat + 1

             !Store the locations and name of the stations                       
             MagmLat_I(nStat)    = iMagmLat
             MagmLon_I(nStat)    = iMagmLon
             MagName_I(nStat)    = iMagName

          end do READPOINTS

       end if
    end do READFILE

    close(iunit)

    if(DoTest)write(*,*) NameSub,': nstat=',nStat

    ! Number of magnetometers                    

    nMagnetometer = nStat

    write(*,*) NameSub, ': Number of Magnetometers: ', nMagnetometer

    ! Save the positions (maglatitude, maglongitude)   

    do iMag=1, nMagnetometer
       PosMagnetometer_II(1,iMag) = MagmLat_I(iMag)
       PosMagnetometer_II(2,iMag) = MagmLon_I(iMag)
    end do

    !   if(nProc>1)then
    !      ! Tell the coordinates to the other processors
    !      call MPI_Bcast(MagInCoord, 3, MPI_CHARACTER, 0, iComm, iError)
    !      ! Tell the number of magnetometers to the other processors
    !      call MPI_Bcast(nMagnetometer, 1, MPI_INTEGER, 0, iComm, iError)
    !      ! Tell the magnetometer name to the other processors
    !      call MPI_Bcast(MagName_I, nMagnetometer*3, MPI_CHARACTER,0,iComm,iError)
    !      ! Tell the other processors the coordinates       
    !      call MPI_Bcast(PosMagnetometer_II, 2*nMagnetometer, MPI_REAL, 0, &
    !           iComm, iError)
    !   end if
  end subroutine read_mag_input_file

  !======================================================================
  subroutine write_iono_magperturb_file

    use CON_axes, ONLY: transform_matrix
    use ModMpi
    use ModUtilities, ONLY: flush_unit
    implicit none

    real, dimension(3,nMagnetometer):: Xyz_DI, MagPerturb_Jh_DI, MagPerturb_Jp_DI, &
         MagVarSum_Jh_DI, MagVarSum_Jp_DI
    real, dimension(3):: Xyz_D
    real, dimension(3,3) :: MagtoGsm_DD, GsmtoSmg_DD
    integer :: iMag, iError, i

    MagtoGsm_DD = transform_matrix(Time_Simulation, &
         MagInCoord, 'GSM')
    GsmtoSmg_DD = transform_matrix(Time_Simulation, 'GSM', 'SMG')

    do iMag = 1 , nMagnetometer
       ! (360,360) is for the station at the center of the planet
       if ( nint(PosMagnetometer_II(1,iMag)) == 360 .and. &
            nint(PosMagnetometer_II(2,iMag)) == 360) then
          Xyz_DI(:,iMag) = 0.0
       else
          call  sph_to_xyz(IONO_Radius,             &
               (90-PosMagnetometer_II(1,iMag))*cDegToRad, &
               PosMagnetometer_II(2,iMag)*cDegToRad,      &
               Xyz_D)
          Xyz_D = matmul(MagtoGsm_DD, Xyz_D)
          Xyz_DI(:,iMag) = matmul(GsmtoSmg_DD, Xyz_D)
       end if

    end do
    ! calculate the magnetic perturbation caused by Hall and Perdersen currents
    call iono_mag_perturb(Xyz_DI, MagPerturb_Jh_DI, MagPerturb_Jp_DI)

     !\
     ! Collect the variables from all the PEs
     !/
      MagVarSum_Jh_DI = 0.0
      MagVarSum_Jp_DI = 0.0
      if(nProc>1)then 
         call MPI_reduce(MagPerturb_Jh_DI, MagVarSum_Jh_DI, 3*nMagnetometer, &
              MPI_REAL, MPI_SUM, 0, iComm, iError)
         call MPI_reduce(MagPerturb_Jp_DI, MagVarSum_Jp_DI, 3*nMagnetometer, &
              MPI_REAL, MPI_SUM, 0, iComm, iError)
         if(iProc==0) then
            MagPerturb_Jh_DI = MagVarSum_Jh_DI
            MagPerturb_Jp_DI = MagVarSum_Jp_DI
         end if
      end if
    
     if(iProc==0) then
       do iMag = 1, nMagnetometer
          ! writing
          write(iUnitMag,'(i5)',ADVANCE='NO') nSolve
          write(iUnitMag,'(i5,5(1x,i2.2),1x,i3.3)',ADVANCE='NO') &
               (Time_array(i),i=1,7)
          write(iUnitMag,'(1X,i2)', ADVANCE='NO')  iMag

          ! Write position of magnetometer in SGM Coords
          write(iUnitMag,'(3es13.5)',ADVANCE='NO') Xyz_DI(:,iMag)
          ! Get the magnetic perturb data and Write out
          write(iUnitMag, '(3es13.5)', ADVANCE='NO') MagPerturb_Jh_DI(:,iMag)
          ! Write the Jp perturbations
          write(iUnitMag, '(3es13.5)') MagPerturb_Jp_DI(:,iMag)

          call flush_unit(iUnitMag)
       end do

    end if

  end subroutine write_iono_magperturb_file

  !=====================================================================
  subroutine close_iono_magperturb_file

    implicit none

    close(iUnit)

  end subroutine close_iono_magperturb_file

end module ModIonoMagPerturb
