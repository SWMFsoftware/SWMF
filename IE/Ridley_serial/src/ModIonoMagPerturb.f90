!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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

  logical :: IsInitiated = .false.
  integer :: nMagnetometer = 0
  real,             allocatable :: PosMagnetometer_II(:,:)
  character(len=3), allocatable :: TypeCoordMag_I(:)

contains 
  !======================================================================
  subroutine iono_mag_init
    ! Initialize ionospheric magnetometers by allocating arrays.
    !--------------------------------------------------------------------

    IsInitiated = .true.
    if(.not.allocated(PosMagnetometer_II)) allocate( &
         PosMagnetometer_II(2, nMagnetometer), TypeCoordMag_I(nMagnetometer))

  end subroutine iono_mag_init

  !======================================================================
  subroutine iono_mag_perturb(nMag, Xyz0_DI, JhMagPerturb_DI, JpMagPerturb_DI)
    ! For a series of nMag virtual observatories at SMG coordinates Xyz0_DI, 
    ! calculate the magnetic pertubation from the ionospheric Pederson currents
    ! (JpMagPerturb_DI) and Hall currents (JhMagPerturb_DI) in three orthogonal
    ! directions.

    use CON_planet_field, ONLY: get_planet_field

    implicit none

    integer,intent(in)                     :: nMag
    real,   intent(in),  dimension(3,nMag) :: Xyz0_DI
    real,   intent(out), dimension(3,nMag) :: JhMagPerturb_DI, JpMagPerturb_DI

    integer, parameter :: nTheta = IONO_nTheta, nPsi = IONO_nPsi

    real, dimension(nTheta*2, nPsi, 3) :: Jh_IID, Jp_IID, Xyz_IID, eIono_IID
    real, dimension(nTheta*2, nPsi)    :: Phi, Theta, Psi, ETh, EPs, &
         SigmaH, SigmaP, SinTheta, SinPhi, CosTheta, CosPhi
    real, dimension(3)                 :: bIono_D,  Xyz0_D, MagJh_D, MagJp_D, &
         XyzIono_D, tempJh_dB, tempJp_dB, TempMagJh_D,TempMagJp_D

    real :: dTheta(nTheta*2),dPsi(nPsi)
    real :: dv
    real :: XyzSph_DD(3,3)
    integer :: i, j, iMag
    !======================================================================
    call timing_start('iono_mag_perturb')

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

    SinTheta = sin(Theta)
    SinPhi   = sin(Psi)
    CosTheta = cos(Theta)
    CosPhi   = cos(Psi)

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
                  (dPsi(j)*Radius*SinTheta(i,j))
          end do
          ETh(nTheta*2,j)   = ETh(nTheta*2-1,j)
          EPs(nTheta*2,j)   = EPs(nTheta*2-1,j)

       else 
          do i = 1, nTheta*2 -1
             ETh(i,j) = -(PHI(i+1,j)-PHI(i,j))/                     &
                  (dTheta(i)*Radius)
             EPs(i,j) = -(PHI(i,1)-PHI(i,j))/                       &
                  (dPsi(j)*Radius*SinTheta(i,j))
          end do

          ETh(nTheta*2,j)   = ETh(nTheta*2-1,j)
          EPs(nTheta*2,j)   = EPs(nTheta*2-1,j)
       end if

    end do

    ! convert to xyz coords
    eIono_IID(:,:,1) =  ETh*CosTheta*CosPhi - EPs*SinPhi
    eIono_IID(:,:,2) =  ETh*CosTheta*SinPhi + EPs*CosPhi
    eIono_IID(:,:,3) = -ETh*SinTheta

    do i = 1, nTheta*2
       do j = 1, nPsi
          call sph_to_xyz(Radius, Theta(i,j), Psi(i,j), XyzIono_D)
          Xyz_IID(i,j,:) = XyzIono_D
          ! get the magnetic field in SMG coords.
          call get_planet_field(Time_simulation, XyzIono_D, 'SMG', bIono_D)
          bIono_D = bIono_D/sqrt(sum(bIono_D**2))

          ! get the Hall and Perdersen currents in xyz coords
          Jh_IID(i,j,:) = cross_product(bIono_D, eIono_IID(i,j,:))*SigmaH(i,j)
          Jp_IID(i,j,:) = eIono_IID(i,j,:) * SigmaP(i,j)
       end do
    end do

    call timing_start('iono_mag_db')

    do iMag=1, nMag

       Xyz0_D = Xyz0_DI(:,iMag)

       MagJh_D = 0.0
       MagJp_D = 0.0
       ! Biot-Savart integral to calculate the magnetic perturbations
       if (Xyz0_D(3) < 0) then           
          ! southern hemisphere
          if (iProc /= nProc-1)CYCLE
          do i = nTheta+1, nTheta*2
             do j = 1, nPsi

                tempJh_dB = &
                     cross_product(Jh_IID(i,j,:), Xyz0_D-Xyz_IID(i,j,:)) &
                     / (sqrt( sum( (Xyz_IID(i,j,:)-Xyz0_D)**2 )) )**3

                tempJp_dB = &
                     cross_product(Jp_IID(i,j,:), Xyz0_D-Xyz_IID(i,j,:)) &
                     / (sqrt( sum( (Xyz_IID(i,j,:)-Xyz0_D)**2 )) )**3

                dv = cMu/(4*cPi) * Radius**2 * dTheta(i)*dPsi(j)*SinTheta(i,j)

                MagJh_D = MagJh_D + tempJh_dB * dv
                MagJp_D = MagJp_D + tempJp_dB * dv                

             end do
          end do

       else
          ! northern hemisphere
          if(iProc /= 0)CYCLE
          do i = 1, nTheta
             do j = 1, nPsi

                tempJh_dB = &
                     cross_product(Jh_IID(i,j,:), Xyz0_D-Xyz_IID(i,j,:)) &
                     / (sqrt(sum((Xyz_IID(i,j,:)-Xyz0_D)**2)))**3

                tempJp_dB = &
                     cross_product(Jp_IID(i,j,:), Xyz0_D-Xyz_IID(i,j,:)) &
                     / (sqrt(sum((Xyz_IID(i,j,:)-Xyz0_D)**2)))**3

                dv = cMu/(4*cPi) * Radius**2 * dTheta(i)*dPsi(j)*SinTheta(i,j)

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
    call timing_stop('iono_mag_db')
    call timing_stop('iono_mag_perturb')

  end subroutine iono_mag_perturb

  !=====================================================================
  subroutine get_iono_magperturb_now(PerturbJh_DI, PerturbJp_DI, Xyz_DI)
    ! For all virtual magnetometers, update magnetometer coordinates in SMG 
    ! coordinates. Then, calculate the perturbation from the Hall and Pederson 
    ! currents and return them to caller as PerturbJhOut_DI, PerturbJpOut_DI.
    ! Updated magnetometer coordinates are also returned as Xyz_DI.

    use CON_axes, ONLY: transform_matrix
    use ModMpi

    implicit none

    real, intent(out), dimension(3,nMagnetometer) :: PerturbJh_DI, PerturbJp_DI
    real, intent(out), dimension(3,nMagnetometer) :: Xyz_DI

    real, dimension(3,nMagnetometer):: MagVarSum_Jh_DI, MagVarSum_Jp_DI
    real, dimension(3):: Xyz_D
    real, dimension(3,3) :: MagtoSmg_DD

    integer :: iMag, iError
    !--------------------------------------------------------------------------
    ! Get current positions of magnetometers in SMG coordinates.
    do iMag = 1 , nMagnetometer
       ! Create rotation matrix.
       MagtoSmg_DD = transform_matrix(Time_Simulation, TypeCoordMag_I(iMag), &
            'SMG')

       ! (360,360) is for the station at the center of the planet
       if ( nint(PosMagnetometer_II(1,iMag)) == 360 .and. &
            nint(PosMagnetometer_II(2,iMag)) == 360) then
          Xyz_DI(:,iMag) = 0.0
       else
          call  sph_to_xyz(IONO_Radius, &
               (90-PosMagnetometer_II(1,iMag))*cDegToRad, &
               PosMagnetometer_II(2,iMag)*cDegToRad, Xyz_D)
          Xyz_DI(:,iMag) = matmul(MagtoSmg_DD, Xyz_D)
       end if

    end do

    ! calculate the magnetic perturbation caused by Hall and Perdersen currents
    call iono_mag_perturb(nMagnetometer, Xyz_DI, PerturbJh_DI, PerturbJp_DI)

     !\
     ! Collect the variables from all the PEs
     !/
    MagVarSum_Jh_DI = 0.0
    MagVarSum_Jp_DI = 0.0
    if(nProc > 1)then 
       call MPI_reduce(PerturbJh_DI, MagVarSum_Jh_DI, 3*nMagnetometer, &
            MPI_REAL, MPI_SUM, 0, iComm, iError)
       call MPI_reduce(PerturbJp_DI, MagVarSum_Jp_DI, 3*nMagnetometer, &
            MPI_REAL, MPI_SUM, 0, iComm, iError)
       if(iProc == 0) then
          PerturbJh_DI = MagVarSum_Jh_DI
          PerturbJp_DI = MagVarSum_Jp_DI
       end if
    end if
    
  end subroutine get_iono_magperturb_now
  !=====================================================================

end module ModIonoMagPerturb
