!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModIonoMagPerturb

  use ModCoordTransform, ONLY: &
       sph_to_xyz, xyz_to_sph, cross_product, rot_xyz_sph
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
    real, dimension(nTheta*2, nPsi)    :: Theta, Psi, ETh, EPs, &
         SigmaH, SigmaP, SinTheta, SinPhi, CosTheta, CosPhi

    ! Potential with ghost cells in Psi direction
    real:: Phi_G(1:2*nTheta, 0:nPsi+1)

    real, dimension(3)                 :: bIono_D,  Xyz0_D, MagJh_D, MagJp_D, &
         XyzIono_D, tempJh_dB, tempJp_dB, TempMagJh_D,TempMagJp_D

    real :: dTheta, dPsi
    real :: dv
    real :: XyzSph_DD(3,3)
    integer :: i, j, iMag

    integer, parameter:: iDebug=10, jDebug=10
    logical, parameter:: DoDebug = .false.

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'iono_mag_perturb'
    !======================================================================
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting with XyzSm(iMag=1)=', &
         Xyz0_DI(:,1)

    call timing_start(NameSub)

    !\
    ! calculate the magnetic perturbations at the location of (SMLat, SMLon)
    ! in the SM coordinates, by integrating over the Hall and Perdersen
    ! current systems.
    !/
    JhMagPerturb_DI = 0.0
    JpMagPerturb_DI = 0.0

    Phi_G(1:nTheta,1:nPsi)            = Iono_North_PHi
    Phi_G(nTheta+1:2*nTheta,1:nPsi)   = Iono_South_PHi

    ! periodic in Phi with 1 and nPsi being the same!
    Phi_G(:,0)               = Phi_G(:,nPsi-1)
    Phi_G(:,nPsi+1)          = Phi_G(:,2)

    Theta(1:nTheta,:)          = Iono_North_Theta
    Theta(nTheta+1:nTheta*2,:) = Iono_South_Theta

    Psi(1:nTheta,:)            = Iono_North_Psi
    Psi(nTheta+1:nTheta*2,:)   = Iono_South_Psi

    SigmaH(1:nTheta,:)         = Iono_North_SigmaH
    SigmaH(nTheta+1:nTheta*2,:)= Iono_South_SigmaH

    SigmaP(1:nTheta,:)         = Iono_North_SigmaP
    SigmaP(nTheta+1:nTheta*2,:)= Iono_South_SigmaP

    ! Uniform grid resolution in Psi
    dPsi = cTwoPi/(IONO_nPsi-1)

    ! Uniform grid resolution in Theta (we ignore the poles)
    dTheta = cHalfPi/(IONO_nTheta-1)

    SinTheta = sin(Theta)
    SinPhi   = sin(Psi)
    CosTheta = cos(Theta)
    CosPhi   = cos(Psi)

    ! Calculate electric field. Ignore the poles
    Eth(1,:) = 0.0; Eth(2*nTheta,:) = 0.0
    EPs(1,:) = 0.0; Eps(2*nTheta,:) = 0.0
    do j = 1, nPsi; do i = 2, 2*nTheta-1
       ETh(i,j) = -(Phi_G(i+1,j) - Phi_G(i-1,j))/(2*dTheta*Radius)
       EPs(i,j) = -(Phi_G(i,j+1) - Phi_G(i,j-1)) &
            /(2*dPsi*Radius*SinTheta(i,j)) 
    end do; end do

    ! convert to xyz coords
    eIono_IID(:,:,1) =  ETh*CosTheta*CosPhi - EPs*SinPhi
    eIono_IID(:,:,2) =  ETh*CosTheta*SinPhi + EPs*CosPhi
    eIono_IID(:,:,3) = -ETh*SinTheta

    do i = 2, nTheta*2-1
       do j = 1, nPsi
          call sph_to_xyz(Radius, Theta(i,j), Psi(i,j), XyzIono_D)
          Xyz_IID(i,j,:) = XyzIono_D
          ! get the magnetic field in SMG coords.
          call get_planet_field(Time_simulation, XyzIono_D, 'SMG', bIono_D)
          bIono_D = bIono_D/sqrt(sum(bIono_D**2))

          ! get the Hall and Perdersen currents in xyz coords
          Jh_IID(i,j,:) = cross_product(bIono_D, eIono_IID(i,j,:))*SigmaH(i,j)
          Jp_IID(i,j,:) = eIono_IID(i,j,:) * SigmaP(i,j)

          if(DoDebug .and. iProc==0.and.i==iDebug.and.j==jDebug)then
             write(*,*)NameSub,': iono_mag_perturb'
             write(*,*)NameSub,': i,j,Theta,Psi=', i,j,Theta(i,j),Psi(i,j)
             write(*,*)NameSub,': rad, dTheta, dPsi, SinTheta=', &
                  Radius, dTheta, dPsi, SinTheta(i,j)
             write(*,*)NameSub,': SigmaH,SigmaP=', SigmaH(i,j),SigmaP(i,j)
             write(*,*)NameSub,': IonoPotential=', Phi_G(i-1:i+1,j-1:j+1)
             write(*,*)NameSub,': ETh, EPs     =', ETh(i,j), Eps(i,j)
             write(*,*)NameSub,': SinTheta,CosTheta,SinPhi,CosPhi=', &
                  SinTheta(i,j), CosTheta(i,j), SinPhi(i,j), CosPhi(i,j)
             write(*,*)NameSub,': bUnit_D      =', bIono_D, eIono_IID(i,j,:)
             write(*,*)NameSub,': eIono_D      =', eIono_IID(i,j,:)
             write(*,*)NameSub,': Jhall        =', Jh_IID(i,j,:)
             write(*,*)NameSub,': JPedersen    =', Jp_IID(i,j,:)
          endif

       end do
    end do


    call timing_start('iono_mag_db')

    do iMag=1, nMag

       Xyz0_D = Xyz0_DI(:,iMag)

       MagJh_D = 0.0
       MagJp_D = 0.0
       ! Biot-Savart integral to calculate the magnetic perturbations

       ! southern hemisphere
       if (iProc == nProc-1)then
          do i = nTheta+1, 2*nTheta-1   ! SKIP i=2*nTheta south pole
             do j = 1, nPsi

                tempJh_dB = &
                     cross_product(Jh_IID(i,j,:), Xyz0_D-Xyz_IID(i,j,:)) &
                     / (sqrt( sum( (Xyz_IID(i,j,:)-Xyz0_D)**2 )) )**3

                tempJp_dB = &
                     cross_product(Jp_IID(i,j,:), Xyz0_D-Xyz_IID(i,j,:)) &
                     / (sqrt( sum( (Xyz_IID(i,j,:)-Xyz0_D)**2 )) )**3

                ! dArea*mu/4pi
                dv = cMu/(4*cPi)*Radius**2*dTheta*dPsi*SinTheta(i,j)

                MagJh_D = MagJh_D + tempJh_dB * dv
                MagJp_D = MagJp_D + tempJp_dB * dv                

             end do
          end do
       end if
       ! northern hemisphere
       if(iProc == 0)then
          do i = 2, nTheta   ! SKIP i=1 north pole
             do j = 1, nPsi

                tempJh_dB = &
                     cross_product(Jh_IID(i,j,:), Xyz0_D-Xyz_IID(i,j,:)) &
                     / (sqrt(sum((Xyz_IID(i,j,:)-Xyz0_D)**2)))**3

                tempJp_dB = &
                     cross_product(Jp_IID(i,j,:), Xyz0_D-Xyz_IID(i,j,:)) &
                     / (sqrt(sum((Xyz_IID(i,j,:)-Xyz0_D)**2)))**3

                ! dArea*mu/4pi
                dv = cMu/(4*cPi)*Radius**2*dTheta*dPsi*SinTheta(i,j)

                MagJh_D = MagJh_D + tempJh_dB * dv
                MagJp_D = MagJp_D + tempJp_dB * dv

                if(DoDebug.and.i==iDebug.and.j==jDebug.and.iMag==1)then
                   write(*,*)NameSub,': Time=', Time_simulation
                   write(*,*)NameSub,': XyzSm,XyzIono=', Xyz0_D, Xyz_IID(i,j,:)
                   write(*,*)NameSub,': Radius**2, dTheta, dPsi, SinTheta=', &
                        Radius**2, dTheta, dPsi, SinTheta(i,j)
                   write(*,*)NameSub,': dArea, r3    =', dv*4*cPi/cMu, &
                        (sqrt(sum((Xyz_IID(i,j,:)-Xyz0_D)**2)))**3
                   write(*,*)NameSub,': dBHall,sum =', tempJh_dB * dv, MagJh_D
                   write(*,*)NameSub,': dBPede,sum =', tempJp_dB * dv, MagJp_D
                end if
             end do
          end do
       end if

       if(DoTestMe .and. iMag==1) write(*,*) &
            NameSub,': dBHall(iMag=1), dBPede(iMag=1)=', &
            MagJh_D, MagJp_D

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
    call timing_stop(NameSub)

    if(DoTestMe) write(*,*) &
         NameSub,': dBHall(iMag=1), dBPede(iMag=1) in NED=', &
         JhMagPerturb_DI(:,1), JpMagPerturb_DI(:,1)

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
