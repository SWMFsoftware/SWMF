subroutine IE_solve

  use IE_ModMain
  use IE_ModIO, ONLY: write_prefix
  use ModProcIE
  use ModIonosphere

  implicit none
  character(len=*), parameter :: NameSub = 'IE_solve'

  real    :: Radius
  integer :: iBlock

  logical DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  SinThetaTilt = sin(ThetaTilt)
  CosThetaTilt = cos(ThetaTilt)

  if (DoTestMe) then
     write(*,*) "=> Computing Ionospheric Potential at Time :"
     write(*,*) "           Year   : ",Time_Array(1)
     write(*,*) "           Month  : ",Time_Array(2)
     write(*,*) "           Day    : ",Time_Array(3)
     write(*,*) "           Hour   : ",Time_Array(4)
     write(*,*) "           Minute : ",Time_Array(5)
     write(*,*) "           Second : ",Time_Array(6)
  endif

  Radius = IONO_Radius + IONO_Height

  do iBlock = 1, 2
     select case(iBlock)
     case(1) ! Northern hemisphere
        if(iProc /= 0) CYCLE
        if (IONO_NORTH_nMagBndPts < 8) CYCLE

        if(DoTest)write(*,*)NameSub,': sum(abs(IONO_NORTH_JR))=',&
             sum(abs(IONO_NORTH_JR))

        if (UseFakeRegion2) then
           call Create_Region2_Currents(1)
           IONO_NORTH_JR = IONO_NORTH_JR + IONO_NORTH_Fake_JR
           if(DoTest)write(*,*)NameSub,': after UseFakeRegion2=',&
                sum(abs(IONO_NORTH_JR))
        end if

        call FACs_to_fluxes(conductance_model, iBlock)

        ! This is tricky, JR is passed as PHI to the ionosphere_solver
        IONO_NORTH_PHI = IONO_NORTH_JR   - &
             Iono_North_Tgcm_Jr

        call ionosphere_conductance(IONO_NORTH_Sigma0,               &
             IONO_NORTH_SigmaH, IONO_NORTH_SigmaP,    &
             IONO_NORTH_SigmaThTh,                    &
             IONO_NORTH_SigmaThPs,                    &
             IONO_NORTH_SigmaPsPs,                    &
             IONO_NORTH_dSigmaThTh_dTheta,            &
             IONO_NORTH_dSigmaThPs_dTheta,            &
             IONO_NORTH_dSigmaPsPs_dTheta,            &
             IONO_NORTH_dSigmaThTh_dPsi,              &
             IONO_NORTH_dSigmaThPs_dPsi,              &
             IONO_NORTH_dSigmaPsPs_dPsi,              &
             IONO_NORTH_EFlux, IONO_NORTH_Ave_E,      &
             IONO_NORTH_Theta, IONO_NORTH_Psi,        &
             IONO_nTheta, IONO_nPsi,                  &
             dTheta_North, dPsi_North,                &
             conductance_model, f107_flux)

        call ionosphere_solver(IONO_NORTH_PHI,             &
             IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs,   &
             IONO_NORTH_SigmaPsPs,                         &
             IONO_NORTH_dSigmaThTh_dTheta,                 &
             IONO_NORTH_dSigmaThPs_dTheta, &
             IONO_NORTH_dSigmaPsPs_dTheta, &
             IONO_NORTH_dSigmaThTh_dPsi,  &
             IONO_NORTH_dSigmaThPs_dPsi, &
             IONO_NORTH_dSigmaPsPs_dPsi, &
             IONO_NORTH_Theta, IONO_NORTH_Psi, &
             Radius, IONO_nTheta, IONO_nPsi,               &
             dTheta_North, dPsi_North,                     &
             4, 0)


        if(DoTest)then
           call write_prefix; 
           write(*,*) "Northern Cross Polar Cap Potential=",cpcp_north," kV"
        end if

        ! Calculate Currents and Boundary Conditions for North

        call ionosphere_currents(IONO_NORTH_Jx,IONO_NORTH_Jy,IONO_NORTH_Jz,&
             IONO_NORTH_Ex,IONO_NORTH_Ey,IONO_NORTH_Ez, &
             IONO_NORTH_ETh,IONO_NORTH_EPs, &
             IONO_NORTH_Ux,IONO_NORTH_Uy,IONO_NORTH_Uz, &
             IONO_NORTH_PHI, &
             IONO_NORTH_SigmaThTh, IONO_NORTH_SigmaThPs, &
             IONO_NORTH_SigmaPsPs, &
             IONO_NORTH_X, IONO_NORTH_Y, IONO_NORTH_Z, &
             IONO_NORTH_Theta, IONO_NORTH_Psi, &
             Radius, IONO_nTheta, IONO_nPsi,              &
             dTheta_North, dPsi_North)

     case(2) ! Southern hemisphere

        if(iProc /= nProc-1) CYCLE

        if (IONO_SOUTH_nMagBndPts < 8) CYCLE

        if(DoTest)write(*,*)NameSub,': sum(abs(IONO_SOUTH_JR))=',&
             sum(abs(IONO_SOUTH_JR))

        if (UseFakeRegion2) then
           call Create_Region2_Currents(iBlock)
           IONO_SOUTH_JR = IONO_SOUTH_JR + IONO_SOUTH_Fake_JR
           if(DoTest)write(*,*)NameSub,': after UseFakeRegion2=',&
                sum(abs(IONO_SOUTH_JR))
        end if

        call FACs_to_fluxes(conductance_model, iBlock)

        ! This is tricky, JR is passed as PHI to the ionosphere_solver
        IONO_SOUTH_PHI = IONO_SOUTH_JR   - &
             Iono_South_Tgcm_Jr

        call ionosphere_conductance(IONO_SOUTH_Sigma0,               &
             IONO_SOUTH_SigmaH, &
             IONO_SOUTH_SigmaP, &
             IONO_SOUTH_SigmaThTh, &
             IONO_SOUTH_SigmaThPs, &
             IONO_SOUTH_SigmaPsPs, &
             IONO_SOUTH_dSigmaThTh_dTheta, &
             IONO_SOUTH_dSigmaThPs_dTheta, &
             IONO_SOUTH_dSigmaPsPs_dTheta, &
             IONO_SOUTH_dSigmaThTh_dPsi, &
             IONO_SOUTH_dSigmaThPs_dPsi, &
             IONO_SOUTH_dSigmaPsPs_dPsi, &
             IONO_SOUTH_EFlux, IONO_SOUTH_Ave_E,  &
             IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
             IONO_nTheta, IONO_nPsi,                  &
             dTheta_South, dPsi_South,                &
             conductance_model, f107_flux)

        call ionosphere_solver(IONO_SOUTH_PHI, &
             IONO_SOUTH_SigmaThTh, &
             IONO_SOUTH_SigmaThPs, &
             IONO_SOUTH_SigmaPsPs, &
             IONO_SOUTH_dSigmaThTh_dTheta, &
             IONO_SOUTH_dSigmaThPs_dTheta, &
             IONO_SOUTH_dSigmaPsPs_dTheta, &
             IONO_SOUTH_dSigmaThTh_dPsi, &
             IONO_SOUTH_dSigmaThPs_dPsi, &
             IONO_SOUTH_dSigmaPsPs_dPsi, &
             IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
             Radius, IONO_nTheta, IONO_nPsi,               &
             dTheta_South, dPsi_South,                     &
             4, 1)

        if(DoTest)then
           call write_prefix; 
           write(*,*) "Southern Cross Polar Cap Potential=",cpcp_south," kV"
        end if

        ! Calculate Currents and Boundary Conditions for South

        call ionosphere_currents(IONO_SOUTH_Jx,IONO_SOUTH_Jy,IONO_SOUTH_Jz,&
             IONO_SOUTH_Ex,IONO_SOUTH_Ey,IONO_SOUTH_Ez, &
             IONO_SOUTH_ETh,IONO_SOUTH_EPs, &
             IONO_SOUTH_Ux,IONO_SOUTH_Uy,IONO_SOUTH_Uz, &
             IONO_SOUTH_PHI, &
             IONO_SOUTH_SigmaThTh, IONO_SOUTH_SigmaThPs, &
             IONO_SOUTH_SigmaPsPs, &
             IONO_SOUTH_X, IONO_SOUTH_Y, IONO_SOUTH_Z, &
             IONO_SOUTH_Theta, IONO_SOUTH_Psi, &
             Radius, IONO_nTheta, IONO_nPsi,              &
             dTheta_South, dPsi_South)

     end select

  end do

end subroutine IE_solve
