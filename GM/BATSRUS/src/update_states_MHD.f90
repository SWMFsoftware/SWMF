!^CFG COPYRIGHT UM
subroutine update_states_MHD(iStage,iBLK)
  use ModProcMH
  use ModMain
  use ModAdvance
  use ModGeometry, ONLY : R_BLK,vInv_CB,RMin_BLK,body_BLK
  use ModPhysics
  use ModNumConst
  use ModKind, ONLY: Real8_
  use ModPointImplicit, ONLY: UsePointImplicit, UsePointImplicit_B, &
       update_point_implicit
  use ModUser, ONLY: user_calc_sources, user_init_point_implicit
  use ModMultiIon, ONLY: multi_ion_source_impl, multi_ion_init_point_impl, &
       multi_ion_set_restrict, multi_ion_update, DoRestrictMultiIon
  use ModEnergy
  use ModWaves, ONLY: nWave, UseWavePressure, update_wave_group_advection
  implicit none

  integer, intent(in) :: iStage, iBLK


  integer :: i,j,k, iVar

  ! These variables have to be double precision for accurate Boris scheme
  real(Real8_) :: FullBx, FullBy, FullBz, fullBB, rhoc2, UdotBc2, gA2_Boris,&
       FullBxOld, FullByOld, FullBzOld, Ux, Uy, Uz, UxOld, UyOld, UzOld,&
       Bx, By, Bz, BxOld, ByOld, BzOld, B0x, B0y, B0z, RhoUx, RhoUy, RhoUz,&
       mBorisMinusRhoUxOld, mBorisMinusRhoUyOld, mBorisMinusRhoUzOld,&
       Rho, RhoInv, eCorr, p
  real:: DtFactor
  real:: DtLocal
  real:: B0_DC(3,nI,nJ,nK)
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub = 'update_states_mhd'
  !--------------------------------------------------------------------------
  if(iBLK==BLKtest .and. iProc==PROCtest)then
     call set_oktest(NameSub,DoTest,DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  endif

  !\
  ! Set variables depending on stage number
  !/
  DtFactor = iStage*(Cfl/nStage)

  !\
  ! Update the new solution state and calc residuals for the mth stage.
  ! Note must copy state to old state only if m is 1.
  !/

  if(iStage==1) then
     StateOld_VCB(1:nVar,1:nI,1:nJ,1:nK,iBLK) = & 
          State_VGB(1:nVar,1:nI,1:nJ,1:nK,iBLK)
     EnergyOld_CBI(:,:,:,iBLK,:) = Energy_GBI(1:nI,1:nJ,1:nK,iBLK,:)
  end if

  !Get Residual.
  do k = 1,nK; do j = 1,nJ; do i = 1,nI
     DtLocal=DtFactor*time_BLK(i,j,k,iBLK)
     Source_VC(:,i,j,k) = &
          DtLocal* (Source_VC(:,i,j,k) + &
          vInv_CB(i,j,k,iBLK) * &
          ( Flux_VX(:,i,j,k) - Flux_VX(:,i+1,j,k) &
          + Flux_VY(:,i,j,k) - Flux_VY(:,i,j+1,k) &
          + Flux_VZ(:,i,j,k) - Flux_VZ(:,i,j,k+1) )) 
  end do; end do; end do

  if(UseMultiIon .and. DoRestrictMultiIon)call multi_ion_set_restrict(iBlk)

  if(DoTestMe)write(*,*) NameSub, ' original testvar and energy=', &
       State_VGB(VarTest, iTest, jTest, kTest, iBlk), &
       Energy_GBI(iTest,jTest,kTest,iBlk,1)

  call update_explicit

  if(UseMultiIon .and. IsMhd)then
     call multi_ion_update(iBlk, IsFinal = .false.)

     if(DoTestMe)write(*,*) NameSub, ' after multiion update1=', &
          State_VGB(VarTest,iTest,jTest,kTest,iBlk), &
          Energy_GBI(iTest,jTest,kTest,iBlk,1)

  end if

  ! The point implicit update and other stuff below are only done in last stage
  if(iStage < nStage) RETURN

  ! Add point implicit user or multi-ion source terms
  if (UsePointImplicit .and. UsePointImplicit_B(iBLK))then
     if(UseMultiIon)then
        call update_point_implicit(iBLK, multi_ion_source_impl, &
             multi_ion_init_point_impl)
     elseif(UseUserSource) then
        call update_point_implicit(iBLK, user_calc_sources, &
             user_init_point_implicit)
     end if

     if(DoTestMe)write(*,*) NameSub, ' after point impl state=', &
          State_VGB(VarTest, iTest,jTest,kTest,iBlk), &
          Energy_GBI(iTest,jTest,kTest,iBlk,1)
  end if

  if(UseMultiIon .and. IsMhd)then
     call multi_ion_update(iBlk, IsFinal = .true.)
     if(DoTestMe)write(*,*) NameSub,' after multiion update2=', &
          State_VGB(VarTest,iTest,jTest,kTest,iBlk), &
          Energy_GBI(iTest,jTest,kTest,iBlk,1)
  end if

  if(UseHyperbolicDivb .and. HypDecay > 0) &
       State_VGB(Hyp_,1:nI,1:nJ,1:nK,iBlk) = &
       State_VGB(Hyp_,1:nI,1:nJ,1:nK,iBlk)*(1 - HypDecay)

  if(DoTestMe)write(*,*) NameSub, ' final state=', &
       State_VGB(VarTest,iTest,jTest,kTest,iBlk), &
       Energy_GBI(iTest,jTest,kTest,iBlk,1)

contains

  subroutine update_explicit

    do k=1,nK; do j=1,nJ; do i=1,nI
       do iVar=1,nVar
          State_VGB(iVar,i,j,k,iBLK) = &
               StateOld_VCB(iVar,i,j,k,iBLK) + &
               Source_VC(iVar,i,j,k)
       end do
       ! Compute energy. Choose which to keep below in the where statement
       Energy_GBI(i,j,k,iBLK,:) = EnergyOld_CBI(i,j,k,iBLK,:) + &
            Source_VC(Energy_:Energy_-1+nFluid,i,j,k)
    end do; end do; end do

    if(DoTestMe)write(*,*) NameSub, ' after flux/source=', &
         State_VGB(VarTest,iTest,jTest,kTest,iBlk), &
         Energy_GBI(iTest,jTest,kTest,iBlk,1)

    if(UseMultiSpecies)then
       ! Fix negative species densities
       State_VGB(SpeciesFirst_:SpeciesLast_,1:nI,1:nJ,1:nK,iBLK) = &
            max(0.0, State_VGB(SpeciesFirst_:SpeciesLast_,1:nI,1:nJ,1:nK,iBLK))

       if(DoReplaceDensity)then
          ! Add up species densities to total density
          do k=1,nK; do j=1,nJ; do i=1,nI
             State_VGB(Rho_,i,j,k,iBLK) = &
                  sum(State_VGB(SpeciesFirst_:SpeciesLast_,i,j,k,iBLK))
          end do; end do; end do
       end if

       if(DoTestMe)write(*,*) NameSub, ' after multispecies correct=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlk)
    
    end if

    if( IsMhd .and. &
         ((nStage==1.and..not.time_accurate).or.(nStage>1.and.iStage==1)))then
       do k=1,nK; do j=1,nJ; do i=1,nI
          Energy_GBI(i,j,k,iBLK,1) = Energy_GBI(i,j,k,iBLK,1) + cHalf*( &
               Source_VC(Bx_,i,j,k)**2 + &
               Source_VC(By_,i,j,k)**2 + &
               Source_VC(Bz_,i,j,k)**2)
       end do; end do; end do

       if(DoTestMe)write(*,*) NameSub, ' after energy dB correct=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlk), &
            Energy_GBI(iTest,jTest,kTest,iBlk,1)

    end if

    if(iStage == nStage.and. UseWavePressure.and.nWave>=2)&
         call update_wave_group_advection(iBLK) 

    if(boris_correction) then                 !^CFG IF BORISCORR BEGIN
       if(UseB0)then
          B0_DC=B0_DGB(:,1:nI,1:nJ,1:nK,iBLK)
       else
          B0_DC=0.00
       end if
       do k=1,nK; do j=1,nJ; do i=1,nI
          B0x= B0_DC(x_,i,j,k)
          B0y= B0_DC(y_,i,j,k)
          B0z= B0_DC(z_,i,j,k)

          BxOld= StateOld_VCB(Bx_,i,j,k,iBLK)
          ByOld= StateOld_VCB(By_,i,j,k,iBLK)
          BzOld= StateOld_VCB(Bz_,i,j,k,iBLK)
          fullBxOld = B0x + BxOld
          fullByOld = B0y + ByOld
          fullBzOld = B0z + BzOld

          Rho = StateOld_VCB(rho_,i,j,k,iBLK)
          rhoc2  = Rho*c2LIGHT

          RhoUx = StateOld_VCB(rhoUx_,i,j,k,iBLK)
          RhoUy = StateOld_VCB(rhoUy_,i,j,k,iBLK)
          RhoUz = StateOld_VCB(rhoUz_,i,j,k,iBLK)

          RhoInv=1/Rho

          UxOld=RhoUx*RhoInv
          UyOld=RhoUy*RhoInv
          UzOld=RhoUz*RhoInv

          UdotBc2= (RhoUx*fullBxOld + &
               RhoUy*fullByOld + &
               RhoUz*fullBzOld)/rhoc2

          gA2_Boris= (fullBxOld**2 + fullByOld**2 + fullBzOld**2)/rhoc2

          ! rhoU_Boris_old - rhoU_Old= (U B^2 - B U.B)/c^2
          !    

          MBorisMinusRhoUxOld = RhoUx*ga2_Boris - fullBxOld*UdotBc2
          MBorisMinusRhoUyOld = RhoUy*ga2_Boris - fullByOld*UdotBc2
          MBorisMinusRhoUzOld = RhoUz*ga2_Boris - fullBzOld*UdotBc2

          Bx= State_VGB(Bx_,i,j,k,iBLK)
          By= State_VGB(By_,i,j,k,iBLK)
          Bz= State_VGB(Bz_,i,j,k,iBLK)

          fullBx = B0x + Bx
          fullBy = B0y + By
          fullBz = B0z + Bz
          fullBB = fullBx**2 + fullBy**2 + fullBz**2

          Rho  = State_VGB(rho_,i,j,k,iBLK)
          rhoc2  = Rho*c2LIGHT
          RhoUx = State_VGB(rhoUx_,i,j,k,iBLK)
          RhoUy = State_VGB(rhoUy_,i,j,k,iBLK)
          RhoUz = State_VGB(rhoUz_,i,j,k,iBLK)

          UdotBc2= (RhoUx*fullBx + &
               MBorisMinusRhoUxOld*Source_VC(Bx_,i,j,k)+ &
               RhoUy*fullBy +      &
               MBorisMinusRhoUyOld*Source_VC(By_,i,j,k)+ &
               RhoUz*fullBz+       &
               MBorisMinusRhoUzOld*Source_VC(Bz_,i,j,k))/rhoc2
          gA2_Boris=rhoc2/(fullBB+rhoc2)

          ! rhoU = 1/[1+BB/(rho c^2)]* (rhoU_Boris 
          !                             + (rhoUBorisdotB/(rho c^2) * B)
          !  ((M_Boris_old-RhoU) /cdot FullB_old)=0! 

          RhoUx = gA2_Boris * &
               (RhoUx+MBorisMinusRhoUxOld+UdotBc2*fullBx)
          RhoUy = gA2_Boris * &
               (RhoUy+MBorisMinusRhoUyOld+UdotBc2*fullBy)
          RhoUz = gA2_Boris * &
               (RhoUz+MBorisMinusRhoUzOld+UdotBc2*fullBz)


          State_VGB(rhoUx_,i,j,k,iBLK) = RhoUx
          State_VGB(rhoUy_,i,j,k,iBLK) = RhoUy
          State_VGB(rhoUz_,i,j,k,iBLK) = RhoUz

          RhoInv=1/Rho

          Ux=RhoUx*RhoInv
          Uy=RhoUy*RhoInv
          Uz=RhoUz*RhoInv

          ! E = E_Boris - (UxB)^2/(2 c^2)
          ECorr = (UyOld*BzOld-Uy*Bz+(UyOld-Uy)*B0z  &
               -(UzOld*ByOld-Uz*By+(UzOld-Uz)*B0y))* &
               (UyOld*fullBzOld+Uy*fullBz            &
               -UzOld*fullByOld-Uz*fullBy)           &
               +(UxOld*BzOld-Ux*Bz+(UxOld-Ux)*B0z    &
               -(UzOld*BxOld-Uz*Bx+(UzOld-Uz)*B0x))* &
               (UxOld*FullBzOld+Ux*FullBz            &
               -UzOld*FullBxOld-Uz*FullBx)           &
               +(UxOld*ByOld-Ux*By+(UxOld-Ux)*B0y    &
               -(UyOld*BxOld-Uy*Bx+(UyOld-Uy)*B0x))* &
               (UxOld*fullByOld+Ux*FullBy            &
               -UyOld*fullBxOld-Uy*FullBx)

          Energy_GBI(i,j,k,iBLK,1) = Energy_GBI(i,j,k,iBLK,1)  &
               + cHalf*inv_c2LIGHT*ECorr

       end do; end do; end do

       if(DoTestMe)write(*,*) NameSub, ' after Boris update=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlk), &
            Energy_GBI(iTest,jTest,kTest,iBlk,1)

    endif                                    !^CFG END BORISCORR
    if(UseBorisSimple) then                  !^CFG IF SIMPLEBORIS BEGIN
       if(UseB0)then
          B0_DC=B0_DGB(:,1:nI,1:nJ,1:nK,iBLK)
       else
          B0_DC=0.00
       end if
       ! Convert simple Boris variables back to MHD variables

       do k=1,nK; do j=1,nJ; do i=1,nI
          fullBx = B0_DC(x_,i,j,k) + StateOld_VCB(Bx_,i,j,k,iBLK)
          fullBy = B0_DC(y_,i,j,k) + StateOld_VCB(By_,i,j,k,iBLK)
          fullBz = B0_DC(z_,i,j,k) + StateOld_VCB(Bz_,i,j,k,iBLK)
          fullBB = fullBx**2 + fullBy**2 + fullBz**2
          rhoc2  = StateOld_VCB(rho_,i,j,k,iBLK)*c2LIGHT
          gA2_Boris=fullBB/rhoc2

          ! rhoU_BorisSimple = rhoU*(1+BB/(rho*c2))
          State_VGB(rhoUx_:rhoU_+nDim,i,j,k,iBLK) = &
               State_VGB(rhoUx_:rhoU_+nDim,i,j,k,iBLK)+&
               StateOld_VCB(rhoUx_:rhoU_+nDim,i,j,k,iBLK)*ga2_Boris


          fullBx = B0_DGB(x_,i,j,k,iBLK) + State_VGB(Bx_,i,j,k,iBLK)
          fullBy = B0_DGB(y_,i,j,k,iBLK) + State_VGB(By_,i,j,k,iBLK)
          fullBz = B0_DGB(z_,i,j,k,iBLK) + State_VGB(Bz_,i,j,k,iBLK)
          fullBB = fullBx**2 + fullBy**2 + fullBz**2
          rhoc2  = State_VGB(rho_,i,j,k,iBLK)*c2LIGHT
          gA2_Boris=rhoc2/(fullBB+rhoc2)

          ! rhoU = 1/[1+BB/(rho c^2)]* rhoU_BorisSimple
          State_VGB(rhoUx_:rhoU_+nDim,i,j,k,iBLK) = gA2_Boris * &
               State_VGB(rhoUx_:rhoU_+nDim,i,j,k,iBLK)

       end do; end do; end do

       if(DoTestMe)write(*,*) NameSub, ' after BorisSimple update=', &
            State_VGB(VarTest,iTest,jTest,kTest,iBlk), &
            Energy_GBI(iTest,jTest,kTest,iBlk,1)

    end if                                   !^CFG END SIMPLEBORIS

    ! Update energy or pressure based on UseConservative and IsConserv_CB
    call calc_energy_or_pressure(iBlk)

    if(DoTestMe)write(*,*) NameSub, ' after pressure/energy update=', &
         State_VGB(VarTest,iTest,jTest,kTest,iBlk), &
         Energy_GBI(iTest,jTest,kTest,iBlk,1)

  end subroutine update_explicit

end subroutine update_states_mhd
