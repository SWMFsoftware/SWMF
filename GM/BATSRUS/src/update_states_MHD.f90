!^CFG COPYRIGHT UM
subroutine update_states_MHD(iStage,iBLK)
  use ModProcMH
  use ModMain
  use ModAdvance
  use ModGeometry, ONLY : &
       iVolumeCounterBLK,iVolumeCounterI,&  !^CFG IF NOT CARTESIAN
       R_BLK,VolumeInverse_I,RMin_BLK,body_BLK
  use ModPhysics
  use ModNumConst
  implicit none

  integer, intent(in) :: iStage,iBLK
  integer :: iVolumeCounter,i,j,k, iVar

  real*8 :: fullBx, fullBy, fullBz, fullBB, rhoc2, UdotBc2, gA2_Boris,&
       FullBxOld,FullByOld,FullBzOld,Ux,Uy,Uz,UxOld,UyOld,UzOld,&
       Bx,By,Bz,BxOld,ByOld,BzOld,B0x,B0y,B0z,RhoUx,RhoUy,RhoUz,&
       MBorisMinusRhoUxOld, MBorisMinusRhoUyOld, MBorisMinusRhoUzOld,&
       Rho,RhoInv,ECorr
  real::cfl_factor
  real:: DtLocal,VInvLocal
  !^CFG IF POINTIMPLICIT BEGIN

  real :: sigma_ptImplicit, &
       TimeCell, ResidCell, &
       ResidualCell(1:NVar), LHSCell(1:Nvar,1:Nvar)
  real :: inv_rho, inv_rho2
  real :: LHSMAX, LHSTEMP, TOTALSUM, SCALING(NVar)

  common /timestepping/ cfl_factor, sigma_ptImplicit, &
       TimeCell, ResidCell, &
       ResidualCell, LHSCell
  common /localvars/ inv_rho, inv_rho2, &
       LHSMAX, LHSTEMP, TOTALSUM, SCALING

  !T3E! !dir$ cache_align /timestepping/
  !T3E! !dir$ cache_align /localvars/
  !^CFG END POINTIMPLICIT
  logical :: oktest, oktest_me
  !--------------------------------------------------------------------------
  if(iBLK==BLKtest .and. iProc==PROCtest)then
     call set_oktest('update_states_MHD',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  endif

  !\
  ! Set variables depending on stage number
  !/

  cfl_factor = iStage*(cfl/nStage)

  !\
  ! Update the new solution state and calc residuals for the mth stage.
  ! Note must copy state to old state only if m is 1.
  !/

  if(iStage==1) then
     StateOld_VCB(1:nVar,1:nI,1:nJ,1:nK,iBLK) = & 
          State_VGB(1:nVar,1:nI,1:nJ,1:nK,iBLK)
     E_o_BLK(1:nI,1:nJ,1:nK,iBLK) = E_BLK(1:nI,1:nJ,1:nK,iBLK)
  end if

  !Get Residual.
  iVolumeCounter=iBLK
  iVolumeCounter=iVolumeCounterBLK*(iVolumeCounter-iVolumeCounterI) !^CFG IF NOT CARTESIAN
  do k = 1,nK; do j = 1,nJ; do i = 1,nI
     iVolumeCounter=iVolumeCounter+iVolumeCounterI   !^CFG IF NOT CARTESIAN
     DtLocal=cfl_factor*time_BLK(i,j,k,iBLK)
     VInvLocal=VolumeInverse_I(iVolumeCounter)
     Source_VC(:,i,j,k) = &
          DtLocal* (Source_VC(:,i,j,k) + &
          VInvLocal * &
          ( Flux_VX(:,i,j,k) - Flux_VX(:,i+1,j,k) &
          + Flux_VY(:,i,j,k) - Flux_VY(:,i,j+1,k) &
          + Flux_VZ(:,i,j,k) - Flux_VZ(:,i,j,k+1) )) 
  end do; end do; end do
  if (UsePointImplicit) then   !^CFG IF POINTIMPLICIT BEGIN
     call update_point_implicit
  else                         !^CFG END POINTIMPLICIT
     call update_explicit
  end if                       !^CFG IF POINTIMPLICIT

contains

  subroutine update_explicit
    do k=1,nK; do j=1,nJ; do i=1,nI
       do iVar=1,nVar
          State_VGB(iVar,i,j,k,iBLK) = &
               StateOld_VCB(iVar,i,j,k,iBLK) + &
               Source_VC(iVar,i,j,k)
       end do
       ! Compute energy. Choose which to keep below in the where statement
       E_BLK(i,j,k,iBLK) = E_o_BLK(i,j,k,iBLK) + Source_VC(Energy_,i,j,k)
    end do; end do; end do

    if((nStage==1.and..not.time_accurate).or.(nStage>1.and.iStage==1))then
       do k=1,nK; do j=1,nJ; do i=1,nI
          E_BLK(i,j,k,iBLK) = E_BLK(i,j,k,iBLK) + cHalf*( &
               Source_VC(Bx_,i,j,k)**2 + &
               Source_VC(By_,i,j,k)**2 + &
               Source_VC(Bz_,i,j,k)**2)
       end do; end do; end do
    end if

    if(boris_correction) then                 !^CFG IF BORISCORR BEGIN

       do k=1,nK; do j=1,nJ; do i=1,nI
          B0x= B0xCell_BLK(i,j,k,iBLK)
          B0y= B0yCell_BLK(i,j,k,iBLK)
          B0z= B0zCell_BLK(i,j,k,iBLK)

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

          E_BLK(i,j,k,iBLK) = E_BLK(i,j,k,iBLK)      &
               + cHalf*inv_c2LIGHT*ECorr

       end do; end do; end do
    endif                                    !^CFG END BORISCORR
    if(UseBorisSimple) then                  !^CFG IF SIMPLEBORIS BEGIN
       ! Convert simple Boris variables back to MHD variables

       do k=1,nK; do j=1,nJ; do i=1,nI
          fullBx = B0xCell_BLK(i,j,k,iBLK) + StateOld_VCB(Bx_,i,j,k,iBLK)
          fullBy = B0yCell_BLK(i,j,k,iBLK) + StateOld_VCB(By_,i,j,k,iBLK)
          fullBz = B0zCell_BLK(i,j,k,iBLK) + StateOld_VCB(Bz_,i,j,k,iBLK)
          fullBB = fullBx**2 + fullBy**2 + fullBz**2
          rhoc2  = StateOld_VCB(rho_,i,j,k,iBLK)*c2LIGHT
          gA2_Boris=fullBB/rhoc2

          ! rhoU_BorisSimple = rhoU*(1+BB/(rho*c2))
          State_VGB(rhoUx_:rhoU_+nDim,i,j,k,iBLK) = &
               State_VGB(rhoUx_:rhoU_+nDim,i,j,k,iBLK)+&
               StateOld_VCB(rhoUx_:rhoU_+nDim,i,j,k,iBLK)*ga2_Boris


          fullBx = B0xCell_BLK(i,j,k,iBLK) + State_VGB(Bx_,i,j,k,iBLK)
          fullBy = B0yCell_BLK(i,j,k,iBLK) + State_VGB(By_,i,j,k,iBLK)
          fullBz = B0zCell_BLK(i,j,k,iBLK) + State_VGB(Bz_,i,j,k,iBLK)
          fullBB = fullBx**2 + fullBy**2 + fullBz**2
          rhoc2  = State_VGB(rho_,i,j,k,iBLK)*c2LIGHT
          gA2_Boris=rhoc2/(fullBB+rhoc2)

          ! rhoU = 1/[1+BB/(rho c^2)]* rhoU_BorisSimple
          State_VGB(rhoUx_:rhoU_+nDim,i,j,k,iBLK) = gA2_Boris * &
               State_VGB(rhoUx_:rhoU_+nDim,i,j,k,iBLK)

       end do; end do; end do
    end if                                   !^CFG END SIMPLEBORIS

    if(UseNonConservative) then
       if(nConservCrit > 0)then
          where(IsConserv_CB(:,:,:,iBLK))
             State_VGB(P_,1:nI,1:nJ,1:nK,iBLK) = &
                  gm1*(E_BLK(1:nI,1:nJ,1:nK,iBLK) &
                  - cHalf*((State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBLK)**2 +&
                  State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBLK)**2 +&
                  State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBLK)**2) &
                  /State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)  &
                  +   State_VGB(Bx_,1:nI,1:nJ,1:nK,iBLK)**2 + &
                  State_VGB(By_,1:nI,1:nJ,1:nK,iBLK)**2 + &
                  State_VGB(Bz_,1:nI,1:nJ,1:nK,iBLK)**2) )
          elsewhere
             E_BLK(1:nI,1:nJ,1:nK,iBLK) = &
                  inv_gm1*State_VGB(P_,1:nI,1:nJ,1:nK,iBLK) &
                  + cHalf*((State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBLK)**2 +&
                  State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBLK)**2 +&
                  State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBLK)**2) &
                  /State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)  &
                  +   State_VGB(Bx_,1:nI,1:nJ,1:nK,iBLK)**2 + &
                  State_VGB(By_,1:nI,1:nJ,1:nK,iBLK)**2 + &
                  State_VGB(Bz_,1:nI,1:nJ,1:nK,iBLK)**2)
          end where
       else
          ! All cells are non-conservative
          E_BLK(1:nI,1:nJ,1:nK,iBLK) = &
               inv_gm1*State_VGB(P_,1:nI,1:nJ,1:nK,iBLK) &
               + cHalf*((State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBLK)**2 +&
               State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBLK)**2 +&
               State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBLK)**2) &
               /State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)  &
               +   State_VGB(Bx_,1:nI,1:nJ,1:nK,iBLK)**2 + &
               State_VGB(By_,1:nI,1:nJ,1:nK,iBLK)**2 + &
               State_VGB(Bz_,1:nI,1:nJ,1:nK,iBLK)**2)
       end if
    else
       ! All cells are conservative
       State_VGB(P_,1:nI,1:nJ,1:nK,iBLK) = gm1*(&
            E_BLK(1:nI,1:nJ,1:nK,iBLK) &
            - cHalf*((State_VGB(rhoUx_,1:nI,1:nJ,1:nK,iBLK)**2 +&
            State_VGB(rhoUy_,1:nI,1:nJ,1:nK,iBLK)**2 +&
            State_VGB(rhoUz_,1:nI,1:nJ,1:nK,iBLK)**2) &
            /State_VGB(rho_,1:nI,1:nJ,1:nK,iBLK)  &
            +   State_VGB(Bx_,1:nI,1:nJ,1:nK,iBLK)**2 + &
            State_VGB(By_,1:nI,1:nJ,1:nK,iBLK)**2 + &
            State_VGB(Bz_,1:nI,1:nJ,1:nK,iBLK)**2) )
    end if
  end subroutine update_explicit

  !^CFG IF POINTIMPLICIT BEGIN
  !================================================ 
  subroutine update_point_implicit
    sigma_ptImplicit = cOne
    do k=1,nK; do j=1,nJ; do i=1,nI 
       TimeCell = cfl_factor*time_BLK(i,j,k,iBLK)*sigma_ptImplicit
       inv_rho = 1.00/State_VGB(rho_,i,j,k,iBLK)
       inv_rho2 = inv_rho*inv_rho
       LHSCell = 0.00

       ResidualCell(1:7) = Source_VC(1:7,i,j,k)
       ResidualCell(8)   = Source_VC(Energy_,i,j,k)

       do iVar=1,7
          LHSCell(iVar,iVar) = 1.00
       end do

       LHSCell(2,1) = -TimeCell*fbody_x_BLK(i,j,k,iBLK)

       if (problem_type == problem_heliosphere) &
            LHSCell(2,3) = -2.00*TimeCell*OMEGAbody
       LHSCell(2,5) = TimeCell*DivB1_GB(i,j,k,iBLK)

       LHSCell(3,1) = -TimeCell*fbody_y_BLK(i,j,k,iBLK)
       if (problem_type == problem_heliosphere) &
            LHSCell(3,2) = 2.00*TimeCell*OMEGAbody
       LHSCell(3,6) = sigma_ptImplicit*TimeCell*DivB1_GB(i,j,k,iBLK)

       LHSCell(4,1) = -TimeCell*fbody_z_BLK(i,j,k,iBLK)
       LHSCell(4,7) = TimeCell*DivB1_GB(i,j,k,iBLK)

       LHSCell(5,1) = -TimeCell*DivB1_GB(i,j,k,iBLK)* &
            State_VGB(rhoUx_,i,j,k,iBLK)*inv_rho2
       LHSCell(5,2) = TimeCell*DivB1_GB(i,j,k,iBLK)*inv_rho


       LHSCell(6,1) = -TimeCell*DivB1_GB(i,j,k,iBLK)* &
            State_VGB(rhoUy_,i,j,k,iBLK)*inv_rho2
       LHSCell(6,3) = TimeCell*DivB1_GB(i,j,k,iBLK)*inv_rho



       LHSCell(7,1) = -TimeCell*DivB1_GB(i,j,k,iBLK)* &
            State_VGB(rhoUz_,i,j,k,iBLK)*inv_rho2
       LHSCell(7,4) = TimeCell*DivB1_GB(i,j,k,iBLK)*inv_rho

       LHSCell(8,1) = -TimeCell*( &
            DivB1_GB(i,j,k,iBLK)* &
            (State_VGB(Bx_,i,j,k,iBLK)* &
            State_VGB(rhoUx_,i,j,k,iBLK) + &
            State_VGB(By_,i,j,k,iBLK)* &
            State_VGB(rhoUy_,i,j,k,iBLK) + &
            State_VGB(Bz_,i,j,k,iBLK)* &
            State_VGB(rhoUz_,i,j,k,iBLK))*inv_rho2+ &
            (Theat0(i,j,k)-0.50*g*gm1* &
            (State_VGB(rhoUx_,i,j,k,iBLK)* &
            State_VGB(rhoUx_,i,j,k,iBLK) + &
            State_VGB(rhoUy_,i,j,k,iBLK)* &
            State_VGB(rhoUy_,i,j,k,iBLK) + &
            State_VGB(rhoUz_,i,j,k,iBLK)* &
            State_VGB(rhoUz_,i,j,k,iBLK))* &
            inv_rho2)*qheat_BLK(i,j,k,iBLK))
       LHSCell(8,2) = TimeCell*( &
            DivB1_GB(i,j,k,iBLK)*State_VGB(Bx_,i,j,k,iBLK)*inv_rho - &
            fbody_x_BLK(i,j,k,iBLK) - &
            g*gm1*qheat_BLK(i,j,k,iBLK)* &
            State_VGB(rhoUx_,i,j,k,iBLK)*inv_rho)
       LHSCell(8,3) = TimeCell*( &
            DivB1_GB(i,j,k,iBLK)*State_VGB(By_,i,j,k,iBLK)*inv_rho - &
            fbody_y_BLK(i,j,k,iBLK) - &
            g*gm1*qheat_BLK(i,j,k,iBLK)* &
            State_VGB(rhoUy_,i,j,k,iBLK)*inv_rho)
       LHSCell(8,4) = TimeCell*( &
            DivB1_GB(i,j,k,iBLK)*State_VGB(Bz_,i,j,k,iBLK)*inv_rho - &
            fbody_z_BLK(i,j,k,iBLK) - &
            g*gm1*qheat_BLK(i,j,k,iBLK)* &
            State_VGB(rhoUz_,i,j,k,iBLK)*inv_rho)
       LHSCell(8,5) = TimeCell*( &
            DivB1_GB(i,j,k,iBLK)*State_VGB(rhoUx_,i,j,k,iBLK)*inv_rho - &
            g*gm1*qheat_BLK(i,j,k,iBLK)* &
            State_VGB(Bx_,i,j,k,iBLK))
       LHSCell(8,6) = TimeCell*( &
            DivB1_GB(i,j,k,iBLK)*State_VGB(rhoUy_,i,j,k,iBLK)*inv_rho - &
            g*gm1*qheat_BLK(i,j,k,iBLK)* &
            State_VGB(By_,i,j,k,iBLK))
       LHSCell(8,7) = TimeCell*( &
            DivB1_GB(i,j,k,iBLK)*State_VGB(rhoUz_,i,j,k,iBLK)*inv_rho - &
            g*gm1*qheat_BLK(i,j,k,iBLK)* &
            State_VGB(Bz_,i,j,k,iBLK))
       LHSCell(8,8) = 1.00+TimeCell*( &
            g*gm1*qheat_BLK(i,j,k,iBLK))

       call linear_equation_solver

       State_VGB(rho_:Bz_,i,j,k,iBLK) =  &
            StateOld_VCB(rho_:Bz_,i,j,k,iBLK) &
            + ResidualCell(rho_:Bz_)

       E_BLK(i,j,k,iBLK) =   E_o_BLK(i,j,k,iBLK) &
            + ResidualCell(8)

       State_VGB(P_,i,j,k,iBLK) = &
            gm1*(E_BLK(i,j,k,iBLK) &
            - 0.5*((State_VGB(rhoUx_,i,j,k,iBLK)*&
            State_VGB(rhoUx_,i,j,k,iBLK) + &
            State_VGB(rhoUy_,i,j,k,iBLK)*&
            State_VGB(rhoUy_,i,j,k,iBLK) + &
            State_VGB(rhoUz_,i,j,k,iBLK)*&
            State_VGB(rhoUz_,i,j,k,iBLK)) &
            /State_VGB(rho_,i,j,k,iBLK))  &
            - 0.5*(State_VGB(Bx_,i,j,k,iBLK)*&
            State_VGB(Bx_,i,j,k,iBLK) + &
            State_VGB(By_,i,j,k,iBLK)*&
            State_VGB(By_,i,j,k,iBLK) + &
            State_VGB(Bz_,i,j,k,iBLK)*&
            State_VGB(Bz_,i,j,k,iBLK)))

    end do; end do; end do

  end subroutine update_point_implicit


  subroutine linear_equation_solver
    ! This routine solves the system of Nvar linear equations:
    ! 
    !               LHSCell*dUCell = ResidualCell.
    ! 
    ! The routine performs a lower-upper (LU) decomposition of the 
    ! square matrix LHSCell of rank Nvar and then uses forward and
    ! backward substitution to obtain the solution vector dUCell.
    ! Crout's method with partial implicit pivoting is used to perform
    ! the decompostion.
    integer :: IL, II, ILMAX, JL, KL, LL, INDX(1:Nvar)
    real, parameter :: TINY=1.0E-20


    !\
    ! Loop through each row to get implicit scaling
    ! information.
    !/
    DO IL=1,Nvar
       LHSMAX=0.00
       DO JL=1,Nvar
          IF (ABS(LHSCell(IL,JL)).GT.LHSMAX) LHSMAX=ABS(LHSCell(IL,JL))
       END DO
       SCALING(IL)=1.00/LHSMAX
    END DO

    !\
    ! Peform the LU decompostion using Crout's method.
    !/
    DO JL=1,Nvar
       DO IL=1,JL-1
          TOTALSUM=LHSCell(IL,JL)
          DO KL=1,IL-1
             TOTALSUM=TOTALSUM-LHSCell(IL,KL)*LHSCell(KL,JL)
          END DO
          LHSCell(IL,JL)=TOTALSUM
       END DO
       LHSMAX=0.00
       DO IL=JL,Nvar
          TOTALSUM=LHSCell(IL,JL)
          DO KL=1,JL-1
             TOTALSUM=TOTALSUM-LHSCell(IL,KL)*LHSCell(KL,JL)
          END DO
          LHSCell(IL,JL)=TOTALSUM
          LHSTEMP=SCALING(IL)*ABS(TOTALSUM)
          IF (LHSTEMP.GE.LHSMAX) THEN
             ILMAX=IL
             LHSMAX=LHSTEMP
          END IF
       END DO
       IF (JL.NE.ILMAX) THEN
          DO KL=1,Nvar
             LHSTEMP=LHSCell(ILMAX,KL)
             LHSCell(ILMAX,KL)=LHSCell(JL,KL)
             LHSCell(JL,KL)=LHSTEMP
          END DO
          SCALING(ILMAX)=SCALING(JL)
       END IF
       INDX(JL)=ILMAX
       IF (abs(LHSCell(JL,JL)).EQ.0.00) LHSCell(JL,JL)=TINY
       IF (JL.NE.Nvar) THEN
          LHSTEMP=1.00/LHSCell(JL,JL)
          DO IL=JL+1,Nvar
             LHSCell(IL,JL)=LHSCell(IL,JL)*LHSTEMP
          END DO
       END IF
    END DO

    !\
    ! Peform the forward and back substitution to obtain
    ! the solution vector.
    !/
    II=0
    DO IL=1,Nvar
       LL=INDX(IL)
       TOTALSUM=ResidualCell(LL)
       ResidualCell(LL)=ResidualCell(IL)
       IF (II.NE.0) THEN
          DO JL=II,IL-1
             TOTALSUM=TOTALSUM-LHSCell(IL,JL)*ResidualCell(JL)
          END DO
       ELSE IF (TOTALSUM.NE.0.00) THEN
          II=IL
       END IF
       ResidualCell(IL)=TOTALSUM
    END DO
    DO IL=Nvar,1,-1
       TOTALSUM=ResidualCell(IL)
       DO JL=IL+1,Nvar
          TOTALSUM=TOTALSUM-LHSCell(IL,JL)*ResidualCell(JL)
       END DO
       ResidualCell(IL)=TOTALSUM/LHSCell(IL,IL)
    END DO

  end subroutine linear_equation_solver

  !^CFG END POINTIMPLICIT

end subroutine update_states_mhd
