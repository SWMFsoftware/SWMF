!^CFG FILE ROEFLUX
!^CFG COPYRIGHT UM

! Use maximum eigenvalue for all (like Rusanov) at inner BC

subroutine option_roeflux(TrueOption,NameOption)

  logical, intent(out) :: TrueOption
  character (len=40), intent(out) :: NameOption

  TrueOption  = .true.
  NameOption  = 'ROE FLUX De Zeeuw, Sokolov 1.1'

end subroutine option_roeflux

subroutine calc_flux_roe(DoResChangeOnly)
  use ModProcMH
  use ModMain
  use ModGeometry, ONLY : fAx_BLK,fAy_BLK,fAz_BLK,true_cell
  use ModPhysics, ONLY : g,gm1,inv_gm1
  use ModNumConst
  use ModAdvance
  use ModParallel, ONLY : neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  implicit none

  logical, intent (in) :: DoResChangeOnly

  integer, parameter :: nFlux=9, nWave=8, MaxStrip=10, E_=8
  integer, parameter :: RhoUn_=RhoUx_, RhoUt1_=RhoUy_, RhoUt2_=RhoUz_, &
                        B1n_=Bx_, B1t1_=By_, B1t2_=Bz_

  integer, parameter :: EntropyW_=1, AlfvenRW_=2, AlfvenLW_=3, &
                        SlowRW_=4, FastRW_=5, SlowLW_=6, FastLW_=7, DivBW_=8 

  real :: cSqrt2, cSqrt2Inv

!!!  real :: zz(68)  ! for equivalence statement
  integer :: i,j,k,iFlux,iWave,iDir,iStart,iStrip,nStrip

  !  integer (KIND=8) :: len

  ! Left face
  real, dimension(MaxStrip) :: RhoL_I,UnL_I,Ut1L_I,Ut2L_I
  real, dimension(MaxStrip) :: BnL_I,Bt1L_I,Bt2L_I,BbL_I
  real, dimension(MaxStrip) :: B1nL_I,B1t1L_I,B1t2L_I,Bb1L_I
  real, dimension(MaxStrip) :: pL_I,eL_I,aL_I,CsL_I,CfL_I

  ! Right face
  real, dimension(MaxStrip) :: RhoR_I,UnR_I,Ut1R_I,Ut2R_I
  real, dimension(MaxStrip) :: BnR_I,Bt1R_I,Bt2R_I,BbR_I
  real, dimension(MaxStrip) :: B1nR_I,B1t1R_I,B1t2R_I,Bb1R_I
  real, dimension(MaxStrip) :: pR_I,eR_I,aR_I,CsR_I,CfR_I

  ! Average (hat)
  real, dimension(MaxStrip) :: RhoH_I,UnH_I,Ut1H_I,Ut2H_I
  real, dimension(MaxStrip) :: BnH_I,Bt1H_I,Bt2H_I,BbH_I
  real, dimension(MaxStrip) :: B1nH_I,B1t1H_I,B1t2H_I,Bb1H_I
  real, dimension(MaxStrip) :: pH_I,eH_I,UuH_I
  real, dimension(MaxStrip) :: aH_I,CsH_I,CfH_I

  ! More face variables
  real, dimension(MaxStrip) :: B0n_I,B0t1_I,B0t2_I

  real, dimension(MaxStrip) :: BetaY_I, BetaZ_I, AlphaS_I, AlphaF_I

  real, dimension(MaxStrip) :: RhoInvL_I,RhoInvR_I,RhoInvH_I
  real, dimension(MaxStrip) :: RhoSqrtH_I,    RhoSqrtL_I,    RhoSqrtR_I, &
                               RhoInvSqrtH_I, RhoInvSqrtL_I, RhoInvSqrtR_I
  real, dimension(MaxStrip) :: Tmp1_I !!! v_tmp2,v_tmp3

  ! Jumps in conservative variables
  real, dimension(nVar) :: Delta_V                      

  ! Eigenvalues and jumps in characteristic variable
  real, dimension(nWave) :: Eigenvalue_V,DeltaWave_V 

  ! Eigenvectors
  real, dimension(nVar ,nWave):: EigenvectorL_VV       ! Left Eigenvectors
  real, dimension(nWave,nFlux):: EigenvectorR_VV       ! Right Eigenvector

  ! Fluxes
  real, dimension(nFlux)      :: FluxL_V, FluxR_V, FluxFull_V  

  ! Logical to use Rusanov flux at inner boundary
  logical, dimension(MaxStrip) :: UseFluxRusanov_I

  ! Misc. scalar variables
  real :: SignBnH,Tmp1,Tmp2,Tmp3,Gamma1A2Inv,DtInvVolume,AreaFace

  common/Strip/ &
       UnL_I,Ut1L_I,Ut2L_I,BnL_I,Bt1L_I,Bt2L_I,             &
       pL_I,eL_I,BbL_I,                                     &
       aL_I,aR_I,aH_I,CsL_I,CsR_I,CsH_I,CfL_I,CfR_I,CfH_I,  &
       UnR_I,Ut1R_I,Ut2R_I,BnR_I,Bt1R_I,Bt2R_I,             &
       pR_I,eR_I,BbR_I,                                     &
       RhoH_I,RhoL_I,RhoR_I,UnH_I,Ut1H_I,Ut2H_I,BnH_I,      &
       Bt1H_I,Bt2H_I,pH_I,eH_I,UuH_I,BbH_I,                 &
       RhoInvL_I,RhoInvR_I,RhoInvH_I,                       &
       BetaY_I,BetaZ_I,AlphaF_I,AlphaS_I,                   &
       RhoSqrtH_I,RhoSqrtL_I,RhoSqrtR_I,RhoInvSqrtH_I,      &
       RhoInvSqrtL_I,RhoInvSqrtR_I,Tmp1_I

  common/Strip2/ &
       B1nL_I,B1t1L_I,B1t2L_I,Bb1L_I,  &
       B1nR_I,B1t1R_I,B1t2R_I,Bb1R_I,  &
       B1nH_I,B1t1H_I,B1t2H_I,Bb1H_I,  &
       B0n_I,B0t1_I,B0t2_I

  common/Scalar/ &
       SignBnH,Tmp1,Tmp2,Tmp3,Gamma1A2Inv,DtInvVolume      

  common/Eig1/ Eigenvalue_V,DeltaWave_V                
  common/Eig2/ EigenvectorL_VV,EigenvectorR_VV 
  common/Flux/ FluxL_V,FluxR_V,FluxFull_V         

!!!  equivalence(zz(1),rho_lf)

  logical :: DoTest, DoTestMe, DoTestRow

  !---------------------------------------------------------------------------

  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('calc_facefluxes',DoTest,DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if

  !T3E! !dir$ cache_align /vface/ 
  !T3E! !dir$ cache_align /vface2/ 
  !T3E! !dir$ cache_align /vface3/ 
  !T3E! !dir$ cache_align /face1/ 
  !T3E! !dir$ cache_align /face2/
  !T3E! !dir$ cache_align /face3/ 
  !T3E! !dir$ cache_align /Eig1/ 
  !T3E! !dir$ cache_align /Eig2/  
  !T3E! !dir$ cache_align /flux/

  cSqrt2 = sqrt(2.)
  cSqrt2Inv = 1./cSqrt2

  !\
  ! x-face fluxes --- x-face fluxes --- x-face fluxes
  !/
  iDir = 1
  AreaFace = fAx_BLK(globalBLK)
  if (.not.DoResChangeOnly) then
     do k=kMinFaceX,kMaxFaceX
        do j=jMinFaceX,jMaxFaceX
           do iStart=1,nIFace,MaxStrip
              nStrip = min(MaxStrip,nIFace-(iStart-1))
              call get_fluxes
           end do
        end do
     end do
  else if (neiLeast(globalBLK)==+1) then
     iStart=1
     nStrip=1
     do k=1,nK
        do j=1,nJ
           call get_fluxes
        end do
     end do
  else if ( neiLwest(globalBLK)==+1) then
     iStart = nIFace
     nStrip = 1
     do k=1,nK
        do j=1,nJ
           call get_fluxes
        end do
     end do
  end if

  ! y-face fluxes --- y-face fluxes --- y-face fluxes
  iDir = 2
  AreaFace = fAy_BLK(globalBLK)
  if (.not.DoResChangeOnly) then
     do k=kMinFaceY,kMaxFaceY
        do j=1,nJFace
           do iStart=iMinFaceY,iMaxFaceY,MaxStrip
              nStrip = min(MaxStrip,iMaxFaceY-(iStart-1))
              call get_fluxes
           end do
        end do
     end do
  else if(neiLsouth(globalBLK)==+1)then
     j = 1
     do k=1,nK
        do iStart=1,nI,MaxStrip
           nStrip = min(MaxStrip,nI-(iStart-1))
           call get_fluxes
        end do
     end do
  else if (neiLnorth(globalBLK)==+1) then
     j = nJFace 
     do k=1,nK
        do iStart=1,nI,MaxStrip
           nStrip = min(MaxStrip,nI-(iStart-1))
           call get_fluxes
        end do
     end do
  end if

  ! z-face fluxes --- z-face fluxes --- z-face fluxes 
  iDir = 3
  AreaFace = fAz_BLK(globalBLK)
  if (.not.DoResChangeOnly) then
     do k=1,nKFace
        do j=jMinFaceZ,jMaxFaceZ
           do iStart=iMinFaceZ,iMaxFaceZ,MaxStrip
              nStrip = min(MaxStrip,iMaxFaceZ-(iStart-1))
              call get_fluxes
           end do
        end do
     end do
  else if (neiLbot(globalBLK)==+1) then
     k = 1
     do j=1,nJ
        do iStart=1,nI,MaxStrip
           nStrip = min(MaxStrip,nI-(iStart-1))
           call get_fluxes
        end do
     end do
  else if ( neiLtop(globalBLK)==+1) then
     k = nKFace            
     do j=1,nJ
        do iStart=1,nI,MaxStrip
           nStrip = min(MaxStrip,nI-(iStart-1))
           call get_fluxes
        end do
     end do
  end if

Contains

  !==========================================================================
  !==========================================================================
  !==========================================================================

  subroutine get_fluxes

    DoTestRow = DoTestMe .and. j==Jtest .and. k==Ktest

    !-------------------------
    select case (iDir)
    case (1)
       !\
       ! x face
       !/
       do i=1,nStrip
          iStrip = i+iStart-1
          !\
          ! B0 on the face
          !/
          B0n_I(i)  = B0xFace_x_BLK(iStrip,j,k,globalBLK)
          B0t1_I(i) = B0yFace_x_BLK(iStrip,j,k,globalBLK)
          B0t2_I(i) = B0zFace_x_BLK(iStrip,j,k,globalBLK)
          !\
          ! Left face
          !/
          RhoL_I(i)  =  LeftState_VX(rho_,iStrip,j,k)
          UnL_I(i)   =  LeftState_VX(Ux_ ,iStrip,j,k)
          Ut1L_I(i)  =  LeftState_VX(Uy_ ,iStrip,j,k)
          Ut2L_I(i)  =  LeftState_VX(Uz_ ,iStrip,j,k)
          B1nL_I(i)  =  LeftState_VX(Bx_ ,iStrip,j,k)
          B1t1L_I(i) =  LeftState_VX(By_ ,iStrip,j,k)
          B1t2L_I(i) =  LeftState_VX(Bz_ ,iStrip,j,k)
          pL_I(i)    =  LeftState_VX(P_  ,iStrip,j,k)
          !\
          ! Right face
          !/
          RhoR_I(i)  =  RightState_VX(rho_,iStrip,j,k)
          UnR_I(i)   =  RightState_VX(Ux_ ,iStrip,j,k)
          Ut1R_I(i)  =  RightState_VX(Uy_ ,iStrip,j,k)
          Ut2R_I(i)  =  RightState_VX(Uz_ ,iStrip,j,k)
          B1nR_I(i)  =  RightState_VX(Bx_ ,iStrip,j,k)
          B1t1R_I(i) =  RightState_VX(By_ ,iStrip,j,k)
          B1t2R_I(i) =  RightState_VX(Bz_ ,iStrip,j,k)
          pR_I(i)    =  RightState_VX(P_  ,iStrip,j,k)

          UseFluxRusanov_I(i)=true_cell(iStrip-1,j,k,globalBLK).neqv.&
               true_cell(iStrip  ,j,k,globalBLK)
       end do
    case (2)
       !\
       ! y face
       !/
       do i=1,nStrip
          iStrip = i+iStart-1
          !\
          ! B0 on the face
          !/
          B0n_I(i)  = B0yFace_y_BLK(iStrip,j,k,globalBLK)
          B0t1_I(i) = B0zFace_y_BLK(iStrip,j,k,globalBLK)
          B0t2_I(i) = B0xFace_y_BLK(iStrip,j,k,globalBLK)
          !\
          ! Left face
          !/
          RhoL_I(i)  =  LeftState_VY(rho_,iStrip,j,k)
          UnL_I(i)   =  LeftState_VY(Uy_ ,iStrip,j,k)
          Ut1L_I(i)  =  LeftState_VY(Uz_ ,iStrip,j,k)
          Ut2L_I(i)  =  LeftState_VY(Ux_ ,iStrip,j,k)
          B1nL_I(i)  =  LeftState_VY(By_ ,iStrip,j,k)
          B1t1L_I(i) =  LeftState_VY(Bz_ ,iStrip,j,k)
          B1t2L_I(i) =  LeftState_VY(Bx_ ,iStrip,j,k)
          pL_I(i)    =  LeftState_VY(P_  ,iStrip,j,k)
          !\
          ! Right face
          !/
          RhoR_I(i)  =  RightState_VY(rho_,iStrip,j,k)
          UnR_I(i)   =  RightState_VY(Uy_ ,iStrip,j,k)
          Ut1R_I(i)  =  RightState_VY(Uz_ ,iStrip,j,k)
          Ut2R_I(i)  =  RightState_VY(Ux_ ,iStrip,j,k)
          B1nR_I(i)  =  RightState_VY(By_ ,iStrip,j,k)
          B1t1R_I(i) =  RightState_VY(Bz_ ,iStrip,j,k)
          B1t2R_I(i) =  RightState_VY(Bx_ ,iStrip,j,k)
          pR_I(i)    =  RightState_VY(P_  ,iStrip,j,k)

          UseFluxRusanov_I(i)=true_cell(iStrip,j-1,k,globalBLK).neqv.&
               true_cell(iStrip,j  ,k,globalBLK)
       end do
    case (3)
       !\
       ! z face
       !/
       do i=1,nStrip
          iStrip = i+iStart-1
          !\
          ! B0 on the face
          !/
          B0n_I(i)  = B0zFace_z_BLK(iStrip,j,k,globalBLK)
          B0t1_I(i) = B0xFace_z_BLK(iStrip,j,k,globalBLK)
          B0t2_I(i) = B0yFace_z_BLK(iStrip,j,k,globalBLK)
          !\
          ! Left face
          !/
          RhoL_I(i)  =  LeftState_VZ(rho_,iStrip,j,k)
          UnL_I(i)   =  LeftState_VZ(Uz_ ,iStrip,j,k)
          Ut1L_I(i)  =  LeftState_VZ(Ux_ ,iStrip,j,k)
          Ut2L_I(i)  =  LeftState_VZ(Uy_ ,iStrip,j,k)
          B1nL_I(i)  =  LeftState_VZ(Bz_ ,iStrip,j,k)
          B1t1L_I(i) =  LeftState_VZ(Bx_ ,iStrip,j,k)
          B1t2L_I(i) =  LeftState_VZ(By_ ,iStrip,j,k)
          pL_I(i)    =  LeftState_VZ(P_  ,iStrip,j,k)
          !\
          ! Right face
          !/
          RhoR_I(i)  =  RightState_VZ(rho_,iStrip,j,k)
          UnR_I(i)   =  RightState_VZ(Uz_ ,iStrip,j,k)
          Ut1R_I(i)  =  RightState_VZ(Ux_ ,iStrip,j,k)
          Ut2R_I(i)  =  RightState_VZ(Uy_ ,iStrip,j,k)
          B1nR_I(i)  =  RightState_VZ(Bz_ ,iStrip,j,k)
          B1t1R_I(i) =  RightState_VZ(Bx_ ,iStrip,j,k)
          B1t2R_I(i) =  RightState_VZ(By_ ,iStrip,j,k)
          pR_I(i)    =  RightState_VZ(P_  ,iStrip,j,k)

          UseFluxRusanov_I(i)=true_cell(iStrip,j,k-1,globalBLK).neqv.&
               true_cell(iStrip,j,k  ,globalBLK)
       end do
    end select

    do i=1,nStrip
       RhoInvL_I(i) = 1./RhoL_I(i)
       BnL_I(i)  = B0n_I(i)+B1nL_I(i)
       Bt1L_I(i) = B0t1_I(i)+B1t1L_I(i)
       Bt2L_I(i) = B0t2_I(i)+B1t2L_I(i)
       BbL_I(i)  = BnL_I(i)**2 + Bt1L_I(i)**2 + Bt2L_I(i)**2
       Bb1L_I(i) = B1nL_I(i)**2 + B1t1L_I(i)**2 + B1t2L_I(i)**2
       aL_I(i)   = g*pL_I(i)*RhoInvL_I(i)

       RhoInvR_I(i) = 1./RhoR_I(i)
       BnR_I(i)  = B0n_I(i)+B1nR_I(i)
       Bt1R_I(i) = B0t1_I(i)+B1t1R_I(i)
       Bt2R_I(i) = B0t2_I(i)+B1t2R_I(i)
       BbR_I(i)  = BnR_I(i)**2 + Bt1R_I(i)**2 + Bt2R_I(i)**2
       Bb1R_I(i) = B1nR_I(i)**2 + B1t1R_I(i)**2 + B1t2R_I(i)**2
       aR_I(i)   = g*pR_I(i)*RhoInvR_I(i)

       !\
       ! Hat face
       !/
       RhoH_I(i) = 0.5*(RhoL_I(i) + RhoR_I(i))
       RhoInvH_I(i) = 1./RhoH_I(i)
       UnH_I(i)  = 0.5*(  UnL_I(i) +   UnR_I(i))
       Ut1H_I(i) = 0.5*( Ut1L_I(i) +  Ut1R_I(i))
       Ut2H_I(i) = 0.5*( Ut2L_I(i) +  Ut2R_I(i))
       BnH_I(i)  = 0.5*(  BnL_I(i) +   BnR_I(i))
       Bt1H_I(i) = 0.5*( Bt1L_I(i) +  Bt1R_I(i))
       Bt2H_I(i) = 0.5*( Bt2L_I(i) +  Bt2R_I(i))
       B1nH_I(i) = 0.5*( B1nL_I(i) +  B1nR_I(i))
       B1t1H_I(i)= 0.5*(B1t1L_I(i) + B1t1R_I(i))
       B1t2H_I(i)= 0.5*(B1t2L_I(i) + B1t2R_I(i))
       pH_I(i)   = 0.5*(   pL_I(i) +    pR_I(i))
       BbH_I(i)  = BnH_I(i)**2  + Bt1H_I(i)**2  + Bt2H_I(i)**2
       Bb1H_I(i) = B1nH_I(i)**2 + B1t1H_I(i)**2 + B1t2H_I(i)**2
       aH_I(i)   = g*pH_I(i)*RhoInvH_I(i)
    end do

    !VSQRT!    do i=nStrip+1,MaxStrip
    !VSQRT!       aL_I(i) = 0.
    !VSQRT!       aR_I(i) = 0.
    !VSQRT!    end do
    !VSQRT!   len = 2*MaxStrip + nStrip
    !VSQRT!    call sqrt_v32(len,aL_I,aL_I)
    do i=1,nStrip
       !if(aL_I(i)<0.0)then
       !   write(*,*)'NEGATIVE aL_I Me, iDir, i, j, k, globalBLK',&
       !        aL_I(i),iProc,iDir,i,j,k,&
       !        x_BLK(i,j,k,globalBLK),&
       !        y_BLK(i,j,k,globalBLK),&
       !        z_BLK(i,j,k,globalBLK)
       !   call stop_mpi
       !end if
       aL_I(i)=sqrt(aL_I(i))
       aR_I(i)=sqrt(aR_I(i))
       aH_I(i)=sqrt(aH_I(i))
    end do

    do i=1,nStrip
       eL_I(i) = aL_I(i)*aL_I(i) + BbL_I(i)*RhoInvL_I(i)
       CfL_I(i) = max(0.,                                &
            (eL_I(i)**2 - 4.*aL_I(i)**2 * BnL_I(i)**2 * RhoInvL_I(i)))
       eR_I(i) = aR_I(i)**2 + BbR_I(i)*RhoInvR_I(i)
       CfR_I(i) = max(0.,                                &
            (eR_I(i)**2 - 4.*aR_I(i)**2 * BnR_I(i)**2 * RhoInvR_I(i)))
       eH_I(i) = aH_I(i)**2 + BbH_I(i)*RhoInvH_I(i)
       CfH_I(i) = max(0.,                                &
            (eH_I(i)**2 - 4.*aH_I(i)**2 * BnH_I(i)**2 * RhoInvH_I(i)))
    end do

    !VSQRT!    do i=nStrip+1,MaxStrip
    !VSQRT!       CfL_I(i) = 0.
    !VSQRT!       CfR_I(i) = 0.
    !VSQRT!    end do
    !VSQRT!    call sqrt_v32(len,CfL_I,CfL_I)
    do i=1,nStrip
       CfL_I(i)=sqrt(CfL_I(i))
       CfR_I(i)=sqrt(CfR_I(i))
       CfH_I(i)=sqrt(CfH_I(i))
    end do

    do i=1,nStrip
       CsL_I(i)  = max(0.,0.5*(eL_I(i)-CfL_I(i)))
       CfL_I(i)  = 0.5*(eL_I(i)+CfL_I(i))
       eL_I(i) = pL_I(i)*inv_gm1 + 0.5*RhoL_I(i)*          &
            (UnL_I(i)**2 + Ut1L_I(i)**2 + Ut2L_I(i)**2) + 0.5*Bb1L_I(i)
       CsR_I(i)  = max(0.,0.5*(eR_I(i)-CfR_I(i)))
       CfR_I(i)  = 0.5*(eR_I(i)+CfR_I(i))
       eR_I(i)   = pR_I(i)*inv_gm1 + 0.5*RhoR_I(i)*      &
            (UnR_I(i)**2 + Ut1R_I(i)**2 + Ut2R_I(i)**2) + 0.5*Bb1R_I(i)
       CsH_I(i)  = max(0.,0.5*(eH_I(i)-CfH_I(i)))
       CfH_I(i)  = 0.5*(eH_I(i)+CfH_I(i))
       UuH_I(i)  = UnH_I(i)**2 + Ut1H_I(i)**2 + Ut2H_I(i)**2
       eH_I(i)   = pH_I(i)*inv_gm1 + 0.5*RhoH_I(i)*UuH_I(i) + 0.5*Bb1H_I(i)
    end do

    !VSQRT!    do i=nStrip+1,MaxStrip
    !VSQRT!       CsL_I(i) = 0.
    !VSQRT!       CsR_I(i) = 0.
    !VSQRT!       CsH_I(i) = 0.
    !VSQRT!    end do
    !VSQRT!    len = 5*MaxStrip + nStrip
    !VSQRT!    call sqrt_v32(len,CsL_I,CsL_I)
    do i=1,nStrip
       CsL_I(i)=sqrt(CsL_I(i))
       CsR_I(i)=sqrt(CsR_I(i))
       CsH_I(i)=sqrt(CsH_I(i))
       CfL_I(i)=sqrt(CfL_I(i))
       CfR_I(i)=sqrt(CfR_I(i))
       CfH_I(i)=sqrt(CfH_I(i))
    end do

    do i=1,nStrip
       CsL_I(i)  = min(CsL_I(i),aL_I(i))
       CfL_I(i)  = max(CfL_I(i),aL_I(i))
       CsR_I(i)  = min(CsR_I(i),aR_I(i))
       CfR_I(i)  = max(CfR_I(i),aR_I(i))
       CsH_I(i)  = min(CsH_I(i),aH_I(i))
       CfH_I(i)  = max(CfH_I(i),aH_I(i))
       !\
       ! Non-dimensional scaling factors
       !/
       Tmp1_I(i) = max(1.00e-08, Bt1H_I(i)**2 + Bt2H_I(i)**2)
    end do

    !VSQRT!    do i=nStrip+1,MaxStrip
    !VSQRT!       Tmp1_I(i) = 1.0
    !VSQRT!    end do
    !VSQRT!    len = MaxStrip + nStrip
    !VSQRT!    call sqrti_v32(len,Tmp1_I,Tmp1_I)
    do i=1,nStrip
       Tmp1_I(i)=sqrt(1./Tmp1_I(i))
    end do

    do i=1,nStrip
       if (Tmp1_I(i) < 1.0e04) then
          BetaY_I(i) = Bt1H_I(i)*Tmp1_I(i)
          BetaZ_I(i) = Bt2H_I(i)*Tmp1_I(i)
       else
          BetaY_I(i) = cSqrt2Inv
          BetaZ_I(i) = cSqrt2Inv
       end if

       Tmp1 = CfH_I(i)**2 - CsH_I(i)**2
       if (Tmp1 > 1.0e-08) then
          AlphaF_I(i) = (aH_I(i)**2  - CsH_I(i)**2)/Tmp1
          AlphaS_I(i) = (CfH_I(i)**2 - aH_I(i)**2 )/Tmp1
       else if (BnH_I(i)**2 * RhoInvH_I(i) <= aH_I(i)**2 ) then
          AlphaF_I(i) = 1.00
          AlphaS_I(i) = 0.00
       else
          AlphaF_I(i) = 0.00
          AlphaS_I(i) = 1.00
       endif
    end do

    !VSQRT!    do i=nStrip+1,MaxStrip
    !VSQRT!       AlphaF_I(i) = 0.
    !VSQRT!    end do
    !VSQRT!    call sqrt_v32(len,AlphaF_I,AlphaF_I)
    do i=1,nStrip
       AlphaF_I(i)=sqrt(AlphaF_I(i))
       AlphaS_I(i)=sqrt(AlphaS_I(i))
    end do

    !\
    ! Set some values that are reused over and over
    !/
    !VSQRT!    do i=nStrip+1,MaxStrip
    !VSQRT!       RhoSqrtH_I(i) = 0.
    !VSQRT!       RhoSqrtL_I(i) = 0.
    !VSQRT!    end do
    !VSQRT!    len = 2*MaxStrip + nStrip
    !VSQRT!    call sqrt_v32(len,RhoH_I,RhoSqrtH_I)
    do i=1,nStrip
       RhoSqrtH_I(i)   =sqrt(RhoH_I(i))
       RhoSqrtL_I(i)   =sqrt(RhoL_I(i))
       RhoSqrtR_I(i)   =sqrt(RhoR_I(i))
       RhoInvSqrtH_I(i)=1./RhoSqrtH_I(i)
       RhoInvSqrtL_I(i)=1./RhoSqrtL_I(i)
       RhoInvSqrtR_I(i)=1./RhoSqrtR_I(i)
    end do

    do i=1,nStrip
       iStrip = i+iStart-1

       SignBnH     = sign(1.,BnH_I(i))
       Gamma1A2Inv = gm1 / aH_I(i)**2

       !\
       ! Left state face fluxes
       !/
       FluxL_V(Rho_)    = RhoL_I(i)*UnL_I(i)
       FluxL_V(RhoUn_)  = RhoL_I(i)*UnL_I(i)* UnL_I(i) - &
            B1nL_I(i)*B1nL_I(i)+0.5*Bb1L_I(i)+pL_I(i) &
            -((B0n_I(i)*B1nL_I(i))+(B1nL_I(i)*B0n_I(i))) &
            +(B1nL_I(i)*B0n_I(i)+B1t1L_I(i)*B0t1_I(i)+B1t2L_I(i)*B0t2_I(i))
       FluxL_V(RhoUt1_) = RhoL_I(i)*UnL_I(i)*Ut1L_I(i) - &
            B1nL_I(i)*B1t1L_I(i) &
            -((B0n_I(i)*B1t1L_I(i))+(B1nL_I(i)*B0t1_I(i)))
       FluxL_V(RhoUt2_) = RhoL_I(i)*UnL_I(i)*Ut2L_I(i)-  &
            B1nL_I(i)*B1t2L_I(i) &
            -((B0n_I(i)*B1t2L_I(i))+(B1nL_I(i)*B0t2_I(i)))
       FluxL_V(B1n_)    = 0.
       FluxL_V(B1t1_)   = UnL_I(i)*B1t1L_I(i)-B1nL_I(i)*Ut1L_I(i) &
            + UnL_I(i)*B0t1_I(i)-B0n_I(i)*Ut1L_I(i)
       FluxL_V(B1t2_)   = UnL_I(i)*B1t2L_I(i)-B1nL_I(i)*Ut2L_I(i) &
            +UnL_I(i)*B0t2_I(i)-B0n_I(i)*Ut2L_I(i)
       FluxL_V(Energy_)      = UnL_I(i)*(eL_I(i)+pL_I(i)+0.5*Bb1L_I(i)) &
            -B1nL_I(i)*( B1nL_I(i) *UnL_I(i) + &
            B1t1L_I(i)*Ut1L_I(i) + &
            B1t2L_I(i)*Ut2L_I(i)) &
            +UnL_I(i)*(B1nL_I(i)*B0n_I(i) + &
            B1t1L_I(i)*B0t1_I(i) + &
            B1t2L_I(i)*B0t2_I(i)) &
            -B0n_I(i)*(B1nL_I(i)*UnL_I(i) + &
            B1t1L_I(i)*Ut1L_I(i) + &
            B1t2L_I(i)*Ut2L_I(i))
       FluxL_V(P_)      = UnL_I(i)*pL_I(i)

       !\
       ! Right state face fluxes
       !/
       FluxR_V(Rho_)    = RhoR_I(i)*UnR_I(i)
       FluxR_V(RhoUn_)  = RhoR_I(i)*UnR_I(i)* UnR_I(i) - &
            B1nR_I(i)*B1nR_I(i)+0.5*Bb1R_I(i)+pR_I(i) &
            -((B0n_I(i)*B1nR_I(i))+(B1nR_I(i)*B0n_I(i))) &
            +(B1nR_I(i)*B0n_I(i)+B1t1R_I(i)*B0t1_I(i)+B1t2R_I(i)*B0t2_I(i))
       FluxR_V(RhoUt1_) = RhoR_I(i)*UnR_I(i)*Ut1R_I(i) - &
            B1nR_I(i)*B1t1R_I(i) &
            -((B0n_I(i)*B1t1R_I(i))+(B1nR_I(i)*B0t1_I(i)))
       FluxR_V(RhoUt2_) = RhoR_I(i)*UnR_I(i)*Ut2R_I(i)- &
            B1nR_I(i)*B1t2R_I(i) &
            -((B0n_I(i)*B1t2R_I(i))+(B1nR_I(i)*B0t2_I(i)))
       FluxR_V(B1n_)    = 0.
       FluxR_V(B1t1_)   = UnR_I(i)*B1t1R_I(i)-B1nR_I(i)*Ut1R_I(i) &
            +UnR_I(i)*B0t1_I(i)-B0n_I(i)*Ut1R_I(i)
       FluxR_V(B1t2_)   = UnR_I(i)*B1t2R_I(i)-B1nR_I(i)*Ut2R_I(i) &
            +UnR_I(i)*B0t2_I(i)-B0n_I(i)*Ut2R_I(i)
       FluxR_V(Energy_)      = UnR_I(i)*(eR_I(i)+pR_I(i)+0.5*Bb1R_I(i)) &
            -B1nR_I(i)*( B1nR_I(i) *UnR_I(i) + &
            B1t1R_I(i)*Ut1R_I(i) + &
            B1t2R_I(i)*Ut2R_I(i)) &
            +UnR_I(i)*(B1nR_I(i)*B0n_I(i) + &
            B1t1R_I(i)*B0t1_I(i)+ &
            B1t2R_I(i)*B0t2_I(i)) &
            -B0n_I(i)*(B1nR_I(i)*UnR_I(i) + &
            B1t1R_I(i)*Ut1R_I(i) + &
            B1t2R_I(i)*Ut2R_I(i))
       FluxR_V(P_)      = UnR_I(i)*pR_I(i)


       !\
       ! Conserved state delta
       !/
       Delta_V(rho_)    = RhoR_I(i)        - RhoL_I(i)
       Delta_V(rhoUn_)  = RhoR_I(i)* UnR_I(i) - RhoL_I(i)* UnL_I(i)
       Delta_V(rhoUt1_) = RhoR_I(i)*Ut1R_I(i) - RhoL_I(i)*Ut1L_I(i)
       Delta_V(rhoUt2_) = RhoR_I(i)*Ut2R_I(i) - RhoL_I(i)*Ut2L_I(i)
       Delta_V(B1n_)    = B1nR_I(i)         - B1nL_I(i)
       Delta_V(B1t1_)   = B1t1R_I(i)        - B1t1L_I(i)
       Delta_V(B1t2_)   = B1t2R_I(i)        - B1t2L_I(i)
       Delta_V(E_)      = eR_I(i)           - eL_I(i)

       !\
       ! Eigenvalues
       !/
       Eigenvalue_V(EntropyW_) = UnH_I(i)
       Eigenvalue_V(AlfvenRW_) = UnH_I(i) + BnH_I(i)*RhoInvSqrtH_I(i)
       Eigenvalue_V(AlfvenLW_) = UnH_I(i) - BnH_I(i)*RhoInvSqrtH_I(i)
       Eigenvalue_V(SlowRW_)   = UnH_I(i) + CsH_I(i)
       Eigenvalue_V(FastRW_)   = UnH_I(i) + CfH_I(i)
       Eigenvalue_V(SlowLW_)   = UnH_I(i) - CsH_I(i)
       Eigenvalue_V(FastLW_)   = UnH_I(i) - CfH_I(i)
       Eigenvalue_V(DivBW_)    = UnH_I(i)

       !\
       ! Entropy fix for Eigenvalues
       !/
       Tmp1 = (UnR_I(i)) - (UnL_I(i))
       Tmp1 = max(0.,4.*Tmp1)
       if (abs(Eigenvalue_V(1)) < Tmp1*0.5) then
          Eigenvalue_V(1) = sign(1.,Eigenvalue_V(1))*   &
               ((Eigenvalue_V(1)*Eigenvalue_V(1)/Tmp1) + Tmp1*0.25)
       end if

       Tmp1 = (UnR_I(i)+BnR_I(i)*RhoInvSqrtR_I(i)) - &
            (UnL_I(i)+BnL_I(i)*RhoInvSqrtL_I(i))
       Tmp1 = max(0.,4.*Tmp1)
       if (abs(Eigenvalue_V(2)) < Tmp1*0.5) then
          Eigenvalue_V(2) = sign(1.,Eigenvalue_V(2))*   &
               ((Eigenvalue_V(2)*Eigenvalue_V(2)/Tmp1) + Tmp1*0.25)
       end if

       Tmp1 = (UnR_I(i)-BnR_I(i)*RhoInvSqrtR_I(i)) - &
            (UnL_I(i)-BnL_I(i)*RhoInvSqrtL_I(i))
       Tmp1 = max(0.,4.*Tmp1)
       if (abs(Eigenvalue_V(3)) < Tmp1*0.5) then
          Eigenvalue_V(3) = sign(1.,Eigenvalue_V(3))*   &
               ((Eigenvalue_V(3)*Eigenvalue_V(3)/Tmp1) + Tmp1*0.25)
       end if

       Tmp1 = (UnR_I(i)+CsR_I(i)) - (UnL_I(i)+CsL_I(i))
       Tmp1 = max(0.,4.*Tmp1)
       if (abs(Eigenvalue_V(4)) < Tmp1*0.5) then
          Eigenvalue_V(4) = sign(1.,Eigenvalue_V(4))*   &
               ((Eigenvalue_V(4)*Eigenvalue_V(4)/Tmp1) + Tmp1*0.25)
       end if

       Tmp1 = (UnR_I(i)+CfR_I(i)) - (UnL_I(i)+CfL_I(i))
       Tmp1 = max(0.,4.*Tmp1)
       if (abs(Eigenvalue_V(5)) < Tmp1*0.5) then
          Eigenvalue_V(5) = sign(1.,Eigenvalue_V(5))*   &
               ((Eigenvalue_V(5)*Eigenvalue_V(5)/Tmp1) + Tmp1*0.25)
       end if

       Tmp1 = (UnR_I(i)-CsR_I(i)) - (UnL_I(i)-CsL_I(i))
       Tmp1 = max(0.,4.*Tmp1)
       if (abs(Eigenvalue_V(6)) < Tmp1*0.5) then
          Eigenvalue_V(6) = sign(1.,Eigenvalue_V(6))*   &
               ((Eigenvalue_V(6)*Eigenvalue_V(6)/Tmp1) + Tmp1*0.25)
       end if

       Tmp1 = (UnR_I(i)-CfR_I(i)) - (UnL_I(i)-CfL_I(i))
       Tmp1 = max(0.,4.*Tmp1)
       if (abs(Eigenvalue_V(7)) < Tmp1*0.5) then
          Eigenvalue_V(7) = sign(1.,Eigenvalue_V(7))*   &
               ((Eigenvalue_V(7)*Eigenvalue_V(7)/Tmp1) + Tmp1*0.25)
       end if

       Tmp1 = (UnR_I(i)) - (UnL_I(i))
       Tmp1 = max(0.,4.*Tmp1)
       if (abs(Eigenvalue_V(8)) < Tmp1*0.5) then
          Eigenvalue_V(8) = sign(1.,Eigenvalue_V(8))*   &
               ((Eigenvalue_V(8)*Eigenvalue_V(8)/Tmp1) + Tmp1*0.25)
       end if

       !\
       ! Timur's divergence wave fix!!!
       !/
       !original  version
       !      Eigenvalue_V(8)=abs(Eigenvalue_V(8))+aH_I(i)
       !
       !Enhanced diffusion, the maximum eigenvalue
       Eigenvalue_V(DivbW_) = max( Eigenvalue_V(FastRW_), &
            -Eigenvalue_V(FastLW_))

       ! The original version was proposed by Timur Linde for heliosphere
       ! simulations. Enhanced version was found to be more robust on 8/20/01
       ! The original version was commented out in versions 6x resulting in 
       ! worse stability for the Roe solver.

       ! At inner BC replace all eigenvalues with the enhanced eigenvalue
       ! of divB, which is the maximum eigenvalue
       if(UseFluxRusanov_I(i)) Eigenvalue_V(1:nWave) = Eigenvalue_V(DivBW_)

       !\
       ! Eigenvectors
       !/
       Tmp1=1./(2.*RhoH_I(i)*aH_I(i)**2)
       Tmp2=RhoInvH_I(i)*cSqrt2Inv
       Tmp3=RhoInvSqrtH_I(i)*cSqrt2Inv

       ! Left eigenvector for Entropy wave
       EigenvectorL_VV(1,1) = 1.-0.5*Gamma1A2Inv*UuH_I(i)
       EigenvectorL_VV(2,1) = Gamma1A2Inv*UnH_I(i)
       EigenvectorL_VV(3,1) = Gamma1A2Inv*Ut1H_I(i)
       EigenvectorL_VV(4,1) = Gamma1A2Inv*Ut2H_I(i)
       EigenvectorL_VV(5,1) = Gamma1A2Inv*B1nH_I(i)
       EigenvectorL_VV(6,1) = Gamma1A2Inv*B1t1H_I(i)
       EigenvectorL_VV(7,1) = Gamma1A2Inv*B1t2H_I(i)
       EigenvectorL_VV(8,1) = -Gamma1A2Inv

       ! Left eigenvector for Alfven wave +
       EigenvectorL_VV(1,2) = (Ut1H_I(i)*BetaZ_I(i)-  &
            Ut2H_I(i)*BetaY_I(i))*Tmp2
       EigenvectorL_VV(2,2) = 0.
       EigenvectorL_VV(3,2) = -(BetaZ_I(i)*Tmp2)
       EigenvectorL_VV(4,2) = (BetaY_I(i)*Tmp2)
       EigenvectorL_VV(5,2) = 0.
       EigenvectorL_VV(6,2) = (BetaZ_I(i)*Tmp3)
       EigenvectorL_VV(7,2) = -(BetaY_I(i)*Tmp3)
       EigenvectorL_VV(8,2) = 0.

       ! Left eigenvector for Alfven wave -
       EigenvectorL_VV(1,3) = (Ut1H_I(i)*BetaZ_I(i)-  &
            Ut2H_I(i)*BetaY_I(i))*Tmp2
       EigenvectorL_VV(2,3) = 0.
       EigenvectorL_VV(3,3) = -(BetaZ_I(i)*Tmp2)
       EigenvectorL_VV(4,3) = (BetaY_I(i)*Tmp2)
       EigenvectorL_VV(5,3) = 0.
       EigenvectorL_VV(6,3) = -(BetaZ_I(i)*Tmp3)
       EigenvectorL_VV(7,3) = (BetaY_I(i)*Tmp3)
       EigenvectorL_VV(8,3) = 0.

       ! Left eigenvector for Slow magnetosonic wave +
       EigenvectorL_VV(1,4) = Tmp1* &
            (AlphaS_I(i)*(gm1*UuH_I(i)/2. - &
            UnH_I(i)*CsH_I(i)) - &
            AlphaF_I(i)*CfH_I(i)*SignBnH* &
            (Ut1H_I(i)*BetaY_I(i)+ &
            Ut2H_I(i)*BetaZ_I(i)))
       EigenvectorL_VV(2,4) = Tmp1* &
            (AlphaS_I(i)*(-UnH_I(i)*gm1+CsH_I(i)))
       EigenvectorL_VV(3,4) = Tmp1* &
            (-gm1*AlphaS_I(i)*Ut1H_I(i) + &
            AlphaF_I(i)*CfH_I(i)*BetaY_I(i)*SignBnH)
       EigenvectorL_VV(4,4) = Tmp1* &
            (-gm1*AlphaS_I(i)*Ut2H_I(i) + &
            AlphaF_I(i)*CfH_I(i)*BetaZ_I(i)*SignBnH)
       EigenvectorL_VV(5,4) = Tmp1* &
            (-gm1*B1nH_I(i)*AlphaS_I(i))
       EigenvectorL_VV(6,4) = Tmp1* &
            (-AlphaF_I(i)*BetaY_I(i)*aH_I(i)* &
            RhoSqrtH_I(i)-gm1*B1t1H_I(i)*AlphaS_I(i))
       EigenvectorL_VV(7,4) = Tmp1* &
            (-AlphaF_I(i)*BetaZ_I(i)*aH_I(i)* &
            RhoSqrtH_I(i)-gm1*B1t2H_I(i)*AlphaS_I(i))
       EigenvectorL_VV(8,4) = Tmp1* &
            (gm1*AlphaS_I(i))

       ! Left eigenvector for Fast magnetosonic wave +
       EigenvectorL_VV(1,5) = Tmp1* &
            (AlphaF_I(i)*(gm1*UuH_I(i)/2. - &
            UnH_I(i)*CfH_I(i))+AlphaS_I(i)* &
            CsH_I(i)*SignBnH* &
            (Ut1H_I(i)*BetaY_I(i) + &
            Ut2H_I(i)*BetaZ_I(i)))
       EigenvectorL_VV(2,5) = Tmp1* &
            (AlphaF_I(i)*(-UnH_I(i)*gm1+CfH_I(i)))
       EigenvectorL_VV(3,5) = Tmp1* &
            (-gm1*AlphaF_I(i)*Ut1H_I(i) - &
            AlphaS_I(i)*CsH_I(i)*BetaY_I(i)*SignBnH)
       EigenvectorL_VV(4,5) = Tmp1* &
            (-gm1*AlphaF_I(i)*Ut2H_I(i) - &
            AlphaS_I(i)*CsH_I(i)*BetaZ_I(i)*SignBnH)
       EigenvectorL_VV(5,5) = Tmp1* &
            (-gm1*B1nH_I(i)*AlphaF_I(i))
       EigenvectorL_VV(6,5) = Tmp1* &
            (AlphaS_I(i)*BetaY_I(i)*aH_I(i)* &
            RhoSqrtH_I(i)-gm1*B1t1H_I(i)*AlphaF_I(i))
       EigenvectorL_VV(7,5) = Tmp1*                   &
            (AlphaS_I(i)*BetaZ_I(i)*aH_I(i)* &
            RhoSqrtH_I(i)-gm1*B1t2H_I(i)*AlphaF_I(i))
       EigenvectorL_VV(8,5) = Tmp1* &
            (gm1*AlphaF_I(i))

       ! Left eigenvector for Slow magnetosonic wave -
       EigenvectorL_VV(1,6) = Tmp1* &
            (AlphaS_I(i)*(gm1*UuH_I(i)/2. + &
            UnH_I(i)*CsH_I(i)) + &
            AlphaF_I(i)*CfH_I(i)*SignBnH* &
            (Ut1H_I(i)*BetaY_I(i) + &
            Ut2H_I(i)*BetaZ_I(i)))
       EigenvectorL_VV(2,6) = Tmp1* &
            (AlphaS_I(i)*(-UnH_I(i)*gm1-CsH_I(i)))
       EigenvectorL_VV(3,6) = Tmp1* &
            (-gm1*AlphaS_I(i)*Ut1H_I(i) - &
            AlphaF_I(i)*CfH_I(i)*BetaY_I(i)*SignBnH)
       EigenvectorL_VV(4,6) = Tmp1* &
            (-gm1*AlphaS_I(i)*Ut2H_I(i) - &
            AlphaF_I(i)*CfH_I(i)*BetaZ_I(i)*SignBnH)
       EigenvectorL_VV(5,6) = Tmp1* &
            (-gm1*B1nH_I(i)*AlphaS_I(i))
       EigenvectorL_VV(6,6) = Tmp1* &
            (-AlphaF_I(i)*BetaY_I(i)*aH_I(i)* &
            RhoSqrtH_I(i)-gm1*B1t1H_I(i)*AlphaS_I(i))
       EigenvectorL_VV(7,6) = Tmp1* &
            (-AlphaF_I(i)*BetaZ_I(i)*aH_I(i)* &
            RhoSqrtH_I(i)-gm1*B1t2H_I(i)*AlphaS_I(i))
       EigenvectorL_VV(8,6) = Tmp1* &
            (gm1*AlphaS_I(i))

       ! Left eigenvector for Fast magnetosonic wave -
       EigenvectorL_VV(1,7) = Tmp1* &
            (AlphaF_I(i)*(gm1*UuH_I(i)/2. + &
            UnH_I(i)*CfH_I(i)) - &
            AlphaS_I(i)*CsH_I(i)*SignBnH* &
            (Ut1H_I(i)*BetaY_I(i) + &
            Ut2H_I(i)*BetaZ_I(i)))
       EigenvectorL_VV(2,7) = Tmp1* &
            (AlphaF_I(i)*(-UnH_I(i)*gm1-CfH_I(i)))
       EigenvectorL_VV(3,7) = Tmp1* &
            (-gm1*AlphaF_I(i)*Ut1H_I(i) + &
            AlphaS_I(i)*CsH_I(i)*BetaY_I(i)*SignBnH)
       EigenvectorL_VV(4,7) = Tmp1* &
            (-gm1*AlphaF_I(i)*Ut2H_I(i) + &
            AlphaS_I(i)*CsH_I(i)*BetaZ_I(i)*SignBnH)
       EigenvectorL_VV(5,7) = Tmp1* &
            (-gm1*B1nH_I(i)*AlphaF_I(i))
       EigenvectorL_VV(6,7) = Tmp1* &
            (AlphaS_I(i)*BetaY_I(i)*aH_I(i)* &
            RhoSqrtH_I(i)-gm1*B1t1H_I(i)*AlphaF_I(i))
       EigenvectorL_VV(7,7) = Tmp1* &
            (AlphaS_I(i)*BetaZ_I(i)*aH_I(i)* &
            RhoSqrtH_I(i)-gm1*B1t2H_I(i)*AlphaF_I(i))
       EigenvectorL_VV(8,7) = Tmp1* &
            (gm1*AlphaF_I(i))

       ! Left eigenvector for Divergence wave
       EigenvectorL_VV(1,8) = 0.
       EigenvectorL_VV(2,8) = 0.
       EigenvectorL_VV(3,8) = 0.
       EigenvectorL_VV(4,8) = 0.
       EigenvectorL_VV(5,8) = 1.
       EigenvectorL_VV(6,8) = 0.
       EigenvectorL_VV(7,8) = 0.
       EigenvectorL_VV(8,8) = 0.

       Tmp1_I(i)=g*max(pL_I(i),pR_I(i)) !coefficient for pressure component of the Right vector

       !Pressure component is not linearly independent and obeys the equation as follows:
       !EigenvectorR_VV(1:8,9)=(0.5*UuH_I(i)*EigenvectorR_VV(1:8,1)-&
       !                     UnH_I(i)*EigenvectorR_VV(1:8,2)-&
       !                     Ut1H_I(i)*EigenvectorR_VV(1:8,3)-&
       !                     Ut2H_I(i)*EigenvectorR_VV(1:8,4)-&
       !                     B1nH_I(i)*EigenvectorR_VV(1:8,5)-&
       !                     B1t1H_I(i)*EigenvectorR_VV(1:8,6)-&
       !                     B1t2H_I(i)*EigenvectorR_VV(1:8,7)+
       !                     EigenvectorR_VV(1:8,8))*inv_gm1         

       ! Right eigenvector for Entropy wave
       EigenvectorR_VV(1,1) = 1.
       EigenvectorR_VV(1,2) = UnH_I(i)
       EigenvectorR_VV(1,3) = Ut1H_I(i)
       EigenvectorR_VV(1,4) = Ut2H_I(i)
       EigenvectorR_VV(1,5) = 0.
       EigenvectorR_VV(1,6) = 0.
       EigenvectorR_VV(1,7) = 0.
       EigenvectorR_VV(1,Energy_) = 0.5*UuH_I(i)
       EigenvectorR_VV(1,P_)=cZero

       ! Right eigenvector for Alfven wave +
       EigenvectorR_VV(2,1) = 0.
       EigenvectorR_VV(2,2) = 0.
       EigenvectorR_VV(2,3) = -BetaZ_I(i)*RhoH_I(i)*cSqrt2Inv
       EigenvectorR_VV(2,4) = BetaY_I(i)*RhoH_I(i)*cSqrt2Inv
       EigenvectorR_VV(2,5) = 0.
       EigenvectorR_VV(2,6) = BetaZ_I(i)*RhoSqrtH_I(i)*cSqrt2Inv
       EigenvectorR_VV(2,7) = -BetaY_I(i)*RhoSqrtH_I(i)*cSqrt2Inv
       EigenvectorR_VV(2,Energy_) = (BetaY_I(i)*Ut2H_I(i) - &
            BetaZ_I(i)*Ut1H_I(i))*RhoH_I(i)*cSqrt2Inv &
            +(B1t1H_I(i)*BetaZ_I(i)-B1t2H_I(i)*BetaY_I(i))* &
            RhoSqrtH_I(i)*cSqrt2Inv
       EigenvectorR_VV(2,P_)=cZero

       ! Right eigenvector for Alfven wave -
       EigenvectorR_VV(3,1) = 0.
       EigenvectorR_VV(3,2) = 0.
       EigenvectorR_VV(3,3) = -BetaZ_I(i)*RhoH_I(i)*cSqrt2Inv
       EigenvectorR_VV(3,4) = BetaY_I(i)*RhoH_I(i)*cSqrt2Inv
       EigenvectorR_VV(3,5) = 0.
       EigenvectorR_VV(3,6) = -BetaZ_I(i)*RhoSqrtH_I(i)*cSqrt2Inv
       EigenvectorR_VV(3,7) = BetaY_I(i)*RhoSqrtH_I(i)*cSqrt2Inv
       EigenvectorR_VV(3,Energy_) = (BetaY_I(i)*Ut2H_I(i) - &
            BetaZ_I(i)*Ut1H_I(i))*RhoH_I(i)*cSqrt2Inv &
            -(B1t1H_I(i)*BetaZ_I(i)-B1t2H_I(i)*BetaY_I(i))* &
            RhoSqrtH_I(i)*cSqrt2Inv
       EigenvectorR_VV(3,P_)=cZero

       ! Right eigenvector for Slow magnetosonic wave +
       EigenvectorR_VV(4,1) = RhoH_I(i)*AlphaS_I(i)
       EigenvectorR_VV(4,2) = RhoH_I(i)*AlphaS_I(i)* &
            (UnH_I(i)+CsH_I(i))
       EigenvectorR_VV(4,3) = RhoH_I(i)* &
            (AlphaS_I(i)*Ut1H_I(i) + &
            AlphaF_I(i)*CfH_I(i)*BetaY_I(i)*SignBnH)
       EigenvectorR_VV(4,4) = RhoH_I(i)* &
            (AlphaS_I(i)*Ut2H_I(i) + &
            AlphaF_I(i)*CfH_I(i)*BetaZ_I(i)*SignBnH)
       EigenvectorR_VV(4,5) = 0.
       EigenvectorR_VV(4,6) = -AlphaF_I(i)*aH_I(i)*BetaY_I(i)* &
            RhoSqrtH_I(i)
       EigenvectorR_VV(4,7) = -AlphaF_I(i)*aH_I(i)*BetaZ_I(i)* &
            RhoSqrtH_I(i)
       EigenvectorR_VV(4,Energy_) = AlphaS_I(i)*(RhoH_I(i)*UuH_I(i)*0.5 + &
            g*pH_I(i)*inv_gm1+RhoH_I(i)*UnH_I(i)* &
            CsH_I(i))-AlphaF_I(i)*(aH_I(i)* &
            RhoSqrtH_I(i)*(BetaY_I(i)*B1t1H_I(i) + &
            BetaZ_I(i)*B1t2H_I(i))-RhoH_I(i)* &
            CfH_I(i)*SignBnH* &
            (Ut1H_I(i)*BetaY_I(i)+Ut2H_I(i)*BetaZ_I(i)))
       EigenvectorR_VV(4,P_)=Tmp1_I(i)*AlphaS_I(i)

       ! Right eigenvector for Fast magnetosonic wave +
       EigenvectorR_VV(5,1) = RhoH_I(i)*AlphaF_I(i)
       EigenvectorR_VV(5,2) = RhoH_I(i)*AlphaF_I(i)* &
            (UnH_I(i)+CfH_I(i))
       EigenvectorR_VV(5,3) = RhoH_I(i)* &
            (AlphaF_I(i)*Ut1H_I(i) - &
            AlphaS_I(i)*CsH_I(i)*BetaY_I(i)*SignBnH)
       EigenvectorR_VV(5,4) = RhoH_I(i)* &
            (AlphaF_I(i)*Ut2H_I(i) - &
            AlphaS_I(i)*CsH_I(i)*BetaZ_I(i)*SignBnH)
       EigenvectorR_VV(5,5) = 0.
       EigenvectorR_VV(5,6) = AlphaS_I(i)*aH_I(i)*BetaY_I(i)* &
            RhoSqrtH_I(i)
       EigenvectorR_VV(5,7) = AlphaS_I(i)*aH_I(i)*BetaZ_I(i)* &
            RhoSqrtH_I(i)
       EigenvectorR_VV(5,Energy_) = AlphaF_I(i)*(RhoH_I(i)*UuH_I(i)*0.5 + &
            g*pH_I(i)*inv_gm1+RhoH_I(i)*UnH_I(i)* &
            CfH_I(i))+AlphaS_I(i)*(aH_I(i)* &
            RhoSqrtH_I(i)*(BetaY_I(i)*B1t1H_I(i)+ &
            BetaZ_I(i)*B1t2H_I(i))-RhoH_I(i)* &
            CsH_I(i)*SignBnH* &
            (Ut1H_I(i)*BetaY_I(i)+Ut2H_I(i)*BetaZ_I(i)))
       EigenvectorR_VV(5,P_)=Tmp1_I(i)*AlphaF_I(i)

       ! Right eigenvector for Slow magnetosonic wave -
       EigenvectorR_VV(6,1) = RhoH_I(i)*AlphaS_I(i)
       EigenvectorR_VV(6,2) = RhoH_I(i)*AlphaS_I(i)* &
            (UnH_I(i)-CsH_I(i))
       EigenvectorR_VV(6,3) = RhoH_I(i)* &
            (AlphaS_I(i)*Ut1H_I(i) - &
            AlphaF_I(i)*CfH_I(i)*BetaY_I(i)*SignBnH)
       EigenvectorR_VV(6,4) = RhoH_I(i)* &
            (AlphaS_I(i)*Ut2H_I(i) - &
            AlphaF_I(i)*CfH_I(i)*BetaZ_I(i)*SignBnH)
       EigenvectorR_VV(6,5) = 0.
       EigenvectorR_VV(6,6) = - AlphaF_I(i)*aH_I(i)*BetaY_I(i)* &
            RhoSqrtH_I(i)
       EigenvectorR_VV(6,7) = - AlphaF_I(i)*aH_I(i)*BetaZ_I(i)* &
            RhoSqrtH_I(i)
       EigenvectorR_VV(6,Energy_) = AlphaS_I(i)*(RhoH_I(i)*UuH_I(i)*0.5 + &
            g*pH_I(i)*inv_gm1-RhoH_I(i)*UnH_I(i)* &
            CsH_I(i))-AlphaF_I(i)*(aH_I(i)* &
            RhoSqrtH_I(i)*(BetaY_I(i)*B1t1H_I(i)+ &
            BetaZ_I(i)*B1t2H_I(i))+RhoH_I(i)* &
            CfH_I(i)*SignBnH* &
            (Ut1H_I(i)*BetaY_I(i)+Ut2H_I(i)*BetaZ_I(i)))
       EigenvectorR_VV(6,P_)=Tmp1_I(i)*AlphaS_I(i)

       ! Right eigenvector for Fast magnetosonic wave -
       EigenvectorR_VV(7,1) = RhoH_I(i)*AlphaF_I(i)
       EigenvectorR_VV(7,2) = RhoH_I(i)*AlphaF_I(i)* &
            (UnH_I(i)-CfH_I(i))
       EigenvectorR_VV(7,3) = RhoH_I(i)* &
            (AlphaF_I(i)*Ut1H_I(i) + &
            AlphaS_I(i)*CsH_I(i)*BetaY_I(i)*SignBnH)
       EigenvectorR_VV(7,4) = RhoH_I(i)* &
            (AlphaF_I(i)*Ut2H_I(i) + &
            AlphaS_I(i)*CsH_I(i)*BetaZ_I(i)*SignBnH)
       EigenvectorR_VV(7,5) = 0.
       EigenvectorR_VV(7,6) = AlphaS_I(i)*aH_I(i)*BetaY_I(i)* &
            RhoSqrtH_I(i)
       EigenvectorR_VV(7,7) = AlphaS_I(i)*aH_I(i)*BetaZ_I(i)* &
            RhoSqrtH_I(i)
       EigenvectorR_VV(7,Energy_) = AlphaF_I(i)*(RhoH_I(i)*UuH_I(i)*0.5 + &
            g*pH_I(i)*inv_gm1-RhoH_I(i)*UnH_I(i)* &
            CfH_I(i))+AlphaS_I(i)*(aH_I(i)* &
            RhoSqrtH_I(i)*(BetaY_I(i)*B1t1H_I(i) + &
            BetaZ_I(i)*B1t2H_I(i))+RhoH_I(i)* &
            CsH_I(i)*SignBnH* &
            (Ut1H_I(i)*BetaY_I(i)+Ut2H_I(i)*BetaZ_I(i)))
       EigenvectorR_VV(7,P_)=Tmp1_I(i)*AlphaF_I(i)

       ! Right eigenvector for Divergence wave
       EigenvectorR_VV(8,1) = 0.
       EigenvectorR_VV(8,2) = 0.
       EigenvectorR_VV(8,3) = 0.
       EigenvectorR_VV(8,4) = 0.
       EigenvectorR_VV(8,5) = 1.
       EigenvectorR_VV(8,6) = 0.
       EigenvectorR_VV(8,7) = 0.
       EigenvectorR_VV(8,Energy_) = B1nH_I(i)
       EigenvectorR_VV(8,P_) = cZero

       !\
       ! Alphas (elemental wave strengths)
       !/
       do iWave=1,nWave
          DeltaWave_V(iWave)=sum(EigenvectorL_VV(1:nVar,iWave)*Delta_V(1:nVar))
       end do

       !\
       ! Calculate the Roe Interface fluxes 
       ! F = A * 0.5 * [ F_L+F_R - sum_k(|lambda_k| * alpha_k * r_k) ]
       !/
       do iFlux=1,nFlux
          FluxFull_V(iFlux) = 0.5*AreaFace*((FluxL_V(iFlux)+FluxR_V(iFlux)) - &
               sum(abs(Eigenvalue_V(1:nWave))*DeltaWave_V(1:nWave)* &
               EigenvectorR_VV(1:nWave,iFlux)))
       end do

       !\
       ! Compute face time step (actually Volume/dt)
       !/
       DtInvVolume = AreaFace * (abs(UnH_I(i))+CfH_I(i))

       !\
       ! Move interface fluxes to face fluxes for use elsewhere
       !   (also DtInvVolume)
       !/
       select case (iDir)
       case (1)
          !\
          ! x face
          !/
          Flux_VX(rho_  ,iStrip,j,k) = FluxFull_V(1)
          Flux_VX(rhoUx_,iStrip,j,k) = FluxFull_V(2)
          Flux_VX(rhoUy_,iStrip,j,k) = FluxFull_V(3)
          Flux_VX(rhoUz_,iStrip,j,k) = FluxFull_V(4)
          Flux_VX(Bx_   ,iStrip,j,k) = FluxFull_V(5)
          Flux_VX(By_   ,iStrip,j,k) = FluxFull_V(6)
          Flux_VX(Bz_   ,iStrip,j,k) = FluxFull_V(7)
          Flux_VX(P_    ,iStrip,j,k) = FluxFull_V(8)
          Flux_VX(Energy_    ,iStrip,j,k) = FluxFull_V(9)

          VdtFace_x(iStrip,j,k)    = DtInvVolume
          UDotFA_X(iStrip,j,k)     = &
               cHalf*AreaFace*(UnL_I(i)+UnR_I(i))
       case (2)
          !\
          ! y face
          !/
          Flux_VY(rho_  ,iStrip,j,k) = FluxFull_V(1)
          Flux_VY(rhoUx_,iStrip,j,k) = FluxFull_V(4)
          Flux_VY(rhoUy_,iStrip,j,k) = FluxFull_V(2)
          Flux_VY(rhoUz_,iStrip,j,k) = FluxFull_V(3)
          Flux_VY(Bx_   ,iStrip,j,k) = FluxFull_V(7)
          Flux_VY(By_   ,iStrip,j,k) = FluxFull_V(5)
          Flux_VY(Bz_   ,iStrip,j,k) = FluxFull_V(6)
          Flux_VY(P_    ,iStrip,j,k) = FluxFull_V(8)
          Flux_VY(Energy_,iStrip,j,k)     = FluxFull_V(9)

          VdtFace_y(iStrip,j,k)    = DtInvVolume
          UDotFA_Y(iStrip,j,k)     = &
               cHalf*AreaFace*(UnL_I(i)+UnR_I(i))
       case (3)
          !\
          ! z face
          !/
          Flux_VZ(rho_  ,iStrip,j,k) = FluxFull_V(1)
          Flux_VZ(rhoUx_,iStrip,j,k) = FluxFull_V(3)
          Flux_VZ(rhoUy_,iStrip,j,k) = FluxFull_V(4)
          Flux_VZ(rhoUz_,iStrip,j,k) = FluxFull_V(2)
          Flux_VZ(Bx_   ,iStrip,j,k) = FluxFull_V(6)
          Flux_VZ(By_   ,iStrip,j,k) = FluxFull_V(7)
          Flux_VZ(Bz_   ,iStrip,j,k) = FluxFull_V(5)
          Flux_VZ(P_    ,iStrip,j,k) = FluxFull_V(8)
          Flux_VZ(Energy_    ,iStrip,j,k) = FluxFull_V(9)

          VdtFace_z(iStrip,j,k)    = DtInvVolume
          UDotFA_z(iStrip,j,k)     = &
               cHalf*AreaFace*(UnL_I(i)+UnR_I(i))
       end select
    end do
  end subroutine get_fluxes

end subroutine calc_flux_Roe


!==========================================================================
!==========================================================================
!==========================================================================

!\\\
!<ddz> Sample print statements for debugging
!
! if (iDir == 1 .and. (iStrip) == 1 .and. j == 1 .and. k == 1) then
!    write (*,'(a,3i2,a,f16.10)') "  OPT flux  I,J,K=",iStrip,j,k,", time=",DtInvVolume
!    write (*,'(4x,2f18.10,2x,f18.10)') FluxL_V(1),FluxR_V(1),FluxFull_V(1)
!    write (*,'(4x,2f18.10,2x,f18.10)') FluxL_V(2),FluxR_V(2),FluxFull_V(2)
!    write (*,'(4x,2f18.10,2x,f18.10)') FluxL_V(3),FluxR_V(3),FluxFull_V(3)
!    write (*,'(4x,2f18.10,2x,f18.10)') FluxL_V(4),FluxR_V(4),FluxFull_V(4)
!    write (*,'(4x,2f18.10,2x,f18.10)') FluxL_V(5),FluxR_V(5),FluxFull_V(5)
!    write (*,'(4x,2f18.10,2x,f18.10)') FluxL_V(6),FluxR_V(6),FluxFull_V(6)
!    write (*,'(4x,2f18.10,2x,f18.10)') FluxL_V(7),FluxR_V(7),FluxFull_V(7)
!    write (*,'(4x,2f18.10,2x,f18.10)') FluxL_V(8),FluxR_V(8),FluxFull_V(8)
!    write (*,*) ' '
!    write (*,*) '  I                LF                RF                 HAT'
!    write (*,*) ' '
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  0,RhoL_I(i),RhoR_I(i),RhoH_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  1,UnL_I(i),UnR_I(i),UnH_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  2,Ut1L_I(i),Ut1R_I(i),Ut1H_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  3,Ut2L_I(i),Ut2R_I(i),Ut2H_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  4,pL_I(i),pR_I(i),pH_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  5,BnL_I(i),BnR_I(i),BnH_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  6,Bt1L_I(i),Bt1R_I(i),Bt1H_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  7,Bt2L_I(i),Bt2R_I(i),Bt2H_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  8,
!    write (*,'(2x,i2,2f18.10,2x,f18.10)')  9,UuH_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)') 10,BbL_I(i),BbR_I(i),BbH_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)') 11,aL_I(i),aR_I(i),aH_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)') 12,CsL_I(i),CsR_I(i),CsH_I(i)
!    write (*,'(2x,i2,2f18.10,2x,f18.10)') 13,CfL_I(i),CfR_I(i),CfH_I(i)
!    write (*,*) ' '
!    write (*,*) '       Eigenvalue_V             DeltaWave_V               DEL'
!    write (*,*) ' '
!    write (*,'(3f18.10)') abs(Eigenvalue_V(1)),DeltaWave_V(1),Delta_V(1)
!    write (*,'(3f18.10)') abs(Eigenvalue_V(2)),DeltaWave_V(2),Delta_V(2)
!    write (*,'(3f18.10)') abs(Eigenvalue_V(3)),DeltaWave_V(3),Delta_V(3)
!    write (*,'(3f18.10)') abs(Eigenvalue_V(4)),DeltaWave_V(4),Delta_V(4)
!    write (*,'(3f18.10)') abs(Eigenvalue_V(5)),DeltaWave_V(5),Delta_V(5)
!    write (*,'(3f18.10)') abs(Eigenvalue_V(6)),DeltaWave_V(6),Delta_V(6)
!    write (*,'(3f18.10)') abs(Eigenvalue_V(7)),DeltaWave_V(7),Delta_V(7)
!    write (*,'(3f18.10)') abs(Eigenvalue_V(8)),DeltaWave_V(8),Delta_V(8)
!    write (*,*) ' '
!    write (*,*) 'J,I     EigenvectorL_VV  I,J     EigenvectorR_VV'
!    write (*,*) ' '
!    do iFlux=1,8
!       write (*,'(" 1,",i1,f18.10,i3,",1",f18.10)') iFlux,EigenvectorL_VV(1,iFlux),iFlux,EigenvectorR_VV(iFlux,1)
!       write (*,'(" 2,",i1,f18.10,i3,",2",f18.10)') iFlux,EigenvectorL_VV(2,iFlux),iFlux,EigenvectorR_VV(iFlux,2)
!       write (*,'(" 3,",i1,f18.10,i3,",3",f18.10)') iFlux,EigenvectorL_VV(3,iFlux),iFlux,EigenvectorR_VV(iFlux,3)
!       write (*,'(" 4,",i1,f18.10,i3,",4",f18.10)') iFlux,EigenvectorL_VV(4,iFlux),iFlux,EigenvectorR_VV(iFlux,4)
!       write (*,'(" 5,",i1,f18.10,i3,",5",f18.10)') iFlux,EigenvectorL_VV(5,iFlux),iFlux,EigenvectorR_VV(iFlux,5)
!       write (*,'(" 6,",i1,f18.10,i3,",6",f18.10)') iFlux,EigenvectorL_VV(6,iFlux),iFlux,EigenvectorR_VV(iFlux,6)
!       write (*,'(" 7,",i1,f18.10,i3,",7",f18.10)') iFlux,EigenvectorL_VV(7,iFlux),iFlux,EigenvectorR_VV(iFlux,7)
!       write (*,'(" 8,",i1,f18.10,i3,",8",f18.10)') iFlux,EigenvectorL_VV(8,iFlux),iFlux,EigenvectorR_VV(iFlux,8)
!       write (*,*) ' '
!    end do
!    write (*,*) ' '
!    write (*,'(a,f18.10)') "betay:  ",BetaY_I(i)
!    write (*,'(a,f18.10)') "betaz:  ",BetaZ_I(i)
!    write (*,'(a,f18.10)') "alphas: ",AlphaS_I(i)
!    write (*,'(a,f18.10)') "alphaf: ",AlphaF_I(i)
!    write (*,*) ' '
! end if
!///    
















