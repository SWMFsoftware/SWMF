!^CFG FILE LINDEFLUX
!^CFG COPYRIGHT UM
subroutine option_lindeflux(TrueOption,NameOption)

  logical, intent(out) :: TrueOption
  character (len=40), intent(out) :: NameOption

  TrueOption  = .true.
  NameOption  = 'LINDE FLUX Linde, Sokolov 1.1'

end subroutine option_lindeflux

subroutine calc_flux_Linde(DoResChangeOnly)
  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModGeometry, ONLY : fAx_BLK,fAy_BLK,fAz_BLK
  use ModPhysics, ONLY : g,inv_gm1,cLIGHT,inv_c2LIGHT
  use ModNumConst
  use ModAdvance
  use ModParallel, ONLY : neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  implicit none

  logical, intent (in) :: DoResChangeOnly

  integer, parameter :: nFlux=Energy_

  integer, parameter :: MaxStrip=10

  integer :: i,j,k,n,iDir,iStart,nStrip,iStrip

  ! B0 at face

  real, dimension(MaxStrip) :: v_B0n, v_B0t1, v_B0t2

  ! Left face

  real, dimension(MaxStrip) :: v_rho_lf,v_Un_lf,v_Ut1_lf,v_Ut2_lf
  real, dimension(MaxStrip) :: v_Bn_lf,v_Bt1_lf,v_Bt2_lf,v_BB_lf
  real, dimension(MaxStrip) :: v_B1n_lf,v_B1t1_lf,v_B1t2_lf,v_BB1_lf
  real, dimension(MaxStrip) :: v_p_lf,v_E_lf
  real, dimension(MaxStrip) :: v_En_lf,v_Et1_lf,v_Et2_lf,v_EE_lf
  real :: UU_lf, UdotB_lf

  ! Right face

  real, dimension(MaxStrip) :: v_rho_rf,v_Un_rf,v_Ut1_rf,v_Ut2_rf
  real, dimension(MaxStrip) :: v_Bn_rf,v_Bt1_rf,v_Bt2_rf,v_BB_rf
  real, dimension(MaxStrip) :: v_B1n_rf,v_B1t1_rf,v_B1t2_rf,v_BB1_rf
  real, dimension(MaxStrip) :: v_p_rf,v_E_rf
  real, dimension(MaxStrip) :: v_En_rf,v_Et1_rf,v_Et2_rf,v_EE_rf
  real :: UU_rf, UdotB_rf

  ! Hat face

  real, dimension(MaxStrip) :: v_rho_hf, v_Un_hf, v_p_hf, v_cmax_hf
  real, dimension(MaxStrip) :: v_Bn_hf, v_Bt1_hf, v_Bt2_hf, v_BB_hf
  real, dimension(MaxStrip) :: v_B1n_hf
  real, dimension(MaxStrip) :: v_max_hf
  real, dimension(MaxStrip) :: v_diffrB1n

  ! Wave speed variables
  real :: inv_rho, a2, vA2, a2Va2, vAn2, discr, discr_slow, discr_fast
  !^CFG IF BORISCORR BEGIN
  ! Boris
  real :: vA2Boris, vAn2Boris, a2Boris, a2Va2Boris, gammaA2, gammaU2
  !^CFG END BORISCORR
  !^CFG IF SIMPLEBORIS BEGIN
  ! Simple Boris
  real :: Gamma2R, Gamma2L
  !^CFG END SIMPLEBORIS
  ! Flux, Jump
  real, dimension(nFlux,MaxStrip) :: L_Flux,R_Flux,delta_U
  real, dimension(nFlux+1&              !+1 for u_n averaged
       +1&                              !+1 for E_n averaged !^CFG IF BORISCORR
       ,MaxStrip) :: Full_Flux
  ! HLLEL variables
  real :: lambda_l, lambda_r, WeightL, WeightR, Diffusion
  real, dimension(MaxStrip) :: v_cleft_lf, v_cleft_hf, v_cright_rf, v_cright_hf


  ! Misc and temporary 
  real :: inv_c, B1dotB0, UdotB1, Vdt, fA, v_diffusion , v_coeffAW, coeff

  logical :: oktest, oktest_me, oktest_row
  !--------------------------------------------------------------------------

  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('calc_facefluxes',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if
  v_En_lf=cZero; v_En_rf=cZero                  !^CFG IF BORISCORR

  !\
  ! x-face fluxes --- x-face fluxes --- x-face fluxes
  !/
  iDir = 1
  fA = fAx_BLK(globalBLK)
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
     do k=1,nK
        do j=1,nJ
           iStart=1
           nStrip=1
           call get_fluxes
        end do
     end do
  else if ( neiLwest(globalBLK)==+1) then
     do k=1,nK
        do j=1,nJ
           iStart=nIFace
           nStrip = 1
           call get_fluxes
        end do
     end do
  end if

  ! y-face fluxes --- y-face fluxes --- y-face fluxes
  iDir = 2
  fA = fAy_BLK(globalBLK)
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
     do k=1,nK
        j=1
        do iStart=1,nI,MaxStrip
           nStrip = min(MaxStrip,nI-(iStart-1))
           call get_fluxes
        end do
     end do
  else if (neiLnorth(globalBLK)==+1) then
     do k=1,nK
        j = nJFace 
        do iStart=1,nI,MaxStrip
           nStrip = min(MaxStrip,nI-(iStart-1))
           call get_fluxes
        end do
     end do
  end if

  ! z-face fluxes --- z-face fluxes --- z-face fluxes 
  iDir = 3
  fA = fAz_BLK(globalBLK)
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
     k=1
     do j=1,nJ
        do iStart=1,nI,MaxStrip
           nStrip = min(MaxStrip,nI-(iStart-1))
           call get_fluxes
        end do
     end do
  else if ( neiLtop(globalBLK)==+1) then
     k=nKFace            
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
    !-------------------------
    oktest_row = oktest_me .and. (j==Jtest .or. (iDir==2 .and. j==Jtest+1)) &
         .and. (k==Ktest .or. (iDir==3 .and. k==Ktest+1))

    select case (iDir)
    case (1)
       !\
       ! x face
       !/
       do i=1,nStrip
          iStrip=i+iStart-1
          !\
          ! B0 on the face
          !/
          v_B0n(i)  = B0xFace_x_BLK(iStrip,j,k,globalBLK)
          v_B0t1(i) = B0yFace_x_BLK(iStrip,j,k,globalBLK)
          v_B0t2(i) = B0zFace_x_BLK(iStrip,j,k,globalBLK)

          !\
          ! Left face
          !/
          v_rho_lf(i)  =  LeftState_VX(rho_,iStrip,j,k)
          v_Un_lf(i)   =  LeftState_VX(Ux_ ,iStrip,j,k)
          v_Ut1_lf(i)  =  LeftState_VX(Uy_ ,iStrip,j,k)
          v_Ut2_lf(i)  =  LeftState_VX(Uz_ ,iStrip,j,k)
          v_B1n_lf(i)  =  LeftState_VX(Bx_ ,iStrip,j,k)
          v_B1t1_lf(i) =  LeftState_VX(By_ ,iStrip,j,k)
          v_B1t2_lf(i) =  LeftState_VX(Bz_ ,iStrip,j,k)
          v_p_lf(i)    =  LeftState_VX(P_  ,iStrip,j,k)

          !\
          ! Right face
          !/
          v_rho_rf(i)  =  RightState_VX(rho_,iStrip,j,k)
          v_Un_rf(i)   =  RightState_VX(Ux_ ,iStrip,j,k)
          v_Ut1_rf(i)  =  RightState_VX(Uy_ ,iStrip,j,k)
          v_Ut2_rf(i)  =  RightState_VX(Uz_ ,iStrip,j,k)
          v_B1n_rf(i)  =  RightState_VX(Bx_ ,iStrip,j,k)
          v_B1t1_rf(i) =  RightState_VX(By_ ,iStrip,j,k)
          v_B1t2_rf(i) =  RightState_VX(Bz_ ,iStrip,j,k)
          v_p_rf(i)    =  RightState_VX(P_  ,iStrip,j,k)




       end do
    case (2)
       !\
       ! y face
       !/
       do i=1,nStrip
          iStrip=i+iStart-1
          !\
          ! B0 on the face
          !/
          v_B0n(i)  = B0yFace_y_BLK(iStrip,j,k,globalBLK)
          v_B0t1(i) = B0zFace_y_BLK(iStrip,j,k,globalBLK)
          v_B0t2(i) = B0xFace_y_BLK(iStrip,j,k,globalBLK)

          !\
          ! Left face
          !/
          v_rho_lf(i)  =  LeftState_VY(rho_,iStrip,j,k)
          v_Un_lf(i)   =  LeftState_VY(Uy_ ,iStrip,j,k)
          v_Ut1_lf(i)  =  LeftState_VY(Uz_ ,iStrip,j,k)
          v_Ut2_lf(i)  =  LeftState_VY(Ux_ ,iStrip,j,k)
          v_B1n_lf(i)  =  LeftState_VY(By_ ,iStrip,j,k)
          v_B1t1_lf(i) =  LeftState_VY(Bz_ ,iStrip,j,k)
          v_B1t2_lf(i) =  LeftState_VY(Bx_ ,iStrip,j,k)
          v_p_lf(i)    =  LeftState_VY(P_  ,iStrip,j,k)

          !\
          ! Right face
          !/
          v_rho_rf(i)  =  RightState_VY(rho_,iStrip,j,k)
          v_Un_rf(i)   =  RightState_VY(Uy_ ,iStrip,j,k)
          v_Ut1_rf(i)  =  RightState_VY(Uz_ ,iStrip,j,k)
          v_Ut2_rf(i)  =  RightState_VY(Ux_ ,iStrip,j,k)
          v_B1n_rf(i)  =  RightState_VY(By_ ,iStrip,j,k)
          v_B1t1_rf(i) =  RightState_VY(Bz_ ,iStrip,j,k)
          v_B1t2_rf(i) =  RightState_VY(Bx_ ,iStrip,j,k)
          v_p_rf(i)    =  RightState_VY(P_  ,iStrip,j,k)
       end do
    case (3)
       !\
       ! z face
       !/
       do i=1,nStrip
          iStrip=i+iStart-1
          !\
          ! B0 on the face
          !/
          v_B0n(i)  = B0zFace_z_BLK(iStrip,j,k,globalBLK)
          v_B0t1(i) = B0xFace_z_BLK(iStrip,j,k,globalBLK)
          v_B0t2(i) = B0yFace_z_BLK(iStrip,j,k,globalBLK)

          !\
          ! Left face
          !/
          v_rho_lf(i)  =  LeftState_VZ(rho_,iStrip,j,k)
          v_Un_lf(i)   =  LeftState_VZ(Uz_ ,iStrip,j,k)
          v_Ut1_lf(i)  =  LeftState_VZ(Ux_ ,iStrip,j,k)
          v_Ut2_lf(i)  =  LeftState_VZ(Uy_ ,iStrip,j,k)
          v_B1n_lf(i)  =  LeftState_VZ(Bz_ ,iStrip,j,k)
          v_B1t1_lf(i) =  LeftState_VZ(Bx_ ,iStrip,j,k)
          v_B1t2_lf(i) =  LeftState_VZ(By_ ,iStrip,j,k)
          v_p_lf(i)    =  LeftState_VZ(P_  ,iStrip,j,k)

          !\
          ! Right face
          !/
          v_rho_rf(i)  =  RightState_VZ(rho_,iStrip,j,k)
          v_Un_rf(i)   =  RightState_VZ(Uz_ ,iStrip,j,k)
          v_Ut1_rf(i)  =  RightState_VZ(Ux_ ,iStrip,j,k)
          v_Ut2_rf(i)  =  RightState_VZ(Uy_ ,iStrip,j,k)
          v_B1n_rf(i)  =  RightState_VZ(Bz_ ,iStrip,j,k)
          v_B1t1_rf(i) =  RightState_VZ(Bx_ ,iStrip,j,k)
          v_B1t2_rf(i) =  RightState_VZ(By_ ,iStrip,j,k)
          v_p_rf(i)    =  RightState_VZ(P_ ,iStrip,j,k)
       end do
    end select

    !
    !jump in the normal component of B1 field is removed
    !
    do i=1,nStrip
       v_diffrB1n(i)=cHalf*(v_B1n_lf(i)-v_B1n_rf(i))	
       v_B1n_lf(i)=v_B1n_rf(i)+ v_diffrB1n(i)
       v_B1n_rf(i)=v_B1n_lf(i)
    end do
    !
    !
    !\
    ! Common variables used for physical fluxes
    !/
    do i=1,nStrip
       !\
       ! Left face
       !/
       ! Total magnetic field
       v_Bn_lf(i)  = v_B0n(i) +v_B1n_lf(i)
       v_Bt1_lf(i) = v_B0t1(i)+v_B1t1_lf(i)
       v_Bt2_lf(i) = v_B0t2(i)+v_B1t2_lf(i)

       ! Magnetic field squared
       v_BB_lf(i)  = v_Bn_lf(i)**2 + v_Bt1_lf(i)**2 + v_Bt2_lf(i)**2
       v_BB1_lf(i) = v_B1n_lf(i)**2 + v_B1t1_lf(i)**2 + v_B1t2_lf(i)**2

       ! Energy density
       UU_lf  = v_Un_lf(i)**2 + v_Ut1_lf(i)**2 + v_Ut2_lf(i)**2
       v_E_lf(i) = v_p_lf(i)*inv_gm1 + cHalf*v_rho_lf(i)*UU_lf + cHalf*v_BB1_lf(i)

       !\
       ! Right face
       !/
       ! Total magnetic field
       v_Bn_rf(i)  = v_B0n(i) +v_B1n_rf(i)
       v_Bt1_rf(i) = v_B0t1(i)+v_B1t1_rf(i)
       v_Bt2_rf(i) = v_B0t2(i)+v_B1t2_rf(i)

       ! Magnetic field squared
       v_BB_rf(i)  = v_Bn_rf(i)**2 + v_Bt1_rf(i)**2 + v_Bt2_rf(i)**2
       v_BB1_rf(i) = v_B1n_rf(i)**2 + v_B1t1_rf(i)**2 + v_B1t2_rf(i)**2

       ! Energy density
       UU_rf  = v_Un_rf(i)**2 + v_Ut1_rf(i)**2 + v_Ut2_rf(i)**2
       v_E_rf(i) = v_p_rf(i)*inv_gm1 + cHalf*v_rho_rf(i)*UU_rf + cHalf*v_BB1_rf(i)

    end do

    ! Average (hat) face primitive 
    do i=1,nStrip
       iStrip=i+iStart-1
       v_rho_hf(i) = cHalf*(v_rho_lf(i) + v_rho_rf(i))
       v_Un_hf(i)  = cHalf*(v_Un_lf(i)  + v_Un_rf(i))
       v_B1n_hf(i) = cHalf*(v_B1n_lf(i) + v_B1n_rf(i))
       v_Bn_hf(i)  = cHalf*(v_Bn_lf(i)  + v_Bn_rf(i))
       v_Bt1_hf(i) = cHalf*(v_Bt1_lf(i) + v_Bt1_rf(i))
       v_Bt2_hf(i) = cHalf*(v_Bt2_lf(i) + v_Bt2_rf(i))
       v_p_hf(i)   = cHalf*(v_p_lf(i)   + v_p_rf(i))
       v_BB_hf(i)  = v_Bn_hf(i)**2 + v_Bt1_hf(i)**2 + v_Bt2_hf(i)**2

       if(oktest_row.and.(iStrip==Itest .or. (iDir==1.and.iStrip==Itest+1)))then
          write(*,*)'Hat state for face=',iDir,' at I=',iStrip,' J=',j,' K=',k
          write(*,*)'rho=',v_rho_hf(i)
          write(*,*)'Un =',v_Un_hf(i)
          write(*,*)'P  =',v_p_hf(i)
          write(*,*)'B  =',v_Bn_hf(i),v_Bt1_hf(i),v_Bt2_hf(i)
          write(*,*)'BB =',v_BB_hf(i)
       end if

    end do


    if(boris_correction) then  !^CFG IF BORISCORR BEGIN
       call get_flux_boris
    else                       !^CFG END BORISCORR
       call get_flux_mhd
    end if                     !^CFG IF BORISCORR

    !\
    ! HLLEL interface fluxes
    !/
    do i=1,nStrip

       lambda_l = min(v_cleft_hf(i),  v_cleft_lf(i), cZero)
       lambda_r = max(v_cright_hf(i), v_cright_rf(i), cZero)
       WeightL = lambda_r/(lambda_r-lambda_l)
       WeightR = cOne-WeightL
       Diffusion=lambda_l*WeightL
       do n=1,nFlux
          Full_Flux(n,i) = fA*(WeightL*L_Flux(n,i) &
               +WeightR*R_Flux(n,i) + Diffusion*delta_U(n,i)) 
       end do
       Full_Flux(nFlux+1,i) = fA*(WeightL*v_Un_lf(i) + WeightR*v_Un_rf(i))
       Full_Flux(nFlux+2,i) = &                         !^CFG IF BORISCORR
            fA*(WeightL*v_En_lf(i) + WeightR*v_En_rf(i))!^CFG IF BORISCORR
    end do

    do i=1,nStrip
       coeff = fa*v_max_hf(i)*v_diffrB1n(i)
       Full_Flux(Bx_,i) = Full_Flux(Bx_,i) + coeff
       Full_Flux(Energy_ ,i) = Full_Flux(Energy_ ,i) + coeff*v_B1n_hf(i)
    end do


    if(oktest_row) then
       do i=1,nStrip
          iStrip=i+iStart-1

          if(iStrip==iTEST.or.(iDir==1.and.iStrip==iTEST+1))then
             write(*,*)'Fluxes for dir=',iDir,' at I=',iStrip,' J=',j,' K=',k
             write(*,*)'Eigenvalue_maxabs=',v_max_hf(i)
             do n=1,nFlux
                write(*,'(a,i2,4(1pe13.5))') 'Var,F,F_L,F_R,dU=',&
                     n,Full_Flux(n,i),L_Flux(n,i),R_Flux(n,i),delta_U(n,i)
             end do
          endif

       end do
    end if

    do i=1,nStrip
       iStrip=i+iStart-1

       !\
       ! Move interface fluxes to face fluxes for use elsewhere
       !   (also Vdt)
       !/
       select case (iDir)
       case (1)
          !\
          ! x face
          !/
          Flux_VX(rho_   ,iStrip,j,k)    = Full_Flux(1,i)
          Flux_VX(rhoUx_ ,iStrip,j,k)    = Full_Flux(2,i)
          Flux_VX(rhoUy_ ,iStrip,j,k)    = Full_Flux(3,i)
          Flux_VX(rhoUz_ ,iStrip,j,k)    = Full_Flux(4,i)
          Flux_VX(Bx_    ,iStrip,j,k)    = Full_Flux(5,i)
          Flux_VX(By_    ,iStrip,j,k)    = Full_Flux(6,i)
          Flux_VX(Bz_    ,iStrip,j,k)    = Full_Flux(7,i)
          Flux_VX(P_ ,iStrip,j,k)= Full_Flux(P_,i)
          Flux_VX(Energy_,iStrip,j,k)= Full_Flux(Energy_,i)

          UDotFA_x(iStrip,j,k)     = Full_Flux(Energy_+1,i)    
          EDotFA_x(iStrip,j,k)     = Full_Flux(Energy_+2,i)!^CFG IF BORISCORR
          ! Compute face time step (actually Volume/dt)
          VdtFace_x(iStrip,j,k)    = fA * v_max_hf(i)
       case (2)
          !\
          ! y face
          !/
          Flux_VY(rho_  ,iStrip,j,k) = Full_Flux(1,i)
          Flux_VY(rhoUx_,iStrip,j,k) = Full_Flux(4,i)
          Flux_VY(rhoUy_,iStrip,j,k) = Full_Flux(2,i)
          Flux_VY(rhoUz_,iStrip,j,k) = Full_Flux(3,i)
          Flux_VY(Bx_   ,iStrip,j,k) = Full_Flux(7,i)
          Flux_VY(By_   ,iStrip,j,k) = Full_Flux(5,i)
          Flux_VY(Bz_   ,iStrip,j,k) = Full_Flux(6,i)
          Flux_VY(P_ ,iStrip,j,k) = Full_Flux(P_,i)
          Flux_VY(Energy_,iStrip,j,k) = Full_Flux(Energy_,i)

          UDotFA_y(iStrip,j,k)     = Full_Flux(10,i)    
          EDotFA_y(iStrip,j,k)     = Full_Flux(11,i)!^CFG IF BORISCORR
          ! Compute face time step (actually Volume/dt)
          VdtFace_y(iStrip,j,k)    = fA * v_max_hf(i)
       case (3)
          !\
          ! z face
          !/
          Flux_VZ(rho_  ,iStrip,j,k) = Full_Flux(1,i)
          Flux_VZ(rhoUx_,iStrip,j,k) = Full_Flux(3,i)
          Flux_VZ(rhoUy_,iStrip,j,k) = Full_Flux(4,i)
          Flux_VZ(rhoUz_,iStrip,j,k) = Full_Flux(2,i)
          Flux_VZ(Bx_   ,iStrip,j,k) = Full_Flux(6,i)
          Flux_VZ(By_   ,iStrip,j,k) = Full_Flux(7,i)
          Flux_VZ(Bz_   ,iStrip,j,k) = Full_Flux(5,i)
          Flux_VZ(P_    ,iStrip,j,k) = Full_Flux(P_,i)
          Flux_VZ(Energy_    ,iStrip,j,k) = Full_Flux(Energy_,i)
          UDotFA_z(iStrip,j,k)     = Full_Flux(Energy_+1,i)    
          EDotFA_z(iStrip,j,k)     = Full_Flux(Energy_+2,i)!^CFG IF BORISCORR
          ! Compute face time step (actually Volume/dt)
          VdtFace_z(iStrip,j,k)    = fA * v_max_hf(i)
       end select
    end do

  end subroutine get_fluxes

  !==========================================================================
  !==========================================================================
  !==========================================================================

  subroutine get_flux_mhd
    !\
    ! Absolute value of fast wave speed for hat face
    !/
    do i=1,nStrip
       inv_rho  = cOne/v_rho_hf(i)
       a2       = g*v_p_hf(i)*inv_rho           ! Sound speed squared
       vA2      = v_BB_hf(i)*inv_rho            ! Maximum Alfven speed squared
       vAn2     = v_Bn_hf(i)**2*inv_rho         ! Normal  Alfven speed squared

       ! Approximate slow and fast wave speeds
       a2Va2 = a2 + vA2
       discr = max(cZero,a2Va2**2-cFour*a2*vAn2)
       discr = sqrt(discr)

       if(UseBorisSimple)then    !^CFG IF SIMPLEBORIS BEGIN
          ! Reduce fast speed with the "Alfven Lorentz" factor
          ! This is exact for u=0, and an upper estimate otherwise
          discr_fast = sqrt( cHalf*( a2Va2 + discr )/(cOne+vA2*inv_c2LIGHT) )
       else                      !^CFG END SIMPLEBORIS
          discr_fast = sqrt( cHalf*( a2Va2 + discr ) )
       end if                    !^CFG IF SIMPLEBORIS

       v_cright_hf(i)= v_Un_hf(i) + discr_fast
       v_cleft_hf(i) = v_Un_hf(i) - discr_fast

       v_max_hf(i)   = max(v_cright_hf(i),-v_cleft_hf(i))

    end do

    do i=1,nStrip

       !\
       ! Left going speed 
       !/
       inv_rho  = cOne/v_rho_lf(i)
       a2       = g*v_p_lf(i)*inv_rho        ! Sound speed squared
       vA2      = v_BB_lf(i)*inv_rho         ! Maximum Alfven speed squared
       vAn2     = v_Bn_lf(i)**2*inv_rho      ! Normal  Alfven speed squared

       ! Approximate slow and fast wave speeds
       a2Va2 = a2 + vA2
       discr = max(cZero,a2Va2**2-cFour*a2*vAn2)
       discr = sqrt(discr)
       if(UseBorisSimple)then    !^CFG IF SIMPLEBORIS BEGIN
          discr_fast = sqrt( cHalf*( a2Va2 + discr )/(cOne+vA2*inv_c2LIGHT) )
       else                      !^CFG END SIMPLEBORIS
          discr_fast = sqrt( cHalf*( a2Va2 + discr ) )
       end if                    !^CFG IF SIMPLEBORIS

       v_cleft_lf(i) = v_Un_lf(i) - discr_fast

       !\
       ! Right going speed 
       !/
       inv_rho  = cOne/v_rho_rf(i)
       a2       = g*v_p_rf(i)*inv_rho         ! Sound speed squared
       vA2      = v_BB_rf(i)*inv_rho          ! Maximum Alfven speed squared
       vAn2     = v_Bn_rf(i)**2*inv_rho       ! Normal  Alfven speed squared

       ! Approximate slow and fast wave speeds
       a2Va2 = a2 + vA2
       discr = max(cZero,a2Va2**2-cFour*a2*vAn2)
       discr = sqrt(discr)

       if(UseBorisSimple)then    !^CFG IF SIMPLEBORIS BEGIN
          ! Reduce fast speed with the "Alfven Lorentz" factor
          ! This is exact for u=0, and an upper estimate otherwise
          discr_fast = sqrt( cHalf*( a2Va2 + discr )/(cOne+vA2*inv_c2LIGHT) )
       else                      !^CFG END SIMPLEBORIS
          discr_fast = sqrt( cHalf*( a2Va2 + discr ) )
       end if                    !^CFG IF SIMPLEBORIS

       v_cright_rf(i)= v_Un_rf(i) + discr_fast

    end do


    !\
    ! Left state physical fluxes
    !/
    do i=1,nStrip

       B1dotB0 = &
            v_B1n_lf(i)*v_B0n(i) + v_B1t1_lf(i)*v_B0t1(i) + v_B1t2_lf(i)*v_B0t2(i)

       UdotB1 = v_B1n_lf(i)*v_Un_lf(i)+ v_B1t1_lf(i)*v_Ut1_lf(i)+ &
            v_B1t2_lf(i)*v_Ut2_lf(i)

       L_Flux(rho_,i) = v_rho_lf(i)*v_Un_lf(i)

       L_Flux(rhoUx_,i) = &
            v_rho_lf(i)*v_Un_lf(i)* v_Un_lf(i) + v_p_lf(i)         & !HD
            - v_B1n_lf(i)**2 + cHalf*v_BB1_lf(i)                       & !MHD B1
            + B1dotB0 - v_B0n(i)*v_B1n_lf(i) - v_B1n_lf(i)*v_B0n(i)    !Split

       L_Flux(rhoUy_,i) = &
            v_rho_lf(i)*v_Un_lf(i)*v_Ut1_lf(i)                     & !HD
            - v_B1n_lf(i)*v_B1t1_lf(i)                               & !MHD B1
            - v_B0n(i)*v_B1t1_lf(i) - v_B1n_lf(i)*v_B0t1(i)            !Split

       L_Flux(rhoUz_,i) = &
            v_rho_lf(i)*v_Un_lf(i)*v_Ut2_lf(i)                     & !HD
            - v_B1n_lf(i)*v_B1t2_lf(i)                               & !MHD B1
            - v_B0n(i)*v_B1t2_lf(i) - v_B1n_lf(i)*v_B0t2(i)            !Split

       L_Flux(Bx_,i) = cZero

       L_Flux(By_,i) = v_Un_lf(i)*v_Bt1_lf(i) - v_Ut1_lf(i)*v_Bn_lf(i)

       L_Flux(Bz_,i) = v_Un_lf(i)*v_Bt2_lf(i) - v_Ut2_lf(i)*v_Bn_lf(i)

       L_Flux(Energy_,i)  = &
            v_Un_lf(i)*(v_E_lf(i) + v_p_lf(i) + cHalf*v_BB1_lf(i) + B1dotB0) &
            - v_Bn_lf(i)*UdotB1

       L_Flux(P_,i) = v_p_lf(i)*v_Un_lf(i)
    end do

    !\
    ! Right state physical fluxes
    !/
    do i=1,nStrip

       B1dotB0 = &
            v_B1n_rf(i)*v_B0n(i) + v_B1t1_rf(i)*v_B0t1(i) + v_B1t2_rf(i)*v_B0t2(i)

       UdotB1 = v_B1n_rf(i)*v_Un_rf(i)+ v_B1t1_rf(i)*v_Ut1_rf(i)+ &
            v_B1t2_rf(i)*v_Ut2_rf(i)

       R_Flux(rho_,i) = v_rho_rf(i)*v_Un_rf(i)

       R_Flux(rhoUx_,i) = &
            v_rho_rf(i)*v_Un_rf(i)* v_Un_rf(i) + v_p_rf(i)         & !HD
            - v_B1n_rf(i)**2 + cHalf*v_BB1_rf(i)                       & !MHD B1
            + B1dotB0 - v_B0n(i)*v_B1n_rf(i) - v_B1n_rf(i)*v_B0n(i)    !Split

       R_Flux(rhoUy_,i) = &
            v_rho_rf(i)*v_Un_rf(i)*v_Ut1_rf(i)                     & !HD
            - v_B1n_rf(i)*v_B1t1_rf(i)                               & !MHD B1
            - v_B0n(i)*v_B1t1_rf(i) - v_B1n_rf(i)*v_B0t1(i)            !Split

       R_Flux(rhoUz_,i) = &
            v_rho_rf(i)*v_Un_rf(i)*v_Ut2_rf(i)                     & !HD
            - v_B1n_rf(i)*v_B1t2_rf(i)                               & !MHD B1
            - v_B0n(i)*v_B1t2_rf(i) - v_B1n_rf(i)*v_B0t2(i)            !Split

       R_Flux(Bx_,i) = cZero

       R_Flux(By_,i) = v_Un_rf(i)*v_Bt1_rf(i) - v_Ut1_rf(i)*v_Bn_rf(i)

       R_Flux(Bz_,i) = v_Un_rf(i)*v_Bt2_rf(i) - v_Ut2_rf(i)*v_Bn_rf(i)

       R_Flux(Energy_,i)  = &
            v_Un_rf(i)*(v_E_rf(i) + v_p_rf(i) + cHalf*v_BB1_rf(i) + B1dotB0) &
            - v_Bn_rf(i)*UdotB1

       R_Flux(P_,i) = v_p_rf(i)*v_Un_rf(i)
    end do

    !\
    ! Conserved state solution jumps
    !/
    do i=1,nStrip

       delta_U(rho_,i)    = v_rho_rf(i) -v_rho_lf(i)

       delta_U(rhoUx_,i)  = v_rho_rf(i)*v_Un_rf(i)  - v_rho_lf(i)*v_Un_lf(i)
       delta_U(rhoUy_,i)  = v_rho_rf(i)*v_Ut1_rf(i) - v_rho_lf(i)*v_Ut1_lf(i)
       delta_U(rhoUz_,i)  = v_rho_rf(i)*v_Ut2_rf(i) - v_rho_lf(i)*v_Ut2_lf(i)

       delta_U(Bx_,i) = v_B1n_rf(i)  -v_B1n_lf(i)
       delta_U(By_,i) = v_B1t1_rf(i) -v_B1t1_lf(i)
       delta_U(Bz_,i) = v_B1t2_rf(i) -v_B1t2_lf(i)

       delta_U(Energy_,i)  = v_E_rf(i) -v_E_lf(i)

       delta_U(P_,i)  = v_p_rf(i) -v_p_lf(i)
    end do

    if(UseBorisSimple)then  !^CFG IF SIMPLEBORIS BEGIN
       ! Correct the jump in the momentum using the (1+VA2/c^2)*rho*U
       do i=1,nStrip
          Gamma2R = 1+v_BB_rf(i)/v_rho_rf(i)*inv_c2LIGHT
          Gamma2L = 1+v_BB_lf(i)/v_rho_lf(i)*inv_c2LIGHT

          delta_U(rhoUx_,i)  = v_rho_rf(i)*v_Un_rf(i) *Gamma2R &
               - v_rho_lf(i)*v_Un_lf(i) *Gamma2L
          delta_U(rhoUy_,i)  = v_rho_rf(i)*v_Ut1_rf(i)*Gamma2R &
               - v_rho_lf(i)*v_Ut1_lf(i)*Gamma2L
          delta_U(rhoUz_,i)  = v_rho_rf(i)*v_Ut2_rf(i)*Gamma2R &
               - v_rho_lf(i)*v_Ut2_lf(i)*Gamma2L
       end do
    end if                  !^CFG END SIMPLEBORIS

  end subroutine get_flux_mhd

  !==========================================================================
  !==========================================================================
  !==========================================================================
  !^CFG IF BORISCORR BEGIN
  subroutine get_flux_boris
    !\
    ! Extra variables used for physical fluxes
    !/
    inv_c = cOne/cLIGHT

    do i=1,nStrip
       !\
       ! Left face
       !/
       ! Electric field divided by speed of light: E= - U x B / c
       v_En_lf(i)  = inv_c*(v_Bt1_lf(i)*v_Ut2_lf(i) - v_Bt2_lf(i)*v_Ut1_lf(i))
       v_Et1_lf(i) = inv_c*(v_Bt2_lf(i)*v_Un_lf(i)  - v_Bn_lf(i) *v_Ut2_lf(i))
       v_Et2_lf(i) = inv_c*(v_Bn_lf(i) *v_Ut1_lf(i) - v_Bt1_lf(i)*v_Un_lf(i) )

       ! Electric field squared/c^2 (with second and first order B0)
       v_EE_lf(i)  = v_En_lf(i)**2 + v_Et1_lf(i)**2 + v_Et2_lf(i)**2

       !\
       ! Right face
       !/
       ! Electric field divided by speed of light: E= - U x B / c = (1/c)B x U
       v_En_rf(i)  = inv_c*(v_Bt1_rf(i)*v_Ut2_rf(i) - v_Bt2_rf(i)*v_Ut1_rf(i))
       v_Et1_rf(i) = inv_c*(v_Bt2_rf(i)*v_Un_rf(i)  - v_Bn_rf(i) *v_Ut2_rf(i))
       v_Et2_rf(i) = inv_c*(v_Bn_rf(i) *v_Ut1_rf(i) - v_Bt1_rf(i)*v_Un_rf(i) )

       ! Electric field squared/c^2 (with second and first order B0)
       v_EE_rf(i)  = v_En_rf(i)**2 + v_Et1_rf(i)**2 + v_Et2_rf(i)**2

    end do

    !\
    ! Absolute value of fast wave speed for hat face
    !/
    do i=1,nStrip
       inv_rho  = cOne/v_rho_hf(i)
       a2       = g*v_p_hf(i)*inv_rho           ! Sound speed squared
       vA2      = v_BB_hf(i)*inv_rho            ! Maximum Alfven speed squared
       vAn2     = v_Bn_hf(i)**2*inv_rho         ! Normal  Alfven speed squared

       gammaA2  = cOne/(cOne+vA2*inv_c2LIGHT)       ! "Alfven Lorentz" factor
       gammaU2  = max(cZero,cOne-gammaA2*v_Un_hf(i)**2*inv_c2LIGHT) !1-gA^2*Un^2/c^2

       ! Modified speeds
       a2Boris  = a2  *gammaA2*(1+vAn2*inv_c2LIGHT)
       vA2Boris = vA2 *gammaA2*gammaU2
       vAn2Boris= vAn2*gammaA2*gammaU2

       ! Approximate slow and fast wave speeds
       a2Va2Boris = a2Boris + vA2Boris
       discr      = max(cZero,a2Va2Boris**2-cFour*a2*vAn2Boris)
       discr      = sqrt(discr)

       discr_fast = sqrt( cHalf*( a2Va2Boris + discr ) )
       discr_slow = sqrt( cHalf*( max(cZero, a2Va2Boris - discr ) ) )

       ! In extreme cases "slow" wave can be faster than "fast" wave
       ! so take maximum of the two
       v_cright_hf(i)= max(v_Un_hf(i)*gammaA2 + discr_fast, v_Un_hf(i) + discr_slow)
       v_cleft_hf(i) = min(v_Un_hf(i)*gammaA2 - discr_fast, v_Un_hf(i) - discr_slow)

       v_max_hf(i)   = max(abs(v_cright_hf(i)), abs(v_cleft_hf(i)))

    end do

    do i=1,nStrip

       !\
       ! Left going speed 
       !/
       inv_rho  = cOne/v_rho_lf(i)
       a2       = g*v_p_lf(i)*inv_rho           ! Sound speed squared
       vA2      = v_BB_lf(i)*inv_rho            ! Maximum Alfven speed squared
       vAn2     = v_Bn_lf(i)**2*inv_rho         ! Normal  Alfven speed squared

       gammaA2  = cOne/(cOne+vA2*inv_c2LIGHT)       ! "Alfven Lorentz" factor
       gammaU2  = max(cZero,cOne-gammaA2*v_Un_lf(i)**2*inv_c2LIGHT) !1-gA^2*Un^2/c^2

       ! Modified speeds
       a2Boris  = a2  *gammaA2*(1+vAn2*inv_c2LIGHT)
       vA2Boris = vA2 *gammaA2*gammaU2
       vAn2Boris= vAn2*gammaA2*gammaU2

       ! Approximate slow and fast wave speeds
       a2Va2Boris = a2Boris + vA2Boris
       discr      = max(cZero,a2Va2Boris**2-cFour*a2*vAn2Boris)
       discr      = sqrt(discr)

       discr_fast = sqrt( cHalf*( a2Va2Boris + discr ) )
       discr_slow = sqrt( cHalf*(max(cZero, a2Va2Boris - discr ) ) )

       ! In extreme cases "slow" wave can be faster than "fast" wave
       ! so take maximum of the two
       v_cleft_lf(i) = min(v_Un_lf(i)*gammaA2 - discr_fast, v_Un_lf(i) - discr_slow)

       !\
       ! Right going speed 
       !/
       inv_rho  = cOne/v_rho_rf(i)
       a2       = g*v_p_rf(i)*inv_rho           ! Sound speed squared
       vA2      = v_BB_rf(i)*inv_rho            ! Maximum Alfven speed squared
       vAn2     = v_Bn_rf(i)**2*inv_rho         ! Normal  Alfven speed squared

       gammaA2  = cOne/(cOne+vA2*inv_c2LIGHT)       ! "Alfven Lorentz" factor
       gammaU2  = max(cZero,cOne-gammaA2*v_Un_rf(i)**2*inv_c2LIGHT) !1-gA^2*Un^2/c^2

       ! Modified speeds
       a2Boris  = a2  *gammaA2*(1+vAn2*inv_c2LIGHT)
       vA2Boris = vA2 *gammaA2*gammaU2
       vAn2Boris= vAn2*gammaA2*gammaU2

       ! Approximate slow and fast wave speeds
       a2Va2Boris = a2Boris + vA2Boris
       discr      = max(cZero,a2Va2Boris**2-cFour*a2*vAn2Boris)
       discr      = sqrt(discr)

       discr_fast = sqrt( cHalf*( a2Va2Boris + discr ) )
       discr_slow = sqrt( cHalf*(max(cZero, a2Va2Boris - discr ) ) )

       ! In extreme cases "slow" wave can be faster than "fast" wave
       ! so take maximum of the two
       v_cright_rf(i)= max(v_Un_rf(i)*gammaA2 + discr_fast, v_Un_rf(i) + discr_slow)

    end do

    !\
    ! Left state physical fluxes
    !/
    do i=1,nStrip

       B1dotB0 = &
            v_B1n_lf(i)*v_B0n(i) + v_B1t1_lf(i)*v_B0t1(i) + v_B1t2_lf(i)*v_B0t2(i)

       UdotB1 = v_B1n_lf(i)*v_Un_lf(i)+ v_B1t1_lf(i)*v_Ut1_lf(i)+ &
            v_B1t2_lf(i)*v_Ut2_lf(i)

       L_Flux(rho_,i) = v_rho_lf(i)*v_Un_lf(i)

       L_Flux(rhoUx_,i) = &
            v_rho_lf(i)*v_Un_lf(i)* v_Un_lf(i) + v_p_lf(i)         & !HD
            - v_B1n_lf(i)**2 + cHalf*v_BB1_lf(i)                       & !MHD B1
            + B1dotB0 - v_B0n(i)*v_B1n_lf(i) - v_B1n_lf(i)*v_B0n(i)  & !Split
            - v_En_lf(i)**2  + cHalf*v_EE_lf(i)                          !Boris

       L_Flux(rhoUy_,i) = &
            v_rho_lf(i)*v_Un_lf(i)*v_Ut1_lf(i)                     & !HD
            - v_B1n_lf(i)*v_B1t1_lf(i)                               & !MHD B1
            - v_B0n(i)*v_B1t1_lf(i) - v_B1n_lf(i)*v_B0t1(i)          & !Split
            - v_En_lf(i)*v_Et1_lf(i)                                   !Boris

       L_Flux(rhoUz_,i) = &
            v_rho_lf(i)*v_Un_lf(i)*v_Ut2_lf(i)                     & !HD
            - v_B1n_lf(i)*v_B1t2_lf(i)                               & !MHD B1
            - v_B0n(i)*v_B1t2_lf(i) - v_B1n_lf(i)*v_B0t2(i)          & !Split
            - v_En_lf(i)*v_Et2_lf(i)                                   !Boris

       L_Flux(Bx_,i) = cZero

       L_Flux(By_,i) = v_Un_lf(i)*v_Bt1_lf(i) - v_Ut1_lf(i)*v_Bn_lf(i)

       L_Flux(Bz_,i) = v_Un_lf(i)*v_Bt2_lf(i) - v_Ut2_lf(i)*v_Bn_lf(i)

       L_Flux(Energy_,i)  = &
            v_Un_lf(i)*(v_E_lf(i) + v_p_lf(i) + cHalf*v_BB1_lf(i) + B1dotB0) &
            - v_Bn_lf(i)*UdotB1

       L_Flux(P_,i) = v_p_lf(i)*v_Un_lf(i)
    end do

    !\
    ! Right state physical fluxes
    !/
    do i=1,nStrip

       B1dotB0 = &
            v_B1n_rf(i)*v_B0n(i) + v_B1t1_rf(i)*v_B0t1(i) + v_B1t2_rf(i)*v_B0t2(i)

       UdotB1 = v_B1n_rf(i)*v_Un_rf(i)+ v_B1t1_rf(i)*v_Ut1_rf(i)+ &
            v_B1t2_rf(i)*v_Ut2_rf(i)

       R_Flux(rho_,i) = v_rho_rf(i)*v_Un_rf(i)

       R_Flux(rhoUx_,i) = &
            v_rho_rf(i)*v_Un_rf(i)* v_Un_rf(i) + v_p_rf(i)         & !HD
            - v_B1n_rf(i)**2 + cHalf*v_BB1_rf(i)                       & !MHD B1
            + B1dotB0 - v_B0n(i)*v_B1n_rf(i) - v_B1n_rf(i)*v_B0n(i)  & !Split
            - v_En_rf(i)**2  + cHalf*v_EE_rf(i)                          !Boris

       R_Flux(rhoUy_,i) = &
            v_rho_rf(i)*v_Un_rf(i)*v_Ut1_rf(i)                     & !HD
            - v_B1n_rf(i)*v_B1t1_rf(i)                               & !MHD B1
            - v_B0n(i)*v_B1t1_rf(i) - v_B1n_rf(i)*v_B0t1(i)          & !Split
            - v_En_rf(i)*v_Et1_rf(i)                                   !Boris

       R_Flux(rhoUz_,i) = &
            v_rho_rf(i)*v_Un_rf(i)*v_Ut2_rf(i)                     & !HD
            - v_B1n_rf(i)*v_B1t2_rf(i)                               & !MHD B1
            - v_B0n(i)*v_B1t2_rf(i) - v_B1n_rf(i)*v_B0t2(i)          & !Split
            - v_En_rf(i)*v_Et2_rf(i)                                   !Boris

       R_Flux(Bx_,i) = cZero

       R_Flux(By_,i) = v_Un_rf(i)*v_Bt1_rf(i) - v_Ut1_rf(i)*v_Bn_rf(i)

       R_Flux(Bz_,i) = v_Un_rf(i)*v_Bt2_rf(i) - v_Ut2_rf(i)*v_Bn_rf(i)

       R_Flux(Energy_,i)  = &
            v_Un_rf(i)*(v_E_rf(i) + v_p_rf(i) + cHalf*v_BB1_rf(i) + B1dotB0) &
            - v_Bn_rf(i)*UdotB1

       R_Flux(P_,i) = v_p_rf(i)*v_Un_rf(i)
    end do

    !\
    ! Conserved state solution jumps
    !/
    do i=1,nStrip

       UdotB_lf=v_Un_lf(i)*v_Bn_lf(i) + v_Ut1_lf(i)*v_Bt1_lf(i) &
            + v_Ut2_lf(i)*v_Bt2_lf(i)
       UdotB_rf=v_Un_rf(i)*v_Bn_rf(i) + v_Ut1_rf(i)*v_Bt1_rf(i) &
            + v_Ut2_rf(i)*v_Bt2_rf(i)

       delta_U(rho_,i)    = v_rho_rf(i) -v_rho_lf(i)

       ! rhoU_Boris = rhoU - ((U x B) x B)/c^2 = rhoU + (U B^2 - B U.B)/c^2
       delta_U(rhoUx_,i)  = v_rho_rf(i)*v_Un_rf(i) - v_rho_lf(i)*v_Un_lf(i) &
            +((v_Un_rf(i) *v_BB_rf(i) - v_Bn_rf(i) *UdotB_rf)    &
            -(v_Un_lf(i) *v_BB_lf(i) - v_Bn_lf(i) *UdotB_lf))*inv_c2LIGHT

       delta_U(rhoUy_,i)  = v_rho_rf(i)*v_Ut1_rf(i) - v_rho_lf(i)*v_Ut1_lf(i) &
            +((v_Ut1_rf(i)*v_BB_rf(i) - v_Bt1_rf(i)*UdotB_rf)    &
            -(v_Ut1_lf(i)*v_BB_lf(i) - v_Bt1_lf(i)*UdotB_lf))*inv_c2LIGHT

       delta_U(rhoUz_,i)  = v_rho_rf(i)*v_Ut2_rf(i) - v_rho_lf(i)*v_Ut2_lf(i) &
            +((v_Ut2_rf(i)*v_BB_rf(i) - v_Bt2_rf(i)*UdotB_rf)    &
            -(v_Ut2_lf(i)*v_BB_lf(i) - v_Bt2_lf(i)*UdotB_lf))*inv_c2LIGHT

       delta_U(Bx_,i) = v_B1n_rf(i)  -v_B1n_lf(i)
       delta_U(By_,i) = v_B1t1_rf(i) -v_B1t1_lf(i)
       delta_U(Bz_,i) = v_B1t2_rf(i) -v_B1t2_lf(i)

       delta_U(Energy_,i)  = v_E_rf(i)    -v_E_lf(i) + cHalf*(v_EE_rf(i)-v_EE_lf(i))

       delta_U(P_,i) = v_p_rf(i) - v_p_lf(i)
    end do

  end subroutine get_flux_boris
  !^CFG END BORISCORR
end subroutine calc_flux_Linde
