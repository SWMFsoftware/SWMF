!^CFG COPYRIGHT UM
!^CFG FILE NOT CARTESIAN

!^CFG FILE LINDEFLUX 
subroutine option_lindeflux(TrueOption,NameOption)

  logical, intent(out) :: TrueOption
  character (len=40), intent(out) :: NameOption

  TrueOption  = .true.
  NameOption  = 'COVARIANT LINDE FLUX  Sokolov'

end subroutine option_lindeflux
!=============================================================================

subroutine calc_flux_Linde(DoResChangeOnly)
  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK
  use ModCovariant

  use ModNumConst
  use ModAdvance
  use ModParallel, ONLY : neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModFlux

  implicit none
  logical, intent (in) :: DoResChangeOnly

  integer :: i,j,k
  integer :: iBegin,iEnd,jBegin,jEnd,kBegin,kEnd
  integer :: nStrip,iStart,iStrip
  logical :: oktest,oktest_me,oktest_row



  if(iProc==PROCtest .and. globalBLK==BLKtest)then
     call set_oktest('calc_facefluxes',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if


  !--------------------------------------------------------------------------
  !\
  ! x-face fluxes --- x-face fluxes --- x-face fluxes
  !/

  if (.not.DoResChangeOnly) then
     kBegin=kMinFaceX ; kEnd=kMaxFaceX
     jBegin=jMinFaceX ; jEnd=jMaxFaceX
     iBegin=1 ;             iEnd=nIFace
  else
     kBegin=1 ;        kEnd=nK
     jBegin=1 ;        jEnd=nJ
     iBegin=1 ;        iEnd=nIFace
     if (neiLeast(globalBLK)/=+1) iBegin=nIFace
     if (neiLwest(globalBLK)/=+1) iEnd=1
  end if


  !\
  ! x-face fluxes --- x-face fluxes --- x-face fluxes
  !/
  !
  FaceArea2Min=FaceArea2MinI_B(globalBLK)
  do iStart=iBegin,iEnd,MaxStrip 
     nStrip=min(MaxStrip,1+iEnd-iStart) 
     do k=kBegin,kEnd
        oktest_row=oktest_me.and.k==Ktest ! for test output
        do j=jBegin,jEnd
           oktest_row=oktest_row.and.j==Jtest !for test output
           do i=1,nStrip
              iStrip=i+iStart-1		  

              call calc_faceareaI(iStrip,j,k,globalBLK,FaceArea_DI(:,i))
             
              ! GET PRIMITIVES
              !\
              ! B0 on the face
              !/
              v_B0x(i) = B0xFace_x_BLK(iStrip,j,k,globalBLK)
              v_B0y(i) = B0yFace_x_BLK(iStrip,j,k,globalBLK)
              v_B0z(i) = B0zFace_x_BLK(iStrip,j,k,globalBLK)
              !\
              ! Left face
              !/
              v_rho_lf(i) = LeftState_VX(rho_,iStrip,j,k)
              v_rhoSp_lf(1:MaxSpecies,i)=&
                   LeftState_VX(rho_+1:rho_+MaxSpecies,iStrip,j,k)
              v_rho_lf(i)=sum( v_rhoSp_lf(1:MaxSpecies,i))

              v_Ux_lf(i)  = LeftState_VX(Ux_ ,iStrip,j,k)
              v_Uy_lf(i)  = LeftState_VX(Uy_ ,iStrip,j,k)
              v_Uz_lf(i)  = LeftState_VX(Uz_ ,iStrip,j,k)
              v_B1x_lf(i) = LeftState_VX(Bx_ ,iStrip,j,k)
              v_B1y_lf(i) = LeftState_VX(By_ ,iStrip,j,k)
              v_B1z_lf(i) = LeftState_VX(Bz_ ,iStrip,j,k)
              v_p_lf(i)   = LeftState_VX(P_  ,iStrip,j,k)
              !\
              ! Right face
              !/
              v_rho_rf(i) = RightState_VX(rho_,iStrip,j,k)
              v_rhoSp_rf(1:MaxSpecies,i)=&
                   RightState_VX(rho_+1:rho_+MaxSpecies,iStrip,j,k)
              v_rho_rf(i)=sum( v_rhoSp_rf(1:MaxSpecies,i))

              v_Ux_rf(i)  = RightState_VX(Ux_ ,iStrip,j,k)
              v_Uy_rf(i)  = RightState_VX(Uy_ ,iStrip,j,k)
              v_Uz_rf(i)  = RightState_VX(Uz_ ,iStrip,j,k)
              v_B1x_rf(i) = RightState_VX(Bx_ ,iStrip,j,k) 
              v_B1y_rf(i) = RightState_VX(By_ ,iStrip,j,k) 
              v_B1z_rf(i) = RightState_VX(Bz_ ,iStrip,j,k) 
              v_p_rf(i)   = RightState_VX(P_  ,iStrip,j,k) 

              !\
              ! Normal components
              !/

           end do
           !\
           ! get fluxes
           !/

           call get_flux_mhdLinde(nStrip)	
           !
           !\
           ! test output
           !/	
           if(oktest_row) then
              do i=iStart,iStart+nStrip-1
                 if(i==Itest.or.i==Itest+1) then
                    write(*,*)'Hat var-s and Fluxes for iDir=1 at I=',i, &
                         ' J=',j,' K=',k  
                    call test_fluxes_wri
                 end if
              end do
           end if
           !\
           ! save fluxes
           !/
           !\
           ! x face
           !/
           do i=1,nStrip
              iStrip=i+iStart-1
              Flux_VX(1:Energy_,iStrip,j,k) = &
                   Full_Flux(1:Energy_,i)
              UDotFA_x(iStrip,j,k)          = Full_Flux(Energy_+2, i)        

              ! Compute face time step (actually Volume/dt)

              VdtFace_x(iStrip,j,k)    = v_max_hf(i)
           end do
        end do
     end do
  end do



  ! y-face fluxes --- y-face fluxes --- y-face fluxes

  if (.not.DoResChangeOnly) then
     kBegin = kMinFaceY ; kEnd = kMaxFaceY
     jBegin = 1 ;             jEnd = nJFace
     iBegin = iMinFaceY ; iEnd = iMaxFaceY
  else
     kBegin = 1 ; kEnd = nK
     jBegin = 1 ; jEnd = nJFace
     iBegin = 1 ; iEnd = nI
     if(neiLsouth(globalBLK)/=+1) jBegin=nJFace
     if (neiLnorth(globalBLK)/=+1) jEnd=1
  end if

  !_________________________________________________________________________
  ! y-face fluxes --- y-face fluxes --- y-face fluxes
  FaceArea2Min=FaceArea2MinJ_B(globalBLK)
  do iStart=iBegin,iEnd,MaxStrip
     nStrip=min(MaxStrip,1+iEnd-iStart) 
     do k=kBegin,kEnd
        oktest_row=oktest_me.and.k==Ktest ! for test output
        do j=jBegin,jEnd
           oktest_row=oktest_row.and.(j==Jtest.or.j==Jtest+1) !for test output

           !for non-cartesian the faceArea vector should be defined here..
           !\
           ! y face
           !/
           !\
           ! B0 on the face
           !/
           do i=1,nStrip
              iStrip=i+iStart-1

              call calc_faceareaJ(iStrip,j,k,globalBLK,FaceArea_DI(:,i))

              v_B0x(i) = B0xFace_y_BLK(iStrip,j,k,globalBLK)
              v_B0y(i) = B0yFace_y_BLK(iStrip,j,k,globalBLK)
              v_B0z(i) = B0zFace_y_BLK(iStrip,j,k,globalBLK)

              !\
              ! Left face
              !/
              v_rho_lf(i) = LeftState_VY(rho_,iStrip,j,k)
              v_rhoSp_lf(1:MaxSpecies,i)=&
                   LeftState_VY(rho_+1:rho_+MaxSpecies,iStrip,j,k)
              v_rho_lf(i)=sum( v_rhoSp_lf(1:MaxSpecies,i))
              
              v_Ux_lf(i)  = LeftState_VY(Ux_ ,iStrip,j,k)
              v_Uy_lf(i)  = LeftState_VY(Uy_ ,iStrip,j,k)
              v_Uz_lf(i)  = LeftState_VY(Uz_ ,iStrip,j,k)
              v_B1x_lf(i) = LeftState_VY(Bx_ ,iStrip,j,k)
              v_B1y_lf(i) = LeftState_VY(By_ ,iStrip,j,k)
              v_B1z_lf(i) = LeftState_VY(Bz_ ,iStrip,j,k)
              v_p_lf(i)   = LeftState_VY(P_  ,iStrip,j,k)
              !\
              ! Right face
              !/
              v_rho_rf(i) = RightState_VY(rho_,iStrip,j,k)
              v_rhoSp_rf(1:MaxSpecies,i)=&
                   RightState_VY(rho_+1:rho_+MaxSpecies,iStrip,j,k)
              v_rho_rf(i)=sum( v_rhoSp_rf(1:MaxSpecies,i))

              v_Ux_rf(i)  = RightState_VY(Ux_ ,iStrip,j,k)
              v_Uy_rf(i)  = RightState_VY(Uy_ ,iStrip,j,k)
              v_Uz_rf(i)  = RightState_VY(Uz_ ,iStrip,j,k)
              v_B1x_rf(i) = RightState_VY(Bx_ ,iStrip,j,k)
              v_B1y_rf(i) = RightState_VY(By_ ,iStrip,j,k) 
              v_B1z_rf(i) = RightState_VY(Bz_ ,iStrip,j,k) 
              v_p_rf(i)   = RightState_VY(P_  ,iStrip,j,k) 
           end do
           !\
           ! get fluxes
           !/

           call get_flux_mhdLinde(nStrip)
           !\
           ! test output
           !/	
           if(oktest_row) then
              do i=iStart,iStart+nStrip-1
                 if(i==Itest) then
                    write(*,*)'Hat var-s and Fluxes for iDir=2 at I=',&
                         i,' J=',j,' K=',k  
                    call test_fluxes_wri
                 end if
              end do
           end if
           !\
           ! save fluxes
           !/
           do i=1,nStrip
              iStrip=i+iStart-1
              Flux_VY(1:Energy_,iStrip,j,k) = &
                   Full_Flux(1:Energy_,i) 
              UDotFA_y(iStrip,j,k)          = Full_Flux(Energy_+2, i)

              ! Compute face time step (actually Volume/dt)

              VdtFace_y(iStrip,j,k) = v_max_hf(i)
           end do
        end do
     end do
  end do

  ! z-face fluxes --- z-face fluxes --- z-face fluxes 

  if (.not.DoResChangeOnly) then
     kBegin = 1 ;             kEnd = nKFace
     jBegin = jMinFaceZ ; jEnd = jMaxFaceZ
     iBegin = iMinFaceZ ; iEnd = iMaxFaceZ
  else
     kBegin = 1 ; kEnd = nKFace
     jBegin = 1 ; jEnd = nJ
     iBegin = 1 ; iEnd = nI
     if (neiLbot(globalBLK)/=+1) kBegin = nKFace
     if ( neiLtop(globalBLK)/=+1) kEnd = 1
  end if



  !_____________________________________________________________________       
  ! z-face fluxes --- z-face fluxes --- z-face fluxes 

  FaceArea2Min=FaceArea2MinK_B(globalBLK)
  do iStart=iBegin,iEnd,MaxStrip
     nStrip=min(MaxStrip,1+iEnd-iStart) 
     do k=kBegin,kEnd
        oktest_row=oktest_me.and.(k==Ktest.or.k==Ktest+1) ! for test output
        do j=jBegin,jEnd
           oktest_row=oktest_row.and.j==Jtest !for test output

           !\
           ! z face
           !/
           do i=1,nStrip
              iStrip=i+iStart-1

              call calc_faceareaK(iStrip,j,k,globalBLK,FaceArea_DI(:,i))
              !\
              ! B0 on the face
              !/
              v_B0x(i)  = B0xFace_z_BLK(iStrip,j,k,globalBLK)
              v_B0y(i)  = B0yFace_z_BLK(iStrip,j,k,globalBLK)
              v_B0z(i)  = B0zFace_z_BLK(iStrip,j,k,globalBLK)
              !\
              ! Left face
              !/
              v_rho_lf(i) = LeftState_VZ(rho_,iStrip,j,k)
              v_rhoSp_lf(1:MaxSpecies,i)=&
                   LeftState_VZ(rho_+1:rho_+MaxSpecies,iStrip,j,k)
              v_rho_lf(i)=sum( v_rhoSp_lf(1:MaxSpecies,i))
              
              v_Ux_lf(i)  = LeftState_VZ(Ux_ ,iStrip,j,k)
              v_Uy_lf(i)  = LeftState_VZ(Uy_ ,iStrip,j,k)
              v_Uz_lf(i)  = LeftState_VZ(Uz_ ,iStrip,j,k)
              v_B1x_lf(i) = LeftState_VZ(Bx_ ,iStrip,j,k)
              v_B1y_lf(i) = LeftState_VZ(By_ ,iStrip,j,k)
              v_B1z_lf(i) = LeftState_VZ(Bz_ ,iStrip,j,k)
              v_p_lf(i)   = LeftState_VZ(P_  ,iStrip,j,k)
              !\
              ! Right face
              !/
              v_rho_rf(i) = RightState_VZ(rho_,iStrip,j,k)
              v_rhoSp_rf(1:MaxSpecies,i)=&
                   RightState_VZ(rho_+1:rho_+MaxSpecies,iStrip,j,k) 
              v_rho_rf(i)=sum( v_rhoSp_rf(1:MaxSpecies,i))
             
              v_Ux_rf(i)  = RightState_VZ(Ux_ ,iStrip,j,k) 
              v_Uy_rf(i)  = RightState_VZ(Uy_ ,iStrip,j,k) 
              v_Uz_rf(i)  = RightState_VZ(Uz_ ,iStrip,j,k) 
              v_B1x_rf(i) = RightState_VZ(Bx_ ,iStrip,j,k) 
              v_B1y_rf(i) = RightState_VZ(By_ ,iStrip,j,k) 
              v_B1z_rf(i) = RightState_VZ(Bz_ ,iStrip,j,k)
              v_p_rf(i)   = RightState_VZ(P_  ,iStrip,j,k)  
           end do
           !\
           ! get fluxes
           !/
           call get_flux_mhdLinde(nStrip)	

           if(oktest_row) then
              do i=iStart,iStart+nStrip-1
                 if(i==Itest) then
                    write(*,*)'Hat var-s and Fluxes for iDir=3 at I=',&
                         i,' J=',j,' K=',k  
                    call test_fluxes_wri
                 end if
              end do
           end if
           !\
           ! save fluxes
           !/
           !\
           ! z face
           !/
           do i=1,nStrip
              iStrip=i+iStart-1

              Flux_VZ(1:Energy_,iStrip,j,k) = &
                   Full_Flux(1:Energy_,i)  
              UDotFA_z(iStrip,j,k)          = Full_Flux(Energy_+2, i)

              ! Compute face time step (actually Volume/dt)

              VdtFace_z(iStrip,j,k)    = v_max_hf(i)
           end do
        end do
     end do
  end do
Contains
  !_____________________________________________________________________
  subroutine test_fluxes_wri
    integer :: n    
    real    :: v_rho_hf, v_p_hf
    !_______________________________________________________________________

    iStrip=i-iStart+1
    v_rho_hf   = cHalf*(v_rho_lf(iStrip) + v_rho_rf(iStrip))
    v_p_hf     = cHalf*(v_p_lf(iStrip)   + v_p_rf(iStrip))
    write(*,*)'rho=',v_rho_hf
    write(*,*)'P  =', v_p_hf
    write(*,*)'Eigenvalue_maxabs*fA=',v_max_hf(iStrip)
    do n=1,nFlux
       write(*,'(a,i2,4(1pe13.5))') 'Var,F,F_L,F_R,dU=',&
            n,Full_Flux(n,iStrip)
    end do

  end subroutine test_fluxes_wri

end subroutine calc_flux_Linde

!==========================================================================
subroutine get_flux_mhdLinde(nStrip)
  use ModFlux
  use ModNumConst
  use ModMain, ONLY :x_,y_,z_
  use ModVarIndexes,ONLY:rho_, rhoUx_, rhoUy_, rhoUz_, Bx_, By_,&
       Bz_, Energy_,P_, MaxSpecies
  use ModPhysics, ONLY : g_half, inv_gm1, g
  implicit none
  integer, intent(in) :: nStrip                   ! the length of "strip", of flux
  real, dimension(MaxStrip) ::Bn_lf,Bn_rf
  real, dimension(MaxStrip) :: a2,Misc1_I,discr, mP_lf, mP_rf,B0B0,BFull2_I
  real, dimension(MaxStrip) :: v_cright_hf, v_cleft_hf
  real, dimension(nFlux,MaxStrip):: Fl_lf, Fl_rf
  real, dimension(MaxStrip)::C_E_lf ,C_Bx_lf,C_By_lf,C_Bz_lf,C_P_lf
  real, dimension(MaxStrip)::C_E_rf ,C_Bx_rf,C_By_rf,C_Bz_rf,C_P_rf

!----variables for HLLEL flux ---------------------------------------
  real, dimension(nFlux,MaxStrip):: U_l, U_r,dFlux,  delta_U
  real, dimension(MaxStrip)::v_cstar_hf, v_Un_hf, v_diff_hf,min_v, &
       lambda_l, lambda_r
  real :: Alpha

  integer :: i,n

  real,dimension(MaxStrip):: FaceArea2_I,GHalfFA2_I
  real,dimension(x_:z_,MaxStrip) :: DiffrB1n_DI
  !__________________________________________________________________________


  do i=1,nStrip
     FaceArea2_I(i) = max(FaceArea_DI(x_,i)**2+FaceArea_DI(y_,i)**2+&
          FaceArea_DI(z_,i)**2,FaceArea2Min)
     v_Un_lf(i)=FaceArea_DI(x_,i)*v_Ux_lf(i)+&
          FaceArea_DI(y_,i)*v_Uy_lf(i)+FaceArea_DI(z_,i)*v_Uz_lf(i)
     v_Un_rf(i)=FaceArea_DI(x_,i)*v_Ux_rf(i)+&
          FaceArea_DI(y_,i)*v_Uy_rf(i)+FaceArea_DI(z_,i)*v_Uz_rf(i)
     B0n(i)=FaceArea_DI(x_,i)*v_B0x(i)+&
          FaceArea_DI(y_,i)*v_B0y(i)+ FaceArea_DI(z_,i)*v_B0z(i)
     B1n_lf(i)=FaceArea_DI(x_,i)*v_B1x_lf(i)+&
          FaceArea_DI(y_,i)*v_B1y_lf(i)+FaceArea_DI(z_,i)*v_B1z_lf(i)
     B1n_rf(i)=FaceArea_DI(x_,i)*v_B1x_rf(i)+&
          FaceArea_DI(y_,i)*v_B1y_rf(i)+FaceArea_DI(z_,i)*v_B1z_rf(i)
  end do

  do i=1,nStrip
     Misc1_I(i)=cHalf*(B1n_lf(i)-B1n_rf(i))/FaceArea2_I(i)
     DiffrB1n_DI(x_,i) = FaceArea_DI(x_,i)*Misc1_I(i)
     DiffrB1n_DI(y_,i) = FaceArea_DI(y_,i)*Misc1_I(i)
     DiffrB1n_DI(z_,i) = FaceArea_DI(z_,i)*Misc1_I(i)
     v_B1x_rf(i) = v_B1x_rf(i) + DiffrB1n_DI(x_,i)
     v_B1x_lf(i) = v_B1x_lf(i) - DiffrB1n_DI(x_,i)
     v_B1y_rf(i) = v_B1y_rf(i) + DiffrB1n_DI(y_,i)
     v_B1y_lf(i) = v_B1y_lf(i) - DiffrB1n_DI(y_,i)
     v_B1z_rf(i) = v_B1z_rf(i) + DiffrB1n_DI(z_,i)
     v_B1z_lf(i) = v_B1z_lf(i) - DiffrB1n_DI(z_,i)
     B1n_rf(i) = cHalf*( B1n_rf(i)+B1n_lf(i))
     B1n_lf(i) = B1n_rf(i)
  end do


  do i = 1,nStrip
     mP_lf(i) = cHalf*(v_B1x_lf(i)**2+v_B1y_lf(i)**2+&
          v_B1z_lf(i)**2)                                     ! B1 squared/2 
     a2(i) = v_Ux_lf(i)**2+v_Uy_lf(i)**2+v_Uz_lf(i)**2        ! U^2
     C_E_lf(i) = inv_gm1*v_p_lf(i)+mP_lf(i)+&
          cHalf*v_rho_lf(i)*a2(i)                             !energy density 

     C_Bx_lf(i) = v_B0x(i)+v_B1x_lf(i)		          !total field	
     C_By_lf(i) = v_B0y(i)+v_B1y_lf(i)
     C_Bz_lf(i) = v_B0z(i)+v_B1z_lf(i)

     Bn_lf(i) = B1n_lf(i)+B0n(i)                              ! total normal MF
     mP_lf(i) = mP_lf(i)+v_B1x_lf(i)*v_B0x(i)+&
          v_B1y_lf(i)*v_B0y(i)+v_B1z_lf(i)*v_B0z(i)       !mP=magnetic pressure
     B0B0(i)=(v_B0x(i)**2+v_B0y(i)**2+v_B0z(i)**2)       !  B0^2
 
  end do


  do i=1,nStrip
     mP_rf(i) = cHalf*(v_B1x_rf(i)**2+v_B1y_rf(i)**2+&
          v_B1z_rf(i)**2) 
     a2(i) = v_Ux_rf(i)**2+v_Uy_rf(i)**2+v_Uz_rf(i)**2 
     C_E_rf(i) = inv_gm1*v_p_rf(i)+mP_rf(i)+&
          cHalf*v_rho_rf(i)*a2(i)

     C_Bx_rf(i) = v_B0x(i)+v_B1x_rf(i)
     C_By_rf(i) = v_B0y(i)+v_B1y_rf(i)
     C_Bz_rf(i) = v_B0z(i)+v_B1z_rf(i)

     Bn_rf(i)=B1n_rf(i)+B0n(i)
     mP_rf(i)=mP_rf(i)+v_B1x_rf(i)*v_B0x(i)+&
          v_B1y_rf(i)*v_B0y(i)+v_B1z_rf(i)*v_B0z(i)      ! magnetic pressure

     BFull2_I(i)=(v_B0x(i)+(v_B1x_rf(i)+v_B1x_lf(i))*cHalf)**2+&
          (v_B0y(i)+ (v_B1y_rf(i)+v_B1y_lf(i))*cHalf)**2+&
          (v_B0z(i)+(v_B1z_rf(i)+v_B1z_lf(i))*cHalf)**2     
  end do

  !all the propagation velocities and fluxes are mutiplied by the faceArea
  !begin the propagation velocities computation

  do i=1,nStrip
     GHalfFA2_I(i)=g_half*FaceArea2_I(i)

     a2(i)   = GHalfFA2_I(i)*(v_p_lf(i)+v_p_rf(i))
     Misc1_I(i)  = BFull2_I(i)*FaceArea2_I(i)+a2(i)
     discr(i)=Misc1_I(i)**2-a2(i)*(Bn_lf(i)+Bn_rf(i))**2     ! va2**2-4*Bn_hf**2*a2
     discr(i)=Misc1_I(i)+sqrt(max(discr(i),cZero))
     discr(i)=discr(i)/(v_rho_lf(i)+v_rho_rf(i))
     discr(i)=sqrt(discr(i))                             !fast sound speed*fA

     v_Un_hf(i) = cHalf*(v_Un_rf(i)+v_Un_lf(i))
!     v_cright_hf(i) =discr(i)+cHalf*(v_Un_lf(i)+v_Un_rf(i))        
!     v_cleft_hf(i)  =discr(i)-cHalf*(v_Un_lf(i)+v_Un_rf(i))          !With the oppsite sign!

     v_cright_hf(i) =discr(i) + v_Un_hf(i)       
     v_cleft_hf(i)  =discr(i) - v_Un_hf(i)

     v_max_hf(i)    =max(v_cright_hf(i), v_cleft_hf(i))              ! for further calculation of CFL


  end do

  do i=1,nStrip

     a2(i)   = GHalfFA2_I(i)*(v_p_lf(i)+v_p_lf(i))
     Misc1_I(i)  = (B0B0(i)+mP_lf(i)+mP_lf(i))&
          *FaceArea2_I(i)+a2(i)
     discr(i)=Misc1_I(i)**2-a2(i)*(Bn_lf(i)+Bn_lf(i))**2 ! va2**2-4*Bn_hf**2*a2
     discr(i)=Misc1_I(i)+sqrt(max(discr(i),cZero))
     discr(i)=discr(i)/(v_rho_lf(i)+v_rho_lf(i))
     discr(i)=sqrt(discr(i))                             !fast sound speed*fA

     v_cleft_hf(i) =max( discr(i)- v_Un_lf(i),v_cleft_hf(i),cZero) ! with oppisite SIGN

  end do

  do i=1,nStrip
     a2(i)   = GHalfFA2_I(i)*(v_p_rf(i)+v_p_rf(i))
     Misc1_I(i)  = (B0B0(i)+mP_rf(i)+mP_rf(i))&
          *FaceArea2_I(i)+a2(i)
     discr(i)=Misc1_I(i)**2-a2(i)*(Bn_rf(i)+Bn_rf(i))**2     ! va2**2-4*Bn_hf**2*a2
     discr(i)=Misc1_I(i)+sqrt(max(discr(i),cZero))
     discr(i)=discr(i)/(v_rho_rf(i)+v_rho_rf(i))
     discr(i)=sqrt(discr(i))                             !fast sound speed*fA

     v_cright_hf(i) =max(v_Un_rf(i)+discr(i),v_cright_hf(i),cZero)         

     a2(i)=v_cleft_hf(i)/(v_cleft_hf(i)+v_cright_hf(i))  !weight coefficient
  end do

!------------------------------------------------------------------
!calculate U_l and U_r and  delta_U
!------------------------------------------------------------------
  do i=1,nStrip
     U_l(rho_,i)=v_rho_lf(i)     
     U_l(rhoUx_,i)=v_Ux_lf(i)*v_rho_lf(i)
     U_l(rhoUy_,i)=v_Uy_lf(i)*v_rho_lf(i)
     U_l(rhoUz_,i)=v_Uz_lf(i)*v_rho_lf(i)
     U_l(Bx_,i)=	C_Bx_lf(i)
     U_l(By_,i)=	C_By_lf(i)
     U_l(Bz_,i)=	C_Bz_lf(i)
     U_l(Energy_,i)=  C_E_lf(i)
     U_l(P_,i)=v_p_lf(i)
     U_l(Energy_+1:nFlux,i) =0.0
     U_l(rho_+1:rho_+MaxSpecies,i)=v_rhoSp_lf(1:MaxSpecies,i)
  end do
  do i=1,nStrip
     U_r(rho_,i)=v_rho_rf(i)     
     U_r(rhoUx_,i)=v_Ux_rf(i)*v_rho_rf(i)
     U_r(rhoUy_,i)=v_Uy_rf(i)*v_rho_rf(i)
     U_r(rhoUz_,i)=v_Uz_rf(i)*v_rho_rf(i)
     U_r(Bx_,i)=	C_Bx_rf(i)
     U_r(By_,i)=	C_By_rf(i)
     U_r(Bz_,i)=	C_Bz_rf(i)
     U_r(Energy_,i)=  C_E_rf(i)
     U_r(P_,i)=v_p_rf(i)
     U_r(Energy_+1:nFlux,i) =0.0
     U_r(rho_+1:rho_+MaxSpecies,i)=v_rhoSp_rf(1:MaxSpecies,i)
   end do
   do i=1,nStrip
!      v_diff_hf(i) = cHalf*(v_cright_hf(i)+v_cleft_hf(i))     
      v_diff_hf(i) = cHalf*(v_cright_hf(i)-v_cleft_hf(i))     
!      if (abs( v_diff_hf(i) - v_Un_hf(i) ) > 1.0e-3 ) &
!       write(*,*)'v_diff, v_Un=', v_diff_hf(i),v_Un_hf(i)
         
!      v_Un_hf(i) = cHalf*(v_cright_hf(i)-v_cleft_hf(i))
      min_v(i)=min(v_cleft_hf(i), v_cright_hf(i))
      delta_U(1:nFlux,i)=U_r(1:nFlux,i)-U_l(1:nFlux,i)
      dFlux(1:nFlux,i)=v_cright_hf(i)*U_r(1:nFlux,i)+ &
                       v_cleft_hf(i)*U_l(1:nFlux,i)-&
                       v_Un_hf(i)*delta_U(1:nFlux,i)
      lambda_r(i)=v_cright_hf(i)
      lambda_l(i)=-v_cleft_hf(i)
 
   end do
!--------------------------------------------------------



  !\
  ! Conserved state solutions, Fluxes
  !/
  do i=1,nStrip
     mP_lf(i)=mP_lf(i)+v_P_lf(i)                             ! full pressure
     Misc1_I(i)=v_Ux_lf(i)*v_B1x_lf(i)+v_Uy_lf(i)*v_B1y_lf(i)+&  !UdotB1
          v_Uz_lf(i)*v_B1z_lf(i)
     v_cleft_hf(i)=v_Un_lf(i)+v_cleft_hf(i)    

     ! advection velocity is shifted
     ! to calculate both the flux and
     ! diffusion: F_L-D_L*U_L is claculated
     ! below

     Fl_lf(rho_,i)=v_cleft_hf(i)*v_rho_lf(i)
     Fl_lf(rho_+1:rho_+MaxSpecies,i)=v_cleft_hf(i)*v_rhoSp_lf(1:MaxSpecies,i)
     Fl_lf(rhoUx_,i)=v_Ux_lf(i)*v_rho_lf(i)*v_cleft_hf(i)-&
          (Bn_lf(i)*v_B1x_lf(i)+B1n_lf(i)*v_B0x(i))
     Fl_lf(rhoUy_,i)=v_Uy_lf(i)*v_rho_lf(i)*v_cleft_hf(i)-&
          (Bn_lf(i)*v_B1y_lf(i)+B1n_lf(i)*v_B0y(i))
     Fl_lf(rhoUz_,i)=v_Uz_lf(i)*v_rho_lf(i)*v_cleft_hf(i)-&
          (Bn_lf(i)*v_B1z_lf(i)+B1n_lf(i)*v_B0z(i))
     Fl_lf(Bx_,i)=	C_Bx_lf(i)*v_cleft_hf(i)-Bn_lf(i)*v_Ux_lf(i)
     Fl_lf(By_,i)=	C_By_lf(i)*v_cleft_hf(i)-Bn_lf(i)*v_Uy_lf(i)
     Fl_lf(Bz_,i)=	C_Bz_lf(i)*v_cleft_hf(i)-Bn_lf(i)*v_Uz_lf(i)
     Fl_lf(P_,i)=v_cleft_hf(i)*v_p_lf(i)
     Fl_lf(Energy_,i)=  C_E_lf(i)*v_cleft_hf(i)+v_Un_lf(i)*mP_lf(i)-&
          Bn_lf(i)*Misc1_I(i)                   ! the last term is B(UdotB1)
     Fl_lf(Energy_+1,i)=mP_lf(i)
     Fl_lf(Energy_+2,i)=v_Un_lf(i)                  
     ! the term mP*{\delta}_ij gives
     !the input mP*\overrightarrow(fA) to the flux
  end do
  do i=1,nStrip
     mP_rf(i)=mP_rf(i)+v_P_rf(i)
     Misc1_I(i)=v_Ux_rf(i)*v_B1x_rf(i)+v_Uy_rf(i)*&
          v_B1y_rf(i)+v_Uz_rf(i)*v_B1z_rf(i)
     v_cright_hf(i)=v_Un_rf(i)-v_cright_hf(i)
     Fl_rf(rho_,i)=v_cright_hf(i)*v_rho_rf(i)
     Fl_rf(rho_+1:rho_+MaxSpecies,i)=v_cright_hf(i)*v_rhoSp_rf(1:MaxSpecies,i)
     Fl_rf(rhoUx_,i)=v_Ux_rf(i)*v_rho_rf(i)*v_cright_hf(i)-&
          (Bn_rf(i)*v_B1x_rf(i)+B1n_rf(i)*v_B0x(i))
     Fl_rf(rhoUy_,i)=v_Uy_rf(i)*v_rho_rf(i)*v_cright_hf(i)-&
          (Bn_rf(i)*v_B1y_rf(i)+B1n_rf(i)*v_B0y(i))
     Fl_rf(rhoUz_,i)=v_Uz_rf(i)*v_rho_rf(i)*v_cright_hf(i)-&
          (Bn_rf(i)*v_B1z_rf(i)+B1n_rf(i)*v_B0z(i))
     Fl_rf(Bx_,i)=C_Bx_rf(i)*v_cright_hf(i)-Bn_rf(i)*v_Ux_rf(i)
     Fl_rf(By_,i)=C_By_rf(i)*v_cright_hf(i)-Bn_rf(i)*v_Uy_rf(i)
     Fl_rf(Bz_,i)=C_Bz_rf(i)*v_cright_hf(i)-Bn_rf(i)*v_Uz_rf(i)
     Fl_rf(P_,i)=v_cright_hf(i)*v_p_rf(i)
     Fl_rf(Energy_,i) =C_E_rf(i)*v_cright_hf(i)+v_Un_rf(i)*mP_rf(i)-&
          Bn_rf(i)*Misc1_I(i)
     Fl_rf(Energy_+1,i)=mP_rf(i)
     Fl_rf(Energy_+2,i)=v_Un_rf(i)
  end do

!----------for HLLEL flux------------------------------------------
    do i=1,nStrip
     Fl_rf(rhoUx_,i) =  Fl_rf(rhoUx_,i) + Fl_rf(Energy_+1,i)*FaceArea_DI(x_,i)
     Fl_rf(rhoUy_,i) =  Fl_rf(rhoUy_,i) + Fl_rf(Energy_+1,i)*FaceArea_DI(y_,i)
     Fl_rf(rhoUz_,i) =  Fl_rf(rhoUz_,i) + Fl_rf(Energy_+1,i)*FaceArea_DI(z_,i)
     Fl_lf(rhoUx_,i) =  Fl_lf(rhoUx_,i) + Fl_lf(Energy_+1,i)*FaceArea_DI(x_,i)
     Fl_lf(rhoUy_,i) =  Fl_lf(rhoUy_,i) + Fl_lf(Energy_+1,i)*FaceArea_DI(y_,i)
     Fl_lf(rhoUz_,i) =  Fl_lf(rhoUz_,i) + Fl_lf(Energy_+1,i)*FaceArea_DI(z_,i)
!     dFlux(1:nFlux,i)=abs(Fl_rf(1:nFlux,i)-Fl_lf(1:nFlux,i)+&
!           v_diff_hf(i)*(U_r(1:nFlux,i)+U_l(1:nFlux,i)))
     dFlux(1:nFlux,i)=abs( Fl_rf(1:nFlux,i)-Fl_lf(1:nFlux,i)+&
           dFlux(1:nFlux,i) )

   end do
!-----------------------------------------------------------


  !\	
  !Full flux=a2F_r+(1-a2)F_L+d(U_l-U_r)  
  !
  do i=1,nStrip
     do n=1,nFlux
        Full_Flux(n,i) = Fl_lf(n,i)+(Fl_rf(n,i)-Fl_lf(n,i))*a2(i) 
     end do
  end do

  do i=1,nStrip
  end do
  do i=1,nStrip
     Fl_rf(Bx_:Bz_,i) =   v_max_hf(i)*DiffrB1n_DI(:,i)
     Full_Flux(Bx_:Bz_,i) = Full_Flux(Bx_:Bz_,i)+  Fl_rf(Bx_:Bz_,i)
     Full_Flux(Energy_,i) = Full_Flux(Energy_,i)+ &
          Fl_rf(Bx_,i)*v_B1x_lf(i)+Fl_rf(By_,i)*&
          v_B1y_lf(i)+Fl_rf(Bz_,i)*v_B1z_lf(i)
end do

!------------------------------------------------------------------
! HLLEL interface fluxes
!------------------------------------------------------------------
  
  do i=1,nStrip

     if( min_V(i)>0.0 ) then

        v_cstar_hf(i) = sqrt( FaceArea2_I(i)*g*(v_p_lf(i)+v_p_rf(i))/&
             (v_rho_lf(i)+v_rho_rf(i)) )
        
!        if(UseMultiSpecies) then
!           Alpha = 1.00 - ( sum(dFlux(rho_+1:E_,i))+sum(dFlux(12:nFlux,i)) )/&
!                ( (sum(abs(delta_U(rho_+1:E_,i))+1.0e-08) + &
!                sum(abs(delta_U(12:nFlux,i))+1.0e-08)) &
!                *v_cstar_hf(i) )

!        else
!           Alpha = 1.00 - (dFlux(rho_,i)+dFlux(Energy_,i)+sum(dFlux(rhoUx_:Bz_,i)))/&
!                ( (abs(delta_U(rho_,i))+abs(delta_U(Energy_,i))+sum(abs(delta_U(rhoUx_:Bz_,i)))+1.0e-08)&
!                *v_cstar_hf(i) )

           Alpha = 1.00 - (dFlux(Energy_,i)+sum(dFlux(rho_+1:Bz_,i)))/&
                ( (abs(delta_U(Energy_,i))+sum(abs(delta_U(rho_+1:Bz_,i)))+1.0e-08)&
                *v_cstar_hf(i) )

!        end if
        Alpha = max( 0.00, Alpha)

        if(v_Un_hf(i)>=0.0) then
           do n=1,nFlux
              Full_Flux(n,i) = Full_Flux(n,i)+ Alpha *&
                   lambda_l(i)*(v_Un_hf(i)-lambda_r(i))*delta_U(n,i)&
                   /(lambda_r(i)-lambda_l(i))
           end do
        else
           do n=1,nFlux
              Full_Flux(n,i) = Full_Flux(n,i)+ Alpha *&
                   lambda_r(i)*(v_Un_hf(i)-lambda_l(i))*delta_U(n,i)&
                   /(lambda_r(i)-lambda_l(i))

           end do
        end if
     end if !v_cleft_hf(i) > 0.0 .and. v_cright_hf(i) > 0.0
  end do
end subroutine get_flux_mhdLinde
