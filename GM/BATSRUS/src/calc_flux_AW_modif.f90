subroutine option_awflux(TrueOption,NameOption)

  logical, intent(out) :: TrueOption
  character (len=40), intent(out) :: NameOption

  TrueOption  = .true.
  NameOption  = 'AW FLUX Sokolov 1.1'

end subroutine option_awflux

Module ModFlux
  use ModVarIndexes,ONLY: Energy_
  implicit none
  integer, parameter :: nFlux=Energy_+2   ! +1 for pressure +2 for u_n averaged 
  integer,parameter  :: MaxStrip=10

  ! symmetric face variables
  real, dimension(3):: fA  ! the faceArea vector


  real, dimension(MaxStrip) :: v_B0x,v_B0y,v_B0z,B0n
  real, dimension(MaxStrip) :: v_max_hf          ! Maximum perturbation velocity
  real, dimension(MaxStrip) :: x_face,y_face,z_face

  ! left face 
  real, dimension(MaxStrip) :: v_rho_lf,v_p_lf,v_er_lf
  real, dimension(MaxStrip) :: v_Ux_lf,v_B1x_lf,v_Uy_lf,v_B1y_lf,v_Uz_lf,v_B1z_lf
  real, dimension(MaxStrip) :: v_Un_lf,B1n_lf

  ! right face	
  real, dimension(MaxStrip) :: v_rho_rf,v_p_rf,v_er_rf
  real, dimension(MaxStrip) :: v_Ux_rf,v_B1x_rf,v_Uy_rf,v_B1y_rf,v_Uz_rf,v_B1z_rf
  real, dimension(MaxStrip) :: v_Un_rf,B1n_rf

  ! flux
  real, dimension(nFlux,MaxStrip) :: Full_Flux     ! Flux

  ! These common blocks improve speed slightly
  common/igors/ v_rho_lf,v_Ux_lf,v_Uy_lf,v_Uz_lf,&
       v_B1x_lf,v_B1y_lf,v_B1z_lf,v_p_lf,v_er_lf,&
       v_B0x,v_B0y,v_B0z,&
       v_rho_rf,v_Ux_rf,v_Uy_rf,v_Uz_rf,&
       v_B1x_rf,v_B1y_rf,v_B1z_rf,v_p_rf,&
       v_er_rf,v_Un_lf,v_Un_rf,B1n_lf,B1n_rf,B0n,&
       Full_Flux,v_max_hf,fA
end module ModFlux

!=============================================================================

subroutine calc_flux_AW(DoResChangeOnly)
  use ModProcMH
  use ModMain
  use ModGeometry, ONLY : fAx_BLK,fAy_BLK,fAz_BLK,&
       x_BLK,y_BLK,z_BLK
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
  real, dimension(MaxStrip)::v_diffrB1n



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

  !\
  fA(1) = fAx_BLK(globalBLK);  fA(2)=cZero ;  fA(3)=cZero
  !/

  do iStart=iBegin,iEnd,MaxStrip 
     nStrip=min(MaxStrip,1+iEnd-iStart) 
     do k=kBegin,kEnd
        oktest_row=oktest_me.and.k==Ktest ! for test output
        do j=jBegin,jEnd
           oktest_row=oktest_row.and.j==Jtest !for test output
           do i=1,nStrip
              iStrip=i+iStart-1		  

              x_face(i) = (x_BLK(iStrip,j,k,globalBLK)+x_BLK(iStrip-1,j,k,globalBLK))*cHalf
              y_face(i) = y_BLK(iStrip,j,k,globalBLK)
              z_face(i) = z_BLK(iStrip,j,k,globalBLK)

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
              v_Ux_lf(i)  = LeftState_VX(Ux_ ,iStrip,j,k)
              v_Uy_lf(i)  = LeftState_VX(Uy_ ,iStrip,j,k)
              v_Uz_lf(i)  = LeftState_VX(Uz_ ,iStrip,j,k)
              v_B1x_lf(i) = LeftState_VX(Bx_ ,iStrip,j,k)
              v_B1y_lf(i) = LeftState_VX(By_ ,iStrip,j,k)
              v_B1z_lf(i) = LeftState_VX(Bz_ ,iStrip,j,k)
              v_er_lf(i)  = LeftState_VX(EnergyRL_,iStrip,j,k)
              v_p_lf(i)   = LeftState_VX(P_   ,iStrip,j,k)
              !\
              ! Right face
              !/
              v_rho_rf(i) = RightState_VX(rho_,iStrip,j,k)
              v_Ux_rf(i)  = RightState_VX(Ux_ ,iStrip,j,k)
              v_Uy_rf(i)  = RightState_VX(Uy_ ,iStrip,j,k)
              v_Uz_rf(i)  = RightState_VX(Uz_ ,iStrip,j,k)
              v_B1x_rf(i) = RightState_VX(Bx_ ,iStrip,j,k) 
              v_B1y_rf(i) = RightState_VX(By_ ,iStrip,j,k) 
              v_B1z_rf(i) = RightState_VX(Bz_ ,iStrip,j,k) 
              v_er_rf(i)  = RightState_VX(EnergyRL_,iStrip,j,k)
              v_p_rf(i)   = RightState_VX(P_   ,iStrip,j,k) 
             
              !\
              ! Normal components
              !/
              v_diffrB1n(i)=cHalf*(v_B1x_lf(i)-v_B1x_rf(i)) 
	      v_B1x_lf(i)=v_B1x_rf(i)+v_diffrB1n(i)
	      v_B1x_rf(i)=v_B1x_lf(i)
              v_Un_lf(i)=fA(1)*v_Ux_lf(i)
              v_Un_rf(i)=fA(1)*v_Ux_rf(i)
              B0n(i)=fA(1)*v_B0x(i)      
              B1n_lf(i)=fA(1)*v_B1x_lf(i)
              B1n_rf(i)=fA(1)*v_B1x_rf(i)

           end do
           !\
           ! get fluxes
           !/

           call get_flux_mhdLAW(nStrip)	
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
              Flux_VX(rho_  ,iStrip,j,k) = Full_Flux(1,i)
              Flux_VX(rhoUx_,iStrip,j,k) = Full_Flux(2,i)+fA(1)*Full_Flux(Energy_+1,i)    
              Flux_VX(rhoUy_,iStrip,j,k) = Full_Flux(3,i) 
              Flux_VX(rhoUz_,iStrip,j,k) = Full_Flux(4,i) 
              v_diffrB1n(i)=v_max_hf(i)*v_diffrB1n(i)                            
              Flux_VX(Bx_   ,iStrip,j,k) = Full_Flux(5,i)+v_diffrB1n(i)               
              Flux_VX(By_   ,iStrip,j,k) = Full_Flux(6,i)
              Flux_VX(Bz_   ,iStrip,j,k) = Full_Flux(7,i)
              Flux_VX(EnergyRL_,iStrip,j,k) = Full_Flux(EnergyRL_,i)   
              Flux_VX(P_     ,iStrip,j,k) = Full_Flux(P_,i)
              Flux_VX(Energy_,iStrip,j,k) = Full_Flux(Energy_,i)+v_diffrB1n(i)*v_B1x_lf(i) 
              UDotFA_x(iStrip,j,k)     = Full_Flux(Energy_+2,i)
 
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
  !\ 
  fA(1)=cZero; fA(2)= fAy_BLK(globalBLK);  fA(3)=cZero
  !/
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

              x_face(i) = x_BLK(iStrip,j,k,globalBLK)
              y_face(i) = (y_BLK(iStrip,j,k,globalBLK)+y_BLK(iStrip,j-1,k,globalBLK))*cHalf
              z_face(i) = z_BLK(iStrip,j,k,globalBLK)

              v_B0x(i) = B0xFace_y_BLK(iStrip,j,k,globalBLK)
              v_B0y(i) = B0yFace_y_BLK(iStrip,j,k,globalBLK)
              v_B0z(i) = B0zFace_y_BLK(iStrip,j,k,globalBLK)

              !\
              ! Left face
              !/
              v_rho_lf(i) = LeftState_VY(rho_,iStrip,j,k)
              v_Ux_lf(i)  = LeftState_VY(Ux_ ,iStrip,j,k)
              v_Uy_lf(i)  = LeftState_VY(Uy_ ,iStrip,j,k)
              v_Uz_lf(i)  = LeftState_VY(Uz_ ,iStrip,j,k)
              v_B1x_lf(i) = LeftState_VY(Bx_ ,iStrip,j,k)
              v_B1y_lf(i) = LeftState_VY(By_ ,iStrip,j,k)
              v_B1z_lf(i) = LeftState_VY(Bz_ ,iStrip,j,k)
              v_er_lf(i)  = LeftState_VY(EnergyRL_,iStrip,j,k)
              v_p_lf(i)   = LeftState_VY(P_   ,iStrip,j,k)
              !\
              ! Right face
              !/
              v_rho_rf(i) = RightState_VY(rho_,iStrip,j,k)
              v_Ux_rf(i)  = RightState_VY(Ux_ ,iStrip,j,k)
              v_Uy_rf(i)  = RightState_VY(Uy_ ,iStrip,j,k)
              v_Uz_rf(i)  = RightState_VY(Uz_ ,iStrip,j,k)
              v_B1x_rf(i) = RightState_VY(Bx_ ,iStrip,j,k)
              v_B1y_rf(i) = RightState_VY(By_ ,iStrip,j,k) 
              v_B1z_rf(i) = RightState_VY(Bz_ ,iStrip,j,k) 
              v_er_rf(i)  = RightState_VY(EnergyRL_,iStrip,j,k)
              v_p_rf(i)   = RightState_VY(P_   ,iStrip,j,k) 

              !\
              ! normal components
              !/
	      v_diffrB1n(i) = cHalf*(v_B1y_lf(i)-v_B1y_rf(i))
	      v_B1y_lf(i) = v_B1y_rf(i)+v_diffrB1n(i)
	      v_B1y_rf(i) = v_B1y_lf(i)
              v_Un_lf(i) = fA(2)*v_Uy_lf(i)
              v_Un_rf(i) = fA(2)*v_Uy_rf(i)
              B0n(i) = fA(2)*v_B0y(i)      
              B1n_lf(i) = fA(2)*v_B1y_lf(i)
              B1n_rf(i) = fA(2)*v_B1y_rf(i)
           end do
           !\
           ! get fluxes
           !/

           call get_flux_mhdLAW(nStrip)
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

              Flux_VY(rho_  ,iStrip,j,k) = Full_Flux(1,i)
              Flux_VY(rhoUx_,iStrip,j,k) = Full_Flux(2,i)  
              Flux_VY(rhoUy_,iStrip,j,k) = Full_Flux(3,i)+fA(2)*Full_Flux(Energy_+1,i)
              Flux_VY(rhoUz_,iStrip,j,k) = Full_Flux(4,i) 
              Flux_VY(Bx_   ,iStrip,j,k) = Full_Flux(5,i)
              v_diffrB1n(i) = v_max_hf(i)*v_diffrB1n(i)                        
              Flux_VY(By_   ,iStrip,j,k) = Full_Flux(6,i)+v_diffrB1n(i)           
              Flux_VY(Bz_   ,iStrip,j,k) = Full_Flux(7,i)
              Flux_VY(EnergyRL_,iStrip,j,k) = Full_Flux(EnergyRL_,i)
              Flux_VY(P_     ,iStrip,j,k) = Full_Flux(P_,i)
              Flux_VY(Energy_,iStrip,j,k) = Full_Flux(Energy_,i)+v_diffrB1n(i)*v_B1y_lf(i)
              UDotFA_Y(iStrip,j,k)     = Full_Flux(Energy_+2,i)

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
  
  !\ The last portion, relating to the cartesian grid only
  fA(1)=0.0;  fA(2)=0.0;  fA(3) = fAz_BLK(globalBLK)
  !/
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

              x_face(i) = x_BLK(iStrip,j,k,globalBLK)
              y_face(i) = y_BLK(iStrip,j,k,globalBLK)
              z_face(i) = (z_BLK(iStrip,j,k,globalBLK)+z_BLK(iStrip,j,k-1,globalBLK))*cHalf

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
              v_Ux_lf(i)  = LeftState_VZ(Ux_ ,iStrip,j,k)
              v_Uy_lf(i)  = LeftState_VZ(Uy_ ,iStrip,j,k)
              v_Uz_lf(i)  = LeftState_VZ(Uz_ ,iStrip,j,k)
              v_B1x_lf(i) = LeftState_VZ(Bx_ ,iStrip,j,k)
              v_B1y_lf(i) = LeftState_VZ(By_ ,iStrip,j,k)
              v_B1z_lf(i) = LeftState_VZ(Bz_ ,iStrip,j,k)
              v_er_lf(i)  = LeftState_VZ(EnergyRL_,iStrip,j,k)
              v_p_lf(i)   = LeftState_VZ(P_   ,iStrip,j,k)

              !\
              ! Right face
              !/
              v_rho_rf(i) = RightState_VZ(rho_,iStrip,j,k)
              v_Ux_rf(i)  = RightState_VZ(Ux_ ,iStrip,j,k) 
              v_Uy_rf(i)  = RightState_VZ(Uy_ ,iStrip,j,k) 
              v_Uz_rf(i)  = RightState_VZ(Uz_ ,iStrip,j,k) 
              v_B1x_rf(i) = RightState_VZ(Bx_ ,iStrip,j,k) 
              v_B1y_rf(i) = RightState_VZ(By_ ,iStrip,j,k) 
              v_B1z_rf(i) = RightState_VZ(Bz_ ,iStrip,j,k)
              v_er_rf(i)  = RightState_VZ(EnergyRL_,iStrip,j,k)
              v_p_rf(i)   = RightState_VZ(P_   ,iStrip,j,k)
              !\
              ! normal components
              !/
              v_diffrB1n(i)=cHalf*(v_B1z_lf(i)-v_B1z_rf(i))
              v_B1z_lf(i)=v_diffrB1n(i)+v_B1z_rf(i)
              v_B1z_rf(i)=v_B1z_lf(i)
              v_Un_lf(i)=fA(3)*v_Uz_lf(i) !normal velocities*fA
              v_Un_rf(i)=fA(3)*v_Uz_rf(i)
              B0n(i)=fA(3)*v_B0z(i)       !normal B-zero field
              B1n_lf(i)=fA(3)*v_B1z_lf(i)
              B1n_rf(i)=fA(3)*v_B1z_rf(i)
           end do
           !\
           ! get fluxes
           !/
           call get_flux_mhdLAW(nStrip)	

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

              Flux_VZ(rho_  ,iStrip,j,k) = Full_Flux(1,i)
              Flux_VZ(rhoUx_,iStrip,j,k) = Full_Flux(2,i)
              Flux_VZ(rhoUy_,iStrip,j,k) = Full_Flux(3,i)
              Flux_VZ(rhoUz_,iStrip,j,k) = Full_Flux(4,i)+fA(3)*Full_Flux(Energy_+1,i)    
              Flux_VZ(Bx_   ,iStrip,j,k) = Full_Flux(5,i)
              Flux_VZ(By_   ,iStrip,j,k) = Full_Flux(6,i)
              v_diffrB1n(i)=v_max_hf(i)*v_diffrB1n(i)
              Flux_VZ(Bz_   ,iStrip,j,k) = Full_Flux(7,i)+v_diffrB1n(i) 
              Flux_VZ(EnergyRL_,iStrip,j,k) = Full_Flux(EnergyRL_,i)  
              Flux_VZ(P_     ,iStrip,j,k) = Full_Flux(P_,i)            
              Flux_VZ(Energy_,iStrip,j,k) = Full_Flux(Energy_,i)+v_diffrB1n(i)*v_B1z_lf(i)
              UDotFA_z(iStrip,j,k)     = Full_Flux(Energy_+2,i)

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
    real    :: v_rho_hf,v_p_hf
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

end subroutine calc_flux_AW

!==========================================================================
subroutine get_flux_mhdLAW(nStrip)
  use ModFlux
  use ModNumConst
  use ModVarIndexes, ONLY : rho_, rhoUx_, rhoUy_, rhoUz_, Bx_, By_, &
                      Bz_, Energy_, P_, EnergyRL_
  use ModPhysics, ONLY : g_half, inv_gm1
  implicit none
  integer, intent(in) :: nStrip                   ! the length of "strip", of flux
  real, dimension(MaxStrip) ::Bn_lf,Bn_rf
  real, dimension(MaxStrip) :: a2,vA2,discr, mP_lf, mP_rf, B0B0
  real, dimension(MaxStrip) :: v_cright_hf, v_cleft_hf,v_Ua_lf,v_Ua_rf
  real, dimension(nFlux,MaxStrip):: Fl_lf, Fl_rf
  real, dimension(MaxStrip)::C_E_lf ,C_Bx_lf,C_By_lf,C_Bz_lf,C_P_lf
  real, dimension(MaxStrip)::C_E_rf ,C_Bx_rf,C_By_rf,C_Bz_rf,C_P_rf
  integer :: i,n
  real :: GHalfFA2, fA2, sign_br

  !__________________________________________________________________________

  fA2 = sum(fA(1:3)*fA(1:3))




  do i = 1,nStrip
     mP_lf(i) = cHalf*(v_B1x_lf(i)**2+v_B1y_lf(i)**2+&
          v_B1z_lf(i)**2)                                 ! B1 squared/2 
     a2(i) = v_Ux_lf(i)**2+v_Uy_lf(i)**2+v_Uz_lf(i)**2    ! U^2
     C_E_lf(i) = inv_gm1*v_p_lf(i)+mP_lf(i)+&
          cHalf*v_rho_lf(i)*a2(i)                         ! energy density 

     C_Bx_lf(i) = v_B0x(i)+v_B1x_lf(i)	     	          ! total field	
     C_By_lf(i) = v_B0y(i)+v_B1y_lf(i)
     C_Bz_lf(i) = v_B0z(i)+v_B1z_lf(i)

     Bn_lf(i) = B1n_lf(i)+B0n(i)                          ! total normal MF
     mP_lf(i) = mP_lf(i)+v_B1x_lf(i)*v_B0x(i)+&
          v_B1y_lf(i)*v_B0y(i)+v_B1z_lf(i)*v_B0z(i)       ! mP=magnetic pressure
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
     B0B0(i)=(v_B0x(i)**2+v_B0y(i)**2+v_B0z(i)**2)       ! B0^2
  end do

  !all the propagation velocities and fluxes are mutiplied by the faceArea
  !begin the propagation velocities computation

  GHalfFA2=g_half*fA2

  do i=1,nStrip



                                                  
     a2(i)   = GHalfFA2*(v_p_lf(i)+v_p_rf(i)) 
     vA2(i)  = (B0B0(i)+mP_lf(i)+mP_rf(i))&
          *fA2+a2(i)	                           

     discr(i)=vA2(i)**2-a2(i)*(Bn_lf(i)+Bn_rf(i))**2     ! va2**2-4*Bn_hf**2*a2
     discr(i)=vA2(i)+sqrt(max(discr(i),cZero))
     discr(i)=discr(i)/(v_rho_lf(i)+v_rho_rf(i))
     discr(i)=sqrt(discr(i))                             ! fast sound speed*fA
     v_cright_hf(i) =max( v_Un_lf(i) + discr(i),&        ! AW right velocity*fA 
          v_Un_rf(i)+discr(i),cZero)         
     v_cleft_hf(i) =max( discr(i)- v_Un_lf(i),&    
          discr(i)-v_Un_rf(i),cZero)                     ! with oppisite SIGN
     v_max_hf(i)=max(v_cright_hf(i),v_cleft_hf(i))       ! for further calculation of CFL
     a2(i)=v_cleft_hf(i)/(v_cleft_hf(i)+v_cright_hf(i))  ! weight coefficient
  end do

  !\
  ! Conserved state solutions, Fluxes
  !/
  do i=1,nStrip
     mP_lf(i)=mP_lf(i)+v_P_lf(i)                             ! full pressure
     vA2(i)=v_Ux_lf(i)*v_B1x_lf(i)+v_Uy_lf(i)*v_B1y_lf(i)+&  ! UdotB1
          v_Uz_lf(i)*v_B1z_lf(i)
     v_cleft_hf(i)=v_Un_lf(i)+v_cleft_hf(i)
     ! advection velocity is shifted
     ! to calculate both the flux and
     ! diffusion: F_L-D_L*U_L is claculated
     ! below

     Fl_lf(rho_,i)=v_cleft_hf(i)*v_rho_lf(i)
     Fl_lf(rhoUx_,i)=v_Ux_lf(i)*v_rho_lf(i)*v_cleft_hf(i)-&
          (Bn_lf(i)*v_B1x_lf(i)+B1n_lf(i)*v_B0x(i))      
     Fl_lf(rhoUy_,i)=v_Uy_lf(i)*v_rho_lf(i)*v_cleft_hf(i)-&
          (Bn_lf(i)*v_B1y_lf(i)+B1n_lf(i)*v_B0y(i))      
     Fl_lf(rhoUz_,i)=v_Uz_lf(i)*v_rho_lf(i)*v_cleft_hf(i)-&
          (Bn_lf(i)*v_B1z_lf(i)+B1n_lf(i)*v_B0z(i))     
     Fl_lf(Bx_,i)   =C_Bx_lf(i)*v_cleft_hf(i)-Bn_lf(i)*v_Ux_lf(i)
     Fl_lf(By_,i)   =C_By_lf(i)*v_cleft_hf(i)-Bn_lf(i)*v_Uy_lf(i)
     Fl_lf(Bz_,i)   =C_Bz_lf(i)*v_cleft_hf(i)-Bn_lf(i)*v_Uz_lf(i)
     Fl_lf(EnergyRL_  ,i)=v_cleft_hf(i)*v_er_lf(i)
     Fl_lf(P_     ,i)=v_cleft_hf(i)*v_p_lf(i)
     Fl_lf(Energy_,i)=C_E_lf(i)*v_cleft_hf(i)+v_Un_lf(i)*mP_lf(i)-&
          Bn_lf(i)*vA2(i)                   ! the last term is B(UdotB1)
     Fl_lf(Energy_+1,i)=mP_lf(i)
     Fl_lf(Energy_+2,i)=v_Un_lf(i)                  
     ! the term mP*{\delta}_ij gives
     !the input mP*\overrightarrow(fA) to the flux
  end do

  !\
  ! Conserved state solutions, Fluxes
  !/
  do i=1,nStrip
     mP_rf(i)=mP_rf(i)+v_P_rf(i)                          ! full pressure
     vA2(i)=v_Ux_rf(i)*v_B1x_rf(i)+v_Uy_rf(i)*&           ! UdotB1
          v_B1y_rf(i)+v_Uz_rf(i)*v_B1z_rf(i)
     v_cright_hf(i)=v_Un_rf(i)-v_cright_hf(i)

     Fl_rf(rho_,i)=v_cright_hf(i)*v_rho_rf(i)
     Fl_rf(rhoUx_,i)=v_Ux_rf(i)*v_rho_rf(i)*v_cright_hf(i)-&
          (Bn_rf(i)*v_B1x_rf(i)+B1n_rf(i)*v_B0x(i)) 
     Fl_rf(rhoUy_,i)=v_Uy_rf(i)*v_rho_rf(i)*v_cright_hf(i)-&
          (Bn_rf(i)*v_B1y_rf(i)+B1n_rf(i)*v_B0y(i))
     Fl_rf(rhoUz_,i)=v_Uz_rf(i)*v_rho_rf(i)*v_cright_hf(i)-&
          (Bn_rf(i)*v_B1z_rf(i)+B1n_rf(i)*v_B0z(i))      
     Fl_rf(Bx_,i)   =C_Bx_rf(i)*v_cright_hf(i)-Bn_rf(i)*v_Ux_rf(i)
     Fl_rf(By_,i)   =C_By_rf(i)*v_cright_hf(i)-Bn_rf(i)*v_Uy_rf(i)
     Fl_rf(Bz_,i)   =C_Bz_rf(i)*v_cright_hf(i)-Bn_rf(i)*v_Uz_rf(i)
     Fl_rf(EnergyRL_  ,i)=v_cright_hf(i)*v_er_rf(i)
     Fl_rf(P_     ,i)=v_cright_hf(i)*v_p_rf(i)
     Fl_rf(Energy_,i)=C_E_rf(i)*v_cright_hf(i)+v_Un_rf(i)*mP_rf(i)-&
          Bn_rf(i)*vA2(i)
     Fl_rf(Energy_+1,i)=mP_rf(i)
     Fl_rf(Energy_+2,i)=v_Un_rf(i)
  end do
  !\	
  ! Full flux=a2F_r+(1-a2)F_L+d(U_l-U_r)  
  !/
  do i=1,nStrip
     do n=1,nFlux
        Full_Flux(n,i) = Fl_lf(n,i)+(Fl_rf(n,i)-Fl_lf(n,i))*a2(i) 
     end do
  end do



end subroutine get_flux_mhdLAW

