!^CFG COPYRIGHT UM
module ModLimiter
  use ModSize
  use ModVarIndexes,ONLY:nVar
  implicit none

  integer,parameter                  :: MaxIJK= max(nI,nJ,nK)
      ! the number of variables

  ! index ranges for optimized slope limiter calculations
  integer, parameter::Lo2_=nVar+1, Hi2_=nVar+nVar, Lo3_=Hi2_+1, Hi3_=nVar+Hi2_

  real, dimension(1:nVar,0:MaxIJK+1) :: dVarLim_VI
  real, dimension(1:nVar,-1:MaxIJK+2):: Primitive_VI
  logical, dimension(-1:MaxIJK+2)    :: IsTrueCell_I
end module ModLimiter

subroutine calc_facevalues(DoResChangeOnly)

  ! The subroutine calculates right and left face values.
  ! If DoResChangeOnly is true, only facevalues next to a coarser 
  ! neighbor block are calculated.
  ! To improve the code stability, ONLY those values are calculated, which
  ! are actually used by the calc_facefluxes.

  use ModMain
  use ModVarIndexes
  use ModLimiter
  use ModGeometry, ONLY : true_cell,body_BLK
  use ModNumConst
  use ModPhysics, ONLY : c2LIGHT,inv_c2LIGHT  !^CFG IF BORISCORR BEGIN
  use ModAdvance, ONLY:   B0xCell_BLK,B0yCell_BLK,B0zCell_BLK,&
       B0xFace_x_BLK,B0yFace_x_BLK,B0zFace_x_BLK,&
       B0xFace_y_BLK,B0yFace_y_BLK,B0zFace_y_BLK,&
       B0xFace_z_BLK,B0yFace_z_BLK,B0zFace_z_BLK    !^CFG END BORISCORR
  use ModAdvance, ONLY: State_VGB,&
       LeftState_VX,      &  ! Face Left  X
       RightState_VX,      &  ! Face Right X
       LeftState_VY,      &  ! Face Left  Y
       RightState_VY,      &  ! Face Right Y
       LeftState_VZ,      &  ! Face Left  Z
       RightState_VZ          ! Face Right Z

  use ModParallel, ONLY : neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  implicit none
  logical, intent(in):: DoResChangeOnly
  logical::DoTest,DoTestMe
  integer:: i,j,k,m,l
  real:: RhoInv

  real:: RhoC2Inv, BxFull, ByFull, BzFull, B2Full,& !^CFG IF BORISCORR
         uBC2Inv,Ga2Boris                           !^CFG IF BORISCORR

  !primitive variables
  real, dimension(nVar,&
       -1:nI+2,-1:nJ+2,-1:nK+2)::Primitive_VG
  !---------------------------------------------------------------------------



  if(globalBLK==BLKtest .and. .not. DoResChangeOnly )then
     call set_oktest('calc_facevalues',DoTest,DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if



  ! first, calculate the CELL values for the variables to be limited
  ! for non-borus corrections they are: density, velocity, pressure
  ! for boris correction momentum is used instead of the velocity
  !
  if(boris_correction)then                     !^CFG IF BORISCORR BEGIN 
     do k=kMinFaceX,kMaxFaceX
        do j=jMinFaceX,jMaxFaceX
           do i=1-nOrder,nI+nOrder
              call calc_primitives_boris       !needed for x-faces
           end do
        end do
     end do
     do k=kMinFaceY,kMaxFaceY
        do i=iMinFaceY,iMaxFaceY
           do j=1-nOrder,jMinFaceX-1
              call calc_primitives_boris       ! additional calculations for 
           end do                              ! y -faces
           do j=jMaxFaceX+1,nJ+nOrder
              call calc_primitives_boris       ! additional calculations for 
           end do	                       ! y-faces
        end do
     end do
     do j=jMinFaceZ,jMaxFaceZ
        do i=iMinFaceZ,iMaxFaceZ
           do k=1-nOrder,kMinFaceX-1
              call calc_primitives_boris       ! additional calculations for 
           end do                              ! z-faces
           do k=kMaxFaceX+1,nK+nOrder
              call calc_primitives_boris       ! additional calculations for 
           end do	                       ! z-faces
        end do
     end do
  else                                         !^CFG END BORISCORR
     do k=kMinFaceX,kMaxFaceX
        do j=jMinFaceX,jMaxFaceX
           do i=1-nOrder,nI+nOrder
              call calc_primitives_MHD         !needed for x-faces
           end do
        end do
     end do
     do k=kMinFaceY,kMaxFaceY
        do i=iMinFaceY,iMaxFaceY
           do j=1-nOrder,jMinFaceX-1
              call calc_primitives_MHD         ! additional calculations for 
           end do                              ! y -faces
           do j=jMaxFaceX+1,nJ+nOrder
              call calc_primitives_MHD         ! additional calculations for 
           end do	                       ! y-faces
        end do
     end do
     do j=jMinFaceZ,jMaxFaceZ
        do i=iMinFaceZ,iMaxFaceZ
           do k=1-nOrder,kMinFaceX-1
              call calc_primitives_MHD         ! additional calculations for 
           end do                              ! z-faces
           do k=kMaxFaceX+1,nK+nOrder
              call calc_primitives_MHD         ! additional calculations for 
           end do	                       ! z-faces
        end do
     end do
  end if                                       !^CFG IF BORISCORR

  !\
  ! write(*,*)' primitive variables have been calculated'
  !/
  !
  ! Now the first or second order face values are calcuted
  !
  select case(nOrder)
  case(1)
     !\
     ! First order reconstruction
     !/
     if (.not.DoResChangeOnly) then
        call get_faceX_first(1,nIFace,jMinFaceX,jMaxFaceX,kMinFaceX,kMaxFaceX)
        call get_faceY_first(iMinFaceY,iMaxFaceY,1,nJFace,kMinFaceY,kMaxFaceY)
        call get_faceZ_first(iMinFaceZ,iMaxFaceZ,jMinFaceZ,jMaxFaceZ,1,nKFace)
     else
        if(neiLeast(globalBLK)==+1)&
             call get_faceX_first(1,1,1,nJ,1,nK)
        if(neiLwest(globalBLK)==+1)&
             call get_faceX_first(nIFace,nIFace,1,nJ,1,nK)
        if(neiLsouth(globalBLK)==+1) &
             call get_faceY_first(1,nI,1,1,1,nK)
        if(neiLnorth(globalBLK)==+1) &
             call get_faceY_first(1,nI,nJFace,nJFace,1,nK)
        if(neiLbot(globalBLK)==+1) &
             call get_faceZ_first(1,nI,1,nJ,1,1)
        if(neiLtop(globalBLK)==+1) &
             call get_faceZ_first(1,nI,1,nJ,nKFace,nKFace)
     end if
  case(2)
     ! For second order scheme (nOrder==2)   
     ! use second order limited reconstruction.
     ! When prolong_order==1 we can use first order reconstruction 
     ! at resolution changes.
     ! However, constrained transport requires facevalues !^CFG IF CONSTRAINB
     ! to be independent of the resolution changes.       !^CFG IF CONSTRAINB

     if (.not.DoResChangeOnly)then
        call get_faceX_second(1,nIFace,jMinFaceX,jMaxFaceX,kMinFaceX,kMaxFaceX)
        call get_faceY_second(iMinFaceY,iMaxFaceY,1,nJFace,kMinFaceY,kMaxFaceY)
        call get_faceZ_second(iMinFaceZ,iMaxFaceZ,jMinFaceZ,jMaxFaceZ,1,nKFace)
     end if
     if(prolong_order==1 &
          .and..not.UseConstrainB & !^CFG IF CONSTRAINB
        )then
        ! First order facevalues at resolution change
        if(neiLeast(globalBLK)==+1)&
             call get_faceX_first(1,1,1,nJ,1,nK)
        if(neiLwest(globalBLK)==+1)&
             call get_faceX_first(nIFace,nIFace,1,nJ,1,nK)
        if(neiLsouth(globalBLK)==+1) &
             call get_faceY_first(1,nI,1,1,1,nK)
        if(neiLnorth(globalBLK)==+1) &
             call get_faceY_first(1,nI,nJFace,nJFace,1,nK)
        if(neiLbot(globalBLK)==+1) &
             call get_faceZ_first(1,nI,1,nJ,1,1)
        if(neiLtop(globalBLK)==+1) &
             call get_faceZ_first(1,nI,1,nJ,nKFace,nKFace)
     else if(DoResChangeOnly) then
        ! Second order face values at resolution changes
        if(neiLeast(globalBLK)==+1)&
             call get_faceX_second(1,1,1,nJ,1,nK)
        if(neiLwest(globalBLK)==+1)&
             call get_faceX_second(nIFace,nIFace,1,nJ,1,nK)
        if(neiLsouth(globalBLK)==+1) &
             call get_faceY_second(1,nI,1,1,1,nK)
        if(neiLnorth(globalBLK)==+1) &
             call get_faceY_second(1,nI,nJFace,nJFace,1,nK)
        if(neiLbot(globalBLK)==+1) &
             call get_faceZ_second(1,nI,1,nJ,1,1)
        if(neiLtop(globalBLK)==+1) &
             call get_faceZ_second(1,nI,1,nJ,nKFace,nKFace)
     endif
  end select  !end second order

  if(DoTestMe)write(*,*)&
       'calc_facevalues_opt: final UyFaceR_y=',RightState_VY(Uy_,Itest,Jtest,Ktest)

contains
  !========================================================
  subroutine calc_primitives_MHD
    Primitive_VG(:,i,j,k) = &
         State_VGB(1:nVar,i,j,k,globalBLK)
    RhoInv=cOne/Primitive_VG(rho_,i,j,k)
    Primitive_VG(Ux_:Uz_,i,j,k)=RhoInv*&
         Primitive_VG(rhoUx_:rhoUz_,i,j,k)
  end subroutine calc_primitives_MHD
  !======================================================
  !^CFG IF BORISCORR BEGIN
  subroutine calc_primitives_boris                                       
    !momentum is limited

    ! rhoU_Boris = rhoU - ((U x B) x B)/c^2
    !            = rhoU + (U B^2 - B U.B)/c^2
    !            = rhoU*(1+BB/(rho*c2)) - B UdotB/c^2
    Primitive_VG(:,i,j,k) = &
         State_VGB(1:nVar,i,j,k,globalBLK)
    BxFull = B0xCell_BLK(i,j,k,globalBLK) + Primitive_VG(Bx_,i,j,k)
    ByFull = B0yCell_BLK(i,j,k,globalBLK) + Primitive_VG(By_,i,j,k)
    BzFull = B0zCell_BLK(i,j,k,globalBLK) + Primitive_VG(Bz_,i,j,k)
    B2Full = BxFull**2 + ByFull**2 + BzFull**2
    RhoC2Inv  =cOne/(Primitive_VG(rho_,i,j,k)*c2LIGHT)
    uBC2Inv= (Primitive_VG(rhoUx_,i,j,k)*BxFull + &
         Primitive_VG(rhoUy_,i,j,k)*ByFull + &
         Primitive_VG(rhoUz_,i,j,k)*BzFull)*RhoC2Inv
    Ga2Boris=cOne+B2Full*RhoC2Inv

    Primitive_VG(Ux_,i,j,k)= Primitive_VG(rhoUx_,i,j,k)*&
         Ga2Boris - BxFull*uBC2Inv
    Primitive_VG(Uy_,i,j,k)= Primitive_VG(rhoUy_,i,j,k)*&
         Ga2Boris - ByFull*uBC2Inv
    Primitive_VG(Uz_,i,j,k)= Primitive_VG(rhoUz_,i,j,k)*&
         Ga2Boris - BzFull*uBC2Inv
  end subroutine calc_primitives_boris
  !======================================================                  
  !^CFG END BORISCORR
  subroutine get_faceX_first(iMin,iMax,jMin,jMax,kMin,kMax)
    integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
       LeftState_VX(:,i,j,k)=Primitive_VG(:,i-1,j,k)
       RightState_VX(:,i,j,k)=Primitive_VG(:,i,j,k)              
    end do; end do; end do
    if(boris_correction)&                                    !^CFG IF BORISCORR
         call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) !^CFG IF BORISCORR
  end subroutine get_faceX_first
  !========================================================================
  subroutine get_faceY_first(iMin,iMax,jMin,jMax,kMin,kMax)
    integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
       LeftState_VY(:,i,j,k)=Primitive_VG(:,i,j-1,k)
       RightState_VY(:,i,j,k)=Primitive_VG(:,i,j,k)              
    end do; end do; end do
    if(boris_correction)&                                    !^CFG IF BORISCORR
         call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) !^CFG IF BORISCORR
  end subroutine get_faceY_first
  !========================================================================
  subroutine get_faceZ_first(iMin,iMax,jMin,jMax,kMin,kMax)
    integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax
       LeftState_VZ(:,i,j,k)=Primitive_VG(:,i,j,k-1)
       RightState_VZ(:,i,j,k)=Primitive_VG(:,i,j,k)              
    end do; end do; end do
    if(boris_correction)&                                    !^CFG IF BORISCORR
         call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) !^CFG IF BORISCORR
  end subroutine get_faceZ_first
  !============================================================================
  !^CFG IF BORISCORR BEGIN
  subroutine BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)

    ! Convert face centered Boris momenta to MHD velocities

    integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
    !--------------------------------------------------------------------------
    ! U_Boris=rhoU_Boris/rho
    ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

       ! Left face values
       RhoInv=cOne/LeftState_VX(rho_,i,j,k)
       BxFull = B0xFace_x_BLK(i,j,k,globalBLK) + LeftState_VX(Bx_,i,j,k)
       ByFull = B0yFace_x_BLK(i,j,k,globalBLK) + LeftState_VX(By_,i,j,k)
       BzFull = B0zFace_x_BLK(i,j,k,globalBLK) + LeftState_VX(Bz_,i,j,k)
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  = inv_c2LIGHT*RhoInv
       LeftState_VX(Ux_,i,j,k)=LeftState_VX(Ux_,i,j,k)*RhoInv
       LeftState_VX(Uy_,i,j,k)=LeftState_VX(Uy_,i,j,k)*RhoInv
       LeftState_VX(Uz_,i,j,k)=LeftState_VX(Uz_,i,j,k)*RhoInv
       uBC2Inv= (LeftState_VX(Ux_,i,j,k)*BxFull + &
                 LeftState_VX(Uy_,i,j,k)*ByFull + &
                 LeftState_VX(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris=cOne/(cOne+ B2Full*RhoC2Inv)

       LeftState_VX(Ux_,i,j,k) = Ga2Boris * (LeftState_VX(Ux_,i,j,k)+uBC2Inv*BxFull)
       LeftState_VX(Uy_,i,j,k) = Ga2Boris * (LeftState_VX(Uy_,i,j,k)+uBC2Inv*ByFull)
       LeftState_VX(Uz_,i,j,k) = Ga2Boris * (LeftState_VX(Uz_,i,j,k)+uBC2Inv*BzFull)

       ! Right face values
       RhoInv=cOne/RightState_VX(rho_,i,j,k)
       BxFull = B0xFace_x_BLK(i,j,k,globalBLK) + RightState_VX(Bx_,i,j,k)
       ByFull = B0yFace_x_BLK(i,j,k,globalBLK) + RightState_VX(By_,i,j,k)
       BzFull = B0zFace_x_BLK(i,j,k,globalBLK) + RightState_VX(Bz_,i,j,k)
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  =inv_c2LIGHT*RhoInv
       RightState_VX(Ux_,i,j,k)=RightState_VX(Ux_,i,j,k)*RhoInv
       RightState_VX(Uy_,i,j,k)=RightState_VX(Uy_,i,j,k)*RhoInv
       RightState_VX(Uz_,i,j,k)=RightState_VX(Uz_,i,j,k)*RhoInv
       uBC2Inv= (RightState_VX(Ux_,i,j,k)*BxFull + &
                 RightState_VX(Uy_,i,j,k)*ByFull + &
                 RightState_VX(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris=cOne/(cOne+B2Full*RhoC2Inv)

       RightState_VX(Ux_,i,j,k) = Ga2Boris * (RightState_VX(Ux_,i,j,k)+uBC2Inv*BxFull)
       RightState_VX(Uy_,i,j,k) = Ga2Boris * (RightState_VX(Uy_,i,j,k)+uBC2Inv*ByFull)
       RightState_VX(Uz_,i,j,k) = Ga2Boris * (RightState_VX(Uz_,i,j,k)+uBC2Inv*BzFull)

    end do; end do; end do

  end subroutine BorisFaceXtoMHD
  !======================================================

  subroutine BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)
    integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax

    ! U_Boris=rhoU_Boris/rho
    ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

       ! Left face values
       RhoInv=cOne/LeftState_VY(rho_,i,j,k)
       BxFull = B0xFace_y_BLK(i,j,k,globalBLK) + LeftState_VY(Bx_,i,j,k)
       ByFull = B0yFace_y_BLK(i,j,k,globalBLK) + LeftState_VY(By_,i,j,k)
       BzFull = B0zFace_y_BLK(i,j,k,globalBLK) + LeftState_VY(Bz_,i,j,k)
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  =inv_c2LIGHT*RhoInv
       LeftState_VY(Ux_,i,j,k)=LeftState_VY(Ux_,i,j,k)*RhoInv
       LeftState_VY(Uy_,i,j,k)=LeftState_VY(Uy_,i,j,k)*RhoInv
       LeftState_VY(Uz_,i,j,k)=LeftState_VY(Uz_,i,j,k)*RhoInv
       uBC2Inv= (LeftState_VY(Ux_,i,j,k)*BxFull + &
                 LeftState_VY(Uy_,i,j,k)*ByFull + &
                 LeftState_VY(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris=cOne/(cOne+ B2Full*RhoC2Inv)

        LeftState_VY(Ux_,i,j,k) = Ga2Boris * (LeftState_VY(Ux_,i,j,k)+uBC2Inv*BxFull)
        LeftState_VY(Uy_,i,j,k) = Ga2Boris * (LeftState_VY(Uy_,i,j,k)+uBC2Inv*ByFull)
        LeftState_VY(Uz_,i,j,k) = Ga2Boris * (LeftState_VY(Uz_,i,j,k)+uBC2Inv*BzFull)

       ! Right face values
       RhoInv=cOne/RightState_VY(rho_,i,j,k) 
       BxFull = B0xFace_y_BLK(i,j,k,globalBLK) + RightState_VY(Bx_,i,j,k)
       ByFull = B0yFace_y_BLK(i,j,k,globalBLK) + RightState_VY(By_,i,j,k)
       BzFull = B0zFace_y_BLK(i,j,k,globalBLK) + RightState_VY(Bz_,i,j,k)
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv=inv_c2LIGHT*RhoInv
       RightState_VY(Ux_,i,j,k)=RightState_VY(Ux_,i,j,k)*RhoInv
       RightState_VY(Uy_,i,j,k)=RightState_VY(Uy_,i,j,k)*RhoInv
       RightState_VY(Uz_,i,j,k)=RightState_VY(Uz_,i,j,k)*RhoInv
       uBC2Inv= (RightState_VY(Ux_,i,j,k)*BxFull + &
                 RightState_VY(Uy_,i,j,k)*ByFull + &
                 RightState_VY(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris=cOne/(cOne + B2Full*RhoC2Inv)

       RightState_VY(Ux_,i,j,k) = Ga2Boris * (RightState_VY(Ux_,i,j,k)+uBC2Inv*BxFull)
       RightState_VY(Uy_,i,j,k) = Ga2Boris * (RightState_VY(Uy_,i,j,k)+uBC2Inv*ByFull)
       RightState_VY(Uz_,i,j,k) = Ga2Boris * (RightState_VY(Uz_,i,j,k)+uBC2Inv*BzFull)
    end do; end do; end do

  end subroutine BorisFaceYtoMHD
  !============================================================================
  subroutine BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax)
    integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax

    ! Convert face centered Boris momenta/rho to MHD velocities
    !--------------------------------------------------------------------------
    ! U_Boris=rhoU_Boris/rho
    ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

       ! Left face values
       RhoInv=cOne/LeftState_VZ(rho_,i,j,k)
       BxFull = B0xFace_z_BLK(i,j,k,globalBLK) + LeftState_VZ(Bx_,i,j,k)
       ByFull = B0yFace_z_BLK(i,j,k,globalBLK) + LeftState_VZ(By_,i,j,k)
       BzFull = B0zFace_z_BLK(i,j,k,globalBLK) + LeftState_VZ(Bz_,i,j,k)
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  =inv_c2LIGHT*RhoInv
       LeftState_VZ(Ux_,i,j,k)=LeftState_VZ(Ux_,i,j,k)*RhoInv
       LeftState_VZ(Uy_,i,j,k)=LeftState_VZ(Uy_,i,j,k)*RhoInv
       LeftState_VZ(Uz_,i,j,k)=LeftState_VZ(Uz_,i,j,k)*RhoInv
       uBC2Inv= (LeftState_VZ(Ux_,i,j,k)*BxFull + &
                 LeftState_VZ(Uy_,i,j,k)*ByFull + &
                 LeftState_VZ(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris=cOne/(cOne + B2Full*RhoC2Inv)

       LeftState_VZ(Ux_,i,j,k) = Ga2Boris * (LeftState_VZ(Ux_,i,j,k)+uBC2Inv*BxFull)
       LeftState_VZ(Uy_,i,j,k) = Ga2Boris * (LeftState_VZ(Uy_,i,j,k)+uBC2Inv*ByFull)
       LeftState_VZ(Uz_,i,j,k) = Ga2Boris * (LeftState_VZ(Uz_,i,j,k)+uBC2Inv*BzFull)

       ! Right face values
       RhoInv=cOne/RightState_VZ(rho_,i,j,k)
       BxFull = B0xFace_z_BLK(i,j,k,globalBLK) + RightState_VZ(Bx_,i,j,k)
       ByFull = B0yFace_z_BLK(i,j,k,globalBLK) + RightState_VZ(By_,i,j,k)
       BzFull = B0zFace_z_BLK(i,j,k,globalBLK) + RightState_VZ(Bz_,i,j,k)
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  =inv_c2LIGHT*RhoInv
       RightState_VZ(Ux_,i,j,k)=RightState_VZ(Ux_,i,j,k)*RhoInv
       RightState_VZ(Uy_,i,j,k)=RightState_VZ(Uy_,i,j,k)*RhoInv
       RightState_VZ(Uz_,i,j,k)=RightState_VZ(Uz_,i,j,k)*RhoInv
       uBC2Inv= (RightState_VZ(Ux_,i,j,k)*BxFull + &
                 RightState_VZ(Uy_,i,j,k)*ByFull + &
                 RightState_VZ(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris=cOne/(cOne + B2Full*RhoC2Inv)

       RightState_VZ(Ux_,i,j,k) = Ga2Boris * (RightState_VZ(Ux_,i,j,k)+uBC2Inv*BxFull)
       RightState_VZ(Uy_,i,j,k) = Ga2Boris * (RightState_VZ(Uy_,i,j,k)+uBC2Inv*ByFull)
       RightState_VZ(Uz_,i,j,k) = Ga2Boris * (RightState_VZ(Uz_,i,j,k)+uBC2Inv*BzFull)
    end do; end do; end do

  end subroutine BorisFaceZtoMHD
  !============================================================================
  !^CFG END BORISCORR

  subroutine get_faceX_second(iMin,iMax,jMin,jMax,kMin,kMax)
    integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
    integer::i1
    do k=kMin, kMax; do j=jMin, jMax; 
       Primitive_VI(:,iMin-2:iMax+1)=Primitive_VG(:,iMin-2:iMax+1,j,k)
       if(body_BLK(globalBLK))then
          IsTrueCell_I(iMin-2:iMax+1)=&
               true_cell(iMin-2:iMax+1,j,k,globalBLK)
          call limiter_body(iMin,iMax)
       else
          call limiter(iMin,iMax)
       end if
       do i=iMin,iMax
          i1=i-1
          LeftState_VX(:,i,j,k) =Primitive_VI(:,i1)+dVarLim_VI(:,i1)
          RightState_VX(:,i,j,k) =Primitive_VI(:,i )-dVarLim_VI(:,i )

       end do
    end do; end do
    if(boris_correction) &                                   !^CFG IF BORISCORR
         call BorisFaceXtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) !^CFG IF BORISCORR
  end subroutine get_faceX_second
  !============================================================================
  subroutine get_faceY_second(iMin,iMax,jMin,jMax,kMin,kMax)
    integer,intent(in):: iMin,iMax,jMin,jMax,kMin,kMax
    integer::j1
    do k=kMin, kMax; do i=iMin,iMax
       Primitive_VI(:,jMin-2:jMax+1)=Primitive_VG(:,i,jMin-2:jMax+1,k)
       if(body_BLK(globalBLK))then
          IsTrueCell_I(jMin-2:jMax+1)=&
               true_cell(i,jMin-2:jMax+1,k,globalBLK)
          call limiter_body(jMin,jMax)
       else
          call limiter(jMin,jMax)
       end if
       do j=jMin, jMax
          j1=j-1
          LeftState_VY(:,i,j,k) =Primitive_VI(:,j1)+dVarLim_VI(:,j1)
          RightState_VY(:,i,j,k) =Primitive_VI(:,j )-dVarLim_VI(:,j )
       end do
    end do; end do
    if(boris_correction) &                                   !^CFG IF BORISCORR
         call BorisFaceYtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) !^CFG IF BORISCORR
  end subroutine get_faceY_second
  !============================================================================
  subroutine get_faceZ_second(iMin,iMax,jMin,jMax,kMin,kMax)
    integer,intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax
    integer::k1
    do j=jMin,jMax; do i=iMin,iMax; 
       Primitive_VI(:,kMin-2:kMax+1)=Primitive_VG(:,i,j,kMin-2:kMax+1)
       if(body_BLK(globalBLK))then
          IsTrueCell_I(kMin-2:kMax+1)=&
               true_cell(i,j,kMin-2:kMax+1,globalBLK)
          call limiter_body(kMin,kMax)
       else
          call limiter(kMin,kMax)
       end if
       do k=kMin,kMax
          k1=k-1
          LeftState_VZ(:,i,j,k) =Primitive_VI(:,k1)+dVarLim_VI(:,k1)
          RightState_VZ(:,i,j,k) =Primitive_VI(:,k )-dVarLim_VI(:,k )
       end do
    end do; end do
    if(boris_correction)&                                    !^CFG IF BORISCORR
         call BorisFaceZtoMHD(iMin,iMax,jMin,jMax,kMin,kMax) !^CFG IF BORISCORR
  end subroutine get_faceZ_second
end subroutine calc_facevalues
!============================================================================
!
! beta-limiter is implemented (see C.Hirsch, Numerical Computation of
! internal and external flows, V2,page.544-545.
! In particular case beta=1.0 (default)
! this limiter coincides with minmod.
!
!=============================================
subroutine limiter_body(lMin,lMax)
  use ModLimiter
  use ModNumConst
  use ModMain, ONLY: limiter_type, BetaLimiter !^CFG IF NOT SIMPLE
  implicit none
  integer, intent(in)::lMin,lMax
  real,dimension(Hi3_):: dVar2_I,dVar1_I
  integer::l

  select case(limiter_type)                   !^CFG IF NOT SIMPLE BEGIN
  case('beta')
     dVar1_I(1:nVar)=Primitive_VI(:,lMax+1)-Primitive_VI(:,lMax)
     dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
     dVar1_I(Lo3_:Hi3_)=BetaLimiter*&
          dVar1_I(Lo2_:Hi2_)
     dVar1_I(1:nVar)=sign(cQuarter,dVar1_I(1:nVar))
     do l=lMax,lMin-1,-1
        dVar2_I=dVar1_I
        dVar1_I(1:nVar)=Primitive_VI(:,l)-Primitive_VI(:,l-1)
        dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
        dVar1_I(Lo3_:Hi3_)=BetaLimiter*&
             dVar1_I(Lo2_:Hi2_)
        dVar1_I(1:nVar)=sign(cQuarter,dVar1_I(1:nVar))
        if(all(IsTrueCell_I(l-1:l+1)))then
           dVar2_I(1:nVar)=dVar2_I(1:nVar)+dVar1_I(1:nVar)
           dVar2_I(Lo2_:Hi2_)=min(dVar2_I(Lo2_:Hi2_),dVar1_I(Lo3_:Hi3_))
           dVar2_I(Lo3_:Hi3_)=min(dVar1_I(Lo2_:Hi2_),dVar2_I(Lo3_:Hi3_))
           dVar2_I(Lo2_:Hi2_)=max(dVar2_I(Lo2_:Hi2_),dVar2_I(Lo3_:Hi3_))
           dVarLim_VI(:,l)=dVar2_I(1:nVar)* dVar2_I(Lo2_:Hi2_)
        else
           dVarLim_VI(:,l)=cZero
        end if
     end do
  case('minmod')                                      !^CFG END SIMPLE
     dVar1_I(1:nVar)=Primitive_VI(:,lMax+1)-Primitive_VI(:,lMax)
     dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
     dVar1_I(1:nVar)=sign(cQuarter,dVar1_I(1:nVar))
     do l=lMax,lMin-1,-1
        dVar2_I(1:Hi2_)=dVar1_I(1:Hi2_)
        dVar1_I(1:nVar)=Primitive_VI(:,l)-Primitive_VI(:,l-1)
        dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
        dVar1_I(1:nVar)=sign(cQuarter,dVar1_I(1:nVar))
        if(all(IsTrueCell_I(l-1:l+1)))then
           dVar2_I(1:nVar)=dVar2_I(1:nVar)+dVar1_I(1:nVar)
           dVar2_I(Lo2_:Hi2_)=min(dVar2_I(Lo2_:Hi2_),dVar1_I(Lo2_:Hi2_))
           dVarLim_VI(:,l)=dVar2_I(1:nVar)* dVar2_I(Lo2_:Hi2_)
        else
           dVarLim_VI(:,l)=cZero
        end if
     end do
  case default
     call stop_mpi('limiter_body: unknown TypeLimiter='//limiter_type)
  end select

end subroutine limiter_body
!===========================================================================
subroutine limiter(lMin,lMax)
  use ModLimiter
  use ModNumConst
  use ModMain, ONLY: limiter_type, BetaLimiter !^CFG IF NOT SIMPLE
  implicit none
  integer, intent(in)::lMin,lMax
  real,dimension(Hi3_):: dVar2_I,dVar1_I
  integer::l

  select case(limiter_type)                  !^CFG IF NOT SIMPLE BEGIN
  case('beta')
     dVar1_I(1:nVar)=Primitive_VI(:,lMax+1)-Primitive_VI(:,lMax)
     dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
     dVar1_I(Lo3_:Hi3_)=BetaLimiter*&
          dVar1_I(Lo2_:Hi2_)
     dVar1_I(1:nVar)=sign(cQuarter,dVar1_I(1:nVar))
     do l=lMax,lMin-1,-1
        dVar2_I=dVar1_I
        dVar1_I(1:nVar)=Primitive_VI(:,l)-Primitive_VI(:,l-1)
        dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
        dVar1_I(Lo3_:Hi3_)=BetaLimiter*dVar1_I(Lo2_:Hi2_)
        dVar1_I(1:nVar)=sign(cQuarter,dVar1_I(1:nVar))
        dVar2_I(1:nVar)=dVar2_I(1:nVar)+dVar1_I(1:nVar)
        dVar2_I(Lo2_:Hi2_)=min(dVar2_I(Lo2_:Hi2_),dVar1_I(Lo3_:Hi3_))
        dVar2_I(Lo3_:Hi3_)=min(dVar1_I(Lo2_:Hi2_),dVar2_I(Lo3_:Hi3_))
        dVar2_I(Lo2_:Hi2_)=max(dVar2_I(Lo2_:Hi2_),dVar2_I(Lo3_:Hi3_))
        dVarLim_VI(:,l)=dVar2_I(1:nVar)*dVar2_I(Lo2_:Hi2_)
     end do
  case('minmod')                             !^CFG END SIMPLE
     dVar1_I(1:nVar)=Primitive_VI(:,lMax+1)-Primitive_VI(:,lMax)
     dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
     dVar1_I(1:nVar)=sign(cQuarter,dVar1_I(1:nVar))
     do l=lMax,lMin-1,-1
        dVar2_I(1:Hi2_)=dVar1_I(1:Hi2_)
        dVar1_I(1:nVar)=Primitive_VI(:,l)-Primitive_VI(:,l-1)
        dVar1_I(Lo2_:Hi2_)=abs(dVar1_I(1:nVar))
        dVar1_I(1:nVar)=sign(cQuarter,dVar1_I(1:nVar))
        dVar2_I(1:nVar)=dVar2_I(1:nVar)+dVar1_I(1:nVar)
        dVar2_I(Lo2_:Hi2_)=min(dVar2_I(Lo2_:Hi2_),dVar1_I(Lo2_:Hi2_))
        dVarLim_VI(:,l)=dVar2_I(1:nVar)* dVar2_I(Lo2_:Hi2_)
     end do
  case default
     call stop_mpi('limiter: unknown TypeLimiter='//limiter_type)
  end select

end subroutine limiter

