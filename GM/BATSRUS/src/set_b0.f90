!^CFG COPYRIGHT UM
!---------------------------------------------------------------------------
!                     B0 module calculates the B0 field
!---------------------------------------------------------------------------

!\
! These routines compute the intrinsic magnetic field B0 for each of 
! the cells on the grid.  
! The cell center as well as the face center values are calculated.
!/

!============================================================================
subroutine set_b0_cell(iBlock)


  ! Cell center values are calculated through all the block

  use ModProcMH, ONLY: iProc
  use ModMain,   ONLY: ProcTest, BlkTest, iTest, jTest, kTest
  use ModB0,     ONLY: B0_DGB
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK
  use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

  implicit none

  integer, intent(in) :: iBlock

  integer :: i, j, k

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  if(iProc==PROCtest.and.iBlock==BLKtest)then
     call set_oktest('set_b0_cell',DoTest,DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  endif

  do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
     call get_b0(x_BLK(i,j,k,iBlock),&
          y_BLK(i,j,k,iBlock),&
          z_BLK(i,j,k,iBlock),&
          B0_DGB(:,i,j,k,iBlock))
  end do; end do; end do

  if(DoTestMe)write(*,*)'B0*Cell_BLK=',&
       B0_DGB(:,Itest,Jtest,Ktest,BLKtest)

end subroutine set_b0_cell

!============================================================================

subroutine get_b0(X0,Y0,Z0,B0)

  ! this interface allows use of various models for B0
  use ModMain,          ONLY : Time_Simulation, NameThisComp, &
       TypeCoordSystem, IsStandAlone
  use ModPhysics,       ONLY : Si2No_V, UnitB_, DipoleStrengthSi
  use CON_planet_field, ONLY : get_planet_field
  use ModMain,          ONLY : UseBody2             !^CFG IF SECONDBODY
  use ModMain,          ONLY : UseUserB0, UseMagnetogram
  use ModUser,          ONLY: user_get_b0
  implicit none

  real, intent(in) :: X0,Y0,Z0
  real, intent(out), dimension(3) :: B0

  if(UseMagnetogram)then
     call get_coronal_b0(X0,Y0,Z0,B0)
  elseif(UseUserB0)then
     call user_get_b0(X0,Y0,Z0,B0)
  elseif(IsStandAlone .and. DipoleStrengthSi==0.0)then
     B0 = 0.0
     RETURN
  elseif(NameThisComp=='GM')then
     call get_planet_field(Time_Simulation,X0,Y0,Z0,TypeCoordSystem//' NORM',&
          B0)
     B0 = B0*Si2No_V(UnitB_)
  else
     call get_b0_multipole(X0,Y0,Z0,B0) 
  end if
  if(UseBody2)call add_b0_body2(X0,Y0,Z0,B0)        !^CFG IF SECONDBODY

end subroutine get_b0

!============================================================================

subroutine get_b0_multipole(X0,Y0,Z0,B0)

  !\
  ! This routine returns the intrinsic magnetic field using a multipole
  ! expansion at the given location (X0,Y0,Z0).
  !
  ! Any multiple expansion has an inner and an outer boundary condition
  ! For many applications an outer boundary of B -> 0 at as r-> infinity
  ! In some sense this multiple expansion considers only constraints
  ! imposed by the inner boundary. 
  ! 
  ! Another option commonly used treating the solar magnetic field is to 
  ! introduce a radial boundary known as the source surface where the imposed 
  ! boundary condition  requires a radial magnetic field. This is 
  ! known as a source surface model and is denoted here as "B0SourceSurface"
  !/

  use ModPhysics
  use ModNumConst
  implicit none

  real, intent(in) :: X0,Y0,Z0
  real, intent(out), dimension(3) :: B0

  integer :: i, j, k, l
  real :: R0, rr, rr_inv, rr2_inv, rr3_inv, rr5_inv, rr7_inv
  real, dimension(3) :: xxt, bb
  real :: Dp, temp1, temp2

  logical :: do_quadrupole, do_octupole
  !-------------------------------------------------------------------------
  !\
  ! Determine radial distance and powers of it
  !/

  R0 = sqrt(X0*X0 + Y0*Y0 + Z0*Z0)

  ! Avoid calculating B0 inside a critical radius = 1.E-6*Rbody
  if(R0 <= cTiny*Rbody)then
     B0=0.0
     RETURN
  end if

  rr     = R0
  rr_inv = cOne/rr
  rr2_inv=rr_inv*rr_inv
  rr3_inv=rr_inv*rr2_inv

  !\
  ! Compute dipole moment of the intrinsic magnetic field B0.
  !/

  Bdpx = -sinTHETAtilt*Bdp       ! 1
  Bdpy = cZero                   ! 2
  Bdpz = cosTHETAtilt*Bdp        ! 3 

  Dp = (Bdpx*X0+Bdpy*Y0+Bdpz*Z0)*3.00*rr2_inv

  B0(1) = (Dp*X0-Bdpx)*rr3_inv
  B0(2) = (Dp*Y0-Bdpy)*rr3_inv
  B0(3) = (Dp*Z0-Bdpz)*rr3_inv

  do_quadrupole=any(abs(Qqp)>cTiny)
  do_octupole  =any(abs(Oop)>cTiny)

  if(do_quadrupole.or.do_octupole)then
     !\
     ! Compute the xx's in the tilted reference frame aligned with
     ! the magnetic field.
     !/
     xxt(1) = cosTHETAtilt*X0 + sinTHETAtilt*Z0
     xxt(2) = Y0
     xxt(3)= -sinTHETAtilt*X0 + cosTHETAtilt*Z0

     rr5_inv = rr3_inv*rr2_inv
     rr7_inv = rr5_inv*rr2_inv
  end if

  if(do_quadrupole)then
     !\
     ! Compute quadrupole moment of the intrinsic 
     ! magnetic field B0.
     !/
     do k=1,3
        temp1 = cZero
        temp2 = cZero
        do i=1,3
           temp1 = temp1 + Qqp(k,i)*xxt(i)
           do j=1,3
              temp2 = temp2 + Qqp(i,j)*xxt(i)*xxt(j)*xxt(k)
           end do
        end do
        bb(k) = 5.0*cHalf*temp2*rr7_inv - temp1*rr5_inv
     end do

     B0(1) = B0(1) + cosTHETAtilt*bb(1) - sinTHETAtilt*bb(3) 
     B0(2) = B0(2) + bb(2)
     B0(3) = B0(3) + sinTHETAtilt*bb(1) + cosTHETAtilt*bb(3)
  end if

  if(do_octupole)then
     !\
     ! Compute octupole moment of the intrinsic 
     ! magnetic field B0.
     !/
     do k = 1, 3
        temp1 = cZero
        temp2 = cZero
        do i = 1, 3
           do j = 1, 3
              temp1 = temp1 + Oop(i,j,k)*xxt(i)*xxt(j)
              do l = 1, 3
                 temp2 = temp2 + Oop(i,j,l)*xxt(i)*xxt(j)*xxt(l)*xxt(k)
              end do
           end do
        end do
        bb(k) = 7.0*temp2*rr7_inv - 3.0*temp1*rr5_inv
     end do

     B0(1) = B0(1) + cosTHETAtilt*bb(1) - sinTHETAtilt*bb(3) 
     B0(2) = B0(2) + bb(2)
     B0(3) = B0(3) + sinTHETAtilt*bb(1) + cosTHETAtilt*bb(3)
  end if

end subroutine get_b0_multipole

!^CFG IF SECONDBODY BEGIN
!=============================================================================
subroutine add_b0_body2(X0,Y0,Z0,B0)

  !\
  ! If there is a second body that has a magnetic field the contribution
  ! to the field from the second body should be computed here (inside the
  ! if block.
  !/
  use ModPhysics
  use ModNumConst
  implicit none

  real, intent(in)   :: X0,Y0,Z0
  real, intent(inout):: B0(3)

  real :: Xyz_D(3),R0,rr_inv,rr2_inv,rr3_inv,Dp
  !--------------------------------------------------------------------------
  !\
  ! Determine normalized relative coordinates and radial distance from body 2
  !/
  Xyz_D(1) = (X0-xBody2)/rBody2
  Xyz_D(2) = (Y0-yBody2)/rBody2
  Xyz_D(3) = (Z0-zBody2)/rBody2
  R0 = sqrt(sum(Xyz_D**2))

  ! Avoid calculating B0 inside a critical normalized radius = cTiny
  if(R0 <= cTiny) RETURN

  rr_inv = cOne/R0
  rr2_inv=rr_inv*rr_inv
  rr3_inv=rr_inv*rr2_inv

  ! Add dipole field of the second body
  Dp = sum(BdpBody2_D*Xyz_D)*3*rr2_inv

  B0 = B0 + (Dp*Xyz_D - BdpBody2_D)*rr3_inv

end subroutine add_b0_body2
!^CFG END SECONDBODY

!==============================================================================

subroutine update_b0

  use ModMain,          ONLY: nBlock, unusedBLK, &
       time_simulation, NameThisComp
  use ModPhysics,       ONLY: ThetaTilt
  use ModAdvance,       ONLY: Bx_, By_, Bz_, State_VGB, B0_DGB
  use ModGeometry,      ONLY: true_cell, body_BLK
  use CON_axes,         ONLY: get_axes
  use ModNumConst,      ONLY: cRadToDeg
  use ModIO,            ONLY: iUnitOut, write_prefix
  use ModEnergy,        ONLY: calc_energy_ghost
  use ModB0,            ONLY: set_b0_reschange

  implicit none

  character(len=*), parameter :: NameSub = 'update_b0'
  logical :: DoTest, DoTestMe
  integer :: iBlock
  !============================================================================

  call set_oktest(NameSub, DoTest, DoTestMe)

  ! Update ThetaTilt
  if(NameThisComp=='GM') &
       call get_axes(Time_Simulation, MagAxisTiltGsmOut=ThetaTilt)

  if (DoTestMe) then
     if(NameThisComp=='GM')then
        call write_prefix; write(iUnitOut,*) &
          "update_b0 at tSimulation, TiltGsm=", &
          Time_Simulation, ThetaTilt*cRadToDeg
     else
        call write_prefix; write(iUnitOut,*) &
             "update_b0 at tSimulation=",Time_Simulation
     end if
  end if
  call timing_start(NameSub)

  do iBlock=1,nBlock
     if(unusedBLK(iBlock)) CYCLE

     ! Save total magnetic field into Bx_BLK,By_BLK,Bz_BLK
     State_VGB(Bx_:Bz_,:,:,:,iBlock) = State_VGB(Bx_:Bz_,:,:,:,iBlock) &
          + B0_DGB(:,:,:,:,iBlock)

     call set_b0_cell(iBlock)

     ! Split total B again using new B0
     State_VGB(Bx_:Bz_,:,:,:,iBlock) = State_VGB(Bx_:Bz_,:,:,:,iBlock) &
          - B0_DGB(:,:,:,:,iBlock)
       
     ! Set B1 to 0 inside bodies
     if(Body_BLK(iBlock))then
        where(.not.true_cell(:,:,:,iBlock))
           State_VGB(Bx_,:,:,:,iBlock)=0.0
           State_VGB(By_,:,:,:,iBlock)=0.0
           State_VGB(Bz_,:,:,:,iBlock)=0.0
        end where
     end if

     ! Recalculate energy
     call calc_energy_ghost(iBlock)

  end do

  ! Recalculate B0 face values at resolution changes
  call set_b0_reschange

  call timing_stop(NameSub)

end subroutine update_b0

!===========================================================================
subroutine get_coronal_b0(xInput,yInput,zInput,B0_D)
  use ModPhysics,     ONLY:Si2No_V,UnitB_
  use ModMagnetogram, ONLY: get_magnetogram_field
  implicit none

  real, intent(in):: xInput,yInput,zInput
  real, intent(out), dimension(3):: B0_D
  !--------------------------------------------------------------------------
  call get_magnetogram_field(xInput,yInput,zInput,B0_D)
  B0_D = B0_D*Si2No_V(UnitB_)

end subroutine get_coronal_b0
