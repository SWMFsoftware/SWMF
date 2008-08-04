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
subroutine set_b0(iBlock)

  implicit none
  integer, intent(in) :: iBlock

  call set_b0_cell(iBlock)
  call set_b0_source(iBlock)

end subroutine set_b0
!============================================================================
subroutine set_b0_cell(iBlock)
  use ModProcMH
  use ModMain
  use ModB0
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK
  use ModNumConst
  implicit none

  integer, intent(in) :: iBlock

  integer :: i,j,k
  real :: B0_D(3)
  real :: x,y,z

  logical :: oktest, oktest_me
  

  !--------------------------------------------------------------------------
  if(iProc==PROCtest.and.iBlock==BLKtest)then
     call set_oktest('set_b0_cell',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  endif



  !\
  ! Cell center values are calculated through all the block
  !/
  do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
     call get_b0(x_BLK(i,j,k,iBlock),&
          y_BLK(i,j,k,iBlock),&
          z_BLK(i,j,k,iBlock),&
          B0_DGB(:,i,j,k,iBlock))
  end do; end do; end do

  if(oktest_me)write(*,*)'B0*Cell_BLK=',&
       B0_DGB(:,Itest,Jtest,Ktest,BLKtest)
end subroutine set_b0_cell

!============================================================================
subroutine set_b0_source(iBlock)
  use ModProcMH
  use ModMain
  use ModB0
  use ModGeometry, ONLY : &       
       dx_BLK,dy_BLK,dz_BLK,XyzStart_BLK
  use ModGeometry, ONLY : UseCovariant                
  use ModParallel, ONLY : neiLtop,neiLbot,neiLeast,neiLwest,neiLnorth,neiLsouth
  use ModNumConst
  implicit none

  integer, intent(in) :: iBlock

  integer :: i,j,k
  real,dimension(3):: B0_D,RefXyzStart_D,RefDXyz_D
  real,dimension(nDim,0:1,0:1,0:1)::RefB0_DIII
  real ::x,y,z
  ! inverse of Dx, Dy, Dz                              
  real:: DxInv, DyInv, DzInv                                   


  logical :: oktest, oktest_me

  !--------------------------------------------------------------------------
 
  if(iProc==PROCtest.and.iBlock==BLKtest)then
     call set_oktest('set_b0_face',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  endif


  if(UseCovariant)then                   
     call calc_b0source_covar(iBlock)         
     return
  end if                                 


  RefDXyz_D(x_)=cHalf*dx_BLK(iBlock)
  RefDXyz_D(y_)=cHalf*dy_BLK(iBlock)
  RefDXyz_D(z_)=cHalf*dz_BLK(iBlock)
  RefXyzStart_D(:)=XyzStart_BLK(:,iBlock)-cHalf*RefDXyz_D(:)
 
  if (neiLeast(iBlock)==-1) then
     if(oktest_me .and. Itest==1)then
        call set_b0_face(iBlock)
          write(*,*)'Before correcting east refinement: B0*Face_x_BLK=',&
          B0_DX(:,Itest,Jtest,Ktest)
     end if
     call correct_b0_face(East_)
  end if

  if(neiLwest(iBlock)==-1) then
     if(oktest_me .and. Itest==nIFace)then
        call set_b0_face(iBlock)
          write(*,*)'Before correcting west refinement: B0*Face_x_BLK=',&
          B0_DX(:,Itest,Jtest,Ktest)
     end if
     call correct_b0_face(West_)
  end if
  if (neiLsouth(iBlock)==-1) then
     if(oktest_me .and. jTest==1)then
        call set_b0_face(iBlock)
          write(*,*)'Before correcting south refinement: B0*Face_y_BLK=',&
          B0_DY(:,Itest,Jtest,Ktest)
     end if
     call correct_b0_face(South_)
  end if

  if(neiLnorth(iBlock)==-1) then
     if(oktest_me .and. jTest==nJFace)then
        call set_b0_face(iBlock)
          write(*,*)'Before correcting north refinement: B0*Face_y_BLK=',&
          B0_DY(:,Itest,Jtest,Ktest)
     end if
     call correct_b0_face(North_)
  end if

  if (neiLbot(iBlock)==-1) then
     if(oktest_me .and. kTest==1)then
        call set_b0_face(iBlock)
          write(*,*)'Before correcting bot refinement: B0*Face_z_BLK=',&
          B0_DZ(:,Itest,Jtest,Ktest)
     end if
     call correct_b0_face(Bot_)
  end if

  if(neiLtop(iBlock)==-1) then 

     if(oktest_me .and. kTest==nKFace)then
        call set_b0_face(iBlock)
          write(*,*)'Before correcting top refinement: B0*Face_z_BLK=',&
          B0_DZ(:,Itest,Jtest,Ktest)
     end if
     call correct_b0_face(Top_)
  end if

  if(oktest_me)then
     write(*,*)'B0*Face_x_BLK=',&
          B0_DX(:,Itest,Jtest,Ktest)
     write(*,*)'B0*Face_y_BLK=',&
          B0_DY(:,Itest,Jtest,Ktest)
     write(*,*)'B0*Face_z_BLK=',&
          B0_DZ(:,Itest,Jtest,Ktest)
  end if

  call set_b0_matrix(iBlock)
 

contains
  !===========================================================================
  subroutine correct_b0_face(iSide)
    implicit none
    integer,intent(in)::iSide
    integer::iFace,jFace,kFace,iDim
    !-------------------------------------------------------------------------
    select case(iSide)
    case(East_,West_)
       iFace=1+nI*(iSide-East_)
       do k=1,nK; do j=1,nJ
          call get_refined_b0(2*iFace-3,2*j-2,2*k-2)
          do iDim=1,nDim
             B0ResChange_DXSB(iDim,j,k,iSide,iBlock) = 0.125 * &
                  sum(RefB0_DIII(iDim,:,:,:))
          end do
       end do; end do
    case(South_,North_)
       jFace=1+nJ*(iSide-South_)
       do k=1,nK; do i=1,nI
          call get_refined_b0(2*i-2,2*jFace-3,2*k-2)
          do iDim=1,nDim
             B0ResChange_DYSB(iDim,i,k,iSide,iBlock) = 0.125 * &
                  sum(RefB0_DIII(iDim,:,:,:))
          end do
       end do; end do
    case(Bot_,Top_)
       kFace=1+nK*(iSide-Bot_)
       do j=1,nJ; do i=1,nI
          call get_refined_b0(2*i-2,2*j-2,2*kFace-3)
          do iDim=1,nDim
             B0ResChange_DZSB(iDim,i,j,iSide,iBlock) = 0.125 * &
                  sum(RefB0_DIII(iDim,:,:,:))
          end do
       end do; end do
    end select
  end subroutine correct_b0_face
  !===========================================================================
  subroutine get_refined_b0(iRef,jRef,kRef)
    integer,intent(in)::iRef,jRef,kRef
    integer::i2,j2,k2
    !-------------------------------------------------------------------------
    do k2=0,1; do j2=0,1; do i2=0,1
       x=RefXyzStart_D(x_)+(iRef+i2)*RefDXyz_D(x_)
       y=RefXyzStart_D(y_)+(jRef+j2)*RefDXyz_D(y_)
       z=RefXyzStart_D(z_)+(kRef+k2)*RefDXyz_D(z_)
       call get_b0(x,y,z,RefB0_DIII(:,i2,j2,k2))
    end do; end do; end do	
  end subroutine get_refined_b0
end subroutine set_b0_source

!============================================================================
subroutine calc_db0_dt(dTime)

  ! Calculate dB0/dt for true cells in used blocks

  use ModMain
  use ModAdvance, ONLY: B0_DGB, Db0Dt_CDB
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK
  use ModPhysics, ONLY: No2Si_V, UnitT_
  implicit none

  real, intent(in) :: dTime

  real :: B0_D(3), TimeSimulationOrig
  integer :: i,j,k,iBlock
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------

  call set_oktest('calc_db0_dt', DoTest, DoTestMe)
  if(DoTestMe)write(*,*)'calc_db0_dt starting, Dtime=',Dtime

  if(.not.allocated(Db0Dt_CDB))allocate(Db0Dt_CDB(nI,nJ,nK,3,MaxBlock))

  ! If time step is zero, so are the source terms and we may return
  if(Dtime <= 0.0)then
     Db0Dt_CDB = 0.0
     RETURN
  end if

  ! Increase simulation time and to what it will be at the end of the time step
  TimeSimulationOrig = Time_Simulation
  Time_Simulation = Time_Simulation + dTime*No2Si_V(UnitT_)

  ! Calculate dB0dt
  do iBlock=1,nBlock
     if(unusedBLK(iBlock))CYCLE
     do k=1,nK; do j=1,nJ; do i=1,nI
        call get_b0(&
             x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock),z_BLK(i,j,k,iBlock),B0_D)
        Db0Dt_CDB(i,j,k,x_,iBlock)=(B0_D(x_)-B0_DGB(x_,i,j,k,iBlock))/dTime
        Db0Dt_CDB(i,j,k,y_,iBlock)=(B0_D(y_)-B0_DGB(y_,i,j,k,iBlock))/dTime
        Db0Dt_CDB(i,j,k,z_,iBlock)=(B0_D(z_)-B0_DGB(z_,i,j,k,iBlock))/dTime
     end do; end do; end do
  end do
  ! Reset time, date and dipole tilt
  Time_Simulation = TimeSimulationOrig

end subroutine calc_db0_dt
!============================================================================
subroutine set_b0_matrix(iBlock)

  ! Calculate the elements of the B0 Source term

  use ModProcMH
  use ModMain
  use ModNumConst
  use ModB0
  use ModAdvance, ONLY  : &
       CurlB0_DCB,DivB0_CB, NormB0_CB                      
  use ModGeometry, ONLY : dx_BLK,dy_BLK,dz_BLK,true_cell
  use ModNumConst
  implicit none

  integer, intent(in) :: iBlock

  integer :: i,j,k,iDir,jDir

  ! inverse of Dx, Dy, Dz
  real :: DxInv, DyInv, DzInv

  logical :: oktest, oktest_me
  real::CurlB02,CurlB02_DD(3,3),B0Nabla_DD(3,3)
  real:: eigenvalue_max
  external eigenvalue_max
  !--------------------------------------------------------------------------

  if(iProc==PROCtest.and.iBlock==BLKtest)then
     call set_oktest('set_b0_matrix',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  endif

!  if(.not.allocated(CurlB0_DCB))then                !^CFG UNCOMMENT IF DYNAMIC
!       allocate(CurlB0_DCB(3,nI,nJ,nK,MaxBlock),&   !^CFG UNCOMMENT IF DYNAMIC
!                DivB0_CB(nI,nJ,nK,MaxBlock))        !^CFG UNCOMMENT IF DYNAMIC
!  if(UseCurlB0)allocate(NormB0_CB(nI,nJ,nK,MaxBlock)!^CFG UNCOMMENT IF DYNAMIC
!  end if                                            !^CFG UNCOMMENT IF DYNAMIC

  
  DxInv = cOne/dx_BLK(iBlock)
  DyInv = cOne/dy_BLK(iBlock)
  DzInv = cOne/dz_BLK(iBlock)
  call set_b0_face(iBlock)
  ! face areas
  do k=1,nK; do j=1,nJ; do i=1,nI
     DivB0_CB(i,j,k,iBlock)= &
          DxInv*(B0_DX(x_,i+1,j,k)-B0_DX(x_,i,j,k))+&
          DyInv*(B0_DY(y_,i,j+1,k)-B0_DY(y_,i,j,k))+&
          DzInv*(B0_DZ(z_,i,j,k+1)-B0_DZ(z_,i,j,k))   
  
     CurlB0_DCB(z_,i,j,k,iBlock) = &
          DxInv*(B0_DX(y_,i+1,j,k)-B0_DX(y_,i,j,k))-&
          DyInv*(B0_DY(x_,i,j+1,k)-B0_DY(x_,i,j,k))

     CurlB0_DCB(y_,i,j,k,iBlock) = &
          -DxInv*(B0_DX(z_,i+1,j,k)-B0_DX(z_,i,j,k))+&
           DzInv*(B0_DZ(x_,i,j,k+1)-B0_DZ(x_,i,j,k))
     CurlB0_DCB(x_,i,j,k,iBlock) = & 
          DyInv*(B0_DY(z_,i,j+1,k)-B0_DY(z_,i,j,k))-&
          DzInv*(B0_DZ(y_,i,j,k+1)-B0_DZ(y_,i,j,k))
     if(.not.UseCurlB0)CYCLE
     CurlB02=sum(CurlB0_DCB(:,i,j,k,iBlock)**2)
     if(CurlB02==cZero.or..not.true_cell(i,j,k,iBlock))then
        NormB0_CB(i,j,k,iBlock)=cZero
     else
        CurlB02_DD=CurlB02*cUnit_DD
        do jDir=1,3
           do iDir=1,3
              CurlB02_DD(iDir,jDir)=CurlB02_DD(iDir,jDir)-&
                   CurlB0_DCB(iDir,i,j,k,iBlock)*&
                   CurlB0_DCB(jDir,i,j,k,iBlock)
           end do
        end do
        B0Nabla_DD(x_,:)=DxInv*(B0_DX(:,i+1,j,k)-B0_DX(:,i,j,k))
        B0Nabla_DD(y_,:)=DyInv*(B0_DY(:,i,j+1,k)-B0_DY(:,i,j,k))
        B0Nabla_DD(z_,:)=DzInv*(B0_DZ(:,i,j,k+1)-B0_DZ(:,i,j,k))
        B0Nabla_DD=B0Nabla_DD-DivB0_CB(i,j,k,iBlock)*cUnit_DD
        NormB0_CB(i,j,k,iBlock)= sqrt(sqrt(eigenvalue_max(&
             matmul(transpose(B0Nabla_DD),matmul(&
                              CurlB02_DD,&
                              B0Nabla_DD))               )&
                                     )    )
     end if
  end do; end do; end do

end subroutine set_b0_matrix

!============================================================================
real function eigenvalue_max(A_DD)
  use ModNumConst
  !Returns the maximum eigenvalue of a SYMMETRIC matrix
  !Good only for express estimates of the spectral radius 
  implicit none
  real,dimension(3,3),intent(in)::A_DD
  real::Lambda, DetF, F_DD(3,3),FPrime_DDD(3,3,3)
  integer::i,iIteration
  real::Tolerance
  !-----------------------------------------------------------
  Lambda=cZero
  !Upper estimate for a spectral radius
  do i=1,3
     Lambda=max(Lambda,sum(abs(A_DD(i,:))))
  end do
  if(Lambda==cZero)then
     eigenvalue_max=cZero
     return
  end if
  Tolerance=(0.1*Lambda)**3

  F_DD=A_DD-Lambda*cUnit_DD
  DetF=det(F_DD) 
  iIteration=0
  do while (abs(DetF).gt.Tolerance)
     iIteration=iIteration+1
     do i=1,3
        FPrime_DDD(:,:,i)=F_DD
        FPrime_DDD(:,i,i)=cZero
        FPrime_DDD(i,:,i)=cZero
        FPrime_DDD(i,i,i)=cOne
     end do
     Lambda=Lambda+DetF/&
          (det(FPrime_DDD(:,:,1))+det(FPrime_DDD(:,:,2))+det(FPrime_DDD(:,:,3)))
     if(iIteration==10)exit
     F_DD=A_DD-Lambda*cUnit_DD
     DetF=det(F_DD)
  end do
  eigenvalue_max=Lambda
contains
  real function det(A_DD)
    implicit none
    real,dimension(3,3),intent(in)::A_DD
    det=A_DD(1,1)*(A_DD(2,2)*A_DD(3,3)-&
         A_DD(3,2)*A_DD(2,3))-&
         A_DD(1,2)*(A_DD(2,1)*A_DD(3,3)-&
         A_DD(2,3)*A_DD(3,1))+&
         A_DD(1,3)*(A_DD(2,1)*A_DD(3,2)-&
         A_DD(2,2)*A_DD(3,1))
  end function det
end function eigenvalue_max


!============================================================================

subroutine get_b0(X0,Y0,Z0,B0)

  ! this interface allows use of various models for B0
  use ModMain,          ONLY : Time_Simulation, NameThisComp, &
       TypeCoordSystem, IsStandAlone
  use ModPhysics,       ONLY : Si2No_V, UnitB_, DipoleStrengthSi
  use CON_planet_field, ONLY : get_planet_field
  use ModMain,          ONLY : UseBody2             !^CFG IF SECONDBODY
  use ModMain,          ONLY : UseUserB0
  use ModUser,          ONLY: user_get_b0
  implicit none

  real, intent(in) :: X0,Y0,Z0
  real, intent(out), dimension(3) :: B0

  if(UseUserB0)then
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
  real :: R0, rr, rr_inv, rr2_inv, rr3_inv, rr5_inv, rr7_inv, rr9_inv
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

  !\
  ! Add dipole field of the second body
  !/

  Dp = sum(BdpBody2_D*Xyz_D)*3*rr2_inv

  B0 = B0 + (Dp*Xyz_D - BdpBody2_D)*rr3_inv

end subroutine add_b0_body2
!^CFG END SECONDBODY

!==============================================================================

subroutine update_b0

  use ModMain,          ONLY: DoSplitDb0Dt, nBlock, unusedBLK, &
       time_simulation, NameThisComp
  use ModPhysics,       ONLY: ThetaTilt
  use ModAdvance,       ONLY: Bx_, By_, Bz_, State_VGB, B0_DGB
  use ModGeometry,      ONLY: true_cell, body_BLK
  use CON_axes,         ONLY: get_axes
  use ModNumConst,      ONLY: cRadToDeg
  use ModIO,            ONLY: iUnitOut, write_prefix
  use ModEnergy,        ONLY: calc_energy_ghost

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

     if(DoSplitDb0Dt)then
        ! Save total magnetic field into Bx_BLK,By_BLK,Bz_BLK
        State_VGB(Bx_:Bz_,:,:,:,iBlock) = State_VGB(Bx_:Bz_,:,:,:,iBlock) &
             + B0_DGB(:,:,:,:,iBlock)
     end if

     call set_b0(iBlock)

     if(DoSplitDb0Dt)then
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

     end if
  end do
  call timing_stop(NameSub)
end subroutine update_b0

!===========================================================================
