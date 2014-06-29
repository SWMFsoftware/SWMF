!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module ModFaceGradient

  use ModSize, ONLY: MaxDim, nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       j0_, j2_, nJp1_, nJm1_, k0_, k2_, nKp1_, nKm1_, &
       jRatio, kRatio, InvIjkRatio

  implicit none
  save

  private ! except

  ! Public methods
  public :: set_block_field2
  public :: set_block_field3
  public :: get_face_gradient
  public :: get_face_gradient_field
  public :: get_face_curl
  public :: set_block_jacobian_face

  ! Jacobian matrix for general grid: Dgencoord/Dcartesian
  real, public :: DcoordDxyz_DDFD(MaxDim,MaxDim,1:nI+1,1:nJ+1,1:nK+1,MaxDim)

contains

  !============================================================================

  subroutine set_block_field2(iBlock, nVar, Field1_VG, Field_VG)

    ! correct the ghostcells of the given scalar/vector field on iBlock
    ! using second order interpolation

    use ModParallel, ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, NOBLK
    use BATL_lib, ONLY: DiLevelNei_IIIB

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: Field1_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(inout) :: Field_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    real, parameter :: c0 = 0.5, p0 = 1./6., F1 = 1./3.

    integer :: i1, j1, k1, i2, j2, k2, iC, jC, kC, iSide, jSide, kSide
    integer :: iL, iR, jL, jR, kL, kR
    integer :: ip, jp, kp

    logical :: IsEqualLevel_G(0:nI+1,0:nJ+1,0:nK+1)
    !--------------------------------------------------------------------------

    Field1_VG = Field_VG

    do kSide = -1,1; 
       if(nK == 1 .and. kSide /= 0) CYCLE
       do jSide = -1,1
          if(nJ == 1 .and. jSide /= 0) CYCLE
          do iSide = -1,1
             if(iSide==0)then
                iL = 1; iR = nI
             elseif(iSide==1)then
                iL = nI+1; iR = iL
             else
                iL = 0; iR = iL
             end if
             if(jSide==0)then
                jL = 1; jR = nJ
             elseif(jSide==1)then
                jL = nJ+1; jR = jL
             else
                jL = 0; jR = jL
             end if
             if(kSide==0)then
                kL = 1; kR = nK
             elseif(kSide==1)then
                kL = nK+1; kR = kL
             else
                kL = 0; kR = kL
             end if
             if( DiLevelNei_IIIB(iSide, jSide, kSide, iBlock) == 0 )then
                IsEqualLevel_G(iL:iR,jL:jR,kL:kR) = .true.
             else
                IsEqualLevel_G(iL:iR,jL:jR,kL:kR) = .false.
             end if
          end do; 
       end do; 
    end do

    ! Do the faces
    if(NeiLeast(iBlock) == 1)then
       if(nJ == 1)then
          ! Interpolate in 1D. Fine cell center is twice further than coarse
          Field_VG(:,0,1,1) = (2*Field1_VG(:,0,1,1) + Field1_VG(:,1,1,1))/3
       else
          ! Interpolate 3 co-planar cell centers 
          ! 0,j2,k2 is the fine ghost cell to be filled in, 
          !         originally it is the coarse cell (1st order prolongation)
          ! 1,j2,k2 is the fine physical cell next to the ghost cell
          ! 0,jp,kp is a fine/coarse nearby ghost cell chosen in the same plane
          !
          ! When 0,jp,kp is a fine cell they are interpolated to the corner
          ! with weights 1/2, when 0,jp,kp is coarse, it is interpolated with
          ! weights 1/3 and 2/3.
          !
          ! The corner and coarse cell center are then interpolated

          do k1=1, nK, 2; do j1=1, nJ, 2; 
             do k2 = k1,k1+min(1,nK-1); do j2 = j1,j1+1
                jp = 3*j2 - 2*j1 -1
                if(nK == 1) kp = 1
                if(nK >  1) kp = 3*k2 - 2*k1 -1
                if(IsEqualLevel_G(0,jp,kp))then
                   Field_VG(:,0,j2,k2) = c0*Field1_VG(:,0,j2,k2) &
                        + 0.25*(Field1_VG(:,0,jp,kp) + Field_VG(:,1,j2,k2))
                else
                   Field_VG(:,0,j2,k2) = c0*Field1_VG(:,0,j2,k2) &
                        + p0*Field1_VG(:,0,jp,kp) + F1*Field_VG(:,1,j2,k2)
                end if
             end do; end do; 
          end do; end do
       end if
    end if

    if(NeiLwest(iBlock) == 1)then
       if(nJ == 1)then
          Field_VG(:,nI+1,1,1) = &
               (2*Field1_VG(:,nI+1,1,1) + Field1_VG(:,nI,1,1))/3
       else
          do k1=1, nK, 2; do j1=1, nJ, 2; 
             do k2 = k1,k1+min(1,nK-1); do j2 = j1,j1+1
                jp = 3*j2 - 2*j1 -1
                if(nK == 1) kp = 1
                if(nK >  1) kp = 3*k2 - 2*k1 -1
                if(IsEqualLevel_G(nI+1,jp,kp))then
                   Field_VG(:,nI+1,j2,k2) = c0*Field1_VG(:,nI+1,j2,k2) &
                        + 0.25*(Field1_VG(:,nI+1,jp,kp) + Field_VG(:,nI,j2,k2))
                else
                   Field_VG(:,nI+1,j2,k2) = c0*Field1_VG(:,nI+1,j2,k2) &
                        + p0*Field1_VG(:,nI+1,jp,kp) + F1*Field_VG(:,nI,j2,k2)
                end if
             end do; end do; 
          end do; end do
       end if
    end if

    if(nJ == 1) RETURN

    if(NeiLsouth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2; 
          do k2 = k1,k1+min(1,nK-1); do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1
             if(nK == 1) kp = 1
             if(nK >  1) kp = 3*k2 - 2*k1 -1
             if(IsEqualLevel_G(ip,j0_,kp))then
                Field_VG(:,i2,j0_,k2) = c0*Field1_VG(:,i2,j0_,k2) &
                     + 0.25*(Field1_VG(:,ip,j0_,kp) + Field_VG(:,i2,1,k2))
             else
                Field_VG(:,i2,j0_,k2) = c0*Field1_VG(:,i2,j0_,k2) &
                     + p0*Field1_VG(:,ip,j0_,kp) + F1*Field_VG(:,i2,1,k2)
             end if
          end do; end do; 
       end do; end do
    end if

    if(NeiLnorth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2; 
          do k2 = k1,k1+min(1,nK-1); do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1
             if(nK == 1) kp = 1
             if(nK >  1) kp = 3*k2 - 2*k1 -1
             if(IsEqualLevel_G(ip,nJp1_,kp))then
                Field_VG(:,i2,nJp1_,k2) = c0*Field1_VG(:,i2,nJp1_,k2) &
                     + 0.25*(Field1_VG(:,ip,nJp1_,kp) + Field_VG(:,i2,nJ,k2))
             else
                Field_VG(:,i2,nJp1_,k2) = c0*Field1_VG(:,i2,nJp1_,k2) &
                     + p0*Field1_VG(:,ip,nJp1_,kp) + F1*Field_VG(:,i2,nJ,k2)
             end if
          end do; end do;
       end do; end do
    end if

    if(nK > 1 .and. NeiLbot(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          jp = 3*j2 - 2*j1 -1
          if(IsEqualLevel_G(ip,jp,k0_))then
             Field_VG(:,i2,j2,k0_) = c0*Field1_VG(:,i2,j2,k0_) &
                  + 0.25*Field1_VG(:,ip,jp,k0_) + 0.25*Field_VG(:,i2,j2,1)
          else
             Field_VG(:,i2,j2,k0_) = c0*Field1_VG(:,i2,j2,k0_) &
                  + p0*Field1_VG(:,ip,jp,k0_) + F1*Field_VG(:,i2,j2,1)
          end if
       end do; end do; end do; end do
    end if

    if(nK > 1 .and. NeiLtop(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1
          jp = 3*j2 - 2*j1 -1
          if(IsEqualLevel_G(ip,jp,nKp1_))then
             Field_VG(:,i2,j2,nKp1_) = c0*Field1_VG(:,i2,j2,nKp1_) &
                  + 0.25*Field1_VG(:,ip,jp,nKp1_) + 0.25*Field_VG(:,i2,j2,nK)
          else
             Field_VG(:,i2,j2,nKp1_) = c0*Field1_VG(:,i2,j2,nKp1_) &
                  + p0*Field1_VG(:,ip,jp,nKp1_) + F1*Field_VG(:,i2,j2,nK)
          end if
       end do; end do; end do; end do
    end if

    ! Do the edges

    ! 4 Z edges
    do jSide = -1,1,2; do iSide = -1,1,2
       if(  DiLevelNei_IIIB(iSide, jSide, 0, iBlock) /= 1.and. .not. ( &
            DiLevelNei_IIIB(iSide, jSide, 0, iBlock) == NOBLK .and. ( &
            DiLevelNei_IIIB(iSide, 0, 0, iBlock) == 1 .or. &
            DiLevelNei_IIIB(0, jSide, 0, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; iC = i1+iSide
       j1=1; if(jSide==1) j1=nJ; jC = j1+jSide
       do k1 = 1, nK, 2 ; do k2 = k1, k1 + min(nK-1,1)
          if(nK == 1) then
             kp = 1
          else
             kp = 3*k2 - 2*k1 -1
          end if
          if(IsEqualLevel_G(iC,jC,kp))then
             Field_VG(:,iC,jC,k2) = c0*Field1_VG(:,iC,jC,k2) &
                  + 0.25*Field1_VG(:,iC,jC,kp) + 0.25*Field_VG(:,i1,j1,k2)
          else
             Field_VG(:,iC,jC,k2) = c0*Field1_VG(:,iC,jC,k2) &
                  + p0*Field1_VG(:,iC,jC,kp) + F1*Field_VG(:,i1,j1,k2)
          end if
       end do; end do         
    end do; end do

    ! The X and Y edges are not needed in 2D
    if(nK == 1) RETURN

    ! 4 X edges
    do kSide = -1,1,2; do jSide = -1,1,2
       if(  DiLevelNei_IIIB(0, jSide, kSide, iBlock) /= 1 .and. .not. ( &
            DiLevelNei_IIIB(0, jSide, kSide, iBlock) == NOBLK .and. ( &
            DiLevelNei_IIIB(0, jSide, 0, iBlock) == 1 .or. &
            DiLevelNei_IIIB(0, 0, kSide, iBlock) == 1))) CYCLE

       j1=1; if(jSide==1) j1=nJ; jC = j1+jSide
       k1=1; if(kSide==1) k1=nK; kC = k1+kSide
       do i1 = 1,nI,2; do i2 = i1, i1+1
          ip = 3*i2 - 2*i1 -1
          if(IsEqualLevel_G(ip,jC,kC))then
             Field_VG(:,i2,jC,kC) = c0*Field1_VG(:,i2,jC,kC) &
                  + 0.25*Field1_VG(:,ip,jC,kC) + 0.25*Field_VG(:,i2,j1,k1)
          else
             Field_VG(:,i2,jC,kC) = c0*Field1_VG(:,i2,jC,kC) &
                  + p0*Field1_VG(:,ip,jC,kC) + F1*Field_VG(:,i2,j1,k1)
          end if
       end do; end do
    end do; end do
    ! 4 Y edges
    do kSide = -1,1,2; do iSide = -1,1,2
       if(  DiLevelNei_IIIB(iSide, 0, kSide, iBlock) /= 1 .and. .not. ( &
            DiLevelNei_IIIB(iSide, 0, kSide, iBlock) == NOBLK .and. ( &
            DiLevelNei_IIIB(iSide, 0, 0, iBlock) == 1 .or. &
            DiLevelNei_IIIB(0, 0, kSide, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; iC = i1+iSide
       k1=1; if(kSide==1) k1=nK; kC = k1+kSide
       do j1 = 1, nJ, 2; do j2 = j1, j1+1
          jp = 3*j2 - 2*j1 -1
          if(IsEqualLevel_G(iC,jp,kC))then
             Field_VG(:,iC,j2,kC) = c0*Field1_VG(:,iC,j2,kC) &
                  + 0.25*Field1_VG(:,iC,jp,kC) + 0.25*Field_VG(:,i1,j2,k1)
          else
             Field_VG(:,iC,j2,kC) = c0*Field1_VG(:,iC,j2,kC) &
                  + p0*Field1_VG(:,iC,jp,kC) + F1*Field_VG(:,i1,j2,k1)
          end if
       end do; end do
    end do; end do


  end subroutine set_block_field2

  !============================================================================

  subroutine set_block_field3(iBlock, nVar, Field1_VG, Field_VG)

    ! correct the ghost cells of the given scalar/vector field on iBlock
    ! using third order interpolation

    use ModParallel, ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, NOBLK
    use BATL_lib, ONLY: DiLevelNei_IIIB

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: Field1_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(inout) :: Field_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    real,parameter :: C1 = 8./15., F1 = 2./3., F2 = -1./5.
    real,parameter :: p0 = 5./32., m0 =-3./32., c0= 15./16.

    integer :: i1, j1, k1, i2, j2, k2, iC, jC, kC, iSide, jSide, kSide
    integer :: iL, iR, jL, jR, kL, kR
    integer :: ip, im, jp, jm, kp, km, iVar

    real :: Fc_V(nVar) ! interpolated coarse cell value
    !-------------------------------------------------------------------------

    Field1_VG = Field_VG

    ! Fix ghost edge and corner ghost cells
    do kSide = -1,1; do jSide = -1,1; do iSide = -1,1
       ! If the corner or edge is not coarser but the face is
       ! then average out the 8 fine cells so that the
       ! general interpolation formulas remain 2nd order
       ! accurate
       if(  DiLevelNei_IIIB(iSide,jSide,kSide,iBlock) /= 1 .and. &
            DiLevelNei_IIIB(iSide,jSide,kSide,iBlock) /= NOBLK .and. ( &
            DiLevelNei_IIIB(iSide,0,0,iBlock) == 1 .or. &
            DiLevelNei_IIIB(0,jSide,0,iBlock) == 1 .or. &
            DiLevelNei_IIIB(0,0,kSide,iBlock) == 1)) then

          if(iSide==0)then
             iL = 1; iR = nI
          elseif(iSide==1)then
             iL = nI+1; iR=nI+2
          else
             iL = -1; iR = 0
          end if
          if(jSide==0)then
             jL = 1; jR = nJ
          elseif(jSide==1)then
             jL = nJ+1; jR=nJ+2
          else
             jL = -1; jR = 0
          end if
          if(kSide==0)then
             kL = 1; kR = nK
          elseif(kSide==1)then
             kL = nK+1; kR=nK+2
          else
             kL = -1; kR = 0
          end if

          do k1=kL,kR,2; do j1=jL,jR,2; do i1=iL,iR,2; do iVar=1,nVar
             Field1_VG(iVar,i1:i1+1,j1:j1+jRatio-1,k1:k1+kRatio-1)= &
                  sum(Field_VG(iVar,i1:i1+1,j1:j1+jRatio-1,k1:k1+kRatio-1)) &
                  *InvIjkRatio
          end do; end do; end do; end do

       end if
    end do; end do; end do

    ! Do six faces
    if(NeiLeast(iBlock) == 1)then
       if(nJ == 1)then
          Field_VG(:,0,1,1) = max(min( C1*Field1_VG(:,0,1,1) &
               + F1*Field_VG(:,1,1,1) + F2*Field_VG(:,2,1,1), &
               max(Field1_VG(:,0,1,1), Field_VG(:,1,1,1))), &
               min(Field1_VG(:,0,1,1), Field_VG(:,1,1,1)))
       else
          do k1=1, nK, 2; do j1=1, nJ, 2;
             do k2 = k1,k1+min(1,nK-1); do j2 = j1,j1+1
                jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
                if(nK == 1)then
                   kp = 1; km = 1
                else
                   kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
                end if
                Fc_V = c0*Field1_VG(:,0,j2,k2) &
                     + p0*Field1_VG(:,0,jp,kp) &
                     + m0*Field1_VG(:,0,jm,km)

                ! Limit the interpolation formula with the max and min
                ! of the 3 points surrounding the interpolated point
                Field_VG(:,0,j2,k2) = max(min( C1*Fc_V &
                     + F1*Field_VG(:,1,j2,k2) + F2*Field_VG(:,2,j2,k2),&
                     max(Field1_VG(:,0,j2,k2), &
                     Field1_VG(:,0,jp,kp), Field_VG(:,1,j2,k2))), &
                     min(Field1_VG(:,0,j2,k2), &
                     Field1_VG(:,0,jp,kp), Field_VG(:,1,j2,k2)))
             end do; end do
          end do; end do
       end if
    end if

    if(NeiLwest(iBlock) == 1)then
       if(nJ == 1)then
          Field_VG(:,nI+1,1,1) = max(min( C1*Field1_VG(:,nI+1,1,1) &
               + F1*Field_VG(:,nI,1,1) + F2*Field_VG(:,nI-1,1,1), &
               max(Field1_VG(:,nI+1,1,1), Field_VG(:,nI,1,1))), &
               min(Field1_VG(:,nI+1,1,1), Field_VG(:,nI,1,1)))
       else
          do k1=1, nK, 2; do j1=1, nJ, 2
             do k2 = k1,k1+min(1,nK-1); do j2 = j1,j1+1
                jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
                if(nK == 1)then
                   kp = 1; km = 1
                else
                   kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
                end if
                Fc_V = c0*Field1_VG(:,nI+1,j2,k2) &
                     + p0*Field1_VG(:,nI+1,jp,kp) &
                     + m0*Field1_VG(:,nI+1,jm,km)
                Field_VG(:,nI+1,j2,k2) = max(min( C1*Fc_V &
                     + F1*Field_VG(:,nI,j2,k2) + F2*Field_VG(:,nI-1,j2,k2), &
                     max(Field1_VG(:,nI+1,j2,k2), &
                     Field1_VG(:,nI+1,jp,kp), Field_VG(:,nI,j2,k2))), &
                     min(Field1_VG(:,nI+1,j2,k2), &
                     Field1_VG(:,nI+1,jp,kp), Field_VG(:,nI,j2,k2)))
             end do; end do
          end do; end do
       end if
    end if

    if(nJ == 1) RETURN

    if(NeiLsouth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2
          do k2 = k1,k1+min(1,nK-1); do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             if(nK == 1)then
                kp = 1; km = 1
             else
                kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
             end if
             Fc_V = c0*Field1_VG(:,i2,j0_,k2) &
                  + p0*Field1_VG(:,ip,j0_,kp) &
                  + m0*Field1_VG(:,im,j0_,km)
             Field_VG(:,i2,j0_,k2) = max(min( &
                  C1*Fc_V + F1*Field_VG(:,i2,1,k2) + F2*Field_VG(:,i2,j2_,k2), &
                  max(Field1_VG(:,i2,j0_,k2), &
                  Field1_VG(:,ip,j0_,kp), Field_VG(:,i2,1,k2))), &
                  min(Field1_VG(:,i2,j0_,k2), &
                  Field1_VG(:,ip,j0_,kp), Field_VG(:,i2,1,k2)))
          end do; end do
       end do; end do
    end if

    if(NeiLnorth(iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2
          do k2 = k1,k1+min(1,nK-1); do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             if(nK == 1)then
                kp = 1; km = 1
             else
                kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
             end if
             Fc_V = c0*Field1_VG(:,i2,nJp1_,k2) &
                  + p0*Field1_VG(:,ip,nJp1_,kp) &
                  + m0*Field1_VG(:,im,nJp1_,km)
             Field_VG(:,i2,nJp1_,k2) = max(min( C1*Fc_V &
                  + F1*Field_VG(:,i2,nJ,k2) + F2*Field_VG(:,i2,nJm1_,k2), &
                  max(Field1_VG(:,i2,nJp1_,k2), &
                  Field1_VG(:,ip,nJp1_,kp), Field_VG(:,i2,nJ,k2))), &
                  min(Field1_VG(:,i2,nJp1_,k2), &
                  Field1_VG(:,ip,nJp1_,kp), Field_VG(:,i2,nJ,k2)))
          end do; end do
       end do; end do
    end if

    if(nK > 1 .and. NeiLbot(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
          jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
          Fc_V = c0*Field1_VG(:,i2,j2,k0_) &
               + p0*Field1_VG(:,ip,jp,k0_) &
               + m0*Field1_VG(:,im,jm,k0_)
          Field_VG(:,i2,j2,k0_) = max(min( &
               C1*Fc_V + F1*Field_VG(:,i2,j2,1) + F2*Field_VG(:,i2,j2,k2_), &
               max(Field1_VG(:,i2,j2,k0_), &
               Field1_VG(:,ip,jp,k0_), Field_VG(:,i2,j2,1))), &
               min(Field1_VG(:,i2,j2,k0_), &
               Field1_VG(:,ip,jp,k0_), Field_VG(:,i2,j2,1)))
       end do; end do; end do; end do
    end if

    if(nK > 1 .and. NeiLtop(iBlock) == 1)then
       do j1=1, nJ, 2; do i1=1, nI, 2; do j2 = j1,j1+1; do i2 = i1,i1+1
          ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
          jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
          Fc_V = c0*Field1_VG(:,i2,j2,nKp1_) &
               + p0*Field1_VG(:,ip,jp,nKp1_) &
               + m0*Field1_VG(:,im,jm,nKp1_)
          Field_VG(:,i2,j2,nKp1_) = max(min( &
               C1*Fc_V + F1*Field_VG(:,i2,j2,nK) + F2*Field_VG(:,i2,j2,nKm1_), &
               max(Field1_VG(:,i2,j2,nKp1_), &
               Field1_VG(:,ip,jp,nKp1_), Field_VG(:,i2,j2,nK))), &
               min(Field1_VG(:,i2,j2,nKp1_), &
               Field1_VG(:,ip,jp,nKp1_), Field_VG(:,i2,j2,nK)))
       end do; end do; end do; end do
    end if

    ! Do the edges

    ! 4 Z edges
    do jSide = -1,1,2; do iSide = -1,1,2
       if(  DiLevelNei_IIIB(iSide, jSide, 0, iBlock) /= 1.and. .not. ( &
            DiLevelNei_IIIB(iSide, jSide, 0, iBlock) == NOBLK .and. ( &
            DiLevelNei_IIIB(iSide, 0, 0, iBlock) == 1 .or. &
            DiLevelNei_IIIB(0, jSide, 0, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; i2 = i1-iSide; iC = i1+iSide
       j1=1; if(jSide==1) j1=nJ; j2 = j1-jSide; jC = j1+jSide
       do k1 = 1, nK, 2 ; do k2 = k1, k1 + min(1,nK-1)
          if(nK == 1)then
             kp = 1; km = 1
          else
             kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
          end if
          Fc_V = c0*Field1_VG(:,iC,jC,k2) &
               + p0*Field1_VG(:,iC,jC,kp) &
               + m0*Field1_VG(:,iC,jC,km)
          Field_VG(:,iC,jC,k2) = max(min( &
               C1*Fc_V + F1*Field_VG(:,i1,j1,k2) + F2*Field_VG(:,i2,j2,k2), &
               max(Field1_VG(:,iC,jC,k2), &
               Field1_VG(:,iC,jC,kp), Field_VG(:,i1,j1,k2))), &
               min(Field1_VG(:,iC,jC,k2), &
               Field1_VG(:,iC,jC,kp), Field_VG(:,i1,j1,k2)))
       end do;end do
    end do;end do

    ! The X and Y edges are not needed in 2D
    if(nK == 1) RETURN

    ! 4 X edges
    do kSide = -1,1,2; do jSide = -1,1,2
       if(  DiLevelNei_IIIB(0, jSide, kSide, iBlock) /= 1 .and. .not. ( &
            DiLevelNei_IIIB(0, jSide, kSide, iBlock) == NOBLK .and. ( &
            DiLevelNei_IIIB(0, jSide, 0, iBlock) == 1 .or. &
            DiLevelNei_IIIB(0, 0, kSide, iBlock) == 1))) CYCLE

       j1=1; if(jSide==1) j1=nJ; j2 = j1-jSide; jC = j1+jSide
       k1=1; if(kSide==1) k1=nK; k2 = k1-kSide; kC = k1+kSide
       do i1 = 1,nI,2; do i2 = i1, i1+1
          ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
          Fc_V = c0*Field1_VG(:,i2,jC,kC) &
               + p0*Field1_VG(:,ip,jC,kC) &
               + m0*Field1_VG(:,im,jC,kC)
          Field_VG(:,i2,jC,kC) = max(min( &
               C1* Fc_V + F1*Field_VG(:,i2,j1,k1) + F2*Field_VG(:,i2,j2,k2), &
               max(Field1_VG(:,i2,jC,kC), &
               Field1_VG(:,ip,jC,kC), Field_VG(:,i2,j1,k1))), &
               min(Field1_VG(:,i2,jC,kC), &
               Field1_VG(:,ip,jC,kC), Field_VG(:,i2,j1,k1)))
       end do;end do
    end do;end do
    ! 4 Y edges
    do kSide = -1,1,2; do iSide = -1,1,2
       if(  DiLevelNei_IIIB(iSide, 0, kSide, iBlock) /= 1 .and. .not. ( &
            DiLevelNei_IIIB(iSide, 0, kSide, iBlock) == NOBLK .and. ( &
            DiLevelNei_IIIB(iSide, 0, 0, iBlock) == 1 .or. &
            DiLevelNei_IIIB(0, 0, kSide, iBlock) == 1))) CYCLE

       i1=1; if(iSide==1) i1=nI; i2 = i1-iSide; iC = i1+iSide
       k1=1; if(kSide==1) k1=nK; k2 = k1-kSide; kC = k1+kSide
       do j1 = 1, nJ, 2; do j2 = j1, j1+1
          jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
          Fc_V = c0*Field1_VG(:,iC,j2,kC) &
               + p0*Field1_VG(:,iC,jp,kC) &
               + m0*Field1_VG(:,iC,jm,kC)
          Field_VG(:,iC,j2,kC) = max(min( &
               C1*Fc_V + F1*Field_VG(:,i1,j2,k1) + F2*Field_VG(:,i2,j2,k2), &
               max(Field1_VG(:,iC,j2,kC), &
               Field1_VG(:,iC,jp,kC), Field_VG(:,i1,j2,k1))), &
               min(Field1_VG(:,iC,j2,kC), &
               Field1_VG(:,iC,jp,kC), Field_VG(:,i1,j2,k1)))
       end do;end do
    end do;end do

  end subroutine set_block_field3

  !============================================================================

  subroutine set_block_jacobian_face(iBlock, UseFirstOrderBcIn)

    use ModMain, ONLY: x_, y_, z_, nDim
    use ModNumConst, ONLY: i_DD
    use ModCoordTransform, ONLY: inverse_matrix
    use BATL_lib, ONLY: Xyz_DGB, CellSize_DB
    use ModParallel,   ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, NoBlk
    integer, intent(in):: iBlock
    logical, optional, intent(in):: UseFirstOrderBcIn

    ! Dxyz/Dcoord matrix for one cell
    real:: DxyzDcoord_DD(MaxDim,MaxDim)

    ! Transverse gradients
    real:: TransGrad_DDG(MaxDim,MaxDim,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    ! Inverse of cell size
    real :: InvDx, InvDy, InvDz

    ! Indexes
    integer:: i, j, k

    !coeff of Ui+2 and Ui+1 to get normal derivative
    real, parameter:: fp2 = -1./24.0, fp1 = 9.0/8.0 
    !coeff of Ui+2 and Ui+1 for transverse derivatives
    real, parameter:: dp2 = -1./12.0, dp1 = 2.0/3.0 
    !coeff to average transverse derivatives
    real, parameter:: ap2 = -1./16.0, ap1 = 9.0/16.
    real :: Dxyz_D(MaxDim) 

    logical :: DoTest, DoTestMe, UseFirstOrderBc
    character(len=*), parameter:: NameSub='set_block_jacobian_face'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Calculate the dGencoord/dCartesian matrix

    InvDx = 1.0/CellSize_DB(x_,iBlock)
    InvDy = 1.0/CellSize_DB(y_,iBlock)
    InvDz = 1.0/CellSize_DB(z_,iBlock)

    if(present(UseFIrstOrderBcIn))then
       UseFirstOrderBc = UseFirstOrderBcIn
    else
       UseFirstOrderBc =  .false.
    end if

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=1,nI
       TransGrad_DDG(:,1,i,j,k)=  &
            ( dp1* (Xyz_DGB(:,i+1,j,k,iBlock) - Xyz_DGB(:,i-1,j,k,iBlock)) &
            + dp2* (Xyz_DGB(:,i+2,j,k,iBlock) - Xyz_DGB(:,i-2,j,k,iBlock)))
    end do; end do; end do

    if(nJ > 1)then
       do k=MinK,MaxK; do j=1,nJ; do i=MinI,MaxI
          TransGrad_DDG(:,2,i,j,k)=  &
               ( dp1* (Xyz_DGB(:,i,j+1,k,iBlock) - Xyz_DGB(:,i,j-1,k,iBlock)) &
               + dp2* (Xyz_DGB(:,i,j+2,k,iBlock) - Xyz_DGB(:,i,j-2,k,iBlock)))
       end do; end do; end do
    end if

    if(nK > 1)then
       do k=1,nK; do j=MinJ,MaxJ; do i=MinI,MaxI
          TransGrad_DDG(:,3,i,j,k)=  &
               ( dp1* (Xyz_DGB(:,i,j,k+1,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock)) &
               + dp2* (Xyz_DGB(:,i,j,k+2,iBlock) - Xyz_DGB(:,i,j,k-2,iBlock)))
       end do; end do; end do
    end if

    ! coord1 face
    do k=1,nK; do j=1,nJ; do i=1,nI+1
       ! DxyzDcoord along coord1 face
       DxyzDcoord_DD(:,1) = InvDx* &
            (  fp1*(Xyz_DGB(:,i  ,j,k,iBlock) - Xyz_DGB(:,i-1,j,k,iBlock)) &
            +  fp2*(Xyz_DGB(:,i+1,j,k,iBlock) - Xyz_DGB(:,i-2,j,k,iBlock)))
       if(nJ > 1)then
          DxyzDcoord_DD(:,2) = InvDy* &
               ( ap1*( TransGrad_DDG(:,2,i  ,j,k)+TransGrad_DDG(:,2,i-1,j,k)) &
               + ap2*( TransGrad_DDG(:,2,i+1,j,k)+TransGrad_DDG(:,2,i-2,j,k)))
       else
          DxyzDcoord_DD(:,2) = i_DD(:,2)
       end if
       if(nK > 1)then
          DxyzDcoord_DD(:,3) = InvDz* &
               ( ap1*( TransGrad_DDG(:,3,i  ,j,k)+TransGrad_DDG(:,3,i-1,j,k)) &
               + ap2*( TransGrad_DDG(:,3,i+1,j,k)+TransGrad_DDG(:,3,i-2,j,k)))
       else
          DxyzDcoord_DD(:,3) = i_DD(:,3)
       end if
       DcoordDxyz_DDFD(:,:,i,j,k,1) = &
            inverse_matrix(DxyzDcoord_DD, DoIgnoreSingular=.true.)
    end do; end do; end do

    ! coord2 face
    if(nJ > 1)then
       do k=1,nK; do j=1,nJ+1; do i=1,nI
          ! DxyzDcoord along coord2 face
          DxyzDcoord_DD(:,1) = InvDx* &
               ( ap1*( TransGrad_DDG(:,1,i,j  ,k)+TransGrad_DDG(:,1,i,j-1,k)) &
               + ap2*( TransGrad_DDG(:,1,i,j+1,k)+TransGrad_DDG(:,1,i,j-2,k)))
          DxyzDcoord_DD(:,2) = InvDy* &
               (  fp1*(Xyz_DGB(:,i,j  ,k,iBlock) - Xyz_DGB(:,i,j-1,k,iBlock)) &
               +  fp2*(Xyz_DGB(:,i,j+1,k,iBlock) - Xyz_DGB(:,i,j-2,k,iBlock)))
          if(nK > 1)then
             DxyzDcoord_DD(:,3) = InvDz* &
                  ( ap1*( TransGrad_DDG(:,3,i,j  ,k)  &
                  +       TransGrad_DDG(:,3,i,j-1,k)) &
                  + ap2*( TransGrad_DDG(:,3,i,j+1,k)  &
                  +       TransGrad_DDG(:,3,i,j-2,k)))
          else
             DxyzDcoord_DD(:,3) = i_DD(:,3)
          end if
          DcoordDxyz_DDFD(:,:,i,j,k,2) = &
               inverse_matrix(DxyzDcoord_DD, DoIgnoreSingular=.true.)
       end do; end do; end do
    end if

    ! coord3 face
    if(nK > 1)then
       do k=1,nK+1; do j=1,nJ; do i=1,nI
          ! DxyzDcoord along coord3 face
          DxyzDcoord_DD(:,1) = InvDx* &
               ( ap1*( TransGrad_DDG(:,1,i,j,k  )+TransGrad_DDG(:,1,i,j,k-1)) &
               + ap2*( TransGrad_DDG(:,1,i,j,k+1)+TransGrad_DDG(:,1,i,j,k-2)))
          DxyzDcoord_DD(:,2) = InvDy* &
               ( ap1*( TransGrad_DDG(:,2,i,j,k  )+TransGrad_DDG(:,2,i,j,k-1)) &
               + ap2*( TransGrad_DDG(:,2,i,j,k+1)+TransGrad_DDG(:,2,i,j,k-2)))
          DxyzDcoord_DD(:,3) = InvDz* &
               (  fp1*(Xyz_DGB(:,i,j,k  ,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock)) &
               +  fp2*(Xyz_DGB(:,i,j,k+1,iBlock) - Xyz_DGB(:,i,j,k-2,iBlock)))
          DcoordDxyz_DDFD(:,:,i,j,k,3) = &
               inverse_matrix(DxyzDcoord_DD, DoIgnoreSingular=.true.)
       end do; end do; end do
    end if
    if(.not.UseFirstOrderBc)RETURN
    if(neiLEast(iBlock)==NoBlk)then
       do k=1,nK; do j=1,nJ
          Dxyz_D = Xyz_DGB(:,1,j,k,iBlock) - Xyz_DGB(:,0,j,k,iBlock)
          Dxyz_D = Dxyz_D*(CellSize_DB(x_,iBlock)/sum(Dxyz_D**2))
          DcoordDxyz_DDFD(x_,:,1,j,k,x_) = Dxyz_D
       end do; end do
    end if
    if(neiLWest(iBlock)==NoBlk)then
       do k=1,nK; do j=1,nJ
          Dxyz_D = Xyz_DGB(:,nI + 1,j,k,iBlock) - Xyz_DGB(:,nI,j,k,iBlock)
          Dxyz_D = Dxyz_D*(CellSize_DB(x_,iBlock)/sum(Dxyz_D**2))
          DcoordDxyz_DDFD(x_,:,nI + 1,j,k,x_) = Dxyz_D
       end do; end do
    end if
    if(nJ==1)RETURN
    if(neiLSouth(iBlock)==NoBlk)then
       do k=1,nK; do i=1,nI
          Dxyz_D = Xyz_DGB(:,i,1,k,iBlock) - Xyz_DGB(:,i,0,k,iBlock)
          Dxyz_D = Dxyz_D*(CellSize_DB(y_,iBlock)/sum(Dxyz_D**2))
          DcoordDxyz_DDFD(y_,:,i,1,k,y_) = Dxyz_D
       end do; end do
    end if
    if(neiLNorth(iBlock)==NoBlk)then
       do k=1,nK; do i=1,nI
          Dxyz_D = Xyz_DGB(:,i,nJ + 1,k,iBlock) - Xyz_DGB(:,i,nJ,k,iBlock)
          Dxyz_D = Dxyz_D*(CellSize_DB(y_,iBlock)/sum(Dxyz_D**2))
          DcoordDxyz_DDFD(y_,:,i,nJ+1,k,y_) = Dxyz_D
       end do; end do
    end if
   if(neiLBot(iBlock)==NoBlk)then
       do j=1,nJ; do i=1,nI
          Dxyz_D = Xyz_DGB(:,i,j,1,iBlock) - Xyz_DGB(:,i,j,0,iBlock)
          Dxyz_D = Dxyz_D*(CellSize_DB(z_,iBlock)/sum(Dxyz_D**2))
          DcoordDxyz_DDFD(z_,:,i,j,1,z_) = Dxyz_D
       end do; end do
    end if
    if(nK==1)RETURN
    if(neiLTop(iBlock)==NoBlk)then
       do j=1,nJ; do i=1,nI
          Dxyz_D = Xyz_DGB(:,i,j,nK + 1,iBlock) - Xyz_DGB(:,i,j,nK,iBlock)
          Dxyz_D = Dxyz_D*(CellSize_DB(z_,iBlock)/sum(Dxyz_D**2))
          DcoordDxyz_DDFD(z_,:,i,j,nK +1,z_) = Dxyz_D
       end do; end do
    end if

  end subroutine set_block_jacobian_face

  !============================================================================
  subroutine get_face_gradient_field(iDir, i, j, k, iBlock, nField, &
       IsNewBlock, Var_IG, FaceGrad_DI)

    ! calculate the gradient FaceGrad_DI of field Var_IG 
    ! on face iDir of cell i, j, k of block iBlock

    use ModParallel,   ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot
    use BATL_lib,      ONLY: IsCartesianGrid, nDim, Dim1_, Dim2_, Dim3_, &
         CellSize_DB, DiLevelNei_IIIB, MaxDim, x_, y_, z_, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    integer, intent(in) :: iDir, i, j, k, iBlock, nField
    logical, intent(inout) :: IsNewBlock
    real, intent(inout) :: Var_IG(nField,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out) :: FaceGrad_DI(nDim,nField)

    integer :: iL, iR, jL, jR, kL, kR, iField
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz
    Real :: InvDx, InvDy, InvDz
    real, allocatable :: Var1_IG(:,:,:,:)

    character(len=*), parameter:: NameSub='get_face_gradient_field'
    !--------------------------------------------------------------------------
    InvDx = 1.0/CellSize_DB(x_,iBlock)
    InvDy = 1.0/CellSize_DB(y_,iBlock)
    InvDz = 1.0/CellSize_DB(z_,iBlock)

    if(IsNewBlock)then
       allocate(Var1_IG(nField,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
       call set_block_field3(iBlock, nField, Var1_IG, Var_IG)
       deallocate(Var1_IG)
       if(.not.IsCartesianGrid) call set_block_jacobian_face(iBlock)

       IsNewBlock = .false.
    end if

    ! Central difference with averaging in orthogonal direction
    iR = i+1; iL = i-1; 
    jR = j+1; jL = j-1; 
    kR = k+1; kL = k-1; 

    Ax = -0.25*InvDx; Bx = 0.0; Cx = +0.25*InvDx
    Ay = -0.25*InvDy; By = 0.0; Cy = +0.25*InvDy
    Az = -0.25*InvDz; Bz = 0.0; Cz = +0.25*InvDz

    if(i==1)then
       if(NeiLeast(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. DiLevelNei_IIIB(-1,-1, 0,iBlock)==-1) .or. &
            (j==nJ+1 .and. DiLevelNei_IIIB(-1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. DiLevelNei_IIIB(-1, 0,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. DiLevelNei_IIIB(-1, 0, 1,iBlock)==-1)) &
            )then
          iL = i+1; iR = i+2; Ax=InvDx; Bx=-0.75*InvDx; Cx=-0.25*InvDx
       end if
    elseif((i==nI+1 .or. i==nI.and.iDir/=Dim1_) .and. NeiLwest(iBlock)==-1 .or. &
         i==nI .and. ((iDir==y_ .and. &
         (j==1    .and. DiLevelNei_IIIB( 1,-1, 0,iBlock)==-1) .or. &
         (j==nJ+1 .and. DiLevelNei_IIIB( 1, 1, 0,iBlock)==-1)) &
         .or.         (iDir==z_ .and. &
         (k==1    .and. DiLevelNei_IIIB( 1, 0,-1,iBlock)==-1) .or. &
         (k==nK+1 .and. DiLevelNei_IIIB( 1, 0, 1,iBlock)==-1))) &
         )then
       iL = i-1; iR = i-2; Ax=-InvDx; Bx=0.75*InvDx; Cx=0.25*InvDx
    end if

    if(j==1)then
       if(NeiLsouth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. DiLevelNei_IIIB(-1,-1, 0,iBlock)==-1) .or. &
            (i==nI+1 .and. DiLevelNei_IIIB( 1,-1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. DiLevelNei_IIIB( 0,-1,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. DiLevelNei_IIIB( 0,-1, 1,iBlock)==-1)) &
            )then
          jL = j+1; jR = j+2; Ay=InvDy; By=-0.75*InvDy; Cy=-0.25*InvDy
       end if
    elseif((j==nJ+1 .or. j==nJ.and.iDir/=y_) .and. NeiLnorth(iBlock)==-1 .or. &
         j==nJ .and. ((iDir==x_ .and. &
         (i==1    .and. DiLevelNei_IIIB(-1, 1, 0,iBlock)==-1) .or. &
         (i==nI+1 .and. DiLevelNei_IIIB( 1, 1, 0,iBlock)==-1)) &
         .or.         (iDir==z_ .and. &
         (k==1    .and. DiLevelNei_IIIB( 0, 1,-1,iBlock)==-1) .or. &
         (k==nK+1 .and. DiLevelNei_IIIB( 0, 1, 1,iBlock)==-1)))&
         )then
       jL = j-1; jR = j-2; Ay=-InvDy; By=0.75*InvDy; Cy=0.25*InvDy
    end if

    if(k==1)then
       if(NeiLbot(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. DiLevelNei_IIIB(-1, 0,-1,iBlock)==-1) .or. &
            (i==nI+1 .and. DiLevelNei_IIIB( 1, 0,-1,iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. DiLevelNei_IIIB( 0,-1,-1,iBlock)==-1) .or. &
            (j==nJ+1 .and. DiLevelNei_IIIB( 0, 1,-1,iBlock)==-1)) &
            )then
          kL = k+1; kR = k+2; Az=InvDz; Bz=-0.75*InvDz; Cz=-0.25*InvDz
       end if
    elseif((k==nK+1 .or. k==nK.and.iDir/=z_) .and. NeiLtop(iBlock)==-1 .or. &
         k==nK .and. ((iDir==x_ .and. &
         (i==1    .and. DiLevelNei_IIIB(-1, 0, 1,iBlock)==-1) .or. &
         (i==nI+1 .and. DiLevelNei_IIIB( 1, 0, 1,iBlock)==-1)) &
         .or.         (iDir==y_ .and. &
         (j==1    .and. DiLevelNei_IIIB( 0,-1, 1,iBlock)==-1) .or. &
         (j==nJ+1 .and. DiLevelNei_IIIB( 0, 1, 1,iBlock)==-1))) &
         )then
       kL = k-1; kR = k-2; Az=-InvDz; Bz=0.75*InvDz; Cz=0.25*InvDz
    end if

    do iField = 1, nField
       ! Use central difference to get gradient at face
       select case(iDir)
       case(1)
          FaceGrad_DI(Dim1_,iField) = &
               InvDx*(Var_IG(iField,i,j,k) - Var_IG(iField,i-1,j,k))

          if(nJ > 1) FaceGrad_DI(Dim2_,iField) = &
               + Ay*(Var_IG(iField,i-1,jL,k) + Var_IG(iField,i,jL,k)) &
               + By*(Var_IG(iField,i-1,j ,k) + Var_IG(iField,i,j ,k)) &
               + Cy*(Var_IG(iField,i-1,jR,k) + Var_IG(iField,i,jR,k))

          if(nK > 1) FaceGrad_DI(Dim3_,iField) = &
                  + Az*(Var_IG(iField,i-1,j,kL) + Var_IG(iField,i,j,kL)) &
                  + Bz*(Var_IG(iField,i-1,j,k ) + Var_IG(iField,i,j,k )) &
                  + Cz*(Var_IG(iField,i-1,j,kR) + Var_IG(iField,i,j,kR))
       case(2)
          FaceGrad_DI(Dim1_,iField) = &
               + Ax*(Var_IG(iField,iL,j-1,k) + Var_IG(iField,iL,j,k)) &
               + Bx*(Var_IG(iField,i ,j-1,k) + Var_IG(iField,i ,j,k)) &
               + Cx*(Var_IG(iField,iR,j-1,k) + Var_IG(iField,iR,j,k))

          FaceGrad_DI(Dim2_,iField) = &
               InvDy*(Var_IG(iField,i,j,k) - Var_IG(iField,i,j-1,k))

          if(nK > 1) FaceGrad_DI(Dim3_,iField) = &
                  + Az*(Var_IG(iField,i,j-1,kL) + Var_IG(iField,i,j,kL)) &
                  + Bz*(Var_IG(iField,i,j-1,k ) + Var_IG(iField,i,j,k )) &
                  + Cz*(Var_IG(iField,i,j-1,kR) + Var_IG(iField,i,j,kR))
       case(3)
          FaceGrad_DI(Dim1_,iField) = &
               + Ax*(Var_IG(iField,iL,j,k-1) + Var_IG(iField,iL,j,k)) &
               + Bx*(Var_IG(iField,i ,j,k-1) + Var_IG(iField,i ,j,k)) &
               + Cx*(Var_IG(iField,iR,j,k-1) + Var_IG(iField,iR,j,k))

          FaceGrad_DI(Dim2_,iField) = &
               + Ay*(Var_IG(iField,i,jL,k-1) + Var_IG(iField,i,jL,k))  &
               + By*(Var_IG(iField,i,j ,k-1) + Var_IG(iField,i,j ,k))  &
               + Cy*(Var_IG(iField,i,jR,k-1) + Var_IG(iField,i,jR,k))

          FaceGrad_DI(Dim3_,iField) = &
               InvDz*(Var_IG(iField,i,j,k) - Var_IG(iField,i,j,k-1))
       case default
          write(*,*)'Error face index iDir=',iDir
          call stop_mpi(NameSub)
       end select

       ! multiply by the coordinate transformation matrix to obtain the
       ! cartesian gradient from the partial derivatives dScalar/dGencoord
       if(.not.IsCartesianGrid) FaceGrad_DI(:,iField) = &
            matmul(FaceGrad_DI(:,iField), &
            DcoordDxyz_DDFD(1:nDim,1:nDim,i,j,k,iDir))

    end do

  end subroutine get_face_gradient_field
  !==========================================================================

  subroutine get_face_gradient(iDir, i, j, k, iBlock, IsNewBlock, Scalar_G, &
       FaceGrad_D, UseFirstOrderBcIn) 

    ! calculate the cell face gradient of Scalar_G

    use BATL_lib,      ONLY: IsCartesianGrid
    use ModMain,       ONLY: x_, y_, z_
    use ModParallel,   ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot, NoBlk
    use BATL_lib,      ONLY: CellSize_DB, DiLevelNei_IIIB

    integer, intent(in) :: iDir, i, j, k, iBlock
    logical, intent(inout) :: IsNewBlock
    real, intent(inout) :: Scalar_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out) :: FaceGrad_D(3)
    logical, optional, intent(in):: UseFirstOrderBcIn
    !\
    !Limits for the cell index for the cells involoved in calculating
    !the vector components of gradient, which are parallel to the face
    !/
    integer :: iL, iR, jL, jR, kL, kR
    !\
    !Limits for the cell index for the cells involoved in calculating
    !the vector component of gradient, which is orthogonal to the face
    !/
    integer :: iD, iU, jD, jU, kD, kU
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz
    Real :: InvDx, InvDy, InvDz
    real :: Scalar1_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    !\
    ! If UseFirstOrderBc, then near the domain boundary the
    ! ghost cell value is only used to calculate the gradient
    ! in the direction across the computational domain boundary
    ! and only calculating  the gradient at the face which is
    ! a part of the computational domain boundary. Otherwise the
    ! physical cells and the ghostcells within the computational 
    ! domain are used, with eliminating the out-of-domain ghost cell 
    ! from the interpolation stencil.
    ! 
    !/
    logical :: UseFirstOrderBc
    !--------------------------------------------------------------------------
    InvDx = 1.0/CellSize_DB(x_,iBlock)
    InvDy = 1.0/CellSize_DB(y_,iBlock)
    InvDz = 1.0/CellSize_DB(z_,iBlock)

    if(IsNewBlock)then
       call set_block_field3(iBlock, 1, Scalar1_G, Scalar_G)
       if(.not.IsCartesianGrid) &
            call set_block_jacobian_face(iBlock,UseFirstOrderBcIn)
       IsNewBlock = .false.
    end if
    if(present(UseFirstOrderBcIn))then
       UseFirstOrderBc = UseFirstOrderBcIn
    else
       UseFirstOrderBc = .false.
    end if

    ! Central difference with averaging in orthogonal direction
    iR = i+1; iL = i-1; iD = i - 1; iU = i;
    jR = j+1; jL = j-1; jD = j - 1; jU = j;
    kR = k+1; kL = k-1; kD = k - 1; kU = k

    Ax = -0.25*InvDx; Bx = 0.0; Cx = +0.25*InvDx
    Ay = -0.25*InvDy; By = 0.0; Cy = +0.25*InvDy
    Az = -0.25*InvDz; Bz = 0.0; Cz = +0.25*InvDz

    if(i==1)then
       if(NeiLeast(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. DiLevelNei_IIIB(-1,-1, 0,iBlock)==-1) .or. &
            (j==nJ+1 .and. DiLevelNei_IIIB(-1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. DiLevelNei_IIIB(-1, 0,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. DiLevelNei_IIIB(-1, 0, 1,iBlock)==-1)) &
            )then
          iL = i+1; iR = i+2; Ax=InvDx; Bx=-0.75*InvDx; Cx=-0.25*InvDx
       elseif(UseFirstOrderBc.and.NeiLeast(iBlock)==NoBlk)then
          iL = i; iD = i; Ax = 0.0; Bx = -0.50*InvDx; Cx = 0.50*InvDx   
       end if
    elseif((i==nI+1 .or. i==nI.and.iDir/=x_) .and. NeiLwest(iBlock)==-1 .or. &
         i==nI .and. ((iDir==y_ .and. &
         (j==1    .and. DiLevelNei_IIIB( 1,-1, 0,iBlock)==-1) .or. &
         (j==nJ+1 .and. DiLevelNei_IIIB( 1, 1, 0,iBlock)==-1)) &
         .or.         (iDir==z_ .and. &
         (k==1    .and. DiLevelNei_IIIB( 1, 0,-1,iBlock)==-1) .or. &
         (k==nK+1 .and. DiLevelNei_IIIB( 1, 0, 1,iBlock)==-1))) &
         )then
       iL = i-1; iR = i-2; Ax=-InvDx; Bx=0.75*InvDx; Cx=0.25*InvDx
    elseif(UseFirstOrderBc.and.(i==nI+1 .or. i==nI.and.iDir/=x_)&
         .and. NeiLwest(iBlock)==NoBlk)then
       iR = i; iU = i-1; Ax =-0.50*InvDx; Bx = 0.50*InvDx; Cx = 0.0
    end if

    if(j==1)then
       if(NeiLsouth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. DiLevelNei_IIIB(-1,-1, 0,iBlock)==-1) .or. &
            (i==nI+1 .and. DiLevelNei_IIIB( 1,-1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. DiLevelNei_IIIB( 0,-1,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. DiLevelNei_IIIB( 0,-1, 1,iBlock)==-1)) &
            )then
          jL = j+1; jR = j+2; Ay=InvDy; By=-0.75*InvDy; Cy=-0.25*InvDy
       elseif(UseFirstOrderBc.and.NeiLsouth(iBlock)==NoBlk)then
          jL = i; jD = j; Ay = 0.0; By = -0.50*InvDy; Cy = 0.50*InvDy 
       end if
    elseif((j==nJ+1 .or. j==nJ.and.iDir/=y_) .and. NeiLnorth(iBlock)==-1 .or. &
         j==nJ .and. ((iDir==x_ .and. &
         (i==1    .and. DiLevelNei_IIIB(-1, 1, 0,iBlock)==-1) .or. &
         (i==nI+1 .and. DiLevelNei_IIIB( 1, 1, 0,iBlock)==-1)) &
         .or.         (iDir==z_ .and. &
         (k==1    .and. DiLevelNei_IIIB( 0, 1,-1,iBlock)==-1) .or. &
         (k==nK+1 .and. DiLevelNei_IIIB( 0, 1, 1,iBlock)==-1)))&
         )then
       jL = j-1; jR = j-2; Ay=-InvDy; By=0.75*InvDy; Cy=0.25*InvDy
    elseif(UseFirstOrderBc.and.(j==nJ+1 .or. j==nJ.and.iDir/=y_)&
         .and. NeiLnorth(iBlock)==NoBlk)then
       jR = j; jU = j-1; Ay =-0.50*InvDy; By = 0.50*InvDy; Cy = 0.0
    end if

    if(k==1)then
       if(NeiLbot(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. DiLevelNei_IIIB(-1, 0,-1,iBlock)==-1) .or. &
            (i==nI+1 .and. DiLevelNei_IIIB( 1, 0,-1,iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. DiLevelNei_IIIB( 0,-1,-1,iBlock)==-1) .or. &
            (j==nJ+1 .and. DiLevelNei_IIIB( 0, 1,-1,iBlock)==-1)) &
            )then
          kL = k+1; kR = k+2; Az=InvDz; Bz=-0.75*InvDz; Cz=-0.25*InvDz
       elseif(UseFirstOrderBc.and.NeiLbot(iBlock)==NoBlk)then
          kL = k; kD = k; Az = 0.0; Bz = -0.50*InvDz; Cz = 0.50*InvDz 
       end if
    elseif((k==nK+1 .or. k==nK.and.iDir/=z_) .and. NeiLtop(iBlock)==-1 .or. &
         k==nK .and. ((iDir==x_ .and. &
         (i==1    .and. DiLevelNei_IIIB(-1, 0, 1,iBlock)==-1) .or. &
         (i==nI+1 .and. DiLevelNei_IIIB( 1, 0, 1,iBlock)==-1)) &
         .or.         (iDir==y_ .and. &
         (j==1    .and. DiLevelNei_IIIB( 0,-1, 1,iBlock)==-1) .or. &
         (j==nJ+1 .and. DiLevelNei_IIIB( 0, 1, 1,iBlock)==-1))) &
         )then
       kL = k-1; kR = k-2; Az=-InvDz; Bz=0.75*InvDz; Cz=0.25*InvDz
    elseif(UseFirstOrderBc.and.(k==nK+1 .or. k==nK.and.iDir/=z_)& 
         .and. NeiLtop(iBlock)==NoBlk)then
       kR = k; kU = k-1; Az =-0.50*InvDz; Bz = 0.50*InvDz; Cz = 0.0
    end if

    ! Use central difference to get gradient at face
    select case(iDir)
    case(x_)
       FaceGrad_D(x_) = InvDx*(Scalar_G(i,j,k) - Scalar_G(i-1,j,k))
       if(nJ > 1)then
          FaceGrad_D(y_) = &
               + Ay*(Scalar_G(iD,jL,k) + Scalar_G(iU,jL,k)) &
               + By*(Scalar_G(iD,j ,k) + Scalar_G(iU,j ,k)) &
               + Cy*(Scalar_G(iD,jR,k) + Scalar_G(iU,jR,k))
       else
          FaceGrad_D(y_) = 0.0
       end if
       if(nK > 1)then
          FaceGrad_D(z_) = &
               + Az*(Scalar_G(iD,j,kL) + Scalar_G(iU,j,kL)) &
               + Bz*(Scalar_G(iD,j,k ) + Scalar_G(iU,j,k )) &
               + Cz*(Scalar_G(iD,j,kR) + Scalar_G(iU,j,kR))
       else
          FaceGrad_D(z_) = 0.0
       end if
    case(y_)
       FaceGrad_D(x_) = &
            + Ax*(Scalar_G(iL,jD,k) + Scalar_G(iL,jU,k)) &
            + Bx*(Scalar_G(i ,jD,k) + Scalar_G(i ,jU,k)) &
            + Cx*(Scalar_G(iR,jD,k) + Scalar_G(iR,jU,k))
       FaceGrad_D(y_) = InvDy*(Scalar_G(i,j,k) - Scalar_G(i,j-1,k))
       if(nK > 1)then
          FaceGrad_D(z_) = &
               + Az*(Scalar_G(i,jD,kL) + Scalar_G(i,jU,kL)) &
               + Bz*(Scalar_G(i,jD,k ) + Scalar_G(i,jU,k )) &
               + Cz*(Scalar_G(i,jD,kR) + Scalar_G(i,jU,kR))
       else
          FaceGrad_D(z_) = 0.0
       end if
    case(z_)
       FaceGrad_D(x_) = &
            + Ax*(Scalar_G(iL,j,kD) + Scalar_G(iL,j,kU)) &
            + Bx*(Scalar_G(i ,j,kD) + Scalar_G(i ,j,kU)) &
            + Cx*(Scalar_G(iR,j,kD) + Scalar_G(iR,j,kU))
       FaceGrad_D(y_) = &
            + Ay*(Scalar_G(i,jL,kD) + Scalar_G(i,jL,kU))  &
            + By*(Scalar_G(i,j ,kD) + Scalar_G(i,j ,kU))  &
            + Cy*(Scalar_G(i,jR,kD) + Scalar_G(i,jR,kU))
       FaceGrad_D(z_) = InvDz*(Scalar_G(i,j,k) - Scalar_G(i,j,k-1))
    case default
       write(*,*)'Error in get_face_gradient: iDir=',iDir
       call stop_mpi('DEBUG')
    end select

    ! multiply by the coordinate transformation matrix to obtain the
    ! cartesian gradient from the partial derivatives dScalar/dGencoord
    if(.not.IsCartesianGrid) then
       FaceGrad_D = matmul(FaceGrad_D, DcoordDxyz_DDFD(:,:,i,j,k,iDir))
   end if
  end subroutine get_face_gradient

  !============================================================================

  subroutine get_face_curl(iDir, i, j, k, iBlock, IsNewBlock, Vector_DG, &
       FaceCurl_D)

    use ModMain,      ONLY: x_, y_, z_
    use BATL_lib,     ONLY: IsCartesianGrid, IsRzGeometry
    use ModParallel,  ONLY: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot
    use BATL_lib,     ONLY: CellSize_DB, DiLevelNei_IIIB

    integer, intent(in) :: iDir, i, j, k, iBlock
    logical, intent(inout) :: IsNewBlock
    real, intent(inout) :: Vector_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out)  :: FaceCurl_D(3)

    integer :: iL, iR, jL, jR, kL, kR
    real :: Ax, Bx, Cx, Ay, By, Cy, Az, Bz, Cz
    Real :: InvDx, InvDy, InvDz
    real :: Vector1_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    !--------------------------------------------------------------------------

    InvDx = 1.0/CellSize_DB(x_,iBlock)
    InvDy = 1.0/CellSize_DB(y_,iBlock)
    InvDz = 1.0/CellSize_DB(z_,iBlock)

    if(IsNewBlock)then
       call set_block_field3(iBlock, 3, Vector1_DG, Vector_DG)
       if(.not. IsCartesianGrid) call set_block_jacobian_face(iBlock)
       IsNewBlock = .false.
    end if

    ! Central difference with averaging in orthogonal direction
    iR = i+1; iL = i-1; 
    jR = j+1; jL = j-1; 
    kR = k+1; kL = k-1; 

    Ax = -0.25*InvDx; Bx = 0.0; Cx = +0.25*InvDx
    Ay = -0.25*InvDy; By = 0.0; Cy = +0.25*InvDy
    Az = -0.25*InvDz; Bz = 0.0; Cz = +0.25*InvDz

    if(i==1)then
       if(NeiLeast(iBlock)==-1 &
            .or. (iDir==y_ .and. &
            (j==1    .and. DiLevelNei_IIIB(-1,-1, 0,iBlock)==-1) .or. &
            (j==nJ+1 .and. DiLevelNei_IIIB(-1, 1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. DiLevelNei_IIIB(-1, 0,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. DiLevelNei_IIIB(-1, 0, 1,iBlock)==-1)) &
            )then
          iL = i+1; iR = i+2; Ax=InvDx; Bx=-0.75*InvDx; Cx=-0.25*InvDx
       end if
    elseif((i==nI+1 .or. i==nI.and.iDir/=x_) .and. NeiLwest(iBlock)==-1 .or. &
         i==nI .and. ((iDir==y_ .and. &
         (j==1    .and. DiLevelNei_IIIB( 1,-1, 0,iBlock)==-1) .or. &
         (j==nJ+1 .and. DiLevelNei_IIIB( 1, 1, 0,iBlock)==-1)) &
         .or.         (iDir==z_ .and. &
         (k==1    .and. DiLevelNei_IIIB( 1, 0,-1,iBlock)==-1) .or. &
         (k==nK+1 .and. DiLevelNei_IIIB( 1, 0, 1,iBlock)==-1))) &
         )then
       iL = i-1; iR = i-2; Ax=-InvDx; Bx=0.75*InvDx; Cx=0.25*InvDx
    end if

    if(j==1)then
       if(NeiLsouth(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. DiLevelNei_IIIB(-1,-1, 0,iBlock)==-1) .or. &
            (i==nI+1 .and. DiLevelNei_IIIB( 1,-1, 0,iBlock)==-1)) &
            .or. (iDir==z_ .and. &
            (k==1    .and. DiLevelNei_IIIB( 0,-1,-1,iBlock)==-1) .or. &
            (k==nK+1 .and. DiLevelNei_IIIB( 0,-1, 1,iBlock)==-1)) &
            )then
          jL = j+1; jR = j+2; Ay=InvDy; By=-0.75*InvDy; Cy=-0.25*InvDy
       end if
    elseif((j==nJ+1 .or. j==nJ.and.iDir/=y_) .and. NeiLnorth(iBlock)==-1 .or. &
         j==nJ .and. ((iDir==x_ .and. &
         (i==1    .and. DiLevelNei_IIIB(-1, 1, 0,iBlock)==-1) .or. &
         (i==nI+1 .and. DiLevelNei_IIIB( 1, 1, 0,iBlock)==-1)) &
         .or.         (iDir==z_ .and. &
         (k==1    .and. DiLevelNei_IIIB( 0, 1,-1,iBlock)==-1) .or. &
         (k==nK+1 .and. DiLevelNei_IIIB( 0, 1, 1,iBlock)==-1)))&
         )then
       jL = j-1; jR = j-2; Ay=-InvDy; By=0.75*InvDy; Cy=0.25*InvDy
    end if

    if(k==1)then
       if(NeiLbot(iBlock)==-1 &
            .or. (iDir==x_ .and. &
            (i==1    .and. DiLevelNei_IIIB(-1, 0,-1,iBlock)==-1) .or. &
            (i==nI+1 .and. DiLevelNei_IIIB( 1, 0,-1,iBlock)==-1)) &
            .or. (iDir==y_ .and. &
            (j==1    .and. DiLevelNei_IIIB( 0,-1,-1,iBlock)==-1) .or. &
            (j==nJ+1 .and. DiLevelNei_IIIB( 0, 1,-1,iBlock)==-1)) &
            )then
          kL = k+1; kR = k+2; Az=InvDz; Bz=-0.75*InvDz; Cz=-0.25*InvDz
       end if
    elseif((k==nK+1 .or. k==nK.and.iDir/=z_) .and. NeiLtop(iBlock)==-1 .or. &
         k==nK .and. ((iDir==x_ .and. &
         (i==1    .and. DiLevelNei_IIIB(-1, 0, 1,iBlock)==-1) .or. &
         (i==nI+1 .and. DiLevelNei_IIIB( 1, 0, 1,iBlock)==-1)) &
         .or.         (iDir==y_ .and. &
         (j==1    .and. DiLevelNei_IIIB( 0,-1, 1,iBlock)==-1) .or. &
         (j==nJ+1 .and. DiLevelNei_IIIB( 0, 1, 1,iBlock)==-1))) &
         )then
       kL = k-1; kR = k-2; Az=-InvDz; Bz=0.75*InvDz; Cz=0.25*InvDz
    end if

    if(IsCartesianGrid)then
       call calc_cartesian_curl
    else
       call calc_gencoord_curl
    end if                             

  contains

    !==========================================================================

    subroutine calc_cartesian_curl

      use BATL_lib, ONLY: Xyz_DGB
      !------------------------------------------------------------------------
      select case(iDir)
      case(x_)
         FaceCurl_D(y_) = -InvDx*(Vector_DG(z_,i,j,k) - Vector_DG(z_,i-1,j,k))
         if(nK > 1) FaceCurl_D(y_) = FaceCurl_D(y_) &
              + Az*(Vector_DG(x_,i-1,j,kL) + Vector_DG(x_,i,j,kL)) &
              + Bz*(Vector_DG(x_,i-1,j,k ) + Vector_DG(x_,i,j,k )) &
              + Cz*(Vector_DG(x_,i-1,j,kR) + Vector_DG(x_,i,j,kR))

         FaceCurl_D(z_) = +InvDx*(Vector_DG(y_,i,j,k) - Vector_DG(y_,i-1,j,k))
         if(nJ > 1) FaceCurl_D(z_) = FaceCurl_D(z_) &
              - Ay*(Vector_DG(x_,i-1,jL,k) + Vector_DG(x_,i,jL,k)) &
              - By*(Vector_DG(x_,i-1,j ,k) + Vector_DG(x_,i,j ,k)) &
              - Cy*(Vector_DG(x_,i-1,jR,k) + Vector_DG(x_,i,jR,k)) 

         if(nJ > 1)then
            FaceCurl_D(x_) = &
                 + Ay*(Vector_DG(z_,i-1,jL,k) + Vector_DG(z_,i,jL,k )) &
                 + By*(Vector_DG(z_,i-1,j ,k) + Vector_DG(z_,i,j ,k )) &
                 + Cy*(Vector_DG(z_,i-1,jR,k) + Vector_DG(z_,i,jR,k ))
         else
            FaceCurl_D(x_) = 0.0
         end if
         if(nK > 1) FaceCurl_D(x_) = FaceCurl_D(x_) &
              - Az*(Vector_DG(y_,i-1,j,kL) + Vector_DG(y_,i,j ,kL)) &
              - Bz*(Vector_DG(y_,i-1,j,k ) + Vector_DG(y_,i,j ,k )) &
              - Cz*(Vector_DG(y_,i-1,j,kR) + Vector_DG(y_,i,j ,kR))

         ! Correct current for rz-geometry: Jz = Jz + Bphi/radius
         if(IsRzGeometry) FaceCurl_D(x_) = FaceCurl_D(x_) &
              + 0.5*(Vector_DG(z_,i,j,k) + Vector_DG(z_,i-1,j,k)) &
              / Xyz_DGB(y_,i,j,k,iBlock)

      case(y_)
         FaceCurl_D(x_) = &
              +InvDy*(Vector_DG(z_,i,j,k) - Vector_DG(z_,i,j-1,k))
         if(nK > 1) FaceCurl_D(x_) = FaceCurl_D(x_) &
              - Az*(Vector_DG(y_,i,j-1,kL) + Vector_DG(y_,i,j,kL)) &
              - Bz*(Vector_DG(y_,i,j-1,k ) + Vector_DG(y_,i,j,k )) &
              - Cz*(Vector_DG(y_,i,j-1,kR) + Vector_DG(y_,i,j,kR))

         FaceCurl_D(z_) = &
              -InvDy*(Vector_DG(x_,i,j,k) - Vector_DG(x_,i,j-1,k)) &
              + Ax*(Vector_DG(y_,iL,j-1,k) + Vector_DG(y_,iL,j,k)) &
              + Bx*(Vector_DG(y_,i ,j-1,k) + Vector_DG(y_,i ,j,k)) &
              + Cx*(Vector_DG(y_,iR,j-1,k) + Vector_DG(y_,iR,j,k))

         FaceCurl_D(y_) = &
              - Ax*(Vector_DG(z_,iL,j-1,k) + Vector_DG(z_,iL,j,k)) &
              - Bx*(Vector_DG(z_,i ,j-1,k) + Vector_DG(z_,i ,j,k)) &
              - Cx*(Vector_DG(z_,iR,j-1,k) + Vector_DG(z_,iR,j,k))
         if(nK > 1) FaceCurl_D(y_) = FaceCurl_D(y_) &
              + Az*(Vector_DG(x_,i,j-1,kL) + Vector_DG(x_,i,j,kL)) &
              + Bz*(Vector_DG(x_,i,j-1,k ) + Vector_DG(x_,i,j,k )) &
              + Cz*(Vector_DG(x_,i,j-1,kR) + Vector_DG(x_,i,j,kR))

         ! Correct current for rz-geometry: Jz = Jz + Bphi/radius
         if(IsRzGeometry)then
            if(Xyz_DGB(y_,i,j-1,k,iBlock)<0.0)then
               ! Just for bookkeeping. It's effect is zeroed by zero face area
               FaceCurl_D(x_) = FaceCurl_D(x_) &
                    + Vector_DG(z_,i,j,k)/Xyz_DGB(y_,i,j,k,iBlock)
            else
               FaceCurl_D(x_) = FaceCurl_D(x_) &
                    + (Vector_DG(z_,i,j,k) + Vector_DG(z_,i,j-1,k)) &
                    / (Xyz_DGB(y_,i,j,k,iBlock) + Xyz_DGB(y_,i,j-1,k,iBlock))
            end if
         end if

      case(z_)
         FaceCurl_D(x_) = &
              -InvDz*(Vector_DG(y_,i,j,k) - Vector_DG(y_,i,j,k-1)) & 
              + Ay*(Vector_DG(z_,i,jL,k-1) + Vector_DG(z_,i,jL,k)) &
              + By*(Vector_DG(z_,i,j ,k-1) + Vector_DG(z_,i,j ,k)) &
              + Cy*(Vector_DG(z_,i,jR,k-1) + Vector_DG(z_,i,jR,k))

         FaceCurl_D(y_) = &
              +InvDz*(Vector_DG(x_,i,j,k) - Vector_DG(x_,i,j,k-1)) &
              - Ax*(Vector_DG(z_,iL,j,k-1) + Vector_DG(z_,iL,j,k)) &
              - Bx*(Vector_DG(z_,i ,j,k-1) + Vector_DG(z_,i ,j,k)) &
              - Cx*(Vector_DG(z_,iR,j,k-1) + Vector_DG(z_,iR,j,k))

         FaceCurl_D(z_) = &
              + Ax*(Vector_DG(y_,iL,j,k-1) + Vector_DG(y_,iL,j,k)) &
              + Bx*(Vector_DG(y_,i ,j,k-1) + Vector_DG(y_,i ,j,k)) &
              + Cx*(Vector_DG(y_,iR,j,k-1) + Vector_DG(y_,iR,j,k)) &
              - Ay*(Vector_DG(x_,i,jL,k-1) + Vector_DG(x_,i,jL,k)) &
              - By*(Vector_DG(x_,i,j ,k-1) + Vector_DG(x_,i,j ,k)) &
              - Cy*(Vector_DG(x_,i,jR,k-1) + Vector_DG(x_,i,jR,k))

      case default
         write(*,*)'Error in calc_cartesian_curl: iDir=',iDir
         call stop_mpi('DEBUG')
      end select

    end subroutine calc_cartesian_curl

    !==========================================================================

    subroutine calc_gencoord_curl

      real :: DvectorDcoord_DD(MaxDim,MaxDim)
      !------------------------------------------------------------------------

      ! Calculate the partial derivatives dVector/dCoord
      select case(iDir)
      case(x_)
         DvectorDcoord_DD(:,1) = &
              InvDx*(Vector_DG(:,i,j,k) - Vector_DG(:,i-1,j,k))
         if(nJ > 1)then
            DvectorDcoord_DD(:,2) = &
                 + Ay*(Vector_DG(:,i-1,jL,k) + Vector_DG(:,i,jL,k)) &
                 + By*(Vector_DG(:,i-1,j ,k) + Vector_DG(:,i,j ,k)) &
                 + Cy*(Vector_DG(:,i-1,jR,k) + Vector_DG(:,i,jR,k))
         else
            DvectorDcoord_DD(:,2) = 0.0
         end if
         if(nK > 1)then
            DvectorDcoord_DD(:,3) = &
                 + Az*(Vector_DG(:,i-1,j,kL) + Vector_DG(:,i,j,kL)) &
                 + Bz*(Vector_DG(:,i-1,j,k ) + Vector_DG(:,i,j,k )) &
                 + Cz*(Vector_DG(:,i-1,j,kR) + Vector_DG(:,i,j,kR))
         else
            DvectorDcoord_DD(:,3) = 0.0
         end if

      case(y_)
         DvectorDcoord_DD(:,1) = &
              + Ax*(Vector_DG(:,iL,j-1,k) + Vector_DG(:,iL,j,k)) &
              + Bx*(Vector_DG(:,i ,j-1,k) + Vector_DG(:,i ,j,k)) &
              + Cx*(Vector_DG(:,iR,j-1,k) + Vector_DG(:,iR,j,k))
         DvectorDcoord_DD(:,2) = &
              InvDy*(Vector_DG(:,i,j,k) - Vector_DG(:,i,j-1,k))
         if(nK > 1)then
            DvectorDcoord_DD(:,3) = &
                 + Az*(Vector_DG(:,i,j-1,kL) + Vector_DG(:,i,j,kL)) &
                 + Bz*(Vector_DG(:,i,j-1,k ) + Vector_DG(:,i,j,k )) &
                 + Cz*(Vector_DG(:,i,j-1,kR) + Vector_DG(:,i,j,kR))
         else
            DvectorDcoord_DD(:,3) = 0.0
         end if

      case(z_)
         DvectorDcoord_DD(:,1) = &
              + Ax*(Vector_DG(:,iL,j,k-1) + Vector_DG(:,iL,j,k)) &
              + Bx*(Vector_DG(:,i ,j,k-1) + Vector_DG(:,i ,j,k)) &
              + Cx*(Vector_DG(:,iR,j,k-1) + Vector_DG(:,iR,j,k))
         DvectorDcoord_DD(:,2) = &
              + Ay*(Vector_DG(:,i,jL,k-1) + Vector_DG(:,i,jL,k)) &
              + By*(Vector_DG(:,i,j ,k-1) + Vector_DG(:,i,j ,k)) &
              + Cy*(Vector_DG(:,i,jR,k-1) + Vector_DG(:,i,jR,k))
         DvectorDcoord_DD(:,3) = &
              InvDz*(Vector_DG(:,i,j,k) - Vector_DG(:,i,j,k-1))
      end select

      ! Curl_x = Dvector_z/Dy - Dvector_y/Dz
      FaceCurl_D(x_) = &
           + sum(DvectorDcoord_DD(z_,:)*DcoordDxyz_DDFD(:,y_,i,j,k,iDir)) &
           - sum(DvectorDcoord_DD(y_,:)*DcoordDxyz_DDFD(:,z_,i,j,k,iDir))

      ! Curl_y = Dvector_x/Dz - Dvector_z/Dx
      FaceCurl_D(y_) = &
           + sum(DvectorDcoord_DD(x_,:)*DcoordDxyz_DDFD(:,z_,i,j,k,iDir)) &
           - sum(DvectorDcoord_DD(z_,:)*DcoordDxyz_DDFD(:,x_,i,j,k,iDir))

      ! Curl_z = Dvector_y/Dx - Dvector_x/Dy
      FaceCurl_D(z_) = &
           + sum(DvectorDcoord_DD(y_,:)*DcoordDxyz_DDFD(:,x_,i,j,k,iDir)) &
           - sum(DvectorDcoord_DD(x_,:)*DcoordDxyz_DDFD(:,y_,i,j,k,iDir))

    end subroutine calc_gencoord_curl

  end subroutine get_face_curl

end module ModFaceGradient
