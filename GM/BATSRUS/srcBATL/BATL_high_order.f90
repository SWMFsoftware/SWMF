!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_high_order
  implicit none

  save

contains

  !======================================================================
  real function correct_face_value(FaceValue, CellValue_I)

    ! FaceValue is at cell face. CellValue_I are cell centered.
    ! Return 6th order approximation

    real, intent(in):: CellValue_I(4), FaceValue
    real:: Der2, Der4
    real, parameter:: c1over6 = 1./6, c1over180 = 1./180
    !--------------------------------------------------------------------
    Der2 = c1over6*(CellValue_I(2) - 2*FaceValue + CellValue_I(3))
    Der4 = c1over180*(16*FaceValue - &
         9*(CellValue_I(2) + CellValue_I(3)) + &
         CellValue_I(1) + CellValue_I(4))
    correct_face_value = FaceValue - Der2 + Der4

  end function correct_face_value
  !======================================================================

  real function calc_center_first_derivate(CellValue_I, DxIn)

    ! Calculate df/dx at x=x_i with f(k), where k = i-3,i-2 ... i+3.
    ! Directly combine CellValue_I to get df/dx can save some 
    ! computation compare to the approach in this subroutine.

    real, intent(in):: CellValue_I(7)
    real, optional, intent(in):: DxIn
    real:: Dx
    real:: FaceL, FaceR ! FaceL = F_{i-1/2}, FaceR = F_{i+1/2}
    real:: CorrectedFaceL, CorrectedFaceR
    real, parameter:: c3over256 = 3./256, c25over256 = 25./256, &
         c150over256 = 150./256
    integer, parameter:: i = 4
    ! ----------------------------------------------------------------------

    Dx = 1.0
    if(present(DxIn)) Dx = DxIn

    FaceL = calc_face_value(CellValue_I(1:6))
    FaceR = calc_face_value(CellValue_I(2:7))

    CorrectedFaceL = correct_face_value(FaceL, CellValue_I(i-2:i+1))
    CorrectedFaceR = correct_face_value(FaceR, CellValue_I(i-1:i+2)) 

    calc_center_first_derivate = &
         (CorrectedFaceR - CorrectedFaceL)/Dx
  end function calc_center_first_derivate
  !======================================================================

  real function calc_face_value(CellValue_I, DoLimitIn)
    ! Calculate f_{i+1/2} with f(k), where k = i-2,i-1 ... i+3
    real, intent(in):: CellValue_I(6)
    logical, optional, intent(in):: DoLimitIn
    logical:: DoLimit
    real, parameter:: c3over256 = 3./256, c25over256 = 25./256, &
         c150over256 = 150./256
    real:: FaceValue
    real:: Distance_I(4) = (/-1.5, -0.5, 0.5, 1.5/)
    !----------------------------------------------------------------------

    DoLimit = .false. 
    if(present(DoLimitIn)) DoLimit = DoLimitIn

    FaceValue = c3over256*(CellValue_I(1) + CellValue_I(6)) - &
         c25over256*(CellValue_I(2) + CellValue_I(5)) + &
         c150over256*(CellValue_I(3) + CellValue_I(4))

    calc_face_value = FaceValue
    if(DoLimit) calc_face_value = &
         limit_interpolation(FaceValue, CellValue_I(2:5), Distance_I)
  end function calc_face_value

  !===========================================================================

  real function limit_interpolation(FaceOrig, CellValue_I, Distance_I)
    ! This Limiter works for ununiform grid interpolation. 
    ! See (2.18) in 'Accurate Monotonicity-Preserving Schemes with 
    ! Runge-Kutta Time Stepping' by A. Suresh & H. T. Huynh (1997)

    real, intent(in):: FaceOrig, CellValue_I(4), Distance_I(4)
    real:: FaceL, FaceR, FaceAV, FaceMD, FaceMin, FaceMax
    !----------------------------------------------------------------------

    FaceAV = two_points_interpolation(CellValue_I(2:3), Distance_I(2:3)) 
    FaceL = two_points_interpolation(CellValue_I(1:2), Distance_I(1:2))
    FaceR = two_points_interpolation(CellValue_I(3:4), Distance_I(3:4))
    
    FaceMD = median(FaceAV, FaceL, FaceR)
    FaceMin = min(FaceMD, CellValue_I(2), CellValue_I(3))
    FaceMax = max(FaceMD, CellValue_I(2), CellValue_I(3))

    limit_interpolation = median(FaceOrig, FaceMin, FaceMax)
    
  end function limit_interpolation

  !======================================================================

  real function median(a,b,c)
    real, intent(in):: a, b, c
    median = max(min(a, max(b,c)), min(b,c))
  end function median
  !======================================================================

  real function two_points_interpolation(Cell_I, Distance_I)
    ! Cell_I(i) is at xi, calculate the value at x=0. 
    ! Distance_I(i) = xi - x0
    real, intent(in):: Cell_I(2), Distance_I(2)
    real:: c1, c2
    !----------------------------------------------------------------------
    c1 = Distance_I(2)/(Distance_I(2) - Distance_I(1))
    c2 = -Distance_I(1)/(Distance_I(2) - Distance_I(1))
    two_points_interpolation = c1*Cell_I(1) + c2*Cell_I(2)
  end function two_points_interpolation

  !======================================================================
  subroutine get_ghost_for_coarse_blk(CoarseCell, FineCell_III, Ghost_I)
    ! For 2D: 
    !         _________________________________
    !         | u1|   |   |   |   |   |   |   |
    !         |___|___|___|___|___|___|___|___|
    !         | u2|   |   |   |   |   |   |   |
    !  _______|___|___|___|___|___|___|___|___|
    ! |       | u3| u7|   |   |   |   |   |   |
    ! |Coarse |__G1___|__G2___|__G3___|___|___|
    ! | u0    | u4| u8|   |   |   |   |   |   |
    ! |_______|___|___|___|___|___|___|___|___|
    !         | u5|   |   |   |   |   |   |   |
    !         |___|___|___|___|___|___|___|___|
    !         | u6|   |   |   |   |   |   |   |
    !         |___|___|___|___|___|___|___|___|
   

    ! First calculate the face value (f1) between u3 and u4 with u1, u2...u6.
    ! Face value between u7 and u8 (f2) and more face value can be got in the
    ! same way. 8 face values are needed. 

    ! Use u0, f1, f2, f3, f4 to calculate the ghost cell G1. 
    ! Use f1, f2, f3, f4, f5, f6 to interpolate G2
    ! Use f3, f4...f7,f8 to interpolate the third ghost cell (G3). 

    use BATL_size, ONLY: nJ, nk

    ! In 2D there is only 1 cell in the K direction
    integer, parameter:: k6_ = min(nK, 6), j6_ = min(nJ,6)
    
    real, intent(in) :: CoarseCell              ! value of coarse neighbor cell
    real, intent(in) :: FineCell_III(8,j6_,k6_) ! value of local fine cells
    real, intent(out):: Ghost_I(3)              ! coarse ghost cells for neighbor

    ! Local variables
    real:: FineCell_I(8), FaceValue_I(6)
    real:: Ghost, Cell_I(4)
    real:: CellInput_I(6), CellInput_II(6,6)
    
    ! Distances of coarse neighbor and fine cells from 1st ghost cell
    real:: Distance_I(4)=(/-4,-1,1,3/)

    ! Interpolation coefficients for G1
    real, parameter:: c1=-1./63, c2=5./12, c3 = 3./4, c4=-5./28, c5=1./36

    logical, parameter:: DoLimit = .true. ! change only for debugging
    integer:: i, j, k
    !----------------------------------------------------------------------

    ! Integerpolate fine cell centers to line connecting coarse cell centers
    if(nK == 1) then  ! 2D resolution change       
       do i = 1, 8          
          ! Use a temporary variable CellInput_I to avoid compile error. 
          CellInput_I(1:j6_) = FineCell_III(i,:,1)
          FineCell_I(i) = calc_face_value(CellInput_I, DoLimit)
       enddo
    else   
       ! 3D resolution change. Need more tests to make sure it works!!
       do i = 1, 8
          CellInput_II(1:j6_,1:k6_) = FineCell_III(I,:,:)
          FineCell_I(i) = &
               calc_edge_value(CellInput_II,Dolimit)
       enddo 
    endif

    ! High order interpolation for first ghost cell using coarse neighbor
    Ghost = c1*CoarseCell + c2*FineCell_I(1) + c3*FineCell_I(2) + &
         c4*FineCell_I(3) + c5*FineCell_I(4)

    ! Limit value if requested
    if(DoLimit) then
       Cell_I(1)   = CoarseCell
       Cell_I(2:4) = FineCell_I(1:3)
       Ghost_I(1)  = limit_interpolation(Ghost, Cell_I, Distance_I)
    else
       Ghost_I(1) = Ghost
    endif

    ! Interpolate G2 and G3 from fine edge values
    Ghost_I(2) = calc_face_value(FineCell_I(1:6), DoLimit)
    Ghost_I(3) = calc_face_value(FineCell_I(3:8), DoLimit)

  end subroutine get_ghost_for_coarse_blk
  !======================================================================
  real function calc_edge_value(CellValue_II,DoLimitIn)
    ! For 3D, need more tests. 
    real, intent(in) :: CellValue_II(6,6)
    logical, optional, intent(in):: DoLimitIn
    logical:: DoLimit
    real:: CellValue_I(6)
    integer:: i, j
    integer:: iBegin=1, iEnd=6
    !----------------------------------------------------------------------
    DoLimit = .true. 
    if(present(DoLimitIn)) DoLimit = DoLimitIn

    do i = iBegin, iEnd       
       CellValue_I(i) = calc_face_value(CellValue_II(i,:),DoLimit)
    enddo
    calc_edge_value = calc_face_value(CellValue_I,DoLimit)
  end function calc_edge_value
  !======================================================================

  subroutine get_ghost_for_fine_blk(CoarseCell_II, FineCell_I, Ghost_I)
    ! 2D: 
    ! __________________________
    ! |        |       |       |
    ! |        |       |  jmm  |
    ! |        |       |       |
    ! |________|_______|_______|
    ! |        |       |       |
    ! |        |       |  jm   |
    ! |        |       |       |
    ! |________|_______|_______|_____________
    ! |        |       |       |   |    |   |
    ! |        |       |  j2   |___|____|___|
    ! |   *    |  *  ? |?  * ? |u0 | u1 |u2 |
    ! |________|_______|_______|___|____|___|
    ! |        |       |       |
    ! |        |       |  jp   |
    ! |        |       |       |
    ! |________|_______|_______|
    ! |        |       |       | 
    ! |        |       | jpp   |
    ! |        |       |       |
    ! | _______|_______|_______|

    ! The 3*5 cells represent CoarseCell_II. u0, u1, and u2 are FineCell_I.
    ! '?' are the ghost values needed. 
    ! First calculate the values represented by '*', and then use these star 
    ! values and u0, u1, u2 to calculate the ghost cells. 

    real, intent(in):: CoarseCell_II(5,3), FineCell_I(3)
    real, intent(out):: Ghost_I(3)

    real:: CoarseCell_I(3)
    integer:: i
    logical:: DoLimit = .true.
    DoLimit = .true. 
    do i = 1, 3       
       CoarseCell_I(i) = &
            interpolate_in_coarse_blk(CoarseCell_II(:,i), DoLimitIn=DoLimit)
    enddo

    call interpolate_ghost_for_fine_blk

  contains
    !----------------------------------------------------------------------
    subroutine interpolate_ghost_for_fine_blk
      real, parameter:: c11=-4./231, c12=4./7,c13=5./7, c14=-1./3, c15=5./77
      real, parameter:: c21=-9./572, c22=1./6,c23=1.05, c24=-3./11, c25=14./195
      real, parameter:: c31=-9./286, c32=5./7,c33=0.5, c34=-20./77, c35=1./13
      real:: Ghost, Cell_I(4), Distance_I(4)

      !-----------------
      Ghost = c11*CoarseCell_I(2) + c12*CoarseCell_I(1) + &
           c13*FineCell_I(1) + c14*FineCell_I(2) + c15*FineCell_I(3)

      if(DoLimit) then
         Distance_I(1) = -5; Distance_I(2) = -1
         Distance_I(3) = 2; Distance_I(4) = 4
         Cell_I(1) = CoarseCell_I(2); Cell_I(2) = CoarseCell_I(1)
         Cell_I(3) = FineCell_I(1); Cell_I(4) = FineCell_I(2)

         Ghost_I(1) = limit_interpolation(Ghost, Cell_I, Distance_I)
      else
         Ghost_I(1) = Ghost
      endif
      !---------------------

      Ghost = c21*CoarseCell_I(3) + c22*CoarseCell_I(2) + &
           c23*CoarseCell_I(1) + c24*FineCell_I(1) + c25*FineCell_I(2)

      if(DoLimit) then
         Distance_I(1) = -7; Distance_I(2) = -3
         Distance_I(3) = 1; Distance_I(4) = 4
         Cell_I(1) = CoarseCell_I(3); Cell_I(2) = CoarseCell_I(2)
         Cell_I(3) = CoarseCell_I(1); Cell_I(4) = FineCell_I(1)

         Ghost_I(2) = limit_interpolation(Ghost, Cell_I, Distance_I)
      else
         Ghost_I(2) = Ghost
      endif
      !----------------------

      Ghost = c31*CoarseCell_I(3) + c32*CoarseCell_I(2) + &
           c33*CoarseCell_I(1) + c34*FineCell_I(1) + c35*FineCell_I(2)

      if(DoLimit) then
         Distance_I(1) = -5; Distance_I(2) = -1
         Distance_I(3) = 3; Distance_I(4) = 6
         Cell_I(1) = CoarseCell_I(3); Cell_I(2) = CoarseCell_I(2)
         Cell_I(3) = CoarseCell_I(1); Cell_I(4) = FineCell_I(1)

         Ghost_I(3) = limit_interpolation(Ghost, Cell_I, Distance_I)
      else
         Ghost_I(3) = Ghost
      endif
    end subroutine interpolate_ghost_for_fine_blk

  end subroutine get_ghost_for_fine_blk
  !======================================================================
  real function interpolate_in_coarse_blk(Cell_I, DoLimitIn)
    real, intent(in):: Cell_I(5)
    logical, optional, intent(in):: DoLimitIn
    logical:: DoLimit
    real, parameter:: cpp=-45./2048, cp=105./512, &
         c0=945./1024, cm=-63./512, cmm=35./2048
    integer, parameter:: Ipp_=1, Ip_=2, I_=3, Im_=4, Imm_=5
    real:: Temp, Distance_I(4)=(/-7,-3,1,5/)
    !----------------------------------------------------------------------

    DoLimit = .false. 
    if(present(DoLimitIn)) DoLimit = DoLimitIn

    Temp = cpp*Cell_I(Ipp_) + cp*Cell_I(Ip_) + &
         c0*Cell_I(I_) + cm*Cell_I(Im_) + cmm*Cell_I(Imm_)        
    if(DoLimit) then
       interpolate_in_coarse_blk = limit_interpolation(Temp, Cell_I(Ipp_:Im_), Distance_I)
    else
       interpolate_in_coarse_blk = Temp
    endif

  end function interpolate_in_coarse_blk
  !======================================================================

  subroutine calc_high_ghost_for_fine_blk(iBlock, nVar, Field1_VG, Field_VG,&
       neiLeast, neiLwest, neiLsouth, neiLnorth, neiLtop, neiLbot)
    ! Works for 2D.
    ! The 3D part is also partially implemented, but need more tests. 
    
    use BATL_tree, ONLY: DiLevelNei_IIIB, iNodeNei_IIIB, iTree_IA, block_, &
         unset_, Proc_
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         nI, nJ, nK, i0_,i2_,j0_, j2_, nJp1_, nJm1_, k0_, k2_, nKp1_, nKm1_, &
         jm2_,jm1_,nJm2_,nJm1_,nJp2_,nJp3_,km2_,km1_,nKm2_,nKm1_,nKp2_,&
         nKp3_,im2_,im1_,nIm2_,nIm1_,nIp1_,nIp2_,nIp3_, i3_,j3_,k3_
    
    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: Field1_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(inout) :: Field_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    integer, intent(in):: neiLeast, neiLwest, neiLsouth, &
         neiLnorth, neiLtop, neiLbot

    integer :: i1, j1, k1, i2, j2, k2
    !integer :: iL, iR, jL, jR, kL, kR
    integer :: ip, im, jp, jm, kp, km, iVar, jpp, jmm, ipp, imm, kpp, kmm

    real :: FieldCoarse_VII(nVar,5,3), FieldFine_VI(nVar,3),Ghost_I(3)
    integer, parameter:: Ipp_=1, Ip_=2, I_=3, Im_=4, Imm_=5

    integer:: iNode1, iNode2, iNode3, iNode0
    !-------------------------------------------------------------------------

    Field1_VG = Field_VG

    ! Do six faces
    if(NeiLeast == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2;
          do k2 = k1,k1+min(1,nK-1); do j2 = j1,j1+1
             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             jpp = 7*j2 - 6*j1 -3; jmm = 8*j1 - 7*j2 + 4

             iNode0 = iNodeNei_IIIB(0,0,1,iBlock)
             iNode1 = iNodeNei_IIIB(0,1,1,iBlock)
             iNode2 = iNodeNei_IIIB(0,2,1,iBlock)
             iNode3 = iNodeNei_IIIB(0,3,1,iBlock)

             if(DiLevelNei_IIIB(-1,-1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,-1,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                ! Be careful what are stored in the corner ghost cells!!
                jpp = max(-1, jpp)
                jmm = max(-1, jmm)
             endif

             if(DiLevelNei_IIIB(-1,1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,1,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                jpp = min(jpp, nJ+2)
                jmm = min(jmm, nJ+2)
             endif

             if(nK == 1)then
                kp  = 1; km  = 1
                kpp = 1; kmm = 1
             else
                kp  = 3*k2 - 2*k1 -1; km  = 4*k1 -3*k2 +2
                kpp = 7*k2 - 6*k1 -3; kmm = 8*k1 -7*k2 +4
             end if

             !      Store the values closest to the block boundary
             !                     \ /
             FieldCoarse_VII(:,Ipp_,1) = Field1_VG(:,i0_,jpp,kpp)
             FieldCoarse_VII(:,Ip_,1)  = Field1_VG(:,i0_,jp,kp)
             FieldCoarse_VII(:,I_,1)   = Field1_VG(:,i0_,j2,k2)
             FieldCoarse_VII(:,Im_,1)  = Field1_VG(:,i0_,jm,km)
             FieldCoarse_VII(:,Imm_,1) = Field1_VG(:,i0_,jmm,kmm)

             FieldCoarse_VII(:,Ipp_,2) = Field1_VG(:,im1_,jpp,kpp)
             FieldCoarse_VII(:,Ip_,2)  = Field1_VG(:,im1_,jp,kp)
             FieldCoarse_VII(:,I_,2)   = Field1_VG(:,im1_,j2,k2)
             FieldCoarse_VII(:,Im_,2)  = Field1_VG(:,im1_,jm,km)
             FieldCoarse_VII(:,Imm_,2) = Field1_VG(:,im1_,jmm,kmm)

             FieldCoarse_VII(:,Ipp_,3) = Field1_VG(:,im2_,jpp,kpp)
             FieldCoarse_VII(:,Ip_,3)  = Field1_VG(:,im2_,jp,kp)
             FieldCoarse_VII(:,I_,3)   = Field1_VG(:,im2_,j2,k2)
             FieldCoarse_VII(:,Im_,3)  = Field1_VG(:,im2_,jm,km)
             FieldCoarse_VII(:,Imm_,3) = Field1_VG(:,im2_,jmm,kmm)

             FieldFine_VI(:,1) = Field1_VG(:,1,  j2,k2)
             FieldFine_VI(:,2) = Field1_VG(:,i2_,j2,k2)
             FieldFine_VI(:,3) = Field1_VG(:,i3_,j2,k2)

             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VII(iVar,:,:), &
                     FieldFine_VI(iVar,:), Ghost_I)
                Field_VG(iVar,i0_ ,j2,k2) = Ghost_I(1)
                Field_VG(iVar,im1_,j2,k2) = Ghost_I(2)
                Field_VG(iVar,im2_,j2,k2) = Ghost_I(3)

             enddo

          end do; end do
       end do; end do
    end if

    if(NeiLwest == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2
          do k2 = k1,k1+min(1,nK-1); do j2 = j1,j1+1
             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             jpp = 7*j2 - 6*j1 -3; jmm = 8*j1 - 7*j2 + 4

             iNode0 = iNodeNei_IIIB(3,0,1,iBlock)
             iNode1 = iNodeNei_IIIB(3,1,1,iBlock)
             iNode2 = iNodeNei_IIIB(3,2,1,iBlock)
             iNode3 = iNodeNei_IIIB(3,3,1,iBlock)

             if(DiLevelNei_IIIB(1,-1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,-1,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                jpp = max(-1, jpp)
                jmm = max(-1, jmm)
             endif

             if(DiLevelNei_IIIB(1,1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,1,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                jpp = min(jpp, nJ+2)
                jmm = min(jmm, nJ+2)
             endif

             if(nK == 1)then
                kp = 1; km = 1
                kpp= 1; kmm= 1
             else
                kp  = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
                kpp = 7*k2 - 6*k1 -3; kmm = 8*k1 -7*k2 +4
             end if

             FieldCoarse_VII(:,Ipp_,1) = Field1_VG(:,nIp1_,jpp,kpp)
             FieldCoarse_VII(:,Ip_,1)  = Field1_VG(:,nIp1_,jp,kp)
             FieldCoarse_VII(:,I_,1)   = Field1_VG(:,nIp1_,j2,k2)
             FieldCoarse_VII(:,Im_,1)  = Field1_VG(:,nIp1_,jm,km)
             FieldCoarse_VII(:,Imm_,1) = Field1_VG(:,nIp1_,jmm,kmm)

             FieldCoarse_VII(:,Ipp_,2) = Field1_VG(:,nIp2_,jpp,kpp)
             FieldCoarse_VII(:,Ip_,2)  = Field1_VG(:,nIp2_,jp,kp)
             FieldCoarse_VII(:,I_,2)   = Field1_VG(:,nIp2_,j2,k2)
             FieldCoarse_VII(:,Im_,2)  = Field1_VG(:,nIp2_,jm,km)
             FieldCoarse_VII(:,Imm_,2) = Field1_VG(:,nIp2_,jmm,kmm)

             FieldCoarse_VII(:,Ipp_,3) = Field1_VG(:,nIp3_,jpp,kpp)
             FieldCoarse_VII(:,Ip_,3)  = Field1_VG(:,nIp3_,jp,kp)
             FieldCoarse_VII(:,I_,3)   = Field1_VG(:,nIp3_,j2,k2)
             FieldCoarse_VII(:,Im_,3)  = Field1_VG(:,nIp3_,jm,km)
             FieldCoarse_VII(:,Imm_,3) = Field1_VG(:,nIp3_,jmm,kmm)

             FieldFine_VI(:,1) = Field1_VG(:,nI,j2,k2)
             FieldFine_VI(:,2) = Field1_VG(:,nIm1_,j2,k2)
             FieldFine_VI(:,3) = Field1_VG(:,nIm2_,j2,k2)

             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VII(iVar,:,:), &
                     FieldFine_VI(iVar,:), Ghost_I)
                Field_VG(iVar,nIp1_,j2,k2) = Ghost_I(1)
                Field_VG(iVar,nIp2_,j2,k2) = Ghost_I(2)
                Field_VG(iVar,nIp3_,j2,k2) = Ghost_I(3)
             enddo

          end do; end do
       end do; end do
    end if

    if(NeiLsouth == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2
          do k2 = k1,k1+min(1,nK-1); do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             ipp = 7*i2 - 6*i1 -3; imm = 8*i1 - 7*i2 + 4

             iNode0 = iNodeNei_IIIB(0,0,1,iBlock)
             iNode1 = iNodeNei_IIIB(1,0,1,iBlock)
             iNode2 = iNodeNei_IIIB(2,0,1,iBlock)
             iNode3 = iNodeNei_IIIB(3,0,1,iBlock)

             if(DiLevelNei_IIIB(-1,-1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(-1,0,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                ipp = max(-1, ipp)
                imm = max(-1, imm)
             endif

             if(DiLevelNei_IIIB(1,-1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(1,0,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                ipp = min(ipp, nI+2)
                imm = min(imm, nI+2)
             endif


             if(nK == 1)then
                kp = 1; km = 1
             else
                kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
             end if

             ! Only works for 2D so far. 
             FieldCoarse_VII(:,Ipp_,1) = Field1_VG(:,ipp,j0_,k1)
             FieldCoarse_VII(:,Ip_,1)  = Field1_VG(:,ip, j0_,k1)
             FieldCoarse_VII(:,I_,1)   = Field1_VG(:,i2, j0_,k1)
             FieldCoarse_VII(:,Im_,1)  = Field1_VG(:,im, j0_,k1)
             FieldCoarse_VII(:,Imm_,1) = Field1_VG(:,imm,j0_,k1)

             FieldCoarse_VII(:,Ipp_,2) = Field1_VG(:,ipp,jm1_,k1)
             FieldCoarse_VII(:,Ip_,2)  = Field1_VG(:,ip, jm1_,k1)
             FieldCoarse_VII(:,I_,2)   = Field1_VG(:,i2, jm1_,k1)
             FieldCoarse_VII(:,Im_,2)  = Field1_VG(:,im, jm1_,k1)
             FieldCoarse_VII(:,Imm_,2) = Field1_VG(:,imm,jm1_,k1)

             FieldCoarse_VII(:,Ipp_,3) = Field1_VG(:,ipp,jm2_,k1)
             FieldCoarse_VII(:,Ip_,3)  = Field1_VG(:,ip, jm2_,k1)
             FieldCoarse_VII(:,I_,3)   = Field1_VG(:,i2, jm2_,k1)
             FieldCoarse_VII(:,Im_,3)  = Field1_VG(:,im, jm2_,k1)
             FieldCoarse_VII(:,Imm_,3) = Field1_VG(:,imm,jm2_,k1)

             FieldFine_VI(:,1) = Field1_VG(:,i2,1,k1)
             FieldFine_VI(:,2) = Field1_VG(:,i2,j2_,k1)
             FieldFine_VI(:,3) = Field1_VG(:,i2,j3_,k1)

             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VII(iVar,:,:), &
                     FieldFine_VI(iVar,:), Ghost_I)
                Field_VG(iVar,i2,j0_ ,k2) = Ghost_I(1)
                Field_VG(iVar,i2,jm1_,k2) = Ghost_I(2)
                Field_VG(iVar,i2,jm2_,k2) = Ghost_I(3)

             enddo

          end do; end do
       end do; end do
    end if

    if(NeiLnorth == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2
          do k2 = k1,k1+min(1,nK-1); do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             ipp = 7*i2 - 6*i1 -3; imm = 8*i1 - 7*i2 + 4

             iNode0 = iNodeNei_IIIB(0,3,1,iBlock)
             iNode1 = iNodeNei_IIIB(2,3,1,iBlock)
             iNode2 = iNodeNei_IIIB(2,3,1,iBlock)
             iNode3 = iNodeNei_IIIB(3,3,1,iBlock)

             if(DiLevelNei_IIIB(-1,1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(-1,0,0,iBlock) == 1 .or. &
                  iNode0 /=iNode1)) then
                ipp = max(-1, ipp)
                imm = max(-1, imm)
             endif

             if(DiLevelNei_IIIB(1,1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(1,0,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                ipp = min(ipp, nI+2)
                imm = min(imm, nI+2)
             endif

             if(nK == 1)then
                kp = 1; km = 1
             else
                kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
             end if

             ! Only works for 2D so far. 
             FieldCoarse_VII(:,Ipp_,1) = Field1_VG(:,ipp,nJp1_,k1)
             FieldCoarse_VII(:,Ip_,1)  = Field1_VG(:,ip, nJp1_,k1)
             FieldCoarse_VII(:,I_,1)   = Field1_VG(:,i2, nJp1_,k1)
             FieldCoarse_VII(:,Im_,1)  = Field1_VG(:,im, nJp1_,k1)
             FieldCoarse_VII(:,Imm_,1) = Field1_VG(:,imm,nJp1_,k1)

             FieldCoarse_VII(:,Ipp_,2) = Field1_VG(:,ipp,nJp2_,k1)
             FieldCoarse_VII(:,Ip_,2)  = Field1_VG(:,ip, nJp2_,k1)
             FieldCoarse_VII(:,I_,2)   = Field1_VG(:,i2, nJp2_,k1)
             FieldCoarse_VII(:,Im_,2)  = Field1_VG(:,im, nJp2_,k1)
             FieldCoarse_VII(:,Imm_,2) = Field1_VG(:,imm,nJp2_,k1)

             FieldCoarse_VII(:,Ipp_,3) = Field1_VG(:,ipp,nJp3_,k1)
             FieldCoarse_VII(:,Ip_,3)  = Field1_VG(:,ip, nJp3_,k1)
             FieldCoarse_VII(:,I_,3)   = Field1_VG(:,i2, nJp3_,k1)
             FieldCoarse_VII(:,Im_,3)  = Field1_VG(:,im, nJp3_,k1)
             FieldCoarse_VII(:,Imm_,3) = Field1_VG(:,imm,nJp3_,k1)

             FieldFine_VI(:,1) = Field1_VG(:,i2,nJ,   k1)
             FieldFine_VI(:,2) = Field1_VG(:,i2,nJm1_,k1)
             FieldFine_VI(:,3) = Field1_VG(:,i2,nJm2_,k1)

             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VII(iVar,:,:), &
                     FieldFine_VI(iVar,:), Ghost_I)
                Field_VG(iVar,i2,nJp1_, k2) = Ghost_I(1)
                Field_VG(iVar,i2,nJp2_, k2) = Ghost_I(2)
                Field_VG(iVar,i2,nJp3_, k2) = Ghost_I(3)
             enddo
          end do; end do
       end do; end do

    end if
  end subroutine calc_high_ghost_for_fine_blk
  !======================================================================

  subroutine correct_ghost_for_fine_blk(iBlock, nVar, Field_VG)
    ! For example: 

    !--------- 
    !|   |    |
    !|8__|_9__|
    !|   |    |
    !|6  | 7  |
    ! -------------------
    !|        | 4  |  5 |
    !|   1    |____|____|
    !|        |    |    |
    !|        | 2  | 3  |
    !--------------------

    ! For block 4, the ghost cells -2 <= i <= 0, nJ - 3 <= j <= nJ will be corredted. 
    ! For block 7, the ghost cells nI-3 <= i <= nI, -2 <= j <= 0 will be corredted. 

    use BATL_tree, ONLY: DiLevelNei_IIIB
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         nI, nJ, nK
    integer, intent(in):: iBlock, nVar
    real, intent(inout):: Field_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real,parameter::  Distance_II(4,4) = reshape((/&
         -2,-1,4,5, &
         -3,-2,3,4, &
         -4,-3,2,3, &
         -5,-4,1,2 &
         /), (/4,4/))

    real, parameter::  Coef_II(6,4) = reshape(([&
         5./21, -15./14, 12./7, 3./7,  -3./7,  5./42, &
         5./14, -10./7,  12./7, 8./7, -15./14, 2./7, &
         2./7,  -15./14, 8./7,  12./7, -10./7, 5./14,&
         5./42, -3./7,   3./7,  12./7, -15./14,5./21 &
         ]),&
         ([6,4]))

    integer:: iDir, jDir, kDir, iDir1, jDir1
    integer:: jBegin, jEnd, iBegin, iEnd
    integer:: Di, Dj, i, j, k
    integer:: DiLevel, DiLevel1, Count, iVar

    real:: Orig, CellValue_I(6)
    !----------------------------------------------------------------------

    kDir = 0
    do iDir = -1, 1; do jDir = -1, 1
       if(iDir*jDir /=0) CYCLE
       DiLevel = DiLevelNei_IIIB(iDir,jDir,kDir,iBlock)
       if(DiLevel /= 1) CYCLE

       ! DiLevel = 1
       if(jDir == 0) then ! Resolution change in x-dir
          if(iDir ==1) then
             iBegin = nI+1; iEnd = nI+3
          else
             iBegin = -2; iEnd = 0
          endif
          Di = 1

          do jDir1 = -1, 1, 2
             DiLevel1 = DiLevelNei_IIIB(iDir,jDir1,kDir,iBlock)
             if(DiLevel1 /=0) CYCLE

             if(jDir1 == -1) then
                jBegin = 1; jEnd = 4; Dj = 1
             else
                jBegin = nJ; jEnd = nJ - 3; Dj = -1
             endif

             k = 1
             do i = iBegin, iEnd, Di
                do iVar = 1, nVar
                   CellValue_I(1) = Field_VG(iVar,i,jBegin-3*Dj,k)
                   CellValue_I(2) = Field_VG(iVar,i,jBegin-2*Dj,k)
                   CellValue_I(3) = Field_VG(iVar,i,jBegin-Dj,k)
                   CellValue_I(4) = Field_VG(iVar,i,jEnd+Dj,k)
                   CellValue_I(5) = Field_VG(iVar,i,jEnd+2*Dj,k)                  
                   CellValue_I(6) = Field_VG(iVar,i,jEnd+3*Dj,k)                  

                   Count = 1
                   do j = jBegin, jEnd, Dj
                      Orig = sum(CellValue_I*Coef_II(:,Count))
                      Field_VG(iVar,i,j,k) = limit_interpolation(Orig,&
                           CellValue_I(2:5), Distance_II(:,Count))
                      Count = Count + 1                      
                   enddo
                enddo ! iVar
             enddo ! i

          enddo ! jDir1
       elseif(iDir == 0) then
          if(jDir == 1) then
             jBegin = nJ+1; jEnd = nJ+3
          else
             jBegin = -2; jEnd = 0
          endif
          Dj = 1

          do iDir1 = -1, 1, 2
             DiLevel1 = DiLevelNei_IIIB(iDir1,jDir,kDir,iBlock)
             if(DiLevel1 /=0) CYCLE

             if(iDir1 == -1) then
                iBegin = 1; iEnd = 4; Di = 1
             else
                iBegin = nI; iEnd = nI - 3; Di = -1
             endif

             k = 1
             do j = jBegin, jEnd, Dj
                do iVar = 1, nVar
                   CellValue_I(1) = Field_VG(iVar,iBegin-3*Di,j,k)
                   CellValue_I(2) = Field_VG(iVar,iBegin-2*Di,j,k)
                   CellValue_I(3) = Field_VG(iVar,iBegin-  Di,j,k)
                   CellValue_I(4) = Field_VG(iVar,iEnd  +  Di,j,k)
                   CellValue_I(5) = Field_VG(iVar,iEnd  +2*Di,j,k)
                   CellValue_I(6) = Field_VG(iVar,iEnd  +3*Di,j,k)

                   Count = 1
                   do i = iBegin, iEnd, Di
                      Orig = sum(CellValue_I*Coef_II(:,Count))
                      Field_VG(iVar,i,j,k) = limit_interpolation(Orig,&
                           CellValue_I(2:5),Distance_II(:,Count))
                      Count = Count + 1
                   enddo
                enddo ! ivar
             enddo ! j
          enddo ! idir1
       endif


    enddo; enddo
  end subroutine correct_ghost_for_fine_blk
  !======================================================================

  subroutine correct_ghost_for_coarse_blk(iBlock, nVar, Field_VG)

    !         |---------
    !         |         |
    !         |   6     | 
    !         |         |
    ! -------------------
    !|        | 4  |  5 |
    !|   1    |____|____|
    !|        |    |    |
    !|        | 2  | 3  |
    !--------------------

    ! The right top ghost cells of block 1 
    ! at j = nJ and i = nI+1,nI+2,nI+3 will be corrected.
    ! At input these cells are set as if block 6 was fine.
    ! Now we make the ghost cells high order accurate.

    use BATL_tree, ONLY: DiLevelNei_IIIB
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         nI, nJ, nK

    integer, intent(in):: iBlock, nVar
    real, intent(inout):: Field_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    integer:: iDir, jDir, kDir, i, j, k, j1
    integer:: jDir1, iDir1
    integer:: iBegin, iEnd, jBegin, jEnd
    integer:: DiYp, DiYm, DiLevel, iVar, DiLevel1
    real:: NewGhost, Orig
    real, parameter:: c1=0.05, c2=-0.3, c3=0.75
    real:: CellValue_I(4), Distance_I(4) = (/-2,-1,1,2/)

    !----------------------------------------------------------------------    
    ! Only works for x-y 2D
    kDir = 0

    ! Find finer face neighbors first (blocks 4 and 2 above)
    do iDir = -1, 1; do jDir = -1, 1
       if(iDir*jDir /=0) CYCLE ! If it is a edge or corner.
       DiLevel = DiLevelNei_IIIB(iDir,jDir,kDir,iBlock)
       if(DiLevel /= -1) CYCLE

       if(jDir == 0) then ! Resolution changes in x-dir.
          ! Find edge neighbor at the same level (block 6 above)
          do jDir1 = -1, 1, 2
             DiLevel1 = DiLevelNei_IIIB(iDir,jDir1, kDir,iBlock)
             if(DiLevel1 /=0) CYCLE ! 
             if(jDir1==-1) then
                j = 1
             else ! jDir1 = 1
                j = nJ
             endif

             if(iDir == 1) then
                iBegin = nI + 1
                iEnd = nI + 3
             elseif(iDir == -1) then
                iBegin = -2
                iEnd = 0
             else            
                call stop_mpi('wrong!! In correct_ghost_for_coarse_blk')
             endif

             k = 1
             do i = iBegin, iEnd
                do iVar = 1, nVar
                   CellValue_I(1) = Field_VG(iVar,i,j-2,k)
                   CellValue_I(2) = Field_VG(iVar,i,j-1,k)
                   CellValue_I(3) = Field_VG(iVar,i,j+1,k)
                   CellValue_I(4) = Field_VG(iVar,i,j+2,k)
                   Orig = &
                        c3*(CellValue_I(2) + CellValue_I(3)) + &
                        c2*(CellValue_I(1) + CellValue_I(4)) + &
                        c1*(Field_VG(iVar,i,j+3,k) + & ! Also copy these two values to
                        Field_VG(iVar,i,j-3,k))        ! CellValue_I in the future.
                   Field_VG(iVar,i,j,k) = limit_interpolation(Orig,&
                        CellValue_I,Distance_I)

                enddo ! ivar
             enddo ! i

          enddo ! jDir1

       elseif(iDir == 0) then ! Resolution changes in y-dir.
          do iDir1 = -1, 1, 2
             DiLevel1 = DiLevelNei_IIIB(iDir1,jDir,kDir,iBlock)
             if(DiLevel1 /=0) CYCLE

             if(iDir1 == -1) then
                i = 1
             else ! iDir1 = 1
                i = nI
             endif

             if(jDir == 1) then
                jBegin = nJ + 1
                jEnd = nJ + 3
             elseif(jDir == -1) then
                jBegin = -2
                jEnd = 0
             endif


             k = 1
             do j = jBegin, jEnd
                do iVar = 1, nVar
                   CellValue_I(1) = Field_VG(iVar,i-2,j,k)
                   CellValue_I(2) = Field_VG(iVar,i-1,j,k)
                   CellValue_I(3) = Field_VG(iVar,i+1,j,k)
                   CellValue_I(4) = Field_VG(iVar,i+2,j,k)
                   Orig = &
                        c3*(CellValue_I(2) + CellValue_I(3)) + &
                        c2*(CellValue_I(1) + CellValue_I(4)) + &
                        c1*(Field_VG(iVar,i+3,j,k) + &
                        Field_VG(iVar,i-3,j,k)) 

                   Field_VG(iVar,i,j,k) = limit_interpolation(Orig,&
                        CellValue_I,Distance_I)                                      

                enddo
             enddo
          enddo ! iDir1
       endif
    enddo; enddo
  end subroutine correct_ghost_for_coarse_blk

end module BATL_high_order
