!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS============================================
module ModTestInterpolateAMR
  use ModInterpolateAMR, ONLY: interpolate_amr, interpolate_amr_gc
  implicit none
  !\
  ! Shift of the iGrid point in the stencil with respect to the
  ! first one
  !/
  integer, dimension(3,8),parameter :: iShift_DI = reshape((/&
       0, 0, 0, &
       1, 0, 0, &
       0, 1, 0, &
       1, 1, 0, &
       0, 0, 1, &
       1, 0, 1, &
       0, 1, 1, &
       1, 1, 1/),(/3,8/))

  !\
  ! For test: arrays of the refinement levels to be passed to
  ! find_test routine
  !/
  integer:: iLevelTest_I(8)
  integer,parameter:: nCell = 2
  integer,parameter:: nG    = 2
contains
  !==================================================================
  subroutine test_interpolate_amr(nDim,nSample, UseGeneric)
    use ModRandomNumber, ONLY: random_real

    integer, intent(in)::nDim, nSample
    logical, intent(in)::UseGeneric

    integer :: iIndexes_II(0:nDim+1,2**nDim)
    logical :: IsSecondOrder, IsPossible, IsOut
    real, dimension(nDim):: DxyzDomain_D, DxyzCoarseBlock_D, &
         DxyzFineBlock_D, DxyzCoarse_D, &
         DxyzFine_D, Xyz_D,             &
         XyzCont_D,                     &
         XyzInterpolated_D, XyzCorner_D,&
         Dxyz_D
    real    ::VarInterpolated, VarContInterpolated
    real, allocatable::Xyz0_DGB(:,:,:,:,:)    
    real, allocatable::Xyz_DGB(:,:,:,:,:)    
    ! index of neigboring block, last index is resolution level of neighbor
    integer, allocatable::DiLevelNei_IIIB(:,:,:,:)
    integer :: DiLevelNei_I(3**nDim)
    real, allocatable, dimension(:,:,:,:) :: Var_GB
    real    :: Weight_I(2**nDim)
    !Loop variables
    integer :: iCase, iSample, iGrid, iSubGrid, i, j, k, iBlock, iDir
    integer :: iProc, iBlockNei, iBoundary
    integer :: nCell_D(3)  ! Cells per block
    integer :: iCellIndex_D(3)
    integer :: nIndexes
    integer:: iMisc , nGridOut

    integer:: iSeed = 1
    !-------------------------------------------------------------------
    nCell_D = 1; nCell_D(1:nDim) = nCell
    DxyzDomain_D      = 2*nCell
    DxyzCoarseBlock_D = nCell
    DxyzFineBlock_D   = 0.5*nCell
    DxyzCoarse_D      = 1
    DxyzFine_D        = 0.5
    allocate(DiLevelNei_IIIB(-1:1,-1:1,-1:1,(2**nDim)*(2**nDim+1)))
    DiLevelNei_IIIB = 0
    allocate(Xyz0_DGB(nDim, 1-nG:nCell_D(1)+nG, &
         1-nG:nCell_D(2)+nG, &
         1-nG:nCell_D(3)+nG, (2**nDim)*(2**nDim+1)))
    Xyz0_DGB = 0
    allocate(Xyz_DGB(nDim, 1-nG:nCell_D(1)+nG, &
         1-nG:nCell_D(2)+nG, &
         1-nG:nCell_D(3)+nG, (2**nDim)*(2**nDim+1)))
    Xyz_DGB = 0
    allocate(Var_GB(1-nG:nCell_D(1)+nG, &
         1-nG:nCell_D(2)+nG, &
         1-nG:nCell_D(3)+nG, (2**nDim)*(2**nDim+1)))
    Var_GB = 0
    do iGrid = 1, 2**nDim
       iBlock = iGrid
       XyzCorner_D = DxyzCoarseBlock_D*iShift_DI(1:nDim,iGrid)
       do k = 1-nG, nCell_D(3)+nG
          do j = 1-nG, nCell_D(2)+nG
             do i = 1-nG, nCell_D(1)+nG
                iCellIndex_D = (/i,j,k/)
                Xyz_DGB(:,i,j,k,iBlock) = XyzCorner_D +&
                     DxyzCoarse_D*(iCellIndex_D(1:nDim) - 0.50)
                ! fill values at physical cells first
                if(all(iCellIndex_D >= 1).and.all(iCellIndex_D <= nCell_D))&
                     Var_GB(i,j,k,iBlock) = random_real(iSeed)                 
             end do
          end do
       end do
       do iSubGrid = 1, 2**nDim
          iBlock = iGrid*(2**nDim)+iSubGrid
          XyzCorner_D = DxyzCoarseBlock_D*iShift_DI(1:nDim,iGrid) + &
               DxyzFineBlock_D*iShift_DI(1:nDim,iSubGrid)
          do k = 1-nG, nCell_D(3)+nG
             do j = 1-nG, nCell_D(2)+nG
                do i = 1-nG, nCell_D(1)+nG
                   iCellIndex_D = (/i,j,k/)
                   Xyz_DGB(:,i,j,k,iBlock) = XyzCorner_D +&
                        DxyzFine_D*(iCellIndex_D(1:nDim) - 0.50)
                   ! fill values at physical cells first
                   if(all(iCellIndex_D >= 1).and.all(iCellIndex_D <= nCell_D))&
                        Var_GB(i,j,k,iBlock) = 0.25 + 0.50 * random_real(iSeed)
                end do
             end do
          end do
       end do
    end do
    Xyz0_DGB = Xyz_DGB

    nIndexes = nDim +1
    CASE:do iCase = 0, 2**(2**nDim) - 2
       iLevelTest_I = 0; iGrid = 0
       iMisc = iCase
       do while(iMisc > 0)
          iGrid = iGrid + 1
          iLevelTest_I(iGrid) = mod(iMisc, 2)
          iMisc = (iMisc - iLevelTest_I(iGrid))/2
       end do
       !write(*,*)'Case=',iLevelTest_I(1:2**nDim)
       call fill_ghost_cells(iLevelTest_I)
       call fill_nei_levels(iLevelTest_I)
       !\
       ! We generated refinement, now sample points
       !/
       SAMPLE:do iSample = 1, nSample
          do iDir = 1, nDim
             Xyz_D(iDir) = (0.01 +0.98*random_real(iSeed))*DxyzDomain_D(iDir)
          end do
          !\
          ! call interpolate_amr
          !/
          if(UseGeneric)then
             call interpolate_amr(&
                  nDim=nDim, &
                  XyzIn_D=Xyz_D, &
                  nIndexes=nDim+1,&
                  find=find_test, &
                  nCell_D=nCell_D(1:nDim),&
                  nGridOut=nGridOut,&
                  Weight_I=Weight_I,&
                  iIndexes_II=iIndexes_II,&
                  IsSecondOrder=IsSecondOrder,&
                  UseGhostCell=.true.)
          else
             call find_test(nDim, Xyz_D, &
                  iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
             Xyz_D = XyzCorner_D + Xyz_D
             call check_interpolate_test(nDim, Xyz_D, iBlock, &
                  IsPossible, iProc, iBlockNei, iBoundary)
             if(.not. IsPossible) then
                iBlock = iBlockNei
                XyzCorner_D = Xyz_DGB(:,1,1,1,iBlock) - 0.5*DxyzFine_D
                Dxyz_D = DxyzFine_D
             end if
             if(nDim==2)then
                DiLevelNei_I=reshape(DiLevelNei_IIIB(:,:,0,iBlock),(/3**nDim/))
             else
                DiLevelNei_I=reshape(DiLevelNei_IIIB(:,:,:,iBlock),(/3**nDim/))
             end if

             call interpolate_amr_gc(&
                  nDim         = nDim, &
                  Xyz_D        = Xyz_D, &
                  XyzMin_D     = XyzCorner_D, &
                  DXyz_D       = DXyz_D, &
                  nCell_D      = nCell_D, &
                  DiLevelNei_I = DiLevelNei_I, &
                  nCellOut     = nGridOut, &
                  iCellOut_II  = iIndexes_II(1:nDim,:), &
                  Weight_i     = Weight_I, &
                  IsSecondOrder= IsSecondOrder)
             iIndexes_II(nIndexes,:) = iBlock
          end if
          !\          
          !Compare with interpolated:
          !/
          XyzInterpolated_D = 0
          VarInterpolated   = 0
          do iGrid = 1, nGridOut
             iCellIndex_D = 1
             iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
             iBlock = iIndexes_II(nIndexes,iGrid)
             XyzInterpolated_D = XyzInterpolated_D + Weight_I(iGrid)*&
                  Xyz_DGB(:,iCellIndex_D(1), iCellIndex_D(2), &
                  iCellIndex_D(3), iBlock)
             VarInterpolated = VarInterpolated + &
                  Weight_I(iGrid)*&
                  Var_GB(iCellIndex_D(1), iCellIndex_D(2), &
                  iCellIndex_D(3), iBlock)
          end do
          if(any(abs(Xyz_D - XyzInterpolated_D) > 1.0e-6).and.&
               IsSecondOrder)then
             write(*,*)'Approximation test failed'
             write(*,*)'Grid:', iLevelTest_I(1:2**nDim)
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'Point=', Xyz_D
             write(*,*)'Cell_D  iBlock XyzGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCellIndex_D = 1
                iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Xyz_DGB(:,iCellIndex_D(1), iCellIndex_D(2), &
                     iCellIndex_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'Xyz_D=',Xyz_D
             write(*,*)'XyzInterpolated_D=',XyzInterpolated_D
             call CON_stop('Correct code and redo test')
          end if
          !\
          ! Test continuity
          !/
          do iDir =1, nDim
             XyzCont_D(iDir) = Xyz_D(iDir) + (0.02*random_real(iSeed) - 0.01)
          end do
          !\
          ! call interpolate_amr
          !/
          if(UseGeneric)then
             call interpolate_amr(&
                  nDim=nDim, &
                  XyzIn_D=XyzCont_D, &
                  nIndexes=nDim+1,&
                  find=find_test, &
                  nCell_D=nCell_D(1:nDim),&
                  nGridOut=nGridOut,&
                  Weight_I=Weight_I,&
                  iIndexes_II=iIndexes_II,&
                  IsSecondOrder=IsSecondOrder,&
                  UseGhostCell=.true.)
          else
             call find_test(nDim, XyzCont_D, &
                  iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
             XyzCont_D = XyzCorner_D + XyzCont_D
             call check_interpolate_test(nDim, XyzCont_D, iBlock, &
                  IsPossible, iProc, iBlockNei, iBoundary)
             if(.not. IsPossible) then
                iBlock = iBlockNei
                XyzCorner_D = Xyz_DGB(:,1,1,1,iBlock) - 0.5*DxyzFine_D
                Dxyz_D = DxyzFine_D
             end if
             if(nDim==2)then
                DiLevelNei_I=reshape(DiLevelNei_IIIB(:,:,0,iBlock),(/3**nDim/))
             else
                DiLevelNei_I=reshape(DiLevelNei_IIIB(:,:,:,iBlock),(/3**nDim/))
             end if
             call interpolate_amr_gc(&
                  nDim         = nDim, &
                  Xyz_D        = XyzCont_D, &
                  XyzMin_D     = XyzCorner_D, &
                  DXyz_D       = DXyz_D, &
                  nCell_D      = nCell_D, &
                  DiLevelNei_I = DiLevelNei_I, &
                  nCellOut     = nGridOut, &
                  iCellOut_II  = iIndexes_II(1:nDim,:), &
                  Weight_i     = Weight_I, &
                  IsSecondOrder= IsSecondOrder)
             iIndexes_II(nIndexes,:) = iBlock
          end if
          !\          
          !Compare interpolated values of Var:
          !/
          VarContInterpolated = 0
          XyzInterpolated_D = 0
          do iGrid = 1, nGridOut
             iCellIndex_D = 1
             iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
             iBlock = iIndexes_II(nIndexes,iGrid)
             VarContInterpolated = VarContInterpolated + &
                  Weight_I(iGrid)*&
                  Var_GB(iCellIndex_D(1), iCellIndex_D(2), &
                  iCellIndex_D(3), iBlock)
             XyzInterpolated_D = XyzInterpolated_D + Weight_I(iGrid)*&
                  Xyz_DGB(:,iCellIndex_D(1), iCellIndex_D(2), &
                  iCellIndex_D(3), iBlock)
          end do
          if(any(abs(XyzCont_D - XyzInterpolated_D) > 1.0e-6).and.&
               IsSecondOrder)then
             write(*,*)'Approximation test failed'
             write(*,*)'Grid:', iLevelTest_I(1:2**nDim)
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'PointCont=', XyzCont_D
             write(*,*)'Cell_D  iBlock XyzGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCellIndex_D = 1
                iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Xyz_DGB(:,iCellIndex_D(1), iCellIndex_D(2), &
                     iCellIndex_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'XyzCont_D=',XyzCont_D
             write(*,*)'XyzInterpolated_D=',XyzInterpolated_D
             call CON_stop('Correct code and redo test')
          end if
          if(abs(VarContInterpolated - VarInterpolated) > nDim*0.01)then
             write(*,*)'Continuity test failed'
             write(*,*)'Grid:', iLevelTest_I
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'XyzCont=', XyzCont_D
             write(*,*)'Cell_D  iBlock XyzGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCellIndex_D = 1
                iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Xyz_DGB(:,iCellIndex_D(1), iCellIndex_D(2), &
                     iCellIndex_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'Xyz_D=',Xyz_D
             call interpolate_amr(&
                  nDim=nDim, &
                  XyzIn_D=Xyz_D, &
                  nIndexes=nDim+1,&
                  find=find_test, &
                  nCell_D=nCell_D(1:nDim),&
                  nGridOut=nGridOut,&
                  Weight_I=Weight_I,&
                  iIndexes_II=iIndexes_II,&
                  IsSecondOrder=IsSecondOrder,&
                  UseGhostCell=.true.)
             write(*,*)'Cell_D  iBlock XyzGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCellIndex_D = 1
                iCellIndex_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Xyz_DGB(:,iCellIndex_D(1), iCellIndex_D(2), &
                     iCellIndex_D(3), iBlock), Weight_I(iGrid)
             end do
             call CON_stop('Correct code and redo test')
          end if
       end do SAMPLE
    end do CASE
    deallocate(Xyz_DGB, Var_GB, DiLevelNei_IIIB)
  contains
    !============================
    subroutine check_interpolate_test(nDim, Xyz_D, iBlockIn, &
         IsPossible, iProcOut, iBlockOut, iBoundary)
      integer, intent(in) :: nDim
      real,    intent(in) :: Xyz_D(nDim)
      integer, intent(in) :: iBlockIn
      logical, intent(out):: IsPossible
      integer, intent(out):: iProcOut
      integer, intent(out):: iBlockOut
      integer, intent(out):: iBoundary

      integer:: iDiscr_D(3)
      integer:: iShift_D(3)
      integer:: iBlock
      real   :: XyzNei_D(nDim), XyzCorner_D(nDim), Dxyz_D(nDim)
      logical:: IsOut
      !------------------------------------------------------------
      IsPossible = .true.
      if(iBlockIn > 2**nDim) RETURN ! block is Fine

      iDiscr_D = 0
      where(Xyz_D < Xyz_DGB(:,1,1,1,iBlockIn))
         iDiscr_D(1:nDim) = -1
      elsewhere(Xyz_D >= Xyz_DGB(:,nCell_D(1),nCell_D(2),nCell_D(3),iBlockIn))
         iDiscr_D(1:nDim) =  1
      end where

      if(all(iDiscr_D==0)) RETURN ! within coarse block and far from boundary

      if(any(DiLevelNei_IIIB(&
           MIN(iDiscr_D(1),0):MAX(0,iDiscr_D(1)), &
           MIN(iDiscr_D(2),0):MAX(0,iDiscr_D(2)), &
           MIN(iDiscr_D(3),0):MAX(0,iDiscr_D(3)), &
           iBlockIn) > 0) )then
         IsPossible = .false.
         do iBlock = 1, 2**nDim
            iShift_D = iShift_DI(:,iBlock) - iShift_DI(:,iBlockIn)
            if(iLevelTest_I(iBlock)==1.and. &
                 all(iShift_D==iDiscr_D.or.iShift_D==0)) EXIT
         end do
         XyzNei_D = Xyz_D + DxyzFine_D*iShift_D(1:nDim)
         call find_test(nDim, XyzNei_D, &
            iProcOut, iBlockOut, XyzCorner_D, Dxyz_D, IsOut)
         iBoundary = sum(iDiscr_D*(/1, 3, 9/))-3**nDim/2
      end if
      
    end subroutine check_interpolate_test
    !===========================
    subroutine get_level(Xyz_D, iLevel_I, iLevelOut) 
      ! gets refinement level at location Xyz_D
      real,    intent(in ):: Xyz_D(nDim)
      integer, intent(in ):: iLevel_I(2**nDim)
      integer, intent(out):: iLevelOut

      integer:: iDim, iGrid
      integer, parameter:: nTwoPower_I(3) = (/1,2,4/)
      !-----------------------------------------------
      if(any(Xyz_D < 0) .or. any(Xyz_D >= DxyzDomain_D(1:nDim)))then
         iLevelOut = -100
         RETURN
      end if
      iGrid = 1 + SUM(nTwoPower_I(1:nDim), &
           MASK = Xyz_D >= 0.5 * DxyzDomain_D(1:nDim))
      iLevelOut = iLevel_I(iGrid)
    end subroutine get_level
    !===========================
    subroutine fill_nei_levels(iLevel_I)
      ! fill DiLevelNei_IIB that store levels of neighbors
      integer, intent(in):: iLevel_I(2**nDim)

      integer:: iGrid, iGridNei, iSubGrid, iSubGridNei
      integer:: iLevelNei
      integer:: iShift_D(3), iShiftGhost_DI(3,2)
      real:: XyzBlockCenterCoarse_D(nDim)
      real:: XyzBlockCenterFine_D(nDim)
      real:: Xyz_D(nDim)
      integer:: i,j,k, iIndex_I(3)
      !-----------------------------------------------
      DiLevelNei_IIIB = 0
      do iGrid = 1, 2**nDim
         XyzBlockCenterCoarse_D = &
              (0.25 + 0.5*iShift_DI(1:nDim, iGrid))*DxyzDomain_D
         if(iLevel_I(iGrid) == 0)then
            ! iGrid is Coarse
            do i = -1, 1; do j = -1, 1; do k = -(nDim-2), (nDim-2)
               iIndex_I = (/i,j,k/)
               Xyz_D = XyzBlockCenterCoarse_D + &
                    0.6 * iIndex_I(1:nDim) * DxyzCoarseBlock_D
               call get_level(Xyz_D(1:nDim), iLevel_I, iLevelNei)
               DiLevelNei_IIIB(i,j,k,iGrid) = iLevelNei
            end do; end do; end do
         else
            ! iGrid is Fine
            do iSubGrid = 1, 2**nDim
               XyzBlockCenterFine_D = XyzBlockCenterCoarse_D + &
                    (2*iShift_DI(1:nDim,iSubGrid) - 1) * 0.5 * DxyzFineBlock_D
               do i = -1, 1; do j = -1, 1; do k = -(nDim-2), (nDim-2)
                  iIndex_I = (/i,j,k/)
                  Xyz_D = XyzBlockCenterFine_D + &
                       0.6 * iIndex_I(1:nDim) * DxyzFineBlock_D
                  call get_level(Xyz_D(1:nDim), iLevel_I, iLevelNei)
                  DiLevelNei_IIIB(i,j,k,iGrid*(2**nDim)+iSubGrid) = iLevelNei-1
               end do; end do; end do
            end do
         end if
      end do
    end subroutine fill_nei_levels
    !============================
    subroutine fill_ghost_cells(iLevel_I)
      ! values at ghost cells depend on the refinement
      integer, intent(in):: iLevel_I(2**nDim)

      logical:: Use_B(2**(2*nDim)+2**nDim), IsOut
      integer:: iDim, iProc, iBlock, iBlockNei, i, j ,k
      integer:: iGrid, iSubGrid
      integer:: iCellIndex_D(3), iCellIndexNei_D(3)
      real, dimension(nDim):: Xyz_D, XyzCorner_D, Dxyz_D
      !-----------------------------------------------
      Xyz_DGB = Xyz0_DGB
      ! find blocks that are used in the current refinement
      Use_B = .false.
      do iGrid = 1, 2**nDim
         if(iLevel_I(iGrid)==0)then
            Use_B(iGrid) = .true.
         else
            do iSubGrid = 1, 2**nDim
               Use_B(iGrid*(2**nDim)+iSubGrid) = .true.
            end do
         end if
      end do
      ! fill values
      do iBlock = 1, 2**(2*nDim)+2**nDim
         if(.not. Use_B(iBlock))CYCLE
         do k = 1-nG, nCell_D(3)+nG 
            do j = 1-nG, nCell_D(2)+nG 
               do i = 1-nG, nCell_D(1)+nG
                  iCellIndex_D = (/i,j,k/)
                  if(all(iCellIndex_D >= 1).and.all(iCellIndex_D <= nCell_D))&
                       CYCLE
                  Xyz_D = Xyz_DGB(:,i,j,k,iBlock)
                  call find_test(nDim, Xyz_D, &
                       iProc, iBlockNei, XyzCorner_D, Dxyz_D, IsOut)
                  if(IsOut) CYCLE
                  if(iBlock <= 2**nDim            &! iBlock    is Coarse
                       .and. iBlockNei > 2**nDim) &! iBlockNei is Fine
                       CYCLE

                  iCellIndexNei_D = 1
                  iCellIndexNei_D(1:nDim) = &
                       nint(0.3+Xyz_D(1:nDim)/Dxyz_D(1:nDim))
                  Var_GB(i,j,k,iBlock) = &
                       Var_GB(iCellIndexNei_D(1),&
                       iCellIndexNei_D(2),&
                       iCellIndexNei_D(3),iBlockNei)
                  Xyz_DGB(:,i,j,k,iBlock) = &
                       Xyz_DGB(:,iCellIndexNei_D(1),&
                       iCellIndexNei_D(2),&
                       iCellIndexNei_D(3),iBlockNei)
               end do
            end do
         end do
      end do
    end subroutine fill_ghost_cells
    !============================
    subroutine find_block(nDim, XyzIn_D, iBlock)
      integer, intent(in) :: nDim
      real,    intent(in) :: XyzIn_D(nDim)
      integer, intent(out):: iBlock 

      integer:: iProc
      logical:: IsOut
      real   :: Xyz_D(nDim), XyzCorner_D(nDim), Dxyz_D(nDim)
    !---------------------------------------------------------
    Xyz_D = XyzIn_D
    call find_test(nDim, Xyz_D, &
            iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
    end subroutine find_block
  end Subroutine test_interpolate_amr
  !============================
  subroutine find_test(nDim, Xyz_D, &
            iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
    integer, intent(in) :: nDim
    !\
    ! "In"- the coordinates of the point, "out" the coordinates of the
    ! point with respect to the block corner. In the most cases
    ! XyzOut_D = XyzIn_D - XyzCorner_D, the important distinction,
    ! however, is the periodic boundary, near which the jump in the
    ! stencil coordinates might occur. To handle the latter problem,
    ! we added the "out" intent. The coordinates for the stencil
    ! and input point are calculated and recalculated below with
    ! respect to the block corner. 
    !/
    real,  intent(inout):: Xyz_D(nDim)
    integer, intent(out):: iProc, iBlock !processor and block number
    !\
    ! Block left corner coordinates and the grid size:
    !/
    real,    intent(out):: XyzCorner_D(nDim), Dxyz_D(nDim)
    logical, intent(out):: IsOut !Point is out of the domain.
    real, dimension(nDim):: DxyzDomain_D, DxyzCoarseBlock_D ,&
         DxyzFineBlock_D, DxyzCoarse_D, DxyzFine_D
    integer:: iShift_D(3), iGrid, iSubGrid
    integer, dimension(0:1,0:1,0:1), parameter:: iGridFromShift_III=reshape(&
         (/1, 2, 3, 4, 5, 6, 7, 8/),(/2, 2, 2/))
    logical, dimension(nDim) :: IsAboveCenter_D
    !------------------- 
    DxyzDomain_D      = 2*nCell
    DxyzCoarseBlock_D = nCell
    DxyzFineBlock_D   = 0.5*nCell 
    DxyzCoarse_D      = 1
    DxyzFine_D        = 0.5
    iProc = 0; iBlock=0; XyzCorner_D=0.0; Dxyz_D = 0.0
    IsOut = any(Xyz_D < 0.0 .or. Xyz_D >= DxyzDomain_D)
    if(IsOut) RETURN
    !\
    ! Find into which coarse block the point fall
    !/ 
    IsAboveCenter_D = Xyz_D >= DxyzCoarseBlock_D
    iShift_D = 0
    where(IsAboveCenter_D)iShift_D(1:nDim) = 1
    XyzCorner_D = XyzCorner_D + DxyzCoarseBlock_D*iShift_D(1:nDim)
    Xyz_D       = Xyz_D       - DxyzCoarseBlock_D*iShift_D(1:nDim)
    iGrid = iGridFromShift_III(iShift_D(1),iShift_D(2),iShift_D(3))
    !\
    ! Check if the coarse block is used
    !/
    if(iLevelTest_I(iGrid)==0)then
       iBlock = iGrid
       Dxyz_D = DxyzCoarse_D
       RETURN
    end if
    !\
    ! The coarser block is refined, find into which fine block 
    ! the point falls
    !/
    IsAboveCenter_D = Xyz_D >= DxyzFineBlock_D
    iShift_D = 0
    where(IsAboveCenter_D)iShift_D(1:nDim) = 1
    XyzCorner_D = XyzCorner_D + DxyzFineBlock_D*iShift_D(1:nDim)
    Xyz_D       = Xyz_D       - DxyzFineBlock_D*iShift_D(1:nDim)
    iSubGrid = iGridFromShift_III(iShift_D(1),iShift_D(2),iShift_D(3))
    iBlock = iGrid*(2**nDim)+iSubGrid
    Dxyz_D = DxyzFine_D
  end subroutine find_test

end module ModTestInterpolateAMR

program test_interpolate_amr

  use ModTestInterpolateAMR, test => test_interpolate_amr

  implicit none

  call test(2,20000, .true.)
  call test(2,20000, .false.)
  call test(3,20000, .true.)
  call test(3,20000, .false.)

end program test_interpolate_amr

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  !----------------------------------------------------------------------------

  write(*,'(a)')StringError
  write(*,'(a)')'!!! SWMF_ABORT !!!'
  stop

end subroutine CON_stop

