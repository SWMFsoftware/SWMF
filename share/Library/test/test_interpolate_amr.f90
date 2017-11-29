!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS============================================
module ModTestInterpolateAMR

  use ModInterpolateAMR, ONLY: &
       interpolate_amr, interpolate_amr_gc, get_reference_block
  use ModUtilities, ONLY: CON_stop

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
  integer,parameter:: Out_  = -100
  integer,parameter:: nCell = 2
  integer,parameter:: nG    = 1
  !\
  ! Ratio of lengths in different dimensions
  !/
  real, parameter:: SizeRatio_D(3) = (/1.0, 100.0, 0.01/)
  !\
  ! which coordinates are periodic
  !/
  logical:: IsPeriodic_D(3) = .false.
  
contains
  !==================================================================
  subroutine test_interpolate_amr(&
       nDim,IsPeriodicIn_D,&
       nSample, UseGeneric, UseGhostCell)
    use ModRandomNumber, ONLY: random_real

    integer, intent(in)::nDim, nSample
    logical, intent(in)::IsPeriodicIn_D(nDim)
    logical, intent(in)::UseGeneric
    logical, optional, intent(in):: UseGhostCell

    integer :: iIndexes_II(0:nDim+1,2**nDim)
    logical :: IsSecondOrder, IsPossible, IsOut
    real, dimension(nDim):: DxyzDomain_D, DxyzCoarseBlock_D, &
         DxyzFineBlock_D, DxyzCoarse_D, &
         DxyzFine_D, Xyz_D, XyzPass_D,  &
         XyzCont_D,                     &
         XyzInterpolated_D, XyzCorner_D,&
         XyzModulo_D,                   &
         Dxyz_D
    real    ::VarInterpolated, VarContInterpolated
    real, allocatable::Xyz0_DGB(:,:,:,:,:)    
    real, allocatable::Xyz_DGB(:,:,:,:,:)    
    ! index of neigboring block, last index is resolution level of neighbor
    integer, allocatable::DiLevelNei_IIIB(:,:,:,:)
    real, allocatable, dimension(:,:,:,:) :: Var_GB
    real    :: Weight_I(2**nDim)

    ! whether to test approximation
    logical:: DoTestApproximation
    !Loop variables
    integer :: iCase, iSample, iGrid, iSubGrid, i, j, k, iBlock, iDir
    integer :: iProc, iBlockNei
    integer :: nCell_D(3)  ! Cells per block
    integer :: iCellIndex_D(3)
    integer :: nIndexes
    integer:: iMisc , nGridOut

    integer:: iSeed = 1
    !-------------------------------------------------------------------
    nCell_D = 1; nCell_D(1:nDim) = nCell
    DxyzDomain_D        = 2*nCell*SizeRatio_D(1:nDim)
    DxyzCoarseBlock_D   = nCell*SizeRatio_D(1:nDim)
    DxyzFineBlock_D     = 0.5*nCell*SizeRatio_D(1:nDim)
    DxyzCoarse_D        = 1*SizeRatio_D(1:nDim)
    DxyzFine_D          = 0.5*SizeRatio_D(1:nDim)
    IsPeriodic_D(1:nDim)= IsPeriodicIn_D
    if(present(UseGhostCell))then
       DoTestApproximation = .not.any(IsPeriodicIn_D) .or. UseGhostCell
    else
       DoTestApproximation = .not.any(IsPeriodicIn_D) .or. .not.UseGeneric
    end if
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
                  UseGhostCell=UseGhostCell)
          else
             call find_test(nDim, Xyz_D, &
                  iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
             Xyz_D = XyzCorner_D + Xyz_D
             call check_interpolate_test(nDim, Xyz_D, iBlock, &
                  iProc, iBlockNei)
             if(iBlockNei /= iBlock) then
                iBlock = iBlockNei
                XyzCorner_D = Xyz_DGB(:,1,1,1,iBlock) - 0.5*DxyzFine_D
                Dxyz_D = DxyzFine_D
             end if
             call fix_coord(nDim, Xyz_D, XyzCorner_D, XyzCorner_D+nCell*Dxyz_D, XyzPass_D)
             call interpolate_amr_gc(&
                  nDim         = nDim, &
                  Xyz_D        = XyzPass_D, &
                  XyzMin_D     = XyzCorner_D, &
                  DXyz_D       = DXyz_D, &
                  nCell_D      = nCell_D, &
                  DiLevelNei_III = DiLevelNei_IIIB(:,:,:,iBlock), &
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
          ! for periodic boundary conditions apply modulo
          XyzModulo_D = XyzInterpolated_D
          where(IsPeriodic_D(1:nDim))
             XyzModulo_D = modulo(XyzModulo_D, DxyzDomain_D)
          end where
          if(  DoTestApproximation.and.&
               any(abs(Xyz_D - XyzModulo_D) > 1.0e-6).and.&
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
             XyzCont_D(iDir) = Xyz_D(iDir) + &
                  SizeRatio_D(iDir)*(0.02*random_real(iSeed) - 0.01)
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
                  UseGhostCell=UseGhostCell)
          else
             call find_test(nDim, XyzCont_D, &
                  iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
             XyzCont_D = XyzCorner_D + XyzCont_D
             call check_interpolate_test(nDim, XyzCont_D, iBlock, &
                  iProc, iBlockNei)
             if(iBlockNei /= iBlock) then
                iBlock = iBlockNei
                XyzCorner_D = Xyz_DGB(:,1,1,1,iBlock) - 0.5*DxyzFine_D
                Dxyz_D = DxyzFine_D
             end if
             call fix_coord(nDim, XyzCont_D, XyzCorner_D, XyzCorner_D+nCell*Dxyz_D, XyzPass_D)
             call interpolate_amr_gc(&
                  nDim         = nDim, &
                  Xyz_D        = XyzPass_D, &
                  XyzMin_D     = XyzCorner_D, &
                  DXyz_D       = DXyz_D, &
                  nCell_D      = nCell_D, &
                  DiLevelNei_III = DiLevelNei_IIIB(:,:,:,iBlock), &
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
          XyzModulo_D = XyzInterpolated_D
          where(IsPeriodic_D(1:nDim))
             XyzModulo_D = modulo(XyzModulo_D, DxyzDomain_D)
          end where
          if(  DoTestApproximation.and.&
               any(abs(XyzCont_D - XyzModulo_D) > 1.0e-6).and.&
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
                  IsSecondOrder=IsSecondOrder)
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
    subroutine fix_coord(nDim, XyzIn_D, XyzBlockMin_D, XyzBlockMax_D, XyzOut_D)
      ! for periodic/flipped coordinate may need to adjust point's coordinates
      ! for calling interpolate_amr_gc
      integer,intent(in) :: nDim
      real,   intent(in) :: XyzIn_D(nDim)
      real,   intent(in) :: XyzBlockMin_D(nDim)
      real,   intent(in) :: XyzBlockMax_D(nDim)
      real,   intent(out):: XyzOut_D(nDim)
      
      real   :: Dxyz
      integer:: iDim
      !------------------------------------------------------------------------
      XyzOut_D = XyzIn_D
      do iDim = 1, nDim
         if(  XyzOut_D(iDim) < XyzBlockMax_D(iDim) .and.&
              XyzOut_D(iDim) >=XyzBlockMin_D(iDim))&
              CYCLE
         Dxyz = (XyzBlockMax_D(iDim) - XyzBlockMin_D(iDim)) / nCell_D(iDim)
         if(IsPeriodic_D(iDim))then
            if(XyzBlockMax_D(iDim)==DxyzDomain_D(iDim).and.&
                 XyzOut_D(iDim) < Dxyz)then
               XyzOut_D(iDim) = XyzOut_D(iDim) + DxyzDomain_D(iDim)
            elseif(XyzBlockMin_D(iDim)==0.0.and.&
                 XyzOut_D(iDim) >= DxyzDomain_D(iDim)-Dxyz)then
               XyzOut_D(iDim) = XyzOut_D(iDim) - DxyzDomain_D(iDim)
            end if
         end if
      end do
    end subroutine fix_coord
    !============================
    subroutine check_interpolate_test(nDim, Xyz_D, iBlockIn, &
         iProcOut, iBlockOut)
      integer, intent(in) :: nDim
      real,    intent(in) :: Xyz_D(nDim)
      integer, intent(in) :: iBlockIn
      integer, intent(out):: iProcOut
      integer, intent(out):: iBlockOut

      integer:: iLevel_I(2**nDim)
      logical:: IsOut_I(2**nDim)
      integer:: iDiscr_D(3)
      real   :: XyzCentral_D(nDim)
      real   :: XyzNei_D(nDim), XyzCorner_D(nDim), Dxyz_D(nDim)
      real   :: XyzGrid_DI(nDim, 2**nDim)
      logical:: IsOut
      integer:: iGridRef, iGrid
      integer:: iShiftRef_D(nDim)
      !------------------------------------------------------------
      ! determine displacement of the point relative to the block's interior
      iDiscr_D = 0
      where(Xyz_D < Xyz_DGB(:,1,1,1,iBlockIn))
         iDiscr_D(1:nDim) = -1
      elsewhere(Xyz_D >= Xyz_DGB(:,nCell_D(1),nCell_D(2),nCell_D(3),iBlockIn))
         iDiscr_D(1:nDim) =  1
      end where

      ! array of levels of neighbors
      iLevel_I = reshape(DiLevelNei_IIIB(&
           (/MIN(iDiscr_D(1),0), MAX(0,iDiscr_D(1))/), &
           (/MIN(iDiscr_D(2),0), MAX(0,iDiscr_D(2))/), &
           (/MIN(iDiscr_D(3),0), MAX(0,iDiscr_D(3))/), iBlockIn), (/2**nDim/))

      ! find those that are outside of the computational domain
      IsOut_I = iLevel_I < -1
      
      ! iLevel_I may contain -1's;
      ! fix so there are only 0's (Coarser) and 1's (Finer)
      if(any(iLevel_I == -1 .and. .not. IsOut_I))&
           iLevel_I = iLevel_I + 1
      
      ! cell size of this block
      if(iBlockIn > 2**nDim)then
         Dxyz_D = DxyzFine_D   ! the block is Finer 
      else
         Dxyz_D = DxyzCoarse_D ! the block is Coarse
      end if

      ! coordinates of block's junction
      where(    iDiscr_D(1:nDim) == 1)
         XyzCentral_D = Xyz_DGB(:,nCell,nCell,nCell,iBlockIn)  + 0.5 * Dxyz_D
      elsewhere(iDiscr_D(1:nDim) ==-1)
         XyzCentral_D = Xyz_DGB(:,1,1,1,iBlockIn)  - 0.5 * Dxyz_D
      elsewhere
         XyzCentral_D = Xyz_D
      end where

      ! supergrid
      do iGrid = 1, 2**nDim
         XyzGrid_DI(:, iGrid) = XyzCentral_D + &
              (iShift_DI(1:nDim, iGrid)-0.5) * DxyzCoarse_D
      end do
     
      ! find the block that has enough information to perform interpolation
      ! using only 1 layer of gc
      call get_reference_block(&
           nDim, Xyz_D, XyzGrid_DI, iLevel_I, IsOut_I, iGridRef)
      ! displacement to this block
      iShiftRef_D = iShift_DI(1:nDim,iGridRef)

      ! check if it is the input block
      if(all((iDiscr_D(1:nDim)-1)/2 == iShiftRef_D*ABS(iDiscr_D(1:nDim))))then
         iBlockOut = iBlockIn
         RETURN
      end if

      ! it is a different block, find it
      XyzNei_D = XyzCentral_D + 0.5*DxyzFine_D*(iShiftRef_D-0.5)
      call find_test(nDim, XyzNei_D, &
           iProcOut, iBlockOut, XyzCorner_D, Dxyz_D, IsOut)

    end subroutine check_interpolate_test
    !===========================s
    subroutine get_level(Xyz_D, iLevel_I, iLevelOut) 
      ! gets refinement level at location Xyz_D
      real,    intent(in ):: Xyz_D(nDim)
      integer, intent(in ):: iLevel_I(2**nDim)
      integer, intent(out):: iLevelOut

      real:: XyzModulo_D(nDim)
      integer:: iDim, iGrid
      integer, parameter:: nTwoPower_I(3) = (/1,2,4/)
      !-----------------------------------------------
      XyzModulo_D = Xyz_D
      where(IsPeriodic_D(1:nDim))
         XyzModulo_D = modulo(Xyz_D, DxyzDomain_D(1:nDim))
      end where
      if(any(XyzModulo_D < 0).or.any(XyzModulo_D >= DxyzDomain_D(1:nDim)))then
         iLevelOut = Out_
         RETURN
      end if
      iGrid = 1 + SUM(nTwoPower_I(1:nDim), &
           MASK = XyzModulo_D >= 0.5 * DxyzDomain_D(1:nDim))
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
      real, dimension(nDim):: Xyz_D, XyzCorner_D, Dxyz_D, DxyzModulo_D
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
                  DxyzModulo_D = Xyz_DGB(:,i,j,k,iBlock) - Xyz_D - XyzCorner_D
                  where(    DxyzModulo_D > Dxyz_D)
                     DxyzModulo_D = DxyzDomain_D
                  elsewhere(DxyzModulo_D <-Dxyz_D)
                     DxyzModulo_D =-DxyzDomain_D
                  elsewhere
                     DxyzModulo_D = 0.0
                  end where
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
                  where(IsPeriodic_D(1:nDim))
                     Xyz_DGB(:,i,j,k,iBlock) = &
                          Xyz_DGB(:,i,j,k,iBlock) + DxyzModulo_D
                  end where
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
    DxyzDomain_D      = 2*nCell*SizeRatio_D(1:nDim)
    DxyzCoarseBlock_D = nCell*SizeRatio_D(1:nDim)
    DxyzFineBlock_D   = 0.5*nCell*SizeRatio_D(1:nDim)
    DxyzCoarse_D      = 1*SizeRatio_D(1:nDim)
    DxyzFine_D        = 0.5*SizeRatio_D(1:nDim)
    iProc = 0; iBlock=0; XyzCorner_D=0.0; Dxyz_D = 0.0
    ! fix periodic coordinates
    where(IsPeriodic_D(1:nDim))
       Xyz_D = modulo(Xyz_D, DxyzDomain_D)
    end where
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

  call test(&
       nDim          = 2,&
       IsPeriodicIn_D= (/.true.,.false./), &
       nSample       = 20000, &
       UseGeneric    = .true., &
       UseGhostCell  = .false.)

  call test(&
       nDim          = 2,&
       IsPeriodicIn_D= (/.false.,.true./), &
       nSample       = 20000, &
       UseGeneric    = .true., &
       UseGhostCell  = .true.)

  call test(&
       nDim          = 2,&
       IsPeriodicIn_D= (/.false.,.true./), &
       nSample       = 20000, &
       UseGeneric    = .false.)


  call test(&
       nDim          = 3,&
       IsPeriodicIn_D= (/.false.,.true.,.true./), &
       nSample       = 20000, &
       UseGeneric    = .true., &
       UseGhostCell  = .false.)

  call test(&
       nDim          = 3,&
       IsPeriodicIn_D= (/.false.,.false.,.true./), &
       nSample       = 20000, &
       UseGeneric    = .true., &
       UseGhostCell  = .true.)

  call test(&
       nDim          = 3,&
       IsPeriodicIn_D= (/.true.,.false.,.true./), &
       nSample       = 20000, &
       UseGeneric    = .false.)

end program test_interpolate_amr


