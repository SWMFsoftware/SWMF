!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS============================================
module ModTestInterpolateAMR

  use ModInterpolateAMR, ONLY: &
       interpolate_amr, interpolate_amr_gc, get_reference_block
  use ModUtilities, ONLY: CON_stop
  use ModNumConst, ONLY: cPi
  
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

  integer, parameter::&
       R_ = 1, Phi_ = 2, Theta_ = 3

  integer, parameter:: Coarse_ = 0, Fine_ = 1

  !\
  ! For test: arrays of the refinement levels to be passed to
  ! find_test routine
  !/
  integer:: iLevelTest_I(16)
  integer,parameter:: Out_  = -100
  integer,parameter:: nCell = 4
  integer,parameter:: nG    = 1
  !\
  ! Ratio of lengths in different dimensions
  !/
  real, parameter:: SizeRatio_D(3) = (/1.0, 100.0, 0.01/)
  !\
  ! which coordinates are periodic
  !/
  logical:: IsPeriodic_D(3) = .false.

  integer :: nCell_D(3)  ! Cells per block

  character(len=*), parameter:: StringType = 'cartesian'

  real,    allocatable::  Coord_DB(:,:)
  real,    allocatable:: DCoord_DB(:,:)
  logical, allocatable:: Used_B(:)

  integer:: nRoot, nBlock
  logical:: IsSpherical

  ! Indices:
  ! 1-8: Children
  ! 9  : Level
  integer,allocatable:: iTree_IB(:,:)
  integer, parameter:: Level_ = 9

  real,allocatable:: CoordMin_D(:),CoordMax_D(:),DomainSize_D(:)
contains
  !==================================================================
  subroutine test_interpolate_amr(&
       nDim,IsCartesian,IsPeriodicIn_D,&
       nSample, UseGeneric, UseGhostCell)
    use ModRandomNumber, ONLY: random_real

    
    integer, intent(in)::nDim, nSample
    logical, intent(in)::IsPeriodicIn_D(nDim)
    logical, intent(in)::IsCartesian
    logical, intent(in)::UseGeneric
    logical, optional, intent(in):: UseGhostCell

    integer :: iIndexes_II(0:nDim+1,2**nDim)
    logical :: IsSecondOrder, IsPossible, IsOut
    real, dimension(nDim):: &
         CellSizeCoarse_D, CellSizeFine_D, &
         Coord_D, CoordPass_D,  &
         CoordCont_D,                     &
         CoordInterpolated_D, CoordCorner_D,&
         CoordModulo_D,                   &
         CellSize_D
    real    ::VarInterpolated, VarContInterpolated
    real, allocatable::Coord0_DGB(:,:,:,:,:)    
    real, allocatable::Coord_DGB(:,:,:,:,:)
    ! index of neigboring block, last index is resolution level of neighbor
    integer, allocatable::DiLevelNei_IIIB(:,:,:,:)
    real, allocatable, dimension(:,:,:,:) :: Var_GB
    real    :: Weight_I(2**nDim)


    ! whether to test approximation
    logical:: DoTestApproximation
    !Loop variables
    integer :: iCase, iSample, iGrid, iSubGrid, i, j, k, iBlock, iDir
    integer :: iProc, iBlockNei
    integer :: iCell_D(3)
    integer :: nIndexes
    integer:: iMisc , nGridOut

    integer:: iSeed = 1
    !-------------------------------------------------------------------
    nCell_D = 1; nCell_D(1:nDim) = nCell
    IsSpherical = nDim==3 .and. .not.IsCartesian
    if(IsCartesian)then
       IsPeriodic_D(1:nDim)= IsPeriodicIn_D
    else
       IsPeriodic_D(1:nDim) = .false.
       IsPeriodic_D(  Phi_) = .true.
    end if
    
    if(present(UseGhostCell))then
       DoTestApproximation = .not.any(IsPeriodicIn_D) .or. UseGhostCell
    else
       DoTestApproximation = .not.any(IsPeriodicIn_D) .or. .not.UseGeneric
    end if
    call generate_grid

    nIndexes = nDim +1
    CASE:do iCase = 0, 2**(2**nDim)*nRoot - 2 + (nRoot-1)
       iLevelTest_I = 0; iGrid = 0
       call set_levels(iCase)
       !write(*,*)            'Case=',iLevelTest_I(        1:2** nDim   )
       !if(nRoot>1) write(*,*)'     ',iLevelTest_I(2**nDim+1:2**(nDim+1))
       call fill_ghost_cells
       call fill_nei_levels
       !\
       ! We generated refinement, now sample points
       !/
       SAMPLE:do iSample = 1, nSample
          do iDir = 1, nDim
             Coord_D(iDir) = CoordMin_D(iDir) + &
                  (0.01 +0.98*random_real(iSeed))*DomainSize_D(iDir)
          end do
          !\
          ! call interpolate_amr
          !/
          if(UseGeneric)then
             call interpolate_amr(&
                  nDim=nDim, &
                  XyzIn_D=Coord_D, &
                  nIndexes=nDim+1,&
                  find=find_test, &
                  nCell_D=nCell_D(1:nDim),&
                  nGridOut=nGridOut,&
                  Weight_I=Weight_I,&
                  iIndexes_II=iIndexes_II,&
                  IsSecondOrder=IsSecondOrder,&
                  UseGhostCell=UseGhostCell)
          else
             call find_test(nDim, Coord_D, &
                  iProc, iBlock, CoordCorner_D, CellSize_D, IsOut)
             Coord_D = CoordCorner_D + Coord_D
             call check_interpolate_test(nDim, Coord_D, iBlock, &
                  iProc, iBlockNei)
             if(iBlockNei /= iBlock) then
                iBlock = iBlockNei
                CoordCorner_D = Coord_DB(:,iBlock)
                CellSize_D = DCoord_DB(:,iBlock) / nCell_D(1:nDim)
             end if
             call fix_coord(nDim, iBlock, Coord_D, CoordPass_D)
             call interpolate_amr_gc(&
                  nDim         = nDim, &
                  Xyz_D        = CoordPass_D, &
                  XyzMin_D     = CoordCorner_D, &
                  DXyz_D       = CellSize_D, &
                  nCell_D      = nCell_D(1:nDim), &
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
          CoordInterpolated_D = 0
          VarInterpolated   = 0
          do iGrid = 1, nGridOut
             iCell_D = 1
             iCell_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
             iBlock = iIndexes_II(nIndexes,iGrid)
             CoordInterpolated_D = CoordInterpolated_D + Weight_I(iGrid)*&
                  Coord_DGB(:,iCell_D(1), iCell_D(2), &
                  iCell_D(3), iBlock)
             VarInterpolated = VarInterpolated + &
                  Weight_I(iGrid)*&
                  Var_GB(iCell_D(1), iCell_D(2), &
                  iCell_D(3), iBlock)
          end do
          ! for periodic boundary conditions apply modulo
          CoordModulo_D = CoordInterpolated_D
          where(IsPeriodic_D(1:nDim))
             CoordModulo_D = modulo(CoordModulo_D, DomainSize_D)
          end where
          if(IsSpherical)then
             if(CoordModulo_D(Theta_) < CoordMin_D(Theta_))then
                CoordModulo_D(Theta_) = 2*CoordMin_D(Theta_) - CoordModulo_D(Theta_)
                CoordModulo_D(Phi_  ) = modulo(CoordModulo_D(Phi_)+0.5*DOmainSize_D(Phi_), DomainSize_D(Phi_))
             elseif(CoordModulo_D(Theta_) > CoordMax_D(Theta_))then
                CoordModulo_D(Theta_) = 2*CoordMax_D(Theta_) - CoordModulo_D(Theta_)
                CoordModulo_D(Phi_  ) = modulo(CoordModulo_D(Phi_)+0.5*DOmainSize_D(Phi_), DomainSize_D(Phi_))
             end if
          end if
          if(  DoTestApproximation.and.&
               any(abs(Coord_D - CoordModulo_D) > 1.0e-6).and.&
               IsSecondOrder)then
             write(*,*)'Approximation test failed'
             write(*,*)'nDim=',nDim
             write(*,*)'IsCartesian=',IsCartesian
             write(*,*)'IsPeriodic_D=',IsPeriodic_D
             write(*,*)'UseGeneric=',UseGeneric
             write(*,*)'UseGhostCell=',UseGhostCell
             write(*,*)'Grid:', iLevelTest_I(1:2**nDim)
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'Point=', Coord_D
             write(*,*)'Cell_D  iBlock CoordGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCell_D = 1
                iCell_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Coord_DGB(:,iCell_D(1), iCell_D(2), &
                     iCell_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'Coord_D=',Coord_D
             write(*,*)'CoordInterpolated_D=',CoordInterpolated_D
             call CON_stop('Correct code and redo test')
          end if
          !\
          ! Test continuity
          !/
          do iDir =1, nDim
             CoordCont_D(iDir) = Coord_D(iDir) + &
                  CellSizeCoarse_D(iDir)*(0.02*random_real(iSeed) - 0.01)
          end do
          !\
          ! call interpolate_amr
          !/
          if(UseGeneric)then
             call interpolate_amr(&
                  nDim=nDim, &
                  XyzIn_D=CoordCont_D, &
                  nIndexes=nDim+1,&
                  find=find_test, &
                  nCell_D=nCell_D(1:nDim),&
                  nGridOut=nGridOut,&
                  Weight_I=Weight_I,&
                  iIndexes_II=iIndexes_II,&
                  IsSecondOrder=IsSecondOrder,&
                  UseGhostCell=UseGhostCell)
          else
             call find_test(nDim, CoordCont_D, &
                  iProc, iBlock, CoordCorner_D, CellSize_D, IsOut)
             CoordCont_D = CoordCorner_D + CoordCont_D
             call check_interpolate_test(nDim, CoordCont_D, iBlock, &
                  iProc, iBlockNei)
             if(iBlockNei /= iBlock) then
                iBlock = iBlockNei
                CoordCorner_D = Coord_DB(:,iBlock)
                CellSize_D = DCoord_DB(:,iBlock) / nCell_D(1:nDim)
             end if
             call fix_coord(nDim, iBlock, CoordCont_D, CoordPass_D)
             call interpolate_amr_gc(&
                  nDim         = nDim, &
                  Xyz_D        = CoordPass_D, &
                  XyzMin_D     = CoordCorner_D, &
                  DXyz_D       = CellSize_D, &
                  nCell_D      = nCell_D(1:nDim), &
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
          CoordInterpolated_D = 0
          do iGrid = 1, nGridOut
             iCell_D = 1
             iCell_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
             iBlock = iIndexes_II(nIndexes,iGrid)
             VarContInterpolated = VarContInterpolated + &
                  Weight_I(iGrid)*&
                  Var_GB(iCell_D(1), iCell_D(2), &
                  iCell_D(3), iBlock)
             CoordInterpolated_D = CoordInterpolated_D + Weight_I(iGrid)*&
                  Coord_DGB(:,iCell_D(1), iCell_D(2), &
                  iCell_D(3), iBlock)
          end do
          CoordModulo_D = CoordInterpolated_D
          where(IsPeriodic_D(1:nDim))
             CoordModulo_D = modulo(CoordModulo_D, DomainSize_D)
          end where
          if(IsSpherical)then
             if(CoordModulo_D(Theta_) < CoordMin_D(Theta_))then
                CoordModulo_D(Theta_) = 2*CoordMin_D(Theta_) - CoordModulo_D(Theta_)
                CoordModulo_D(Phi_  ) = modulo(CoordModulo_D(Phi_)+0.5*DOmainSize_D(Phi_), DomainSize_D(Phi_))
             elseif(CoordModulo_D(Theta_) > CoordMax_D(Theta_))then
                CoordModulo_D(Theta_) = 2*CoordMax_D(Theta_) - CoordModulo_D(Theta_)
                CoordModulo_D(Phi_  ) = modulo(CoordModulo_D(Phi_)+0.5*DOmainSize_D(Phi_), DomainSize_D(Phi_))
             end if
          end if
          if(  DoTestApproximation.and.&
               any(abs(CoordCont_D - CoordModulo_D) > 1.0e-6).and.&
               IsSecondOrder)then
             write(*,*)'Approximation test failed'
             write(*,*)'nDim=',nDim
             write(*,*)'IsCartesian=',IsCartesian
             write(*,*)'IsPeriodic_D=',IsPeriodic_D
             write(*,*)'UseGeneric=',UseGeneric
             write(*,*)'UseGhostCell=',UseGhostCell
             write(*,*)'Grid:', iLevelTest_I(1:2**nDim)
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'PointCont=', CoordCont_D
             write(*,*)'Cell_D  iBlock CoordGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCell_D = 1
                iCell_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Coord_DGB(:,iCell_D(1), iCell_D(2), &
                     iCell_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'CoordCont_D=',CoordCont_D
             write(*,*)'CoordInterpolated_D=',CoordInterpolated_D
             call CON_stop('Correct code and redo test')
          end if
          if(abs(VarContInterpolated - VarInterpolated) > nDim*0.01)then
             write(*,*)'Continuity test failed'
             write(*,*)'nDim=',nDim
             write(*,*)'IsCartesian=',IsCartesian
             write(*,*)'IsPeriodic_D=',IsPeriodic_D
             write(*,*)'UseGeneric=',UseGeneric
             write(*,*)'UseGhostCell=',UseGhostCell
             write(*,*)'Grid:', iLevelTest_I
             write(*,*)'nGridOut=',nGridOut
             write(*,*)'CoordCont=', CoordCont_D
             write(*,*)'Cell_D  iBlock CoordGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCell_D = 1
                iCell_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Coord_DGB(:,iCell_D(1), iCell_D(2), &
                     iCell_D(3), iBlock), Weight_I(iGrid)
             end do
             write(*,*)'Coord_D=',Coord_D
             call interpolate_amr(&
                  nDim=nDim, &
                  XyzIn_D=Coord_D, &
                  nIndexes=nDim+1,&
                  find=find_test, &
                  nCell_D=nCell_D(1:nDim),&
                  nGridOut=nGridOut,&
                  Weight_I=Weight_I,&
                  iIndexes_II=iIndexes_II,&
                  IsSecondOrder=IsSecondOrder)
             write(*,*)'Cell_D  iBlock CoordGrid_D Weight_I(iGrid)'
             do iGrid = 1, nGridOut
                iCell_D = 1
                iCell_D(1:nDim) = iIndexes_II(1:nDim,iGrid)
                iBlock = iIndexes_II(nIndexes,iGrid)
                write(*,*)iIndexes_II(1:nDim,iGrid), iBlock ,&
                     Coord_DGB(:,iCell_D(1), iCell_D(2), &
                     iCell_D(3), iBlock), Weight_I(iGrid)
             end do
             call CON_stop('Correct code and redo test')
          end if
       end do SAMPLE
    end do CASE
    call clean_grid
  contains
    !============================
    subroutine clean_grid
      deallocate(CoordMin_D)
      deallocate(CoordMax_D)
      deallocate(DomainSize_D)
      deallocate(Coord_DB)
      deallocate(DCoord_DB)
      deallocate(Used_B)  
      deallocate(iTree_IB)
      deallocate(Coord_DGB, Var_GB, DiLevelNei_IIIB)
    end subroutine clean_grid
    !============================
    subroutine generate_grid
      ! generates the grid used in the test
      integer:: nBlockInit, iBlockChild, iChild
      !------------------------------------------
      allocate(CoordMin_D(1:nDim))
      allocate(CoordMax_D(1:nDim))
      allocate(DomainSize_D(1:nDim))
      if(IsCartesian)then
         ! number of root blocks
         nRoot  = 1
         ! total number of blocks
         nBlock = (2**nDim) * (2**nDim+1) + 1
         ! domains boundaries and size
         DomainSize_D = 2*nCell_D(1:nDim)*SizeRatio_D(1:nDim)
         CoordMin_D   = 0.0
         CoordMax_D   = DomainSize_D


         ! allocate and initialize (where needed) grid storage
         allocate( Coord_DB(1:nDim,1:nBlock))
         allocate(DCoord_DB(1:nDim,1:nBlock))
         allocate(Used_B(1:nBlock))
         Used_B = .false.
         allocate(iTree_IB(1:Level_,1:nBlock))
         iTree_IB = -1
         CellSizeCoarse_D = 1*SizeRatio_D(1:nDim)
         CellSizeFine_D   = 0.5*SizeRatio_D(1:nDim)
         ! initialize the root block
         Coord_DB( :,1)= CoordMin_D
         DCoord_DB(:,1)= DomainSize_D

         ! at the moment only one block is initialized
         nBlockInit = 1
      else ! spherical or polar
         ! number of root blocks
         nRoot  = 2
         ! total number of blocks
         nBlock = 2*((2**nDim) * (2**nDim+1) + 1)
         ! domains boundaries and size
         DomainSize_D(1:2) = (/2.0, 2*cPi/)
         CoordMin_D(1)     = 1.0
         CoordMin_D(2:nDim)= 0.0
         if(nDim==3) DomainSize_D(nDim) = cPi
         CoordMax_D   = CoordMin_D + DomainSize_D

         ! allocate and initialize (where needed) grid storage
         allocate( Coord_DB(1:nDim,1:nBlock))
         allocate(DCoord_DB(1:nDim,1:nBlock))
         allocate(Used_B(1:nBlock))
         Used_B = .false.
         allocate(iTree_IB(1:Level_,1:nBlock))
         iTree_IB = -1
         CellSizeCoarse_D(1)     = 1.0 / nCell_D(1)
         CellSizeCoarse_D(2:nDim)= 0.5 * cPi / nCell_D(2:nDim)
         CellSizeFine_D          = 0.5*CellSizeCoarse_D
         ! initialize the root blocks
         Coord_DB( :,1)= CoordMin_D
         DCoord_DB(:,1)= DomainSize_D
         DCoord_DB(2,1)= cPi
         Coord_DB( :,2)= CoordMin_D
         Coord_DB( 2,2)= CoordMin_D(2) + cPi
         DCoord_DB(:,2)= DCoord_DB(:,1)

         ! at the moment only one block is initialized
         nBlockInit = 2
      end if

      !\
      ! Refine block up to level Fine_
      ! ----------------------------------------
      ! the index of the block to be refined
      iBlock = 1
      do while(nBlockInit < nBlock)
         ! don't refine Fine_ blocks
         if(iTree_IB(Level_,iBlock)==Fine_)&
              CYCLE
         do iChild = 1, 2**nDim
            iBlockChild = nBlockInit + iChild
            ! store indax of iChild for iBlock
            iTree_IB(iChild,iBlock) = nBlockInit + iChild
            ! initialize iChild
            ! tree data
            iTree_IB(Level_, iBlockChild) = iTree_IB(Level_,iBlock) + 1
            ! coordinate data
            DCoord_DB(1:nDim,iBlockChild) = 0.5 * DCoord_DB(1:nDim,iBlock)
            Coord_DB(1:nDim, iBlockChild) = Coord_DB(1:nDim, iBlock) + &
                 iShift_DI(1:nDim, iChild) * DCoord_DB(1:nDim,iBlockChild)
         end do
         nBlockInit = nBlockInit + 2**nDim
         iBlock     = iBlock + 1
      end do

      ! initialize additional arrays for all blocks
      allocate(DiLevelNei_IIIB(-1:1,-1:1,-1:1,nBlock))
      DiLevelNei_IIIB = 0
      allocate(Coord0_DGB(nDim, 1-nG:nCell_D(1)+nG, &
           1-nG:nCell_D(2)+nG, &
           1-nG:nCell_D(3)+nG, nBlock))
      Coord0_DGB = 0
      allocate(Coord_DGB(nDim, 1-nG:nCell_D(1)+nG, &
           1-nG:nCell_D(2)+nG, &
           1-nG:nCell_D(3)+nG, nBlock))
      Coord_DGB = 0
      allocate(Var_GB(1-nG:nCell_D(1)+nG, &
           1-nG:nCell_D(2)+nG, &
           1-nG:nCell_D(3)+nG, nBlock))
      Var_GB = 0

      !\
      ! fills cells
      !/
      do iBlock = 1, nBlock
         do k = 1-nG, nCell_D(3)+nG 
            do j = 1-nG, nCell_D(2)+nG 
               do i = 1-nG, nCell_D(1)+nG
                  iCell_D = (/i,j,k/)
                  ! coords of cell centers
                  Coord_DGB(:,i,j,k,iBlock) = Coord_DB(:,iBlock) +&
                       DCoord_DB(:,iBlock) * (iCell_D(1:nDim) - 0.50) / nCell_D(1:nDim)
                  ! for now, fill values at physical cells only
                  if(all(iCell_D >= 1).and.all(iCell_D <= nCell_D))then
                     select case(iTree_IB(Level_,iBlock))
                     case(Fine_)
                        Var_GB(i,j,k,iBlock) = 0.25 + 0.50*random_real(iSeed) 
                     case(Coarse_)
                        Var_GB(i,j,k,iBlock) = random_real(iSeed) 
                     case default
                        ! do nothing   
                     end select
                  end if
               end do
            end do
         end do
      end do
      Coord0_DGB = Coord_DGB
    end subroutine generate_grid
    !============================================================
    subroutine set_levels(iCase)
      integer, intent(in) :: iCase

      integer:: iMisc, iChild, iUseFine, iRoot, iBlockChild
      !------------------------------------
      ! reset
      Used_B = .false.
      do iRoot = 1, nRoot
         Used_B(iTree_IB(1:2**nDim,iRoot)) = .true.
      end do

      iMisc = iCase
      iRoot = 1 + iMisc / 2**(2**nDim)
      iMisc = mod(iMisc, 2**(2**nDim))
      iChild = 0
      do while(iMisc > 0)
         iChild = iChild + 1
         iBlockChild = iTree_IB(iChild,iRoot)
         iUseFine = mod(iMisc, 2)
         iMisc = (iMisc - iUseFine)/2
         iLevelTest_I((iRoot-1)*2**nDim+iChild) = iUseFine
         if(iUseFine == 1)then
            Used_B(iBlockChild) = .true.
         else
            Used_B(iBlockChild) = .false.
            Used_B(iTree_IB(1:2**nDim,iBlockChild)) = .true.
         end if
      end do
    end subroutine set_levels
    !============================
    subroutine fix_coord(nDim, iBlock, CoordIn_D, CoordOut_D)
      ! for periodic/flipped coordinate may need to adjust point's coordinates
      ! for calling interpolate_amr_gc; 
      ! fix depends on reference block used in interpolate_amr_gc
      integer,intent(in) :: nDim
      integer,intent(in) :: iBlock
      real,   intent(in) :: CoordIn_D(nDim)
      real,   intent(out):: CoordOut_D(nDim)

      real, dimension(nDim) :: CellSize_D, CoordBlockMin_D, CoordBlockMax_D, iDiscr_D
      logical:: DoFix_D(nDim)
      real:: PhiAux
      real:: PhiAuxLo, PhiAuxHi
      integer:: iDim
      real, parameter:: cTol = 1E-8
      !------------------------------------------------------------------------
      CoordOut_D = CoordIn_D
      CoordBlockMin_D = Coord_DB(:,iBlock)
      CoordBlockMax_D = Coord_DB(:,iBlock) + DCoord_DB(:, iBlock)
      if(  all(CoordOut_D < CoordBlockMax_D) .and. &
           all(CoordOut_D >=CoordBlockMin_D))&
           RETURN
      CellSize_D = (CoordBlockMax_D - CoordBlockMin_D)/ nCell_D(1:nDim)
      do iDim = 1, nDim
         if(.not.IsPeriodic_D(iDim)) CYCLE
         if(CoordBlockMax_D(iDim)==CoordMax_D(iDim).and.&
              CoordOut_D(iDim) - CoordMin_D(iDim) < CellSize_D(iDim)*(1+cTol))then
            CoordOut_D(iDim) = CoordOut_D(iDim) + DomainSize_D(iDim)
         elseif(CoordBlockMin_D(iDim)==CoordMin_D(iDim).and.&
              CoordMax_D(iDim) - CoordOut_D(iDim) < CellSize_D(iDim)*(1+cTol))then
            CoordOut_D(iDim) = CoordOut_D(iDim) - DomainSize_D(iDim)
         end if
      end do
      if(IsSpherical)then
         PhiAuxLo =  CoordOut_D(Phi_) - 0.5*DomainSize_D(Phi_)
         PhiAuxHi =  CoordOut_D(Phi_) + 0.5*DomainSize_D(Phi_)
         if(  PhiAuxLo >= CoordBlockMin_D(Phi_) - CellSize_D(Phi_)*(1+cTol) .and.&
              PhiAuxLo <  CoordBlockMax_D(Phi_) + CellSize_D(Phi_)*(1+cTol))then
            PhiAux = PhiAuxLo
         elseif(PhiAuxHi >= CoordBlockMin_D(Phi_) - CellSize_D(Phi_)*(1+cTol) .and.&
              PhiAuxHi <  CoordBlockMax_D(Phi_) + CellSize_D(Phi_)*(1+cTol))then
            PhiAux = PhiAuxHi
         else
            RETURN
         end if
         if(CoordBlockMax_D(Theta_)==CoordMax_D(Theta_))then
            CoordOut_D(Theta_) = 2*CoordMax_D(Theta_) - CoordOut_D(Theta_)
            CoordOut_D(Phi_  ) = PhiAux
         elseif(CoordBlockMin_D(Theta_)==CoordMin_D(Theta_))then
            CoordOut_D(Theta_) = 2*CoordMin_D(Theta_) - CoordOut_D(Theta_)
            CoordOut_D(Phi_  ) = PhiAux
         end if
      end if
    end subroutine fix_coord
    !============================
    subroutine reverse_fix(nDim, iBlock, CoordIn_D, CoordOut_D)
      ! for periodic/flipped coordinate may need to adjust point's coordinates
      ! for calling interpolate_amr_gc; 
      ! fix depends on reference block used in interpolate_amr_gc
      integer,intent(in) :: nDim
      integer,intent(in) :: iBlock
      real,   intent(in) :: CoordIn_D(nDim)
      real,   intent(out):: CoordOut_D(nDim)

      real   :: CellSize, CoordBlockMin_D(nDim), CoordBlockMax_D(nDim)
      integer:: iDim
      !------------------------------------------------------------------------
      CoordOut_D = CoordIn_D
      CoordBlockMin_D = Coord_DB(:,iBlock)
      CoordBlockMax_D = Coord_DB(:,iBlock) + DCoord_DB(:, iBlock)
      if(  all(CoordOut_D < CoordBlockMax_D) .and. &
           all(CoordOut_D >=CoordBlockMin_D))&
           RETURN

      do iDim = 1, nDim
         if(.not.IsPeriodic_D(iDim)) CYCLE
         CellSize = &
              (CoordBlockMax_D(iDim) - CoordBlockMin_D(iDim))/ nCell_D(iDim)
         if(CoordBlockMax_D(iDim)==CoordMax_D(iDim).and.&
              CoordOut_D(iDim) - CoordMin_D(iDim) < CellSize)then
            CoordOut_D(iDim) = CoordOut_D(iDim) + DomainSize_D(iDim)
         elseif(CoordBlockMin_D(iDim)==CoordMin_D(iDim).and.&
              CoordMax_D(iDim) - CoordOut_D(iDim) < CellSize)then
            CoordOut_D(iDim) = CoordOut_D(iDim) - DomainSize_D(iDim)
         end if
      end do

      if(IsSpherical)then
         if(CoordBlockMax_D(Theta_)==CoordMax_D(Theta_).and.&
              abs(CoordOut_D(Phi_)-CoordBlockMax_D(Phi_)) > &
              0.25*DomainSize_D(Phi_))then
            CoordOut_D(Theta_) = 2*CoordMax_D(Theta_) - CoordOut_D(Theta_)
            CoordOut_D(Phi_  ) = CoordMin_D(Phi_) + modulo(&
                 CoordOut_D(Phi_) - CoordMin_D(Phi_) + 0.5*DomainSize_D(Phi_),&
                 DomainSize_D(Phi_))
         elseif(CoordBlockMin_D(Theta_)==CoordMin_D(Theta_).and.&
              abs(CoordOut_D(Phi_)-CoordBlockMin_D(Phi_)) > &
              0.25*DomainSize_D(Phi_))then
            CoordOut_D(Theta_) = 2*CoordMin_D(Theta_) - CoordOut_D(Theta_)
            CoordOut_D(Phi_  ) = CoordMin_D(Phi_) + modulo(&
                 CoordOut_D(Phi_) - CoordMin_D(Phi_) + 0.5*DomainSize_D(Phi_),&
                 DomainSize_D(Phi_))
         end if
      end if
    end subroutine reverse_fix
    !============================
    subroutine check_interpolate_test(nDim, Coord_D, iBlockIn, &
         iProcOut, iBlockOut)
      integer, intent(in) :: nDim
      real,    intent(in) :: Coord_D(nDim)
      integer, intent(in) :: iBlockIn
      integer, intent(out):: iProcOut
      integer, intent(out):: iBlockOut

      integer:: iLevel_I(2**nDim)
      logical:: IsOut_I(2**nDim)
      integer:: iDiscr_D(3)
      real   :: CoordCentral_D(nDim)
      real   :: CoordNei_D(nDim), CoordCorner_D(nDim), CellSize_D(nDim)
      real   :: CoordGrid_DI(nDim, 2**nDim)
      logical:: IsOut
      integer:: iGridRef, iGrid
      integer:: iShiftRef_D(nDim)
      !------------------------------------------------------------
      ! determine displacement of the point relative to the block's interior
      iDiscr_D = 0
      where(Coord_D < Coord_DGB(:,1,1,1,iBlockIn))
         iDiscr_D(1:nDim) = -1
      elsewhere(Coord_D >= Coord_DGB(:,nCell_D(1),nCell_D(2),nCell_D(3),iBlockIn))
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
      if(iTree_IB(Level_,iBlockIn)==Fine_)then
         CellSize_D = CellSizeFine_D   ! the block is Finer 
      else
         CellSize_D = CellSizeCoarse_D ! the block is Coarse
      end if

      ! coordinates of block's junction
      where(    iDiscr_D(1:nDim) == 1)
         CoordCentral_D = &
              Coord_DGB(:,nCell_D(1),nCell_D(2),nCell_D(3),iBlockIn) + 0.5*CellSize_D
      elsewhere(iDiscr_D(1:nDim) ==-1)
         CoordCentral_D = Coord_DGB(:,1,1,1,iBlockIn) - 0.5*CellSize_D
      elsewhere
         CoordCentral_D = Coord_D
      end where

      ! supergrid
      do iGrid = 1, 2**nDim
         CoordGrid_DI(:, iGrid) = CoordCentral_D + &
              (iShift_DI(1:nDim, iGrid)-0.5) * CellSizeCoarse_D
      end do

      ! find the block that has enough information to perform interpolation
      ! using only 1 layer of gc
      call get_reference_block(&
           nDim, Coord_D, CoordGrid_DI, iLevel_I, IsOut_I, iGridRef)
      ! displacement to this block
      iShiftRef_D = iShift_DI(1:nDim,iGridRef)

      ! check if it is the input block
      if(all((iDiscr_D(1:nDim)-1)/2 == iShiftRef_D*ABS(iDiscr_D(1:nDim))))then
         iBlockOut = iBlockIn
         RETURN
      end if

      ! it is a different block, find it
      CoordNei_D = CoordCentral_D + 0.5*CellSizeFine_D*(iShiftRef_D-0.5)
      call find_test(nDim, CoordNei_D, &
           iProcOut, iBlockOut, CoordCorner_D, CellSize_D, IsOut)
    end subroutine check_interpolate_test
    !===========================
    subroutine fill_nei_levels
      ! fill DiLevelNei_IIIB that store levels of neighbors
      integer:: iLevelNei, iLevel
      real, dimension(nDim):: Coord_D, CellSize_D, CoordCenter_D
      integer:: i,j,k, iIndex_I(3)
      !-----------------------------------------------
      DiLevelNei_IIIB = 0
      do iBlock = 1, nBlock
         ! skip the unused
         if(.not.Used_B(iBlock))&
              CYCLE

         iLevel = iTree_IB(Level_, iBlock)
         ! search neigbors
         CoordCenter_D = Coord_DB(:,iBlock) + 0.5 * DCoord_DB(:,iBlock)
         do i = -1, 1; do j = -1, 1; do k = -(nDim-2), (nDim-2)
            iIndex_I = (/i,j,k/)
            Coord_D = CoordCenter_D + 0.6*iIndex_I(1:nDim) *Dcoord_DB(:,iBlock)
            call find_test(nDim, Coord_D, &
                 iProc, iBlockNei, CoordCorner_D, CellSize_D, IsOut)
            if(IsOut)then
               iLevelNei = Out_
            else
               iLevelNei = iTree_IB(Level_, iBlockNei)
            end if
            DiLevelNei_IIIB(i,j,k,iBlock) = iLevelNei - iLevel
         end do; end do; end do
      end do
    end subroutine fill_nei_levels
    !============================
    subroutine fill_ghost_cells
      ! values at ghost cells depend on the refinement
      logical:: IsOut
      integer:: iProc, iBlock, iBlockNei, i, j ,k
      integer:: iCell_D(3), iCellIndexNei_D(3)
      real, dimension(nDim):: Coord_D, CoordCorner_D, CellSize_D, DcoordModulo_D
      !-----------------------------------------------
      Coord_DGB = Coord0_DGB
      ! fill values based on whether bocks are used
      do iBlock = 1, nBlock
         if(.not.Used_B(iBlock))CYCLE
         ! loop over cells
         do k = 1-nG, nCell_D(3)+nG
            do j = 1-nG, nCell_D(2)+nG
               do i = 1-nG, nCell_D(1)+nG
                  iCell_D = (/i,j,k/)
                  ! skip physical cells
                  if(all(iCell_D >= 1).and.all(iCell_D <= nCell_D))&
                       CYCLE
                  Coord_D = Coord_DGB(:,i,j,k,iBlock)
                  call find_test(nDim, Coord_D, &
                       iProc, iBlockNei, CoordCorner_D, CellSize_D, IsOut)
                  DcoordModulo_D = Coord_DGB(:,i,j,k,iBlock) - Coord_D - CoordCorner_D
                  where(    DcoordModulo_D > CellSize_D.and.IsPeriodic_D(1:nDim))
                     DcoordModulo_D = DomainSize_D
                  elsewhere(DcoordModulo_D <-CellSize_D.and.IsPeriodic_D(1:nDim))
                     DcoordModulo_D =-DomainSize_D
                  elsewhere
                     !DcoordModulo_D = 0.0
                  end where
                  if(IsOut) CYCLE
                  ! don't fill ghost cells for Coarse_ block if nei is Fine_
                  if(  iTree_IB(Level_,iBlock   )==Coarse_ .and.&
                       iTree_IB(Level_,iBlockNei)==Fine_)&
                       CYCLE
                  iCellIndexNei_D = 1
                  iCellIndexNei_D(1:nDim) = &
                       nint(0.3+Coord_D(1:nDim)/CellSize_D(1:nDim))
                  Var_GB(i,j,k,iBlock) = &
                       Var_GB(iCellIndexNei_D(1),&
                       iCellIndexNei_D(2),&
                       iCellIndexNei_D(3),iBlockNei)
                  Coord_DGB(:,i,j,k,iBlock) = &
                       Coord_DGB(:,iCellIndexNei_D(1),&
                       iCellIndexNei_D(2),&
                       iCellIndexNei_D(3),iBlockNei)
                  call fix_coord(nDim, iBlock, Coord_DGB(:,i,j,k,iBlock), Coord_D)
                  Coord_DGB(:,i,j,k,iBlock) = Coord_D
               end do
            end do
         end do
      end do
    end subroutine fill_ghost_cells
  end subroutine test_interpolate_amr
  !============================
  subroutine find_test(nDim, Coord_D, &
       iProc, iBlock, CoordCorner_D, CellSize_D, IsOut)
    integer, intent(in) :: nDim
    !\
    ! "In"- the coordinates of the point, "out" the coordinates of the
    ! point with respect to the block corner. In the most cases
    ! CoordOut_D = CoordIn_D - CoordCorner_D, the important distinction,
    ! however, is the periodic boundary, near which the jump in the
    ! stencil coordinates might occur. To handle the latter problem,
    ! we added the "out" intent. The coordinates for the stencil
    ! and input point are calculated and recalculated below with
    ! respect to the block corner. 
    !/
    real,  intent(inout):: Coord_D(nDim)
    integer, intent(out):: iProc, iBlock !processor and block number
    !\
    ! Block left corner coordinates and the grid size:
    !/
    real,    intent(out):: CoordCorner_D(nDim), CellSize_D(nDim)
    logical, intent(out):: IsOut !Point is out of the domain.
    real, dimension(nDim):: CoordCenter_D
    integer:: iShift_D(3)
    integer, dimension(0:1,0:1,0:1), parameter:: iChildFromShift_III=reshape(&
         (/1, 2, 3, 4, 5, 6, 7, 8/),(/2, 2, 2/))
    logical, dimension(nDim) :: IsAboveCenter_D
    integer:: iRoot, iChild
    !------------------- 
    iProc = 0
    iBlock= 0
    CoordCorner_D = CoordMin_D
    CellSize_D    = 0.0
    ! fix periodic coordinates
    where(IsPeriodic_D(1:nDim))
       Coord_D = modulo(Coord_D, DomainSize_D)
    end where
    ! fix "flipped" coordinate
    if(IsSpherical)then
       if(Coord_D(Theta_) < CoordMin_D(Theta_))then
          Coord_D(Theta_) = 2*CoordMin_D(Theta_) - Coord_D(Theta_)
          Coord_D(Phi_  ) = &
               modulo(Coord_D(Phi_)+0.5*DomainSize_D(Phi_),DomainSize_D(Phi_))
       elseif(Coord_D(Theta_) > CoordMax_D(Theta_))then
          Coord_D(Theta_) = 2*CoordMax_D(Theta_) - Coord_D(Theta_)
          Coord_D(Phi_  ) = &
               modulo(Coord_D(Phi_)+0.5*DomainSize_D(Phi_),DomainSize_D(Phi_))
       end if
    end if
    ! check whether the point is out if the domain
    IsOut = any(Coord_D < CoordMin_D .or. Coord_D >= CoordMax_D)
    if(IsOut) RETURN
    !\
    ! determine which root block the point falls into
    !/
    do iRoot = 1, nRoot
       if(  all(Coord_D >= Coord_DB(:,iRoot)) .and. &
            all(Coord_D <  Coord_DB(:,iRoot)+DCoord_DB(:,iRoot)))&
            EXIT
    end do
    !\
    ! Find into which coarse block the point fall
    !/ 
    CoordCenter_D = Coord_DB(:,iRoot) + 0.5 * DCoord_DB(:,iRoot)
    IsAboveCenter_D = Coord_D >= CoordCenter_D
    iShift_D = 0
    where(IsAboveCenter_D)iShift_D(1:nDim) = 1
    iChild = iChildFromShift_III(iShift_D(1),iShift_D(2),iShift_D(3))
    iBlock = iTree_IB(iChild,iRoot)
    !\
    ! Check if the coarse block is used
    !/
    if(Used_B(iBlock))then
       CoordCorner_D = Coord_DB(:,iBlock)
       Coord_D       = Coord_D - CoordCorner_D
       CellSize_D    = DCoord_DB(:,iBlock) / nCell_D(1:nDim)
       RETURN
    end if
    !\
    ! The coarser block is refined, find into which fine block 
    ! the point falls
    !/
    CoordCenter_D = Coord_DB(:,iBlock) + 0.5 * DCoord_DB(:,iBlock)
    IsAboveCenter_D = Coord_D >= CoordCenter_D
    iShift_D = 0
    where(IsAboveCenter_D)iShift_D(1:nDim) = 1
    iChild = iChildFromShift_III(iShift_D(1),iShift_D(2),iShift_D(3))
    iBlock = iTree_IB(iChild,iBlock)
    CoordCorner_D = Coord_DB(:,iBlock)
    Coord_D       = Coord_D - CoordCorner_D
    CellSize_D    = DCoord_DB(:,iBlock) / nCell_D(1:nDim)
  end subroutine find_test
end module ModTestInterpolateAMR
!=============================================================================
!=============================================================================
program test_interpolate_amr

  use ModTestInterpolateAMR, test => test_interpolate_amr

  implicit none

  integer :: nSampleTest = 1000

  call test(&
       nDim          = 2,&
       IsCartesian   = .false.,&
       IsPeriodicIn_D= (/.true.,.false./), &
       nSample       = nSampleTest, &
       UseGeneric    = .true., &
       UseGhostCell  = .false.)

  call test(&
       nDim          = 2,&
       IsCartesian   = .false.,&
       IsPeriodicIn_D= (/.false.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .true., &
       UseGhostCell  = .true.)

  call test(&
       nDim          = 2,&
       IsCartesian   = .false.,&
       IsPeriodicIn_D= (/.false.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .false.)

  call test(&
       nDim          = 3,&
       IsCartesian   = .false.,&
       IsPeriodicIn_D= (/.false.,.true.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .true., &
       UseGhostCell  = .false.)

  call test(&
       nDim          = 3,&
       IsCartesian   = .false.,&
       IsPeriodicIn_D= (/.false.,.false.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .true., &
       UseGhostCell  = .true.)

  call test(&
       nDim          = 3,&
       IsCartesian   = .false.,&
       IsPeriodicIn_D= (/.true.,.false.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .false.)



  call test(&
       nDim          = 2,&
       IsCartesian   = .true.,&
       IsPeriodicIn_D= (/.true.,.false./), &
       nSample       = nSampleTest, &
       UseGeneric    = .true., &
       UseGhostCell  = .false.)

  call test(&
       nDim          = 2,&
       IsCartesian   = .true.,&
       IsPeriodicIn_D= (/.false.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .true., &
       UseGhostCell  = .true.)

  call test(&
       nDim          = 2,&
       IsCartesian   = .true.,&
       IsPeriodicIn_D= (/.false.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .false.)


  call test(&
       nDim          = 3,&
       IsCartesian   = .true.,&
       IsPeriodicIn_D= (/.false.,.true.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .true., &
       UseGhostCell  = .false.)

  call test(&
       nDim          = 3,&
       IsCartesian   = .true.,&
       IsPeriodicIn_D= (/.false.,.false.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .true., &
       UseGhostCell  = .true.)

  call test(&
       nDim          = 3,&
       IsCartesian   = .true.,&
       IsPeriodicIn_D= (/.true.,.false.,.true./), &
       nSample       = nSampleTest, &
       UseGeneric    = .false.)

end program test_interpolate_amr
