!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE PT

!BOP
!MODULE: CON_couple_gm_pt - couple GM and PT components
!
!DESCRIPTION:
! Couple GM and PT components both ways. 
!
!INTERFACE:
module CON_couple_gm_pt

  !USES:
  use CON_coupler

  use CON_couple_points

  implicit none
  save

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_pt_init ! initialize both couplings
  public :: couple_gm_pt      ! couple GM to PT

  !REVISION HISTORY:
  ! 09/24/2013 G.Toth <gtoth@umich.edu> - initial version
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  logical       :: UseMe = .true.

  ! Check if the GM and PT processors coincide
  logical:: IsSameLayout

  ! Proc index inside SWMF
  integer:: iProcWorld

  ! Router communicator info
  integer:: iCommGmPt, nProcGmPt, iProcGmPt, iProc0Gm, iProc0Pt, nProcCommon

  type(CouplePointsType) :: Coupler
  
contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt_init - initialize GM-PT couplings
  !INTERFACE:
  subroutine couple_gm_pt_init

    integer:: iError
    !integer:: nDimPt

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='couple_gm_pt_init'

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s
    !EOP
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    Coupler%iCompTarget = PT_
    Coupler%iCompSource = GM_
    
    ! GM sends all its variables to PT. Take information from Grid_C
    Coupler%NameVar = Grid_C(Coupler%iCompSource)%NameVar
    Coupler%nVar    = Grid_C(Coupler%iCompSource)%nVar

    call couple_points_init(Coupler)

  end subroutine couple_gm_pt_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt - couple GM to PT
  !INTERFACE:
  subroutine couple_gm_pt(tSimulation)
    
    interface

       subroutine PT_put_from_gm( &
            NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

         implicit none
         ! List of variables
         character(len=*), intent(inout):: NameVar 

         ! Number of variables in Data_VI
         integer,          intent(inout):: nVar    

         ! Number of points in Pos_DI
         integer,          intent(inout):: nPoint  

         ! Position vectors
         real, pointer:: Pos_DI(:,:)               

         ! Recv data array
         real,    intent(in), optional:: Data_VI(:,:)

         ! Order of data
         integer, intent(in), optional:: iPoint_I(nPoint)    
       end subroutine PT_put_from_gm
       subroutine GM_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

         implicit none

         ! dimension of position vectors
         integer, intent(in) :: nDimIn             

         ! number of positions
         integer, intent(in) :: nPoint                

         ! positions
         real,    intent(in) :: Xyz_DI(nDimIn,nPoint) 

         ! processor owning position
         integer, intent(out):: iProc_I(nPoint)       

       end subroutine GM_find_points


        subroutine GM_get_for_pt(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
             Data_VI)
          logical,          intent(in):: IsNew   ! true for new point array
          character(len=*), intent(in):: NameVar ! List of variables
          integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
          integer,          intent(in):: nDimIn  ! Dimensionality of positions
          integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

          real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
          real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array
        end subroutine GM_get_for_pt

        subroutine GM_get_grid_info(nDim, iGrid, iDecomp)
          integer, intent(out) :: nDim, iGrid, iDecomp
        end subroutine GM_get_grid_info

        subroutine PT_get_grid_info(nDim, iGrid, iDecomp)
          integer, intent(out) :: nDim, iGrid, iDecomp
        end subroutine PT_get_grid_info
     end interface
    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere       (GM) source\\
    !    Particle Tracker           (PT) target
    !
    ! Send information from GM to PT. 
    !EOP

    ! Stored variables for GM->PT coupler
    !------------------------------------

    ! List of variables to pass
    character(len=lNameVar):: NameVar

    ! Grid index
    integer:: iDecompLastGm = -1, iDecompLastPt = -1
!    integer:: iGridGm, iDecompGm, iGridPt, iDecompPt

    ! Number of local points for PT and GM cores
!    integer:: nPointPt = 0, nPointGm = 0
    
    ! Number of data points found on GM component
!    integer:: nData = 0

    ! Number of points that belong to a given processor of the OTHER component
!    integer, allocatable, save:: nPointGm_P(:), nPointPt_P(:)

    ! Permutation of PT points after data is returned
!    integer, allocatable, save:: iPointPt_I(:)

    ! Point positions local on a GM processor
!    real, allocatable, save:: PosGm_DI(:,:)
    !---------------------------------------------------------------
    ! Temporary variables

    ! Is there a need to recalculate the data transfer route?
!    logical:: IsNewRoute

    ! Storage for original PT point positions. 
!    real, pointer:: PosPt_DI(:,:)

    ! Positions sorted according to the correspongin GM processors
!    real, allocatable:: PosSortPt_DI(:,:)

    ! GM processor index for PT points
!    integer, allocatable:: iProcPt_I(:)
!    integer, allocatable:: iProcGm_I(:)

    ! Buffers for data on GM and PT
!    real, allocatable:: DataGm_VI(:,:), DataPt_VI(:,:)

    integer:: iPoint, i, iError

    logical :: DoTest, DoTestMe

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_pt'

    ! Variables for general-case coupling

    ! number of processors of the OTHER component to communicate with
    integer, save:: nCoupleGm, nCouplePt

    ! processors of the OTHER component to communicate with
    integer, allocatable, save:: iCoupleProcGm_I(:), iCoupleProcPt_I(:)

    ! number of entries received/sent by a processor during rendezvous
    integer, allocatable, save:: nCouplePointGm_I(:), nCouplePointPt_I(:)
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',Coupler%iProcWorld

    call couple_points(Coupler, GM_get_for_pt, PT_put_from_gm, PT_get_grid_info, GM_get_grid_info, GM_find_points)

    if(DoTest) write(*,*) NameSub,' finished, iProc=',Coupler%iProcWorld
  end subroutine couple_gm_pt

end module CON_couple_gm_pt
