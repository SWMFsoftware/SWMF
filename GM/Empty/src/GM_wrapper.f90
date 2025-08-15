!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! Wrapper for the empty Global Magnetosphere (GM) component

module GM_wrapper

  use CON_coupler

  implicit none

contains
  !============================================================================

  subroutine GM_set_param(CompInfo, TypeAction)

    use CON_comp_info

    ! Arguments
    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*), intent(in)    :: TypeAction ! What to do
    character(len=*), parameter:: NameSub = 'GM_set_param'
    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.false., &
            NameVersion='Empty', &
            Version    =0.0)

    case default
       call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
    end select

  end subroutine GM_set_param
  !============================================================================

  subroutine GM_init_session(iSession, TimeSimulation)

    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'GM_init_session'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_init_session
  !============================================================================

  subroutine GM_finalize(TimeSimulation)

    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'GM_finalize'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_finalize
  !============================================================================

  subroutine GM_save_restart(TimeSimulation)

    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter:: NameSub = 'GM_save_restart'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_save_restart
  !============================================================================

  subroutine GM_run(TimeSimulation,TimeSimulationLimit)

    real, intent(inout):: TimeSimulation   ! current time of component

    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter:: NameSub = 'GM_run'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_run
  !============================================================================

  subroutine GM_get_grid_info(nDimOut, iGridOut, iDecompOut)

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index (increases with AMR)
    integer, intent(out):: iDecompOut ! decomposition index

    character(len=*), parameter:: NameSub = 'GM_get_grid_info'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_grid_info
  !============================================================================

  subroutine GM_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    integer, intent(in) :: nDimIn ! dimension of position vectors
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position

    character(len=*), parameter:: NameSub = 'GM_find_points'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_find_points
  !============================================================================

  subroutine GM_use_pointer(iComp, tSimulation)

    integer, intent(in):: iComp
    real,    intent(in):: tSimulation

    character(len=*), parameter:: NameSub = 'GM_use_pointer'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_use_pointer
  !============================================================================

  subroutine GM_synchronize_refinement(iProc0,iCommUnion)

    integer,intent(in) :: iProc0,iCommUnion

    character(len=*), parameter:: NameSub = 'GM_synchronize_refinement'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_synchronize_refinement
  !============================================================================

  subroutine GM_get_for_im(Buffer_IIV,BufferKp,iSize,jSize,nVar,NameVar)

    integer, intent(in) :: iSize,jSize,nVar
    real, intent(out)   :: BufferKp
    real, intent(out), dimension(iSize,jSize,nVar) :: Buffer_IIV
    character (len=*), intent(in) :: NameVar

    character(len=*), parameter:: NameSub = 'GM_get_for_im'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
  end subroutine GM_get_for_im
  !============================================================================

  subroutine GM_get_for_im_trace(nRadius, nLon, nVarLine, nPointLine, NameVar)

    ! Ray tracing for RAM type codes
    ! Provides total number of points along rays
    ! and the number of variables to pass to IM

    integer, intent(in)           :: nRadius, nLon
    integer, intent(out)          :: nVarLine, nPointLine
    character (len=*), intent(in) :: NameVar

    character(len=*), parameter:: NameSub = 'GM_get_for_im_trace'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_im_trace
  !============================================================================

  subroutine GM_get_for_im_line(nRadius, nLon, MapOut_DSII, &
       nVarLine, nPointLine, BufferLine_VI)

    integer, intent(in) :: nRadius, nLon
    real,    intent(out):: MapOut_DSII(3,2,nRadius,nLon)
    integer, intent(in) :: nPointLine, nVarLine
    real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)

    character(len=*), parameter:: NameSub = 'GM_get_for_im_line'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_im_line
  !============================================================================

  subroutine GM_put_from_im(Buffer_II,iSizeIn,jSizeIn,nVar,NameVar)

    integer, intent(in) :: iSizeIn,jSizeIn,nVar
    real, intent(in) :: Buffer_II(iSizeIn,jSizeIn)
    character(len=*), intent(in) :: NameVar

    character(len=*), parameter:: NameSub = 'GM_put_from_im'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
  end subroutine GM_put_from_im
  !============================================================================
  subroutine GM_put_from_im_cimi(Buffer_IIV,iSizeIn,jSizeIn,nVarIm,NameVarIm)

    integer, intent(in) :: iSizeIn,jSizeIn,nVarIm
    real, intent(in) :: Buffer_IIV(iSizeIn,jSizeIn,nVarIm)
    character(len=*), intent(in) :: NameVarIm

    character(len=*), parameter:: NameSub = 'GM_put_from_im_cimi'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
  end subroutine GM_put_from_im_cimi
  !============================================================================

  subroutine GM_satinit_for_im(nSats)

    integer, intent(out) :: nSats

    character(len=*), parameter:: NameSub = 'GM_satinit_for_im'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

  end subroutine GM_satinit_for_im
  !============================================================================

  subroutine GM_get_sat_for_im(Buffer_III, Buffer_I, nSats)

    integer, intent(in)               :: nSats
    real, intent(out)                 :: Buffer_III(3,2,nSats)
    character (len=100), intent(out)  :: Buffer_I(nSats)

    !--------------------------------------------------------------------------
  end subroutine GM_get_sat_for_im
  !============================================================================

  subroutine GM_get_sat_for_im_crcm(Buffer_III, Buffer_I, nSats)

    integer, intent(in)               :: nSats
    real, intent(out)                 :: Buffer_III(4,2,nSats)
    character (len=100), intent(out)  :: Buffer_I(nSats)

    !--------------------------------------------------------------------------
  end subroutine GM_get_sat_for_im_crcm
  !============================================================================

  subroutine GM_get_for_im_trace_crcm(iSizeIn, jSizeIn, nDensityIn, &
       NameVar, nVarLine, nPointLine)

    integer, intent(in)           :: iSizeIn, jSizeIn, nDensityIn
    character (len=*), intent(in) :: NameVar
    integer, intent(out)          :: nVarLine, nPointLine

    character(len=*), parameter:: NameSub = 'GM_get_for_im_trace_crcm'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_im_trace_crcm
  !============================================================================

  subroutine GM_get_for_im_crcm( &
       Buffer_IIV, KpOut,AeOut,iSizeIn, jSizeIn, nDensity,&
       nVarIn, BufferLine_VI, nVarLine, nPointLine, BufferSolarWind_V, NameVar)

    integer, intent(in) :: iSizeIn, jSizeIn, nDensity, nVarIn
    real, intent(out)   :: Buffer_IIV(iSizeIn,jSizeIn,nVarIn), KpOut, AeOut

    integer, intent(in) :: nPointLine, nVarLine
    real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)
    real,    intent(out):: BufferSolarWind_V(8)
    character(len=*), intent(in):: NameVar

    character(len=*), parameter:: NameSub = 'GM_get_for_im_crcm'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_im_crcm
  !============================================================================

  subroutine GM_put_from_ps(Buffer_II,iSizeIn,jSizeIn,nVar,NameVar)

    integer, intent(in) :: iSizeIn,jSizeIn,nVar
    real, intent(in) :: Buffer_II(iSizeIn,jSizeIn)
    character(len=*), intent(in) :: NameVar

    character(len=*), parameter:: NameSub = 'GM_put_from_ps'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
  end subroutine GM_put_from_ps
  !============================================================================

  subroutine GM_get_for_rb_trace(iSize,jSize,NameVar,nVarLine,nPointLine)

    integer, intent(in) :: iSize,jSize
    character (len=*), intent(in) :: NameVar
    integer, intent(out):: nVarLine,nPointLine

    character(len=*), parameter:: NameSub = 'GM_get_for_rb_trace'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
  end subroutine GM_get_for_rb_trace
  !============================================================================

  subroutine GM_get_for_rb(Buffer_IIV,iSize,jSize,nVar, &
       BufferLine_VI, nVarLine, nPointLine, NameVar)

    integer, intent(in) :: iSize,jSize,nVar
    real, intent(out)   :: Buffer_IIV(iSize,jSize,nVar)
    integer, intent(in) :: nVarLine, nPointLine
    real, intent(out)   :: BufferLine_VI(nVarLine, nPointLine)
    character (len=*), intent(in):: NameVar

    character(len=*), parameter:: NameSub = 'GM_get_for_rb'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
  end subroutine GM_get_for_rb
  !============================================================================

  subroutine GM_satinit_for_rb(nSats)
    integer :: nSats
    !--------------------------------------------------------------------------
  end subroutine GM_satinit_for_rb
  !============================================================================

  subroutine GM_get_sat_for_rb(Buffer_III, Buffer_I, nSats)

    integer, intent(in)               :: nSats
    real, intent(out)                 :: Buffer_III(4,2,nSats)
    character (len=100), intent(out)  :: Buffer_I(nSats)
    !--------------------------------------------------------------------------
  end subroutine GM_get_sat_for_rb
  !============================================================================

  subroutine GM_get_for_ie(Buffer_IIV,iSize,jSize,nVar)

    integer, intent(in) :: iSize,jSize,nVar
    real, intent(out), dimension(iSize,jSize,nVar) :: Buffer_IIV

    character(len=*), parameter:: NameSub = 'GM_get_for_ie'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_ie
  !============================================================================

  subroutine GM_get_info_for_ie(nVarIeGm, nVarGmIe, NameVar_I)

    integer, intent(out) :: nVarIeGM
    integer, intent(out), optional :: nVarGmIe
    character(len=*), intent(out), optional:: NameVar_I(:)

    character(len=*), parameter:: NameSub = 'GM_get_info_for_ie'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_info_for_ie
  !============================================================================

  subroutine GM_put_from_ie(Buffer_IIV, iSize, jSize, nVar, NameVar_I)

    integer,          intent(in):: iSize, jSize, nVar
    real,             intent(in):: Buffer_IIV(iSize,jSize,nVar)
    character(len=*), intent(in):: NameVar_I(nVar)

    character(len=*), parameter:: NameSub = 'GM_put_from_ie'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_put_from_ie
  !============================================================================

  subroutine GM_put_from_mh(nPartial,iPutStart,Put,Weight,DoAdd,StateSI_V,&
       nVar)
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::Weight
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::StateSI_V

    ! Derived type arguments, it is easier not to declare them

    character(len=*), parameter:: NameSub = 'GM_put_from_mh'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_put_from_mh
  !============================================================================

  subroutine GM_put_from_ih_buffer( &
       NameCoord, nY, nZ, yMin, yMax, zMin, zMax, Buffer_VII)

    character(len=*), intent(in) :: NameCoord
    integer,          intent(in) :: nY, nZ
    real,             intent(in) :: yMin, yMax, zMin, zMax
    real,             intent(in) :: Buffer_VII(8, nY, nZ)

    character(len=*), parameter:: NameSub = 'GM_put_from_ih_buffer'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_put_from_ih_buffer
  !============================================================================
  subroutine GM_get_for_global_buffer(&
       nR, nLon, nLat, BufferMinMax_DI, Buffer_VG)
    ! Buffer size and limits
    integer,intent(in) :: nR, nLon, nLat
    real, intent(in)   :: BufferMinMax_DI(3,2)

    ! OUTPUT ARGUMENTS
    ! State variables to be fiiled in all buffer grid points
    real,dimension(nVarCouple, nR, nLon, nLat), intent(out):: &
         Buffer_VG
    character(len=*), parameter:: NameSub = 'GM_get_for_global_buffer'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_global_buffer
  !============================================================================
  subroutine GM_put_from_pw(Buffer_VI, nVar, nFieldLine, Name_V)

    integer, intent(in)           :: nVar, nFieldLine
    real, intent(out)             :: Buffer_VI(nVar, nFieldLine)
    character (len=*), intent(in) :: Name_V(nVar)

    character(len=*), parameter:: NameSub = 'GM_put_from_pw'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_put_from_pw
  !============================================================================

  subroutine GM_get_for_pw(nTotalLine,p_I)

    integer, intent(in)           :: nTotalLine
    real, intent(out)             :: p_I(nTotalLine)

    character(len=*), parameter:: NameSub = 'GM_get_for_pw'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_pw
  !============================================================================

  subroutine GM_get_for_pt(IsNew, NameVar, nVarIn, nDimIn, nPoint, Pos_DI, &
       Data_VI)

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Pos_DI

    real, intent(in) :: Pos_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    character(len=*), parameter:: NameSub = 'GM_get_for_pt'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_pt
  !============================================================================

  subroutine GM_get_for_pc_init(nParamInt, nParamReal, iParam_I, Param_I)

    integer, intent(inout) :: nParamInt, nParamReal
    integer, optional, intent(out):: iParam_I(nParamInt)
    real,    optional, intent(out) :: Param_I(nParamReal)

    character(len=*), parameter:: NameSub = 'GM_get_for_pc_init'
    !--------------------------------------------------------------------------

    call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_pc_init
  !============================================================================

  subroutine GM_get_for_pc_grid_info(nInt, nPicGrid, AccumulatedSize_I, Int_I)
    integer, intent(inout) :: nInt, nPicGrid
    integer, optional, intent(out):: Int_I(nInt), AccumulatedSize_I(nPicGrid)

    character(len=*), parameter:: NameSub = 'GM_get_for_pc_grid_info'
    !--------------------------------------------------------------------------

    call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_pc_grid_info
  !============================================================================

  subroutine GM_get_for_pc_dt(DtSi)

    real, intent(out) ::  DtSi

    character(len=*), parameter:: NameSub = 'GM_get_for_pc_dt'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_pc_dt
  !============================================================================

  subroutine GM_get_for_pc(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    character(len=*), parameter:: NameSub = 'GM_get_for_pc'
    !--------------------------------------------------------------------------

    call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

  end subroutine GM_get_for_pc
  !============================================================================

  subroutine GM_put_from_pc( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI

    real,    intent(in), optional:: Data_VI(:,:)           ! Recv data array
    integer, intent(in), optional:: iPoint_I(nPoint)       ! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:) ! Position vectors

    character(len=*), parameter:: NameSub = 'GM_put_from_pc'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

  end subroutine GM_put_from_pc
  !============================================================================
  function GM_is_right_boundary_d(iBlock) RESULT(IsRightBoundary_D)
    integer, intent(in) :: iBlock
    logical :: IsRightBoundary_D(3)
    character(len=*), parameter:: NameSub = 'GM_is_right_boundary_d'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
  end function GM_is_right_boundary_d
  !============================================================================

  subroutine GM_put_from_ua(NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar
    integer,          intent(inout):: nVar
    integer,          intent(inout):: nPoint

    real,    intent(in), optional:: Data_VI(:,:)
    integer, intent(in), optional:: iPoint_I(nPoint)
    real, intent(out), allocatable, optional:: Pos_DI(:,:)

    character(len=*), parameter:: NameSub = 'GM_put_from_ua'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
    
  end subroutine GM_put_from_ua

end module GM_wrapper
!==============================================================================

! This subroutine is only needed because of SC|IH/BATSRUS/src/ModCellBoundary
subroutine read_ih_buffer(y, z, State_V)

  use ModUtilities, ONLY: CON_stop

  real, intent(in) :: y, z
  real, intent(out):: State_V(8)

  character(len=*), parameter:: NameSub = 'read_ih_buffer'
  !----------------------------------------------------------------------------
  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine read_ih_buffer
!==============================================================================

! This subroutine is only needed because of SC|IH/BATSRUS/src/ModFaceBoundary
subroutine read_pw_buffer(FaceCoords_D,nVar,FaceState_V)

  use ModUtilities, ONLY: CON_stop
  
  implicit none

  real, intent(in) :: FaceCoords_D(3)
  integer, intent(in) :: nVar
  real, intent(inout) :: FaceState_V(nVar)

  character(len=*), parameter:: NameSub = 'read_pw_buffer'
  !----------------------------------------------------------------------------
  call CON_stop(NameSub//'GM_ERROR: empty version cannot be used!')

end subroutine read_pw_buffer
!==============================================================================

