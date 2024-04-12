!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE IM

module CON_couple_gm_im

  ! Couple GM and IM components both ways.

  use CON_coupler
  use CON_transfer_data, ONLY: transfer_integer, transfer_real_array, &
       transfer_real, transfer_string, transfer_string_array

  use GM_wrapper, ONLY: GM_get_for_im, GM_get_for_im_line, &
       GM_get_for_im_trace, GM_satinit_for_im, GM_get_sat_for_im, &
       GM_get_for_im_crcm, GM_get_for_im_trace_crcm, GM_get_sat_for_im_crcm, &
       GM_put_from_im, GM_put_from_im_cimi

  use IM_wrapper, ONLY: IM_get_for_gm, IM_put_from_gm, &
       IM_put_from_gm_line, IM_put_from_gm_crcm, IM_put_sat_from_gm

  implicit none

  private ! except

  public :: couple_gm_im_init ! initialize both couplings
  public :: couple_gm_im      ! couple GM to IM
  public :: couple_im_gm      ! couple IM to GM

  ! revision history:
  ! 07/25/2003 G.Toth <gtoth@umich.edu> - initial version
  !            O.Volberg and D.DeZeeuw
  !
  ! 08/27/2003 G.Toth - external subroutines combined into a module
  ! 01/01/2007 D.Welling - added satellite info tranfer

  logical :: IsInitialized = .false.

  ! Size of the 2D spherical structured (possibly non-uniform) IM grid
  integer:: iSize = 0, jSize = 0

  ! Number of satellites in GM that will also be traced in IM
  integer:: nShareSats = 0

  ! Number of coupled density and pressure variables
  integer :: nRhoPCoupled = 0

  logical, save :: DoMultiFluidIMCoupling, DoAnisoPressureIMCoupling

  ! Number of fluids in GM (useful as multifluid can have more than 2 fluids)
  integer:: nDensityGM = 0

contains
  !============================================================================
  subroutine couple_gm_im_init

    ! Store IM grid size.

    use ModProcessVarName,  ONLY: process_var_name

    integer :: nSpeedGm, nPGm, nPparGm, nWaveGm, nMaterialGm, nChargeStateAllGm
    integer :: nDensityIm, nSpeedIm, nPIm, nPparIm, nWaveIm, nMaterialIm, &
         nChargeStateAllIm

    ! General error code
    !--------------------------------------------------------------------------

    if(IsInitialized) RETURN
    IsInitialized = .true.

    if(.not. (is_proc(IM_) .or. is_proc(GM_))) RETURN

    ! This works for a regular IM grid only
    iSize = Grid_C(IM_) % nCoord_D(1)
    jSize = Grid_C(IM_) % nCoord_D(2)

    call set_couple_var_info(GM_,IM_)
    nRhoPCoupled = nVarBuffer
    ! if(is_proc0(GM_))then
    !   write(*,*)'(GM_)%NameVar    =', Grid_C(GM_)%NameVar
    !   write(*,*)'(IM_)%NameVar    =', Grid_C(IM_)%NameVar
    !   write(*,*)'nVarBuffer       =', nVarBuffer
    !   write(*,*)'iVarSource_V     =', iVarSource_V(1:12)
    !   write(*,*)'iVarTarget_V     =', iVarTarget_V(1:12)
    !   !    call con_stop('')
    ! end if

    ! this will likely be removed when coupling generalization if done
    call process_var_name(Grid_C(GM_)%NameVar, nDensityGm, nSpeedGm, &
         nPGm, nPparGm, nWaveGm, nMaterialGm, nChargeStateAllGm)
    call process_var_name(Grid_C(IM_)%NameVar, nDensityIm, nSpeedIm, &
         nPIm, nPparIm, nWaveIm, nMaterialIm, nChargeStateAllIm)

    DoMultiFluidIMCoupling = nDensityGm > 1 .and. nDensityIm > 1

    DoAnisoPressureIMCoupling = nPparGm > 0 .and. nPparIm > 0

    ! Set number of satellites shared between GM and IM for tracing.
    if(is_proc(GM_)) call GM_satinit_for_im(nShareSats)
    call transfer_integer(GM_, IM_, nShareSats, UseSourceRootOnly=.false.)

  end subroutine couple_gm_im_init
  !============================================================================
  subroutine couple_gm_im(tSimulation)

    use CON_world, ONLY: get_comp_info
    use CON_comp_param, ONLY: lNameVersion

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Couple between two components:
    !    Global Magnetosphere (GM) source
    !    Inner Magnetosphere  (IM) target
    !
    ! Send field line volumes, average density and pressure and
    ! geometrical information.

    ! Which IM model is used?
    character(len=lNameVersion):: NameVersionIm

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_gm_im'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    call get_comp_info(IM_,NameVersion=NameVersionIm)
    select case(NameVersionIm(1:3))
    case('RCM')
       call couple_rcm
    case('RAM')
       call couple_ram
    case('CRC', 'CIM')
       call couple_crcm
    case default
       call CON_stop(NameSub//': unknown IM version='//NameVersionIm)
    end select

    if(DoTest)write(*,*)NameSub,' finished, iProc=', i_proc()

  contains
    !==========================================================================
    subroutine couple_rcm

      ! Number of variables to pass
      integer :: nVarGmIm

      character(len=100) :: NameVar

      ! Buffer for the variables on the 2D IM grid
      real, allocatable :: Buffer_IIV(:,:,:)

      ! Buffer for satellite locations
      real, allocatable :: SatPos_DII(:,:,:)

      ! Buffer for scalar Kp
      real :: BufferKp

      ! Buffer for satellite names
      character(len=100), allocatable:: NameFileSat_I(:)

      character(len=*), parameter :: NameSubSub=NameSub//'.couple_rcm'
      !------------------------------------------------------------------------
      if(DoTest)write(*,*)NameSubSub,' starting, iProc=', i_proc()

      if(DoMultiFluidIMCoupling) then
         NameVar='vol:z0x:z0y:bmin:Hprho:Oprho:Hpp:Opp:pe'
         nVarGmIm = 9
      else
         NameVar='vol:z0x:z0y:bmin:rho:p:pe'
         nVarGmIm = 7
      endif

      ! Get field line integrals from GM to IM
      allocate(Buffer_IIV(iSize,jSize,nVarGmIm))
      ! All GM processors participate in calculating the integrals
      if(is_proc(GM_)) &
           call GM_get_for_im(Buffer_IIV, BufferKp, &
           iSize, jSize, nVarGmIm, NameVar)
      ! The integrals are available on GM root only
      call transfer_real_array(GM_, IM_, size(Buffer_IIV), Buffer_IIV)
      call transfer_real(GM_,IM_, BufferKp)
      if(is_proc0(IM_)) &
           call IM_put_from_gm(Buffer_IIV, BufferKp, &
           iSize, jSize, nVarGmIm, NameVar)
      deallocate(Buffer_IIV)

      ! If IM sat tracing is enabled, get sat locations from GM
      if (nShareSats > 0) then
         allocate(SatPos_DII(3,2,nShareSats), NameFileSat_I(nShareSats))
         if(is_proc(GM_)) then
            call GM_get_sat_for_im(SatPos_DII, NameFileSat_I, nShareSats)
         end if
         call transfer_string_array(GM_, IM_, nShareSats, NameFileSat_I, &
              UseSourceRootOnly = .false.)
         call transfer_real_array(GM_, IM_, size(SatPos_DII), SatPos_DII)
         if(is_proc(IM_)) then
            call IM_put_sat_from_gm(nShareSats, NameFileSat_I, SatPos_DII)
         end if
         deallocate(SatPos_DII, NameFileSat_I)
      end if

    end subroutine couple_rcm
    !==========================================================================
    subroutine couple_ram

      ! Some variables do not change during a run, so transfer once only
      logical:: IsFirstTime = .true.

      ! List of variable names
      character(len=100), save :: NameVar

      ! Number of variables saved into the line data
      integer:: nVarLine

      ! Number of points saved into the line data
      integer:: nPointLine

      ! Buffers for the line data and mapping from equator to ionosphere
      real, allocatable:: BufferLine_VI(:,:), Map_DSII(:,:,:,:)

      character(len=*), parameter :: NameSubSub=NameSub//'.couple_ram'
      !------------------------------------------------------------------------
      if(DoTest)write(*,*)NameSubSub,' starting, iProc=', i_proc()

      ! Get field line trace sizes from GM and transfer it to IM
      if(is_proc(GM_)) call GM_get_for_im_trace( &
           iSize, jSize, nVarLine, nPointLine, NameVar)

      if(DoTest .and. is_proc0(GM_)) write(*,*) NameSubSub, &
           ' GM root has nVarLine, nPointLine=', nVarLine, nPointLine

      if(IsFirstTime)then
         ! This cannot change during a run, so only transfer them once
         call transfer_string(GM_, IM_, NameVar)
         IsFirstTime = .false.
      end if
      call transfer_integer(GM_, IM_, nVarLine, nPointLine)

      if(DoTest .and. is_proc0(IM_)) write(*,*) NameSubSub, &
           ' IM root has nVarLine, nPointLine=', nVarLine, nPointLine

      ! Now the size is known on GM root and IM so allocate arrays
      if(is_proc0(GM_) .or. is_proc(IM_))then
         allocate(Map_DSII(3,2,iSize,jSize), &
              BufferLine_VI(nVarLine,nPointLine))
      else
         allocate(Map_DSII(1,1,1,1), BufferLine_VI(1,1))
      end if

      ! Transfer the trace data and mapping data from GM root to IM
      if(is_proc0(GM_)) call GM_get_for_im_line( &
           iSize, jSize, Map_DSII, nVarLine, nPointLine, BufferLine_VI)
      call transfer_real_array(GM_, IM_, size(Map_DSII), Map_DSII)
      call transfer_real_array(GM_, IM_, size(BufferLine_VI), BufferLine_VI)

      if(is_proc(IM_)) call IM_put_from_gm_line(iSize, jSize, Map_DSII, &
           nVarLine, nPointLine, BufferLine_VI, NameVar)
      deallocate(BufferLine_VI, Map_DSII)

    end subroutine couple_ram
    !==========================================================================
    subroutine couple_crcm

      ! Couple between two components:
      !    Global Magnetosphere (GM) source
      !    Inner  Magnetosphere (IM) target

      ! Coupling variables

      ! Number of varibles at minimum B to pass
      integer  :: nVarBmin

      ! Number of variables saved into the line data
      integer:: nVarLine

      ! Number of points saved into the line data
      integer:: nPointLine

      ! Names of variables to pass
      character (len=100) :: NameVar

      ! Buffer for the variables on the 2D IM grid and line data
      real, allocatable :: Buffer_IIV(:,:,:), BufferLine_VI(:,:)

      ! Buffer for passing the solar wind variables
      real :: BufferSolarWind_V(8)

      ! Buffer for satellite locations
      real, allocatable :: SatPos_DII(:,:,:)

      ! Buffer for scalar Kp
      real :: BufferKp, BufferAe

      ! Buffer for satellite names
      character(len=100), allocatable:: NameFileSat_I(:)

      logical :: DoTest, DoTestMe

      character(len=*), parameter:: NameSub = 'couple_crcm'
      !------------------------------------------------------------------------
      call CON_set_do_test(NameSub, DoTest, DoTestMe)

      if(DoTest)write(*,*)NameSub,' starting, iProc=', i_proc()

      ! Number of variables in Buffer_IIV indexed by Lon, Lat
      ! x, y, B, densities and pressures at the min B surface
      nVarBmin = nRhoPCoupled + 5

      NameVar = 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p'

      ! Get size of field line traces
      if(is_proc(GM_)) call GM_get_for_im_trace_crcm( &
           iSize, jSize, nDensityGM, NameVar, nVarLine, nPointLine)

      call transfer_integer(GM_, IM_, nVarLine, nPointLine)

      ! Now the size is known on GM root and IM so allocate arrays
      if(is_proc0(GM_) .or. is_proc(IM_))then
         allocate(BufferLine_VI(nVarLine, nPointLine), &
              Buffer_IIV(iSize,jSize,nVarBmin))
      else
         allocate(BufferLine_VI(1,1), Buffer_IIV(1,1,1))
      end if

      ! Only GM root returns useful info but all processors should be called
      ! so they can deallocate ray tracing
      if(is_proc0(GM_)) call GM_get_for_im_crcm( &
           Buffer_IIV, BufferKp, BufferAe, iSize, jSize, nDensityGM, nVarBmin, &
           BufferLine_VI, nVarLine, nPointLine, BufferSolarWind_V, NameVar)

      call transfer_real_array(GM_, IM_, size(Buffer_IIV), Buffer_IIV)
      call transfer_real_array(GM_, IM_, size(BufferLine_VI), BufferLine_VI)
      call transfer_real_array(GM_, IM_, size(BufferSolarWind_V), BufferSolarWind_V)
      call transfer_real(GM_,IM_, BufferKp)
      call transfer_real(GM_,IM_, BufferAe)

      if(is_proc(IM_)) call IM_put_from_gm_crcm(&
           Buffer_IIV, BufferKp, BufferAe, iSize, jSize, nVarBmin, &
           BufferLine_VI, nVarLine, nPointLine, NameVar, &
           BufferSolarWind_V, tSimulation)

      deallocate(Buffer_IIV, BufferLine_VI)

      ! If IM sat tracing is enabled, get sat locations from GM
      if (nShareSats > 0) then
         allocate(SatPos_DII(4,2,nShareSats), NameFileSat_I(nShareSats))
         if(is_proc(GM_)) &
              call GM_get_sat_for_im_crcm(SatPos_DII, NameFileSat_I, nShareSats)

         ! Transfer satellite names from GM to IM
         call transfer_string_array(GM_, IM_, nShareSats, NameFileSat_I, &
              UseSourceRootOnly = .false.)
         ! Transfer satellite locations from GM to IM
         call transfer_real_array(GM_, IM_, size(SatPos_DII), SatPos_DII)

         if(is_proc(IM_)) &
              call IM_put_sat_from_gm(nShareSats, NameFileSat_I, SatPos_DII)
         deallocate(SatPos_DII, NameFileSat_I)
      end if

    end subroutine couple_crcm
    !==========================================================================
  end subroutine couple_gm_im
  !============================================================================
  subroutine couple_im_gm(tSimulation)

    use CON_world, ONLY: get_comp_info
    use CON_comp_param, ONLY: lNameVersion

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Couple between two components:
    !    Inner Magnetosphere  (IM) source
    !    Global Magnetosphere (GM) target
    !
    ! Send pressure from IM to GM.

    ! Which IM model is used?
    character(len=lNameVersion):: NameVersionIm

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_im_gm'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    call get_comp_info(IM_,NameVersion=NameVersionIm)

    select case(NameVersionIm(1:3))
    case('RCM', 'RAM','CRC')
       call couple_im_gm_default
    case('CIM')
       call couple_im_gm_cimi
    case default
       call CON_stop(NameSub//': unknown IM version='//NameVersionIm)
    end select

    if(DoTest)write(*,*)NameSub,' finished, iProc=', i_proc()

  contains
    !==========================================================================
    subroutine couple_im_gm_default

      ! Number of variables to pass
      integer:: nVarImGm

      ! Names of variables to pass
      character(len=100) :: NameVar

      ! Buffer for the variables on the 2D IM grid
      real, allocatable:: Buffer_IIV(:,:,:)

      logical:: DoTest, DoTestMe

      character(len=*), parameter:: NameSub = 'couple_im_gm_default'
      !------------------------------------------------------------------------
      call CON_set_do_test(NameSub,DoTest,DoTestMe)

      if(DoTest)write(*,*)NameSub,' starting, iProc=',i_proc()

      if(DoMultiFluidIMCoupling)then
         NameVar='pe:p:rho:Hpp:Opp:Hprho:Oprho'
         nVarImGm=7
      else if(DoAnisoPressureIMCoupling)then
         NameVar='p:rho:ppar:bmin'
         nVarImGm=4
      else
         NameVar='pe:p:rho'
         nVarImGm=3
      end if

      allocate(Buffer_IIV(iSize,jSize,nVarImGm))
      if(is_proc(IM_)) &
           call IM_get_for_gm(Buffer_IIV, iSize, jSize, nVarImGm, NameVar)
      call transfer_real_array(IM_, GM_, size(Buffer_IIV), Buffer_IIV)
      if(is_proc(GM_)) &
           call GM_put_from_im(Buffer_IIV, iSize, jSize, nVarImGm, NameVar)
      deallocate(Buffer_IIV)

      if(DoTest)write(*,*)NameSub,': finished iProc=', i_proc()

    end subroutine couple_im_gm_default
    !==========================================================================
    subroutine couple_im_gm_cimi

      ! Number of variables to pass
      integer:: nVarImGm

      ! Names of variables to pass
      character(len=100) :: NameVar

      ! Buffer for the variables on the 2D IM grid
      real, allocatable:: Buffer_IIV(:,:,:)

      logical:: DoTest, DoTestMe

      character(len=*), parameter:: NameSub = 'couple_im_gm_cimi'
      !------------------------------------------------------------------------
      call CON_set_do_test(NameSub,DoTest,DoTestMe)

      if(DoTest)write(*,*)NameSub,' starting, iProc=',i_proc()

      ! get number of variable coupled from IM to GM from the grid descriptor
      ! add 1 since we need also the minimum B which we tack onto the end
      nVarImGm = Grid_C(IM_)%nVar+1

      ! for now pass in the IM namevar although this is not really needed
      NameVar = Grid_C(IM_)%NameVar

      allocate(Buffer_IIV(iSize,jSize,nVarImGm))
      if(is_proc0(IM_)) &
           call IM_get_for_gm(Buffer_IIV, iSize, jSize, nVarImGm, NameVar)
      call transfer_real_array(IM_, GM_, size(Buffer_IIV), Buffer_IIV)
      if(is_proc(GM_)) &
           call GM_put_from_im_cimi(Buffer_IIV, iSize, jSize, nVarImGm, &
           NameVar)
      deallocate(Buffer_IIV)

      if(DoTest)write(*,*)NameSub,': finished iProc=', i_proc()
    end subroutine couple_im_gm_cimi
    !==========================================================================
  end subroutine couple_im_gm
  !============================================================================
end module CON_couple_gm_im
!==============================================================================
