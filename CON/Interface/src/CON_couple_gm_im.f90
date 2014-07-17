!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE IM

!BOP
!MODULE: CON_couple_gm_im - couple GM and IM components
!
!DESCRIPTION:
! Couple GM and IM components both ways. 
!
!INTERFACE:
module CON_couple_gm_im

  !USES:
  use CON_coupler
  use CON_transfer_data, ONLY: transfer_integer, transfer_real_array, &
       transfer_string, transfer_string_array

  use GM_wrapper, ONLY: GM_get_for_im, GM_get_for_im_line, &
       GM_get_for_im_trace, GM_satinit_for_im, GM_get_sat_for_im, &
       GM_get_for_im_crcm, GM_get_for_im_trace_crcm, GM_get_sat_for_im_crcm, &
       GM_put_from_im, GM_print_variables

  use IM_wrapper, ONLY: IM_get_for_gm, IM_put_from_gm, &
       IM_put_from_gm_line, IM_put_from_gm_crcm, IM_put_sat_from_gm

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_im_init ! initialize both couplings
  public :: couple_gm_im      ! couple GM to IM
  public :: couple_im_gm      ! couple IM to GM

  !REVISION HISTORY:
  ! 07/25/2003 G.Toth <gtoth@umich.edu> - initial version
  !            O.Volberg and D.DeZeeuw
  !
  ! 08/27/2003 G.Toth - external subroutines combined into a module
  ! 01/01/2007 D.Welling - added satellite info tranfer
  !EOP

  logical :: IsInitialized = .false.

  ! Size of the 2D spherical structured (possibly non-uniform) IM grid
  integer, save :: iSize, jSize

  ! Number of satellites in GM that will also be traced in IM
  integer, save :: nShareSats

  logical, save :: DoMultiFluidIMCoupling, DoAnisoPressureIMCoupling

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_im_init - initialize GM-IM coupling
  !INTERFACE:
  subroutine couple_gm_im_init
    !DESCRIPTION:
    ! Store IM grid size.
    !EOP
    
    use ModProcessVarName,  ONLY: process_var_name

    integer :: nDensityGm, nSpeedGm, nPGm, nPparGm, nWaveGm, nMaterialGm
    integer :: nDensityIm, nSpeedIm, nPIm, nPparIm, nWaveIm, nMaterialIm 

    ! General error code
    !------------------------------------------------------------------------

    if(IsInitialized) RETURN
    IsInitialized = .true.

    if(.not. (is_proc(IM_) .or. is_proc(GM_))) RETURN

    ! This works for a regular IM grid only
    iSize = Grid_C(IM_) % nCoord_D(1)
    jSize = Grid_C(IM_) % nCoord_D(2)

    call process_var_name(Grid_C(GM_)%NameVar, nDensityGm, nSpeedGm, &
         nPGm, nPparGm, nWaveGm, nMaterialGm)
    call process_var_name(Grid_C(IM_)%NameVar, nDensityIm, nSpeedIm, &
         nPIm, nPparIm, nWaveIm, nMaterialIm)

    DoMultiFluidIMCoupling = nDensityGm > 1 .and. nDensityIm > 1
    
    DoAnisoPressureIMCoupling = nPparGm > 0 .and. nPparIm > 0 

    ! Set number of satellites shared between GM and IM for tracing.
    if(is_proc(GM_)) call GM_satinit_for_im(nShareSats)
    call transfer_integer(GM_, IM_, nShareSats, UseSourceRootOnly=.false.)
    
  end subroutine couple_gm_im_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_im - couple GM to IM component
  !INTERFACE:
  subroutine couple_gm_im(tSimulation)

    use CON_world, ONLY: get_comp_info
    use CON_comp_param, ONLY: lNameVersion

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere (GM) source\\
    !    Inner Magnetosphere  (IM) target
    !
    ! Send field line volumes, average density and pressure and
    ! geometrical information.
    !EOP

    ! Which IM model is used?
    character(len=lNameVersion):: NameVersionIm

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='couple_gm_im'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    call get_comp_info(IM_,NameVersion=NameVersionIm)
    select case(NameVersionIm(1:3))
    case('RCM')
       call couple_rcm
    case('RAM')
       call couple_ram
    case('CRC')
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

      ! Buffer for satellite names
      character(len=100), allocatable:: NameSat_I(:)

      character(len=*), parameter :: NameSubSub=NameSub//'.couple_rcm'
      !------------------------------------------------------------------------
      if(DoTest)write(*,*)NameSubSub,' starting, iProc=', i_proc()

      if(DoMultiFluidIMCoupling) then
         NameVar='vol:z0x:z0y:bmin:rho:p:Hprho:Oprho:Hpp:Opp'
         nVarGmIm = 10
      else
         NameVar='vol:z0x:z0y:bmin:rho:p'
         nVarGmIm = 6
      endif

      !\
      ! Get field line integrals from GM to IM
      !/
      allocate(Buffer_IIV(iSize,jSize,nVarGmIm))
      ! All GM processors participate in calculating the integrals
      if(is_proc(GM_)) &
           call GM_get_for_im(Buffer_IIV, iSize, jSize, nVarGmIm, NameVar)
      ! The integrals are available on GM root only
      call transfer_real_array(GM_, IM_, size(Buffer_IIV), Buffer_IIV)
      if(is_proc0(IM_)) &
           call IM_put_from_gm(Buffer_IIV, iSize, jSize, nVarGmIm, NameVar)
      deallocate(Buffer_IIV)

      !\
      ! If IM sat tracing is enabled, get sat locations from GM
      !/
      if (nShareSats > 0) then
         allocate(SatPos_DII(3,2,nShareSats), NameSat_I(nShareSats))
         if(is_proc(GM_)) &
              call GM_get_sat_for_im(SatPos_DII, NameSat_I, nShareSats)
         call transfer_string_array(GM_, IM_, nShareSats, NameSat_I, &
              UseSourceRootOnly = .false.)
         call transfer_real_array(GM_, IM_, size(SatPos_DII), SatPos_DII)
         if(is_proc(IM_)) &
              call IM_put_sat_from_gm(nShareSats, NameSat_I, SatPos_DII)
         deallocate(SatPos_DII, NameSat_I)
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

      !DESCRIPTION:
      ! Couple between two components:\\
      !    Global Magnetosphere (GM) source\\
      !    Inner  Magnetosphere (IM) target
      !EOP

      !\
      ! Coupling variables
      !/

      ! Number of varibles at minimum B to pass
      integer  :: nVarBmin

      ! Number of variables saved into the line data
      integer:: nVarLine

      ! Number of points saved into the line data      
      integer:: nPointLine

      ! Names of variables to pass
      character (len=100) :: NameVar

      ! Buffer for the variables on the 2D IM grid and line data
      real, allocatable:: Buffer_IIV(:,:,:), BufferLine_VI(:,:)

      ! Buffer for satellite locations   
      real, allocatable :: SatPos_DII(:,:,:)

      ! Buffer for satellite names   
      character(len=100), allocatable:: NameSat_I(:)

      logical :: DoTest, DoTestMe
      character(len=*), parameter :: NameSub='couple_gm_im_crcm'
      !------------------------------------------------------------------------
      call CON_set_do_test(NameSub, DoTest, DoTestMe)

      if(DoTest)write(*,*)NameSub,' starting, iProc=', i_proc()

      if(DoMultiFluidIMCoupling) then
         NameVar='x:y:bmin:I_I:S_I:R_I:B_I:rho:p:Hprho:Oprho:Hpp:Opp'
         nVarBmin = 10
      else if(DoAnisoPressureIMCoupling)then
         NameVar='x:y:bmin:I_I:S_I:R_I:B_I:rho:p:ppar'
         nVarBmin = 7
      else
         NameVar='x:y:bmin:I_I:S_I:R_I:B_I:rho:p'
         nVarBmin = 6
      endif

      !\
      ! Get size of field line traces
      !/
      if(is_proc(GM_)) call GM_get_for_im_trace_crcm( &
           iSize, jSize, NameVar, nVarLine, nPointLine)

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
           Buffer_IIV, iSize, jSize, nVarBmin, &
           BufferLine_VI, nVarLine, nPointLine, NameVar)

      call transfer_real_array(GM_, IM_, size(Buffer_IIV), Buffer_IIV)
      call transfer_real_array(GM_, IM_, size(BufferLine_VI), BufferLine_VI)

      if(is_proc(IM_)) call IM_put_from_gm_crcm(&
           Buffer_IIV, iSize, jSize, nVarBmin,&
           BufferLine_VI, nVarLine, nPointLine, NameVar, tSimulation)

      deallocate(Buffer_IIV, BufferLine_VI)

      ! If IM sat tracing is enabled, get sat locations from GM
      if (nShareSats > 0) then
         allocate(SatPos_DII(4,2,nShareSats), NameSat_I(nShareSats))
         if(is_proc(GM_)) &
              call GM_get_sat_for_im_crcm(SatPos_DII, NameSat_I, nShareSats)
         
         ! Transfer satellite names from GM to IM
         call transfer_string_array(GM_, IM_, nShareSats, NameSat_I, &
              UseSourceRootOnly = .false.)
         ! Transfer satellite locations from GM to IM
         call transfer_real_array(GM_, IM_, size(SatPos_DII), SatPos_DII)

         if(is_proc(IM_)) &
              call IM_put_sat_from_gm(nShareSats, NameSat_I, SatPos_DII)
         deallocate(SatPos_DII, NameSat_I)
      end if

    end subroutine couple_crcm

  end subroutine couple_gm_im

  !BOP =======================================================================
  !IROUTINE: couple_im_gm - couple IM to GM component
  !INTERFACE:
  subroutine couple_im_gm(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Inner Magnetosphere  (IM) source\\
    !    Global Magnetosphere (GM) target
    !
    ! Send pressure from IM to GM.
    !EOP

    ! Number of variables to pass
    integer:: nVarImGm

    ! Names of variables to pass
    character(len=100) :: NameVar

    ! Buffer for the variables on the 2D IM grid
    real, allocatable:: Buffer_IIV(:,:,:)

    logical:: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='couple_im_gm'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting, iProc=',i_proc()

    if(DoMultiFluidIMCoupling)then
       NameVar='p:rho:Hpp:Opp:Hprho:Oprho'
       nVarImGm=6
    else if(DoAnisoPressureIMCoupling)then 
       NameVar='p:rho:ppar:bmin'
       nVarImGm=4
    else
       NameVar='p:rho'
       nVarImGm=2
    end if

    allocate(Buffer_IIV(iSize,jSize,nVarImGm))
    if(is_proc(IM_)) &
         call IM_get_for_gm(Buffer_IIV, iSize, jSize, nVarImGm, NameVar)
    call transfer_real_array(IM_, GM_, size(Buffer_IIV), Buffer_IIV)
    if(is_proc(GM_)) &
         call GM_put_from_im(Buffer_IIV, iSize, jSize, nVarImGm, NameVar)
    deallocate(Buffer_IIV)

    if(DoTest)write(*,*)NameSub,': finished iProc=', i_proc()
    if(DoTest.and.is_proc0(GM_)) call GM_print_variables('IM')

  end subroutine couple_im_gm

end module CON_couple_gm_im
