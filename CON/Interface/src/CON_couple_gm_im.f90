!^CMP COPYRIGHT UM
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

  ! Communicator and logicals to simplify message passing and execution
  integer, save :: iProc0Im, iProc0Gm, iCommWorld

  logical :: IsInitialized = .false., UseMe=.true.

  ! Size of the 2D spherical structured (possibly non-uniform) IM grid
  integer, save :: iSize, jSize, nCells_D(2)

  ! Number of satellites in GM that will also be traced in IM
  integer, save :: nShareSats

  logical, save :: DoMultiFluidIMCoupling

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_im_init - initialize GM-IM coupling
  !INTERFACE:
  subroutine couple_gm_im_init
    !DESCRIPTION:
    ! Store IM grid size.
    !EOP
    !------------------------------------------------------------------------

    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! General error code
    integer :: iError
    integer :: nFluid

    if(IsInitialized) RETURN
    IsInitialized = .true.

    UseMe = is_proc(IM_) .or. is_proc(GM_)
    if(.not.UseMe) RETURN

    ! Store these for root-to-root communication
    iProc0Im   = i_proc0(IM_)
    iProc0Gm   = i_proc0(GM_)
    iCommWorld = i_comm()

    ! This works for a regular IM grid only
    nCells_D = ncells_decomposition_d(IM_)
    iSize = nCells_D(1); jSize = nCells_D(2)

    ! Set number of satellites shared between GM and IM for tracing.
    call GM_satinit_for_im(nShareSats)

    ! Send number of satellites GM to IM
    if(iProc0Im /= iProc0Gm)then
       if(is_proc0(GM_)) &
            call MPI_send(nShareSats,1,MPI_INTEGER,iProc0Im,&
            1,iCommWorld,iError)
       if(is_proc0(IM_)) &
            call MPI_recv(nShareSats,1,MPI_INTEGER,iProc0Gm,&
            1,iCommWorld,iStatus_I,iError)
    end if

    ! Get the logical variable between GM and IM coupling: if multifluid coupling
    call GM_get_multi_for_im(DoMultiFluidIMCoupling)

    ! Send DoMultiFluidIMCoupling from GM to IM
    if(iProc0Im /= iProc0Gm)then
       if(is_proc0(GM_)) then
          if (DoMultiFluidIMCoupling) then
             call MPI_send(2,1,MPI_INTEGER,iProc0Im,&
                  2,iCommWorld,iError)
          else
             call MPI_send(1,1,MPI_INTEGER,iProc0Im,&
                  2,iCommWorld,iError)
          end if
       end if
       if(is_proc0(IM_)) then
          call MPI_recv(nFluid,1,MPI_INTEGER,iProc0Gm,&
               2,iCommWorld,iStatus_I,iError)
          if (nFluid == 2) then
             DoMultiFluidIMCoupling = .true.
          else
             DoMultiFluidIMCoupling = .false.
          endif
       end if
    end if


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

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_im'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld
    character(len=lNameVersion):: NameVersionIm
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()

    call get_comp_info(IM_,NameVersion=NameVersionIm)
    select case(NameVersionIm(1:3))
    case('RCM')
       call couple_rcm
    case('RAM')
       call couple_ram
    case('CRC')
       call couple_crcm
    end select

    if(DoTest)write(*,*)NameSub,': finished iProc=',iProcWorld

  contains

    !==========================================================================
    subroutine couple_rcm

      character (len=*), parameter :: NameSubSub=NameSub//'.couple_rcm'

      ! Number of variables to pass
      integer :: nVarGmIm
 
      character (len=100) :: NameVar

      ! Buffer for the variables on the 2D IM grid
      real, dimension(:,:,:), allocatable :: Buffer_IIV

      ! Buffer for satellite locations
      real, dimension(:,:,:), allocatable :: SatPos_DII

      ! Buffer for satellite names
      character (len=100), dimension(:), allocatable:: NameSat_I

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, GMi_iProc0, iProc0Im=', &
           iProcWorld,iProc0Gm,iProc0Im

       if(DoMultiFluidIMCoupling) then
          NameVar='vol:z0x:z0y:bmin:rho:p:Hprho:Oprho:Hpp:Opp'
          nVarGmIm = 10
       else
          NameVar='vol:z0x:z0y:bmin:rho:p'
          nVarGmIm = 6
       endif

      !\
      ! Allocate buffers both in GM and IM
      !/

      allocate(Buffer_IIV(iSize,jSize,nVarGmIm), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_IIV")
      if (nShareSats > 0) then
         allocate(SatPos_DII(3,2,nShareSats), stat=iError)
         call check_allocate(iError,NameSubSub//": SatPos_DII")
         allocate(NameSat_I(nShareSats),       stat=iError)
         call check_allocate(iError,NameSubSub//": NameSat_I")
      end if

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld

      !\
      ! Get field line integrals from GM
      !/
      if(is_proc(GM_)) &           
           call GM_get_for_im(Buffer_IIV,iSize,jSize,nVarGmIm,NameVar)

      !\
      ! If IM sat tracing is enabled, get sat locations from GM
      !/
      if(is_proc(GM_).AND.(nShareSats > 0)) &
           call GM_get_sat_for_im(SatPos_DII, NameSat_I, nShareSats)

      !\
      ! Transfer physical variables from GM to IM
      !/
      if(iProc0Im /= iProc0Gm)then
         nSize = iSize*jSize*nVarGmIm
         if(is_proc0(GM_)) then 
            call MPI_send(Buffer_IIV,nSize,MPI_REAL,iProc0Im,&
                 1,iCommWorld,iError)
         endif
         if(is_proc0(IM_)) then
            call MPI_recv(Buffer_IIV,nSize,MPI_REAL,iProc0Gm,&
                 1,iCommWorld,iStatus_I,iError)
         endif
      end if

      !\
      ! Transfer satellite names from GM to IM
      !/   

      if(nShareSats > 0 .and. iProc0Im /= iProc0Gm)then
         nSize = nShareSats*100
         if(is_proc0(GM_)) then
            call MPI_send(NameSat_I,nSize,MPI_BYTE,iProc0Im,&
                 1,iCommWorld,iError)
         endif
         if(is_proc0(IM_)) then
            call MPI_recv(NameSat_I,nSize,MPI_BYTE,iProc0Gm,&
                 1,iCommWorld,iStatus_I,iError)
         endif

         ! Transfer satellite locations from GM to IM

         nSize = 3*2*nShareSats
         if(is_proc0(GM_)) then
            call MPI_send(SatPos_DII,nSize,MPI_REAL,iProc0Im,&
                 1,iCommWorld,iError)
         endif
         if(is_proc0(IM_)) then
            call MPI_recv(SatPos_DII,nSize,MPI_REAL,iProc0Gm,&
                 1,iCommWorld,iStatus_I,iError)
         endif
      end if

      if(DoTest)write(*,*)NameSubSub,', variables transferred',&
           ', iProc:',iProcWorld

      !\
      ! Put variables into IM
      !/
      if(is_proc0(IM_))then
         call IM_put_from_gm(Buffer_IIV,iSize,jSize,nVarGmIm,NameVar)
         if(nShareSats > 0) &
              call IM_put_sat_from_gm(nShareSats, NameSat_I, SatPos_DII)
         if(DoTest) &
              write(*,*)'get_fieldline_volume: IM iProc, Buffer(1,1)=',&
              iProcWorld,Buffer_IIV(1,1,:)
      end if

      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_IIV)

      if (nShareSats > 0) then
         deallocate(NameSat_I)
         deallocate(SatPos_DII)
      end if

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_rcm

    !==========================================================================
    subroutine couple_ram

      character (len=*), parameter :: NameSubSub=NameSub//'.couple_ram'

      ! Number of variables and points saved into the line data
      integer :: nVarLine, nPointLine

      ! Buffer for the line data
      real, allocatable :: BufferLine_VI(:,:), Map_DSII(:,:,:,:)
 
      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      ! List of variables
      integer, parameter :: lNameVar=100
      character (len=lNameVar) :: NameVar
      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, GMi_iProc0, iProc0Im=', &
           iProcWorld,iProc0Gm,iProc0Im

      !\
      ! Get field line info from GM
      !/
      if(is_proc(GM_)) then
         call GM_get_for_im_trace(iSize, jSize, nVarLine, nPointLine, NameVar)
         if(is_proc0(GM_))then
            allocate(Map_DSII(3,2,iSize,jSize), &
                 BufferLine_VI(nVarLine, nPointLine))
            call GM_get_for_im_line(iSize, jSize, Map_DSII, &
                 nVarLine, nPointLine, BufferLine_VI)
         end if
      end if

      !\
      ! Transfer variables from GM to IM
      !/
      if(iProc0Im /= iProc0Gm)then

         if(is_proc0(GM_)) then
            call MPI_send(NameVar,100,MPI_CHARACTER,iProc0Im,&
                 1,iCommWorld,iError)
            call MPI_send(nVarLine,1,MPI_INTEGER,iProc0Im,&
                 2,iCommWorld,iError)
            call MPI_send(nPointLine,1,MPI_INTEGER,iProc0Im,&
                 3,iCommWorld,iError)
            call MPI_send(Map_DSII,size(Map_DSII),MPI_REAL,iProc0Im,&
                 4,iCommWorld,iError)
            call MPI_send(BufferLine_VI,size(BufferLine_VI),MPI_REAL,iProc0Im,&
                 5,iCommWorld,iError)
         end if
         if(is_proc0(IM_))then
            ! setup BufferLine in IM when not sharing proc with GM
            call MPI_recv(NameVar,100,MPI_CHARACTER,iProc0Gm,&
                 1,iCommWorld,iStatus_I,iError)
            call MPI_recv(nVarLine,1,MPI_INTEGER,iProc0Gm,&
                 2,iCommWorld,iStatus_I,iError)
            call MPI_recv(nPointLine,1,MPI_INTEGER,iProc0Gm,&
                 3,iCommWorld,iStatus_I,iError)

            ! Allocate buffer on the IM root processor
            allocate(Map_DSII(3,2,iSize,jSize), &
                 BufferLine_VI(nVarLine, nPointLine))

            ! recieve variables from GM
            call MPI_recv(Map_DSII,size(Map_DSII),MPI_REAL,iProc0Gm,&
                 4,iCommWorld,iStatus_I,iError)
            call MPI_recv(BufferLine_VI,size(BufferLine_VI),MPI_REAL,iProc0Gm,&
                 5,iCommWorld,iStatus_I,iError)
         end if
      end if

      !\
      ! Put variables into IM
      !/
      if(is_proc0(IM_))then
         call IM_put_from_gm_line( &
              iSize, jSize, Map_DSII, &
              nVarLine, nPointLine, BufferLine_VI, NameVar)
         if(DoTest) &
              write(*,*)'IM got from GM: IM iProc, Buffer(1,1)=', &
              iProcWorld,BufferLine_VI(:,1)
      end if

      ! Deallocate buffer on root processors of GM and IM
      if(allocated(BufferLine_VI)) deallocate(BufferLine_VI)
      if(allocated(Map_DSII)) deallocate(Map_DSII)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_ram
    
    !==========================================================================
    subroutine couple_crcm

 
    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere (GM) source\\
    !    Inner  Magnetosphere (IM) target
    !
    ! Send field line volumes, average density and pressure and
    ! geometrical information.
    !EOP

    !\
    ! Coupling variables
    !/

    ! Number of integrals to pass
    integer  :: nIntegral
!    integer, parameter :: nIntegral=6

    ! Names of variables to pass
    character (len=100) :: NameVar
!    character (len=*), parameter :: NameVar='Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:IMF'

    ! Number of variables and points saved into the line data
    integer :: nVarLine, nPointLine

    ! Buffer for the variables on the 2D IM grid and line data
    real, allocatable :: Integral_IIV(:,:,:), BufferLine_VI(:,:)

    ! Buffer for satellite locations   
    real, dimension(:,:,:), allocatable :: SatPos_DII
    
    ! Buffer for satellite names   
    character (len=100), dimension(:), allocatable:: NameSat_I
    
    ! MPI related variables

    ! MPI status variable
    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! General error code
    integer :: iError

    ! Message size
    integer :: nSize

    ! Name of this method
    character (len=*), parameter :: NameSub='couple_gm_im_crcm'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld

    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()

    ! After everything is initialized exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld
    if(DoTest)write(*,*)NameSub,', iProc, GMi_iProc0, i_proc0(IM_)=', &
         iProcWorld,i_proc0(GM_),i_proc0(IM_)

    if(DoMultiFluidIMCoupling) then
       NameVar='vol:Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:rho:p:Hprho:Oprho:Hpp:Opp'
       nIntegral = 10
    else
       NameVar='vol:Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:rho:p'
       nIntegral = 6
    endif


    !\
    ! Allocate buffers both in GM and IM
    !/
    allocate(Integral_IIV(iSize,jSize,nIntegral), stat=iError)
    call check_allocate(iError,NameSub//": Integral_IIV")

    if (nShareSats > 0) then
       allocate(SatPos_DII(4,2,nShareSats), stat=iError)
       call check_allocate(iError,NameSub//": SatPos_DII")
       allocate(NameSat_I(nShareSats),       stat=iError)
       call check_allocate(iError,NameSub//": NameSat_I")
    end if

    if(DoTest)write(*,*)NameSub,', variables allocated',&
         ', iProc:',iProcWorld

    !\
    ! Get field line integrals from GM
    !/
    if(is_proc(GM_)) then
       call GM_get_for_im_trace_crcm(iSize, jSize, NameVar, nVarLine, nPointLine)
       allocate(BufferLine_VI(nVarLine, nPointLine))
       call GM_get_for_im_crcm(Integral_IIV, iSize, jSize, nIntegral, &
            BufferLine_VI, nVarLine, nPointLine, NameVar)
    end if
    !\
    ! If IM sat tracing is enabled, get sat locations from GM
    !/
    if(is_proc(GM_).AND.(nShareSats > 0)) &
         call GM_get_sat_for_im_crcm(SatPos_DII, NameSat_I, nShareSats)

    !\
    ! Transfer variables from GM to IM
    !/
    if(i_proc0(IM_) /= i_proc0(GM_))then
       nSize = iSize*jSize*nIntegral
       if(is_proc0(GM_)) then
          call MPI_send(nVarLine,1,MPI_INTEGER,i_proc0(IM_),&
               1,i_comm(),iError)
          call MPI_send(nPointLine,1,MPI_INTEGER,i_proc0(IM_),&
               1,i_comm(),iError)
          call MPI_send(Integral_IIV,nSize,MPI_REAL,&
               i_proc0(IM_),1,i_comm(),iError)
          call MPI_send(BufferLine_VI,nVarLine*nPointLine,MPI_REAL,&
               i_proc0(IM_),2,i_comm(),iError)
       end if
       if(is_proc0(IM_))then
          ! get nVarLine and nPointLine from GM 0 processor
          call MPI_recv(nVarLine,1,MPI_INTEGER,i_proc0(GM_),&
               1,i_comm(),iStatus_I,iError)
          call MPI_recv(nPointLine,1,MPI_INTEGER,i_proc0(GM_),&
               1,i_comm(),iStatus_I,iError)
          ! setup BufferLine in IM when not sharing proc with GM
          ! If IM and GM share procs not equal to zero proc then 
          ! BufferLine_VI must be first dealloacted and then reallocated 
          ! with proper size
          if (allocated(BufferLine_VI)) deallocate(BufferLine_VI)
          allocate(BufferLine_VI(nVarLine, nPointLine))

          !recieve variables from GM
          call MPI_recv(Integral_IIV,nSize,MPI_REAL,&
               i_proc0(GM_),1,i_comm(),iStatus_I,iError)
          call MPI_recv(BufferLine_VI,nVarLine*nPointLine,MPI_REAL,&
               i_proc0(GM_),2,i_comm(),iStatus_I,iError)
       end if
    end if


        !Broadcast variables inside IM
    if(n_proc(IM_)>1 .and. is_proc(IM_)) then
       nSize = iSize*jSize*nIntegral
       call MPI_bcast(nVarLine,1,MPI_INTEGER,0,i_comm(IM_),iError)
       call MPI_bcast(nPointLine,1,MPI_INTEGER,0,i_comm(IM_),iError)
       call MPI_bcast(Integral_IIV,nSize,MPI_REAL,0,i_comm(IM_),iError)
       if (.not. allocated(BufferLine_VI))&
            allocate(BufferLine_VI(nVarLine, nPointLine))
       call MPI_bcast(BufferLine_VI,nVarLine*nPointLine,MPI_REAL,0,&
            i_comm(IM_),iError)
    endif
      

    !\
    ! Transfer satellite names from GM to IM
    !/   
    
    if(nShareSats > 0 .and. i_proc0(IM_) /= i_proc0(GM_))then
       nSize = nShareSats*100
       if(is_proc0(GM_)) then
          call MPI_send(NameSat_I,nSize,MPI_BYTE,i_proc0(IM_),&
               1,i_comm(),iError)
       endif
       if(is_proc0(IM_)) then
          call MPI_recv(NameSat_I,nSize,MPI_BYTE,i_proc0(GM_),&
               1,i_comm(),iStatus_I,iError)
       endif
    
   ! Transfer satellite locations from GM to IM
       
       nSize = 3*2*nShareSats
       if(is_proc0(GM_)) then
          call MPI_send(SatPos_DII,nSize,MPI_REAL,i_proc0(IM_),&
               1,i_comm(),iError)
       endif
       if(is_proc0(IM_)) then
          call MPI_recv(SatPos_DII,nSize,MPI_REAL,i_proc0(GM_),&
               1,i_comm(),iStatus_I,iError)
       endif
    end if
 
!    ! Broadcast SatPos info in IM
!    if(nShareSats>0 .and. n_proc(IM_)>1 .and. is_proc(IM_)) then
!       call MPI_bcast(nShareSats,1,MPI_INTEGER,0,i_comm(IM_),iError)
!       nSize = nShareSats*100
!       call MPI_bcast(NameSat_I,nSize,MPI_BYTE,0,i_comm(IM_),iError)
!       nSize = 3*2*nShareSats
!       call MPI_bcast(SatPos_DII,nSize,MPI_REAL,0,i_comm(IM_),iError)
!    endif

    if(DoTest)write(*,*)NameSub,', variables transferred',&
         ', iProc:',iProcWorld



       
    !\
    ! Put variables into IM
    !/
    if(is_proc(IM_)) then
       call IM_put_from_gm_crcm(Integral_IIV,iSize,jSize,nIntegral,&
            BufferLine_VI,nVarLine,nPointLine,NameVar,tSimulation)
    endif

    if(is_proc0(IM_))then
       if(nShareSats > 0) &
            call IM_put_sat_from_gm(nShareSats, NameSat_I, SatPos_DII)
       if(DoTest) &
            write(*,*)'IM got from GM: IM iProc, Buffer(1,1)=',&
            iProcWorld,Integral_IIV(1,1,:)
    end if
    

    !\
    ! Deallocate buffer to save memory
    !/
    deallocate(Integral_IIV, BufferLine_VI)

    if (nShareSats > 0) then
       deallocate(NameSat_I)
       deallocate(SatPos_DII)
    end if

    if(DoTest)write(*,*)NameSub,', variables deallocated',&
         ', iProc:',iProcWorld

    if(DoTest)write(*,*)NameSub,' finished, iProc=',iProcWorld

    ! if(DoTest.and.is_proc0(IM_)) call IM_print_variables('GM',iError)      
      
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

    !\
    ! General coupling variables
    !/

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_im_gm'

    logical :: DoTest, DoTestMe
    integer :: iProcWorld

    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    iProcWorld = i_proc()
    call couple_mpi

    if(DoTest)write(*,*)NameSub,': finished iProc=',iProcWorld

    if(DoTest.and.is_proc0(GM_)) call GM_print_variables('IM')

  contains

    !=========================================================================
    subroutine couple_mpi

      character (len=*), parameter :: NameSubSub=NameSub//'.couple_mpi'

      ! Number of variables to pass
      integer :: nVarImGm
     ! Number of variables for Multi-Fluid to pass                                              

      ! Variable to pass is pressure
      character (len=100) :: NameVar

      ! Buffer for the variables on the 2D IM grid
      real, dimension(:,:,:), allocatable :: Buffer_IIV

      ! MPI related variables

      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)

      ! General error code
      integer :: iError

      ! Message size
      integer :: nSize

      ! Communicator and logicals to simplify message passing and execution
      logical :: IsUninitialized = .true., UseMe=.true.
      !-----------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,' starting, iProc=',iProcWorld
      if(DoTest)write(*,*)NameSubSub,', iProc, GMi_iProc0, IMi_iProc0=', &
           iProcWorld,iProc0Gm,iProc0Im

     if(DoMultiFluidIMCoupling)then
        NameVar='p:rho:Hpp:Opp:Hprho:Oprho'
        nVarImGm=6
     else
        NameVar='p:rho'
        nVarImGm=2
     end if

      !\
      ! Allocate buffers both in GM and IM
      !/
      allocate(Buffer_IIV(iSize,jSize,nVarImGm), stat=iError)
      call check_allocate(iError,NameSubSub//": Buffer_IIV")

      if(DoTest)write(*,*)NameSubSub,', variables allocated',&
           ', iProc:',iProcWorld

      !\
      ! Get pressure from IM
      !/
      if(is_proc(IM_)) &
           call IM_get_for_gm(Buffer_IIV,iSize,jSize,nVarImGm,NameVar)
      !\
      ! Transfer variables from IM to GM
      !/ 
      nSize = iSize*jSize*nVarImGm
      if(iProc0Im /= iProc0Gm)then
         if(is_proc0(IM_)) &
              call MPI_send(Buffer_IIV,nSize,MPI_REAL,iProc0Gm,&
              1,iCommWorld,iError)
         if(is_proc0(GM_)) &
              call MPI_recv(Buffer_IIV,nSize,MPI_REAL,iProc0Im,&
              1,iCommWorld,iStatus_I,iError)
      end if

      ! Broadcast variables inside GM
      if(n_proc(GM_)>1 .and. is_proc(GM_)) &
           call MPI_bcast(Buffer_IIV,nSize,MPI_REAL,0,i_comm(GM_),iError)

      if(DoTest)write(*,*)NameSubSub,', variables transferred',&
           ', iProc:',iProcWorld

      !\
      ! Put variables into GM
      !/
      if(is_proc(GM_))then
         call GM_put_from_im(Buffer_IIV,iSize,jSize,nVarImGm,NameVar)
         if(DoTest) &
              write(*,*)NameSubSub//' iProc, Buffer(1,1,1)=',&
              iProcWorld,Buffer_IIV(1,1,1)
      end if

      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Buffer_IIV)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_im_gm


end module CON_couple_gm_im
