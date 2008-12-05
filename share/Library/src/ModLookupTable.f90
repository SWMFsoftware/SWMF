module ModLookupTable

  ! Use lookup tables to calculate lookup properties, such as
  ! equation of state, opacities, ionization level etc. 
  ! For example interpolate pressure as a function of the logarithm of 
  ! density and logarithm of internal energy. All variables are in SI units.

  use ModReadParam, ONLY: read_var
  use ModPlotFile,  ONLY: read_plot_file, save_plot_file
  use ModUtilities, ONLY: split_string
  use ModInterpolate, ONLY: bilinear

  implicit none
  SAVE

  private ! except

  public:: read_lookup_table_param  ! read parameters of the lookup table(s)
  public:: make_lookup_table        ! create table from calculations (and save)
  public:: load_lookup_table        ! load lookup table from file
  public:: interpolate_lookup_table ! interpolate from lookup table
  public:: test_lookup_table        ! unit test

  integer, public, parameter:: MaxTable = 2, TableRhoE_ = 1, TableRhoP_ = 2
  character(len=8), public, parameter:: NameTable_I(MaxTable) = &
       (/ "TableRhoE", "TableRhoP" /)

  ! private variables

  type TableType
     integer:: nValue                     ! number of values in each element
     integer:: nIndex_I(2)                ! number of columns and rows
     real   :: IndexMin_I(2)              ! minimum values for indexes
     real   :: IndexMax_I(2)              ! maximum values for indexes
     real   :: dIndex_I(2)                ! increment
     logical:: IsLogIndex_I(2)            ! true if arguments are logarithmic
     real, allocatable:: Value_VII(:,:,:) ! array of actual values
     character(len=100):: NameTable       ! description of table
     character(len=100):: NameVar         ! name of indexes and values
     character(len=100):: NameFile        ! file name containing the table
     character(len=100):: NameCommand     ! action for the table
  end type

  ! The array of tables
  type(TableType), target :: Table_I(MaxTable)

  ! Array for variable names
  integer, parameter:: MaxVar=20
  character(len=20):: NameVar_I(MaxVar)

contains

  !==========================================================================
  subroutine read_lookup_table_param(NameCommand)

    character(len=100), intent(in):: NameCommand

    integer :: iTable, iIndex
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub = 'read_lookup_table_param'
    !-----------------------------------------------------------------------
    select case(NameCommand)
    case("#TABLEMAKE", "#TABLESAVE")
       call read_var('iTable', iTable)
       Ptr => Table_I(iTable)
       Ptr%NameCommand = NameCommand
       if(NameCommand == "#TABLESAVE") &
            call read_var('NameFile for '//NameTable_I(iTable), Ptr%NameFile)
       call read_var('NameTable', Ptr%NameTable)
       call read_var('NameVar',   Ptr%NameVar)
       call split_string(Ptr%NameVar, MaxVar, NameVar_I, Ptr%nValue)
       ! Do not count the names of the indexes
       Ptr%nValue = Ptr%nValue - 2
       ! Figure out which index is logarithmic
       Ptr%IsLogIndex_I = index(NameVar_I(1:2), "log") == 1
       
       do iIndex = 1, 2
          call read_var('nIndex',     Ptr%nIndex_I(iIndex))
          call read_var('IndexMin',   Ptr%IndexMin_I(iIndex))
          call read_var('IndexMax',   Ptr%IndexMax_I(iIndex))
          
          ! Take logarithm of the ranges if logarithmic
          if(Ptr%IsLogIndex_I(iIndex)) then
             Ptr%IndexMin_I(iIndex) = log(Ptr%IndexMin_I(iIndex))
             Ptr%IndexMax_I(iIndex) = log(Ptr%IndexMax_I(iIndex))
          end if
       end do
       ! Calculate increments
       Ptr%dIndex_I = &
            (Ptr%IndexMax_I - Ptr%IndexMin_I)/(Ptr % nIndex_I - 1)

    case("#TABLELOAD")
       call read_var('iTable', iTable)
       Ptr => Table_I(iTable)
       Ptr%NameCommand = NameCommand
       call read_var('NameFile for '//NameTable_I(iTable), Ptr%NameFile)
       call load_lookup_table(iTable)

    case default
       call CON_stop(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine read_lookup_table_param

  !===========================================================================

  subroutine load_lookup_table(iTable)

    integer, intent(in) :: iTable

    type(TableType), pointer:: Ptr
    integer :: nVar

    character(len=*), parameter:: NameSub = 'load_lookup_table'
    !------------------------------------------------------------------------

    Ptr => Table_I(iTable)
    if(Ptr%NameCommand /= "#TABLELOAD") RETURN

    ! Make sure it is not loaded again
    Ptr%NameCommand = ""

    call read_plot_file( Ptr%NameFile,      &
         StringHeaderOut = Ptr%Nametable,   &
         n1Out           = Ptr%nIndex_I(1), &
         n2Out           = Ptr%nIndex_I(2), &
         nVarOut         = Ptr%nValue,      &
         NameVarOut      = Ptr%NameVar)

    ! Figure out which index is logarithmic
    call split_string(Ptr%NameVar, MaxVar, NameVar_I, nVar)
    Ptr%IsLogIndex_I = index(NameVar_I(1:2), "log") == 1

    if(allocated(Ptr%Value_VII)) deallocate(Ptr%Value_VII)
    allocate(Ptr%Value_VII(Ptr%nValue, Ptr%nIndex_I(1), Ptr%nIndex_I(2)))
    
    call read_plot_file( Ptr%NameFile,    &
         CoordMinOut_D = Ptr%IndexMin_I,  &
         CoordMaxOut_D = Ptr%IndexMax_I,  &
         VarOut_VII    = Ptr%Value_VII)

    ! Calculate increments
    Ptr%dIndex_I = (Ptr%IndexMax_I - Ptr%IndexMin_I)/(Ptr % nIndex_I - 1)

  end subroutine load_lookup_table

  !===========================================================================

  subroutine make_lookup_table(iTable, calc_table_var)

    ! Fill in tables using the subroutine calc_table_var

    integer, intent(in):: iTable
    interface
       subroutine calc_table_var(Arg1, Arg2, Value_V)
         real,    intent(in) :: Arg1, Arg2
         real,    intent(out):: Value_V(:)
       end subroutine calc_table_var
    end interface

    integer :: i1, i2, n1, n2, nValue
    logical :: IsLog1, IsLog2
    real :: Index1Min, Index2Min, dIndex1, dIndex2, Index1, Index2
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub='make_lookup_table'
    !------------------------------------------------------------------------

    Ptr => Table_I(iTable)

    if(  Ptr%NameCommand /= "#TABLEMAKE" .and. &
         Ptr%NameCommand /= "#TABLESAVE") RETURN

    ! Use simple scalars for sake of legibility
    n1     = Ptr%nIndex_I(1)
    n2     = Ptr%nIndex_I(2)
    nValue = Ptr%nValue
    IsLog1 = Ptr%IsLogIndex_I(1)
    IsLog2 = Ptr%IsLogIndex_I(2)

    Index1Min = Ptr%IndexMin_I(1)
    Index2Min = Ptr%IndexMin_I(2)
    dIndex1   = Ptr%dIndex_I(1)
    dIndex2   = Ptr%dIndex_I(2)

    ! Allocate Value_VII array
    if(allocated(Ptr%Value_VII)) deallocate(Ptr%Value_VII)
    allocate(Ptr%Value_VII(nValue, n1, n2))

    ! Fill up lookup table, This could be done in parallel !!!
    do i2 = 1, n2
       Index2 = Index2Min + (i2 - 1)*dIndex2
       if(IsLog2) Index2 = exp(Index2)
       do i1 = 1, n1
          Index1 = Index1Min + (i1 - 1)*dIndex1
          if(IsLog1) Index1 = exp(Index1)
          call calc_table_var( Index1, Index2, Ptr%Value_VII(:,i1,i2))
       end do
    end do

    if(Ptr%NameCommand == "#TABLESAVE") call save_plot_file( &
         Ptr%NameFile,                    &
         StringHeaderIn = Ptr%NameTable,  &
         NameVarIn      = Ptr%NameVar,    &
         CoordMinIn_D   = Ptr%IndexMin_I, &
         CoordMaxIn_D   = Ptr%IndexMax_I, &
         VarIn_VII      = Ptr%Value_VII)

    ! Make sure it is not saved again
    Ptr%NameCommand = ""

  end subroutine make_lookup_table

  !===========================================================================

  subroutine interpolate_lookup_table(iTable, Arg1In, Arg2In, Value_V)

    integer, intent(in) :: iTable            ! table
    real,    intent(in) :: Arg1In, Arg2In    ! input arguments
    real,    intent(out):: Value_V(:)        ! output values

    real :: Arg_I(2)
    type(TableType), pointer:: Ptr
    
    character(len=*), parameter:: NameSub='interpolate_lookup_table'
    !--------------------------------------------------------------------------
    Arg_I = (/Arg1In, Arg2In/)
    Ptr => Table_I(iTable)

    ! This line is broken so that emacs does not get confused
    where(Ptr%IsLogIndex_I) &
         Arg_I = log(Arg_I)

    Value_V = bilinear(Ptr%Value_VII, Ptr%nValue, &
         1, Ptr%nIndex_I(1), 1, Ptr%nIndex_I(2), &
         (Arg_I - Ptr%IndexMin_I)/Ptr%dIndex_I  + 1)

  end subroutine interpolate_lookup_table

  !===========================================================================

  subroutine test_lookup_table

    ! testing the read_lookup_table_param is left for the functionality tests

    type(TableType), pointer :: Ptr, Ptr2
    real :: p_I(3), pGood_I(3)

    character(len=*), parameter:: NameSub = 'test_lookup_table'
    !------------------------------------------------------------------------
    Ptr => Table_I(TableRhoE_)

    Ptr%NameCommand = "#TABLESAVE"
    Ptr%NameFile    = "test_lookup_table1.out"
    Ptr%NameTable   = "eos: p_i(rho,e) for i=0,1,2 materials"
    Ptr%NameVar     = "logrho e pXe pBe pPl"
    Ptr%nValue      = 3
    Ptr%IsLogIndex_I= (/.true., .false./)
    Ptr%nIndex_I    = (/15, 10/)
    Ptr%IndexMin_I  = (/log(0.001),   1.0/)
    Ptr%IndexMax_I  = (/log(1000.0), 10.0/)
    Ptr%dIndex_I    = (Ptr%IndexMax_I - Ptr%IndexMin_I)/(Ptr%nIndex_I - 1)

    write(*,*)'testing make_lookup_table'
    call make_lookup_table(TableRhoE_, eos_rho_e)

    write(*,*)'testing interpolate_lookup_table'    
    call interpolate_lookup_table(TableRhoE_, 1.0, 2.0, p_I)
    ! rho=1.0 is exactly in the middle, e=2.0 is also an index, so exact result

    pGood_I = (/ 4./3., 4./5., 3. /)
    if(any(abs(p_I - pGood_I) > 1e-5))then
       write(*,*)'p_I=',p_I,' is different from pGood_I=',pGood_I
       call CON_stop(NameSub)
    end if

    write(*,*)'testing load_lookup_table'

    ! Load the saved file into the second table
    Ptr2 => Table_I(TableRhoP_)
    Ptr2%NameCommand = "#TABLELOAD"
    Ptr2%NameFile    = "test_lookup_table1.out"

    call load_lookup_table(TableRhoP_)

    if(Ptr2%NameTable /= Ptr%NameTable) call CON_stop(NameSub // &
         ' NameTable='//trim(Ptr2%NameTable)//' is different from '// &
         trim(Ptr2%NameTable))

    if(Ptr2%NameVar /= Ptr%NameVar) call CON_stop(NameSub // &
         ' NameVar='//trim(Ptr2%NameVar)//' is different from '// &
         trim(Ptr2%NameVar))

    if(Ptr%nValue /= Ptr2%nValue)then
       write(*,*)'nValue=',Ptr2%nValue,' is different from ',Ptr%nValue
       call CON_stop(NameSub)
    end if

    write(*,*)'testing interpolate_lookup_table for loaded table'
    call interpolate_lookup_table(TableRhoE_, 1.0, 2.0, p_I)
    if(any(abs(p_I - pGood_I) > 1e-5))then
       write(*,*)'p_I=',p_I,' is different from pGood_I=',pGood_I
       call CON_stop(NameSub)
    end if


  end subroutine test_lookup_table

  !===========================================================================

  subroutine eos_rho_e(rho, e, p_I)

    real, intent(in):: rho, e
    real, intent(out):: p_I(:)

    p_I(1) = (2./3.)*e
    p_I(2) = (2./5.)*e
    p_I(3) = e + rho

  end subroutine eos_rho_e

end module ModLookupTable
