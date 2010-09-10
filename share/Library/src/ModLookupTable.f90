module ModLookupTable

  ! Use lookup tables to calculate lookup properties, such as
  ! equation of state, opacities, ionization level etc. 
  ! For example interpolate pressure as a function of the logarithm of 
  ! density and logarithm of internal energy. All variables are in SI units.

  use ModReadParam, ONLY: read_var
  use ModPlotFile,  ONLY: read_plot_file, save_plot_file
  use ModUtilities, ONLY: split_string, lower_case
  use ModInterpolate, ONLY: bilinear, find_cell
  use ModMpi

  implicit none
  SAVE

  private ! except

  public:: init_lookup_table        ! set parameters of  the lookup table(s)
  public:: read_lookup_table_param  ! read parameters of the lookup table(s)
  public:: i_lookup_table           ! function returning the index of table
  public:: make_lookup_table        ! create table from calculations (and save)
  public:: interpolate_lookup_table ! interpolate from lookup table
  public:: test_lookup_table        ! unit test

  integer, public, parameter:: MaxTable = 20 ! maximum number of tables
  integer, public :: nTable = 0     ! actual number of tables

  ! private variables

  interface interpolate_lookup_table
     module procedure interpolate_with_known_arg   !Both arguments are known
     module procedure interpolate_with_known_val   !Table value is given
  end interface

  type TableType
     character(len=100):: NameTable        ! unique name for identification
     character(len=4)  :: NameCommand      ! command: load, make, save
     character(len=100):: NameFile         ! file name containing the table
     character(len=10) :: TypeFile         ! file type (ascii, real4, real8)
     character(len=100):: StringDescription! description of table
     character(len=500):: NameVar          ! name of indexes and values
     integer:: nValue                      ! number of values in each element
     integer:: nIndex_I(2)                 ! number of columns and rows
     real   :: IndexMin_I(2)               ! minimum values for indexes
     real   :: IndexMax_I(2)               ! maximum values for indexes
     real   :: dIndex_I(2)                 ! increment of indexes
     logical:: IsLogIndex_I(2)             ! true if arguments are logarithmic
     real, allocatable :: Value_VII(:,:,:) ! array of actual values
  end type

  ! The array of tables
  type(TableType), target :: Table_I(MaxTable)

  ! Array for variable names
  integer, parameter:: MaxVar = 200
  character(len=20):: NameVar_I(MaxVar)

contains
  !==========================================================================
  subroutine init_lookup_table(NameTable, NameCommand, NameVar, &
    nIndex_I, IndexMin_I, IndexMax_I, &
    NameFile, TypeFile, StringDescription)

    character(len=*), intent(in):: NameTable, NameCommand 

    character(len=*),   optional, intent(in):: NameVar
    integer,            optional, intent(in):: nIndex_I(2)
    real, dimension(2), optional, intent(in):: IndexMin_I, IndexMax_I
    character(len=*),   optional, intent(in):: &
                                     NameFile, TypeFile, StringDescription


    integer :: iTable, iIndex, nTable2Read, nFile2Read
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub = 'init_lookup_table'
    !---------------------------------------

    !If NameTable includes the combination like Prefix{Xe Be Pl}Suffix
    !the tables PrefixXeSuffix, PrefixBeSuffix and PrefixPlSuffix will be
    !created. This may make sense only for loading everal similar tables
    !like: Xe_eos, Be_eos, Pl_eos, therefore the loop over the tables
    !ends just after reading the file names
    

    ! Check if the table has been set already (say in a previous session)
       
       
    iTable = i_lookup_table(NameTable)
    if(iTable < 0)then
       ! new table
       nTable = nTable + 1
       
       if(nTable > MaxTable)then
          write(*,*)NameSub,' MaxTable =',MaxTable
          call CON_stop(NameSub//': number of tables exceeded MaxTable')
       end if
       
       iTable = nTable
    end if
       
    ! For sake of more concise source code, use a pointer to the table
    Ptr => Table_I(iTable)
    Ptr%NameTable = NameTable
    
    Ptr%NameCommand = NameCommand
    call lower_case(Ptr%NameCommand)
    
    select case(Ptr%NameCommand)
    case("load","save")
       Ptr%NameFile = NameFile
       Ptr%TypeFile = TypeFile
       
    case("make")
       ! will be done below
    case default
       call CON_stop(NameSub//': unknown command='//Ptr%NameCommand)
    end select
    
    
    if(NameCommand == "load")&
         call load_lookup_table(iTable)

   
    if(NameCommand == "load") RETURN
   

    if(present(StringDescription))&
         Ptr%StringDescription = StringDescription
 
    Ptr%NameVar = NameVar

    call split_string(Ptr%NameVar, MaxVar, NameVar_I, Ptr%nValue, &
         UseArraySyntaxIn=.true.)

    ! Do not count the names of the indexes
    Ptr%nValue = Ptr%nValue - 2
   
    ! Figure out which index is logarithmic
   
    Ptr%nIndex_I   = nIndex_I
    Ptr%IndexMin_I = IndexMin_I
    Ptr%IndexMax_I = IndexMax_I
    
    Ptr%IsLogIndex_I = index(NameVar_I(1:2), "log") == 1
  
    ! Take logarithm of the ranges if logarithmic
    do iIndex = 1, 2
       if(Ptr%IsLogIndex_I(iIndex)) then
          Ptr%IndexMin_I(iIndex) = log10(Ptr%IndexMin_I(iIndex))
          Ptr%IndexMax_I(iIndex) = log10(Ptr%IndexMax_I(iIndex))
       end if
    end do
    ! Calculate increments
    Ptr%dIndex_I = (Ptr%IndexMax_I - Ptr%IndexMin_I)/(Ptr % nIndex_I - 1)
  end subroutine init_lookup_table
  !==========================================================================
  subroutine read_lookup_table_param

    ! Read parameters for one table. The table is identified by a name string

    character(len=100):: NameTable, NameFile, NameCommand, TypeFile
    integer :: iTable, iIndex, nTable2Read, iTableLoop, nFile2Read
    integer, parameter:: MaxString = 200
    character(LEN=100), dimension(MaxString):: NameTable_I, NameFile_I
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub = 'read_lookup_table_param'
    !-----------------------------------------------------------------------
    call read_var('NameTable', NameTable)

    call read_var('NameCommand', NameCommand)
    call lower_case(NameCommand)

    !If NameTable includes the combination like Prefix{Xe Be Pl}Suffix
    !the tables PrefixXeSuffix, PrefixBeSuffix and PrefixPlSuffix will be
    !created. This may make sense only for loading several similar tables
    !like: Xe_eos, Be_eos, Pl_eos, therefore the loop over the tables
    !ends just after reading the file names

    call check_braces(NameTable, NameTable_I, nTable2Read)

    if(nTable2read /=1.and.NameCommand/="load")call CON_stop(&
         'Multiple names can be used only for tables to be loaded')

    do iTableLoop = 1, nTable2Read

       NameTable = NameTable_I(iTableLoop)

       ! Check if the table has been set already (say in a previous session)
       
       
       iTable = i_lookup_table(NameTable)
       if(iTable < 0)then
          ! new table
          nTable = nTable + 1
          
          if(nTable > MaxTable)then
             write(*,*)NameSub,' MaxTable =',MaxTable
             call CON_stop(NameSub//': number of tables exceeded MaxTable')
          end if
       
          iTable = nTable
       end if

       ! For sake of more concise source code, use a pointer to the table
       Ptr => Table_I(iTable)
       Ptr%NameTable = NameTable
       if(iTableLoop == 1)then

          ! Read the parameters for this table
         
          select case(NameCommand)
          case("save")
             call read_var('NameFile', NameFile_I(1))
             call read_var('TypeFile', TypeFile)
          case("load")
             call read_var('NameFile', NameFile)
             call check_braces(NameFile, NameFile_I, nFile2Read)
             if(nTable2read /= nFile2Read)call CON_stop(&
                  'The number of tables to load is not equal to the number of files to read')

             call read_var('TypeFile', TypeFile)
             
          case("make")
             ! will be done below
          case default
             call CON_stop(NameSub//': unknown command='//Ptr%NameCommand)
          end select
       end if

       Ptr%NameCommand = NameCommand
       Ptr%NameFile = NameFile_I(iTableLoop)
       Ptr%TypeFile = TypeFile
       if(NameCommand == "load")&
            call load_lookup_table(iTable)

 
    end do
   
    if(NameCommand == "load") RETURN
    
    
    call read_var('StringDescription', Ptr%StringDescription)
    call read_var('NameVar',           Ptr%NameVar)
    
    call split_string(Ptr%NameVar, MaxVar, NameVar_I, Ptr%nValue, &
         UseArraySyntaxIn=.true.)
    
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
          Ptr%IndexMin_I(iIndex) = log10(Ptr%IndexMin_I(iIndex))
          Ptr%IndexMax_I(iIndex) = log10(Ptr%IndexMax_I(iIndex))
       end if
    end do
    ! Calculate increments
    Ptr%dIndex_I = (Ptr%IndexMax_I - Ptr%IndexMin_I)/(Ptr % nIndex_I - 1)
  contains
    
    subroutine check_braces(Name, Name_I, nString)
      character(LEN=*), intent(in) :: Name
      character(LEN=100), intent(out) :: Name_I(MaxString)
      integer, intent(out) :: nString

      integer:: iBracePosition1, iBracePosition2, iString

      !-----------------------------------------

      iBracePosition1 = index(Name,'{')
      iBracePosition2 = index(Name,'}')

      if(iBracePosition1 < 1 .or. iBracePosition2 < 1 .or.&
           iBracePosition2< iBracePosition1 )then
         nString = 1 
         Name_I(1) = Name
         return
      end if

      call split_string(Name(iBracePosition1 + 1:iBracePosition2 - 1),&
           MaxString, Name_I, nString)

      if(iBracePosition1 > 1)then
         do iString = 1, nString
            Name_I(iString) = Name(1:iBracePosition1-1)//trim(Name_I(iString))
         end do
      end if

      if(iBracePosition2 < len_trim(Name))then
         do iString = 1, nString
            Name_I(iString) =trim(Name_I(iString))//Name(iBracePosition2 + 1: len_trim(Name))
         end do
      end if

    end subroutine check_braces
  end subroutine read_lookup_table_param

  !===========================================================================

  integer function i_lookup_table(Name)

    ! return the index of the lookup table based on its name
    ! return -1 if the table was not found

    character(len=*), intent(in):: Name
    integer :: iTable
    !------------------------------------------------------------------------

    do iTable = 1, nTable
       if(Table_I(iTable)%NameTable == Name) then
          i_lookup_table = iTable
          RETURN
       end if
    end do
    i_lookup_table = -1

  end function i_lookup_table

  !===========================================================================

  subroutine load_lookup_table(iTable)

    integer, intent(in) :: iTable

    type(TableType), pointer:: Ptr
    integer :: nVar

    character(len=*), parameter:: NameSub = 'load_lookup_table'
    !------------------------------------------------------------------------

    if(iTable > nTable) call CON_stop(NameSub//' iTable larger than nTable')

    Ptr => Table_I(iTable)

    call read_plot_file( Ptr%NameFile,            &
         TypeFileIn      = Ptr%TypeFile,          &
         StringHeaderOut = Ptr%StringDescription, &
         n1Out           = Ptr%nIndex_I(1),       &
         n2Out           = Ptr%nIndex_I(2),       &
         nVarOut         = Ptr%nValue,            &
         NameVarOut      = Ptr%NameVar)

    ! Figure out which index is logarithmic
    call split_string(Ptr%NameVar, MaxVar, NameVar_I, nVar)
    Ptr%IsLogIndex_I = index(NameVar_I(1:2), "log") == 1

    if(allocated(Ptr%Value_VII)) deallocate(Ptr%Value_VII)
    allocate(Ptr%Value_VII(Ptr%nValue, Ptr%nIndex_I(1), Ptr%nIndex_I(2)))
    
    call read_plot_file( Ptr%NameFile,    &
         TypeFileIn      = Ptr%TypeFile,  &
         CoordMinOut_D = Ptr%IndexMin_I,  &
         CoordMaxOut_D = Ptr%IndexMax_I,  &
         VarOut_VII    = Ptr%Value_VII)

    ! Calculate increments
    Ptr%dIndex_I = (Ptr%IndexMax_I - Ptr%IndexMin_I)/(Ptr % nIndex_I - 1)

  end subroutine load_lookup_table

  !===========================================================================

  subroutine make_lookup_table(iTable, calc_table_var, iComm)

    ! Fill in table iTable using the subroutine calc_table_var
    ! The optional communicator allows for parallel execution

    integer, intent(in):: iTable  ! table index
    interface
       subroutine calc_table_var(iTable, Arg1, Arg2, Value_V)
         integer, intent(in) :: iTable
         real,    intent(in) :: Arg1, Arg2
         real,    intent(out):: Value_V(:)
       end subroutine calc_table_var
    end interface
    integer, optional, intent(in):: iComm

    integer:: iProc, nProc, iError
    integer:: i1, i2, n1, n2, nValue
    logical:: IsLog1, IsLog2
    real   :: Index1Min, Index2Min, dIndex1, dIndex2, Index1, Index2
    real, allocatable:: Value_VII(:,:,:)
    type(TableType), pointer:: Ptr

    character(len=*), parameter:: NameSub='make_lookup_table'
    !------------------------------------------------------------------------

    Ptr => Table_I(iTable)

    if(Ptr%NameCommand /= "make" .and. Ptr%NameCommand /= "save") &
         RETURN

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

    ! Get processor index and total number of processors
    if(present(iComm))then
       call MPI_comm_rank(iComm,iProc,iError)
       call MPI_comm_size(iComm,nProc,iError)
    else
       iProc = 0
       nProc = 1
    end if

    ! Allocate Value_VII array
    if(allocated(Ptr%Value_VII)) deallocate(Ptr%Value_VII)
    allocate(Value_VII(nValue, n1, n2), Ptr%Value_VII(nValue, n1, n2))
    Value_VII = 0.0

    ! Fill up lookup table in parallel
    do i2 = iProc+1, n2, nProc
       Index2 = Index2Min + (i2 - 1)*dIndex2
       if(IsLog2) Index2 = 10**Index2
       do i1 = 1, n1
          Index1 = Index1Min + (i1 - 1)*dIndex1
          if(IsLog1) Index1 = 10**Index1
          call calc_table_var(iTable, Index1, Index2, Value_VII(:,i1,i2))
       end do
    end do

    ! Collect (or copy) result into table
    if(nProc > 1)then
       call MPI_allreduce(Value_VII, Ptr%Value_VII, n1*n2*nValue, MPI_REAL, &
            MPI_SUM, iComm, iError)
    else
       Ptr%Value_VII = Value_VII
    end if

    deallocate(Value_VII)

    if(Ptr%NameCommand == "save" .and. iProc == 0) call save_plot_file( &
         Ptr%NameFile,                                 &
         TypeFileIn     = Ptr%TypeFile,                &
         StringHeaderIn = Ptr%StringDescription,       &
         NameVarIn      = Ptr%NameVar,                 &
         CoordMinIn_D   = Ptr%IndexMin_I,              &
         CoordMaxIn_D   = Ptr%IndexMax_I,              &
         VarIn_VII      = Ptr%Value_VII)

    ! Make sure that all processors are done
    call MPI_barrier(iComm, iError)

    ! Make sure it is not saved again
    Ptr%NameCommand = "done"

  end subroutine make_lookup_table

  !===========================================================================

  subroutine interpolate_with_known_arg(iTable, Arg1In, Arg2In, Value_V, &
       DoExtrapolate)

    ! Return the array of values Value_V corresponding to arguments
    ! Arg1In and Arg2In in iTable. Use a bilinear interpolation.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a 
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    real,    intent(in) :: Arg1In, Arg2In    ! input arguments
    real,    intent(out):: Value_V(:)        ! output values

    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    real :: Arg_I(2)
    type(TableType), pointer:: Ptr
    
    character(len=*), parameter:: NameSub='interpolate_lookup_table'
    !--------------------------------------------------------------------------
    Arg_I = (/Arg1In, Arg2In/)
    Ptr => Table_I(iTable)

    ! This line is broken so that emacs does not get confused
    where(Ptr%IsLogIndex_I) &
         Arg_I = log10(Arg_I)

    ! If value is outside table, use the last value (works well for constant)
    Value_V = bilinear(Ptr%Value_VII, Ptr%nValue, &
         1, Ptr%nIndex_I(1), 1, Ptr%nIndex_I(2), &
         (Arg_I - Ptr%IndexMin_I)/Ptr%dIndex_I  + 1, &
         DoExtrapolate = DoExtrapolate)

  end subroutine interpolate_with_known_arg

  !===========================================================================

  subroutine interpolate_with_known_val(&
       iTable, iVal, ValIn, Arg2In, Value_V, &
       Arg1Out, DoExtrapolate)

    ! Return the array of values Value_V corresponding to the argument
    ! Arg2In in iTable, with Arg1 being calculated from the condition, that 
    ! Value_V(iVal) equals the given value, ValIn.
    !
    ! Use a bilinear interpolation.
    ! If DoExtrapolate is not present, stop with an error if the arguments
    ! are out of range. If it is present and false, return the value of
    ! the closest element in the table. If it is present and true, do a 
    ! linear extrapolation.

    integer, intent(in) :: iTable            ! table
    integer, intent(in) :: iVal              ! which value is known
    real,    intent(in) :: ValIn             ! known table value
    real,    intent(in) :: Arg2In            ! second input argument
    real,    intent(out):: Value_V(:)        ! output values
    
    real, optional, intent(out) :: Arg1Out        ! optional calculated Arg
    logical, optional, intent(in):: DoExtrapolate ! optional extrapolation

    real :: Arg1, Arg2
    type(TableType), pointer:: Ptr

    real    :: Dx1, Dx2, Dy1, Dy2
    integer :: i1, i2, j1, j2
    
    character(len=*), parameter:: NameSub='interpolate_lookup_table'
    !--------------------------------------------------------------------------
    
    Ptr => Table_I(iTable)

    Arg2 = Arg2In
    
    If(Ptr%IsLogIndex_I(2)) Arg2 = log10(Arg2)

    call find_cell(1, Ptr%nIndex_I(2), &
         (Arg2- Ptr%IndexMin_I(2))/Ptr%dIndex_I(2)+ 1 , &
         j1, Dy1, &
         DoExtrapolate=DoExtrapolate, &
         StringError = 'Called from '//NameSub)

    j2 = j1 + 1; Dy2 = 1.0 - Dy1


    call find_cell(1, Ptr%nIndex_I(1), &
         ValIn, &
         i1, Dx1, &
         Dy2*Ptr%Value_VII(iVal,:,j1) + &
         Dy1*Ptr%Value_VII(iVal,:,j2), &
         DoExtrapolate, &
         'Called from '//NameSub)
    i2 = i1 + 1; Dx2 = 1.0 - Dx1

    ! If value is outside table, use the last value (works well for constant)
    Value_V = Dy2*( Dx2*Ptr%Value_VII(:,i1,j1)   &
         +          Dx1*Ptr%Value_VII(:,i2,j1))  &
         +    Dy1*( Dx2*Ptr%Value_VII(:,i1,j2)   &
         +          Dx1*Ptr%Value_VII(:,i2,j2))
    if(present(Arg1Out))then
       Arg1Out = (i1 - 1 + Dx1)*Ptr%dIndex_I(1) + Ptr%IndexMin_I(1)
       if(Ptr%IsLogIndex_I(1)) Arg1Out = 10**Arg1Out
    end if
       
  end subroutine interpolate_with_known_val

  !===========================================================================

  subroutine test_lookup_table

    ! testing the read_lookup_table_param is left for the functionality tests

    type(TableType), pointer :: Ptr, Ptr2
    integer :: iTable, iProc, iError
    real :: p_I(3), pGood_I(3), Arg

    character(len=*), parameter:: NameSub = 'test_lookup_table'
    !------------------------------------------------------------------------
    call MPI_comm_rank(MPI_COMM_WORLD,iProc,iError)
    call init_lookup_table(&
         NameTable   = "RhoE",                  &
         NameCommand = "save",                  &
         NameVar     = "logrho e pXe pBe pPl",  &
         NameFile    = "test_lookup_table1.out",& 
         TypeFile    = "ascii",                 &
         nIndex_I    = (/15, 10/),              &
         IndexMin_I  = (/0.001,   1.0/),        &
         IndexMax_I  = (/1000.0, 10.0/))
    

    if(iProc==0) write(*,*)'testing i_lookup_table'
   
    iTable = i_lookup_table("xxx")
    if(iTable /= -1)then
       write(*,*)'iTable = ',iTable,' should be -1'
       call CON_stop(NameSub)
    end if
    iTable = i_lookup_table("RhoE")
    if(iTable /= 1)then
       write(*,*)'iTable = ',iTable,' should be 1'
       call CON_stop(NameSub)
    end if
    Ptr=>Table_I(1)

    if(iProc==0) write(*,*)'testing make_lookup_table'
    call make_lookup_table(1, eos_rho_e, MPI_COMM_WORLD)
    
    if(iProc==0) write(*,*)'testing interpolate_lookup_table'    
    call interpolate_lookup_table(1, 1.0, 2.0, p_I)
    ! rho=1.0 is exactly in the middle, e=2.0 is also an index, so exact result

    pGood_I = (/ 4./3., 4./5., 3. /)
    if(any(abs(p_I - pGood_I) > 1e-5))then
       write(*,*)'p_I=',p_I,' is different from pGood_I=',pGood_I
       call CON_stop(NameSub)
    end if

    if(iProc==0) write(*,*)'testing load_lookup_table'

    ! Load the saved file into the second table
    call init_lookup_table(&
         NameTable   = "RhoE2"                 ,&
         NameCommand = "load"                  ,&
         NameFile    = "test_lookup_table1.out",&
         TypeFile    = "ascii")

    

    if(iProc==0) write(*,*)'testing i_lookup_table for table 2'
    iTable = i_lookup_table("RhoE2")
    if(iTable /= 2)then
       write(*,*)'iTable = ',iTable,' should be 2'
       call CON_stop(NameSub)
    end if
    Ptr2=>Table_I(iTable)

    if(Ptr2%StringDescription /= Ptr%StringDescription) &
         call CON_stop(NameSub // &
         ' Description='//trim(Ptr2%StringDescription)// &
         ' is different from '// trim(Ptr%StringDescription))

    if(Ptr2%NameVar /= Ptr%NameVar) call CON_stop(NameSub // &
         ' NameVar='//trim(Ptr2%NameVar)//' is different from '// &
         trim(Ptr2%NameVar))

    if(Ptr%nValue /= Ptr2%nValue)then
       write(*,*)'nValue=',Ptr2%nValue,' is different from ',Ptr%nValue
       call CON_stop(NameSub)
    end if

    if(iProc==0) write(*,*)'testing interpolate_lookup_table for loaded table'
    call interpolate_lookup_table(2, 1.0, 2.0, p_I)
    if(any(abs(p_I - pGood_I) > 1e-5))then
       write(*,*)'p_I=',p_I,' is different from pGood_I=',pGood_I
       call CON_stop(NameSub)
    end if
    
    call interpolate_lookup_table(2, 3, 3.0, 2.0, p_I, Arg) 
    if(any(abs(p_I - pGood_I) > 1e-5) .or. abs(Arg - 1.0) >  1e-5)then
       write(*,*)'p_I=',p_I, ' Arg =', Arg,&
            ' are different from pGood_I=',pGood_I, ' ArgGood = 1.0'
       call CON_stop(NameSub)
    end if

  end subroutine test_lookup_table

  !===========================================================================

  subroutine eos_rho_e(iTable, rho, e, p_I)
    ! This is an example for the subroutine passed to make_lookup_table

    integer, intent(in):: iTable
    real, intent(in)   :: rho, e
    real, intent(out)  :: p_I(:)
    !-----------------------------------------------------------------------
    p_I(1) = (2./3.)*e
    p_I(2) = (2./5.)*e
    p_I(3) = e + rho

  end subroutine eos_rho_e

end module ModLookupTable
