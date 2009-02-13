module ModPlotFile

  ! Save or read VAC/IDL type plotfiles from 1 up to 3 dimensions.
  ! ASCII, single (real4), or double (real8) precision binary file formats 
  ! can be used.
  ! The plot file contains 5 header lines: 
  !
  !    Header
  !    nStep, Time, nDim, nParam, nVar
  !    n1 .. nNDim
  !    Param1 .. ParamNParam
  !    NameVar
  !
  ! Header   (string) describes the plot. It is up to 500 characters.
  ! nStep    (integer) is the number of time steps/iterations etc.
  ! Time     (real) is the simulation time.
  ! nDim     (integer) number of dimensions. Negative for non-Cartesian grids.
  ! nParam   (integer) number of parameters (for example adiabatic index)
  ! n1 ..    (nDim integers) grid sizes in the nDim dimensions
  ! Param1.. (nParam reals) parameter values
  ! NameVar  (string) space separated list of names for the 
  !                   nDim coordinates, nVar variables and nParam parameters
  !
  ! The header is followed by the coordinate/variable data. 
  ! In ASCII files each line contains the coordinates+variables for one cell.
  ! The cells are ordered such that the first coordinate index changes fastest.
  ! In binary files the coordinates are saved in a single array, followed by
  ! the variables, saved as one record per variable.
  !
  ! For save_plot_file the coordinates can be either given as full arrays,
  ! or for Cartesian grids they can be given as 1D arrays for each dimension,
  ! or for uniform Cartesian grids they can be given with min and max values.
  ! The number of dimensions, the size of the grid and the number of the 
  ! variables is determined from the size of the variable array.
  !
  ! For read_plot_file the number of dimensions, variables, parameters, grid
  ! size are optionaal parameters.

  use ModIoUnit,    ONLY: UnitTmp_
  use ModKind,  ONLY: Real4_

  implicit none

  private ! except

  public:: read_plot_file
  public:: save_plot_file
  public:: test_plot_file

  integer, parameter :: MaxDim = 3

contains

  !=========================================================================

  subroutine save_plot_file(NameFile, TypeFileIn, &
       StringHeaderIn, nStepIn, TimeIn, &
       ParamIn_I, NameVarIn, &
       IsCartesianIn, &
       nDimIn,&
       CoordMinIn_D, CoordMaxIn_D, &
       Coord1In_I, Coord2In_I, Coord3In_I, &
       CoordIn_I, CoordIn_DII, CoordIn_DIII, &
       VarIn_VI, VarIn_VII, VarIn_VIII)

    character(len=*),           intent(in):: NameFile       ! Name of plot file
    character(len=*), optional, intent(in):: TypeFileIn     ! ascii/real8/real4
    character(len=*), optional, intent(in):: StringHeaderIn ! header line
    integer,          optional, intent(in):: nStepIn        ! number of steps
    real,             optional, intent(in):: TimeIn         ! simulation time  
    real,             optional, intent(in):: ParamIn_I(:)   ! parameters
    character(len=*), optional, intent(in):: NameVarIn      ! list of names
    logical,          optional, intent(in):: IsCartesianIn  ! Cartesian grid?
    integer,          optional, intent(in):: nDimIn         ! grid dimensions
    real,             optional, intent(in):: CoordIn_I(:)   ! coords in 1D
    real,             optional, intent(in):: CoordIn_DII(:,:,:)       ! 2D
    real,             optional, intent(in):: CoordIn_DIII(:,:,:,:)    ! 3D
    real,             optional, intent(in):: Coord1In_I(:)  ! coords for axis 1
    real,             optional, intent(in):: Coord2In_I(:)  ! coords for axis 2
    real,             optional, intent(in):: Coord3In_I(:)  ! coords for axis 3
    real,             optional, intent(in):: CoordMinIn_D(:)! min coordinates
    real,             optional, intent(in):: CoordMaxIn_D(:)! max coordinates
    real,             optional, intent(in):: VarIn_VI(:,:)  ! variables in 1D
    real,             optional, intent(in):: VarIn_VII(:,:,:)            ! 2D
    real,             optional, intent(in):: VarIn_VIII(:,:,:,:)         ! 3D

    character(len=20)  :: TypeFile
    character(len=500) :: StringHeader
    character(len=500) :: NameVar
    integer            :: nStep, nDim, nParam, nVar, n1, n2, n3
    real               :: Time, Coord
    logical            :: IsCartesian
    real,         allocatable :: Param_I(:), Coord_ID(:,:), Var_IV(:,:)

    integer :: n_D(0:MaxDim)
    integer :: i, j, k, i_D(3), iDim, iVar, n, nDimOut, iError

    character(len=*), parameter:: NameSub = 'save_plot_file'
    !---------------------------------------------------------------------

    TypeFile = 'ascii'
    if(present(TypeFileIn)) TypeFile = TypeFileIn

    StringHeader = 'No header info'
    if(present(StringHeaderIn)) StringHeader = StringHeaderIn

    nStep = 0
    if(present(nStepIn)) nStep = nStepIn

    Time = 0.0
    if(present(TimeIn)) Time = TimeIn

    if(present(ParamIn_I))then
       nParam = size(ParamIn_I)
       allocate(Param_I(nParam))
       Param_I = ParamIn_I
    else
       nParam = 1
       allocate(Param_I(1))
       Param_I(1) = 0.0
    end if

    ! Figure out grid dimensions and number of variables
    n_D = 1
    if(present(VarIn_VI))then
       nDim = 1
       n_D(0:1) = shape(VarIn_VI)
    elseif(present(VarIn_VII))then
       nDim = 2
       n_D(0:2) = shape(VarIn_VII)
    elseif(present(VarIn_VIII))then
       nDim = 3
       n_D(0:3) = shape(VarIn_VIII)
    else
       call CON_stop(NameSub // &
            ' none of VarIn_VI, VarIn_VII, VarIn_VIII are present')
    endif

    ! Extract information
    nVar = n_D(0)
    n1   = n_D(1)
    n2   = n_D(2)
    n3   = n_D(3)

    ! The plot dimension may be different from the dimensionality of VarIn
    if(present(nDimIn))then
       nDim = nDimIn
       if(n1 == 1 .and. n2 == 1)then
          n_D(1:3) = (/ n3, 1, 1/)
       elseif(n1 == 1)then
          n_D(1:3) = (/ n2, n3, 1/)
       elseif(n2 == 1)then
          n_D(1:3) = (/ n1, n3, 1/)
       end if
    end if

    IsCartesian = .true.
    if(present(IsCartesianIn)) IsCartesian = IsCartesianIn

    ! nDim is saved with a negative sign for non-Cartesian grid
    nDimOut = nDim
    if(.not. IsCartesian) nDimOut = -nDim
    if(present(NameVarIn))then
       NameVar = NameVarIn
    else
       ! Create some arbitrary variable names
       NameVar = 'x1'
       do i = 2, nDim
          write(NameVar, "(a, i1)") trim(NameVar) // ' x', i
       end do
       do i = 1, nVar
          write(NameVar, "(a, i2.2)") trim(NameVar) // ' v', i
       end do
       do i = 1, nParam
          write(NameVar, "(a, i2.2)") trim(NameVar) // ' p', i
       end do
    end if
    ! Allocate arrays with a shape that is convenient for saving data
    allocate(Coord_ID(n1*n2*n3,nDim), Var_IV(n1*n2*n3,nVar))

    ! Fill in the 2D coordinate array using the available information
    do iDim = 1, nDim
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1, n1
          n = n + 1
          Coord = huge(1.0)
          if(present(CoordIn_I))    Coord = CoordIn_I(i)
          if(present(CoordIn_DII))  Coord = CoordIn_DII(iDim,i,j)
          if(present(CoordIn_DIII)) Coord = CoordIn_DIII(iDim,i,j,k)
          if(present(Coord1In_I) .and. iDim==1) Coord = Coord1In_I(i)
          if(present(Coord2In_I) .and. iDim==2) Coord = Coord2In_I(j)
          if(present(Coord3In_I) .and. iDim==3) Coord = Coord3In_I(k)
          if(present(CoordMinIn_D)) then
             i_D = (/i, j, k/)
             Coord = CoordMinIn_D(iDim) + (i_D(iDim)-1)* &
                  ((CoordMaxIn_D(iDim) - CoordMinIn_D(iDim))/(n_D(iDim)-1))
          end if
          Coord_ID(n, iDim) = Coord
       end do; end do; end do; 
    end do

    ! Check if all coordinates were set
    if(any(Coord_ID == huge(1.0))) call CON_stop(NameSub // & 
        ' coordinates were not defined')

    ! Fill in the 2D variable array using the available information
    Var_IV = huge(1.0)
    do iVar = 1, nVar
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1,n1
          n = n + 1
          if(present(VarIn_VI))   Var_IV(n,iVar) = VarIn_VI(iVar,i)
          if(present(VarIn_VII))  Var_IV(n,iVar) = VarIn_VII(iVar,i,j)
          if(present(VarIn_VIII)) Var_IV(n,iVar) = VarIn_VIII(iVar,i,j,k)
       end do; end do; end do; 
    end do
   
    ! Check if all variables were set
    if(any(Var_IV == huge(1.0))) call CON_stop(NameSub // & 
         ' variables were not defined')

    select case(TypeFile)
    case('formatted', 'ascii')
       open(UnitTmp_, file=NameFile, iostat=iError)
       if(iError /= 0)call CON_stop(NameSub // &
            ' could not open ascii file=' // trim(NameFile))

       write(UnitTmp_, "(a)")             trim(StringHeader)
       write(UnitTmp_, "(i7,es13.5,3i3)") nStep, Time, nDimOut, nParam, nVar
       write(UnitTmp_, "(3i8)")           n_D(1:nDim)
       write(UnitTmp_, "(100es13.5)")     Param_I
       write(UnitTmp_, "(a)")             trim(NameVar)

       ! write out coordinates and variables line by line
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1, n1
          n = n + 1
          write(UnitTmp_, "(100es18.10)") Coord_ID(n,:), Var_IV(n, :) 
       end do; end do; end do

    case('real8')
       open(UnitTmp_, file=NameFile, form='unformatted', iostat=iError)
       if(iError /= 0)call CON_stop(NameSub // &
            ' could not open real8 file=' // trim(NameFile))
       write(UnitTmp_) StringHeader
       write(UnitTmp_) nStep, Time, nDimOut, nParam, nVar
       write(UnitTmp_) n_D(1:nDim)
       write(UnitTmp_) Param_I
       write(UnitTmp_) NameVar
       write(UnitTmp_) Coord_ID
       ! write out variables 1 by 1 to avoid segmentation fault 
       ! for very large Var_IV array
       do iVar = 1, nVar
          write(UnitTmp_) Var_IV(:,iVar)
       end do
    case('real4')
       open(UnitTmp_, file=NameFile, form='unformatted', iostat=iError)
       if(iError /= 0)call CON_stop(NameSub // &
            ' could not open real4 file=' // trim(NameFile))
       write(UnitTmp_) StringHeader
       write(UnitTmp_) nStep, real(Time, Real4_), nDimOut, nParam, nVar
       write(UnitTmp_) n_D(1:nDim)
       write(UnitTmp_) real(Param_I, Real4_)
       write(UnitTmp_) NameVar
       write(UnitTmp_) real(Coord_ID, Real4_)
       ! write out variables 1 by 1 to avoid segmentation fault 
       ! for very large Var_IV array
       do iVar = 1, nVar
          write(UnitTmp_) real(Var_IV(:,iVar), Real4_)
       end do

    case default
       call CON_stop(NameSub // ' unknown TypeFile =' // trim(TypeFile))
    end select
    close(UnitTmp_)

    deallocate(Param_I, Coord_ID, Var_IV)

  end subroutine save_plot_file

  !=========================================================================

  subroutine read_plot_file(NameFile, TypeFileIn, &
       StringHeaderOut, &
       nStepOut, TimeOut, nDimOut, nParamOut, nVarOut, &
       IsCartesianOut, &
       n1Out, n2Out, n3Out, nOut_D, &
       ParamOut_I, NameVarOut, &
       CoordMinOut_D, CoordMaxOut_D, &
       Coord1Out_I, Coord2Out_I, Coord3Out_I, &
       CoordOut_I, CoordOut_DII, CoordOut_DIII, &
       VarOut_VI, VarOut_VII, VarOut_VIII)

    character(len=*),            intent(in) :: NameFile
    character(len=*), optional, intent(in) :: TypeFileIn
    character(len=*), optional, intent(out):: StringHeaderOut
    character(len=*), optional, intent(out):: NameVarOut
    real,             optional, intent(out):: TimeOut  
    integer,          optional, intent(out):: nStepOut
    integer,          optional, intent(out):: nDimOut   ! number of dimensions
    integer,          optional, intent(out):: nParamOut ! number of parameters
    integer,          optional, intent(out):: nVarOut   ! number of variables
    integer,          optional, intent(out):: n1Out, n2Out, n3Out ! grid size
    integer,          optional, intent(out):: nOut_D(:) ! grid size array
    logical,          optional, intent(out):: IsCartesianOut ! Cartesian grid?
    real,             optional, intent(out):: ParamOut_I(:)  ! parameters
    real,             optional, intent(out):: CoordMinOut_D(:)
    real,             optional, intent(out):: CoordMaxOut_D(:)
    real,             optional, intent(out):: CoordOut_I(:)
    real,             optional, intent(out):: Coord1Out_I(:)
    real,             optional, intent(out):: Coord2Out_I(:)
    real,             optional, intent(out):: Coord3Out_I(:)
    real,             optional, intent(out):: CoordOut_DII(:,:,:)
    real,             optional, intent(out):: CoordOut_DIII(:,:,:,:)
    real,             optional, intent(out):: VarOut_VI(:,:)
    real,             optional, intent(out):: VarOut_VII(:,:,:)
    real,             optional, intent(out):: VarOut_VIII(:,:,:,:)

    character(len=20)  :: TypeFile
    character(len=500) :: StringHeader
    character(len=500) :: NameVar
    integer            :: nStep, nDim, nParam, nVar, n1, n2, n3, n_D(MaxDim)
    real               :: Time, Coord
    real(Real4_)       :: Time4
    logical            :: IsCartesian
    real(Real4_), allocatable:: Param4_I(:), Coord4_ID(:,:), Var4_IV(:,:)
    real,         allocatable:: Param_I(:),  Coord_ID(:,:),  Var_IV(:,:)

    integer :: i, j, k, iDim, iVar, n, iError

    character(len=*), parameter:: NameSub = 'read_plot_file'
    !---------------------------------------------------------------------

    TypeFile = 'ascii'
    if(present(TypeFileIn)) TypeFile = TypeFileIn

    n_D = 1
    select case(TypeFile)
    case('ascii', 'formatted')
       open(UnitTmp_, file=NameFile, status='old', iostat=iError)
       if(iError /= 0) call CON_stop(NameSub // &
            ' could not open ascii file=' // trim(NameFile))

       read(UnitTmp_, '(a)') StringHeader
       read(UnitTmp_, *) nStep, Time, nDim, nParam, nVar
       read(UnitTmp_, *) n_D(1:abs(nDim))
       allocate(Param_I(nParam))
       read(UnitTmp_, *) Param_I
       read(UnitTmp_, '(a)') NameVar

    case('real8')
       open(UnitTmp_, file=NameFile, status='old', form='unformatted', &
            iostat=iError)
       if(iError /= 0) call CON_stop(NameSub // &
            ' could not open ascii file=' // trim(NameFile))

       read(UnitTmp_) StringHeader
       read(UnitTmp_) nStep, Time, nDim, nParam, nVar
       read(UnitTmp_) n_D(1:abs(nDim))
       allocate(Param_I(nParam))
       read(UnitTmp_) Param_I
       read(UnitTmp_) NameVar
    case('real4')
       open(UnitTmp_, file=NameFile, status='old', form='unformatted', &
            iostat=iError)
       if(iError /= 0) call CON_stop(NameSub // &
            ' could not open ascii file=' // trim(NameFile))

       read(UnitTmp_) StringHeader
       read(UnitTmp_) nStep, Time4, nDim, nParam, nVar
       Time = Time4
       read(UnitTmp_) n_D(1:abs(nDim))
       allocate(Param_I(nParam), Param4_I(nParam))
       read(UnitTmp_) Param4_I
       Param_I = Param4_I
       deallocate(Param4_I)
       read(UnitTmp_) NameVar
    case default
       call CON_stop(NameSub // ' unknown TypeFile =' // trim(TypeFile))
    end select

    IsCartesian = nDim > 0
    nDim = abs(nDim)
    n1 = n_D(1); n2 = n_D(2); n3 = n_D(3) 
   
    ! Read coordinates and variables into suitable 2D arrays
    allocate(Coord_ID(n1*n2*n3, nDim), Var_IV(n1*n2*n3, nVar))
    select case(TypeFile)
    case('ascii', 'formatted')
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1, n1
          n = n + 1
          read(UnitTmp_, *) Coord_ID(n, :), Var_IV(n, :)
       end do; end do; end do

    case('real8')
       read(UnitTmp_) Coord_ID
       do iVar = 1, nVar
          read(UnitTmp_) Var_IV(:, iVar)
       end do

    case('real4')
       allocate(Coord4_ID(n1*n2*n3, nDim), Var4_IV(n1*n2*n3, nVar))
       read(UnitTmp_) Coord4_ID
       Coord_ID = Coord4_ID
       do iVar = 1, nVar
          read(UnitTmp_) Var4_IV(:, iVar)
       end do
       Var_IV = Var4_IV
       deallocate(Coord4_ID, Var4_IV)
    end select
    close(UnitTmp_)

    if(present(StringHeaderOut)) StringHeaderOut = trim(StringHeader)
    if(present(NameVarOut))      NameVarOut      = trim(NameVar)
    if(present(TimeOut))         TimeOut         = Time
    if(present(nStepOut))        nStepOut        = nStep
    if(present(nDimOut))         nDimOut         = nDim
    if(present(nParamOut))       nParamOut       = nParam
    if(present(nVarOut))         nVarOut         = nVar
    if(present(n1Out))           n1Out           = n1
    if(present(n2Out))           n2Out           = n2
    if(present(n3Out))           n3Out           = n3
    if(present(nOut_D))          nOut_D(1:nDim)  = n_D(1:nDim)
    if(present(IsCartesianOut))  IsCartesianOut  = IsCartesian
    if(present(ParamOut_I))      ParamOut_I(1:nParam) = Param_I
    if(present(CoordMinOut_D)) CoordMinOut_D(1:nDim) = minval(Coord_ID, DIM=1)
    if(present(CoordMaxOut_D)) CoordMaxOut_D(1:nDim) = maxval(Coord_ID, DIM=1)

    ! Fill in output coordinate arrays
    do iDim = 1, nDim
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1, n1
          n = n + 1
          Coord = Coord_ID(n, iDim)
          if(present(CoordOut_I))    CoordOut_I(i)             = Coord
          if(present(CoordOut_DII))  CoordOut_DII(iDim,i,j)    = Coord
          if(present(CoordOut_DIII)) CoordOut_DIII(iDim,i,j,k) = Coord
          if(present(Coord1Out_I) .and. iDim==1 .and. j==1 .and. k==1) &
               Coord1Out_I(i) = Coord
          if(present(Coord2Out_I) .and. iDim==2 .and. i==1 .and. k==1) &
               Coord2Out_I(j) = Coord
          if(present(Coord3Out_I) .and. iDim==3 .and. i==1 .and. j==1) &
               Coord3Out_I(k) = Coord
       end do; end do; end do
    end do

    ! Fill in output variable arrays
    do iVar = 1, nVar
       n = 0
       do k = 1, n3; do j = 1, n2; do i = 1, n1
          n = n + 1
          if(present(VarOut_VI))   VarOut_VI(iVar, i)      = Var_IV(n, iVar)
          if(present(VarOut_VII))  VarOut_VII(iVar,i,j)    = Var_IV(n, iVar)
          if(present(VarOut_VIII)) VarOut_VIII(iVar,i,j,k) = Var_IV(n, iVar)
       end do; end do; end do
    end do

    deallocate(Param_I, Coord_ID, Var_IV)

  end subroutine read_plot_file

  !=========================================================================

  subroutine test_plot_file

    ! Set up a hydro shock tube initial condition on a 2D Cartesian grid
    ! Save plot file then read it and check consistency
    ! Do this multiple times with various settings

    character(len=*), parameter:: StringHeaderIn = "test_hd22"
    real,    parameter :: TimeIn = 25.0
    integer, parameter :: nStepIn = 10, nDimIn = 2, nParamIn = 2, nVarIn = 4
    integer, parameter :: n1In= 10, n2In = 2
    real,    parameter :: CoordMinIn_D(nDimIn) = (/ 0.5, -0.5 /)
    real,    parameter :: CoordMaxIn_D(nDimIn) = (/ 9.5,  0.5 /)
    real,    parameter :: ParamIn_I(nParamIn) = (/ 1.667, 2.5 /)
    character(len=*), parameter:: NameVarIn = "x y rho ux uy p gamma rbody"
    real    :: CoordIn_DII(nDimIn, n1In, n2In), VarIn_VII(nVarIn, n1In, n2In)
    real    :: CoordIn_DIII(nDimIn, n1In, 1, n2In), VarIn_VIII(nVarIn, n1In, 1, n2In)

    ! Do tests with ascii/real8/real4 files, 
    ! Cartesian/non-Cartesian coordinates
    ! 2D/3D input arrays
    integer, parameter:: nTest = 12
    character(len=5)  :: TypeFileIn_I(nTest) = &
         (/ 'ascii', 'real8', 'real4', 'ascii', 'real8', 'real4', &
            'ascii', 'real8', 'real4', 'ascii', 'real8', 'real4' /)
    logical           :: IsCartesianIn_I(nTest) = &
         (/ .true.,   .true., .true.,  .false.,   .false., .false.,&
            .true.,   .true., .true.,  .false.,   .false., .false. /)

    ! Input and output of tests
    character(len=80)    :: NameFile
    character(len=20)    :: TypeFileIn
    character(len=100)   :: StringHeaderOut
    real                 :: TimeOut
    integer              :: nStepOut, nDimOut, nParamOut, nVarOut
    integer              :: n1Out, n2Out
    integer              :: nOut_D(3)
    real                 :: ParamOut_I(100)
    character(len=100)   :: NameVarOut
    logical              :: IsCartesianIn, IsCartesianOut
    real                 :: CoordMinOut_D(nDimIn), CoordMaxOut_D(nDimIn)
    real                 :: Coord1Out_I(n1In), Coord2Out_I(n2In)
    real                 :: CoordOut_DII(nDimIn, n1In, n2In)
    real                 :: CoordOut_DIII(nDimIn, n1In, n2In, 1)
    real                 :: VarOut_VII(nVarIn, n1In, n2In)
    real                 :: VarOut_VIII(nVarIn, n1In, n2In, 1)

    ! Tolerance for errors
    real :: Eps

    ! Indexes
    integer :: i, j, iTest

    character(len=*), parameter:: NameSub = 'test_plot_file'
    !----------------------------------------------------------------------

    ! Initialize coordinates an variables: shock tube on a 2D uniform grid
    do j = 1, n2In; do i = 1, n1In
       CoordIn_DII(1, i, j) = CoordMinIn_D(1) &
            + (i-1)*((CoordMaxIn_D(1)-CoordMinIn_D(1))/(n1In - 1))
       CoordIn_DII(2, i, j) = CoordMinIn_D(2) &
            + (j-1)*((CoordMaxIn_D(2)-CoordMinIn_D(2))/(n2In - 1))
       CoordIn_DIII(1, i, 1, j) = CoordIn_DII(1,i,j)
       CoordIn_DIII(2, i, 1, j) = CoordIn_DII(2,i,j)

       if(i <= n1In/2)then
          VarIn_VII(:, i, j) = (/ 1.0, 0.0, 0.0, 1.0 /)
          VarIn_VIII(:, i, 1, j) = (/ 1.0, 0.0, 0.0, 1.0 /)
       else
          VarIn_VII(:,i,j) = (/ 0.1, 0.0, 0.0, 0.125 /)
          VarIn_VIII(:, i, 1, j) = (/ 0.1, 0.0, 0.0, 0.125 /)  
       end if
    end do; end do

    ! Test ascii, real8 and real4 files
    do iTest = 1, nTest 
       write(NameFile, '(a,i2.2,a)') 'test_plot_file',iTest,'.out'
       write(*,*) NameSub, ' writing file=', trim(NameFile)

       TypeFileIn    = TypeFileIn_I(iTest)
       IsCartesianIn = IsCartesianIn_I(iTest)

       if(TypeFileIn == 'real4')then
          Eps = 1e-5
       else
          Eps = 1e-12
       end if

       ! Test saving it
       select case(iTest)
       case(1)
          ! Use coordinate ranges
          call save_plot_file(NameFile,      &
               TypeFileIn     = TypeFileIn,     &
               StringHeaderIn = StringHeaderIn, &
               nStepIn        = nStepIn,        &
               TimeIn         = TimeIn,         &
               ParamIn_I      = ParamIn_I,      &
               NameVarIn      = NameVarIn,      &
               IsCartesianIn  = IsCartesianIn,  &
               CoordMinIn_D   = CoordMinIn_D,   &
               CoordMaxIn_D   = CoordMaxIn_D,   &
               VarIn_VII      = VarIn_VII)
       case(2)
          ! Use 1D coordinate arrays
          call save_plot_file(NameFile,      &
               TypeFileIn     = TypeFileIn,     &
               StringHeaderIn = StringHeaderIn, &
               nStepIn        = nStepIn,        &
               TimeIn         = TimeIn,         &
               ParamIn_I      = ParamIn_I,      &
               NameVarIn      = NameVarIn,      &
               IsCartesianIn  = IsCartesianIn,  &
               Coord1In_I     = CoordIn_DII(1,:,1), &
               Coord2In_I     = CoordIn_DII(2,1,:), &
               VarIn_VII      = VarIn_VII)
       case(3:6)
          ! Use full coordinate array
          call save_plot_file(NameFile,      &
               TypeFileIn     = TypeFileIn,     &
               StringHeaderIn = StringHeaderIn, &
               nStepIn        = nStepIn,        &
               TimeIn         = TimeIn,         &
               ParamIn_I      = ParamIn_I,      &
               NameVarIn      = NameVarIn,      &
               IsCartesianIn  = IsCartesianIn,  &
               CoordIn_DII    = CoordIn_DII,    &
               VarIn_VII      = VarIn_VII)
       case default
          ! Test 3D input array
          ! Use full coordinate array
          call save_plot_file(NameFile,      &
               TypeFileIn     = TypeFileIn,     &
               StringHeaderIn = StringHeaderIn, &
               nStepIn        = nStepIn,        &
               TimeIn         = TimeIn,         &
               ParamIn_I      = ParamIn_I,      &
               NameVarIn      = NameVarIn,      &
               IsCartesianIn  = IsCartesianIn,  &
               nDimIn         = nDimIn,         &
               CoordIn_DIII   = CoordIn_DIII,   &
               VarIn_VIII     = VarIn_VIII)
       
       end select

          call read_plot_file(NameFile,        &
            TypeFileIn      = TypeFileIn,      &
            StringHeaderOut = StringHeaderOut, &
            nStepOut        = nStepOut,        &
            TimeOut         = TimeOut,         &
            nDimOut         = nDimOut,         &
            nParamOut       = nParamOut,       &
            nVarOut         = nVarOut,         &
            ParamOut_I      = ParamOut_I,      &
            NameVarOut      = NameVarOut,      &
            IsCartesianOut  = IsCartesianOut,  &
            CoordOut_DII    = CoordOut_DII,    &
            Coord1Out_I     = Coord1Out_I,     &
            Coord2Out_I     = Coord2Out_I,     &
            CoordMinOut_D   = CoordMinOut_D,   &
            CoordMaxOut_D   = CoordMaxOut_D,   &
            VarOut_VII      = VarOut_VII)



       if(nStepOut /= nStepIn)then
          write(*,*)'nStepIn=', nStepIn,' nStepOut=', nStepOut
          call CON_stop(NameSub)
       end if

       if(abs(TimeOut - TimeIn) > Eps)then
          write(*,*)'TimeIn=', TimeIn,' TimeOut=', TimeOut
          call CON_stop(NameSub)
       end if

       if(nDimOut /= nDimIn)then
          write(*,*)'nDimIn=', nDimIn,' nDimOut=', nDimOut
          call CON_stop(NameSub)
       end if

       if(nParamOut /= nParamIn)then
          write(*,*)'nParamIn=', nParamIn,' nParamOut=', nParamOut
          call CON_stop(NameSub)
       end if

       if(nVarOut /= nVarIn)then
          write(*,*)'nVarIn=', nVarIn,' nVarOut=', nVarOut
          call CON_stop(NameSub)
       end if

       if(any(abs(ParamOut_I(1:nParamIn) - ParamIn_I) > Eps))then
          write(*,*)'ParamIn=', ParamIn_I,' ParamOut=', ParamOut_I(1:nParamIn)
          call CON_stop(NameSub)
       end if
       
       if(IsCartesianOut .neqv. IsCartesianIn)then
          write(*,*)'IsCartesianIn, Out=', IsCartesianIn, IsCartesianOut
          call CON_stop(NameSub)
       end if
 
       if(NameVarOut /= NameVarIn)then
             write(*,*)'NameVarIn=', NameVarIn,' NameVarOut=', NameVarOut
             call CON_stop(NameSub)
          end if
       
       !To simpify, replace the 3D input array with 2D 
       if(iTest > 6)then
          CoordIn_DII = CoordIn_DIII(:,:,1,:)
          VarIn_VII = VarIn_VIII(:,:,1,:)
       end if 
       do j = 1, n2In; do i = 1, n1In
          if(any(abs(CoordIn_DII(:,i,j) - CoordOut_DII(:,i,j)) > Eps))then
              write(*,*)'i,j=', i, j
              write(*,*)'CoordIn =', CoordIn_DII(:,i,j)
              write(*,*)'CoordOut=', CoordOut_DII(:,i,j)
              call CON_stop(NameSub)
           end if
           if(abs(CoordIn_DII(1,i,j) - Coord1Out_I(i)) > Eps )then
              write(*,*)'i,j=', i, j
              write(*,*)'CoordIn(1)=', CoordIn_DII(1,i,j)
              write(*,*)'Coord1Out =', Coord1Out_I(i)
              call CON_stop(NameSub)
           end if
           if(abs(CoordIn_DII(2,i,j) - Coord2Out_I(j)) > Eps )then
              write(*,*)'i,j=', i, j
              write(*,*)'CoordIn(2)=', CoordIn_DII(2,i,j)
              write(*,*)'Coord2Out =', Coord2Out_I(j)
              call CON_stop(NameSub)
           end if

           if(any(abs(VarIn_VII(:,i,j) - VarOut_VII(:,i,j)) > Eps))then
              write(*,*)'i,j=', i, j
              write(*,*)'VarIn =', VarIn_VII(:,i,j)
              write(*,*)'VarOut=', VarOut_VII(:,i,j)
              call CON_stop(NameSub)
           end if
        end do; end do

        if(abs(CoordMinOut_D(1) -  minval(CoordIn_DII(1,:,:))) >Eps)then
           write(*,*)'CoordMinOut_D(1)     =',CoordMinOut_D(1)
           write(*,*)'minval(CoordIn_DII(1)=',minval(CoordIn_DII(1,:,:))
           call CON_stop(NameSub)
        end if
        if(abs(CoordMinOut_D(2) -  minval(CoordIn_DII(2,:,:))) > Eps)then
           write(*,*)'CoordMinOut_D(2)     =',CoordMinOut_D(2)
           write(*,*)'minval(CoordIn_DII(2)=',minval(CoordIn_DII(2,:,:))
           call CON_stop(NameSub)
        end if
        if(abs(CoordMaxOut_D(1) -  maxval(CoordIn_DII(1,:,:))) > Eps )then
           write(*,*)'CoordMaxOut_D(1)     =',CoordMaxOut_D(1)
           write(*,*)'maxval(CoordIn_DII(1)=',maxval(CoordIn_DII(1,:,:))
           call CON_stop(NameSub)
        end if
        if(abs(CoordMaxOut_D(2) - maxval(CoordIn_DII(2,:,:))) >Eps)then
           write(*,*)'CoordMaxOut_D(2)     =',CoordMaxOut_D(2)
           write(*,*)'maxval(CoordIn_DII(2)=',maxval(CoordIn_DII(2,:,:))
           call CON_stop(NameSub)
        end if
      
    end do

    ! Test using defaults for 2D input array
    NameFile = 'test_plot_file13.out'       
    call save_plot_file(NameFile, VarIn_VII=VarIn_VII, CoordIn_DII=CoordIn_DII)

    ! Read header info
    call read_plot_file(NameFile, &
         StringHeaderOut=StringHeaderOut, NameVarOut=NameVarOut, &         
         nDimOut=nDimOut, nVarOut=nVarOut, nParamOut=nParamOut, &
         IsCartesianOut=IsCartesianOut, nOut_D=nOut_D)

    if( StringHeaderOut /= 'No header info')call CON_stop(NameSub // &
         ' incorrect value for default StringHeaderOut='//StringHeaderOut)
    if( NameVarOut /= 'x1 x2 v01 v02 v03 v04 p01') call CON_stop(NameSub // &
         ' incorrect value for default NameVarOut='//NameVarOut)
    if(.not.IsCartesianOut) call CON_stop(NameSub // &
         ' incorrect value for IsCartesianOut: should be true')
    if( any(nOut_D(1:nDimOut) /= (/ n1In, n2In /)) ) then
       write(*,*) 'n1In, n2In=', n1In, n2In
       write(*,*) 'nOut_D    =',nOut_D
       call CON_stop(NameSub)
    end if

    ! Test using defaults for 3D input array
    NameFile = 'test_plot_file14.out'       
    call save_plot_file(NameFile,nDimIn = nDimIn, VarIn_VIII = VarIn_VIII,&
                        CoordIn_DIII = CoordIn_DIII)

    ! Read header info
    call read_plot_file(NameFile, &
         StringHeaderOut = StringHeaderOut, NameVarOut = NameVarOut, &         
         nDimOut = nDimOut, nVarOut = nVarOut, nParamOut = nParamOut, &
         IsCartesianOut = IsCartesianOut, nOut_D = nOut_D)

    if( StringHeaderOut /= 'No header info')call CON_stop(NameSub // &
         ' incorrect value for default StringHeaderOut='//StringHeaderOut)
    if( NameVarOut /= 'x1 x2 v01 v02 v03 v04 p01') call CON_stop(NameSub // &
         ' incorrect value for default NameVarOut='//NameVarOut)
    if(.not.IsCartesianOut) call CON_stop(NameSub // &
         ' incorrect value for IsCartesianOut: should be true')
    if( any(nOut_D(1:nDimOut) /= (/ n1In, n2In /)) ) then
       write(*,*) 'n1In, n2In=', n1In, n2In
       write(*,*) 'nOut_D    =',nOut_D
       call CON_stop(NameSub)
    end if
    

    ! Now that we have the dimensions, we could allocate coordinate and
    ! variable arrays and read them

  end subroutine test_plot_file

end module ModPlotFile
