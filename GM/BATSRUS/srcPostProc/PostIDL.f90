!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================================
program PostIDL

  ! Read a .h file from STDIN, then read all the files based on the .h file
  ! and put structured Coord_DC and State_VC together and save into a VAC file

  use ModPlotFile, ONLY: save_plot_file
  use ModNumConst, ONLY: cHalfPi
  use ModInterpolate, ONLY: linear
  use ModNumConst,       ONLY: cTwoPi, cRadToDeg
  use ModCoordTransform, ONLY: rot_matrix_z
  use ModUtilities,      ONLY: lower_case, split_string, join_string
  use ModSort,     ONLY: sort_quick
  use ModReadParam
  
  implicit none

  ! This is copied from ModKind, because PostIDL.exe may be compiled with
  ! different precision then the rest of the codes.
  integer, parameter :: Real4_=selected_real_kind(6,30)
  integer, parameter :: Real8_=selected_real_kind(12,100)
  integer, parameter :: nByteReal = 4 + (1.00000000041 - 1.0)*10000000000.0

  ! Global variables

  integer, parameter :: unit_tmp=99
  character(len=20) :: TypeFile='real4'
  !TypeFile = 'ascii', 'real8', 'real4'

  integer :: nxyz(3), icell, countcell, ncell, nw, nEqPar, nProc, it
  real :: t
  real, allocatable :: Coord_DC(:,:,:,:), State_VC(:,:,:,:)
  real, allocatable :: GenCoord_DI(:,:) ! gen. coords for unstructured grids
  real, allocatable :: w1(:), EqPar_I(:)
  real, allocatable :: Param_I(:)

  real(Real4_)              :: DxCell4, Xyz4_D(3)
  real(Real4_), allocatable :: State4_V(:)
  real(Real8_)              :: DxCell8, Xyz8_D(3)
  real(Real8_), allocatable :: State8_V(:)

  ! Coordinates, sizes, indices
  real, dimension(3) :: Xyz_D, xyzmin, xyzmax
  real, dimension(3) :: CellSizePlot_D, dxyz, dxyzmin
  real, dimension(3) :: XyzGen_D
  real ::    x, y, z, xmin, ymin, zmin
  real ::    dx, dy, dz, dyperdx, dzperdx, dxcell, dycell, dzcell
  real ::    frac
  real(selected_real_kind(12))  :: total, volume
  integer :: i, j, k, imin, imax, jmin, jmax, kmin, kmax, nx, ny, nz, iw

  ! Variables related to sorting and averaging unstructured data
  integer:: nSum
  integer, allocatable :: iSort_I(:)
  real, allocatable :: Sort_I(:)
  real, allocatable :: StateSum_V(:), CellSizeMin_D(:)

  integer :: idim, icutdim(3), nDim, nspecialpar
  real    :: specialpar(3)
  character(len=5):: NameCoord_D(3) = (/'x    ','y    ','z    '/)
  character(len=5):: NameCoordPlot_D(3)

  logical :: IsStructured, DoReadBinary=.false.
  character (len=100) :: filename, filenamehead, coordnames
  character (len=500) :: varnames, unitnames
  integer :: l, me

  ! Variables for the 2D lookup table
  integer :: nx1,nx2,idim1,idim2
  integer :: idim0 ! the ignored dimension
  integer :: iError
  real    :: xmin1, xmax1, xmin2, xmax2, dx1, dx2

  ! Variables for checking binary compatibility
  integer            :: nByteRealRead

  ! Variables for generalized coordinates
  character (len=79) :: TypeGeometry='cartesian'
  logical            :: UseDoubleCut = .false.

  ! Toroidal geometry
  integer:: iPoint, nPoint
  real :: rTorusSmall, rTorusLarge
  real, allocatable:: TorusSurface_I(:)

  ! Logarithmic radial coordinate
  logical:: IsLogRadius = .false.

  ! Stretched radial coordinates
  logical:: IsGenRadius = .false.
  integer:: nRgen = -1
  real, allocatable:: LogRgen_I(:)

  ! Roundcube parameters
  real:: rRound0    ! fully Cartesian distance
  real:: rRound1    ! fully round distance
  real:: SqrtNDim   ! sqrt 2 or sqrt 3

  integer :: nDimSim ! Dimension of simulation domain.

  ! Parameters read from head file, but not used for this subroutine.
  real:: cLight, ThetaTilt, rBody
  logical:: IsPeriodic_D(3)
  
  character (len=lStringLine) :: NameCommand, StringLine
  !---------------------------------------------------------------------------

  write(*,'(a)')'PostIDL (G.Toth 2000-) starting'


  call read_file()
  call read_init()

  
  call read_echo_set(.true.)

  READPARAM: do
     if(.not.read_line(StringLine) )then
        EXIT READPARAM
     end if

     if(.not.read_command(NameCommand)) CYCLE READPARAM

     select case(NameCommand)
     case('#HEADFILE')
        call read_var('HeadFileName', filenamehead)
        call read_var('nProc',nProc)
        call read_var('SaveBinary', DoReadBinary)
        nByteRealRead = -1
        if(DoReadBinary) then
           call read_var('nByteReal', nByteRealRead)
           if(nByteRealRead==nByteReal)then
              write(*,*)'nByteReal=',nByteReal
           else if(nByteRealRead==-1)then
              write(*,*)'!!! Warning: PostIDL was compiled with ',&
                   nByteReal,' byte reals but nByteReal is not given in file !!'
           else if(nByteRealRead < nByteReal)then
              write(*,*)'!!! Warning: PostIDL was compiled with ',&
                   nByteReal,' byte reals but file contains nByteReal=', &
                   nByteRealRead
           end if
        endif
        
        ! Get rid of the directory part
        filenamehead = filenamehead( &
             index(filenamehead,'/',BACK=.true.)+1:len(filenamehead))

     case('#NDIM')
        call read_var('nDimSim',nDimSim)

     case('#NSTEP')
        call read_var('nStep',it)

     case('#TIMESIMULATION')
        call read_var('TimeSimulation',t)

     case('#PLOTRANGE')
        xyzmin = 0; xyzmax = 0        
        do i = 1, nDimSim
           call read_var('CoordMin',xyzmin(i))
           call read_var('CoordMax',xyzmax(i))
        enddo

     case('#PLOTRESOLUTION')
        do i = 1, nDimSim
           call read_var('DxSavePlot',CellSizePlot_D(i))
        enddo

     case('#CELLSIZE')
        dxyzmin = 1
        do i = 1, nDimSim
           call read_var('CellSizeMin',dxyzmin(i))
        enddo

     case('#NCELL')
        call read_var('nCellPlot',ncell)

     case('#PLOTVARIABLE')
        call read_var('nPlotVar', nw)
        call read_var('VarNames',varnames)
        call read_var('Unit',unitnames)

     case('#SCALARPARAM')
        call read_var('nParam', nEqPar)
        allocate(EqPar_I(nEqPar))
        do i = 1, nEqPar
           call read_var('Param',EqPar_I(i))
        enddo
        call read_var('cLight',cLight)
        call read_var('ThetaTild',ThetaTilt)
        call read_var('rBody',rBody)

     case('#GRIDGEOMETRYLIMIT')
        call read_var('TypeGeometry', TypeGeometry)
        IsLogRadius = index(TypeGeometry,'lnr')  > 0
        IsGenRadius = index(TypeGeometry,'genr') > 0
        if(IsGenRadius)then
           read(*,*) nRgen
           allocate(LogRgen_I(nRgen))
           do i = 1, nRgen
              read(*,*) LogRgen_I(i)
           end do
        end if

        if(TypeGeometry == 'roundcube')then
           call read_var('rRound0', rRound0)
           call read_var('rRound1', rRound1)
           call read_var('SqrtNDim',SqrtNDim)
        endif

     case('#PERIODIC')
        do i = 1, nDimSim
           call read_var('IsPeriodic',IsPeriodic_D(i))
        enddo
        
     case('#OUTPUTFORMAT')
        call read_var('OutPutFormat',TypeFile)
        
     case default
        write(*,*) 'WARNING: unknow command ', NameCommand
     end select

  enddo READPARAM
  
  if(filenamehead(1:3) == 'sph')then
     NameCoord_D    = (/'r    ','theta','phi  '/)
  elseif(TypeGeometry == 'rz' .or. TypeGeometry == 'xr')then
     TypeGeometry = 'cartesian'
     NameCoord_D    = (/'x    ','r    ','phi  '/)
  elseif(filenamehead(1:3) == 'cut')then
     select case(TypeGeometry(1:5))
     case('spher')
        ! In reality this is the r-lon-lat coordinate system
        NameCoord_D = (/'r    ','phi  ','lat  '/)
     case('cylin')
        NameCoord_D = (/'r    ','phi  ','z    '/)
     case('round')
        NameCoord_D = (/'xi   ','eta  ','zeta '/)
     end select
  end if

  ! Save input CellSizePlot_D into dxyz that may get overwritten
  dxyz = CellSizePlot_D

  ! Unstructured grid has negative dx
  IsStructured = dxyz(1) >= 0.0

  ! If dx<=0. use the smallest cell as resolution
  if(dxyz(1)<1.e-6) dxyz = dxyzmin

  ! Calculate structured grid size
  nxyz = max(1, nint((xyzmax - xyzmin)/dxyz))

  write(*,*)'plot area size=', nxyz

  ! Calculate dimensionality of the cut and add specialparameters if needed
  nDim=0
  nspecialpar=0
  icutdim=0
  do i = 1, 3
     if(nxyz(i)>1)then
        nDim=nDim+1
        icutdim(nDim)=i
     else
        icutdim(3)=i
        nspecialpar = nspecialpar + 1
        specialpar(nspecialpar) = 0.5*(xyzmax(i) + xyzmin(i))
        varnames=trim(varnames)//' cut'//trim(NameCoord_D(i))
     end if
  end do
  
  if(nDim==2)then
     ! Make a lookup table to check coinciding cells
     idim1=icutdim(1)
     idim2=icutdim(2)
     idim0=icutdim(3)
     xmin1=xyzmin(idim1)
     xmin2=xyzmin(idim2)
     xmax1=xyzmax(idim1)
     xmax2=xyzmax(idim2)
     dx1=dxyzmin(idim1)  ! Note that we use smallest cell size
     dx2=dxyzmin(idim2)
     nx1=nint((xmax1-xmin1)/dx1)
     nx2=nint((xmax2-xmin2)/dx2)

     ! Sph/cyl. X=0 and Y=0 cuts require doubled lookup table (+/- r)
     if(idim0==2)then
        if(TypeGeometry(1:9)=='spherical' .and. xmax2 > cHalfPi) then
           ! Use LatMin < Lat' < 2*LatMax-LatMin as generalized coordinate
           UseDoubleCut = .true.; nx2 = 2*nx2; 
           ! nxyz(3) = 2*nxyz(3)
        elseif(TypeGeometry=='cylindrical' .and. xmax2 > cHalfPi)then
           ! Use rMin < r' < 2*rMax - rMin as generalized coordinate
           UseDoubleCut = .true.; nx1 = 2*nx1; 
           ! nxyz(1) = 2*nxyz(1)
        end if

        ! Do not attempt to use structured grid for double cuts
        if(UseDoubleCut)then
           IsStructured = .false.
           CellSizePlot_D = -1.0
        end if
     end if

  endif

  ! For unstructured grid make the Coord_DC and State_VC arrays linear
  if(.not.IsStructured)then
     nxyz(1)=ncell
     nxyz(2:3)=1
  end if

  ! Get components for sake of efficiency
  xmin=xyzmin(1); ymin=xyzmin(2); zmin=xyzmin(3)
  dx = dxyz(1);   dy = dxyz(2);   dz = dxyz(3)
  nx = nxyz(1);   ny = nxyz(2);   nz = nxyz(3)

  ! Cell aspect ratios
  dyperdx=dxyzmin(2)/dxyzmin(1); dzperdx=dxyzmin(3)/dxyzmin(1)

  ! Allocate State_VC and Coord_DC, the arrays of variables and coordinates
  allocate(w1(nw),State_VC(nw,nx,ny,nz), Coord_DC(nDim,nx,ny,nz), STAT=iError)
  if(iError /= 0) stop 'PostIDL.exe ERROR: could not allocate arrays'

  if(.not.IsStructured)then     
     allocate(GenCoord_DI(nDim,nCell), STAT=iError)
     if(iError /= 0) stop 'PostIDL.exe ERROR: could not allocate enCoord_DI'
  end if

  if(DoReadBinary.and.nByteRealRead==4) allocate(State4_V(nw))
  if(DoReadBinary.and.nByteRealRead==8) allocate(State8_V(nw))

  
  ! Initialize State_VC and Coord_DC
  State_VC = 0.0
  Coord_DC = 0.0

  call set_strings

  ! Collect info from all files and put it into State_VC and Coord_DC
  total=0.0
  icell=0
  countcell=0
  l=len_trim(filenamehead)
  do me=0,nProc-1
     if(    nProc > 100000)then
        write(filename,'(a,i6.6,a)')filenamehead(1:l-2)//"_pe",me,'.idl'
     elseif(nProc > 10000)then
        write(filename,'(a,i5.5,a)')filenamehead(1:l-2)//"_pe",me,'.idl'
     else
        write(filename,'(a,i4.4,a)')filenamehead(1:l-2)//"_pe",me,'.idl'
     end if

     if(me==0)write(*,*)'reading files=',trim(filename),&
          '...',nProc-1,'.idl'

     if(DoReadBinary)then
        open(unit_tmp, file=filename, status='old', form='unformatted', &
             iostat=iError)
     else
        open(unit_tmp, file=filename, status='old', iostat=iError)
     end if

     ! Assume that missing files were empty. 
     if(iError /=0) CYCLE

     ! Read file
     do
        !Debug
        !write(*,*)'START READING'
        if(DoReadBinary)then
           if(nByteRealRead == 4)then
              read(unit_tmp,ERR=999,END=999) DxCell4, Xyz4_D, State4_V
              DxCell = DxCell4; Xyz_D = Xyz4_D; w1 = State4_V
           else
              read(unit_tmp,ERR=999,END=999) DxCell8, Xyz8_D, State8_V
              DxCell = DxCell8; Xyz_D = Xyz8_D; w1 = State8_V
           end if
        else
           read(unit_tmp,*,ERR=999,END=999) DxCell, Xyz_D, w1
        end if

        countcell = countcell + 1
        dycell=dxcell*dyperdx; dzcell=dxcell*dzperdx

        ! Set XyzGen_D, possibly modify Xyz_D too
        call set_gen_coord

        if(.not.IsStructured)then

           ! Simply put data into array. 
           ! Sorting and averaging will be done at the end.

           iCell = iCell + 1
           State_VC(:,iCell,1,1) = w1
           do iDim = 1, nDim
              Coord_DC(iDim,iCell,1,1) = Xyz_D(icutdim(iDim))
              GenCoord_DI(iDim,iCell)  = XyzGen_D(iCutDim(iDim))
           end do

           ! We are finished with unstructured
           CYCLE

        endif

        x = XyzGen_D(1); y = XyzGen_D(2); z = XyzGen_D(3)

        if(dxcell < dx+1.e-6)then
           ! Cell has the correct size or finer
           i = max(1, nint((x-xmin)/dx+0.5))
           j = max(1, nint((y-ymin)/dy+0.5))
           k = max(1, nint((z-zmin)/dz+0.5))

           if(dxcell<dx-1.e-6)then
              ! Cell is finer, calculate volume fraction
              frac=(dxcell/dx)**ndim
           else
              frac=1.0
           end if

           State_VC(:,i,j,k) = State_VC(:,i,j,k) + frac*w1
           do iDim = 1, nDim
              Coord_DC(iDim,i,j,k) = Coord_DC(iDim,i,j,k) &
                   + frac*Xyz_D(iCutDim(iDim))
           end do
           total=total + frac
        else
           ! Cell is coarser than required resolution
           imin=min(nx,max(1,nint((x-0.5*dxcell-xmin)/dx+1)))
           imax=min(nx,max(1,nint((x+0.5*dxcell-xmin)/dx)))
           jmin=min(ny,max(1,nint((y-0.5*dycell-ymin)/dy+1)))
           jmax=min(ny,max(1,nint((y+0.5*dycell-ymin)/dy)))
           kmin=min(nz,max(1,nint((z-0.5*dzcell-zmin)/dz+1)))
           kmax=min(nz,max(1,nint((z+0.5*dzcell-zmin)/dz)))

           ! First order prolongation
           do iw = 1, nw
              State_VC(iw,imin:imax,jmin:jmax,kmin:kmax)= &
                   State_VC(iw,imin:imax,jmin:jmax,kmin:kmax) &
                   + w1(iw)
           end do
           do iDim = 1, nDim
              Coord_DC(iDim,imin:imax,jmin:jmax,kmin:kmax) = &
                   Coord_DC(iDim,imin:imax,jmin:jmax,kmin:kmax) &
                   + Xyz_D(iCutDim(iDim))
           end do

           if(imax<imin.or.jmax<jmin.or.kmax<kmin)&
                write(*,*)'!!! Empty box for cell dx,x,y,z=',dxcell,x,y,z

           total = total + (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
        end if
     end do ! read file

999  continue

     close(unit_tmp)
  end do ! me
  
  if(countcell/=ncell)&
       write(*,*)'!!! Discrepancy: countcell=',countcell,' ncell=',ncell,' !!!'

  if(IsStructured)then
     volume=product(real(nxyz))
     if(ndim==1 .and. abs(total/volume-4.0)<0.0001)then
        State_VC = 0.25*State_VC
        Coord_DC = 0.25*Coord_DC
        write(*,*)'Averaged 1D structured file everywhere'
     elseif(abs(total/volume-2.0)<0.0001)then
        State_VC = 0.5*State_VC
        Coord_DC = 0.5*Coord_DC
        write(*,*)'Averaged structured file everywhere'
     elseif(abs(total/volume-1.0)>0.0001)then
        write(*,*)'!!! Discrepancy in structured file:',&
             'filled total=',total,' volume=',volume,' !!!'
     end if
  else
     if(icell /= ncell) &
          write(*,*)'!!! Error: ncell=',ncell,' /= icell=',icell,' !!!'

     nx=icell
     nxyz(1)=icell
  end if

  if(filenamehead(1:3) == 'cut')then

     ! Convert radians to degrees
     do iDim = 1, nDim
        select case(NameCoordPlot_D(iDim))
        case('r')
           ! Convert generalized coordinates to radius and degrees
           if(IsLogRadius)then
              Coord_DC(iDim,1:nx,:,:) = exp(Coord_DC(iDim,1:nx,:,:))
           elseif(IsGenRadius)then
              do k = 1, nz; do j = 1, ny; do i = 1, nx
                 Coord_DC(iDim,i,j,k) = exp(linear(LogRgen_I, 0, nRgen-1, &
                      Coord_DC(iDim,i,j,k)*(nRgen-1), DoExtrapolate=.true.) )
              end do; end do; end do
           end if

        case('phi', 'theta', 'lat')
           Coord_DC(iDim,1:nx,:,:) = Coord_DC(iDim,1:nx,:,:)*cRadToDeg
        end select
     end do

  endif

  filename = filenamehead(1:l-2)//'.out'
  write(*,*)'writing file =',trim(filename)

  ! Param_I is the combination of eqpar and specialpar
  allocate(Param_I(nEqPar+nspecialpar))
  do i = 1, nEqPar
     Param_I(i)=EqPar_I(i)
  end do
  do i = 1, nSpecialPar
     Param_I(i+nEqPar) = SpecialPar(i)
  end do
  
  if(.not.IsStructured)then
     ! Sort points based on (generalized) coordinates

     allocate(Sort_I(nx), iSort_I(nx), STAT=iError)
     if(iError /= 0) stop 'PostIDL.exe ERROR: could not allocate sort arrays'

     ! Form sorting function from the generalized coordinates
     Sort_I = GenCoord_DI(1,:)
     if(nDim > 1) Sort_I = Sort_I + exp(1.0)*GenCoord_DI(2,:)
     if(nDim > 2) Sort_I = Sort_I + exp(2.0)*GenCoord_DI(3,:)

     ! Sort points according to the sorting function
     call sort_quick(nx, Sort_I, iSort_I)
     Coord_DC(:,:,1,1) = Coord_DC(:,iSort_I,1,1)
     State_VC(:,:,1,1) = State_VC(:,iSort_I,1,1)

     ! Average out coinciding points
     if(nDim < 3) then
        GenCoord_DI = GenCoord_DI(:,iSort_I)

        allocate(StateSum_V(nw), CellSizeMin_D(nDim))
        CellSizeMin_D = dxyzmin(icutdim(1:nDim))
        i = 1
        k = 1
        do while(i < nx)
           StateSum_V = State_VC(:,i,1,1)
           nSum       = 1
           j = i + 1
           do while( sum(abs(GenCoord_DI(:,j) - GenCoord_DI(:,i)) &
                /CellSizeMin_D) < 0.01)
              StateSum_V = StateSum_V + State_VC(:,j,1,1)
              nSum = nSum + 1
              j = j + 1
              if(j > nx) EXIT
           end do
           if(j > i+1) then
              ! Put average value into i-th element
              State_VC(:,i,1,1) = StateSum_V/nSum              
           end if
           ! Save the index for the unique coordinates 
           iSort_I(k) = i
           k = k + 1
           i = j 
        end do
        deallocate(StateSum_V, CellSizeMin_D, GenCoord_DI)

        ! Special Judgement for the last point
        if(j == nx) then
           iSort_I(k) = j
           nx = k
        else
           nx = k - 1
        end if
        
        ! move the elements after finding out all the coinciding ones 
        Coord_DC(:,1:nx,1,1) = Coord_DC(:,iSort_I(1:nx),1,1)
        State_VC(:,1:nx,1,1) = State_VC(:,iSort_I(1:nx),1,1)
     
     end if
  
     deallocate(Sort_I, iSort_I)
  end if

  ! the sizes of Coord_DC and State_VC may be modified by cell averaging 
  ! in unstructured grids. Only the first dimension (1:nx) needs to be set
  call save_plot_file(filename,&
       TypeFileIn = TypeFile, &
       StringHeaderIn = unitnames, &
       nStepIn = it, TimeIn = t, &
       ParamIn_I = Param_I, &
       NameVarIn = varnames, &
       IsCartesianIn = TypeGeometry=='cartesian' .and. IsStructured,&
       nDimIn = nDim,&
       CoordIn_DIII = Coord_DC(:,1:nx,:,:), & 
       VarIn_VIII = State_VC(:,1:nx,:,:))

  deallocate(Coord_DC, State_VC, Param_I)
  deallocate(w1, EqPar_I)
  if(allocated(State4_V)) deallocate(State4_V)
  if(allocated(State8_V)) deallocate(State8_V)

  write(*,'(a)')'PostIDL finished'

contains

  !===========================================================================

  subroutine set_strings

    ! Produce coordinate names 
    !         ('x y z', 'x y', 'x z', 'y z' or 'r theta', 'r phi' ...)

    ! Reorder coordinate names for cuts
    NameCoordPlot_D(1:nDim) = NameCoord_D(iCutDim(1:nDim))

    ! Fix first coordinate name for non-cartesian cut along phi=90,270 deg
    ! that is a cut in the SECOND generalized coordinate.
    if(filenamehead(1:3) == 'x=0') NameCoordPlot_D(1) = 'y'

    call join_string(NameCoordPlot_D(1:nDim), coordnames)

    varnames = trim(coordnames)//' '//trim(varnames)

  end subroutine set_strings

  !===========================================================================

  subroutine set_gen_coord

    ! Calculate the generalized coordinates XyzGen_D from Xyz_D
    ! mostly for lookup

    real:: rCyl ! distance from axis Z

    ! Toroidal variables
    real:: PoloidalAngle, r, z, StretchCoef, dAngle, Residual, WallRadius

    ! Rotation matrix for rotated Cartesian grid
    real, allocatable, save:: GridRot_DD(:,:)

    ! List of vector variables for rotated Cartesian coordinates
    integer:: nVector = 0
    integer, allocatable, save:: iVarVector_I(:)

    ! Temporary variables to process the variable name string
    integer:: iVector, iVar, MaxWord, nWord, l
    character(len=10), allocatable:: NameVar_V(:)
    character(len=10)             :: NameVar, NameVector

    ! Variables for roundcube coordinate transformation
    real :: r2, Dist1, Dist2, Coef1, Coef2
    !---------------------------------------------------------------------
    if(TypeGeometry == 'cartesian' .or. filenamehead(1:3) == 'cut')then
       XyzGen_D = Xyz_D
       
       RETURN
    end if

    if(TypeGeometry == 'roundcube')then

       r2 = sum(Xyz_D**2)
       if (r2 > 0.0) then
          ! L1 and L2 distance
          Dist1 = maxval(abs(Xyz_D))
          Dist2 = sqrt(r2)
          if (rRound1 > rRound0 ) then
             ! The rounded (distorted) grid is outside of the non-distorted part
             if (Dist1 > rRound0) then
                ! Outside the undistorted region
                ! Assume Coord = w * Xyz and Replace Xyz in transformation Coord_to_Xyz
                ! We have to solve a quadratic equation of w. 
                ! w^2 - Coef1*w- Coef2 = 0
                Coef1 = -1 + rRound0/(rRound1-rRound0)*(dist1*SqrtNDim/Dist2 - 1)
                Coef2 = Dist1/(rRound1-rRound0)*(SqrtNDim*Dist1/Dist2 - 1)
                XyzGen_D = Xyz_D/(-Coef1 + sqrt(Coef1**2 + 4*Coef2))*2
             else
                ! No distortion
                XyzGen_D = Xyz_D
             end if
         
          else
             ! The rounded (distorted) grid is inside of the non-distorted part
             if (Dist2 < rRound1) then
                ! Solving w^2 - w + Coef1 = 0
                Coef1 = Dist1/rRound1*(1 - Dist1/Dist2)
                XyzGen_D = Xyz_D / (1 + sqrt(1-4*Coef1))*2
             else
                ! Solving w^2 + Coef1*w + Coef2 = 0
                Coef1 = -1 + (1 - Dist1/Dist2)/(rRound0-rRound1)*rRound0
                Coef2 = -(1 - Dist1/Dist2)/(rRound0 - rRound1)*Dist1
                Coef2 = (-Coef1 + sqrt(Coef1**2 - 4*Coef2))*0.5
                XyzGen_D = Xyz_D / Coef2
             end if
          end if
           
       else
          XyzGen_D = 0.0
       end if
       RETURN
    end if

    if(TypeGeometry(1:7)=='rotated')then

       if(.not.allocated(GridRot_DD))then
          ! Setup rotation matrix
          allocate(GridRot_DD(3,3))

          ! The rotation matrix should be the same as in BATL_geometry
          GridRot_DD = rot_matrix_z(0.6,0.8)

          ! Find vectors
          MaxWord = nDim + nW + nEqPar
          allocate(NameVar_V(MaxWord))
          call split_string(varnames, MaxWord, NameVar_V, nWord)

          ! Make this array large enough
          allocate(iVarVector_I(nW/3))

          NameVector = ' '
          do iVar = 1, nW - 2
             ! Add nDim to skip the coordinate names
             NameVar = NameVar_V(nDim+iVar)
             call lower_case(NameVar)

             l = len_trim(NameVar)

             ! Identify vectors as 3 strings ending with x, y, z
             if(NameVar(l:l) /= 'x') CYCLE

             ! Prospective vector component
             NameVector = NameVar(1:l-1)

             ! Check the next two names
             NameVar = NameVar_V(nDim+iVar+1)
             call lower_case(NameVar)
             if(NameVar /= trim(NameVector)//'y') CYCLE

             NameVar = NameVar_V(nDim+iVar+2)
             call lower_case(NameVar)
             if(NameVar /= trim(NameVector)//'z') CYCLE

             nVector = nVector + 1
             iVarVector_I(nVector) = iVar
          end do
          deallocate(NameVar_V)

          !write(*,*)'nVector, iVarVector_I=', nVector, iVarVector_I(1:nVector)

       end if

       ! Unrotate the coordinates for comparison with Cartesian runs
       Xyz_D = matmul(Xyz_D, GridRot_DD)
       XyzGen_D = Xyz_D

       ! Unrotate vectors
       do iVector = 1, nVector
          iVar = iVarVector_I(iVector)
          w1(iVar:iVar+2) = matmul(w1(iVar:iVar+2), GridRot_DD)
       end do

       RETURN
    end if

    rCyl = sqrt(Xyz_D(1)**2 + Xyz_D(2)**2)

    ! Calculate phi
    if (rCyl == 0.0) then
       XyzGen_D(2) = 0.0
    else
       XyzGen_D(2) = &
            modulo(atan2(Xyz_D(2), Xyz_D(1)) - XyzMin(2), cTwoPi) + XyzMin(2)
    end if

    select case(TypeGeometry)
    case('cylindrical', 'cylindrical_lnr', 'cylindrical_genr')
       XyzGen_D(1) = rCyl
       XyzGen_D(3) = Xyz_D(3)

       if(nDim==2)then
          ! Set the 'X-Y' coordinates for plotting a 2D cut
          select case(idim0)
          case(1)
             ! This is R=const slice, use longitude [deg] vs height
             Xyz_D(2)=XyzGen_D(2)*cRadToDeg
          case(2)
             ! This is x=0 or y=0 plane, use signed radius vs Z
             Xyz_D(1) = sign(1.0, Xyz_D(1)+Xyz_D(2))*rCyl
             ! Radial distance
             XyzGen_D(1) = rCyl
             ! The generalized coordinate runs from rMin to 2*rMax-rMin
             if(UseDoubleCut .and. Xyz_D(1) < 0.0) &
                  XyzGen_D(1) = XyzGen_D(1) + xmax1 - xmin1
          end select
       end if
    case('spherical', 'spherical_lnr', 'spherical_genr')
       XyzGen_D(1) = sqrt(rCyl**2 + Xyz_D(3)**2)

       if(nDim==2)then
          ! Set the 'X-Y' coordinates for plotting a 2D cut
          select case(idim0)
          case(1)
             ! This is R=const slice, use longitude vs latitude in degs.
             Xyz_D(2:3)=XyzGen_D(2:3)*cRadToDeg
          case(2)
             ! This is x=0 or y=0 plane, use axial radius vs Z
             Xyz_D(1) = sign(1.00,Xyz_D(1)+Xyz_D(2))*rCyl
             XyzGen_D(2) = xmin2 ! could be useful for structured grid
          case(3)
             ! This is the z=0 plane
             ! Stretch X and Y with rSph/rCyl instead of simply
             ! projecting the points down to the X-Y plane
             Xyz_D(1:2) = Xyz_D(1:2)*XyzGen_D(1)/rCyl
          end select
       end if
       ! Latitude
       XyzGen_D(3) = asin(Xyz_D(3)/XyzGen_D(1))
       ! Shift by width of latitude range for the left half 
       if(UseDoubleCut .and. Xyz_D(1) < 0.0) & 
            XyzGen_D(3) = XyzGen_D(3) + xmax2 - xmin2

    case('axial_torus')

       !This is x=0 or y=0 plane
       if(nDim==2 .and. idim0==2) Xyz_D(1) = rCyl

       r = rCyl - rTorusLarge
       z = Xyz_D(3)
       if(.not.(z==0.0 .and. r==0.0))then
          PoloidalAngle = modulo(atan2(Z,R), cTwoPi)

          ! Use linear interpolation to obtain the distance to the wall
          dAngle = cTwoPi/nPoint
          iPoint = modulo(int(PoloidalAngle/dAngle), nPoint)
          Residual = PoloidalAngle - iPoint*dAngle
          WallRadius = &
               (1-Residual)*TorusSurface_I(iPoint  ) &
               +  Residual*TorusSurface_I(iPoint+1)

          StretchCoef = rTorusSmall/ &
               (WallRadius * max(abs(r), abs(z))/sqrt(r**2 + z**2))

          r = r*StretchCoef
          z = z*StretchCoef
       end if
       XyzGen_D(3) = z
       XyzGen_D(1) = r + rTorusLarge

    case default
       write(*,*)'Unknown TypeGeometry='//TypeGeometry
       stop
    end select

    if(IsLogRadius .or. IsGenRadius) XyzGen_D(1) = log(XyzGen_D(1))

    if(IsGenRadius)then
       i = min(nRgen-1, count(LogRgen_I < XyzGen_D(1)))
       XyzGen_D(1) = ( i -1  &
            + (XyzGen_D(1) - LogRgen_I(i)) / (LogRgen_I(i+1) - LogRgen_I(i)) )&
            / (nRgen - 1)
    end if

  end subroutine set_gen_coord

end program PostIDL

!=============================================================================

subroutine CON_stop(String)

  ! This routine is needed for ModPlotFile

  implicit none

  character(len=*), intent(in):: String
  write(*,*) 'ERROR in PostIDL: '//String
  stop

end subroutine CON_stop

!=============================================================================

double precision function MPI_WTIME()

  ! This is needed for ModUtilities.F90

  integer:: clock,clockrate,count_max

  call system_clock(clock,clockrate,count_max)
  MPI_WTIME=dble(clock)/clockrate

end function MPI_WTIME

!=============================================================================

subroutine MPI_comm_rank(iComm, iProc, iError)
  integer, intent(in)  :: iComm
  integer, intent(out) :: iProc, iError
  iProc = 0
  iError = 0
  write(*,*) 'I am MPI_comm_rank'
end subroutine MPI_comm_rank

!=============================================================================

subroutine MPI_Bcast(Buffer, Count, Datatype, Root, iComm, iError)
  integer, intent(in) :: Buffer ! Do not allow changing buffer!!
  integer, intent(in) :: Count
  integer, intent(in) :: Datatype
  integer, intent(in) :: Root
  integer, intent(in) :: iComm
  integer, intent(out) :: iError

  write(*,*) 'I am MPI_Bcast'
  
end subroutine MPI_Bcast

!=============================================================================

! subroutine MPI_Bcast(Buffer, Count, Datatype, Root, iComm, iError)
!   character(len=*), intent(in) :: Buffer(:) ! Do not allow changing buffer!!
!   integer, intent(in) :: Count
!   integer, intent(in) :: Datatype
!   integer, intent(in) :: Root
!   integer, intent(in) :: iComm
!   integer, intent(out) :: iError

!   write(*,*) 'I am MPI_Bcast'
  
! end subroutine MPI_Bcast
