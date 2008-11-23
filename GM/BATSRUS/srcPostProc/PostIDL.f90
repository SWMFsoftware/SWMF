!^CFG COPYRIGHT UM
!=============================================================================
program PostIDL

  ! Read a .h file from STDIN, then read all the files based on the .h file
  ! and put structured x and w together and save into a VAC file

  implicit none

  ! This is copied from ModKind, because PostIDL.exe may be compiled with
  ! different precision then the rest of the codes.
  integer, parameter :: Real4_=selected_real_kind(6,30)
  integer, parameter :: Real8_=selected_real_kind(12,100)
  integer, parameter :: nByteReal = 4 + (1.00000000041 - 1.0)*10000000000.0


  ! Global variables

  integer, parameter :: unit_tmp=99
  logical, parameter :: write_binary=.true.
  real, parameter :: halfeps=0.6

  integer :: nxyz(3), icell, countcell, ncell, nw, neqpar, numprocs, it
  real :: t
  real, dimension(:,:,:,:), allocatable :: xx, w
  real, dimension(:), allocatable :: w1, eqpar, dxdoubled

  real(Real4_)              :: DxCell4, Xyz4_D(3)
  real(Real4_), allocatable :: State4_V(:)
  real(Real8_)              :: DxCell8, Xyz8_D(3)
  real(Real8_), allocatable :: State8_V(:)

  ! Coordinates, sizes, indices
  real, dimension(3) :: Xyz_D, xyzmin, xyzmax, dxyz, dxyzmin, dxyzcell
  real, dimension(3) :: XyzGen_D
  real ::    x, y, z, xmin, ymin, zmin
  real ::    dx, dy, dz, dyperdx, dzperdx, dxcell, dycell, dzcell
  real ::    frac
  real(selected_real_kind(12))  :: total, volume
  integer :: i, j, k, imin, imax, jmin, jmax, kmin, kmax, nx, ny, nz, iw

  integer :: idim, icutdim(3), ndim, nspecialpar
  real    :: xcut(3), specialpar(3)
  character (LEN=5), dimension(3)            :: coord
  character (LEN=5), dimension(3), parameter :: &
       coord_xyz=(/'x    ','y    ','z    '/), &
       coord_sph=(/'r    ','theta','phi  '/)     !Theta is colatitude

  logical :: structured, read_binary=.false., UseLookup=.false.
  character (len=100) :: filename, filenamehead, coordnames
  character (len=500) :: varnames, unitnames, fileheadout
  integer :: l, ll, me

  ! Variables for the 2D lookup table
  integer :: ix1,ix2,ixmin1,ixmax1,ixmin2,ixmax2,nx1,nx2,idim1,idim2,jcell=0
  integer :: idim0 ! the ignored dimension
  integer :: iError
  real    :: xmin1, xmax1, xmin2, xmax2, dx1, dx2, dx1cell, dx2cell
  integer, dimension(:,:), allocatable :: lookup

  ! Variables for checking binary compatibility
  integer            :: nByteRealRead

  ! Variables for generalized coordinates
  character (len=79) :: TypeGeometry='cartesian', TypeGeometryRead
  logical            :: UseDoubleCut = .false.

  !Toroidal geometry
  integer:: iPoint, nPoint
  real :: rTorusSmall, rTorusLarge
  real, allocatable:: TorusSurface_I(:)
  !---------------------------------------------------------------------------

  write(*,'(a)')'PostIDL (G.Toth 2000-2002) starting'

  ! Read information from STDIN
  read(*,'(a)')filenamehead

  ! Get rid of the directory part
  filenamehead = filenamehead( &
       index(filenamehead,'/',BACK=.true.)+1:len(filenamehead))

  read(*,*)numprocs
  write(*,*)trim(filenamehead),', numprocs=',numprocs
  if(filenamehead(1:2) == 'sp')then
     coord=coord_sph
  else
     coord=coord_xyz
  end if
  read(*,*)it
  read(*,*)t
  write(*,*)'n_step=',it,' time_simulation=',t
  read(*,*)(xyzmin(i),xyzmax(i),i=1,3)
  write(*,*)'xyzmin=',xyzmin
  write(*,*)'xyzmax=',xyzmax
  read(*,*)dxyz,dxyzmin,ncell
  write(*,*)'dxyz,dxyzmin,ncell=',dxyz,dxyzmin,ncell
  read(*,*)nw
  read(*,*)neqpar
  allocate(eqpar(neqpar))
  read(*,*)eqpar
  write(*,*)'nvar=',nw,' neqpar=',neqpar,' eqpar=',eqpar
  read(*,'(a)')varnames
  write(*,*)'varnames =',trim(varnames)
  read(*,'(a)')unitnames
  write(*,*)'unitnames=',trim(unitnames)
  if(unitnames=='')unitnames='normalized units'

  read_binary = .false.
  read(*,'(l8)',err=1,end=1)read_binary
1 continue
  write(*,*)'binary   =',read_binary

  if(read_binary)then
     nByteRealRead = -1
     read(*,'(i8)',err=2,end=2)nByteRealRead
2    continue
     if(nByteRealRead==nByteReal)then
        write(*,*)'nByteReal=',nByteReal
     else if(nByteRealRead==-1)then
        write(*,*)'!!! Warning: PostIDL was compiled with ',&
             nByteReal,' byte reals but nByteReal is not given in file !!!'
     else if(nByteRealRead < nByteReal)then
        write(*,*)'!!! Warning: PostIDL was compiled with ',&
             nByteReal,' byte reals but file contains nByteReal=',nByteRealRead
     end if
     write(*,*)'nByteReal=',nByteRealRead
  end if

  !Read TypeGeometry, if possible
  read(*,'(a)',err=3,end=3) TypeGeometryRead
  TypeGeometry = TypeGeometryRead
3 continue
  write(*,*)'TypeGeometry = ',trim(TypeGeometry)
  if(index(TypeGeometry,'torus')>0)then
     open(unit_tmp,file='torus.dat',status='old')
     read(unit_tmp,*)nPoint, rTorusSmall, rTorusLarge
     allocate(TorusSurface_I(0:nPoint))
     do i=0,nPoint
        read(unit_tmp,*)iPoint,TorusSurface_I(iPoint)
     end do
  end if
  ! Unstructured grid has dx=-1.
  structured = dxyz(1) > -0.9

  ! If dx<=0. use the smallest cell as resolution
  if(dxyz(1)<1.e-6)dxyz=dxyzmin

  ! Calculate structured grid size
  nxyz=max(1,nint((xyzmax-xyzmin)/dxyz))

  write(*,*)'plot area size=', nxyz

  ! Calculate dimensionality of the cut and add specialparameters if needed
  ndim=0
  nspecialpar=0
  icutdim=0
  do i=1,3
     if(nxyz(i)>1)then
        ndim=ndim+1
        icutdim(ndim)=i
     else
        icutdim(3)=i
        nspecialpar=nspecialpar+1
        specialpar(nspecialpar)=0.5*(xyzmax(i)+xyzmin(i))
        varnames=trim(varnames)//' cut'//trim(coord(i))
     end if
  end do

  if(ndim==2)then
     !Make a lookup table to check coinciding cells
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
        if(TypeGeometry(1:9)=='spherical') then
           ! Use -pi/2 < theta' < 3/2* pi as generalized coordinate
           UseDoubleCut = .true.; nx2 = 2*nx2
        elseif(TypeGeometry=='cylindrical')then
           ! Use 0 < r' < 2*rmax as generalized coordinate
           UseDoubleCut = .true.; nx1 = 2*nx1
        end if
     end if

     if(.not.structured)then
        if(real(nx1)*real(nx2) > 1e8)then
           write(*,*)'PostIDL WARNING: very fine grid, no averaging is done!'
        else
           allocate(lookup(nx1,nx2),stat=iError)

           if(iError/=0 .or. size(lookup) < real(nx1)*real(nx2)-0.9 )then
              write(*,*)'Allocating lookup table was not successful!'
              write(*,*)'iError,size(lookup)=',iError,size(lookup)
              write(*,*)'No averaging is done!'
           else
              UseLookup=.true.
              lookup=0
              ! Cell sizes have to be stored for unstructured 2D grid
              allocate(dxdoubled(ncell))
           end if
           write(*,*)'allocate done'
        end if
     end if
  endif

  ! For unstructured grid make the xx and w arrays linear
  if(.not.structured)then
     nxyz(1)=ncell
     nxyz(2:3)=1
  end if

  ! Get components for sake of efficiency
  xmin=xyzmin(1); ymin=xyzmin(2); zmin=xyzmin(3)
  dx  =dxyz(1);   dy  =dxyz(2);   dz  =dxyz(3)
  nx  =nxyz(1);   ny  =nxyz(2);   nz  =nxyz(3)

  ! Cell aspect ratios
  dyperdx=dxyzmin(2)/dxyzmin(1); dzperdx=dxyzmin(3)/dxyzmin(1)

  ! Allocate w and xx, the arrays of variables and coordinates
  allocate(w1(nw),w(nx,ny,nz,nw),xx(nx,ny,nz,ndim),STAT=iError)
  if(iError /= 0) stop 'PostIDL.exe ERROR: could not allocate arrays'

  if(read_binary.and.nByteRealRead==4) allocate(State4_V(nw))
  if(read_binary.and.nByteRealRead==8) allocate(State8_V(nw))


  !Initialize w
  w=0.0

  !Calculate xx for structured grid
  if(structured)then
     do k=1,nz
        xcut(3)=zmin+(k-0.5)*dz
        do j=1,ny
           xcut(2)=ymin+(j-0.5)*dy
           do i=1,nx
              xcut(1)=xmin+(i-0.5)*dx
              do idim=1,ndim
                 xx(i,j,k,idim)=xcut(icutdim(idim))
              end do
           end do
        end do
     end do
  endif

  call set_strings

  ! Collect info from all files and put it into w and xx
  total=0.0
  icell=0
  countcell=0
  l=len_trim(filenamehead)
  do me=0,numprocs-1
     write(filename,'(a,i4,a)')filenamehead(1:l-2)//"_pe",me,'.idl'
     do ll=l+1,l+6
        if(filename(ll:ll)==' ')filename(ll:ll)='0'
     end do
     if(me==0)write(*,*)'reading files=',trim(filename),&
          '...',numprocs-1,'.idl'

     if(read_binary)then
        open(unit_tmp,file=filename,status='old',form='unformatted')
     else
        open(unit_tmp,file=filename,status='old')
     end if

     ! Read file
     do
        !Debug
        !write(*,*)'START READING'
        if(read_binary)then
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

        countcell=countcell+1
        dycell=dxcell*dyperdx; dzcell=dxcell*dzperdx

        if(TypeGeometry == 'cartesian')then
           XyzGen_D = Xyz_D
        else
           call set_gen_coord
        end if

        if(.not.structured)then

           if(.not.UseLookup)then
              ! In unstructured 3D grid or a very fine 2D grid
              ! no averaging is possible
              icell=icell+1
              call weighted_average(1.,0.,-1,icell)
              CYCLE
           endif

           ! Calculate indices for lookup table
           ix1=nint((XyzGen_D(idim1)-xmin1)/dx1+halfeps)
           ix2=nint((XyzGen_D(idim2)-xmin2)/dx2+halfeps)

           call unstructured_2D

           ! We are finished with unstructured
           CYCLE
        endif

        x = XyzGen_D(1); y = XyzGen_D(2); z = XyzGen_D(3)

        if(dxcell<dx+1.e-6)then
           ! Cell has the correct size or finer
           i=max(1,nint((x-xmin)/dx+0.5))
           j=max(1,nint((y-ymin)/dy+0.5))
           k=max(1,nint((z-zmin)/dz+0.5))

           if(dxcell<dx-1.e-6)then
              ! Cell is finer, calculate volume fraction
              frac=(dxcell/dx)**ndim
           else
              frac=1.0
           end if
           w(i,j,k,:)=w(i,j,k,:)+frac*w1
           total=total+frac
        else
           ! Cell is coarser than required resolution
           imin=min(nx,max(1,nint((x-0.5*dxcell-xmin)/dx+1)))
           imax=min(nx,max(1,nint((x+0.5*dxcell-xmin)/dx)))
           jmin=min(ny,max(1,nint((y-0.5*dycell-ymin)/dy+1)))
           jmax=min(ny,max(1,nint((y+0.5*dycell-ymin)/dy)))
           kmin=min(nz,max(1,nint((z-0.5*dzcell-zmin)/dz+1)))
           kmax=min(nz,max(1,nint((z+0.5*dzcell-zmin)/dz)))

           ! First order prolongation
           do iw=1,nw
              w(imin:imax,jmin:jmax,kmin:kmax,iw)= &
                   w(imin:imax,jmin:jmax,kmin:kmax,iw)+w1(iw)
           end do

           if(imax<imin.or.jmax<jmin.or.kmax<kmin)&
                write(*,*)'!!! Empty box for cell dx,x,y,z=',dxcell,x,y,z

           total=total+(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
        end if
     end do ! read file

999  continue

     close(unit_tmp)
  end do ! me

  if(countcell/=ncell)&
       write(*,*)'!!! Discrepancy: countcell=',countcell,' ncell=',ncell,' !!!'

  if(structured)then
     volume=product(real(nxyz))
     if(ndim==1 .and. abs(total/volume-4.0)<0.0001)then
        w=0.25*w
        write(*,*)'Averaged 1D structured file everywhere'
     elseif(abs(total/volume-2.0)<0.0001)then
        w=0.5*w
        write(*,*)'Averaged structured file everywhere'
     elseif(abs(total/volume-1.0)>0.0001)then
        write(*,*)'!!! Discrepancy in structured file:',&
             'filled total=',total,' volume=',volume,' !!!'
     end if
  else
     if(UseLookup)then
        volume=(xmax1-xmin1)*(xmax2-xmin2)
        ! For axysimmetric cut planes with phi being the negligible coordinate
        ! we plot both phi=cut and phi=cut+pi, so the volume is doubled
        if(UseDoubleCut) volume = 2*volume

        if(abs(total/volume-1.0)<0.0001)then
           write(*,*)'Averaged 2D unstructured file'
        else
           write(*,*)'!!! Discrepancy in averaging 2D unstructured file:',&
                'filled total=',total,' volume=',volume,' !!!'
        end if
     else
        if(ndim/=2.and.icell /= ncell) &
             write(*,*)'!!! Error: ncell=',ncell,' /= icell=',icell,' !!!'
     end if
     nx=icell
     nxyz(1)=icell
  end if

  filename=filenamehead(1:l-2)//'.out'
  write(*,*)'writing file =',trim(filename)
  if(write_binary)then
     open(unit_tmp,file=filename,status='unknown',form='unformatted')
     call save_vacfile_bin
  else
     open(unit_tmp,file=filename,status='unknown')
     call save_vacfile_ascii
  end if

  deallocate(w1, w, xx, eqpar)
  if(read_binary.and.nByteRealRead==4) deallocate(State4_V)
  if(read_binary.and.nByteRealRead==8) deallocate(State8_V)
  if(UseLookup) deallocate(lookup, dxdoubled)

  write(*,'(a)')'PostIDL finished'

contains
  !===========================================================================

  subroutine unstructured_2D

    ! Cell size
    dxyzcell(1)=dxcell; dxyzcell(2)=dycell; dxyzcell(3)=dzcell
    dx1cell=dxyzcell(idim1); dx2cell=dxyzcell(idim2)

    if(dx1cell>1.9*dx1)then
       ! Lookup indices of possible finer pairs
       ixmin1=nint((XyzGen_D(idim1)-0.25*dx1cell-xmin1)/dx1+halfeps)
       ixmax1=nint((XyzGen_D(idim1)+0.25*dx1cell-xmin1)/dx1+halfeps)
       ixmin2=nint((XyzGen_D(idim2)-0.25*dx2cell-xmin2)/dx2+halfeps)
       ixmax2=nint((XyzGen_D(idim2)+0.25*dx2cell-xmin2)/dx2+halfeps)

    endif

    jcell=lookup(ix1,ix2)

    if(jcell>0)then
       ! A cell has already been found for this projected location

       ! Check relative size of current cell with respect to the pair
       select case(nint(dxdoubled(jcell)/dxcell))
       case(1)
          ! Finer neighbor, check the four corners
          call check_corners
       case(2)
          ! Same size pair, use simple average
          call weighted_average(0.5,0.5,jcell,jcell)

          ! Negate lookup for safety check 
          lookup(ix1,ix2)=-lookup(ix1,ix2)
       case(4,8,16)
          ! Coarse pair but this is NOT the last fine neighbor yet
          ! Create new cell with weighted average
          icell=icell+1
          call weighted_average(2./3.,1./3.,jcell,icell)

          ! Increase dxdoubled by another factor of 2 to count fine neighbors
          dxdoubled(jcell)=dxdoubled(jcell)*2

          ! Negate lookup for safety check
          lookup(ix1,ix2)=-lookup(ix1,ix2)
       case(32)
          ! Coarse pair and this is the LAST fine neighbor
          ! Use weighted average and overwrite jcell
          call weighted_average(2./3.,1./3.,jcell,jcell)

          total = total - 3*dx1cell*dx2cell

          ! Negate lookup for safety check
          lookup(ix1,ix2)=-lookup(ix1,ix2)
       case default
          write(*,*)''!!! Error: Impossible dx ratio !!!'
          write(*,*)'ix1,ix2,icell,dxcell,xyz=',ix1,ix2,icell,dxcell,XyzGen_D
          write(*,*)'jcell,dxdoubled,xx=',jcell,dxdoubled(jcell),&
               xx(jcell,1,1,:)
          stop
       end select

    elseif(jcell==0)then
       ! No same size pair found yet
       ! Check the corners for finer neighbors

       call check_corners
    else
       ! Negative lookup value means an error
       write(*,*)'!!! Error: 3rd data for same projected position in ',filename
       write(*,*)'ix1,ix2,icell,jcell,dx,xyz=', &
            ix1,ix2,icell,jcell,dxcell,XyzGen_D
       stop
    end if

  end subroutine unstructured_2D

  !===========================================================================

  subroutine check_corners

    integer :: i1,i2,count

    ! Check four corners for finer neighbors
    count=0
    if(dx1cell>1.9*dx1)then
       do i1=ixmin1,ixmax1,ixmax1-ixmin1
          do i2=ixmin2,ixmax2,ixmax2-ixmin2
             jcell=lookup(i1,i2)
             if(jcell==0)then
                ! Mark lookup table for possible fine neighbor
                lookup(i1,i2)=icell+1
                CYCLE
             endif

             if(jcell<0)then
                write(*,*) '!!! Error: negative jcell when looking for finer'
                write(*,*)'ix1,ix2,icell,jcell,dx,xyz=',&
                     ix1,ix2,icell,jcell,dxcell,XyzGen_D
                stop
             endif

             count=count+1
             if(nint(dxdoubled(jcell)/dxcell)/=1)then
                write(*,*) '!!! Error: incorrect finer cell size !!!'
                write(*,*)'ix1,ix2,icell,xyz=',ix1,ix2,icell,XyzGen_D
                write(*,*)'i1,i2,jcell,xx(j)=',i1,i2,jcell,xx(jcell,1,1,:)
                write(*,*)'dxdoubled, dxcell=',dxdoubled(jcell),dxcell
                stop
             end if

             ! Average current data with finer neighbor
             call weighted_average(1./3.,2./3.,jcell,jcell)

             ! Negate lookup for safety check
             lookup(i1,i2)=-lookup(i1,i2)
          end do
       end do
    endif
    if(count<4)then
       ! Not all (if any) fine neighbors were found
       ! Store data in new cell

       icell=icell+1
       call weighted_average(1.,0.,-1,icell)

       ! Save cell size and index in lookup position(s)
       dxdoubled(icell)=2**(count+1)*dxcell
       lookup(ix1,ix2)=icell
    endif

  end subroutine check_corners

  !===========================================================================

  subroutine weighted_average(new_weight,from_weight,from_cell,to_cell)

    ! Average current cell coordinates and cell values with 
    ! an already stored twice bigger cell indexed by from_cell
    ! and put the result into cell to_cell

    real,    intent(in) :: new_weight, from_weight
    integer, intent(in) :: from_cell,to_cell
    !---------------------------------------------------------------
    if(from_cell<0)then
       w(to_cell,1,1,:)=w1
       do idim=1,ndim
          xx(to_cell,1,1,idim)=Xyz_D(icutdim(idim))
       end do
    else
       w(to_cell,1,1,:)=new_weight*w1+from_weight*w(from_cell,1,1,:)
       do idim=1,ndim
          xx(to_cell,1,1,idim)=new_weight*Xyz_D(icutdim(idim)) + &
               from_weight*xx(from_cell,1,1,idim)
       enddo
    end if
    if(to_cell/=from_cell)total=total+dx1cell*dx2cell

  end subroutine weighted_average

  !===========================================================================

  subroutine set_strings

    ! Produce fileheadout for the VAC file based on the name of the headerfile
    fileheadout=unitnames(1:494)

    ! Length of the headerline
    l=len_trim(fileheadout)

    !Replace underscores with dashes 
    do ll=1,l
       if(fileheadout(ll:ll)=='_')fileheadout(ll:ll)='-'
    enddo
    ! Add _xxx13, _xxx23 or _xxx33 to fileheadout based on ndim
    ! The _xxx comes from filenamehead (e.g. y=0_var_... --> _var 
    ! where var is 2 or 3 character-long)
    ll = 7; if(filenamehead(7:7) == '_') ll = 6
    write(fileheadout,'(a,i1,i1)') fileheadout(1:l)//filenamehead(4:ll),ndim,3

    ! Produce coordinate names 
    !         ('x y z', 'x y', 'x z', 'y z' or 'r theta', 'r phi' ...)

    coordnames=coord(icutdim(1))
    ! Fix coordinate name for non-cartesian cut along phi=90 deg:
    if(index(filenamehead,'x=0')>0) coordnames='y'
    do idim = 2, ndim
       coordnames=trim(coordnames)//' '//trim(coord(icutdim(idim)))
    end do
    varnames=trim(coordnames)//' '//trim(varnames)

  end subroutine set_strings

  !==========================================================================

  subroutine save_vacfile_ascii

    write(unit_tmp,"(a)")trim(fileheadout)
    if(structured)then
       write(unit_tmp,"(i7,1pe13.5,3i3)")it,t,ndim,neqpar+nspecialpar,nw
       write(unit_tmp,"(3i4)") (nxyz(icutdim(idim)),idim=1,ndim)
    else
       write(unit_tmp,"(i7,1pe13.5,3i3)")it,t,-ndim,neqpar+nspecialpar,nw
       write(unit_tmp,"(i8,i2,i2)") nxyz(1:ndim)
    endif
    write(unit_tmp,"(100(1pe13.5))")eqpar,specialpar(1:nspecialpar)
    write(unit_tmp,"(a)")trim(varnames)

    do k= 1,nz
       do j= 1,ny
          do i= 1,nx
             write(unit_tmp,"(100(1pe18.10))")xx(i,j,k,:),w(i,j,k,:)
          enddo
       enddo
    enddo

  end subroutine save_vacfile_ascii

  !==========================================================================

  subroutine save_vacfile_bin

    integer :: iw, i, j, k, idim, n

    real, allocatable :: Buffer_I(:)
    !----------------------------------------------------------------------
    write(unit_tmp)fileheadout

    if(structured)then
       write(unit_tmp)it,t,ndim,neqpar+nspecialpar,nw
       write(unit_tmp)(nxyz(icutdim(idim)),idim=1,ndim)
    else
       write(unit_tmp)it,t,-ndim,neqpar+nspecialpar,nw
       write(unit_tmp)nxyz(1:ndim)
    endif
    write(unit_tmp)eqpar,specialpar(1:nspecialpar)
    write(unit_tmp)varnames

    ! Use a 1D buffer, because some versions of the NAG compiler are
    ! not able to write out a 4D array into an unformatted file correctly :-(
    allocate(Buffer_I(nx*ny*nz*ndim))
    n = 0
    do idim=1,ndim; do k=1,nz; do j=1,ny; do i=1,nx; n = n + 1;
       Buffer_I(n) = xx(i,j,k,idim)
    end do; end do; end do; end do
    write(unit_tmp) Buffer_I
    deallocate(Buffer_I)
    allocate(Buffer_I(nx*ny*nz))
    do iw=1,nw
       n = 0
       do k=1,nz; do j=1,ny; do i=1,nx; n = n + 1;
          Buffer_I(n) = w(i,j,k,iw)
       end do; end do; end do
       write(unit_tmp) Buffer_I
    end do
    deallocate(Buffer_I)

  end subroutine save_vacfile_bin

  !===========================================================================
  subroutine set_gen_coord

    ! Calculate the generalized coordinates mostly for lookup

    real, parameter:: cPi= 3.1415926535897932384626433832795
    real, parameter:: cTwoPi = 2*cPi, cHalfPi = cPi/2, cRadToDeg=180/cPi

    real:: rCyl ! distance from axis Z

    ! Toroidal variables
    real:: PoloidalAngle, r, z, StretchCoef, dAngle, Residual, WallRadius

    !---------------------------------------------------------------------
    rCyl = sqrt(Xyz_D(1)**2 + Xyz_D(2)**2)

    ! Calculate phi
    if (rCyl == 0.0) then
       XyzGen_D(2) = 0.0
    else
       XyzGen_D(2) = modulo(atan2(Xyz_D(2), Xyz_D(1)), cTwoPi)
    end if

    select case(TypeGeometry)
    case('cylindrical')
       XyzGen_D(1) = rCyl
       XyzGen_D(3) = Xyz_D(3)

       if(ndim==2)then
          ! Set the 'X-Y' coordinates for plotting a 2D cut
          select case(idim0)
          case(1)
             ! This is R=const slice, use longitude [deg] vs height
             Xyz_D(2)=XyzGen_D(2)*cRadToDeg
          case(2)
             ! This is x=0 or y=0 plane, use signed radius vs Z
             Xyz_D(1) = sign(1.0, Xyz_D(1)+Xyz_D(2))*rCyl
             ! The generalized coordinate will run from 0 to 2r?
             XyzGen_D(1) = Xyz_D(1) + XyzMax(1)
          end select
       end if

    case('spherical','spherical_lnr')
       XyzGen_D(1) = sqrt(rCyl**2 + Xyz_D(3)**2)

       if(ndim==2)then
          ! Set the 'X-Y' coordinates for plotting a 2D cut
          select case(idim0)
          case(1)
             ! This is R=const slice, use longitude vs latitude in degs.
             Xyz_D(2:3)=XyzGen_D(2:3)*cRadToDeg
          case(2)
             ! This is x=0 or y=0 plane, use axial radius vs Z
             Xyz_D(1) = sign(1.00,Xyz_D(1)+Xyz_D(2))*rCyl
          case(3)
             ! Stretch X and Y with rSph/rCyl since instead of simply
             ! projecting the points down to the X-Y plane
             Xyz_D(1:2)=Xyz_D(1:2)*XyzGen_D(1)/rCyl
          end select
       end if
       if(ndim==2 .and. idim0==2) then
          ! Use 360 degrees for theta to distinguish left and right halves
          XyzGen_D(3) = cHalfPi - atan2(Xyz_D(3),Xyz_D(1))
          if (XyzGen_D(3) < -cHalfPi) XyzGen_D(3) = XyzGen_D(3) + cTwoPi
       else
          ! Colatitude
          XyzGen_D(3) = cHalfPi - acos(Xyz_D(3)/XyzGen_D(1))
       end if
       if(TypeGeometry=='spherical_lnr') XyzGen_D(1) = log(XyzGen_D(1))

    case('axial_torus')

       !This is x=0 or y=0 plane
       if(ndim==2 .and. idim0==2) Xyz_D(1) = rCyl

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
  end subroutine set_gen_coord

end program PostIDL
