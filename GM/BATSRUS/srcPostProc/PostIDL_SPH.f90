!^CFG COPYRIGHT UM
!^CFG FILE COVARIANT
!=============================================================================
program PostIDL

  ! Read a .h file from STDIN, then read all the files based on the .h file
  ! and put structured x and w together and save into a VAC file

  implicit none

  ! Global variables

  integer, parameter :: unit_tmp=99
  logical, parameter :: write_binary=.false.
  real, parameter :: halfeps=0.6

  integer :: nxyz(3), icell, countcell, ncell, nw, neqpar, numprocs, it
  real :: t
  real, dimension(:,:,:,:), allocatable :: xx, w
  real, dimension(:), allocatable :: w1, eqpar, dxdoubled

  ! Coordinates, sizes, indices
  real, dimension(3) :: XyzCart_D,XyzGen_D, xyzmin, xyzmax,&
       dxyz, dxyzmin, dxyzcell
  real ::    x, y, z, xmin, xmax, ymin, ymax, zmin, zmax
  real ::    dx, dy, dz, dyperdx, dzperdx, dxcell, dycell, dzcell
  real ::    frac, total, volume
  integer :: i, j, k, imin, imax, jmin, jmax, kmin, kmax, nx, ny, nz, iw

  integer :: idim, icutdim(3), ndim, nspecialpar
  real    :: xcut(3), specialpar(3)
  character (LEN=5), dimension(3), parameter :: &
       coord_xyz=(/'x    ','y    ','z    '/), &
       coord_sph=(/'r    ','phi  ','theta'/)

  logical :: fileopen, structured, read_binary=.false.,&
       IsSpherical=.true., UseLnr=.false.
  character (len=79) :: fileheadout, varnames, filename, filenamehead, &
       coordnames,TypeGeometry
  character (len=500) :: unitnames
  integer :: l, ll, me

  ! Variables for the 2D lookup table
  integer :: ix1,ix2,ixmin1,ixmax1,ixmin2,ixmax2,nx1,nx2,idim1,idim2,idim0
  integer:: jcell=0
  real    :: xmin1, xmax1, xmin2, xmax2, dx1, dx2, dx1cell, dx2cell,RCyl
  integer, dimension(:,:), allocatable :: lookup

  ! Variables for checking binary compatibility
  integer :: nByteReal, nByteRealRead
  real,parameter::cPi= 3.1415926535897932384626433832795
  !---------------------------------------------------------------------------

  write(*,'(a)')'PostIDL (G.Toth 2000-2002) starting'

  ! Read information from STDIN
  read(*,'(a)')filenamehead
  read(*,*)numprocs
  write(*,*)filenamehead(1:len_trim(filenamehead)),', numprocs=',numprocs


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
  write(*,*)'varnames =',varnames(1:len_trim(varnames))
  read(*,'(a)')unitnames
  write(*,*)'unitnames=',unitnames(1:len_trim(unitnames))

  read_binary = .false.
  read(*,'(l8)',err=1,end=1)read_binary
1 continue
  write(*,*)'binary   =',read_binary

  if(read_binary)then
     inquire(iolength=nByteReal)1.0
     nByteRealRead = -1
     read(*,'(i8)',err=2,end=2)nByteRealRead
2    continue
     if(nByteRealRead==nByteReal)then
        write(*,*)'nByteReal=',nByteReal
     else if(nByteRealRead==-1)then
        write(*,*)'!!! Warning: PostIDL was compiled with ',&
            nByteReal,' byte reals but nByteReal is not given in file !!!'
     else
        write(*,*)'!!! Error: PostIDL was compiled with ',&
            nByteReal,' byte reals but file contains nByteReal=',nByteRealRead
        stop '!!! Change PRECISION in Makefile.${OS} and make PIDL !!!'
     end if
  end if
  
  read(*,'(a)',err=3,end=3)TypeGeometry
  IsSpherical=index(TypeGeometry,'spher')>0
  UseLnr=IsSpherical.and.index(TypeGeometry,'lnr')>0
3 continue
  ! Unstructured grid has dx=-1.
  structured = dxyz(1) > 0.00

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
        varnames=varnames(1:len_trim(varnames))//' cut'//coord_sph(i)(1:len_trim(coord_sph(i)))
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
   !  if(idim0==2)nx2=2*nx2                  !^CFG UNCOMMENT IF COVARIANT
     allocate(lookup(nx1,nx2))
     lookup=0
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
  allocate(w1(nw),w(nx,ny,nz,nw),xx(nx,ny,nz,ndim))

  ! Cell sizes have to be stored for unstructured 2D grid
  if(.not.structured.and.ndim==2)allocate(dxdoubled(ncell))

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
     if(me==0)write(*,*)'reading files=',filename(1:len_trim(filename)),&
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
           read(unit_tmp,ERR=999,END=999)dxcell,XyzCart_D,w1
        else
           read(unit_tmp,*,ERR=999,END=999)dxcell,XyzCart_D,w1
        end if
        
        countcell=countcell+1
        dycell=dxcell*dyperdx; dzcell=dxcell*dzperdx
        
        XyzGen_D=XyzCart_D          !^CFG IF NOT COVARIANT

! call cart_to_spher(XyzCart_D,XyzGen_D)               !^CFG UNCOMMENT IF COVARIANT

        if(ndim==2)then
           ! Calculate indices for lookup table
           ix1=nint((XyzGen_D(idim1)-xmin1)/dx1+halfeps)
           ix2=nint((XyzGen_D(idim2)-xmin2)/dx2+halfeps)
   !        select case(idim0)                                       !^CFG UNCOMMENT IF COVARIANT
   !        case(1)                                                  !^CFG UNCOMMENT IF COVARIANT
   !            XyzCart_D(2)=mod(XyzGen_D(2)*12./cPi+12.0,24.0)      !^CFG UNCOMMENT IF COVARIANT
   !            XyzCart_D(3)=90.0-XyzGen_D(3)*180./cPi               !^CFG UNCOMMENT IF COVARIANT
   !        case(2)                                                  !^CFG UNCOMMENT IF COVARIANT
   !            XyzCart_D(1)=&                                       !^CFG UNCOMMENT IF COVARIANT
   !            sign(1.00,XyzCart_D(1)+XyzCart_D(2))*RCyl            !^CFG UNCOMMENT IF COVARIANT
   !        case(3)                                                  !^CFG UNCOMMENT IF COVARIANT
   !            if(UseLnr) then                                      !^CFG UNCOMMENT IF COVARIANT
   !               XyzCart_D(1)=XyzCart_D(1)*exp(XyzGen_D(1))/Rcyl   !^CFG UNCOMMENT IF COVARIANT
   !               XyzCart_D(2)=XyzCart_D(2)*exp(XyzGen_D(1))/Rcyl   !^CFG UNCOMMENT IF COVARIANT
   !            else                                                 !^CFG UNCOMMENT IF COVARIANT
   !               XyzCart_D(1)=XyzCart_D(1)*XyzGen_D(1)/Rcyl        !^CFG UNCOMMENT IF COVARIANT
   !               XyzCart_D(2)=XyzCart_D(2)*XyzGen_D(1)/Rcyl        !^CFG UNCOMMENT IF COVARIANT
   !            end if                                               !^CFG UNCOMMENT IF COVARIANT
   !        end select                                               !^CFG UNCOMMENT IF COVARIANT
        end if

        if(.not.structured)then

           if(ndim/=2)then
              ! In unstructured 3D grid no averaging is possible
              icell=icell+1
              call weighted_average(1.,0.,-1,icell)
              CYCLE
           endif

           call unstructured_2D

           ! We are finished with unstructured
           CYCLE
        endif
        x=XyzGen_D(1); y=XyzGen_D(2); z=XyzGen_D(3)

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

     999 continue

     close(unit_tmp)
  end do ! me

  if(countcell/=ncell)&
       write(*,*)'!!! Discrepancy: countcell=',countcell,' ncell=',ncell,' !!!'

  if(structured)then
     volume=product(real(nxyz))
     if(abs(total/volume-2.0)<0.0001)then
        w=0.5*w
        write(*,*)'Averaged 2D structured file everywhere'
     elseif(abs(total/volume-1.0)>0.0001)then
        write(*,*)'!!! Discrepancy in averaging 2D structured file:',&
             'filled total=',total,' volume=',volume,' !!!'
     end if
  else
     if(ndim/=2.and.icell /= ncell) &
           write(*,*)'!!! Error: ncell=',ncell,' /= icell=',icell,' !!!'
     nx=icell
     nxyz(1)=icell
  end if

  filename=filenamehead(1:l-2)//'.out'
  write(*,*)'writing file =',filename(1:len_trim(filename))
  if(write_binary)then
     open(unit_tmp,file=filename,status='unknown',form='unformatted')
     call save_vacfile_bin
  else
     open(unit_tmp,file=filename,status='unknown')
     call save_vacfile_ascii
  end if

  deallocate(w1,w,xx,eqpar)
  if(ndim==2)deallocate(lookup)
  if(ndim==2.and. .not.structured)deallocate(dxdoubled)

  write(*,'(a)')'PostIDL finished'

contains

  subroutine unstructured_2D

    ! Cell size
    dxyzcell(1)=dxcell; dxyzcell(2)=dycell; dxyzcell(3)=dzcell
    dx1cell=dxyzcell(idim1); dx2cell=dxyzcell(idim2)

    if(dx1cell>dx1)then
       ! Lookup indices of possible finer pairs
       ixmin1=nint((XyzGen_D(idim1)-0.25*dx1cell-xmin1)/dx1+halfeps)
       ixmax1=nint((XyzGen_D(idim1)+0.25*dx1cell-xmin1)/dx1+halfeps)
       ixmin2=nint((XyzGen_D(idim2)-0.25*dx2cell-xmin2)/dx2+halfeps)
       ixmax2=nint((XyzGen_D(idim2)+0.25*dx2cell-xmin2)/dx2+halfeps)

    endif

    jcell=lookup(ix1,ix2)
    if(jcell>0)then
       ! A cell has already been found for this projected location

       ! Check relative size of curret cell with respect to the pair
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

          ! Increase dxdoubled by another factor of 2 to count fine neigbors
          dxdoubled(jcell)=dxdoubled(jcell)*2
 
          ! Negate lookup for safety check
          lookup(ix1,ix2)=-lookup(ix1,ix2)
       case(32)
          ! Coarse pair and this is the LAST fine neigbor
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
       write(*,*)'!!! Error: 3rd data for same projected position'
       write(*,*)'ix1,ix2,icell,jcell,dx,xyz=',ix1,ix2,icell,jcell,dxcell,XyzGen_D
       stop
    end if

  end subroutine unstructured_2D

  subroutine check_corners

    integer :: i1,i2,count

    ! Check four corners for finer neighbors
    count=0
    if(dx1cell>dx1)then
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

  subroutine weighted_average(new_weight,from_weight,from_cell,to_cell)

    ! Average current cell coordinates and cell values with 
    ! an already stored twice bigger cell indexed by from_cell
    ! and put the result into cell to_cell

    real,    intent(in) :: new_weight, from_weight
    integer, intent(in) :: from_cell,to_cell
       if(from_cell<0)then
          w(to_cell,1,1,:)=w1
          do idim=1,ndim
             xx(to_cell,1,1,idim)=XyzCart_D(icutdim(idim))
          end do
       else
          w(to_cell,1,1,:)=new_weight*w1+from_weight*w(from_cell,1,1,:)
          do idim=1,ndim
             xx(to_cell,1,1,idim)=new_weight*XyzCart_D(icutdim(idim)) + &
                  from_weight*xx(from_cell,1,1,idim)
          enddo
       end if
    if(to_cell/=from_cell)total=total+dx1cell*dx2cell

  end subroutine weighted_average

  subroutine set_strings

    ! Produce fileheadout for the VAC file based on the name of the headerfile
    fileheadout=unitnames(1:73)

    ! Length of the headerline
    l=len_trim(fileheadout)

    !Replace underscores with dashes 
    do ll=1,l
       if(fileheadout(ll:ll)=='_')fileheadout(ll:ll)='-'
    enddo
    ! Add _mhd13, _mhd23 or _mhd33 to fileheadout based on ndim
    write(fileheadout,'(a,i1,i1)') fileheadout(1:l)//filenamehead(10:13),ndim,3

    ! Produce coordinate names ('x y z ', 'x y ', 'x z ', 'y z ' or 'theta','phi')
    if(index(filenamehead,'x=0')>0) then
       coordnames=coord_xyz(2)
    elseif(index(filenamehead,'x=y')>0)then
       coordnames=coord_sph(1)
    else
       coordnames=coord_xyz(icutdim(1))
    end if
    do idim=2,ndim
       coordnames=coordnames(1:len_trim(coordnames)+1)//&
            coord_xyz(icutdim(idim))(1:len_trim(coord_xyz(icutdim(idim))))
    end do
    varnames=coordnames(1:len_trim(coordnames)+1)//varnames

  end subroutine set_strings

  !==========================================================================
  subroutine save_vacfile_ascii

    write(unit_tmp,"(a)")fileheadout
    if(structured)then
       write(unit_tmp,"(i7,1pe13.5,3i3)")it,t,ndim,neqpar+nspecialpar,nw
       write(unit_tmp,"(3i4)") (nxyz(icutdim(idim)),idim=1,ndim)
    else
       write(unit_tmp,"(i7,1pe13.5,3i3)")it,t,-ndim,neqpar+nspecialpar,nw
       write(unit_tmp,"(i8,i2,i2)") nxyz(1:ndim)
    endif
    write(unit_tmp,"(100(1pe13.5))")eqpar,specialpar(1:nspecialpar)
    write(unit_tmp,"(a)")varnames

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

    integer :: iw, i, j, k, idim

    !Debug
    !write(*,*)'save_vacfile_bin start nx,ny,nz,nw w(98,1,58,:)=',&
    !     nx,ny,nz,nw,w(98,1,58,:)

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
    write(unit_tmp)((((xx(i,j,k,idim),i=1,nx),j=1,ny),k=1,nz),idim=1,ndim)
    do iw=1,nw
       write(unit_tmp)(((w(i,j,k,iw),i=1,nx),j=1,ny),k=1,nz)
    end do

  end subroutine save_vacfile_bin
  subroutine cart_to_spher(Cart_D,Spher_D)
    implicit none
    real,dimension(3),intent(in)::Cart_D
    real,dimension(3),intent(out)::Spher_D
    RCyl=sqrt(Cart_D(1)**2+Cart_D(2)**2)
    Spher_D(1)= sqrt(Rcyl**2+Cart_D(3)**2)
    if (Cart_D(1) == 0.00 .and. Cart_D(2) == 0.00) then
       Spher_D(2) = 0.00
    else
       Spher_D(2) = atan2(Cart_D(2), Cart_D(1))
    end if
    if (Spher_D(2) < 0.00) Spher_D(2) = Spher_D(2) + cPi*2
    if(ndim==2.and.idim0==2) then
       Spher_D(3) = atan2(Cart_D(3),sign(1.00,Cart_D(1)+Cart_D(2))*RCyl)
       if (Spher_D(3) < 0.00) Spher_D(3) = Spher_D(3) + cPi*2
    else
       Spher_D(3) = acos(Cart_D(3)/Spher_D(1))
    end if
    if(UseLnr)Spher_D(1)=log(Spher_D(1))
  end subroutine cart_to_spher
end program PostIDL

