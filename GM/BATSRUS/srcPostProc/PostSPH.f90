!This code is a copyright protected software (c) 2002- University of Michigan
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                                                                |
!/////////////////////////////////////////////////////////////////
program PostSPH
  implicit none
  
  !\
  ! Some array definitions.
  !/
  integer, parameter :: linesize = 500
  integer, parameter :: maxtheta = 360     ! 90 degrees * 4
  integer, parameter :: maxphi   = 1440    !360 degrees * 4
  integer, parameter :: maxerror = 50000   !~maxtheta*maxphi*.01
  integer :: overcounts(maxtheta,maxphi)
  integer :: out_of_range_i(maxerror),out_of_range_j(maxerror)

  character(LEN=linesize), allocatable :: lines_out(:,:)
  integer :: iError

  logical :: is_missing(maxtheta,maxphi),is_overcounted(maxtheta,maxphi)
  integer, dimension(7) :: date_time_array

  !\
  ! Some IO stuff
  !/
  integer, parameter :: unit_in=30, unit_out=31
  integer :: l
  character (len=80) :: filename_in, filename_out, filename_original
  character (LEN=4) :: TimeH4
  character (LEN=2) :: TimeM2,TimeS2

  !\
  ! Some other variables
  !/
  integer :: i,j,ii
  integer :: nnodes, nnodesmax, nlines, ndouble, nmissing, nout_of_range
  integer :: ntheta, nphi
  integer :: imin,imax,jmin,jmax
  real :: dtheta, dphi
  real :: x_out,y_out,z_out,theta_out,phi_out,thetatilt,phitilt
  logical :: doDebug

  !\
  ! Read in variables
  !/
  integer :: numprocs, it
  character (LEN=linesize) :: varnames, hemisphere, line
  real :: t, rplot


  !\
  ! initialize some variables
  !/
  doDebug=.false.
  is_missing = .true.
  is_overcounted = .false.
  overcounts = 0
  nout_of_range = 0
  out_of_range_i = 0
  out_of_range_j = 0

  nmissing  = 0
  ndouble   = 0
  nnodes    = 0
  nnodesmax = 0
  nlines    = 0

  dtheta    = 0
  dphi      = 0
  ntheta    = 0
  nphi      = 0

  filename_in  = 'inputfile.DAT'
  filename_out = 'outputfileN.DAT'

  ! Begin by reading the headerfile

  write(*,'(a)') '-------------------------------------------------------------------'
  write(*,*) ' '
  write(*,'(a)') 'PostSPH: Post-processing package for converting'
  write(*,'(a)') '         concatenated spherical plot files into'
  write(*,'(a)') '         tecplot format.'
  write(*,'(a)') '    -Written by K.C. Hansen, 2001,2003'

  ! Read information from STDIN
  read(*,'(a)') filename_original
  read(*,*) numprocs

  read(*,*) it
  read(*,*) t
  read(*,'(a)') varnames
  read(*,*) dtheta, dphi
  read(*,*) date_time_array(:) 
  read(*,*) thetatilt, phitilt
  read(*,*) rplot
  read(*,'(a)') hemisphere

  write(TimeH4,'(i4.4)') &
       int(                            T/3600.)
  write(TimeM2,'(i2.2)') &
       int((T-(3600.*int(T/3600.)))/60.)
  write(TimeS2,'(i2.2)') &
       int( T-(  60.*int(T/  60.)))


  ! Compute the number of nodes that we expect

  ntheta = 1 + 90.0/dtheta
  nphi   = 360.0/dphi

  allocate(lines_out(ntheta, nphi), STAT=iError)
  if (iError/=0) stop 'PostSPH could not allocate lines_out array'

  nnodesmax = ntheta*nphi
  nmissing = nnodesmax


  ! write output to the screen

  write(*,*)                        ' '
  write(*,*)                        'Writing plot for ', hemisphere(1:len_trim(hemisphere))
  write(*,*)                        '   - R=',rplot
  write(*,*)                        ' '
  write(*,*)                        '    n_step = ',it
  write(*,'(" ",a,1pe12.5)')        '  sim_time = ',t
  write(*,*)                        ' date_time = ',date_time_array
  write(*,*)                        ' '
  write(*,'(" ",a,i5,a,1pe12.5,a)') '    nTheta = ',nTheta, '     dTheta = ',dtheta,' deg'
  write(*,'(" ",a,i5,a,1pe12.5,a)') '      nPhi = ',nPhi  , '       dPhi = ',dphi,  ' deg'
  write(*,*)                        ' '
  write(*,*)                        'theta_tilt = ',thetatilt,' deg'
  write(*,*)                        '  phi_tilt = ',phitilt,  ' deg'

  if (doDebug) then
     write(*,*) ' '
     write(*,*) filename_original(1:len_trim(filename_original))
     write(*,*) 'varnames =',varnames(1:len_trim(varnames))
  end if


  ! Start the calculation.  Keep track of missing nodes.

  imin = ntheta
  imax = 0
  jmin = nphi
  jmax = 0

  if (doDebug) write(*,'(a,i9)') 'Expect -> nnodes:', nnodesmax


  ! Now read in the data file line by line and store the line in an array
  ! This will allow the array to be written out in the correct order for tecplot
  ! and also will allow statistics to be built

  open(unit_in,file=filename_in,status='old')

  do 

     read(unit_in,'(i7,i7,a)',END=999)i,j,line

     nlines = nlines + 1

     imin = min(imin,i)
     imax = max(imax,i)
     jmin = min(jmin,j)
     jmax = max(jmax,j)

     if (i>maxtheta .or. j>maxphi) then
        nout_of_range = nout_of_range+1
        out_of_range_i(nout_of_range) = i
        out_of_range_j(nout_of_range) = j
     else

        if (is_missing(i,j)) then
           is_missing(i,j) = .false.
           lines_out(i,j) = line
           nnodes = nnodes + 1
           nmissing = nmissing-1
        else
           is_overcounted(i,j) = .true.
           ndouble = ndouble + 1
           overcounts(i,j) = overcounts(i,j) + 1
        end if

     end if

  end do

999 continue

  ! Now do some error reporting
  if (doDebug) then
     write(*,*) 'Total number of lines read:           ', nlines
     write(*,*) 'Total number of nodes expected:       ', nnodesmax
     write(*,*) 'Total number of unique nodes read:    ', nnodes
     write(*,*) 'Total number of redundant nodes read: ', ndouble
     write(*,*) 'Total number of nodes missing:        ', nmissing
     write(*,*) 'Total number of nodes out of range:   ', nout_of_range
     write(*,*) ''
     write(*,*) 'imin,imax: ', imin,imax
     write(*,*) 'jmin,jmax: ', jmin,jmax
   
     write(*,*) ''
     write(*,*) '****************List of double counted nodes follows*****************'
  end if

  ii = 0  
  ndouble = 0
  do i=1,ntheta; do j=1,nphi
     if(is_overcounted(i,j))  then
        ii = ii+1 
        ndouble = ndouble + overcounts(i,j)
        if (doDebug) write(*,'(a,i7,i7,i7,i7,i7)') 'list#,overcount#,i,j,ncounts: ', &
               ii,ndouble,i,j,overcounts(i,j)+1
     end if
  end do; end do

  if (doDebug) write(*,*) ''
  if (doDebug) write(*,*) '****************List of missing nodes follows************************'

  ii = 0  
  do i=1,ntheta; do j=1,nphi
     if(is_missing(i,j))  then
        ii = ii+1 
        if (doDebug) write(*,'(a,i7,i7,i7)') 'list#,i,j: ', &
               ii,i,j
     end if
  end do; end do

  if (doDebug) write(*,*) ''
  if (doDebug) write(*,*) '****************List of out of range nodes follows*******************'

  do i=1,nout_of_range
     if (doDebug) write(*,'(a,i7,i7,i7)') 'list#,i,j: ', &
              i,out_of_range_i(i),out_of_range_j(i)
  end do
  
  write(*,*) ''

  ! Write out the tecplot file header.
  open(unit_out,file=filename_out,status='unknown')

  write(unit_out,'(a,i4.4,"-",5(i2.2,"-"),i3.3,a,2f6.2,a,1pe11.3,a)') & 
     'TITLE="BATSRUS: Spherical Shell, ',date_time_array(:),', BTiltDeg=',thetatilt,phitilt, &
     ', R= ',rplot, ', '//hemisphere(1:len_trim(hemisphere))//'"'
  write(unit_out,'(a)') varnames(1:len_trim(varnames))
  write(unit_out,'(a,i7.7,a,i8,a,i8,a)') 'ZONE T="Sph'//hemisphere(1:1)//' N=',it, &
     ' T='//TimeH4//':'//TimeM2//':'//TimeS2//'"'//', I=',ntheta,', J=',nphi+1,', K=1, F=POINT'

  ! Write out all of the data to the file
  do j=1,nphi; do i=1,ntheta
     write(unit_out,'(a)') lines_out(i,j)(1:len_trim(lines_out(i,j)))
  end do; end do



  ! write the first phi point again at the end to close the grid in tecplot
  do i=1,ntheta
     read(lines_out(i,1),'(5(E14.6),a)') x_out,y_out,z_out,theta_out,phi_out,line
     phi_out = phi_out+360.0
     write(unit_out,'(5(E14.6),a)') x_out,y_out,z_out,theta_out,phi_out, &
                 line(1:len_trim(line))
  end do





  
end program PostSPH



