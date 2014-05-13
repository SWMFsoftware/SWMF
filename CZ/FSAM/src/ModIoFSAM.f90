module ModIoFSAM
  implicit none

  private

  public :: readrst_mpi

  public :: writedata_mpi

  public :: writegrid

contains
  
  !===================================================================================

  subroutine writegrid
    use ModPar,      ONLY: inmax, jnmax, knmax, myid
    use ModOutfile,  ONLY: gridfile
    use ModGrid
    implicit none
    
    integer i,j,k
    real :: vect1(inmax), vect2(jnmax)
    !---------------------------------------------------------------------------------

    if(myid .eq. 0) then
       open(unit=13, file=trim(gridfile)//'.dat', form='unformatted',access='stream')
       write(13) inmax,jnmax,knmax
    endif

    ! write x1a
    call combine_nproc1(x1a,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write x1b
    call combine_nproc1(x1b,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write x2a
    call combine_nproc2(x2a,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax)
    ! write x2b
    call combine_nproc2(x2b,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax-1)
    ! write x3a
    if(myid .eq. 0) write(13) (x3a(k),k=1,knmax)
    ! write x3b
    if(myid .eq. 0) write(13) (x3b(k),k=1,knmax-1)
    ! write dx1a
    call combine_nproc1(dx1a,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write dx1b
    call combine_nproc1(dx1b,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write dx2a
    call combine_nproc2(dx2a,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax-1)
    ! write dx2b
    call combine_nproc2(dx2b,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax)
    ! write dx3a
    if(myid .eq. 0) write(13) (dx3a(k),k=1,knmax-1)
    ! write dx3b
    if(myid .eq. 0) write(13) (dx3b(k),k=1,knmax)
    ! write g2a
    call combine_nproc1(g2a,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write g2b
    call combine_nproc1(g2b,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write g31a
    call combine_nproc1(g31a,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write g31b
    call combine_nproc1(g31b,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write g32a
    call combine_nproc2(g32a,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax)
    ! write g32b
    call combine_nproc2(g32b,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax-1)
    ! write dg2bd1
    call combine_nproc1(dg2bd1,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write dg2ad1
    call combine_nproc1(dg2ad1,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write dg31bd1
    call combine_nproc1(dg31bd1,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write dg31ad1
    call combine_nproc1(dg31ad1,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write dg32bd2
    call combine_nproc2(dg32bd2,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax)
    ! write dg32ad2
    call combine_nproc2(dg32ad2,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax-1)
    ! write dvl1a
    call combine_nproc1(dvl1a,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write dvl1b
    call combine_nproc1(dvl1b,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write dvl2a
    call combine_nproc2(dvl2a,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax-1)
    ! write dvl2b
    call combine_nproc2(dvl2b,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax)
    ! write dvl3a
    if(myid .eq. 0) write(13) (dvl3a(k),k=1,knmax-1)
    ! write dvl3b
    if(myid .eq. 0) write(13) (dvl3b(k),k=1,knmax)
    ! write dx1ai
    call combine_nproc1(dx1ai,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write dx1bi
    call combine_nproc1(dx1bi,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write dx2ai
    call combine_nproc2(dx2ai,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax-1)
    ! write dx2bi
    call combine_nproc2(dx2bi,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax)
    ! write dx3ai
    if(myid .eq. 0) write(13) (dx3ai(k),k=1,knmax-1)
    ! write dx3bi
    if(myid .eq. 0) write(13) (dvl3bi(k),k=1,knmax)
    ! write g2ai
    call combine_nproc1(g2ai,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write g2bi
    call combine_nproc1(g2bi,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write g31ai
    call combine_nproc1(g31ai,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write g31bi
    call combine_nproc1(g31bi,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write g32ai
    call combine_nproc2(g32ai,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax)
    ! write g32bi
    call combine_nproc2(g32bi,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax-1)
    ! write dvl1ai
    call combine_nproc1(dvl1ai,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax-1)
    ! write dvl1bi
    call combine_nproc1(dvl1bi,vect1)
    if(myid .eq. 0) write(13) (vect1(i),i=1,inmax)
    ! write dvl2ai
    call combine_nproc2(dvl2ai,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax-1)
    ! write dvl2bi
    call combine_nproc2(dvl2bi,vect2)
    if(myid .eq. 0) write(13) (vect2(j),j=1,jnmax)
    ! write dvl3ai
    if(myid .eq. 0) write(13) (dvl3ai(k),k=1,knmax-1)
    ! write dvl3bi
    if(myid .eq. 0) write(13) (dvl3bi(k),k=1,knmax)

    if(myid .eq. 0) close(13)
    
  end subroutine writegrid
  
  !===================================================================================
  
  subroutine combine_nproc1(subvect,vect)
    use ModPar,  ONLY: in, inmax, myid, nproc1, myid1, myid2
    use ModMpi
    use ModFSAM, ONLY: iComm
    implicit none
    
    real :: subvect(in),vect(inmax), sendbuf1(in),recvbuf1(in)
    integer :: cputo,cpufrom,icpu, istatus(MPI_STATUS_SIZE), ierr
    !--------------------------------------------------------------------------------
    
    if(myid==0) then       
       vect(1:in) = subvect
       ! if necessary receive data
       if(nproc1 > 1) then
          do icpu=1,nproc1-1
             cpufrom = icpu
             call MPI_Recv(recvbuf1, in, MPI_DOUBLE_PRECISION, cpufrom, &
                  icpu, iComm, istatus, ierr)
             vect(icpu*(in-5)+1:icpu*(in-5)+in) = recvbuf1(1:in)
          enddo
       endif
    else
       if(myid2==0) then
          cputo = 0
          sendbuf1 = subvect
          call MPI_Send(sendbuf1, in, MPI_DOUBLE_PRECISION, cputo, &
               myid1, iComm, ierr)
       endif
    endif
    
    call MPI_BARRIER(iComm, ierr)
    call MPI_BCAST(vect, inmax, MPI_DOUBLE_PRECISION, 0, iComm, ierr)
    
  end subroutine combine_nproc1
  
  !====================================================================================
  
  subroutine combine_nproc2(subvect,vect)
    use ModPar,   ONLY: myid, myid1, myid2, nproc1, nproc2, jn, jnmax
    use ModMpi
    use ModFSAM,  ONLY: iComm
    implicit none
    
    real :: subvect(jn), vect(jnmax), sendbuf2(jn), recvbuf2(jn)
    integer :: cputo, cpufrom, jcpu, istatus(MPI_STATUS_SIZE), ierr
    !----------------------------------------------------------------------------------
    
    if(myid == 0) then       
       vect(1:jn) = subvect
       ! if necessary receive data
       if(nproc2 > 1) then
          do jcpu=1,nproc2-1
             cpufrom=jcpu*nproc1
             call MPI_Recv(recvbuf2,jn,MPI_DOUBLE_PRECISION, cpufrom, &
                  jcpu, iComm, istatus, ierr)
             vect(1+jcpu*(jn-5):jn+jcpu*(jn-5)) = recvbuf2
          enddo
       endif
    else
       if(myid1 == 0) then
          cputo=0
          sendbuf2 = subvect
          call MPI_Send(sendbuf2,jn,MPI_DOUBLE_PRECISION, cputo, myid2, &
               iComm, ierr)
       endif
    endif
    
    call MPI_BARRIER(iComm, ierr)
    call MPI_BCAST(vect, jnmax, MPI_DOUBLE_PRECISION, 0, iComm, ierr)
    
  end subroutine combine_nproc2
  
  !===================================================================================
  
  subroutine readrst_mpi
    use ModSundry,  ONLY: time
    use ModPar,     ONLY: myid, myid1, myid2, in, jn, kn, inmax, jnmax, knmax
    use ModOutfile, ONLY: b1file, b2file, b3file, v1file, v2file, v3file, sfile, pfile
    use ModField,   ONLY: v1, v2, v3, b1, b2, b3, s, p
    use ModMpi
    use ModFSAM,    ONLY: iComm
    implicit none
    
    character(len=20) :: varfile
    integer(kind=MPI_OFFSET_KIND) :: disp
    real, allocatable :: buf(:,:,:)
    integer :: filetype, fileinfo, fhandl, mysize, ierr, starts(1:3)
    integer, parameter, dimension(1:3) :: sizes = (/ inmax-1, jnmax-1, knmax-1 /), &
         subsizes = (/ in-1, jn-1, kn-1 /)
    !--------------------------------------------------------------------------------
    
    starts = (/ myid1*(in-5), myid2*(jn-5), 0 /)
    
    call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
         MPI_DOUBLE_PRECISION, filetype, ierr)
    call MPI_Type_commit(filetype, ierr)
    call MPI_Info_create(fileinfo, ierr)
    
    varfile = trim(v1file)//'rst.dat'
    allocate(buf(1:in-1,1:jn-1,1:kn-1))
    call readvar_mpi(varfile, buf)
    v1(1:in-1, 1:jn-1, 1:kn-1) = buf
    
    varfile = trim(v2file)//'rst.dat'
    call readvar_mpi(varfile, buf)
    v2(1:in-1, 1:jn-1, 1:kn-1) = buf
    
    varfile = trim(v3file)//'rst.dat'
    call readvar_mpi(varfile, buf)
    v3(1:in-1, 1:jn-1, 1:kn-1) = buf
    
    varfile = trim(b1file)//'rst.dat'
    call readvar_mpi(varfile, buf)
    b1(1:in-1, 1:jn-1, 1:kn-1) = buf
    
    varfile = trim(b2file)//'rst.dat'
    call readvar_mpi(varfile, buf)
    b2(1:in-1, 1:jn-1, 1:kn-1) = buf
    
    varfile = trim(b3file)//'rst.dat'
    call readvar_mpi(varfile, buf)
    b3(1:in-1, 1:jn-1, 1:kn-1) = buf
    
    varfile = trim(sfile)//'rst.dat'
    call readvar_mpi(varfile, buf)
    s(1:in-1, 1:jn-1, 1:kn-1) = buf
    
    varfile = trim(pfile)//'rst.dat'
    call readvar_mpi(varfile, buf)
    p(1:in-1, 1:jn-1, 1:kn-1) = buf
    
    deallocate(buf)
    
  contains
    !---------------------------------------------------------------------------------
    subroutine readvar_mpi(varfile, buf)
      implicit none
      
      character(len=20), intent(in)  :: varfile
      real,              intent(out) :: buf(1:in-1,1:jn-1,1:kn-1)
      
      integer :: mpi_status(MPI_STATUS_SIZE)
      integer :: inmax_in, jnmax_in, knmax_in, i, j, k
      real    :: time_in
      !------------------------------------------------------------------------------

      fhandl=13
      call MPI_Barrier(iComm,ierr)
      call MPI_File_open(iComm, varfile, MPI_MODE_RDONLY, fileinfo, &
           fhandl, ierr)
      if(myid .eq. 0) then
         call MPI_File_read(fhandl, time_in, 1, MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_read(fhandl, inmax_in, 1, MPI_INTEGER, mpi_status, ierr)
         call MPI_File_read(fhandl, jnmax_in, 1, MPI_INTEGER, mpi_status, ierr)
         call MPI_File_read(fhandl, knmax_in, 1, MPI_INTEGER, mpi_status, ierr)
      endif
      call MPI_BCAST(time_in, 1, MPI_DOUBLE_PRECISION, 0, iComm, ierr)
      time = time_in
      call MPI_BCAST(inmax_in, 1, MPI_INTEGER, 0, iComm, ierr)
      call MPI_BCAST(jnmax_in, 1, MPI_INTEGER, 0, iComm, ierr)
      call MPI_BCAST(knmax_in, 1, MPI_INTEGER, 0, iComm, ierr)
      if(inmax_in .ne. inmax .or. jnmax_in .ne. jnmax .or. knmax_in .ne. knmax) then
         write(6,*) 'inmax,jnmax,knmax in ', varfile,' :', inmax_in, jnmax_in, &
              knmax_in, ' do not match inmax, jnmax, knmax in ModPar.f90'
         call MPI_ABORT(iComm, 1,ierr)
      endif
      !inquire(IOLENGTH = disp) time_in, inmax, jnmax, knmax
      disp = sizeof(time_in) + 3*sizeof(inmax)
      call MPI_FILE_SET_VIEW(fhandl, disp, MPI_DOUBLE_PRECISION, filetype, &
           'native', fileinfo, ierr)
      mysize = subsizes(1)*subsizes(2)*subsizes(3)
      call MPI_File_read_all(fhandl, buf, mysize, MPI_DOUBLE_PRECISION, &
           mpi_status, ierr)
      call MPI_Barrier(iComm, ierr)
      call MPI_File_close(fhandl, ierr)
      
    end subroutine readvar_mpi
    
  end subroutine readrst_mpi
  
  !===================================================================================

  subroutine writedata_mpi(DoWriteRestart)
    use ModPar,      ONLY: inmax, jnmax, knmax, nproc1, nproc2, myid1, myid2, myid
    use ModGrid,     ONLY: in, jn, kn, is, js
    use ModOutfile,  ONLY: v1file, v2file, v3file, b1file, b2file, b3file, sfile, &
         pfile, ifile
    use ModSundry,   ONLY: time
    use ModField,    ONLY: v1, v2, v3, b1, b2, b3, s, p
    use ModBack,     ONLY: fact
    use ModMpi
    use ModFSAM,     ONLY: iComm
    implicit none
    
    logical, intent(in), optional :: DoWriteRestart
    integer, parameter, dimension(1:3) :: sizes = (/ inmax-1, jnmax-1, knmax-1 /)
    integer, dimension(1:3) :: subsizes, starts
    integer :: filetype, fileinfo, fhandl, ierr
    character(len=20) :: idout
    character(len=30):: filename
    !---------------------------------------------------------------------------------
    
    subsizes = (/    in-5,    jn-5,    kn-1 /)
    starts   = (/ is-1+myid1*(in-5), js-1+myid2*(jn-5), 0 /)
    
    if(myid1 .eq. 0) then
       subsizes(1) = subsizes(1) + 2
       starts(1)   = 0
    endif
    if(myid1 .eq. nproc1-1) subsizes(1) = subsizes(1) + 2
    if(myid2 .eq. 0) then
       subsizes(2) = subsizes(2) + 2
       starts(2)   = 0
    endif
    if(myid2 .eq. nproc2-1) subsizes(2) = subsizes(2) + 2
    
    call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
         MPI_DOUBLE_PRECISION, filetype, ierr)
    call MPI_Type_commit(filetype, ierr)
    call MPI_Info_create(fileinfo, ierr)
    
    fhandl = 13
    write(idout,'(I4.4)') ifile
    if(present(DoWriteRestart).and.DoWriteRestart) write(idout,'(a3)') 'rst'
    call write_var_file(trim(v1file)//trim(idout)//'.dat', v1)
    call write_var_file(trim(v2file)//trim(idout)//'.dat', v2)
    call write_var_file(trim(v3file)//trim(idout)//'.dat', v3)
    call write_var_file(trim(b1file)//trim(idout)//'.dat', b1)
    call write_var_file(trim(b2file)//trim(idout)//'.dat', b2)
    call write_var_file(trim(b3file)//trim(idout)//'.dat', b3)
    call write_var_file(trim( sfile)//trim(idout)//'.dat',  s)
    call write_var_file(trim( pfile)//trim(idout)//'.dat',  p, DoPressure=.true.)
    
  contains
    !----------------------------------------------------------------------------------
    subroutine write_var_file(filename, var, DoPressure)
      implicit none
      
      character(len=*), intent(in) :: filename
      real,             intent(in) :: var(1:in,1:jn,1:kn)
      logical, optional,intent(in) :: DoPressure
      
      integer(kind=MPI_OFFSET_KIND) :: disp
      integer :: mpi_status(MPI_STATUS_SIZE), i, j, k, counter, mysize
      real, allocatable :: buf(:)
      !--------------------------------------------------------------------------------
      
      call MPI_Barrier(iComm,ierr)
      call MPI_File_open(iComm, filename, &
           IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), fileinfo, fhandl, ierr)
      if(myid .eq. 0) then
         call MPI_File_write(fhandl, time, 1, MPI_DOUBLE_PRECISION, mpi_status, ierr)
         call MPI_File_write(fhandl, inmax, 1, MPI_INTEGER, mpi_status, ierr)
         call MPI_File_write(fhandl, jnmax, 1, MPI_INTEGER, mpi_status, ierr)
         call MPI_File_write(fhandl, knmax, 1, MPI_INTEGER, mpi_status, ierr)
      endif
      !inquire(IOLENGTH = disp) time, inmax, jnmax, knmax
      disp = sizeof(time) + 3*sizeof(inmax)
      call MPI_FILE_SET_VIEW(fhandl, disp, MPI_DOUBLE_PRECISION, filetype, 'native', &
           fileinfo, ierr)
      mysize = subsizes(1)*subsizes(2)*subsizes(3)
      allocate(buf(in*jn*kn))
      do k=starts(3)+1,starts(3)+subsizes(3)
         do j=starts(2)+1,starts(2)+subsizes(2)
            do i=starts(1)+1,starts(1)+subsizes(1)
               counter = 1+i-starts(1)-1 +(j-starts(2)-1)*subsizes(1) + &
                    (k-starts(3)-1)*subsizes(1)*subsizes(2)
               buf(counter) = var(i-myid1*(in-5),j-myid2*(jn-5),k)
               if(present(DoPressure).and.DoPressure) buf(counter) = &
                    buf(counter)/fact(i)
            enddo
         enddo
      enddo
      call MPI_File_write_all(fhandl, buf, mysize, MPI_DOUBLE_PRECISION, &
           mpi_status, ierr)
      deallocate(buf)
      call MPI_Barrier(iComm,ierr)
      call MPI_File_close(fhandl, ierr)
      
    end subroutine write_var_file
    
  end subroutine writedata_mpi

end module ModIoFSAM
