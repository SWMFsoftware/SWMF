PROGRAM Extract_arrays
  USE Rcm_variables
  USE Rcm_IO
  IMPLICIT NONE
  CHARACTER (LEN=10) :: time_string
  CHARACTER (LEN=80) :: filename
!     
  INTEGER :: rec_beg,rec_end, n_rec, ii,jj,kk, i0
  INTEGER, DIMENSION (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc) :: Iarray_2d
  REAL, DIMENSION (1-n_gc:isize+n_gc, 1-n_gc:jsize+n_gc) :: Rarray_2d
!
!
!
1 WRITE (*,'(/A)') 'PROGRAM TO PRINT OUT RCM ARRAYS'
  WRITE (*,'(A)') 'PRINTING OUT:  POTENTIAL, BIRKLAND CURRENTS'
  WRITE (*,'(A)') 'enter file name for output'
  READ (*,*) filename
  filename = TRIM (filename)
  OPEN (UNIT=50, FILE=trim(NameRcmDir)//filename, STATUS='REPLACE')
  write(50,'(a)') 'TITLE="RCM"'
  write(50,'(a)') 'VARIABLES='
  write(50,*) '"I", "J"'
  write(50,*) '"colat", "aloct"'
  write(50,*) '"X [GSM]", "Y [GSM]"'
  write(50,*) '"bndloc"'
  write(50,*) '"Potential [kV]"'
  write(50,*) '"Birkland Currents"'
  do kk=1,kcsize
     write(50,'(a,i3.3,a)') ' "Eta ',kk,'"'
  end do
!
  WRITE (6,'(A)',ADVANCE = 'NO') 'ENTER RECORD # range: '
  READ (5,*) rec_beg,rec_end
!
  do n_rec = rec_beg,rec_end
!
     CALL Read_grid ()
     CALL Read_plasma ()
!
     CALL Read_array ('rcmbirk',   n_rec, label, ARRAY_2D = birk)
     CALL Read_array ('rcmv'   ,   n_rec, label, ARRAY_2D = v)
     CALL Read_array ('rcmbndloc', n_rec, label, ARRAY_1D = bndloc)
     CALL Read_array ('rcmxmin',   n_rec, label, ARRAY_2D = xmin)
     CALL Read_array ('rcmymin',   n_rec, label, ARRAY_2D = ymin)
!     CALL Read_array ('rcmbi',     n_rec, label, ARRAY_1D = bi)
!     CALL Read_array ('rcmbj',     n_rec, label, ARRAY_1D = bj)
!     CALL Read_array ('rcmitrack', n_rec, label, ARRAY_1D = itrack)
!     CALL Read_array ('rcmmpoint', n_rec, label, ARRAY_1D = mpoint)
!     CALL Read_array ('rcmnpoint', n_rec, label, ARRAY_1D = npoint)
     CALL Read_array ('rcmeeta',   n_rec, label, ARRAY_3D = eeta)
!
     time_string = TRIM(Get_time_char_string(label))
     WRITE (6,'(A,I6.6,A,A)') ' RECORD ',n_rec,' TIME IS ',time_string
!
!
     write(50,'(a,a,a,i6.6,a,i4,a,i4,a)') &
          'ZONE T="T=',time_string,' REC=',n_rec,'", I=',isize, &
          ', J=',jsize+1,', K=1, F=BLOCK'
!
     do jj = 1,jsize+1; do ii = 1,isize
        Iarray_2d(ii,jj) = ii
     end do; end do
     call writeInteger
!
     do jj = 1,jsize+1; do ii = 1,isize
        Iarray_2d(ii,jj) = jj
     end do; end do
     call writeInteger
!
     Rarray_2d = colat
     call writeReal
!
     Rarray_2d = aloct
!!;     Rarray_2d(:,jsize) = 2.*Rarray_2d(:,jsize-1) - Rarray_2d(:,jsize-2)
     call writeReal
!
     Rarray_2d = xmin
     call writeReal
!
     Rarray_2d = ymin
     call writeReal
!
     do jj=1,jsize+1
        Rarray_2d(:,jj) = bndloc(jj)
     end do
     call writeReal
!
     Rarray_2d = v/1000.
     call writeReal
!
     Rarray_2d = birk
     call writeReal
!
     do kk=1,kcsize
        Rarray_2d = eeta(:,:,kk)
        call writeReal
     end do
!
  end do
  CLOSE (50)
!
  STOP
  
Contains
  
  subroutine writeInteger
    do jj = 1,jsize+1
       do ii = 1,isize,16
          write(50,'(16i4)') (Iarray_2d(i0,jj), i0=ii,min(ii+15,isize))
       end do
    end do
    write(50,*) ' '
  end subroutine writeInteger
  
  subroutine writeReal
    do jj = 1,jsize+1
       do ii = 1,isize,5
          write(50,'(5G14.6)') (Rarray_2d(i0,jj), i0=ii,min(ii+4,isize))
       end do
    end do
    write(50,*) ' '
  end subroutine writeReal
  
END PROGRAM Extract_arrays
