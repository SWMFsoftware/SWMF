MODULE Rcm_io
    USE Rcm_variables, ONLY : iprec, rprec, label_def, lun, NameRcmDir
    IMPLICIT NONE
!
!    use ModMpiOrig
!
    INTERFACE Read_array
       MODULE PROCEDURE Read_real_1d_array, Read_real_2d_array, Read_real_3d_array,&
                        Read_intg_1d_array, Read_intg_2d_array, Read_intg_3d_array
    END INTERFACE
!
    INTERFACE Write_array
       MODULE PROCEDURE Write_real_1d_array, Write_real_2d_array, Write_real_3d_array,&
                        Write_intg_1d_array, Write_intg_2d_array, Write_intg_3d_array
    END INTERFACE
!
    INTERFACE Outp
       MODULE PROCEDURE Outp_real, Outp_integer, Outp_logical
    END INTERFACE
!
!
CONTAINS
!
!
    SUBROUTINE Read_real_1d_array (filename, rec_num, label, &
                                   array_1d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    REAL (rprec),      INTENT (OUT) :: array_1d (:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: length, istat
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
!
!    write(*,*) '!!! Read_real_1d_array:',filename, rec_num, asci_format

    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_1d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_1d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_1d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_REAL_1D_ARRAY'
       WRITE (*,'(T2,A)') 'FILE NAME: '//filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF
!
    OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN       
       WRITE (*,*) 'ERROR OPENING FILE: ',trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90::read_real_1d_array')
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, &
             ERR = 1, FMT = form_string)        label, array_1d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1) label, array_1d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file ', filename, rec_num, Istat
       write (*,*) label
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90::Read_real_1d_array')
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_real_1d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_real_2d_array (filename, rec_num, label, &
                                   array_2d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    REAL (rprec),      INTENT (OUT) :: array_2d (:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    INTEGER (iprec)    :: length, istat
    CHARACTER (LEN=80) :: form_string
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
!
!    write(*,*) '!!! Read_real_2d_array:',filename, rec_num, asci_format

    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_2d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_2d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_2d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_REAL_2D_ARRAY'
       WRITE (*,'(T2,A)') 'FILE NAME: '//filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF
!
    OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90::Read_real_2d_array')
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1, FMT = form_string) label, array_2d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1 ) label, array_2d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file is: ', filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90::Read_real_2d_array')
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_real_2d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_real_3d_array (filename, rec_num, label, array_3d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    REAL (rprec),      INTENT (OUT) :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: k, rec_num_mod
    LOGICAL            :: flag, asci_format
!
!
!    write(*,*) '!!! Read_real_1d_array:',filename, rec_num

    DO k = 1, SIZE (array_3d, DIM = 3)
       rec_num_mod = (rec_num - 1) * SIZE(array_3d, DIM = 3) + k
       CALL Read_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                        error_flag, asci)
    END DO
!
    RETURN
    END SUBROUTINE Read_real_3d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_intg_1d_array (filename, rec_num, label, &
                                   array_1d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    INTEGER (iprec),   INTENT (OUT) :: array_1d (:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: length, istat
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
!
!    write(*,*) '!!! Read_intg_1d_array:',filename, rec_num, asci_format

    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_1d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_1d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,I10))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_1d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_INTG_1D_ARRAY'
       WRITE (*,'(T2,A)') 'FILE NAME: '//filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF
!
    OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90::Read_intg_1d_array')
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1, FMT = form_string) label, array_1d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1) label, array_1d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file is: ', filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90::Read_intg_1d_array')
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_intg_1d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_intg_2d_array (filename, rec_num, label, &
                                   array_2d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    INTEGER (iprec),   INTENT (OUT) :: array_2d (:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: length, istat
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF

!    write(*,*) '!!! Read_intg_2d_array:',filename, rec_num, asci_format
!
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_2d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_2d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,I10))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_2d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_INTG_2D_ARRAY'
       WRITE (*,'(T2,A)') 'FILE NAME: '//filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF
!
    OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90::Read_intg_2d_array')
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1, FMT = form_string) label, array_2d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1                   ) label, array_2d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file is: ', filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90::Read_intg_2d_array')
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_intg_2d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_intg_3d_array (filename, rec_num, label, &
                                   array_3d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    INTEGER (iprec),   INTENT (OUT) :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    INTEGER (iprec)    :: k, rec_num_mod
!
!    write(*,*) '!!! Read_intg_3d_array:',filename, rec_num
!
    DO k = 1, SIZE (array_3d, DIM = 3)
       rec_num_mod = (rec_num - 1) * SIZE(array_3d, DIM = 3) + k
       CALL Read_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                        error_flag, asci)
    END DO
!
    RETURN
    END SUBROUTINE Read_intg_3d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_real_1d_array (filename, rec_num, label, array_1d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    REAL (rprec),      INTENT (IN)      :: array_1d (:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    CHARACTER (LEN=11):: form_type_char
    CHARACTER (LEN=80) :: form_string
    LOGICAL           :: flag, asci_format
!
!    write(*,*) '!!! Write_real_1d_array:',filename, rec_num, present(setup), present(asci)

    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       write(*,'(T2,A)')  'ERROR in Write_real_1D array'
       WRITE (*,'(T2,A)') 'FILE NAME: '//filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90:UNIT ALREADY OPEN')
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF

!    write(*,*) '!!! Write_real_1d_array:',filename, rec_num, status_char, asci_format


    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_1d
!
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_1d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxxxx(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_1d)
       form_type_char = 'FORMATTED  '
    END IF

    OPEN (UNIT=LUN, FILE = trim(NameRcmDir)//filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = form_type_char,&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90::Write_real_1d_array')
    END IF

    IF (asci_format) THEN
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat, FMT = form_string) label, array_1d
    ELSE
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat) label, array_1d
    END IF
    IF (istat /= 0) then
       write(*,*) 'ERROR READING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_real_1d_array')
    END IF
    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) THEN
       write(*,*) 'ERROR CLOSING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_real_1d_array')
    END IF
    RETURN
    END SUBROUTINE Write_real_1d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_real_2d_array (filename, rec_num, label, array_2d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    REAL (rprec),      INTENT (IN)      :: array_2d (:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    CHARACTER (LEN=11):: form_type_char
    CHARACTER (LEN=80) :: form_string
    LOGICAL           :: flag, asci_format
!

!    write(*,*) '!!! Write_real_2d_array:',filename, rec_num, present(setup), present(asci)

    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       write(*,'(T2,A)')  'ERROR in Write_real_2D array'
       write(*,*)'FILENAME=',filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90:UNIT ALREADY OPEN')
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF

!    write(*,*) '!!! Write_real_2d_array:',filename, rec_num, status_char, asci_format

    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_2d
!
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_2d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxxxx(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_2d)
       form_type_char = 'FORMATTED  '
    END IF

    OPEN (UNIT=LUN, FILE = trim(NameRcmDir)//filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = form_type_char,&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF

    IF (asci_format) THEN
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat, FMT = form_string) label, array_2d
    ELSE
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat                   ) label, array_2d
    END IF
    IF (istat /= 0) THEN
       write(*,*) 'ERROR READING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_real_2d_array')
    END IF

    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) THEN
       write(*,*) 'ERROR CLOSING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_real_2d_array')
    END IF

    RETURN
    END SUBROUTINE Write_real_2d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_real_3d_array (filename, rec_num, label, array_3d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    REAL (rprec),      INTENT (IN)      :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: k, rec_num_mod
    CHARACTER (LEN=7) :: status_char
!
!    write(*,*) '!!! Write_real_3d_array:',filename, rec_num, present(setup), present(asci)

    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
!
    DO k = 1, SIZE (array_3d, DIM = 3)
       rec_num_mod = (rec_num - 1) * SIZE (array_3d, DIM = 3) + k
       IF (k == 1) THEN
          CALL Write_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                            SETUP = setup, ASCI = asci)
       ELSE
          CALL Write_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                            ASCI = asci)
       END IF
    END DO
!
    RETURN
    END SUBROUTINE Write_real_3d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_intg_1d_array (filename, rec_num, label, array_1d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    INTEGER (iprec),   INTENT (IN)      :: array_1d (:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    CHARACTER (LEN=11):: form_type_char
    CHARACTER (LEN=80) :: form_string
    LOGICAL           :: flag, asci_format
!
!    write(*,*) '!!! Write_intg_1d_array:',filename, rec_num, present(setup), present(asci)

    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       write(*,*)'ERROR in Write_intg_1D array, filename=',filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90:UNIT ALREADY OPEN')
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF

!    write(*,*) '!!! Write_intg_1d_array:',filename, rec_num, status_char, asci_format

    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_1d
!
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_1d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxxxx(TR2,I10))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_1d)
       form_type_char = 'FORMATTED  '
    END IF

    OPEN (UNIT=LUN, FILE = trim(NameRcmDir)//filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = form_type_char,&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE in write_intg_1d_array: ',&
            trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF

    IF (asci_format) THEN
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat, FMT = form_string) label, array_1d
    ELSE
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat                   ) label, array_1d
    END IF
    IF (istat /= 0) THEN
       write(*,*) 'ERROR WRITING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_intg_1d_array')
    END IF

    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) THEN
       write(*,*) 'ERROR CLOSING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_intg_1d_array')
    END IF

    RETURN
    END SUBROUTINE Write_intg_1d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_intg_2d_array (filename, rec_num, label, array_2d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    INTEGER (iprec),   INTENT (IN)      :: array_2d (:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    CHARACTER (LEN=11):: form_type_char
    CHARACTER (LEN=80) :: form_string
    LOGICAL           :: flag, asci_format
!
!    write(*,*) '!!! Write_intg_2d_array:',filename, rec_num, present(setup), present(asci)

    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       write(*,*)'ERROR in Write_intg_2D array, filename=',filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90:UNIT ALREADY OPEN')
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF

!    write(*,*) '!!! Write_intg_2d_array:',filename, rec_num, status_char, asci_format
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_2d
!
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_2d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxxxx(TR2,I10))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_2d)
       form_type_char = 'FORMATTED  '
    END IF

    OPEN (UNIT=LUN, FILE = trim(NameRcmDir)//filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = form_type_char,&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF

    IF (asci_format) THEN
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat, FMT = form_string) label, array_2d
    ELSE
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat                   ) label, array_2d
    END IF
    IF (istat /= 0) THEN
       write(*,*) 'ERROR WRITING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_intg_2d_array')
    END IF

    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) THEN
       write(*,*) 'ERROR CLOSING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_intg_2d_array')
    END IF

    RETURN
    END SUBROUTINE Write_intg_2d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_intg_3d_array (filename, rec_num, label, array_3d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    INTEGER (iprec),   INTENT (IN)      :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: k, rec_num_mod
    CHARACTER (LEN=7) :: status_char
!
!    write(*,*) '!!! Write_intg_3d_array:',filename, rec_num, present(setup), present(asci)

    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
!
    DO k = 1, SIZE (array_3d, DIM = 3)
       rec_num_mod = (rec_num - 1) * SIZE (array_3d, DIM = 3) + k
       IF (k == 1) THEN
          CALL Write_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                            SETUP = setup, ASCI = asci)
       ELSE
          CALL Write_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                            ASCI = asci)
       END IF
    END DO
!
    RETURN
    END SUBROUTINE Write_intg_3d_array
!
!
    SUBROUTINE Outp_real (r, isize, jsize, ibeg, iend, iinc, jbeg, jend, jinc, &
                          xscale, ilabel, title, ntp, ncol)
    USE Rcm_variables, ONLY : n_gc
    IMPLICIT NONE
    INTEGER (iprec), INTENT(IN) :: ibeg, iend, iinc, jbeg, jend, jinc, &
                                   ntp, ncol, ilabel(20), isize, jsize
    CHARACTER (LEN = * ), INTENT (IN) :: title
    REAL (rprec), INTENT (IN) :: r (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc), xscale
!                                                                       
!------------------------------------------------------------
!
! author: r.w. spiro        last modified: 06-25-85                     
!                                                                       
! DESCRIPTION OF INPUT PARAMETERS:                          
!      r(isize,jsize)=array to be output                                
!      this subroutine can output selected elements of array     
!      ibeg= initial i value to be output                               
!      iend= final i value to be output                                 
!      iinc= i value increment for output                               
!      (note: (iend-ibeg)/iinc should be an integer)                    
!      jbeg,jend, and jinc are defined similarly                        
!      scale=scale factor                                               
!         if(scale.eq.0.) scale is calculated to give best 
!            display      
!         all elements or array are divided by scale before
!            being output
!      ilabel= vector of length 20 that gives the label                 
!      title=character string that identifies the array 
!             being output    
!      ntp= output unit number                                          
!      ncol= number of columns for output device (80 or 132)            
!                                                                       
!------------------------------------------------------------
!                                                                       
    INTEGER (iprec), PARAMETER :: jcol = 16
    REAL (rprec) ::  y (jcol), scale_tmp, sum0, sum1, sum2, test, ave, sd
    INTEGER (iprec) :: mxline, jjcol, isclmx, i, j, itest, &
                       ipower, jstart, istart, ifinal, iutt, lcount, &
                       jfinal, jj, ii, jcount
!
!   initialization
!
    scale_tmp = xscale
    isclmx = - 12
    sum0 = 0.0_rprec
    sum1 = 0.0_rprec
    sum2 = 0.0_rprec
    mxline = 57
    jjcol = (ncol - 4) / 8
    IF (ibeg > isize+n_gc .OR. iend > isize+n_gc .OR. &
        jbeg > jsize+n_gc .OR. jend > jsize+n_gc .OR. &
        ibeg < 1-n_gc .OR. iend < 1-n_gc .OR. jbeg < 1-n_gc .OR. jend < 1-n_gc) THEN
        call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90:INDICES WRONG IN OUTP')
    END IF
!                                                                       
!   if scale=0. then compute auto scale factor
!                                                                       
    IF (ABS(scale_tmp) < TINY(1.0_rprec)) THEN
!
!      Determine maximum # of places to left of decimal pt
       DO i = ibeg, iend, iinc
       DO j = jbeg, jend, jinc
          IF (ABS(r(i,j)) > TINY(1.0_rprec)) THEN
             test = LOG10 (ABS (r (i, j) ) )
             IF (test >= 0.0_rprec) then
                test = test + .000001
             ELSE
                test = test - .999999
             END IF
             itest = test
          ELSE
             itest = 0
          END IF
          IF (r (i, j) < 0.0_rprec) itest = itest + 1
          isclmx = MAX (isclmx, itest)
       END DO
       END DO
!
!      determine scale factor such that max # of places
!      to left of decimal point is 3:
!
       ipower = isclmx - 2
       scale_tmp = 10.0_rprec**ipower
    ELSE
       ipower = NINT (LOG10 (scale_tmp) )
    END IF
!                                                                       
!   Compute mean and std dev. of outputted array elements:
!
    DO i = 1, isize
    DO j = 1, jsize
       sum0 = sum0 + 1.0_rprec
       sum1 = sum1 + r (i,j)
!       sum2 = sum2 + r (i,j) **2
    END DO
    END DO
    ave = sum1 / sum0
!    sd  = SQRT ((sum2-sum0*ave**2)/(sum0-1.0_rprec))
    sd = sum2
    sd = 0.0_rprec
!                                                                       
!   Output data:
!                                                                       
    jstart = jbeg
    istart = ibeg
    ifinal = iend
!
!   Below ilabel(6) is ITIME variable, ilabel(1) is IRDW  variable,
!                              
    WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6), ilabel(3),ilabel(4), &
                     ilabel(5)
    iutt = ilabel (2)
!
    WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
    lcount = 0
!                                                                       
!   Begin output loop:
!
    main_loop: DO
         jfinal = (jjcol - 1) * jinc + jstart 
         IF (jfinal > jend) jfinal = jend 
         IF (jstart > jend) EXIT main_loop
         IF (mxline-lcount <= 5) then 
!                                 start a new page                      
!
            WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                            ,ilabel(3),ilabel(4),ilabel(5)
            iutt = ilabel (2) 
            WRITE (ntp, 810) TRIM (title), ipower, iutt, ave, sd 
            lcount = 0 
         END IF
!
         WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc) 
         lcount = lcount + 2
         DO ii = istart, ifinal, iinc 
            lcount = lcount + 1 
            IF (lcount > mxline) then 
!                                    start a new page                   
!
               WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                               ,ilabel(3),ilabel(4),ilabel(5)
               iutt = ilabel (2) 
               WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd 
               WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc) 
               lcount = 2 
            END IF 
            jcount = 0 
            DO jj = jstart, jfinal, jinc 
               jcount = jcount + 1 
               y (jcount) = r (ii, jj) / scale_tmp
            END DO 
            WRITE (ntp, 840) ii, (y (jj), jj = 1, jcount) 
         END DO
         jstart = jfinal + jinc 
         istart = ibeg 
    END DO main_loop
!
    RETURN
800 FORMAT (//,TR3,'# ', I6.6, TR2, 'ut=', I6.6, TR2,     &
            'ITIME=',I5.5,'(', I2.2,':',I2.2,':',  I2.2,')' )
810 FORMAT (/,T4, A, '/(1.E', I2, ')', TR4, 'ut=', I6.6, TR2, &
            'ave=', ES11.4, TR2, 'sd=', ES11.4)
830 FORMAT ( / ,T3, 16(TR5,I3))
840 FORMAT (T2, I3, 16(F8.3))
    END SUBROUTINE Outp_real
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Outp_integer (r, ibeg, iend, iinc, jbeg, jend, jinc, &
                             xscale, ilabel, title, ntp, ncol)
    IMPLICIT NONE
    INTEGER (iprec), INTENT(IN) :: ibeg, iend, iinc, jbeg, jend, jinc, &
                                   ntp, ncol, ilabel(20)
    CHARACTER (LEN = * ), INTENT (IN) :: title
    REAL (rprec),INTENT (IN) :: xscale
    INTEGER (iprec), INTENT (IN) :: r (:,:)
!                                                                       
!------------------------------------------------------------
!
! author: r.w. spiro        last modified: 06-25-85                     
!                                                                       
! DESCRIPTION OF INPUT PARAMETERS:                          
!      r(isize,jsize)=array to be output                                
!      this subroutine can output selected elements of array     
!      ibeg= initial i value to be output                               
!      iend= final i value to be output                                 
!      iinc= i value increment for output                               
!      (note: (iend-ibeg)/iinc should be an integer)                    
!      jbeg,jend, and jinc are defined similarly                        
!      scale=scale factor                                               
!         if(scale.eq.0.) scale is calculated to give best 
!            display      
!         all elements or array are divided by scale before
!            being output
!      ilabel= vector of length 20 that gives the label                 
!      title=character string that identifies the array 
!             being output    
!      ntp= output unit number                                          
!      ncol= number of columns for output device (80 or 132)            
!                                                                       
!------------------------------------------------------------
!                                                                       
    INTEGER (iprec), PARAMETER :: jcol = 16
    REAL (rprec)::  y (jcol), scale_tmp, sum0, sum1, sum2, test, ave, sd
    INTEGER (iprec) :: mxline, jjcol, isclmx, i, j, itest, ipower, &
                       jstart, istart, ifinal, iutt, lcount, jfinal, &
                       jj, ii, jcount, isize, jsize
!
!   initialization
!
    isize = SIZE (r, DIM = 1)
    jsize = SIZE (r, DIM = 2)
    scale_tmp = xscale
    isclmx = - 12
    sum0 = 0.0_rprec
    sum1 = 0.0_rprec
    sum2 = 0.0_rprec
    mxline = 57
    jjcol = (ncol - 4) / 8
    IF (ibeg > isize .OR. iend > isize .OR. &
        jbeg > jsize .OR. jend > jsize .OR. &
        ibeg < 1 .OR. iend < 1 .OR. jbeg < 1 .OR. jend < 1) THEN
        call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90:INDICES WRONG IN OUTP')
    END IF
!
!   if scale=0. then compute auto scale factor
!                                                                       
    IF (ABS(scale_tmp) < TINY(1.0_rprec)) then
!
!      Determine maximum # of places to left of decimal pt
       DO i = ibeg, iend, iinc
       DO j = jbeg, jend, jinc
          IF (ABS(r(i,j)) > TINY(1.0_rprec)) THEN
             test = LOG10 (ABS (REAL(r (i, j)) ) )
             IF (test >= 0.0_rprec) THEN
                test = test + .000001
             ELSE
                test = test - .999999
             END IF
             itest = test
          ELSE
             itest = 0
          END IF
          IF (r (i, j) < 0.0_rprec) itest = itest + 1
          isclmx = MAX (isclmx, itest)
       END DO
       END DO
!
!      determine scale factor such that max # of places
!      to left of decimal point is 3:
!
       ipower = isclmx - 2
       scale_tmp = 10.0_rprec**ipower
    ELSE
       ipower = NINT (LOG10 (scale_tmp) )
    END IF
!                                                                       
!   Compute mean and std dev. of outputted array elements:
!
    DO i = 1, isize
    DO j = 1, jsize
       sum0 = sum0 + 1.0_rprec
       sum1 = sum1 + REAL (r (i,j))
!       sum2 = sum2 + REAL (r (i,j))**2
    END DO
    END DO
    ave = SUM1 / sum0
!    sd  = SQRT ((sum2-sum0*ave**2)/(sum0-1.0_rprec))
    sd = sum2
    sd = 0.0_rprec
!                                                                       
!   Output data:
!                                                                       
    jstart = jbeg
    istart = ibeg
    ifinal = iend
!
!                   Below ilabel(6) is ITIME variable,
!                         ilabel(1) is IRDW  variable,
!                              
    WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6), &
                     ilabel(3),ilabel(4),ilabel(5)
    iutt = ilabel (2)
!
    WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
    lcount = 0
!                                                                       
!   Begin output loop:
!
    main_loop: DO
       jfinal = (jjcol - 1) * jinc + jstart
       IF (jfinal > jend) jfinal = jend
       IF (jstart > jend) EXIT main_loop
       IF (mxline-lcount <= 5) then
!                               start a new page
!
          WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                          ,ilabel(3),ilabel(4),ilabel(5)
          iutt = ilabel (2)
          WRITE (ntp, 810) TRIM (title), ipower, iutt, ave, sd
          lcount = 0
       END IF
!
       WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc)
       lcount = lcount + 2
       DO ii = istart, ifinal, iinc
          lcount = lcount + 1
          IF (lcount > mxline) then
!                                  start a new page
!
             WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                             ,ilabel(3),ilabel(4),ilabel(5)
             iutt = ilabel (2)
             WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
             WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc)
             lcount = 2
          END IF
          jcount = 0
          DO jj = jstart, jfinal, jinc
             jcount = jcount + 1
             y (jcount) = r (ii, jj) / scale_tmp
          END DO
          WRITE (ntp, 840) ii, (y (jj), jj = 1, jcount)
       END DO
       jstart = jfinal + jinc
       istart = ibeg
    END DO main_loop
!
800 FORMAT (//, T3, '# ', I6.6, TR2, 'ut=', I6.6, TR2,     &
            'ITIME=', I5.5, '(', I2.2, ':', I2.2, ':', I2.2, ')' )
810 FORMAT (/, T4, A, '/(1.E', I2, ')', TR4, 'ut=', I6.6, TR2, &
            'ave=', ES11.4, TR2, 'sd=', ES11.4)
830 FORMAT ( / , T3, 16(TR5,i3))
840 FORMAT (T2, I3, 16(G8.3))
    RETURN
    END SUBROUTINE Outp_integer
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Outp_logical (r, ibeg, iend, iinc, jbeg, jend, jinc, &
                             xscale, ilabel, title, ntp, ncol)
    IMPLICIT NONE
    INTEGER (iprec), INTENT(IN) :: ibeg, iend, iinc, jbeg, jend, jinc, &
                                   ntp, ncol, ilabel(20)
    CHARACTER (LEN = * ), INTENT (IN) :: title
    REAL (rprec),INTENT(IN) :: xscale
    LOGICAL, INTENT(IN) :: r (:,:)
!                                                                       
!------------------------------------------------------------
!
! author: r.w. spiro        last modified: 06-25-85                     
!                                                                       
! DESCRIPTION OF INPUT PARAMETERS:                          
!      r(isize,jsize)=array to be output                                
!      this subroutine can output selected elements of array     
!      ibeg= initial i value to be output                               
!      iend= final i value to be output                                 
!      iinc= i value increment for output                               
!      (note: (iend-ibeg)/iinc should be an integer)                    
!      jbeg,jend, and jinc are defined similarly                        
!      scale=scale factor                                               
!         if(scale.eq.0.) scale is calculated to give best 
!            display      
!         all elements or array are divided by scale before
!            being output
!      ilabel= vector of length 20 that gives the label                 
!      title=character string that identifies the array 
!             being output    
!      ntp= output unit number                                          
!      ncol= number of columns for output device (80 or 132)            
!                                                                       
!------------------------------------------------------------
!                                                                       
    INTEGER (iprec), PARAMETER :: jcol = 16
    LOGICAL ::  y (jcol)
!
    REAL (rprec) :: scale_tmp, ave, sd
    INTEGER (iprec) :: mxline, jjcol, ipower, &
                       jstart, istart, ifinal, iutt, lcount, jfinal, &
                       jj, ii, jcount, isize, jsize
!
!
    isize = SIZE (r, DIM = 1)
    jsize = SIZE (r, DIM = 2)
    mxline = 57
    jjcol = (ncol - 4) / 8
    iutt   = 0
    ave    = 0
    sd     = 0
    scale_tmp = xscale
    ipower = 0
    IF (ABS(scale_tmp) > TINY(1.0_rprec)) THEN
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90:'// &
            'IF CALL OUTP_LOGICAL, SET SCALE = 0')
    END IF
    IF (ibeg > isize .OR. iend > isize .OR. &
        jbeg > jsize .OR. jend > jsize .OR. &
        ibeg < 1 .OR. iend < 1 .OR. jbeg < 1 .OR. jend < 1) THEN
        call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90:INDICES WRONG IN OUTP')
    END IF
!
!
!   Output data:
!                                                                       
    jstart = jbeg
    istart = ibeg
    ifinal = iend
!
!                 Below ilabel(6) is ITIME variable, ilabel(1) is IRDW
!                              
    WRITE (ntp,800) &
          ilabel(1),ilabel(2),ilabel(6), ilabel(3),ilabel(4),ilabel(5)
    iutt = ilabel (2)
!
    WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
    lcount = 0
!                                                                       
!   Begin output loop:
!
    main_loop: DO
       jfinal = (jjcol - 1) * jinc + jstart
       IF (jfinal > jend) jfinal = jend
       IF (jstart > jend) EXIT main_loop
       IF (mxline-lcount <= 5) then
!                               start a new page
!
          WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                          ,ilabel(3),ilabel(4),ilabel(5)
          iutt = ilabel (2)
          WRITE (ntp, 810) TRIM (title), ipower, iutt, ave, sd
          lcount = 0
       END IF
!
       WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc)
       lcount = lcount + 2
       DO ii = istart, ifinal, iinc
          lcount = lcount + 1
          IF (lcount > mxline) then
!                                  start a new page
!
             WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                             ,ilabel(3),ilabel(4),ilabel(5)
             iutt = ilabel (2)
             WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
             WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc)
             lcount = 2
          END IF
          jcount = 0
          DO jj = jstart, jfinal, jinc
             jcount = jcount + 1
             y (jcount) = r (ii, jj)
          END DO
          WRITE (ntp, 840) ii, (y (jj), jj = 1, jcount)
       END DO
       jstart = jfinal + jinc
       istart = ibeg
    END DO main_loop
!
!
800 FORMAT (//,T3,'# ', I6.6, TR2, 'ut=', I6.6, TR2,     &
            'ITIME=', I5.5, '(', I2.2, ':', I2.2,':', I2.2, ')' )
810 FORMAT (/,T4, A, '/(1.E', I2, ')', TR4,'ut=', I6.6, TR2, &
            'ave=', ES11.4, TR2, 'sd=', ES11.4)
830 FORMAT ( /, T3, 16(TR5,I3))
840 FORMAT (T2, I3, 16L8)
!
    RETURN
    END SUBROUTINE Outp_logical
!
!
!
!
    SUBROUTINE Write_real_3d_array_old (filename, rec_num, label, array_3d, setup)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    REAL (rprec),      INTENT (IN)      :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    LOGICAL           :: flag
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90:UNIT ALREADY OPEN')
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
    INQUIRE (IOLENGTH = length ) label, array_3d

    OPEN (UNIT=LUN, FILE = trim(NameRcmDir)//filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = 'UNFORMATTED',&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF

    WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat) label, array_3d
    IF (istat /= 0) THEN
       write(*,*) 'ERROR WRITING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_real_3d_array_old')
    END IF

    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) THEN
       write(*,*) 'ERROR CLOSING FILE: ',trim(NameRcmDir),filename
       call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Write_real_3d_array_old')
    END IF

    RETURN
    END SUBROUTINE Write_real_3d_array_old
!
    SUBROUTINE Read_real_3d_array_old (filename, rec_num, label, &
                                   array_3d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    REAL (rprec),      INTENT (OUT) :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: length, istat
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
!
    INQUIRE (IOLENGTH = length ) label, array_3d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_3d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_3d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_REAL_3D_ARRAY'
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF
!
    OPEN (UNIT = LUN, FILE = trim(NameRcmDir)//filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',trim(NameRcmDir),filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1, FMT = form_string) label, array_3d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1                   ) label, array_3d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file is: ', filename
       call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_real_3d_array_old
!
!
!
!
      FUNCTION Check_logical_units ()
      USE Rcm_variables, ONLY : LUN, LUN_2, LUN_3
      IMPLICIT NONE
      LOGICAL :: Check_logical_units, L1, L2
        INQUIRE (UNIT = LUN, EXIST = L1, OPENED = L2)
        IF (.NOT.L1) THEN
           Check_logical_units = .FALSE.
           RETURN
        ELSE IF (L2) THEN
           Check_logical_units = .FALSE.
           RETURN
        END IF
        INQUIRE (UNIT = LUN_2, EXIST = L1, OPENED = L2)
        IF (.NOT.L1) THEN
           Check_logical_units = .FALSE.
           RETURN
        ELSE IF (L2) THEN
           Check_logical_units = .FALSE.
           RETURN
        END IF
        INQUIRE (UNIT = LUN_3, EXIST = L1, OPENED = L2)      
        IF (.NOT.L1) THEN
           Check_logical_units = .FALSE.
           RETURN
        ELSE IF (L2) THEN
           Check_logical_units = .FALSE.
           RETURN
        END IF
        Check_logical_units = .TRUE.
      RETURN
      END FUNCTION Check_logical_units
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Read_grid ( )
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec) :: istat, isize_pr, jsize_pr
      CHARACTER (LEN=80) :: form_length

      OPEN (UNIT = LUN, STATUS = 'OLD', FORM = 'FORMATTED', &
            FILE = trim(NameRcmDir)//'input/rcmcrd11', IOSTAT = istat)
      IF (istat /= 0) THEN
         write(*,*) 'ERROR OPENING',trim(NameRcmDir)//'input/rcmcrd11'
         call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Read_grid')
      END IF
      READ (LUN, '(A80)') form_length
      READ (UNIT = LUN, FMT = form_length) isize_pr, jsize_pr, dlam, dpsi, ri, re
      CLOSE (UNIT = LUN)
!
!
      IF (isize /= isize_pr .OR. jsize /= jsize_pr ) THEN
         WRITE (*,*) ' GRID SIZES IN rcmcrd11 DO NOT MATCH THOSE IN THE CODE'
         call CON_stop('ERROR in IM/RCM2/src/rcm_io.f90')
      END IF
!
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FORM = 'FORMATTED', &
            FILE = trim(NameRcmDir)//'input/rcmcrd21', IOSTAT = istat)
      IF (istat /= 0) THEN
         WRITE (*,*)  'ERROR OPENING RCMCRD21'
         call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90::Read_grid')            
      END IF
      READ (LUN, '(A80)') form_length
      READ (UNIT=LUN, FMT = form_length) alpha
      READ (UNIT=LUN, FMT = form_length) beta
      READ (UNIT=LUN, FMT = form_length) colat
      READ (UNIT=LUN, FMT = form_length) aloct
      READ (UNIT=LUN, FMT = form_length) vcorot
      READ (UNIT=LUN, FMT = form_length) bir
      READ (UNIT=LUN, FMT = form_length) sini

      CLOSE (UNIT = LUN)
      RETURN
      END SUBROUTINE Read_grid
!
!
!
      SUBROUTINE Write_grid (isize_pr, jsize_pr )
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec) :: istat, isize_pr, jsize_pr
      CHARACTER (LEN=80) :: form_string
!
      IF (isize_pr /= isize) call CON_stop( &
           'Error in IM/RCM2/rcm_io::Write_grid: isizes do not match')

      IF (jsize_pr /= jsize) call CON_stop( &
           'Error in IM/RCM2/rcm_io::Write_grid: jsizes do not match')

!
      OPEN (UNIT = LUN, STATUS = 'REPLACE', FORM = 'FORMATTED', &
            FILE = trim(NameRcmDir)//'input/rcmcrd11', IOSTAT = istat)
      IF (istat /= 0) call  CON_stop( &
           'Error in IM/RCM2/rcm_io::Write_grid: ERROR OPENING RCMCRD11')
      form_string = '(2(TR2,I10),4(TR2,ES23.15))'
      WRITE (UNIT = LUN, FMT = '(A80)') form_string
      WRITE (UNIT = LUN, FMT=form_string) &
           isize, jsize, dlam, dpsi, Re, Ri
      CLOSE (UNIT = LUN)
!
!
      OPEN (UNIT = LUN, STATUS = 'REPLACE', FORM = 'FORMATTED', &
            FILE = trim(NameRcmDir)//'input/rcmcrd21', IOSTAT = istat)
      IF (istat /= 0) call  CON_stop( &
           'Error in IM/RCM2/rcm_io::Write_grid: ERROR OPENING RCMCRD21')

      form_string = '(3(TR2,ES23.15))'
      WRITE (LUN, '(A80)') form_string
      WRITE (LUN, form_string) alpha
      WRITE (LUN, form_string) beta
      WRITE (LUN, form_string) colat
      WRITE (LUN, form_string) aloct
      WRITE (LUN, form_string) vcorot
      WRITE (LUN, form_string) bir
      WRITE (LUN, form_string) sini
      CLOSE (UNIT = LUN)
      RETURN
      END SUBROUTINE Write_grid
!
!
!
      SUBROUTINE Read_plasma ()
      USE Rcm_variables
      INTEGER (iprec) :: n,k
      CHARACTER (LEN=80) :: form_string
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FORM = 'FORMATTED', &
           FILE = trim(NameRcmDir)//'input/rcmlas1')
!
         READ (LUN,'(TR2,I10.10)') n
         IF (n /= ksize) call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90: '// &
              'Incorrect ksize in IM/inputrcmlas1')
         READ (LUN,'(A80)') form_string
         DO k = 1, n
            READ (LUN, form_string) alam(k), eta(k), ikflav(k), fudge(k)
         END DO
!
         READ (LUN,'(TR2,I10.10)') n
         IF (n /= kcsize) call CON_STOP('ERROR in IM/RCM2/src/rcm_io.f90: '//&
              'Incorrect kcsize in IM/inputrcmlas1')
         READ (LUN,'(A80)') form_string
         DO k = 1, n
            READ (LUN, form_string) alamc(k), etac(k), ikflavc(k), fudgec(k)
         END DO
!
      CLOSE (UNIT = LUN)
!
      RETURN
      END SUBROUTINE Read_plasma
!
!
!
      SUBROUTINE Write_plasma ()
      USE Rcm_variables
      INTEGER (iprec) :: n,k
      CHARACTER (LEN=80) :: form_string
!
!
      OPEN (LUN, FORM = 'FORMATTED', STATUS='REPLACE', &
           FILE = trim(NameRcmDir)//'input/rcmlas1')
!
        form_string = '(2(TR2,ES15.7), TR2,I3, TR2,ES15.7)'
        WRITE (LUN, '(I10.10)') SIZE(alam)
        WRITE (LUN, '(A80)') form_string
        DO k = 1, SIZE (alam)
           WRITE (LUN, form_string) alam(k), eta(k), ikflav(k), fudge(k)
        END DO
!
        form_string = '(2(TR2,ES15.7), TR2,I3, TR2,ES15.7)'
        WRITE (LUN, '(I5.5)') SIZE (alamc)
        WRITE (LUN, '(A80)') form_string
        DO k = 1, SIZE (alamc)
           WRITE (LUN, form_string) alamc(k), etac(k), ikflavc(k),fudgec(k)
        END DO
!
      CLOSE (LUN)
!
      RETURN
      END SUBROUTINE Write_plasma
!
!
!
!
      FUNCTION Get_time_char_string (label) RESULT (time_string)
      USE Rcm_variables, ONLY : label_def
      IMPLICIT NONE
      TYPE (label_def), INTENT (IN) :: label
      CHARACTER (LEN=8) :: time_string
      WRITE (time_string,'(I2.2,A1,I2.2,A1,I2.2)') &
             label%intg(3),':',label%intg(4),':',label%intg(5)
      RETURN
      END FUNCTION Get_time_char_string
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Read_qtcond ()
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec) :: n, i, j, istat
      CHARACTER (LEN=80) :: form_string
      LOGICAL, SAVE :: called_already = .FALSE.
!
      IF (called_already) RETURN
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FORM = 'FORMATTED', &
            FILE = trim(NameRcmDir)//'input/rcmcond', ACTION = 'READ', IOSTAT=istat) 
      IF (istat /= 0) call CON_STOP( &
           'ERROR in IM/RCM2/rcm_io::read_qtcond: ERROR OPENING rcmcond')
!
      READ (LUN, '(I10.10)') n
      IF (n /= isize*jsize) call CON_STOP( &
           'ERROR in IM/RCM2/rcm_io::read_qtcond: sizes do not match')
      ! IF(n/=(isize)*(jsize+3)) ...
      READ (LUN,'(A80)') form_string
      DO j = 1, jsize
         DO i = 1, isize
           READ (LUN,form_string) qtplam(i,j), qthall(i,j), qtped(i,j)
        END DO
     END DO
!
!
     READ (LUN, '(I10.10)') n
     IF (n /= (jsize)) call CON_STOP(&
          'ERROR in IM/RCM2/rcm_io::read_qtcond: jsize does not match')
     READ (LUN,'(A80)') form_string
     DO j = 1, jsize
        READ (LUN,form_string) ss(j)
     END DO
!
     CALL Wrap_around_ghostcells (qtplam, isize, jsize, n_gc)
     CALL Wrap_around_ghostcells (qtped, isize, jsize, n_gc)
     CALL Wrap_around_ghostcells (qthall, isize, jsize, n_gc)
     ss (0) = ss (jsize)
     ss (-1) = ss (jsize-1)
     ss (jsize+1) = ss (1)
     ss (jsize+2) = ss (2)
!
     CLOSE (LUN)
!ss = 0.0
!qtplam = 10.
!qtped  = 10.
!qthall = 10.
!print*,'reset cond'
      called_already = .TRUE.
      RETURN
      END SUBROUTINE Read_qtcond
!
!
      SUBROUTINE Write_qtcond ()
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec) :: n, i, j, istat
      CHARACTER (LEN=80) :: form_string
!
      OPEN (UNIT = LUN, STATUS = 'REPLACE', FORM = 'FORMATTED', &
            FILE = trim(NameRcmDir)//'input/rcmcond', ACTION = 'WRITE', IOSTAT=istat) 
      IF (istat /= 0) call CON_STOP( &
           'ERROR in IM/RCM2/rcm_io::Write_qtcond: ERROR OPENING rcmcond')
!
        WRITE (LUN, '(I10.10)') isize*jsize
        form_string = '(3(TR2,ES23.15))'
        WRITE (LUN,'(A80)') form_string
        DO j = 1, jsize
        DO i = 1, isize
           WRITE (LUN,form_string) qtplam(i,j), qthall(i,j), qtped(i,j)
        END DO
        END DO
!
!
        WRITE (LUN, '(I10.10)') jsize
        form_string = '(1(TR2,ES23.15))'
        WRITE (LUN,'(A80)') form_string
        DO j = 1 , jsize
           WRITE (LUN,form_string) ss(j)
        END DO
!
      CLOSE (LUN)
      RETURN
      END SUBROUTINE Write_qtcond
!
      SUBROUTINE Wrap_around_ghostcells (array, isize, jsize, n_gc)
      USE rcm_variables, only : iprec, rprec
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: n_gc, isize, jsize
      REAL (rprec), INTENT (IN OUT) :: array (1-n_gc:isize+n_gc,1-n_gc:jsize+n_gc)
!
      array (1:isize, -1) = array (1:isize, jsize - 1)
      array (1:isize, 0)  = array (1:isize, jsize )
      array (1:isize, jsize+1) = array (1:isize, 1)
      array (1:isize, jsize+2) = array (1:isize, 2)
!
      array (-1, :) = array (1, :)
      array (0,  :) = array (1, :)
      array (isize+1, :) = array (isize,:)
      array (isize+2, :) = array (isize,:)
!
      RETURN
      END SUBROUTINE Wrap_around_ghostcells




    SUBROUTINE Read_dktime (L_dktime)
    USE Rcm_variables, ONLY : NameRcmDir, dktime, LUN
    IMPLICIT NONE
    LOGICAL, INTENT (IN) :: L_dktime
!
!
    IF (L_dktime) THEN
       OPEN (UNIT=LUN, FILE=trim(NameRcmDir)//'input/dktable', STATUS='OLD', ACTION='READ')
       READ (LUN,800) dktime
 800   FORMAT (8(1X,1PE9.3))
       CLOSE (LUN)
    ELSE
        dktime = 0.0
    END IF
    RETURN
    END SUBROUTINE Read_dktime


      SUBROUTINE Read_trf ()
      USE Rcm_variables
      IMPLICIT NONE
      INTEGER (iprec) :: i
      REAL (rprec) :: dtau
!
      OPEN (LUN, FILE=trim(NameRcmDir)//'input/trf.dat', STATUS='OLD')
      DO i = 1, SIZE(trf,1)
          READ (LUN,*) trf(i,1), trf(i,2), trf(i,3), trf(i,4)
      END DO
      CLOSE (LUN)
!
!     Now interpolate for given DOY and sunspot_number
!
      IF (doy >= 172 .AND. doy <= 355) THEN
         dtau = doy-172.0
         DO i = 1, 19
            trf(i,5) = trf(i,2) + &
             ((trf(i,3)-trf(i,2))/(355.-172.0))*dtau
         END DO
      ELSE 
         IF (DOY < 172) THEN
             dtau = doy + 10.0
         ELSE
             dtau = doy - 355.
         END IF
         DO i = 1, 19
            trf(i,5) = trf(i,3) + &
               ((trf(i,2)-trf(i,3))/183.)*dtau
         END DO
      END IF
!
      IF (sunspot_number >= 15.0 .AND. sunspot_number <= 165.) THEN
          DO i = 1, 19
             trf(i,5) = trf (i,4) + ((trf(i,5)-trf(i,4))/150.)*(sunspot_number-15.0)
          END DO
      ELSE
          IF (sunspot_number < 15.0) THEN
              trf(:,5) = trf(:,4)
          END IF
      END IF
!   
      RETURN
      END SUBROUTINE Read_trf


END MODULE Rcm_io
