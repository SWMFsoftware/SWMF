MODULE CSR_module
  !
  !  This module is part of SMLIB v. 1.1.  It contain the data-structure
  !  and subroutines for the CSR sparse matrix format.  For a detailed
  !  account of it's functionality see the documentation.
  !
  ! ------------------------------------------------------------------------
  !
  !  Copyright (C) 1996 Ernst A. Meese
  !    Refer to the file copyright.doc for details and important disclaimer.
  !
  ! Version 1.0b created 6.9.95
  ! Version 1.0.1b, validated 1.2.96 by TESTS/CSR/tester
  ! Version 1.0.2b, validated 26.2.96 by TESTS/CSR/tester
  !
  !
  !  Get the working precision
  !
  USE Precision_Module, ONLY: prec

  IMPLICIT NONE

  PRIVATE
  !
  !  Define the block size for memory allocations
  !
  INTEGER, PARAMETER :: Block_Size = 50
  !
  !  Define the derived data type for the CSR matrices
  !
  TYPE CSR
     INTEGER :: N
     REAL(prec), DIMENSION (:), POINTER :: A
     INTEGER, DIMENSION (:), POINTER :: IA, JA
  END TYPE CSR
  !
  !  Define the generic interfaces
  !
  INTERFACE operator (*)
     MODULE PROCEDURE CSR_vector_product
  END INTERFACE
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE CSR_to_array
  END INTERFACE
  INTERFACE allocate_matrix
     MODULE PROCEDURE allocate_CSR
  END INTERFACE
  INTERFACE reallocate_matrix
     MODULE PROCEDURE reallocate_CSR
  END INTERFACE
  INTERFACE deallocate_matrix
     MODULE PROCEDURE deallocate_CSR
  END INTERFACE
  INTERFACE nullify_matrix
     MODULE PROCEDURE nullify_CSR
  END INTERFACE
  INTERFACE ALLOCATED_MATRIX
     MODULE PROCEDURE allocated_CSR
  END INTERFACE
  INTERFACE ENTRY
     MODULE PROCEDURE entry_CSR
  END INTERFACE
  INTERFACE setrow
     MODULE PROCEDURE set_CSR_row
  END INTERFACE
  INTERFACE SAVE
     MODULE PROCEDURE Save_CSR
  END INTERFACE
  INTERFACE LOAD
     MODULE PROCEDURE Load_CSR
  END INTERFACE
  INTERFACE is_strictly_lower
     MODULE PROCEDURE is_strictly_lower_CSR
  END INTERFACE
  INTERFACE is_lower
     MODULE PROCEDURE is_lower_CSR
  END INTERFACE
  INTERFACE is_strictly_upper
     MODULE PROCEDURE is_strictly_upper_CSR
  END INTERFACE
  INTERFACE is_upper
     MODULE PROCEDURE is_upper_CSR
  END INTERFACE
  INTERFACE is_ok
     MODULE PROCEDURE is_ok_CSR
  END INTERFACE
  INTERFACE SetVal
     MODULE PROCEDURE SetVal_CSR
  END INTERFACE
  !
  !  Define what should be publicly available
  !
  PUBLIC :: prec, CSR, OPERATOR (*), ASSIGNMENT (=), Nullify_Matrix,&
       & Allocate_Matrix, Reallocate_Matrix, Deallocate_Matrix,&
       & Allocated_MATRIX, Entry, SetRow, Save, Load, Is_Strictly_Lower,&
       & Is_Strictly_Upper, Is_Lower, Is_Upper, Is_OK, SetVal

CONTAINS

  SUBROUTINE nullify_CSR (A)
    TYPE(CSR) :: A

    A % N = 0
    NULLIFY ( A % A, A % IA, A % JA )

  END SUBROUTINE nullify_CSR

  SUBROUTINE allocate_CSR (A, N, NZMAX)
    TYPE (CSR)                    :: A
    INTEGER,          INTENT (in) :: N
    INTEGER, OPTIONAL, INTENT(in) :: NZMAX
    INTEGER :: SIZE

    IF (N <= 0) STOP "allocate_CSR: matrix must have positive size"
    IF (PRESENT(NZMAX)) THEN
       SIZE = MAX(NZMAX,1)    ! One element is the least possible
    ELSE
       SIZE = Block_Size      ! Default memory
    END IF
    ALLOCATE ( A % A (SIZE), A % JA (SIZE), A % IA (N+1) )
    !
    !  IMPORTANT: Initially all row start pointers point to the first element
    !             of JA.  This defines all rows to be of length zero, and makes
    !             row insertion easier.
    A % IA = 1
    A % N = N

  END SUBROUTINE allocate_CSR

  SUBROUTINE reallocate_CSR (S, DELTA)
    TYPE (CSR)                         :: S         ! Matrix to reallocate
    INTEGER, OPTIONAL, INTENT (in)     :: DELTA     ! Additional memory size
    REAL(prec), DIMENSION (:), POINTER :: A         ! New storage
    INTEGER,    DIMENSION (:), POINTER :: JA        ! New storage
    INTEGER                            :: NZMAX, N
    !
    !  Note that this routine assumes the start of row pointers to be
    !  in increasing order.  This should be a matter of concern only
    !  if this routine is called at an intermediate state in the
    !  construction of a matrix where the matrix data structure is not
    !  well defined as defined in is_ok (S).  If you use only routines from 
    !  this library to construct the matrix, it should pose no problem.
    !
    IF (.NOT.ALLOCATED_CSR(S)) STOP "reallocate_CSR: matrix not allocated"
    N = S % N
    IF (PRESENT(DELTA)) THEN
       NZMAX = MAX(SIZE(S % A) + DELTA, S % IA (N+1) - 1)
    ELSE
       NZMAX = S % IA (N+1) - 1
    END IF
    ALLOCATE ( A (NZMAX), JA (NZMAX) )
    !
    !  The reallocation involves a memory copy
    !
    A (1:S % IA (N+1)-1) = S % A (1:S % IA (N+1)-1)
    JA (1:S % IA (N+1)-1) = S % JA (1:S % IA (N+1)-1)
    DEALLOCATE ( S % A, S % JA )
    S % A => A
    S % JA => JA

  END SUBROUTINE reallocate_CSR

  SUBROUTINE deallocate_CSR (A)
    TYPE(CSR) :: A

    A % N = 0
    DEALLOCATE (A % A, A % IA, A % JA)

  END SUBROUTINE deallocate_CSR

  FUNCTION allocated_CSR ( A ) result (res)
    TYPE(CSR), INTENT (in) :: A
    LOGICAL                :: res

    res = ASSOCIATED ( A % A )
    IF ((res.AND.(.NOT.ASSOCIATED(A % IA).OR..NOT.ASSOCIATED(A % JA))) .OR.&
         & (.NOT.res.AND.(ASSOCIATED(A % IA).OR.ASSOCIATED(A % JA))) .OR.&
         & (.NOT.ASSOCIATED (A % IA).AND.ASSOCIATED(A % JA)) .OR.&
         & (ASSOCIATED(A % IA).AND..NOT.ASSOCIATED(A % JA)) .OR.&
         & (res.AND.A % N <=0)) &
         & STOP "Illegal associacion in CSR, probably not NULLIFIED"

  END FUNCTION allocated_CSR

  FUNCTION entry_CSR ( A, i, j ) Result (TheEntry)
    TYPE(CSR), INTENT(inout) :: A
    INTEGER,             INTENT(in)    :: i, j
    REAL(prec)                               :: TheEntry
    INTEGER :: k
    LOGICAL :: found

    IF (A % N <= 0) STOP "entry_CSR: matrix not allocated"
    IF (i<1.OR.i>A % N) STOP "entry_CSR: row out range"
    IF (j<1.OR.j>A % N) STOP "entry_CSR: column out range"
    found = .FALSE.
    DO k = A % IA (i), A % IA (i+1) - 1
       IF (A % JA (k) == j) THEN
          found = .TRUE.
          EXIT
       END IF
    END DO
    IF (found) THEN
       TheEntry = A % A (k)
    ELSE
       TheEntry = 0.0
    END IF

  END FUNCTION entry_CSR

  SUBROUTINE set_CSR_row (A, i, cols, vals)
    TYPE (CSR),                           INTENT(inout) :: A
    INTEGER,                              INTENT(in)    :: i
    INTEGER, DIMENSION(:),                INTENT(in)    :: cols
    REAL(prec),    DIMENSION(SIZE(cols)), INTENT(in)    :: vals
    INTEGER :: old_len, new_len, diff_len, N
    !
    ! Your fair warning:  set_CSR_row do not check for multiple
    ! entries in cols due to the cost of this check.  Multiple entries
    ! in cols will give unpredictable results since the matrix will no
    ! longer be uniqely defined.
    !
    IF (.NOT.allocated_CSR(A)) STOP "set_CSR_row: matrix not allocated"
    N = A % N
    IF (i<1.OR.i>A % N) STOP "set_CSR_row: row out range"
    !
    old_len = A % IA(i+1) - A % IA(i)
    new_len = SIZE(cols)
    diff_len = new_len - old_len
    !
    !  Check to see if there is room for the new elements in A
    !
    IF ( A % IA (N+1) + diff_len > SIZE( A % A ) ) CALL reallocate_CSR (A, &
         & MAX(Block_Size, A % IA (N+1) + diff_len - SIZE(A % A)))
    !
    !  Move old values in A to fit in the new row.  Note that if the row is
    !  appended (that is; there is no row with number > i inserted) the range
    !  of the indices is empty and the next line is a 'no operation'.  Hence
    !  the fastest way to construct a matrix is to append rows successively. 
    !  This functionality is heavily dependent upon the initialisation of all
    !  start of row pointers (AJ(1:N+1)) to the base of the storeage
    !  vector (N+2)  (see the allocate_CSR routine)
    !
    A % A (A % IA (i+1)+diff_len:A % IA (N+1)+diff_len) = &
         & A % A (A % IA (i+1):A % IA (N+1))
    A % JA (A % IA (i+1)+diff_len:A % IA (N+1)+diff_len) = &
         & A % JA (A % IA (i+1):A % IA (N+1))
    !
    !  Adjust start of row pointers
    !
    A % IA (i+1:N+1) = A % IA (i+1:N+1) + diff_len
    !
    !  Insert elements
    !
    A % A (A % IA(i):A % IA(i+1)-1) = vals
    A % JA (A % IA(i):A % IA(i+1)-1) = cols

  END SUBROUTINE set_CSR_row

  FUNCTION CSR_vector_product ( A, x ) Result (y)
    TYPE(CSR),                          INTENT(in) :: A
    REAL(prec), DIMENSION (:),          INTENT(in) :: x
    REAL(prec), DIMENSION (SIZE(x))                :: y
    REAL(prec), DIMENSION (A % IA (A % N + 1) - 1) :: t
    INTEGER                                        :: i, nnz
    
    nnz = A % IA (A % N + 1) - 1
    t = A % A (1:nnz) * x ( A % JA (1:nnz) )
    y = (/ (SUM( t (A % IA (i):A % IA (i+1)-1) ), i = 1, A % N) /)

  END FUNCTION CSR_vector_product

  SUBROUTINE CSR_to_array (A, S)
    TYPE (CSR), INTENT (in) :: S
    REAL(prec), DIMENSION(:,:), INTENT(out) :: A
    INTEGER :: i, j

    A = 0.0_prec
    DO i = 1, S % N
       DO j = S % IA (i), S % IA(i+1)-1
          A(i,S % JA(j)) = S % A (j)
       END DO
    END DO

  END SUBROUTINE CSR_to_array

  FUNCTION is_strictly_lower_CSR ( S ) Result (Res)
    TYPE (CSR) :: S
    LOGICAL :: Res
    INTEGER :: i

    Res = ALL ( (/ (ALL ( S % JA ( S % IA (i):S % IA (i+1)-1 ) < i ), i = 1, S % N) /) )

  END FUNCTION is_strictly_lower_CSR

  FUNCTION is_lower_CSR ( S ) Result (Res)
    TYPE (CSR) :: S
    LOGICAL :: Res
    INTEGER :: i

    Res = ALL ( (/ (ALL ( S % JA ( S % IA (i):S % IA (i+1)-1 ) <= i ), i = 1, S % N) /) )

  END FUNCTION is_lower_CSR

  FUNCTION is_strictly_upper_CSR ( S ) Result (Res)
    TYPE (CSR) :: S
    LOGICAL :: Res
    INTEGER :: i

    Res = ALL ( (/ (ALL ( S % JA ( S % IA (i):S % IA (i+1)-1 ) > i ), i = 1, S % N) /) )

  END FUNCTION is_strictly_upper_CSR

  FUNCTION is_upper_CSR ( S ) Result (Res)
    TYPE (CSR) :: S
    LOGICAL :: Res
    INTEGER :: i

    Res = ALL ( (/ (ALL ( S % JA ( S % IA (i):S % IA (i+1)-1 ) >= i ), i = 1, S % N) /) )

  END FUNCTION is_upper_CSR

  FUNCTION is_ok_CSR ( S ) RESULT (is_ok)
    TYPE (CSR) :: S
    LOGICAL :: is_ok
    LOGICAL :: has_col ( S % N )
    INTEGER :: i, k

    IF (.not.Allocated_CSR (S)) THEN
       is_ok = .FALSE.
       RETURN
    END IF
    IF ( S % N == size ( S % IA ) - 1 .AND. size (S % A) == size (S%JA)) THEN
       is_ok = .TRUE.
    ELSE
       is_ok = .FALSE.
       RETURN
    END IF
    has_col = .FALSE.
    DO i = 1, S % N
       !
       !  Check for rows sorted by increasing order
       !
       IF ( S % IA (i) > S % IA (i+1) ) THEN
          is_ok = .FALSE.
          EXIT
       END IF
       !
       !  Check for coloumn indices out of range
       !
       IF ( ANY ( S % JA ( S % IA(i) : S % IA(i+1)-1 ) < 1) .OR. &
            ANY ( S % JA ( S % IA(i) : S % IA(i+1)-1 ) > S % N ) ) THEN
          is_ok = .FALSE.
          EXIT
       END IF
       !
       !  Check for multiple coloumn indices
       !
       DO k = S % IA(i), S % IA(i+1)-1
          IF (.NOT.has_col ( S % JA (k) )) THEN
             has_col ( S % JA (k) ) = .TRUE.
          ELSE
             is_ok = .FALSE.
          END IF
       END DO
       IF (.NOT.is_ok) EXIT
       has_col ( S % JA (S % IA(i) : S % IA(i+1)-1) ) = .FALSE.
    END DO

  END FUNCTION is_ok_CSR

  SUBROUTINE save_CSR ( A, file_name, form )

    USE IO

    TYPE (CSR),                 INTENT (in) :: A
    CHARACTER(LEN=*),           INTENT (in) :: file_name
    CHARACTER(LEN=*), OPTIONAL, INTENT (in) :: form
    CHARACTER(LEN=11)                       :: the_form
    INTEGER                                 :: iounit, ios, i
    INTEGER                                 :: Field_Len, Rec_Len, Num_Recs
    CHARACTER(LEN=80)                       :: format

    iounit = new_unit ()
    IF (iounit==-1) STOP "BIG TROUBLE: No more Fortran UNITs available."

    IF (PRESENT(form)) THEN
       the_form = ADJUSTL(form)
    ELSE
       the_form = "formatted"
    END IF

    OPEN (iounit, FILE=file_name, STATUS = "UNKNOWN", IOSTAT = ios,&
         & ACTION = "write", form = the_form)

    IF (ios/=0) THEN
       WRITE(*,"(A)") "save_CSR: Could not open write file '"// file_name //"'"
       STOP
    END IF

    SELECT CASE (the_form(1:1))
    CASE ("f", "F")
       write (iounit,*) A % N
       Field_Len = 2 + INT(LOG10(REAL(A % IA (A % N + 1)))+EPSILON(1.0E0))
       Rec_Len   = 80 / Field_Len
       Num_Recs  = CEILING (REAL(A % N) / REAL(Rec_Len))
       write (format,*) '(', Num_Recs, '(', Rec_Len, 'I', Field_Len, ':/))'
       write (iounit,format) A % IA
       Field_Len = 2 + INT(LOG10(REAL(A % N))+EPSILON(1.0E0))
       Rec_Len   = 80 / Field_Len
       Num_Recs  = CEILING (REAL(A % IA (A % N + 1)-1) / REAL(Rec_Len))
       write (format,*) '(', Num_Recs, '(', Rec_Len, 'I', Field_Len, ':/))'
       write (iounit,format) A % JA (1:A % IA (A % N + 1) - 1)
       Field_Len = 8 + Precision (1.0_prec)
       Rec_Len   = 80 / Field_Len
       Num_Recs  = CEILING (REAL(A % IA (A % N + 1)-1) / REAL(Rec_Len))
       write (format,*) '(', Num_Recs, '(', Rec_Len, 'ES', Field_Len,&
            & '.', Precision (1.0_prec), 'E2:/))'
       write (iounit,format) A % A (1:A % IA (A % N + 1) - 1)
    CASE ("u", "U")
       WRITE (iounit) A % N
       DO i = 1, A % N + 1
          WRITE (iounit) A % IA (i)
       END DO
       DO i = 1, A % IA (A % N + 1) - 1
          WRITE (iounit) A % JA (i)
       END DO
       DO i = 1, A % IA (A % N + 1) - 1
          WRITE (iounit) A % A (i)
       END DO

    END SELECT

    CLOSE(iounit, STATUS = "keep")

  END SUBROUTINE save_CSR


  SUBROUTINE load_CSR (A, file_name, form) 

    USE IO

    TYPE (CSR),                 INTENT (inout) :: A
    CHARACTER(LEN=*),           INTENT(in)     :: file_name
    CHARACTER(LEN=*), OPTIONAL, INTENT (in)    :: form
    CHARACTER(LEN=11)                          :: the_form
    INTEGER                                    :: iounit, ios, nnz, i 

    IF (ALLOCATED_CSR(A)) CALL Deallocate_Matrix ( A )

    IF (PRESENT(form)) THEN
       the_form = ADJUSTL(form)
    ELSE
       the_form = "formatted"
    END IF

    iounit = new_unit ()
    OPEN (iounit, FILE=file_name, STATUS = "OLD", IOSTAT = ios,&
         & ACTION = "READ", form = the_form)

    IF (ios/=0) THEN
       WRITE(*,"(A)") "load_CSR: Could not open read file '" // file_name //"'"
       RETURN
    END IF
    
    SELECT CASE (the_form(1:1))
    CASE ("f","F")
       READ (iounit, *) A % N 
       IF (A % N <= 0) STOP "Load_CSR: TRYING TO LOAD EMPTY MATRIX"
       ALLOCATE ( A % IA (A % N+1) )
       READ (iounit, *) A % IA
       nnz = A % IA( A % N + 1) - 1
       IF (nnz <= 0) STOP "Load_CSR: MATRIX CONTAINS NO ELEMENTS" 
       ALLOCATE ( A % JA (nnz), A % A (nnz) )
       READ (iounit, *) A % JA
       READ (iounit, *) A % A
    CASE ("u","U")
       READ (iounit) A % N
       ALLOCATE ( A % IA (A % N + 1) )
       DO i = 1, A % N + 1
          READ (iounit) A % IA (i)
       END DO
       nnz = A % IA (A % N + 1) - 1; ALLOCATE ( A % JA(nnz), A % A(nnz) )
       DO i = 1, A % IA (A % N + 1) - 1
          READ (iounit) A % JA (i)
       END DO
       DO i = 1, A % IA (A % N + 1) - 1
          READ (iounit) A % A (i)
       END DO
    END SELECT

    CLOSE (iounit) 
    
  END SUBROUTINE load_CSR

  SUBROUTINE setVal_CSR ( S, i, j, val ) 
    TYPE(CSR),  INTENT(INOUT) :: S
    INTEGER,    INTENT(IN)    :: i, j
    REAL(prec), INTENT(IN)    :: val
    INTEGER                   :: k
    
    DO k = S % IA ( i ), S % IA ( i + 1 ) - 1
       if ( S % JA ( k ) == j ) EXIT
    END DO
    IF ( k < S % IA ( i + 1 ) ) THEN
       S % A ( k ) = val
    ELSE
       STOP "Setting value violating CSR matrix sparsity pattern."
    END IF

  END SUBROUTINE setVal_CSR

END MODULE CSR_module
