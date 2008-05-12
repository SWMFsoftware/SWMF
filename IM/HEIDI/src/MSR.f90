MODULE MSR_module
  !
  !  This module is part of SMLIB v. 1.1.  It contain the data-structure
  !  and subroutines for the MSR sparse matrix format.  For a detailed
  !  account of it's functionality see the documentation.
  !
  ! ------------------------------------------------------------------------
  !
  !  Copyright (C) 1996 Ernst A. Meese
  !    Refer to the file copyright.doc for details and important disclaimer.
  !
  ! Created: 30.1.96 by Ernst A. Meese
  ! Version 1.0b, validated by TESTS/MSR/tester 1.2.96
  ! Version 1.0.1b, validated by TESTS/MSR/tester 3.9.97
  !
  USE Precision_Module, ONLY : prec
  !DIR$ MODINLINE

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: Block_Size = 50

  TYPE MSR
     INTEGER                            :: N
     REAL(prec), DIMENSION (:), POINTER :: A
     INTEGER,    DIMENSION (:), POINTER :: JA
  END TYPE MSR

  INTERFACE operator (*)
     MODULE PROCEDURE MSR_vector_product
  END INTERFACE
  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE MSR_to_array
  END INTERFACE
  INTERFACE allocate_matrix
     MODULE PROCEDURE allocate_MSR
  END INTERFACE
  INTERFACE reallocate_matrix
     MODULE PROCEDURE reallocate_MSR
  END INTERFACE
  INTERFACE deallocate_matrix
     MODULE PROCEDURE deallocate_MSR
  END INTERFACE
  INTERFACE nullify_matrix
     MODULE PROCEDURE nullify_MSR
  END INTERFACE
  INTERFACE allocated_matrix
     MODULE PROCEDURE allocated_MSR
  END INTERFACE
  INTERFACE ENTRY
     MODULE PROCEDURE entry_MSR
  END INTERFACE
  INTERFACE setrow
     MODULE PROCEDURE set_MSR_row
  END INTERFACE
  INTERFACE is_lower
     MODULE PROCEDURE is_lower_MSR
  END INTERFACE
  INTERFACE is_upper
     MODULE PROCEDURE is_upper_MSR
  END INTERFACE
  INTERFACE is_ok
     MODULE PROCEDURE is_ok_MSR
  END INTERFACE
  INTERFACE setVal
     MODULE PROCEDURE setVal_MSR
  END INTERFACE
  INTERFACE addVal
     MODULE PROCEDURE addVal_MSR
  END INTERFACE
  INTERFACE rowSum
     MODULE PROCEDURE rowSum_MSR
  END INTERFACE
  INTERFACE SAVE
     MODULE PROCEDURE Save_MSR
  END INTERFACE

  PUBLIC :: MSR, operator (*), ASSIGNMENT (=), allocate_matrix, reallocate_matrix, deallocate_matrix, nullify_matrix, &
       &    allocated_matrix, ENTRY, setrow, is_lower, is_upper, is_ok, setVal, addVal, rowSum, RowNNZ, SAVE

CONTAINS

  SUBROUTINE nullify_MSR (A)
    TYPE(MSR) :: A

    A % N = 0
    NULLIFY ( A % A, A % JA )

  END SUBROUTINE nullify_MSR

  SUBROUTINE allocate_MSR (A, N, NZMAX)
    TYPE (MSR)                     :: A
    INTEGER,           INTENT (in) :: N
    INTEGER, OPTIONAL, INTENT (in) :: NZMAX
    INTEGER :: SIZE

    IF (N <= 0) STOP "allocate_MSR: matrix must have positive size"
    IF (PRESENT (NZMAX)) THEN
       SIZE = MAX(NZMAX+1,N+2)
    ELSE
       SIZE = N + Block_Size + 1
    END IF
    ALLOCATE ( A % A (SIZE), A % JA (SIZE) )
    A % JA (1:N+1) = N+2
    A % A (N+1:N+2) = 0.0_prec
    A % N = N

  END SUBROUTINE allocate_MSR

  SUBROUTINE reallocate_MSR (S, DELTA)
    TYPE (MSR)                     :: S
    INTEGER, OPTIONAL, INTENT (in) :: DELTA
    REAL(prec), DIMENSION (:), POINTER :: A
    INTEGER,    DIMENSION (:), POINTER :: JA
    INTEGER                            :: NZMAX, N
    !
    !  Note that this routine assumes the start of row pointers to be
    !  in increasing order.  This should be a matter of concern only
    !  if this routine is called at an intermediate state in the
    !  construction of a matrix where the matrix data structure is not
    !  well defined as defined in is_ok (S). If you use only routines from 
    !  this library to construct the matrix, it should pose no problem.
    !
    IF (.NOT.ALLOCATED_MSR(S)) STOP "reallocate_MSR: matrix not allocated"
    N = S % N
    IF (PRESENT(DELTA)) THEN
       NZMAX = MAX(SIZE(S % A) + DELTA, S % N + 2, S % JA (N+1) )
    ELSE
       NZMAX = S % JA (N+1) 
    END IF
    ALLOCATE ( A (NZMAX), JA (NZMAX) )
    A (1:S % JA (N+1)-1) = S % A (1:S % JA (N+1)-1)
    JA (1:S % JA (N+1)-1) = S % JA (1:S % JA (N+1)-1)
    DEALLOCATE ( S % A, S % JA )
    S % A => A
    S % JA => JA

  END SUBROUTINE reallocate_MSR

  SUBROUTINE deallocate_MSR (A)
    TYPE(MSR) :: A

    A % N = 0
    DEALLOCATE (A % A, A % JA)

  END SUBROUTINE deallocate_MSR

  FUNCTION allocated_MSR ( A ) result (res)
    TYPE(MSR), INTENT (in) :: A
    LOGICAL                :: res

    res = ASSOCIATED ( A % A ) .AND. ASSOCIATED ( A % JA )

  END FUNCTION allocated_MSR

  FUNCTION entry_MSR ( A, i, j ) Result (TheEntry)
    TYPE(MSR), INTENT(inout) :: A
    INTEGER,             INTENT(in)    :: i, j
    REAL(prec)                               :: TheEntry
    INTEGER :: k

    IF (A % N <= 0) STOP "entry_MSR: matrix not allocated"
    IF (i<1.OR.i>A % N) STOP "entry_MSR: row out range"
    IF (j<1.OR.j>A % N) STOP "entry_MSR: column out range"
    IF ( i == j) THEN
       TheEntry = A % A (i)
    ELSE
       TheEntry = 0.0
       DO k = A % JA (i), A % JA (i+1) - 1
          IF (A % JA (k) == j) THEN
             TheEntry = A % A (k)
             EXIT
          END IF
       END DO
    END IF

  END FUNCTION entry_MSR

  SUBROUTINE set_MSR_row (A, i, cols, vals)
    TYPE (MSR),                           INTENT(inout) :: A
    INTEGER,                              INTENT(in)    :: i
    INTEGER, DIMENSION(:),                INTENT(in)    :: cols
    REAL(prec),    DIMENSION(SIZE(cols)), INTENT(in)    :: vals
    INTEGER :: old_len, new_len, diff_len, N, j, k
    !
    ! Your fair warning:  set_MSR_row do not check for multiple
    ! entries in cols due to the cost of this check.  Multiple entries
    ! in cols will give unpredictable results since the matrix will no
    ! longer be uniqely defined.
    !
    IF (.NOT.allocated_MSR(A)) STOP "set_MSR_row: matrix not allocated"
    N = A % N
    IF (i<1.OR.i>N) STOP "set_MSR_row: row out range"
    !
    !  Treat diagonal element by itself.  k = position in cols with 
    !  diagonal element of A.  k = 0 implies no diagonal element.
    !
    k = 0
    DO j = 1, SIZE(cols)
       IF (i == cols(j)) THEN
          k = j
          EXIT
       END IF
    END DO
    IF ( k > 0) THEN
       A % A (i) = vals(k)
    ELSE
       A % A (i) = 0.0_prec
    END IF
    !
    old_len = A % JA(i+1) - A % JA(i)
    new_len = SIZE(cols); IF (k > 0) new_len = new_len - 1
    diff_len = new_len - old_len
    !
    !  Check to see if there is room for the new elements in A
    !
    IF ( A % JA (N+1) + diff_len > SIZE( A % A ) ) CALL reallocate_MSR (A, &
         & MAX(Block_Size, A % JA (N+1) + diff_len - SIZE(A % A)))
    !
    !  Move old values in A to fit in the new row.  Note that if the row is
    !  appended (that is; there is no row with number > i inserted) the range
    !  of the indices is empty and the next line is a 'no operation'.  Hence
    !  the fastest way to construct a matrix is to append rows successively. 
    !  This functionality is heavily dependent upon the initialisation of all
    !  start of row pointers (AJ(1:N+1)) to the base of the storeage
    !  vector (N+2)  (see the allocate_MSR routine)
    !
    A % A (A % JA (i+1)+diff_len:A % JA (N+1)+diff_len) =&
         & A % A (A % JA (i+1):A % JA (N+1))
    A % JA (A % JA (i+1)+diff_len:A % JA (N+1)+diff_len) =&
         & A % JA (A % JA (i+1):A % JA (N+1))
    !
    !  Adjust start of row pointers
    !
    A % JA (i+1:N+1) = A % JA (i+1:N+1) + diff_len
    !
    !  Insert off-diagonal elements
    !
    A % A (A % JA(i):A % JA(i+1)-1) = (/ vals(1:k-1), vals(k+1:SIZE(vals)) /)
    A % JA (A % JA(i):A % JA(i+1)-1) = (/ cols(1:k-1), cols(k+1:SIZE(cols)) /)

  END SUBROUTINE set_MSR_row

!!$  FUNCTION MSR_vector_product ( A, x ) Result (y)
!!$    TYPE(MSR),                         INTENT(in) :: A
!!$    REAL(prec), DIMENSION (:),         INTENT(in) :: x
!!$    REAL(prec), DIMENSION (SIZE(x))               :: y
!!$    REAL(prec), DIMENSION (A%JA(1):A%JA(A%N+1)-1) :: t
!!$    INTEGER                                       :: i, last
!!$    
!!$    last = A % JA (A % N + 1) - 1
!!$    t = A % A (A%JA(1):last) * x ( A % JA (A%JA(1):last) )
!!$    y = A % A (1:A%N) * x + (/ ( SUM( t (A % JA (i):A % JA (i+1)-1) ), i = 1, A%N) /)
!!$
!!$  END FUNCTION MSR_vector_product

  FUNCTION MSR_vector_product ( A, x ) Result (y)
    TYPE(MSR),                          INTENT(in) :: A
    REAL(prec), DIMENSION (A%N),        INTENT(in) :: x
    REAL(prec), DIMENSION (A%N)                    :: y
    INTEGER                                        :: i, k, k_max, idx
    !
    ! This routine is optimized on a Cray-J90
    !
    k_max = MAXVAL( A % JA (2:A % N + 1) - A % JA (1:A%N) ) - 1
    y = A % A(1:A%N) * x

    DO k = 0, k_max
       !DIR$ IVDEP
       DO i = 1, A % N
          IF ( k >= A % JA (i+1) - A % JA (i)) CYCLE
          idx = A % JA (i) + k
          y(i) = y(i) + A % A ( idx ) * x( A % JA (idx) )
       END DO
    END DO

  END FUNCTION MSR_vector_product


  SUBROUTINE MSR_to_array (A, S)
    TYPE (MSR), INTENT (in) :: S
    REAL(prec), DIMENSION(:,:), INTENT(out) :: A
    INTEGER :: i, j

    A = 0.0_prec
    DO i = 1, S % N
       A (i,i) = S % A (i)
       DO j = S % JA (i), S % JA(i+1)-1
          A(i,S % JA(j)) = S % A (j)
       END DO
    END DO

  END SUBROUTINE MSR_to_array

  FUNCTION is_lower_MSR ( S ) Result (Res)
    TYPE (MSR) :: S
    LOGICAL :: Res
    INTEGER :: i

    Res = ALL ( (/ (ALL ( S % JA ( S % JA (i):S % JA (i+1)-1 ) <= i ), i = 1, S % N) /) )

  END FUNCTION is_lower_MSR

  FUNCTION is_upper_MSR ( S ) Result (Res)
    TYPE (MSR) :: S
    LOGICAL :: Res
    INTEGER :: i

    Res = ALL ( (/ (ALL ( S % JA ( S % JA (i):S % JA (i+1)-1 ) >= i ), i = 1, S % N) /) )

  END FUNCTION is_upper_MSR

  FUNCTION is_ok_MSR ( S ) RESULT (is_ok)
    TYPE (MSR) :: S
    LOGICAL :: is_ok
    LOGICAL :: has_col ( S % N )
    INTEGER :: i, k

    IF (.not.Allocated_MSR (S)) THEN
       is_ok = .FALSE.
       RETURN
    END IF
    IF ( S % N == S % JA (1) - 2 .AND. size (S % A) == size (S % JA) ) THEN
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
       IF ( S % JA (i) > S % JA (i+1) ) THEN
          is_ok = .FALSE.
          EXIT
       END IF
       !
       !  Check for coloumn indices out of range
       !
       IF ( ANY ( S % JA ( S % JA(i) : S % JA(i+1)-1 ) < 1) .OR. &
            ANY ( S % JA ( S % JA(i) : S % JA(i+1)-1 ) > S % N ) ) THEN
          is_ok = .FALSE.
          EXIT
       END IF
       !
       !  Check for multiple coloumn indices
       !
       DO k = S % JA(i), S % JA(i+1)-1
          IF (.NOT.has_col ( S % JA (k) )) THEN
             has_col ( S % JA (k) ) = .TRUE.
          ELSE
             is_ok = .FALSE.
          END IF
       END DO
       IF (.NOT.is_ok) EXIT
       has_col ( S % JA (S % JA(i) : S % JA(i+1)-1) ) = .FALSE.
    END DO

  END FUNCTION is_ok_MSR

  SUBROUTINE setVal_MSR ( S, i, j, val ) 
    TYPE(MSR),  INTENT(INOUT) :: S
    INTEGER,    INTENT(IN)    :: i, j
    REAL(prec), INTENT(IN)    :: val
    INTEGER                   :: k
    
    IF ( i == j ) THEN
       S % A ( i ) = val
    ELSE
       DO k = S % JA (i), S % JA ( i + 1 ) - 1
          if ( S % JA ( k ) == j ) EXIT
       END DO
       IF ( k < S % JA ( i + 1 ) ) THEN
          S % A ( k ) = val
       ELSE
          STOP "Setting value violating MSR matrix sparsity pattern."
       END IF
    END IF

  END SUBROUTINE setVal_MSR
  
  SUBROUTINE addVal_MSR ( S, i, j, val ) 
    TYPE(MSR),  INTENT(INOUT) :: S
    INTEGER,    INTENT(IN)    :: i, j
    REAL(prec), INTENT(IN)    :: val
    INTEGER                   :: k
    
    IF ( i == j ) THEN
       S % A ( i ) = S % A ( i ) + val
    ELSE
       DO k = S % JA (i), S % JA ( i + 1 ) - 1
          if ( S % JA ( k ) == j ) EXIT
       END DO
       IF ( k < S % JA ( i + 1 ) ) THEN
          S % A ( k ) = S % A ( k ) + val
       ELSE
          STOP "Setting value violating MSR matrix sparsity pattern."
       END IF
    END IF

  END SUBROUTINE addVal_MSR

  FUNCTION rowSum_MSR ( S, i ) RESULT ( sum )
    TYPE(MSR),  INTENT(IN) :: S
    INTEGER,    INTENT(IN) :: i
    REAL(prec)             :: sum
    INTEGER                :: k
    
    sum = S % A ( i )
    DO k = S % JA (i), S % JA ( i + 1 ) - 1
       sum = sum + S % A ( k )
    END DO
     
  END FUNCTION rowSum_MSR

  FUNCTION RowNNZ( S, i )
    TYPE(MSR),  INTENT(IN) :: S
    INTEGER,    INTENT(IN) :: i
    INTEGER                :: RowNNZ

    RowNNZ =  S % JA (i+1) - S % JA (i) + 1

  END FUNCTION RowNNZ

  SUBROUTINE save_MSR ( A, file_name, form )
    !
    !  Writes MSR to a file in CSR format
    !
    USE IO
    !
    TYPE (MSR),                 INTENT (in) :: A
    CHARACTER(LEN=*),           INTENT (in) :: file_name
    CHARACTER(LEN=*), OPTIONAL, INTENT (in) :: form
    CHARACTER(LEN=11)                       :: the_form
    INTEGER                                 :: iounit, ios, i, k
    INTEGER                                 :: Field_Len, Rec_Len, Num_Recs
    CHARACTER(LEN=80)                       :: format

    iounit = new_unit ()
    IF (iounit==-1) STOP "BIG TROUBLE: No more Fortran UNITs available."

    IF (PRESENT(form)) THEN
       the_form = ADJUSTL(form)
    ELSE
       the_form = "formatted"
    END IF

    OPEN (iounit, FILE=file_name, STATUS = "UNKNOWN", IOSTAT = ios, ACTION = "write", form = the_form)

    IF (ios/=0) THEN
       WRITE(*,"(A)") "save_MSR: Could not open write file '"// file_name //"'"
       STOP
    END IF

    SELECT CASE (the_form(1:1))
    CASE ("f", "F")
       write (iounit,*) A % N
       Field_Len = 2 + INT(LOG10(REAL(A % JA (A % N + 1) + A % N))+EPSILON(1.0E0))
       Rec_Len   = 80 / Field_Len
       Num_Recs  = CEILING (REAL(A % N) / REAL(Rec_Len))
       write (format,*) '(', Num_Recs, '(', Rec_Len, 'I', Field_Len, ':/))'
       write (iounit,format) (/ (A % JA(i) - A % N - 1 + i - 1, i = 1, A % N + 1) /)
       Field_Len = 2 + INT(LOG10(REAL(A % N))+EPSILON(1.0E0))
       Rec_Len   = 80 / Field_Len
       Num_Recs  = CEILING (REAL(A % JA (A % N + 1)-2) / REAL(Rec_Len))
       write (format,*) '(', Num_Recs, '(', Rec_Len, 'I', Field_Len, ':/))'
       write (iounit,format) (/ ( (/ i, A % JA ( A % JA(i):A % JA (i + 1) - 1) /), i = 1, A % N ) /)
       Field_Len = 8 + Precision (1.0_prec)
       Rec_Len   = 80 / Field_Len
       Num_Recs  = CEILING (REAL(A % JA (A % N + 1)-2) / REAL(Rec_Len))
       write (format,*) '(', Num_Recs, '(', Rec_Len, 'ES', Field_Len,&
            & '.', Precision (1.0_prec), 'E2:/))'
       write (iounit,format) (/ ( (/ A % A(i), A % A ( A % JA(i):A % JA (i + 1) - 1) /), i = 1, A % N ) /)
    CASE ("u", "U")
       WRITE (iounit) A % N
       DO i = 1, A % N + 1
          WRITE (iounit) A % JA(i) - A % N - 1 + i - 1
       END DO
       DO i = 1, A % N
          WRITE(iounit) i
          DO k = A % JA(i), A % JA (i + 1) - 1
             WRITE (iounit) A % JA (k)
          END DO
       END DO
       DO i = 1, A % N
          WRITE(iounit) A % A(i)
          DO k = A % JA(i), A % JA (i + 1) - 1
             WRITE (iounit) A % A (k)
          END DO
       END DO

    END SELECT

    CLOSE(iounit, STATUS = "keep")

  END SUBROUTINE save_MSR

END MODULE MSR_module
