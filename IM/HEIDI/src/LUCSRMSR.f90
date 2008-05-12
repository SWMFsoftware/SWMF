MODULE LU_CSR_MSR_module
  !
  !  This module is part of SMLIB v. 1.1.  It contain the data-structure
  !  and subroutines for the LUCSRMSR sparse matrix factorisation format.
  !  For a detailed account of it's functionality, see the documentation.
  !
  ! ------------------------------------------------------------------------
  !
  !  Copyright (C) 1996 Ernst A. Meese
  !    Refer to the file copyright.doc for details and important disclaimer.
  !
  USE CSR_Module
  USE MSR_Module

  IMPLICIT NONE

  PRIVATE

  TYPE LU_CSR_MSR
     TYPE (CSR) :: L
     TYPE (MSR) :: U
  END TYPE LU_CSR_MSR

  INTERFACE Nullify_Matrix
     MODULE PROCEDURE Nullify_LU_CSR_MSR
  END INTERFACE
  INTERFACE Allocate_Matrix
     MODULE PROCEDURE Allocate_LU_CSR_MSR
  END INTERFACE
  INTERFACE Reallocate_Matrix
     MODULE PROCEDURE Reallocate_LU_CSR_MSR
  END INTERFACE
  INTERFACE Deallocate_Matrix
     MODULE PROCEDURE Deallocate_LU_CSR_MSR
  END INTERFACE
  INTERFACE ALLOCATED_MATRIX
     MODULE PROCEDURE Allocated_LU_CSR_MSR
  END INTERFACE
  INTERFACE Entry_L
     MODULE PROCEDURE Entry_sparse_L
  END INTERFACE
  INTERFACE Entry_U
     MODULE PROCEDURE Entry_sparse_U
  END INTERFACE
  INTERFACE Solve
     MODULE PROCEDURE LU_CSR_MSR_Solve
  END INTERFACE
  INTERFACE is_ok
     MODULE PROCEDURE is_ok_LU_CSR_MSR
  END INTERFACE

  PUBLIC :: LU_CSR_MSR, Nullify_Matrix, Allocate_Matrix, Reallocate_Matrix, Deallocate_Matrix, ALLOCATED_MATRIX, &
       &    Entry_L, Entry_U, Solve, is_ok

CONTAINS

  SUBROUTINE Nullify_LU_CSR_MSR (LU)
    TYPE(LU_CSR_MSR), INTENT(inout) :: LU

    CALL Nullify_Matrix ( LU % L )
    CALL Nullify_Matrix ( LU % U )

  END SUBROUTINE Nullify_LU_CSR_MSR

  SUBROUTINE Allocate_LU_CSR_MSR (LU, N, LNZMAX, UNZMAX)
    TYPE (LU_CSR_MSR) :: LU
    INTEGER, INTENT (in) :: N
    INTEGER, OPTIONAL, INTENT(in) :: LNZMAX, UNZMAX

    IF (PRESENT (LNZMAX)) THEN
       CALL Allocate_Matrix ( LU % L, N, LNZMAX )
    ELSE
       CALL Allocate_Matrix ( LU % L, N )
    END IF
    IF (PRESENT (UNZMAX)) THEN
       CALL Allocate_Matrix ( LU % U, N, UNZMAX )
    ELSE
       CALL Allocate_Matrix ( LU % U, N )
    END IF

  END SUBROUTINE Allocate_LU_CSR_MSR

  SUBROUTINE Reallocate_LU_CSR_MSR (LU, LDELTA, UDELTA)
    TYPE (LU_CSR_MSR) :: LU
    INTEGER, OPTIONAL, INTENT(in) :: LDELTA, UDELTA

    IF (PRESENT (LDELTA)) THEN
       CALL Reallocate_Matrix ( LU % L, LDELTA )
    ELSE
       CALL Reallocate_Matrix ( LU % L )
    END IF
    IF (PRESENT (UDELTA)) THEN
       CALL Reallocate_Matrix ( LU % U, UDELTA )
    ELSE
       CALL Reallocate_Matrix ( LU % U )
    END IF

  END SUBROUTINE Reallocate_LU_CSR_MSR

  SUBROUTINE Deallocate_LU_CSR_MSR (LU)
    TYPE(LU_CSR_MSR), INTENT(inout) :: LU

    CALL Deallocate_Matrix ( LU % L )
    CALL Deallocate_Matrix ( LU % U )

  END SUBROUTINE Deallocate_LU_CSR_MSR

  FUNCTION Allocated_LU_CSR_MSR (LU) Result (ALL)
    TYPE(LU_CSR_MSR), INTENT(in) :: LU
    LOGICAL :: ALL

    ALL = ALLOCATED_MATRIX ( LU % L) .AND. ALLOCATED_MATRIX ( LU % U )

  END FUNCTION Allocated_LU_CSR_MSR

  FUNCTION Entry_sparse_L ( LU, i, j ) Result (TheEntry)
    TYPE(LU_CSR_MSR), INTENT(inout) :: LU
    INTEGER,          INTENT(in)    :: i, j
    REAL(prec)                      :: TheEntry

    IF ( i == j) THEN
       TheEntry = 1.0_prec
    ELSE
       TheEntry = ENTRY ( LU % L, i, j )
    END IF

  END FUNCTION Entry_sparse_L

  FUNCTION Entry_sparse_U ( LU, i, j ) Result (TheEntry)
    TYPE(LU_CSR_MSR), INTENT(inout) :: LU
    INTEGER,          INTENT(in)    :: i, j
    REAL(prec)                      :: TheEntry

    TheEntry = ENTRY ( LU % U, i, j)

  END FUNCTION entry_sparse_U

  SUBROUTINE LU_CSR_MSR_solve (LU, b, x)
    TYPE (LU_CSR_MSR),        INTENT (in)  :: LU
    REAL(prec), DIMENSION(LU % L % N), INTENT (in)  :: b
    REAL(prec), DIMENSION(LU % L % N), INTENT (out) :: x
    REAL(prec), DIMENSION(LU % L % N) :: y
    INTEGER :: i, i1, i2
    !
    ! Your fair warning: This routine do not make any checks as to L being
    ! strictly lower diagonal and U being upper diagonal.  Too query for this
    ! use the function is_ok ( LU ).
    !
    ! Forward substitution
    !
    i2 = 0
    DO i = 1, LU % L % N
       i1 = i2 + 1; i2 = LU % L % IA (i+1)-1
       y(i) = b(i) - DOT_PRODUCT(LU % L % A (i1:i2), y(LU % L % JA (i1:i2)))
    END DO
    !
    ! Backward substitution
    !
    DO i = LU % U % N, 1, -1
       i1 = LU % U % JA (i); i2 = LU % U % JA (i+1) - 1
       x(i) = ( y(i) & 
                - DOT_PRODUCT (LU % U % A(i1:i2), x(LU % U % JA(i1:i2))) &
              ) / LU % U % A (i)
    END DO

  END SUBROUTINE LU_CSR_MSR_solve

  FUNCTION is_ok_LU_CSR_MSR ( LU ) Result (is_ok_res)
    TYPE (LU_CSR_MSR) :: LU
    LOGICAL :: is_ok_res

    IF (.not.Allocated_MATRIX (LU)) THEN
       is_ok_res = .FALSE.
       RETURN
    END IF
    is_ok_res = is_strictly_lower ( LU % L ) .AND. is_upper ( LU % U ) &
         .AND. is_ok ( LU % L ) .AND. is_ok ( LU % U ) &
         .AND. LU % L % N == LU % U % N

  END FUNCTION is_ok_LU_CSR_MSR

END MODULE LU_CSR_MSR_module
