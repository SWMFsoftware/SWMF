MODULE ILU_MODULE
  !
  !  This module is part of SMLIB v. 1.1.  It contains procedures needed
  !  to perform ILU(p) factorisations on sparse matrices stored in the MSR
  !  format.
  !
  ! ------------------------------------------------------------------------
  !
  !  Copyright (C) 1996 Ernst A. Meese
  !    Refer to the file copyright.doc for details and important disclaimer.
  !
  !  Created: January 1996
  !  Development version
  !
  ! ------------------------------------------------------------------------
  !  Module contents
  !
  !   Procedures:
  !      Symbolic_Cancellation   - Create a symbolic cancellation fitting
  !                                a given sparse matrix.
  !      ILU                     - Calculate the sparse matrix factorisation.
  !
  ! -----------------------------------------------------------------------
  !
  !  We use the modules defining abstract datatypes for the matrices
  !
  USE Matrix_Arithmetic_Module

  IMPLICIT NONE

  SAVE

  INTERFACE Symbolic_Cancellation
     MODULE PROCEDURE Symbolic_Cancellation_ilu
  END INTERFACE

CONTAINS

  SUBROUTINE Symbolic_Cancellation_ilu (p, S, LU)
    !
    !  Purpose: Perform a symbolic cancellation on the matrix S to 
    !           create the data structure for the incomplete
    !           factorisation LU.
    !
    ! ---------------------------------------------------------------
    !
    !  Argument list:
    !
    !    p  (in)  - Maximum level of fill-in (integer)
    !    S  (in)  - A matrix (TYPE(MSR))
    !    LU (out) - The symbolic factors (TYPE(LU_CSR_MSR))
    !
    !  Local Variables:
    !
    !    i          - Matrix row index.
    !    k          - Element number to cancel on row i.
    !    ik         - LU % L % JA (ik) corresponds to matrix element ik
    !                 (cf. the algorithm in the documentation).
    !    kj         - LU % L % JA (kj) or LU % U % JA (kj) corresponds to matrix
    !                 element kj (cf. the algorithm in the documentation).
    !    l          - Dummy index.
    !    level      - Level of fill-in calculated for the latest matrix element.
    !    CurrCol    - Current matrix coloumn index.
    !    Sorter     - Array of indices to LU % L % JA on the current row
    !                 under cancellation so that 
    !                 LU % L % JA (Sorter (1:Sorter_Len)) is sorted in 
    !                 increasing order.
    !    Sorter_Len - Number of elements in Sorter.
    !    NonZero    - Array of 1-1 correspondence with a matrix row (lenght S % N).
    !                 If on the current row in the matrix, element j is non-zero,
    !                 NonZero(j) points to the corresponding element in LU % L % JA
    !                 of LU % U % JA.  If element j is zero, NonZero(j) = 0.
    !    Mem_Block_Size - Everytime memory allocated to either LU % L or LU % U is
    !                 extended, it is extended in blocks of this size.
    !
    !  Internal procedures:
    !    Reallocate_L - Increases the memory of LU % L by Mem_Block_Size
    !    Reallocate_U - Increases the memory of LU % U by Mem_Block_Size
    !    Sort_Sorter  - Sorts the level zero elements in LU % L before the
    !                   cancellation begins.
    !
    ! ---------------------------------------------------------------
    !
    !  Argument list variables
    !
    INTEGER, INTENT(in) :: p
    TYPE (MSR)          :: S
    TYPE (LU_CSR_MSR)   :: LU
    !
    !  Local variables
    !
    INTEGER                    :: i, k, l, ik, kj, level, CurrCol, Sorter_Len, k2,k3
    INTEGER, DIMENSION (S % N) :: Sorter, NonZero
    INTEGER, PARAMETER         :: Mem_Block_Size = 100
    integer :: test

    !
    !  Checking allocation status for LU.  The fitness of S should
    !  be checked with is_ok(S) before this procedure is called.
    !
    IF (.NOT.ALLOCATED_MATRIX(LU))&
         & CALL Allocate_Matrix (LU, S % N, S % N, S % JA (2) )
    IF (LU % U % N /= S % N .OR. SIZE( LU % U % A ) < S % JA (2)) THEN
       CALL Deallocate_Matrix ( LU )
       CALL Allocate_Matrix (LU, S % N, S % N, S % JA (2) )
    END IF
    !
    !  Initiate arrays.  We use the % A arrays to store the element levels.
    !
    LU % L % A = 0.0_prec; LU % U % A = 0.0_prec; NonZero = 0
    !
    !  The first line of L is empty, and the first line of U will equal the
    !  first line of S
    !
    LU % L % IA (1) = 1; LU % L % IA (2) = 1
    LU % U % JA (1) = S % JA (1); LU % U % JA (2) = S % JA (2)
    LU % U % JA ( S % JA (1):S % JA (2) - 1 ) = S % JA (S % JA (1):S % JA (2) - 1)
    !
    !  Main loop, cancellation on matrix row i
    !
    DO i = 2, S % N
       !
       ! Insert level zero elements into L and set up Sorter
       !
       k = COUNT (S % JA (S % JA (i):S % JA (i+1)-1) < i)
       LU % L % IA (i+1) = LU % L % IA (i) + k
       IF (LU % L % IA (i+1) >= SIZE (LU % L % A)) CALL Reallocate_L
!\\\
       k3 = LU % L % IA (i)
       do k2 = S % JA (i), S % JA (i+1) - 1
          if (S % JA (k2) < i) then
             LU % L % JA (k3) = S % JA (k2)
             k3=k3+1
          end if
       end do
!---
!!$       LU % L % JA (LU % L % IA (i):LU % L % IA (i+1)-1) = &
!!$            & PACK(S % JA (S % JA (i):S % JA (i+1) - 1), S % JA (S % JA (i):S % JA (i+1) - 1) < i)
!///
       Sorter_Len = k
       Sorter(1:Sorter_Len) = (/ (k, k = LU % L % IA (i), LU % L % IA (i+1) - 1) /)
       CALL Sort_Sorter
       !
       ! Insert level zero elements into U
       !
       k = COUNT (S % JA (S % JA (i):S % JA (i+1) - 1) > i)
       LU % U % JA (i+1) = LU % U % JA (i) + k

!       write(6,*) 'i : ',i, s%n,S % JA (i),S % JA (i+1)

       IF (LU % U % JA (i+1) >= SIZE (LU % U % A)) CALL Reallocate_U

!       write(6,*) 'After'

       if (LU % U % JA (i+1) > LU % U % JA (i)) then

!\\\
          k3 = LU % U % JA (i)
          do k2 = S % JA (i), S % JA (i+1) - 1
             if (S % JA (k2) > i) then
                LU % U % JA (k3) = S % JA (k2)
                k3=k3+1
             end if
          end do
!---
!!$          LU % U % JA (LU % U % JA (i):LU % U % JA (i+1)-1) = &
!!$               & PACK(S % JA (S % JA (i):S % JA (i+1) - 1), S % JA (S % JA (i):S % JA (i+1) - 1) > i)
!///

       endif

!       write(6,*) 'after after'
       !
       ! Set the NonZero pointer for the level zero elements
       !
       NonZero (i) = i
!\\\
       do k=LU % L % IA (i), LU % L % IA (i+1) - 1
          NonZero (LU % L % JA (k)) = k
       end do
       do k=LU % U % JA (i), LU % U % JA (i+1) - 1
          NonZero (LU % U % JA (k)) = k
       end do
!---
!!$       NonZero ((/ (LU % L % JA (k), k = LU % L % IA (i), LU % L % IA (i+1) - 1) /)) = &
!!$            (/ (k, k = LU % L % IA (i), LU % L % IA (i+1) - 1) /)
!!$       NonZero ((/ (LU % U % JA (k), k = LU % U % JA (i), LU % U % JA (i+1) - 1) /)) = &
!!$            (/ (k, k = LU % U % JA (i), LU % U % JA (i+1) - 1) /)
!///
       !
       ! Cancel elements up to level p in current row i
       ! 
       k = 1
       DO
          IF (k > Sorter_Len) EXIT
          ik = Sorter (k)
          kj = LU % U % JA ( LU % L % JA (ik) )
          !
          ! Create fill-in elements
          !
          DO
             IF (kj >= LU % U % JA ( LU % L % JA (ik) + 1)) EXIT 
             CurrCol = LU % U % JA (kj)
             level = NINT(LU % L % A (ik) + LU % U % A (kj)) + 1
             IF (CurrCol < i) THEN
                IF (NonZero (CurrCol) == 0) THEN
                   IF (level <= p) THEN
                      LU % L % JA ( LU % L % IA (i+1) ) = CurrCol
                      LU % L % A  ( LU % L % IA (i+1) ) = level
                      NonZero (CurrCol) =  LU % L % IA (i+1) 
                      DO l = k + 1, Sorter_Len
                         IF ( CurrCol < LU % L % JA (Sorter (l)) ) EXIT
                      END DO
                      Sorter (l+1:Sorter_Len+1) = Sorter (l:Sorter_Len)
                      Sorter (l) = LU % L % IA (i+1)
                      Sorter_Len = Sorter_Len + 1
                      LU % L % IA (i+1) = LU % L % IA (i+1) + 1
                      IF (LU % L % IA (i+1) >= SIZE (LU % L % A)) CALL Reallocate_L
                   END IF
                ELSE
                   LU % L % A (NonZero (CurrCol)) = &
                        & MIN (REAL(level,prec), LU % L % A (NonZero (CurrCol)))
                END IF
             ELSE
                IF (NonZero (CurrCol) == 0) THEN
                   IF (level <= p) THEN
                      LU % U % JA (LU % U % JA (i+1)) = CurrCol
                      LU % U % A  (LU % U % JA (i+1)) = level
                      NonZero (CurrCol) = LU % U % JA (i+1)
                      LU % U % JA (i+1) = LU  % U % JA (i+1) + 1
                      IF (LU % U % JA (i+1) >= SIZE (LU % U % A)) CALL Reallocate_U
                   END IF
                ELSE
                   LU % U % A (NonZero (CurrCol)) = &
                        & MIN (REAL(level,prec), LU % U % A (NonZero (CurrCol)))
                END IF
             END IF
             kj = kj + 1
          END DO
          k = k + 1
       END DO
       !
       !  Sparse set to zero
       !
       NonZero (i) = 0
!\\\
       do k=LU % L % IA (i), LU % L % IA (i+1) - 1
          NonZero (LU % L % JA (k)) = 0
       end do
       do k=LU % U % JA (i), LU % U % JA (i+1) - 1
          NonZero (LU % U % JA (k)) = 0
       end do
!---
!!$       NonZero ((/ (LU % L % JA (k), k = LU % L % IA (i), LU % L % IA (i+1) - 1) /)) = 0
!!$       NonZero ((/ (LU % U % JA (k), k = LU % U % JA (i), LU % U % JA (i+1) - 1) /)) = 0
!///
    END DO

  CONTAINS
    SUBROUTINE Reallocate_L 
      REAL(prec), DIMENSION(:), POINTER :: A
      INTEGER,    DIMENSION(:), POINTER :: JA
      
      ALLOCATE (  A ( SIZE(LU % L % A ) + Mem_Block_Size ) )
      ALLOCATE ( JA ( SIZE(LU % L % JA) + Mem_Block_Size ) )
      A (1:SIZE(LU % L % A )) = LU % L % A  ;  A( SIZE(LU % L %  A ) + 1:SIZE( A) ) = 0.0_prec
      JA(1:SIZE(LU % L % JA)) = LU % L % JA
      DEALLOCATE ( LU % L % A,  LU % L % JA )
      LU % L % A  => A
      LU % L % JA => JA
    END SUBROUTINE Reallocate_L

    SUBROUTINE Reallocate_U 
      REAL(prec), DIMENSION(:), POINTER :: A
      INTEGER,    DIMENSION(:), POINTER :: JA
      
      ALLOCATE (  A ( SIZE(LU % U % A ) + Mem_Block_Size ) )
      ALLOCATE ( JA ( SIZE(LU % U % JA) + Mem_Block_Size ) )
      A (1:SIZE(LU % U % A )) = LU % U % A  ;  A( SIZE(LU % U %  A ) + 1:SIZE( A) ) = 0.0_prec
      JA(1:SIZE(LU % U % JA)) = LU % U % JA
      DEALLOCATE ( LU % U % A,  LU % U % JA )
      LU % U % A  => A
      LU % U % JA => JA
   END SUBROUTINE Reallocate_U

    SUBROUTINE Sort_Sorter
      INTEGER :: i, j, TEMP

      DO i = 1, Sorter_Len - 1
         DO j = i + 1, Sorter_Len
            IF (LU % L % JA (Sorter(i)) > LU % L % JA (Sorter(j))) THEN
               TEMP = Sorter (i); Sorter (i) = Sorter (j); Sorter (j) = TEMP
            END IF
         END DO
      END DO

    END SUBROUTINE Sort_Sorter

  END SUBROUTINE Symbolic_Cancellation_ilu

  SUBROUTINE ILU (S, LU)
    !
    !  Purpose: Factorises a sparse matrix into L and U factors by the
    !           level of fill-in version of incomplete LU factorisation.
    !           The LU datastructure must previously have been 
    !           constructed by the Symbolic_Cancellation procedure.
    !
    ! ---------------------------------------------------------------
    !
    !  Argument list:
    !
    !    S  (in)  - A matrix (TYPE(MSR))
    !    LU (out) - The symbolic factors (TYPE(LU_CSR_MSR))
    !
    !  Local Variables:
    !
    !    i          - Matrix row index.
    !    ik         - LU % L % JA (ik) corresponds to matrix element ik
    !                 (cf. the algorithm in the documentation).
    !    kj         - LU % L % JA (kj) or LU % U % JA (kj) corresponds to matrix
    !                 element kj (cf. the algorithm in the documentation).
    !    ij         - LU % L % JA (ij) or LU % U % JA (ij) corresponds to matrix
    !                 element ij (cf. the algorithm in the documentation).
    !    k          - Dummy index.
    !    Sorter     - Array of indices to LU % L % JA on the current row
    !                 under cancellation so that 
    !                 LU % L % JA (Sorter (1:Sorter_Len)) is sorted in 
    !                 increasing order.
    !    Sorter_Len - Number of elements in Sorter.
    !    NonZero    - Array of 1-1 correspondence with a matrix row (lenght S % N).
    !                 If on the current row in the matrix, element j is non-zero,
    !                 NonZero(j) points to the corresponding element in LU % L % JA
    !                 of LU % U % JA.  If element j is zero, NonZero(j) = 0.
    !
    ! ---------------------------------------------------------------
    !
    !  Argument list variables
    !
    TYPE (MSR)          :: S
    TYPE (LU_CSR_MSR)   :: LU
    !
    !  Local variables
    !
    INTEGER                    :: i, ik, kj, ij, k, Sorter_Len, k2,k3
    INTEGER, DIMENSION( S % N) :: Sorter, NonZero

    NonZero = 0
    LU % L % A ( 1 : LU % L % IA ( LU % L % N + 1 ) - 1) = Real(0, prec)
    LU % U % A ( 1 : LU % U % JA ( LU % U % N + 1 ) - 1) = Real(0, prec)
    LU % U % A (1) = S % A (1)
    LU % U % A (LU % U % JA(1):LU % U % JA(2) - 1) = S % A (S % JA(1):S % JA(2) - 1) 
    DO i = 2, S % N

!\\\
       k=0
       do k2=LU % L % IA(i),LU % L % IA(i+1)-1
          NonZero (LU % L % JA (k2)) = LU % L % IA(i)+k
          k=k+1
       end do
       NonZero(i) = i
       k=0
       do k2=LU % U % JA(i),LU % U % JA(i+1)-1
          NonZero (LU % U % JA (k2)) = LU % U % JA(i)+k
          k=k+1
       end do
!---
!!$       NonZero ((/ LU % L % JA (LU % L % IA(i):LU % L % IA(i+1)-1), i, &
!!$            &      LU % U % JA (LU % U % JA(i):LU % U % JA(i+1)-1) /)) = &
!!$            &   (/ (LU % L % IA(i)+k, k=0, LU % L % IA (i+1) - LU % L % IA (i) - 1), i, &
!!$            &      (LU % U % JA(i)+k, k=0, LU % U % JA (i+1) - LU % U % JA (i) - 1) /)
!///
       LU % U % A (i) = S % A (i)
       k = COUNT (S % JA (S % JA (i):S % JA (i+1)-1) < i)
!\\\
       k3 = LU % L % IA (i)
       do k2 = S % JA (i), S % JA (i+1) - 1
          if (S % JA (k2) < i) then
             LU % L % A (k3) = S % A (k2)
             k3=k3+1
          end if
       end do
!---
!!$       LU % L % A (LU % L % IA (i):LU % L % IA (i)+k-1) = &
!!$            & PACK(S % A (S % JA (i):S % JA (i+1) - 1), S % JA (S % JA (i):S % JA (i+1) - 1) < i)
!///
       Sorter_Len = LU % L % IA (i+1) - LU % L % IA (i)
       Sorter(1:Sorter_Len) = (/ (k, k = LU % L % IA (i), LU % L % IA (i+1) - 1) /)
       CALL Sort_Sorter
       k = COUNT (S % JA (S % JA (i):S % JA (i+1)-1) > i)
!\\\
       if (k>=1) then
          k3 = LU % U % JA (i)
          do k2 = S % JA (i), S % JA (i+1) - 1
             if (S % JA (k2) > i ) then
                LU % U % A (k3) = S % A (k2)
                k3=k3+1
             end if
          end do
       endif
!---
!!$       if (LU % U % JA (i)+k-1 >= LU % U % JA (i)) then
!!$          LU % U % A (LU % U % JA (i):LU % U % JA (i)+k-1) = &
!!$               & PACK(S % A (S % JA (i):S % JA (i+1) - 1), S % JA (S % JA (i):S % JA (i+1) - 1) > i)
!!$       endif
!///
       DO k = 1, Sorter_Len
          ik = Sorter(k)
          LU % L % A (ik) = LU % L % A (ik) / LU % U % A (LU % L % JA (ik))
          DO kj = LU % U % JA (LU % L % JA (ik)), LU % U % JA (LU % L % JA (ik)+1) - 1
             ij = NonZero(LU % U % JA (kj))
             IF ( ij /= 0) THEN
               IF (LU % U % JA (kj) < i) THEN
                   LU % L % A (ij) = LU % L % A (ij) - LU % L % A (ik) * LU % U % A (kj)
                ELSE
                   LU % U % A (ij) = LU % U % A (ij) - LU % L % A (ik) * LU % U % A (kj)
                END IF
             END IF
          END DO
       END DO

!\\\
       k=0
       do k2=LU % L % IA(i),LU % L % IA(i+1)-1
          NonZero (LU % L % JA (k2)) = 0
          k=k+1
       end do
       NonZero(i) = 0
       k=0
       do k2=LU % U % JA(i),LU % U % JA(i+1)-1
          NonZero (LU % U % JA (k2)) = 0
          k=k+1
       end do
!---
!!$       NonZero ((/ LU % L % JA (LU % L % IA(i):LU % L % IA(i+1)-1), i, &
!!$            &      LU % U % JA (LU % U % JA(i):LU % U % JA(i+1)-1) /)) = 0
!///
    END DO

  CONTAINS
    SUBROUTINE Sort_Sorter
      INTEGER :: i, j, TEMP
      
      DO i = 1, Sorter_Len - 1
         DO j = i + 1, Sorter_Len
            IF (LU % L % JA (Sorter(i)) > LU % L % JA (Sorter(j))) THEN
               TEMP = Sorter (i); Sorter (i) = Sorter (j); Sorter (j) = TEMP
            END IF
         END DO
      END DO
    END SUBROUTINE Sort_Sorter

  END SUBROUTINE ILU

END MODULE ILU_MODULE
