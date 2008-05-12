MODULE Band_LU_Module
  !
  !  This module is part of SMLIB v. 1.1.  It contain the data-structure
  !  and subroutines for the Band_LU banded matrix LU format.  For a detailed
  !  account of it's functionality, see the documentation.
  !
  ! ------------------------------------------------------------------------
  !
  !  Copyright (C) 1998 Ernst A. Meese
  !    Refer to the file copyright.doc for details and important disclaimer.
  !
  ! Created: 28.5.98 by Ernst A. Meese
  !
  USE Precision_module

  IMPLICIT NONE

  PRIVATE

  TYPE Band_LU
     INTEGER                             :: N, Upper_band_width, Lower_band_width
     REAL(prec), DIMENSION(:,:), POINTER :: abd
     INTEGER, DIMENSION(:),      POINTER :: ipvt
  END TYPE Band_LU

  INTERFACE Allocate_Matrix
     MODULE PROCEDURE Allocate_Band_LU
  END INTERFACE

  INTERFACE Deallocate_Matrix
     MODULE PROCEDURE Deallocate_Band_LU
  END INTERFACE

  INTERFACE Nullify_Matrix
     MODULE PROCEDURE Nullify_Band_LU
  END INTERFACE

  INTERFACE ALLOCATED_Matrix
     MODULE PROCEDURE ALLOCATED_Band_LU
  END INTERFACE

  PUBLIC :: Band_LU, Allocate_Matrix, Deallocate_Matrix, Nullify_Matrix, ALLOCATED_Matrix

CONTAINS

  SUBROUTINE Allocate_Band_LU (A, Lower_band_width, Upper_band_width, N)
    TYPE(Band_LU), INTENT(INOUT) :: A
    INTEGER           :: Lower_band_width, Upper_band_width, N

    A % N = N
    A % Lower_band_width = Lower_band_width
    A % Upper_band_width = Upper_band_width
    ALLOCATE ( A % abd (2*Lower_band_width+Upper_band_width+1,N), A % ipvt (N)  )

  END SUBROUTINE Allocate_Band_LU

  SUBROUTINE Deallocate_Band_LU ( A )
    TYPE(Band_LU), INTENT(INOUT) :: A

    A % N = 0; A % Lower_band_width = 0; A % Upper_band_width = 0
    DEALLOCATE ( A % abd, A % ipvt )

  END SUBROUTINE Deallocate_Band_LU

  SUBROUTINE Nullify_Band_LU ( A )
    TYPE(Band_LU), INTENT(INOUT) :: A

    A % N = 0; A % Lower_band_width = 0; A % Upper_band_width = 0
    NULLIFY ( A % abd); NULLIFY ( A % ipvt )

  END SUBROUTINE Nullify_Band_LU

  FUNCTION ALLOCATED_Band_LU ( A ) RESULT ( ALLOC )
    TYPE(Band_LU) :: A
    LOGICAL :: ALLOC

    ALLOC = ASSOCIATED(A%abd) .and. ASSOCIATED(A%ipvt)

  END FUNCTION ALLOCATED_Band_LU

END MODULE Band_LU_Module
