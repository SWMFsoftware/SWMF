MODULE Flop_Module
  !
  !  This module is part of SMLIB v. 1.1.  It defines dummy data and procedures for computers
  !  having no floating point operation counting mechanism (that is; most computers).
  !
  ! ---------------------------------------------------------------------------------------------
  !
  !  Copyright (C) 1996 Ernst A. Meese
  !    Refer to the file copyright.doc for details and important disclaimer.
  !
  !  Created: January 1996
  !  Development version
  !
  ! ---------------------------------------------------------------------------------------------
  !  Module contents
  !
  !   Parameters:
  !      Flop_Counter_Exists (Logical) - Set to .FALSE. to indicate the non-existence of such a
  !                                      counter.
  !   Data Types:
  !      Flop_Int - Smalles possible integer kind to save memory.
  !
  !   Procedures:
  !      init_flop - Returns a warning about the non-existence of a flop counter.
  !      dflop     - Returns -1. 
  !      flop      - Returns -1.
  !
  ! --------------------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  LOGICAL, PARAMETER :: Flop_Counter_Exists = .FALSE.
  INTEGER, PARAMETER :: flop_int = KIND(1)

CONTAINS
  SUBROUTINE init_flop () 

    WRITE (*,*) 'Warning init_flop: A flop counter is not implemented on this system.'

  END SUBROUTINE init_flop

  FUNCTION flop ()
    INTEGER(flop_int) :: flop

    flop = -1

  END FUNCTION flop

  FUNCTION dflop ()
    INTEGER(flop_int) :: dflop

    dflop = -1

  END FUNCTION dflop

END MODULE Flop_Module
