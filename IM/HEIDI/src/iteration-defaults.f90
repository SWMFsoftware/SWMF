MODULE Iteration_Defaults_Module
  !
  !  This module is part of SMLIB v. 1.1.  It contain data, data structures, 
  !  and procedures common to all iterative solvers.  This includes default
  !  values and convergence history maintainance.
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
  !   Parameters:
  !       NONE, SOME, ALL - Defining the amount of information the iterative
  !                         solvers should print out.
  !       Default_tolerance - Default stop residual thereshold
  !       Default_minimum_iterations - Default least number of iterations
  !       Default_maximum_iterations - Default largest number of iterations
  !       Default_infolevel - Default amout of information to print out
  !       Default_use_init_value - Use initial value or not
  !       F1 - Formatting string for iteration output
  !       F2 - Formatting string for iteration output
  !       F3 - Formatting string for iteration output
  ! 
  !   Data types:
  !       ConvHistItem - Storing information about one iteration
  !       ConvHist     - Linked list of ConvHistItem's
  !
  !   Interfaces:
  !       Write_ConvHist - Generic interface for  write_history_to_screen, 
  !                        write_history_to_file, and write_history_to_string
  !
  !   Subroutines:
  !       Nullify_ConvHist        - Nullify pointers and initate variables
  !       Initiate_ConvHist       - Reset the counters
  !       Reallocate_ConvHist     - Create more memory if needed
  !       Deallocate_ConvHist     - Release memory used by the linked list
  !       SetIterationHistory     - Set data for some iteration number
  !       write_history_to_screen - Write the linked list to screen
  !       write_history_to_file   - Write the linked list to file
  !       write_history_to_string - Write the linked list to a string
  !
  !   Functions:
  !       No_of_its         - Returns the total number of iterations
  !       Time_ConvHist     - Returns the CPU-time used at some iteration
  !       Flop_ConvHist     - Returns the Flops used at some iteration
  !       Residual_ConvHist - Returns the residual at some iteration
  !
  ! -----------------------------------------------------------------------
  !
  USE Precision_Module, ONLY: prec    !  Get the library working real precision.
  USE Flop_Module                     !  On some Crays, it is possible to count
                                      !  float point operations.

  IMPLICIT NONE

  SAVE

  PRIVATE

  INTEGER, PARAMETER :: NONE = 0, SOME = 1, ALL = 2
  REAL(prec)         :: Default_tolerance = 0.000001_prec
  INTEGER            :: Default_minimum_iterations = 0, Default_maximum_iterations = 50
  INTEGER            :: Default_infolevel = NONE
  LOGICAL            :: Default_use_init_value = .FALSE.
  !
  ! Some common format strings
  !
  CHARACTER(LEN=*), PARAMETER :: &
       F1 = '(A,"; tolerance:", 1PE9.2, ", min it.:", I4, ", max it.:", I4)', &
       F2 = '(2X,"No of its.:", I5, ", 2-norm of calc. residual:", 1PE11.4)', &
       F3 = '(2X,"It. no.:", I4,"; 2-norm of calc. residual:",1PE10.3,"; 2-norm of res.:",1PE10.3)'
  !
  ! Data type for iteration information
  !
  TYPE ConvHistItem
     PRIVATE
     REAL              :: Time
     INTEGER(Flop_Int) :: Flop
     REAL(prec)        :: Residual
  END TYPE ConvHistItem
  !
  ! Data type for the convergence history list
  !
  TYPE ConvHist
     PRIVATE
     INTEGER           :: N, Last_Iter
     INTEGER(Flop_Int) :: Init_Flop
     REAL              :: Init_Time
     TYPE(ConvHistItem), DIMENSION(:), POINTER :: Item
  END TYPE ConvHist
  !
  ! Generic procedure interface for output for history data
  !
  INTERFACE Write_ConvHist
     MODULE PROCEDURE write_history_to_screen
     MODULE PROCEDURE write_history_to_file
     MODULE PROCEDURE write_history_to_string
  END INTERFACE
  !
  ! Declaring PUBLIC entries
  !
  PUBLIC :: ConvHist, Write_ConvHist, Initiate_ConvHist, SetIterationHistory, Nullify_ConvHist, &
       & NONE, SOME, ALL, F1, F2, F3, &
       & Default_tolerance, Default_minimum_iterations, Default_maximum_iterations, &
       & Default_infolevel, Default_use_init_value

CONTAINS
  
  Subroutine Nullify_ConvHist (History)
    TYPE(ConvHist) :: History

    History % N = -1; History % Last_Iter = -1
    NULLIFY ( History % Item )

  END Subroutine Nullify_ConvHist

  SUBROUTINE Reallocate_ConvHist (History)
    TYPE (ConvHist)                            :: History
    TYPE (ConvHistItem), DIMENSION(:), POINTER :: New_Item
    REAL, PARAMETER :: Scale = 1.5
    
    ALLOCATE ( New_Item (INT(Scale*History % N)))
    New_Item (1:History % N) = History % Item
    DEALLOCATE ( History % Item )
    History % Item => New_Item
    History % N = INT(Scale*History % N)

  END SUBROUTINE Reallocate_ConvHist

  SUBROUTINE Deallocate_ConvHist (History)
    TYPE(ConvHist) :: History
 
    History % N = -1; History % Last_Iter = -1
    DEALLOCATE (History % Item)

  END SUBROUTINE Deallocate_ConvHist

  SUBROUTINE Initiate_ConvHist (History, MAX)
    TYPE(ConvHist)      :: History
    INTEGER, INTENT(in) :: MAX
    INTERFACE CPU_TIME
       SUBROUTINE CPU_TIME (TIME)
         REAL :: TIME
       END SUBROUTINE CPU_TIME
    END INTERFACE

    IF (History % N /= MAX .and. History % N /= -1) CALL deallocate_ConvHist (History)
    History % Last_Iter = -1; History % N = MAX
    IF (.NOT.ASSOCIATED(History % Item)) ALLOCATE ( History % Item(0:MAX) )

    History % Init_Flop = dflop ()
!    CALL CPU_TIME ( History % Init_Time )    

  END SUBROUTINE Initiate_ConvHist

  SUBROUTINE SetIterationHistory ( history, residual )
    REAL(prec),     INTENT(in)    :: residual
    TYPE(ConvHist), INTENT(inout) :: history
    REAL :: Time
    INTERFACE CPU_TIME
       SUBROUTINE CPU_TIME (TIME)
         REAL :: TIME
       END SUBROUTINE CPU_TIME
    END INTERFACE

    History % Last_Iter = History % Last_Iter + 1
    IF (History % Last_Iter > History % N) Call Reallocate_ConvHist ( History )
!    CALL CPU_TIME (Time)
    History % Item (History % Last_Iter) = ConvHistItem (Time, dflop(), residual)

  END SUBROUTINE SetIterationHistory

  SUBROUTINE write_history_to_screen (History)
    TYPE(ConvHist), INTENT(in) :: History
    INTEGER :: i

    IF (Flop_Counter_Exists) THEN
       WRITE (*,"(A4,A8,A15,A11)") "  IC", "    t[s]", "           Flop", "   Residual"
    ELSE
       WRITE (*,"(A4,A8,A11)") "  IC", "    t[s]", "   Residual"
    END IF
    DO i = 0, History % Last_Iter
       IF (Flop_Counter_Exists) THEN
          WRITE (*,"(I4,F8.2,I15,1PE11.3)") i, &
               History % Item (i) % Time - History % Init_Time, &
               & History % Item (i) % Flop, History % Item (i) % Residual
       ELSE
          WRITE (*,"(I4,F8.2,1PE11.3)") i, &
               History % Item (i) % Time - History % Init_Time, &
               History % Item (i) % Residual
       END IF
    END DO
    
  END SUBROUTINE write_history_to_screen

  SUBROUTINE write_history_to_file (iounit, History)
    TYPE(ConvHist), INTENT(in) :: History
    INTEGER,        INTENT(in) :: iounit
    INTEGER :: i

    IF (Flop_Counter_Exists) THEN
       WRITE (iounit,"(A4,A8,A15,A11)") "  IC", "    t[s]", "           Flop", "   Residual"
    ELSE
       WRITE (iounit,"(A4,A8,A11)") "  IC", "    t[s]", "   Residual"
    END IF
    DO i = 0, History % Last_Iter
       IF (Flop_Counter_Exists) THEN
          WRITE (iounit,"(I4,F8.2,I15,1PE11.3)") i, &
               History % Item (i) % Time - History % Init_Time,&
               & History % Item (i) % Flop, History % Item (i) % Residual
       ELSE
          WRITE (iounit,"(I4,F8.2,1PE11.3)") i, &
               History % Item (i) % Time - History % Init_Time, &
               History % Item (i) % Residual
       END IF
    END DO
    
  END SUBROUTINE write_history_to_file

  SUBROUTINE write_history_to_string ( ic, History, Line )
    INTEGER,        INTENT (in)  :: ic
    TYPE(ConvHist), INTENT (in)  :: History
    CHARACTER(*),   INTENT (out) :: Line

    IF (ic>History % Last_Iter .or. ic<0) THEN
       Line = " **** Out of History Bounds ****"
    ELSE
       IF (Flop_Counter_Exists) THEN
          WRITE (Line,"(I4,F8.2,I15,1PE11.3)") ic, &
               History % Item (ic) % Time - History % Init_Time,&
               & History % Item (ic) % Flop, History % Item (ic) % Residual
       ELSE
          WRITE (Line,"(I4,F8.2,1PE11.3)") ic, &
               History % Item (ic) % Time - History % Init_Time, &
               History % Item (ic) % Residual
       END IF
    END IF

  END SUBROUTINE Write_history_to_string

  INTEGER FUNCTION No_of_Its ( History )
    TYPE (ConvHist) :: History

    No_of_Its = History % Last_Iter

  END FUNCTION No_of_Its

  FUNCTION Time_ConvHist (ic, History )
    TYPE (Convhist), INTENT (in) :: History
    INTEGER,         INTENT (in) :: ic
    REAL(prec) :: Time_ConvHist

    Time_ConvHist = REAL(History % Item (ic) % Time - History % Init_Time, prec)

  END FUNCTION Time_ConvHist

  FUNCTION Flop_ConvHist ( ic, History )
    TYPE (Convhist), INTENT (in) :: History
    INTEGER,         INTENT (in) :: ic
    INTEGER(flop_int) :: Flop_ConvHist

    IF (Flop_Counter_Exists) THEN
       Flop_ConvHist = History % Item (ic) % Flop - History % Init_Flop
    ELSE
       Flop_ConvHist = -1
    END IF

  END FUNCTION Flop_ConvHist

  FUNCTION Residual_ConvHist (ic, History )
    TYPE (Convhist), INTENT (in) :: History
    INTEGER,         INTENT (in) :: ic
    REAL(prec) :: Residual_ConvHist

    Residual_ConvHist = History % Item (ic) % Residual

  END FUNCTION Residual_ConvHist

END MODULE Iteration_Defaults_Module
