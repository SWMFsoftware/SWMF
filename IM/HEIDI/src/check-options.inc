IF (PRESENT(no_of_its)) no_of_its = ic
IF (PRESENT(residual_norm)) residual_norm = norm

CONTAINS

SUBROUTINE Check_Options ()

  IF ( PRESENT ( init ) ) THEN
     IF (.not.init) X = 0.0_prec
  ELSE
     IF (.not.Default_use_init_value) X = 0.0_prec
  END IF

  IF (PRESENT(reduction_factor) ) toler=reduction_factor * SQRT(SUM((b - A*x)**2))

  IF ( PRESENT ( tol ) ) THEN
     toler = tol
  ELSE
     IF(.NOT.PRESENT(reduction_factor)) toler = Default_tolerance
  END IF

  IF ( PRESENT ( min_it ) ) THEN
     min_iter = min_it
  ELSE
     min_iter = Default_minimum_iterations
  END IF

  IF ( PRESENT ( max_it ) ) THEN
     max_iter = max_it
  ELSE
     max_iter = Default_maximum_iterations
  END IF

  IF ( PRESENT ( INFO ) ) THEN
     infolevel = INFO
  ELSE
     infolevel = Default_infolevel
  END IF

  IF ( PRESENT ( HISTORY ) ) Call Initiate_ConvHist ( History, max_iter )

END SUBROUTINE Check_Options
