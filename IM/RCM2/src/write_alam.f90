   PROGRAM Write_alam_values
   IMPLICIT NONE
!
   INTEGER :: k, kmin, kmax
   REAL :: kc, alam_0, dx, y
   REAL, ALLOCATABLE :: alam(:)
!
   kmin = 1
   kmax = 100
   alam_0 = 202.48
   y = 6.409
   dx = 2.0*y/(kmax-kmin)
   kc = (REAL(kmin)+REAL(kmax))*0.5
   ALLOCATE (alam(kmin:kmax))
!
   OPEN (UNIT=1, FILE='alam_values.dat', STATUS='REPLACE')
!
   DO k = kmin, kmax
      IF (k < kc) THEN
         alam(k) = alam_0*EXP(dx*(k-kc))
      ELSE
         alam(k) = alam_0*(1.0+dx*(k-kc))
      END IF
      WRITE (1,*) k, alam(k), 2, 0.0
   END DO
   DEALLOCATE (alam)
!
!
   kmin = 101
   kmax = 150
   alam_0 = -202.48/7.0
   y = 6.409
   dx = 2.0*y/(kmax-kmin)
   kc = (REAL(kmin)+REAL(kmax))*0.5
   ALLOCATE (alam(kmin:kmax))
   DO k = kmin, kmax
      IF (k < kc) THEN
         alam(k) = alam_0*EXP(dx*(k-kc))
      ELSE
         alam(k) = alam_0*(1.0+dx*(k-kc))
      END IF
      WRITE (1,*) k, alam(k), 1, 0.333
   END DO

!
!
   CLOSE (1)
   DEALLOCATE (alam)
!
   STOP
   END PROGRAM Write_alam_values
