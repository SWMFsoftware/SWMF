      SUBROUTINE Bc2 (maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &               dx,dy,q,maux,aux,t,dt,mthbc)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: maxmx, maxmy, meqn, mbc, mx, my, maux
      INTEGER, INTENT (IN) :: mthbc (4)
      DOUBLE PRECISION, INTENT (IN) :: xlower, ylower, dx, dy, t, dt
      DOUBLE PRECISION, INTENT (IN OUT) ::
     &     q   (1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn),
     &     aux (1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
!
!
c     =====================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my)
c     # to the ghost cells outside the region:
c     #   (i, 1-jbc)   for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (i, my+jbc)  for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (1-ibc, j)   for ibc = 1,mbc,  j = 1-mbc, my+mbc
c     #   (mx+ibc, j)  for ibc = 1,mbc,  j = 1-mbc, my+mbc
c
      INTEGER :: ibc, jbc, i, j    
!
!
      SELECT CASE (mthbc(1)+1)   ! # left boundary:
!
      CASE (1)   !  # user-specified boundary conditions 
!
      CASE (2)   !  # zero-order extrapolation:
!
         DO ibc=1,mbc
            q   (1-ibc,:,:) = q   (1,:,:)
            aux (1-ibc,:,:) = aux (1,:,:)
         END DO
!
      CASE (3)   !  # periodic:  
!
         DO ibc=1,mbc
            q   (1-ibc,:,:) = q   (mx+1-ibc,:,:)
            aux (1-ibc,:,:) = aux (mx+1-ibc,:,:)
         END DO
!
      CASE (4) ! # solid wall (assumes 2nd comp. is vel. or momentum in x):
!
         DO ibc=1,mbc
         DO j = 1-mbc, my+mbc
            q   (1-ibc,j,:) = q   (ibc,j,:)
         END DO
         END DO
!
!       # negate the normal velocity:
         DO ibc=1,mbc
         DO j = 1-mbc, my+mbc
            q(1-ibc,j,2) = -q(ibc,j,2)
         END DO
         END DO
!
      CASE DEFAULT
         STOP 'ILLEGAL VALUE OF MTHBC(1)'
      END SELECT
! 
!
      SELECT CASE (mthbc(2)+1) ! # right boundary:
!
      CASE (1)
!
!        # user-specified boundary conditions
!
      CASE (2)
!
!       # zero-order extrapolation:
!
         DO ibc=1,mbc
            q   (mx+ibc,:,:) = q   (mx,:,:)
            aux (mx+ibc,:,:) = aux (mx,:,:)
         END DO
!
      CASE (3)
!
!        # periodic:  
!
         DO ibc=1,mbc
            q   (mx+ibc,:,:) = q   (ibc,:,:)
            aux (mx+ibc,:,:) = aux (ibc,:,:)
         END DO
!
      CASE (4)
!
!         # solid wall (assumes 2'nd component is velocity or momentum in x):
!
         DO ibc=1,mbc
         DO j = 1-mbc, my+mbc
            q(mx+ibc,j,:) = q(mx+1-ibc,j,:)
         END DO
         END DO
!
!        # negate the normal velocity:
         do ibc=1,mbc
         do j = 1-mbc, my+mbc
            q(mx+ibc,j,2) = -q(mx+1-ibc,j,2)
         END DO
         END DO
!      
      CASE DEFAULT
         STOP 'ILLEGAL VALUE FOR MTHBC(2)'
      END SELECT
!
!
!
      SELECT CASE (mthbc(3)+1) ! # bottom boundary:
!
      CASE (1)   ! # user-specified boundary conditions
!
      CASE (2)   ! # zero-order extrapolation:
!
         DO jbc=1,mbc
         DO i = 1-mbc, mx+mbc
            q   (i,1-jbc,:) = q   (i,1,:)
            aux (i,1-jbc,:) = aux (i,1,:)
         END DO
         END DO
!
      CASE (3)   ! # periodic:  
!
         DO jbc=1,mbc
         DO i = 1-mbc, mx+mbc
            q   (i,1-jbc,:) = q   (i,my+1-jbc,:)
            aux (i,1-jbc,:) = aux (i,my+1-jbc,:)
         END DO
         END DO
!
      CASE (4) !# solid wall (assumes 3rd comp. is vel. or momentum in y):
!
         DO jbc=1,mbc
         DO i = 1-mbc, mx+mbc
            q(i,1-jbc,:) = q(i,jbc,:)
         END DO
         END DO
!
!        # negate the normal velocity:
         DO jbc=1,mbc
         DO i = 1-mbc, mx+mbc
            q(i,1-jbc,3) = -q(i,jbc,3)
         END DO
         END DO
!
      CASE DEFAULT 
         STOP 'ILLEGAL VALUE FOR MTHBC(3)'
      END SELECT
!
!
!
      SELECT CASE (mthbc(4)+1)   ! # top boundary:
!
      CASE (1)    ! # user-specified boundary conditions
!
      CASE (2)    ! # zero-order extrapolation:
!
         DO jbc=1,mbc
         DO i = 1-mbc, mx+mbc
            q   (i,my+jbc,:) = q   (i,my,:)
            aux (i,my+jbc,:) = aux (i,my,:)
         END DO
         END DO
!
      CASE (3)    ! # periodic:  
!
         DO jbc=1,mbc
         DO i = 1-mbc, mx+mbc
            q   (i,my+jbc,:) = q   (i,jbc,:)
            aux (i,my+jbc,:) = aux (i,jbc,:)
         END DO
         END DO
!
      CASE (4) ! # solid wall (assumes 3rd comp. is vel. or momentum in y):
!
         DO jbc=1,mbc
         DO i = 1-mbc, mx+mbc
            q(i,my+jbc,:) = q(i,my+1-jbc,:)
         END DO
         END DO
!
!        # negate the normal velocity:
         DO jbc=1,mbc
         DO i = 1-mbc, mx+mbc
            q(i,my+jbc,3) = -q(i,my+1-jbc,3)
         END DO
         END DO
!
      CASE DEFAULT
        STOP 'ILLEGAL VALUE FOR MTHBC(4)'
      END SELECT
!
      RETURN
      END SUBROUTINE BC2
