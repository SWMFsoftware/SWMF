c      =======================================================
       subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,t,dt)
c      =======================================================
c
       implicit double precision (a-h,o-z)
       dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn),
     &            aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
c
c      # dummy subroutine for use when equation has no source term.
c      # If method(5)=0 then this routine is never called, but its
c      # existence may be required by some compilers.
c
       DO i = 1-mbc, mx+mbc
       DO j = 1-mbc, my+mbc
       DO m = 1, meqn
          q(i,j,m) = q(i,j,m)*EXP(-aux(i,j,3)*dt)
       END DO
       END DO
       END DO
c
       RETURN
       END SUBROUTINE Src2
