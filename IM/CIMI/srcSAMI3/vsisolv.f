*******************************************
*******************************************

!            vsisolv

*******************************************
*******************************************

      subroutine vsisolv ( vii,vid,viold,snuj,nfl,nll,cs )

      include 'param3_mpi-1.98.inc'
      include 'com3_mpi-1.98.inc'

      dimension a(nz), b(nz), c(nz), d(nz)
      real vii(nz),vid(nz),viold(nz),snuj(nz),cs(nz)

! initialize

      do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
      enddo

      do j = 2,nz-1

        ujm1 = vii(j-1)
        uj   = vii(j)
        ujp1 = vii(j+1)
        ur = .25*( uj +ujp1)
        ul = .25*( uj +ujm1)

        if (ur .ge. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 =  ur
          c0 =  0.
        endif
        if (ur .le. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = -ul
          c0 = ur
        endif
        if (ur .ge. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = ur - ul
          c0 = 0.
        endif
        if (ur .le. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 = 0.
          c0 = ur
        endif

! anomalous drag to prevent ions from going supersonic

          delcs        = 0.1 * cs(j)
          alpha_drag_p = ( vii(j) - 0.9 * cs(j) )  / delcs
          alpha_drag_n = ( vii(j) + 0.9 * cs(j) )  / delcs
          anu_drag    = 0.5 * anu_drag0 * ( 1. + tanh(alpha_drag_p)) + 
     .                  0.5 * anu_drag0 * ( 1. - tanh(alpha_drag_n)) 

         
        a(j) = a0 / d22s(j,nfl,nll) * bms(j,nfl,nll)
        b(j) = 1/dt + snuj(j) + b0 / d22s(j,nfl,nll) *
     .                          bms(j,nfl,nll) + anu_drag
        c(j) = c0 / d22s(j,nfl,nll) * bms(j,nfl,nll)
        d(j) = viold(j)/dt + vid(j)

      enddo

! we will assume that the bc's are zero
! at both ends of the field line

! lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = 0.

! upper bc
 
      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = 0.

      call rtryds(a,b,c,d,vii,nz)

      return
      end
