
!
!      INTERP-3D-U1P.F
!

!      takes u1pu in sami3 coordinates (zaltu,glatu,glonu):(nz,nf,nl,nt)
!      and interpolates to a regular grid (zalt0,glat0,glon0):(nx,ny,nl,nt)
!        - the grid in the z-direction can be nonuniform (gamy)

       include "param_diag.inc"

!      arrays of 'real' data

       real glatu(nz,nf,nl),zaltu(nz,nf,nl),glonu(nz,nf,nl)
       real u1pu(nz,nf,nl,nt),u1ptmp(nz,nf,nl)
       real xu(nz,nf,nl),zu(nz,nf,nl)
       real zualtmax(nf,nl)    !gj


!      arrays of interpolated data

       real x0(nx,ny,nl),z0(nx,ny,nl)
       integer i0(nx,ny,nl),j0(nx,ny,nl)

       real glat0i(nx,ny,nl),glon0i(nx,ny,nl),zalt0i(nx,ny,nl)
       real u1p0i(nx,ny,nl,nt)

       real glat0(nx,ny,nl),glon0(nx,ny,nl),zalt0(nx,ny,nl)
       real u1p0(nx,ny,nl,nt)

!      open 'real' data files

       open(unit= 9,file='glonu.dat', form='unformatted')
       open(unit=10,file='glatu.dat', form='unformatted')
       open(unit=11,file='zaltu.dat', form='unformatted')
       open(unit=12,file='u1pu.dat' , form='unformatted')

!      open interpolated data files

       open(unit=13,file='glat0i.dat', form='unformatted')
       open(unit=23,file='glon0i.dat', form='unformatted')
       open(unit=14,file='zalt0i.dat', form='unformatted')
       open(unit=15,file='u1p0i.dat', form='unformatted')

       open(unit=31,file='glat0.dat', form='unformatted')
       open(unit=32,file='glon0.dat', form='unformatted')
       open(unit=33,file='zalt0.dat', form='unformatted')
       open(unit=34,file='u1p0.dat', form='unformatted')

!      some parameters

       re   = 6370.
       pi   = 4. * atan(1.)
       dtor = pi / 180.
       gamy = 3.

!      read in 'real' data

       read( 9) glonu
       read(10) glatu
       read(11) zaltu

       do l = 1,nt
         print *,' reading data: l = ',l
         read(12) u1ptmp
         do k = 1,nl
           do j = 1,nf
             do i = 1,nz
               u1pu(i,j,k,l) = u1ptmp(i,j,k)
             enddo
           enddo
         enddo
       enddo

!      get min and max of glatu and zaltu

       print *,'getting min/max glatu/zaltu - setting grid0'

       do k = 1,nl

         glatmin =   90.
         glatmax =  -90.
         zaltmin =   85.
         zaltmax =    0.

!      zaltmin is maximum altitude of lowest field line 
 
         j = 1 
         do i = 1,nz 
           if ( zaltu(i,j,k) .ge. zaltmin ) zaltmin = zaltu(i,j,k) 
         enddo 

         do j = 1,nf
           do i = 1,nz
             if ( glatu(i,j,k) .le. glatmin ) glatmin = glatu(i,j,k)
             if ( glatu(i,j,k) .ge. glatmax ) glatmax = glatu(i,j,k)
             if ( zaltu(i,j,k) .ge. zaltmax ) zaltmax = zaltu(i,j,k)
           enddo
           zualtmax(j,k) = zaltmax
         enddo

!        set up grid for glat0i (uniform)

         delg0  = ( glatmax - glatmin ) / float(nx-1)

         do iy = 1,ny 
           do ix = 1,nx
             glat0i(ix,iy,k) = glatmin + ( ix - 1 ) * delg0
             glat0i(ix,iy,k) = -90. + (ix - 1) * 180. / (nx - 1)
           enddo
         enddo

!        set up grid for zalt0i (non-uniform)

         zaltmax0 = 2000.

         do iy = 1,ny
           do ix = 1,nx
             iy0  = ny + 1 - iy
             dy   = float(iy-1) / float(ny-1)
             y2   = dy * sinh ( gamy )
             y3   = alog ( y2 + sqrt ( y2**2 + 1. ) ) / gamy
             zalt0i(ix,iy0,k) = zaltmin + ( zaltmax0 - zaltmin ) 
     .                                           * ( 1. - y3 )
           enddo
         enddo

       enddo


!      obtain xu and zu 
!      (cartesian coordinates of glatu and zaltu)

       print *,'setting up xu/zu'

       do k = 1,nl
         do j = 1,nf
           do i = 1,nz
             xu(i,j,k) = ( zaltu(i,j,k) + re ) * 
     .                 cos ( glatu(i,j,k) * dtor )
             zu(i,j,k) = ( zaltu(i,j,k) + re ) * 
     .                 sin ( glatu(i,j,k) * dtor )
           enddo
         enddo
       enddo

!        obtain x0 and z0 
!        (cartesian coordinates of glat0i and zalt0i, i.e.,
!        the interpolated grid)

       print *,'setting up glat0i/zalt0i'

       do k = 1,nl
         do j = 1,ny
           do i = 1,nx
             x0(i,j,k) = ( zalt0i(i,j,k) + re ) 
     .                   * cos ( glat0i(i,j,k) * dtor )
             z0(i,j,k) = ( zalt0i(i,j,k) + re ) 
     .                    * sin ( glat0i(i,j,k) * dtor )
           enddo
         enddo
       enddo

!      determine interior point 
!      using cross-products from opposing vertices
!      (determines where interpolated grid data is
!      in terms of the 'read' grid data)

       print *,'getting interior points'

       do k = 1,nl
         print *,'k = ',k
         do iy = 1,ny
! limit the number of field lines searched.
! if zu is greater than any point on the "real" field line, we don't 
!  need to search it  -- gj
             jtmp = 1
!  Remember zalt0i is constant in ix
             do while(zalt0i(1,iy,k) .gt. zualtmax(jtmp,k) .and. 
     .            jtmp .le. nf-1)
                jtmp = jtmp +1
             enddo
             j00 = max(jtmp-1,1)
!             print *, ' iy k j00 ',iy,k,j00
!             print *,  ' jtmp zalt0i zualtmax ',jtmp,zalt0i(1,iy,k)
!     .               ,zualtmax(jtmp,k) 
             istart = 1
           do ix = 1,nx
             i0(ix,iy,k) = 0
             j0(ix,iy,k) = 0
             ifind   = 0

             do j = j00,nf-1   !-- gj

               if ( ifind .eq. 0 ) then

                 do i = istart,nz-1
                   if ( ifind .eq. 0 ) then
                     ax  = xu(i+1,j,k) - xu(i,j,k)
                     ay  = zu(i+1,j,k) - zu(i,j,k)
                     bx  = xu(i,j+1,k) - xu(i,j,k)
                     by  = zu(i,j+1,k) - zu(i,j,k)
                     cx  = x0(ix,iy,k)   - xu(i,j,k)
                     cy  = z0(ix,iy,k)   - zu(i,j,k)
                     axc = ax * cy - ay * cx
                     bxc = bx * cy - by * cx
                     if ( axc * bxc .le. 0 .and. 
     .                    axc .le. 0       .and.
     .                    ifind .eq. 0           ) then
                       ax  = xu(i+1,j,k) - xu(i+1,j+1,k)
                       ay  = zu(i+1,j,k) - zu(i+1,j+1,k)
                       bx  = xu(i,j+1,k) - xu(i+1,j+1,k)
                       by  = zu(i,j+1,k) - zu(i+1,j+1,k)
                       cx  = x0(ix,iy,k)   - xu(i+1,j+1,k)
                       cy  = z0(ix,iy,k)   - zu(i+1,j+1,k)
                       axc = ax * cy - ay * cx
                       bxc = bx * cy - by * cx
                       if ( axc * bxc .le. 0 .and.
     .                      axc .ge. 0             ) then
                         i0(ix,iy,k) = i
                         j0(ix,iy,k) = j
                         ifind   = 1
                         istart = max(1,i-2)
!                         istart  = 1
                       endif
                     endif
                   endif
                 enddo
               endif
             enddo
           enddo
         enddo
       enddo     

!      set default value of interpolated data
!      (in this case the electron density)

       print *,'initializing u1p0i'

       do n = 1,nt
         do k = 1,nl
           do j = 1,ny
             do i = 1,nx
                u1p0i(i,j,k,n) = 1.
             enddo
           enddo
         enddo
       enddo

!      lay down interpolated data on interpolated grid
!      using an area weighted lay down scheme

       print *,'setting u1p0i'

       do n = 1,nt
         do k = 1,nl
           do iy = 1,ny
             do ix = 1,nx
                i           = i0(ix,iy,k)
                j           = j0(ix,iy,k)
                if ( i .ne. 0 .and. j .ne. 0 ) then
                  xp             = x0(ix,iy,k)
                  yp             = z0(ix,iy,k)
                  call area_subxz(i,j,k,xu,zu,xp,yp,as1,as2,as3,as4)
                  a_tot = as1 + as2 + as3 + as4
                  u1p0i(ix,iy,k,n) = 
     .            ( as1 * u1pu(i,j,k,n)   + as2 * u1pu(i+1,j,k,n) +
     .              as4 * u1pu(i,j+1,k,n) + as3 * u1pu(i+1,j+1,k,n) ) 
     .            / a_tot
                endif
             enddo
           enddo
         enddo
       enddo

!      initialize glon0i 

       do k = 1,nl 
         do j = 1,ny
           do i = 1,nx
             glon0i(i,j,k) = 0.
           enddo
         enddo
       enddo

!      calculate interpolated value of glon

       do k = 1,nl
         do iy = 1,ny
           do ix = 1,nx
              i              = i0(ix,iy,k)
              j              = j0(ix,iy,k)
              if ( i .ne. 0 .and. j .ne. 0 ) then
                xp             = x0(ix,iy,k)
                yp             = z0(ix,iy,k)
                call area_subxz(i,j,k,xu,zu,xp,yp,as1,as2,as3,as4)
                a_tot = as1 + as2 + as3 + as4
                glon1 = glonu(i,j,k)
                glon2 = glonu(i+1,j,k)
                glon3 = glonu(i,j+1,k)
                glon4 = glonu(i+1,j+1,k)
                glont = glon1 + glon2 + glon3 + glon4
                if ( (glon1 .gt. 355. .or.
     .                glon2 .gt. 355. .or.
     .                glon3 .gt. 355. .or.
     .                glon4 .gt. 355.     ) .and.
     .                glont .lt. 1200.            ) then
                  glon2 = glon1
                  glon3 = glon1
                  glon4 = glon1
                endif
                glon0i(ix,iy,k) = 
     .          (   as1 * glon1
     .            + as2 * glon2
     .            + as4 * glon3
     ,            + as3 * glon4  )
     .          / a_tot 
              endif
           enddo
         enddo
       enddo

!      define longitude 'outside' of data

       j00 = 1
       do k = 1,nl 
         do i = 1,nx
           do j = 1,ny
             if ( glon0i(i,j,k) .ne. 0 ) j00 = j
             if ( glon0i(i,j,k) .eq. 0 ) 
     .            glon0i(i,j,k) = glon0i(i,j00,k)
           enddo
         enddo
       enddo

!      fix end points

       do k = 1,nl
         do j = 1,ny
           glon0i(1,j,k)  = glon0i(2,1,k)
           glon0i(nx,j,k) = glon0i(nx-1,1,k)
         enddo
       enddo

!      write out interpolated data files
  
       write(13) glat0i
       write(23) glon0i
       write(14) zalt0i
       write(15) u1p0i

!      set longitude 'straight' at longitude at
!      peak of the magnetic flux tube
!      set glat0 to be glat0i
!      set zalt0 to be zalt0i
!      set u1p0 to be u1p0i

       do k = 1,nl
         do j = 1,ny
           do i = 1,nx
             glon0(i,j,k) = glon0i(nx/2,ny,k)
             glat0(i,j,k) = glat0i(i,j,k)
             zalt0(i,j,k) = zalt0i(i,j,k)
           enddo
         enddo
         print *,'k,glon0',
     .            k,glon0(nx/2,ny,k)
       enddo

       do n = 1,nt
         do k = 1,nl
           do j = 1,ny
             do i = 1,nx
               u1p0(i,j,k,n) = u1p0i(i,j,k,n)
             enddo
           enddo
         enddo
       enddo

       do n = 1,nt
         do k = 2,nl-1
           do j = 2,ny-1
             do i = 2,nx-1
               if ( u1p0(i,j,k,n) .eq. 1. ) 
     .            u1p0(i,j,k,n) = u1p0(i,j,k-1,n)
c$$$     .                           ( te0i(i-1,j-1,k-1,n) +
c$$$     .                             te0i(i+1,j-1,k-1,n) +
c$$$     .                             te0i(i-1,j+1,k-1,n) +
c$$$     .                             te0i(i+1,j+1,k-1,n) +
c$$$     .                             te0i(i-1,j-1,k+1,n) +
c$$$     .                             te0i(i+1,j-1,k+1,n) +
c$$$     .                             te0i(i-1,j+1,k+1,n) +
c$$$     .                             te0i(i+1,j+1,k+1,n)  ) / 8.
             enddo
           enddo
         enddo
       enddo


       write(31) glat0
       write(32) glon0
       write(33) zalt0
       write(34) u1p0


       stop
       end

       
*******************************************
*******************************************

!            area_subxz

*******************************************
*******************************************


        subroutine area_subxz(i,j,k,x,z,x0,z0,as1,as2,as3,as4)

!       calculate areas of cell sides
!       break each quadrilateral side into
!       two triangles and use the formula: 
!           A = (1/2)|a x b|
!       where
!           a: vector from A to B
!           b: vector from A to C

        include "param_diag.inc"

        real x(nz,nf,nl),z(nz,nf,nl)

!       as1

        ax1 = 0.5 * ( x(i+1,j,k) + x(i+1,j+1,k) ) - x0
        az1 = 0.5 * ( z(i+1,j,k) + z(i+1,j+1,k) ) - z0

        bx1 = x(i+1,j+1,k) - x0
        bz1 = z(i+1,j+1,k) - z0

        cz1 =    ax1 * bz1 - az1 * bx1

        a1 = 0.5 * sqrt ( cz1 * cz1 )

        ax2 = 0.5 * ( x(i,j+1,k) + x(i+1,j+1,k) ) - x0
        az2 = 0.5 * ( z(i,j+1,k) + z(i+1,j+1,k) ) - z0

        bx2 = x(i+1,j+1,k) - x0
        bz2 = z(i+1,j+1,k) - z0

        cz2 =    ax2 * bz2 - az2 * bx2

        a2 = 0.5 * sqrt ( cz2*cz2 )

        as1 =  a1 + a2  

!       as2

        ax1 = 0.5 * ( x(i,j+1,k) + x(i+1,j+1,k) ) - x0
        az1 = 0.5 * ( z(i,j+1,k) + z(i+1,j+1,k) ) - z0

        bx1 = x(i,j+1,k) - x0
        bz1 = z(i,j+1,k) - z0

        cz1 =    ax1 * bz1 - az1 * bx1

        a1 = 0.5 * sqrt ( cz1 * cz1 )

        ax2 = 0.5 * ( x(i,j,k) + x(i,j+1,k) ) - x0
        az2 = 0.5 * ( z(i,j,k) + z(i,j+1,k) ) - z0

        bx2 = x(i,j+1,k) - x0
        bz2 = z(i,j+1,k) - z0

        cz2 =    ax2 * bz2 - az2 * bx2

        a2 = 0.5 * sqrt ( cz2*cz2 )

        as2 =  a1 + a2  

!       as3

        ax1 = 0.5 * ( x(i,j,k) + x(i,j+1,k) ) - x0
        az1 = 0.5 * ( z(i,j,k) + z(i,j+1,k) ) - z0

        bx1 = x(i,j,k) - x0
        bz1 = z(i,j,k) - z0

        cz1 =    ax1 * bz1 - az1 * bx1

        a1 = 0.5 * sqrt ( cz1 * cz1 )

        ax2 = 0.5 * ( x(i,j,k) + x(i+1,j,k) ) - x0
        az2 = 0.5 * ( z(i,j,k) + z(i+1,j,k) ) - z0

        bx2 = x(i,j,k) - x0
        bz2 = z(i,j,k) - z0

        cz2 =    ax2 * bz2 - az2 * bx2

        a2 = 0.5 * sqrt ( cz2*cz2 )

        as3 =  a1 + a2  

!       as4

        ax1 = 0.5 * ( x(i,j,k) + x(i+1,j,k) ) - x0
        az1 = 0.5 * ( z(i,j,k) + z(i+1,j,k) ) - z0

        bx1 = x(i+1,j,k) - x0
        bz1 = z(i+1,j,k) - z0

        cz1 =    ax1 * bz1 - az1 * bx1

        a1 = 0.5 * sqrt ( cz1 * cz1 )

        ax2 = 0.5 * ( x(i+1,j+1,k) + x(i+1,j,k) ) - x0
        az2 = 0.5 * ( z(i+1,j+1,k) + z(i+1,j,k) ) - z0

        bx2 = x(i+1,j,k) - x0
        bz2 = z(i+1,j,k) - z0

        cz2 =    ax2 * bz2 - az2 * bx2

        a2 = 0.5 * sqrt ( cz2*cz2 )

        as4 =  a1 + a2  

        return
        end
