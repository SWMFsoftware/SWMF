
!
!      INTERP-3D-DENE-DB.F
!

!      takes deneu in sami3 coordinates (baltu,blatu,blonu):(nz,nf,nl,nt)
!      and interpolates to a regular grid (balt0,blat0,blon0):(nx,ny,nl,nt)
!        - the grid in the z-direction can be nonuniform (gamy)

       include "param_diagB.inc"

!      arrays of 'real' data

       real blatu(nz,nf,nl),baltu(nz,nf,nl),blonu(nz,nf,nl)
       real deneu(nz,nf,nl,nt),denetmp(nz,nf,nl)
       real den1u(nz,nf,nl,nt),den1tmp(nz,nf,nl)
       real den2u(nz,nf,nl,nt),den2tmp(nz,nf,nl)
       real den5u(nz,nf,nl,nt),den5tmp(nz,nf,nl)
       real vs1u(nz,nf,nl,nt),vs1tmp(nz,nf,nl)
       real vs2u(nz,nf,nl,nt),vs2tmp(nz,nf,nl)
       real vs5u(nz,nf,nl,nt),vs5tmp(nz,nf,nl)
       real xu(nz,nf,nl),zu(nz,nf,nl)
       real zualtmax(nf,nl)    !gj
       real time(5,nt),tm(nt)

!      arrays of interpolated data

       real x0(nx,ny,nl),z0(nx,ny,nl)
       integer i0(nx,ny,nl),j0(nx,ny,nl)

       real blat0i(nx,ny,nl),blon0i(nx,ny,nl),balt0i(nx,ny,nl)
       real dene0i(nx,ny,nl,nt)
       real den10i(nx,ny,nl,nt),den20i(nx,ny,nl,nt),den50i(nx,ny,nl,nt)
       real vs10i(nx,ny,nl,nt),vs20i(nx,ny,nl,nt),vs50i(nx,ny,nl,nt)

       real blat0(nx,ny,nl),blon0(nx,ny,nl),balt0(nx,ny,nl)
       real dene0(nx,ny,nl,nt)
       real den10(nx,ny,nl,nt),den20(nx,ny,nl,nt),den50(nx,ny,nl,nt)
       real den50flash(nx,ny,nl,nt)

       real dene0derot(nx,ny,nl,nt),den50derot(nx,ny,nl,nt)
       real den10derot(nx,ny,nl,nt),den20derot(nx,ny,nl,nt)
       real den50derotflash(nx,ny,nl,nt)

       real vs10(nx,ny,nl,nt),vs20(nx,ny,nl,nt),vs50(nx,ny,nl,nt)

!      open 'real' data files

       open(unit= 8,file='time.dat', form='formatted')
       open(unit= 9,file='blonu.dat', form='unformatted')
       open(unit=10,file='blatu.dat', form='unformatted')
       open(unit=11,file='baltu.dat', form='unformatted')
       open(unit=12,file='deneu.dat' , form='unformatted')
!      1 H+, 2 O+, 5 He+
       open(unit=41,file='deniu1.dat' , form='unformatted')
       open(unit=42,file='deniu2.dat' , form='unformatted')
       open(unit=45,file='deniu5.dat' , form='unformatted')
       open(unit=51,file='vsiu1.dat' , form='unformatted')
       open(unit=52,file='vsiu2.dat' , form='unformatted')
       open(unit=55,file='vsiu5.dat' , form='unformatted')

!      open interpolated data files

       open(unit=13,file='blat0i.dat', form='unformatted')
       open(unit=23,file='blon0i.dat', form='unformatted')
       open(unit=14,file='balt0i.dat', form='unformatted')
       open(unit=15,file='dene0Bi.dat', form='unformatted')

       open(unit=31,file='blat0.dat', form='unformatted')
       open(unit=32,file='blon0.dat', form='unformatted')
       open(unit=33,file='balt0.dat', form='unformatted')
       open(unit=34,file='dene0B.dat', form='unformatted')
       open(unit=35,file='deni10B.dat', form='unformatted')
       open(unit=36,file='deni20B.dat', form='unformatted')
       open(unit=37,file='deni50B.dat', form='unformatted')
       open(unit=38,file='vsi10B.dat', form='unformatted')
       open(unit=39,file='vsi20B.dat', form='unformatted')
       open(unit=30,file='vsi50B.dat', form='unformatted')

       open(unit=64,file='dene0B-derot.dat', form='unformatted')
       open(unit=65,file='deni10B-derot.dat', form='unformatted')
       open(unit=66,file='deni20B-derot.dat', form='unformatted')
       open(unit=67,file='deni50B-derot.dat', form='unformatted')

       open(unit=77,file='deni50B-derot-flash.dat', form='unformatted')
       open(unit=78,file='deni50B-flash.dat', form='unformatted')

!      some parameters

       re   = 6370.
       pi   = 4. * atan(1.)
       dtor = pi / 180.
       gamy = 3.

!      read in 'real' data

       read( 8,*) time
       read( 9) blonu
       read(10) blatu
       read(11) baltu

       tm = time(2,:) + time(3,:)/60.0 + time(4,:)/3600.0

       do k = 1,nl
         do j = 1,nf
           do i = 1,nz
             baltu(i,j,k) = baltu(i,j,k) - re
           enddo
         enddo
       enddo

       do l = 1,nt
         print *,' reading dene,i1,i2,i5,vi1,vi2,vi5: l = ',l
         read(12) denetmp
         read(41) den1tmp
         read(42) den2tmp
         read(45) den5tmp
         read(51) vs1tmp
         read(52) vs2tmp
         read(55) vs5tmp
         do k = 1,nl
           do j = 1,nf
             do i = 1,nz
               deneu(i,j,k,l) = denetmp(i,j,k)
               den1u(i,j,k,l) = den1tmp(i,j,k)
               den2u(i,j,k,l) = den2tmp(i,j,k)
               den5u(i,j,k,l) = den5tmp(i,j,k)
               vs1u(i,j,k,l) = vs1tmp(i,j,k)
               vs2u(i,j,k,l) = vs2tmp(i,j,k)
               vs5u(i,j,k,l) = vs5tmp(i,j,k)
             enddo
           enddo
         enddo
       enddo

!      get min and max of blatu and baltu

       print *,'getting min/max blatu/baltu - setting grid0'

       do k = 1,nl

         blatmin =   90.
         blatmax =  -90.
         baltmin =   85.
         baltmax =    0.

!      baltmin is maximum altitude of lowest field line 
 
         j = 1 
         do i = 1,nz 
           if ( baltu(i,j,k) .ge. baltmin ) baltmin = baltu(i,j,k) 
         enddo 

         do j = 1,nf
           do i = 1,nz
             if ( blatu(i,j,k) .le. blatmin ) blatmin = blatu(i,j,k)
             if ( blatu(i,j,k) .ge. blatmax ) blatmax = blatu(i,j,k)
             if ( baltu(i,j,k) .ge. baltmax ) baltmax = baltu(i,j,k)
           enddo
             zualtmax(j,k) = baltmax
         enddo

!        set up grid for blat0i (uniform)

         delg0  = ( blatmax - blatmin ) / float(nx-1)

         do iy = 1,ny 
           do ix = 1,nx
             blat0i(ix,iy,k) = blatmin + ( ix - 1 ) * delg0
           enddo
         enddo

!        set up grid for balt0i (non-uniform)

         baltmax0 = 40000.
 
         do iy = 1,ny
           do ix = 1,nx
             iy0  = ny + 1 - iy
             dy   = float(iy-1) / float(ny-1)
             y2   = dy * sinh ( gamy )
             y3   = alog ( y2 + sqrt ( y2**2 + 1. ) ) / gamy
             balt0i(ix,iy0,k) = baltmin + ( baltmax0 - baltmin ) 
     .                                           * ( 1. - y3 )
           enddo
         enddo
       enddo


!      obtain xu and zu 
!      (cartesian coordinates of blatu and baltu)

       print *,'setting up xu/zu'

       do k = 1,nl
         do j = 1,nf
           do i = 1,nz
             xu(i,j,k) = ( baltu(i,j,k) + re ) * 
     .                 cos ( blatu(i,j,k) * dtor )
             zu(i,j,k) = ( baltu(i,j,k) + re ) * 
     .                 sin ( blatu(i,j,k) * dtor )
           enddo
         enddo
       enddo

!        obtain x0 and z0 
!        (cartesian coordinates of blat0i and balt0i, i.e.,
!        the interpolated grid)

       print *,'setting up blat0i/balt0i'

       do k = 1,nl
         do j = 1,ny
           do i = 1,nx
             x0(i,j,k) = ( balt0i(i,j,k) + re ) 
     .                   * cos ( blat0i(i,j,k) * dtor )
             z0(i,j,k) = ( balt0i(i,j,k) + re ) 
     .                    * sin ( blat0i(i,j,k) * dtor )
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
!  Remember balt0i is constant in ix
             do while(balt0i(1,iy,k) .gt. zualtmax(jtmp,k) .and. 
     .            jtmp .le. nf-1)
                jtmp = jtmp +1
             enddo
             j00 = max(jtmp-1,1)
!             print *, ' iy k j00 ',iy,k,j00
!             print *,  ' jtmp balt0i zualtmax ',jtmp,balt0i(1,iy,k)
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

       print *,'initializing dene0i,den1/2/50i,vs1/2/50i'

       do n = 1,nt
         do k = 1,nl
           do j = 1,ny
             do i = 1,nx
                dene0i(i,j,k,n) = 1.
                den10i(i,j,k,n) = 1.
                den20i(i,j,k,n) = 1.
                den50i(i,j,k,n) = 1.
                vs10i(i,j,k,n) = 0.01
                vs20i(i,j,k,n) = 0.01
                vs50i(i,j,k,n) = 0.01
             enddo
           enddo
         enddo
       enddo

!      lay down interpolated data on interpolated grid
!      using an area weighted lay down scheme

       print *,'setting dene0i,den10i,den20i,den50i,vs10i,vs20i,vs50i'

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
                  dene0i(ix,iy,k,n) = 
     .            ( as1 * deneu(i,j,k,n)   + as2 * deneu(i+1,j,k,n) +
     .              as4 * deneu(i,j+1,k,n) + as3 * deneu(i+1,j+1,k,n) ) 
     .            / a_tot
                   den10i(ix,iy,k,n) = 
     .            ( as1 * den1u(i,j,k,n)   + as2 * den1u(i+1,j,k,n) +
     .              as4 * den1u(i,j+1,k,n) + as3 * den1u(i+1,j+1,k,n) ) 
     .            / a_tot
                   den20i(ix,iy,k,n) = 
     .            ( as1 * den2u(i,j,k,n)   + as2 * den2u(i+1,j,k,n) +
     .              as4 * den2u(i,j+1,k,n) + as3 * den2u(i+1,j+1,k,n) ) 
     .            / a_tot
                  den50i(ix,iy,k,n) = 
     .            ( as1 * den5u(i,j,k,n)   + as2 * den5u(i+1,j,k,n) +
     .              as4 * den5u(i,j+1,k,n) + as3 * den5u(i+1,j+1,k,n) ) 
     .            / a_tot
                   vs10i(ix,iy,k,n) = 
     .            ( as1 * vs1u(i,j,k,n)   + as2 * vs1u(i+1,j,k,n) +
     .              as4 * vs1u(i,j+1,k,n) + as3 * vs1u(i+1,j+1,k,n) ) 
     .            / a_tot
                   vs20i(ix,iy,k,n) = 
     .            ( as1 * vs2u(i,j,k,n)   + as2 * vs2u(i+1,j,k,n) +
     .              as4 * vs2u(i,j+1,k,n) + as3 * vs2u(i+1,j+1,k,n) ) 
     .            / a_tot
                   vs50i(ix,iy,k,n) = 
     .            ( as1 * vs5u(i,j,k,n)   + as2 * vs5u(i+1,j,k,n) +
     .              as4 * vs5u(i,j+1,k,n) + as3 * vs5u(i+1,j+1,k,n) ) 
     .            / a_tot
                endif
             enddo
           enddo
         enddo
       enddo

!      initialize blon0i 

       do k = 1,nl 
         do j = 1,ny
           do i = 1,nx
             blon0i(i,j,k) = 0.
           enddo
         enddo
       enddo

!      calculate interpolated value of blon

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
                blon1 = blonu(i,j,k)
                blon2 = blonu(i+1,j,k)
                blon3 = blonu(i,j+1,k)
                blon4 = blonu(i+1,j+1,k)
                blont = blon1 + blon2 + blon3 + blon4
                if ( (blon1 .gt. 355. .or.
     .                blon2 .gt. 355. .or.
     .                blon3 .gt. 355. .or.
     .                blon4 .gt. 355.     ) .and.
     .                blont .lt. 1200.            ) then
                  blon2 = blon1
                  blon3 = blon1
                  blon4 = blon1
                endif
                blon0i(ix,iy,k) = 
     .          (   as1 * blon1
     .            + as2 * blon2
     .            + as4 * blon3
     ,            + as3 * blon4  )
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
             if ( blon0i(i,j,k) .ne. 0 ) j00 = j
             if ( blon0i(i,j,k) .eq. 0 ) 
     .            blon0i(i,j,k) = blon0i(i,j00,k)
           enddo
         enddo
       enddo

!      fix end points

       do k = 1,nl
         do j = 1,ny
           blon0i(1,j,k)  = blon0i(2,1,k)
           blon0i(nx,j,k) = blon0i(nx-1,1,k)
         enddo
       enddo

!      write out interpolated data files
  
       write(13) blat0i
       write(23) blon0i
       write(14) balt0i
       write(15) dene0i

!      set longitude 'straight' at longitude at
!      peak of the magnetic flux tube
!      set blat0 to be blat0i
!      set balt0 to be balt0i
!      set dene0 to be dene0i
!      set den10 to be den10i
!      set den20 to be den20i, etc

       do k = 1,nl
         do j = 1,ny
           do i = 1,nx
             blon0(i,j,k) = blon0i(nx/2,ny,k)
             blat0(i,j,k) = blat0i(i,j,k)
             balt0(i,j,k) = balt0i(i,j,k)
           enddo
         enddo
         print *,'k,blon0',
     .            k,blon0(nx/2,ny,k)
       enddo

       do n = 1,nt
         kshift = nint(tm(n)*nl/24.0)
         do k = 1,nl
           do j = 1,ny
             do i = 1,nx
               dene0(i,j,k,n) = dene0i(i,j,k,n)
               den10(i,j,k,n) = den10i(i,j,k,n)
               den20(i,j,k,n) = den20i(i,j,k,n)
               den50(i,j,k,n) = den50i(i,j,k,n)
               vs10(i,j,k,n) = vs10i(i,j,k,n)
               vs20(i,j,k,n) = vs20i(i,j,k,n)
               vs50(i,j,k,n) = vs50i(i,j,k,n)
!              These are de-corotated for plotting
               kde = k - kshift
               if (kde .lt. 1) kde = kde + nl
               dene0derot(i,j,kde,n) = dene0i(i,j,k,n)
               den10derot(i,j,kde,n) = den10i(i,j,k,n)
               den20derot(i,j,kde,n) = den20i(i,j,k,n)
               den50derot(i,j,kde,n) = den50i(i,j,k,n)
!              These have extra density at local time noon for diagnostics
               den50derotflash(i,j,kde,n) = den50i(i,j,k,n)
               den50flash(i,j,k,n) = den50i(i,j,k,n)
               if (k .eq. 49) then 
                 den50derotflash(i,j,kde,n)= 2.0*den50i(i,j,k,n)
                 den50flash(i,j,k,n)= 2.0*den50i(i,j,k,n)
               endif
             enddo
           enddo
         enddo
       enddo

       write(31) blat0
       write(32) blon0
       write(33) balt0
       write(34) dene0
       write(35) den10
       write(36) den20
       write(37) den50
       write(38) vs10
       write(39) vs20
       write(30) vs50

       write(64) dene0derot
       write(65) den10derot
       write(66) den20derot
       write(67) den50derot
       write(77) den50derotflash
       write(78) den50flash

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

        include "param_diagB.inc"

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
