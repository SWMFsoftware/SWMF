!
!      PSPHERE.F
!

!      Takes deniu1 (H+) and deniu5(He+) and computes integrated densities 
!      in the equatorial plane for H and He to approximate IMAGE EUV 
!      images (He+) and their H+ counterparts

       include "param_diagB.inc"

!       parameter ( nt = 7 ) 
!       parameter ( nx = 100,
!     .             ny = 200,
!     .             nr =  60,
!     .             nl =  96  )
       parameter ( nr =  60)

!      arrays of interpolated SAMI3 output on an interpolated grid
!      x = lat (deg), y = alt (km)

       real blat0(nx,ny,nl),balt0(nx,ny,nl),blon0(nx,ny,nl)
       real den10B(nx,ny,nl,nt),den20B(nx,ny,nl,nt),den50B(nx,ny,nl,nt)

!      arrays of integrated/projected densities

       real psprad(nr,nl),psplonb(nr,nl)
       real psp1nb(nr,nl,nt),psp2nb(nr,nl,nt),psp5nb(nr,nl,nt)
       real psp1sb(nr,nl,nt),psp2sb(nr,nl,nt),psp5sb(nr,nl,nt)

!      Open interpolated SAMI3 data files

       open(unit=71,file='deni10B.dat' , form='unformatted')
       open(unit=72,file='deni20B.dat' , form='unformatted')
       open(unit=73,file='deni50B.dat' , form='unformatted')

       open(unit=61,file='blon0.dat', form='unformatted')
       open(unit=62,file='blat0.dat', form='unformatted')
       open(unit=63,file='balt0.dat', form='unformatted')

!      open output data files

       open(unit=41,file='psprad0b.dat', form='unformatted')
       open(unit=42,file='psplon0b.dat', form='unformatted')
       open(unit=45,file='psp1n0b.dat', form='unformatted')
       open(unit=46,file='psp2n0b.dat', form='unformatted')
       open(unit=47,file='psp5n0b.dat', form='unformatted')
       open(unit=48,file='psp1s0b.dat', form='unformatted')
       open(unit=49,file='psp2s0b.dat', form='unformatted')
       open(unit=50,file='psp5s0b.dat', form='unformatted')

!      some parameters

       re   = 6370.
       pi   = 4. * atan(1.)
       dtor = pi / 180.

!      read in interpolated data (geo amd mag)

       read(61) blon0
       read(62) blat0
       read(63) balt0
       read(71) den10B
       read(72) den20B
       read(73) den50B
       close(61)
       close(62)
       close(63)
       close(71)
       close(72)
       close(73)

!      Set up grid for plasmasphere
!      We have nr radii in the plasmsphere (equatorial plane) grid

       ix = nx/2
       jy = ny/2
       delr = 8.0*re/float(nr-1)
       do k = 1,nl
         do j = 1,nr
           psplonb(j,k) = blon0(ix,jy,k)
           psprad(j,k)  = delr*float(j-1)
         enddo
       enddo

!      Zero out plasmasphere "images"

       do l = 1,nt
         do k = 1,nl
           do j = 1,nr
             psp1nb(j,k,l) = 0.0
             psp2nb(j,k,l) = 0.0
             psp5nb(j,k,l) = 0.0
             psp1sb(j,k,l) = 0.0
             psp2sb(j,k,l) = 0.0
             psp5sb(j,k,l) = 0.0
           enddo
         enddo
       enddo

!      Fill in plasmsphere images based on ions per cm^2

       dellon = 360.0/float(nl)
       dellatb = blat0(2,1,1) - blat0(1,1,1)
       do l = 1,nt
         print *,' fill in psp1,2,5: l = ',l
         do k = 1,nl
           do j = 1,ny-1

!            Magnetic coordinates             
             do i = 1,nx

               dlon = blon0(i,j,k)
               drad = (balt0(i,j,k) + re)*cos(dtor*blat0(i,j,k))
               dxlat = dtor*dellatb*drad
               dyalt = balt0(i,j+1,k) - balt0(i,j,k)
               dzlon = dtor*dellon*drad
               vol = dxlat*dyalt*dzlon*1.0e15

!              Get number of H+, He+ ions at this grid point
!              We will later convert this to a number per unit area.
               dn1 = den10B(i,j,k,l)*vol
               dn2 = den20B(i,j,k,l)*vol
               dn5 = den50B(i,j,k,l)*vol

!              All of these ions will be assigned to a pair of 
!              grids in the equatorial plane
!              Get radius index, interpolation factors, area
               jrad = int(drad/delr) + 1
               if (jrad .lt. 1) jrad = 1
               if (jrad .gt. nr-2) then
                 jrad = nr-2
                 dn1 = 0.
                 dn5 = 0.
               endif
               fjp = (drad - psprad(jrad,k))/delr
               fj = 1.0 - fjp

!              if ((fjp .gt. 1.0 .or. fjp .le. 0.) .and. i .ne. 101) then
!                  print *, " drad/re, jrad = ",drad/re,jrad 
!                  print *, " dlon, drad = ",dlon, drad
!                  print *, " balt, blat = ",baltu(i,j,k),blatu(i,j,k)
!                  print *, " i,j,k = ", i,j,k
!                  print *, " "
!               endif

               area = dzlon*delr
               dn1 = dn1/(area*1.e10)
               dn2 = dn2/(area*1.e10)
               dn5 = dn5/(area*1.e10)

!               print *," dn1,2,5,vol,area", dn1,dn2,dn5,vol,area 

!              Add this contribution to the north image, 
!              the south image, or both
               if (drad .le. re) then
                 if (blat0(i,j,k) .ge. 0.) then
                   psp1nb(jrad,k,l)   = psp1nb(jrad,k,l)   + fj*dn1
                   psp1nb(jrad+1,k,l) = psp1nb(jrad+1,k,l) + fjp*dn1
                   psp2nb(jrad,k,l)   = psp2nb(jrad,k,l)   + fj*dn2
                   psp2nb(jrad+1,k,l) = psp2nb(jrad+1,k,l) + fjp*dn2
                   psp5nb(jrad,k,l)   = psp5nb(jrad,k,l)   + fj*dn5
                   psp5nb(jrad+1,k,l) = psp5nb(jrad+1,k,l) + fjp*dn5
                 else
                   psp1sb(jrad,k,l)   = psp1sb(jrad,k,l)   + fj*dn1
                   psp1sb(jrad+1,k,l) = psp1sb(jrad+1,k,l) + fjp*dn1
                   psp2sb(jrad,k,l)   = psp2sb(jrad,k,l)   + fj*dn2
                   psp2sb(jrad+1,k,l) = psp2sb(jrad+1,k,l) + fjp*dn2
                   psp5sb(jrad,k,l)   = psp5sb(jrad,k,l)   + fj*dn5
                   psp5sb(jrad+1,k,l) = psp5sb(jrad+1,k,l) + fjp*dn5
                 endif
               else if (drad .gt. re .and. drad .lt. 6.0*re) then
                 psp1nb(jrad,k,l)   = psp1nb(jrad,k,l)   + fj*dn1
                 psp1nb(jrad+1,k,l) = psp1nb(jrad+1,k,l) + fjp*dn1

                 psp2nb(jrad,k,l)   = psp2nb(jrad,k,l)   + fj*dn2
                 psp2nb(jrad+1,k,l) = psp2nb(jrad+1,k,l) + fjp*dn2

                 psp5nb(jrad,k,l)   = psp5nb(jrad,k,l)   + fj*dn5
                 psp5nb(jrad+1,k,l) = psp5nb(jrad+1,k,l) + fjp*dn5

                 psp1sb(jrad,k,l)   = psp1sb(jrad,k,l)   + fj*dn1
                 psp1sb(jrad+1,k,l) = psp1sb(jrad+1,k,l) + fjp*dn1

                 psp2sb(jrad,k,l)   = psp2sb(jrad,k,l)   + fj*dn2
                 psp2sb(jrad+1,k,l) = psp2sb(jrad+1,k,l) + fjp*dn2

                 psp5sb(jrad,k,l)   = psp5sb(jrad,k,l)   + fj*dn5
                 psp5sb(jrad+1,k,l) = psp5sb(jrad+1,k,l) + fjp*dn5
               else
                 fsum = 0.5*fj + 0.25*fjp
                 fsump = 0.25*fj + 0.5*fjp
                 psp1nb(jrad-1,k,l) = psp1nb(jrad-1,k,l) + 0.25*fj*dn1
                 psp1nb(jrad,k,l)   = psp1nb(jrad,k,l)   + fsum*dn1
                 psp1nb(jrad+1,k,l) = psp1nb(jrad+1,k,l) + fsump*dn1
                 psp1nb(jrad+2,k,l) = psp1nb(jrad+2,k,l) + 0.25*fjp*dn1

                 psp2nb(jrad-1,k,l) = psp2nb(jrad-1,k,l) + 0.25*fj*dn2
                 psp2nb(jrad,k,l)   = psp2nb(jrad,k,l)   + fsum*dn2
                 psp2nb(jrad+1,k,l) = psp2nb(jrad+1,k,l) + fsump*dn2
                 psp2nb(jrad+2,k,l) = psp2nb(jrad+2,k,l) + 0.25*fjp*dn2

                 psp5nb(jrad-1,k,l) = psp5nb(jrad-1,k,l) + 0.25*fj*dn5
                 psp5nb(jrad,k,l)   = psp5nb(jrad,k,l)   + fsum*dn5
                 psp5nb(jrad+1,k,l) = psp5nb(jrad+1,k,l) + fsump*dn5
                 psp5nb(jrad+2,k,l) = psp5nb(jrad+2,k,l) + 0.25*fjp*dn5

                 psp1sb(jrad-1,k,l) = psp1sb(jrad-1,k,l) + 0.25*fj*dn1
                 psp1sb(jrad,k,l)   = psp1sb(jrad,k,l)   + fsum*dn1
                 psp1sb(jrad+1,k,l) = psp1sb(jrad+1,k,l) + fsump*dn1
                 psp1sb(jrad+2,k,l) = psp1sb(jrad+2,k,l) + 0.25*fjp*dn1

                 psp2sb(jrad-1,k,l) = psp2sb(jrad-1,k,l) + 0.25*fj*dn2
                 psp2sb(jrad,k,l)   = psp2sb(jrad,k,l)   + fsum*dn2
                 psp2sb(jrad+1,k,l) = psp2sb(jrad+1,k,l) + fsump*dn2
                 psp2sb(jrad+2,k,l) = psp2sb(jrad+2,k,l) + 0.25*fjp*dn2

                 psp5sb(jrad-1,k,l) = psp5sb(jrad-1,k,l) + 0.25*fj*dn5
                 psp5sb(jrad,k,l)   = psp5sb(jrad,k,l)   + fsum*dn5
                 psp5sb(jrad+1,k,l) = psp5sb(jrad+1,k,l) + fsump*dn5
                 psp5sb(jrad+2,k,l) = psp5sb(jrad+2,k,l) + 0.25*fjp*dn5
               endif
             enddo
           enddo
         enddo
       enddo

!      write out interpolated data files

       write(41) psprad
       write(42) psplonb
       write(45) psp1nb
       write(46) psp2nb
       write(47) psp5nb
       write(48) psp1sb
       write(49) psp2sb
       write(50) psp5sb

       close(41)
       close(42)
       close(43)
       close(44)
       close(45)
       close(46)
       close(47)
       close(48)
       close(49)
       close(50)

       stop
       end

       
