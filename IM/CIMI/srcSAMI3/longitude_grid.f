! LONGITUDE_GRID.F

       parameter ( ntotal = 96 )

       real*8 theta(ntotal),theta_smooth(ntotal),tmp(ntotal)
       real*8 dtheta_smooth(ntotal),theta_smooth_shift(ntotal)
       real*8 dtheta_smooth_shift(ntotal)
       real*8 theta_min,theta_max,dtheta
       real*8 dtheta_left,dtheta_right,dtheta_big,dtheta_small

       logical uniform_grid

       uniform_grid = .true.
       

       if ( uniform_grid ) then
       
         theta_min = 0.
         theta_max = 360.
         dtheta    = theta_max/float(ntotal)
                
         do i = 1,ntotal
           theta(i) = float(i-1) * dtheta
         enddo

         open(20,file='longitude.inp',form='unformatted')
         write(20) theta
         close(20)

       else

         theta_min     = 168.
         theta_max     = 192.
         dtheta_small  = 0.1

         nfine  = int ( ( theta_max - theta_min ) / dtheta_small ) + 1
         nlarge = ntotal - nfine

         dtheta_left   = theta_min
         dtheta_right  = ( 360. - theta_max )
         dtheta_big    = ( dtheta_left + dtheta_right ) /
     .                    float(nlarge)

         nleft         = (dtheta_left / dtheta_big) 
         nright        = nleft + nfine

         do i = 1,nleft
           theta(i) = float(i-1) * dtheta_big
         enddo      
         
         do i = nleft+1,nleft+nfine+1
           theta(i) = theta_min + float(i-(nleft+1)) * dtheta_small
         enddo

         do i = nleft+nfine+2,ntotal
           theta(i) = theta_max + float(i-(nleft+nfine+1)) * dtheta_big
         enddo

         do i = 1,ntotal
           theta_smooth(i) = theta(i)
           tmp(i)          = theta(i)
         enddo

         msmooth = 124

         do j = 1,msmooth
           do i = 3,ntotal-4 
             im1 = i - 1
             ip1 = i + 1
             theta_smooth(i) = 0.25 * ( tmp(im1)    +
     .                                  2. * tmp(i) +
     .                                  tmp(ip1)      )
           enddo
           do i = 1,ntotal
             tmp(i) = theta_smooth(i)
           enddo
         enddo

         do i = 2,ntotal-1
           print *,i,theta(i),0.5*(theta(i+1)-theta(i-1))
         enddo

         do i = 2,ntotal-1
           print *,i,theta_smooth(i),
     .               0.5*(theta_smooth(i+1)-theta_smooth(i-1))
           dtheta_smooth(i) = 0.5*(theta_smooth(i+1)-theta_smooth(i-1))
         enddo

         open(20,file='theta_smooth.inp', form='unformatted')
         open(21,file='dtheta_smooth.inp',form='unformatted')

         write(20) theta_smooth
         write(21) dtheta_smooth

         close(20)
         close(21)

         istart = ntotal - 23
         idiff  = ntotal - istart

         do i = 1,idiff
           i0 = i + istart
           theta_smooth_shift(i) = theta_smooth(i0) -
     .                             theta_smooth(istart+1)
         enddo       

         do i = idiff+1,ntotal
           i0 = i - idiff
           theta_smooth_shift(i) = theta_smooth(i0) +
     .                             theta_smooth(idiff+1)
         enddo

         do i = 2,ntotal-1
           dtheta_smooth_shift(i) = 0.5*( theta_smooth_shift(i+1)
     .                                   - theta_smooth_shift(i-1))
         enddo

         open(20,file='theta_smooth_shift.inp', form='unformatted')
         open(21,file='dtheta_smooth_shift.inp',form='unformatted')

         write(20) theta_smooth_shift
         write(21) dtheta_smooth_shift

         close(20)
         close(21)

       endif

       end

