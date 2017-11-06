! RW_SAMI3.F

       parameter ( nz = 160,
     .             nf = 200,
     .             nl =  90,
     .             nt = 80  )

       parameter ( re = 6370. )

       real blat(nz,nf,nl),
     .      blon(nz,nf,nl),
     .      balt(nz,nf,nl)

       real dene(nz,nf,nl,nt)

       real x(nz,nf,nl),
     .      y(nz,nf,nl),
     .      z(nz,nf,nl)

       real dene0(nz,nf,nl)

       open(101, file='blatu.dat', form='unformatted')
       open(102, file='blonu.dat', form='unformatted')
       open(103, file='baltu.dat', form='unformatted')
       open(104, file='deneu.dat', form='unformatted')
       open(106, file='dene_dx.dat', form='unformatted')

       iskip = 1
       pi    = 4. * atan(1.)
       dtor  = pi / 180.

       read(101) blat
       read(102) blon
       read(103) balt
    
       print *,'reading dene'
       do l = 1,nt
         print *,'l = ',l
         read(104) dene0
         do k = 1,nl
           do j = 1,nf
             do i = 1,nz
               dene(i,j,k,l) = dene0(i,j,k)
             enddo
           enddo
         enddo
       enddo

       do k = 1,nl
         do j = 1,nf
           do i = 1,nz
             r        = balt(i,j,k)
             x(i,j,k) = r * cos(blat(i,j,k)*dtor) 
     .                    * cos(blon(i,j,k)*dtor)
             y(i,j,k) = r * cos(blat(i,j,k)*dtor) 
     .                    * sin(blon(i,j,k)*dtor)
             z(i,j,k) = r * sin(blat(i,j,k)*dtor)
           enddo
         enddo
       enddo

       j0 = 1

       do l = 1,nt,iskip
         print *,'j0',j0
         do k = 1,nl
           do j = 1,nf
             do i = 1,nz
               dene0(i,j,k) = dene(i,j,k,l)
             enddo
           enddo
         enddo
         write(106) x,y,z,dene0
         j0 = j0 + 1
       enddo

       stop
       end

       
