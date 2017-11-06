
          parameter ( nextra = 4, nseg = 6, nze = nseg * nextra )
          parameter ( nz0 = 200, nz = nz0+nze)

          real qpnew(nz0+1)
          real qpextra(nze+nseg+1)
          real qp(nz+1)

          do i = 1,nz0+1
            qpnew(i) = i - 101.
            print *,i,qpnew(i)
          enddo

          nz0h   = nz0/2
          iextra = 0
          do ie = 1,nseg
            do iee = 1,nextra+1
              iextra = iextra + 1
!              print *,'ie,iee,iextra',ie,iee,iextra
               delq   = qpnew(nz0h+(1-nseg/2)+ie) - 
     .                  qpnew(nz0h+(1-nseg/2)+(ie-1))
             qpextra(iextra) = qpnew(nz0h+(1-nseg/2)+(ie-1)) 
     .                         + delq/float(nextra+1) * (iee-1)
!             print *,'and ...',delq,qpnew(nz0h-2+(ie-1)),qpextra(iextra)
            enddo
          enddo

          qpextra(nze+nseg+1) = qpnew(nz0h+(1+nseg/2))

          do i = 1,nze+nseg+1
            print *,'i,qpextra',i,qpextra(i)
          enddo

          do i = 1,nz0h+(1-nseg/2)-1
            qp(i) = qpnew(i)
!            print *,'i1',i,qp(i)
          enddo

          do i = nz0h+(1-nseg/2),nz0h+(1-nseg/2)+(nze+nseg+1)-1
            qp(i) = qpextra(i+1-(nz0h+(1-nseg/2)))
!            print *,'i2',i,qp(i)
          enddo

          do i = nz0h+(1-nseg/2)+(nze+nseg+1),nz+1
            ii    = i + 1 - (nz0h+(1-nseg/2)+(nze+nseg+1))
     .             + nz0h + (1+nseg/2)
            qp(i) = qpnew(ii)
!            print *,'i3',i,ii,qpnew(ii)
          enddo
!               stop
          do i = 1,nz+1
            print *,'i,qp',i,qp(i)
          enddo

          end

