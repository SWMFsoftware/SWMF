!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine read_ring_current(iteration)

  implicit none

  use ModIonosphere
  use ModHeidiSize
  use ModHeidiCurrents

  integer, intent(in)          :: iteration

  real, dimension(1:100,1:100) :: FAC, FAC_tmp
  real, dimension(1:100)       :: lats, mlts, dlat, dmlt, dummy

  character*100 :: filein, input, format
  integer       :: ierr, inunit, nlats, nmlts, i, j, k, l
  real          :: Re, T, P

  logical       :: done

  Re = 6372.0 * 1000.0

  do i=1,LonFac
     mlts(i) = LonFac(i) * IONO_PI / 12.0 
  enddo
  mlts(LonFac+1) = mlts(1) + 2.0 * IONO_PI

  dmlt(1) = mlts(2) - mlts(1)
  do i=2,IlFac
     dmlt(i) = (mlts(i+1) - mlts(i-1))/2.0
  enddo

  read(inunit,'(a)') input
  read(inunit,format) (lats(i),i=1,nlats)

  do i=1,Ifac
     lats(i) = lats(i) * 3.141592 / 180.0 
  enddo

  dlat(1) = lats(2) - lats(1)
  do i=2,nlats-1
     dlat(i) = (lats(i+1) - lats(i-1))/2.0
  enddo
  dlat(nlats) = lats(nlats) - lats(nlats-1)

! Turn dlat into a distance. The sin(lat) is only good if we are going to
! be multiplying by the dmlt also, which we are...

  do i=1,nlats
     dlat(i) = dlat(i) * Re * Re * sin(lats(i))
  enddo

  read(inunit,'(a)') input

  write(format, '(A1,I2,A6)') "(",IlFac,"E11.3)"

! this is the time string

  do l=0,iteration

     read(inunit,'(a)') input

     if (l == iteration) write(6,*) input

     ! the / 2.0 is for seperating into northern and southern hemispheres

     do i=1,nlats
        read(inunit,format) (FAC(i,j),j=1,IlFac)
        do j=1,IlFac
           FAC(i,j) = FAC(i,j) / (dlat(i) * dmlt(j)) / 2.0
!           if (i >= nlats-2) FAC(i,j) = 0.0
        enddo
        FAC(i,IlFac+1) = FAC(i,1)
     enddo

     read(inunit,format) (dummy(j),j=1,IlFac)

  enddo

  FAC_tmp = FAC

  do j=3,IlFac-1
     i = nlats-1
     done = .false.
     do while (.not.done)
        if (((FAC(i,j).ne.0).and.(FAC(i+1,j).eq.0)).or. &
             ((FAC(i,j).ne.0).and.(FAC(i+1,j-1).eq.0)).or. &
             ((FAC(i,j).ne.0).and.(FAC(i+1,j+1).eq.0))) then
            FAC_tmp(i,j) = 0.0
!           done = .true.
        else
           FAC_tmp(i,j) = FAC_tmp(i,j)
        endif
        i = i - 1
        if (i.eq.1) done = .true.
     enddo
  enddo

  FAC = FAC_tmp

  do i = 1, IONO_nTheta

     do j = 1, IONO_nPsi

        T = IONO_PI/2.0 - IONO_NORTH_Theta(i,j)
        P = mod(IONO_NORTH_Psi(i,j) + IONO_PI, IONO_PI*2)

        if ((T < lats(1)).or.(T > lats(nlats))) then
           IONO_NORTH_RCM_JR(i,j) = 0.0
        else 

           k = 1
           do while (T > lats(k))
              k = k + 1
           enddo

           l = 1
           do while (P > mlts(l))
              l = l + 1
           enddo

           IONO_NORTH_RCM_JR(i,j) = FAC(k,l)

           T = IONO_NORTH_Theta(i,j)

           IONO_NORTH_RCM_JR(i,j) = -2.0 * IONO_NORTH_RCM_JR(i,j) *       &
                cos(T) / sqrt(1.0 + 3.0 * cos(T) ** 2)

        endif

     enddo

  enddo

  do i = 1, IONO_nTheta

     do j = 1, IONO_nPsi

        T = IONO_SOUTH_Theta(i,j) - IONO_PI/2
        P = mod(IONO_SOUTH_Psi(i,j) + IONO_PI, IONO_PI*2)

        if ((T < lats(1)).or.(T > lats(nlats))) then
           IONO_SOUTH_RCM_JR(i,j) = 0.0
        else 

           k = 1
           do while (T > lats(k))
              k = k + 1
           enddo

           l = 1
           do while (P > mlts(l))
              l = l + 1
           enddo

           if (l > 1) l = l - 1

           IONO_SOUTH_RCM_JR(i,j) = -1.0*FAC(k,l)

        endif

        T = IONO_SOUTH_Theta(i,j)

        IONO_SOUTH_RCM_JR(i,j) = -2.0 * IONO_SOUTH_RCM_JR(i,j) *       &
             cos(T) / sqrt(1.0 + 3.0 * cos(T) ** 2)

     enddo

  enddo

  return

end subroutine read_ring_current

