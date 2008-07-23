
subroutine read_ring_current

  use ModIonosphere
  use ModHeidiSize
  use ModHeidiCurrents

  implicit none

  real, dimension(1:100,1:100) :: FAC, FAC_tmp
  real, dimension(1:100)       :: lat, mlts, dlat, dmlt

  integer       :: nlats, nmlts, i, j, k, l
  real          :: Re, T, P

  logical       :: done

  Re = 6372.0 * 1000.0

  nmlts = IlFac
  nlats = IrFac

!!  write(*,*) "nmlts : ",nmlts, nlats

  do i=1,nmlts
     mlts(i) = LonFac(i) * IONO_PI / 12.0 
  enddo
  mlts(nmlts+1) = mlts(1) + 2.0 * IONO_PI

!!  write(*,*) "Lon : ",LonFac

  dmlt(1) = mlts(2) - mlts(1)
  do i=2,nmlts
     dmlt(i) = (mlts(i+1) - mlts(i-1))/2.0
  enddo

  do i=1,nlats
     lat(i) = Latfac(i) * 3.141592 / 180.0 
  enddo

  dlat(1) = lat(2) - lat(1)
  do i=2,nlats-1
     dlat(i) = (lat(i+1) - lat(i-1))/2.0
  enddo
  dlat(nlats) = lat(nlats) - lat(nlats-1)

!!  write(*,*) "Lat : ",LatFac

! Turn dlat into a distance. The sin(lat) is only good if we are going to
! be multiplying by the dmlt also, which we are...

  do i=1,nlats
     dlat(i) = dlat(i) * Re * Re * sin(lat(i))
  enddo

  do i=1,nlats
     do j=1,nmlts
        FAC(i,j) = Jfac(i,j) !/ (dlat(i) * dmlt(j))
     enddo
     FAC(i,nmlts+1) = FAC(i,1)
  enddo

!!  write(*,*) "FAC : ",JFac

! increment nmlts to account for adding a wrapping point:

  nmlts = nmlts + 1

  FAC_tmp = FAC

  do j=3,nmlts-2
     i = nlats-1
     done = .false.
     do while (.not.done)
        if (((FAC(i,j).ne.0).and.(FAC(i+1,j).eq.0)).or. &
             ((FAC(i,j).ne.0).and.(FAC(i+1,j-1).eq.0)).or. &
             ((FAC(i,j).ne.0).and.(FAC(i+1,j+1).eq.0))) then
            FAC_tmp(i,j) = 0.0
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

        if ((T < lat(1)).or.(T > lat(nlats))) then
           IONO_NORTH_RCM_JR(i,j) = 0.0
        else 

           k = 1
           do while (T > lat(k))
              k = k + 1
           enddo

           l = 1
           do while (P > mlts(l))
              l = l + 1
           enddo

           IONO_NORTH_RCM_JR(i,j) = FAC(k,l)

           T = IONO_NORTH_Theta(i,j)

!!           write(*,*) "Computing jr : ",IONO_NORTH_RCM_JR(i,j), t

           IONO_NORTH_RCM_JR(i,j) = -2.0 * IONO_NORTH_RCM_JR(i,j) *       &
                cos(T) / sqrt(1.0 + 3.0 * cos(T) ** 2)

        endif

     enddo

  enddo

!!  write(*,*) "Done with north"

  do i = 1, IONO_nTheta

     do j = 1, IONO_nPsi

        T = IONO_SOUTH_Theta(i,j) - IONO_PI/2

        P = mod(IONO_SOUTH_Psi(i,j) + IONO_PI, IONO_PI*2)

        if ((T < lat(1)).or.(T > lat(nlats))) then
           IONO_SOUTH_RCM_JR(i,j) = 0.0
        else 

           k = 1
           do while (T > lat(k))
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

  write(*,*) "Done with read_ring_current", &
       minval(IONO_NORTH_RCM_JR), maxval(IONO_NORTH_RCM_JR)

  write(*,*) "slice : ", IONO_NORTH_RCM_JR(:,97)

  return

end subroutine read_ring_current

