
subroutine read_amie(iter)

  use ModIonosphere
  implicit none

  integer, intent(in)          :: iter

  real, dimension(0:100,0:100) :: FAC, tmp, cond_eflx, cond_avee
  real, dimension(1:100)       :: lats, mlts, dlat, dmlt, dummy

  character*100 :: filein, input, format
  character*30  :: field
  integer       :: ierr, inunit, nlats, nmlts, i, j, k, l, ntime,ntimes, n
  integer       :: iyear,imon,ida,ihr,imin, nfields, ibinary_unit
  real          :: Re, T, P, dx, dy

  real          :: SWV,BXGSM, BYGSM,BZGSM,AENDX
  real          :: AECALC,AUCALC,ALCALC,DSTRD,DSTNDX
  real          :: HEMPWR,SMPLJL,potdif

  Re = 6372.0 * 1000.0

  nmlts = 25
  nlats = 23

  inunit = 73

  filein = "b910604-07.mo"

  open(unit = inunit, file = filein, iostat = ierr, status = 'old', &
       form='UNFORMATTED')

  if (ierr.ne.0) then 
     write(6,*) "An error occured while trying to read ",filein
     write(6,*) "in subroutine read_amie"
     stop
  endif

  read(inunit) nlats, nmlts, ntimes
  read(inunit) (lats(i),i=1,nlats)
  read(inunit) (mlts(i),i=1,nmlts)

  do i=1,nmlts
     mlts(i) = mlts(i) * IONO_PI / 12.0 
  enddo

  dmlt(1) = mlts(2) - mlts(1)
  do i=2,nmlts-1
     dmlt(i) = (mlts(i+1) - mlts(i-1))/2.0
  enddo
  dmlt(nmlts) = mlts(nmlts) - mlts(nmlts-1)

  do i=1,nlats
     lats(i) = lats(i) * 3.141592 / 180.0 
  enddo

  dlat(1) = lats(2) - lats(1)
  do i=2,nlats-1
     dlat(i) = (lats(i+1) - lats(i-1))/2.0
  enddo
  dlat(nlats) = lats(nlats) - lats(nlats-1)

  read(inunit) nfields

  do i=1,nfields
     read(inunit) field
  enddo

  do n=0,iter*12

     read(inunit) ntime,iyear,imon,ida,ihr,imin

     read(inunit) SWV,BXGSM, BYGSM,BZGSM,AENDX,     &
          AECALC,AUCALC,ALCALC,DSTRD,DSTNDX,        &
          HEMPWR,SMPLJL,potdif

! potential
     read(inunit) tmp(1:nmlts,1:nlats)

! pedersen conductance
     read(inunit) tmp(1:nmlts,1:nlats)

! pedersen conductance - aurora only
     read(inunit) tmp(1:nmlts,1:nlats)

! hall conductance
     read(inunit) tmp(1:nmlts,1:nlats)

! hall conductance - aurora only
     read(inunit) tmp(1:nmlts,1:nlats)

! Auroral Mean Energy (keV)
     read(inunit) cond_avee(1:nmlts,1:nlats)

! Auroral Energy Flux (W/m2)
     read(inunit) cond_eflx(1:nmlts,1:nlats)

! Electric Field (East) 
     read(inunit) tmp(1:nmlts,1:nlats)

! Electric Field (North)
     read(inunit) tmp(1:nmlts,1:nlats)

! Field Aligned Current
     read(inunit) FAC(1:nmlts,1:nlats)

! other stuff...

     do i=1,7
        read(inunit) tmp(1:nmlts,1:nlats)
     enddo

  enddo

  write(6,*) "reading time : ",ntime, iyear, imon, ida, ihr, imin

  close(inunit)

  do i = 1, IONO_nTheta

     do j = 1, IONO_nPsi

        T = IONO_NORTH_Theta(i,j)
        P = mod(IONO_NORTH_Psi(i,j) + IONO_PI, IONO_PI*2)

        if ((T < lats(1)).or.(T > lats(nlats))) then
           IONO_NORTH_AMIE_JR(i,j) = 0.0
        else 

           k = 1
           do while (T > lats(k))
              k = k + 1
           enddo

           l = 1
           do while (P > mlts(l))
              l = l + 1
           enddo

           dx = 1.0 - (lats(k)-T)/dlat(k)
           dy = 1.0 - (mlts(l)-P)/dmlt(l)

           IONO_NORTH_AMIE_JR(i,j)  =(    dx)*(    dy)*FAC(l  ,k  ) + &
                                     (1.0-dx)*(1.0-dy)*FAC(l-1,k-1) + &
                                     (    dx)*(1.0-dy)*FAC(l-1,k  ) + &
                                     (1.0-dx)*(    dy)*FAC(l  ,k-1)

           IONO_NORTH_EFlux(i,j) = (    dx)*(    dy)*cond_eflx(l  ,k  ) + &
                                   (1.0-dx)*(1.0-dy)*cond_eflx(l-1,k-1) + &
                                   (    dx)*(1.0-dy)*cond_eflx(l-1,k  ) + &
                                   (1.0-dx)*(    dy)*cond_eflx(l  ,k-1)

           IONO_NORTH_Ave_E(i,j) = (    dx)*(    dy)*cond_avee(l  ,k  ) + &
                                   (1.0-dx)*(1.0-dy)*cond_avee(l-1,k-1) + &
                                   (    dx)*(1.0-dy)*cond_avee(l-1,k  ) + &
                                   (1.0-dx)*(    dy)*cond_avee(l  ,k-1)

           T = IONO_NORTH_Theta(i,j)

           IONO_NORTH_AMIE_JR(i,j) = -2.0 * IONO_NORTH_AMIE_JR(i,j) *       &
                cos(T) / sqrt(1.0 + 3.0 * cos(T) ** 2)

        endif

     enddo

  enddo

  do i = 1, IONO_nTheta

     do j = 1, IONO_nPsi

        T = IONO_PI - IONO_SOUTH_Theta(i,j) 
        P = mod(IONO_SOUTH_Psi(i,j) + IONO_PI, IONO_PI*2)

        if ((T < lats(1)).or.(T > lats(nlats))) then
           IONO_SOUTH_AMIE_JR(i,j) = 0.0
        else 

           k = 1
           do while (T > lats(k))
              k = k + 1
           enddo

           l = 1
           do while (P > mlts(l))
              l = l + 1
           enddo

           dx = 1.0 - (lats(k)-T)/dlat(k)
           dy = 1.0 - (mlts(l)-P)/dmlt(l)

           IONO_SOUTH_AMIE_JR(i,j)  =(    dx)*(    dy)*FAC(l  ,k  ) + &
                                     (1.0-dx)*(1.0-dy)*FAC(l-1,k-1) + &
                                     (    dx)*(1.0-dy)*FAC(l-1,k  ) + &
                                     (1.0-dx)*(    dy)*FAC(l  ,k-1)

           IONO_SOUTH_EFlux(i,j) = (    dx)*(    dy)*cond_eflx(l  ,k  ) + &
                                   (1.0-dx)*(1.0-dy)*cond_eflx(l-1,k-1) + &
                                   (    dx)*(1.0-dy)*cond_eflx(l-1,k  ) + &
                                   (1.0-dx)*(    dy)*cond_eflx(l  ,k-1)

           IONO_SOUTH_Ave_E(i,j) = (    dx)*(    dy)*cond_avee(l  ,k  ) + &
                                   (1.0-dx)*(1.0-dy)*cond_avee(l-1,k-1) + &
                                   (    dx)*(1.0-dy)*cond_avee(l-1,k  ) + &
                                   (1.0-dx)*(    dy)*cond_avee(l  ,k-1)

        endif

        T = IONO_SOUTH_Theta(i,j)

        IONO_SOUTH_AMIE_JR(i,j) = -2.0 * IONO_SOUTH_AMIE_JR(i,j) *       &
             cos(T) / sqrt(1.0 + 3.0 * cos(T) ** 2)

     enddo

  enddo

  return

end subroutine read_amie

