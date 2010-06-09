subroutine write_ring_current

  use ModIonoHeidi
  use ModHeidiSize
  use ModHeidiCurrents
  use ModHeidiDGCPM
  use ModNumConst, ONLY: cPi

  implicit none

  real, dimension(1:100,1:100) :: FAC
  real, dimension(1:100)       :: lat, mlts
  real          :: T, P, kfac, lfac, potfacn, potfacs, potfac
  real          :: potmax, potmin
  integer       :: nlats, nmlts, i, j, k, l, kk, l1
  integer       :: kkmax, lmax, kmin, lmin
  !----------------------------------------------------------------------------

  nmlts = IlFac
  nlats = IrFac

  do i=1,nmlts
     mlts(i) = LonFac(i) * cPi / 12.0 
  enddo
  mlts(nmlts+1) = mlts(1) + 2.0 * cPi

  do i=1,nlats
     lat(i) = Latfac(i) * cPi / 180.0 
  enddo

  kkmax=1
  lmax=1
  kmin=1
  lmin=1
  potmax=0.
  potmin=0.

  do i = 1, nlats

     do j = 1, nmlts

        T = cPi/2.0 - lat(i)
        P = mod(mlts(j) + cPi, cPi*2)

        k = 1
        do while (T > IONO_NORTH_THETA(k,1))
           k = k + 1
        enddo
	kk=k-1
	If (k.eq.1) Then
           kk=k
           kfac=0.
	else
           kfac=(IONO_NORTH_THETA(k,1)-T)/(IONO_NORTH_THETA(k,1)-IONO_NORTH_THETA(k-1,1))
	endif

        l = 1
        do while (P > IONO_NORTH_PSI(1,l))
           l = l + 1
        enddo
	l1=l-1
	If (l.eq.1) Then
           l1=l
           lfac=0.
	else
           lfac=(IONO_NORTH_PSI(1,l)-P)/(IONO_NORTH_PSI(1,l)-IONO_NORTH_PSI(1,l-1))
	endif

!!!...Use the next 4 lines to average N and S hemispheres
 !        FPOT(i,j) = (1.-kfac)*(1.-lfac)*(IONO_NORTH_PHI(k,l) + IONO_SOUTH_PHI(IONO_nTheta-k,l))/2.0
 !        FPOT(i,j) = FPOT(i,j) + kfac*(1.-lfac)*(IONO_NORTH_PHI(kk,l) + IONO_SOUTH_PHI(IONO_nTheta-kk,l))/2.0
 !        FPOT(i,j) = FPOT(i,j) + (1.-kfac)*lfac*(IONO_NORTH_PHI(k,l1) + IONO_SOUTH_PHI(IONO_nTheta-k,l1))/2.0
 !        FPOT(i,j) = FPOT(i,j) + kfac*lfac*(IONO_NORTH_PHI(kk,l1) + IONO_SOUTH_PHI(IONO_nTheta-kk,l1))/2.0

!!!...OR

!!!...Use the next lines to choose the smaller potential of the 2
	potfacn=(1.-kfac)*(1.-lfac)*IONO_NORTH_PHI(k,l)
	potfacn=potfacn+kfac*(1.-lfac)*IONO_NORTH_PHI(kk,l)
	potfacn=potfacn+(1.-kfac)*lfac*IONO_NORTH_PHI(k,l1)
	potfacn=potfacn+kfac*lfac*IONO_NORTH_PHI(kk,l1)
!!!	potfacs=(1.-kfac)*(1.-lfac)*IONO_SOUTH_PHI(IONO_nTheta-k,l)
!!!	potfacs=potfacs+kfac*(1.-lfac)*IONO_SOUTH_PHI(IONO_nTheta-kk,l)
!!!	potfacs=potfacs+(1.-kfac)*lfac*IONO_SOUTH_PHI(IONO_nTheta-k,l1)
!!!	potfacs=potfacs+kfac*lfac*IONO_SOUTH_PHI(IONO_nTheta-kk,l1)
!!!	If (abs(potfacn).le.abs(potfacs)) then
        FPOT(i,j) = potfacn
!!!	else
!!!           FPOT(i,j) = potfacs
!!!	endif

	if (FPOT(i,j).gt.potmax) then
           potmax=FPOT(i,j)
           kkmax=i
           lmax=j
	else if (FPOT(i,j).lt.potmin) then
           potmin=FPOT(i,j)
           kmin=i
           lmin=j
	endif

     enddo

  enddo

  print *, 'FPOT max: ',potmax,kkmax,lmax,1./cos(lat(kkmax))**2,mlts(lmax)*12./cPi
  print *, 'FPOT min: ',potmin,kmin,lmin,1./cos(lat(kmin))**2,mlts(lmin)*12./cPi

!!!...Now fill in the potentials on the DGCPM grid

  kkmax=1
  lmax=1
  kmin=1
  lmin=1
  potmax=0.
  potmin=0.

  do i = 1, nthetacells

     do j = 1, nphicells

        T = vthetacells(i)*cPi/180.  !vthetacells is colatitude
        P = mod(vphicells(j)*cPi/180. + cPi, cPi*2)

        k = 1
        do while (T > IONO_NORTH_THETA(k,1))
           k = k + 1
        enddo
	kk=k-1
	If (k.eq.1) Then
           kk=k
           kfac=0.
	else
           kfac=(IONO_NORTH_THETA(k,1)-T)/(IONO_NORTH_THETA(k,1)-IONO_NORTH_THETA(k-1,1))
	endif

        l = 1
        do while (P > IONO_NORTH_PSI(1,l))
           l = l + 1
        enddo
	l1=l-1
	If (l.eq.1) Then
           l1=l
           lfac=0.
	else
           lfac=(IONO_NORTH_PSI(1,l)-P)/(IONO_NORTH_PSI(1,l)-IONO_NORTH_PSI(1,l-1))
	endif

!!!...Use the next 4 lines to average N and S hemispheres
 !        potdgcpm(i,j) = potdgcpm(i,j) + (1.-kfac)*(1.-lfac)*(IONO_NORTH_PHI(k,l) + IONO_SOUTH_PHI(IONO_nTheta-k,l))/2.0
 !        potdgcpm(i,j) = potdgcpm(i,j) + kfac*(1.-lfac)*(IONO_NORTH_PHI(kk,l) + IONO_SOUTH_PHI(IONO_nTheta-kk,l))/2.0
 !        potdgcpm(i,j) = potdgcpm(i,j) + (1.-kfac)*lfac*(IONO_NORTH_PHI(k,l1) + IONO_SOUTH_PHI(IONO_nTheta-k,l1))/2.0
 !        potdgcpm(i,j) = potdgcpm(i,j) + kfac*lfac*(IONO_NORTH_PHI(kk,l1) + IONO_SOUTH_PHI(IONO_nTheta-kk,l1))/2.0

!!!...OR
	potfacn=(1.-kfac)*(1.-lfac)*IONO_NORTH_PHI(k,l)
	potfacn=potfacn+kfac*(1.-lfac)*IONO_NORTH_PHI(kk,l)
	potfacn=potfacn+(1.-kfac)*lfac*IONO_NORTH_PHI(k,l1)
	potfacn=potfacn+kfac*lfac*IONO_NORTH_PHI(kk,l1)
 !	potfacs=(1.-kfac)*(1.-lfac)*IONO_SOUTH_PHI(IONO_nTheta-k,l)
 !	potfacs=potfacs+kfac*(1.-lfac)*IONO_SOUTH_PHI(IONO_nTheta-kk,l)
 !	potfacs=potfacs+(1.-kfac)*lfac*IONO_SOUTH_PHI(IONO_nTheta-k,l1)
 !	potfacs=potfacs+kfac*lfac*IONO_SOUTH_PHI(IONO_nTheta-kk,l1)
 !	If (abs(potfacn).le.abs(potfacs)) then
        potfac = potfacn
        !	else
        !           potfac = potfacs
        !	endif

	potdgcpm(i,j) = potdgcpm(i,j) + potfac

	if (potfac.gt.potmax) then
           potmax=potfac
           kkmax=i
           lmax=j
	else if (potfac.lt.potmin) then
           potmin=potfac
           kmin=i
           lmin=j
	endif

     enddo

  enddo

  print *, 'potfac max: ',potmax,kkmax,lmax,vlzcells(kkmax),vphicells(lmax)/15.
  print *, 'potfac min: ',potmin,kmin,lmin,vlzcells(kmin),vphicells(lmin)/15.

  return

end subroutine write_ring_current

