! File name: heidi_currents.f90
!
! Contains: pressure and current calculation routines for HEIDI
!	PRESSURES
!	CURRENTSETUP
!	CURRENTCALC
!	Funcpf
!	Funcpc
!
! Last Modified: March 2006, Mike Liemohn
!
!=======================================================================
!			 	PRESSURES
!     Routine calculates the pressures from the distribution function
!=======================================================================
subroutine PRESSURES

  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiWaves
  use ModHeidiCurrents

  implicit none

  integer :: I,J,K,L,I_1,I_2
  real    :: RFAC,SUME,SUMA,SUMN,SUMTE,SUMTN,J_FAC
  !---------------------------------------------------------------------   
  
!open (unit = 3, file = "PPAR_S.dat")
!   write (3,*)'Numerical values '
!    write (3,*)'PPAR PPER'


  


  !...Start the main species loop
  do S=1,NS
     if (SCALC(S).eq.1) then

        !   RFAC is 4*PI*SQRT(2)*(kg->amu * keV->J)^1.5 *(m->km)^6 *(cm->m)^3
        !   and should be ~5.3E-7/M(amu)^1.5

        RFAC=4.*PI*sqrt(2.)*1.E-24*(Q*1.E3/MP/M1(S))**1.5

        !   J_FAC is (keV->J)*(T->G)*(m->cm)^3 /RE(in m)
        !   and should be ~2.5E-13
	J_FAC=Q*1.E13/RE
	SUMTN=0.
	SUMTE=0.
	do J=1,JO
           do I=2,ILMP(J)
              PPER(I,J,S)=0.
              PPAR(I,J,S)=0.
              RNHT(I,J,S)=0.
              EDEN(I,J,S)=0.
              Nspace(I,J,S)=0.
              Espace(I,J,S)=0.
              do K=2,KO
                 F2(I,J,K,1,S)=F2(I,J,K,2,S)
                 SUME=0.
                 SUMA=0.
                 SUMN=0.
                 do L=1,LO
                    SUME=SUME+F2(I,J,K,L,S)*EPME(I,j,K,L)
                    SUMA=SUMA+F2(I,J,K,L,S)*EPMA(I,j,K,L)
                    SUMN=SUMN+F2(I,J,K,L,S)*ERNM(I,j,K,L)

                   ! write(3,*) l,ERNM(i,j,k,l), EPMA(i,j,k,l), EPME(i,j,k,l) 
                    

                    Nspace(I,J,S)=Nspace(I,J,S)+F2(I,J,K,L,S)*WE(K)*DR*DPHI*WMU(L)*CONSL(K,S)
                    Espace(I,J,S)=Espace(I,J,S)+F2(I,J,K,L,S)*WE(K)*DR*DPHI*WMU(L)*EKEV(K)*CONSL(K,S)
                 enddo  ! L loop
                 
                 PPER(I,J,S)=PPER(I,J,S)+EPP(K,S)*SUME
                 PPAR(I,J,S)=PPAR(I,J,S)+EPP(K,S)*SUMA
                 RNHT(I,J,S)=RNHT(I,J,S)+ERNH(K,S)*SUMN
                 EDEN(I,J,S)=EDEN(I,J,S)+ERNH(K,S)*EKEV(K)*SUMN
                 
                 !write(3,*) 'SUME,SUMA',SUME,SUMA
                 
                 !write(3,*) 'PPER, PPAR',PPER(i,j,1), PPAR(i,j,1)
                 
                 !write(3,*) 'SUME,SUMA,EPP(K,S)*SUME,EPP(K,S)*SUMA',&
                 !SUME,SUMA,EPP(K,S)*SUME,EPP(K,S)*SUMA
                 
                 !write(3,*) 'i,j,s,PPER, PPAR',i,j,s,PPER(i,j,s), PPAR(i,j,s)


              enddo   ! K loop
              !...These parameters are equatorial plane values (not bounce-integrated)
              ANIS(I,J,S)=PPER(I,J,S)/2./PPAR(I,J,S)-1.
              EPAR(I,J,S)=2*PPAR(I,J,S)/RNHT(I,J,S)  ! kT parallel [keV] 
              RNHT(I,J,S)=RNHT(I,J,S)*RFAC	! RC dens [1/cm3]
              EDEN(I,J,S)=EDEN(I,J,S)*RFAC	! RC ener dens [keV/cm3]
              PPAR(I,J,S)=2*RFAC*PPAR(I,J,S) ! P parallel [keV/cm3]
              PPER(I,J,S)=RFAC*PPER(I,J,S)	! P perpend [keV/cm3]
              !...These parameters are globally integrated values
              SUMTN=SUMTN+Nspace(I,J,S)
              SUMTE=SUMTE+Espace(I,J,S)
           end do   ! I loop
           do I=ILMP(J)+1,IO
              ANIS(I,J,S)=0.
              EPAR(I,J,S)=0.
              RNHT(I,J,S)=0.
              EDEN(I,J,S)=0.
              PPAR(I,J,S)=0.
              PPER(I,J,S)=0.
              JPER(I,J,S)=0.
           end do   ! I loop
           do I=2,ILMP(J)
              I_1=I+1
              I_2=I
              if (I.eq.ILMP(J)) then
                 I_1=ILMP(J)
                 I_2=ILMP(J)-1
              end if
              JPER(I,J,S)=-J_FAC*((PPER(I_1,J,S)-PPER(I_2,J,S))/DL1+   &
                   3.*(PPER(I,J,S)-PPAR(I,J,S))/LZ(I))/BE(I,1) !J perp [A/m2]
           end do   ! I loop
	end do    ! J loop
        !stop
        !...These values are bounce-integrated
        !	NTOT(S)=2.E6*RFAC*RE*RE*SUMTN	! Ntotal for RC species [parts]
        !	ETOT(S)=2.E6*RFAC*RE*RE*SUMTE	! Etotal for RC species [keV]
	NTOT(S)=SUMTN
	ETOT(S)=SUMTE
	Dst(S)=-3.98E-30*ETOT(S)	! Dst* [nT] 
        !	print *, 'Dst: ',S,Dst(S),ETOT(S),NTOT(S)

     end if    ! SCALC check
  end do    ! S loop
end subroutine PRESSURES
!=======================================================================
!			 	CURRENTSETUP
!     Routine sets up the arrays for the current calculation
!=======================================================================
subroutine CURRENTSETUP

  use ModHeidiSize
  use ModHeidiMain
  use ModHeidiCurrents

  implicit none

  real    :: Lambda(Slen),R22(NR+3,Slen),Rionos,Lmax,dlam,lam1,lam2
  integer :: kk,i,j,k
  !---------------------------------------------------------------------   
  !  Setup calculations

  Ir=IO+3   ! add 1 to the bottom end, 2 to the top end
  Lsh(2:IO+1)=LZ(1:IO)
  Lsh(1)=Lsh(2)-DL1
  Lsh(IO+2)=Lsh(IO+1)+DL1
  Lsh(IO+3)=Lsh(IO+2)+DL1
  sp(1:JO)=sin(phi(1:JO))
  cp(1:JO)=cos(phi(1:JO))
  Lmax=70.
  Ko2=15
  dlam=Lmax/real(Ko2)
  do k=1,Ko2 
     Lambda(k)=dlam*real(k-1)
  enddo
  do k=1,Ir 
     Kmax(k) = Ko2
  enddo
  Rionos=1.02    ! in RE
  Lats(1:Ir)=acos(sqrt(Rionos/Lsh(1:Ir)))
  cl(1:Ir)=cos(Lats(1:Ir))
  sl(1:Ir)=sin(Lats(1:Ir))
  Lats(1:Ir)=Lats(1:Ir)*180.0/pi
  do i=1,Ir
     do k=Ko2,1,-1
        if (Lambda(k).gt.Lats(i)) Kmax(i)=k-1
     end do
  end do
  do i=1,Ir
     rl(i,1:Ko2)=Lambda(1:Ko2)*pi/180.
     rl(i,Kmax(i))=Lats(i)*pi/180.
  end do
  do i=1,Ir
     drl(i,2:Ko2-1)=.5*(rl(i,3:Ko2)-rl(i,1:Ko2-2))
     drl(i,Kmax(i))=rl(i,Kmax(i))-rl(i,Kmax(i)-1)
     drl(i,1)=rl(i,2)-rl(i,1)
  end do
  do k=1,Ko2
     sr(1:Ir,k)=sin(rl(1:Ir,k))
     cr(1:Ir,k)=cos(rl(1:Ir,k))
     sr3(1:Ir,k)=1.+3.*sr(1:Ir,k)**2
     BBr(1:Ir,k)=cr(1:Ir,k)**6/sqrt(sr3(1:Ir,k))
  end do
  do k=1,Ko2
     do j=1,Jo
        r2(1,1:Ir,j,k)=Re*Lsh(1:Ir)*cr(1:Ir,k)**3*cp(j)
        r2(2,1:Ir,j,k)=Re*Lsh(1:Ir)*cr(1:Ir,k)**3*sp(j)
        r2(3,1:Ir,j,k)=Re*Lsh(1:Ir)*cr(1:Ir,k)**2*sr(1:Ir,k)
     end do
  end do
  do j=1,Jo
     do i=1,Ir
        r2(1,i,j,Kmax(i))=Rionos*Re*cl(i)*cp(j)
        r2(2,i,j,Kmax(i))=Rionos*Re*cl(i)*sp(j)
        r2(3,i,j,Kmax(i))=Rionos*Re*sl(i)
     end do
  end do
  do k=1,Ko2 
     R22(1:Ir,k)=(Re*Lsh(1:Ir))**2*cr(1:Ir,k)**4
     Rxy(1:Ir,k)=Re*Lsh(1:Ir)*cr(1:Ir,k)**3
     do j=1,Jo
        do i=1,Ir
           dBdrB(1,i,j,k)=-r2(1,i,j,k)/R22(i,k)*(3.+sr(i,k)**2/sr3(i,k))
           dBdrB(2,i,j,k)=-r2(2,i,j,k)/R22(i,k)*(3.+sr(i,k)**2/sr3(i,k))
           dBdrB(3,i,j,k)=-1./R22(i,k)*(-3.*r2(3,i,j,k)    &
     		+Rxy(i,k)*sr(i,k)*cr(i,k)/sr3(i,k))
        end do
     end do
  end do
  do k=1,Ko2
     Bz(1:Ir,k)=ME*(1.-3.*sr(1:Ir,k)**2)/R22(1:Ir,k)**1.5
     Bxy(1:Ir,k)=3.*abs(cr(1:Ir,k)*sr(1:Ir,k))*ME/R22(1:Ir,k)**1.5
     Bf2(1:Ir,k)=Bz(1:Ir,k)**2+Bxy(1:Ir,k)**2
  end do
  j=1
  do i=1,Ir 
     ds1(i,1)=0.
     ds2(i,Kmax(i))=0.
     ds1(i,2:Kmax(i))=sqrt((r2(1,i,j,2:Kmax(i))-r2(1,i,j,1:Kmax(i)-1))**2   &
          +(r2(2,i,j,2:Kmax(i))-r2(2,i,j,1:Kmax(i)-1))**2   &
          +(r2(3,i,j,2:Kmax(i))-r2(3,i,j,1:Kmax(i)-1))**2)
     ds2(i,1:Kmax(i)-1)=sqrt((r2(1,i,j,1:Kmax(i)-1)-r2(1,i,j,2:Kmax(i)))**2   &
          +(r2(2,i,j,1:Kmax(i)-1)-r2(2,i,j,2:Kmax(i)))**2   &
          +(r2(3,i,j,1:Kmax(i)-1)-r2(3,i,j,2:Kmax(i)))**2)
     ds(i,1:Kmax(i))=.5*(ds1(i,1:Kmax(i))+ds2(i,1:Kmax(i)))
     ds1(i,1)=ds1(i,2)
     ds2(i,Kmax(i))=1.E20
  end do
  do j=1,Jo
     j1(j)=j-1
     j2(j)=j+1
     if (j.eq.Jo) j2(j)=1 
     if (j.eq.1) j1(j)=Jo
     gam1(j)=phi(j)-.5*dphi
     gam2(j)=phi(j)+.5*dphi
     sg1(j)=sin(gam1(j))
     sg2(j)=sin(gam2(j))
     cg1(j)=cos(gam1(j))
     cg2(j)=cos(gam2(j))
  end do

  do i=1,Ir
     i1(i)=i-1
     i2(i)=i+1
     if (i.eq.Ir) i2(i)=i2(i)-1
     if (i.eq.1) i1(i)=1
     dRm(i)=.5
     if (i.eq.1 .or. i.eq.Ir) then 
        dRm(i)=1.
     endif
  end do
  j=1  ! Pick one, r2(3) is independent of j
  do i=1,Ir
     do k=1,Kmax(i)
        k1(i,k)=k-1
        k2(i,k)=k+1
        if (k.eq.Kmax(i)) k2(i,k)=k2(i,k)-1
        if (k.eq.1) k1(i,k)=1
        beta1(i,k)=asin((Bz(i,k)+Bz(i,k1(i,k)))/(sqrt(Bf2(i,k))+sqrt(Bf2(i,k1(i,k)))))
        beta2(i,k)=asin((Bz(i,k)+Bz(i,k2(i,k)))/(sqrt(Bf2(i,k))+sqrt(Bf2(i,k2(i,k)))))
        sb1(i,k)=sin(beta1(i,k))
        sb2(i,k)=sin(beta2(i,k))
        cb1(i,k)=cos(beta1(i,k))
        cb2(i,k)=cos(beta2(i,k))
        alpha1(i,k)=.5*pi - asin(Bz(i,k)/sqrt(Bf2(i,k)))
        sa1(i,k)=sin(alpha1(i,k))
        ca1(i,k)=cos(alpha1(i,k))
        delR(i,k)=dR*BBr(i,k)/cr(i,k)**3
        ikk1(i,k)=k
        ikk2(i,k)=k
        ik1(i,k)=0
        ik2(i,k)=0
        if (i.eq.1 .or. k.gt.Kmax(i)) then 
           fac1(i,k)=0. 
        else
           lam1=atan((r2(3,i,j,k)-delR(i,k)*sa1(i,k))/   &
                (Rxy(i,k)-delR(i,k)*ca1(i,k)))
           ikk1(i,k)=-1
           kk=1
           do while (ikk1(i,k).eq.-1)
              kk=kk+1
              if (rl(i1(i),kk).ge.lam1) then 
                 ikk1(i,k)=kk 
              else
                 if (kk.eq.Kmax(i1(i))) ikk1(i,k)=kk
              endif
           end do
           if (ikk1(i,k).gt.1) ik1(i,k)=1
           fac1(i,k)=amax1(0.,amin1(1.,(rl(i1(i),ikk1(i,k))-lam1)/   &
                (rl(i1(i),ikk1(i,k))-rl(i1(i),ikk1(i,k)-1))))
        end if
        if (i.eq.Ir .or. k.gt.Kmax(i)) then 
           fac2(i,k)=0. 
        else
           lam2=atan((r2(3,i,j,k)+delR(i,k)*sa1(i,k))/    &
                (Rxy(i,k)+delR(i,k)*ca1(i,k)))
           ikk2(i,k)=-1
           kk=1
           do while (ikk2(i,k).eq.-1)
              kk=kk+1
              if (rl(i,kk).ge.lam2) then 
                 ikk2(i,k)=kk 
              else 
                 if (kk.eq.Kmax(i2(i))) ikk2(i,k)=kk
              endif
           end do
           if (ikk2(i,k).gt.1) ik2(i,k)=1
           fac2(i,k)=amax1(0.,amin1(1.,(rl(i,ikk2(i,k))-lam2)/   &
                (rl(i,ikk2(i,k))-rl(i,ikk2(i,k)-1))))
        end if
     end do
  end do
  do i=1,Ir
     do k=1,Kmax(i) 
        dR1(i,k)=.5*(delR(i,k)+dR*((1.-fac1(i,k))*   &
      	     sqrt(BBr(i1(i),ikk1(i,k)))+fac1(i,k)*   &
      	     sqrt(BBr(i1(i),ikk1(i,k)-ik1(i,k)))))
        dR2(i,k)=.5*(delR(i,k)+dR*(fac2(i,k)*   &
      	     sqrt(BBr(i2(i),ikk2(i,k)-ik2(i,k)))    &
      	     +(1.-fac2(i,k))*sqrt(BBr(i2(i),ikk2(i,k)))))
     end do
  end do

  Irfac=0
  Ilfac=0
  do j=1,Jo
     do I=1,Ir
        FPOT(I,J)=0.
     end do
  end do
  do j=1,Jo
     do I=1,Ir
        BASEPOT(I,J)=0.
     end do
  end do
  do S=1,NS
     do j=1,Jo
        do I=1,Ir
           Iphi(I,J,S)=0.
        end do
     end do
  end do
  do S=1,NS
     do j=1,Jo
        do I=1,Ir
           Irad(I,J,S)=0.
        end do
     end do
  end do
end subroutine CURRENTSETUP
!=======================================================================
!			 	CURRENTCALC
!     Routine calculates the perpendicular "ring current" and the
!      field-aligned "partial ring current" closure currents
!=======================================================================
subroutine CURRENTCALC

  use ModHeidiSize
  use ModHeidiCurrents
  use ModHeidiMain
  use ModHeidiIO

  implicit none

  real :: As1(NR+3,NT),Tf(NT),Tc(NT),Jc(3,NR+3,NT,Slen),Pzero,Nzero,   &
       dPi1,dPi2,dPj1,dPj2,dPk1,dPk2,dPx,dPy,dPz,Jx1,Jx2,Jy1,Jy2,c0,   &
       Jz1,Jz2,PsubP,Funcpc,Funcpf,Ji1(3),Ji2(3),Jj1(3),Jj2(3),kfac,   &
       Js1,Js2,Js3,Js4,Js,J3

  real     :: pf1(NR+3,NT),pc1(NR+3,NT),nr1(NR+3,NT)
  integer  :: i,j,k,m,Imax
  external :: Funcpc,Funcpf

  !---------------------------------------------------------------------   

  !  Initialize a few numbers
  Pzero=1.E-8
  Nzero=1.E-8
  !Irad(1:IO,1:JO,1:NS)=0.   ! Radially outward current totals (A)
  !Iphi(1:IO,1:JO,1:NS)=0.   ! Eastward current totals (A)
  !Jion1(1:IO,1:JO,1:NS)=0.  ! Current into ionosphere (A/m2)

  Irad=0.   ! Radially outward current totals (A)
  Iphi=0.   ! Eastward current totals (A)
  Jion1=0.  ! Current into ionosphere (A/m2)


  !  Start main species loop
  do S=1,NS
     if (SCALC(S).eq.1) then

        !  Calculate bulk values for an empty loss cone distribution
	do j=1,Jo
           Pf1(1:Ir,j)=Pzero
           Pc1(1:Ir,j)=Pzero
           Nr1(1:Ir,j)=Nzero
           Pf1(2:IO+1,j)=PPAR(1:IO,j,S)
           Pc1(2:IO+1,j)=PPER(1:IO,j,S)
           Nr1(2:IO+1,j)=RNHT(1:IO,j,S)
           do i=2,1,-1
              Pf1(i,j)=Pf1(i+1,j)
              Pc1(i,j)=Pc1(i+1,j)*amin1(.25,Pc1(i+1,j)/(Pc1(i+2,j)+1e-30))
              Nr1(i,j)=Nr1(i+1,j)*amin1(.25,Nr1(i+1,j)/(Nr1(i+2,j)+1e-30))
           end do
           do i=ILMP(J)+1,Ir
              Pf1(i,j)=Pf1(i-1,j)*(Lsh(i-1)/Lsh(i))**3
              PC1(i,j)=PC1(i-1,j)*(Lsh(i-1)/Lsh(i))**6
              Nr1(i,j)=Nr1(i-1,j)*(Lsh(i-1)/Lsh(i))**3
           end do
           do i=1,Ir
              Pf1(i,j)=amax1(Pf1(i,j),Pzero)  ! in case model result is 0
              Pc1(i,j)=amax1(Pc1(i,j),Pzero)
              Nr1(i,j)=amax1(Nr1(i,j),Nzero)
              c0=1./(1./BBr(i,Kmax(i))-1.)
              Tf(j)=Pf1(i,j)/Nr1(i,j)*(1.+Pf1(i,j)/(Pc1(i,j)+1e-30)*c0)
              Tc(j)=(Pc1(i,j)+Pf1(i,j)*c0)/(Nr1(i,j)+1e-30)
              Pf1(i,j)=Nr1(i,j)*Tf(j)
              Pc1(i,j)=Nr1(i,j)*Tc(j)
              As1(i,j)=Pf1(i,j)/(Pc1(i,j)+1e-30) - 1.
           end do
	end do

        ! Jperp calculation
	do j=1,Jo
           do i=1,Ir
              do k=1,Kmax(i)
                 dPi1=Pf1(i2(i),j)*(fac2(i,k)*Funcpc(As1(i2(i),j),   &
                      BBr(i2(i),ikk2(i,k)-ik2(i,k)),BBr(i2(i),Kmax(i2(i))))   &
                      +(1.-fac2(i,k))*Funcpc(As1(i2(i),j),BBr(i2(i),   &
                      ikk2(i,k)),BBr(i2(i),Kmax(i2(i)))))-Pf1(i,j)*Funcpc(As1(i,j),   &
                      BBr(i,k),BBr(i,Kmax(i)))
                 dPi2=Pf1(i,j)*Funcpc(As1(i,j),BBr(i,k),BBr(i,Kmax(i)))   &
                      -Pf1(i1(i),j)*(fac1(i,k)*Funcpc(As1(i1(i),j),   &
                      BBr(i1(i),ikk1(i,k)-ik1(i,k)),BBr(i1(i),Kmax(i1(i))))   &
                      +(1.-fac1(i,k))*Funcpc(As1(i1(i),j),BBr(i1(i),ikk1(i,k))   &
                      ,BBr(i1(i),Kmax(i1(i)))))
                 dPj1=Pf1(i,j)*Funcpc(As1(i,j),BBr(i,k),BBr(i,Kmax(i)))   &
                      -Pf1(i,j1(j))*Funcpc(As1(i,j1(j)),BBr(i,k),BBr(i,Kmax(i)))
                 dPj2=Pf1(i,j2(j))*Funcpc(As1(i,j2(j)),BBr(i,k),BBr(i,Kmax(i)))   &
                      -Pf1(i,j)*Funcpc(As1(i,j),BBr(i,k),BBr(i,Kmax(i)))
                 dPk1=Pf1(i,j)*Funcpc(As1(i,j),BBr(i,k),BBr(i,Kmax(i)))   &
                      -Pf1(i,j)*Funcpc(As1(i,j),BBr(i,k1(i,k)),BBr(i,Kmax(i)))
                 dPk2=Pf1(i,j)*Funcpc(As1(i,j),BBr(i,k2(i,k)),BBr(i,Kmax(i)))   &
                      -Pf1(i,j)*Funcpc(As1(i,j),BBr(i,k),BBr(i,Kmax(i)))
                 dPk1=0.   ! no influence B x dP since they are along B!
                 dPk2=0.
                 dPx=dRm(i)*ca1(i,k)*cp(j)*(dPi1/dR1(i,k)+dPi2/dR2(i,k))    &
                      -.5*(dPj1*sg1(j)+dPj2*sg2(j))/(Rxy(i,k)*dphi)   !&
                 !      	      -.5*dPk1*cb1(i,k)*cp(j)/ds1(i,k)   &
                 !      	      -.5*dPk2*cb2(i,k)*cp(j)/ds2(i,k)
                 dPy=dRm(i)*ca1(i,k)*sp(j)*(dPi1/dR1(i,k)+dPi2/dR2(i,k))    &
                      +.5*(dPj1*cg1(j)+dPj2*cg2(j))/(Rxy(i,k)*dphi)   !&
                 !      	      -.5*dPk1*cb1(i,k)*sp(j)/ds1(i,k)    &
                 !      	      -.5*dPk2*cb2(i,k)*sp(j)/ds2(i,k)
                 dPz=dRm(i)*sa1(i,k)*(dPi1/dR1(i,k)+dPi2/dR2(i,k))   !&
                 !      	      +.5*dPk1*sb1(i,k)/ds1(i,k)    !&
                 !      	      +.5*dPk2*sb2(i,k)/ds2(i,k) 
                 Jx1=-Bxy(i,k)*sp(j)/Bf2(i,k)*dPz - Bz(i,k)/Bf2(i,k)*dPy
                 Jy1=Bz(i,k)/Bf2(i,k)*dPx + Bxy(i,k)*cp(j)/Bf2(i,k)*dPz
                 Jz1=-Bxy(i,k)*cp(j)/Bf2(i,k)*dPy   &
                      +Bxy(i,k)*sp(j)/Bf2(i,k)*dPx
                 PsubP=Pf1(i,j)*Funcpf(As1(i,j),BBr(i,k),BBr(i,Kmax(i)))   &
                      -Pf1(i,j)*Funcpc(As1(i,j),BBr(i,k),BBr(i,Kmax(i)))
                 Jx2=-Bxy(i,k)*sp(j)/Bf2(i,k)*PsubP*dBdrB(3,i,j,k)   &
                      -Bz(i,k)/Bf2(i,k)*PsubP*dBdrB(2,i,j,k)
                 Jy2=Bz(i,k)/Bf2(i,k)*PsubP*dBdrB(1,i,j,k)   &
                      +Bxy(i,k)*cp(j)/Bf2(i,k)*PsubP*dBdrB(3,i,j,k)
                 Jz2=-Bxy(i,k)*cp(j)/Bf2(i,k)*PsubP*dBdrB(2,i,j,k)   &
                      +Bxy(i,k)*sp(j)/Bf2(i,k)*PsubP*dBdrB(1,i,j,k)
                 Jc(1,i,j,k)=Jx1+Jx2   ! save Jperp for Jpara calculation
                 Jc(2,i,j,k)=Jy1+Jy2
                 Jc(3,i,j,k)=Jz1+Jz2
                 !* Integrate for total current output  (A) (doubled: N&S hemispheres)
                 if (i.ge.3 .and. i.le.ILMP(J)) then
                    Irad(i,j,s)=Irad(i,j,s)+2.*.16E-9*Re*cr(i,k)**3*Lsh(i)*dphi   &
                         *ds(i,k)*((Jc(1,i,j,k)*cp(j)+Jc(2,i,j,k)*sp(j))   &
                         *ca1(i,k)+Jc(3,i,j,k)*sa1(i,k))
                    Iphi(i,j,s)=Iphi(i,j,s)+2.*.16E-9*delR(i,k)*ds(i,k)*   &
                         (-Jc(1,i,j,k)*sp(j)+Jc(2,i,j,k)*cp(j))
                    if (Iphi(i,j,s)-Iphi(i,j,s).ne.0.) then
                       !print *,'Iphi:',Iphi(i,j,s),s,i,j,k,Kmax(i),Jc(1,i,j,k),   &
                       !     Jc(2,i,j,k),Jc(3,i,j,k),Jx1,Jx2,Jy1,Jy2,Jz1,Jz2
                       ! print *,'dP#:',dPx,dPy,dPz,dPi1,dPi2,dPj1,dPj2,dPk1,dPk2
                       ! print *,'dPi:',i1(i),i2(i),Pf1(i2(i),j),Pf1(i1(i),j),fac2(i,k),fac1(i,k)
                       !print *,'Funcpc:',Funcpc(As1(i,j),BBr(i,k),BBr(i,Kmax(i))),BBr(i,k),   &
                       !     BBr(i,Kmax(i)),As1(i,j)
                       !print *,'FuncPi1:',Pf1(i2(i),j),fac2(i,k),Funcpc(As1(i2(i),j),   &
                       !     BBr(i2(i),ikk2(i,k)-ik2(i,k)),BBr(i2(i),Kmax(i2(i)))),   &
                       !     1.-fac2(i,k),Funcpc(As1(i2(i),j),BBr(i2(i),   &
                       !     ikk2(i,k)),BBr(i2(i),Kmax(i2(i)))),Pf1(i,j),   &
                       !     Funcpc(As1(i,j),BBr(i,k),BBr(i,Kmax(i)))
                       call CON_stop('ERROR in heidi_currents.f90')
                    end if
                 end if  !! Only save where we have pressures
              end do  ! k loop
           end do   ! i loop
	end do    ! j loop

        !  Begin Jparallel calculation
	do j=1,Jo
           Imax=ILMP(J)-1
           if (J.ge.2 .and.J.le.Jo-1) then
              Imax=min(Imax,min(ILMP(J-1)-1,ILMP(J+1)-1))
           end if
           do i=2,Imax
              J3=0.  ! reset for each hemispherical calculation
              do k=1,Kmax(i)
                 kfac=1.
                 if (k.eq.1) kfac=.5
                 do m=1,3   ! 3 components of Jpara 
                    Ji1(m)=.5*(fac1(i,k)*Jc(m,i1(i),j,ikk1(i,k)-ik1(i,k))   &
                         +(1.-fac1(i,k))*Jc(m,i1(i),j,ikk1(i,k))+Jc(m,i,j,k))
                    Ji2(m)=.5*(fac2(i,k)*Jc(m,i2(i),j,ikk2(i,k)-ik2(i,k))   &
                         +(1.-fac2(i,k))*Jc(m,i2(i),j,ikk2(i,k))+Jc(m,i,j,k))
                    Jj1(m)=.5*(Jc(m,i,j,k)+Jc(m,i,j1(j),k))
                    Jj2(m)=.5*(Jc(m,i,j2(j),k)+Jc(m,i,j,k))
                 end do  ! m loop
                 Js1=Ji1(1)*ca1(i,k)*cp(j)+Ji1(2)*ca1(i,k)*sp(j)   &
                      +Ji1(3)*sa1(i,k)  ! positive for flow into volume
                 Js2=-Ji2(1)*ca1(i,k)*cp(j)-Ji2(2)*ca1(i,k)*sp(j)   &
                      -Ji2(3)*sa1(i,k)
                 Js3=-Jj1(1)*sg1(j)+Jj1(2)*cg1(j)
                 Js4=Jj2(1)*sg2(j)-Jj2(2)*cg2(j)
                 Js=Js1*.5*(ds(i,k)+fac1(i,k)*ds(i1(i),ikk1(i,k)-ik1(i,k))+   &
                      (1.-fac1(i,k))*ds(i1(i),ikk1(i,k)))   &
                      *(Lsh(i)*Re-.5*delR(i,k))*cr(i,k)**3*dphi   &
                      +Js2*.5*(ds(i,k)+fac2(i,k)*ds(i2(i),ikk2(i,k)-ik2(i,k))+   &
                      (1.-fac2(i,k))*ds(i2(i),ikk2(i,k)))   &
                      *(Lsh(i)*Re+.5*delR(i,k))*cr(i,k)**3*dphi+(Js3+Js4)   &
                      *ds(i,k)*delR(i,k)
                 Js=Js/(delR(i,k)*Lsh(i)*Re*cr(i,k)**3*dphi)
                 J3=J3*BBr(i,k1(i,k))/BBr(i,k)+kfac*Js
              end do ! k loop
              Jion1(i,j,s)=0.16E-9*J3    ! A/m2 into one hemisphere
           end do  ! i loop
	end do   ! j loop

     end if   ! SCALC check
  end do   ! S loop

end subroutine CURRENTCALC
!=======================================================================
!                        Function Funcpf
!  This function calculates the parallel pressure of the distribution 
!  function at some point along the field line
!=======================================================================
real function Funcpf(As,BBr,BBm)

  implicit none

  real:: c1, c2, BBr, BBm, As
  !---------------------------------------------------------------------  
  c1=sqrt(amax1(0.,(BBr/BBm-1.)/(BBr/BBm-BBr)))
  c2=(1.+As*BBr)/(BBr/BBm+As*BBr)
  FuncPf=(As+1.)/(1.+As*BBr)*c1*(1.-c2)

end function Funcpf
!=======================================================================
!                        Function Funcpc 
!  This function calculates the perpendicular pressure of the 
!  distribution function at some point along the field line 
!=======================================================================
real function Funcpc(As,BBr,BBm)

  implicit none

  real:: c1, c2, BBr, BBm, As
  !--------------------------------------------------------------------    
  c1=sqrt(amax1(0.,(BBr/BBm-1.)/(BBr/BBm-BBr)))
  c2=(1.+As*BBr)/(BBr/BBm+As*BBr)
  FuncPc=(As+1.)/(1.+As*BBr)**2*c1*(1.-c2)

end function Funcpc


