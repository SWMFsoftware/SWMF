!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!======================================================================
!			       heidi_cepara
! Calculates the cross-section of charge exchange of ring current
! ion species with the neutral hydrogen and the charge exchange rate
! Also calculates the atmospheric loss rate (out the loss cone)
!======================================================================
subroutine heidi_cepara

  use ModHeidiSize,   ONLY: nR, nT, nPa, nS, s, io, jo, lo, ko, scalc, dt
  use ModHeidiMain,   ONLY: eKev, upa, mu, lz, funt, m1, NeutralHydrogen, v, IsBFieldNew
  use ModHeidiDrifts, ONLY: acharge, atlos
  use ModIoUnit,      ONLY: UNITTMP_
  use ModHeidiIO,     ONLY: NameInputDirectory
  use ModHeidiInput,  ONLY: TypeBCalc

  implicit none

  integer            :: I, K, L, LUP, ier, J
  real               :: x, y, fac, cosd 
  real               :: RLAMBDA(71), PA(71)
  real               :: HDNS(NR,NT,NPA) = 0.0
  real               :: LH(20), HDNSIN(20,71)
  character(len=80)  :: TITLE
  external           :: cosd

  integer :: iR
  !---------------------------------------------------------------------
  !\
  ! Open a file with the bounce-averaged H dens [m-3] in geocorona
  ! obtained through bounce-averaging the Chamberlain model
  ! (Mei-Ching's fit by a poly of order 4 good for L from 1 to 6.5)
  !/

  s = 0
  do i = 2, nS
     s = s + sCalc(I)   ! S=0, no ions in calc, S>0, ions in calc
  end do
  if (s .ge. 1) then	! Only needed for ion charge exchange

     select case(TypeBCalc)

     case('analytic')
        open(UNITTMP_,file=NameInputDirectory//'hgeo71.in',STATUS='OLD')
        LUP=71
        do I = 1,20
           read(UNITTMP_,3) TITLE,LH(i)
           do L = 1,LUP
              read(UNITTMP_,4) PA(L),RLAMBDA(L),HDNSIN(I,L)
           end do	! L loop
        end do	        ! I loop
        close(UNITTMP_)
        do L = 1, LUP
           PA(L) = COSD(PA(L))
        end do
        do L=1,LO
           do j = 1, jo
              do I=1,IO
                 call heidi_lintp2(LH,PA,HDNSIN,20,71,LZ(i),MU(L),FAC,IER)
                 HDNS(i,j,l)=FAC
              enddo
           end do
        end do

     case('numeric')

        if (IsBFieldNew) then 
           write(*,*) 'heidi_coefficients: heidi_cepara--->get_neutral_hydrogen'
           call get_neutral_hydrogen(NeutralHydrogen)
        end if

        do L = 1, LO
           do j = 1, jo
              do I = 1, IO
                 HDNS(i,j,l) = NeutralHydrogen(i,j,l)
              enddo
           end do
        end do
     end select

  end if



3 format(A61,F5.2)
4 format(2F7.3,2X,1PE12.5)
  !\
  ! Calculate charge exchange cross-section of species S with H
  ! and then the charge exchange decay rate ACHARGE
  !/

  !\
  ! H+ charge exchange
  !/

  if (SCALC(2).eq.1) then
     do K = 2,KO
        X=ALOG10(EKEV(K))
        if (X.lt.-2.) X=-2.
        Y=-18.767-0.11017*X-3.8173e-2*X**2-0.1232*X**3-5.0488e-2*X**4
        !............10**Y is cross section of H+ in m2
        
        do j = 1, jo
           do L=2,Lo
              acharge(2:io,j,k,l,2)=exp(-(10.**Y*V(K,2)*HDNS(2:io,j,L)*DT))
           end do ! L loop
        end do    ! j loop
     end do	  ! K loop

!!$     open(unit=3, file='H+ChargeDecayRate_3PA_Phi90_E10keV_LAll.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to H+ with H'
!!$     write(3,*) '  Y'
!!$
!!$     do K = 27, 27
!!$        do iR = 1, nR
!!$           write(3,*) LZ(iR), acharge(iR,6,k,2,2),  acharge(iR,1,6,20,2),  acharge(iR,6,k,42,2)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$     open(unit=3, file='H+ChargeDecayRate_3PA_Phi0_E10keV_LAll.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to H+ with H'
!!$     write(3,*) '  Y'
!!$
!!$     do K = 27, 27
!!$        do iR = 1, nR
!!$           write(3,*) LZ(iR), acharge(iR,1,k,2,2),  acharge(iR,1,k,20,2),  acharge(iR,1,k,42,2)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$
!!$ open(unit=3, file='H+ChargeDecayRate_3PA_Phi90_E107keV_LAll.dat', status='unknown')
!!$  write(3,*) 'Decay Rate due to H+ with H'
!!$  write(3,*) '  Y'
!!$  
!!$  do K = 37, 37
!!$     do iR = 1, nR
!!$        write(3,*) LZ(iR), acharge(iR,6,k,2,2),  acharge(iR,1,6,20,2),  acharge(iR,6,k,42,2)
!!$     end do
!!$  end do	  ! K loop
!!$close(3)
!!$
!!$
!!$  open(unit=3, file='H+ChargeDecayRate_3PA_Phi0_E107keV_LAll.dat', status='unknown')
!!$  write(3,*) 'Decay Rate due to H+ with H'
!!$  write(3,*) '  Y'
!!$  
!!$  do K = 37, 37
!!$     do iR = 1, nR
!!$        write(3,*) LZ(iR), acharge(iR,1,k,2,2),  acharge(iR,1,k,20,2),  acharge(iR,1,k,42,2)
!!$     end do
!!$  end do	  ! K loop
!!$close(3)
!!$
!!$
!!$open(unit=3, file='H+ChargeDecayRate_3PA_Phi180_E107keV_LAll.dat', status='unknown')
!!$  write(3,*) 'Decay Rate due to H+ with H'
!!$  write(3,*) '  Y'
!!$  
!!$  do K = 37, 37
!!$     do iR = 1, nR
!!$        write(3,*) LZ(iR), acharge(iR,12,k,2,2),  acharge(iR,12,k,20,2),  acharge(iR,12,k,42,2)
!!$     end do
!!$  end do	  ! K loop
!!$close(3)
!!$
!!$
!!$open(unit=3, file='H+ChargeDecayRate_3PA_Phi180_Eall_L5Re.dat', status='unknown')
!!$write(3,*) 'Decay Rate due to H+ with H'
!!$write(3,*) '  Y'
!!$
!!$do K = 1, ko
!!$     do iR = 14, 14
!!$        write(3,*) Ekev(k), acharge(iR,12,k,2,2),  acharge(iR,12,k,20,2),  acharge(iR,12,k,42,2)
!!$     end do
!!$  end do	  ! K loop
!!$close(3)
!!$
!!$
!!$open(unit=3, file='H+ChargeDecayRate_3PA_Phi0_Eall_L5Re.dat', status='unknown')
!!$write(3,*) 'Decay Rate due to H+ with H'
!!$write(3,*) '  Y'
!!$
!!$do K = 1, ko
!!$     do iR = 14, 14
!!$        write(3,*) Ekev(k), acharge(iR,1,k,2,2),  acharge(iR,1,k,20,2),  acharge(iR,1,k,42,2)
!!$     end do
!!$  end do	  ! K loop
!!$close(3)
!!$
!!$
!!$
!!$open(unit=3, file='H+ChargeDecayRate_3PA_Phi90_Eall_L5Re.dat', status='unknown')
!!$write(3,*) 'Decay Rate due to H+ with H'
!!$write(3,*) '  Y'
!!$
!!$do K = 1, ko
!!$     do iR = 14, 14
!!$        write(3,*) Ekev(k), acharge(iR,6,k,2,2),  acharge(iR,6,k,20,2),  acharge(iR,6,k,42,2)
!!$     end do
!!$  end do	  ! K loop
!!$close(3)
!!$

  end if

  !\
  ! He+ charge exchange
  !/

  if (SCALC(3).eq.1) then
     do K=2,KO
        X=ALOG10(EKEV(K))
        if (X.lt.-2.) X=-2.
        Y=-20.789+0.92316*X-0.68017*X**2+0.66153*X**3-0.20998*X**4
        !..........10**Y is cross sect of He+ in m2
        do j = 1, jo
           do L=2,Lo
              acharge(2:io,j,K,L,3)=exp(-(10.**Y*V(K,3)*HDNS(2:io,j,L)*DT))
           end do ! L loop
        end do    ! j loop
     end do	  ! K loop
  end if

  !\
  ! O+ charge exchange
  !/

  if (SCALC(4).eq.1) then
     do K=2,KO
        X=ALOG10(EKEV(K))
        if (X.lt.-2.) X=-2.
        Y=-18.987-0.10613*X-5.4841E-3*X**2-1.6262E-2*X**3   &
             -7.0554E-3*X**4
        ! 10**Y is cross sect of O+ in m2
        do j =1, jo
           do L=2,Lo
              acharge(2:io,j,k,l,4)=exp(-(10.**Y*V(K,4)*HDNS(2:io,j,L)*DT))
           end do  ! L loop
        end do     ! j loop
     end do	   ! K loop

  
!!$     open(unit=3, file='O+ChargeDecayRate_3PA_Phi0_E10keV_LAll.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to O+ with H'
!!$     write(3,*) '  Y'
!!$
!!$     do K = 27, 27
!!$        do iR = 1, nR
!!$           write(3,*) LZ(iR), acharge(iR,1,k,2,4),  acharge(iR,1,k,20,4),  acharge(iR,1,k,42,4)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$
!!$     open(unit=3, file='O+ChargeDecayRate_3PA_Phi90_E10keV_LAll.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to O+ with H'
!!$     write(3,*) '  Y'
!!$
!!$     do K = 27, 27
!!$        do iR = 1, nR
!!$           write(3,*) LZ(iR), acharge(iR,6,k,2,4),  acharge(iR,6,k,20,4),  acharge(iR,6,k,42,4)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$
!!$
!!$     open(unit=3, file='O+ChargeDecayRate_3PA_Phi0_E107keV_LAll.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to O+ with H'
!!$     write(3,*) '  Y'
!!$     
!!$     do K = 37, 37
!!$        do iR = 1, nR
!!$           write(3,*) LZ(iR), acharge(iR,1,k,2,4),  acharge(iR,1,k,20,4),  acharge(iR,1,k,42,4)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$     
!!$     
!!$     open(unit=3, file='O+ChargeDecayRate_3PA_Phi180_E107keV_LAll.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to O+ with H'
!!$     write(3,*) '  Y'
!!$     
!!$     do K = 37, 37
!!$        do iR = 1, nR
!!$           write(3,*) LZ(iR), acharge(iR,12,k,2,4),  acharge(iR,12,k,20,4),  acharge(iR,12,k,42,4)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$    
!!$     
!!$     open(unit=3, file='O+ChargeDecayRate_3PA_Phi90_E107keV_LAll.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to O+ with H'
!!$     write(3,*) '  Y'
!!$     
!!$     do K = 37, 37
!!$        do iR = 1, nR
!!$           write(3,*) LZ(iR), acharge(iR,6,k,2,4),  acharge(iR,6,k,20,4),  acharge(iR,6,k,42,4)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$     
!!$     open(unit=3, file='O+ChargeDecayRate_3PA_Phi180_Eall_L5Re.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to O+ with H'
!!$     write(3,*) '  Y'
!!$     
!!$     do K = 1, ko
!!$        do iR = 14, 14
!!$           write(3,*) Ekev(k), acharge(iR,12,k,2,4),  acharge(iR,12,k,20,4),  acharge(iR,12,k,42,4)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$     
!!$     
!!$     open(unit=3, file='O+ChargeDecayRate_3PA_Phi0_Eall_L5Re.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to O+ with H'
!!$     write(3,*) '  Y'
!!$     
!!$     do K = 1, ko
!!$        do iR = 14, 14
!!$           write(3,*) Ekev(k), acharge(iR,1,k,2,4),  acharge(iR,1,k,20,4),  acharge(iR,1,k,42,4)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$     
!!$     open(unit=3, file='O+ChargeDecayRate_3PA_Phi90_Eall_L5Re.dat', status='unknown')
!!$     write(3,*) 'Decay Rate due to O+ with H'
!!$     write(3,*) '  Y'
!!$     
!!$     do K = 1, ko
!!$        do iR = 14, 14
!!$           write(3,*) Ekev(k), acharge(iR,6,k,2,4),  acharge(iR,6,k,20,4),  acharge(iR,6,k,42,4)
!!$        end do
!!$     end do	  ! K loop
!!$     close(3)
!!$

  end if

  !\
  ! Calculate the losses due to the collis with atmosphere
  !/

  do s = 1, NS
     do I=2,IO  
        do L = UPA(I),LO
           do K=2,KO
              do j = 1, jo
                 ATLOS(I,J,K,L,S)=exp(-DT/(29.1347*LZ(I)*FUNT(L,I,J)*sqrt(M1(S)/EKEV(K))))
              end do	! I loop
           end do	! J loop
        end do		! K loop
     end do		! L loop
  end do                ! S loop

end subroutine Heidi_cepara
!======================================================================
!				OTHERPARA
!	  Calculate parameters used in drifts and Coulomb drags
!======================================================================
subroutine OTHERPARA

  use ModHeidiSize,   ONLY: nE, nPA, nS, nSth, s, io, jo, ko, lo, dT, scalc
  use ModHeidiIO,     ONLY: rns, rnl, res, rel, esn, eln, ele, ecn, ece, aln, &
       ale, cen, cee, ilmp, lmp, ese
  use ModHeidiMain,   ONLY: IsBFieldNew, Re, dPhi, Me, dl1, mp, q, cedr, cidr, &
       mas, dE, vbnd, ebnd, facmu, VPhi_IIII, VR_IIII, funt, funi, ekev, upa,  &
       lz, z, mlt, mu, wmu, dmu, v
  use ModHeidiDrifts, ONLY: VrConv, Vr, P2, j6, j18, coule, couli, &
       gtae, gtai, atae, atai, btae, btai
  use ModHeidiInput,  ONLY: TypeBCalc
  use ModNumConst,    ONLY: cPi
  implicit none

  integer :: i,j,k,l,is,iss,ier
  real    :: eps,dln,sqm,qe,gama,cco,ccd,edrco,x,xd,g,cci,cdi,edri,ccde,ccdi
  real    :: bane,badif,c,gpa,cce,cde,edre
  real    :: erf
  real    :: COULDE(NE,NPA),COULDI(NE,NPA),AFIR,ASEC
  real    :: VF(NSTH),RA(NSTH),MUBOUN,TMAS(NSTH),TM1(NSTH)

  external:: erf


  data TM1/5.4462E-4,1.0,4.0,16.0/  ! Thermal mass in AMU
  data RA/1.,.77,.2,.03/            ! Thermal proportions in plasmasphere
  !---------------------------------------------------------------------  

  !\
  ! Parameters used in calculating drifts at boundaries
  !/


  C=1.44E-2*RE**2		! Constant of corotation, [C]=V*m
  ISS=-1			! sign of specie's charge
  if (S.ge.2) ISS=1

  !\
  ! Determine the dayside
  !/

  J6=0
  J=0
  do while (J6.eq.0)
     J=J+1
     if (MLT(J).ge.5.99) J6=J
  end do	! While J6 loop
  J18=0
  J=0
  do while (J18.eq.0)
     J=J+1
     if (MLT(J).ge.17.99) J18=J
  end do	! While J18 loop


  if (TypeBCalc == 'numeric') then
     if (IsBFieldNew) then
        write(*,*) 'heidi_coefficients: OTHERPARA---> get_grad_curv_drift'
        call  get_grad_curv_drift(VPhi_IIII,VR_IIII)
     end if
  end if

  do i = 1, io      
     do L = 1, UPA(I)    ! Kp independent part of azimuthal drift
        do K = 1, KO
           do j = 1,  jo
              select case(TypeBCalc)
              case('analytic')
                 GPA=1.-FUNI(L,I,J)/6./FUNT(L,I,J)
                 P2(I,J,K,L)=(C-ISS*3.*EKEV(K)*1000.*Z(I)*GPA)*DT/DPHI/ME
              case('numeric')
                 P2(I,J,K,L)=(C/Me-ISS*VPhi_IIII(i,j,k,l))*DT/DPHI
              end select
           end do
        end do
     end do
  end do

  do i = 1, io 
     do L = UPA(I)+1, LO
        do K = 1, KO
           do j = 1,  jo
              select case(TypeBCalc)
              case('analytic')
                 P2(I,J,K,L)=(C-ISS*3.*EKEV(K)*1000.*Z(I)*GPA)*DT/DPHI/ME
              case('numeric')
                 P2(I,J,K,L)=(C/Me - ISS*VPhi_IIII(i,j,k,l))*DT/DPHI
              end select
           end do
        end do
     end do
  end do

  do l = 1, lo
     do k = 1, ko
        do i = 1, IO
           VrConv(I,J6,k,l)=0.
           VrConv(I,J18,k,l)=0.
           select case(TypeBCalc)
           case('numeric')
              VR(i,j6,k,l)  = VrConv(i,j6,k,l) + ISS* VR_IIII(i,j6,k,l)*DT/DL1
              VR(i,j18,k,l) = VrConv(i,j18,k,l)+ ISS* VR_IIII(i,j18,k,l)*DT/DL1
           case('analytic')
              VR(i,j6,k,l) = VrConv(I,J6,k,l)
              VR(i,j6,k,l) = VrConv(I,J6,k,l)
           end select
        end do	! I loop processing VR
     end do
  end do

  !\
  ! We assume Te=Ti=1eV (kT=1eV)
  !/

  do I = 1, NSTH
     TMAS(I) = TM1(I)*MP
     VF(I) =sqrt(2.*Q/TMAS(I))	! [m/s]
  end do	! I loop
  EPS=8.854E-12                ! [F/m]
  DLN=21.5                     ! Coulomb logarithm
  QE=(Q**2/EPS)

  do S = 1, NS
     if (SCALC(S).eq.1) then
        SQM=sqrt(MAS(S))
        ! GAMA and EDRCO divided by an extra 2 to match Khazanov defn of Q
        GAMA=DLN/8./cPi*QE*1E6/SQM*QE ! *1E6 to transform ne into m-3
        CCO=GAMA/Q*sqrt(2.)*DT
        CCD=GAMA*2*DT/16./sqrt(2.)   ! 2*DT due to double the time step
        EDRCO=2.*DLN*QE*RE/MAS(S)*QE/MAS(S)*(1.E-10) ! to invert in eV/cm2/s
        ! There was an autoscope directive here, but I removed it
        do j =1, jo
           do i =1, io
              do K=1,KO
                 X=VBND(K,S)/VF(1)
                 XD=V(K,S)/VF(1)
                 CCE=RA(1)*G(X)		     ! electrons
                 CDE=RA(1)*(ERF(XD,IER)-G(XD))
                 EDRE=RA(1)*G(XD)
                 CCI=0.
                 CDI=0.
                 EDRI=0.
                 do IS=2,NSTH
                    X=VBND(K,S)/VF(IS)
                    XD=V(K,S)/VF(IS)
                    CCI=CCI+RA(IS)*G(X)	     ! ions
                    EDRI=EDRI+RA(IS)*G(XD)
                    CDI=CDI+RA(IS)*(ERF(XD,IER)-G(XD))
                 end do	! IS loop
                 !\
                 ! Collisions with plasmaspheric e- (minus: dF/dt-d(cF)/dt=0)
                 !/
                 COULE(1,J,K,1,S)=-CCE*sqrt(EBND(K)/Q/1000.)/DE(K)*CCO
                 !\
                 ! Collisions with plasmaspheric ions
                 !/
                 COULI(1,J,K,1,S)=-CCI*sqrt(EBND(K)/Q/1000.)/DE(K)*CCO
                 !\
                 ! Energy deposition rate [eV/cm2/s]
                 !/
                 CEDR(I,J,K,1,S)=EDRCO*EKEV(K)*EDRE
                 CIDR(I,J,K,1,S)=EDRCO*EKEV(K)*EDRI
                 !\
                 ! Coulomb scattering
                 !/
                 CCDE=CCD*CDE/(EKEV(K)*Q*1000)**(1.5)
                 CCDI=CCD*CDI/(EKEV(K)*Q*1000)**(1.5)
                 COULDE(K,1)=0.
                 GTAE(I,J,K,2,S)=0.
                 COULDI(K,1)=0.
                 GTAI(I,J,K,2,S)=0.
                 do L=2,LO-1
                    BANE=(1.-FUNI(L,I,J)/2./FUNT(L,I,J))/(1.-MU(L)*MU(L))
                    COULE(1,J,K,L,S)=COULE(1,J,K,1,S)*BANE
                    COULI(1,J,K,L,S)=COULI(1,J,K,1,S)*BANE
                    CEDR(I,J,K,L,S)=CEDR(I,J,K,1,S)*BANE*FACMU(L,I,J)*WMU(L)
                    CIDR(I,J,K,L,S)=CIDR(I,J,K,1,S)*BANE*FACMU(L,I,J)*WMU(L)
                    MUBOUN=MU(L)+0.5*WMU(L)
                    BADIF=(1.-MUBOUN*MUBOUN)*(0.5*(FUNI(L+1,I,J)+FUNI(L,I,J)))/MUBOUN
                    !\
                    ! Collisions with plasmaspheric electrons
                    !/
                    COULDE(K,L)=CCDE*BADIF

                    ATAE(I,J,K,L,S)=COULDE(K,L)/FACMU(L+1,I,J)/DMU(L)/WMU(L)
                    AFIR=COULDE(K,L)/FACMU(L,I,J)/DMU(L)/WMU(L)
                    ASEC=COULDE(K,L-1)/FACMU(L,I,J)/DMU(L-1)/WMU(L)
                    BTAE(I,J,K,L,S)=AFIR+ASEC
                    GTAE(I,J,K,L+1,S)=COULDE(K,L)/FACMU(L,I,J)/WMU(L+1)/DMU(L)
                    !\
                    ! Collisions with plasmaspheric ions
                    !/
                    COULDI(K,L)=CCDI*BADIF
                    ATAI(I,J,K,L,S)=COULDI(K,L)/FACMU(L+1,I,J)/DMU(L)/WMU(L)
                    AFIR=COULDI(K,L)/FACMU(L,I,J)/DMU(L)/WMU(L)
                    ASEC=COULDI(K,L-1)/FACMU(L,I,J)/DMU(L-1)/WMU(L)
                    BTAI(I,J,K,L,S)=AFIR+ASEC
                    GTAI(I,J,K,L+1,S)=COULDI(K,L)/FACMU(L,I,J)/WMU(L+1)/DMU(L)
                 end do	! L loop
                 COULDE(K,LO)=0.
                 ATAE(I,J,K,LO,S)=0.
                 COULDI(K,LO)=0.
                 ATAI(I,J,K,LO,S)=0.
              end do	! K loop
           end do
        end do

        do K=1,KO
           do j=1, jo        
              do I=IO,1,-1		! Set to edge of loss cone when
                 do L=1,UPA(I)-1	! in the loss cone
                    COULE(I,J,K,L,S)=COULE(1,J,K,L,S)
                    COULI(I,J,K,L,S)=COULI(1,J,K,L,S)
                 end do	! L loop	
                 do L=UPA(I),LO
                    COULE(I,J,K,L,S)=COULE(I,J,K,L-1,S)
                    COULI(I,J,K,L,S)=COULI(I,J,K,L-1,S)
                 end do	! 2nd L loop
              end do	! I loop
           end do	! J loop
        end do          ! K loop
     end if	        ! SCALC check
  end do		! S loop

  !\
  ! Set continuous source/loss integrals to zero
  !/

  RNS=0.	! Radial drift particle source
  RNL=0.	! Radial drift particle loss
  RES=0.	! Radial drift energy source
  REL=0.	! Radial drift energy loss
  ESN=0.	! Energy drift particle gain
  ELN=0.	! Energy drift particle loss
  ESE=0.	! Energy drift energy gain
  ELE=0.	! Energy drift energy loss
  ECN=0.	! Coulomb collision particle gain/loss
  ECE=0.	! Coulomb collision energy gain/loss
  ALN=0.	! Atmospheric loss particle gain/loss
  ALE=0.	! Atmospheric loss energy gain/loss
  CEN=0.	! Charge exchange particle gain/loss
  CEE=0.	! Charge exchange energy gain/loss
  !\
  ! Set "magnetopause" at outer boundary
  !/

  ILMP(1:JO)=IO
  LMP(1:JO)=LZ(IO)+DL1



end subroutine OTHERPARA
!======================================================================
!                             MAGCONV
!	Calculates magnetospheric convection strength and puts it
!	into the drift terms 
!	(added to already-calculated corotation and gradient-curvature)
!  Calculation options:
!	IA=0  No convection field
!	IA=1  Kp driven V-S with Maynard and Chen activity dependence
!	IA=2  MBI driven V-S, acitivity from force balance at dusk
!	IA=3  McIlwain field
!	IA=4  McIlwain field, Eo driven by Cross Polar Cap Potential
!	IA=5  Kp V-S field plus Burke-Wygant penetration field
!	IA=6  McIlwain field plus Burke-Wygant penetration field
!	IA=7  McIlwain(CPCP) field plus B-W penetration field
!	IA=8  Kp V-S field plus self-consistent penetration field
!	IA=9  McIlwain field plus self-consistent penetration field
!	IA=10 McIlwain(CPCP) field plus S-C penetration field
!	IA=11 Unshielded V-S(CPCP) field plus penetration field
!	IA=12 Unshielded VS(CPCP) field (no Epen)
!       IA=13 W96 field plus S-C penetration field
!       IA=14 W96 field everywhere (done in S-C field block)
!       IA=15 AMIE field plus S-C penetration field
!       IA=16 AMIE field everywhere (done in S-C field block)
!       IA=20 Read in electric field from a file
!       IA=21 AMIE potentials everywhere
!       IA=22 Weimer-96 potentials everywhere
!       IA=23 Foster potentials everywhere
!       IA=?? More field models in Aaron's new subroutine? Put them in...
!  Outputs (both are unitless advection coefficients):
!	VR => Radial drift
!	P1 => Activity dependent azimuthal drift
!======================================================================

subroutine MAGCONV(I3,NST)

  use ModHeidiSize,     ONLY: nS, nR, nT, nPhiCells, nThetaCells, s, io, jo, ko, lo, &
       scalc, dT
  use ModHeidiMain,     ONLY: LZ, VPhi_IIII, VR_IIII, mlt, dPhi, Phi, IsBFieldNew, dl1, nParallelSpecies, &
       t, Me, Re, A
  use ModHeidiIO,       ONLY: NameInputDirectory, NameRun, iUnitPot, iUnitStdOut, f107, day, year, &
       BySW, BzSW, uSW, Kp, IA, iPPC, lmp, tppc, lamgam, nStep, ppc, time, write_prefix
  use ModHeidiDrifts,   ONLY: Vr, VrConv, P1
  use ModHeidiCurrents, ONLY: Jfac, Jion1, Lats, Lsh, LatFac, fpot, basepot, fpot, ilfac, irfac, &
       lonfac, iR
  use ModHeidiDGCPM,    ONLY: potdgcpm, vthetacells, vmltcells, vphicells, vlzcells
  use ModProcIM,        ONLY: iComm, iError, iProc, nProc
  use ModIoUnit,        ONLY: io_unit_new
  use ModHeidiInput,    ONLY: TypeBCalc
  use ModNumConst,      ONLY: cPi
  use ModMpi,           ONLY: MPI_Real
  implicit none

  integer :: ABASE(24),AEPEN(24),IER,j,i,ii,I4,I3,NST,jj,bc_choice
  integer :: L,IOpot,JOpot   
  real    :: KR,SPJ,CPJ,PJ,DLMAG,Ppc0,Eom,EOJ,EOJ1,NY(4),FAC,CRo
  real    :: KGAM,KS(4),KALP,KBETA,KPHI,KEX,KEY,LP,DLP,PCO,PHIPOFF
  real    :: Ro,RR,RoRg,KEPS,KSD,SJ,CJ,PHIP,EPR,EPP,DP1,FACP
  real    :: DPP,DP2,sLP,dLPdphi,Jreducer
  real    :: sind,cosd
  real    :: jfac_temp(NR+3,NT)
  real    :: FPOT1(3*NR,NT),FPOT2(3*NR,NT),FPOT12(3*NR,NT)
  real    :: TP1,TP2,LZpot(3*NR),data(3*NR),Fmin,Fmax,MLTpot(NT)
  real    :: t_sawtooth(9),tdiff,dt_saw
  character(len=80):: header
  data ABASE/1,1,1,2,3,1,2,3,1,2,3,4,4,5,5,6,6,1,1,1,10,11,12,13/
  data AEPEN/0,0,0,0,0,1,1,1,2,2,2,2,0,2,3,2,3,0,0,0, 0, 0, 0, 0/
  data KGAM/8./, KS/9.8,-1.4,-0.9,-0.3/
  data KALP/0.3/, KBETA/0.8/, KPHI/3000./
  data KEX/3.E-5/, KEY/1.2E-4/, PCO/9.17E4/
  data t_sawtooth/9.6e4,9.84e4,1.062e5,1.152e5,1.278e5,   &
       1.374e5,1.458e5,1.554e5,1.626e5/
  ! Use the next t_sawtooth definition for 4.5 h of 2*Vsw early on 4/18
  ! DATA t_sawtooth/8.64e4,8.82e4,9.0e4,9.18e4,9.36e4,
  !  &                  9.54e4,9.72e4,9.9e4,1.008e5/
  data dt_saw/600./

  external :: sind,cosd

  ! Variables to send to Aaron's subroutine
  real :: eyear, eday, ehour, eminute, esecond, eby, ebz, evsw
  integer :: nemlts1, nelats1, nemlts2, nelats2, nechoice,k
  real :: emlts1(NT,NR+3), elats1(NT,NR+3), epots1(NT,NR+3)
  real :: emlts2(nphicells,nthetacells), elats2(nphicells,nthetacells),  &
       epots2(nphicells,nthetacells)

  save I4,TP1,TP2,FPOT1,FPOT2,IOpot,LZpot,JOpot,MLTpot

  integer :: edayplus
  real    :: univ_time

  real :: ISS

  !---------------------------------------------------------------------  

  ! Initialize the field aligned current
  Jfac = 0.0


  PHIPOFF=0.
  DP1=.4*cPi
  DP2=2.*cPi/3.
  FACP=AMIN1((.25*KP)**2,1.5)
  Kr=KP/(1.+.1*KP)
  ! The next line is to cut down the size of the field-aligned currents...
  Jreducer=1.0

  ! Find polar cap potential value (if needed)
  if (ABASE(IA+1).eq.3 .or. ABASE(IA+1).eq.4) then  
     if (I3.eq.NST) then
        I4=0
        ii=2
        do while (I4.eq.0)
           if (TPPC(ii).ge.T) then
              I4=ii-1
           else if (ii.eq.IPPC-1) then
              I4=IPPC-1
           else
              ii=ii+1
           end if
        end do
     else if (T+2.*DT.ge.TPPC(IPPC)) then
        I4=IPPC-1
     else if (T.ge.TPPC(I4+1)) then
        I4=I4+1
     end if
     if (T.le.TPPC(1)) then
        PPC0=PPC(1)
     else if (T.ge.TPPC(IPPC)) then
        Ppc0=PPC(IPPC)
     else
        Ppc0=PPC(I4)+(PPC(I4+1)-PPC(I4))*(T-TPPC(I4))   &
             /AMAX1(.1*DT,TPPC(I4+1)-TPPC(I4))
     end if
  end if

  DLMAG=LMP(7)+LMP(19)     ! Width across magnetosphere at X=0
  Eom=sqrt(KEX**2+KEY**2)  ! Eo factor for McIlwain field

  if (IA.eq.0) A=0.        ! No convection field

  !\
  ! Calculate base convection electric field
  !/

  if (TypeBCalc=='numeric') then
     if (IsBFieldNew) then
        write(*,*) 'heidi_coefficients: MAGCONV---> get_grad_curv_drift'
        call  get_grad_curv_drift(VPhi_IIII,VR_IIII)
     end if
  end if

  ISS=-1			! sign of specie's charge
  if (S.ge.2) ISS=1

  if (ABASE(IA+1).eq.1) then
     write(*,*) 'This is the Volland-Stern Potential'
     do J=1,JO   ! Fill in drift values
        do i=1,io
           do l =1, lo
              do k = 1, ko
                 VrConv(i,j,k,l) = -A*cos(PHI(J))*(LZ(I)+0.5*DL1)**(LAMGAM+2.)*   &
                      DT/DL1*(RE*RE/ME)
                 select case(TypeBCalc)
                 case('analytic')
                    VR(i,j,k,l) = VrConv(i,j,k,l) 
                 case('numeric')
                    VR(i,j,k,l) = VrConv(i,j,k,l) + ISS*VR_IIII(i,j,k,l)*DT/DL1
                 end select
              end do
           end do

           P1(I,J)=A*LAMGAM*LZ(I)**(LAMGAM+1.)*sin(PHI(J)+0.5*DPHI)*   &
                DT/DPHI*(RE*RE/ME)
           BASEPOT(i+1,j)=A*RE*LZ(i)**(LAMGAM)*sin(phi(j))
        end do
        !\
        ! Compute the V-S potential outside the edges of the radial grid. 
        ! Needed for the current calculation.
        !/
        BASEPOT(1,j)    = A*RE*(LZ(1)-DL1)**(LAMGAM)*sin(phi(j))
        BASEPOT(io+2,j) = A*RE*(LZ(io)+dl1)**(LAMGAM)*sin(phi(j))
     end do

     do J = 1, nphicells   ! Fill in DGCPM potentials
        do I = 1, nthetacells
           potdgcpm(i,j)=A*RE*vlzcells(I)**(LAMGAM)*SIND(vphicells(J))
        enddo
     enddo

     !\
     ! McIlwain or modified McIlwain field
     !/

  else if (ABASE(IA+1).eq.2 .or. ABASE(IA+1).eq.3) then

     KEPS=1.+KALP*Kr   ! McIlwain E5D potential structure
     if (ABASE(IA+1).eq.3) KEPS=Ppc0/(DLMAG*Re*Eom)
     do J=1,JO    ! Fill in radial drift values
        SJ=sin(PHI(J))
        CJ=cos(PHI(J))
        Ro=KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
        if (LAMGAM.eq.0.) Ro=0.
        EOJ=RE*(KEY*SJ+KEX*CJ)
        EOJ1=RE*(KEY*CJ-KEX*SJ)
        do i=1,io
           RR=LZ(i)+0.5*DL1
           RoRg=(Ro/RR)**(KGAM-1.)
           KSD=1.+RoRg*(Ro/RR)
           do l = 1, lo
              do k = 1, ko
                 VrConv(i,j,k,l)=-DT/DL1*RE/ME*RR**3*KEPS/KSD*(EOJ1+KGAM*SJ/KSD   &
                      *(EOJ+KPHI/RR)*RoRg*KBETA/RR*(KS(2)+Kr*KS(4)))
                 select case(TypeBCalc)
                 case('analytic')
                    VR(i,j,k,l) = VrConv(i,j,k,l)
                 case('numeric')
                    VR(i,j,k,l) = VrConv(i,j,k,l) +ISS* VR_IIII(i,j,k,l)*DT/DL1
                 end select
              end do
           end do
        end do
     end do

     do J = 1, JO   ! Fill in azimuthal drift values
        SJ=sin(PHI(J)+0.5*DPHI)
        CJ=cos(PHI(J)+0.5*DPHI)
        Ro=KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
        if (LAMGAM.eq.0.) Ro=0.
        EOJ=RE*(KEY*SJ+KEX*CJ)
        do i=1,io
           RR=LZ(i)
           RoRg=(Ro/RR)**KGAM
           KSD=1.+RoRg
           P1(I,J)=DT/DPHI*RE/ME*RR**2*KEPS*EOJ/KSD*(1.+KGAM/KSD   &
                *(1.+KPHI/EOJ/RR)*RoRg)
        end do
     end do
     do j=1,jo          ! Fill in potentials for output
        SJ=sin(phi(j))
        CJ=cos(phi(j))
        EOJ=RE*(KEY*sj+KEX*cj)
        Ro=KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
        do i=1,io+3
           RoRg=(Ro/Lsh(i))**KGAM
           KSD=1.+RoRg
           BASEPOT(i,j)=KEPS*(Lsh(i)*EOJ+KPHI)/KSD
        enddo
     enddo
     do j=1,nphicells   ! Fill in DGCPM potentials
        SJ=sind(vphicells(j))
        CJ=cosd(vphicells(j))
        EOJ=RE*(KEY*sj+KEX*cj)
        Ro=KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
        if (LAMGAM.eq.0) Ro=0.
        do i=1,nthetacells
           RoRg=(Ro/vlzcells(i))**KGAM
           KSD=1.+RoRg
           potdgcpm(i,j)=KEPS*(vlzcells(i)*EOJ+KPHI)/KSD
        enddo
     enddo

     !\
     ! Unshielded Volland-Stern field
     !/
  else if (ABASE(IA+1).eq.4) then
     KEPS=Ppc0/DLMAG
     do J=1,JO   !  Fill in drift values
        do i=1,io
           do l =1, lo
              do k = 1, ko
                 VrConv(i,j,k,l)=-KEPS*cos(PHI(J))*(LZ(I)+0.5*DL1)**3   &
                      *DT/DL1*(RE/ME)
                 select case(TypeBCalc)
                 case('analytic')
                    VR(i,j,k,l) = VrConv(i,j,k,l)
                 case('numeric')
                    VR(i,j,k,l) = VrConv(i,j,k,l) +ISS* VR_IIII(i,j,k,l)*DT/DL1
                 end select
              end do
           end do
           P1(I,J)=KEPS*LZ(I)**2*sin(PHI(J)+0.5*DPHI)*   &
                DT/DPHI*(RE/ME)
        end do
        do i=1,io+3
           BASEPOT(i,j)=KEPS*Lsh(i)*sin(phi(j))
        end do
     end do
     do J=1,nphicells   ! Fill in DGCPM potentials
        do I=1,nthetacells
           potdgcpm(i,j)=KEPS*vlzcells(I)*SIND(vphicells(J))
        enddo
     enddo

     !\
     ! Base field is specified during the self-consistent calculation
     !/

  else if (ABASE(IA+1).ge.5) then
     do J = 1, JO   !  Zero out drift values
        do i = 1, io
           do l =1, lo
              do k =1, ko
                 VrConv(i,j,k,l) = 0.0
                 select case(TypeBCalc)
                 case('analytic')
                    VR(i,j,k,l) = VrConv(i,j,k,l)
                 case('numeric')
                    VR(i,j,k,l) = VrConv(i,j,k,l) + ISS*VR_IIII(i,j,k,l)*DT/DL1
                 end select
              end do
           end do
           P1(I,J) = 0.0
        end do
     end do
     do J = 1, nphicells   ! Zero out DGCPM potentials
        do I = 1, nthetacells
           potdgcpm(i,j) = 0.0
        enddo
     enddo

     !\
     ! Read in potentials from a file
     !/

  else if (ABASE(IA+1).eq.10) then

     iUnitPot = io_unit_new()

     if (T.eq.TIME) then
        TP2=TIME-1.
        TP1=TP2
        open(UNIT=iUnitPot,FILE=NameInputDirectory//trim(NameRun)//'_pot.in',status='old')
        do i=1,5                   ! lines of header material
           read (iUnitPot,*) HEADER
        end do
        read (iUnitPot,*) JOpot
        read (iUnitPot,*) (MLTpot(j),j=1,JOpot)
        read (iUnitPot,*) header
        read (iUnitPot,*) IOpot
        read (iUnitPot,*) (LZpot(I),i=1,IOpot)
        do j=JOpot+1,NT
           MLTpot(j)=MLTpot(j-1)+(MLTpot(JOpot)-MLTpot(JOpot-1))
        enddo
        do i=IOpot+1,3*NR
           LZpot(i)=LZpot(i-1)+(LZpot(IOpot)-LZpot(IOpot-1))
        enddo
        FPOT2(1:3*NR,1:NT)=0.
     end if
     if (TP2.lt.T) then			! Read in potentials
        do while (TP2.le.T)		! Best if final TP2 > final T
           FPOT1(1:IOpot,1:JOpot)=FPOT2(1:IOpot,1:JOpot)
           read (iUnitPot,*,IOSTAT=L) TP2
           do j=1,JOpot
              read (iUnitPot,*,IOSTAT=L) (data(I),I=1,IOpot)
              FPOT2(1:IOpot,J)=data(1:IOpot)  ! Potential in Volts
           end do
           FPOT2(1:IOpot,JOpot+1)=FPOT2(1:IOpot,1)
           if (L.lt.0) TP2=TIME+2*DT*(NSTEP+1)
           if (T.eq.TIME) then		! In case T2>T already
              TP1=TIME
              FPOT1(1:IOpot,1:JOpot+1)=FPOT2(1:IOpot,1:JOpot+1)
           end if
        end do
     end if
     close(iUnitPot)

     FAC=(T-TP1)/(TP2-TP1)			! Linearly interpolate
     Fmax=0.
     Fmin=0.
     do j=1,JOpot+1
        FPOT12(1:IOpot,J)=FAC*FPOT2(1:IOpot,J)+(1.-FAC)*FPOT1(1:IOpot,J)
        BASEPOT(1:IO+3,J)=FPOT12(1:IO+3,J)
        do i=1,IOpot
           if (FPOT12(i,j).gt.Fmax) Fmax=FPOT12(i,j)
           if (FPOT12(i,j).lt.Fmin) Fmin=FPOT12(i,j)
        enddo
     enddo

     call write_prefix; write(iUnitStdOut,*) 'FPOT12 Max, Min:',Fmin,Fmax,Fmax-Fmin
     !\
     ! transform FPOT to equatorial plane drifts
     !/
     do j=1,jo
        SJ=sin(PHI(J))
        CJ=cos(PHI(J))
        jj=j+1
        !	    IF (j.eq.jo) jj=1
        do i=1,io   ! Note this is shifted from FPOT12's I grid by 1
           RR=LZ(I)+0.5*DL1
           EPP=-(FPOT12(i+1,jj)-FPOT12(i+1,j))/RR/DPHI
           do l = 1, lo
              do k = 1, ko
                 VrConv(i,j,k,l) = VrConv(i,j,k,l)+EPP*DT/DL1*RR**3*RE/ME
                 select case(TypeBCalc)
                 case('analytic')
                    VR(i,j,k,l) = VrConv(i,j,k,l)
                 case('numeric')
                    VR(i,j,k,l) = VrConv(i,j,k,l) + ISS*VR_IIII(i,j,k,l)*DT/DL1
                 end select
              end do
           end do
        end do
        SJ=sin(PHI(J)+.5*DPHI)
        CJ=cos(PHI(J)+.5*DPHI)
        do i=1,io
           RR=LZ(I)
           EPR=-(FPOT12(i+2,j)-FPOT12(i+1,j))/DL1
           P1(i,j)=P1(i,j)-EPR*DT/DPHI*RR**2*RE/ME
        end do
     end do

     !\
     ! potdgcpm filled in from FPOT12
     !/

     Fmax=0.
     Fmin=0.

     do j = 1, nphicells   !   Fill in DGCPM potentials
        do i = 1, nthetacells
           call heidi_lintp2(LZpot,MLT,FPOT12,3*NR,NT,   &
                vlzcells(i),vmltcells(J),FAC,IER)
           potdgcpm(i,j)=FAC
        enddo
        do i = 1, nthetacells
           if (potdgcpm(i,j).gt.Fmax) Fmax=potdgcpm(i,j)
           if (potdgcpm(i,j).lt.Fmin) Fmin=potdgcpm(i,j)
        enddo
     enddo

     call write_prefix; write(iUnitStdOut,*)  'potdgcpm Max, Min:',Fmin,Fmax,Fmax-Fmin

     !\
     ! Calculate potentials using an ionospheric model (AMIE, Weimer, etc.)
     !/

  else if (ABASE(IA+1).ge.11) then

     ! Transfer date and time to new variables
     eyear=year
     eday=day
     edayplus=aint(t/86400.)
     eday=eday+edayplus
     univ_time=t-edayplus*86400.
     ehour=aint(univ_time/3600.)
     univ_time=univ_time-ehour*3600.
     eminute=aint(univ_time/60.)
     univ_time=univ_time-eminute*60.
     esecond=univ_time

     ! Transfer solar wind values to new variables
     eby=BYSW
     ebz=BZSW
     evsw=USW
     ! Transfer grids to new variables
     ! Note that we have 2 sets of grids, for RAM and for DGCPM
     nemlts1=NT
     nelats1=NR+3
     do i=1,nelats1
        emlts1(1:nemlts1,i)=MLT(1:nemlts1)
     end do
     do j=1,nemlts1
        elats1(j,1:nelats1)=lats(1:nelats1)
     end do
     nemlts2=nphicells
     nelats2=nthetacells
     do i=1,nelats2
        emlts2(1:nemlts2,i)=vmltcells(1:nemlts2)
     end do
     do j=1,nemlts2
        elats2(j,1:nelats2)=90.-vthetacells(1:nelats2)
     end do
     !\
     ! Transfer the potential model choice to a new variable
     ! nechoice = 1 means AMIE potentials
     !          = 2 means Weimer-96 potentials
     !          = 3 means Foster potentials
     !          = ?? we have to add in more to the DATA defs above
     !/

     nechoice=ABASE(IA+1) - 10

     ! Aaron: insert your call here, using the "e" variables

     call get_potential_heidi(nechoice,eby,ebz,evsw,eyear,eday,ehour,   &
          eminute,esecond,elats1,emlts1,epots1,elats2,emlts2,epots2)

     ! Use epots1 to update the advection coefficients:
     do j=1,jo  ! transform pots1 to equatorial plane drifts
        jj=j+1
        if (j.eq.jo) jj=1
        BASEPOT(1:io+3,j)=epots1(j,1:io+3)
        do i=1,io   ! Note this is shifted from epots1's I grid
           RR=LZ(I)+0.5*DL1
           EPP=-(epots1(jj,i+1)-epots1(j,i+1))/RR/DPHI
           do l = 1, lo
              do k = 1, ko
                 VrConv(i,j,k,l) = EPP*DT/DL1*RR**3*RE/ME
                 select case(TypeBCalc)
                 case('analytic')
                    VR(i,j,k,l) = VrConv(i,j,k,l)
                 case('numeric')
                    VR(i,j,k,l) = VrConv(i,j,k,l) + ISS*VR_IIII(i,j,k,l)*DT/DL1
                 end select
              end do
           end do
        end do
        do i=1,io
           RR=LZ(I)
           EPR=-(epots1(j,i+2)-epots1(j,i+1))/DL1
           P1(i,j)=-EPR*DT/DPHI*RR**2*RE/ME
        end do
     end do
     ! Write pots2 into potdgcpm (for Dan Ober's code)
     do j=1,nphicells
        potdgcpm(1:nthetacells,j)=epots2(j,1:nthetacells)
     end do

     !\
     ! Done with base field calculations
     !/

  end if

  !\
  ! Calculate Burke-Wygant penetration electric field
  !/

  if (AEPEN(IA+1).eq.1) then
     if (ABASE(IA+1).eq.1) then
        PHIP=sqrt(PCO*RE*A*(1.5*DLMAG)**(LAMGAM-1))
     else if (ABASE(IA+1).ge.2) then
        PHIP=sqrt(PCO*RE*Eom*KEPS)
     end if
     PHIP=AMAX1(0.,1.4*PHIP-3200.)  ! From Burke's analysis, fig12
     PHIP=-PHIP*FACP                ! To have the Burke value at MLT=18
     do J=1,JO                      ! Fill in radial drift values
        CJ=cos(PHI(J))
        SJ=sin(PHI(J))
        PJ=PHI(J)-PHIPOFF
        CPJ=abs(PJ)
        SPJ=1.
        if (exp(-(CPJ/DP1)**3) .lt. exp(-((2.*cPi-CPJ)/DP2)**3)) then
           CPJ=2.* cPi-CPJ
           if (PJ.gt.0.) SPJ=-1.
           DPP=DP2
        else
           if (PJ.lt.0.) SPJ=-1.
           DPP=DP1
        endif
        CPJ=CPJ/DPP
        LP=0.5+KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
        dLPdphi=-KBETA*SJ*(KS(2)+Kr*KS(4))
        sLP=sqrt(LP)
        do i=1,io
           RR=LZ(I)+0.5*DL1
           DLP=10.
           if (RR.lt.LP) DLP=sLP-1.
           CRo=(RR-LP)/DLP
           RoRg=exp(-CRo**2-CPJ**3)
           KSD=3.*SPJ*CPJ**2-2.*CRo*dLPdphi/DLP*(1.+CRo/sLP)
           EPP=PHIP*RoRg/RR*KSD
           FAC=EPP*DT/DL1*RR**3*RE/ME
           do l = 1, lo
              do k = 1, ko
                 VrConv(i,j,k,l)=VrConv(i,j,k,l)+FAC
                 select case(TypeBCalc)
                 case('analytic')
                    VR(i,j,k,l) = VrConv(i,j,k,l)
                 case('numeric')
                    VR(i,j,k,l) = VrConv(i,j,k,l) + ISS*VR_IIII(i,j,k,l)*DT/DL1
                 end select
              end do
           end do
        end do
     end do
     do J=1,JO   ! Fill in azimuthal drift values
        CJ=cos(PHI(J)+.5*DPHI)
        PJ=PHI(J)-PHIPOFF+.5*DPHI
        CPJ=abs(PJ)
        if (exp(-(CPJ/DP1)**3) .lt. exp(-((2.* cPi-CPJ)/DP2)**3)) then
           CPJ=2.* cPi-CPJ
           DPP=DP2
        else
           DPP=DP1
        endif
        CPJ=CPJ/DPP
        LP=0.5+KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
        sLP=sqrt(LP)
        do i=1,io
           RR=LZ(I)
           DLP=10.
           if (RR.lt.LP) DLP=sLP-1.
           CRo=(RR-LP)/DLP
           RoRg=exp(-CRo**2-CPJ**3)
           KSD=2.*CRo/DLP
           EPR=PHIP*RoRg*KSD
           FAC=-EPR*DT/DPHI*RR**2*RE/ME
           P1(I,J)=P1(I,J)+FAC
        end do
     end do
     do j=1,nphicells   ! Fill in DGCPM potentials
        CJ=cosd(vphicells(j))
        PJ=vphicells(j)*cPi/180.+phipoff
        CPJ=abs(PJ)
        if (exp(-(CPJ/DP1)**3) .lt. exp(-((2.* cPi-CPJ)/DP2)**3)) then
           CPJ=2.* cPi-CPJ
           DPP=DP2
        else
           DPP=DP1
        endif
        CPJ=CPJ/DPP
        LP=0.5+KBETA*(KS(1)+KS(2)*CJ+Kr*(KS(3)+KS(4)*CJ))
        sLP=sqrt(LP)
        do i=1,nthetacells
           RR=vlzcells(i)
           DLP=10.
           if (RR.lt.LP) DLP=sLP-1.
           CRo=(RR-LP)/DLP
           RoRg=exp(-CRo**2-CPJ**3)
           potdgcpm(i,j)=potdgcpm(i,j)+PHIP*RoRg
        enddo
     enddo
  end if
51 format(f6.3,f6.2,1p,10e11.3)

  !\
  ! Calculate self-consistent penetration electric field
  !/

  if (AEPEN(IA+1).ge.2) then
     call write_prefix; write(iUnitStdOut,*)  '...Calling PRESSURES'
     call PRESSURES    ! Updates the bulk ring current parameters
     call write_prefix; write(iUnitStdOut,*) '...Calling CURRENTCALC' 
     call CURRENTCALC  ! Finds the perpendicular "ring current"
     ! and the field-aligned closure currents

    write(*,*) 'IA ===========', IA

!!$STOP

     Irfac=Ir
     Latfac(1:Ir)=Lats(1:Ir)

     Ilfac=JO
     Lonfac(1:JO)=MLT(1:JO)


     bc_choice=0       ! High-lat boundary has pot=0.
     if (ABASE(IA+1).eq.5) bc_choice=1  ! W-96 high-lat BC
     if (ABASE(IA+1).eq.6) bc_choice=2  ! W-96 high-lat BC
     if (bc_choice.ne.0 .and. AEPEN(IA+1).eq.3) bc_choice=-bc_choice
     !             !! W-96/AMIE everywhere, done in S-C potential call

     ! The Jfac adjustment that's commented out is if you only do 1
     ! species and want to proportionally increase the Jfac value
     ! to reflect the current from all species
     !	  NY(2)=0.34*EXP(0.054*KP)	! From Young et al 1982
     !	  NY(3)=5.1E-3*EXP(6.6E-3*F107)
     !	  NY(4)=0.011*EXP(0.24*KP+0.011*F107)
     !	  FAC=(NY(2)+NY(3)+NY(4))/NY(S)
     !	  DO j=1,JO
     !	    Jfac(1:Ir,j)=Jion1(1:Ir,j)*FAC
     !	  end do
     ! Here the Jfac from the various species are summed together

     if (nProc > 1) then
        call MPI_BARRIER(iComm,iError)
        Jfac = 0.0
        do i = 1, nProc
           if (i-1.eq.iProc) then
              jfac_temp(:,:) = Jion1(:,:,nParallelSpecies(i))
           endif
           call MPI_Bcast(jfac_temp,(NR+3)*NT,MPI_Real,   &
                i-1,iComm,iError)
           Jion1(:,:,nParallelSpecies(i)) = jfac_temp(:,:)
           JFAC(:,:)=Jfac(:,:)+Jreducer*Jion1(:,:,nParallelSpecies(i))
        enddo
     else
        do J  = 1, JO
           Jfac(1:Ir,J)=0.
           do S = 1, NS   
              if (SCALC(S).eq.1) JFAC(1:Ir,J)=Jfac(1:Ir,j)+Jreducer*Jion1(1:Ir,j,s)
           end do
        end do
     endif

     if (iProc==0) then
        call write_prefix; write(iUnitStdOut,*) 'FACs being sent to potential solver' 
     end if
50   format(F5.1,1P,25E12.4)

     !\
     ! The next few lines are for doubling convection during sawtooth "snaps"
     !/

     evsw=USW
     do ii=1,9
        tdiff = t - t_sawtooth(ii)
        if (tdiff.ge.0. .and. tdiff.le.dt_saw) evsw=2.*evsw
     end do
     ! End the sawtooth snap block
     if (iProc==0) then
        call write_prefix; write(iUnitStdOut,*) '...Calling EPENCALC'
        call EPENCALC(t,f107,bc_choice,BYSW,BZSW,evsw) 
        ! Aaron Ridley's solver 
        ! for the potential
        ! from the field-aligned currents
        ! Note: By and Bz in nT, Usw in m/s

        call write_prefix; write(iUnitStdOut,*) &
             'Potentials returned from the solver'
     end if

     call MPI_Bcast(fpot,(NR+3)*NT,MPI_Real,0,iComm,iError)

     !\
     ! While EPENCALC is commented out, use these lines to zero out FPOT
     !	do j=1,jo
     !	 do i=1,ir
     !	  FPOT(i,j)=0.
     !	 enddo
     !	end do  ! End FPOT zeroing
     !/


     if (AEPEN(IA+1).gt.1) then  !! >2 Prohibits inclusion of Epen
        do j=1,jo  ! transform FPOT to equatorial plane drifts
           SJ=sin(PHI(J))
           CJ=cos(PHI(J))
           jj=j+1
           if (j.eq.jo) jj=1
           do i=1,io   ! Note this is shifted from FPOT's I grid
              RR=LZ(I)+0.5*DL1
              EPP=-(FPOT(i+1,jj)-FPOT(i+1,j))/RR/DPHI
              do l = 1, lo
                 do k = 1, ko
                    VrConv(i,j,k,l) = VrConv(i,j,k,l)+EPP*DT/DL1*RR**3*RE/ME
                    select case(TypeBCalc)
                    case('analytic')
                       VR(i,j,k,l) = VrConv(i,j,k,l)
                    case('numeric')
                       VR(i,j,k,l) = VrConv(i,j,k,l) + ISS*VR_IIII(i,j,k,l)*DT/DL1

                    end select
                 end do
              end do
           end do
           SJ=sin(PHI(J)+.5*DPHI)
           CJ=cos(PHI(J)+.5*DPHI)
           do i=1,io
              RR=LZ(I)
              EPR=-(FPOT(i+2,j)-FPOT(i+1,j))/DL1
              P1(i,j)=P1(i,j)-EPR*DT/DPHI*RR**2*RE/ME
           end do
        end do

        !\
        !potdgcpm filled in write_ring_current from epencalc
        !	  do j=1,nphicells   !   Fill in DGCPM potentials
        !	   do i=1,nthetacells
        !	    CALL heidi_lintp2(LZ,MLT,FPOT,NR,NT,
        !     &	      vlzcells(I),vmltcells(J),FAC,IER)
        !	    potdgcpm(i,j)=potdgcpm(i,j)+FAC
        !	   enddo
        !	  enddo
        !/

     end if  !! Prohibits inclusion of Epen

  end if

  call write_prefix; write(iUnitStdOut,*)'Done with MAGCONV'


  write(*,*) 'ISS*DT/DL1', ISS*DT/DL1
  

!!$  open(unit=3, file='DriftVsPA.dat')
!!$  do i =14,14
!!$     do j = 1,nT
!!$        do k =37, 37
!!$           do l =1, nPa
!!$              write(3,*) VPhi_IIII(i,j,k,l), P2(i,j,k,l), P1(i,j), Vr(i,j,k,l),&
!!$                   VR_IIII(i,j,k,l), VrConv(i,j,k,l), VR_IIII(i,j,k,l)+VrConv(i,j,k,l),&
!!$                   i, j, k, l, LZ(i),  Phi(j),  acos(mu(l))
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  close(3)
!!$
!!$  open(unit=3, file='DriftVsL.dat')
!!$
!!$  do j =1, nT
!!$     do k =37, 37
!!$        do l =2, 2
!!$           do i = 1, nR
!!$              write(3,*) VPhi_IIII(i,j,k,l), P2(i,j,k,l), P1(i,j), Vr(i,j,k,l),&
!!$                   VR_IIII(i,j,k,l), VrConv(i,j,k,l), VR_IIII(i,j,k,l)+VrConv(i,j,k,l),&
!!$                   i, j, k, l, LZ(i),  Phi(j),  acos(mu(l))
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  close(3)
!!$
!!$
!!$
!!$  open(unit=3, file='DriftVsPhi.dat')
!!$  write(3,*) 'drifts'
!!$  write(3,*) 'VPhi_IIII       P2       P1   VR  VrConv   VR_IIII+VrConv     i   j   k   l   lz      Phi pa'
!!$
!!$  do i =14, 14
!!$     do j =1, nT
!!$        do k =37, 37
!!$           do l =2, 2
!!$              write(3,*) VPhi_IIII(i,j,k,l), P2(i,j,k,l), P1(i,j),Vr(i,j,k,l),&
!!$                   VR_IIII(i,j,k,l), VrConv(i,j,k,l), VR_IIII(i,j,k,l)+VrConv(i,j,k,l),&
!!$                   i, j, k, l, LZ(i),  Phi(j),  acos(mu(l))
!!$
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$  close(3)
!!$

end subroutine MAGCONV

!===================================================================

