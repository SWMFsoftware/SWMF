subroutine get_E_mu_dot
  
  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts
  use ModHeidiInput, ONLY: TypeBField
  
  implicit none
  
  integer :: i,j,k,l!,is,iss,ier
  real    :: MUBOUN,MULC
  real    :: gpa

  !real    :: eps,dln,sqm,qe,gama,cco,ccd,edrco,x,xd,g,cci,cdi,edri,ccde,ccdi
  !real    :: bane,badif,c,gpa,cce,cde,edre
  !real    :: erf
  !real    :: COULDE(NE,NPA),COULDI(NE,NPA),AFIR,ASEC
  !real    :: VF(NSTH),RA(NSTH),MUBOUN,MULC,TMAS(NSTH),TM1(NSTH)
  
  real, dimension(nR,nT,nE,nPA) :: dRdt_IIII ,dPhiDt_IIII, InvRdRdt_IIII,dEdt_IIII
  real, dimension(nR,nT,nPA)    :: dMudt_III 

 ! external:: erf

!-----------------------------------------------------------------------------

  open (unit = 2, file = "EDOT.dat")
  write (2,*)'Values '
  write (2,*)'i,j,k,l,EDOT' 
  
  
  if (TypeBField == 'analytic') then 
     
     do I=1,IO
        do J=1,JO
           do L=1,UPA(I)
              MUBOUN=MU(L)+0.5*WMU(L)           ! MU at boundary of grid
              MUDOT(I,J,L)=(1.-MUBOUN**2)*(0.5*(FUNI(L+1,I,J)+FUNI(L,I,J)))/LZ(I)  &
                   /MUBOUN/4./(0.5*(FUNT(L+1,I,J)+FUNT(L,I,J)))*DL1/DMU(L)
              GPA=1.-FUNI(L,I,J)/6./FUNT(L,I,J)
              do K=1,KO
                 EDOT(I,J,K,L)=-3.*EBND(K)/LZ(I)*GPA*DL1/DE(K)
                 write(2,*) 'i,j,k,l,EDOT',i,j,k,l,EBND(k),EDOT(I,J,K,L)
              end do	! K loop
           end do 	! L loop
           MULC=MU(UPA(I))+0.5*WMU(UPA(I))
           do L=UPA(I)+1,LO-1
              MUBOUN=MU(L)+0.5*WMU(L)
              if (l== LO-1) MU(L+1) = MU(L)
              MUDOT(I,J,L)=(1.-MUBOUN**2)*(0.5*(FUNI(L+1,I,J)+FUNI(L,I,J)))/LZ(I)  &
                   /MUBOUN/4./(0.5*( FUNT(L+1,I,J)+FUNT(L,I,J)))*DL1/DMU(L)
              do K=1,KO
                 EDOT(I,J,K,L)=-3.*EBND(K)/LZ(I)*GPA*DL1/DE(K)
              end do	! K loop
           end do	! L loop
           MUDOT(I,J,LO)=0.
           do K=1,KO
              EDOT(I,J,K,LO)=-3.*EBND(K)/LZ(I)*GPA*DL1/DE(K)
           end do	! K loop
        end do 	! J loop
     end do	! I loop
 

     write(*,*) 'EDOT_analytic', EDOT(20,24,42,67)
    
  end if

  
  if (TypeBField == 'numeric') then 
     call get_coef(dRdt_IIII ,dPhiDt_IIII, InvRdRdt_IIII,dEdt_IIII,dMudt_III)

      do I=1,IO
        do J=1,JO
           do L=1,UPA(I)
              MUBOUN=MU(L)+0.5*WMU(L)           ! MU at boundary of grid
              MUDOT(I,J,L) = -((1.-MUBOUN**2)/MUBOUN)*dMudt_III(i,j,l)*DL1/DMU(L)
              do K=1,KO
                 EDOT(I,J,K,L)= dEdt_IIII(i,j,k,l)*DL1/DE(K)
                 write(2,*) 'i,j,k,l,EDOT',i,j,k,l,EBND(k),EDOT(I,J,K,L)
              end do	! K loop
           end do 	! L loop
           MULC=MU(UPA(I))+0.5*WMU(UPA(I))
           do L=UPA(I)+1,LO-1
              MUBOUN=MU(L)+0.5*WMU(L)
              if (l== LO-1) MU(L+1) = MU(L)
              MUDOT(I,J,L)= -((1.-MUBOUN**2)/MUBOUN)*dMudt_III(i,j,l)*DL1/DMU(L)
              do K=1,KO
                 EDOT(I,J,K,L)= dEdt_IIII(i,j,k,l)*DL1/DE(K)
              end do	! K loop
           end do	! L loop
           MUDOT(I,J,LO)=0.
           do K=1,KO
              EDOT(I,J,K,L)= dEdt_IIII(i,j,k,l)*DL1/DE(K)
           
           end do	! K loop
        end do 	! J loop
     end do	! I loop
 
     write(*,*) 'EDOT_numeric', EDOT(20,24,42,67)   
     
  end if

  close(2)
stop

end subroutine get_E_mu_dot
