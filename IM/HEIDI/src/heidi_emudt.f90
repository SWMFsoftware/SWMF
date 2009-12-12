
subroutine get_E_mu_dot
  
  use ModHeidiSize
  use ModHeidiIO
  use ModHeidiMain
  use ModHeidiDrifts
  use ModHeidiInput, ONLY: TypeBField
  
  implicit none
  
  
  real,dimension(nR,nT,nE,nPA)  :: dMudt_III 
  real,dimension(nR,nT,nE,nPa)  :: dEdt_IIII,VPhi_IIII,VR_IIII
  real                          :: MUBOUN,MULC
  real                          :: gpa
  integer                       :: i,j,k,l
!-----------------------------------------------------------------------------

  if (TypeBField == 'analytic') then 
     do I=1,IO
        do J=1,JO
           do L=1,UPA(I)
              do K=1,KO
                 MUBOUN=MU(L)+0.5*WMU(L)           ! MU at boundary of grid
                 MUDOT(I,J,K,L)=(1.-MUBOUN**2)*(0.5*(FUNI(L+1,I,J)+FUNI(L,I,J)))/LZ(I)  &
                      /MUBOUN/4./(0.5*(FUNT(L+1,I,J)+FUNT(L,I,J)))*DL1/DMU(L) * VR(i,j,k,l)

                 GPA=1.-FUNI(L,I,J)/6./FUNT(L,I,J)
                 EDOT(I,J,K,L)=-3.*EBND(K)/LZ(I)*GPA*DL1/DE(K)*VR(i,j,k,l)
              end do	! K loop
           end do 	! L loop
           MULC=MU(UPA(I))+0.5*WMU(UPA(I))
           do L=UPA(I)+1,LO-1
              do K=1,KO
                 MUBOUN=MU(L)+0.5*WMU(L)
                 if (l== LO-1) MU(L+1) = MU(L)
                 MUDOT(I,J,K,L)=(1.-MUBOUN**2)*(0.5*(FUNI(L+1,I,J)+FUNI(L,I,J)))/LZ(I)  &
                      /MUBOUN/4./(0.5*( FUNT(L+1,I,J)+FUNT(L,I,J)))*DL1/DMU(L) * VR(i,j,k,l)
                 EDOT(I,J,K,L)=-3.*EBND(K)/LZ(I)*GPA*DL1/DE(K)*VR(i,j,k,l)
              end do	! K loop
           end do	! L loop
           do K = 1, KO
              MUDOT(I,J,K,LO)=0.
              EDOT(I,J,K,LO)=-3.*EBND(K)/LZ(I)*GPA*DL1/DE(K)*VR(i,j,k,l)
           end do	! K loop
        end do 	! J loop
     end do	! I loop
     
  end if


  
  if (TypeBField == 'numeric') then 
     call get_coef(dEdt_IIII,dMudt_III)

      do I=1,IO
        do J=1,JO
           do L=1,UPA(I)
              do K=1,KO
                 MUBOUN=MU(L)+0.5*WMU(L)           ! MU at boundary of grid
                 MUDOT(I,J,K,L) = ((1.-MUBOUN**2)/MUBOUN)*dMudt_III(i,j,k,l)*DL1/DMU(L)
                 EDOT(I,J,K,L)= dEdt_IIII(i,j,k,l)*DL1/DE(K)
              end do	! K loop
           end do 	! L loop
           MULC=MU(UPA(I))+0.5*WMU(UPA(I))
           do L=UPA(I)+1,LO-1
              do K=1,KO
                 MUBOUN=MU(L)+0.5*WMU(L)
                 if (l== LO-1) MU(L+1) = MU(L)
                 MUDOT(I,J,K,L)= ((1.-MUBOUN**2)/MUBOUN)*dMudt_III(i,j,k,l)*DL1/DMU(L)
                 EDOT(I,J,K,L)= dEdt_IIII(i,j,k,l)*DL1/DE(K)
              end do	! K loop
           end do	! L loop
           do K=1,KO
              MUDOT(I,J,K,LO)=0.
              EDOT(I,J,K,L)= dEdt_IIII(i,j,k,l)*DL1/DE(K)
           end do	! K loop
        end do 	! J loop
     end do	! I loop
     
  end if

end subroutine get_E_mu_dot
