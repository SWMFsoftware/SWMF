!============================================================================!
! This routine solves three-diagonal system of equations:                    !
!  ||m_1 u_1  0....        || ||w_1|| ||r_1||                                !
!  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||                                !
!  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||                                !
!  ||...                   || ||...|| ||...||                                !
!  ||.............0 l_n m_n|| ||w_n|| ||r_n||                                !
! From: Numerical Recipes, Chapter 2.6, p.40.                                !
!============================================================================!
subroutine tridag(n,L_I,M_I,U_I,R_I,W_I)
  use SP_ModMain,ONLY:iStdOut,prefix
  implicit none
  !--------------------------------------------------------------------------!
  integer, intent(in):: n
  real, intent(in):: L_I(n),M_I(n),U_I(n),R_I(n)
  real, intent(out):: W_I(n)
  !--------------------------------------------------------------------------!
  integer:: j
  real:: Aux,Aux_I(2:n)
  !--------------------------------------------------------------------------!
  if (M_I(1).eq.0.0) then
       write(iStdOut,*)'Error in tridag: M_I(1)=0'
       stop
    end if
  Aux = M_I(1)
  W_I(1) = R_I(1)/Aux
  do j=2,n
     Aux_I(j) = U_I(j-1)/Aux
     Aux = M_I(j)-L_I(j)*Aux_I(j)
     if (Aux.eq.0.0) then
        write(iStdout,*)prefix,'M_I(j), L_I(j), Aux_I(j) = ',M_I(j),L_I(j),Aux_I(j)
        write(iStdOut,*)'Tridag failed: for j=',j
        stop
     end if
     W_I(j) = (R_I(j)-L_I(j)*W_I(j-1))/Aux
  end do
  do j=n-1,1,-1
     W_I(j) = W_I(j)-Aux_I(j+1)*W_I(j+1)
  end do
  !------------------------------------ DONE --------------------------------!
end subroutine tridag
!============================================================================!
! This routine solves the diffusion equation:                                !
!         f_t-D_outer(D_inner*f_x)_x=0,                                      !
! with zero Neumann boundary condition. The solution is advanced in time     !
! using fully implicit scheme.                                               !
!============================================================================!
subroutine advance_diffusion(Dt,n,X_DI,F_I,DOuter_I,DInner_I)
  use ModNumConst
  implicit none
  !--------------------------------------------------------------------------!
  real,intent(in):: Dt        !Time step.                                    !
  integer,intent(in):: n      !Number of meshes along the x-coordinate.      !
  real,dimension(3,n),intent(in):: &
       X_DI                   !Three Cartesian coordinates of each mesh.     !
  real,dimension(n),intent(inout):: &
       F_I                    !In: sol. to be advanced; Out: advanced sol.   !
  real,dimension(n),intent(in):: &
       DOuter_I,DInner_I      !Laplace multiplier and diffusion coefficient. !
  !--------------------------------------------------------------------------1
  real:: DsMesh_I(2:n), &
       DsFace_I(2:n-1)        !Mesh spacing and face spacing.                !
  real, dimension(n):: &
       Main_I,Upper_I,Lower_I !Main, upper, and lower diagonals.             !
  integer:: i
  real:: Aux1,Aux2
  !--------------------------------------------------------------------------!
  do i=2,n
     DsMesh_I(i) = max(sqrt(sum((X_DI(:,  i)-X_DI(:,i-1))**2)),cTiny)
  end do
  do i=2,n-1
     DsFace_I(i) = max(sqrt(sum((X_DI(:,i+1)-X_DI(:,i-1))**2))*cHalf,cTiny)
  end do
  !--------------------------------------------------------------------------!
  ! f^(n+1)_i-Dt*DOuter_I/DsFace_I*(&                                        !
  !          DInner_(i+1/2)*(f^(n+1)_(i+1)-f^(n+1)_i)/DsMesh_(i+1/2)-&       !
  !          DInner_(i-1/2)*(f^(n+1)_i -f^(n+1)_(i-1)/DsMesh_(i-1/2))=f^n_i  !
  !--------------------------------------------------------------------------!
  Main_I = cOne
  !\
  ! For i=1:
  !/
  Aux1 = Dt*DOuter_I(1)*cHalf*(DInner_I(1)+DInner_I(2))/&
       DsMesh_I(2)**2
  Main_I(1) = Main_I(1)+Aux1
  Upper_I(1) = -Aux1
  !\
  ! For i=2,n-1:
  !/
  do i=2,n-1
     Aux1 = Dt*DOuter_I(i)*cHalf*(DInner_I(i  )+DInner_I(i+1))/&
          (DsMesh_I(i+1)*DsFace_I(i))
     Aux2 = Dt*DOuter_I(i)*cHalf*(DInner_I(i-1)+DInner_I(i  ))/&
          (DsMesh_I(i  )*DsFace_I(i))
     Main_I(i) = Main_I(i)+Aux1+Aux2
     Upper_I(i) = -Aux1
     Lower_I(i) = -Aux2
  end do
  !\
  ! For i=n:
  !/
  Aux2 = Dt*DOuter_I(n)*cHalf*(DInner_I(n-1)+DInner_I(n))/&
       DsMesh_I(n)**2
  Main_I(n) = Main_I(n)+Aux2
  Lower_I(n) = -Aux2
  !\
  ! Update the solution from f^(n) to f^(n+1):
  !/
  call tridag(n,Lower_I,Main_I,Upper_I,F_I,F_I)
  !------------------------------------ DONE --------------------------------!
end subroutine advance_diffusion
!============================================================================!
subroutine enhance_diffusion(Din_I,Dout_I)
  use SP_ModMain
  implicit none
  real,dimension(nX),intent(inout):: Din_I  
  real,dimension(nX),intent(in):: Dout_I
  Din_I = max(Din_I,DiffCoeffMin/Dout_I) 
end subroutine enhance_diffusion
!============================================================================!
