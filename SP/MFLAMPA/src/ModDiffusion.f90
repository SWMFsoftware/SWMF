!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModDiffusion
  !Solve the diffusion term in the Parker equation
  !Version from FLAMPA, I.Sokolov and I.Roussev
  !Updated:  I.Sokolov, Dec.17, 2017
  !Revision history:
  implicit none
  PRIVATE
  public:: advance_diffusion
contains
  !================================================================
  ! This routine solves three-diagonal system of equations: 
  !  ||m_1 u_1  0....        || ||w_1|| ||r_1||
  !  ||l_2 m_2 u_2...        || ||w_2|| ||r_2||
  !  || 0  l_3 m_3 u_3       ||.||w_3||=||r_3||
  !  ||...                   || ||...|| ||...||
  !  ||.............0 l_n m_n|| ||w_n|| ||r_n||
  ! From: Numerical Recipes, Chapter 2.6, p.40.                               
  subroutine tridiag(n, L_I, M_I, U_I, R_I, W_I)
    !\
    ! input parameters
    !/
    integer,            intent(in):: n
    real, dimension(n), intent(in):: L_I, M_I ,U_I ,R_I
    !\
    ! Output parameters
    !/
    real, intent(out):: W_I(n)
    !\
    ! Misc
    !/
    integer:: j
    real:: Aux,Aux_I(2:n)
    !-------------------!
    if (M_I(1)==0.0)&
         call CON_stop(' Error in tridiag: M_I(1)=0')
    Aux = M_I(1)
    W_I(1) = R_I(1)/Aux
    do j=2,n
       Aux_I(j) = U_I(j-1)/Aux
       Aux = M_I(j)-L_I(j)*Aux_I(j)
       if (Aux.eq.0.0) then
          write(*,*)'M_I(j), L_I(j), Aux_I(j) = ',&
               M_I(j),L_I(j),Aux_I(j)
          write(*,*)'  For j=',j
          call CON_stop('Tridiag failed')
       end if
       W_I(j) = (R_I(j)-L_I(j)*W_I(j-1))/Aux
    end do
    do j=n-1,1,-1
       W_I(j) = W_I(j)-Aux_I(j+1)*W_I(j+1)
    end do
  end subroutine tridiag
  !=========================================================================
  ! This routine solves the diffusion equation:
  !         f_t-D_outer(D_inner*f_x)_x=0,
  ! with zero Neumann boundary condition. The solution is advanced in time
  ! using fully implicit scheme.
  subroutine advance_diffusion(Dt,n,D_I,F_I,DOuter_I,DInner_I)
    use ModNumConst, ONLY: cTiny

    real,   intent(in   ):: Dt     !Time step                       
    integer,intent(in   ):: n      !Number of meshes along the x-coordinate
    real,   intent(in   ):: D_I(n) !Distance to the next mesh
    real,   intent(inout):: F_I(n) !In:sol.to be advanced; Out:advanced sol
    !Laplace multiplier and diffusion coefficient.
    real,   intent(in   ):: DOuter_I(n), DInner_I(n)
   
    !Mesh spacing and face spacing.
    real                 :: DsMesh_I(2:n), DsFace_I(2:n-1)
    !Main, upper, and lower diagonals.
    real, dimension(n)   :: Main_I,Upper_I,Lower_I, R_I 
    integer:: i
    real:: Aux1,Aux2
    !-----------------------------------------------------------------
    !D_I(i)      is the distance between meshes i   and i+1
    !DsMesh_I(i) is the distance between meshes i-1 and i 
    do i=2,n
       DsMesh_I(i) = max(D_I(i-1),cTiny)
    end do
    !Distance between the faces is enumerated with the cell index
    do i=2,n-1
       DsFace_I(i) = max(0.5*(D_I(i) + D_I(i-1)),cTiny)
    end do
    !\
    ! f^(n+1)_i-Dt*DOuter_I/DsFace_I*(&
    !     DInner_(i+1/2)*(f^(n+1)_(i+1)-f^(n+1)_i)/DsMesh_(i+1)-&
    !     DInner_(i-1/2)*(f^(n+1)_i -f^(n+1)_(i-1)/DsMesh_(i ))=f^n_i
    !/
    Main_I = 1.0
    !\
    ! For i=1:
    !/
    Aux1 = Dt*DOuter_I(1)*0.50*(DInner_I(1)+DInner_I(2))/&
         DsMesh_I(2)**2
    Main_I( 1) = Main_I(1)+Aux1
    Upper_I(1) = -Aux1
    !\
    ! For i=2,n-1:
    !/
    do i=2,n-1
       Aux1 = Dt*DOuter_I(i)*0.50*(DInner_I(i  ) + DInner_I(i+1))/&
            (DsMesh_I(i+1)*DsFace_I(i))
       Aux2 = Dt*DOuter_I(i)*0.50*(DInner_I(i-1) + DInner_I(i  ))/&
            (DsMesh_I(i  )*DsFace_I(i))
       Main_I(i)  = Main_I(i) + Aux1 + Aux2
       Upper_I(i) = -Aux1
       Lower_I(i) = -Aux2
    end do
    !\
    ! For i=n:
    !/
    Aux2 = Dt*DOuter_I(n)*0.50*(DInner_I(n-1) + DInner_I(n))/&
         DsMesh_I(n)**2
    Main_I( n) = Main_I(n) + Aux2
    Lower_I(n) = -Aux2
    !\
    ! Update the solution from f^(n) to f^(n+1):
    !/
    R_I = F_I
    call tridiag(n,Lower_I,Main_I,Upper_I,R_I,F_I)
  end subroutine advance_diffusion
end module SP_ModDiffusion
