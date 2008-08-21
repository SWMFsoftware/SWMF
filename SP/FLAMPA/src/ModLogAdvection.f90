module ModLogAdvection
  implicit none
  !Revision history
  !Prototype: Sokolov&Roussev, FLAMPA code, 2004
  !Present version: Sokolov& Roussev, Jan,2008

  !DESCRIPTION
!============================================================================!
! In space physics applications one often needs to solve the "log-advection" !
! equations of the kind of:                                                  !
!         f_t+A*f_lnp=0             (NC)                                     !
! or                                                                         !
!         f_t+A*(pf)_p=0 .          (C)                                      !
! The latter equation which is a CONSERVATION LAW can be also written as     !
!         f_t+A*f_lnp+ A*f=0,                                                !
! which is the total of the NONCONSERVATIVE equation plus Af                 !
! Both equations propagate the initial values of f along the lines           !
!         p=p_0*exp(At),                                                     !
! which results in a uniform expansion (A>0) or shrinking (A<0) of the       !
! distribution, keeping its shape unchanged. For the conservative equation   !
! the solution also scales accordingly, keeping unchanged the integral       !
!\int{f dp}                                                                  !
!                                                                            !
! EXAMPLES: Non-conservative formulation: first-order Fermi acceleration,    !
!           Alfven-wave turbulence evolution in the solar wind               !
!           Conservative formulation: Kolmogorov-type cascade, wave-particle !
!           interactions                                                     !
! In all these applications, A does not depend on the phase coordinate, p    !
!                                                                            !
!============================================================================!

contains
  !===========advance_log_advection==========================================!
  !DESCRIPTION: the procedure integrates the log-advection equation, in the  !
  !conservative or non-conservative formulation, at a logarithmic grid, using!
  !a one-state second order scheme                                           !
  !==========================================================================!
  subroutine advance_log_advection(&
       CFLIn,         & 
       nP,            &
       nGCLeft,       &
       nGCRight,      &
       FInOut_I,      &
       IsConservative,&
       DeltaLnP)      

    !--------------------------------------------------------------------------!
    real,intent(in):: CFLIn           !Time step * acceleration rate/(Dlnp).   !
    integer,intent(in):: nP           !Number of meshes along lnp-coordinate.  !
    integer,intent(in):: nGCLeft      !The solution in the ghost cells is not  !
    integer,intent(in):: nGCRight     !advanced in time but used as the        ! 
    !boundary condition, if any              !
    real,dimension(1-nGCLeft:nP+nGCRight),intent(inout):: &
         FInOut_I                     !In: sol. to be advanced; Out: advanced  !
    !sol. index 0 is to set boundary         !
    !condition at the injection energy.      !
    logical,intent(in)::IsConservative!Solve (C) if .true. (NC) otherwise      !
    real,optional,intent(in)::DeltaLnP!Used only for NonConservative=.true.    !
    !--------------------------------------------------------------------------!
    integer:: iP,iStep,nStep
    real,dimension(1-max(nGCLeft,2):nP+max(nGCRight,2))::F_I
    real,dimension(0:nP+1):: FSemiintUp_I,FSemiintDown_I
    real::CFL,HalfADtIfNeeded
    !-------------------------NonConservative----------------------------------!
    ! This is a one-stage second-order scheme for the advection equation:      !
    !             f_t+A f_{ln p}=0,         .                                  !
    !                                                                          !
    ! Herewith CFL=A*Delta t/Delta (ln p):                                     !
    ! f^(n+1)_i-CFL*(f^n_(i-1/2)-f^n_(i+1/2)=f^n_i, where                      !
    ! For CFL>0:                                                               !
    ! f^n_(i-1/2)=f^n_(i-1)+cHalf*(cOne-CFL)*df_lim^n_(i-1)                    !
    ! For CFL<0:                                                               !
    ! f^n_(i-1/2)=f^n_(i  )-cHalf*(cOne+CFL)*df_lim^n_(i  )                    !
    !--------------------------------------------------------------------------!
    !-------------------------   Conservative-----------------------------------!
    ! This is a one-stage second-order scheme for the advection equation:      !
    !             f_t+A (f p}_p=0,         .                                   ! 
    ! We now need to know both                                                 !
    ! A*Delta t=CFLIn*Delta (LnP)                                              !
    ! and                                                                      !
    ! CFL=A*Delta t/(2 tanh(Delta (ln p)/2)                                    !
    ! f^(n+1)_i=CFL*(f^n_(i-1/2)-f^n_(i+1/2))-  (1/2)(f^n_(i-1/2)+f^n_(i+1/2))* !
    !                                   *A*Delta t+f^n_i,                      !
    !where                                                                     !
    ! For CFL>0:                                                               !
    ! f^n_(i-1/2)=f^n_(i-1)*(1-(1/2)*A*Delta t)+cHalf*(cOne-CFL)*df_lim^n_(i-1)!
    ! For CFL<0:                                                               !
    ! f^n_(i-1/2)=f^n_(i  )*(1-(1/2)*A*Delta t)-cHalf*(cOne+CFL)*df_lim^n_(i  )!
    !--------------------------------------------------------------------------!

    

    if(IsConservative)then
       HalfADtIfNeeded=0.50*CFLIn*DeltaLnP
       CFL=HalfADtIfneeded/tanh(0.50 * DeltaLnP)
       HalfADtIfNeeded=-HalfADtIfNeeded
    else
       HalfADtIfNeeded=0.0; CFL=CFLIn
    end if

    nStep=1+int(CFL); CFL=CFL/real(nStep)
    HalfADtIfNeeded=HalfADtIfNeeded/real(nStep)
    F_I(1-nGCLeft:nP+nGCRight)=FInOut_I(1-nGCLeft:nP+nGCRight)

    !Check for positivity
    if(any(F_I(1-nGCLeft:nP+nGCRight)<=0.0))then
       write(*,*)'Before advection F_I <=0'
       write(*,*)IsConservative
       write(*,*)F_I
       stop
    end if

    !One stage second order upwind scheme

    if (CFL>0.0) then
       do iStep=1,nStep

          !Boundary condition at the left boundary
          if(nGCLeft<2)F_I(            -1:0-nGCLeft) = F_I( 1-nGCLeft )
          !Boundary condition at the right boundary
          if(nGCRight<2)F_I(nP+1-nGCRight:nP+2     ) = F_I(nP+nGCRight)

          do iP=0,nP

             ! f_(i+1/2):
             FSemiintUp_I(iP) = F_I(iP)*(1.00 + HalfADtIfNeeded)+&
                  0.50*(1.00 - CFL)*df_lim(iP)
          end do
          ! f_(i-1/2): 
          FSemiintDown_I(1:nP) = FSemiintUp_I(0:nP-1)
          !\
          ! Update the solution from f^(n) to f^(n+1):
          !/
          
          F_I(1:nP) = F_I(1:nP)+CFL*(FSemiintDown_I(1:nP)-FSemiintUp_I(1:nP))+&
               HalfADtIfNeeded*(FSemiintDown_I(1:nP)+FSemiintUp_I(1:nP))
       end do
    else 
       do iStep=1,nStep
          !Boundary condition at the left boundary
          if(nGCLeft<2)F_I(            -1:0-nGCLeft) = F_I( 1-nGCLeft)
          !Boundary condition at the right boundary
          if(nGCRight<2)F_I(nP+1-nGCRight:nP+2     ) = F_I(nP+nGCRight)
          do iP=1,nP+1
             ! f_(i-1/2):
             FSemiintDown_I(iP) = F_I(iP)*(1.00 + HalfADtIfNeeded)&
                  -0.50 * (1.00 + CFL)*df_lim(iP)
          end do
          ! f_(i+1/2):
          FSemiintUp_I(1:nP) = FSemiintDown_I(2:nP+1)
          !\
          ! Update the solution from f^(n) to f^(n+1):
          !/
          F_I(1:nP) = F_I(1:nP)+CFL*(FSemiintDown_I(1:nP)-FSemiintUp_I(1:nP))+&
               HalfADtIfNeeded*(FSemiintDown_I(1:nP)+FSemiintUp_I(1:nP))
       end do
    end if

    FInOut_I(1:nP)=F_I(1:nP)
    if(any(FInOut_I(1:nP)<=0.0))then
       write(*,*)'After advection F_I <=0, for CFLFermi= ',CFL
       write(*,*)F_I
       stop
    end if
    !------------------------------------ DONE --------------------------------!
  Contains
    real function df_lim(i)
      integer, intent(in):: i
      !------------------------------------------------------------------------!
      real:: dF1,dF2
      !------------------------------------------------------------------------!
      dF1 = F_I(i+1)-F_I(i)
      dF2 = F_I(i)-F_I(i-1)

      !df_lim=0 if dF1*dF2<0, sign(dF1) otherwise:

      df_lim = sign(0.50,dF1)+sign(0.50,dF2)
      dF1 = abs(dF1)
      dF2 = abs(dF2)
      df_lim = df_lim*min(max(dF1,dF2),2.0*dF1,2.0*dF2)
      !---------------------------------- DONE --------------------------------!
    end function df_lim
  end subroutine advance_log_advection
!============================================================================!
end module ModLogAdvection
