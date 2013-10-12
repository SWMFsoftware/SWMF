!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! File name: dgcpm_fields_011.f9
!
! Contains: Volland-Stern KP Driven Fields, Simple Potential Model, 
!           Sonjka Potential Model
!
! Last Modified: April 2011, Aron Dodger
!
! **********************************************************************
!                             VOLLAND_STERN
!	Calculates magnetospheric convection strength and puts it
!	into the drift terms (added to already-calculated 
!       corotation and gradient-curvature)
! **********************************************************************
SUBROUTINE VOLLAND_STERN

  use ModSizeDGCPM
  use ModMainDGCPM
  use ModIoDGCPM
  use ModConstants

  implicit none

  integer j,i

  !  Calculate base convection electric field
  !  CC Shielded Volland-Stern

  do J=1,nphicells   ! Fill in DGCPM potentials
!     write(*,*) j, vphicells(j)
     do I=1,nthetacells
!        if (j == nphicells) write(*,*) i, vlzcells(i)
        potdgcpm(i,j)=a*re*vlzcells(i)**(lamgam)*sin(dtor*vphicells(J))
     enddo
  enddo

  if (debug .gt. 0) write(*,*) "cpcp : ", maxval(potdgcpm)-minval(potdgcpm)
  return

end subroutine volland_stern

! **********************************************************************
!                             EPOTSIMPLE
! **********************************************************************
subroutine epotsimple()

    use ModSizeDGCPM
    use ModMainDGCPM
    use ModTimeDGCPM

! input: dtheta, and dphi in degrees
! phi is zero at 24 MLT positive towards dawn
! theta is zero at the pole

!      Assuming a uniform dawn-dusk electric field
!      and an assumed kp relationship
!
!     Stagnation point (Re)     Electric Field (V/Re)   Kp
!           10                         919               1
!            9                        1134
!            8                        1436
!            7                        1875
!            6                        2552
!            5                        3675               7

      real kp
      real pot
      real pi, rad
      real chi,dr,sint
      real dx
      integer i, j

    pi = 3.14159          ! rad
    rad = pi / 180.0      ! rad/degree

    i3 = tSimulation/(2.0*dt)
    kp = rkph(i3)

    chi = 7350.0 / (9.0 - kp)

    do i=1, nthetacells
        do j=1, nphicells
        
            sint = sin(vthetacells(i) * rad)
            dr = 1.0 / (sint * sint)
            dx = dr * sin((180.0+vphicells(j))*rad)

            mgridpot(i,j) = -chi * dx
        enddo
    enddo

    return
end subroutine EPOTSIMPLE

! **********************************************************************
!                             EPOTSOJKA
! **********************************************************************

subroutine epotsojka()

    use ModSizeDGCPM
    use ModMainDGCPM
    use ModTimeDGCPM

! input: dphi,dtheta (spherical coordinates) in degrees
! theta is zero at the pole and 
! Input: Kp index
      real kp

! output: pot potential in volts
      real pot

      real kp
      real pi, rad
      real theta, phi
      real theta_eq, theta_pc, theta_max, chi_0
      real phix, theta_3, chi_pc, chi
      real theta_1, theta_2

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree

      i3 = tSimulation/(2.0*dt)
      kp = rkph(i3)

    do i=1, nthetacells
    do j=1, nphicells

        theta=vthetacells(i)*rad
        phi=vphicells(j)*rad

        theta_eq=(25.0+(2.0*kp))*rad
        theta_pc=(15.0+(0.3*kp))*rad
        theta_max=(theta_pc)+(0.3*(theta_eq-theta_pc)*(abs(sin(phi))))
        chi_0=(10.0+(6.5*kp))*1000.0

      if ((theta.lt.theta_pc).and.(theta.ge.0.0)) then

       phix=asin(sqrt(1.0-((sin(theta)*sin(theta)*
     *    cos(phi)*cos(phi))/(sin(theta_pc)*sin(theta_pc)))))
       theta_3=((sin(theta)*sin(phi))/(sqrt((sin(theta_pc)*
     *    sin(theta_pc))-(sin(theta)*sin(theta)*cos(phi)*cos(phi)))))

       if  ((phix.ge.0.0).and.(phix.le.(pi/3.0))) then
        chi_pc=chi_0*sin((3.0*phix)/(2.0))
        chi=chi_pc*theta_3

       else if ((phix.gt.(pi/3.0)).and.
     *    (phix.le.((2.0*pi)/3.0))) then
        chi_pc=chi_0
        chi=chi_pc*theta_3

       else if ((phix.gt.((2.0*pi)/3.0)).and.
     *    (phix.le.((4.0*pi)/3.0))) then
        chi_pc=-chi_0*cos((3.0*phix)/(2.0))
        chi=chi_pc*theta_3

       else if ((phix.gt.((4.0*pi)/3.0)).and.
     *    (phix.le.((5.0*pi)/3.0))) then
        chi_pc=-chi_0
        chi=chi_pc*theta_3

       else if ((phix.gt.((5.0*pi)/3.0)).and.
     *    (phix.le.((6.0*pi)/3.0))) then
        chi_pc=-chi_0*sin((3.0*phix)/(2.0))
        chi=chi_pc*theta_3

       end if

      else if ((theta.lt.theta_max).and.
     *   (theta.ge.theta_pc)) then

       theta_1=(1.0-(((theta-theta_pc)*(theta-theta_pc))/
     *    ((theta_max-theta_pc)*(theta_eq-theta_pc))))

       if ((phi.ge.0.0).and.(phi.le.(pi/3.0))) then
        chi_pc=chi_0*sin((3.0*phi)/(2.0))
        chi=chi_pc*theta_1

       else if ((phi.gt.(pi/3.0)).and.
     *    (phi.le.((2.0*pi)/3.0))) then
        chi_pc=chi_0
        chi=chi_pc*theta_1

       else if ((phi.gt.((2.0*pi)/3.0)).and.
     *    (phi.le.((4.0*pi)/3.0))) then
        chi_pc=-chi_0*cos((3.0*phi)/(2.0))
        chi=chi_pc*theta_1

       else if ((phi.gt.((4.0*pi)/3.0)).and.
     *    (phi.le.((5.0*pi)/3.0))) then
        chi_pc=-chi_0
        chi=chi_pc*theta_1

       else if ((phi.gt.((5.0*pi)/3.0)).and.
     *    (phi.le.((6.0*pi)/3.0))) then
        chi_pc=-chi_0*sin((3.0*phi)/(2.0))
        chi=chi_pc*theta_1

       end if

      else if ((theta.lt.theta_eq).and.
     *   (theta.ge.theta_max)) then

       theta_2=((((theta-theta_eq)*(theta-theta_eq))/
     *    ((theta_eq-theta_max)*(theta_eq-theta_pc))))

       if  ((phi.ge.0.0).and.(phi.le.(pi/3.0))) then
        chi_pc=chi_0*sin((3.0*phi)/(2.0))
        chi=chi_pc*theta_2

       else if ((phi.gt.(pi/3.0)).and.
     *    (phi.le.((2.0*pi)/3.0))) then
        chi_pc=chi_0
        chi=chi_pc*theta_2

       else if ((phi.gt.((2.0*pi)/3.0)).and.
     *    (phi.le.((4.0*pi)/3.0))) then
        chi_pc=-chi_0*cos((3.0*phi)/(2.0))
        chi=chi_pc*theta_2

       else if ((phi.gt.((4.0*pi)/3.0)).and.
     *    (phi.le.((5.0*pi)/3.0))) then
        chi_pc=-chi_0
        chi=chi_pc*theta_2

       else if ((phi.gt.((5.0*pi)/3.0)).and.
     *    (phi.le.((6.0*pi)/3.0))) then
        chi_pc=-chi_0*sin((3.0*phi)/(2.0))
        chi=chi_pc*theta_2

       end if

      else if ((theta.le.(pi/2.0)).and.
     *   (theta.ge.theta_eq)) then
       chi=0.0
      end if

      mgridpot(i,j)=chi


    enddo
    enddo
      return
end subroutine epotsojka
