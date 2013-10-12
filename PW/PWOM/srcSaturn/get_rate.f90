!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!******************************************************************************
! get_rate takes Alt as input and returns the rate coef for the 
! H+ + H2 --> H2+ + H reaction. Data from Moses and Bass, 2000
! Created by Alex Glocer, 07/06
!******************************************************************************

Subroutine get_rate( Alt,Rate)

  real, intent(in)  :: Alt
  real, intent(out) :: Rate

  
  real,dimension(15)  :: MeasuredAlt, MeasuredRate
  
  
  MeasuredRate(1)=2.02711E-18
  MeasuredAlt(1) =8.00E+2
  
  MeasuredRate(2)=8.06354E-18
  MeasuredAlt(2) =8.73E+2
  
  MeasuredRate(3)=6.24330E-17
  MeasuredAlt(3) =9.91E+2
  
  MeasuredRate(4)=1.65461E-16
  MeasuredAlt(4) =1.091E+3
  
  MeasuredRate(5)=3.01799E-16
  MeasuredAlt(5) =1.191E+3

  MeasuredRate(6)=5.32882E-16
  MeasuredAlt(6) =1.310E+3
  
  MeasuredRate(7)=1.00407E-15
  MeasuredAlt(7) =1.414E+3
  
  MeasuredRate(8)=2.01890E-15
  MeasuredAlt(8) =1.524E+3
  
  MeasuredRate(9)=3.80406E-15
  MeasuredAlt(9) =1.615E+3
  
  MeasuredRate(10)=9.60204E-15
  MeasuredAlt(10) =1.706E+3
  
  MeasuredRate(11)=1.86898E-14
  MeasuredAlt(11) =1.8828E+3
  
  MeasuredRate(12)=4.56680E-14
  MeasuredAlt(12) =2.020E+3
  
  MeasuredRate(13)=6.63544E-14
  MeasuredAlt(13) =2.224E+3
  
  MeasuredRate(14)=7.55626E-14
  MeasuredAlt(14) =2.402E+3
  
  MeasuredRate(15)=7.80578E-14
  MeasuredAlt(15) =3.076E+3
  
  if (Alt .ge. MeasuredAlt(15)) then
     Rate = MeasuredRate(15)
     return
  endif
  
  if (Alt .lt. MeasuredAlt(1)) then
     Rate = MeasuredRate(1)
     return
  endif
  
  
  do i=1,14
     if (Alt .ge. MeasuredAlt(i) .and. Alt .lt. MeasuredAlt(i+1)) then
        Rate = &
          (MeasuredRate(i+1)-MeasuredRate(i))  &
         /(MeasuredAlt(i+1)-MeasuredAlt(i)) &
         *(Alt-MeasuredAlt(i)) + MeasuredRate(i)
        return
     endif
  enddo

end Subroutine get_rate
