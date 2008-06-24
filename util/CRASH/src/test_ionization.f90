!^CFG COPYRIGHT UM

!***********************************************************************
!    calculation of ionization equilibrium, single material
!    for given - concentration of  atomic particles Co+Ciz [cm^{-3}]
!              - electron temperature Te[eV]
!***********************************************************************
program saha
  use ModStatSum
  implicit NONE
  real :: &
       Nao = 1.00e18,  &  ! cm-3
       vTe=5., NaTrial, vU 

  integer,parameter :: nN=5 , nT=120, nU=100
  real    :: dTe, dLogN, dU
  integer :: iT, nT1=1000000, iN, iU, iLoop
  
  real    :: z_I(0:nN),z2_I(0:nN),Uav_I(0:nN),Cv_I(0:nN), Te_I(0:nN), iIter_I(0:nN)
  !character(LEN=*),parameter,dimension(0:nN) :: Separator_I='|'
  !character(LEN=*),parameter,dimension(0:nN) :: Separator1_I='/'
  logical :: IsDegenerated
  !-------------------------------------------
  dTe = 5.0; dLogN=log(10.0)

  
  call set_element( 54 )

!  tm_1 = diff_sec()!
!  nT1 =  (nN +1)*nT/1000000 
!  write(*,*)"Start,", nT1 , " million iterations"



  open(24,file='../doc/Table1.tex',status='unknown')
  write(24,'(a)')'\begin{table}'
  write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(24,'(a)')'\hline'



  write(*,'(a)')' \ Na |   10^18 |   10^19 |   10^20 |'//& 
'   10^21 |   10^22 |   10^23 |'

  write(24,'(a)')'Na/cm3 & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Te[eV] & z | $<z^2>/z$ & z | $<z^2>/z$ & z | $<z^2>/z$ & z | $<z^2>/z$ &'//&
                 'z | $<z^2>/z$ & z | $<z^2>/z$\tabularnewline'
  write(24,'(a)')'& Uav | Cv &  Uav | Cv &  Uav | Cv &  Uav | Cv &  Uav | Cv &'//&
                 'Uav | Cv\tabularnewline'
  write(24,'(a)')'\hline' 
  write(24,'(a)')'\hline'
  


  do iT  = 1,nT
     if (count((/(iT == iLoop*25+1, iLoop = 1,10)/))>0) then
        write(24,'(a)')'\end{tabular}'
        write(24,'(a)')'\end{table}', char(10)
        !--------------
        write(24,'(a)')'\begin{table}'
        write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Na/cm3 & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Te[eV] & $z | <z^2>/z$ & $z | <z^2>/z$ & $z | <z^2>/z$ & $z | <z^2>/z$ &'//&
                        '$z | <z^2>/z$ & $z | <z^2>/z$\tabularnewline'
        write(24,'(a)')'& Uav | Cv &  Uav | Cv &  Uav | Cv &  Uav | Cv &  Uav | Cv &'//&
                 'Uav | Cv\tabularnewline'
        write(24,'(a)')'\hline' 
        write(24,'(a)')'\hline'
     end if
     vTe = dTe * iT 

     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call set_ionization_equilibrium(vTe,NaTrial*1000000.0,IsDegenerated)
        Z_I(iN) = z_averaged() 
        Z2_I(iN)= z2_averaged()/Z_I(iN)
        Uav_I(iN)=internal_energy()
        Cv_I(iN)=heat_capacity()
        if(IsDegenerated)then
           Z_I(iN) = -1.0
           Z2_I(iN)= -1.0
           Uav_I(iN)= -1.0
           Cv_I(iN)= -1.0
        end if
     end do
     write(* ,'(a,f5.0,a,6(f4.1,a,f4.1,a))')'|',vTe,'|',&
               (Z_I(iN), '/', Z2_I(iN), '|', iN=0,nN )
     write(24,'(f5.0,6(a,f7.1,a,f7.1),a)') vTe,&
               (' & ', Z_I(iN), ' | ', Z2_I(iN), iN=0,nN ),'\tabularnewline'
     write(24,'(6(a,f7.1,a,f7.1),a)') &
               (' & ', Uav_I(iN), ' | ', Cv_I(iN), iN=0,nN ),'\tabularnewline'
     write(24,'(a)')'\hline'


  end do
  



  write(24,'(a)')'\end{tabular}'
  write(24,'(a)')'\end{table}'

  close(24)

!_____________________________________

 open(25,file='../doc/Table2.tex',status='unknown')
  write(25,'(a)')'\begin{table}'
  write(25,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(25,'(a)')'\hline'
  write(25,'(a)')'Na/cm3 & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(25,'(a)')'\hline'
  write(25,'(a)')'U[eV] & Te | Iterations &  Te | Iterations &  Te | Iterations & '//&
                 ' Te | Iterations &  Te | Iterations &  Te | Iterations \tabularnewline'
  write(25,'(a)')'\hline' 
  write(25,'(a)')'\hline'
  
  do iU  = 1,nU
     if (count((/(iU == iLoop*25+1, iLoop = 1,10)/))>0) then
        write(25,'(a)')'\end{tabular}'
        write(25,'(a)')'\end{table}', char(10)
        !------------
        write(25,'(a)')'\begin{table}'
        write(25,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(25,'(a)')'\hline'
        write(25,'(a)')'Na/cm3 & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(25,'(a)')'\hline'
        write(25,'(a)')'U[eV] & Te | Iterations &  Te | Iterations &  Te | Iterations & '//&
                 ' Te | Iterations &  Te | Iterations &  Te | Iterations \tabularnewline'
        write(25,'(a)')'\hline' 
        write(25,'(a)')'\hline'
     end if
     vU = dU * iU 
     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call set_temperature(vU,NaTrial*1000000.0)
        Te_I(iN) = Te 
        iIter_I(iN) = iIterTe
        !if(IsDegenerated)then
        !   Z_I(iN) = -1.0
        !   Z2_I(iN)= -1.0
        !end if
     end do

     write(25,'(f5.0,6(a,f7.1,a,f7.1),a)') vU,&
               (' & ', Te_I(iN), ' | ', iIter_I(iN), iN=0,nN ),'\tabularnewline'
     write(25,'(a)')'\hline'

  end do


  write(25,'(a)')'\end{tabular}'
  write(25,'(a)')'\end{table}'
  close(25)



contains
  !================================== current time, in SECONDs 
  function  diff_sec  ( ) result (sec) 
    integer,dimension (8) :: val
    integer               :: sec
    call date_and_time(VALUES=val)    
    write (*,'("time ",i2.2,":",i2.2,":",i2.2))') val(5:7)
    sec = val(7) +(val(6) +val(5)*60)*60                
  end function diff_sec 

end program saha

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test


