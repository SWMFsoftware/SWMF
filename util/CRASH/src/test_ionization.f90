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

  integer,parameter :: nN=5 , nT=10, nU=10
  real    :: dTe, dLogN, dU
  integer :: iT, nT1=1000000, iN, iU, iLoop
  
  real    :: z_I(0:nN),z2_I(0:nN),Uav_I(0:nN),Cv_I(0:nN), Te_I(0:nN) 
  integer :: iIter_I(0:nN)
  !character(LEN=*),parameter,dimension(0:nN) :: Separator_I='|'
  !character(LEN=*),parameter,dimension(0:nN) :: Separator1_I='/'
  logical :: IsDegenerated
  !-------------------------------------------

  dTe = 50.0; dLogN=log(10.0); dU=1000.0

  
  call set_element( 54 )

!  tm_1 = diff_sec()!
!  nT1 =  (nN +1)*nT/1000000 
!  write(*,*)"Start,", nT1 , " million iterations"



  open(24,file='../doc/Table1.tex')
        write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Te[eV]\textbackslash \textbackslash Na[$1/cm^3$] & $10^{18}$'//&
                      ' & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(24,'(a)')'\hline' 
        write(24,'(a)')'\hline'
  

  do iT  = 1,nT
     if (((iT-1)/25)*25==(iT-1).and.iT>25) then
        
        write(24,'(a)')'\end{tabular}', char(10)
        !--------------
        write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(24,'(a)')'\hline'
        write(24,'(a)')'Te[eV]\textbackslash \textbackslash Na[$1/cm^3$] & $10^{18}$ & $10^{19}$'//&
                       ' & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
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
     write(24,'(f5.0,6(a,f7.1),a)') vTe,&
               (' & ', Z_I(iN), iN=0,nN ),'\tabularnewline'
     write(24,'(a)')'\hline'


  end do
  
  write(24,'(a)')'\end{tabular}'

  close(24)
!_______________________________________________
open(24,file='../doc/Table2.tex')
  write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(24,'(a)')'\hline'
  write(24,'(a)')'Te[eV] & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ '//&
                 '& $U_{av} | C_v$ & $U_{av} | C_v$\tabularnewline'
  write(24,'(a)')'\hline' 
  write(24,'(a)')'\hline'
  

  do iT  = 1,nT
     if (((iT-1)/25)*25==(iT-1).and.iT>25) then
        
        write(24,'(a)')'\end{tabular}', char(10)
        !--------------
		write(24,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
		write(24,'(a)')'\hline'
		write(24,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
		write(24,'(a)')'\hline'
		write(24,'(a)')'Te[eV] & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ & $U_{av} | C_v$ '//&
					'& $U_{av} | C_v$ & $U_{av} | C_v$\tabularnewline'
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
     write(24,'(f5.0,6(a,f8.1,a,f7.1),a)') vTe,&
               (' & ', Uav_I(iN), ' | ', Cv_I(iN), iN=0,nN ),'\tabularnewline'
     write(24,'(a)')'\hline'


  end do
  
  write(24,'(a)')'\end{tabular}'

  close(24)

!_____________________________________

  open(25,file='../doc/Table3.tex')
  write(25,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
  write(25,'(a)')'\hline'
  write(25,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
  write(25,'(a)')'\hline'
  write(25,'(a)')'U[eV] & Te (Iterations) &  Te (Iterations) &  Te (Iterations) & '//&
        ' Te (Iterations) &  Te (Iterations) &  Te (Iterations) \tabularnewline'
  write(25,'(a)')'\hline' 
  write(25,'(a)')'\hline'
  
  do iU  = 1,nU
     if (((iU-1)/50)*50==(iU-1).and.iU>50) then
        write(25,'(a)')'\end{tabular}', char(10)
        !------------
        write(25,'(a)')'\begin{tabular}{|c||c|c|c|c|c|c|}'
        write(25,'(a)')'\hline'
        write(25,'(a)')'Na[$1/cm^3$] & $10^{18}$ & $10^{19}$ & $10^{20}$ & $10^{21}$ & $10^{22}$ & $10^{23}$\tabularnewline'
        write(25,'(a)')'\hline'
        write(25,'(a)')'U[eV] & Te (Iterations) &  Te (Iterations) &  Te (Iterations) & '//&
             ' Te (Iterations) &  Te (Iterations) &  Te (Iterations) \tabularnewline'
        write(25,'(a)')'\hline' 
        write(25,'(a)')'\hline'
     end if
     vU = dU * iU 
     do iN = 0,nN
        NaTrial = Nao*exp(iN*dLogN)
        call set_temperature(vU,NaTrial*1000000.0,IsDegenerated)
        Te_I(iN) = Te 
        iIter_I(iN) = iIterTe
        if(IsDegenerated)Te_I(iN)=-1.0
     end do

     write(25,'(f6.0,6(a,f7.1,a,i7,a),a)') vU,&
               (' & ', Te_I(iN), ' (', iIter_I(iN),')', iN=0,nN ),'\tabularnewline'
     write(25,'(a)')'\hline'

  end do


  write(25,'(a)')'\end{tabular}'
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


