
subroutine init_msis

  implicit none

  call report("init_msis",1)

  write(*,*) "You can not use MSIS with any planet except Earth!!!"
  write(*,*) "If you ARE running Earth, then make the code again, using"
  write(*,*) "configure Earth ; make"
  call stop_gitm("I can not continue...")

end subroutine init_msis

subroutine msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS, Temp, LogRho)

  write(*,*) "You can not use MSIS with any planet except Earth!!!"
  write(*,*) "If you ARE running Earth, then make the code again, using"
  write(*,*) "configure Earth ; make"
  call stop_gitm("I can not continue...")

end subroutine msis_bcs
