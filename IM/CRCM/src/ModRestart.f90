Module ModCrcmRestart
  implicit none
  
  SAVE

  logical :: IsRestart = .false.
contains
  !============================================================================

  subroutine crcm_read_restart
    use ModCrcmPlanet,ONLY: nspec
    use ModCrcmGrid,  ONLY: np,nt,nm,nk
    use ModCrcm,      ONLY: f2!, phot, Pressure_IC, FAC_C
    use ModIoUnit,    ONLY: UnitTmp_
    !--------------------------------------------------------------------------
    open(unit=UnitTmp_,file='IM/restartIN/data.restart',status='old',form='unformatted')
    read(UnitTmp_) f2             
  end subroutine crcm_read_restart
  

  !============================================================================
  subroutine crcm_write_restart
    use ModCrcmPlanet,ONLY: nspec
    use ModCrcmGrid,  ONLY: np,nt,nm,nk
    use ModCrcm,      ONLY: f2!, phot, Pressure_IC, FAC_C
    use ModIoUnit,    ONLY: UnitTmp_
    !--------------------------------------------------------------------------
    open(unit=UnitTmp_,file='IM/restartOUT/data.restart',form='unformatted')
    write(UnitTmp_) f2

  end subroutine crcm_write_restart

end Module ModCrcmRestart
