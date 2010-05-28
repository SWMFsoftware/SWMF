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
    close(UnitTmp_)
  end subroutine crcm_read_restart
  

  !============================================================================
  subroutine crcm_write_restart
    use ModCrcmPlanet,ONLY: nspec
    use ModCrcmGrid,  ONLY: np,nt,nm,nk
    use ModCrcm,      ONLY: f2,time!, phot, Pressure_IC, FAC_C
    use ModIoUnit,    ONLY: UnitTmp_
    !--------------------------------------------------------------------------
    open(unit=UnitTmp_,file='IM/restartOUT/data.restart',form='unformatted')
    write(UnitTmp_) f2
    close(UnitTmp_)

    open(unit=UnitTmp_,file='IM/restartOUT/restart.H')
    write(UnitTmp_,'(a)') '#TIMESIMULATION'
    write(UnitTmp_,'(es15.8,a25)') time,'tSimulation'
    close(UnitTmp_)
  end subroutine crcm_write_restart

end Module ModCrcmRestart
