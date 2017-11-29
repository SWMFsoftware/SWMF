!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program test
  use ModHydrogenicModel
  implicit none
  character(LEN=120)::NameLong
  real:: Energy0, Energy
  call read_sigma_coef
  !=====================HYDROGEN============
  write(*,*)' HYDROGEN '
  write(*,*)
  nZ=1
  nElectron_I=0
  nElectron_I(1)=1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy0 = configuration_energy()
  write(*,*)Energy0
  write(*,*)
  nElectron_I(1)=0;nElectron_I(2)=1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  nElectron_I(2)=0;nElectron_I(3)=1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  nElectron_I(3)=0;nElectron_I(4)=1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
 !=====================HELIUMN============
  write(*,*)' HELIUM '
  write(*,*)
  nZ=2
  nElectron_I=0
  nElectron_I(1)=2
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy0 = configuration_energy()
  write(*,*)Energy0
  write(*,*)
  nElectron_I(1)=1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  nElectron_I(2)=1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  nElectron_I(2)=0;nElectron_I(3)=1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  !=========================================
  write(*,*)' OXYGEN '
  write(*,*)
  nZ=8
  nElectron_I=0
  nElectron_I(1:4) =2 
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy0 = configuration_energy()
  write(*,*)Energy0
  write(*,*)
  nElectron_I(4) = 1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  Energy0=Energy
  nElectron_I(4) = 0
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  Energy0=Energy
  nElectron_I(3) = 1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  Energy0=Energy
  nElectron_I(3) = 0
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  Energy0=Energy
  nElectron_I(2) = 1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  Energy0=Energy
  nElectron_I(2) = 0
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  Energy0=Energy
  nElectron_I(1) = 1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  write(*,*)
  Energy0=Energy
  write(*,*)-Energy0
  !=========================================
  
  nZ=54
  nElectron_I(1:19) = (/2,&
                    2,2,4,&
                    2,2,4,4,6,&
                    2,2,4,4,6,0,0,&
                    2,2,4/)
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy0 = configuration_energy()
  write(*,*)Energy0
  nElectron_I(19) = 3; nElectron_I(26) = 1
  NameLong = TypeConf()
  write(*,'(a)')NameLong
  Energy = configuration_energy()
  write(*,*)Energy-Energy0
  
end program test
!============================================================================
