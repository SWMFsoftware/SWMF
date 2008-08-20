!============================================================================!
subroutine write_logfile_SP(TypeActionLogFile)
  use SP_ModMain
  use ModIOUnit
  implicit none
  character(LEN=*),intent(in):: TypeActionLogFile
  !--------------------------------------------------------------------------!
  character(LEN=*),parameter:: NameSub=' write_logfile_SP'
  character(LEN=30):: NameLogFile
  integer,parameter:: nEChannel=6
  integer:: iError,iEChannel,iLnP
  integer,save:: iFile
  real,allocatable,save,dimension(:):: P_I,E_I
  real,dimension(nEChannel),save:: EnergyLowLimit_I
  real,dimension(nEChannel):: TotalDEFlux_I
  !--------------------------------------------------------------------------!
 
  select case(TypeActionLogFile)
  case('OPEN','open')
     !-----------------------------------------------------------------------!
     write(iStdOut,*)prefix//NameSub//' '//TypeActionLogFile
     allocate(P_I(1:nP),stat=iError)
     call check_allocate(iError,NameSub//'P_I')
     allocate(E_I(1:nP),stat=iError)
     call check_allocate(iError,NameSub//'E_I')
     !-----------------------------------------------------------------------!
     iFile=io_unit_new()
     write(NameLogFile,'(a)')trim(SP_DirOut)//&
          'SP_logfile.dat'
     write(iStdOut,*)prefix//'Open '//trim(NameLogFile)
     open(iFile,file=trim(NameLogFile),status='unknown',&
          form='formatted')
     write(iFile,'(a)')'  SP_Time,   Flux(E>5MeV),'//&
          '   Flux(E>10MeV),   Flux(E>30MeV),'//&
          '   Flux(E>50MeV),   Flux(E>60MeV),'//&
          '   Flux(E>100MeV)'   
     !-----------------------------------------------------------------------!
     EnergyLowLimit_I(1) =   5.0           !5   In units of [MeV].
     EnergyLowLimit_I(2) =  10.0           !10  In units of [MeV].
     EnergyLowLimit_I(3) =  30.0           !30  In units of [MeV].
     EnergyLowLimit_I(4) =  50.0           !50  In units of [MeV].
     EnergyLowLimit_I(5) =  60.0           !60  In units of [MeV].
     EnergyLowLimit_I(6) = 100.0           !100 In units of [MeV].
     do iEChannel=1,nEChannel
        EnergyLowLimit_I(iEChannel) = &
             EnergyLowLimit_I(iEChannel)*&
             energy_in('MeV')            !In units of [kg*m^2/s^2].
     end do
     !-----------------------------------------------------------------------!
     do iLnP=1,nP
        P_I(iLnP) = PInjection*exp(real(iLnP)* &
             DeltaLnP)                   !This gives P in [kg*m^2/s^2].
        E_I(iLnP) = momentum_to_kinetic_energy(&
             P_I(iLnP),NameParticle)     !This gives E in [J].
     end do
     !-----------------------------------------------------------------------!
  case('WRITE','write')
     !-----------------------------------------------------------------------!
     do iEChannel=1,nEChannel
        call integrate_diff_energy_flux(P_I,E_I,&
             EnergyLowLimit_I(iEChannel),&
             TotalDEFlux_I(iEChannel))
     end do
     !-----------------------------------------------------------------------!
     write(iFile,'(e11.4,30(1X,e13.5))')SP_Time,&
          TotalDEFlux_I(1:nEChannel)
     !-----------------------------------------------------------------------!
  case('CLOSE','close')
     write(iStdOut,*)prefix//NameSub//' '//TypeActionLogFile
     close(iFile)
     deallocate(P_I,E_I)
  end select
end subroutine write_logfile_SP
!============================================================================!
subroutine integrate_diff_energy_flux(P_I,E_I,ELow,TotalDEFlux)
  use SP_ModMain
  implicit none
  !--------------------------------------------------------------------------!
  real,intent(in):: ELow
  real,dimension(nP),intent(in):: P_I,E_I
  real,intent(out):: TotalDEFlux
  !--------------------------------------------------------------------------!
  logical:: IsFound
  integer:: iLnP
  !--------------------------------------------------------------------------!
  IsFound=.false.
  do iLnP = 1,nP-1
     if ((E_I(iLnP)>ELow).and.(.not.IsFound)) then
        IsFound = .true.
        TotalDEFlux =(&
             (0.50*(E_I(iLnP+1)+ELow)-E_I(iLnP))*&
                     F_II(iLnP+1,nX)*P_I(iLnP+1)**2+ &
              0.50*(E_I(iLnP+1)-ELow)*&
                     F_II(iLnP  ,nX)*P_I(iLnP  )**2)*&
                    (E_I(iLnP+1)-ELow)/(E_I(iLnP+1)-E_I(iLnP))
     end if
     if (IsFound) then
        TotalDEFlux = TotalDEFlux+&
             0.50*(E_I (iLnP+1) -  E_I(iLnP  ))*  &
                   (F_II(iLnP  ,nX)*P_I(iLnP  )**2+&
                    F_II(iLnP+1,nX)*P_I(iLnP+1)**2)
     end if
  end do
end subroutine integrate_diff_energy_flux
!============================================================================!
