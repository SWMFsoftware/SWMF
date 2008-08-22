!=============================================================================
subroutine write_plotfile_SP(DoPlot,TypeOutput)
  use ModIOUnit
  use SP_ModMain
  use ModTurbulence
  implicit none
  !--------------------------------------------------------------------------!
  integer:: iFile
  logical,intent(in):: DoPlot
  character(LEN=50),intent(in):: TypeOutput
  !--------------------------------------------------------------------------!
  integer,parameter:: nVarMax=10
  integer,parameter:: nTecplotVersion=10
  integer,parameter:: nPSkip=2,nXSkip=2
  integer:: iX,iLnP,iVar,nVar,iError,iLength
  character(LEN=3)  :: TypePlot
  character(LEN=36) :: NameFile
  character(len=50) :: SP_NameAllVars
  character(LEN=3)  :: SP_NamePlotVar
  character(len=10) :: SP_NamePlotVar_I(nVarMax)
  character(LEN=4)  :: SP_NameExtension
  character(LEN=500):: SP_NameUnitTec='',SP_NameUnitIdl=''
  character(LEN=20) :: SP_NameUnitTec_G
  character(LEN=20) :: SP_NameUnitTec_X,SP_NameUnitTec_E
  character(LEN=70) :: SP_NameUnitTec_DEF,SP_NameUnitTec_CDF
  character(LEN=20) :: SP_NameUnitIdl_G
  character(LEN=20) :: SP_NameUnitIdl_X,SP_NameUnitIdl_E
  character(LEN=70) :: SP_NameUnitIdl_DEF,SP_NameUnitIdl_CDF
  real:: Location,Momentum,Energy
  real:: SP_Unit_X,SP_Unit_E
  real,allocatable:: PlotVar_III(:,:,:)
  integer::iLnK
  !--------------------------------------------------------------------------!
  if (.not.DoPlot) return
  iFile=io_unit_new()
  SP_iPlot = int(SP_Time/SP_TimePlot)
  !\
  ! Set units for spatial and energy coordinates.
  !/
  SP_Unit_X = 1.0 !/Rsun
  SP_Unit_E = 1.0/energy_in('KeV')
  !\
  ! Set strings of units for Tecplot & IDL variables.
  !/
  if (nTecplotVersion==10) then
     SP_NameUnitTec_X   = '[R<sub>s</sub>]'
     SP_NameUnitTec_E   = '[keV]'
     SP_NameUnitTec_G   = '[   ]'
     SP_NameUnitTec_CDF = &
          '[s<sup>3</sup> kg<sup>-3</sup> m<sup>-6</sup>]'
     SP_NameUnitTec_DEF = &
          '[s<sup>-1</sup> m<sup>-2</sup> str<sup>-1</sup> keV<sup>-1</sup>]'
  else
     SP_NameUnitTec_X   = '[R_s]'
     SP_NameUnitTec_E   = '[keV]'
     SP_NameUnitTec_CDF = '[s^3 kg^-^3 m^-^6]'
     SP_NameUnitTec_DEF = '[s^-^1 m^-^2 str^-^1 keV^-^1]'
  end if
  SP_NameUnitIdl_X   = SP_NameUnitTec_X
  SP_NameUnitIdl_E   = SP_NameUnitTec_E
  SP_NameUnitIdl_G   = SP_NameUnitTec_G
  SP_NameUnitIdl_DEF = SP_NameUnitTec_DEF
  SP_NameUnitIdl_CDF = SP_NameUnitTec_CDF
  !\
  ! Resolve the test string for data output.
  !/
  call sp_split_str(TypeOutput,nVarMax,SP_NamePlotVar_I,nVar)
  TypePlot = SP_NamePlotVar_I(nVar); nVar = nVar-1
  !\
  ! Allocate array for plot variables.
  !/
  allocate(PlotVar_III(nP,nX,nVar),stat=iError)
  call check_allocate(iError,'SP: PlotVar_III')
  !\
  ! Construct a string of variable names for IDL output.
  !/
  SP_NameAllVars = 'X E'
  iLength=3
  do iVar=1,nVar
     SP_NameAllVars(iLength+1:iLength+1+&
          len_trim(SP_NamePlotVar_I(iVar)))=' '//&
          trim(SP_NamePlotVar_I(iVar))
     iLength = iLength+1+len_trim(SP_NamePlotVar_I(iVar))
  end do
  !\
  ! Select the type of plot file: IDL or Tecplot
  ! and also the type of file header & units.
  !/
  select case(TypePlot)
  case('TEC','tec')
     SP_NameExtension = '.dat'
     call sp_get_tec_units
  case('IDL','idl')
     SP_NameExtension = '.out'
     call sp_get_idl_units
  end select
  !\
  ! Construct the name of output file.
  !/
  write(NameFile,'(a,i4.4,a)')trim(SP_DirOut)//'SP_data_n',SP_iPlot,&
       SP_NameExtension
  write(iStdout,*)prefix,' '
  write(iStdout,*)prefix,'Opening File '//trim(NameFile)//' to Save SP Data'//&
       ' at SP_Time =',SP_Time
  write(iStdout,*)prefix,' '
  !\
  ! Open file unit and write plot file.
  !/
  open(iFile,file=trim(NameFile),status='replace',form='formatted')
  select case(TypePlot)
  case('TEC','tec')
     write(iFile,*)'TITLE = "BATSRUS: Solar Particle Data at SP_Time = ',&
          SP_Time,'" '
     write(iFile,*)trim(SP_NameUnitTec)
     write(iFile,*)'ZONE T = "SP Image"',', I=',(nP-1)/nPSkip+1,', J=',(nX-1)/nXSkip+1,', K=1, F=POINT'
     !\
     ! Write data values in the file.
     !/
     do iX=1,nX,nXSkip
        !        Location = real(iX)    !in MKS units
        Location = sqrt(sum(X_DI(:,iX)**2))       !in MKS units
        Location = Location/Rsun                  !in units of Rsun
        do iLnP=1,nP,nPSkip
           Momentum = PInjection*exp(real(iLnP)* &
                DeltaLnP)                         !in MKS units.
           Energy   = momentum_to_kinetic_energy(&
                Momentum,NameParticle)            !in MKS units.
           call sp_set_plotvar
           write(iFile,'(50(E14.6))') Location*SP_Unit_X,&
                Energy*SP_Unit_E,PlotVar_III(iLnP,iX,1:nVar)
        end do
     end do
  case('IDL','idl')
     !\
     ! Description of file contains units, physics and
     ! dimension.
     !/
     write(iFile,'(a)')'BATSRUS: Solar Particle Data'
     !\
     ! *2* in the next line means two-dimensional plot, *1* in
     ! the next line is for a 2D-cut; in other words,
     ! one-dimension is left out.
     !/
     write(iFile,'(i7,1pe13.5,3i3)')&
          SP_iPlot,SP_Time,2,nVar+2,nVar+2
     !\
     ! Write grid-size.
     !/
     write(iFile,'(2i4)')nP,nX
     !\
     ! Write coordinate, variable, and equation parameter
     ! names.
     !/
     write(iFile,'(a)')SP_NameAllVars
     !\
     ! Write data values in the file.
     !/
     do iX=1,nX
        !        Location = real(iX)    !in MKS units
        Location = sqrt(sum(X_DI(:,iX)**2))       !in MKS units
        Location = Location/Rsun                  !in units of Rsun
        do iLnP=1,nP
           Momentum = PInjection*exp(real(iLnP)* &
                DeltaLnP)                         !in MKS units.
           Energy   = momentum_to_kinetic_energy(&
                Momentum,NameParticle)            !in MKS units.
           call sp_set_plotvar
           write(iFile,'(50(1pe13.5))')& 
                Location*SP_Unit_X,&
                Energy*SP_Unit_E,PlotVar_III(iLnP,iX,1:nVar)
        end do
     end do
  end select
  close(iFile)
  !\
  ! Deallocate array for plot variables.
  !/
  deallocate(PlotVar_III)
  if(DoOutputGamma)then
     iFile=io_unit_new()
     !\
     ! Construct the name of output file.
     !/
     write(NameFile,'(a,i4.4,a,i4.4,a)')trim(SP_DirOut)//'Wave_ix',&
          iXOutputGamma,'_n',SP_iPlot,SP_NameExtension
     write(iStdout,*)prefix,' '
     write(iStdout,*)prefix,'Opening File '//trim(NameFile)//' to Save wave spectra'//&
       ' at SP_Time =',SP_Time
     write(iStdout,*)prefix,' '
     !\
     ! Open file unit and write plot file.
     !/
     open(iFile,file=trim(NameFile),status='replace',form='formatted')
     write(iFile,*)'Waves at SP_Time = ',&
          SP_Time,'" '
     write(iFile,*) 'VARIABLES = "k", "I+", "I-", "gamma" '

      do iLnK=1,nP
         write(iFile,*) B_I(iXOutputGamma)*Kmin*exp(real(iLnK-1)*DeltaLnK),&
              IPlus_IX(iLnK,iXOutputGamma),IMinus_IX(iLnK,iXOutputGamma),&
              Gamma_I(iLnK)
      end do
      close(iFile)
   end if
        
Contains
  !============================================================================!
  subroutine sp_set_plotvar
    implicit none
    !--------------------------------------------------------------------------!
    do iVar=1,nVar
       SP_NamePlotVar = SP_NamePlotVar_I(iVar)
       select case(SP_NamePlotVar)
       case('CDF','cdf')
          ! Canonical distribution function: f in units of
          ! [s^3 kg^-^3 m^-^6].
          PlotVar_III(iLnP,iX,iVar) = F_II(iLnP,iX)
       case('DEF','def')
          ! Differential energy flux: F = f*p^2 in units of
          ! [s^-^1 m^-^2 str^-^1 keV^-^1].
          PlotVar_III(iLnP,iX,iVar) = (F_II(iLnP,iX)*&
               Momentum)*(Momentum*energy_in(NameEnergyUnit))
       case('IND','ind')
          if(F_II(iLnP-1,iX)<=0.0)then
             write(*,*)'In plot writing: negative value ',F_II(iLnP-1,iX),&
                  ' at iLnP-1=',iLnP-1,' ix=', iX
             stop
          end if
          PlotVar_III(iLnP,iX,iVar) = &
               log(F_II(iLnP  ,iX)/   &
                   F_II(iLnP-1,iX))/  &
                   DeltaLnP+2.0
       case default
          write(iStdOut,*)prefix//' unknown plot var '//SP_NamePlotVar
          stop
       end select
    end do
  end subroutine sp_set_plotvar
  !============================================================================!
  subroutine sp_get_tec_units
    implicit none
    !--------------------------------------------------------------------------!
    !\
    ! Construct the string with variable names.
    !/
    write(SP_NameUnitTec,'(a)') 'VARIABLES = '
    write(SP_NameUnitTec,'(a)') &
         trim(SP_NameUnitTec)//'"X '//   &
         trim(SP_NameUnitTec_X)
    write(SP_NameUnitTec,'(a)') &
         trim(SP_NameUnitTec)//'", "E '//&
         trim(SP_NameUnitTec_E)
    do iVar=1,nVar
       write(SP_NameUnitTec,'(a)') &
            trim(SP_NameUnitTec)//'", "'
       SP_NamePlotVar = SP_NamePlotVar_I(iVar)
       select case(SP_NamePlotVar)
       case('CDF','cdf')
          write(SP_NameUnitTec,'(a)') & 
               trim(SP_NameUnitTec)//'f'//' '//&
               trim(SP_NameUnitTec_CDF)
       case('DEF','def')
          write(SP_NameUnitTec,'(a)') & 
               trim(SP_NameUnitTec)//'F'//' '//&
               trim(SP_NameUnitTec_DEF)
       case('IND','ind')
          write(SP_NameUnitTec,'(a)') & 
               trim(SP_NameUnitTec)//'Gamma'//' '//&
               trim(SP_NameUnitTec_G)
       case default
          write(SP_NameUnitTec,'(a)') & 
               trim(SP_NameUnitTec)//'Default'
       end select
    end do
    write(SP_NameUnitTec,'(a)') &
         trim(SP_NameUnitTec)//'" '
  end subroutine sp_get_tec_units
  !============================================================================!
  subroutine sp_get_idl_units
    implicit none
    !--------------------------------------------------------------------------!
    !\
    ! Construct the string with variable names.
    !/
    write(SP_NameUnitIdl,'(a)')       & 
         trim(SP_NameUnitIdl_X)//' '//&
         trim(SP_NameUnitIdl_E)
    do iVar=1,nVar
       SP_NamePlotVar = SP_NamePlotVar_I(iVar)
       select case(SP_NamePlotVar)
       case('CDF','cdf')
          write(SP_NameUnitIdl,'(a)') & 
               trim(SP_NameUnitIdl)//' '//&
               trim(SP_NameUnitIdl_CDF)
       case('DEF','def')
          write(SP_NameUnitIdl,'(a)') & 
               trim(SP_NameUnitIdl)//' '//&
               trim(SP_NameUnitIdl_DEF)  
       case('IND','ind')
          write(SP_NameUnitIdl,'(a)') & 
               trim(SP_NameUnitIdl)//' '//&
               trim(SP_NameUnitIdl_G)  
       case default
          write(SP_NameUnitIdl,'(a)') & 
               trim(SP_NameUnitIdl)//'" Default"'
       end select
    end do
  end subroutine sp_get_idl_units
  
end subroutine write_plotfile_SP
!=============================================================================!
subroutine sp_split_str(str,nmax,strarr,n)
  implicit none
  !---------------------------------------------------------------------------!
  character(len=*),intent(in):: str
  integer,intent(in):: nmax
  character(len=10),intent(out):: strarr(nmax)
  integer, intent(out):: n
  !---------------------------------------------------------------------------!
  character(len=100):: s
  integer:: i,l
  !---------------------------------------------------------------------------!
  n=0
  l=len_trim(str)
  s=str(1:l)
  !  s(l+1:l+1)=' '
  do
     ! Check leading spaces
     i=1
     do while(s(i:i)==' ' .and. i<=l)
        i=i+1
     end do
     if(i>l) EXIT       ! All spaces.
     if(i>1)s=s(i:l)    ! Delete leading spaces.
     i=index(s,' ')     ! Find end of first word.
     n=n+1              ! Put word into strarr.
     strarr(n)=s(1:i-1)
     s=s(i+1:l)         ! Delete word and 1 space from string.
     if(n==nmax)exit
  end do
end subroutine sp_split_str
!=============================================================================!
