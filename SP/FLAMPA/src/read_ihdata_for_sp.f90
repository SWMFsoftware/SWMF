
!============================================================================!
subroutine read_ihdata_for_sp(nLineHeader,nLineBad)
  use SP_ModMain
  use ModIOUnit
  implicit none
  !--------------------------------------------------------------------------!
  integer,intent(in):: nLineHeader
  integer,intent(in):: nLineBad
  !--------------------------------------------------------------------------!
  character(LEN=*),parameter:: IH_dir='./SP/MHDATA/'
  character(LEN=80):: NameFile
  character(LEN=80):: TextHeader
  integer:: iLine,nPoint,iFile,iError=0
  real,dimension(11):: State_V
  !--------------------------------------------------------------------------!
  write(NameFile,'(a,i4.4,a)')IH_dir//'mhdata_',iDataSet,'.dat'
  write(iStdout,*)prefix,'Reading IH data from ',trim(NameFile)
  !\
  ! First reading of file sets the value of nX::
  !/
  iFile=io_unit_new()
  open(iFile,file=NameFile,status='old',iostat=iError)
  if(iError>0)call CON_stop('File does not exist')
  if (nLineHeader>0) then
     do iLine=1,nLineHeader
        read(iFile,'(a)') TextHeader
        read(TextHeader(20:40),*) DataInputTime
     end do
  end if
  nPoint=0
  do
     read(iFile,*,iostat=iError) State_V(:)
     if(iError/=0)exit
     nPoint=nPoint+1
  end do
  close(iFile)
  nX=nPoint-nLineBad
  write(iStdOut,*)prefix,' nX=',nX
  ! Last nLineBad lines of data that are bad!
  if(.not.allocated(DInner_I))call SP_allocate
  nX=min(nX,ubound(DInner_I,1))
  
  !\
  ! Second reading of file with data from IH_::
  !/
  iFile=io_unit_new()
  open(iFile,file=NameFile,status='old',iostat=iError)
  if (nLineHeader>0) then
     do iLine=1,nLineHeader
        read(iFile,'(a)') TextHeader
        write(iStdout,*)prefix,'at ',TextHeader
     end do
  end if

  do iLine=1,nX
     read(iFile,*)X_DI(:,iLine), State_VI(:,iLine)
  end do
  close(iFile)
  call mh_transform_for_flampa
end subroutine read_ihdata_for_sp
subroutine mh_transform_for_flampa
  use ModConst
  use SP_ModMain
  implicit none
  integer::iPoint

  RhoSmoothOld_I(1:nX)= RhoSmooth_I(1:nX)
  BSmoothOld_I(1:nX)= BSmooth_I(1:nX)
  iShockOld=iShock

  do iPoint=1,nX
     X_DI(:,iPoint)=X_DI(1:3,iPoint)*Rsun           ! in [m]
     RhoSmooth_I(iPoint) =State_VI(1,iPoint)/cProtonMass      ! in [m^-3]
     T_I(iPoint)   =State_VI(8,iPoint)/RhoSmooth_I(iPoint)/&    ! in NameEnergyUnit [KeV]
          energy_in(NameEnergyUnit)
     BSmooth_I(iPoint)   =sqrt(sum(State_VI(5:7,iPoint)**2)) ! in [Tesla]
     U_I(iPoint)   =sqrt(sum(State_VI(2:4,iPoint)**2))  ! in [m/s]
  end do
  !\
  ! Write min/max values of MHD variables being read from file::
  !/
  write(iStdout,*)prefix//'nX is set to: ',nX 
  write(iStdout,*)prefix,'min/max value of number density in [m^-3]    = ',&
       minval(RhoSmooth_I),maxval(RhoSmooth_I)
  write(iStdout,*)prefix,'min/max value of temperature in ['//             &
       trim(NameEnergyUnit)//']       = ',minval(T_I),maxval(T_I)
  write(iStdout,*)prefix,'min/max value of radial distance in [Rsun]   = ',&
       sqrt(minval(X_DI(1,:)**2+X_DI(2,:)**2+X_DI(3,:)**2))/Rsun,&
       sqrt(maxval(X_DI(1,:)**2+X_DI(2,:)**2+X_DI(3,:)**2))/Rsun
  write(iStdout,*)prefix,'min/max value of field strength in [Tesla]   = ',&
       minval(BSmooth_I),maxval(BSmooth_I)
  iShock=i_shock(RhoSmooth_I(1:nX))
end subroutine mh_transform_for_flampa
!=============================================================================
