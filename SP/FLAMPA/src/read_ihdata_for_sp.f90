
!============================================================================!
subroutine read_ihdata_for_sp(NameFileIn,nLineHeader,nLineBad)
  use SP_ModMain
  use CON_world, ONLY: CON_stop
  use ModConst,  ONLY: cProtonMass,Rsun
  implicit none
  include 'stdout.h'
  !--------------------------------------------------------------------------!
  character(LEN=*),intent(in):: NameFileIn
  integer,intent(in):: nLineHeader
  integer,intent(in):: nLineBad
  !--------------------------------------------------------------------------!
  character(LEN=*),parameter:: IH_dir='./IO_IH/IHDATA/'
  character(LEN=80):: NameFile
  character(LEN=80):: TextHeader
  integer:: iLine,nPoint,iFile,iError=0
  real,dimension(11):: State_V
  !--------------------------------------------------------------------------!
  NameFile=IH_dir//trim(NameFileIn)
  write(iStdout,*)prefix,'Reading IH data from ',NameFile
  !\
  ! First reading of file sets the value of nX::
  !/
  call get_io_unit_new(iFile)
  open(iFile,file=NameFile,status='old',iostat=iError)
  if (nLineHeader>0) then
     do iLine=1,nLineHeader
        read(iFile,'(a)') TextHeader
        read(TextHeader(20:40),*) TimeToRead
     end do
  end if
  nPoint=1
  do
     read(iFile,*,iostat=iError) State_V(:)
     if(iError/=0)exit
     nPoint=nPoint+1
  end do
  close(iFile)
  nX=nPoint-nLineBad
  ! Last nLineBad lines of data that are bad!
  if(.not.allocated(DInner_I))call SP_allocate
  nX=min(nX,ubound(DInner_I,1))
  
  !\
  ! Second reading of file with data from IH_::
  !/
  call get_io_unit_new(iFile)
  open(iFile,file=NameFile,status='old',iostat=iError)
  if (nLineHeader>0) then
     do iLine=1,nLineHeader
        read(iFile,'(a)') TextHeader
        write(iStdout,*)prefix,'at ',TextHeader
     end do
  end if
  RhoSmoothOld_I(1:nX)= RhoSmooth_I(1:nX)
  iShockOld=iShock
  do iLine=1,nX
     read(iFile,*,iostat=iError) State_V(:)
     X_DI(:,iLine)=State_V(1:3)*Rsun           ! in [m]
     RhoSmooth_I(iLine) =State_V(4)/cProtonMass      ! in [m^-3]
     T_I(iLine)   =State_V(11)/RhoSmooth_I(iLine)/&  ! in NameEnergyUnit [KeV]
          energy_in(NameEnergyUnit)
     B_I(iLine)   =sqrt(sum(State_V(8:10)**2)) ! in [Tesla]
     U_I(iLine)   =sqrt(sum(State_V(5:7)**2))  ! in [m/s]
  end do
  close(iFile)
  !\
  ! Write min/max values of MHD variables being read from file::
  !/
  write(iStdout,*)prefix//'nX is set to: ',nX 
  write(iStdout,*)prefix,'min/max value of number density in [m^-3]    = ',&
       minval(Rho_I),maxval(Rho_I)
  write(iStdout,*)prefix,'min/max value of temperature in ['//             &
       trim(NameEnergyUnit)//']       = ',minval(T_I),maxval(T_I)
  write(iStdout,*)prefix,'min/max value of radial distance in [Rsun]   = ',&
       sqrt(minval(X_DI(1,:)**2+X_DI(2,:)**2+X_DI(3,:)**2))/Rsun,&
       sqrt(maxval(X_DI(1,:)**2+X_DI(2,:)**2+X_DI(3,:)**2))/Rsun
  write(iStdout,*)prefix,'min/max value of field strength in [Tesla]   = ',&
       minval(B_I),maxval(B_I)
  iShock=i_shock(RhoSmooth_I(1:nX))
end subroutine read_ihdata_for_sp
!=============================================================================
