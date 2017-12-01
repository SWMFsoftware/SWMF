!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                                                                |
!/////////////////////////////////////////////////////////////////
program PostIONO
  implicit none
  
  !\
  ! PostProcessing variables.
  !/
  integer :: ptLength,pvLength
  integer :: nvars,nTheta,nPhi
  integer :: nYear,nMonth,nDay,nHour,nMinute,nSecond,nMillisecond
  integer :: i,j,n,n_step
  
  real :: Theta,Phi, X,Y,Z, DegToRad, TimeSimulation, ttilt,ptilt, Radius

  character (len=1000) :: PlotTitle, PlotVars, ZoneTitle
  character (len=1000) :: text,text2
  character (len=22)  :: textNandT
  character (len=80) :: stmp

  character (LEN=4) :: TimeH4
  character (LEN=2) :: TimeM2,TimeS2
  !--------------------------------------------------------------------------
  Radius = (6378.+100.)/6378.

  ! Setup Degree to Radian conversion factor
  DegToRad = (2.*asin(1.))/180.
  
  !\
  ! Begin processing
  !/
  write(*,*) ' '
  write(*,*) 'PostIONO: BATSRUS Post-Processing Package for converting'
  write(*,*) '          ionosphere files from idl format to tecplot format.'
  write(*,*) '    -Written by  Darren De Zeeuw          '
  write(*,*) ' '
  
  ! Open file
  text='inputfile.DAT'
  call FixChar(n)
  open(unit=50,file=text(1:n),status='old',form='formatted')
  rewind(50)

  ! Read title
  read(50,'(a)') text
  call FixChar(n)
  PlotTitle=text(1:n)
  PlotTitle(n+1:n+1)='='
  ptLength=n+1
  read(50,'(a)') text
  call FixChar(n)
  PlotTitle(ptLength+1:ptLength+n)=text(1:n)
  ptLength=ptLength+n

  ! Skip 2 lines
  read(50,'(a)') text
  read(50,'(a)') text

  ! Read nvars, nTheta, and nPhi
  read(50,'(i6,a)') nvars,text
  read(50,'(i6,a)') nTheta,text
  read(50,'(i6,a)') nPhi,text
  write(*,*) '    nvars =',nvars,' nTheta =',nTheta,' nPhi =',nPhi

  ! Skip 2 lines
  read(50,'(a)') text
  read(50,'(a)') text

  ! Read variable names
  PlotVars='VARIABLES = "X [R]", "Y [R]", "Z [R]", '
  pvLength=39
  do i=1,nvars
     read(50,'(a6,a)') text2,text
     call FixChar(n)
     if(i>1)then
        PlotVars(pvLength+1:pvLength+2)=', '
        pvLength=pvLength+2
     end if

     PlotVars(pvLength+1:pvLength+1)='"'
     pvLength=pvLength+1

     PlotVars(pvLength+1:pvLength+n)=text(1:n)
     pvLength=pvLength+n

     PlotVars(pvLength+1:pvLength+1)='"'
     pvLength=pvLength+1
  end do

  ! Skip 2 lines
  read(50,'(a)') text
  read(50,'(a)') text

  ! Read date
  read(50,'(i6,a)') nYear,text
  read(50,'(i6,a)') nMonth,text
  read(50,'(i6,a)') nDay,text
  read(50,'(i6,a)') nHour,text
  read(50,'(i6,a)') nMinute,text
  read(50,'(i6,a)') nSecond,text
  read(50,'(i6,a)') nMillisecond,text
  write(*,*) '     time =',nYear,nMonth,nDay,nHour,nMinute,nSecond,nMillisecond

  ! Skip 2 lines
  read(50,'(a)') text
  read(50,'(a)') text

  ! Read simulation time
  read(50, '(i13,a)')     n_step, text
  read(50, '(1pe13.5,a)') TimeSimulation, text

  ! Skip 2 lines
  read(50,'(a)') text
  read(50,'(a)') text

  ! Read tilts
  read(50, '(f7.2,a)') ttilt, text
  read(50, '(f7.2,a)') ptilt, text

  ! Close IDL file
  close(50)

  ! Create zone name
  write(TimeH4,'(i4.4)') &
     int(                           TimeSimulation/3600.)
  write(TimeM2,'(i2.2)') &
     int((TimeSimulation-(3600.*int(TimeSimulation/3600.)))/60.)
  write(TimeS2,'(i2.2)') &
     int( TimeSimulation-(  60.*int(TimeSimulation/  60.)))
  write(textNandT,'(a,i7.7,a)') &
       "N=",n_step," T="//TimeH4//":"//TimeM2//":"//TimeS2

  write(ZoneTitle,'(a,i8,a,i8,a)') 'ZONE T="IonN '//textNandT//'"'//', I=',nTheta,', J=',nPhi,', K=1, F=POINT'

  ! Open and write the North datafile header
  text='outputfileN.DAT'
  call FixChar(n)
  open (unit=51,file=text(1:n),status='unknown',form='formatted')
  rewind(51)
  write(51,'(a)') PlotTitle(1:ptLength)
  write(51,'(a)') PlotVars(1:pvLength)
  write(51,'(a)') ZoneTitle(1:len_trim(ZoneTitle))

  ! Write AUXDATA
  write(51,'(a)') 'AUXDATA HEMISPHERE="North"'
  write(stmp,'(i4.4,a,5(i2.2,a),i3.3)') &
       nYear,'-',nMonth,'-',nDay,' ',nHour,':',nMinute,':',nSecond,'.',nMillisecond
  write(51,'(a,a,a)') 'AUXDATA TIMEEVENT="',trim(adjustl(stmp)),'"'
  write(stmp,'(a)') " T="//TimeH4//":"//TimeM2//":"//TimeS2
  write(51,'(a,a,a)') 'AUXDATA TIMESIM="',trim(adjustl(stmp)),'"'
  write(stmp,'(a)') " T="//TimeH4//":"//TimeM2
  write(51,'(a,a,a)') 'AUXDATA TIMESIMSHORT="',trim(adjustl(stmp)),'"'
  write(stmp,'(f7.2)') ttilt
  write(51,'(a,a,a)') 'AUXDATA THETATILTDEG="',trim(adjustl(stmp)),'"'
  write(stmp,'(f7.2)') ptilt
  write(51,'(a,a,a)') 'AUXDATA PHITILTDEG="',trim(adjustl(stmp)),'"'

  ! Open and read until 'BEGIN NORTHERN HEMISPHERE'
  text='inputfile.DAT'
  call FixChar(n)
  open(unit=50,file=text(1:n),status='old',form='formatted')
  rewind(50)
  do
     read(50,'(a)') text
     if(index(text,'BEGIN NORTHERN HEMISPHERE') > 0) EXIT
  end do

  ! Read data and write to North datafile
  do i=1,nTheta
     do j=1,nPhi
        read(50,'(a)') text
        call FixCharEnd(n)
        read(text,*) Theta,Phi
        Phi=Phi*DegToRad
        Theta=Theta*DegToRad
        X=sin(Theta)*cos(Phi)
        Y=sin(Theta)*sin(Phi)
        Z=cos(Theta)
        write(51,'(3E14.6,a)') X,Y,Z,text(1:n)
     end do
  end do

  ! Close IDL file and North datafile
  close(50)
  close(51)

  ! Open and write the South datafile header
  write(ZoneTitle,'(a,i8,a,i8,a)') 'ZONE T="IonS '//textNandT//'"'//', I=',nTheta,', J=',nPhi,', K=1, F=POINT'
  text='outputfileS.DAT'
  call FixChar(n)
  open (unit=51,file=text(1:n),status='unknown',form='formatted')
  rewind(51)
  write(51,'(a)') PlotTitle(1:ptLength)
  write(51,'(a)') PlotVars(1:pvLength)
  write(51,'(a)') ZoneTitle(1:len_trim(ZoneTitle))

  ! Write AUXDATA
  write(51,'(a)') 'AUXDATA HEMISPHERE="South"'
  write(stmp,'(i4.4,a,5(i2.2,a),i3.3)') &
       nYear,'-',nMonth,'-',nDay,' ',nHour,':',nMinute,':',nSecond,'.',nMillisecond
  write(51,'(a,a,a)') 'AUXDATA TIMEEVENT="',trim(adjustl(stmp)),'"'
  write(stmp,'(a)') " T="//TimeH4//":"//TimeM2//":"//TimeS2
  write(51,'(a,a,a)') 'AUXDATA TIMESIM="',trim(adjustl(stmp)),'"'
  write(stmp,'(a)') " T="//TimeH4//":"//TimeM2
  write(51,'(a,a,a)') 'AUXDATA TIMESIMSHORT="',trim(adjustl(stmp)),'"'
  write(stmp,'(f7.2)') ttilt
  write(51,'(a,a,a)') 'AUXDATA THETATILTDEG="',trim(adjustl(stmp)),'"'
  write(stmp,'(f7.2)') ptilt
  write(51,'(a,a,a)') 'AUXDATA PHITILTDEG="',trim(adjustl(stmp)),'"'

  ! Open and read until 'BEGIN SOUTHERN HEMISPHERE'
  text='inputfile.DAT'
  call FixChar(n)
  open(unit=50,file=text(1:n),status='old',form='formatted')
  rewind(50)
  do
     read(50,'(a)') text
     if(index(text,'BEGIN SOUTHERN HEMISPHERE') > 0) EXIT
  end do

  ! Read data and write to South datafile
  do i=1,nTheta
     do j=1,nPhi
        read(50,'(a)') text
        call FixCharEnd(n)
        read(text,*) Theta,Phi
        Phi=Phi*DegToRad
        Theta=(180.-Theta)*DegToRad
        X=Radius*sin(Theta)*cos(Phi)
        Y=Radius*sin(Theta)*sin(Phi)
        Z=Radius*cos(Theta)
        write(51,'(3E14.6,a)') X,Y,Z,text(1:n)
     end do
  end do

  ! Close IDL file and South datafile
  close(50)
  close(51)

  ! Execution completed.
  write(*,*) ' '
  write(*,*) 'PostIONO: Execution completed.'
  write(*,*) ' '

contains
  
  !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  ! Character String Fix Subroutine                                |
  !/////////////////////////////////////////////////////////////////
  
  subroutine FixChar(LenTot)
    integer, intent(inout) :: LenTot
    integer :: k
    
    if (text.eq.' ') then
       LenTot=0
       text=' '
    else
       LenTot=len(text)
       k=0
       do while (text(1:1).eq.' '.and.k.le.LenTot)
          k=k+1
          text=text(2:LenTot)
          text(LenTot:LenTot)=' '
       end do
       k=LenTot
       do while (text(k:k).eq.' '.and.k.ge.1)
          k=k-1
          LenTot=LenTot-1
       end do
    end if
    
  end subroutine FixChar
  
  subroutine FixCharEnd(LenTot)
    integer, intent(inout) :: LenTot
    integer :: k
    
    if (text.eq.' ') then
       LenTot=0
       text=' '
    else
       LenTot=len(text)
       k=LenTot
       do while (text(k:k).eq.' '.and.k.ge.1)
          k=k-1
          LenTot=LenTot-1
       end do
    end if
    
  end subroutine FixCharEnd
  
end program PostIONO
