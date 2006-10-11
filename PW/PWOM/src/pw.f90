program pw

  use Mod_PW
  use ModPass

!******************************************************************************
! Initiallize MPI and get number of processors and rank of given processor
!******************************************************************************
  include 'mpif.h'
  integer errcode,comm,iProc,nProc
  call MPI_INIT(errcode)
  comm = MPI_COMM_WORLD
  
  call MPI_COMM_RANK(comm,iProc,errcode)
  call MPI_COMM_SIZE(comm,nProc,errcode)
  write(*,*) iProc,nProc

!******************************************************************************
!  Set the number of fieldlines that each processor solves for
!******************************************************************************
  if (iproc .lt. mod(maxLine,nproc)) then
     nLine=int(ceiling(real(maxline)/real(nproc)))
  else
     nLine=int(floor(real(maxline)/real(nproc)))
  endif
  nLine=int(nLine)
  
  

!******************************************************************************
!  Define file names and unit numbers, and open for reading and writing.
!******************************************************************************
  NameInput          = 'pw.input'
  NameOutput         = 'log.out'
  !NameGraphics       = 'plots.out'
  NameSourceGraphics = 'plot_sources.out'
!  NameRestart        = 'restart.out'
  NameCollision      = 'plots_collision.out'
!  NameRestartIn      = 'restart.in'
  NamePhiNorth  = 'North.dat'
  NamePhiSouth  = 'South.dat'




  iUnitInput         = 11 
  iUnitOutput        = 16  
  iUnitSourceGraphics= 18
!  iUnitRestart       = 7
  iUnitCollision     = 17  
!  iUnitRestartIn     = 8
  iUnitSouth         = 101
  iUnitNorth         = 102

  do iLine=1,nLine
     if (iproc .lt. mod(maxLine,nproc)) then
        iLineGlobal(iLine)=&
             iproc*ceiling(real(maxline)/real(nproc))+iLine
     else
        iLineGlobal(iLine)=&
             (mod(maxLine,nproc))*ceiling(real(maxline)/real(nproc)) &
             + ((iproc)-mod(maxLine,nproc))                        &
             *floor(real(maxline)/real(nproc))+iLine
     endif
     write(*,*) 'iLineGlobal',iLineGlobal(iLine)
     
     write(NameRestartIn(iLine),"(a,i4.4,a)") &
          'PW/restartIN/restart_iline',iLineGlobal(iLine),'.dat'
     write(NameRestart(iLine),"(a,i4.4,a)") &
          'PW/restartOUT/restart_iline',iLineGlobal(iLine),'.dat'

     write(NameGraphics(iLine),"(a,i4.4,a)") &
          'PW/plots/plots_iline',iLineGlobal(iLine),'.out'

     iUnitRestartIn(iLine) = 200+iLineGlobal(iLine)
     iUnitRestart(iLine)   = 300+iLineGlobal(iLine)
     iUnitGraphics(iLine)  = 400+iLineGlobal(iLine)

     OPEN(iUnitRestart(iLine) ,FILE=NameRestart(iLine))
     OPEN(iUnitGraphics(iLine),FILE=NameGraphics(iLine))
  enddo
  
   

  
  OPEN(iUnitInput,         FILE=NameInput)
  OPEN(UNIT=iUnitOutput,   FILE=NameOutput)

!  OPEN(iUnitSourceGraphics,FILE=NameSourceGraphics)

!  OPEN(iUnitCollision,     FILE=NameCollision)
  open(iUnitNorth,file=NamePhiNorth)  


!******************************************************************************
! Read the input file
!******************************************************************************
 
  READ(iUnitInput,*) TMAX
  WRITE(iUnitOutput,*) TMAX
  
  READ(iUnitInput,*) DToutput
  WRITE(iUnitOutput,*) DToutput
  
  READ(iUnitInput,*) TypeSolver
  WRITE(iUnitOutput,*) TypeSolver
  
  READ(iUnitInput,*) IsImplicit
  WRITE(iUnitOutput,*) IsImplicit
  
  read(iUnitInput,*) IsRestart
  if (IsRestart) then
     write(*,*) 'Is Restart', IsRestart
  endif
  read(iUnitInput,*) IsVariableDt
  if (IsVariableDt) then
     write(*,*) 'IsVariableDT', IsVariableDt
  endif
  READ(iUnitInput,*)   DTpolarwind
  WRITE(iUnitOutput,*) DTpolarwind
  
  READ(iUnitInput,*)   IsMoveFluxTube 
  READ(iUnitInput,*)   IsUseJr
  READ(iUnitInput,*)   IsCentrifugal

!******************************************************************************
!  Read the restart file
!******************************************************************************
  nDim= 390
  if(IsRestart)then
         
     do iLine=1,nLine
        OPEN(UNIT=iUnitRestartIn(iLine),FILE=NameRestartIn(iLine),STATUS='OLD')
        READ (iUnitRestartIn(iLine),2001) TIME,DDT1,NS
        READ (iUnitRestartIn(iLine),*)    GeoMagLat(iLine),GeoMagLon(iLine)
        
        READ (iUnitRestartIn(iLine),2002) &
             (XXX,uOxyg(i,iLine),pOxyg(i,iLine),dOxyg(i,iLine),TOxyg(i,iLine),&
             i=1,NDIM)
        READ (iUnitRestartIn(iLine),2002) &
             (XXX,uHel(i,iLine),pHel(i,iLine),dHel(i,iLine),THel(i,iLine),    &
             i=1,NDIM)
        READ (iUnitRestartIn(iLine),2002) &
             (XXX,uHyd(i,iLine),pHyd(i,iLine),dHyd(i,iLine),THyd(i,iLine),    &
             i=1,NDIM)
        READ (iUnitRestartIn(iLine),2002) &
             (XXX,uElect(i,iLine),pElect(i,iLine),dElect(i,iLine),            &
             TElect(i,iLine),i=1,NDIM)
        CLOSE(UNIT=iUnitRestartIn(iLine))
     enddo

2001 FORMAT(2(1PE16.6),I10)
     ! 2002    FORMAT(5(1PE16.6))
2002 FORMAT(5(1PE25.16))
  Else
     Time=0.0
  endif
  
!****************************************************************************
! Use Get_GITM to bring in neutral atmosphere from GITM
!****************************************************************************  
  !call GetNeutralData

!******************************************************************************
!  Set parameters for reading in potential and time of simulation
!******************************************************************************

!  IsMoveFluxTube = .false.


!  IsUseJr = .false.

!  IsCentrifugal =.true.


  Dtheta  = 0.0242
  Dphi    = 0.0245
  
  Dt      =    50.0
!  Dt      =    0.2
!  Dt      =    Tmax
  !maxTime = 10000.0
  maxTime = Tmax
  Time    =     0.0

  
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        read(unit=iUnitNorth,fmt='(6(1PE13.5))') &
             Theta_G(iPhi,iTheta),Phi_G(iPhi,iTheta),SigmaH_G(iPhi,iTheta),&
             SigmaP_G(iPhi,iTheta),Jr_G(iPhi,iTheta),Potential_G(iPhi,iTheta)

     enddo
  enddo

  

!******************************************************************************
!  Change angles to radians
!******************************************************************************
  Theta_G(:,:) = Theta_G(:,:)*3.1416/180.0
  Phi_G  (:,:) = Phi_G  (:,:)*3.1416/180.0

!******************************************************************************
!  Convert potential from kilovolts to Volts
!******************************************************************************
  Potential_G(:,:) = Potential_G(:,:)*1.0e3
!******************************************************************************
!  Calc Bfield components
!******************************************************************************
  rPlanet        = 6378000.0
  rLowerBoundary = rPlanet+110.0e3
  MagMoment      = 7.84e15
  Bcoef          = MagMoment/(rLowerBoundary)**3 
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        Br_G(iPhi,iTheta)     = -Bcoef*2.0*sin(Theta_G(iPhi,iTheta))
        Btheta_G(iPhi,iTheta) = Bcoef*cos(Theta_G(iPhi,iTheta))
     enddo
  enddo
  
  BmagnitudeSquared_G(:,:) = Br_G(:,:)**2 + Btheta_G(:,:)**2
!******************************************************************************
!  Fill Ghost cells
!******************************************************************************
  Theta_G    (nPhi+1,:) = Theta_G    (2,:) 
  Phi_G      (nPhi+1,:) = Phi_G      (2,:)
  SigmaH_G   (nPhi+1,:) = SigmaH_G   (2,:)
  SigmaP_G   (nPhi+1,:) = SigmaP_G   (2,:) 
  Jr_G       (nPhi+1,:) = Jr_G       (2,:)
  Potential_G(nPhi+1,:) = Potential_G(2,:)

  Theta_G    (0,:) = Theta_G    (nPhi-1,:) 
  Phi_G      (0,:) = Phi_G      (nPhi-1,:)
  SigmaH_G   (0,:) = SigmaH_G   (nPhi-1,:)
  SigmaP_G   (0,:) = SigmaP_G   (nPhi-1,:) 
  Jr_G       (0,:) = Jr_G       (nPhi-1,:)
  Potential_G(0,:) = Potential_G(nPhi-1,:)

  Theta_G    (:,nTheta+1) = Theta_G    (:,nTheta) 
  Phi_G      (:,nTheta+1) = Phi_G      (:,nTheta)
  SigmaH_G   (:,nTheta+1) = SigmaH_G   (:,nTheta)
  SigmaP_G   (:,nTheta+1) = SigmaP_G   (:,nTheta) 
  Jr_G       (:,nTheta+1) = Jr_G       (:,nTheta)
  Potential_G(:,nTheta+1) = Potential_G(:,nTheta)

  do iPhi=1,nPhi
     Theta_G    (iPhi,0) = Theta_G    (mod(iPhi+128,nPhi-1),2) 
     Phi_G      (iPhi,0) = Phi_G      (mod(iPhi+128,nPhi-1),2)
     SigmaH_G   (iPhi,0) = SigmaH_G   (mod(iPhi+128,nPhi-1),2)
     SigmaP_G   (iPhi,0) = SigmaP_G   (mod(iPhi+128,nPhi-1),2) 
     Jr_G       (iPhi,0) = Jr_G       (mod(iPhi+128,nPhi-1),2)
     Potential_G(iPhi,0) = Potential_G(mod(iPhi+128,nPhi-1),2)
  enddo



!******************************************************************************
!  Calc electric field from E=-grad Potential
!******************************************************************************
 
  
  do iPhi=1,nPhi
     do iTheta=1,nTheta
        Etheta(iPhi,iTheta) = &
             (Potential_G(iPhi,iTheta)-Potential_G(iPhi,iTheta+1))&
             /((rLowerBoundary)*Dtheta)
        if (iTheta == 1) then
           Ephi(iPhi,iTheta) = 0.0
        else
           Ephi(iPhi,iTheta)   = &
                (Potential_G(iPhi-1,iTheta)-Potential_G(iPhi+1,iTheta))&
                / ((rLowerBoundary)*2.0*Dphi  &
                * sin(Theta_G(iPhi,iTheta)))
        endif
    enddo
  enddo
  Er(:,:) = 0.0
!******************************************************************************
!  Calc VelocityExB drift from E and B
!******************************************************************************
  

  do iPhi=1,nPhi
     do iTheta=1,nTheta
        VelocityExBtheta(iPhi,iTheta) = &
             (Ephi(iPhi,iTheta)*Br_G(iPhi,iTheta)) &
             / BmagnitudeSquared_G(iPhi,iTheta)
        
        VelocityExBphi(iPhi,iTheta)   = &
             (-Etheta(iPhi,iTheta)*Br_G(iPhi,iTheta)) &
             / BmagnitudeSquared_G(iPhi,iTheta)
        
        VelocityExBr(iPhi,iTheta)     = &
             (-Ephi(iPhi,iTheta)*Btheta_G(iPhi,iTheta)) &
             / BmagnitudeSquared_G(iPhi,iTheta)


     enddo
  enddo
  VelocityExBtheta(:,1) = VelocityExBtheta(:,2)
  VelocityExBphi  (:,1) = VelocityExBphi  (:,2)
  
  if (.not.IsMoveFluxTube) Then 
     VelocityExBtheta(:,:) = 0.0
     VelocityExBPhi  (:,:) = 0.0
     VelocityExBr    (:,:) = 0.0
  endif

  if (.not.IsUseJr) Then 
     Jr_G(:,:) = 0.0
  endif


!******************************************************************************
!  Move flux tube around
!******************************************************************************
  
  !initialize field line locations
!  FieldLineTheta (1) = 10.0 * 3.141592653589/180.0
  do iLine=1,nLine
     if (.not. IsRestart) then
        FieldLineTheta (iLine) = 10.0 * 3.141592653589/180.0
        FieldLinePhi   (iLine) = 0.0  * 3.141592653589/180.0
     else
        FieldLineTheta (iLine) = 1.57079632679-GeoMagLat(iLine)&
             * 3.141592653589/180.0
        FieldLinePhi   (iLine) = GeoMagLon(iLine)              &
             * 3.141592653589/180.0
     endif
     
     iFieldLineTheta(iLine) = &
          ceiling(  FieldLineTheta(iLine)  / Dtheta)
     
     iFieldLinePhi(iLine)   = &
          mod(ceiling(FieldLinePhi(iLine)  / Dphi),nPhi-1)
     
     if (FieldLinePhi(iLine) .eq. 0.0)  iFieldLinePhi(iLine) = 1
     
     
     FieldLineX(iLine)      = &
          rLowerBoundary*sin(FieldLineTheta(iLine))*cos(FieldLinePhi(iLine))
     
     FieldLineY(iLine)      = &
          rLowerBoundary*sin(FieldLineTheta(iLine))*sin(FieldLinePhi(iLine))
     
     FieldLineZ(iLine)      = &
          rLowerBoundary*cos(FieldLineTheta(iLine))
     
     OldFieldLineX(iLine)   = FieldLineX(iLine)
     OldFieldLineY(iLine)   = FieldLineY(iLine)
     OldFieldLineZ(iLine)   = FieldLineZ(iLine)
  enddo

!Call New_Interpolate_Velocity( &
!              Theta_G(iFieldLinePhi(1),iFieldLineTheta(1)),           &
!              Theta_G(iFieldLinePhi(1),iFieldLineTheta(1)-1),         &
!              FieldLineTheta(1),                                          &
!              VelocityExBTheta(iFieldLinePhi(1),iFieldLineTheta(1)),  &
!              VelocityExBTheta(iFieldLinePhi(1),iFieldLineTheta(1)-1),&
!              VelocityExBPhi(iFieldLinePhi(1),iFieldLineTheta(1)),    &
!              VelocityExBPhi(iFieldLinePhi(1),iFieldLineTheta(1)-1),  &
!              FieldLineVelTheta(1),FieldLineVelPhi(1))

  AccelerationTheta(:)=0.0
  AccelerationPhi  (:)=0.0
  


  !write(13,*) 'VARIABLES = "X", "Y", "Z", "Ux", "Uy", "Uz"'
  !write(13,*) 'Zone I=', int(maxTime/DtOutput),' ,',' DATAPACKING=POINT'  
  !
  !write(14,*) 'VARIABLES = "X", "Y", "Z", "Time", "Jr"'
  !write(14,*) 'Zone I=', int(maxTime/DtOutput),' ,',' DATAPACKING=POINT'  
  
  do
     if (Time .ge. Tmax) exit
     !if (GeoMagLon(1) .ge. 170.0 .and. GeoMagLon(1) .le. 172) exit
     do iLine=1,nLine
        





        ! if the field line is in the polar circle then set velocity equal
        ! to that of the nearest point. Otherwise, get the velocity from an
        !interpolation.
 
        if (iFieldLineTheta(iLine) .le. 1) then
           write(*,*) 'TT',FieldLineTheta(iLine)/Dtheta,FieldLineZ(iLine),OldFieldLineZ(iLine)
           


           !stop
           FieldLineVelTheta(iLine)= &
                VelocityExBTheta(iFieldLinePhi(iLine),iFieldLineTheta(iLine))
           FieldLineVelPhi(iLine)= &
                VelocityExBPhi(iFieldLinePhi(iLine),iFieldLineTheta(iLine))
           
           FieldLineJr(iLine)    = &
                Jr_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine))
           write(*,*) iFieldLineTheta(iline),time,FieldLineVelX(iLine),FieldLineX(iLine)
        else
           Call New_Interpolate_Velocity( &
              Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),           &
              Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),         &
              FieldLineTheta(iLine),                                          &
              VelocityExBTheta(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),  &
              VelocityExBTheta(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),&
              VelocityExBPhi(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),    &
              VelocityExBPhi(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),  &
              FieldLineVelTheta(iLine),FieldLineVelPhi(iLine))

           Call New_Interpolate_Jr( &
              Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),           &
              Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),         &
              FieldLineTheta(iLine),                                          &
              Jr_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),  &
              Jr_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)-1),&
              FieldLineJr(iLine))
        endif


        
        OldFieldLineTheta(iLine) = FieldLineTheta(iLine)
        OldFieldLinePhi(iLine)   = FieldLinePhi(iLine)
        
        
        OldFieldLineX(iLine) = FieldLineX(iLine)
        OldFieldLineY(iLine) = FieldLineY(iLine)
        OldFieldLineZ(iLine) = FieldLineZ(iLine)

        
        FieldLineVelx(iLine)   = &
             FieldLineVelTheta(iLine)*cos(OldFieldLineTheta(iLine)) &
             * cos(OldFieldLinePhi(iLine)) &
             - FieldLineVelPhi(iLine)*sin(OldFieldLinePhi(iLine)) 
        
        FieldLineVely(iLine)   = &
             FieldLineVelTheta(iLine)*cos(OldFieldLineTheta(iLine)) &
             * sin(OldFieldLinePhi(iLine)) & 
            + FieldLineVelPhi(iLine)*cos(OldFieldLinePhi(iLine))
        
        FieldLineVelz(iLine)   = &
             -FieldLineVelTheta(iLine)*sin(OldFieldLineTheta(iLine))
        
        FieldLineX(iLine) =  &
             (FieldLineVelx(iLine)) *Dt     &
             + OldFieldLineX(iLine)                                   
        
        FieldLineY(iLine) =  &
             (FieldLineVely(iLine)) *Dt     &
             + OldFieldLineY(iLine)                                   

        FieldLineZ(iLine) =  &
             (FieldLineVelz(iLine)) *Dt     &
             + OldFieldLineZ(iLine) 
        
        if (iFieldLineTheta(iLine) .le. 1) write(*,*) FieldLineZ(iLine)/rLowerBoundary,Time

        ! Get new angles from new XYZ positions.Use Theta=arccos(Z/R) and
        ! use Phi = arctan(y/x). make sure to use case statements for phi
        ! depending on what quadrant the position lies in. 

        if (FieldLineZ(iLine)/rLowerBoundary .ge.1) then
           FieldLineTheta(iLine) = 1.57079632679 &
                - acos(sqrt(FieldLineX(iLine)**2+FieldLineY(iLine)**2) & 
                / rLowerBoundary)
        else
           FieldLineTheta(iLine) = &
                abs(acos(FieldLineZ(iLine)/rLowerBoundary))
        endif




        if (FieldLineX(iLine) .gt. 0.0 .and. FieldLineY(iLine) .ge. 0.0) then
           
           FieldLinePhi(iLine)   = &
                (atan(FieldLineY(iLine)/FieldLineX(iLine)))
        else if(FieldLineX(iLine) .lt. 0.0 .and. FieldLineY(iLine).ge.0.0) then
           FieldLinePhi(iLine)   = &
                3.14159-(atan(abs(FieldLineY(iLine)/FieldLineX(iLine))))
        else if(FieldLineX(iLine) .lt.0.0.and. FieldLineY(iLine) .le. 0.0) then
            FieldLinePhi(iLine)   = &
                 3.14159+(atan(abs(FieldLineY(iLine)/FieldLineX(iLine))))
        else if(FieldLineX(iLine) .gt.0.0.and. FieldLineY(iLine) .le. 0.0) then
           FieldLinePhi(iLine)   = &
                 6.283185-(atan(abs(FieldLineY(iLine)/FieldLineX(iLine))))
        else if(FieldLineX(iLine) .eq.0.0.and. FieldLineY(iLine) .ge.0.0) then
           FieldLinePhi(iLine)   = &
                3.14159/2.0
        else if(FieldLineX(iLine) .eq.0.0.and. FieldLineY(iLine) .le.0.0) then
           FieldLinePhi(iLine)   = &
                3.0*3.14159/2.0
        endif


        ! Deal with posibility that phi is negative
        if (FieldLinePhi(iLine) .lt. 0.0) then
           write(*,*) 'TTT',FieldLinePhi(iLine)
           FieldLinePhi(iLine) = FieldLinePhi(iLine) + 6.283185
           write(*,*) 'TTTT',FieldLinePhi(iLine)
           write(*,*) FieldLineX(iLine),FieldLineY(iLine),FieldLineZ(iLine)
           stop
        endif

             

        ! Extract the nearest grid point greater then or equal to the 
        ! actual location. The result is used to get the velocity. 
        iFieldLineTheta(iLine) = &
             ceiling(  FieldLineTheta(iLine)  / Dtheta)
        
        iFieldLinePhi(iLine)   = &
             mod(ceiling(FieldLinePhi(iLine)  / Dphi),nPhi-1)
        if (iFieldLinePhi(iLine) .eq. 0)  iFieldLinePhi(iLine) = 1
        if (iFieldLineTheta(iLine) .eq. 0)  iFieldLineTheta(iLine) = 1
        
        ! Put field line locations back on the sphere

        if (FieldLineTheta(iLine) .ne.0.0 ) then
           FieldLineX(iLine)      = &
                rLowerBoundary*sin(FieldLineTheta(iLine))     &
                * cos(FieldLinePhi(iLine))
           
           FieldLineY(iLine)      = &
                rLowerBoundary*sin(FieldLineTheta(iLine))     &
                * sin(FieldLinePhi(iLine))
    
           FieldLineZ(iLine)      = &
                rLowerBoundary*cos(FieldLineTheta(iLine))
           
        endif

!        if (floor(Time/DToutput).NE.floor((Time-DT)/DToutput) ) Then
!           write(13,*) &
!              OldFieldLineX(iLine),OldFieldLineY(iLine),OldFieldLineZ(iLine), &
!              FieldLineVelx(iLine),FieldLineVely(iLine),FieldLineVelz(iLine)
!           
!           write(14,*) FieldLineX(iLine),FieldLineY(iLine),FieldLineZ(iLine), &
!              Time,FieldLineJr(iLine)
!        endif


   !***************************************************************************
   !  Call the flux tube to be solved
   !***************************************************************************
        ! Get the GeoMagnetic latitude and longitude 
        GeoMagLat(iLine) = (1.57079632679 - OldFieldLineTheta(iLine)) &
             * 180.0 / 3.141592653589
        
        GeoMagLon(iLine) = (OldFieldLinePhi(iLine)) &
             * 180.0 / 3.141592653589
        
!        write(*,*) 'lat',GeoMagLat(iLine),'lon',GeoMagLon(iLine),&
!             'theta',FieldLineTheta(iLine)* 180.0 / 3.141592653589,&
!             'phi',FieldLinePhi(iLine)* 180.0 / 3.141592653589
        
        
        !GeoMagLat(iLine) = 80.0
        !GeoMagLon(iLine) = 0.0
        
        !Calculate the horizontal fieldline velocity in cm/s to be used for 
        !centrifugal force in model
        if (IsCentrifugal) then
           OmegaHorFieldLine(iLine)= sqrt(( &
                FieldLineVelTheta(iLine)**2.0+FieldLineVelPhi(iLine)**2.0) &
                /rLowerBoundary**2.0)
        else
           OmegaHorFieldLine(iLine) = 0.0
        endif


        Call Put_Mod_PW(&
          dOxyg(:,iLine), uOxyg(:,iLine), pOxyg(:,iLine), TOxyg(:,iLine),     &
          dHel(:,iLine), uHel(:,iLine), pHel(:,iLine), THel(:,iLine),         &
          dHyd(:,iLine), uHyd(:,iLine), pHyd(:,iLine), THyd(:,iLine),         &
          dElect(:,iLine), uElect(:,iLine), pElect(:,iLine), TElect(:,iLine), &
          IsRestart,IsVariableDt, Time,DT,DToutput,DTpolarwind,TypeSolver,    &
          GeoMagLat(iLine),GeoMagLon(iLine),FieldLineJr(iLine),               &
          OmegaHorFieldLine(iLine),iUnitInput,iUnitOutput,                    &
          iUnitGraphics(iLine),                                               &
          iUnitSourceGraphics,iUnitRestart(iLine),iUnitCollision,             &
          iUnitRestartIn(iLine),iLine,int(nLine))   
        
        write(*,*) nLine
        

        Call polar_wind

        if (iLine == nLine) then
           Call Get_Mod_PW( &
           dOxyg(:,iLine), uOxyg(:,iLine), pOxyg(:,iLine), TOxyg(:,iLine),    &
           dHel(:,iLine), uHel(:,iLine), pHel(:,iLine), THel(:,iLine),        &
           dHyd(:,iLine), uHyd(:,iLine), pHyd(:,iLine), THyd(:,iLine),        &
           dElect(:,iLine), uElect(:,iLine), pElect(:,iLine), TElect(:,iLine),&
           IsRestart,IsVariableDt,Time,XXX,XXX,XXX,TypeSolver,XXX,XXX,XXX,XXX,&
           XXX,XXX,XXX,XXX,XXX,XXX,XXX,XXX,XXX )
        else
           Call Get_Mod_PW( &
           dOxyg(:,iLine), uOxyg(:,iLine), pOxyg(:,iLine), TOxyg(:,iLine),    &
           dHel(:,iLine), uHel(:,iLine), pHel(:,iLine), THel(:,iLine),        &
           dHyd(:,iLine), uHyd(:,iLine), pHyd(:,iLine), THyd(:,iLine),        &
           dElect(:,iLine), uElect(:,iLine), pElect(:,iLine), TElect(:,iLine),&
           IsRestart,IsVariableDt,XXX,XXX,XXX,XXX,TypeSolver,XXX,XXX,XXX,XXX,&
           XXX,XXX,XXX,XXX,XXX,XXX,XXX,XXX,XXX )
        endif
        





        
        
     enddo
!     Time = Time + Dt
  enddo

  

!******************************************************************************
!  Write output, use cartesian coords for output
!******************************************************************************
 
! do iPhi=1,nPhi
!     do iTheta=1,nTheta
!        ux(iPhi,iTheta) =  & 
!             VelocityExBtheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
!             * cos(Phi_G(iPhi,iTheta)) &
!             - VelocityExBphi(iPhi,iTheta)*sin(Phi_G(iPhi,iTheta))   &
!             + VelocityExBr(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
!             * cos(Phi_G(iPhi,iTheta))
!        
!        uy(iPhi,iTheta) =  &
!             VelocityExBtheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
!             * sin(Phi_G(iPhi,iTheta)) &
!             + VelocityExBphi(iPhi,iTheta)*cos(Phi_G(iPhi,iTheta))   &
!             + VelocityExBr(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
!             * sin(Phi_G(iPhi,iTheta))
!        
!        uz(iPhi,iTheta) =  &
!             -VelocityExBtheta(iPhi,iTheta)*sin(Theta_G(iPhi,iTheta)) &
!             + VelocityExBr(iPhi,iTheta)  *cos(Theta_G(iPhi,iTheta)) 
!        
!        x(iPhi,iTheta)  =  &
!             rLowerBoundary*sin(Theta_G(iPhi,iTheta))*cos(Phi_G(iPhi,iTheta))
!        
!        y(iPhi,iTheta)  =  &
!             rLowerBoundary*sin(Theta_G(iPhi,iTheta))*sin(Phi_G(iPhi,iTheta))
!        
!        z(iPhi,iTheta)  =  &
!             rLowerBoundary*cos(Theta_G(iPhi,iTheta))
!        
!        Ex(iPhi,iTheta) =  & 
!             Etheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
!             * cos(Phi_G(iPhi,iTheta))                     &
!             - Ephi(iPhi,iTheta)*sin(Phi_G(iPhi,iTheta))   &
!             + Er(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
!             * cos(Phi_G(iPhi,iTheta))
!        
!        Ey(iPhi,iTheta) =  &
!             Etheta(iPhi,iTheta)*cos(Theta_G(iPhi,iTheta)) &
!             * sin(Phi_G(iPhi,iTheta))                     &
!             + Ephi(iPhi,iTheta)*cos(Phi_G(iPhi,iTheta))   &
!             + Er(iPhi,iTheta)  *sin(Theta_G(iPhi,iTheta)) &
!             * sin(Phi_G(iPhi,iTheta))
!        
!        Ez(iPhi,iTheta) =  &
!             -Etheta(iPhi,iTheta)*sin(Theta_G(iPhi,iTheta))&
!             + Er(iPhi,iTheta)  *cos(Theta_G(iPhi,iTheta)) 
!     enddo
!  enddo
!  
!  
!  write(12,*) &
!     'VARIABLES = "X", "Y", "Z", "Ux", "Uy", "Uz", "V", "Ex", "Ey", "Ez", "Jr"'
!
!  write(12,*) 'Zone I=', nPhi, ', J=', nTheta,', DATAPACKING=POINT'
!
!  do iTheta=1,nTheta
!     do iPhi=1,nPhi
!        write(12,*) &
!             x(iPhi,iTheta),y(iPhi,iTheta),z(iPhi,iTheta),     &
!             ux(iPhi,iTheta),uy(iPhi,iTheta),uz(iPhi,iTheta),  &
!             Potential_G(iPhi,iTheta),                         &
!             Ex(iPhi,iTheta), Ey(iPhi,iTheta), Ez(iPhi,iTheta),&
!             Jr_G(iPhi,iTheta)
!        
!     enddo
!  enddo
  
  do iLine=1,nLine
     CLOSE(UNIT=iUnitRestart(iLine))
     CLOSE(UNIT=iUnitGraphics(iLine))

  enddo
  close(iUnitNorth)

  CLOSE(UNIT=iUnitOutput)
  CLOSE(UNIT=iUnitInput)
  !CLOSE(UNIT=iUnitCollision)
  !CLOSE(UNIT=iUnitSourceGraphics)
  
  call MPI_FINALIZE(errcode)
  
end program pw



!******************************************************************************
!  Determine the interpolated velocity at point theta3, and phi3. 
!  Uses Shepard Interpolation given on page 430 of 
!  Numerical Analysis by Kincaid
!******************************************************************************


Subroutine interpolate_velocity(  Theta1,Theta2,Theta3, &
                                  Phi1,  Phi2,  Phi3,   &
                                  VelocityTheta11, VelocityTheta12, &
                                  VelocityTheta21, VelocityTheta22, &
                                  VelocityPhi11,   VelocityPhi12,   &
                                  VelocityPhi21,   VelocityPhi22,   &
                                  VelocityThetaOut,VelocityPhiOut)


  real, intent(in)  ::           Theta1,Theta2,Theta3, &
                                Phi1,  Phi2,  Phi3,   &
                                VelocityTheta11, VelocityTheta12, &
                                VelocityTheta21, VelocityTheta22, &
                                VelocityPhi11,   VelocityPhi12,   &
                                VelocityPhi21,   VelocityPhi22   
  
  real, intent(out)  ::          VelocityThetaOut,VelocityPhiOut
  

  real,dimension(4) ::          Theta, Phi,VelocityTheta,VelocityPhi,numerator
  
  real              ::          multiplyer
  real, dimension(4,4) ::       d

  integer           ::          i,j

  Theta(1)         = Theta1
  Phi(1)           = Phi1
  VelocityTheta(1) = VelocityTheta11
  VelocityPhi(1)   = VelocityPhi11

  Theta(2)         = Theta1
  Phi(2)           = Phi2
  VelocityTheta(2) = VelocityTheta21
  VelocityPhi(2)   = VelocityPhi21

  Theta(3)         = Theta2
  Phi(3)           = Phi1
  VelocityTheta(3) = VelocityTheta12
  VelocityPhi(3)   = VelocityPhi12

  Theta(4)         = Theta2
  Phi(4)           = Phi2
  VelocityTheta(4) = VelocityTheta22
  VelocityPhi(4)   = VelocityPhi22


  do i=1,4
     
     Numerator(i) = (Theta3-Theta(i))**2.0 + (Phi3-Phi(i))**2.0

     do j=1,4
        d(i,j)= (Theta(i)-Theta(j))**2.0 + (Phi(i)-Phi(j))**2.0
     enddo
  enddo
     
  VelocityThetaOut = 0.0
  VelocityPhiOut   = 0.0
  
  do i=1,4
     multiplyer = 1.0
     do j=1,4
        if (j .ne. i) then
           multiplyer=multiplyer*numerator(j)/(d(i,j))
        endif
     enddo
     VelocityThetaOut = VelocityThetaOut+VelocityTheta(i)*multiplyer
     VelocityPhiOut   = VelocityPhiOut  +VelocityPhi(i)  *multiplyer
     
  enddo
     
end Subroutine interpolate_velocity


Subroutine New_Interpolate_Velocity(Theta1,Theta2,Theta3,uTheta1,&
                                    uTheta2,uPhi1,uPhi2,uTheta3,uPhi3)

real, intent(in)            ::   Theta1,Theta2,Theta3,uTheta1,&
                                 uTheta2,uPhi1,uPhi2

real, intent(out)           ::   uTheta3,uPhi3
  
  
uTheta3 = (uTheta2 - uTheta1) / (Theta2-Theta1) * (Theta3-Theta1) + uTheta1

uPhi3   = (uPhi2   - uPhi1)   / (Theta2-Theta1) * (Theta3-Theta1) + uPhi1
  
  
end Subroutine New_Interpolate_Velocity


Subroutine New_Interpolate_Jr(Theta1,Theta2,Theta3,Jr1,&
                                    Jr2,Jr3)

real, intent(in)            ::   Theta1,Theta2,Theta3,Jr1,&
                                 Jr2

real, intent(out)           ::   Jr3
  
  
Jr3 = (Jr2 - Jr1) / (Theta2-Theta1) * (Theta3-Theta1) + Jr1


  
  
end Subroutine New_Interpolate_Jr



!******************************************************************************
!  Put polarwind variables into Mod_PW for passing 
!*****************************************************************************


Subroutine Put_Mod_PW(dOxyg, uOxyg, pOxyg, TOxyg,     &
                      dHel, uHel, pHel, THel,         &
                      dHyd, uHyd, pHyd, THyd,         &
                      dElect, uElect, pElect, TElect, &
                      IsRestart,IsVariableDt,Time,DT,DToutput,DTpolarwind,&
                      TypeSolver,GeoMagLat,GeoMagLon,Jr,wHorizontal,      &
                      iUnitInput,iUnitOutput,iUnitGraphics,               &
                      iUnitSourceGraphics,iUnitRestart,iUnitCollision,    &
                      iUnitRestartIn,iLine,nLine     )

  use ModParameters
  use ModPass
  
  real, intent(in),dimension(maxGrid):: dOxyg, uOxyg, pOxyg, TOxyg,     &
                                        dHel, uHel, pHel, THel,         &
                                        dHyd, uHyd, pHyd, THyd,         &
                                        dElect, uElect, pElect, TElect
 
  real,    intent(in)     :: Time, DT, DToutput, DTpolarwind,GeoMagLat, &
                             GeoMagLon,Jr, wHorizontal                  
  integer, intent(in)     :: iUnitInput,iUnitOutput,iUnitGraphics,      &
                             iUnitSourceGraphics,iUnitRestart,          &
                             iUnitCollision,iUnitRestartIn,iLine,nLine     
  logical, intent(in)     :: IsRestart,IsVariableDt
  CHARACTER(7), intent(in):: TypeSolver

  dOxygPW (:) = dOxyg(:)
  uOxygPW (:) = uOxyg(:)
  pOxygPW (:) = pOxyg(:)
  TOxygPW (:) = TOxyg(:)
  dHelPW  (:) = dHel(:)
  uHelPW  (:) = uHel(:)
  pHelPW  (:) = pHel(:)
  THelPW  (:) = THel(:)
  dHydPW  (:) = dHyd(:)
  uHydPW  (:) = uHyd(:)
  pHydPW  (:) = pHyd(:)
  THydPW  (:) = THyd(:)
  dElectPW(:) = dElect(:)
  uElectPW(:) = uElect(:)
  pElectPW(:) = pElect(:)
  TElectPW(:) = TElect(:) 
  IsRestartPW = IsRestart
  IsVariableDtPW=IsVariableDT
  TimePW      = Time
  TmaxPW      = Time+DT
  DToutputPW  = DToutput
  DTpolarwindPW = DTpolarwind
  TypeSolverPW  = TypeSolver
  GeoMagLatPW = GeoMagLat
  GeoMagLonPW = GeoMagLon
  JrPW        = Jr
  wHorizontalPW   = wHorizontal
  iUnitInputPW    =iUnitInput
  iUnitOutputPW   =iUnitOutput
  iUnitGraphicsPW =iUnitGraphics
  iUnitSourceGraphicsPW=iUnitSourceGraphics
  iUnitRestartPW  =iUnitRestart
  iUnitCollisionPW=iUnitCollision
  iUnitRestartInPW=iUnitRestartIn
  iLinePW = iLine
  nLinePW = nLine
end Subroutine Put_Mod_PW




!******************************************************************************
!  Get polarwind variables from Mod_PW 
!*****************************************************************************


Subroutine Get_Mod_PW(dOxyg, uOxyg, pOxyg, TOxyg,     &
                      dHel, uHel, pHel, THel,         &
                      dHyd, uHyd, pHyd, THyd,         &
                      dElect, uElect, pElect, TElect, &
                      IsRestart,IsVariableDt,Time,Tmax,DToutput,DTpolarwind,&
                      TypeSolver,GeoMagLat,GeoMagLon,Jr,wHorizontal,        &
                      iUnitInput,iUnitOutput,iUnitGraphics,               &
                      iUnitSourceGraphics,iUnitRestart,iUnitCollision,    &
                      iUnitRestartIn,iLine,nLine     )

  use ModParameters
  use ModPass
  
  real, intent(out),dimension(maxGrid):: dOxyg, uOxyg, pOxyg, TOxyg,     &
                                        dHel, uHel, pHel, THel,         &
                                        dHyd, uHyd, pHyd, THyd,         &
                                        dElect, uElect, pElect, TElect
 
  real,    intent(out)     :: Time,DToutput, DTpolarwind,Tmax,&
                              GeoMagLat,GeoMagLon, Jr, wHorizontal

  integer,intent(out)      :: iUnitInput,iUnitOutput,iUnitGraphics,      &
                              iUnitSourceGraphics,iUnitRestart,          &
                              iUnitCollision,iUnitRestartIn,iLine,nLine
  logical, intent(out)     :: IsRestart,IsVariableDt
  CHARACTER(7),intent(out) :: TypeSolver
  

  dOxyg (:) = dOxygPW(:)
  uOxyg (:) = uOxygPW(:)
  pOxyg (:) = pOxygPW(:)
  TOxyg (:) = TOxygPW(:)
  dHel  (:) = dHelPW(:)
  uHel  (:) = uHelPW(:)
  pHel  (:) = pHelPW(:)
  THel  (:) = THelPW(:)
  dHyd  (:) = dHydPW(:)
  uHyd  (:) = uHydPW(:)
  pHyd  (:) = pHydPW(:)
  THyd  (:) = THydPW(:)
  dElect(:) = dElectPW(:)
  uElect(:) = uElectPW(:)
  pElect(:) = pElectPW(:)
  TElect(:) = TElectPW(:) 
  IsRestart = IsRestartPW
  IsVariableDT = IsVariableDtPW
  Time      = TimePW
  Tmax      = TmaxPW
  DToutput  = DToutputPW
  DTpolarwind = DTpolarwindPW
  TypeSolver  = TypeSolverPW
  GeoMagLat = GeoMagLatPW
  GeoMagLon = GeoMagLonPW
  Jr        = JrPW
  wHorizontal = wHorizontalPW
  iUnitInput    =iUnitInputPW
  iUnitOutput   =iUnitOutputPW
  iUnitGraphics =iUnitGraphicsPW
  iUnitSourceGraphics=iUnitSourceGraphicsPW
  iUnitRestart  =iUnitRestartPW
  iUnitCollision=iUnitCollisionPW
  iUnitRestartIn=iUnitRestartInPW
  iLine=iLinePW
  nLine=nLinePW
end Subroutine Get_Mod_PW







!        FieldLineVelTheta(iLine) = & 
!             VelocityExBTheta(iFieldLinePhi(iLine),iFieldLineTheta(iLine))

!        FieldLineVelPhi(iLine) = & 
!             VelocityExBPhi(iFieldLinePhi(iLine),iFieldLineTheta(iLine))


!        Call interpolate_velocity(  &
!            Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)), &
!            Theta_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)+1),&
!            FieldLineTheta(iLine),           &
!            Phi_G(iFieldLinePhi(iLine),iFieldLineTheta(iLine)), &
!            Phi_G(iFieldLinePhi(iLine)+1,iFieldLineTheta(iLine)),&
!            FieldLinePhi(iLine), &
!            VelocityExBTheta(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),&
!            VelocityExBTheta(iFieldLinePhi(iLine),iFieldLineTheta(iLine)+1),&
!            VelocityExBTheta(iFieldLinePhi(iLine)+1,iFieldLineTheta(iLine)),&
!            VelocityExBTheta(iFieldLinePhi(iLine)+1,iFieldLineTheta(iLine)+1),&
!            VelocityExBPhi(iFieldLinePhi(iLine),iFieldLineTheta(iLine)),&
!            VelocityExBPhi(iFieldLinePhi(iLine),iFieldLineTheta(iLine)+1),&
!            VelocityExBPhi(iFieldLinePhi(iLine)+1,iFieldLineTheta(iLine)),&
!            VelocityExBPhi(iFieldLinePhi(iLine)+1,iFieldLineTheta(iLine)+1),&
!            FieldLineVelTheta(iLine),FieldLineVelPhi(iLine))
