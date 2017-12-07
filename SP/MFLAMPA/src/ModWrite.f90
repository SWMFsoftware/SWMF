module SP_ModWrite

  ! This module contains methods for writing output files

  use SP_ModSize, ONLY: &
       nDim, nLat, nLon, nNode, &
       iParticleMin, iParticleMax, nParticle,&
       nMomentumBin, &
       Particle_, OriginLat_, OriginLon_

  use SP_ModGrid, ONLY: &
       get_node_indexes, &
       iComm, nProc, &
       nVar, nVarRead, nBlock, State_VIB, iGridLocal_IB, iNode_B, &
       Distribution_IIB, LogEnergyScale_I, LogMomentumScale_I, &
       DMomentumOverDEnergy_I, &
       Proc_, Begin_, End_, Shock_, X_, Y_, Z_, Bx_, By_, Bz_, Wave1_,Wave2_,&
       B_, Ux_, Uy_, Uz_, U_, Rho_, T_, S_, LagrID_, DLogRho_,  &
       EFlux_, Flux0_, Flux1_, Flux2_, Flux3_, Flux4_, Flux5_, Flux6_, &
       NameVar_V

  use SP_ModAdvance, ONLY: TimeGlobal, iIterGlobal, DoTraceShock

  use ModPlotFile, ONLY: save_plot_file, read_plot_file

  use ModUtilities, ONLY: split_string

  implicit none

  SAVE

  private ! except
 
  public:: set_write_param, write_output, NamePlotDir

  type TypeOutputFile
     !\
     ! General information
     !/
     ! kind of data printed to a file
     integer:: iKindData
     ! format of a file
     integer:: iFormat
     character(len=4 ):: NameFormat
     character(len=20):: TypeFile
     ! whether it is the first call
     ! USED ONLY IN write_mh_time  FOR NOW!!!
     logical:: IsFirstCall
     ! names of variables to be written
     character(len=300):: NameVarPlot
     ! output buffer
     real, pointer:: Buffer_II(:,:)
     !\
     ! MH data
     !/
     ! variables from the state vector to be written
     logical:: DoPlot_V(nVar)
     ! total number of variables to be written
     integer:: nVarPlot
     ! their indices in the state vector
     integer, pointer:: iVarPlot_V(:)
     !\
     ! Distribution
     !/
     ! scale of distribution function (momentum or energy)
     integer:: iScale
     !\
     ! Data on the sphere
     !/
     ! radius of the sphere the data to be written at
     real:: Radius
  end type TypeOutputFile

  !\
  !----------------------------------------------------------------------------
  ! Number of output files
  integer:: nFile = 0
  ! The output files
  type(TypeOutputFile), allocatable:: File_I(:)
  ! Types of output files in terms of output dataa
  integer, parameter:: &
       ! Background mhd data
       MH1D_      = 0, & ! along each line
       MH2D_      = 1, & ! at given radius as Lon-Lat plot
       MHTime_    = 2, & ! at given radius for each line as time series 
       ! Distribution
       Distr1D_   = 3, & ! along each line
       Distr2D_   = 4, & !  at given radius as Lon-Lat plot
       DistrTime_ = 5    ! at given radius for each line as time series 
  ! Format of output files
  integer, parameter:: &
       Tec_ = 0, &
       Idl_ = 1
  ! Scale for writing distribution
  integer, parameter:: &
       MomentumScale_= 1, &
       EnergyScale_  = 2
  !----------------------------------------------------------------------------
  ! the output directory
  character (len=100) :: NamePlotDir="SP/IO2/"
  !----------------------------------------------------------------------------
  ! auxilary array, used to write data on a sphere
  integer:: iNodeIndex_I(nNode)
  !/
contains
  
  subroutine set_write_param
    use ModReadParam, ONLY: read_var
    ! set parameters of output files: file format, kind of output etc.
    character(len=300):: StringPlot
    ! loop variables
    integer:: iFile, iNode
    character(len=*), parameter :: NameSub='SP:set_write_param'
    !--------------------------------------------------------------------------
    ! initialize auxilary array
    do iNode = 1, nNode
       iNodeIndex_I(iNode) = iNode
    end do
    ! number of output files
    call read_var('nFile', nFile)
    ! check correctness
    if(nFile == 0) RETURN ! no output file requested
    if(nFile  < 0) call CON_stop(NameSub//': incorrect SAVEPLOT section')

    ! allocate the storage for file info
    allocate(File_I(nFile))

    ! read info about each file
    do iFile = 1, nFile
       ! reset and read the file info
       StringPlot = ''
       call read_var('StringPlot', StringPlot)

       ! the format of output file must be set
       if(    index(StringPlot,'tec') > 0)then
          File_I(iFile) % iFormat   = Tec_
          File_I(iFile) % NameFormat='.dat'
          File_I(iFile) % TypeFile  ='tec'
       elseif(index(StringPlot,'idl') > 0)then
          File_I(iFile) % iFormat   = Idl_
          File_I(iFile) % NameFormat='.out'
          File_I(iFile) % TypeFile  ='ascii'
       else
          call CON_stop(NameSub//': output format was not set in PARAM.in')
       end if

       ! the kind of output data must be set
       if(    index(StringPlot,'mh1d' )    > 0)then
          File_I(iFile) % iKindData = MH1D_
       elseif(index(StringPlot,'mh2d')    > 0)then
          File_I(iFile) % iKindData = MH2d_
       elseif(index(StringPlot,'mhtime')    > 0)then
          File_I(iFile) % iKindData = MHTime_
       elseif(index(StringPlot,'distr1d')  > 0)then
          File_I(iFile) % iKindData = Distr1D_
       elseif(index(StringPlot,'distr2d') > 0)then
          File_I(iFile) % iKindData = Distr2d_
          call CON_stop(NameSub//&
               ': spherical output for distribution is not implemented yet')
       elseif(index(StringPlot,'distrtime')    > 0)then
          File_I(iFile) % iKindData = DistrTime_
          call CON_stop(NameSub//&
               ': time series output for distribution is not implemented yet')
       else
          call CON_stop(NameSub//': kind of data was not set in PARAM.in')
       end if

       ! reset variables' names
       File_I(iFile) % NameVarPlot = ''

       ! based on kind of data process the requested output
       select case(File_I(iFile) % iKindData)
          case(MH1D_)
             call process_mh
             ! prepare the output data container
             allocate(File_I(iFile) % Buffer_II(&
                  File_I(iFile)%nVarPlot, iParticleMin:iParticleMax))
             ! add particle index to variable names
             File_I(iFile) % NameVarPlot = &
                  'ParticleID '//trim(File_I(iFile) % NameVarPlot)
             if(DoTraceShock)&
                  File_I(iFile) % NameVarPlot = &
                  trim(File_I(iFile) % NameVarPlot)//' iShock RShock'
          case(MH2D_)
             call process_mh
             ! prepare the output data container
             allocate(File_I(iFile) % Buffer_II(&
                  File_I(iFile)%nVarPlot, nNode))
             ! add line index to variable names
             File_I(iFile) % NameVarPlot = &
                  'LineIndex '//trim(File_I(iFile) % NameVarPlot)
             ! get radius
             call read_var('Radius [Rs]', File_I(iFile) % Radius)
          case(MHTime_)
             call process_mh
             ! prepare the output data container
             allocate(File_I(iFile) % Buffer_II(&
                  File_I(iFile)%nVarPlot, 1))
             ! add time interval index to variable names
             File_I(iFile) % NameVarPlot = &
                  'TimeIntervalIndex '//trim(File_I(iFile) % NameVarPlot)
             ! get radius
             call read_var('Radius [Rs]', File_I(iFile) % Radius)
             ! reset indicator of the first call
             File_I(iFile) % IsFirstCall = .true.
          case(Distr1D_)
             call process_distr
             ! prepare the output data container
             allocate(File_I(iFile) % &
                  Buffer_II(nMomentumBin,iParticleMin:iParticleMax))
          case(Distr2D_)
             call process_distr
             ! prepare the output data container
             allocate(File_I(iFile) % Buffer_II(nMomentumBin, nNode))
             ! get radius
             call read_var('Radius [Rs]', File_I(iFile) % Radius)
       end select       
    end do

  contains
    subroutine process_mh
      ! process variables to plot
      ! NOTE: for iKindData == MH1D_ certain variables are always printed:
      !       Rho_, T_, Ux_:Uz_, Bx_:Bz_, Wave1_, Wave2_
      integer:: iVar, iVarPlot, nStringPlot, iStringPlot
      character(len=10) :: StringPlot_I(2*nVar)
      !-----------------------------------------------
      ! reset
      File_I(iFile) % DoPlot_V = .false.
      ! for MH1D_ minimal set of variables is printed
      if(File_I(iFile)%iKindData == MH1D_)then
         ! for MH1D_ minimal set of variables is printed
         File_I(iFile) % DoPlot_V(1:nVarRead) = .true.
      else
         ! coordinates are always printed
         File_I(iFile) % DoPlot_V(1:nDim) = .true.
      end if

      ! get individual variables' names
      call split_string(StringPlot, StringPlot_I, nStringPlot)
      do iStringPlot = 1, nStringPlot
         select case(StringPlot_I(iStringPlot))
         case('dist')
            ! distance along line -----
            File_I(iFile) % DoPlot_V(S_) = .true.
         case('lagr')
            ! lagrangian coord along line -----
            File_I(iFile) % DoPlot_V(LagrID_) = .true.
         case('rho')
            ! plasma density ----------
            File_I(iFile) % DoPlot_V(Rho_) = .true.
         case('drho')
            ! lagrangian derivative of plasma density ----------
            File_I(iFile) % DoPlot_V(DLogRho_) = .true.
         case('temp')
            ! temperature -------------
            File_I(iFile) % DoPlot_V(T_) = .true.
         case('ux','uy','uz')
            ! velocity ----------------
            File_I(iFile) % DoPlot_V((/Ux_, Uy_, Uz_/)) = .true.
         case('|u|')
            ! velocity ----------------
            File_I(iFile) % DoPlot_V(U_) = .true.
         case('bx','by','bz')
            ! magnetic field ----------
            File_I(iFile) % DoPlot_V((/Bx_, By_, Bz_/)) = .true.
         case('|b|')
            ! magnetic field ----------
            File_I(iFile) % DoPlot_V(B_) = .true.
         case('wave', 'wave1', 'wave2')
            ! turbulence intensity ----
            File_I(iFile) % DoPlot_V((/Wave1_,Wave2_/)) = .true.
         case('flux')
            ! energy flux -------------
            File_I(iFile) % DoPlot_V(EFlux_) = .true.
            File_I(iFile) % DoPlot_V(Flux0_:Flux6_) = .true.
         end select
      end do
      File_I(iFile) % nVarPlot = count(File_I(iFile) % DoPlot_V)
      ! indices in the state vector
      allocate(File_I(iFile) % iVarPlot_V(File_I(iFile)%nVarPlot))
      ! determine indices and names of variables
      iVarPlot = 1
      do iVar = 1, nVar
         if(.not.File_I(iFile) % DoPlot_V(iVar)) CYCLE
         File_I(iFile) % NameVarPlot = &
              trim(File_I(iFile) % NameVarPlot)//' '//&
              trim(NameVar_V(iVar))
         File_I(iFile) % iVarPlot_V(iVarPlot) = iVar
         iVarPlot = iVarPlot + 1
      end do
    end subroutine process_mh
    !------------------------------------------------------------------------
    subroutine process_distr
      ! process output parameters for distribution output
      !------------------------
      File_I(iFile) % nVarPlot = 1
      if(    index(StringPlot, 'momentum') > 0)then
         File_I(iFile) % iScale = MomentumScale_
         File_I(iFile) % NameVarPlot = &
              'Log10Momentum Distance Log10Distribution'
      elseif(index(StringPlot, 'energy') > 0)then
         File_I(iFile) % iScale = EnergyScale_
         File_I(iFile) % NameVarPlot = &
              'Log10Energy Distance Log10Distribution'
      else
         call CON_stop(NameSub//&
              ': type of scale for distribution output wasnot set in PARAM.in')
      end if
    end subroutine process_distr
  end subroutine set_write_param

  !============================================================================

  subroutine write_output(IsInitialOutput)
    ! write the output data
    logical, intent(in), optional:: IsInitialOutput

    ! loop variables
    integer:: iFile
    integer:: iKindData

    logical:: IsInitialOutputLocal

    character(len=*), parameter:: NameSub = 'SP:write_output'
    !--------------------------------------------------------------------------
    ! check whether this is a call for initial output
    if(present(IsInitialOutput))then
       IsInitialOutputLocal = IsInitialOutput
    else
       IsInitialOutputLocal = .false.
    end if

    if(nFile == 0) RETURN
    
    do iFile = 1, nFile
       iKindData = File_I(iFile) % iKindData

       ! during initial call only background 1D data is printed
       if(IsInitialOutputLocal .and. iKindData /= MH1D_)&
            iKindData = -1
           
       select case(iKindData)
       case(MH1D_)
          call write_mh_1d
       case(MH2D_)
          call write_mh_2d
       case(MHTime_)
          call write_mh_time
       case(Distr1D_)
          call write_distr_1d
       case(Distr2D_)
          !call write_distr_2d
       end select
    end do

  contains

    subroutine write_mh_1d
      ! write output with 1D MH data in the format to be read by IDL/TECPLOT;
      ! separate file is created for each field line, name format is
      ! MH_data_<iLon>_<iLat>_n<ddhhmmss>_n<iIter>.{out/dat}
      !------------------------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      integer:: iFirst, iLast
      ! for better readability
      integer:: nVarPlot
      ! shock location
      integer:: iShock
      real   :: RShock
      ! timestamp
      character(len=8):: StringTime
      !------------------------------------------------------------------------
      nVarPlot = File_I(iFile) % nVarPlot
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! set the file name
         call get_time_string(TimeGlobal, StringTime)
         write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
              trim(NamePlotDir)//'MH_data_',iLon,'_',iLat,&
              '_t'//StringTime//'_n',iIterGlobal,&
              File_I(iFile) % NameFormat

         ! get min and max particle indexes on this field line
         iFirst = iGridLocal_IB(Begin_, iBlock)
         iLast  = iGridLocal_IB(End_,   iBlock)

         do iParticle = iFirst, iLast
            ! fill the output buffer
            do iVarPlot = 1, nVarPlot
               File_I(iFile) % Buffer_II(iVarPlot, iParticle) = &
                    State_VIB(File_I(iFile) % iVarPlot_V(iVarPlot), &
                    iParticle, iBlock)
            end do
         end do

         ! shock location
         iShock = iGridLocal_IB(Shock_,iBlock)
         RShock = sqrt(sum(State_VIB(X_:Z_,iShock,iBlock)**2))
         ! print data to file
         call save_plot_file(&
              NameFile     = NameFile, &
              TypeFileIn   = File_I(iFile) % TypeFile, &
              nDimIn       = 1, &
              TimeIn       = TimeGlobal, &
              nStepIn      = iIterGlobal, &
              CoordMinIn_D = (/real(iFirst)/), &
              CoordMaxIn_D = (/real(iLast)/), &
              NameVarIn    = File_I(iFile) % NameVarPlot, &
              VarIn_VI     = &
              File_I(iFile) % Buffer_II(1:nVarPlot,iFirst:iLast),&
              ! additionally print the shock location both as
              ! index and radial distance
              ParamIn_I    = pack(&
              (/real(iShock), RShock/),DoTraceShock)&
              )
      end do
    end subroutine write_mh_1d

    !=========================================================================

    subroutine write_mh_2d
      use ModMpi
      use CON_world, ONLY: is_proc0, i_proc0
      use CON_comp_param, ONLY: SP_
      ! write output with 2D MH data in the format to be read by IDL/TECPLOT;
      ! single file is created for all field lines, name format is
      ! MH_data_R=<Radius [AU]>_t<ddhhmmss>_n<iIter>.{out/dat}
      !------------------------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot, iVarIndex
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      integer:: iFirst, iLast
      ! index of particle just above the radius
      integer:: iAbove
      ! radii of particles, added for readability
      real:: Radius0, Radius1
      ! interpolation weight
      real:: Weight
      ! for better readability
      integer:: nVarPlot
      ! MPI error
      integer:: iError
      ! skip a field line if it fails to reach radius of output sphere
      logical:: DoPrint_I(nNode)
      ! timestamp
      character(len=8):: StringTime
      !------------------------------------------------------------------------
      nVarPlot = File_I(iFile) % nVarPlot
      ! reset the output buffer
      File_I(iFile) % Buffer_II = 0
      
      ! reset, all field lines are printed unless fail to reach output sphere
      DoPrint_I = .true.

      ! go over all lines on the processor and find the point of intersection
      ! with output sphere if present
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)

         ! set the file name
         call get_time_string(TimeGlobal, StringTime)
         write(NameFile,'(a,i4.4,f0.2,a,i6.6,a)') &
              trim(NamePlotDir)//'MH_data_R=', int(File_I(iFile) % Radius), &
              File_I(iFile) % Radius - int(File_I(iFile) % Radius), &
              '_t'//StringTime//'_n', iIterGlobal, File_I(iFile) % NameFormat

         ! get min and max particle indexes on this field line
         iFirst = iGridLocal_IB(Begin_, iBlock)
         iLast  = iGridLocal_IB(End_,   iBlock)

         ! find the particle just above the given radius
         do iParticle = iFirst , iLast
            Radius0 = sum(State_VIB(X_:Z_, iParticle, iBlock)**2)**0.5
            if( Radius0 > File_I(iFile) % Radius) then
               iAbove = iParticle
               !check if line started above output sphere, i.e. no intersection
               DoPrint_I(iNode) = iAbove == iFirst
               EXIT
            end if
            ! check if reached the end, i.e. there is no intersection
            if(iParticle == iLast) &
                 DoPrint_I(iNode) = .false.
         end do

         ! if no intersection found -> proceed to the next line
         if(.not.DoPrint_I(iNode)) CYCLE

         ! intersection is found -> get data at that location;
         ! interpolate data and fill buffer
         Radius0 = sum(State_VIB(X_:Z_, iAbove-1, iBlock)**2)**0.5
         Radius1 = sum(State_VIB(X_:Z_, iAbove,   iBlock)**2)**0.5
         Weight  = (File_I(iFile)%Radius - Radius0) / (Radius1 - Radius0)
         ! interpolate each requested variable
         do iVarPlot = 1, nVarPlot
            iVarIndex = File_I(iFile) % iVarPlot_V(iVarPlot)
            File_I(iFile) % Buffer_II(iVarPlot, iNode) = &
                 State_VIB(iVarIndex, iAbove-1, iBlock) * (1-Weight) + &
                 State_VIB(iVarIndex, iAbove,   iBlock) *    Weight
         end do

      end do
      
      ! gather interpolated data on the source processor
      if(nProc > 1)then
         if(is_proc0(SP_))then
            call MPI_Reduce(MPI_IN_PLACE, File_I(iFile) % Buffer_II, &
                 nNode * File_I(iFile) % nVarPlot, MPI_REAL, MPI_Sum, &
                 i_proc0(SP_), iComm, iError)
            call MPI_Reduce(MPI_IN_PLACE, DoPrint_I, &
                 nNode, MPI_Logical, MPI_Land, &
                 i_proc0(SP_), iComm, iError)
         else
            call MPI_Reduce(File_I(iFile) % Buffer_II, File_I(iFile) % Buffer_II,&
                 nNode * File_I(iFile) % nVarPlot, MPI_REAL, MPI_Sum, &
                 i_proc0(SP_), iComm, iError)
            call MPI_Reduce(DoPrint_I, DoPrint_I, &
                 nNode, MPI_Logical, MPI_Land, &
                 i_proc0(SP_), iComm, iError)
         end if
      end if

      if(is_proc0(SP_))&
           ! print data to file
           call save_plot_file(&
           NameFile     = NameFile, &
           TypeFileIn   = File_I(iFile) % TypeFile, &
           nDimIn       = 1, &
           TimeIn       = TimeGlobal, &
           nStepIn      = iIterGlobal, &
           Coord1In_I   = real(pack(iNodeIndex_I, MASK=DoPrint_I)), &
           NameVarIn    = File_I(iFile) % NameVarPlot, &
           VarIn_VI     = &
           reshape(&
           pack(File_I(iFile) % Buffer_II(1:nVarPlot,1:nNode),&
           MASK = spread(DoPrint_I, 1, nVarPlot)), &
           (/nVarPlot, count(DoPrint_I)/))&
           )
    end subroutine write_mh_2d

    !=========================================================================

    subroutine write_mh_time
      ! write output w/time series MH data in format to be read by IDL/TECPLOT;
      ! a file is created for each field lines, name format is
      ! MH_data_R=<Radius [AU]>_<iLon>_<iLat>.{out/dat}
      ! the file has no timestamp as it is updated during the run
      !------------------------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot, iVarIndex
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      integer:: iFirst, iLast
      ! index of particle just above the radius
      integer:: iAbove
      ! radii of particles, added for readability
      real:: Radius0, Radius1
      ! interpolation weight
      real:: Weight
      ! for better readability
      integer:: nVarPlot
      ! skip a field line if it fails to reach radius of output sphere
      logical:: DoPrint
      ! size of the already written data
      integer:: nDataLine
      ! current size of the buffer
      integer:: nBufferSize
      !------------------------------------------------------------------------
      nVarPlot = File_I(iFile) % nVarPlot
      ! reset the output buffer
      File_I(iFile) % Buffer_II = 0


      ! go over all lines on the processor and find the point of intersection
      ! with output sphere if present
      do iBlock = 1, nBlock

         ! reset, the field line is printed unless fail to reach output sphere
         DoPrint = .true.

         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! get min and max particle indexes on this field line
         iFirst = iGridLocal_IB(Begin_, iBlock)
         iLast  = iGridLocal_IB(End_,   iBlock)

         ! find the particle just above the given radius
         do iParticle = iFirst , iLast
            Radius0 = sum(State_VIB(X_:Z_, iParticle, iBlock)**2)**0.5
            if( Radius0 > File_I(iFile) % Radius)then
               iAbove = iParticle
               !check if line started above output sphere, i.e. no intersection
               DoPrint = iAbove == iFirst
               EXIT
            end if
            ! check if reached the end, i.e. there is no intersection
            if(iParticle == iLast) &
                 DoPrint = .false.
         end do
         ! if no intersection found -> proceed to the next line
         if(.not.DoPrint) CYCLE


         ! set the file name
         write(NameFile,'(a,i4.4,f0.2,a,i3.3,a,i3.3,a)') &
              trim(NamePlotDir)//'MH_data_R=', int(File_I(iFile) % Radius), &
              File_I(iFile) % Radius - int(File_I(iFile) % Radius), &
              '_', iLon, '_', iLat, File_I(iFile) % NameFormat

         !\
         ! if file already exists -> read its content
         nDataLine = 0
         if(.not.File_I(iFile)%IsFirstCall)then
            
            ! first, determine its size
            call read_plot_file(&
                 NameFile   = NameFile, &
                 TypeFileIn = File_I(iFile) % TypeFile, &
                 n1Out      = nDataLine)
            ! if buffer is too small then reallocate it
            nBufferSize = ubound(File_I(iFile)%Buffer_II, 2)
            if(nBufferSize < nDataLine + 1)then
               deallocate(File_I(iFile) % Buffer_II)
               allocate(File_I(iFile) % Buffer_II(&
                    File_I(iFile)%nVarPlot, 2*nBufferSize))
            end if

            ! read the data itself
            call read_plot_file(&
                 NameFile   = NameFile, &
                 TypeFileIn = File_I(iFile) % TypeFile,&
                 VarOut_VI  = File_I(iFile) % Buffer_II)
         end if

         !\
         ! add new data
         nDataLine = nDataLine+1

         ! interpolate data and fill buffer
         Radius0 = sum(State_VIB(X_:Z_, iAbove-1, iBlock)**2)**0.5
         Radius1 = sum(State_VIB(X_:Z_, iAbove,   iBlock)**2)**0.5
         Weight  = (File_I(iFile)%Radius - Radius0) / (Radius1 - Radius0)
         ! interpolate each requested variable
         do iVarPlot = 1, nVarPlot
            iVarIndex = File_I(iFile) % iVarPlot_V(iVarPlot)
            File_I(iFile) % Buffer_II(iVarPlot, nDataLine) = &
                 State_VIB(iVarIndex, iAbove-1, iBlock) * (1-Weight) + &
                 State_VIB(iVarIndex, iAbove,   iBlock) *    Weight
         end do

         ! reprint data to file
         call save_plot_file(&
              NameFile     = NameFile, &
              TypeFileIn   = File_I(iFile) % TypeFile, &
              nDimIn       = 1, &
              TimeIn       = TimeGlobal, &
              nStepIn      = iIterGlobal, &
              CoordMinIn_D = (/real(iIterGlobal - nDataLine + 1)/), &
              CoordMaxIn_D = (/real(iIterGlobal)/), &
              NameVarIn    = File_I(iFile) % NameVarPlot, &
              VarIn_VI     = &
              File_I(iFile) % Buffer_II(1:nVarPlot,1:nDataLine)&
              )
      end do

      ! mark that the first call is done
      if(File_I(iFile)%IsFirstCall)&
           File_I(iFile)%IsFirstCall = .false.

    end subroutine write_mh_time

    !=========================================================================
    
    subroutine write_distr_1d
      ! write file with distribution in the format to be read by IDL/TECPLOT;
      ! separate file is created for each field line, name format is
      ! Distribution_<iLon>_<iLat>_t<ddhhmmss>_n<iIter>.{out/dat}
      !------------------------------------------------------------------------
      ! name of the output file
      character(len=100):: NameFile
      ! loop variables
      integer:: iBlock, iParticle, iVarPlot
      ! indexes of corresponding node, latitude and longitude
      integer:: iNode, iLat, iLon
      ! index of first/last particle on the field line
      integer:: iFirst, iLast
      ! scale and conversion factor
      real:: Scale_I(nMomentumBin), Factor_I(nMomentumBin)
      real:: Unity_I(nMomentumBin) = 1.0
      ! timestamp
      character(len=8):: StringTime
      !------------------------------------------------------------------------
      select case(File_I(iFile) % iScale)
      case(MomentumScale_)
         Scale_I = LogMomentumScale_I / log(10.)
         Factor_I= Unity_I
      case(EnergyScale_)
         Scale_I = LogEnergyScale_I / log(10.)
         Factor_I= DMomentumOverDEnergy_I
      end select
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! set the file name
         call get_time_string(TimeGlobal, StringTime)
         write(NameFile,'(a,i3.3,a,i3.3,a,i6.6,a)') &
              trim(NamePlotDir)//'Distribution_',iLon,'_',iLat,&
              '_t'//StringTime//'_n',iIterGlobal,&
              File_I(iFile) % NameFormat

         ! get min and max particle indexes on this field line
         iFirst = iGridLocal_IB(Begin_, iBlock)
         iLast  = iGridLocal_IB(End_,   iBlock)
         
         do iParticle = iParticleMin, iParticleMax
            ! reset values outside the line's range
            if(iParticle < iFirst .or. iParticle > iLast)then
               File_I(iFile) % Buffer_II(:,iParticle) = 0.0
               CYCLE
            end if
            ! the actual distribution
            File_I(iFile) % Buffer_II(:,iParticle) = &
                 log10(Distribution_IIB(:,iParticle,iBlock) * Factor_I(:))
         end do

         ! print data to file
         call save_plot_file(&
              NameFile   = NameFile, &
              TypeFileIn = File_I(iFile) % TypeFile, &
              nDimIn     = 2, &
              TimeIn     = TimeGlobal, &
              nStepIn    = iIterGlobal, &
              Coord1In_I = Scale_I, &
              Coord2In_I = State_VIB(S_,iFirst:iLast,iBlock), &
              NameVarIn  = File_I(iFile) % NameVarPlot, &
              VarIn_II   = File_I(iFile) % Buffer_II(:,iFirst:iLast) &
              )
      end do     
    end subroutine write_distr_1d

  end subroutine write_output

  !==========================================================================

  subroutine get_time_string(Time, StringTime)
    ! the subroutine converts real variable Time into a string,
    ! the structure of the string is 'ddhhmmss', 
    ! i.e shows number of days, hours, minutes and seconds 
    ! after the beginning of the simulation
    real,             intent(in) :: Time
    character(len=8), intent(out):: StringTime
    !--------------------------------------------------------------------------
    ! This is the value if the time is too large
    StringTime = '99999999'
    if(Time < 100.0*86400) &
         write(StringTime,'(i2.2,i2.2,i2.2,i2.2)') &
         int(                  Time          /86400.), & ! # days
         int((Time-(86400.*int(Time/86400.)))/ 3600.), & ! # hours
         int((Time-( 3600.*int(Time/ 3600.)))/   60.), & ! # minutes
         int( Time-(   60.*int(Time/   60.)))            ! # seconds
  end subroutine get_time_string

end module SP_ModWrite
