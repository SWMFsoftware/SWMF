module SP_ModWrite

  ! This module contains methods for writing output files

  use SP_ModSize, ONLY: &
       nDim, nLat, nLon, nNode, &
       iParticleMin, iParticleMax, nParticle,&
       nMomentumBin, &
       Particle_, OriginLat_, OriginLon_

  use SP_ModGrid, ONLY: &
       get_node_indexes, &
       iComm, &
       nVar, nBlock, State_VIB, iGridLocal_IB, iNode_B, &
       Distribution_IIB, LogEnergyScale_I, LogMomentumScale_I, &
       DMomentumOverDEnergy_I, &
       Proc_, Begin_, End_, X_, Y_, Z_, Bx_, By_, Bz_, &
       B_, Ux_, Uy_, Uz_, U_, Rho_, T_, S_, EFlux_, &
       NameVar_V

  use ModPlotFile, ONLY: save_plot_file

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
       MH3D_    = 0, & ! Background mhd data along all lines
       MHSph_   = 1, & ! Background mhd data on a given sphere
       Distr3D_ = 2, & ! Distribution along all lines
       DistrSph_= 3    ! Distribution on a given sphere
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
       if(    index(StringPlot,'mh3d' )    > 0)then
          File_I(iFile) % iKindData = MH3D_
       elseif(index(StringPlot,'mhsph')    > 0)then
          File_I(iFile) % iKindData = MHSph_
       elseif(index(StringPlot,'distr3d')  > 0)then
          File_I(iFile) % iKindData = Distr3D_
       elseif(index(StringPlot,'distrsph') > 0)then
          File_I(iFile) % iKindData = DistrSph_
          call CON_stop(NameSub//&
               ': spherical output for distribution is not implemented yet')
       else
          call CON_stop(NameSub//': kind of data was not set in PARAM.in')
       end if

       ! reset variables' names
       File_I(iFile) % NameVarPlot = ''

       ! based on kind of data process the requested output
       select case(File_I(iFile) % iKindData)
          case(MH3D_)
             call process_mh
             ! prepare the output data container
             allocate(File_I(iFile) % Buffer_II(&
                  File_I(iFile)%nVarPlot, iParticleMin:iParticleMax))
             ! add particle index to variable names
             File_I(iFile) % NameVarPlot = &
                  'ParticleIndex '//trim(File_I(iFile) % NameVarPlot)
          case(MHSph_)
             call process_mh
             ! prepare the output data container
             allocate(File_I(iFile) % Buffer_II(&
                  File_I(iFile)%nVarPlot, nNode))
             ! add line index to variable names
             File_I(iFile) % NameVarPlot = &
                  'LineIndex '//trim(File_I(iFile) % NameVarPlot)
             ! get radius
             call read_var('Radius [Rs]', File_I(iFile) % Radius)
          case(Distr3D_)
             call process_distr
             ! prepare the output data container
             allocate(File_I(iFile) % &
                  Buffer_II(nMomentumBin,iParticleMin:iParticleMax))
          case(DistrSph_)
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
      integer:: iVar, iVarPlot
      !----------------------
      ! reset
      File_I(iFile) % nVarPlot = 0
      File_I(iFile) % DoPlot_V = .false.
      ! coordinates are always printed
      File_I(iFile) % DoPlot_V((/X_, Y_, Z_/)) = .true.
      File_I(iFile) % nVarPlot = File_I(iFile) % nVarPlot + 3
      ! distance along line -----
      if(index(StringPlot,' dist ') > 0)then
         File_I(iFile) % DoPlot_V(S_) = .true.
         File_I(iFile) % nVarPlot = File_I(iFile) % nVarPlot + 1
      end if
      ! plasma density ----------
      if(index(StringPlot,' rho ') > 0)then
         File_I(iFile) % DoPlot_V(Rho_) = .true.
         File_I(iFile) % nVarPlot = File_I(iFile) % nVarPlot + 1
      end if
      ! temperature -------------
      if(index(StringPlot,' temp ')> 0)then
         File_I(iFile) % DoPlot_V(T_) = .true.
         File_I(iFile) % nVarPlot = File_I(iFile) % nVarPlot + 1
      end if
      ! velocity ----------------
      if(index(StringPlot,' ux ')> 0 .or. index(StringPlot,' uy ')> 0 .or. &
           index(StringPlot,' uz ')> 0)then
         File_I(iFile) % DoPlot_V((/Ux_, Uy_, Uz_/)) = .true.
         File_I(iFile) % nVarPlot = File_I(iFile) % nVarPlot + 3
      end if
      if(index(StringPlot,' |u| ')> 0)then
         File_I(iFile) % DoPlot_V(U_) = .true.
         File_I(iFile) % nVarPlot = File_I(iFile) % nVarPlot + 1
      end if
      ! magnetic field ----------
      if(index(StringPlot,' bx ')> 0 .or. index(StringPlot,' by ')> 0 .or. &
           index(StringPlot,' bz ')> 0)then
         File_I(iFile) % DoPlot_V((/Bx_, By_, Bz_/)) = .true.
         File_I(iFile) % nVarPlot = File_I(iFile) % nVarPlot + 3
      end if
      if(index(StringPlot,' |b| ')> 0)then
         File_I(iFile) % DoPlot_V(B_) = .true.
         File_I(iFile) % nVarPlot = File_I(iFile) % nVarPlot + 1
      end if
      ! energy flux -------------
      if(index(StringPlot,' eflux ')> 0)then
         File_I(iFile) % DoPlot_V(EFlux_) = .true.
         File_I(iFile) % nVarPlot = File_I(iFile) % nVarPlot + 1
      end if
      !--------------------------
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
              'LogMomentum Distance LogDistribution'
      elseif(index(StringPlot, 'energy') > 0)then
         File_I(iFile) % iScale = EnergyScale_
         File_I(iFile) % NameVarPlot = &
              'LogEnergy Distance LogDistribution'
      else
         call CON_stop(NameSub//&
              ': type of scale for distribution output wasnot set in PARAM.in')
      end if
    end subroutine process_distr
  end subroutine set_write_param

  !============================================================================

  subroutine write_output(Time, iIter)
    ! write the output data
    real,    intent(in):: Time ! current time
    integer, intent(in):: iIter! current iteration

    ! loop variables
    integer:: iFile

    character(len=*), parameter:: NameSub = 'SP:write_output'
    !--------------------------------------------------------------------------
    if(nFile == 0) RETURN

    do iFile = 1, nFile
       select case(File_I(iFile) % iKindData)
       case(MH3D_)
          call write_mh_3d
       case(MHSph_)
          call write_mh_sph
       case(Distr3D_)
          call write_distr_3d
       case(DistrSph_)
          !call write_distr_sph
       end select
    end do

  contains

    subroutine write_mh_3d
      ! write output with 3D MH data in the format to be read by IDL/TECPLOT;
      ! separate file is created for each field line, name format is
      ! MH_data_<iLon>_<iLat>_n<iIter>.{out/dat}
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
      !------------------------------------------------------------------------
      nVarPlot = File_I(iFile) % nVarPlot
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! set the file name
         write(NameFile,'(a,i3.3,a,i3.3,a,i8.8,a,i6.6,a)') &
              trim(NamePlotDir)//'MH_data_',iLon,'_',iLat,&
              '_t',floor(Time),'_n',iIter,&
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

         ! print data to file
         call save_plot_file(&
              NameFile     = NameFile, &
              TypeFileIn   = File_I(iFile) % TypeFile, &
              nDimIn       = 1, &
              TimeIn       = Time, &
              nStepIn      = iIter, &
              CoordMinIn_D = (/real(iFirst)/), &
              CoordMaxIn_D = (/real(iLast)/), &
              NameVarIn    = File_I(iFile) % NameVarPlot, &
              VarIn_VI     = &
              File_I(iFile) % Buffer_II(1:nVarPlot,iFirst:iLast)&
              )
      end do
    end subroutine write_mh_3d

    !=========================================================================

    subroutine write_mh_sph
      use ModMpi
      use CON_world, ONLY: is_proc0, i_proc0
      use CON_comp_param, ONLY: SP_
      ! write output with 3D MH data in the format to be read by IDL/TECPLOT;
      ! single file is created for all field lines, name format is
      ! MH_data_R=<Radius [AU]>_n<iIter>.{out/dat}
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
         write(NameFile,'(a,i4.4,f0.2,a,i8.8,a,i6.6,a)') &
              trim(NamePlotDir)//'MH_data_R=', int(File_I(iFile) % Radius), &
              File_I(iFile) % Radius - int(File_I(iFile) % Radius), &
              '_t',floor(Time),'_n', iIter, File_I(iFile) % NameFormat

         ! get min and max particle indexes on this field line
         iFirst = iGridLocal_IB(Begin_, iBlock)
         iLast  = iGridLocal_IB(End_,   iBlock)

         ! find the particle just above the given radius
         do iParticle = iFirst , iLast
            Radius0 = sum(State_VIB(X_:Z_, iParticle, iBlock)**2)**0.5
            if( Radius0 > File_I(iFile) % Radius**2) EXIT
            ! check if reached the end, i.e. there is no intersection
            if(iParticle == iLast) &
                 DoPrint_I(iNode) = .false.
         end do
         !check if field line started above output sphere, i.e. no intersection
         if(iParticle == iFirst) &
              DoPrint_I(iNode) = .false.

         ! if no intersection found -> proceed to the next line
         if(.not.DoPrint_I(iNode)) CYCLE

         !intersection is found -> get data at that location
         iAbove = iParticle
         
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

      if(is_proc0(SP_))&
           ! print data to file
           call save_plot_file(&
           NameFile     = NameFile, &
           TypeFileIn   = File_I(iFile) % TypeFile, &
           nDimIn       = 1, &
           TimeIn       = Time, &
           nStepIn      = iIter, &
           Coord1In_I   = real(pack(iNodeIndex_I, MASK=DoPrint_I)), &
           NameVarIn    = File_I(iFile) % NameVarPlot, &
           VarIn_VI     = &
           reshape(&
           pack(File_I(iFile) % Buffer_II(1:nVarPlot,1:nNode),&
           MASK = spread(DoPrint_I, 1, nVarPlot)), &
           (/nVarPlot, count(DoPrint_I)/))&
           )
    end subroutine write_mh_sph

    !=========================================================================
    
    subroutine write_distr_3d
      ! write file with distribution in the format to be read by IDL/TECPLOT;
      ! separate file is created for each field line, name format is
      ! Distribution_<iLon>_<iLat>_n<iIter>.{out/dat}
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
      real, pointer:: Scale_I(:), Factor_I(:)
      real, target:: Unity_I(nMomentumBin) = 1.0
      !------------------------------------------------------------------------
      select case(File_I(iFile) % iScale)
      case(MomentumScale_)
         Scale_I => LogMomentumScale_I
         Factor_I=> Unity_I
      case(EnergyScale_)
         Scale_I => LogEnergyScale_I
         Factor_I=> DMomentumOverDEnergy_I
      end select
      do iBlock = 1, nBlock
         iNode = iNode_B(iBlock)
         call get_node_indexes(iNode, iLon, iLat)

         ! set the file name
         write(NameFile,'(a,i3.3,a,i3.3,a,i8.8,a,i6.6,a)') &
              trim(NamePlotDir)//'Distribution_',iLon,'_',iLat,&
              '_t',floor(Time),'_n',iIter,&
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
                 log(Distribution_IIB(:,iParticle,iBlock) * Factor_I(:))
         end do

         ! print data to file
         call save_plot_file(&
              NameFile   = NameFile, &
              TypeFileIn = File_I(iFile) % TypeFile, &
              nDimIn     = 2, &
              TimeIn     = Time, &
              nStepIn    = iIter, &
              Coord1In_I = Scale_I, &
              Coord2In_I = State_VIB(S_,iFirst:iLast,iBlock), &
              NameVarIn  = File_I(iFile) % NameVarPlot, &
              VarIn_II   = File_I(iFile) % Buffer_II(:,iFirst:iLast) &
              )
      end do     
    end subroutine write_distr_3d

  end subroutine write_output

end module SP_ModWrite
