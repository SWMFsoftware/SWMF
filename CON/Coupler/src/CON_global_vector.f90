!\documentstyle[12pt]{letter}      !^CFG UNCOMMENT IF PRINTPS   !
!\setlength{\textwidth}{6.5 in}    !^CFG UNCOMMENT IF PRINTPS   !
!\setlength{\oddsidemargin}{0.5 in}!^CFG UNCOMMENT IF PRINTPS   !
!\begin{document}                  !^CFG UNCOMMENT IF PRINTPS   !
!\begin{verbatim}                  !^CFG UNCOMMENT IF PRINTPS   !
!^CFG COPYRIGHT UM                                              !
!BOP
!MODULE: CON_global_vector - pointwise state vectors
!INTERFACE:
Module CON_global_vector
  !DESCRIPTION:
  !The toolkit for coupling the different data sets within a      
  !single component code, or within a framework as it is now (see 
  !the date below) consists of: 
  !                                  
  !ModMpi (iclude for mpif90.h)
  !                                   
  !ModNumConst
  !                                                    
  !CON\_world
  !                                                                   
  !CON\_domain\_decomposition
  !                                     
  !CON\_grid\_descriptor
  !                            
  !CON\_global_vector
  !                                          
  !CON\_router
  !                                                    
  !CON\_global\_message\_pass
  !
  !CON_global_vector handles with the poitntwise state vectors, 
  !first index of the array enumerates the state component,
  !second is implied to be the global cell index in a grid 
  !descriptor 


  !USES:

  use CON_grid_descriptor
  !REVISION HISTORY:
  ! Sokolov I.V.                                                  
  ! 6.18.04                                              
  ! igorsok@umich.edu                                             
  ! phone(734)647-4705                                            
  !EOP
  implicit none
  !BOP
  !DESCRIPTION:
  !
  !\begin{verbatim}                                
  !=====================DERIVED TYPE==========================!
  private !except
  type GlobalVectorType
     real,dimension(:,:),pointer::I
  end type GlobalVectorType
  type GlobalMaskType
     logical,dimension(:),pointer::I
  end type GlobalMaskType
  !EOP
  integer,parameter::nVectorMax=20,nLength=10
  integer,save::nVector=0,nMask=0

  type(GlobalVectorType),dimension(nVectorMax),save::&
       Vector_I
  character(LEN=nLength),dimension(nVectorMax),save::&
       NameVector_I


  type(GlobalMaskType),dimension(nVectorMax),save::&
       Mask_I
  character(LEN=nLength),dimension(nVectorMax),save::&
       NameMask_I

  public::allocate_vector,allocate_mask
  interface allocate_vector
     module procedure allocate_vector_ilength
     module procedure allocate_vector_for_gd
     module procedure allocate_vector_for_mask
  end interface
  interface allocate_mask
     module procedure allocate_mask_ilength
     module procedure allocate_mask_for_gd
     module procedure allocate_mask_for_vector
  end interface
  public::deallocate_vector,deallocate_mask
  public::associate_with_global_vector
  public::associate_with_global_mask
  interface bcast_global_vector
     module procedure bcast_in_comm
     module procedure bcast_in_grid
     module procedure bcast_in_union
  end interface
  public::bcast_global_vector
  public::used_vector,used_mask,ubound_vector,ubound_mask
  public::set_mask,count_mask
  public::point_state_v
  interface point_state_v
     module procedure point_state_vi
     module procedure point_state_vx
  end interface
  public:: save_global_vector,read_global_vector
contains
  !===========================================================!
  subroutine allocate_vector_ilength(NameVector,nVar,nI)
    character(LEN=*),intent(in)::NameVector
    integer,intent(in)::nVar,nI
    integer:: iError
    nVector=nVector+1
    if(nVector>nVectorMax)call CON_stop(&
         'Insufficient nVectorMax in CON_global_vector')
    NameVector_I(nVector)=trim(NameVector)
    nullify(Vector_I(nVector)%I)
    allocate(Vector_I(nVector)%I(nVar,nI),stat=iError)
    call check_allocate(iError,'Global vector '//NameVector)
    Vector_I(nVector)%I=cZero
  end subroutine allocate_vector_ilength
  !===========================================================!
  subroutine allocate_vector_for_gd(NameVector,nVar,GD)
    character(LEN=*),intent(in)::NameVector
    integer,intent(in)::nVar
    type(GridDescriptorType),intent(in)::GD
    call allocate_vector_ilength(NameVector,nVar,&
         n_block_total(GD%DD%Ptr)*n_grid_points_per_block(GD))
  end subroutine allocate_vector_for_gd
  !===========================================================!
  subroutine allocate_vector_for_mask(NameVector,nVar,NameMask)
    character(LEN=*),intent(in)::NameVector,NameMask
    integer,intent(in)::nVar
    call allocate_vector_ilength(NameVector,nVar,&
         ubound(Mask_I(i_mask(NameMask))%I,1))
  end subroutine allocate_vector_for_mask
  !===========================================================!
  subroutine allocate_mask_ilength(NameMask,nI)
    character(LEN=*),intent(in)::NameMask
    integer,intent(in)::nI
    integer:: iError
    nMask=nMask+1
    NameMask_I(nMask)=trim(NameMask)
    nullify(Mask_I(nMask)%I)
    allocate(Mask_I(nMask)%I(nI),stat=iError)
    call check_allocate(iError,'Global mask '//NameMask)
    Mask_I(nMask)%I=.true.
  end subroutine allocate_mask_ilength
  !===========================================================!
  subroutine allocate_mask_for_gd(NameMask,GD)
    character(LEN=*),intent(in)::NameMask
    type(GridDescriptorType),intent(in)::GD
    call allocate_mask_ilength(NameMask,&
         n_block_total(GD%DD%Ptr)*n_grid_points_per_block(GD))
  end subroutine allocate_mask_for_gd
  !==========================================================!
  subroutine allocate_mask_for_vector(NameMask,NameVector)
    character(LEN=*),intent(in)::NameMask,NameVector
    call allocate_mask_ilength(NameMask,&
         ubound(Vector_I(i_vector(NameVector))%I,2))
  end subroutine allocate_mask_for_vector
  !==========================================================!
  integer function i_vector(Name)
    character(LEN=*),intent(in)::Name
    integer::iLoop
    do iLoop=1,nVector
       i_vector=iLoop
       if(trim(Name)==trim(NameVector_I(iLoop)))return
    end do
    call CON_stop('Global vector '//Name//' is not allocated')
  end function i_vector
  !==========================================================!
  integer function i_mask(Name)
    character(LEN=*),intent(in)::Name
    integer::iLoop
    do iLoop=1,nMask
       i_mask=iLoop
       if(trim(Name)==trim(NameMask_I(iLoop)))return
    end do
    call CON_stop('Global mask '//Name//' is not allocated')
  end function i_mask
  !==========================================================!
  subroutine deallocate_vector(NameVector)
    character(LEN=*),intent(in)::NameVector
    integer::iVector,iLoop
    iVector=i_vector(NameVector)
    deallocate (Vector_I(iVector)%I)
    do iLoop=iVector,nVector-1
       nullify(Vector_I(iLoop)%I)
       Vector_I(iLoop)%I=>Vector_I(iLoop+1)%I
       NameVector_I(iLoop)=NameVector_I(iLoop+1)
    end do
    nullify(Vector_I(nVector)%I)
    NameVector_I(nVector)=''
    nVector=nVector-1
  end subroutine deallocate_vector
  !===========================================================!
  subroutine deallocate_mask(NameMask)
    character(LEN=*),intent(in)::NameMask
    integer::iMask,iLoop
    iMask=i_mask(NameMask)
    deallocate (Mask_I(iMask)%I)
    do iLoop=iMask,nMask-1
       nullify(Mask_I(iLoop)%I)
       Mask_I(iLoop)%I=>Mask_I(iLoop+1)%I
       NameMask_I(iLoop)=NameMask_I(iLoop+1)
    end do
    nullify(Mask_I(nMask)%I)
    NameMask_I(nMask)=''
    nMask=nMask-1
  end subroutine deallocate_mask
  !===========================================================!
  subroutine associate_with_global_vector(Ptr,Name)
    real,dimension(:,:),pointer::Ptr
    character(LEN=*),intent(in)::Name
    integer :: iTemp
    iTemp=i_vector(Name)
    nullify(Ptr)
    Ptr=>Vector_I(iTemp)%I
  end subroutine associate_with_global_vector
  !===========================================================!
  subroutine associate_with_global_mask(Ptr,Name)
    logical,dimension(:),pointer::Ptr
    character(LEN=*),intent(in)::Name
    integer :: iTemp
    iTemp=i_mask(Name)
    nullify(Ptr)
    Ptr=>Mask_I(iTemp)%I   
  end subroutine associate_with_global_mask
  !===========================================================!
  logical function used_vector(Name)
    character(LEN=*),intent(in)::Name
    integer::iVector
    used_vector=.false.
    do iVector=1,nVector
       used_vector=&
            trim(NameVector_I(iVector))==trim(Name)
       if(used_vector)return
    end do
  end function used_vector
  !===========================================================!
  logical function used_mask(Name)
    character(LEN=*),intent(in)::Name
    integer::iMask
    used_mask=.false.
    do iMask=1,nMask
       used_mask=trim(NameMask_I(iMask))==trim(Name)
       if(used_mask)return
    end do
  end function used_mask
  !===========================================================!
  subroutine bcast_in_comm(&
       NameVector,iProc,iComm,NameMask)
    character(LEN=*),intent(in)::NameVector,NameMask
    optional::NameMask
    integer,intent(in)::iProc,iComm
    integer::iError,iVector,iMask,iLoop,iStart,nU_I(2)
    iVector=i_vector(NameVector)
    nU_I=ubound(Vector_I(iVector)%I)
    if(present(NameMask))then
       iMask=i_mask(NameMask)
       iLoop=0
       FINDTRUE:do while (iLoop<nU_I(2))
          iLoop=iLoop+1
          if(.not.Mask_I(iMask)%I(iLoop))CYCLE FINDTRUE
          iStart=iLoop
          FINDFALSE:do while (iLoop<nU_I(2))
             if(.not.Mask_I(iMask)%I(iLoop+1))EXIT FINDFALSE
             iLoop=iLoop+1
          end do FINDFALSE
          call MPI_BCAST(Vector_I(iVector)%I(1,iStart),&
               (iLoop-iStart+1)*nU_I(1),&
               MPI_REAL,&
               iProc,iComm,iError)
       end do FINDTRUE
    else
       call MPI_BCAST(Vector_I(iVector)%I(1,1),&
            nU_I(1)*nU_I(2),&
            MPI_REAL,&
            iProc,iComm,iError)
    end if
  end subroutine bcast_in_comm
  !=======================================================!
  subroutine bcast_in_grid(&
       NameVector,GD,NameMask)
    character(LEN=*),intent(in)::NameVector,NameMask
    optional::NameMask
    type(GridDescriptorType),intent(in)::GD
    integer::iError,iVector,iMask,iLoop
    integer::iStart,iFinal
    integer::iBlockAll,nBlockAll
    integer::nGridPointsPerBlock
    integer::iProc,iComm,iProcNew
    integer::nU_I(2)
    logical::UseMask
    !----------------
    if(.not.GD%DD%Ptr%IsLocal)&
         call CON_stop(&
         'Use bcast_in_router if a global vector has a global GD')
    nGridPointsPerBlock=n_grid_points_per_block(GD)
    nBlockAll=n_block_total(GD%DD%Ptr)
    iComm=i_comm(compid_grid(GD%DD%Ptr))
    iVector=i_vector(NameVector)
    UseMask=present(NameMask)
    if(UseMask)iMask=i_mask(NameMask)
    nU_I=ubound(Vector_I(iVector)%I)
    if(nU_I(2)/=nBlockAll*nGridPointsPerBlock)&
         call CON_stop(&
         'The gird descriptor does not fit the global vector '//&
         NameVector)
    iBlockAll=0
    iFinal=0

    !Start PE is used for the first block
    iProcNew=pe_decomposition(GD%DD%Ptr,&
         i_global_node_a(GD%DD%Ptr,1))

    DIFFBLOCK:do while(iBlockAll<nBlockAll)
      
       iBlockAll=iBlockAll+1
       iStart=iFinal+1
       iFinal=iFinal+nGridPointsPerBlock
       iProc= iProcNew
       SAMEBLOCK:do while(iBlockAll<nBlockAll)
          !Find a sequence of block set at the same
          !PE = iProc

          iProcNew=pe_decomposition(&
               GD%DD%Ptr,&
               i_global_node_a(GD%DD%Ptr,iBlockAll+1))
          if(iProc/=iProcNew)EXIT SAMEBLOCK
          iBlockAll=iBlockAll+1
          iFinal=iFinal+nGridPointsPerBlock
       end do SAMEBLOCK

       if(UseMask)then
          iLoop=iStart-1

          FINDTRUE:do while (iLoop<iFinal)

             !Find the start of a sequence of true points

             iLoop=iLoop+1
             if(.not.Mask_I(iMask)%I(iLoop))CYCLE FINDTRUE
             iStart=iLoop
             FINDFALSE:do while (iLoop<iFinal)

                !Find the end of a sequence of true points

                if(.not.Mask_I(iMask)%I(iLoop+1))EXIT FINDFALSE
                iLoop=iLoop+1
             end do FINDFALSE

             call MPI_BCAST(Vector_I(iVector)%I(1,iStart),&
                  (iLoop-iStart+1)*nU_I(1),&
                  MPI_REAL,&
                  iProc,iComm,iError)
          end do FINDTRUE
       else
          call MPI_BCAST(Vector_I(iVector)%I(1,iStart),&
               (iFinal-iStart+1)*nU_I(1),&
               MPI_REAL,&
               iProc,iComm,iError)
       end if
    end do DIFFBLOCK
  end subroutine bcast_in_grid
  !=======================================================!
  subroutine bcast_in_union(&
       NameVector,GD,iTranslated_P,iCommUnion,NameMask)
    character(LEN=*),intent(in)::NameVector,NameMask
    optional::NameMask
    type(GridDescriptorType),intent(in)::GD
    integer,dimension(:),intent(in)::&
         iTranslated_P
    integer,intent(in)::iCommUnion
    integer::iError,iVector,iMask,iLoop
    integer::iStart,iFinal
    integer::iBlockAll,nBlockAll
    integer::nGridPointsPerBlock
    integer::iProc,iProcNew
    integer::nU_I(2)
    logical::UseMask
    integer::iPETransFrom_I(1),iPETransTo_I(1)
    nGridPointsPerBlock=n_grid_points_per_block(GD)
    nBlockAll=n_block_total(GD%DD%Ptr)

    iVector=i_vector(NameVector)
    UseMask=present(NameMask)
    if(UseMask)iMask=i_mask(NameMask)
    nU_I=ubound(Vector_I(iVector)%I)
    if(nU_I(2)/=nBlockAll*nGridPointsPerBlock)&
         call CON_stop(&
         'The gird descriptor does not fit the global vector '//&
         NameVector)
    iBlockAll=0
    iFinal=0
    iProcNew=pe_decomposition(GD%DD%Ptr,&
         i_global_node_a(GD%DD%Ptr,1))
    DIFFPROC:do while(iBlockAll<nBlockAll)
       iBlockAll=iBlockAll+1
       iStart=iFinal+1
       iFinal=iFinal+nGridPointsPerBlock
       iProc= iProcNew
       SAMEPROC:do while(iBlockAll<nBlockAll)
          iProcNew=pe_decomposition(&
               GD%DD%Ptr,&
               i_global_node_a(GD%DD%Ptr,iBlockAll+1))
          if(iProc/=iProcNew)EXIT SAMEPROC
          iBlockAll=iBlockAll+1
          iFinal=iFinal+nGridPointsPerBlock
       end do SAMEPROC
       if(UseMask)then
          iLoop=iStart-1
          FINDTRUE:do while (iLoop<iFinal)
             iLoop=iLoop+1
             if(.not.Mask_I(iMask)%I(iLoop))CYCLE FINDTRUE
             iStart=iLoop
             FINDFALSE:do while (iLoop<iFinal)
                if(Mask_I(iMask)%I(iLoop+1))EXIT FINDFALSE
                iLoop=iLoop+1
             end do FINDFALSE
             call MPI_BCAST(Vector_I(iVector)%I(1,iStart),&
                  (iLoop-iStart+1)*nU_I(1),&
                  MPI_REAL,&
                  iTranslated_P(iProc),iCommUnion,iError)
          end do FINDTRUE
       else
          call MPI_BCAST(Vector_I(iVector)%I(1,iStart),&
               (iFinal-iStart+1)*nU_I(1),&
               MPI_REAL,&
               iTranslated_P(iProc),iCommUnion,iError)
       end if
    end do DIFFPROC
  end subroutine bcast_in_union
  !=======================================================
  function ubound_vector(NameVector)
    integer,dimension(2)::ubound_vector
    character(LEN=*),intent(in)::NameVector
    ubound_vector=ubound(Vector_I(i_vector(NameVector))%I)
  end function ubound_vector
  !=======================================================
  integer function ubound_mask(NameMask)
    character(LEN=*),intent(in)::NameMask
    ubound_mask=ubound(Mask_I(i_mask(NameMask))%I,1)
  end function ubound_mask
  !=======================================================
  subroutine set_mask(NameMask,NameVector,used)
    character(LEN=*),intent(in)::NameMask,NameVector
    interface
       logical function used(State_V)
         implicit none
         real,dimension(:),intent(in)::State_V
       end function used
    end interface
    integer::iPoint,nPoint,iVector,iMask
    iVector=i_vector(NameVector)
    iMask=i_mask(NameMask)
    nPoint=ubound(Mask_I(iMask)%I,1)
    do iPoint=1,nPoint
       Mask_I(iMask)%I(iPoint)=&
            used(Vector_I(iVector)%I(:,iPoint))
    end do
  end subroutine set_mask
  !=======================================================
  integer function count_mask(NameMask)
    character(LEN=*),intent(in)::NameMask
    integer :: iTemp
    iTemp=i_mask(NameMask)
    count_mask=count(Mask_I(iTemp)%I)
  end function count_mask
  !======================================================!
  function point_state_vi(NameVector,nVar,iPoint)
    character(LEN=*),intent(in)::NameVector
    integer,intent(in)::iPoint,nVar
    real,dimension(nVar)::point_state_vi
    integer,save::iVector=0
    character(LEN=nLength),save::NameSaved=''
    if(trim(NameSaved)/=trim(NameVector))then
       iVector=i_vector(NameVector)
       NameSaved=NameVector
    end if
    point_state_vi=Vector_I(iVector)%I(:,iPoint)
  end function point_state_vi
  !======================================================!
  function point_state_vx(&
       NameVector, &
       nVar,       &
       nDim,       &
       Xyz_D,      &
       GD,         &
       interpolate)
    type(GridDescriptorType),intent(in):: GD
    integer,intent(in)::nDim
    real,intent(in):: Xyz_D(nDim)
    integer,intent(in)::nVar
    character(LEN=*),intent(in)::NameVector
    optional::interpolate
    interface
       subroutine interpolate(&
            nDim,&
            Xyz_D,&
            GridDescriptor,&
            nIndexes,&
            Index_II,&
            nImages,Weight_I)
         use CON_grid_descriptor
         implicit none
         integer,intent(in)::nDim
         real,intent(inout)::Xyz_D(nDim) 
         type(GridDescriptorType)::GridDescriptor     
         integer,intent(in)::nIndexes
         integer           ::Index_II(0:nIndexes,2**nDim)
         integer,intent(out)::nImages
         real,intent(out)::Weight_I(2**nDim)
       end subroutine interpolate
    end interface
    !Return value:
    real,dimension(nVar)::point_state_vx
    !Local variables:
    real,dimension(nDim)::XyzMisc_D 
    integer::nIndexes,iPoint
    integer,dimension(0:nDim+1,2**nDim)::Index_II
    integer::nImages,iImages,lGlobalNode
    real,dimension(2**nDim)::Weight_I
    real,dimension(nVar)::State_V
    !---------------------------------------------------------
    nIndexes=nDim+1;XyzMisc_D=Xyz_D
    if(present(interpolate))then
       call interpolate(&
            nDim,&
            XyzMisc_D,&
            GD,&
            nIndexes,&
            Index_II,&
            nImages,Weight_I)
    else
       call nearest_grid_points(&
            nDim,&
            XyzMisc_D,&
            GD,&
            nIndexes,&
            Index_II,&
            nImages,Weight_I)
    end if
    point_state_vx=cZero
    do iImages=1,nImages
       lGlobalNode=i_global_node_bp(GD%DD%Ptr,&
            Index_II(nIndexes,iImages),&
            Index_II(0,iImages))
       iPoint=i_grid_point_global(GD,&
            lGlobalNode,&
            Index_II(1:nDim,iImages))
       point_state_vx=point_state_vx+&
            point_state_vi(NameVector,nVar,iPoint)&
            *Weight_I(iImages)
    end do
  end function point_state_vx
  subroutine save_global_vector(NameVector,NameMask,iFileIn)
    use ModIOUnit
    character(LEN=*),intent(in)::NameVector,NameMask
    integer,intent(in)::iFileIn
    optional::NameMask,iFileIn
    character(LEN=nLength+11)::NameFile
    character(LEN=2)::NameComp
    logical::UseMask
    integer::lLengthMask,iFile,nU_I(2),iPoint,iVector,iMask
    NameComp=NameVector(1:2)
    NameFile='./'//NameComp//'/'//NameVector(&
         4:len_trim(NameVector))
    if(present(iFileIn))&
         write(NameFile,'(a,i5.5)')trim(NameFile)//'_',iFileIn
    iVector=i_vector(NameVector)
    iFile=io_unit_new()
    nU_I=ubound(Vector_I(iVector)%I)
  
    UseMask=present(NameMask)
    if(UseMask)then
       lLengthMask=len_trim(NameMask)
       NameComp=NameMask(-1+lLengthMask:lLengthMask)
       iMask=i_mask(NameMask)
    end if
    NameFile=trim(NameFile)//'_'//NameComp
    open(iFile,file=trim(NameFile),status='replace')

    if(UseMask)then
       write(iFile,*)nU_I,count_mask(NameMask)
       do iPoint=1,nU_I(2)
          if(.not.Mask_I(iMask)%I(iPoint))CYCLE
          write(iFile,*)iPoint,Vector_I(iVector)%I(:,iPoint)
       end do
    else
       write(iFile,*)nU_I
       do iPoint=1,nU_I(2)
          write(iFile,*)iPoint,Vector_I(iVector)%I(:,iPoint)
       end do
    end if
    close(iFile)
  end subroutine save_global_vector
  subroutine read_global_vector(NameVector,iFileIn)
    use ModIOUnit
    use CON_comp_param,ONLY:NameComp_I
    character(LEN=*),intent(in)::NameVector
    integer,intent(in)::iFileIn
    optional::iFileIn
    integer::iFile,lComp,nU_I(2),nPoint,iVector,iError
    integer::iPoint,jPoint,nPointHere
    character(LEN=nLength+11)::NameFile,NameFilePreffix
    NameFilePreffix='./'//NameVector(1:2)//'/'//NameVector(&
         4:len_trim(NameVector))//'_'
    if(present(iFileIn))&
         write(NameFilePreffix,'(a,i5.5,a)')&
         trim(NameFilePreffix),iFileIn,'_'
    nPoint=0
    iFile=io_unit_new()
    do lComp=1,n_comp()
       NameFile=trim(NameFilePreffix)//NameComp_I(i_comp(lComp))
       open(iFile,FILE=NameFile,STATUS='old',IOSTAT=iError)
       if(iError>0)CYCLE
       if(NameFile(1:2)==NameComp_I(i_comp(lComp)))then
          read(iFile,*)nU_I
          nPointHere=nU_I(2)
       else
          read(iFile,*)nU_I,nPointHere
       end if
       if(.not.used_vector(NameVector))then
          call allocate_vector(NameVector,nU_I(1),nU_I(2))
          iVector=i_vector(NameVector)
       end if
       do jPoint=1,nPointHere
          read(iFile,*)iPoint,Vector_I(iVector)%I(:,iPoint)
          nPoint=nPoint+1
       end do
       close(iFile)
    end do
    if(nPoint/=nU_I(2))then
       write(*,*)'Vector '//NameVector//' with dimensions ',nU_I
       write(*,*)'nPoint=',nPoint
       call CON_stop('For vector '//NameVector//' fewer points are read')
    end if
  end subroutine read_global_vector
end Module CON_global_vector
      
