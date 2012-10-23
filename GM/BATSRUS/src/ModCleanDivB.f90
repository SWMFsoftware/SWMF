!This code is a copyright protected software (c) 2002- University of Michigan

Module ModDivbCleanup
  use ModSize
  real::OneDotMDotOneInv
  real, dimension(1:nI,1:nJ,1:nK,nBLK):: Prec_CB
  real:: BoundaryCoef = 1.0
  integer:: nCleanDivb = 0
end Module ModDivbCleanup
!==========================================================================
subroutine clean_divb
  use ModProcMH
  use ModSize
  use ModNumConst, ONLY: cTiny
  use ModDivbCleanup
  use ModMain,ONLY: iNewGrid, iNewDecomposition, nBlock, Unused_B
  use ModAdvance,ONLY: nVar,State_VGB, Bx_, By_, Bz_, P_,tmp1_BLK,tmp2_BLK,&
       Residual_GB=>tmp1_blk,Dir_GB=>tmp2_blk
  use ModAdvance,ONLY:tmp3_blk=>divB1_GB
  use ModGeometry, ONLY: true_blk, body_blk, true_cell
  use ModParallel, ONLY : NOBLK, neiLEV
  use ModPhysics,ONLY:gm1
  use ModMpi
  use BATL_lib, ONLY: message_pass_cell, CellFace_DB, CellVolume_GB

  implicit none

  integer::i,j,k,iBlock
  logical::DoConservative=.false.
  real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)::DivBV_G
  real:: DivBInt(2),DivBTemp(2)
  integer, parameter::ResDotM1DotRes_=2,ResDotOne_=1

  !Conjugated gradients algorithm, see
  !http://netlib2.cs.utk.edu/linalg/old_html_templates/subsection2.6.3.1.html
  !Notations
  !M^(-1)=diag(Prec_CB) - the Jacobi preconditioner
  !Residual_GB=M^{-1}r^(i)/rho_(i-1) (temporal)
  !Dir_GB=p^(i)/rho_(i-1)

  real::DirDotDir,DirDotDirInv
  !DirDotDirInv=1/(Dir.A.Dir)

  real::Tolerance=cTiny
  integer:: iError
  integer::Iteration                        

  real,dimension(nI,nJ,nK)::vDotGradX_C,vDotGradY_C,vDotgradZ_C
  !For convenience, what is stored in these arrays are the gradients
  !multiplied by the cell volume
  !While calculating the magnetic field, we use vInv*vDotGrad Phi
  !While calculating the sum V_i grad^2 Phi, we use vInv vDotGrad^2 Phi 


  logical::oktest,oktest_me

  integer :: iLastGrid=-100, iLastDecomposition=-100

  !---------------------------------------------------------------------------

  call set_oktest('clean_divb',oktest,oktest_me)

  if(iLastGrid /= iNewGrid .or. iLastDecomposition /= iNewDecomposition)then
     call init_divb_cleanup
     iLastGrid          = iNewGrid
     iLastDecomposition = iNewDecomposition
  endif
  call timing_start('clean_divb')

  Iteration=1                                

  do 
     call message_pass_cell(nVar, State_VGB, nWidthIn=1, &
          DoSendCornerIn=.false., nProlongOrderIn=1, DoRestrictFaceIn=.true.)

     !Get the ghostcell values for MF 
     DivBInt = 0.0;DivBTemp = 0.0
     do iBlock=1,nBlock
        Residual_GB(:,:,:,iBlock) = 0.0
        if(Iteration==1) Dir_GB(:,:,:,iBlock) = 0.0
        !Initialize the vector "Dir"

        if(Unused_B(iBlock))CYCLE
        call div_3d_b1(iBlock,&
             State_VGB(Bx_,:,:,:,iBlock),&
             State_VGB(By_,:,:,:,iBlock),&
             State_VGB(Bz_,:,:,:,iBlock),&
             DivBV_G)
        !Multiply residual by M^{-1}
        Residual_GB(1:nI,1:nJ,1:nK,iBlock) = &
             DivBV_G(1:nI,1:nJ,1:nK)*Prec_CB(:,:,:,iBlock)

        !DivBAbsMax=max(DivBAbsMax,vInv_CB(iBlock)*abs(&
        !               DivBV_G(1:nI,1:nJ,1:nK))))

        !Integrate divB by the whole volume
        DivBInt(ResDotOne_)=DivBInt(ResDotOne_)&
             +sum(DivBV_G(1:nI,1:nJ,1:nK))

        !Calculate Resudual.M.Residual
        DivBInt(ResDotM1DotRes_)=DivBInt(ResDotM1DotRes_)&
             +sum(DivBV_G(1:nI,1:nJ,1:nK)*Residual_GB(1:nI,1:nJ,1:nK,iBlock))
     end do

     if(nProc>1)then
        call MPI_allreduce(DivBInt,DivBTemp, 2,  MPI_REAL, MPI_SUM, &
             iComm, iError)
     else
        DivBTemp=DivBInt
     end if
     !Now evaluate the norm only for that part of the residual,
     !which is orthogobal to 1: 
     !(Res.M^{-1}.Res)_{mod}+(Res.M^{-1}.(1/sqrt(1.M^{-1}.1)))^2=Res.M^{-1}.Res
     DivBInt(ResDotM1DotRes_)=DivBTemp(ResDotM1DotRes_)&
          -OneDotMDotOneInv*DivBTemp(ResDotOne_)**2
     if(iProc==0.and.oktest)&
          write(*,*) ' Iteration=',Iteration,& 
          'SigmaRes,SigmaRes^2,integral error:',&
          DivBTemp,DivBInt(ResDotM1DotRes_)
     if(DivBInt(ResDotM1DotRes_)<Tolerance)EXIT

     !This is what we should add to the Residual to make it ortogonal to the 
     !constant solution
     DivBInt(ResDotOne_)=-OneDotMDotOneInv*DivBTemp(ResDotOne_)

     !Below we divide per Res.M.Res, so now we get inverse of it
     DivBInt(ResDotM1DotRes_)=1.0/DivBInt(ResDotM1DotRes_)
     do iBlock=1,nBlock
        if(Unused_B(iBlock))CYCLE
        if(true_blk(iBlock))then
           Dir_GB(1:nI,1:nJ,1:nK,iBlock)=Dir_GB(1:nI,1:nJ,1:nK,iBlock)+&
                DivBInt(ResDotM1DotRes_)*(Residual_GB(1:nI,1:nJ,1:nK,iBlock)+&
                DivBInt(ResDotOne_))
        else
           do k=1,nK;do j=1,nJ;do i=1,nI
              if(.not.true_cell(i,j,k,iBlock))CYCLE
              tmp2_blk(i,j,k,iBlock)=tmp2_blk(i,j,k,iBlock)+&
                   DivBInt(ResDotM1DotRes_)*(tmp1_blk(i,j,k,iBlock)+&
                   DivBInt(ResDotOne_))
           end do;end do;end do
        end if
     end do

     call message_pass_cell(tmp2_blk, nWidthIn=1, DoSendCornerIn=.false. ,&
          nProlongOrderIn=1, DoRestrictFaceIn=.true.)

     !Calculate Dir.A.Dir
     DirDotDir = 0.0
     do iBlock=1,nBlock
        if(Unused_B(iBlock))CYCLE
        call v_grad_phi(tmp2_blk,iBlock)
        DirDotDir = DirDotDir + &
             sum( (vDotGradX_C**2+vDotGradY_C**2+vDotGradZ_C**2) &
             /CellVolume_GB(1:nI,1:nJ,1:nK,iBlock) )
     end do
     if(nProc>1)then
        call MPI_allreduce(DirDotDir,DirDotDirInv, 1,  MPI_REAL, MPI_SUM, &
             iComm, iError)
        DirDotDirInv = 1.0/DirDotDirInv
     else
        DirDotDirInv = 1.0/DirDotDir
     end if

     if(oktest .and. iProc==0) write(*,*)'Effective diffusion coefficient = ',&
          DirDotDirInv*DivBInt(ResDotM1DotRes_)
     ! If iterations are not used, the diffusion coefficient 
     ! should be less than 1, to ensure the convergence in the r.M^{-1}.r norm
     if(nCleanDivb==0)DirDotDirInv = &
          min(DirDotDirInv,0.99/DivBInt(ResDotM1DotRes_))

     do iBlock=1,nBlock
        if (Unused_B(iBlock)) CYCLE
        call v_grad_phi(Dir_GB,iBlock)

        State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock) = &
             State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock)-&
             vDotGradX_C*DirDotDirInv/CellVolume_GB(1:nI,1:nJ,1:nK,iBlock)
        State_VGB(By_,1:nI,1:nJ,1:nK,iBlock) = &
             State_VGB(By_,1:nI,1:nJ,1:nK,iBlock)-&
             vDotGradY_C*DirDotDirInv/CellVolume_GB(1:nI,1:nJ,1:nK,iBlock)
        State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)=&
             State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)-&
             vDotGradZ_C*DirDotDirInv/CellVolume_GB(1:nI,1:nJ,1:nK,iBlock)
        if(DoConservative)&
             State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)=& 
             State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)+&
             0.5*gm1*DirDotDirInv*DirDotDirInv*&
             (vDotGradX_C**2+vDotGradY_C**2+vDotGradZ_C**2) &
             /CellVolume_GB(1:nI,1:nJ,1:nK,iBlock)**2
     end do
     Iteration=Iteration+1   
     if(Iteration>nCleanDivb)EXIT    
  end do

  call timing_stop('clean_divb')

contains
  !============================================================================
  subroutine init_divb_cleanup
    !    use ModAdvance,ONLY: DivB1_GB
    implicit none

    integer::i,j,k,iBlock,iError,iLimit
    real,dimension(0:nI+1,0:nJ+1,0:nK+1)::Q_G
    real::EstimateForMAMNorm,OneDotMDotOne 
    real:: divb_diffcoeff, VInvHalf
    !--------------------------------------------------------------------------
    tmp1_blk = 0.0; Prec_CB = 0.0

    do iBlock=1,nBlock
       if(Unused_B(iBlock))CYCLE
       tmp1_blk(1:nI,1:nJ,1:nK,iBlock) = 1/CellVolume_GB(1:nI,1:nJ,1:nK,iBlock)
    end do

    !tmp1 is equal to the inverse of the volume, including the ghostcells 
    call message_pass_cell(tmp1_blk,DoSendCornerIn=.false. ,&
         nProlongOrderIn=1, DoRestrictFaceIn=.true.)

    do iBlock=1,nBlock
       if (Unused_B(iBlock)) CYCLE

       Q_G = 1.0

       if(any(NeiLev(:,iBlock)/=0))then
          if(NeiLev(1,iBlock)==NoBLK)then
             tmp1_blk(0,1:nJ,1:nK,iBlock) = &
                  1.0/CellVolume_GB(1,1:nJ,1:nK,iBlock)
             !to define somehow tmp1 in the outer ghostcells

          elseif(abs(NeiLev(1,iBlock))==1)then           
             Q_G(0,:,:)=4.0**NeiLev(1,iBlock)
             !if the neighboring block is coarser,
             !  the input from FA^2 should be multipled by four
             !If the neighborig block is finer,
             !  the input from FA^2  should be multiplied by 1/4
          end if
          if(NeiLev(2,iBlock)==NoBLK)then
             tmp1_blk(nI+1,1:nJ,1:nK,iBlock) = &
                  1.0/CellVolume_GB(nI,1:nJ,1:nK,iBlock)
          elseif(abs(NeiLev(2,iBlock))==1)then
             Q_G(nI+1,:,:)=4.0**NeiLev(2,iBlock)
          end if
          if(NeiLev(3,iBlock)==NoBLK)then
             tmp1_blk(1:nI,0,1:nK,iBlock) = &
                  1.0/CellVolume_GB(1:nI,1,1:nK,iBlock)
          elseif(abs(NeiLev(3,iBlock))==1)then
             Q_G(:,0,:)=4.0**NeiLev(3,iBlock)
          end if
          if(NeiLev(4,iBlock)==NoBLK)then
             tmp1_blk(1:nI,nJ+1,1:nK,iBlock) = &
                  1.0/CellVolume_GB(1:nI,nJ,1:nK,iBlock)
          elseif(abs(NeiLev(4,iBlock))==1)then
             Q_G(:,nJ+1,:)=4.0**NeiLev(4,iBlock)
          end if
          if(NeiLev(5,iBlock)==NoBLK)then
             tmp1_blk(1:nI,1:nJ,0,iBlock) = &
                  1.0/CellVolume_GB(1:nI,1:nJ,1,iBlock)
          elseif(abs(NeiLev(5,iBlock))==1)then
             Q_G(:,:,0)=4.0**NeiLev(5,iBlock)
          end if
          if(NeiLev(6,iBlock)==NoBLK)then
             tmp1_blk(1:nI,1:nJ,nK+1,iBlock) = &
                  1.0/CellVolume_GB(1:nI,1:nJ,nK,iBlock)
          elseif(abs(NeiLev(6,iBlock))==1)then
             Q_G(:,:,nK+1)=4.0**NeiLev(6,iBlock)
          end if
       end if
       Prec_CB(:,:,:,iBlock)=4.0/(&
            (Q_G(2:nI+1,1:nJ,1:nK)*tmp1_blk(2:nI+1,1:nJ,1:nK,iBlock)+& 
            Q_G(0:nI-1,1:nJ,1:nK)*tmp1_blk(0:nI-1,1:nJ,1:nK,iBlock))& 
            *CellFace_DB(x_,iBlock)**2+&
            (Q_G(1:nI,2:nJ+1,1:nK)*tmp1_blk(1:nI,2:nJ+1,1:nK,iBlock)+&
            Q_G(1:nI,0:nJ-1,1:nK)*tmp1_blk(1:nI,0:nJ-1,1:nK,iBlock))&
            *CellFace_DB(y_,iBlock)**2+&
            (Q_G(1:nI,1:nJ,2:nK+1)*tmp1_blk(1:nI,1:nJ,2:nK+1,iBlock)+&
            Q_G(1:nI,1:nJ,0:nK-1)*tmp1_blk(1:nI,1:nJ,0:nK-1,iBlock))&
            *CellFace_DB(z_,iBlock)**2)
    end do
    do iLimit=1,2
       do iBlock=1,nBlock
          if(Unused_B(iBlock))CYCLE
          if(any(NeiLev(:,iBlock)==NoBLK))then
             if(NeiLev(1,iBlock)==NoBLK)&
                  tmp1_blk(0,1:nJ,1:nK,iBlock)=Prec_CB(1,:,:,iBlock)
             if(NeiLev(2,iBlock)==NoBLK)&
                  tmp1_blk(nI+1,1:nJ,1:nK,iBlock)=Prec_CB(nI,:,:,iBlock)
             if(NeiLev(3,iBlock)==NoBLK)&
                  tmp1_blk(1:nI,0,1:nK,iBlock)=Prec_CB(:,1,:,iBlock)
             if(NeiLev(4,iBlock)==NoBLK)&
                  tmp1_blk(1:nI,nJ+1,1:nK,iBlock)=Prec_CB(:,nJ,:,iBlock)
             if(NeiLev(5,iBlock)==NoBLK)&
                  tmp1_blk(1:nI,1:nJ,0,iBlock)=Prec_CB(:,:,1,iBlock)
             if(NeiLev(6,iBlock)==NoBLK)&
                  tmp1_blk(1:nI,1:nJ,nK+1,iBlock)=Prec_CB(:,:,nK,iBlock)
          end if
          tmp1_blk(1:nI,1:nJ,1:nK,iBlock)=Prec_CB(:,:,:,iBlock)
       end do

       call message_pass_cell(tmp1_blk,DoSendCornerIn=.false. ,&
            nProlongOrderIn=1, DoRestrictFaceIn=.true.)

       do iBlock=1,nBlock
          if (Unused_B(iBlock)) CYCLE

          do k=1,nK;do j=1,nJ;do i=1,nI
             Prec_CB(i,j,k,iBlock)=min(&
                  minval(tmp1_blk(i-1:i+1,j,k,iBlock)),&
                  minval(tmp1_blk(i,j-1:j+1:2,k,iBlock)),& 
                  minval(tmp1_blk(i,j,k-1:k+1:2,iBlock)))
          end do;end do; end do
       end do
    end do

    !With the preconditioner, defined in the manner like this
    !the upper estimate of the spectral radius fo the matrix
    !||Prec^{1/2}.A.Prec^1/2 is 2, according to the Frobenius inequality 

    !Herewith Prec is also denoted as M^{-1}
    !Calculate the divB diffusion coefficient =2/||M^(-1/2).A.M^(-1/2)||
    !For the grid, which is at least partially uniform cartesian,
    !this coefficient equals 1
    do iBlock=1,nBlock
       if (Unused_B(iBlock)) CYCLE
       if(any(NeiLev(:,iBlock)==NoBLK))then
          if(NeiLev(1,iBlock)==NoBLK)&
               tmp1_blk(0,1:nJ,1:nK,iBlock)=sqrt(Prec_CB(1,:,:,iBlock))
          if(NeiLev(2,iBlock)==NoBLK)&
               tmp1_blk(nI+1,1:nJ,1:nK,iBlock)=sqrt(Prec_CB(nI,:,:,iBlock))
          if(NeiLev(3,iBlock)==NoBLK)&
               tmp1_blk(1:nI,0,1:nK,iBlock)=sqrt(Prec_CB(:,1,:,iBlock))
          if(NeiLev(4,iBlock)==NoBLK)&
               tmp1_blk(1:nI,nJ+1,1:nK,iBlock)=sqrt(Prec_CB(:,nJ,:,iBlock))
          if(NeiLev(5,iBlock)==NoBLK)&
               tmp1_blk(1:nI,1:nJ,0,iBlock)=sqrt(Prec_CB(:,:,1,iBlock))
          if(NeiLev(6,iBlock)==NoBLK)&
               tmp1_blk(1:nI,1:nJ,nK+1,iBlock)=sqrt(Prec_CB(:,:,nK,iBlock))
       end if
       tmp1_blk(1:nI,1:nJ,1:nK,iBlock)=sqrt(Prec_CB(:,:,:,iBlock))
    end do
    !    if(iProc==0)write(*,*)' Cleanup Initialization second message pass'

    call message_pass_cell(tmp1_blk,DoSendCornerIn=.false. ,&
         nProlongOrderIn=1, DoRestrictFaceIn=.true.)

    !Now the elements of diag(Prec_CB)^{1/2} are in tmp1_blk

    do iBlock=1,nBlock
       if (Unused_B(iBlock)) CYCLE

       Q_G=tmp1_blk(0:nI+1,0:nJ+1,0:nK+1,iBlock)
       tmp1_BLK(:,:,:,iBlock) = 0.0
       tmp2_BLK(:,:,:,iBlock) = 0.0
       tmp3_BLK(:,:,:,iBlock) = 0.0

       do k=1,nK;do j=1,nJ;do i=1,nI
          VInvHalf = 0.5/CellVolume_GB(i,j,k,iBlock)
          tmp1_BLK(i,j,k,iBlock) = &
               CellFace_DB(x_,iBlock)*VInvHalf*(&
               Q_G(i+1,j,k)+&
               Q_G(i-1,j,k))
          tmp2_BLK(i,j,k,iBlock) = &
               CellFace_DB(y_,iBlock)*VInvHalf*(&
               Q_G(i,j+1,k)+&
               Q_G(i,j-1,k))
          tmp3_BLK(i,j,k,iBlock) = &
               CellFace_DB(z_,iBlock)*VInvHalf*(&
               Q_G(i,j,k+1)+&
               Q_G(i,j,k-1))
       end do;end do; end do
    end do
    !In tmp1,tmp2 and divB1 are the estimates of gradX, gradY,gradZ correspondenyly

    call message_pass_cell(tmp1_blk,DoSendCornerIn=.false. ,&
         nProlongOrderIn=1, DoRestrictFaceIn=.true.)
    call message_pass_cell(tmp2_blk,DoSendCornerIn=.false. ,&
         nProlongOrderIn=1, DoRestrictFaceIn=.true.)
    call message_pass_cell(tmp3_blk,DoSendCornerIn=.false. ,&
         nProlongOrderIn=1, DoRestrictFaceIn=.true.)

    EstimateForMAMNorm = 0.0
    do iBlock=1,nBlock
       if(Unused_B(iBlock))CYCLE
       EstimateForMAMNorm=max(EstimateForMAMNorm,&
            maxval(sqrt(Prec_CB(:,:,:,iBlock))*&
            0.5*( &
            CellFace_DB(x_,iBlock)*&
            (tmp1_blk(2:nI+1, 1:nJ, 1:nK,iBlock)+&
            tmp1_blk(0:nI-1, 1:nJ, 1:nK,iBlock))&
            +CellFace_DB(y_,iBlock)*&
            (tmp2_blk(1:nI, 2:nJ+1, 1:nK,iBlock)+&
            tmp2_blk(1:nI, 0:nJ-1, 1:nK,iBlock))&
            +CellFace_DB(z_,iBlock)*&
            (tmp3_blk(1:nI, 1:nJ, 2:nK+1,iBlock)+&
            tmp3_blk(1:nI, 1:nJ, 0:nK-1,iBlock)) ) ))
    end do

    if(nProc>1)then
       call MPI_allreduce(EstimateForMAMNorm,divb_diffcoeff,  &
            1, MPI_REAL, MPI_MAX, iComm, iError)
       divb_diffcoeff = 2/divb_diffcoeff
    else
       divb_diffcoeff = 2/EstimateForMAMNorm
    end if
    if(iProc==0)write(*,*)"Divb diffusion coefficient is: ",divb_diffcoeff

    !Compute 1/sum(M_i)
    OneDotMDotOne = 0.0
    do iBlock=1,nBlock
       if (Unused_B(iBlock)) CYCLE
       !     Prec_CB(:,:,:,iBlock)=Prec_CB(:,:,:,iBlock)*divb_diffcoeff
       if(true_blk(iBlock))then
          OneDotMDotOne=OneDotMDotOne+sum(1.0/Prec_CB(:,:,:,iBlock))
       elseif(any(true_cell(1:nI,1:nJ,1:nK,iBlock)))then
          OneDotMDotOne=OneDotMDotOne+sum(1.0/Prec_CB(:,:,:,iBlock)&
               ,MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))
       end if
    end do
    if(nProc>1)then
       call MPI_allreduce(OneDotMDotOne,OneDotMDotOneInv, 1,  MPI_REAL, MPI_SUM, &
            iComm, iError)
       OneDotMDotOneInv = 1.0/OneDotMDotOneInv
    else
       OneDotMDotOneInv = 1.0/OneDotMDotOne
    end if
    if(iProc==0)write(*,*)' init_divb_cleanup finishes with OneDotMDotOneInv=',OneDotMDotOneInv
    !    write(*,*)'Maxval loc and minval loc of Prec_CB=', &
    !         maxval(Prec_CB(:,:,:,1:nBlock)), maxloc(Prec_CB(:,:,:,1:nBlock)), &
    !         minval(Prec_CB(:,:,:,1:nBlock)), minloc(Prec_CB(:,:,:,1:nBlock)) 
  end subroutine init_divb_cleanup
  !===========================================================================
  subroutine v_grad_phi(Phi_GB,iBlock)

    integer,intent(in)::iBlock
    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK),intent(inout)::Phi_GB
    !------------------------------------------------------------------------

    vDotGradX_C = 0.0;vDotGradY_C = 0.0;vDotGradZ_C = 0.0
!!! Apply continuous solution at east and west
    !    if (NeiLev(1,iBlock)==NOBLK)&
    !         Phi_GB(0   ,1:nJ,1:nK,iBlock) = Phi_GB(1 ,1:nJ,1:nK,iBlock)
    !    if (NeiLev(2,iBlock)==NOBLK)&
    !         Phi_GB(nI+1,1:nJ,1:nK,iBlock) = Phi_GB(nI,1:nJ,1:nK,iBlock)
!!! Apply shearing at north and south
    !    if (NeiLev(3,iBlock)==NOBLK)&
    !         Phi_GB(1:nI,0   ,1:nK,iBlock) = Phi_GB(0:nI-1,2 ,1:nK,iBlock)
    !    if (NeiLev(4,iBlock)==NOBLK)&
    !         Phi_GB(1:nI,nJ+1,1:nK,iBlock) = Phi_GB(2:nI+1,nJ-1,1:nK,iBlock)
!!! Apply translation invariant solution at bottom and top
    !    if (NeiLev(5,iBlock)==NOBLK)&
    !         Phi_GB(1:nI,1:nJ,0   ,iBlock) = Phi_GB(1:nI,1:nJ,1 ,iBlock)
    !    if (NeiLev(6,iBlock)==NOBLK)&
    !         Phi_GB(1:nI,1:nJ,nK+1,iBlock) = Phi_GB(1:nI,1:nJ,nK,iBlock)

    if (NeiLev(1,iBlock) == NOBLK) Phi_GB(0,1:nJ,1:nK,iBlock) = &
         -BoundaryCoef*Phi_GB(1 ,1:nJ,1:nK,iBlock)
    if (NeiLev(2,iBlock) == NOBLK) Phi_GB(nI+1,1:nJ,1:nK,iBlock) = &
         -BoundaryCoef*Phi_GB(nI,1:nJ,1:nK,iBlock)
    if (NeiLev(3,iBlock) == NOBLK) Phi_GB(1:nI,j0_ ,1:nK,iBlock) = &
         -BoundaryCoef*Phi_GB(1:nI,1 ,1:nK,iBlock)
    if (NeiLev(4,iBlock) == NOBLK) Phi_GB(1:nI,nJp1_,1:nK,iBlock) = &
         -BoundaryCoef*Phi_GB(1:nI,nJ,1:nK,iBlock)
    if (NeiLev(5,iBlock) == NOBLK .and. nK>1)Phi_GB(1:nI,1:nJ,k0_,iBlock) = &
         -BoundaryCoef*Phi_GB(1:nI,1:nJ,1 ,iBlock)
    if (NeiLev(6,iBlock) == NOBLK .and. nK>1)Phi_GB(1:nI,1:nJ,nKp1_,iBlock) = &
         -BoundaryCoef*Phi_GB(1:nI,1:nJ,nK,iBlock)

    if(body_blk(iBlock))then
       do k = 1, nK; do j = 1, nJ;do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock))CYCLE
          where(.not.true_cell(i-1:i+1:2,j,k,iBlock))&
               Phi_GB(i-1:i+1:2,j,k,iBlock)=-BoundaryCoef*Phi_GB(i,j,k,iBlock)
          where(.not.true_cell(i,j-1:j+1:2,k,iBlock))&
               Phi_GB(i,j-1:j+1:2,k,iBlock)=-BoundaryCoef*Phi_GB(i,j,k,iBlock)
          where(.not.true_cell(i,j,k-1:k+1:2,iBlock))&
               Phi_GB(i,j,k-1:k+1:2,iBlock)=-BoundaryCoef*Phi_GB(i,j,k,iBlock)
          vDotGradX_C(i,j,k)=&
               0.5*CellFace_DB(x_,iBlock)*(&
               Phi_GB(i+1,j,k,iBlock)-&
               Phi_GB(i-1,j,k,iBlock))
          vDotGradY_C(i,j,k)=&
               0.5*CellFace_DB(y_,iBlock)*(&
               Phi_GB(i,j+1,k,iBlock)-&
               Phi_GB(i,j-1,k,iBlock))
          vDotGradZ_C(i,j,k)=&
               0.5*CellFace_DB(z_,iBlock)*(&
               Phi_GB(i,j,k+1,iBlock)-&
               Phi_GB(i,j,k-1,iBlock))
       end do;end do;end do
    else
       do k=1,nK;do j=1,nJ;do i=1,nI
          vDotGradX_C(i,j,k)=&
               0.5*CellFace_DB(x_,iBlock)*(&
               Phi_GB(i+1,j,k,iBlock)-&
               Phi_GB(i-1,j,k,iBlock))
          vDotGradY_C(i,j,k)=&
               0.5*CellFace_DB(y_,iBlock)*(&
               Phi_GB(i,j+1,k,iBlock)-&
               Phi_GB(i,j-1,k,iBlock))
          vDotGradZ_C(i,j,k)=&
               0.5*CellFace_DB(z_,iBlock)*(&
               Phi_GB(i,j,k+1,iBlock)-&
               Phi_GB(i,j,k-1,iBlock))
       end do;end do;end do
    end if
  end subroutine v_grad_phi
end subroutine clean_divb
!===================================================================
subroutine div_3d_b1(iBlock,VecX_G,VecY_G,VecZ_G,Out_G)     
  use ModSize
  use ModGeometry, ONLY: body_blk, true_cell
  use ModParallel, ONLY: neilev, NOBLK
  use ModDivbCleanup, ONLY: BoundaryCoef
  use BATL_lib, ONLY: CellFace_DB
  implicit none

  integer,intent(in) :: iBlock
  real,dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK),intent(in)::&
       VecX_G,VecY_G,VecZ_G
  real,dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK),intent(out)::Out_G

  real, dimension(0:nI+1, 0:nJ+1, 0:nK+1) :: OneTrue_G

  integer :: i, j, k

  !\
  ! Can only be used for divB diffusion and projection scheme!!!!
  ! DivB is multiplied by -V_Cell!!!!
  ! With this modification DivB[grad Phi] is a symmetric positive definite 
  ! operator!
  !/
  Out_G = 0.0
  if(.not.(body_blk(iBlock).or.any(neilev(:,iBlock)==NOBLK)))then 
     do k=1,nK; do j=1,nJ; do i=1,nI
        Out_G(i,j,k) = - 0.5*(&
             CellFace_DB(x_,iBlock)*&
             (VecX_G(i+1, j, k)-VecX_G(i-1,j,k))&
             +CellFace_DB(y_,iBlock)*&
             (VecY_G(i ,j+1, k)-VecY_G(i,j-1,k))&
             +CellFace_DB(z_,iBlock)*&
             (VecZ_G(i, j, k+1)-VecZ_G(i,j,k-1)) )
     end do; end do; end do
  else
     where(true_cell(0:nI+1, 0:nJ+1, 0:nK+1,iBlock)) 
        OneTrue_G = 1.0
     elsewhere
        OneTrue_G = 0.0
     end where
     if(neilev(1 ,iBlock)==NOBLK) OneTrue_G(0   ,:,:) = 0.0
     if(neilev(2 ,iBlock)==NOBLK) OneTrue_G(nI+1,:,:) = 0.0
     if(neilev(3,iBlock)==NOBLK) OneTrue_G(:,0   ,:) = 0.0
     if(neilev(4,iBlock)==NOBLK) OneTrue_G(:,nJ+1,:) = 0.0
     if(neilev(5  ,iBlock)==NOBLK) OneTrue_G(:,:,0   ) = 0.0
     if(neilev(6  ,iBlock)==NOBLK) OneTrue_G(:,:,nK+1) = 0.0
     !
     !\
     ! Where .not.true_cell, all the gradients are zero
     ! In true_cell the input to gradient from the face neighbor
     ! is ignored, if the face neighbor is .not.true_cell
     !/
     !
     do k=1,nK; do j=1,nJ; do i=1,nI
        Out_G(i,j,k) = - 0.5*OneTrue_G(i,j,k)*(+&
             CellFace_DB(x_,iBlock)*&
             (VecX_G(i+1,j,k)-&
             BoundaryCoef*(VecX_G(i+1,j,k)-VecX_G(i,j,k))*&
             (1.0-OneTrue_G(i+1,j,k))-&
             VecX_G(i-1,j,k)-&
             BoundaryCoef*(VecX_G(i,j,k)-VecX_G(i-1,j,k))*&
             (1.0-OneTrue_G(i-1,j,k)))+&
             CellFace_DB(y_,iBlock)*&
             (VecY_G(i,j+1,k)-&
             BoundaryCoef*(VecY_G(i,j+1,k)-VecY_G(i,j,k))*&
             (1.0-OneTrue_G(i,j+1,k))-&
             VecY_G(i,j-1,k)-&
             BoundaryCoef*(VecY_G(i,j,k)-VecY_G(i,j-1,k))*&
             (1.0-OneTrue_G(i,j-1,k)))+&
             CellFace_DB(z_,iBlock)*&
             (VecZ_G(i,j,k+1)-&
             BoundaryCoef*(VecZ_G(i,j,k+1)-VecZ_G(i,j,k))*&
             (1.0-OneTrue_G(i,j,k+1))-&
             VecZ_G(i,j,k-1)-&
             BoundaryCoef*(VecZ_G(i,j,k)-VecZ_G(i,j,k-1))*&
             (1.0-OneTrue_G(i,j,k-1))) &
             )
     end do; end do; end do
  end if

end subroutine div_3d_b1
