!^CFG COPYRIGHT UM
!^CFG FILE DIVBDIFFUSE
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
  use ModNumConst
  use ModDivbCleanup
  use ModMain,ONLY: iNewGrid, iNewDecomposition, nBlock, unusedblk
  use ModAdvance,ONLY:State_VGB, Bx_, By_, Bz_, P_,tmp1_BLK,tmp2_BLK,&
       Residual_GB=>tmp1_blk,Dir_GB=>tmp2_blk
  use ModAdvance,ONLY:tmp3_blk=>divB1_GB
  use ModGeometry,ONLY:vInv_CB,true_blk,&
       FaX_BLK,FaY_BLK,FaZ_BLK,body_blk,true_cell
  use ModParallel, ONLY : NOBLK, neiLEV
  use ModPhysics,ONLY:gm1
  use ModMpi
  implicit none

  integer::i,j,k,iBlock
  logical::DoConservative=.false.
  real,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn)::DivBV_G
  real:: DivBInt(2),DivBTemp(2)
  integer,parameter::ResDotM1DotRes_=2,ResDotOne_=1

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
     call message_pass_cells_8state(.true.,.true.,.true.)    
     !Get the ghostcell values for MF 
     DivBInt=cZero;DivBTemp=cZero
     do iBlock=1,nBlock
        Residual_GB(:,:,:,iBlock)=cZero
        if(Iteration==1) Dir_GB(:,:,:,iBlock)=cZero 
        !Initialize the vector "Dir"
        
        if(unusedBLK(iBlock))CYCLE
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
     DivBInt(ResDotM1DotRes_)=cOne/DivBInt(ResDotM1DotRes_)
     do iBlock=1,nBlock
        if(unusedBLK(iBlock))CYCLE
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
     call message_pass_cells(.true.,.true.,.true.,tmp2_blk)
     !Calculate Dir.A.Dir
     DirDotDir=cZero
     do iBlock=1,nBlock
        if(unusedBLK(iBlock))CYCLE
        call v_grad_phi(tmp2_blk,iBlock)
        DirDotDir=DirDotDir+&
             sum(vInv_CB(:,:,:,iBlock)*&
             (vDotGradX_C**2+vDotGradY_C**2+vDotGradZ_C**2))
     end do
     if(nProc>1)then
        call MPI_allreduce(DirDotDir,DirDotDirInv, 1,  MPI_REAL, MPI_SUM, &
             iComm, iError)
        DirDotDirInv=cOne/DirDotDirInv
     else
        DirDotDirInv=cOne/DirDotDir
     end if

     if(oktest .and. iProc==0) write(*,*)'Effective diffusion coefficient = ',&
                DirDotDirInv*DivBInt(ResDotM1DotRes_)
     ! If iterations are not used, the diffusion coefficient 
     ! should be less than 1, to ensure the convergence in the r.M^{-1}.r norm
     if(nCleanDivb==0)DirDotDirInv = &
          min(DirDotDirInv,0.99/DivBInt(ResDotM1DotRes_))

     do iBlock=1,nBlock
        if (unusedBLK(iBlock)) CYCLE
        call v_grad_phi(Dir_GB,iBlock)
       
        State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock) = &
             State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock)-&
             vDotGradX_C*vInv_CB(:,:,:,iBlock)*DirDotDirInv
        State_VGB(By_,1:nI,1:nJ,1:nK,iBlock) = &
             State_VGB(By_,1:nI,1:nJ,1:nK,iBlock)-&
             vDotGradY_C*vInv_CB(:,:,:,iBlock)*DirDotDirInv
        State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)=&
             State_VGB(Bz_,1:nI,1:nJ,1:nK,iBlock)-&
             vDotGradZ_C*vInv_CB(:,:,:,iBlock)*DirDotDirInv
        if(DoConservative)&
             State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)=& 
             State_VGB(P_,1:nI,1:nJ,1:nK,iBlock)+&
             cHalf*gm1*DirDotDirInv*DirDotDirInv*&
             (vDotGradX_C**2+vDotGradY_C**2+vDotGradZ_C**2)*&
             vInv_CB(:,:,:,iBlock)**2
     end do
     Iteration=Iteration+1   
     if(Iteration>nCleanDivb)EXIT    
  end do

  call timing_stop('clean_divb')

contains
  !=============================================================================
  subroutine init_divb_cleanup
!    use ModAdvance,ONLY: DivB1_GB
    implicit none

    integer::i,j,k,iBlock,iError,iLimit
    real,dimension(0:nI+1,0:nJ+1,0:nK+1)::Q_G
    real::EstimateForMAMNorm,OneDotMDotOne 
    real:: divb_diffcoeff
    tmp1_blk=cZero;Prec_CB=cZero
    do iBlock=1,nBlock
       if(unusedBLK(iBlock))CYCLE
       tmp1_blk(1:nI,1:nJ,1:nK,iBlock)=vInv_CB(:,:,:,iBlock)
    end do



    call message_pass_cells(.false.,.true.,.true.,tmp1_blk)
!tmp1 is equal to the inverse of the volume, including the ghostcells 

    do iBlock=1,nBlock
       if (unusedBLK(iBlock)) CYCLE

       Q_G=cOne

       if(any(NeiLev(:,iBlock)/=0))then
          if(NeiLev(East_,iBlock)==NoBLK)then
             tmp1_blk(0,1:nJ,1:nK,iBlock)=vInv_CB(1,:,:,iBlock) 
             !to define somehow tmp1 in the outer ghostcells

          elseif(abs(NeiLev(East_,iBlock))==1)then           
             Q_G(0,:,:)=4.0**NeiLev(East_,iBlock)
             !if the neighboring block is coarser, the input from FA^2 should be multipled by four
             !If the neighborig block is finer, the input from FA^2  should be multiplied by 1/4
          end if
          if(NeiLev(West_,iBlock)==NoBLK)then
             tmp1_blk(nI+1,1:nJ,1:nK,iBlock)=vInv_CB(nI,:,:,iBlock)
          elseif(abs(NeiLev(West_,iBlock))==1)then
             Q_G(nI+1,:,:)=4.0**NeiLev(West_,iBlock)
          end if
          if(NeiLev(South_,iBlock)==NoBLK)then
             tmp1_blk(1:nI,0,1:nK,iBlock)=vInv_CB(:,1,:,iBlock)
          elseif(abs(NeiLev(South_,iBlock))==1)then
             Q_G(:,0,:)=4.0**NeiLev(South_,iBlock)
          end if
          if(NeiLev(North_,iBlock)==NoBLK)then
             tmp1_blk(1:nI,nJ+1,1:nK,iBlock)=vInv_CB(:,nJ,:,iBlock)
          elseif(abs(NeiLev(North_,iBlock))==1)then
             Q_G(:,nJ+1,:)=4.0**NeiLev(North_,iBlock)
          end if
          if(NeiLev(Bot_,iBlock)==NoBLK)then
             tmp1_blk(1:nI,1:nJ,0,iBlock)=vInv_CB(:,:,1,iBlock)
          elseif(abs(NeiLev(Bot_,iBlock))==1)then
             Q_G(:,:,0)=4.0**NeiLev(Bot_,iBlock)
          end if
          if(NeiLev(Top_,iBlock)==NoBLK)then
             tmp1_blk(1:nI,1:nJ,nK+1,iBlock)=vInv_CB(:,:,nK,iBlock)
          elseif(abs(NeiLev(Top_,iBlock))==1)then
             Q_G(:,:,nK+1)=4.0**NeiLev(Top_,iBlock)
          end if
       end if
       Prec_CB(:,:,:,iBlock)=4.0/(&
            (Q_G(2:nI+1,1:nJ,1:nK)*tmp1_blk(2:nI+1,1:nJ,1:nK,iBlock)+& 
            Q_G(0:nI-1,1:nJ,1:nK)*tmp1_blk(0:nI-1,1:nJ,1:nK,iBlock))& 
            *FaX_BLK(iBlock)**2+&
            (Q_G(1:nI,2:nJ+1,1:nK)*tmp1_blk(1:nI,2:nJ+1,1:nK,iBlock)+&
            Q_G(1:nI,0:nJ-1,1:nK)*tmp1_blk(1:nI,0:nJ-1,1:nK,iBlock))&
            *FaY_BLK(iBlock)**2+&
            (Q_G(1:nI,1:nJ,2:nK+1)*tmp1_blk(1:nI,1:nJ,2:nK+1,iBlock)+&
            Q_G(1:nI,1:nJ,0:nK-1)*tmp1_blk(1:nI,1:nJ,0:nK-1,iBlock))&
            *FaZ_BLK(iBlock)**2)
    end do
    do iLimit=1,2
       do iBlock=1,nBlock
          if(unusedBLK(iBlock))CYCLE
          if(any(NeiLev(:,iBlock)==NoBLK))then
             if(NeiLev(East_,iBlock)==NoBLK)&
                  tmp1_blk(0,1:nJ,1:nK,iBlock)=Prec_CB(1,:,:,iBlock)
             if(NeiLev(West_,iBlock)==NoBLK)&
                  tmp1_blk(nI+1,1:nJ,1:nK,iBlock)=Prec_CB(nI,:,:,iBlock)
             if(NeiLev(South_,iBlock)==NoBLK)&
                  tmp1_blk(1:nI,0,1:nK,iBlock)=Prec_CB(:,1,:,iBlock)
             if(NeiLev(North_,iBlock)==NoBLK)&
                  tmp1_blk(1:nI,nJ+1,1:nK,iBlock)=Prec_CB(:,nJ,:,iBlock)
             if(NeiLev(Bot_,iBlock)==NoBLK)&
                  tmp1_blk(1:nI,1:nJ,0,iBlock)=Prec_CB(:,:,1,iBlock)
             if(NeiLev(Top_,iBlock)==NoBLK)&
                  tmp1_blk(1:nI,1:nJ,nK+1,iBlock)=Prec_CB(:,:,nK,iBlock)
          end if
          tmp1_blk(1:nI,1:nJ,1:nK,iBlock)=Prec_CB(:,:,:,iBlock)
       end do
       call message_pass_cells(.false.,.true.,.true.,tmp1_blk)
       do iBlock=1,nBlock
          if (unusedBLK(iBlock)) CYCLE

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
       if (unusedBLK(iBlock)) CYCLE
       if(any(NeiLev(:,iBlock)==NoBLK))then
          if(NeiLev(East_,iBlock)==NoBLK)&
               tmp1_blk(0,1:nJ,1:nK,iBlock)=sqrt(Prec_CB(1,:,:,iBlock))
          if(NeiLev(West_,iBlock)==NoBLK)&
               tmp1_blk(nI+1,1:nJ,1:nK,iBlock)=sqrt(Prec_CB(nI,:,:,iBlock))
          if(NeiLev(South_,iBlock)==NoBLK)&
               tmp1_blk(1:nI,0,1:nK,iBlock)=sqrt(Prec_CB(:,1,:,iBlock))
          if(NeiLev(North_,iBlock)==NoBLK)&
               tmp1_blk(1:nI,nJ+1,1:nK,iBlock)=sqrt(Prec_CB(:,nJ,:,iBlock))
          if(NeiLev(Bot_,iBlock)==NoBLK)&
               tmp1_blk(1:nI,1:nJ,0,iBlock)=sqrt(Prec_CB(:,:,1,iBlock))
          if(NeiLev(Top_,iBlock)==NoBLK)&
               tmp1_blk(1:nI,1:nJ,nK+1,iBlock)=sqrt(Prec_CB(:,:,nK,iBlock))
       end if
       tmp1_blk(1:nI,1:nJ,1:nK,iBlock)=sqrt(Prec_CB(:,:,:,iBlock))
    end do
!    if(iProc==0)write(*,*)' Cleanup Initialization second message pass'
    call message_pass_cells(.false.,.true.,.true.,tmp1_blk)
    !Now the elements of diag(Prec_CB)^{1/2} are in tmp1_blk

    do iBlock=1,nBlock
       if (unusedBLK(iBlock)) CYCLE
       Q_G=tmp1_blk(0:nI+1,0:nJ+1,0:nK+1,iBlock)
       tmp1_BLK(:,:,:,iBlock)=cZero
       tmp2_BLK(:,:,:,iBlock)=cZero
       tmp3_BLK(:,:,:,iBlock)=cZero
       do k=1,nK;do j=1,nJ;do i=1,nI
          tmp1_BLK(i,j,k,iBlock)= &
               cHalf*FaX_BLK(iBlock)*vInv_CB(i,j,k,iBlock)*(&
               Q_G(i+1,j,k)+&
               Q_G(i-1,j,k))
          tmp2_BLK(i,j,k,iBlock)=&
               cHalf*FaY_BLK(iBlock)*vInv_CB(i,j,k,iBlock)*(&
               Q_G(i,j+1,k)+&
               Q_G(i,j-1,k))
          tmp3_BLK(i,j,k,iBlock)=&
               cHalf*FaZ_BLK(iBlock)*vInv_CB(i,j,k,iBlock)*(&
               Q_G(i,j,k+1)+&
               Q_G(i,j,k-1))
       end do;end do; end do
    end do
    !In tmp1,tmp2 and divB1 are the estimates of gradX, gradY,gradZ correspondenyly
    call message_pass_cells(.false.,.true.,.true.,tmp1_BLK)
    call message_pass_cells(.false.,.true.,.true.,tmp2_BLK)
    call message_pass_cells(.false.,.true.,.true.,tmp3_BLK)
    EstimateForMAMNorm=cZero
    do iBlock=1,nBlock
       if(unusedBLK(iBlock))CYCLE
       EstimateForMAMNorm=max(EstimateForMAMNorm,&
            maxval(sqrt(Prec_CB(:,:,:,iBlock))*&
            cHalf*(&
            fAX_BLK(iBlock)*&
            (tmp1_blk(2:nI+1, 1:nJ, 1:nK,iBlock)+&
            tmp1_blk(0:nI-1, 1:nJ, 1:nK,iBlock))&
            +fAY_BLK(iBlock)*&
            (tmp2_blk(1:nI, 2:nJ+1, 1:nK,iBlock)+&
            tmp2_blk(1:nI, 0:nJ-1, 1:nK,iBlock))&
            +fAZ_BLK(iBlock)*&
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
    OneDotMDotOne=cZero
    do iBlock=1,nBlock
       if (unusedBLK(iBlock)) CYCLE
       !     Prec_CB(:,:,:,iBlock)=Prec_CB(:,:,:,iBlock)*divb_diffcoeff
       if(true_blk(iBlock))then
          OneDotMDotOne=OneDotMDotOne+sum(cOne/Prec_CB(:,:,:,iBlock))
       elseif(any(true_cell(1:nI,1:nJ,1:nK,iBlock)))then
          OneDotMDotOne=OneDotMDotOne+sum(cOne/Prec_CB(:,:,:,iBlock)&
               ,MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))
       end if
    end do
    if(nProc>1)then
       call MPI_allreduce(OneDotMDotOne,OneDotMDotOneInv, 1,  MPI_REAL, MPI_SUM, &
            iComm, iError)
       OneDotMDotOneInv=cOne/OneDotMDotOneInv
    else
       OneDotMDotOneInv=cOne/OneDotMDotOne
    end if
    if(iProc==0)write(*,*)' init_divb_cleanup finishes with OneDotMDotOneInv=',OneDotMDotOneInv
!    write(*,*)'Maxval loc and minval loc of Prec_CB=', &
!         maxval(Prec_CB(:,:,:,1:nBlock)), maxloc(Prec_CB(:,:,:,1:nBlock)), &
!         minval(Prec_CB(:,:,:,1:nBlock)), minloc(Prec_CB(:,:,:,1:nBlock)) 
  end subroutine init_divb_cleanup
  !============================================================================                         
  subroutine v_grad_phi(Phi_GB,iBlock)
    integer,intent(in)::iBlock
    real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,nBLK),intent(inout)::Phi_GB
    vDotGradX_C=cZero;vDotGradY_C=cZero;vDotGradZ_C=cZero
!!! Apply continuous solution at east and west
!    if (NeiLev(East_,iBlock)==NOBLK)&
!         Phi_GB(0   ,1:nJ,1:nK,iBlock) = Phi_GB(1 ,1:nJ,1:nK,iBlock)
!    if (NeiLev(West_,iBlock)==NOBLK)&
!         Phi_GB(nI+1,1:nJ,1:nK,iBlock) = Phi_GB(nI,1:nJ,1:nK,iBlock)
!!! Apply shearing at north and south
!    if (NeiLev(South_,iBlock)==NOBLK)&
!         Phi_GB(1:nI,0   ,1:nK,iBlock) = Phi_GB(0:nI-1,2 ,1:nK,iBlock)
!    if (NeiLev(North_,iBlock)==NOBLK)&
!         Phi_GB(1:nI,nJ+1,1:nK,iBlock) = Phi_GB(2:nI+1,nJ-1,1:nK,iBlock)
!!! Apply translation invariant solution at bottom and top
!    if (NeiLev(Bot_,iBlock)==NOBLK)&
!         Phi_GB(1:nI,1:nJ,0   ,iBlock) = Phi_GB(1:nI,1:nJ,1 ,iBlock)
!    if (NeiLev(Top_,iBlock)==NOBLK)&
!         Phi_GB(1:nI,1:nJ,nK+1,iBlock) = Phi_GB(1:nI,1:nJ,nK,iBlock)
     if (NeiLev(East_,iBlock)==NOBLK)&
          Phi_GB(0   ,1:nJ,1:nK,iBlock)=-BoundaryCoef*Phi_GB(1 ,1:nJ,1:nK,iBlock)
     if (NeiLev(West_,iBlock)==NOBLK)&
          Phi_GB(nI+1,1:nJ,1:nK,iBlock)=-BoundaryCoef*Phi_GB(nI,1:nJ,1:nK,iBlock)
     if (NeiLev(South_,iBlock)==NOBLK)&
          Phi_GB(1:nI,0   ,1:nK,iBlock)=-BoundaryCoef*Phi_GB(1:nI,1 ,1:nK,iBlock)
     if (NeiLev(North_,iBlock)==NOBLK)&
          Phi_GB(1:nI,nJ+1,1:nK,iBlock)=-BoundaryCoef*Phi_GB(1:nI,nJ,1:nK,iBlock)
     if (NeiLev(Bot_,iBlock)==NOBLK)&
          Phi_GB(1:nI,1:nJ,0   ,iBlock)=-BoundaryCoef*Phi_GB(1:nI,1:nJ,1 ,iBlock)
     if (NeiLev(Top_,iBlock)==NOBLK)&
          Phi_GB(1:nI,1:nJ,nK+1,iBlock)=-BoundaryCoef*Phi_GB(1:nI,1:nJ,nK,iBlock)
    if(body_blk(iBlock))then
       do k=1,nK;do j=1,nJ;do i=1,nI
          if(.not.true_cell(i,j,k,iBlock))CYCLE
          where(.not.true_cell(i-1:i+1:2,j,k,iBlock))&
               Phi_GB(i-1:i+1:2,j,k,iBlock)=-BoundaryCoef*Phi_GB(i,j,k,iBlock)
          where(.not.true_cell(i,j-1:j+1:2,k,iBlock))&
               Phi_GB(i,j-1:j+1:2,k,iBlock)=-BoundaryCoef*Phi_GB(i,j,k,iBlock)
          where(.not.true_cell(i,j,k-1:k+1:2,iBlock))&
               Phi_GB(i,j,k-1:k+1:2,iBlock)=-BoundaryCoef*Phi_GB(i,j,k,iBlock)
          vDotGradX_C(i,j,k)=&
               cHalf*FaX_BLK(iBlock)*(&
               Phi_GB(i+1,j,k,iBlock)-&
               Phi_GB(i-1,j,k,iBlock))
          vDotGradY_C(i,j,k)=&
               cHalf*FaY_BLK(iBlock)*(&
               Phi_GB(i,j+1,k,iBlock)-&
               Phi_GB(i,j-1,k,iBlock))
          vDotGradZ_C(i,j,k)=&
               cHalf*FaZ_BLK(iBlock)*(&
               Phi_GB(i,j,k+1,iBlock)-&
               Phi_GB(i,j,k-1,iBlock))
       end do;end do;end do
    else
       do k=1,nK;do j=1,nJ;do i=1,nI
          vDotGradX_C(i,j,k)=&
               cHalf*FaX_BLK(iBlock)*(&
               Phi_GB(i+1,j,k,iBlock)-&
               Phi_GB(i-1,j,k,iBlock))
          vDotGradY_C(i,j,k)=&
               cHalf*FaY_BLK(iBlock)*(&
               Phi_GB(i,j+1,k,iBlock)-&
               Phi_GB(i,j-1,k,iBlock))
          vDotGradZ_C(i,j,k)=&
               cHalf*FaZ_BLK(iBlock)*(&
               Phi_GB(i,j,k+1,iBlock)-&
               Phi_GB(i,j,k-1,iBlock))
       end do;end do;end do
    end if
  end subroutine v_grad_phi
end subroutine clean_divb
!===================================================================
subroutine div_3d_b1(iBlock,VecX_G,VecY_G,VecZ_G,Out_G)     
use ModSize
use ModGeometry,ONLY:body_blk, true_cell, &
    fAX_BLK, fAY_BLK, fAZ_BLK
use ModParallel,ONLY:neilev,NOBLK
use ModDivbCleanup, ONLY: BoundaryCoef
use ModNumConst
implicit none

integer,intent(in) :: iBlock
real,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn),intent(in)::&
    VecX_G,VecY_G,VecZ_G
real,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn),intent(out)::Out_G

real, dimension(0:nI+1, 0:nJ+1, 0:nK+1) :: OneTrue_G

integer :: i, j, k

!\
! Can only be used for divB diffusion and projection scheme!!!!
! DivB is multiplied by -V_Cell!!!!
! With this modification DivB[grad Phi] is a symmetric positive definite 
! operator!
!/
Out_G=cZero
if(.not.(body_blk(iBlock).or.any(neilev(:,iBlock)==NOBLK)))then 
  do k=1,nK; do j=1,nJ; do i=1,nI
     Out_G(i,j,k) = - cHalf*(&
          fAX_BLK(iBlock)*&
          (VecX_G(i+1, j, k)-VecX_G(i-1,j,k))&
          +fAY_BLK(iBlock)*&
          (VecY_G(i ,j+1, k)-VecY_G(i,j-1,k))&
          +fAZ_BLK(iBlock)*&
          (VecZ_G(i, j, k+1)-VecZ_G(i,j,k-1)) )
  end do; end do; end do
else
  where(true_cell(0:nI+1, 0:nJ+1, 0:nK+1,iBlock)) 
     OneTrue_G=cOne
  elsewhere
     OneTrue_G=cZero
  end where
  if(neilev(East_ ,iBlock)==NOBLK) OneTrue_G(0   ,:,:)=cZero
  if(neilev(West_ ,iBlock)==NOBLK) OneTrue_G(nI+1,:,:)=cZero
  if(neilev(South_,iBlock)==NOBLK) OneTrue_G(:,0   ,:)=cZero
  if(neilev(North_,iBlock)==NOBLK) OneTrue_G(:,nJ+1,:)=cZero
  if(neilev(Bot_  ,iBlock)==NOBLK) OneTrue_G(:,:,0   )=cZero
  if(neilev(Top_  ,iBlock)==NOBLK) OneTrue_G(:,:,nK+1)=cZero
  !
  !\
  ! Where .not.true_cell, all the gradients are zero
  ! In true_cell the input to gradient from the face neighbor
  ! is ignored, if the face neighbor is .not.true_cell
  !/
  !
  do k=1,nK; do j=1,nJ; do i=1,nI
     Out_G(i,j,k) = - cHalf*OneTrue_G(i,j,k)*&
                                  (+&
          fAX_BLK(iBlock)*&
          (VecX_G(i+1,j,k)-&
          BoundaryCoef*(VecX_G(i+1,j,k)-VecX_G(i,j,k))*&
          (cOne-OneTrue_G(i+1,j,k))-&
          VecX_G(i-1,j,k)-&
          BoundaryCoef*(VecX_G(i,j,k)-VecX_G(i-1,j,k))*&
          (cOne-OneTrue_G(i-1,j,k)))+&
          fAY_BLK(iBlock)*&
          (VecY_G(i,j+1,k)-&
          BoundaryCoef*(VecY_G(i,j+1,k)-VecY_G(i,j,k))*&
          (cOne-OneTrue_G(i,j+1,k))-&
          VecY_G(i,j-1,k)-&
          BoundaryCoef*(VecY_G(i,j,k)-VecY_G(i,j-1,k))*&
          (cOne-OneTrue_G(i,j-1,k)))+&
          fAZ_BLK(iBlock)*&
          (VecZ_G(i,j,k+1)-&
          BoundaryCoef*(VecZ_G(i,j,k+1)-VecZ_G(i,j,k))*&
          (cOne-OneTrue_G(i,j,k+1))-&
          VecZ_G(i,j,k-1)-&
          BoundaryCoef*(VecZ_G(i,j,k)-VecZ_G(i,j,k-1))*&
          (cOne-OneTrue_G(i,j,k-1))) &
          )
  end do; end do; end do
end if

end subroutine div_3d_b1
