module ModGetQtys
  private
  
  public :: get_global

  public :: get_mean_entropy

contains
  
  !===================================================================================

  subroutine get_global(em, ek, ek1, ek2, ek3, eth, anglm, anglmnorm)
    use ModPar,      ONLY: myid1, in
    use ModGrid,     ONLY: is, ie, js, je, ks, ke, dx1a, g2b, dx2a, g31b, g32b, dx3a
    use ModField,    ONLY: v1, v2, v3, b1, b2, b3, s
    use ModBack,     ONLY: temp0, d0
    use ModMpi
    use ModFSAM,     ONLY: iComm
    implicit none
    
    real, intent(out) :: em, ek, ek1, ek2, ek3, eth, anglm, anglmnorm
    
    integer :: i, j, k, ierr
    real :: dVol, Density, energylc(1:8), energy(1:8)
    !--------------------------------------------------------------------------------
    
    energy   = 0.d0
    energylc = 0.d0
    
    do k=ks,ke; do j=js,je; do i=is,ie
       Density =  d0(i+myid1*(in-5))
       dVol = dx1a(i)*g2b(i)*dx2a(j)*g31b(i)*g32b(j)*dx3a(k)
       energylc(1) = energylc(1) + 0.125d0*dVol*((b1(i,j,k) + b1(i+1,j,k))**2 + &
            (b2(i,j,k) + b2(i,j+1,k))**2 + (b3(i,j,k) + b3(i,j,k+1))**2)
       energylc(3) = energylc(3) + 0.125d0*dVol*Density*(v1(i,j,k)+v1(i+1,j,k))**2
       energylc(4) = energylc(4) + 0.125d0*dVol*Density*(v2(i,j,k)+v2(i,j+1,k))**2
       energylc(5) = energylc(5) + 0.125d0*dVol*Density*(v3(i,j,k)+v3(i,j,k+1))**2
       energylc(6) = energylc(6) + Density*temp0(i+myid1*(in-5))*s(i,j,k)*dVol
       energylc(7) = energylc(7) + 0.5d0*dVol*Density* &
            (v3(i,j,k) + v3(i,j,k+1))*g31b(i)*g32b(j)
       energylc(8) = energylc(8) + 0.5d0*dVol*Density* &
            abs(v3(i,j,k) + v3(i,j,k+1))*g31b(i)*g32b(j)
    enddo; enddo; enddo
    energylc(2) = sum(energylc(3:5))
    call MPI_ALLREDUCE(energylc, energy, 8, MPI_DOUBLE_PRECISION, MPI_SUM, &
         iComm, ierr)
    
    em  = energy(1)
    ek  = energy(2)
    ek1 = energy(3) 
    ek2 = energy(4)
    ek3 = energy(5)
    eth = energy(6)
    anglm     = energy(7)
    anglmnorm = energy(8)
    
  end subroutine get_global

  !===================================================================================

  subroutine get_mean_entropy(smeanbot, smeantop)
    use ModPar,    ONLY: myid1, nproc1
    use ModGrid,   ONLY: g32b, dx2a, dx3a, js, je, ks, ke, is, ie
    use ModField,  ONLY: s
    use ModMpi
    use ModFSAM,   ONLY: iComm
    implicit none
    
    real, intent(out) :: smeanbot, smeantop
    
    real :: smeanbotlc, smeantoplc, arealc, area
    integer :: ierr, j, k
    !--------------------------------------------------------------------------------
    
    smeanbotlc = 0.d0
    smeantoplc = 0.d0
    arealc     = 0.d0
    if(myid1==0)then
       do k=ks,ke; do j=js,je
          smeanbotlc = smeanbotlc + s(is,j,k)*g32b(j)*dx2a(j)*dx3a(k)
          arealc     = arealc + g32b(j)*dx2a(j)*dx3a(k)
       enddo; enddo
    end if
    
    if(myid1==(nproc1-1))then
       do k=ks,ke; do j=js,je
          smeantoplc = smeantoplc + s(ie,j,k)*g32b(j)*dx2a(j)*dx3a(k) 
       enddo; enddo
    end if
    
    call MPI_ALLREDUCE(smeanbotlc, smeanbot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         iComm, ierr)
    call MPI_ALLREDUCE(smeantoplc, smeantop, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         iComm, ierr)
    call MPI_ALLREDUCE(arealc, area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         iComm, ierr)
    
    smeanbot = smeanbot/area
    smeantop = smeantop/area
    
  end subroutine get_mean_entropy

end module ModGetQtys
