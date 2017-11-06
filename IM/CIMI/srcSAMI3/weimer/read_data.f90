!-----------------------------------------------------------------------
module read_data
  implicit none
!
! Data read from W05scEpot.dat or W05scBpot.dat:
!
  integer,parameter :: csize=28, d1_pot=15, d2_pot=18
  integer :: ab(csize), ls(csize), ms(csize)
  real :: alschfits(d2_pot,csize), schfits(d1_pot,csize), ex_pot(2)
  integer :: maxl_pot,maxm_pot
!
! Data read from SCHAtable.dat:
!
  integer,parameter :: d1_scha=19, d2_scha=7, d3_scha=68
  real :: allnkm(d1_scha,d2_scha,d3_scha)
  integer :: maxk_scha,maxm_scha
  real :: th0s(d3_scha)
!
! Data read from W05scBndy.dat:
!
  integer,parameter :: na=6, nb=7
  real :: bndya(na),bndyb(nb),ex_bndy(2)
!
  contains 
!-----------------------------------------------------------------------
subroutine read_potential(infile)
!
! Read ascii data file W05scEpot.dat or W05scBpot.dat, written by 
!   pro write_potential (write_data.pro)
!
  implicit none
!
! Args:
  character(len=*),intent(in) :: infile 
!
! Local:
!
  character(len=16) :: fname
  integer :: i,lu=20
  integer :: csize_rd,d1_rd,d2_rd
!
!!  PRINT *,infile
  open(lu,file=infile,status='old')
  read(lu,"(a)") fname
  read(lu,"(28i3)") ab
  read(lu,"(3i3)") csize_rd,d1_rd,d2_rd
  if (csize_rd /= csize) then 
    write(6,"('>>> read_potential: file ',a,': incompatable csize: ',&
      &'csize_rd=',i4,' csize=',i4)") fname,csize_rd,csize
    stop 'csize'
  endif
  if (d1_rd /= d1_pot) then 
    write(6,"('>>> read_potential: file ',a,': incompatable d1: ',&
      &'d1_rd=',i4,' d1_pot=',i4)") fname,d1_rd,d1_pot
    stop 'd1'
  endif
  if (d2_rd /= d2_pot) then 
    write(6,"('>>> read_potential: file ',a,': incompatable d2: ',&
      &'d2_rd=',i4,' d2_pot=',i4)") fname,d2_rd,d2_pot
    stop 'd2'
  endif
  do i=1,csize
    read(lu,"(6e20.9)") alschfits(:,i)
  enddo
  read(lu,"(2f10.3)") ex_pot
  read(lu,"(28i3)") ls
  read(lu,"(2i3)") maxl_pot,maxm_pot
  read(lu,"(28i3)") ms
  do i=1,csize 
    read(lu,"(6e20.9)") schfits(:,i)
  enddo

!!  write(6,"(/,'read_potential: Opened file ',a)") infile
!!  write(6,"('ab=',28i3)") ab
!!  write(6,"('csize=',i4,' d1=',i4,' d2=',i4)") csize,d1_pot,d2_pot
!!  write(6,"('alschfits min,max=',2e12.4)") minval(alschfits),maxval(alschfits)
! do i=1,csize
!   write(6,"('i=',i3,' alschfits(:,i)=',/,(6e20.9))") i,alschfits(:,i)
! enddo
! write(6,"('ex_pot=',2f10.3)") ex_pot
!!  write(6,"('ls=',28i3)") ls
!!  write(6,"('maxl_pot,maxm_pot=',2i4)") maxl_pot,maxm_pot
!!  write(6,"('ms=',28i3)") ms
!!  write(6,"('schfits min,max=',2e12.4)") minval(schfits),maxval(schfits)
! do i=1,csize 
!   write(6,"('i=',i3,' schfits(:,i)=',/,(6e20.9))") i,schfits(:,i)
! enddo
  close(lu)
!!  write(6,"('read_potential: Closed file ',a)") infile
end subroutine read_potential
!-----------------------------------------------------------------------
subroutine read_schatable(infile)
!
! Read ascii data file SCHAtable.dat, written by pro write_scha
!   (write_data.pro)
!
  implicit none
!
! Args:
  character(len=*),intent(in) :: infile 
!
! Local:
!
  character(len=16) :: fname
  integer :: i,j,lu=20
!
  open(lu,file=infile,status='old')
  read(lu,"(a)") fname
  read(lu,"(2i3)") maxk_scha,maxm_scha
  do i=1,d3_scha
    do j=1,d2_scha
      read(lu,"(6e20.9)") allnkm(:,j,i)
    enddo
  enddo
  read(lu,"(8f10.4)") th0s
!
! write(6,"(/,'read_schatable: Opened file ',a)") infile
! write(6,"('maxk_scha=',i4,' maxm_scha=',i4)") maxk_scha,maxm_scha
! do i=1,d3_scha
!   do j=1,d2_scha
!     write(6,"('i=',i3,' j=',i3,' allnkm(:,j,i)=',/,(6e20.9))") &
!       i,j,allnkm(:,j,i)
!   enddo
! enddo
! write(6,"('th0s=',/,(8f10.4))") th0s
  close(lu)
! write(6,"('read_schatable: Closed file ',a)") infile
end subroutine read_schatable
!-----------------------------------------------------------------------
subroutine read_bndy(infile)
!
! Read ascii data file W05scBndy.dat, written by pro write_bndy
!   (write_data.pro)
!
  implicit none
!
! Args:
  character(len=*),intent(in) :: infile 
!
! Local:
!
  character(len=16) :: fname
  integer :: rd_na,rd_nb,lu=20
!
  open(lu,file=infile,status='old')
  read(lu,"(a)") fname
  read(lu,"(2i3)") rd_na,rd_nb
  if (rd_na /= na) then 
    write(6,"('>>> read_potential: file ',a,': incompatable na: ',&
      &'rd_na=',i4,' na=',i4)") fname,rd_na,na
    stop 'na'
  endif
  if (rd_nb /= nb) then 
    write(6,"('>>> read_potential: file ',a,': incompatable nb: ',&
      &'rd_nb=',i4,' nb=',i4)") fname,rd_nb,nb
    stop 'nb'
  endif
  read(lu,"(8e20.9)") bndya
  read(lu,"(8e20.9)") bndyb
  read(lu,"(8e20.9)") ex_bndy
!
! write(6,"(/,'read_bndy: Opened file ',a)") infile
! write(6,"('na=',i3,' nb=',i3)") na,nb
! write(6,"('bndya=',/,(8e20.9))") bndya
! write(6,"('bndyb=',/,(8e20.9))") bndyb
! write(6,"('ex_bndy=',(8e20.9))") ex_bndy
  close(lu)
! write(6,"('read_bndy: Closed file ',a)") infile
end subroutine read_bndy
!-----------------------------------------------------------------------
end module read_data
