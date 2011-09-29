module M_EOS
contains
  !------
  subroutine set_kbR(ro,Natom)
    use M_localProperties,only : kBr_E,kBr_P ,ERGperEV,DYNEperEV &
         ,avogadro,EVperK,Atomass,kB_ztf_E,kB_ztf_P	! ,ro,NI
    implicit none
    real,optional,intent(IN) :: ro,Natom
    
    if(present(ro) ) then
       kBr_E= ro * kB_ztf_E
       kBr_P= ro * kB_ztf_P
    else
       kBr_E= Atomass * (Natom/avogadro) * kB_ztf_E
       kBr_P= Atomass * (Natom/avogadro) * kB_ztf_P
    end if
    return
  end subroutine set_kBr
  !=======================
  subroutine setOptions(brent,EElog,caleos)
    use M_NLTE,only : useZbrent,useEElog
    implicit none
    logical,optional,intent(IN) :: brent,EElog,caleos
    
    if(present(brent)) useZbrent=brent
    if(present(EElog)) useEElog=EElog
    ! write(0,*) '- - useEElog,useZbrent=',useEElog,useZbrent
    return
  end subroutine setOptions
end module M_EOS

 !------
subroutine setMethod
  use M_EOS, ONLY: setOptions

  implicit none
  integer ib,ie

  write(0,*)'input 2 values : 1 for useZbrent, 1 for useEElog'
  read(*,*) ib,ie
  call setOptions(ib.gt.0,ie.gt.0,.false.)

  		write(0,*)'calls setZTF(au)'
  		
  ! initialize Thomas_Fermi EOS (analytical fit)
  call setZTF('au')
  !===================================================
  !Check the consitency of the input parameters
  call verify()
return
end




 ! voids function, replacing interface to CALEOS
 
subroutine inital
 implicit none
 stop '"inital" not implemented'
end

subroutine ltesta(e,p,r,t,z,cv,mode1)
 implicit none
 integer,intent(IN) :: mode1
 real,intent(IN) :: e
 real,intent(OUT) :: p,r,t,z,cv
 stop '"ltesta" not implemented'
end

subroutine CALEOS_dir(te,ro,Etot,Ptot,Zbar,Cv)
 implicit none
 real,intent(IN) :: te,ro
 real,intent(OUT) :: Ptot,Zbar,Cv
 real,intent(INOUT) :: Etot
 stop '"CALEOS_dir" not implemented'
end

subroutine CALEOS_inv(te,ro,Etot,Ptot,Zbar,Cv)
 implicit none
 real,intent(IN) :: te,ro
 real,intent(OUT) :: Ptot,Zbar,Cv
 real,intent(INOUT) :: Etot
 stop '"CALEOS_inv" not implemented'
end

