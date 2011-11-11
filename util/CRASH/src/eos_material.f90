!=====================Above is for ZTF
!HERE CAL means CACculated and tabulated EOS in contrast with
!ZTF which is the analytical fit for the Thomas-Fermi model
module CRASH_M_EOS
  !The CRASH meterial number
  SAVE
  integer:: iMaterial = -1
  logical:: UseCrashEos = .false.
contains
  !------
  subroutine set_kbR(ro,Natom)
    use CRASH_M_localProperties,only : kBr_E,kBr_P ,ERGperEV,DYNEperEV &
         ,avogadro,EVperK,Atomass,kB_ztf_E,kB_ztf_P	! ,ro,NI
    use ModConst,ONLY: cBoltzmann,cEV,cKToEv,cEVToK, cAtomicMass
    use CRASH_ModEos,ONLY:cAtomicMassCRASH_I
    implicit none
    real,optional,intent(IN) :: ro,Natom
    if(useCrashEos)then
       if(present(ro) ) then
          kBr_E= ro/(cAtomicMass * cAtomicMassCRASH_I(iMaterial))& !Natomic
               * cBoltzmann * cEVToK 
          kBr_P= kBr_E
       else
          kBr_E= Natom * cBoltzmann * cEVToK 
          kBr_P= Natom * cBoltzmann * cEVToK 
       end if
    else
       if(present(ro) ) then
          kBr_E= ro * kB_ztf_E
          kBr_P= ro * kB_ztf_P
       else
          kBr_E= Atomass * (Natom/avogadro) * kB_ztf_E
          kBr_P= Atomass * (Natom/avogadro) * kB_ztf_P
       end if
    end if
    return
  end subroutine set_kBr
  !=======================
  subroutine setOptions(brent,EElog,caleos)
    use CRASH_M_NLTE,only : useZbrent,useEElog
    implicit none
    logical,optional,intent(IN) :: brent,EElog,caleos
    
    if(present(brent)) useZbrent=brent
    if(present(EElog)) useEElog=EElog
    ! write(0,*) '- - useEElog,useZbrent=',useEElog,useZbrent
    return
  end subroutine setOptions
end module CRASH_M_EOS
!===============
!  to develop and test "radiom.f90"  we use this simple EOS based on
!  Thomas-Fermi  fit  by  F.Perrot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! \
!!   PROBABLY IT SHOULD BE BETTER TO PASS "ro" through argument all along the hierachy of calls
!! /
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module CRASH_M_ZTF
  use CRASH_M_EOS, ONLY: useCALEOS =>useCrashEos
  !This module is used in this file only
  implicit none
  integer,save :: nbCold=0
  real,save,allocatable,dimension(:) :: ROcolds,EEcolds,TEcolds
  real,save :: atoNumCold=0,atoMassCold=0  ,refden1=0,refden2=0
  character*(80),save :: fileCold=' '
  real,save :: zt		! shared variable between "ZTF_EOS_..." and ZTFdiff
  logical,save :: setCALEOS=.false.
  character(16),save :: endianCALEOS='LITTLE_ENDIAN'
  !===================Added October, 3rd, 2011=====================================
  logical,save:: UseCRASHEos = .false.
  !------
contains
  !------
  subroutine getEcold_table(AtoNum)
    implicit none
    real,intent(IN) :: atoNum
    integer :: ier,i
    real :: TEcold
    character*(4) :: symboldCold=' '
    !
    if(fileCold.eq.' ') then
       call CON_stop('-E- getEcold_table : file not assigned')
    end if
    open(1,file=fileCold,status='OLD',form='FORMATTED',iostat=ier)
    if(ier.ne.0) then
       write(*,*) '-E- cannot open "cold material file " : ',fileCold &
            ,', I/O error number=',ier
       stop '-I/O cannot open file'
    end if
    read(1,*)
    read(1,*) atoNumCold,atoMassCold,symboldCold	
    read(1,*) TEcold
    read(1,*) nbCold
    read(1,*)
    if(allocated(ROcolds)) deallocate(ROcolds,EEcolds,TEcolds)
    allocate(ROcolds(nbCold),EEcolds(nbCold),TEcolds(nbCold),stat=ier)
    if(ier.ne.0) then
       stop '-E- cannot allocate ROcolds,EEcolds,TEcolds'
    end if
    do i=1,nbCold
       read(1,*) ROcolds(i),EEcolds(i)
    end do
    TEcolds(:)=TEcold
    close(1,iostat=ier)
    !
    return
  end subroutine getEcold_table
  !------
  subroutine getEcoldZTF(atoNum,ro,EEcold,TEcold)
    implicit none
    real,intent(IN) :: atoNum,ro
    real,intent(OUT) :: EEcold,TEcold
    real :: r,rn
    integer n,i,ilo
    !
    if(atoNum.ne.atoNumCold) call getEcold_table(atoNum)
    if(ro.le.ROcolds(1)) then
       EEcold=EEcolds(1)
       TEcold=TEcolds(1)
    elseif(ro.gt.ROcolds(nbCold))then
       EEcold=EEcolds(nbCold)
       TEcold=TEcolds(nbCold)
    else
       r=log(ROcolds(2)/ROcolds(1))
       n=1+int(log(ro/ROcolds(1)) / r )
       if(ROcolds(n).le.ro) then
          do i=n,nbCold-1
             ilo=i
             if(ROcolds(i+1).gt.ro) exit
          end do
       else
          do i=n-1,1,-1
             ilo=i
             if(ROcolds(i).le.ro) exit
          end do
       end if
       r =ro/ROcolds(ilo)
       rn=ROcolds(ilo+1)/ROcolds(ilo)
       r=log(r)/log(rn)
       EEcold=EEcolds(ilo) + r * (EEcolds(ilo+1)-EEcolds(ilo))
       TEcold=TEcolds(ilo) + r * (TEcolds(ilo+1)-TEcolds(ilo))
    end if
    !
    return
  end subroutine getEcoldZTF
  !------
  subroutine ZTF_EOS_dir(te,Etot,Ptot,Zbar,Cv)	! ro : in module M_localProperties

    use CRASH_M_localProperties,only : ro,atoNum,Atomass,kBr_E,kBr_P &
         ,ERGperEV,DYNEperEV ,Pcold,Efloor,zion
    implicit none
    real,intent(IN) :: te
    real,intent(OUT) :: Etot,Ptot,Zbar,Cv
    real,parameter :: x_3o2=1.5d0,x_1=1d0
    real,parameter :: smallZ=0.05
    !
    ! args. passed through MODULE:
    !
    call ZTF_Perrot(atoNum,AtoMass,te,ro,zbar)
    zbar=max(zbar,smallZ)
    Cv  = x_3o2 * (zbar+zion)      * kBr_E
    zt  = (zbar+zion) * te
    Etot= x_3o2 * zt * kBr_E + Efloor*ro
    Ptot=         zt * kBr_P + Pcold
    !
    return
  end subroutine ZTF_EOS_dir
  !------
  subroutine ZTF_EOS_inv(te,Etot,Ptot,Zbar,Cv)	! ro : in module M_localProperties

    use CRASH_M_localProperties,only : ro,atoNum,Atomass,kBr_E,kBr_P &
         ,ERGperEV,DYNEperEV ,Pcold,Efloor,zion	! ,zt
    implicit none
    real,intent(IN) :: Etot
    real,intent(OUT) :: te,Ptot,Zbar,Cv
    real,parameter :: x_0=0,x_3o2=1.5d0,x_1=1d0,x_1o2=0.5d0,two=2d0,three=3d0
    real,parameter :: zmin=1e-4 , tol=1d-4
    real :: z1,z2 ,Ecold,Tcold
    !
    !
    ! args. passed through MODULE:
    !
    call getEcold(ro,Ecold,Tcold)		! mb.110827
    !
    zt=(Etot-Efloor*ro)/(x_3o2*kBr_E)  !  == (zbar+zion) * te
    if(zt.le.0 .or. Etot.le.Ecold) then		! mb.110827
       call ZTF_More_cold(atoNum,AtoMass,ro,zbar)
       te=Tcold
       Ptot= Pcold
       return
    end if

    z1=zmin
    z2=atoNum
    zbar=zbrentZM(z1,z2,zmin,tol,te)		! cvg until |Znew-Zold| < tol
    te=zt/(zbar+zion)
    Cv  = x_3o2 * zt/te      * kBr_E
    Ptot=         zt * kBr_P + Pcold
    !
    return
  end subroutine ZTF_EOS_inv

  !------
  function ZTF_dif(z_est,t_est)
    use CRASH_M_localProperties,only : atoNum,AtoMass,ro,zion
    implicit none
    real :: z_est,ZTF_dif,z_new,t_est 

    t_est=zt/(z_est+zion)
    call ZTF_Perrot(atoNum,AtoMass,t_est,ro,z_new)
    ZTF_dif=z_new - z_est

    return
  end function ZTF_dif
  !------
  function zbrentZM(x1,x2,xmin,tol,te)
    !
    !- this routine was taken from Numerical Recipes. it uses the Brent
    !- method to find a zero of a function without given derivative.
    !- we use here this function to iterate on the chemical potential
    !- till convergence  : ZTF_dif  variation < tol  (in abs.value, not relative)
    !-    when not bracketed, extension of the range has been added	(mb.031206)
    !-    using a minimal bound  "xmin"				(mb.031206)
    !
    implicit none
    !
    real :: x1,x2,xmin,tol,te
    !
    real :: zbrentZM
    !
    integer,parameter ::  itmax=100
    integer ::iter
    real,parameter :: epsZbrent=3.D-8 &
         , tenth=0.1d0,zero=0,half=0.4d0,one=1d0,two=2d0,three=3d0,x1o2=0.5d0
    !
    real :: a,b,c,d,e,fa,fb,fc,tenn
    real :: p,q,r,s,xm,tol1
    integer,save :: nbZbrent=0
    logical :: dbg=.false.
    !
    e=0  ! as this may not be initialized

    zbrentZM=(x1+x2)*half	! 080715
    !
    a=x1
    b=x2
    if(a.eq.b) then	
       a=min(a,b)-0.1	
       b=a+0.1	
    end if
21  fa=ZTF_dif(a,te)
    if(fa.eq.zero) then	! 100629
       zbrentZM=a		! 100629
       return			! 100629
    end if			! 100629
22  fb=ZTF_dif(b,te)
    if(fb.eq.zero) then	! 100629
       zbrentZM=b		! 100629
       return			! 100629
    end if			! 100629
1010 iter=0			! 100629
2   if((fb*fa).gt.zero)then	! 100629
       if(b.gt.10.) then
          b=1.d0
          fb=ZTF_dif(b,te)
          goto 2
       end if
       write(*,*)'-R- zbrentZM: root must be bracketed'
       stop 'stop  zbrentZM .not. bracketed'
    end if
    !
100 format(5x,'n = ',i4,/,5x,'x1 = ',e12.5,/,5x,'x2 = ',e12.5 &
         ,/,5x,'f1 = ',e12.5,/,5x,'f2 = ',e12.5)
    !
    iter=0	! 051025
    !-  extend range if required
10  tenn=max(tenth*abs(a),tenth)
    if((fa*fb).gt.0) then
       if(iter.ge.itmax)then	! 051025
          goto 25
       end if	! 051025
       iter=iter+1
       if(a.lt.b)then
          a=max(xmin,a-tenn)
          b=b+tenn
       else
          b=max(xmin,b-tenn)
          a=a+tenn
       end if
       fb=ZTF_dif(b,te)
       fa=ZTF_dif(a,te)
       goto 10
    end if	! (fb*fa).gt.zero
    !
    !-  use  Brent algorithm
20  c=a
    fc=fa
    d = b - a
    e = d
    LOOP11: do iter=1,itmax
       if((fb*fc).gt.zero) then	! at iter#1 fb=fc
          c=a
          fc=fa
          d=b-a
          e=d
       end if
       if(abs(fc).lt.abs(fb)) then	! at iter#1 fb=fc
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       end if
       tol1=epsZbrent*(abs(b)+abs(a))		! 110902
       xm=half*(c-b)
       if(abs(xm).le.tol1 .or. fb.eq.zero)then
8         zbrentZM=b
          return
       elseif(fa.eq.zero) then		! mb
          zbrentZM=a			! mb
          return
       elseif(fb.eq.zero) then		! mb
          zbrentZM=b			! mb
          return
       end if
       if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
	     p=2*xm*s
	     q=one-s
          else
	     q=fa/fc
	     r=fb/fc
	     p=s*(2*xm*q*(q-r)-(b-a)*(r-one))
	     q=(q-one)*(r-one)*(s-one)
          end if
          if(p.gt.zero) then
	     q=-q
          elseif(p.le.zero) then
	     p=-p
          end if
          if(2*p .lt. min(3*xm*q-abs(tol1*q),abs(e*q))) then
             e=d
             d=p/q
          else
             d=xm
             e=d
          end if
       else	! (abs(e).ge.tol1 ...
          d=xm
          e=d
       end if	! (abs(e).ge.tol1 ...
       a=b
       fa=fb
       if(abs(d) .gt. tol1) then
          b=b+d
       else
          b=b+sign(tol1,xm)
       end if
       b=max(b,xmin)		! mb.080715
       fb=ZTF_dif(b,te)
    end do LOOP11
    !
    write(* ,111) itmax
111 format('-R- zbrentZM exceeding maximum iterations' &
         ,',(=',i5,')')
    call CON_stop('-R- no convergence in zbrentZM')
    zbrentZM=b
    return
    !- out of range
25  write(* ,112) nbZbrent,x1,x2,tol,a,b,fa,fb,c,fc
112 format('-R- from ',i8,'# zbrentZM' &
         ,' not bracketed ' &
         ,/,'x1,x2,tol,lfix=',1p,3e12.3,l2 &
         ,/,' a,b,fa,fb=',1p,4e12.3,' c,fc=',2e12.3)
    call CON_stop('-R- no bracketing in zbrentZM')
  end function zbrentZM

  SUBROUTINE ZTF_Perrot(atoNum,AtoMass,Te_eV,RO_g_cm3, Zbar)
    implicit none
    real, intent(IN) :: AtoMass,Te_eV,RO_g_cm3
    real, intent(IN) :: atoNum
    real, intent(out) :: Zbar
    real,parameter :: x_1=1d0,x_2=2d0,x_1o3=1.d0/3.d0,x_4o3=4.d0/3.d0,x_1o2=0.5d0
    real,parameter :: c1=1.31,c2=-0.64,g1=17.96,g2=18.37  &
         ,y1=2.22e-3,y2=1.3,y3=0.7
    real,parameter :: Alpha=33.12,Beta=0.6847
    real :: t,t2,t3,y,c,d,g,gc,gd,zz,zb,x,r
    real,parameter :: xmaxx=1d0-1d-4		! 110819
    r=RO_g_cm3/(atoNum*atoMass)
    t=Te_eV/atoNum**x_4o3
    t2=t/(x_1+t)
    t3= t**x_1o3
    y=y1*t*(y2*t+t3)/(x_1+y3*t3)
    c=c1 + c2*t2
    g=g1/sqrt(x_1 + g2*t2)
    gc= g**c
    gd=gc-x_1
    d=c*gc/gd
    zz=r*(x_1 + gd*(y/r)**d)**(x_1/c)
    zb=zz**Beta
    x=Alpha *zb*(x_1 + zb/x_2)
    !
    Zbar= atoNum*min(xmaxx, x/(x_1+x+sqrt(x_1 + 2*x)))
    !
    RETURN
  END SUBROUTINE ZTF_Perrot

  !------
  SUBROUTINE ZTF_More_Cold(atoNum,AtoMass,RO_g_cm3, Zbar)
    implicit none
    real, intent(IN) :: AtoMass,RO_g_cm3
    real, intent(IN) :: atoNum
    real, intent(out) :: Zbar
    real,parameter :: x_1=1d0,x_2=2d0,x_1o3=1.d0/3.d0,x_4o3=4.d0/3.d0
    real,parameter :: c1=1.31,c2=-0.64,g1=17.96,g2=18.37  &
         ,y1=2.22e-3,y2=1.3,y3=0.7
    real,parameter :: b0=-1.763,b1=1.43715,b2=0.31546
    real,parameter :: Alpha=14.3139,Beta=0.6624
    real :: R,x
    !
    R=RO_g_cm3/(atoNum*AtoMass)
    x=alpha* R**Beta

    Zbar=atoNum * x / (x_1+x+SQRT(x_1+2*x))
    !
    return
  end SUBROUTINE ZTF_More_Cold
  !------
  SUBROUTINE ZTF_More(atoNum,AtoMass,Te_eV,RO_g_cm3, Zbar)
    implicit none
    real, intent(IN) :: AtoMass,Te_eV,RO_g_cm3
    real, intent(IN) :: atoNum
    real, intent(out) :: Zbar
    real,parameter :: x_1=1d0,x_2=2d0,x_1o3=1.d0/3.d0,x_4o3=4.d0/3.d0
    real,parameter :: a1=3.323d-3,a2=0.9718d0,a3=9.26148d-5,a4=3.10165d0
    real,parameter :: b0=-1.763,b1=1.43715,b2=0.31546
    real,parameter :: Alpha=14.3139,Beta=0.6624
    real,parameter :: c1=-1.1d0/3.d0,c2=2.95d0/3.d0
    real :: T,R,A,TF,B,C,Q,Q1,x
    !
    T=Te_eV/ atoNum**x_4o3
    R=RO_g_cm3/(atoNum*AtoMass)

    A=a1 * T**a1 + a3 * T**a4
    TF=T/(x_1+T)

    B=-EXP(b0+b1*TF+b2*TF**7)
    C=C1*TF+C2
    Q1=A * R**B
    Q=( R**C + Q1*C) ** (x_1/C)

    x=alpha* Q**Beta

    Zbar=atoNum * x / (x_1+x+SQRT(x_1+2*x))
    !
    return
  end SUBROUTINE ZTF_More

  !------
end module CRASH_M_ZTF

!===below are the interface routines 
!------
subroutine verify()	! check transmission of Z,A,...
  use CRASH_M_localProperties
  implicit none
  write(*,*)'verify: Z,A=',atonum,atomass,' roS=',roSolid
  if(atonum.gt.0  .and. atomass.gt.0) return
  stop '-P-  atoNum or atoMass not set'
end subroutine verify
!------
subroutine setZTF(symb)
  use CRASH_M_localProperties,only : atoNum,Atomass,roSolid &
       ,Efloor,TEfloor,Pcold
  use CRASH_M_ZTF
  implicit none
  character*(*),optional,intent(IN) :: symb

  if(.not.present(symb)) goto 100
  if(symb.eq.'au' .or. symb.eq.'Au') then
     atoMass=196.966
     atoNum=79.
     roSolid=19.3
     Efloor=1e-15
     TEfloor=0.025d0
     Pcold=0
  else
     write(*,*)'-E- only "au" accepted in ZTF at this time'
     goto 100
  end if
  fileCold='ztf_'//symb//'.Ecold'
  call getEcold_table(atoNum)
  return
100 call CON_stop('-E- setZTF requires "atom=au"')
end subroutine setZTF

subroutine LTE_EOS_dir(te,Etot,Ptot,Zbar,Cv)	! ro : in module M_localProperties
  use CRASH_M_ZTF,only : ZTF_EOS_dir,useCALEOS,setCALEOS
  use CRASH_M_localProperties,only : atoNum,ro,zion
  use CRASH_ModEos, ONLY: eos
  use CRASH_M_EOS,ONLY:iMaterial
  implicit none
  real,intent(IN) :: te
  real,intent(OUT) :: Etot,Ptot,Zbar,Cv
  real,external :: zbrent_EE_inv

  if(useCALEOS) then
      if(zion==0) then !The electron contribution only
        call eos(iMaterial=iMaterial,&
                 Rho=ro, &
                 TeIn=Te,&
                 eElectronOut=ETot,&
                 pElectronOut=PTot,&
                 zAverageOut=zBar,&
                 CvElectronOut=Cv)
     else
        call eos(iMaterial=iMaterial,&
                 Rho=ro, &
                 TeIn=Te,&
                 eTotalOut=ETot,&
                 pTotalOut=PTot,&
                 zAverageOut=zBar,&
                 CvTotalOut=Cv)
     end if 
  else
     call ZTF_EOS_dir(te,Etot,Ptot,Zbar,Cv)
  end if

end subroutine LTE_EOS_dir	! ro : in module M_localProperties

!------
subroutine LTE_EOS_inv(te,Etot,Ptot,Zbar,Cv)	! ro : in module M_localProperties
  use CRASH_M_localProperties,only : ro, zion
  use CRASH_M_ZTF,only : ZTF_EOS_inv,useCALEOS,setCALEOS
  use CRASH_ModEos, ONLY: eos
  use CRASH_M_EOS,ONLY:iMaterial
  implicit none
  real,intent(IN) :: Etot
  real,intent(OUT) :: te,Ptot,Zbar,Cv

  if(useCALEOS) then
     if(zion==0) then !The electron contribution only
        call eos(iMaterial=iMaterial,&
                 Rho=ro, &
                 eElectronIn=ETot,&
                 TeOut=Te,&
                 pElectronOut=PTot,&
                 zAverageOut=zBar,&
                 CvElectronOut=Cv)
     else
        call eos(iMaterial=iMaterial,&
                 Rho=ro, &
                 eTotalIn=ETot,&
                 TeOut=Te,&
                 pTotalOut=PTot,&
                 zAverageOut=zBar,&
                 CvTotalOut=Cv)
     end if
  else
     call ZTF_EOS_inv(te,Etot,Ptot,Zbar,Cv)
  end if

end subroutine LTE_EOS_inv	! ro : in module M_localProperties

!------
subroutine getEcold(ro,EEcold,TEcold)
  use CRASH_M_ZTF,only : getEcoldZTF,useCALEOS,setCALEOS
  use CRASH_M_localProperties,only : atoNum

  implicit none
  real,intent(IN) :: ro
  real,intent(OUT) :: EEcold,TEcold
  real :: Zcold,CVcold,Pcold

  if(useCALEOS) then
     !At the time being, Ecold=0 and for
     !any positive energy density the LTE regime
     !is not applied. In the future, if we know
     !some good estimate, for Ecold(rho) below which
     !the LTE approximation can be employed, we can save
     !the computation time at low energy densities (or
     !temperatures)
     EEcold=0
     TEcold=0
  else
     call getEcoldZTF(atoNum,ro,EEcold,TEcold)
  end if
end subroutine getEcold
!============================================

