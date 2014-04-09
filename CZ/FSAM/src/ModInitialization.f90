module ModInitialization

  implicit none

  private

  public :: initialize
  public :: input_param
  public :: atmoin
  public :: blk_init
contains
  
  !=================================================================================

  subroutine initialize
    use ModPar,        ONLY: myid
    use ModSundry,     ONLY: time, ctime
    use ModBval,       ONLY: bvalv
    use ModRHS,        ONLY: p_init
    use ModIteration,  ONLY: nudt
    use ModIO,         ONLY: readrst_mpi
    use ModUserSetup
    use ModFSAM,       ONLY: DoRestart
    use ModControl,    ONLY: itnow
    implicit none
    !-------------------------------------------------------------------------------

    ! input solar atmosphere and define units
    call atmoin
    
    ! input run parameters
    call input_param

    ! define grid
    call grid
       
    ! define blk matrix
    call blk_init
    
    ! initialize sundry.h and initial fields
    call field_init
    if(.not. DoRestart) then
       itnow = 0
       time  = 0.D0
       ctime = time
    else
       call readrst_mpi
       ctime = time
       if(myid==0) write(6,*) 'Restart from restart files at time t = ', time
       call bvalv
    endif
    call p_init
    call nudt
  
  end subroutine initialize
  
  !================================================================================
  
  subroutine atmoin
    ! setup solar model table and define units
    use ModPar, ONLY: gg, rgas, pi, myid
    use ModSolar
    use ModPhysunits
    use ModBack,   ONLY: rtab, te, pe, re, gacc, bvfsq, gamaad, gradad, delta, &
         gradrad, vc, ross_kappa, nptjcd
    use ModInterp,     ONLY: lint
    
    implicit none
    
    character(len=80) :: cdata(4)
    integer :: nmod, iform, nn, nrdtmd, nidtmd, ndtgng, nvar, nbccf
    integer :: i, n, ip
    real    :: fracm, yp
    real, allocatable :: datmod(:), datgng(:), bccoef(:), yvar(:,:)
    integer, allocatable :: idatmd(:)
    !------------------------------------------------------------------------------
    
    allocate(idatmd(1:30))
    allocate(datmod(1:120), datgng(1:60), bccoef(1:48), yvar(1:30,1:nptjcd+1))
    allocate(rtab(1:nptjcd), te(1:nptjcd), pe(1:nptjcd), re(1:nptjcd), &
         gacc(1:nptjcd), bvfsq(1:nptjcd), gamaad(1:nptjcd), gradad(1:nptjcd), &
         delta(1:nptjcd), gradrad(1:nptjcd), vc(1:nptjcd), ross_kappa(1:nptjcd))
    
    open(unit=17, file='gong.l4b.14', form='unformatted',status='old')
    read(17) (cdata(i),i=1,4), nmod, iform, nn, nrdtmd, nidtmd, ndtgng, nvar, &
         nbccf, (datmod(i),i=1,120), (idatmd(i),i=1,30), (datgng(i),i=1,60), &
        (bccoef(i),i=1,48), ((yvar(i,n),i=1,30),n=1,nptjcd+1)
    close(17)
    
    ! checking
    if(myid==0) then
       write(6,*) '------------------ read from gong.14b.14 --------------------'
       write(6,*) cdata(1)
       write(6,*) cdata(2)
       write(6,*) cdata(3)
       write(6,*) cdata(4)
       write(6,*) 'nmod=',nmod
       write(6,*) 'iform=',iform
       write(6,*) 'nn=',nn
       write(6,*) 'nrdtmd=',nrdtmd
       write(6,*) 'nidtmd=',nidtmd
       write(6,*) 'ndtgng=',ndtgng
       write(6,*) 'nvar=',nvar
       write(6,*) 'nbccf=',nbccf
    endif

    msol  = datmod(23)
    r_sun = datmod(24)
    r_czb = r_sun - datgng(5)*r_sun
    sigma = 5.67D-05
    lsol  = datgng(4)
    
    ! checking
    if(myid==0) then
       write(6,*) 'gg=',gg
       write(6,*) 'msol=',msol
       write(6,*) 'r_sun=',r_sun
       write(6,*) 'r_czb=',r_czb
       write(6,*) 'rgas=',rgas
       write(6,*) 'sigma=',sigma
       write(6,*) 'lsol=',lsol
    endif
    
    ! set up table
    if(myid==0) then
       open(unit=16, file='JCD.table')
       write(16,*) '1-r/rs,r,fm,T,P,rho,g,N2,gamaad,gradad,delta,gradrad,vc,kappa'
    endif
    do i=1,nptjcd
       ip = nptjcd + 1 - i
       rtab(i) = 10.d0**(yvar(2,ip))*1.d11
       te(i)   = yvar(8,ip)
       pe(i)   = yvar(9,ip)
       re(i)   = yvar(10,ip)
       fracm   = 10.d0**(yvar(1,ip))
       gacc(i) = gg*fracm*msol/rtab(i)**2
       bvfsq(i) = gacc(i)/rtab(i)*yvar(18,ip)
       gamaad(i) = yvar(14,ip)
       gradad(i) = yvar(15,ip)
       delta(i) = yvar(27,ip)
       gradrad(i) = yvar(26,ip)+gradad(i)
       vc(i) = yvar(30,ip)
       ross_kappa(i) = yvar(12,ip)
       if(myid==0) &
            write(16,'(14e20.10)') 1.d0 - rtab(i)/r_sun, rtab(i), fracm, te(i), &
            pe(i), re(i), gacc(i), bvfsq(i), gamaad(i), gradad(i), delta(i), &
            gradrad(i), vc(i), ross_kappa(i)
    enddo
    if(myid==0)  close(16)
    deallocate(datmod, idatmd, datgng, bccoef, yvar)

    ! calculate parameters in ModSolar
    call lint(rtab,re,nptjcd,r_czb,yp)
    d_czb = yp
    call lint(rtab,gacc,nptjcd,r_czb,yp)
    gacc_czb = yp
    call lint(rtab,pe,nptjcd,r_czb,yp)
    hp_czb = yp/d_czb/gacc_czb
    call lint(rtab,te,nptjcd,r_czb,yp)
    th_czb = yp
    ! checking
    if(myid==0) then
       write(6,*) 'hp_czb=',hp_czb
       write(6,*) 'd_czb=',d_czb
       write(6,*) 'th_czb=',th_czb
       write(6,*) 'gacc_czb=',gacc_czb
       write(6,*) '-----------------------------------------------------------'
    endif

    ! define units
    unit_l = hp_czb
    unit_b = 1.D5
    unit_d = d_czb
    unit_temp = th_czb
    unit_v = unit_b/sqrt(4.D0*pi*unit_d)
    unit_t = unit_l/unit_v
    unit_s = unit_v**2/unit_temp
    unit_p = unit_d*unit_v**2
    unit_flux = unit_d*unit_v**3*unit_l**2
    unit_eng  = unit_d*unit_v**2*unit_l**3

  end subroutine atmoin
  
  !================================================================================

  subroutine input_param
    use ModPar, ONLY: gg, rgas, kboltz, mproton, pi, myid
    use ModSolar
    use ModPhysunits
    use ModDomain
    use ModCoeff
    use ModOutfile
    use ModControl
    use ModIoUnit,   ONLY: UnitTmp_
    use ModReadParam
    use ModFSAM
    implicit none

    real :: omega_in, nu_in, kdiff_in, eta_in
    real :: tout_in, tend_in
    character (len=100) :: NameCommand, NameItNow
    !-----------------------------------------------------------------------------
    
    if(IsStandAlone)then
       call read_file('PARAM.in',iComm)
       call read_init(' ', iSessionIn=1, iLineIn=0)
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case('#STOP')
          call read_var('itmax', itmax)
          call read_var('tend',tend_in)
          call read_var('wtLimit', wtLimit)
       case('#GRID')
          call read_var('rmin', rmin)
          call read_var('rmax', rmax)
          call read_var('thmin', thmin)
          call read_var('thmax', thmax)
          call read_var('phmin', phmin)
          call read_var('phmax', phmax)
       case('#ROTATION')
          call read_var('omega', omega_in)
       case('#PARAMETER')
          call read_var('nu',nu_in)
          call read_var('kdiff',kdiff_in)
          call read_var('eta',eta_in)
       case('#RESTART')
          call read_var('DoRestart', DoRestart)
          call read_var('ifile', ifile)
          call read_var('itnow', itnow)
       case('#HEATCONDUCTION')
          call read_var('UseImplCond', UseImplCond)
       case('#ERRPSOLV')
          call read_var('ErrPSolv', ErrPSolv)
       case('#SAVEPLOT')
          call read_var('itintv', itintv)
          call read_var('itout', itout)
          call read_var('tout', tout_in)
       case('#FILENAME')
          call read_var('rstfile', rstfile)
          call read_var('b1file', b1file)
          call read_var('b2file', b2file)
          call read_var('b3file', b3file)
          call read_var('v1file', v1file)
          call read_var('v2file', v2file)
          call read_var('v3file', v3file)
          call read_var('sfile',  sfile)
          call read_var('pfile', pfile)
          call read_var('gridfile',gridfile)
          call read_var('engfile',engfile)
       end select
    enddo
    
    write(NameItNow,'("_n",I8.8)') itnow
    engfile = trim(engfile)//trim(NameItNow)
    tout = tout_in/unit_t
    tend = tend_in/unit_t
    ovrb = omega_in*unit_l/unit_v
    ovre = nu_in/(unit_l*unit_v)
    ovrt = kdiff_in/(unit_l*unit_v)
    ovrm = eta_in/(unit_l*unit_v)
    
    ! output physical parameters and units
    if(myid==0) then
       open(UnitTmp_,file='physparams.dat', form='unformatted',access='stream',&
            status='replace')
       write(UnitTmp_) gg,rgas,kboltz,mproton,pi
       write(UnitTmp_) msol,lsol,sigma,r_sun,r_czb, hp_czb,d_czb,th_czb,gacc_czb
       write(UnitTmp_) unit_l,unit_b,unit_d,unit_temp, &
            unit_v,unit_t,unit_s,unit_p,unit_flux,unit_eng
       write(UnitTmp_) omega_in, nu_in, kdiff_in, eta_in
       close(UnitTmp_)
    endif
    
  end subroutine input_param
  
  !===============================================================================
  
  subroutine blk_init
    use ModPar,       ONLY: inmax, jnmax
    use ModGrid
    use ModBack,      ONLY: fact
    use ModBlk
    implicit none
    
    integer :: i, j
    !-----------------------------------------------------------------------------
    
    do i=is,inmax-3
       cr(i-is+1) = fact(i)/dxxa(i)*0.5D0*(1.D0/fact(i+1)+1.D0/fact(i)) &
            *g2xxa(i+1)*g31xxa(i+1)/dxxb(i+1)
       ar(i-is+1) = fact(i)/dxxa(i)*0.5D0*(1.D0/fact(i)+1.D0/fact(i-1)) &
            *g2xxa(i)*g31xxa(i)/dxxb(i)
       br(i-is+1) = - ar(i-is+1) - cr(i-is+1)
    enddo
    
    do j=js,jnmax-3
       cth(j-js+1) = 1.D0/(g32yyb(j)*dyya(j))*g32yya(j+1)/dyyb(j+1)
       ath(j-js+1) = 1.D0/(g32yyb(j)*dyya(j))*g32yya(j)/dyyb(j)
       bth(j-js+1) = - cth(j-js+1) - ath(j-js+1)
       bthk(j-js+1) = dx3ai(ks)*dx3ai(ks)/g32yyb(j)/g32yyb(j)
    enddo
    
  end subroutine blk_init
  
end module ModInitialization
