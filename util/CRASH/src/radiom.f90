!
!  code will use :
!
! initialisation:
!	call exp_tab8()		! prepare tabulated exponentials
!	call prepCorrUbar()	! or it will be called automatically
!	call prep_projE(hnuGr(0:nbGr),nbGR)
!	call calTZ(Te_eV,Ne_cn3, Tz, hnuGr,Erad,Brad,ng)
!  or
!	call calTZ0(Te_eV,Ne_cn3, Tz, Erad,Brad,ng)
!  or
!	call calTE(Te_estimate,Tz,Ne, Tenew, hnug,eg,bg,ng)
!

!------- radiom.f90
!
!   check array dim (0:n) / (1+n+1) : grBound(0:ngr), *ctrb
!
!------- radiom.f90
!-------
! \
MODULE CRASH_M_projE
  ! /
  !  {U=hnu/Kte,Erad,Brad} will be projected on a working grid , with even log. spacing of U
  !
  implicit none
  !\
  ! Minimal and maximal values for E/T_e ratio
  !/
  real,parameter :: Umin=1d-2, Umax=1d2
  !\
  ! Number of grid points per a decimal order
  !/ 
  real,parameter :: nbPerDec =25	! -> 100
  !\
  ! log10(uMax/uMin)
  !/
  real,parameter :: nbDec =4
  !\
  ! Number of grid points for u=E/T_e
  !/
  real,parameter :: nbUout =nbPerDec*nbDec	! =100

  integer,parameter :: mxOut= 300	
  !\
  ! log(UMax/UMin)/nBOut  - set in prep_proj
  !/
  real,save         :: lgdu
  !\
  ! exp(lgdu)
  !/
  real,save         :: rdu 

  real,save         :: Efirst,Elast
  real,save,dimension(mxOut+1) :: Eout,Uout,SPout

  !\
  ! Number of energy grouprs in the input 
  ! file of radiation group energies
  !/
  integer,save :: nbIn=0

  !\
  ! Number of grid points in the internal grid
  ! = log(Elast/Efirst)/lgdu+1
  !/
  integer,save :: nbOut=0

  integer,save :: nbContrib=0,out1,out2
  integer,parameter :: mxIn=300,mxContrib=150000
  integer,save :: to_ctrb(mxContrib),fr_ctrb(mxContrib)
  integer,save :: last_ctrb(0:mxOut),first_ctrb(1:mxOut+1)
  real,save :: coef_ctrb(mxContrib)
  integer,save :: ctrb_to1(mxContrib),ctrb_to2(mxContrib)
  logical,save :: dbgProj =.false.

  !-------
contains
  !========
  !\
  !Set debug flags
  !/
  subroutine setDbgPrep(flag)
    logical,intent(IN) :: flag
    dbgProj=flag
    !write(*,*)'- - flag  dbgProj  set to ',dbgProj
  end subroutine setDbgPrep
  !========
  subroutine getEout(E,nb)
    ! returns  boundaries of projection groups
    integer,intent(OUT) :: nb
    real,intent(OUT),dimension(:) :: E
    nb=nbOut
    E(1:nb)=Eout(1:nbOut)
  end subroutine getEout
  !========
  subroutine prep_projE(Ein,nbE)	
    !  Initialization of the projection coefficients 
    !\
    !Inputs
    !/
    !Number of energy groups
    integer,intent(In) :: nbE
    !Photon eergies
    real,intent(In),dimension(0:nbE) :: Ein

    real,parameter :: ten=10.0
    real :: u1,u2,du
    integer :: n,nIn,nOut 
    logical :: b1,b2
    integer :: fr,to
    real :: c
    logical :: dbg
    !------------
    dbg=dbgProj
    ! 
    if(dbg)write(*,*)'___ prep_projE ____  nbE=',nbE &
         ,'  E(0,',nbE,')=',Ein(0),Ein(nbE)
    nbIn=nbE
    lgdu=log(Umax/Umin)/nbUout
    rdu=exp(lgdu)
    Efirst=Ein(0)
    Elast=Ein(nbIn)
    nbOut=1+int(log(Elast/Efirst)/lgdu)
    if(dbg)write(*,*)'Elast,Efirst,lgdu,nbIn=',Elast,Efirst,lgdu,nbIn
    if(nbOut.gt.mxOut)then
       write(*,125) nbUout,nbOut
125    format(// '-E- nbUout(=',i5,') should be LESS THAN nbout(=',i5,')',//)
       call CON_stop('-R- : prep_projE  nbUout .GT. nbOUT')
    end if
    if(dbg)write(*,*)'Ef,El=',Efirst,Elast,' nbO=',nbOut
    Eout(1)=Efirst
    first_ctrb(1)=0		!!
    last_ctrb(0)=0
    do n=2,nbOut+1
       first_ctrb(n)=0	!!
       last_ctrb(n-1)=0
       Eout(n)=Eout(n-1)*rdu
    end do
    ! 
    nOut=0
    out1=0
    out2=0
    nbContrib=0
    du=Eout(nbOut+1)-Eout(nbOut)
    LOOP10:	do nIn=1,nbIn	! U(nIn-1)-U(nIn)
       last_ctrb(nIn-1)=nbContrib
       if(dbg) write(*,*)'E(',nIn,',:)='  &
            ,REAL(Ein(nIn-1)),REAL(Ein(nIn))  &
            ,' Eout(',nOut,'+1)=',REAL(Eout(nOut+1))

       do while (Ein(nIn-1).ge.Eout(nOut+1)) 
          nOut = nOut+1
          first_ctrb(nOut)=nIn	
          du=Eout(nOut+1)-Eout(nOut)
       end do
       if(Ein(nIn).le.Eout(nOut)) then
          cycle LOOP10
       end if
2      b1=Ein(nIn-1).gt.Eout(nOut)
       if(b1) then
          u1=Ein(nIn-1)
       else
          u1=Eout(nout)
          if(out1.eq.0) out1=nOut
       end if
       b2=Ein(nIn).lt.Eout(nOut+1)
       if(b2) then
          u2=Ein(nIn)
       else
          u2=Eout(nOut+1)
          out2=nOut
       end if
       nbContrib=nbContrib+1
       if(nbContrib.le.mxContrib) then
          to_ctrb(nbContrib)=nOut
          fr_ctrb(nbContrib)=nIn	! in group  [ Ein(nIn-1) : Ein(nIn) ]
          coef_ctrb(nbContrib)=(u2-u1)/du
          if(dbg) write(*,*)'#',nbContrib,' fr,to=',nIn,nOut &
               ,' u1,u2=',REAL(u1),REAL(u2) &
               ,' du=',REAL(du),' C=',REAL(coef_ctrb(nbContrib))
       end if
       !
       !   Ein :	  |      |
       !   Eout:	    |..a   b
       !   u1,u2           |..A B
       !
       !   Ein :	    |       |
       !   Eout:	  a    b..|
       !   u1,u2           A  B..|
       !
       if(.not.b2)then	! input segment overlaps 
          if(nOut.ge.nbOut) exit LOOP10	! more than 1 output segment
          nOut=nOut+1
          first_ctrb(nOut)=nIn		!!
          du=Eout(nOut+1)-Eout(nOut)
          goto 2
       end if
    end do LOOP10
    ! 
20  last_ctrb(nbIn)=nbContrib
    first_ctrb(nOut+1)=nIn		!!

    if(nbContrib.ge.mxContrib) then
       write(*,*)'-D- needs to increase "mxContrib" to ',nbContrib
       call CON_stop('')
    else
       do n=1,nbContrib
          ctrb_to1(n)=0
          ctrb_to2(n)=0
       end do
       do n=1,nbContrib
          ctrb_to2(to_ctrb(n))=n
       end do
       do n=nbContrib,1,-1
          ctrb_to1(to_ctrb(n))=n
       end do
    end if
    ! 
    if(.not.dbg) return
    !\
    !Print out the debug info
    !/
    write(*,*)'- - found ',nbContrib,' contrib.'
    WRITE(*,401) 'to_ctrb=' &
         ,(to_ctrb(nIn),nIn=1,nbContrib)
    WRITE(*,401) ' FR_ctrb=' &
         ,(fr_ctrb(nIn),nIn=1,nbContrib)
    WRITE(*,401) 'ctrb_to1=' &
         ,(ctrb_to1(nIn),nIn=1,nbContrib)
    WRITE(*,401) 'CTRB_to2=' &
         ,(ctrb_to2(nIn),nIn=1,nbContrib)
    WRITE(*,401) 'first_ctrb=' &
         ,(first_ctrb(nIn),nIn=1,nbOut+1)
    WRITE(*,401) ' LAST_ctrb=' &
         ,(last_ctrb(nIn),nIn=0,nbOut)
401 format(a,/,(20i4))
    write(*,*)'prep_projE -> nIn,nOut=',nIn,nOut &
         ,' nbContrib=',nbContrib,'  out1,out2=',out1,out2
    call wr5('Ein=',Ein(0:nbIn),nbIn+1)
    call wr5('Eout=',Eout,nbOut)
101 format('in#',i4,' E=',1p,2e13.3)
    do nIn=1,nbIn
       write(*,101) nIn,Ein(nIn-1),Ein(nIn)
       do nOut=last_ctrb(nIn-1)+1,last_ctrb(nIn)
  	  fr=fr_ctrb(nOut)
  	  to=to_ctrb(nOut)
  	  c=coef_ctrb(nOut)
  	  write(*,102) to,fr &
               ,c &
               ,Ein(fr-1),Ein(fr),Eout(to),Eout(to+1) &
               ,nOut,nIn
       end do
102    format(' to,fr=',2i4,' c=',1p,e13.3 &
            ,' Ein=',2e13.3,'  Eout=',2e13.3,' #',i4,i5)
    end do
    ! 
  end subroutine prep_projE
  !======
  subroutine projSP_E(Te,SPin,nOut,gotoLTE)
    real,intent(In) :: Te,SPin(nbIn)
    integer,intent(Out) :: nOut
    real :: Ufirst,Ulast ,uBef,uAft,c,r,u
    integer :: n1,n2,n3,ctr,fr,to,nBef,nAft,n,tt1,tt2
    real, parameter :: zero=0.d0,one=1.d0
    logical :: gotoLTE,dbg

    dbg=dbgProj
    ! 
    if(nbIn.eq.0) then
       write(*,*)'-P- subroutine "prep_projE" has not been called'
       call CON_stop('')
    end if
    Ufirst=Efirst/Te
    Ulast=Elast/Te
    gotoLTE=.false.
    if(dbg) then
       write(*,*)'-- projSP_E --','  nbIn=',nbIn
       write(*,*)'Uf,Ul=',REAL(Ufirst),REAL(Ulast) &
            ,' Umin,max=',REAL(Umin),REAL(Umax)
    end if
    !
    ! special cases:  range of interest and group span dont overlap:
    if(Ulast.le.Umin) then
       if(dbg) then
          write(*,*)'  --- |  |'
       end if
       n1=2
       n2=1
       r=zero
       goto 110
    elseif(Ufirst.gt.Umax) then
       gotoLTE=.true.
       if(dbg) then
          write(*,*)' |  |  ---'
       end if
       n1=nbContrib
       n2=nbContrib-1
       r=one
       goto 110
    end if
    !
    Uaft=max(Umin,Ulast)
    Uaft=min(Uaft,Umax)
    Naft=log(Umax/Uaft)/lgdu
    Naft=max(Naft,0)
    !
    Ubef=min(Umax,Ufirst)
    Ubef=min(Ubef,Ulast)
    Nbef=log(Ubef/Umin)/lgdu
    Nbef=max(Nbef,0)
    !
    n1=log(Umin/Ufirst)/lgdu
    n1=n1+1
    n2=log(Umax/Ufirst)/lgdu
    n2=n2+1
    n1=max(n1,0)
    n2=max(n1,min(n2,nbOut))
    IF(dbg) THEN
       write(*,*)'Umin,Umax,lgdu,Te,Ufirst,Ulast=',Umin,Umax,lgdu,Te,Ufirst,Ulast
       write(*,*)'nBef,nAft=',nBef,nAft
       write(*,*)'uBef,uAft=',uBef,uAft
       write(*,*)'n1,n2 {Eout}=',n1,n2,' nbOut=',nbOut
       if(n1.eq.0) then
          write(*,*)' uOut(1)=',Eout(1)/te
       elseif(n1.eq.1) then
          write(*,*)' uOut(n1:)=',Eout(n1)/te,Eout(n1+1)/te
       else
          write(*,*)' uOut(:n1:)=',Eout(n1-1)/te,Eout(n1)/te &
               ,Eout(n1+1)/te
       end if
       if(n2==0) then
          write(*,*)' uOut(1)=',Eout(1)/te
       elseif(n2.eq.1) then
          write(*,*)' uOut(n2:)=',Eout(n2)/te,Eout(n2+1)/te
       else
          write(*,*)' uOut(:n2:)=',Eout(n2-1)/te,Eout(n2)/te &
               ,Eout(n2+1)/te
       end if
    ENDIF
    ! 
    if(nBef.eq.0 .and. nAft.eq.0) then
       n1=ctrb_to1(n1+1)
       n3=ctrb_to2(n2)				! 111004
       if(n3.eq.0) n3=ctrb_to2(max(1,n2-1))		! 111004
       n2=n3						! 111004
       if(dbg) write(*,*)'420 : n1,n2=',n1,n2,'   n3=',n3
    elseif(nBef.eq.0)then
       n1=ctrb_to1(n1+1)
       n1=max(n1,1)		! 110827
       n2=nbContrib
    elseif(nAft.eq.0) then
       n1=1
       n2=ctrb_to2(n2)
    else
       n1=1
       n2=nbContrib
    end if
    !	
    IF(dbg) THEN
       write(*,*)'Ufirst,Ulast=',Ufirst,Ulast
       write(*,*) &
            ' Ubef,Uaft=',Ubef,Uaft
       write(*,*) &
            ' lgdu,nBef,nAft=',REAL(lgdu),nBef,nAft
       write(*,*)'n1,n2 {ctrb}=',n1,n2
       write(*,*)'fr,to1=',fr_ctrb(n1),to_ctrb(n1) &
            ,' c=',coef_ctrb(n1)
       write(*,*)'fr,to2=',fr_ctrb(n2),to_ctrb(n2) &
            ,' c=',coef_ctrb(n2)
    ENDIF
    !c
    !c     Umin    Umax
    if(Ulast.le.Umin) then	! very high temp.
       !c	|	|
       !c  |--|			   {Emin:Emax}/Te
       write(*,*)'  --- |  |'
    elseif(Ufirst.gt.Umax) then	! very low temp, go to LTE ?
       !c	|	|
       !c 		   |--|   {Emin:Emax}/Te

       gotoLTE=.true.
       write(*,*)' |  |  ---'
    elseif(Ufirst.ge.Umin .and. Ulast.le.Umax) then
       !c	|	|
       !c 	  |--|  	  {Emin:Emax}/Te : nBef>0, nAft>0
       if(dbg)write(*,*)'  | --- |'
    elseif(Ufirst.le.Umin .and. Ulast.ge.Umax) then
       !c	|	|
       !c     |-----------|	  {Emin:Emax}/Te
       if(dbg)write(*,*)' --|---|--'
    elseif(Ufirst.le.Umin) then
       !c	|	|
       !c     |-----|	     {Emin:Emax}/Te
       if(dbg)write(*,*)' --|-  |  '
    elseif(Ulast.ge.Umax) then
       !c	|	|
       !c            |-----|	  {Emin:Emax}/Te
       if(dbg)write(*,*)'  |  -|--'
    else
       write(*,*)'-P- can not distribute U over Umin..Umax'
       call CON_stop('-P- 473')
    end if
    ! 
    IF(dbg)THEN
       write(*,*)'_______ projSP_E, w/ Te=',te
       write(*,*)'Ef,El=',Efirst,Elast,' nbO=',nbOut
       write(*,*)'Ufirst,Umin=',Ufirst,Umin,' nBef,n1=',nBef,n1
       write(*,*)'Umax,Ulast=',Umax,Ulast,' nAft,n2=',nAft,n2
       write(*,*)'nbIn,nbOut=',nbIn,nbOut,'  n2=',n2
    ENDIF	! (dbg) THEN
    ! 
    ! 
    nOut=0
    if(nBef.ne.0) then
       u=Ubef
       r=SPin(1)
       if(nBef.gt.mxout) goto 1101
       do n=nBef,1,-1
          u=u/rdu
          SPout(n)=r
          Uout(n)=u
       end do
       nOut=nOut+nBef
       if(nout.gt.mxout) goto 1100
    else
       Uout(1)=u
       r=SPin(1)			! 111003
       SPout(1)=r			! 111003
       if(dbg)write(*,*)'506. Uout(1)=',Uout(1),' SPout(1:2)=',SPout(1:2)
    end if
    ! 
    do n=nOut+1,min(mxOut,nbOut)
       SPout(n)=zero
    end do
    u=Ubef
    Uout(nOut+1)=u
    if(n2.ge.n1) then
       tt1=n1
       tt2=to_ctrb(n1)
       do ctr=tt1,n2
          to=to_ctrb(ctr)
          u=Eout(to)/te
          if(u.gt.Umin .or. to.ne.tt2) goto 201
          n1=ctr+1
           if(dbg)write(*,*)'563 : n1,n2=',n1,n2
       end do
201    tt1=to_ctrb(n2)
       tt2=n2
       do ctr=tt2,n1,-1
          to=to_ctrb(ctr)
          u=Eout(to)/te
          if(u.gt.Umin .or. to.ne.tt1) goto 202
          n2=ctr
       end do
       ! 
202    if(nBef.eq.0) then
          to=to_ctrb(n1)
          Uout(1)=Eout(to)/te
       end if
       tt1=nOut+1-to_ctrb(n1)
       do ctr=n1,n2
          fr=fr_ctrb(ctr)
          to=to_ctrb(ctr)
          u=Eout(to+1)/te
          to=to+tt1
          if(to.ge.size(Uout)) then
             write(*,*)'to,size(Uout)=',to,size(Uout)
             call CON_stop('to > size(Uout)')
          end if
          if(to.lt.size(Uout)) &
               Uout(to+1)=u
          c=coef_ctrb(ctr)
          SPout(to)=SPout(to)+c*SPin(fr)
       end do
       nOut=to
    end if
    ! 

    !! ?? que doit-on faire quand  nOut=0 ???

    if(nAft.ne.0) then
       do n=1,nAft
          if(nOut.ge.nbOut) exit		! 110806
          u=u*rdu
          nOut=nOut+1
          SPout(nOut)=zero
          Uout(nOut+1)=u
       end do
    end if
    ! 
    return
    ! 
110 nOut=min(mxOut,int(nbUout))
    u=Umin
    do n=1,nOut
       SPout(n)=r
       Uout(n)=u
       u=u*rdu
    end do
    Uout(nOut+1)=u
    ! 
    return
1100 write(*,*)'-R- nOut=',nOut,'  > mxOut=',mxOUt
    call CON_stop('- projSP_E, error=1100 -')
1101 write(*,*)'-R- nbef=',nbef,'  > mxOut=',mxOUt
    call CON_stop('- projSP_E, error=1101 -')
  end subroutine projSP_E
  !-------
  subroutine projSP(Uin,SPin,nbIn)	
    !
    !  Uin(i)=eg(i)/Te, for i=1,nbIn+1
    !
    implicit none
    integer, intent(in) :: nbIn
    integer             :: out1,out2	
    real,intent(in) :: Uin(nbIn+1),SPin(nbIn)	
    real :: du,u1,u2,r
    !
    integer :: nIN,nOut
    real,parameter :: zero=0,one=1
    !
    do nOut=1,nbOut
       SPout(nOut)=0
    end do
    nOut=1
    du=Uout(nOut+1)-Uout(nOut)
    out1=0
    out2=0
    LOOP10:	do nIn=1,nbIn	
1      if(Uin(nIn).ge.Uout(nOut+1)) then	! input segment does not
          if(nOut.ge.nbOut) exit LOOP10 !  goto 20		!  overlap actual output segment
          nOut=nOut+1
          du=Uout(nOut+1)-Uout(nOut)
          goto 1
       elseif(Uin(nIn+1).le.Uout(nOut)) then
          cycle LOOP10 
       end if
2      u1=max(Uin(nIn),Uout(nOut))
       u2=min(Uin(nIn+1),Uout(nOut+1))
       SPout(nOut)=SPout(nOut) + (u2-u1)/du * SPin(nIn)
       !   Uin :	  |      |
       !   Uout:	    |..a   b
       !   u1,u2           |  A B
       if(out1.eq.0 .and. Uout(nOut).ge.u1) out1=nOut
       !   Uin :	    |      |
       !   Uout:	  a   b ..|
       !   u1,u2           A B   |
       if(u2.eq.Uout(nOut+1)) out2=nOut 
       if(u2.ne.Uin(nIn+1))then	! input segment overlaps 
          if(nOut.ge.nbOut) exit LOOP10	! more than 1 output segment
          nOut=nOut+1
          du=Uout(nOut+1)-Uout(nOut)
          goto 2
       end if
    end do LOOP10
    !
    IF(out1.ne.0 .and. out1.lt.nOut) THEN		! 080831
       if(SPout(out1).gt.0 .and. SPout(out1+1).gt.0) then
          r=( log(SPout(out1)) -log(SPout(out1+1)) ) &
               /( log( Uout(out1)) -log( Uout(out1+1)) )
          r=min(zero,r)
          do nOut=out1-1,1,-1
             SPout(nOut)=SPout(nOut+1)*exp(r*  &
                  (log(Uout(nOut))-log(Uout(nOut+1))) )
          end do
       end if
    ENDIF						! 080831
    ! 
    IF(out2.gt.1 .and. out2.le.nOut) THEN		! 080831
       if(SPout(out2).gt.0 .and. SPout(out2-1).gt.0) then
          r=( log(SPout(out2)) -log(SPout(out2-1)) ) &
               /( log( Uout(out2)) -log( Uout(out2-1)) )
          r=max(zero,r)
          do nOut=out2+1,nbOut
             SPout(nOut)=SPout(nOut-1)*exp(r*  &
                  (log(Uout(nOut))-log(Uout(nOut-1))) )
          end do
       end if
    ENDIF						! 080831
    ! 
    return
  end subroutine projSP
  !====================
  subroutine wr5(nom,val,nb)
    character(LEN=*),intent(in)::nom
    integer,intent(in)::nb
    real,intent(in)::val(nb)

    write(*,'(a,/,1p,(5e12.3))') nom,val(1:nb)

  end subroutine wr5
  !===================
  ! \
end MODULE CRASH_M_projE
! /
!------- radiom.f90

