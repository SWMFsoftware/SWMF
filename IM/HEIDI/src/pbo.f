c 02/16/2001 2:22:28 PM

cccccccccccccccccccc
ccc file pbo.for ccc
cccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccc
ccc subroutine initmain, plasmasphere ccc
ccccccccccccccccccccccccccccccccccccccccc

      subroutine initmain()

c Set nthetacells, nphicells array index
      integer nthetacells, nphicells
ccc Use this line for Lmax=15
c      parameter (nthetacells = 100, nphicells = 120)
ccc Use this line for Lmax=10
      parameter (nthetacells = 60, nphicells = 120)
c Set thetamin, thetamax in Degrees
      real thetamin, thetamax
ccc Use this line for Lmax=15
c      parameter (thetamin = 14.963217, thetamax = 60.0)
ccc Use this line for Lmax=10
      parameter (thetamin = 18.434949, thetamax = 60.0)

c Input for entry getgrid
c Input: nt, np size of thetagrid, phigrid arrays
c Output: vthetacells, vphicells are put into thetagrid, phigrid

      integer nt, np
      real thetagrid(nt),phigrid(np)

c inputs for entry setfluxtubevol
c Input: nt, np size of thetagrid, phigrid arrays
c input: fluxtubevol of sixe nt, np

      real fluxtubevol(nt,np)

c inputs for entry setxygrid
c Input: nt, np size of thetagrid, phigrid arrays
c input: gridx,gridy,gridoc of sixe nt, np

      real gridx(nt,np)
      real gridy(nt,np)
      real gridoc(nt,np)

c inputs for getdensity
c Input: nt, np size of thetagrid, phigrid arrays
c input: density of sixe nt, np

      real density(nt,np)

c inputs for getpot
c Input: nt, np size of thetagrid, phigrid arrays
c input: pot of sixe nt, np
c input: wc=0 corotating; wc=1 inertial 

      real pot(nt,np)
      integer wc

c Inputs for entry plasmasphere
c Input: Delt in second
      real delt
c Input: par
      real par(2)

c Inputs for entry saveplasmasphere
c Input: filename
      character filename*80

c Internal variables

c delr in Re, delphi in degrees
      real delr, delphi
c vrcells in Re, vthetacells, vphicells in degrees 
      real vrcells(nthetacells), vthetacells(nthetacells)
      real vphicells(nphicells)
c mgridb in tesla
      real mgridb(nthetacells,nphicells)
c mgridbi in tesla
      real mgridbi(nthetacells,nphicells)
c mgridpot in volts
      real mgridpot(nthetacells,nphicells)
c mgrider, mgridep in volts/meter
      real mgrider(nthetacells,nphicells)
      real mgridep(nthetacells,nphicells)
c mgridvr in meter/sec, and mgridvp in degree /sec
      real mgridvr(nthetacells,nphicells)
      real mgridvp(nthetacells,nphicells)
c mgridn in particles / weber
      real mgridn(nthetacells,nphicells)
c mgridhalf in particles / weber (Work space for upwind)
      real mgridhalf(nthetacells,nphicells)
c mgridden in particles / m**3
      real mgridden(nthetacells,nphicells)
c mgridvol in m**3 / weber
      real mgridvol(nthetacells,nphicells)
c mgridx, mgridy, in Re
      real mgridx(nthetacells,nphicells)
      real mgridy(nthetacells,nphicells)
c mgridoc, open(0) or closed(1) table
      real mgridoc(nthetacells,nphicells)

      real pari(2)
      data pari /-1,-1/

      real pi, rad, re
      integer i, j
      real rmin, rmax
      integer nrcells
      real maxvr, maxvp, deltr, deltp, deltmax, time

c variables for getgrid
      real gdelr, gvrcell

      save nrcells,rmin,rmax,vrcells,vthetacells,
     *   vphicells,delr,delphi,mgridb,mgridbi,
     *   mgridpot,mgrider,mgridep,mgridvr,mgridvp,
     *   mgridn,mgridden,mgridvol,mgridx,mgridy,
     *   mgridoc,pari,deltmax

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters

      nrcells = nthetacells
      rmin = 1.0/(sin(thetamax*rad)*sin(thetamax*rad))
      rmax = 1.0/(sin(thetamin*rad)*sin(thetamin*rad))

      print*,'rmin = ',rmin
      print*,'rmax = ',rmax

      delr = ((rmax - rmin) / float(nrcells-1))
      delphi  = (360.0 / float(nphicells))

      print*,'delr = ',delr
      print*,'delphi = ',delphi

      do i = 1, nrcells
       vrcells(i) = rmin + (float(i-1) * delr)
       vthetacells(i) = (asin(sqrt(1.0/vrcells(i))))/rad
      enddo

      do j = 1, nphicells
       vphicells(j) = (j-1)*360.0/nphicells
      enddo

      print*, 'Number of middle grid cells = ',nrcells*nphicells

      print*, 'Getting equatorial B field on middle grid'
      call getmgridb(nrcells,nphicells,vrcells,mgridb)

      print*, 'Getting ionospheric B field on middle grid'
      call getmgridbi(nrcells,nphicells,vthetacells,mgridbi)

c get flux tube volumes
       print*, 'Getting volume of flux tubes on middle grid'
       call getdipolevol(nthetacells,nphicells,vthetacells,mgridvol)

c get equatorial locations of flux tubes
      print*, 'Getting x, y values for flux tubes'
      call getxydipole(nthetacells,nphicells,vthetacells,
     *   vphicells,mgridx,mgridy,mgridoc)

c set initial particle distribution
      print*, 'Setting initial content of flux tubes on middle grid'
      call initmgridn(nrcells,nphicells,vrcells,mgridn,mgridden,
     *   mgridvol,mgridoc)

      return

ccccccccccccccccccccc
ccc entry getgrid ccc
ccccccccccccccccccccc

      entry getgrid(thetagrid,nt,phigrid,np)

      gdelr = ((rmax - rmin) / (float(nt-1)))
      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree

      do i = 1, nt
       gvrcell = rmin + (float(i-1) * gdelr)
       thetagrid(i) = (asin(sqrt(1.0/gvrcell)))/rad
      enddo

      do j = 1, np
       phigrid(j) = (j-1)*360.0/np
      enddo

      return

cccccccccccccccccccccccccccc
ccc entry setfluxtubevol ccc
cccccccccccccccccccccccccccc

      entry setfluxtubevol(thetagrid,nt,phigrid,np,fluxtubevol)

      call interpol2dpolar(thetagrid,nt,phigrid,np,fluxtubevol,
     *   vthetacells,nthetacells,vphicells,nphicells,mgridvol)

      return

ccccccccccccccccccccccc
ccc entry setxygrid ccc
ccccccccccccccccccccccc

      entry setxygrid(thetagrid,nt,phigrid,np,gridx,gridy,gridoc)

      call interpol2dpolar(thetagrid,nt,phigrid,np,gridx,
     *   vthetacells,nthetacells,vphicells,nphicells,mgridx)

      call interpol2dpolar(thetagrid,nt,phigrid,np,gridy,
     *   vthetacells,nthetacells,vphicells,nphicells,mgridy)

      call interpol2dpolar(thetagrid,nt,phigrid,np,gridoc,
     *   vthetacells,nthetacells,vphicells,nphicells,mgridoc)

      return

ccccccccccccccccccccccccc
ccc entry initdensity ccc
ccccccccccccccccccccccccc

      entry initdensity()

c set initial particle distribution
      print*, 'Setting initial content of flux tubes on middle grid'
      call initmgridn(nrcells,nphicells,vrcells,mgridn,mgridden,
     *   mgridvol,mgridoc)

      return

cccccccccccccccccccccccc
ccc entry getdensity ccc
cccccccccccccccccccccccc

      entry getdensity(thetagrid,nt,phigrid,np,density)

      call interpol2dpolar(vthetacells,nthetacells,vphicells,
     *   nphicells,mgridden,thetagrid,nt,phigrid,np,density)

      return

cccccccccccccccccccc
ccc entry setpot ccc
cccccccccccccccccccc

      entry setpot(thetagrid,nt,phigrid,np,pot)

      call interpol2dpolar(thetagrid,nt,phigrid,np,pot,
     *   vthetacells,nthetacells,vphicells,nphicells,mgridpot)

c      print*, 'Differencing setpot on middle grid to get electric ',
c     *   'field'
      call gradpot(nrcells,nphicells,vrcells,delr,delphi,mgridpot,
     *   mgrider,mgridep)

      maxvr = 0.0
      maxvp = 0.0

c      print*, 'Calculating the E cross B drift velocity on middle grid'
      call ecrossb(nrcells,nphicells,vrcells,mgrider,mgridep,mgridb,
     *   mgridvr,mgridvp,maxvr,maxvp)

c deltr,deltp in units of seconds
      re    = 6.378E6   ! in meters
      deltr = (delr*re)/maxvr
      deltp = delphi/maxvp

c      print*, 'Maximum radial velocity = ',maxvr,' meters/sec'
c      print*, 'Maximum radial time step = ',deltr,' seconds'
c      print*, 'Maximum azimuthal velocity = ',maxvp,' degrees/sec'
c      print*, 'Maximum azimuthal time step = ',deltp,' seconds'

c calculate maximum time step to use
c delt in units of seconds
      deltmax = 10.0*aint(amin1(deltp,deltr)/10.0)
      if (deltmax.eq.0.0) deltmax = aint(amin1(deltp,deltr))
      if (deltmax.eq.0.0) then
       print*,'Use a coarser grid.'
       print*, 'Maximum radial velocity = ',maxvr,' meters/sec'
       print*, 'Maximum radial time step = ',deltr,' seconds'
       print*, 'Maximum azimuthal velocity = ',maxvp,' degrees/sec'
       print*, 'Maximum azimuthal time step = ',deltp,' seconds'
       print*, 'deltmax,deltp,deltr = ',deltmax,deltp,deltr
       call CON_stop('ERROR in pbo.f')
      end if

      return

cccccccccccccccccccc
ccc entry getpot ccc
cccccccccccccccccccc

      entry getpot(thetagrid,nt,phigrid,np,pot,par,wc)

c      print*, 'Getting electric potential on middle grid'
      call getmgridpot(nrcells,nphicells,vrcells,vthetacells,
     *   vphicells,mgridpot,par)

      if (wc.eq.1) then
       call addcorotpot(nrcells,nphicells,vrcells,vthetacells,
     *    vphicells,mgridpot)
      endif

c      print*,'Got pot'

      call interpol2dpolar(vthetacells,nthetacells,vphicells,
     *   nphicells,mgridpot,thetagrid,nt,phigrid,np,pot)

c      print*,'finished'

      return

cccccccccccccccccccccccccc
ccc entry plasmasphere ccc
cccccccccccccccccccccccccc

      entry plasmasphere(delt,par)

      if ((par(1).ne.pari(1)).or.(par(2).ne.pari(2))) then 

c       print*, 'Getting electric potential on middle grid'
       call getmgridpot(nrcells,nphicells,vrcells,vthetacells,
     *    vphicells,mgridpot,par)

c       print*, 'Differencing potential on middle grid to get electric ',
c     *    'field'
       call gradpot(nrcells,nphicells,vrcells,delr,delphi,mgridpot,
     *    mgrider,mgridep)

       maxvr = 0.0
       maxvp = 0.0

c       print*, 'Calculating the E cross B drift velocity on middle grid'
       call ecrossb(nrcells,nphicells,vrcells,mgrider,mgridep,mgridb,
     *    mgridvr,mgridvp,maxvr,maxvp)

c deltr,deltp in units of seconds
       deltr = (delr*re)/maxvr
       deltp = delphi/maxvp

c       print*, 'Maximum radial velocity = ',maxvr,' meters/sec'
c       print*, 'Maximum radial time step = ',deltr,' seconds'
c       print*, 'Maximum azimuthal velocity = ',maxvp,' degrees/sec'
c       print*, 'Maximum azimuthal time step = ',deltp,' seconds'

c calculate maximum time step to use
c delt in units of seconds
       deltmax = 10.0*aint(amin1(deltp,deltr)/10.0)
       if (deltmax.eq.0.0) deltmax = aint(amin1(deltp,deltr))
       if (deltmax.eq.0.0) then
        print*,'Use a coarser grid.'
        print*, 'Maximum radial velocity = ',maxvr,' meters/sec'
        print*, 'Maximum radial time step = ',deltr,' seconds'
        print*, 'Maximum azimuthal velocity = ',maxvp,' degrees/sec'
        print*, 'Maximum azimuthal time step = ',deltp,' seconds'
        call CON_stop('ERROR in pbo.f')
       end if

       pari(1) = par(1)
       pari(2) = par(2)

      end if

c      print*, 'Upwind differencing on middle grid to advance ',
c     *   'solution in time'

      if (deltmax.ge.delt) deltmax = delt

c      print*, 'Using a time step of ',deltmax,' seconds'

      time = 0.0

      do while (time.lt.delt)

      if ((delt-time).lt.deltmax) deltmax = (delt-time)

c      print*,'Time = ',time
c      print*,'Time step = ',deltmax

c Upwind differencing on middle grid to advance solution in time
c       call upwind(nrcells,nphicells,vrcells,delr,delphi,mgridn,
c     *    mgridvr,mgridvp,deltmax,mgridhalf)

c Upwind/LaxWendroff-superbee differencing on middle grid to 
c advance solution in time
       call superbee(nrcells,nphicells,vrcells,delr,delphi,mgridn,
     *    mgridvr,mgridvp,deltmax,mgridhalf)

c calculate filling and draining of flux tubes
       call filling(nrcells,nphicells,vrcells,vthetacells,vphicells,
     *    mgridn,mgridden,mgridvol,mgridoc,mgridbi,deltmax)

       time = time + deltmax

      end do

c      print*,'Time = ',time
c      print*,'Time step = ',deltmax

c      print*,'finished'

      return

cccccccccccccccccccccccccccccc
ccc entry saveplasmasphere ccc
cccccccccccccccccccccccccccccc

      entry saveplasmasphere(filename)

      call saveit(vthetacells,nthetacells,vphicells,nphicells,
     *   mgridden,mgridx,mgridy,mgridoc,mgridpot,filename)

      return

cccccccccccccccccccccccccccccc
ccc entry loadplasmasphere ccc
cccccccccccccccccccccccccccccc

      entry loadplasmasphere(filename)

      call loadit(vthetacells,nthetacells,vphicells,nphicells,
     *   mgridden,mgridx,mgridy,mgridoc,filename)

      call denton(nthetacells,nphicells,mgridden,mgridvol,
     *    mgridoc,mgridn)

      return
      end

ccccccccccccccccccccccccc
ccc subroutine saveit ccc
ccccccccccccccccccccccccc

      subroutine saveit(vthetacells,nthetacells,vphicells,nphicells,
     *   mgridden,mgridx,mgridy,mgridoc,mgridpot,filename)

      use ModIoUnit, ONLY : UNITTMP_
      
c Input: nthetacells, nphicells array index
      integer nthetacells, nphicells
c Input: vthetacells, vphicells in degrees 
      real vthetacells(nthetacells), vphicells(nphicells)
c Input: mgridden in particles / m**3
      real mgridden(nthetacells,nphicells)
c Input: filename
      character filename*80
c mgridx, mgridy, in Re
      real mgridx(nthetacells,nphicells)
      real mgridy(nthetacells,nphicells)
c mgridoc, open(0) or closed(1) table
      real mgridoc(nthetacells,nphicells)
c mgridpot in volts
      real mgridpot(nthetacells,nphicells)

            
      open(unit = UNITTMP_, file=filename, status = 'unknown',
     *   form = 'formatted')
	print *, 'Saving the DGCPM output'
      write(UNITTMP_,*) nthetacells, nphicells
      write(UNITTMP_,*) vthetacells
      write(UNITTMP_,*) vphicells
      write(UNITTMP_,*) mgridden
      write(UNITTMP_,*) mgridx
      write(UNITTMP_,*) mgridy
      write(UNITTMP_,*) mgridoc
      write(UNITTMP_,*) mgridpot
      close(unit = UNITTMP_)

      return
      end

ccccccccccccccccccccccccc
ccc subroutine loadit ccc
ccccccccccccccccccccccccc

      subroutine loadit(vthetacells,nthetacells,vphicells,nphicells,
     *   mgridden,mgridx,mgridy,mgridoc,filename)

      use ModIoUnit, ONLY : io_unit_new


c Input: filename
      character filename*80

c Output: nthetacells, nphicells array index
      integer nthetacells, nphicells
c Output: vthetacells, vphicells in degrees 
      real vthetacells(nthetacells), vphicells(nphicells)
c Output: mgridden in particles / m**3
      real mgridden(nthetacells,nphicells)
c Output: mgridx, mgridy, in Re
      real mgridx(nthetacells,nphicells)
      real mgridy(nthetacells,nphicells)
c Output: mgridoc, open(0) or closed(1) table
      real mgridoc(nthetacells,nphicells)

c Internal: nthetacells1, nphicells1 array index
      integer nthetacells1, nphicells1
      integer nthetacells2, nphicells2
c      parameter (nthetacells2 = 200, nphicells2 = 720)
ccc Set the next line to match the input file you are using:
      parameter (nthetacells2 = 60, nphicells2 = 120)
c Internal: vthetacells1, vphicells1 in degrees 
      real vthetacells1(nthetacells2), vphicells1(nphicells2)
c Internal: mgridden1 in particles / m**3
      real mgridden1(nthetacells2,nphicells2)
c Internal: mgridx1, mgridy1, in Re
      real mgridx1(nthetacells2,nphicells2)
      real mgridy1(nthetacells2,nphicells2)
c Internal: mgridoc1, open(0) or closed(1) table
      real mgridoc1(nthetacells2,nphicells2)

      integer :: iUnitIn != 33
      
      iUnitIn = io_unit_new() 
      
      open(unit =iUnitIn, file=filename, status = 'old',
     *   form = 'formatted')
      read(iUnitIn,*) nthetacells1, nphicells1

      if (nthetacells1.ne.nthetacells2.or.nphicells1.ne.nphicells2) then
       print*,'File size mismatch in subroutine Loadit'
       call CON_stop('ERROR in pbo.f')
      endif

      read(iUnitIn,*) vthetacells1
      read(iUnitIn,*) vphicells1
      read(iUnitIn,*) mgridden1
      read(iUnitIn,*) mgridx1
      read(iUnitIn,*) mgridy1
      read(iUnitIn,*) mgridoc1
      close(unit = iUnitIn)

      call interpol2dpolar(vthetacells1,nthetacells1,vphicells1,
     *   nphicells1,mgridden1,vthetacells,nthetacells,vphicells,
     *   nphicells,mgridden)

      call interpol2dpolar(vthetacells1,nthetacells1,vphicells1,
     *   nphicells1,mgridx1,vthetacells,nthetacells,vphicells,
     *   nphicells,mgridx)

      call interpol2dpolar(vthetacells1,nthetacells1,vphicells1,
     *   nphicells1,mgridy1,vthetacells,nthetacells,vphicells,
     *   nphicells,mgridy)

      call interpol2dpolar(vthetacells1,nthetacells1,vphicells1,
     *   nphicells1,mgridoc1,vthetacells,nthetacells,vphicells,
     *   nphicells,mgridoc)

      return
      end

ccccccccccccccccccccccccccc
ccc function saturation ccc
ccccccccccccccccccccccccccc

      real function saturation(l)

c Carpenter and Anderson's saturation density in units of
c particles / m**3 
c Carpenter and Anderson, JGR, p. 1097, 1992.

c input: l in re
      real l

c output: saturation in particles / m**3

      saturation = (1.0e6) * 10.0**((-0.3145*l)+3.9043)

      return
      end

ccccccccccccccccccccccc
ccc function trough ccc
ccccccccccccccccccccccc

      real function trough(l)

c Carpenter and Anderson's trough density in units of
c particles / m**3 
c Carpenter and Anderson, JGR, p. 1097, 1992.

c input: l in re
      real l

c output: trough in particles / m**3

      trough = (1.0e6) * 0.5 * ((10.0/l)**4.0)

      return
      end

cccccccccccccccccccccccccccccccccc
ccc function DipoleFluxTubeVol ccc
cccccccccccccccccccccccccccccccccc

      real function dipoleFluxTubeVol(l)

c calculates the unit volume of a dipole magnetic field flux tube
c (the volume in m**3 per unit of magnetic flux(weber))

c 1 tesla = newton/(ampere-meter) or (volt-sec)/meter**2
c 1 weber = Tesla-m**2 or joule/ampere or volt-sec

c input: l in re
      real l

c output: dipoleFluxTubevol in m/Tesla or m**3/weber

      real pi, re, mu, m

      pi = 3.14159          ! rad
      re = 6.378e6          ! radius of Earth in meters
      mu = 4.0*pi*1.0e-7    ! newtons/amps**2
      m = 8.05e22           ! amps*meter**2

      dipoleFluxTubeVol = ((4.0*pi)/(mu*m)) * (32.0/35.0) * (l**4) * 
     *    sqrt(1.0-(1.0/l)) * (1.0+(1.0/(2.0*l))+(3.0/(8.0*l*l))+
     *    (5.0/(16.0*l*l*l))) * (re**4.0)

      return
      end

ccccccccccccccccccccccccccc
ccc subroutine mydipole ccc
ccccccccccccccccccccccccccc

      subroutine mydipole(r,theta,br,btheta)

c calculates the two components of a dipole magnetic field

c input: r in re, theta and phi in degrees
      real r, theta

c output: br,btheta,bphi in tesla
c 1 tesla = 1 newton/(ampere-meter)
      real br, btheta

      real pi, rad, re, mu, m
      real thetarad, mum

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters
      mu = 4.0*pi*1.0e-7    ! newtons/amps**2
      m = 8.05e22           ! amps*meter**2

      thetarad = theta * rad
      mum = (mu*m) / (2.0*pi*r*r*r*re*re*re)

      br = - mum * cos(thetarad)
      btheta = - (mum*sin(thetarad)) / 2.0

      return
      end

cccccccccccccccccccccccccc
ccc subroutine dipoleb ccc
cccccccccccccccccccccccccc

      subroutine dipoleb(r,theta,eb)

c calculates the magnitude of a dipole magnetic field
c 1 tesla = 1 newton/(ampere-meter)

c input: r in re, theta in degrees
      real r, theta

c output: eb in tesla
      real eb

      real pi, rad, re, mu, m
      real cost

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters
      mu = 4.0*pi*1.0e-7    ! newtons/amps**2
      m = 8.05e22           ! amps*meter**2

      cost = cos(theta*rad)
      eb = ((mu*m)/(4.0*pi*r*r*r*re*re*re))*sqrt(1.0+(3.0*cost*cost))

      return
      end

ccccccccccccccccccccccccccccccc
ccc subroutine DipoleLshell ccc
ccccccccccccccccccccccccccccccc

      subroutine dipoleLshell(r,theta,l)

c calculates the L parameter of a dipole magnetic field line

c input: r in re, and theta in degress
      real r, theta

c output: l parameter in re
      real l

      real pi, rad
      real sint

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree

      sint = sin(theta*rad)
      l = r / (sint*sint)

      return
      end

ccccccccccccccccccccccc
ccc subroutine coro ccc
ccccccccccccccccccccccc

      subroutine coro(r,theta,vphi)

c calculates the corotation velocity in meter / sec

c input: r in Re, theta in degrees
      real r, theta

c output: vphi meters/second
      real vphi

c 1 volt = 1 joule/coulomb
c 1 coulomb = 1 amp sec

      real pi, rad, re, w

      pi = 3.14159                    ! rad
      rad = pi / 180.0                ! rad/degree
      re = 6.378e6                    ! radius of Earth in meters
      w = (2.0*pi) / (24.0*3600.0)    ! rad/sec

      vphi = w*r*re*sin(theta*rad)    ! meters/sec

      return
      end

ccccccccccccccccccccccccc
ccc subroutine dipsph ccc
ccccccccccccccccccccccccc

      subroutine dipsph(s,q,r,theta,j)

c converts dipole coords into spherical ones and vica versa
c (theta in degrees).

c              j>0            j<0
c input:       s,q           r,theta
c output:     r,teta          s,q

      real s, q, r, theta
      integer j

      real pi, rad
      real f, err, step, g, ft, thetar, sint

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree

      if (j.gt.0) then
       f = s / (q*q)
       err = 1.0
       step = 22.5 * rad
       g = 45.0 * rad
       do while (abs(err).gt.1e-6)
        ft = cos(g) / ((sin(g))**4)
        err = (f-ft) / f
        if (err.gt.0.0) then
         g = g - step
        else
         g = g + step
        end if
        step = step / 2.0
        if (step.eq.0.0) then
         print*,'dipsph failed'
         print*,f,ft,err,g/rad,step/rad
         call CON_stop('ERROR in pbo.f')
        end if
       end do
       theta = g / rad
       sint = sin(theta*rad)
       r = sint * sint / q
      else
       thetar = theta * rad
       sint = sin(thetar)
       s = cos(thetar) / (r*r)
       q = (sint*sint) / (r)
      end if

      return
      end

ccccccccccccccccccccccccc
ccc subroutine vpocar ccc
ccccccccccccccccccccccccc

      subroutine vpocar(theta,vr,vtheta,vx,vy)

c Calculates cartesian vector components from polar

c Input: theta in degrees, and vr,vtheta vector components
      real theta, vr, vtheta

c Output: vx,vy vector components
      real vx, vy

      real pi, rad
      real tr, sf, cf

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree

      tr = theta * rad
      sf = sin(tr)
      cf = cos(tr)
      vx = vr*cf - vtheta*sf
      vy = vr*sf + vtheta*cf

      return
      end

cccccccccccccccccccccccccccc
ccc subroutine getmgridb ccc
cccccccccccccccccccccccccccc

      subroutine getmgridb(nrcells,nphicells,vrcells,mgridb)

c input: nrcells, nphicells array index
      integer nrcells, nphicells
c input: vrcells in Re
      real vrcells(nrcells)

c output: mgridb in tesla
      real mgridb(nrcells,nphicells)

      integer i, j
      real bfield

      do i = 1, nrcells
       call dipoleb(vrcells(i),90.0,bfield)
       do j = 1, nphicells
        mgridb(i,j) = bfield
       enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccc
ccc subroutine getmgridbi ccc
ccccccccccccccccccccccccccccc

      subroutine getmgridbi(nrcells,nphicells,vthetacells,mgridbi)

c input: nrcells, nphicells array index
      integer nrcells, nphicells
c input: vthetacells in degrees
      real vthetacells(nrcells)

c output: mgridbi in tesla
      real mgridbi(nrcells,nphicells)

      integer i, j
      real bfield

      do i = 1, nrcells
       call dipoleb(1.0,vthetacells(i),bfield)
       do j = 1, nphicells
        mgridbi(i,j) = bfield
       enddo
      enddo

      return
      end

cccccccccccccccccccccccccc
ccc subroutine ecrossb ccc
cccccccccccccccccccccccccc

      subroutine ecrossb(nrcells,nphicells,vrcells,mgrider,mgridep,
     *              mgridb,mgridvr,mgridvp,maxvr,maxvp)

c Input: nrcells, nphicells array index
      integer nrcells, nphicells
c Input: vrcells in Re 
      real vrcells(nrcells)
c Input: mgridb in tesla
      real mgridb(nrcells,nphicells)
c Input: mgrider, mgridep in volts/meter
      real mgrider(nrcells,nphicells), mgridep(nrcells,nphicells)

c Output: mgridvr in meter/sec, and mgridvp in degree /sec
      real mgridvr(nrcells,nphicells), mgridvp(nrcells,nphicells)
c Output: maxvr, maxvp in meter / sec
      real maxvr, maxvp

      real pi, rad, re
      integer i, j, ivr, jvr, ivp, jvp
      real vt, vc

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters

      maxvr = 0.0
      maxvp = 0.0
      ivr=0
      jvr=0
      ivp=0
      jvp=0
      do i = 1, nrcells
       do j = 1, nphicells
        mgridvr(i,j) = mgridep(i, j) / mgridb(i,j)
        if (abs(mgridvr(i,j)).gt.maxvr) then
	  maxvr = abs(mgridvr(i,j))
	  ivr = i
          jvr = j
	endif
        vt = - mgrider(i, j) / mgridb(i,j)
        call coro(vrcells(i),90.0,vc)
        mgridvp(i,j) = (vt + vc)/(vrcells(i)*re*rad)
        if (abs(mgridvp(i,j)).gt.maxvp) then
	  maxvp = abs(mgridvp(i,j))
	  ivp  = i
	  jvp = j
	endif
       enddo
      enddo
	print *, 'Max Vr:',ivr,jvr,maxvr,vrcells(ivr)
	print *, 'Max Vp:',ivp,jvp,maxvp,vrcells(ivp)

      return
      end

ccccccccccccccccccccccccccccccc
ccc subroutine getdipolevol ccc
ccccccccccccccccccccccccccccccc

      subroutine getdipolevol(nthetacells,nphicells,
     *    vthetacells,mgridvol)

c Input: nthetacells, nphicells array index
      integer nthetacells, nphicells
c Input: vthetacells in degrees 
      real vthetacells(nthetacells)

c Output: mgridvol in m**3 / weber
      real mgridvol(nthetacells,nphicells)

      real pi, rad, re
      integer i, j
      real dvol, vrcell, st

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters

      do i = 1, nthetacells
        st = sin(vthetacells(i)*rad)
        vrcell = 1.0/(st*st)
        dvol = dipolefluxtubevol(vrcell)
       do j = 1, nphicells
        mgridvol(i,j) = dvol
       enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccc
ccc subroutine getxydipole ccc
cccccccccccccccccccccccccccccc

      subroutine getxydipole(nthetacells,nphicells,
     *   vthetacells,vphicells,mgridx,mgridy,mgridoc)

c Theta is zero at the north pole, positive towards the equator 
c phi is zero at 24 MLT (antisunward), positive rotation towards dawn
c x is positive towards dusk (phi = 270.0)
c y is positive towards the sun (phi = 180.0)

c Set nthetacells, nphicells array index
      integer nthetacells, nphicells
c vthetacells, vphicells in degrees 
      real vthetacells(nthetacells)
      real vphicells(nphicells)
c mgridx, mgridy, in Re
      real mgridx(nthetacells,nphicells)
      real mgridy(nthetacells,nphicells)
c mgridoc, open(0) or closed(1) table
      real mgridoc(nthetacells,nphicells)

      real pi, rad, angle, ca, sa, st, r
      integer i, j

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree

      do i = 1, nthetacells
       st = sin(vthetacells(i)*rad)
       r = 1.0/(st*st)
       do j = 1, nphicells
        angle = (180.0+vphicells(j))*rad
        ca = cos(angle)
        sa = sin(angle)
        mgridx(i,j) = r * ca
        mgridy(i,j) = r * sa
        mgridoc(i,j) = 1
       enddo
      enddo

      return
      end

ccccccccccccccccccccccccccccc
ccc subroutine initmgridn ccc
ccccccccccccccccccccccccccccc

      subroutine initmgridn(nrcells,nphicells,vrcells,mgridn,mgridden,
     *   mgridvol,mgridoc)

c Input: nrcells, nphicells array index
      integer nrcells, nphicells
c Input: vrcells in Re 
      real vrcells(nrcells)
c Input: mgridvol in m**3 / weber
      real mgridvol(nrcells,nphicells)
c mgridoc, open(0) or closed(1) table
      real mgridoc(nrcells,nphicells)

c Output: mgridn in particles / weber
      real mgridn(nrcells,nphicells)
c Output: mgridden in particles / m**3
      real mgridden(nrcells,nphicells)

c mgridn in units of particles per unit of magnetic flux
c particles/weber or particles/(Tesla-meter**2)
c mgridden in units of particles per m**3

      integer i, j
      real dn

      do i = 1, nrcells-1
       dn = saturation(vrcells(i))
       do j = 1, nphicells
        if (mgridoc(i,j).gt.0.999) then
         mgridden(i,j) = dn
         mgridn(i,j) = mgridden(i,j) * mgridvol(i,j)
        else
         mgridn(i,j) = 100.0
         mgridden(i,j) = 0.0
        end if
       enddo
      enddo
      i = nrcells
      dn = trough(vrcells(i))
      do j = 1, nphicells
       if (mgridoc(i,j).gt.0.999) then
        mgridden(i,j) = dn
        mgridn(i,j) = mgridden(i,j) * mgridvol(i,j)
       else
        mgridn(i,j) = 100.0
        mgridden(i,j) = 0.0
       end if
      enddo

      return
      end

ccccccccccccccccccccccccc
ccc subroutine denton ccc
ccccccccccccccccccccccccc

      subroutine denton(nthetacells,nphicells,mgridden,mgridvol,
     *    mgridoc,mgridn)

c Input: nthetacells, nphicells array index
      integer nthetacells, nphicells
c Input: mgridden in particles / m**3
      real mgridden(nthetacells,nphicells)
c Input: mgridvol in m**3 / weber
      real mgridvol(nthetacells,nphicells)
c Imput: mgridoc, open(0) or closed(1) table
      real mgridoc(nthetacells,nphicells)

c Output: mgridn in particles / weber
      real mgridn(nthetacells,nphicells)

      integer i, j

      do i = 1, nthetacells
       do j = 1, nphicells
        if (mgridoc(i,j).gt.0.999) then
         mgridn(i,j) = mgridden(i,j) * mgridvol(i,j)
        endif
       enddo
      enddo

      return
      end

cccccccccccccccccccccccccc
ccc subroutine filling ccc
cccccccccccccccccccccccccc

      subroutine filling(nrcells,nphicells,vrcells,vthetacells,
     *    vphicells,mgridn,mgridden,mgridvol,mgridoc,mgridbi,delt)

c Input: nrcells, nphicells array index
      integer nrcells, nphicells
c Input: vrcells in Re, vthetacells, vphicells in degrees 
      real vrcells(nrcells), vthetacells(nrcells), vphicells(nphicells)
c Input: delt in seconds
      real delt
c Input: mgridbi in tesla
      real mgridbi(nrcells,nphicells)
c Input: mgridoc, open(0) or closed(1) table
      real mgridoc(nrcells,nphicells)

c Output: mgridn in particles / weber
      real mgridn(nrcells,nphicells)
c Output: mgridden in particles / m**3
      real mgridden(nrcells,nphicells)
c Output: mgridvol in m**3 / weber
      real mgridvol(nrcells,nphicells)

      real pi, rad, re
      integer i, j
      real fmax, dsat, br, f, tden, tn

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters

c fmax is the upward flow of particles in units of 
c particles/m**2/sec
      fmax  = 2.0e12

      do i = 1, nrcells
       dsat = saturation(vrcells(i))
       do j = 1, nphicells
        if (mgridoc(i,j).gt.0.999) then
         if ((vphicells(j).ge.90.0).and.(vphicells(j).le.270.0)) then
          mgridden(i,j) = mgridn(i,j) / mgridvol(i,j)
          tn = mgridn(i,j)
          tden = mgridden(i,j)
          if (tden.lt.dsat) then
           f = ((dsat-mgridden(i,j))/dsat)*fmax
           br = mgridbi(i,j)
           mgridn(i,j) = mgridn(i,j) + ((f*delt)/br)
           mgridden(i,j) = mgridn(i,j) / mgridvol(i,j)
          else
           mgridn(i,j) = dsat * mgridvol(i,j)
           mgridden(i,j) = mgridn(i,j) / mgridvol(i,j)
          endif
         else
          mgridn(i,j) = mgridn(i,j) - 
     *       (mgridn(i,j)*(delt/(10.0*24.0*3600.0)))
          mgridden(i,j) = mgridn(i,j) / mgridvol(i,j)
         end if
         if (mgridden(i,j).le.0.0) then
          print*,'subroutine: filling'
          print*,'i,j,mgridden = ',i,j,mgridden(i,j)
          print*,'mgridoc,mgridvol',mgridoc(i,j),mgridvol(i,j)
          print*,'mgridden',mgridden(i,j)
          print*,'mgridn',mgridn(i,j)
          print*,'tden,br,f',tden,br,f
          print*,'tn',tn
          print*,'vr,dsat',vrcells(i),dsat
          print*,'delt',delt
          print*,'deln',((f*delt)/br)
          call CON_stop('ERROR in pbo.f')
         end if
        else
c         print*,mgridoc(i,j),i,j
         mgridn(i,j) = mgridn(i,j) - mgridn(i,j)*(delt/(24.0*3600.0))
         mgridden(i,j) = 0.0
        end if
       enddo
      enddo

      return
      end

ccccccccccccccccccccccccc
ccc subroutine upwind ccc
ccccccccccccccccccccccccc

      subroutine upwind(nrcells,nphicells,vrcells,delr,delphi,mgridn,
     *   mgridvr,mgridvp,delt,mgridhalf)

C First order upwind differencing 

c Input: nrcells, nphicells array index
      integer nrcells, nphicells
c Input: vrcells in Re
      real vrcells(nrcells)
c Input: delr in Re, delphi in degrees
      real delr, delphi
c Input: delt in seconds
      real delt
c Input: mgridvr in meter/sec, and mgridvp in degree /sec
      real mgridvr(nrcells,nphicells), mgridvp(nrcells,nphicells)
c Input: mgridn in particles / weber
      real mgridn(nrcells,nphicells),mgridhalf(nrcells,nphicells)

C Output: mgridn in particles / weber

      real pi, rad, re
      integer i, j, ip, im, jp, jm
      real delret,delphit
      real small

      small = 0.0001

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters

      delret = delt / (delr * re)

c Advect the radial component first
c For BC on the outer cells, the two rows on the edge 
c are calculated separately outside of the r loop

c do the first r cell
      i = 1
      ip = 2
      do j = 1, nphicells
       if (mgridvr(i,j).lt.-small) then
        mgridhalf(i,j) = mgridn(i,j)  - 
     *     (mgridvr(i,j)*delret) *
     *     (mgridn(ip,j) - mgridn(i,j))
       else
        mgridhalf(i,j) = mgridn(i,j)
       end if
      enddo
c do the middle r cells in a loop
      do i = 2, nrcells - 1
       ip = i + 1
       im = i - 1
       do j = 1, nphicells
        if (abs(mgridvr(i,j)).gt.small) then
         if (mgridvr(i,j).gt.0.0) then
          mgridhalf(i,j) = mgridn(i,j)  - 
     *       (mgridvr(i,j)*delret) *
     *       (mgridn(i,j) - mgridn(im,j))
         else
          mgridhalf(i,j) = mgridn(i,j)  - 
     *       (mgridvr(i,j)*delret) *
     *       (mgridn(ip,j) - mgridn(i,j))
         end if
        else
         mgridhalf(i,j) = mgridn(i,j)
        end if
       enddo
      enddo
c do the last r cell
      i = nrcells
      im = i - 1
      do j = 1, nphicells
       if (mgridvr(i,j).gt.small) then
        mgridhalf(i,j) = mgridn(i,j)  - 
     *     (mgridvr(i,j)*delret) *
     *     (mgridn(i,j) - mgridn(im,j))
       else
        mgridhalf(i,j) = mgridn(i,j)
       end if
      enddo

c Advect the azimuthal component next
c Due to rap around of the phi cells, the two columns on the edge 
c are calculated separately outside of the phi loop

      delphit = delt / delphi

      do i = 1, nrcells
c do the first phi cell
       j = 1
       jp = 2
       jm = nphicells
       if (abs(mgridvp(i,j)).gt.small) then
        if (mgridvp(i,j).gt.0.0) then
         mgridn(i,j) = mgridhalf(i,j)  - 
     *      (mgridvp(i,j)*delphit) * 
     *      (mgridhalf(i,j) - mgridhalf(i,jm))
        else
         mgridn(i,j) = mgridhalf(i,j)  - 
     *      (mgridvp(i,j)*delphit) * 
     *      (mgridhalf(i,jp) - mgridhalf(i,j))
        end if
       else
        mgridn(i,j) = mgridhalf(i,j)
       end if
c Do the rest of the phi cells in a loop
       do j = 2, nphicells - 1
        jp = j + 1
        jm = j - 1
        if (abs(mgridvp(i,j)).gt.small) then
         if (mgridvp(i,j).gt.0.0) then
          mgridn(i,j) = mgridhalf(i,j)  - 
     *       (mgridvp(i,j)*delphit) * 
     *       (mgridhalf(i,j) - mgridhalf(i,jm))
         else
          mgridn(i,j) = mgridhalf(i,j)  - 
     *       (mgridvp(i,j)*delphit) * 
     *       (mgridhalf(i,jp) - mgridhalf(i,j))
         end if
        else
         mgridn(i,j) = mgridhalf(i,j)
        end if
       enddo
c Do the last phi cell
       j = nphicells
       jp = 1
       jm = nphicells - 1
       if (abs(mgridvp(i,j)).gt.small) then
        if (mgridvp(i,j).gt.0.0) then
         mgridn(i,j) = mgridhalf(i,j)  - 
     *      (mgridvp(i,j)*delphit) * 
     *      (mgridhalf(i,j) - mgridhalf(i,jm))
        else
         mgridn(i,j) = mgridhalf(i,j)  - 
     *      (mgridvp(i,j)*delphit) * 
     *      (mgridhalf(i,jp) - mgridhalf(i,j))
        end if
       else
        mgridn(i,j) = mgridhalf(i,j)
       end if
      enddo

      return
      end

ccccccccccccccccccccccccccc
ccc subroutine superbee ccc
ccccccccccccccccccccccccccc

      subroutine superbee(nrcells,nphicells,vrcells,delr,delphi,mgridn,
     *   mgridvr,mgridvp,delt,mgridhalf)

C Mixed first order upwind and
c second order Lax-Wendroff 
c differencing with the 
c Superbee limiter function

c Input: nrcells, nphicells array index
      integer nrcells, nphicells
c Input: vrcells in Re
      real vrcells(nrcells)
c Input: delr in Re, delphi in degrees
      real delr, delphi
c Input: delt in seconds
      real delt
c Input: mgridvr in meter/sec, and mgridvp in degree /sec
      real mgridvr(nrcells,nphicells), mgridvp(nrcells,nphicells)
c Input: mgridn in particles / weber
      real mgridn(nrcells,nphicells),mgridhalf(nrcells,nphicells)

C Output: mgridn in particles / weber

      real pi, rad, re
      integer i, j, ip, ipp, im, imm, jp, jpp, jm, jmm
      real delret,delphit
      real small
      real cou, scou
      real fphup, fmhup, fphlw, fmhlw
      real rjp1, rjp2, rjp, sip
      real rjm1, rjm2, rjm, sim
      real fph, fmh
      real mscale, mdiv

      small = 0.0001

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters
      mscale = 1.E15        ! scaling for mgridn if numbers are too big

      delret = delt / (delr * re)

c Advect the radial component first
c For BC on the outer cells, the two rows on the edge 
c are calculated separately outside of the r loop

c do the first r cell
      i = 1
      ip = 2
      do j = 1, nphicells
       if (mgridvr(i,j).lt.-small) then
        mgridhalf(i,j) = mgridn(i,j)  - 
     *     (mgridvr(i,j)*delret) *
     *     (mgridn(ip,j) - mgridn(i,j))
       else
        mgridhalf(i,j) = mgridn(i,j)
       end if
      enddo

c do the second r cell
      i = 2
      ip = 3
      im = 1
      do j = 1, nphicells
       if (abs(mgridvr(i,j)).gt.small) then
        if (mgridvr(i,j).gt.0.0) then
         mgridhalf(i,j) = mgridn(i,j)  - 
     *      (mgridvr(i,j)*delret) *
     *      (mgridn(i,j) - mgridn(im,j))
        else
         mgridhalf(i,j) = mgridn(i,j)  - 
     *      (mgridvr(i,j)*delret) *
     *      (mgridn(ip,j) - mgridn(i,j))
        end if
       else
        mgridhalf(i,j) = mgridn(i,j)
       end if
      enddo

c do the middle r cells in a loop
      do i = 3, nrcells - 2
       ip = i + 1
       ipp = i + 2
       im = i - 1
       imm = i - 2
       do j = 1, nphicells
        if (abs(mgridvr(i,j)).gt.small) then
c courant number
         cou = delret*mgridvr(i,j)
c sign of courant number
         scou = cou / (abs(cou))
c upwind
         fphup = 0.5*((1+scou)*mgridn(i,j) + (1-scou)*mgridn(ip,j))
         fmhup = 0.5*((1+scou)*mgridn(im,j) + (1-scou)*mgridn(i,j))
c lax-Wendroff
         fphlw = 0.5*((1+cou)*mgridn(i,j) + (1-cou)*mgridn(ip,j))
         fmhlw = 0.5*((1+cou)*mgridn(im,j) + (1-cou)*mgridn(i,j))
c calculate the limiter function
         if (cou.gt.0.0) then

          mdiv=1.
          if (abs(mgridn(i,j)).gt.mscale) mdiv=mscale

          rjp1 = (mgridn(i,j)/mdiv - mgridn(im,j)/mdiv)
          rjp2 = (mgridn(ip,j)/mdiv - mgridn(i,j)/mdiv)
          if (abs(rjp2).gt.abs(0.5*rjp1)) then 
           rjp = (rjp1 / rjp2)
          else
           rjp = 2.0
          end if
          if (rjp1*rjp2.lt.0.0) then 
           sip = 0.0
          else 
           sip = max(min(2*rjp,1.),min(rjp,2.))
          end if

          rjm1 = (mgridn(im,j)/mdiv - mgridn(imm,j)/mdiv)
          rjm2 = (mgridn(i,j)/mdiv - mgridn(im,j)/mdiv)
          if (abs(rjm2).gt.abs(0.5*rjm1)) then 
           rjm = (rjm1 / rjm2)
          else
           rjm = 2.0
          end if
          if (rjm1*rjm2.lt.0.0) then
           sim = 0.0
          else 
           sim = max(min(2*rjm,1.),min(rjm,2.))
          end if

         else

          mdiv=1.
          if (abs(mgridn(i,j)).gt.mscale) mdiv=mscale 

          rjp1 = (mgridn(ipp,j)/mdiv - mgridn(ip,j)/mdiv)
          rjp2 = (mgridn(ip,j)/mdiv - mgridn(i,j)/mdiv)
c          print *,'rjp #1:',rjp1,rjp2,mgridn(ipp,j),mgridn(ip,j)
          if (abs(rjp2).gt.abs(0.5*rjp1)) then 
           rjp = (rjp1 / rjp2)
          else
           rjp = 2.0
          end if
c          print *,'rjp #2:',rjp1,rjp2
          if (rjp1*rjp2.lt.0.0) then
           sip = 0.0
          else 
           sip = max(min(2*rjp,1.),min(rjp,2.))
          end if

          rjm1 = (mgridn(ip,j)/mdiv - mgridn(i,j)/mdiv)
          rjm2 = (mgridn(i,j)/mdiv - mgridn(im,j)/mdiv)
          if (abs(rjm2).gt.abs(0.5*rjm1)) then 
           rjm = (rjm1 / rjm2)
          else
           rjm = 2.0
          end if
          if (rjm1*rjm2.lt.0.0) then 
           sim = 0.0
          else 
           sim = max(min(2*rjm,1.),min(rjm,2.))
          end if

         end if

c difference
         fph = fphup + (fphlw - fphup)*sip
         fmh = fmhup + (fmhlw - fmhup)*sim
         mgridhalf(i,j) = mgridn(i,j) - (cou*(fph - fmh))
        else
         mgridhalf(i,j) = mgridn(i,j)
        end if
       enddo
      enddo

c do the second to last r cell
      i = nrcells - 1
      ip = i + 1
      im = i - 1
      do j = 1, nphicells
       if (abs(mgridvr(i,j)).gt.small) then
        if (mgridvr(i,j).gt.0.0) then
         mgridhalf(i,j) = mgridn(i,j)  - 
     *      (mgridvr(i,j)*delret) *
     *      (mgridn(i,j) - mgridn(im,j))
        else
         mgridhalf(i,j) = mgridn(i,j)  - 
     *      (mgridvr(i,j)*delret) *
     *      (mgridn(ip,j) - mgridn(i,j))
        end if
       else
        mgridhalf(i,j) = mgridn(i,j)
       end if
      enddo

c do the last r cell
      i = nrcells
      im = i - 1
      do j = 1, nphicells
       if (mgridvr(i,j).gt.small) then
        mgridhalf(i,j) = mgridn(i,j)  - 
     *     (mgridvr(i,j)*delret) *
     *     (mgridn(i,j) - mgridn(im,j))
       else
        mgridhalf(i,j) = mgridn(i,j)
       end if
      enddo

c Advect the azimuthal component next
c Due to rap around of the phi cells, the two columns on the edge 
c are calculated separately outside of the phi loop

      delphit = delt / delphi

c Do the phi cells in a loop
      do j = 1, nphicells

      if (j.eq.1) then
       jp = j + 1
       jpp = j + 2
       jm = nphicells
       jmm  = nphicells - 1
      else if (j.eq.2) then
       jp = j + 1
       jpp = j + 2
       jm = j - 1
       jmm  = nphicells
      else if (j.eq.nphicells-1) then
       jp = j + 1
       jpp = 1
       jm = j - 1
       jmm  = j - 2
      else if (j.eq.nphicells) then
       jp = 1
       jpp = 2
       jm = j - 1
       jmm  = j - 2
      else
       jp = j + 1
       jpp = j + 2
       jm = j - 1
       jmm  = j - 2
      end if

       do i = 1, nrcells

        if (abs(mgridvp(i,j)).gt.small) then
c courant number
         cou = delphit*mgridvp(i,j)
c sign of courant number
         scou = cou/abs(cou)
c upwind
         fphup = 0.5*((1+scou)*mgridhalf(i,j) + 
     *      (1-scou)*mgridhalf(i,jp))
         fmhup = 0.5*((1+scou)*mgridhalf(i,jm) + 
     *      (1-scou)*mgridhalf(i,j))
c lax-Wendroff
         fphlw = 0.5*((1+cou)*mgridhalf(i,j) + 
     *      (1-cou)*mgridhalf(i,jp))
         fmhlw = 0.5*((1+cou)*mgridhalf(i,jm) + 
     *      (1-cou)*mgridhalf(i,j))
c calculate the limiter function
         if (cou.gt.0.0) then

          mdiv=1.
          If (abs(mgridhalf(i,j)).gt.mscale) mdiv=mscale

          rjp1 = (mgridhalf(i,j)/mdiv - mgridhalf(i,jm)/mdiv)
          rjp2 = (mgridhalf(i,jp)/mdiv - mgridhalf(i,j)/mdiv)
          if (abs(rjp2).gt.abs(0.5*rjp1)) then 
           rjp = (rjp1 / rjp2)
          else
           rjp = 2.0
          end if
          if (rjp1*rjp2.lt.0.0) then 
           sip = 0.0
          else 
           sip = max(min(2*rjp,1.),min(rjp,2.))
          end if

          rjm1 = (mgridhalf(i,jm)/mdiv - mgridhalf(i,jmm)/mdiv)
          rjm2 = (mgridhalf(i,j)/mdiv - mgridhalf(i,jm)/mdiv)
          if (abs(rjm2).gt.abs(0.5*rjm1)) then 
           rjm = (rjm1 / rjm2)
          else
           rjm = 2.0
          end if
          if (rjm1*rjm2.lt.0.0) then
           sim = 0.0
          else 
           sim = max(min(2*rjm,1.),min(rjm,2.))
          end if

         else

          mdiv=1.
          If (abs(mgridhalf(i,j)).gt.mscale) mdiv=mscale
          rjp1 = (mgridhalf(i,jpp)/mdiv - mgridhalf(i,jp)/mdiv)
          rjp2 = (mgridhalf(i,jp)/mdiv - mgridhalf(i,j)/mdiv)
          if (abs(rjp2).gt.abs(0.5*rjp1)) then 
           rjp = (rjp1 / rjp2)
          else
           rjp = 2.0
          end if
          if (rjp1*rjp2.lt.0.0) then
           sip = 0.0
          else 
           sip = max(min(2*rjp,1.),min(rjp,2.))
          end if

          rjm1 = (mgridhalf(i,jp)/mdiv - mgridhalf(i,j)/mdiv)
          rjm2 = (mgridhalf(i,j)/mdiv - mgridhalf(i,jm)/mdiv)
          if (abs(rjm2).gt.abs(0.5*rjm1)) then 
           rjm = (rjm1 / rjm2)
          else
           rjm = 2.0
          end if
          if (rjm1*rjm2.lt.0.0) then 
           sim = 0.0
          else 
           sim = max(min(2*rjm,1.),min(rjm,2.))
          end if

         end if

c difference
         fph = fphup + (fphlw - fphup)*sip
         fmh = fmhup + (fmhlw - fmhup)*sim
         mgridn(i,j) = mgridhalf(i,j) - (cou*(fph - fmh))
        else
         mgridn(i,j) = mgridhalf(i,j)
        end if
       enddo

      enddo

      return
      end

ccccccccccccccccccccccccccccc
ccc subroutine epotsimple ccc
ccccccccccccccccccccccccccccc

      subroutine epotsimple(dtheta,dphi,kp,pot)

c input: dtheta, and dphi in degrees
c phi is zero at 24 MLT positive towards dawn
c theta is zero at the pole
      real dtheta, dphi
c Input: Kp index
      real kp

c output: pot potential in volts
      real pot

c   Assuming a uniform dawn-dusk electric field
c      and an assumed kp relationship
c
c     Stagnation point (Re)     Electric Field (V/Re)   Kp
c           10                         919               1
c            9                        1134
c            8                        1436
c            7                        1875
c            6                        2552
c            5                        3675               7

      real pi, rad
      real chi,dr,sint
      real dx

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree

      chi = 7350.0 / (9.0 - kp)

      sint = sin(dtheta * rad)
      dr = 1.0 / (sint * sint)
      dx = dr * sin((180.0+dphi)*rad)

      pot = -chi * dx

      return
      end

cccccccccccccccccccccccccccc
ccc subroutine epotsojka ccc
cccccccccccccccccccccccccccc

      subroutine epotsojka(dtheta,dphi,kp,pot)

c input: dphi,dtheta (spherical coordinates) in degrees
c theta is zero at the pole and 
c phi is zero at 24 MLT positive towards dawn
      real dtheta, dphi
c Input: Kp index
      real kp

c output: pot potential in volts
      real pot

      real pi, rad
      real theta, phi
      real theta_eq, theta_pc, theta_max, chi_0
      real phix, theta_3, chi_pc, chi
      real theta_1, theta_2

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree

      theta=dtheta*rad
      phi=dphi*rad

      theta_eq=(25.0+(2.0*kp))*rad
      theta_pc=(15.0+(0.3*kp))*rad
      theta_max=(theta_pc)+(0.3*(theta_eq-theta_pc)*(abs(sin(phi))))
      chi_0=(10.0+(6.5*kp))*1000.0

      if ((theta.lt.theta_pc).and.(theta.ge.0.0)) then

       phix=asin(sqrt(1.0-((sin(theta)*sin(theta)*
     *    cos(phi)*cos(phi))/(sin(theta_pc)*sin(theta_pc)))))
       theta_3=((sin(theta)*sin(phi))/(sqrt((sin(theta_pc)*
     *    sin(theta_pc))-(sin(theta)*sin(theta)*cos(phi)*cos(phi)))))

       if  ((phix.ge.0.0).and.(phix.le.(pi/3.0))) then
        chi_pc=chi_0*sin((3.0*phix)/(2.0))
        chi=chi_pc*theta_3

       else if ((phix.gt.(pi/3.0)).and.
     *    (phix.le.((2.0*pi)/3.0))) then
        chi_pc=chi_0
        chi=chi_pc*theta_3

       else if ((phix.gt.((2.0*pi)/3.0)).and.
     *    (phix.le.((4.0*pi)/3.0))) then
        chi_pc=-chi_0*cos((3.0*phix)/(2.0))
        chi=chi_pc*theta_3

       else if ((phix.gt.((4.0*pi)/3.0)).and.
     *    (phix.le.((5.0*pi)/3.0))) then
        chi_pc=-chi_0
        chi=chi_pc*theta_3

       else if ((phix.gt.((5.0*pi)/3.0)).and.
     *    (phix.le.((6.0*pi)/3.0))) then
        chi_pc=-chi_0*sin((3.0*phix)/(2.0))
        chi=chi_pc*theta_3

       end if

      else if ((theta.lt.theta_max).and.
     *   (theta.ge.theta_pc)) then

       theta_1=(1.0-(((theta-theta_pc)*(theta-theta_pc))/
     *    ((theta_max-theta_pc)*(theta_eq-theta_pc))))

       if ((phi.ge.0.0).and.(phi.le.(pi/3.0))) then
        chi_pc=chi_0*sin((3.0*phi)/(2.0))
        chi=chi_pc*theta_1

       else if ((phi.gt.(pi/3.0)).and.
     *    (phi.le.((2.0*pi)/3.0))) then
        chi_pc=chi_0
        chi=chi_pc*theta_1

       else if ((phi.gt.((2.0*pi)/3.0)).and.
     *    (phi.le.((4.0*pi)/3.0))) then
        chi_pc=-chi_0*cos((3.0*phi)/(2.0))
        chi=chi_pc*theta_1

       else if ((phi.gt.((4.0*pi)/3.0)).and.
     *    (phi.le.((5.0*pi)/3.0))) then
        chi_pc=-chi_0
        chi=chi_pc*theta_1

       else if ((phi.gt.((5.0*pi)/3.0)).and.
     *    (phi.le.((6.0*pi)/3.0))) then
        chi_pc=-chi_0*sin((3.0*phi)/(2.0))
        chi=chi_pc*theta_1

       end if

      else if ((theta.lt.theta_eq).and.
     *   (theta.ge.theta_max)) then

       theta_2=((((theta-theta_eq)*(theta-theta_eq))/
     *    ((theta_eq-theta_max)*(theta_eq-theta_pc))))

       if  ((phi.ge.0.0).and.(phi.le.(pi/3.0))) then
        chi_pc=chi_0*sin((3.0*phi)/(2.0))
        chi=chi_pc*theta_2

       else if ((phi.gt.(pi/3.0)).and.
     *    (phi.le.((2.0*pi)/3.0))) then
        chi_pc=chi_0
        chi=chi_pc*theta_2

       else if ((phi.gt.((2.0*pi)/3.0)).and.
     *    (phi.le.((4.0*pi)/3.0))) then
        chi_pc=-chi_0*cos((3.0*phi)/(2.0))
        chi=chi_pc*theta_2

       else if ((phi.gt.((4.0*pi)/3.0)).and.
     *    (phi.le.((5.0*pi)/3.0))) then
        chi_pc=-chi_0
        chi=chi_pc*theta_2

       else if ((phi.gt.((5.0*pi)/3.0)).and.
     *    (phi.le.((6.0*pi)/3.0))) then
        chi_pc=-chi_0*sin((3.0*phi)/(2.0))
        chi=chi_pc*theta_2

       end if

      else if ((theta.le.(pi/2.0)).and.
     *   (theta.ge.theta_eq)) then
       chi=0.0
      end if

      pot=chi

      return
      end

cccccccccccccccccccccccccccccc
ccc subroutine getmgridpot ccc
cccccccccccccccccccccccccccccc

      subroutine getmgridpot(nrcells,nphicells,vrcells,vthetacells,
     *   vphicells,mgridpot,par)

c Get the electric potential on the grid

c Input: nrcells, nphicells array index
      integer nrcells, nphicells
c Input: vrcells in Re, vthetacells, vphicells in degrees 
      real vrcells(nrcells), vthetacells(nrcells), vphicells(nphicells)
c Input: par index
      real par(2)

c Output: mgridpot in volts
      real mgridpot(nrcells,nphicells)

      real pi, rad, re
      integer i, j
      real pot

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters

      do i = 1, nrcells
       do j = 1, nphicells
        mgridpot(i,j) = 0.0
        if (par(1).eq.1.0) then
         call epotsimple(vthetacells(i),vphicells(j),par(2),pot)
         mgridpot(i,j) = mgridpot(i,j) + pot
        else if (par(1).eq.2.0) then
         call epotsojka(vthetacells(i),vphicells(j),par(2),pot)
         mgridpot(i,j) = mgridpot(i,j) + pot
        else
         print*,'par(1) ne 1 or 2'
         call CON_stop('ERROR in pbo.f')
        endif
       enddo
      enddo

      return
      end

cccccccccccccccccccccccccccccc
ccc subroutine addcorotpot ccc
cccccccccccccccccccccccccccccc

      subroutine addcorotpot(nrcells,nphicells,vrcells,vthetacells,
     *   vphicells,mgridpot)

c Get the corotation electric potential on the grid

c Input: nrcells, nphicells array index
      integer nrcells, nphicells
c Input: vrcells in Re, vthetacells, vphicells in degrees 
      real vrcells(nrcells), vthetacells(nrcells), vphicells(nphicells)

c Output: mgridpot in volts
      real mgridpot(nrcells,nphicells)

      real pi, re, w, mu, m
      integer i, j
      real coro

      pi = 3.14159                    ! rad
      re = 6.378e6                    ! radius of Earth in meters
      w = (2.0*pi) / (24.0*3600.0)    ! rad/sec
      mu = 4.0*pi*1.0e-7              ! newtons/amps**2
      m = 8.05e22                     ! amps*meter**2

      do i = 1, nrcells
       coro = - (w*mu*m) / (4.0*pi*vrcells(i)*re)
       do j = 1, nphicells
        mgridpot(i,j) = mgridpot(i,j) + coro
       enddo
      enddo

      return
      end

cccccccccccccccccccccccccc
ccc subroutine gradpot ccc
cccccccccccccccccccccccccc

      subroutine gradpot(nrcells,nphicells,vrcells,delr,delphi,
     *   mgridpot,mgrider,mgridep)

c Calculates the two components of the electric field
c from the gradient of the electric potential on the grid

c Input: nrcells, nphicells array index
      integer nrcells, nphicells
c Input: vrcells in Re
      real vrcells(nrcells)
c Input: delr in Re, delphi in degrees
      real delr, delphi
c Input: mgridpot in volts
      real mgridpot(nrcells,nphicells)

c Output: mgrider, mgridep in volts/meter
      real mgrider(nrcells,nphicells), mgridep(nrcells,nphicells)

      real pi, rad, re
      integer i, j, ip, im, jp, jm

      pi = 3.14159          ! rad
      rad = pi / 180.0      ! rad/degree
      re = 6.378e6          ! radius of Earth in meters

      i = 1
      ip = i + 1
      im = i
      j = 1
      jp = j + 1
      jm = nphicells
      mgrider(i,j) = - (1.0/re) * 
     *   ((mgridpot(ip,j)-mgridpot(im,j))/(delr))
      mgridep(i,j) = - (1.0/(vrcells(i)*re)) *
     *   ((mgridpot(i,jp)-mgridpot(i,jm))/(2.0*delphi*rad))
      do j = 2, nphicells-1
       jp = j + 1
       jm = j - 1
       mgrider(i,j) = - (1.0/re) * 
     *    ((mgridpot(ip,j)-mgridpot(im,j))/(delr))
       mgridep(i,j) = - (1.0/(vrcells(i)*re)) *
     *    ((mgridpot(i,jp)-mgridpot(i,jm))/(2.0*delphi*rad))
      enddo
      j = nphicells
      jp = 1
      jm = nphicells - 1
      mgrider(i,j) = - (1.0/re) * 
     *   ((mgridpot(ip,j)-mgridpot(im,j))/(delr))
      mgridep(i,j) = - (1.0/(vrcells(i)*re)) *
     *   ((mgridpot(i,jp)-mgridpot(i,jm))/(2.0*delphi*rad))

      do i = 2, nrcells-1
       ip = i + 1
       im = i - 1
       j = 1
       jp = j + 1
       jm = nphicells
       mgrider(i,j) = - (1.0/re) * 
     *    ((mgridpot(ip,j)-mgridpot(im,j))/(2.0*delr))
       mgridep(i,j) = - (1.0/(vrcells(i)*re)) *
     *    ((mgridpot(i,jp)-mgridpot(i,jm))/(2.0*delphi*rad))
       do j = 2, nphicells-1
        jp = j + 1
        jm = j - 1
        mgrider(i,j) = - (1.0/re) * 
     *     ((mgridpot(ip,j)-mgridpot(im,j))/(2.0*delr))
        mgridep(i,j) = - (1.0/(vrcells(i)*re)) *
     *     ((mgridpot(i,jp)-mgridpot(i,jm))/(2.0*delphi*rad))
       enddo
       j = nphicells
       jp = 1
       jm = nphicells - 1
       mgrider(i,j) = - (1.0/re) * 
     *    ((mgridpot(ip,j)-mgridpot(im,j))/(2.0*delr))
       mgridep(i,j) = - (1.0/(vrcells(i)*re)) *
     *    ((mgridpot(i,jp)-mgridpot(i,jm))/(2.0*delphi*rad))
      enddo

      i = nrcells
      ip = i
      im = i - 1
      j = 1
      jp = j + 1
      jm = nphicells
      mgrider(i,j) = - (1.0/re) * 
     *   ((mgridpot(ip,j)-mgridpot(im,j))/(delr))
      mgridep(i,j) = - (1.0/(vrcells(i)*re)) *
     *   ((mgridpot(i,jp)-mgridpot(i,jm))/(2.0*delphi*rad))
      do j = 2, nphicells-1
       jp = j + 1
       jm = j - 1
       mgrider(i,j) = - (1.0/re) * 
     *    ((mgridpot(ip,j)-mgridpot(im,j))/(delr))
       mgridep(i,j) = - (1.0/(vrcells(i)*re)) *
     *    ((mgridpot(i,jp)-mgridpot(i,jm))/(2.0*delphi*rad))
      enddo
      j = nphicells
      jp = 1
      jm = nphicells - 1
      mgrider(i,j) = - (1.0/re) * 
     *   ((mgridpot(ip,j)-mgridpot(im,j))/(delr))
      mgridep(i,j) = - (1.0/(vrcells(i)*re)) *
     *   ((mgridpot(i,jp)-mgridpot(i,jm))/(2.0*delphi*rad))

      return
      end

ccccccccccccccccccccccc
ccc subroutine hunt ccc
ccccccccccccccccccccccc

      SUBROUTINE HUNT(XX,N,X,JLO)
c
c if x.le.min(xx) then jlo = 0
c if x.gt.max(xx) then jlo = n
c else xx(jlo) < x < xx(jlo+1)
c
      integer n, jlo
      real XX(N), x
      integer inc, jhi, jm
      LOGICAL ASCND
      ASCND=XX(N).GT.XX(1)

      jlo = 1

      IF(JLO.LE.0.OR.JLO.GT.N)THEN
        JLO=0
        JHI=N+1
        GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
        IF(JHI.GT.N)THEN
          JHI=N+1
        ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
          JLO=JHI
          INC=INC+INC
          GO TO 1
        ENDIF
      ELSE
        JHI=JLO
2       JLO=JHI-INC
        IF(JLO.LT.1)THEN
          JLO=0
        ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
          JHI=JLO
          INC=INC+INC
          GO TO 2
        ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)RETURN
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
        JLO=JM
      ELSE
        JHI=JM
      ENDIF
      GO TO 3
      END

ccccccccccccccccccccccccc
ccc subroutine locate ccc
ccccccccccccccccccccccccc

      SUBROUTINE LOCATE(XX,N,X,J)
      DIMENSION XX(N)
      JL=0
      JU=N+1
10    IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
      GO TO 10
      ENDIF
      J=JL
      RETURN
      END

cccccccccccccccccccccccccccccccccc
ccc subroutine interpol2dpolar ccc
cccccccccccccccccccccccccccccccccc

      subroutine interpol2dpolar(theta1,ntheta1,phi1, 
     *   nphi1,data1,theta2,ntheta2,phi2,nphi2,data2)

c interpolate/extrapolate values from data1 into data2

      integer ntheta1, nphi1, ntheta2, nphi2
      real theta1(ntheta1),phi1(nphi1),data1(ntheta1,nphi1)
      real theta2(ntheta2),phi2(nphi2),data2(ntheta2,nphi2)

      integer i, ii, j, jj, jjp
      real v1, v2, v3, v4, stheta, sphi

      ii = 1
      jj = 1

      do i = 1, ntheta2
       if (theta2(i).lt.0.0) then
        print*,'interpol2dpolar: theta2(',i,') is less than zero'
        call CON_stop('ERROR in pbo.f')
       endif
       if (theta2(i).gt.90.0) then
        print*,'interpol2dpolar: theta2(',i,') is greater than 90'
        call CON_stop('ERROR in pbo.f')
       endif
       call hunt(theta1,ntheta1,theta2(i),ii)
       if (ii.eq.0) ii = 1
       if (ii.eq.ntheta1) ii = ntheta1 - 1
       stheta = (theta2(i)-theta1(ii))/(theta1(ii+1)-theta1(ii))
       do j = 1, nphi2
        if (phi2(j).lt.0.0) then
         print*,'interpol2dpolar: phi2(',j,') is less than zero'
         call CON_stop('ERROR in pbo.f')
        endif
        if (phi2(j).gt.360.0) then
         print*,'interpol2dpolar: phi2(',j,') is greater than 360'
         call CON_stop('ERROR in pbo.f')
        endif
        call hunt(phi1,nphi1,phi2(j),jj)
        if (jj.eq.0) then
         jj = nphi1
         jjp = 1
         if (phi1(1).lt.phi1(nphi1)) then
          sphi = (phi2(j) - (phi1(jj)-360.0))/(phi1(jjp) - 
     *       (phi1(jj)-360.0))
         else
          sphi = (phi2(j) - (phi1(jj)+360.0))/(phi1(jjp) - 
     *       (phi1(jj)+360.0))
         endif
        else if (jj.eq.nphi1) then
         jjp = 1
         if (phi1(1).lt.phi1(nphi1)) then
          sphi = (phi2(j) - phi1(jj))/((phi1(jjp)+360.0) - 
     *       phi1(jj))
         else
          sphi = (phi2(j) - phi1(jj))/((phi1(jjp)-360.0) - 
     *       phi1(jj))
         endif
        else
         jjp = jj + 1
         sphi = (phi2(j) - phi1(jj))/(phi1(jjp) - phi1(jj))
        endif
        v1 = data1(ii,jj)
        v2 = data1(ii+1,jj)
        v3 = data1(ii+1,jjp)
        v4 = data1(ii,jjp)
        data2(i,j) = ((1-stheta)*(1-sphi)*v1) + 
     *     (stheta*(1-sphi)*v2) + (stheta*sphi*v3) + 
     *     ((1-stheta)*sphi*v4)

       enddo
      enddo

      return
      end

cccccccccccccccccccccccccc
ccc subroutine savet96 ccc
cccccccccccccccccccccccccc

      subroutine savet96(vthetacells,nthetacells,vphicells,nphicells,
     *   mgridvol,mgridx,mgridy,mgridoc,parmod)

      use ModIoUnit, ONLY : UNITTMP_
      
c Input: nthetacells, nphicells array index
      integer nthetacells, nphicells
c Input: vthetacells, vphicells in degrees 
      real vthetacells(nthetacells), vphicells(nphicells)
c Output: mgridvol in m**3 / weber
      real mgridvol(nthetacells,nphicells)
c mgridx, mgridy, in Re
      real mgridx(nthetacells,nphicells)
      real mgridy(nthetacells,nphicells)
c mgridoc, open(0) or closed(1) table
      real mgridoc(nthetacells,nphicells)
      real parmod(10)

      open(unit = UNITTMP_, file='t96_2.dato', status = 'new',
     *   form = 'formatted')
      write(UNITTMP_,*) parmod
      write(UNITTMP_,*) nthetacells, nphicells
      write(UNITTMP_,*) vthetacells
      write(UNITTMP_,*) vphicells
      write(UNITTMP_,*) mgridx
      write(UNITTMP_,*) mgridy
      write(UNITTMP_,*) mgridvol
      write(UNITTMP_,*) mgridoc
      close(unit = UNITTMP_)

      return
      end

cccccccccccccccccccccccccc
ccc subroutine readt96 ccc
cccccccccccccccccccccccccc

      subroutine readt96(vthetacells,nthetacells,vphicells,nphicells,
     *   mgridvol,mgridx,mgridy,mgridoc,parmod)

      use ModIoUnit, ONLY : UNITTMP_

c Input: nthetacells, nphicells array index
      integer nthetacells, nphicells
c Input: vthetacells, vphicells in degrees 
      real vthetacells(nthetacells), vphicells(nphicells)
c Output: mgridvol in m**3 / weber
      real mgridvol(nthetacells,nphicells)
c mgridx, mgridy, in Re
      real mgridx(nthetacells,nphicells)
      real mgridy(nthetacells,nphicells)
c mgridoc, open(0) or closed(1) table
      real mgridoc(nthetacells,nphicells)
      real parmod(10)

      open(unit = UNITTMP_, file='t96_2.dato', status = 'old',
     *   form = 'formatted')
      read(UNITTMP_,*) parmod
      read(UNITTMP_,*) nthetacells, nphicells
      read(UNITTMP_,*) vthetacells
      read(UNITTMP_,*) vphicells
      read(UNITTMP_,*) mgridx
      read(UNITTMP_,*) mgridy
      read(UNITTMP_,*) mgridvol
      read(UNITTMP_,*) mgridoc
      close(unit = UNITTMP_)

      return
      end

c Scratch space

c23456789012345678901234567890123456789012345678901234567890123456789012

c Input: nrcells, nphicells array index
c      integer nrcells, nphicells
c Input: vrcells in Re, vthetacells, vphicells in degrees 
c      real vrcells(nrcells), vthetacells(nrcells), vphicells(nphicells)
c Input: delr in Re, delphi in degrees
c      real delr, delphi
c Output: mgridb in tesla
c      real mgridb(nrcells,nphicells)
c Output: mgridbi in tesla
c      real mgridbi(nrcells,nphicells)
c Output: mgridpot in volts
c      real mgridpot(nrcells,nphicells)
c Output: mgrider, mgridep in volts/meter
c      real mgrider(nrcells,nphicells), mgridep(nrcells,nphicells)
c Output: mgridvr in meter/sec, and mgridvp in degree /sec
c      real mgridvr(nrcells,nphicells), mgridvp(nrcells,nphicells)
c Output: mgridn in particles / weber
c      real mgridn(nrcells,nphicells)
c Output: mgridden in particles / m**3
c      real mgridden(nrcells,nphicells)
c Output: mgridvol in m**3 / weber
c      real mgridvol(nrcells,nphicells)
c mgridx, mgridy, in Re
c      real mgridx(nthetacells,nphicells)
c      real mgridy(nthetacells,nphicells)
c mgridoc, open(0) or closed(1) table
c      real mgridoc(nthetacells,nphicells)

c      real pi, rad, re
c      integer i, j, ip, im, jp, jm

c      pi = 3.14159          ! rad
c      rad = pi / 180.0      ! rad/degree
c      re = 6.378e6          ! radius of Earth in meters

