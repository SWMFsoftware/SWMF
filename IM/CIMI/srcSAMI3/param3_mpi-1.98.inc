!*******************************************
!*******************************************

!       PARAM3_MPI-1.98.INC

!*******************************************
!*******************************************

!      total number of workers
       integer,parameter :: numwork = 23

!      number of altitudes (e.g., field lines)

!       parameter ( nf   = 120,
      	integer, parameter :: nf   = 124
      	integer, parameter :: nfp1 = nf + 1
      	integer, parameter :: nfm1 = nf - 1
       	integer, parameter :: nfm2 = nf - 2  

!      number of 'banded latitude'

      	integer, parameter :: nf1   = 20
      	integer, parameter :: nf2   = 66
      	integer, parameter :: nf3   = 112
      	integer, parameter :: nf4   = nf  


!      number of grid cells along field line
!      nseg needs to be even
     
      integer, parameter ::nz0 = 200, nextra = 2, nseg = 2
      integer, parameter ::nze = nseg * nextra, nz = nz0+nze
      integer, parameter ::nzp1 = nz + 1
      integer, parameter ::nzm1 = nz - 1 

!      number of grid cells in longitudinal direction per worker

      integer, parameter ::  nl   = 6
      integer, parameter ::nlp1 = nl + 1
      integer, parameter ::nlm1 = nl - 1  

!      Total number of grid cells in longitudinal direction 

      integer, parameter ::nlt   = numwork*(nl - 2)
      integer, parameter ::nltp1 = nlt + 1
      integer, parameter ::nltm1 = nlt - 1  

!      grid for madala potential solver 

      integer, parameter :: nnx =  nlt + 1, nny = nf - 1 
      integer, parameter :: nyextra = 5, nnyt = nny + nyextra 

!      ion densities

       integer nion,pthp,pthep,ptnp,ptop,ptn2p,ptnop,pto2p

       parameter ( nion  = 7 )   ! number of ions
       parameter ( pthp  = 1 )   ! h+
       parameter ( pthep = 5 )   ! he+
       parameter ( ptnp  = 7 )   ! n+
       parameter ( ptop  = 2 )   ! o+
       parameter ( ptn2p = 6 )   ! n2+
       parameter ( ptnop = 3 )   ! no+
       parameter ( pto2p = 4 )   ! o2+

!      neutrals 

       integer nneut,pth,pthe,ptn,pto,ptn2,ptno,pto2

       parameter ( nneut = 7 )   ! number of neutrals
       parameter ( pth   = 1 )   ! h
       parameter ( pthe  = 5 )   ! he
       parameter ( ptn   = 7 )   ! n
       parameter ( pto   = 2 )   ! o
       parameter ( ptn2  = 6 )   ! n2
       parameter ( ptno  = 3 )   ! no
       parameter ( pto2  = 4 )   ! o2

!      number of chemical reactions

      integer, parameter ::nchem = 21 

!      various constants

!      ftnchek is giving some meaningless errors about precision, 
!      but i am going to lower the precision of some of these
!      variables to keep down the error messages


       real, parameter ::pie    = 3.1415927   
       real, parameter ::po180  = 1.745329e-02 
       real, parameter :: rtod   = 57.295780    
       real, parameter :: tm18   = 1.e-18        
       real, parameter :: charge = 4.8032e-10   

       real, parameter :: spday  = 86400., sphr = 3600.  
       real, parameter :: sol    = 3.e10    

       real, parameter :: gzero  = 980.665, re = 6370., bmag = 0.31 

       real, parameter :: bolt   = 1.38044e-16 
       real, parameter :: amu    = 1.67252e-24 
       real, parameter :: evtok  = 1.1604e+04   


!      (plat,plon) =  coords of north magnetic pole.  
!      the geomagnetic longitude is measured east from the  
!      geomagnetic prime meridian at 291 degrees east geographic. 

!      now in namelist 

!       parameter ( plat = 1.375, plon = 5.079  )

!       parameter ( linesuv = 105, linesnt = 4 )  ! fism

       integer, parameter ::linesuv = 37, linesnt = 4     ! euvac

       real, parameter :: dayve = 80., sidyr = 365.4, solinc = 23.5 

!      these are for the error function

       real, parameter :: pas =   .3275911, z1 =  .2548295
       real, parameter :: z2  = - .284496 , z3 = 1.421413
       real, parameter :: z4  = -1.453152 , z5 = 1.0614054 


