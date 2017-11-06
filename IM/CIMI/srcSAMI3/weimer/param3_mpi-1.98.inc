*******************************************
*******************************************

!       PARAM3_MPI-1.98.INC

*******************************************
*******************************************

!      total number of workers
       parameter ( numwork = 16 )

!      number of altitudes (e.g., field lines)

       parameter ( nf   = 120,
!       parameter ( nf   = 90,
     .             nfp1 = nf + 1,
     .             nfm1 = nf - 1,
     .             nfm2 = nf - 2  )

!      number of 'banded latitude'

       parameter ( nf1   = 20,
     .             nf2   = 66,
     .             nf3   = 112,
     .             nf4   = nf  )

!!       parameter ( nf1   = 20,
!!     .             nf2   = 60,
!!     .             nf3   = 106,
!!     .             nf4   = nf  )

!       parameter ( nf1   = 10,
!     .             nf2   = 30,
!     .             nf3   = 62,
!     .             nf4   = nf  )

!      number of grid cells along field line
!      nseg needs to be even
     
       parameter ( nz0 = 200, nextra = 2, nseg = 2,
!       parameter ( nz0 = 200, nextra = 0, nseg = 0,
     .             nze = nseg * nextra, nz = nz0+nze,
     .             nzp1 = nz + 1,
     .             nzm1 = nz - 1  )

!      number of grid cells in longitudinal direction per worker

       parameter ( nl   = 8,
     .             nlp1 = nl + 1,
     .             nlm1 = nl - 1  )

!      Total number of grid cells in longitudinal direction 

       parameter ( nlt   = numwork*(nl - 2),
     .             nltp1 = nlt + 1,
     .             nltm1 = nlt - 1  )

!      grid for madala potential solver 

       parameter ( nnx =  nlt + 1, nny = nf - 1 )
!       parameter ( nnx =  nlt + 1, nny = nf2 - 1 )
       parameter ( nyextra = 5, nnyt = nny + nyextra )

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

       parameter ( nchem = 21 )

!      various constants

!      ftnchek is giving some meaningless errors about precision, 
!      but i am going to lower the precision of some of these
!      variables to keep down the error messages


       parameter ( pie    = 3.1415927   )
       parameter ( po180  = 1.745329e-02 )
       parameter ( rtod   = 57.295780    )
       parameter ( tm18   = 1.e-18       ) 
       parameter ( charge = 4.8032e-10   )

       parameter ( spday  = 86400., sphr = 3600.   )
       parameter ( sol    = 3.e10        )

       parameter ( gzero  = 980.665, re = 6370., bmag = 0.31 )  

       parameter ( bolt   = 1.38044e-16 )
       parameter ( amu    = 1.67252e-24 )
       parameter ( evtok  = 1.1604e+04  ) 


!      (plat,plon) =  coords of north magnetic pole.  
!      the geomagnetic longitude is measured east from the  
!      geomagnetic prime meridian at 291 degrees east geographic. 

!      now in namelist 

!       parameter ( plat = 1.375, plon = 5.079  )

!       parameter ( linesuv = 105, linesnt = 4 )  ! fism

       parameter ( linesuv = 37, linesnt = 4 )    ! euvac

       parameter ( dayve = 80., sidyr = 365.4, solinc = 23.5 )

!      these are for the error function

       parameter ( pas =   .3275911, z1 =  .2548295,
     .             z2  = - .284496 , z3 = 1.421413,
     .             z4  = -1.453152 , z5 = 1.0614054  )


