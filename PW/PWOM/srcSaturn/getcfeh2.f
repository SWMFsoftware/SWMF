
C This function numerically integrates the momentum crossesction to 
C get the average crossection for a given temperature. The function returns
C the electron collission frequency with H2

      real FUNCTION getcfeh2(TELECT,XH2,XMSE,UELECT)
      use ModConst, ONLY: cBoltzmann

      REAL sigmaD(29)
      REAL vsigma(29)
      REAL esigma(29)

      REAL c2
      REAL QD

C Data for calculating crossection      
      DATA (sigmaD(i),i=1,29)/0,7.2E-20,8.1E-20,8.5E-20,9E-20,9.2E-20,
     &9.6E-20,9.8E-20,10E-20,10.2E-20,10.5E-20,12E-20,13E-20,13.8E-20,
     &14.5E-20,15.2E-20,16E-20,16.5E-20,17E-20,17.4E-20,18E-20,15.5E-20,
     &13.5E-20,10.8E-20,6.5E-20,5.5E-20,4.8E-20,4E-20,3.4E-20/

C Velocity of incident electrons (m/s)      
      DATA (vsigma(i),i=1,29)/0,59309,83876,1.02727E5,1.18619E5,
     &1.32620E5,
     &1.45278E5,1.56918E5,1.67753E5,1.77929E5,1.87553E5,2.65241E5, 
     &3.24852E5,3.75107E5,4.19382E5,4.59411E5,4.96220E5,5.30482E5,
     &5.62661E5,5.93097E5,8.38766E5,1.027274E6,1.186194E6,1.326205E6,
     &1.452785E6,1.567187E6,1.677531E6,1.779291E6,1.875537E6/ 

C Energy of incident electrons (eV)
      DATA (esigma(i),i=1,29)/0,.01,.02,.03,.04,
     &.05,
     &.06,.07,.08,.09,.1,.2, 
     &.3,.4,.5,.6,.7,.8,
     &.9,1.,2.,3.,4.,5.,6.,7.,8.,9.,10./ 
      
      QD=0.
      c2=2.*cBoltzmann*TELECT/XMSE

C calculate the kinetic energy in eV of the electrons.      
      kinetic=.5*XMSE*(UELECT/(1.E2))**2.*1.60217733e-19
      
      if (kinetic .lt. esigma(2)) then
         QD=sigmaD(2)
      endif

      do i=3,29
         if (kinetic .eq. esigma(i-1)) then 
           QD=sigmaD(i-1) 
         endif
         if ( kinetic .gt. esigma(i-1)) then
            if (kinetic .lt. esigma(i))then
               QD=(sigmaD(i-1)+sigmaD(i))/2.
            endif
         endif
         
         if (kinetic .eq. esigma(i)) then
            QD=sigmaD(i)
         endif
      enddo
      
      if (kinetic .gt. esigma(29)) then
         QD=sigmaD(29)
      endif

      getcfeh2=(1.E6)* XH2*abs(UELECT)*QD

      
      return
      end
