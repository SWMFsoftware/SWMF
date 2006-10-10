
C ============================================================================

      SUBROUTINE prntCollision
      use ModCommonVariables

      write (iUnitCollision,"(a79)") 'Collision Frequencies'
      write (iUnitCollision,"(i8,1pe13.5,3i3)") nint(time/dt),time,1,1,12
      write (iUnitCollision,"(3i4)") nDim
      write (iUnitCollision,"(100(1pe13.5))") Gamma
      write (iUnitCollision,"(a79)") 
     &    'alt HpH3p ELHp ELH3p HpH ELH2 H3pHp HpEL H3pEL HpH2 H3pH H3pH2 ELH gamma'
      
      DO K=1,NDIM
         QS1=CFHpH3p(K)
         QS2=CFELHp(K)
         QS3=CFELH3p(K)
         QS4=CFHpH(K)
         QS5=CFELH2(K)
         QS6=CFH3pHp(K)
         QS7=CFHpEL(K)
         QS8=CFH3pEL(K)
         QS9=CFHpH2(K)
         QS10=CFH3pH(K)
         QS11=CFH3pH2(K)
         QS12=CFELH(K)
CALEX helium not considered at saturn so we need a selector
         
         WRITE (iUnitCollision,"(100(1pe18.10))") 
     &        ALTD(K),QS1,QS2,QS3,QS4,QS5,QS6,QS7,QS8,QS9,QS10,QS11,QS12,gamma
      enddo

      RETURN
      END
