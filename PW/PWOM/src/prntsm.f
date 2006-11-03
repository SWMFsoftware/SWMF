C ============================================================================

      SUBROUTINE PRNTSM
      use ModCommonVariables
      real MO,MH,MHe,Me

      write (iUnitGraphics,"(a79)") 'Polarwind output_var11'
      write (iUnitGraphics,"(i8,1pe13.5,3i3)") nint(time/dt),time,1,1,20
      write (iUnitGraphics,"(3i4)") nDim+2
      write (iUnitGraphics,"(100(1pe13.5))") Gamma
      write (iUnitGraphics,"(a79)")
     &     'r Lat Lon uO uHe uH ue lgnO lgnHe lgnH lgne TO THe TH Te MO MH MHe Me Ef Pe g'
         

C write out lower ghost cell values
         QS1=UsurfO/1.E5
         QS2=UsurHe/1.E5
         QS3=UsurfH/1.E5
         QS4=Usurfe/1.E5
         QS5=alog10_check(DsurfO/XMSO)
CALEX helium not considered at saturn so we need a selector
         if (namePlanet .eq. 'Earth ') then
            QS6=alog10_check(DsurHe/XMSHE)
            MHe=UsurHe/sqrt(gamma*PsurHe/DsurHe)
         endif

         if (namePlanet .eq. 'Saturn') then
            QS6=0.
            MHe=0.
         endif

         QS7=alog10_check(DsurfH/XMSH)
         QS8=alog10_check(Dsurfe/XMSE)
         MO=UsurfO/sqrt(gamma*PsurfO/DsurfO)
         MH=UsurfH/sqrt(gamma*PsurfH/DsurfH)
         Me=Usurfe/sqrt(gamma*Psurfe/Dsurfe)
        
                  
         WRITE (iUnitGraphics,"(100(1pe18.10))") 
     &        AltMin,GmLat, GmLong,QS1,QS2,QS3,QS4,QS5,QS6,QS7,QS8,
     &        TsurfO,TsurHe,TsurfH,Tsurfe,
     &        MO,MH,MHe,Me,Efield(1),Psurfe
         


C Write out for the middle of the grid

      DO K=1,NDIM
         QS1=UOXYG(K)/1.E5
         QS2=UHEL(K)/1.E5
         QS3=UHYD(K)/1.E5
         QS4=UELECT(K)/1.E5
         QS5=alog10_check(DOXYG(K)/XMSO)
CALEX helium not considered at saturn so we need a selector
         if (namePlanet .eq. 'Earth ') then
            QS6=alog10_check(DHEL(K)/XMSHE)
            MHe=UHEL(K)/sqrt(gamma*PHEL(K)/DHEL(K))
         endif

         if (namePlanet .eq. 'Saturn') then
            QS6=0.
            MHe=0.
         endif

         QS7=alog10_check(DHYD(K)/XMSH)
         QS8=alog10_check(DELECT(K)/XMSE)
         MO =UOXYG(K)/sqrt(gamma*POXYG(K)/DOXYG(K))
         MH =UHYD(K)/sqrt(gamma*PHYD(K)/DHYD(K))
         Me =UELECT(K)/sqrt(gamma*PELECT(K)/DELECT(K))
         

         WRITE (iUnitGraphics,"(100(1pe18.10))") 
     &        ALTD(K),GmLat, GmLong,QS1,QS2,QS3,QS4,QS5,QS6,QS7,QS8,
     &        TOXYG(K),THEL(K),THYD(K),TELECT(K),
     &        MO,MH,MHe,Me,Efield(k),Pelect(k)
      enddo


C Write out the upper ghost cell

         QS1=UbgndO/1.E5
         QS2=UbgnHe/1.E5
         QS3=UbgndH/1.E5
         QS4=Ubgnde/1.E5
         QS5=alog10_check(DbgndO/XMSO)
CALEX helium not considered at saturn so we need a selector
         if (namePlanet .eq. 'Earth ') then
            QS6=alog10_check(DbgnHe/XMSHE)
            MHe=UbgnHe/sqrt(gamma*PbgnHe/DbgnHe)
         endif

         if (namePlanet .eq. 'Saturn') then
            QS6=0.
            MHe=0.
         endif

         QS7=alog10_check(DbgndH/XMSH)
         QS8=alog10_check(Dbgnde/XMSE)
         MO=UbgndO/sqrt(gamma*PbgndO/DbgndO)
         MH=UbgndH/sqrt(gamma*PbgndH/DbgndH)
         Me=Ubgnde/sqrt(gamma*Pbgnde/Dbgnde)
                  
         
         WRITE (iUnitGraphics,"(100(1pe18.10))") 
     &        AltMax,GmLat, GmLong, QS1,QS2,QS3,QS4,QS5,QS6,QS7,QS8,
     &        TbgndO,TbgnHe,TbgndH,Tbgnde,
     &        MO,MH,MHe,Me,Efield(nDim),Pbgnde
             



      RETURN
      END

C     ========================================================================
      real function alog10_check(x)
      implicit none
      real, intent(in) :: x

      if(x < 0.0) then
         write(*,*)'negative argument for alog10_check:',x
         call CON_stop('PWOM ERROR: negative argument for alog10')
      endif
      alog10_check=alog10(x)

      end function alog10_check


