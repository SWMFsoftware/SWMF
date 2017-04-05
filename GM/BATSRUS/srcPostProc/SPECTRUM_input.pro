pro SPECTRUM_input,abund_unity=abund_unity,notstandard=notstandard,photoexc=photoexc

; This program calculates the emissivities logG(logT,logNe) for all the ions in
; the CHIANTI database and prints out into one file to be used for BATS-R-US.
;
; This program needs in input the following things:
;
;               ion fraction file
;               abundance file
;               masterlist file
;               if /notstandard keyword is set, the CHIANTI directory
;
;           Keywords:
;               /notstandard : User chooses CHIANTI version to use
;               /photoexc    : Include photoexcitation in calculation
;               /abund_unity : All element abundances are set to 1.
;
; NOTES
;
; 1 - Output file name is set to SPECTRUM_chianti_tbl.dat
;     (containing ion.erg files data)
;
; 2 - It does not deal with dielectronic satellite lines
;
; 3 - Delta logT = 0.1
;
; 4 - logG(logT,logNe) are calculated at temperatures 1e3-1e9 K
;     and densities 1e6-1e14 cm^3.
;
; VERSION
;
;  This is a modified version of crea_files.pro.  
;
; 30-Mar-2016 - Output into one file only for BATSRUS
  version='8.0'
; Close lun=10 because we are going to use that for writing the output file
  close,10

; Density grid 
  LogN=[6,7,8,9,10,11,12,13,14]*1.0
  nLogN=n_elements(LogN)
  LogNmin=min(LogN)
  LogNmax=max(LogN)
  dLogN=(max(LogN)-min(LogN))/(nLogN-1.0)

; Temperature grid
  LogTmin=4.0
  LogTmax=9.0
  dLogT=0.05
  nLogT=(LogTmax-LogTmin)/dLogT
  LogT=LogTmin+dLogT*findgen(nLogT+1.0)
  nLogT=n_elements(LogT)

; List of elements
  Element=['h' ,'he','li','be','b' ,'c' ,'n' ,'o' ,'f' ,'ne','na','mg','al', $
           'si','p' ,'s' ,'cl','ar','k' ,'ca','sc','ti','v' ,'cr','mn','fe', $
           'co','ni','cu','zn']

  nElement=n_elements(Element)

; Atomic weight of elements
  Aelement=[ 1.008,  4.003,  6.94 ,  9.012, 10.81 , 12.011, 14.007, 15.999, $
             18.998, 20.18 , 22.99 , 24.305, 26.982, $
             28.085, 30.974, 32.06 , 35.45 , 39.948, 39.098, 40.078, 44.956, $
             47.867, 50.942, 51.996, 54.938 ,55.845, $
             58.933, 58.693, 63.546, 65.38]
  
; Set CHIANTI version (if keyword 'notstandard' is set)
  if keyword_set(notstandard) then begin 
     version=''
     versionpath=''
     read,'CHIANTI version (e.g. 6, 5.2.1, set 0 to current version) : '$
          ,version
     read,'CHIANTI version path : ',versionpath
     version=strcompress(versionpath,/remove_all) + $
             strcompress(version,/remove_all)
     use_chianti,version
  endif

; Input the ion fraction file
  read_ioneq,'',t_ioneq,ioneq,ref
;  help,t_ioneq,ioneq

; Input the abundance file
  if keyword_set(abund_unity) then abundance=fltarr(nElement)+1.0d0
  if not keyword_set(abund_unity) then read_abund,'',abundance,ref
;  print,'!!!!!!!! read_abund in SPECTRUM = ', ref
; Provide photoexcitation parameters (if keyword 'photoexc' is set)
  if keyword_set(photoexc) then begin
     read,'rPhot = ',rPhot
     read,'RadTemp = ',RadTemp
  endif

; Input the masterlist file
  read_masterlist,'',MasterList
  nIon=n_elements(MasterList)
  Dielectronic=100.

; Remove dielectronic list elements
  for i=0L ,nIon-1 do begin
     position=strpos(MasterList(i),'d')
     if position ge 2 then Dielectronic=[Dielectronic,i]
  endfor

  if n_elements(Dielectronic) gt 1 then begin
     Dielectronic=Dielectronic(1:*)
     remove,Dielectronic,MasterList
     nIon = n_elements(MasterList)
  endif

; Set output file name
  NameFileOut='SPECTRUM_chianti_tbl.dat'

; Write header into the output file
  openw,10,NameFileOut
  
  printf,10,'log(G(LogN,LogT)) table created with CHIANTI : ',version, $
         format='(a42)'
  printf,10,'#GRID',format='(a)'
  printf,10,LogNmin,LogNmax,dLogN
  printf,10,LogTmin,LogTmax,dLogT
  printf,10,'number of ions = ',nIon,format='(a17,i3)'
  printf,10,'ion ','atomic mass ','lvl1 ','lvl2 ','wavelength ','logn ','logT ',$
         'logG(n,T) ',format='(a4,a12,2a5,a11,2a5,a9)'
  printf,10,'in units [A] [cm^-3] [K] [erg cm^3 s^-1]'
  printf,10,'#START',format='(a)'

;  iNe6 = where(MasterList eq 'ne_6')
;  iNe6 = iNe6(0)
;  print,'!!!!! i, MasterList = ',iNe6, MasterList(iNe6)

  
;  print,'!!!! Element=', Element


; START main loop f calculation and writing
  for i=0L,nIon-1 do begin
;;;for i=iNe6,iNe6 do begin
; Select ion
     position=strpos(MasterList(i),'_')
     LocalElement=strmid(MasterList(i),0,position)
     zElem=where(Element eq LocalElement)+1
     Aion = Aelement(where(Element eq LocalElement)) 
     zElem=zElem(0)

     Ion=MasterList(i)
     zIon=float(strmid(Ion,position+1,strlen(Ion)))

;     print,'!!!! i, LocalElement, zElem=', i, LocalElement, zElem
     print,'doing ',MasterList(i),' which is the ',i+1,'th ion out of ',nIon

; Select ioneq

     ion_eq=ioneq(*,zElem-1,zIon-1)
;     help,ion_eq,zElem,zIon
;     print,'!!!!!!!!!!!! t_ioneq, ioneq = '

;      for it = 0,100 do print, t_ioneq(it),ion_eq(it)

; Limit LogT range to fraction > limit_ioneq
;     limit_ioneq = 0. ;;;;1e-6

;     if zElem-1 eq 0 and zIon-1 eq 0 then begin
;        location = where(ion_eq gt limit_ioneq/100.)
;     endif else begin
;        location = where(ion_eq gt limit_ioneq)
;     endelse
;     LogTmin=t_ioneq(min(location))
;     LogTmax=t_ioneq(max(location))
;     dLogT=0.1
;     nLogT=(LogTmax-LogTmin)/dLogT
;     LogT=LogTmin+dLogT*findgen(nLogT+1)
;     nLogT=n_elements(LogT)
     
;     dt1=t_ioneq(min(location)+1)-t_ioneq(min(location))
;     if dt1 lt 0.09 then begin
;        q1=location(2*findgen(nLogT))
;        location=q1
;     endif

;print,location
;stop


; Calculate the emissivities
     data=emiss_calc(zElem,zIon,dens=LogN,temp=LogT,/quiet)
     
     if keyword_set(photoexc) then data=emiss_calc(zElem,zIon,dens=LogN,$
                                                   temp=LogT,rphot=rPhot,$
                                                   radtemp=RadTemp,/quiet)

; Select the element abundance
     ElementAbundance=abundance(zElem-1)

; Calculate log(G(T,Ne))
     LevelFrom=data.level1
     LevelTo=data.level2
     Flag=data.flag
     LogG=data.em
     Wavelength=data.lambda
     Missing=where(LogG lt 1.0d-99)
     LogG(Missing)=1.0d-45
     Missing=where(strcompress(string(LogG),/remove_all) eq 'Infinity')
     LogG(Missing)=1.0d-45
     Missing=where(strcompress(string(LogG),/remove_all) eq '-Infinity')
     LogG(Missing)=1.0d-45
     Missing=where(strcompress(string(LogG),/remove_all) eq 'NaN')
     LogG(Missing)=1.0d-45
     Missing=where(strcompress(string(LogG),/remove_all) eq '-NaN')
     LogG(Missing)=1.0d-45


     if 1 eq 0 then begin
        print,'i = ', i
        print,' abundance(zElem-1) = ', abundance(zElem-1)
        print, 'ElementAbundance = ', ElementAbundance
        print, 'ElementAbundance*0.83 = ',ElementAbundance*0.83
        j = 2
        l = 4
        k = 2765
   ;;; k = where(abs(wavelength-300.562) lt 0.001)
        print, 'j, LogN = ', j, LogN(j)
        print, 'l, LogT = ', l, LogT(l)
        print, 'k, WaveL= ', k, Wavelength(k)
        g = LogG(l,j,k)
        print, 'LogG=', g
        g = g*ElementAbundance*0.83
        print, 'LogG*abund*0.83=', g
        print,'10^LogN(j)=',10^LogN(j)
        print,'location(k)=', location(l)
        ioneq_ = ion_eq(location(l))/(10^LogN(j))
        print,'ion_eq=', ioneq_
        g = g*ioneq_
        print, 'LogG*abund*0.83*ioneq=', g
        
     endif

;print, 'LogG = ',LogG

     LogG=LogG*ElementAbundance ;;;; *0.83
     for l=0,nLogT-1 do begin
        for j=0,nLogN-1 do begin
           LogG(l,j,*)=LogG(l,j,*)*ion_eq(l)/(10^LogN(j))
           position=abs(logT(l)-t_ioneq)
        endfor
     endfor
     LogG=alog10(LogG)
     Missing=where(strcompress(string(LogG),/remove_all) eq '-Infinity')
     LogG(Missing)=-100.
; Write into output file
     nTransition=n_elements(data)   
     help,nTransition
     for k=0L,nTransition-1 do begin
        if max(LogG(*,*,k)) gt -26.0 then begin
           for j=0,nLogN-1 do begin
              for l=0,nLogT-1 do begin
                 if Flag(k) eq 0 then begin
                    printf,10,MasterList(i),Aion,LevelFrom(k),LevelTo(k),$
                           Wavelength(k),LogN(j),$
                           LogT(l),LogG(l,j,k),$
                           format='(a6,f13.3,2i5,f13.3,2f13.4,f13.3)'
                    
                 endif
                 if Flag(k) lt 0 then begin
                    printf,10,MasterList(i),Aion,LevelFrom(k),LevelTo(k),$
                           Wavelength(k)*Flag(k),$
                           LogN(j),LogT(l),LogG(l,j,k),$
                           format='(a6,f13.3,2i5,f13.3,2f13.4,f13.3)'
                 endif
              endfor
           endfor
        endif
     endfor
  endfor
  close,10
end
