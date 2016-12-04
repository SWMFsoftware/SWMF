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
; 1 - Output file name is set to SPECTRUM.input.erg
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
  LogTmin=3.0
  LogTmax=9.0
  dLogT=0.1
  nLogT=(LogTmax-LogTmin)/dLogT
  LogT=LogTmin+dLogT*findgen(nLogT+1.0)
  nLogT=n_elements(LogT)

; List of elements
  Element=['h' ,'he','li','be','b' ,'c' ,'n' ,'o' ,'f' ,'ne','na','mg','al', $
           'si','p' ,'s' ,'cl','ar','k' ,'ca','sc','ti','v' ,'cr','mn','fe', $
           'co','ni','cu','zn']

  nElement=n_elements(Element)
  
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

; Input the abundance file
  if keyword_set(abund_unity) then abundance=fltarr(nElement)+1.0d0
  if not keyword_set(abund_unity) then read_abund,'',abundance,ref

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
  for i=0,nIon-1 do begin
     position=strpos(MasterList(i),'d')
     if position ge 2 then Dielectronic=[Dielectronic,i]
  endfor

  if n_elements(Dielectronic) gt 1 then begin
     Dielectronic=Dielectronic(1:*)
     remove,Dielectronic,MasterList
     nIon = n_elements(MasterList)
  endif

; Set output file name
  NameFileOut='SPECTRUM.input.erg'

; Write header into the output file
  openw,10,NameFileOut
  
  printf,10,'log(G(LogN,LogT)) table created with CHIANTI : ',version, $
         format='(a42)'
  printf,10,'#GRID',format='(a)'
  printf,10,LogNmin,LogNmax,dLogN
  printf,10,LogTmin,LogTmax,dLogT
  printf,10,'number of ions = ',nIon,format='(a17,i3)'
  printf,10,'ion ','lvl1 ','lvl2 ','wavelength ','logn ','logT ',$
         'logG(n,T) ',format='(a4,2a5,a11,2a5,a9)'
  printf,10,'in units [A] [cm^-3] [K] [erg cm^3 s^-1]'
  printf,10,'#START',format='(a)'

; START main loop f calculation and writing
  for i=0,nIon-1 do begin
     
; Select ion
     position=strpos(MasterList(i),'_')
     LocalElement=strmid(MasterList(i),0,position)
     zElem=where(Element eq LocalElement)+1
     zElem=zElem(0)
     Ion=MasterList(i)
     zIon=float(strmid(Ion,position+1,strlen(Ion)))
     print,'doing ion ',MasterList(i),' which is the ',i,'th ion out of ',nIon

; Select ioneq 
     ion_eq=ioneq(*,zElem-1,zIon-1)

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
     LogG=LogG*ElementAbundance*0.83
     for k=0,nLogT-1 do begin
        for j=0,nLogN-1 do begin
           if not keyword_set(ioneq_unity) then begin
              LogG(k,j,*)=LogG(k,j,*)*ion_eq(k)/(10^LogN(j))
              position=abs(logT(k)-t_ioneq)
           endif
           if keyword_set(ioneq_unity) then LogG(k,j,*)=LogG(k,j,*)/(10^LogN(j))
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
                    printf,10,MasterList(i),LevelFrom(k),LevelTo(k),$
                           Wavelength(k),LogN(j),$
                           LogT(l),LogG(l,j,k),$
                           format='(a6,2i5,f13.3,2f13.4,f13.3)'
                 endif
                 if Flag(k) lt 0 then begin
                    printf,10,MasterList(i),LevelFrom(k),LevelTo(k),$
                           Wavelength(k)*Flag(k),$
                           LogN(j),LogT(l),LogG(l,j,k),$
                           format='(a6,2i5,f13.3,2f13.4,f13.3)'
                 endif
              endfor
           endfor
        endif
     endfor
  endfor
  close,10
end
