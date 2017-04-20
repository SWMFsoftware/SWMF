pro SPECTRUM_input,abund_unity=abund_unity,notstandard=notstandard,photoexc=photoexc,interval=interval

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
; 3 - Delta logT = as in ioneq file
;
; 4 - logG(logT,logNe) are calculated at temperatures defined in
;     ionequation file and densities 1e6-1e14 cm^3.
;
; VERSION
;
;  This is a modified version of crea_files.pro.  
;
; 06-Apr-2017 - Generating full tables into one file.
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

; Wavelengthrange
  wmin = 10
  wmax = 1000
  if keyword_set(interval) then begin
     read,'minimum wavelength [A]= ',wmin
     read,'minimum wavelength [A]= ',wmax
  endif

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
  ioneqfile = '/Applications/ssw/packages/chianti/dbase/ioneq/chianti.ioneq'
  read_ioneq,ioneqfile,t_ioneq,ioneq,ref

; Input the abundance file
  abundancefile= '/Applications/ssw/packages/chianti/dbase/abundance/sun_coronal_1992_feldman.abund'

  if keyword_set(abund_unity) then abundance=fltarr(nElement)+1.0d0
  if not keyword_set(abund_unity) then read_abund,abundancefile,abundance,ref

; Provide photoexcitation parameters (if keyword 'photoexc' is set)
  if keyword_set(photoexc) then begin
     read,'rPhot = ',rPhot
     read,'RadTemp = ',RadTemp
  endif

; Input the masterlist file
  masterfile='/Applications/ssw/packages/chianti/dbase/masterlist/masterlist.ions'
  read_masterlist,masterfile,MasterList
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

; START main loop f calculation and writing
  for i=0L,nIon-1 do begin
     
; Select ion
     position=strpos(MasterList(i),'_')
     LocalElement=strmid(MasterList(i),0,position)
     zElem=where(Element eq LocalElement)+1
     Aion = Aelement(where(Element eq LocalElement)) 
     zElem=zElem(0)
     Ion=MasterList(i)
     zIon=float(strmid(Ion,position+1,strlen(Ion)))
     print,'doing ',MasterList(i),' which is the ',i+1,'th ion out of ',nIon

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
     G=data.em
     Wavelength=data.lambda

     Missing1=where(G lt 1.0d-99)
     Missing2=where(strcompress(string(G),/remove_all) eq 'Infinity')
     Missing3=where(strcompress(string(G),/remove_all) eq '-Infinity')
     Missing4=where(strcompress(string(G),/remove_all) eq 'NaN')
     Missing5=where(strcompress(string(G),/remove_all) eq '-NaN')

     G(Missing1)=1.0d-100
     G(Missing2)=1.0d-100
     G(Missing3)=1.0d-100
     G(Missing4)=1.0d-100
     G(Missing5)=1.0d-100
     
     G=G*ElementAbundance 
     for l=0,nLogT-1 do begin
        for j=0,nLogN-1 do begin
           G(l,j,*)=G(l,j,*)*ion_eq(l)/(10^LogN(j))
           position=abs(LogT(l)-t_ioneq)
        endfor
     endfor
     
     LogG=alog10(G)
     Missing6=where(strcompress(string(LogG),/remove_all) eq '-Infinity')
     
     LogG(Missing1)=-100
     LogG(Missing2)=-100
     LogG(Missing3)=-100
     LogG(Missing4)=-100
     LogG(Missing5)=-100
     LogG(Missing6)=-100

; Write into output file
     nTransition=n_elements(data)   
     help,nTransition
     for k=0L,nTransition-1 do begin
        if Wavelength(k) gt wmin and $
           Wavelength(k) lt wmax and $
           max(LogG(*,*,k)) gt -26 then begin 
           for j=0,nLogN-1 do begin
              for l=0,nLogT-1 do begin
                 if Flag(k) eq 0 then begin
                    printf,10,strtrim(MasterList(i),2)+STRING(9B)+ $
                           strtrim(string(Aion,format='(f13.3)'),2)+STRING(9B)+ $
                           strtrim(string(LevelFrom(k)),2)+STRING(9B)+ $
                           strtrim(string(LevelTo(k)),2)+STRING(9B)+ $
                           strtrim(string(Wavelength(k),format='(f13.3)'),2)+STRING(9B)+ $
                           strtrim(string(LogN(j),format='(f13.1)'),2)+STRING(9B)+ $
                           strtrim(string(LogT(l),format='(f13.2)'),2)+STRING(9B)+ $
                           strtrim(string(LogG(l,j,k),format='(f13.3)'),2)
                 endif
                 if Flag(k) lt 0 then begin
                    printf,10,strtrim(MasterList(i),2)+STRING(9B)+ $
                           strtrim(string(Aion,format='(f13.3)'),2)+STRING(9B)+ $
                           strtrim(string(LevelFrom(k)),2)+STRING(9B)+ $
                           strtrim(string(LevelTo(k)),2)+STRING(9B)+ $
                           strtrim(string(Wavelength(k)*Flag(k),format='(f13.3)'),2)+STRING(9B)+ $
                           strtrim(string(LogN(j),format='(f13.1)'),2)+STRING(9B)+ $
                           strtrim(string(LogT(l),format='(f13.2)'),2)+STRING(9B)+ $
                           strtrim(string(LogG(l,j,k),format='(f13.3)'),2)
                 endif
              endfor
           endfor
        endif
     endfor
  endfor
  close,10
end
