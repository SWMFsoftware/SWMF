; JSzente 2017
; Procedure for synthetic DEM calculation
pro dem

  common getpict_param
  common plot_data
  common file_head

; Constants
  Mp = 1.6726e-27 * 1e3         ; proton mass in g
  KB = 1.38064852e-16           ; Boltzmann Constant in cm2 g s-2 K-1
  RSun = 7e10                   ; solar radius in cm

; Set IsDebug flag
  read, IsDebug, prompt = 'Debug mode? (1 means yes ) ='

; IsDebug setup
  if IsDebug eq 1 then begin
     IsTest = 0
     IsVerbose = 0

     NameFile = 'box8.out'

     iRho = 0
     IsTe = 0
     iTe = 12
     iPe = 10

     Tmin = 1e4
     Tmax = 1e9
     DT = 1e4
     NameFileOut = 'IsDebug-box8'

     !x.margin=[10,3]
  endif else begin
; Command-line setup
     read, IsTest, prompt = ' Testing? (1 means yes ) = '
     read, IsVerbose, prompt = ' Verbose? (1 means yes ) = '
; String holders
     filename = ''
     NameFileOut=''
; Non-testing setup
     if IsTest ne 1 then begin 
        read, filename, prompt = 'data file = '
        read, NameFileOut, prompt = 'name of output file (.pdf)= '
; Get indexes of variables
        read,iRho, prompt = 'index of electron density variable = '
        read,IsTe, prompt = $
             'use electron pressure (0) or electron temperature (1) ='
        if IsTe eq 0 then read,iPe, prompt = $
                               'index of electron pressure variable = '
        if IsTe eq 1 then read,iTe, prompt = $
                               'index of electron temperature variable  = '
        iRho = fix(iRho)
        if IsTe eq 0 then iPe = fix(iPe)
        if IsTe eq 1 then iTe = fix(iTe)
; Get temperature grid information
        read, Tmin, prompt = 'Tmin = '
        read, Tmax, prompt = 'Tmax = '
        read, DT, prompt = 'DT = '
     endif
  endelse

; Generate test-data
  if IsTest eq 1 then begin 
     filename = 'dem-test-box.out'
     NameFileOut = 'dem-test-out.eps'
     read_data

     Tmin = 1
     Tmax = 101
     DT = 1

     IsTe = 1

     iRho = 0
     iTe = 11

     w[0:10,*,*,iTe] = 10
     w[11:1000,*,*,iTe] = 99
     
     w[0:10,*,*,iRho] = Mp
     w[11:1000,*,*,iRho] = 100*Mp

     DX = 1./RSun

  endif

; Data grid size
  sizex = size(x)
  nX = sizex[1]
  nY = sizex[2]
  nZ = sizex[3]
  DX = (max(x[*,0,0,0]) - min(x[*,0,0,0]))/max([1,nX-1])
  DY = (max(x[0,*,0,1]) - min(x[0,*,0,1]))/max([1,nY-1])
  DZ = (max(x[0,0,*,2]) - min(x[0,0,*,2]))/max([1,nZ-1])

; Set up temperature grid
  NT = round((Tmax-Tmin)/DT)
  Tarray = Tmin+DT*FINDGEN(NT+1)
  DEMarray = Tarray*0. + 1e-10

; In case of verbose print grid information
  if IsVerbose eq 1 then begin
     print, 'Tarray = ', Tarray
     print, 'nX, nY, nZ = ',nX,nY,nZ
     print, 'DX, DY, DZ = ',DX,DY,DZ
  endif
  
; Main loop
; Width of observed plasma
  Dh = DX * RSun             
  for k = 0, nZ-1 do begin
     for j = 0, nY-1 do begin
        for i = 0, nX-1 do begin
; Avoid body cells
           if w[i,j,k,0] eq 0 then continue 
; Locate on temperature grid
           for l = 0, NT-2 do begin         
              if IsTe eq 0 then Te = w[i,j,k,iPe]/w[i,j,k,iRho] * Mp/KB
              if IsTe eq 1 then Te = w[i,j,k,iTe]
              if (Te lt Tarray[l+1]) and (Te ge Tarray[l])then begin 
                 DEMarray[l] = DEMarray[l] + (w[i,j,k,iRho]/Mp)^2*Dh/DT
                 break
              endif
           endfor
        endfor
     endfor
  endfor
  
; Print if IsVerbose is on                                
  if IsVerbose eq 1 then begin
     print, 'DEMarray = ', DEMarray
  endif

; Set plotting environment
  plot_spacex = 10
  !p.charsize=2
  !p.color=0
  !p.background=255
  !p.thick = 4
  bottomline=0
  !p.title = 'Synthetic Differential Emission Measure'
  !x.title = 'Temperature [K] '
  !y.title = 'DEM [K!U-1 !Ncm!U-5!N]'

; Save plot when not debugging
  if IsDebug ne 1 then begin
     set_device,NameFileOut,/eps,/land
     !x.margin=[4,2]
  endif

  plot,Tarray,DEMarray,/xlog,/ylog, psym=2, symsize=1, linestyle=2, thick=2,$
       yrange=[1e18,1e28], xrange=[1e4, 1e7]

  if IsDebug ne 1 then  close_device, /pdf

; Plot on screen
  set_plot,'X'
  window, 1
  plot,Tarray,DEMarray,/xlog,/ylog, psym=2, symsize=1, linestyle=2, thick=2

end
