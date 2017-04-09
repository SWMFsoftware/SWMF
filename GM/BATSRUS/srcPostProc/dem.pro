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

; Set Debug flag
  read, IsDebug, prompt = 'Debug mode? (1 means yes ) ='

; Debug setup
  if IsDebug eq 1 then begin
     IsTest = 0
     IsVerbose = 1

     filename = 'debug-box.out'
     read_data

     iRho = 0
     IsTe = 0
     iPe = 10

     LogTmin = 5
     LogTmax = 7
     DLogT = .05

  endif else begin

; Test setup 
     read, IsTest, prompt = ' Testing? (1 means yes ) = '
     if IsTest eq 1 then begin 
        IsDebug = 0
        IsVerbose = 1
        NameFileOut = 'test-dem.eps'

        iRho = 0
        IsTe = 1
        iTe = 1

        LogTmin = 1
        LogTmax = 3
        DLogT = 1

        w =  MAKE_ARRAY(1000, 1, 1, 3, /DOUBLE, VALUE = 0.)
        x = w
        x[0,*,*] = 1.

        w[0:10,*,*,iTe] = 10
        w[11:999,*,*,iTe] = 99
        
        w[0:10,*,*,iRho] = Mp
        w[11:999,*,*,iRho] = 100*Mp

        DX = 1./RSun

     endif else begin

; Command-line setup
        IsDebug = 0
        IsTest = 0
        read, IsVerbose, prompt = ' Verbose? (1 means yes ) = '
        filename = ''
        NameFileOut=''
        read, NameFileOut, prompt = 'output file = '
        read, filename, prompt = 'data file = '
        read_data

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

        read, LogTmin, prompt = 'LogTmin = '
        read, LogTmax, prompt = 'LogTmax = '
        read, DLogT, prompt = 'DLogT = '
     endelse
  endelse

; Data grid size
  sizex = size(x)
  nX = sizex[1]
  nY = sizex[2]
  nZ = sizex[3]
  DX = (max(x[*,0,0,0]) - min(x[*,0,0,0]))/max([1,nX-1])
  DY = (max(x[0,*,0,1]) - min(x[0,*,0,1]))/max([1,nY-1])
  DZ = (max(x[0,0,*,2]) - min(x[0,0,*,2]))/max([1,nZ-1])

; Set up temperature grid
  NT = round((LogTmax - LogTmin) / DLogT)
  Tarray = MAKE_ARRAY(NT+1, /DOUBLE, VALUE = 0.)
  for i = 0,NT do begin
     Tarray[i] = 10.^LogTmin * 10.^(DLogT * i)
  endfor
  DEMarray = MAKE_ARRAY(NT+1, /DOUBLE, VALUE = 1e-10)

; In case of verbose print grid information
  if IsVerbose eq 1 then begin
     print, 'LogTmin, LogTmax, DLogT, NT, Tarray = ', $
            LogTmin, LogTmax, DLogT, NT, Tarray
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
                 DEMarray[l] = DEMarray[l] + (w[i,j,k,iRho]/Mp)^2*Dh/10.^DLogT
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

  if IsDebug ne 1 then begin
     close_device, /pdf
  endif else begin
; Plot on screen
     set_plot,'X'
     window, 1
     plot,Tarray,DEMarray,/xlog,/ylog, psym=2, symsize=1, linestyle=2, thick=2
  endelse
end
