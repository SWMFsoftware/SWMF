; procedure for synthetic DEM calculation

pro dem

  common getpict_param
  common plot_data
  common file_head

; Constants
  Mp = 1.6726e-27

; Testing or not
  read, test, prompt = ' testing ? (1 means yes ) = '

; String holders
  filename = ''
  filenameout=''

; If test is on do not use prompt
  if test eq 1 then begin 
     filename = 'dem-test-box.out'
     filenameout = 'dem-test-out.eps'
     read_data

     tmin = 1
     tmax = 101
     dt = 1

     iRho = 0
     iTe = 11

     w[0:10,*,*,iTe] = 10
     w[11:1000,*,*,iTe] = 99
     
     w[0:10,*,*,iRho] = 1* 1e-3*Mp
     w[11:1000,*,*,iRho] = 100* 1e-3*Mp

  endif else begin 
; Get data file
     read, filename, prompt = 'filename = '
     read_data
; Indexes of variables
     read,iRho, prompt = 'index of electron density variable = '
     read,iTe, prompt = 'index of electron temperature variable = '
     iRho = fix(iRho)
     iTe = fix(iTe)

; Set up temperature grid
     read, tmin, prompt = 'Tmin = '
     read, tmax, prompt = 'Tmax = '
     read, dt, prompt = 'dT = '
     read, filenameout, prompt = 'name of output pdf file = '
  endelse

; Set up temperature and dem arrays
  nt = round((tmax-tmin)/dt)
  tarray = tmin+dt*FINDGEN(nt+1)
  demarray = tarray*0.

; Grid size
  sizex = size(x)
  nx = sizex[1]
  ny = sizex[2]
  nz = sizex[3]
  dx = (max(x[*,0,0,0]) - min(x[*,0,0,0]))/max([1,nx-1])
  dy = (max(x[*,0,0,1]) - min(x[*,0,0,1]))/max([1,ny-1])
  dz = (max(x[*,0,0,2]) - min(x[*,0,0,2]))/max([1,nz-1])

; In case of testing print temperatre array and overwrite dx
  if test eq 1 then begin
     dx = 1./7e8
     print, 'tarray = ', tarray
  endif

; Main loop
  dh = dx * 7e8 ; width of observed plasma
  for k = 0, nz-1 do begin
     for j = 0, ny-1 do begin
        for i = 0, nx-1 do begin
           if w[i,j,k,0] eq 0 then continue ;avoid body cells
           for l = 0, nt-2 do begin ; check temperature
              if (w[i,j,k,iTe] lt tarray[l+1]) and (w[i,j,k,iTe] ge tarray[l])then begin 
                 demarray[l] = demarray[l] + (w[i,j,k,iRho]*1e3/Mp)^2*dh/dt
                 break
              endif
           endfor
        endfor
     endfor
  endfor
  
; Print dem array if testing
  if test eq 1 then begin
     print, 'demarray = ', demarray
  endif

; Set plotting environment
  plot_spacex = 10
  !p.charsize=2
  !p.color=0
  !p.background=255
  !p.thick = 4
  bottomline=0

; Save plot
  set_device,filenameout,/eps,/land
  !x.margin=[4,2]
  !p.title = 'Synthetic Differential Emission Measure'
  !x.title = 'Temperature [K] '
  !y.title = 'DEM [K!U-1 !Ncm!U-5!N]'
  plot,tarray,demarray,/ylog, psym=2, symsize=2, linestyle=2, thick=3
  close_device, /pdf

end
