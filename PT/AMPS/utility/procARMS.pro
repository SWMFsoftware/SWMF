pro procARMS, infile, outfile
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; procedure to read ARMS (Adaptively Refined MHD Solver) data 
; from IDL .sav file and to write it to output files 'outfile.t=0.dat',...
; USAGE example:
;   procARMS, './structure_data.sav', 'dataARMS' 
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  RESTORE, infile
;data is restored to struc (all variables are in cgs units):
;variable name     description                                 units
;---------------------------------------------------------------------------
;STEP      = time-step counter in ARMS simulation   ;#
;TIME      = simulation time in seconds             ;s
;NX        = number of points in the x-direction    ;#
;NZ        = number of points in the z-direction    ;#
;XPOS      = x-position of the grid point (cm)      ;cm
;ZPOS      = z-position of the grid point (cm)      ;cm
;N_E       = electron mass density                  ;1/cm^3
;PRE       = thermal pressure                       ;dyn/cm^2
;TEM       = temperature                            ;K
;VEL_RAD   = velocity in the radial direction       ;cm/s
;VEL_THETA = velocity in the theta direction        ;cm/s 
;VEL_PHI   = velocity in the phi direction          ;cm/s
;B_RAD     = magnetic field in the radial direction ;Gauss    
;B_THETA   = magnetic field in the theta direcion   ;Gauss 
;B_PHI     = magnetic field in the phi direcion     ;Gauss 
;J_RAD     = current density in the radial direction;gr^(1/2)/(cm^(1/2) s^2)=esu/(cm^2 s)
;J_THETA   = current density in the theta direction ;gr^(1/2)/(cm^(1/2) s^2)=esu/(cm^2 s)
;J_PHI     = current density in the phi direction   ;gr^(1/2)/(cm^(1/2) s^2)=esu/(cm^2 s)
;ALF       = local Alfven speed                     ;cm/s
;BETA      = plasma beta                            ;
;FLUX      = poloidal flux function                 ;arbitrary units
;---------------------------------------------------------------------------
;e.g., to get the value of the magnetic field in the radial direction 
;at time = t0, and the grid point [x=x0,z=z0], use: struc(t0).B_RAD(x0,y0)
;To get the x-position of the grid point x=x0, use: struc(t0).XPOS(x0,y0)
;To visualize the radial component of the field: tvscl,
;struc(t0).B_RAD(*,*)
;---------------------------------------------------------------------------

; write data to files
nTime = n_elements(struc.TIME)
header = 'Data from ARMS simulation; data is 2.5D, sampled on uniform X-Z grid; processed from IDL structure'
for iTime=0, nTime-1 do begin
   curfile = outfile + '.t=' + strcompress(string(iTime),/remove_all) + '.dat'
   nX = struc(iTime).NX
   nZ = struc(iTime).NZ
   openw,  lun, curfile, /get_lun
   printf, lun, header
; current time ----------------------------------------------------------------
   printf, lun, '#TIME'
   printf, lun, struc(iTime).TIME - struc(0).TIME,format='(E13.5)'
; mext time -------------------------------------------------------------------
   printf, lun, '#TIME_NEXT'
   if(iTime eq nTime-1)then begin
      printf, lun, '-1'
   endif else begin
      printf, lun, struc(iTime+1).TIME - struc(0).TIME,format='(E13.5)'
   endelse
   printf, lun, '#NX'
   printf, lun, nX, format='(I4)'
   printf, lun, '#NZ'
   printf, lun, nZ, format='(I4)'
; grid ------------------------------------------------------------------------
   printf, lun, '#XPOS'
   for iX = 0, nX-1 do begin
      printf, lun, struc(iTime).XPOS(iX),format='(E13.5,A1,$)'
   endfor
   printf, lun, ''
   printf, lun, '#ZPOS'
   for iZ = 0, nZ-1 do begin
      printf, lun, struc(iTime).ZPOS(iZ),format='(E13.5,A1,$)'
   endfor
   printf, lun, ''
; electron mass density -------------------------------------------------------
   printf, lun, '#N_E'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).N_E(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; thermal pressure ------------------------------------------------------------
   printf, lun, '#PRE'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).PRE(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; temperature -----------------------------------------------------------------
   printf, lun, '#TEM'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).TEM(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; velocity in radial direction ------------------------------------------------
   printf, lun, '#VEL_RAD'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).VEL_RAD(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; velocity in theta direction -------------------------------------------------
   printf, lun, '#VEL_THETA'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).VEL_THETA(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; velocity in phi direction ---------------------------------------------------
   printf, lun, '#VEL_PHI'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).VEL_PHI(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; magnetic field in radial direction ------------------------------------------
   printf, lun, '#B_RAD'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).B_RAD(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; magnetic field in theta direction -------------------------------------------
   printf, lun, '#B_THETA'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).B_THETA(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; magnetic field in phi direction ---------------------------------------------
   printf, lun, '#B_PHI'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).B_PHI(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; current density in radial direction -----------------------------------------
   printf, lun, '#J_RAD'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).J_RAD(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; current density in theta direction ------------------------------------------
   printf, lun, '#J_THETA'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).J_THETA(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; current density in phi direction --------------------------------------------
   printf, lun, '#J_PHI'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).J_PHI(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; Alfven speed ----------------------------------------------------------------
   printf, lun, '#ALF'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).ALF(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; plasma beta -----------------------------------------------------------------
   printf, lun, '#BETA'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).BETA(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
; poloidal flux function ------------------------------------------------------
   printf, lun, '#FLUX'
   for iZ = 0, nZ-1 do begin
      for iX = 0, nX-1 do begin
         printf, lun, struc(iTime).FLUX(iX,iZ),format='(E13.5,A1,$)'
      endfor
      printf, lun, ''
   endfor
   printf, lun, ''
   printf, lun, '#END'
   free_lun, lun
endfor

end
