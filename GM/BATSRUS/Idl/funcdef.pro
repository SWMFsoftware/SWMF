;^CFG COPYRIGHT VAC_UM
;===========================================================================
function funcdef,xx,w,func,physics,eqpar,wnames
;
; Originally developed for the Versatile Advection Code by G. Toth (1996-99).
; Rewritten for the BATSRUS code by G. Toth (2000).
;
; This is the function called by ".r animate" and ".r plotfunc".
;
; "xx" array contains the "ndim" components of the coordinates for the grid.
; "w" array contains the "nw" conservative variables for the grid.
; "func" string describes the function to be returned.
; "physics" string defines the meaning of the variables in "w".
; "eqpar" array contains the equation parameters.
; "wnames" string array contains the names for the conservative variables.
;
; The "func" string is interpreted by the following rules:
;
;   If the first character is '-' it is stripped off and the function
;       is multiplied by -1 before returned.
;
;   A single number between 0 and nw-1 returns the variable indexed by
;       that number, e.g. '2' in 3D returns w(*,*,*,2). 
;
;   Variable names in the "wnames" array mean the appropriate variable.
;
;   Function names listed below are calculated and returned.
;
;   Expressions formed from the standard conservative variable names,
;     (i.e. rho ux uy uz p bx by bz), coordinate names (x, y, z),
;     and equation parameter names (gamma) and any IDL function
;     and operator are evaluated and returned. 
;
; Examples for valid strings: '3', 'rho', '-T', 'rho*ux', ...
;
; One can use "funcdef" interactively too, e.g.
; 
; ekin =funcdef(x,w,'(ux^2+uy^2)*rho','mhd23')
; cfast=funcdef(x,w,'cfast','mhd33',eqpar)
;
; Note that "eqpar" is needed for the pressure but not for the kinetic energy, 
; while "wnames" is not needed in either case.
;===========================================================================

; In 1D xx(n1), in 2D xx(n1,n2,2), in 3D xx(n1,n2,n3,3)
siz=size(xx)
ndim=siz(0)-1
if ndim eq 0 then ndim=1
n1=siz(1)
if ndim gt 1 then n2=siz(2)
if ndim gt 2 then n3=siz(3)
; For 1 variable: w(n1),   w(n1,n2),   w(n1,n2,n3)
; for more      : w(n1,nw),w(n1,n2,nw),w(n1,n2,n3,nw)
siz=size(w)
if siz(0) eq ndim then nw=1 else nw=siz(ndim+1)

if n_elements(wnames) eq 0 then wnames=strarr(nw)

; Check for a negative sign in func
if strmid(func,0,1) eq '-' then begin 
   f=strmid(func,1,strlen(func)-1)
   sign=-1
endif else begin
   f=func
   sign=1
endelse

; Check if f is among the variable names listed in wnames or if it is a number
for iw=0,nw-1 do $
   if f eq strtrim(string(iw),2) or f eq wnames(iw) then $
   case ndim of
      1:result=w(*,iw)
      2:result=w(*,*,iw)
      3:result=w(*,*,*,iw)
   endcase

if n_elements(result) gt 0 then return,sign*result

; Extract variables from w based on name
if not keyword_set(physics) then begin
   print,'Error in funcdef: physics=',physics,' is missing'
   retall
endif

; Extract phys (physics without the numbers) and ndir from physics
l=strlen(physics)
phys=strmid(physics,0,l-2)
ndir=long(strmid(physics,l-1,1))
if ndim ne long(strmid(physics,l-2,1)) then begin
   print,'Error in funcdef: physics=',physics,$
         ' is inconsistent with ndim=',ndim
   retall
endif

; Put elements of eqpar into simple variables based on phys
neqpar = n_elements(eqpar)
if neqpar ne 0 then case 1 of
   phys eq 'mhd' or phys eq 'var' or phys eq 'raw': begin 
       gamma=eqpar(0)
       if neqpar gt 1 then clight=eqpar(1) else clight=59.958
   end
   phys eq 'ray': R_ray=eqpar(0)
   phys eq 'vol': begin
      dtheta=eqpar(0) & dphi=eqpar(1)
   end
   else :
endcase

if phys eq 'mhd' or phys eq 'raw' or phys eq 'ful' or phys eq 'var' then begin
   ; Extract primitive variables
   ux=0 & uy=0 & uz=0 & bx=0 & by=0 & bz=0
   for iw=0,nw-1 do case ndim of
   1: case wnames(iw) of
      'rho': rho=w(*,iw)
      'ux' : ux=w(*,iw)
      'uy' : uy=w(*,iw)
      'uz' : uz=w(*,iw)
      'bx' : bx=w(*,iw)
      'by' : by=w(*,iw)
      'bz' : bz=w(*,iw)
      'b1' : bx=w(*,iw)
      'b2' : by=w(*,iw)
      'b3' : bz=w(*,iw)
      'b1x': b1x=w(*,iw)
      'b1y': b1y=w(*,iw)
      'b1z': b1z=w(*,iw)
      'p'  : p=w(*,iw)
      'pth': p=w(*,iw)
      'mx' : ux=w(*,iw)/rho
      'my' : uy=w(*,iw)/rho
      'mz' : uz=w(*,iw)/rho
      'm1' : ux=w(*,iw)/rho
      'm2' : uy=w(*,iw)/rho
      'm3' : uz=w(*,iw)/rho
      'e'  : e=w(*,iw)
      'e1' : e1=w(*,iw)
      else :
      endcase
   2: case wnames(iw) of
      'RHO': rho=w(*,*,iw)
      'UX' : ux=w(*,*,iw)
      'UY' : uy=w(*,*,iw)
      'UZ' : uz=w(*,*,iw)
      'BX' : bx=w(*,*,iw)
      'BY' : by=w(*,*,iw)
      'BZ' : bz=w(*,*,iw)
      'P'  : p=w(*,*,iw)
      'rho': rho=w(*,*,iw)
      'ux' : ux=w(*,*,iw)
      'uy' : uy=w(*,*,iw)
      'uz' : uz=w(*,*,iw)
      'bx' : bx=w(*,*,iw)
      'by' : by=w(*,*,iw)
      'bz' : bz=w(*,*,iw)
      'b1' : bx=w(*,*,iw)
      'b2' : by=w(*,*,iw)
      'b3' : bz=w(*,*,iw)
      'b1x': b1x=w(*,*,iw)
      'b1y': b1y=w(*,*,iw)
      'b1z': b1z=w(*,*,iw)
      'p'  : p=w(*,*,iw)
      'pth': p=w(*,*,iw)
      'mx' : ux=w(*,*,iw)/rho
      'my' : uy=w(*,*,iw)/rho
      'mz' : uz=w(*,*,iw)/rho
      'm1' : ux=w(*,*,iw)/rho
      'm2' : uy=w(*,*,iw)/rho
      'm3' : uz=w(*,*,iw)/rho
      'e'  : e=w(*,*,iw)
      'e1' : e1=w(*,*,iw)
      else :
      endcase
   3: case wnames(iw) of
      'RHO': rho=w(*,*,*,iw)
      'UX' : ux=w(*,*,*,iw)
      'UY' : uy=w(*,*,*,iw)
      'UZ' : uz=w(*,*,*,iw)
      'BX' : bx=w(*,*,*,iw)
      'BY' : by=w(*,*,*,iw)
      'BZ' : bz=w(*,*,*,iw)
      'P'  : p=w(*,*,*,iw)
      'rho': rho=w(*,*,*,iw)
      'ux' : ux=w(*,*,*,iw)
      'uy' : uy=w(*,*,*,iw)
      'uz' : uz=w(*,*,*,iw)
      'bx' : bx=w(*,*,*,iw)
      'by' : by=w(*,*,*,iw)
      'bz' : bz=w(*,*,*,iw)
      'b1' : bx=w(*,*,*,iw)
      'b2' : by=w(*,*,*,iw)
      'b3' : bz=w(*,*,*,iw)
      'b1x': b1x=w(*,*,*,iw)
      'b1y': b1y=w(*,*,*,iw)
      'b1z': b1z=w(*,*,*,iw)
      'p'  : p=w(*,*,*,iw)
      'pth': p=w(*,*,*,iw)
      'mx' : ux=w(*,*,*,iw)/rho
      'my' : uy=w(*,*,*,iw)/rho
      'mz' : uz=w(*,*,*,iw)/rho
      'm1' : ux=w(*,*,*,iw)/rho
      'm2' : uy=w(*,*,*,iw)/rho
      'm3' : uz=w(*,*,*,iw)/rho
      'e'  : e=w(*,*,*,iw)
      'e1' : e1=w(*,*,*,iw)
      else :
      endcase
   endcase
   ; Some extra variables
   uu=ux^2+uy^2+uz^2
   bb=bx^2+by^2+bz^2
   ; Change energy into pressure
   if n_elements(p) eq 0 and n_elements(e) gt 0 then $
      p=(gamma-1)*(e-0.5*(rho*uu+bb))
   ; Calculate gamma*p+bb if that is needed for the function
   if strmid(f,1,4) eq 'fast' or strmid(f,1,4) eq 'slow' then cc=gamma*p+bb
endif else begin
   case physics of
   'ray13':begin
      bx =w(*,0) & by =w(*,1) & bz =w(*,2)
      theta1=w(*,3) & phi1=w(*,4) & theta2=w(*,5) & phi2=w(*,6) 
      status=w(*,7) & blk=w(*,8)
      f1x=w(*,9)  & f1y=w(*,10) & f1z=w(*,11)
      f2x=w(*,12) & f2y=w(*,13) & f2z=w(*,14)
   end
   'ray23':begin
      bx =w(*,*,0) & by =w(*,*,1) & bz =w(*,*,2)
      theta1=w(*,*,3) & phi1=w(*,*,4) & theta2=w(*,*,5) & phi2=w(*,*,6) 
      status=w(*,*,7) & blk=w(*,*,8)
      f1x=w(*,*,9)  & f1y=w(*,*,10) & f1z=w(*,*,11)
      f2x=w(*,*,12) & f2y=w(*,*,13) & f2z=w(*,*,14)
   end
   'ray33':begin
      bx =w(*,*,*,0) & by =w(*,*,*,1) & bz =w(*,*,*,2)
      theta1=w(*,*,*,3) & phi1=w(*,*,*,4)
      theta2=w(*,*,*,5) & phi2=w(*,*,*,6) 
      status=w(*,*,*,7) & blk=w(*,*,*,8)
      f1x=w(*,*,*,9)  & f1y=w(*,*,*,10) & f1z=w(*,*,*,11)
      f2x=w(*,*,*,12) & f2y=w(*,*,*,13) & f2z=w(*,*,*,14)
   end
   'vol22':begin
      vol1=w(*,*,0) & vol2=w(*,*,1)
   end
   else: begin
     print,'Error in funcdef: physics=',physics,' is not known'
     retall
   end
   endcase
endelse

; Extract coordinate variables from x:
case ndim of
   1:begin
      x=xx
   end
   2:begin
      x=xx(*,*,0) & y=xx(*,*,1)
   end
   3:begin
      x=xx(*,*,*,0) & y=xx(*,*,*,1) & z=xx(*,*,*,2)
   end
endcase

;==== Put your function definition(s) below using the variables and eq. params
;     extracted from x, w and eqpar, and select cases by the function name "f"
;     and the "phys, ndim, ndir" variables extracted from "physics".

case f of
  'Btheta'   : result=atan(by,sqrt(bx^2+bz^2))
  ; Electric field
  'Ex'       : result=(by*uz-uy*bz)
  'Ey'       : result=(bz*ux-uz*bx)
  'Ez'       : result=(bx*uy-ux*by)
  ; momenta
  'mx'       : result=rho*ux
  'my'       : result=rho*uy
  'mz'       : result=rho*uz
  ; Boris momenta
  'mBx'      : result=rho*ux + (bb*ux - (ux*bx+uy*by+uz*bz)*bx)/clight^2
  'mBy'      : result=rho*uy + (bb*uy - (ux*bx+uy*by+uz*bz)*by)/clight^2
  'mBz'      : result=rho*uz + (bb*uz - (ux*bx+uy*by+uz*bz)*bz)/clight^2
  ; Total energy
  'e'        : result=p/(gamma-1)+0.5*(rho*uu + bb)
  ; Alfven speed and Mach number
  'calfvenx' : result=bx/sqrt(rho)
  'calfveny' : result=by/sqrt(rho)
  'calfvenz' : result=bz/sqrt(rho)
  'calfven'  : result=sqrt(bb/rho)
  'Malfvenx' : result=ux/bx*sqrt(rho)
  'Malfveny' : result=uy/by*sqrt(rho)
  'Malfvenz' : result=uz/bz*sqrt(rho)
  'Malfven'  : result=sqrt(uu/bb*rho)
  ; pressure, plazma beta, temperature, entropy
  'pbeta'    : result=2*p/bb
  'T'        : result=p/rho
  's'        : result=p/rho^gamma
  ; sound speed, Mach number
  'csound'   :result=sqrt(gamma*p/rho)
  'mach'  :result=sqrt(uu/(gamma*p/rho))
  'machx' :result=ux/sqrt(gamma*p/rho)
  'machy' :result=uy/sqrt(gamma*p/rho)
  'machz' :result=uz/sqrt(gamma*p/rho)
  ; fast and slow speeds and Mach numbers
  'cfast' :result=sqrt(cc/rho)
  'cfastx':result=sqrt((cc+sqrt(cc^2-4*gamma*p*bx^2))/2/rho)
  'cfasty':result=sqrt((cc+sqrt(cc^2-4*gamma*p*by^2))/2/rho)
  'cfastz':result=sqrt((cc+sqrt(cc^2-4*gamma*p*bz^2))/2/rho)
  'cslowx':result=sqrt((cc-sqrt(cc^2-4*gamma*p*bx^2))/2/rho)
  'cslowy':result=sqrt((cc-sqrt(cc^2-4*gamma*p*by^2))/2/rho)
  'cslowz':result=sqrt((cc-sqrt(cc^2-4*gamma*p*bz^2))/2/rho)
  'Mfast' :result=sqrt(mm/cc/rho)
  'Mfastx':result=ux/sqrt((cc+sqrt(cc^2-4*gamma*p*bx^2))/2/rho)
  'Mfasty':result=uy/sqrt((cc+sqrt(cc^2-4*gamma*p*by^2))/2/rho)
  'Mfastz':result=uz/sqrt((cc+sqrt(cc^2-4*gamma*p*bz^2))/2/rho)
  'Mslowx':result=ux/sqrt((cc-sqrt(cc^2-4*gamma*p*bx^2))/2/rho)
  'Mslowy':result=uy/sqrt((cc-sqrt(cc^2-4*gamma*p*by^2))/2/rho)
  'Mslowz':result=uz/sqrt((cc-sqrt(cc^2-4*gamma*p*bz^2))/2/rho)
  'jz': case ndim of
     1: result=deriv(x,by)
     2: result=curl(bx,by,x,y)
     3: begin
          print,'Error in funcdef: jz is not implemented for 3D yet'
          retall
        end
  endcase

  ; Magnetic vector potential A in slab symmetry
  ; The density of contourlines is proportional to B
  'Ax': begin
     result=dblarr(n1,n2)
     ; Integrate along the first row
     for i1=1,n1-1 do result(i1,0)=result(i1-1,0) $
               +(bz(i1,0)+bz(i1-1,0))*(x(i1,0)-x(i1-1,0))*0.5 $
               -(by(i1,0)+by(i1-1,0))*(y(i1,0)-y(i1-1,0))*0.5
     ; Integrate all columns vertically
     for i2=1,n2-1 do result(*,i2)=result(*,i2-1) $
               +(bz(*,i2)+bz(*,i2-1))*(x(*,i2)-x(*,i2-1))*0.5 $
               -(by(*,i2)+by(*,i2-1))*(y(*,i2)-y(*,i2-1))*0.5
  end
  'Ay': begin
     result=dblarr(n1,n2)
     ; Integrate along the first row
     for i1=1,n1-1 do result(i1,0)=result(i1-1,0) $
               -(bz(i1,0)+bz(i1-1,0))*(x(i1,0)-x(i1-1,0))*0.5 $
               +(bx(i1,0)+bx(i1-1,0))*(y(i1,0)-y(i1-1,0))*0.5
     ; Integrate all columns vertically
     for i2=1,n2-1 do result(*,i2)=result(*,i2-1) $
               -(bz(*,i2)+bz(*,i2-1))*(x(*,i2)-x(*,i2-1))*0.5 $
               +(bx(*,i2)+bx(*,i2-1))*(y(*,i2)-y(*,i2-1))*0.5
  end
  'Az': begin
     result=dblarr(n1,n2)
     ; Integrate along the first row
     for i1=1,n1-1 do result(i1,0)=result(i1-1,0) $
               +(by(i1,0)+by(i1-1,0))*(x(i1,0)-x(i1-1,0))*0.5 $
               -(bx(i1,0)+bx(i1-1,0))*(y(i1,0)-y(i1-1,0))*0.5
     ; Integrate all columns vertically
     for i2=1,n2-1 do result(*,i2)=result(*,i2-1) $
               +(by(*,i2)+by(*,i2-1))*(x(*,i2)-x(*,i2-1))*0.5 $
               -(bx(*,i2)+bx(*,i2-1))*(y(*,i2)-y(*,i2-1))*0.5
  end

  ; Some functions used in specific problems

  'e1':result=p/(gamma-1)+0.5*(rho*uu+$
       w(*,*,*,10)^2+w(*,*,*,11)^2+w(*,*,*,12)^2)

  'e1B':result=p/(gamma-1)+0.5*(rho*uu + $
       w(*,*,*,10)^2+w(*,*,*,11)^2+w(*,*,*,12)^2 + $
       ((by*uz-uy*bz)^2+(bz*ux-uz*bx)^2+(bx*uy-ux*by)^2)/clight^2)

  'bb1':result=w(*,*,*,10)^2+w(*,*,*,11)^2+w(*,*,*,12)^2

  'Jy_cut': begin
	result=diff2(2,w(*,*,9),z)-diff2(1,w(*,*,11),x)
	loc=where(x^2+z^2 lt 16.,count)
	if count gt 0 then result(loc)=0.0
            end

  'logrho':result=alog10(rho)

  ;If "f" has not matched any function yet, try evaluating it as an expression
  else:begin
     if sign eq -1 then begin
        f='-'+f
        sign=1
     endif
     if not execute('result='+f) then begin
        print,'Error in funcdef: cannot evaluate function=',func
        retall
     endif
  end
endcase

if n_elements(result) gt 0 then begin
   return,sign*result
endif else begin
   print,'Error in funcdef: function=',func,' was not calculated ?!'
   retall
endelse

end
