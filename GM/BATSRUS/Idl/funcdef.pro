;^CFG COPYRIGHT VAC_UM
;===========================================================================
function funcdef,xx,w,func,time,eqpar,variables,rcut
;
; Originally developed for the Versatile Advection Code by G. Toth (1996-99).
; Rewritten for the BATSRUS code by G. Toth (2000).
;
; This is the function called by ".r animate" and ".r plotfunc".
;
; "xx"        array contains the "ndim" components of the coordinates.
; "w"         array contains the "nw" variables.
; "func"      string describes the function to be returned.
; "time"      simulatuion time
; "eqpar"     array contains the equation parameters.
; "variables" string array contains the names of the coordinates in xx
;             the variables in w and the equation parameters in eqpar
; "rcut"      radius for central cut out circle
;
; The "func" string is interpreted by the following rules:
;
;   If the first character is '-' it is stripped off and the function
;       is multiplied by -1 before returned.
;
;   A single number between 0 and nw-1 returns the variable indexed by
;       that number, e.g. '2' in 3D returns w(*,*,*,2). 
;
;   Variable names in the "variables" array mean the appropriate variable.
;
;   Function names listed below are calculated and returned.
;
;   Expressions formed from the standard conservative variable names,
;     (i.e. rho ux uy uz p bx by bz), coordinate names (x, y, z),
;     and equation parameter names (gamma) and any IDL function
;     and operator are evaluated and returned. 
;
;   Expressions can also contain coordinate, variable or equation
;   parameter names enclosed in curly brackets. For example
;   {r} or {bx1} or {rbody}. These will be replaced with strings like
;   xx(*,*,0) or w(*,*,10) or eqpar(2) before evaluation.


; Examples for valid strings: 
;   '3', 'rho', 'Mfast', '-T', 'rho*ux', '{bx1}^2+{by1}^2'...
;
; One can use "funcdef" interactively too, e.g.
; 
; ekin =funcdef(x,w,'(ux^2+uy^2)*rho','mhd23',eqpar,variables)
; cfast=funcdef(x,w,'cfast','mhd33',eqpar,variables)
;
;===========================================================================

; In 1D xx(n1), in 2D xx(n1,n2,2), in 3D xx(n1,n2,n3,3)
siz=size(xx)
ndim=siz(0)-1
if ndim eq 0 then ndim=1
n1=siz(1)
if ndim gt 1 then n2=siz(2)
if ndim gt 2 then n3=siz(3)

; For 1 variable: w(n1), w(n1*n2), w(n1,n2),  w(n1,n2,n3)
; for more      : w(n1,nw),w(n1,n2,nw),w(n1,n2,n3,nw)
siz=size(w)
if siz(0) le ndim then nw=1 else nw=siz(ndim+1)

; Define some constants for "planetary" units
amu        = 1.6726e-24         ; [g]
kboltzmann = 1.3807e-23         ; [SI]
tunit      = 1e-15/kboltzmann   ; tunit*p/rho       [nPa] [amu/cc] -> [K]
mu0        = 4*!pi*1e-7*1e9     ; p/(bb/2/mu0)      [nPa] [nT]
mu0A       = 4*!pi*1e-7*amu*1e27; sqrt(bb/mu0A/rho) [nT] [amu/cc] -> [km/s]

; Number of equation parameters
nEqpar = n_elements(eqpar)

; Extract equation parameters
gamma  =  5./3.
clight =  1.0
rbody  = -1.0
if nEqpar gt 0 then begin
    for i=nDim+nW,n_elements(variables)-1 do begin
        iEqpar = i-nDim-nW
        case variables(i) of
            'g'     : gamma  = eqpar(iEqpar)
            'gamma' : gamma  = eqpar(iEqpar)
            'c'     : clight = eqpar(iEqpar)
            'r'     : rbody  = eqpar(iEqpar)
            'rbody' : rbody  = eqpar(iEqpar)
            else:
        endcase
    endfor
endif

; Variable names
if n_elements(variables) eq 0 then variables=strarr(ndim+nw+neqpar)
wnames = strlowcase(variables(ndim:ndim+nw-1))

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
  if f eq strtrim(string(iw),2) or strlowcase(f) eq wnames(iw) then $
  case ndim of
    1:result = w(*,iw)
    2:result = w(*,*,iw)
    3:result = w(*,*,*,iw)
endcase

if n_elements(result) gt 0 and rcut le 0 then return, sign*result

; Extract coordinate variables from xx:
case ndim of
    1:begin
        x=xx & y=0 & z=0
    end
    2:begin
        x=xx(*,*,0) & y=xx(*,*,1) & z=0
    end
    3:begin
        x=xx(*,*,*,0) & y=xx(*,*,*,1) & z=xx(*,*,*,2)
    end
endcase
r = sqrt(x^2+y^2+z^2)

; Return result after cutting out the body
if n_elements(result) gt 0 then begin
    loc = where(r le rcut, count) 
    if count gt 0 then result(loc)=min(result)
    return, sign*result
endif

; Extract primitive variables for calculating MHD type functions

; initialize all the variables as scalars 
rho=0 & ux=0 & uy=0 & uz=0 & bx=0 & by=0 & bz=0 & p=0 & e=0

; set the variables from the w
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
        'p'  : p=w(*,iw)
        'pth': p=w(*,iw)
        'mx' : ux=w(*,iw)/rho
        'my' : uy=w(*,iw)/rho
        'mz' : uz=w(*,iw)/rho
        'm1' : ux=w(*,iw)/rho
        'm2' : uy=w(*,iw)/rho
        'm3' : uz=w(*,iw)/rho
        'e'  : e=w(*,iw)
        else :
    endcase
    2: case wnames(iw) of
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
        'p'  : p=w(*,*,iw)
        'pth': p=w(*,*,iw)
        'mx' : ux=w(*,*,iw)/rho
        'my' : uy=w(*,*,iw)/rho
        'mz' : uz=w(*,*,iw)/rho
        'm1' : ux=w(*,*,iw)/rho
        'm2' : uy=w(*,*,iw)/rho
        'm3' : uz=w(*,*,iw)/rho
        'e'  : e=w(*,*,iw)
        else :
    endcase
    3: case wnames(iw) of
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
        'p'  : p=w(*,*,*,iw)
        'pth': p=w(*,*,*,iw)
        'mx' : ux=w(*,*,*,iw)/rho
        'my' : uy=w(*,*,*,iw)/rho
        'mz' : uz=w(*,*,*,iw)/rho
        'm1' : ux=w(*,*,*,iw)/rho
        'm2' : uy=w(*,*,*,iw)/rho
        'm3' : uz=w(*,*,*,iw)/rho
        'e'  : e=w(*,*,*,iw)
        else :
    endcase
endcase
                                ; Some extra variables
uu=ux^2+uy^2+uz^2
bb=bx^2+by^2+bz^2
u =sqrt(uu)
b =sqrt(bb)

                                ; Change energy into pressure
if n_elements(p) le 1 and n_elements(e) gt 1 then $
  p=(gamma-1)*(e-0.5*(rho*uu+bb))

                                ; Calculate gamma*p+bb if needed
if strmid(f,1,4) eq 'fast' or strmid(f,1,4) eq 'slow' then $
  cc=gamma*p+bb


;==== Put your function definition(s) below using the variables and eq. params
;     extracted from x, w and eqpar, and select cases by the function name "f"

case f of
    'Btheta'   : result=atan(by,sqrt(bx^2+bz^2))
                                ; momenta
    'mx'       : result=rho*ux
    'my'       : result=rho*uy
    'mz'       : result=rho*uz
                                ; Boris momenta
    'mBx'      : result=rho*ux + (bb*ux - (ux*bx+uy*by+uz*bz)*bx)/clight^2
    'mBy'      : result=rho*uy + (bb*uy - (ux*bx+uy*by+uz*bz)*by)/clight^2
    'mBz'      : result=rho*uz + (bb*uz - (ux*bx+uy*by+uz*bz)*bz)/clight^2
                                ; Electric field
    'Ex'       : result=(by*uz-uy*bz)
    'Ey'       : result=(bz*ux-uz*bx)
    'Ez'       : result=(bx*uy-ux*by)
                                ; Total energy
    'e'        : result=p/(gamma-1)+0.5*(rho*uu + bb)
                                ; Alfven speed and Mach number
    'calfvenx' : result=bx/sqrt(rho)
    'calfveny' : result=by/sqrt(rho)
    'calfvenz' : result=bz/sqrt(rho)
    'calfven'  : result=b /sqrt(rho)
    'Malfvenx' : result=ux/bx*sqrt(rho)
    'Malfveny' : result=uy/by*sqrt(rho)
    'Malfvenz' : result=uz/bz*sqrt(rho)
    'Malfven'  : result=u /b *sqrt(rho)
                                ; pressure, plasma beta, temperature, entropy
    'pbeta'    : result=2*p/bb
    's'        : result=p/rho^gamma
    'T'        : result=p/rho
    'T_SC'     : result=1.211E-8*p/rho ; [K] amu/kB*(dyne/cm^2) / (g/cm^3)
    'T_IH'     : result=1.211E-8*p/rho ; [K] amu/kB*(dyne/cm^2) / (g/cm^3)
    'T_GM'     : result=7.243E+7*p/rho ; [K] amu/kB*nPa / (amu/cm^3)
    'n_IH'     : result=rho/amu ; [/cc] n = rho/amu in CGS
    'n_SC'     : result=rho/amu ; [/cc] n = rho/amu in CGS
                                ; sound speed, Mach number
    'csound'   :result=sqrt(gamma*p/rho)
    'mach'  :result=u /sqrt(gamma*p/rho)
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
    'Mfast' :result=sqrt(rho*uu/cc)
    'Mfastx':result=ux/sqrt((cc+sqrt(cc^2-4*gamma*p*bx^2))/2/rho)
    'Mfasty':result=uy/sqrt((cc+sqrt(cc^2-4*gamma*p*by^2))/2/rho)
    'Mfastz':result=uz/sqrt((cc+sqrt(cc^2-4*gamma*p*bz^2))/2/rho)
    'Mslowx':result=ux/sqrt((cc-sqrt(cc^2-4*gamma*p*bx^2))/2/rho)
    'Mslowy':result=uy/sqrt((cc-sqrt(cc^2-4*gamma*p*by^2))/2/rho)
    'Mslowz':result=uz/sqrt((cc-sqrt(cc^2-4*gamma*p*bz^2))/2/rho)
    'divbxy':result=div(bx,by,x,y)
    'absdivbxy':result=abs(div(bx,by,x,y))
    'jz': case ndim of
        1: result=deriv(x,by)
        2: result=curl(bx,by,x,y)
        3: begin
            print,'Error in funcdef: jz is not implemented for 3D yet'
            retall
        end
    endcase

                                ; Magnetic vector potential A in slab symmetry
                                ; Density of contourlines is proportional to B
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

                                ; If "f" has not matched any function,
                                ; try evaluating it as an expression
    else:begin
        if sign eq -1 then begin
            f='-'+f
            sign=1
        endif

                                ; Create * or *,* or *,*,* for 
                                ; {var} --> w(*,*,iVar) replacements
        case nDim of
            1: Stars = '*,'
            2: Stars = '*,*,'
            3: Stars = '*,*,*,'
        endcase

                                ; replace {name} with appropriate variable
        while strpos(f,'{') ge 0 do begin
            iStart = strpos(f,'{')
            iEnd   = strpos(f,'}')
            Name   = strmid(f,iStart+1,iEnd-iStart-1)
            iVar   = where(variables eq Name)
            iVar   = iVar(0)
            if iVar lt 0 then begin
                print,'Error in funcdef: cannot find variable "{',Name, $
                  '}" among variables'
                print,variables
                retall
            endif
            fStart = strmid(f,0,iStart)
            fEnd   = strmid(f,iEnd+1,strlen(f)-iEnd-1)
            if iVar lt ndim then $
              f = fStart + 'xx(' + Stars + strtrim(iVar,2) + ')' + fEnd $
            else if iVar lt ndim+nw then $
              f = fStart + 'w(' + Stars + strtrim(iVar-nDim,2) + ')' + fEnd $
            else $
              f = fStart + 'eqpar(' + strtrim(iVar-nDim-nW,2) + ')' + fEnd

        endwhile


        if not execute('result='+f) then begin
            print,'Error in funcdef: cannot evaluate function=',func
            retall
        endif
    end
endcase

if n_elements(result) gt 0 then begin
                                ; exclude r < rcut
    loc=where(r le rcut, count) 
    if count gt 0 then result(loc)=min(result)

    return,sign*result
endif else begin
    print,'Error in funcdef: function=',func,' was not calculated ?!'
    retall
endelse

end
