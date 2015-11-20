;  Copyright (C) 2002 Regents of the University of Michigan, 
;  portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
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
;   Function names defined in the functiondef array or listed in the
;   case statement below are calculated and returned.
;
;   Expressions formed from the 
;   1. standard variables:  rho, mx, my, mz, ux, uy, uz, uu, u, p, bx, by, bz, bb, b
;   2. standard coordinates: x, y, z, r
;   3. standard equation parameters: gamma, rbody, c0
;   4. standard IDL functions and operators
;   5. {names of coordinates}, {names of variables}, {names of equation parameters}
;      Examples: {r}, {bx1}, {rbody} will be replaced with
;                xx(*,*,0), w(*,*,10), eqpar(2)
;   6. {names of functions defined in the functiondef array}
;      Examples: {T_GM}, {Mfast}
;
; Examples for valid strings: 
;   '3', 'rho', '{Mfast}+ux', '-T', 'rho*ux', '{bx1}^2+{by1}^2'...
;
; One can use "funcdef" interactively too, e.g.
; 
; ekin =funcdef(x,w,'(ux^2+uy^2)*rho,time,eqpar,variables)
; cfast=funcdef(x,w,'cfast',time,eqpar,variables)
;
;===========================================================================

  ;; Define various functions of the basic MHD variables
  ;; The functions names are evaluated in lower case
  functiondef = strlowcase(transpose([ $
                ['mxB'    , 'rho*ux+(bb*ux-(ux*bx+uy*by+uz*bz)*bx)/c0^2'], $ ; Boris momenta
                ['myB'    , 'rho*uy+(bb*uy-(ux*bx+uy*by+uz*bz)*by)/c0^2'], $
                ['mzB'    , 'rho*uz+(bb*uz-(ux*bx+uy*by+uz*bz)*bz)/c0^2'], $
                ['mx'       , 'rho*ux'                                  ], $ ; momenta
                ['my'       , 'rho*uy'                                  ], $
                ['mz'       , 'rho*uz'                                  ], $
                ['uH'       , 'uH0*sqrt({jx}^2+{jy}^2+{jz}^2)/rho'      ], $ ; Hall velocity
                ['uHx'      , 'uH0*{jx}/rho'                            ], $ 
                ['uHy'      , 'uH0*{jy}/rho'                            ], $ 
                ['uHz'      , 'uH0*{jz}/rho'                            ], $ 
                ['uex'      , 'ux-uH0*{jx}/rho'                         ], $ ; electron velocity
                ['uey'      , 'uy-uH0*{jy}/rho'                         ], $ 
                ['uez'      , 'uz-uH0*{jz}/rho'                         ], $ 
                ['ue'       , 'sqrt({uex}^2+{uey}^2+{uez}^2)'           ], $
                ['j'        , 'sqrt({jx}^2+{jy}^2+{jz}^2)'              ], $ ; current density
                ['divbxy'   , 'div(bx,by,x,y)'                          ], $ ; div(B) in 2D
                ['Ex'       , 'by*uz-uy*bz'                             ], $ ; electric field
                ['Ey'       , 'bz*ux-uz*bx'                             ], $
                ['Ez'       , 'bx*uy-ux*by'                             ], $
                ['e'        , 'p/(gamma-1)+0.5*(rho*uu + bb)'           ], $ ; energy density
                ['pbeta'    , '2*mu0*p/bb'                              ], $ ; plasma beta
                ['s'        , 'p/rho^gamma'                             ], $ ; entropy
                ['Ti'       , 'ti0*p/rho'                               ], $ ; ion temperature [K]
                ['Te'       , 'ti0*{pe}/rho'                            ], $ ; electron temp. [K]
                ['calfvenx' , 'bx/sqrt(rho*mu0A)'                       ], $ ; Alfven velocity
                ['calfveny' , 'by/sqrt(rho*mu0A)'                       ], $
                ['calfvenz' , 'bz/sqrt(rho*mu0A)'                       ], $
                ['calfven'  , 'b /sqrt(rho*mu0A)'                       ], $
                ['Malfvenx' , 'ux/bx*sqrt(rho*mu0A)'                    ], $ ; Alfven Mach number
                ['Malfveny' , 'uy/by*sqrt(rho*mu0A)'                    ], $
                ['Malfvenz' , 'uz/bz*sqrt(rho*mu0A)'                    ], $
                ['Malfven'  , 'u /b *sqrt(rho*mu0A)'                    ], $
                ['csound'   , 'sqrt(gs*p/rho)'                          ], $ ; ion sound speed
                ['csounde'  , 'sqrt(gs*pe/rho*mi/me)'                   ], $ ; electron sound speed
                ['mach'     , 'u /sqrt(gs*p/rho)'                       ], $ ; Mach number
                ['machx'    , 'ux/sqrt(gs*p/rho)'                       ], $
                ['machy'    , 'uy/sqrt(gs*p/rho)'                       ], $
                ['machz'    , 'uz/sqrt(gs*p/rho)'                       ], $
                ['cfast'    , 'sqrt(cc/rho)'                            ], $ ; fast magnetosonic speed
                ['cfastx'   , 'sqrt((cc+sqrt(cc^2-c4*p*bx^2))/2/rho)'   ], $
                ['cfasty'   , 'sqrt((cc+sqrt(cc^2-c4*p*by^2))/2/rho)'   ], $
                ['cfastz'   , 'sqrt((cc+sqrt(cc^2-c4*p*bz^2))/2/rho)'   ], $
                ['cslowx'   , 'sqrt((cc-sqrt(cc^2-c4*p*bx^2))/2/rho)'   ], $ ; slow speed
                ['cslowy'   , 'sqrt((cc-sqrt(cc^2-c4*p*by^2))/2/rho)'   ], $
                ['cslowz'   , 'sqrt((cc-sqrt(cc^2-c4*p*bz^2))/2/rho)'   ], $
                ['Mfast'    , 'sqrt(rho*uu/cc)'                         ], $ ; fast Mach number
                ['Mfastx'   , 'ux/sqrt((cc+sqrt(cc^2-c4*p*bx^2))/2/rho)'], $
                ['Mfasty'   , 'uy/sqrt((cc+sqrt(cc^2-c4*p*by^2))/2/rho)'], $
                ['Mfastz'   , 'uz/sqrt((cc+sqrt(cc^2-c4*p*bz^2))/2/rho)'], $
                ['Mslowx'   , 'ux/sqrt((cc-sqrt(cc^2-c4*p*bx^2))/2/rho)'], $ ; slow Mach number
                ['Mslowy'   , 'uy/sqrt((cc-sqrt(cc^2-c4*p*by^2))/2/rho)'], $
                ['Mslowz'   , 'uz/sqrt((cc-sqrt(cc^2-c4*p*bz^2))/2/rho)'], $
                ['uth'      , 'sqrt(cs0*p/rho)'                         ], $ ; ion thermal speed
                ['uthe'     , 'sqrt(cs0*{pe}/rho*mi/me)'                ], $ ; electron thermal speed
                ['omegapi'  , 'op0*sqrt(rho)'                           ], $ ; ion plasma frequency
                ['omegape'  , 'op0*sqrt(rho)*mi/me'                     ], $ ; electron plasma freq.
                ['omegaci'  , 'oc0*b'                                   ], $ ; ion gyro frequency
                ['omegace'  , 'oc0*b*mi/me'                             ], $ ; electron gyro freq.
                ['rgyro'    , 'rg0*sqrt(p/rho)/(b>1e-30)'               ], $ ; gyro radius  
                ['rgSI'     , 'rg0*sqrt(p/rho)/(b>1e-30)*xSI'           ], $ ; gyro radius in SI
                ['rgyroe'   , 'rg0*sqrt(p/rho*me/mi)/(b>1e-30)'         ], $ ; electron gyro radius  
                ['rgeSI'    , 'rg0*sqrt(p/rho*me/mi)/(b>1e-30)*xSI'     ], $ ; electron gyro radius in SI
                ['dinertial', 'di0/sqrt(rho)'                           ], $ ; inertial length
                ['diSI'     ,' di0/sqrt(rho)*xSI'                       ], $ ; ion inertial length in SI
                ['skindepth',' di0/sqrt(rho)*me/mi'                     ], $ ; electron skin depth
                ['deSI'     ,' di0/sqrt(rho)*me/mi*xSI'                 ], $ ; electron skin depth in SI
                ['ldebye'   , 'ld0/c0*sqrt(p)/rho'                      ], $ ; Debye length
                ['ldSI'     , 'ld0/c0*sqrt(p)/rho*xSI'                  ]  $ ; Debye length in SI
                                     ]))

  common phys_units, $
     fixunits, typeunit, xSI, tSI, rhoSI, uSI, pSI, bSI, jSI, Mi, Me
  common phys_convert, ti0, cs0, mu0A, mu0, c0, uH0, op0, oc0, rg0, di0, ld0
  common phys_const, kbSI, mpSI, mu0SI, eSI, ReSI, RsSI, cSI, e0SI

  if n_elements(xx) eq 0 or n_elements(w) eq 0 then begin
     print,'ERROR in funcdef: xx or w are not defined'
     help,xx,w
     retall
  endif

  if n_elements(rcut) eq 0 then rcut = -1

  ;; In 1D xx(n1), in 2D xx(n1,n2,2), in 3D xx(n1,n2,n3,3)
  siz=size(xx)
  ndim=siz(0)-1
  if ndim eq 0 then ndim=1
  n1=siz(1)
  if ndim gt 1 then n2=siz(2)
  if ndim gt 2 then n3=siz(3)

  ;; For 1 variable: w(n1), w(n1*n2), w(n1,n2),  w(n1,n2,n3)
  ;; for more      : w(n1,nw),w(n1,n2,nw),w(n1,n2,n3,nw)
  siz=size(w)
  if siz(0) le ndim then nw=1 else nw=siz(ndim+1)

  ;; Number of equation parameters
  nEqpar = n_elements(eqpar)

  ;; Extract equation parameters
  gamma  =  5./3.
  c0     =  1.0
  rbody  = -1.0

  if nEqpar gt 0 then begin
     for i=nDim+nW,n_elements(variables)-1 do begin
        iEqpar = i-nDim-nW
        case variables(i) of
           'g'     : gamma  = eqpar(iEqpar)
           'gamma' : gamma  = eqpar(iEqpar)
           'c'     : c0     = eqpar(iEqpar)
           'r'     : rbody  = eqpar(iEqpar)
           'rbody' : rbody  = eqpar(iEqpar)
           else:
        endcase
     endfor
  endif

  ;; Define gamma for the sound speed = sqrt(gs*p/rho) with units
  gs = gamma*cs0

  ;; Variable names
  if n_elements(variables) eq 0 then variables=strarr(ndim+nw+neqpar)
  wnames = strlowcase(variables(ndim:ndim+nw-1))

  ;; Check for a negative sign in func
  if strmid(func,0,1) eq '-' then begin 
     f=strmid(func,1,strlen(func)-1)
     sign=-1
  endif else begin
     f=func
     sign=1
  endelse

  ;; Check if f is among the variable names listed in wnames or if it is a number
  for iw=0,nw-1 do $
     if f eq strtrim(string(iw),2) or strlowcase(f) eq wnames(iw) then $
        case ndim of
     1:result = w(*,iw)
     2:result = w(*,*,iw)
     3:result = w(*,*,*,iw)
  endcase

  if n_elements(result) gt 0 and rcut le 0 then return, sign*result

  ;; set radial distance (assuming cartesian coordinates)
  case ndim of
     1: r = abs(xx)
     2: r = sqrt(xx(*,*,0)^2 + xx(*,*,1)^2)
     3: r = sqrt(xx(*,*,*,0)^2 + xx(*,*,*,1)^2 + xx(*,*,*,2)^2)
  endcase

  if n_elements(result) gt 0 then begin
                                ; Return result after cutting out at rcut 
                                ; set value to the minimum of the remaining values
     loc = where(r le rcut, count)
     loc1= where(r gt rcut)
     if count gt 0 then result(loc) = min(result(loc1))
     return, sign*result
  endif

  ;; set the coordinate arrays x, y, z if they occur in the variable names
  x = 0 & y = 0 & z = 0
  for idim = 0, ndim-1 do case ndim of
     1: case variables(idim) of
        'x': x = xx
        'y': y = xx
        'z': z = xx
        else: x = xx
     endcase
     2: case variables(idim) of
        'x': x = xx(*,*,idim)
        'y': y = xx(*,*,idim)
        'z': z = xx(*,*,idim)
        else: case idim of
           0: x = xx(*,*,idim)
           1: y = xx(*,*,idim)
        endcase
     end
     3: case variables(idim) of
        'x': x = xx(*,*,*,idim)
        'y': y = xx(*,*,*,idim)
        'z': z = xx(*,*,*,idim)
        else: case idim of
           0: x = xx(*,*,*,idim)
           1: y = xx(*,*,*,idim)
           2: z = xx(*,*,*,idim)
        endcase
     end
  endcase

  ;; Extract primitive variables for calculating MHD type functions

  ;; initialize all the variables as scalars 
  rho=0 & ux=0 & uy=0 & uz=0 & bx=0 & by=0 & bz=0 & p=0 & e=0

  ;; set the variables from the w
  for iw = 0, nw-1 do case ndim of
     1: case wnames(iw) of
        'rho': rho=w(*,iw)
        'ux' : ux=w(*,iw)
        'uy' : uy=w(*,iw)
        'uz' : uz=w(*,iw)
        'bx' : bx=w(*,iw)
        'by' : by=w(*,iw)
        'bz' : bz=w(*,iw)
        'p'  : p=w(*,iw)
        'pth': p=w(*,iw)
        'mx' : ux=w(*,iw)/rho
        'my' : uy=w(*,iw)/rho
        'mz' : uz=w(*,iw)/rho
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
        'p'  : p=w(*,*,iw)
        'pth': p=w(*,*,iw)
        'mx' : ux=w(*,*,iw)/rho
        'my' : uy=w(*,*,iw)/rho
        'mz' : uz=w(*,*,iw)/rho
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
        'p'  : p=w(*,*,*,iw)
        'pth': p=w(*,*,*,iw)
        'mx' : ux=w(*,*,*,iw)/rho
        'my' : uy=w(*,*,*,iw)/rho
        'mz' : uz=w(*,*,*,iw)/rho
        'e'  : e=w(*,*,*,iw)
        else :
     endcase
  endcase

  ;; Extra variables
  uu = ux^2 + uy^2 + uz^2       ; velocity squared, 
  bb = bx^2 + by^2 + bz^2       ; magnetic field squared
  u  = sqrt(uu)                 ; speed
  b  = sqrt(bb)                 ; magnetic field strength

  ;; Change energy into pressure
  if n_elements(p) le 1 and n_elements(e) gt 1 then $
     p = (gamma-1)*(e - 0.5*(rho*uu + bb))

  ;; Calculate gamma*p+bb if needed
  if strpos(f,'fast') ge 0 or strpos(f, 'slow') ge 0 then begin
     c4 = 4*gs/mu0A
     cc = gs*p + bb/mu0A
  end

  ;; Add functions to the basic variable list
  functions = [strlowcase(variables), functiondef(*,0)]

  ;;==== Functions that cannot be expressed from the basic variables
  case f of
     'Btheta'   : result=atan(by,sqrt(bx^2+bz^2))
                                ; sound speed, Mach number
     'jz': case ndim of
        1: result=deriv(x,by)
        2: result=curl(bx,by,x,y)
        3: begin
           print,'Error in funcdef: jz is not implemented for 3D yet'
           retall
        end
     endcase

     ;; Magnetic vector potential A in slab symmetry
     ;; Density of contourlines is proportional to B
     'Ax': begin
        result=dblarr(n1,n2)
                                ; Integrate along the first row
        for i1=1,n1-1 do result(i1,0)=result(i1-1,0) $
           +(bz(i1,0)+bz(i1-1,0))*(y(i1,0)-y(i1-1,0))*0.5 $
           -(by(i1,0)+by(i1-1,0))*(z(i1,0)-z(i1-1,0))*0.5
                                ; Integrate all columns vertically
        for i2=1,n2-1 do result(*,i2)=result(*,i2-1) $
           +(bz(*,i2)+bz(*,i2-1))*(y(*,i2)-y(*,i2-1))*0.5 $
           -(by(*,i2)+by(*,i2-1))*(z(*,i2)-z(*,i2-1))*0.5
     end
     'Ay': begin
        result=dblarr(n1,n2)
                                ; Integrate along the first row
        for i1=1,n1-1 do result(i1,0)=result(i1-1,0) $
           -(bz(i1,0)+bz(i1-1,0))*(x(i1,0)-x(i1-1,0))*0.5 $
           +(bx(i1,0)+bx(i1-1,0))*(z(i1,0)-z(i1-1,0))*0.5
                                ; Integrate all columns vertically
        for i2=1,n2-1 do result(*,i2)=result(*,i2-1) $
           -(bz(*,i2)+bz(*,i2-1))*(x(*,i2)-x(*,i2-1))*0.5 $
           +(bx(*,i2)+bx(*,i2-1))*(z(*,i2)-z(*,i2-1))*0.5
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
                                ; check if f is the name of a function
        iFunc = where(functiondef(*,0) eq STRLOWCASE(f))
        iFunc = iFunc(0)

        if iFunc ge 0 then f = functiondef(iFunc,1)

        if sign eq -1 then begin
           f='-'+f
           sign=1
        endif

        ;; Create * or *,* or *,*,* for {var} --> w(*,*,iVar) replacements
        case nDim of
           1: Stars = '*,'
           2: Stars = '*,*,'
           3: Stars = '*,*,*,'
        endcase

        ;; replace {name} with appropriate variable
        while strpos(f,'{') ge 0 do begin
           iStart = strpos(f,'{')
           iEnd   = strpos(f,'}')
           Name   = strmid(f,iStart+1,iEnd-iStart-1)
           iVar   = where(functions eq STRLOWCASE(Name))
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
           else if iVar lt ndim+nw+neqpar then $
              f = fStart + 'eqpar(' + strtrim(iVar-nDim-nW,2) + ')' + fEnd $
           else $
              f = fStart + '(' + functiondef(iVar-nDim-nW-nEqpar,1) + ')' + fEnd

        endwhile

        if not execute('result='+f) then begin
           print,'Error in funcdef: cannot evaluate function=',func
           retall
        endif

     end


  endcase

  if n_elements(result) gt 0 then begin
     
     if rcut gt 0 then begin
                                ; exclude r < rcut
        loc = where(r le rcut, count) 
        loc1= where(r gt rcut)
        if count gt 0 then result(loc) = min(result(loc1))
     endif
     return,sign*result
  endif else begin
     print,'Error in funcdef: function=',func,' was not calculated ?!'
     retall
  endelse

end
