;^CFG COPYRIGHT UM
;------------------------------------------------------------------------
;------------------------------------------------------------------------

;
; function mm
;
;  computes min and max of an array and returns them in an array
;

function mm, array
return, [min(array),max(array)]
end

;------------------------------------------------------------------------
;------------------------------------------------------------------------

function tostr,value
  return, strcompress(string(fix(value)),/remove_all)
end

;------------------------------------------------------------------------
;------------------------------------------------------------------------

;*****************************************************************************
 
pro readct, ncolors, ctname
 
;*****************************************************************************
 
  openr,11,ctname
  ncolors=0
  readf,11,ncolors
  color = fltarr(3,255)
  readf,11,color
  close,11

  tvlct,color(0,*),color(1,*),color(2,*)
 
  return
 
end

;------------------------------------------------------------------------
;------------------------------------------------------------------------

pro plotdumb

  plot, [0,1], 			$
	xstyle = 5,		$
	ystyle = 5,		$
	pos = [0,0,1,1],	$
	/nodata

  return

end

;------------------------------------------------------------------------
;------------------------------------------------------------------------


function ask, what, orig_ans, set_orig = set_orig

  if n_elements(orig_ans) eq 0 then orig_ans = ''

  answer = ''

  read, 'Enter '+what+' ['+orig_ans+'] : ', answer

  if strlen(answer) eq 0 then answer = orig_ans

  if n_elements(set_orig) gt 0 then orig_ans = answer

  return, answer

end

;------------------------------------------------------------------------
;------------------------------------------------------------------------

PRO c_a_to_r, timearray, timereal

  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]
  if ((timearray(0) mod 4) eq 0) then dayofmon(1) = dayofmon(1) + 1

  timereal = double(0.0)

  if timearray(0) lt 65 then timearray(0) = timearray(0) + 2000
  if timearray(0) gt 1900 then numofyears = timearray(0)-1965 		      $
  else numofyears = timearray(0)-65	
  numofleap = floor(float(numofyears)/4.0)
  numofmonths = timearray(1) - 1
  numofdays = 0

  for i = 0, numofmonths-1 do begin

    numofdays = numofdays + dayofmon(i)

  endfor

  numofdays = numofdays + timearray(2) - 1
  numofhours = timearray(3)
  numofminutes = timearray(4)
  numofseconds = timearray(5)

  timereal = double(numofseconds*1.0) +       $
	     double(numofminutes*60.0) +             $
	     double(numofhours*60.0*60.0) +          $
	     double(numofdays*24.0*60.0*60.0) +      $
	     double(numofleap*24.0*60.0*60.0) +      $
	     double(numofyears*365.0*24.0*60.0*60.0)

  RETURN

END

;------------------------------------------------------------------------
;------------------------------------------------------------------------

pro c_a_to_s, timearray, strtime

  mon='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC' 

  sd = '0'+tostr(timearray(2))
  sd = strmid(sd,strlen(sd)-2,2)
  sm = strmid(mon,(timearray(1)-1)*3,3)
  if timearray(0) lt 1900 then year = timearray(0) 		$
  else year = timearray(0)-1900
  if (year ge 100) then year = year - 100
  sy = chopr('0'+tostr(year),2)
  sh = '0'+tostr(timearray(3))
  sh = strmid(sh,strlen(sh)-2,2)
  si = '0'+tostr(timearray(4))
  si = strmid(si,strlen(si)-2,2)
  ss = '0'+tostr(timearray(5))
  ss = strmid(ss,strlen(ss)-2,2)

  strtime = sd+'-'+sm+'-'+sy+' '+sh+':'+si+':'+ss+'.000'

  RETURN

END

;------------------------------------------------------------------------
;------------------------------------------------------------------------

PRO c_r_to_a, timearray, timereal

  dayofmon = [31,28,31,30,31,30,31,31,30,31,30,31]

  timearray = intarr(6)

  speryear = double(31536000.0)
  sperday  = double(86400.0)
  sperhour = double(3600.0)
  spermin  = double(60.0)

  numofyears = floor(timereal/speryear)
  if (numofyears+65) mod 4 eq 0 then dayofmon(1) = dayofmon(1) + 1
  numofdays = floor((timereal mod speryear)/sperday)
  numofleap = floor(numofyears / 4)
  numofdays = numofdays - numofleap
  if numofdays lt 0 then begin
    if (numofyears+65) mod 4 eq 0 then dayofmon(1) = dayofmon(1) - 1
    numofyears = numofyears - 1
    numofdays = numofdays + numofleap + 365
    if (numofyears+65) mod 4 eq 0 then dayofmon(1) = dayofmon(1) + 1
    numofleap = floor(numofyears / 4)
    numofdays = numofdays - numofleap
  endif
  numofhours = floor((timereal mod sperday)/sperhour)
  numofmin = floor((timereal mod sperhour)/spermin)
  numofsec = floor(timereal mod spermin)

  numofmon = 0

  while numofdays ge dayofmon(numofmon) do begin

    numofdays = numofdays - dayofmon(numofmon)
    numofmon = numofmon + 1

  endwhile

  timearray(0) = (numofyears + 65) mod 100
  timearray(1) = numofmon + 1
  timearray(2) = numofdays + 1
  timearray(3) = numofhours
  timearray(4) = numofmin
  timearray(5) = numofsec

  RETURN

END

;------------------------------------------------------------------------
;------------------------------------------------------------------------

pro c_s_to_a, timearray, strtime

  timearray = intarr(6)

  mon='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'

  timearray(2) = fix(strmid(strtime,0,2))
  smon = strmid(strtime,3,3)
  bs = byte(smon)
  loc = where((bs ge 97) and (bs lt 122), count)
  if count gt 0 then bs(loc) = bs(loc)-byte(32)
  smon = string(bs)

  for j=0,11 do 							      $
    if strmid(mon,j*3,3) eq smon then timearray(1)=j+1
  timearray(0)=fix(strmid(strtime,7,2))+1900

  if (strlen(strmid(strtime,10,2)) gt 0) and			$
     (strmid(strtime,10,2) ne '  ') then			$
    timearray(3)=fix(strmid(strtime,10,2))			$
  else timearray(3) = 0
  if (strlen(strmid(strtime,13,2)) gt 0) and			$
     (strmid(strtime,13,2) ne '  ') then			$
    timearray(4)=fix(strmid(strtime,13,2))			$
  else timearray(4) = 0
  if (strlen(strmid(strtime,16,6)) gt 0) and			$
     (strmid(strtime,16,2) ne '  ') then			$
    timearray(5)=fix(strmid(strtime,16,6))			$
  else timearray(5) = 0

  RETURN

END

;------------------------------------------------------------------------
;------------------------------------------------------------------------

function chopr, svalue, n
  if strlen(svalue) lt n then n = strlen(svalue)
  return, strmid(svalue, strlen(svalue)-n,n)
end

;------------------------------------------------------------------------
;------------------------------------------------------------------------

pro closedevice

  if !d.name eq 'PS' then begin
    device, /close
    set_plot, 'X'
  endif

  return

end

;------------------------------------------------------------------------
;------------------------------------------------------------------------


;
; get_position
;
; used in conjunction with pos_space. Determines the position of the current
; plotting region, given the output parameters from pos_space.
;
; Input parameters:
; nb, space, bs, nbx, nby, xoff, yoff, xf, yf - Outputs from pos_space
; pos_num - the number of the plot, ranges from 0 : bs-1
;
; Output parameters:
;
; pos - the position of the plot, used in the plot command
;
; modified to make rectangles on Jan 2, 1998

pro get_position, nb, space, sizes, pos_num, pos, rect = rect,		$
		  xmargin = xmargin, ymargin = ymargin

  xipos = fix(pos_num) mod sizes.nbx
  yipos = fix(pos_num)/sizes.nbx

  yf2 = sizes.yf
  yf = sizes.yf*(1.0-space)
  xf2 = sizes.xf
  xf = sizes.xf*(1.0-space)

  if n_elements(rect) gt 0 then begin

    if n_elements(xmargin) gt 0 then xmar = xmargin(0) 			$
    else xmar = space/2.0

    if n_elements(ymargin) gt 0 then ymar = ymargin(0) 			$
    else ymar = space/2.0

    xtotal = 1.0 - (space*float(sizes.nbx-1) + xmar + xf2*space/2.0)
    xbs = xtotal/(float(sizes.nbx)*xf)

    xoff = xmar - xf2*space/2.0

    ytotal = 1.0 - (space*float(sizes.nby-1) + ymar + yf2*space/2.0)
    ybs = ytotal/(float(sizes.nby)*yf)

    yoff = 0.0

  endif else begin

    xbs  = sizes.bs
    xoff = sizes.xoff
    ybs  = sizes.bs
    yoff = sizes.yoff

  endelse

  xpos0 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0
  xpos1 = float(xipos) * (xbs+space)*xf + xoff + xf2*space/2.0 + xbs*xf

  ypos0 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf + ybs*yf) - yoff
  ypos1 = (1.0-yf2*space/2) - (yipos * (ybs+space)*yf) - yoff

  pos= [xpos0,ypos0,xpos1,ypos1]

  RETURN

END

;------------------------------------------------------------------------
;------------------------------------------------------------------------

function mklower, string

  temp = byte(string)

  loc = where((temp ge 65) and (temp le 90), count)

  if count ne 0 then temp(loc) = temp(loc)+32

  return, string(temp)

end

;------------------------------------------------------------------------
;------------------------------------------------------------------------

pro plotmlt, maxran, white = white, black = black, 		$
      no00 = no00, no06 = no06, no12 = no12, no18 = no18,	$
      longs = longs

  if n_elements(white) gt 0 then color = 255
  if n_elements(black) gt 0 then color = 0
  if n_elements(color) eq 0 then begin
    if !d.name eq 'PS' then color = 0 else color = 255
  endif

  if n_elements(no00) eq 0 then no00 = 0
  if n_elements(no06) eq 0 then no06 = 0
  if n_elements(no12) eq 0 then no12 = 0
  if n_elements(no18) eq 0 then no18 = 0

  if n_elements(longs) eq 0 then begin
    p00 = '00'
    p06 = '06'
    p12 = '12'
    p18 = '18'
  endif else begin
    p00 = '00'
    p06 = '90'
    p12 = '180'
    p18 = '270'
  endelse

  t = findgen(182.0)*2.0*!pi/180.0
  xp = cos(t)
  yp = sin(t)

  plots, maxran*xp, maxran*yp, color = color
  for i=10,maxran, 10 do					$
    oplot, float(i)*xp, float(i)*yp,linestyle=1, color = color

  oplot, [-maxran,maxran],[0.0,0.0], linestyle =1, color = color
  oplot, [0.0,0.0], [-maxran,maxran], linestyle = 1, color = color

  xs  = float(!d.x_size)
  ys  = float(!d.y_size)
  pos = float(!p.clip(0:3))
  pos([0,2]) = pos([0,2])/xs
  pos([1,3]) = pos([1,3])/ys

  mid_x = (pos(2) + pos(0))/2.0
  mid_y = (pos(3) + pos(1))/2.0

  ch = 0.8

  y_ch = ch*float(!d.y_ch_size)/ys
  x_ch = ch*float(!d.x_ch_size)/xs

  if no00 eq 0 then 							$
    xyouts, mid_x, pos(1) - y_ch*1.1, p00, alignment=0.5, 		$
      charsize=ch, /norm

  if no06 eq 0 then 							$
    xyouts, pos(2)+x_ch*0.15, mid_y - y_ch/2.0, p06, 			$
      charsize=ch, /norm

  if no12 eq 0 then 							$
    xyouts, mid_x, pos(3) + y_ch*0.15, p12, alignment=0.5, 		$
      charsize=ch, /norm

  if no18 eq 0 then 							$
    xyouts, pos(0)-x_ch*0.15, mid_y - y_ch/2.0, p18, 			$
      charsize=ch, /norm, alignment = 1.0

  return

end

;------------------------------------------------------------------------
;------------------------------------------------------------------------

;
; pos_space
;
; Determines the size and multiplication factors for plotting perfect circles
; or squares. This routine is used to simply find the parameters, another
; procedure, get_position, is used to actually find the position of the 
; circle or square.
; This routine maxamizes the area used by the plots, determinining the best
; positions for the number of plots that the user has selected.
;
; input parameters:
; nb - number of plots on a page
; space - amount of space in between each of the plots in normalized
;	  coordinates
;
; output parameters:
; bs - box size (size of the plotting region)
; nbx, nby - number of plots in the x and y directions
; xoff, yoff - x and y offsets for positions
; xf, yf - x and y multiplication factors for making perfect squares
;
; This has been adapted to allow the user to define how many objects
;   are in the x and y direction on Jan 2, 1998

pro pos_space, nb, space, sizes, nx = nx, ny = ny

  sizes = {bs:0.0, nbx:0, nby:0, xoff:0.0, yoff:0.0, xf:0.0, yf:0.0}

  xsi = float(!d.x_size)
  ysi = float(!d.y_size)

  xs = xsi - 5.0*space*xsi
  ys = ysi - 5.0*space*ysi

  if nb eq 1 then begin

    sizes.nbx = 1
    sizes.nby = 1
    sizes.bs = 1.0 - space

    if xs gt ys then begin

       sizes.yf = 1.0
       sizes.xf = ys/xs

    endif else begin

       sizes.xf = 1.0
       sizes.yf = xs/ys

     endelse

  endif else begin

    if (n_elements(nx) gt 0) then begin
      sizes.nbx = nx(0)
      if n_elements(ny) eq 0 then sizes.nby = nb/nx(0) else sizes.nby = ny(0)
    endif else begin
      if (n_elements(ny) gt 0) then begin
        sizes.nby = ny(0)
        sizes.nbx = nb/ny(0)
      endif else begin
        if xs gt ys then begin
          sizes.nbx = round(sqrt(nb))
          sizes.nby = fix(nb/sizes.nbx)
        endif else begin
          sizes.nby = round(sqrt(nb))
          sizes.nbx = fix(nb/sizes.nby)
        endelse
      endelse
    endelse

    if xs gt ys then begin

      if (sizes.nbx*sizes.nby lt nb) then 				$
	if (sizes.nbx le sizes.nby) then sizes.nbx = sizes.nbx + 1	$
	else sizes.nby = sizes.nby + 1					$
      else								$
	if (sizes.nbx lt sizes.nby) and					$
	   (n_elements(nx) eq 0) and					$
	   (n_elements(ny) eq 0) then begin
	  temp = sizes.nby
	  sizes.nby = sizes.nbx
	  sizes.nbx = temp
	endif

      sizes.yf = 1.0
      sizes.xf = ys/xs
      sizes.bs = ((1.0-space*(sizes.nbx-1))/sizes.nbx )/sizes.xf
      if sizes.nby*sizes.bs+space*(sizes.nby-1) gt 1.0 then 		$
	sizes.bs = (1.0- space*(sizes.nby-1))/sizes.nby 

    endif else begin

      if (sizes.nbx*sizes.nby lt nb) then				$
	if (sizes.nby le sizes.nbx) then sizes.nby = sizes.nby + 1	$
	else sizes.nbx = sizes.nbx + 1					$
      else								$
	if (sizes.nby lt sizes.nbx) and					$
	   (n_elements(nx) eq 0) and					$
	   (n_elements(ny) eq 0) then begin
	  temp = sizes.nby
	  sizes.nby = sizes.nbx
	  sizes.nbx = temp
	endif

      sizes.xf = 1.0
      sizes.yf = xs/ys
      sizes.bs = ((1.0 - space*(sizes.nby-1))/sizes.nby)/sizes.yf
      if sizes.nbx*sizes.bs+space*(sizes.nbx-1) gt 1.0 then 		$
	sizes.bs = (1.0 - space*(sizes.nbx-1))/sizes.nbx

    endelse

  endelse

  sizes.xoff = (1.0 - sizes.xf*(sizes.bs*sizes.nbx + space*(sizes.nbx-1)))/2.0
  sizes.yoff = (1.0 - sizes.yf*(sizes.bs*sizes.nby + space*(sizes.nby-1)))/2.0

  RETURN

END

;------------------------------------------------------------------------
;------------------------------------------------------------------------

pro setdevice, psfile, orient, psfont, percent, eps=eps, 	$
	       psname_inq = psname_inq

  if n_elements(psfile) eq 0 then begin

    psfile = ''
    if n_elements(psname_inq) gt 0 then begin
      read, 'Enter ps filename : ',psfile
    endif
    if strlen(psfile) eq 0 then psfile = 'idl.ps'

  endif

  if n_elements(percent) eq 0 then percent = 1.0		$
  else if percent gt 1.0 then percent = float(percent)/100.0
  if n_elements(orient) eq 0 then orient = 'landscape'
  if n_elements(psfont) eq 0 then psfont = 28
  if n_elements(eps) eq 0 then eps = 0 else eps = 1
  set_plot, 'ps', /copy, /interpolate

  !p.font = 0

  if (strmid(orient,0,1) eq 'p') or (strmid(orient,0,1) eq 'P') then begin

    changep = percent
    xs = 7.5
    ys = 9.5

    if eps eq 0 then begin

      case (psfont) of

	0  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Courier 
	1  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Courier, /Bold 
    	2  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Courier, /Oblique 
	3  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Courier, /Bold, /Oblique
       	4  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Helvetica
      	5  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Helvetica, /Bold
    	6  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Helvetica, /Oblique
       	8  : device, file = psfile, /color, bits=8,      	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Helvetica, /Bold, /Oblique 
    	12 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Avantgarde, /Book 
     	13 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Avantgarde, /Book, /Oblique
	 14 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Avantgarde, /Demi 
      	15 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Avantgarde, /Demi, /Oblique
       	20 : device, file = psfile, /color, bits=8, $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Schoolbook
   	21 : device, file = psfile, /color, bits=8,$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Schoolbook, /Bold
      	22 : device, file = psfile, /color, bits=8,$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Schoolbook, /Italic
       	23 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Schoolbook, /Bold, /Italic 
	28 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Times
	29 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Times, /Bold
	 30 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Times, /Italic
	31 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		xoff = (8.5-xs*changep)/2.0, yoff = (11.0-ys*changep)/2.0,  $
		/Times, /Bold, /Italic

      endcase

    endif else begin

      case (psfont) of

	0  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Courier, /encapsulated 
	1  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Courier, /Bold, /encapsulated 
    	2  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Courier, /Oblique, /encapsulated 
	3  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Courier, /Bold, /Oblique, /encapsulated
       	4  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Helvetica, /encapsulated
      	5  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Helvetica, /Bold, /encapsulated
    	6  : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Helvetica, /Oblique, /encapsulated
       	8  : device, file = psfile, /color, bits=8,      	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Helvetica, /Bold, /Oblique, /encapsulated 
    	12 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Avantgarde, /Book, /encapsulated 
     	13 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Avantgarde, /Book, /Oblique, /encapsulated
	 14 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Avantgarde, /Demi, /encapsulated 
      	15 : device, file = psfile, /color, bits=8,	      $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Avantgarde, /Demi, /Oblique, /encapsulated
       	20 : device, file = psfile, /color, bits=8, $
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Schoolbook, /encapsulated
   	21 : device, file = psfile, /color, bits=8,$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Schoolbook, /Bold, /encapsulated
      	22 : device, file = psfile, /color, bits=8,$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Schoolbook, /Italic, /encapsulated
       	23 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Schoolbook, /Bold, /Italic, /encapsulated 
	28 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Times, /encapsulated
	29 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Times, /Bold, /encapsulated
	 30 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Times, /Italic, /encapsulated
	31 : device, file = psfile, /color, bits=8,	$
		/inches, /portrait, xsize = xs*changep, ysize = ys*changep, $
		/Times, /Bold, /Italic, /encapsulated 

      endcase

    endelse

  endif else begin

    xs = 10.0
    ys = 7.0
    change = percent

    if eps eq 0 then begin

      case (psfont) of

	0  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Courier 
	1  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Courier, /Bold 
    	2  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Courier, /Oblique 
	3  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Courier, /Bold, /Oblique
       	4  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Helvetica
      	5  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Helvetica, /Bold
    	6  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Helvetica, /Oblique
       	8  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Helvetica, /Bold, /Oblique 
    	12 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Avantgarde, /Book 
     	13 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Avantgarde, /Book, /Oblique
	14 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
	 	/Avantgarde, /Demi 
      	15 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Avantgarde, /Demi, /Oblique
       	20 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Schoolbook
   	21 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Schoolbook, /Bold
      	22 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Schoolbook, /Italic
       	23 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Schoolbook, /Bold, /Italic 
	28 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Times
	29 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Times, /Bold
	30 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Times, /Italic
	31 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		yoff=11.0-(11.0-xs*change)/2.0, 		 $
		xoff=(8.5-ys*change)/2.0,			 $
		/inches,					 $
		/Times, /Bold, /Italic

      endcase

    endif else begin

      case (psfont) of

	0  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Courier, /encapsulated  
	1  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Courier, /Bold, /encapsulated  
    	2  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Courier, /Oblique, /encapsulated  
	3  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Courier, /Bold, /Oblique, /encapsulated 
       	4  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Helvetica, /encapsulated 
      	5  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Helvetica, /Bold, /encapsulated 
    	6  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Helvetica, /Oblique, /encapsulated 
       	8  : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Helvetica, /Bold, /Oblique, /encapsulated  
    	12 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Avantgarde, /Book, /encapsulated  
     	13 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Avantgarde, /Book, /Oblique, /encapsulated 
	14 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
	 	/Avantgarde, /Demi, /encapsulated  
      	15 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Avantgarde, /Demi, /Oblique, /encapsulated 
       	20 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Schoolbook, /encapsulated 
   	21 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Schoolbook, /Bold, /encapsulated 
      	22 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Schoolbook, /Italic, /encapsulated 
       	23 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Schoolbook, /Bold, /Italic, /encapsulated  
	28 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Times, /encapsulated 
	29 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Times, /Bold, /encapsulated 
	30 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Times, /Italic, /encapsulated 
	31 : device, file = psfile, /color, bits=8, /landscape,  $
		xsize=xs*change, ysize=ys*change,		 $
		/inches,					 $
		/Times, /Bold, /Italic, /encapsulated 

      endcase

    endelse

  endelse

  return

end

;------------------------------------------------------------------------
;------------------------------------------------------------------------

list = findfile("-t i*.idl")
if strlen(list(0)) gt 0 then filein = list(0) $
else filein = 'ionosphere_n000000.dat'

filein = ask('filename',filein)

list = findfile(filein)
nfiles = n_elements(list)
if nfiles eq 1 and strlen(list(0)) eq 0 then begin
  print, "I can't seem to find that file."
  stop
endif

line = ''

mr = float(ask('maximum range','40'))

time = dblarr(nfiles)
itime = intarr(6)

for n = 0, nfiles-1 do begin

  filein = list(n)

  if strpos(filein,'save') gt -1 then begin
    restore, filein
  endif else begin

    openr,1,filein

    done = 0

    while (not done) do begin

      readf,1, line

      if (strpos(mklower(line),"numerical") gt -1) then begin

        readf,1, nvars
        readf,1, nlats
        readf,1, nlons

        tmp = fltarr(nvars)
        data = fltarr(2,nvars,nlons,nlats)

        vars = strarr(nvars)

      endif

      if (strpos(mklower(line),"variable") gt -1) then begin

        for i=0,nvars-1 do begin
          readf,1,line
          vars(i) = strmid(line,6,strlen(line)-6)
	  if n eq 0 then print, chopr(' '+tostr(i),2),'. ',vars(i)
        endfor

      endif

      if (strpos(mklower(line),"time") gt -1) then begin

        int_tmp = 0
        for i=0,5 do begin
          readf, 1, int_tmp
          itime(i) = int_tmp
        endfor

        c_a_to_r, itime, rtime
        time(n) = rtime

      endif

      if (strpos(mklower(line),"northern") gt -1) then begin

        for j=0,nlons-1 do for i=0,nlats-1 do begin
          readf,1,tmp
          data(0,*,j,i) = tmp
        endfor

      endif

      if (strpos(mklower(line),"southern") gt -1) then begin

        for j=0,nlons-1 do for i=0,nlats-1 do begin
          readf,1,tmp
          data(1,*,j,i) = tmp
        endfor

      endif

      if eof(1) then done = 1

    endwhile

    close,1

    if (n eq 0) then begin

      var = 0
      nvars_to_plot = 0
      while (var ge 0) do begin
        var = fix(ask('Variable Number to plot (-1 to exit)','-1'))
        if (var ge 0) and (var lt nvars) then begin
          if nvars_to_plot eq 0 then $
            varlist = [var]          $
          else varlist = [varlist,var]
          nvars_to_plot = nvars_to_plot + 1
        endif
      endwhile

      print, "You have selected "+tostr(nvars_to_plot)+" variables to plot." 

      if (nvars_to_plot eq 0) then begin
	print, "Can not continue!"
	stop
      endif

      data_to_plot = fltarr(2,nfiles,nvars_to_plot,nlons,nlats)
      theta = fltarr(2,nlons,nlats)
      phi   = fltarr(2,nlons,nlats)

      nt = -1
      for i=0,nvars-1 do if strpos(mklower(vars(i)),'theta') gt -1 then nt = i
      np = -1
      for i=0,2 do if strpos(mklower(vars(i)),'psi') gt -1 then np = i
;I have it named wrong for a number of runs...
      for i=0,2 do if strpos(mklower(vars(i)),'phi') gt -1 then np = i
      if (nt eq -1) or (np eq -1) then begin
        print, "Can't file Theta or Phi variable. Please check file."
        stop
      endif

      theta(*,*,*) = reform(data(*,nt,*,*))
      phi(*,*,*)   = reform(data(*,np,*,*))

    endif

    data_to_plot(*,n,*,*,*) = reform(data(*,varlist,*,*))

  endelse

  print, 'Finished Reading File '+filein

endfor

;
; We want to write an indices file, so figure out maximum & minimum for all
; Selected variables
;

indices = fltarr(2,nfiles,nvars_to_plot,2)
indices_vars = strarr(nvars_to_plot)

for i = 0, nvars_to_plot-1 do indices_vars(i) = vars(varlist(i))

for hem = 0, 1 do for i = 0, nvars_to_plot-1 do for j = 0, nfiles-1 do begin
  indices(hem,j,i,*) = mm(data_to_plot(hem,j,i,*,*))
endfor

;save, indices, indices_vars, file = list(0)+'.ind.save'

;
; Figure out plots per page...
;

space = 0.05

if nfiles eq 1 then begin
  ; if we have just 1 or two variables, just plot northern and southern
  ; hemisphere on the same page. If not, plot them on different pages:
  if nvars_to_plot le 2 then ppp = nvars_to_plot*2 else ppp = nvars_to_plot
endif else begin
  ; if we have less than 15 total images (per hemisphere), plot them all
  ; on the same page. If not, then plot each variable on a seperate page,
  ; with a limit of 9 plots per page.
  if (nfiles le 5) and (nvars_to_plot le 3) then begin
    ppp = nfiles * nvars_to_plot
    space = 0.01
  endif else begin
    if nfiles le 9 then ppp = nfiles else ppp = 9
  endelse
endelse

hems = ["Northern Hemisphere", "Southern Hemisphere"]

for hem = 0, 1 do begin

  print, "Working on ",hems(hem)

  if hem eq 0 then begin
    loc = where(reform(theta(0,0,*)) le mr)
    rang = reform(theta(0,*,loc))
    lons = reform(phi(0,*,loc))*!pi/180.0 + !pi/2
    xpos = rang*cos(lons)
    ypos = rang*sin(lons)
  endif else begin
    loc = where(reform(theta(1,0,*)) gt 180.0-mr)
    rang = 180.0-reform(theta(1,*,loc))
    lons = reform(phi(0,*,loc))*!pi/180.0 + !pi/2
    xpos = rang*cos(lons)
    ypos = rang*sin(lons)
  endelse

  for i = 0, nvars_to_plot-1 do begin

    if (nfiles eq 1) and (nvars_to_plot le 2) then begin
      if hem eq 0 and i eq 0 then begin
        pn = -1
        setdevice, list(0)+'.ps','p',4
	pos_space, ppp, space, sizes
      endif
    endif else begin

      if (nfiles le 5) and (nvars_to_plot le 3) then begin
        if hem eq 1 and i eq 0 then closedevice
        if i eq 0 then begin
          setdevice, list(0)+'_'+strmid(hems(hem),0,5)+'.ps','p',4
          pos_space, ppp, space, sizes, nx = nvars_to_plot
        endif
      endif else begin

        if ((i eq 0) and (nfiles eq 1)) or (nfiles gt 1) then begin
          pn = -1
          if i gt 0 then closedevice
          setdevice, list(0)+'_Var'+tostr(i)+'_'+strmid(hems(hem),0,5)+'.ps',$
                     'p',4
          pos_space, ppp, space, sizes
        endif

      endelse

    endelse

    print, "Working on variable ",vars(varlist(i))

    mini = min(data_to_plot(hem,*,i,*,loc))
    maxi = max(data_to_plot(hem,*,i,*,loc))

    if n_elements(ct_path) eq 0 then begin
      path = !path
      Idl_loc = strpos(path,'Idl')
      if (Idl_loc lt 0) then begin
        ct_path = '.'
      endif else begin
	list = str_sep(path,':')
        nlist_eles = n_elements(list)
	for nl = 0,nlist_eles-1 do begin
          if (strpos(list(nl),'Idl') gt -1) then ct_path = list(nl)
        endfor
      endelse
    endif

    if (mini ge 0.0) then begin
      levels = maxi*findgen(9)/8.0 + mini
      readct,ncolors, ct_path+"/white_red.ct"
    endif else begin
      maxi = max([abs(mini),maxi])
      mini = -maxi
      levels = fltarr(8)
      levels(0:3) = -maxi*(4-findgen(4))/4.0
      levels(4:7) = maxi*(findgen(4)+1)/4.0
      readct,ncolors, ct_path+"/blue_white_red.ct"
    endelse

    c_levels = (maxi-mini)*findgen(30)/29.0 + mini
    c_colors = (ncolors-1)*findgen(30)/29.0 + 1

    for n = 0, nfiles - 1 do begin

      if (nfiles le 5) and (nvars_to_plot le 3) then begin
        if (nfiles eq 1) and (nvars_to_plot le 2) then pn = (pn + 1) mod ppp $
        else pn = n*nvars_to_plot + i
      endif else begin 
        pn = (pn + 1) mod ppp
      endelse

      if pn eq 0 then plotdumb
      get_position, ppp, space, sizes, pn, pos

      contour, data_to_plot(hem,n,i,*,loc), xpos, ypos,        $
        /follow, nlevels=30, /noerase, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], levels = c_levels, c_colors = c_colors,/cell_fill

      contour, data_to_plot(hem,n,i,*,loc), xpos, ypos,        $
        /follow, levels=levels, $
	pos = pos, xstyle = 5, ystyle = 5, xrange = [-mr,mr], 		$
	yrange = [-mr,mr], /noerase


;--------------------------------------------------------------
; Figure out where we are on the page, and whether we need to
; labels or not for the MLT grid
;--------------------------------------------------------------

      no00 = 1
      no06 = 1
      no12 = 1
      no18 = 1

      if (pn+1 gt ppp-sizes.nbx) then no00 = 0
      if (pn mod sizes.nbx eq sizes.nbx-1) then no06 = 0
      if (pn lt sizes.nbx) then no12 = 0
      if (pn mod sizes.nbx eq 0) then no18 = 0

;--------------------------------------------------------------
; Draw the MLT grid
;--------------------------------------------------------------

      plotmlt, mr, no00 = no00, no06 = no06, no12 = no12, no18 = no18

      mini_tmp = min(data_to_plot(hem,n,i,*,loc))
      maxi_tmp = max(data_to_plot(hem,n,i,*,loc))

      if (abs(mini_tmp) gt 1000.0) or (abs(maxi_tmp) gt 1000.0) or        $
         (abs(maxi_tmp) lt 0.1) then begin
        maxs = string(maxi_tmp,format="(e8.2)")
        mins = string(mini_tmp,format="(e9.2)")
      endif else begin
        maxs = string(maxi_tmp,format="(f5.2)")
        mins = string(mini_tmp,format="(f7.2)")
      endelse

      charsize = 1.0
      if ppp gt 9 then charsize = 0.75

      xyouts, pos(0),pos(1), mins, /norm, charsize = charsize
      xyouts, pos(2),pos(1), maxs, /norm, align=1.0, charsize = charsize

      if (nfiles eq 1) then begin
        xyouts, pos(2),pos(3),vars(varlist(i)), 	$
	  align=1.0, /norm
      endif

      if (nfiles le 5) and (nvars_to_plot le 3) then begin
        if i eq 0 then begin
          c_r_to_a, itime, time(n)
          c_a_to_s, itime, stime
          xyouts, pos(0)-0.03, (pos(1)+pos(3))/2.0, strmid(stime,0,15), $
                  alignment = 0.5, /norm, orient = 90
        endif
        if n eq 0 then begin
          xyouts, (pos(0)+pos(2))/2.0, 1.01, vars(varlist(i)), $
                  alignment = 0.5, /norm
          if pn eq 0 then begin
            xyouts, 0.5, 1.05, hems(hem), alignment = 0.5, /norm
            xyouts, -0.01, -0.01, filein, /norm, charsize = 0.5, orient = 90
          endif
        endif
      endif else begin

        if pn eq 0 then begin
          if (nfiles gt 1) then begin
            xyouts, 0.0, 1.01, hems(hem), alignment = 0.0, /norm
            xyouts, 1.0, 1.01, vars(varlist(i)), alignment = 1.0, /norm, $
                  charsize = 1.25
          endif else begin
            xyouts, 0.5, 1.01, hems(hem), alignment = 0.5, /norm
          endelse
          xyouts, -0.01, -0.01, filein, /norm, charsize = 0.5
        endif

      endelse

    endfor

  endfor

endfor

closedevice

end


