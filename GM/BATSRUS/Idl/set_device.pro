;^CFG COPYRIGHT UM
pro set_device, psfile, port=port, land=land, eps=eps, $
                psfont=psfont, percent=percent

  ; Parameter defaults and conversions

  if not keyword_set(psfile) then psfile = 'idl.ps'
  orientation = 'normal'
  if keyword_set(port) then orientation='port'
  if keyword_set(land) then orientation='land'
  if not keyword_set(percent) then percent = 1.0		$
  else if percent gt 1.0 then percent = float(percent)/100.0

  if n_elements(psfont) eq 0  then psfont = 28

  common SETDEVICE, NameFile
  NameFile = psfile

  ; Set sizes and offsets
  case (orientation) of
     'normal': begin
        xs   = 10.0*percent
        ys   = 7.0 *percent
        xoff = 11.0-(11.0-xs)/1.5
        yoff = (8.5-ys)/2.0
        land=0
     end
     'land': begin
        xs   = 10.0*percent
        ys   = 7.0 *percent
        xoff = (8.5-ys)/2.0
        yoff = 11.0-(11.0-xs)/1.5
        land=1
     end
     'port': begin
        xs = 7.5*percent
        ys = 9.5*percent
        xoff = (8.5-xs)/2.0
        yoff = (11.0-ys)/2.0
        land=0
     end
  endcase

  set_plot, 'PS', /copy, /interpolate

  !p.font = 0

  case (psfont) of
       -1  : device, file = psfile, encapsulated=eps, /color, bits=8,         $
                /inches, landscape=land, xsize = xs, ysize = ys, $
               xoff = xoff, yoff = yoff

	0  : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Courier 
	1  : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Courier, /Bold 
    	2  : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Courier, /Oblique 
	3  : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Courier, /Bold, /Oblique
       	4  : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Helvetica
      	5  : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Helvetica, /Bold
    	6  : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Helvetica, /Oblique
       	8  : device, file = psfile, encapsulated=eps, /color, bits=8,         $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Helvetica, /Bold, /Oblique 
    	12 : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Avantgarde, /Book 
     	13 : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Avantgarde, /Book, /Oblique
	14 : device, file = psfile, encapsulated=eps, /color, bits=8,	$
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Avantgarde, /Demi 
      	15 : device, file = psfile, encapsulated=eps, /color, bits=8,	      $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Avantgarde, /Demi, /Oblique
       	20 : device, file = psfile, encapsulated=eps, /color, bits=8, $
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Schoolbook
   	21 : device, file = psfile, encapsulated=eps, /color, bits=8,$
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Schoolbook, /Bold
      	22 : device, file = psfile, encapsulated=eps, /color, bits=8,$
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Schoolbook, /Italic
       	23 : device, file = psfile, encapsulated=eps, /color, bits=8,	$
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Schoolbook, /Bold, /Italic 
	28 : device, file = psfile, encapsulated=eps, /color, bits=8,	$
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Times
	29 : device, file = psfile, encapsulated=eps, /color, bits=8,	$
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Times, /Bold
	30 : device, file = psfile, encapsulated=eps, /color, bits=8,	$
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Times, /Italic
	31 : device, file = psfile, encapsulated=eps, /color, bits=8,	$
		/inches, landscape=land, xsize = xs, ysize = ys, $
		xoff = xoff, yoff = yoff,  $
		/Times, /Bold, /Italic

    endcase

  return

end

