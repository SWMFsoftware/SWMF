;------------------------------------------------
; rough IDL routine to take EUV / SXR los files 
; and convert them into PNG images
;------------------------------------------------
;Inputs
; FileInEuv: input EUV file
; FileInSxr: input SXR file
; PngOut171: ouput PNG of EUV 171
; PngOut195: ouput PNG of EUV 195
; PngOut284: ouput PNG of EUV 284
; PngOutSXR: ouput PNG of SXR
; nPix: number of pixels of input image 
;
;Comments
; 
; This is designed to read the ascii tecplot output format
; from the SWMF los euv tec command
; (NOT THE IDL formatted saves!)
;
; Right now to get the color tables you need
; the SSW framework on the system  with eit
; installed in order to use 'eit_colors' routine
; (Ideally eit_colors should become part of the framework)

pro euv, FileInEuv, FileInSxr, PngOut171, PngOut195, PngOut284, PngOutSXR, nPix

; load the EUV file
hdr=strarr(3)
npts=(1L*npix)^2
datavector=fltarr(5,npts)
openr,1,FileInEuv
readf,1,hdr
readf,1,datavector
close,1
i1=fltarr(npix,npix)
i2=fltarr(npix,npix)
i3=fltarr(npix,npix)
for i=0,npix-1 do begin
	for j=0,npix-1 do begin
		i1[i,j]=datavector[2,npix*1L*i+j]
	endfor
endfor
for i=0,npix-1 do begin
	for j=0,npix-1 do begin
		i2[i,j]=datavector[3,npix*1L*i+j]
	endfor
endfor
for i=0,npix-1 do begin
	for j=0,npix-1 do begin
		i3[i,j]=datavector[4,npix*1L*i+j]
	endfor
endfor

; load the SXR file
hdr=strarr(3)
npts=(1L*npix)^2
datavector=fltarr(3,npts)
openr,1,FileInSxr
readf,1,hdr
readf,1,datavector
close,1
i4=fltarr(npix,npix)
for i=0,npix-1 do begin
	for j=0,npix-1 do begin
		i4[i,j]=datavector[2,npix*1L*i+j]
	endfor
endfor

; set some resizing parameters since going to cut into square
ntot=(size(i1))[1]
smult=2
pix=floor(ntot*(1-1.4/2.0)/2)
imx=smult*(ntot-2*pix)
imy=smult*(ntot-2*pix)

; this is final image size
size=imx/2

; congrid the images to right size
it1=congrid(i1[pix:ntot-pix-1,pix:ntot-pix-1],size,size)
it2=congrid(i2[pix:ntot-pix-1,pix:ntot-pix-1],size,size)
it3=congrid(i3[pix:ntot-pix-1,pix:ntot-pix-1],size,size)
it4=congrid(i4[pix:ntot-pix-1,pix:ntot-pix-1],size,size)

; save 171
eit_colors,42,/silent
pmin=0.02 & pmax = 3d3
imscl=bytscl(alog10(it1>pmin<pmax),min=alog10(pmin),max=alog10(pmax))
tvlct,r,g,b,/get
imrgb=bytarr(3,size,size)
imrgb[0,*,*] = r[imscl]
imrgb[1,*,*] = g[imscl]
imrgb[2,*,*] = b[imscl]
write_png,PngOut171,imrgb

; save 195
eit_colors,43,/silent
pmin=0.02 & pmax = 3d3
imscl=bytscl(alog10(it2>pmin<pmax),min=alog10(pmin),max=alog10(pmax))
tvlct,r,g,b,/get
imrgb=bytarr(3,size,size)
imrgb[0,*,*] = r[imscl]
imrgb[1,*,*] = g[imscl]
imrgb[2,*,*] = b[imscl]
write_png,PngOut195,imrgb

; save 284
eit_colors,44,/silent
pmin=0.02 & pmax = 1d3
imscl=bytscl(alog10(it3>pmin<pmax),min=alog10(pmin),max=alog10(pmax))
tvlct,r,g,b,/get
imrgb=bytarr(3,size,size)
imrgb[0,*,*] = r[imscl]
imrgb[1,*,*] = g[imscl]
imrgb[2,*,*] = b[imscl]
write_png,PngOut284,imrgb

; save sxr
loadct,3,/silent
pmin=0.01 & pmax = 1d3
imscl=bytscl(alog10(it4>pmin<pmax),min=alog10(pmin),max=alog10(pmax))
tvlct,r,g,b,/get
imrgb=bytarr(3,size,size)
imrgb[0,*,*] = r[imscl]
imrgb[1,*,*] = g[imscl]
imrgb[2,*,*] = b[imscl]
write_png,PngOutSXR,imrgb

end
