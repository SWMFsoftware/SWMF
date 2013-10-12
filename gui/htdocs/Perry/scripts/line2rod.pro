;  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
;  For more information, see http://csem.engin.umich.edu/tools/swmf
; This macro will take a line trace, convert it into a rod with
; variable radius, and write it out to tmp.dat.

nPtsInCircle = 21
angle=double(2.0*3.141592654/(nPtsInCircle-1))

; Read file. nLines is number of data lines to read
file = "satellite5min.dat"
file = ask("file name to process",file)
nLines=file_lines(file)

; Read options.
rVariesString = "no"
rVariesString = ask("variable rod radius (yes/no)",rVariesString)
rVaries = strpos(rVariesString,'yes')

radius=double(0.1)
answer = ''
read, 'Enter rod radius ['+strcompress(string(radius),/remove_all)+'] : ', answer
if strlen(answer) gt 0 then radius=double(answer)

if rVaries ne -1 then begin
  print, 'Using variable radius ',radius
endif else begin
  print, 'Using fixed radius ',radius
endelse

; Find out how many zones there are
sarr = strarr(nLines)
openr,unit,file,/get_lun
readf,unit,sarr
free_lun,unit
nZones=0L
for i=0L,nLines-1 do begin
  result=strcmp('ZONE',sarr(i),4,/fold_case)
  if result eq 1 then begin
    nZones=nZones+1
  endif
endfor
print, 'Zones found:',nZones

; Get number of points to read
nStart=lonarr(nZones)
nEnd=lonarr(nZones)
k=0L
for i=1L,nLines do begin
  result=strcmp('ZONE',sarr(i-1),4,/fold_case)
  if result eq 1 then begin
    k=k+1
    nStart(k-1)=i+4
    if k ne 1 then nEnd(k-2)=i-1
  endif
endfor
nEnd(k-1)=nLines

; open files to write incremental data
openw,unitX,'tmpX',/get_lun
openw,unitP,'tmpP',/get_lun

shift=0L
; Loop over zones
for z=0L,nZones-1 do begin

  ; Read xyzR
  n=nEnd(z)-nStart(z)+1
  xyzR = dblarr(3,n)
  dummy = ""
  openr,1,file
  for i=1L,nStart(z)-1 do readf,1,dummy
  readf,1,xyzR
  close,1

  ; Find points to actually use.  Skip if step less than radius.
  lR = dblarr(n)
  iR = lonarr(n)
  nUsed = 0L
  for i=0L,n-1 do begin
    if i eq 0 then begin
      lR(i)=0.
      nUsed=nUsed+1
      iR(nUsed-1)=i
    endif else begin
      lR(i) = lR(i-1)+sqrt((xyzR(0,i)-xyzR(0,i-1))^2+(xyzR(1,i)-xyzR(1,i-1))^2+(xyzR(2,i)-xyzR(2,i-1))^2)
      if (lR(i)-lR(iR(nUsed-1))) gt radius then begin
        nUsed=nUsed+1
        iR(nUsed-1)=i
      endif
    endelse
  endfor
  iR(nUsed-1)=n-1

  print, 'zone:',z+1,'  n:',n,'  nUsed:',nUsed,'  length:',lR(n-1)

  ; Create rod of points around read in values.  Result to xyz array.
  Nv = dblarr(3)
  Av = dblarr(3)
  Bv = dblarr(3)
  xyz = dblarr(3,nPtsInCircle*nUsed)
  for i=0L,nUsed-1 do begin

    ; Unit normal vector from point before and after
    if i eq 0 then begin
      Nv(*)=xyzR(*,iR(i+1))-xyzR(*,iR(i))
    endif else begin
      if i eq (nUsed-1) then begin
        Nv(*)=xyzR(*,iR(i))-xyzR(*,iR(i-1))
      endif else begin
        Nv(*)=xyzR(*,iR(i+1))-xyzR(*,iR(i-1))
      endelse
    endelse
    norm = sqrt(mean(Nv^2)*n_elements(Nv))
    Nv(*)=Nv(*)/norm

    ; Unit vector Av (Av = (1,0,0) x Nv)
    Av(0)=0.
    Av(1)=-Nv(2)
    Av(2)=Nv(1)
    norm = sqrt(mean(Av^2)*n_elements(Av))
    Av(*)=Av(*)/norm

    ; Unit vector Bv (Bv = Av x Nv)
    Bv(0)=Av(1)*Nv(2)-Nv(1)*Av(2)
    Bv(1)=Av(2)*Nv(0)
    Bv(2)=-Nv(0)*Av(1)
    norm = sqrt(mean(Bv^2)*n_elements(Bv))
    Bv(*)=Bv(*)/norm

    r=radius
    if rVaries ne -1 then begin
      norm = sqrt(xyzR(0,iR(i))^2+xyzR(1,iR(i))^2+xyzR(2,iR(i))^2)
      if norm gt 25. then norm=25.
      r=radius*2.*norm
    endif
    for j=0,nPtsInCircle-1 do begin
      ; Untransformed values
      xyz(0,nPtsInCircle*i+j)=xyzR(0,iR(i))+r*(cos(j*angle)*Av(0)+sin(j*angle)*Bv(0))
      xyz(1,nPtsInCircle*i+j)=xyzR(1,iR(i))+r*(cos(j*angle)*Av(1)+sin(j*angle)*Bv(1))
      xyz(2,nPtsInCircle*i+j)=xyzR(2,iR(i))+r*(cos(j*angle)*Av(2)+sin(j*angle)*Bv(2))
    endfor

  endfor

  ; 
  e = (nUsed-1)*(nPtsInCircle-1)
  p = ulonarr(5,e)
  p(0,*) = 4
  j=0L-1
  for i=0L,((nUsed-1)*nPtsInCircle)-1 do begin
    if (i mod nPtsInCircle eq nPtsInCircle-1) then begin
    endif else begin
      j=j+1
      p(1,j) = shift+i
      p(2,j) = shift+i  +nPtsInCircle
      p(3,j) = shift+i+1+nPtsInCircle
      p(4,j) = shift+i+1
    endelse
  endfor

  ; Write temporary file
  printf,unitX,xyz
  printf,unitP,p

  ; Save shift
  shift=shift+nPtsInCircle*nUsed

endfor

; Close up temporary files
free_lun,unitX
free_lun,unitP

; Reload completed files
nLinesX=file_lines('tmpX')
xyz=dblarr(3,nLinesX)
openr,unitX,'tmpX',/get_lun
readf,unitX,xyz
free_lun,unitX

nLinesP=file_lines('tmpP')
p=ulonarr(5,nLinesP)
openr,unitP,'tmpP',/get_lun
readf,unitP,p
free_lun,unitP

; Write tecplot file
openw,unitT,'tmp.dat',/get_lun
printf,unitT,'TITLE = "IDL script output"'
printf,unitT,'VARIABLES = "X [R]"'
printf,unitT,' "Y [R]"'
printf,unitT,' "Z [R]"'
printf,unitT,'ZONE T="surface"'
printf,unitT,' N=',nLinesX,', E=',nLinesP,', ZONETYPE=FEQuadrilateral'
printf,unitT,' DATAPACKING=POINT'
printf,unitT,' DT=(SINGLE SINGLE SINGLE)'
printf,unitT,xyz
for i=0L,nLinesP-1 do begin
  printf,unitT,p(1,i)+1,p(2,i)+1,p(3,i)+1,p(4,i)+1
endfor
free_lun,unitT

end

