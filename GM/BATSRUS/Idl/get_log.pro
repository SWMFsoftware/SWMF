;^CFG COPYRIGHT VAC_UM

;=============================================================================
pro get_log,unit=unit,file=file,wlog,wlognames,verbose=verbose
;
; Read the log data from a unit or a file into wlog and the variable 
; names into wlognames. If verbose is set show verbose information.
; If versbose is a string, attach it to 'wlog' in the verbose info.
; 

if not keyword_set(unit) and not keyword_set(file) then begin
    print,'get_log error: neither unit nor file is set!'
    retall
end

if not keyword_set(unit) then begin
    ; open file
    unit = 1
    data = fstat(unit)
    if data.open ne 0 then close,1
    openr,unit,file
endif

if not keyword_set(verbose) then verbose = 0
; If verbose is a string set the index string to it
if size(verbose,/type) eq 7 then index=verbose else index=''

headline=''
readf,unit,headline
wlognames=''
readf,unit,wlognames

wlognames=str_sep(strtrim(strcompress(wlognames),2),' ')
nwlog=n_elements(wlognames)

if verbose then begin
  if keyword_set(file) then print,'logfile',index,'  =',file
  print,'headline',index,' =',strtrim(headline,2)
  for i=0,nwlog-1 do $
    print,FORMAT='("  wlog",A,"(*,",I2,")= ",A)',index,i,wlognames(i)
endif

; Use buffers for efficient reading
buf   = long(10000)
dbuf  = long(10000)
wlog  = dblarr(nwlog,buf)
wlog_ = dblarr(nwlog)
nt    = long(0)
while not eof(unit) do begin
    on_ioerror,close_file
    readf,unit,wlog_
    wlog(*,nt)=wlog_
    nt=nt+1
    if nt ge buf then begin
        buf=buf+dbuf
        wlog=[[wlog],[dblarr(nwlog,buf)]]
    endif
endwhile
close_file: if not keyword_set(unit) then close,unit

if verbose then print,'Number of recorded timesteps: nt=',nt

wlog = transpose(wlog(*,0:nt-1))

end
