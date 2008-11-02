;^CFG COPYRIGHT VAC_UM
;=============================================================================
;
; Read the log data from 1, 2, or 3 files into the wlog, wlog1, wlog2 arrays
; Store the names of the variables into wlognames, wlognames1 and wlognames2
; Calculate and store log time into the logtime, logtime1 and logtime2
; arrays. The time units are set by timeunit.
;

askstr,'logfilename(s) ',logfilename,doask

nlogfile=0
if stregex(logfilename, '[?*[]', /boolean) then begin
    spawn,'ls '+logfilename, logfilenames
    nlogfile = n_elements(logfilenames)
endif else $
  str2arr,logfilename,logfilenames,nlogfile

if nlogfile gt 3 then begin
   print,'Error in getlog: cannot handle more than 3 files.'
   retall
endif

get_log,   logfilenames(0), wlog,  wlognames,  logtime, timeunit, /verbose
if nlogfile ge 2 then $
  get_log, logfilenames(1), wlog1, wlognames1, logtime1, timeunit, verbose='1'
if nlogfile eq 3 then $
  get_log, logfilenames(2), wlog2, wlognames2, logtime2, timeunit, verbose='2'

end
