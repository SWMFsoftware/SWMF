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

if nlogfile gt 10 then begin
   print,'Error in getlog: cannot handle more than 10 files.'
   retall
endif

get_log,   logfilenames(0), wlog,  wlognames,  logtime, timeunit, /verbose
if nlogfile ge 2 then $
  get_log, logfilenames(1), wlog1, wlognames1, logtime1, timeunit, verbose='1'
if nlogfile ge 3 then $
  get_log, logfilenames(2), wlog2, wlognames2, logtime2, timeunit, verbose='2'
if nlogfile ge 4 then $
  get_log, logfilenames(3), wlog3, wlognames3, logtime3, timeunit, verbose='3'
if nlogfile ge 5 then $
  get_log, logfilenames(4), wlog4, wlognames4, logtime4, timeunit, verbose='4'
if nlogfile ge 6 then $
  get_log, logfilenames(5), wlog5, wlognames5, logtime5, timeunit, verbose='5'
if nlogfile ge 7 then $
  get_log, logfilenames(6), wlog6, wlognames6, logtime6, timeunit, verbose='6'
if nlogfile ge 8 then $
  get_log, logfilenames(7), wlog7, wlognames7, logtime7, timeunit, verbose='7'
if nlogfile ge 9 then $
  get_log, logfilenames(8), wlog8, wlognames8, logtime8, timeunit, verbose='8'
if nlogfile gt 9 then $
  get_log, logfilenames(9), wlog9, wlognames9, logtime9, timeunit, verbose='9'

end
