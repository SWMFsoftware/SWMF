pro close_device, pdf=pdf, delete=delete

  closedevice

  if not keyword_set(pdf) then return

  ; Convert PS/EPS file to PDF and remove original if required
  common SETDEVICE, NameFile
  PdfFile = NameFile
  i = strpos(PdfFile,'.',/reverse_search)
  PdfFile = strmid(PdfFile, 0, i) + '.pdf'
  Command = 'convert ' + NameFile + ' ' + PdfFile
  if keyword_set(delete) then Command = Command + '; /bin/rm ' + NameFile
  spawn, Command

end
