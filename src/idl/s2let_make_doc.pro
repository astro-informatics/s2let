pro s2let_make_doc

loc = GETENV('S2LET')

if loc eq '' then begin
   print,'You must define environment variable S2LET (base directory of s2let)'
endif else begin
   codeloc = loc + '/src/main/idl/'
   docfile = loc + '/doc/idl/index.html'
   MK_HTML_HELP, codeloc, docfile, /verbose, title='S2LET library : documentation for IDL routines and interfaces'
endelse

end
