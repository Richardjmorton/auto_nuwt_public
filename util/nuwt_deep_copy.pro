;+
;NAME: NUWT_DEEP_COPY
;
;PURPOSE:
;   Makes a 'deep copy' of the input variable
;
;INPUTS:
;   [to come later]
;
;OPTIONAL INPUTS:
;   [to come later]
;
;OUTPUTS:
;   [TO COME LATER]
;
;HISTORY: Name---------Date---------Description
;         M Weberg  31 AUG, 2017  Initial coding
;
;TO-DO / LIMITATIONS:
;   - WRITE DOC STRING
;-
FUNCTION NUWT_DEEP_COPY, input_var

COMPILE_OPT IDL2

temp_filename = 'nuwt_deep_copy_temp_file.sav'
var = input_var
SAVE, var, filename=temp_filename
var = 0 ;This is apparently a key step to cut the chain of references
RESTORE, temp_filename
FILE_DELETE, temp_filename
RETURN, var
END
