;+
;NAME: NUWT_CONVERT_MONTH_NAME
;
;PURPOSE:
;   Converts month names to month numbers. Works fine for truncated names.
;
;INPUTS:
;   [TO BE UPDATED]
;
;OPTIONAL OUTPUTS:
;   [TO BE UPDATED]

;
;HISTORY: Name---------Date---------Description
;         M Weberg  23 MAR, 2017   Initial coding (as part of "FETCH_AIA_DATA.pro")
;         M Weberg  15 SEPT, 2017  Split function out into its own file
;
;TO DO:
;   - Add options for converting entire arrays of input values.
;-

FUNCTION NUWT_CONVERT_MONTH_NAME, input_str

COMPILE_OPT IDL2

;Converts month names to month numbers. Works fine for truncated names
month_names = ['January', 'February', 'March', 'April', $
               'May', 'June', 'July', 'August', $
               'September', 'October', 'November', 'December']

output_str = input_str ;defaults to returning the input string
FOR m=0, 11 DO BEGIN
    test_month = month_names[m]
    IF test_month.StartsWith(input_str, /FOLD_CASE) THEN output_str = str(m+1)
ENDFOR

RETURN, output_str
END
