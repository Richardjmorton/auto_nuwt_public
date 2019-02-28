;+
;NAME: GET_CURRENT_TIMESTAMP
;
;PURPOSE:
;   Create a string with the current time. This is useful for creating unique
;   filenames. By default, the local time is given in the format YYYYMMDD_hhmm.
;
;INPUTS:
;   No required inputs
;
;OPTIONAL INPUTS
;   /utc - if set, will output the time in UTC instead of the local time zone
;   format - string containing the desired output format. Can accept a wide
;            range of formats as long as they contain one or more of the
;            following format codes:
;               YYYY - four-digit year
;               MM - two-digit month number (01 = Jan, ..., 12 = Dec)
;               DD - two-digit day number (01 - 31)
;               hh - two-digit hour in 24-hour format (00 - 23)
;               mm - two-digit minutes (00 - 59)
;               ss - two-digit second (00 - 59)
;            Format codes ARE case sensitive and all numbers are zero-padded to
;            have constant widths. Any extra characters in the input format
;            string will be copied to the same relative location in the the output
;            string. The default format is 'YYYYMMDD_hhmm'
;
;OUTPUTS:
;   timestamp_str - string with the current time in the requested format
;
;HISTORY: Name---------Date---------Description
;         M Weberg  25 JAN, 2018  Initial coding
;
;-

FUNCTION GET_CURRENT_TIMESTAMP, utc=utc, format=format
COMPILE_OPT IDL2

IF KEYWORD_SET(UTC) THEN utc = 1 ELSE UTC = 0
IF NOT KEYWORD_SET(format) THEN format = 'YYYYMMDD_hhmm'

timestamp_str = format[0] ;if given an array, only use the first format

CALDAT, SYSTIME(/julian, utc=utc), month_int, day_int, year_int, hr_int, min_int, sec_flt

regex_arr = ['YYYY', 'MM', 'DD', 'hh', 'mm', 'ss']
time_arr = [year_int, month_int, day_int, hr_int, min_int, sec_flt]
fmt_arr = ['(I04)', '(I02)', '(I02)', '(I02)', '(I02)', '(I02)']
FOR i=0, 5 DO BEGIN
    timestamp_str = timestamp_str.Replace(regex_arr[i], STRING(time_arr[i], format=fmt_arr[i]))
ENDFOR

RETURN, timestamp_str
END
