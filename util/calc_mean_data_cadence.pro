;+
;NAME: CALC_MEAN_DATA_CADENCE
;
;PURPOSE:
;   Read an array of time stamps and estimate the mean cadence
;
;INPUTS:
;   fits_index - index array from the fits file headers. MUST include a structures
;                tag named "date_obs" (standard at least for SDO / AIA data).
;
;OPTIONAL INPUTS
;   /trim - if set, will attempt to trim off any fractional seconds from the end of the
;           timestamp
;
;OUTPUTS:
;   mean_cadence - mean data cadence found in units of [s]
;   time_diff - array of time differnces (in units of [s]) between adjacent
;               timestamps
;
;HISTORY: Name---------Date---------Description
;         M Weberg  12 JUNE, 2017  Initial coding
;         M Weberg  12 SEPT, 2017  Added "time_diff" optional output
;-

FUNCTION CALC_MEAN_DATA_CADENCE, fits_index, trim=trim, time_diff=time_diff
;WARNING! currently only works if all of the data is from the same day
COMPILE_OPT IDL2

num_timesteps = N_ELEMENTS(fits_index.date_obs)

IF KEYWORD_SET(trim) THEN BEGIN
    trim_locs = STRSPLIT(fits_index[0].date_obs, '.') ;find where the "." is
    trim_fits_date = STRMID(fits_index.date_obs, 0, LONG(trim_locs[-1]-1)) ;cuts off the fractional seconds
ENDIF ELSE BEGIN
    trim_fits_date = fits_index.date_obs ;keeps the full timestamp
ENDELSE

split_date = STRSPLIT(trim_fits_date, 'T:', /EXTRACT) ; YYYY-MM-DD, Hr, Min, Sec

datetime_arr = split_date.ToArray() ;array of form [timestep_col, time_element_row]
time_in_sec = float(datetime_arr[*,1])*3600.0+float(datetime_arr[*,2])*60.0+float(datetime_arr[*,3])

sec_diff = time_in_sec - SHIFT(time_in_sec, 1)

mean_cadence = MEAN(sec_diff[1:-1]) ;ignore 0th element since SHIFT wraps around
time_diff = sec_diff

RETURN, mean_cadence
END
