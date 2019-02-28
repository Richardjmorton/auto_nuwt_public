;+
;NAME: CONVERT_NUWT_UNITS
;
;PURPOSE:
;   Calculates the conversion factor needed to convert the units of one of the
;   primary NUWT output paramters. Note, this code does not actully apply the
;   unit conversion in the output parameter array.
;
;INPUTS:
;   selected_parameter - Parameter or units to convert FROM. May select from either
;                        one of the following NUWT parameters:
;                        ['amp', 'freq', 'period', td_dist', 'td_time', 'enbw']
;                        or one of the valid output units
;   requested_units - desired output units. May select from one of the following:
;                     ['pixels', 'arcsec', 'm', 'km', 'Mm', 'timesteps', 's',
;                      'min', 'hr', 'ks', 'timesteps^{-1}', 'Hz', mHz']
;                     The user is responsible for ensuring the input and output
;                     unit types are compatible
;
;OPTIONAL INPUTS:
;   res - spatial resolution of the data (in [arcsec]). If not given by the user,
;         the program will load the 'res' value stored by NUWT for the current
;         set of results. If NUWT does not have a valid res value, then all
;         distance values will revert to units of [pixels]
;   cad - temporal cadence of the data (in [s]). If not given by the user, the
;         program will load the 'cad' value stored by NUWT for the current set of
;         results. If NUWT does not have a valid cad value, then all time values
;         will revert to units of [timesteps]
;   km_per_arcsec - ratio of [km]/[arcsec]. Defaults to 725.27 which is the
;                   mean scale of the solar surface as seen from 1 AU.
;   /no_autoload - If set, will NOT autoload unit, res, or cad values from NUWT.
;                  This is good for wehn you just want to perform simple unit
;                  conversions without  using NUWT information.
;   /quiet - If set, will not NOT print a warning when the output unit is reset
;            due to incomplete "res" or "cad" information.
;
;OUTPUTS:
;   conversion_factor - multiplicative factor needed to convert from the units
;                       stored by NUWT to the requested units
;   output_units - units actually outputted by the given conversion factor. In
;                  some cases, the requested unit conversion may not be possible
;                  using the given information.
;
;HISTORY: Name---------Date---------Description
;         M Weberg  16  FEB, 2018  Initial coding
;         M Weberg  04  APR, 2018  Added the /quiet keyword
;
;TO-DO / LIMITATIONS:
;   - add a check to see if the input and output units are compatable
;-

FUNCTION CONVERT_NUWT_UNITS, selected_parameter, requested_units, $
                             res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                             output_units=output_units, $
                             no_autoload=no_autoload, $
                             quiet=quiet, $
                             use_temp_common_blocks=use_temp_common_blocks

COMPILE_OPT IDL2

IF KEYWORD_SET(no_autoload) THEN BEGIN
    ;Ignore NUWT information and just do simple unit conversions

    user_res_override = 0
    IF KEYWORD_SET(res) THEN BEGIN
        ;user input res
        IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27
        res_val = res
        res_unknown = 0
    ENDIF ELSE BEGIN
        ;no res value given
        IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27
        res_val = 1.0
        res_unknown = 1
    ENDELSE

    user_cad_override = 0
    IF KEYWORD_SET(cad) THEN BEGIN
        ;user input cad
        cad_val = cad
        nuwt_freq_units = 'timesteps^{-1}'
        user_cad_override = 1
        cad_unknown = 0
    ENDIF ELSE BEGIN
        ;no cad value given
        cad_val = 1.0
        cad_unknown = 1
    ENDELSE
ENDIF ELSE BEGIN
    ;@load_nuwt#####################################################################
    ;LOAD IN NUWT DATA CONTAINING UNIT INFORMATION
    ;###############################################################################
    COMMON all_nuwt_dat, nuwt_located, nuwt_threads, nuwt_fft_spec, nuwt_fft_peaks

    IF NOT KEYWORD_SET(use_temp_common_blocks) THEN BEGIN
        ;[Default] For running after the main NUWT program as been run
        slit_located = nuwt_located[0]
        slit_fft_peaks = nuwt_fft_peaks[0]

        IF KEYWORD_SET(use_refit) THEN BEGIN
            COMMON refit_dat, refit_spec, refit_stats
            ;Quick and dirty swap in structure. use with care!!!
            slit_fft_peaks = refit_stats[0]
        ENDIF
    ENDIF ELSE BEGIN
        ;Used when running as part of the main NUWT data processor (use with care!)
        COMMON located_dat, located
        COMMON fft_results_dat, fft_spec, fft_peaks

        slit_located = located
        slit_fft_peaks = fft_peaks
    ENDELSE

    ;@read_units####################################################################
    ;READ IN UNIT INFORMATION FROM NUWT
    ;###############################################################################
    temp_fft_peaks = slit_fft_peaks[0]

    nuwt_amp_units = temp_fft_peaks.amp_units
    amp_to_pxls = temp_fft_peaks.amp_to_pxls

    nuwt_freq_units = temp_fft_peaks.freq_units
    freq_to_timesteps = temp_fft_peaks.freq_to_timesteps

    IF KEYWORD_SET(res) THEN BEGIN
        ;user override
        IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27
        res_val = res
        nuwt_amp_units = 'pixels'
        user_res_override = 1
        res_unknown = 0
    ENDIF ELSE BEGIN
        ;in NUWT we trust
        IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = slit_located.km_per_arcsec
        res_val = slit_located.res
        user_res_override = 0
        IF slit_located.units_dx EQ 'pixels' THEN res_unknown = 1 ELSE res_unknown = 0
    ENDELSE

    IF KEYWORD_SET(cad) THEN BEGIN
        ;user override
        cad_val = cad
        nuwt_freq_units = 'timesteps^{-1}'
        user_cad_override = 1
        cad_unknown = 0
    ENDIF ELSE BEGIN
        ;in NUWT we trust
        cad_val = slit_located.cad
        nuwt_dt_units = slit_located.units_dt
        user_cad_override = 0
        IF slit_located.units_dt EQ 'timesteps' THEN cad_unknown = 1 ELSE cad_unknown = 0
    ENDELSE

    freq_to_period_units = ORDEREDHASH('timesteps^{-1}', 'timesteps', $
                                       'Hz', 's', $
                                       's^{-1}', 's', $
                                       'mHz', 'ks')
    nuwt_period_units = freq_to_period_units[nuwt_freq_units] ;inverse of the freq units
ENDELSE

;@tables########################################################################
;CREATE UNIT CONVERSION TABLES
;###############################################################################
;Note: if you expand these conversion tables, don't forget to add the new units
;      to the lists in the next section!

index_dist_units = ORDEREDHASH('pixels', 0, 'arcsec', 1, 'm', 2, 'km', 3, 'Mm', 4)
a2km = km_per_arcsec
;FROM [pixels],                  [arcsec],           [m],                 [km],               [Mm]
convert_dist_table = [$
    ; [1.0,             amp_to_pxls,  amp_to_pxls,          amp_to_pxls,        amp_to_pxls], $       ;TO [pixels]
    [1.0,               1.0/res_val,  1e-3/(res_val*a2km),  1.0/(res_val*a2km),  1e3/(res_val*a2km)], $ ;TO [pixels]
    [res_val,           1.0,          1e-3/a2km,            1.0/a2km,           1e3/a2km], $           ;TO [arcsec]
    [res_val*a2km*1e3,  a2km*1e3,     1.0,                  1e3,                1e6], $                ;TO [m]
    [res_val*a2km,      a2km,         1e-3,                 1.0,                1e3], $                ;TO [km]
    [res_val*a2km/1e3,  a2km/1e3,     1e-6,                 1e-3,               1.0]]                  ;TO [Mm]

index_time_units = ORDEREDHASH('timesteps', 0, 's', 1, 'min', 2, 'hr', 3, 'ks', 4)
;FROM [timesteps],     [s],          [min],      [s],          [ks]
convert_time_table = [$
    [1.0,             1.0/cad_val,  60.0/cad_val,  3600.0/cad_val,  1e3/cad_val], $    ;TO [timesteps]
    [cad_val,         1.0,          60.0,          3600.0,          1e3], $        ;TO [s]
    [cad_val/60.0,    1.0/60.0,     1.0,           60.0,            1e3/60.0], $   ;TO [min]
    [cad_val/3600.0,  1.0/3600.0,   1.0/60.0,      1.0,             1e3/3600.0], $ ;TO [hr]
    [cad_val/1e3,     1e-3,         60.0/1e3,      3600.0/1e3,      1.0]]          ;TO [ks]

index_freq_units = ORDEREDHASH('timesteps^{-1}', 0, '1/timesteps', 1, 'Hz', 2, 's^{-1}', 3, '1/s', 4, 'mHz', 5)
;FROM [timesteps^{-1}], [1/timesteps],  [Hz],            [s^{-1}],            [1/s],             [mHz]
convert_freq_table = [$
    ; [1.0,           1.0,        freq_to_timesteps, freq_to_timesteps, freq_to_timesteps, freq_to_timesteps], $ ;TO [timesteps^{-1}]
    ; [1.0,           1.0,        freq_to_timesteps, freq_to_timesteps, freq_to_timesteps, freq_to_timesteps], $ ;TO [1/timesteps]
    [1.0,           1.0,              cad_val,           cad_val,           cad_val,           cad_val/1e3], $ ;TO [timesteps^{-1}]
    [1.0,           1.0,              cad_val,           cad_val,           cad_val,           cad_val/1e3], $ ;TO [1/timesteps]
    [1.0/cad_val,   1.0/cad_val,        1.0,               1.0,               1.0,               1e-3], $       ;TO [Hz]
    [1.0/cad_val,   1.0/cad_val,        1.0,               1.0,               1.0,               1e-3], $       ;TO [s^{-1}]
    [1.0/cad_val,   1.0/cad_val,        1.0,               1.0,               1.0,               1e-3], $       ;TO [1/s]
    [1e3/cad_val,   1e3/cad_val,        1e3,               1e3,               1e3,               1.0]]          ;TO [mHz]

;@check#########################################################################
;VALIDATE REQUESTED UNITS
;###############################################################################
;Note: if you add to these lists, don't forget to expand the conversion tables
;      above as needed!
valid_parameters = ['amplitude', 'amp', 'frequency', 'freq', 'period', $
                    'td_distance', 'td_dist', 'td_timestep', 'td_time', $
                    'enbw']

valid_units = ['pixels', 'arcsec', 'm', 'km', 'Mm', $
               'timesteps', 's', 'min', 'hr', 'ks', $
               'timesteps^{-1}', '1/timesteps', 'Hz', 's^{-1}', '1/s', 'mHz']

IF KEYWORD_SET(no_autoload) THEN BEGIN
    IF TOTAL(STRMATCH(valid_units, selected_parameter)) EQ 0 THEN BEGIN
        PRINT, 'ERROR: Invaild input units!'
        PRINT, 'Please select from one of the following units:'
        FOR i=0, (N_ELEMENTS(valid_units) -1) DO PRINT, '   '+valid_units[i]
        MESSAGE, selected_parameter+' cannot be used in the /no_autoload mode.'
    ENDIF
ENDIF

IF TOTAL(STRMATCH(valid_parameters, selected_parameter, /FOLD_CASE)) EQ 0 THEN BEGIN
    IF TOTAL(STRMATCH(valid_units, selected_parameter)) EQ 0 THEN BEGIN
        PRINT, 'ERROR: Invaild parameter selected!'
        PRINT, 'Please select from one of the following parameters:'
        FOR i=0, (N_ELEMENTS(valid_parameters) -1) DO PRINT, '   '+valid_parameters[i]
        PRINT, 'Alternatively, you can select one of the valid units to convert FROM.'
        MESSAGE, selected_parameter+' is unknown'
    ENDIF
ENDIF

IF TOTAL(STRMATCH(valid_units, requested_units)) EQ 0 THEN BEGIN
    PRINT, 'ERROR: Invaild units requested!'
    PRINT, 'Please select from one of the following output units:'
    FOR i=0, (N_ELEMENTS(valid_units) -1) DO PRINT, '   '+valid_units[i]
    MESSAGE, requested_units+' is unknown'
ENDIF

;####################################
;Selecting which set of tables to use
;####################################
;Converting from the NUWT reported units for a given parameter
IF TOTAL(STRMATCH(['amplitude', 'amp'], selected_parameter, /FOLD_CASE)) EQ 1 THEN BEGIN
    input_units = nuwt_amp_units
    lookup_table = convert_dist_table
    use_index = index_dist_units
    IF user_res_override THEN rescale = amp_to_pxls ELSE rescale = 1.0
ENDIF ELSE IF TOTAL(STRMATCH(['frequency', 'freq'], selected_parameter, /FOLD_CASE)) EQ 1 THEN BEGIN
    input_units = nuwt_freq_units
    lookup_table = convert_freq_table
    use_index = index_freq_units
    IF user_cad_override THEN rescale = freq_to_timesteps ELSE rescale = 1.0
ENDIF ELSE IF TOTAL(STRMATCH(['period'], selected_parameter, /FOLD_CASE)) EQ 1 THEN BEGIN
    input_units = nuwt_period_units
    lookup_table = convert_time_table
    use_index = index_time_units
    IF user_cad_override THEN rescale = 1.0/freq_to_timesteps ELSE rescale = 1.0
ENDIF ELSE IF TOTAL(STRMATCH(['td_distance', 'td_dist'], selected_parameter, /FOLD_CASE)) EQ 1 THEN BEGIN
    input_units = 'pixels'
    lookup_table = convert_dist_table
    use_index = index_dist_units
    rescale = 1.0
ENDIF ELSE IF TOTAL(STRMATCH(['td_timestep', 'td_time'], selected_parameter, /FOLD_CASE)) EQ 1 THEN BEGIN
    input_units = 'timesteps'
    lookup_table = convert_time_table
    use_index = index_time_units
    rescale = 1.0
ENDIF ELSE IF TOTAL(STRMATCH(['enbw'], selected_parameter, /FOLD_CASE)) EQ 1 THEN BEGIN
    IF cad_unknown THEN input_units = 'timesteps^{-1}' ELSE input_units = 'Hz'
    IF requested_units EQ 'mHz' THEN requested_units = 'Hz'
    lookup_table = convert_freq_table
    use_index = index_freq_units
    IF user_cad_override THEN rescale = slit_located.cad ELSE rescale = 1.0
ENDIF

;Converting between units of the same type
IF TOTAL(STRMATCH(['pixels', 'arcsec', 'm', 'km', 'Mm'], selected_parameter, /FOLD_CASE)) EQ 1 THEN BEGIN
    input_units = selected_parameter
    lookup_table = convert_dist_table
    use_index = index_dist_units
    rescale = 1.0
ENDIF ELSE IF TOTAL(STRMATCH(['timesteps^{-1}', '1/timesteps', 'Hz', 's^{-1}', '1/s', 'mHz'], selected_parameter, /FOLD_CASE)) EQ 1 THEN BEGIN
    input_units = selected_parameter
    lookup_table = convert_freq_table
    use_index = index_freq_units
    rescale = 1.0
ENDIF ELSE IF TOTAL(STRMATCH(['timesteps', 's', 'min', 'hr', 'ks'], selected_parameter, /FOLD_CASE)) EQ 1 THEN BEGIN
    input_units = selected_parameter
    lookup_table = convert_time_table
    use_index = index_time_units
    rescale = 1.0
ENDIF

output_units = requested_units

;Enforcing only conversions we have the correct information to perform
IF res_unknown THEN BEGIN
    IF input_units EQ 'pixels' AND requested_units NE 'pixels' THEN BEGIN
        w_txt = 'WARNING: "res" is unknown! Conversion FROM units of "pixels" is disabled. Output units set to pixels'
        IF NOT KEYWORD_SET(quiet) THEN PRINT, w_txt
        output_units = 'pixels'
    ENDIF ELSE IF input_units NE 'pixels' AND requested_units EQ 'pixels' THEN BEGIN
        w_txt = 'WARNING: "res" is unknown! Conversion TO units of "pixels" is disabled. Output units set to '+input_units
        IF NOT KEYWORD_SET(quiet) THEN PRINT, w_txt
        output_units = input_units
    ENDIF
ENDIF

IF cad_unknown THEN BEGIN
    ;converting freq
    IF input_units EQ 'timesteps^{-1}' AND requested_units NE 'timesteps^{-1}' THEN BEGIN
        w_txt = 'WARNING: "cad" is unknown! Converting FROM units of "timesteps^{-1}" is disabled. Output units set to timesteps^{-1}'
        IF NOT KEYWORD_SET(quiet) THEN PRINT, w_txt
        output_units = 'timesteps^{-1}'
    ENDIF ELSE IF input_units NE 'timesteps^{-1}' AND requested_units EQ 'timesteps^{-1}' THEN BEGIN
        w_txt = 'WARNING: "cad" is unknown! Converting TO units of "timesteps^{-1}" is disabled. Output units set to '+input_units
        IF NOT KEYWORD_SET(quiet) THEN PRINT, w_txt
        output_units = input_units
    ENDIF

    ;converting periods
    IF input_units EQ 'timesteps' AND requested_units NE 'timesteps' THEN BEGIN
        w_txt = 'WARNING: "cad" is unknown! Converting FROM units of "timesteps" is disabled. Output units set to timesteps'
        IF NOT KEYWORD_SET(quiet) THEN PRINT, w_txt
        output_units = 'timesteps'
    ENDIF ELSE IF input_units NE 'timesteps' AND requested_units EQ 'timesteps' THEN BEGIN
        w_txt = 'WARNING: "cad" is unknown! Converting TO units of "timesteps" is disabled. Output units set to '+input_units
        IF NOT KEYWORD_SET(quiet) THEN PRINT, w_txt
        output_units = input_units
    ENDIF
ENDIF

;@convert#######################################################################
;PERFORM THE ACTUAL UNIT CONVERSION
;###############################################################################
conversion_factor = rescale*lookup_table[use_index[input_units], use_index[output_units]]

RETURN, conversion_factor
END
