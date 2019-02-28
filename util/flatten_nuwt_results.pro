;+
;NAME: FLATTEN_NUWT_RESULTS
;
;PURPOSE:
;   Flattens the wave parameters found by NUWT into a set of 1D arrays. This
;   makes certain kinds of plots easier to make (and others perhaps harder).
;
;INPUTS:
;   savfile - [Optional] filename and path for the nuwt_results file to be flattened
;
;OPTIONAL INPUTS:
;   slitnum - virtual slit number to plot. By default, will show slit number "0".
;             Note: can also set to 'all' to plot a combined histogram of all slits.
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
;   /final_units - [DEPRECIATED] (will be removed in the near future)
;   amp_units - desired output units for the wave amplitudes. May select from
;               'pixels', 'arcsec', 'm', 'km', and 'Mm'. Default is 'km'
;   freq_units - desired output units for the wave frequencies. May select from
;                'timesteps^{-1}', 'Hz', and mHz'. Default is 'Hz'
;   period_units - desired output units for the wave periods. May select from
;                  'timesteps', 's', 'min', and 'hr'. Default is 's'
;   min_user_flag - minimum user quality flag to plot. These flags are set with
;                   'set_nuwt_qual_flags.pro'. Default is -1000 (all waves)
;   max_user_flag - maximum user quality flag to plot (see above). Default is 1000
;   min_auto_flag - minimum auto quality flag to plot. These flags are set with
;                   'set_nuwt_qual_flags.pro'. Default is -1000 (all waves)
;   max_auto_flag - maximum auto quality flag to plot (see above). Default is 1000
;   /use_temp_common_blocks - if set, will use the temporary common blocks rather
;                             than the all_nuwt_dat block. Use with care!
;
;OUTPUTS:
;   params_out - stucture with 1D arrays containing all of the wave parameters
;
;HISTORY: Name---------Date---------Description
;         M Weberg  ?? AUG, 2017 - Initial coding
;         M Weberg  27 FEB, 2018 - Updated to include more unit conversion options
;
;TO-DO / LIMITATIONS:
;   - none at the moment
;-

FUNCTION FLATTEN_NUWT_RESULTS, savfile, slitnum=slitnum, $
                               use_temp_common_blocks=use_temp_common_blocks, $
                               res=res, cad=cad, km_per_arcsec=km_per_arcsec, final_units=final_units, $
                               amp_units=amp_units, $
                               freq_units=freq_units, $
                               period_units=period_units, $
                               min_user_flag=min_user_flag, max_user_flag=max_user_flag,$
                               min_auto_flag=min_auto_flag, max_auto_flag=max_auto_flag, $
                               vel_amp_mode=vel_amp_mode, use_refit=use_refit

COMPILE_OPT IDL2
;###############################################################################
;Setting default values
;###############################################################################
IF N_ELEMENTS(slitnum) EQ 0 THEN slitnum = 0
IF N_ELEMENTS(slitnum) GT 1 THEN slitnum = slitnum[0]

IF NOT KEYWORD_SET(min_user_flag) THEN min_user_flag = -1000
IF NOT KEYWORD_SET(max_user_flag) THEN max_user_flag = 1000
IF min_user_flag GT max_user_flag THEN BEGIN
    print, 'WARNING: min_user_flag > max_user_flag!'
    print, '   max_user_flag has been reset to equal min_user_flag'
    max_user_flag = min_user_flag
ENDIF

IF NOT KEYWORD_SET(min_auto_flag) THEN min_auto_flag = -1000
IF NOT KEYWORD_SET(max_auto_flag) THEN max_auto_flag = 1000
IF min_user_flag GT max_user_flag THEN BEGIN
    print, 'WARNING: min_auto_flag > max_auto_flag!'
    print, '   max_auto_flag has been reset to equal min_auto_flag'
    max_auto_flag = min_auto_flag
ENDIF

;###############################################################################
;LOADING DATA AND SELECTING CORRECT SET OF RESULTS
;###############################################################################
IF N_ELEMENTS(savfile) EQ 0 THEN BEGIN
    ;[default] load the data from NUWT COMMON blocks
    COMMON all_nuwt_dat, nuwt_located, nuwt_threads, nuwt_fft_spec, nuwt_fft_peaks
    COMMON nuwt_meta_dat, nuwt_meta, nuwt_bulk_stats
ENDIF ELSE BEGIN
    RESTORE, savfile
    print, 'NUWT results loaded from file but have NOT been restore to the COMMON blocks'
    IF KEYWORD_SET(nuwt_fft_stats) THEN BEGIN
        ;needed for backwards compatibility with older NUWT runs (before Feb 2018)
        nuwt_fft_peaks = nuwt_fft_stats
    ENDIF
    use_temp_common_blocks = 0
ENDELSE

IF NOT KEYWORD_SET(use_temp_common_blocks) THEN BEGIN
    ;[Default] For running after the main NUWT program as been run
    IF typename(slitnum) EQ 'STRING' THEN BEGIN
        ;Concat results from all data slits and plot them together
        total_num_slits = n_elements(nuwt_fft_peaks)

        slit_located = nuwt_located[0]
        slit_fft_peaks = nuwt_fft_peaks[0]

        slit_threads = LIST()
        FOR s=0,(total_num_slits-1) DO BEGIN
            IF s GT 0 THEN BEGIN
                slit_fft_peaks = [temporary(slit_fft_peaks), nuwt_fft_peaks[s]]
            ENDIF

            temp_threads = nuwt_threads[s]
            load_th_from_temp = N_ELEMENTS(nuwt_threads[s])
            FOR h=0, (load_th_from_temp-1) DO BEGIN
                slit_threads.ADD, temp_threads[h]
            ENDFOR
        ENDFOR
    ENDIF ELSE BEGIN
        ;Typical use, plotting just one slit of data
        IF slitnum GT n_elements(nuwt_fft_peaks) THEN BEGIN
            last_slit = n_elements(nuwt_fft_peaks)-1
            print, 'WARNING! slit number '+strtrim(slitnum, 2)+' does not exist!'
            print, '   slitnum set to '+strtrim(last_slit, 2)+' (last set of results)'
            slitnum = last_slit
        ENDIF

        slit_located = nuwt_located[slitnum]
        slit_threads = nuwt_threads[slitnum]
        slit_fft_peaks = nuwt_fft_peaks[slitnum]

        IF KEYWORD_SET(use_refit) THEN BEGIN
            COMMON refit_dat, refit_spec, refit_stats
            ;Quick and dirty swap in structure. use with care!!!
            slit_fft_peaks = refit_stats[slitnum]
        ENDIF

    ENDELSE
ENDIF ELSE BEGIN
    ;Used when running as part of the main NUWT data processor (use with care!)
    COMMON located_dat, located
    COMMON threads_dat, threads
    COMMON fft_results_dat, fft_spec, fft_peaks

    slit_located = located
    slit_threads = threads
    slit_fft_peaks = fft_peaks
ENDELSE

n_threads = n_elements(slit_threads)

;@units#########################################################################
;UNIT CONVERSIONS AND LABELS (AUTOMATED METHODS)
;###############################################################################
;If not given, these parameter will be loaded from the NUWT results where needed
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 0.0
IF NOT KEYWORD_SET(res) THEN res = 0.0
IF NOT KEYWORD_SET(cad) THEN cad = 0.0

IF NOT KEYWORD_SET(amp_units) THEN amp_units = 'km'
IF NOT KEYWORD_SET(freq_units) THEN freq_units = 'Hz'
IF NOT KEYWORD_SET(period_units) THEN period_units = 's'

IF KEYWORD_SET(use_temp_common_blocks) THEN use_temp_common_blocks = 1 ELSE use_temp_common_blocks = 0
amp_dx = convert_nuwt_units('amp', amp_units, res=res, km_per_arcsec=km_per_arcsec, output_units=output_amp_units, $
                            use_temp_common_blocks=use_temp_common_blocks)
freq_dt = convert_nuwt_units('freq', freq_units, cad=cad, output_units=output_freq_units, $
                             use_temp_common_blocks=use_temp_common_blocks)
period_dt = convert_nuwt_units('period', period_units, cad=cad, output_units=output_period_units, $
                               use_temp_common_blocks=use_temp_common_blocks)
enbw_dt = convert_nuwt_units('enbw', 'Hz', cad=cad, output_units=output_enbw_units, $
                             use_temp_common_blocks=use_temp_common_blocks)

;ALWAYS gives psd in time units of [s] or [timesteps] (NEVER [min], [hr], or [ks])
;However, the distance component changes with the selected amplitude units
pow_dp = (amp_dx^2)/enbw_dt
psd_time_unit = 's'
IF period_units EQ 'timesteps' THEN psd_time_unit = 'timestep'
IF output_period_units EQ 'timesteps' OR output_period_units EQ 's' THEN BEGIN
    vel_pow_dp = 1.0
ENDIF ELSE IF output_period_units EQ 'min' THEN BEGIN
    vel_pow_dp = (1.0/60.0)^2
ENDIF ELSE IF output_period_units EQ 'hr' THEN BEGIN
    vel_pow_dp = (1.0/3600.0)^2
ENDIF ELSE IF output_period_units EQ 'ks' THEN BEGIN
    vel_pow_dp = (1.0/1000.0)^2
ENDIF

;If NUWT res & cad values were used, load them for inclusion in the output structure
output_duration_units = 's'
IF NOT KEYWORD_SET(cad) THEN BEGIN
    ;Needed for calculating the thread durations (ALWAYS in [s] or [timesteps])
    cad = slit_located.cad
    IF slit_located.units_dt EQ 'timesteps' THEN output_duration_units = 'timesteps'
ENDIF
IF NOT KEYWORD_SET(res) THEN res = slit_located.res
IF NOT KEYWORD_SET(cad) THEN km_per_arcsec = slit_located.km_per_arcsec

output_vel_amp_units = amp_units+' '+period_units+'^{-1}''
output_pow_units = amp_units+'^{2} '+psd_time_unit
output_vel_pow_units = amp_units+'^{2} '+psd_time_unit+'^{-1}'

IF KEYWORD_SET(final_units) THEN BEGIN
    ;remove once code is fully refactored
    print, "Notice: the /final_units keyword option is unnecessary and has been depreciated!"
ENDIF

;###############################################################################
;FIND VALID PEAK WAVES AND CONVERT TO PROPER UNITS
;###############################################################################
; threads_with_waves = where(slit_fft_peaks.num_signif_peaks GT 0, /NULL)
threads_with_waves = where(slit_fft_peaks.num_saved_waves GT 0, /NULL)
IF N_ELEMENTS(threads_with_waves) GT 0 THEN BEGIN
    num_th_waves = slit_fft_peaks[threads_with_waves].num_saved_waves
ENDIF ELSE BEGIN
    num_th_waves = intarr(n_threads)
ENDELSE

;Cap the number of allowed waves to four (or whatever is the limit of the current results)
max_num_waves_per_thread = n_elements(slit_fft_peaks[0].peak_amplitude)
loc_too_many_waves = where(num_th_waves GT max_num_waves_per_thread, /NULL)
IF n_elements(loc_too_many_waves) GT 0 THEN BEGIN
    num_th_waves[loc_too_many_waves] = max_num_waves_per_thread
ENDIF

;Loop over all valid wave results (and filtering as we go)
found_first_good_wave = 0
FOR tt=0L, n_elements(threads_with_waves)-1 DO BEGIN
    temp_thread = slit_threads[threads_with_waves[tt]]
    temp_fft_peaks = slit_fft_peaks[threads_with_waves[tt]]
    th_user_flag = temp_fft_peaks.user_qual_flag
    th_auto_flag = temp_fft_peaks.auto_qual_flag
    th_bin_flags = temp_thread.bin_flags
    th_len = temp_thread.length
    num_good_pnts = N_ELEMENTS(WHERE(th_bin_flags EQ 2, /NULL))
    num_gap_pnts = N_ELEMENTS(WHERE(th_bin_flags EQ 0, /NULL))
    num_filled_pnts = th_len - num_good_pnts
    th_percent_gaps = (100.0*num_gap_pnts)/th_len
    IF (th_user_flag GE min_user_flag) AND (th_user_flag LE max_user_flag) AND $
       (th_auto_flag GE min_auto_flag) AND (th_auto_flag LE max_auto_flag) THEN BEGIN
        FOR ww=0, num_th_waves[tt]-1 DO BEGIN
            IF NOT found_first_good_wave THEN BEGIN
                peak_amps = [temp_fft_peaks.peak_amplitude[0]]
                peak_periods = [1.0/temp_fft_peaks.peak_freq[0]]
                peak_freqs = [temp_fft_peaks.peak_freq[0]]
                peak_vel_amps = [temp_fft_peaks.peak_vel_amp[0]]
                peak_powers = [temp_fft_peaks.peak_power[0]]
                peak_enbw = [temp_fft_peaks.enbw]
                peak_th_lens = [temp_thread.length]
                peak_percent_gaps = [th_percent_gaps]
                peak_wave_orders = [1]
                peak_th_index = [threads_with_waves[tt]]
                peak_user_flags = [th_user_flag]
                peak_auto_flags = [th_auto_flag]
                found_first_good_wave = 1
            ENDIF ELSE BEGIN
                peak_amps = [peak_amps, temp_fft_peaks.peak_amplitude[ww]]
                peak_periods = [peak_periods, 1.0/temp_fft_peaks.peak_freq[ww]]
                peak_freqs = [peak_freqs, temp_fft_peaks.peak_freq[ww]]
                peak_vel_amps = [peak_vel_amps, temp_fft_peaks.peak_vel_amp[ww]]
                peak_powers = [peak_powers, temp_fft_peaks.peak_power[ww]]
                peak_enbw = [peak_enbw, temp_fft_peaks.enbw]
                peak_th_lens = [peak_th_lens, temp_thread.length]
                peak_percent_gaps = [peak_percent_gaps, th_percent_gaps]
                peak_wave_orders = [peak_wave_orders, ww+1]
                peak_th_index = [peak_th_index, threads_with_waves[tt]]
                peak_user_flags = [peak_user_flags, th_user_flag]
                peak_auto_flags = [peak_auto_flags, th_auto_flag]
            ENDELSE
        ENDFOR
    ENDIF
ENDFOR

IF n_elements(threads_with_waves) GT 0 THEN BEGIN
    ;Apply the unit scaling (easier to read and track if it is all in one place)
    peak_amps = peak_amps*amp_dx
    peak_periods = peak_periods*period_dt
    peak_freqs = peak_freqs*freq_dt
    peak_vel_amps = peak_vel_amps*amp_dx ;CURRENTLY OVERWRITTEN IN THE SECTION BELOW!
    peak_powers = peak_powers*pow_dp
    peak_enbw = peak_enbw*enbw_dt

    peak_vel_amps = 2*!PI*peak_amps/peak_periods
    peak_durations = peak_th_lens*cad

    IF KEYWORD_SET(vel_amp_mode) THEN BEGIN
        peak_vel_amps = peak_amps/period_dt
        peak_amps = peak_vel_amps*peak_periods/(2.0*!PI)
    ENDIF

    peak_vel_pows = (peak_vel_amps^2)*vel_pow_dp/peak_enbw ;Velocity Power Spectral Density (PSD)
    num_all_waves = n_elements(peak_amps)
ENDIF ELSE BEGIN
    ;for the case of no waves found
    peak_amps = [1.0]
    peak_periods = [1.0]
    peak_freqs = [1.0]
    peak_vel_amps = [1.0]
    peak_powers = [1.0]
    peak_enbw = [1.0]
    peak_th_lens = [1.0]
    peak_percent_gaps = [100.0]
    peak_wave_orders = [0]
    peak_th_index = [0]
    peak_user_flags = [-1000]
    peak_auto_flags = [-1000]

    peak_durations = [1.0]
    peak_vel_pows = [1.0]
    num_all_waves = 0
ENDELSE

;##### OUTPUTTING ALL PARAMETER VALUES
;outputs values to structure for manual inspection
params_out = {nuwt_run_date:nuwt_meta.run_date, $
              slitnum:slitnum, $
              num_threads:n_threads, $
              num_waves:num_all_waves, $
              res:res, cad:cad, km_per_arcsec:km_per_arcsec, $
              amps:peak_amps, amp_units:output_amp_units, $
              amp_dx:amp_dx, $
              periods:peak_periods, period_units:output_period_units, $
              period_dt:period_dt, $
              freqs:peak_freqs, freq_units:output_freq_units, $
              freq_dt:freq_dt, $
              vel_amps:peak_vel_amps, vel_amp_units:output_vel_amp_units, $
              pows:peak_powers, pow_units:output_pow_units, $
              pow_dp:pow_dp, $
              vel_pows:peak_vel_pows, vel_pow_units:output_vel_pow_units, $
              vel_pow_dp:vel_pow_dp, $
              enbw:peak_enbw, enbw_units:output_enbw_units, $
              durations:peak_durations, duration_units:output_duration_units, $
              wave_orders:peak_wave_orders, $
              th_len:peak_th_lens, $
              percent_gaps:peak_percent_gaps, $
              th_indices:peak_th_index, $
              user_flags:peak_user_flags, $
              auto_flags:peak_auto_flags}

RETURN, params_out
END
