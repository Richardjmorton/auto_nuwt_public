;+
;NAME: PLOT_NUWT_FFT_RESULTS
;
;PURPOSE:
;   Plots the FFT wave results from 'nuwt_apply_fft.pro' for each thread as a
;   seperate multi-panel plot in a single PDF. Also prints various parameters
;   for diagnostic purposes.
;
;INPUTS:
;   None directly. Loads in the common blocks containing the results from
;   'locate_things.pro' (not currently used), 'follow_thread.pro', and
;   'nuwt_apply_fft.pro'
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
;   dist_units - string indicating what units to use on the distance axis of plot
;                Defaults to units of 'arcsec' if 'res' is known or 'pixels' if
;                'res' is unknown.
;   time_units - string indicating what units to use on the time axis of plots.
;                Defaults to units of 's' if 'cad' is known or 'timesteps' if
;                'cad' is unknown.
;   VVV_units - string with the plot units for amp, freq, or period (vel_amp and
;               power units are automatically determined). The defaults are
;               [km] for amp_units, [Hz] for freq_units, and [s] for period_units
;   /final_units - [DEPRECIATED] (will be removed in the near future)
;   first_thread - index of the first thread to plot. Default is 0
;   last_thread - index of the last thread to plot. Default is (num_threads-1)
;   plot_indices - indices of the threads to plot (overides first_thread and
;                  last_thread). This allows more complicated plotting
;   slitnum - virtual slit number for header text. By defualt assumes 1
;   header - custom header text. Useful for identifying the source data.
;            Defaults to 'NUWT FFT Results'
;   run_tag - custom string that will be appended to the end of the header.
;             Normally used to keep track of differnt runs of the same dataset.
;             There is no default run_tag string.
;   save_folder - folder in which to save the plots. Defaults to the user's home folder.
;                 Will also append a '/' to the end if not included.
;   filename - name for the output PDF file. Default is 'NUWT_FFT_results'
;
;OUTPUTS:
;   PDF_file - multi-page PDF containing FFT results and diagnostics for each
;              thread / set of waves found. If there are multiple significant
;              wave results, the program will plot the combined waveform of
;              all of the waves with power above the significance threshold.
;
;HISTORY: Name---------Date---------Description
;         M Weberg  30 JUNE, 2016  Initial coding
;         M Weberg  26  AUG, 2016  Expanded options and updated plot format
;                                  - Added basic unit conversions
;                                  - Will now plot / list the four largest waves
;                                  - Can not be passed an arbitary header text
;                                  - Other visual tweaks and small improvements
;         M Weberg  14 SEPT, 2016  reconfigured to use NUWT master COMMON block
;         M Weberg  ??  MAR, 2017  Now can load unit information directly from
;                                  NUWT common block structures
;
;TO-DO / LIMITATIONS:
;   - Automatic metadata handling to denote dataset and source of observations.
;     (currently managed via user defined text strings)
;-

FUNCTION TEST_FFT_WAVE_FUNCTION, amplitude_in, freq_in, phase_in, total_len, dt=dt
    IF NOT KEYWORD_SET(dt) THEN dt = 1.0
    ;print, 'input wave parameters: amp =', amplitude_in, ' freq =', freq_in, ' phase =', phase_in, ' total len', total_len
    wave_func = amplitude_in*cos(2.0*!PI*freq_in*findgen(total_len)*dt+phase_in)
    ;print, 'num vals returned =', n_elements(wave_func)
    ;if n_elements(wave_func) EQ 1 THEN print, 'wave val =', wave_func
RETURN, wave_func
END

PRO PLOT_NUWT_FFT_RESULTS, res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                           dist_units=dist_units, time_units=time_units, $
                           amp_units=amp_units, freq_units=freq_units, $
                           period_units=period_units, $
                           final_units=final_units, $

                           show_power_errors=show_power_errors, $
                           vel_amp_mode=vel_amp_mode, amp_scale=amp_scale, $
                           ref_waves=ref_waves, use_refit=use_refit, $
                           sim_res=sim_res, sim_cad=sim_cad, sim_km_per_arcsec=sim_km_per_arcsec, $

                           first_thread=first_thread, last_thread=last_thread, $
                           plot_indices=plot_indices, $
                           slitnum=slitnum, header=header, run_tag=run_tag, $
                           save_folder=save_folder, filename=filename

;###############################################################################
;Setting default values
;###############################################################################
IF NOT KEYWORD_SET(amp_scale) THEN amp_scale = 1.0
IF n_elements(slitnum) EQ 0 THEN slitnum = 0
IF NOT KEYWORD_SET(header) THEN header = 'NUWT FFT Results'
IF KEYWORD_SET(run_tag) THEN header = header + ' - '+run_tag
IF NOT KEYWORD_SET(filename) THEN filename = 'NUWT_FFT_results'
IF NOT KEYWORD_SET(save_folder) THEN CD, CURRENT=save_folder ;i.e. defaults to current folder
IF strlen(save_folder) GT 0 AND NOT save_folder.endswith('/') THEN save_folder = save_folder+'/'

wave_order_colors = ['red', 'orange', 'purple', 'violet']

update_progress_interval = 5.0 ;[percent]

;###############################################################################
;LOADING DATA AND SELECTING CORRECT SET OF RESULTS
;###############################################################################
COMMON all_nuwt_dat, nuwt_located, nuwt_threads, nuwt_fft_spec, nuwt_fft_peaks

IF slitnum GT n_elements(nuwt_fft_peaks) THEN BEGIN
    last_slit = n_elements(nuwt_fft_peaks)-1
    print, 'WARNING! slit number '+strtrim(slitnum, 2)+' does not exist!'
    print, '   slitnum set to '+strtrim(last_slit, 2)+' (last set of results)'
    slitnum = last_slit
ENDIF

slit_located = nuwt_located[slitnum]
slit_threads = nuwt_threads[slitnum]
slit_fft_spec = nuwt_fft_spec[slitnum]
slit_fft_peaks = nuwt_fft_peaks[slitnum]

IF KEYWORD_SET(use_refit) THEN BEGIN
    COMMON refit_dat, refit_spec, refit_stats
    ;Quick and dirty swap in structure. use with care!!!
    slit_fft_spec = refit_spec[slitnum]
    slit_fft_peaks = refit_stats[slitnum]
ENDIF

old_fft_peaks_struc = TAG_EXIST(slit_fft_peaks, 'fft_length')

n_threads = n_elements(slit_threads)

IF NOT KEYWORD_SET(first_thread) THEN first_thread = 0
IF NOT KEYWORD_SET(last_thread) THEN last_thread = n_threads-1
IF last_thread GE n_threads THEN last_thread = n_threads-1
IF NOT KEYWORD_SET(plot_indices) THEN BEGIN
    plot_indices = indgen((last_thread - first_thread) + 1) + first_thread
ENDIF

;@units#########################################################################
;UNIT CONVERSIONS AND LABELS (AUTOMATED METHODS)
;###############################################################################
;If not given, these parameter will be loaded from the NUWT results where needed
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 0.0
IF NOT KEYWORD_SET(res) THEN res = 0.0
IF NOT KEYWORD_SET(cad) THEN cad = 0.0

;Requested units for plotting (not always available)
IF NOT KEYWORD_SET(dist_units) THEN dist_units = 'arcsec'
IF NOT KEYWORD_SET(time_units) THEN time_units = 's'
IF NOT KEYWORD_SET(amp_units) THEN amp_units = 'km'
IF NOT KEYWORD_SET(freq_units) THEN freq_units = 'Hz'
IF NOT KEYWORD_SET(period_units) THEN period_units = 's'

dx = convert_nuwt_units('td_dist', dist_units, res=res, km_per_arcsec=km_per_arcsec, output_units=plot_dist_units)
dt = convert_nuwt_units('td_time', time_units, cad=cad, output_units=plot_time_units)
amp_dx = convert_nuwt_units('amp', amp_units, res=res, km_per_arcsec=km_per_arcsec, output_units=plot_amp_units)
freq_dt = convert_nuwt_units('freq', freq_units, cad=cad, output_units=plot_freq_units)
period_dt = convert_nuwt_units('period', period_units, cad=cad, output_units=plot_period_units)
enbw_dt = convert_nuwt_units('enbw', 'Hz', cad=cad, output_units=output_enbw_units)

;ALWAYS gives psd in time units of [s] or [timesteps] (NEVER [min], [hr], or [ks])
;However, the distance component changes with the selected amplitude units
pow_dp = (amp_dx^2)/enbw_dt
psd_time_unit = 's'
IF period_units EQ 'timesteps' THEN psd_time_unit = 'timestep'
IF plot_period_units EQ 'timesteps' OR plot_period_units EQ 's' THEN BEGIN
    vel_pow_dp = 1.0
ENDIF ELSE IF plot_period_units EQ 'min' THEN BEGIN
    vel_pow_dp = (1.0/60.0)^2
ENDIF ELSE IF plot_period_units EQ 'hr' THEN BEGIN
    vel_pow_dp = (1.0/3600.0)^2
ENDIF ELSE IF plot_period_units EQ 'ks' THEN BEGIN
    vel_pow_dp = (1.0/1000.0)^2
ENDIF

plot_pow_units = amp_units+'^{2} '+psd_time_unit
plot_vel_pow_units = amp_units+'^{2} '+psd_time_unit+'^{-1}'

;Needed for plotting correctly
temp_fft_peaks = slit_fft_peaks[0]
amp_to_pxls = temp_fft_peaks.amp_to_pxls
freq_to_timesteps = temp_fft_peaks.freq_to_timesteps

IF KEYWORD_SET(vel_amp_mode) THEN BEGIN
    amp_dx = amp_dx*freq_dt
    amp_units = '('+ amp_units + ' ' + freq_units +')'
    amp_to_pxls = amp_to_pxls/dt
ENDIF

;###############################################################################
;LOADING SIM DATA AS REF WAVES - CURRENTLY CANNOT CROSS-COMPARE TWO NUWT RUNS!
;###############################################################################
IF KEYWORD_SET(ref_waves) THEN BEGIN
    ;DOUBLE CHECK ALL UNITS!!!!!!!!
    IF NOT keyword_set(sim_res) THEN sim_res = 0.6
    IF NOT keyword_set(sim_cad) THEN sim_cad = 12.0
    IF NOT keyword_set(sim_km_per_arcsec) THEN sim_km_per_arcsec = 725.27

    sim_amp_scale = convert_nuwt_units(ref_waves.amp_units, plot_amp_units, $
                                       res=sim_res, km_per_arcsec=sim_km_per_arcsec)
    sim_period_scale = convert_nuwt_units(ref_waves.period_units, plot_period_units, cad=sim_cad)

    ; lookup_sim_amp_scale = ORDEREDHASH('pixels', 1.0/sim_res, $
    ;                                    'arcsec', 1.0, $
    ;                                    'm', sim_km_per_arcsec*1e3, $
    ;                                    'km', sim_km_per_arcsec, $
    ;                                    'Mm', sim_km_per_arcsec/1e3)
    ; sim_amp_scale = lookup_sim_amp_scale[plot_amp_units]
    ;
    ; lookup_sim_period_scale = ORDEREDHASH('timesteps', 1.0/sim_cad, $
    ;                                       's', 1.0, $
    ;                                       'min', 1.0/60.0, $
    ;                                       'hr', 1.0/3600.0, $
    ;                                       'ks', 1e-3)
    ; sim_period_scale = lookup_sim_period_scale[plot_period_units]

    ; IF n_elements(max_pos_offset) EQ 0 THEN max_pos_offset = 5
    ; IF n_elements(max_bin_offset) EQ 0 THEN max_bin_offset = 5
    max_pos_offset = 5
    max_bin_offset = 5

    ;Compares the NUWT reported values to simulated values ONLY (no comparing different NUWT runs)
    n_ref = n_elements(ref_waves.phase)
    ref_threads = ref_waves

    ;EXTRACTING REFERENCE DATA TO STANDARDIZED STRUCTURES ######################
    ;Extract reference wave start positions and bins (in pixels)
    ref_coords = {start_pos:fltarr(n_ref), end_pos:fltarr(n_ref), $
                  min_pos:fltarr(n_ref), max_pos:fltarr(n_ref), $
                  start_bin:intarr(n_ref), end_bin:intarr(n_ref), $
                  length:intarr(n_ref), paired_index:intarr(n_ref), $
                  paired_flag:intarr(n_ref)}
    ref_coords.paired_index[0:-1] = -1
    ref_coords.start_pos = ref_threads.start_dist/sim_res
    ref_coords.min_pos = reform(ref_threads.box_coords[0,*])
    ref_coords.max_pos = reform(ref_threads.box_coords[2,*])
    ref_coords.end_pos = ref_coords.min_pos + (ref_coords.max_pos - ref_coords.min_pos)/2
    ref_coords.start_bin = reform(ref_threads.box_coords[1,*])
    ref_coords.end_bin = reform(ref_threads.box_coords[3,*])
    ref_coords.length = ref_coords.end_bin - ref_coords.start_bin + 1

    ;Extract wave parameters
    ref_params = {amplitude:fltarr(n_ref), period:fltarr(n_ref), phase:fltarr(n_ref), $
                  velocity_amp:fltarr(n_ref), mean:fltarr(n_ref), slope:fltarr(n_ref), $
                  num_cyc:fltarr(n_ref), num_waves:fltarr(n_ref)}
    ref_params.amplitude = ref_waves.amp*sim_amp_scale ;[km]
    ref_params.period = ref_waves.period*sim_period_scale ;[s]
    ref_params.velocity_amp = 2*!PI*ref_params.amplitude/ref_params.period ;[km/s]
    ref_params.num_cyc = (ref_coords.length*sim_cad*sim_period_scale)/ref_params.period
    ref_params.phase = ref_waves.phase
    loc_shift_phase = where(ref_params.phase GT !PI, /NULL)
    IF n_elements(loc_shift_phase) GT 0 THEN BEGIN
        ;Adjust phase to be in the range of -pi to +pi (as returned by fft)
        ref_params.phase[loc_shift_phase] = ref_params.phase[loc_shift_phase] - 2*!PI
    ENDIF
    ref_params.num_waves[0,-1] = 1
    FOREACH TAG, TAG_NAMES(ref_waves) DO BEGIN
        IF STRLOWCASE(TAG) EQ 'num_waves' THEN BEGIN
            ref_params.num_waves = ref_waves.num_waves
        ENDIF
    ENDFOREACH

    ;EXTRACTING TEST (NUWT) COORDS TO A STANDARDIZED STRUCTURE #################
    ;Extract start positions and bins of each test thread
    test_coords = {start_pos:fltarr(n_threads), end_pos:fltarr(n_threads), $
                   min_pos:fltarr(n_threads), max_pos:fltarr(n_threads), $
                   start_bin:intarr(n_threads), end_bin:intarr(n_threads), $
                   length:intarr(n_threads), paired_index:intarr(n_threads), $
                   paired_flag:intarr(n_threads)}
    test_coords.paired_index[0:-1] = -1
    test_coords.start_bin = slit_threads.start_bin
    test_coords.end_bin = slit_threads.end_bin
    test_coords.length = slit_threads.length
    FOR h=0L, (n_threads-1) DO BEGIN
        test_coords.start_pos[h] = slit_threads[h].pos[test_coords.start_bin[h]]
        test_coords.end_pos[h] = slit_threads[h].pos[test_coords.end_bin[h]]
    ENDFOR
    min_start_pos = test_coords.start_pos - slit_fft_peaks.peak_amplitude[0]*amp_to_pxls
    max_start_pos = test_coords.start_pos + slit_fft_peaks.peak_amplitude[0]*amp_to_pxls
    min_end_pos = test_coords.end_pos - slit_fft_peaks.peak_amplitude[0]*amp_to_pxls
    max_end_pos = test_coords.end_pos + slit_fft_peaks.peak_amplitude[0]*amp_to_pxls
    test_coords.min_pos = min_start_pos < min_end_pos
    test_coords.max_pos = max_start_pos > max_end_pos

    ;PAIRING NUWT LOCATED THREADS WITH SIMULATED DATA ##########################
    paired_indices = {test:indgen(n_threads), ref:intarr(n_threads), $
                      degen:intarr(n_threads), num_ref:intarr(n_threads)}
    FOR h=0L, (n_threads-1) DO BEGIN
        ;[NEW METHOD] Search for overlapping sim and wave boxes
        loc_overlapped_ref = where((((ref_coords.min_pos GE test_coords.min_pos[h]) AND $
                                     (ref_coords.min_pos LE test_coords.max_pos[h])) OR $
                                    ((ref_coords.max_pos GE test_coords.min_pos[h]) AND $
                                     (ref_coords.max_pos LE test_coords.max_pos[h])) OR $
                                    ((ref_coords.min_pos LE test_coords.min_pos[h]) AND $
                                     (ref_coords.max_pos GE test_coords.max_pos[h]))) $
                                   AND $
                                   (((ref_coords.start_bin GE test_coords.start_bin[h]) AND $
                                     (ref_coords.start_bin LE test_coords.end_bin[h])) OR $
                                    ((ref_coords.end_bin GE test_coords.start_bin[h]) AND $
                                     (ref_coords.end_bin LE test_coords.end_bin[h])) OR $
                                    ((ref_coords.start_bin LE test_coords.start_bin[h]) AND $
                                     (ref_coords.end_bin GE test_coords.end_bin[h]))))

        ;count number of ref threads that might be stitched together to make the one test thread
        paired_indices.ref[h] = loc_overlapped_ref[0]
        IF loc_overlapped_ref[0] NE -1 THEN BEGIN
            paired_indices.num_ref[h] = n_elements(loc_overlapped_ref)
            ref_coords.paired_index[loc_overlapped_ref] = h
        ENDIF
    ENDFOR

    test_coords.paired_index = paired_indices.ref

    ;Determine degeneracy (i.e. number of TEST threads paired to the same REF thread)
    FOR h=0L, (n_threads-1) DO BEGIN
        IF paired_indices.ref[h] EQ -1 THEN BEGIN
            ;unpaired TEST threads
            paired_indices.degen[h] = -1
        ENDIF ELSE BEGIN
            IF paired_indices.num_ref[h] GT 1 THEN BEGIN
                ;multiple REF threads connected to same TEST thread
                paired_indices.degen[h] = 0
            ENDIF ELSE BEGIN
                ;one (or more) TEST threads connected to same REF thread
                loc_same_ref = where(paired_indices.ref EQ paired_indices.ref[h])
                num_same_ref = n_elements(loc_same_ref)
                paired_indices.degen[h] = num_same_ref
            ENDELSE
        ENDELSE
    ENDFOR

    loc_ref_paired = where(ref_coords.paired_index GE 0, /NULL)
    loc_ref_unpaired = where(ref_coords.paired_index EQ -1, /NULL)
    ref_coords.paired_flag[loc_ref_paired] = 1

    loc_test_paired = where(test_coords.paired_index GE 0, /NULL)
    test_coords.paired_flag[loc_test_paired] = 1
ENDIF

;###############################################################################
;PLOTTING THE DATA
;###############################################################################
; num_plot_threads = (last_thread - first_thread) + 1
num_plot_threads = N_ELEMENTS(plot_indices)
print, 'Plotting FFT results for '+strtrim(num_plot_threads, 2)+' out of '+strtrim(n_threads, 2)+' threads ...'
print, '--- Please wait, this may take a few minutes ---'
print_progress_threshold = update_progress_interval

fig_fft = window(name='fig_fft', dimensions=[1200, 725], /buffer)

num_plots_done = 0
start_time = systime(/seconds)
FOREACH h, plot_indices DO BEGIN
    ;print, 'Thread number', h
    active_th = slit_threads[h]
    tpos = active_th.pos ;dummy array with postions of the selected thread
    terr = active_th.err_pos ;errors on postions of thread
    th_flags = active_th.bin_flags ;flags for thread data quality

    spec = slit_fft_spec[h] ;results from 'nuwt_apply_fft.pro'
    stats = slit_fft_peaks[h] ;results from 'nuwt_apply_fft.pro'
    tpos = tpos[active_th.start_bin:active_th.end_bin]
    terr = terr[active_th.start_bin:active_th.end_bin]
    trend = spec.trend
    s_len = active_th.length
    IF s_len LE 0 THEN s_len = 1

    IF NOT old_fft_peaks_struc THEN BEGIN
        ;[default] current structure format
        fft_len = spec.fft_length
        signif_test_name = spec.signif_test
        signif_level_val = spec.signif_level
        win_func_name = spec.window_func
        win_param_val = spec.window_param
        cpg_val = spec.cpg
    ENDIF ELSE BEGIN
        ;for plotting old NUWT fft_results (from before 2017-SEPT-19)
        fft_len = stats.fft_length
        signif_test_name = stats.signif_test
        signif_level_val = stats.signif_level
        win_func_name = stats.window_func
        win_param_val = stats.window_param
        cpg_val = stats.cpg
    ENDELSE

    IF KEYWORD_SET(vel_amp_mode) THEN BEGIN
        ;[EXPERIMENTAL] calculate the velocity time series as a means
        ;to directly determine velocity amplitude independently of freq
        ;Central difference method (with forwards/backwards on the ends)
        ; tpos = smooth(tpos, 3, /EDGE_TRUNCATE)
        ; terr = smooth(terr, 3, /EDGE_TRUNCATE)
        dv = (SHIFT(tpos, -1) - SHIFT(tpos, 1))*1.0/(2.0*dt)
        dv[0] = (tpos[1]-tpos[0])*1.0/dt ;forwards
        dv[-1] = (tpos[-1]-tpos[-2])*1.0/dt ;backwards
        dv_err = (SHIFT(terr, -1) - SHIFT(terr, 1))*1.0/(2.0*dt)
        dv_err[0] = (terr[1]-terr[0])*1.0/dt
        dv_err[-1] = (terr[-1]-terr[-2])*1.0/dt
        tpos = dv
        terr = dv_err
        trend = trend/dt
        nuwt_amp_type = 'vel. amp.'
        label_time_series = 'Velocity [$'+plot_dist_units+' '+plot_time_units+'^{-1}$]'
        label_residuals = 'Residuals [$'+plot_dist_units+' '+plot_time_units+'^{-1}$]'
    ENDIF ELSE BEGIN
        nuwt_amp_type = 'amplitude'
        label_time_series = 'Distance ['+plot_dist_units+']'
        label_residuals = 'Residuals ['+plot_dist_units+']'
    ENDELSE

    ;Finding wave properties (including any secondary waves)
    total_fit_params = 0
    total_num_waves = 0
    wave_vals = fltarr(s_len)
    tag_color = strarr(n_elements(stats.peak_bin))
    tag_color[*] = 'dark grey'
    wave_amp_ratios = fltarr(4)
    wave_freq_ratios = fltarr(4)
    FOR w=0, 3 DO BEGIN ;Note: only consider the first four waves for now [HARDCODED!]
        IF stats.peak_bin[w] GT -1 THEN BEGIN
            add_wave = TEST_FFT_WAVE_FUNCTION(stats.peak_amplitude[w]*amp_to_pxls*amp_scale, $
                                              stats.peak_freq[w]*freq_to_timesteps, stats.peak_phase[w], s_len)
            wave_vals = wave_vals + add_wave
            total_num_waves += 1
            total_fit_params += 3
            tag_color[w] = 'black'
            wave_amp_ratios[w] = stats.peak_amplitude[w]/stats.peak_amplitude[0]
            wave_freq_ratios[w] = stats.peak_freq[w]/stats.peak_freq[0]
        ENDIF
    ENDFOR

    apod_pos = (tpos-trend)*spec.apod_window + trend
    ts_xvals = (findgen(s_len) + active_th.start_bin)*dt
    residuals = tpos - (wave_vals + trend)

    ;Extracting various bits of information for printing
    num_good_pnts = n_elements(th_flags[where(th_flags EQ 2, /NULL)])
    num_filled_pnts = s_len - num_good_pnts
    percent_filled = (100.0*num_filled_pnts)/s_len
    num_freq_bins = n_elements(spec.freq)
    num_fft_pad_zeros = fft_len - s_len
    IF num_fft_pad_zeros GT 0 THEN tag_pad = 'YES' ELSE tag_pad = 'no'

    mean_trend = mean(trend)*dx
    slope_trend = ((trend[-1] - trend[0])/(ts_xvals[-1] - ts_xvals[0]))*(dx/dt)

    ;##### UPPER LEFT SUBPLOT #####
    ;Time series
    ;##############################
    ;Plot the time series of thread positions along with apodized values and estimated wave
    plt_pos = errorplot(ts_xvals, tpos*dx, terr*dx, name='Observed', $
                        color='blue', symbol='o', /sym_filled, linestyle='none', errorbar_color='blue',$
                        position=[0.06, 0.52, 0.35, 0.9], /current)

    plt_apod_pos = plot(ts_xvals, apod_pos*dx, name='Windowed', $
                        color='sky blue', symbol='x', linestyle='none', /overplot)

    plt_trend = plot(ts_xvals, trend*dx, name='Trend', color='black', symbol='none', linestyle=':', /overplot)

    plt_wave = plot(ts_xvals, (wave_vals + trend)*dx, name='Fit Wave', $
                    color='black', symbol='none', linestyle='-', /overplot)

    leg_pos = legend(target=[plt_pos, plt_apod_pos, plt_trend, plt_wave], /normal, $
                     position=[0.485, 0.9], linestyle='none', transparency=100, font_size=12)
    plt_pos.title = 'Time Series'
    plt_pos.ytitle = label_time_series

    ;##### LOWER LEFT SUBPLOT #####
    ;Residuals
    ;##############################
    ;Plot the residuals of the time series (i.e. what is left after subtracting the 'fit' function)
    plt_residuals = plot(ts_xvals, residuals*dx, color='blue', symbol='none', linestyle='-', $
                         position=[0.06, 0.1, 0.35, 0.48], /current)
    plt_residuals.ytitle = label_residuals
    plt_residuals.xtitle = 'Time ['+plot_time_units+']'

    ;##### UPPER RIGHT SUBPLOT ####
    ;FFT power spectrum with the significance limit
    ;##############################
    plt_pow = plot(spec.freq*freq_dt, spec.power*amp_dx^2, name='FFT results', $
                   color='black', symbol='none', linestyle='-', $
                   position=[0.56, 0.52, 0.85, 0.9], /current, /stairstep)

    IF KEYWORD_SET(show_power_errors) THEN BEGIN
        plt_errs_pow = errorplot(spec.freq*freq_dt, spec.power*amp_dx^2, spec.err_power*amp_dx^2, $
                                 color='black', symbol='none', linestyle='none', $
                                 errorbar_color='grey', /overplot)
    ENDIF

    plt_signif_pow = plot(spec.freq*freq_dt, spec.signif_vals*amp_dx^2, name='Signif. limit', $
                          color='black', symbol='none', linestyle='--', /overplot)

    plt_primary_peak = plot([stats.peak_freq[0]*freq_dt, stats.peak_freq[0]*freq_dt], $
                            [stats.peak_power[0]*amp_dx^2, stats.peak_power[0]*amp_dx^2], $
                            name='Primary Peak', color=wave_order_colors[0], symbol='D', linestyle='none', /overplot)

    plt_secondary_peak = plot([stats.peak_freq[1]*freq_dt, stats.peak_freq[1]*freq_dt], $
                              [stats.peak_power[1]*amp_dx^2, stats.peak_power[1]*amp_dx^2], $
                              name='Secondary', color=wave_order_colors[1], symbol='*', linestyle='none', /overplot)

    plt_tertiary_peak = plot([stats.peak_freq[2]*freq_dt, stats.peak_freq[2]*freq_dt], $
                             [stats.peak_power[2]*amp_dx^2, stats.peak_power[2]*amp_dx^2], $
                             name='Tertiary', color=wave_order_colors[2], symbol='+', linestyle='none', /overplot)

    plt_quaternary_peak = plot([stats.peak_freq[3]*freq_dt, stats.peak_freq[3]*freq_dt], $
                               [stats.peak_power[3]*amp_dx^2, stats.peak_power[3]*amp_dx^2], $
                               name='Quaternary', color=wave_order_colors[3], symbol='Td', linestyle='none', /overplot)

    leg_fft = legend(target=[plt_pow, plt_signif_pow], /normal, position=[0.985, 0.9], $
                     linestyle='none', transparency=100, font_size=12)
    plt_pow.title = 'FFT'
    plt_pow.ytitle = 'PSD [$'+plot_pow_units+'$]'
    plt_pow.xtitle = 'Frequency [$'+plot_freq_units+'$]'

    ;##### PRINTED TEXT #####
    ;Header and ID text (top of plot window)
    txt_slit_ID = text(0.01, 0.97, 'Slit # '+strtrim(slitnum, 2), font_size=14, font_style=1)
    txt_thread_ID = text(0.08, 0.97, 'Thread # '+strtrim(h, 2), font_size=14, font_style=1)
    txt_header = text(0.225, 0.97, header, font_size=14, font_style=1) ;OLD px-pos = 0.35

    ;Time series text
    txt_trend = text(0.36, 0.67, 'Linear Trend:$\n$'+$
                     '   mean = '+strtrim(string(mean_trend, format='(e10.3)'), 2)+'$\n$'+$
                     '   slope = '+strtrim(string(slope_trend, format='(e10.3)'), 2), font_size=12)
    txt_percent_filled = text(0.36, 0.60, strtrim(percent_filled, 2)+'% filled', font_size=12)
    txt_thread_len = text(0.36, 0.52, strtrim(active_th.length, 2)+' data points:$\n$   '+$
                          strtrim(num_good_pnts, 2)+' good$\n$   '+$
                          strtrim(num_filled_pnts, 2)+' filled', font_size=12)

    ;FFT text
    txt_signif_test = text(0.86, 0.812, '('+signif_test_name+')', font_size=9.5)
    txt_signif_level = text(0.86, 0.782, 'Signif. level: '+string(signif_level_val, format='(f5.3)'), font_size=12)
    txt_window_func = text(0.86, 0.705, 'Window Function:$\n$   '+win_func_name, font_size=12)
    txt_window_param = text(0.86, 0.675, 'Parameter: '+strtrim(string(win_param_val, format='(f6.3)'), 2), font_size=12)
    txt_cpg = text(0.86, 0.645, 'CPG: '+strtrim(string(cpg_val, format='(f7.4)'), 2), font_size=12)
    txt_freq_len = text(0.86, 0.60, strtrim(num_freq_bins, 2)+' freq. bins', font_size=12)
    txt_fft_len = text(0.86, 0.52, strtrim(fft_len, 2)+' FFT input points:$\n$   '+$
                          strtrim(active_th.length, 2)+' data (windowed)$\n$   '+$
                          strtrim(num_fft_pad_zeros, 2)+' zeros (padding)', font_size=12)

    ;Goodness-of-fit (GOF) test
    txt_GOF_header = text(0.36, 0.43, 'GOF Tests', font_size=14, font_style=1)
    txt_KS_test = text(0.36, 0.385, 'Kolmogorov-Smirnov', font_size=12)
    txt_KS_stat = text(0.36, 0.355, '   stat = '+strtrim(stats.KS_stat[total_num_waves], 2), font_size=12)
    txt_KS_prob = text(0.36, 0.325, '   prob = '+strtrim(stats.KS_prob[total_num_waves], 2), font_size=12)
    txt_AD_test = text(0.36, 0.275, 'Anderson-Darling', font_size=12)
    txt_AD_stat = text(0.36, 0.245, '   $A^2$ = '+strtrim(stats.AD_stat[total_num_waves], 2), font_size=12)
    txt_AD_crit = text(0.36, 0.215, '   crit val = '+strtrim(stats.AD_crit[total_num_waves], 2), font_size=12)
    txt_LB_test = text(0.36, 0.165, 'Ljung-Box', font_size=12)
    txt_LB_stat = text(0.36, 0.135, '   $Q$ = '+strtrim(stats.LB_stat[total_num_waves], 2), font_size=12)
    txt_LB_crit = text(0.36, 0.105, '   $\chi^2$ = '+strtrim(stats.LB_chisqrd[total_num_waves], 2), font_size=12)

    txt_result_header = text(0.56, 0.43, 'Results', font_size=14, font_style=1)
    txt_primary_tag = symbol(0.505, 0.375, 'D', sym_size=1.0, sym_color=wave_order_colors[0],  $
                             label_string='Primary$\n$Wave', label_font_size=13.5, $
                             label_font_style=1, label_color=wave_order_colors[0])
    txt_secondary_tag = symbol(0.505, 0.27, '*', sym_size=1.0, sym_color=wave_order_colors[1],  $
                               label_string='Secondary$\n$Wave', label_font_size=13.5, $
                               label_font_style=1, label_color=wave_order_colors[1])
    txt_tertiary_tag = symbol(0.505, 0.165, '+', sym_size=1.0, sym_color=wave_order_colors[2],  $
                              label_string='Tertiary$\n$Wave', label_font_size=13.5, $
                              label_font_style=1, label_color=wave_order_colors[2])
    txt_quaternary_tag = symbol(0.505, 0.06, 'Td', sym_size=1.0, sym_color=wave_order_colors[3],  $
                                label_string='Quaternary$\n$Wave', label_font_size=13.5, $
                                label_font_style=1, label_color=wave_order_colors[3])

    result_txt_y_locs = [0.33, 0.225, 0.12, 0.015]
    FOR i=0, 3 DO BEGIN
        IF NOT KEYWORD_SET(ref_waves) THEN BEGIN
            ;[DEFAULT] displays the peak power along with whichever amp. mode is used
            txt_wave_results = text(0.59, result_txt_y_locs[i], $
                                    'peak PSD = '+strtrim(string(stats.peak_power[i]*amp_dx^2, '(e10.3)'), 2)+' $'+plot_pow_units+'\n$'+$
                                     nuwt_amp_type+' = '+strtrim(stats.peak_amplitude[i]*amp_dx*amp_scale, 2)+' $'+plot_amp_units+'\n$'+$
                                    'frequency = '+strtrim(string(stats.peak_freq[i]*freq_dt, '(e10.3)'), 2)+' $'+plot_freq_units+'\n$'+$
                                    'phase = '+strtrim(stats.peak_phase[i], 2)+' rad', font_size=12, color=tag_color[i])
            IF i GT 0 THEN BEGIN
                txt_param_ratio = text(0.80, result_txt_y_locs[i]+0.042, $
                                       nuwt_amp_type+' = '+strtrim(string(wave_amp_ratios[i], '(f8.4)'), 2)+'$\n$'+$
                                       'frequency = '+strtrim(string(wave_freq_ratios[i], '(f8.4)'), 2),$
                                       font_size=12, color=tag_color[i])
            ENDIF
        ENDIF ELSE BEGIN
            ;Easier to compare to input simulated waves
            nuwt_amp = stats.peak_amplitude[i]*amp_dx
            nuwt_freq = stats.peak_freq[i]*freq_dt
            nuwt_vel_amp = 2*!PI*nuwt_amp*nuwt_freq
            nuwt_phase = stats.peak_phase[i]
            IF KEYWORD_SET(vel_amp_mode) THEN BEGIN
                nuwt_vel_amp = nuwt_amp ;note: the conversion to units of [km/s]] is included in amp_+dx at this time
                IF nuwt_freq LE 0.0 THEN nuwt_amp = 0.0 ELSE nuwt_amp = nuwt_vel_amp/(2*!PI*nuwt_freq)
            ENDIF
            txt_wave_results = text(0.59, result_txt_y_locs[i], $
                                    'amplitude = '+strtrim(nuwt_amp, 2)+' $'+plot_amp_units+'\n$'+$
                                    'frequency = '+strtrim(string(nuwt_freq, '(e10.3)'), 2)+' $'+plot_freq_units+'\n$'+$
                                    'vel. amp. = '+strtrim(nuwt_vel_amp, 2)+' $'+plot_amp_units+' '+plot_period_units+'^{-1}\n$'+$
                                    'phase = '+strtrim(nuwt_phase, 2)+' rad', font_size=12, color=tag_color[i])
        ENDELSE
    ENDFOR

    IF NOT KEYWORD_SET(ref_waves) THEN BEGIN
        ;[Default] Print FFT window information
        txt_fft_header = text(0.80, 0.43, 'FFT Details', font_size=14, font_style=1)
        txt_adj_peaks = text(0.80, 0.39, 'Adjacent peaks allowed? '+stats.adjacent_peaks, font_size=12)
        txt_ratio_header = text(0.80, 0.33, 'Parameter Ratios', font_size=14, font_style=1) ;see above for ratios

        txt_timestamp = text(0.85, 0.01, 'Date & Time Created$\n$'+systime(), font_size=10, color='dark grey')
    ENDIF ELSE BEGIN
        ;Print the simulation input values
        ;Simulated data from 'generate_kink_waves.pro'

        sim_index = test_coords.paired_index[h]
        num_sim_waves = ref_params.num_waves[sim_index]
        txt_sim_header = text(0.80, 0.43, 'Simulation Inputs', font_size=14, font_style=1)
        result_txt_y_locs = [0.33, 0.225, 0.12, 0.015]
        FOR i=0, 3 DO BEGIN
            sim_amp = 0.0
            sim_period = 0.0
            sim_vel_amp = 0.0
            sim_freq = 0.0
            sim_phase = 0.0
            tag_wave_color = 'dark_grey'
            IF (sim_index GE 0) AND (i LT num_sim_waves) THEN BEGIN
                sim_amp = ref_waves.amp[sim_index,i]*sim_amp_scale;[km]
                sim_period = ref_waves.period[sim_index,i]*sim_period_scale ;[s]
                sim_vel_amp = 2*!PI*sim_amp/sim_period ;[km/s]
                sim_freq = 1.0/sim_period
                sim_phase = ref_waves.phase[sim_index,i]
                IF sim_phase GT !PI THEN sim_phase = sim_phase - 2*!PI ;Adjust phase to be in the range of -pi to +pi
                tag_wave_color = 'black'
            ENDIF
            txt_wave_results = text(0.80, result_txt_y_locs[i], $
                                    'amplitude = '+strtrim(sim_amp, 2)+' $'+plot_amp_units+'\n$'+$
                                    'frequency = '+strtrim(string(sim_freq, '(e10.3)'), 2)+' $'+plot_freq_units+'\n$'+$
                                    'vel. amp. = '+strtrim(sim_vel_amp, 2)+' $'+plot_amp_units+' '+plot_period_units+'^{-1}\n$'+$
                                    'phase = '+strtrim(sim_phase, 2)+' rad', font_size=12, color=tag_wave_color)
        ENDFOR
        txt_timestamp = text(0.85, 0.965, 'Date & Time Created$\n$'+systime(), font_size=10, color='dark grey')
    ENDELSE

    ; IF h LT (num_plot_threads-1) THEN BEGIN
    IF num_plots_done LT (num_plot_threads-1) THEN BEGIN
        fig_fft.save, save_folder+filename+'.pdf', page_size=[12.0, 7.25], /append
        fig_fft.erase
    ENDIF ELSE BEGIN
        fig_fft.save, save_folder+filename+'.pdf', page_size=[12.0, 7.25], /append, /close
    ENDELSE

    num_plots_done= num_plots_done + 1

    ;##### PRINTING STATUS AND ESTIMATED TIME REMAINING #####
    ; percent_done = (float(h-first_thread+1)/float(num_plot_threads))*100.0
    percent_done = (float(num_plots_done)/float(num_plot_threads))*100.0
    IF percent_done GE print_progress_threshold THEN BEGIN
        current_time = systime(/seconds)
        ; avg_sec_per_loop = (current_time - start_time)/double(h-first_thread+1)
        ; est_min_remaining = (last_thread-h)*avg_sec_per_loop/60.0
        ; print, strtrim(string(percent_done, format='(F6.2)'), 2)+'% completed ('+strtrim(h-first_thread+1,2)+'/'+strtrim(num_plot_threads, 2)+'), '+$
        ;        'estimated '+strtrim(string(est_min_remaining, format='(F8.2)'), 2)+' min remaining'
        avg_sec_per_loop = (current_time - start_time)/double(num_plots_done)
        est_min_remaining = (num_plot_threads-num_plots_done)*avg_sec_per_loop/60.0
        print, strtrim(string(percent_done, format='(F6.2)'), 2)+'% completed ('+strtrim(num_plots_done,2)+'/'+strtrim(num_plot_threads, 2)+'), '+$
               'estimated '+strtrim(string(est_min_remaining, format='(F8.2)'), 2)+' min remaining'
        print_progress_threshold = print_progress_threshold + update_progress_interval
    ENDIF
; ENDFOR
ENDFOREACH
print, 'Finished plotting FFT results to a multi-page PDF!'
print, 'Save Folder: ', save_folder
print, 'Filename: ', filename+'.pdf'
END
