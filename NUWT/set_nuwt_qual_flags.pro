;+
;NAME: SET_NUWT_QUAL_FLAGS
;
;PURPOSE:
;   Interactively plots the FFT wave results from 'nuwt_apply_fft.pro' for each thread
;   and asks the user to rate the quality of the fit on a three step scale of:
;   0 (bad fit), 1 (questionable fit), & 2 (good fit)
;   Various GOF statistics are given as an aid for selection
;
;INPUTS:
;   None directly. Loads in the NUWT common blocks containing the wave results
;
;OPTIONAL INPUTS:
;   res - spatial resolution of data (in arcsec). Default is 1
;   cad - temporal cadence of the data (in s). Default is 1
;   dist_units - string indicating what units to use for distance in the plot
;                     labels. Defaults to units of 'pixels' unless 'res' is set,
;                     in which case the default is 'arcsec'
;   time_units - string indicating what units to use for time in the plot
;                     labels. Defaults to units of 'steps' unless 'cad' is set,
;                     in which case the default is 's'
;   /final_units - [DEPRECIATED] if set, will assume wave parameters are already in the correct
;                  units of [km] for amplitude and [s^-1] for frequency
;   slitnum - virtual slit number for header text. By defualt assumes 1
;   header - custom header text. Useful for identifying the source data.
;            Defaults to 'NUWT FFT Wave Result'
;
;OUTPUTS:
;   fit_qual - interger array with the user's selected quality flags
;
;HISTORY: Name---------Date---------Description
;         M Weberg  19 SEPT, 2016  Initial coding
;         M Weberg  ?? JULY, 2017  Added automatic flagging based on size of
;                                  data gaps
;
;TO-DO / LIMITATIONS:
;   - Automatic metadata handling to denote dataset and source of observations.
;     (currently managed via user defined text strings)
;   - Better interfacing with COMMON block data to store the results
;-

FUNCTION TEST_FFT_WAVE_FUNCTION, amplitude_in, freq_in, phase_in, total_len, dt=dt
    IF NOT keyword_set(dt) THEN dt = 1.0
    ;print, 'input wave parameters: amp =', amplitude_in, ' freq =', freq_in, ' phase =', phase_in, ' total len', total_len
    wave_func = amplitude_in*cos(2.0*!PI*freq_in*findgen(total_len)*dt+phase_in)
    ;print, 'num vals returned =', n_elements(wave_func)
    ;if n_elements(wave_func) EQ 1 THEN print, 'wave val =', wave_func
RETURN, wave_func
END

PRO SET_NUWT_QUAL_FLAGS, res=res, cad=cad, km_per_arcsec=km_per_arcsec, final_units=final_units, $
                         dist_units=dist_units, time_units=time_units, $
                         amp_units=amp_units, freq_units=freq_units, period_units=period_units, $
                         interactive=interactive, vel_amp_mode=vel_amp_mode, $
                         max_percent_gaps=max_percent_gaps, $
                         quiet=quiet, $
                         slitnum=slitnum, header=header, $
                         use_temp_common_blocks=use_temp_common_blocks, $
                         user_flags=user_flags, auto_flags=auto_flags

COMPILE_OPT IDL2

;Seting default values
IF n_elements(slitnum) EQ 0 THEN slitnum = 0
IF NOT KEYWORD_SET(header) THEN header = 'NUWT FFT Wave Result'
IF KEYWORD_SET(vel_amp_mode) THEN vel_amp_mode = 1 ELSE vel_amp_mode = 0

IF KEYWORD_SET(quiet) THEN quiet = 1 ELSE quiet = 0
IF NOT KEYWORD_SEt(max_percent_gaps) THEN max_percent_gaps = 50.0

;@load_data#####################################################################
;LOADING DATA AND SELECTING CORRECT SET OF RESULTS
;###############################################################################
COMMON all_nuwt_dat, nuwt_located, nuwt_threads, nuwt_fft_spec, nuwt_fft_peaks

IF NOT KEYWORD_SET(use_temp_common_blocks) THEN BEGIN
    use_temp_common_blocks = 0
    ;[Default] For running after the main NUWT program as been run
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
ENDIF ELSE BEGIN
    use_temp_common_blocks = 1
    ;Used when running as part of the main NUWT data processor (use with care!)
    COMMON located_dat, located
    COMMON threads_dat, threads
    COMMON fft_results_dat, fft_spec, fft_peaks

    slit_located = located
    slit_threads = threads
    slit_fft_spec = fft_spec
    slit_fft_peaks = fft_peaks
ENDELSE

n_threads = n_elements(slit_threads)

user_flags = intarr(n_threads) ;optional output
auto_flags = intarr(n_threads) ;optional output

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

IF KEYWORD_SET(use_temp_common_blocks) THEN use_temp_common_blocks = 1 ELSE use_temp_common_blocks = 0
dx = convert_nuwt_units('td_dist', dist_units, res=res, km_per_arcsec=km_per_arcsec, $
                        output_units=plot_dist_units, $
                        use_temp_common_blocks=use_temp_common_blocks, $
                        quiet=quiet)
dt = convert_nuwt_units('td_time', time_units, cad=cad, $
                        output_units=plot_time_units, $
                        use_temp_common_blocks=use_temp_common_blocks, $
                        quiet=quiet)
amp_dx = convert_nuwt_units('amp', amp_units, res=res, km_per_arcsec=km_per_arcsec, $
                            output_units=plot_amp_units, $
                            use_temp_common_blocks=use_temp_common_blocks, $
                            quiet=quiet)
freq_dt = convert_nuwt_units('freq', freq_units, cad=cad, $
                             output_units=plot_freq_units, $
                             use_temp_common_blocks=use_temp_common_blocks, $
                             quiet=quiet)
period_dt = convert_nuwt_units('period', period_units, cad=cad, $
                               output_units=plot_period_units, $
                               use_temp_common_blocks=use_temp_common_blocks, $
                               quiet=quiet)

;Needed for plotting correctly
temp_fft_peaks = slit_fft_peaks[0]
amp_to_pxls = temp_fft_peaks.amp_to_pxls
freq_to_timesteps = temp_fft_peaks.freq_to_timesteps

IF KEYWORD_SET(vel_amp_mode) THEN BEGIN
    amp_dx = amp_dx*freq_dt
    amp_units = '('+ amp_units + ' ' + freq_units +')'
    amp_to_pxls = amp_to_pxls/dt
ENDIF

;@calc_stats####################################################################
;CALCULATES BULK STATS FOR FOR STAT-BASED FLAGGING
;###############################################################################
; print, '##### Setting auto quality flags for all '+strtrim(n_threads, 2)+' threads #####'
CALC_NUWT_BULK_STATS, slitnum=slitnum, vel_amp_mode=vel_amp_mode, $
                 amp_units=amp_units, freq_units=freq_units, period_units=period_units, $
                 bulk_stats_out=slit_bulk_stats, params_out=slit_wave_params, $
                 use_temp_common_blocks=use_temp_common_blocks

;side note: this automatically filters out threads with NO GAPS but also NO SIGNIF. WAVES!
;           Is this really what we should do here...
loc_first_order = WHERE(slit_wave_params.wave_orders EQ 1)
IF N_ELEMENTS(loc_first_order) GT 0 THEN BEGIN
    first_th_indices = slit_wave_params.th_indices[loc_first_order]
    th_percent_gaps = slit_wave_params.percent_gaps[loc_first_order]
    loc_good_fits = WHERE(th_percent_gaps LT max_percent_gaps)
ENDIF

;flag data and update slit_fft_peaks array
IF N_ELEMENTS(loc_good_fits) GT 0 THEN BEGIN
    auto_flags[first_th_indices[loc_good_fits]] = 2 ;"good" fit quality
    slit_fft_peaks.auto_qual_flag = auto_flags
ENDIF

;@plot##########################################################################
;PLOTTING THE DATA AND GETTING FLAGS
;###############################################################################
IF KEYWORD_SET(interactive) THEN BEGIN
    IF n_threads GT 20 THEN BEGIN
        ask_cont = 'yes'
        print, ' ' ;empty line for cleaner printing
        print, 'WARNING! There are '+strtrim(n_threads, 2)+' total threads to be flagged!'
        READ, ask_cont, PROMPT='Are you sure you want to continue? (y/n) '
        ask_cont = strlowcase(ask_cont)
        IF NOT ask_cont.startswith('y') THEN BEGIN
            print, 'Exiting program. Quality flags NOT set.'
            STOP
        ENDIF
    ENDIF

    print, ' ' ;empty line for cleaner printing
    print, '##### Setting user quality flags for all '+strtrim(n_threads, 2)+' threads #####'
    print, '--- Type "q" at any time to quit early ---'

    fig_wave = window(name='fig_wave', window_title='Setting Wave Quality Flags', dimensions=[800, 800], /no_toolbar)

    end_early = 0
    FOR h=0L, (n_threads-1) DO BEGIN
        fig_wave.erase ;clearing the plot window

        active_th = slit_threads[h]
        tpos = active_th.pos ;dummy array with postions of the selected thread
        terr = active_th.err_pos ;errors on postions of thread
        th_flags = active_th.bin_flags ;flags for thread data quality

        spec = slit_fft_spec[h] ;results from 'nuwt_apply_fft.pro'
        stats = slit_fft_peaks[h] ;results from 'nuwt_apply_fft.pro'
        tpos = tpos[active_th.start_bin:active_th.end_bin]
        terr = terr[active_th.start_bin:active_th.end_bin]
        trend = spec.trend
        apod_pos = (tpos-trend)*spec.apod_window + trend

        s_len = active_th.length
        ts_xvals = (findgen(s_len) + active_th.start_bin)*dt

        ;Finding wave properties (including any secondary waves)
        total_fit_params = 0
        total_num_waves = 0
        wave_vals = fltarr(s_len)
        tag_color = strarr(n_elements(stats.peak_bin))
        tag_color[*] = 'dark grey'
        FOR w=0, 3 DO BEGIN ;only consider the first four waves for now
            IF stats.peak_bin[w] GT -1 THEN BEGIN
                add_wave = TEST_FFT_WAVE_FUNCTION(stats.peak_amplitude[w]*amp_to_pxls, $
                                                  stats.peak_freq[w]*freq_to_timesteps, stats.peak_phase[w], s_len)
                wave_vals = wave_vals + add_wave
                total_num_waves += 1
                total_fit_params += 3
                tag_color[w] = 'black'
            ENDIF
        ENDFOR

        residuals = tpos - (wave_vals + trend)

        ;Extracting various bits of information for printing
        num_good_pnts = n_elements(th_flags[where(th_flags EQ 2, /NULL)])
        num_filled_pnts = s_len - num_good_pnts
        percent_filled = (100.0*num_filled_pnts)/s_len

        mean_trend = mean(trend)*dx
        slope_trend = ((trend[-1] - trend[0])/(ts_xvals[-1] - ts_xvals[0]))*(dx/dt)

        ;##### UPPER LEFT SUBPLOT #####
        ;Plot the time series of thread positions along with apodized values and estimated wave
        plt_pos = errorplot(ts_xvals, tpos*dx, terr, name='Observed', $
                            color='blue', symbol='o', /sym_filled, linestyle='none', errorbar_color='blue',$
                            position=[0.12, 0.52, 0.74, 0.9], current=fig_wave)

        plt_apod_pos = plot(ts_xvals, apod_pos*dx, name='Windowed', $
                            color='sky blue', symbol='x', linestyle='none', /overplot)

        plt_trend = plot(ts_xvals, trend*dx, name='Trend', color='black', symbol='none', linestyle=':', /overplot)

        plt_wave = plot(ts_xvals, (wave_vals + trend)*dx, name='Fit Wave(s)', $
                        color='black', symbol='none', linestyle='-', /overplot)

        leg_pos = legend(target=[plt_pos, plt_apod_pos, plt_trend, plt_wave], /normal, $
                         position=[0.95, 0.9], linestyle='none', transparency=100, font_size=12)
        plt_pos.title = 'Time Series'
        plt_pos.ytitle = 'Distance ['+plot_dist_units+']'

        ;##### LOWER LEFT SUBPLOT #####
        ;Plot the residuals of the time series (i.e. what is left after subtracting the 'fit' function)
        plt_residuals = plot(ts_xvals, residuals*dx, color='blue', symbol='none', linestyle='-', $
                             position=[0.12, 0.1, 0.74, 0.48], current=fig_wave)
        plt_residuals.ytitle = 'Residuals ['+plot_dist_units+']'
        plt_residuals.xtitle = 'Time ['+plot_time_units+']'

        ;##### PRINTED TEXT #####
        ;Printing useful information in the plot
        txt_slit_ID = text(0.01, 0.97, 'Slit # '+strtrim(slitnum, 2), font_size=14, font_style=1)
        txt_thread_ID = text(0.12, 0.97, 'Thread # '+strtrim(h, 2), font_size=14, font_style=1)
        txt_header = text(0.35, 0.97, header, font_size=14, font_style=1)

        tx_num_waves = text(0.755, 0.75, strtrim(total_num_waves, 2) + ' total waves', font_size=12)
        txt_trend = text(0.755, 0.67, 'Linear Trend:$\n$'+$
                         '   mean = '+strtrim(mean_trend, 2)+'$\n$'+$
                         '   slope = '+strtrim(string(slope_trend, format='(e10.3)'), 2), font_size=12)
        txt_thread_len = text(0.755, 0.57, strtrim(active_th.length, 2)+' data points:$\n$   '+$
                              strtrim(num_good_pnts, 2)+' good$\n$   '+$
                              strtrim(num_filled_pnts, 2)+' filled', font_size=12)
        txt_percent_filled = text(0.755, 0.52, strtrim(percent_filled, 2)+'% filled', font_size=12)

        txt_GOF_header = text(0.755, 0.43, 'GOF Tests', font_size=14, font_style=1)
        txt_KS_test = text(0.755, 0.38, 'Kolmogorov-Smirnov', font_size=12)
        txt_KS_stat = text(0.755, 0.35, '   stat = '+strtrim(stats.KS_stat[total_num_waves], 2), font_size=12)
        txt_KS_prob = text(0.755, 0.32, '   prob = '+strtrim(stats.KS_prob[total_num_waves], 2), font_size=12)
        txt_AD_test = text(0.755, 0.27, 'Anderson-Darling', font_size=12)
        txt_AD_stat = text(0.755, 0.24, '   $A^2$ = '+strtrim(stats.AD_stat[total_num_waves], 2), font_size=12)
        txt_AD_crit = text(0.755, 0.21, '   crit val = '+strtrim(stats.AD_crit[total_num_waves], 2), font_size=12)
        txt_LB_test = text(0.755, 0.16, 'Ljung-Box', font_size=12)
        txt_LB_stat = text(0.755, 0.13, '   $Q$ = '+strtrim(stats.LB_stat[total_num_waves], 2), font_size=12)
        txt_LB_crit = text(0.755, 0.10, '   $\chi^2$ = '+strtrim(stats.LB_chisqrd[total_num_waves], 2), font_size=12)

        ;##### GETTING USER INPUT #####
        print, ' ' ;empty line for cleaner printing
        print, 'Plotting thread '+strtrim(h+1, 2)+' out of '+strtrim(n_threads, 2)
        input_flag = '0' ;initializing as a string allows for more robust options
        keep_asking = 1
        WHILE keep_asking DO BEGIN
            READ, input_flag, PROMPT='Please enter fit quality flag (0=bad, 1=questionable, 2=good): '
            test_for_int = stregex(input_flag, '^[-+]?[0-9][0-9]*$')
            IF test_for_int NE -1 THEN BEGIN
                keep_asking = 0 ; quit loop if actually given an integer
            ENDIF ELSE BEGIN
                input_flag = strlowcase(input_flag)
                IF input_flag.startswith('q') THEN BEGIN
                    print, 'Exiting program early. Remaining user_qual_flags set to 0 (bad fit)'
                    input_flag = '0'
                    keep_asking = 0
                    end_early = 1
                ENDIF ELSE BEGIN
                    print, '   You must enter an integer value or type "q" to quit!'
                ENDELSE
            ENDELSE
        ENDWHILE
        user_flags[h] = fix(input_flag)
        slit_fft_peaks[h].user_qual_flag = input_flag
        IF end_early THEN h = n_threads ; should break out of the loop cleanly
    ENDFOR

    fig_wave.close
ENDIF

IF NOT KEYWORD_SET(use_temp_common_blocks) THEN BEGIN
    ; Updating the NUWT master COMMON block
    nuwt_fft_peaks[slitnum] = slit_fft_peaks
ENDIF ELSE BEGIN
    fft_peaks = slit_fft_peaks
ENDELSE

END
