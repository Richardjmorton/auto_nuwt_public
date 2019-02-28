;+
;NAME: PLOT_NUWT_WAVE_HIST
;
;PURPOSE:
;   Plots histograms of the peak FFT wave values from 'nuwt_apply_fft.pro' for all
;   significant waves. Returns one set of plots for each order (primary, secondary,
;   etc.) of waves as well as the combined distibution of all waves. Additionally,
;   a calculated log-normal normal distribution of each parameter will be overplotted
;   The vertical scale of the log-nromal distributions will be adjusted to match
;   the y-axis units of the histograms
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
;   /final_units - [DEPRECIATED] (will be removed in the near future)
;   min_user_flag - minimum user quality flag to plot. These flags are set with
;                   'set_nuwt_qual_flags.pro'. Default is -1 (all waves)
;   max_user_flag - maximum user quality flag to plot (see above). Default is 1000
;   VVV_range - range to plot for amp, period, or vel amp. Defaults are [0,2000] km,
;               [0,2000] s, & [0,50] km/s respectively
;   VVV_units - string with the plot units for amp or period (vel_amp units are
;               automatically determined from the other two). Default is [km]
;               for amp_units and [s] for period_units
;   VVV_binsize - size of histogram bins for amp, period, or vel amp in the selected
;                 units. Defaults are 25 km, 25 s, & 1 km/s for amp, period,
;                 and vel_amp respectively.
;   VVV_nbins - number of histogram bins for amp, period, or vel amp. By default,
;               VVV_nbins =  1 + (VVV_range[1] - VVV_range[0]) / VVV_binsize
;   /normalized - if set, will normalize the histograms by the number of waves
;   /hide_log_norm_dist - if set, will hide the lines showing the calculated
;                         log-normal distributions for each parameter.
;   /clean - if set, will suppress plotting symbols and text showing the mean
;            and median parameter values
;   ref_waves - stucture containing either input simulated wave parameters (from
;               'generate_kink_waves.pro') or the 'fft_peaks' output from a different
;               run of NUWT. Will be histogrammed and compared to the current results.
;               Used for testing and validation.
;   /simulated - if set, will assume the ref_waves structure contains input
;                simulated parameters (default)
;   /nuwt - if set, will assume the ref_waves structure contains output from
;           a seperate NUWT run. Note: if both /simulated and /nuwt are set, the
;           program will give precedence to the /simulated option.
;   slitnum - virtual slit number to plot. By default, will show slit number "0".
;             Note: can also set to 'all' to plot a combined histogram of all slits.
;   header - custom header text. Useful for identifying the source data.
;            Defaults to 'NUWT FFT stats'
;   run_tag - custom string that will be appended to the end of the header.
;             Normally used to keep track of differnt runs of the same dataset.
;             There is no default run_tag string.
;   save_folder - folder in which to save the plots. Defaults to the user's home folder.
;                 Will also append a '/' to the end if not included.
;   filename - name for the output PDF file. Default is "NUWT_fft_peaks"
;
;OUTPUTS:
;   PDF_file - multi-page PDF containing FFT stats and diagnostics for each
;              ordered group of waves (primary, secondary, etc.) The first page
;              will show the distributions of all waves combined.
;   bulk_stats_out - stucture with the calculated summary statistics for the
;                    selected NUWT waves
;   ref_stats_out - stucture with the calculated summary statistics for the
;                   reference data
;
;HISTORY: Name---------Date---------Description
;         M Weberg  8 Sept, 2016  Initial coding
;         M Weberg 14 SEPT, 2016  reconfigured to use NUWT master COMMON block
;         M Weberg 18  OCT, 2016   added output structure for summary stats
;         M Weberg 09  JAN, 2017   added the "all" option to "slitnum"
;         M Weberg ??  MAR, 2017  Now can load unit information directly from
;                                 NUWT common block structures
;         M Weberg ?? JUNE, 2017  corrected how log-normal values are calculated
;         M Weberg 25 JULY, 2017  Modified slitnum="all" mode to allow for merged
;                                 nuwt results from runs with different number of
;                                 timesteps.
;         M Weberg     FEB, 2017  Reworked how NUWT data is loaded and unit
;                                 calculations performed (cleaner code with more options)
;                                 Also changed the default to always plot the calculated
;                                 log-normal distributions (can be hidden with
;                                 /hide_log_norm_dist)
;
;TO-DO / LIMITATIONS:
;   - more options for the 'compare_waves' distibution that would allow comparisions
;     to other NUWT results and not just the simulated waves input for testing
;   - Automatic metadata handling to denote dataset and source of observations.
;     (currently managed via user defined text strings)
;-

PRO PLOT_NUWT_WAVE_HIST, res=res, cad=cad, km_per_arcsec=km_per_arcsec, final_units=final_units, $
                         min_user_flag=min_user_flag, max_user_flag=max_user_flag,$
                         min_auto_flag=min_auto_flag, max_auto_flag=max_auto_flag, $

                         amp_range=amp_range, amp_units=amp_units, $
                         amp_binsize=amp_binsize, amp_nbins=amp_nbins, $
                         period_range=period_range, period_units=period_units, $
                         period_binsize=period_binsize, period_nbins=period_nbins, $
                         vel_amp_range=vel_amp_range, $ ;vel_amp units are determined by the amp and period units
                         vel_amp_binsize=vel_amp_binsize, vel_amp_nbins=vel_amp_nbins, $
                         normalized=normalized, hide_log_norm_dist=hide_log_norm_dist, $
                         nonzero_only=nonzero_only, clean=clean, $

                         apply_rcf=apply_rcf, vel_amp_mode=vel_amp_mode, use_refit=use_refit, $
                         ref_waves=ref_waves, simulated=simulated, nuwt=nuwt, $
                         sim_res=sim_res, sim_cad=sim_cad, sim_km_per_arcsec=sim_km_per_arcsec, $
                         bulk_stats_out=bulk_stats_out, ref_stats_out=ref_stats_out, $
                         slitnum=slitnum, header=header, run_tag=run_tag, $
                         save_folder=save_folder, filename=filename, debug=debug

COMPILE_OPT IDL2
IF NOT KEYWORD_SET(debug) THEN BEGIN
    ;hides math underflow errors
    CurrentExceptVal = !Except
    !Except = 0
    void = CHECK_MATH()
ENDIF
;###############################################################################
;Setting default values
;###############################################################################
IF NOT KEYWORD_SET(min_user_flag) THEN min_user_flag = -1
IF NOT KEYWORD_SET(max_user_flag) THEN max_user_flag = 1000
IF min_user_flag GT max_user_flag THEN BEGIN
    print, 'WARNING: min_user_flag > max_user_flag!'
    print, '   max_user_flag has been reset to equal min_user_flag'
    max_user_flag = min_user_flag
ENDIF

IF NOT KEYWORD_SET(min_auto_flag) THEN min_auto_flag = 2
IF NOT KEYWORD_SET(max_auto_flag) THEN max_auto_flag = 1000
IF min_user_flag GT max_user_flag THEN BEGIN
    print, 'WARNING: min_auto_flag > max_auto_flag!'
    print, '   max_auto_flag has been reset to equal min_auto_flag'
    max_auto_flag = min_auto_flag
ENDIF

;Selecting which type (if any) of reference waves to show
IF KEYWORD_SET(simulated) AND KEYWORD_SET(nuwt) THEN nuwt = 0 ;simulated takes precedence

IF n_elements(slitnum) EQ 0 THEN slitnum = 0
IF N_ELEMENTS(slitnum) GT 1 THEN slitnum = slitnum[0]

IF NOT KEYWORD_SET(header) THEN header = 'NUWT Wave Histograms'
IF KEYWORD_SET(run_tag) THEN header = header + ' - '+run_tag
IF NOT KEYWORD_SET(filename) THEN filename = 'NUWT_wave_hist'
IF NOT KEYWORD_SET(save_folder) THEN CD, CURRENT=save_folder ;i.e. defaults to current folder
IF strlen(save_folder) GT 0 AND NOT save_folder.endswith('/') THEN save_folder = save_folder+'/'

wave_group_titles = ['All waves', 'Primary waves only', 'Secondary waves only', $
                     'Tertiary waves only', 'Quaternary waves only']
wave_group_colors = ['black', 'red', 'orange', 'purple', 'violet']

;@load_data#####################################################################
;LOAD FLATTENED NUWT WAVE PARAMETERS (ALSO INCLUDES UNIT CONVERSIONS)
;###############################################################################
;If not given, these parameter will be loaded from the NUWT results where needed
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 0.0
IF NOT KEYWORD_SET(res) THEN res = 0.0
IF NOT KEYWORD_SET(cad) THEN cad = 0.0
IF KEYWORD_SET(amp_units) THEN plot_amp_units = amp_units ELSE plot_amp_units = 'km'
;freq_units is NOT currently a user-changable parameter since freq is not plotted in this code
IF KEYWORD_SET(freq_units) THEN plot_freq_units = freq_units ELSE plot_freq_units = 'Hz'
IF KEYWORD_SET(period_units) THEN plot_period_units = period_units ELSE plot_period_units = 's'
IF KEYWORD_SET(use_refit) THEN use_refit = 1.0 ELSE use_refit = 0.0
IF KEYWORD_SET(vel_amp_mode) THEN vel_amp_mode = 1.0 ELSE vel_amp_mode = 0.0
all_params = FLATTEN_NUWT_RESULTS(slitnum=slitnum, res=res, cad=cad, $
                                  km_per_arcsec=km_per_arcsec, $
                                  amp_units=plot_amp_units, $
                                  freq_units=plot_freq_units, $
                                  period_units=plot_period_units, $
                                  min_user_flag=min_user_flag, $
                                  max_user_flag=max_user_flag, $
                                  min_auto_flag=min_auto_flag, $
                                  max_auto_flag=max_auto_flag, $
                                  vel_amp_mode=vel_amp_mode, $
                                  use_refit=use_refit)

;update plot units (in case something defaulted somewhere)
plot_amp_units = all_params.amp_units
plot_period_units = all_params.period_units

IF KEYWORD_SET(apply_rcf) THEN BEGIN
    rcf = SQRT(2) ;rotation correction factor
    all_params.amps = all_params.amps*rcf
    all_params.vel_amps = all_params.vel_amps*rcf
ENDIF

;@set_ranges####################################################################
;SETTING THE DEFAULT PLOT RANGES
;###############################################################################
;Note: these defaults were adjusted for waves in coronal polar plumes visable
;      in the 171 channel of SDO/AIA. Units conversions are included.
;      Different data sources and/or waves may require modified plot ranges

scale_amp_range = convert_nuwt_units('km', plot_amp_units, res=res, km_per_arcsec=km_per_arcsec)
IF NOT KEYWORD_SET(amp_range) THEN amp_range = [0.0, 2000.0]*scale_amp_range
IF NOT KEYWORD_SET(amp_binsize) THEN amp_binsize = 25.0*scale_amp_range
IF NOT KEYWORD_SET(amp_nbins) THEN amp_nbins = FIX(1 + (amp_range[1] - amp_range[0])/amp_binsize)

scale_period_range = convert_nuwt_units('s', plot_period_units, cad=cad)
IF NOT KEYWORD_SET(period_range) THEN period_range = [0.0, 2000.0]*scale_period_range
IF NOT KEYWORD_SET(period_binsize) THEN period_binsize = 25.0*scale_period_range
IF NOT KEYWORD_SET(period_nbins) THEN period_nbins = FIX(1 + (period_range[1] - period_range[0])/period_binsize)

IF NOT KEYWORD_SET(vel_amp_range) THEN vel_amp_range = [0.0, 50.0]*(scale_amp_range/scale_period_range)
IF NOT KEYWORD_SET(vel_amp_binsize) THEN vel_amp_binsize = 1.0*scale_amp_range/scale_period_range
IF NOT KEYWORD_SET(vel_amp_nbins) THEN vel_amp_nbins = FIX(1 + (vel_amp_range[1] - vel_amp_range[0])/vel_amp_binsize)

;@ref_params####################################################################
;EXTRACT AND CONVERT REFERENCE PARAMETERS
;###############################################################################
ref_count_outliers = {amp:0, period:0, vel_amp:0}
IF n_elements(ref_waves) GT 0 THEN BEGIN
    ;Extract reference wave parameters
    n_ref = n_elements(ref_waves.phase)
    ref_params = {amplitude:fltarr(n_ref), period:fltarr(n_ref), phase:fltarr(n_ref), $
                  velocity_amp:fltarr(n_ref), mean:fltarr(n_ref), slope:fltarr(n_ref), $
                  num_cyc:fltarr(n_ref), num_waves:fltarr(n_ref)}
    IF KEYWORD_SET(nuwt) THEN BEGIN
        ;Results from an alternate run of 'run_nuwt.pro'
        ;WARNING: THIS ASSUMES THE SAME UNITS WERE USED BY NUWT IN BOTH RUNS!!!
        ref_params.amplitude = ref_waves.peak_amplitude[0]*all_params.amp_dx
        ref_params.period = 1.0*all_params.period_dt/(ref_waves.peak_freq[0])
        ref_params.velocity_amp = 2*!PI*ref_params.amplitude/ref_params.period
        ; ref_params.num_cyc = (ref_coords.length*cad)/ref_params.period
        ref_params.phase = ref_waves.peak_phase[0]
        ;fixing periods of INF and NaN
        loc_nan = where(~finite(ref_params.period))
        ref_params.period[loc_nan] = 0.0
        ref_params.velocity_amp[loc_nan] = 0.0
        ref_params.num_waves = ref_waves.num_signif_peaks
    ENDIF ELSE BEGIN
        ;Simulated data from 'generate_kink_waves.pro'
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

        ref_params.amplitude = ref_waves.amp*sim_amp_scale
        ref_params.period = ref_waves.period*sim_period_scale
        ref_params.velocity_amp = 2*!PI*ref_params.amplitude/ref_params.period
        ; ref_params.num_cyc = (ref_coords.length*cad)/ref_params.period
        ref_params.phase = ref_waves.phase
        loc_shift_phase = where(ref_params.phase GT !PI, /NULL)
        IF n_elements(loc_shift_phase) GT 0 THEN BEGIN
            ;Adjust phase to be in the range of -pi to +pi (as returned by fft)
            ref_params.phase[loc_shift_phase] = ref_params.phase[loc_shift_phase] - 2*!PI
        ENDIF
        ref_params.num_waves[0,-1] = 1
    ENDELSE

    IF KEYWORD_SET(nonzero_only) THEN BEGIN
        hist_ref_ind = WHERE(ref_params.amplitude GT 0.0, /NULL)
    ENDIF ELSE BEGIN
        hist_ref_ind = INDGEN(n_ref)
    ENDELSE

    ref_amp_hist = histogram(ref_params.amplitude[hist_ref_ind], LOCATIONS=ref_amp_bins, $
                             nbins=amp_nbins, min=amp_range[0], max=amp_range[1])
    ref_period_hist = histogram(ref_params.period[hist_ref_ind], LOCATIONS=ref_period_bins, $
                                nbins=period_nbins, min=period_range[0], max=period_range[1])
    ref_vel_amp_hist = histogram(ref_params.velocity_amp[hist_ref_ind], LOCATIONS=ref_vel_amp_bins, $
                                 nbins=vel_amp_nbins, min=vel_amp_range[0], max=vel_amp_range[1])
ENDIF

IF n_elements(ref_waves) GT 0 THEN BEGIN
    ;Calculate summary statistics of ref distribution
    ref_summary_stats = ORDEREDHASH('amp_mean', 0.0, 'amp_stddev', 0.0,$
                                    'amp_median', 0.0, 'amp_MAD', 0.0,$
                                    'log_norm_amp_mean', 0.0, 'log_norm_amp_stddev', 0.0,$
                                    'log_norm_amp_mode', 0.0, $
                                    'period_mean', 0.0, 'period_stddev', 0.0,$
                                    'period_median', 0.0, 'period_MAD', 0.0,$
                                    'log_norm_period_mean', 0.0, 'log_norm_period_stddev', 0.0,$
                                    'log_norm_period_mode', 0.0, $
                                    'vel_amp_mean', 0.0, 'vel_amp_stddev', 0.0,$
                                    'vel_amp_median', 0.0, 'vel_amp_MAD', 0.0, $
                                    'log_norm_vel_amp_mean', 0.0, 'log_norm_vel_amp_stddev', 0.0,$
                                    'log_norm_vel_amp_mode', 0.0, /FOLD_CASE, /LOWERCASE)

    FOREACH VAR, LIST(LIST('amp', ref_params.amplitude[hist_ref_ind]), $
                      LIST('period', ref_params.period[hist_ref_ind]), $
                      LIST('vel_amp', ref_params.velocity_amp[hist_ref_ind])) DO BEGIN
        ref_summary_stats[VAR[0]+'_mean'] = MEAN(VAR[1])
        ref_summary_stats[VAR[0]+'_median'] = MEDIAN(VAR[1])
        ref_summary_stats[VAR[0]+'_stddev'] = STDDEV(VAR[1])
        ref_summary_stats[VAR[0]+'_MAD'] = MEANABSDEV(VAR[1], /median)
        ref_summary_stats['log_norm_'+VAR[0]+'_mean'] = EXP(MEAN(alog(VAR[1])))
        ref_summary_stats['log_norm_'+VAR[0]+'_stddev'] = EXP(STDDEV(alog(VAR[1])))
        ref_summary_stats['log_norm_'+VAR[0]+'_mode'] =  EXP(MEAN(alog(VAR[1])) - STDDEV(alog(VAR[1]))^2)
    ENDFOREACH

    ; Counting datapoints outside the histogram range
    ref_count_outliers = {amp:0, period:0, vel_amp:0}
    ref_count_outliers.amp = n_elements(where(ref_params.amplitude[hist_ref_ind] GT amp_range[1], /NULL))
    ref_count_outliers.period = n_elements(where(ref_params.period[hist_ref_ind] GT period_range[1], /NULL))
    ref_count_outliers.vel_amp = n_elements(where(ref_params.velocity_amp[hist_ref_ind] GT vel_amp_range[1], /NULL))

    ref_stats_out = ref_summary_stats
ENDIF

;@plot##########################################################################
;MAKE THE PLOTS
;###############################################################################
FOR page_num=1, 5 DO BEGIN
    ;##### INITIALIZE THE FIGURE AND FILTER WAVES
    fig_fft_peaks = window(name='fig_fft_peaks', dimensions=[600, 720], /buffer)

    IF page_num EQ 1 THEN BEGIN ;All waves
        wave_order_loc = where(all_params.wave_orders GT 0, /NULL)
    ENDIF ELSE BEGIN ;Groups of ONLY a certain wave order (primary, secondary, etc.)
        wave_order_loc = where(all_params.wave_orders EQ page_num-1, /NULL)
    ENDELSE

    total_num_waves = n_elements(wave_order_loc)

    ;##### CALCULATE HISTOGRAMS AND VARIOUS SIMPLE STATISTICS
    summary_stats = ORDEREDHASH('amp_mean', 0.0, 'amp_stddev', 0.0,$
                                'amp_median', 0.0, 'amp_MAD', 0.0,$
                                'log_norm_amp_mean', 0.0, 'log_norm_amp_stddev', 0.0,$
                                'log_norm_amp_mode', 0.0, $
                                'period_mean', 0.0, 'period_stddev', 0.0,$
                                'period_median', 0.0, 'period_MAD', 0.0,$
                                'log_norm_period_mean', 0.0, 'log_norm_period_stddev', 0.0,$
                                'log_norm_period_mode', 0.0, $
                                'vel_amp_mean', 0.0, 'vel_amp_stddev', 0.0,$
                                'vel_amp_median', 0.0, 'vel_amp_MAD', 0.0, $
                                'log_norm_vel_amp_mean', 0.0, 'log_norm_vel_amp_stddev', 0.0,$
                                'log_norm_vel_amp_mode', 0.0, /FOLD_CASE, /LOWERCASE)

    IF total_num_waves GT 0 THEN BEGIN
        filtered_amps = all_params.amps[wave_order_loc]
        filtered_periods = all_params.periods[wave_order_loc]
        filtered_vel_amps = all_params.vel_amps[wave_order_loc]

        IF KEYWORD_SET(nonzero_only) THEN BEGIN
            hist_nuwt_ind = WHERE(filtered_amps GT 0.0, /NULL)
        ENDIF ELSE BEGIN
            hist_nuwt_ind = INDGEN(total_num_waves)
        ENDELSE

        plot_amp_hist = histogram(filtered_amps[hist_nuwt_ind], LOCATIONS=plot_amp_bins, $
                                  nbins=amp_nbins, min=amp_range[0], max=amp_range[1])
        plot_period_hist = histogram(filtered_periods[hist_nuwt_ind], LOCATIONS=plot_period_bins, $
                                     nbins=period_nbins, min=period_range[0], max=period_range[1])
        plot_vel_amp_hist = histogram(filtered_vel_amps[hist_nuwt_ind], LOCATIONS=plot_vel_amp_bins, $
                                      nbins=vel_amp_nbins, min=vel_amp_range[0], max=vel_amp_range[1])

        ;Looping over a list of lists with variable:array pairs allows for more compact code
        color_summary_stats = 'black'
        FOREACH VAR, LIST(LIST('amp', filtered_amps[hist_nuwt_ind]), $
                          LIST('period', filtered_periods[hist_nuwt_ind]), $
                          LIST('vel_amp', filtered_vel_amps[hist_nuwt_ind])) DO BEGIN
            summary_stats[VAR[0]+'_mean'] = MEAN(VAR[1], /NAN)
            summary_stats[VAR[0]+'_median'] = MEDIAN(VAR[1])
            summary_stats[VAR[0]+'_stddev'] = STDDEV(VAR[1], /NAN)
            summary_stats[VAR[0]+'_MAD'] = MEANABSDEV(VAR[1], /median, /NAN)
            N = float(N_ELEMENTS(VAR[1]))
            mu = MEAN(alog(VAR[1]), /NAN)
            sigma = ((N-1)/N)*STDDEV(alog(VAR[1]), /NAN)
            summary_stats['log_norm_'+VAR[0]+'_mean'] = EXP(mu + (sigma^2)/2.0)
            summary_stats['log_norm_'+VAR[0]+'_stddev'] = SQRT((EXP(sigma^2) - 1)*EXP(2*mu+sigma^2))
            summary_stats['log_norm_'+VAR[0]+'_mode'] =  EXP(mu - sigma^2)
        ENDFOREACH
    ENDIF ELSE BEGIN
        ;If there if is no waves of the given order, create dummy arrays
        color_summary_stats = 'dark grey'
        filtered_amps = fltarr(10)
        filtered_periods = fltarr(10)
        filtered_vel_amps = fltarr(10)
        plot_amp_hist = fltarr(10)
        plot_amp_bins = findgen(10)
        plot_period_hist = fltarr(10)
        plot_period_bins = findgen(10)
        plot_vel_amp_hist = fltarr(10)
        plot_vel_amp_bins = findgen(10)/10.0
    ENDELSE

    ; Counting datapoints outside the histogram range
    test_count_outliers = {amp:0, period:0, vel_amp:0}
    test_count_outliers.amp = n_elements(where(filtered_amps GT amp_range[1], /NULL))
    test_count_outliers.period = n_elements(where(filtered_periods GT period_range[1], /NULL))
    test_count_outliers.vel_amp = n_elements(where(filtered_vel_amps GT vel_amp_range[1], /NULL))

    IF page_num EQ 1 THEN bulk_stats_out = summary_stats

    ;Normalize the histogram by the total number of waves
    IF KEYWORD_SET(normalized) AND total_num_waves GT 0 THEN BEGIN
        plot_amp_hist = float(plot_amp_hist)/float(total_num_waves)
        plot_period_hist = float(plot_period_hist)/float(total_num_waves)
        plot_vel_amp_hist = float(plot_vel_amp_hist)/float(total_num_waves)
        y_axis_label = 'Fraction of waves'
    ENDIF ELSE y_axis_label = 'Number of waves'

    ;##### TOP SUBPLOT #####
    plt_amp = barplot(plot_amp_bins, plot_amp_hist, name='Amplitude', $
                      color='black', fill_color='light grey', linestyle='-', $
                      position=[0.10, 0.7, 0.72, 0.92], /current, histogram=1, xrange=amp_range)
    plt_amp.xtitle = 'Displacement Amplitude [$'+plot_amp_units+'$]'
    plt_amp.ytitle = y_axis_label

    ;##### MIDDLE SUBPLOT #####
    plt_period = barplot(plot_period_bins, plot_period_hist , name='Period', $
                         color='black', fill_color='light grey', linestyle='-', $
                         position=[0.10, 0.39, 0.72, 0.61], /current, histogram=1, xrange=period_range)
    plt_period.xtitle = 'Period [$'+plot_period_units+'$]'
    plt_period.ytitle = y_axis_label

    ;##### BOTTOM SUBPLOT #####
    plt_vel_amp = barplot(plot_vel_amp_bins, plot_vel_amp_hist, name='Velocity Amplitude', $
                          color='black', fill_color='light grey', linestyle='-', $
                          position=[0.10, 0.08, 0.72, 0.30], /current, histogram=1, xrange=vel_amp_range)
    plt_vel_amp.xtitle = 'Velocity Amplitude [$'+plot_amp_units+' '+plot_period_units+'^{-1}$]'
    plt_vel_amp.ytitle = y_axis_label

    IF n_elements(ref_waves) GT 0 AND page_num LT 3 THEN BEGIN
        ; Plotting Reference distributions
        IF KEYWORD_SET(normalized) AND n_ref GT 0 THEN BEGIN
            ref_amp_hist = float(ref_amp_hist)/float(n_ref)
            ref_period_hist = float(ref_period_hist)/float(n_ref)
            ref_vel_amp_hist = float(ref_vel_amp_hist)/float(n_ref)
        ENDIF
        plt_ref_amp = plot(ref_amp_bins, ref_amp_hist, /histogram, $
                            color='dark red', linestyle='-', overplot=plt_amp)
        plt_ref_period = plot(ref_period_bins, ref_period_hist, /histogram, $
                               color='dark red', linestyle='-', overplot=plt_period)
        plt_ref_vel_amp = plot(ref_vel_amp_bins, ref_vel_amp_hist, /histogram, $
                                color='dark red', linestyle='-', overplot=plt_vel_amp)
    ENDIF

    IF NOT KEYWORD_SET(hide_log_norm_dist) AND total_num_waves GT 0 THEN BEGIN
        ; Overplotting sample log-normal distributions
        ; amp_kde = CALC_KDE(filtered_amps[hist_nuwt_ind], range=amp_range, num_vals=amp_nbins*5, xvals_out=amp_kde_xvals)
        ; scaled_amp_kde = amp_kde*TOTAL(plot_amp_hist*amp_binsize)
        ; plt_amp_kde = plot(amp_kde_xvals, scaled_amp_kde, $
        ;                    color='blue', linestyle='-', overplot=plt_amp)

        amp_dist_x_vals = (amp_range[1] - amp_range[0])*findgen(amp_nbins*2)/(amp_nbins*2) + amp_range[0]
        ; amp_dist_y_vals = LOG_NORMAL_PDF(amp_dist_x_vals, arith_mean=summary_stats['amp_mean'], $
        ;                                  arith_stddev=summary_stats['amp_stddev'])
        amp_dist_y_vals = LOG_NORMAL_PDF(amp_dist_x_vals, arith_mean=summary_stats['log_norm_amp_mean'], $
                                         arith_stddev=summary_stats['log_norm_amp_stddev'])
        scaled_amp_dist_y_vals = amp_dist_y_vals*TOTAL(plot_amp_hist*amp_binsize)
        plt_log_norm_amp = plot(amp_dist_x_vals, scaled_amp_dist_y_vals, $
                            color='green', linestyle='-', overplot=plt_amp)

        period_dist_x_vals = (period_range[1] - period_range[0])*findgen(period_nbins*2)/(period_nbins*2) + period_range[0]
        ; period_dist_y_vals = LOG_NORMAL_PDF(period_dist_x_vals, arith_mean=summary_stats['period_mean'], $
        ;                                     arith_stddev=summary_stats['period_stddev'])
        period_dist_y_vals = LOG_NORMAL_PDF(period_dist_x_vals, arith_mean=summary_stats['log_norm_period_mean'], $
                                            arith_stddev=summary_stats['log_norm_period_stddev'])
        scaled_period_dist_y_vals = period_dist_y_vals*TOTAL(plot_period_hist*period_binsize)
        plt_log_norm_period = plot(period_dist_x_vals, scaled_period_dist_y_vals, $
                               color='green', linestyle='-', overplot=plt_period)

        vel_amp_dist_x_vals = (vel_amp_range[1] - vel_amp_range[0])*findgen(vel_amp_nbins*2)/(vel_amp_nbins*2) + vel_amp_range[0]
        ; vel_amp_dist_y_vals = LOG_NORMAL_PDF(vel_amp_dist_x_vals, arith_mean=summary_stats['vel_amp_mean'], $
        ;                                      arith_stddev=summary_stats['vel_amp_stddev'])
        vel_amp_dist_y_vals = LOG_NORMAL_PDF(vel_amp_dist_x_vals, arith_mean=summary_stats['log_norm_vel_amp_mean'], $
                                             arith_stddev=summary_stats['log_norm_vel_amp_stddev'])
        scaled_vel_amp_dist_y_vals = vel_amp_dist_y_vals*TOTAL(plot_vel_amp_hist*vel_amp_binsize)
        plt_log_norm_period = plot(vel_amp_dist_x_vals, scaled_vel_amp_dist_y_vals, $
                               color='green', linestyle='-', overplot=plt_vel_amp)
    ENDIF

    ; Plotting arrows indicating outliers
    upper_outlier_xloc = 0.69
    FOREACH VAR, LIST(LIST(test_count_outliers.amp, 0.895, 'black'), $
                      LIST(ref_count_outliers.amp, 0.875, 'dark red'), $
                      LIST(test_count_outliers.period, 0.585, 'black'), $
                      LIST(ref_count_outliers.period, 0.565, 'dark red'), $
                      LIST(test_count_outliers.vel_amp, 0.275, 'black'), $
                      LIST(ref_count_outliers.vel_amp, 0.255, 'dark red')) DO BEGIN
        ; List of lists of [test_var, x_coord, y_coord, sym_color]
        IF VAR[0] GT 0 THEN BEGIN
            sym_outliers = symbol(upper_outlier_xloc, VAR[1], sym_text='$\rightarrow$', sym_color=VAR[2], $
                                  label_color=VAR[2], label_position='left', label_shift=[0.0, -0.005], $
                                  label_string=strtrim(VAR[0], 2), label_font_size=10)
        ENDIF
    ENDFOREACH

    IF NOT KEYWORD_SET(clean) THEN BEGIN
        IF total_num_waves GT 0 THEN BEGIN
            ;add symbols along the edge showing the mean, median, log-norm mean, and log-norm mode
            rel_sym_ylocs = [-0.04, 1.04]

            FOREACH VAR, LIST(LIST('amp', amp_range[1], plt_amp), $
                              LIST('period', period_range[1], plt_period), $
                              LIST('vel_amp', vel_amp_range[1], plt_vel_amp)) DO BEGIN

                rel_mean_xloc = summary_stats[VAR[0]+'_mean']/VAR[1]
                arith_mean_sym = symbol(rel_mean_xloc, -0.035, $
                                        'triangle', sym_color='blue', $
                                        sym_size=0.9, /relative, target=VAR[2])
                arith_mean_sym = symbol(rel_mean_xloc, 1.035, $
                                        'triangle_down', sym_color='blue', $
                                        sym_size=0.9, /relative, target=VAR[2])

                rel_median_xloc = summary_stats[VAR[0]+'_median']/VAR[1]
                arith_median_sym = symbol([rel_median_xloc, rel_median_xloc],rel_sym_ylocs, $
                                          'vline', sym_color='blue', $
                                          sym_size=0.9, /relative, target=VAR[2])

                rel_log_norm_mean_xloc = summary_stats['log_norm_'+VAR[0]+'_mean']/VAR[1]
                log_norm_mean_sym = symbol(rel_log_norm_mean_xloc, -0.035, $
                                           'triangle', sym_color='green', $
                                           sym_size=0.9, /relative, target=VAR[2])
                log_norm_mean_sym = symbol(rel_log_norm_mean_xloc, 1.035, $
                                           'triangle_down', sym_color='green', $
                                           sym_size=0.9, /relative, target=VAR[2])

                rel_log_norm_mode_xloc = summary_stats['log_norm_'+VAR[0]+'_mode']/VAR[1]
                log_norm_mode_sym = symbol([rel_log_norm_mode_xloc, rel_log_norm_mode_xloc], rel_sym_ylocs, $
                                           'vline', sym_color='green', $
                                           sym_size=0.9, /relative, target=VAR[2])
            ENDFOREACH
        ENDIF

        ;##### PRINTED STATISTICS #####
        fake_legend_xloc = 0.74
        txt_stat_xloc = 0.76
        FOREACH VAR, LIST(LIST('amp', 0.90), $
                          LIST('period', 0.59), $
                          LIST('vel_amp', 0.28)) DO BEGIN

            fake_legend_arith_mean = symbol(fake_legend_xloc, VAR[1], 'triangle_right', $
                                            sym_color='blue', sym_size=0.9, /normal)
            fake_legend_arith_median = symbol(fake_legend_xloc, VAR[1]-0.05, 'vline', $
                                              sym_color='blue', sym_size=0.9, /normal)
            fake_legend_log_norm_mean = symbol(fake_legend_xloc, VAR[1]-0.13, 'triangle_right', $
                                               sym_color='green', sym_size=0.9, /normal)
            fake_legend_log_norm_mode = symbol(fake_legend_xloc, VAR[1]-0.18, 'vline', $
                                               sym_color='green', sym_size=0.9, /normal)

            tag_arithmetic = text(fake_legend_xloc, VAR[1]+0.025, $
                                  '--- Arithmetic ---', font_size=11, color='blue')
            txt_amp_mean = text(txt_stat_xloc, VAR[1]-0.02, 'Mean $\pm$ Std. Dev.$\n$'+$
                                strtrim(summary_stats[VAR[0]+'_mean'], 2)+' $\pm$ '+$
                                strtrim(summary_stats[VAR[0]+'_stddev'], 2), $
                                font_size=11, color=color_summary_stats)
            txt_amp_median = text(txt_stat_xloc, VAR[1]-0.07, 'Median $\pm$ MAD$\n$'+$
                                  strtrim(summary_stats[VAR[0]+'_median'], 2)+' $\pm$ '+$
                                  strtrim(summary_stats[VAR[0]+'_mad'], 2), $
                                  font_size=11, color=color_summary_stats)

            tag_log_normal = text(fake_legend_xloc, VAR[1]-0.105, $
                                  '--- Log-Normal ---', font_size=11, color='green')
            txt_ln_amp_mean = text(txt_stat_xloc, VAR[1]-0.15, 'Mean $\pm$ Std. Dev.$\n$'+$
                                   strtrim(summary_stats['log_norm_'+VAR[0]+'_mean'], 2)+' $\pm$ '+$
                                   strtrim(summary_stats['log_norm_'+VAR[0]+'_stddev'], 2), $
                                   font_size=11, color=color_summary_stats)
            txt_ln_amp_mode = text(txt_stat_xloc, VAR[1]-0.20, 'Mode$\n$'+$
                                     strtrim(summary_stats['log_norm_'+VAR[0]+'_mode'], 2), $
                                     font_size=11, color=color_summary_stats)
        ENDFOREACH
    ENDIF

    ;Header text
    txt_slit_ID = text(0.01, 0.97, 'Slit # '+strtrim(slitnum, 2), font_size=14, font_style=1)
    txt_header = text(0.20, 0.97, header, font_size=14, font_style=1)
    txt_num_threads = text(0.01, 0.94, strtrim(all_params.num_threads, 2)+' Threads', font_size=12)
    txt_wave_group = text(0.25, 0.94, wave_group_titles[page_num-1], font_size=12, $
                          font_style=1, color=wave_group_colors[page_num-1])
    txt_num_waves = text(0.59, 0.94, strtrim(total_num_waves, 2)+' Waves', font_size=12)

    txt_timestamp = text(0.76, 0.01, 'Date & Time Created$\n$'+systime(), font_size=9, color='dark grey')

    IF page_num LT 5 THEN BEGIN
        fig_fft_peaks.save, save_folder+filename+'.pdf', page_size=[6.0, 7.2], /append
    ENDIF ELSE BEGIN
        fig_fft_peaks.save, save_folder+filename+'.pdf', page_size=[6.0, 7.2], /append, /close
    ENDELSE
ENDFOR

IF NOT KEYWORD_SET(debug) THEN BEGIN
    ;manually check for other types of math errors
    floating_point_underflow = 32
    status = Check_Math() ; Get status and reset accumulated math error register.
    IF(status AND NOT floating_point_underflow) NE 0 THEN BEGIN
        print, 'IDL Check_Math() error: ' + StrTrim(status, 2)
    ENDIF

   !Except = CurrentExceptVal
ENDIF

print, 'Finished plotting wave histograms!'
print, 'Save Folder: ', save_folder
print, 'Filename: ', filename+'.pdf'
END
