;+
;NAME: CALC_NUWT_BULK_STATS
;
;PURPOSE:
;   Calcuates the bulk statistics for ALL good waves found by NUWT.
;
;INPUTS:
;   None directly. Loads in the common blocks containing the results from
;   'locate_things.pro', 'follow_thread.pro', and
;   'nuwt_apply_fft.pro'
;
;OPTIONAL INPUTS:
;   res - spatial resolution of data (in arcsec). Default is 1
;   cad - temporal cadence of the data (in s). Default is 1
;   km_per_arcsec - ratio of [km]/[arcsec]. Defaults to 725.27 which is the
;                   mean scale of the solar surface.
;   /final_units - [DEPRECIATED] if set, will assume wave parameters are already in the correct
;                  units of [km] for amplitude and [s^-1] for frequency
;   min_user_flag - minimum user quality flag to plot. These flags are set with
;                   'set_nuwt_qual_flags.pro'. Default is -1 (all waves)
;   max_user_flag - maximum user quality flag to plot (see above). Default is 1000
;   slitnum - virtual slit number to plot. By default, will show slit number "0".
;             Note: can also set to 'all' to plot a combined histogram of all slits.
;
;OUTPUTS:
;   bulk_stats_out - stucture with the calculated summary statistics for the
;                    bulk (ALL good waves) NUWT results
;
;HISTORY: Name---------Date---------Description
;         M Weberg 20 JUNE, 2017  Initial coding
;         M Weberg 27 JULY, 2017  Major code update
;                                 - renamed to procedure to 'CALC_NUWT_BULK_STATS'
;                                 - updated load methods to allow for merged data
;                                 - now returns the flattened parameters for ALL
;                                   non-zero nuwt waves BUT only calculates the
;                                   statistics based on the selected flag values.
;
;TO-DO / LIMITATIONS:
;   - update doc string
;-

PRO CALC_NUWT_BULK_STATS, res=res, cad=cad, km_per_arcsec=km_per_arcsec, final_units=final_units, $
                     amp_units=amp_units, freq_units=freq_units, period_units=period_units, $
                     slitnum=slitnum, use_temp_common_blocks=use_temp_common_blocks, $
                     min_user_flag=min_user_flag, max_user_flag=max_user_flag, $
                     min_auto_flag=min_auto_flag, max_auto_flag=max_auto_flag, $
                     vel_amp_mode=vel_amp_mode, $
                     bulk_stats_out=bulk_stats_out, params_out=params_out

COMPILE_OPT IDL2
;###############################################################################
;Setting default values
;###############################################################################
IF N_ELEMENTS(slitnum) EQ 0 THEN slitnum = 0
IF N_ELEMENTS(slitnum) GT 1 THEN slitnum = slitnum[0]
IF NOT KEYWORD_SET(min_user_flag) THEN min_user_flag = -1
IF NOT KEYWORD_SET(max_user_flag) THEN max_user_flag = 1000
IF min_user_flag GT max_user_flag THEN BEGIN
    print, 'WARNING: min_user_flag > max_user_flag!'
    print, '   max_user_flag has been reset to equal min_user_flag'
    max_user_flag = min_user_flag
ENDIF

IF NOT KEYWORD_SET(min_auto_flag) THEN min_auto_flag = -1
IF NOT KEYWORD_SET(max_auto_flag) THEN max_auto_flag = 1000
IF min_user_flag GT max_user_flag THEN BEGIN
    print, 'WARNING: min_auto_flag > max_auto_flag!'
    print, '   max_auto_flag has been reset to equal min_auto_flag'
    max_auto_flag = min_auto_flag
ENDIF

;@load_data#####################################################################
;LOAD FLATTENED NUWT WAVE PARAMETERS (ALSO INCLUDES UNIT CONVERSIONS)
;###############################################################################
;If not given, these parameter will be loaded from the NUWT results where needed
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 0.0
IF NOT KEYWORD_SET(res) THEN res = 0.0
IF NOT KEYWORD_SET(cad) THEN cad = 0.0
IF KEYWORD_SET(amp_units) THEN output_amp_units = amp_units ELSE output_amp_units = 'km'
;freq_units is NOT currently a user-changable parameter since freq is not plotted in this code
IF KEYWORD_SET(freq_units) THEN output_freq_units = freq_units ELSE output_freq_units = 'Hz'
IF KEYWORD_SET(period_units) THEN output_period_units = period_units ELSE output_period_units = 's'
IF KEYWORD_SET(use_temp_common_blocks) THEN use_temp_common_blocks = 1 ELSE use_temp_common_blocks = 0
IF KEYWORD_SET(use_refit) THEN use_refit = 1.0 ELSE use_refit = 0.0
IF KEYWORD_SET(vel_amp_mode) THEN vel_amp_mode = 1.0 ELSE vel_amp_mode = 0.0
peak_params = FLATTEN_NUWT_RESULTS(slitnum=slitnum, res=res, cad=cad, $
                                  km_per_arcsec=km_per_arcsec, $
                                  amp_units=output_amp_units, $
                                  freq_units=output_freq_units, $
                                  period_units=output_period_units, $
                                  min_user_flag=min_user_flag, $
                                  max_user_flag=max_user_flag, $
                                  min_auto_flag=min_auto_flag, $
                                  max_auto_flag=max_auto_flag, $
                                  vel_amp_mode=vel_amp_mode, $
                                  use_refit=use_refit, $
                                  use_temp_common_blocks=use_temp_common_blocks)

;update plot units (in case something defaulted somewhere)
output_amp_units = peak_params.amp_units
output_freq_units = peak_params.freq_units
output_period_units = peak_params.period_units
output_vel_amp_units = peak_params.vel_amp_units
output_duration_units = peak_params.duration_units

max_num_waves_per_thread = MAX(peak_params.wave_orders)
;@calc_stats####################################################################
;CALCULATE THE BULK STATS
;###############################################################################
bulk_stats = ORDEREDHASH('num_threads', 0, 'num_waves', 0, 'filtered_num_waves', 0,$
                         'wave_counts', intarr(max_num_waves_per_thread+1), $
                         'user_flag_range', [0.0, 0.0], $
                         'auto_flag_range', [0.0, 0.0], $
                         'wave_order_range', [0.0, 0.0], $
                         'amp_mean', 0.0, 'amp_stddev', 0.0, $
                         'amp_median', 0.0, 'amp_MAD', 0.0, $
                         'log_norm_amp_mean', 0.0, 'log_norm_amp_stddev', 0.0, $
                         'log_norm_amp_mode', 0.0, $
                         'amp_units', output_amp_units, $
                         'period_mean', 0.0, 'period_stddev', 0.0, $
                         'period_median', 0.0, 'period_MAD', 0.0,$
                         'log_norm_period_mean', 0.0, 'log_norm_period_stddev', 0.0,$
                         'log_norm_period_mode', 0.0, $
                         'period_units', output_period_units, $
                         'freq_mean', 0.0, 'freq_stddev', 0.0, $
                         'freq_median', 0.0, 'freq_MAD', 0.0,$
                         'log_norm_freq_mean', 0.0, 'log_norm_freq_stddev', 0.0,$
                         'log_norm_freq_mode', 0.0, $
                         'freq_units', output_freq_units, $
                         'vel_amp_mean', 0.0, 'vel_amp_stddev', 0.0, $
                         'vel_amp_median', 0.0, 'vel_amp_MAD', 0.0, $
                         'log_norm_vel_amp_mean', 0.0, 'log_norm_vel_amp_stddev', 0.0, $
                         'log_norm_vel_amp_mode', 0.0, $
                         'vel_amp_units', output_vel_amp_units, $
                         'duration_mean', 0.0, 'duration_stddev', 0.0, $
                         'duration_median', 0.0, 'duration_MAD', 0.0, $
                         'duration_units', output_duration_units, $
                         'log_norm_duration_mean', 0.0, 'log_norm_duration_stddev', 0.0, $
                         'log_norm_duration_mode', 0.0, /FOLD_CASE, /LOWERCASE)

filter_waves = WHERE((peak_params.user_flags GE min_user_flag) AND (peak_params.user_flags LE max_user_flag) AND $
                     (peak_params.auto_flags GE min_auto_flag) AND (peak_params.auto_flags LE max_auto_flag), /NULL)

total_num_waves = N_ELEMENTS(peak_params.amps)
filtered_num_waves = N_ELEMENTS(filter_waves)

bulk_stats['num_threads'] = peak_params.num_threads
bulk_stats['num_waves'] = total_num_waves
bulk_stats['filtered_num_waves'] = filtered_num_waves

IF filtered_num_waves GT 0 THEN BEGIN
    sub_wave_orders = peak_params.wave_orders[filter_waves]
    bulk_stats['wave_order_range'] = [min(sub_wave_orders), max(sub_wave_orders)]

    sub_user_flags = peak_params.user_flags[filter_waves]
    sub_auto_flags = peak_params.auto_flags[filter_waves]
    bulk_stats['user_flag_range'] = [min(sub_user_flags), max(sub_user_flags)]
    bulk_stats['auto_flag_range'] = [min(sub_auto_flags), max(sub_auto_flags)]

    FOR oo=1, max_num_waves_per_thread DO BEGIN
        ;find number of each order of waves
        loc_sub_order_waves = where(sub_wave_orders EQ oo)
        IF loc_sub_order_waves[0] NE -1 THEN BEGIN
            bulk_stats['wave_counts',oo] = n_elements(loc_sub_order_waves)
        ENDIF
    ENDFOR
    bulk_stats['wave_counts',0] = peak_params.num_threads - bulk_stats['wave_counts',1] ;num_threads without waves

    sub_peak_amps = peak_params.amps[filter_waves]
    sub_peak_periods = peak_params.periods[filter_waves]
    sub_peak_freqs = peak_params.freqs[filter_waves]
    sub_peak_vel_amps = peak_params.vel_amps[filter_waves]
    sub_peak_durations = peak_params.durations[filter_waves]

    ;Looping over a list of lists with variable:array pairs allows for more compact code
    FOREACH VAR, LIST(LIST('amp', sub_peak_amps), $
                      LIST('period', sub_peak_periods), $
                      LIST('freq', sub_peak_freqs), $
                      LIST('vel_amp', sub_peak_vel_amps), $
                      LIST('duration', sub_peak_durations)) DO BEGIN
        bulk_stats[VAR[0]+'_mean'] = MEAN(VAR[1])
        bulk_stats[VAR[0]+'_median'] = MEDIAN(VAR[1])
        bulk_stats[VAR[0]+'_stddev'] = STDDEV(VAR[1])
        bulk_stats[VAR[0]+'_MAD'] = MEANABSDEV(VAR[1], /median)
        ;ALT method using equations on Wikipedia
        N = FLOAT(N_ELEMENTS(VAR[1]))
        ; m = MEAN(VAR[1])
        ; sd = STDDEV(VAR[1])
        ; mu = alog(m/SQRT(1 + (sd/m)^2))
        ; sigma = SQRT(alog(1 + (sd/m)^2))
        mu = MEAN(alog(VAR[1]), /NAN)
        sigma = ((N-1)/N)*STDDEV(alog(VAR[1]), /NAN)
        bulk_stats['log_norm_'+VAR[0]+'_mean'] = EXP(mu + (sigma^2)/2.0)
        bulk_stats['log_norm_'+VAR[0]+'_stddev'] = SQRT((EXP(sigma^2) - 1)*EXP(2*mu+sigma^2))
        bulk_stats['log_norm_'+VAR[0]+'_mode'] =  EXP(mu - sigma^2)
    ENDFOREACH
ENDIF

;@output########################################################################
;OUTPUT STRUCTURES AND ARRAYS
;###############################################################################
bulk_stats_out = bulk_stats
params_out = peak_params

END
