;+
;NAME: NUWT_APPLY_FFT
;
;PURPOSE:
;   Performs a Fast Fourier Transform (FFT) on the positions of each thread
;   found by 'follow_thread.pro' and then calculates the Power Spectral Density
;   (PSD). Results are saved to a list of structures, one for each thread. Also
;   computes and saves a number of other variables and reference values for
;   later plotting. Based on the program 'pwr_frm_fft.pro'.
;
;INPUTS:
;   None directly. Loads in the common block data from 'follow_thread.pro'.
;
;OPTIONAL INPUTS:
;   res - spatial resolution of data in [arcsec]. Defaults to a value of 1 [pixel]
;   cad - temporal cadence of the data in [s]. Defaults to a value of 1 [timestep]
;   km_per_arcsec - ratio of [km]/[arcsec]. Defaults to 725.27 which is the
;                   mean scale of the solar surface.
;   detrend - if keyword set will remove linear trend from each thread. If not
;            set, the mean value will be removed instead (needed for a good FFT)
;   /pad_fft - if set, will pad the array with zeros to improve the precision
;              of peak frequencies returned by the fft. Note, this does NOT
;              actually increase the resolution and can results in extra spectral
;              leakage. Use with care.
;   pad_length - number of zeros to pad to the end of the array when /pad_fft is set.
;                If /fill_pad is ALSO set, will instead pad the array until the
;                total length equals pad_length (or just pad with 1 zero if the
;                array is already longer than pad_length). Default is 1000.
;   /fill_pad - if set, will fill the array with extra zeros until the length equals
;               "pad_length" instead of appending a set number zeros (default).
;               Note: if the pad_length is set too low, using /fill_pad may
;               introducing "banding" in the possible frequency values returned.
;   /bootstrap - if set, will perform bootstrapping to resample the data and obtaining
;                error estimates on the output values
;   num_bootstrap - number of bootstrap resamples to use. Default is 1000
;   include_nyquist - [EXPERIMENTAL!] For even number of data points, will
;                     include the results at the nyquist frequency.
;   min_freq_cutoff - minimum allowable frequency for a significant peak. needed
;                     to suppress peaks of noise that incorrectly get selected
;                     when using a large amount of FFT padding. Default is 1e-4.
;                     Note: this is only used when "min_cyc_cutoff" is set to 0.0
;   min_cyc_cutoff - minimum number of cycles to allow for a valid wave result.
;                    Will not return any waves in which there is not enough actually
;                    data to observed the given number of cycles. Default is 0.75
;                    Note: this overides the "min_freq_cutoff" variable above.
;   window_func - select which windowing (also known as tapering or apodization)
;                 function to apply to the data. Choose from 'split_cosine_bell'
;                 (default) or 'hann'. See FFT tutorials / textbooks for more
;                 information.
;   window_param - set the parameter, p, to use in the generation of the window
;                  function. For the split cosine bell window, p indicates the
;                  total fraction of data that will be tappered (half on each end
;                  of the data series). At p = 1.0, the split cosine bell window
;                  is identical to the hann window. Default is 0.4 (i.e. 40%)
;   signif_levels - sets significance levels for FFT results (default=0.95)
;   /vaughn - if set, will use the significance test from Vaughan 2005. Otherwise
;             will default to the method of Torrence & Compo 1998.
;   /recover_power - [EXPERIMENTAL] attempts to recover leaked PSD by summing adjacent
;                    bins and performing a weighted mean to get the frequency
;   /adjacent_peaks - [EXPERIMENTAL]  Allows adjacent FFT bins to be classified as a
;                     seperate waves. WARNING! May not work well with FFT padding
;   /iterate_padding - [EXPERIMENTAL] Mode used for running the FFT a twice, the
;                      first time normally, and the second time attemping to pad the
;                      FFT spectrum such that the primary peak lands in the center
;                      of an FFT bin (hopefully minimizing spectral leakage)
;   /vel_amp_mode - [EXPERIMENTAL] alternate mode that attempts to directly calculate
;                   the velocity ampitude by first determining the "velocity time series"
;                   using a central differnce method and then running the FFT on it.
;                   The desire if for a more accurate vel. amp. value that does not
;                   depend as much on our frequencey resolution.
;
;OUTPUTS (in the 'fft_results_dat' COMMON block):
;   fft_spec - LIST of structures with the FFT results for each thread. Since
;                 the length of each thread is different, the FFT PSD spectra will
;                 have different array sizes (hence the list of structures)
;   fft_peaks - ARRAY of structures containing parameters for the four most
;               significant waves, GOF test statistics, and various other bits
;               of information used for tracking and plotting
;
;   The format of each structure is as follows:
;
;   spec_th_fft = {power:fltarr(nf), $
;                    amplitude:fltarr(nf), $
;                    freq:fltarr(nf), $
;                    phase:fltarr(nf), $
;                    err_power:fltarr(nf), $
;                    err_amplitude:fltarr(nf), $
;                    err_phase:fltarr(nf), $
;                    trend:fltarr(len_thread), $
;                    apod_window:fltarr(len_thread), $
;                    bin_flags:intarr(nf), $
;                    signif_vals:fltarr(nf), $
;                    power_to_pxls:convert_pow_to_pxls, $
;                    power_units:pow_output_units, $
;                    amp_to_pxls:convert_amp_to_pxls, $
;                    amp_units:amp_output_units, $
;                    freq_to_timesteps:convert_freq_to_timesteps, $
;                    freq_units:freq_output_units}
;
;   peaks_th_fft =  {peak_power:[0.0, 0.0, 0.0, 0.0], $
;                    peak_amplitude:[0.0, 0.0, 0.0, 0.0], $
;                    peak_freq:[0.0, 0.0, 0.0, 0.0], $
;                    peak_vel_amp:[0.0, 0.0, 0.0, 0.0], $
;                    peak_phase:[0.0, 0.0, 0.0, 0.0], $
;                    err_peak_power:[0.0, 0.0, 0.0, 0.0], $
;                    err_peak_amplitude:[0.0, 0.0, 0.0, 0.0], $
;                    err_peak_phase:[0.0, 0.0, 0.0, 0.0], $
;                    peak_bin:[-1, -1, -1, -1], $
;                    KS_stat:[-1.0, -1.0, -1.0, -1.0, -1.0], $
;                    KS_prob:[-1.0, -1.0, -1.0, -1.0, -1.0], $
;                    AD_stat:[-1.0, -1.0, -1.0, -1.0, -1.0], $
;                    AD_crit:[-1.0, -1.0, -1.0, -1.0, -1.0], $
;                    LB_stat:[-1.0, -1.0, -1.0, -1.0, -1.0], $
;                    LB_chisqrd:[-1.0, -1.0, -1.0, -1.0, -1.0], $
;                    num_signif_peaks:0, $
;                    num_saved_waves:0, $
;                    fft_length:-1, $
;                    window_func:'not_applicable', $
;                    window_param:-1.0, $
;                    enbw:-1.0, $
;                    cpg:-1.0, $
;                    signif_level:0.0, $
;                    signif_test:'not_applicable', $
;                    adjacent_peaks:allow_adj_pks, $
;                    user_qual_flag:-1, $
;                    auto_qual_flag:1, $
;                    power_to_pxls:convert_pow_to_pxls, $
;                    power_units:pow_output_units, $
;                    amp_to_pxls:convert_amp_to_pxls, $
;                    amp_units:amp_output_units, $
;                    freq_to_timesteps:convert_freq_to_timesteps, $
;                    freq_units:freq_output_units}

;   Quick key for values in 'bin_flags' :
;       -2 : invalid thread (too litle real data)
;       -1 : empty or invalid frequency bin (not used in fft of particular thread)
;        0 : fft result below selected significance level
;        1 : significant value in fft result
;        2 : local maximum in nearby significant fft results
;
;EXTERNAL CALLS - locate_things_min.pro, locate_things.pro, follow_thread.pro, cospd.pro
;
;HISTORY: Name---------Date---------Description
;         R Morton   ???, 2014 - Initial coding of 'pwr_frm_fft.pro'
;         R Morton   MAR, 2015 - Updated to include structures. [note: this older
;                                format had much larger, mostly empty, arrays].
;         M Weberg  JUNE, 2016 - Reworked the code. Seperated the FFT code (this
;                                program) from the rest of the NUWT code, cleaned
;                                up the logic a little bit, and added the new
;                                output format (with extra information for plotting)
;         M Weberg  SEPT, 2016 - Merged with the experimental bootstrapping fork of
;                                of this code. Also split the output to 'fft_results'
;                                & 'fft_peaks', added calculation of GOF test stats,
;                                and coded some experimental options (/adjacent_peaks
;                                and /recover_power)
;         M Weberg   OCT, 2016 - Reworked the FFT padding variable names and added
;                                the /fill_pad option as an alternative to just
;                                appending a set number of zeros. Also added min
;                                freq and cycle cutoffs to avoid bad padding results
;                                caused by peaks in the noise.
;         M Weberg   FEB, 2017 - Added an experimental "vel_amp_mode" to try directly
;                                calculating the velocity amplitude from a velocity
;                                time-series. Inital tests give mixed results.
;         M Weberg   MAR, 2017 - Moderate update.
;                              - Added in automatic conversion to physical units,
;                              - Changed default power & amp. scaling to those given
;                                in the Heinzel FFT notes (peak values should be the same;
;                                only the power spectrum units are different).
;                              - Renamed the "ctrend" option to "detrend"
;         M Weberg  SEPT, 2017 - Renamed "fft_results" structures to "fft_spec".
;                              - Moved (or just copied) some variables from
;                                "fft_peaks" to "fft_spec" (keep things bundled better)
;
;RESTRICTIONS: LOTS AT THE MOMENT. ROUTINE IS VERY MUCH DEPENDENT UPON USERS
;
;TO DO / RESTRICTIONS:
;   - Currently only saves the wave parameters for the 4 largest, significant
;     peaks. 5 or more peaks should not be all that common but we may be wrong...
;     The 'num_signif_peaks' variable will still count correctly so, it may not
;     be too hard to find the 'missing' wave parameters when needed (if ever).
;   - Add calculation for the power averaged over all threads. Use the following?
;       "Average is calculated by creating frequency bins, taking the ln of power
;        in frequency bins, plot the ln power PDF and fitting a Gaussian. Centroid
;        of Gaussian (mu) gives median value of distribution, i.e. exp(mu)=median."
;   - Add method to transfer results to a structure of arrays? This would allow
;     cleaner & faster indexing later but will require the output array to be
;     scaled for the largest results (with fill values for the smaller arrays).
;-

;bootstrap resampling of the data for obtaining error estimations
FUNCTION BOOTSTRAP_1D_DATA, data_in, errors_in, num_resamples=num_resamples
    ;Notes: (1) this is designed only for 1D input arrays
    ;       (2) if you are zero padding (for FFT), it is recommended to pad
    ;           AFTER performing bootstrapping on the real data (saves time)
    COMMON bootstrap_seed_num, bootseed

    IF NOT keyword_set(num_resamples) THEN num_resamples = 1000

    len_data = n_elements(data_in)
    output_array = fltarr(len_data, num_resamples)
    mean_array = rebin(data_in, len_data, num_resamples)
    sigma_array = rebin(errors_in, len_data, num_resamples)

    ;generate normally distributed numbers for each data point
    FOR i=0L, (len_data-1) DO BEGIN
        output_array[i,*] = RANDOMN(bootseed, num_resamples)
    ENDFOR

    ;adjust the mean and width to match the data points and their errors
    ;we should then have noramally (i.e. gaussian) distributed numbers
    ;centered on each point and sampled from a dist. with sigma = error
    output_array = temporary(output_array)*sigma_array + mean_array

RETURN, output_array
END

;temporal apodization using a Split Cosine Bell (or Hann for p = 1.0) window
FUNCTION GET_APOD_WINDOW, tn, cpg=cpg, parameter=p, win_function=win_function
    w_t = fltarr(tn) + 1

    IF n_elements(p) EQ 0 THEN p = 0.4

    IF win_function EQ 'split_cosine_bell' OR win_function EQ 'hann' THEN BEGIN
        num_taper = tn * (p / 2.0)
        ;w_t[0] = 0.5*(1.0 - cos(2*!PI*findgen(num_taper)/(tn*p)))
        w_t[0] = (sin(!PI*findgen(num_taper)/(tn*p)))^2
        w_t = w_t*shift(rotate(w_t, 2), 1) ;the rotate simply reverses the values
        cpg = total(w_t)/n_elements(w_t)
    ENDIF ELSE BEGIN ;defaults to a rectangular window
        cpg = 1.0
    ENDELSE
RETURN, w_t
END

FUNCTION TEST_FFT_WAVE_FUNCTION, amplitude_in, freq_in, phase_in, total_len, dt=dt
    IF NOT keyword_set(dt) THEN dt = 1.0
    wave_func = amplitude_in*cos(2.0*!PI*freq_in*findgen(total_len)*dt+phase_in)
RETURN, wave_func
END

PRO NUWT_APPLY_FFT, res=res, cad=cad, km_per_arcsec=km_per_arcsec, detrend=detrend,$
                   pad_fft=pad_fft, pad_length=pad_length, fill_pad=fill_pad, $
                   bootstrap=bootstrap, num_bootstrap=num_bootstrap, $

                   include_nyquist=include_nyquist, $
                   min_freq_cutoff=min_freq_cutoff, $
                   min_cyc_cutoff=min_cyc_cutoff, $
                   window_func=window_func, $
                   window_param=window_param, $
                   signif_levels=signif_levels, vaughan=vaughan, $

                   old_pow_scale=old_pow_scale, $
                   recover_power=recover_power, adjacent_peaks=adjacent_peaks, $
                   iterate_padding=iterate_padding, vel_amp_mode=vel_amp_mode

;@defaults######################################################################
;LOADING COMMON DATA BLOCKS AND SETTING DEFAULT OPTIONS
;###############################################################################
COMMON threads_dat, threads
COMMON fft_results_dat, fft_spec, fft_peaks

;Setting default values
IF NOT keyword_set(dx) THEN dx = 1.0
IF NOT keyword_set(dt) THEN dt = 1.0
IF NOT keyword_set(pad_length) THEN pad_length = 1000
IF NOT KEYWORD_SET(fill_pad) THEN num_pad_zeros = pad_length
IF KEYWORD_SET(fill_pad) AND NOT KEYWORD_SET(pad_fft) THEN pad_fft = 1
IF NOT keyword_set(num_bootstrap) THEN num_bootstrap = 1000
IF NOT keyword_set(window_func) THEN window_func = 'split_cosine_bell'
IF NOT keyword_set(window_param) THEN window_param = 0.4
IF strlowcase(window_func) EQ 'hann' THEN window_param = 1.0
IF window_param EQ 1.0 THEN window_func = 'hann'
;note: a split cosine bell window with p = 1.0 is the same as the hann window
IF strlowcase(window_func) EQ 'rectangular' THEN window_param = 0.0
IF NOT keyword_set(signif_levels) THEN signif_levels = 0.95
IF keyword_set(adjacent_peaks) THEN allow_adj_pks = 'yes' ELSE allow_adj_pks = 'no'
IF KEYWORD_SET(pad_fft) THEN allow_adj_pks = 'no' ;padding spreads power out and would give incorrect results

If n_elements(min_freq_cutoff) LE 0 THEN min_freq_cutoff = 1e-4
If n_elements(min_cyc_cutoff) LE 0 THEN min_cyc_cutoff = 0.75 ;units of [cycles]

n_threads = N_ELEMENTS(threads)

;catch the case of zero threads found
IF n_threads EQ 1 THEN BEGIN
    only_tlen = threads[0].length
    IF only_tlen LE 0 THEN BEGIN
        PRINT, 'WARNING: no threads were found in the given td-diagram. Returning empty FFT results.'
    ENDIF
ENDIF

;create list to contain the results
fft_spec = LIST()
;each element in fft_spec will contain a structure with the fft results
;for a single thread. See the header for reference of the format.

IF KEYWORD_SET(iterate_padding) THEN BEGIN
    old_fft_peaks = fft_peaks
    pad_fft = 1
    fill_pad = 1
ENDIF

;@units#########################################################################
;COMPUTE UNIT CONVERSION FACTORS
;###############################################################################
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27

IF NOT KEYWORD_SET(res) THEN BEGIN
    res = 1.0
    amp_dx = 1.0
    amp_output_units = 'pixels'
ENDIF ELSE BEGIN
    ;Input resolution must be in [arcsec]!
    amp_dx = res*km_per_arcsec
    amp_output_units = 'km'
ENDELSE

IF NOT KEYWORD_SET(cad) THEN BEGIN
    cad = 1.0
    freq_dt = 1.0
    freq_output_units = 'timesteps^{-1}'
    period_output_units = 'timesteps'
ENDIF ELSE BEGIN
    ;Input cadence must be in [s]!
    freq_dt = 1.0/cad
    freq_output_units = 'Hz'
    period_output_units = 's'
ENDELSE

vel_amp_output_units = amp_output_units+'/'+period_output_units

convert_amp_to_pxls = double(1.0/amp_dx)
convert_freq_to_timesteps = double(1.0/freq_dt)
pow_output_units = amp_output_units+'^{2} '+period_output_units
convert_pow_to_pxls = double(convert_amp_to_pxls^2)

;###############################################
;PERFORM FOURIER ANALYSIS ON EACH GOOD THREAD
;###############################################
FOR h=0L, (n_threads-1) DO BEGIN

    tpos = threads[h].pos[*] ;temp array with postions of the selected thread
    terr = threads[h].err_pos[*] ;errors on postions of thread

    ;Locate the start and end of thread
    t_start = -1000
    t_end = -1000 ;reset each time

    ; ;OLD CODE THAT FOUND THE VALUES FROM THE DATA
    ; loc_th_vals = where(tpos ne -1.0)
    ; len_thread = n_elements(loc_th_vals)
    ; t_start = loc_th_vals[0]
    ; t_end = loc_th_vals[len_thread-1]

    t_start = threads[h].start_bin
    t_end = threads[h].end_bin
    len_thread = threads[h].length

    IF len_thread LE 0 THEN len_thread = 1 ;prevents crashing on a empty thread result

    ; IF KEYWORD_SET(vel_amp_mode) THEN len_thread = threads[h].length -1

    IF KEYWORD_SET(iterate_padding) THEN BEGIN
        ;EXPERIMENTAL iterative method for centering peaks in freq bins
        target_peak_bin = 20
        guess_freq = old_fft_peaks[h].peak_freq[0]*freq_dt ; [Hz]
        IF guess_freq EQ 0.0 THEN guess_freq = 600.0 ;default to prevent crashes
        pad_length = round(target_peak_bin/(guess_freq*cad))
    ENDIF

    IF KEYWORD_SET(fill_pad) THEN num_pad_zeros = pad_length - len_thread
    IF num_pad_zeros LE 0 THEN num_pad_zeros = 1

    ;Determine output array sizes
    IF keyword_set(pad_fft) THEN len_fft_array = len_thread + num_pad_zeros $
        ELSE len_fft_array = len_thread

    IF keyword_set(include_nyquist) THEN BEGIN
        nf = ceil((len_fft_array-1)/2.0) ;INCLUDING Nyquist freq for even N
    ENDIF ELSE BEGIN
        nf = floor((len_fft_array-1)/2.0) ;EXCLUDING Nyquist freq for even N
    ENDELSE

    IF nf LE 0 THEN nf = 1 ;for when there is not enough data to perform an fft

    IF min_cyc_cutoff GT 0.0 THEN BEGIN
        ;Enforcing a hard lower limit on the frequencies returned for peaks
        max_possible_period = float(len_thread)/min_cyc_cutoff ;in units of [pixels]
        min_freq_cutoff = 1.0/max_possible_period ;in units of [pixels^-1]
    ENDIF

    ;Initialize structures to contain the fft results for this thread
    spec_th_fft = {power:fltarr(nf), $
                   err_power:fltarr(nf), $
                   power_units:pow_output_units, $
                   power_to_pxls:convert_pow_to_pxls, $
                   amplitude:fltarr(nf), $
                   err_amplitude:fltarr(nf), $
                   amp_units:amp_output_units, $
                   amp_to_pxls:convert_amp_to_pxls, $
                   freq:fltarr(nf), $
                   freq_units:freq_output_units, $
                   freq_to_timesteps:convert_freq_to_timesteps, $
                   phase:fltarr(nf), $
                   err_phase:fltarr(nf), $
                   trend:fltarr(len_thread), $
                   trend_poly_degree:0, $ ;NEW!
                   apod_window:fltarr(len_thread), $
                   window_func:'not_applicable', $ ;MOVED!
                   window_param:-1.0, $ ;MOVED!
                   signif_vals:fltarr(nf), $
                   fft_length:-1, $ ;MOVED!
                   signif_test:'not_applicable', $ ;MOVED!
                   signif_level:0.0, $ ;COPIED!
                   bin_flags:intarr(nf), $
                   enbw:-1.0, $ ;COPIED!
                   cpg:-1.0} ;MOVED!

    peaks_th_fft =  {analysis_method:'FFT', $ ;NEW!
                     peak_power:[0.0, 0.0, 0.0, 0.0], $
                     err_peak_power:[0.0, 0.0, 0.0, 0.0], $
                     power_units:pow_output_units, $
                     power_to_pxls:convert_pow_to_pxls, $
                     peak_amplitude:[0.0, 0.0, 0.0, 0.0], $
                     err_peak_amplitude:[0.0, 0.0, 0.0, 0.0], $
                     amp_units:amp_output_units, $
                     amp_to_pxls:convert_amp_to_pxls, $
                     peak_freq:[0.0, 0.0, 0.0, 0.0], $
                     freq_units:freq_output_units, $
                     freq_to_timesteps:convert_freq_to_timesteps, $
                     peak_vel_amp:[0.0, 0.0, 0.0, 0.0], $
                     vel_amp_units:vel_amp_output_units, $ ;NEW!
                     peak_phase:[0.0, 0.0, 0.0, 0.0], $
                     err_peak_phase:[0.0, 0.0, 0.0, 0.0], $
                     peak_bin:[-1, -1, -1, -1], $
                     adjacent_peaks:allow_adj_pks, $
                     num_signif_peaks:0, $
                     num_saved_waves:0, $
                     signif_level:0.0, $
                     KS_stat:[-1.0, -1.0, -1.0, -1.0, -1.0], $
                     KS_prob:[-1.0, -1.0, -1.0, -1.0, -1.0], $
                     AD_stat:[-1.0, -1.0, -1.0, -1.0, -1.0], $
                     AD_crit:[-1.0, -1.0, -1.0, -1.0, -1.0], $
                     LB_stat:[-1.0, -1.0, -1.0, -1.0, -1.0], $
                     LB_chisqrd:[-1.0, -1.0, -1.0, -1.0, -1.0], $
                     enbw:-1.0, $
                     user_qual_flag:-1, $
                     auto_qual_flag:1}

    ;Initializing some array values
    spec_th_fft.apod_window[0:-1] = 1.0
    spec_th_fft.bin_flags[0:-1] = -2

    ;Recording some basic information
    spec_th_fft.fft_length = len_fft_array
    spec_th_fft.signif_level = signif_levels
    peaks_th_fft.signif_level = signif_levels

    ;###############################################
    ;SORTING THREADS AND SKIPPING THE BAD ONES
    ;###############################################
    ;skips enteries with less than 2 positive values
    IF (n_elements(where(tpos gt 0.0))) GT 2.0 THEN BEGIN

        ;skips enteries where half the data points are set to zero, i.e. no
        ;value was obtained at fitting stage.
        IF (n_elements(where(tpos EQ 0.0))) LT 0.5*(n_elements(where(tpos GE 0.0))) THEN BEGIN

            ;if value missing set to same as last pixel and set
            ;note: this should NOT be needed if 'patch_up.pro' is used prior to calling this procedure
            zer = where(tpos EQ 0.0)
            IF zer[0] NE -1 THEN BEGIN
                FOR ii=0,n_elements(zer)-1 DO BEGIN
                    tpos[zer[ii]] = tpos[zer[ii]-1]
                    terr[zer[ii]] = 1.0
                ENDFOR
            ENDIF

            ;###############################################
            ;FIT THE GOOD THREADS
            ;###############################################
            IF t_start ge 0 THEN BEGIN

                tpos = tpos[t_start:t_end]
                terr = terr[t_start:t_end]

                IF KEYWORD_SET(vel_amp_mode) THEN BEGIN
                    ;[EXPERIMENTAL] calculate the velocity time series as a means
                    ;to directly determine velocity amplitude independently of freq
                    ;Central difference method (with forwards/backwards on the ends)
                    ; tpos = smooth(tpos, 3, /EDGE_TRUNCATE)
                    ; terr = smooth(terr, 3, /EDGE_TRUNCATE)
                    dv = (SHIFT(tpos, -1) - SHIFT(tpos, 1))*1.0/(2.0)
                    dv[0] = (tpos[1]-tpos[0])*1.0 ;forwards
                    dv[-1] = (tpos[-1]-tpos[-2])*1.0 ;backwards
                    dv_err = (SHIFT(terr, -1) - SHIFT(terr, 1))*1.0/(2.0)
                    dv_err[0] = (terr[1]-terr[0])*1.0
                    dv_err[-1] = (terr[-1]-terr[-2])*1.0
                    tpos = dv
                    terr = dv_err
                ENDIF

                ; calculates linear trend or mean
                IF keyword_set(detrend) THEN BEGIN
                    line_fit_result = poly_fit(findgen(len_thread), tpos, 1, yfit=trend) ;fit linear trend line
                    spec_th_fft.trend_poly_degree = 1
                ENDIF ELSE BEGIN
                    trend = fltarr(len_thread) + mean(tpos) ;creates an array with the mean of the data
                    spec_th_fft.trend_poly_degree = 0
                ENDELSE

                spec_th_fft.trend = trend

                ;###############################################
                ;Perform the FFT and calculate the power and phase
                ;###############################################

                oscil = tpos-trend ;detrend the time series (or just subtract the mean from it)

       		    s_len = n_elements(oscil) ;length of the SOURCE (S), unpadded time-series array
       		    apodt = get_apod_window(s_len, cpg=cpg, parameter=window_param, win_function=window_func)
                ;note: cpg is Coherent Power Gain - correction factor needed for power after apodisation

                spec_th_fft.apod_window = apodt
                spec_th_fft.window_func = window_func
                spec_th_fft.window_param = window_param

                ;applying window and padding the array (if asked to)
                IF KEYWORD_SET(bootstrap) THEN BEGIN
                    ;bootstrapping the input data. Returns 2D array with dimensions of [s_len, num_bootstrap]
                    boot_arr = bootstrap_1d_data(oscil, terr, num_resamples = num_bootstrap)

           		    IF NOT keyword_set(pad_fft) THEN BEGIN
                        oscil = oscil*apodt
                        boot_arr = boot_arr*rebin(apodt, s_len, num_bootstrap)
                  	    n_len = s_len
                    ENDIF ELSE BEGIN
                        ; oscil = smooth(boot_arr, [3, 1], /edge_truncate)
                        oscil = [oscil*apodt , fltarr(num_pad_zeros)]
                        boot_arr = smooth(boot_arr, [3, 1], /edge_truncate)
                      	boot_arr = [boot_arr*rebin(apodt, s_len, num_bootstrap), fltarr(num_pad_zeros, num_bootstrap)]
                      	n_len = s_len + num_pad_zeros
           		    ENDELSE
                ENDIF ELSE BEGIN
           		    IF NOT keyword_set(pad_fft) THEN BEGIN
                        oscil = oscil*apodt
                  	    n_len = n_elements(oscil)
                    ENDIF ELSE BEGIN
                        ; oscil = smooth(oscil,3, /edge_truncate)
                      	oscil = [oscil*apodt, fltarr(num_pad_zeros)]
                      	n_len = n_elements(oscil)
           		    ENDELSE
                ENDELSE

                IF keyword_set(include_nyquist) THEN BEGIN
                    end_fft_index = ceil((n_len-1)/2.0) ;last pos freq result in fft (INCLUDING Nyquist freq for even N)
                ENDIF ELSE BEGIN
                    end_fft_index = floor((n_len-1)/2.0) ;last pos freq result in fft (EXCLUDING Nyquist freq for even N)
                ENDELSE
                df = 1.0/(n_len);*dt)
                f = findgen(n_len)*df
                f = f[1:end_fft_index]

                ;Take FFT and then calculate power (PSD) and phase
                IF KEYWORD_SET(bootstrap) THEN BEGIN
                    fft_of_boot = fft(boot_arr, -1, dimension=1)
                    boot_pow = (2.0*(abs(fft_of_boot))^2)[1:end_fft_index,*]
                    boot_phase = atan((fft_of_boot)[1:end_fft_index,*], /phase)

                    pow = mean(boot_pow, dimension=2)
                    phase = mean(boot_phase, dimension=2)
                    ; pow = median(boot_pow, dimension=2)
                    ; phase = median(boot_phase, dimension=2)
                    pow_err = stddev(boot_pow, dimension=2)
                    phase_err = stddev(boot_phase, dimension=2)
                ENDIF ELSE BEGIN
                    fft_of_oscill = fft(oscil, -1)
                    pow = (2.0*(abs(fft_of_oscill))^2)[1:end_fft_index]
                    phase = atan((fft_of_oscill)[1:end_fft_index], /phase)
                ENDELSE

                IF NOT KEYWORD_SET(old_pow_scale) THEN BEGIN
                    ;[DEFAULT] Power and Amp. scaling as given by Heinzel FFT Notes
                    S1 = total(apodt)
                    S2 = total(apodt^2)
                    power_correction = (cad)*(n_len^2)/S2 ;the 12.0 is due to a factor of 1/fs where fs = 1/cad
                    amp_correction = n_len/S1;*sqrt(2.0)/2.0 ;do we need the 2.0 factors to correct for an error?
                    signif_correction = ((cad)*(S1^2))/S2
                ENDIF ELSE BEGIN
                    ;[OLD METHOD] - Assumes power is evenly distributed in apod window (not always true)
                    ; note: CPG below is equal to S1/s_len (see above)
                    S1 = total(apodt)
                    S2 = total(apodt^2)
                    power_correction = (n_len/(s_len*CPG))^2
                    amp_correction = n_len/(s_len*CPG)
                    signif_correction = 1.0
                ENDELSE

                IF end_fft_index LT (nf - 1) THEN BEGIN
                    spec_th_fft.bin_flags[end_fft_index-1:-1] = -1 ;flag fft bins outside results with a value of -1
                ENDIF

                ;###############################################
                ;Calculates significance level of the FFT signal
                ; (see e.g., Torrence & Compo 1998 and Vaughan 2005)
                ;###############################################
        	    ;Doesn't let straight lines through (i.e. fft results with very low varience)
       		    IF ((moment(pow))[1]) GT 1e-30 THEN BEGIN
                ; IF ((moment(pow))[1]) GT 1e-9 THEN BEGIN ;OLD threashold that does not work well for large fft padding

                    ;Computing various results and saving to the output structure
                    ;UPDATE [2017-03-09]: now includes unit conversions to output physical units by default
                    spec_th_fft.power = pow*power_correction*(amp_dx^2)
                    spec_th_fft.amplitude = 2.0*sqrt(pow/2.0)*amp_correction*amp_dx
           			spec_th_fft.freq = f*freq_dt
                    spec_th_fft.phase = phase
                    spec_th_fft.enbw = (1.0/cad)*(S2/S1^2)
                    spec_th_fft.cpg = CPG
                    peaks_th_fft.enbw = (1.0/cad)*(S2/S1^2)

                    IF KEYWORD_SET(bootstrap) THEN BEGIN
                        spec_th_fft.err_power = pow_err*power_correction*(amp_dx^2)
                        spec_th_fft.err_amplitude = 2.0*sqrt(pow_err/2.0)*amp_correction*amp_dx
                        spec_th_fft.err_phase = phase_err
                    ENDIF

                    ;Significant tests for power spectra
         			IF keyword_set(vaughan) THEN BEGIN
                        ; Alternative from Vaughan 2005, Astron. & Astophys. (A&A)
                        ; Note: this method assumes a red-noise spectrum
                        ; TEMP note: Need to still add unit conversions (and test signif corrections)
                        spec_th_fft.signif_test = 'Vaughan_2005'
           			    ;Normalise power
         				npw = s_len*pow/(moment(oscil))[1]
            			prob = 1.0 - chisqr_pdf(2.0*npw,2)
            			nprob = MAKE_ARRAY(SIZE(prob, /DIM))
            			log_pp = s_len * ALOG(DOUBLE(1.0 - prob))
            			indx = WHERE(log_pp gt -30.0, count, COMPLEMENT=indx_c)

                        IF (count gt 0) THEN nprob[indx] = exp(log_pp[indx])
            			IF (count lt s_len) THEN nprob[indx_c] = 0.0

                        spec_th_fft.signif_vals = nprob
                        loc_signif_pow = where(nprob gt signif_levels, /NULL)
         			ENDIF ELSE BEGIN
                        ; [DEFAULT] from Torrence & Compo 1998, Bul. Amer. Met. Soc. (BAMS)
                        ; Note: can select from white- or red-noise (normally just use white)
                        spec_th_fft.signif_test = 'Torrence_&_Compo_1998'
                        sig_vals = SIGNIF_NOISE_SPEC((tpos-trend), signif_levels, n_elements(f), color='white', /bonferroni)
                        spec_th_fft.signif_vals = sig_vals*signif_correction*(amp_dx^2)
                        loc_signif_pow = where(spec_th_fft.power GT spec_th_fft.signif_vals, /NULL)

                        ;[OLD] Method from before changing the PSD scaling (Dec. 2016)
            			; sig_vals = SIGNIF_CONF(oscil, signif_levels)
                        ; sig_vals = SIGNIF_NOISE_SPEC(oscil, signif_levels, n_elements(f), color='white')
                        ; sig_vals = SIGNIF_NOISE_SPEC_BOOTSTRAP(boot_arr, signif_levels, n_elements(f))
                        ; boot_sig_vals = SIGNIF_NOISE_SPEC(mean(boot_arr, dimension=2), signif_levels, n_elements(f))
                        ; print, 'thread #', h, ' sig val diff.'
                        ; print, boot_sig_vals - sig_vals
                        ; spec_th_fft.signif_vals = sig_vals*power_correction
            			; loc_signif_pow = where(pow gt sig_vals, /NULL)
        			ENDELSE

                    ; Flagging significant power bins
                    spec_th_fft.bin_flags[0:end_fft_index-1] = 0 ;fill all valid fft bin flags with a value of 0
                    spec_th_fft.bin_flags[loc_signif_pow] = 1 ;flag significant results with a quality value of 1

                    ;Finding peaks in the significant power results (NEW METHOD AS OF JUNE 2016)
                    temp_pow = spec_th_fft.power[0:-1] ;temporary array for finding peaks
                    temp_pow[where(spec_th_fft.bin_flags LT 1, /NULL)] = 0 ;clear out values that are not "significant"
                    compare_right = temp_pow - SHIFT(temp_pow, -1)
                    compare_right[-1] = 0 ;ignores the value to the "right" of the last value in the array
                    compare_left = temp_pow - SHIFT(temp_pow, 1)
                    compare_left[0] = 0 ;ignores the value to the "left" of the last value in the array
                    IF allow_adj_pks EQ 'no' THEN BEGIN
                        loc_signif_peaks = where(temp_pow GT 0 and (compare_right GE 0 and compare_left GE 0), /NULL)
                    ENDIF ELSE BEGIN
                        loc_signif_peaks = where(temp_pow GT 0, /NULL)
                    ENDELSE

                    spec_th_fft.bin_flags[loc_signif_peaks] = 2 ;flag local peaks in the significant results with a value of 2

                    ;Saving the wave parameters of the largest significant peaks
                    ;note: will only save the FOUR largest peaks, even if there are more
                    peaks_th_fft.num_signif_peaks = n_elements(loc_signif_peaks)
                    sav_p = 0 ;counter & index for Saved_Peaks
                    IF peaks_th_fft.num_signif_peaks GT 0.0 THEN BEGIN
                        ;sort based on the power at each peak
         			    pow_at_signif_peaks = spec_th_fft.power[loc_signif_peaks]
                        sorted_indices = REVERSE(SORT(pow_at_signif_peaks)) ;sort values high to low
                        FOR jj=0, peaks_th_fft.num_signif_peaks-1 DO BEGIN
                            bin_of_peak = loc_signif_peaks[sorted_indices[jj]]
                            IF (sav_p LT 4) AND (f[bin_of_peak] GT min_freq_cutoff) THEN BEGIN
                                ;save only realistic values  (i.e. greater than a min frequencey)
                                peaks_th_fft.peak_power[sav_p] = spec_th_fft.power[bin_of_peak]
                                peaks_th_fft.peak_amplitude[sav_p] = spec_th_fft.amplitude[bin_of_peak]
                                peaks_th_fft.peak_freq[sav_p] = spec_th_fft.freq[bin_of_peak]
                                peaks_th_fft.peak_vel_amp[sav_p] = 2*!PI*spec_th_fft.amplitude[bin_of_peak]*spec_th_fft.freq[bin_of_peak]
                                peaks_th_fft.peak_phase[sav_p] = phase[bin_of_peak]
                                peaks_th_fft.peak_bin[sav_p] = bin_of_peak
                                IF KEYWORD_SET(bootstrap) THEN BEGIN
                                    peaks_th_fft.err_peak_power[sav_p] = spec_th_fft.err_power[bin_of_peak]
                                    peaks_th_fft.err_peak_amplitude[sav_p] = spec_th_fft.err_amplitude[bin_of_peak]
                                    peaks_th_fft.err_peak_phase[sav_p] = phase_err[bin_of_peak]
                                ENDIF
                                sav_p = sav_p + 1 ;count number of wave results saved so far
                            ENDIF
                        ENDFOR
                    ENDIF

                    peaks_th_fft.num_saved_waves = sav_p
                    ; IF peaks_th_fft.num_signif_peaks GT sav_p THEN peaks_th_fft.num_signif_peaks = sav_p

                    ; [EXPERIMENTAL] Recovering leaked PSD
                    ; note: this will ONLY work if adjacent peaks are disallowed
                    ;       currently only set to recover power from directly adjacent bins
                    IF KEYWORD_SET(recover_power) AND allow_adj_pks EQ 'no' THEN BEGIN
                        print, 'Attempting to recover leaked power in FFT!'
                        FOR w=0, 3 DO BEGIN
                            center_bin = peaks_th_fft.peak_bin[w]
                            IF center_bin GT -1 THEN BEGIN
                                left_bin = center_bin - 1
                                right_bin = center_bin + 1
                                ; cap bins to end of array
                                IF left_bin LT 0 THEN left_bin = 0
                                IF right_bin GT (nf-1) THEN right_bin = nf-1
                                num_recover_bins = 1 + right_bin - left_bin

                                ; Get weights for sum (will share power between peaks adjacent to same bin)
                                sum_weights = fltarr(num_recover_bins) + 1.0
                                IF center_bin GE 2 THEN BEGIN
                                    IF spec_th_fft.bin_flags[left_bin-1] EQ 2 THEN sum_weights[0] = 0.5
                                ENDIF
                                IF center_bin LE (nf-3) THEN BEGIN
                                    IF spec_th_fft.bin_flags[right_bin+1] EQ 2 THEN sum_weights[-1] = 0.5
                                ENDIF

                                ; sum power, average freq, and average phase
                                tot_power = total(spec_th_fft.power[left_bin:right_bin]*sum_weights)
                                mean_weights = (spec_th_fft.power[left_bin:right_bin]/tot_power)*sum_weights
                                mean_freq = total(f[left_bin:right_bin]*mean_weights)/total(mean_weights)
                                mean_phase= total(phase[left_bin:right_bin]*mean_weights)/total(mean_weights)

                                ; recalc amplitude (should already have the correct power correction)
                                tot_amp = 2.0*sqrt(tot_power/2.0)

                                ; overwrite the saved results with these adjusted values
                                peaks_th_fft.peak_power[w] = tot_power
                                peaks_th_fft.peak_amplitude[w] = tot_amp
                                peaks_th_fft.peak_freq[w] = mean_freq
                                peaks_th_fft.peak_phase[w] = mean_phase
                            ENDIF
                        ENDFOR
                    ENDIF

                    ; GOF tests of residuals
                    combined_wave_vals = fltarr(s_len)
                    total_fit_params = 0
                    ; note: the results for the first N waves combined is stored in the i = N position.
                    ;       In other words, the 1st result (index 0) is the stats for NO waves,
                    ;       the 2nd (index 1) results is for the first wave, the 3rd (index 2)
                    ;       is for the first TWO waves inclusive, and so on and so forth...
                    FOR w=0, peaks_th_fft.num_saved_waves DO BEGIN
                        IF w GT 0 AND peaks_th_fft.peak_bin[w-1] GT -1 THEN BEGIN
                            add_wave = TEST_FFT_WAVE_FUNCTION(peaks_th_fft.peak_amplitude[w-1], peaks_th_fft.peak_freq[w-1], $
                                                              peaks_th_fft.peak_phase[w-1], s_len)
                            combined_wave_vals = combined_wave_vals + add_wave
                            total_fit_params += 3
                        ENDIF
                        residuals = tpos - (combined_wave_vals + trend)
                        resid_stats = test_residuals(residuals, num_fit_params=total_fit_params, $
                                                     signif_level=peaks_th_fft.signif_level)
                        peaks_th_fft.KS_stat[w] = resid_stats.KS_stat
                        peaks_th_fft.KS_prob[w] = resid_stats.KS_prob
                        peaks_th_fft.AD_stat[w] = resid_stats.AD_stat
                        peaks_th_fft.AD_crit[w] = resid_stats.AD_crit
                        peaks_th_fft.LB_stat[w] = resid_stats.LB_stat
                        peaks_th_fft.LB_chisqrd[w] = resid_stats.LB_chisqrd
                    ENDFOR

       		    ENDIF
            ENDIF
        ENDIF
    ENDIF

    ;Add the fft results to output LIST of structures (with variable size array tags)
    fft_spec.ADD, spec_th_fft

    ;Create or concat fft stats to output ARRAY of structures (with fixed size and shape)
    IF h EQ 0 THEN BEGIN
        fft_peaks = [peaks_th_fft]
    ENDIF ELSE BEGIN
        fft_peaks = [temporary(fft_peaks), peaks_th_fft]
    ENDELSE
ENDFOR
END
