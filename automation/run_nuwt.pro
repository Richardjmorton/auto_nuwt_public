;+
;NAME: RUN_NUWT
;
;PURPOSE:
;   One of the top level Northumbria University Wave Tracking (NUWT) run scripts.
;   Takes a time distance diagram, picks out peaks in the intensity (also called
;   'threads'), then feeds them into an FFT to calculate the power and frequencies
;   of oscillation (if any). See the user guide (once available) and the code of
;   each sub-program for more information.
;
;INPUTS:
;   input_data - an array of time-distance (t-d) diagrams of the form (x,t,z)
;
;OPTIONAL INPUTS:
;   errors - estimated errors on intensity values, supplied to gaussian fitting routine.
;            if not supplied a default value of 10% of intensity values is used
;   /invert - If set, will invert the intensity values of the td-diagram before
;             finding peak locations. This results in the code finding the local
;             MINIMA in the original image.
;   grad - gradient limit used by "locate_things.pro". Default is 0.5 (good for
;          unsharp masked data)
;   min_tlen - minimum thread length used by "follow_threads.pro". Default is 20
;   max_dist_jump - maximum allowable distance (in pixels) between two peaks in
;                   the same thread. Used by "follow_threads.pro". Default is 3 pixels
;   max_time_skip - maximum allowable timesteps beweewn two peaks in the same
;                   thread. Used BY "follow_threads.pro". Default is 4 time steps.
;   res - spatial resolution of data in [arcsec]. Defaults to a value of 1 [pixel]
;   cad - temporal cadence of the data in [s]. Defaults to a value of 1 [timestep]
;   km_per_arcsec - ratio of [km]/[arcsec]. Defaults to 725.27 which is the
;                   mean scale of the solar surface.
;   /gauss - uses a gaussian fit to locate thread centres; provides sub-pixel measurements
;   /full_gauss - similar to the above. Will also return the fitted gaussian widths (for debugging)
;   /pad_fft - if set, will pad the array with zeros to improve the precision
;              of peak frequencies returned by the fft. Note, this does NOT
;              actually increase the resolution and can results in extra spectral
;              leakage. Use with care.
;   /fill_pad - if set, will fill the array with extra zeros until the length equals
;               "pad_length" instead of appending a set number zeros (default).
;               Note: if the pad_length is set too low, using /fill_pad may
;               introducing "banding" in the possible frequency values returned.
;   pad_length - number of zeros to pad to the end of the array when /pad_fft is set.
;                If /fill_pad is ALSO set, will instead pad the array until the
;                total length equals pad_length (or just pad with 1 zero if the
;                array is already longer than pad_length). Default is 1000.
;   /bootstrap - if set, will perform bootstrapping to resample the thread positions
;                before running the FFT. This gives a rough error estimates on the
;                FFT spectrum and output wave paramters
;   num_bootstrap - number of bootstrap resamples to use. Default is 1000
;   /vel_amp_mode - [EXPERIMENTAL] alternate mode for "nuwt_apply_fft.pro" (see
;                   that code for more information)
;   slit_meta - metadata structure outputted by "nuwt_arc_slit.pro" or "nuwt_diag_slit.pro"
;   /aia - if set, will assume input data came from SDO / AIA and proceed to use
;          a few convenience functions.
;   wavelength - wavelength of input AIA data from calculating intensity errors.
;   slit_norm_num - width of data slit. Alos used for calcuating AIA intensity errors
;   /interactive - if set, will allow for interatively picking gradient and min thread length
;   save_folder - folder in which to save the results. Defaults to current folder.
;                 Will also append a '/' to the end if not included.
;   filetag - unique tag to append to the START of the output file
;   default_plot_header - string that can be optioanlly defined an used when
;                         plotting the data using the NUWT plot procedures
;
;OUTPUTS:
;   common blocks contain the results and are used for all internal passing of
;   data between subprograms.
;
;EXTERNAL CALLS - locate_things.pro, follow_thread.pro, patch_up.pro, nuwt_apply_fft.pro
;
;HISTORY: Name---------Date---------Description
;         M Weberg  JUNE, 2016  Initial coding.
;         M weberg  SEPT, 2016  Made some quality-of-life improvements
;                               - Added a master COMMON block that compiles the
;                                 results from all of the data slits into lists.
;                                 Can also output this block to a '.save' file
;                               - Added an interactive mode for selecting gradients
;                                 and minimum thread lengths
;         M Weberg   FEB, 2017  Added vel_amp_mode options
;         M Weberg   MAR, 2017  Updated all NUWT procedures to automatically
;                               calculate parameters in physical units when given
;                               the "res" and "cad" for the input data. Also
;                               added some convenient defaults when using aia data.
;         M Weberg   MAY, 2017  Added slit_meta functions
;         M Weberg  JUNE, 2017  Added the "nuwt_meta_dat" COMMON block and structure
;         M Weberg  JULY, 2017  Exposed the "max_dist_jump" & "max_time_skip" options
;                               from "follow_thread.pro" to be changable via inputs
;                               to this script
;         M Weberg  SEPT, 2017  Renamed "nuwt_fft_results" to "nuwt_fft_spec" in
;                               all NUWT programs (more appropriate name)
;         M Weberg   JAN, 2018  Renamed a number of subprograms and output structure
;                               to have more intuitive names (most importanly,
;                               "nuwt_fft_stats" is now "nuwt_fft_peaks" and
;                               "plot_nuwt_fft_peaks" is now "plot_nuwt_wave_hist")
;
;TO DO / RESTRICTIONS:
;   - Add funcationality for passing a structure with the desired keyword
;     arguments (kwargs) for each subprogram.
;   - Add options for plotting the data after processing
;-

PRO RUN_NUWT, input_data, errors=errors, $

              invert=invert, grad=grad, $
              smooth=smooth, sm_width=sm_width, $

              min_tlen=min_tlen, $
              max_dist_jump=max_dist_jump, $
              max_time_skip=max_time_skip, $

              gauss=gauss, full_gauss=full_gauss, $
              pad_fft=pad_fft, fill_pad=fill_pad, pad_length=pad_length, $
              bootstrap=bootstrap, num_bootstrap=num_bootstrap, $
              vel_amp_mode=vel_amp_mode, $

              res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
              slit_meta=slit_meta, aia=aia, wavelength=wavelength, $
              slit_norm_num=slit_norm_num, $
              interactive=interactive, $
              save_folder=save_folder, $
              filetag=filetag, $
              default_plot_header=default_plot_header

COMPILE_OPT IDL2
nuwt_code_version_num = 2.21
PRINT, ' ' ;just an empty line to make console output easier to read
PRINT, '##### Running Northumbria Wave Tracking Code (NUWT) #####'
;@check_data####################################################################
; CHECK THE INPUT DATA AND EXTRACT VARIABLES FROM SLIT STRUCTURES
;###############################################################################
IF N_TAGS(input_data) GT 0 THEN BEGIN
    ;Input is a structure containing both slit(s) and meta data
    input_tags = TAG_NAMES(input_data)
    IF TOTAL(STRCMP(input_tags, 'slit', /fold_case)) GE 1 THEN BEGIN
        img_data = input_data.slit
        slit_meta = input_data.meta
    ENDIF ELSE BEGIN
        PRINT, 'ERROR: invalid structure format!'
        MESSAGE, 'Input data strutures must have both "slit" and "meta" tags'
    ENDELSE
ENDIF ELSE BEGIN
    ;Input is an array with one (or more) data slits
    img_data = input_data
ENDELSE

;@int_common####################################################################
; INITIALIZATION OF COMMON BLOCKS
;###############################################################################
;Structures, arrays, and list containing the results for a single data slit
COMMON located_dat, located
COMMON threads_dat, threads
COMMON fft_results_dat, fft_spec, fft_peaks

;Create a master COMMON block that holds the results for ALL slits
COMMON all_nuwt_dat, nuwt_located, nuwt_threads, nuwt_fft_spec, nuwt_fft_peaks
COMMON nuwt_meta_dat, nuwt_meta, nuwt_bulk_stats ;assorted useful meta data about the nuwt run
nuwt_located = LIST()
nuwt_threads = LIST()
nuwt_fft_spec = LIST()
nuwt_fft_peaks = LIST()
nuwt_bulk_stats = LIST()

;@set_defaults##################################################################
; SETTING DEFAULT PARAMETER VALUES
;###############################################################################
;Setting variables which control the search box size
IF KEYWORD_SET(invert) THEN invert = 1 ELSE invert = 0
IF KEYWORD_SET(smooth) THEN smooth = 1 ELSE smooth = 0
IF N_ELEMENTS(sm_width) EQ 0 THEN sm_width = [1, 3]
IF N_ELEMENTS(sm_width) EQ 1 THEN sm_width = [0, sm_width]
IF N_ELEMENTS(sm_width) GT 2 THEN sm_width = [sm_width[0], sm_width[1]]
IF N_ELEMENTS(grad) EQ 0 THEN grad = 0.5

IF N_ELEMENTS(min_tlen) EQ 0 THEN min_tlen = 20
IF NOT KEYWORD_SET(max_dist_jump) THEN max_dist_jump = 3
IF NOT KEYWORD_SET(max_time_skip) THEN max_time_skip = 4

IF NOT KEYWORD_SET(res) THEN res = 0
IF NOT KEYWORD_SET(cad) THEN cad = 0
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27

IF KEYWORD_SET(gauss) OR KEYWORD_SET(full_gauss) THEN gauss = 1 ELSE gauss = 0
IF KEYWORD_SET(full_gauss) THEN full_gauss = 1 ELSE full_gauss = 0

IF KEYWORD_SET(pad_fft) THEN pad_fft = 1 ELSE pad_fft = 0
IF KEYWORD_SET(fill_pad) THEN fill_pad = 1 ELSE fill_pad = 0
IF KEYWORD_SET(fill_pad) AND NOT KEYWORD_SET(pad_fft) THEN pad_fft = 1
IF NOT KEYWORD_SET(pad_length) THEN pad_length = 1000
IF KEYWORD_SET(bootstrap) THEN bootstrap = 1 ELSE bootstrap = 0
IF NOT KEYWORD_SET(num_bootstrap) THEN num_bootstrap = 1000

IF KEYWORD_SET(vel_amp_mode) THEN vel_amp_mode = 1 ELSE vel_amp_mode = 0

IF KEYWORD_SET(filetag) THEN BEGIN
    save_filename = filetag+'_nuwt_results.sav'
ENDIF ELSE BEGIN
    save_filename='nuwt_results_run_'+GET_CURRENT_TIMESTAMP()+'.sav' ;YYYYMMDD_hhmm format
ENDELSE
IF NOT KEYWORD_SET(save_folder) THEN CD, CURRENT=save_folder ;i.e. defaults to current folder
IF STRLEN(save_folder) GT 0 AND NOT save_folder.endswith('/') THEN save_folder = save_folder+'/'
IF NOT KEYWORD_SET(default_plot_header) THEN default_plot_header = 0

img_size = SIZE(img_data)
IF img_size[0] EQ 3 THEN nslits = img_size[3] ELSE nslits = 1

;check the shape and size of the input error values
IF N_ELEMENTS(errors) EQ 0 THEN BEGIN
    PRINT, 'WARNING: No errors input. NUWT will assume errors = 10% of the input data.'
ENDIF ELSE BEGIN
    err_size = SIZE(errors)
    IF err_size[0] NE img_size[0] OR TOTAL(err_size NE img_size) GE 1 THEN BEGIN
        MESSAGE, 'ERROR: "errors" array must have the same dimensions as the input data!'
    ENDIF
ENDELSE

IF NOT KEYWORD_SET(gauss) THEN BEGIN
    PRINT, 'Finding whole pixel maxima locations'
    maxima_mode = 'nearest whole pixel'
ENDIF ELSE BEGIN
    PRINT, 'Gaussian location enabled. Finding sub-pixel maxima locations'
    maxima_mode = 'sub-pixel Gaussian fit'
ENDELSE

IF KEYWORD_SET(slit_meta) THEN BEGIN
    IF NOT KEYWORD_SET(res) THEN res = slit_meta.res
    IF NOT KEYWORD_SET(res) THEN cad = slit_meta.cad
    IF NOT KEYWORD_SET(res) THEN km_per_arcsec = slit_meta.km_per_arcsec
    IF NOT KEYWORD_SET(res) THEN wavelength = slit_meta.wavelength
    ; IF NOT KEYWORD_SET(res) THEN slit_norm_num = slit_meta.width
ENDIF

IF KEYWORD_SET(aia) THEN BEGIN
    IF NOT KEYWORD_SET(res) THEN res = 0.6
    IF NOT KEYWORD_SET(cad) THEN cad = 12.0
    IF NOT KEYWORD_SET(wavelength) THEN wavelength = 171
    ; IF NOT KEYWORD_SET(slit_norm_num) THEN slit_norm_num = 1
    IF N_ELEMENTS(errors) EQ 0 THEN BEGIN
        PRINT, 'NOTICE: AIA errors are no longer calculated in run_nuwt.pro'
        PRINT, '        Please use the function CALC_AIA_ERRORS() before running NUWT.'
    ENDIF
    ; errors = calc_aia_errors(img_data, wavelength, norm_num=slit_norm_num)
ENDIF

;@query_inputs##################################################################
; QUERY USER FOR GRAD AND MIN THREAD LENGTH PARAMETERS (IF NOT GIVEN)
;###############################################################################
IF KEYWORD_SET(interactive) THEN BEGIN
    ;Initialize the threads structure (prevents odd behavior in interactive mode)
    IF NOT KEYWORD_SET(full_gauss) THEN BEGIN
        ;standard mode (either sub-pixel fitting OR nearest-pixel)
        threads = {pos:FLTARR(10), err_pos:FLTARR(10), $
                   bin_flags:INTARR(10), $
                   start_bin:-1, end_bin:-1, length:-1}
    ENDIF ELSE BEGIN
        ;/full_gauss mode
        threads = {pos:fltarr(10), err_pos:fltarr(10), $
                   inten:fltarr(10), err_inten:fltarr(10), $
                   wid:fltarr(10), err_wid:fltarr(10), $
                   bin_flags:intarr(10), $
                   start_bin:-1, end_bin:-1, length:-1}
    ENDELSE

    PRINT, 'Interactive mode activated!'
    PRINT, 'NOTE: plots will automatically close after each command prompt is answered'
    PRINT, ' ' ;and another...
    grad = 1d
    READ, grad, PROMPT='Please enter initial gradient to try: '

    min_tlen = 1d
    READ, min_tlen, PROMPT='Please enter initial minimum thread length to try: '
ENDIF

;@int_nuwt_meta#################################################################
; RECORD IMPORTANT INFORMATION TO THE "NUWT_META" STRUCTURE
;###############################################################################
IF NOT KEYWORD_SET(slit_meta) THEN slit_meta = 'unknown'
nuwt_meta = {run_date:SYSTIME(), nuwt_version:nuwt_code_version_num, $
             res:res, cad:cad, km_per_arcsec:km_per_arcsec, $
             num_slits:nslits, $
             grad:FLTARR(nslits), $
             inverted_img:invert, $
             smoothed_img:smooth, $
             min_tlen:INTARR(nslits), $
             max_dist_jump:INTARR(nslits), max_time_skip:INTARR(nslits), $
             num_threads:INTARR(nslits), num_waves:INTARR(nslits), $
             maxima_fitting:maxima_mode, $
             pad_fft:pad_fft, fill_pad:fill_pad, pad_length:pad_length, $
             bootstrap:bootstrap, num_bootstrap:num_bootstrap, $
             vel_amp_mode:vel_amp_mode, $
             default_plot_header:default_plot_header, $
             slit_meta:slit_meta}

;@loop_slits####################################################################
; PROCESS EACH TIME-DISTANCE DIAGRAM
;###############################################################################
PRINT, '   ' ;pentultimate empty line printed
FOR k=0, nslits-1 DO BEGIN
    PRINT, 'Processing Slit Number: '+STRTRIM(k, 2)

    ;Finding peak locations
    PRINT, '   Finding local intensity maxima ...'
    ;load or roughly estimate the intensity errors
    IF N_ELEMENTS(errors) GT 0 THEN slit_inten_err = errors[*,*,k]

    ;If not user-supplied, errors are assumed to be 10% of the intensity
    IF N_ELEMENTS(errors) EQ 0 THEN BEGIN
        slit_inten_err = FLTARR(img_size[1],img_size[2])
        slit_inten_err = ABS(img_data[*,*,k])*0.1
    ENDIF

    loop_grad = 1
    ask_grad = 'unknown'
    WHILE loop_grad DO BEGIN
        IF (NOT KEYWORD_SET(gauss)) OR KEYWORD_SET(interactive) THEN BEGIN
            ;Whole pixel resolution (no gaussian fitting)
            ;Also used as part of the interactive mode since the sub-pixel fitting
            ;takes a long time to run.
            NUWT_LOCATE_THINGS, img_data[*,*,k], grad=grad, res=res, cad=cad, $
                                km_per_arcsec=km_per_arcsec, errors=slit_inten_err, $
                                /nearest_pixel, invert=invert, /despike
        ENDIF ELSE BEGIN
            ;Sub-pixel resolution using gaussin fitting of intensity maxima
            NUWT_LOCATE_THINGS, img_data[*,*,k], grad=grad, res=res, cad=cad, $
                                km_per_arcsec=km_per_arcsec, errors=slit_inten_err, $
                                invert=invert, /despike, $
                                smooth=smooth, sm_width=sm_width, $
                                full_gauss=full_gauss, cut_chisq=1e12 ;removed chi-sqrd
        ENDELSE

        IF KEYWORD_SET(interactive) THEN BEGIN
            PLOT_NUWT_PEAKS, img_data[*,*,k], slitnum=k, grad=grad, $
                             /plot_peaks, /screen, /use_temp_common_blocks
            PRINT, ''
            PRINT, 'Current gradient = '+STRTRIM(grad, 2)
            READ, ask_grad, PROMPT='Would you like to select a different gradient? (y/n) '
            ask_grad = STRLOWCASE(ask_grad)
            IF ask_grad.startswith('n') THEN BEGIN
                loop_grad = 0
            ENDIF ELSE BEGIN
                READ, grad, PROMPT='Please enter a new gradient to try: '
            ENDELSE
            ;close the figures afterwards
            fig_td = GETWINDOWS('fig_td')
            fig_peaks = GETWINDOWS('fig_peaks')
            fig_grads = GETWINDOWS('fig_grads')
            IF fig_td THEN fig_td.close
            IF fig_peaks THEN fig_peaks.close
            IF fig_grads THEN fig_grads.close
        ENDIF ELSE BEGIN
            ;end 'loop' if not in interactive mode
            loop_grad = 0
        ENDELSE
    ENDWHILE

    IF KEYWORD_SET(interactive) AND KEYWORD_SET(gauss) THEN BEGIN
        ;Sub-pixel resolution using gaussin fitting of intensity maxima
        PRINT, '   Running the sub-pixel resolution method with the final gradient ...'
        NUWT_LOCATE_THINGS, img_data[*,*,k], grad=grad, res=res, cad=cad, $
                            km_per_arcsec=km_per_arcsec, errors=slit_inten_err, $
                            invert=invert, /despike, $
                            smooth=smooth, sm_width=sm_width, $
                            full_gauss=full_gauss, cut_chisq=1e12 ;removed chi-sqrd
    ENDIF

    ;Save final used values to the meta structure
    nuwt_meta.grad[k] = grad

    ;Identify connected 'threads' of peaks and track them over time
    PRINT, '   Identifing and tracking threads ...'
    loop_tlen = 1
    ask_tlen = 'unknown'
    WHILE loop_tlen DO BEGIN
        NUWT_FOLLOW_THREADS, min_tlen=min_tlen, max_dist_jump=max_dist_jump, $
                             max_time_skip=max_time_skip

        n_threads = N_ELEMENTS(threads)

        ;Catch the case when no threads are found
        IF n_threads EQ 1 THEN BEGIN
            only_tlen = threads[0].length
            IF only_tlen LE 0 THEN BEGIN
                PRINT, '   No threads found for the selected parameters!'
                PRINT, '   Please try running NUWT with different grad or min_tlen values.'
                n_threads = 0
            ENDIF
        ENDIF ELSE BEGIN
            PRINT, '   Number of threads found: '+STRTRIM(n_threads, 2)
        ENDELSE

        IF KEYWORD_SET(interactive) THEN BEGIN
            PLOT_NUWT_PEAKS, img_data[*,*,k], slitnum=k, grad=grad, $
                             /plot_threads, /screen, /use_temp_common_blocks
            PRINT, ''
            PRINT, 'Current min threads length = '+STRTRIM(min_tlen, 2)
            READ, ask_tlen, PROMPT='Would you like to select a different min thread length? (y/n) '
            ask_tlen = STRLOWCASE(ask_tlen)
            IF ask_tlen.startswith('n') THEN BEGIN
                loop_tlen = 0
            ENDIF ELSE BEGIN
                READ, min_tlen, PROMPT='Please enter a new min thread length to try: '
            ENDELSE
            ;close the figures afterwards
            fig_peaks = GETWINDOWS('fig_peaks')
            fig_threads = GETWINDOWS('fig_threads')
            IF fig_peaks THEN fig_peaks.close
            IF fig_threads THEN fig_threads.close
        ENDIF ELSE BEGIN
            ;end 'loop' if not in interactive mode
            loop_tlen = 0
        ENDELSE
    ENDWHILE

    ;save the values used to the meta structure
    nuwt_meta.min_tlen[k] = min_tlen
    nuwt_meta.max_dist_jump[k] = max_dist_jump
    nuwt_meta.max_time_skip[k] = max_time_skip

    ;Fill in data gaps and then find waves using FFT
    PRINT, '   Filling data gaps ...'
    NUWT_PATCH_UP_THREADS

    IF KEYWORD_SET(bootstrap) THEN BEGIN
        IF KEYWORD_SET(pad_fft) THEN BEGIN
            PRINT, '   Running FFT with both padding and bootstrapping ...'
            NUWT_APPLY_FFT, res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                           /pad_fft, pad_length=pad_length, fill_pad=fill_pad, $
                           /bootstrap, num_bootstrap=num_bootstrap
        ENDIF ELSE BEGIN
            PRINT, '   Running FFT with bootstrapping ...'
            NUWT_APPLY_FFT, res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                           /bootstrap, num_bootstrap=num_bootstrap
        ENDELSE
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(pad_fft) THEN BEGIN
            PRINT, '   Running FFT with padding ...'
            NUWT_APPLY_FFT, res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                           /pad_fft, pad_length=pad_length, fill_pad=fill_pad, $
                           vel_amp_mode=vel_amp_mode
        ENDIF ELSE BEGIN
            PRINT, '   Running FFT ...'
            NUWT_APPLY_FFT, res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                            vel_amp_mode=vel_amp_mode
        ENDELSE
    ENDELSE

    n_waves = TOTAL(fft_peaks.num_saved_waves)
    PRINT, '   Number of waves identified: '+STRTRIM(n_waves, 2)

    ;Setting the auto_qual_flag values
    SET_NUWT_QUAL_FLAGS, /use_temp_common_blocks, $
                         max_percent_gaps=35.0

    ;Calculate the bulk statistics
    CALC_NUWT_BULK_STATS, /use_temp_common_blocks, $
                          min_auto_flag=2, bulk_stats_out=slit_bulk_stats

    ;Append all results to master lists and output structures
    nuwt_located.ADD, located
    nuwt_threads.ADD, threads
    nuwt_fft_spec.ADD, fft_spec
    nuwt_fft_peaks.ADD, fft_peaks

    nuwt_meta.num_threads[k] = n_threads
    nuwt_meta.num_waves[k] = n_waves
    nuwt_bulk_stats.ADD, slit_bulk_stats
ENDFOR

;@save_results##################################################################
; SAVE MASTER LISTS AND STRUCTURES WITH THE RESULTS FOR ALL DATA SLITS
;###############################################################################

PRINT, 'Saving COMMON block data for all slits in:'
PRINT, '   Save folder: '+save_folder
PRINT, '   filename: '+save_filename
SAVE, nuwt_meta, nuwt_bulk_stats, nuwt_located, nuwt_threads, $
      nuwt_fft_spec, nuwt_fft_peaks, FILENAME=save_folder+save_filename


PRINT, '##### Finished running NUWT #####'
PRINT, '   ' ;last empty line printed
END
