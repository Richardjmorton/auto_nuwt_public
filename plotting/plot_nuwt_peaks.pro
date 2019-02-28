;+
;NAME: PLOT_NUWT_PEAKS
;
;PURPOSE:
;   Plots the intensity peaks found by 'locate_things.pro' and the threads found
;   by 'follow_thread.pro'. Will also plot the input td-diagram as well as histograms
;   of the gradients on both the left and right sides (as an aid for selecting
;   a good gradient. Can plot either to the screen or a PDF (for later refference)
;
;INPUTS:
;   input_data - [optional positional input] image that was processed by NUWT
;               (typically a td-diagram). If no image is given, the program will
;               default to a copy of the image stored within the 'located' structures
;               produced by "locate_things.pro" (note: this image may be despiked
;               or have inverted values).
;
;   All other required inputs are loaded from the 'all_nuwt_dat' COMMON block
;   containing the results from 'locate_things.pro' and 'follow_threads.pro'
;
;OPTIONAL INPUTS:
;   res - spatial resolution of data (in [arcsec]). Default is 1
;   cad - temporal cadence of the data (in [s]). Default is 1
;   km_per_arcsec - ratio of [km]/[arcsec]. Defaults to 725.27 which is the
;                   mean scale of the solar surface.
;   grad - gradient used to filter peaks. If not set, will estimate from the data
;   /final_units - [DEPRECIATED] if set, will load "res", "cad", and "_units" values from the
;                  NUWT "located" structure.
;   dist_units - string indicating what units to use for distance in the plot
;                labels. Defaults to units of [pixels] unless "res" is set,
;                in which case the default is [arcsec]
;   time_units - string indicating what units to use for time in the plot
;                labels. Defaults to units of [steps] unless "cad" is set,
;                in which case the default is [s]
;   /underplot_td - if set, will underplot the td-diagram beneath the "Peak locations"
;                   and "Threads found" images. Note: plots may become harder to read.
;   /show_rejected_peaks - if set, will plot all of the peaks rejected by the gradient
;   /normalized - if set, will normalize the gradient histograms by the number of peaks
;   /plot_td_only - if set, will only plot the input td-diagram (useful mainly for
;                   presentations or debugging)
;   /plot_peaks - if set, will plot input td, located peaks, and histograms of gradients
;   /plot_threads - if set, will plot located peaks and found threads
;   /screen - if set, will plot to the screen rather than a file
;       Note: Defaults to plotting everything if neither /plot_peaks or /plot_threads is set
;             For consistancy, will also plot everything if saving to a file.
;   /use_temp_common_blocks - if set, will use the temporary common blocks rather
;                             than the all_nuwt_dat block. Use with care!
;   slitnum - virtual slit number for header text. By defualt assumes 1
;   header - custom header text. Useful for identifying the source data.
;            Defaults to 'NUWT FFT stats'.
;   run_tag - custom string that will be appended to the end of the header.
;             Normally used to keep track of differnt runs of the same dataset.
;             There is no default run_tag string.
;   save_folder - folder in which to save the plots. Defaults to the user's home folder.
;                 Will also append a '/' to the end if not included.
;   filename - name for the output PDF file. Default is 'NUWT_selected_peaks'
;
;OUTPUTS:
;   PDF_file - multi-page PDF containing (1) the input td-diagram, (2) selected
;              peaks, (3) threads found, and (4) histograms of gradients
;
;HISTORY: Name---------Date---------Description
;         M Weberg   MAY, 2016  Initial coding. Was part of a modified version
;                               of 'locate_things.pro'
;         M Weberg  JULY, 2016  Moved to seperate program and extended options
;         M Weberg   AUG, 2016  Fixed gradient estimation & added basic unit conversions
;         M Weberg  SEPT, 2016  Large update:
;                                - more plots (better PDF output, can also show threads)
;                                - more diagnostic information
;                                - more options (for interatively exploring parameters)
;                                - reconfigured to use NUWT master COMMON block
;         M Weberg   JAN, 2017  Added the "underplot_td" option
;         M Weberg   MAR, 2017  Now can load unit information from the NUWT COMMON
;                               block structures
;
;TO-DO / LIMITATIONS:
;   - Add inputs for metadata indicating the data and source of observations.
;-
FUNCTION FIND_PERCENTILE, input_vals, percentile, nonzero=nonzero, invert=invert
    IF KEYWORD_SET(invert) THEN percentile = 100 - percentile
    IF KEYWORD_SET(nonzero) THEN BEGIN
        IF NOT KEYWORD_SET(invert) THEN BEGIN
            input_vals = input_vals[where(input_vals GT 0.0)]
        ENDIF ELSE BEGIN
            input_vals = input_vals[where(input_vals LT 0.0)]
        ENDELSE
    ENDIF
    len_data = double(N_ELEMENTS(input_vals))
    sorted_data = input_vals[SORT(input_vals)]
    percent_index = FIX(len_data*(percentile/100.0))
    output_val = sorted_data[percent_index]
RETURN, output_val
END

PRO PLOT_NUWT_PEAKS, input_data, res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                     grad=grad, final_units=final_units, $
                     dist_units=dist_units, time_units=time_units, $
                     underplot_td=underplot_td, show_rejected_peaks=show_rejected_peaks, $
                     normalized=normalized,  $

                     plot_td_only=plot_td_only, plot_peaks=plot_peaks, $
                     plot_threads=plot_threads, screen=screen, $

                     use_temp_common_blocks=use_temp_common_blocks, $
                     slitnum=slitnum, header=header, run_tag=run_tag, $
                     save_folder=save_folder, filename=filename

;@defaults######################################################################
;Seting default values
;###############################################################################
IF N_ELEMENTS(slitnum) EQ 0 THEN slitnum = 0
IF NOT KEYWORD_SET(header) THEN header = 'NUWT Peaks'
IF KEYWORD_SET(run_tag) THEN header = header + ' - '+run_tag
IF NOT KEYWORD_SET(filename) THEN filename = 'NUWT_peaks'
IF NOT KEYWORD_SET(save_folder) THEN CD, CURRENT=save_folder ;i.e. defaults to current folder
IF strlen(save_folder) GT 0 AND NOT save_folder.endswith('/') THEN save_folder = save_folder+'/'
IF NOT KEYWORD_SET(plot_peaks) AND NOT KEYWORD_SET(plot_threads) THEN BEGIN
    ;Defaults to plotting everything
    plot_td = 1
    plot_peaks = 1
    plot_threads = 1
ENDIF
IF NOT KEYWORD_SET(screen) THEN BEGIN
    ;Also plots everything if saving to file
    plot_td = 1
    plot_peaks = 1
    plot_threads = 1
ENDIF
IF KEYWORD_SET(plot_td_only) THEN BEGIN
    plot_td = 1
    plot_peaks = 0
    plot_threads = 0
ENDIF

;Colors to help differentiate threads
th_color_list = ['white', 'dark grey', 'light coral', 'crimson', 'gold', 'light green', $
                 'forest green', 'aqua', 'dodger blue', 'violet', 'dark violet']

;@load_common###################################################################
;LOADING DATA AND SELECTING CORRECT SET OF RESULTS
;###############################################################################
IF NOT KEYWORD_SET(use_temp_common_blocks) THEN BEGIN
    ;Default and preferred data source
    COMMON all_nuwt_dat, nuwt_located, nuwt_threads, nuwt_fft_spec, nuwt_fft_peaks

    IF slitnum GT N_ELEMENTS(nuwt_located) THEN BEGIN
        last_slit = N_ELEMENTS(nuwt_located)-1
        print, 'WARNING! slit number '+strtrim(slitnum, 2)+' does not exist!'
        print, '   slitnum set to '+strtrim(last_slit, 2)+' (last set of results)'
        slitnum = last_slit
    ENDIF

    slit_located = nuwt_located[slitnum]
    slit_threads = nuwt_threads[slitnum]
ENDIF ELSE BEGIN
    ;Used for running NUWT in interactive mode.
    ;WARNING: NOT recommended unless you know what you are doing
    COMMON located_dat, located
    COMMON threads_dat, threads

    slit_located = located
    slit_threads = threads
ENDELSE

n_threads = N_ELEMENTS(slit_threads)
num_th_colors = N_ELEMENTS(th_color_list)

;catch the case of zero threads found
IF n_threads EQ 1 THEN BEGIN
    only_tlen = threads[0].length
    IF only_tlen LE 0 THEN n_threads = 0
ENDIF

;@check_data####################################################################
; CHECK THE INPUT DATA AND EXTRACT VARIABLES FROM SLIT STRUCTURES
;###############################################################################
IF N_ELEMENTS(input_data) EQ 0 THEN BEGIN
    ;Set default image to plot if the user does not provide one
    IF TOTAL(STRMATCH(tag_names(slit_located), 'TD_IMG', /FOLD_CASE)) GE 1 THEN BEGIN
        ;This image may be despiked or have inverted values
        img_data = slit_located.td_img[*,*]
    ENDIF ELSE BEGIN
        ;Needed for backwards compatibility with older NUWT runs (before Dec. 2017)
        print, 'WARNING: incomplete intensity values! Plots show peak intensities only.'
        img_data = slit_located.peaks[*,*,1] ;intensity values at peaks ONLY
    ENDELSE
ENDIF ELSE BEGIN
    IF N_TAGS(input_data) GT 0 THEN BEGIN
        ;input_data is a structure
        IF TOTAL(STRCMP(TAg_NAMES(input_data), 'slit', /fold_case)) GE 1 THEN BEGIN
            img_data = input_data.slit
        ENDIF ELSE BEGIN
            PRINT, 'ERROR: invalid structure format!'
            MESSAGE, 'Input data strutures must have both "slit" and "meta" tags'
        ENDELSE
    ENDIF ELSE BEGIN
        img_data = input_data
    ENDELSE
ENDELSE

sz = size(img_data)
nx = sz[1]
nt = sz[2]
IF sz[0] EQ 3 THEN BEGIN
    ;select the correct td-diagram for the requested slitnum
    slit_img = img_data[*,*,slitnum]
ENDIF ELSE BEGIN
    ;[DEFAULT] Only one 2D image was given
    slit_img = img_data[*,*,0]
ENDELSE

;@units#########################################################################
;UNIT CONVERSIONS AND LABELS (AUTOMATED METHODS)
;###############################################################################
;Doing it this way is less wordy and prevents leakage of variables
IF KEYWORD_SET(dist_units) THEN plot_dist_units = dist_units
IF KEYWORD_SET(time_units) THEN plot_time_units = time_units

IF KEYWORD_SET(res) THEN BEGIN
    IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27
    IF NOT KEYWORD_SET(plot_dist_units) THEN plot_dist_units = 'arcsec'
ENDIF ELSE BEGIN
    ;[DEFAULT] Load dist unit information automatically from the NUWT output
    IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = slit_located.km_per_arcsec
    res = slit_located.res
    dx = slit_located.dx
    nuwt_dx_units = slit_located.units_dx
    IF KEYWORD_SET(plot_dist_units) AND nuwt_dx_units EQ 'pixels' THEN BEGIN
        print, 'WARNING: "res" is unknown! Setting plot_dist_units = "pixels"'
        plot_dist_units = 'pixels'
    ENDIF ELSE IF NOT KEYWORD_SET(plot_dist_units) THEN BEGIN
        plot_dist_units = nuwt_dx_units
    ENDIF
ENDELSE

IF KEYWORD_SET(cad) THEN BEGIN
    IF NOT KEYWORD_SET(plot_time_units) THEN plot_time_units = 'min'
ENDIF ELSE BEGIN
    ;[DEFAULT] Load time unit information automatically from the NUWT output
    cad = slit_located.cad
    dt = slit_located.dt
    nuwt_dt_units = slit_located.units_dt
    IF KEYWORD_SET(plot_time_units) AND nuwt_dt_units EQ 'timesteps' THEN BEGIN
        print, 'WARNING: "cad" is unknown! Setting plot_time_units = "timesteps"'
        plot_time_units = 'timesteps'
    ENDIF ELSE IF NOT KEYWORD_SET(plot_time_units) THEN BEGIN
        plot_time_units = nuwt_dt_units
    ENDIF
ENDELSE

;Unit conversions
IF plot_dist_units EQ 'pixels' THEN dx = 1.0
IF plot_dist_units EQ 'arcsec' THEN dx = res
IF plot_dist_units EQ 'm' THEN dx = res*km_per_arcsec*1000.0
IF plot_dist_units EQ 'km' THEN dx = res*km_per_arcsec
IF plot_dist_units EQ 'Mm' THEN dx = res*km_per_arcsec/1000.0

IF plot_time_units EQ 'timesteps' THEN dt = 1.0
IF plot_time_units EQ 's' THEN dt = cad
IF plot_time_units EQ 'min' THEN dt = cad/60.0
IF plot_time_units EQ 'hr' THEN dt = cad/3600.0

IF KEYWORD_SET(final_units) THEN BEGIN
    ;remove once code is fully refactored
    print, "Notice: the /final_units keyword option is unnecessary and has been depreciated!"
ENDIF

;@calc_stats####################################################################
;FILTERING AND BASIC STATISTICS
;###############################################################################
; Finding the locations of various types of peaks
all_pk_locs = where(slit_located.allpeaks[*,*] GT 0, /NULL)
valid_pk_locs = where(slit_located.allpeaks[*,*]  GE 2, /NULL) ;includes both whole pixel (default) and sub-pixel (gaussin fit) peaks
rejected_pk_locs = where(slit_located.allpeaks[*,*] EQ 1, /NULL)
failed_gauss_locs = where(slit_located.allpeaks[*,*] EQ 2, /NULL)
gauss_fit_locs = where(slit_located.allpeaks[*,*] EQ 3, /NULL)

IF N_ELEMENTS(grad) EQ 0 THEN BEGIN
    ;tries to look at the results and infer what the cutoff gradient was
    abs_min_valid_grad_left = 0.0
    abs_min_valid_grad_right = 0.0
    IF N_ELEMENTS(all_pk_locs) NE N_ELEMENTS(valid_pk_locs) THEN BEGIN
        nonzero_valid_gl_locs = where(slit_located.peaks[*,*,1] NE 0 AND $
                                      slit_located.allpeaks[*,*] GE 2 AND $
                                      slit_located.grad_left GT 0, /NULL)
        nonzero_valid_gr_locs = where(slit_located.peaks[*,*,1] NE 0 AND $
                                      slit_located.allpeaks[*,*] GE 2 AND $
                                      slit_located.grad_right LT 0, /NULL)
        IF N_ELEMENTS(nonzero_valid_gl_locs) GT 0 AND N_ELEMENTS(nonzero_valid_gl_locs) GT 0 THEN BEGIN
            abs_min_valid_grad_left = min(abs(slit_located.grad_left[nonzero_valid_gl_locs]))
            abs_min_valid_grad_right = min(abs(slit_located.grad_right[nonzero_valid_gr_locs]))
        ENDIF
    ENDIF
    grad = mean([abs_min_valid_grad_left, abs_min_valid_grad_right])
    tag_grad_calc = 'estimated'
ENDIF ELSE tag_grad_calc = 'selected'
; IF grad GE 0 THEN print, tag_grad_calc+' gradient cutoff = ', grad

; Histogramming the gradients
gl_max = FIND_PERCENTILE(slit_located.grad_left[all_pk_locs], 95, /nonzero)
gr_min = FIND_PERCENTILE(slit_located.grad_right[all_pk_locs], 95, /nonzero, /invert)

abs_max_filtered_grad = MAX([ABS(gl_max), ABS(gr_min)])
; gl_hist = histogram(slit_located.grad_left[all_pk_locs], nbins=100, min=1e-20, max=gl_max, LOCATIONS=gl_bins)
gl_hist = histogram(slit_located.grad_left[all_pk_locs], nbins=100, min=1e-20, max=abs_max_filtered_grad, LOCATIONS=gl_bins)

;note: the inverting is needed to fix some buggy IDL histogramming
; gr_hist = histogram(-slit_located.grad_right[all_pk_locs], nbins=100, min=1e-20, max=-gr_min, LOCATIONS=gr_bins)
gr_hist = histogram(-slit_located.grad_right[all_pk_locs], nbins=100, min=1e-20, max=abs_max_filtered_grad, LOCATIONS=gr_bins)
gr_bins = -gr_bins

; Calculate some basic statistics
num_all_peaks = N_ELEMENTS(all_pk_locs)
num_rejected_peaks = N_ELEMENTS(rejected_pk_locs)
num_selected_peaks = N_ELEMENTS(valid_pk_locs)
num_sub_pixel_peaks = N_ELEMENTS(gauss_fit_locs)
num_whole_pixel_peaks = N_ELEMENTS(failed_gauss_locs)

percent_rejected = 100.0 - 100.0*FLOAT(num_selected_peaks)/FLOAT(num_all_peaks)

IF KEYWORD_SET(normalized) AND num_all_peaks GT 0 THEN BEGIN
    gl_hist = float(gl_hist)/(float(num_all_peaks))
    gr_hist = float(gr_hist)/(float(num_all_peaks))
    grad_y_axis_label = 'Fraction of peaks'
ENDIF ELSE grad_y_axis_label = 'Number of peaks'

;@plot##########################################################################
;PLOTTING THE DATA
;###############################################################################
IF KEYWORD_SET(plot_td) OR KEYWORD_SET(plot_peaks) THEN BEGIN
    ; ##### PAGE 1 (or 1st plot) - INPUT TD-DIAGRAM #####
    IF KEYWORD_SET(screen) THEN BEGIN
        fig_td = window(name='fig_td', dimensions=[600, 720])
    ENDIF ELSE BEGIN
        fig_td = window(name='fig_td', dimensions=[600, 720], /buffer)
    ENDELSE
    txt_slit_ID = text(0.01, 0.97, 'Slit # '+strtrim(slitnum, 2), font_size=14, font_style=1)
    txt_im_size = text(0.05, 0.94, strtrim(nx, 2)+' x '+strtrim(nt, 2)+' pixels', font_size=12)
    txt_header = text(0.20, 0.97, header, font_size=14, font_style=1)
    txt_timestamp = text(0.76, 0.01, 'Date & Time Created$\n$'+systime(), font_size=9, color='dark grey')

    text_im_title = text(0.485, 0.94, 'Image Data', font_size=12, font_style=1)
    plt_im = IMAGE(slit_img, axis_style=0, position=[0.125, 0.10, 0.995, 0.92], /current)
    im_x_ax = AXIS('X', location=0.0, target=plt_im, coord_transform=[0.0, dx], tickdir=1, $
                   title='Distance ['+plot_dist_units+']', ticklen=0.02)
    im_t_ax = AXIS('Y', location=0.0, target=plt_im, coord_transform=[0.0, dt], tickdir=1, $
                   title='Time ['+plot_time_units+']', ticklen=0.02)
ENDIF

IF KEYWORD_SET(plot_peaks) OR KEYWORD_SET(plot_threads) THEN BEGIN
    ; ##### PAGE 2 (or 2nd plot) - DATA AND PEAK LOCS #####
    IF KEYWORD_SET(screen) THEN BEGIN
       fig_peaks = window(name='fig_peaks', dimensions=[600, 720])
    ENDIF ELSE BEGIN
       fig_peaks = window(name='fig_peaks', dimensions=[600, 720], /buffer)
    ENDELSE
    txt_slit_ID = text(0.01, 0.97, 'Slit # '+strtrim(slitnum, 2), font_size=14, font_style=1)
    txt_num_selected = text(0.05, 0.94, strtrim(num_selected_peaks, 2)+' selected peaks', font_size=12)
    txt_grad_cutoff = text(0.71, 0.94, '~gradient: '+strtrim(grad, 2), font_size=12)
    txt_header = text(0.20, 0.97, header, font_size=14, font_style=1)
    txt_timestamp = text(0.76, 0.01, 'Date & Time Created$\n$'+systime(), font_size=9, color='dark grey')

    txt_pk_title  = text(0.465, 0.94, 'Peak Locations', font_size=12, font_style=1)
    IF KEYWORD_SET(underplot_td) THEN BEGIN
        plt_pk = IMAGE(slit_img, axis_style=0, position=[0.125, 0.10, 0.995, 0.92], /current)
    ENDIF ELSE BEGIN
        plt_pk = IMAGE(fltarr(nx,nt), axis_style=0, position=[0.125, 0.10, 0.995, 0.92], /current)
    ENDELSE

    ;Note: the code below plots gaussin fit peaks in white (255, 255, 255) and
    ;      peaks defaulting to nearest whole pixels in green (0, 255, 0)
    valid_pk_image = fltarr(nx,nt,4) ;RGBA image (i.e. has an alpha channel)
    temp_peaks = fltarr(nx,nt)
    ; sub-pixel, gaussian fit peaks only
    temp_peaks[gauss_fit_locs] = 255
    valid_pk_image[*,*,0] = temp_peaks
    valid_pk_image[*,*,2] = temp_peaks
    ; both nearest pixel AND gaussian fit peaks
    temp_peaks[valid_pk_locs] = 255
    valid_pk_image[*,*,1] = temp_peaks
    valid_pk_image[*,*,3] = temp_peaks ;full opacity for ALL valid peaks
    plt_pk = IMAGE(valid_pk_image, axis_style=0, overplot=plt_pk, background_transparency=100)
    fail_gauss_fit_tag = text(0.05, 0.03, 'peaks defaulted to nearest whole pixel', color='green', font_size=10)

    IF KEYWORD_SET(show_rejected_peaks) THEN BEGIN
        ;plot invalid peaks in red (255, 0, 0)
        invalid_pk_with_alpha = fltarr(nx,nt,4) ;RGBA image (i.e. has an alpha channel)
        temp_peaks = fltarr(nx,nt)
        temp_peaks[rejected_pk_locs] = 255
        invalid_pk_with_alpha[*,*,0] = temp_peaks ;full red
        invalid_pk_with_alpha[*,*,3] = temp_peaks ;full opacity
        plt_other_pk = IMAGE(invalid_pk_with_alpha, axis_style=0, overplot=plt_pk, background_transparency=100)
        rejected_tag = text(0.05, 0.01, 'peaks rejected by gradient', color='red', font_size=10)
    ENDIF

    pk_x_ax = AXIS('X', location=0.0, target=plt_pk, coord_transform=[0.0, dx], tickdir=1, $
                   title='Distance ['+plot_dist_units+']', ticklen=0.02)
    pk_t_ax = AXIS('Y', location=0.0, target=plt_pk, coord_transform=[0.0, dt], tickdir=1, $
                   title='Time ['+plot_time_units+']', ticklen=0.02)
ENDIF

IF KEYWORD_SET(plot_threads) THEN BEGIN
    ; ##### PAGE 3 (or 3rd plot) - found threads #####
    IF KEYWORD_SET(screen) THEN BEGIN
        fig_threads = window(name='fig_threads', dimensions=[600, 720])
    ENDIF ELSE BEGIN
        fig_threads = window(name='fig_threads', dimensions=[600, 720], /buffer)
    ENDELSE
    txt_slit_ID = text(0.01, 0.97, 'Slit # '+strtrim(slitnum, 2), font_size=14, font_style=1)
    txt_num_threads = text(0.05, 0.94, strtrim(n_threads, 2)+' threads', font_size=12)
    txt_header = text(0.20, 0.97, header, font_size=14, font_style=1)
    txt_timestamp = text(0.76, 0.01, 'Date & Time Created$\n$'+systime(), font_size=9, color='dark grey')
    th_color_tag = text(0.05, 0.02, 'colors indicate different threads', color='black', font_size=10)

    txt_th_title = text(0.465, 0.94, 'Threads Found', font_size=12, font_style=1)
    IF KEYWORD_SET(underplot_td) THEN BEGIN
        plt_th_im = IMAGE(slit_img, axis_style=0, position=[0.125, 0.10, 0.995, 0.92], /current)
    ENDIF ELSE BEGIN
        plt_th_im = IMAGE(fltarr(nx,nt), axis_style=0, position=[0.125, 0.10, 0.995, 0.92], /current)
    ENDELSE

    IF n_threads GT 0 THEN BEGIN
        all_tsteps = findgen(nt)
        min_th_len = 1e10 ; initialize with unreasonably large number
        FOR h=0L, n_threads-1 DO BEGIN
            c_num = FIX(h MOD num_th_colors)
            tpos = slit_threads[h].pos
            tpos = tpos[slit_threads[h].start_bin:slit_threads[h].end_bin]
            tstep = all_tsteps[slit_threads[h].start_bin:slit_threads[h].end_bin]
            th_len = slit_threads[h].length
            IF th_len GT 0 THEN BEGIN
                plt_th = PLOT(tpos, tstep, linestyle='none', symbol='.', color=th_color_list[c_num], overplot=plt_th_im)
            ENDIF
            IF th_len LT min_th_len THEN min_th_len = th_len
        ENDFOR
        txt_min_tlen = text(0.74, 0.94, 'min length: '+strtrim(min_th_len, 2), font_size=12)
    ENDIF ELSE BEGIN
        txt_no_threads = text(0.45, 0.5, 'No threads found!', font_size=12, color='white')
    ENDELSE

    im_x_ax = AXIS('X', location=0.0, target=plt_th_im, coord_transform=[0.0, dx], tickdir=1, $
                   title='Distance ['+plot_dist_units+']', ticklen=0.02)
    im_t_ax = AXIS('Y', location=0.0, target=plt_th_im, coord_transform=[0.0, dt], tickdir=1, $
                   title='Time ['+plot_time_units+']', ticklen=0.02)
ENDIF

IF KEYWORD_SET(plot_peaks) THEN BEGIN
    ; ##### PAGE 4 (or 4th plot) - GRADIENTS #####
    IF KEYWORD_SET(screen) THEN BEGIN
        fig_grads = window(name='fig_grads', dimensions=[600, 720])
    ENDIF ELSE BEGIN
        fig_grads = window(name='fig_grads', dimensions=[600, 720], /buffer)
    ENDELSE

    ; Left-hand side gradients (TOP SUBPLOT)
    ; plt_left = barplot(gl_bins, gl_hist, histogram=1, xrange=[0.0, gl_max], $abs_max_filtered_grad
    plt_left = barplot(gl_bins, gl_hist, histogram=1, xrange=[0.0, abs_max_filtered_grad], $
                      color='black', fill_color='light grey', linestyle='-', $
                      position=[0.12, 0.51, 0.92, 0.84], /current)
    plt_left.title = 'Left-hand side peak gradients'
    plt_pos_grad_line = plot([grad, grad], [0.0, plt_left.yrange[1]], symbol='none', linestyle='-', color='red', overplot=plt_left)
    plt_left.xtitle = 'Intensity Gradient'
    plt_left.ytitle = grad_y_axis_label

    ; Right-hand side gradients (BOTTOM SUBPLOT)
    ; plt_right = barplot(gr_bins, gr_hist, histogram=1, xrange=[gr_min, 0.0], $
    plt_right = barplot(gr_bins, gr_hist, histogram=1, xrange=[-abs_max_filtered_grad, 0.0], $
                      color='black', fill_color='light grey', linestyle='-', $
                      position=[0.12, 0.08, 0.92, 0.41], /current)
    plt_right.title = 'Right-hand side peak gradients'
    plt_neg_grad_line = plot([-grad, -grad], [0.0, plt_right.yrange[1]], symbol='none', linestyle='-', color='red', overplot=plt_right)
    plt_right.xtitle = 'Intensity Gradient'
    plt_right.ytitle = grad_y_axis_label

    ; ##### PRINTED TEXT (on page 4)#####
    txt_slit_ID = text(0.01, 0.97, 'Slit # '+strtrim(slitnum, 2), font_size=14, font_style=1)
    txt_header = text(0.20, 0.97, header, font_size=14, font_style=1)

    txt_total_num_peaks = text(0.05, 0.94, strtrim(num_all_peaks, 2)+' total local peaks', font_size=12)
    txt_num_rejected = text(0.05, 0.915, strtrim(num_rejected_peaks, 2)+' rejected by gradient'+ $
                                         ' ('+STRTRIM(string(percent_rejected, format='(f6.2)'), 2)+'%)', font_size=12)
    txt_num_selected = text(0.05, 0.89, strtrim(num_selected_peaks, 2)+' selected ('+ $
                            strtrim(num_sub_pixel_peaks, 2)+' sub-pixel accuracy / '+$
                            strtrim(num_whole_pixel_peaks, 2)+' whole pixel)', font_size=12)
    txt_grad_cutoff = text(0.48, 0.94, tag_grad_calc+' gradient cutoff: '+strtrim(grad, 2), font_size=12)

    txt_timestamp = text(0.76, 0.01, 'Date & Time Created$\n$'+systime(), font_size=9, color='dark grey')
ENDIF

; Saving the plot to PDF (default)
IF NOT KEYWORD_SET(screen) THEN BEGIN

    IF KEYWORD_SET(plot_td) THEN fig_td.save, save_folder+filename+'.pdf', page_size=[6.0, 7.2], /append
    IF KEYWORD_SET(plot_peaks) THEN fig_peaks.save, save_folder+filename+'.pdf', page_size=[6.0, 7.2], /append
    IF KEYWORD_SET(plot_threads) THEN fig_threads.save, save_folder+filename+'.pdf', page_size=[6.0, 7.2], /append
    IF KEYWORD_SET(plot_peaks) THEN fig_grads.save, save_folder+filename+'.pdf', page_size=[6.0, 7.2], /append, /close
    print, 'Finished plotting NUWT located peaks!'
    print, 'Save Folder: ', save_folder
    print, 'Filename: ', filename+'.pdf'
ENDIF

END
