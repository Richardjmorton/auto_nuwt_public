;+
;NAME: MERGE_NUWT_RESULTS
;
;PURPOSE:
;   Combine results from seperate NUWT runs into single set of result lists.
;   This intended to make plotting the combined results easier. Care should be
;   when merging data since most plotting routines assume that all data slits
;   within a single "run" share the same spatial resolution and sampling cadence.
;   Note: not all metadata will be properly copied over. In particular, the NUWT
;   run parameters are assumed to be the same for all slits.
;
;INPUTS:
;   None required. By default, will attept to merge all files in the the current
;   directory that contain "nuwt_results" in their name and end in ".sav".
;
;OPTIONAL INPUTS:
;   load_folder - folder to search of valid "nuwt_results" files. Defaults to the
;                 current working directory. Will automatically append a '/' to
;                 the end if not included.
;   timestamps - string or array of strings containing the timestamps of the
;                input files. May be used by itself or in conjuntion with the
;                "labels" varible. If set, the timestamps are assumed to be at
;                the BEGINNING of the filenames with nothing before them.
;   labels - string or array of strings containing the labels of the input files.
;            May be used by itself or in conjuntion with the "timestamps" varible.
;            If set, labels are assumed to directly DIRECTLY PRECEDE "nuwt_results.sav"
;            in the filenames. If both "timestamps" and "labels" are set, then
;            the filenames are constructed as "TIMESTAMP_LABEL_nuwt_results.sav"
;   save_folder - folder in which to save the results. Defaults to current folder.
;                 Will automatically append a '/' to the end if not included.
;   filetag - unique tag to append to the START of the output file
;   default_plot_header - string that can be optioanlly defined and used when
;                         plotting the data using the NUWT plot procedures.
;                         Defaults to "Merged NUWT Results".
;
;OUTPUTS:
;   An IDL .sav file named "merged_nuwt_results.sav" contains the merged lists
;   Also loads the lists into the "nuwt_meta_dat" & "all_nuwt_dat" COMMON blocks.
;
;HISTORY: Name---------Date---------Description
;         M Weberg  24 JULY, 2017  Initial coding.
;
;TO-DO / LIMITATIONS:
;   - Assumes that all data slits share the same resolution and cadence and where
;     processed by NUWT using the same options (some parameters are allowed to be different)
;-


PRO MERGE_NUWT_RESULTS, load_folder=load_folder, timestamps=timestamps, labels=labels, $
                        save_folder=save_folder, filetag=filetag, default_plot_header=default_plot_header

COMPILE_OPT IDL2
;###############################################################################
;SETTING DEFAULT VALUES AND VALIDATING INPUTS
;###############################################################################
COMMON all_nuwt_dat, nuwt_located, nuwt_threads, nuwt_fft_results, nuwt_fft_peaks
COMMON nuwt_meta_dat, nuwt_meta, nuwt_bulk_stats

IF NOT KEYWORD_SET(load_folder) THEN CD, CURRENT=load_folder ;i.e. defaults to current folder
IF NOT load_folder.endswith('/') THEN load_folder = load_folder+'/'
IF NOT KEYWORD_SET(save_folder) THEN CD, CURRENT=save_folder ;i.e. defaults to current folder
IF NOT save_folder.endswith('/') THEN save_folder = save_folder+'/'
IF NOT KEYWORD_SET(filetag) THEN filetag='' ELSE filetag = filetag+'_'
save_filename = filetag+'merged_nuwt_results.sav'

IF NOT KEYWORD_SET(default_plot_header) THEN default_plot_header = 'Merged NUWT Results'
;###############################################################################
;ASSEMBLING FILE PATHS (OR SEARCHING FOR FILES IN THE GIVEN FOLDER)
;###############################################################################
IF KEYWORD_SET(timestamps) OR KEYWORD_SET(labels) THEN BEGIN
    ;Contruct the filenames using input values. Good for when there are other
    ;"nuwt_results" files in the same folder that you DO NOT want to merge
    separator = ''
    IF KEYWORD_SET(timestamps) AND KEYWORD_SET(labels) THEN separator = '_'
    IF NOT KEYWORD_SET(timestamps) THEN timestamps = ''
    IF NOT KEYWORD_SET(labels) THEN labels = ''

    num_timestamps = N_ELEMENTS(timestamps)
    num_labels = N_ELEMENTS(labels)

    ;ensure that single values are NOT contained in a single element array
    IF num_timestamps EQ 1 THEN timestamps = timestamps[0]
    IF num_labels EQ 1 THEN labels = labels[0]

    IF (num_timestamps NE num_labels) AND $
       (num_timestamps GT 1 AND num_labels GT 1)THEN BEGIN
        PRINT, 'ERROR: number of input timestamps and labels are incompatible!'
        MESSAGE, 'Please input either single values or arrays with the same number of elements'
    ENDIF

    nuwt_files = load_folder+timestamps+separator+labels+'_nuwt_results.sav'

ENDIF ELSE BEGIN
    ;Search "load_folder" for all .sav files containing "nuwt_results" in the filename
    nuwt_files = findfile(load_folder+'*'+'nuwt_results'+'*'+'.sav')
ENDELSE

num_files = N_ELEMENTS(nuwt_files)
IF num_files GT 1 THEN BEGIN
    PRINT, ' '
    PRINT, 'Found '+STRTRIM(num_files, 2)+' "nuwt_results" files in '+load_folder
ENDIF ELSE IF num_files GT 1 THEN BEGIN
    MESSAGE, 'NOTICE: Only one "nuwt_results" file found! No merging is needed.'
ENDIF ELSE BEGIN
    MESSAGE, 'ERROR: No "nuwt_results" files found!'
ENDELSE

;###############################################################################
;MERGING DATA TOGETHER INTO A "SINGLE" NUWT RUN
;###############################################################################
temp_nuwt_meta = LIST()
temp_slit_meta = LIST()
total_num_merged_slits = 0

merged_located = LIST()
merged_threads = LIST()
merged_fft_results = LIST()
merged_fft_peaks = LIST()
merged_bulk_stats = LIST()

PRINT, ' '
PRINT, 'Merging results from the following files:'
FOR f=0, (num_files-1) DO BEGIN
    PRINT, '   '+nuwt_files[f]

    RESTORE, nuwt_files[f]

    num_slits_in_file = nuwt_meta.num_slits
    FOR s=0, (num_slits_in_file-1) DO BEGIN
        merged_located.ADD, nuwt_located[s]
        merged_threads.ADD, nuwt_threads[s]
        merged_fft_results.ADD, nuwt_fft_results[s]
        merged_fft_peaks.ADD, nuwt_fft_peaks[s]
        merged_bulk_stats.ADD, nuwt_bulk_stats[s]
    ENDFOR

    total_num_merged_slits = total_num_merged_slits + num_slits_in_file
    temp_nuwt_meta.ADD, nuwt_meta
    temp_slit_meta.ADD, nuwt_meta.slit_meta
ENDFOR

;###############################################################################
;COMBINING NUWT_META DATA
;###############################################################################
nuwt_meta = {run_date:strarr(total_num_merged_slits), merged_date:systime(), $
             res:nuwt_meta.res, cad:nuwt_meta.cad, km_per_arcsec:nuwt_meta.km_per_arcsec, $
             num_slits:total_num_merged_slits, $
             grad:fltarr(total_num_merged_slits), min_tlen:intarr(total_num_merged_slits), $
             max_dist_jump:intarr(total_num_merged_slits), max_time_skip:intarr(total_num_merged_slits), $
             num_threads:intarr(total_num_merged_slits), num_waves:intarr(total_num_merged_slits), $
             maxima_fitting:'please see the "nuwt_meta" structure for each run', $
             pad_fft:nuwt_meta.pad_fft, fill_pad:nuwt_meta.fill_pad, pad_length:nuwt_meta.pad_length, $
             bootstrap:nuwt_meta.bootstrap, num_bootstrap:nuwt_meta.num_bootstrap, $
             vel_amp_mode:nuwt_meta.vel_amp_mode, $
             default_plot_header:default_plot_header, $
             slit_meta:temp_slit_meta}

FOR f=0, (num_files-1) DO BEGIN
    temp = temp_nuwt_meta[f]
    nuwt_meta.run_date[f] = temp.run_date
    nuwt_meta.grad[f] = temp.grad
    nuwt_meta.min_tlen[f] = temp.min_tlen
    nuwt_meta.max_dist_jump[f] = -1 ;temp.max_dist_jump
    nuwt_meta.max_time_skip[f] = -1 ;temp.max_time_skip
    nuwt_meta.num_threads[f] = temp.num_threads
    nuwt_meta.num_waves[f] = temp.num_waves
ENDFOR

;###############################################################################
;TRANSFERING VALUES TO VARIABLE NAMES EXPECTED BY OTHER NUWT PROGRAMS
;###############################################################################
nuwt_located = merged_located
nuwt_threads = merged_threads
nuwt_fft_results = merged_fft_results
nuwt_fft_peaks = merged_fft_peaks
nuwt_bulk_stats = merged_bulk_stats

SAVE, nuwt_meta, nuwt_bulk_stats, nuwt_located, nuwt_threads, nuwt_fft_results, nuwt_fft_peaks, FILENAME=save_folder+save_filename

PRINT, 'Finished merging NUWT results!'
PRINT, 'Save filepath: '+save_folder+save_filename
PRINT, 'Results have also been preloaded into the "nuwt_meta_dat", & "all_nuwt_dat" COMMON blocks'

END
