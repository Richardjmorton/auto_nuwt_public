;+
;NAME: PREPROCESS_AIA_DATA_FOR_NUWT
;
;PURPOSE:
;   Loads in a set of SDO / AIA .fits files, selects images within a given time
;   range (using the header timestamps), and applies preprocesssing such as
;   alignment, unsharp masking, and smoothing.
;
;INPUTS:
;   start_date - string containing the start time and date. Recommended format:
;                'YYYY-MM-DD hh:mm'
;   end_date - string containing the end time and date.
;
;OPTIONAL INPUTS:
;   wavelength - integer value of a valid AIA wavelength (in units of [Angstroms])
;   /cutout - if set, will extract and process only a cutout window of data
;   x0, y0 - lower left pixel coords of the cutout window. Used only if "/cutout"
;            is set. Default values are x0 = 1025 &  y0 = 0
;   xsize, ysize - pixel width and height of the cutout window. Used only if
;                  "/cutout" is set. Default values are xsize = 2048 & ysize = 1024
;                  When conbined with the defauts for x0 & y0, this extracts a
;                  high-res, 1/8 size image centered on the solar south pole
;   /unsharp - if set, will apply an unsharp filter to the data in the x & y dimensions
;   us_width - pixel width of the unsharp mask in both the x & y dimensions.
;              Default is 10 pixels (~6"x6" for AIA).
;   /smooth - if set, will apply boxcar smoothing along the time dimensions
;   sm_width - number of timesteps to apply the baxcar smoothing over. Default
;              is 3 timesteps (i.e. image frames).
;   filetag - string that will be inserted in the output filename to identify
;             the processed image region. Default is "full_disk" or "cutout" (if
;             /cutout is set). Other common tags include "north", "east", & the like.
;   base_data_folder - local top level folder where ALL input data is stored.
;   full_load_path - full path to the folder containing the requested data.
;                    If not set, will default to the following folder:
;                       base_data_folder+'SDO/AIA/YYYY/MM/DD/WAVELENGTH/level1/'
;   full_save_path - full path where the output files should be saved.
;                    If not set, will default to the following folder:
;                       base_data_folder+'SDO/processed/'
;   input_prefix - file prefix string shared in common by all of the input files.
;                  Used to prefilter for only valid fits files. Default is "aia".
;   fits_ext - file extension of the input fits files. Default is ".fits" which
;              appears to be the most commonly used extension. Less common
;              extensions include .fit and .fts
;   /save_lv_1_5 - if set, will save the processed lv 1.5 data (produced by AIA_PREP.pro)
;                  By default, the lv 1.5 data is NOT saved unless neither /unsharp
;                  or /smooth are set.
;
;OUTPUTS:
;   aia_index - structure containing the header information for the processed data
;   sm_aia_data - smoothed data array (or the string "none" if /smooth is not set)
;   us_aia_data - unsharp masked data array (or the string "none" if /unsharp is not set)
;
;   The above structure and arrays are all saved within the same .sav file:
;   full_save_path+YYYYMMDD_hhmm-hhmm_aia_WAVELENGTH_FILETAG_sm_and_us.sav'
;
;   If /save_lv_1_5 is set (or neither /smooth or /unsharp are set), the lv 1.5
;   data is stored in the array "aia_data" and it is (long with aia_header) saved
;   to the file: full_save_path+YYYYMMDD_hhmm-hhmm_aia_WAVELENGTH_FILETAG_lv_1_5.sav'
;
;HISTORY: Name---------Date---------Description
;         M Weberg  18 MAY, 2017   Initial coding
;         M Weberg  18 SEPT, 2017  Updated documentation and split off a helper function
;-

FUNCTION SELECT_FILES_IN_TIME_RANGE, test_fits_files, start_time, end_time
    ;Reads the file header timesteps and selects only the requested subrange
    num_fits_files = N_ELEMENTS(test_fits_files)
    fits_time_in_sec = FLTARR(num_fits_files)

    read_sdo, test_fits_files, test_fits_index, empty_fits_data, /nodata, /uncomp_delete

    split_fits_date = STRSPLIT(test_fits_index.date_obs, 'T:', /EXTRACT)
    FOR i=0,(num_fits_files-1) DO BEGIN
        fits_time_in_sec[i] = FLOAT(split_fits_date[i,1])*3600.0+FLOAT(split_fits_date[i,2])*60.0+FLOAT(split_fits_date[i,3])
    ENDFOR

    ;note: currently ignores input seconds (if any). Selects time based only on hour and min
    split_start_time = STRSPLIT(start_time, ':', /EXTRACT)
    split_end_time = STRSPLIT(end_time, ':', /EXTRACT)
    start_sec = FLOAT(split_start_time[0])*3600.0+FLOAT(split_start_time[1])*60
    end_sec = FLOAT(split_end_time[0])*3600.0+FLOAT(split_end_time[1])*60

    ;finds the data within the given time range
    selected_indices = WHERE((fits_time_in_sec GE start_sec) AND (fits_time_in_sec LE end_sec))

    IF selected_indices[-1] EQ -1 THEN MESSAGE, 'ERROR: no data found within the given time range!'

    RETURN, selected_indices
END

;###############################################################################
;@main##########################################################################
;###############################################################################
PRO PREPROCESS_AIA_DATA_FOR_NUWT, start_date, end_date, wavelength=wavelength, $
                                  cutout=cutout, x0=x0, y0=y0, xsize=xsize, ysize=ysize, $
                                  unsharp=unsharp, us_width=us_width, $
                                  smooth=smooth, sm_width=sm_width, $
                                  filetag=filetag, $
                                  base_data_folder=base_data_folder, $
                                  full_load_path=full_load_path, $
                                  full_save_path=full_save_path, $
                                  input_prefix=input_prefix, $
                                  fits_ext=fits_ext, $
                                  save_lv_1_5=save_lv_1_5

COMPILE_OPT IDL2
;@set_defaults##################################################################
; SET DEFAULT PARAMETER VALUES
;###############################################################################
IF NOT KEYWORD_SET(wavelength) THEN wavelength = 171
IF KEYWORD_SET(cutout) AND NOT KEYWORD_SET(filetag) THEN filetag = 'cutout'
IF NOT KEYWORD_SET(filetag) THEN filetag = 'full_disk'
IF NOT KEYWORD_SET(base_data_folder) THEN base_data_folder = '/DATADRIVE1/data/'
IF NOT KEYWORD_SET(full_save_path) THEN BEGIN
    full_save_path = base_data_folder+'SDO/processed/'
ENDIF
IF STRLEN(base_data_folder) GT 0 AND NOT base_data_folder.endswith('/') THEN base_data_folder = base_data_folder+'/'
IF STRLEN(full_save_path) GT 0 AND NOT full_save_path.endswith('/') THEN full_save_path = full_save_path+'/'

IF NOT KEYWORD_SET(input_prefix) THEN input_prefix = 'aia'
IF NOT KEYWORD_SET(fits_ext) THEN fits_ext = '.fits'
IF (NOT KEYWORD_SET(unsharp) AND NOT KEYWORD_SET(smooth)) THEN save_lv_1_5 = 1

;Default coords select the south polar region (reminder: AIA coords are 1 indexed)
;These values are only used if /cutout is set!
IF N_ELEMENTS(x0) EQ 0 THEN x0 = 1025
IF N_ELEMENTS(y0) EQ 0 THEN y0 = 1
IF N_ELEMENTS(xsize) EQ 0 THEN xsize = 2048
IF N_ELEMENTS(ysize) EQ 0 THEN ysize = 1024

IF NOT KEYWORD_SET(us_width) THEN us_width = 10
IF NOT KEYWORD_SET(sm_width) THEN sm_width = 3

PRINT, '  ' ;just an empty line to make console output easier to read
PRINT, '##### RUNNING PREPROCESS_AIA_DATA_FOR_NUWT.pro #####'
start_time = SYSTIME(/seconds)

;@read_dates####################################################################
; READING AND CONVERTING INPUT TIMESTAMPS
;###############################################################################
PRINT, '   Parsing input dates ...'
split_start_date = STRSPLIT(start_date, '/- ', /extract) ;note: splits on a space too
start_year = split_start_date[0]
start_month = NUWT_CONVERT_MONTH_NAME(split_start_date[1])
start_day = split_start_date[2]
IF STRLEN(start_month) EQ 1 THEN start_month = '0'+start_month
IF STRLEN(start_day) EQ 1 THEN start_day = '0'+start_day

split_end_date = STRSPLIT(end_date, '/- ', /extract);note: splits on a space too
end_year = split_end_date[0]
end_month = NUWT_CONVERT_MONTH_NAME(split_end_date[1])
end_day = split_end_date[2]
IF STRLEN(end_month) EQ 1 THEN end_month = '0'+end_month
IF STRLEN(end_day) EQ 1 THEN end_day = '0'+end_day

IF start_day NE end_day THEN MESSAGE, 'ERROR: please select a time range within a single day!'

split_start_time = STRSPLIT(split_start_date[3], ':', /EXTRACT)
start_hour = split_start_time[0]
start_min = split_start_time[1]
IF STRLEN(start_hour) EQ 1 THEN start_hour = '0'+start_hour
IF STRLEN(start_min) EQ 1 THEN start_min = '0'+start_min

split_end_time = STRSPLIT(split_end_date[3], ':', /EXTRACT)
end_hour = split_end_time[0]
end_min = split_end_time[1]
IF STRLEN(end_hour) EQ 1 THEN end_hour = '0'+end_hour
IF STRLEN(end_min) EQ 1 THEN end_min = '0'+end_min

;timestamp for output files
save_timestamp = start_year+start_month+start_day+'_'+start_hour+start_min+'-'+end_hour+end_min

;@check_dir#####################################################################
; CHECKING AND MAKING SAVE FOLDERS AS NEEDED
;###############################################################################
PRINT, '   Checking input and output folders ...'
IF NOT KEYWORD_SET(full_load_path) THEN BEGIN
    full_load_path = base_data_folder+'SDO/AIA/'+start_year+'/'+start_month+'/'+start_day+'/'+str(wavelength)+'/level1/'
ENDIF
IF STRLEN(full_load_path) GT 0 AND NOT full_load_path.endswith('/') THEN full_load_path = full_load_path+'/'
IF NOT FILE_TEST(full_load_path) THEN MESSAGE, 'ERROR: no local data is available for the selected day!'

IF NOT FILE_TEST(full_save_path) THEN BEGIN
    FILE_MKDIR, full_save_path
    PRINT, '   Created folder: '+full_save_path
ENDIF

;@load_data#####################################################################
; LOADING DATA WITHIN THE REQUESTED CUTOUT AND RUNING AIA_PREP
;###############################################################################
PRINT, '   Loading in data and running aia_prep ...'
fits_files = findfile(full_load_path+input_prefix+'*'+fits_ext)

select_indices = select_files_in_time_range(fits_files, split_start_date[3], split_end_date[3])

IF KEYWORD_SET(cutout) THEN BEGIN
    READ_SDO, fits_files[select_indices], lv_1_index, lv_1_data, x0, y0, xsize, ysize, /uncomp_delete
    AIA_PREP, lv_1_index, lv_1_data, aia_index, aia_data, /cutout
ENDIF ELSE BEGIN
    READ_SDO, fits_files[select_indices], lv_1_index, lv_1_data, /uncomp_delete
    AIA_PREP, lv_1_index, lv_1_data, aia_index, aia_data
ENDELSE

DELVAR, lv_1_index, lv_1_data

;@apply_filters#################################################################
; APPLYING UNSHARP MASK AND SMOOTHING IF REQUESTED
;###############################################################################
us_aia_data = 'none' ; placeholders for saving method
sm_aia_data = 'none'
IF KEYWORD_SET(unsharp) THEN BEGIN
    print, '   Applying unsharp mask of width '+STRTRIM(us_width, 2)+' ...'
    us_aia_data = aia_data - SMOOTH(aia_data, [us_width, us_width, 1])
ENDIF

IF KEYWORD_SET(smooth) THEN BEGIN
    print, '   Smoothing data over '+STRTRIM(sm_width, 2)+' timesteps ...'
    sm_aia_data = SMOOTH(aia_data, [1, 1, sm_width])
    IF KEYWORD_SET(unsharp) THEN us_aia_data = SMOOTH(us_aia_data, [1, 1, sm_width])
ENDIF

;@save_data#####################################################################
; SAVING DATA TO OUTPUT FOLDER
;###############################################################################
print, '   Saving processed data ...'
IF KEYWORD_SET(save_lv_1_5) THEN BEGIN
    lv_1_5_filename = save_timestamp+'_aia_'+str(wavelength)+'_'+filetag+'_lv_1_5.sav'
    SAVE, aia_index, aia_data, filename=full_save_path+lv_1_5_filename
    PRINT, 'lv 1.5 filename: '+lv_1_5_filename
ENDIF
IF KEYWORD_SET(unsharp) OR KEYWORD_SET(smooth) THEN BEGIN
    sm_and_us_filename = save_timestamp+'_aia_'+str(wavelength)+'_'+filetag+'_sm_and_us.sav'
    SAVE, aia_index, sm_aia_data, us_aia_data, filename=full_save_path+sm_and_us_filename
    PRINT, 'smoothed & unsharp filename: '+sm_and_us_filename
ENDIF

;@print_status##################################################################
; PRINTING FINAL STATUS AND ASSORTED INFORMATION
;###############################################################################
finished_time = SYSTIME(/seconds)
total_run_time = (finished_time - start_time)/60.0

PRINT, '  ' ;just an empty line to make console output easier to read
PRINT, 'Finished preprocessing data!'
PRINT, 'Total time elapsed: '+STRTRIM(STRING(total_run_time, format='(F8.2)'), 2)+' min'
PRINT, 'All files saved to the following directory:'
PRINT, '   '+full_save_path

END
