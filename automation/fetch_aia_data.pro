;+
;NAME: FETCH_AIA_DATA
;
;PURPOSE:
;   Searches for SDO / AIA in the offical VSO database and attempts to
;   automatically download the compressed .fits files to the local drive. Will
;   create local folders as needed.
;
;INPUTS:
;   start_date - string containing the start time and date. Recommended format:
;                'YYYY-MM-DD hh:mm'
;   end_date - string containing the end time and date.
;
;OPTIONAL INPUTS:
;   wavelength - integer value of a valid AIA wavelength (in units of [Angstroms])
;   base_data_folder - top level folder where the data should be downloaded to.
;
;OUTPUTS:
;   If .fits files are found, they will be downloaded to the following folder
;   base_data_folder/SDO/AIA/YYYY/MM/DD/str(wavelength)/level1/
;   Files that already exist on the local computer will NOT be download
;
;NOTES:
;   Often, there will be files that fail to download for some reason or another.
;   The 'fetch_aia_dat' COMMON block contains more information about which files
;   have been downloaded. The common block arrays are also saved to to a file
;   named: YYYYMMDD_AIA_WWW_fetch_status.sav'. You can use the FUNCTION
;   'RETRY_AIA_DATA' to later try to download the missing files.
;
;HISTORY: Name---------Date---------Description
;         M Weberg  23 MAR, 2017  Initial coding
;-

FUNCTION FLAG_MISSING_VSO_FILES, full_vso_index, test_dir, missing_indices=missing_indices
    ;Compares a index of VSO files available vs the headers of files stored locally
    num_vso_files = N_ELEMENTS(full_vso_index)
    file_flags = INTARR(num_vso_files)

    fits_files = findfile(test_dir+'aia*.fits')
    read_sdo, fits_files, fits_index, fits_data, /nodata, /uncomp_delete

    full_vso_index_date = full_vso_index[*].time.start
    trim_fits_date = STRMID(fits_index.date_obs, 0, 19)

    FOR i=0, (num_vso_files-1) DO BEGIN
        file_flags[i] = TOTAL(STRMATCH(trim_fits_date, full_vso_index_date[i]))
    ENDFOR

    missing_indices = WHERE(file_flags EQ 0)

    RETURN, file_flags
END

;---

FUNCTION RETRY_AIA_DATA, filename=filename, use_pseudo_vso_get=use_pseudo_vso_get
    ;Retry downloading VSO data that are missing locally
    ;JSOC ref num for missing single image in May 2014: JSOC_20170613_236
    COMMON fetch_aia_dat, vso_index, download_flag, missing_files, full_save_path

    IF KEYWORD_SET(filename) THEN RESTORE, filename=filename
    PRINT, 'Retrying to download missing AIA files ...'

    IF KEYWORD_SET(use_pseudo_vso_get) THEN BEGIN
        ;(Hopefully) Temporary workaround to downloading issues
        PSEUDO_VSO_GET, vso_index[missing_files], out_dir=full_save_path
    ENDIF ELSE BEGIN
        ; [DEFAULT]
        vso_status = vso_get(vso_index[missing_files], out_dir=full_save_path, /rice)
    ENDELSE

    download_flag = flag_missing_vso_files(vso_index, full_save_path, missing_indices=missing_files)
    SAVE, vso_index, download_flag, missing_files, full_save_path, filename=filename

    RETURN, download_flag
END

;###############################################################################
;@main##########################################################################
;###############################################################################
PRO FETCH_AIA_DATA, start_date, end_date, wavelength=wavelength, sample=sample, $
                    base_data_folder=base_data_folder, use_pseudo_vso_get=use_pseudo_vso_get

COMPILE_OPT IDL2
COMMON fetch_aia_dat, vso_index, download_flag, missing_files, full_save_path

IF NOT KEYWORD_SET(wavelength) THEN wavelength = 171
IF NOT KEYWORD_SET(base_data_folder) THEN base_data_folder = '/DATADRIVE1/data/'

PRINT, '  ' ;just an empty line to make console output easier to read
PRINT, '##### RUNNING FETCH_AIA_DATA.pro #####'
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

split_end_date = STRSPLIT(end_date, '/- ', /extract) ;note: splits on a space too
end_year = split_end_date[0]
end_month = NUWT_CONVERT_MONTH_NAME(split_end_date[1])
end_day = split_end_date[2]
IF STRLEN(end_month) EQ 1 THEN end_month = '0'+end_month
IF STRLEN(end_day) EQ 1 THEN end_day = '0'+end_day

IF start_day NE end_day THEN MESSAGE, 'ERROR: please select a time range with a single day!'

;@check_dir#####################################################################
; CHECKING AND MAKING SAVE FOLDERS AS NEEDED
;###############################################################################
PRINT, '   Checking output folder ...'
full_save_path = base_data_folder+'SDO/AIA/'+start_year+'/'+start_month+'/'+start_day+'/'+str(wavelength)+'/level1/'
IF NOT FILE_TEST(full_save_path) THEN BEGIN
    FILE_MKDIR, full_save_path
    PRINT, '   Created folder: '+full_save_path
ENDIF

;@vso_download##################################################################
; FINDING AND DOWNLOADING DATA FROM THE VSO
;###############################################################################
PRINT, '   Downloading full resolution images from VSO (this may take a while) ...'
PRINT, '  ' ;just an empty line to make console output easier to read
IF KEYWORD_SET(sample) THEN BEGIN
    vso_index = VSO_SEARCH(start_date, end_date, inst='aia', wave=wavelength, pixels=4096, sample=sample)
ENDIF ELSE BEGIN
    vso_index = VSO_SEARCH(start_date, end_date, inst='aia', wave=wavelength, pixels=4096)
ENDELSE
IF STRCMP(TYPENAME(vso_index), 'STRING', /FOLD_CASE) THEN MESSAGE, 'ERROR: no VSO data found!'


IF KEYWORD_SET(use_pseudo_vso_get) THEN BEGIN
    ;(Hopefully) Temporary workaround to downloading issues
    PSEUDO_VSO_GET, vso_index, out_dir=full_save_path
ENDIF ELSE BEGIN
    ; [DEFAULT]
    vso_status = VSO_GET(vso_index, out_dir=full_save_path, /rice)
ENDELSE

download_flag = flag_missing_vso_files(vso_index, full_save_path, missing_indices=missing_files)

num_vso_files = N_ELEMENTS(vso_index)
num_missing = N_ELEMENTS(missing_files)
IF missing_files[-1] EQ -1 THEN num_missing = 0
num_downloaded = num_vso_files - num_missing

;@print_status##################################################################
; PRINTING FINAL STATUS AND ASSORTED INFORMATION
;###############################################################################
finished_time = SYSTIME(/seconds)
total_download_time = (finished_time - start_time)/60.0

PRINT, '  ' ;just an empty line to make console output easier to read
PRINT, 'Finished! Number of files downloaded: '+STRTRIM(num_downloaded, 2)+'/'+STRTRIM(num_vso_files, 2)
PRINT, 'Total time elapsed: '+STRTRIM(STRING(total_download_time, format='(F8.2)'), 2)+' min'
PRINT, 'Files were downloaded to the following directory:'
PRINT, '   '+full_save_path
PRINT, 'See the "fetch_aia_dat" COMMON block for more status information'
IF num_missing GT 0 THEN BEGIN
    fetch_status_filename = start_year+start_month+start_day+'_AIA_'+str(wavelength)+'_fetch_status.sav'
    SAVE, vso_index, download_flag, missing_files, full_save_path, filename=fetch_status_filename
    PRINT, 'Status filename: '+fetch_status_filename
ENDIF ELSE BEGIN
    PRINT, 'All files download. Fetch status file NOT created.'
ENDELSE
END
