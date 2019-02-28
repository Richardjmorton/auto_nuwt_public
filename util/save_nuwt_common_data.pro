;+
;NAME: SAVE_NUWT_COMMON_DATA
;
;PURPOSE:
;   Saving the master "all_nuwt_dat" COMMON block data to a '.sav' file
;
;INPUTS:
;   save_folder - folder in which to save the results. Defaults to the user's home folder.
;                 Will also append a '/' to the end if not included.
;   filename - name of the file to save the results in
;   filetag - unique tag to append to the START of the output file
;
;OUTPUTS:
;   A '.sav' file containing the with the lists of structures currently stored
;   "all_nuwt_dat" COMMON block.
;
;HISTORY: Name---------Date---------Description
;         M Weberg  19 SEPT, 2016  Initial coding
;         M Weberg  26 JULY, 2017  Updated code to save "nuwt_meta_dat" too
;         M Weberg  18 SEPT, 2017  Renamed "nuwt_fft_results" to "nuwt_fft_spec"
;-

PRO SAVE_NUWT_COMMON_DATA, save_folder=save_folder, filename=filename, filetag=filetag

COMPILE_OPT IDL2

COMMON all_nuwt_dat, nuwt_located, nuwt_threads, nuwt_fft_spec, nuwt_fft_peaks
COMMON nuwt_meta_dat, nuwt_meta, nuwt_bulk_stats

;###############################################
;SETTING DEFAULT VALUES
;###############################################
IF NOT KEYWORD_SET(filename) THEN filename = 'nuwt_results'
IF NOT KEYWORD_SET(filetag) THEN filetag='' ELSE filetag = filetag+'_'
IF NOT KEYWORD_SET(save_folder) THEN save_folder = '' ;i.e. defaults to home folder
IF strlen(save_folder) GT 0 AND NOT save_folder.endswith('/') THEN save_folder = save_folder+'/'

save_filename = filetag+filename+'.sav'
print, 'Saving NUWT "all_nuwt_dat" & "nuwt_meta_dat" COMMON blocks to file:'
print, '   Save folder: '+save_folder
print, '   filename: '+save_filename
SAVE, nuwt_meta, nuwt_bulk_stats, nuwt_located, nuwt_threads, nuwt_fft_spec, nuwt_fft_peaks, FILENAME=save_folder+save_filename

END
