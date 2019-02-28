;+
;NAME: RESTORE_NUWT_COMMON_DATA
;
;PURPOSE:
;   Loading in results stored in a '.sav' file created by 'run_nuwt.pro'
;   and placing it back in the 'all_nuwt_dat' master COMMON block
;
;INPUTS:
;   savfile - filename and path for the '.sav' file to be restored
;
;OPTIONAL OUTPUTS:
;   meta_out - structure containing the metadata about the given NUWT run
;   bulk_stats_out - list of structures containing the bulk wave statistics
;                    for each data slit.
;   located_out - list containing the peaks located by NUWT. Each element in the
;                 list contains the stucture for a single dataslit.
;   threads_out - list containing the threads strung together. Each element
;                 in the list contains an array for the given slit. Each array
;                 contains multiple structures, one for each thread in the slit.
;   fft_spec_out - list of the calculated fft spectra. Each element in the list
;                     is, itself, a list containing the results for a single slit.
;                     each sub_structure contains the spectra for a slingle thread.
;   fft_peaks_out - list of the various calculated and selected values. Similar to
;                   "threads_out", each element is an array of structures. One
;                   array each slit & one structure for each thread in the slit
;   /old_run - if set, will not try to export "meta_out" or "bulk_stats_out"
;              prevents the code for crashing on results from very old NUWT runs.
;
;HISTORY: Name---------Date---------Description
;         M Weberg  14 SEPT, 2016  Initial coding
;         M Weberg  18 SEPT, 2017  Renamed "nuwt_fft_results" to "nuwt_fft_spec"
;-

PRO RESTORE_NUWT_COMMON_DATA, savfile, $
                              meta_out=meta_out, $
                              bulk_stats_out=bulk_stats_out, $
                              located_out=located_out, $
                              threads_out=threads_out, $
                              fft_spec_out=fft_spec_out, $
                              fft_peaks_out=fft_peaks_out, $
                              old_run=old_run

COMPILE_OPT IDL2

COMMON all_nuwt_dat, nuwt_located, nuwt_threads, nuwt_fft_spec, nuwt_fft_peaks
COMMON nuwt_meta_dat, nuwt_meta, nuwt_bulk_stats
nuwt_fft_results = 0 ;placeholder for old "fft_spec" list name
nuwt_fft_stats = 0 ;placeholder for old "fft_peaks" list name

RESTORE, savfile
print, 'NUWT results restored to the "nuwt_meta_dat", & "all_nuwt_dat" COMMON blocks'

;update list names (used for loading old NUWT runs)
IF KEYWORD_SET(nuwt_fft_results) THEN nuwt_fft_spec = nuwt_fft_results

;exporting requested lists
IF NOT KEYWORD_SET(old_run) THEN BEGIN
    meta_out = nuwt_meta
    bulk_stats_out = nuwt_bulk_stats
ENDIF

IF KEYWORD_SET(nuwt_fft_stats) THEN BEGIN
    ;needed for backwards compatibility with older NUWT runs (before Feb 2018)
    nuwt_fft_peaks = nuwt_fft_stats
ENDIF

located_out = nuwt_located
threads_out = nuwt_threads
fft_spec_out = nuwt_fft_spec
fft_peaks_out = nuwt_fft_peaks

END
