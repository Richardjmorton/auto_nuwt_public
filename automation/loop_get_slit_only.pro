;+
;NAME: LOOP_GET_SLIT_ONLY
;
;PURPOSE:
;   [to write later]
;
;INPUTS:
;   [to write later]
;
;OPTIONAL INPUTS:
;   [to write later]
;
;OUTPUTS:
;   various
;
;HISTORY: Name---------Date---------Description
;         M Weberg 13 JUNE, 2017  Initial coding.
;
;TO-DO / LIMITATIONS:
;   - GENERALIZE THE PROCEDURE SO THERE ARE NO HARD CODED VALUES!
;   - Give more explict control over how multiple slits from the same datacube are
;     extracted
;-

PRO LOOP_GET_SLIT_ONLY, timestamps=timestamps, img_locs=img_locs, $
                                regions=regions, $
                                slit_type=slit_type, $
                                radius=radius, start_ang=start_ang, $
                                subtend=subtend, width=width, $
                                clockwise=clockwise, $
                                x1=x1, x2=x2, y1=y1, y2=y2, $
                                load_folder=load_folder, $
                                save_folder=save_folder, $
                                rerun_nuwt=rerun_nuwt, $
                                old_slit_folder=old_slit_folder

COMPILE_OPT IDL2
;@defaults######################################################################
;SETTING DEFAULT VALUES AND VALIDATING INPUTS
;###############################################################################
IF NOT KEYWORD_SET(slit_type) THEN slit_type = 'arc'
us_width = 10

; IF NOT KEYWORD_SET(load_folder) THEN CD, CURRENT=load_folder ;i.e. defaults to current folder
IF NOT KEYWORD_SET(load_folder) THEN load_folder = '/DATADRIVE1/data/SDO/processed/'
IF NOT load_folder.endswith('/') THEN load_folder = load_folder+'/'
IF NOT KEYWORD_SET(save_folder) THEN CD, CURRENT=save_folder ;i.e. defaults to current folder
IF NOT save_folder.endswith('/') THEN save_folder = save_folder+'/'
IF NOT KEYWORD_SET(old_slit_folder) THEN CD, CURRENT=old_slit_folder ;i.e. defaults to current folder
IF NOT old_slit_folder.endswith('/') THEN old_slit_folder = old_slit_folder+'/'

; yearly run list ['20100523_0400-0800', '20110523_0000-0400', '20120523_0000-0400', '20130523_0000-0400', '20140523_0000-0400', '20150523_0000-0400', '20160523_0000-0400', '20170523_0000-0400']
; p1_times = ['20100523_0400-0800', '20100523_0400-0800', '20101006_0000-0400', '20120327_1600-2000', '20120327_1600-2000']
; p1_img_locs = ['south', 'south', 'north', 'north', 'north']
; p1_regions = ['ch', 'qs', 'ch', 'ch', 'ch']
; p1_radii = [15, 15, 15, 7, 15]
; p1_start_ang = [265, 250, 95, 95, 95]

; p3_ts = ['20100523_0400-0800', '20110523_0000-0400', '20120527_0000-0400', '20130424_1600-2000', '20140501_0000-0400', '20150523_0000-0400', '20160523_0000-0400', '20170523_0000-0400']
;
; p3_locs = ['south']
; p3_regions = ['ch']
; p3_rads = [5.0, 7.5, 10.0 ,12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0]
; p3_start_ang = [265.0, 265.0, 265.0, 268.0, 225.0, 265.0, 265.0, 265.0]

start_time = SYSTIME(/seconds)
;@find_files####################################################################
;ASSEMBLING FILE PATHS (OR SEARCHING FOR FILES IN THE GIVEN FOLDER)
;###############################################################################
IF KEYWORD_SET(timestamps) OR KEYWORD_SET(img_locs) THEN BEGIN

    ; restore, filename='/DATADRIVE1/data/SDO/processed/'+FILE_TS+'_aia_171_south_sm_and_us.sav'

    ;Contruct the filenames using input values.
    img_loc_sep = ''
    IF KEYWORD_SET(timestamps) AND KEYWORD_SET(img_locs) THEN img_loc_sep = '_'
    IF NOT KEYWORD_SET(timestamps) THEN timestamps = ''
    IF NOT KEYWORD_SET(img_locs) THEN img_locs = ''

    num_timestamps = N_ELEMENTS(timestamps)
    num_img_locs = N_ELEMENTS(img_locs)

    ;cast single values into arrays with the correct sizes
    IF num_timestamps EQ 1 THEN timestamps = strarr(num_img_locs) + timestamps[0]
    IF num_img_locs EQ 1 THEN img_locs = strarr(num_timestamps) + img_locs[0]

    IF (num_timestamps NE num_img_locs) AND $
       (num_timestamps GT 1 AND num_img_locs GT 1)THEN BEGIN
        PRINT, 'ERROR: number of input timestamps and img_locs are incompatible!'
        MESSAGE, 'Please input either single values or arrays with the same number of elements'
    ENDIF

    ; input_files = load_folder+timestamps+'_aia_171'+img_loc_sep+img_locs+'_sm_and_us.sav'
    input_files = load_folder+timestamps+'_aia_171'+img_loc_sep+img_locs+'_lv_1_5.sav'

ENDIF ELSE BEGIN
    ;Search "load_folder" for ALL .sav files (WARNING: USE WITH CARE!)
    input_files = findfile(load_folder+'*.sav')
ENDELSE

num_files = N_ELEMENTS(input_files)
IF num_files GE 1 THEN BEGIN
    PRINT, ' '
    PRINT, 'Found '+STRTRIM(num_files, 2)+' input file(s) in '+load_folder
ENDIF ELSE BEGIN
    MESSAGE, 'ERROR: No input files found!'
ENDELSE

;@slit_params###################################################################
;VALIDATING ARRAYS OF SLIT PARAMETERS
;###############################################################################
IF KEYWORD_SET(regions) THEN BEGIN
    IF N_ELEMENTS(regions) NE num_files THEN BEGIN
        slit_regions = strarr(num_files) + regions[0]
    ENDIF ELSE BEGIN
        slit_regions = regions
    ENDELSE
ENDIF ELSE BEGIN
    slit_regions = strarr(num_files) + slit_type
ENDELSE

IF STRLOWCASE(slit_type) EQ 'arc' THEN BEGIN

    ;[2017-09-19 update] slit radii has been given special funcationality!
    IF NOT KEYWORD_SET(radius) THEN BEGIN
        ;[default]
        slit_radii = fltarr(num_files) + 15.0
        num_slits_per_file = 1
    ENDIF ELSE BEGIN
        IF N_ELEMENTS(radius) EQ num_files THEN BEGIN
            ;each file has its own defined slit radius
            slit_radii = radius
            num_slits_per_file = 1
        ENDIF ELSE IF N_ELEMENTS(radius) EQ 1 THEN BEGIN
            ;single radius used for ALL data files
            slit_radii = fltarr(num_files) + radius[0]
            num_slits_per_file = 1
        ENDIF ELSE BEGIN
            ;standard set of multiple slit radii are applied to ALL data files
            slit_radii = radius
            num_slits_per_file = N_ELEMENTS(radius)
        ENDELSE
    ENDELSE

    ; IF KEYWORD_SET(radius) THEN BEGIN ;OLD CODE KEPT FOR NOW...
    ;     IF N_ELEMENTS(radius) NE num_files THEN BEGIN
    ;         slit_radii = fltarr(num_files) + radius[0]
    ;     ENDIF ELSE BEGIN
    ;         slit_radii = radius
    ;     ENDELSE
    ; ENDIF ELSE BEGIN
    ;     slit_radii = fltarr(num_files) + 15.0
    ; ENDELSE

    ;Each of the folowing varibles behave similarly
    IF KEYWORD_SET(start_ang) THEN BEGIN
        IF N_ELEMENTS(start_ang) NE num_files THEN BEGIN
            slit_start_ang = fltarr(num_files) + start_ang[0]
        ENDIF ELSE BEGIN
            slit_start_ang = start_ang
        ENDELSE
    ENDIF ELSE BEGIN
        slit_start_ang = fltarr(num_files) + 95.0
    ENDELSE

    IF KEYWORD_SET(subtend) THEN BEGIN
        IF N_ELEMENTS(subtend) NE num_files THEN BEGIN
            slit_subtend = fltarr(num_files) + subtend[0]
        ENDIF ELSE BEGIN
            slit_subtend = subtend
        ENDELSE
    ENDIF ELSE BEGIN
        slit_subtend = fltarr(num_files) + 10.0
    ENDELSE

    IF KEYWORD_SET(width) THEN BEGIN
        IF N_ELEMENTS(width) NE num_files THEN BEGIN
            slit_width = fltarr(num_files) + width[0]
        ENDIF ELSE BEGIN
            slit_width = width
        ENDELSE
    ENDIF ELSE BEGIN
        slit_width = intarr(num_files) + 3
    ENDELSE

    ; IF KEYWORD_SET(clockwise) THEN BEGIN
    ;     IF N_ELEMENTS(clockwise) NE num_files THEN BEGIN
    ;         slit_direction = fltarr(num_files) + clockwise[0]
    ;     ENDIF ELSE IF num_files EQ 1 THEN BEGIN
    ;         slit_direction = [1]
    ;     ENDIF ELSE BEGIN
    ;         slit_direction = clockwise
    ;     ENDELSE
    ; ENDIF ELSE BEGIN
    ;     loc_north_hemi = WHERE(slit_start_ang LE 180.0)
    ;     slit_direction = intarr(num_files)
    ;     slit_direction[loc_north_hemi] = 1
    ; ENDELSE

    slit_direction = intarr(num_files)

ENDIF ELSE IF STRLOWCASE(slit_type) EQ 'diag' THEN BEGIN
    print, 'diag slit currently not enabled (will be added in the future)'
ENDIF

FOR f=0, (num_files-1) DO BEGIN
    print, ' '
    print, ' '
    print, ' '
    print, 'Now selecting slit(s) and running NUWT for: '+input_files[f]

    ;LOAD the the AIA data cube
    IF NOT KEYWORD_SET(rerun_nuwt) THEN BEGIN
        RESTORE, filename=input_files[f]

        ;Apply unsharp mask to the data
        us_aia_data = aia_data - SMOOTH(aia_data, [us_width, us_width, 1])
        us_notes = 'Unsharp masked using smooth widths of [10, 10, 1]'
    ENDIF

    FOR slitnum=0, (num_slits_per_file-1) DO BEGIN
        IF num_slits_per_file EQ 1 THEN BEGIN
            current_slit_radius = slit_radii[f]
        ENDIF ELSE BEGIN
            current_slit_radius = slit_radii[slitnum]
        ENDELSE

        lowcase_img_loc = strlowcase(img_locs[f])
        slit_filename = timestamps[f]+img_loc_sep+lowcase_img_loc+'_'+strlowcase(slit_regions[f])+$
                        '_'+strtrim(string(current_slit_radius, format='(f4.1)'), 2)+'Mm_td.sav'

        nuwt_tag = timestamps[f]+img_loc_sep+lowcase_img_loc+'_'+strlowcase(slit_regions[f])+$
                   '_'+strtrim(string(current_slit_radius, format='(f4.1)'), 2)+'Mm_us'

        nuwt_plot_header = timestamps[f]+', '+lowcase_img_loc.CapWords()+' '+STRUPCASE(slit_regions[f])+', '$
                           +strtrim(string(current_slit_radius, format='(f4.1)'), 2)+'Mm, Unsharp'

        IF NOT KEYWORD_SET(rerun_nuwt) THEN BEGIN
            IF STRLOWCASE(slit_type) EQ 'arc' THEN BEGIN
                print, 'Slit parameters:'
                print, '   radius = '+STRTRIM(current_slit_radius, 2)
                print, '   start_ang = '+STRTRIM(slit_start_ang[f], 2)
                print, '   subtend = '+STRTRIM(slit_subtend[f], 2)
                print, '   width = '+STRTRIM(slit_width[f], 2)
                IF slit_direction[f] EQ 1 THEN print, '   clockwise sampling direction'
                ; NUWT_ARC_SLIT, sm_aia_data, sm_slit_td, index=aia_index, radius=current_slit_radius, units_radius='Mm', /rel_limb, $
                NUWT_ARC_SLIT, aia_data, slit_td, index=aia_index, radius=current_slit_radius, units_radius='Mm', /rel_limb, $
                          start_ang=slit_start_ang[f], subtend=slit_subtend[f], width=slit_width[f], $
                          clockwise=slit_direction[f], /noopen, meta_out=slit_meta
                NUWT_ARC_SLIT, us_aia_data, us_slit_td, index=aia_index, radius=current_slit_radius, units_radius='Mm', /rel_limb, $
                          start_ang=slit_start_ang[f], subtend=slit_subtend[f], width=slit_width[f], $
                          clockwise=slit_direction[f], /noopen, notes=us_notes
            ENDIF ELSE IF STRLOWCASE(slit_type) EQ 'diag' THEN BEGIN
                print, 'diag slit currently not enabled (will be added in the future)'
            ENDIF

            slit_td_errs = CALC_AIA_ERRORS(slit_td)
            ; save, slit_meta, sm_slit_td, us_slit_td, filename=save_folder+slit_filename ;OUTDATED!
            save, slit_td, us_slit_td, slit_td_errs, filename=save_folder+slit_filename
        ENDIF ELSE BEGIN
            ;reloading slits previously extracted
            PRINT, 'Using previously extracted data slit!'
            RESTORE, filename=old_slit_folder+slit_filename
        ENDELSE


        ; PLOT_NUWT_SLIT_LOC, sm_aia_data, slit_meta=slit_meta, title=timestamps[f]+', SDO / AIA 171, 15Mm', filetag=FILE_TS+'_15Mm'

        RUN_NUWT, us_slit_td, errors=slit_td_errs, /aia, grad=0.5, min_tlen=20, /gauss, /pad_fft, $
                  /smooth, sm_width = [1, 3], $
                  default_plot_header=nuwt_plot_header, save_folder=save_folder, filetag=nuwt_tag

    ENDFOR
ENDFOR

;@print_status##################################################################
; PRINTING FINAL STATUS AND ASSORTED INFORMATION
;###############################################################################
finished_time = SYSTIME(/seconds)
total_run_time = (finished_time - start_time)/60.0

PRINT, '  ' ;just an empty line to make console output easier to read
PRINT, 'Finished extracting slits and running NUWT!'
PRINT, 'Total time elapsed: '+STRTRIM(STRING(total_run_time, format='(F8.2)'), 2)+' min'
PRINT, 'All files saved to the following directory:'
PRINT, '   '+save_folder

END
