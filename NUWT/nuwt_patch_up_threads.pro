;+
;NAME: NUWT_PATCH_UP_THREADS
;
;FUNCTION:
;   'Patchs up' the threads found from 'locate_things.pro', followed by
;   'follow_thread.pro'. Follow thread can skip forward time frames leaving
;   zero results. This uses interpolation (or simple value filling) to fix
;   these gaps and weights these results with a larger error.
;
;PROCEDURE OUTLINE:
;   The first step is to remove any examples where the number of pixels in the
;   thread is less than 2 and any threads where more than half the pixels could
;   not be found with locate things (i.e. given the value 0). The routine then
;   locates any points that have the value 0. These values are then replaced via
;   linear interpolations given a large error (1 pixel). This is to ensure the
;   weighted fitting routines give limited significance to this value.
;
;OPTIONAL INPUTS:
;   fit_flag - {0-3} defines which time series to work on (see follow_thread_fg.pro/moscill.pro)
;   /simp_fill - if set, will fill missing values with the last non-zero value instead of interpolating
;
;HISTORY:
;   R Morton   May 2014  - Created initial program
;   M Weberg  June 2016  - Cleaned-up code, fixed the call to LINFILL,
;                          and added a simple fill option
;   M Weberg   Mar 2018  - Removed and replaced calls to LINFILL.PRO and CONSEC.PRO
;                          This should make the code more portable and easier to
;                          extend  or debug in the future
;
;TO DO OR THINK ABOUT:
;   - Is linear interpolation the best option?
;   - Should weighting be calculated using traditional error analysis?
;   - Does funny things if first elements in array is a 0 - need to look into this
;-

PRO NUWT_PATCH_UP_THREADS, fit_flag=fit_flag, simp_fill=simp_fill, debug=debug

COMPILE_OPT IDL2
;Loads common data generated from follow_thread_fg
COMMON threads_dat, threads

IF n_elements(fit_flag) EQ 0 THEN fit_flag=0

sz = size(threads)
n_threads = sz[1]

;Loading in the correct variables to be patched
CASE fit_flag OF
    ;Option for simple fitting structures
    0: BEGIN
        tval = threads.pos
        terr = threads.err_pos
        END
    ;Options for advanced fitting structures (_FG procedures)
    1: BEGIN
        tval = threads.pos
        terr = threads.err_pos
        END
    2: BEGIN
        tval = threads.inten
        terr = threads.err_inten
        END
    3: BEGIN
        tval = threads.wid
        terr = threads.err_wid
        END
ENDCASE

IF NOT keyword_set(simp_fill) THEN BEGIN
    ;Filling gaps with linear interpolation
    FOR j=0,(n_threads-1) DO BEGIN
        loc_data = WHERE(tval[*,j] GT 0.0)
        loc_zeroes = WHERE(tval[*,j] EQ 0.0)
        len_thread = N_ELEMENTS(WHERE(tval[*,j] GE 0.0))

        ;skips threads with less than 2 positive values
        IF n_elements(loc_data) GT 2 THEN BEGIN

            ;skips threads where half the data points are set to zero
            ;(i.e. follow_thread.pro skipped too many timesteps)
            IF n_elements(loc_zeroes) LT 0.5*len_thread THEN BEGIN

                index_diff = loc_data - SHIFT(loc_data, 1)
                ;note, a differneces in index number GT 1 indicae the presence of a data gap
                gap_ends = WHERE(index_diff GT 1) ;reminder, index location in LOC_DATA array!
                IF gap_ends[0] GT -1 THEN BEGIN
                    num_gaps = N_ELEMENTS(gap_ends)
                    FOR i=0,(num_gaps-1) DO BEGIN
                        start_ind = loc_data[gap_ends[i]-1]
                        end_ind = loc_data[gap_ends[i]]

                        start_val = tval[start_ind, j]
                        end_val = tval[end_ind, j]

                        num_fill = index_diff[gap_ends[i]] - 1
                        fill_ind_arr = INDGEN(num_fill) + 1 + start_ind ;indices of the fill values
                        m = FLOAT(end_val - start_val)/FLOAT(end_ind - start_ind)
                        tval[start_ind+1:end_ind-1, j] = start_val + m*(fill_ind_arr - start_ind)
                        terr[start_ind+1:end_ind-1, j] = 1.0
                    ENDFOR
                ENDIF

                IF KEYWORD_SET(debug) THEN BEGIN
                    ;Double check to make sure all zeros have been filled
                    ;(This SHOULD never happen. If it does, then something else is wrong...)
                    still_zero = WHERE(tval[*,j] EQ 0.0)
                    IF still_zero[0] GT -1 THEN BEGIN
                        print, 'thread index '+STRTRIM(j, 2)+' still has zeroes!'
                    ENDIF
                ENDIF

                ; OLD VERSION OF THE CODE FOR REFERENCE
                ; ;FIND CONSECUTIVE ZEROS IN THREAD
                ; b = where(tval[*,j] EQ 0)
                ; IF b[0] GT -1 THEN CONSEC, b, lo, hi, num
                ;
                ; ;FILL IN USING LINEAR INTERPOLATION
                ; IF b[0] GT -1 THEN BEGIN
                ;     FOR i=0,num-1 DO BEGIN
                ;         ;note: the messy variable juggling is needed to pass data to LINFILL
                ;         ;using pass-by-reference (so we actually get the results "passed" back out)
                ;         temp_sub_arr = tval[*,j]
                ;         LINFILL, temp_sub_arr, b[lo[i]]-1, b[hi[i]]+1
                ;         tval[*,j] = temp_sub_arr
                ;         terr[lo[i]:hi[i],j] = 1.0
                ;     ENDFOR
                ; ENDIF
                ; IF b[0] GT -1 THEN FOR i=0,num-1 DO terr[lo[i]:hi[i],j] = 1.
                ;
                ; ;FIND SINGLE ZEROS AND FILL VIA LINEAR INTERPOLATION
                ; b = where(tval[*,j] EQ 0)
                ; IF b[0] GT -1 THEN BEGIN
                ;     FOR i=0,(n_elements(b)-2) DO BEGIN
                ;         tval[b[i],j,0] = 0.5*tval[b[i]-1,j] + 0.5*tval[b[i]+1,j]
                ;         ; terr[b[i],j] = sqrt(tval[b[i]-1,j]^2+tval[b[i]+1,j]^2)
                ;         terr[b[i],j] = 1.0
                ;     ENDFOR
                ; ENDIF

            ENDIF ELSE BEGIN
                ;Erases all entries where number of 0 elements gt 1/2 positive elements
                tval[*,j] = -1.0
            ENDELSE

        ENDIF ELSE BEGIN
            ;Erases all entries where number of positive elements lt 2
            tval[*,j] = -1.0
        ENDELSE
    ENDFOR
ENDIF ELSE BEGIN
    ;Simple filling of gaps with the previous positive value
    FOR j=0,(n_threads-1) DO BEGIN
        loc_data = WHERE(tval[*,j] GT 0.0)
        loc_zeroes = WHERE(tval[*,j] EQ 0.0)
        len_thread = N_ELEMENTS(WHERE(tval[*,j] GE 0.0))

        ;skips threads with less than 2 positive values
        IF n_elements(loc_data) GT 2 THEN BEGIN

            ;skips threads where half the data points are set to zero
            ;(i.e. follow_thread.pro skipped too many timesteps)
            IF n_elements(loc_zeroes) LT 0.5*len_thread THEN BEGIN

                ;if value is missing, set to the same value as the last timestep
                IF loc_zeroes[0] NE -1 THEN BEGIN
                    FOR ii=0,(n_elements(loc_zeroes)-1) DO BEGIN
                        tval[loc_zeroes[ii]] = tval[loc_zeroes[ii]-1]
                        terr[loc_zeroes[ii]] = 1.0
                    ENDFOR
                ENDIF

            ENDIF ELSE BEGIN
                ;Erases all entries where number of 0 elements gt 1/2 positive elements
                tval[*,j] = -1.0
            ENDELSE

        ENDIF ELSE BEGIN
            ;Erases all entries where number of positive elements lt 2
            tval[*,j] = -1.0
        ENDELSE
    ENDFOR
ENDELSE

;Saving the patched values back to the input variable arrays
CASE fit_flag OF
    ;Option for simple fitting structures
    0: BEGIN
        threads.pos = tval
        threads.err_pos = terr
        END
    ;Options for advanced fitting structures (_FG procedures)
    1: BEGIN
        threads.pos = tval
        threads.err_pos = terr
        END
    2: BEGIN
        threads.inten = tval
        threads.err_inten = terr
        END
    3: BEGIN
        threads.wid = tval
        threads.err_wid = terr
        END
ENDCASE

END
