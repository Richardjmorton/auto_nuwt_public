;+
;NAME: NUWT_SUBPLOT_GRID
;
;PURPOSE:
;   Generate the position values for a grid of subplots. Options are included
;   for controlling the amount of whitespace between the subplots and
;   specifying custom width and height ratios. This function is inspired by the
;   GridSpec() function in the matplotlib module of Python.
;
;   IMPORTANT! By default, both the input and output coordinates are assumed to
;   be in normalized units ranging from 0.0 to 1.0 and are measured relative to
;   LOWER LEFT corner of the plotting window. HOWEVER, if the "window_dimensions"
;   parameter is set to a two-element array with the width and height of the plot
;   window (in pixels), then any input values in normalized units (values <= 1.0)
;   will be automatically converted to pixels and all output coordinates will be
;   given in units of pixels. As a reminder, pixel values can only be given to
;   the POSITION keyword of an IDL plot function if the /DEVICE keyword is also set.
;
;REQUIRED INPUTS:
;   ncol - number of columns
;   nrow - number of rows
;
;OPTIONAL INPUTS
;   window_dimensions - Two-element array of [width, height] giving the
;                       dimensions of the full plot window in units of pixels.
;                       If set, input values in normalized units (values <= 1.0)
;                       will be automatically converted and all output coordinates
;                       values be given in units of pixels.
;
;   left - Left edge of the entire subplot grid. Default is 0.10
;   bottom - Bottom edge of the entire subplot grid. Default is  0.10
;   right - Right edge of the entire subplot grid. Default is 0.90
;   top - Top edge of the entire subplot grid. Default is 0.90
;
;   hspace - Horizontal whitespace between subplots. The user can set different
;            spacing for each gap between plots by passing an array of ncol - 1
;            values. Horizontal gaps are ordered from left to right. Default is 0.05
;
;   vspace - Vertical whitespace between subplots. The user can set different
;            spacing for each gap between plots by passing an array of nrow - 1
;            values. Vertical gaps are ordered from top to bottom. Default is 0.05
;
;   width_ratios - Ratios of the widths of each column relative to each other.
;                  If used, MUST be given an array with ncol total values.
;                  Ratios may given as integers or flaoting point values.
;                  By default, all columns will have the same width.
;
;   height_ratios - Ratios of the heights of each row relative to each other.
;                   If used, MUST be given an array with nrow total values.
;                   Ratios may given as integers or flaoting point values.
;                   By default, all rows will have the same height.
;
;OUTPUT:
;   Returns an array of shape [ncol, nrow, 4]. The z-axis of the array contains
;   the the full set of [left, bottom, right, top] postion values in the order
;   expected by the POSITION keyword of all IDL plotting functions.
;
;   Please note, the rows and columns are zero indexed and start with the
;   UPPER LEFT subplot. Example: subplot_grid[1,3,*] would be the position
;   values the subplot located in the second column from the left and the fourth
;   row from the top.
;
;HISTORY: Name---------Date---------Description
;         M Weberg  22 MAY, 2018  Initial coding
;         M Weberg  25 MAY, 2018  Added WINDOW_DIMENSIONS keyword
;-

FUNCTION NUWT_SUBPLOT_GRID, ncol, nrow, window_dimensions=window_dimensions, $
                            left=left, bottom=bottom, right=right, top=top, $
                            hspace=hspace, vspace=vspace, $
                            width_ratios=width_ratios, height_ratios=height_ratios

COMPILE_OPT IDL2
;@defaults######################################################################
; VALITDATE INPUT VALUES AND SET DEFAULTS
;###############################################################################
num_cols = MAX([FIX(ncol), 1])
num_rows = MAX([FIX(nrow), 1])
num_subplots = num_cols*num_rows

IF N_ELEMENTS(left) NE 1 THEN left = 0.10 ELSE left = FLOAT(left[0])
IF N_ELEMENTS(bottom) NE 1 THEN bottom = 0.10 ELSE bottom = FLOAT(bottom[0])
IF N_ELEMENTS(right) NE 1 THEN right = 0.90 ELSE right = FLOAT(right[0])
IF N_ELEMENTS(top) NE 1 THEN top = 0.90 ELSE top = FLOAT(top[0])

IF N_ELEMENTS(hspace) LT 1 THEN hspace = 0.05 ELSE hspace = FLOAT(hspace)
IF num_cols GT 1 THEN BEGIN
    IF N_ELEMENTS(hspace) EQ 1 THEN BEGIN
        hspace_arr = FLTARR(num_cols-1) + hspace
    ENDIF ELSE IF N_ELEMENTS(hspace) EQ (num_cols-1) THEN BEGIN
        hspace_arr = hspace
    ENDIF ELSE BEGIN
        PRINT, 'ERROR: Invalid number of elements in hspace!'
        MESSAGE, 'Please input either a single value or an array of ncol - 1 values.'
    ENDELSE
ENDIF ELSE BEGIN
    ;If there is only one column, ignore hspace.
    hspace_arr = 0.0
ENDELSE

IF N_ELEMENTS(vspace) LT 1 THEN vspace = 0.05 ELSE vspace = FLOAT(vspace)
IF num_rows GT 1 THEN BEGIN
    IF N_ELEMENTS(vspace) EQ 1 THEN BEGIN
        vspace_arr = FLTARR(num_rows-1) + vspace
    ENDIF ELSE IF N_ELEMENTS(vspace) EQ (num_rows-1) THEN BEGIN
        vspace_arr = vspace
    ENDIF ELSE BEGIN
        PRINT, 'ERROR: Invalid number of elements in vspace!'
        MESSAGE, ' Please input either a single value or an array of nrow - 1 values.'
    ENDELSE
ENDIF ELSE BEGIN
    ;If there is only one row, ignore vspace.
    vspace_arr = 0.0
ENDELSE

IF N_ELEMENTS(width_ratios) LE 1 THEN BEGIN
    width_ratio_arr = FLTARR(num_cols) + 1.0
ENDIF ELSE IF N_ELEMENTS(width_ratios) EQ num_cols THEN BEGIN
    width_ratio_arr = FLOAT(width_ratios)
ENDIF ELSE BEGIN
    PRINT, 'ERROR: Invalid number of elements in width_ratios!'
    MESSAGE, ' Please input either a single value or an array of ncol values.'
ENDELSE

IF N_ELEMENTS(height_ratios) LE 1 THEN BEGIN
    height_ratio_arr = FLTARR(num_rows) + 1.0
ENDIF ELSE IF N_ELEMENTS(height_ratios) EQ num_rows THEN BEGIN
    height_ratio_arr = FLOAT(height_ratios)
ENDIF ELSE BEGIN
    PRINT, 'ERROR: Invalid number of elements in width_ratios!'
    MESSAGE, ' Please input either a single value or an array of nrow values.'
ENDELSE

;@convert_to_pxls###############################################################
; IF WINDOW_DIMENSIONS IS SET, CONVERT TO UNITS OF PIXELS
;###############################################################################
IF N_ELEMENTS(window_dimensions) EQ 1 OR N_ELEMENTS(window_dimensions) GT 2 THEN BEGIN
    PRINT, 'ERROR: Invalid number of elements in window_dimensions!'
    MESSAGE, 'Please either input a two-element array or omit the keyword'
ENDIF ELSE IF N_ELEMENTS(window_dimensions) EQ 2 THEN BEGIN
    win_width = FIX(window_dimensions[0])
    win_height = FIX(window_dimensions[1])
    ;Values < 1.0 are assumed to be in normalized units
    IF left LT 1.0 THEN left = left*win_width
    IF bottom LT 1.0 THEN bottom = bottom*win_height
    IF right LT 1.0 THEN right = right*win_width
    IF top LT 1.0 THEN top = top*win_height
    IF hspace_arr[0] LT 1.0 THEN hspace_arr = hspace_arr*win_width
    IF vspace_arr[0] LT 1.0 THEN vspace_arr = vspace_arr*win_height
ENDIF

;@calc_coords###################################################################
; CALCULATE THE SUPLOT COORDINATES
;###############################################################################

;Left and right subplot edges
available_width = right - left - TOTAL(hspace_arr)
unit_width = available_width/TOTAL(width_ratio_arr)
col_widths = width_ratio_arr*unit_width
subplot_lefts = FLTARR(num_cols, num_rows)
subplot_rights = FLTARR(num_cols, num_rows)
FOR c=0, (num_cols-1) DO BEGIN
    IF c EQ 0 THEN BEGIN
        subplot_lefts[0,*] = left
        subplot_rights[0,*] = left + col_widths[0]
    ENDIF ELSE BEGIN
        subplot_lefts[c,*] = subplot_rights[c-1,0] + hspace_arr[c-1]
        subplot_rights[c,*] = subplot_lefts[c,0] + col_widths[c]
    ENDELSE
ENDFOR

;Bottom and top subplot edges
available_height = top - bottom - TOTAL(vspace_arr)
unit_height = available_height/TOTAL(height_ratio_arr)
row_heights = height_ratio_arr*unit_height
subplot_bots = FLTARR(num_cols, num_rows)
subplot_tops = FLTARR(num_cols, num_rows)
FOR r=0, (num_rows-1) DO BEGIN
    IF r EQ 0 THEN BEGIN
        subplot_tops[*,0] = top
        subplot_bots[*,0] = top - row_heights[0]
    ENDIF ELSE BEGIN
        subplot_tops[*,r] = subplot_bots[0,r-1] - vspace_arr[r-1]
        subplot_bots[*,r] = subplot_tops[0,r] - row_heights[r]
    ENDELSE
ENDFOR

;Combine all the postion values into a single output array
subplot_grid = FLTARR(num_cols, num_rows, 4)
subplot_grid[*,*,0] = subplot_lefts
subplot_grid[*,*,1] = subplot_bots
subplot_grid[*,*,2] = subplot_rights
subplot_grid[*,*,3] = subplot_tops

RETURN, subplot_grid
END
