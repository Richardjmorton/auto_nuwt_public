;+
;NAME: NUWT_REGION_SELECT
;
;PURPOSE:
;   Simple interactive code that extracts a user selected region from an input
;   data cube (or 2D image)
;
;INPUTS:
;   in_data - 2D or 3D array containing image data
;
;OPTIONAL INPUTS:
;   fn - frame number to plot
;   /log_plot - if set, will log-scale the example image intensity values
;               OUTPUT VALUES WILL BE UNEFFECTED!
;   xsize, ysize - size of graphics window (command for widgets)
;   win_id - set to value of open graphics window you want to use. If not set,
;            opens next free window.
;   /noopen - if set, suppresses the opening of plot windows
;
;OUTPUTS:
;   out_data - output array containing the region selected
;   outvec - [OPTIONAL] array of [x_left, y_bot, x_right, y_top] coordinates used for slit
;
;HISTORY: Name---------Date---------Description
;         M Weberg  8 JUNE, 2016  Initial coding
;         M Weberg  2  MAR, 2018  Refined the region selection process to show a
;                                 preview of currently selected box
;
;TO DO:
;   - Add optional inputs for SDO/AIA header data. Then, output a copy of the
;     header data with updated coordinates and cutout size values.
;   - Related to the above, allow a user to input full-disk AIA fits files
;-

PRO NUWT_REGION_SELECT, in_data, out_data, outvec=outvec, $
                        fn=fn, log_plot=log_plot, $
                        xsize=xsize, ysize=ysize, win_id=win_id, noopen=noopen

COMPILE_OPT IDL2
;###############################################################################
;VALIDATING INPUT DATA
;###############################################################################
;Check dimensions of input data
sz = size(in_data)
IF sz[0] EQ 3 THEN BEGIN
    nx_data = sz[1]
    ny_data = sz[2]
    nt_data = sz[3]
ENDIF ELSE BEGIN
    IF sz[0] EQ 2 THEN BEGIN
        print, '2D input data detected'
        nx_data = sz[1]
        ny_data = sz[2]
        nt_data = 1
        fn = 0
    ENDIF ELSE BEGIN
        MESSAGE, 'ERROR: please input a valid 2D or 3D array'
    ENDELSE
ENDELSE

;###############################################################################
;PLOT EXAMPLE IMAGE FOR REFERENCE AND INTERACTIVE SELECTION
;###############################################################################
IF NOT keyword_set(noopen) THEN BEGIN
    IF n_elements(fn) EQ 0 THEN fn = 0
    IF fn GT (nt_data-1) THEN fn = nt_data-1

    ex_img = in_data[*,*,fn]
    IF KEYWORD_SET(log_plot) THEN BEGIN
        min_ex = MIN(ex_img)
        ex_img = ex_img + ABS(min_ex) + 1
        ex_img = ALOG10(ex_img)
        print, min_ex
    ENDIF

    ;Set plotting elements
    ;writing to buffer pixmap reduces on-screen jitter!
    IF n_elements(win_id) EQ 0 THEN BEGIN
        determine_window_size, xsize=xsize, ysize=ysize, openwin=openwin, /new_win
        window, openwin[0], xs=xsize, ys=ysize
        window_var = openwin[0]
    ENDIF ELSE BEGIN
        wset, win_id
        window_var = win_id
    ENDELSE
    window, xs=xsize, ys=ysize, /pixmap, /free
    pixid = !d.window
    tvim, ex_img
    wset, window_var
    device, copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
ENDIF

;###############################################################################
;QUERY USER FOR SLIT PARAMETERS AS NEEDED
;###############################################################################
;Ask for user to select the lower left and upper right corners

loop_select = 1
ask_select = 'y'
WHILE loop_select DO BEGIN

    IF KEYWORD_SET(noopen) THEN keyboard_select = 1 ELSE keyboard_select = 0

    ;Selecting / setting coroner #1
    PRINT, ''
    IF NOT keyword_set(noopen) THEN BEGIN
        PRINT,'Please left click on the first corner with the cursor'
        PRINT,'[or right click to enter the values using the keyboard]'
        window, xs=xsize, ys=ysize, /pixmap, /free
        pixid = !d.window
        tvim, ex_img
        WSET, window_var
        DEVICE, copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
        CURSOR, x1, y1, /Data, /Down
        IF !MOUSE.button EQ 1 THEN BEGIN
            x1 = FIX(ROUND(x1))
            y1 = FIX(ROUND(y1))
            PRINT, 'x1 = '+STRTRIM(x1, 2)+',   y1 = '+STRTRIM(y1, 2)
        ENDIF ELSE BEGIN
            keyboard_select = 1
        ENDELSE
    ENDIF

    IF keyboard_select THEN BEGIN
        PRINT, 'Please enter the x-coord for the first corner: '
        READ, x1, PROMPT='[waiting for input]> '
        PRINT, 'Please enter the y-coord for the first corner: '
        READ, y1, PROMPT='[waiting for input]> '
        x1 = FIX(ROUND(x1))
        y1 = FIX(ROUND(y1))
    ENDIF

    ;reset the value for the selecting the next corner
    IF KEYWORD_SET(noopen) THEN keyboard_select = 1 ELSE keyboard_select = 0

    ;Selecting / setting corner #2
    PRINT, ''
    IF NOT keyword_set(noopen) THEN BEGIN
        PRINT,'Please left click on the opposite corner with the cursor'
        PRINT,'[or right click to enter the values using the keyboard]'
        IF !MOUSE.button EQ 1 THEN BEGIN
            CURSOR, do_nothing_x, do_nothing_y, /UP, /WAIT
        ENDIF
        !MOUSE.button = 0 ;flushes all clicks
        WHILE !MOUSE.button EQ 0 DO BEGIN
            CURSOR, x2, y2, /change
            window, xs=xsize, ys=ysize, /pixmap, /free
            pixid = !d.window
            tvim, ex_img
            WSET, window_var
            DEVICE, copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
            PLOTS, [x1,x1,x2,x2,x1], [y1,y2,y2,y1,y1]
        ENDWHILE

        IF !MOUSE.button EQ 1 THEN BEGIN
            x2 = FIX(ROUND(x2))
            y2 = FIX(ROUND(y2))
            print, 'x2 = '+strtrim(x2, 2)+',   y2 = '+strtrim(y2, 2)
        ENDIF ELSE BEGIN
            ;revert to keyboard select mode but first clear the preview box from the image
            window, xs=xsize, ys=ysize, /pixmap, /free
            pixid = !d.window
            tvim, ex_img
            WSET, window_var
            DEVICE, copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
            PLOTS, [x1], [y1], psym=2
            keyboard_select = 1
        ENDELSE
    ENDIF

    IF keyboard_select THEN BEGIN
        PRINT, 'Please enter the x-coord for the opposite corner: '
        READ, x2, PROMPT='[waiting for input]> '
        PRINT, 'Please enter the y-coord for the opposite corner: '
        READ, y2, PROMPT='[waiting for input]> '
        x2 = FIX(ROUND(x2))
        y2 = FIX(ROUND(y2))

        IF NOT keyword_set(noopen) THEN BEGIN
            window, xs=xsize, ys=ysize, /pixmap, /free
            pixid = !d.window
            tvim, ex_img
            WSET, window_var
            DEVICE, copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
            PLOTS, [x1,x1,x2,x2,x1], [y1,y2,y2,y1,y1]
        ENDIF
    ENDIF

    ;Cacluate and print region size
    x_diff = ABS(x2-x1)
    y_diff = ABS(y2-y1)
    PRINT, ''
    PRINT, 'Region dimensions = ['+strtrim(x_diff, 2)+', '+strtrim(y_diff, 2)+']'

    ;-----------------------------------
    ;plot preview of the selected region
    IF NOT keyword_set(noopen) THEN BEGIN
        x_left = MAX([MIN([x1, x2]), 0])
        x_right = MIN([MAX([x1, x2]), nx_data - 1])
        y_bot = MAX([MIN([y1, y2]), 0])
        y_top = MIN([MAX([y1, y2]), ny_data - 1])

        ;Slicing example image (already includes scaling)
        ex_region = ex_img[x_left:x_right, y_bot:y_top]

        ;Set plotting elements
        ;writing to buffer pixmap reduces on-screen jitter!
        IF n_elements(preview_win_id) EQ 0 THEN BEGIN
            determine_window_size, xsize=xsize, ysize=ysize, openwin=openwin, /new_win
            window, openwin[0], xs=xsize, ys=ysize
            preview_win_id = openwin[0]
        ENDIF ELSE BEGIN
            wset, preview_win_id
        ENDELSE
        window, xs=xsize, ys=ysize, /pixmap, /free
        pixid = !d.window
        tvim, ex_region
        wset, preview_win_id
        device, copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
    ENDIF

    IF NOT KEYWORD_SET(noopen) THEN BEGIN
        ;Ask user if they want to try again or stop
        PRINT, ''
        PRINT, 'Type "e" to extract the current region'
        READ, ask_select, PROMPT='[or hit ENTER to select a different region]: '
        ask_select = strlowcase(ask_select)
        IF ask_select.startswith('e') THEN loop_select = 0 ;ends the loop
    ENDIF ELSE BEGIN
        loop_select = 0
    ENDELSE
ENDWHILE

;###############################################################################
;VALIDATING REGION COORDS AND RESETTING IF NEEDED (prevents unexpected behavior)
;###############################################################################
IF x1 LT 0 THEN x1 = 0
IF x1 GT nx_data THEN x1 = nx_data - 1
IF x2 LT 0 THEN x2 = 0
IF x2 GT nx_data THEN x2 = nx_data - 1

IF y1 LT 0 THEN y1 = 0
IF y1 GT ny_data THEN y1 = ny_data - 1
IF y2 LT 0 THEN y2 = 0
IF y2 GT ny_data THEN y2 = ny_data - 1

;###############################################################################
;EXTRACT THE DATA REGION
;###############################################################################
x_left = MIN([x1, x2])
x_right = MAX([x1, x2])
y_bot = MIN([y1, y2])
y_top = MAX([y1, y2])

;Slicing the input data cube
out_data = in_data[x_left:x_right, y_bot:y_top, *]

;Outputting slit coordinates (might be used by other scripts?)
outvec = [x_left, y_bot, x_right, y_top]

END
