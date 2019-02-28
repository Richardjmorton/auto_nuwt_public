;+
;NAME: NUWT_ARC_SLIT
;
;PURPOSE:
;   Extracts a time-distance slice along a circular arc in a series of images.
;   Can operate in an interactive mode. Basic code framework was adapted from
;   diag_slit.pro
;
;EXAMPLE CALL:
;   nuwt_arc_slit, in_data, out_slit, center=[500, 500], radius=200
;
;INPUTS:
;   in_data - 2D or 3D array containing image data
;
;OPTIONAL INPUTS:
;   index - index structure containing the .fits file header information
;           (currently only accepts SDO/AIA index structures). If provided,
;           many parameters (solar radius, center, image res., ect...) will be
;           extracted from the index structure.
;   center - [x,y] center pixel coordinates of the circle that defines the arc
;   radius - radius of the circle from which the arc slit will be selected.
;            By default, is assumed to be in units of [pixels]
;   start_ang - start angle in degrees from the positive x-direction
;   subtend - total angle (in degrees) for the slection arc to subtend. If set to
;             360 [deg] or greater, the program will return all points on the
;             circle that fall within the given image. Note: a negative subtend
;             anglue will extract data in a clockwise direction
;   ang_res - angular resolution (in degrees) for sampling the image.
;             WARNING! if not careful, it is easy to under- or over-sample the
;             image. The default method is scaled based on the circumference of
;             the given circle and should sample each point only once.
;   /clockwise - if set, will sample points in a clockwise direction.
;   /rel_limb - if set, will assume the requested radius is given relative to the
;               solar limb. The user MUST also supply either the solar radius or
;               an index structure containing the solar radius.
;   solar_radius - radius of the Sun Used in conjunction with /rel_limb. By defualt
;                  is assumed to be in units of [pixels].
;   units_radius - defained what units the input radii (arc and solar) are given
;                  in. Choose from "pixels" (default), "arcsec", "km", or "Mm".
;   res - spatial resolution of data (in [arcsec]). Default is 0
;   cad - temporal cadence of the data (in [s]). Default is 0
;   km_per_arcsec - ratio of [km]/[arcsec]. Defaults to 725.27 which is the
;                   mean scale of the solar surface.
;   fn - frame number to plot
;   /log_plot - if set, will log-scale the example image intensity values
;               OUTPUT VALUES WILL BE UNEFFECTED!
;   xsize, ysize - size of graphics window (command for widgets)
;   win_id - set to value of open graphics window you want to use. If not set,
;            opens next free window.
;   /noopen - if set, suppresses the opening of plot windows
;
;OUTPUTS:
;   out_slit - output array where data along the arc is stored
;   meta_out - structure containing various metadata about the extracted slit
;   outvec - array of [h, k, radius, start_ang, subtend] parameters used for the circle
;   slit_coords - array of shape [T, 2] containing the x and y coordinates of
;                 virtual slit used for interpolating and extracting the data
;
;HISTORY: Name---------Date---------Description
;         M Weberg  8 JUNE, 2016  Initial coding
;         M Weberg  9 JUNE, 2016  Large upgrade, added options for:
;                                  - variable angular resolution
;                                  - start and subtending angle
;                                  - direction of data sampling
;         M Weberg   1 DEC, 2016  Added option to average over multiple slices
;                                 in order to improve signal/noise ratio
;         M Weberg  29 MAR, 2017  Now stores basic metadata about slit coords and
;                                 input parameters in a COMMON block. It can even
;                                 read and record basic image details if given an
;                                 associated "index" array of AIA FTIS headers.
;         M Weberg  1  MAR, 2018  Updated the 'outvec' values to give more information
;
;TO DO:
;   - NONE (at the moment)
;-

PRO NUWT_ARC_SLIT, in_data, out_slit, index=index, $
                   center=center, radius=radius, start_ang=start_ang, $
                   subtend=subtend, width=width, ang_res=ang_res, clockwise=clockwise, $
                   rel_limb=rel_limb, solar_radius=solar_radius, units_radius=units_radius, $

                   res=res, cad=cad, km_per_arcsec=km_per_arcsec, $
                   wavelength=wavelength, notes=notes, $

                   meta_out=meta_out, slit_only=slit_only, $
                   slit_coords=slit_coords, outvec=outvec, $
                   fn=fn, log_plot=log_plot, $
                   xsize=xsize, ysize=ysize, win_id=win_id, noopen=noopen

COMPILE_OPT IDL2
;@int_meta######################################################################
;INITIALIZE METADATA STRUCTURE AND COMMON BLOCK
;###############################################################################
COMMON slit_meta_dat, slit_meta

IF NOT KEYWORD_SET(wavelength) THEN wavelength = -1
IF KEYWORD_SET(notes) THEN user_notes = STRING(notes)+', ' ELSE user_notes = ''
slit_meta = {slit_type:'arc', date_extracted:systime(), $
             dataset:'unknown', wavelength:wavelength, $
             source_size:[0,0,0], $
             start_date:'unknown', end_date:'unknown', $
             notes:user_notes, $
             res:0.0d, cad:0.0d, km_per_arcsec:0.0d, $
             ref_center:[0.0,0.0], r_sun:0.0, $
             center:[0.0,0.0], radius:0.0, $
             start_ang:0.0, subtend:0.0, $
             ang_res:0.0, width:1, clockwise:0, $
             start_coord:[0.0,0.0], end_coord:[0.0,0.0]}

IF KEYWORD_SET(rel_limb) THEN name_of_center = 'Solar limb' ELSE name_of_center = 'user inputed center'
IF KEYWORD_SET(clockwise) THEN clockwise = 1 ELSE clockwise = 0

;@set_units#####################################################################
;UNIT CONVERSIONS AND TAGGING
;###############################################################################
IF NOT KEYWORD_SET(res) THEN res = 0d
IF NOT KEYWORD_SET(cad) THEN cad = 0d
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27

IF NOT KEYWORD_SET(solar_radius) THEN solar_radius = 0.0
IF NOT KEYWORD_SET(units_radius) THEN units_radius = 'pixels'

;convert user inputted solar radius to units of [pixels]
zero_res = 0 ;this is a flag used to prevent math errors while retaining a value of 0 when no res is given
IF res EQ 0d THEN BEGIN
    ;prevents math errors
    res = 1d
    zero_res = 1
ENDIF
IF STRLOWCASE(units_radius) EQ 'arcsec' THEN BEGIN
    pxl_solar_radius = solar_radius/res
ENDIF ELSE IF STRLOWCASE(units_radius) EQ 'km' THEN BEGIN
    pxl_solar_radius = (solar_radius)/(km_per_arcsec*res)
ENDIF ELSE IF STRLOWCASE(units_radius) EQ 'Mm' THEN BEGIN
    pxl_solar_radius = (solar_radius*1000.0)/(km_per_arcsec*res)
ENDIF ELSE BEGIN
    pxl_solar_radius = solar_radius[0]
ENDELSE

IF zero_res EQ 1 THEN res = 0d

;@get_meta######################################################################
;IF GIVEN AN INDEX STRUCTURE WITH IMAGE METADATA, LOOKUP VALUES AS NEEDED
;###############################################################################
IF KEYWORD_SET(index) THEN BEGIN
    index_tag_names = TAG_NAMES(index[0])
    test_for_telescope_tag = WHERE(STRUPCASE(index_tag_names) EQ 'TELESCOP')
    IF test_for_telescope_tag GE 0 THEN BEGIN
        IF index[0].TELESCOP.StartsWith('SDO', /FOLD_CASE) THEN BEGIN
            ;SDO / AIA data dectected!
            slit_meta.dataset = index[0].TELESCOP+'_'+index[0].INSTRUME+'_lvl_'+string(index[0].LVL_NUM, format='(f4.2)')
            slit_meta.wavelength = index[0].WAVELNTH
            slit_meta.start_date = index[0].DATE_OBS
            slit_meta.end_date = index[-1].DATE_OBS
            slit_meta.ref_center = [(index[0].CRPIX1), (index[0].CRPIX2)]
            res = index[0].CDELT1
            cad = CALC_MEAN_DATA_CADENCE(index)
            km_per_arcsec = (index[0].RSUN_REF/index[0].RSUN_OBS)/1000.0
            pxl_solar_radius = index[0].R_SUN ;units of [pixels]
            IF N_ELEMENTS(center) LE 1 THEN BEGIN
                ;note: the AIA FITS coordinates system labels the lower left corner
                ;of a full-disk image as [1,1] not [0,0] as in IDL (hence the shift)
                center = [(index[0].CRPIX1 - 1), (index[0].CRPIX2 - 1)]
                IF NOT KEYWORD_SET(rel_limb) THEN name_of_center = 'solar disk center'
            ENDIF
        ENDIF ELSE BEGIN
            PRINT, 'WARNING: unknown telescope. Index structure will be ignored.'
        ENDELSE
    ENDIF ELSE BEGIN
        PRINT, 'WARNING: unknown index format'
    ENDELSE
ENDIF

;@check_inputs##################################################################
;CHECKING DIMENSIONS OF INPUT DATA
;###############################################################################
sz = size(in_data)
IF sz[0] EQ 3 THEN BEGIN
    nx_data = sz[1]
    ny_data = sz[2]
    nt_data = sz[3]
ENDIF ELSE BEGIN
    IF sz[0] EQ 2 THEN BEGIN
        PRINT, 'WARNING: input data is only 2D! Output array will be 1D.'
        nx_data = sz[1]
        ny_data = sz[2]
        nt_data = 1
        fn = 0
    ENDIF ELSE BEGIN
        MESSAGE, 'ERROR: please input a valid 2D or 3D array'
    ENDELSE
ENDELSE

; IF KEYWORD_SET(errors) THEN BEGIN
;     err_sz = SIZE(errors)
;     IF err_sz[0] NE sz[0] OR TOTAL(err_sz NE sz) GE 1 THEN BEGIN
;         MESSAGE, 'ERROR: "errors" array must have the same dimensions as the input data!'
;     ENDIF
; ENDIF


;@plot_ex########################################################################
;PLOT EXAMPLE IMAGE FOR REFERENCE AND/OR INTERACTIVE SELECTION
;###############################################################################
IF NOT keyword_set(noopen) THEN BEGIN
    IF n_elements(fn) eq 0 THEN fn=0
    IF fn GT (nt_data-1) THEN fn = nt_data-1

    ex_img = in_data[*,*,fn]
    IF KEYWORD_SET(log_plot) THEN BEGIN
        min_ex = MIN(ex_img)
        ex_img = ex_img + ABS(min_ex) + 1
        ex_img = ALOG10(ex_img)
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

;@query#########################################################################
;QUERY USER FOR SLIT PARAMETERS AS NEEDED
;###############################################################################
IF n_elements(center) LE 1 THEN BEGIN
    ;Ask for the center of circle and radius (if not set)
    PRINT, ''
    PRINT, 'Please enter x-coord of the circle center: '
    READ, cx, PROMPT='[waiting for input]> '
    PRINT, PROMPT='Please enter y-coord of the circle center: '
    READ, cy, PROMPT='[waiting for input]> '
    center = [cx,cy]
ENDIF

IF NOT keyword_set(noopen) THEN plots, center[0], center[1], psym=2

IF n_elements(radius) EQ 0 OR n_elements(start_ang) EQ 0 THEN BEGIN
    ;Ask for the radius and start angle of the slit (if not set)
    PRINT, ''
    IF NOT keyword_set(noopen) THEN BEGIN
        PRINT, 'Please select the start point of the arc slit on the image'
        CURSOR, x, y, /Data, /Down
        plots, x, y, psym=7
        rel_x = x - center[0]
        rel_y = y - center[1]
        radius = SQRT(rel_x^2 + rel_y^2)
        IF rel_x EQ 0 THEN BEGIN
            IF rel_y GE 0 THEN start_ang = 90.0d
            IF rel_y LT 0 THEN start_ang = 270.0d
        ENDIF ELSE BEGIN
            quad_corr = 0.0d ;QI
            IF rel_x LT 0 THEN quad_corr = 180.0 ;QII or QIII
            IF rel_y LT 0 AND rel_x GT 0 THEN quad_corr = 360.0 ;QIV
            start_ang = !CONST.RtoD*ATAN(rel_y/rel_x) + quad_corr
        ENDELSE
        PRINT, 'x = '+strtrim(x, 2)+',  y = '+strtrim(y, 2)
        PRINT, 'radius = '+strtrim(radius, 2)+' [pixels],  start_ang = '+strtrim(start_ang, 2)+' [deg]'
        units_radius = 'pixels'
    ENDIF ELSE BEGIN
        IF n_elements(radius) EQ 0 THEN BEGIN
            PRINT, 'Please enter the radius of the arc slit: '
            READ, radius, PROMPT='[waiting for input]> '
        ENDIF
        IF n_elements(start_ang) EQ 0 THEN BEGIN
            PRINt, PROMPT='Please enter the start angle of the arc slit: '
            READ, radius, PROMPT='[waiting for input]> '
        ENDIF
    ENDELSE
ENDIF

IF n_elements(subtend) EQ 0 THEN BEGIN
    ;Ask for end point if subtend is not set
    PRINT, ''
    IF NOT keyword_set(noopen) THEN BEGIN
        PRINT, 'Please select a point along the same radial line as your desired end point'
        CURSOR, x, y, /Data, /Down
        plots, x, y, psym=7
        rel_x = x - center[0]
        rel_y = y - center[1]
        IF rel_x EQ 0 THEN BEGIN
            IF rel_y GE 0 THEN end_ang = 90.0d
            IF rel_y LT 0 THEN end_ang = 270.0d
        ENDIF ELSE BEGIN
            quad_corr = 0.0d ;QI
            IF rel_x LT 0 THEN quad_corr = 180.0 ;QII or QIII
            IF rel_y LT 0 AND rel_x GT 0 THEN quad_corr = 360.0 ;QIV
            end_ang = !CONST.RtoD*ATAN(rel_y/rel_x) + quad_corr

            ;always select the shortest arc, even if it spans the 0 deg axis
            IF (end_ang - start_ang) GT 180.0 THEN end_ang = end_ang - 360.0d
            IF (start_ang - end_ang) GT 180.0 THEN end_ang = end_ang + 360.0d
        ENDELSE
        subtend = (end_ang - start_ang)
        PRINT, 'x = '+strtrim(x, 2), '  y = '+strtrim(y, 2)
        PRINT, 'subtend = '+strtrim(subtend, 2)+' [deg]'
        IF subtend LT 0 THEN BEGIN
            PRINT, 'Note: the slit will be extracted in the clockwise direction'
        ENDIF ELSE BEGIN
            PRINT, 'Note: the slit will be extracted in the anticlockwise direction'
        ENDELSE
    ENDIF ELSE BEGIN
        PRINT, 'Please enter the number of degrees for the the arc to subtend: '
        PRINT, '(note: negative angles will be extracted in a clockwise direction) '
        READ, subtend, PROMPT='[waiting for input]> '
    ENDELSE
ENDIF

IF NOT keyword_set(width) THEN BEGIN
    ;Ask for slit width (if not set)
    PRINT, ''
    PRINT, 'Please enter the slit pixel width to average over: '
    READ, width, PROMPT='[waiting for input]> '
ENDIF

input_radius = radius ;for saving the input radius to the slit_meta.notes string

;@validate######################################################################
;VALIDATING SLIT VALUES AND RESETTING IF NEEDED (prevents unexpected behavior)
;###############################################################################

;Convert input radius to units of [pixels] (unless selected from the image)
zero_res = 0 ;this is a flag used to prevent math errors while retaining a value of 0 when no res is given
IF res EQ 0d THEN BEGIN
    ;prevents math errors
    res = 1d
    zero_res = 1
ENDIF
IF STRLOWCASE(units_radius) EQ 'arcsec' THEN BEGIN
    pxl_radius = radius/res
ENDIF ELSE IF units_radius EQ 'km' THEN BEGIN
    pxl_radius = (radius)/(km_per_arcsec*res)
ENDIF ELSE IF units_radius EQ 'Mm' THEN BEGIN
    pxl_radius = (radius*1000.0)/(km_per_arcsec*res)
ENDIF ELSE BEGIN
    ;units of [pixels]
    pxl_radius = radius[0] ; ensures no backward leakage into the user inputs
ENDELSE

IF zero_res EQ 1 THEN res = 0d

IF STRLOWCASE(units_radius) NE 'pixels' THEN BEGIN
    PRINT, ''
    PRINT, 'Note: input radius of '+STRTRIM(input_radius, 2)+' '+units_radius+' equates to '+STRTRIM(pxl_radius, 2)+' pixels'
ENDIF

IF KEYWORD_SET(rel_limb) THEN BEGIN
    ;input radius was given relative to the solar limb
    pxl_radius = pxl_solar_radius + pxl_radius
ENDIF

IF pxl_radius LT 0.0 THEN pxl_radius = -1*pxl_radius
IF width LT 1 THEN width = 1
IF width GT ny_data THEN width = ny_data
half_width = FIX(FLOOR(width/2.0)) ;i.e. number of data points to select on either side
true_width = 2*half_width + 1

IF subtend LT 0 THEN clockwise = 1 ELSE clockwise = 0
IF subtend EQ 0.0 THEN subtend = 360.0
IF subtend GT 360.0 THEN subtend = 360.0
IF KEYWORD_SET(clockwise) AND subtend GT 0 THEN subtend = -1.0*subtend

;@gen_coords####################################################################
;GENERATE SLIT COORDINATES
;###############################################################################
;Determine angular resolution needed for a critical sampling rate
IF n_elements(ang_res) EQ 0 THEN BEGIN
    circ = 2.0*!CONST.PI*pxl_radius ;circumference in units of pixels
    ang_res = 360.0/circ
ENDIF

;Generate the array of angles (double precision in units of degrees)
;note: these angles are always positive and are measured from the pos x-direction
IF keyword_set(clockwise) THEN BEGIN
    ;note: the extra +1 ensures the angles START at the correct ang after reversal
    theta = reverse((dindgen(360d/ang_res) + 1)*ang_res + start_ang)
    rotation_dir= -1.0
ENDIF ELSE BEGIN
    theta = dindgen(360d/ang_res)*ang_res + start_ang
    rotation_dir = 1.0
ENDELSE

;Trim the angles down to the desired arc
IF abs(subtend) LT 360.0 THEN BEGIN
    begin_angle = theta[0]
    end_angle = begin_angle + rotation_dir*double(abs(subtend))
    IF keyword_set(clockwise) THEN BEGIN
        loc_arc = where(theta LE begin_angle and theta GE end_angle)
    ENDIF ELSE BEGIN
        loc_arc = where(theta GE begin_angle and theta LE end_angle)
    ENDELSE
    theta = theta[loc_arc]
ENDIF

num_ang = n_elements(theta)
offset_vals = indgen(true_width) - half_width
x_vals = fltarr(num_ang, true_width)
y_vals = fltarr(num_ang, true_width)

;Generate the x and y coordinates on the circle
FOR s=0, true_width-1 DO BEGIN
    x_vals[*,s] = center[0] + (pxl_radius + offset_vals[s])*cos(theta*!CONST.DtoR)
    y_vals[*,s] = center[1] + (pxl_radius + offset_vals[s])*sin(theta*!CONST.DtoR)
ENDFOR

;@clip_coords###################################################################
;FILTER AND CLIP SLIT COORDS
;###############################################################################
;Filters slit coords for only the region inside the image (based on the center of the slit)
loc_valid = where(x_vals[*,half_width] GE 0 and x_vals[*,half_width] LE nx_data-1 and $
                  y_vals[*,half_width] GE 0 and y_vals[*,half_width] LE ny_data-1, /NULL)
num_valid = n_elements(loc_valid)
IF num_valid GT 0 THEN BEGIN
    IF num_valid LT num_ang THEN BEGIN
        ;only bother slicing data if there is a section that falls outside
        x_vals = x_vals[loc_valid[0]:loc_valid[-1], *]
        y_vals = y_vals[loc_valid[0]:loc_valid[-1], *]
    ENDIF
ENDIF ELSE BEGIN
    MESSAGE, 'ERROR: the selected arc lies completely outside of the image!'
ENDELSE
slit_lnth = n_elements(x_vals[*,0])

;Clips edges of the slit where the center is still valid
;i.e. any part of a slit that runs off the image will reset to the nearest edge
;Hmm... may need to reconsider the need for this code bit... the filtering above
;suggests that this should never actually be used...at the smae time, there is
;still an issue with the finite width of the slit running off the image...
x_vals = x_vals > 0
x_vals = x_vals < (nx_data-1)
y_vals = y_vals > 0
y_vals = y_vals < (ny_data-1)

;Plot the selected (and validated) slit on the image
IF NOT keyword_set(noopen) THEN plots, x_vals[*,half_width], y_vals[*,half_width]

;@get_slit######################################################################
;EXTRACT THE DATA SLIT
;###############################################################################
mindat = MIN(in_data)
IF mindat GT 0 THEN mindat = -1.0 ELSE mindat = mindat - 0.1*STDDEV(in_data)

;Slicing the input data cube
slit_data = fltarr(slit_lnth, nt_data, true_width)
FOR s=0, (true_width-1) DO BEGIN
    FOR i=0L, (nt_data-1) DO BEGIN
        slit_data[*,i,s] = interpolate(reform(in_data[*,*,i]), x_vals[*,s], y_vals[*,s], cubic=-0.5, missing=mindat)
    ENDFOR
ENDFOR

;Averaging over the width (if needed)
IF true_width GT 1 THEN BEGIN
    slit_data = MEAN(slit_data, dimension=3)
ENDIF ELSE BEGIN
    slit_data = slit_data[*,*,0]
ENDELSE

;@save_meta#####################################################################
;SAVING INFORMATION TO THE SLIT_META STRUCTURE (AS WELL AS OPTIONAL OUTPUTS)
;###############################################################################
slit_meta.source_size = [nx_data, ny_data, nt_data]
slit_meta.notes = slit_meta.notes+STRTRIM(input_radius, 2)+' '+units_radius+' from the '+name_of_center
slit_meta.res = res
slit_meta.cad = cad
slit_meta.km_per_arcsec = km_per_arcsec
slit_meta.r_sun = pxl_solar_radius
slit_meta.center = center
slit_meta.radius = pxl_radius
slit_meta.start_ang = start_ang
slit_meta.subtend = subtend
slit_meta.width = width
slit_meta.ang_res = ang_res
slit_meta.clockwise = clockwise ;0 for false, 1 for true
slit_meta.start_coord = [x_vals[0,half_width],y_vals[0,half_width]]
slit_meta.end_coord = [x_vals[1,half_width],y_vals[1,half_width]]

;MAIN OUTPUT
IF KEYWORD_SET(slit_only) THEN BEGIN
    out_slit = slit_data
ENDIF ELSE BEGIN
    out_slit = {slit:slit_data, meta:slit_meta}
ENDELSE

;Outputting slit_meta as a seperate variable
meta_out = slit_meta

;Outputting circle parameters (might be used by other scripts?)
outvec = [center[0], center[1], pxl_radius, start_ang, subtend]

;Output the coords of the slit center
arc_coords = fltarr(slit_lnth, 2)
arc_coords[*,0] = x_vals[*,half_width]
arc_coords[*,1] = y_vals[*,half_width]
slit_coords = arc_coords
END
