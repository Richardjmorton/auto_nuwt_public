;+
;NAME: PLOT_NUWT_SLIT_LOC
;
;PURPOSE:
;   Plots an example solar image with the specified dataslit(s) overplotted.
;   For ease of plotting, can accept slit_meta structures outputted by
;   "nuwt_arc_slit.pro" or "nuwt_diag_slit.pro". Can even plot a combination of slit types.
;
;INPUTS:
;   in_data - 2D or 3D array containing image data
;
;OPTIONAL INPUTS:
;   res - spatial resolution of the data (in [arcsec]). If not given by the user,
;         the program will load the 'res' value stored by NUWT for the current
;         set of results. If NUWT does not have a valid res value, then all
;         distance values will revert to units of [pixels]
;   cad - temporal cadence of the data (in [s]). If not given by the user, the
;         program will load the 'cad' value stored by NUWT for the current set of
;         results. If NUWT does not have a valid cad value, then all time values
;         will revert to units of [timesteps]
;   km_per_arcsec - ratio of [km]/[arcsec]. Defaults to 725.27 which is the
;                   mean scale of the solar surface as seen from 1 AU.
;   dist_units - string indicating what units to use on the distance axis of plot
;                Defaults to units of 'arcsec' if 'res' is known or 'pixels' if
;                'res' is unknown.
;   time_units - string indicating what units to use on the time axis of plots.
;                Defaults to units of 's' if 'cad' is known or 'timesteps' if
;                'cad' is unknown.
;   wavelength - wavelength (in [angstrom]) of input SDO / AIA images. Used to
;                load the correct color table from SolarSoft.
;   slit_meta - metadata structure outputed by "nuwt_arc_slit.pro" or "nuwt_diag_slit.pro".
;               Should contain all the information required to draw the slit.
;   center - (arc slits) [x,y] center pixel coordinates of the circle
;   radius - (arc slits) radius of the circle from which the arc slit will be
;            selected. By default, is assumed to be in units of [pixels]
;   start_ang - (arc slits) start angle in degrees from the positive x-direction
;   subtend - (arc slits) total angle (in degrees) for the slection arc to subtend.
;             If not defined, the program will return all points on the circle that
;             fall within the given image. Note: sign is ignored - i.e. the
;             actual value used is abs(subtend)
;   ang_res - (arc slits) angular resolution (in degrees) for sampling the image.
;             WARNING! if not careful, it is easy to under- or over-sample the
;             image. The default method is scaled based on the circumference of
;             the given circle and should sample each point only once.
;   /clockwise - (arc slits) if set, will sample points in a clockwise direction.
;   x1,x2,y1,y2 - (diag slits) X & Y coordinates at either end of the diagonal line
;   fn - frame number to plot
;
;OUTPUTS:
;   out_data - output array where data along diagonal line is stored
;   outvec - array of [h, k, r] parameters used for the circle
;   slit_coords - array of shape [T, 2] containing the x and y coordinates of
;                 virtual slit used for interpolating and extracting the data
;
;HISTORY: Name---------Date---------Description
;         M Weberg  8 MAY, 2017  Initial coding
;
;TO DO:
;   - NONE (at the moment)
;-

PRO PLOT_NUWT_SLIT_LOC, in_data, res=res, km_per_arcsec=km_per_arcsec, $
                        wavelength=wavelength, slit_meta=slit_meta, $
                        slit_type=slit_type, center=center, radius=radius, $
                        start_ang=start_ang, subtend=subtend, $
                        width=width, ang_res=ang_res, clockwise=clockwise, $
                        x1=x1, x2=x2, y1=y1, y2=y2, $
                        slit_coords=slit_coords, $
                        close_region=close_region, end_to_start=end_to_start, $
                        title=title, filetag=filetag, $
                        cutout_coords=cutout_coords, thick=thick, $
                        fn=fn, log_plot=log_plot, screen=screen, scale=scale, $
                        dist_units=dist_units, plot_size=plot_size

IF NOT KEYWORD_SET(title) THEN title = ''
IF NOT KEYWORD_SET(filetag) THEN filetag='' ELSE filetag = filetag+'_'

;@check_inputs##################################################################
;CHECKING DIMENSIONS OF INPUT DATA
;###############################################################################
input_sz = size(in_data)
IF input_sz[0] EQ 3 THEN BEGIN
    nx_data = input_sz[1]
    ny_data = input_sz[2]
    nt_data = input_sz[3]
ENDIF ELSE BEGIN
    IF input_sz[0] EQ 2 THEN BEGIN
        nx_data = input_sz[1]
        ny_data = input_sz[2]
        nt_data = 1
        fn = 0
    ENDIF ELSE BEGIN
        MESSAGE, 'ERROR: please input a valid 2D or 3D image array'
    ENDELSE
ENDELSE

IF N_ELEMENTS(fn) EQ 0 THEN fn = 0
IF N_ELEMENTS(fn) GT 1 THEN fn = fn[0]
IF fn GT (nt_data-1) THEN fn = nt_data-1

IF N_ELEMENTS(cutout_coords) EQ 0 THEN BEGIN
    ;[default] No cutout
    ex_img = in_data[*,*,fn]
    x_offset = 0
    y_offset = 0
ENDIF ELSE IF N_ELEMENTS(cutout_coords) NE 4 THEN BEGIN
    MESSAGE, 'ERROR: "cutout_coords" array is invalid. Please input an array of the form [X, Y, xsize, ysize]'
ENDIF ELSE BEGIN
    x_offset = cutout_coords[0]
    cutout_end_x = cutout_coords[0] + (cutout_coords[2] - 1)
    y_offset = cutout_coords[1]
    cutout_end_y = cutout_coords[1] + (cutout_coords[3] - 1)
    ex_img = in_data[x_offset:cutout_end_x, y_offset:cutout_end_y, fn]
    cutout_sz = size(ex_img)
    nx_data = cutout_sz[1]
    ny_data = cutout_sz[2]
ENDELSE

IF KEYWORD_SET(log_plot) THEN BEGIN
    min_ex = MIN(ex_img)
    ex_img = ex_img + ABS(min_ex) + 1
    ex_img = ALOG10(ex_img)
ENDIF

; A handful of control and ref values
IF NOT KEYWORD_SET(scale) THEN scale = 1.5;2.0
IF KEYWORD_SET(thick) THEN slit_line_thickness = thick ELSE slit_line_thickness = 5;3
slit_color = 'white'
im_center = [nx_data/2.0, ny_data/2.0] ;already adjusted for cutouts

;@read_meta#####################################################################
;IF GIVEN A STRUCTURE WITH SLIT METADATA, LOOKUP NEEDED VALUES
;###############################################################################
IF KEYWORD_SET(slit_meta) THEN BEGIN
    num_slits = n_elements(slit_meta)
    res = slit_meta[0].res
    km_per_arcsec = slit_meta[0].km_per_arcsec
    IF N_ELEMENTS(wavelength) EQ 0 THEN wavelength = slit_meta[0].wavelength ;can now override the plot color scale

    ;Initialize array of slit coords
    slit_type = strarr(num_slits)
    all_centers = fltarr(2, num_slits)
    all_radii = fltarr(num_slits)
    all_start_ang = fltarr(num_slits) + 180.0
    all_subtend = fltarr(num_slits)
    all_x1 = fltarr(num_slits)
    all_y1 = fltarr(num_slits)
    all_x2 = fltarr(num_slits)
    all_y2 = fltarr(num_slits)
    FOR s=0, (num_slits-1) DO BEGIN
        slit_type[s] = slit_meta[s].slit_type
        IF strlowcase(slit_type[s]) EQ 'arc' THEN BEGIN
            all_centers[*,s] = slit_meta[s].center - [x_offset, y_offset]
            all_radii[s] = slit_meta[s].radius
            all_start_ang[s] = slit_meta[s].start_ang
            all_subtend[s] = slit_meta[s].subtend
            ; ang_res = slit_meta[s].ang_res
            clockwise = slit_meta[s].clockwise
            ; clockwise = 0
        ENDIF ELSE IF strlowcase(slit_type[s]) EQ 'diag' THEN BEGIN
            all_centers = slit_meta.ref_center - [x_offset, y_offset]
            all_x1 = slit_meta[s].start_coord[0] - x_offset
            all_y1 = slit_meta[s].start_coord[1] - y_offset
            all_x2 = slit_meta[s].end_coord[0] - x_offset
            all_y2 = slit_meta[s].end_coord[1] - y_offset
        ENDIF
    ENDFOR
ENDIF

;@units#########################################################################
;UNIT CONVERSIONS AND LABELS
;###############################################################################
IF KEYWORD_SET(dist_units) THEN plot_dist_units = dist_units ELSE plot_dist_units = 'arcsec'

dx = 1.0
IF NOT KEYWORD_SET(km_per_arcsec) THEN km_per_arcsec = 725.27
IF NOT KEYWORD_SET(res) THEN BEGIN
    res = 1.0
    plot_dist_units = 'pixels'
ENDIF

IF plot_dist_units EQ 'arcsec' THEN dx = res
IF plot_dist_units EQ 'm' THEN dx = res*km_per_arcsec*1000.0
IF plot_dist_units EQ 'km' THEN dx = res*km_per_arcsec
IF plot_dist_units EQ 'Mm' THEN dx = res*km_per_arcsec/1000.0

;@read_other####################################################################
;PARSE OTHER INPUTS IF SLIT METADATA IS NOT PROVIDED (AND SET DEFAULTS)
;###############################################################################
IF NOT KEYWORD_SET(slit_meta) THEN BEGIN
    num_arc = MAX([n_elements(radius), n_elements(start_ang)])
    num_diag = MAX([n_elements(x1), n_elements(y1), n_elements(x2), n_elements(y2)])
    num_slits = num_arc + num_diag

    ;Initialize array of slit coords
    slit_type = strarr(num_slits) + 'diag'
    all_centers = fltarr(2, num_slits)
    all_radii = fltarr(num_slits)
    all_start_ang = fltarr(num_slits)
    all_subtend = fltarr(num_slits)
    all_x1 = fltarr(num_slits)
    all_y1 = fltarr(num_slits)
    all_x2 = fltarr(num_slits)
    all_y2 = fltarr(num_slits)
    IF num_arc GT 0 THEN BEGIN
        slit_type[0:num_arc-1] = 'arc'
        IF n_elements(center) LT 2 THEN center = im_center
        FOR s=0, (num_arc-1) DO BEGIN
            IF n_elements(center) NE num_arc/2 THEN BEGIN
                ;same center for all arc slits
                all_centers[*,s] = center[*,0] - [x_offset, y_offset]
            ENDIF ELSE BEGIN
                ;different centers
                all_centers[*,s] = center[*,s] - [x_offset, y_offset]
            ENDELSE
        ENDFOR

        IF n_elements(radius) EQ 0 THEN radius = MIN(im_center)/2.0
        IF n_elements(radius) NE num_arc THEN BEGIN
            ;same radius for all arc slits
            all_radii[0:num_arc-1] = radius[0]
        ENDIF ELSE BEGIN
            ;different radii
            all_radii[0:num_arc-1] = radius
        ENDELSE

        IF n_elements(start_ang) EQ 0 THEN start_ang = 180.0
        IF n_elements(start_ang) NE num_arc THEN BEGIN
            ;same start angles for all arc slits
            all_start_ang[0:num_arc-1] = start_ang[0]
        ENDIF ELSE BEGIN
            ;differnt start angles
            all_start_ang[0:num_arc-1] = start_ang
        ENDELSE

        IF n_elements(subtend) EQ 0 THEN subtend = 0.0
        IF n_elements(subtend) NE num_arc THEN BEGIN
            ;same angular length for all arc slits
            all_subtend[0:num_arc-1] = subtend[0]
        ENDIF ELSE BEGIN
            ;differnt angular length
            all_subtend[0:num_arc-1] = subtend
        ENDELSE
    ENDIF
    ;DIAG SLITS CURRENTLY CANNOT BE LOADED VIA DIRECT INPUTS! USE SLIT_META INSTEAD!
ENDIF

;@gen_coords####################################################################
;GENERATE SLIT COORDINATES
;###############################################################################
all_x_vals = LIST()
all_y_vals = LIST()
FOR s=0, (num_slits-1) DO BEGIN
    IF strlowcase(slit_type[s]) EQ 'arc' THEN BEGIN
        ; ##### ARC SLIT
        slit_center = all_centers[*,s]
        slit_radius = all_radii[s]
        slit_start_ang = all_start_ang[s]
        slit_subtend = all_subtend[s]

        ;Determine angular resolution
        IF n_elements(ang_res) EQ 0 THEN BEGIN
            circ = 2.0*!CONST.PI*slit_radius ;circumference in units of pixels
            ang_res = 360.0/circ
        ENDIF

        ;Generate the array of angles (double precision in units of degrees)
        ;note: these angles are always positive and are measured from the pos x-direction
        IF KEYWORD_SET(clockwise) THEN BEGIN
            ;note: the extra +1 ensures the angles START at the correct ang after reversal
            theta = reverse((dindgen(360d/ang_res) + 1)*ang_res + slit_start_ang)
            rotation_dir= -1.0
        ENDIF ELSE BEGIN
            theta = dindgen(360d/ang_res)*ang_res + slit_start_ang
            rotation_dir = 1.0
        ENDELSE

        ;Trim the angles down to the desired arc
        IF slit_subtend GT 0 THEN BEGIN
            slit_subtend = double(abs(slit_subtend))
            ;Ignore invalid subtending angles
            IF slit_subtend GT 0.0 AND slit_subtend LT 360.0 THEN BEGIN
                begin_angle = theta[0]
                end_angle = begin_angle + rotation_dir*slit_subtend
                IF KEYWORD_SET(clockwise) THEN BEGIN
                    loc_arc = where(theta LE begin_angle and theta GE end_angle)
                ENDIF ELSE BEGIN
                    loc_arc = where(theta GE begin_angle and theta LE end_angle)
                ENDELSE
                theta = theta[loc_arc]
            ENDIF
        ENDIF

        num_ang = n_elements(theta)
        slit_x_vals = fltarr(num_ang)
        slit_y_vals = fltarr(num_ang)

        ;Generate the x and y coordinates on the circle
        slit_x_vals = slit_center[0] + slit_radius*cos(theta*!CONST.DtoR)
        slit_y_vals = slit_center[1] + slit_radius*sin(theta*!CONST.DtoR)

        ;Filters slit coords for only the region inside the image (based on the center of the slit)
        loc_valid = WHERE(slit_x_vals GE 0 and slit_x_vals LE (nx_data-1) AND $
                          slit_y_vals GE 0 and slit_y_vals LE (ny_data-1), /NULL)
        num_valid = n_elements(loc_valid)
        IF num_valid GT 0 THEN BEGIN
            IF num_valid LT num_ang THEN BEGIN
                ;only bother slicing data if there is a section that falls outside
                slit_x_vals = slit_x_vals[loc_valid[0]:loc_valid[-1], *]
                slit_y_vals = slit_y_vals[loc_valid[0]:loc_valid[-1], *]
            ENDIF
        ENDIF ELSE BEGIN
            MESSAGE, 'ERROR: the selected arc lies completely outside of the image!'
        ENDELSE
        length = n_elements(slit_x_vals[*,0])

        ;Clips edges of the slit where the center is still valid
        ;i.e. any part of a slit that runs off the image will reset to the nearest edge
        slit_x_vals = slit_x_vals > 0
        slit_x_vals = slit_x_vals < (nx_data-1)
        slit_y_vals = slit_y_vals > 0
        slit_y_vals = slit_y_vals < (ny_data-1)

        ;Transfer to plotting array
        all_x_vals.add, slit_x_vals
        all_y_vals.add, slit_y_vals
    ENDIF ELSE IF strlowcase(slit_type[s]) EQ 'diag' THEN BEGIN
        ;##### DIAG SLIT
        ;Set default values
        slit_center = all_centers[*,s]
        slit_x1 = all_x1[s]
        slit_y1 = all_y1[s]
        slit_x2 = all_x2[s]
        slit_y2 = all_y2[s]

        ;resets invalid values to edges of the image
        IF slit_x1 LT 0 THEN slit_x1 = 0.0
        IF slit_x1 GE (nx_data) THEN slit_x1 = 1.0*nx_data - 1
        IF slit_y1 LT 0 THEN slit_y1 = 0.0
        IF slit_y1 GE (ny_data) THEN slit_y1 = 1.0*ny_data - 1
        IF slit_x2 LT 0 THEN slit_x2 = 0.0
        IF slit_x2 GE (nx_data) THEN slit_x2 = 1.0*nx_data - 1
        IF slit_y2 LT 0 THEN slit_y2 = 0.0
        IF slit_y2 GE (ny_data) THEN slit_y2 = 1.0*ny_data - 1

        ;Transfer to plotting array
        all_x_vals.add, [slit_x1,slit_x2]
        all_y_vals.add, [slit_y1,slit_y2]
    ENDIF
ENDFOR

;@close_region##################################################################
;ADD EXTRA SLITS TO CLOSE THE REGION BOUNDED BY THE FIRST AND LAST SLITS
;###############################################################################
IF KEYWORD_SET(close_region) THEN BEGIN
    num_slits = num_slits + 2
    first_x_start = all_x_vals[0, 0]
    first_y_start = all_y_vals[0, 0]
    first_x_end = all_x_vals[0, -1]
    first_y_end = all_y_vals[0, -1]
    last_x_start = all_x_vals[-1, 0]
    last_y_start = all_y_vals[-1, 0]
    last_x_end = all_x_vals[-1, -1]
    last_y_end = all_y_vals[-1, -1]
    IF KEYWORD_SET(end_to_start) THEN BEGIN
        ;Connect the START of one slit to the END of the other (and vice versa)
        all_x_vals.add, [first_x_start, last_x_end]
        all_y_vals.add, [first_y_start, last_y_end]
        all_x_vals.add, [last_x_start, first_x_end]
        all_y_vals.add, [last_y_start, first_y_end]
    ENDIF ELSE BEGIN
        ;[Default]
        ;Connected the like-ends (START to START and END to END)
        all_x_vals.add, [first_x_start, last_x_start]
        all_y_vals.add, [first_y_start, last_y_start]
        all_x_vals.add, [last_x_end, first_x_end]
        all_y_vals.add, [last_y_end, first_y_end]
    ENDELSE
ENDIF

;@find_origin###################################################################
;DETERMINE WHERE THE AXIS ORIGIN IS RELATIVE TO THE LOWER-LEFT CORNER OF THE PLOT
;###############################################################################
origin = [all_centers[0,0], all_centers[1,0]] ;units of [pixels]

;@color#########################################################################
;LOADING COLOR TABLES
;###############################################################################
IF KEYWORD_SET(wavelength) THEN BEGIN
    aia_lct, rr, gg, bb, wave=wavelength, /load
    img_color_table = [[rr], [gg], [bb]]
ENDIF ELSE BEGIN
    img_color_table = 0
ENDELSE

;@plot_slit#####################################################################
;PLOT EXAMPLE IMAGE AND DATA SLIT(S)
;###############################################################################

IF N_ELEMENTS(plot_size) EQ 2 THEN BEGIN
    fig_dimensions = plot_size
    plt_positions = [50, 50, plot_size[0]-50, plot_size[1]-50]
ENDIF ELSE IF KEYWORD_SET(screen) THEN BEGIN
    fig_dimensions = [(170+(nx_data/2))*scale, (110+(ny_data/2))*scale]
    plt_positions = [115*scale, 55*scale, (115+(nx_data/2))*scale, (55+(ny_data/2))*scale]
ENDIF ELSE BEGIN
    fig_dimensions = [(285+(nx_data/2))*scale, (165+(ny_data/2))*scale]
    plt_positions = [230*scale, 110*scale, (230+(nx_data/2))*scale, (110+(ny_data/2))*scale]
ENDELSE

IF KEYWORD_SET(screen) THEN BEGIN
    fig_im = window(name='fig_im', dimensions=fig_dimensions)
ENDIF ELSE BEGIN
    fig_im = window(name='fig_im', dimensions=fig_dimensions, /buffer)
ENDELSE

plt_im = IMAGE(ex_img, axis_style=0, position=plt_positions, $
               RGB_table=img_color_table, /DEVICE, /CURRENT, title=title, font_size=24*scale)
im_x_ax = AXIS('X', location=0.0, target=plt_im, coord_transform=[-origin[0]*dx, dx], tickdir=1, $
               title='Solar X ['+plot_dist_units+']', ticklen=0.02, tickfont_size=24*scale)
im_y_ax = AXIS('Y', location=0.0, target=plt_im, coord_transform=[-origin[1]*dx, dx], tickdir=1, $
               title='Solar Y ['+plot_dist_units+']', ticklen=0.02, tickfont_size=24*scale)

;plotting the data slit
FOR s=0, (num_slits-1) DO BEGIN
    plt_slit = plot(all_x_vals[s], all_y_vals[s], color=slit_color, symbol='none', linestyle='-', thick=slit_line_thickness, overplot=plt_im)
ENDFOR

; Saving the plot to png (default)
IF NOT KEYWORD_SET(screen) THEN BEGIN
    save_filename = filetag+'td_slit.png'
    ; fig_im.save, save_filename, width=900, bit_depth=24, resolution=300
    fig_im.save, save_filename, bit_depth=24, resolution=300
    print, 'Finished plotting source data and slit!'
ENDIF

END
