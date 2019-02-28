;PURPOSE:
;;;; ------------------------------------------------------------------
;;;| This file contains all the required routines to track a solar    |
;;;| feature by fitting a curve to the general position of the feature|
;;;| and checking the value around the fitted curve in time.          |
;;; ------------------------------------------------------------------
;
; PRO NUWT_SPLINE_SLIT, in_data, slits=slits, frames=frames,$
;                       pickn=pickn, distance=distance, anchor=anchor,$
;                       gpoints=gpoints, spacing=spacing, $
;                       oversample=oversample, plot_slit=plot_slit, $
;                       vector=vector, curvepoints=curvepoints, intim=intim, $
;                       slit_coords=slit_coords, outvec=outvec
;
;:POSITIONAL INPUTS:
;   in_data - input datacube. Can be a 2D or 3D array.
;
;:KEYWORD INPUTS:
;   distance - Half length of each slit in pixels. Default is 10 pixels
;              (i.e. each slit is 21 total pixels long)
;
;   frames - Start and end image frames to use for slit selection (e.g. [t1, t2]).
;            The program will sum together all images from t1 to t2 and plot the
;            resulting image. This gives the user an idea of how the structures
;            in the image may change over time. Default is to sum the entire dataset.
;
;   /pickn - If set, the user will be asked to pick the guide points for the
;            desired spline by left clicking on an example image. The user may pick
;            as many guide points as they would like and must right click to end
;            picking process. /pickn will be automatiaclly enabled if "gpoints"
;            is unspecified.
;
;   gpoints  - Guide points for/from spline fit. Can be used as either an INPUT
;              or an OUTPUT,
;              INPUT - [N,2] array of guide points that will be used for the
;                      spline fit. If /pickn is also set, then any values in
;                      "gpoints" will be ignored and overwirtten with the values
;                      selected by the user.
;              OUTPUT - Returns the values seleted by the user
;
;   vector - Overplots vectors normal to the centre spline curve at each slit location.
;            In other words, plots vectors indicating the direction in which each
;            slit is extracted.
;
;   spacing - Scalar value that gives the spacing (in pixels) between points
;             on the fitted spline curve (i.e. the spacing between splits).
;             If not specified, the default spacing will be used in accordance
;             with the IDL routine 'spline_p', which gives ~8 interpolated points
;             between each input guide point. Note: spline_p intervals typically
;             do not give the desired spacing value (normally smaller values)
;             To get intervals close to desired spacing, use oversample
;
;   oversample - oversamples the curve in spline_p by 10^oversample. Obviously,
;                large values of oversample will take a long time to process
;
;   anchor - will use first and last selected points as anchor points for defining a curve, not used in
;            slit extraction (HMM... THIS LAST BIT MY NOT BE TRUE!).
;
;   intim - provide an image to use for picking the guide points, e.g., unsharp mask image
;
;OUTPUTS - curvepoints - spline data points
;          slits - cube of time distance diagrams
;          gpoints - selected points from pickn
;
;
;:CALL PROCEDURE: spline_slit,data,pickn=6,slits=slits,gpoints=gpoints - if picking points from image, returns picked
;                   values to gpoints
;                spline_slit,data,gpoints=gpoints - if data points for spline-ing known
;
;:EXTERNAL CALLS: tvim.pro
;
;:HISTORY: Author:
;   Oliver Scott   06/2015  - Initial coding
;   R. Morton      11/2015  - Removed redundant code, added common variables,
;                             and removed/added some keywords
;   M. Weberg      03/2018  - Upgraded the pick_n function to keep picking until
;                             the user decides to stop.
;                           - Removed the call to partvelvec.pro
;                           - Also cleaned up the formatting of the code and
;                             renamed internal variables to make it easier to read
;                             and understand
;
;:TO DO:
;   - [DONE!] Use mouse to bring point picking to end rather than entering
;     number at command line
;   - Upgrade the code to offer more output information and play nice with the
;     automated version of NUWT
;
;


;########################################################################################
;;;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;;;
;; The procedure takes an image, displays it,
;; and allows the user to select as many points
;; on the image as they would like (using the cursor).
;; The coordinates of these points are output.
;;;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;;;

FUNCTION PICK_GUIDE_POINTS, n

COMMON SOLVAR, imagesum, times

sz_img = SIZE(imagesum)
nx = sz_img[1]
ny = sz_img[2]

; Plot the summed image
WINDOW, /free
tvim, imagesum, title='Sum of the images from time ' + STRTRIM(times[0], 2) $
                       + ' to time ' + STRTRIM(times[1], 2)

temp_xpoints = LIST()
temp_ypoints = LIST()

PRINT, ''
PRINT, 'Please pick the guide points for the spline by left-clicking on the image'
PRINT, '---Right click at anytime to stop picking points---'

;Let the user pick as many guide points as they would like
num_gpoints = 0
loop_and_pick = 1
WHILE loop_and_pick EQ 1 DO BEGIN
    ;Removed the needlessly opaque call to "PICK.PRO"
    ;(more self-contained and only added 1 extra line of code!)
    CURSOR, xpoint, ypoint, /Data, /Down
    IF !MOUSE.BUTTON EQ 1 THEN BEGIN
        IF (xpoint GE 0) AND (xpoint LE nx-1) $
        AND (ypoint GE 0) AND (ypoint LE nx-1) THEN BEGIN
            ;Only save points that lie inside the actual image
            temp_xpoints.ADD, xpoint
            temp_ypoints.ADD, ypoint
            num_gpoints = num_gpoints + 1
            PLOTS, [xpoint,ypoint], psym=1
        ENDIF ELSE BEGIN
            PRINT, 'ERROR: Please select a point inside the image!'
        ENDELSE
    ENDIF ELSE IF !MOUSE.BUTTON EQ 4 THEN BEGIN
        IF num_gpoints GE 2 THEN BEGIN
            loop_and_pick = 0
        ENDIF ELSE BEGIN
            ;Don't stop if there are not enough points selected
            PRINT, 'Please select more guide points! At least two points are needed'
        ENDELSE
    ENDIF
ENDWHILE

PRINT, ''
PRINT, STRTRIM(num_gpoints, 2)+' guide points selected'
;Transfer guide points to the output array
gpoints = MAKE_ARRAY(num_gpoints,2)
FOR g=0,(num_gpoints-1) DO gpoints[g,*] = [temp_xpoints[g], temp_ypoints[g]]

RETURN, gpoints
END



FUNCTION EXTRACT_SPLINE_SLITS, timeslice, normals, slit_coords=slit_coords, $
                               d=d, plot_slit=plot_slit, missing=missing
;+
;:DESCRIPTION:
;   Takes an input timeslice of image data and extracts (using interpolation)
;   intensity values along slits described by a center point and a normal vector
;   (as given by the output of the procedure SOLCURVEFIT procedure.
;   This function is a revised version of the procedure IMAGE_INTERPOLATE originally
;   coded by Oliver Scott in 2015. Redundant outputs have been removed and the array
;   shapes have been modified for more efficient memory access.
;
;:POSITIONAL INPUTS:
;   timeslice - single 2D image from the input 3D datacube
;   normals - array of slit centers, gradients, and normal unit vector components
;
;:KEYWORD INPUTS:
;   d - half-length of a single slit. Default is 10
;   /plot_slit - If set, will overplot the slits on the example image
;   missing - Value to assign to missing data points. Default is -999 (in this function)
;             but the main function uses -1 if MIN(datacube) > 0) OR
;             MIN(datacube) - STDDEV(datacube)
;
;:RETURN VALUE:
;   interp_mat - Array of interpolated intensity values taken along each slit.
;                The dimensions of the array is [2d+1, S] where S is the number
;                of slits extracted
;
;:KEYWORD OUTPUTS:
;   slit_coords - x & y coordinates for all points along each slit.
;                 The dimensions of the array is [2d+1, 2, S] where S is the
;                 number of slits extracted. Row 0 contains the X values and
;                 row 1 contains the y values
;-

IF NOT KEYWORD_SET(d) THEN d = 10
IF N_ELEMENTS(missing) EQ 0 THEN missing = -999

;quick key to the normals array:
;   normals[s,0] & normals[s,1] - x & y coordinates of the slit center
;   normals[s,2]                - differential of the cerve at the point [x,y]
;   normals[s,3] & normals[s,4] - x & y components of the normal vector
;where s indicates which slit is onder consideration
sz_normals = SIZE(normals)
num_slits = sz_normals[1]
slit_coords = MAKE_ARRAY(2*d+1, 2, num_slits)
interp_mat = MAKE_ARRAY(2*d+1, num_slits)

FOR i=0,(num_slits-1) DO BEGIN
    slit_coords[*,0,i] = (normals[i,0]-d*normals[i,3])+2*d*normals[i,3]*FINDGEN(2*d+1)/(2*d)
    slit_coords[*,1,i] = (normals[i,1]-d*normals[i,4])+2*d*normals[i,4]*FINDGEN(2*d+1)/(2*d)
    interp_mat[*,i] = INTERPOLATE(timeslice, slit_coords[*,0,i], slit_coords[*,1,i], cubic=-0.5, missing=missing)
    IF KEYWORD_SET(plot_slit) THEN OPLOT, slit_coords[*,0,i], slit_coords[*,1,i]
ENDFOR

RETURN, interp_mat
END

;########################################################################################
;;;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;;;
;;
;; The procedure performs parametric spline fitting to a data set that
;; is input either through inputting an array of points, or by using a
;; procedure to allow the user to place points onto an image directly.
;;
;; Once a curve is fit to the data, the procedure calculates the
;; differential of the curve (dy/dx) using centred finite differences
;; (forward finite differences and backwards finite differences used
;; for the first and last points respectively). The differential at
;; each point is used to calculate the gradient of the normal vector
;; at each point in the fitted curve. Subsequent unit normal vectors
;; are then calculated for each point.
;;
;;; Procedure Arguments
;;
;;    spline_arr - A two dimensional array containing x,y values of
;;                 the fitted curve
;;
;;    normal_arr  - A two dimensional array containing x,y values of
;;                 points, the differential to the curve at these
;;                 points, and the unit normal vector to the
;;                 curve at these points
;;
;;                 These values are stored in rows, thus the first two
;;                 rows give the x,y values respectively; the next row
;;                 gives the differential of the curve at the
;;                 corresponding point;and the next two rows give the
;;                 x and y components of the unit normal vector at the
;;                 corresponding point.
;;
;;;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;;;

PRO SOLCURVEFIT, spline_arr, normal_arr, diffvec, $
                 spacing=spacing, vector=vector, pickn=pickn, $
                 gpoints=gpoints, oversample=oversample, half_length=half_length


COMMON SOLVAR, imagesum, times

IF KEYWORD_SET(pickn) THEN BEGIN
    gpoints = pick_guide_points(pickn)
    ; data vectors extracted from data array
    xvec = gpoints[*,0]
    yvec = gpoints[*,1]

ENDIF ELSE BEGIN

    ;; Check for valid array size and orientation
    ; Do note the edge case of when N=2 : it is possible for the data to
    ; be inputted incorrectly but there is no method to check this.
    sz_gpoints = SIZE(gpoints)

    IF sz_gpoints[0] NE 2 THEN BEGIN
        ;1D and 3D (or larger) input arrays
        MESSAGE, 'gpoints array has the wrong dimensions. Please input a 2D array of dimensions [N,2].'
    ENDIF ELSE IF (sz_gpoints[1] EQ 2) AND (sz_gpoints[2] NE 2) THEN BEGIN
        ;2D array with rows and columns reversed
        gpoints = TRANSPOSE(gpoints)
    ENDIF ELSE IF sz_gpoints[2] NE 2 THEN BEGIN
        ;2D array with the wrong shape
        MESSAGE, 'gpoints array has the wrong shape. Please input a 2D array of dimensions [N,2].'
    ENDIF

    ; Data vectors extracted from data array
    xvec = gpoints[*,0]
    yvec = gpoints[*,1]

    ;Question: What is the point of running the slit extraction program without
    ;          any input data!?! Debugging?

    ; ;; Finds data plot range for x and y if using supplied data array
    ; ;; without an image supplied
    ; ; ranges found by extending the x and y ranges by 50% centred around
    ; ; central value for x and y respectively.
    ;
    ; ; for x
    ; xmin = MIN(xvec,max=xmax)
    ; oxrange = [xmin,xmax] ; original x range stored
    ; xrange = xmax-xmin
    ; xmid = 0.5*(xmin+xmax)
    ; xmin = xmid-0.75*xrange
    ; xmax = xmid+0.75*xrange
    ;
    ; ; for y
    ; ymin = MIN(yvec,max=ymax)
    ; oyrange = [ymin,ymax] ; original y range stored
    ; yrange = ymax-ymin
    ; ymid = 0.5*(ymin+ymax)
    ; ymin = ymid-0.75*yrange
    ; ymax = ymid+0.75*yrange

    tvim, imagesum
    OPLOT, xvec, yvec, psym=1

ENDELSE

;;;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;;;
;;; Parametric Splining

IF N_ELEMENTS(spacing) GT 0 THEN BEGIN
    IF N_ELEMENTS(oversample) GT 0 THEN BEGIN
        ;; Oversampling works by setting the spacing to be smaller to make the
        ;; point density larger. Then by summing the lengths between the
        ;; points together to make up the original desired length the correct
        ;; point spacing is achieved while reducing the impact of any error in
        ;; the splining routine.
        ;;
        ;; The check for when the inter-point distance (from oversampled
        ;; curve) has reached the desired spacing is achieved by finding the
        ;; cumulative distance along the (oversampled) curve. Then by dividing
        ;; this cumulative distance vector by the original sampling it is
        ;; apparent that every time the desired original spacing is found the
        ;; floor of the cumulative distance will increase by 1. Flooring has
        ;; been achieved by changing variable type to integer. Creating a
        ;; logical vector that is true when the respective element in the
        ;; floored cumulative distance vector increases by 1. This logical
        ;; vector is used to sample points in the oversampled curve, achieving
        ;; the original interpoint spacing as desired, with reduced errors.

        ospacing = spacing
        spacing = DOUBLE(spacing)/(10^oversample)

        SPLINE_P, xvec, yvec, xsplines, ysplines, interval=spacing

        point_num = N_ELEMENTS(xsplines)
        spline_arr = [[xsplines],[ysplines]]

        ;linear distance between spline points
        diff_vec = SQRT( (SHIFT(spline_arr[*,0],-1)-spline_arr[*,0])^2 + $
                         (SHIFT(spline_arr[*,1],-1)-spline_arr[*,1])^2)
        diff_vec = [0, diff_vec[0:-2]]

        cum_diff_vec = FLTARR(point_num)
        FOR i=0,(point_num-1) DO cum_diff_vec[i] = TOTAL(diff_vec[0:i])

        cum_diff_vec = cum_diff_vec/ospacing
        cum_diff_vec = ULONG64(cum_diff_vec) ;i.e. remove the fraction bits

        logical_vec = cum_diff_vec
        FOR i=1,(point_num-1) DO logical_vec[i] = (cum_diff_vec[i] NE cum_diff_vec[i-1])
        logical_vec[0] = 1

        xsplines = xsplines[WHERE(logical_vec EQ 1)]
        ysplines = ysplines[WHERE(logical_vec EQ 1)]

    ENDIF ELSE BEGIN
        ;User-defined slit spacing
        ; use spline_p routine to generate spline interpolation points
        SPLINE_P, xvec, yvec, xsplines, ysplines, interval=spacing
    ENDELSE
ENDIF ELSE BEGIN
    ;Slit spacing automatically determined by the spline_p procudure
    ;Note: returns ~8 interpolated points for each segment (sometimes less)
    SPLINE_P, xvec, yvec, xsplines, ysplines
ENDELSE

;Plot the interpolated spline curve
OPLOT, xsplines, ysplines, psym=0

spline_arr = [[xsplines],[ysplines]]

;;;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~;;;
;;; Calculation of Unit Normal Vectors
;; find finite differences

num_slits = N_ELEMENTS(xsplines)

; initalise vector to store differences
tan_arr = MAKE_ARRAY(num_slits, 2, /float) ; array of tangent vectors between points
norm_arr = MAKE_ARRAY(num_slits, 2, /float) ; array of normals at each point

; use forward and backward differences on first and last elements where centred difference is ill defined
tan_arr[0,*] = [[xsplines[1]-xsplines[0]],[ysplines[1]-ysplines[0]]]
tan_arr[-1,*] = [[xsplines[-1]-xsplines[-2]],[ysplines[-1]-ysplines[-2]]]

; array manipulation for centred finite differences
; [2018-03-08]: corrected the equation by adding the missing factors of 1/2 (thankfully the error canceled out in the old code)
tan_arr[1:-2,*] = [[xsplines[2:-1]-xsplines[0:-3]]/2.0,[ysplines[2:-1]-ysplines[0:-3]]/2.0]


;; create normals by finding gradients and creating vectors
; initalisation of array of normal vectors

diff_arr = tan_arr[*,1]/tan_arr[*,0]

; this loop creates the normals
FOR i=0,(num_slits-1) DO BEGIN
    IF (diff_arr[i] eq 0) THEN BEGIN
        norm_arr[i,*] = [[0],[1]]
    ENDIF ELSE BEGIN
        norm_arr[i,*] = (1/SQRT(1+(-1/diff_arr[i])^2))*[[1],[-1/diff_arr[i]]]
    ENDELSE

    ; Test and use the cross product to ensure correct orientation of normals
    IF tan_arr[i,0]*norm_arr[i,1] LT tan_arr[i,1]*norm_arr[i,0] THEN BEGIN
        norm_arr[i,*] = -norm_arr[i,*]
    ENDIF
ENDFOR

;One of the key outputs!
normal_arr = [[spline_arr],[diff_arr],[norm_arr]]

;Plot the vectors
;[2018-03-08]: Removed the call to partvelvec.pro since it is non-standard IDL
;              and had too many other dependencies (namely, the coyote graphics programs)
IF KEYWORD_SET(vector) THEN BEGIN
    d = half_length
    vec_color_str = 'lime'
    IF TAG_EXIST(!COLOR, vec_color_str, index=color_index) THEN BEGIN
        vc_color_rgb = !COLOR.(color_index)
    ENDIF ELSE BEGIN
        vc_color_rgb = !COLOR.WHITE ;default vector color
    ENDELSE
    TVLCT, current_colortable, /GET ;save a copy of the current colortable
    TVLCT, REFORM(vc_color_rgb, 1, 3), 255 ;modify the color table to have the desired color at the end
    ARROW, spline_arr[*,0], spline_arr[*,1], $
           spline_arr[*,0]+norm_arr[*,0]*d, spline_arr[*,1]+norm_arr[*,1]*d, $
           /data, hsize=-0.25, color=255
    TVLCT, current_colortable ;reset the color table to whatever it was before
ENDIF

;diffvec is the linear differnce between the centre of each slit taken along the spline
;(this does not seem to be used anywhere, perhaps it is just for debugging...)
diffvec = SQRT(((spline_arr[1:-1,0]-spline_arr[0:-2,0])^2)+((spline_arr[1:-1,1]-spline_arr[0:-2,1])^2))

END





;########################################################################################
;MAIN ROUTINE

PRO NUWT_SPLINE_SLIT, in_data, slits=slits, frames=frames,$
                      pickn=pickn, distance=distance, anchor=anchor,$
                      gpoints=gpoints, spacing=spacing, $
                      oversample=oversample, plot_slit=plot_slit, $
                      vector=vector, curvepoints=curvepoints, intim=intim, $
                      slit_coords=slit_coords, outvec=outvec

COMMON SOLVAR, imagesum, times
on_error,2

;Set Default values
IF N_ELEMENTS(distance) EQ 0 THEN distance = 10
IF N_ELEMENTS(gpoints) EQ 0 AND NOT KEYWORD_SET(pickn) THEN pickn = 1
IF KEYWORD_SET(pickn) AND n_elements(gpoints) gt 0 THEN BEGIN
    PRINT, 'WARNING: the input array of gpoints will be overwritten with the points selected using the cursor.'
ENDIF
IF N_ELEMENTS(spacing) GT 1 THEN BEGIN
    PRINT, 'WARNING: multiple spacing values given! Only the first value will be used.'
    spacing = spacing[0]
ENDIF

;@check_dims####################################################################
;CHECKING DIMENSIONS OF INPUT DATA
;###############################################################################
sz_data = size(in_data)
IF sz_data[0] EQ 3 THEN BEGIN
    nx = sz_data[1]
    ny = sz_data[2]
    nt = sz_data[3]
ENDIF ELSE BEGIN
    IF sz_data[0] EQ 2 THEN BEGIN
        print, 'WARNING: input data is only 2D! Each output slit will be 1D.'
        nx = sz_data[1]
        ny = sz_data[2]
        nt = 1
        intim = in_data[*,*] ;will be passed to imagesum (prevents a crash)
    ENDIF ELSE BEGIN
        MESSAGE, 'ERROR: please input a valid 2D or 3D array of images'
    ENDELSE
ENDELSE

IF N_ELEMENTS(frames) NE 2 THEN times = [0,nt-1] ELSE times = [frames[0],frames[1]]

IF N_ELEMENTS(intim) GT 0 THEN imagesum = intim ELSE imagesum = total(in_data[*,*,times[0]:times[1]],3)


; uses the solcurvefit procedure to fit a curve to the feature
solcurvefit, spline_arr, normal_arr, diffvec, $
             spacing=spacing, vector=vector, pickn=pickn, $
             oversample=oversample, gpoints=gpoints, $
             half_length=distance

curvepoints = spline_arr

IF keyword_set(anchor) THEN BEGIN
    xvec = gpoints[*,0]
    yvec = gpoints[*,1]

    ;Define circle that passes through loop footpoints
    vecsize = n_elements(xvec)
    cutoff = [xvec[1],yvec[1]]
    cutoff2 = [xvec[vecsize-2],yvec[vecsize-2]]

    r = sqrt((cutoff2[0]-cutoff[0])^2+(cutoff2[1]-cutoff[1])^2)
    mid = [cutoff2[0]+cutoff[0],cutoff2[1]+cutoff[1]]/2 ;midpoint between footpoints


    cx = mid[0]-sqrt(r^2-(r/2)^2)*(cutoff[1]-cutoff2[1])/r
    cy = mid[1]-sqrt(r^2-(r/2)^2)*(cutoff2[0]-cutoff[0])/r


    ;define polar array
    pol = fltarr(nx,ny,2)
    x = findgen(nx)
    y = findgen(ny)
    x = rebin(x,nx,ny)
    y = transpose(rebin(y,ny,nx))

    pol[*,*,0] = sqrt((x-cx)^2+(y-cy)^2)
    pol[*,*,1] = atan((y-cy)/(x-cx))

    sp_pol = sqrt((curvepoints[*,0]-cx)^2+(curvepoints[*,1]-cy)^2)
    in = where(sp_pol GE r, complement=in2)
    ;   IF n_elements(in) GT in2 THEN

    curvepoints = curvepoints[in,*]
    normal_arr = normal_arr[in,*]

ENDIF

;Count how many output slits there are
sz_curvepoints = SIZE(curvepoints)
num_slits = sz_curvepoints[1]

mindat = MIN(in_data)
IF mindat GT 0 THEN mindat = -1.0 ELSE mindat = mindat - 0.1*STDDEV(in_data)

; initialize the array for the extracted slits
interp_arr = MAKE_ARRAY(2*distance+1, nt, num_slits)

FOR i=0,nt-1 DO BEGIN
   interp_arr[*,i,*] = extract_spline_slits(in_data[*,*,i], normal_arr, slit_coords=s_coords, d=distance, plot_slit=plot_slit)
ENDFOR

outvec = fltarr(4, num_slits)
outvec[0,*] = s_coords[0,0,*] ;start x vals for each slit
outvec[1,*] = s_coords[0,1,*] ;start y vals
outvec[2,*] = s_coords[-1,0,*] ;end x vals
outvec[3,*] = s_coords[-1,1,*] ;end y vals

;Finally, output the extracted slits
slits = interp_arr
slit_coords = s_coords
END
