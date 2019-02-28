;+
;NAME: CALC_KDE
;
;PURPOSE:
;   Calcuate the Kernal Density Estimate (KDE) of the input univariate data.
;   Currently limited to using a gaussian kernal
;
;INPUTS:
;   input_data - data values to compute the KDE of
;
;OPTIONAL INPUTS:
;   [to-be-written]
;
;OUTPUTS:
;   [to-be-written]
;
;HISTORY: Name---------Date---------Description
;         M Weberg   Nov, 2017 - Initial coding
;-

;@1d_kde################
; Univariate (1D) KDE
;#######################
FUNCTION CALC_UNIVARIATE_KDE, data_in, x_vals, bandwidth=bandwidth, h_out=h_out, print_h=print_h
    COMPILE_OPT IDL2, HIDDEN

    n_data = N_ELEMENTS(data_in)
    n_x = N_ELEMENTS(x_vals)

    data_mad = MEANABSDEV(data_in, /MEDIAN)
    data_mean = MEAN(data_in, /NAN)
    data_stddev = STDDEV(data_in, /NAN)

    sorted = data_in[SORT(data_in)]
    ;Type 8 quartiles in R (recommended by Hyndman & Fan, 1996, American Statistician, 50, 361-365)
    q1_ind = 0.25*(n_data+1.0/3.0)+1.0/3.0 - 1
    q2_ind = 0.5*(n_data+1.0/3.0)+1.0/3.0 - 1
    q3_ind = 0.75*(n_data+1.0/3.0)+1.0/3.0 - 1
    q1_floor_ind = FLOOR(q1_ind)
    q1_val = sorted[q1_floor_ind] + (q1_ind-q1_floor_ind)*(sorted[q1_floor_ind+1]-sorted[q1_floor_ind])
    q2_floor_ind = FLOOR(q2_ind)
    data_median = sorted[q2_floor_ind] + (q2_ind-q2_floor_ind)*(sorted[q2_floor_ind+1]-sorted[q2_floor_ind])
    q3_floor_ind = FLOOR(q3_ind)
    q3_val = sorted[q3_floor_ind] + (q3_ind-q3_floor_ind)*(sorted[q3_floor_ind+1]-sorted[q3_floor_ind])
    IQR = q3_val - q1_val

    IF NOT KEYWORD_SET(bandwidth) THEN BEGIN
        ;The normal distribution approximation (aka "Silverman's rule of thumb")
        ;Note: this is a very rough estimate which works best when the data is
        ;      normally distributed
        ; h = data_stddev*(4.0/(3.0*N))^(0.2) ;"Classical"
        h = 0.9*MIN([data_stddev, IQR/1.34])*n_data^(-0.2) ;"Reduced"
        IF KEYWORD_SET(print_h) THEN print, 'h =', h
    ENDIF ELSE BEGIN
        h = bandwidth[0]
        n_h = N_ELEMENTS(h)
        IF N_ELEMENTs(h) GT 1 THEN h = REBIN(h, n_data, n_X)
    ENDELSE

    ; ;Simple loop method (just to get it up and running)
    ; FOR i=0, (N-1) DO BEGIN
    ;     K_xi = EXP(-((x_vals-input_data[i])^(2.0))/(2.0*h^(2.0)))
    ;     sum_K = sum_K + k_xi
    ; ENDFOR

    ;Vectorized method
    data_matrix = REBIN(data_in, n_data, n_X)
    X_matrix = REBIN(REFORM(x_vals, 1, n_X), n_data, n_X)

    K_matrix = (1.0/h)*EXP(-((X_matrix-data_matrix)^(2.0))/(2.0*h^(2.0)))
    sum_K = TOTAL(K_matrix, 1)

    KDE = (1.0/(n_data*SQRT(2.0*!PI)))*sum_K
    h_out = h

    RETURN, KDE
END

;@main##########################################################################
;MAIN FUNCTION #################################################################
;###############################################################################
FUNCTION CALC_KDE, input_data, test_x=test_x, range=range, num_vals=num_vals, xvals_out=xvals_out, $
                   bandwidth=bandwidth, adaptive=adaptive, print_h=print_h, h_out=h_out, debug=debug

COMPILE_OPT IDL2
IF NOT KEYWORD_SET(debug) THEN BEGIN
    ;hides math underflow errors
    CurrentExceptVal = !Except
    !Except = 0
    void = CHECK_MATH()
ENDIF

;Calc basic properties of the input data
n_data = N_ELEMENTS(input_data)
data_median = MEDIAN(input_data)
data_mad = MEANABSDEV(input_data, /MEDIAN)
data_mean = MEAN(input_data, /NAN)
data_stddev = STDDEV(input_data, /NAN)

IF NOT KEYWORD_SET(num_vals) THEN num_vals = 2000
IF NOT KEYWORD_SET(range) THEN range = [data_median-3*data_mad, data_median+3*data_mad]
IF NOT KEYWORD_SET(bandwidth) THEN bandwidth = 0
IF NOT KEYWORD_SET(print_h) THEN print_h = 0

IF KEYWORD_SET(test_x) THEN BEGIN
    x_vals = test_x[*]
ENDIF ELSE BEGIN
    x_vals = FINDGEN(num_vals)*range[1]/(FLOAT(num_vals)) + range[0]
ENDELSE

IF KEYWORD_SET(adaptive) THEN BEGIN
    ;Adaptive KDE (i.e. the bandwidth, h, depends on location)
    KDE = calc_univariate_kde(input_data, input_data, bandwidth=bandwidth, h_out=h_fixed)
    ; h_i = ((4.0/(3.0*n_data))^(0.2))*SQRT(data_stddev/KDE)
    ; h_i = ((0.9/n_data)^(0.2))*SQRT(data_stddev/KDE) ;"reduced?"

    geometric_mean_KDE = EXP(MEAN(ALOG(KDE)))
    h_i = h_fixed*(geometric_mean_KDE/KDE)^(0.5) ;stolen from another code...
    ;note: all of the adpative methods seem to be oversmoothed. The basic one are
    ;      worse than the standard rule of thumb and the stolen version is
    ;      halfway between the standard and reduced rules of thumb.
    ;      At least for my test dataset, that is. Maybe it works better on more
    ;      complicated data....
    KDE = calc_univariate_kde(input_data, x_vals, bandwidth=h_i, h_out=h_val_out)
ENDIF ELSE BEGIN
    ;[DEFAULT] Simple KDE calculation
    KDE = calc_univariate_kde(input_data, x_vals, bandwidth=bandwidth, h_out=h_val_out)
ENDELSE

h_out = h_val_out
xvals_out = x_vals

IF NOT KEYWORD_SET(debug) THEN BEGIN
    ;manually check for other types of math errors
    floating_point_underflow = 32
    status = Check_Math() ; Get status and reset accumulated math error register.
    IF(status AND NOT floating_point_underflow) NE 0 THEN BEGIN
        print, 'IDL Check_Math() error: ' + StrTrim(status, 2)
    ENDIF

   !Except = CurrentExceptVal
ENDIF

RETURN, KDE
END
