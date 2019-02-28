;+
;NAME: ANDERSON_DARLING
;
;PURPOSE:
;   Determine if a given set of data is normally distributed by application
;   of the Anderson Darling test. This test is generally seen as more powerful
;   than the
;
;INPUTS:
;   x_vals - sample data to be tested
;
;OPTIONAL INPUTS:
;   adjusted - use the adjusted AD test statistic as given by D'Agostino, 1986
;   signif_level - significance level to use for determining the critical value
;   alt_2d_interp - if set, will use two 1D interpolations instead of one 2D
;                   interpolation. MUST be set when running on a mac.
;
;OUTPUTS:
;   A_squared - the Anderson Darling test statistic. Should be compared to
;               a table of critical values to determine quality of fit.
;   crit_val - critial value at the 5% significance level
;
;HISTORY: Name---------Date---------Description
;         M Weberg  July, 2016 - Initial coding
;         M Weberg  Sept, 2016 - Added tables of critical values with interpolation
;                                to the selected signif_level with the given
;                                number of data points
;         M Weberg   Nov, 2016 - Added alt 2D interpolation option since there is
;                                an IDL bug in TRIGRID on mac computers.
;
;TO DO:
;   - Fix the 2D TRIGRID interpolation on mac computers!
;-

FUNCTION ANDERSON_DARLING, x_vals, adjusted=adjusted, signif_level=signif_level, $
                                   crit_val=crit_val, alt_2d_interp=alt_2d_interp

;Set default significance level
IF n_elements(signif_level) EQ 0 THEN signif_level = 0.95

N = n_elements(x_vals)

mu = total(x_vals)/N
sigma = sqrt((1.0/(N-1))*total((x_vals-mu)^2))

sort_ind = SORT(x_vals)
Y = (x_vals[sort_ind]-mu)/sigma

norm_CDF_vals = standard_normal_cdf(Y)

i_arr = indgen(N)+1

A_sqrd = -N -1.0/N*total((2*i_arr - 1)*(ALOG(norm_CDF_vals) + ALOG(1 - REVERSE(norm_CDF_vals))))

; Find the critical value for the given significance level
IF KEYWORD_SET(adjusted) AND N GE 8 THEN BEGIN
    ;from D'Agostino, 1986, (book)
    ;Crit vals table (N => 8)
    ; 10%,   5%,   2.5%,   1%,   0.5%
    ;0.631, 0.752, 0.873, 1.035, 1.159
    A_sqrd = A_sqrd*(1 + 0.75/N + 2.25/N^2)
    crit_val = interpol([0.631, 0.752, 0.873, 1.035, 1.159], [0.90, 0.95, 0.975, 0.99, 0.995], signif_level)
ENDIF ELSE BEGIN
    ;Crit val table (from Stephens, 1974, JASA - pg 733, Table 2, Case 3)
    ; n,    15%,    10%,    5%,     2.5%,   1%
    ; 10    0.514   0.578   0.683   0.779   0.926
    ; 20    0.528   0.591   0.704   0.815   0.969
    ; 50    0.546   0.616   0.735   0.861   1.021
    ; 100   0.559   0.631   0.754   0.884   1.047
    ; Inf   0.576   0.656   0.787   0.918   1.092

    ;Interpolation of critical values on [signif_vals, N_vals] grid
    N_coords = [10, 20, 50, 100, 1e37] ;note: 1e37 will serve as a proxy for infinity (better results than simple extrapolation)
    S_coords = [0.85, 0.90, 0.95, 0.975, 0.99]
    crits = fltarr(5, 5)
    ; sig_level,  0.85,  0.90,  0.95,  0.975, 0.99 -- N=
    crits[*,0] = [0.514, 0.578, 0.683, 0.779, 0.926] ;10
    crits[*,1] = [0.528, 0.591, 0.704, 0.815, 0.969] ;20
    crits[*,2] = [0.546, 0.616, 0.735, 0.861, 1.021] ;50
    crits[*,3] = [0.559, 0.631, 0.754, 0.884, 1.047] ;100
    crits[*,4] = [0.576, 0.656, 0.787, 0.918, 1.092] ;1e37 (proxy for infinity)

    IF KEYWORD_SET(alt_2d_interp) THEN BEGIN
        ;Alternative method which uses two 1D interpolations as a proxy for one 2D
        ;This is the method that should be used on macs for the time being
        interp_crit_table = fltarr(5)
        FOR s_index=0, 4 DO BEGIN
            interp_crit_table[s_index] = interpol(crits[s_index, *], N_coords, N)
        ENDFOR
        crit_val = interpol(interp_crit_table, S_coords, signif_level)
    ENDIF ELSE BEGIN
        ;Trinagular interpolation
        ;Note: this seems to have issues on mac computers
        NN = rebin(reform(N_coords, 1, 5), 5, 5)
        SS = rebin(S_coords, 5, 5)
        TRIANGULATE, SS, NN, triangles, b_nodes
        output_grid = TRIGRID(SS, NN, crits, triangles, xout=[signif_level, 0.995], yout=[N, 1000], extrapolate=b_nodes)
        crit_val = output_grid[0,0]
    ENDELSE
ENDELSE

RETURN, A_sqrd
END
