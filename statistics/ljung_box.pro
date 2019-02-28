;+
;NAME: LJUNG_BOX
;
;PURPOSE:
;   Calculate the Ljung-Box statistic, Q, to test for autocorelation in the
;   residuals of a model fit. If Q > the chi^2 cutoff value, then the null
;   hypothesis that the data is randomly distributed with NO autocorrelation
;   is rejected.
;
;INPUTS:
;   residuals - values obtained by subtracting a model results from the data
;
;OPTIONAL INPUTS:
;   signif_level - significance level to use for determining the chi squared
;                  cut-off limit for the Ljung-Box test. Default is 0.95
;   fit_params - number of parameters fit by the model (or combination of models).
;                Used to apply a correction to the degrees of freedom in determining
;                the chi squared cut-off value. Default is 0 (no model)
;   max_lag - maximum autocorrelation lag to include in Ljung-Box test.
;             Default is num_residuals/2 (this is also the maximum value allowed)
;
;
;OUTPUTS:
;   test_stats - structure containing the test statistics as well as their
;                threshold values and probabilities
;
;
;HISTORY: Name---------Date---------Description
;         M Weberg  July, 2016 - Initial coding
;         M Weberg   OCT, 2016 - Fixed an issue with negative DOF for very short
;                                input data and many fit parameters
;-

FUNCTION LJUNG_BOX, residuals, crit_val=crit_val, signif_level=signif_level, $
                    num_fit_params=num_fit_params, max_lag=max_lag

N = n_elements(residuals)
half_n = floor(N/2.0)

;Set default significance level
IF n_elements(signif_level) EQ 0 THEN signif_level = 0.95

;Set default number of parameters that have been fit by the data model being tested
IF n_elements(num_fit_params) EQ 0 THEN num_fit_params = 0

;Set default number of lags to test in Ljung-Box test
;Note: will cap the maximum allowable lag at half the number of input data points
IF NOT KEYWORD_SET(max_lag) THEN BEGIN
    max_lag = half_n
ENDIF ELSE BEGIN
    IF max_lag GT half_n THEN max_lag = half_n
ENDELSE

;Performing the Ljung-Box test to detect autocorrelation
;note: we use the equation, AR(k) = sum_t=1_to_n-k[(x_t - mu)(X_{t+k} - mu)]/sum[(X - mu)^2]
lags = indgen(max_lag) + 1
AR = fltarr(max_lag)
mu = mean(residuals)
mean_diff = residuals - mu
FOR k = 0L, max_lag-1 DO BEGIN
    AR[k] = total(mean_diff[0:-lags[k]-1]*mean_diff[lags[k]:*])
ENDFOR

AR = AR/total(mean_diff^2)

Q = N*(N+2)*total(AR^2/(N-lags))

IF max_lag GE num_fit_params THEN BEGIN
    crit_val = chisqr_CVF(1.0-signif_level, max_lag-num_fit_params)
ENDIF ELSE BEGIN
    ;Prevents crashes and indicates odd results
    crit_val = 0.0
ENDELSE

RETURN, Q
END
