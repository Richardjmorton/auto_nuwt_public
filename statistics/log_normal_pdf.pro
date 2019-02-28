;+
;NAME: LOG_NORMAL_PDF
;
;PURPOSE:
;   Compute the probability distribution function (PDF) for a generalized
;   log-normal distribution (i.e. a gaussian with a given log_norm_mean and log_norm_stddev).
;
;INPUTS:
;   x_vals - x locations to sample the log-normal PDF at
;
;OPTIONAL INPUTS:
;   mu - mean of the distribution in log space
;   sigma - standard deviation of the distibution in log space
;   arith_mean - arithmetic mean of the distribution
;   arit_stddev - arithmetic stddev of the distribution
;
;OUTPUTS:
;   gen_norm_PDF - general normal PDF values at the given x locations
;
;HISTORY: Name---------Date---------Description
;         M Weberg  Mar, 2017 - Initial coding
;         M Weberg  Jun, 2017 - Corrected the equation (just scales peak height)
;-

FUNCTION LOG_NORMAL_PDF, x_vals, mu=mu, sigma=sigma, arith_mean=arith_mean, arith_stddev=arith_stddev
IF N_ELEMENTS(mu) EQ 0 THEN mu = 1.0
IF N_ELEMENTS(sigma) EQ 0 THEN sigma = 1.0

IF KEYWORD_SET(arith_mean) GT 0 AND N_ELEMENTS(arith_stddev) GT 0 THEN BEGIN
    mu = ALOG((arith_mean^2)/SQRT(arith_stddev^2 + arith_mean^2))
    sigma = SQRT(ALOG(1 + (arith_stddev^2)/(arith_mean^2)))
ENDIF

gen_log_norm_PDF = (1.0/(x_vals*sigma*SQRT(2*!PI)))*EXP(-((ALOG(x_vals) - mu)^2)/(2*sigma^2))
loc_nan = WHERE(FINITE(gen_log_norm_PDF, /NAN))

IF N_ELEMENTS(loc_nan) GT 0 THEN gen_log_norm_PDF[loc_nan] = 0.0

RETURN, gen_log_norm_PDF

END
