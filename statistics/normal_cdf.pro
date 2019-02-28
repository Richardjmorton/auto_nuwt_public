;+
;NAME: NORMAL_CDF
;
;PURPOSE:
;   Compute the cumulative distribution function (CDF) for a generalized normal
;   distribution (i.e. a gaussian with a given mean and stddev).
;   note: this is a modified version of "standard_normal_cdf.pro"
;
;INPUTS:
;   x_vals - x locations to sample the normal CDF at
;
;OPTIONAL INPUTS:
;   mu - mean of the normal distribution
;   sigma - standard deviation of the distibution
;
;OUTPUTS:
;   gen_norm_CDF - general normal CDF values at the given x locations
;
;HISTORY: Name---------Date---------Description
;         M Weberg  Mar, 2017 - Initial coding
;-

FUNCTION NORMAL_CDF, x_vals, mu=mu, sigma=sigma
IF NOT KEYWORD_SET(mu) THEN mu = 1.0
IF NOT KEYWORD_SET(sigma) THEN sigma = 1.0

gen_norm_CDF = 0.5*(1 + ERF((x_vals - mu)/(sigma*SQRT(2))))
RETURN, gen_norm_CDF

END
