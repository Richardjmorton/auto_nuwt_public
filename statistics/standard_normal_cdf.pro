;+
;NAME: STANDARD_NORMAL_CDF
;
;PURPOSE:
;   Compute the cumulative distribution function (CDF) for the standard normal
;   distribution (i.e. a gaussian with mean = 0 and stddev = 1).
;
;INPUTS:
;   x_vals - x locations to sample the standard normal CDF at
;
;OPTIONAL INPUTS:
;   NONE (at the moment)
;
;OUTPUTS:
;   norm_CDF - standard normal CDF values at the given x locations
;
;HISTORY: Name---------Date---------Description
;         M Weberg  July, 2016 - Initial coding
;-

FUNCTION STANDARD_NORMAL_CDF, x_vals

norm_CDF = 0.5*(1 + ERF(x_vals/SQRT(2)))
RETURN, norm_CDF

END
