;+
;NAME: TEST_RESIDUALS
;
;PURPOSE:
;   Quantify the goodness-of-fit of a given data model by appliying two statistical
;   tests on the residuals. Tests for both a normal distribution, using the
;   Kolmogorov–Smirnov (KS) test, and autocovarience, using the Ljung–Box test.
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
;         M Weberg   Nov, 2016 - Cleaned up code a little and added option for AD test
;-

FUNCTION TEST_RESIDUALS, residuals, signif_level=signif_level, num_fit_params=num_fit_params, $
                         max_lag=max_lag, debug_hist=debug_hist

;Set default significance level
IF n_elements(signif_level) EQ 0 THEN signif_level = 0.95

;Performing the Kolmogorov–Smirnov one sample test vs a gaussian CDF to determine normality
KS_stat = kolmogorov_smirnov(residuals, prob=ks_prob)

;Performing the Anderson Darling test for normality
;note: the /alt_2d_interp option is needed on macs to avoid an IDL bug
IF !version.OS_NAME EQ 'Mac OS X' THEN $ 
A_sqrd = anderson_darling(residuals, signif_level=signif_level, crit_val=AD_crit_val, /alt_2d_interp) $
ELSE A_sqrd = anderson_darling(residuals, signif_level=signif_level, crit_val=AD_crit_val);, /alt_2d_interp)

;Performing the Ljung-Box test to detect autocorrelation
Q = ljung_box(residuals, signif_level=signif_level, max_lag=max_lag, $
              num_fit_params=num_fit_params, crit_val=chisqr_cutoff)

;Initialize output array
test_stats = {signif_level:signif_level, $
              KS_stat:KS_stat, $
              KS_prob:KS_prob, $
              AD_stat:A_sqrd, $
              AD_crit:AD_crit_val, $
              LB_stat:Q, $
              LB_chisqrd:chisqr_cutoff}

; debug_hist = {residuals:residuals, $
;               norm_resid:standard_residuals, $
;               hist_vals:resid_hist, $
;               hist_bins:resid_bins}

RETURN, test_stats
END
