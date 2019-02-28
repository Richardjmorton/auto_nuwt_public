;+
;NAME: KOLMOGOROV_SMIRNOV
;
;PURPOSE:
;   Calculate the one-sided Kolmogorov-Smirnov (ks) test statistic of given data
;   relative to a normal (gaussian) distribution. Adapted from 'ksone.pro' and
;   'prob_ks.pro' from the IDL astro library. This version is less robust however,
;   it also operates more quietly (no extra 'COMPILED MODULE: XYZ' statements
;   with every single use).
;
;INPUTS:
;   x_vals - sample data to be tested
;
;OUTPUTS:
;   D - the KS test statistic which corresponds to the maximum difference
;       between the data cumulative distribution functions (CDF) and
;       the standard normal CDF.
;   prob - significance level (between 0 and 1) of the results. A very LOW
;          value indicates the data is NOT normally distributed
;
;HISTORY: Name---------Date---------Description
;         M Weberg  Sept, 2016 - Initial adaptation
;
;TO DO:
;   - Include options for other test distributions
;-

FUNCTION KOLMOGOROV_SMIRNOV, x_vals, prob=prob

N = n_elements(x_vals)
mu = total(x_vals)/N
sigma = sqrt((1.0/(N-1))*total((x_vals-mu)^2))

sort_ind = SORT(x_vals)

;Sorted (ascending) and standardized data values
Y = (x_vals[sort_ind]-mu)/sigma

;##### Calculate CDFs and the KS test statistic (ADAPTED FROM 'ksone.pro')
f0 = findgen(N)/ N
fn = (findgen(N) +1.0) / N
ff = standard_normal_cdf(Y)
D = max([max(abs(f0-ff), sub0), max(abs(fn-ff), subn)], msub)

;##### Calculate the significance level (ADAPTED FROM 'prob_ks.pro')
eps1 = 0.001    ;Stop if current term less than EPS1 times previous term
eps2 = 1.e-8    ;Stop if current term changes output by factor less than EPS2

en = sqrt(N)
lambda = (en + 0.12 + 0.11/en)*D

a2 = -2.0*lambda^2
probks = 0.0
termbf = 0.0
sign = 1.0

j = 1 ; loop variable
converged = 0 ; Flag for breaking out of the loop
WHILE j LE 100 AND converged EQ 0 DO BEGIN
    term = sign*2*exp(a2*j^2)
    probks = probks + term

    IF ( abs(term) LE eps1*termbf ) OR $
       ( abs(term) LE eps2*probks ) THEN converged = 1

    sign = -sign                  ;Series alternates in sign
    termbf = abs(term)
    j += 1
ENDWHILE

IF converged EQ 0 THEN probks = 1.0 ;Sum did not converge after 100 iterations

prob = probks

RETURN, D
END
