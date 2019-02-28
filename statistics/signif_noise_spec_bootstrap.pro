;+
;NAME: SIGNIF_NOISE_SPEC_BOOTSTRAP
;
;PURPOSE: 
;   Calculate and return the significant FFT power spectrum for either white or
;   red (aka brown) noise at a given significance level. Assumes two degrees of
;   freedom (DOF) but may be modified as desired. This program is based on the
;   'signif_conf.pro' code and implements the equations found in Torrence and 
;   Compo, 1998 (BAMS). 
;
;INPUTS:
;   data_in - FFT power spectrum of the data
;   p - desired significance level. Example, 0.95 corresponds to the 95% 
;       significance level which means there is only a 5% chance that a random
;       process will return a value larger than the test value. This is also 
;       sometime called the "5% level of significance" or the "95% confidence level". 
;       (latter is supposedly something of a misnomer).
;   num_FFT - number of FFT bins for required in the output significance spectrum.
;             WARNING! currently will not work properly for zero padded FFTs
;
;OPTIONAL INPUTS:
;   color - noise model to use for generating the significance values. choose from 
;           'white' or 'red'. 'White noise' assumes a Gaussian distribution with a
;           flat FFT power spectrum while 'red noise' assumes more power in the
;           low frequency bins than the high frequency bins. Default is 'white'. 
;
;OUTPUTS:
;   signif_vals - array of same size as the input data containing the significance
;                 value in each frequency bin
;   AR_coeff - average autoregression coeffient calculated from (a_1 + sqrt(a_2))/2.
;              White noise has an AR coeffient of 0.0
;
; 
;HISTORY: Name---------Date---------Description
;         M Weberg  July, 2016 - Initial "coding" (modification of'signif_conf.pro')
;
;ORIGINAL HEADER TEXT:
;   commented using http://cires.colorado.edu/geo_data_anal/assign/assign4.html
;   Assuming a Gaussian distribution for each point of your NINO3 timeseries, then each
;   point of the FFT Power Spectrum should have a chi-square
;   distribution with 2 degrees of freedom (DOF) (see Figure).
;   Looking at the curve, one can see that as one goes further to
;   the right (larger x), the probability of being greater than x
;   decreases. One can define the 95% level as that value of x
;   where there is only a 5% chance of being greater than x. This is referred to as either the
;   "95% significance level" or the "5% level of significance". Occasionally it is referred to as
;   the "95% confidence level", but the term "confidence" should really be reserved for
;   reference to "confidence intervals", which are the error bars seen on spectral plots. 
;
;   For the FFT Power Spectrum, this means that if one were to choose 20 frequencies at
;   random, then only 1 of these frequencies would be expected to have FFT Power greater
;   than this value. Therefore, if you look at your power spectrum, and several peaks are
;   above the 95% level, then you can be reasonably "confident" that these are "real" peaks,
;   not just random noise. 
;-

FUNCTION SIGNIF_NOISE_SPEC_BOOTSTRAP, data_in, p, num_FFT, color=color, AR_coeff=AR_coeff

;Setting default values
IF NOT keyword_set(color) THEN color = 'white'

size_data = size(data_in)
N = size_data[1] ;number of data points in each sample
N_boot = size_data[2] ;number of bootstrapping resamples

s = rebin(reform(stddev(data_in, dimension=1), 1, N_boot), num_fft, N_boot)

DOF = 2

chi = chisqr_CVF(1.0-p, DOF)

;   The value of chi-square with 2 DOF at the 5% level of significance (95% significance) is
;   5.99. To convert this into a Power that can be plotted on your FFT Power Spectrum,
;   one needs to divide by the DOF (=2), multiply by the timeseries variance s^2, and divide
;   by the number of points in the spectrum N/2. 

;   In general, the formula for significance levels is:

;                  s^2 chi_sqr(1-p,DOF)
;        signif = --------------------  x P(k)
;                      0.5 N DOF

;   where s^2 is the variance, p is the desired significance level (such as 0.95), DOF is the
;   degrees of freedom (usually 2), and N is the number of points in the timeseries. 

IF color EQ 'white' THEN BEGIN
    P_k = 1.0 + fltarr(num_FFT, N_boot)
    alpha = 0.0
ENDIF

IF color EQ 'red' THEN BEGIN
    mu = rebin(reform(mean(data_in, dimension=1), 1, N_boot), N, N_boot)
    sigma = stddev(data_in, dimension=1)
    data_i_plus_1 = SHIFT(data_in, -1, 0)
    data_i_plus_2 = SHIFT(data_in, -2, 0)
    ;calculate the lag-1 (a1) and lag-2 (a2) coeffients
    a1 = 1.0/((N - 1)*sigma^2)*total(((data_in - mu)*(data_i_plus_1 - mu))[0:-2,*], 1)
    a2 = 1.0/((N - 2)*sigma^2)*total(((data_in - mu)*(data_i_plus_2 - mu))[0:-3,*], 1)
    a1 = rebin(reform(a1, 1, N_boot), num_FFT, N_boot)
    a2 = rebin(reform(a2, 1, N_boot), num_FFT, N_boot)
    alpha = (a1 + sqrt(abs(a2)))/2
    ;print, 'a1, a2, alpha =', a1, a2, alpha
    P_k = (1 - alpha^2)/(1 + alpha^2 - 2*alpha*cos(2*!PI*rebin(indgen(num_fft), num_fft, N_boot)/N))
ENDIF

signif = ((s^2)*(chi))/((N/2)*DOF)*P_k

signif = mean(signif, dimension=2)

AR_coeff = alpha

RETURN, signif
END


