;+
;NAME: CALC_AIA_ERRORS
;
;PURPOSE:
;   Calculate the intensity errors (i.e. noise) for the main SDO / AIA
;   channels following the method of Yuan & Nakariakov, 2012 (A&A) and using
;   the calibration values given by Boerner et al., 2012 (Solar Physics).
;
;INPUTS:
;   input_data - array of unaltered AIA intensity values
;   wavelength - central wavelength of the AIA channel. MUST select from 94, 131,
;                171, 193, 211, 304, or 335
;
;OPTIONAL INPUTS
;   norm_num - Number of timesteps used in smoothing the td diagram. Used to
;                apply a normalization factor. Default is 1 (assumes no smoothing)
;   /spike_noise - If set, will also include the error due to the despiking method.
;                  Normally this is a very small correction and can be ignored.
;   percent_spikes - Estimated percent of the pixels affected of energetic particles.
;                    Only used if /spike_noise is set. Default is 0.6%
;
;OUTPUTS:
;   output_data - array of intensity errors. Will have the same size and shape
;                 as "input_data"
;NOTES:
;   For the 171 channel, the default case is equivalent to,
;       output_data = SQRT(2.3 + 0.06*input_data) / SQRT(norm_num) ;(see Y&N 2012)
;
;HISTORY: Name---------Date---------Description
;         M Weberg  7 JUNE, 2016  Initial coding (for 171 only)
;         M Weberg 16  MAR, 2017  Expanded to include the other AIA wavelengths
;         M Weberg 27  MAR, 2018  Added the
;-

FUNCTION CALC_AIA_ERRORS, input_data, wavelength, norm_num=norm_num, $
                          spike_noise=spike_noise, percent_spikes=percent_spikes

; Setting defaults and channel indepenent variables
IF NOT keyword_set(percent_spikes) THEN percent_spikes = 0.6 ;[%]

;@check_data####################################################################
;CHECK THE INPUT DATA AND EXTRACT VARIABLES IF GIVEN A SLIT STRUCTURE
;###############################################################################
;Note: passing values to internal variables prevents the unintentional
;      modification of the inputs
IF N_TAGS(input_data) GT 0 THEN BEGIN
    ;Input is a structure containing both slit(s) and meta data
    input_tags = TAG_NAMES(input_data)
    IF TOTAL(STRCMP(input_tags, 'slit', /fold_case)) GE 1 THEN BEGIN
        DN = input_data.slit
        slit_meta = input_data.meta
        IF N_ELEMENTS(wavelength) EQ 0 THEN channel = FIX(slit_meta.wavelength) $
                                       ELSE channel = FIX(wavelength)
        IF N_ELEMENTS(norm_num) EQ 0 THEN norm_num = FIX(slit_meta.width)
    ENDIF ELSE BEGIN
        PRINT, 'ERROR: invalid structure format!'
        MESSAGE, 'Input data strutures must have both "slit" and "meta" tags'
    ENDELSE
ENDIF ELSE BEGIN
    ;Input is an array with one (or more) data slits
    DN = input_data
    channel = FIX(wavelength)
    IF NOT KEYWORD_SET(norm_num) THEN norm_num = 1
ENDELSE

;@lookup########################################################################
;LOOKUP THE CALIBRATED NOISE COMPONENTS FOR EACH DATA CHANNEL
;###############################################################################
;The values below are from Table 6 of Boerner et al., 2012 (Solar Physics)
IF channel EQ 131 OR channel EQ 335 THEN BEGIN
    ;AIA camera 1
    gain = 17.6
    read_noise = 1.18
ENDIF ELSE IF channel EQ 193 OR channel EQ 211 THEN BEGIN
    ;AIA camera 2
    gain = 18.3
    read_noise = 1.20
ENDIF ELSE IF channel EQ 171 THEN BEGIN
    ;AIA camera 3
    gain = 17.7
    read_noise = 1.15
ENDIF ELSE IF channel EQ 94 OR channel EQ 304 THEN BEGIN
    ;AIA camera 4
    gain = 18.3
    read_noise = 1.14
ENDIF ELSE BEGIN
    MESSAGE, 'ERROR: "wavelength" MUST be 94, 131, 171, 193, 211, 304, or 335'
ENDELSE

photon_factor = 1.0625/gain ;combines photon & compression noise which dep. on data number (DN)
const_factor = read_noise^2 + 1.0 ;readout, rounding, dark current, and subtraction noise
spike_factor = (percent_spikes/100.0)*0.15 ;noise from despiking method

;@validate######################################################################
;VALIDATE THE INPUT DATA VALUES
;###############################################################################
;Test for values that would cause a math error (i.e. negative value in SQRT)
min_valid = -1.0*const_factor/(photon_factor + spike_factor^2)
IF MIN(DN) LT min_valid THEN BEGIN
    loc_invalid = WHERE(DN LT min_valid)
    DN[loc_invalid] = min_valid
ENDIF

IF KEYWORD_SET(spike_noise) THEN BEGIN
    ;Includes the noise caused by the despiking process
    output_data = SQRT(photon_factor*DN + const_factor + (DN*spike_factor^2)) / SQRT(norm_num)
ENDIF ELSE BEGIN
    ;[DEFAULT] ignores despiking noise
    output_data = SQRT(photon_factor*DN + const_factor) / SQRT(norm_num)
ENDELSE

RETURN, output_data

END
