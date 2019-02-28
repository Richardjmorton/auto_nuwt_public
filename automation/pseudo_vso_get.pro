;+
;NAME: PSEUDO_VSO_GET
;
;PURPOSE:
;   This program attempts to download SDO / AIA files
;   The standard VSO_GET() function sometimes claims that a file URL is not
;   accessible even though the file can be manually downloaded just fine using
;   the very same URL that was attempted by VSO_GET!
;
;   This program attempts to workaround the issue by downloading the files using
;   very basic calls to the .get method of an 'IDLnetUrl' object.
;
;   WARNING: This program should only be used as a last resort. Use with care!
;
;INPUTS:
;   vso_index - index of files found using the VSO_SEARCH() function in SolarSoft
;
;OPTIONAL INPUTS:
;   out_dir - full_path to the folder where you wish to save the files.
;                 Defaults to the current working directory.
;
;OUTPUTS:
;   If .fits files are downloaded and saved to the specifiedfolder
;   Files that already exist on the local computer WILL BE OVERWRITTEN!
;
;NOTES:
;   Currently only downloads compressed files using the /rice keyword in VSO_GET()
;
;   The saved filesnames try to mimic the standard format of the filenames
;   downloaded using VSO_GET(). However, the timestamps are only given to the
;   nearest second, rather than ms. This is limited by the timestamps reported by
;   VSO_SEARCH().
;
;HISTORY: Name---------Date---------Description
;         M Weberg  18 MAY, 2018  Initial coding
;-


;THE CODE BELOW IS A MODIFIED VERSION OF THE EXAMPLE CODE IN THE IDLnetURL::Get DOCS ONLINE
; FUNCTION Url_Callback, status, progress, data
;     ; print the info msgs from the url object
;     PRINT, status
;
;     ; return 1 to continue, return 0 to cancel
;     RETURN, 1
; END

;-----------------------------------------------------------------

PRO PSEUDO_VSO_GET, vso_index, out_dir=out_dir

    IF NOT KEYWORD_SET(out_dir) THEN CD, CURRENT=out_dir ;i.e. defaults to current folder
    IF STRLEN(out_dir) GT 0 AND NOT out_dir.endswith('/') THEN out_dir = out_dir+'/'

    ; If the url object throws an error it will be caught here
    CATCH, errorStatus
    IF (errorStatus NE 0) THEN BEGIN
        CATCH, /CANCEL

        ; Display the error msg in a dialog and in the IDL output log
        r = DIALOG_MESSAGE(!ERROR_STATE.msg, TITLE='URL Error', /ERROR)

        PRINT, !ERROR_STATE.msg

        ; Get the properties that will tell us more about the error.
        oUrl->GetProperty, RESPONSE_CODE=rspCode, $
                           RESPONSE_HEADER=rspHdr, $
                           RESPONSE_FILENAME=rspFn

        PRINT, 'rspCode = ', rspCode
        PRINT, 'rspHdr= ', rspHdr
        PRINT, 'rspFn= ', rspFn

        ; Destroy the url object
        OBJ_DESTROY, oUrl
        RETURN
    ENDIF

    ;This is just used to get the URLs of the files we want to download
    vso_get_status = VSO_GET(vso_index, out_dir=full_save_path, /rice, /NODOWNLOAD)

    num_urls = N_ELEMENTS(vso_get_status)
    num_ind = N_ELEMENTS(vso_index)

    ; create a new IDLnetURL object
    oUrl = OBJ_NEW('IDLnetUrl')

    ; ; Specify the callback function
    ; oUrl->SetProperty, CALLBACK_FUNCTION ='Url_Callback'

    ; ; Set verbose to 1 to see more info on the transacton
    ; oUrl->SetProperty, VERBOSE = 0 ;WARNING: PRINTS A TON OF OUTPUT LINES!

    IF num_ind EQ num_urls THEN BEGIN
        FOR uu=0, (num_urls-1) DO BEGIN

            split_ID = STRSPLIT(vso_index[uu].fileid, '_:', /EXTRACT)
            inst_str = split_ID[0]
            lv_str = split_ID[1]
            wave_str = split_ID[2]
            ; wave_str = strtrim(string(vso_index[uu].wave.min, format='(i5.0)'), 2)
            ts_str = vso_index[uu].time.start.Replace(':', '_')
            save_filename = inst_str+'.'+lv_str+'.'+wave_str+'A_'+ts_str+'Z.image.'+lv_str+'.fits'

            data_url = vso_get_status[uu].URL ;url of the file to be downloaded

            ; Acutally download the file
            file_start_time = SYSTIME(/seconds)
            downloaded_file = oUrl->Get(FILENAME=out_dir+save_filename, URL=data_url)
            file_end_time = SYSTIME(/seconds)

            ; Report on the progress
            file_download_time = file_end_time - file_start_time
            download_time_str = STRTRIM(STRING(file_download_time, format='(f6.1)'), 2)
            PRINT, STRTRIM(uu, 2)+' : downloaded '+save_filename+' in '+download_time_str+' seconds'
        ENDFOR
    ENDIF ELSE BEGIN
        MESSAGE, 'ERROR: Number of URLs does not match the number of files found by VSO_SEARCH()'
    ENDELSE

    ; Destroy the url object
    OBJ_DESTROY, oUrl

END
