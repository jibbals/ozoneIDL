; NAME: TemperatureTropopause
; 
; Purpose: 
;         Return altitude of tropopause using temperature and altitude inputs
; 
; Inputs:
;         temperature: ozone ppbv array
;         Z: altitudes(km) array corresponding to the ozone ppbv array
; Outputs:
;         Index: index of the tropopause height in the Z input
;         
; Keywords:
;         Polar: if we are in the polar regions there is a slightly altered
;         tropopause definition (only for ozone?)
; 
; Returns:
;         Tropopause altitude(float)
;         third argument will match the index of the tropopause
;         
; Example:
;         TP = TemperatureTropopause(data.temperature, data.altitude, tpIndex)
;         ; TP will be the tropopause height, tpIndex will be it's index in the ozone.altitude array
;
FUNCTION temperaturetropopause, temp, Z, Index=Index
  
  ; lapse rate tropopause
  rate = -2.0 ; minimum dt/dz 
  minh = 2.0 ; look above this height only(avoid ground effects)
  dh = 2.0 ; average dt/dz over what height needs to beat the rate?

  ; first find dt/dz
  lapse=deriv(Z,temp)
  if max(Z,/nan) lt 5 then return, !values.f_nan
  
  ; find where lapse rate gt -2K/km
  testinds=where(lapse gt rate and Z gt minh)
  alt=0.
  foreach ind, testinds, i do begin
    ; if the lapse rate for the next 2k is over -2K/km then we are done
    checks=where(Z ge Z[ind] and Z le Z[ind]+dh)
    if mean(lapse[checks]) gt rate then begin
      alt=Z[ind]
      break
    endif   
  endforeach
  if alt eq 0. then return, !values.f_nan
  Index=ind
  return, alt
  
END
