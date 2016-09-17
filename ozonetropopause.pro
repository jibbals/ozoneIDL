; NAME: OzoneTropopause
; 
; Purpose: 
;         Return altitude of tropopause using ozone ppbv and altitude inputs
; 
; Inputs:
;         ppbv: ozone ppbv array
;         Z: altitudes(km) array corresponding to the ozone ppbv array
; Outputs:
;         Index: index of the tropopause height in the Z input
;         
; Keywords:
;         Polar: if we are in the polar regions there is a slightly altered
;         ozone tropopause definition
; 
; Returns:
;         Tropopause altitude(float)
;         third argument will match the index of the tropopause
;         
; Example:
;         TP = OzoneTropopause(ozone.ppbv, ozone.altitude, tpIndex)
;         ; TP will be the tropopause height, tpIndex will be it's index in the ozone.altitude array
;
FUNCTION ozonetropopause, ppbv, Z, Index=Index, POLAR=POLAR
  
  ; we want the lowest altitude with dmr/dz > 60ppbv/km where altitude+.5 up to +1.5 is > 110ppbv
  
  ; first find indices of dmr/dz > 60
  ppbv1=[ppbv, 0]
  z1 = [Z, 500]
  dmr=ppbv1[1:-1]-ppbv
  dz=z1[1:-1] - z
  gt60 = where(dmr/dz gt 60)
  gt80 = where(ppbv GT 80)
  gt110 = where(ppbv GT 110)
  testinds=cgSetIntersection(gt60,gt80)
  
  upper=2.0
  IF KEYWORD_SET(POLAR) THEN $
    upper = 1.5
  
  ; for each vmr gt 60 check the .5 to 1.5 k higher values are gt 110
  FOREACH ind, testinds DO BEGIN
    alt=Z[ind]
    checks = where(Z GE alt+0.5 AND Z LE alt+upper)
    inter = cgSetIntersection(checks, gt110)
    
    ; if all the indices in our check range are gt 110 ppbv we are finished
    IF n_elements(checks) EQ n_elements(inter) THEN BREAK
    
  ENDFOREACH
  
  Index=ind
  return, alt
  
END
