;+
; NAME:
;    GETINTERPOLATEDDAVIS
;
; Purpose:
;    Restore the davis sondes data file and return the ppmv ozone interpolated at 
;    any pressure(s) 
;
;
; Calling Sequence:
;
;    OzonePPMV = GetInterpolatedData([pressure1,p2,p3,...],[julian start, julian end]) 
;
; Inputs:
;    Pressures: you want to interpolate to,
;    Timelimits: [start, end] date (Julian) which you want to limit your data set to
;
; Optional Parameters:
;    DAVISFILE: where is the ozonesondes.dat 
;
; Notes:
;    This should be updated to read a netcdf rather than the saved variable output
;    
; Example:
;
;          
;
; MODIFICATION HISTORY:
;    7-Jan-2015: Created by Jesse Greenslade
;-

FUNCTION GetInterpolatedDavis, Pressures, TimeLimits=TimeLimits, DAVISFILE=DAVISFILE

  ; Set Defaults
  if n_elements(DAVISFILE) eq 0 then $
    DAVISFILE="C:\Users\asp_transfer\Documents\Data\Davis\ozonesondes.dat"
  if n_elements(TimeLimits) eq 0 then $
    TimeLimits=[-1, 1e10] ; should cover all time
  M=N_ELEMENTS(Pressures)

  ; Read the davis data
  ;
  RESTORE, DAVISFILE   
  
  ; Limit time returned using the TimeLimits parameter
  ;
  dateinds=where(OZONE_JTIME gt TimeLimits[0] and OZONE_JTIME lt TimeLimits[1])
  dates=OZONE_JTIME[dateinds]
  N=N_ELEMENTS(dates)
  
  ; Get the interpolated ozone ppmv for each pressure at every time within the time limits
  ;
  vmr=MAKE_ARRAY(N,M,/DOUBLE,VALUE=0d)
  foreach P, Pressures, Pind do begin
    foreach i,dateinds, iind do begin
      ; last element with more pressure
      below=max(where(OZONE_PRESSURE[i,*] ge P))
      ; first element with less pressure 
      above=min(where(OZONE_PRESSURE[i,*] le P))
      ; Pressure below 200hpa(will be more than 200hpa)
      a=OZONE_PRESSURE[i,below]
      ; Pressure above 200hpa(will be less than 200hpa)
      b=OZONE_PRESSURE[i,above]
      vmr[iind,Pind] = OZONE_MR[i,below] * ( 1 - (a-P)/double(a-b) ) + $
        OZONE_MR[i,above] * ( 1 - (P-b)/double(a-b) )
    endforeach
  endforeach
  
  return, {ppmv:vmr, dates:dates}
END