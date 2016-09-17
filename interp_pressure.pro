; array: float[profiles, altitudes]
; pressure: float[profiles, altitudes]
; pvals: float[ pressures to interpolate to for each profile ]
function interp_pressure, array, pressure, pvals
  arraysize=size(array)
  N=arraysize[1]
  M=n_elements(pvals)
  
  retarr=MAKE_ARRAY(N,M,/DOUBLE,VALUE=0d)
  foreach P, pvals, Pind do begin
    for i=0, N-1 do begin
      ; last element with more pressure
      below=max(where(pressure[i,*] ge P))
      ; first element with less pressure 
      above=min(where(pressure[i,*] le P))
      ; Pressure below 200hpa(will be more than 200hpa)
      a=pressure[i,below]
      ; Pressure above 200hpa(will be less than 200hpa)
      b=pressure[i,above]
      
      ; don't interpolate too far 
      if a-b gt 75 then begin
        retarr[i,Pind]=!values.f_nan
        continue
      endif
      
      if a eq b then $
        retarr[i, Pind] = array[i,below]
      if a ne b then $
        retarr[i,Pind] = array[i,below] * ( 1 - (a-P)/double(a-b) ) + $
        array[i,above] * ( 1 - (P-b)/double(a-b) )
    endfor
  endforeach

  return, retarr

end