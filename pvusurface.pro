PRO pvusurface, jday, $
  MELBOURNE=MELBOURNE, DAVIS=DAVIS, MACQUARIE=MACQUARIE

  ; get the data, we need high res for the surface plot
  era=eradata(jday, MELBOURNE=MELBOURNE, DAVIS=DAVIS, MACQUARIE=MACQUARIE, /HIGHRES)

  ; get the 2pvu surface to plot 
  help, era
  pvu = era.potential_vorticity * 1.0e6
  s=size(pvu, /dimension)
  pvu2 = fltarr([s[0],s[1]]) + !values.f_nan 
  
  ; couldn't get polyshade etc to work, just looking at top 2pvu surface for now
  foreach p, era.pressure, i do begin
    
    twos=where(pvu[*,*,i] gt -2.0 and not finite(pvu2)) 
    if twos[0] ne -1 then $
      pvu2[twos] = p
    
  endforeach

  ; set up graphics for colored surface
  setdecomposedstate, 0, currentstate=currentstate
  set_shading, values=[0,249]
  cgloadct, 5, /brewer, /reverse, ncolors=250
  
  cgsurface, pvu2, era.longitude, era.latitude, $
    zrange=[800, 100], $
    ztitle='pressure level', $
    xtitle='longitude', $
    ytitle='latitude', $
    ctable=5, $
    /shaded ; same as style=2

  setdecomposedstate, currentstate
  ;stop
END


;  X = era.longitude
;  Y = era.latitude
;  X2D = make_array(n_elements(X),n_elements(Y))
;  Y2D = make_array(n_elements(X),n_elements(Y))
;  for i=0,n_elements(Y)-1 do X2D[*,i] = X
;  for i=0,n_elements(X)-1 do Y2D[i,*] = Y
  

