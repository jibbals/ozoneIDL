;   Procedure: 
;       plot the longitude( or latitude) slice of an event using high res ERA data
PRO event_lon, jday, $
  MELBOURNE=MELBOURNE, DAVIS=DAVIS, MACQUARIE=MACQUARIE, $
  LATITUDE=LATITUDE

  ; need lat/lon of chosen station
  if keyword_set(davis) then begin 
    loc=[77.967,-68.58]
  endif else $
  if keyword_set(macquarie) then begin
    loc=[158.967, -54.5]
  endif else $ 
    loc=[144.95, -37.69]
  iloc=1  ;latitude
  
  era=eradata(jday, MELBOURNE=MELBOURNE, DAVIS=DAVIS, MACQUARIE=MACQUARIE, /HIGHRES)

  n_levs=n_elements(era.pressure)
  ;n_levs=n_elements(era.geopotential)
  if keyword_set(LATITUDE) then begin 
    X=era.latitude
    iloc=0  ;longitude
    tmp=min(abs(era.longitude - loc[iloc]),Zind)
  endif else begin
    X=era.longitude
    tmp=min(abs(era.latitude - loc[iloc]),Zind)
  endelse 
  n_X=n_elements(X)
  Y=era.pressure
  Y2D=make_array([n_X,n_levs], /FLOAT, value=!values.f_nan)
  X2D=make_array([n_X,n_levs], /FLOAT, value=!values.f_nan)
  for i=0,n_levs-1 do X2D[*,i]=X
  for i=0, n_X-1 do Y2D[i,*]=Y
  xrange=[loc[1-iloc]-10.,loc[1-iloc]+10.]
  yrange=[1000,100]
  
  ;stop
  ;ct_load, 'calipso'
  loadct, 33
  cgdisplay, 800, 600, wid=0, title='base lon, lat:[ '+string(loc[0],format='(f6.2)')+', '+string(loc[1],format='(f6.2)')+' ]'
  ;!p.multi=[0,2,1]
  
  ;for i=2,3 do begin
  i=3
    CASE i OF
      0:BEGIN
          ARRAY=era.geopotential / 1.0e3            ;height(km)
          levels=findgen(18)/2.
          title='GPH'
        END
      1:BEGIN
          ARRAY=era.temperature                     ;temp(K)
          levels=findgen(10)*12.5+200.
          title='Temp(K)'
        END
      2:BEGIN
          ARRAY=era.potential_vorticity * 1.0E6     ;pvu
          levels=REVERSE(-1.0*findgen(8)/2.-0.5)
          title='PVU'
        END
      3:BEGIN
          ARRAY=era.ozone_mmr* 1.0E9                ;ppbv
          levels=findgen(12)*12.5+25
          title='PPBV O3'
        END
    ENDCASE
    IF keyword_set(LATITUDE) THEN $
      ARRAY=reform(ARRAY[Zind,*,*]) $
    ELSE $
      ARRAY=reform(ARRAY[*,Zind,*])
      
    cgcontour, ARRAY, X2D, Y2D, $
      /YLOG, $
      title=title, $
      xrange=xrange, yrange=yrange, $
      levels=levels, $
      /fill
    CGcontour, ARRAY, X2D, Y2D, $
      /YLOG, $
      levels=levels[0:-1:3], $
      /FOLLOW, /OVER
    
    ; overplot the 2pvu contour
    if i eq 2 then $
      CGcontour, ARRAY, X2D, Y2D, $
        /YLOG, $
        levels=[-2], $
        thick=2, $
        /OVER
    
    ; white line where station is located
    ;
    cgoplot, fltarr(2)+loc[1-iloc], yrange, color=cgcolor('white')
      
  ;endfor
  !p.multi=0
  
END
