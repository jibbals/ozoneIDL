;   Procedure: 
;       ExamineMap
;
;   Purpose:
;       Look at the wind and gph map at one jday and pressure level
;       
;   Example:
;       windmap, julday(2,3,2005,0),/highres, /melbourne
;           shows GPH and wind map
;       windmap, julday(1,13,2010), /highres, /melb, maptype=2
;           shows PVU map
;
;   Prerequisites:
;       coyotegraphics
;       eradata script

PRO windmap, jday, $
  imageprefix=imageprefix, plev=plev, maptype=maptype,$
  melbourne=melbourne, davis=davis, macquarie=macquarie, highres=highres

  ; get the ERA Data
  era=eradata(jday, $
    melbourne=melbourne, davis=davis, macquarie=macquarie, highres=highres)

  ; We want to look at a (single time and pressure) map of winds and geopotential height
  ; 
  if n_elements(plev) eq 0 then plev = 500
  if n_elements(maptype) eq 0 then maptype=0

  ret = Min(Abs(era.jultime - jday), jindex)
  pindex= where(era.pressure eq plev)
  if ret gt 1. then print, ret, "days from requested date"

  ; set up map coordinates
  ;
  X = era.longitude
  Y = era.latitude
  Lims=[Y[0],X[0],Y[-1],X[-1]]
  Center=[(Y[0]+Y[-1]) / 2., (X[0]+X[-1])/2.]
  X2D = make_array(n_elements(X),n_elements(Y))
  Y2D = make_array(n_elements(X),n_elements(Y))
  for i=0,n_elements(Y)-1 do X2D[*,i] = X
  for i=0,n_elements(X)-1 do Y2D[i,*] = Y
  

  if keyword_set(melbourne) then $
    plot_title='Melbourne ' $
  else if keyword_set(davis) then $
    plot_title='Davis ' $
  else $  ; macca
    plot_title='Macquarie '


  ; velocities
  ;
  xvel = era.zonal_velocity[*,*,pindex]
  yvel = era.meridional_velocity[*,*,pindex]
  gph  = ERA.GEOPOTENTIAL[*,*,pindex]/1000.0 ; height(km)
  pvu  = abs(era.potential_vorticity[*,*,pindex]) * 1e6         ; abs pvu

  
  ; Plot stuff
  ;
  ytitle = 'Latitude'
  caldat, jday, MM,DD,YY,HH
  YMDH= string(YY,format='(i5)')+string(MM,format='(i3)')+$
    string(DD,format='(i3)')+string(HH,format='(i3)')
  xtitle ='Longitude'
  plot_title=plot_title+YMDH+'(YMDH)'
  ;!Y.Margin=[2,8]
  
  ; load the color tables(cg library script)
  ;CGLOADCT, 33
  loadct, 33
  
  ;===================================================
  ;Data achieved, commence plotting
  ;===================================================
  cbarposition= [0.09, 0.05, 0.11, 0.4]
  ctitle='GPH(km)'

  
  ; save the image output
  YMDHstr = string(YY,format='(i4)')+string(MM,format='(i02)')+$
    string(DD,format='(i02)')+string(HH,format='(i02)')


  cgDisplay, 1200, 800, title="Julian Day: "+string(jday,format='(f10.1)'), $
    wid=0
  
  ; Set up the map for our contour plots
  CGMAP_SET, Center[0], Center[1],/MERCATOR, LIMIT=Lims, $
    charsize=2, /noborder
  
  ; plot the filled in contour map
  CGcontour, gph, X, Y, title=plot_title, $
    nlev=50, /CELL, /OVER
        
  ; Add a line at pvu=1,2,3
  CGcontour, pvu, X, Y, levels=[1,2], /FOLLOW, /OVER, $
    c_thick=[2,3], c_colors=[cgcolor('purple'),cgcolor('white')]
  
  ; vector lines for horizontal winds
  ; draw wind vectors for gph map
  if keyword_set(highres) then frac=0.015 else frac= 0.1
  ;cgdrawvectors, xvel, yvel, $
  ;  X2D, Y2D, fraction=frac, /ORDERED, /OVERPLOT, HSIZE=!D.X_SIZE / 150, $
  ;  color=0, thick=0.5, length=.05
  velocity_field, xvel,yvel, x, y, $
    HSIZE=!D.X_SIZE / 150
      
        
  ; add the continents
  CGMAP_CONTINENTS, thick=2,color=cgcolor('black')
        
  ; put a star where the station is
  if keyword_set(melbourne) then $
    cgplots, 144.95, -37.69, psym=2, symsize=2, thick=2 $
  else if keyword_set(davis) then $
    cgplots, 77.9674, -68.5766, psym=2, symsize=2, thick=2 $
  else $  ; macca
    cgplots, 159.0, 54.5, psym=2, symsize=2, thick=2
 
  ; add in a colorbar
  cgColorbar, Divisions=3, $
    range=[min(gph,/nan), max(gph,/nan)], $
    Title=ctitle, /vertical, $
    TLocation='right', Format='(f5.1)', Position=cbarposition, $
    tcharsize=1., charsize=1.5

END
