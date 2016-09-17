;   Procedure: 
;       makemap
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
;       eradata2 script...

PRO makemap, jday, $
  imageprefix=imageprefix, plev=plev, maptype=maptype,$
  melbourne=melbourne, davis=davis, macquarie=macquarie, highres=highres

  ; get the ERA Data
  caldat, jday,themonth,theday,theyear
  if theyear lt 2010 then $
    era=eradata(jday, $
      melbourne=melbourne, davis=davis, macquarie=macquarie, highres=highres) $
  else $
    era=eradata2(jday, melbourne=melbourne, davis=davis, macquarie=macquarie)
  ; We want to look at a (single time and pressure) map of winds and geopotential height
  ; 
  if n_elements(plev) eq 0 then plev = 500
  if n_elements(maptype) eq 0 then maptype=0
  
  ; Get closest matching day to the input jday
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
  xvel = era.zonal_velocity[*,*,pindex,jindex]
  yvel = era.meridional_velocity[*,*,pindex,jindex]
  zvel = era.vertical_velocity[*,*,pindex,jindex]
  
  ; Plot stuff
  ;
  ytitle = 'Latitude'
  caldat, jday, MM,DD,YY,HH

  YMDHstr = string(YY,format='(i4)')+string(MM,format='(i02)')+$
    string(DD,format='(i02)')+string(HH,format='(i02)')
  if long(total(finite(era.jultime))) eq 0L then begin
    ; no matches..
    get_lun,lun
    print, "Missing data for ", YMDHstr
    OpenW, lun, imageprefix+"Missing_"+YMDHstr
    printf, lun, "missing"
    Close, lun
    free_lun, lun
    return
  endif
  
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
  ctitles=['GPH(km)', 'Temp(K)','PVU', 'Ozone(ppbv)']
  ctitles2=['GPH','Temp','PVU','Ozone']
  ctitle=ctitles[maptype]
  ctitle2=ctitles2[maptype]
  
  ; save the image output

  ; if image prefix is set then save the image  
  if n_elements(imageprefix) gt 0 then $
    cgps_open, filename=imageprefix+ctitle2+YMDHstr+'.ps'

  CASE maptype OF
    0: ARRAY=ERA.GEOPOTENTIAL[*,*,pindex,jindex]/1000.0            ;height(km)
    1: ARRAY=ERA.TEMPERATURE[*,*,pindex,jindex]                   ;temp(K)
    2: ARRAY=abs(ERA.POTENTIAL_VORTICITY[*,*,pindex,jindex]*1.0E6) ;abs pvu
    3: ARRAY=o3mmrtovmr(ERA.OZONE_MMR[*,*,pindex,jindex]) * 1.0E9  ;ppbv
  ENDCASE
    
  cgDisplay, 1200, 800, title="Julian Day: "+string(jday,format='(f10.1)'), $
    wid=0
  
  ; Set up the map for our contour plots
  CGMAP_SET, Center[0], Center[1],/MERCATOR, LIMIT=Lims, $
    charsize=2, /noborder
  
  ; plot the filled in contour map
  CGcontour, ARRAY, X, Y, title=plot_title, $
    nlev=100, /CELL, /OVER
        
  ; for pvu's add a line at pvu=1,2,3
  if maptype eq 2 then $
    CGcontour, ARRAY, X, Y, levels=2, /FOLLOW, /OVER
        
  
  ; vector lines for horizontal winds
  ; draw wind vectors for gph map
  if keyword_set(highres) then frac=0.015 else frac= 0.1
  if maptype eq 0 then begin
  ;  cgdrawvectors, xvel, yvel, $
  ;    X2D, Y2D, fraction=frac, /ORDERED, /OVERPLOT, HSIZE=!D.X_SIZE / 150, $
  ;    color=0, thick=0.5, length=.05
    xvelnew=xvel
    yvelnew=yvel
    xnew=X
    ynew=Y
    if theyear gt 2009 then begin
      xvelnew=xvel[0:*:4,0:*:3]
      yvelnew=yvel[0:*:4,0:*:3]
      xnew=X[0:*:4]
      ynew=Y[0:*:3]
    endif
    velocity_field, xvelnew,yvelnew, xnew, ynew, $
      HSIZE=!D.X_SIZE / 150
      
    ; add contour for vertical wind velocity
    ;CGcontour, zvel, X, Y, /OVERPLOT, $
    ;  levels=[-.4,-.2,0.,.2,.4], c_linestyle=[1,1,0,2,2], $
    ;  c_thick=[1.5,1.,1.,1.,1.5], c_color=[2,2,1,3,3]
  endif
                  
  ; add the continents
  CGMAP_CONTINENTS, thick=2,color=cgcolor('black')
  ; put a star where the station is
  if keyword_set(melbourne) then $
    cgplots, 144.95, -37.69, psym=2, symsize=2, thick=2 $
  else if keyword_set(davis) then $
    cgplots, 77.9674, -68.5766, psym=2, symsize=2, thick=2 $
  else $  ; macca
    cgplots, 159.0, -54.5, psym=2, symsize=2, thick=2
 
  ; add in a colorbar
  cgColorbar, Divisions=3, $
    range=[min(ARRAY,/nan), max(ARRAY,/nan)], $
    Title=ctitles[maptype], /vertical, $
    TLocation='right', Format='(f5.1)', Position=cbarposition, $
    tcharsize=1., charsize=1.5

  if n_elements(imageprefix) gt 0 then $
    cgps_close
  
END
