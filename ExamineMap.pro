;   Procedure: 
;       ExamineMap
;       
;   Purpose:
;       Look at the weather at and +- 6 hours of a particular julian time. weather maps are at 200, 300, 500hpa
;       
;   Example:
;       ExamineMap, jday=julday(5,15,2010), /Davis
;
PRO ExamineMap, jday, $
  imageprefix=imageprefix,$
  melbourne=melbourne, davis=davis, macquarie=macquarie, highres=highres

  ; get the ERA Data
  era=eradata(jday + [-.25, 0, .25], melbourne=melbourne, davis=davis, macquarie=macquarie, highres=highres)

  ; We want to look at a (single time and pressure) map of winds and geopotential height
  ; 
  ret = Min(Abs(era.jultime - jday), jindex)
  pinds= where(era.pressure eq 200 or era.pressure eq 300 or era.pressure eq 500)

  if ret gt 1. then print, ret, "days from requested date"
  X = era.longitude
  Y = era.latitude
  Lims=[Y[0],X[0],Y[-1],X[-1]]
  Center=[(Y[0]+Y[-1]) / 2., (X[0]+X[-1])/2.]
  X2D = make_array(n_elements(X),n_elements(Y))
  Y2D = make_array(n_elements(X),n_elements(Y))
  for i=0,n_elements(Y)-1 do X2D[*,i] = X
  for i=0,n_elements(X)-1 do Y2D[i,*] = Y
  ; ERA arrays we will examine
  jindex=jindex+[-1,0,1]
  
  ERAppbv = era.ozone_mmr[*,*,*,jindex]*1.0e9  

  zonal = era.zonal_velocity
  meridional = era.meridional_velocity
  vertical = era.vertical_velocity
  
  ; Plot stuff
  ;
  ytitles = ['200 hpa','300 hpa','500 hpa']
  caldat, jday, MM,DD,YY,HH
  YMDH= string(YY,format='(i5)')+string(MM,format='(i3)')+$
    string(DD,format='(i3)')+string(HH,format='(i3)')
  xtitles = [' - 6 hours','Y M D H:'+YMDH , '+6 hours']
  !Y.Margin=[2,8]
  
  ;CT_LOAD, 'CALIPSO'
  cgloadct, 33
  ;===================================================
  ;Data achieved, commence plotting
  ;===================================================
  cx=0.025 & cxd=0.01 & cy=0.015 & cyd=.15 & cyn=.33
  cbarpositions=fltarr(4,3)
  foreach k,[2,1,0],ki do $
    cbarpositions[*,ki]= [cx, cy, cx, cy] + k*[0, cyn, 0, cyn] + [0,0,cxd, cyd]
  ctitles=['GPH', 'Temp', 'PVU', 'ppbv']
  
  ;plotting gph, temp, pvu, and ozone
  FOR plotcount=0, 3 DO BEGIN
  
    ; save the image output
    YMDH = string(YY,format='(i4)')+string(MM,format='(i02)')+$
      string(DD,format='(i02)')+string(HH,format='(i02)')
  
    if n_elements(imageprefix) gt 0 then $
      cgps_open, filename=imageprefix+ctitles[plotcount]+YMDH+'.ps'
  
    CASE plotcount OF
      0: ARRAY=ERA.GEOPOTENTIAL/1000.0            ;height(km)
      1: ARRAY=ERA.TEMPERATURE                    ;temp(K)
      2: ARRAY=abs(ERA.POTENTIAL_VORTICITY*1.0E6) ;abs pvu
      3: ARRAY=o3mmrtovmr(ERA.OZONE_MMR) * 1.0E9  ;ppbv
    ENDCASE
    
    !P.multi = [0,3,3] ; 3 rows 3 columns
    
    cgDisplay, 1400, 960, title="Julian Day: "+string(jday,format='(f10.1)')+"+[-6,0,6] hours", $
      wid=plotcount
    
    ; for each pressure level
    FOREACH P, pinds,Pi DO BEGIN
      ; plot the -6 hours map, and the +6 hours map on the left and right respectively
      FOREACH J, jindex,Ji DO BEGIN
        
        ; Titles along top and down the left
        title=''
        if Pi eq 0 then title=xtitles[Ji] 
        if Ji eq 0 then title=ytitles[Pi]
        if Pi+Ji eq 0 then title=xtitles[0]+' at '+ytitles[0]
        
        ; Set up the map for our contour plots
        CGMAP_SET, Center[0], Center[1],/MERCATOR, LIMIT=Lims, advance=Ji+Pi, $
          title=title, charsize=2
        
        ; plot the filled in contour map
        CGcontour, ARRAY[*,*,P, J], X, Y, nlev=50, /CELL, /OVER
        
        ; for pvu's add a line at pvu=1,2,3
        if plotcount eq 2 then $
          CGcontour, ARRAY[*,*,P, J], X, Y, levels=[2], /FOLLOW, /OVER
        
        ; add in a colorbar
        if Ji eq 0 then $
          cgColorbar, Divisions=3, range=[min(array[*,*,P,J],/nan),max(array[*,*,P,J],/nan)], $
            Title=ctitles[plotcount], /vertical, $
            TLocation='right', Format='(f5.1)', Position=cbarpositions[*,Pi], $
            tcharsize=1., charsize=1.5
        
        ; vector lines for horizontal winds
        if keyword_set(highres) then frac=0.015 else frac= 0.1
        if plotcount eq 0 then begin
          cgdrawvectors, zonal[*,*,P,J],meridional[*,*,P,J], X2D,Y2D, $
            fraction=frac, /ORDERED, /OVERPLOT, HSIZE=(!D.X_SIZE / 150), $
            color=cgcolor('black'), thick=0.5, length=.05
          
          ; add contour for vertical wind velocity
          CGcontour, vertical[*,*,P,J], X, Y, /OVERPLOT, $
            levels=[-.4,-.2,0.,.2,.4], c_linestyle=[1,1,0,2,2], $
            c_thick=[1.5,1.,1.,1.,1.5], c_color=[2,2,1,3,3]
        endif
                  
        ; add the continents
        CGMAP_CONTINENTS, thick=2,color=cgcolor('black')
        
        ; put a star where the station is
        if keyword_set(melbourne) then $
          cgplots, 144.95, -37.69, psym=2, symsize=2, thick=2 $
        else if keyword_set(davis) then $
          cgplots, 77.9674, -68.5766, psym=2, symsize=2, thick=2 $
        else $  ; macca
          cgplots, 159.0, 54.5, psym=2, symsize=2, thick=2
        
      ENDFOREACH
    ENDFOREACH
    !P.multi = 0
    
    if n_elements(imageprefix) gt 0 then $
      cgps_close
    
    
    
  ENDFOR
END
