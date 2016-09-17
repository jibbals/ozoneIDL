;   Procedure: 
;       EventMaps
;       
;   Purpose:
;       Look at the weather maps of the 9 or so most perturbed events
;       
;   Example:
;       EventMaps, /melbourne
;
PRO EventMaps, plevel=plevel, $
  imageprefix=imageprefix,$
  melbourne=melbourne, davis=davis, macquarie=macquarie

  if n_elements(plevel) eq 0 then plevel=300

  ; Get events from ozonesondes
  ; 
  data_sonde=sondedata(melbourne=melbourne, davis=davis, macquarie=macquarie)
  jdates=data_sonde.jtime
  
  getevents, data_sonde, events=eventstruct
  
  ; dates of events
  eventinds=eventstruct.indices
  eventdates=jdates[eventinds]

  ; get the ERA Data
  ;
  era=eradata(eventdates, davis=davis,melbourne=melbourne,macquarie=macquarie)
    
  ; We want to look at a single time and pressure map of winds
  ; 
  pind= where(era.pressure eq plevel)

  X = era.longitude
  Y = era.latitude
  Lims=[Y[0],X[0],Y[-1],X[-1]]
  Center=[(Y[0]+Y[-1]) / 2., (X[0]+X[-1])/2.]
  X2D = make_array(n_elements(X),n_elements(Y))
  Y2D = make_array(n_elements(X),n_elements(Y))
  for i=0,n_elements(Y)-1 do X2D[*,i] = X
  for i=0,n_elements(X)-1 do Y2D[i,*] = Y
  
  
  ERAgph =ERA.GEOPOTENTIAL/1000.0            ;height(km)
  ERApvu = -1* ERA.POTENTIAL_VORTICITY * 1.0e6   ;pvu
  zonal = era.zonal_velocity
  meridional = era.meridional_velocity
  vertical = era.vertical_velocity
  
  ; 'biggest' n events from ozonesonde data
  ;
  n_events=9
  mags=eventstruct.magnitudes[eventinds]
  bigeventdays=eventdates[ reverse(sort(mags)) ]
  bigsondeinds=float(eventinds[ reverse(sort(mags)) ])
  bigeventdays=bigeventdays[0:n_events-1] 
  bigsondeinds=bigsondeinds[0:n_events-1]
  bigerainds=fltarr(n_elements(bigeventdays)) + !values.f_nan
  for i=0, n_elements(bigeventdays)-1 do begin
    tmp=min(abs(era.jultime-bigeventdays[i]),tmpind)
    if tmp lt 1 then $
      bigerainds[i]=tmpind $
    else $
      bigsondeinds[i]=!values.f_nan
  endfor
  bigerainds=bigerainds[where(finite(bigerainds))]
  bigsondeinds=bigsondeinds[where(finite(bigsondeinds))]

  ;===================================================
  ;Data achieved, commence plotting
  ;===================================================
  CT_LOAD, 'CALIPSO'
  if n_elements(imageprefix) gt 0 then $
    cgps_open, filename=imageprefix+'eventmaps.ps'
  
  cgDisplay, 1400, 960
  !P.multi = [0,3,3] ; 3 rows 3 columns
  
  FOREACH J, bigerainds, Ji DO BEGIN

    caldat, era.jultime[J], MM,DD,YY,HH  
    YMDH = string(YY,format='(i4)')+string(MM,format='(i02)')+$
      string(DD,format='(i02)')+string(HH,format='(i02)')

    print, Ji,' : ', YMDH
  
    ; Set up the map for our contour plots
    CGMAP_SET, Center[0], Center[1],/MERCATOR, LIMIT=Lims, advance=Ji, $
        title=YMDH, charsize=2
    
    ; plot the filled in contour map
    CGcontour, ERAGPH[*,*,Pind, J], X, Y, nlev=50, /CELL, /OVER
      
    ; vector lines for horizontal winds
    cgdrawvectors, zonal[*,*,Pind,J], meridional[*,*,Pind,J], X2D,Y2D, $
      fraction=0.1, /ORDERED, /OVERPLOT, HSIZE=(!D.X_SIZE / 150), $
      color=cgcolor('black'), thick=0.5, length=.05
      
    ; add contour for vertical wind velocity
    CGcontour, ERApvu[*,*,Pind,J], X, Y, /OVERPLOT, $
      levels=[1,2,3], /FOLLOW, $
      c_thick=[1,2,1]
                
    ; add the continents
    CGMAP_CONTINENTS, thick=3,color=cgcolor('dark grey')
      
    ; put a star where the station is
    if keyword_set(davis) then $
      cgplots, 77.9674, -68.5766, psym=2, symsize=2, thick=2 $
    else if keyword_set(macquarie) then $
      cgplots, 158.967, -54.5, psym=2, symsize=2, thick=2 $
    else $; melbourne
      cgplots, 144.95, -37.69, psym=2, symsize=2, thick=2
    
  ENDFOREACH
  !P.multi = 0
  
  if n_elements(imageprefix) gt 0 then $
    cgps_close
  
  ;==================================================
  ;             Plot N biggest profiles
  ;==================================================
  cgdisplay, 1600, 800, wid=1, title='SomeEvents'
  N=n_elements(bigerainds)
  
  !P.multi=[0,N,1]
  foreach i,bigsondeinds,ii do begin
 
    caldat, data_sonde.jtime[i], M, D, Y, H
    fi2='(i02)' & fi4='(i04)'
    title=string(Y, format=fi4)+string(M,format=fi2)+string(D,format=fi2)+$
      string(H,format=fi2)
    print, ii, ' : ', title
    
    ppbv=reform(data_sonde.o3ppbv[i,*])
    temp=reform(data_sonde.temperature[i,*])+273.15
    RH = reform(data_sonde.RH[i,*])
    gph=reform(data_sonde.gph[i,*]/1000.0)
    tpo3=data_sonde.tpo3[i]
    tplr=data_sonde.tplr[i]
    
    plotprofile, ppbv, gph, tpo3=tpo3, tplr=tplr, title=title, $
      /altitude, top=tpo3+5, temp=temp, Humid=RH
    
  endforeach
  !p.multi=0
    
END