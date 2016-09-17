;   PROCEDURE: eventweather
;
;   PURPOSE:  
;           Map the composite weather from the 'event days'
;
;   KEYWORDS:
;           MELBOURNE/MACQUARIE/DAVIS
;           
;   OPTIONAL:
;           TOPN=int  : only use top N events to make average picture 
PRO eventweather, $
  melbourne=melbourne, macquarie=macquarie, davis=davis, $
  topn=topn
  
  ; file for shortcut reading of era data events
  ; will eventually need options for mac and davis
  file='C:\Users\asp_transfer\Documents\Data\Broadmeadows\MelbEvents.dat'
  
  ; Get events from ozonesondes
  ; 
  data_sonde=sondedata(melbourne=melbourne, davis=davis, macquarie=macquarie)
  
  centre=[data_sonde.lon, data_sonde.lat]
  jdates=data_sonde.jtime
  years=JesseDate(jdates)
  
  getevents, data_sonde.o3ppbv, data_sonde.gph/1000.0, $
    data_sonde.tropopause, events=eventstruct
  
  ; dates of events
  eventinds=eventstruct.indices
  eventdates=jdates[eventinds]
  ; 'biggest' 20 events
  mags=eventstruct.magnitudes[eventinds]
  bigeventdays=eventdates[ reverse(sort(mags)) ]
  if n_elements(topn) gt 0 then bigeventdays=bigeventdays[0:topn-1] 
  bigeventinds=fltarr(n_elements(bigeventdays)) + !values.f_nan
  
  ; Get weather from the ERA datasets
  ;
;  if n_elements(topn) gt 0 then eventdates=jdates[bigeventinds[0:topn-1]]
;  data_era=eradata(eventdates,melbourne=melbourne, macquarie=macquarie, davis=davis)
  data_era=eradata(file=file)
  for i=0, n_elements(bigeventdays)-1 do begin
    tmp=min(abs(data_era.jultime-bigeventdays[i]),tmpind)
    if tmp lt 1 then $
    bigeventinds[i]=tmpind
  endfor
  bigeventinds=bigeventinds[where(finite(bigeventinds))]
  print, 'averaging ',total(finite(bigeventinds)), ' event weather maps'
  X = data_era.longitude
  Y = data_era.latitude
  Lims=[Y[0],X[0],Y[-1],X[-1]]
  
  ; get the averages over the event space
  ;
  avg_gph   =   mean(data_era.geopotential[*,*,*,bigeventinds], dimension=4, /nan) / 1.0e3
  avg_temp  =   mean(data_era.temperature[*,*,*,bigeventinds], dimension=4, /nan)
  avg_pv    =   -1 * mean(data_era.potential_vorticity[*,*,*,bigeventinds], dimension=4, /nan) * 1.0e6
  avg_ozone =   mean(data_era.ozone_mmr[*,*,*,bigeventinds], dimension=4, /nan) * 1.0e9
  titles=['gph(km)','temp(K)','pvu','ppbv o3']
  
  
  
  ; plot stuff like a boss
  ;
  CT_load, 'Calipso'
  cbarloc=[.045,.01, .055, .2]
  cbarlocdx=[.5, 0, .5, 0]
  cbarlocdy=[0,.5,0,.5]
  
  foreach p, data_era.pressure, pind do begin
    cgDisplay, 1200, 900, title="event mean weather for "+string(p)+"hPa", wid=pind
    !p.multi = [0,2,2]
    for i=0,3 do begin
      CASE i OF
        0: ARRAY=avg_gph                     ;height(km)
        1: ARRAY=avg_temp                    ;temp(K)
        2: ARRAY=avg_pv                      ;pvu
        3: ARRAY=avg_ozone                   ;ppbv
      ENDCASE
      
      ; Set up the map
      CGMAP_SET, Centre[1], Centre[0], /MERCATOR, LIMIT=Lims, advance=i, $
        title=titles[i]
      
      ; plot the filled in contour map
      CGcontour, ARRAY[*,*,pind], X, Y, nlev=50, /CELL, /OVER
      
      ; add some lines
      CGcontour, ARRAY[*,*,pind], X, Y, nlev=4, /FOLLOW, /OVER
      
      ; add in a colorbar
      cgColorbar, Divisions=3, range=[min(array[*,*,Pind],/nan),max(array[*,*,Pind],/nan)], $
        Position = cbarloc + ( i mod 2 ) * cbarlocdx + ( i lt 2 ) * cbarlocdy, $
        Title=titles[i], /vertical, $
        TLocation='right', Format='(f5.1)', $
        tcharsize=1., charsize=1.5
      
      ; add the continents
      CGMAP_CONTINENTS, thick=2,color=cgcolor('black')
    endfor
  endforeach
  !p.multi=0
END