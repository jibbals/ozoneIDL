;   PROCEDURE: CompositeWeather
;
;   PURPOSE:  
;           Map the composite weather of N events sorted by type
;
;   KEYWORDS:
;           MELBOURNE/MACQUARIE/DAVIS
;           
;   OPTIONAL:
;           N=int  : only use N events to make average picture 
PRO CompositeWeather, $
  melbourne=melbourne, macquarie=macquarie, davis=davis, $
  n=n
   
  ntypes=3
  
  ; Get events from ozonesondes
  ; 
  data_sonde=sondedata(melbourne=melbourne, davis=davis, macquarie=macquarie)
  
  centre=[data_sonde.lon, data_sonde.lat]
  jdates=data_sonde.jtime
  years=JesseDate(jdates)
  
  getevents, data_sonde.o3ppbv, data_sonde.gph/1000.0, $
    data_sonde.tpo3, events=ev
  eventdates=jdates[ev.indices]
  
  ; Get weather from the ERA datasets
  ;
  data_era = ptrarr(3, /allocate_heap)
  for i=0,ntypes-1 do begin
    dates=jdates[where(ev.type eq i)] ; dates of events of type i
    data_era[i]=ptr_new(eradata(dates,melbourne=melbourne, macquarie=macquarie, davis=davis))
  endfor
  
  print, 'count of ERA maps to be averaged: '
  count=fltarr(ntypes)
  for i=0,ntypes-1 do begin
    count[i]=total(finite((*data_era[i]).jultime))
    print, count[i]
  endfor
    
  ; plot stuff like a boss
  ;
  CT_load, 'Calipso'
  cbarloc=[.045,.01, .055, .2]
  cbarlocdx=[.5, 0, .5, 0]
  cbarlocdy=[0,.5,0,.5]
  titles=['gph(km)','temp(K)','pvu','ppbv o3']
  X = (*data_era[0]).longitude
  Y = (*data_era[0]).latitude
  Lims=[Y[0],X[0],Y[-1],X[-1]]
  
  pind = where((*data_era[0]).pressure eq 300.0)
  ; for each type
  ; 
  for type=0,ntypes-1 do begin
    
    ; plot 4 averages unless there is nothing to average
    ;
    if count[type] le 1 then continue
    
    cgDisplay, 1200, 900, title="event mean weather for 300 hPa, type="+string(type), wid=type
    !p.multi = [0,2,2]
    
    for i=0,3 do begin
      CASE i OF
        ;height(km)
        0: ARRAY=mean((*data_era[type]).geopotential, dimension=4, /nan) / 1.0e3
        ;temp(K)
        1: ARRAY=mean((*data_era[type]).temperature, dimension=4, /nan)
        ;pvu
        2: ARRAY=-1 * mean((*data_era[type]).potential_vorticity, dimension=4, /nan) * 1.0e6
        ;ppbv
        3: ARRAY=o3mmrtovmr(mean((*data_era[type]).ozone_mmr, dimension=4, /nan)) * 1.0e9
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
  endfor
  !p.multi=0
END