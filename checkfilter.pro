;
PRO checkfilter ,melbourne=melbourne, davis=davis, macquarie=macquarie
  
  ; read in the data
  data=sondedata(melbourne=melbourne, davis=davis, macquarie=macquarie)  
  n_sondes=n_elements(data.jtime)
  sondetps=min([[data.tplr], [data.tpo3]], dimension=2,/nan)
  sondetps[where(sondetps lt 3)] = !values.f_nan
  
  events=getevents(data)
  jtimes=events.jtime
  dates=jessedate(jtimes)
  ppbvs=events.o3ppbv
  n_events=n_elements(jtimes)
  
  print, n_sondes, ' ozonesondes'
  print, long(total(finite(data.gph[*,0]))), ' ozonesondes with gph'
  print, n_events, ' events'
  
  locs=events.locations
  gphs=events.gph / 1.0e3
  tps=events.tp
  caldat, jtimes, Months, Days, Years, Hours
   
  stationstr='Melbourne'
  if keyword_set(macquarie) then stationstr='Macquarie'
  if keyword_set(davis) then stationstr='Davis'
  summarytitle=stationstr+'('+string(n_events, format='(i4)')+' events)'

  ;==================================================
  ;             Plot 8 example profiles
  ;==================================================
  cgdisplay, 1500, 800, wid=0, title='SomeEvents'
  
  !P.multi=[0,8,1]
  fi2='(i02)' & fi4='(i04)'
  for i=0,7 do begin
    title='YMDH: ' + string(Years[i], format=fi4)+string(Months[i],format=fi2)+string(Days[i],format=fi2)+$
      string(Hours[i],format=fi2)
    temp=reform(events.temperature[i,*])+273.15
    RH = reform(events.RH[i,*])
    ppbv=reform(ppbvs[i,*])
    gph= reform(gphs[i,*])
    
    plotprofile, ppbv, gph, tpo3=events.tpo3[i], title=title, $
      /altitude, top=14, temp=temp, Humid=RH, tplr=events.tplr[i]
    
  endfor
  !p.multi=0
  
  ;==================================================
  ;             Event Time Plots
  ;==================================================
  cgdisplay, 1000, 800, wid=1, title='EventSummary'
  
  cgplot, dates, tps, $
    POSITION=[0.10, 0.55, 0.65, 0.90], $
    psym=1, $
    XTITLE='Time (years)', $
    YTITLE='Altitude (kms)', $
    TITLE='Event TP Heights'
  cgoplot, [min(dates,/nan),max(dates,/nan)], fltarr(2)+mean(tps,/nan), $
    linestyle=2, color=cgcolor('blue')
  cgoplot, [min(dates,/nan),max(dates,/nan)], fltarr(2)+mean(sondetps,/nan), $
    linestyle=2, color=cgcolor('black')
  cglegend, titles=['events','overall average'], colors=['black','blue'], $
    location=[.05,.95], linestyles=2
  
  histo = HISTOGRAM(months, BINSIZE=1, LOCATIONS=binvals)
  cgBARPLOT, histo, /NOERASE, $ ;binvals, histo, $
    POSITION=[0.70, 0.55, 0.95, 0.90], $
    XTITLE='Frequency', $
    TITLE='Events by month', $
    COLORS='GOLD', $
    BARNAMES=['J','F','M','A','M','J','J','A','S','O','N','D'], $ 
    /ROTATE
    
  ;==================================================
  ;             Event Altitude Plots
  ;==================================================
  bs=0.5
  cghistoplot, locs, BINSIZE=bs, $
    title='Event Altitudes',$
    POSITION=[0.10, 0.10, 0.45, 0.45], $
    XTITLE='Altitude(km)', $
    YTITLE='frequency', $
    POLYcolor = 'dark red', $
    /NOERASE, /FILL, /NAN
  cghistoplot, tps-locs, BINSIZE=bs, $
    MININPUT=1, $
    title='Tropopause to Event peak',$
    POSITION=[0.50, 0.10, 0.95, 0.45], $ 
    XTITLE='km from tropopause', $
    POLYcolor= 'turquoise', $
    /NOERASE, /FILL, /NAN

  ;add the title
  cgtext, .35, .95, summarytitle, $
    charsize=2, $
    /normal

  ;==================================================
  ;             Event Pressures Plots
  ;==================================================
;  !p.multi=0
;  bs=25
;  cgdisplay, 800, 500, wid=2, title='EventSummary(pressure)'
;  !P.multi=[0,2,1]
;  cghistoplot, pressures, BINSIZE=bs, $
;    title='Event Altitudes',$
;    XTITLE='Pressure(hPa)', $
;    XRANGE=[700,100], $
;    YTITLE='frequency', $
;    POLYcolor = 'dark red', $
;    /FILL, /NAN
;  cghistoplot, pressures-tps_press, BINSIZE=bs, $
;    title='(Event - TP) Pressures',$
;    XTITLE='Altitude(hPa)', $
;    POLYcolor= 'turquoise', $
;    /FILL, /NAN
  
  !p.multi=0
END



;  ; temporary loop to look at tp adjacent events:
;  cgdisplay, 500, 700, wid=5
;  toexamine=where(tps-locs lt 1)
;  print, toexamine
;  foreach i, toexamine, ei do begin
;    plotprofile, reform(ppbv[i,*]), reform(gph[i,*]), tp=tps[i], title=string(i), $
;    /altitude, top=tps[i]+5, temp=reform(data.temperature[i,*])+273.15, $
;    humid=reform(data.RH[i,*])
;    stop
;  endforeach
