;
PRO show_min_summary, melbourne=melbourne, davis=davis, macquarie=macquarie
  
  ; read in the data
  data=sondedata(melbourne=melbourne, davis=davis, macquarie=macquarie)  
  n_sondes=n_elements(data.jtime)
  sondetps=min([[data.tplr], [data.tpo3]], dimension=2,/nan)
  sondetps[where(sondetps lt 3)] = !values.f_nan
  
  events=getevents(data)
  jtimes=events.jtime
  n_events=n_elements(jtimes)
  
  print, n_sondes, ' ozonesondes'
  print, long(total(finite(data.gph[*,0]))), ' ozonesondes with gph'
  print, n_events, ' events'
  
  locs=events.locations
  gphs=events.gph / 1.0e3
  caldat, jtimes, Months, Days, Years, Hours
   
  stationstr='Melbourne'
  if keyword_set(macquarie) then stationstr='Macquarie'
  if keyword_set(davis) then stationstr='Davis'
  summarytitle=stationstr+' ('+string(n_events, format='(i2)')+' events)'

  !p.multi=0
  !p.charsize=2.25
  cgdisplay,  1000, 800 
  ;===================================================
  ;             Yearly Occurrences
  ;===================================================
  histo = HISTOGRAM(months, BINSIZE=1, LOCATIONS=binvals)
  cgBARPLOT, histo, /NOERASE, $ ;binvals, histo, $
    POSITION=[0.10, 0.6, 0.95, 0.90], $
    TITLE=summarytitle, $
    COLORS='GOLD', $
    BARNAMES=['J','F','M','A','M','J','J','A','S','O','N','D'], $ 
    XMINOR=1.0, YMINOR=1.0
    
  ;==================================================
  ;             Event Altitude Plots
  ;==================================================
  bs=0.5
  cghistoplot, locs, BINSIZE=bs, $
    title='Event Altitudes',$
    POSITION=[0.10, 0.10, 0.95, 0.45], $
    XTITLE='Altitude(km)', $
    POLYcolor = 'dark red', $
    /NOERASE, /FILL, /NAN

  !p.multi=0
END

