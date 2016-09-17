;
; non_transport_summary.pro
; Purpose
;   Summary of event occurrence, event peak altitude, and event depth
;   for a particular site
;
; HISTORY
;   Created at some point in 2016 by Jesse Greenslade
;   Updated 20160616 to make x axes the same between stations
;
PRO non_transport_summary
  
; read in the data
for i=0,2 do begin
  if i eq 0 then data=sondedata(/melbourne)
  if i eq 1 then data=sondedata(/macquarie)
  if i eq 2 then data=sondedata(/davis)
  stationstr=(['Melbourne','Macquarie','Davis'])[i]
  
  n_sondes=n_elements(data.jtime)
  sondetps=min([[data.tplr], [data.tpo3]], dimension=2,/nan)
  sondetps[where(sondetps lt 3)] = !values.f_nan
  
  events=getevents(data)
  jtimes=events.jtime
  n_events=n_elements(jtimes)
  n_trans=n_elements(events.transportinds)  

  print, n_sondes, ' ozonesondes'
  print, long(total(finite(data.gph[*,0]))), ' ozonesondes with gph'
  print, n_events, ' events'
  if n_trans eq 1 then n_trans=0 ; handle davis
  print, n_trans, ' possibly caused by transport'
  allinds=indgen(n_events)

  n_nt_events=n_events-n_trans 
  summarytitle=stationstr+' ('+string(n_nt_events, format='(i2)')+' events)'
 
  tinds=long(events.transportinds)
  if n_trans eq 0 then tinds=-1
  non_trans_inds=cgsetdifference(allinds,tinds)
  times=jtimes[non_trans_inds]
  trans_times=jtimes[tinds]

  locs=events.locations
  tps =events.tp
  gphs=events.gph / 1.0e3
  caldat, times, Months, Days, Years, Hours
  caldat, trans_times, tmonths, tdays, tyears, thours


  !p.multi=0
  !p.font=0
  !p.charsize=2.25
  cgPS_Open, stationstr+'.ps'
  cgdisplay,  1000, 1000 
  ;===================================================
  ;             Yearly Occurrences
  ;===================================================
  histo = HISTOGRAM(months, BINSIZE=1, LOCATIONS=binvals)
  thisto = HISTOGRAM(tmonths, binsize=1, nbins=12, min=1)
  cgBARPLOT, histo, /NOERASE, $ ;binvals, histo, $
    POSITION=[0.10, 0.75, 0.95, 0.95], $
    TITLE=summarytitle, $
    ytitle='Ozone Events', $
    COLORS='GOLD', $
    BARNAMES=['J','F','M','A','M','J','J','A','S','O','N','D'], $ 
    XMINOR=1.0, YMINOR=1.0
  
  ; overplot the transport events (only if not davis)
  if not keyword_set(davis) then $
    cgbarplot,  thisto, /overplot, baselines=histo, Colors='dark red'


  ;==================================================
  ;             Event Altitude Plots
  ;==================================================
  ; event altitudes and distance from tropopause
  bs=0.5
  cghistoplot, locs[non_trans_inds], BINSIZE=bs, $
    title='Event Altitudes',$
    POSITION=[0.10, 0.42, 0.95, 0.65], $
    XTITLE='Altitude(km)', $
    YTITLE='Occurences', $
    POLYcolor = 'gold', $
    mininput= 4.0, $
    xrange=[3.5,14], $
    /NOERASE, /FILL, /NAN
  ; overplot the transport events (only if not davis)
  if i ne 2 then $
    cghistoplot, locs[tinds], binsize=bs, $
      /fill, /oplot, polycolor='dark red', mininput=4.0, $
      /line_fill, orientation=45, line_thick=2.0
  
  ;==================================================
  ;             Event Depth Plots
  ;==================================================
  ; event peaks distance from tropopauses
  cghistoplot, (tps-locs)[non_trans_inds], BINSIZE=bs, $
    title='Event Depths',$
    POSITION=[0.10, 0.10, 0.95, 0.30], $ 
    XTITLE='km below tropopause', $
    POLYcolor= 'gold', $
    XRANGE=[0,9], $
    mininput=0, $
    /NOERASE, /FILL, /NAN, /OPROBABILITY
  ; overplot the transport events (only if not davis)
  if not keyword_set(davis) then $
    cghistoplot, (tps-locs)[tinds], binsize=bs, $
      /fill, /oplot, polycolor='dark red', $
      /line_fill, orientation=45, line_thick=2, mininput=0
  ; save the output as .ps
  CGPS_Close
  ;Output=stationstr+'.ps'
endfor 
END

