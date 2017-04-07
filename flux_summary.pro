; PURPOSE:
;   Plot small summary of relative and absolute ozone flux at three sonde sites
;
; UPDATES:
;   Now just look at non transport events ~2016/Feb
;   Box and whisker plot instead of error bars 2016/07/26
;
PRO flux_summary
  
  ; read in the data
  davi = sondedata(/davis)
  macq = sondedata(/macquarie)
  melb = sondedata(/melbourne)
  ptrlist = [ptr_new(davi), ptr_new(macq), ptr_new(melb)]
  sitelist = ['Davis','Macquarie','Melbourne']
  fluxes = ptrarr(3)
  relfluxes= ptrarr(3)
  n_events = intarr(3)
  n_profiles = intarr(3)
  relbox=dblarr(5,3)
  absbox=dblarr(5,3)
  FOREACH DATAptr, ptrlist, ii DO BEGIN
    DATA = *DATAptr
    n_sondes=n_elements(DATA.jtime)
    sondetps=min([[DATA.tplr], [DATA.tpo3]], dimension=2,/nan)
    sondetps[where(sondetps lt 3)] = !values.f_nan
  
    ; EventFlux is stored in events returned by getevents(DATA)
    events=getevents(DATA)
    jtimes=events.jtime
    dates=jessedate(jtimes)
    ppbvs=events.o3ppbv
    ev_count=n_elements(jtimes)
    flux=events.flux    ;molecules/cm2
    tropozone=events.tropozone
    relflux=flux/tropozone * 100.
    
    n_gph = long(total(finite(DATA.gph[*,0])))
    ;print, n_sondes, ' ozonesondes'
  
    locs=events.locations
    gphs=events.gph / 1.0e3
    tps =events.tp
    caldat, jtimes, Months, Days, Years, Hours
    
    ; UPDATED: also store data with transport removed
    nontpinds = cgsetdifference(indgen(ev_count), long(events.transportinds))
      
    n_events[ii] = n_elements(nontpinds)
    
    fluxes[ii]=ptr_new(dblarr(n_events[ii]))
    (*fluxes[ii])[*]=flux[nontpinds]
    relfluxes[ii]=ptr_new(dblarr(n_events[ii]))
    (*relfluxes[ii])[*]=relflux[nontpinds]
    
    print, 'Tropospheric ozone in molecules/cm3 at ', sitelist[ii]
    print, 'Mean ozone on event days: ', mean(tropozone)
    print, 'Max ozone on event days : ', max(tropozone)
    print, 'Mean event flux: ', mean(flux)
    print, 'Max event flux : ', max(flux)
    print, 'Mean relative flux: ', mean(relflux)
    print, 'Max relative flux : ', max(relflux)
    print, n_gph, ' ozonesondes with gph'  
    print, ev_count, ' events'
    print, 'relative mean sans transport: ', mean(*relfluxes[ii])
    print, 'absolute mean sans transport: ', mean(*fluxes[ii])
    
    ; UPDATE box plots made from saving (min, lower quartile, median, upper q, max)
    relbox[*,ii] = createboxplotdata(relflux[nontpinds])
    absbox[*,ii] = createboxplotdata(flux[nontpinds])
  ENDFOREACH
  ; BoxPlots:
  cgLoadct, 33
  width=0.8
  colr=cgcolor("burlywood")
  xticknames=['','Davis','Macquarie','Melbourne','']
  cgps_open, "flux_relative.ps"
  cgPlot, [0.,15.], /nodata, ytitle='%', YStyle=8, $
    title="Portion of ozone from STTs", xtitle="", XStyle=8, xrange=[-0.5, 2.5], $
    xtickname=xticknames, xtickv=[0,1,2]
  for i=0,2 do $
    cgboxplot, *relfluxes[i], xlocation=i, boxcolor=colr, $
      thick=2, /overplot, width=width ;, /fillbox ; can't see mean
  cgps_close
  
  cgps_open, "flux_absolute.ps"
  cgPlot, [min(absbox),max(absbox)], /nodata, $
    ytitle='Ozone [molec/cm2]', YStyle=8, xrange=[-0.5,2.5], xtickv=[0,1,2],$
    title="STT Ozone flux", xtitle="", XStyle=8, xtickname=xticknames
  for i=0,2 do $
    cgboxplot, *fluxes[i], xlocation=i, boxcolor=colr, $
      thick=2, /overplot, width=width ;, /fillbox ; can't see mean
  cgps_close

  stop
END
