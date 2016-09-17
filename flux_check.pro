;
PRO flux_check, date, davis=davis, melbourne=melbourne, macquarie=macquarie, cgps=cgps
    
    ; read in the data
    site = sondedata(davis=davis, melbourne=melbourne, macquarie=macquarie)
  
    ; EventFlux is stored in events returned by getevents(DATA)
    events=getevents(site)
    jtimes=events.jtime
    ppbvs=events.o3ppbv
    ev_count=n_elements(jtimes)
    flux=events.flux
    tropozone=events.tropozone
    relflux=flux/tropozone * 100.
    ;print, n_sondes, ' ozonesondes'
  
    locs=events.locations
    gphs=events.gph / 1.0e3
    tps =events.tp
    tmp=min(abs(jtimes-date), mindex)
    caldat, jtimes, Months, Days, Years, Hours
    print,'ymd:', Years[mindex],Months[mindex],Days[mindex]
    print,tmp
    print,'flux: ',flux[mindex]
    print,'relflux: ',relflux[mindex]

    show_profile, jday=jtimes[mindex], davis=davis, melbourne=melbourne, macquarie=macquarie, cgps=cgps
    stop
END
