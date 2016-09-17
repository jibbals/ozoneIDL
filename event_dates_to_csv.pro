;
; Put all event dates into .csv files
PRO event_dates_to_csv

  melb=ptr_new(sondedata(melbourne=1))
  mac =ptr_new(sondedata(macquarie=1))
  dav =ptr_new(sondedata(davis=1))
  ;dataptr = ptrarr(3)
  dataptr = [melb, mac, dav]
  filenames = ['melb.csv','mac.csv','dav.csv']
  ;header=['Year','Month', 'Day', 'Hour']
  for i=0,2 do begin
    data=*(dataptr[i])
    ; for this station, pull out all the dates
    events=getevents(data)
    jtimes=events.jtime
    caldat, jtimes, Months, Days, Years, Hours
    
    ; also pull out the flux amount and percentage
    flux=events.flux
    tropozone=events.tropozone
    ; flag the fire transport events
    fire_influence=strarr(n_elements(jtimes))
    fire_influence[*]=0
    if TOTAL(finite(events.transportinds)) ne 0 then $
      fire_influence[events.transportinds] = 1
    event_heights=events.locations
    event_tps=events.tp
    data=transpose([[Years], [Months], [Days], [Hours], [tropozone], [flux], $
      [event_heights], [event_tps], [fire_influence]])
    header=['YYYY','MM','DD','HH','tropozone','flux','peak','tp','fire']
    write_csv, filenames[i], data, header=header
  endfor

end
   
  
