;
; Put all event dates into .csv files
PRO event_dates_to_csv

  melb=ptr_new(sondedata(melbourne=1))
  mac =ptr_new(sondedata(macquarie=1))
  dav =ptr_new(sondedata(davis=1))
  ;dataptr = ptrarr(3)
  dataptr = [melb, mac, dav]
  sitestr = ['melb','mac','dav']
  ;header=['Year','Month', 'Day', 'Hour']
  for i=0,2 do begin
    data=*(dataptr[i])
    for ii=0,3 do begin
      trop_def=([2,0,0,0])[ii]
      cutoff=([0.99,0.99,0.985,0.98])[ii]
      fname=sitestr[i]+(['_orig', '_tpo3','_tpo3_co985','_tpo3_co980'])[ii]+'.csv'
      ; for this station, pull out all the dates
        events=getevents(data,trop_def=trop_def,cutoff_percentile=cutoff)
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
      write_csv, fname, data, header=header
    endfor
  endfor

end
   
  
