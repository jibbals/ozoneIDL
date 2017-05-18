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
    co_arr=[0.99, 0.99, 0.985, 0.98, 0.97, 0.96, 0.95]
    trop_def_arr=intarr(n_elements(co_arr))
    trop_def_arr[0]=2 ; only use mintp for 'orig'
    name_arr=['_orig', '_tpo3','_co985','_co98','_co97','_co96','_co95']
    for ii=0,n_elements(co_arr)-1 do begin
      trop_def=trop_def_arr[ii]
      cutoff=co_arr[ii]
      fname='Data/'+sitestr[i]+name_arr[ii]+'.csv'
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
      filedata=transpose([[Years], [Months], [Days], [Hours], [tropozone], [flux], $
        [event_heights], [event_tps], [fire_influence]])
      header=['YYYY','MM','DD','HH','tropozone','flux','peak','tp','fire']
      write_csv, fname, filedata, header=header
    endfor
  endfor

end
   
  
