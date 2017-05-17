
;   NAME:
;       getevent
;       
;   Purpose:
;       loop through columns and return 
;       columns which are classed as events:
;       
;   Inputs:
;       sondes: structure returned from sondedata() function
;   
;   Outputs: events
;       events: structure{Profiles, Locations, Magnitudes, indices, type}
;       Profiles are event profiles, locations are height of maximum ozone perturbance, 
;   Example:
;       getevents(sondedata(/melbourne), /example)
;           this will return melbourne events, as well as save an example picture
;           of the filtering process
;

FUNCTION regrid_profile, T, F, delta, top, $
  TGrid=TGrid
;   Purpose:
;     Regrid profile from 2km up to tropopause with interpolated regular grid spacing
;   Inputs:
;     T: altitude of profile
;     F: ppbv of profile
;   Returns:
;     Regridded ppbv
;   Outputs:
;     TGrid: regridded altitude
;
  ; set up standard uniform grid:
  ;
  gridmax=top
  gridlen=gridmax/delta
  Tgrid=findgen(gridlen)*delta 
  FGrid=interpol(F,T,TGrid,/nan)
  return, FGrid
END

PRO show_profile_filter, ppbv, filt, gph, evtp, cutoff, $
    name=name, $
    title=title, peak=peak
;
;  Show the profile and filter passed in against the gph
;  set name to save out using cgps_open, name and cgps_close
    
    if n_elements(name) ne 0 then begin
        cgps_open, name
        !p.font=0
    endif
    ; show the profile
    !p.multi=[0,2,1]
    !p.charsize=1
    !X.OMargin = [2, 6]
    !Y.OMargin = [2, 6]
    cgdisplay, 900, 900, wid=6
    xlims=[20,100]
    ylims=[2,15]
    ytickvalues=[5,7,9,10,11,12,14]
    cgplot, ppbv, gph, title='Profile', /ylog,$
        xrange=xlims, yrange=ylims, ytickv=ytickvalues, $
        ytitle='Altitude (km)', xtitle='Ozone (ppbv)'
    c0=cgcolor('purple')
    c1=cgcolor('orange')
    ; line at tropopause: dashed, vertical, red
    cgoplot, [xlims[0], xlims[-1]], [evtp, evtp], $
        linestyle=2, color=c0
    
    ; if peak set then plot an x on the peak
    if n_elements(peak) gt 0 then begin
      peakind=where(ppbv eq peak)
      if n_elements(peakind) gt 1 then peakind=peakind[0]
      cgoplot, peak, gph[peakind], psym=1
    endif
    
    ; Add plot of transformed profile, against zeroline
    cgplot, filt, gph, /ylog, $
      title='Transformed', ytitle='Altitude (km)', $
      yrange=ylims, xtitle='perturbation (ppbv)', ytickv=ytickvalues
    cgoplot, [0,0], ylims, color=c0, linestyle=2
    cgoplot, [cutoff, cutoff], ylims, color=c1, linestyle=2
    if n_elements(title) gt 0 then begin
      cgText, 0.5, 0.95, ALIGNMENT=0.5, CHARSIZE=1.25, /NORMAL, $
        title
    endif
    if n_elements(name) ne 0 then cgps_close
    !p.multi=0
END

;
; THE FUNCTION:
;
function getevents, sondes, verbose=verbose, CUTOFF_PERCENTILE=CUTOFF_PERCENTILE, $
  analyse=analyse, example=example, show_storm=show_storm, PLOT_GRADIENTS=PLOT_GRADIENTS, $
    FOURIER_SCALE=FOURIER_SCALE

  ; SET SOME THINGS HERE
  ;
  ; Lower boundary to clip off of our filtered profile
  LOWERCLIP=2.0
  ; LOWER EVENT BOUNDARY  
  EVENT_BOUND_LOW=4.0
  ; Upper boundary to clip off of our filtered profile
  UPPERCLIP=0.5
  ; Upper boundary for determining cutoff region
  UPPERCLIP_CUTOFF=1.0
  ; Percentile to use determining cutoff
  IF N_ELEMENTS(CUTOFF_PERCENTILE) EQ 0 THEN CUTOFF_PERCENTILE=0.99
  IF N_ELEMENTS(FOURIER_SCALE) EQ 0 THEN FOURIER_SCALE=[0.5,5.0]
  ; Drop required before we say the event is seperate from the stratosphere
  GRADIENT_DROP=20
  ; value of ozone ppb in the 3km above the event peak must drop below this
  GRADIENT_THRESH=80

  if float(sondes.lat) lt -37 and float(sondes.lat) gt -38 then sitename='Melbourne'
  if float(sondes.lat) lt -54 and float(sondes.lat) gt -55 then sitename='Macquarie'
  if float(sondes.lat) lt -68 and float(sondes.lat) gt -69 then sitename='Davis'
  imglocation='./images/'+sitename+'/'
  
  ppbv=sondes.o3ppbv
  vmrs=1.0d-9*ppbv
  alt=sondes.gph/1.0e3  ; m to km
  alts=alt ; THIS IS kms duuuude
  alts[where(alts eq 0)] = !values.f_nan
  kels = sondes.temperature + 273.15 ; Celcius to Kelvins
  hpas = sondes.pressure; pressure profile in hPa
  tpo3 = sondes.tpo3    ; ozone tropopause
  tplr = sondes.tplr    ; lapse rate tropopause
  jtimes = sondes.jtime ; julian time axis
  years = jessedate(jtimes)   ; julian time array to array of 'years from start of array'
  
  ; minimum tropopause = 0 if either tp is NAN
  mintp=min( [[tplr],[tpo3]], dimension=2, /nan )
  if n_elements(mintp) lt n_elements(tpo3) then stop
  
  N_profs=n_elements(ppbv[*,0]) ; how many profiles
  N_alts=n_elements(ppbv[0,*])  ; how many altitudes
  
  ; Make arrays to hold redefined ozone troposphere profiles
  ;
  delta=0.02 ; 20m new profile resolution
  top=14.0   ; up to 14km
  regridlen=top/delta
  TGrid=findgen(regridlen)*delta ; ALSO IN kms
  troposphere = make_array(N_profs, regridlen, /FLOAT, VALUE=!values.f_nan )
  SFTo3 = make_array(N_profs, regridlen, /FLOAT, VALUE=!values.f_nan )
  SFTo3_for_percentiles=make_array(N_profs, regridlen, /FLOAT, VALUE=!values.f_nan )
  SFTo3_for_detection=make_array(N_profs, regridlen, /FLOAT, VALUE=!values.f_nan )
  locs = make_array(N_profs, /FLOAT, VALUE=!values.f_nan)
  flux = make_array(N_profs, /FLOAT, VALUE=!values.f_nan)
  tropozone=make_array(N_profs, /FLOAT, VALUE=!values.f_nan)
  ; Make arrays to hold ozone bubble upper and lower bounds
  ubs=fltarr(N_profs)+!values.f_nan
  ubinds=fltarr(N_profs)+!values.f_nan
  lbinds=fltarr(N_profs)+!values.f_nan
  lbs=fltarr(N_profs)+!values.f_nan
  ; flags for edge case removals
  tpflags=intarr(n_profs)
  gradflags=intarr(n_profs)
  lbflags=intarr(n_profs)
  
  ; work out density columns!
  ; 
  density = make_array(n_profs, regridlen, /double, value=!values.d_nan)
  ; need regridded Pressure and Temperature
  reP = make_array(n_profs, regridlen, /float, value=!values.f_nan)
  reT = make_array(n_profs, regridlen, /float, value=!values.f_nan)
  reVMR = make_array(n_profs, regridlen, /float, value=!values.f_nan)
  for i=0,n_profs-1 do begin
    reP[i,*]   = interpol(reform(hpas[i,*]),reform(alts[i,*]),TGrid, /nan) ; hPa
    reT[i,*]   = interpol(reform(kels[i,*]),reform(alts[i,*]),TGrid, /nan) ; Kelvins
    reVMR[i,*] = interpol(reform(vmrs[i,*]),reform(alts[i,*]),TGrid, /nan)
  endfor
  
  ; gas constant R
  R = 8.3144621 ; [cm3 MPa k-1 mol-1]
  R = R * 1.0d4 ; [cm3 hPa k-1 mol-1]
  ; number density(n_i / V): [ mol_{O_3} / cm^3 ] 
  ; from P_i = n_i R T / V, and P_i = P vmr
  ;   n_i / V = P vmr / ( R T ) = mol_{O_3} / cm^3 
  ;  mol/cm3 * Avegadro's = molecules/cm3
  density = reVMR * reP / R / reT * 6.02214129e23 ; in molecules/cm3
  
  low=FOURIER_SCALE[0]/delta           ;km freq min
  high=FOURIER_SCALE[1]/delta           ;km frequency max
  
  ; for each profile
  for i=0,N_profs-1 do begin
    tpbound=mintp[i]
    if not finite(tpbound) then continue  ; ignore columns with no tp
    if tpbound lt 2 then continue ; ignore columns with tp of 2 or less
    ; Regrid the first 14km
    ;
    troposphere[i,*] = regrid_profile(reform(alt[i,*]),reform(ppbv[i,*]), delta, top)
 
    ; remove ppbv above the tropopause before filtering
    ;
    tpup=where(TGrid ge tpbound)
    troposphere[i,tpup] = !values.f_nan
    
    ; Here we skip profiles which start with NAN value
    inds=where(finite(troposphere[i,*]))
    if inds[0] eq -1 then continue
    
    ; run the bandpass filter on the column up to the tropopause
    ;
    sa_filter_1d, reform(troposphere[i,inds]), low, high, filtered
    SFTo3[i,inds] = filtered  ; Simon's Fourier Transformed O3
    SFTo3_for_percentiles[i,inds] = filtered  ; three copies
    SFTo3_for_detection[i,inds] = filtered

    ; set range on basic filter, used in plots and event boundary determination
    clip=where(TGrid ge tpbound-UPPERCLIP or TGrid lt LOWERCLIP)
    SFTo3[i,clip]=!values.f_nan
    
    ; set range on filter for event detection
    clip_det=where(TGrid ge tpbound-UPPERCLIP or TGrid lt EVENT_BOUND_LOW)
    SFTo3_for_detection[i,clip_det]=!values.f_nan
    
    ; set range on filter for percentiles
    clip_perc=where(TGrid ge tpbound-UPPERCLIP_CUTOFF or TGrid lt LOWERCLIP)
    SFTo3_for_percentiles[i,clip_perc] = !values.f_nan

    ; location of ozone peak in profile
    ; OLD LOCATION was set using maximum within filtered profile
    tmp=max(SFTo3_for_detection[i,*],maxind,/nan)
    locs[i]=tgrid[maxind]
    ; UPDATE, LOWEST THRESH EXCEEDENCE (only done for events)
  endfor
  
  ; Percentiles of filtered ozone profile datapoints
  pcs=[CUTOFF_PERCENTILE]
  percentiles=cgpercentiles(SFTo3_for_percentiles[where(finite(SFTo3_for_percentiles))], percentiles=pcs)
  cutoff=percentiles[-1]  ; we use the percentile set by cutoff_percentile variable.
  
  ; remove values below the cutoff
  testarr=SFTo3_for_detection
  testarr[where(testarr lt cutoff)] = !values.f_nan
  testprof=total(testarr,2,/NAN)
  eventinds=where(testprof gt 0.0)
  n_events=n_elements(eventinds)
  if keyword_set(verbose) then $
    print, n_events,' events detected'

  ; DETERMINE FLUX CAPACITOR
  ; 
  ; for each event
  tpcount=0 & gradcount=0
  foreach ind, eventinds, ii do begin
    evtrop=reform(troposphere[ind,*])
    evfilt=reform(SFTo3[ind,*])
    evdens=reform(density[ind,*])
    evlocold=locs[ind]
    evlocindold=where(TGrid eq evlocold)
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; we determine event peak using max ozone within lowest threshhold exceedence range
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    peak_filtered=reform(SFTo3_for_detection[ind,*])
    over_thresh= where(peak_filtered ge cutoff)
    below_thresh=where(peak_filtered lt cutoff)
    peak_index_low=over_thresh[0]
    first_below_thresh= (where(below_thresh gt peak_index_low))[0]
    peak_index_high=below_thresh[first_below_thresh]
    ; handle when the profile does not drop back below the threshold
    if first_below_thresh eq -1 then peak_index_high=over_thresh[-1]
    if peak_index_high lt peak_index_low then begin
      print, "trouble determining peak:" 
      stop
    endif
    peak_o3=max(evtrop[peak_index_low:peak_index_high])
    evlocind=(where(evtrop eq peak_o3))[0]
    evloc=(TGrid[evlocind])
    locs[ind]=evloc
    ; lower bound is first zero below 'event' ozone values in the filtered profile
    lbinds[ind] = max(where(evfilt le 0 and tgrid lt evloc))
    if lbinds[ind] eq -1 then begin
      lbinds[ind] = where(tgrid eq 2.0)
      print, 'event ', ind, ' at ', evloc 
      caldat, jtimes[ind], tempm, tempd, tempy
      print, 'date: ', tempy, tempm, tempd
      lbflags[ind]=1
      ;stop
    endif
    lbs[ind] = tgrid[lbinds[ind]]
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;; Gradients flag set if pertubed does not pass zero ;;;;;;
    ;;;;;;; And (peak-min <20 OR min < 80)                    ;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ; upper bound is first zero above 'event' or else minimum above event location
    ubinds[ind] = min(where(evfilt le 0 and tgrid gt evloc))
    ; if pertubed does not drop to zero between event and TP:
    if ubinds[ind] eq -1 then begin
      minppbv=min(evtrop[evlocind:*],minind,/nan)
      ubinds[ind]=minind+evlocind
      tpflags[ind]=1
      tpcount=tpcount+1
      ;;;; UPDATE: separation is OK if peak to min distance is gt thresh1 and min < thresh2
      if (evtrop[evlocind]-minppbv lt GRADIENT_DROP) OR (minppbv gt GRADIENT_THRESH) then begin 
        gradflags[ind]=1
        gradcount = gradcount+1
        ; plot all these so I can see what's going down in Sydney town
        if KEYWORD_SET(PLOT_GRADIENTS) gt 0 then begin
          yaxi = where(tgrid ge 2 and tgrid le 14)
          yax  = tgrid[yaxi]
          evtime = jtimes[ind]
          evtp= mintp[ind]
          evFT = reform(SFTo3[ind,*])
          caldat, evtime, mm,dd,yy
          datestr=string(yy,format='(i04)')+string(mm,format='(i02)')+string(dd,format='(i02)')
          title="Filter on "+ datestr+ ' flag = '+string(gradflags[ind], '(I2)')
          name=imglocation+'filters/flagged'+datestr+'.ps'
          show_profile_filter, evtrop[yaxi], evfilt[yaxi], yax, evtp, cutoff, title=title, $
            name=name, peak=locs[ind]
        endif
      endif
    endif
    ubs[ind] = tgrid[ubinds[ind]]
    
    ;;;;;;;;
    ;;;; OLD WAY
    ;;;;;;;
    dx=TGrid[ubinds[ind]]-TGrid[evlocind] ; change in altitude
    dy=evtrop[ubinds[ind]]-evtrop[evlocind] ; change in ppbv
    grad = dy/dx
    if dx lt .1 then grad = 0

    ; These should already be caught by other method... otherwise plot them
    if (grad gt -20) and (gradflags[ind] ne 1) then begin
    
      if KEYWORD_SET(PLOT_GRADIENTS) gt 0 then begin
        yaxi = where(tgrid ge 2 and tgrid le 14)
        yax  = tgrid[yaxi]
        evtime = jtimes[ind]
        evtp= mintp[ind]
        evFT = reform(SFTo3[ind,*])
        
        caldat, evtime, mm,dd,yy
        datestr=string(yy,format='(i04)')+string(mm,format='(i02)')+string(dd,format='(i02)')
        title="Filter on "+ datestr+ '\n gradflag = '+string(gradflags[ind], '(I2)')
        name=imglocation+'filters/oldway_flagged'+datestr+'.ps'
        show_profile_filter, evtrop[yaxi], evfilt[yaxi], yax, evtp, cutoff, title=title, $
          name=name, peak=locs[ind]
      endif
    endif

    if keyword_set(show_storm) then begin 
      if jtimes[ind] gt 2453403 and jtimes[ind] lt 2453407 then begin
        print,'storm gradient=', grad
        print, 'using dx= ',dx,' dy=',dy
        print, 'from ub alt= ',ubs[ind], ' and ub o3 = ', evtrop[ubinds[ind]], ' and peak o3= ',evtrop[evlocind]
        print, 'event ', ind, ' at ', evloc 
        print, 'gradflag=', gradflags[ind]
      endif
    endif
    ; line between where filtered zeros occur above and below the ozone peak(using density profile)
    baseline=interpol(evdens[[lbinds[ind],ubinds[ind]]],[lbs[ind],ubs[ind]], tgrid[lbinds[ind]:ubinds[ind]] )
    area=evdens[lbinds[ind]:ubinds[ind]] - baseline
    ; conservative estimate for event flux is integral over baseline 
    ; of 'event'
    ; Do this using density column!
    ;
    flux[ind]=int_tabulated(tgrid[lbinds[ind]:ubinds[ind]],area)
    ; integrate density column UP TO TROPOPAUSE 
    evtpind=max(where(tgrid le mintp[ind]))
    tropozone[ind] = int_tabulated(tgrid[0:evtpind], evdens[0:evtpind])
    ; since we integrate over the vertical (which is in kms) this is now in
    ; molecules kilometres / cm3  so we just multiply by cms/kilometre
    tropozone[ind] = tropozone[ind] * 1e5
    flux[ind] = flux[ind] * 1e5
  endforeach

  if keyword_set(verbose) then $
    print,'jtimes with gradient and tp flags=',n_elements(jtimes[where(gradflags gt 0 and tpflags gt 0)])
  
  ; remove gradient shallow flagged events
  ;
  removes=where(gradflags eq 1 and tpflags eq 1)
  finalevents = cgsetdifference(eventinds, removes)
  
  ; plot all the events and their filters
  if keyword_set(PLOT_GRADIENTS) then begin
    foreach evind, finalevents, evindi do begin
      yaxi = where(tgrid ge 2 and tgrid le 14)
      yax  = tgrid[yaxi]
      evtime = jtimes[evind]
      evppbv = reform(troposphere[evind, *])
      evtp= mintp[evind]
      evFT = reform(SFTo3[evind,*])
      
      caldat, evtime, mm,dd,yy
      datestr=string(yy,format='(i04)')+string(mm,format='(i02)')+string(dd,format='(i02)')
      title=sitename+' '+datestr
      name=imglocation+'filters/EVENTS/'+datestr+'.ps'
      ; show the profile
      show_profile_filter, evppbv[yaxi], evFT[yaxi], yax, evtp, cutoff, title=title, $
        name=name, peak=locs[evind]
    endforeach
  endif
  
  ; Determine event types based on distance from tropopause
  ; 0-2km = shallow, 2-4 = medium, 4+ = deep
  ;
  distance=mintp-locs
  bins=[2,4]  ; 0 to 2 is shallow, 2 to 4 is mediumn, 4+ is deep event
  type = (distance[finalevents] gt bins[0]) + (distance[finalevents] ge bins[1])
  
  ; List of Events probably caused by BB transport
  ;
  
  ; Melbourne Transport Events
  if strmatch(sitename,'Melbourne') then begin 
    ; 20041007 20051026 20051103 20051108 20051117 20061016 (Singapore fires)
    ;   20071017 20071031 20090916 20100902 20101027 20101103 20111025 20121003
    ;   20121030 20131004
    ; julday(9,2,2010) removed as below 4km
    transportDays=[ julday(8,11,2004), julday(10,26,2005), $
      julday(11,3,2005), $
      julday(11,8,2005),julday(11,17,2005),julday(10,16,2006),julday(10,17,2007),$
      julday(10,31,2007),julday(9,16,2009),julday(10,27,2010),$
      julday(10,25,2011),julday(10,3,2012),julday(10,30,2012),$
      julday(10,4,2013)] - 0.5
    n_transport_excluded=n_elements(transportDays)
    ; subtract 0.5 to use midnight of julian days instead of midday, this is
    ;  so the UTC+0 event times are more closely matched

    ; get indices of transport events by closest date match
    transportInds=make_array(n_transport_excluded)
    ii=0
    foreach tday, transportDays do begin
      tmp=min(abs(jtimes[finalevents]-tday), tind)
      ; if more than a day between transport event and matched event then stop
      if tmp gt 1 then begin
        print, "bad transport event date: ", string(tday, format='(F10.1)'), tmp
        caldat, tday, mm,dd,yy
        print, string(yy,format='(i04)')+string(mm,format='(i02)')+string(dd,format='(i02)')
      endif else begin
        transportInds[ii]=tind
        ii=ii+1
      endelse
    endforeach

    if ii ne n_transport_excluded then begin
      print, "missing transport event?"
    endif
  ; Macquarie Transport Events
  endif else if strmatch(sitename,'Macquarie') then begin
    ;  20041020 20050825 20051019 20061012 20070117 20081022 20110921 20120926
    ;
    transportDays=[  julday(10,20,2004),julday(8,25,2005),julday(10,12,2005), $
        julday(10,12,2006), julday(1,17,2007), $
        julday(10,22,2008), julday(9,21,2011), julday(9,26,2012) ] - 0.5
    n_transport_excluded=n_elements(transportDays)
    transportInds=make_array(n_transport_excluded)
    ii=0
    foreach tday,transportDays do begin
      tmp=min(abs(jtimes[finalevents]-tday), tind)
      if tmp gt 1 then begin
        print, "bad transport event date: ", string(tday, format='(F10.1)'), tmp
        caldat, tday, mm,dd,yy
        print, string(yy,format='(i04)')+string(mm,format='(i02)')+string(dd,format='(i02)')
      endif else begin
        transportInds[ii]=tind
        ii=ii+1
      endelse
    endforeach

    if ii ne n_transport_excluded then begin
      print, "missing transport event?"
    endif
  ; not melb or mac
  endif else transportInds=!values.f_nan

  events={jtime:jtimes[finalevents], $    ; julian time stamps
    o3ppbv:ppbv[finalevents, *], $          ; ozone ppbv
    gph:sondes.gph[finalevents, *], $       ; geopotential height
    temperature:sondes.temperature[finalevents,*], $
    pressure:sondes.pressure[finalevents,*], $
    tpo3:sondes.tpo3[finalevents], $        ; ozone defined tp
    tplr:sondes.tplr[finalevents], $        ; lapse rate tp
    rh:sondes.rh[finalevents, *], $         ; relative humidity
    lat:sondes.lat,$                        ;
    lon:sondes.lon,$                        ;
    height:sondes.height,$                  ; altitude in kms
    tp:mintp[finalevents], $                ; minimum tropopause
    locations:locs[finalevents], $          ; peak ozone height in km
    density:density[finalevents,*], $       ; column ozone in molecs/cm3
    tropozone:tropozone[finalevents], $     ; tropospheric column ozone in molecs/cm2
    flux:flux[finalevents], $               ; event estimated ozone flux in molecs/cm2
    tpflag:tpflags[finalevents], $          ; flag for upper bound of event hitting tropopause-1
    gradflag:gradflags[finalevents], $      ; flag for gradiend above event not falling fast enough
    lbflag:lbflags[finalevents], $          ; flag for lower bound of event hitting 2km limit
    indices:finalevents, $                  ; indices of events for comparison with initial dataset
    eventfilters:SFTo3[finalevents,*], $    ; filter used determining events 
    type:type, $                            ; class based on ozone peak distance to tp
    transportInds:transportInds }           ; events caused by BB transport

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;; EXAMPLE PLOT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if keyword_set(example) then begin
    ; Put out an example of the bandfilter
    ;
    evind = finalevents[0]
    yaxi = where(tgrid ge 2 and tgrid le 14)
    yax  = tgrid[yaxi]
    evtime = jtimes[evind]
    evppbv = reform(troposphere[evind, *])
    ; updatE: use event density rather than ppbv to show flux
    ; update2: use ppbv for both plots
    evdens = reform(density[evind,*])
    evP = reform(sondes.pressure[evind,*])
    evtp= mintp[evind]
    evFT = reform(SFTo3[evind,*])
    
    caldat, evtime, mm,dd,yy
    title1="Ozone at Melbourne on "+string(yy,format='(i04)')+'/'+ $
        string(mm,format='(i02)')+'/'+string(dd,format='(i02)')

    ; show the profile
    name='images/filtereg.ps'
    cgps_open, name
    !p.multi=[0,2,1]
    !p.charsize=2
    !p.font=0
    cgdisplay, 1200, 1200, wid=6
    xlims=[20,100]
    ;xlims=[min(evdens[yaxi],/nan), max(evdens[yaxi],/nan)]
    ylims=[1,10]
    cgplot, evppbv[yaxi], yax, title=title1, $
        xrange=xlims, yrange=ylims, $
        ytitle='Altitude (km)', xtitle='Ozone (ppbv)'
    c0=cgcolor('purple')
    c1=cgcolor('orange')
    ; line at tropopause: dashed, horizontal, red
    cgoplot, [xlims[0], xlims[-1]], [evtp, evtp], $
        linestyle=2, color=c0

    ; flux calculations
    evlbind = lbinds[evind]
    evubind = ubinds[evind]
    evlb = lbs[evind]
    evub = ubs[evind]
    evex = tgrid[evlbind:evubind]
    evbaseline=interpol(evppbv[[evlbind,evubind]],[evlb,evub], evex )
    evtopline =evppbv[evlbind:evubind]
    
    ; overlay flux outline
    cgoplot, evtopline, tgrid[evlbind:evubind], linestyle=2, color=c1, thick=2
    cgoplot, evbaseline, tgrid[evlbind:evubind], linestyle=2, color=c1, thick=2
    
    ; Add plot of transformed profile, against zeroline
    cgplot, evFT, tgrid, $
      title='Fourier transformed profile', ytitle='Altitude (km)', $
      yrange=ylims, xtitle='perturbation (ppbv)'
    cgoplot, [0,0], ylims, color=c0, linestyle=2
    cgoplot, [cutoff, cutoff], ylims, color=c1, linestyle=2
    ; line at tropopause: dashed, horizontal, red
    cgoplot, [-30, 30], [evtp, evtp], $
        linestyle=2, color=c0
    cgps_close
    !p.multi=0
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;; ANALYSIS PLOTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if keyword_set(analyse) then begin
    ; printables
    ;
    stationstr='[' + sondes.lon + 'E, ' + sondes.lat + 'N]'
    print, 'station: ', stationstr
    print, 'percentiles: ', [percentiles]
    print, 'Using ',cutoff,' as cutoff'
    print, N_profs, ' Sonde profiles'
    print, n_elements(eventinds), ' Events found.'
    print, tpcount, ' Events occured with upper bound set to the tropopause - 1 km'
    print, gradcount, ' Events have descent gradients greater than -20ppbv/km'
    ; plotables
    ;
    !p.charsize=1.0
    titles='averaged '+['', 'filtered '] + 'troposphere ppbv ' + stationstr 
    cgdisplay, 1100, 700, wid=3
    ; Plot avg+- 1std for both the tropospheric ozone and the filtered perturbations
    ;
    for i=0,1 do begin
      case i of
        0:array=troposphere
        1:array=SFTo3
      endcase 
      moments = moment(array, dimension=1, /nan)
      mn=moments[0]
      std=sqrt(moments[1])
      cgplot, TGrid, mn, title=titles[i], xtitle='altitude(km)', $
        position=[.1,.55,.65,.95] - i*[0,.5,0,.5], $
        NOERASE=i
      cgoplot, TGrid, mn+std, color='grey'
      cgoplot, TGrid, mn-std, color='grey'
      ; store the mean tropopause for next graphic
      if i eq 0 then meantropo=mn
    endfor
    cglegend, length=0.0, tcolor='grey', titles='+- 1 std dev', location=[.1,.89]
    cglegend, length=0.02, tcolor='black', titles='cutoff', location=[.84,.89]
    
    ; plot distribution of tropospheric points
    ;
    pdf=histogram(SFTo3, nbins=20,locations=ybins, /NAN)
    bs=ybins[1]-ybins[0]
    cghistoplot, SFTo3, binsize=bs, $
      /ROTATE, /FILL, /NAN, /NOERASE, $
      title='Filtered distribution', $
      position=[.70,.05,.96,.95] 
    cgoplot, !X.crange, fltarr(2)+cutoff
  
    ; plot average profiles, split by type
    ;
    evtrop=troposphere[eventinds,*]
    titles='Average '+['shallow','medium','deep'] + ' event profiles'
    cgdisplay, 900, 900, wid=4
    !p.multi=[0,1,3]
    for i=0, 2 do begin
      class=where(type eq i)
      titles[i] = titles[i] + '(' + string(n_elements(class),format='(i3)') + ' events) ' +stationstr 
      tropi=evtrop[class,*]
      moments = moment(tropi, dimension=1, /nan)
      mn=moments[0]
      std=sqrt(moments[1])
      cgplot, TGrid, mn, title=titles[i], xtitle='altitude(km)', ytitle='Ozone ppbv', $
        yrange=[10.0, 110.0]
      cgoplot, TGrid, mn+std, color='grey'
      cgoplot, TGrid, mn-std, color='grey'
      cgoplot, TGrid, meantropo, color='turquoise'
    endfor
    cglegend, length=0.0, tcolor=['grey','turquoise'],$
      titles=['+- 1 std dev','Average profile'], location=[.1,.89]
    !p.multi=0
    !p.charsize=2.0
    
    ; plot the ones which are tpflagged but not gradflagged:
    ; This should be good edge cases, we remove the edgecases with gradflags eq 1
    ;jj=where(tpflags eq 1 and gradflags eq 0)
    jj=[finalevents[0], removes[0]]
    egtitles=['good example','removed example']
    print, n_elements(removes), ' Events removed due to gradient being too shallow'
    foreach j, jj, ji do begin
      ; Plot one good and one bad example,
      ; plot is of event, density, filter, and selection of flux area
      ;
      ind = j
      evtrop=reform(troposphere[ind,*])
      evdens=reform(density[ind,*])
      evfilt=reform(SFTo3[ind,*])
      evtest=reform(testarr[ind,*])
      evloc=locs[ind]
      tpind=max(where(tgrid le mintp[ind]))
      ; bounds
      lbind = lbinds[ind]
      lb = lbs[ind] 
      ubind = ubinds[ind]
      ub = ubs[ind]
      
      cgdisplay, 900, 700, wid=ji, title=egtitles[ji]+string(years[ind])
      !p.multi=[0,1,3]
      !p.charsize=1.5
      !X.margin=[10,2]
      !Y.margin=[2,2]
      xrange=[0,mintp[ind]+1]
      
      ; plot profile
      ;
      cgplot, TGrid, evtrop, title='tropospheric profile', $
        xrange=xrange
      ; baseline of event
      cgoplot, [lb,ub], evtrop[[lbind,ubind]], color='maroon', $
        linestyle=2
      ; shade area of event
      
      ; molar density profile
      ;
      !Y.margin=[4,3]
      cgplot, TGrid[0:tpind], evdens[0:tpind], $
        xrange=xrange, $
        title='mol/cm3, flux='+string(flux[ind] / tropozone[ind] * 100.0) + '%'
      ; baseline of event
      cgoplot, [lb,ub], evdens[[lbind,ubind]], color='maroon', $
        linestyle=2, thick=2
      cgoplot, TGrid[lbind:ubind], evdens[lbind:ubind], color='maroon', $
        linestyle=2, thick=2
      
      ; bandfiltered profile
      ;
      cgplot, TGrid, evfilt, title='bandfiltered profile, tpflag+2gradflag='+string(tpflags[ind]+2*gradflags[ind]), $
        xtitle='altitude(km)', $
        xrange=xrange
      ; red line along threshold
      cgoplot, [2, mintp[ind]-1], fltarr(2)+cutoff, color='red'
      ; dashed line along zero axis
      cgoplot, xrange, fltarr(2), linestyle=2, color='grey'
      ; marks where event intersects zeros
      cgoplot, [lb,ub], evfilt[[lbind,ubind]], color='maroon', psym=1

      ;if j eq 3 then stop
      
      !p.multi=0
      !p.charsize=1.0
    endforeach
  endif
  
  
  
  return, events
END


