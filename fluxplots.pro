FUNCTION regrid_profile, T, F, delta, top, $
  TGrid=TGrid
;   Purpose:
;     Regrid profile from 2km up to tropopause with interpolated regular grid spacing
;     
;   Inputs:
;     T: altitude of profile
;     F: ppbv of profile
;     
;   Returns:
;     Regridded ppbv
;     
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

PRO fluxplots, melbourne=melbourne, davis=davis, macquarie=macquarie
;   NAME:
;       fluxplots
;       Copy of getevents focussed on plotting all the flux calculations
;       
;   Inputs:
;       sondes: structure returned from sondedata() function
;   
;   Outputs: none
;   
;

  sondes = sondedata(melbourne=melbourne, davis=davis, macquarie=macquarie)

  ppbv=sondes.o3ppbv
  vmrs=1.0d-9*ppbv
  alt=sondes.gph/1.0e3  ; m to km
  alts=alt
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
  TGrid=findgen(regridlen)*delta
  troposphere = make_array(N_profs, regridlen, /FLOAT, VALUE=!values.f_nan )
  SFTo3 = make_array(N_profs, regridlen, /FLOAT, VALUE=!values.f_nan )
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
  density = reVMR * reP / R / reT * 6.02214129d23 ; in molecules/cm3
  
  low=.5/delta           ;km freq min
  high=5/delta           ;km frequency max
  
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
    ;SFTo3[i,tpup] = !values.f_nan
    troposphere[i,tpup] = !values.f_nan
    inds=where(finite(troposphere[i,*]))
    if inds[0] eq -1 then continue
    
    ; run the bandpass filter on the column up to the tropopause
    ;
    sa_filter_1d, reform(troposphere[i,inds]), low, high, filtered
    SFTo3[i,inds] = filtered  ; Simon's Fourier Transformed O3
 
    ; remove the ppbv above the tropopause - 1km
    ; and then Determine peak of filtered trimmed data
    tp1up=where(TGrid ge tpbound-1)
    filtered[tp1up] = !values.f_nan
    filtered[where(TGrid lt 2.0) ] = !values.f_nan
    SFTo3[i,tp1up] = !values.f_nan
    
    ; location of ozone peak in profile
    tmp=max(filtered,maxind,/nan)
    locs[i]=tgrid[maxind]
 
  endfor
  
  ; Clip the bottom 2km of the filtered data(transform noise)
  SFTo3[*,where(TGrid lt 2.0)] = !values.f_nan
  
  
  ; Percentiles of filtered ozone profile datapoints
  pcs=[.5,.75,.9,.95, .96, .97, .98, .99]
  percentiles=cgpercentiles(SFTo3[where(finite(SFTo3))], percentiles=pcs)
  cutoff=percentiles[-1]; we use the 99th percentile
  
  ; remove values below the cutoff
  testarr=SFTo3
  testarr[where(testarr lt cutoff)] = !values.f_nan
  ;testcounts=total(finite(testarr), 2)  ; profiles containing points above cutoff
  ; find the profile indices where values exist above the cutoff
  testprof=total(testarr,2,/NAN)
  eventinds=where(testprof gt 0.0)
  n_events=n_elements(eventinds)
  
  ; DETERMINE FLUX CAPACITOR
  ; 
  ; for each event
  tpcount=0 & gradcount=0
  foreach ind, eventinds, ii do begin
    evtrop=reform(troposphere[ind,*])
    evfilt=reform(SFTo3[ind,*])
    evdens=reform(density[ind,*])
    evloc=locs[ind]
    evlocind=where(TGrid eq evloc)
    ; lower zero bound is first zero below 'event' ozone values
    lbinds[ind] = max(where(evfilt le 0 and tgrid lt evloc))
    if lbinds[ind] eq -1 then begin
      lbinds[ind] = where(tgrid eq 2.0)
      print, 'event ', ind, ' at ', evloc 
      lbflags[ind]=1
    endif
    lbs[ind] = tgrid[lbinds[ind]]
    ; upper bound is first zero above 'event' or else tropopause - 1 
    ubinds[ind] = min(where(evfilt le 0 and tgrid gt evloc))
    if ubinds[ind] eq -1 then begin
      tmp=min(evtrop[evlocind:*],minind,/nan)
      ubinds[ind]=minind+evlocind
      tpflags[ind]=1
      tpcount=tpcount+1
    endif
    ubs[ind] = tgrid[ubinds[ind]]

    ; remove events where the gradient from peak to the min before the tp is not steep enough
    ; NB: events removed after foreach loop
    ;
    dx=TGrid[ubinds[ind]]-TGrid[evlocind]
    dy=evtrop[ubinds[ind]]-evtrop[evlocind]
    grad = dy/dx
    if dx lt .1 then grad = 0
    if grad gt -20 then begin ; if profile has gt -20 ppbv/km gradient flag it
      gradcount = gradcount+1
      gradflags[ind]=1
    endif   

    ; line between where filtered zeros occur above and below the ozone peak(using density profile)
    baseline=interpol(evdens[[lbinds[ind],ubinds[ind]]],[lbs[ind],ubs[ind]], tgrid[lbinds[ind]:ubinds[ind]] )
    area=evdens[lbinds[ind]:ubinds[ind]] - baseline
    ; conservative estimate for event flux is integral over baseline 
    ; of 'event'
    ; Do this using density column!
    ;
    flux[ind]=int_tabulated(tgrid[lbinds[ind]:ubinds[ind]],area)
    tropozone[ind] = int_tabulated(tgrid, evdens)
  endforeach
  
  ; remove gradient shallow flagged events
  ;
  removes=where(gradflags eq 1 and tpflags eq 1)
  finalevents = cgsetdifference(eventinds, removes)
  
  ; Determine event types based on distance from tropopause
  ; 0-2km = shallow, 2-4 = medium, 4+ = deep
  ;
  distance=mintp-locs
  bins=[2,4]  ; 0 to 2 is shallow, 2 to 4 is mediumn, 4+ is deep event
  type = (distance[finalevents] gt bins[0]) + (distance[finalevents] ge bins[1])
  
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; PLOTS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; show the profile
  ylims=[20,200] ; ppbv limits
  xlims=[1,13.5]
  ylims2=[1.0e11, 1.4e12] ; density limits

  kap = ceil(sqrt(n_elements(finalevents)))
  !P.multi = [0, kap, kap]
  ; Density profile plot
  cgdisplay, 1800, 1000, wid=0
 
  FOREACH evind, finalevents, ii DO BEGIN
    ; Plot out bandfilter and flux image for all events
    ;
    xaxi = where(tgrid ge 2 and tgrid le 14)
    xax  = tgrid[xaxi]
    evtime = reform(jtimes[evind, *])
    evppbv = reform(troposphere[evind, *])
    evP = reform(sondes.pressure[evind,*])
    evtp= mintp[evind]
    evFT = reform(SFTo3[evind,*])
    evdens=reform(density[evind,*])

    cgplot, xax, evdens[xaxi], $
      xrange=xlims, yrange=ylims2, $
      xstyle = 4, ystyle = 4 ; Removes axes
   
    ; line at tropopause: dashed, vertical, red
    cgoplot, [evtp, evtp], [min(evdens[xaxi]), max(evdens[xaxi])], $
      linestyle=2, color=2

    ; flux calculations
    evlbind = lbinds[evind]
    evubind = ubinds[evind]
    evlb = lbs[evind]
    evub = ubs[evind]
    evex = tgrid[evlbind:evubind]
    evbaseline=interpol(evdens[[evlbind,evubind]],[evlb,evub], evex )
    evtopline =evdens[evlbind:evubind]
    
    cgoplot, tgrid[evlbind:evubind], evtopline, $
      linestyle=2, color=3, thick=2
    cgoplot, tgrid[evlbind:evubind], evbaseline, $
      linestyle=2, color=3, thick=2

  

  ENDFOREACH

  !p.multi=0

  ; NOW LOOK AT A SUBSET IN MORE DETAIL

  !P.multi = [0, 1, 3]
  !P.charsize=2
  ; Density profile plot
  kap=9
  kapoff=5
  FOREACH evind, finalevents[kap:kap+kapoff], ii DO BEGIN

    cgdisplay, 900, 600, wid=ii+1
    ; Put out an example of the bandfilter
    ;
    xaxi = where(tgrid ge 2 and tgrid le 14)
    xax  = tgrid[xaxi]
    evtime = jessedate(jtimes[evind])
    evppbv = reform(troposphere[evind, *])
    evP = reform(sondes.pressure[evind,*])
    evtp= mintp[evind]
    evFT = reform(SFTo3[evind,*])
    evdens=reform(density[evind,*])

    ; event bounds
    ;
    evlbind = lbinds[evind]
    evubind = ubinds[evind]
    evlb = lbs[evind]
    evub = ubs[evind]
    evex = tgrid[evlbind:evubind]
    evbaseline=interpol(evdens[[evlbind,evubind]],[evlb,evub], evex )
    evtopline =evdens[evlbind:evubind]
    
    ; conservative estimate for event flux is integral over baseline 
    ; of 'event'
    ; Do this using density column!
    ;
    area=evdens[evlbind:evubind] - evbaseline
    evflux=int_tabulated(tgrid[evlbind:evubind],area)
    evtpind = max(where(tgrid lt evtp))
    evtropozone = int_tabulated(tgrid[0:evtpind], evdens[0:evtpind])  
    

    ; show the profile
    cgplot, xax, evppbv[xaxi], title='ozone ppbv', $
        xrange=xlims, yrange=ylims

    ; Density profile plot
    cgplot, xax, evdens[xaxi], $
      title='molecules/cm3 ' + $
      string(evflux/evtropozone*100) + "% stratospheric", $
      xrange=xlims, yrange=ylims2
    
    ; line at tropopause: dashed, vertical, red
    cgoplot, [evtp, evtp], [min(evdens[xaxi]), max(evdens[xaxi])],$
      linestyle=2, color=2

    cgoplot, tgrid[evlbind:evubind], evtopline, $
      linestyle=2, color=3, thick=2
    cgoplot, tgrid[evlbind:evubind], evbaseline, $
      linestyle=2, color=3, thick=2

    ; Add plot of transformed profile, against zeroline
    cgplot, tgrid, evFT, $
      title='Fourier transformed profile', xtitle='Altitude(km)', $
      xrange=xlims
    cgoplot, xlims, [0,0], color=2, linestyle=2
    cgoplot, xlims, [cutoff, cutoff], color=3, linestyle=2
    

  ENDFOREACH

  !p.multi=0

END


