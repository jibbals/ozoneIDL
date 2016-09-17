
PRO examineprofile, jday, $
  davis=davis, macquarie=macquarie, melbourne=melbourne, $
  highres=highres

  ; get the ERA Data
  ERA=eradata(jday+[-.25,0,.25], $
    davis=davis, macquarie=macquarie, melbourne=melbourne, $
    highres=highres)
  if not finite(era.jultime[1]) then begin
    print, 'No matching ERA data'
    return
  endif 
  sondes=sondedata(davis=davis, macquarie=macquarie, melbourne=melbourne)
  tmp = min(abs(sondes.jtime - jday), sondeind)
  baseloc=[sondes.lon, sondes.lat]
  
  caldat, sondes.jtime[sondeind], M, D, Y, H
  
  if tmp gt 0.01 then begin 
    print, "sonde event missing?, tmp=",tmp
    help, jday
    help, sondes.jtime
    print, sondes.jtime[sondeind]
  endif
  
  title="sonde:"+string(Y,format='(i04)')+string(M,format='(i02)')+$
      string(D,format='(i02)')+string(H,format='(i02)')
  print, title
  tmp = min(abs(era.jultime - jday), eraind)
  if tmp gt 0.01 then print, "era profile days from sonde event:",tmp
  
  ; sonde arrays we will examine
  sondeppbv=sondes.o3ppbv    ; mPa to hPa and vmr to ppbv
  sondealt=sondes.gph/1.0e3   ; m to km
  sonderh=sondes.rh
  sondetemp=sondes.temperature+273.15   ; C to K
  sondeppbv=reform(sondeppbv[sondeind,*])
  sondealt=reform(sondealt[sondeind,*])
  sondepressure=reform(sondes.pressure[sondeind,*])
  sondetemp=reform(sondetemp[sondeind,*])
  sonderh=reform(sonderh[sondeind,*])
  sondetp_alt=sondes.tpo3[sondeind]
  sondetplr_alt=sondes.tplr[sondeind]
  
  ; ERA arrays we will examine
  ;
  tmp = Min(Abs(ERA.latitude - baseloc[1]), Y)
  tmp = Min(Abs(ERA.longitude - baseloc[0]), X)
  ERAppbv = o3mmrtovmr(reform(era.ozone_mmr[X,Y,*,*]))*1.0e9
  ERApressure = ERA.pressure
  ERAalt = reform(ERA.geopotential[X,Y,*,*] / 1.0e3)
  ERArh = reform(ERA.relative_humidity[X,Y,*,*])
  ERAtemp=reform(ERA.temperature[X,Y,*,*])
  ERApv=reform(ERA.potential_vorticity[X,Y,*,*])

;  Not enough resolution for ozonetp definition to work
;  ERAtp = fltarr(3)
;  for ii=0,2 do $
;    ERAtp[ii] = OzoneTropopause(ERAppbv[*,ii], ERAalt[*,ii])
  
  ; Plotting
  cgdisplay, 800, 1000, wid=0
  !P.Multi = [0,2,2]
  
  top=sondetp_alt+5
  PlotProfile, sondeppbv, sondealt, tpo3=sondetp_alt, tplr=sondetplr_alt, $
    title=title, top=top, temp=sondetemp, humid=sonderh, /altitude
  
  cgplot, erapv[*,1]*1.0e6, eraalt[*,1], title='ERA PVU', $
    xrange=[0, -5], yrange=[0,sondetp_alt+5]

  for i=1,2 do begin
    caldat, ERA.jultime[i], M,D,Y,H
    titles="ERA: "+string(Y,format='(i04)')+string(M,format='(i02)')+$
      string(D,format='(i02)')+string(H,format='(i02)')

    PlotProfile, ERAppbv[*,i], ERAalt[*,i], $
      title=titles, temp=ERAtemp[*,i], $
      humid=ERArh[*,i], top=top, $
      /altitude
    
  endfor
  !P.Multi=0
  
  ; can't remember what I was doing here:
  ;medians=monthlymean(sondes,/medians)
  ;thresh=0.45*medians[M-1]
  ;cgdisplay, 800, 800, wid=1
  ;print, isevent(sondealt, sondeppbv, sondeTP_alt, thresh=thresh,/analyse)
  
END
