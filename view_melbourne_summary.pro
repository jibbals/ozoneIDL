;
;
pro view_melbourne_summary, plotcontour=plotcontour
  ; read the melbourne data
  melb=readmelbournecsv(file='../Data/Broadmeadows/Sondes/melbsondes.dat')

  ; lat = -37.68
  ; lon = 144.95
  ; height = 110m

  dates=JesseDate(melb.jtime)
  dates[where(dates lt 1980)] = !values.f_nan
  N_t = n_elements(dates)
  N_p = n_elements(melb.pressure[0,*])
  dates2d = fltarr(N_t, N_p)
  for j=0, N_p-1 do $
    dates2d[*,j] = dates
  
  levels = indgen(17)
  xtickv=1999+indgen(15)
  ytickv = indgen(11)*100
  press = melb.pressure
  press[where(press le 0.0001)] = !values.f_nan
  
  CT_LOAD, 'CALIPSO'

  if keyword_set(plotcontour) then begin
    cgdisplay, 1200, 800, wid=0
    
    cgcontour, melb.o3pp, dates2d, press, /YLOG, $
      /cell_fill, position=[.08,.05, .98,.8], $
      levels=levels, missingvalue=0.0, $
      xtickv=xtickv, xticklen=1, xtickinterval=1, $
      ytitle='Pressure(hpa)', xtitle='Year', $
      yrange=[1000,5]
    
  
    cgColorbar, divisions=9, range=[levels[0],levels[-1]], /DISCRETE,  $
            tickinterval=2, $
            Title='O3 PP', $
            TLocation='top', Format='(i2)', $
            tcharsize=1.5, charsize=1, $
            position = [.1, .85, .9, .9]
  endif
  ; Now let's look at the 200, 300, 400, 500 hpa ozone time series
  
  pressures=[200,300,400,500]
  ppseries=interp_pressure(melb.o3pp, melb.pressure,pressures)
  
  
  cgdisplay, 1000, 800, wid=1
  K=n_elements(pressures)
  !p.multi=[0,(K+1)/2,(K+1)/2]
  maxinds=fltarr(K)
  titles=strarr(K)
  for i=0,K-1 do begin
    cgplot, dates, ppseries[*,i], xtitle='years',ytitle='o3 PP', $
      title='o3 PP at '+string(pressures[i], format='(i3)')+'hpa', $
      psym=1
      
    tmp=max(ppseries[*,i], maxind)
    maxinds[i]=maxind
    caldat, melb.jtime[maxind], M, D, Y, H
    fi2='(i2)'
    fi5='(i5)'
    titles[i]='(D M Y H)' + string(D,format=fi2)+string(M,format=fi2)+$
      string(Y, format=fi5)+string(H,format=fi2)
    print, string(pressures[i])+'hpa has max o3 pp on '+titles[i]
  endfor

  cgdisplay, 1200, 800, wid=2
  !p.multi=[0,K,1]
  for i=0,K-1 do begin
    press=melb.pressure[maxinds[i],*]
    press[where(press lt 0.2)] = !values.f_nan
    PlotProfile, melb.o3pp[maxinds[i],*], press, $
      title=titles[i]
  endfor
  
  !p.multi=0
  
  ;one more to check gph profile
  if total(melb.gph[maxinds[1],0:3]) gt 100 then begin
    cgdisplay, 700, 700, wid=3
    gph=reform(melb.gph[maxinds[1],*]/1000.0)
    gph[where(gph lt 0.001)] = !values.f_nan
    ; vmr is PP(mPa) / P(hPa) * 1/100 000
    vmr=melb.o3pp[maxinds[1],*] / melb.pressure[maxinds[1],*] * 1.0e-5
    ppbv=reform(vmr*1.0e9)
    tp=OzoneTropopause(ppbv, gph)
    temp=reform(melb.temperature[maxinds[1],*])
      PlotProfile, ppbv, gph, tp=tp, $
        /ALTITUDE, title=titles[1], top=tp+5, $
        temp=temp
    
    ; TODO: Remove
    stop
  endif
  
end
