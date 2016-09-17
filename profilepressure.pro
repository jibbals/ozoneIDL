;
;   Read ppbv, pressure, and more arguments and plot a standard profile
;
;   ARGUMENTS: ppbv, pressure
;
;   OPTIONALS: 
;       title: string setting plot title
;       top: float setting upper y axis
;       tplr: lapse rate tropopause
;       tp03: ozone tropopause
;       temp: temperature profile(Kelvin)
;       humid: relative humidity profile
;
;   NOTES:
;       
PRO profilepressure, ppbv, pressure, $
  tpo3=tpo3, title=title, $
  color=color, top=top, $
  tplr=tplr, temp=temp, Humid=Humid
  
  if N_ELEMENTS(color) eq 0 then color = cgcolor('black')
  if N_ELEMENTS(color) gt 0 then color = cgcolor(color)
  if N_ELEMENTS(top) eq 0 then top=300 ; 200hPa is pretty high

  !y.margin=[4,4]
  !x.margin=[8,3]
  xcs=1.0 ;xcharsize

  ; Set up yaxis (because idl farts on log axes)
  ticks = (indgen(11-top/100)+top/100)*100
  if ticks[0] eq 0 then ticks[0] = top
  nticks = N_Elements(ticks) 
  

  inds=where(pressure gt top)
  ppbv1 = ppbv[inds]
  pressure1 = pressure[inds]
  cgplot, ppbv1, pressure1, charsize=2, /YLOG,$
      title=title, $
      xtitle="O3 ppbv",ytitle="Pressure(hPa)", $
      yrange=[1000, top], xrange=[20,200], $
      color=color, background='white', $
      thick=2, $
      yticks=nticks-1, ytickv=Reverse(ticks), $
      XSTYLE = 8, YSTYLE = 8 ; top and right axes removed
        
  ; add the tropopause dashed line
  IF N_ELEMENTS(tpo3) gt 0 then $
    oplot, !X.crange, [tpo3, tpo3], linestyle=2, color=cgcolor('green')
  ; add the lapserate tropopause dashed line
  IF N_ELEMENTS(tplr) gt 0 then $
    oplot, !X.crange, [tplr, tplr], linestyle=2, color=cgcolor('red')
  

  ; optionally plot temperature profile
  if N_ELEMENTS(temp) gt 0 then begin
    temp1 = temp[inds]
    axis, 0,top+40, 0, /save, xstyle=1, xrange=[200, 300], $
      xtitle='Temp(K)', xcharsize=xcs, color=cgcolor('red')
    oplot, temp1, pressure1, color=cgcolor('red')
  endif
      
  ; optionally plot humidity profile
  if N_ELEMENTS(Humid) gt 0 then begin
    humid1 = Humid[inds]
    axis, 0,top+10, 0, /save, xstyle=1, xrange=[0, 100], $
      xtitle='RH', xcharsize=xcs, color=cgcolor('blue')
    oplot, humid1, pressure1, color=cgcolor('blue')
  endif
  

END
