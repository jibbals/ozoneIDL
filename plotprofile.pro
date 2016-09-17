;
;
;   If Altitude keyword set, assume km's for height
;     otherwise assume height is in hpa 
PRO plotprofile, pp, height, $
  tpo3=tpo3, title=title, OVERPLOT=OVERPLOT, $
  color=color, altitude=altitude, top=top, $
  tplr=tplr, temp=temp, Humid=Humid
  
  if N_ELEMENTS(color) eq 0 then color = cgcolor('black')
  if N_ELEMENTS(color) gt 0 then color = cgcolor(color)
  
  !y.margin=[4,4]
  !x.margin=[7.5,3]
  xcs=1.0 ;xcharsize
  if KEYWORD_SET(OVERPLOT) then begin
    oplot, pp, height, color=color
  endif else begin
    IF KEYWORD_SET(ALTITUDE) THEN BEGIN
      if N_ELEMENTS(top) eq 0 then top=max(height,/nan)
      pheight=height[where(height gt 0)] ; exclude 0 as it is the missingvalue.
      inds=where(pheight le top)
      cgplot, pp[inds], pheight[inds], charsize=2,$
        title=title, $
        xtitle="ppbv",ytitle="height(km)", $
        yrange=[0.,top], xrange=[1,150], $
        xticks=5, $
        color=color, background='white', $
        thick=2
       
      ; optionally plot temperature profile
      if N_ELEMENTS(temp) gt 0  then begin
        
        axis, 0,0.9*top, 0, /save, xstyle=1, xrange=[200, 300], $
          xtitle='Temp(K)', $
          charsize=xcs, xcharsize=xcs, $
          color=cgcolor('red')
        oplot, temp[inds], pheight[inds], color=cgcolor('red')
      endif
      
      ; optionally plot humidity profile
      if N_ELEMENTS(Humid) gt 0 then begin
        axis, 0,0.95*top, 0, /save, xstyle=1, xrange=[0, 100], $
          xtitle='RH', charsize=xcs, $
          xcharsize=xcs, color=cgcolor('blue')
        oplot, Humid[inds], pheight[inds], color=cgcolor('blue')
      endif
        
    ENDIF ELSE BEGIN
      if N_ELEMENTS(top) eq 0 then top=min(height,/nan)
      inds=where(height gt top)
      cgplot, pp, height , charsize=2, /YLOG,$
        title=title, $
        xtitle="pp",ytitle="height(hpa)", $
        yrange=[1000, top+5], xrange=[1,max(pp[inds],/nan)], $
        color=color, background='white', $
        thick=2
        
      ; optionally plot temperature profile
      if N_ELEMENTS(temp) gt 0 then begin
        offtemp = temp[inds]
        axis, 0,top+70, 0, /save, xstyle=1, xrange=[200, 300], $
          xtitle='Temp(K)', xcharsize=xcs, color=cgcolor('red')
        oplot, offtemp, height[inds], color=cgcolor('red')
      endif
      
      ; optionally plot humidity profile
      if N_ELEMENTS(Humid) gt 0 then begin
        hum = Humid[inds]
        axis, 0,top+20, 0, /save, xstyle=1, xrange=[0, 100], $
          xtitle='RH', xcharsize=xcs, color=cgcolor('blue')
        oplot, hum, height[inds], color=cgcolor('blue')
      endif
    ENDELSE
  endelse
  
  ; add the tropopause dashed line
  IF N_ELEMENTS(tpo3) gt 0 then oplot, !X.crange, [tpo3, tpo3], linestyle=2, color=cgcolor('black')
  ; add the lapserate tropopause dashed line
  IF N_ELEMENTS(tplr) gt 0 then oplot, !X.crange, [tplr, tplr], linestyle=2, color=cgcolor('red')
END
