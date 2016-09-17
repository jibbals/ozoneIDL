PRO sondeozonecontours, startyear=startyear
  ; get the Sondes data
  data=ptrarr(3,/allocate_heap)
  *data[0]=sondedata(/davis)
  *data[1]=sondedata(/macq)
  *data[2]=sondedata(/melb)
  
  titles=['Davis','Macquarie Island','Melbourne']
  ; contour pressure bounds:
  lbound=100 & ubound= 600
  missing=0.0
  
  ; contour colour bounds, just look at 25-200 o3 ppbv levels
  nlevels = 15
  interval = 15
  oob_low = 0
  oob_high = nlevels
  levels = Indgen(nlevels)*interval + interval ; ppbv colours
  range = [Min(levels), Max(levels)]

  ;load fresh blue - red colour table:
  cgLoadCT, 0 ; Start fresh
  cgLoadCT, 33, NColors=nlevels+1, Bottom=0

  ; for our three datasets
  for i=0,2 do begin
    cgdisplay, 1400, 500, wid=i, title=titles[i]
    datai=*data[i]
    ; convert julian to years
    dates=JesseDate(datai.JTIME)
    press=datai.pressure
    nz=n_elements(press[0,*])
    ; mPa/hPa to hPa/hPa
    vmr=datai.o3pp/press*1.0e-5
    ppbv=vmr*1.0e9
    
    ; remove data before start year(optional)
    if n_elements(startyear) gt 0 then begin
      keeps=where(dates ge startyear)
      dates = dates[keeps]
      ppbv = ppbv[keeps, *]
      press = press[keeps, *]
    endif
    
    ; make dates into 2d array for contouring
    dates2d=make_array(n_elements(dates),nz)
    for j=0,nz-1 do $
      dates2d[*,j] = dates
    ; clip bounds to reduce saturation
    ppbv[where(press lt lbound or press gt ubound)] = missing
    ;interval=25
    ;levels=indgen(20)*interval
    n_years=ceil(dates[-1] - dates[0])
    xtickv=indgen(n_years+1)+floor(dates[0])
    yinterval=100
    n_press= (ubound-lbound) / yinterval
    ytickv=indgen(n_press)*yinterval + lbound
    
    ; ================================================ 
    ; Plot contours of ppbv ozone against altitude 
    ; ================================================
    !p.charsize=1.8
    cgcontour, ppbv, dates2d, press, $
      position=[.07,.07, .93,.8], $
      levels=levels, missingvalue=missing, $
      xtitle='Year', $
      ;xtickv=xtickv, $
      ;xticklen=1, $
      xtickinterval=1, $
      xrange=[xtickv[0],xtickv[-1]], $ 
      ytitle='Pressure(hPa)', $
      yticks=n_elements(ytickv)-1, $
      ytickv=reverse(ytickv), $
      yrange=[ubound,lbound], $
      C_colors=indgen(nlevels+1), $
      /YLOG, /fill ;, /cell_fill
      
  
    ;cgColorbar, divisions=9, range=[levels[0],levels[-1]], /DISCRETE,  $
    ;  tickinterval=100, $
    cgColorBar, NColors=nlevels-1, Bottom=1, $
      OOB_Low=fix(oob_low), OOB_High=fix(oob_high), $
      /discrete, Range=range, $  
      Title=titles[i]+' Ozone (ppbv)', $
      TLocation='top', Format='(i3)', $
      tcharsize=2, charsize=1.6, $
      position = [.1, .86, .9, .91]
  endfor

  ; why not compare pressure and altitude example
  melb=*data[2]
  pressure=melb.pressure[-1,*]
  keep=where(pressure gt 0.0001)
  pressure=pressure[keep]
  altitude=melb.gph[-1,keep]/1.0e3 ; in km
  cgdisplay, 600, 600, wid=3, title='height vs pressure'

  cgplot, altitude, pressure, title='Melbourne example', $
    XTicklen=1.0, YTicklen=1.0, XGridStyle=1, YGridStyle=1, $
    xtitle='Altitude(km)', ytitle='Pressure(hPa)'
END
