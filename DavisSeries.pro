;
;
PRO DavisSeries

  start=JULDAY(1,1,2006)
  finish=JULDAY(3,1,2014)
  ret=GetInterpolatedDavis([200,300,500],[start,finish])
  vmr=ret.ppmv
  dates=ret.dates
  
  ; For the x axis I like using floating point years
  ; caldat converts jultime to month,day,year
  ;
  CALDAT, dates[0], M, D, Y
  CALDAT, dates[-1], Me, De, Ye
  print, "Start date (D,M,Y): ", D, M, Y
  print, "End date (D,M,Y): ", De, Me, Ye 
  years = JesseDate(dates)
  
  ; plot the data
  ;
  plot1=PLOT(years, vmr[*,0], 'g', NAME='200hpa',xtitle='Year',ytitle='ozone PPMV', title='Davis OzoneSondes')
  plot2=PLOT(years, vmr[*,1], /OVERPLOT, 'b', NAME='300hpa')
  plot3=PLOT(years, vmr[*,2], /OVERPLOT, 'm', NAME='500hpa')
  
  ; add the legend
  lege = LEGEND(target=[plot1,plot2,plot3], /DATA, /AUTO_TEXT_COLOR)

END
