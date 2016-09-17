;
;
PRO OzoneSeries, DATELIMITS=DATELIMITS, ERAFILE=ERAFILE

  ; Set Defaults
  if n_elements(ERAFILE) eq 0 then $
    ERAFILE="C:\Users\asp_transfer\Documents\Data\Pressure\pressure_20120101-20120430.nc"

  ; Read ERA-Interim Pressure Data:
  ; 
  sa_read_era_interim, pdata, FILE=ERAFILE
  
  if n_elements(DATELIMITS) eq 0 then $
    DATELIMITS=[ pdata.jultime[0], pdata.jultime[-1] ]
  
  ; Read the davis data
  ;
  ret=GetInterpolatedDavis([200,300,500],DATELIMITS)
  d_ppmv=ret.ppmv
  d_dates=ret.dates
  d_years=JesseDate(d_dates)
  
  ; index of data over Davis site ( 69S, 80E ) 
  ;
  tmp = Min(Abs(pdata.latitude + 69), latdex)
  tmp = Min(Abs(pdata.longitude - 80), londex)
  
  ; pressure is 200, 300, 500
  ; let's look at a timeline of the ozone at these pressures over davis
  ; 
  p200 = pdata.ozone_mmr[londex,latdex,0,*] * 1.0E6 ; let's use ppmv
  p300 = pdata.ozone_mmr[londex,latdex,1,*] * 1.0E6
  p500 = pdata.ozone_mmr[londex,latdex,2,*] * 1.0E6
  
  ; For the x axis I like using floating point years
  ; caldat converts jultime to month,day,year
  ;
  CALDAT, pdata.jultime[0], M, D, Y
  CALDAT, pdata.jultime[-1], Me, De, Ye
  print, "Start date (D,M,Y): ", D, M, Y
  print, "End date (D,M,Y): ", De, Me, Ye 
  date=JesseDate(pdata.jultime)
  
  ; plot the ERA data
  ;
  plot1=PLOT(date, p200, 'g', NAME='200hpa',xtitle='Year',ytitle='ozone PPMV')
  plot4=PLOT(d_years, d_ppmv[*,0], '--2g+', /OVERPLOT)
  plot2=PLOT(date, p300, /OVERPLOT, 'b', NAME='300hpa')
  plot5=PLOT(d_years, d_ppmv[*,1], '--2b+', /OVERPLOT)
  plot3=PLOT(date, p500, /OVERPLOT, 'm', NAME='500hpa')
  plot6=PLOT(d_years, d_ppmv[*,2], '--2m+', /OVERPLOT)
  
  ; add the legend
  lege = LEGEND(target=[plot1,plot2,plot3], /DATA, /AUTO_TEXT_COLOR)

END
