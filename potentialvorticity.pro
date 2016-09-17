

PRO potentialvorticity, ERAFILE=ERAFILE

    ; Set Defaults
    IF n_elements(ERAFILE) EQ 0 THEN $
        ERAFILE="../Data/Davis/ERA/Pressure/pressure_20040101-20040430.nc"

    ; Read ERA-Interim Pressure Data:
    ; 
    sa_read_era_interim, pdata, FILE=ERAFILE
  
    ; index of data over Davis site ( 69S, 80E ) 
    ;
    tmp = Min(Abs(pdata.latitude + 69), latdex)
    tmp = Min(Abs(pdata.longitude - 80), londex)
    
    ; pressure is 200, 300, 500
    ; let's look at a timeline of the potential vorticity at these pressures over davis
    ; 
    pv200 = pdata.potential_vorticity[londex,latdex,0,*] * 1.0E6 ; in pvu's
    pv300 = pdata.potential_vorticity[londex,latdex,1,*] * 1.0E6
    pv500 = pdata.potential_vorticity[londex,latdex,2,*] * 1.0E6
   
    ; For the x axis I like using floating point years
    ; caldat converts jultime to month,day,year
    ;
    CALDAT, pdata.jultime[0], M, D, Y
    CALDAT, pdata.jultime[-1], Me, De, Ye
    print, "Start date (D,M,Y): ", D, M, Y
    print, "End date (D,M,Y): ", De, Me, Ye 
    date=JesseDate(pdata.jultime)
    ; Grab the location string from the ERAFILE string:
    loc = locfromstr(ERAFILE)
  
    ; plot the ERA data
    ;
    yaxis=[max([pv200,pv300,pv500]), min(pv200)-0.4]
    plot1=PLOT(date, pv200, 'g', NAME='200hpa',xtitle='Year',ytitle='potential vorticity',$
      yrange=yaxis, title="potential vorticity over "+loc)
    plot2=PLOT(date, pv300, /OVERPLOT, 'b', NAME='300hpa')
    plot3=PLOT(date, pv500, /OVERPLOT, 'm', NAME='500hpa')
    
    ; add the legend
    ;lege = LEGEND(target=[plot1,plot2,plot3], /DATA, /AUTO_TEXT_COLOR)
    ;cglegend, alignment=1, colors=['g','b','m'], $
    ;    titles=['200hPa','300hPa','500hPa']
END

