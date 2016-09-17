PRO FoldEvents, ERAFile=ERAFile

  ; get the ERA Data
  if N_ELEMENTS(ERAFile) eq 0 then $
    ERAFile="C:\Users\asp_transfer\Documents\Data\Pressure\pressure_20120101-20120430.nc"
  sa_read_era_interim, ERA, FILE=ERAFile
  
  ; index of data over Davis site ( 69S, 78E ) 
  tmp = Min(Abs(ERA.latitude + 69), latdex)
  tmp = Min(Abs(ERA.longitude - 78), londex)
  
  ; get the Davis Sondes data
  sondes=DavisSondes() ; all the data
  ret=GetInterpolatedDavis([300,400,500]) ; data interpolated at pressures 200, 300, 400
  d_ppmv=ret.ppmv
  d_years=JesseDate(ret.dates)
  
  ; First plot the series for our interpolated davis data
  WINDOW, 0, XSIZE=1200, YSIZE=800
  !P.multi = [0,3,2] ; 3 cols 2 rows
  !P.charsize=2.0
  PLOT, d_years, d_ppmv[*,0], title='300hpa', psym=1
  PLOT, d_years, d_ppmv[*,1], title='400hpa', psym=1
  PLOT, d_years, d_ppmv[*,2], title='500hpa', psym=1
  
  ; for the maximum ppmv ozone at each pressure let's check the profile
  ;
  tmp=Max(d_ppmv[*,0], ind1)
  print, "max 1 occurs on :", ret.dates[ind1]
  tmp=Max(d_ppmv[*,1], ind2)
  print, "max 2 occurs on :", ret.dates[ind2]
  tmp=Max(d_ppmv[*,2], ind3)
  print, "max 3 occurs on :", ret.dates[ind3]
  inds=[ind1,ind2,ind3]
  
  pressures=sondes.pressure[inds,*]
  ppbvs=sondes.ppmv[inds,*] * 1.0e3
  alts = sondes.altitude[inds,*]
  
  press1=pressures[0,*] & ppbv1=REFORM(ppbvs[0,*])
  tp1 = OzoneTropopause(ppbv1, REFORM(alts[0,*]),tp1ind) 
  press2=pressures[1,*] & ppbv2=REFORM(ppbvs[1,*])
  tp2 = OzoneTropopause(ppbv2, REFORM(alts[1,*]),tp2ind) 
  press3=pressures[2,*] & ppbv3=REFORM(ppbvs[2,*])
  tp3 = OzoneTropopause(ppbv3, REFORM(alts[2,*]),tp3ind)
    
  plot, ppbv1, press1 , /XLOG,$
    title="Sonde O3 Profile: "+string(d_years[ind1],format='(f7.2)'), $
    xtitle="ppbv",ytitle="hpa", $
    yrange=[1000.,100.], xrange=[1,max(ppbv1)], $
    ystyle=8 ; no right y axis
  ;yaxis = AXIS(YAXIS=1, YTICKV=reform(alts[0,*]), LOCATION='right', $
  ;  TITLE='Altitude(km)', YRANGE=[0, max(alts[0,*])])
  oplot, [1,max(ppbv1)], [press1[tp1ind], press1[tp1ind]], linestyle=2
  plot, ppbv2, press2 , /XLOG,$
    title="Sonde O3 Profile: "+string(d_years[ind2],format='(f7.2)'), $
    xtitle="ppbv",ytitle="hpa", $
    yrange=[1000.,100.], xrange=[1,max(ppbv2)]
  oplot, [1, max(ppbv2)], [press2[tp2ind], press2[tp2ind]], linestyle=2
  plot, ppbv3, press3 , /XLOG,$
    title="Sonde O3 Profile: "+string(d_years[ind3],format='(f7.2)'), $
    xtitle="ppbv",ytitle="hpa", $
    yrange=[1000.,100.], xrange=[1,max(ppbv3)]
  oplot, [1, max(ppbv3)], [press3[tp3ind], press3[tp3ind]], linestyle=2
 
  ; overplot horizontal line at tropopause
  ; TODO: use def from davis paper
  
  ; Total ozone below 400 hpa
  ; TODO: total column ozone function.
  
  !P.multi=0
  
END