; Seasonal frequency distribution of O3 at a particular plevel
; default = 300hPa
PRO frequencyDist, plevel, xmax=xmax, binsize=binsize
  
  ; default plevel
  if n_elements(plevel) eq 0 then plevel = 300
  if n_elements(binsize) eq 0 then binsize = 15
  if n_elements(xmax) eq 0 then xmax = 650 - plevel
  
  ; grab all the ozonesonde data
  ;
  data=ptrarr(3, /allocate_heap)
  data[0] = ptr_new(sondedata(/davis))
  data[1] = ptr_new(sondedata(/macquarie))
  data[2] = ptr_new(sondedata(/melbourne)) 

  ; for each location we want the freq distr.
  ; array [ times, season, station ]
  o3dat=make_array( [ 1000, 4, 3], /FLOAT, VALUE=!values.f_nan)
  foreach ptr, data, ptri do begin
    dat=*ptr
    
    ;interpolate out the plevel we want:
    o3=interp_pressure(dat.o3ppbv, dat.pressure, plevel)
    n=n_elements(o3)
    
    ; split seasonally
    caldat, dat.jtime, M, D, Y
    season=ptrarr(4,/allocate_heap)
    season[0]=ptr_new(where(M eq 1 or M eq 2 or M eq 12)) ; summer
    season[1]=ptr_new(where(M eq 3 or M eq 4 or M eq 5))  ; autumn
    season[2]=ptr_new(where(M eq 6 or M eq 7 or M eq 8))  ; winter
    season[3]=ptr_new(where(M eq 9 or M eq 10 or M eq 11)); spring
   
    ; store seasonally split data
    for i=0,3 do $
      o3dat[*season[i], i, ptri] = o3[*season[i]]
    
  endforeach
  pstring=string(plevel,format='(i4)')+'hPa'
  bstring='(bin size='+string(binsize, format='(i3)')+')'
  cgdisplay, 1200, 800, wid=0, title='Relative frequencies of ozone ppbv' 
  !p.multi=[0,2,2]
  
  ; go through all data for each season and add to histoplots
  titles=['Summer','Autumn','Winter','Spring'] + ' ' + pstring + ' ' + bstring
  colors= ['dark red', 'dark green', 'royal blue']
  bases=['Davis','Macquarie','Melbourne']
  counts=fltarr(3)  
  for j=0,2 do $
    counts[j] = total(finite(o3dat[*,*,j]))
  print, 'plotting davis(red), macca(green), melbourne(blue)'
  bs=binsize
  for i=0,3 do begin
    !p.multi=[4-i,2,2]
    
    ; each season has 3 stations
    for j=0,2 do begin
       
      o3dati = HISTOGRAM(o3dat[*,i,j], BINSIZE=bs, MIN=0,locations=xbin,/nan)
      n_points=total(finite(o3dat[*,i,j]))
      cgplot, xbin, o3dati/n_points, $ ; BINSIZE=bs, MINinput=10, $
        title=titles[i],$
        XTITLE='ozone ppbv', $
        YTITLE='relative frequency', $
        YRANGE=[0,0.5], $
        XRANGE=[0,xmax], $
        color = colors[j], $
        thick = 2, $
        /NOERASE, OVERPLOT=J
    endfor
    
  endfor
  
  ; add legend!
  cglegend, tcolors=colors, $
    titles=bases+' (N='+string(counts,format='(i3)')+')', $
    length=0.0, $ ; remove legend lines
    charthick=1.5, $
    charsize=2., $
    vspace=2.4, $
    location = [0.24, 0.88]
  
  !p.multi=0

;  ; temporary checking of outlier data point
;  tmp=max(o3dat[*,1,1], maxind, /nan)
;  print, 'max:', max(o3dat), ' at ', maxind
;  maxp=(*data[1]).pressure[maxind,0:100]
;  maxppbv=(*data[1]).o3ppbv[maxind,0:100]
;  print, reform(maxp)
;  print, reform(maxppbv)
;
;  cgdisplay, 400, 600, wid=1
;  plotprofile, maxppbv, maxp, top=200
END 