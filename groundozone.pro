
;   Function: groundozone
;
;   Purpose: compare first layer of ozone between event and nonevent days
;   
;   KEYWORDS: trend
;       trend: plot trend instead of molecs/cm2 
;   Returns, Outputs: None
;   
;   NOTES:
;       I once ate a pocket watch, it was time consuming.
 
PRO groundozone, trend=trend

  sites=ptrarr(3)
  sites[0]=ptr_new(sondedata(/davis))
  sites[1]=ptr_new(sondedata(/macquarie))
  sites[2]=ptr_new(sondedata(/melbourne))
  titles=['Davis','Macquarie','Melbourne']

  ; Plot setup
  cgdisplay, 1200,900, wid=1
  !p.charsize=2.5
  !p.multi=[0,2,3]

  means=fltarr([2,3,12]) ; monthly meaned molecules per cm2 and pp
  gstds =fltarr([2,3,12]) ; 
  emeans=fltarr([2,3,12]) ; event monthly meaned molecs per cm2 and pp
  foreach site, sites, si do begin
  
    sondes=(*site)
    N=n_elements(sondes.jtime)
    sondejtime=sondes.jtime
    sondealts=sondes.gph / 1.0e3
    sondepps=sondes.o3pp
    sondeppbvs=sondes.o3ppbv
    sondepress = sondes.pressure  
  
    ; set up array of NaNs for holding tropo molecs/cm2
    groundmolecs = make_array(N, /float, value=!values.f_nan)
    groundpps = make_array(N, /float, value=!values.f_nan)
  
    ; sort ppbv's into monthly bins
    caldat, sondejtime, M, D, Y
  
    ;for each month
    for i=0,11 do begin
      ; finite indexes of month i
      mi = where(M eq (i+1) and $
                finite(sondealts[*,0]) and $
                sondealts[*,0] ne 0.0)
      ; For each tropopause in this month, get molecs/cm2
      ;
      foreach j, mi do begin
        ; data from 'ground' level
        ; pressure edges and ppbv edges are direct from sonde, interpolate mids
        pedges=sondepress[j,0:1]
        ppbvmids= (sondeppbvs[j, 0]+sondeppbvs[j,1]) / 2.0
        groundmolecs[j] = (ppbv_to_molecs_per_cm2( pedges, ppbvmids ))[0]

      endforeach ; end monthly loop
      ; also get pp 
      groundpps[mi] = sondepps[mi,0]
      means[0,si,i] = mean(groundmolecs[mi])
      means[1,si,i] = mean(groundpps[mi])
      gstds[0,si,i] = stddev(groundmolecs[mi])
      gstds[1,si,i] = stddev(groundpps[mi])
    endfor

    
    ; first column either plot trend or molecs/cm2
    if keyword_set(trend) then begin
      dating=where(sondejtime gt 2.40e6); ignore pre 1000BC
      ;dummy = LABEL_DATE(DATE_FORMAT=['%M','%Y'])
      dummy = LABEL_DATE(DATE_FORMAT='%Y')
      cgplot, sondejtime[dating], sondepps[dating,0], $ 
        xtitle='date', title=titles[si], $
        xminor=1, yminor=1, xstyle=8, ystyle=8, $
        XTICKUNITS = 'Time', XTICKINTERVAL=1, $
        XTICKFORMAT='LABEL_DATE'
    endif else begin
      title=titles[si]+ ' ground level ozone'
      ppbvscale=1.0e-15
      meansi=means[0,si,*]*ppbvscale
      gstdsi = gstds[0,si,*]*ppbvscale
      cgplot, meansi, color=cgcolor('black'), $
        yrange=[min(meansi)*0.99,max(meansi+gstdsi)*1.01], $
        xtickname=[' ','J','F','M','A','M','J','J','A','S','O','N','D',' '], $
        xticks=13, $
        xtitle='months', $
        ytitle='molecs/cm2 (/ 10^15)', $
        title=title, $
        xminor=1, yminor=1, xstyle=8, ystyle=8
      ; add grey stddevs
      cgoplot, meansi+gstdsi, color='grey'
      cgoplot, meansi-gstdsi, color='grey'
    endelse
    ;if si eq 0 then $
    ;  cglegend, titles=['nonevent','event', 'Mean +- 1 Std.'], $
    ;    tcolors=['black','red','grey'], $
    ;    length=0.0, $
    ;    location=[.12,.95], $
    ;    charsize=1.5

    meansi=means[1,si,*]
    gstdsi=gstds[1,si,*]
    title=titles[si] + ' ground level ozone partial pressure(mPa)'
    cgplot, meansi, $
      yrange=[min(meansi-gstdsi)*0.99,max(meansi+gstdsi)*1.01], $
      xtickname=[' ','J','F','M','A','M','J','J','A','S','O','N','D',' '], $
      xticks=13, $
      xtitle='months', $
      ytitle='pp', $
      title=title, $
      xminor=1, yminor=1, xstyle=8, ystyle=8
    cgoplot, meansi+gstdsi, color='grey'
    cgoplot, meansi-gstdsi, color='grey'

    ; overplot event pp
    ; color=cgcolor('red')         

  endforeach
  !p.multi=0
 
END
