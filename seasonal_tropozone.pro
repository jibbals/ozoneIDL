
;   Function: check_sondes.pro
;
;   Purpose: create some plots looking at sonde release dataset
;
;   Inputs:
;
;   Returns:
;   
;   Outputs:
;   
;   KEYWORDS:
;
;   NOTES:
;       Salmon is the colour of desire
 
FUNCTION seasonal_tropozone

  sites=ptrarr(3)
  sites[0]=ptr_new(sondedata(/davis))
  sites[1]=ptr_new(sondedata(/macquarie))
  sites[2]=ptr_new(sondedata(/melbourne))
  titles=['Davis','Macquarie Island','Melbourne']

  cgdisplay, 800,1000, wid=1
  !p.charsize=2.5
  !p.multi=[0,1,3]

  ; split the sondes monthly, in order to count them.
  ; This is done in python hopefully(!)
  seasons=make_array([3,12],/float,value=!values.f_nan)
  means=fltarr([3,12])
  stds=fltarr([3,12])
  meds=fltarr([3,12])
  tps=fltarr([3,12])
  ppbvs = fltarr([3,12, 200]) ; ppbv average seasonal up to 20km
  hppbvs = findgen(200)/10.0  ; ppbvs heights

  foreach site, sites, si do begin
  
    sondes=(*site)
    N=n_elements(sondes.jtime)
    sondejtime=sondes.jtime
    sondealts=sondes.gph / 1.0e3
    sondeppbvs=sondes.o3ppbv
    
    sondetplrs = sondes.tplr
    sondetpo3s = sondes.tpo3
    sondetps = min([[sondetplrs],[sondetpo3s]],dimension=2,/nan)
    sondetpps = max([[sondes.tpo3p], [sondes.tplrp]],dimension=2 ,/nan)
    sondepress = sondes.pressure  
  
    ; set up array of NaNs for holding tropo molecs/cm2
    sondemolecs = make_array(N, /float, value=!values.f_nan)
  
    ; sort ppbv's into monthly bins
    caldat, sondejtime, M, D, Y
    
    ;===================================================
    ; at each site plot monthly sonde releases
    ;===================================================
    histo = HISTOGRAM(M, BINSIZE=1, LOCATIONS=binvals)
    !p.multi=[si,1,3]
    cgBARPLOT, histo, /NOERASE, $
        TITLE=titles[si], $
        COLORS='black', $
        BARNAMES=['J','F','M','A','M','J','J','A','S','O','N','D'], $ 
        XMINOR=1.0, YMINOR=1.0
    
    
    ;for each month
    for i=0,11 do begin
      ; finite indexes of month i
      mi = where(M eq (i+1) and $
                sondetpps gt 2 and $
                finite(sondetpps) )
      ; For each tropopause in this month, get molecs/cm2
      ;
      foreach j, mi do begin
        ; data from below tropopause
        tpi = Value_Locate(sondepress[j,*], sondetpps[j])
        pedges = sondepress[j, 0:tpi]
        ; pressure edges and ppbv edges are direct from sonde, interpolate mids
        ppbvmids = (sondeppbvs[j, 0:tpi-1] + sondeppbvs[j,1:tpi]) / 2.0
        sondemolecs[j] = total(ppbv_to_molecs_per_cm2( pedges, ppbvmids ))
      endforeach

      ; for each height level get average ppbv
      ;
      monthalts= sondealts[mi,*]
      monthppbvs=sondeppbvs[mi,*]
      foreach ub, hppbvs, ui do begin
        if ui eq 0 then continue    ; skip first row(zeros anyway)
        lb = hppbvs[ui-1] 
        sinds = where(monthalts le ub and monthalts gt lb)
        if sinds[0] eq -1 then ppbvs[si,i,ui] = !values.f_nan $
          else ppbvs[si,i,ui] = mean(monthppbvs[sinds],/nan)
      endforeach

      tps[si,i] = mean(sondetps[mi],/nan)
      means[si,i] = mean(sondemolecs[mi],/nan)
      stds[si,i] = stddev(sondemolecs[mi],/nan)
      meds[si,i] = median(sondemolecs[mi])
    endfor

    if keyword_set(analyse) then begin
        title='seasonal tropospheric ozone at '+titles[si]
        meansi=means[si,*]
        stdsi =stds[si,*]
        
        cgplot, meansi, color=cgcolor('black'), $
          yrange=[min(meansi-stdsi)*0.99,max(meansi+stdsi)*1.01], $
          xtickname=[' ','J','F','M','A','M','J','J','A','S','O','N','D',' '], $
          xticks=13, $
          xtitle='months', $
          ytitle='molecs/cm2', $
          title=title, $
          xminor=1, xstyle=8, ystyle=8
        ; add grey stddevs
        cgoplot, meansi+stdsi, color='grey'
        cgoplot, meansi-stdsi, color='grey'
    
        cgoplot, meds[si,*], color=cgcolor('red')
        if si eq 0 then $
          cglegend, titles=['Mean','Median', 'Mean +- 1 Std.'], $
            tcolors=['black','red','grey'], $
            length=0.0, $
            location=[.12,.95], $
            charsize=1.5
      
    endif
  endforeach
  if keyword_set(contourppbv) then begin
    !p.multi=[0,1,3]
    !p.charsize=3.5
    !X.margin=[8,22]
    ; set up colourbar  
    ; contour colour bounds, just look at base to top o3 ppbv lvls
    range=[10,100]
    nlevels=91
    levels=indgen(nlevels)+10
    cgloadct, 0 ; start fresh
    cgloadct, 33, NColors=nlevels-1, Bottom=1   ; load colours
    oob_low=0
    oob_high=nlevels
    TVLCT, cgcolor('white',/triple), oob_low ; load white for oob
    TVLCT, cgcolor('white',/triple), oob_high

    cgdisplay, 900,900, wid=2
    
    xtitles=['','','months']
    for si = 0,2 do begin
      title='Seasonal ozone at '+titles[si]
      ppbvsi=reform(ppbvs[si,*,*])
      cgcontour, ppbvsi, findgen(12)+1., hppbvs, $
        yrange=[0.0, 15.0], $
        xtickname=[' ','J','F','M','A','M','J','J','A','S','O','N','D',' '], $
        xticks=13, $
        xtitle=xtitles[si], $
        ytitle='Altitude(kms)', $
        title=title, $
        yminor=1, xminor=1, xstyle=8, ystyle=8, $
        /fill , levels = [Min(ppbvsi), levels], c_colors=indgen(nlevels+1)

      ; add tropopause lines
      cgoplot, findgen(12)+1.0, tps[si,*], color='black'
      
    endfor
    cgColorBar, NColors=nlevels-1, bottom=1,$ ;nlevels-base, 
      OOB_Low=fix(oob_low), OOB_High=fix(oob_low), $; high bound out of range gives white color
      Range=range, $
      Title='Ozone (ppbv)', $
      TLocation='left', Format='(i3)', $
      tcharsize=3.0, charsize=2.5, $
      /vertical, position=[.9,.25, .95,.75]
  endif
  !p.multi=0
  return, means
END
