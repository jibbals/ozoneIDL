; tptest.pro
;   Determine whether the event tropopause is significantly different 
;   to the non event tropopause
;
;   Use Welch's T Test before and after deseasonalising the data
;
pro tptest, showseason=showseason

    ; First get events and the minimum tropopause sets
    data = ptrarr(3, /allocate_heap)
    data[0] = ptr_new(sondedata(/davis))
    data[1] = ptr_new(sondedata(/macquarie))
    data[2] = ptr_new(sondedata(/melbourne)) 
    sites=['Davis','Macquarie','Melbourne']
    results=fltarr([3,2])

    foreach ptr, data, ptri do begin
        dat=*ptr
        ; get events and the event tropopauses
        ev=getevents(dat)
        evi = ev.indices
        etps = ev.tp
        ; get overall tropopauses
        tpo3 = dat.tpo3
        tplr = dat.tplr
        tps = min( [[tplr],[tpo3]], dimension=2, /nan )
        if n_elements(tps) ne n_elements(tpo3) then stop
        ; remove event and NAN tps from all tps set
        tps[evi] = !values.f_nan
        tps = tps[where(finite(tps))]

        ; Now we have the event and non event minimum tropopauses, do t test
        ;print "T-Test for "+sites[ptri]
        ;print 
        results[ptri, *] = TM_TEST(tps, etps,/UNEQUAL)
    endforeach
    print, "========= Not deseasonalised ========="
    print, sites
    print, 'T values:'
    print, string(results[*,0])
    print, 'p values:' 
    print, string(results[*,1])
    diffs = where(results[*,1] lt 0.05)
    print, 'event tropopauses are significantly different(alpha = 0.05) in'
    print, sites[diffs]

    if keyword_set(showseason) then begin
        !p.charsize=2.0
        ;set up plot
        title='Tropopause Heights'
        titles=['Davis','Macquarie','Melbourne']
        colours=['black','orange','purple']
        cgdisplay, 800,500, wid=0
        cgplot, intarr(12) ,/NODATA, color=cgcolor('black'), $
            yrange=[6, 15], $
            xtickname=[' ','J','F','M','A','M','J','J','A','S','O','N','D',' '], $
            xticks=13, xminor=1, xstyle=8, ystyle=8, $
            xtitle='Month', ytitle='Height(km)', $
            title=title
    endif
    ; Now deseasonalise, and then show seasonal pattern 
    ; and deseasonalised t-scores
    ;
    foreach ptr, data, ptri do begin
        dat=*ptr
        n = (size(dat.jtime))[1] ; how many profiles
        caldat, dat.jtime, m,d,y ; m contains month of datapoint
        ev = getevents(dat)
        caldat, ev.jtime, em, ed, ey
        evi = ev.indices
        etps = ev.tp
        newetps = etps
        tpo3 = dat.tpo3
        tplr = dat.tplr
        tps = min( [ [tplr],[tpo3] ], dimension=2, /nan )
        newtps = tps
        avg = fltarr(12) ; seasonal average
        eavg = fltarr(12) + !values.f_nan ; event seasonal avg
        for i = 1, 12 do begin
            mi = where(m eq i)
            if mi[0] eq -1 then stop ; every month should have some data
            emi = where(em eq i)
            avg[i-1] = mean(tps[mi],/nan)
            newtps[mi] = newtps[mi] - avg[i-1]
            if emi[0] eq -1 then continue
            newetps[emi] = newetps[emi] - avg[i-1]
            eavg[i-1] = mean(etps[emi],/nan)
        endfor
        if keyword_set(showseason) then begin
            ; plot the seasonal tpheights
            cgoplot, avg, color=colours[ptri], thick=3
            ; plot the avg event tp height
            cgoplot, eavg, color=colours[ptri], $
                linestyle=2, thick=2
        endif
        
        ; run t-test
        newtps[evi] = !values.f_nan
        newtps = newtps[where(finite(newtps))]
        results[ptri, *] = TM_TEST(newtps,newetps,/UNEQUAL)
    endforeach
    if keyword_set(showseason) then begin
        cglegend, titles=titles, charsize=2.5,$
            tcolors=colours, $
            length=0.0, vspace=2.5,$
            location=[0.60, 0.83]
        
    endif

    print, "=========  deseasonalised ========="
    print, sites
    print, 'T values:'
    print, string(results[*,0])
    print, 'p values:' 
    print, string(results[*,1])
    diffs = where(results[*,1] lt 0.05)
    print, 'event tropopauses are significantly different(alpha = 0.05) in'
    print, sites[diffs]

END
