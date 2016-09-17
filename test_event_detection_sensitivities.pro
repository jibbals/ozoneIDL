
PRO test_event_detection_sensitivities
sites=ptrarr(3)
sites[0]=ptr_new(sondedata(/davis))
sites[1]=ptr_new(sondedata(/macquarie))
sites[2]=ptr_new(sondedata(/melbourne))
titles=['Davis','Macca','Melb']
CPS=[0.95,0.98,0.985,0.99,0.995]
SCALES=[ [ 0.5,5 ], [0.4, 5.1], [0.6,4.9], [0.,5.5] ]
foreach site, sites, si do begin
    ; sensitivity to cutoff percentile:
    foreach CP, CPS do begin
        print, "======================================"
        print, titles[si], " Cutoff percentile = ", CP
        tmp=getevents((*site), /verbose, CUTOFF_PERCENTILE=CP)
        
        print, "======================================"
        print, "======================================"
        print, n_elements(tmp.jtime), " Events after filtering with cutoff_percentile = ", CP
        print, '======================================'
    
    endforeach
endforeach
foreach site, sites, si do begin
    ; sensitivity to upper and lower filter bounds
    FOR i=0,2 DO BEGIN
        print, "======================================"
        print, titles[si], " Fourier scale = ", SCALES[*,i]
        tmp=getevents((*site), /verbose, FOURIER_SCALE=SCALES[*,i])
        
        print, "======================================"
        print, "======================================"
        print, n_elements(tmp.jtime), " Events after filtering using SCALES bounds = ", SCALES[*,i]
        print, '======================================'
    ENDFOR
endforeach

end

