
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
 
PRO check_sondes

  sites=ptrarr(3)
  sites[0]=ptr_new(sondedata(/davis))
  sites[1]=ptr_new(sondedata(/macquarie))
  sites[2]=ptr_new(sondedata(/melbourne))
  titles=['Davis','Macquarie Island','Melbourne']
  
  cgps_open,'release_count.ps'
  cgdisplay, 800,1000, wid=1
  !p.charsize=2.5
  !p.font=0
  !p.multi=[0,1,3]

  ; split the sondes monthly, in order to count them.
  ; This is done in python hopefully(!)

  foreach site, sites, si do begin
  
    sondes=(*site)
    N=n_elements(sondes.jtime)
    sondejtime=sondes.jtime
  
    ; sort ppbv's into monthly bins
    caldat, sondejtime, M, D, Y
    
    ;===================================================
    ; at each site plot monthly sonde releases
    ;===================================================
    histo = HISTOGRAM(M, BINSIZE=1, LOCATIONS=binvals)
    cgBARPLOT, histo, /NOERASE, $
        TITLE=titles[si], $
        COLORS='black', $
        BARNAMES=['J','F','M','A','M','J','J','A','S','O','N','D'], $ 
        XMINOR=1.0, YMINOR=1.0
    !p.multi=[2-si,1,3]
    
  endforeach
  cgps_close
  !p.multi=0
END
