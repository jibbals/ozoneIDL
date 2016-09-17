;
;   PROCEDURE: eventdepth
;
;   PURPOSE:  
;           Look at distribution of event ozone peak depths below the tropopause
;
;   NOTES: 
;           no species of millipede has 1000 legs...
PRO eventdepth
  
  
  ; grab all the ozonesonde data
  ;
  data=ptrarr(3, /allocate_heap)
  data[0] = ptr_new(sondedata(/davis))
  data[1] = ptr_new(sondedata(/macquarie))
  data[2] = ptr_new(sondedata(/melbourne)) 

  bs=0.5
  bstring='(bin size='+string(bs, format='(i3)')+')'
  cgdisplay, 900, 700, wid=0, title='Event peak ozone depths' 
  
  colors= ['dark red', 'dark green', 'royal blue']
  sites=['Davis','Macquarie','Melbourne']
  ytitles=['','Frequency of occurrence','']
  xtitles=['','','Event Depth(km from tropopause)']
  thk=2.0
  
  !p.multi = [0,1,3]  
  ; for each station
  foreach ptr, data, ptri do begin
    dat=*ptr
    !Y.MARGIN=[4,2]

    ev=getevents(dat)
    jtimes=ev.jtime
    n=n_elements(jtimes)
    caldat, jtimes, Month, D, Y
    type = ev.type 
    depths=ev.tp-ev.locations
    ;pdf = HISTOGRAM(depths, BINSIZE=bs, locations=sbin, /nan) 
    
    ; set up plot
    
    cghistoplot, depths,$
      xtitle=xtitles[ptri], $
      title=sites[ptri]+'('+ string(n,format='(i3)')+' events)', $
      xstyle=8,ystyle=8, yminor=1, xminor=1, $
      charsize=3, xrange=[0,10], $
      binsize=bs, /fillpolygon, polycolor=colors[ptri]
    
  endforeach
   
   
!p.multi=0
!y.margin=0
END
