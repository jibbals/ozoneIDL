;
;   PROCEDURE: eventtypes
;
;   PURPOSE:  
;           Look at distribution of event types: shallow, medium, and deep
;
PRO eventtypes
  
  
  ; grab all the ozonesonde data
  ;
  data=ptrarr(3, /allocate_heap)
  data[0] = ptr_new(sondedata(/davis))
  data[1] = ptr_new(sondedata(/macquarie))
  data[2] = ptr_new(sondedata(/melbourne)) 

  bs=1
  bstring='(bin size='+string(bs, format='(i3)')+')'
  cgdisplay, 900, 700, wid=0, title='Relative frequencies of event types' 
  
  colors= ['dark red', 'dark green', 'royal blue']
  styles=[1,2,0]
  bases=['Davis','Macquarie','Melbourne']
  ytitles=['','Frequency of occurrence','']
  xtitles=['','','Month']
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
  
    ; for each type
    shallowi=where(type eq 0.0)
    mediumi=where(type eq 1)
    deepi = where(type eq 2)
    
    shallow=[month[shallowi]]
    shallowpdf = HISTOGRAM(shallow, BINSIZE=bs, locations=sbin, /nan) 
    sn=total(finite(shallow))
    
    medium=[month[mediumi]]
    mediumpdf = HISTOGRAM(medium, BINSIZE=bs, locations=mbin, /nan)
    mn=total(finite(medium))
    
    deep=[month[deepi]]
    deeppdf = HISTOGRAM(deep, BINSIZE=bs, locations=dbin, /nan)
    dn=total(finite(deep))
    
    ; if there is one or less deep events we need to handle that
    if n_elements(deep) eq 1 then begin
      deeppdf = [0.,(deepi ne -1),0.]
      dbin = [dbin-bs, dbin, dbin+bs]
    endif
 
    xrange=[0,13]
    yrange=[0, max([shallowpdf, mediumpdf, deeppdf])]
  
    print, 'types:', sn, mn, dn
    
    ; set up plot
    xlabels=['x','J','F','M','A','M','J','J','A','S','O','N','D','x']
    cgplot, xrange, yrange, /nodata, $
      xtitle=xtitles[ptri], ytitle=ytitles[ptri], $
      title=bases[ptri]+'('+ string(n,format='(i4)')+' events)', $
      xtickname=xlabels, xtickv=indgen(14), $
      xticks=13, yminor=1,$
      charsize=3
    
    ; for each type
    if sn gt 0 then $
      cgoplot, sbin, shallowpdf, $
        color = colors[ptri], $
        thick = thk, $
        linestyle=styles[0]
    
    if mn gt 0 then $
      cgoplot, mbin, mediumpdf, $
        color = colors[ptri], $
        thick = thk, $
        linestyle=styles[1]
    
    if dn gt 0 then $  
      cgoplot, dbin, deeppdf, $
        color = colors[ptri], $
        thick = thk, $
        linestyle=styles[2]
    
  endforeach
   
  
;  ; add legend!
;  cglegend, tcolors=colors, $
;    titles=bases +'('+ string(n,format='(I3)')+' events)', $
;    length=0.0, $ ; remove legend lines
;    charthick=1.5, $
;    charsize=2., $
;    vspace=2.4, $
;    location = [0.20, 0.88]
  
  ; add another legend!
  cglegend, $;colors=0, $
    titles = ['shallow', 'medium', 'deep'], $
    linestyles=styles, $
    length = 0.05, $
    charsize=2, $
    thick=1.5, $
    location = [0.74, 0.59]
  
!p.multi=0
!y.margin=0
END
