;   Procedure: 
;       get_all_maps
;
;   Purpose:
;       Look at the wind and gph map for each event
;       
;   Example:
;       get_all_maps
;
;   Prerequisites:
;       coyotegraphics
;       eradata script
;       makemap script
;       sondedata, getevents scripts

PRO get_all_maps, after2010=after2010, parts=parts
  
  if n_elements(parts) eq 0 then parts=0
  ; Read each site
  sites=ptrarr(3)
  sites[0]=ptr_new(sondedata(/davis))
  sites[1]=ptr_new(sondedata(/macqu))
  sites[2]=ptr_new(sondedata(/melbo))
  titles= ['Davis','Macquarie','Melbourne']
  ; Save each map into the map folder
  foreach siteptr, sites, sind do begin
    print, "Running for ", titles[sind]
    ; grab site data, and determine events for that site
    site=*siteptr
    events=getevents(site)
    
    ; for each event
    foreach jday, events.jtime, evind do begin
      if parts eq 1 and evind gt 20 then continue
      if parts eq 2 and (evind lt 21 or evind gt 40) then continue
      if parts eq 3 and (evind lt 41 or evind gt 60) then continue
      if parts eq 4 and (evind lt 61 or evind gt 80) then continue
      if parts eq 5 and (evind lt 81 or evind gt 100) then continue
      if parts eq 6 and evind lt 101 then continue
      
      caldat, jday, cmon, cday, cyr
      ; skip everything before 2010 ( save time, already plotted )
      if keyword_set(after2010) and cyr lt 2010 then continue
      datestr=string(cyr,format='(I4)')+$
        string(cmon,format='(I02)')+$
        string(cday,format='(I02)')
      print, "mapping ",datestr
      
      ; map the era data
      case sind of
        0:makemap, jday, imageprefix="images/maps/Davis/", /davis
        1:makemap, jday, imageprefix="images/maps/Macquarie/", /macquarie
        2:makemap, jday, imageprefix="images/maps/Melbourne/", /melbourne
      endcase
      ; check memory
      mem = MEMORY(/CURRENT)
      ; Prepare dialog message:
      message = 'Current amount of dynamic memory used is '
      print, message + STRTRIM(mem,2)+' bytes.'
    endforeach
  endforeach
  print, "Finished get_all_maps, part="+string(part,format='(i1)')
END
