pro check_updated_algorithm

  newdav=ptr_new(getevents(sondedata(/dav)))
  newmac=ptr_new(getevents(sondedata(/mac)))
  newmelb=ptr_new(getevents(sondedata(/melb)))
  olddav=ptr_new(old_getevents(sondedata(/dav)))
  oldmac=ptr_new(old_getevents(sondedata(/mac)))
  oldmelb=ptr_new(old_getevents(sondedata(/melb)))

  oldptrs= [ olddav, oldmac, oldmelb ]
  newptrs= [ newdav, newmac, newmelb ]
  sites=['Davis','Macquarie','Melbourne']
  for i=0,2 do begin
    old=*(oldptrs[i])
    new=*(newptrs[i])
    oldint=long(10*old.jtime)
    newint=long(10*new.jtime)
    print, '================================='
    print, sites[i], ' old: ', n_elements(oldint), ' new: ',n_elements(newint)
    matches=cgsetintersection(oldint, newint)
    newdates= cgsetdifference(newint, matches) * 0.1
    removals= cgsetdifference(oldint, matches) * 0.1
    caldat, newdates, mm,dd,yy
    datelist=strarr(n_elements(mm))
    for j=0,n_elements(mm)-1 do $
        datelist[j]=string(yy[j],format='(i04)')+string(mm[j],format='(i02)')+string(dd[j],format='(i02)')
    ;print, string(newdates,format='(f9.1)')
    print, 'new dates:', datelist
    print, 'dates removed:', string(removals,format='(f9.1)')
  endfor
end
