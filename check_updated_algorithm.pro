function jtime_to_str, jtimes
  caldat, jtimes, mm,dd,yy
  list=strarr(n_elements(mm))
  for j=0,n_elements(mm)-1 do $
        list[j]=string(yy[j],format='(i04)')+string(mm[j],format='(i02)')+string(dd[j],format='(i02)')
  return, list
end

pro check_updated_algorithm, plot_profiles=plot_profiles

  newdav=ptr_new(getevents(sondedata(/dav), TROP_DEF=0))
  newmac=ptr_new(getevents(sondedata(/mac), TROP_DEF=0))
  newmelb=ptr_new(getevents(sondedata(/melb), TROP_DEF=0))
  olddav=ptr_new(getevents(sondedata(/dav), TROP_DEF=2))
  oldmac=ptr_new(getevents(sondedata(/mac), TROP_DEF=2))
  oldmelb=ptr_new(getevents(sondedata(/melb), TROP_DEF=2))

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
    newdatelist=strarr(n_elements(mm))
    caldat, removals, mmr,ddr,yyr
    removedlist=strarr(n_elements(mmr))
    for j=0,n_elements(mmr)-1 do $
        removedlist[j]=string(yyr[j],format='(i04)')+string(mmr[j],format='(i02)')+string(ddr[j],format='(i02)')
    for j=0,n_elements(mm)-1 do $
        newdatelist[j]=string(yy[j],format='(i04)')+string(mm[j],format='(i02)')+string(dd[j],format='(i02)')
    ;print, string(newdates,format='(f9.1)')
    print, 'new dates:', newdatelist
    print, 'dates removed:', removedlist

    if keyword_set(plot_profiles) then begin
      olddates=jtime_to_str(old.jtime)
      newdates=jtime_to_str(new.jtime)
      locstr=(['Dav','Mac','Melb'])[i]
      for ii=0,n_elements(old.jtime)-1 do begin
        name='prof_'+olddates[ii]
        cgps_open, 'images/profiles/old/'+locstr+'/'+name+'.png', $
            XSize=12, YSize=18, inches=0
        profilepressure,reform(old.o3ppbv[ii,*]),reform(old.pressure[ii,*]),$
        tpo3=old.tpo3[ii], humid=reform(old.rh[ii]), top=200, $
        title=locstr+name,tplr=old.tplr[ii],temp=reform(old.temperature[ii])
        cgps_close, /png
      endfor
      for ii=0,n_elements(new.jtime)-1 do begin
        name='prof_'+newdates[ii]
        cgps_open, 'images/profiles/new/'+locstr+'/'+name+'.png', $
            XSize=12, YSize=18, inches=0
        profilepressure,reform(new.o3ppbv[ii,*]),reform(new.pressure[ii,*]),$
        tpo3=new.tpo3[ii], humid=reform(new.rh[ii]), top=200,$
        title=locstr+name,tplr=new.tplr[ii],temp=reform(new.temperature[ii])
        cgps_close, /png
      endfor
    endif
  endfor
end
