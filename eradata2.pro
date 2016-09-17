;   Function: 
;           eradata2
;   
;   Purpose:
;           Read the ERA netcdf files for specific dates
;           Uses updated ERAI files
;   
;   Parameters:
;           dates: array of julian dates to extract from the ERA datafiles
;           
;   Keywords:
;           Melbourne/Davis/Macquarie: set one of these(Melbourne default)
;   
;   Requirements:
;           Knife and fork for dinner.
;           sa_read_era_interim.pro - jesselib/ozone/...

FUNCTION eradata2, dates, $
  MELBOURNE=MELBOURNE, DAVIS=DAVIS, MACQUARIE=MACQUARIE, $
  FILE=FILE

  ; shortcut file:
  if n_elements(FILE) gt 0 then begin
    restore, file=file
    return, data
  endif

  caldat, dates,mnths,dys,yrs
  ystr=string(yrs[0],'(i4)')
    
  ; new datafiles look like this 
  ;OzoneWork/Data/ERAI/ERAI_2004_pl.nc 
  path='../Data/ERAI/'
  ptrn='ERAI_*_pl.nc'
  ; if only looking for one day, match the year
  if n_elements(dates) eq 1 then ptrn='ERAI_'+ystr+'_pl.nc'
  files=file_search(path,ptrn)
  sa_read_era_interim, tmp, file=files[0]
  N_times=n_elements(dates)
  N_lons=n_elements(tmp.longitude)
  N_lats=n_elements(tmp.latitude)
  N_levs=n_elements(tmp.pressure)
  
  ; create the arrays we will pass back to caller
  ;
  jultime=fltarr(N_times) + !values.f_nan
  longitude=tmp.longitude
  latitude=tmp.latitude
  pressure=tmp.pressure
  temp = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
  o3vmr = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
  gpot = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
  potv = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
  rh = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
  vvel = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
  zonalvel = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
  meridvel = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
  
  savind=0
  dates2=dates[sort(dates)]
  
  ; for each file
  foreach file, files, fi do begin
    sa_read_era_interim, era, file=file
    
    jdays=era.jultime
    getind=where(dates2 ge jdays[0] and dates2 le jdays[-1])
    if getind[0] eq -1 then continue
     
    ;grab the weather maps closest to the dates
    foreach i, getind, gi do begin
      tmp=min(abs(jdays-dates2[i]), imatch)
      if tmp gt 1 then $
        print, 'WARNING: closest date is ',tmp,' days from sonde event'
      
      ; skip if we are more than three days from a match
      if tmp gt 3 then continue
      
      jultime[savind] = jdays[imatch]
      if n_elements(era.temperature) gt 1 then $
        temp[*,*,*,savind] =  era.temperature[*,*,*,imatch]
      if n_elements(era.ozone_mmr) gt 1 then $
        o3vmr[*,*,*,savind] = o3mmrtovmr(era.ozone_mmr[*,*,*,imatch])
      if n_elements(era.geopotential) gt 1 then $
        gpot[*,*,*,savind] =  era.geopotential[*,*,*,imatch]
      if n_elements(era.potential_vorticity) gt 1 then $
        potv[*,*,*,savind] =  era.potential_vorticity[*,*,*,imatch]
      if n_elements(era.relative_humidity) gt 1 then $
        rh[*,*,*,savind] = era.relative_humidity[*,*,*,imatch]
      if n_elements(era.specific_humidity) gt 1 then $
        shumid[*,*,*,savind] = era.specific_humidity[*,*,*,imatch]
      if n_elements(era.vertical_velocity) gt 1 then $
        vvel[*,*,*,savind] =  era.vertical_velocity[*,*,*,imatch]
      if n_elements(era.zonal_velocity) gt 1 then $
        zonalvel[*,*,*,savind] = era.zonal_velocity[*,*,*,imatch]
      if n_elements(era.meridional_velocity) gt 1 then $
        meridvel[*,*,*,savind] = era.meridional_velocity[*,*,*,imatch]
      savind = savind + 1
    endforeach
    ; leave early if we have all the dates
    if savind eq n_elements(dates) then break
  endforeach
  
  data={jultime:jultime,longitude:longitude,latitude:latitude, $
    pressure:pressure, temperature:temp, ozone_mmr:o3vmr, $
    geopotential:gpot, potential_vorticity:potv, $
    relative_humidity:rh, $
    vertical_velocity:vvel, $
    zonal_velocity:zonalvel, meridional_velocity:meridvel }

  print, long(total(finite(jultime))), ' date matches for ', n_times, ' events' 
  
  return, data
END
