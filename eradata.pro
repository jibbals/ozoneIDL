;   Function: 
;           eradata
;   
;   Purpose:
;           Read the ERA netcdf files for specific dates
;   
;   Parameters:
;           dates: array of julian dates to extract from the ERA datafiles
;           
;   Keywords:
;           Melbourne/Davis/Macquarie: set one of these(Melbourne default)
;           HIGHRES: set to look only at the high resolution ERA data
;   
;   Requirements:
;           sa_read_era_interim: Dr Alexander's era reading script
FUNCTION eradata, dates, $
  MELBOURNE=MELBOURNE, DAVIS=DAVIS, MACQUARIE=MACQUARIE, $
  HIGHRES=HIGHRES, $
  FILE=FILE

  ; shortcut file:
  if n_elements(FILE) gt 0 then begin
    restore, file=file
    return, data
  endif

  filelen=29
  n_pressures=4
  if keyword_set(HIGHRES) then begin
    filelen=20
    n_pressures=34
  endif

  if keyword_set(davis) then begin
    path='../Data/Davis/ERA/Pressure/'
  endif else $  
  if keyword_set(macquarie) then begin
    path='../Data/Macquarie/ERA/Pressure/'
  endif else begin   ; melbourne default
  ;if keyword_set(melbourne) then begin
    path='../Data/Broadmeadows/ERA/Pressure/'
  endelse
  
  ; just get the files which are synoptic weather format
  ; eg: 'pressure_20050101-20050430.nc'
  files=file_search(path,'*.nc')
  
  files=files[where(strlen(file_basename(files)) eq filelen)]
  
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
     
    ; highres files will not be uniform
    if keyword_set(highres) then begin
        N_lons=n_elements(era.longitude)
        N_lats=n_elements(era.latitude)
        N_levs=n_elements(era.pressure)
        longitude=era.longitude
        latitude=era.latitude
        pressure=era.pressure
        temp = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
        o3vmr = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
        gpot = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
        potv = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
        rh = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
        vvel = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
        zonalvel = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
        meridvel = fltarr([N_lons,N_lats,N_levs,N_times])+!values.f_nan
    endif
     
    if n_elements(era.pressure) ne n_pressures then begin
        print, 'WARNING: non standard levels'
        caldat, jdays[0], cmon,cday,cyear
        print, 'Skipping file matching date: ', cyear, cmon, cday
        continue
    endif
     
    ;grab the weather maps closest to the dates
    foreach i, getind, gi do begin
      tmp=min(abs(jdays-dates2[i]), imatch)
      if tmp gt 1 then $
        print, 'WARNING: closest date is ',tmp,' days from sonde event'
      
      
      ; skip if we are more than three days from a match
      if tmp gt 3 then continue
      
      jultime[savind] = jdays[imatch]
      temp[*,*,*,savind] =  era.temperature[*,*,*,imatch]
      o3vmr[*,*,*,savind] = o3mmrtovmr(era.ozone_mmr[*,*,*,imatch])
      gpot[*,*,*,savind] =  era.geopotential[*,*,*,imatch]
      potv[*,*,*,savind] =  era.potential_vorticity[*,*,*,imatch]
      rh[*,*,*,savind] = era.relative_humidity[*,*,*,imatch]
;      if n_elements(era.specific_humidity) gt 1 then $
;        shumid[*,*,*,savind] = era.specific_humidity[*,*,*,imatch]
      vvel[*,*,*,savind] =  era.vertical_velocity[*,*,*,imatch]
      zonalvel[*,*,*,savind] = era.zonal_velocity[*,*,*,imatch]
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
