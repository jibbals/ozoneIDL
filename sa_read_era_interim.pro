;+
;
;  Name                 SA_READ_ERA_INTERIM.pro
;
;  Purpose              Read in NetCDF ECMWF data from Interim Re-analysis
;
;  Inputs               None
;
;  Keywords             FILE - Full filepath of the NetCDF file to read in
;
;  Optional
;   Keywords
;                       TIME_RANGE - Double 2-element array of start & end Jultime
;                       
;                       LON_RANGE - Double 2-element array of start & end longitudes
;
;                       LAT_RANGE - Double 2-element array of start & end latitudes
;
;                       PRESSURE_RANGE - Double 2-element array of start & end pressure
;
;
;
;-

PRO SA_READ_ERA_INTERIM,ERA_Interim,FILE=FILE $
       ; Optional - if reading part of the array
     ,LON_RANGE=LON_RANGE,LAT_RANGE=LAT_RANGE $
     ,TIME_RANGE=TIME_RANGE,PRESSURE_RANGE=PRESSURE_RANGE

    ; Missing Value
    dmiss = !values.f_nan;999.

    ; File Open
    fd=NCDF_OPEN(FILE,/NOWRITE)

    ; Get Variable ID
    lonid=NCDF_VARID(fd,'longitude')
    latid=NCDF_VARID(fd,'latitude')
    timeid=NCDF_VARID(fd,'time')
    ;levelid=NCDF_VARID(fd,'levelist')
    levelid=NCDF_VARID(fd,'level')
    IF levelid EQ -1 THEN levelid=NCDF_VARID(fd,'levelist')  ; try different version of ERA data

    ; Read Data
    NCDF_VARGET,fd,lonid,longitude
    NCDF_VARGET,fd,latid,latitude
    NCDF_VARGET,fd,timeid,time
    IF levelid NE -1 THEN NCDF_VARGET,fd,levelid,pressure

    ; Time is hours since 1900-01-01 00:00:00
    start=sa_conv_time([0,0,0,1,1,1900],/to)
    jultime=start+time/24d

    ; May only have some of these variables
    Tid=NCDF_VARID(fd,'t') & O3id=NCDF_VARID(fd,'o3')
    Zid=NCDF_VARID(fd,'z') & void=NCDF_VARID(fd,'vo')
    PVid=NCDF_VARID(fd,'pv') & Qid=NCDF_VARID(fd,'q')
    Wid=NCDF_VARID(fd,'w') & Rid=NCDF_VARID(fd,'r')
    CLWCid=NCDF_VARID(fd,'clwc') & CIWCid=NCDF_VARID(fd,'ciwc')
    Uid=NCDF_VARID(fd,'u') & Vid=NCDF_VARID(fd,'v')
    Did=NCDF_VARID(fd,'d')
    params=[Tid,O3id,Zid,Void,PVid,Qid,Wid,Rid,CLWCid,CIWCid,Uid,Vid,Did]
    n_possible_params=N_ELEMENTS(params)

    ; Zero them initially for putting in structure
    temperature=0  & ozone_mmr=0 & geopot=0  &  zeta=0
    potential_vorticity=0 & specific_humidity=0
    vertical_velocity=0 & relative_humidity=0
    cloud_LWC=0 & cloud_IWC=0
    zonal_velocity=0 & meridional_velocity=0 & divergence=0

    ; Defaults all data if only setting some of below 
    lon_idx=WHERE(longitude GE -180. AND longitude LE 360.)
    time_idx=WHERE(jultime GT 0 AND FINITE(jultime) EQ 1)
    IF KEYWORD_SET(pressure) THEN pressure_idx=WHERE(pressure LT 1020.)
    lat_idx=WHERE(latitude GE -90. AND latitude LE 90.)

    ; Set Data Range - See OLR_REAR_EAR.PRO
    IF KEYWORD_SET(LON_RANGE) THEN $
      lon_idx = WHERE(longitude ge lon_range(0) and longitude le lon_range(1))
    IF KEYWORD_SET(LAT_RANGE) THEN $
      lat_idx = WHERE(latitude ge lat_range(0) and latitude le lat_range(1))  
    IF KEYWORD_SET(TIME_RANGE) THEN $
      time_idx = WHERE(jultime ge time_range(0) and $
                       jultime le time_range(1))
    IF KEYWORD_SET(PRESSURE_RANGE) THEN $
      pressure_idx = WHERE(pressure GE MIN(pressure_range) and $
                           pressure LE MAX(pressure_range))
 
     IF KEYWORD_SET(LON_RANGE) OR KEYWORD_SET(LAT_RANGE) OR KEYWORD_SET(TIME_RANGE) THEN $
      read_part=1 $
    ELSE $
      read_part=0

  dummy=0.
  FOR i=0,n_possible_params-1 DO $
    BEGIN
    pp=params[i]
    IF pp NE -1 THEN $
      BEGIN
        IF KEYWORD_SET(read_part) THEN $
          BEGIN
            IF levelid NE -1 THEN $
              NCDF_VARGET,fd,pp,dummy,$
                count =[N_ELEMENTS(lon_idx),N_ELEMENTS(lat_idx),N_ELEMENTS(pressure_idx),N_ELEMENTS(time_idx)],$
                OFFSET=[lon_idx(0),lat_idx(0),pressure_idx[0],time_idx(0)] $
            ELSE $
              NCDF_VARGET,fd,pp,dummy,$
                count =[N_ELEMENTS(lon_idx),N_ELEMENTS(lat_idx),N_ELEMENTS(time_idx)],$
                OFFSET=[lon_idx(0),lat_idx(0),time_idx(0)]
          ENDIF ELSE $
            NCDF_VARGET,fd,pp,dummy
        dummy=FLOAT(dummy)
        NCDF_ATTGET,fd,pp,'scale_factor',scale
        NCDF_ATTGET,fd,pp,'add_offset',add_offset
        NCDF_ATTGET,fd,pp,'missing_value',missing_value
      ; To save space, comment out these 2 lines
        idx = where(dummy ne missing_value)
        IF idx(0) NE -1 THEN dummy(idx) = FLOAT(dummy(idx))*scale + add_offset
        idx = WHERE(dummy EQ missing_value)
        IF idx(0) NE -1 THEN dummy(idx) = dmiss
    ENDIF ELSE dummy=0
    IF i EQ 0 THEN temperature=dummy
    IF i EQ 1 THEN ozone_mmr=dummy
    IF i EQ 2 THEN geopot=dummy/9.80665  ; divide by gravity to get into metres (see ERA website FAQ)
    IF i EQ 3 THEN zeta=dummy
    IF i EQ 4 THEN potential_vorticity=dummy
    IF i EQ 5 THEN specific_humidity=dummy
    IF i EQ 6 THEN vertical_velocity=dummy
    IF i EQ 7 THEN relative_humidity=dummy
    IF i EQ 8 THEN cloud_LWC=dummy
    IF i EQ 9 THEN cloud_IWC=dummy
    IF i EQ 10 THEN zonal_velocity=dummy
    IF i EQ 11 THEN meridional_velocity=dummy
    IF i EQ 12 THEN divergence=dummy
  ENDFOR

  IF KEYWORD_SET(read_part) THEN $
    BEGIN
      longitude=longitude(lon_idx)
      latitude=latitude(lat_idx)
      jultime=jultime(time_idx)
      IF KEYWORD_SET(pressure) THEN pressure=pressure[pressure_idx]
  ENDIF


  IF levelid EQ -1 THEN pressure=-1  ; dummy - probably only 1 pressure altitude to read in

  ERA_Interim={ $
        ; Always include the basic parameters
        jultime:jultime,longitude:longitude $
        ,latitude:latitude,pressure:pressure $
        ; Define all possible variables: zero if not in file
        ,temperature:temperature,ozone_mmr:ozone_mmr $
        ,geopotential:geopot $  
        ,zeta:zeta $
        ,potential_vorticity:potential_vorticity $
        ,specific_humidity:specific_humidity $
        ,vertical_velocity:vertical_velocity $
        ,relative_humidity:relative_humidity $
        ,cloud_LWC:cloud_LWC,cloud_IWC:cloud_IWC $
        ,zonal_velocity:zonal_velocity $
        ,meridional_velocity:meridional_velocity $
        ,divergence:divergence}

  RETURN

END