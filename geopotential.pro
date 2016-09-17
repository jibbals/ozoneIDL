
PRO geopotential, file=file

  if n_elements(file) eq 0 then $
    file="../Data/Davis/ERA/Pressure/pressure_20110101-20110430.nc"
  sa_read_era_interim, pdata, FILE=file

  ; plot out the geopotential for our three pressures
  ; 
  
  X = pdata.longitude   ; lons
  Y = pdata.latitude    ; lats
  Z = MEAN(pdata.geopotential, 4, /NAN)*1.0e-3 ; [lons,lats,levels,time] in kms
  loc= locfromstr(file)
  titles = 'GPH(kms) at '+['200','300','500'] +'hPa over ' + loc
  
  ; colortable script from AAD
  LOADCT, 33
  
  ; note: can't do loops in a script, it has to be in a procedure or function
  FOR J = 0, 2 DO BEGIN
    cgdisplay, 600,500, wid=J
    ; Set up the map for our contour plots
    cgMAP_SET, -50, 120, /ROBINSON, LIMIT=[-85, 40, -40, 160], title=titles[J]
    ; plot the filled in contour map
    cgcontour, Z[*,*,J], X, Y, nlev=100, /CELL, /OVER
    ; add in some labeled contours
    cgcontour, Z[*,*,J], X, Y, /OVER, nlev=7, /FOLLOW
    ; add the continents
    cgMAP_CONTINENTS
  ENDFOR
  
lEND
