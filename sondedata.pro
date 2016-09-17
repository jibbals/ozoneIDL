; read the sonde dataset
FUNCTION SondeData, $
  davis=davis, macquarie=macquarie, melbourne=melbourne
  
  ; Get the stored ozonesonde data(compiled from csvs)
  ;
  if keyword_set(davis) then begin
    davisflag=davis
    file='/home/jwg366/OzoneWork/Data/Davis/Sondes/davsondes.dat'
    restore, file
    data=davis
  endif else $  
  if keyword_set(macquarie) then begin
    macquarieflag=macquarie
    file='/home/jwg366/OzoneWork/Data/Macquarie/Sondes/macsondes.dat'
    restore, file
    data=Macca
    ; one anomolously high ozone day seems wrong, unless it was released into a fire plume?
    data.o3pp[182,*] = data.o3pp[182,*]+!values.f_nan
    data.pressure[182,*] = data.pressure[182,*]+!values.f_nan
  endif else begin   ; melbourne default
    melbourneflag=1
    file='/home/jwg366/OzoneWork/Data/Broadmeadows/Sondes/melbsondes.dat'
    restore, file
    data=Melb
  endelse
  
  ; Add some usefull stuff to sonde data structure:
  ;
  
  ; vmr is PP(mPa) / P(hPa) * 1/100 000
  vmr=data.o3pp / data.pressure * 1.0e-5
  ppbv=vmr*1.0e9
  
  ; tropopause height
  n=n_elements(data.jtime)
  tpo3=fltarr(n) + !values.f_nan
  tplr=fltarr(n) + !values.f_nan
  tpo3p=fltarr(n) + !values.f_nan
  tplrp=fltarr(n) + !values.f_nan
  ; loop through non missing columns
  gphinds=where(data.gph[*,0] gt 0.0)
  foreach i, gphinds do begin
    tpo3[i]=OzoneTropopause(reform(ppbv[i,*]), $
      reform(data.gph[i,*]/1000.0), POLAR=DAVISFLAG)
    tplr[i]=TemperatureTropopause(reform(data.temperature[i,*]), reform(data.gph[i,*]/1000.0))
    ; get pressures as well as km heights for our tropopauses
    o3ind = where(data.gph[i,*]/1000.0 eq tpo3[i], count)
    if count eq 1 then tpo3p[i] = data.pressure[i, o3ind]
    lrind = where(data.gph[i,*]/1000.0 eq tplr[i], count)
    if count eq 1 then tplrp[i] = data.pressure[i, lrind]
  end  
  
  ; list of transport caused events (using AIRS analysis)
  

  ; idlutils function to add tags
  struct={o3ppbv:ppbv, tpo3:tpo3, tplr:tplr, tpo3p:tpo3p, tplrp:tplrp}; , tppvu:tppvu}
  datanew=struct_addtags(data, struct)
  
  
  return, datanew
END
