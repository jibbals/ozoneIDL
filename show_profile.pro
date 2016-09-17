;Procedure: show_profile
;
;Purpose: create some plots looking at sonde release dataset
;
;Inputs:
;
;Returns:
;   
;Outputs:
;   
;KEYWORDS:
;    jday=day you want profiled(julian), eg julday(2,3,2005) ;= 3/feb/2005
;    top= height you want to plot up to.
;    cgps= turn on if you want to store output using cgps_open...cgps_close
;NOTES:
;    Outputs to a .ps with !P.font=0 so that text objects are created for inkscape editting

PRO show_profile, jday=jday, $
  melbourne=melbourne, davis=davis, macquarie=macquarie, $
  top=top, cgps=cgps

  if keyword_set(melbourne) then $
    title="Melbourne"
  if keyword_set(davis) then $
    title="Davis"
  if keyword_set(macquarie) then $
    title="Macquarie"

  ; read in the data
  data=sondedata(melbourne=melbourne, davis=davis, macquarie=macquarie)
    
  dates=JesseDate(data.jtime)
  
  press=data.pressure
  ppbv=data.o3ppbv

  n=n_elements(dates)
  gph=data.gph/1000.0
  
  tmp=min(abs(data.jtime - jday), jind)
  
  tps= [data.tpo3p[jind], data.tplrp[jind]] 
  
  ;==================================================
  ;             Plot profile
  ;==================================================
  if keyword_set(cgps) then begin
    !P.font=0
    CGPS_OPEN, title+'profile.ps'
  endif
  
  ; date stuff  
  caldat, data.jtime[jind], M, D, Y, H
  fi2='(i02)'
  dmy = string(D,format=fi2)+'/'+string(M,format=fi2)+'/'+string(Y MOD 100, format=fi2)
  title=title+' ' + dmy
  
  cgdisplay, 500, 800, wid=0, title='Profile'+string(dates[jind])
  
  temp=reform(data.temperature[jind,*])+273.15
  RH = reform(data.RH[jind,*]) 
    
  profilepressure, reform(ppbv[jind,*]), reform(press[jind,*]), tpo3=tps[0], title=title, $
    top=top, temp=temp, Humid=RH, tplr=tps[1]

  if keyword_set(cgps) then CGPS_CLOSE
end 
