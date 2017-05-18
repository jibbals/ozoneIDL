;+
;
;  Name                    sa_conv_time.pro
;
;  Purpose                 Convert [s,min,hr,day,month, year]
;							to Julian day format or vice versa
;
;  Inputs                  in - Julian date
;
;  Output				   jtime - in Julian seconds
;
;
;  Return Value:           Julian time in seconds if /to_jul set,
;							else a six element array of form:
;							[second, minute, hour, day, month, year]
;
;  Notes:                  Note that the IDL Julday procedure is formatted
; 							[month,day,year,hour,min,sec] which is
;							inconvenient... Hence this program...
;-

FUNCTION SA_CONV_TIME,in,to_jul=to_jul


    IF KEYWORD_SET(to_jul) THEN $
      BEGIN
        ; going from [s,min,hr,day,month, year] to Jultime
        input=lonarr(6)
        input(0)=in(4)  &  input(1)=in(3)  &  input(2)=in(5)
        input(3)=in(2)  &  input(4)=in(1)  &  input(5)=in(0)
        out=JULDAY(input(0),input(1),input(2),input(3),input(4),input(5))
        RETURN,out
    ENDIF ELSE $
      BEGIN
        ; going from Jultime to [s,min,hr,day,month, year]
        CALDAT,in,mon,day,year,hr,minu,sec
        out=[sec,minu,hr,day,mon,year]
        RETURN,LONG(out)
    ENDELSE

    RETURN,-1

END
