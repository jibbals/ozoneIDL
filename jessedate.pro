; function returning floating year date from a julian array
;

FUNCTION jessedate, julian
  CALDAT, julian[0], M, D, Y
  RETURN, (julian-julian[0]) / 365.25 + Y + (M-1)/12.0 + (D-1)/365.25
END
