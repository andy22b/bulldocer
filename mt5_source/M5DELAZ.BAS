sub m5delaz(slat,slon,eplat,eplon,delta,azim) static
'------------ get delta and azimuth
  F=1.745329E-02 
 ' --- distaz
        TSLAT = slat*F
        TSLON = slon*F
  GLAT=ATN(.993277*TAN(TSLAT))
        DCA = COS(GLAT)*COS(TSLON)
        DCB = COS(GLAT)*SIN(TSLON)
        DCC = SIN(GLAT)
        SLAT = TSLAT/F
  GEPLAT = ATN(.993277*TAN(EPLAT*F))
  EPDCA = COS(GEPLAT)*COS(EPLON*F)
  EPDCB = COS(GEPLAT)*SIN(EPLON*F)
  EPDCC = SIN(GEPLAT)
 'calculate distance and azimuth from earthquake to station
  PAR1=(1-.5*((DCA -EPDCA)^2 + (DCB -EPDCB)^2+(DCC  - EPDCC)^2))
        DELTA  = 57.29578*ATN(SQR((1/(PAR1*PAR1)) - 1))
        IF (PAR1<0) THEN DELTA  = 180 - DELTA
        IF(DELTA=0) THEN DELTA = .001
        T1 = DCA  - SIN(EPLON*F)
        T2 = DCB  + COS(EPLON*F)
        T3 = SIN(DELTA *F)
        PAR2 = (.5*(T1*T1 + T2*T2 + DCC *DCC ) - 1)/T3
        T1 = DCA  - SIN(GEPLAT)*COS(EPLON*F)
        T2 = DCB  - SIN(GEPLAT)*SIN(EPLON*F)
        T4 = DCC  + COS(GEPLAT)
        APAR2 = ABS(PAR2)
        IF (APAR2 -1>=0) THEN azim  = 90 : GOTO daz1
        azim  = 57.29578*ATN(1/SQR((1/PAR2^2)-1))
daz1:   PAR3 = (.5*(T1*T1+T2*T2+T4*T4) - 1)/T3
        IF (PAR2>=0) THEN goto daz2
      IF (PAR3<0) THEN azim  = 180 + azim :  exit sub
        azim  = 360 - azim :  exit sub
daz2:  IF(PAR3<0) THEN azim  = 180 - azim
end sub

