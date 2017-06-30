sub m5auxpl (str, dip, slip, str2, dip2, slip2) static
'---- get P-wave auxiliary plane
'  shared pii
     aa=cos(slip)*cos(str)+cos(dip)*sin(slip)*sin(str)
     ba=cos(slip)*sin(str)-cos(dip)*sin(slip)*cos(str)
     ca=sin(slip)*sin(dip)
     xlen=sqr(aa*aa+ba*ba+ca*ca)
     dip2=fnacos(ca/xlen)
     str2=fnatan2(ba,aa)
     str2=str2-piihalf
aux1: if str2<0 then str2=str2+2*pii
     if dip2>=0 and dip2 <=piihalf then goto aux2
     if dip2<0 then dip2=dip2+pii
     if dip2>piihalf then dip2=fnamin(piihalf,pii-dip2): str2=str2-pii: goto aux1
aux2: '--- get slip direction
      az1=str-piihalf : pl1=piihalf-dip   'normal to fault plane
      az2=str2 : pl2=0                'strike of aux plane

'-- does not return slip angle 2
slip2=0

' slip2=az2-az1
'      x1=cos(az1) : y1=sin(az1) : z1=sin(pl1)
'      x2=cos(az2) : y2=sin(az2) : z2=sin(pl2)
'    xm = sqr(x1*x1+y1*y1+z1*z1) * sqr(x2*x2+y2*y2+z2*z2)
'    xprod= (x1*x2 + y1*y2 + z1*z2)
'    t= fnacos (xprod/xm)
'    slip2=t
'    if slip>pii then slip2=pii+slip2

end sub

