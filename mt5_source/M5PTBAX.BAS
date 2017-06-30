sub m5ptbax(str,dip,slip,paz,ppl,taz,tpl,baz,bpl) static
  dim xp(2),yp(2),zp(2),pax(2),sax(2)
  const dpr=57.2957795

 call m5auxpl(str,dip,slip,str2,dip2,slip2)  'get auxiliary plane
    ' get P and T axes
    xp(1)=sin(str)*sin(dip)+sin(str2)*sin(dip2)
    xp(2)=sin(str)*sin(dip)-sin(str2)*sin(dip2)
    yp(1)=-cos(str)*sin(dip)-cos(str2)*sin(dip2)
    yp(2)=-cos(str)*sin(dip)+cos(str2)*sin(dip2)
    zp(1)=cos(dip)+cos(dip2)
    zp(2)=cos(dip)-cos(dip2)

  for i=1 to 2
    pax(i)=fnatan2(zp(i),sqr(xp(i)^2+yp(i)^2))
    sax(i)=fnatan2(yp(i),xp(i))
    if pax(i)<0 then pax(i)=-pax(i): sax(i)=sax(i)+pii
    if pax(i)>piihalf then pax(i)=pii-pax(i): sax(i)=sax(i)+pii
    s=sax(i)*dpr
    if not fninlim(s,0,360) then s= s mod 360
    if s<0 then s=s+360
    sax(i)=s/dpr                  'put strike between 0 and 360
  next

  paz=sax(2) : ppl=pax(2) : taz=sax(1) : tpl=pax(1)

  if not fninlim(slip,0,pii) then 
    swap paz,taz :  swap ppl,tpl
  end if 

'--- get b axis

    p1=paz : p2=ppl
    py=cos(p2)*cos(p1) 'X,Y,Z COORDS
    px=cos(p2)*sin(p1) '+X=e, +Y=n, +Z=DOWN
    pz=sin(p2)

    p1=taz : p2=tpl
    ty=cos(p2)*cos(p1) 
    tx=cos(p2)*sin(p1) 
    tz=sin(p2)

'-- cross product
    bx=py*tz-pz*ty
    by=pz*tx-px*tz
    bz=px*ty-py*tx

    baz=fnatan2(bx,by)
    bpl=atn(bz/sqr(bx*bx+by*by))
    if bpl>piihalf then bpl=pii-bpl : baz=baz+pii
    if bpl<0. then bpl=-bpl : baz=baz+pii
    bazd=dpr*baz+720 : baz=(bazd mod 360)/dpr

end sub

