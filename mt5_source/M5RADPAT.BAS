sub m5radpat (kgetder,kpwave,strike,dip,rake,az,rp,alph2,beta2,rho2,f(1),df(2),vslow(1)) static

'--- subroutine for radiation patterns and derivatives
'      Input: KGETDER - flag to get derivatives
'             KPWAVE - true for P wave, false for SH
'      Output: F(i) - radiation pattern amplitudes for down and up rays
'              DF(i,j) - derivatives of F(i) w.r.t. parameter j
'              VSLOW(i) - vertical slownesses
'
'   pii=3.141592654
   a=az-strike: sa=sin(a): ca=cos(a): s2a=sin(2.*a): c2a=cos(2.*a): sd=sin(dip)
   sa2=sa*sa
   cd=cos(dip): s2d=sin(2.*dip): c2d=cos(2.*dip): sr=sin(rake): cr=cos(rake)
   rp2=rp*rp: vsb2=1/(beta2*beta2)-rp2 : vsb=sqr(vsb2)

 if kpwave then 
'---- P-SV radiation:  --------------
    vsa2=1/(alph2*alph2)-rp2 : vsa=sqr(vsa2)
    vslow(1)=vsa : vslow(2)=vsb

'---- amplitudes: F(i), i= 1-P down, 2-P up, 3-SV up
      pcor=4.*pii*rho2*alph2: scor=4*pii*rho2*beta2*(vsb/vsa)
      a1=sd*cr*s2a-s2d*sr*sa2: a2=s2d*sr: a3=-cd*cr*ca+c2d*sr*sa
      e1=a1*rp2+a2*vsa2: e2=2*a3*rp*vsa
      e3=(vsb2-rp2)*a3: e4=rp*vsb*(a2-a1)
    f(1)=(e1+e2)/pcor   'P down
    f(2)=(e1-e2)/pcor   'P up
    f(3)=-(e3+e4)/scor  'SV up

'---- derivatives; i indexing 1-P down, 2-P up, 3-SV up
'                  j indexing 1-strike, 2-dip, 3-rake
  if kgetder then
      b3=-cd*cr*sa-c2d*sr*ca: b1=-2.*sd*cr*c2a+s2d*sr*s2a
      e1=b1*rp2: e2=2.*b3*rp*vsa
      e3=(vsb2-rp2)*b3: e4=-b1*rp*vsb
    df(1,1)=(e1+e2)/pcor
    df(2,1)=(e1-e2)/pcor
    df(3,1)=-(e3+e4)/scor  'STRIKE

      b2=2.*c2d*sr: b3=sd*cr*ca-2.*s2d*sr*sa
      b1=cd*cr*s2a-2.*c2d*sr*sa2
      e1=b1*rp2+b2*vsa2: e2=2.*b3*rp*vsa
      e3=(vsb2-rp2)*b3: e4=rp*vsb*(b2-b1)
    df(1,2)=(e1+e2)/pcor
    df(2,2)=(e1-e2)/pcor
    df(3,2)=-(e3+e4)/scor  'DIP

      b2=s2d*cr: b3=cd*sr*ca+c2d*cr*sa
      b1=-sd*sr*s2a-s2d*cr*sa2
      e1=b1*rp2+b2*vsa2: e2=2.*b3*rp*vsa
      e3=(vsb2-rp2)*b3: e4=rp*vsb*(b2-b1)
    df(1,3)=(e1+e2)/pcor
    df(2,3)=(e1-e2)/pcor
    df(3,3)=-(e3+e4)/scor 'RAKE
  end if

else      '---- SH radiation pattern:  ------------
 vslow(1)=vsb

'---- amplitudes for SH.  F(1)-SH down, F(2)-SH up
      scor=4.0*pii*rho2*beta2*beta2
      a2=cd*cr*sa+c2d*sr*ca
      a1=sd*cr*c2a-0.5*s2d*sr*s2a
      e1=a1*rp : e2=a2*vsb
    f(1)=(e1+e2)/scor
    f(2)=(e1-e2)/scor
    f(3)=0


'--- derivatives; i indexing 1-SH down, 2-SH up
'                 j indexing 1-strike, 2-dip, 3-rake
 if kgetder then 
      b1=2*sd*cr*s2a+s2d*sr*c2a
      b2=-cd*cr*ca+c2d*sr*sa
      e1=b1*rp : e2=b2*vsb
    df(1,1)=(e1+e2)/scor
    df(2,1)=(e1-e2)/scor
    df(3,1)=0                 'STRIKE
      b1=cd*cr*c2a-c2d*sr*s2a 
      b2=-sd*cr*sa-2*s2d*sr*ca
      e1=b1*rp : e2=b2*vsb
    df(1,2)=(e1+e2)/scor
    df(2,2)=(e1-e2)/scor
    df(3,2)=0                 'DIP
      b1=-sd*sr*c2a-.5*s2d*sr*s2a
      b2=-cd*sr*sa+c2d*cr*ca
     e1=b1*rp :e2=b2*vsb
   df(1,3)=(e1+e2)/scor
   df(2,3)=(e1-e2)/scor
   df(3,3)=0                  'RAKE
 end if

end if

end sub

