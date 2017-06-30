sub m5radmt (kgetder,kpwave,amnn,amee,amdd,amne,amnd,amed,az,rp,alph2,beta2,rho2,f(1),df(2),vslow(1)) static

'--- subroutine for radiation patterns and derivatives
'      Input: KGETDER - flag to get derivatives
'             KPWAVE - true for P wave, false for SH
'      Output: F(i) - radiation pattern amplitudes for down and up rays
'              DF(i,j) - derivatives of F(i) w.r.t. parameter j
'                Units of f() and df() are sec^3/g*1.0e-15 for
'                amnn-amed elements as dimensionless normalized moments
'              VSLOW(i) - vertical slownesses
'      Note that SVup amplitude and derivatives have been sign reversed to be 
'       consistent with Aki and Richards designation of + in horizontal
'       direction of ray path

   if kpwave then 
     toang=fnasin(rp*alph2) 
   else 
     toang=fnasin(rp*beta2)
   end if
   sa=sin(az): ca=cos(az): s2a=sin(2.*az): c2a=cos(2.*az) 
   sa2=sa*sa : ca2=ca*ca
   st=sin(toang) : ct=cos(toang): s2t=sin(2.*toang): c2t=cos(2.*toang)
   st2=st*st : ct2=ct*ct
   rp2=rp*rp: vsb=sqr(1/(beta2*beta2)-rp2) 

 if kpwave then 
'---- P-SV radiation:  --------------
    toangs=fnasin(rp*beta2) 'negative of take-off angle of SVup
    s2ts=sin(2.*toangs) : c2ts=cos(2.*toangs)
    vsa=sqr(1/(alph2*alph2)-rp2)
    vslow(1)=vsa : vslow(2)=vsb

'---- amplitudes: F(i), i= 1-P down, 2-P up, 3-SV up
      pcor=4.*pii*rho2*alph2^3: scor=4*pii*rho2*beta2^3*(vsb/vsa)

    e1=st2*(amnn*ca2 + amee*sa2 + amne*s2a) + amdd*ct2
    e2=s2t*(amnd*ca + amed*sa)
    e3=0.5*s2ts*(amnn*ca2+amee*sa2-amdd+amne*s2a)
    e4=c2ts*(amnd*ca+amed*sa)

    f(1)= (e1+e2)/pcor   'P down
    f(2)= (e1-e2)/pcor   'P up
    f(3)= -(-e3+e4)/scor  'SV up

'---- derivatives; i indexing 1-P down, 2-P up, 3-SV up
'                  j indexing 1-Mnn, 2-Mee, 3-Mdd, 4-Mne, 5-Mnd, 6-Med
  if kgetder then
    df(1,1)= st2*ca2/pcor
    df(2,1)= df(1,1)
    df(3,1)= 0.5*s2ts*ca2/scor  'Mnn

    df(1,2)= st2*sa2/pcor
    df(2,2)= df(1,2)
    df(3,2)= 0.5*s2ts*sa2/scor  'Mee

    df(1,3)= ct2/pcor
    df(2,3)= df(1,3)
    df(3,3)= -0.5*s2ts/scor  'Mdd

    df(1,4)= st2*s2a/pcor
    df(2,4)= df(1,4)
    df(3,4)= 0.5*s2ts*s2a/scor  'Mne

    df(1,5)= s2t*ca/pcor
    df(2,5)= -df(1,5)
    df(3,5)= -c2ts*ca/scor  'Mnd

    df(1,6)= s2t*sa/pcor
    df(2,6)= -df(1,6)
    df(3,6)= -c2ts*sa/scor  'Med

  end if

else      '---- SH radiation pattern:  ------------
 vslow(1)=vsb

'---- amplitudes for SH.  F(1)-SH down, F(2)-SH up
      scor=4.0*pii*rho2*beta2^3
      e1=st*(s2a*0.5*(-amnn+amee)+amne*c2a) : e2=ct*(-amnd*sa+amed*ca)
    f(1)=(e1+e2)/scor
    f(2)=(e1-e2)/scor
    f(3)=0


'--- derivatives; i indexing 1-SH down, 2-SH up
'                 j indexing 1-Mnn, 2-Mee, 3-Mdd, 4-Mne, 5-Mnd, 6-Med
 if kgetder then 
    for i=1 to 6 : df(3,i)= 0 : next  'because f(3)=0

'Mnn
    df(1,1)= -0.5*st*s2a/scor
    df(2,1)= df(1,1)
'Mee
    df(1,2)= -df(1,1)
    df(2,2)= df(1,2)
'Mdd
    df(1,3)= 0
    df(2,3)= 0
'Mne
    df(1,4)= st*c2a/scor
    df(2,4)= df(1,4)
'Mnd
    df(1,5)= -ct*sa/scor
    df(2,5)= -df(1,5)
'Med
    df(1,6)= ct*ca/scor
    df(2,6)= -df(1,6)

 end if

end if

end sub

