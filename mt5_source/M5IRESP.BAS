sub m5iresp (ss(1),inst,instinp,tstar,dt,npts,hpfilt,gain) static

shared nzero,npole,zero(),pole()
shared ts, ds, tg, dg, tgain,drive$


' Instruments:  (1) WSSN 15-100 (2) WSSN 30-100 (3) WSSN SP (4) GDSN LP
'               (5) GDSN SP (7) GDSN BB
  
  n2= int( log(npts)/log(2) ) +1
  nt=2^n2
  nfreq=nt/2+1

redim fresp(nfreq,2)
  
  dfreq=1/(nt*dt)

instin=instinp AND 7  'first 3 bits are instrument type code
instva=instinp AND 24 'next 2 bits for vel & accel-type response
'-- defaults for instrument types
if instin<>1 then
 if instin=0 then
  if inst=1 then  ts=15: ds=.93: tg=100: dg=1.0: tgain=15
  if inst=2 then  ts=30: ds=1.5: tg=100: dg=1.0: tgain=30
  if inst=3 then  ts=1.: ds=1.0: tg=1.33: dg=1.0: tgain=1
 end if

 if inst<3 then call wssnlp (fresp(),nfreq,dfreq,gain,ts,ds,tg,dg,tgain)
 if inst=3 then call wssnsp (fresp(),nfreq,dfreq,gain,ts,tg,tgain)
 if inst>3 then call gdsn (fresp(),nfreq,dfreq,gain,nzero,npole,zero(),pole())
elseif instin=1 then
 call forfft(nt,dt,fresp(),ss())
end if

' convert acceleration response type to velocity type
if (instva AND 16)>0 then call resp2d (fresp(),nfreq,dfreq)
' convert velocity response type to displacement type
if (instva AND 24)>0 then call resp2d (fresp(),nfreq,dfreq)

'-- apply attenuation, filter, and inverse fft
  call tstars (fresp(),tstar, dt, nt)
  call but3p (nfreq, dfreq, fresp(), hpfilt)
  call invfft (nt,dt,fresp(),ss())
end sub

'===================================================================
sub tstars(fresp(2),tstar,dt,nt) static
' shared pii
  defint i-n
 '
 '   mantle attenuation correction
 '   
 '   Nabelek
 '
 '  complex catn,ai,c1,c2,fresp
 '
 '
       air=0 : aii=1
       tempr=-0.5 : tempi=0
       twpi=pii+pii
       nfreq=nt/2+1
       if(tstar=0.) then exit sub
      dom=twpi/(nt*dt)
      ommax=twpi/(2.*dt)
      omega=dom

      for i=2 to nfreq
      ot=omega*tstar
'----------- carpenter/nabelek -------------
     t=log(omega/ommax)/pii
      c2r=air*t
      c2i=aii*t
    call comadd (tempr,tempi,c2r,c2i,c1r,c1i)
      c1r=c1r*ot: c1i=c1i*ot
    call comexp(c1r,c1i,catnr,catni)
      f1=fresp(i,1) : f2=fresp(i,2)
    call comult(f1,f2,catnr,catni,f3,f4)
      fresp(i,1)=f3 : fresp(i,2)=f4
      omega=omega+dom
    next 
end sub
'===================================================================
sub wssnlp(fresp(2),nfreq,dfreq,gain,ts,ds,tg,dg,tgain) static
'  shared pii
'
'   Hagiwara's formula for the 30-100 WWSSN instrument

  dim sig15(5),sig30(5),fmag(5)
  sig30(1)=.003:sig30(2)=.012:sig30(3)=.044:sig30(4)=.195:sig30(5)=.767
  sig15(1)=.003:sig15(2)=.013:sig15(3)=.047:sig15(4)=.204:sig15(5)=.805
  fmag(1)=375.:fmag(2)=750.:fmag(3)=1500.:fmag(4)=3000.:fmag(5)=6000.

  bmag=gain
      ai1=0 : ai2=1
      pii2=pii+pii
      for i=2 to 5
        if(fmag(i)>=bmag) then goto w2
      next

w2:     fmi=fmag(i) : fmi1=fmag(i-1)
      dm=fmi-fmi1
   if ts=15 then si=sig15(i) : si1=sig15(i-1)
   if ts=30 then si=sig30(i) : si1=sig30(i-1)
      xds=si-si1
      xdsdm=xds/dm
      dm=bmag-fmi1
      sigma2=si1+xdsdm*dm

      u0=tg/ts
      c0=1.0
      c1=-2.0*(ds+dg/u0)
      c2=1.0+1.0/(u0*u0)+4.0*ds*dg*(1.0-sigma2)/u0
      c3=2.0*(ds/u0+dg)/u0
      c4=1.0/(u0*u0)
      f=dfreq
      df=dfreq
      fresp(1,1)=0. : fresp(1,2)=0. 
    for i=2 to nfreq
      u1=1.0/(ts*f)
      fim=c0-c2*u1^2+c4*u1^4
      fre1=c1*u1+c3*u1^3
      amag=u1/sqr(fre1*fre1+fim*fim)
      phase= fnatan2(fim,fre1)+pii
      fresp(i,1) = amag * exp(ai1*phase) * cos(ai2*phase)
      fresp(i,2) = amag * exp(ai1*phase) * sin(ai2*phase)
      f=f+df
    next
      u1=tgain/ts
      fim=c0-c2*u1^2+c4*u1^4
      fre1=c1*u1+c3*u1^3
      amax=u1/sqr(fre1*fre1+fim*fim)
      ratio=bmag/amax
    for if1=2 to nfreq
      for k=1 to 2 : frr=fresp(if1,k): fresp(if1,k)=frr*ratio : next
  next 
  end sub
'===================================================================
sub wssnsp(fresp(2),nfreq,dfreq,gain,ts,tg,tgain) static
'  shared pii
'
'     short-period instrument response of the WWSSN
'     following L. Burdick's thesis, App. 1 with coupling
'     constant=0.
'     Wang-Ping Chen
'
    ts=1.0
    tg=1.33
    tgain=1.0

'      complex ai,c1,fresp
      air=0. : aii=1.
      fsi=1./ts*pii*2.
      fsi2=fsi*fsi
      fgi=1./tg*pii*2.
      fgi2=fgi*fgi
      esi=0.49*fsi
      egi=(4./3.)*pii*2.
'  circuit parameters
      cl=6.8
      cr=152.5
      r=170.
      s=18.
      s2=s*s
      cqs=r*cr+r*s+s*cr
      fresp(1,1)=0. : fresp(1,2)=0.
      freq=0.
      dfr=dfreq*2.*pii
    for if1=2 to nfreq
      freq=freq+dfr
      freq2=freq*freq
      u=cl*(r+s)*freq/cqs
      cf=1./sqr(1.+u*u)
      cf2=cf*cf
      fssq=fsi2+2.*esi*freq*u*cf2
      ax=freq*cl*s2/cqs/(cr+s)
      fgsq=fgi2+2.*egi*ax*freq*cf2
      es=esi*cf2
      ay=1.+cl*freq*u/(cr+s)
      eg=egi*ay*cf2
      bxs=fssq-freq2
      bxg=fgsq-freq2
      cb=bxs*bxg-4.*es*eg*freq2
      cc=(es*bxg+eg*bxs)*2.*freq
      amp=-1.*freq2*freq*cf/sqr(cc*cc+cb*cb)
      ph=fnatan2(cb,cc)-atn(u)
      c1r=air*ph : c1i=aii*ph  
      call comexp(c1r,c1i,f1,f2)
      fresp(if1,1)=f1*amp : fresp(if1,2)=f2*amp
  next ' if1 continue
      freq=1./tgain*2.*pii
      freq2=freq*freq
      u=cl*(r+s)*freq/cqs
      cf=1./sqr(1.+u*u)
      cf2=cf*cf
      fssq=fsi2+2.*esi*freq*u*cf2
      ax=freq*cl*s2/cqs/(cr+s)
      fgsq=fgi2+2.*egi*ax*freq*cf2
      es=esi*cf2
      ay=1.+cl*freq*u/(cr+s)
      eg=egi*ay*cf2
      bxs=fssq-freq2
      bxg=fgsq-freq2
      cb=bxs*bxg-4.*es*eg*freq2
      cc=(es*bxg+eg*bxs)*2.*freq
      amp=freq2*freq*cf/sqr(cc*cc+cb*cb)
      ratio=abs(gain/amp)
    for if1=2 to nfreq
      for k=1 to 2 : fresp(if1,k)=fresp(if1,k)*ratio : next
    next        'continue
 end sub
'===================================================================
sub gdsn(fresp(2),nfreq,dfreq,gain,nzero,npole,zero(1),pole(1)) static

'  shared pii

    dom=2*pii*dfreq
    om=dom

    for i=2 to nfreq
      sr= 0.0 : si=om
      vr1=gain : vi1=0.0

      '---- zeroes
      for jf=1 to nzero
         zr=zero(2*jf-1)  : zi=zero(2*jf)
         vr2=sr-zr : vi2=si-zi
         call comult (vr1,vi1,vr2,vi2,vr3,vi3)
         vr1=vr3 : vi1=vi3
      next

      '---- poles
      for jf=1 to npole
         pr=pole(2*jf-1) : pi=pole(2*jf)
         vr2=sr-pr : vi2=si-pi
         call cominv(vr2,vi2,vr,vi)
         call comult (vr1,vi1,vr,vi,vr3,vi3)
         vr1=vr3 : vi1=vi3
      next

      fresp(i,1)=vr1
      fresp(i,2)=vi1
      om=om+dom
   next
   fresp(1,1)=0. : fresp(1,2)=0. 
end sub
'===================================================================
sub resp2d (fresp(2),nfreq,dfreq)
' convert velocity-type response function to displacement-type
' by Zwick -94/5
'
' fresp = a*cos(omega*t) + b*sin(omega*t)
'
' need to do conversion frespx() = frespv() * F(f)
' Known are: v=cos(omega*t)
' and output= G(f)= Ar*cos(omega*t)+Ai*sin(omega*t),
' where Ar, Ai are real and imag coef of the velocity response
'
' to get this same output with x=integral(v)=1/omega*sin(omega*t),
' we need to multiply coef by omega and add pi/2 phase shift
' the new transfer function is: G'(f)=(-omega*Ai,omega*Ar).

  dom=2*pii*dfreq
  om=0.
  for i=1 to nfreq
    a = fresp(i,1)
    fresp(i,1)= -fresp(i,2)*om
    fresp(i,2)= a*om
    om=om+dom
  next
end sub
'===================================================================
rem $INCLUDE: 'm5comfun'
rem $INCLUDE: 'm5fft'
'===================================================================

