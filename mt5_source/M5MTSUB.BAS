'-- M5MTSUB.BAS  Subroutines for Moment Tensor Plotting in MT5PLOT
'   By Peter Zwick        Created Nov, 1991
'   Last Modified:  Dec 2,1991

DEFINT I-N: OPTION BASE 1
const npltpt=1000 'points to plot

'following are shared with MT5PLOT module
COMMON /PLOT/ fisx,fisy,gscale,io$,irotate,khpcol,kopen7,_
  lcode,letset,prad,xmax,xscale,yscale

rem $INCLUDE: 'M5COMMON'
rem $INCLUDE: 'M5FUNCTN'

def fnisx(x) = cint(x*fisx)
def fnisy(y) = cint(kscry-(y-scrorig)*fisy)

redim xtmp1(npltpt),ytmp1(npltpt),xtmp2(npltpt),ytmp2(npltpt)
redim n(2),kl(2),pinc(2),ndplt(2)

'========================================================================
sub PLOTPMTRAD(x0,y0,prad) static
'-------------- plot P-wave nodal plane for MT solution
shared scratio,sorpar(),nsor
shared xtmp1(),ytmp1(),xtmp2(),ytmp2(),n(),i,ndplt()

redim ixhld(2),iyhld(2),ndplt(2),newln(2)
redim xtmp1(npltpt),ytmp1(npltpt),xtmp2(npltpt),ytmp2(npltpt),n(2),kl(2),pinc(2)

azjincmax=5 : azjincmin=.5
for isor=1 to nsor
  redim phi(5)
  amnn=sorpar(1+noffmt,isor)
  amee=sorpar(2+noffmt,isor)
  amdd=sorpar(3+noffmt,isor)
  amne=sorpar(4+noffmt,isor)
  amnd=sorpar(5+noffmt,isor)
  amed=sorpar(6+noffmt,isor)
  call linetype(isor,kline,pincr)

'-- initialize variables
'-- identify azimuths (phi) with radially directed nodal plane segments
  phitest=amne*amne-amee*amnn
  if phitest>=0 and abs(amee)>1e-5 then
    phitest=sqr(phitest)
    k=1
    for i=1 to -1 step -2
      phi(k)=fndeg(atn((-amne+i*phitest)/amee))
      phi(k+1)=180.+phi(k)
      for j=0 to 1 'make positive
        if phi(k+j)<0 then phi(k+j)=360.+phi(k+j)
      next j
      k=k+2
    next i
    for i=1 to 3  'sort in increasing order
      for j=i+1 to 4
        if phi(i)>phi(j) then SWAP phi(i),phi(j)
    next j,i
  end if
  azj=phi(1)
  phi(5)=azj+360.
  k=1  'phitest index
  do
    az=fnrad(azj) : ca=cos(az) : sa=sin(az)
    a=amnn*ca*ca+amee*sa*sa+amne*sin(2*az)
    b=2*(amnd*ca+amed*sa)
    c=amdd
    c2=b*b-4*a*c
    plast=piihalf
    p=piihalf
    if c2>0 or (c2=0 and abs(b)>=1e-7) then 
      sqrc=sqr(c2)
      imin=1 : imax=2
      if abs(a)<1e-7 then
        if b>0 then imax=1 else imin=2
      end if
      for i=imin to imax
        newln(i)=kfalse
        if i=1 then ix=1 else ix=-1
        if abs(a)>=1e-7 then 
          p=(-b+ix*sqrc)/(2*a) 
          p=ATN(p)
        else
          p=ATN(-c/b)
        end if
        if fninlim(p,0,piihalf) then
          n(i)=n(i)+1 
          newln(i)=(n(i)=1)
          call eaproj2(x0,y0,prad,p,az,x,y)
          if i=1 then 
            xtmp1(n(i))=x : ytmp1(n(i))=y
          else
            if azj>phi(1) then
'if two roots just found, connect them
              if newln(1) and newln(2) then
                xtmp2(1)=xtmp1(1) : ytmp2(1)=ytmp1(1)
                n(i)=n(i)+1
'if one root just appeared, connect to circle
              elseif newln(1) then
                xtmp1(2)=xtmp1(1) : ytmp1(2)=ytmp1(1) : p=piihalf
                call eaproj2(x0,y0,prad,p,az,xtmp1(1),ytmp1(1))
                n(1)=2
              elseif newln(2) then
                p=piihalf
                call eaproj2(x0,y0,prad,p,az,xtmp2(1),ytmp2(1))
                n(2)=2
              end if
            end if
            xtmp2(n(i))=x : ytmp2(n(i))=y
          end if
        else 
          ndplt(i)=(n(i)>0) 'plot it
        end if
      next i
    else 
      for i=1 to 2 : ndplt(i)=(n(i)>0) : next i
    end if
'-- plot it if needed
    if ndplt(2) then
      n(2)=n(2)+1
      if ndplt(1)  then 'connect lines
        xtmp2(n(2))=xtmp1(n(1)) : ytmp2(n(2))=ytmp1(n(1))
      else   'draw line to circle
        p=piihalf
        call eaproj2(x0,y0,prad,p,az,xtmp2(n(2)),ytmp2(n(2)))
      end if
    elseif ndplt(1) then   'draw line to circle
      n(1)=n(1)+1
      p=piihalf
      call eaproj2(x0,y0,prad,p,az,xtmp1(n(1)),ytmp1(n(1)))
    end if
    for i=1 to 2
      if ndplt(i) or n(i)>=npltpt then call plotnow(kline,pincr,kcoldat)
    next i
'-- decide on a step size that will get closest to the outer radius
    azjinc=azjincmin
    if azj>=phi(k+1)-.1 then 
      k=k+1
    elseif azj-phi(k)>=azjincmax and phi(k+1)-azj>=azjincmax then
      azjinc=azjincmax
    elseif phi(k+1)-azj<=azjincmin then
      azjinc=phi(k+1)-azj
    end if
    if abs(c2)<.04 then azjinc=fnamin(azjinc,azjincmin)
    azj=azj+azjinc
  loop until azj>phi(5)+.1
  for i=1 to 2 'finish up lines
    if n(i)>0 then call plotnow(kline,pincr,kcoldat)
  next i
next isor
end sub
'========================================================================
sub PLOTSMTRAD(x0,y0,prad) static
'-------------- plot SH-wave nodal plane for MT solution
shared scratio,sorpar(),nsor
shared xtmp1(),ytmp1(),xtmp2(),ytmp2(),n(),i,ndplt()

redim ixhld(2),iyhld(2),ndplt(2)
redim xtmp1(npltpt),ytmp1(npltpt),xtmp2(npltpt),ytmp2(npltpt),n(2),kl(2),pinc(2)

azjincmax=3 : azjincmin=1
for isor=1 to nsor
  amnn=sorpar(1+noffmt,isor)
  amee=sorpar(2+noffmt,isor)
  amdd=sorpar(3+noffmt,isor)
  amne=sorpar(4+noffmt,isor)
  amnd=sorpar(5+noffmt,isor)
  amed=sorpar(6+noffmt,isor)
  call linetype(isor,kline,pincr)

'-- initialize variables
  azj=0
  azjinclast=10
  do
    az=fnrad(azj)
    a=sin(2*az)/2*(-amnn+amee)+amne*cos(2*az)
    b=-amnd*sin(az)+amed*cos(az)
    plast=piihalf
    p=piihalf
    i=1
    if abs(a)>1e-7 then 
      p=ATN(-b/a)
    elseif abs(b)>1e-7 then 
      p=-sgn(b)*piihalf
    end if
    if fninlim(p,0,piihalf) then
      n(i)=n(i)+1
      call eaproj2(x0,y0,prad,p,az,x,y)
      xtmp1(n(i))=x : ytmp1(n(i))=y
      plast=fnamin(plast,p)
    else 
      if n(i)>0 then call plotnow(kline,pincr,kcoldat)
    end if
    if n(i)>=npltpt then call plotnow(kline,pincr,kcoldat)
'-- decide on a step size that will get closest to the outer radius
    dp=fnamax(abs(plast-p),1e-5)
    azjinc=fnamin(fnamin((piihalf-p)*azjinclast/dp-.15,azjincmax),360-azj)
    azjinc=fnamax(azjinc,azjincmin)
    azjinclast=azjinc
    azj=azj+azjinc
  loop until azj>360
  if n(i)>0 then call plotnow(kline,pincr,kcoldat)
next isor
end sub
'========================================================================
sub eaproj2 (x0,y0,rad0,ang,az,x,y) static
    r2=rad0*sqr2*sin(ang/2)
    x=x0 + r2*cos(az-piihalf)
    y=y0 - r2*sin(az-piihalf)
end sub
'============================================================================
sub scrscale(fisx,fisy) static
'---- convert plot coords (250x200mm) to undistorted screen coordinates
shared kscry,kscrx,scratio
  if (kscry-200*scratio*kscrx/250)<0 then  'if screen too wide
     fisx=kscry/(scratio*200)  '           'pixels/(mm x dimension)
     fisy=kscry/200                        'pixels/(mm y dimension)
  else
     fisx=kscrx/250
     fisy=kscrx*scratio/250
  end if
end sub
'-----------------------------------------------------------------------------
sub PLOTPRAD(x0,y0,prad) static '--- plot nodal planes
shared sorpar(),nsor,kletter,evct(),mtflag,fisx,fisy,gsrot,kgrcol()
redim xtmp(npltpt), ytmp(npltpt)
call scrscale(fisx,fisy)
for isor=1 to nsor
  str=sorpar(1,isor)
  if str>-998 then
    call linetype(isor,kline,pincr)
    dip=sorpar(2,isor)
    rk=sorpar(3,isor)

    for npl=1 to 2
      dip=fnamin(dip,piihalf-.001)
      n=1
      xtmp(1)=x0+prad*cos(str-piihalf)
      ytmp(1)=y0-prad*sin(str-piihalf)

'-- plot plane
      for i=0 to 180 step 3
        n=n+1
        d1=fnrad(i) + str
        p=piihalf-atn(tan(dip)*sin(fnrad(i)))
        call eaproj2(x0,y0,prad,p,d1,x,y)
        xtmp(n)=x : ytmp(n)=y
      next i

      call gsline (xtmp(), ytmp(), n, kline, pincr, kcolbak)
      if npl=1 then  'plot P and T axes
        if mtflag then
'---- get P and T axes
          taz=fnatan2(evct(1,2,isor),evct(1,1,isor))
          tpl=fnatan2(evct(1,3,isor),sqr(evct(1,2,isor)*evct(1,2,isor)+_
              evct(1,1,isor)*evct(1,1,isor)))
          if evct(1,3,isor)<0 then taz=taz+piihalf:tpl=-tpl
          while taz>pii : taz=taz-pii : wend
          paz=fnatan2(evct(2,2,isor),evct(2,1,isor))
          ppl=fnatan2(evct(2,3,isor),sqr(evct(2,2,isor)*evct(2,2,isor)+_
              evct(2,1,isor)*evct(2,1,isor)))
          if evct(2,3,isor)<0 then paz=paz+piihalf:ppl=-ppl
          while paz>pii  : paz=paz-pii : wend
          kc=kgrcol(kcoldat)
        else
          call m5ptbax (str,dip,rk,paz,ppl,taz,tpl,baz,bpl)
          kc=kgrcol(kcolbak)
        end if
'--- plot T axis as open circle
        tpl=piihalf-tpl : call eaproj2 (x0,y0,prad,tpl,taz,x,y)
        ix=fnisx(x):iy=fnisy(y)
        axsz=prad/8
        call gssym (kfalse,x,y, 33, gsrot, axsz)
        call sccircle(ix,iy,3,kc)
        if not kletter then call gssym (kfalse, x, y, 40, gsrot, axsz)
'--- plot P axis as filled 
        ppl=piihalf-ppl : call eaproj2 (x0,y0,prad,ppl,paz,x,y)
        ix=fnisx(x):iy=fnisy(y)
        call gssym (ktrue, x, y, 33, gsrot, axsz)
        call sccircle(ix,iy,3,kc)
        call sccircle(ix,iy,1,kc)
'--- get the auxiliary plane
        call m5auxpl(str,dip,rk,strt,dipt,rkt)
        str=strt : dip=dipt : rk=rkt
      end if

    next npl
  end if 'str>-998
next isor
end sub
'return
'============================================================================
sub plotnow(kline,pincr,kcol) static
shared xtmp1(),ytmp1(),xtmp2(),ytmp2(),n(),i,ndplt()

if i=1 then
  call gsline (xtmp1(), ytmp1(), n(1), kline, pincr, kcol)
else
  call gsline (xtmp2(), ytmp2(), n(2), kline, pincr, kcol)
end if
n(i)=0 : ndplt(i)=kfalse
end sub

'-----------------------------------------------------------------------------
sub PLOTSRAD(x0,y0,prad) static'----- subroutine to plot S radiation pattern (fixed 4/30/86)
shared sorpar(),nsor,kletter
shared cd,cs,c2d,ss,sd,s2d,dz,str,f,rp,vsb,a1,a2,sa,ca,s2a,c2a
shared ain,az,x,y    

redim xtmp(npltpt), ytmp(npltpt)

for isor=1 to nsor
  call linetype(isor,kline,pincr)
  n=0

  str=sorpar(1,isor)
  if str>-998 then

    dip=sorpar(2,isor) : rk=sorpar(3,isor) 
    sd=sin(dip): cd=cos(dip): s2d=sin(2.*dip): c2d=cos(2.*dip)
    ss=sin(rk): cs=cos(rk)
    rp=0 : vsb=0 

    '--- singular cases:
    if dip<1e-7 then       'dip=0 deg
       ain=90 : az=fndeg(str-rk) : call PROJNET(x0,y0,prad) : xtmp(1)=x : ytmp(1)=y
       ain=-ain : call PROJNET(x0,y0,prad) : xtmp(2)=x : ytmp(2)=y
       n=2 : go to 3180
    elseif s2d<1e-7 and abs(cs)<1e-7 then   'dip=0 or 90 deg and rk=90 or 270 deg
       ain=90 : az=fndeg(str)+90 : call PROJNET(x0,y0,prad) : xtmp(1)=x : ytmp(1)=y
       ain=-ain : call PROJNET(x0,y0,prad) : xtmp(2)=x : ytmp(2)=y
       n=2 : go to 3180
    elseif cd<1e-7 and abs(ss)<1e-7 then  'dip=90 deg and rk=0 or 180 deg
       ain=90 : az=fndeg(str)+45 : call PROJNET(x0,y0,prad) : xtmp(1)=x : ytmp(1)=y
       ain=-ain : call PROJNET(x0,y0,prad) : xtmp(2)=x : ytmp(2)=y
       n=2 : call gsline(xtmp(),ytmp(),n,kline,pincr,kcolbak)
       az=fndeg(str)-45 : call PROJNET(x0,y0,prad) : xtmp(1)=x : ytmp(1)=y
       ain=-ain : call PROJNET(x0,y0,prad) : xtmp(2)=x : ytmp(2)=y
       go to 3180
    elseif abs(c2d)<1e-7 and abs(cs)<1e-7 then     'dip=45 deg and rk=90 or 270 deg
       ain=90 : az=fndeg(str) : call PROJNET(x0,y0,prad) : xtmp(1)=x : ytmp(1)=y
       ain=-ain : call PROJNET(x0,y0,prad) : xtmp(2)=x : ytmp(2)=y
       n=2 : call gsline(xtmp(),ytmp(),n,kline,pincr,kcolbak)
       az=fndeg(str)+90 : call PROJNET(x0,y0,prad) : xtmp(1)=x : ytmp(1)=y
       ain=-ain : call PROJNET(x0,y0,prad) : xtmp(2)=x : ytmp(2)=y
       go to 3180
    end if
    
    notyet=ktrue : sgnold=999 : n=0
    azstep=1.0 : dz=-azstep

3080 '--- start loop for radiation pattern
    dz=dz+azstep : if dz > 180. then 3180
    az=dz
    call AMPLSH
    if abs(a1)<1e-7 then  '--- special case : ain = 90
       if notyet then ain=90. : call PROJNET(x0,y0,prad) : goto 3120
       if a2 <>0 then ain=90.*sgnold : call PROJNET(x0,y0,prad) : goto 3120
       ain=-90. : call PROJNET(x0,y0,prad) : xx=x : yy=y 
       ain=90. : call PROJNET(x0,y0,prad) 
       if n>0 then call gsline(xtmp(),ytmp(),n,kline,pincr,kcolbak)
       n=2
       xtmp(1)=xx : ytmp(1)=yy
       xtmp(2)=x : ytmp(2)=y 
       call gsline(xtmp(),ytmp(),n,kline,pincr,kcolbak) 
       n=0 
       goto 3130 
    end if
     
3100 ain=fndeg(atn(-a2/a1))
    call PROJNET(x0,y0,prad)
    if  notyet then 3115           '----  first time
    fcn=-a2/a1
    sdfdaz=sgn(((-2*sd*cs*s2a-s2d*ss*c2a)*a2-(cd*cs*ca-c2d*ss*sa)*a1)/a1/a1)
    sdshdrp=sgn(a1*cos(fnrad(ain))-a2*sin(fnrad(ain)))
    if sdfdaz=sgn(ain-ainold) then 3120        '--- check for discontinuity
       if sdfdazold<>sdfdaz and sdshdrpold=sdshdrp then 3120
         xx1=x : yy1=y : ain=90.*sdfdaz : call PROJNET(x0,y0,prad) 
         n=n+1 : xtmp(n)=x : ytmp(n)=y
         call gsline (xtmp(),ytmp(),n,kline,pincr,kcolbak)
         n=0
         ain=-ain : call PROJNET(x0,y0,prad)
         xtmp(1)=x : ytmp(1)=y : xtmp(2)=xx1 : ytmp(2)=yy1 : n=2
         goto 3130

3115 '--- first point, just move
         n=n+1 : xtmp(1)=x : ytmp(1)=y
         x1=x : y1=y
         goto 3130

3120 n=n+1 : xtmp(n)=x : ytmp(n)=y

3130 sgnold=sgn(ain) : ainold=ain : sdfdazold=sdfdaz : sdshdrpold=sdshdrp

    notyet=kfalse
    goto 3080

3180 call gsline (xtmp(),ytmp(),n,kline,pincr,kcolbak)
  end if 'str>-998
next isor
end sub
'-----------------------------------------------------------------------------
sub PROJNET(x0,y0,prad) static ' find x,y for projection
shared ain,az,x,y    
    r2=prad*sqr2*sin(fnrad(ain)/2)
    x=x0 + r2*cos(fnrad(az)-piihalf)
    y=y0 - r2*sin(fnrad(az)-piihalf)
end sub
'-----------------------------------------------------------------------------
sub AMPLSH static '---- amplitudes for SH
shared cd,cs,c2d,ss,sd,s2d,dz,str,f,rp,vsb,a1,a2,sa,ca,s2a,c2a
 a=fnrad(dz)-str : sa=sin(a): ca=cos(a): s2a=sin(2.*a): c2a=cos(2.*a)
 a2=cd*cs*sa+c2d*ss*ca : a1=sd*cs*c2a-.5*s2d*ss*s2a
 f=a1*rp + a2*vsb
end sub
'-----------------------------------------------------------------------------
sub linetype (nsor,kline,ds) static
'-- get line type KLINE and increment DS for source NSOR
  if nsor=1 then kline=2 : ds=0.00
  if nsor=2 then kline=0 : ds=1.30
  if nsor=3 then kline=1 : ds=1.30
  if nsor=4 then kline=0 : ds=2.60
  if nsor=5 then kline=1 : ds=2.60
  if nsor=6 then kline=0 : ds=3.90
end sub
'-----------------------------------------------------------------------------
sub getplotnm(gsin$,gsf$,gss$)
shared iplform,i2page,iphase,kgsnul
  gsf$=alltrim$(lcase$(gsin$))
  gss$=fnsubext$(gsf$)
  kgsnul=(len(gss$)=0 or gss$="nul")
  if kgsnul then 
     gsf$="nul"
     gss$="nul"
     exit sub
  elseif iplform=1 then
    if i2page then 
      if iphase=1 then gsext$=".p" else gsext$=".sh"
    else 
      gsext$=".psh"
    end if
  else
    gsext$=".phm"
  end if
  gsf$=gss$+gsext$
end sub
'-----------------------------------------------------------------------------
rem $INCLUDE: 'M5GRSUBS'

