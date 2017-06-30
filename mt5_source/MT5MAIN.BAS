' main MT5MAIN (MT5.00)

' teleseismic body waveform inversion program for PC's written by
'   Rob McCaffrey and Geoff Abers
'     SYN3.01 finished April 1988,
'     revision SYN3.02 June 1989
'     revision SYN3.04 June 1990
'     revision SYN4.00 September 1990 (with Peter Zwick)
'   Revised Jan 1991 to be 2nd module in MT4 chaining to conserve memory use -PZ
'   Revised Jan 1993 to check for .INV/.PRM file mismatch -PZ
'   Revised Mar 1993 to null variable parameter list when reading parameters -PZ
'   Revised Nov 1993 to use long integer in calculating maxsyn, and
'     changed calculation of memory required  -PZ
'   Link with m5sub and specify QB as library during link

DEFINT I-N: OPTION BASE 1

' additional common block needed by m5iresp
dim zero(100), pole(100)
COMMON /INSTRESP/ nzero,npole,zero(),pole(),ts,ds,tg,dg,tgain

rem $INCLUDE: 'M5COMMON'
rem $INCLUDE: 'M5FUNCTN'

if not imt5 then call notalone

nraymax=200
dim refcs(4,4), aincs(4)
dim sor(nraymax)
dim work(nraymax)

DEF FNINSAME(i,j)
'-- decide if the instruments, phase, t*, etc. are the same to
'    cut down on recalculating instrument responses

  static ins
  ins= (iph(j)=iph(i)) and (inst(j)=inst(i)) and (tstar(j)=tstar(i))
  ins= (ins) and (dt(j)=dt(i)) and (tau(1,j)=tau(1,i)) and (tau(2,j)=tau(2,i)) and (tau(3,j)=tau(3,i))
  ins= (ins) and (ndat(j)>=ndat(i)) and (hpfilt(j)=hpfilt(i))
  ins=(ins) and ((inst(i)>2) or (xmag(j)=xmag(i))) and (instin(i)=instin(j))

  if fningdsn(i) and fningdsn(j) then
    ins= (ins) and nzero(i)=nzero(j)
    ins= (ins) and npole(i)=npole(j)
    ins= (ins) and calcons(i)=calcons(j)
    if ins then
      ip1=indexzp(i)-1 : jp1=indexzp(j)-1
      for k=1 to 2*nzero(i)+2*npole(i) 
        ip1=ip1+1 : jp1=jp1+1 
        ins= (ins) and zeropole(ip1)=zeropole(jp1) 
      next
    end if
  end if

  fninsame=ins
end def

'-- start

redim tmp(nesmax)

w$="i"
101 gosub MENU ' menu

if mkgrn then
  gosub SUBINST ' get instrument responses
end if

if mkgrn or mktau then
  gosub CRESP ' get Green's functions
  gosub TAUCONV ' convolve with TAU
end if
mkgrn=kfalse: mktau=kfalse: mkdat=kfalse
chain syndir$+"MT5MEN1"
end

MENU: '-- setup menu
call setscreen(0)
call scclear
klst$="igtvhwlcrexqsam"
call scbox(1,21,ktrue,"SETUP MENU","Press C to Continue to MAIN MENU")
call scprint(2,2,"KEY  OPTION")         : call scprint(2,55,"CURRENT VALUE")
call scprompt(3,3,"I - Input/output file prefix")
call scprompt(4,3,"G - crustal response functions to be Generated")
call scprompt(5,3,"T - half-duration of Triangular STF element (sec)")
call scprompt( 6,3,"V - Velocity, azimuth and plunge of line source")
call scprompt( 7,3,"H - source depth for calculating ray parameters")
call scprompt( 8,3,"W - Water depth in km")
call scprompt( 9,3, "L - crustal structure  Layer #:(Ztop, Vp, Vs, Rho)")
if nlyr>1 then
  call scprompt(11,3,"N - Number of complete reverberations in layer 1")
  call scprintv(11,55, str$(nrevcom))
  klst$=klst$+"n"
end if
call scprompt(12,3, "C - Continue to Main Menu")
call scprompt(13,3, "R - Read parameter (.PRM) and data (.INV) files")
call scprompt(14,3, "E - change default directories and screen settings")
call scprompt(15,3,"X/"):call scprompt(15,5,"Q - Shell/Quit")

call scprint(17,7,"PHASE-DATA defaults    ")
  for i=1 to nphtype
    print "  ";fnkphlab$(i);" ";
  next
call scprompt(18,3,"S - Seismogram length, sec")
call scprompt(19,3,"A - Attenuation t*, sec")
call scprompt(20,3,"M - Magnification")

'-- now write variables
call setvarcol
call scprint(3,56,fhdr$)
call scprint(4,56,fnyesno$(mkgrn))
call scprint(5,55,str$(tau0))
locate 6,55:print using fnfform$(6,1); vrupt; azrupt; plrupt
call scprint(7,55,str$(deptot))
call scprint(8,55,str$(wthk))
for i=1 to nlyr
  locate  8+i,55: print using " #:(";i;
  print using fnfform$(5,1); ztop(i)+wthk; alph(i); beta(i); rho(i);
  print ")";
next
locate 18,29
for i=1 to nphtype
  print using fnfform$(7,1);phrlen(i);
next
locate 19,29
for i=1 to nphtype
  print using fnfform$(7,2);phtstar(i);
next
locate 20,29
for i=1 to nphtype
  print using fniform$(7);phmag(i);
next
call settxtcol


'-- end of printing menu
prompt1:
if input1(22,1,"Enter KEY for option ",w$,kfalse,kfalse,ktrue,klst$)=3 then w$="q"
'Ctrl-C checked

if fniquit(w$) then call ENDCHECK:w$="":goto prompt1
if w$="c" then goto MENUEND

if w$="a" then
  cls: print crin$
  irow=0:icol=0
  call scprintc(irow,icol,"Enter default t* (ESC to exit)") : mkgrn=ktrue

  for i=1 to nphtype
    t$=str$(phtstar(i))
    if inputn(irow+i,1,fnkphlab$(i),t$,ktrue,ktrue,kfalse,"f",8)=27 then exit for
    if len(t$)>0 then
      t=fnxlimit(val(t$),0,99999)
      if phtstar(i)<>t then phtstar(i)=t : mkinst=ktrue
    end if
  next

  if not nodat then
    if igetyesno(14,2,"Change t* for individual stations","y",ktrue) then ich=ktrue : gosub CHGVALUES
  end if

elseif w$="n" and nlyr>1 then
  a$=str$(nrevcom)
  if inputn((0),2,"Number of reverberations to calculate in layer 1 ",_
   a$,ktrue,ktrue,ktrue,"i",3)=27 then goto prompt1
  n=val(a$)
  if n<>nrevcom then mkgrn=ktrue
  nrevcom=fnamax(1,n)

elseif w$="e" then
  ichoice=asc(w$): chain syndir$+"MT5"

elseif w$="g" then
  if isfile((drive$+"gfxi.tmp"),(""),1)=0 and isfile((drive$+"inst.tmp"),(""),1)=0 then
    print "If the .TMP files exist from a previous session, you can use the stored response
    print "  functions. Be sure that the parameters are identical.":y$="y"
    if input1((0),(0),"Generate crustal/instrument responses (y/n) ",y$,kfalse,ktrue,ktrue,"yn")<>27 then
      mkgrn=fnkyes(y$) : mkinst=mkgrn
    end if
  else
    a$=displaymsg$(drive$+"   .TMP files not found. You must generate new ones.","",kfalse)
    mkgrn=ktrue : mkinst=ktrue
  end if

elseif w$="h" then
  print "*** Source must be in lowest layer (halfspace) ***"
  a$=str$(deptot)
  if inputn((0),(0)," Enter source depth estimate for ray parameters",_
   a$,ktrue,ktrue,ktrue,"f",6)<>27 then deptot=val(a$)

elseif w$="i" then
  f$=fhdr$:irow=0:call scprintc(irow,1,"I/O file prefix (include pathname but no extension)")
  if inputn(irow+1,2,"--",f$,ktrue,ktrue,ktrue,"a",40)<>27 then
    if len(f$)>0 then
      f$=fnsubext$(f$)
      fhdr$=f$ : mkdat=ktrue : iparflg=kfalse:nsynmax=15000  'reset nsynmax
    end if
    w$="r"
  end if

elseif w$="l" then
  cls:d$=str$(nlyr)
  print "*** Source must be in lowest layer, which is a halfspace ***":irow=0
  if input1(irow,(0),"Number of layers(1-"+fnstrnum$(nlyrmax)+")",d$,kfalse,ktrue,ktrue,numlist$(1,nlyrmax))=27 then
    n=0
  else
    n=val(d$)
  endif
  if fninlim(n,1,nlyrmax) then
    nlyr=n : mkgrn=ktrue
    for i=1 to nlyr
      irow=irow+5:call scprint(irow,4,"For Layer"+str$(i)+" - ") :icol=0
      s$=str$(ztop(i))
      if inputn(irow,icol,"Enter depth of top of layer (km)",s$,ktrue,ktrue,kfalse,"f",6)<>27 then z=val(s$): if z>=0 then ztop(i)=z
      s$=str$(alph(i))
      if inputn(irow+1,icol,"Enter P wave velocity  (km/s)   ",s$,ktrue,ktrue,kfalse,"f",6)<>27 then z=val(s$): if z>=0 then alph(i)=z
      s$=str$(beta(i))
      if inputn(irow+2,icol,"Enter S wave velocity  (km/s)   ",s$,ktrue,ktrue,kfalse,"f",6)<>27 then z=val(s$): if z>=0 then beta(i)=z
      s$=str$(rho(i))
      if inputn(irow+3,icol,"Enter density rho (g/cm^3)      ",s$,ktrue,ktrue,kfalse,"f",6)<>27 then z=val(s$): if z>=0 then rho(i)=z
      if i=1 then z1=ztop(i)
      ztop(i)=ztop(i)-z1
    next
  end if

elseif w$="m" then
  cls:print crin$
  print "Enter default reference magnification for:"
  for i=1 to nphtype
    print fnkphlab$(i);" ";
    phmag(i)=fncrin(phmag(i))
  next
  if not nodat then
    if igetyesno(14,2,"Change Reference Magnification for individual stations","y",ktrue) then ich=kfalse : gosub CHGVALUES
  end if

elseif w$="r" then
  cls
  gosub M5RDPARM
  isorm=1 'start on source #1
  deptot=sorpar(4,1)+cthk
  mkinst=ktrue : mkgrn=ktrue: itf=0
  for i=1 to nsormax: vrp$(i)="": next
  if mkdat then gosub RDDATA 'get data
  w$="c"
elseif w$="s" then
  cls: print crin$
  print "Enter length of calculated seismograms (in seconds) for:"
  for i=1 to nphtype
    print fnkphlab$(i);" ";
    phrlen(i)=fncrin(phrlen(i))
  next
  mkinst=ktrue
  if not nodat then
    for ista=1 to nsta
      nes= phrlen(kphtype(ista))/dt(ista)
      nes(ista)=fnamin(nes,ndat(ista))
      call chkbound(ista)
    next
  end if

elseif w$="t" then
  s$=str$(tau0 )
  if inputn((0),(0),"Enter TAU, the duration of the STF element (in seconds)",s$,ktrue,ktrue,ktrue,"f",6)<>27 then
    s=val(s$)
    if (s<>tau0) and (s>0) then mktau=ktrue : tau0=s
  end if

elseif w$="v" then
  irow=0:x=vrupt:y=azrupt:z=plrupt
  s$=str$(x)
  if inputn(irow,1,"Rupture velocity (V=0 for no directivity)",s$,ktrue,ktrue,kfalse,"f",6)<>27 then x=val(s$)
  s$=str$(y)
  if inputn(irow+1,1,"Rupture azimuth",s$,ktrue,ktrue,kfalse,"f",6)<>27 then y=val(s$)
  s$=str$(z)
  if inputn(irow+2,1,"Rupture plunge",s$,ktrue,ktrue,ktrue,"f",6)<>27 then z=val(s$)
  if x<>vrupt or y<>azrupt or z<>plrupt then
    mktau=ktrue : vrupt=fnamax(x,0)
    azrupt=fnxlimit(y,-360,360) :if azrupt<0 then azrupt=360+azrupt
    plrupt=fnxlimit(z,-90,90)
  end if

elseif w$="w" then
  s=wthk:s$=fnstrnum$(s)
  if inputn((0),(0),"Water depth in kilometers",s$,ktrue,ktrue,ktrue,"f",6)<>27 then s=val(s$)
  if s<>wthk then wthk=fnamax(s,0) : mkgrn=ktrue
  if wthk>20 then wthk=wthk/1000.

elseif w$="x" then call SHELLCOM

end if

'-- depth of top of source half-space
if nlyr>1 then
  sthk=ztop(2)
else
  sthk=0
end if
cthk=wthk+sthk
deptot = fnamax(deptot, cthk+.1)
goto MENU

MENUEND:
if nodat then
  locate 23,1 :  print " *** No data file read in, use option R ***"
  call scpause(24,1) : goto MENU
end if
alph2=alph(nlyr) : beta2=beta(nlyr) : rho2=rho(nlyr)
return

'======================================================================

RDDATA: '--Enter data from file
print:print: d$=ucase$(fhdr$)+".INV":irow=0
call scprint(irow,1,"Default seismogram amplitude data file is "+d$)
call scprintc(irow+2,2,"Enter data (.INV) file (Esc or Q if none)")
if inputn(irow+3,2,"--",d$,kfalse,ktrue,kfalse,"a",40)=27 then return
if fniquit(d$) then return
print:

'-- assign input file name
if len(d$)=0 then d$=fhdr$+".inv"

df$=fnaddext$(d$,"inv")

if ifopen("I",1,df$,0,-1)=0 then go to RDDATA

fhdr$=fnsubext$(df$)

'open "i",1,df$

'-- estimate nessum
mbytesta&=nesmax*11.25+83 'bytes/station max
lenfile&=lof(1)
nestrec&=int(lenfile&/mbytesta&)*nesmax + (lenfile& mod mbytesta&)/11.25  'conservative estimate of # of data points
if nestrec&>32767& then
  call scprintc(irow+3,2,fnstrnum$(nestrec&)+" estimated records in "+df$+".Truncating to 32767.")
  nestrec&=32767&
endif
nessum=nestrec&
  
'nessum=lof(1)/11.25  'very conservative estimate

'-- estimate nraymax
nraymax =200
redim itpos(40,2)
if iparflg then
  if nlyr<2 then if wthk>0 then nraymax=25 else nraymax=2
else
  nsta=nstamax
  stacode$="" 
  itpos(31,1)=ktrue 'default STF pos'n
end if

'-- recover memory
redim es(1), igrflag(1),wtim(1),wamp(1),sp(1)
redim nrays(1)
redim timout(1), ampout(1,1)

'-- try to second guess max array storage required
memleft&=fre(-1)
'****memreq&=nraymax*(12.*nsta+20.)+20000.+4.*nessum 'not counting es() and seistmp()

'memory requirements are based on dynamic array needs: 
' wtim(nraymax),wamp(nraymax),timout(nraymax),ampout(nraymax,2)
' igrflag(nsta),nrays(nsta),sp(nessum) - 45000 for difference in smaller
' code requirements of MT5IT + adjusted somewhat empirically by nessum

' and in MT5IT: tmp(2*nesmax),es0(2*nes),d(nesp,np),es3(nesp,nph1),es2(nesp)
' where nes & nesp can be up to nesmax, and np <= nparmax, nph1<=3

'This estimate does not include es(maxsyn) and seistmp(maxsyn) which can get stored in tmp files

' get estimate of safe limit of inversion parameters
ntmp=(memleft&-nraymax*10.-nsta*4.-nessum*4.)/(1.5*nesmax) -8.
if ntmp<nparmax then
  if ntmp<5 then
    call schit("Data file "+df$+" too big for memory.")
    GOTO EXITRD
  end if
  call schit("Suggested limit of inversion parameters is"+str$(ntmp)+".")
  locate irow+4
else
  ntmp=nparmax
end if

memreq&=nraymax*10.+nesmax*(8.+1.5*ntmp)+nsta*4.+nessum*4.'-5000 adjusted empirically somewhat

maxsyn=fnamin((memleft&-memreq&)/4.,nsynmax)
maxsyn=fnxlimit(maxsyn,1,nessum*2&) 'could spool all data to temp files
redim seistmp(maxsyn)

'-- read event header line
print : print "Press Q to interrupt" : print
print : print "Reading file "; df$;

line input #1, hypoc$

ista=0 : nphs=1 : nessum=0 : nessumes&=0 : nscheck=ktrue : maxnes=0  
nzp=0
okprm=ktrue

erase kphflag,nodat,ndat,dext 'reset static arrays to 0

'-- start reading seismograms
while not eof(1)
  if fninquit(kfalse) then goto EXITRD
  ista=ista+1
  print "."; 'added PZ -5/10/93

'-- read and decipher header line for seismogram
  line input #1, a$
  s$=left$(a$,4)
  if iparflg then
    okprm= (ucase$(s$)=left$(fnstanm$(ista),4))
  else
    stacode$=stacode$+s$
  end if
  iphn=fnvalue(5,1)
  instin(ista)=fnvalue(6,1)
  instn=fnvalue(7,1)
  if iparflg and inst(ista)>0 then
    okprm=(iph(ista)=iphn) and (inst(ista)=instn)
  else
    iph(ista)=iphn: inst(ista)=instn
  end if
  if not okprm then
    print "ERROR in following line:";
    print a$ 'debug line
    call schit(".PRM and .INV files do not match on station entry #"+str$(ista))
    goto EXITRD
  end if
  del(ista)=fnvalue(8,8)
  azd(ista)=fnvalue(16,8)
  xmag=fnvalue(24,8)
  tstime(ista)=fnvalue(32,8)
  wt=fnvalue(40,5)
  ipf1=fnvalue(45,5)
  ndata=fnvalue(50,5)
  dt=fnvalue(63,6)
  thtime(ista)=fnvalue(69,6)
  hpfilt(ista)=fnvalue(75,6)

  inst=inst(ista)
  iph=iph(ista)

'-- check for valid PDT
  if (not fninlim(iph,1,2)) or (not fninlim(inst,1,7)) then
    print "ERROR in following line:"
    print a$ '*****debug line
    call schit("Station "+s$+" PHASE or PDT not recognized   ")
    goto EXITRD
  end if

  instva= (instin(ista) AND 6)*4 'set flag for vel(=8) & accel(=16) instr response calib
  if instin(ista) AND 1 then
    line input #1, i$
    insfile$(ista)=i$
'-- check for instrument type
    if ifopen("i",2, insfile$(ista),(0),-1)=0 then goto EXITRD
    line input #2, i$ : i$=mid$(lcase$(i$),1,2)
    if i$="ts" then
      instin(ista)=1
    elseif i$="co" then
      instin(ista)=2
    elseif i$="pz" then
      instin(ista)=3
    else
      call schit ("STA "+s$+" PHASE "+str$(iph)+" Instrument file incorrect")
      close 2 : goto EXITRD
    end if
    close 2
    if (inst>3 and instin(ista)<3) or (inst<4 and instin(ista)=3) then
      print "STA ";s$;" PHASE ";iph;" PDT designation incorrect; changing"
      select case inst
        case 4,6,7
          inst(ista)=1
        case 5
          inst(ista)=3
        case 1,2
          inst(ista)=4
        case 3
          inst(ista)=5
      end select
    end if
  else
    instin(ista)=0
  end if
  instin(ista)=instin(ista) OR instva 'add 8 for vel or 16 for accel calib

'-- assign phase type and set flag
  if iph(ista)=1 and inst(ista)<3 then kphflag(1)=ktrue : kphtype(ista)=1
  if iph(ista)=2 and inst(ista)<3 then kphflag(2)=ktrue : kphtype(ista)=2
  if inst(ista)=3 then kphflag(3)=ktrue : kphtype(ista)=3
  if iph(ista)=1 and inst(ista)=4 then kphflag(4)=ktrue : kphtype(ista)=4
  if iph(ista)=2 and inst(ista)=4 then kphflag(5)=ktrue : kphtype(ista)=5
  if inst(ista)=5 then kphflag(6)=ktrue : kphtype(ista)=6
  if inst(ista)=7 then kphflag(7)=ktrue : kphtype(ista)=7
'-- use default DT if DT not in file
  if dt=0 then dt=phdt(kphtype(ista))

'-- read zeroes and poles for GDSN seismograms and store in ZEROPOLE
  if fningdsn(ista) and (instin(ista) AND 7)=0 then
    indexzp(ista)=nzp+1
    input #1, nzero(ista), npole(ista), calcons(ista)
    for j=1 to 2*nzero(ista) + 2*npole(ista)
      nzp=nzp+1 : input #1, z
      zeropole(nzp)=z
    next
  end if

  xmag(ista)=xmag : dt(ista)=dt
  if refmag(ista)=0 or (not iparflag) then refmag(ista)=phmag(kphtype(ista))

'-- sum phases for indexing positions of elementary seismograms
  iphes(ista)=nphs
  if fnpwave(ista) then nph=2 else nph=1
  nphs=nphs+nph

'-- number of points NES in synthetic seismograms is the lesser of data and
'    default seismogram length
  ndat(ista)=fnamin(ndata, nesmax)
  maxnes=fnamax(ndat(ista),maxnes)
  nes=phrlen(kphtype(ista))/dt(ista)

'-- if no data is read in set flag NODAT
  if ndat(ista)=0 then nodat(ista)=ktrue : ndat(ista)=nes

'-- NES is the number of points in the elementary seismograms and synthetics
  nes(ista)=fnamin (ndat(ista), nes)

'-- INDEX is the starting record -1 in the array SEISTMP storing the data
  index(ista)=nessum
  nessum=nessum + ndat(ista)

'-- INPTS is the number of points in the instrument responses
  inpts(ista)=ndat(ista)

'-- INDEXES is the starting record -1 in the array ES storing the elementary seismograms
  indexes&(ista)=nessumes&
  nessumes&=nessumes& + ndat(ista)*nph

'-- check which seismogram fills up array SEISTMP, the rest is kept in random
'    access file SEISTMP$
  if nscheck then
    nsamax=ista
    if 2&*nessum>maxsyn then nsamax=ista-1 : nscheck=kfalse
  end if

'-- set values if no parameter file was read in
  if not iparflg then
    if (not nodat(ista)) then
      wt(ista)=wt: ipf1(ista)=ipf1
      n=twindow(kphtype(ista))/dt(ista)
      ninv(ista)=fnamin(n,ndat(ista))
    end if
    ipf2(ista)=ipf1 : iflag(ista)=ktrue
    j=int(azd(ista)/24)+1
    j2=int(azd(ista)/15)+17
    while itpos(j,iph) and itpos(j2,iph) 'itpos=-1 means that plot posn taken
      j=j+1
      j2=j2+1
      if j>40 then j=1
      if j2>40 then j2=1
    wend
    if itpos(j,iph) then j=j2
    ispos(ista,1)=j
    itpos(j,iph)=ktrue
    tstar(ista)=phtstar(kphtype(ista))
  end if

  call chkbound(ista)
  if thtime(ista)<=0 then thtime(ista)=(ipf1(ista)-1)*dt(ista)

320
'-- data is input in microns, convert to cm at magnification =1
  x=1.0/(xmag*10000.)

'-- limits of data window
  i1=ipf2(ista): i2=i1+ninv(ista)-1

'-- read amplitude data and find maximum value DEXT
'  erase odat 
  'redim odat(NESMAX) 'not needed
  dext=0 : s2=0
  if not nodat(ista) then
    for j=1 to ndat(ista)
      input #1, odat: odat=odat*x
      if fninlim(j,i1,i2) then dext=fnamax(abs(odat),dext)
      s2=s2+odat*odat
      odat(j)=odat
    next

'-- maximum value of data but DEXT >= 1 micron (=1.0e-5 cm)
    dext=fnamax(dext,1.0e-5)
    dext(ista)=dext*refmag(ista)

'-- RMS of data, >=1.0e-5 cm
    n=i2-i1+1
    if n>0 then datrms(ista)=fnamax(sqr(s2*dt)*refmag(ista),1.0e-5)

'-- put data in array
    call m5synio (0,0,ista,ndat(ista),0,odat(),amax,dummy)

'-- skip extra data
    if ndata>ndat(ista) then
      for j=ndat(ista)+1 to ndata: input #1, odat : next
    end if
  end if

'-- check for end-of-file if maximum seismograms reached
  if ista=nstamax and not eof(1) then d$=displaymsg("WARNING: limiting data to first"+str$(nstamax)+" seismograms","",kfalse) : goto 5490

wend

5490  nodat=kfalse : nsta=ista

'-- get rid of old screen save file
if isfile((drive$+"screen.tmp"),(""),1)=0 then kill drive$+"screen.tmp"

'-- set all phase flags
gosub RESETFLG

EXITRD: 
erase itpos
close 1
return

'----------------

SUBINST: '-- generate instrument responses
'-- WWSSN LP and SP response
'-- instype =1 WWSSN 15-100, =2 WWSSN 30-100, =3 WWSSN SP, =4 GDSN LP, =5 GDSN SP, -7 BBP
redim sp(nessum)
call memcheck
call cline(22,1) : call cline(24,1)
call scprintc(22,1,"Press Q to interrupt")

if mkinst then
  call scprint(24,1,"Calculating instrument responses")

  for ista=1 to nsta
    if fninquit(kfalse) then return 101
    instin=instin(ista) AND 7
    i=instin*10+kphtype(ista)
    call scprint(24,40,fnstanm$(ista)+fnkphlab$(i))
'    erase ss
    redim ss(2*nesmax)

'-- read file with instrument constants or time series
    if instin>0 then
      open "i",1, insfile$(ista)
      line input #1, a$  ' a$=mid$(lcase$(a$),1,2)

'-- time series
      if instin=1 then
        input #1, dt, ximag
        if dt<>dt(ista) then
          call schit("Invalid sampling interval (DT) in file "+insfile$(ista)+"  ")
          close 1: return 101
        end if
        i=0
        while (not eof(1)) and (i<256)
          i=i+1
          input #1, s : ss(i)=s/ximag
        wend
        if i<255 then
          for j=i+1 to 255 : ss(j)=0 : next
        end if

'-- instrument constants
      elseif instin=2 then
        input #1, ts,ds,tg,dg,tgain
     
'-- poles and zeros
      elseif instin=3 then
        if indexzp(ista)=0 then indexzp(ista)=nzp+1 else nzp=indexzp(ista)
        input #1, nzero(ista), npole(ista), calcons(ista)
        for j=1 to 2*nzero(ista) + 2*npole(ista)
          nzp=nzp+1 : input #1, z
          zeropole(nzp)=z 
        next
      end if

      close 1
    end if

    iph=iph(ista) : inst=inst(ista) : tstar=tstar(ista) : dt=dt(ista)
    ndat=ndat(ista) : hpfilt=hpfilt(ista)
     insame=0

'-- check for previous station with same response and steal it
     if instin=0 then
       if ista>1 then
         for j=1 to ista-1
           if fninsame(ista,j) then insame=j : goto SUBINST1
         next
       end if
     end if

SUBINST1:

'-- take previous instrument which is the same or else calculate response

    if insame > 0 then 
       inpts=fnamin(inpts(insame),ndat(ista))
       nrec1is=index(ista) : nrec1ins=index(insame)
       for i=1 to inpts : sp(i+nrec1is)=sp(i+nrec1ins) : next
    else
      instype=inst : npinst=255 
      xmag0=xmag(ista)

'-- get instrument constants  
      if fningdsn(ista) then
        nzero=nzero(ista) : npole=npole(ista)
        kp1=indexzp(ista)-1
        for i=1 to 2*nzero : kp1=kp1+1 : zero(i)=zeropole(kp1) : next
        for i=1 to 2*npole : kp1=kp1+1 : pole(i)=zeropole(kp1) : next
        xmag0=xmag0*calcons(ista)
      end if

'-- calculate instrument response

      call m5iresp (ss(),instype,instin(ista),tstar,dt,npinst,hpfilt,xmag0)

      if fningdsn(ista) then xmag0=xmag0/calcons(ista)

'-- normalize instrument to magnification of 1
      for i=1 to inpts(ista) : ss(i)=ss(i)/xmag0 : next

'-- set i=1 at time=0 and store instrument response in array SP
      nrec1is=index(ista)
      for i=2 to inpts(ista) : sp(i+nrec1is)=ss(i-1) : next
      sp(1+nrec1is)=0
    end if
  next

'-- write responses of all instruments to RAF INST.TMP
  open drive$+"inst.tmp" as 3 len=4 : field #3, 4 as io$
  lset io$=mks$(sp(1)) : put #3, 1
  for i=2 to nessum : lset io$=mks$(sp(i)) : put #3 : next
  close 3

else
  locate 24,1 : print "Reading instrument responses ";

'-- read responses of all instruments from INST.TMP
  open drive$+"inst.tmp" as 3 len=4 : field #3, 4 as io$
  get #3, 1 : sp(1)=cvs(io$)
  for i=2 to nessum : get #3 : sp(i)=cvs(io$) : next
  close 3

end if
mkinst=kfalse
return

'--------

TAUCONV: ' convolve elementary seismograms with TAU

'-- get instrument responses first
if not mkgrn then gosub SUBINST

call cline (24,1)
call scprint(24,1,"Convolving elementary seismogram with STF element")

nrec=0

for ista=1 to nsta       '------ STATION LOOP

  dt=dt(ista) : azsta=azd(ista) : p=rapar(ista)
  if fnpwave(ista) then nph=2 else nph=1

TAUCONVR: '--- REFLECTED PHASES:
          '--- convolve file GFXI$ with unit source time function element
          '      and output results to ES

603 for j=1 to nph 

    instin=instin(ista) AND 7
    i=instin*10+kphtype(ista)
    call scprint(24,55,fnstanm$(ista)+fnkphlab$(i)) 'added PZ-5/10/94
'    erase ss
    redim ss(2*nesmax)
    nrec=nrec+1

    if (fnpwave(ista) and j=1) then sorvel=alph2 else sorvel=beta2
    toa=pii-fnasin(p*sorvel)

'-- get TAU for this phase
    call taucurve ( dt, tau0, azsta, sorvel, toa, tau, nt, sor(), ierror )
    if ierror then goto TAUCONV
    tau(j+1,ista)=tau

' --- get elementary seismogram, put in TMP
    call m5esio (1, ista, j, ndat(ista), tmp() )

' --- convolve TMP with SOR, put in SS, then move to ES
    for i=1 to ndat(ista)
      for it=1 to nt
        jj=i+it-1: ss(jj)=ss(jj)+ tmp(i)*sor(it)
    next it, i

    call m5esio (0, ista, j, ndat(ista), ss() )
  next j

604 '-- DIRECT PHASE:
' -- convolve TAU with SP

'-- check for previous station with same SP
  insame=0
  if ista>1 then
     for j=1 to ista-1
       if fninsame(j,ista) then insame=j : goto TAUCONVD
     next j
  end if

TAUCONVD:
  inpts=inpts(ista)
  if insame > 0 then 
     inpts=fnamin(inpts(insame), ndat(ista) )
     nrec1is=index(ista) : nrec1ins=index(insame)
     for i=1 to inpts : sp(i+nrec1is)=sp(i+nrec1ins) : next
  else
    redim ss(2*nesmax)

    if (fnpwave(ista)) then sorvel=alph2 else sorvel=beta2
    toa=fnasin(p*sorvel)

'-- get TAU for this phase
    call taucurve ( dt, tau0, azsta, sorvel, toa, tau, nt, sor(), ierror )
    if ierror then goto TAUCONV
    tau(1,ista)=tau

' -- convolve SP with SOR, put back in SP
    nrec1is=index(ista)
    for i=1 to inpts
      for it=1 to nt
        jj=i+it-1: ss(jj)=ss(jj)+sp(i+nrec1is)*sor(it)
    next it,i

    for i=1 to inpts : sp(i+nrec1is)=ss(i) : next
  end if

next ista
return

'-------

CRESP: '--subroutine to get crustal response and elementary seismograms

'-- check memory again
redim es(1), igrflag(1),wtim(1),wamp(1)
redim nrays(1)
redim timout(1), ampout(1,1)
nraymax=2
if nlyr>1 then
  nraymax=200
elseif wthk>0 then
  nraymax=25
end if
memleft&=fre(-1)
'-- get max array size actually used
if nsamax<nsta then maxsyn=2*(index(nsamax+1)) else maxsyn=nessum*2
maxsyn=fnamax(maxsyn,1)
if maxsyn<ubound(seistmp) and memleft&>maxsyn*2+2000 then 'compress seistmp()
  dim tmp1(maxsyn)
  for i=1 to maxsyn: tmp1(i)=seistmp(i) :next
  redim seistmp(maxsyn)
  for i=1 to maxsyn: seistmp(i)=tmp1(i) :next
  erase tmp1
end if

'****memreq&=nraymax*(12.*nsta+20.)+13000.+(nessum-ubound(seistmp))*4. 'w/all data spooled to disk

' memory requirements are same as estimated above except that we now know size of
' seistmp() (array which can be reduced to 0 if needed) and sp(nessum) is allocated

memreq&=nraymax*10.+nesmax*(8.+1.5*ntmp)+nsta*4.-maxsyn*2.' -10000. not counting es() & seistmp()
d2$="Kb required, only"+str$(memleft&\1024.)+" remaining."
if memleft&<memreq& then
  d$=str$(memreq&\1024.)
  d$=displaymsg$("No memory left for data -"+d$+d2$,"",kfalse):return 101
else
  if memleft&<memreq&+(maxsyn)*4& then
    d$=str$((memreq&+(maxsyn)*4&)\1024.)
    for i=22 to 25:call cline(i,1):next :call setmsgcol
    call scprint(22,1,"Too much data for memory -"+d$+d2$)
    call scprint(23,1,"Memory can be conserved by spooling more data to .TMP files (slower)")
    call settxtcol
    d$=displaymsg$("OK to spool more data to .TMP files","YN",kfalse)
    if d$="Y" then
      nsynmax=fnamax((memleft&-memreq&)/8.,1.)
      d$=displaymsg$("Select R to read data again","",kfalse):w$="R"
    end if
    return 101
  end if
  d$=""
end if

redim es(maxsyn)
redim wtim(nraymax), wamp(nraymax), igrflag(nsta), nrays(nsta)
call memcheck
if nsta>nsamax then call cline(23,1): call scprint(23,1,"Data from"+str$(nsta-nsamax)+_
 " of"+str$(nsta)+" stations being stored in .TMP files")

nrec=0

'-- get velocities of halfspace(2) and layer just above halfspace(1)
if nlyr=1 then
  alph2=alph(1): beta2=beta(1): rho2=rho(1)
  if wthk>0 then
    alph1=alphw: beta1=betaw: rho1=rhow   'water
  else
    alph1=alphf: beta1=betaf: rho1=rhof   'air
  end if
  sthk=0
else
  alph1=alph(1): beta1=beta(1): rho1=rho(1)
  alph2=alph(2): beta2=beta(2): rho2=rho(2)
  sthk=ztop(2)-ztop(1)
end if
vps=alph2/beta2

'-- get elementary seismogram, put in ES
call cline(24,1)
if not mkgrn then  'load ES from file

  call scprint(24,1,"Reading elementary seismograms ")
  for ista=1 to nsta
    if fninquit(kfalse) then return 101
    if fnpwave(ista) then nph=2 else nph=1

'-- get ray parameter and take-off angle
    call m5rapars (deptot, del(ista), iph(ista),alph2,beta2,rapar(ista),ptoa(ista))

    for j=1 to nph
      nrec=nrec+1
      call m5tmpio (1, ista, j, ndat(ista), tmp(), gfxi$ )
      call m5esio (0, ista, j, ndat(ista), tmp() )
    next
  next
  return
end if

call scprint(24,1,"Computing elementary seismograms ")

for ista=1 to nsta
  if fninquit(kfalse) then return 101

  redim timout(nraymax), ampout(nraymax,2)
  nrays=0

  iph=iph(ista) : inst=inst(ista) : dt=dt(ista) : nes=nes(ista)
  locate 24,40: print "Station ";fnstanm$(ista);" phase ";iph;" inst ";inst(ista);
  if fnpwave(ista) then nph=2 else nph=1        'SH has one upgoing phase

'-- get ray parameter and take-off angle
  call m5rapars (deptot, del(ista), iph(ista),alph2,beta2,rapar(ista),ptoa(ista))

  rp=rapar(ista)           'ray parameter
  igrflag(ista)=ktrue      'flag that the response is calculated, not borrowed
  timax=(ndat(ista)-1)*dt  'maximum arrival time for any phase

 if fnpwave(ista) then

'-- get reflection and transmission coefficients, and angles at top of halfspace

  call psvcoef (rp, alph1, beta1, rho1, alph2, beta2, rho2, refcs(), aincs(),r0)

  nrays=1
  timout(nrays)=0.0
  ampout(nrays,1)=refcs(3,3)   'pP arrival from base
  ampout(nrays,2)=refcs(4,3)   'sP arrival

  if nlyr=1 and wthk>0 then     'get water reverbs
    call water (rp, alphw, wthk, refcs(1,1), r0, timax, wtim(), wamp(), nwater)

    pupu=refcs(3,1) : pdpd=refcs(1,3) : supu=refcs(4,1)
    for i=1 to nwater
      nrays=fnamin(nrays+1,nraymax)
      timout(nrays)=wtim(i)
      ampout(nrays,1)= pupu*pdpd*wamp(i)
      ampout(nrays,2)= supu*pdpd*wamp(i)
    next
  end if

  if nlyr>1 then

'-- save coefficients at bottom of layer 1
    pupu1=refcs(3,1) : pusu1=refcs(3,2): supu1=refcs(4,1)
    susu1=refcs(4,2) : pdsu1=refcs(1,2): pdpd1=refcs(1,3)
    sdpd1=refcs(2,3) : pdpu1=refcs(1,1): sdsu1=refcs(2,2)
    sdpu1=refcs(2,1)

'-- get coefficients at top of layer 1
    if wthk>0 then
      vpu=alphw : vsu=betaw : rhou=rhow
    else
      vpu=alphf : vsu=betaf : rhou=rhof
    end if

    call psvcoef (rp, vpu, vsu, rhou, alph1, beta1, rho1, refcs(), aincs(), r0)

    ai2=aincs(3) : aj2=aincs(4)

    pdpu=refcs(1,1) : pupd=refcs(3,3) : pupu=refcs(3,1) : pdpd=refcs(1,3)
    pusd=refcs(3,4) : supd=refcs(4,3) : pdsd=refcs(1,4) : supu=refcs(4,1)
    susd=refcs(4,4)

    if wthk>0 then call water (rp, alphw, wthk, pdpu, r0, timax, wtim(), wamp(), nwater)

'-- time correction to wavefront for each bounce in upper layer

    rxs=sthk*tan(aj2) : rxp=sthk*tan(ai2)
    ppc= rxp*rp : ssc=rxs*rp
    tp1=sthk/(alph1*cos(ai2))-ppc  'delays in upper layer, one-way
    ts1=sthk/(beta1*cos(aj2))-ssc

'-- reverberations in upper layer, for NREV between 2 and NREVCOM call the
'    routine M5CRESP that gets all the rays; Otherwise, use this approximation
'   Note: the following routine is exact for 1 bounce and it calculates
'     water reverberations

    nrevmax=fnamin(int(timax/(2*tp1))+1, 20)      'maximum reverbs
    nrev1 = fnamin(nrevmax,nrevcom)

    if nrev1 >1 then _
     call m5cresp (sthk, rp, alph1, beta1, rho1, vpu, vsu, rhou, alph2, beta2, rho2, nrev1, timout(), ampout(), nrays)

    for n=1 to nrevmax
      if fninquit(kfalse) then return 101
      if n=1 or n>nrev1 then

        tpp= n*2*tp1  'PP path
        if tpp>timax then goto 7018
        nrays=fnamin(nrays+1,nraymax)
        timout(nrays) = tpp
        ppn=pupd^n*pdpu1^(n-1)
        ppp=pupu1*ppn*pdpd1 : ampout(nrays,1)= ppp
        spp=supu1*ppn*pdpd1 : ampout(nrays,2)= spp
'water
        if n=1 and wthk>0 then
          p1=ppp/pupd * pupu*pdpd
          s1=spp/pupd * pupu*pdpd
          for i=1 to nwater
            nrays=fnamin(nrays+1,nraymax)
            timout(nrays)=tpp + wtim(i)
            ampout(nrays,1)= p1 * wamp(i)
            ampout(nrays,2)= s1 * wamp(i)
          next
        end if

        tsp=n*(ts1+tp1)  'SP and PS paths
        if tsp>timax then goto 7018
        nrays=fnamin(nrays+1,nraymax)
        timout(nrays) = tsp
        psn=pusd^n*sdpu1^(n-1): spn=supd^n*pdsu1^(n-1)
        ppsps=pupu1*psn*sdpd1 : pspsp=pusu1*spn*pdpd1
        sspsp=susu1*spn*pdpd1 : spsps=supu1*psn*sdpd1
        ampout(nrays,1) = ppsps + pspsp
        ampout(nrays,2) = sspsp + spsps

'water
        if n=1 and wthk>0 then
          p1=(ppsps/pusd*pupu*pdsd+pspsp/supd*supu*pdpd)
          s1=(sspsp/supd*supu*pdpd+spsps/pusd*pupu*pdsd)
          for i=1 to nwater
            nrays=fnamin(nrays+1,nraymax)
            timout(nrays)=tsp + wtim(i)
            ampout(nrays,1)=p1*wamp(i)
            ampout(nrays,2)=s1*wamp(i)
          next
        end if

        tss=n*2*ts1  'SS path
        if tss>timax then goto 7018
        nrays=fnamin(nrays+1,nraymax)
        timout(nrays) = tss
        ssn=susd^n*sdsu1^(n-1)
        pss=pusu1*ssn*sdpd1 : ampout(nrays,1) = pss
        sss=susu1*ssn*sdpd1 : ampout(nrays,2) = sss
'water
        if n=1 and wthk>0 then
          p1=pss/susd * supu*pdpd
          s1=sss/susd * supu*pdpd
          for i=1 to nwater
            nrays=fnamin(nrays+1,nraymax)
            timout(nrays)=tss + wtim(i)
            ampout(nrays,1)= p1* wamp(i)
            ampout(nrays,2)= s1* wamp(i)
          next
        end if
      end if
7018  next

  end if

 elseif (not fnpwave(ista)) then            'get ES for SH phase

  call shcoef (rp, alph1, beta1, rho1, alph2, beta2, rho2, refcs(), aincs())
  nrec1=nrec+1
  nrays=1
  timout(nrays) = 0.0
  ampout(nrays,1) = refcs(2,2)    'sS phase

  if nlyr>1 then

    aj1=aincs(1)
    susu=refcs(2,1) : sdsu=refcs(1,1) : sdsd=refcs(1,2)

    ts1=2.*sthk/(beta1*cos(aj1))       'delay in upper layer
    tc=2.*sthk*tan(aj1)*rp             'wavefront correction
    ts1=ts1-tc

    sref=1   'reflection from the seafloor
    for n=1 to ndat(ista)
      tss = n*ts1
      if tss>timax then 7330
      nrays=fnamin(nrays+1,nraymax)
      timout(nrays) = tss
      ampout(nrays,1) = susu*sref^n*sdsu^(n-1)*sdsd
    next
  end if

 end if

7330   '--- convolve series TIMOUT,AMPOUT with instrument SP and store in ES
 if fninquit(kfalse) then return 101

'put SP in TMP
 nrec1is=index(ista) : for i=1 to inpts(ista) : tmp(i)=sp(i+nrec1is) : next

 for j=1 to nph
  nrec=nrec+1
'  erase ss 
  redim ss(2*nesmax) : esm=0   'max of ES
  for i=1 to nrays: esm=fnamax(abs(ampout(i,j)),esm) : next
  esm=0.01*esm     ' ignore if <1% of max of ES

  for i=1 to nrays
    esamp=ampout(i,j)
    if abs(esamp) >= esm then
      t=timout(i)/dt +1            '+1 added to be the same as before
      itt=int(t)                   'so that i=1 at t=0
      lim=fnamin(ndat(ista)-itt-1,inpts(ista))
      if lim>0 then

'-- interpolate ES so it arrives at exact time
        dd=t-itt : omdd=1.0-dd

        for it=1 to lim
          k=itt+it
          ss(k)=ss(k)+ esamp* ( tmp(it)*dd + tmp(it+1)*omdd )
        next it

      end if
    end if
  next i

 ' -- output response ES to files
  call m5tmpio (0, ista, j, ndat(ista), ss(), gfxi$ )
  call m5esio (0, ista, j, ndat(ista), ss() )

 next j

'-- write times and amplitudes to AMPS.TMP
 if ista=1 then open "o",9,drive$+"amps.tmp" else open "a",9,drive$+"amps.tmp"
 print #9, "Station "; fnstanm$(ista); "  Phase "; iph(ista); " Inst ";inst(ista)
 print #9, nrays
 nrays(ista)=nrays
 for i=1 to nrays
  print #9, using fnfform$(12,3)+"  "; timout(i);
  print #9, using fndform$(8,5)+" "; ampout(i,1);ampout(i,2)
 next
 close 9

7420 next ista
return

'---------------------------------------------------------------------------
sub psvcoef (rp, alph1, beta1, rho1, alph2, beta2, rho2, refcs(2), aincs(1), r0) static

'-- get reflection and transmission coefficients for P-SV
'--  (Aki & Richards p. 149-151)
'
' notation in REF(i,j) and AIN(i)
'                i
'       1     2     3     4
'   1  pdpu, sdpu, pupu, supu
'j  2  pdsu, sdsu, pusu, susu
'   3  pdpd, sdpd, pupd, supd
'   4  pdsd, sdsd, pusd, susd
'
' aincs() ai1,  aj1,  ai2,  aj2

p=rp : p2=rp*rp    'ray parameter

ai1=fnasin(p*alph1)  : aincs(1)=ai1
aj1=fnasin(p*beta1)  : aincs(2)=aj1
ai2=fnasin(p*alph2)  : aincs(3)=ai2
aj2=fnasin(p*beta2)  : aincs(4)=aj2

ci1=cos(ai1)/alph1 : ci2=cos(ai2)/alph2
cj1=cos(aj1)/beta1 : cj2=cos(aj2)/beta2
a=rho2*(1-2*beta2^2*p2)-rho1*(1-2*beta1^2*p2)
b=rho2*(1-2*beta2^2*p2)+2*rho1*beta1^2*p2
c=rho1*(1-2*beta1^2*p2)+2*rho2*beta2^2*p2
d=2*(rho2*beta2^2-rho1*beta1^2)
e=b*ci1+c*ci2
f=b*cj1+c*cj2
g=a-d*ci1*cj2
h=a-d*ci2*cj1
dd=e*f+g*h*p2

a1=(a*c+b*d*ci1*cj1)
a2=(b*ci1-c*ci2)*f
a3=(a+d*ci2*cj1)*p2*g
a4=(a+d*ci1*cj2)*h*p2
a5=(a*b+c*d*ci2*cj2)
a6=(b*cj1-c*cj2)*e

pdpu=(a2-a4)/dd
pdsu=-2*ci1*a5*p*alph1/(beta1*dd)
pdpd=2*rho1*ci1*f*alph1/(alph2*dd)
pdsd=2*rho1*ci1*h*p*alph1/(beta2*dd)
sdpu=-2*cj1*a5*p*beta1/(alph1*dd)
sdsu=-(a6-a3)/dd
sdpd=-2*rho1*cj1*g*p*beta1/(alph2*dd)
sdsd=2*rho1*cj1*e*beta1/(beta2*dd)
pupu=2*rho2*ci2*f*alph2/(alph1*dd)
pusu=-2*rho2*ci2*g*p*alph2/(beta1*dd)
pupd=-(a2+a3)/dd
pusd=2*ci2*a1*p*alph2/(beta2*dd)
supu=2*rho2*cj2*h*p*beta2/(alph1*dd)
susu=2*rho2*cj2*e*beta2/(beta1*dd)
supd=2*cj2*a1*p*beta2/(alph2*dd)
susd=(a6+a4)/dd

rs0=4*p2*(ci1)*(cj1)
rs1=(1/beta1^2-2*p2)^2
r0=(rs0-rs1)/(rs0+rs1)

refcs(1,1)=pdpu : refcs(2,1)=sdpu : refcs(3,1)=pupu : refcs(4,1)=supu
refcs(1,2)=pdsu : refcs(2,2)=sdsu : refcs(3,2)=pusu : refcs(4,2)=susu
refcs(1,3)=pdpd : refcs(2,3)=sdpd : refcs(3,3)=pupd : refcs(4,3)=supd
refcs(1,4)=pdsd : refcs(2,4)=sdsd : refcs(3,4)=pusd : refcs(4,4)=susd

end sub
'---------------------------------------------------------------------------
sub shcoef (rp, alph1, beta1, rho1, alph2, beta2, rho2, refcs(2), aincs(1)) static

'-- get reflection and transmission coefficients for SH waves  (A&R p. 144)

'
' notation in REF(i,j) and AIN(i)
'          i
'       1     2
'   1  sdsu, susu
'j  2  sdsd, susd
'
' aincs() aj1,  aj2

aj1=fnasin(rp*beta1)  : aincs(1)=aj1
aj2=fnasin(rp*beta2)  : aincs(2)=aj2

dd=rho1*beta1*cos(aj1)+rho2*beta2*cos(aj2)
sdsu=(rho1*beta1*cos(aj1)-rho2*beta2*cos(aj2))/dd
susu=2*rho2*beta2*cos(aj2)/dd
sdsd=2*rho1*beta1*cos(aj1)/dd
susd=-sdsu

refcs(1,1)=sdsu : refcs(2,1)=susu
refcs(1,2)=sdsd : refcs(2,2)=susd

end sub
'===============================================================
sub m5cresp (sthk, p, vp, vs, rho, vpu, vsu, rhou, vpd, vsd, rhod, nrev, timout(1), ampout(2), nrays) static

'-- calculate complete P-SV response of a single layer for NREV bounces

dim ampin0(2,2), ampout0(2), a(4,2), refcs(4,4), aincs(4)
shared nraymax

'-- reflection coefficients at upper surface
call psvcoef (p, vpu, vsu, rhou, vp, vs, rho, refcs(), aincs(), r0)
a(1,1)=refcs(3,3) : a(2,1)=refcs(3,4) : a(3,1)=refcs(4,3) : a(4,1)=refcs(4,4)

'-- reflection coefficients at lower surface
call psvcoef (p, vp, vs, rho, vpd, vsd, rhod, refcs(), aincs(), r0)
a(1,2)=refcs(1,1) : a(2,2)=refcs(1,2) : a(3,2)=refcs(2,1) : a(4,2)=refcs(2,2)

ai1=aincs(1) : aj1=aincs(2)
xp=sthk*tan(ai1) : xs=sthk*tan(aj1)

tp=sthk/(vp*cos(ai1))-xp*p : ts=sthk/(vs*cos(aj1))-xs*p

'--- incoming transmission coefficients at lower surface
ampin0(1,1)=refcs(3,1) : ampin0(1,2)=refcs(3,2)
ampin0(2,1)=refcs(4,1) : ampin0(2,2)=refcs(4,2)

'--- outgoing transmission coefficients at lower surface
ampout0(1)=refcs(1,3)
ampout0(2)=refcs(2,3)

for irev=2 to nrev

  nrays=fnamin(nrays+1,nraymax)
  nlegs=2*irev  'number of legs, times, and amplitudes
  nray=2^nlegs

  if (nrays+nlegs)>nraymax then exit sub

  for k=0 to nlegs : timout(k+nrays)=k*ts + (nlegs-k)*tp : next

  for jray=0 to nray-1
    call binconv (jray, nlegs, bin$)  'get binary value for ray

'-- KSUM is the number of S phases, Time = ksum*Ts+(nlegs-ksum)*Tp 
    amp=1 : mm= ktrue   
    ksum=0 
    for i=1 to nlegs: ksum=ksum+val(mid$(bin$,i,1)) : next

' -- look at each bounce, KRAY1 is the phase going into the bounce, KRAY2 is 
'       the phase coming out

    for i=1 to nlegs-1
      kray1=val(mid$(bin$,i,1))

      if (i=1) then kr1=kray1+1         'phase of first leg in layer

      kray2=val(mid$(bin$,i+1,1))

      if (i=nlegs-1) then kr9=kray2+1    'phase of last leg in layer

      kray=2*kray1+kray2               'KRAY=0(PP) =1(PS) =2(SP) =3(SS)
      amp=amp*a(kray+1,mm+2)           'find coefficient
      mm=not mm                        'MM=true for top interface, false for bottom
    next i
    k2ray = ksum+nrays

    ampout(k2ray,1)=ampout(k2ray,1) + amp * ampin0(1,kr1) * ampout0(kr9)   'pP
    ampout(k2ray,2)=ampout(k2ray,2) + amp * ampin0(2,kr1) * ampout0(kr9)   'sP
  next jray 

  nrays=nrays+nlegs
next irev
end sub
'=====
sub binconv (nint, nlen, bin$) static
'--- convert integer NINT into binary string BIN$ of length NLEN

binoct$="000001010011100101110111"

bin$="" : o$=oct$(nint) : l=len(o$)

for k=1 to l
  bin$=bin$+mid$(binoct$, 1+3*val(mid$(o$,k,1)), 3)
next

lbin=len(bin$)-nlen

if lbin<0 then 
  for k=1 to abs(lbin) : bin$="0"+bin$ : next
elseif lbin>0 then
  bin$=mid$(bin$,lbin+1,nlen)
end if
end sub
'=======
sub water(rp, alphw, wthk, pdpu, r0, timax, t(1), amp(1), nrev ) static
'--------- get reverberations in water layer
shared nraymax

if wthk=0 then exit sub
cutoff = 1.0e-5
ai1=fnasin(rp*alphw)  'angle in water

wdly = 2*wthk/(alphw*cos(ai1))     'two-way delay in water
twc =  2*wthk * tan(ai1)*rp        'correction to wavefront
wdly=wdly-twc

nrev = int(timax/wdly)
nrev=fnamin(nrev,nraymax)
pw=1 : n=0

while abs(pw) > cutoff and n < nrev
  n=n+1
  pw=r0^n*pdpu^(n-1)
  t(n)= n*wdly
  amp(n)=pw
wend

nrev=n
end sub
'======
sub CHKBOUND (ista) static
'--------- check bounds on NINV, IPF1 and IPF2
shared ipf1(),ipf2(),ninv(),ndat(),nes()

ipf1(ista)=fnamax(1,ipf1(ista))
ipf1(ista)=fnamin(ipf1(ista),ndat(ista))
ipf2(ista)=fnamax(1,ipf2(ista))
ipf2(ista)=fnamin(ipf2(ista),ndat(ista))
ninv(ista)=fnamin(ninv(ista),ndat(ista)-ipf2(ista))
ninv(ista)=fnamax(ninv(ista),0)
ninv(ista)=fnamin(ninv(ista),nes(ista))
end sub
'======
sub taucurve ( dt, tau0, azsta, sorvel, toa, tau, nt, sor(1), ierror ) static
shared vrupt, azrupt, plrupt
ierror=0

'-- get TAU from directivity
if vrupt=0 then
  tau=tau0
else

  ar=fnrad(azrupt) : dr=fnrad(plrupt) : az=fnrad(azsta)
'-- cos of angle between ray and rupture
  cosomega= cos(dr)*sin(ar)*sin(toa)*sin(az) +_
  cos(dr)*cos(ar)*sin(toa)*cos(az) +_
  sin(dr)*cos(toa)

  tau=tau0*(1.0-vrupt/sorvel*cosomega)

  if tau <=0 then
    call cline(23,1) : call cline(24,1)
    locate 23,1: print "VRUPT is too large, TAU is negative"
    d$=str$(vrupt)+str$(azrupt)+str$(plrupt)
    do
      if inputn(24,1,"Enter new rupture velocity,azimuth,plunge ",d$,ktrue,ktrue,ktrue,"a",24)=27 then exit do
      s1$=lftvar$(d$):s2$=lftvar$(d$):s3$=lftvar$(d$)
      if istype(s1$,"f") and istype(s2$,"f") and istype(s3$,"f") then
        vrupt=val(s1$): azrupt=val(s2$): plrupt=val(s3$): exit do
      else
        d$=displaymsg$("Entry error","",kfalse)
      end if
    loop
    ierror=-1 : exit sub
  end if

end if

'-- set up array SOR with the TAU triangle amplitudes at DT increments
'-- if TAU<DT then a delta-function is used for the triangle
'-- if isvel flag is set, convert SOR to approximate a box function 
'     to provide source velocity = dx/dt where d(t) is triangle function
'     (added to allow use of broadband inst response for SEED data
'     PZ - 94/5/13)


nt=2*tau/dt + 1 : h=dt/tau
if tau<=dt then nt=3 : h=1.0
if ubound(sor)<nt then redim sor(nt) 'added PZ -5/10/94
for i=1 to nt
  v=(i-1)*h
  if v<=1 then sor=v else sor=2-v
  sor(i)=fnamax(sor,0)*h
next i
end sub

'---------------------------------------------------------------------------
CHGVALUES:
'-- change t-star or ref mag for individual stations

if ich then
  m1$="Change Attenuation t*" : m2$="      t*"
  for i=1 to nsta : sor(i)=tstar(i) : next
  frm$=fnfform$(8,2) : dval=0.1
else
  m1$="Change Reference Magnification" : m2$="  RefMag"
  for i=1 to nsta : sor(i)=refmag(i) : next
  frm$=fniform$(8) : dval=1000.
end if

call scbox(1,22,ktrue,m1$,"Press Esc to exit")
i=0

for icol=2 to 58 step 31
  call scprint(2,icol,"Sta Flag Del Phase"+m2$)
  for k=1 to 19
    i=i+1 : if i>nsta then 402
    call settxtcol
    irow=k+2
    call sclocate(irow, icol) 
    print using fnaform$(5); fnstanm$(i);
    print using fniform$(2)+" "; -iflag(i);
    print using fniform$(4)+" "; del(i);
    instin=instin(i) AND 7
    it=instin*10+kphtype(i)
    print using " "+fnaform$(4); fnkphlab$(it);
    call setvarcol
    print using frm$; sor(i)
  next
next

402 irow=1 : icol=1
403 'get new value
call scprintc(irow+2, 20+(icol-1)*31, ">")
call scprintc(23,1,"KEYS:   Up   Down (or arrows) C-change (or +/-)")
call inkey(d$)
call scprintv(irow+2, 20+(icol-1)*31," ")

if fniesc(d$) then 
 if ich then
  for i=1 to nsta : tstar(i)=sor(i) : next
 else
  for i=1 to nsta : refmag(i)=sor(i) : next
 end if
   return
 end if

if fnmoveup(d$) or d$="u" then irow=irow-1 
if fnmovedn(d$) or d$="d" then irow=irow+1
if d$="c" or fniplus(d$) or fnminus(d$) then
  i= (icol-1)*19 + irow : valin=sor(i)
  if d$="c" then
    call cline(23,1):d$=str$(sor(i))
    if inputn(23,1,"New value ",d$,ktrue,ktrue,ktrue,"f",6)=27 then 403
    valin=val(d$):if valin<0 then 403
  elseif fniplus(d$) then 
    valin=valin+dval
  elseif fnminus(d$) then 
    valin=valin-dval
  end if
  valin = fnamax(valin,0.0)
  if sor(i)<>valin then 
    sor(i)=valin  : if ich then mkinst=ktrue
    call sclocate(irow+2, 20+(icol-1)*31) 
    call setvarcol
    print using frm$; valin; 
    call settxtcol
  end if
end if
408  irow=fnamax(irow,1)
if irow>19 then irow=1 : icol=icol+1
if (icol-1)*19+irow > nsta then irow=1 : icol=1
goto 403 
return

'======
rem $INCLUDE: 'm5rdparm'
'======
