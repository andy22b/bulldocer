' main MT5INT
'  interpolation and windowing to set up data file for waveform inversion
'  data is input in seconds and microns from the original seismogram
'  predicted arrival times are based on J-B tables

' Written by Rob McCaffrey and Geoff Abers, MIT, 1988
' Rewritten Fall, 1990 by Peter Zwick for SYN4
' modified 28/Mar/1992 by Peter Zwick for MT4
' corrected npts(sta) when nn truncated and made indexs long, 5/Nov/1993  PZ
' increased NINTMAX0 9/Nov/1993 PZ
'
' IDATYPE = 1 for WWSSN, 2=GDSN

DEFINT I-N: OPTION BASE 1
CONST NSTAMAX=1 'dummy variable

rem $INCLUDE: 'M5COMSCR'
rem $INCLUDE: 'm5fndecl'
rem $INCLUDE: 'm5functn'

'--- welcome
prgnm$="MT5INT"
authors$="P. Zwick, R. McCaffrey, G. Abers"
version$="1.0"
call welcome
call scprint(13,3,"Setting up..")

'-- NSEISMAX=maximum number of seismograms
nseismax=100
dim npts(nseismax), tfirstpt(nseismax), tlastpt(nseismax)
dim del(nseismax), az(nseismax), baz(nseismax), flag$(nseismax)
dim iord(nseismax), inst(nseismax), indexdat&(nseismax), idatype(nseismax)
dim dt(nseismax), ttime(nseismax), xmag(nseismax), kphtyp(nseismax)
dim npoles(nseismax), nzeros(nseismax), calconst(nseismax)

'dim zeroi(100,2), polei(100,2)
dim zeroi(nseismax,2), polei(nseismax,2)
dim staphase$(5), stacode$(5)
dim instva(nseismax) '8=vel type, 16=accel type calibration


nphtype=7
dim tkept(7), tlead(7), dtdef(7), fcuts(7)

dim phin$(5), datin$(2), smag(2), istanum(2)
dim tmin(2), tmax(2), iwindow(2)

'-- maximum number of stations in file is MAXWWSTA
maxwwsta=200
dim wwscode$(6), adlat(maxwwsta), adlon(maxwwsta)

'-- maximum number of interpolated points is NINTMAX0
nintmax0=3072  'revised from 1028 Apr 15,92, and from 1024 Nov 9,1993
dim yout(nintmax0), yout2(nintmax0,2)
dim work(nintmax0)

'-- maximum number of digitized time-amplitude points is INMAX
'inmax=500
inmax=nintmax0
dim txdat(inmax), ydat(inmax)

'-- travel time tables, set up for 2 phases, 120 deltas, 15 depths
dim ttimtab(120,15,2), dep(15,2), vsurf(2)
dim del0(2), ddelt(2), ndel(2), ndep(2)

'-- functions for efficient storage of string data
DEF FNWWSNM$(I)
  ijt=int(i/50.1) +1
  ijs= ((i-1) mod 50)*4 +1
  fnwwsnm$=mid$(WWSCODE$(ijt),ijs,4)
end def

DEF FNSTANM$(I)
  ijt=int(i/50.1) +1
  ijs= ( (i-1) mod 50 )*4 +1
  fnstanm$=mid$(STACODE$(ijt),ijs,4)
end def

DEF FNSTAPH$(I)
  ijt=int(i/50.1) +1
  ijs= ( (i-1) mod 50 )*3 +1
  fnstaph$=mid$(STAPHASE$(ijt),ijs,3)
end def

'-- plot functions
def fntime(i)= tmin + (i-1)*dt
def fnix(t)=kscrx*(t-tmin)/(tmax-tmin)
def fniy(y)=(iyzero-0.8*yamp*y/ymax)
def fnia(a)=(ix-ssc*a/ymax)
def fnit(t)=t*xtsc+noff

'-- function keys
def fnkeyf1(x$) = (x$=chr$(59))
def fnkeyf2(x$) = (x$=chr$(60))
def fnkeyf3(x$) = (x$=chr$(61))
def fnkeyf4(x$) = (x$=chr$(62))
def fnkeyf5(x$) = (x$=chr$(63))
def fnkeyf6(x$) = (x$=chr$(64))
def fnkeyf7(x$) = (x$=chr$(65))

'--- phase and data types
DEF FNPWAVE = (kphtype=1) OR (kphtype=3) OR (kphtype=4) OR (kphtype=6) OR (kphtype=7)
DEF FNILONG = (kphtype=1) OR (kphtype=2) OR (kphtype=4) OR (kphtype=5)
DEF FNISHORT = (kphtype=3) OR (kphtype=6) OR (kphtype=7)
DEF FNINGDSN = (kphtype>3)
DEF FNIWWSSN = (kphtype<4)
DEF FNKPHLAB$(I) = mid$("WWLPWWLSWWSPGDLPGDLSGDSPGDBP",i*4-3,4)

'--- time lengths to be kept
  tkept(1)=60. : tkept(2)=80. : tkept(3)=15. : tkept(4)=100.
  tkept(5)=120. : tkept(6)=15. : tkept(7)=60.
'--- lead times
  tlead(1)=10. : tlead(2)=15. : tlead(3)=3. : tlead(4)=15.
  tlead(5)=20. : tlead(6)=3. : tlead(7)=10.
'--- default DT
  dtdef(1)=0.5 : dtdef(2)=1.0 : dtdef(3)=0.1 : dtdef(4)=1.
  dtdef(5)=1.0 : dtdef(6)=0.05 : dtdef(7)=0.2
'--- cutoff period
  fcuts(1)=70. : fcuts(2)=70. : fcuts(3)=20. : fcuts(4)=70.
  fcuts(5)=70. : fcuts(6)=20. : fcuts(7) = 66.67

'--- default fileheader
fhdr$="No-File"
data "lpz","lpn","lpe","spz","bbz"
for i=1 to 5 : read phin$(i): next

data "WWSSN","GDSN"
for i=1 to 2 : read datin$(i) : next

littic=kscry/10 : nstep=1

'-- read in stations
   gosub READSTA

'-- read travel time table
numttype=2
dim ttfilenm$(numttype),ttname$(numttype)
data "m5j-b.dat","m5iasp91.dat"
data "Jeffreys-Bullen (used for EQ loc'n by ISC)","IASPEI 91"
for i=1 to 2 : read ttfilenm$(i) : next
for i=1 to 2 : read ttname$(i) : next
  ittype=1
  ttfile$=syndir$+ttfilenm$(ittype)
  call readtime (ttfile$)

  call scpause(18,3)

100 call setscreen(0)
    inord=0

call scbox(1,9,ktrue,"STARTUP MENU","Press Q to quit")
call scprint(2,2,"KEY   OPTION")
call scprompt(4,3,"1 - Process new data from .DIG and .DSN files (create .INV file)")
call scprompt(5,3,"2 - Plot interpolated data on screen from .INV file")
call scprompt(6,3,"3 - Plot digitized data on screen from .DIG and .DSN files")
call scprompt(7,3,"4 - Change travel time table from "+ttname$(ittype))
call scprompt(8,3,"X/"):call scprompt(8,5,"Q - Shell/Quit")

a$=""
101   ik=input1(10,1,"Enter KEY for option ",a$,kfalse,ktrue,ktrue,numlist$(1,4)+"xq")
     if a$="q" then 
       call endcheck:goto 101
     elseif a$="x" then 
       call shellcom
     elseif a$="1" then
       goto 250
     elseif a$="2" then
       gosub PLOTSEIS1
     elseif a$="3" then
       gosub PLOTRAW
     elseif a$="4" then
       call scprint(10,3,"KEY  TRAVEL TIME TABLE SELECTION")
       b$=fnstrnum$(ittype)
       for i=1 to numttype
         call scprompt(10+i,3,fnstrnum$(i)+" - "+ttname$(i))
       next
       if input1(10+i,1,"Enter KEY for selection ",b$,kfalse,ktrue,ktrue,numlist$(1,numttype))=27 then goto 100
       ia=val(b$)
       if ia<>ittype then
         ittype=ia
         ttfile$=syndir$+ttfilenm$(ittype)
         call readtime (ttfile$)
       end if
     end if
     goto 100

250 call scclear
    call setscreen(0)
    call scbox(1,17,ktrue,"PARAMETER SETUP MENU","Press Esc to exit")

call scprint(2,2,"KEY   OPTION")
call scprompt(4,3,"I - Input file header   "):call scprintv((0),(0),"<"+ fhdr$+">")
call scprompt(5,3,"O - Output file         "):call scprintv((0),(0),invfile$)
    if newopen then call scprintv((0),(0)," as NEW file") else call scprintv((0),(0)," APPENDING")

call scprint(7,7, "parameter / data type ")
   for i=1 to nphtype
      print "  ";fnkphlab$(i);" ";
   next

call scprompt(8,3,"S - Seismogram length, s ")
   call setvarcol
   for i=1 to nphtype
      print using fnfform$(7,1);tkept(i);
   next
   call settxtcol

call scprompt(9,3,"L - Lead time, s         ")
   call setvarcol
   for i=1 to nphtype
      print using fnfform$(7,1); tlead(i);
   next
   call settxtcol

call scprompt(10,3,"D - DT, s                ")
   call setvarcol
   for i=1 to nphtype
      print using fnfform$(7,2);dtdef(i);
   next
   call settxtcol

call scprompt(11,3,"F - Filter cutoff period ")
   call setvarcol
   for i=1 to nphtype
      print using fnfform$(7,1);fcuts(i);
   next
   call settxtcol

call scprompt(13,3,"X - Shell to DOS")
call scprompt(14,3,"B - Begin")

    a$=""
    if input1(18,1,"Enter KEY for option ",a$,kfalse,ktrue,ktrue,"iosldfxb")=27 then 100

if a$="i" then
    gosub GETNAMES
     if invfile$="" and fhdr$<>"" then invfile$=fhdr$+".inv" : newopen=ktrue

elseif a$="d" then
  print "Default time increment DT for "
   for i=1 to nphtype
     print fnkphlab$(i);" ";
     dtdef(i)=fncrin(dtdef(i))
   next

elseif a$="s" then 
  print "Default seismogram length for " 
   for i=1 to nphtype
     print fnkphlab$(i);" ";
     tkept(i)=fncrin(tkept(i))
   next

elseif a$="l" then 
  print "Default lead time for " 
   for i=1 to nphtype
     print fnkphlab$(i);" ";
     tlead(i)=fncrin(tlead(i))
   next

elseif a$="f" then 
  print "Default Filter cutoff for " 
   for i=1 to nphtype
     print fnkphlab$(i);" ";
     fcuts(i)=fncrin(fcuts(i))
   next

elseif a$="o" then 
    input "Enter output file (.INV assumed, CR for default) "; invfile$
     if len(invfile$)=0 then invfile$=fhdr$
     invfile$=fnaddext$(invfile$,"INV")
     y$="n"  ' this and following 2 lines added 12 Jan,92; revised 4/25/1992
     if isfile(invfile$,f$,1)=0 then input " Append mode (y/n) ";  y$ : invfile$=f$
     newopen=not fnkyes(y$)

elseif a$="x" then call shellcom

elseif a$="b" then goto 260

end if : goto 250

260 '--- print header line to file
if newopen then
    open "o",2,invfile$
    print #2, using fnaform$(6); dt$;
    print #2, using fniform$(2); ohr; omin;
    print #2, using fniform$(3); osec*10;
    print #2, using fniform$(5); eplat*100;
    print #2, using fniform$(6); eplon*100; 
    print #2, using fniform$(3); dep;
    print #2, hend$
     newopen = kfalse
    close 2
end if

319 call setscreen(0)
320 call scclear
321 call scbox(1,17,ktrue,"SEISMOGRAM MENU","Enter M for SETUP MENU")
locate 3,3
 call scprint(3,3,"Input file header = "):call scprintv((0),(0),fhdr$)
call scprint((0),(0),"   Output file = "):call scprintv((0),(0),invfile$)

call scprint(5,3,"KEY   OPTION")
call scprompt(7,4,"L - list seismograms")
call scprompt(8,4,"M - return to PARAMETER SETUP MENU")
call scprompt(9,4,"N - proceed to next seismogram in list")
call scprompt(10,4,"R - reset COUNTER for seismograms, now at "):call scprintv((0),(0),str$(inord+1))
call scprompt(11,4,"F - reset FLAGS for seismograms")
call scprintc(12,5,"-- Enter KEY for option, OR --")
call scprintc(13,4,"enter station NAME, PHASE, and data TYPE (use format a4,i1,i1)")
call scprint(14,5,"PHASE 1:LP P    2:LP SH    3:SP P   4:BB P")
call scprint(15,5,"TYPE  1:WWSSN   2:GDSN")

 call scprintc(18,3,"-->") : line input a$

 if len(a$)=0 then 320
 if len(a$)=6 then 330

 st$=fnvalue$(1,1) : st$=lcase$(st$)

 if fniesc(st$) then 1140
 if st$="l" then gosub STALIST : goto 320
 if st$="m" then goto 250
 if st$="n" then goto 340
 if st$="r" then
   input " New seismogram start number "; ino
   if ino>0 then inord=ino-1
   inord=fnxlimit(inord,0,nsta-1)
 elseif st$="f" then
   for i=1 to nsta : flag$(i)="  ": next
 end if

 goto 320

330 '------ seismogram input by user
  ss$=lcase$(fnvalue$(1,4)): kph=fnvalue(5,1): ktype=fnvalue(6,1)
  if (not fninlim(kph,1,4)) or (not fninlim(ktype,1,2)) then 320 else 350

340 '------ the next seismogram in order
   locate 20,1

   inord=inord+1
      if inord>nsta then print "End of list" : call scpause(19,1) : goto 321
      ino=iord(inord)   
      if flag$(ino)="w " or flag$(ino)="s " then 340
      ss$=fnstanm$(ino)
      ktype=idatype(ino)

      if fnstaph$(ino)="lpz" then kph=1
      if fnstaph$(ino)="lpn" or fnstaph$(ino)="lpe" then kph=2
      if fnstaph$(ino)="spz" then kph=3
      if fnstaph$(ino)="bbz" then kph=4

350 '----- set up for new seismogram
     locate 20,1
     nout=0 : hpfilt=0
     m=1 : kcutleft=ktrue : noright=ktrue
     ksave=kfalse
     iphasetim=1 : iwindow(1)=1 : iwindow(2)=1

     if kph>2 then kph=kph+1
                     
     if kph=1 or kph=4 or kph=5 then ph$="P" else ph$="S"
     if kph=1 or kph=4 or kph=5 then nph=1 else nph=2

390  if nph=1 then

'--------- get P seismogram and interpolate -------------------- 

  for ista=1 to nsta
   if ss$=fnstanm$(ista) and fnstaph$(ista)=phin$(kph) and ktype=idatype(ista) then 460
  next
  print "Station ";ss$;" ";datin$(ktype);" ";phin$(kph);" not found" : call scpause(19,1) : goto 320

460 print "Processing station ";fnstanm$(ista);" ";datin$(ktype);" ";phin$(kph)
    kphtype=kphtyp(ista)

    if fninquit(kfalse) then inord=inord-1 : goto 320

    tmin=int(tfirstpt(ista))+1: tmax=int(tlastpt(ista))-1
    if idatype(ista)=1 then dt(ista)=dtdef(kphtype)
    dt=dt(ista)

    nintmax = (tmax-tmin)/dt + 1
    nintmax = fnamin(nintmax,nintmax0)
    nyoutmax=2^(int(fnalog2(nintmax))+1)
    nyoutmax = fnamin(nyoutmax,nintmax0)
    redim yout(nyoutmax)

    gosub GETDATA
    if ttime(ista)<>0 then tcalc=ttime(ista) else tcalc=0
    baz=baz(ista): az=az(ista): del=del(ista)

else

'-------- get 2 S components, interpolate, and rotate------------

570 nfnd=0

590 for iphz=1 to 2
600   istanum(iphz)=0
      tmin(iphz)=0 : tmax(iphz)=0

   for ista=1 to nsta
     if ss$=fnstanm$(ista) and fnstaph$(ista)=phin$(iphz+1) and ktype=idatype(ista) then 650
   next
   print "Station ";ss$;" ";datin$(ktype);" ";phin$(iphz+1);" not found": call scpause(19,1) : GOTO 730

650  tmin(iphz)=tfirstpt(ista) : tmax(iphz)=tlastpt(ista)
670  tmin=tmin(iphz) : tmax=tmax(iphz)
     kphtype=kphtyp(ista)
     baz=baz(ista): az=az(ista): del=del(ista)

 print "Processing station ";fnstanm$(ista);" ";datin$(ktype);" ";phin$(iphz+1)

     nfnd=nfnd+1 : istanum(iphz)=ista 
     if ttime(ista)<>0 then tcalc=ttime(ista) else tcalc=0
     smag(iphz)=xmag(ista)
     xmag0=xmag(ista)

730 next iphz

735 if nfnd=0 then d$=input$(1) : goto 320

    ista=fnamax(istanum(1),istanum(2))

740  if nfnd=2 then 
      tmin=fnamax(int(tmin(1)),int(tmin(2))) : tmin=tmin+1
      tmax=fnamin(int(tmax(1)),int(tmax(2))) : tmax=tmax-1
      xmag0=fnamax(smag(1),smag(2))             'use maximum xmag
      if tmax<=tmin then
        print "LPN, LPE do not overlap; skipped"
        for kj=1 to 2
          if istanum(kj)>0 then flag$(istanum(kj))="s "
        next
        goto 320
      end if 
    end if

   if fninquit(kfalse) then inord=inord-1 : goto 320

  for iphz=1 to 2
   ista=istanum(iphz)
    if ista>0 then
     if idatype(ista)=1 then dt(ista)=dtdef(kphtype)
     dt=dt(ista)
    end if
  next

    nintmax = (tmax-tmin)/dt + 1
    nintmax = fnamin(nintmax,nintmax0)
    nyoutmax=2^(int(fnalog2(nintmax))+1)
    nyoutmax = fnamin(nyoutmax,nintmax0)
    redim yout(nyoutmax), yout2(nintmax,2)

' get the data starting at TMIN
  for iphz=1 to 2
   ista=istanum(iphz)
    if ista>0 then
     if idatype(ista)=1 then dt(ista)=dtdef(kphtype)
     dt=dt(ista)
     gosub GETDATA
       for j=1 to nout
         yout2(j,iphz)=yout(j)*xmag0/smag(iphz)
       next
   end if
  next

960  call setscreen(kega)

  for iphz=1 to 2
   if istanum(iphz)>0 then
     ista=istanum(iphz)
     for i=1 to nout : yout(i)=yout2(i,iphz) : next
    if iphz=1 then
      title$=ss$+" LPN"
      iyzero=cint(0.22*kscry) : yamp=0.14*kscry : iyl=1
    else
      title$=ss$+" LPE"
      iyzero=cint(kscry/2) : yamp=0.14*kscry : iyl=0.35*nrowchar
    end if
    call PLOTSEIS (title$, iyzero, yamp, iyl, nout, yout())
   end if
  next

'--- rotate to get transverse SH
  bazr=fnrad(baz)
  for i=1 to nout
    yout(i)=sin(bazr)*yout2(i,1) - cos(bazr)*yout2(i,2)
  next

   if fninquit(kfalse) then inord=inord-1 : goto 319

   title$=ss$+" S" : iyzero=0.78*kscry : iyl=0.64*nrowchar
   call PLOTSEIS (title$, iyzero, yamp, iyl, nout, yout())

  d$=input$(1)

 if d$="n" or d$="e" then
   if d$="n" then ikp=1 else ikp=2
   for i=1 to nout : yout2(i,ikp)=yout2(i+1,ikp) : next
   goto 960
 end if
  
end if

'------------- edit seismogram on screen -----------------------------

1035 nkept=tkept(kphtype)/dt-1 : nlead=tlead(kphtype)/dt-1
     kplot=ktrue
     tpx=tcalc-tmin       ' theoretical time beyond start
     if abs(tpx/dt)>32000 then
       iphasetim = 1
     else
       iphasetim= cint(tpx/dt)+1
     end if
     if not fninlim(iphasetim,1,nout) then iphasetim=1
     iwindow(1)= fnamax(iphasetim-nlead,1)
     iwindow(2)= fnamin(iwindow(1)+nkept,nout)

1040 if not kplot then 1080
     title$=ss$+" "+ph$+" "+datin$(ktype)
1050 call setscreen(kega)
     kplot=kfalse : iyzero=kscry/2 : yamp=(nrowchar-4)/nrowchar*kscry/2 : iyl=1
     call PLOTSEIS (title$,iyzero,yamp,iyl, nout, yout())

1080 call scprintc(nrowchar-2,1,"OPTIONS: F1-Filter F2-New  F3-Write F4-Cut    F5-Phase  F6-DC shift")
     call scprintc(nrowchar-1,1," (arrows) 1-End     4-Left  6-Right  7-Begin 2,8-Step:") : print nstep;

1110 locate 2,1 : call inkey(d$)

if fnkeyf1(d$) or d$="f" then
    gosub FILTER

elseif fnkeyf2(d$) or d$="n" then
    if ksave then goto 319
    call cline(2,1)
    call scprint(2,1,"")
    input "Data not saved, skip (y/n)"; y$
    call cline(2,1)
      if fnkyes(y$) then 
       if nph=1 then 
         flag$(ista)="s "
       else
         for kj=1 to 2
          if istanum(kj)>0 then flag$(istanum(kj))="s "
         next
       end if
       goto 319
      end if

elseif fnkeyf3(d$) or d$="w" then
   gosub OUTPUTS
elseif fnkeyf4(d$) or d$="c" then
   gosub 1850
elseif fnkeyf5(d$) or d$="p" then
   gosub 1810
elseif fnkeyf6(d$) or d$="d" then
   gosub DCSHIFT

elseif d$="s" or fnmoveup(d$) then   'cursor speed
   nstep=nstep+1 
   call setstep (nstep)
   goto 1080

elseif fnmovedn(d$) then   'cursor speed
   nstep=nstep-1 
   call setstep (nstep)
   goto 1080

elseif fnmovelt(d$) or d$="l" then
   mleft=ktrue : gosub MCURSOR
elseif fnmovert(d$) or d$="r" then
   mleft=kfalse: gosub MCURSOR
elseif fnmovehm(d$) or d$="b" then
   kcutleft=ktrue : m=1
elseif fnmovend(d$) or d$="e" then
   kcutleft=kfalse: m=nout
end if

goto 1040  ' while moving cursor

1140 close
1155 goto 100

'======================================================================
OUTPUTS: '---- write to .INV file

'--- first check flag for write status
      if nph=1 and flag$(ista)<>"w " then 1180

      if nph=2 then
       for kj=1 to 2
          if istanum(kj)>0 and flag$(istanum(kj))<>"w " then 1180
       next
      end if

      call scprint (2,1,"")
        input "Data already saved, save again (y/n) "; y$
      call cline(2,1)
        if fnkyes(y$) then 1180 else return 319

1180  open "a",2,invfile$
      ksave=ktrue

     npout=iwindow(2)-iwindow(1)+1
     ymax=0
      for j=iwindow(1) to iwindow(2)
         ymax=fnamax(ymax,abs(yout(j)))
       next

'-- set status flags
     if kph=4 then nph=1
     if nph=1 then flag$(ista)="w "
  for kj=1 to 2
     if nph=2 and istanum(kj)>0 then flag$(istanum(kj))="w "
  next

'-- travel time to start of waveform
  tipxo=tmin + (iwindow(1)-1)*dt

'-- theoretical time of phase arrival after start of seismogram
  ttheor=fnamax(0.0,tcalc-tipxo)
  wt=1

'--- output station header line
     print #2, using fnaform$(4); fnstanm$(ista);
     print #2, using fniform$(1); nph;
     if instva(ista) then print #2,fnstrnum$(instva(ista)\4); : else print #2," ";
     print #2, using fniform$(1); inst(ista);
     print #2, using fnfform$(8,2); del;az;
     print #2, using fniform$(8); xmag(ista);
     print #2, using fnfform$(8,2); tipxo;
     print #2, using fnfform$(5,1); wt;
     print #2, using fniform$(5); iphasetim-iwindow(1)+1; npout;
     print #2, using fniform$(8); ymax;
     print #2, using fnfform$(6,2); dt;
     print #2, using fnfform$(6,2); ttheor;
     print #2, using fnfform$(6,2); hpfilt

 dform$=" "+fndform$(6,4)

 if idatype(ista)=2 then
    print #2, using fniform$(5); nzeros(ista);npoles(ista);
    print #2, using dform$; calconst(ista)

    for i=1 to nzeros(ista)
      print #2, using dform$; zeroi(i,1);zeroi(i,2);
      if (i mod 4=0) or i=nzeros(ista) then print #2,""
    next

    for i=1 to npoles(ista)
      print #2, using dform$; polei(i,1);polei(i,2);
      if (i mod 4=0) or i=npoles(ista) then print #2,""
    next
 end if

'---- output amplitude data
      kout=0
      for j=iwindow(1) to iwindow(2)
         kout=kout+1
         print #2, using dform$; yout(j);
         if (kout mod 8 = 0) or j=iwindow(2) then print #2, ""
       next

     close 2
1290 return 319

'======================================================================
MCURSOR: '-----subroutine to move cursor left & right
1670 call scline2(fnix(fntime(m)),fniy(yout(m))+littic,fnix(fntime(m)),fniy(yout(m))+1,0,kfalse)
1680 if mleft then m=fnamax(1,m-nstep) else m=fnamin(m+nstep,nout)
1690 call scline2(fnix(fntime(m)),fniy(yout(m))+littic,fnix(fntime(m)),fniy(yout(m))+1,kgrcol(kcolsyn),kfalse)
     call scprint(2,1,"Current point-  T= ")
     t1=fntime(m)+origin : call hrminsec(t1,ih,im,s)
     print fnstrnum$(ih)+":"+fnstrnum$(im)+":";
     print using fnfform$(5,2); s;
     print " A=";
     print using fniform$(7); yout(m)
 return
'======================================================================
1810 '----- IPHASETIM, the point which is the beginning of the phase
 if m<iwindow(1) then
   call cline(2,1): call scprint (2,1,"Warning-- this point is before start; reset start time")
   d$=input$(1) : call cline(2,1)
 elseif m>iwindow(2) then
   call cline(2,1): call scprint (2,1,"Warning-- this point is beyond end; reset end time")
   d$=input$(1) : call cline(2,1)
 else
   x1=fntime(iphasetim) : y1=yout(iphasetim)
   call scline2(fnix(x1),fniy(y1)-2*littic,fnix(x1),fniy(y1),0,kfalse)
   iphasetim=m               'starting point for phase
   x1=fntime(iphasetim) : y1=yout(iphasetim)
   call scline2(fnix(x1),fniy(y1)-2*littic,fnix(x1),fniy(y1),kgrcol(kcolsyn),kfalse)
 end if
 return

'======================================================================
1850 '-----get points IWINDOW 1 and 2, start & end points for time series

'--- erase tics
     for ij=1 to 2
       ix1=fnix(fntime(iwindow(ij))) : iy1=fniy(yout(iwindow(ij)))
       call scline2(ix1,iy1-littic,ix1,iy1,0,kfalse)
     next

'-- change values
    if kcutleft then
        if m>iphasetim then
           call cline(2,1) : call scprint (2,1,"Warning -- illegal left cut")
           d$=input$(1)
        else
           iwindow(1)=m
           if noright then iwindow(2)=fnamin(iwindow(1)+nkept,nout)
        end if
     else 
        if m<iphasetim or m<iwindow(1) then
           call cline(2,1) : call scprint (2,1,"Warning -- illegal right cut")
           d$=input$(1)
        else
           iwindow(2)=m
           noright=kfalse
        end if
     end if

'-- redraw tics
     for ij=1 to 2
      ix1=fnix(fntime(iwindow(ij))) : iy1=fniy(yout(iwindow(ij)))
       call scline2(ix1,iy1-littic,ix1,iy1,kgrcol(kcolsyn),kfalse)
     next
 return
'======================================================================
FILTER: '---------- filtering subroutine
2192 call cline(2,1)
     call scprint(2,1,"") : INPUT "Filter (1) or detrend (2)"; INC
     call cline(2,1)

     if not fninlim(inc,1,2) then return

     dt=dt(ista)

2225 if inc=2 then
      call m5detrnd (nout,dt,yout(),a,b) : kplot=ktrue 'detrend only

     elseif inc=1 then
      call scprint(2,1,"Filtering")
      nout = fnamin(nout,nintmax-1) 'revised from ..,nintmax) Apr 15,92
      n4=2^(int(fnalog2(nout-1))+1)
      redim yout2(n4,2)
      yout(nout+1)=yout(nout)/2.
      yout(1)=yout(1)/2.
      for i=nout+2 to n4 : yout(i)=0 : next    'pad with zeroes
      nf=n4/2+1
      dfreq=1./(n4*dt)
      hpfilt=fcuts(kphtype)

      call FORFFT(n4,dt,yout2(),yout())       'forward transform
      call BUT3P(nf, dfreq, yout2(), hpfilt)   'filter
      if fninquit(kfalse) then return
      call INVFFT(n4,dt,yout2(),yout())       'inverse transform

      for i=1 to nout
       yout(i)=yout(i)*dt
      next

     call scprint(2,1,"Filtered seismogram shown by dots; Press any key to continue")
2320 for i=1 to nout: call sccircle(fnix(fntime(i)),fniy(yout(i)),1,kgrcol(kcolsyn)):next
     kplot=ktrue
2330 d$=input$(1)
     end if
     return

'======================================================================
DCSHIFT: ' DC shift amplitude
     call cline(2,1)
3370 CALL SCPRINT(2,1,"") : INPUT "Amount to shift";DC
3380 for i=1 to nout
3390 yout(i)=yout(i)+dc
3400 next
     kplot=ktrue
3410 return 

'======================================================================
READSTA: '-------- read in file of stations and tables
     call scprint(18,3,"Reading station list ")
     if not ifopen("i",4,(syndir$+"m5statio.dat"),0,kfalse) then return
     for i=1 to 6 : wwscode$(i)="" : next
     for i=1 to maxwwsta
5731  line input #4, a$
      if fnvalue$(1,1)="'" then 5731
      if fnvalue$(1,3)="zzz" then 5780
      ijt=int(i/50.1)+1
      wwscode$(ijt)=wwscode$(ijt) + fnvalue$(1,4)
      adlat(i)=fnvalue(5,7) : adlon(i)=fnvalue(12,9)
     next
5780 NWST=I-1
     CLOSE 4
     call cline(18,3)
     return

'======================================================================
PLOTSEIS1: '----- routine to plot seismograms
    cls
    input "Enter length of seismogram in seconds to be plotted "; tlen
       if tlen=0 then return
'*
PSIN:print : print " Default name is ";fhdr$+".INV"
     input "Enter input file (.INV assumed, CR for default, q to quit) "; ifl1$
     if fniquit(ifl1$) then return
'*
     if len(ifl1$)=0 then ifl1$=fhdr$
     invfile$=fnaddext$(ifl1$,"INV")
     if isfile(invfile$, f$, -1)<>0 then goto PSIN
     invfile$=f$
'*
     fhdr$=fnsubext$(invfile$)
'*
3010 'if not ifopen("i",1,invfile$,0,ktrue) then return
     'if (instr(fhdr$,"*")>0 or len(fhdr$)=0) then fhdr$=left$(invfile$,instr(invfile$,".")-1)
'*
     open "i",1,invfile$
     line input #1, a$
     call setscreen(kega)
     nplots=7 : if kega=2 then nplots=4
     plotfile$=fhdr$
'*
3040 for iplot=1 to nplots
       if fninesc then close 1 : return
'*
3050 if eof(1) then 3340
3080 line input #1, a$
     gosub RDHEADER
     title$=s$+" "+fnkphlab$(kphtype)+left$(fncaltyp$(instva),1)
     ishort=fnishort
     if instin>0 then line input #1,d$ 'throw away the extra line
'*
     ymax=0
3150 for j=1 to nout
      input #1, yout(j)
      ymax=fnamax(ymax,abs(yout(j)))
     next
     if ymax>0 then gosub PLOT7
'*
3250 next
'*
3340 locate nrowchar-1,60: print "Press any key... ";: locate 1,1 : a$=input$(1)
     if fniesc(a$) then close 1: return
     if not eof(1) then 
          goto 3040 
       else 
          close 1 
          return 
       end if

'======================================================================
PLOTRAW: '----- routine to plot raw digitized data
     cls : locate 3,1
     input "Enter length of seismogram in seconds to be plotted "; tlen
       if tlen=0 then return

     gosub GETNAMES
     plotfile$=fhdr$

     call setscreen(kega)

     inord=0 : nplots=7 : if kega=2 then nplots=4

7040 for  iplot=1 to nplots

      if fninesc then return

      inord=inord+1
        if inord>nsta then 7340

      ista=iord(inord) 

      iph=1
      if (fnstaph$(ista)=phin$(2) OR fnstaph$(ista)=phin$(3)) then iph=2
      if fnstaph$(ista)=phin$(4) then iph=3

      title$=fnstanm$(ista)+" "+fnstaph$(ista)+" "+datin$(idatype(ista))+left$(fncaltyp$(instva(ista)),1)

      if fnstaph$(ista)="spz" then ishort=ktrue else ishort=false
      tmin=tfirstpt(ista) : tmax=tlastpt(ista)
      if not ishort then dt=1.0 else dt=0.3
      if idatype(ista)=2 then dt=dt(ista)

      gosub GETDATA

      xtsc=kscrx/180.*dt
      del=del(ista) : az=az(ista) : ttime=ttime(ista)

     ymax=0
     for j=1 to nout : ymax=fnamax(ymax,abs(yout(j))) : next

     tpx=ttime-tfirstpt(ista)      ' theoretical time beyond start
     iphasetim= cint(tpx/dt)+1
     nintmax = (tmax-tmin)/dt + 1  'added 8/Aug/91
     nintmax = fnamin(nintmax,nintmax0)  'added 8/Aug/91
     iphasetim=fnxlimit(iphasetim,1,nintmax-1)
     if ymax>0 then gosub PLOT7

7290  next iplot

7340 locate nrowchar-1,ncolchar-20: print "Press any key... ";: locate 1,1 : a$=input$(1)
     if fniesc(a$) then return
     if inord<nsta then goto 7040 else return

'======================================================================
PLOT7: ' plot seismogram on screen
'*
      if iplot=1 then
           call scclear
           call scprint(nrowchar-1,1,"Plot file: "+plotfile$+"    Tics: 5s")
           call scprintc((0),ncolchar-20,"Press Esc to quit")
      end if
'*
'*'-- set up location on screen; use top 21 lines of screen, bottom for text
     ssc=1.5/25. * kscry
     ix= ssc + (iplot-1)*3./25.*kscry
     nc=cint(1+(iplot-1)*3*(nrowchar/25))
'*
'*'-- offset to right 10 characters
     noff=10*kscrx/ncolchar  
     xtsc=(kscrx-noff)/tlen
'*
     call scprint(nc,1,title$)
     call scprint(nc+1,1,"") : print using fniform$(3)+"/"+fniform$(3); del; az;
'*
     it=fnit(iphasetim*dt)
     if fninlim(it,1,kscrx) then call scline2 (it,fnia(yout(iphasetim)+ymax/3),it,fnia(yout(iphasetim)-ymax/3),kgrcol(kcolsyn),kfalse)
'*
     n=0
     for j=1 to nout
         t=j*dt 
       if fninlim(fnit(t),1,kscrx) then 
         n=n+1
        if n=1 then  
          call scpoint(fnit(t),fnia(yout(j)),kgrcol(kcoldat))
        else 
          call scline1(fnit(t),fnia(yout(j)),kgrcol(kcoldat),kfalse)
        end if
       end if
     next
'*
'*'-- time tics every 5 secs
   for t=0 to tlen step 5 : call scpoint (fnit(t),fnia(0),kgrcol(kcolsyn)) : next
'*
return
'======================================================================
STALIST: '--- list stations by azimuth
  ksta=nsta : kcol=1
  icount=0
  b$="  "

     call scclear
     call scprint(nrowchar-1,1,"FLAG: w=written, s=skipped")
     call scprintc(nrowchar-1,40,"Press any key...")

9010 nl=fnamin(ksta,22)
     call scprint(1,kcol,"N F STA   PHS  TYPE DEL  AZ")
     call setvarcol

 for j=1 to nl
  icount=icount+1
  i=iord(icount)
  locate j+1,kcol
  print using fniform$(2);icount;
  print flag$(i);fnstanm$(i);b$;fnstaph$(i);b$;
  print using fnaform$(4); datin$(idatype(i));
  print using fniform$(4); del(i);az(i)
 next
 call settxtcol

ksta=ksta-nl
if ksta<1 then 9020

if kcol=1 then
  kcol = kcol+40 : goto 9010
else
  kcol=1 : d$=input$(1) : call scclear : goto 9010
end if

9020 d$=input$(1)
     return
'======================================================================
REORDER: ' re-order stations by azimuth and in order SP P, BB P, LP P, LP N, LP E
im=1
    for i=1 to nsta
    p=10
     for k=1 to 5
      if fnstaph$(i)=phin$(k) then p=k
     next
       if p=1 then p=1.00002  'lpz
       if p=3 then p=2.00001  'lpe (lpn=2.0000)
       if p=4 then p=1.00000  'spz
       if p=5 then p=1.00001  'bbz
       txdat(i)=p*1000+az(i)
    next

    for i=1 to nsta
     azim=12000.
      for j=1 to nsta
       if txdat(j)<azim then azim=txdat(j) : im=j
      next
     txdat(im)=9999.0
      iord(i)=im
    next
return
'======================================================================
RDHEADER: ' read station header line
       s$=fnvalue$(1,4)
       iph=fnvalue(5,1)
       instin=fnvalue(6,1) AND 1
       instva= (fnvalue(6,1) AND 6)*4
       inst=fnvalue(7,1)
       del=fnvalue(8,8)
       az=fnvalue(16,8)
       xmag=fnvalue(24,8)
       xt=fnvalue(32,8)
       px=fnvalue(40,5)
       iphasetim=fnvalue(45,5)
       nout=fnvalue(50,5)
       ymax=fnvalue(55,8)
       dt=fnvalue(63,5)

'--- assign phase type
     if iph=1 and inst<3 then kphtype=1
     if iph=2 and inst<3 then kphtype=2
     if inst=3 then kphtype=3
     if iph=1 and inst=4 then kphtype=4
     if iph=2 and inst=4 then kphtype=5
     if inst=5 then kphtype=6
     if inst=7 then kphtype=7

  if dt=0 and fnilong then  dt=.5
  if dt=0 and fnishort then dt=.1

'--- read zeroes and poles for GDSN seismograms 
if fningdsn then
   input #1, nzero, npole, calcons
   for j=1 to 2*nzero : input #1, zeroi : next j
   for j=1 to 2*npole : input #1, polei : next j
end if
return
'======================================================================
sub slatlon (ss$,stalat,stalon,ier) static
'------ get lat and lon for station ss$
  shared nwst, wwscode$(), adlat(), adlon()
   ier=0
 for i=1 to nwst
  if fnwwsnm$(i)=ss$ then stalat=adlat(i) : stalon=adlon(i) : exit sub
 next
   ier=-1
end sub
'======================================================================
 sub hrminsec ( sectime, ihr, imin, secs ) static
  secs=sectime
  ihr = int(secs/3600.)
   secs=secs - ihr*3600.0
  imin=int(secs/60.)
   secs=secs - imin*60.
 end sub
'==================================================================
sub PLOTSEIS (title$, iyzero, yamp, iyl, nout, y(1)) static
'--- plot interpolated data

shared dt, dt(), kscrx, kscry, del, az, baz, tmin, tmax, origin
shared nrowchar, ncolchar, iphasetim, littic, iwindow()
shared nkept, tcalc, ymax, m, ista ,kgrcol()

   t=tmin : ixt=fnix(t) : sc=5.0
   while ixt<kscrx   ' tic marks every 5 seconds
     call scpoint(ixt,iyzero,kgrcol(kcolsyn))
     t=t+sc
     ixt=fnix(t)
   wend

     ymax=0
   for i=1 to nout
     ymax=fnamax(abs(y(i)),ymax)
   next

  call scprint(iyl,1,title$)
  if iyl=1 then 
      print " Del:";fnfix(del,1);"Az:"; fnfix(az,1);
      print "Baz:";fnfix(baz,1);"Amax:"; int(ymax);

      t1=tmin+origin : call hrminsec(t1,ih,im,s)
      print " Start time: "; fnstrnum$(ih)+":"+fnstrnum$(im)+":";fnfix(s,1);
   end if

     call scpoint(fnix(fntime(1)),fniy(y(1)),kgrcol(kcoldat))
     for i=2 to nout : call scline1(fnix(fntime(i)),fniy(y(i)),kgrcol(kcoldat),kfalse): next

'Phase arrival
     x1=fntime(iphasetim) : y1=y(iphasetim)
     call scline2(fnix(x1),fniy(y1)-2*littic,fnix(x1),fniy(y1),kgrcol(kcolsyn),kfalse)
'Current point
     x1=fntime(m) : y1=y(m)
     call scline2(fnix(x1),fniy(y1)+littic,fnix(x1),fniy(y1),kgrcol(kcolsyn),kfalse)
'Left cut and Right cut
   for ji=1 to 2
     ix1=fnix(fntime(iwindow(ji))) : iy1=fniy(y(iwindow(ji)))
     call scline2(ix1,iy1-littic,ix1,iy1,kgrcol(kcolsyn),kfalse)
   next
'theoretical time
     kx=fnxlimit(fnix(tcalc),1,kscrx)
     'call scline2(kx,iyzero-littic,kx,iyzero,1,kfalse)
     call sccircle (kx,iyzero,3,kgrcol(kcolsyn))
end sub

'======================================================================
READINPUT:

open "i",1,infile$

'---input event header
line input #1, a$

'-- if the first character is a number, it is assumed to be the header
  iyr1$=fnvalue$(1,1)
  if fninlim(asc(iyr1$),48,57) then line input #1, a$ 
 
  goto RD3

RD1: if eof(1) then goto INPEND
     if fninquit(kfalse) then close 1 : return
 
     line input #1, a$

RD3:
'---- station name must start with a letter
     k=asc(fnvalue$(1,1)) : if not fninlim(k,65,122) then goto RD1

    nsta=nsta+1 : if nsta>nseismax then nsta=nsta-1 : goto RD2

     ss$=lcase$(fnvalue$(1,4))
     inst(nsta)=fnvalue(5,1)
     ph$=lcase$(fnvalue$(6,3))
     instva(nsta)=instr("va",mid$(ph$,2,1))*8
     if instva(nsta) then
       if left$(ph$,1)="b" then mid$(ph$,2,1)="b" : else mid$(ph$,2,1)="p"
     end if

'-- flag WWSSN or GDSN data
      if inst(nsta)<=3 then idatype(nsta)=1 else idatype(nsta)=2
      if idatype(nsta)=2 then igds=ktrue else igds=kfalse

   print using fniform$(2);nsta;
   print " Reading station ";ss$;"  Phase ";ph$;" Type ";datin$(idatype(nsta));left$(fncaltyp$(instva(nsta)),1)

'--- assign phase-data type
   if igds then
     if ph$="lpz" then kphtype=4
     if ph$="lpn" or ph$="lpe" then kphtype=5
     if ph$="spz" then kphtype=6
     if ph$="bbz" then kphtype=7
   else
     if ph$="lpz" then kphtype=1
     if ph$="lpn" or ph$="lpe" then kphtype=2
     if ph$="spz" then kphtype=3
   end if

   kphtyp(nsta)=kphtype

   if igds then
     avscale=fnvalue(9,10)    'New Form
     dt(nsta)=fnvalue(40,5)
     xmag(nsta)=1.0
     xm=xmag(nsta)/avscale
   else
     xmag(nsta)=fnvalue(9,10)
   end if

   if fnpwave then iph=1 else iph=2


     flag$(nsta)="  "

'--- TFIRSTPT= time from origin time to reference time
     s = 60.* fnvalue(33,2) +3600.*fnvalue(31,2)
     tfirstpt= fnvalue(35,5) + s - origin

'--- get station coordinates
    ier=0
    slat=fnvalue(45,7) : slon=fnvalue(52,7)
    if slat=0 and slon=0 then call slatlon(ss$,slat,slon,ier)
      if ier then print " *** station coordinates not found " : goto 1970

     ijt=int(nsta/50.1)+1
     stacode$(ijt)=stacode$(ijt)+ss$
     staphase$(ijt)=staphase$(ijt)+ph$

'--- get delta, azimuth and back-azimuth
      call m5delaz(slat,slon,eplat,eplon,delta,daz)
        del(nsta)=delta : az(nsta)=daz
      call m5delaz(eplat,eplon,slat,slon,delta,daz)
        baz(nsta)=daz

'--- get theoretical travel time
      call travtime (iph, delta, dep, ttime, dtdh, dtdz)
      call ellipcor (iph, delta, eplat, slat, ecor) 
        ttime(nsta)=ttime+ecor
   
1970 if igds then

  input #1, nn, nzeros(nsta), npoles(nsta), calconst(nsta)
    if calconst(nsta)<0 then xm=-xm : calconst(nsta)=-calconst(nsta)

    for i=1 to nzeros(nsta)
        input #1, zeroi(i,1), zeroi(i,2)
     next

     for i=1 to npoles(nsta)
       input #1, polei(i,1), polei(i,2)
     next

'      npts(nsta)=nn
     for jj=1 to nn
      input #1, y : if jj<=nintmax0 then yout(jj)=y*xm
     next

     if nn>nintmax0 then print "More than ";nintmax0;" points; truncating" : nn=nintmax0
      npts(nsta)=nn

'--- time from origin time to first and last points in seismogram
     tfirstpt(nsta)=tfirstpt
     tlastpt(nsta)=tfirstpt+(nn-1)*dt(nsta)
else
    xm=1.0
     if xmag(nsta)<0 then xm=-1.0 : xmag(nsta)=-xmag(nsta)

     input #1, nn : 'npts(nsta)=nn

     for jj=1 to nn
      input #1, x,y
      if jj<=inmax then txdat(jj)=tfirstpt+x/100: ydat(jj)=y*xm
     next

     if nn>inmax then print "More than ";inmax;" points; truncating" : nn=inmax
      npts(nsta)=nn

' time from origin time to first and last points in seismogram
     tfirstpt(nsta)=txdat(1) 
     tlastpt(nsta)=txdat(nn) 

end if

    if ier then 
       nsta=nsta-1  'previous station had no coordinates so overwrite
    else
       gosub SETDATA
    end if

   goto RD1

RD2: print " Too many seismograms in file: MAX = ";nseismax
     print "   Only that many will be read in"

INPEND: close 1
return

'======================================================================
SETDATA:  ' print data to random access file
  open "r",3,drive$+"data.tmp",4 : field #3, 4 as tran$

   indexdat&(nsta)=indexs&
   kst&=indexs&

 if idatype(nsta)=1 then 
'------- WWSSN data out
    lset tran$ = mks$(txdat(1)) : put #3, kst&
    lset tran$ = mks$(ydat(1))  : put #3
   for j=2 to npts(nsta)
    if txdat(j)=txdat(j-1) then txdat(j)=txdat(j)+0.01
    lset tran$ = mks$(txdat(j)) : put #3
    lset tran$ = mks$(ydat(j))  : put #3
   next
  npoints=2*npts(nsta)

 else

'------ GDSN data output
      lset tran$ = mks$(zeroi(1,1)) : put #3, kst&
      lset tran$ = mks$(zeroi(1,2)) : put #3
    for j=2 to nzeros(nsta)
      lset tran$ = mks$(zeroi(j,1)) : put #3
      lset tran$ = mks$(zeroi(j,2)) : put #3
    next
    for j=1 to npoles(nsta)
      lset tran$ = mks$(polei(j,1)) : put #3
      lset tran$ = mks$(polei(j,2)) : put #3
    next
    for j=1 to npts(nsta)
      lset tran$ = mks$(yout(j)) : put #3
    next
    npoints=npts(nsta)+ 2* (nzeros(nsta)+npoles(nsta))
 end if

  indexs&=indexs&+npoints
  close 3
return
'======================================================================
GETDATA:
'--- read data from random access file for station number ISTA
'-- get seismogram from TMIN to TMAX, put NOUT points in YOUT()
'-- interpolate WWSSN data

 xmag=xmag(ista)

  open "r",3,drive$+"data.tmp",4 : field #3, 4 as tran$

 if idatype(ista)=1 then
    kst&=indexdat&(ista)
    npts=npts(ista)
     get #3, kst& : txdat(1)=cvs(tran$)
     get #3      : ydat(1)=cvs(tran$)
   for j=2 to npts
     get #3 : txdat(j)=cvs(tran$)
     get #3 : ydat(j)=cvs(tran$)
   next

   call m5wigint (tmin,tmax,txdat(),ydat(),npts,dt,yout(),nout)

 else
'---   get GDSN data starting at TMIN until TMAX

   nout=fnamin(npts(ista),(tmax-tmin)/dt+1)
   kst&=indexdat&(ista)

'-- first get poles and zeros
      get #3, kst& : zeroi(1,1)=cvs(tran$)
      get #3      : zeroi(1,2)=cvs(tran$)
    for j=2 to nzeros(ista)
      get #3      : zeroi(j,1)=cvs(tran$)
      get #3      : zeroi(j,2)=cvs(tran$)
    next
    for j=1 to npoles(ista)
      get #3      : polei(j,1)=cvs(tran$)
      get #3      : polei(j,2)=cvs(tran$)
    next

    nleader=(tmin-tfirstpt(ista))/dt
    for j=1 to nleader : get #3 : next

   for j=1 to nout
     get #3 : yout(j)=cvs(tran$)
   next

 end if

   call m5detrnd (nout, dt, yout(), a,b)

  close 3
return
'======================================================================
GETNAMES:
    cls
    nsta=0 : indexs&=1 : khypo=kfalse
    for i=1 to 5 :stacode$(i)="" :staphase$(i)="" : next

GETWWSN:
     print
     input "Do you have a WWSSN (.DIG) input file (y/n) "; y$
     inwwssn=fnkyes(y$)

      if inwwssn then
         print : print " Default WWSSN input file is "; ucase$(fhdr$)+".DIG"
         input " Enter file name prefix (.DIG assumed) (CR for default) "; ifl1$

         if len(ifl1$)=0 then ifl1$=fhdr$
         digfile$=fnaddext$(ifl1$,"DIG")

         if isfile(digfile$, f$, -1)<>0 then goto GETWWSN
         fhdr$=fnsubext$(f$)  'revised 4/25/1992
         digfile$=f$
       end if

GETGDSN:
     print : input "Do you have a GDSN (.DSN) input file (y/n) "; y$
     ingdsn=fnkyes(y$)

     if ingdsn then
       print : print " Default GDSN input file is "; ucase$(fhdr$)+".DSN"
       input " Enter file name prefix (.DSN assumed) (CR for default) "; ifl1$
       if len(ifl1$)=0 then ifl1$=fhdr$
       dsnfile$=fnaddext$(ifl1$,"DSN")

       if isfile(dsnfile$, f$, -1)<>0 then goto GETGDSN
       if not inwwssn then fhdr$=fnsubext$(f$)  'revised 4/25/1992
       dsnfile$=f$
     end if

'-- get hypocenter
     if inwwssn then infile$=digfile$ : gosub READHYPO
     if ingdsn  then infile$=dsnfile$ : gosub READHYPO
     if not khypo then gosub GETHYPO

'-- origin time in seconds
         origin=ohr*3600 + omin*60 + osec

'-- read in data, put in RAF DATA.TMP
         if inwwssn then infile$=digfile$ : gosub READINPUT
         if ingdsn  then infile$=dsnfile$ : gosub READINPUT

'-- re-order by azimuth
         if inwwssn or ingdsn then gosub REORDER

return

'======================================================================
READHYPO:
'-- check files for hypocenter information

open "i",1,infile$

'---input event header
line input #1, a$

'-- if the first character is a number, it is assumed to be the year
  iyr1$=fnvalue$(1,1)

  if fninlim(asc(iyr1$),48,57) then 

  if (not khypo) then
    dt$=fnvalue$(1,6)
    ohr=fnvalue(7,2)
    omin=fnvalue(9,2)
    osec=fnvalue(11,3)/10
    eplat=fnvalue(14,5)/100
    eplon=fnvalue(19,6)/100 
    dep=fnvalue(25,3)
    hend$=fnvalue$(28,len(a$)-27)
    khypo=ktrue
  end if

  end if
close 1
return
'======================================================================
GETHYPO:
'-- read in hypocenter data

  if khypo then return
    cls
    print " Enter hypocentral data " : print
    input " Date of event (YYMMDD) "; dt$
    input " Hour "; ohr
    input " Minute "; omin
    input " Second "; osec
    input " Latitude (+North) "; eplat
    input " Longitude (+East) "; eplon
    input " Depth (km) "; dep

return
'==================================================================  
sub setstep (imove) static
'-- set step for cursor motion
  if imove= 0 or imove=31 then imove=60
  if imove=21 or imove=59 then imove=30
  if imove=11 or imove=29 then imove=20
  if imove=19 then imove=10
  if imove=61 then imove=1
end sub
'======================================================================
rem $INCLUDE: 'm5wigint'
rem $INCLUDE: 'm5ttimes'
rem $INCLUDE: 'm5rdset'
'======================================================================

