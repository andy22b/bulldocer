' MT5GDSN.BAS                        QuickBASIC 4.0
' 1/10/87, GA Abers.
' Program to read in USGS-format binary GDSN seismic data in RETRV files
'    based on Compact-Disk version
'    This version reformats to ASCII.
'    Data is stored in 500-point "buffers":  one record for the header,
'      one record for each buffer, and a "last" record
'    The last record has one data point of value 1E+10
'    Reads a FORTRAN binary file, so:
'       (1)  Characters, Integers, and Real Numbers are 4 bytes
'       (2)  Real Numbers are in IEEE format
'       (3)  Each "record" is composed of as many 128-byte "lines" as are
'           necessary for the data
'       (4)  Each line (128 or less bytes) is preceeded by a 2 byte marker
'       "CH2" is used to strip the marker.  Integers are read in LONG.
'
'    modified for SYN3 .DSN output format on 4/15/89 by McCaffrey
'    modified for S4FNOPEN routine on 10/17/90 by Zwick

DEFINT I-N: OPTION BASE 1
CONST NSTAMAX=1 'dummy variable

rem $INCLUDE: 'M5COMSCR'

DECLARE FUNCTION ch4$ ()
DECLARE SUB jdtommdd (iyr&, idoy&, imo, idy)
DECLARE SUB rdseis (nflg%, ndat&, ss() AS SINGLE)

rem $INCLUDE: 'm5fndecl'
rem $INCLUDE: 'm5functn'

'--- welcome
prgnm$="MT5GDSN"
authors$="R. McCaffrey, G. Abers, & P. Zwick"
call welcome

DIM code(4) AS STRING * 4   '-- station
DIM chn(4) AS STRING * 4    '-- instrument
DIM iyr AS LONG, idoy AS LONG, ihr AS LONG, imin AS LONG
DIM np AS LONG, nz AS LONG, nflag(3) AS LONG, ndat AS LONG
DIM poles(2, 30), zeroes(2, 20), fresp(2)

ndatmax = 15000
nesmax=3072
DIM seis(ndatmax), ss(nesmax), iml(12)

def fnix(x)=(kscrx-40)*(1/20 + (x-n1plot+1)/nplotted)
def fniy(a)=kscry/2  - a*0.4*kscry/ampmax

DATA 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
FOR i = 1 TO 12: READ iml(i): NEXT i

nfiles=7
dim  kflag$(nfiles)
DIM  files$(nfiles), infile(nfiles)

DATA "LPZ","LPN","LPE","SPZ","SPN","SPE","IPZ","IPN","IPE"
FOR i = 1 TO nfiles :  READ files$(i) : kflag$(i)=" " :infile(i)=-1 : NEXT i
infile(8)=0 : infile(9)=0

'-- function keys
def fnkeyf1(x$) = (x$=chr$(59))
def fnkeyf2(x$) = (x$=chr$(60))
def fnkeyf3(x$) = (x$=chr$(61))
def fnkeyf4(x$) = (x$=chr$(62))
def fnkeyf5(x$) = (x$=chr$(63))
def fnkeyf6(x$) = (x$=chr$(64))
def fnkeyf7(x$) = (x$=chr$(65))

retrv$="RETRV"
ascfil$="RETRV.DSN"
nstep=1
khypo=0
gosub CHKFILES

call scpause(18,3)

'--- menu
MENU:
call setscreen(0)
call scbox(1,18,ktrue,"SETUP MENU","Press Q to quit")
call scprint(2,2,"KEY  OPTION                                VALUE")
call scprompt(4,3,"N - Name for input RETRV files"):call scprintv(4,40,retrv$)
call scprompt(5,3,"O - Output ASCII file name"):call scprintv(5,40,ascfil$)
call scprompt(7,3,"C - Components to read ( *=component file not found)")
nl=0

FOR i = 1 TO nfiles
    if (i mod 3) =1 then nl=nl+1 : locate 7+nl,8
    print files$(i);" ";
    call scprintv((0),(0),kflag$(i)+" "+fnyesno$(infile(i))):print "    ";
NEXT i

call scprompt(12,3,"H - Hypocenter   "):call setvarcol
   print using fnfform$(7,3); eplat;
    if eplat>=0 then print "N"; else print "S";
   print using fnfform$(8,3); eplon;
    if eplon>=0 then print "E"; else print "W";
   call settxtcol
   print "      Depth:";: call scprintvfu((0),(0)," ###.#",dep):print " km"
   call scprint(13,10,"Origin time: "): call setvarcol: print dt$;
   print using "  ##:##:##.##";ohr;omin;osec

   call scprompt(15,3,"X/"):call scprompt(15,5,"Q - Shell/Quit")
   call scprompt(16,3,"B - Begin")

mnprompt:
   a$=""
   ik = input1(19,1,"Enter KEY for option ",a$,kfalse,kfalse,ktrue,"nochxbq")

if a$="q" then call endcheck:goto mnprompt

if a$="b" then goto BEGIN

if a$="n" then
   p$=retrv$
   irow=0
   call scprintc(irow,1,"Enter new prefix for RETRV files including path names (blank for RETRV)")
   if inputn(irow+1,2,"--",p$,ktrue,ktrue,ktrue,"a",40)<>27 then
      if len(p$)=0 then retrv$="RETRV" else retrv$=p$
      retrv$=fnsubext$(retrv$)
      if ascfil$="RETRV.DSN" then ascfil$=retrv$+".DSN"
      gosub CHKFILES
   end if

elseif a$="h" then
   gosub GETHYPO

elseif a$="o" then
   if len(ascfil$)=0 then p$=fnsubext$(retrv$)+".DSN" else p$=ascfil$
   irow=0
   call scprintc(irow,1,"Enter output ASCII file name (.DSN) (CR for default)")
   if inputn(irow+1,2,"--",p$,ktrue,ktrue,ktrue,"a",40)<>27 then
      if len(p$)=0 then
         ascfil$=fnsubext$(retrv$)
      else
         ascfil$=p$
      endif
      ascfil$=fnaddext$(ascfil$,"DSN")
   end if

elseif a$="c" then
   irow=0
   call scprintc(irow,(0),"Select components to read (y/n for each):")
   FOR i = 1 TO nfiles
     if (i mod 3) =1 then irow=irow+1 : icol=8
     if infile(i) then y$="Y" else y$="N"
     if input1(irow,icol,files$(i)+" "+kflag$(i),y$,kfalse,ktrue,kfalse,"yn")=27 then exit for
     infile(i)=fnkyes(y$)
     icol=icol+13
   NEXT i
   gosub CHKFILES

elseif a$="x" then 
   call shellcom

end if

goto MENU

BEGIN:

OPEN "o",9, ascfil$
'--- print header line to file
if khypo then
    print #9, using fnaform$(6); dt$;
    print #9, using fniform$(2); ohr; omin;
    print #9, using fniform$(3); osec*10;
    print #9, using fniform$(5); eplat*100;
    print #9, using fniform$(6); eplon*100;
    print #9, using fniform$(3); dep
end if
close 9

for icompo=1 to nfiles
   IF infile(icompo) THEN

  finam$ = retrv$ + "." + files$(icompo)
  kflag$(icompo)="*"


if isfile(finam$, f$, 1)=0 then

   kflag$(icompo)=" "

OPEN finam$ FOR BINARY AS #1

  if fninlim(icompo,4,6) then inst=5 else inst=4

40 WHILE NOT EOF(1)

   '---- read in header
    GOSUB READHDR
    secstart = ihr*3600. + imin*60. + sec

    npts = 0
    cls
    print " Reading data for ";stn$;" ";phs$

'--- read data
    amax = -1e20 : amin = 1e20 : sumsi=0

rddata:
    CALL rdseis(nflg%, ndat, ss())
    IF nflg% THEN GOTO PLTDATA0
    IF (ndat + npts) > ndatmax THEN GOTO PLTDATA0
    FOR i = 1 TO ndat
      si=ss(i)
      seis(npts + i) = si 
      sumsi=sumsi+si
      amax=fnamax(si,amax)
      amin=fnamin(si,amin)
    NEXT i
    npts = npts + ndat
    GOTO rddata

PLTDATA0:    '--- plot
    
'-- set initial window to all points
    nstart=1 : nlast=npts

    avg= sumsi/(nlast-nstart+1)
    ampmax = fnamax(abs(amin-avg), abs(amax-avg))
    IF not nflg% then
       PRINT "** Number of amplitude points exceeds array size (";ndatmax;")"
       PRINT "** Decrease the number of points and start over " 
       call scpause(10,3): go to 220
    END IF

    IF npts < 2 THEN 
      PRINT "*** Nothing to Plot ***"
      call scpause(10,3): go to 210
    end if

    GOSUB STATINFO
    call cline(nrowchar-2,1)

    allmax=ampmax : allavg=avg
    call scpause (nrowchar-2,1) : goto 55

PLTDATA:    '--- plot
    gosub GETMAX

55   call setscreen(kega)
     if ampmax<1 then ampmax=1

'-- controls plotting scales
    n1plot=nstart
    nplotted= nlast-nstart+1

'--- plot data
      ix=fnix(nstart) : iy=fniy(0) : call scpoint (ix,iy,kgrcol(kcoldat))

    FOR i = nstart TO nlast
      a=seis(i)-avg 
      x1=fnix(i) : y1=fniy(a)
      call scline1 (cint(x1), cint(y1), kgrcol(kcoldat), kfalse)
    NEXT i

'-- plot time axis
tincr=60
if nplotted*dt < 150 then tincr=10
if nplotted*dt < 60 then tincr=5
if nplotted*dt < 20 then tincr=1
nincr=tincr/dt
    FOR i = 0 TO nlast step nincr
      a=.05*amax
      x1=fnix(i) : y1=kscry*(1-5/nrowchar) :y2=kscry*(1-4/nrowchar)
      call scline2 (cint(x1),cint(y1),cint(x1), cint(y2), kgrcol(kcoldat), kfalse)
    NEXT i

    gosub HEADER
    gosub PLOTTICS

PRINTOPT:
    locate nrowchar-3,1 : print "OPTIONS: F1-Station info.   F3-Write to file  F5-Cursor window limits F7-Replot"
    locate nrowchar-2,1 : print "         F2-New seismogram  F4-All data       F6-Enter window limits   Q-exit"

 call scprintc(nrowchar-1,1,"Enter KEY for option -> ")
 call inkey (d$)
 call cline(nrowchar-1,1) : locate nrowchar-1,1

 if fniquit(d$) then 220

 if fnkeyf3(d$) or d$="w" then
   GOSUB OUTPUT9 : goto 210

 elseif fnkeyf4(d$) or d$="a" then
   nstart=1 : nlast=npts
   ampmax=allmax : avg=allavg
   goto 55

 elseif fnkeyf6(d$) or d$="e" then
   call cline(nrowchar-3,1) : call cline(nrowchar-2,1) : call cline(nrowchar-1,1) :
   s$=fnstrnum$(nstart):n=0
   if inputn(nrowchar-3,1,"Enter new STARTING point for window (1-"+_
       fnstrnum$(npts)+")",s$,ktrue,ktrue,kfalse,"i",5)<>27 then n=val(s$)
    if fninlim(n,1,npts) then nstart=n
   s$=fnstrnum$(nlast):n=0
   if inputn(nrowchar-2,1,"Enter new ENDING point for window ("+fnstrnum$(nstart)+"-"+_
       fnstrnum$(npts)+")",s$,ktrue,ktrue,kfalse,"i",5)<>27 then n=val(s$)
    if fninlim(n,nstart,npts) then nlast=n
   goto PLTDATA

 elseif fnkeyf7(d$) or d$="r" then
    goto 55

 elseif fnkeyf2(d$) or d$="n" then
    call cline(nrowchar-3,1) : call cline(nrowchar-2,1) : locate nrowchar-3,1
    print "Are you sure (y/n) ";
    call inlet (kfalse, y$, klen, kval, ichoice)
    if fnkyes(y$) then 210

 elseif fnkeyf1(d$) or d$="s" then
    gosub STATINFO :  call scpause (nrowchar-2,1) : goto 55

 elseif fnkeyf5(d$) or d$="c" then
    gosub WINDOW1 : goto PLTDATA
 end if
   goto PRINTOPT

210 WEND
999 CLOSE #1

end if
end if
next

220   close : goto MENU

'-------------------------------------------------------------------------
HEADER:  ' print screen header
     locate 1,1 : PRINT stn$; " "; phs$; " Num. points:";npts;
     PRINT " Amin:"; int(amin); " Amax:";int(amax); " Start:";nstart;" End:";nlast;
        secs0=secstart+(nstart-1)*dt
        ihr1 = int(secs0/3600.) : secs0=secs0 - ihr1*3600.0
        imin1=int(secs0/60.)    : secs0=secs0 - imin1*60.
    locate 2,1 : PRINT " Window start time: ";ihr1;"h ";imin1;"m ";fnfix(secs0,2);"s"
        secs0=secstart+(nlast-1)*dt
        ihr1 = int(secs0/3600.) : secs0=secs0 - ihr1*3600.0
        imin1=int(secs0/60.)    : secs0=secs0 - imin1*60.
    locate 3,1 : PRINT "          end time: ";ihr1;"h ";imin1;"m ";fnfix(secs0,2);"s";
     print "   Tics:";tincr; " secs";
 return

'-------------------------------------------------------------------------
WINDOW1: '-- allow windowing of seismogram

      kbegin=ktrue : kend=kfalse

WINDOW2: gosub HEADER
WINDOW3: gosub PLOTTICS

   call cline(nrowchar-3,1) : call cline(nrowchar-2,1)

locate nrowchar-2,1 
print "OPTIONS: (arrows) 4-Left 6-Right 7-Begin 1-End  2,8-Step";
print using fniform$(3); nstep; : print "  Esc to exit";

   call inkey(d$)
   if fniesc(d$) then return

'-- erase tics
  call scline2 (cint(x1), cint(y1), cint(x2), cint(y2), 0, kfalse)
  call scline2 (cint(x3), cint(y3), cint(x4), cint(y4), 0, kfalse)

  call cline(nrowchar-3,1)

   if fnmovelt(d$) or d$="l" then
    if kbegin then nstart=nstart-nstep : nstart=fnamax(1,nstart)
    if kend   then nlast=nlast-nstep : nlast=fnamax(nstart,nlast)
   elseif fnmovert(d$) or d$="r" then
    if kbegin then nstart=nstart+nstep : nstart=fnamin(nstart,nlast)
    if kend   then nlast=nlast+nstep : nlast=fnamin(nlast,npts)
   end if

'-- cursor speed
   if fnmoveup(d$) or d$="s" then 
     nstep=nstep+1 
     call setstep(nstep)
     goto WINDOW3
   elseif fnmovedn(d$) then
     nstep=nstep-1 
     call setstep(nstep)
     goto WINDOW3
   end if

   if fnmovehm(d$) or d$="b" then kbegin=ktrue : kend=kfalse
   if fnmovend(d$) or d$="e" then kbegin=kfalse: kend=ktrue

   goto WINDOW2
return

'-------------------------------------------------------------------------
PLOTTICS:
'-- plot tics
      snstart=seis(nstart)-avg: snlast=seis(nlast)-avg
      x1=fnix(nstart) : y1=fniy(snstart+0.2*ampmax)
      x2=fnix(nstart) : y2=fniy(snstart+0.02*ampmax)
      x3=fnix(nlast) : y3=fniy(snlast-0.2*ampmax)
      x4=fnix(nlast) : y4=fniy(snlast-0.02*ampmax)

      call scline2 (cint(x1), cint(y1), cint(x2), cint(y2), kgrcol(kcolsyn), kfalse)
      call scline2 (cint(x3), cint(y3), cint(x4), cint(y4), kgrcol(kcolsyn), kfalse)
return

'-------------------------------------------------------------------------
STATINFO:
call setscreen(0)
call scbox (1,18,ktrue,"STATION DATA",stn$+" "+phs$)
  locate 3,3
    PRINT "Station: ";stn$; " Component: "; phs$;
  locate 5,3: print "Location: "; alat; 
    if alat>=0 then print "N"; else print "S";
  locate 6,3: print "          "; alon; 
    if alon>=0 then print "E"; else print "W";
  locate 7,3: print "Elevation:"; elev; " meters"
  locate 9,3
    PRINT "Year"; iyr; " Julian Day:"; idoy;

  locate 9,30
    PRINT "Record start time: HH:MM:SS.SS ";
        secs0=secstart
        ihr1 = int(secs0/3600.) : secs0=secs0 - ihr1*3600.0
        imin1=int(secs0/60.)    : secs0=secs0 - imin1*60.
    print using fniform$(2)+":"; ihr1;imin1; 
    print using fnfform$(5,2); secs0
        secs0=secstart+(npts-1)*dt
        ihr1 = int(secs0/3600.) : secs0=secs0 - ihr1*3600.0
        imin1=int(secs0/60.)    : secs0=secs0 - imin1*60.
    locate 10,30: PRINT "Record end                     ";
    print using fniform$(2)+":"; ihr1;imin1; 
    print using fnfform$(5,2); secs0

  locate 12,3 :  PRINT "Sampling rate: ";sampersec;"    interval DT: "; dt;
  locate 12,50:  print "Number of points: ";npts
  locate 13,3 :  print "Instrument Constant:  "; a0; 
  locate 14,3 :  print "Flags: "; nflag(1); nflag(2); nflag(3)
  locate 16,3 :  PRINT "Number of poles: "; np; "   Number of zeroes: "; nz
return
'-------------------------------------------------------------------------
READHDR:
    nrec = 32   '---- start new record

    stn$ = ""
    FOR i = 1 TO 4: code(i) = ch4$: stn$ = stn$ + LEFT$(code(i), 1): NEXT i
    IF UCASE$(stn$) = "END " THEN GOTO 999

    phs$ = ""
    FOR i = 1 TO 4: chn(i) = ch4$: phs$ = phs$ + LEFT$(chn(i), 1): NEXT i

    alat = CVS(ch4$): alon = CVS(ch4$): elev = CVS(ch4$)
    iyr = CVL(ch4$): idoy = CVL(ch4$): ihr = CVL(ch4$): imin = CVL(ch4$)
    sec = CVS(ch4$)
    sampersec = CVS(ch4$): a0 = CVS(ch4$)
    dt=1.0/sampersec

    np = CVL(ch4$)
    FOR i = 1 TO np: FOR j = 1 TO 2: poles(j, i) = CVS(ch4$): NEXT j, i
    nz = CVL(ch4$)
    FOR i = 1 TO nz: FOR j = 1 TO 2: zeroes(j, i) = CVS(ch4$): NEXT j, i
    FOR i = 1 TO 3: nf& = CVL(ch4$): nflag(i) = nf&: NEXT i
RETURN
'-----------------------------------------------------------------------------
OUTPUT9:   '------- ASCII format output routine

'--------------- GDSN format for each station ----------
'    
OPEN "a",9, ascfil$
     CALL jdtommdd(iyr, idoy, imo, idy)
     label$ = str$(iyr)+str$(imo)+str$(idy)

'---- calculate expected reference responses and remove
     scal = 10 ^ nflag(3)
     cons = a0 * scal
     inty$=left$(phs$,2)

     if inty$="LP" or inty$="lp" then freq = 0.04 else freq = 1.0
     call gdsn(fresp(),freq,cons,nz,np,zeroes(),poles()) 

     avscal = sqr(fresp(1)*fresp(1) + fresp(2)*fresp(2))

     CalCons = cons/avscal

     stn$=LCASE$(stn$)
     phs$=LCASE$(phs$)

     secs = ihr*3600. + imin*60. + sec + (nstart-1)*dt
     ihr = int(secs/3600.) : secs=secs - ihr*3600.0
     imin=int(secs/60.)    : sec=secs - imin*60.

     p1$=lcase$(left$(phs$,3))
     if mid$(p1$,1,2)="ip" then p1$="bb"+mid$(p1$,3,1) : inst=7

     PRINT #9, USING fnaform$(4); stn$;
     PRINT #9, USING fniform$(1); inst;
     PRINT #9, USING fnaform$(3); p1$;
     PRINT #9, USING fndform$(6,4); avscal; 
     PRINT #9, USING fnaform$(12); label$;
     PRINT #9, USING fniform$(2); ihr;
     PRINT #9, USING fniform$(2); imin;
     PRINT #9, USING fnfform$(5,2); sec;
     PRINT #9, USING fnfform$(5,2); dt;
     PRINT #9, USING fnfform$(7,2); alat;
     PRINT #9, USING fnfform$(7,2); alon

     PRINT #9, USING fniform$(5); nlast-nstart+1;
     PRINT #9, USING fniform$(5); nz;
     PRINT #9, USING fniform$(5); np;
     PRINT #9, USING " "+fndform$(6,4); CalCons

     FOR i = 1 TO nz
       FOR j = 1 TO 2: PRINT #9, USING fndform$(6,4)+" "; zeroes(j, i); : NEXT j
       IF (i MOD 4) = 0 OR i = nz THEN PRINT #9,
     NEXT i

     FOR i = 1 TO np
       FOR j = 1 TO 2: PRINT #9, USING fndform$(6,4)+" "; poles(j, i); : NEXT j
       IF (i MOD 4) = 0 OR i = np THEN PRINT #9,
     NEXT i

     k=0
     FOR i = nstart TO nlast
       k=k+1
       PRINT #9, USING fndform$(6,4)+" "; seis(i);
       IF (k MOD 8) = 0 OR i = nlast THEN PRINT #9,
     NEXT i
  close 9
RETURN
'-----------------------------------------------------------------------------
FUNCTION ch4$ STATIC
SHARED nrec, a1$
DIM ch2 AS INTEGER, a4 AS STRING * 4
  IF nrec = 32 THEN GET #1, , ch2: nrec = 0
  GET #1, , a4
  nrec = nrec + 1
  ch4$ = a4
END FUNCTION
'==================================================================
SUB jdtommdd (iyr&, idoy&, imo, idy) STATIC
  shared iml()
  idy = idoy&
  kyr = iyr&
  IF (kyr MOD 4) = 0 THEN iml(2) = 29
  FOR i = 1 TO 12
    IF idy <= iml(i) THEN GOTO JDTO1
    idy = idy - iml(i)
  NEXT i
  PRINT "Error:  Julian Day too large"
  imo = 0
  EXIT SUB
JDTO1:
  imo = i
END SUB
'===================================================================
SUB rdseis (nflg%, ndat&, ss() AS SINGLE) STATIC
SHARED nrec, a1$
DIM nd2 AS INTEGER, ch2 AS STRING * 2
    '--- data is on new record
    GET #1, , ch2
    GET #1, , ndat&
    id = 1
    FOR i = 1 TO ndat&
      IF id = 32 THEN GET #1, , ch2: id = 0
      GET #1, , ss(i)
      id = id + 1
    NEXT i
    nflg% = (ndat& = 1) AND (ss(1) > 1E+09)
END SUB
'===================================================================
SUB gdsn(fresp(1),freq,cons,nz&,np&,zero(2),pole(2)) static
'---- calculates response at one circular frequency
'pii = 3.141592653
    nzero = nz&
    npole = np&
    dom=2*pii*freq
    om=dom
      sr= 0.0 : si=om
      vr1=cons : vi1=0.0

      '---- zeroes
      for jf=1 to nzero
         zr=zero(1,jf)  : zi=zero(2,jf)
         vr2=sr-zr : vi2=si-zi
         call comult (vr1,vi1,vr2,vi2,vr3,vi3)
         vr1=vr3 : vi1=vi3
      next

      '---- poles
      for jf=1 to npole
         pr=pole(1,jf) : pi=pole(2,jf)
         vr2=sr-pr : vi2=si-pi
         call cominv(vr2,vi2,vr,vi)
         call comult (vr1,vi1,vr,vi,vr3,vi3)
         vr1=vr3 : vi1=vi3
      next

      fresp(1)=vr1
      fresp(2)=vi1
'     print "GDSN: freq, fresp =";freq;fresp(1);fresp(2)
end sub
'===================================================================
GETMAX: 
'-- find maximum and average of seismogram
    seis0 = seis(nstart)
    amax = seis0: amin = seis0

    sumsi=0

    FOR i = nstart to nlast
      si = seis(i)
      sumsi=sumsi+si
      amax=fnamax(si,amax) 
      amin=fnamin(si,amin)
    NEXT i

    avg= sumsi/(nlast-nstart+1)
    ampmax = fnamax(abs(amin-avg), abs(amax-avg))
return
'======================================================================
GETHYPO:
'-- read in hypocenter data
'    print :
  khypo=kfalse:irow=0
    if inputn(irow,2,"Date of event (YYMMDD)",dt$,ktrue,ktrue,kfalse,"i",6)=27 then return
    s$=fnstrnum$(ohr):if inputn(irow+1,2,"Hour",s$,ktrue,ktrue,kfalse,"i",2)=27 then return else ohr=fnxlimit(val(s$),0,24)
    s$=fnstrnum$(omin):if inputn(irow+2,2,"Minute",s$,ktrue,ktrue,kfalse,"i",2)=27 then return else omin=fnxlimit(val(s$),0,59)
    s$=fnstrnum$(osec):if inputn(irow+3,2,"Second",s$,ktrue,ktrue,kfalse,"f",5)=27 then return else osec=fnxlimit(val(s$),0,59.99)
    s$=fnstrnum$(eplat):if inputn(irow+4,2,"Latitude (+North)",s$,ktrue,ktrue,kfalse,"f",7)=27 then return else eplat=fnxlimit(val(s$),-90,+90)
    s$=fnstrnum$(eplon):if inputn(irow+5,2,"Longitude (+East)",s$,ktrue,ktrue,kfalse,"f",8)=27 then return else eplon=fnxlimit(val(s$),-180,180)
    s$=fnstrnum$(dep):if inputn(irow+6,2,"Depth (km) ",s$,ktrue,ktrue,kfalse,"f",5)=27 then return else dep=fnxlimit(val(s$),0,999.9)
  khypo=ktrue
return
'======================================================================
sub setstep (imove) static
'-- set step for cursor motion
  if imove= 0 or imove=31 then imove=60
  if imove=21 or imove=59 then imove=30
  if imove=11 or imove=29 then imove=20
  if imove=19 then imove=10
  if imove=61 then imove=1
end sub
'==================================================================  
CHKFILES:

'-- check for existence of files

for icompo=1 to nfiles
  if infile(icompo) THEN
   finam$ = retrv$ + "." + files$(icompo)
   kflag$(icompo)="*"
   if isfile(finam$, f$, 1)=0 then kflag$(icompo)=" "
  end if
next
return

'======================================================================
rem $INCLUDE: 'M5COMFUN'
rem $INCLUDE: 'M5RDSET'
'===================================================================

