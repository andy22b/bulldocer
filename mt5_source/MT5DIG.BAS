' main MT5DIG
' program to digitize waveforms and make corrections
' written by Rob McCaffrey, MIT, 1988 as SYN3DIG
' modified by Peter Zwick, May 1990
' Revision History:
'   May,1990  Included test for COM port being avail
'     and digitizing board responding  -PZ
'   Nov,1990  Converted to SYN4DIG  -PZ
'   Jan,1991  Added user configurable digitizer boards  -PZ
'   Jan,1993  Corrected error in detecting input file existence  -PZ
'   Mar,1994  created workaround to allow parity w/8 data bits, nvar>3
'     and improve entry checking for delimited ascii format setup  -PZ
DEFINT I-N: OPTION BASE 1
CONST NSTAMAX=1 'dummy variable

rem $INCLUDE: 'M5COMSCR'
rem $INCLUDE: 'm5fndecl'
declare function isconfig%()
rem $INCLUDE: 'm5functn'
rem $INCLUDE: 'm5wigint'
rem $INCLUDE: 'm5rdset'

'-- opening screen
prgnm$="MT5DIG"
version$="1.0"
authors$="P. Zwick, R. McCaffrey, G. Abers"
call welcome

'=================================================================

 dim x1(5), y1(5), x2(5), y2(5), phz$(6), xmags(10)
 dim minpertr(3), seisscale(3), alphas(3)

'-- number of digitized points per seismogram
 ndigmax=2000
 dim t(ndigmax), a(ndigmax)
 dim t0(ndigmax), a0(ndigmax)

'-- number of interpolated points per seismogram
 nintmax=5000
 dim t2(nintmax), a2(nintmax)

 data "spz","spn","spe","lpz","lpn","lpe"
 for i=1 to 6: read phz$(i): next

nmags=8
 data 375,750,1500,3000,6000,25000,50000,100000
 for i=1 to nmags : read xmags(i) : next

'-- defaults for WWSSN seismograms
data 60, 15, 0.011142
data 30, 30, 0.005571
data 15, 60, 0.002786
 for i=1 to 3 : read minpertr(i),seisscale(i),alphas(i) : next

'-- screen plotting functions
def fnix(x) = kscrx*(x-xmin)/(xmax-xmin)
def fniy(y) = kscry/2*(1.0-0.8*y/ymax)

'-- string translation for board codes
def fnascxlat$(ain$)
static i,b$,c$,a$
  a$=ain$
  ist=1
  do
    i=instr(ist,a$,"\")
    if i>0 then
      c$=mid$(a$,i+1,1)
      if c$="." then
        a$=left$(a$,i-1)+chr$(13)+mid$(a$,i+2)
      elseif c$="\" then  ' two \ together taken as one literal \ character
        a$=left$(a$,i-1)+mid$(a$,i+1)
        ist=ist+i
      else
        b$=""
        for j=1 to 3
          c$=mid$(a$,i+j,1)
          if instr("0123456789",c$)=0 then exit for
          b$=b$+c$
        next
        a$=left$(a$,i-1)+chr$(val(b$))+mid$(a$,i+j)
      end if
    end if
  loop while i>0
  fnascxlat$=a$
end def

'-- defaults
xm2=1500 : icomp=4 : c$=phz$(icomp) :iseistype=1 :ccor=0.0
penrad=0 : iinst=1 : iorchk=ktrue
filout$="digit.dig" : kstart=ktrue
littic=0.05*kscry
indigit=ktrue

dim var$(8)
call digfile(kfalse) 'get default or saved digitizer board settings

call scpause(13,3)
101    call setscreen(0)

call scbox(1,7,ktrue,"DATA INPUT SOURCE","Press Q to QUIT")
call scprint(2,2,"KEY   OPTION")
call scprompt(3,3,"D - input from Digitizer board")
call scprompt(4,3,"F - input from File")
call scprompt(5,3,"C - Change board config from "+boardname$+" on COM"+kpt$+":")
call scprompt(6,3,"X/"):call scprompt(6,5,"Q - Shell/Quit")

      ai$=""
ik=input1(8,1,"Enter KEY for OPTION ",ai$,kfalse,kfalse,ktrue,"cdfxqz") 'include z for debugging
if ai$="q" then call endcheck: goto 101
if ai$="z" then goto 190  ' available for debugging w/o board or file
if ai$="x" then call shellcom:goto 101
if ai$="c" then ik=isconfig:goto 101

'--- setup input type, INDIGIT is the flag for type of input. If true the
'--    digitizer is being used, if false input is from file
indigit = (ai$="d")

if indigit then
  view print 8 TO 25:cls
  GOSUB diginit      'initialize digitizer
  if i=0 then close 3:view print:goto 101
else
  if filin$="" then filin$="DIGIT.TMP"
102  call scprintc(8,1," Enter name of file containing digitizer input (in DIGIT.TMP format, Esc to ABORT)")
  if inputn(9,1,"--",filin$,ktrue,ktrue,ktrue,"a",40)=27 then view print:goto 101
  if filin$="" then filin$="DIGIT.TMP"
  filin$=fnaddext$(filin$,"tmp")
  if not ifopen("I",3,filin$,0,-1) then 102
end if

'--- read station data from input file
180 if not indigit then

'-- read header line if not already read
  if not kain then
    if eof(3) then end
    line input #3, a$
  end if
  a$=lcase$(a$)
  s$=fnvalue$(1,4)
  c$=fnvalue$(5,3)
  iinst=fnvalue(8,1)
  iseistype=fnvalue(9,1)
  inithr=fnvalue(10,2)
  initmin=fnvalue(12,2)
  sec1=fnvalue(14,6)
  xm2=fnvalue(20,10)
  ccor=fnvalue(30,10)
  penrad=fnvalue(40,10)
  minpertr=fnvalue(50,10)
  seisscale=fnvalue(60,10)
  newccor=kfalse
end if

190 view print
call setscreen(0)
if kopen4 then close 4 : kopen4=kfalse
call scbox(1,22,ktrue,"STATION MENU "+string$(4,205)+" WAVEFORM DIGITIZING","Press Q to QUIT")
call scprint(2,2,"KEY  OPTION"):call scprint(2,50,"VALUE")
call scprompt(3,3,"F - File for output")
call scprompt(4,3,"H - Header")
call scprompt(5,3,"O - Order checking (y/n)")
call scprompt(6,3,"S - Station code")
call scprompt(7,3,"M - Magnification")
call scprompt(8,3,"C - Component (1spz 2spn 3spe 4lpz 5lpn 6lpe)")
call scprompt(9,3,"I - Initial point")
call scprompt(10,3,"P - Pen radius in mm")
call scprompt(11,3,"A - clock correction in milliseconds")
call scprompt(12,3,"R - instrument    1- WWSSN 15-100")
call scprint(13,8,"type         2- WWSSN 30-100")
call scprint(14,21,"3- WWSSN SP")
call scprint(15,21,"4- other")
call scprompt(16,3,"T - seismogram    1- WWSSN LP 15mm/m")
call scprint(17,8,"Type         2- WWSSN LP 30mm/m")
call scprint(18,21,"3- WWSSN SP 60mm/m")
call scprint(19,21,"4/5/6- other (dis/vel/acc calib)")
call scprompt(20,3,"D - Digitize")
call scprompt(21,3,"X/"):call scprompt(21,5,"Q - Shell/Quit")

'-- now write variables
call setvarcol
call scprint(3,50,filout$)
call scprint(4,50,header$)
call scprint(5,50,fnyesno$(iorchk))
call scprint(6,50,alltrim$(s$))
call scprint(7,50,fnstrnum$(xm2))
call scprint(8,50,phz$(icomp))
call sclocate(9,50): print using "##h##m##.##s";inithr;initmin;sec1
call scprint(10,50,fnstrnum$(penrad))
call scprint(11,50,fnstrnum$(ccor))
call scprint(12,50,fnstrnum$(iinst))
call scprint(16,50,fnstrnum$(iseistype))
call settxtcol

view print 23 TO 25:cls
w$=""
ik=input1(23,1,"Enter KEY for OPTION ",w$,kfalse,kfalse,ktrue,"fhosmcipartdxq")
if w$="q" then
  call endcheck
elseif w$="a" then
  d$=fnstrnum$(ccor)
  if inputn(23,1,"Clock correction in milliseconds (added to initial time)",_
    d$,ktrue,ktrue,ktrue,"i",6)<>27 then ccor=val(d$):newccor=kfalse
elseif w$="c" then
  d$=fnstrnum$(icomp)
  if input1(23,1,"Component",d$,kfalse,ktrue,ktrue,"123456")<>27 then
    icomp=val(d$)
    if fninlim(icomp,1,6) then c$ = phz$(icomp)
    if fninlim(icomp,1,3) then iinst=3 : iseistype=3
    if fninlim(icomp,4,6) and iinst=3 then iinst=1 : iseistype=1
  end if
elseif w$="d" then
  view print:cls
  newccor=ktrue
  goto 600
elseif w$="f" then
  f$=filout$
  if inputn(23,1,"Output file (.DIG default)",f$,ktrue,ktrue,kfalse,"a",40)<>27 then
    if f$<>"" then
      f$=fnaddext$(f$,"dig")
      if not keepold(f$) then
        ik=isfile(f$,d$,1)
        if ik=0 or ik=2 then filout$=d$
      end if
    end if
  end if
elseif w$="h" then
  ik=inputn(23,1,"Header",header$,ktrue,ktrue,ktrue,"a",12)
elseif w$="i" then
  do
    notdone=kfalse
    d$=fnstrnum$(inithr)+","+fnstrnum$(initmin)+","+fnstrnum$(sec1)
    if inputn(23,1,"Time of first point hh,mm,ss",d$,ktrue,ktrue,ktrue,"a",12)=27 then exit do
    b$=lftvar$(d$):notdone=not (istype(b$,"i") and fninlim(val(b$),0,24))
    ci$=lftvar$(d$):notdone=not (istype(ci$,"i") and fninlim(val(ci$),0,59))
    notdone=not (istype(d$,"f") and fninlim(val(d$),0,60))
    if not notdone then
      inithr=val(b$)
      initmin=val(ci$)
      sec1=val(d$)
      exit do
    end if
  loop
elseif w$="m" then
  view print:view print 9 TO 25
  call scbox(9,23,kfalse,"MAGNIFICATION MENU","")
  call scprint(10,2,"INDEX  MAGNIFICATION")
  for i=1 to nmags
    locate 11+i,2:print i;"-";
    print using fniform$(10); xmags(i);
  next
  call scprint(13+nmags,2,"Or enter mag value")
  xm2$=fnstrnum$(xm2)
  if inputn(24,2,"Magnification index or mag ",xm2$,ktrue,ktrue,ktrue,"f",8)<>27 then
    xm2=val(xm2$)
    if fninlim(xm2,1,nmags) then xm2=xmags(int(xm2))
  end if
elseif w$="o" then
  iorchk=igetyesno(23,1,"Check order of points as they are digitized?","Y",kfalse)
elseif w$="p" then
  d$=fnstrnum$(penrad)
  if inputn(23,1,"Pen radius in mm",d$,ktrue,ktrue,ktrue,"f",6)<>27 then penrad=val(d$)
elseif w$="r" then
  d$=fnstrnum$(iinst)
  if input1(23,1,"Instrument type (1-4)",d$,kfalse,ktrue,ktrue,"1234")<>27 then iinst=val(d$)
elseif w$="s" then
  if inputn(23,1,"Station code ",s$,ktrue,ktrue,ktrue,"a",4)<>27 then 
    s$=lcase$(left$(s$+"    ",4))
    if newccor then ccor=0
  end if
elseif w$="t" then
  d$=fnstrnum$(iseistype)
  if input1(23,1,"Seismogram type (1-6)",d$,kfalse,ktrue,ktrue,"123456")<>27 then iseistype=val(d$)
elseif w$="x" then
  view print:locate 22,80
  call shellcom
end if
goto 190

'-- start digitizing routine
600

'-- get origin time in seconds
osec=inithr*3600. + initmin*60. + sec1

'-- get scaling parameters
if fninlim(iseistype,1,3) then
  minpertr=minpertr(iseistype) : seisscale=seisscale(iseistype) : alphas=alphas(iseistype)
elseif iseistype>3 and indigit then
  do
    d$=fnstrnum$(minpertr)+","+str$(seisscale)
    call scprintc(1,1,"Enter number of minutes per trace and original seismogram scale in mm/min")
    if inputn(2,1,"--",d$,ktrue,ktrue,kfalse,"a",25)=27 then exit do
    b$=lftvar$(d$)
    if not istype(b$,"i") then
      b$=displaymsg$("ERROR: minutes per trace not an integer","",kfalse) 
    elseif not istype(d$,"f") then 
      b$=displaymsg$("ERROR: seismogram scale not a floating point value","",kfalse) 
    else
      minpertr=val(b$):seisscale=val(d$):exit do
    end if
  loop
  print
end if

'-- write station stuff to DIGIT.TMP file, but only if the digitizer is being
'--  used
if not indigit then 610
open "a",4,"digit.tmp"
kopen4=ktrue
if not kstart then print #4,""         'blank line between seismograms
print #4, using fnaform$(4); s$;
print #4, using fnaform$(3); c$;
print #4, using fniform$(1); iinst; iseistype;
print #4, using fniform$(2); inithr;initmin;
print #4, using fnfform$(6,2); sec1;
print #4, using fnfform$(10,1); xm2; ccor; penrad;
print #4, using fnfform$(10,3); minpertr;seisscale
kstart=kfalse

print "Digitize 6 seismogram scaling coordinates as follow:" : print
print "1 - 2     Two points far apart on a single trace"
print "  3       Point on trace below the trace containing pts. 1 & 2"
print "4 - 5     Two minute marks - one minute apart"
print "  6       Minute mark farther along same trace from pts. 4 & 5"
print : print

610 ksave=kfalse
for i=1 to ndigmax : t0(i)=0 : a0(i)=0 : next i
for i=1 to ndigmax : t(i)=0 : a(i)=0 : next i
skewalpha=0.0

'-- read in six points for orientation and scaling
if indigit then
  print "Digitize point 1 ";
  call digitxy (x0,y0,1,iend)
else
  call filexy(x0,y0,1,iend)
end if
if iend then 180

for i=1 to 5
  if indigit then
    print "Digitize point"; i+1;
    if i<>3 then print
    call digitxy(x,y,i+1,iend)
  else
    call filexy(x,y,n,iend)
  end if
  if iend then 180

'-- make the first point the origin and translate the other 5
  x1(i)=x-x0
  y1(i)=y-y0
next : print :print

'-- THETA is the angle the trace makes with the table X-axis
if x1(1)=0 then
  beep
  d$=displaymsg$("Digitizing Error","",kfalse)
  goto 180
end if
theta=atn(y1(1)/x1(1))

'-- rotate points 2-6 about theta
for i=1 to 5
  x2(i)=x1(i)*cos(theta)+y1(i)*sin(theta)
  y2(i)=-x1(i)*sin(theta)+y1(i)*cos(theta)
next

' UNITPERMIN is table units per minute
'-- first get the number of minutes between points 5 and 3
if (x2(4)=x2(3)) or (x2(5)=x2(3)) then
  beep
  d$=displaymsg$("Digitizing Error 1 on points 4-6","",kfalse)
  goto 180
end if

minutes=int((x2(5)-x2(3))/(x2(4)-x2(3))+.5)
if minutes=0 then
  beep
  d$=displaymsg$("Digitizing Error 2 on points 4-6","",kfalse)
  goto 180
end if

'-- divide distance between pts 5 and 3 by MINUTES to get scale
unitpermin=(x2(5)-x2(3))/minutes

' SEISUNITS is the length of the seismogram in table units
seisunits=minpertr*unitpermin

' SSSCALE is the scale in table units per mm
ssscale=unitpermin/seisscale

' ALPHAS - angle of trace wrt true horizontal on `gram
'  if WWSSN record use the theoretical value, otherwise calculate
if iseistype>3 then alphas=atn(abs(y2(2))/seisunits)

'-- digitize the starting point for the waveform, all successive times and
'   amplitudes are determined relative to this point
if indigit then
  print "Digitize Initial Point  at time= ";
  print using fniform$(2); inithr;initmin;
  print using fnfform$(6,2); sec1
  call digitxy(xinit,yinit,7,iend)
else
  call filexy(xinit,yinit,n,iend)
end if
if iend then 180

'--- start digitizing waveform
npoint=0 : nptmp=0

'-- If IEXCHK is on, a warning sounds while digitizing the LPE component when
'   its time exceeds the maximum time of the previous LPN component

if (indigit and icomp=6 and lpnprev and s$=sprev$) then iexchk=ktrue else iexchk=kfalse

if indigit then beep : print "Digitize waveform "

900  ' start loop to digitize waveform
while npoint<ndigmax

  if indigit then
    npoint=npoint+1
    call digitxy(x,y,npoint,iend)
    nptmp=npoint-1
  else
    call filexy(x,y,npoint,iend)
    if npoint=0 and iend then 180
    nptmp=fnamax(nptmp,npoint)
  end if
  if iend then 1110

  j=npoint

960  GOSUB ROTATE  'change coordinates for point J

  if indigit and iexchk and t0(npoint)>prevmax then
    beep
    print " Exceeded previous record length"
    iexchk=kfalse
  end if

'-- check for time going backward
  if indigit then
    t=t0(npoint)+osec
    call hrminsec( t, ih, im, s )
'-- print time and amplitude of this point on screen
    print "Time =";
    call setvarcol
    print using fniform$(3); ih; im;
    print using fnfform$(6,2);s;
    call settxtcol
    print " Amplitude =";
    call setvarcol
    print using "#######"; a0(npoint);
    call settxtcol
    if npoint>1 and iorchk and (t0(npoint)<t0(npoint-1)) then
      beep
      print "  <--":irow=0
      if not igetyesno(irow,1,"Time going backward, keep above point?","N",ktrue) then
        npoint=npoint-1
        call cline((irow-1),1)
      end if
    else
      print
    end if
  end if

wend

1080 print chr$(7);chr$(7);"******WARNING******"
print "Maximum number of data points reached"

1110 npoint=nptmp
call SORTDAT(t0(),a0(),npoint)  'sort data chronologically
for i=1 to npoint : a(i)=a0(i) : t(i)=t0(i) : next

'-- print plot/edit menu
1112 call setscreen(0)
klst$="perscn"
call scbox(1,11,ktrue,"PLOT/EDIT MENU","Press Esc to ABORT")
call scprint(2,10,header$+" "+s$+" "+c$)
call scprint(3,2,"KEY  OPTION")
call scprompt(4,3,"P - Plot seismogram")
call scprompt(5,3,"E - Edit seismogram")
call scprompt(6,3,"R - Rotate, angle ="+str$(skewalpha))
call scprint(7,3,"A - Add new points")
if indigit then call scprintc(7,3,"A"):klst$=klst$+"a"
call scprompt(8,3,"S - Store data in file")
call scprompt(9,3,"C - Change interpolation interval DT; current value ="+str$(dtplot))
call scprompt(10,3,"N - New seismogram (return to STATION MENU)")

w$=""
ik=input1(12,1,"Enter KEY for OPTION ",w$,kfalse,ktrue,ktrue,klst$)
if w$="p" then
  IW1=1 : dt=dtplot : goto 1240
elseif w$="e" then
  IW1=2 : dt=0 : m=1 : goto 1240
elseif w$="a" then
  goto 1235
elseif w$="r" then
  gosub RMSKEW
elseif w$="s" then
  goto 1720
elseif w$="n" or ik=27 then
  goto 1180
elseif w$="c" then
  d$=fnstrnum$(dplot)
  if inputn(12,1,"New DT",d$,ktrue,ktrue,ktrue,"f",6)<>27 then dtplot=val(d$)
end if
goto 1112

1180 if ksave THEN 180
beep : if igetyesno((0),(0),"Data not saved ****** Do you want new seismogram?","N",kfalse) then 180
GOTO 1112

1235 '--- add points
print: print "Time of last point ";
t=t(npoint)+osec
call hrminsec( t, ih, im, s )
print using fniform$(2); ih;im; : print using fnfform$(6,2);s
print " Amp.= ";  : print using fniform$(6); a(npoint)
print "Digitize new points"
goto 900

1240  ' ---- Plotting section
if npoint<2 then print "Not enough points to plot": d$=input$(1) : goto 1112

call setscreen (kega)

1250 xmin = t(1) : xmax = t(npoint)
ymax = 0 : for i=1 to npoint : ymax=fnamax(abs(a(i)),ymax) : next

1330 call scclear
call scprint (25,33,"Tic=5s")
call scprint (1,1,"") : PRINT Header$;" ";S$;" ";C$;" NPTS=";npoint;
call scprint (2,1,"") : PRINT "Tmin=";int(xmin);" Tmax=";int(xmax);

'--- make time axis with dots every 5s
for i=0 to t(npoint) step 5 : call scpoint( fnix(i), fniy(0),kgrcol(kcolsyn)) : next

'--- make spline curve, detrend, and plot
if dt>0 then
  call m5wigint (t(1),t(npoint),t(),a(),npoint,dt,a2(),nout)
  call m5detrnd (nout,dt,a2(),a,b)
  ymax = 0
  for i=1 to nout
    ymax=fnamax(abs(a2(i)),ymax)
    t2(i)=t(1)+(i-1)*dt
  next

  print " Amax=";int(ymax);
  call scpoint (fnix(t2(1)),fniy(a2(1)),kgrcol(kcoldat))
  for k=2 to nout : call scline1 (fnix(t2(k)),fniy(a2(k)),kgrcol(kcoldat),kfalse) : next
  y$=input$(1) : goto 1112
end if

'--- if dt=0 then just connect points
for i=1 to npoint : t2(i)=t(i) : a2(i)=a(i) : next
ymax = 0 : for i=1 to npoint : ymax=fnamax(abs(a2(i)),ymax) : next
print " Amax=";int(ymax);
call scpoint (fnix(t2(1)),fniy(a2(1)),kgrcol(kcoldat))
for i=2 to npoint : call scline1 (fnix(t2(i)),fniy(a2(i)),kgrcol(kcoldat),kfalse): next
if iw1=1 then y$=input$(1) : goto 1112

1544 call scprintc (24,1,"EDIT: 1-Quit 2-reDigitize")
call scprintc (25,1,"3-reMove 4-Left 6-Right")
inc=0
d$=input$(1) : k=asc(d$)
if fninlim(k,49,56) then
  inc=val(d$)
else
  d$=lcase$(d$)
  if d$="q" then inc=1
  if d$="d" then inc=2
  if d$="m" then inc=3
  if d$="l" then inc=4
  if d$="r" then inc=6
end if

if inc=0 then 1544
if inc=1 then
'       gosub rmskew
  go to 1112
end if
on inc gosub 2100,7000,6700,6000,2100,6000
if inc=2 then 1250
if inc=3 then 1240
goto 1544  ' while moving cursor
goto 1112

1720  ' ----- Output section
1725 if not ksave then 1735
1730 beep
if igetyesno((0),(0),"Data already stored * store again?","N",kfalse) then 1735
goto 1112

1735 open "a",1,filout$
ksave=ktrue

'-- set up name and type of previous seismogram
prevmax=t(npoint)
if c$="lpn" then lpnprev=ktrue else lpnprev=kfalse
if lpnprev then sprev$=s$ else sprev$="    "

'-- put out a instrument=9 for other than WWSSN
if iinst=4 then inst=9 else inst=iinst

print #1, using fnaform$(4); s$;
print #1, using fniform$(1); inst;
c2$=c$
if iseistype>4 then mid$(c2$,2,1)=mid$("va",iseistype-4,1)
print #1, using fnaform$(3); c2$;
print #1, using fniform$(10); xm2;
print #1, using fnaform$(12); header$;
print #1, using fniform$(2); inithr;initmin;
print #1, using fnfform$(5,2); sec1+ccor/1000.
print #1, using fniform$(5); npoint

for i=1 to npoint
  print #1, using fniform$(8); t(i)*100; a(i);
  if (i mod 5 =0) or (i=npoint) then print #1, ""
  if i=npoint then 1965
next
1965 close 1 : goto 1112

' ********************************
'   SUBROUTINES
' ********************************
'==============================================================
2100 return  'dummy return
'==============================================================
ROTATE:
'---Rotate point about Theta, correct for slant of trace and curvature of pen
'-- return t0(j), a0(j) = corrected time and amplitude

'--- rotate about THETA
x3=(x-xinit)*cos(theta)+(y-yinit)*sin(theta)
y3=-(x-xinit)*sin(theta)+(y-yinit)*cos(theta)

'--- correct time and amplitude for non-orthogonality between these two axes.
'--   ALPHAS is the angle between the amplitude axis and the normal to the
'--   time axis.
d=y3*tan(alphas)
a0(j)=sqr(y3*y3+d*d)/ssscale*sgn(y3)
t0(j)=60*(x3+d)/unitpermin

'--- correct for arcing of pen
if penrad>0 then
  asign=sgn(a0(j))
  beta=asign*fnasin(abs(a0(j))/penrad)   'deflection angle
  t0(j) = t0(j)- penrad * (1.0 - cos(beta)) * (60./seisscale)
  a0(j) = penrad * beta  ' = r*beta
end if

a0(j)=a0(j)*1000.  'convert amplitude to microns from mm

return
'==============================================================
sub SORTDAT (t(1), a(1), npoint) static
'--   sort data chronologically

for i=1 to npoint
  tmin=100000.0
  for j=i to npoint
    if t(j) < tmin then tmin=t(j) : imin=j
  next
  swap t(i),t(imin) : swap a(i),a(imin)
next

'-- make sure 2 points do not have the same time
for i= 1 to npoint-1
  if t(i)=t(i+1) then t(i+1)=t(i+1)+.01
next

end sub
'==============================================================
6000 '-----subroutine to move cursor
'-- erase old marker
call scline2 (fnix(t2(m)),fniy(a2(m))+littic,fnix(t2(m)),fniy(a2(m))+1,0,kfalse)
'-- find new position
if inc=4 then m=m-1 else m=m+1
m=fnxlimit(m,1,npoint)
'-- draw new mark
call scline2 (fnix(t2(m)),fniy(a2(m))+littic,fnix(t2(m)),fniy(a2(m))+1,kgrcol(kcoldat),kfalse)
'-- print time and amplitude
call cline(2,1)
call scprint (2,1,"Current point: ")
t=t2(m)+osec
call hrminsec( t, ih, im, s )
print "T=";
print using fniform$(2); ih;im;
print using fnfform$(6,2);s;
print " A=";
print int(a2(m));
return
'==============================================================
6700 '----- remove point M
print
npoint=npoint-1
for i=m to npoint
  t(i)=t(i+1)
  a(i)=a(i+1)
  t0(i)=t0(i+1)
  a0(i)=a0(i+1)
next
return
'==============================================================
7000 ' ----- redigitize a point M
call cline (25,1)
if not indigit then
  call scprint (25,1,"Cannot replace point in file input mode")
  d$=input$(1)
else
  call scprint (25,1,"Digitize point")
  call digitxy(x,y,m,iend) :j=m: gosub ROTATE
  if skewalpha<>0. then
    skew=fnrad(skewalpha)
    secpermic=60./(seisscale*1000.)   'seconds per micron scale
    d=a0(j)*tan(skew)
    a(j)=sqr(a0(j)*a0(j)+d*d)*sgn(a0(j))
    t(j)=t0(j)+d*secpermic
  else
    t(j) = t0(j) : a(j) = a0(j)
  end if
end if
return
'==============================================================
sub filexy (x,y,npt,iend) static
shared a$, kain

'----- read X,Y from file
'--    NPT is the point number being read

'-- set IEND flag to false
iend=0 : kain=0

'-- if at end-of-file set IEND flag to true and exit
F1:  if eof(3) then iend=-1 : exit sub

'-- read a line from file
line input #3, a$

'-- if the line is blank set IEND flag to true and exit
if len(a$)=0 then iend =-1 : exit sub

'-- if the line starts with a letter set IEND flag and flag KAIN that the
'--  header has been read in
k=asc(fnvalue$(1,1))
if fninlim(k,65,90) or fninlim(k,97,122) then iend =-1 : kain=-1 : exit sub

'-- get X, Y and the point number
x=fnvalue(1,8) : y=fnvalue(9,8) : npt=fnvalue(17,5)

'-- if the point number is 0 skip this point
if npt=0 then goto F1

end sub
'=================================================================
sub hrminsec ( sectime, ihr, imin, secs ) static
secs=sectime
ihr = int(secs/3600.)
secs=secs - ihr*3600.0
imin=int(secs/60.)
secs=secs - imin*60.
end sub
'======================================================================
RMSKEW:
'--- correct time and amplitude for skew angle SKEWALPHA
d$=fnstrnum$(skewalpha)
if inputn((0),1,"Enter angle for rotation (+ clockwise)",_
    d$,ktrue,ktrue,ktrue,"f",12)<>27 then skewalpha=val(d$)

skew=fnrad(skewalpha)
secpermic=60./(seisscale*1000.)   'seconds per micron scale

for i=1 to npoint
  d=a0(i)*tan(skew)
  a(i)=sqr(a0(i)*a0(i)+d*d)*sgn(a0(i))
  t(i)=t0(i)+d*secpermic
next
call sortdat( t(), a(), npoint)
return
'======================================================================
sub digitxy (x,y,npt,iend) static
shared kopen4,nbend,nvar,kx,ky,kb,kxl,kyl,kbl,ibrdfmt,var$()

'----- read point NPT, coords X,Y from digitizer
'----- this is set up for the Summagraphics MM1812

'-- set IEND flag to false
iend = 0
'-- clear digitizer buffer first
WHILE LOC(3)>0
  cm$=INPUT$(loc(3),3)
WEND
'-- read X and Y and the button number from digitizer
WHILE inkey$<>CHR$(27) AND EOF(3):WEND
IF EOF(3) THEN
  nb=nbend
ELSE
  if ibrdfmt=1 then
    for j=1 to nvar
      input #3, var$(j)  'x,y,nb
    next
    x=val(var$(kx)):y=val(var$(ky)):nb=val(var$(kb))
  else
    line input #3, ai$
    x=val(mid$(ai$,kx,kxl)):y=val(mid$(ai$,ky,kyl)):nb=val(mid$(ai$,kb,kbl))
  end if
  sound 1000,.15
END IF

'-- if the button is 4 then set the flag IEND to true (-1)
if nb=nbend then iend=-1

'-- write raw digitizer data to DIGIT.TMP
if kopen4 and (not iend) then
  print #4, using fniform$(8); x;y;
  print #4, using fniform$(5); npt
end if

end sub
'=================================================================
DIGINIT:
'--- this was originally set up for the Summagraphics MM1812 on serial port

irow=0
if msg$<>"" then print fnascxlat$(msg$)
call scpause(irow,1)

9002 kport = -1
close 3      '-- make sure filenumber 3 not already open
i = 0
kport=val(kpt$)
kp$ = "COM" + kpt$
DEF SEG = 0
i = PEEK(&H400 + 2 * (kport - 1)) + &H100 * PEEK(&H401 + 2 * (kport - 1))
DEF SEG
IF i = 0 THEN
  if input1(irow,1,"Serial port "+ kp$+ " not available. Enter new COM port (1-4, Esc to ABORT)",_
    kpt$,kfalse,ktrue,ktrue,"1234")=27 then return else goto 9002
END IF
print "testing ";kp$;"..."
call opencom(i,kp$,baud$,spar$) 'moved open command into sub -PZ 3/3/94

'-- send Reset command
call digitcmd(reset$)

'-- set digitizer mode
call digitcmd(init$)

'-- send self test command
call digitcmd(selftest$)

if selftest$+response$<>"" then
'Proper response to self test is "O". Min delay for test is about .3 seconds
  i=0
  cm$ = ""
  asec = TIMER
  WHILE ABS(TIMER - asec) < 2 AND cm$ <> fnascxlat$(response$)
    j=LOC(3)
    IF j>0 THEN       'get response
      if ibrdfmt=1 then
        cm$=cm$+INPUT$(j,3)
      else
        line input #3,ci$
        if len(cm$)>0 then cm$=cm$+chr$(13)+ci$ else cm$=ci$
      end if
      i=j
    end if
  WEND

'-- Keep going if self test OK
  IF selftest$<>"" and i = 0 THEN
    b$=displaymsg$("Board did not respond to self test command using "+kp$,"",kfalse)
    return
  ELSEIF i>0 and cm$ <> fnascxlat$(response$) THEN
    i = ASC(cm$)
    PRINT
    IF INSTR(ucase$(boardname$),"SUMMAGRAPHICS")>0 then
      PRINT "Digitizer board failed the self test (Returned byte = &H"; HEX$(i); ")"
      print
      PRINT "The tests that failed are:"
      PRINT "(valid only for Summagraphics MM Series digitizing tablets)"
      PRINT
      IF (i AND 1)=0 THEN PRINT "analog circuitry"
      IF (i AND 2)=0 THEN PRINT "cursor connection &/or coil operation"
      IF (i AND 4)=0 THEN PRINT "digital circuitry & ROM validity"
      IF (i AND 8)=0 THEN PRINT "cursor detection (put cursor on board)"
      IF (i AND 256)<>0 THEN PRINT "parity"
      PRINT
      PRINT "Press I to Ignore test results, or any other key to try again..."
      DO
        ik$=INKEY$
      loop WHILE ik$ = ""
      IF not ucase$(ik$)="I" then i=0:return
    else
      if selftest$<>"" then_
        PRINT "Digitizer board response to self test not same as answer expected."
      print "In the following response from board, \xxx=unprintable ASCII xxx character"
      print "(where \.=CR, \\=\, and \32 =leading or trailing space):"
      ik=len(cm$)
      for j=1 to ik
        ic=asc(mid$(cm$,j,1))
        SELECT CASE ic
          CASE 33 TO 91,93 to 126
            print chr$(ic);
          CASE 92
            print "\\";
          CASE 32
            if j=1 or j=ik then
              print "\32";
            else
              print " ";
            end if
          CASE 13
            print "\.";
          CASE else
            print "\"+fnstrnum$(ic);
        END SELECT
      next
      print
    end if
  END IF
end if
i=0
'--- test the board by entering points and printing them on the screen
call setchccol
PRINT "**** Always use Button";nbend;"or press ESC to end input ****"
print "Test digitizer board position data using other buttons:"
x=0:y=0
print "X","Y"
call setvarcol
DO
  CALL digitxy(x, y, n, iend)
  i=i or x<>0 or y<>0
  if i then PRINT x, y
LOOP UNTIL iend
if i=0 then d$=displaymsg$("No position data received from board","",kfalse)
RETURN
'=================================================================

function isconfig
shared kx,ky,kb,kxl,kyl,kbl,boardname$,nbend,ibrdfmt,nvar
shared reset$,init$,selftest$,response$,spar$,baud$,irow,msg$,kpt$

do
  isconfig=kfalse
  klst$="ncbfritamsu"
  call scbox(1,16,ktrue,"BOARD CONFIGURATION","Press Esc to ABORT")
  call scprint(2,3,"KEY  OPTION"):call scprint(2,55,"VALUE")
  call scprompt(3,4,"N - Name of digitizing board")
  call scprompt(4,4,"C - COM port number")
  call scprint(5,10,"COM port settings")
  call scprompt(6,4,"B - Button number for ending input")
  call scprompt(7,4,"F - Format of board output string")
  if ibrdfmt=1 then
    call scprint(8,10,"Number of variables in board output string")
    call scprint(9,10,"Positions of X,Y,& Button# in output string")
  elseif ibrdfmt=2 then
    call scprint(8,10,"Output string fields (start/length):")
    call scprint(9,35,"X="):call scprint(9,45,"Y="):call scprint(9,55,"Button#=")
  endif
  call scprompt(10,4,"R - Reset command")
  call scprompt(11,4,"I - Initialize command")
  call scprompt(12,4,"T - self Test command")
  call scprompt(13,4,"A - Answer from self test")
  call scprompt(14,4,"M - Message to display")
  call scprompt(15,4,"S/"):call scprompt(15,6,"U - Save/Use changes")

  call setvarcol
  call scprint(3,55,boardname$)
  call scprint(4,55,"COM"+kpt$+":")
  call scprint(5,55,baud$+","+spar$)
  call scprint(6,55,fnstrnum$(nbend))
  if ibrdfmt=1 then
    call scprint(7,55,"delimited ASCII")
    call scprint(8,55,fnstrnum$(nvar))
    kfmt$=fnstrnum$(kx)+","+fnstrnum$(ky)+","+fnstrnum$(kb)
    call scprint(9,55,kfmt$)
  elseif ibrdfmt=2 then
    call scprint(7,55,"fixed field ASCII")
    call scprint(9,37,fnstrnum$(kx)+"/"+fnstrnum$(kxl))
    call scprint(9,47,fnstrnum$(ky)+"/"+fnstrnum$(kyl))
    call scprint(9,63,fnstrnum$(kb)+"/"+fnstrnum$(kbl))
  endif
  call scprint(10,55,reset$)
  call scprint(11,55,init$)
  call scprint(12,55,selftest$)
  call scprint(13,55,response$)
  call scprint(14,55,msg$)

  ai$="":irow=17
  if input1(irow,1,"Enter KEY for OPTION ",ai$,kfalse,ktrue,ktrue,klst$)=27 then
    if igetyesno(irow,1,"ABANDON changes and return to startup configuration?","N",ktrue) then
      call digfile(kfalse)
      exit do
    end if

  elseif ai$="n" then
    ik=inputn(irow,1,"Enter board name",boardname$,ktrue,ktrue,ktrue,"a",30)
    if boardname$="" then boardname$="noname"
  elseif ai$="b" then
    d$=fnstrnum$(nbend)
    if inputn(irow,1,"Enter number transmitted by button used to end input",d$,ktrue,ktrue,ktrue,"i",3)<>27 then nbend=val(d$)  'ending button number

  elseif ai$="f" then
    d$=fnstrnum$(ibrdfmt)
    if input1(irow,1,"Select (1) delimited ASCII, or (2) fixed field position ASCII",d$,kfalse,ktrue,kfalse,"12")<>27 then
      ibrdfmt=val(d$)  'ending button number
      if ibrdfmt=1 then
  '--number of variables from digitizing board for one button push
        d$=fnstrnum$(nvar)
        if input1(irow+1,1,"Enter number of variables in transmitted string (3-8)",d$,ktrue,ktrue,kfalse,numlist$(3,8))<>27 then nvar=val(d$)
        chcs$=numlist$(1,nvar)
        d$=fnstrnum$(kx)
        if input1(irow+2,1,"Enter position of X in string (1-"+fnstrnum$(nvar)+")",d$,ktrue,ktrue,kfalse,chcs$)<>27 then kx=val(d$)
        mid$(chcs$,kx,1)=chr$(0) 'blank out kx entry from choices
        d$=fnstrnum$(ky)
        if input1(irow+3,1,"Enter position of Y in string (1-"+fnstrnum$(nvar)+")",d$,ktrue,ktrue,kfalse,chcs$)<>27 then ky=val(d$)
        mid$(chcs$,ky,1)=chr$(0) 'blank out ky entry from choices
        d$=fnstrnum$(kb)
        if input1(irow+4,1,"Enter position of Button# in string (1-"+fnstrnum$(nvar)+")",d$,ktrue,ktrue,kfalse,chcs$)<>27 then kb=val(d$)
        if ky=kx or kb=kx or kb=ky then
          b$=displaymsg$("ERROR - Positions must be unique","",kfalse)
        end if
      elseif ibrdfmt=2 then
        do
          d$=fnstrnum$(kx)
          if inputn(irow+1,1,"Character start position of X field in string",d$,ktrue,ktrue,kfalse,"i",3)<>27 then kx=val(d$)
          d$=fnstrnum$(kxl)
          if inputn(irow+2,1,"Length of X field",d$,ktrue,ktrue,kfalse,"i",3)<>27 then kxl=val(d$)
          d$=fnstrnum$(ky)
          if inputn(irow+3,1,"Character start position of Y field in string",d$,ktrue,ktrue,kfalse,"i",3)<>27 then ky=val(d$)
          d$=fnstrnum$(kyl)
          if inputn(irow+4,1,"Length of Y field",d$,ktrue,ktrue,kfalse,"i",3)<>27 then kyl=val(d$)
          d$=fnstrnum$(kb)
          if inputn(irow+5,1,"Character start position of Button# field in string",d$,ktrue,ktrue,kfalse,"i",3)<>27 then kb=val(d$)
          d$=fnstrnum$(kbl)
          if inputn(irow+6,1,"Length of Button# field",d$,ktrue,ktrue,kfalse,"i",3)<>27 then kbl=val(d$)
          ng=kfalse
          for i=kx to kx+kxl-1
            for j=ky to ky+kyl-1
              if i=j then ng=ktrue
              for k=kb to kb+kbl-1
                if k=i or k=j then ng=ktrue
                if ng then exit for
          next k,j,i
          if ng then b$=displaymsg$("ERROR - fields overlap","",kfalse)
        loop while ng
      end if
    end if

  elseif instr("rita",ai$)>(0) then
    call scprint(irow,1,"Use \xxx to represent ascii code xxx for untypeable characters")
    call scprint(irow+1,1,"(\.=CR, \\=\, \32 =leading or trailing space")
    if ai$="a" then print ")"; else print ",  \dxx=xx/10 seconds delay time)";
    if ai$="r" then
      ik=inputn(irow+2,1,"Enter board reset command",reset$,ktrue,ktrue,ktrue,"a",40)
    elseif ai$="i" then
      ik=inputn(irow+2,1,"Enter board initialization command",init$,ktrue,ktrue,ktrue,"a",40)
    elseif ai$="t" then
      ik=inputn(irow+2,1,"Enter board self test command",selftest$,ktrue,ktrue,ktrue,"a",40)
    elseif ai$="a" then
      ik=inputn(irow+2,1,"Board should answer self test with:",response$,ktrue,ktrue,ktrue,"a",40)
    end if

  elseif ai$="c" then
        do
          if input1(irow,1,"COM port number (1 or 2)",kpt$,kfalse,ktrue,kfalse,"12")=27 then exit do
          if inputn(irow+1,1,"COM port baud rate (300,1200,2400,4800,9600,19200)",baud$,ktrue,kfalse,kfalse,"i",5)=27 then exit do
          if left$(baud$,2)="19" then  'only allow specified baud  -PZ 3/3/94
            baud$="19200"
          else
            select case left$(baud$,1)
              case "1"
                baud$="1200"
              case "2"
                baud$="2400"
              case "3"
                baud$="300"
              case "4"
                baud$="4800"
              case else
                baud$="9600"  'default
            end select
          end if
          d$=mid$(spar$,1,1)
          if input1(irow+2,1,"Parity (N)one, (O)dd, (E)ven, (S)pace, (M)ark",d$,ktrue,ktrue,kfalse,"noesm")=27 then exit do
          mid$(spar$,1)=ucase$(d$)
          d$=mid$(spar$,3,1)
          if input1(irow+3,1,"Number of data bits (5-8)",d$,ktrue,ktrue,kfalse,"5678")=27 then exit do
          mid$(spar$,3)=d$
          b$=mid$(spar$,5)
stpbit:
          d$=lftvar$(b$)
          if inputn(irow+4,1,"Number of stop bits (1,1.5,2)",d$,ktrue,ktrue,kfalse,"f",3)=27 then exit do
          if d$="1" or d$="1.5" or d$="2" then spar$=left$(spar$,4)+d$ else goto stpbit
          d$=b$
          call scprint(irow+5,1,"Options: ASC (incl XON/XOFF), BIN (binary), OPx (x=ms wait for active line)")
          call scprint(irow+6,10,"CDx,CSx,DSx (x=ms timeout on DCD,CTS,DSR lines), RS (no RTS detect)")
          call scprint(irow+7,10,"(default x is 1000 ms except for OP where default is infinite wait)")
          if inputn(irow+8,1,"Enter Options as comma-delimited string",d$,ktrue,ktrue,kfalse,"a",35)=27 then
            if len(b$)>0 then spar$=spar$+","+b$
          else
            if len(d$)>0 then spar$=spar$+","+ucase$(d$)
          end if
          exit do
        loop
  elseif ai$="m" then
    print "Current message:"
    print"-->"+msg$
    call scprintc((0),1,"Enter message to display to user on startup ("+string$(2,34)+"=above msg, \.=CR, \\=\):"):print
    line input;"-->";d$ 
    do
      ik=instr(d$,string$(2,34))
      if ik>0 then d$=left$(d$,ik-1)+msg$+mid$(d$,ik+2)
    loop while ik>0
    msg$=d$

  elseif ai$="s" then
    if igetyesno(irow,1,"Save all changes to \DIG.SET now?","Y",ktrue) then
      if isfile("\dig.set",d$,1)=0 then
        d$="C"
        ik=input1(irow,1,"\DIG.SET exists - (R)eplace,(B)ackup,or (C)ancel?",d$,kfalse,ktrue,ktrue,"rbc")
        if ik=27 then d$="c"
        if d$="b" then shell "copy \dig.set \dig.sbk >nul"
      end if
      if d$<>"c" then call digfile(ktrue):isconfig=ktrue:exit do
    end if

  elseif ai$="u" then
    if igetyesno(irow,1,"Use changes for this session only?","Y",ktrue) then
      isconfig=ktrue
      exit do
    end if

  end if
loop
cls:irow=0
end function

sub digfile(isave)
shared kx,ky,kb,kxl,kyl,kbl,boardname$,nbend,ibrdfmt,nvar
shared reset$,init$,selftest$,response$,spar$,baud$,msg$,kpt$

if isave then
  open "o",1,"\dig.set"
  print #1,boardname$
  print #1,kpt$
  print #1,baud$+","+spar$
  print #1,nbend
  print #1,ibrdfmt
  write #1,nvar,kx,ky,kb
  if ibrdfmt=2 then write #1,kxl,kyl,kbl
  print #1,reset$
  print #1,init$
  print #1,selftest$
  print #1,response$
  print #1,msg$
  close 1
else
'-- default digitizer constants
  boardname$="Summagraphics MM1812"
  kpt$="1"
  baud$="9600"
  spar$ = "E,7,2,CS,DS"
  nbend=4  'ending button number
  nvar=3   'number of variables from digitizing board for one button push
  kx=1:ky=2:kb=3  'position in nvar for x,y, and button number values
  kxl=5:kyl=5:kbl=5  'default lengths of fixed fields
  reset$= "\32\d2\0\.\d5" ' "+ CHR$(0)   ' Autobaud detection followed by Reset command
  init$="Bj"  ' B = point mode, j= 1000 lpi resolution (highest)
  selftest$="tw\.\d3"  ' t = self test, w = send test results
  response$="O" 'proper response to init
  ibrdfmt=1
  msg$="Make sure cursor is placed on board"
  if isfile("\dig.set",d$,1)=0 then
    open "i",1,"\dig.set"
    input #1,boardname$:if boardname$="" then boardname$="noname"
    input #1,kpt$
    line input #1,spar$
    baud$=lftvar$(spar$)
    input #1,nbend
    input #1,ibrdfmt
    input #1,nvar,kx,ky,kb
    if ibrdfmt=2 then input #1,kxl,kyl,kbl
    input #1,reset$
    input #1,init$
    input #1,selftest$
    input #1,response$
    line input #1,msg$
    close 1
  end if
end if
end sub

sub digitcmd(cmdin$)
cmd$=cmdin$
do
  idelay=1 '0.1 sec default delay
  i=1
  do
    notdone=kfalse
    ik=instr(i,lcase$(cmd$),"\d")
    if ik>0 then
      do
        isl=instr(i,cmd$,"\\")
'--eliminate literal \
        if isl>0 and isl<ik then i=isl+2 :notdone=ktrue else exit do
      loop
    end if
  loop while notdone
  if ik>0 then 'delay specified
    c$=left$(cmd$,ik-1):cmd$=mid$(cmd$,ik+2)
    b$=""
    for j=1 to 2
      b1$=mid$(cmd$,j,1)
      if instr("0123456789",b1$)=0 then exit for
      b$=b$+b1$
    next
    cmd$=mid$(cmd$,j)
    idelay=val(b$)
  else
    c$=cmd$:cmd$=""
  end if
  c$=fnascxlat$(c$)
  PRINT #3, c$;
  asec = TIMER
  WHILE ABS(TIMER - asec) < (idelay/10.): WEND
loop while len(cmd$)>0
end sub

sub opencom(iprt,kp$,baud$,spar$) 'includes parity enabled w/8 data bits
'--- open serial port as unit 3
spa$=spar$
par$=left$(spa$,1)
ndpar=(par$<>"N" and mid$(spa$,3,1)="8") 'QB4 only allows parity N
if ndpar then mid$(spa$,1,1)="N" 'make compatible with QB4
ser$ = kp$+":"+baud$+"," + spa$
if instr(ser$,"RB")=0 then ser$=ser$+",RB2048" 'increase receive buffer
OPEN ser$ FOR RANDOM AS #3
if ndpar then 'change COM port setting
  select case par$
    case "E" 
      icmd=27             'even parity
    case "M" 
      icmd=43             'mark parity (1)
    case "S" 
      icmd=59             'space parity (0)
    case else
      icmd=11             'odd parity
  end select
  if mid$(spa$,5,2)<>"1," then icmd=icmd+4 '1.5 or 2 stop bits
  OUT iprt+3,icmd                           'send to line control register
end if
end sub
