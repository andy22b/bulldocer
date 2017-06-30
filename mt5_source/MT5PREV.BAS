'   program MT5PREV
'   Utility program to quickly plot digitized, interpolated seismograms
'   Output is Golden Software, Postscript, and HPGL
'   SYN3.02 compatible    6/6/89.  R. McCaffrey and G. Abers
'    Adapted for MT4 in 1990 by Peter Zwick
'    modified 9/25/92 by Peter Zwick
'    made max input array size explicit 11/9/93 PZ

DEFINT I-N:OPTION BASE 1
CONST NSTAMAX=1 'dummy variable
rem $INCLUDE: 'M5COMSCR'
rem $INCLUDE: 'm5fndecl'
rem $INCLUDE: 'm5functn'


DEF FNPWAVE = (kphtype=1) OR (kphtype=3) OR (kphtype=4) OR (kphtype=6) OR (kphtype=7)
DEF FNILONG = (kphtype=1) OR (kphtype=2) OR (kphtype=4) OR (kphtype=5)
DEF FNISHORT = (kphtype=3) OR (kphtype=6) OR (kphtype=7)
DEF FNINGDSN = (kphtype>3)
DEF FNIWWSSN = (kphtype<4)
DEF FNKPHLAB$(I) = mid$("WWLPWWLSWWSPGDLPGDLSGDSPGDBP",1+(i-1)*4,4)

' def fnisx(x) = cint(gs2scr*x*scratio)+200  'pixel posn for x in mm
 def fnisy(y) = cint(kscry - gs2scr*y)
  DEF FNHIA(A)=PY + (A*HAMP/DEXT)
  DEF FNHIT(I)=PX + i*dt *tscale + 5*szlet

'-- plot functions
def fntime(i)= tmin + (i-1)*dt
def fnix(t)=kscrx*(t-tmin)/(tmax-tmin)
def fniy(y)=(iyzero-0.8*yamp*y/ymax)
def fnia(a)=(ix-scc*a/ymax)
def fnit(t)=t*xtsc+noff
def fnisx(x) = cint(x*fisx)
'def fnisy(y) = cint(kscry-(y-scrorig)*fisy)
   HAMP=.45*25.4 : SCL=1

  prgnm$="MT5PREV"
  authors$="R. McCaffrey, G. Abers, & P. Zwick"

  CALL welcome
  call scprint(13,3,"Setting up..")

 nesmax=3072
 DIM PXPOS(24),PYPOS(24),ODAT(nesmax)
 DATA .5,9.6, .5,8.7, .5,7.8, .5,6.9, .5,6.0, .5,5.1, .5,4.2, .5,3.3, .5,2.4, .5,1.5
 DATA 4.5,9.6, 4.5,8.7, 4.5,7.8, 4.5,6.9, 4.5,6.0, 4.5,5.1, 4.5,4.2, 4.5,3.3, 4.5,2.4,  4.5,1.5
 nplots=20
 nsuf=97
 FOR I=1 TO nplots
   READ PXPOS(I),PYPOS(I)
   pxpos(i)=pxpos(i)*25.4
   pypos(i)=pypos(i)*25.4
 NEXT

'-- initialize plot parameters
  xmax=9.84:ymax=7.87:xscale=1:yscale=.9907:lcode=0:khpcol=0:scrorig=0.
  kpenhp=kfalse
  letset=20
  szlet=2.  'mm 0.08 in.
  anglet=0
  gscale = 1. '25.4  '--- all commands in inches
  irotate = false
  igraph = 1     '--- 1 for plotcall, 2 for Postscript,etc

'---- convert plot coords (250x200mm) to undistorted screen coordinates
  if (kscry-200*scratio*kscrx/250)<0 then  'if screen too wide
     fisx=kscry/(scratio*200)  '           'pixels/(mm x dimension)
     fisy=kscry/200                        'pixels/(mm y dimension)
  else
     fisx=kscrx/250
     fisy=kscrx*scratio/250
  end if


  gs2scr = kscry/250  '11
'-- initialize variables
  K=0
  inf$="MT5TEST.INV"
  hin$=fnsubext$(inf$)
  nplots=20
  tscale1=50
  plang$="PlotCall2+PlotCall<2PostScriptLaserJet3 Plot8x11  Plot11x17 "
  sizemax=3.9*25.4-5*szlet
  call scpause(13,3)
  ans$="I"

paramenu:
nsuf=97
call setscreen(0)

call scbox(1,13,ktrue,"Menu","Press Q to Quit")
call scprint(3,3,"KEY  OPTION"):call scprint(3,55,"VALUE")
call scprompt(5,4,"I - Input file")
call scprompt(6,4,"O - Output file")
call scprompt(7,4,"N - Number plots per page")
call scprompt(8,4,"T - Time scale of plot (mm/min)")
call scprompt(9,4,"F - plot language")
call scprompt(10,4,"X/"):call scprompt(10,6,"Q - Shell/Quit")
call scprompt(11,4,"P - run plot")
'-- now say variables
call setvarcol
call scprint(5,55,inf$)
call scprint(6,55,hin$)
call scprint(7,55,fnstrnum$(nplots))
call scprint(8,55,fnstrnum$(tscale1))
call scprint(9,55,mid$(plang$,1+igraph*10,10))

210  ik= input1(14,1,"Enter KEY for OPTION ",ans$,kfalse,kfalse,ktrue,"iontfxqp")
if ans$="q" then
    call endcheck:goto 210
elseif ans$="o" then
  do
    h$=hin$
    if inputn(14,1,"Output prefix for Plot (Esc to abort)",H$,ktrue,ktrue,kfalse,"a",40)=27 then exit do
    h$=ucase$(fnsubext$(h$)):hl$=left$(h$,7)+"A.PLT"
    if not keepold(hl$) then hin$=h$:exit do
  loop
elseif ans$="n" then
     if nplots>10 then d$="2" else d$="1"
     if input1((0),(0),"10 plots per page (1) or 20 (2)",d$,kfalse,ktrue,ktrue,"12")<>27 then
       if d$="1" then sizemax=7.9*25.4-5*szlet else sizemax=3.9*25.4-5*szlet
       nplots = 10*val(d$)
     end if
elseif ans$="t" then
     ts$=fnstrnum$(tscale1)
     if inputn((0),(0),"Plotting time scale in mm per minute (5-200)",ts$,ktrue,ktrue,ktrue,"i",3)<>27 then
       if fninlim(val(ts$),5,200) then tscale1=val(ts$)
     end if
elseif ans$="i" then
    d$=inf$
    if inputn((0),(0),"Input file (.INV assumed)",d$,ktrue,ktrue,ktrue,"a",40)<>27 then
      d$=fnaddext$(d$,"inv")
      f$=""
      if isfile(d$,f$,ktrue)=0 then
        inf$=ucase$(f$) 'f$ is full filename
        if hin$="" or hin$="MT5TEST" then hin$=fnsubext$(f$)
      end if
    end if
elseif ans$="f" then
  call scprompt(14,2,"0-PLOTCALL (Golden Software, Version 2 and higher)")
  call scprompt(15,2,"1-PLOTCALL (Golden Software, before Version 2)")
    call scprompt(16,2,"2-PostScript")
    call scprompt(17,2,"3-HP LaserJet3")
    call scprompt(18,2,"4-HP Plotter,8.5"+chr$(34)+"x11"+chr$(34))
    call scprompt(19,2,"5-HP Plotter,11"+chr$(34)+"x17"+chr$(34))
    d$=fnstrnum$(igraph)
    if input1(21,2,"Plot language ",d$,kfalse,ktrue,kfalse,"012345")<>27 then igraph=val(d$)

elseif ans$="x" then call shellcom

elseif ans$="p" then if isfile(inf$,(""),kfalse)=0 then gosub MKPLOT

end if
goto paramenu

'-- start of plot routine
MKPLOT:
 tscale=tscale1/60.   'convert to mm per sec
 hi$=left$(hin$,7)    'leave room for suffix
 OPEN INF$ FOR INPUT AS 1
    open syndir$+"m5letset.raf" as 4 len=2: field #4, 2 as io$

' read header line of data file but don't keep it
    line input #1, a$
    ishdr=len(a$)<50

     call setscreen(kega)
270  view: cls : gosub 700
     view (fnisx(70),0)-(kscrx,kscry)
280 IF EOF(1) THEN 410
    if ishdr then line input #1, a$ else ishdr=ktrue
297   IF LEN(a$)<3 THEN K=K+1 : GOTO 280
       s$=mid$(a$,1,4)
       iphs=val(mid$(a$,5,1))
       instin=fnvalue(6,1) AND 1
       instva= fnvalue(6,1) AND 6 '0=disp,2=vel,4,6=accel
       caltyp$= mid$("dvaa",instva\2+1,1)
       inst=val(mid$(a$,7,1))
       del=val(mid$(a$,8,8))
       azd=val(mid$(a$,16,8))
       xmag=val(mid$(a$,24,8))
       xt0=val(mid$(a$,32,8))
       wt=val(mid$(a$,40,5))
       if0=val(mid$(a$,45,5))
       npt=val(mid$(a$,50,5))
       dext=val(mid$(a$,55,8))
       dt = val(mid$(a$,63,6))
       atime = val(mid$(a$,69,6))
       if instin then
         line input #1, d$
       elseif inst>3 then
         input #1, nz,np,a0
         for i=1 to nz :input #1, zr,zi : next i
         for i=1 to np :input #1, pr,pi : next i
         if xmag=0 then xmag=1
       end if

  if dt=0 and fnilong then  dt=.5
  if dt=0 and fnishort then dt=.1

   DEXT=0
    nplsg = int ( sizemax/tscale /dt )
   NPTS=fnamin(NPT,nesmax) 'limit to ODAT max array size 11/9/93 PZ
   FOR I=1 TO NPTS
     INPUT #1, ODAT(I)
     DEXT=fnamax(DEXT,ABS(ODAT(I)))
   NEXT
   FOR I=NPTS+1 to NPT:INPUT #1, tmp: NEXT 'dump extra data 11/9/93 PZ
   if dext=0 then dext=1

   K=K+1
   if iphs=2 or iphs=4 then PHASE$="S" else phase$="P"
   xcor=1500/xmag : ddext=int(dext*xcor)
      DEXT$ =fnstrnum$(ddext)
      del$ = fnstrnum$(del)
      azd$ = fnstrnum$(azd)
   if NPT>nesmax then print "record ";s$;phase$;" truncated"

'--- new page if changing wave type
'      if k=1 then ph1$=phase$
'      if ph1$<>phase$ then call gsclose(hif$): gosub 700 : k=1 : ph1$=phase$

   locate k+1,1 : PRINT S$;K;NPTS;phase$;xmag;
   X0=PXPOS(K) : Y0=PYPOS(K)
   GOSUB 420 '--- draw seismogram
      if k=1 then gosub 1000 '--- draw tics

 IF K=nplots THEN call gsclose(hif$) : goto 270
    goto 280
410    if kopen7 then call gsclose(hif$)
    CLOSE
    RETURN


420 '--- plot one seismogram  IS at x0,y0
      s$= ucase$(s$)
'--- print station name
     s$=s$+phase$+caltyp$
     call gstext (1,X0-10,Y0+1.4*szlet,s$,anglet,szlet)
     call gstext (1,X0-10,Y0,del$,anglet,szlet)
     call gstext (1,X0-10,Y0-1.4*szlet,azd$,anglet,szlet)
     call gstext (1,X0-10,Y0-2.8*szlet,dext$,anglet,szlet)
     call gstroke

    PX=X0-if0*dt*tscale : PY=Y0
    IF NPTS>NPLSG THEN N=NPLSG ELSE N=NPTS
     'if n<>npts then print "record ";st$;" truncated"
     if px>0 then p=px else p=0

 redim x(n), y(n)
  nin=0
  FOR J=1 TO N
      if fnhit(j)<0 then 590
      nin=nin+1
      x(nin)=fnhit(j) : y(nin)=fnhia(odat(j)) 
590  NEXT

   lsyn=2 : sincr=.1 : kcolor=1
   if nin>0 then call gsline ( x(), y(), nin, lsyn, sincr, kcolor )

 ' tic marks
    call gsmove (fnhit(IF0-1),fnhia(.1*dext),kcolor)
    call gsdraw (fnhit(IF0-1),fnhia(-.1*dext),kcolor)
    call gstroke
  RETURN


700 '--- opens new output file
    hif$=hi$+chr$(nsuf)+".plt"
    print "Output file ";hif$
    k=0
    call gsopen (HIF$)
    call gsinit (hif$)
    xoff=5:yoff=5   'move origin a bit
    call gsorig(xoff,yoff)
    nsuf=nsuf+1
    return


1000   '---- draw time tics every 5 second
   ytop=10
   ybot=6
   for xt=0 to 60 step 5
      i= xt/dt
      itim=i+if0-1
      if fnhit(itim)<7.5*25.4 then
        x=fnhit(itim)
        call gsmove (x,ybot,1)
        call gsdraw (x,ytop,1)
      end if
        call gstroke
   next
      ybot = 8
      call gstext (1,x+15,ybot+2*szlet,PRGNM$+" output on "+date$,anglet,szlet)  'write file name
      call gstext (1,x+2.5,ybot,"60 s",anglet,szlet)  'write file name
      call gstext (1,x+15,ybot,"Input file: "+inf$+"     Page"+str$(nsuf-97),anglet,szlet)  'write file name
   return

rem $INCLUDE: 'm5grsubs'
rem $INCLUDE: 'M5RDSET'

