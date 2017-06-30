' main MT5LIST
' make list of WWSSN and GDSN stations and their distances from event
'  this is to aid in deciding which seismograms to use
' input format for events is in EVENT HEADER format given in manual

DEFINT I-N: OPTION BASE 1
CONST NSTAMAX=1  'dummy variable not used

rem $INCLUDE: 'M5COMSCR'
rem $INCLUDE: 'm5fndecl'
rem $INCLUDE: 'm5functn'

'--- welcome
prgnm$="MT5LIST"
authors$="R. McCaffrey, G. Abers, & P. Zwick"
call welcome

'-- maximum number of stations in file is 250
maxwwsta=250
dim wws$(maxwwsta), adlat(maxwwsta), adlon(maxwwsta),isnum(maxwwsta)
dim network(maxwwsta), azimdel(maxwwsta)
'dim timeon$(maxwwsta), timeoff$(maxwwsta)

dim iah(2),iam(2),asec(2), isflag(10)

'----------- travel time tables ---------
'-- set up for 2 phases, 120 deltas, 15 depths
dim ttimtab(120,15,2), dep(15,2), vsurf(2)
dim del0(2), ddelt(2), ndel(2), ndep(2)

'-- read in stations
   gosub READSTA

'-- read travel time table
numttype=2
dim ttfilenm$(numttype),ttname$(numttype)
data "m5j-b.dat","m5iasp91.dat"
data "Jeffreys-Bullen (used for EQ location by ISC)","IASPEI 91"
for i=1 to 2 : read ttfilenm$(i) : next
for i=1 to 2 : read ttname$(i) : next
  ittype=1      'default travel time table
  ttfile$=syndir$+ttfilenm$(ittype)
  call readtime (ttfile$)

'-- set up sort order options
numisort=5
dim sortname$(numisort)
data "unsorted","by number","by name","by azimuth","by delta"
for i=1 to numisort : read sortname$(i) : next
isort=2         'default sort order

  call scpause(18,3)

   islist$="01" : isflag(1)=ktrue: isflag(2)=ktrue
   neq=0
   fin$="KEYBOARD" : va2=6.0 : dmin=30: dmax=90
   neqs=1: fout$="STALIST.OUT"

MENU:
a$=""
call setscreen(0)
call scbox(1,17,ktrue,"MAIN MENU","Press Q to quit")
call scprint(2,2,"KEY  OPTION                                     VALUE")
call scprompt(3,3,"I - Input file")
call scprompt(4,3,"O - Output file")
call scprompt(6,3,"P - P wave velocity at source (km/sec)")
call scprompt(7,3,"M - Minimum and maximum delta for list")
call scprompt(9,3,"N - Networks to include in list")
call scprint(10,7,"(0)WWSSN (1)DWWSSN, (2)SRO/ASRO, (3-9)others (e.g. 01 for 0 and 1)")
call scprompt(11,3,"S - Station sort order")
call scprompt(12,3,"T - Travel time table")
call scprompt(14,2,"X/"): call scprompt(14,4,"Q - Shell/Quit")
call scprompt(15,3,"B - Begin")

'-- now write variables
call setvarcol
call scprint(3,50,fin$)
call scprint(4,50,fout$)
call sclocate(6,49):print va2;
call sclocate(7,49):print dmin;dmax;
call scprint(9,50,islist$)
call scprint(11,50,sortname$(isort))
call scprint(12,78-fnamax(28,len(ttname$(ittype))),ttname$(ittype))

PRMPT:
ik= input1(18,1,"Enter KEY for option ",a$,kfalse,kfalse,ktrue,"piomnxbqts")

if a$="q" then call endcheck:goto PRMPT

if a$="p" then
    b$=fnstrnum$(va2)
    if inputn(18,2,"Source P-velocity (km/sec)",b$,ktrue,ktrue,ktrue,"f",4)<>27 then
      v=val(b$): if v>0 then va2=v
    end if

elseif a$="b" then
     goto BEGIN

elseif a$="x" then
     call shellcom

elseif a$="i" then
    f$=fin$
    if inputn(18,2,"Input file (K for KEYBOARD)",f$,ktrue,ktrue,ktrue,"a",45)<>27 then
      if len(f$)>0 then
        f$=ucase$(f$)
        if f$="K" or f$="KEYBOARD" then
          fin$="KEYBOARD"
        else
          if isfile(f$, ftmp$, -1)=0 then fin$=ucase$(ftmp$)
        end if
      end if
    end if

elseif a$="m" then
    b$=fnstrnum$(dmin)
    if inputn(18,2,"Minimum delta for list (0-180)",b$,_
         ktrue,ktrue,kfalse,"f",6)<>27 then
      dmin=fnxlimit(val(b$),0,180)
      dmax= fnxlimit(dmax,dmin,180):b$=fnstrnum$(dmax)
      if inputn(19,2,"Maximum delta for list ("+fnstrnum$(dmin)+"-180)",b$,_
          ktrue,ktrue,ktrue,"f",6)<>27 then
        dmax= fnxlimit(val(b$),dmin,180)
      end if
    end if

elseif a$="n" then
    print "Enter networks to include in list
    if inputn((0),2,"(0)WWSSN, (1)DWWSSN, (2)SRO/ASRO, (3-9)others (e.g. 01 for 0 and 1)",_
       islist$,ktrue,ktrue,ktrue,"i",10)=27 goto menu
    for i=0 to 9
      l$=fnstrnum$(i)
      if instr(islist$,l$) then isflag(i+1)=ktrue else isflag(i+1)=kfalse
    next
    if isflag(1) then isflag(2)=ktrue

elseif a$="o" then
    f$=fout$
    if inputn(18,2,"Output file (CON for screen)",f$,ktrue,ktrue,ktrue,"a",45)<>27 then
      f$=LTRIM$(ucase$(f$))
      if (instr(f$,"*")+instr(f$,"?"))=0 then
        if len(f$)=0 then goto prmpt
        if f$<>"CON" then if keepold(f$) goto prmpt
        fout$=f$
      else if isfile(f$,ftmp$,-1)<>0 goto menu
        fout$=ftmp$
      end if
    end if

elseif a$="t" then
   call scprint(18,2,"KEY  TRAVEL TIME TABLE SELECTION")
   b$=fnstrnum$(ittype)
   for i=1 to numttype
     call scprompt(18+i,3,fnstrnum$(i)+" - "+ttname$(i))
   next
   if input1(18+i,1,"Enter KEY for selection ",b$,kfalse,ktrue,ktrue,numlist$(1,numttype))=27 then goto MENU
   ia=val(b$)
   if ia<>ittype then
     ittype=ia
     ttfile$=syndir$+ttfilenm$(ittype)
     for i=18 to 19 + numttype : call cline(i,1) : next
     call readtime (ttfile$)
   end if
elseif a$="s" then
   call scprint(18,2,"KEY  SORT ORDER")
   b$=fnstrnum$(isort)
   for i=1 to numisort
     call scprompt(18+i,3,fnstrnum$(i)+" - "+sortname$(i))
   next
   if input1(18+i,1,"Enter KEY for sort order",b$,kfalse,ktrue,ktrue,numlist$(1,numisort))=27 then goto MENU
   ia=val(b$)
   if ia<>isort then
     isort=ia
     for i=18 to 19 + numisort : call cline(i,1) : next
     if isort=1 then
       gosub readsta
     elseif isort<4 then
       gosub sortsta
     end if
   end if
end if
goto MENU

'========================================================================


BEGIN:
     ikey=(fin$="KEYBOARD")
     if (not ikey) then if not ifopen("i",6,fin$,0,ktrue) then goto MENU

     open "o",2,fout$
     print : print " Press Q or Esc to interrupt "

NEXTEQ:

    if fninquit(0) or fninesc then goto 1100
    if ikey then
      kheader=kfalse
    else
      kheader=ktrue
      if eof(6) then 1100
      line input #6, a$
    end if

      gosub GETHYPO
      if iend then goto 1100
      if ierrin then print "Not a valid event " : print a$ : goto NEXTEQ

      ipg=1 : nout=0 : neq=neq+1
      print : print "Processing ";dat$
478  npp=53

560 gosub 2000

580 for ks=1 to nwsta

       if fninquit(0) or fninesc then goto 1100
       nw=network(ks)
        if not isflag(nw+1) then 1080
'       if dat$<timeon$(ks) then 1080
'       if dat$>timeoff$(ks) then 1080

600    stalat=adlat(ks) : stalon=adlon(ks)
      call m5delaz(stalat,stalon,eplat,eplon,del,az)
       if del <dmin or del >dmax then 1080
680   call m5delaz(eplat,eplon,stalat,stalon,del,baz)

'--- get theoretical arrival time, and derivatives dT/dh, dT/dz
     toa=15
760   for iph=1 to 2
       iah(iph)=0: iam(iph)=0 : asec(iph)=0
       call travtime (iph, del, dep, ttime, dtdh, dtdz)
       if ttime=0 then 805
       call ellipcor (iph, del, eplat, stalat, ecor)
       ttime=ttime+ecor
770   ttime=ttime + origin
780   ah=int(ttime/3600) : iah(iph)=int(ah)
790   ttime=ttime-ah*3600 : am=int(ttime/60) :iam(iph)=int(am)
800   asec(iph)=ttime-am*60
      if iph=1 then toa=fndeg(fnasin(dtdh*va2))
805  next

950   nout=nout+1 : if nout =npp then gosub 2000
980   PRINT #2, USING fnaform$(4); ucase$(WWS$(KS));
1000  PRINT #2, USING fniform$(4);ISNUM(KS);
1010  PRINT #2, USING fnfform$(10,2);DEL; AZ; BAZ; TOA;
1020  print #2," ";
      for iph=1 to 2
      print #2, "  ";
1025  PRINT #2, USING fniform$(2)+":"; IAH(iph); IAM(iph);
1030  PRINT #2, USING fnfform$(4,1); ASEC(iph);
      next : print #2, ""

1080 NEXT

1090 goto NEXTEQ

1100 if (fout$<>"KEYBOARD") then call schit("Saving "+fout$+" -")
     close
     goto menu

2000 '----- subroutine to print header
2010 if neq=1 and ipg=1 then 2020
2015  print #2, chr$(12)  'form feed
2020 print #2,"Event";tab(10);
2022 print #2,"  Latitude Longitude     Depth  Origin Time";
2023 print #2, tab(75);"Page"
2024 print #2, dat$;tab(10);
2028 print #2, using fnfform$(10,3);eplat;eplon;dep;
     print #2, "   ";
2030 print #2, using fniform$(2)+":";ohr;omin;
2040 print #2, using fnfform$(4,1); osec;
2042 print #2, tab(76);ipg
2045 print #2, " "
2050 PRINT #2,"Stat.  N     Delta   Azimuth    Backaz  Take-off     P-time     S-time"
2055 PRINT #2, STRING$(79,45)
2060 nout=0: ipg=ipg+1
2100 return

READSTA: '---- read tables and stations
      call cline(18,1)
      call scprint(18,3 ,"Reading station list ")

10100 iff$=syndir$+"m5statio.dat"
10220 if not ifopen("i",1,iff$,0,ktrue) then end
      print iff$;"...";
      i=0

10240 line input #1, a$
        if fnvalue$(1,1)="'" goto 10240
        if fnvalue$(1,3)="zzz" goto 10360
        if eof(1) then goto 10360
        i=i+1
      wws$(i)=fnvalue$(1,4)  : adlat(i)=fnvalue(5,7)
      adlon(i)=fnvalue(12,9) : isnum(i)=fnvalue(21,4)
      network(i)=fnvalue(27,1)
'     timeon$(i)=fnvalue$(30,6) : timeoff$(i)=fnvalue$(37,6)
      if i=maxwwsta then 10360
      goto 10240

10360 nwsta=i
      close 1

SORTSTA: 'sort stations
'---  re-order stations by number,name,azimuth,or delta
    if isort=1 then return
    print " Sorting..."
    for i=1 to nwsta
      min=i
      if isort=2 then
        minn=isnum(i)
      elseif isort=3 then
        minn$=wws$(i)
      else
        minn=azimdel(i)
      end if
      for j=i+1 to nwsta
        if isort=2 then
          if isnum(j)<minn then
            minn=isnum(j): min=j
          end if
        elseif isort=3 then
          if wws$(j)<minn$ then
            minn$=wws$(j): min=j
          end if
        else
          if azimdel(j)<minn then
            minn=azimdel(j): min=j
          end if
        end if
      next
      if min<>i then
        swap wws$(i),wws$(min)
        swap adlat(i),adlat(min)
        swap adlon(i),adlon(min)
        swap isnum(i),isnum(min)
        swap network(i),network(min)
        if isort>3 then swap azimdel(i),azimdel(min)
'       swap timeon$(i),timeon$(min)
'       swap timeoff$(i),timeoff$(min)
      end if
    next
  return
'======================================================================
GETHYPO:
'-- read in hypocenter data
  iend=0 : ierrin=0

  if kheader then
    if len(a$)<20 then ierrin=ktrue : return
    i=fnvalue(1,2)
      if not fninlim(i,0,99) then ierrin=ktrue : return
    i=fnvalue(3,2)
      if not fninlim(i,1,12) then ierrin=ktrue : return
    i=fnvalue(5,2)
      if not fninlim(i,1,31) then ierrin=ktrue : return
    i=fnvalue(7,2)
      if not fninlim(i,0,24) then ierrin=ktrue : return
    i=fnvalue(9,2)
      if not fninlim(i,0,60) then ierrin=ktrue : return
    dat$=fnvalue$(1,6)
    ohr=fnvalue(7,2)
    omin=fnvalue(9,2)
    osec=fnvalue(11,3)/10
    eplat=fnvalue(14,5)/100
    eplon=fnvalue(19,6)/100
    dep=fnvalue(25,3)
    hend$=fnvalue$(28,len(a$)-27)
    xmag=fnvalue(28,2)/10
 else
    cls
    iend=ktrue

    do
      b$=dat$
      if inputn(1,2,"Date of event (YYMMDD, 000101-991231) (Esc for MENU)",_
        b$,ktrue,ktrue,kfalse,"i",6)=27 then return
    loop until fninlim(val(mid$(b$,3,2)),1,12) and_
       fninlim(val(mid$(b$,5,2)),1,31)
    dat$=b$
    do
      b$=fnstrnum$(ohr)
      if inputn(2,2,"Hour (0-24)",b$,ktrue,ktrue,kfalse,"i",2)=27 then return
    loop until fninlim(val(b$),0,24)
    ohr=val(b$)
    do
      b$=fnstrnum$(omin)
      if inputn(3,2,"Minute (0-60)",b$,ktrue,ktrue,kfalse,"i",2)=27 then return
    loop until fninlim(val(b$),0,60)
    omin=val(b$)
    do
      b$=fnstrnum$(osec)
      if inputn(4,2,"Second (0.00-60.00)",b$,ktrue,ktrue,kfalse,"f",5)=27 then return
    loop until fninlim(val(b$),0,60)
    osec=val(b$)
    do
      b$=fnstrnum$(eplat)
      if inputn(5,2,"Latitude (+North,-90.00 to 90.00)",b$,ktrue,ktrue,kfalse,"f",6)=27 then return
    loop until fninlim(val(b$),-90,90)
    eplat=val(b$)
    do
      b$=fnstrnum$(eplon)
      if inputn(6,2,"Longitude (+East,-180.00 to 180.00)",b$,ktrue,ktrue,kfalse,"f",7)=27 then return
    loop until fninlim(val(b$),-180,180)
    eplon=val(b$)
    do
      b$=fnstrnum$(dep)
      if inputn(7,2,"Depth (km, 0.0-800.0)",b$,ktrue,ktrue,kfalse,"f",5)=27 then return
    loop until fninlim(val(b$),0,800)
    dep=val(b$)
    print
    iend=kfalse

 end if
 if isort>3 then  'need to sort stations by azimuth or delta
   for i=1 to nwsta
     call m5delaz(adlat(i),adlon(i),eplat,eplon,del,az)
     if isort=4 then azimdel(i)=az else azimdel(i)=del
   next
   gosub sortsta
 end if

'-- origin time in seconds
    origin=ohr*3600 + omin*60 + osec
return
'======================================================================
rem $INCLUDE: 'm5delaz'
rem $INCLUDE: 'm5ttimes'
rem $INCLUDE: 'm5rdset'
'======================================================================
