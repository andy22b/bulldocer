'-- M5RDSET.BAS  reads the \MT5.SET setup info and allows user to customize it
'   By Peter Zwick of RPI  Created Nov, 1990
'-- Latest revision for MT4 on 3/28/1992

SUB wiggle(kbeginx,krow,kfinishx,kdely,kdelx,lne$,kcolr)
'-- draws simulated waveform with summation of two sine waves of diff frequencies
'    for color setting purposes

shared ncolchar
kyrange = 1.5*kdely     '+/- pixels high =1.5 columns high for wiggly line
numcycle = 4
tnumlcycle = .5
kstepsize = kdelx*ncolchar / 100 'pixels per step
LOCATE  krow, kbeginx: PRINT lne$; 'name of this line

kbegx=pos(0)*kdelx
kfinx=kfinishx*kdelx
onecycle = (kfinx - kbegx) / numcycle
longcycle = (kfinx - kbegx) / tnumlcycle
krowp = (krow-.5)*kdely
PSET (kbegx, krowp),kcolr
FOR iw = kbegx TO kfinx STEP kstepsize
   theta1 = 6.28318 * ((iw - kbegx) MOD onecycle) / onecycle
   theta2 = 6.28318 * ((iw - kbegx) MOD longcycle) / longcycle
   LINE -(iw, krowp + kyrange * SIN(theta1) * SIN(theta2)),kcolr
NEXT
END SUB

FUNCTION getdataline$(nfile) 'gets line of non-comment data from a setup file
do 
  inln$=""
  if eof(nfile) then exit do
  line input #nfile, inln$
loop while left$(ltrim$(inln$),1)="'"    'line is all comments - ignore
getdataline$=inln$
END FUNCTION

sub READSET
' read graphics setup file
shared syndir$,drive$,kega,kgrcol(),ktxtcol() AS colortype
shared scratio,nrowchar,ncolchar,prgnm$,kscry,kokcol

'-- set defaults
  kgrcol(0)=0
  do
    if kega=0 then 'set defaults for mono screen
      call setdfltcol
      nrowchar=25:ncolchar=80:kscry=199:kokcol=ktrue
    end if
    noset=ktrue
    while not ifopen("i",9,("\mt5.set"),(0),1)
      d$=displaymsg$(prgnm$+" config file not found. Entering installation menu","",kfalse)
      call getsetup:noset=kfalse
    wend
'-- read in as text lines to allow older versions to work
' comments and lines preceded by single quote "'" will be ignored
    inln$=getdataline$(9)
    syndir$=lftvar$(inln$)
    if instr(":\",right$(syndir$,1))=0 then syndir$=syndir$+"\"
    inln$=getdataline$(9)
    drive$=lftvar$(inln$)
    if instr(":\",right$(drive$,1))=0 then drive$=drive$+"\"
    inln$=getdataline$(9)
    kega=val(lftvar$(inln$))
    call setdfltcol
    kgrcol(kcoldat)=val(lftvar$(inln$))
    kgrcol(kcolsyn)=val(lftvar$(inln$))
    kgrcol(kcolbak)=val(lftvar$(inln$))
    for i=1 to 6
      a$=lftvar$(inln$)
      j=instr(a$,"/") '-- uses / as separator between fore and background colors
      if j>0 then
        ktxtcol(i).ifore=val(left$(a$,j-1))
        ktxtcol(i).iback=val(mid$(a$,j+1))
      else
        exit for
      end if
    next
    inln$=getdataline$(9)
    scratio=val(lftvar$(inln$))
    close 9

'------ KEGA = 1 for CGA B&W graphics display (default)
'------ KEGA = 2 for enhanced graphics display
'------ KEGA = 3 for Hercules graphics monochrome display
'------ KEGA = 4 for VGA high resolution display
'------ KEGA = 5 for MCGA high resolution display

    call getscratio
    tmp$=drive$+"screen.tmp"
    ik=isfile(tmp$,(""),1)
    if ik=0 then kill(tmp$)
    if noset and prgnm$="MT5" then
      IF isfile((syndir$+prgnm$+".EXE"),(""),1)<>0 then
        d$=displaymsg$(syndir$+" is invalid path or drive for program files","",kfalse)
      ELSE
        SELECT CASE ik
          CASE 3,&HF,&H15
            d$=displayerr$(ik,("temp files path "+drive$),kfalse)
          CASE ELSE
            exit do
        END SELECT
      END IF
    else
      exit do
    end if
    call getsetup
  loop
end sub
'========================================================================
sub GETSETUP
shared syndir$,drive$,kega,kgrcol(),ktxtcol() AS colortype
shared kscrx,kscry,nrowchar,ncolchar,scratio,prgnm$

dim kcolra(3)            'color for 3 boxes
dim lntyp$(3)          'Labels for 3 line types
kgrcol(0)=0
w$="P"
GETSET:
cls
klist$="ptsxq"
if fninlim(kega,1,5) then bottom$="Press Esc to abort" else bottom$="Must set SCREEN TYPE to continue":w$="S"
call scbox(1,14,kfalse,"\MT5.SET - MENU FOR FILE DIRECTORIES AND SCREEN PARAMETERS",bottom$)
call scprint(3,4,"KEY  OPTION")         : call scprint(3,55,"CURRENT VALUE")
call scprompt(4,5,"P - Program files directory")
call scprompt(5,5,"T - Temporary files directory")
call scprompt(6,5,"S - Screen type (1-5)")
call scprint(7,5,"A - Aspect ratio for graphics")
call scprint(8,5,"C - Colors for text")
call scprint(9,5,"G - colors for Graphics")
call scprompt(10,5,"X - Shell")
call scprompt(11,5,"Q - Quit "+prgnm$)
call scprint(12,5,"D - Done with changes. Save \MT5.SET and return")
if fninlim(kega,1,5) then
  klist$=klist$+"ac"
  call setchccol
  call scprint(7,5,"A")
  call scprint(8,5,"C")
  if kega>1 and kega<>3 and kega<>5 then call scprint(9,5,"G"):klist$=klist$+"g"
  klist$=klist$+"d"
  call scprint(12,5,"D")
  call settxtcol
end if

'-- now write variables
call setvarcol
call scprint(4,55,syndir$)
call scprint(5,55,drive$)
curkega$=fnstrnum$(kega)+" - "+mid$("NONECGA EGA HercVGA MCGA",kega*4+1,4)
call scprint(6,55,curkega$)
call scprint(7,55,str$(scratio))

ik=input1(15,1,"Enter KEY for option ",w$,kfalse,fninlim(kega,1,5),ktrue,klist$)
if ik=3 then end  'Ctrl-Brk
if ik=27 then     'Esc
  if igetyesno(15,1,"Abandon any changes?","n",ktrue) then exit sub

elseif w$="p" then
  d$=syndir$
  do
   if inputn(15,1,"Enter new program files directory ",d$,ktrue,ktrue,kfalse,"a",40)=27 then exit do
   nddir=len(d$)>0
   if nddir then
      d$=lcase$(alltrim$(d$))
      if instr(":\",right$(d$,1))=0 then d$=d$+"\"
      if isfile((d$+prgnm$+".exe"),(""),1)=0 then
         syndir$=d$
         nddir=kfalse
         w$="T"
      else
         call cline(16,1)
         call scprint(16,1,prgnm$+" program files not found in "+d$)
      end if
   end if
  loop while nddir

elseif w$="t" then
  d$=drive$
  do
   if inputn(15,1,"Enter new temporary files directory ",d$,ktrue,ktrue,kfalse,"a",40)=27 then exit do
   nddir=len(d$)>0
   if nddir then
      d$=lcase$(alltrim$(d$))
      if instr(":\",right$(d$,1))=0 then d$=d$+"\"
      ierrnum=2
      ierrnum= isfile((d$+"nul"),(""),1)
      if ierrnum=0 then
            drive$=d$
            nddir=kfalse
            if kega=0 then w$="S" else w$="A"
      else
         a$=displayerr$(ierrnum,d$,kfalse)
      end if
   end if
  loop while nddir

elseif w$="s" then
  call scprint(15,2,"Screen types")
  call scprompt(16,5,"1 - CGA B&W graphics display (640x200)")
  call scprompt(17,5,"2 - enhanced color graphics display (640x350)")
  call scprompt(18,5,"3 - Hercules graphics monochrome display (720x350)")
  call scprompt(19,5,"4 - VGA high resolution color display (640x480)")
  call scprompt(20,5,"5 - MCGA high resolution B&W graphics display (640x480)")

  i=kega
  if kega=0 then
   w$="P"
   kega=1':kgrcol(kcoldat)=1:kgrcol(kcolsyn)=1:kgrcol(kcolbak)=1
  else
   w$="A"
  end if
  k$=fnstrnum$(kega)
  call scprint(22,1,"Be sure to enter the correct screen type or the screen graphics will bomb")
  if input1(23,1,"Enter new screen type (1-5)",k$,kfalse,ktrue,kfalse,"12345")<>27 then
    kega=val(k$)
    if kega=3 then
      call scprint(24,1,"NOTE: To display graphics on the Hercules monitor, you must run MSHERC prior")
      call scprint(25,1," to running this program. - ")
      if igetyesno((0),(0),"QUIT "+prgnm$+" TO RUN MSHERC?","",kfalse) then end
    end if
    call setscreen(kega)
    call setscreen(0)
    if i<>kega then
      scratio=0: call getscratio
      j=instr("124",k$)
      k=instr("124",fnstrnum$(i))
      if not (j>0 and k>0) then call setdfltcol
    end if
  end if

elseif w$="a" then
 call setscreen(kega)
 call getscratio
 s$=" " : idxc=kscrx*.4 : idy=kscry*(.8-4./nrowchar)/2 : ids=kscry*.1
 scrathld=scratio
 DO
   CLS
   call scline2 (idxc-idy/scratio, ids, idxc+idy/scratio, ids + 2*idy, kgrcol(kcolsyn),ktrue)
   call scline2 (idxc-idy/scratio, ids, idxc+idy/scratio, ids + 2*idy, kgrcol(kcolbak),kfalse)
   call scline2 (idxc+idy/scratio, ids, idxc-idy/scratio, ids + 2*idy , kgrcol(kcolbak),kfalse)
   call sccircle(idxc, ids + idy, idy/scratio, kgrcol(kcoldat))
   call scprint(nrowchar-3,1, "Use this display to set screen aspect ratio (currently ")
   call scprintvfu(nrowchar-3,(0),fnfform$(4,2)+")",scratio)
   call scprint(nrowchar-2,1,"Enter "): call scprompt((0),(0),"W to make the box/circle Wider and ")
     call scprompt((0),(0),"N to make them Narrower")
   call scprint(nrowchar-1,1,"Enter "): call scprompt((0),(0),"Q when the box is square and the circle is round")
   while inkey$<>"":wend  'clear buffer
   s$ = INPUT$(1) : s$=lcase$(s$)
   IF s$ = "n" THEN scratio = scratio + .01
   IF s$ = "w" THEN scratio = scratio - .01
 LOOP UNTIL s$ = "q" or s$=chr$(27) 'quit or escape
 if s$=chr$(27) then scratio=scrathld else w$="C"
 call setscreen(0)

elseif w$="c" then
 d$="N"
 DO
  cls
  if kega < 2 or kega = 3 or kega = 5 then
   dfltcol$=" 7 0 1 0 7 015 0 0 7 0 7"
  else
   dfltcol$=" 3 014 0 1 310 0 0 4 0 3"
  end if
  call scbox(1,14,kfalse,"\MT5.SET - TEXT COLOR SELECTION MENU","Press Esc to exit")
  call scprint(3,4,"KEY  OPTION")         : call scprint(3,48,"CURRENT COLOR SETTING")
  call scprompt(4,5,"N - Normal text")
  color 15,0:print " <--(color N)";
  call scprompt(5,5,"V - Value or variable")
  color 15,0:call scprint(5,59,"<--(color V)")
  call scprint(6,2,"<--")
  call scprompt(6,5,"B - Border")
  call scprompt(7,5,"K"):color 15,0:print "<--";
  call settxtcol:print "Key choices available and prompts"
  call scprompt(8,5,"M - Messages (e.g, ")
  call setmsgcol:print "Error on entry...";:call settxtcol:print ")";
  call scprompt(9,5,"H - Highlighted area for your response")
  for i=1 to 6
   call scprintv(3+i,55,fnstrnum$(ktxtcol(i).ifore)+"/"+fnstrnum$(ktxtcol(i).iback))
  next
  color 7,0
  call scprint(11,2,"Color choices:")
  call scprint(12,17,"foreground->")
  call scprint(13,17,"background->"):icol=pos(0)
  for i=0 to 15
   color 7,0
   call scprint(11,icol+(i*3),str$(i))
   color i,0:call scprint(12,icol+(i*3),str$(i))
  next
  for i=0 to 7
   color 7,0
   color 0,i:call scprint(13,icol+(i*3),str$(i))
  next

  klst$="nvbkmh"
  color 15,0
  call scprint(15,26,"<--(color H)")
  if input1(15,1,"Enter KEY for option ",d$,kfalse,ktrue,ktrue,klst$)=27 then exit do
  i=instr(klst$,d$)
  if i>0 then
    c$=fnstrnum$(ktxtcol(i).ifore)+"/"+fnstrnum$(ktxtcol(i).iback)
    d$=ltrim$(mid$(dfltcol$,i*4-3,2))+"/"+ltrim$(mid$(dfltcol$,i*4-1,2))
    print "Color ranges: foreground(0-15), background(0-7). Default is "+d$
    do
     isok=inputn((0),(0),"Enter colors as foreground/background (or D for default)",c$,ktrue,ktrue,ktrue,"a",4)=27
     if not isok then
      if lcase$(c$)="d" then c$=d$
      isl=instr(c$,"/"):if isl then mid$(c$,isl)=" "
      k1$=lftvar$(c$)
      isok=istype(k1$,"i")
      if isok then k2$=lftvar$(c$):isok=istype(k2$,"i")
      if isok then isok = fninlim(val(k1$),0,15) and fninlim(val(k2$),0,7)
      if isok then
        ktxtcol(i).ifore=val(k1$)
        ktxtcol(i).iback=val(k2$)
      end if
     end if
    loop while not isok
    d$=ucase$(mid$(klst$+" ",i+1,1))
  end if
 LOOP
 w$=ucase$(mid$(klist$,8,1))

elseif w$="g" then
 call setscreen(kega)
 call getscratio
 s$=""  : idy=kscry/2
   DO
      CLS
      kdely=(kscry+1)/nrowchar      'pixels per line
      kdelx=(kscrx+1)/ncolchar      'pixels per column
      kulx=ncolchar/4:kuly=2      'upper left corner of line type box display (column & row)
      kystep=4  ' step values for box corner
      lntyp$(1)="DATA"       :kcolra(1)=kgrcol(kcoldat)
      lntyp$(2)="SYNTHETICS" :kcolra(2)=kgrcol(kcolsyn)
      lntyp$(3)="MISC"       :kcolra(3)=kgrcol(kcolbak)
      for i=1 to 3 'boxes 1 to 3
         if i<3 then
            kulyn=kuly+i*kystep
            call wiggle(kulx,kulyn,kulx+ncolchar\2,kdely,kdelx,lntyp$(i),kcolra(i))
         else
            call scline2 ((kulx-2)*kdelx, (kuly-.5)*kdely, (kulx+1+ncolchar\2)*kdelx, (kuly+i*kystep-.5)*kdely,  kcolra(i),ktrue)
            LOCATE  kuly + idy/kdely,kulx: PRINT lntyp$(i);
         end if
      next
      locate nrowchar-6,1:print "Color number choices: ";
      kxstep=3
      kuly=csrlin: kulyp=(kuly-1)*kdely  ' char and pixel rows to start on
      kulx=pos(0)             ' column to start in
      for i=1 to 15
         kulxn=kulx+i*kxstep
         line (kulxn*kdelx,kulyp)-STEP(kxstep*kdelx,-kdely),i,BF  'draw a filled box
         LOCATE kuly,kulxn+1:print i;
      next
      locate nrowchar-4,1
      d$=fnstrnum$(kgrcol(kcoldat))+","+fnstrnum$(kgrcol(kcolsyn))+","+fnstrnum$(kgrcol(kcolbak))
      print "Colors used for graphics of ";lntyp$(1);", ";lntyp$(2);", and ";lntyp$(3);_
        " are currently: "+d$
      kchk=igetyesno((0),(0),"Do you want to change them?","Y",ktrue)
      if kchk then
         if inputn(nrowchar-2,1,"Enter new numbers for colors",d$,ktrue,ktrue,ktrue,"a",8)<>27 then
            k$=lftvar$(d$):if istype(k$,"i") then k=val(k$): if fninlim(k,1,15) then kgrcol(kcoldat)=k
            k$=lftvar$(d$):if istype(k$,"i") then k=val(k$): if fninlim(k,1,15) then kgrcol(kcolsyn)=k
            k$=lftvar$(d$):if istype(k$,"i") then k=val(k$): if fninlim(k,1,15) then kgrcol(kcolbak)=k
         end if
      end if
   loop while kchk
 call setscreen(0)
 w$="D"

elseif w$="x" then
 call shellcom

elseif w$="q" then
 call endcheck

elseif w$="d" then
 if igetyesno(15,1,"Update \MT5.SET and use changes?","Y",kfalse) then
  open "o",9,"\mt5.set"
  print #9, syndir$
  print #9, drive$
  c$=""
for i=1 to 6
   c$=c$+str$(ktxtcol(i).ifore)+"/"+fnstrnum$(ktxtcol(i).iback)
next
print #9, kega;kgrcol(kcoldat);kgrcol(kcolsyn);kgrcol(kcolbak);c$
  write #9, scratio
  close 9
  exit sub
 end if
end if
goto GETSET
end sub
'========================================================================
sub getscratio
shared nrowchar,ncolchar,kscrx,kscry,scratio,kega
'nrowchar=25 : ncolchar=80
' Note that graphics x,y dimensions are from (0) to (dimension-1)
scasp=4./3
if not fninlim(kega,2,5) then
   kega=1 : kscrx=639 : kscry=199: scasp=23/15.3
elseif kega=2 then
   kscrx=639 : kscry=349
elseif kega=3 then
   kscrx=719 : kscry=347
elseif kega>=4 then
   kscrx=639 : kscry=479 : 'nrowchar=30 only change in setscreen()
end if
if scratio=0 then scratio= scasp*(kscry+1)/(kscrx+1)

end sub
'=========================================================================

'--- welcome message
sub welcome
shared prgnm$,version$,authors$,ipmi

ipmi=len(ENVIRON$("IPMI"))>0 ' debug info is printed to drive$+"DEBUG.LOG"

call setscreen(0)
print prgnm$+" - Reading config file...";
call READSET
call scbox(1,10,0,"IASPEI Software Library Volume 3","")
call sccenter("Welcome to "+prgnm$,3)
if len(version$)>0 then call sccenter("<Version "+version$+">",4)
call sccenter("by",5)                                       
call sccenter(authors$,7)
end sub
'=========================================================================
sub setdfltcol
shared kega,kgrcol(),ktxtcol() AS colortype
  kgrcol(0)=0
'   2 char each for fore & background of txt,var,box,chc,msg,hil
  if (kega = 0 or kega = 3 or kega = 5) then
    dfltcol$=" 7 0 1 0 7 015 0 0 7 0 7"
  else
    dfltcol$=" 3 014 0 1 310 0 0 4 0 3"
  end if
  for j=1 to 6  'read defaults
    ktxtcol(j).ifore=val(mid$(dfltcol$,j*4-3,2))
    ktxtcol(j).iback=val(mid$(dfltcol$,j*4-1,2))
  next
if kega < 2 or kega = 3 or kega = 5 then  'can't display colors
  kgrcol(kcoldat)=1 : kgrcol(kcolsyn)=1 : kgrcol(kcolbak)=1
else
  kgrcol(kcoldat)=4 : kgrcol(kcolsyn)=5 : kgrcol(kcolbak)=6
end if
end sub
