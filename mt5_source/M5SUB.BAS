'-- M5SUB.BAS
' Needs to link with QB.QLB which includes INTERRUPT code
' to get this, start QuickBasic with command "QB /l"

DEFINT I-N: OPTION BASE 1
CONST NSTAMAX = 50 'limit on seismograms

rem $INCLUDE: 'QB.BI'

rem $INCLUDE: 'M5COMSCR'
rem $include: 'm5fndecl'
rem $include: 'm5functn'

sub shellcom
'-- run DOS commands with shell escape
 call scprintc((0),(0),"Type EXIT at the DOS prompt to re-enter program")
 shell
end sub
'=======================================================================
sub inlet (kwait, a$, klen, kval, kascii)
' read a letter, convert to lower case, and get its value, and ASCII code
 kval=0 : kascii=0
 if kwait then input a$ else a$=input$(1)
  a$=lcase$(a$) :  klen=len(a$)
  if klen=0 then exit sub
    a$=left$(a$,1) : klen=1
    kascii=asc(a$) : kval=val(a$)
end sub
'=======================================================================
'-- screen subroutines

sub setscreen(k)
'-- Set screen attribute
  shared nrowchar, ncolchar, kscreen, kokcol
  cls

  nrowchar=25:ncolchar=80
  if k=1 then
    kscreen = 2          'medium resolution, no color
  elseif k=2 then
    kscreen = 9          'EGA card, high res with color
  elseif k=3 then
    kscreen = 3          'Hercules card, high res monochrome
  elseif k=4 then
    kscreen = 12         'VGA high resolution graphics with color
  elseif k=5 then
    kscreen = 11         'MCGA high resolution graphics, no color
  else                   'text mode
    kscreen = 0
  end if
  screen kscreen : if kscreen>10 then nrowchar=30: width 80,30
  kokcol = kscreen<>2 and kscreen<>3 and kscreen<>11  ' ok to use color
end sub

sub setcolor (kcolor)
'-- set color for text if screen mode allows it
'-- ktxtcol is array for normal,variable,box,choice,msg,hlite colors
'    (foreground & background for each)
shared kokcol,ktxtcol() AS colortype,kscreen
  if kokcol then  'if ok to issue color command
    if kscreen<12 then
      color ktxtcol(kcolor).ifore,ktxtcol(kcolor).iback
    else
      if ktxtcol(kcolor).ifore=0 then color 15 else color ktxtcol(kcolor).ifore
    end if
  end if
end sub

sub settxtcol
  call setcolor(1)
end sub

sub setvarcol
  call setcolor(2)
end sub

sub setboxcol
  call setcolor(3)
end sub

sub setchccol
  call setcolor(4)
end sub

sub setmsgcol
  call setcolor(5)
end sub

sub sethilcol
  call setcolor(6)
end sub

sub scline1 (ix, iy, kcolor, kbox)
'-- draw line to X,Y; if KBOX make a box
  if kbox then
   line -(ix,iy),kcolor,b
  else
   line -(ix,iy),kcolor
  end if
end sub

sub scline2 (ix1, iy1, ix2, iy2, kcolor, kbox)
'-- draw line from X1,Y1 to X2,Y2; if KBOX make a box
  if kbox then
   line (ix1,iy1)-(ix2,iy2),kcolor,b
  else
   line (ix1,iy1)-(ix2,iy2),kcolor
  end if
end sub

sub scdash2 (ix1, iy1, ix2, iy2, kcolor, kbox)
'-- draw dashed line from X1,Y1 to X2,Y2; if KBOX make a box
  if kbox then
   line (ix1,iy1)-(ix2,iy2),kcolor,b,&Hf0f0
  else
   line (ix1,iy1)-(ix2,iy2),kcolor,,&Hf0f0
  end if
end sub

sub sccircle (ix, iy, irad, kcolor)
'-- draw circle at IX,IY with radius IRAD
shared scratio
  if irad<=0 then iradn=1 else iradn=irad
  if scratio>1 then iradn=iradn*scratio
  circle(ix,iy),iradn,kcolor,,,scratio
end sub

sub scclear
'-- clear the screen
  cls
end sub

sub cline (irow, icol)
'-- clear a text line starting at IROW,ICOL
 call sclocate(irow,icol)
 print string$(81-icol,32);
 locate irow,icol
end sub

sub cnchar (irow, icol, n)
'-- clear n characters of a text line starting at IROW,ICOL
 call settxtcol
 call scprint(irow,icol,string$(n,32))
 locate irow,icol
end sub

sub scpoint(ix,iy,kcolor)
'-- put a point at X,Y
  pset (ix,iy),kcolor
end sub

sub sclocate(irow,icol)
'-- locate with bounds checking turned on and cursor off
shared nrowchar, ncolchar
  if irow<1 then irow = csrlin
  if icol<1 then icol = pos(0)
  if irow>nrowchar then locate nrowchar,ncolchar:print :irow=nrowchar
  if icol>ncolchar then icol=ncolchar
  locate irow,icol,0
end sub

sub scprint(irow,icol,text$)
'-- print text TEXT$ at IX,IY (max to end of a line)
shared ncolchar
call sclocate(irow,icol)
print left$(text$,ncolchar-icol+1);
end sub

sub scprintv(irow,icol,text$)
'-- print text TEXT$ at IX,IY in color for variable
call setvarcol
call scprint(irow,icol,text$)
call settxtcol
end sub

sub scprintb(irow,icol,text$)
'-- print text TEXT$ at IX,IY in color for box
call setboxcol
call scprint(irow,icol,text$)
call settxtcol
end sub

sub scprintc(irow,icol,text$)
'-- print text TEXT$ at IX,IY in color for choice
call setchccol
call scprint(irow,icol,text$)
call settxtcol
end sub

sub scprintvfu(irow,icol,textuse$,a)
'-- print floating point variable A at IX,IY in color for variable using textuse$
call setvarcol
call sclocate(irow,icol)
  print using textuse$;a;
call settxtcol
end sub

sub scprompt(irow,icol,text$)
'-- print 1st char of text$ at ix,iy in choice text color,
'   followed by rest of text$ in normal color
call scprintc(irow,icol,left$(text$,1))
call scprint((0),(0),mid$(text$,2))
end sub

sub schit(a$)
'-- print text a$ on bottom line followed by Press any key...
shared nrowchar
call scprintc(nrowchar,1,a$):call scpause((0),(0))
end sub

sub sccenter(a$,n)
'-- center a text string on screen
  call scprint(n,(40-len(a$)/2),a$)
end sub

sub scpause (irow,icol)
'-- wait for a key
  call scprintc(irow,icol," Press any key to continue... ") :  d$=input$(1)
end sub

sub scbox(ntop,nbottom,ndprgnm,top$,bottom$)
'-- put menu box and label on screen
shared prgnm$
cls
call setboxcol
locate ntop,1
topln$=top$
if ndprgnm then
  if len(topln$)>0 then topln$=prgnm$+"   "+topln$ else topln$=prgnm$
end if
if len(topln$)>0 then
  print chr$(201);string$(74-len(topln$),205);" "+topln$+" ";string$(2,205);chr$(187);
else
  print chr$(201);string$(78,205);chr$(187);
end if

d&=fre(-1)   'check remaining memory
if d&<&h20000 then call scprint(1,4,"<"+fnstrnum$(d&\1024)+"Kb>")

for i=ntop+1 to nbottom-1
 call scprint(i,1,chr$(186))
 call scprint(i,80,chr$(186))
next

locate nbottom,1
call scprintb(nbottom,1,chr$(200)+string$(78,205)+chr$(188))
if len(bottom$)>0 then call scprintc(nbottom,75-len(bottom$)," "+bottom$+" ")
end sub

sub notalone
  print " This program cannot run alone; start with MT5."
  end
end sub

sub inkey (d$)
'-- get a key from the buffer
   do
    d$=inkey$
   loop while len(d$)=0
   if len(d$)=2 then d$=right$(d$,1) else d$=lcase$(d$)
end sub
'=======================================================================

'-- following subroutines for getting user inputs into defined input areas
function keyget (b$, dtype$,inson)
'-- gets a key b$ and returns ascii code of key hit (minus code if extended code)
'   includes type checking, and toggling insert on/off
' also displays time message if itimerow>0
' dtype$ must be lower case i,f,e or anything else for ascii

shared itimerow,itimecol
itype=instr("efi",dtype$)
DO
  DO
    if itimerow>0 and t$<>time$ then
      t$=time$
      irow=csrlin:icol=pos(0)  'gets posn
      call settxtcol:call scprint(itimerow,itimecol," Time "+t$+" ")
      call sethilcol:locate irow,icol,1
    end if
    b$ = INKEY$
  LOOP WHILE LEN(b$) = 0
  ikey=asc(b$)
  isdone=ktrue
  SELECT CASE ikey
    CASE 0  'extended character ok
      ikey = -asc(mid$(b$,2))
      b$=""
      if ikey=-82 then inson = (not inson)  'insert key
    CASE 3, 8, 9, 13, 25, 27   'Ctrl-C, Bksp, Tab, Enter, Ctrl-Y, Esc
      b$=""
    CASE 43, 45, 48 TO 57  ' "+", "-", "0" to "9" always ok
    CASE 46  ' decimal point
      isdone = (itype<3)
    CASE 69, 101   ' "e" or "E"
      isdone = (itype<2)
    CASE 32 TO 126  'all other ascii
      isdone = (itype=0)
    CASE ELSE
      isdone=kfalse
  END SELECT
LOOP until isdone
keyget=ikey
end function

FUNCTION iget% (irow, icol, kwait, a$, dtype$, nin)
'-- at irow & icol get response a$, n characters long
'    return ascii code of last key struck and inext=up/down code for next move
' if kwait then CR must be hit to confirm entry
  SHARED kscry, nrowchar, ncolchar

'-- dtype$ must be lower case
  isnum= (instr("ife",dtype$)>0)
7777  call sclocate(irow,icol)
  n=fnxlimit(nin,1,ncolchar-icol+1) ' adjust max string length if too long
  c$ = LEFT$(ltrim$(a$) + STRING$(n, 32), n)
  call sethilcol
  PRINT c$;

'-- initialize variables
  istrt = ktrue
  b$ = ""
  itab = 3
  kdely=(kscry+8)/nrowchar      'pixels per line
  icurnorm=kdely-1  'normal cursor is 2 scan lines high
  icurbig=kdely-3 'expanded cursor for insert on
  i=1 ':ibeg=1 'i is col of cursor, ibeg is index into string for 1st char displayed
  DO
    IF inson THEN icur = icurbig ELSE icur = icurnorm
    LOCATE , icol + i - 1, 1, icur, kdely
    ikey= keyget(b$,dtype$,inson)
    isasc=fninlim(ikey,32,126)
    IF (istrt and isasc) or ikey=25 then
' a regular key or Ctrl-Y hit - blank the rest of the line
      mid$(c$,i) = STRING$(n, 32): LOCATE , icol: PRINT c$;
    ELSEIF not isasc then
     SELECT CASE ikey
       CASE -81,-80,-73,-72,3,13,27 '-- exit for CR,ESC,Up,Pgup,Down,Pgdn,Ctrl-C
         exit do
       CASE 8   'bksp
         imov = -1
       CASE 9   'tab
         imov = itab
       CASE -15  'shift tab
         imov = -itab
       CASE -71 'home
         imov=-i+1
       CASE -75 'left
         imov=-1
       CASE -77 'right
         imov=1
       CASE -79 'end; just goes to end of text
         imov = LEN(RTRIM$(c$)) - i + 1
       CASE ELSE
         imov=0
     END SELECT
     if isnum then imax=fnamin(LEN(rtrim$(c$))-i+1,n-i) else imax=n-i
     imov = fnxlimit (imov,-i+1,imax)
     i = i + imov
' delete or backspc
     IF ikey=-83 or ikey=8 THEN c$ = LEFT$(c$, i - 1) + MID$(c$, i + 1) + " ": LOCATE , icol: PRINT c$;
    END IF
    istrt = 0
    LOCATE , icol + i - 1
    IF LEN(b$) > 0 THEN
      IF inson THEN
        mid$(c$,i) = b$ + MID$(c$, i): PRINT MID$(c$, i);
      ELSE
        PRINT b$; : MID$(c$, i, 1) = b$
      END IF
      if not (kwait and i=n) then i=i+1
    END IF
  LOOP UNTIL i > n  'here is where field scrolling could be allowed
  IF ikey=27 THEN
  ELSEIF istype(c$,dtype$) then  'do rigorous type checking
    a$=c$
  ELSE
    if displaymsg$(alltrim$(c$)+" is invalid entry (Esc to abort)","",kfalse)<>chr$(27) then 7777
  END IF
  call scprintv(irow,icol, a$)
  iget=ikey
END FUNCTION

function istype(a$,dtype$)
'-- dtype$ must be lcase
istype=ktrue
if instr("fie",dtype$)>0 then
  b$=lcase$(alltrim$(a$))
  l$=numlist$(0,9)
  c$=left$(b$,1)
  if instr(l$+"+-.",c$)=0  then
    istype=kfalse
  elseif dtype$="i" and c$="." then
    istype=kfalse
  else
    e$="e"
    d$="."
    f$=""
    for i=2 to len(b$)
      c$=mid$(b$,i,1)
      if dtype$="i" then
        if instr(l$,c$)=0 then istype=kfalse:exit for
      elseif dtype$="f" then
        if instr(l$+d$,c$)=0 then istype=kfalse:exit for
      else  'dtype$="e"
        if instr(l$+d$+e$+f$,c$)=0 then istype=kfalse:exit for
      end if
      if c$="." then d$=""
      if c$="e" then f$="+-":e$="":d$=""
      if f$="+-" then f$=""
    next
  end if
end if
end function

function input1 (irow,icol,prompt$,a$,kwait,kokesc,kerase,okans$)
'-- displays prompt at irow,icol,and gets a single key (requires CR if kwait)
'    key hit must be in okans$ list (or esc if kokesc)
'    returns last key hit and lowercase of key hit in a$
'-- okans$ must be lcase
lp=len(prompt$)
arrow$="->"
IF lp>0 THEN
  call scprintc(irow,icol,prompt$+arrow$)
  la=len(arrow$)
ELSE
  call sclocate(irow,icol)
  la=0
END IF
nrow=0:ncol=0
DO
  ik=iget(nrow,ncol,kwait,a$,"a",1):a$=lcase$(a$)
loop until (ik<>27 and instr(okans$,a$)>0) or (kokesc and ik=27)
if kerase then
  call cnchar(irow,icol,lp+la+1) 'erase prompt and answer
elseif lp>0 then
  call cnchar(irow,icol+lp,la)  'just erase arrow
endif
input1=ik
end function

function inputn (irow,icol,prompt$,a$,kwait,kokesc,kerase,dtype$,n)
'-- displays prompt at irow,icol,and gets up to n key answer
'    requires CR if kwait or allows esc if kokesc
'    dtype$ is data type:a=ascii,f=float,i=integer,e=exponential
'    returns last key hit and answer in a$

lp=len(prompt$)
arrow$="->"
IF lp>0 THEN
  call scprintc(irow,icol,prompt$+arrow$)
  la=len(arrow$)
ELSE
  call sclocate(irow,icol)
  la=0
END IF
nrow=0:ncol=0
DO
  ik=iget(nrow,ncol,kwait,a$,dtype$,n)
loop until istype(a$,dtype$) or (ik=27 and kokesc)
if kerase then
  call cnchar(irow,icol,lp+la+n) 'erase prompt and answer
else
  call cnchar(irow,icol+lp,la)  'just erase arrow
endif
a$=rtrim$(a$)
inputn=ik
end function

function alltrim$(a$)
  alltrim$=ltrim$(rtrim$(a$))
end function

function numlist$(nl,nu)
  a$=""
  for nn = nl to nu
    a$=a$+alltrim$(str$(nn))
  next
  numlist$=a$
end function

'=======================================================================
'-- following subroutines to provide more expanded fileopening capabilities

SUB parsefilenm (inputname$, fldrv$, flpth$, flnm$,flext$)
inputname$ = LTRIM$(RTRIM$(UCASE$(inputname$)))
istrt = INSTR(inputname$, ":")
fldrv$ = LEFT$(inputname$, istrt)
flpth$ = MID$(inputname$, istrt + 1)
istrthld = INSTR(flpth$, "\")
istrt = istrthld
DO WHILE istrthld > 0
   istrt = istrthld
   istrthld = INSTR(istrthld + 1, flpth$, "\")
LOOP
flnm$ = MID$(flpth$, istrt + 1)
flpth$ = LEFT$(flpth$, istrt)
istrt=instr(flnm$,".")
if istrt=0 then istrt=len(flnm$)+1
flext$=mid$(flnm$,istrt)
flnm$=left$(flnm$,istrt-1)
END SUB

SUB sort (numitem%, listitem$())
limit = numitem
DO
   iswaps = 0
   FOR i = 1 TO (limit - 1)
   ' if two adjacent elements are out of order, swap them
      IF listitem$(i) > listitem$(i + 1) THEN
         SWAP listitem$(i), listitem$(i + 1)
         iswaps = i
      END IF
      IF INKEY$ = CHR$(27) THEN EXIT DO
   NEXT
   IF iswaps>0 THEN limit = iswaps
LOOP WHILE iswaps>0
END SUB

SUB buildname(openfilename$,drv$,filepath$,nm$,fex$,dfltfex$)
'assembles a complete filename w/wildcards if needed

CALL parsefilenm(openfilename$, drv$, filepath$, nm$,fex$)
IF LEN(nm$) = 0 THEN nm$="*"
IF INSTR(fex$, ".") = 0 THEN
   IF dfltfex$ = "" then
      fex$=".*"
   else
      fex$=dfltfex$
   END IF
END IF
openfilename$ = drv$ + filepath$ + nm$ + fex$

END SUB


FUNCTION GetCurDrv$

DIM InRegs AS RegType, OutRegs AS RegType
InRegs.ax = &H1900
CALL INTERRUPT(&H21, InRegs, OutRegs)
GetCurDrv$ = CHR$(65 + (&HFF AND OutRegs.ax)) + ":"

END FUNCTION

FUNCTION displaymsg$ (msgin$, okans$, longans%)
' Input Parameters: msgin$=an error message to display (? is tagged on)
'                   okans$ =acceptable keys to hit
'                           ("" gives "Press any key" prompt)
'                   longans%: true=alternate form where answer may be input
' Function returns ucase$(character of the key hit)

SHARED NROWCHAR, NCOLCHAR

isbeepon = 0

' Save current cursor position
irowsave = CSRLIN
icolsave = POS(0)

IF longans THEN
   msg$ = SPACE$(NCOLCHAR - 2)
ELSE
   IF isbeepon THEN BEEP
   IF LEN(okans$) = 0 THEN
      msg$ = " " + msgin$ + ". Press any key... "
   ELSE
      msg$ = " " + msgin$ + "? ("+left$(okans$,1)
      FOR i = 2 TO LEN(okans$)
         msg$ = msg$ + "/"+MID$(okans$, i, 1)
      NEXT
      msg$ = msg$ + ") >"
   END IF
END IF
nln = 1
ibrk = NCOLCHAR - 2
IF LEN(msg$) > NCOLCHAR - 2 THEN        'need a two line display
   FOR ibrk = NCOLCHAR - 2 TO 1 STEP -1  'find nearest space to end of line
      IF INSTR(" .?)", MID$(msg$, ibrk, 1)) <> 0 THEN EXIT FOR
   NEXT
   nln = 2
END IF
REDIM linesave$(nln)
call sethilcol
' save line and print error msg
FOR ln = 1 TO nln
   linesave$(ln) = ""
   FOR i = 1 TO 80    ' saves entire row
      linesave$(ln) = linesave$(ln) + CHR$(SCREEN(NROWCHAR + ln - 2, i))
   NEXT i
   LOCATE NROWCHAR + ln - 2, 1
   PRINT MID$(msg$, 1 + (ln - 1) * ibrk, ibrk);
NEXT ln
msgans$ = ""
IF longans THEN
   LOCATE NROWCHAR - 1, 1: PRINT " " + msgin$ + " >";
   INPUT ; "", msgans$
ELSE
   DO
      msgans$ = UCASE$(INPUT$(1))
      IF ((LEN(okans$) > 0) AND (INSTR(UCASE$(okans$), msgans$) = 0)) THEN msgans$ = ""
   LOOP UNTIL LEN(msgans$) > 0
END IF
call settxtcol
FOR ln = 1 TO nln
   LOCATE NROWCHAR + ln - 2, 1
   PRINT linesave$(ln);
NEXT
' Go back to saved cursor position
LOCATE irowsave, icolsave
displaymsg$ = msgans$

END FUNCTION

FUNCTION displayerr$ (ierrnum%, openfilename$, igetans%)
SELECT CASE ierrnum
   CASE 2
      errstring$ = "File not found"
      ierrtype = 2
   CASE 3
      errstring$ = "Path not found"
      ierrtype = 1
   CASE 5
      errstring$ = "Access denied"
      ierrtype = 2
   CASE &HF
      errstring$ = "Invalid drive"
      ierrtype = 1
' Max &h12 error codes returned by find file DOS function
   CASE &H15
      errstring$ = "Drive not ready"
      ierrtype = 2
   CASE ELSE
      errstring$ = "Unknown"
      ierrtype = 2
END SELECT
errstring$ = errstring$ + " error with " + UCASE$(openfilename$) + ". "
chc$ = ""
   IF ierrtype > 1 THEN
      errstr2$ = "Retry/"
      chc$ = "R"
   END IF
   IF igetans = ktrue THEN
      errstr2$ = errstr2$ + "Filename change/"
      chc$ = chc$ + "F"
   END IF
   IF chc$<>"" THEN
      chc$ = chc$ + "A"
      errstring$ = errstring$ + errstr2$ + "Abort"
   ENDIF
ans$ = displaymsg$(errstring$, chc$, 0)
IF ans$ = "F" THEN CALL getflnm(openfilename$)
displayerr$ = ans$

END FUNCTION

FUNCTION IChkDrv% (drv$, openfilename$, isoknmchg%)  'returns true if drive is removable
'and ready, if a fixed drive, or if a network device; gives "drive not ready"
'message if not true (unless isoknmchg = 1 meaning "no messages allowed")

DIM InRegs AS RegType, OutRegs AS RegType
numreread = 1      'number of retries for floppy drive being ready
DO
   CALL parsefilenm(openfilename$, drv$, pth$, flnm$,flex$)
' get default drive if no drive specified
   IF LEN(drv$) = 0 THEN drv$ = GetCurDrv$
   numdrv = ASC(drv$) - 65
' Now have drive letter and number for openfilename

' See if drive is removable
   InRegs.ax = &H4408
   InRegs.bx = numdrv + 1
   CALL INTERRUPT(&H21, InRegs, OutRegs)
   itest = OutRegs.ax
'itest: 0=removable,1=fixed drive or network device,&HF=invalid drive specified

   IF itest = 0 THEN       'drive removable, itest it for being ready
      FOR i = 1 TO numreread
         InRegs.dx = numdrv          'drive A=0,B=1,C=2,...; head 0
         InRegs.cx = 1               'sector 1, track 0
         InRegs.ax = &H401           ' Verify Sector to see if drive ready
         CALL INTERRUPT(&H13, InRegs, OutRegs)  ' BIOS call
         itest = (&H1 AND OutRegs.flags) = 0
         IF itest THEN EXIT FOR ELSE itest = &H15 'drive not ready
         InRegs.ax = 0           ' Reset controller
         CALL INTERRUPT(&H13, InRegs, OutRegs)  ' BIOS call
      NEXT
   ELSEIF itest = 1 THEN  'drive is fixed or network device; don't check it
      itest = -1
   END IF

   IF (itest <> -1 and isoknmchg < 1) THEN
      msgans$ = displayerr$(itest, openfilename$, isoknmchg)
      IF msgans$ = "A" THEN itest = 1 ELSE itest = 0
   END IF
LOOP WHILE itest = 0
IChkDrv = (itest = -1)

END FUNCTION

SUB getflnm (oldname$)
msg$ = "Enter new filename"
ans$ = displaymsg$(msg$, "", -1)
IF ans$ <> "" THEN
   CALL parsefilenm(oldname$, drv$, flpth$, nm$,dfltfex$) 'gets the extension
   CALL buildname(ans$,drv$,flpth$,nm$,fex$,dfltex$)
   oldname$ = ans$
ENDIF
END SUB

FUNCTION listselect% (listcnt%, listitem$(), topmsgin$, issorted%)
SHARED NROWCHAR, NCOLCHAR

isbeepon=kfalse
topmsg$=topmsgin$  'makes a modifiable message

' Following constants may be changed to suit
linebottom = NROWCHAR - 1
ix1 = 1
iy1 = 1
icolwidth = 30  'width of file selection box
vl$ = CHR$(179)
hl$ = CHR$(205)
ul$ = CHR$(213)
ur$ = CHR$(184)
ll$ = CHR$(212)
lr$ = CHR$(190)

ix2 = ix1 + icolwidth - 1
lnavail = linebottom - iy1 + 1
IF listcnt < lnavail - 1 THEN
   iy2 = listcnt + iy1 + 1
   ishowcnt = listcnt
ELSE
   iy2 = linebottom
   ishowcnt = lnavail - 2
END IF
irowsave = CSRLIN
icolsave = POS(0)
lncnt = iy2 - iy1 + 1
icolcnt = ix2 - ix1 + 1
DIM linesave$(lncnt)
FOR ln = 1 TO lncnt
   linesave$(ln) = ""
   icurln = iy1 + ln - 1
   FOR icol = ix1 TO ix2
      linesave$(ln) = linesave$(ln) + CHR$(SCREEN(icurln, icol))
   NEXT icol
   LOCATE icurln, ix1: PRINT SPACE$(icolcnt);
NEXT ln
' draw box
call sethilcol
LOCATE iy1, ix1: PRINT ul$; STRING$(icolcnt - 2, hl$); ur$;
FOR ln = 2 TO lncnt - 1
   LOCATE iy1 + ln - 1, ix1: PRINT vl$;
   LOCATE , ix2: PRINT vl$;
NEXT
LOCATE iy2, ix1: PRINT ll$; STRING$(icolcnt - 2, hl$); lr$;

IF LEN(topmsg$) > icolcnt - 2 THEN topmsg$ = ".." + RIGHT$(topmsg$, icolcnt - 4)
LOCATE iy1, (ix2 + ix1 - LEN(topmsg$) - 1) / 2
PRINT topmsg$;
bttmmsg$ = "<Esc=QUIT"
IF ishowcnt < listcnt THEN
   bttmmsg$ = bttmmsg$ + ";Pgdn=MORE"
END IF
LOCATE iy2, (ix2 + ix1 - LEN(bttmmsg$) - 1) / 2
PRINT bttmmsg$ + ">";
call settxtcol

IF issorted THEN
   LOCATE iy1 + 1, ix1 + 1
   PRINT " Sorting...";
   CALL sort(listcnt, listitem$())
END IF
LOCATE iy1 + 1, ix1 + 1
PRINT SPACE$(icolcnt - 2);

FOR ln = 1 TO ishowcnt
   LOCATE iy1 + ln, ix1 + 1
   PRINT listitem$(ln);
NEXT

ishowelement = 1
ihldelement = ishowelement
ichoose = 0
iredraw = 0
itopelement = 1
maxtop = listcnt - ishowcnt + 1
foundstr$=""
DO WHILE ichoose = 0
  LOCATE iy1 + ishowelement - itopelement + 1, ix1 + 1
  call setmsgcol
  PRINT listitem$(ishowelement);
  call settxtcol
  ihldelement = ishowelement
  iaction = kfalse
  WHILE iaction = kfalse
    do
      ans$ = INKEY$
    loop WHILE LEN(ans$) = 0
'    IF ASC(ans$) = 0 THEN ans$ = RIGHT$(ans$, 1)
    ia=ASC(ans$)
    IF ia = 0 THEN ia = -ASC(RIGHT$(ans$, 1))
    iaction = ktrue
    lndelta = 0
    ipgdelta = 0
'     SELECT CASE ASC(ans$)
    ichtmp=kfalse
    SELECT CASE ia
      CASE 13  'selected
        ichoose = ishowelement
      CASE 27  'esc = abort
        EXIT DO
      CASE -72  'up
        lndelta = -1
      CASE -80  'down
        lndelta = 1
      CASE -71  'home
        lndelta = -listcnt
      CASE -79  'end
        lndelta = listcnt
      CASE -73  'pgup
        ipgdelta = -1
        lndelta = -ishowcnt
      CASE -81  'pgdn
        ipgdelta = 1
        lndelta = ishowcnt
      CASE ELSE
        iaction = kfalse
        IF ia>31 THEN 'try to find next matching character
          foundstr$=foundstr$+ucase$(ans$)
          lenfndstr=len(foundstr$)
          for i=1 to listcnt
            if ucase$(left$(listitem$(i),lenfndstr))=foundstr$ then exit for
          next
          if i<=listcnt then  'found a match
            lndelta=i-ishowelement
            iaction = ktrue
            ichtmp=ktrue
          else
            if isbeepon then beep
          end if
        END IF
    END SELECT
    if not ichtmp then foundstr$=""
    ishowelement = ishowelement + lndelta
    IF ishowelement < 1 THEN
      ishowelement = 1
    ELSEIF ishowelement > listcnt THEN
      ishowelement = listcnt
    END IF
    IF ishowelement < itopelement THEN
      IF ipgdelta < 0 THEN
        itopelement = itopelement - ishowcnt
      ELSE
        itopelement = ishowelement
      END IF
      IF itopelement < 1 THEN itopelement = 1
      iredraw = -1
    ELSEIF ishowelement > itopelement + ishowcnt - 1 THEN
      IF ipgdelta > 0 THEN
        itopelement = itopelement + ishowcnt
      ELSE
        itopelement = ishowelement - ishowcnt + 1
      END IF
      IF itopelement > maxtop THEN itopelement = maxtop
      iredraw = -1
    END IF
  WEND
  IF iredraw THEN
    FOR ln = 1 TO ishowcnt
      LOCATE iy1 + ln, ix1 + 1
      PRINT LEFT$(listitem$(ln + itopelement - 1) + SPACE$(icolcnt), icolcnt - 2);
    NEXT
    iredraw = 0
  ELSE
    LOCATE iy1 + ihldelement - itopelement + 1, ix1 + 1
    PRINT listitem$(ihldelement);
  END IF
LOOP
FOR ln = 1 TO lncnt
  LOCATE iy1 + ln - 1, ix1
  PRINT linesave$(ln);
NEXT ln
listselect = ichoose
LOCATE irowsave, icolsave

END FUNCTION

FUNCTION ifilefind% (flnm$, nextfile%, fullflnm$, isoknmchg%)
' Fairly low level routine that does BIOS calls "findfirst" and "findnext"
' Input Parameter:   flnm$ = filename (including path and wildcard characters)
'                    nextfile<>0 if find next matching file (after first calling
'                               function with nextfile=0
'                    isoknmchg%= -1 if name change allowed
'                              =  0 if not, but error messages given
'                              = +1 if no error messages to be given
'                              = +2 if unambiguous file name not allowed
' function returns 0 if file exists, ierrnum if not exist
' Output Parameters: fullflnm$ = full filename (w/drive and path)

'shared ierrnum as integer
DIM InRegsX AS RegTypeX, OutRegsX AS RegTypeX
DO
' Set ds and es values so that they will default to current values
   InRegsX.ds = -1
   InRegsX.es = -1
   itest = -1
   msgans$ = ""

' Put null byte on end of filename for DOS interrupt call
   openfile$=flnm$
   call parsefilenm(openfile$,fldrv$,flpth$,(""),(""))
   openfile$ = flnm$ + CHR$(0)

   IF nextfile THEN
      InRegsX.ax = &H4F00
   ELSEIF IChkDrv(drv$, flnm$, isoknmchg) THEN 'Also returns drive letter if not previously specified
      InRegsX.ax = &H4E00
      InRegsX.dx = SADD(openfile$)
      InRegsX.cx = 0
   ELSE
      itest = 0
      ierrnum=&H15
   END IF

   IF itest THEN
      CALL INTERRUPTX(&H21, InRegsX, OutRegsX)
      ierrnum = (OutRegsX.ax AND &HF)
      itest = ((&H1 AND OutRegsX.flags) = 0)
      IF itest THEN
' what is actual full filename found?
' Set ds and es values so that they will default to current values
         ierrnum=0
         InRegsX.ds = -1
         InRegsX.es = -1
         InRegsX.ax = &H2F00  ' Get Disk Transfer Address
         CALL INTERRUPTX(&H21, InRegsX, OutRegsX)
' es:bx now point to current DTA
         DEF SEG = OutRegsX.es
         fullflnm$ = ""
         FOR noffset = OutRegsX.bx + 30 TO OutRegsX.bx + 41
            nxtasc = PEEK(noffset)
            IF nxtasc = 0 THEN EXIT FOR
            fullflnm$ = fullflnm$ + CHR$(nxtasc)
         NEXT
         DEF SEG
         fullflnm$=fldrv$+flpth$+fullflnm$  'add drv and path
      ELSEIF (nextfile = 0 and isoknmchg<1) THEN
         msgans$ = displayerr$(ierrnum, flnm$, isoknmchg)
      END IF
   END IF
LOOP WHILE msgans$ = "R" OR msgans$ = "F"
ifilefind=ierrnum

END FUNCTION

FUNCTION isfile% (flnm$, fullflnm$, isoknmchg%)

' Input Parameter:   flnm$ = filename (including path and wildcard characters)
'                    isoknmchg%= -1 if name change allowed
'                              =  0 if not, but error messages given
'                              = +1 if no error messages to be given
'                              = +2 if ambiguous file name not allowed
' function returns 0 if file exists, ierrnum if not exist
' ierrnum =-2 if isoknmchg=2 and ambiguous file name given
' Output Parameters: fullflnm$ = full filename (w/drive and path)
fullflnm$=flnm$
isfile= ifileselect(fullflnm$,isoknmchg)

END FUNCTION

FUNCTION ifileselect% (openfilename$, isoknmchg%)
' Checks for file existing (including select from list if wildcards in
'  openfilename$ and isoknmchg=-1
' Returns 0 if ok to open file, otherwise returns errornumber

ierrnum=ifilefind(openfilename$, 0, fullname$, isoknmchg)
IF ierrnum=0 THEN
   IF INSTR(openfilename$, "*") + INSTR(openfilename$, "?") > 0 THEN
     idspy=ktrue
     IF isoknmchg=2 THEN ierrnum=-2 'new error code indicating ambiguous file name not allowed
     IF isoknmchg<>-1 THEN_
      idspy=(displaymsg$("* or ? not allowed here. Show list of "+openfilename$+" files","YN",kfalse)="Y")
     IF idspy THEN
' openfilename$ contains a DOS "wildcard" character
       listcnt = 1
       WHILE ifilefind(openfilename$, -1, fullname$, isoknmchg)=0
         listcnt = listcnt + 1
       WEND
       DIM filelist$(listcnt)
       FOR i = 1 TO listcnt
         itest = (ifilefind(openfilename$, i>1, fullname$, isoknmchg)=0)
         call parsefilenm(fullname$,fldrv$,flpth$,flnm$,flext$)
         filelist$(i) = flnm$+flext$
       NEXT
' now have a (sorted) array of filenames to choose from
       drvfilepath$=fldrv$+flpth$
       ichc = listselect(listcnt, filelist$(), openfilename$, -1)
       IF ichc > 0 AND isoknmchg=-1 THEN
         openfilename$ = drvfilepath$ + filelist$(ichc)
       ELSEIF isoknmchg<>2 THEN
         ierrnum=2  'file not found
       END IF
     END IF
'  ELSE  ' unambiguous file found
   END IF
END IF
ifileselect = ierrnum

END FUNCTION

FUNCTION ifopen% (openmode$, numopenfile%, openfilename$, lenopenrec%, isoknmchg%)
' Input parameters same as OPEN command where MODE is "O","I","R","B",or "A"
' except that wildcard characters in filespec allow selection from list
' returns 0 if not open and nonzero if open

' parse filename and expand with wildcards if not specified
   CALL buildname(openfilename$,drv$,filepath$,nm$,fex$,"")
   ifileopen = ifileselect(openfilename$, isoknmchg)=0 'true means everything checks out

' Only check not done is for write protected disk;
'    this will cause DOS error message and probable program bomb unless
'    Retry is successful.

   IF ifileopen THEN
      IF numopenfile = 0 THEN numopenfile = FREEFILE
      IF lenopenrec > 0 THEN
        OPEN openmode$, numopenfile, openfilename$, lenopenrec
      ELSE
        OPEN openmode$, numopenfile, openfilename$
      END IF
   END IF
ifopen = ifileopen

END FUNCTION

function igetyesno(irow,icol,prompt$,dflt$,isokescape) 'gets yes or no, allows esc, erases when done
 a$=dflt$
 if input1(irow,icol,prompt$+" (Y/N)",a$,kfalse,isokescape,ktrue,"yn")=27 then
   igetyesno=kfalse
 else
   igetyesno=(a$="y")
 end if
end function

'--- check the quit command ---
sub endcheck
shared prgnm$
 y$=""
 if input1((0),(0),"Quit "+prgnm$+". Are you sure? (Y/N required)",y$,kfalse,ktrue,ktrue,"yn")<>27 then
  if y$<>"y" then exit sub
  END
 end if
end sub

function lftvar$(strgvar$)  'returns leftmost variable from strgvar$; strgvar$ is shortened
' separators checked for are comma, space, and equal sign; no check for quote-delimited strings yet
' throw away everything to right of single quote '
ls=instr(strgvar$,"'")
if ls>0 then strgvar$=left$(strgvar$,ls-1)
strgvar$=alltrim$(strgvar$)
'-- discard leading "=" sign
if left$(strgvar$,1)="=" then strgvar$=ltrim$(mid$(strgvar$,2))
ls=len(strgvar$)
if ls>0 then
   for i=1 to ls
      if instr(", =",mid$(strgvar$,i,1))>0 then exit for
   next
   lft$=left$(strgvar$,i-1)
   strgvar$=ltrim$(mid$(strgvar$,i))
   if instr(",",left$(strgvar$,1))>0 then strgvar$=mid$(strgvar$,2) 'take off leading comma
else
   lft$=""
end if
lftvar$=lft$
end function

function keepold(flnm$) 'returns true if existing file should be kept
' requires unambiguous file name
' revised Jan 24, 1993 to explicitly check for ambiguous file name
   kpold=isfile(flnm$,(""),2)
   if kpold = 0 then
     keepold = (displaymsg$("Overwrite existing file "+flnm$,"YN",kfalse)<>"Y")
   else
     keepold= (kpold = -2)  'ambiguous file name given
   end if
end function

sub memcheck
d&=fre(-1)   'check remaining memory
irow=csrlin:icol=pos(0):call setboxcol
if d&<&h20000 then
  call scprint(1,4,"<"+fnstrnum$(d&\1024)+"Kb>อออ")
else
  call scprint(1,4,string$(8,205))
end if
locate irow,icol
call settxtcol
end sub

rem $INCLUDE: 'm5synio'

