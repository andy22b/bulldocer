'---  M5FUNCTN.BAS
'DEFINT I-N

'-- conversions
 DEF FNDEG(x) = x*57.295779
 DEF FNRAD(x) = x*0.017453293

'-- inverse trig functions

 DEF FNASIN(x)
     if abs(x)<1  then fnasin=atn(x/sqr(1-x*x))
     if abs(x)>=1 then fnasin=piihalf
  end def

 DEF FNACOS(x)
     if abs(x)<1  then fnacos=piihalf-atn(x/sqr(1-x*x))
     if abs(x)>=1 then fnacos=0
  end def

 DEF FNATAN2(Y,X)
    if X>0            then fnatan2=ATN(Y/X)
    if (Y>0 AND X<0)  then fnatan2=(pii+ATN(Y/X))
    if (Y<=0 AND X<0) then fnatan2=(ATN(Y/X)-pii)
    if (X=0 AND Y>=0) then fnatan2=piihalf
    if (X=0 AND Y<0)  then fnatan2=-piihalf
  end def

'-- logarithms
 DEF FNALOG10(x) = log(x) * 0.4342945
 DEF FNALOG2(x)  = log(x) * 1.442695

'-- min/max
 DEF FNAMAX(x,y) = -(x>=y)*x -(y>x)*y
 DEF FNAMIN(x,y) = -(x<=y)*x -(y<x)*y

'-- responses to queries
 DEF FNKYES(Y$) = ucase$(left$(Y$,1)) = "Y"
 DEF FNKANSR(Y$) = instr("YN",ucase$(left$(Y$,1))) > 0

 DEF FNYESNO$(k)
    fnyesno$=" NO"
    if k then fnyesno$="YES"
 end def

DEF FNIGETYESNO(isokescape)  'wait for user to hit Y,N (or optional ESC),true if Y hit
   static a$
   do
      a$=ucase$(input$(1))
      if isokescape then
         if a$=chr$(27) then exit do
      end if
   loop until instr("YN",a$)
   fnigetyesno=(a$="Y")
end def

CRIN$=" (value <0 to set to 0, 0 for no change)"
DEF FNCRIN(x)
'-- input and change if >0, set to 0 if <0, no change if 0
static y
   input y
   if y<0 then fncrin=0
   if y=0 then fncrin=x
   if y>0 then fncrin=y
end def

'-- formatting output
DEF FNDFORM$(I,J)=STRING$(I-J-1,35)+"."+STRING$(J,35)+"^^^^"
DEF FNFFORM$(I,J)=STRING$(I-J-1,35)+"."+STRING$(J,35)
DEF FNIFORM$(I)=STRING$(I,35)
DEF FNAFORM$(I)
  if i=1 then
    fnaform$ = "!"
  else
    fnaform$=  "\" + SPACE$(I-2) + "\"
  end if
end def
DEF FNFIX(x,n)= int(x*10^n+0.5)/(10^n)

'-- convert number to string
DEF FNSTRNUM$(x) =ltrim$(str$(x))
'-- get length of number as string
DEF FNLENS(x) = len(ltrim$(str$(x)))
'-- is a between b and c?
DEF FNINLIM(a,b,c)=(a>=b) and (a<=c)
'-- distance between X1,Y1 and X2,Y2
DEF FNDIST(x1,x2,y1,y2) = sqr((x1-x2)^2 + (y1-y2)^2)
'-- value of number in string A$ starting in column I and extending J spaces
DEF FNVALUE(I,J) = val(mid$(a$,i,j))
'-- string in A$ starting in column I and extending J spaces
DEF FNVALUE$(I,J) = mid$(a$,i,j)

'-- force a to be between b and c
DEF FNXLIMIT(a,b,c)
   fnxlimit=a
   if a<b then fnxlimit=b
   if a>c then fnxlimit=c
end def

'=======================================================================
'--- moving the cursor around, number pad keys
'-- must be used with INKEY subroutine below to handle extended characters

def fnmovend(x$) = (x$="1") OR (x$=chr$(79))
def fnmovedn(x$) = (x$="2") OR (x$=chr$(80))
def fnipagedn(x$) = (x$="3") OR (x$=chr$(81))
def fnmovelt(x$) = (x$="4") OR (x$=chr$(75))
def fnmiddle(x$) = (x$="5")
def fnmovert(x$) = (x$="6") OR (x$=chr$(77))
def fnmovehm(x$) = (x$="7") OR (x$=chr$(71))
def fnmoveup(x$) = (x$="8") OR (x$=chr$(72))
def fnipageup(x$) = (x$="9") OR (x$=chr$(73))
def fniplus(x$)   = (x$="+") OR (x$=chr$(78))
def fnminus(x$)  = (x$="-") OR (x$=chr$(74))

'--- function keys
'def fnkeyf1(x$) = (x$=chr$(59))
'def fnkeyf2(x$) = (x$=chr$(60))
'def fnkeyf3(x$) = (x$=chr$(61))
'def fnkeyf4(x$) = (x$=chr$(62))
'def fnkeyf5(x$) = (x$=chr$(63))
'def fnkeyf6(x$) = (x$=chr$(64))
'def fnkeyf7(x$) = (x$=chr$(65))
'def fnkeyf8(x$) = (x$=chr$(66))
'def fnkeyf9(x$) = (x$=chr$(67))
'def fnkeyf10(x$)= (x$=chr$(68))
'=======================================================================
'-- checking for the Esc key
def FNIESC(x$) = (x$=chr$(27))
def FNINESC = (inkey$=chr$(27))
'=======================================================================
'-- checking for the q key to quit

def FNIQUIT(x$) = (x$="q") OR (x$="Q")

DEF FNINQUIT(kwait)
'-- check for the q key, kwait tells whether to pause or not
static d$
  if kwait then d$=input$(1) else d$=inkey$
  fninquit=fniquit(d$)
end def
'=======================================================================
'def fnokcol% = kscreen<>2 and kscreen<>3 and kscreen<>11  ' ok to use color
'=======================================================================
' Add an extension to the end of a filename if no extension present
DEF FNADDEXT$(flnm$,extn$)
  if instr(flnm$,".")=0 then fnaddext$=flnm$+"."+extn$ else fnaddext$=flnm$
end def

' Remove extension from end of a filename
DEF FNSUBEXT$(flnm$)
  static k
  k=instr(flnm$,".")
  if k>0 then fnsubext$=left$(flnm$,k-1) else fnsubext$=flnm$
end def
'=======================================================================
' Get string of instrument calibration type
DEF FNCALTYP$(instin)
  if instin<0 OR (instin AND &H7FE0) then 'anything over bit 4
    fncaltyp$="unknown"
  elseif instin AND 16 then 'bit 4
    fncaltyp$="acceleration"
  elseif instin AND 8 then  'bit 3
    fncaltyp$="velocity"
  else
    fncaltyp$="displacement"
  end if
end def

