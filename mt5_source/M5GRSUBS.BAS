'=============================================================================
' PLOTTING SUBROUTINES
'-- these routines create a plot file in the Golden Software PLOTCALL
'    or PostScript format and also plot on the screen.
'   if IROTATE is true then the x-axis is along the long edge of the paper
'    else it is along the short edge
'   KOPEN7 is a flag so program knows file 7 is open
'   GSCALE is a scale factor that accommodates a smaller page; the 25.4 factor
'    converts from mm to inches
'   started modifying for HPGL on 10/15/90 by Peter Zwick, see igraph=3 entries
'   version completed 1-23-91 used in SYN4
'   modified to correct anamolies in DOT macro and gstext for PS on 3-10-91
'   modified hpcolor to lift pen before selecting a new one - 4-10-91
'   modified to include igraph=0 for PlotCall Version 2.0 and higher Dec,91
'   1-13-94 added % to list of special characters in sub clnpar
'-----------------------------------------------------------------------------

sub gsopen (f$) static

'-- open F$ as file number 7 for output of the plot commands
  shared kopen7
  kopen7=ktrue
  open "o",7,f$
end sub
'-----------------------------------------------------------------------------
sub gsinit (f$) static
shared xmagp,lcode,khpcol,igraph,xoffhp,yoffhp
'-- initialize F$ for output of the plot commands
  if igraph=2 then
    print #7, "%!"
    print #7, "/M { moveto } def"
    print #7, "/D { lineto } def"
    print #7, "/T { translate } def"
    print #7, "/RO { rotate } def"
    print #7, "/S { show } def"
    print #7, "/SN { stroke newpath } def"
    print #7, "/SK { stroke } def"
    print #7, "/SL { setlinewidth } def"
    print #7, "/NP { newpath } def"
    print #7, "/CP { closepath } def"
    print #7, "/OPENLET {true charpath gsave 1 setgray fill grestore SK} def"
    print #7, "/DOT { M gsave [] 0 setdash currentpoint NP 1.5 currentlinewidth mul 0 360 arc fill grestore} def"
    print #7, "/CF {/Courier findfont exch scalefont setfont} def
    print #7, "/TRF {/Times-Roman findfont exch scalefont setfont} def
    print #7, "0.5 SL"
  elseif igraph>2 then
    xoffhp=0:yoffhp=0
    xmagp=1021
    if igraph=3 then print #7,chr$(27);"E";chr$(27);"%1B"  ' 6bytes - HP-LJ III "enter HPGL mode"
    if igraph<5 then
      print #7, "IN;PS4;PA;FT2;PT.2";  '18 bytes+CR/LF, 8-1/2 x 11
    else
      xmagp=1329
      print #7, "IN;PS0;PA;FT2;PT.3";  '18 bytes+CR/LF, 11 x 17
    end if
    lcode=0:khpcol=0
  end if
end sub

'-----------------------------------------------------------------------------
sub gsappend (f$) static  'not used in MT5PLOT

'-- open F$ as file number 7 in append mode for output of the plot commands
  shared kopen7
  open "a",7,f$
  kopen7=-1
end sub
'-----------------------------------------------------------------------------

sub gsclose (f$) static
  shared igraph, kopen7
   if igraph<2 then
     print #7, CHR$(26)  'write EOF
   elseif igraph=2 then
     print #7, "SK"
     print #7, "showpage"
   elseif igraph>2 then
     print #7,";":print #7,"PU;IP;PA1979,0;SP0;"  '19 bytes + CR/LF
     if igraph=3 then print #7,chr$(27);"%0A";chr$(12)  '7 bytes incl formfeed & CR/LF
   end if
  close 7 :  kopen7=0
end sub
'-----------------------------------------------------------------------------

sub gsorig (x,y) static
'-- change the origin by shifting it by X and Y inches

shared gscale, irotate, igraph
  g=gscale/25.4

  if irotate then
     xx=y*g : yy=-x*g
  else
     xx=x*g : yy=y*g
  end if
  if igraph=0 then
    print #7, "TR"; xx; yy
  elseif igraph=1 then
    print #7,  xx; yy; 1
  elseif igraph=2 then
    print #7,  72*xx; 72*yy; "T"
  else   'if igraph>2 then
    call hpcmd(xx,yy,1,0)
  end if
end sub
'-----------------------------------------------------------------------------

sub gsdraw (x,y,kcolor) static

'-- draw a line from the current point to X, Y in color KCOLOR on the screen

  shared gscale, irotate, igraph,kgrcol()
  g=gscale/25.4

  call scline1 (fnisx(x),fnisy(y),kgrcol(kcolor),0)

  if irotate then
     xx=y*g : yy=10-x*g
  else
     xx=x*g : yy=y*g
  end if

  if igraph=0 then
    print #7,"PA"; xx; yy
  elseif igraph=1 then
    write #7,  xx, yy, 2
  elseif igraph=2 then
    print #7,  72*xx; 72*yy; "D"
  else  'if igraph>2 then
    call hpcmd(xx,yy,2,kcolor)
  end if

end sub
'-----------------------------------------------------------------------------
sub gsmove (x,y,kcolor) static
'-- move from the current point to X, Y

  shared gscale, irotate, igraph, kgrcol()
  g=gscale/25.4
  call scpoint (fnisx(x),fnisy(y),kgrcol(kcolor))

  if irotate then
     xx=y*g : yy=10-x*g
  else
     xx=x*g : yy=y*g
  end if

  if igraph=0 then
    print #7, "MA"; xx; yy
  elseif igraph=1 then
    write #7,  xx, yy, 3
  elseif igraph=2 then
    print #7,  72*xx; 72*yy; "M"
  else 'if igraph>2 then
    call hpcmd(xx,yy,3,kcolor)
  end if

end sub

'-----------------------------------------------------------------------------
sub gstext (n, x, y, t$, an, sz) static
'--- write text T$ at X, Y in character set N
'-- letters are plotted on screen

 shared gscale, irotate, letset, igraph

  dq$=chr$(34) 'double quote
  g=gscale/25.4
  rpd=0.0174533

'---- plot on screen
    ix=fnisx(x):iy=fnisy(y):size=fnisx(sz):ian=an:k=0
    call scletter(ix,iy,t$,ian,size,k)

'-- check if letter set already on, otherwise change it
  if letset<>n then
    if n>0 then a$=str$(n) : l$="set"+mid$(a$,2,len(a$)-1)+".sym"
    if n=0 then l$= "centered.sym"
    if igraph=0 then
      print #7, "SS ";dq$+l$+dq$
    elseif igraph=1 then
      write #7, 0,0,5,l$
    end if
    letset=n
  end if

  if irotate then
     ang=an-90 : xx=y*g : yy=10-x*g : size=sz*g
  else
     ang=an : xx=x*g : yy=y*g : size=sz*g
  end if

if igraph=0 then
  print #7,"PS"; xx; yy; size; ang; dq$+t$+dq$

elseif igraph=1 then
  write #7, xx, yy, 4, size, ang, t$

elseif igraph=2 then
  size=size*108 ' ( = 72 * 1.5 )
  if n=1 then print #7, size;" CF"
  if n=7 then print #7, size;" TRF"

  print #7, xx*72;yy*72; "M"
  if ang <> 0 then print #7, ang;" RO"
  t2$=t$: call clnpar(t2$)  '*** added 3-10-91
  print #7, "("; t2$ ; ") S"  '*** revised 3-10-91
  if ang <> 0 then print #7, -ang;" RO"

else  'if igraph>2 then
  call hptext(xx,yy,size,ang,t$)
end if
end sub

sub hptext(xx,yy,size,ang,t$)
shared xscale,yscale,xmagp,lcode
    '========================================================================
    '=          Plot string of characters, op-code 4                        =
    '========================================================================
  call hpcmd(xx,yy,3,kcolbak)
'-- size in centimeters
  print #7, ";SI";fnstrnum$(fnfix(yscale*size*.0015604*xmagp,3));",";fnstrnum$(fnfix(xscale*size*.0024878*xmagp,3));
'-- angle
  a=(90.+ang)*pii/180.
  irun=100*cos(a) : irise=100*sin(a)     'conv to %
  print #7, ";DR";fnstrnum$(irun);",";fnstrnum$(irise);
'-- text
  print #7, ";LB";t$;chr$(3);
  lcode=4
end sub

'-----------------------------------------------------------------------------
sub gsline ( x(1), y(1), n, lsyn, sincr, kcolor ) static
'----- plots a curve defined by X(),Y() with N values
' LSYN is the line type 0-dashed, 1-dotted, 2-solid
' SINCR is the increment for the dashes and dots
shared irotate, igraph

if n=0 then exit sub

if lsyn>=2 then                'solid line
  call gsmove (x(1),y(1),kcolor)
  for i=2 to n
    call gsdraw (x(i),y(i),kcolor)
    if igraph=2 and (i mod 500 =0) then print #7, "SN" : call gsmove (x(i),y(i),kcolor)
  next

else                           'dashed or dotted line

  dx0=x(1) : dy0=y(1) : call gsmove (dx0,dy0,kcolor)

  j=1 : ds=sincr : i=0
  if lsyn=1 then ln=0 else ln=-1

  while j<n    'loop through points
    j=j+1
    x=x(j) : y=y(j)

SAMEPT: s=fndist(x,dx0,y,dy0)

    if s<ds then
      ds=ds-s
      if ln then call gsdraw (x,y,kcolor)
      dx0=x : dy0=y

    else           'finish segment
      delx=ds*(x-dx0)/s : dely=ds*(y-dy0)/s
      dx0=dx0+delx : dy0=dy0+dely

      if lsyn=1 then                 'dotted line
        call gsdot (dx0,dy0) : i=i+1
        if igraph=2 and (i mod 500 =0) then print #7, "SN"
      else                           'dashed line
        if ln then call gsdraw (dx0,dy0,kcolor) else call gsmove (dx0,dy0,kcolor)
        i=i+1
        if igraph=2 and (i mod 500 =0) then print #7, "SN" : call gsmove (dx0,dy0,kcolor)
        ln=not ln
      end if

      if not ln then ds=0.7*sincr
      if ln or lsyn=1 then ds=sincr
      goto SAMEPT

    end if

  wend
end if
if igraph=2 then print #7, "SN"
end sub
'-----------------------------------------------------------------------------
sub gscircle (x, y, rd, dthet, kcol) static
' draw a circle

shared irotate, igraph, gscale, kscrx, kscry, prad, lcode, xmagp,xscale
  g=gscale/25.4
  rpd=0.0174533
  rad=rd*g

  if irotate then
     xx=y*g : yy=10-x*g
  else
     xx=x*g : yy=y*g
  end if

  if igraph<2 then
    xxx=xx+rad:yyy=yy
    if igraph=1 then write #7, xxx, yyy, 3  else print #7, "MA"; xxx; yyy
    for xangle=dthet to 360. step dthet
      xxx= xx + rad * cos(xangle*rpd)
      yyy= yy - rad * sin(xangle*rpd)
      if igraph=1 then write #7, xxx, yyy, 2 else print #7, "PA"; xxx; yyy
    next
  elseif igraph=2 then
    print #7, "SN"
    print #7, 72*xx; 72*yy; 72*rad;" 0 360 arc"
    print #7, "SN"
  elseif igraph>2 then
    call hpcmd(xx,yy,3,kcol)
    print #7,";CI";fnstrnum$(int(rad*xmagp*xscale));",";fnstrnum$(fnxlimit(int(dthet),1,30));
    lcode=0
  end if

end sub
'-----------------------------------------------------------------------------

sub gssym (kfill, x, y, ksymb, an, sz) static
'--- write symbol T$ at X, Y
'-- letters are not plotted on screen
'-- kfill is false for no fill, true for fill

  shared gscale, irotate, letset, igraph,lcode,xmagp,xscale

  dim n1(10), n2(10), n3(10), newch(14)
  dq$=chr$(34) 'double quote
  g=gscale/25.4
  rpd=0.0174533

  ksym=ksymb

  if irotate then
     ang=an-90 : xx=y*g : yy=10-x*g : size=sz*g
  else
     ang=an : xx=x*g : yy=y*g : size=sz*g
  end if

  if igraph<2 then

'-- check if letter set already on, otherwise change it
    if letset<>0 then
      if igraph=1 then write #7, 0,0,5, "centered.sym"_ 
        else print #7, "SS "; dq$+"centered.sym"+dq$
      letset=0
    end if

    if ksym<32 then ksym=ksym+32
    if kfill then kn=4 else kn=1

    for i=1 to kn
      s= size/kn * i
      if igraph=1 then write #7, xx, yy, 4, s, ang, chr$(ksym)_
        else print #7, "PS"; xx; yy; s; ang; dq$+chr$(ksym)+dq$
    next

  else  'if igraph>1 then

'  routine for symbols
'  KSYM is the symbol number
    n1(1)=0: n1(2)=0: n1(3)=45: n1(4)=0: n1(5)=0: n1(6)=0: n1(7)=60
    n2(1)=360: n2(2)=360: n2(3)=405: n2(4)=360: n2(5)=360: n2(6)=360: n2(7)=420
    n3(1)=30: n3(2)=90: n3(3)=90: n3(4)=60: n3(5)=45: n3(6)=120: n3(7)=120

    newch(1)=1 : newch(2)=12: newch(3)=5 : newch(4)=11: newch(5)=3 : newch(6)=10
    newch(7)=4 : newch(8)=6 : newch(9)=8 : newch(10)=2 : newch(11)=14
    newch(12)=9: newch(13)=7 : newch(14)=13

    if ksym>31 then ksym=ksym-32
    if ksym<0 then ksym=1
    if ksym>14 then ksym=14

    ksym=newch(ksym)

    ips=(igraph=2)
    rad=size/2

    if ips then print #7, "SN"
    if ips and ksym=1 then
      print #7, xx*72;yy*72;size*36;" 0 360 arc"
      nfill=kfill
    else

      nc1= ksym mod 2
      ksym = int( (ksym-1)/2 )+1

'-- only odd numbered symbols can be filled in
      if nc1=1 and kfill then nfill=ktrue else nfill=kfalse

      xl= xx + rad*sin(n1(ksym)*rpd)
      yl= yy - rad*cos(n1(ksym)*rpd)
      if ips and nfill then print #7, "NP"
      if nc1=1 then
         if ips then
           print #7, xl*72;yl*72;" M"
         else
           if nfill then 'shade wedge
             call hpcmd(xx,yy,3,kcolbak)
             print #7,";WG";fnstrnum$(int(rad*xmagp*xscale));",";fnstrnum$(n1(ksym));",360,";fnstrnum$(n3(ksym));
           else          'absolute arc
             call hpcmd(xl,yl,3,kcolbak)
             call hpcnvrt(xx,yy)
             print #7,";PD;AA";fnstrnum$(int(xx));",";fnstrnum$(int(yy));",360,";fnstrnum$(n3(ksym));
           end if
           lcode=0
           exit sub
         END IF
      end if
      if nc1=0 then nlast=n2(ksym)-1 else nlast=n2(ksym)

      for i=n1(ksym) to nlast step n3(ksym)
        if nc1=0 then
          if ips then print #7, xx*72;yy*72;" M" else call hpcmd(xx,yy,3,kcolbak)
        end if
        xl=rad*sin(i*rpd)
        yl=rad*cos(i*rpd)
        if ips then print #7, (xx+xl)*72;(yy-yl)*72;" D" else call hpcmd(xx+xl,yy-yl,2,kcolbak)
      next
    end if

    if ips then
      if nfill then
        print #7, "CP "
        print #7, "fill"
      end if
      print #7, "SN"
    end if
  end if
end sub

'-----------------------------------------------------------------------------
sub gsdot (x, y) static
'--- put dot at X, Y
   shared gscale, irotate, kscrx, kscry, igraph,fisx,fisy, scrorig,kgrcol()
   g=gscale/25.4

  call scpoint (fnisx(x),fnisy(y),kgrcol(kcolbak))

 if irotate then
     xx=y*g : yy=10-x*g
 else
     xx=x*g : yy=y*g
 end if

 if igraph=0 then
    print #7, "MA"; xx; yy : print #7, "PA"; xx; yy
 elseif igraph=1 then
    write #7, xx, yy, 3 : write #7, xx, yy, 2
 elseif igraph=2 then
    print #7, 72*xx; 72*yy; " DOT"
 else  'if igraph>2 then
    call hpcmd(xx,yy,3,kcolbak)
    call hpcmd(xx,yy,2,kcolbak)
 end if

end sub
'-----------------------------------------------------------------------------
sub gstroke static

'-- stroke a segment
  shared igraph
   if igraph=2 then print #7, "SK"
end sub

'-----------------------------------------------------------------------------

sub scletter(ix, iy, c$, iangin, size, isize) static
' --- put letters on screen, size in mm
'-- if isize is true then it returns the letter size in the ISIZE variable but
'   does not plot


shared syndir$, scratio, io$

sz=size/100. : ang=-iangin*0.0174533
ca=cos(ang) : sa=sin(ang)

x01=ix : y01=iy
x0=0 : y0=0

'-- opened in main routine MT5PLOT
'open syndir$+"m5letset.raf" as 4 len=2: field #4, 2 as io$

for i=1 to len(c$)
 j=asc(mid$(c$,i,1))

 if j<1 or j>122 then xx=25*sz : goto GL99
 if j=32 then xx=50*sz : goto GL99

  xx=0

km=36*3*(j-1)+1
get #4, km : nm=cvi(io$)

if nm=0 then xx=25*sz : goto GL99

 for k=1 to nm
   get #4 : lm=cvi(io$)
   get #4 : ixi=cvi(io$)
   get #4 : iyi=cvi(io$)

   xi=ixi*sz
   yi=iyi*sz

    x=xi : y=yi

    xl= x*ca + y*sa +x0
    yl=(-x*sa + y*ca)*scratio +y0

     x=x01+xl :y=y01-yl
     if xi>xx then xx=xi

   if not isize then
     if lm=1 then pset (x,y) else line-(x,y)
   end if

 next

GL99:  x0=x0+(xx+3*sz)*ca : y0=y0-(xx+3*sz)*sa*scratio
next
' if isize then isize = fndist(x01, x0, y01, y0)
'close 4
end sub
'-----------------------------------------------------------------------------
sub hpcnvrt(x1,y1) static
' converts x and y coordinates to proper HPGL coordinates
shared xscale,yscale,xmagp,xmax  'ymax not used
shared xoffhp,yoffhp

x1=x1+xoffhp
y1=y1+yoffhp
SWAP x1,y1      'correct to HP plotter orientation
x1=xmax+.8-x1   'full 11" sheet
x1=x1*xscale*xmagp
y1=(y1-.05)*yscale*xmagp  'make it fit
end sub
'-----------------------------------------------------------------------------
sub hpcmd(x1,y1,icode,kcolor) static
shared xscale,yscale,xoffhp,yoffhp,lcode,igraph

' command codes for HP7475a plotter
CONST mov1$="PA",trm$=";":up$="PU":dn$="PD":gpen$="SP1":ppen$="SP0"

lastcol=(igraph=4 and y1>y) or (igraph=5 and x1<x)
if not lastcol then call hpcolor(kcolor) 'get pen when closest to it
x=x1:y=y1
DO
  if icode<>lcode then
    if lcode<>4 then print #7,trm$;
    SELECT CASE icode
      CASE 1   'change origin
        xoffhp=x+xoffhp:yoffhp=y+yoffhp:lcode=4:exit do
      CASE 2   'move with pen down
        print #7,dn$;
      CASE 3   'move with pen up
        print #7,: print #7,up$;
      'CASE 4   'plot string of characters
      'CASE 5   'load symbol set
      CASE 6   'change scaling
        xscale=xscale*y:yscale=yscale*x:lcode=4:exit do
      CASE ELSE
        lcode=4:exit do
    END SELECT
  else
    print #7,",";
  end if
  call hpcnvrt(x,y)
  xs$=fnstrnum$(INT(x)):ys$=fnstrnum$(INT(y))
  if y>0 then print #7,xs$;",";ys$;
  lcode=icode
  exit do
LOOP
if lastcol then call hpcolor(kcolor)
x=x1:y=y1 'save as reference
end sub

sub hpcolor(kcolor)
shared khpcol,kpenhp,lcode
CONST trm$=";"
if kcolor>0 then  'change pens
  IF kpenhp then
    SELECT CASE kcolor
      CASE kcolbak
        kpen=1  'defaults to pen 1 if all the same
      CASE kcoldat
        kpen=2
      CASE kcolsyn
        kpen=3
      CASE ELSE
        kpen=1
    END SELECT
  ELSE
    kpen=1
  END IF
  if kpen<>khpcol then
    khpcol=kpen
    if lcode<>4 then print #7,trm$: else print #7,
    print #7,"PU;SP";fnstrnum$(kpen);:lcode=0 'revised to lift pen 1st - 4/10/1991
  end if
end if
end sub

'*** following added 3-10-91, modified 1/13/1994
sub clnpar(s$) 'to prevent parenthesis in string from interfering with PS
for i=len(s$) to 1 step -1
  if instr("\()%",mid$(s$,i,1))>0 then s$=left$(s$,i-1)+"\"+mid$(s$,i)
next
end sub
