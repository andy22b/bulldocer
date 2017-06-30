' main MT5PLOT (MT5.00)
' plotting routines for MT5
'-- MUST BE COMPILED WITH /AH OPTION FOR HUGE ARRAYS
DEFINT I-N: OPTION BASE 1

'following are shared with M5MTSUB module
COMMON /PLOT/ fisx,fisy,gscale,io$,irotate,khpcol,kopen7,_
  lcode,letset,prad,xmax,xscale,yscale

rem $INCLUDE: 'M5COMMON'
rem $INCLUDE: 'M5FUNCTN'

if not imt5 then call notalone

ntext=0 : kcol0=0
n=2*maxnes  'changed from 1000  11/10/93 PZ 

'-- initialize plot parameters
xmax=9.84:ymax=7.87:xscale=1:yscale=.9907:lcode=0:khpcol=0

dim xtmp(n), ytmp(n), letsym(nsta), phslead(7)
dim mech$(nsor), ks(nsta)

'---- convert plot coords (250x200mm) to undistorted screen coordinates
  if (kscry-200*scratio*kscrx/250)<0 then  'if screen too wide
     fisx=kscry/(scratio*200)  '           'pixels/(mm x dimension)
     fisy=kscry/200.                       'pixels/(mm y dimension)
  else
     fisx=kscrx/250.
     fisy=kscrx*scratio/250
  end if

  def fnisx(x) = cint(x*fisx)
  def fnisy(y) = cint(kscry-(y-scrorig)*fisy)
  def fnyn$(iyn) = mid$("YN",2+(iyn<>0),1)

'---- amp and time functions for seismograms
  def fnamp(a)=ybox + a*scl*seisamp/(-(pext=0.0)-(pext<>0.0)*pext)
  def fntim(t)=xbox + t*spc0

'----- format setup
'---positions of seismograms in mm
'-- 2 configurations can be used, IPLFORM selects the set
'-- set 1 is predefined, set 2 is for sweep
'--  the plot area is 250x200 mm

  pagewide=250 : pagehigh=200
  nforms=2 : nboxmax=40
 dim ipxpos(nboxmax,nforms), ipypos(nboxmax,nforms), nboxes(nforms)
 dim nbrow(nforms), nbcol(nforms)

'-- format 1 defined but can be modified by user
nboxes(1)=40 : boxhigh0(1)=20 : boxlen0(1)=40
data 125,130, 165,130, 165,110, 165, 90, 165, 70, 165, 50, 165, 30, 125, 30
data  85, 30,  45, 30,  45, 50,  45, 70,  45, 90,  45,110,  45,130,  85,130
data 125,150, 165,150, 205,150, 205,130, 205,110, 205, 90, 205, 70, 205, 50
data 205, 30, 205, 10, 165, 10, 125, 10,  85, 10,  45, 10,   5, 10,   5, 30
data   5, 50,   5, 70,   5, 90,   5,110,   5,130,   5,150,  45,150,  85,150
 for i=1 to nboxes(1) 
  read ipxpos(i,1), ipypos(i,1)
 next

'--- lead times for plotting seismograms of each type
 data 7.0, 7.0, 2.0, 7.0, 7.0, 2.0, 7.0
  for i=1 to 7 : read phslead(i) : next i

'-- graphics 
  plang$="PlotCall2+PlotCall<2PostScriptLaserJet3 Plot8x11  Plot11x17 "
  ltypes$="DASHEDDOTTEDSOLID "

'--- sizes of lettering in mm
biglet=4.0      'big letters in mm
smlet= 3.5      'small letters
titlsz=7.0      'size of title letters
letset=20       'current letter set
ldat=2          'data are solid lines

'---- file name initializing
if len(gs$)=0 or gs$="nul" then gs$="nul": kgsnul=ktrue

'---- get the current format parameters
if ifocrad0(1)=0 then gosub GETFORMS else gosub GETSCALES

if boxlen0(2)=0 then
   boxlen0(2)=(pagewide-10.)/10. : boxhigh0(2)=20
end if

ibmake=2 : gosub MKBOXES

'-- check for default plot positions
 for ista=1 to nsta
   if ispos(ista,1)>0 and ixfine(ista)=0 and iyfine(ista)=0 then
    ixfine(ista)=ipxpos(ispos(ista,1),1) : iyfine(ista)=ipypos(ispos(ista,1),1)
   end if
 next

gosub SETGSCALE

'-- open RAF file containing letters
open syndir$+"m5letset.raf" as 4 len=2: field #4, 2 as io$

'-- make gsave array to store graphics image for format 2
SELECT CASE kega
   CASE 1,3,5
     iplane=1    '# planes
   CASE ELSE
     iplane=4
END SELECT
kdely=(kscry+1)/nrowchar 'pixel height of characters
ibpp=1 'bits per pixel per plane: isgsave=kfalse
'igsave = (11+INT(((kscrx+1)*ibpp+7)/8)*iplane*kdely)/8.
'dim glnsave(igsave) as double
igsave = 1+(11+INT(((kscrx+1)*ibpp+7)/8.)*iplane*(kscry-kdely+1))/8. 
mksave = (fre(-1)/8. > igsave+2000)  'leave an extra 16000 bytes to play in
if mksave then
  dim gsave(igsave) as double
  if isfile((drive$+"screen.tmp"),(""),1)=0 then 
    open drive$+"screen.tmp" as 7 len = 8
    i=lof(7)-igsave*8.
    if fninlim(i,0,8) then
      print "loading saved graphics..."
      if igsave>&H2000 then 'more than a segment
        for i=1 to igsave
          get 7,,gsave(i)
        next
        close 7
      else
        close 7
        def seg=varseg(gsave(1))
        bload drive$+"screen.tmp",varptr(gsave(1))
        def seg
      end if
      scrorig=gsave(1) 'last screen orig
      isgsave=ktrue
    else
      close 7
      kill drive$+"screen.tmp"
    endif
  end if
'else
'  erase glnsave
end if

'------------ end setup

200  '--- Plot menu
call setscreen(0)
call scclear

if kopen7 then call gsclose(filopen$)
itimerow=2:itimecol=61
call scbox(1,13,ktrue,"PLOTTING MENU","Esc to return to Main Menu")
call scprint(2, 2,"KEY  OPTION")

call scprompt(4,3,"G - plot file"):call scprintv((0),47,"<"+gs$+">")
if kappend then 
  call scprintv((0),(0)," in APPEND mode ") 
else
  if not kgsnul then call scprintv((0),(0)," as NEW file")
end if
call scprompt(5,3,"C - Change plot positions for seismograms")
call scprompt(6,3,"E - Editing plot positions while plotting")
call scprintv((0),47,fnyesno$(kedit))
call scprompt(7,3,"F - plot language")
call scprintv((0),47,mid$(plang$,1+igraph*10,10))
call scprompt(8,3,"O - go to PLOT OPTION MENU")
call scprompt(9,3,"P - run Plot")
call scprompt(10,3,"W - Write parameter file")
call scprompt(11,3,"X/"):call scprompt(11,5,"Q - Shell/Quit")


if input1(14,2,"Enter KEY for OPTION ",as$,kfalse,ktrue,ktrue,"gcefopwxq")=27 then goto 9999
itimerow=0
if fniquit(as$) then call endcheck

locate 16,1
if as$="c" then 
  gosub CHPLOCS

elseif as$="e" then
  if iplform=1 then
    if kedit then i$="n" else i$="y"
    if input1(14,2,"Edit plot positions while plotting (y/n)",i$,kfalse,ktrue,ktrue,"yn")<>27 then kedit=fnkyes(i$)
  else
    print "Editing not allowed for Format 2. Press any key..."
    while inkey$="":wend
  end if

elseif as$="o" then 
  as$="":gosub OPMENU

elseif as$="p" then 
  goto 720

elseif as$="w" then 
  gosub M5WPARM

elseif as$="x" then
  locate 14,1
  call shellcom

elseif as$="f" then 
  call scprompt(14,2,"0-PLOTCALL (Golden Software, Version 2 and higher)")
  call scprompt(15,2,"1-PLOTCALL (Golden Software, before Version 2)")
  call scprompt(16,2,"2-PostScript")
  call scprompt(17,2,"3-HP LaserJet3") 
  call scprompt(18,2,"4-HP Plotter,8.5"+chr$(34)+"x11"+chr$(34)) 
  call scprompt(19,2,"5-HP Plotter,11"+chr$(34)+"x17"+chr$(34)) 
  d$=fnstrnum$(igraph)
  if input1(21,2,"Plot language ",d$,kfalse,ktrue,kfalse,"012345")<>27 then 
    igraph=val(d$)
    if igraph>3 then 
      call scprint(22,2,"Separate pens are: text=pen 1, data=pen2, synthetics=pen3")
      kpenhp=igetyesno(23,2,"Use separate pens for text, data, and synthetics?","Y",ktrue)
    else
      kpenhp=kfalse
    end if
  end if

elseif as$="g" then
  kapp=kappend
  DO
    notdone=kfalse
    if gs$="nul" then d$=fhdr$ else d$=fnsubext$(gs$)
    call scprintc(14,1,"Enter plot file prefix (I/O file header is "+fhdr$+")")
    if inputn(15,2,"--",d$,ktrue,ktrue,kfalse,"a",40)<>27 then
      igsfn=0 'start out with a new .PHM plot
      d$=lcase$(fnsubext$(d$)): kgsnul=(len(d$)=0 or d$="nul")
      if kgsnul then
        gs$="nul"
      else
        print:print
        i=0
        if isfile(d$+".psh",f$,1)=0 then print f$;" - (format 1, single page)":i=1
        if isfile(d$+".p",f$,1)=0 then print f$;"   - (format 1, page 1 of 2)":i=1
        if isfile(d$+".sh",f$,1)=0 then print f$;"  - (format 1, page 2 of 2)":i=1
        if isfile(d$+".phm",f$,1)=0 then print f$;" - (format 2)":i=2 
        if i>0 then
          notdone= not igetyesno((0),1,"CAUTION - existing plot files listed above can be overwritten. OK?","",ktrue) 
        end if
        kapp=kfalse
        if not notdone then 
          if iplform=2 then
            if i=2 then kapp=igetyesno(22,2,f$+" - Append mode?","N",ktrue)
          end if
          gs$=d$
        end if
      end if
    end if
  LOOP while notdone
  kappend=kapp
end if
goto 200

OPMENU:
'--- plot controls menu
call scclear
call scbox(1,21,ktrue,"PLOT OPTION MENU","Esc to return to Plotting Menu")
call scprint(2, 2,"KEY  OPTION                     VALUE")
call scprint(2,43,"KEY  OPTION                   VALUE")
for i=2 to 11 : locate i,42: print chr$(179); : next

call scprompt(4,3,"A - Amplitudes normalized"):call scprintv(4,35,fnyesno$(normamp))
call scprompt(5,3,"B - letters in focal spheres"):call scprintv(5,35,fnyesno$(kletter))
call scprompt(6,3,"D - synthetics line"):call scprintv(6,34,mid$(ltypes$,1+(lsyn)*6,6))
call scprompt(7,3,"I - dash/dot Increment in mm"):call scprintv(7,35,str$(sincr))
call scprompt(8,3,"E - FOCAL SPHERE CONTROL MENU")
call scprompt(9,3,"J - write page label"):call scprintv(9,35,fnyesno$(kplabel))
call scprompt(10,3,"N - plot position format"):call scprintv(10,35,str$(iplform))
call scprompt(11,3,"P - run Plot")

call scprompt(4,44,"O - synthetics & data Offset"):call scprintv(4,74,fnyesno$(iploff))
call scprompt(5,44,"T - plot Time & amp scales"):call scprintv(5,74,fnyesno$(ktsc))
call scprompt(6,44,"U - S/D/R/H/M control"):call scprintv(6,74,str$(kpmech))
call scprompt(7,44,"V - View plot positions")
call scprompt(8,44,"W - Wait between plots"):call scprintv(8,74,fnyesno$(kdly))
call scprompt(9,44,"1 - plot P waves"):call scprintv(9,74,fnyesno$(kplot(1)))
call scprompt(10,44,"2 - plot SH waves"):call scprintv(10,74,fnyesno$(kplot(2)))
call scprompt(11,44,"3 - new page for SH"):call scprintv(11,74,fnyesno$(i2page))

call scprompt(13, 3,"L - title"):call scprintv(13,14,"<"+title$+">")
print "  position=";:call scprintv(13,(0),str$(labpos))
call scprompt(14, 3,"F - position of time Function "):call scprintv(14,(0),str$(isposit))
if ktsc and not (iplform=2) then print "  and amplitude scale ";:call scprintv(14,(0),str$(iampos))

call scprompt(15,3,"M - Magnification factors:")
locate 16,6
for i=1 to nphtype
  if kphflag(i) then 
    print using fnaform$(4);fnkphlab$(i);
    call scprintvfu(16,(0),"###.# ",phscl(i))
  end if
next i
call scprompt(17,3,"S - Seismogram lengths to plot:")
locate 18,6
for i=1 to nphtype 
  if kphflag(i) then 
    print using fnaform$(4);fnkphlab$(i);
    call scprintvfu(18,(0),"###.# ",phtime(i))
  end if
next i

call scprompt(19,3,"H - time scale (mm/min) : LP/BB")
call scprintv(19,(0),str$(fnfix(60.*spclp*gscale,1)))
if ispflag then print " SP";:call scprintv(19,(0),str$(fnfix(60.*spcsp*gscale,1)))

if input1(22,1,"Enter KEY for OPTION ",as$,kfalse,ktrue,ktrue,"abdiejnpotuvw123lfmshq")=27 then return
if fniquit(as$) then return

if as$="1" then kplot(1)=igetyesno(22,1,"Plot P waves",fnyn$(not kplot(1)),ktrue)
if as$="2" then kplot(2)=igetyesno(22,1,"Plot SH waves",fnyn$(not kplot(2)),ktrue)

if as$="3" then 
  if iplform=2 then 
    print " Plot Format 2 requires single page plot. Press any key.."
    while inkey$="":wend
    i2page=kfalse 
  else
    i=igetyesno(22,1,"SH waves on new page",fnyn$(not i2page),ktrue) 
    if i<>i2page then i2page=i : gosub TIMESCALES
  end if

elseif as$="a" then
  if invflag then
    normamp=igetyesno(22,1,"Normalizing amplitudes",fnyn$(not normamp),ktrue)          
  else
    print "Matching shape requires amplitude normalization. Press any key.."
    while inkey$="":wend
    normamp=ktrue
  end if

elseif as$="b" then 
  kletter=igetyesno(22,1,"Write station letters in focal spheres",fnyn$(not kletter),ktrue) 

elseif as$="d" then
  d$=fnstrnum$(lsyn)
  if input1(22,1,"Synthetics line type: (0)dashed (1)dotted (2)solid",d$,kfalse,ktrue,ktrue,"012")<>27 then lsyn=val(d$)

elseif as$="e" then 
  gosub FSMENU

elseif as$="h" then 
  gosub TIMESCALES

elseif as$="i" then
  d$=fnstrnum$(sincr)
  if inputn(22,1,"Enter increment in mm for dashes or dots in seismograms",d$,ktrue,ktrue,ktrue,"f",5)<>27 then sincr=fnamax(.1,val(d$))

elseif as$="j" then 
  kplabel=igetyesno(22,1,"Write long page label",fnyn$(not kplabel),ktrue)

elseif as$="l" then
  input "Position for title (0 to center on top)"; l
  line input " Title (type 'none' for no title) "; t$
  if len(t$)>0 then title$=t$
  if t$="none" then title$=""
  if fninlim(l,1,npos) then labpos=l else labpos=0

elseif as$="m" then
  print "Scale factor for  (<0 to skip)"
  for i=1 to nphtype
    if kphflag(i) then 
      print fnkphlab$(i);
      input s
      if s>0 then phscl(i)=s
    end if
  next i

elseif as$="n" then 
  d$=fnstrnum$(iplform):i=iplform
  if input1((0),2,"Format for seismogram positions (1-"+fnstrnum$(nforms)+")",d$,kfalse,ktrue,kfalse,numlist$(1,nforms))<>27 then
    iplform=val(d$):print
    call getplotnm(gs$,gsf$,gs$) 
    if iplform=2 then 
      if i<>iplform then  'if this is a mode change
        if isfile(gsf$,f$,1)=0 then kappend=igetyesno((0),2,gsf$+" - Append mode?","Y",ktrue)
      end if
      i2page=kfalse: kedit=kfalse
      input "Number of plot positions across the page "; n
      if n>0 then boxlen0(2)=(pagewide-10)/n : boxhigh0(2)=20 :ibmake=2 : gosub MKBOXES
      input " Vertical offset between lines of seismograms (in mm; 30 mm default) ";y
      if y>0 then yappoff=y
      input " Station codes printed at top(0) or at left(1) of seismograms ";l
      labcode=not (-l)
    else
      kappend=kfalse
    end if
    gosub TIMESCALES
  end if

elseif as$="o" then 
  iploff=igetyesno(22,1,"Seismograms offset (otherwise overlapping)",fnyn$(not iploff),ktrue)

elseif as$="p" then 
  goto 720

elseif as$="f" then
380  input "New position for time function (0=default, <0 for none) "; isp
  if fninlim(isp,0,npos) then 
    if isp=0 then 
      isposit=nboxes(iplform)
    else 
      isposit=isp
    end if
  elseif isp<0 then
    isposit=0
  end if

381  input "New position for amplitude scale (0 default) "; i
  if fninlim(i,0,nboxes(iplform)) then iampos=i else _
   print "position does not exist" : goto 381

elseif as$="s" then
  print "Time length to plot for    ";
  print crin$
  for i=1 to nphtype
    if kphflag(i) then 
      print fnkphlab$(i)+" ";
      phtime(i)=fncrin(phtime(i))
    end if
  next i
  gosub TIMESCALES

elseif as$="t" then 
  ktsc=igetyesno(22,1,"Plot time and amplitude scales",fnyn$(not ktsc),ktrue)

elseif as$="u" then
  print "Position of mechanism label STRIKE/DIP/RAKE/DEPTH/MOMENT"
  input " 0-don't plot, 1-center under title, 2-center over focal sphere";kpmech

elseif as$="v" then 
  gosub PLOTBOX

elseif as$="w" then 
  kdly=igetyesno(22,1,"Pause between plots",fnyn$(not kdly),ktrue)

end if

gosub SETGSCALE

goto OPMENU


720 '----------- start plot --------------------------------------
cls
  if itf=0 then print "Synthetics not yet calculated. Choose Iterate in Main Menu" : d$=input$(1) : goto 200

'-- check that some are to be plotted
 npl=0
 for ista=1 to nsta
   if iflag(ista) and (ispos(ista,iplform)>0) and (not nodat(ista)) then npl=npl+1
 next
 if npl=0 then print "No positions of seismograms entered. Select C in Plotting Menu" : d$=input$(1) : goto 200

gosub SETGSCALE
  if ilpflag and spclp=0 then gosub TIMESCALES
  if ispflag and spcsp=0 then gosub TIMESCALES

  gsrot=0
  if not (iplform=2) then kappend=kfalse else iampos=0
  
 datext=0      'find max of all data 
 for ista=1 to nsta
   if iflag(ista) and (ispos(ista,iplform)>0) and (not nodat(ista)) then
    datext=fnamax(datext,abs(dext(ista)*gcor(ista)))
   end if
 next

 if datext=0 then 
 for ista=1 to nsta
   if iflag(ista) and (not nodat(ista)) then
    datext=fnamax(datext,abs(dext(ista)*gcor(ista)))
   end if
 next
 end if
 pext=datext
  
'-- open file
 if not (iplform=2) and not i2page then
      gsf$=gs$+".psh"
      call gsopen(gsf$) : filopen$=gsf$
      call gsinit(gsf$)
      call gsorig(xoff,yoff) 
      scrorig=0.0
 end if

 if (iplform=2) then
   gsfext$=".phm" 
   if kappend and scrorig+2*yappoff>200 then 'need new file extension for next plot file
     kappend=kfalse
     igsfn=igsfn+1
   end if
   if igsfn>0 then gsfn$=fnstrnum$(igsfn):gsfext$=left$(gsfext$,4-len(gsfn$))+gsfn$
   gsf$=gs$+gsfext$
   appfil$=drive$+"plot.tmp"
   call gsopen (appfil$) : filopen$=appfil$
   if not kappend then 
     call gsinit (appfil$) 
     call gsorig (xoff,yoff) 
     scrorig=0.0
'*** don't reuse the screen image
     isgsave=kfalse
   else
     if mksave then scrorig=gsave(1)+yappoff
     call gsorig (0.0,-yappoff)
   end if
 end if

'----- LOOP through P and SH

for iphase=1 to 2
    if iphase=1 or not (iplform=2) then 
      call setscreen(kega)
'*** display the saved screen image
      if iplform=2 and isgsave then put (0,0),gsave(2) ,PSET
    end if
    if not kplot(iphase) then 850

'-- outline plot area with box
    call scline2(fnisx(0), fnisy(scrorig), fnisx(pagewide), fnisy(pagehigh+scrorig), kgrcol(kcolbak), ktrue)

'-- focal sphere position and radius
    focx0=ifocx0(iphase) : focy0=ifocy0(iphase) : focrad0=ifocrad0(iphase)
    symsiz=fnamin(smlet,focrad0*.05)

'-- see if any seismograms have non-zero scale
       knonzero=kfalse
       for i=1 to nphtype
          if iphase=1 and ipwave(i) and phscl(i)>0 then knonzero=ktrue
          if iphase=2 and (not ipwave(i)) and phscl(i)>0 then knonzero=ktrue
       next
 if (not knonzero) then 850

'-- make file names and open
 if iphase=1 then
     if (not kgsnul and i2page) then gsf$=gs$+".p"
 else
     if (not kgsnul and i2page) then gsf$=gs$+".sh" 
     if not (iplform=2) and (not i2page) then call gsorig(0,-(yoff-10)) 
'     if not (iplform=2) and (not i2page) then call gsorig(xoff,-(yoff-10)) 
 end if

 if i2page then 
   call gsopen(gsf$) : filopen$=gsf$
   call gsinit(gsf$)
   call gsorig(xoff,yoff) 
 end if

'--- write title and MECH$ (= STRIKE/DIP/RAKE/DEPTH/MOMENT)
 if iphase=1 or (iphase=2 and not kplot(1)) then
   gosub PLOTLABL
   if kplabel and not (iplform=2) then
        xlaboff=0.5 : ylaboff=190
        call gstext( 1, xlaboff, ylaboff, art$, gsrot-90., smlet) 
   end if

'--- write title
   if labpos>0 then 
     x0=ipxpos(labpos,iplform) '+ boxlen0(iplform)/2
     y0=ipypos(labpos,iplform) '+ boxhigh0(iplform)
     x0t=x0
   else
     x0=ifocx0(1)
     y0=190    'to center title over P 
     x0t=x0 - (titlsz * len(title$))/(2-(igraph=2))
   end if

   if len(title$)>0 and not (iplform=2) then 
     call gstext(7,x0t, y0, title$, gsrot, titlsz)
   end if

'--- write S/D/R/H/M label
   if fninlim(kpmech,1,2) then
    for i=1 to nsor

      if kpmech=1 and not (iplform=2) then
        s=0.7*titlsz
        tlen= s * len(mech$(i)) 
        x0t=x0 - tlen/(2-(igraph=2))
        y0=y0-1.2*titlsz
        call gstext(7,x0t, y0, mech$(i), gsrot, s) 
      elseif kpmech=2 then
        s=focrad0/3
        tlen= s * len(mech$(i)) 
        x0t=fnamax(ifocx0(1)-tlen/3,0)
        y0=ifocy0(1)+1.1*ifocrad0(1)+1.2*s*(nsor-i)
        call gstext(1, x0t, y0, mech$(i), gsrot, s) 
      end if

    next
   end if

 end if

 if kfoc(iphase) then 
 str=sorpar(1,isor)
    gosub FOCSPH   'make focal sphere
    prad=ifocrad0(iphase)
    if iphase=1 then 
      call PLOTPRAD(x0,y0,prad)
      if mtflag then call PLOTPMTRAD(x0,y0,prad)
    elseif iphase=2 then 
      call PLOTSRAD(x0,y0,prad)
      if mtflag then call PLOTSMTRAD(x0,y0,prad)
    end if
 end if

'-- plot time scale and time function
      if ktsc then gosub PLOTAMP : gosub PLOTSTF

for ktype = 0 to 1
' ktype=0 for data, =1 for synthetics

  ltrplot=64   'ASCII code for letters denoting stations
  kmove=5      'increment for moving

'+++++++++++ Start station loop ++++++++++++++
 for ista=1 to nsta
  if fninquit(kfalse) then 200

    kpos=ispos(ista,iplform)
     if kpos=0 then 780
    scl=phscl(kphtype(ista))
    if scl<=0 or iph(ista)<>iphase or not iflag(ista) then 780
    if not fnilong(ista) then spc0=spcsp else spc0=spclp
    dt = dt(ista)

'---- edit plot position
  if kedit and ktype=0 then

'-- position on focal sphere
  a=fnrad(azd(ista))
  call eaproj2(focx0,focy0,focrad0,ptoa(ista),a,x,y)
  call scline2(fnisx(x+smlet/2),fnisy(y+smlet/2),fnisx(x-smlet/2),fnisy(y-smlet/2),kgrcol(kcolbak),ktrue)
KED1:
     kp=kpos

'-- get x,y position
     px=-(iplform=2)*ipxpos(kp,iplform)-(iplform=1)*ixfine(ista) 
     py=-(iplform=2)*ipypos(kp,iplform)-(iplform=1)*iyfine(ista)

       slen=phtime(kphtype(ista))*spc0
       shigh=boxhigh0(iplform)

     ix=fnisx(px) : iy=fnisy(py) : ix2=fnisx(px+slen) : iy2=fnisy(py+shigh)
     call scline2(ix,iy,ix2,iy2,kgrcol(kcolbak),ktrue)

KED2: call sclocate(1,1) 
     print fnstanm$(ista);" ";fnkphlab$(kphtype(ista));" Az";azd(ista);" Position X:";px;" Y:";py;

    call scprint(nrowchar-2,1,"KEYS: 8-Up 2-Dn 4-Lft 6-Rgt (or arrows) F-fast S-slow C-chg (or +/-) Q-quit")
    call inkey(e$)

        if fniquit(e$) then goto 200
        if e$="p" then kpos=kp:goto KEDOUT

        if e$="c" then 
           call cline (nrowchar-2,1)
           input "New Position index (0 no change, <0 sets to 0) "; newpos
           call cline(2,1)
           if newpos=0 then goto KED3
           if newpos<0 then kpos=0 : ixfine(ista)=0: iyfine(ista)=0:goto 780
           if fninlim(newpos,1,nboxes(iplform)) then 
            kpos=newpos
            ixfine(ista)=ipxpos(newpos,1) : iyfine(ista)=ipypos(newpos,1)
           end if
            goto KED3
         end if

      if iplform=1 then
        if e$="s" then kmove=1 : goto KED2
        if e$="f" then kmove=5 : goto KED2
        if e$="r" or fnmovert(e$) then ixfine(ista)=ixfine(ista)+kmove:kpos=kp
        if e$="l" or fnmovelt(e$) then ixfine(ista)=ixfine(ista)-kmove:kpos=kp
        if e$="u" or fnmoveup(e$) then iyfine(ista)=iyfine(ista)+kmove:kpos=kp
        if e$="d" or fnmovedn(e$) then iyfine(ista)=iyfine(ista)-kmove:kpos=kp
        ixfine(ista)=fnxlimit(ixfine(ista),1,pagewide)
        iyfine(ista)=fnxlimit(iyfine(ista),1,pagehigh)
      end if

KED3: call scline2(ix,iy,ix2,iy2,kcol0,ktrue)   'erase box
      goto KED1

KEDOUT: call scline2(ix,iy,ix2,iy2,kcol0,ktrue)   'erase box
        ispos(ista,iplform)=kpos

end if

    ltr=0
    if ista>1 then
      sista$=mid$(fnstanm$(ista),1,4)
     for i=1 to ista-1         ' see if same station already plotted
       if iph(i)=iphase then
          si$=mid$(fnstanm$(i),1,4)
        if si$=sista$ then ltr=letsym(i)
        if del(i)=del(ista) and azd(i)=azd(ista) then ltr=letsym(i)
       end if
     next
    end if

    if ltr=0 then ltrplot=ltrplot+1 : ltr=ltrplot
    letsym(ista)=ltr

       gosub GETSEIS  'get seismograms and max values

         x0seis=-(iplform=2)*ipxpos(kpos,iplform)-(iplform=1)*ixfine(ista)  
         y0seis=-(iplform=2)*ipypos(kpos,iplform)-(iplform=1)*iyfine(ista)
       if kletter then lbl$=chr$(ltr) else lbl$=""

       gosub PLOTSEIS   '--- plot seismogram

780  next ista
 '++++++++++ END STATION LOOP +++++++++++++++++

   if kdly then 
      if fninquit(ktrue) then 200
   end if

   next ktype

850 if i2page and kopen7 then call gsclose(filopen$)
    next iphase

999 if kopen7 then call gsclose(filopen$)
    if (iplform=2) and (not kgsnul) then
'      if mksave then GET (0,kscry-(kdely-1)) - (kscrx,kscry),glnsave
      if kappend then d$="Append" else d$="Save "
      i=igetyesno(nrowchar,40,d$+" this plot","Y",ktrue)
'      if mksave then PUT (0,kscry-(kdely-1)),glnsave, PSET

      if i then
'*** save the screen image
        if iplform=2 and mksave then 
          GET (0,0)-(kscrx,kscry-kdely),gsave(2)
          gsave(1)=scrorig
          isgsave=ktrue
        end if
        if kappend then
          if mksave then 'already written to
            if igraph<>1  then   '-- postscript or HPGL
              iflnum=FREEFILE
              open "b",#iflnum,gsf$,1
              if igraph=2 then 
                i=12 
              elseif igraph=3 then
                i=30
              else  'if igraph>3 then
                i=23
              end if
              if lof(iflnum)>i then
                seek #iflnum,lof(iflnum)-i+1  'go to start of closing cmd 
                a$=string$(i-2,32)
                put #iflnum,,a$        'write over it
              end if
              close #iflnum
            end if
          endif
          shell "copy "+gsf$+"+"+appfil$+ "  >nul"
        else
          shell "copy "+appfil$+" "+gsf$+" >nul"
          kappend=ktrue 
        end if
      else
        if mksave and scrorig>0 then scrorig=scrorig-yappoff
      end if
    end if
    goto 200
9999 close 4
    if mksave and isgsave then 
      print "saving graphics..."
      if igsave>&h2000 then 
        open drive$+"screen.tmp" as 7 len=8
        for i=1 to igsave
          put 7,,gsave(i)
        next
        close 7
      else
        def seg=varseg(gsave(1))
        bsave drive$+"screen.tmp",varptr(gsave(1)),igsave*8.
        def seg
      end if
    elseif isfile((drive$+"screen.tmp"),(""),1)=0 then
      kill drive$+"screen.tmp"
    end if
    chain SYNDIR$+"MT5MEN1"
    end
'-----------------------------------------------------------------------------
GETSEIS:  ' get data or syn for station ista, return SS, DEXT
  if ktype=0 and nodat(ista) then return
  t0=phtime(kphtype(ista))

  if ktype=1 then nin=fnamin(nes(ista),t0/dt(ista)+1) : jo=0   'synthetics
  if ktype=0 then nin=fnamin(ipf2(ista)+nes(ista),ndat(ista)) : jo=ipf2(ista)-1 'data

  call m5synio (1,ktype,ista,nin,0,ss(),dext,asq)

'-- get scaling normalization factor PEXT and scale seismograms
   if not normamp then pext=datext
   if ktype=0 then for i=1 to nin : ss(i)=ss(i)*gcor(ista) : next
   if normamp and ktype=1 then pext=synrms(ista)
   if normamp and ktype=0 then pext=datrms(ista)*gcor(ista)
   if pext=0 then pext=1
return

'-----------------------------------------------------------------------------
PLOTAMP: 
'--- plot amplitude scales
 if (normamp) or (iplform=2) then return

 nscl=0
 for iphtype=1 to nphtype
   if kphflag(iphtype) then
       if (not ipwave(iphtype)) and (iphase=2) then goto AMP1
       if (ipwave(iphtype)) and (iphase=1) then goto AMP1
       goto 2010

AMP1:
    nscl=nscl+1
    scl=phscl(iphtype) : if scl=0 then scl=1
 '-- A1CM= true seismogram amplitude to produce 10 mm amplitude on plot
     a1cm= 10*pext/(scl*seisamp) : xmicron=a1cm*10000.

2003 if xmicron>=10 then 
2004         y$=str$(xmicron) : k=instr(y$,".")
2005         if k>2 then amp$=mid$(y$,2,k-2) else amp$=y$
2006     else 
2007        print xmicron
         y$=str$(xmicron) : k=instr(y$,".")-1 : n=fnamin(len(y$)-k,5)
             amp$=mid$(y$,k,n)
         end if
2008    amp$=amp$+fnphmark$(iphtype)

'------------- amplitude axis
  if iampos=0 then 
     xbox= ifocx0(iphase)-0.75*ifocrad0(iphase)-7*nscl  : ybox= ifocy0(iphase)-0.5*ifocrad0(iphase)-10
  else 
     xbox=ipxpos(iampos,iplform)+3 + 7*nscl : ybox=ipypos(iampos,iplform)+.75*seisamp
  end if

  spc0=spclp
  call gsmove (xbox+1,fnamp(0),kcolbak)
  call gsdraw (xbox,fnamp(0),kcolbak)
  call gsdraw (xbox,fnamp(a1cm),kcolbak)
  call gsdraw (xbox+1,fnamp(a1cm),kcolbak)
  call gstext (1,xbox-1,fnamp(a1cm/2)-len(amp$)*smlet/2.7,amp$,gsrot+90,smlet)

  end if

2010 next iphtype
  return
'-----------------------------------------------------------------------------
PLOTSTF:
'------- plot time scale and source time function
2020 if isposit=0 or (iphase=2 and kplot(1)) then return   'no time axis for SH

  ayoff=7.0
  xbox=ipxpos(isposit,iplform)+7
  ybox=ipypos(isposit,iplform)+5 - ayoff
  gosub MAXTIMES

if ilpflag and not kappend then
    ybox=ybox+ayoff
    if ispflag then l$="s LP" else l$="s"
    in1=10
    in01=int((boxlen0(iplform)-len(l$)*smlet)/spclp/in1) * in1
    in02=int(timelp/in1)*in1
    in0=fnamin(in01,in02)
    in0=fnamax(30,in0)
    spc0=spclp
    gosub PLOTAXIS
end if

if ispflag and not kappend then 
    ybox=ybox+ayoff
    if ilpflag then l$="s SP" else l$="s"
    in1=2
    in01=int((boxlen0(iplform)-len(l$)*smlet)/spcsp/in1) * in1
    in02=int(timesp/in1)*in1
    in0=fnamin(in01,in02)
    in0=fnamax(8,in0)
    spc0=spcsp
    gosub PLOTAXIS
end if

'--- plot source time function on separate axis
'--- find max of stf's and total moment

  ybox=ybox+ayoff
  totm=0 : taumax=0 : timax=3.0
  for isor=1 to nsor
    for i=1 to ntau(isor) : taumax= fnamax(tam(isor)*atau(i,isor),taumax) :next
     timax=fnamax(timax,sorpar(5,isor)+tau*(1+ntau(isor)))
     totm = totm + tam(isor)
  next isor

       if timax<6 then in1=1 else in1=2
       if timax>=20 then in1=5
       in0= int(timax/in1 +1) * in1
       spc0=fnamax(1.6*smlet,boxlen0(iplform)-4*smlet)/in0
       l$="s STF"
    gosub PLOTAXIS

  for isor=1 to nsor
    n=1
    xtmp(1)=fntim(sorpar(5,isor))
    ytmp(1)=ybox
'-- max amplitude of STF is 10mm
    for i=1 to ntau(isor)
      n=n+1
      xtmp(n)=fntim(i*tau0 + sorpar(5,isor))
      ytmp(n)=ybox+tam(isor)*atau(i,isor)*10/taumax
    next
      n=n+1
      xtmp(n)=fntim((ntau(isor)+1)*tau0 + sorpar(5,isor))
      ytmp(n)=ybox

    call linetype(isor,kline,pincr)
    call gsline (xtmp(), ytmp(), n, kline, pincr, kcolbak)
  next isor

 return
'-----------------------------------------------------------------------------
PLOTAXIS:
'--- plot the time axis
  lb2$=mid$(str$(in0),2)+l$
  call gstext (1,fntim(0)-smlet/2,ybox-5,"0",gsrot,smlet)
  call gstext (1,fntim(in0)-smlet,ybox-5,lb2$,gsrot,smlet)
  call gsmove (xbox,ybox,kcolbak) 

  for i= 0 to in0 step in1   'time scale tics every IN1 secs
    call gsdraw (fntim(i),ybox,kcolbak)
    call gsdraw (fntim(i),ybox-1.5,kcolbak)
    call gsdraw (fntim(i),ybox,kcolbak)
  next
return

'-----------------------------------------------------------------------------
PLOTSEIS:
4000 '--- plot one seismogram  IS at x0seis,y0seis,seisamp

   if ktype=0 and iploff then y0seis=y0seis+1.5*seisamp   'if offset
   if ktype=1 and iploff then y0seis=y0seis+ .5*seisamp
   if (not iploff) then y0seis=y0seis+seisamp    'overlapping

   stat$=fnstanm$(ista) : st$=ucase$(stat$)

   lbl$=lbl$+fnphmark$(kphtype(ista))

'--- print station name and label
 if ktype=0 then
     if kfoc(iphase) then gosub FOCLET

'--- write station name and label
 if (not iplform=2) or (iplform=2 and not labcode) then
   lensta=len(st$) 
   for k=1 to lensta
     a$=mid$(st$,k,1)
     ys=y0seis+ seisamp-biglet-(biglet+1.0)*(k-1)
     call gstext (1,x0seis,ys,a$,gsrot,biglet) 
   next
   lenlab=len(lbl$)
   for k=1 to lenlab
     a$=mid$(lbl$,k,1)
     ys=y0seis+ seisamp-smlet-(smlet+0.5)*(k-1)
     call gstext (1,x0seis+biglet+1,ys,a$,gsrot,smlet) 
   next

 elseif (not kappend) and iplform=2 and labcode then

     if iphase=1 then l$=st$+" P" else l$=st$+" SH"
     l$=l$+fnphmark$(kphtype(ista))
     ts=boxlen0(iplform)/10
     xs=x0seis +2*biglet : ys=y0seis + boxhigh0(iplform)
     call gstext (7, xs, ys, l$, gsrot, ts)

 end if

end if

4020 xbox=x0seis+biglet+smlet+1 : ybox=y0seis
   dt=dt(ista)
   seclead0=phslead(kphtype(ista))
   ipf=fnamax(1, ipf2(ista)-int(seclead0/dt))  ' first data point to plot
   ipoff = ipf2(ista)-ipf
   toff=ipoff*dt
   npls=phtime(kphtype(ista))/dt
   if ktype=1 then ndata=nes(ista) +ipoff                           'synthetics
   if ktype=0 then ndata=fnamin(ipf2(ista)+nes(ista),ndat(ista))-ipf+1  'data
   n = fnamin(npls,ndata)  'number of points to plot
   if ktype=0 and (not nodat(ista)) then

   '------------------- PLOT DATA
       for j = 1 to n
         j1=j+ipf-1 : xtmp(j)=fntim((j-1)*dt(ista)) : ytmp(j)=fnamp(ss(j1))
       next

     call gsline(xtmp(), ytmp(), n, ldat, sincr, kcoldat) 

     '--------------- tic marks for inversion window
       t=fntim(toff)
     call gsmove (t,ybox-1.5,kcolbak)
     call gsdraw (t,ybox+1.5,kcolbak)
       t=fntim(toff+dt(ista)*(ninv(ista)-1))
     call gsmove (t,ybox-1.5,kcolbak)
     call gsdraw (t,ybox+1.5,kcolbak)

   elseif ktype=1 then

     '---------------- plot SYNTHETICS

     n=n-ipoff
       for j=1 to n
         xtmp(j)=fntim(dt(ista)*(j+ipoff-1)) : ytmp(j)=fnamp(ss(j))
       next

     call gsline (xtmp(), ytmp(), n, lsyn, sincr, kcolsyn)

   end if
return

'-----------------------------------------------------------------------------
FOCSPH: '---------- plot equal area net and planes
 x0=ifocx0(iphase) : y0=ifocy0(iphase) : npl=0

  xbox=x0 : ybox=y0 : prad=focrad0
  call gsmove (xbox,ybox+prad-1,kcolbak)
  call gsdraw (xbox,ybox+prad,kcolbak)

'-- draw circle
    dthet=3.0
    call gscircle (x0,y0,prad,dthet,kcolbak) 
    call sccircle (fnisx(x0),fnisy(y0),fnisx(prad),kgrcol(kcolbak))

'-- make + in center
  call scline2(fnisx(x0-1),fnisy(y0),fnisx(x0+1),fnisy(y0),kgrcol(kcolbak),kfalse)
  call scline2(fnisx(x0),fnisy(y0-1),fnisx(x0),fnisy(y0+1),kgrcol(kcolbak),kfalse)
  call gssym (kfalse, xbox, ybox, 39, gsrot, 3.0)

 return

'=============================================================================
CHPLOCS:
 cls
    input " Editing P (0) or SH (1) positions "; i
      if not fninlim(i,0,1) then return
      if i=0 then kpw= ktrue else kpw=false
      if i=0 then kpw$= "P" else kpw$="SH"
    input " Editing position indices (0) or X,Y (1) positions "; i
      if not fninlim(i,0,1) then return
      if i=0 then kpind= ktrue else kpind=false
      if i=0 then kp$= " Indices" else kp$=" position X,Y "

    irow=1 : icol=1 : ip=iplform
    nrows=19 : nwide=26

401 call setscreen(0)
  call scbox(1,22,ktrue,kpw$+kp$+" SEISMOGRAM PLOT LOCATIONS","Press Esc to exit")
  l$="Sta Phase F  Az   P1   P2 "  
  if not kpind then l$="Sta Phase F  Az    X    Y "  
  jsta=0
  call scprint(2,2,l$+l$+l$)
  for i=2 to 21 : call scprint(i,nwide+1,chr$(179)): call scprint(i,2*nwide+1,chr$(179)): next

  for i=1 to nsta
      k=kfalse
      if kpw and fnpwave(i) then jsta=jsta+1 :k=ktrue
      if not kpw and not fnpwave(i) then jsta=jsta+1 : k=ktrue
      if k then
        ks(jsta)=i
        ir= (jsta-1) mod nrows
        ic= int((jsta-1)/nrows)+1
        call sclocate(ir+3,(ic-1)*nwide+2)
       print using fnaform$(5); fnstanm$(i); fnkphlab$(kphtype(i));
       print using fniform$(1); -iflag(i);
       print using fniform$(4); azd(i); 
      call setvarcol
      if kpind then
       print using fniform$(5); ispos(i,1); ispos(i,2); 
      else
       print using fniform$(5); ixfine(i); iyfine(i); 
      end if
      call settxtcol
     end if
  next 

'--- get new value
403 ic=17+(icol-1)*nwide
    call scprintc(irow+2, ic+(ip-1)*5,"->")

    call scprintc(23,1,"KEYS: 8-Up 2-Dn 4-Lft 6-Rgt (or arrows) C-chg (or +/-) V-view plot locations")
    call inkey(d$)

     call scprintv(irow+2, ic+(ip-1)*5, "  ")

   if fniquit(d$) or fniesc(d$) then return
   if d$="v" then gosub PLOTBOX : goto 401
   if d$="u" or fnmoveup(d$) then irow=irow-1 
   if d$="d" or fnmovedn(d$) then irow=irow+1
   if d$="l" or fnmovelt(d$) then ip=ip-1 : ip=fnamax(ip,1)
   if d$="r" or fnmovert(d$) then ip=ip+1 : ip=fnamin(ip,2)

   if d$="c" or fniplus(d$) or fnminus(d$) then
     mc=(d$="c") : mp=fniplus(d$) : mm=fnminus(d$)
     i= (icol-1)*nrows + irow

'-- get new value of item
      if mc then
       call cline(23,1)
        input "New value "; ivalin
      elseif mp or mm then
       if kpind then 
         if mp then ivalin = ispos(ks(i),ip) +1 
         if mm then ivalin = ispos(ks(i),ip) -1 
       else
        if ip=1 then 
         if mp then ivalin = ixfine(i) +1 
         if mm then ivalin = ixfine(i) -1 
        else 
         if mp then ivalin = iyfine(i) +1 
         if mm then ivalin = iyfine(i) -1 
        end if
       end if
      end if

'-- modifying indices
       if kpind then
        if not fninlim(ivalin,0,nboxes(ip)) then 403
        ispos(ks(i),ip)=ivalin 
        if ivalin=0 and ip=1 then ixfine(ks(i))=0 : iyfine(ks(i))=0
        if ivalin>0 and ip=1 then ixfine(ks(i))=ipxpos(ivalin,1) : iyfine(ks(i))=ipypos(ivalin,1)

'-- modifying X,Y values
       else
        if ip=1 then
         if not fninlim(ivalin,1,pagewide) then 403
         ixfine(i)=ivalin : ispos(ks(i),ip)=1
        elseif ip=2 then
         if not fninlim(ivalin,1,pagehigh) then 403
         iyfine(i)=ivalin: ispos(ks(i),ip)=1
        end if
       end if
      call setvarcol:locate irow+2,ic+(ip-1)*5
      print using fniform$(5);ivalin;
      call settxtcol
   end if

408    if irow<1 and icol>1 then irow=nrows : icol=icol-1
       irow=fnamax(irow,1)
       if irow>nrows then irow=1 : icol=icol+1
       if (icol-1)*nrows+irow > jsta then irow=1 : icol=1
     goto 403

'=============================================================================
PLOTBOX:
420 ' show plot locations on screen
  cls :n=1
  input "Show plot format template(0) or current seismogram positions(1) ";i
    kstaf=-i

  if kstaf then
   input "P(0) or S(1) seismogram positions shown"; n
   if not fninlim(n,0,1) then goto BOXEND
    kpw=(n=0): n=n+1
  end if

 input " Plot format to show: 1 or 2 (CR for current format) "; ip
    if not fninlim(ip,0,nforms) then goto BOXEND
    if ip=0 then ip=iplform

  call setscreen(kega)
  print " Format = ";ip

  if iplform=1 then
   if kfoc(n) then call sccircle (fnisx(ifocx0(n)),fnisy(ifocy0(n)),fnisx(ifocrad0(n)),kgrcol(kcolbak))
  else
   for n=1 to 2
    if kfoc(n) then call sccircle (fnisx(ifocx0(n)),fnisy(ifocy0(n)),fnisx(ifocrad0(n)),kgrcol(kcolbak))
   next
  end if

 if not kstaf then
  for i=1 to nboxes(ip)
    ix=fnisx(ipxpos(i,ip)) : iy=fnisy(ipypos(i,ip))
    ix2=fnisx(ipxpos(i,ip)+boxlen0(ip)) 
    iy2=fnisy(ipypos(i,ip)+boxhigh0(ip))
    call scline2(ix,iy,ix2,iy2,kgrcol(kcolbak),ktrue)

    i$=str$(i) :if i=isposit then i$=i$+" STF"
    if i=iampos and invflag then i$=i$+" AMP"   
    ix=fnisx(ipxpos(i,ip)+2) : iy=fnisy(ipypos(i,ip)+2) : sz=fnisx(5.):k=0
    call scletter(ix,iy,i$,0,sz,k)
  next
 else
  for ista=1 to nsta+2
   kpos=0
   if ista<=nsta then
      if (kpw and fnpwave(ista)) or (not kpw and not fnpwave(ista)) or (ip=2) then kpos=ispos(ista,ip)
      a$=fnstanm$(ista)
      if kpos>0 then ix=-(ip=2)*ipxpos(kpos,ip)-(ip=1)*ixfine(ista):iy=-(ip=2)*ipypos(kpos,ip)-(ip=1)*iyfine(ista)
   else
      if ista=nsta+1 and n=1 then kpos=isposit : a$="STF"
      if ista=nsta+2 then kpos=iampos : a$="AMP"
      if kpos>0 then ix=ipxpos(kpos,ip): iy=ipypos(kpos,ip)
   end if

    if kpos>0 then
     ix1=fnisx(ix) : iy1=fnisy(iy)
     ix2=fnisx(ix+boxlen0(ip)) : iy2=fnisy(iy+boxhigh0(ip))
     call scline2 (ix1,iy1,ix2,iy2,kgrcol(kcolbak),ktrue)
     ix1=fnisx(ix+2) : iy1=fnisy(iy+2) : sz=fnisx(5.):k=0
     call scletter(ix1,iy1,a$,0,sz,k)
    end if

  next ista
 end if
d$=input$(1)

BOXEND: call setscreen(0) : return
'=============================================================================
FOCLET:
430 '-- put letters at station positions in focal sphere
    a=fnrad(azd(ista))
    call eaproj2(focx0,focy0,focrad0,ptoa((ista)),a,x,y)
     if kletter and iplform=1 then 
       call gstext (1,x-smlet/2, y-smlet/2, chr$(ltr), gsrot,smlet)
     else 
       call gssym (kfalse, x, y, 33, gsrot, symsiz)
       call scline2(fnisx(x+symsiz/2),fnisy(y+symsiz/2),fnisx(x-symsiz/2),fnisy(y-symsiz/2),kgrcol(kcolbak),ktrue)
     end if
return
'=============================================================================
GETFORMS:
'--- get the format for plotting seismograms

 for ifoc=1 to 2
  if ifocpos(ifoc)=0 then 
   if iplform=1 then ifocx0(ifoc)=125 : ifocy0(ifoc)=90 : ifocrad0(ifoc)=30
   if iplform=2 then 
    ifocx0(1)=ipxpos(1,2)+boxlen0(2)/4 :ifocy0(1)=ipypos(1,2)+boxhigh0(2)/2 : ifocrad0(1)=boxlen0(2)/4.2
    ifocx0(2)=ipxpos(1,2)+3*boxlen0(2)/4 :ifocy0(2)=ipypos(1,2)+boxhigh0(2)/2 : ifocrad0(2)=boxlen0(2)/4.2
   end if
  end if
 next
 gosub GETSCALES
return
'=============================================================================
GETSCALES:
'-- get scaling factors

  seisamp=boxhigh0(iplform)/2

'-- length of seismogram box in mm
  if iplform=2 and labcode then
    boxlen=boxlen0(iplform)
  else 
    boxlen=boxlen0(iplform)-(biglet+smlet+1)
  end if

  npos=nboxes(iplform)
return
'=============================================================================
MKBOXES:
'-- make the boxes for format IBMAKE

   nbcol(ibmake)=int(pagewide/boxlen0(ibmake))
   nbrow(ibmake)=int((pagehigh-20)/boxhigh0(ibmake))
   if ibmake=2 then nbrow(ibmake)=1

   k=0 : nboxes(ibmake)=fnamin(nboxmax,nbcol(ibmake)*nbrow(ibmake))
   for i=1 to nbrow(ibmake)
    for j=1 to nbcol(ibmake)
     k=k+1
     if k<=nboxmax then
      ipxpos(k,ibmake)= 5+(j-1)*boxlen0(ibmake)
      ipypos(k,ibmake)= 180-i*boxhigh0(ibmake)
     end if
   next j,i
return
'=============================================================================
TIMESCALES:
     gosub GETSCALES
     gosub MAXTIMES
     cls : call scprint(3,3,"ADJUSTING PLOT TIME SCALES ")
     locate 5,1
452  if ilpflag then
      print
      print "Default LP time scale is";60.*boxlen/timelp*gscale;"mm/min":print
      print "Current LP time scale is";60.*spclp*gscale: print
        input "Enter new long-period (LP) time scale (0 to keep, <0 for default)"; s
        if s<0 then 
          spclp=boxlen/timelp
        elseif s>0 then
          spclp=s/60/gscale                ' mm per second
          if timelp > boxlen/spclp then
            print : print " This time scale is longer than default"
            input " Some plots may overlap with this scale. Do yo want to keep it (y/n)"; y$
            if not fnkyes(y$) then 452
          end if
        end if
     end if

453  if ispflag then
      print
      print "Default SP time scale is";60.*boxlen/timesp*gscale;"mm/min":print
      print "Current SP time scale is";60.*spcsp*gscale:print
      input "Enter new short-period (SP) time scale (0 to keep, <0 for default)"; s
        if s<0 then 
          spcsp=boxlen/timesp
        elseif s>0 then
          spcsp=s/60/gscale                ' mm per second
          if timesp > boxlen/spcsp then
            print : print " This time scale is longer than default"
            input " Some plots may overlap with this scale. Do yo want to keep it (y/n)"; y$
            if not fnkyes(y$) then 453
          end if
        end if
     end if

return
'=============================================================================
MAXTIMES:
'--- get maximum times for long and short-period data
      timelp=0 : timesp=0
      for i=1 to nphtype
        if iphlong(i) and kphflag(i) then timelp=fnamax(phtime(i),timelp)
        if (not iphlong(i)) and kphflag(i) then timesp=fnamax(phtime(i),timesp)
      next i
return
'=============================================================================
PLOTLABL:  '-----------  get labels
title1$=fnsubext$(df$)
art$=title1$ + "  " + date$+" "

for i=1 to nphtype
  if kphflag(i) then art$=art$ + fnkphlab$(i) + "x" + fnstrnum$(phscl(i))+" "
next i

for isor=1 to nsor
  mech$(isor)=""
  if nsor>1 then mech$(isor)=mech$(isor) + fnstrnum$(isor) +":"
  for i=1 to 4
    if i<4 then 
      if sorpar(1,isor)>-998 then
        x=fndeg(sorpar(i,isor)) 
        a$=fnstrnum$(cint(x))
      else
        a$="NA"
      end if
    else 
      x=sorpar(i,isor)+cthk
      if x>=1 then a$=fnstrnum$(cint(x))
      if x<1 then x=x*10 : a$="0."+fnstrnum$(cint(x))
    end if
    mech$(isor)=mech$(isor) + a$ +"/"
  next 
  logmom=int(fnalog10(tam(isor)))
  ximom=tam(isor)/10^logmom
  mom$=fnstrnum$(fnfix(ximom,3))+"E"+fnstrnum$(logmom)
  mech$(isor)=mech$(isor) + mom$
  art$=art$+" "+mech$(isor)
next isor
return
'=============================================================================
SETGSCALE:
'--- get scaling for plot first
'--- rotated, two-page format

  if i2page or iplform=2 then 
       irotate=ktrue : xoff=-15 : yoff=10 : gscale=1.0
       if igraph<>2 then yoff=0
'--- unrotated, single-page format, a gscale of 0.706 and yoff=160 mm shrinks 
'--- the full plot into half a page, but don't use the top row of SH
  else
       irotate=kfalse: xoff=15 : yoff=170 : gscale=0.706
  end if
return
'=============================================================================
FSMENU:
'-- menu to adjust focal spheres
    call scclear
      call scbox(1,12,ktrue,"FOCAL SPHERE PLOT CONTROL","Press Esc to exit")
      call scprint(3,2,"KEY"):call scprint(3,35,"P"):call scprint(3,43,"SH")
      call scprompt(5,3,"F - Focal sphere plotted")
      call scprompt(6,3,"P - Position number")
      call scprompt(7,3,"X - X position (mm)")
      call scprompt(8,3,"Y - Y position (mm)")
      call scprompt(9,3,"R - Radius - (mm)")
      call scprompt(10,3,"D - use Default values for current format")

'-- now write variables
      call setvarcol
      call scprint(5,33,fnyesno$(kfoc(1))):call scprint(5,41,fnyesno$(kfoc(2)))
      locate 6,28:print using fniform$(8);ifocpos(1);ifocpos(2);
      locate 7,28:print using fniform$(8);ifocx0(1);ifocx0(2);
      locate 8,28:print using fniform$(8);ifocy0(1);ifocy0(2);
      locate 9,28:print using fniform$(8);ifocrad0(1);ifocrad0(2);

      call scprintc(13,1,"Enter KEY for OPTION -->")
       call inlet(kfalse, e$, klen, kval, kasc)
       if fniesc(e$) then return
       locate 15,1
       for k=1 to 2
         if k=1 then p$="  P " else p$="  SH"
               
         if e$="p" then 
           print p$;
           input " position number (0 to center) "; isp
             if fninlim(isp,0,npos) then ifocpos(k)=isp
             if ifocpos(k)>0 then 
                ifocx0(k)=ipxpos(ifocpos(k),iplform)+boxlen0(iplform)/2
                ifocy0(k)=ipypos(ifocpos(k),iplform)+boxhigh0(iplform)/2
             end if
         elseif e$="d" then
            ifocpos(k)=0
            gosub GETFORMS
         elseif e$="f" then
            print p$; : input " plotting focal sphere (y/n) "; y$
            kfoc(k)=fnkyes(y$)
         elseif e$="x" then
           print p$; : input  " x position (10-240 mm) ";x
             if fninlim(x,10,240) then ifocx0(k)=x
         elseif e$="y" then
           print p$; : input  " y position (10-190 mm) ";x
             if fninlim(x,10,190) then ifocy0(k)=x
         elseif e$="r" then
           print p$; : input  " radius (0-100 mm) ";x
             if fninlim(x,0,100) then ifocrad0(k)=x
         end if
        next
   goto FSMENU
return
'=============================================================================
rem $INCLUDE: 'm5auxpl'
rem $INCLUDE: 'm5ptbax'
rem $INCLUDE: 'm5wparm'
'=============================================================================

