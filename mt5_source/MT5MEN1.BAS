' main MT5MEN1  (MT5.00)
' menu 1 for MT5

DEFINT I-N: OPTION BASE 1

rem $INCLUDE: 'M5COMMON'
rem $INCLUDE: 'M5FUNCTN'

if not imt5 then call notalone

kone=1
call setscreen(0)

dim nw(nphtype),sc1(nphtype),sc2(nphtype),am(3,3),eigv(3)

'---  functions for plotting
def fnip(y)=ipy-y*kscry/4*scl/(-(pext=0.0)-(pext<>0.0)*pext)
def fnit(i)=ipx+i*dt*xpixel
def fnix(x)=x+xctr
def fniy(y)=yctr+y*scratio

'--- some formats
varfm$=" "+fndform$(4,3)
fmyn$=fnaform$(40)+fnaform$(3)
fmm$=fnfform$(7,2)+fnaform$(3)+fnfform$(7,2)
fm5m$=fnfform$(5,2)+fnaform$(3)+fnfform$(5,2)
fmmm$=fndform$(6,3)+fnaform$(3)+fndform$(6,3)
fmtau$=fniform$(5)+"("+fniform$(3)+")"

'-- setup first time except if parameter file exists
if not icont then 
  isorm=1             'set source number
  if not iparflg then
    gosub MODELIN       'starting model
    ichoice=99          'get inversion parameters
    chain SYNDIR$+"MT5MEN2"
'    gosub CHINVPARS     'inversion parameters
  end if
end if

410 '---------- Main Menu
menu2on=kfalse:lnprompt=24
if mtflag then mechnm$="MOMENT TENSOR" else mechnm$="DOUBLE COUPLE"
call scbox(1,19,ktrue,"MAIN MENU FOR "+mechnm$+" MODEL","Press O for full OPTION list")
call scprint(3,2,"I/O file prefix  ")
call scprint(4,2,"Source (subevent) number ")
call scprint(4,32,"of")
call scprint(5,2,"Number of time function elements (NTAU) =")
call scprint(6,2,"Length of time function elements (TAU0) =")

call m5stf (kfalse,tau0,isorm,ntau(),atau(),tam(), tcent, t95) ' get duration

if mtflag then
  nparmech=6 : noffmech=noffmt
  call scprint(8,2,"         Mnn   Mee   Mdd   Mne   Mnd   Med Depth  Dur  dTo Dist   Az    Moment")
else
  nparmech=3 : noffmech=0
  call scprint(8,2,"      Strike   Dip  Rake Depth  Dur  dTo Dist   Az    Moment")
end if

call scprint(9,2,"Model ")
call scprint(10,2,"Errors")
if mtflag then
  call scprint(11,59,"Deviatoric=")
  call scprint(12,60,"Isotropic=")
end if
call scprint(12,2,"MODEL VARIANCES: Unwtd=          Wtd=")
call scprint(14,2,"VARIANCES:Total")
k=0
for i=1 to nphtype
  if kphflag(i) then k=k+1: call scprint(14,22+(k-1)*9,fnkphlab$(i))
next

call scprint(15,2,"Data ")
call scprint(16,2,"R/D %")
call scprint(18,2,"T-test:")
if istflag then
  call tstats(ns,avd,sdd,t)
  call scprint(18,10,"DOF="):call scprintv(18,14,str$(ns-1))
  call scprint(18,(0)," Mean Dif= "):call scprintvfu(18,(0),fndform$(5,3),avd)
  call scprint(18,(0),"  SD Dif= "):call scprintvfu(18,(0),fndform$(5,3),sdd)
  call scprint(18,(0),"  t="):call scprintvfu(18,(0),fnfform$(5,2),t)
else
  call scprintv(18,10,"OFF")
end if

'-- now write variables

call setvarcol
call scprint(3,19,"<"+FHDR$+">")
call scprint(4,27,str$(isorm)): call scprint(4,35,str$(nsor))
call scprint(5,43,str$(NTAU(isorm)))
call scprint(6,43,str$(TAU0))

locate 9,8

iok=mtflag or sorpar(1,isorm)>-998

for i=1 to nparmech
  if iok then
    x=sorpar(i+noffmech,isorm)
    if mtflag then
      print using fnfform$(6,2); x;
    else
      x=fndeg(x)
      print using fnfform$(6,1); x;
    end if
  else
    print "  N/A ";
  end if
next

print using fnfform$(6,1); sorpar(4,isorm)+cthk;
print using fnfform$(5,1); t95;
for i= 5 to 6 : print using fnfform$(5,1); sorpar(i,isorm); : next
print using fnfform$(5,1); fndeg(sorpar(7,isorm));
print using fndform$(6,3); tam(isorm);
if mtflag then
  call mtiso(sorpar(),isorm,amdev,amiso)
  tamiso(isorm)=tam(isorm)*amiso
  locate 11,70 : print using fndform$(6,3); tam(isorm)*amdev;
  locate 12,70 : print using fndform$(6,3); tamiso(isorm);
end if

if itf>0 and lastit=mtflag then
  locate 10,8
  for i=1 to nparmech
    if damp(i+noffmech,isorm)>0 then
      call prtform(6,2,ermech(i,isorm))
    else
      locate ,pos(0)+6
    end if
  next
  if damp(4,isorm)>0 then
    call prtform(6,1,erpar(1,isorm))
  else
    locate ,pos(0)+6
  end if
  if stfdamp(isorm)>0 then
    call prtform(5,1,ert95(isorm))
  else
    locate ,pos(0)+5
  end if
  for i=2 to 4
    if damp(3+i,isorm)>0 then
      call prtform(5,1,erpar(i,isorm))
    else
      locate ,pos(0)+5
    end if
  next
  print using fndform$(6,3); ermom(isorm);

  locate 12,25 : print using fndform$(4,3); varnowt*varcor;
  locate 12,39 : print using fndform$(4,3); varres*varcor;

  locate 15,8
  print using varfm$; vardat*varcor;
  for i=1 to nphtype
    if kphflag(i) then print using varfm$; datvar(i)*varcor;
  next

  locate 16,8
  if vardat<>0 then t=varres/vardat*100. else t=0
  t=fnamin(t,999)
  print using fnfform$(9,1); t;
  for i=1 to nphtype
    if datvar(i)>0 then t=resvar(i)/datvar(i)*100. else t=0
    t=fnamin(t,999)
    if kphflag(i) then print using fnfform$(9,1); t;
  next
end if

klst$="abcdflmopqrstvwxz"
call scprintc(20,1,"OPTIONS:")
call scprompt(20,10,"A-summary")
call scprompt(20,22,"M-model menu")
call scprompt(20,38,"C-inversion menu")
call scprompt(20,55,"S-screen plot")

call scprint(21,10,"I-iterate")
if iok then klst$=klst$+"i": call scprintc(21,10,"I")
call scprompt(21,22,"F-focal spheres")
call scprompt(21,38,"L-station menu")
call scprompt(21,55,"W-write .PRM file")

call scprompt(22,2,"Q-quit")
call scprompt(22,10,"P-plot menu")
call scprompt(22,22,"X-DOS shell")
call scprompt(22,38,"Z-setup menu")
if mtflag then
  call scprompt(22,55,"T-DC analysis")
else
  call scprompt(22,55,"T-MT analysis")
end if
if nsor>1 then
  call scprompt(23,55,"#-event shown(1-"+fnstrnum$(nsor)+")")
  klst$=klst$+numlist$(1,nsor)
end if

it$=""

600  itimerow=2:itimecol=61

ik=input1(lnprompt,1,"Enter KEY for option ",it$,kfalse,menu2on,ktrue,klst$)
if ik=3 then it$="q" else if ik=27 then 410
kval=val(it$)
itimerow=0 ' turn off display of time
ichoice=asc(it$)
if fninlim(kval,1,nsor) then isorm = kval : goto 410
if not fninlim(ichoice,97,122) then 600

if ichoice =  97  then chain SYNDIR$+"MT5MEN2"        '-- A
if ichoice =  98  then chain SYNDIR$+"MT5MEN2"        '-- B
if ichoice =  99  then chain SYNDIR$+"MT5MEN2"        '-- C
'if ichoice =  99  then gosub CHINVPARS                 '-- C
if ichoice = 100  then chain SYNDIR$+"MT5MEN2"        '-- D
if ichoice = 102  then chain SYNDIR$+"MT5MEN2"        '-- F
'if ichoice = 104  then gosub HELP                      '-- H
if ichoice = 105  then chain SYNDIR$+"MT5IT"          '-- I
if ichoice = 108  then chain SYNDIR$+"MT5MEN2"        '-- L
'if ichoice = 108  then gosub STADMENU                  '-- L
if ichoice = 109  then gosub MODELIN                   '-- M
if ichoice = 111  then gosub MENU                      '-- O
if ichoice = 112  then chain SYNDIR$+"MT5PLOT"        '-- P
if ichoice = 113  then call ENDCHECK:                  '-- Q
if ichoice = 114  then gosub REALIGN                   '-- R
if ichoice = 115  then gosub SEISPLOT                  '-- S
if ichoice = 116  then gosub TOGGLE                    '-- T
if ichoice = 118  then chain SYNDIR$+"MT5MEN2"        '-- V
if ichoice = 119  then gosub M5WPARM                   '-- W
if ichoice = 120  then call  SHELLCOM                  '-- X
if ichoice = 122  then chain SYNDIR$+"MT5MAIN"        '-- Z
GOTO 410

'********************************
'********  SUBROUTINES  *********
'********************************
'-------------------------------------------------------------------------
MENU:
call scclear
call scbox(1,22,ktrue,"MAIN MENU OPTIONS","Press Esc to return to MAIN MENU")
menu2on=ktrue:lnprompt=23
klst$="abcdflmopqrstvwxz"
call scprint(3,3,"KEY    OPTION DESCRIPTION")
call scprompt(4,4,"A - Send solution summary information to file or printer")
call scprompt(5,4,"B - Generate arrival time data file")
call scprompt(6,4,"C - Change inversion parameters (INVERSION PARAMETER MENU)")
call scprompt(7,4,"D - Generate output .INV file of synthetics or residuals")
call scprompt(8,4,"F - Plot focal spheres and source time functions on screen")
call scprompt(9,4,"I - Iterate and generate synthetic seismograms")
if iok then klst$=klst$+"i" else call scprint(9,4,"I")
call scprompt(10,4,"L - List and change seismogram weights, windows, polarity")
call scprompt(11,4,"M - Change model parameters (MODEL MENU)")
call scprompt(12,4,"P - Generate plot file for hardcopy (PLOT MENU)")
call scprompt(13,4,"Q - Quit program")
call scprompt(14,4,"R - Realign and adjust data window of seismograms")
call scprompt(15,4,"S - Plot seismograms on screen")
if mtflag then
  call scprompt(16,4,"T - Toggle to Double Couple Analysis")
else
  call scprompt(16,4,"T - Toggle to Moment Tensor Analysis")
end if
call scprompt(17,4,"V - View response functions and elementary seismograms")
call scprompt(18,4,"W - Write parameter file")
call scprompt(19,4,"X - Shell (run DOS commands)")
call scprompt(20,4,"Z - Change structure, TAU, data file, t* (return to SETUP MENU)")
call scprint(21,4,"# - Change source (subevent) displayed in menus")
if nsor>1 then
  call scprint((0),(0)," (range 1-"+fnstrnum$(nsor)+")")
  call scprintc(21,4,"#")
  klst$=klst$+numlist$(1,nsor)
end if
return 600

'-------------------------------------------------------------------------

MODELIN: '--- subroutine to change model

call scclear
klst3$="#evhoxymntq"+numlist$(1,nsor)
a$="":itimerow = 2:itimecol = 61
call scbox(1,21,ktrue,mechnm$+" MODEL MENU","Press Esc to exit")
call scprint(3,2,"KEY OPTION                                CURRENT VALUE")
call scprompt(4,3,"# - event number")
call scprompt(5,3,"E - number of events")
if mtflag then
  klst3$=klst3$+"a"
  call scprompt(6,3,"A - moment tensor elements")
  call scprint(6,40,"Mnn") : call scprint(6,60,"Mne")
  call scprint(7,40,"Mee") : call scprint(7,60,"Mnd")
  call scprint(8,40,"Mdd") : call scprint(8,60,"Med")
else
  call scprompt(6,3,"S - strike angle (cwhse from N)")
  call scprompt(7,3,"D - dip angle")
  call scprompt(8,3,"R - rake angle")
  klst3$=klst3$+"s"
  iok=sorpar(1,isorm)>-998
  if iok then klst3$=klst3$+"dr" else call scprint(7,3,"D"):call scprint(8,3,"R")
end if
call scprompt(9,3,"V - Variable range parameters                                BEG   END   INC")
call scprompt(10,3,"H - depth of source in km")
call scprompt(11,3,"O - time offset from event 1")
call scprompt(12,3,"X - distance offset from event 1")
call scprompt(13,3,"Y - azimuth from event 1")
call scprompt(14,3,"M - moment (N-m)")
call scprompt(15,3,"N - # of STF elements")
call scprompt(16,3,"T - amplitudes of STF elements (errors in parentheses):")

call setvarcol
locate 4,43:print using fniform$(7);isorm;
locate 5,43:print using fniform$(7);nsor;
if mtflag then
  for i=1 to 6
    locate ((i-1) mod 3)+6,(i>3)*-20+43:print using fm5m$; sorpar(i+noffmt,isorm);" +-";ermech(i,isorm);
  next
else
  for i=1 to 3
    locate i+5,43
    if iok then 
      print using fmm$;fndeg(sorpar(i,isorm));" +-";ermech(i,isorm);
    else
      print " N/A";
    end if
  next
end if
locate 9,43: print vrp$(isorm);
locate 10,43: print using fmm$;sorpar(4,isorm)+cthk;" +-";erpar(1,isorm)
locate 11,43: print using fmm$;sorpar(5,isorm);" +-";erpar(2,isorm)
locate 12,43: print using fmm$;sorpar(6,isorm);" +-";erpar(3,isorm)
locate 13,43: print using fmm$;fndeg(sorpar(7,isorm));" +-";erpar(4,isorm)
locate 14,40: print using fmmm$;tam(isorm);" +-";ermom(isorm)
locate 15,43: print using fniform$(7);ntau(isorm)
for i=1 to 4
  if instr(vrp$(isorm),mid$("HOXY",i,1)) then
    locate 9+i,62
    vb=varbeg(i,isorm):ve=varend(i,isorm):vi=varinc(i,isorm)
    if i=1 then vb=vb+cthk:ve=ve+cthk
    if i=4 then vb=fndeg(vb):ve=fndeg(ve):vi=fndeg(vi)
    print using "###.# ";vb; ve; vi;
  end if
next
nl=0
for i=1 to ntau(isorm)
  if (i mod 7) =1 then nl=nl+1 : locate 16+nl,3
  print using fmtau$; atau(i,isorm)*100;ertau(i,isorm)*100;
next

if input1(22,1,"Enter KEY for option ",a$,kfalse,ktrue,ktrue,klst3$)=27 then return
kval=val(a$)

if fninlim(kval,1,nsor) then isorm = kval : goto MODELIN

if fniquit(a$) then return

if a$="#" then
  input "Source number"; n
  if fninlim(n,1,nsor) then isorm=n

elseif a$="t" then
  print"Enter relative amplitudes for time function"
  for i=1 to ntau(isorm): input atau(i,isorm) : next
  call m5stf(kfalse,tau0,isorm,ntau(),atau(),tam(),tcent,t95)

elseif a$="v" then
  a$=vrp$(isorm)
  print "Parameter(s) will be varied from BEG to END with INC step size";
  if inputn(23,1,"List parameters to vary (any of: HOXY)",a$,kfalse,ktrue,ktrue,"a",4)<>27 then
     a$=RTRIM$(UCASE$(a$))
     vrp$(isorm)=""
     for i=4 to 7
       damp(i,isorm) = abs(damp(i,isorm)) 'set multiplier on
     next
     for i=1 to len(a$)
       j=instr("HOXY",mid$(a$,i,1))
       if j>0 then
         call cline(24,1): call cline(23,1): call cline(22,1)
         vrp$(isorm)=vrp$(isorm)+mid$("HOXY",j,1)
         vnm$=" value for "+mid$("dpthtimedistazim",j*4-3,4)
         if j>1 then vnm$=vnm$+" offset"
         print "BEGinning"+vnm$;: input ""; vb
         print "ENDing"+vnm$;: input ""; ve
         print "INC"+vnm$;: input ; "" ; vi 
         if vb>ve then vi= -abs(vi) else vi= abs(vi)
         if j=1 then vb=fnamax(vb-cthk,0.1): ve=fnamax(ve-cthk,0.1)
         if j=4 then vb=fnrad(vb): ve=fnrad(ve): vi=fnrad(vi)
         varbeg(j,isorm)=vb: varend(j,isorm)=ve: varinc(j,isorm)=vi
         damp(3+j,isorm) = -abs(damp(3+j,isorm)) 'set multiplier off
       end if
     next i
  end if
elseif a$="a" then
  call scprint(22,3,"Moment tensor elements: 1-Mnn,2-Mee,3-Mdd,4-Mne,5-Mnd,6-Med,A-all")
  ij=0
  do
    if input1(23,3,"Moment tensor element to change (1-6,A,Esc) ",a$,kfalse,ktrue,kfalse,"123456a") <>27 then
      ij=ij+1
      locate 24,3 : input "MOMENT TENSOR value (as ratio of Total Moment)"; d
      if a$="a" then 
        for i=1 to nparmech 
          sorpar(i+noffmt,isorm)=d
          locate ((i-1) mod 3)+6,(i>3)*-20+43:print using fm5m$; d;
        next
      else
        i=val(a$)
        sorpar(i+noffmt,isorm)=d
        locate ((i-1) mod 3)+6,(i>3)*-20+43:print using fm5m$; d;
      end if
    else
      exit do
    end if
  loop
  if ij>0 then
    if ij>1 then i=7
    tamhld=tam(isorm) 'do not change tam()
    call mtadj(sorpar(),isorm,tam(),i) 'normalize mt elements
    tam(isorm)=tamhld
  end if

elseif a$="d" then
  input "DIP         "; d: sorpar(2,isorm)=fnrad(d)

elseif a$="e" then
  input "Number of sources "; n
  if fninlim(n,1,nsormax) then nsor=n
  isorm=fnamin(isorm,nsor)

elseif a$="h" then
  input "SOURCE DEPTH in km"; d: d=d-cthk
  sorpar(4,isorm)=fnamax(d,0.1)

elseif a$="m" then
  input "MOMENT (N-m)"; t
  if t>0 then tam(isorm)=t

elseif a$="n" then
  input "NTAU"; n
  if not fninlim(n,1,ntaumax) then
    print "Invalid number"
  else
    ntau(isorm)=n
    call m5stf(kfalse,tau0,isorm,ntau(),atau(),tam(),tcent,t95)
  end if

elseif a$="o" then
  input "ORIGIN TIME OFFSET in s ";d : sorpar(5,isorm)=d

elseif a$="r" then
  input "RAKE        "; d: sorpar(3,isorm)=fnrad(d)

elseif a$="s" then
  input "STRIKE      "; d: sorpar(1,isorm)=fnrad(d)

elseif a$="x" then
  input "DISTANCE OFFSET in km ";d : sorpar(6,isorm)=d

elseif a$="y" then
  input "AZIMUTH OFFSET in degrees ";d : sorpar(7,isorm)=fnrad(d)

end if
goto MODELIN

'-------------------------------------------------------------------------
REALIGN:  '-------- realign or window ----------------
call scclear : input "Realign(0) or set window(1)"; k
if not fninlim(k,0,1) then return
if k=0 then kreal=ktrue else kreal=kfalse
input "Realign/window automatically(0) or by hand(1)"; kr
if kr=1 and itf=0 then print "Synthetics not yet calculated": d$=input$(1): return

if kreal and kr=0 then 
  gosub AREALIGN
elseif not kreal and kr=0 then 
  gosub AWINDOW
elseif kr=1 then 
  gosub HREALIGN
end if
call setscreen (0)
return

'-------------------------------------------------------------------------
HREALIGN:  '-------- realign/window by hand -------
input "Starting station number to realign/window"; ista
if ista=0 then return
imove=1

if invflag then 
  input "Seismogram amplitudes normalized (y/n)"; y$ 
  if fnkansr(y$) then normamp=fnkyes(y$)
end if

call setscreen (kega)
kd=1:ia1=0:ia2=0:it1=0:it2=0

4150 ' plot and realign
if ista<1 then ista=1
if ista> nsta then call scclear : return

if not iflag(ista) or nodat(ista) then
  ista=ista+kd
  if ista<1 then kd=+1
  goto 4150
end if

ndatin=fnamin(2*nes(ista), ndat(ista))
call m5synio(1,0,ista,ndatin,0,odat(),dext,datsq)      'get data
call m5synio(1,1,ista,nes(ista)-1,0,ss(),spext,synsq)  'get synthetics

ipx=ksrcx/16 : ipy=kscry/2
xnd=(kscrx-ipx)/ndatin

'-- correct data for geometrical spreading
for i=1 to ndatin : odat(i)=odat(i)*gcor(ista) : next

if not normamp then
  pext=fnamax(spext,dext*gcor(ista))
  scl=1.3
else
  scl=2.0
end if

4151
call scclear
call scprintb(1,1,string$(80,"="))
call scprintb(1,20," MT5 ")
call scprintc(1,58," (Press Esc to exit) ")

ip2=ipf2(ista)

'-- plot synthetics
if normamp then pext=synrms(ista)

if kega>1 and (kgrcol(kcolsyn)<>kgrcol(kcoldat)) then
  ix=ipx+xnd*ip2 : iy=fnip(ss(1)) : call scpoint (ix,iy,kgrcol(kcolsyn))
  for i=2 to nes(ista):ix=ipx+xnd*(i+ip2-1):iy=fnip(ss(i)):call scline1(ix,iy,kgrcol(kcolsyn),kfalse):next
else
  for i=1 to nes(ista):it=ipx+xnd*(i+ip2-1):ia=fnip(ss(i)):call sccircle (it,ia,kone,kgrcol(kcolbak)) : next
end if

'-- plot data
if normamp then pext=datrms(ista)*gcor(ista)

call scpoint(ipx,fnip(0),kgrcol(kcoldat))
for i=1 to ndatin
  it=ipx+xnd*i : ia=fnip(odat(i))
  call scline1 (it,ia,kgrcol(kcoldat),kfalse)
next
a$=hypoc$
osec=fnvalue(7,2)*3600+ fnvalue(9,2)*60 + fnvalue(11,3)/10.+ tstime(ista)

4152 '-- erase old window tics and put in new ones 
call scline2(it1,ia1,it1,ia2,0,kfalse)
call scline2(it2,ia1,it2,ia2,0,kfalse)

'-- get variance for alignment
call chkbound(ista)
gosub VARCALC     'get new VAR for station

it1=ipx+xnd*ip2:ia1=ipy+5:ia2=ipy-5 : call scline2(it1,ia1,it1,ia2,kgrcol(kcolbak),kfalse)
it2=ipx+xnd*(ip2+ninv(ista))        : call scline2(it2,ia1,it2,ia2,kgrcol(kcolbak),kfalse)

4155 '-- write out seismogram information
call scprint(2,2,fnstanm$(ista))
call scprint(2,6,fnkphlab$(kphtype(ista)))
call scprint(2,10," Var=") : print using varfm$;stavar(ista)*varcor;
print " T res=";
print using fnfform$(5,2); (ipf2(ista)-1)*dt(ista)-thtime(ista);
asec=osec+ (ipf2(ista)-1)*dt(ista) : call hrminsec(asec,i,j,p)
print " Time=";i;":";j;":"; : print using fnfform$(5,2);p;
locate 1,30 : if kreal then print " REALIGNING "; else print " WINDOWING ";

4160 '-- set cursor speed
call setstep(imove)
d$=str$(imove) : if imove<10 then d$=" "+d$
call scprintc(nrowchar-1,1,"OPTIONS:   Left   Right   Fast   Slow   Change   Next   Previous  step="+d$)
while inkey$<>"": wend  'empty buffer
call inkey(d$)

if fniquit(d$) or fniesc(d$) then call scclear :return

if d$="n" then kd=+1 : ista=ista+1 : goto 4150
if d$="p" then kd=-1 : ista=ista-1 : goto 4150
if d$="c" then kreal=not kreal : goto 4155

if d$="f" or fnmoveup(d$) then imove=imove+1 : goto 4160
if d$="s" or fnmovedn(d$) then imove=imove-1 : goto 4160

mleft=fnmovelt(d$) or d$="l"
mright=fnmovert(d$) or d$="r"
if mleft and kreal then ipf2(ista)=ipf2(ista)-imove : goto 4161
if mright and kreal then ipf2(ista)=ipf2(ista)+imove : goto 4161
if mleft and (not kreal) then ninv(ista)=ninv(ista)-imove : goto 4161
if mright and (not kreal) then ninv(ista)=ninv(ista)+imove : goto 4161

goto 4160

4161
if not kreal then goto 4152
goto 4151

'-------------------------------------------------------------------------
AWINDOW: '---- window automatically
call scclear
call scbox(1,14,ktrue,"AUTOMATIC WINDOWING","Press Esc to exit")
call scprint(5,3,"For each Phase-Data Type (PDT) enter TIME WINDOW in seconds to invert")
call scprint(7,7,"TIME WINDOW <0 will set window to 0")
call scprint(8,7,"TIME WINDOW =0 will skip PDT")

for i=1 to nphtype
  if kphflag(i) then
    call scprint(10,5,"For "+fnkphlab$(i))
    locate 12,5: input " Enter TIME WINDOW ";a$
    if fniesc(a$) then return
    sc1(i)=val(a$)
    call scprint(12,25,"       ")
  end if
next

for ista=1 to nsta
  if wt(ista)>0 and (not nodat(ista)) and ninv(ista)>0 then
    n=sc1(kphtype(ista))/dt(ista)
    if n>0 then
      ninv(ista)=n
      call chkbound(ista)
    elseif n<0 then
      ninv(ista)=0
    end if
  end if
next

return
'-------------------------------------------------------------------------
AREALIGN: '----automatic realignment by cross-correlation 

call scclear
call scbox(1,15,ktrue,"AUTOMATIC REALIGNMENT","Press Q to exit")
call scprint(3,4,"Realignment by cross-correlation")
call scprint(5,3,"For each Phase-Data Type (PDT) enter TIME WINDOW for realignment")
call scprint(6,4,"and TOLERANCE, both in seconds")
call scprint(8,6,"0 TIME WINDOW will skip PDT")
ir=0

for i=1 to nphtype
  if kphflag(i) then
    call scprint(11,3,"For "+ fnkphlab$(i))
    call scprint(12,3,string$(20,32))
    call scprint(13,3,string$(20,32))
    locate 12,3 : input "TIME WINDOW ";a$
    if fniquit(a$) then return
    sc1(i)=val(a$)
    locate 13,3: input "TOLERANCE   ";a$
    if fniquit(a$) then return
    sc2(i)=val(a$)
  end if
  if (not ir) then ir= (sc1(i)>0 and sc2(i)>0)
next

'-- ir is the flag to see if any realignments to be done
if (not ir) then
  locate 17,3: 
  input "No realigning done; do you want to reset alignments to original (y/n) "; y$
  if fnkyes(y$) then
    for ista=1 to nsta : ipf2(ista)=ipf1(ista) :  next
    return
  end if
end if

call scprint(17,3,"  Realigning ... ")
'-- realign
for ista=1 to nsta
  if (not iflag(ista)) or nodat(ista) then 4380

  dt=dt(ista)
  k=kphtype(ista)
  if sc1(k)=0 or sc2(k)=0 then 4380

  ntol=cint(sc2(k)/dt)
  nreal=fnamin(cint(sc1(k)/dt),ndat(ista)-ipf1(ista)-ntol)

  if nreal<1 then 4380 
'  erase odat,ss
  redim odat(NESMAX),ss(2*NESMAX)
  ip=ipf2(ista) 
  ipmax=fnamin(ipf1(ista)+ntol,ndat(ista)) 
  ipmin=fnamax(ipf1(ista)-ntol,1)

  call m5synio (1,0,ista,ndat(ista),0,odat(),amax,datsq)  'get data
  call m5synio (1,1,ista,nes(ista),0,ss(),amax,synsq)     'get synthetics

  summax=0
  for ir1=ipmin to ipmax
    ireal=ir1
    gosub XCOREL
    if sum>summax then summax=sum : ip=ir1
  next

4370  ipf2(ista)=ip
4375  call chkbound(ista)
  if istflag then gosub VARCALC

4380 next ista

return
'========================================================================
XCOREL: 'cross-correlate NREAL points starting at data point IREAL, return SUM
sum=0  
ireal=ireal-1
nr=fnamin(nreal,ndat(ista)-ireal)
if ireal<0 or nr<1 then return
for i=1 to nr
  sum=sum+ss(i)*odat(i+ireal)
next
return
'========================================================================
SEISPLOT: '-----------plotting seismograms on screen
cls
if itf=0 then kpsyn=kfalse

call scprintc(2,1," Change screen plotting parameters (y/n)>")
call inkey(y$)
if fnkyes(y$) then gosub SCREENPLOT
ista=0 : scl=1.0
ioff=-iploff*kscry/20

'-- find max of data
datext=0 : kdat=kfalse : tmax=0
for i=1 to nsta
  if not nodat(i) then kdat=ktrue
  if itf=0 then gcor(i)=1.0    'GCOR not yet calculated
  if iflag(i) then 
    datext=fnamax(datext,abs(dext(i)*gcor(i)))
    tmax=fnamax(tmax,nes(i)*dt(i))
  end if
next
pext=datext

'-- if no data
if not kdat then normamp=ktrue : kpdat=kfalse
if not kpdat and not kpsyn then goto ENDSEIS

call setscreen(kega)

8070 nplt=0 : nplotx=npps/2 : nploty=2
if npps=1 then nplotx=1: nploty=1
if npix=0 then xpix=0.9*(kscrx)/(nplotx*tmax) else xpix=npix

for iipy = 1 to nploty
  ipy= kscry/(2*nploty) + (iipy-1)*kscry/2   'plots NPPS grams per screen
  iys= 1 +(iipy-1)*(nrowchar/2)

  for iipx = 1 to nplotx             'in 2 rows of npps/2
    ipx=(iipx-1)*int(kscrx/nplotx)
    ixs=1+ int((iipx-1)*(ncolchar-1)/nplotx)
    ixs=fnxlimit(ixs,1,ncolchar-10)

    nplt=nplt+1

8135 ista=ista+1
    if ista>nsta then k=fninquit(ktrue) : goto ENDSEIS
    if not iflag(ista) then 8135
    dt=dt(ista) : k=kphtype(ista) : scl=phscl(k) : p$=fnkphlab$(k)
    xpixel=xpix
    if (ilpflag and fnishort(ista)) then xpixel=xpix*5.0 
    if scl<0 then 8135

    datsq=1.0 : synsq=1.0 : spext=0 : dext=0.

'-- get data
    nin=1
    if kpdat or (not nodat(ista)) then
      nin=fnamin(nes(ista)+ipf2(ista)+1,ndat(ista))
      call m5synio(1,0,ista,nin,0,odat(),dext,datsq)
    end if

'-- get synthetics
    if kpsyn then
      call m5synio(1,1,ista,nes(ista),0,ss(),spext,synsq)
    end if

'-- scale STF to 1st seismogram on screen
    if nplt=1 then 
       call scclear : npd=xpix : timax=60: gosub PLOTSTF
    end if

8270 s$=fnstanm$(ista)
    call scprint(iys,ixs, s$)
    call scprint(iys,ixs+5, p$)

'-- start of inversion window
    ip2=ipf2(ista)

'-- scale seismograms
    if kpdat then for i=1 to nin : odat(i)=odat(i)*gcor(ista) : next

'--- plot synthetics with line if color (KEGA) or dots if not
    if kpsyn then 
      if normamp then pext=synrms(ista)

      if kega>1 and (kgrcol(kcolsyn)<>kgrcol(kcoldat)) then
        for i=1 to nes(ista)
          ix=fnit(i+ip2-1) : iy=fnip(ss(i))+ioff
          if i=1 then call scpoint (ix,iy,kgrcol(kcolsyn)) else call scline1(ix,iy,kgrcol(kcolsyn),kfalse)
        next
      else
        kcrad= 1 'fnamax(1,cint(kscry/300))
        for i=1 to nes(ista)
          ix0=fnit(i+ip2-1) : iy0=fnip(ss(i))+ioff
          if abs(ix0-ix1)>=4 or abs(iy0-iy1)>=2 then ix1=ix0 : iy1=iy0 : call sccircle (ix1,iy1,kone,kgrcol(kcoldat))
        next
      end if
    end if

'--- plot data with solid line
    if kpdat and (not nodat(ista)) then
      if normamp then pext=datrms(ista)*gcor(ista)

      call scpoint(fnit(kzero),fnip(kzero)-ioff,kgrcol(kcoldat))
      for i=1 to nin
        it=fnit(i): ia=fnip(odat(i))-ioff
        call scline1 (it,ia,kgrcol(kcoldat),kfalse)
      next
    end if

'-- make tic marks for inversion window
    it=fnit(ip2):ia1=ipy+3:ia2=ipy-3
    call scline2(it,ia1,it,ia2,kgrcol(kcolbak),kfalse)
    it=fnit(ip2+ninv(ista))
    call scline2(it,ia1,it,ia2,kgrcol(kcolbak),kfalse)

    if fninquit(kfalse) then goto ENDSEIS
next iipx, iipy

if ista<nsta and not fninquit(ktrue) then 8070
ENDSEIS: call setscreen(0) : return

'-------------------------------------------------------------------------
VARCALC: '-- calculate Variance for station ISTA

stavar=0
ip=ipf2(ista)-1
ninv=ninv(ista)
if ninv=0 then stavar(ista)=stavar : return

if not invflag then synorm=synrms(ista) : datnorm=datrms(ista)*gcor(ista)

for i=1 to ninv
  if invflag then
    resq=odat(i+ip)-ss(i)
  else
    resq=odat(i+ip)/datnorm-ss(i)/synorm
  end if
  resq=resq*resq
  stavar=stavar+resq
next
stavar(ista)=stavar/ninv
return
'========================================================================
sub CHKBOUND (ista) static
'--------- check bounds on NINV, IPF1 and IPF2
shared ipf1(),ipf2(),ninv(),ndat(),nes()

ipf1(ista)=fnamax(1,ipf1(ista))
ipf1(ista)=fnamin(ipf1(ista),ndat(ista))
ipf2(ista)=fnamax(1,ipf2(ista))
ipf2(ista)=fnamin(ipf2(ista),ndat(ista))
ninv(ista)=fnamin(ninv(ista),ndat(ista)-ipf2(ista))
ninv(ista)=fnamax(ninv(ista),0)
ninv(ista)=fnamin(ninv(ista),nes(ista))

end sub
'========================================================================
sub TSTATS (ns,avd,sdd,t) static
'-- statistics of Huang et al. JGR, 91, 579-598, 1986. (T test)
shared nsta, stavar(), refvar(), iflag(), nodat(), varcor

t=0 : ns=0 : sumy=0 : sumyd=0 : sumyd2=0 : avd=0 : sdd=0

for i=1 to nsta
  if iflag(i) and (not nodat(i)) then
    yd= (stavar(i) - refvar(i)) * varcor
    sumy=sumy+stavar(i)*varcor
    ns=ns+1 
    sumyd=sumyd+yd 
    sumyd2=sumyd2+yd*yd
  end if
next

if ns>1 then avd=sumyd/ns: sdd=sqr( (sumyd2 - sumyd*sumyd/ns)/(ns-1) )
if sdd<>0 and ns>0 then t=avd/(sdd/sqr(ns)) 
end sub
'========================================================================
sub hrminsec ( sectime, ihr, imin, secs ) static
secs=sectime
ihr = int(secs/3600.)
secs=secs - ihr*3600.0
imin=int(secs/60.)
secs=secs - imin*60.
end sub
'=================================================================
SCREENPLOT:

call scbox(1,15,ktrue,"SCREEN PLOTTING MENU","Press Esc to exit")
call scprint(2,2,"KEY OPTION")
call scprint(2,60,"CURRENT VALUE")

call scprompt(4,3,"P - plots per screen")
call scprintv(4,60, str$(npps))

call scprompt(5,3,"N - number of pixels per second (0 default)")
call scprintv(5,60, str$(npix))

call scprompt(6,3,"O - seismograms offset (otherwise overlapping)")
call scprintv(6,60,fnyesno$(iploff))

call scprompt(7,3,"D - plotting Data")
call scprintv(7,60,fnyesno$(kpdat))

call scprompt(8,3,"S - plotting Synthetics")
call scprintv(8,60,fnyesno$(kpsyn))

call scprompt(9,3,"A - amplitudes normalized")
call scprintv(9,60,fnyesno$(normamp))

call scprompt(11,3, "M - scaling factors")
for i=1 to nphtype
  call scprint(12,6+(i-1)*8,fnkphlab$(i))
  call setvarcol
  locate 13, 6+(i-1)*8
  print using fnfform$(4,1); phscl(i);
  call settxtcol
next

a$=""
if input1(16,1,"Enter KEY for option ",a$,kfalse,ktrue,ktrue,"pnodsamq")=27 then return
if fniquit(a$) then return

' print
if a$="p" then
  input "Number of seismograms plotted per screen (1 or even number; 2,4,6,..)";n
  if (int(n/2)=n/2 or n=1) then npps=n
  npps=fnamax(1,npps)

elseif a$="n" then
  input "Seismogram plotting scale on screen, pixels per second (0 for default)"; n
  npix=fnamax(0,n)

elseif a$="o" then
  y$="y":if iploff then y$="n"
  iploff= igetyesno((0),(0),"Seismograms offset in plotting",y$,kfalse)

elseif a$="d" then
  y$="y":if kpdat then y$="n"
  kpdat=igetyesno((0),(0),"Observed seismograms plotted",y$,kfalse)

elseif a$="s" then
  y$="y":if kpsyn then y$="n"
  kpsyn=igetyesno((0),(0),"Synthetic seismograms plotted",y$,kfalse)

elseif a$="a" then
  if invflag then
    y$="y":if normamp then y$="n"
    normamp=igetyesno((0),(0),"Amplitudes normalized",y$,kfalse)
  else
    print "Inverting for shape requires amplitude normalization"
    input "Press any key "; d$
    normamp=ktrue
  end if

elseif a$="m" then
  print "Scale factors  (<0 to skip phase, 0 for no change)"
  for i=1 to nphtype
    if kphflag(i) then
      print "   For "; fnkphlab$(i);
      input s
      if s<>0 then phscl(i)=s
    end if
  next

end if
goto SCREENPLOT
return
'==================================================================  
'HELP:
'call scbox(1,22,ktrue,"HELP SCREEN","")
'
'
'call scpause(25,3)
'return

'==================================================================
TOGGLE:
mtflag= (not mtflag)
if mtflag then  'transform strike,dip,rake into moment tensor elements
'-- from Aki & Richards, p 117 (normalized to tam(isor)
  for isor=1 to nsor
    if sorpar(1,isor)>-998 then 'do not convert if indeterminate solution
      x=sorpar(1,isor) : ss=sin(x) : cs=cos(x) : s2s=sin(2.*x) : c2s=cos(2.*x)
      x=sorpar(2,isor) : sd=sin(x) : cd=cos(x) : s2d=sin(2.*x) : c2d=cos(2.*x)
      x=sorpar(3,isor) : sr=sin(x) : cr=cos(x)

      amiso=tamiso(isor)/tam(isor)
      sorpar(1+noffmt,isor) = -(sd*cr*s2s + s2d*sr*ss*ss) + amiso 'Mnn=Mxx
      sorpar(2+noffmt,isor) = sd*cr*s2s - s2d*sr*cs*cs + amiso    'Mee=Myy
      sorpar(3+noffmt,isor)= s2d*sr + amiso                       'Mdd=Mzz
      sorpar(4+noffmt,isor) = sd*cr*c2s + 0.5*s2d*sr*s2s 'Mne=Mxy
      sorpar(5+noffmt,isor)= -(cd*cr*cs + c2d*sr*ss)     'Mnd=Mxz
      sorpar(6+noffmt,isor)= -(cd*cr*ss - c2d*sr*cs)     'Med=Myz
      call mtadj(sorpar(),isor,tam(),7)
    end if
  next isor
else  'transform moment tensor elements into strike,dip,rake
  for isor=1 to nsor
    call mteigen(ktrue,sorpar(),isor,tam(),tamiso(),ierr)
    if sorpar(1,isor)<-998 then call schit("Indeterminate double couple for source"+str$(isor))
  next isor
end if
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
'==================================================================  
sub prtform (ilen,idec,var) static
if var<10^(ilen-idec-1) then
  print using fnfform$(ilen,idec); var;
elseif ilen>5 then
  print using fndform$(ilen-4,ilen-6); var;
else
  tmp$=fnstrnum$(ilen-idec-1)
  print space$(fnamax(0,ilen-3-len(tmp$)));">1e";tmp$;
end if
end sub
'==================================================================
rem $INCLUDE: 'm5auxpl'
rem $INCLUDE: 'm5ptbax'
rem $INCLUDE: 'm5stf'
rem $INCLUDE: 'm5wparm'
' rem $INCLUDE: 'm5gjinv' 'no longer used in this routine
' rem $INCLUDE: 'm5eigen' 'now linked separately
'========================================================================
