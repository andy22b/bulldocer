' main MT5MEN2  (MT5.00)
' menu 2 for MT5

DEFINT I-N: OPTION BASE 1

rem $INCLUDE: 'M5COMMON'
rem $INCLUDE: 'M5FUNCTN'

if not imt5 then call notalone

dim nw(7),sc1(nphtype)
def fnix(x)=x+xctr  'duplicate of fn in MT5MEN1
def fniy(y)=yctr+y*scratio     'ditto
varfm$=" "+fndform$(4,3)

if mtflag then
  nparmech=6 : noffmech=noffmt
else
  nparmech=3 : noffmech=0
end if

SELECT CASE ichoice
  CASE 97
    gosub PRSUM                     '-- A
  CASE 98
    gosub ATDFILE                   '-- B
  CASE 99
    gosub CHINVPARS                 '-- C
    icont=ktrue
  CASE 100
    gosub SYNFILE                   '-- D
  CASE 102
    gosub PLOTFPS                   '-- F
  CASE 108
    gosub STADMENU                  '-- L
  CASE 118
    call  ESPLOT                    '-- V
  CASE ELSE
END SELECT

chain SYNDIR$+"MT5MEN1"
'-------------------------------------------------------------------------

PRSUM:  '--------- print summary

if outf$="" then d$=fhdr$+".sum" else d$=outf$
if inputn((0),(0),"Enter file name (Esc to abort)",d$,ktrue,ktrue,ktrue,"a",25)=27 then return

if len(d$)=0 then
  outf$=fhdr$+".sum"
else
  outf$=fnaddext$(lcase$(d$),"sum")
end if
if keepold(outf$) then goto PRSUM
call cline(24,1)
call scprint(24,1,"Writing file "+outf$+" and summary lines to MT5FPS.HDR")

open "o",8,outf$
width #8, 100
nt=0

ad$=date$+"   "+time$
fo1$=fniform$(5)+"("+fniform$(3)+")"
fo2$=fnfform$(8,2)+" ("+fnfform$(5,2)+")"
fo3$=fndform$(6,3)+" ("+fndform$(6,3)+")"
fo4$=fnfform$(6,2)+" ("+fnfform$(4,2)+")"
 
print #8, "Earthquake file: ";FHDR$;"   ";ad$ : print #8,
print #8, "Hypocenter Header" : print #8, tab(5); hypoc$
print #8,
                                    
print #8, "-- Source structure --" 
print #8, "   Ztop      Vp      Vs     Rho"
x=0: if wthk>0 then print #8, using fnfform$(8,2);x;alphw;betaw;rhow; : print #8,"w"
for i=1 to nlyr
  print #8, using fnfform$(8,2); ztop(i)+wthk;alph(i);beta(i);rho(i);
  if i=1 and nlyr>1 then print #8, "c" else print #8, ""
next
if nlyr>1 then print #8, "c- complete response calculated for";nrevcom;"bounces"

print #8,
print #8, "-- Station Data --"
print #8," N Sta  FP-I PDT Po1 Po2 +/-   Win   Del    Az   TOA Tres   Wt   Var     DT   t* Filter  RefMag Calib"

for ista=1 to nsta
  dt=dt(ista)
  print #8, using fniform$(2); ista; 
  print #8, using fnaform$(6); " "+fnstanm$(ista);
  print #8, using fniform$(1);-iflag(ista); iph(ista); instin(ista) AND 7; inst(ista);
  print #8, using fnaform$(4); fnkphlab$(kphtype(ista));
  print #8, using fniform$(4); ipf1(ista); ipf2(ista); ipol(ista);
  print #8, using fnfform$(6,1); ninv(ista)*dt; del(ista); azd(ista); fndeg(ptoa(ista));
  print #8, using fnfform$(5,1); dt*(ipf2(ista)-1)-thtime(ista); wt(ista);
  if nodat(ista) then print #8, " No Data "; else print #8, using varfm$; stavar(ista)*varcor;
  print #8, using fnfform$(5,2); dt; tstar(ista);
  print #8, using fnfform$(6,2); hpfilt(ista);
  print #8, using fniform$(8); refmag(ista);
  print #8, using fnaform$(5); " "+fncaltyp$(instin(ista))
  if wt(ista)>0 and iflag(ista) and ninv(ista)>0 then nw(kphtype(ista))=nw(kphtype(ista))+1
next

print #8, 
print #8, "-- Phase Data --"
print #8, "Type  Num    Mag.   Time    ResVar DataVar    R/D %" 
for i=1 to nphtype 
  if kphflag(i) then
    print #8, using fnaform$(5); fnkphlab$(i); 
    print #8, using fniform$(4); nw(i); : nt=nt+nw(i)
    print #8, using fniform$(8); phmag(i);
    print #8, using fnfform$(8,2); phrlen(i);
    if datvar(i)>0 then r=resvar(i)/datvar(i) else r=0
    print #8, using varfm$; resvar(i)*varcor;datvar(i)*varcor;
    print #8, using fnfform$(7,2); r*100
  end if
next
print #8,

if vardat=0 then vrvd=0 else vrvd=varres/vardat
print #8, "Total";
print #8, using fniform$(4);nt;
print #8, "                ";
print #8, using varfm$; varres*varcor; vardat*varcor; 
print #8, using fnfform$(7,2); vrvd*100
print #8,

print #8,"Rupture Velocity, Azimuth, Plunge = ";
print #8, using fnfform$(8,2);vrupt;azrupt;plrupt
print #8, 

for isor=1 to nsor
  print #8, "-- Source number";isor;" --"
  call m5stf (kfalse,tau0,isor,ntau(),atau(),tam(), tcent, t95)

  str=sorpar(1,isor) : iok=str>-998
  if iok then
    dip=sorpar(2,isor) : rk=sorpar(3,isor)
    call m5auxpl (str,dip,rk,strt,dipt,raket)
    call m5ptbax (str,dip,rk,paz,ppl,taz,tpl,baz,bpl)
  end if
  d0=sorpar(4,isor)+cthk

  a$=hypoc$ 
  ts = fnvalue(7,2)*3600 + fnvalue(9,2)*60 + fnvalue(11,3)/10 + sorpar(5,isor)
  ih=int(ts/3600.0)
  ts=ts - 3600.0*ih
  im=int(ts/60.0)
  ts = ts - 60.0*im
  print #8, "Time ";ih;"h";im;"m";fnfix(ts,1);"s    ";
  xlat = fnvalue(14,5)/100 + sorpar(6,isor)*cos(sorpar(7,isor))/111.19
  xlon = fnvalue(19,6)/100 + sorpar(6,isor)*sin(sorpar(7,isor))/111.19
  print #8, " Latitude ";fnfix(xlat,2);"  Longitude "; fnfix(xlon,2)

  print #8,"Source function elements:  TAU = ";TAU0;" NTAU = ";NTAU(ISOR)

  for i=1 to ntau(isor)
    print #8, using fo1$; atau(i,isor)*100; 2*ertau(i,isor)*100;
    if i mod 8 = 0 then print #8,
  next

  print #8,
  print #8, " Half-duration = ";fnfix(tcent,1);"   95% duration = ";
  print #8, using fo2$; t95; ert95(isor)
  if mtflag then print #8,"Best double couple:" else print #8,
  print #8, "      Strike           Dip            Rake           Depth"
  if iok then
    for i=1 to 3
      if mtflag then 
        print #8, using "    "+fnfform$(8,2)+"    "; fndeg(sorpar(i,isor));
      else  
        print #8, using fo2$; fndeg(sorpar(i,isor)); 2*ermech(i,isor);
      endif
    next i
  else
    print #8,"       N/A              N/A             N/A     ";
  end if
  print #8, using fo2$; d0;2*erpar(1,isor)
  if mtflag then
    print #8, "Moment tensor elements (normalized by Moment):"
    print #8, "     Mnn          Mee          Mdd          Mne          Mnd          Med
    for i=8 to 13
      print #8, using fo4$; sorpar(i,isor); 2*ermech(i-7,isor);
    next
    print #8,
    call mtiso(sorpar(),isor,amdev,amiso)
    print #8," Deviatoric Moment = ";: print #8, using fndform$(6,3);tam(isor)*amdev;
    print #8,"  Isotropic Moment = ";: print #8, using fndform$(6,3);tam(isor)*amiso
  end if
  print #8, "        Moment             Or. Time        Distance        Azimuth
  print #8, using fo3$; tam(isor); 2*ermom(isor);
  print #8, using fo2$; sorpar(5,isor); 2*erpar(2,isor);_
                       sorpar(6,isor); 2*erpar(3,isor);_
                       fndeg(sorpar(7,isor)); 2*erpar(4,isor) 
  print #8, "Aux Pl Str,Dip    Paz    Ppl    Taz    Tpl    Baz    Bpl"
  if iok then
    print #8, using fnfform$(7,1); fndeg(strt);fndeg(dipt);fndeg(paz);_
             fndeg(ppl);fndeg(taz);fndeg(tpl);fndeg(baz);fndeg(bpl)
  else 
    print #8," N/A ...(indeterminate double couple)...."
  end if
  print #8,
  gosub SUMLINE        'output summary lines
next isor

close 8

return
'------------------------------------------------------------
SUMLINE:
'-------- print summary lines

open "a",9,"MT5FPS.HDR"
if isor=1 then print #9, fhdr$+" "+ad$

'--- output for ISC format
f1$=mid$(hypoc$,1,6) : f2$=mid$(hypoc$,28,5)
if len(f2$)<5 then f2$=f2$+string$(5-len(f2$),32)
print #9, f1$;
print #9, using fniform$(2);ih;im;
print #9, using fniform$(3);ts*10;
print #9, using fniform$(5);xlat*100;
print #9, using fniform$(6);xlon*100;
print #9, using fniform$(3);d0;
print #9, f2$;
iexp= int(fnalog10(tam(isor))) 
vm  = int(tam(isor)/(10^(iexp-2)))
print #9, using fniform$(3); vm;
print #9, using fniform$(2); iexp;
if iok then
  print #9, using fniform$(4); fndeg(str);
  print #9, using fniform$(3); fndeg(dip);
  print #9, using fniform$(4); fndeg(rk);
else
  print #9,SPC(11);
end if
if mtflag then
  for i=1 to 6
    print #9, using fnfform$(5,2); sorpar(i+noffmt,isor);
  next
end if
print #9,

'---output table format
print #9, f1$;
print #9, using fniform$(3);ih;im;
print #9, using fnfform$(5,1);ts;
print #9, using fnfform$(8,2);xlat;
print #9, using fnfform$(8,2);xlon;
print #9, using fnfform$(6,1);d0; 2*erpar(4,isor);
a$=f2$
print #9, using fnfform$(4,1); fnvalue(1,2)/10;
print #9, using fniform$(4); fnvalue(3,3);
iexpr = int(fnalog10(2.0*ermom(isor)+0.1))
vmr   = int(2.0*ermom(isor)/(10^(iexpr-2)))
print #9, using fniform$(3); vm;
print #9, using fniform$(2); iexp;
print #9, using fniform$(3); vmr;
print #9, using fniform$(2); iexpr;
print #9, using fnfform$(5,1); t95;

'***** Following values only for double couple, hence blank used for invalid
if iok then
  print #9, using fniform$(4); fndeg(str); 
  if mtflag then print #9,"    "; else print #9, using fniform$(4); 2*ermech(1,isor);
  print #9, using fniform$(3); fndeg(dip); 
  if mtflag then print #9,"   "; else print #9, using fniform$(3); 2*ermech(2,isor);
  print #9, using fniform$(4); fndeg(rk); 
  if mtflag then print #9,"    "; else print #9, using fniform$(4); 2*ermech(3,isor);
  print #9, using fniform$(4); fndeg(strt); fndeg(dipt);
  print #9, using fniform$(4); fndeg(paz); fndeg(ppl); fndeg(taz); fndeg(tpl); fndeg(baz); fndeg(bpl);
end if
for i=1 to nphtype :print #9, using fniform$(3); nw(i); : next
print #9, using fnfform$(4,1); wthk;
print #9, using varfm$; varres*varcor; 
print #9, using fniform$(4); vrvd*100;
if mtflag then
  for i=1 to 6
    print #9, using fnfform$(5,2); sorpar(i+noffmt,isor); 2*ermech(i,isor);
  next
end if
print #9,

close 9

return

'-------------------------------------------------------------------------
SYNFILE: ' make data file with synthetics or residuals
cls

if itf=0 then call schit("Synthetics not yet calculated."): return
input "Enter prefix for .INV file to be written (CR to quit)";pf$
if len(pf$)=0 then return
if pf$=fhdr$ then
  print
  input " The input file will be overwritten, OK? (y/n) "; y$
  if not fnkyes(y$) then return
end if
pf$=fnaddext$(pf$,"inv")

print
input "    Writing Synthetics(0) or Residuals(1) "; kres
if kres=1 then kres=ktrue else kres=kfalse

print "Writing file ";pf$
open "o",9,pf$  : print #9, hypoc$
dform$=" "+fndform$(6,4)

for ista=1 to nsta
  if iflag(ista) and not (nodat(ista) and kres) then

9500 
    xmag0=refmag(ista)   '''  phmag(kphtype(ista))
    xc=1.0e4 * xmag(ista) / (xmag0*gcor(ista))

'--- setup for synthetics
    call m5synio(1,1,ista,nes(ista),0,ss(),sext,synsq)  'get synthetics
    kp=9 : nin=nes(ista) : n2=nin+8

'--- set up for residuals
    if kres then 
      sext=1.1/xc : nin=ndat(ista) : n2=nin

      call m5synio(1,0,ista,nin,0,odat(),amax,datsq)  'get data
      kp=ipf2(ista) : ipo=ipf2(ista)-1 

      for j=1 to nin
        jss=j-ipo
        if jss>0 then subtr=ss(jss) else subtr=0.0
        odat(j)=odat(j)*gcor(ista)-subtr
        sext=fnamax(sext,abs(odat(j)))
      next

    end if
  
'--- output station header line

    xo=0 : xt=tstime(ista)+thtime(ista)

    print #9, using fnaform$(4); mid$(fnstanm$(ista),1,4);
    print #9, using fniform$(1); iph(ista); -((instin(ista) AND 7)>0)+(instin(ista) AND 24)\4; inst(ista);
    print #9, using fnfform$(8,2); del(ista);azd(ista);
    print #9, using fniform$(8); xmag(ista);
    print #9, using fnfform$(8,2); xt;
    print #9, using fnfform$(5,1); wt(ista);
    print #9, using fniform$(5); kp;n2;
    print #9, using fniform$(8); sext*xc;
    print #9, using fnfform$(6,2); dt(ista); xo; hpfilt(ista)

    if instin(ista) AND 7 then 
      print #9, insfile$(ista)
    elseif fningdsn(ista) then
      nz=nzero(ista) : np=npole(ista)
      print #9, using fniform$(5); nz; np;
      print #9, using dform$; calcons(ista)

      nz1=indexzp(ista)-1
      for j=1 to 2*nz
        nz1=nz1+1
        print #9, using dform$; zeropole(nz1);
        if (j mod 8 =0) or j=2*nz then print #9,""
      next j

      for j=1 to 2*np
        nz1=nz1+1
        print #9, using dform$; zeropole(nz1);
        if (j mod 8 =0) or j=2*np then print #9,""
      next j
    end if

    if not kres then
      s=0  'print row of zeroes
      for j=1 to 8 : print #9, using dform$; s; : next : print #9,
    end if

    for j3=1 to nin
      if kres then s=odat(j3) else s=ss(j3)
      print #9, using dform$; s*xc;
      if (j3 mod 8 = 0) or j3=nin then print #9, ""
    next

  end if

next  ista
close 9
return

'========================================================================
ATDFILE: ' write arrival time file
input "File name for arrival time data (q to quit, CR for default)"; d$

if len(d$)=0 then
  outf$=fhdr$+".atd"
else
  if fniquit(d$) then return
  outf$=d$
end if

print "Writing file ";outf$
open "o",1,outf$

print #1, hypoc$

'-- origin time
a$=hypoc$
osec=fnvalue(7,2)*3600.+ fnvalue(9,2)*60. + fnvalue(11,3)/10.

for ista=1 to nsta
  if (not nodat(ista)) then 
    print #1, using fnaform$(4); mid$(fnstanm$(ista),1,4);
    print #1, using fniform$(1); -iflag(ista);iph(ista);

'--- change format for output to print arrival time rather than travel time
    asec=osec+ + tstime(ista) + (ipf2(ista)-1)*dt(ista) 
    call hrminsec(asec,i,j,p)
    print #1, using fniform$(2); i;j;
    print #1, using fnfform$(5,2);p

'-- old format= travel time
'     print #1, using fnfform$(8,2); tstime(ista)+dt(ista)*(ipf2(ista)-1)

  end if
next
close 1
return
'========================================================================
CHINVPARS:  '--- change inversion parameters

call scclear
k2lst$=numlist$(1,nsor)+"#fhiklmopqtwvxy" 'also c if istflag
a$="":itimerow = 2:itimecol = 62
call scbox(1,21,ktrue,"INVERSION PARAMETER MENU","Press Esc to exit")
call scprint(3,2,"KEY OPTION                             CURRENT VALUE")

call scprompt(4,3,"M - Minimum % decrease in variance")

call scprompt(5,3,"P - Positivity of STF enforced")
if ipcflag then
  call scprint(6,7,"%tolerance")
  call scprint(7,7,"damping increment")
end if

call scprompt(8,3,"I - matching")
call scprompt(9,3,"T - station t statistics calculated")
if istflag then
  call scprompt(10,3,"C - Change reference station variance values for t-stats")
  k2lst$=k2lst$+"c"
end if
call scprompt(11,3,"K - calculate variance for 1 iter.")
call scprompt(12,3,"W - change seismogram Weights")
call scprompt(13,3,"L - Least squares damper")

call scprompt(15,3,"# - Source number")
call scprompt(16,3,"F - time Function adjustment multiplier")
call scprint(17,3,"-------------- adjustment multipliers (0=fixed) ----------------------")

if mtflag then
  call scprompt(18,3,"A-Mnn  Mee  Mdd  Mne  Mnd  Med")
  k2lst$=k2lst$+"a"
else
  call scprompt(18,5,"S-strike"):call scprompt((0),18,"D-dip"):call scprompt((0),27,"R-rake")
  k2lst$=k2lst$+"sdr"
end if
call scprompt((0),36,"H-depth"):call scprompt((0),47,"O-Time ")
call scprompt((0),54,"X-distnce"):call scprompt((0),64,"Y-azimuth")

call setvarcol
locate 4,43:print using fniform$(5); varmin
locate 5,43:print using fnaform$(3);fnyesno$(ipcflag)
locate 6,43:print using fniform$(5);cint(stfmin*100)
locate 7,43:print using fnfform$(5,1);dampfac
locate 8,43: if invflag then print "AMPLITUDE" else print "SHAPE"
locate 9,43:print using fnaform$(3);fnyesno$(istflag)
locate 11,43:print using fnaform$(3);fnyesno$(krmsc)
locate 13,43:print using fnfform$(5,1);theta
locate 15,43:print using fniform$(5);isorm
locate 16,43:print using fnfform$(5,1);stfdamp(isorm)

locate 19,3
for i=1 to nparmech :print using fnfform$(10+mtflag*5,1);damp(i+noffmech,isorm); : next
for i=4 to 7 :print using fnfform$(10,1);fnamax(0,damp(i,isorm)); : next

if input1(22,1,"Enter KEY for option ",a$,kfalse,ktrue,ktrue,k2lst$)=27 then return
kval=val(a$)
itimerow=0
if kval>0 then isorm = kval : goto CHINVPARS

if fniquit(a$) then return

if a$="s" then input "  Strike adjustment multiplier (0 to fix) ";          damp(1,isorm)
if a$="d" then input "  Dip adjustment multiplier    (0 to fix) ";          damp(2,isorm)
if a$="r" then input "  Rake adjustment multiplier   (0 to fix) ";          damp(3,isorm)
if a$="h" then input "  Depth adjustment multiplier  (0 to fix) ";          damp(4,isorm)
if a$="o" then input "  Time offset adjustment multiplier  (0 to fix) ";    damp(5,isorm)
if a$="x" then input "  Distance offset adjustment multiplier  (0 to fix) ";damp(6,isorm)
if a$="y" then input "  Azimuth offset adjustment multiplier  (0 to fix) "; damp(7,isorm)

i=instr(vrp$(isorm),ucase$(a$)) 'turn off variable range on this parameter
if i>0 then vrp$(isorm)=left$(vrp$(isorm),i-1)+mid$(vrp$(isorm),i+1)

if a$="f" then input "  STF adjustment multiplier ";                       stfdamp(isorm)
if a$="m" then input "  Minimum percent decrease in variance "; varmin

if a$="#" then
  input "Source number"; n
  if fninlim(n,1,nsor) then isorm=n

elseif a$="a" then
  call scprint(22,3,"Moment tensor elements: 1-Mnn,2-Mee,3-Mdd,4-Mne,5-Mnd,6-Med,A-all")
  if input1(23,3,"Moment tensor element to change (1-6,A)) ",a$,kfalse,ktrue,kfalse,"123456a") <>27 then
    locate 24,3
    input " Adjustment multiplier for moment tensor element (0 to fix) "; d
    if a$="a" then 
      for i=1 to nparmech : damp(i+noffmech,isorm)=d : next
    else
      damp(val(a$)+noffmech,isorm)=d
    end if
  end if

elseif a$="i" then
  input "Inversion matching AMPLITUDE(1) or SHAPE(2)"; n
  if n=1 then invflag=ktrue else if n=2 then invflag=kfalse
  if invflag then varcor=1.0e8 else varcor=1.0
  if (not invflag) then normamp=ktrue

elseif a$="c" then
  input "Set reference model statistics to current model (1) or set to zero (2)"; n
  for i=1 to nsta
    if n=1 then refvar(i)=stavar(i) else refvar(i)=0
  next

elseif a$="k" then
  input "Calculate error statistics for single iteration (y/n) "; y$
  krmsc=fnkyes(y$)

elseif a$="l" then
  input "Enter least squares damper "; theta

elseif a$="p" then
  input "Source time function Positivity enforced (y/n) "; y$
  ipcflag=fnkyes(y$)
  if ipcflag then
    input "  Change tolerance and damping increment (y/n) "; y$
    if fnkyes(y$) then
      input "   positivity tolerance for ATAU in percent ";d
      if abs(d)<>0 then stfmin=-abs(d)/100.
      input "   damping increase increment (>1) ";d
      if d>1 then dampfac=d
    end if
  end if

elseif a$="t" then
  input "Calculate station t-test statistics ON (y/n) "; y$
  istflag=fnkyes(y$)

elseif a$="w" then
  gosub RESETWT

end if

goto CHINVPARS

'-------------------------------------------------------------------------
RESETWT:   '--- adjust station weights

call scclear
call scbox(1,13,ktrue,"SEISMOGRAM WEIGHTING MENU","Press Esc to exit")
call scprint(2,2, "KEY  WEIGHT OPTION")
call scprompt(4,3, "1 - Reset ALL weights to 1.0")
call scprompt(5,3, "2 - Reset only NON-ZERO weights to 1.0")
call scprompt(6,3, "3 - Store current weights in temporary array")
call scprompt(7,3, "4 - Reset weights to those stored in temporary array")
call scprompt(8,3, "5 - Multiply weights of all seismograms of same type by same factor")
call scprompt(9,3, "6 - Scale weights to station azimuthal density")
call scprompt(10,3, "7 - Normalize so average weight = 1.0")
call scprompt(11,3, "8 - View and change weights on screen")

a$=""
if input1(14,1,"Enter KEY for option ",a$,kfalse,ktrue,ktrue,"12345678")=27 then return
ic=val(a$)

if (not fninlim(ic,1,8)) then return

if ic=8 then gosub STADMENU

'--- get weighting factors for phases
if ic=5 then
  cls
  for i=1 to nphtype
    if kphflag(i) then
      call scprintc(21,1,"Weight factor for "+ fnkphlab$(i)+" (0 to skip, <0 to set to zero) ")
      input sc
      sc1(i)=sc
      if sc=0 then sc1(i)=1
      if sc<0 then sc1(i)=0
    end if
  next
end if

'-- sum all current weights
if ic=7 then
  sumwt=0.0 : numwt=0
  for ista=1 to nsta
    if iflag(ista) and wt(ista)>0 and (not nodat(ista)) then 
      sumwt=sumwt+wt(ista) : numwt=numwt+1
    end if
  next
  if sumwt>0 then wfac=numwt/sumwt else wfac=1.0
end if

'-- adjust weights according to choice
for ista=1 to nsta

  if iflag(ista) and (not nodat(ista)) then 
    if ic=1 then wt(ista)=1.0
    if ic=3 then wt0(ista)=wt(ista)
    if ic=4 then wt(ista)=wt0(ista)

    if wt(ista)>0 then 
      if ic=2 then wt(ista)=1.0
      if ic=5 then wt(ista)=wt(ista)*sc1(kphtype(ista))
      if ic=7 then wt(ista)=wt(ista)*wfac

'--- adjust weight according to distribution of seismograms of a particular 
'--     phase in azimuth; 
'--  weight adjustment is inversely proportional to 
'--     sum of exp(-((difference in delta)/20.0) ^2)
'--  with the 20.0 in the denominator a station at 35 degrees away in azimuth
'--     will contribute <5% to weight factor 

      if ic=6 then
        azd=azd(ista) : iph=iph(ista)
        swt=0
        for j=1 to nsta
          if iflag(j) and (iph(j)=iph) then
            dis=abs(azd(j)-azd)
            dis = fnamin (dis, abs(360-dis) )
            swt=swt+exp(-(dis/20.0)^2)
          end if
        next j
        wt(ista)=wt(ista)/sqr(swt)
      end if

    end if
  end if
next ista

goto RESETWT
'-------------------------------------------------------------------------
STADMENU: '-- print station info on screen and make changes

4715 ipg=1 : npages=int(nsta/20.1)+1

4720 call scclear
call scbox(1,22,ktrue,"STATION DATA MENU","Press Esc to exit")
call scprint(2,2,"No Sta Phase Cal Del    Az   TOA   Var      Flag   +/-  Po1   Po2 Window  Wt")

for k=1 to 19
  ista=(ipg-1)*19+k : nk=k
  locate k+2,2 : gosub P1LINE
  if ista=nsta then exit for
next

4795 irow=1 : icol=1
4797 call scprintc(irow+2, 43+(icol-1)*6, "->")
call scprintc(23,2,"3-Next page 9-Prev page 6-Rght 4-Lft 8-Up 2-Dn (or arrows) C-change (or +/-)")
while inkey$<>"":wend 'empty buffer
call inkey(d$)
if fniesc(d$) then gosub RESETFLG : return
call scprintv(irow+2, 43+(icol-1)*6, "  ")
if fnipagedn(d$) or d$="n" then
  if ipg<npages then ' gosub RESETFLG : return
    ipg=ipg+1 :goto 4720
  end if

elseif fnipageup(d$) or d$="p" then
  ipg=ipg-1 : ipg=fnamax(1,ipg) : goto 4720

elseif fnmovert(d$) or d$="r" then icol=icol+1
elseif fnmovelt(d$) or d$="l" then icol=icol-1
elseif fnmoveup(d$) or d$="u" then irow=irow-1
elseif fnmovedn(d$) or d$="d" then irow=irow+1
elseif fnmovehm(d$) then icol=1
elseif fnmovend(d$) then icol=6

elseif d$="c" or fniplus(d$) or fnminus(d$) then
  mc=(d$="c") : mp=fniplus(d$) : mm=fnminus(d$)

  if mc then call cline (23,1) : locate 23,1 : input " New value "; valin
  ista= (ipg-1)*19 + irow

  if icol=1 then
    if mc then iflag(ista)= (valin<>0)
    if mm or mp then iflag(ista)=not iflag(ista)
  elseif icol=2 then
    if mc then if abs(valin)=1 then ipol(ista)=valin
    if mm then ipol(ista)=-1
    if mp then ipol(ista)=1
  elseif icol=3 then
    if mc then ipf1(ista)=valin
    if mm then ipf1(ista)=ipf1(ista)-1
    if mp then ipf1(ista)=ipf1(ista)+1
    call chkbound(ista)
  elseif icol=4 then 
    if mc then ipf2(ista)=valin
    if mm then ipf2(ista)=ipf2(ista)-1
    if mp then ipf2(ista)=ipf2(ista)+1
    call chkbound(ista)
  elseif icol=5 then 
    if mc then ninv(ista)=valin/dt(ista)
    if mm then ninv(ista)=ninv(ista)-1
    if mp then ninv(ista)=ninv(ista)+1
    call chkbound(ista)
  elseif icol=6 then 
    if mc then wt(ista)=valin 
    if mm then wt(ista)=wt(ista)-.1
    if mp then wt(ista)=wt(ista)+.1
    wt(ista)=fnamax(wt(ista),0.0)
  end if
  locate irow+2,2 : gosub P1LINE
end if

4800  irow=fnamax(irow,1)
if irow>nk then irow=1 : icol=icol+1
icol=fnxlimit(icol,1,6)
goto 4797

'-------------------------------------------------------------------------
P1LINE:
'--- writes one line for seismogram data menu
print using fniform$(2);ista;
print using fnaform$(7); " "+fnstanm$(ista);
print using fnaform$(5); fnkphlab$(kphtype(ista));
print using fnaform$(1); fncaltyp$(instin(ista));
print using fnfform$(5,1);del(ista);
print using fnfform$(6,1);azd(ista);fndeg(ptoa(ista));
if nodat(ista) then print " No Data "; else print using varfm$; stavar(ista)*varcor;
call setvarcol
print using fniform$(6); -iflag(ista); ipol(ista); ipf1(ista); ipf2(ista);
print using fnfform$(6,1); ninv(ista)*dt(ista);
print using fnfform$(6,2); wt(ista)
call settxtcol
return
'-------------------------------------------------------------------------
RESETFLG:
for i=1 to 7 : kphflag(i)=0 : next
for i=1 to nsta
  if iflag(i) then
    kphflag(kphtype(i))=ktrue
  end if
next
'--- set all phase flags
kplot(1) = (kphflag(1) or kphflag(3) or kphflag(4) or kphflag(6) or kphflag(7))
kplot(2) = (kphflag(2) or kphflag(5))
ilpflag  = (kphflag(1) or kphflag(2) or kphflag(4) or kphflag(5) or kphflag(7))
ispflag  = (kphflag(3) or kphflag(6))
igdsnflag = (kphflag(4) or kphflag(5) or kphflag(6) or kphflag(7))
return
'==================================================================
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
sub ESPLOT static
shared kphtype(), dt(), nes(), nsta, gfxi$, drive$
shared kega, kscrx, kscry, ncolchar, nrowchar, scratio,kgrcol()
shared a$,maxnes, instin()

'nmax=500
nmax=maxnes 'changed PZ- 5/10/94
redim tin(nmax), appin(nmax,2), tes(nmax), aes(nmax)
famp$=drive$+"amps.tmp"

if not ifopen("i",7,famp$,0,-1) then return
nraf=0 : n0=1

for ista=1 to nsta
  call setscreen(kega)
  if fnpwave(ista) then nph=2 else nph=1
  if fnilong(ista) then dtinc=5 else dtinc=1
  dt=dt(ista)

  t0=(nes(ista)+1)*dt : app0=0
  tscale=((kscrx-10)/t0)/2.5

  line input #7, a$
  s$=fnvalue$(9,4)
  input #7, n
  for i=1 to n
    input #7, t, a1, a2
    if i<=nmax then
      tin(i)=t : appin(i,1)=a1 : appin(i,2)=a2
      app0=fnamax(app0, abs(a1))
      app0=fnamax(app0, abs(a2))
    end if
  next
  n=fnamin(n,nmax)
  if app0=0 then app0=1

  kyoff=3/nrowchar*kscry

'--- plot ES
  lbl$=s$+fnkphlab$(kphtype(ista))+left$(fncaltyp$(instin(ista)),1)+"  Crustal response functions and seismograms"
  call scprint(1,1,lbl$)
  locate 2,4:print "Maximum amplitude of response functions:";app0;
  call scprintc(nrowchar,ncolchar-37,"Enter Q to quit; any key to continue..")
  t$="Tics:"+str$(dtinc)+"sec"
  locate nrowchar-1,3: print t$;

  for j=1 to nph
    if nph=1 then l$="sS" else if j=1 then l$="pP" else l$="sP"
    nc=ncolchar/2*(j-1)+1
    call scprint(4,nc+1,l$)

    iy0=((kscry-kyoff)/5)+kyoff : iy01=iy0
    yamp=(kscry-kyoff)/5
    ix0=(j-1)*kscrx/2+15

    call scline2 (ix0-10, kscry, ix0-10, kyoff, kgrcol(kcoldat), kfalse)

    for i=1 to n
      ix=ix0+tin(i)*tscale
      iy=iy0-yamp*appin(i,j)/app0
      call scline2 (ix, iy0, ix, iy, kgrcol(kcoldat), kfalse)
    next
    call m5tmpio (1, ista, j, nes(ista), aes(), gfxi$ )

    amax=0
    for k=1 to nes(ista)
      tes(k)=(k-1)*dt
      amax=fnamax(amax, abs(aes(k)))
    next
    if amax=0 then amax=1

    iy0=3*(kscry-kyoff)/5+kyoff
    ix=ix0+tes(1)*tscale : iy=iy0-yamp*aes(1)/amax

'-- plot seismogram
    call scpoint (ix, iy, kgrcol(kcolsyn))

    for i=2 to nes(ista)
      ix=ix0+tes(i)*tscale : iy=iy0-yamp*aes(i)/amax
      call scline1 (ix, iy, kgrcol(kcolsyn), kfalse)
    next i

    for i=0 to t0 step dtinc
      ix=ix0+i*tscale
      call scline2 (ix, iy0-kscry/80, ix, iy0+kscry/80, kgrcol(kcolbak), kfalse)
    next i

  next j
  d$=input$(1) : if fniquit(d$) then goto ESPEND

next ista

ESPEND: close 7 : call setscreen(0)
end sub
'========================================================================
PLOTFPS: '---------- plot equal area net and planes
d$="1"
if input1((0),(0),"P (1), SH (2) or SV (3) plotted?",d$,kfalse,ktrue,ktrue,"123")=27 then return
kplot=val(d$)
ktype=1
if mtflag then
  d$="3"
  if input1((0),(0),"Double Couple (1), Moment Tensor (2) or Both (3) radiation patterns?",d$,kfalse,ktrue,ktrue,"123")=27 then return
  ktype=val(d$)
end if
call setscreen(kega)
call scprintb(1,1,string$(80,"="))
if kplot=1 then 
  d$="P"
elseif kplot=2 then
  d$="SH"
else
  d$="SV"
end if
call scprintb(1,36," MT5   Focal spheres ("+d$+" Wave) ")

for isorc=1 to nsor
  isor= (isorc-1) mod 3 +1

  irad=kscrx/8 : npl=0 : xctr=kscrx/6+(isor-1)*kscrx/3 : yctr=kscry/2+20
  ix=fnix(0) : iy=fniy(0) : iyr=fniy(-irad): iyr1=fniy(-irad+10)
  call sccircle(ix,iy,irad,kgrcol(kcolsyn))
  call scline2(ix,iyr,ix,iyr1,kgrcol(kcolsyn),kfalse)

  str=sorpar(1,isorc): dip=sorpar(2,isorc): rk=sorpar(3,isorc)
  strt=str : dipt=dip : raket=rk

  if ktype<>1 then
    amnn=sorpar(1+noffmt,isorc) : amee=sorpar(2+noffmt,isorc) : amdd=sorpar(3+noffmt,isorc)
    amne=sorpar(4+noffmt,isorc) : amnd=sorpar(5+noffmt,isorc) : amed=sorpar(6+noffmt,isorc)
  end if

  for ista=1 to nsta 'plot station symbols
    if fninesc then goto ENDFPS
    if iflag(ista) then
      a=fnrad(azd(ista))
      call eaproj (irad, ptoa(ista), a, x, y) : y=-y
      ix4=fnix(x)+4 : iy2=fniy(y)+2
      ixm4=fnix(x)-4 : iym2=fniy(y)-2
      if fnpwave(ista) then
        call scline2(ix4,iy2,ixm4,iym2,kgrcol(kcolbak),ktrue)  'box for P
      else
        call scline2(ix4,iy2-2,ixm4,iy2-2,kgrcol(kcolbak),kfalse)  'plus sign for S
        call scline2(ix4-4,iy2,ix4-4,iym2,kgrcol(kcolbak),kfalse)
      end if
    end if
  next ista
  nc=(isor-1)*ncolchar/3+3 : nc=fnxlimit(nc,1,ncolchar)
  if ispflag and not ilpflag then dt=phdt(3) else dt=phdt(1)
  call scprint(2,nc,"Source") : print isorc;
  if ktype<>2 then
    iok=strt>-998
    call scprint(3,nc,"Str ") : if iok then print cint(fndeg(strt)); else print " N/A";
    call scprint(4,nc,"Dip ") : if iok then print cint(fndeg(dipt)); else print " N/A";
    call scprint(5,nc,"Rake") : if iok then print cint(fndeg(raket)); else print " N/A";

    if iok then
      if kplot=1 then
        gosub PLOTP
        call scprint(6,nc,"NP 2")
        call scprint(7,nc,"Str ") : print cint(fndeg(strt));
        call scprint(8,nc,"Dip ") : print cint(fndeg(dipt));
      elseif kplot=2 then
        gosub PLOTSRAD
      end if
    end if
  end if
  if ktype<>1 then 
    if kplot=1 then 
      call PLOTPMT(amnn,amee,amdd,amne,amnd,amed,xctr,yctr,irad)
    elseif kplot=2 then 
      call PLOTSMT(amnn,amee,amdd,amne,amnd,amed,xctr,yctr,irad)
    else
      call PLOTVMT(amnn,amee,amdd,amne,amnd,amed,xctr,yctr,irad)
    end if
  end if
  gosub PLOTS

  if isor = 3 or isorc=nsor then
    timax=sorpar(5,nsor)+(ntau(nsor)+1)*tau0 
    npd=fnamin(.04*kscrx,(kscrx*.38)/timax)  'changed from npd=20 PZ -94/5/18
    gosub PLOTSTF
    if isorc=nsor then 
      '--write a legend for line types
      np$="Nodal Planes:"
      ir=nrowchar-6 : ic=ncolchar-len(np$)
      call scprint(ir,ic,np$)
      if ktype<>2 then
        ir=ir+1
        if kokcol then color(kgrcol(kcolsyn))
        call scprint(ir,ic,string$(3,196)+" DC")
      end if
      if ktype<>1 then
        ir=ir+1
        if kokcol then color(kgrcol(kcoldat))
        call scprint(ir,ic,chr$(196)+" "+chr$(196)+" MT")
      end if
      call scpause(nrowchar-1,icol): goto ENDFPS
    end if
    if fninquit(ktrue) then goto ENDFPS
    call scclear
  end if
next isorc

ENDFPS: call setscreen (0): return
'-------------------------------------------------------------------------
PLOTP:
call pltnodpl(strt,dipt,xctr,yctr,irad)
call m5auxpl (str,dip,rk,strt,dipt,raket)
call pltnodpl(strt,dipt,xctr,yctr,irad)
return
'-------------------------------------------------------------------------
PLOTS: '---------- subroutine to plot radiation patterns

if kplot=1 then v=alph2 else v=beta2
if iok then jlo=1 else jlo=2
jhi=2
if ktype=1 then jhi=1 else if ktype=2 then jlo=2
for jtype=jlo to jhi  'for double couple and moment tensor radiation
  if kplot=3 and not jtype=2 then iangstep=3 else iangstep=10
  if jtype=1 then
    sd=sin(dip): cd=cos(dip): s2d=sin(2*dip): c2d=cos(2*dip)
    ss=sin(rk): cs=cos(rk)
    cdtcs=cd*cs : c2dtss=c2d*ss : sdtcs=sd*cs : s2dtss=s2d*ss
    kc=kgrcol(kcolsyn): ianglo=1
  else
    kc=kgrcol(kcoldat): ianglo=1+iangstep/2.
  end if

  for iang=ianglo to 89 step iangstep
    if fninquit(kfalse) then return
    ang=fnrad(iang)
    dthet=iangstep/(sqr2*sin(ang/2))
    if jtype=1 then
      rp=sin(ang)/v : rp2=rp*rp: vslo2=1/(v*v)-rp2 : vslo=sqr(vslo2)
      vslo2mrp2=vslo2-rp2 : rptvslo=rp*vslo
    else
      st=sin(ang) : ct=cos(ang) : s2t=sin(2*ang) : c2t=cos(2*ang)
      st2=st*st : ct2=ct*ct
      if kplot=1 then
        st2amnn=st2*amnn : st2amee=st2*amee : st2amne=st2*amne : ct2amdd=ct2*amdd
        s2tamnd=s2t*amnd : s2tamed=s2t*amed 
      elseif kplot=2 then
        stamnnyy=.5*st*(-amnn+amee) : stamne=st*amne
        ctamnd=ct*amnd : ctamed=ct*amed
      else
        s2tsamnn=.5*s2t*amnn : s2tsamee=.5*s2t*amee : s2tsamdd=.5*s2t*amdd 
        s2tsamne=.5*s2t*amne : c2tsamnd=c2t*amnd : c2tsamed=c2t*amed
      end if
    end if

    for az= 0 to 360 step dthet
      a=fnrad(az): if jtype=1 then a=a-str
      sa=sin(a): ca=cos(a): s2a=sin(2*a): c2a=cos(2*a)

      if jtype=1 then
        if kplot<>2 then
          a1=sdtcs*s2a-s2dtss*sa*sa: a2=s2dtss : a3=-cdtcs*ca+c2dtss*sa
          if kplot=1 then                    '--- sign of P down
            e1=a1*rp2+a2*vslo2 : e2=2*a3*rptvslo
            isgn=sgn(e1+e2)
          else                               '---- sign  of SV up
            e3=vslo2mrp2*a3 : e4=rptvslo*(a2-a1)
            isgn=sgn(-(e3+e4)) 
          end if
        else                                 '--- sign of SH down
          a2=cdtcs*sa+c2dtss*ca : a1=sdtcs*c2a-.5*s2dtss*s2a
          isgn=sgn(a1*rp + a2*vslo)
        endif
      else           'plot moment tensor rad pattern
        if kplot<>2 then                     '--- sign of P down
          if kplot=1 then
            e1=st2amnn*ca*ca+st2amee*sa*sa+st2amne*s2a+ct2amdd
            e2=s2tamnd*ca+s2tamed*sa
            isgn=sgn(e1+e2)
          else                               '---- sign  of SV up
            e3=s2tsamnn*ca*ca+s2tsamee*sa*sa-s2tsamdd+s2tsamne*s2a
            e4=c2tsamnd*ca+c2tsamed*sa
            isgn=sgn(-(-e3+e4)) 
          end if
        else                                 '--- sign of SH down
          e1=s2a*stamnnyy+stamne*c2a
          e2=-ctamnd*sa+ctamed*ca
          isgn=sgn(e1+e2)
        end if
      end if

      if isgn>1e-5 then 
        azrad=fnrad(az)
        call eaproj(irad,ang,azrad,x,y) : y=-y
        call scpoint(fnix(x),fniy(y),kc)
      end if

next az, iang, jtype
return

'========================================================================
sub PLTNODPL (strt,dipt,xctr,yctr,irad) static
'-------------- plot P-wave nodal plane
shared scratio,kgrcol()

dipt=dipt-fnrad(.1)
x=irad*cos(strt-piihalf)
y=irad*sin(strt-piihalf)
call scpoint(fnix(x),fniy(y),kgrcol(kcolsyn))

for i=0 to 180 step 10
  d1=fnrad(i)+strt
  p=piihalf-atn(tan(dipt)*sin(fnrad(i)))
  call eaproj(irad,p,d1,x,y)
  call scline1 (fnix(x),fniy(-y),kgrcol(kcolsyn),kfalse)
next
end sub
'========================================================================
sub PLOTPMT(amnn,amee,amdd,amne,amnd,amed,xctr,yctr,irad) static
'-------------- plot P-wave nodal plane for MT solution
shared scratio,kgrcol()
redim ix(2),iy(2),ixhld(2),iyhld(2),ndplt(2),newln(2),phi(5)

'-- initialize variables
azjincmax=5 : azjincmin=.5

'-- identify azimuths (phi) with radially directed nodal plane segments
phitest=amne*amne-amee*amnn
if phitest>=0 and abs(amee)>1e-5 then
  phitest=sqr(phitest)
  k=1
  for i=1 to -1 step -2
    phi(k)=fndeg(atn((-amne+i*phitest)/amee))
    phi(k+1)=180.+phi(k)
    for j=0 to 1 'make positive
      if phi(k+j)<0 then phi(k+j)=360.+phi(k+j)
    next j
    k=k+2
  next i
  for i=1 to 3  'sort in increasing order
    for j=i+1 to 4
      if phi(i)>phi(j) then SWAP phi(i),phi(j)
  next j,i
end if
azj=phi(1)
phi(5)=azj+360.
k=1  'phitest index
do
  az=fnrad(azj) : ca=cos(az) : sa=sin(az)
  a=amnn*ca*ca+amee*sa*sa+amne*sin(2*az)
  b=2*(amnd*ca+amed*sa)
  c=amdd
  c2=b*b-4*a*c
  if c2>0 or (c2=0 and abs(b)>=1e-7) then 
    sqrc=sqr(c2)
    imin=1 : imax=2
    if abs(a)<1e-7 then
      if b>0 then imax=1 else imin=2
    end if
    for i=imin to imax
      newln(i)=kfalse
      if i=1 then ix=1 else ix=-1
      if abs(a)>1e-7 then 
        p=(-b+ix*sqrc)/(2*a) 
        p=ATN(p)
      else 
        p=ATN(-c/b)
      end if
      if fninlim(p,0,piihalf) then
        newln(i)=(not ndplt(i))
        call eaproj(irad,p,az,x,y)
        ix(i)=fnix(x): iy(i)=fniy(-y)
        ndplt(i)=ktrue
      else 
        ndplt(i)=kfalse
      end if
    next i
    if azj>phi(1) then
'if two roots just found, connect them
      if newln(1) and newln(2) then
        ixhld(2)=ix(1): iyhld(2)=iy(2)
        newln(2)=kfalse
'if one root just appeared, connect to circle
      elseif newln(1) then
        p=piihalf
        call eaproj(irad,p,az,x,y)
        ixhld(1)=fnix(x) : iyhld(1)=fniy(-y)
        newln(1)=kfalse
      elseif newln(2) then
        p=piihalf
        call eaproj(irad,p,az,x,y)
        ixhld(2)=fnix(x) : iyhld(2)=fniy(-y)
        newln(2)=kfalse
      end if
    end if
    for i=imin to imax
      if ndplt(i) then
        if not newln(i) then_
           call scdash2 (ixhld(i),iyhld(i),ix(i),iy(i),kgrcol(kcoldat),kfalse)
        ixhld(i)=ix(i):iyhld(i)=iy(i)
      end if
    next
  else 
    ndplt(1)=kfalse:ndplt(2)=kfalse
  end if
'-- decide on a step size that will get closest to the outer radius
  azjinc=azjincmin
  if azj>=phi(k+1)-.1 then 
    k=k+1
  elseif azj-phi(k)>=azjincmax and phi(k+1)-azj>=azjincmax then
    azjinc=azjincmax
  elseif phi(k+1)-azj<=azjincmin then
    azjinc=phi(k+1)-azj
  end if
  if abs(c2)<.04 then azjinc=fnamin(azjinc,azjincmin)
  azj=azj+azjinc
loop until azj>phi(5)+.1
end sub
'========================================================================
sub PLOTSMT(amnn,amee,amdd,amne,amnd,amed,xctr,yctr,irad) static
'-------------- plot SH-wave nodal plane for MT solution
shared scratio,kgrcol()
redim ixhld(2),iyhld(2),ndplt(2)
azjincmax=3 : azjincmin=1
'-- initialize variables
azj=0
azjinclast=10
do
  az=fnrad(azj)
  a=sin(2*az)/2*(-amnn+amee)+amne*cos(2*az)
  b=-amnd*sin(az)+amed*cos(az)
  plast=piihalf
  p=piihalf
  i=1
  if abs(a)>1e-7 then 
    p=-b/a 
    p=ATN(p)
  else 
    p=-sgn(b)*piihalf
  end if
  if fninlim(p,0,piihalf) then
    call eaproj(irad,p,az,x,y): ix=fnix(x): iy=fniy(-y)
    if ndplt(i) then_
            call scdash2 (ixhld(i),iyhld(i),ix,iy,kgrcol(kcoldat),kfalse)
    ixhld(i)=ix:iyhld(i)=iy
    ndplt(i)=ktrue
    plast=fnamin(plast,p)
  else 
    ndplt(i)=kfalse
  end if
'-- decide on a step size that will get closest to the outer radius
  dp=fnamax(abs(plast-p),1e-5)
  azjinc=fnamin(fnamin((piihalf-p)*azjinclast/dp-.15,azjincmax),360-azj)
  azjinc=fnamax(azjinc,azjincmin)
  azjinclast=azjinc
  azj=azj+azjinc
loop until azj>360
end sub
'========================================================================
sub PLOTVMT(amnn,amee,amdd,amne,amnd,amed,xctr,yctr,irad) static
'-------------- plot SV-wave nodal plane for MT solution
shared scratio,kgrcol()
redim ixhld(2),iyhld(2),ndplt(2)
azjincmax=10 : azjincmin=1
'-- initialize variables
azj=0
azjinclast=10
clast=1
do
  az=fnrad(azj) : ca=cos(az) : sa=sin(az)
  a=-(amnd*ca+amed*sa)
  b=amnn*ca*ca+amee*sa*sa-amdd+amne*sin(2*az)
  c=-a
  c2=b*b-4*a*c
  plast=piihalf
  p=piihalf
  if c2>0 or (c2=0 and abs(b)>=1e-7) then 
    sqrc=sqr(c2)
    imin=1 : imax=2
    if abs(a)<1e-7 then
      if b>0 then imax=1 else imin=2
    end if
    for i=imin to imax
      if i=1 then ix=1 else ix=-1
      if abs(a)>1e-7 then 
        p=(-b+ix*sqrc)/(2*a) 
        p=ATN(-p) 'note sign reversal
      else 
        p=ATN(c/b) 'note sign reversal
      end if

      if fninlim(p,0,piihalf) then
        call eaproj(irad,p,az,x,y): ix=fnix(x): iy=fniy(-y)
        if ndplt(i) then_
            call scdash2 (ixhld(i),iyhld(i),ix,iy,kgrcol(kcoldat),kfalse)
        ixhld(i)=ix:iyhld(i)=iy
        ndplt(i)=ktrue
        plast=fnamin(plast,p)
      else 
        ndplt(i)=kfalse
      end if
    next i
  else
    ndplt(1)=kfalse : ndplt(2)=kfalse
  end if
'-- decide on a step size that will get closest to the outer radius
  dc=fnamax(abs(clast-c2),1e-5):dp=fnamax(abs(plast-p),1e-5)
  azjinc=fnamin(c2*azjinclast/dc,(piihalf-p)*azjinclast/dp)-.1
  azjinc=fnamin(fnamin(azjinc,azjincmax),360-azj)
  azjinc=fnamax(azjinc,azjincmin)
  clast=c2:azjinclast=azjinc
  azj=azj+azjinc
loop until azj>360
end sub
'========================================================================
PLOTSRAD:'-- subroutine to plot SH nodal pattern
n=0 : sd=sin(dip): cd=cos(dip): s2d=sin(2.*dip):
c2d=cos(2.*dip) : ss=sin(rk): cs=cos(rk)

'-- singular cases:
if dip<1E-7 then  'dip=0
  ain=90 : az=fndeg(str-rk) : gosub PROJNET : ixt1=ix : iyt1=iy
  ain=-ain : gosub PROJNET
  call scline2(ixt1,iyt1,ix,iy,kgrcol(kcolsyn),kfalse)
  go to 3180
elseif s2d<1E-7 and ABS(cs)<1E-7 then  'dip=0 or 90 deg and rk=90 or 270 deg
  ain=90 : az=fndeg(str)+90 : gosub PROJNET : ixt1=ix : iyt1=iy
  ain=-ain : gosub PROJNET
  call scline2(ixt1,iyt1,ix,iy,kgrcol(kcolsyn),kfalse)
  go to 3180
elseif cd<1E-7 and ABS(ss)<1E-7 then  'dip =90 deg and rk=0 or 180 deg
  ain=90 : az=fndeg(str)+45 : gosub PROJNET : ixt1=ix : iyt1=iy
  ain=-ain : gosub PROJNET
  call scline2(ixt1,iyt1,ix,iy,kgrcol(kcolsyn),kfalse)
  az=fndeg(str)-45 : gosub PROJNET : ixt1=ix : iyt1=iy
  ain=-ain : gosub PROJNET
  call scline2(ixt1,iyt1,ix,iy,kgrcol(kcolsyn),kfalse)
  go to 3180
elseif ABS(c2d)<1E-7 and ABS(cs)<1E-7 then 'dip=45 deg and rk=90 or 270 deg
  ain=90 : az=fndeg(str) : gosub PROJNET : ixt1=ix : iyt1=iy
  ain=-ain : gosub PROJNET
  call scline2(ixt1,iyt1,ix,iy,kgrcol(kcolsyn),kfalse)
  az=fndeg(str)+90 : gosub PROJNET : ixt1=ix : iyt1=iy
  ain=-ain : gosub PROJNET
  call scline2(ixt1,iyt1,ix,iy,kgrcol(kcolsyn),kfalse)
  go to 3180
end if
    
notyet=ktrue : sgnold=999 : n=0
azstep=1.0 : dz=-azstep

STARTSH: '--- start loop for radiation pattern
dz=dz+azstep : if dz > 180. then 3180
az=dz : gosub AMPLSH
if ABS(a1)<1E-7 then  '--- special case : ain = 90
  if notyet then ain=90. : gosub PROJNET : goto 3120
  if a2 <>0 then ain=90.*sgnold : gosub PROJNET : goto 3120
  ain=-90. : gosub PROJNET : ixx=ix : iyy=iy
  ain=90. : gosub PROJNET
  call scline2(ixx,iyy,ix,iy,kgrcol(kcolsyn),kfalse)
  n=0 
  goto 3130
end if
     
ain=fndeg(atn(-a2/a1))
gosub PROJNET
fcn=-a2/a1
sdfdaz=sgn(((-2*sd*cs*s2a-s2d*ss*c2a)*a2-(cd*cs*ca-c2d*ss*sa)*a1)/a1/a1)
sdshdrp=sgn(a1*cos(fnrad(ain))-a2*sin(fnrad(ain)))
if  notyet then 3115           '----  first time
if sdfdaz=sgn(ain-ainold) then 3120        '--- check for discontinuity
if sdfdazold<>sdfdaz and sdshdrpold=sdshdrp then 3120
ixx1=ix : iyy1=iy : ain=90.*sdfdaz : gosub PROJNET
n=n+1 
if n=1 then call scpoint(ix,iy,kgrcol(kcolsyn)) else call scline1(ix,iy,kgrcol(kcolsyn),kfalse)
n=0
ain=-ain : gosub PROJNET
call scline2(ix,iy,ixx1,iyy1,kgrcol(kcolsyn),kfalse)
n=2
goto 3130

3115 '--- first point, just move
n=n+1 : call scpoint(ix,iy,kgrcol(kcolsyn))
ix1=ix : iy1=iy
goto 3130

3120 n=n+1 : if n=1 then call scpoint(ix,iy,kgrcol(kcolsyn)) else_
call scline1(ix,iy,kgrcol(kcolsyn),kfalse)

3130 sgnold=sgn(ain) : ainold=ain : sdfdazold=sdfdaz : sdshdrpold=sdshdrp

notyet=kfalse
goto STARTSH

3180 return
'=================================================================
sub EAPROJ ( irad, toang, az, x, y) static
r2=irad*sqr2*sin(toang/2)
x= +r2*cos(az-piihalf)
y= -r2*sin(az-piihalf)
end sub
'=================================================================
PROJNET: ' find x,y for projection
da1=fnrad(ain) : da2=fnrad(az)
call eaproj(irad,da1,da2,x,y) : iy=fniy(-y): ix=fnix(x)
return
'-----------------------------------------------------------------------------
AMPLSH: '---- amplitudes for SH
a=fnrad(dz)-str : sa=sin(a): ca=cos(a): s2a=sin(2.*a): c2a=cos(2.*a)
a2=cd*cs*sa+c2d*ss*ca : a1=sd*cs*c2a-.5*s2d*ss*s2a
f=a1*rp + a2*vsb
return
'========================================================================
sub hrminsec ( sectime, ihr, imin, secs ) static
secs=sectime
ihr = int(secs/3600.)
secs=secs - ihr*3600.0
imin=int(secs/60.)
secs=secs - imin*60.
end sub
'========================================================================
rem $INCLUDE: 'm5auxpl'
rem $INCLUDE: 'm5stf'
rem $INCLUDE: 'm5ptbax'
' rem $include: 'm5eigen' 'now linked separately
'========================================================================

