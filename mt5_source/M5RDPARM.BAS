M5RDPARM:
'
'---- Read file with model parameters
print:print: d$=ucase$(fhdr$)+".PRM":irow=0
call scprint(irow,1,"Default parameter file is "+d$)
call scprintc(irow+2,2," Enter name of parameter file (Esc or Q if none)")
if inputn(irow+3,2," --",d$,kfalse,ktrue,kfalse,"a",40)=27 then return
if fniquit(d$) then return
print:


if len(d$)=0 then d$=fhdr$+".prm"

pf$=fnaddext$(d$,"prm")
kfnum=1


if ifopen("I",kfnum,pf$,0,-1)=0 then goto M5RDPARM
if lcase$(fhdr$)= "mt5test" then fhdr$=fnsubext$(pf$)

line input #kfnum, a$
npps=abs(fnvalue(1,8))
npix=fnvalue(9,8) : nsor=fnvalue(25,8)
nlyr=fnvalue(41,8) : iplform=fnamax(1,fnvalue(49,8))
nphtmp=fnvalue(57,8) : if nphtmp=0 then nphtmp=7
mtplform=fnamax(1,fnvalue(65,8))

line input #kfnum, a$
tau0=fnvalue(1,8) : dep=fnvalue(9,8) : wthk=fnvalue(17,8)
mktau=ktrue

line input #kfnum, a$
theta=fnvalue(1,8)

'--- flags
line input #kfnum, a$
ipcflag=fnvalue(1,8) : invflag=fnvalue(9,8)
istflag=fnvalue(17,8): normamp=fnvalue(25,8)
mtflag=fnvalue(33,8): nrevcom=fnvalue(41,8)
if invflag then varcor=1.0e8 else varcor=1.0


'--- rupture
line input #kfnum, a$
vrupt=fnvalue(1,8) : azrupt=fnvalue(9,8) : plrupt=fnvalue(17,8)

for i=1 to nphtmp
  line input #kfnum, a$
  phscl(i)=fnvalue(6,8): phrlen(i)=fnvalue(14,8)
  phtime(i)=fnvalue(22,8): phdt(i)=fnvalue(30,8)
  phtstar(i)=fnvalue(38,8): phmag(i)=fnvalue(46,8)
next

'--- source structure
for j=1 to nlyr
  line input #kfnum, a$
  ztop(j)=fnvalue(1,8)  : alph(j)=fnvalue(9,8)
  beta(j)=fnvalue(17,8) : rho(j)=fnvalue(25,8)
  cthk=wthk-(j>1)*(ztop(2)-ztop(1))
next

'--- Station parameters

'--- no. stations in data file
input #kfnum, nsta
okprm=ktrue

if mkdat then stacode$=""

for j=1 to nsta
  line input #kfnum, a$
  s$=left$(a$,4)
  if mkdat then
    stacode$=stacode$+s$
  else
    okprm=(ucase$(s$)=left$(fnstanm$(j),4))
  end if
  iflag(j)=fnvalue(5,4)
  ipf1(j)=fnvalue(9,5)
  ipf2(j)=fnvalue(14,5)
  ninv(j)=fnvalue(19,5)
  ispos(j,1)=fnvalue(24,5)
  ipol(j)=fnvalue(29,5)
  wt(j) = fnvalue(34,7)
  tstar(j)=fnvalue(41,7)
  wt0(j) = fnvalue(48,7)
  ispos(j,2)=fnvalue(55,5)
  ixfine(j)=fnvalue(65,5)
  iyfine(j)=fnvalue(70,5)
  refmag(j)=fnvalue(75,8)
  iphn=fnvalue(84,1)
  instn=fnvalue(85,1)
  if instn=0 then iphn=0
  if mkdat then
    iph(j)=iphn
    inst(j)=instn
  elseif instn>0 then
    okprm=((iph(j)=iphn) and (inst(j)=instn))
  end if
  if not okprm then
    call schit(".PRM and .INV files do not match on station #"+str$(j))
    close #kfnum : return
  end if
next j

for isor=1 to nsor
  line input #kfnum, a$
  for i=1 to 7 : sorpar(i,isor)=fnvalue(10*i-9,10) : next
  sorpar(1,isor) = fnrad(sorpar(1,isor))
  sorpar(2,isor) = fnrad(sorpar(2,isor))
  sorpar(3,isor) = fnrad(sorpar(3,isor))
  sorpar(4,isor) = sorpar(4,isor)-cthk
  sorpar(7,isor) = fnrad(sorpar(7,isor))
  tamiso(isor) = fnvalue(71,10)
  tam(isor)=fnvalue(121,11)
  for i=8 to 13 : sorpar(i,isor)=fnvalue(10*i+52,10) : next

  input #kfnum, ntau(isor)                         '--- no. STF elems

'-- STF amplitudes
  line input #kfnum, a$
  for j=1 to ntau(isor) : atau(j,isor)=fnvalue(7*j-6,7) : next

'-- adjustment limiters
  line input #kfnum, a$ ' : print a$
  for j=1 to 7 : damp(j,isor)=fnvalue(6*j-5,6) : next
  stfdamp(isor)=fnvalue(43,6)
  if len(rtrim$(a$))<49 then
    for j=8 to 13 : damp(j,isor)=1. : next
  else
    for j=8 to 13 : damp(j,isor)=fnvalue(6*j+1,6) : next
  end if
next isor

for ista=1 to nsta
  if abs(ipol(ista))<>1 then ipol(ista)=1
next

for j=1 to nlyr
  if j=1 then zt=ztop(j)
  ztop(j)=ztop(j)-zt
next j

line input #kfnum, a$
title$=a$

line input #kfnum, a$
kletter=fnvalue(1,3) : kplabel=fnvalue(4,3) :labpos=fnvalue(7,3)
kpmech=fnvalue(10,3) :isposit=fnvalue(13,3) :iampos=fnvalue(16,3)
ifocpos(1)=fnvalue(19,3): kfoc(1)=fnvalue(22,3)
ifocpos(2)=fnvalue(25,3): kfoc(2)=fnvalue(28,3)
labcode=fnvalue(31,3) : igraph=fnvalue(34,3)
kpenhp=(igraph<0): igraph=fnxlimit(abs(igraph),1,5)

line input #kfnum, a$
ifocx0(1)=fnvalue(1,4) : ifocy0(1)=fnvalue(5,4) : ifocrad0(1)=fnvalue(9,4)
ifocx0(2)=fnvalue(13,4) : ifocy0(2)=fnvalue(17,4) : ifocrad0(2)=fnvalue(21,4)
spcsp=fnvalue(25,8) :spclp=fnvalue(33,8)
boxlen0(2)=fnvalue(41,8): boxhigh0(2)=fnvalue(49,8):

if not nodat then gosub RESETFLG

mkgrn=ktrue
iparflg = ktrue

close #kfnum
return
'=====
RESETFLG:
for i=1 to 7 : kphflag(i)=0 : next
for i=1 to nsta
  if iflag(i) then
    kphflag(kphtype(i))=ktrue
  end if
next
'-- set all phase flags
kplot(1) = (kphflag(1) or kphflag(3) or kphflag(4) or kphflag(6) or kphflag(7))
kplot(2) = (kphflag(2) or kphflag(5))
ilpflag  = (kphflag(1) or kphflag(2) or kphflag(4) or kphflag(5) or kphflag(7))
ispflag  = (kphflag(3) or kphflag(6))
igdsnflag = (kphflag(4) or kphflag(5) or kphflag(6) or kphflag(7))
return
'=====
