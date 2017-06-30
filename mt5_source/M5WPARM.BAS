  '-------------------------------------------------------------------------
M5WPARM: '---- Write file with model and input parameters
if pf$="" then pf$=fhdr$
irow=0
call scprintc(irow,1,"Enter prefix for .PRM file to write (I/O file="+fhdr$+")")
if inputn(irow+1,2,"--",pf$,ktrue,ktrue,ktrue,"a",40)=27 then return
if len(pf$)=0 then return
pf$=fnsubext$(pf$)
if keepold(pf$+".prm") then return

8510 open "o",1,pf$+".prm"

x=0 :ix=0 '-- for unused spots

'--- various controls
print #1, using fniform$(8); npps; npix; ix; nsor; ix; nlyr; iplform; ix; mtplform
print #1, using fnfform$(8,2); tau0; sorpar(4,1)+cthk; wthk
print #1, using fnfform$(8,2); theta

'--- flags
print #1, using fniform$(8); ipcflag; invflag; istflag; normamp; mtflag; nrevcom

'--- rupture parameters
print #1, using fnfform$(8,2); vrupt; azrupt; plrupt

'--- phase types
for i=1 to nphtype
  print #1, using fnaform$(5); fnkphlab$(i);
  print #1, using fnfform$(8,2); phscl(i);phrlen(i);phtime(i);phdt(i);phtstar(i);
  print #1, using fniform$(8); phmag(i)
next

'--- structure model
for j=1 to nlyr
  print #1, using fnfform$(8,2); ztop(j);alph(j);beta(j);rho(j)
next

'--- station parameters

print #1, using fniform$(5); nsta
for j=1 to nsta
  print #1, using fnaform$(4); fnstanm$(j);
  print #1, using fniform$(4); iflag(j);
  print #1, using fniform$(5); ipf1(j); ipf2(j); ninv(j); ispos(j,1); ipol(j);
  print #1, using fnfform$(7,2); wt(j);tstar(j);wt0(j);
  print #1, using fniform$(5); ispos(j,2), ix; ixfine(j); iyfine(j);
  print #1, using fniform$(8); refmag(j);
  print #1, " ";
  print #1, using fniform$(1); iph(j); inst(j)
next j

'--- source parameters

for isor=1 to nsor
  str=fndeg(sorpar(1,isor))
  dip=fndeg(sorpar(2,isor))
  rk=fndeg(sorpar(3,isor))
  dep=sorpar(4,isor)+cthk
  deltim=sorpar(5,isor)
  delpos=sorpar(6,isor)
  delphi=fndeg(sorpar(7,isor))

  print #1, using fnfform$(10,2); str;dip;rk;dep;deltim;delpos;delphi;
  print #1, using fndform$(6,4); tamiso(isor);
  for i=1 to 4 : print #1, using fnfform$(10,2); x; : next
  print #1, using fndform$(7,5); tam(isor);
'-- print moment tensor normalized elements
  for i=8 to 13 : print #1, using fnfform$(10,4); sorpar(i,isor); : next
  print #1,

  write #1, ntau(isor)                         '--- no. STF elems
  for j=1 to ntau(isor)
    print #1, using fnfform$(7,3); atau(j,isor);
  next j
  print #1,

  for j=1 to 7 : print #1, using fnfform$(6,2); abs(damp(j,isor)); : next
  print #1, using fnfform$(6,2); stfdamp(isor);
'-- print moment tensor adjustment multipliers
  for j=8 to 13 : print #1, using fnfform$(6,2); damp(j,isor); :next
  print #1,

next isor

'--- plotting controls
print #1, title$
print #1, using fniform$(3); kletter;kplabel;labpos;kpmech;isposit;iampos;ifocpos(1);kfoc(1);ifocpos(2);kfoc(2);labcode;igraph*(kpenhp*2+1)
print #1, using fniform$(4); ifocx0(1);ifocy0(1);ifocrad0(1);ifocx0(2);ifocy0(2);ifocrad0(2);
print #1, using fnfform$(8,2); spcsp;spclp;boxlen0(2);boxhigh0(2)
close 1
return
