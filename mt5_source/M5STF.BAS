sub m5stf (inflag,tau,isor,ntau(1),atau(2),tam(1), tcent, t95) static
 
'--- Subroutine to normalize source function and calculate duration
'      INVFLAG flags whether the moment is updated or not
'      Output:  ATAU(i) - normalized STF, and TAM - moment corrected
'               TCENT - centroid time, T95- 95% duration
'
  shared ntaumax

  imflag=inflag
  if ntau(isor)=0 then ntau(isor)=1
  ntau=ntau(isor)
'-- set rest of atau=0
  if ntau < ntaumax then for i=ntau+1 to ntaumax : atau(i,isor)=0 : next

'-- normalize ATAU to total amplitude, and correct moment by that amplitude
stf1:  
   x=0 : for i=1 to ntau : x=x+atau(i,isor) : next
   if abs(x)>1e-6 then goto stf2
     atau(1,isor)=1 : print "    NOTE: STF=0; resetting ATAU(1,";fnstrnum$(isor)+")";
     imflag=kfalse
     goto stf1
stf2:  
   for i=1 to ntau : a=atau(i,isor) : atau(i,isor)=a/x: next  'set area =1
   if imflag then tam(isor)=tam(isor)*x

'-- calculate half duration and 95% duration

  acent=0 
  for i=1 to ntau : acent=acent+fnamax(0., atau(i,isor) )*i : next
  tcent=acent*tau         '---centroid time

'--- 95% duration calculated by second moment method
'         (2-sigma for gaussian STF):
'    t95 = 2 * sqr ( (second moment)/(zeroth moment))  

   sum=0
    for i=1 to ntau
      sum=sum+fnamax(0., atau(i,isor) )*(i*i+1/6)
    next i
     sss= (sum - (tcent/tau)^2)
     sss= fnamax(0.0, sss) 
     t95= 4 * tau * sqr (sss)

end sub
'-------------------------------------------------------------------------
PLOTSTF: '--------- plot source time function

    amax=0
6210 for isor=1 to nsor
       for i=1 to ntau(isor)
        amax=fnamax(amax, abs(tam(isor)*atau(i,isor)))
   next i,isor

   for isor =1 to nsor
     ntau=ntau(isor)
     if ntau<ntaumax then atau(ntau+1,isor)=0
     toff=sorpar(5,isor)  'offset for event

6240 ix=10+npd*toff : iy=.95*kscry : call scpoint (ix,iy,kgrcol(kcolbak))
  for i=1 to ntau
   it=10+(i*tau0+toff)*npd
   ia=iy-(kscry/8)*atau(i,isor)*tam(isor)/amax
   call scline1(it,ia,kgrcol(kcolbak),kfalse)
  next
   it=10+((ntau+1)*tau0+toff)*npd
   ia=0.95*kscry
   call scline1(it,ia,kgrcol(kcolbak),kfalse)
 next isor

   ia1=kscry : ia2=.95*kscry : imax=5*int((timax+4.99)/5)
   call scline2 (10,ia1,10,ia2,kgrcol(kcolbak),kfalse)
   for i=5 to imax step 5
     it=i*npd+10 : call scline1 (it,ia2,kgrcol(kcolbak),kfalse)
     call scline2 (it,ia1,it,ia2,kgrcol(kcolbak),kfalse)
   next
     l$="Tics: 5s"
     if ispflag and ichoice=115 then l$=l$+"-LP 1s-SP"
     it= fnamin(int(it/kscrx*ncolchar+2),ncolchar-1-len(l$))
     call scprint(nrowchar-1,it,l$+"   "):icol=0
     call scprintc(nrowchar-1,icol,"Press Q to quit")
return
