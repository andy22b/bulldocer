' main MT5IT (MT5.00)
' performs waveform synthesis and inversion for MT5

DEFINT I-N: OPTION BASE 1

'following are defined in m5rapars
DECLARE sub m5rapars(dep!,delta!,iph%,alph2!,beta2!,rapar!,ptoa!)
DECLARE sub rayp(i%,j%,iph%,rpjb!)

rem $INCLUDE: 'M5COMMON'
rem $INCLUDE: 'M5FUNCTN'

if not imt5 then call notalone

dim z2(nparmax),diag(nparmax), sumwtph(7)
dim atau2(ntaumax,nsormax), taudamp(ntaumax,nsormax)
dim f(3),df(3,6),vslow(2),pdly0(2),pdly(3),sorparel(10,nsormax)
dim tmp(2*maxnes) 'changed from 2*nesmax    PZ-5/10/94
dim nw(7)

'-- internal seismogram scaling factor UCOR so that velocities are in km/s,
'   densities in g/cm3, amplitudes will be in cm
ucor=1.0e-13
blank$=string$(80,32)
npage=2

'--- iterate
call scclear ': call scprint (1,1,"") : INPUT "Number of iterations (1 to 8) "; nit
if input1(1,1,"Number of iterations (1 to 8)",nit$,kfalse,ktrue,ktrue,numlist$(0,8))=27 then goto 2525
if nit$="0" then goto 2525 else nit=val(nit$)

if mtflag then    'set indexing parameters if moment tensor analysis
   nparmech=6 : noffmech=noffmt : na=3
else
   nparmech=3 : noffmech=0 : na=1
end if

redim DVAR(nit)
redim KPAR(4,nsormax,na), KTAU(ntaumax,nsormax,na), IEVFLAG(nsormax,na)
redim KMECH(nparmech,nsormax,na),npar(na)

'--- find number of parameters NPAR and assign places in parameter list
'    KPAR(i,isor,ia) = which free parameter is SORPAR(i+3,isor)
'    KTAU(i,isor,ia) = which free parameter is ATAU(i,isor)
'    KMECH(i,isor,ia) = which free parameter is SORPAR(i+noffmech,isor)
'                 for i=1 to 3,fault plane analysis
'                  or for i=1 to 6,moment tensor analysis
'                 for isor=1 to number of sources
'                 for ia=1,fault plane analysis
'                   or ia=2,tau analysis only
'                   or ia=3,moment tensor analysis w/out tau analysis

for ia=1 to na
  npar(ia)=0
  needmt=kfalse
  if ia<3 then needtau=kfalse
  for isor=1 to nsor
    nparsor=0
'-- mechanism
    if ia=1 or ia=3 then     'invert for mechanism,depth and rupture
      for i=1 to nparmech
        if damp(i+noffmech,isor)>0 then
          npar(ia)=npar(ia)+1
          kmech(i,isor,ia)=npar(ia)
          nparsor=nparsor+1
          needmt= ktrue
        end if
      next i
'-- depth and rupture multiple source zones
      for i=1 to 4
        if damp(i+3,isor)>0 then
          npar(ia)=npar(ia)+1
          kpar(i,isor,ia)=npar(ia)
          nparsor=nparsor+1
        end if
      next
    end if
'-- source time function
    if ia=1 or (ia=2 and ntau(isor)>1) then    'invert for source time function
      if stfdamp(isor)>0 then
        needtau= ktrue
        for i=1 to ntau(isor)
          npar(ia)=npar(ia)+1
          ktau(i,isor,ia)=npar(ia)
          nparsor=nparsor+1
        next
      end if
    end if

'-- flag the event if it has no parameters to be determined
    if nparsor=0 then ievflag(isor,ia)=kfalse else ievflag(isor,ia)=ktrue

  next isor
next ia

'--- check that sum weights > 0  and number of data > number of parameters
if nit>1 then
  swt=0.0 : ninv=0
  for ista=1 to nsta
    if iflag(ista) then swt=swt+wt(ista)
    if iflag(ista) and wt(ista)>0 then ninv=ninv+ninv(ista)
  next
  if swt=0 then call scprint (22,10,"*** Sum of weights is zero *** ") : goto 2520
  for ia=1 to na
    if ninv<npar(ia) then call scprint (22,10,"*** Number of data < number of parameters ") : goto 2520
  next
end if

'-- one screen is used for each subevent, the top 21 lines of the screen
'    is stored in character string LINER$
'--  lines 22-25 are common to all displays

'-- permanent
call scprintc(24,1,"OPTIONS:")
call scprompt((0),10,"Q-quit"):call scprompt((0),18,"#-display event #")
call scprompt((0),37,"H-help")
call scprint(22,1,"MESSAGE:")

'--clear all screens
redim liner$(nsor+1)
ncharm=2000

'--- write help screen to LINER$(1)
iscreen=1
liner$(iscreen)=string$(ncharm,32)
call scprint1(1,1,string$(80,"-"))
call scprint1(21,1,string$(80,"-"))
call scprint1(1,58," MT5    HELP SCREEN ")
call scprint1(21,37," Enter Source number for ITERATION SCREEN ")

call scprint1(3,1,"The iteration screen displays information about the model and statistics for ")
call scprint1(4,1,"each iteration. ")

call scprint1(5,1,"It - iteration number")
call scprint1(6,1,"Err - standard error for the parameter above it for the most recent model")
call scprint1(7,1,"UnwtdVAR - the unweighted variance of the amplitude residuals")
call scprint1(8,1,"WtdVAR   - the weighted variance of the amplitude residuals")
call scprint1(9,1,"R/D%  - the ratio in percent of the residual variance over the data variance")

call scprint1(11,1,"If there are multiple sources (subevents) one iteration screen is kept for ")
call scprint1(12,1,"each. Press the subevent number (1-6) to see the screen for that event.")

call scprint1(14,1,"The MESSAGE line will describe conditions that cause the iterating to ")
call scprint1(15,1,"terminate. Also if the variance reduction goes below the minimum, iteration ")
call scprint1(16,1,"will cease.")


'-- put labels on all screens
for isor= 1 to nsor
  iscreen=isor+1
  liner$(iscreen)=string$(ncharm,32)
  call scprint1 (1,1,"Source")
  call scprint1 (1,14,"iterations     parameters")
  call scprint1 (1,47,"==== MT5   ITERATION SCREEN ==")
  call scvint (1,10,3,nit)
  call scvint (1,25,3,npar(1))
  if mtflag then
    call scprint1 (2,1,"It   Mnn  Mee  Mdd  Mne  Mnd  Med  Dpt dTo Dist Az  Moment UnwtdVAR WtdVAR  R/D%")
  else
    call scprint1 (2,1,"It Strike   Dip  Rake Depth   dTo  Dist    Az   Moment  UnwtdVAR  WtdVAR   R/D%")
  end if
  call scprint1 (12, 1, "It   STF element amplitudes in percent of total")
next
if nit>1 then call scprint2(npage) ELSE call scprint(1,1,"Source"+str$(npage-1)+"   1 iteration")

' -- initialize flags and variable range parameters for loops
nvarp=0
ivarp=0
for isordp=1 to nsor
  for i=1 to 4
    if instr(vrp$(isordp),mid$("HOXY",i,1))>0 then
      if varinc(i,isordp)<>0 then
        sorpar(3+i,isordp)=varbeg(i,isordp)
        v=abs((varend(i,isordp)-varbeg(i,isordp))/varinc(i,isordp))
        j=int(v+1.99)
        if i=4 then 
          if instr(vrp$(isordp),"X")>0 then
            j=j+(abs(fnamin(varbeg(3,isordp),varend(3,isordp)))<.1)
          else
            j=j+(sorpar(6,isordp)<.1)
          end if
        end if
        nvarp=fnamax(1,nvarp)*fnamax(1,j)
      end if
    end if
  next i
next isordp
ndhdrdp=nvarp>0

'*************** ITERATION LOOP
519 if mtflag and needmt then
  if needtau then ias=2 else ias=3
  na=3
else 
  ias=1
  na=1
end if
if nvarp then 
  ivarp=ivarp+1
  call scprint (22,10,"On variable range step "+fnstrnum$(ivarp)+" of "+fnstrnum$(nvarp))
end if
for itno=1 to nit
'-- set timer
 time1=timer

'--- flag ITLAST for last iteration
 itlast = ((itno=nit) or (ivarstop))
'--- flag MKSTAT to enable statistics
 mkstat= (nit>1) or (nit=1 and krmsc)

'-- for 1 iteration do not redimension arrays or check bounds
 if mkstat then
  if not fninlim(npar(1),1,nparmax) then _
    call scprint (22,10,"*** Invalid number of parameters") : goto 2520
 end if

 depth1=sorpar(4,1)
 nrec=0: wres2=0: res2=0: sumwt=0: neq2=0: npph=0
 sumdswt=0 : sumdos=0

 erase datvar, resvar, sumwtph

'******************* INNER LOOP
 for ia=ias to na
  redim b(npar(ia),npar(ia)), b2(npar(ia),npar(ia)), z(npar(ia))

'--- flag IAFIRST for first inner loop (calc STF elements)
  iafirst= (ia=ias)
'--- flag IALAST for last inner loop of last iteration (calc everything else)
  ialast = (ia=na)

'******************* STATION LOOP

  for ista=1 to nsta

    if (not iflag(ista)) then goto 2090  'skip station

'--- Set bounding indices:
'      NPH= # of upgoing phases
'      IMNO = # pts in Inst*Q*Tau curve (SP), for this station.
'      NINV = # pts to use in inversion from seismogram, starting at IPF2(ista)
'      NES  = # pts total in seismograms
'      NESP = # pts of ES to add; =NINV (or =NES on last iteration)
'      NSYN = # pts of direct phase to add
'
    nes=nes(ista)
    noninv=(ninv(ista)<=1)
    iph=iph(ista)
    kpwave=fnpwave(ista)
    xmag= refmag(ista)   ''phmag(kphtype(ista))
    if kpwave then nph=2 else nph=1
    inst=inst(ista)
    imno=inpts(ista)
    ninv=ninv(ista) : n1=ndat(ista)-ipf2(ista) : ninv=fnamin(ninv,n1)
    neq2=neq2+ninv
    nesp=nes                    'number of ES points to read in
    nsyn=fnamin(imno,nes)       'number of points in direct phase to add
    nphes=iphes(ista)
    isprec1 = index(ista)       'first record in SP array
    dt=dt(ista)
    nv=ninv+1
    nph1=nph+1
    wt=wt(ista)

'--- except for last iteration, only do NINV pts in synthetics
    if (not itlast) then
      if (wt(ista)=0) or noninv or nodat(ista) then goto 2090  'skip this station
      nesp=fnamin(ninv+1,nes) : nsyn=fnamin(ninv+1,imno)  'cut to NINV if smaller
    end if

    np=fnamax(1,npar(ia))
    redim d(nesp,np), es3(nesp,nph1), es2(nesp), r(nv)
'    erase ss
    redim ss(2*maxnes)

'-- get the observed seismograms starting at point IPF2
    ipof=ipf2(ista)-1
    if not noninv and mkstat then
      call m5synio(1,0,ista,ninv,ipof,odat(),dext(ista),datsq)
    end if

'-- calculate ray parameter at depth of first event
    d=depth1+cthk
    call m5rapars (d, del(ista),iph,alph2,beta2,rapar(ista),ptoa(ista))
    az=fnrad(azd(ista))
    rp=rapar(ista)                         'ray parameter for phase
    if kpwave then rpp=rp else rpp=rp/vps  'P ray parameter

'-- Get ACOR, geometrical spreading and receiver response
'--   First, at reference dist DEL0, get ref. correction GSP0
    call geomspr(del0,rpp,rho2,rhor,alph2,alphr,gspr)
    call recvresp (kpwave,rpp,alphr,betar,acor)
    gsp0=acor*gspr

'------ Now, at actual delta
    delta=del(ista)
    call geomspr(delta,rpp,rho2,rhor,alph2,alphr,gspr)

'--- GCORSTA is amplitude correction to data to move station to delta=DEL0
    gcorsta=gsp0/(acor*gspr)
    for i=1 to ninv : odat(i)=odat(i)*gcorsta : next
    gcor(ista)=gcorsta

'******************* SOURCE LOOP

    for isor = 1 to nsor

      lastsor= (isor=nsor)
      dephs =sorpar(4,isor)
      delort=sorpar(5,isor)
      delpos=sorpar(6,isor)
      delphi=sorpar(7,isor)

'-- Calculate radiation pattern F(3) and its derivatives DF(3,6) w.r.t. mechanism
'  Indexing for F(i) and DF(i,j):  i=1,2,3 for P down,P up,SV up; or SH down,SH up;
'   j=1,2,3 for deriv's w.r.t strike, dip, rake; or, j=1 to 6 for deriv's w.r.t.
'   Mnn,Mee,Mdd,Mne,Mnd,Med where n=North,e=East,d=Down

      if mtflag then
        isgn = (tam(isor)<0)*2+1 'occassionally need to correct for neg moment
        tam(isor)=isgn*tam(isor)
        amnn = isgn*sorpar(1+noffmt,isor)
        amee = isgn*sorpar(2+noffmt,isor)
        amdd = isgn*sorpar(3+noffmt,isor)
        amne = isgn*sorpar(4+noffmt,isor)
        amnd = isgn*sorpar(5+noffmt,isor)
        amed = isgn*sorpar(6+noffmt,isor)
        call m5radmt (mkstat,kpwave,amnn,amee,amdd,amne,amnd,amed,az,rp,alph2,beta2,rho2,f(),df(),vslow())
      else
        strike= sorpar(1,isor)
        dip   = sorpar(2,isor)
        rake  = sorpar(3,isor)
        if tam(isor)<0 then 'correct for negative moment
          tam(isor)=-tam(isor)
          if rake <180 then rake=rake+180 else rake=rake-180
        end if
        call m5radpat (mkstat,kpwave,strike,dip,rake,az,rp,alph2,beta2,rho2,f(),df(),vslow())
      end if

'--- Adjust amplitudes and derivatives to moment and reference geom. spreading
      xamp=tam(isor)*ucor*gsp0
      
      for i=1 to 3
        f(i)=f(i)*xamp
        if mkstat then for j=1 to nparmech : df(i,j)=df(i,j)*xamp : next
      next

'--- Get delay times for pP and sP, or sS
'--  PDLY0(1) = delay for pP,sS,  PDLY0(2) = delay for sP
      call deltime (kpwave,rpp,alph2,beta2,dephs,pdly0())

'-- event delay and position delay relative to event 1
'       time offset         position                   depth difference
      evdly = delort - rp * delpos * cos(delphi-az) + (depth1-dephs)*vslow(1)

      evdly = evdly/dt                  'change to units of DT
      pdly(1)=evdly                     'direct phase
      pdly(2)=pdly0(1)/dt+evdly         'pP or sS
      pdly(3)=pdly0(2)/dt+evdly         'sP

      if pdly(1)<0 or pdly(2)<0 or pdly(3)<0 then
        call scprint (22,10,"*** Phase arrival time is negative ***")
        itf=0
        goto 2520
      end if

'-- loop for each phase:  calculate full seismograms and derivatives

      for ip=1 to nph+1
        kflag= (ip=1)  ' flag direct phase
        tau=tau(ip,ista)
        kinteg=(tau/dt=int(tau/dt)) 'if TAU is integer multiple of DT

' --- get elementary seismogram
        redim es0(2*nes)
        ipm1=ip-1
        if not kflag then call m5esio (1, ista, ipm1, nes, es0() )

'------ generate synthetic: --------

'--- For each STF element set ITT, delay for start of elem. seis., and ATAU(i)
'---  interpolate for each delay time: source+phase+STF delays

        for idt=1 to ntau(isor)

          gosub CHECKIN

          atau=atau(idt,isor)
          delay= pdly(ip)*dt + (idt-1)*tau    '--- delay in seconds for start of ES
          itt=int(delay/dt)                   '--- point before start of ES(tau,ph)
          itdly=itt-int(pdly(ip))             '--- start point after start of phase
          dd=delay/dt-itt : omdd=1-dd

'--- loop through the seismogram and interpolate values
          lim=nesp-itt-1
          if lim<1 then goto 1350
          f=f(ip)
          if not kflag then nn=nphes + ip-2

          if not kinteg or idt=1 then
            if kflag then
              for i=1 to lim : es2(i)=sp(i+isprec1)*dd + sp(i+1+isprec1)*omdd : next i
            else
              for i=1 to lim : es2(i)=es0(i)*dd + es0(i+1)*omdd : next i
            end if
          end if

'-- ES2 is interpolated so should be added to SS such that ES2(1) adds to
'--  SS( int(delay/dt) +1) or ITT+1

          if atau<>0 then
            ii=itt
            iid=itdly
            for i=1 to lim
              ii=ii+1
              iid=iid+1
              ess=es2(i)*atau
              es3(iid,ip)=es3(iid,ip) + ess
              ss(ii) = ss(ii) + ess*f
            next i
          end if
997
'--- calculate STF derivatives
'    SS is linear in these, so just do SS without STF elem
          if stfdamp(isor)>0 and mkstat and iafirst and needtau and not noninv then
            ii=itt
            kp=ktau(idt,isor,ia)
            for i=1 to lim
              ii=ii+1
              d(ii,kp)=d(ii,kp) + f*es2(i)
            next i
          end if

1350    next idt

      next ip

'-- SS, the synthetic seismogram, is now formed

'-- Setup and do an inversion
'-- set up matrix of derivatives D

      if (not mkstat) or (not ievflag(isor,ia)) then goto 1920

'-- Accumulate derivatives:
'    IT = index for point before first point of phase
'    ES3(i,ip) = complete elementary seismogram for phase IP
'    KPAR(kp,isor,ia) = index in jacobian for parameter kp & source isor

'-- derivatives for strike, dip, rake or for Mnn,Mee,Mdd,Mne,Mnd,Med
'     as stored in DF array, *ATAU*SP [or EP], first for direct, then refl.
      if ialast then
      for ip=1 to nph+1
        gosub CHECKIN
        it=int(pdly(ip)) : lim=ninv-it
        if lim> 1 then

          for kvar=1 to nparmech
            if damp(kvar+noffmech,isor)>0 then

              kp=kmech(kvar,isor,ia)
              df=df(ip,kvar)
              for i=1 to lim
                k=it+i
                d(k,kp)=d(k,kp) +df*es3(i,ip)
              next i

            end if   '-- damp
          next kvar
        end if    '-- lim
      next ip

'--depth derivatives for pP and sP
' IT = start of phase, LIM = # pts to do in phase, NEQ = current pt in SS
' Calculated by finite differencing w/ time and multiplying by vertical slowness

      if damp(4,isor) > 0 then
        kp=kpar(1,isor,ia)

        for ip=1 to nph
          jp=ip+1             'the PDLY and F index for reflected phases
          it=int(pdly(jp)) : lim=ninv-it-1
          if lim<1 then goto 1610
          neq=it+1 : f=f(jp)*vslow(ip) : f1=f/(2*dt)
          d(neq,kp)=d(neq,kp) - f*(es3(2,jp)-es3(1,jp))/dt
          if lim<2 then goto 1610
          for i=2 to lim
            neq=it+i
            d(neq,kp)=d(neq,kp) - (es3(i+1,jp)-es3(i-1,jp))*f1   'DEPTH
          next i
1610    next ip
      end if


'----- derivatives for origin time and position of secondary events

      for kvar= 2 to 4
        gosub CHECKIN
        if damp(kvar+3,isor) > 0 then
          kp=kpar(kvar,isor,ia)

          for ip=1 to nph+1
            it=int(pdly(ip)) : lim=ninv-it-1
            if lim<1 then goto 1620
'----- du/dDELORT = -du/dt
            if kvar=2 then f = -f(ip)
'----- du/dDELPOS = du/dt rp cos(DELPHI- Az)
            if kvar=3 then f =  f(ip)*rp*cos(delphi-az)
'----- du/dDELPHI = -du/dt rp DELPOS sin(DELPHI-Az)
            if kvar=4 then f = -f(ip)*rp*delpos*sin(delphi-az)

            neq=it+1 : f1=f/(2*dt)
            d(neq,kp)=d(neq,kp) + f*(es3(2,ip)-es3(1,ip))/dt
            if lim<2 then goto 1620
            for i=2 to lim
              neq=neq+1
              d(neq,kp)=d(neq,kp) + (es3(i+1,ip)-es3(i-1,ip))*f1
            next i
1620      next ip
        end if
      next kvar
      end if

'=====  Derivatives are now formed in matrix D ======

1920  next isor
'=== end of loop for sources ==========

'====  Do Statistics :  ========

    if (not noninv) and mkstat then
'-- do statistics for this seismogram
      for i=1 to ninv : tmp(i)=ss(i)*xmag : next
      call makestat(tmp(),odat(),wt,ninv,invflag,r(),stavar(ista),synorm)

'--- weight derivatives with station weights, and magnification
'     and normalize by synthetic rms, if just doing shape inversion

      x= xmag*wt/synorm
      for  j=1 to npar(ia)
        for i=1 to ninv : d(i,j)=d(i,j)*x : next
      next

'---  Add this station to matrix B=DtD, and vector Z=DtR
      call buildb (npar(ia), ninv, d(),r(),z(),b())
    end if    '--- doing statistics

'--- calculate rms and write synthetics to file or array on last iteration
    if itlast and ialast then
      ipof=ipf2(ista)-1
      if ninv(ista)=0 then n=nes(ista) else n=ninv(ista)
      if ipof+n>ndat(ista) then n=ndat(ista)-ipof
      call m5synio (1,0,ista,n,ipof,odat(),dext(ista),datrms(ista))
      s2=0 : for i=1 to n : s2=s2+ss(i)*ss(i) : next
      synrms(ista)= sqr(s2)*xmag
      call m5synio (0,1,ista,nes(ista),0,ss(),amax,dummy)
    end if

2090 next ista   '------ End of station loop -------

'-- set up statistics
  if mkstat then
    varnowt=res2/neq2   'unweighted variance
    if iafirst then 
      totsum0=0 : vardat0=0 : varres0=0
    end if

'-- get variances of data and residuals for each phase-data type
    for i=1 to nphtype
      if sumwtph(i)>0 then
        totsum0=totsum0+sumwtph(i)
        vardat0=vardat0+datvar(i)
        varres0=varres0+resvar(i)
        datvar(i)= datvar(i)/sumwtph(i)
        resvar(i)= resvar(i)/sumwtph(i)
      elseif iafirst then
        datvar(i)=0
        resvar(i)=0
      end if
    next

'-- VARDAT and VARRES are the variances of the data and residuals
    if totsum0>0 then
      vardat=vardat0/totsum0
      varres=varres0/totsum0
    end if
    dvar(itno)=varres

'-- DATERR is estimated unweighted standard error in residuals
    daterr=sqr(res2/(neq2-npar(ia)))

    nrow=neq2 : ncol=npar(ia)  '-- total # data, params

'--- Damping: set least squares dampers scaled to data variances on
'      first iteration.  THETA = damping factor, DIAGWT has parm variance
    if (theta>0) or ipcflag then gosub ADDLSDAMP

'----- solve for adjustments:  Gauss-Jordan inversion for Binv and (Binv)Z
'      ISING flags a singular matrix


    call m5gjinv(npar(ia),b(),z(),d,ising)
    if ising then call scprint (22,10,"*** Singular Matrix *** ") : itf=0 : goto 2520

'---- check STF positivity, and damp out negative STF elements
    if ipcflag and needtau and iafirst then gosub POSITIVT

'----- calculate standard errors of parameters
'-- for orientation or MT elements
    for isor=1 to nsor
    if ialast then
      for i= 1 to nparmech
        kp=kmech(i,isor,ia)
        if damp(i+noffmech,isor)>0 then
          ermech(i,isor)=daterr * sqr(abs(b(kp,kp)))
          if not mtflag then ermech(i,isor)=fndeg(ermech(i,isor))
        else
          ermech(i,isor)=0.
        end if
      next

'-- for depth and multiple source zones
      for i= 1 to 4
        kp=kpar(i,isor,ia)
        if damp(i+3,isor)>0 then
          erpar(i,isor)=daterr * sqr(abs(b(kp,kp)))
          if i=4 then erpar(i,isor)=fndeg(erpar(i,isor))
        else
          erpar(i,isor)=0.
        end if
      next
    end if
'-- for source time function
    if iafirst then
      x=0 : x2=0
      for i=1 to ntau(isor)
        if stfdamp(isor)>0 and needtau then
          kp=ktau(i,isor,ia)
          ertau(i,isor) = daterr * sqr(abs(b(kp,kp)))
          for j=1 to ntau(isor)
            b=b(kp,ktau(j,isor,ia))
            at2 = atau(i,isor)+z(kp)*stfdamp(isor)
            x=x+b
            if at2>0 then x2=x2+i*j*b
          next j
        else
          ertau(i,isor)=0.
        end if
      next
      x=fnamax(x,0.0)
'-- standard error in duration
      ert95(isor)=4.*tau0*daterr*sqr(abs(x2))
      ermom(isor)=sqr(x)*tam(isor)*daterr
    end if
  next isor

'-----  write events, one to a screen
  if ialast then
    for isor=1 to nsor
      iscreen=isor+1
      itline= ((itno-1) mod 8) +1
      call scprint1 (2+itline,1,blank$)
      call scvint   (2+itline,1,2,itno)
      call scprint1 (3+itline,1,blank$)
      call scprint1 (3+itline,1,"Err")
      if mtflag then lvar=5 : ldec=2 else lvar=6 : ldec=1
      il=4
      for i=1 to nparmech
        x=sorpar(i+noffmech,isor)
        if not mtflag then x=fndeg(x)
        x=fnfix(x,2)
        call scvreal( 2+itline, il, lvar,ldec, x)
        if damp(i+noffmech,isor)>0 then
          x=fnfix(ermech(i,isor),2)
          call scvreal( 3+itline, il, lvar,ldec, x)
        end if
        il=il+lvar
      next

      ldec=1
      for i=1 to 4
        x=sorpar(i+3,isor)
        if i=1 then
          x=x+cthk
        elseif i=2 and mtflag then
          lvar=4
        elseif i=4 then
          if mtflag then ldec=0 : lvar=5
          x=fndeg(x)
        end if
        x=fnfix(x,2)
        call scvreal( 2+itline, il, lvar,ldec, x)
        if damp(i+3,isor)>0 then
          x=fnfix(erpar(i,isor),2)
          call scvreal( 3+itline, il, lvar,ldec, x)
        end if
        il=il+lvar
      next

      if mtflag then il=51 : lvar=8 else il=47 : lvar=9
      call scvexp( 2+itline, il, 5,2, tam(isor))
      call scvexp (3+itline, il, 5,2, ermom(isor)) : il=il+lvar+1
      x1=varnowt*varcor : x2=varres*varcor
      call scvexp (2+itline, il, 5,3, x1) :  il=il+lvar
      call scvexp (2+itline, il, 5,3, x2) :  il=il+lvar
      if vardat=0 then x=0 else x=varres/vardat*100
      if x<999 then call scvreal (2+itline, il, 5,1, x) else call scprint1(2+itline,il," >999")

'-- time function elements
      ialine=12+itline : ieline=13+itline : il=4
      if ntau(isor)>18 then ialine=13 : ieline=16

      call scprint1 (ieline-1, 1, "   ")
      call scvint (ialine,1,2,itno)
      for i=1 to ntau(isor)
        ix=atau(i,isor)*100
        call scvint(ialine, il, 4, ix)
        if stfdamp(isor)>0 then
          if ertau(i,isor)<=99.995 then
            ix=ertau(i,isor)*100
            call scvint(ieline, il, 4, ix)
          else
            call scprint1 (ieline,il," Big")
          end if
        end if
        il=il+4
      next
      if stfdamp(isor)>0 then call scprint1 (ieline, 1, "Err")

      next isor
      call scprint2 (npage)
    end if
  end if

'-- check on percent variance reduction
  if itno>1 and ialast then
    varred = (dvar(itno-1) - dvar(itno))/dvar(itno-1) * 100

    if varred < varmin then
      ivarstop=ktrue
      call scprint (23,21,"Variance reduction below cutoff")
      if (varred<-5) then goto 2520 'trap for runaway variance
    else
      if itno>2 and itno<nit then 'allow continuing iteration
        if dvar(itno)<dvar(itno-2) then itlast=kfalse
      end if
      ivarstop=kfalse
      call scprint (23,21,"Variance reduction             ")
      call scprint(23,40,str$(fnfix(varred,1))+"%")
    end if

  end if

'-- if not last iteration
  if itlast and ialast then exit for

'---- adjust model: scale to dampers -------------
  for isor=1 to nsor
    if ialast then
      for i=1 to nparmech
        ik=i+noffmech
        if damp(ik,isor)>0 then
          kp=kmech(i,isor,ia)
          sorpar(ik,isor)=sorpar(ik,isor)+z(kp)*damp(ik,isor)
        end if
      next

      for i=1 to 4
        ik=i+3
        if damp(ik,isor)>0 then
          kp=kpar(i,isor,ia)
          sorpar(ik,isor)=sorpar(ik,isor)+z(kp)*damp(ik,isor)
        end if
      next

'--- reset depth if above halfspace
      sorpar(4,isor)=fnamax(0.1, sorpar(4,isor))

'--- adjust strike, dip, rake to keep in bounds
      if not mtflag then
        strike=fndeg(sorpar(1,isor))
        dip=fndeg(sorpar(2,isor)) mod 360
        rk=fndeg(sorpar(3,isor))

        if dip<0 then 'correct for negative dip
          if dip >= -180 then
            dip= -dip
            strike=strike+180.
            rk=rk+180.
          else
            dip=360+dip
          end if
        end if

        if dip>=90 then  'max dip is 90 degrees
          dip=180-dip
          strike=strike + 180.
          rk = 360-rk
        end if

        if not fninlim(strike,-360,360) then strike = strike mod 360
        if strike<0 then strike=strike+360.
        if not fninlim(rk,-360,360) then rk = rk mod 360
        if rk<0 then rk=rk+360.
        sorpar(1,isor)=fnrad(strike)
        sorpar(2,isor)=fnrad(dip)
        sorpar(3,isor)=fnrad(rk)
      end if

      delphi=fndeg(sorpar(7,isor))
      if not fninlim(delphi,-360,360) then delphi = delphi mod 360
      if delphi<0 then delphi=delphi+360.
      sorpar(7,isor)=fnrad(delphi)
    end if 'if ialast

'--- STF adjustments
    if iafirst then
      if needtau and stfdamp(isor)>0 then
        for i=1 to ntau(isor)
          kp=ktau(i,isor,ia)
          atau(i,isor)=atau(i,isor)+z(kp)*stfdamp(isor)
        next
      end if
'---- normalize source time function and adjust moment
      d=tam(isor)
      call m5stf (invflag,tau0,isor,ntau(),atau(),tam(),tcent,t95)
      if mtflag and needmt then tam(isor)=d
    end if 'iafirst

    if mtflag and needmt and ialast then
'---- normalize moment tensor elements, source time function, and adjust moment
      call mtadj(sorpar(),isor,tam(),7)
    end if

'-- adjust moment for shape inversion by the average of the RMS ratios
    if not invflag then
      if sumdswt>0.0 then tam(isor)=tam(isor)*sumdos/sumdswt
    end if

  next isor

 next ia 
'---- End INNER LOOP

'-- check for input
 gosub CHECKIN
 totime = int(timer-time1)

next itno
itf=itf+1 : lastit=mtflag 'flags another successful iteration
'---- End ITERATION LOOP

if mtflag then
  for isor=1 to nsor 'adj moment and double couple for mt solution
    call mteigen(kfalse,sorpar(),isor,tam(),tamiso(),ierr)
  next
end if

' -- check for variable range parameter looping
if nvarp>0 then
  gosub SUMLINE
  for isordp=1 to nsor
    for i=1 to len(vrp$(isordp))
      j=instr("HOXY",mid$(vrp$(isordp),i,1))
      if j<>4 or sorpar(6,isordp)>.1 then
        js=3+j
        if varinc(j,isordp)<>0 then
    'if within 1% of endpoint, go on to new parameter
          if abs((sorpar(js,isordp)-varend(j,isordp))/varinc(j,isordp))<=.01 then
            sorpar(js,isordp)=varbeg(j,isordp)
          else
            v=sorpar(js,isordp)+varinc(j,isordp)
            if varinc(j,isordp)>0 then
              sorpar(js,isordp)= fnamin(varend(j,isordp),v)
            else
              sorpar(js,isordp)= fnamax(varend(j,isordp),v)
            end if
            goto 519
          end if
        end if
      end if
    next i
  next isordp
end if

2520 call cline (23,55)
call scprintc(23,55," Press any key...   ")
call inlet(kfalse, d$, klen, n, kascii)
if fninlim(n,1,nsor) then n=n+1 : gosub NEWPAGE : goto 2520
if d$="h" then n=1 : gosub NEWPAGE : goto 2520

itf=itf+1 : lastit=mtflag
if open2 then close 2

2525 chain SYNDIR$+"MT5MEN1"

'++++++++++++++++++++ END PROGRAM ++++++++++++++++++++

'******** SUBROUTINES **************
'========================================================================
POSITIVT:
'--- Subroutine to impose positivity on STF by increasing damping
'      The dampers on STF elements that went below STFMIN are incresed
'      by DAMPFAC, and the inversion is redone.  This process is repeated
'      until the constraint is achieved
 
  kdflag=kfalse : dampmax=1e37

'--- TAUDAMP is total damping on STF elements, initialized to theta or 1
  for isor=1 to nsor: for i=1 to ntau(isor) : taudamp(i,isor)=1 : next i,isor

6600
'--- start by checking if any negative elements exist, and setting them to
'     zero. Also check # of increases in damping
   kneg=kfalse
    for isor=1 to nsor
       for i=1 to ntau(isor)
         kp=ktau(i,isor,ia)
          if kp>0 then atau2(i,isor)=atau(i,isor)+z(kp)*stfdamp(isor) else atau2(i,isor)=atau(i,isor)
          if atau2(i,isor)<stfmin then kneg=ktrue : atau(i,isor)=0.
       next i, isor

'--- if STF has no negative elements, or too much damping, quit subroutine
     if not kneg then return

     if itry>ntrys then
       call scprint (22,10,"STF negative, max damping reached")
       return
     end if

'    call cline(22,10)
'    call scprint (22,10,"STF elements negative: ")

'-- get original matrices and add least square dampers to parameters

       for i=1 to ncol
         z(i)=z2(i)
         for j=1 to i
           b(i,j)=b2(i,j)
         next j,i

 for i=1 to npar(ia) : b(i,i)=b(i,i)+theta*diagwt : next

'--- check STF and add dampers where needed
'      TAUDAMP is increased for negative STF elements. Limit damping to prevent 
'      overflow errors

  for isor=1 to nsor
   if stfdamp(isor)>0 then
     for i=1 to ntau(isor)
      kp=ktau(i,isor,ia)
       if (atau2(i,isor) < stfmin) then 
'        print "(";i;",";isor;")";
        if (dampmax/taudamp(i,isor)) =< dampfac then kdflag=ktrue _
           else taudamp(i,isor)=taudamp(i,isor)*dampfac 
       end if 
      if taudamp(i,isor)>1 then b(kp,kp)=b(kp,kp)+taudamp(i,isor)
     next i
  end if
 next isor

'--- re-invert matrix
   call m5gjinv(npar(ia),b(),z(),d,ising)

'--- quit if too much damping
  if kdflag then call scprint (22,10,"STF negative, max damping reached"): goto 6900

'--- restart subroutine to check if STF is now in acceptable bounds
     goto 6600
     
6900 return
'========================================================================
ADDLSDAMP: '---- Subroutine to add least squares damping
'  Scales damping to estimates of residual variance from first iteration
'  Input:  THETA, damping factor
'          B, the matrix being damped
'
'--- store the B (AtA) matrix and the Z (AtR) data array aside
     for i=1 to ncol
       z2(i)=z(i)
       for j=1 to i
         b2(i,j)=b(i,j)
       next
     next

'-- on first iteration, store residual variance
     if itno=1 then diagwt=varres

'--- add damping to diagonal of B
     if theta>0 then  
         for i=1 to ncol : b(i,i)=b(i,i)+theta*diagwt : next
      end if

return
'========================================================================
PRTIME:
  locate 23,1: print left$("Time = "+fnstrnum$(int(timer-time1))+"/"+fnstrnum$(totime)+"s        ",20)
  call scprint (23,56,"Station "+fnstanm$(ista)+" Phase "+fnkphlab$(kphtype(ista)))
 return
'========================================================================
NEWPAGE:
      if n<>npage then call scprint2(n) : npage=n
 return
'========================================================================
sub buildb (npar, ninv, d(2),r(1),z(1),b(2) ) static

'---- Subroutine to build B=DtD and Z=DtR.
'          Only lower-triangle of B stored.
 for k=1 to npar
  for l=1 to k
   for i=1 to ninv
     b(k,l)=b(k,l) + d(i,l)*d(i,k)
   next
  next
   for i=1 to ninv
    z(k)=z(k) + r(i)*d(i,k)
   next
 next
end sub
'========================================================================
sub deltime (kpwave,rpp,alph,beta,dephs,pdly0(1)) static
'---------- Subroutine to get upgoing P and S delay times in the halfspace
'               ALPH, BETA - source structure
'               DEPHS - source depth in halfspace
'        Output:  PDLY0(1) = delay for pP or sS
'                 PDLY0(2) = delay for sP

    if kpwave then vel=alph else vel=beta
    z=fnasin(rpp*alph)  'P or SH angle
    y=fnasin(rpp*beta)  'SV angle

    putime=dephs/(vel*cos(z))        'upgoing P or SH time
    sutime=dephs/(beta*cos(y))       'upgoing SV time
    rx=dephs*tan(z) : svx=dephs*tan(y)    'horizontal distances traveled
   
    pdtime=(putime-2*sin(z)*rx/vel)      'P or SH down corrected to wavefront
    sdtime=(putime-(rx+svx)*sin(z)/vel)  'P down in sP phase

    pdly0(1)=(putime+pdtime)
    pdly0(2)=(sutime+sdtime)
end sub
'========================================================================
sub geomspr(delta,rpp,rho2,rhor,alph2,alphr,gspr) static

'--------- Subroutine for geometrical spreading and receiver crust response
'     Output: GSPR, just the geometrical spreading factor
'
     shared dgs(), cgs()

 if fninlim(delta,30,90) then
'--- get gspr from delta
     for i=2 to 16  'get the distance index
       if dgs(i)>delta then 
        delx=delta-dgs(i-1)
        dg=cgs(i)-cgs(i-1): dd=dgs(i)-dgs(i-1)
        gspr=cgs(i-1)+delx*dg/dd 
        exit for
       end if
     next
 elseif delta < 30. then
     delx=delta-30: dg=cgs(2)-cgs(1)
     gspr=cgs(1)+delx*dg/2.5
 elseif delta > 90. then
     delx=delta-90: dg=cgs(16)-cgs(15)
     gspr=cgs(16)+delx*dg/5.0
 end if

'--- adjust GSPR for appropriate source and receiver velocities
     rpt=rpp: rpt2=rpt^2
     a2=1/(6.5)^2: as2=1/alph2^2: ar2=1/alphr^2
     a=(rho2/rhor)*alph2^2/alphr^2
     c=a*(a2-rpt2)/(sqr(as2-rpt2)*sqr(ar2-rpt2))
     gspr=sqr(c)*gspr
end sub
'========================================================================
sub recvresp (kpwave, rpp, alphr, betar, acor) static
'---Receiver half-space response:  include free-surface effect at receiver
' RPP= ray parameter       ALPHR, BETAR = Vp, Vs of receiver structure
' ACOR = amplitude correction

     acor=2: if not kpwave then exit sub
     aipr=fnasin(rpp*alphr)
     aisr=fnasin(rpp*betar)
     dd=cos(2*aisr)^2+betar^2/alphr^2*sin(2*aipr)*sin(2*aisr)
     acor=2*cos(aipr)*cos(2*aisr)/dd
end sub
'========================================================================
sub makestat(ss(1),odat(1),wt,ninv,invflag,r(1),stavar,synorm) static
'---- routine to calculate residuals and update statistics for an individual 
'     seismogram

  shared wres2, res2, sumwt, ista, gcor()
  shared datvar(), resvar(), sumwtph(), kphtype()
  shared sumdos, sumdswt

'-- sums of amplitudes squared
      sums=0 : sumd=0
      for i=1 to ninv 
        sums=sums + ss(i)*ss(i)
        sumd=sumd + odat(i)*odat(i)
      next

'-- if doing shape inversion calculate normalizations SYNORM and DATNORM 
    synorm=1.0 : datnorm=1.0
    if not invflag then 
      synorm=sqr(sums) : datnorm=sqr(sumd)
      if synorm<>0 then sumdos=sumdos+wt*datnorm/synorm : sumdswt=sumdswt+wt
    end if

'-- Form residuals with normalizations, weighted residuals WRES,
'    and various sums of squares and variances

     wsrms=0 : srms=0 : sumd=0
      for i=1 to ninv
         odat=odat(i)/datnorm : sumd=sumd+odat*odat*wt*wt
         res=odat-ss(i)/synorm
         wres=res*wt : r(i) = wres
         res2=res2+res*res
         wres2=wres2+wres*wres
         srms=srms+res*res         'unweighted station res
         wsrms=wsrms+wres*wres     'weighted station res
      next

      sumwt=sumwt+ninv*wt  

'-- weighted data and residual variancess by PDT
      sumwtph(kphtype(ista))=sumwtph(kphtype(ista))+ninv*wt 
      resvar(kphtype(ista))=resvar(kphtype(ista))+wsrms
      datvar(kphtype(ista))=datvar(kphtype(ista))+sumd
      stavar=srms/ninv
end sub

'========================================================================
CHECKIN:
'-- check for input of new source to be displayed on screen
   d$=lcase$(inkey$) : n=val(d$)
   if fniquit(d$) then goto 2520
   if fninlim(n,1,nsor) then n=n+1 : gosub NEWPAGE
   if d$="h" then n=1 : gosub NEWPAGE
   gosub PRTIME
return

'========================================================================
'-- screen subroutines

sub scprint1(ix,iy,text$) static
 shared iscreen, liner$(), ncharm
  ichar = (ix-1)*80 + iy
  if ichar<= ncharm then mid$(liner$(iscreen),ichar)=text$
end sub
'-----------------------------------------------------------------------
sub scprint2 (k) static
 shared iscreen, liner$(), ncharm
 for i=1 to 21
  locate i,1 : print mid$(liner$(k),80*(i-1)+1,80);
 next
  if k>1 then locate 1,7 : print k-1;
end sub
'-----------------------------------------------------------------------
sub scvint( ix, iy, klen, kvar) static
 shared iscreen, liner$(), ncharm
  ichar = (ix-1)*80 + iy
  text$ = string$(klen,32)+str$(kvar)
  text$=right$(text$,klen)
  if ichar<= ncharm then mid$(liner$(iscreen),ichar)=text$
end sub
'-----------------------------------------------------------------------
sub scvreal( ix, iy, klen, kfrac, varin) static
 shared iscreen, liner$(), ncharm
  ichar = (ix-1)*80 + iy
  var=varin

  v$=str$(var) : kE=instr(v$,"E")
    if kE>0 then
     iex=val(mid$(v$,kE+1,3))
     if -iex>kfrac then var=0 else var=val(mid$(v$,1,kE-1)) * 10^iex
    end if

  text$ = "              "+str$(var) : kdot=instr(text$,".")
  if kdot=0 then text$=text$+"."     : kdot=instr(text$,".")
  text$=text$+"000000000"
  k1=kdot-(klen-kfrac-1)
  text$=mid$(text$,k1,klen)

  if ichar<=ncharm then mid$(liner$(iscreen),ichar)=text$
end sub
'-----------------------------------------------------------------------
sub scvexp ( ix, iy, klen, kfrac, var) static
 shared iscreen, liner$(), ncharm

  kend=klen+iy-1
  ichar = (ix-1)*80 + iy

  if var<>0. then
  v=fnalog10(abs(var))
  iexp=int(v)
  else
   iexp = 0
  end if

  v=var/(10^iexp)
  exp$=str$(iexp)
  if iexp>=0 then exp$=right$(exp$,len(exp$)-1)
  exp$="E"+exp$

  text$ = "            "+str$(v)
  kdot=instr(text$,".")
  if kdot=0 then text$=text$+"." : kdot=instr(text$,".")
  text$=text$+"000000000"
  k1=kdot-(klen-kfrac-1)
  text$=mid$(text$,k1,klen)+exp$

  if ichar<= ncharm then mid$(liner$(iscreen),ichar)=text$
end sub
'========================================================================
SUMLINE:
'-------- print summary lines

open "a",9,"MT5VARY.LOG"
width #9,100
if ndhdrdp then 
  print #9, fhdr$;" (";fnstrnum$(nvarp);" steps of ";fnstrnum$(nit);" iterations/step)    ";date$;"   ";time$
  if mtflag then
    print #9, "Src Mnn  Mee  Mdd  Mne  Mnd  Med Depth   dTo  Dist    Az    Moment  UnwtdVAR    WtdVAR  R/D%"
  else
    print #9, "Src Strike Dip  Rake Depth  dTo   Dist    Az    Moment  UnwtdVAR    WtdVAR  R/D%"
  end if
  xvarmin=1e30
  ivarpmin=0
  ndhdrdp=kfalse
endif

for isorsl=1 to nsor
  print #9,using fniform$(2);isorsl;
  for i=1 to nparmech
    x=sorpar(i+noffmech,isorsl)
    if not mtflag then x=fndeg(x)
    print #9, using fnfform$(6+mtflag,1-mtflag); fnfix(x,2);
  next
    
  for i=4 to 7
    x=sorpar(i,isorsl)
    if i=4 then x=x+cthk
    if i=7 then x=fndeg(x)
    print #9, using fnfform$(6,1); fnfix(x,2);
  next

  print #9, using fndform$(6,3);tam(isorsl);
  if isorsl=nsor then
    print #9,using fndform$(6,3);x1;x2;
    if vardat=0 then x=0 else x=varres/vardat*100
    if x<999 then print #9,using fnfform$(6,1); x; else print #9," >999";
    if x<xvarmin then xvarmin=x:ivarpmin=ivarp
  end if
  print #9,
next isorsl
if ivarp=nvarp then 
  a$="Minimum R/D% of "+fnstrnum$(xvarmin)+" on step "+fnstrnum$(ivarpmin)
  print #9,a$
  call scprint(22,10,a$)
end if
close 9
return
'========================================================================
' rem $INCLUDE: 'm5rapars' 'now linked separately
rem $INCLUDE: 'm5stf'
rem $INCLUDE: 'm5gjinv'
rem $INCLUDE: 'm5radpat'
rem $INCLUDE: 'm5radmt'
' rem $INCLUDE: 'm5eigen'  'now linked separately
'========================================================================
