'-- M5EIGEN.BAS   EIGENVALUE/EIGENVECTOR ROUTINES
'   INCLUDES QR ALGORITHM USED W/HOUSEHOLDER'S METHOD OF SIMILARITY TRANSFORMATION
'   By Peter Zwick    May 4, 1990
'   Converted from C to QUICKBASIC for MT4 update -  9/11/1991
'  (compile and link separately to MT5MEN1, MT5MEN2, MT5IT)
'  if debug mode desired, uncomment following: 
'     "shared ipmi"   (in sub getfp)
'     "call prt...()" (in all subs)
'     "rem $include: 'm5prtsub'" (at end)

DEFINT I-N: OPTION BASE 1

rem $INCLUDE: 'M5COMMON'
rem $INCLUDE: 'M5FUNCTN'

'========================================================================
'-- SUB MMULT - matrix multiplication for n x n matrix
SUB mmult(n, A(2), B(2), C(2)) static
'  returns C() as product of matrices A() and B()
  for i=1 to n
    for j=1 to n
      C(i,j)=0
      for k=1 to n
        C(i,j)=C(i,j)+A(i,k)*B(k,j)
    next k,j
  next
END SUB
'========================================================================
'-- SUB HOUSEHOLDER - returns A as a tridiagonal matrix
SUB householder(n, A(2)) static 'simplified for 3x3 symmetric matrix
'  Inputs n,A() - rank of input matrix and matrix
'  Output A - tridiagonal matrix

  if A(3,1)<>0 then 'A already symmetric and tridiagonal if A(3,1)=0
    redim C(n,n),P(n,n),w(n)
    suma=0
    for k=2 to n : suma=suma+A(k,1)*A(k,1) : next k
    if A(2,1)>0 then alpha#= -sqr(suma) else alpha#= sqr(suma)
    w(1)=0
    w(2)=sqr(.5-A(2,1)/(2*alpha#))
    for j=3 to n : w(j)= -A(j,1)/(2*w(2)*alpha#) : next j
    for j=1 to n : for k=1 to n : P(j,k)= -(j=k) -2*w(j)*w(k) : next k,j
    call mmult(n,P(),A(),C())
    call mmult(n,C(),P(),A()) ' A() is now transformed submatrix
  end if
END SUB
'========================================================================
'-- FUNCTION IQR - QR algorithm to get eigenvalues ----
FUNCTION iqr(irank, a(1),b(1),eigval(1),tol,niter) static
'-- customized for 3x3 symmetric tridiagonal matrix input
  redim d(irank),x(irank),y(irank),z(irank),c(irank)
  redim s(irank),q(irank),r(irank)
  n=irank
  shift#=0.
  ilam=1
  for k=1 to niter
    if abs(b(n))<=tol then
      eigvalue=a(n)+shift#
      eigval(ilam)=eigvalue : ilam=ilam+1
      n=n-1
    elseif abs(b(2))<=tol then
      eigvalue=a(1)+shift#
      eigval(ilam)=eigvalue : ilam=ilam+1
      n=n-1
      a(1)=a(2)
      for j=2 to n
        a(j)=a(j+1)
        b(j)=b(j+1)
      next
    end if
    bb#= -(a(n-1)+a(n))
    cc#=a(n)*a(n-1)-b(n)*b(n)
    dd#=sqr(bb#*bb#-4*cc#)
    if bb#>0 then
      amu1= -2.*cc#/(bb#+dd#)
      amu2= -0.5*(bb#+dd#)
    elseif dd#-bb#<>0 then
      amu1= 0.5*(dd#-bb#)
      amu2= 2.*cc#/(dd#-bb#)
    else
      exit for  'no solution
    end if
    if n=2 then
      eigvalue=amu1+shift#
      eigval(ilam)=eigvalue : ilam=ilam+1
      eigvalue=amu2+shift#
      eigval(ilam)=eigvalue : ilam=ilam+1
      iqr= 0 : exit function
    end if
    if abs(amu1-a(n))<abs(amu2-a(n)) then ss#=amu1 else ss#=amu2
    shift#=shift#+ss#
    for j=1 to n : d(j)=a(j)-ss# : next
    x(1)=d(1)
    y(1)=b(2)
    for j=2 to n
      z(j-1)=sqr(x(j-1)*x(j-1)+b(j)*b(j))
      if z(j-1)<1e-5 then iqr=k : exit function
      c(j)=x(j-1)/z(j-1)
      s(j)=b(j)/z(j-1)
      q(j-1)=c(j)*y(j-1)+s(j)*d(j)
      x(j)= -s(j)*y(j-1)+c(j)*d(j)
      if j<>n then
        r(j-1)=s(j)*b(j+1)
        y(j)  =c(j)*b(j+1)
      end if
    next j
    z(n)=x(n)
    a(1)=s(2)*q(1)+c(2)*z(1)
    b(2)=s(2)*z(2)
    for j=2 to n-1
      a(j)=s(j+1)*q(j)+c(j)*c(j+1)*z(j)
      b(j+1)=s(j+1)*z(j+1)
    next j
    a(n)=c(n)*z(n)
  next k
'  call prtivar(k,"Error on iteration",ktrue)
  iqr= k
END FUNCTION
'========================================================================
'-- SUB EIGEN - gets eigenvalues for lower diagonal matrix
SUB eigen(n,acc,niter,E(2),eigval(1),ierr) static
'-- Inputs :  n,acc,niter = Rank N, Accuracy, Max iterations
'             E(2) = n x n lower diagonal matrix
'-- Outputs:  eigval() = vector of n eigenvalues
'             ierr = 0 if successful, or # of iteration at failure

  redim a(n), b(n)
  tol=acc
  maxiter=niter
'-- make E() the full matrix
  for i=1 to n-1 : for j=i+1 to n : E(i,j)=E(j,i) : next j,i
'  call prtmtrx(n,n,E(),a(),"Initial","a()")
  call householder(n,E())
  a(1)=E(1,1)
  for i=2 to n
    a(i)=E(i,i)
    b(i)=E(i,i-1)
  next
  ierr = iqr(n,a(),b(),eigval(),tol,maxiter) 'get eigenvalues
'  if ierr then call prtivar(ierr,"Failed in iqr on step",ktrue)
END SUB
'========================================================================
'-- SUB MTISO - gets isotropic and deviatoric parts of moment tensor
SUB mtiso(sorpar(2),isor,amdev,amiso) static
  amdev=0 : amiso=0
  for i=1 to 3
    amiso=amiso+sorpar(i+noffmt,isor)
  next i
  amiso=amiso/3
  for i=1 to 3
    am1=sorpar(i+noffmt,isor)-amiso
    amdev=amdev+0.5*am1*am1
  next i
  for i=4 to 6
    am1=sorpar(i+noffmt,isor)
    amdev=amdev+am1*am1
  next i
  amdev=sqr(amdev)
END SUB
'========================================================================
'-- SUB MTADJ - normalizes mt elements and adjusts seismic moment
SUB mtadj(sorpar(2),isor,tam(1),ihld) static
  if ihld<7 then  'put this single element within bounds
    if ihld<4 then h=sqr2 else h=1.
    sorpar(ihld+noffmt,isor)=fnxlimit(sorpar(ihld+noffmt,isor),-h,h)
  end if
4111  am0=0
  amhld=0
  for i=1 to 6
    hld=sorpar(i+noffmt,isor)*sorpar(i+noffmt,isor)
    if i<4 then hld=0.5*hld
    if i<> ihld then am0=am0+hld else amhld=hld
  next i
  am1=sqr(am0)
  am0=am0+amhld
'  call prtdvar(am0,"am0 in mtadj (ihld="+fnstrnum$(ihld)+")",ktrue)
  if am0>1e-6 then 'if some non-zero mt elements
    if ihld<7 then
      am1=am1/sqr(fnamax(1e-6,1.-amhld))
    end if
    am1=fnamax(1e-6,am1) 'if some non-zero and changeable mt elements 
    for i=1 to 6
      if i<> ihld then sorpar(i+noffmt,isor)=sorpar(i+noffmt,isor)/am1
    next i
    tam(isor)=tam(isor)*sqr(am0)
  else
    if ihld<>4 then 
      sorpar(4+noffmt,isor)=1. 'vert strike-slip fault
    else
      sorpar(5+noffmt,isor)=-1. 'vert dip-slip fault
    end if 
    goto 4111
  end if
END SUB
'========================================================================
'-- SUB MTEIGEN - gets eigenvalues (& optionally eigenvectors) of MT elements
'  returns str,dp,rk of fault planes and adjusted total & isentropic moment:tam(),tamiso()
SUB mteigen(isdblcpl,sorpar(2),isor,tam(1),tamiso(1),ierr) static
shared evct()
redim ev(3,3)

  redim am(3,3),eigv(3),str(2),dp(2),rk(2)
  call mtadj(sorpar(),isor,tam(),7) 'normalize mt matrix and adjust seismic moment
'-- fill am() matrix of deviatoric component of moment tensors
  if isdblcpl then  'remove isotropic component
    call mtiso(sorpar(),isor,amdev,amiso) 'get isotropic part of mt
    tamiso(isor)=tam(isor)*amiso  'and save it
    for i=1 to 3 : am(i,i)=sorpar(i+noffmt,isor)-amiso : next i
  else
    for i=1 to 3 : am(i,i)=sorpar(i+noffmt,isor) : next i
  end if
  am(2,1)=sorpar(4+noffmt,isor): am(1,2)=am(2,1)
  am(3,1)=sorpar(5+noffmt,isor): am(1,3)=am(3,1)
  am(3,2)=sorpar(6+noffmt,isor): am(2,3)=am(3,2)
  for i=1 to 3 : for j=1 to 3 : ev(i,j)=am(i,j) : next j,i
  call eigen(3,.0001,50,ev(),eigv(),ierr)
  if ierr=0 then
    call asort(3,eigv(),kfalse)  'put in order T,B,P axes
    SWAP eigv(2),eigv(3)         'chg to order T,P,B
'    call prtmtrx(3,3,am(),eigv(),"deviatoric am()","eigv() sorted as T,P,B axes")
    am0=0.5*(abs(eigv(1))+abs(eigv(2)))
'    call prtdvar(am0,"am0 in mteigen",ktrue) 'for comparison w/am0 from mtadj
    if isdblcpl then tam(isor)=tam(isor)*amdev   'correct to deviatoric portion only
    call evect(am(),eigv(),ev(),ierr)
'    call prtmtrx(3,3,ev(),eigv(),"adjusted eigenvectors for T,P,B axes","eigenvalues")
    if ierr=0 then 
      for i=1 to 3 : for j=1 to 3 : evct(i,j,isor)=ev(i,j) : next j,i
      call getfp(ev(),str(),dp(),rk())  'get fault plane orientation from eigenvectors
' organize into same orientation as previous strike, dip, rake
      diff=1e10
      for i=1 to 2
        for j= 0 to 1
          t=sorpar(1,isor)+j*180.
          if t>360 then t=t-360.
          if abs(t-str(i))<diff then diff=abs(t-str(i)) : k=i
        next j
      next i
      sorpar(1,isor)=str(k)
      sorpar(2,isor)=dp(k)
      sorpar(3,isor)=rk(k)
    end if
  end if
  if ierr<>0 then sorpar(1,isor)=-999
END SUB
'========================================================================
FUNCTION vlen(j,v(2)) static
  vl=0
  for i=1 to 3 : vl=vl+v(j,i)*v(j,i) : next i
  vlen=sqr(vl)
END FUNCTION
'========================================================================
SUB getfp(ev(2),str(1),dip(1),rk(1)) static 'returns fault plane orientation for eivenvector matrix
  redim vrk(2,3),vstr(2,3),vdip(2,3),vt(3),vp(3),vb(3),vfp(2,3),vtmp(2,3)
'  shared ipmi
  for i=1 to 3 : vt(i)=ev(1,i) : vp(i)=ev(2,i) : vb(i)=ev(3,i) : next i
  for i=1 to 3   'get rake and auxiliary rake vectors
    vrk(1,i)=(vt(i)+vp(i))/sqr2 : vfp(2,i)=vrk(1,i)
    vrk(2,i)=(vt(i)-vp(i))/sqr2 : vfp(1,i)=vrk(2,i)
  next i
'-- adjust fault plane normals and rake directions
  for j=1 to 2   'for both fault plane and auxiliary fault plane
    if vfp(j,3)>0 or (vfp(j,3)=0 and (vfp(j,1)>0 or (vfp(j,1)=0 and vfp(j,2)<0))) then for i=1 to 3 : vrk(j,i)= -vrk(j,i)  : vfp(j,i)= -vfp(j,i): next i
    vstr(j,1)= vfp(j,2) : vstr(j,2)= -vfp(j,1) : vstr(j,3)=0  'get strike vector
    vstrl=vlen(j,vstr()) : if vstrl<1e-5 then vstr(j,1)=1 : vstrl=1
    for i=1 to 3 'get strike unit vector
      vstr(j,i)=vstr(j,i)/vstrl
    next i
    call vxprod(j,vstr(),vfp(),vdip()) 'get dip vector
    str(j)=fnatan2(vstr(j,2),vstr(j,1)) : if str(j)<0 then str(j)=str(j)+2*pii
    dip(j)=abs(fnasin(vdip(j,3)))
    call vxprod(j,vstr(),vrk(),vtmp()) 'get |vstr|*|vrk|*sin(rk)=sin(rk)
    rk(j)=fnasin(vlen(j,vtmp()))
    vm=0 : im=1
    for i=1 to 3  'get max length sin vector
      if abs(vtmp(j,i))>vm then vm=abs(vtmp(j,i)):im=i
    next i
    if sgn(vtmp(j,im))<>sgn(vfp(j,im)) then rk(j)= -rk(j)
    call vdotprod(j,vstr(),vrk(),tmp) 'get |vstr|*|vrk|*cos(rk)=cos(rk)
    if tmp<0 then rk(j)=pii-rk(j)
    if rk(j)<0 then rk(j)=rk(j)+2*pii
'    if ipmi then
'      for i=1 to 3
'        i$="("+fnstrnum$(j)+","+fnstrnum$(i)+")"
'        call prtdvar(vrk(j,i),"vrk"+i$,kfalse)
'        call prtdvar(vstr(j,i),"vstr"+i$,kfalse)
'        call prtdvar(vdip(j,i),"vdip"+i$,ktrue)
'      next i
'    end if
  next j
END SUB
'========================================================================
SUB vxprod(j,v1(2),v2(2),v3(2)) static 'returns v3 as cross product of v1,v2
  redim vt(3)
  vt(1)= v1(j,2)*v2(j,3)-v1(j,3)*v2(j,2)
  vt(2)= -v1(j,1)*v2(j,3)+v1(j,3)*v2(j,1)
  vt(3)= v1(j,1)*v2(j,2)-v1(j,2)*v2(j,1)
  for i=1 to 3 : v3(j,i)=vt(i) : next
END SUB
'========================================================================
SUB vdotprod(j,v1(2),v2(2),dotprod) static 'returns dotprod as dot product of v1,v2
  dotprod=0
  for i=1 to 3 : dotprod=dotprod+v1(j,i)*v2(j,i) : next i
END SUB
'========================================================================
SUB asort(n,a(1),isabs) static 'rearranges array in max to min order (slow for large n)
  for i=1 to n-1 : for j=i+1 to n
    if isabs then
      if abs(a(j))>abs(a(i)) then SWAP a(i),a(j)
    else
      if a(j)>a(i) then SWAP a(i),a(j)
    end if
  next j,i
END SUB
'========================================================================
SUB evect(am(2),eigval(1),ev(2),ierr)
  redim a1(3,3),a2r(3,3),a2(3),nrm(3),ncm(3),asum(3)
  for igv=1 to 3 'for each eigenvalue
    ia=0
    a2max=0 : nrm(1)=0 : ncm(1)=0
    for i=1 to 3
      asum(i)=0
      for j=1 to 3
        a1(i,j)=am(i,j)+(i=j)*eigval(igv) 'create matrix of [A-lambda*I]
        asum(i)=asum(i)+abs(a1(i,j))
      next j
      for j=1 to 3    'build matrix of absolute relative values
        if asum(i)>0 then
          a2r(i,j)=abs(a1(i,j))/asum(i)
          if a2r(i,j)>a2max then   'find row & column of largest element
            a2max=a2r(i,j) : nrm(1)=i : ncm(1)=j
          end if
        else
          a2r(i,j)=0
        end if
      next j
    next i
    ierr=(nrm(1)=0) : if ierr then exit sub 'a row with all zeroes
    a2max=0 : nrm(2)=0 : ncm(2)=0
    for i=1 to 3  'find row & column of next largest element
      if i<>nrm(1) then
        for j=1 to 3
          if j<>ncm(1) and a2r(i,j)>a2max then
            a2max=a2r(i,j) : nrm(2)=i : ncm(2)=j 'find largest element
          end if
        next j
      end if
    next i
    ncm(3)=6-ncm(1)-ncm(2)   'remaining column (move to RHS)
    nrm(3)=6-nrm(1)-nrm(2) 'row that won't be used except to check
    ev(igv,ncm(3))=1
    piv=a1(nrm(2),ncm(1))/a1(nrm(1),ncm(1)) 'get row reduction factors
    for j=1 to 3
      a2(ncm(j))=a1(nrm(2),ncm(j))-piv*a1(nrm(1),ncm(j))
    next j
    k=0
    for i=1 to 2      'do row reduction
      ev(igv,ncm(2))= -a2(ncm(3))*ev(igv,ncm(3))/a2(ncm(2))
      ev(igv,ncm(1))= (-a1(nrm(1),ncm(3))*ev(igv,ncm(3))-a1(nrm(1),ncm(2))*ev(igv,ncm(2)))/a1(nrm(1),ncm(1))
      chksum=0
      for j=1 to 3
         chksum=chksum+a1(nrm(3),ncm(j))*ev(igv,ncm(j))
      next j
      if asum(nrm(3))>0 then
        ierr=abs(chksum)/asum(nrm(3)) > 1e-5
      else
        ierr=abs(chksum) > 1e-5
      end if
      if not ierr then exit for
      ev(igv,ncm(3))=0   'else try zeroing out 3rd column
    next i
    if ierr then exit sub
    vnorm=fnamax(1e-5,sqr(ev(igv,1)*ev(igv,1)+ev(igv,2)*ev(igv,2)+ev(igv,3)*ev(igv,3)))
    if ev(igv,3)>=0 then ix=1 else ix=-1  'to move into lower hemisphere
    for i=1 to 3 : ev(igv,i)=ix*ev(igv,i)/vnorm : next i
  next igv
END SUB
'========================================================================

