'-- subroutines to handle reading travel-time tables, calculate travel times
'   and derivatives and ellipticity corrections
'======================================================================

sub readtime (ttfile$) static
'-------- read in file of stations and tables
'  # of depths, first delta, delta increment, # of deltas
'     NDEP        DEL0            DDELT         NDEL
'  DEP() - depths
'  TTIMTAB(i,j,k) - i=delta, j=depth, k=phase

shared nph, ttimtab(), ndep(), del0(), ddelt(), ndel(), dep()

   call cline(18,1)
   call scprint(18,3,"Reading travel-time tables ")
   if not ifopen("i",4,ttfile$,0,-1) then end
   print ttfile$+"..."
  nph=0
    while not eof(4)

'first line is header - not used
     line input #4, h$
     nph=nph+1

     input #4, ndep(nph), del0(nph), ddelt(nph), ndel(nph)
     for j=1 to ndep(nph)
      input #4, dep(j,nph)
      for i=1 to ndel(nph): input #4, ttimtab(i,j,nph) : next
     next

    wend
   close 4
   call cline(18,3)
end sub

'======================================================================

sub travtime ( iph, delta, depth, ptime, dtdh, dtdz ) static
'----------get travel time from table

shared nph, ttimtab(), ndep(), del0(), ddelt(), ndel(), dep()

 deltatmp=delta
 textra=0.0

 idel=int( (deltatmp-del0(iph)) /ddelt(iph) )+1
    if idel<1 then ptime=0 :dtdh=0 :dtdz=0: exit sub

    if idel>=ndel(iph) then 
      if iph=1 then pdiff=4.5 else pdiff=8.3
      idel=ndel(iph)-1
      textra = (deltatmp - idel*ddelt(iph)) * pdiff   'for Pdiff
      deltatmp=idel*ddelt(iph)
    end if

     H=ddelt(iph)
    dH=deltatmp/H-int(deltatmp/H)
     H=H*111.19   'convert to km for derivatives

 idep=ndep(iph)-1
    for i=1 to ndep(iph)-1
      if depth>=dep(i,iph) and depth < dep(i+1,iph) then idep=i
    next
    Z= dep(idep+1,iph)-dep(idep,iph)
   dZ=(depth-dep(idep,iph))/Z

 t0=ttimtab(idel,idep,  iph)    : t2=ttimtab(idel+1,idep,  iph)-t0
 t3=ttimtab(idel,idep+1,iph)-t0 : t4=ttimtab(idel+1,idep+1,iph)-t0

'-- travel time in seconds
 ptime=t0+ t2*dH + t3*dZ + (t4-t2-t3)*dH*dZ + textra

'-- horizontal slowness in seconds per km
 dtdh= t2/H + (t4-t3-t2)/H * dZ

'-- vertical slowness in seconds per km
 dtdz= t3/Z + (t4-t3-t2)/Z * dH

end sub

'======================================================================

sub ellipcor ( nph, delta, eplat, stalat, ecor ) static
'--------- ellipticity correction
height$="  7  7  7  6  5  3  2  0 -2 -4 -5 -7 -9-10-12-13-14-14-14"
dist$="0012233334444556677777788899999999999"

      l1 = int(abs(eplat)/5.0)
      hep = val(mid$(height$, l1*3+1,3))
      l2 = int(abs(stalat/5.0))
      h = val(mid$(height$, l2*3+1,3))
      l3 = int(delta/5.0)
      xfnd = val(mid$(dist$, l3 + 1,1))/100.
      ellip = xfnd*(hep + h)

  if nph=1 then adj=1.0 else adj=1.75
  ecor=ellip*adj
end sub
'======================================================================

