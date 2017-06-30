sub m5synio (io, idat, ista, nin, n1, a(1), amax, asq) static
  shared seistmp$, seistmp(), nsamax, ndat(), ipol(), index()
  shared kphtype(), phmag(), invflag, refmag()

'--- seismogram I/O
' IO = 0 seismogram output, =1 for input
' IDAT = 0 for data, =1 for synthetics
' First NSAMAX records are stored in array SEISTMP, 
'    rest are in external file SEISTMP$

  if ista<=nsamax then 
     nrec1&= 2&*index(ista) + idat*ndat(ista) + n1
  else
     nrec1&= 2&*(index(ista) - index(nsamax+1)) +idat*ndat(ista) + n1
  end if

if io=1 then
   if nin=0 then exit sub

'-- read NIN points starting at point N1 from file or array

'--from array
   if ista<=nsamax then 
      for i=1 to nin :a(i)=seistmp(i+nrec1&) : next 

'--from random access file
   else
      open seistmp$ as 3 len=4 : field #3, 4 as io$
       get #3, nrec1&+1 : a(1)=cvs(io$)
       for i=2 to nin : get #3: a(i)=cvs(io$) : next 
      close 3 
   end if

'-- flip polarity of data, scale to reference magnification, and get max
   if (idat=0) then
     if (ipol(ista)=-1) then pol=-1 else pol=1
   else 
     pol=1
   end if

'-- get maximum and SQR(sum of squares/nin)
 amax=0 : asq=0 :  x=refmag(ista)*pol    '''''' x=phmag(kphtype(ista))*pol
  for i=1 to nin 
    a(i)=a(i)*x
    amax=fnamax(amax,abs(a(i)))
    asq=asq + a(i)*a(i)
  next
    amax=fnamax(1.0e-5,amax)
    asq=fnamax(1.0e-5,sqr(asq))

else  '--- write NIN points starting at point N1 to file or array, 
   if nin=0 then exit sub

'--to array
   if ista<=nsamax then 
      for i=1 to nin : seistmp(i+nrec1&)=a(i) : next 

   else
'--to random access file
      open seistmp$ as 3 len=4 : field #3, 4 as io$
       lset io$=mks$(a(1)) : put #3, nrec1&+1
       for i=2 to nin : lset io$=mks$(a(i)) : put #3 : next 
      close 3
   end if
end if
end sub

'==============================================================================
sub m5esio (io,ista,nph,nin,a(1)) static
  shared gfxixt$, es(), nsamax, ndat(), indexes&()

'--- Elementary seismogram I/O
' IO = 0 seismogram output, =1 for input
' First NSAMAX records are stored in array ES,
'    rest are in external file GFXIXT$

  if ista<=nsamax then
     nrec1&= indexes&(ista) + (nph-1)*ndat(ista)
  else
     nrec1&= indexes&(ista) - indexes&(nsamax+1) +(nph-1)*ndat(ista)
  end if

  if io=1 then

'---- read NIN points from file or array,

'--from array
   if ista<=nsamax then
      for i=1 to nin :a(i)=es(i+nrec1&) : next

'--from random access file
   else
          open gfxixt$ as 3 len=4 : field #3, 4 as io$
           get #3, nrec1&+1 : a(1)=cvs(io$)
           for i=2 to nin : get #3: a(i)=cvs(io$) : next
          close 3
   end if

else  '--- write NIN points starting at point N1 to file or array,


'--to array
   if ista<=nsamax then
      for i=1 to nin : es(i+nrec1&)=a(i) : next

   else
'--to random access file
          open gfxixt$ as 3 len=4 : field #3, 4 as io$
           lset io$=mks$(a(1)) : put #3, nrec1&+1
           for i=2 to nin : lset io$=mks$(a(i)) : put #3 : next
          close 3
   end if
end if

end sub
'==============================================================================

sub m5tmpio (io, ista, nph, nin, a(1), files$) static
  shared indexes&(), ndat()

'--- time series I/O for elementary seismograms
' IO = 0 series output, =1 for input
' phase NPH from station ISTA, read NIN points from FILES$

 nrec1&= indexes&(ista) + (nph-1)*ndat(ista) +1

  open files$ as 2 len=4 : field #2, 4 as IO$

  if io=1 then
'---- read NIN points from random access file
         get #2, nrec1& : a(1)=cvs(io$)
         for i=2 to nin : get #2: a(i)=cvs(io$) : next
        close 2
  else
   '--- write NIN points to file
         lset io$=mks$(a(1)) : put #2, nrec1&
         for i=2 to nin : lset io$=mks$(a(i)) : put #2 : next
        close 2
  end if
end sub
