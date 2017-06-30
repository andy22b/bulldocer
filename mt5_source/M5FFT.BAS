     sub m5fft( work(1),n,isign) static
'
'    The Cooley-Tukey fast fourier transform in Microsoft Basic. 
'   transform(k) = sum(work(j)*exp(isign*2*pi*sqrt(-1)*(j-1)*(k-1) /n)), 
'   summed over all j and k from 1 to n.  work is a one-dimensional 
'   complex array (i.e., the real and imaginary parts are adjacent in 
'   storage) whose length n=2^k, k.ge.0 (if necessary, append zeroes to 
'   the data).  isign is +1 or -1.  if a -1 transform is followed by a +1 
'   one (or vice versa) the original data reappear, multiplied by n.  
'   transform values are returned in array work, replacing the input.  the 
'   time is proportional to n*log2(n), rather than the naive n^2. accuracy 
'   is also greatly improved, the rms relative error being bounded by 
'   6*sqrt(2)*log2(n)*2^(-b), where b is the number of bits in the 
'   floating point fraction.  written by norman brenner of mit lincoln 
'   laboratory.  july 1967.  this is the shortest version of the fft known 
'   to the author.  faster programs four2 and fourt exist that operate on 
'   arbitrarily sized multidimensional arrays. see--ieee audio 
'   transactions (june 1967), special issue on fft. 
'   adapted to Microsoft Basic by Rob McCaffrey, 1986
'   double precision theta,sinth,wstpr,wstpi,wr,wi,tempr 

  defint i-n
      ip0=2
      ip3=ip0*n
      i3rev=1
    for i3=1 to ip3 step ip0
      if(i3-i3rev) < 0 then
        tempr=work(i3)
        tempi=work(i3+1)
        wi3r=work(i3rev)
        wi3i=work(i3rev+1)
        work(i3)=wi3r
        work(i3+1)=wi3i
        work(i3rev)=tempr
        work(i3rev+1)=tempi
      end if

      ip1=ip3/2
f2:   if(i3rev-ip1)<=0 then goto f3
      i3rev=i3rev-ip1
      ip1=ip1/2
      if(ip1-ip0)>=0 then goto f2
f3:   i3rev=i3rev+ip1
    next i3
      ip1=ip0
f4:   if(ip1-ip3)>=0 then goto f5
      ip2=ip1*2
      theta=6.2831853071795865d0/(isign*ip2/ip0)
      sinth=sin(theta/2.)
      wstpr=-2.*sinth*sinth
      wstpi=sin(theta)
      wr=1.
      wi=0.
    for i1=1 to ip1 step ip0
    for i3=i1 to ip3 step ip2
      i2a=i3
      i2b=i2a+ip1
      wi2r=work(i2b) : wi2i=work(i2b+1)
      tempr=wr*wi2r-wi*wi2i
      tempi=wr*wi2i+wi*wi2r
      wi2r=work(i2a)
      wi2i=work(i2a+1)
      work(i2b)=wi2r-tempr
      work(i2b+1)=wi2i-tempi
      work(i2a)=wi2r+tempr
      work(i2a+1)=wi2i+tempi
    next i3
      tempr=wr
      wr=wr*wstpr-wi*wstpi+wr
      wi=wi*wstpr+tempr*wstpi+wi
    next i1
      ip1=ip2
      go to f4
f5:  end sub 
'======================================================================
sub invfft(nt,dt,fresp(2),ss(1)) static
'
'   inverse fourier transform
'   Nabelek
'
  nwork=2055
  redim work(nwork)
  isign=1
      n3=nt+nt+3
      nf=nt/2+1
      df=1./(nt*dt)
    for i=1 to nf-1
      i2=i+i
      i1=i2-1
      work(i1)=fresp(i,1)
      work(i2)=fresp(i,2) 
      work(n3-i1)=-work(i2)
      work(n3-i2)=work(i1)
    next 
      work(nt+1)=fresp(nf,1)
      work(nt+2)=0.
  call m5fft(work(),nt,isign)
    for i=1 to nt
      w=work(i+i-1)
      ss(i)=w*df
    next
end sub
'======================================================================
sub forfft(n4,dt,fresp(2),ss(1)) static
'
'   forward fourier transform
'   Nabelek

  nwork=2055
  redim work(nwork)
  isign=-1

    for i=1 to n4
      i2=i+i
      i1=i2-1
      work(i1)=ss(i)
      work(i2)=0.0
    next 

  call m5fft(work(),n4,isign)
    for i=1 to n4
      i2=i+i
      i1=i2-1
      fresp(i,1)=work(i1)
      fresp(i,2)=work(i2)
    next
end sub
'======================================================================
sub BUT3P(nfreq,dfreq,fresp(2),fcuts) static
'------- filtering subroutine 3 pole Butterworth

  if fcuts<=0 then exit sub

  fcut=1.0/fcuts
  c60=.5 : s60=.866025
  pre=c60*fcut : pim=s60*fcut
  p1r=pre  : p1i=pim : p2r=pre : p2i=-pim : p3r=fcut : p3i=0

 for i=1 to nfreq
  freq=i*dfreq
  sr=0 : si=freq
     call comadd(sr,si,p1r,p1i,sp1r,sp1i)    ' s+p1
     call comadd(sr,si,p2r,p2i,sp2r,sp2i)    ' s+p2
     call comadd(sr,si,p3r,p3i,sp3r,sp3i)    ' s+p3

     call comult(sp1r,sp1i,sp2r,sp2i,sp12r,sp12i)     '(s+p1)*(s+p2)
     call comult(sp12r,sp12i,sp3r,sp3i,spallr,spalli) '             *(s+p3)
     call cominv(spallr,spalli, spinvr,spinvi)        ' invert

     call compwr(sr,si,3,s3r,s3i)                     's^3
     call comult(s3r,s3i,spinvr,spinvi,freal,fim)
     frr=fresp(i,1) : fri=fresp(i,2)
     call comult(frr,fri,freal,fim,forr,fori)
     fresp(i,1)=forr : fresp(i,2)=fori
 next
 end sub
'======================================================================

