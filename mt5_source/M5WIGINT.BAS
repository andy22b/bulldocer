sub m5wigint (tmin,tmax,t(1),a(1),nin,dt,a2(1),nout) static

'this subroutine uses a weighted-average-slope interpolator to interpolate 
'unevenly digitized data points stored in arrays T() and A() with NIN values. 
'Finds the amplitude Y() every DT seconds starting at TMIN until TMAX. NOUT is 
'the number of output points in Y(). 
'see Wiggins, RA, Interpolation of Digitized Curves, Bull. Seismol. Soc. 
'   America, 2077-2081 ,1976

shared nintmax

   epsilon = .0001 : nout=0
   j5=1 : t1=tmin : tnin=tmax

   for x=t1 to tnin step dt

     for j=j5 to nin
       if t(j)>x then goto WIG1
     next j

WIG1: j=j-1 : j5=j
     tj=t(j) : aj=a(j) 
     if j>1 then tjm1=t(j-1) : ajm1=a(j-1)
     if j<nin then tjp1=t(j+1) : ajp1=a(j+1)
     if j+1<nin then tjp2=t(j+2) : ajp2=a(j+2) 
     d=x-tj

     if d=0 then y=a(j) : goto WIG4
       h=tjp1-tj : d1=x-tjp1 : d4=(ajp1-aj)/h
       u4=d4 : a4=d4
       if j<>1 then d4=(aj-ajm1)/(tj-tjm1)
       if nin=j+1 then goto WIG3
       u4=1e5 : if tjp2=tjp1 then goto WIG3
       u4=(ajp2-ajp1)/(tjp2-tjp1)

WIG3: w1=10000. : w=10000. : w2=10000.
     if epsilon<=abs(d4) then w1=1/abs(d4)
     if epsilon<=abs(a4) then w=1/abs(a4)
     if epsilon<=abs(u4) then w2=1/abs(u4)
     r1=aj*(d1^2/h^2+2*d*d1^2/h^3)+ajp1*(d^2/h^2-2*d1*d^2/h^3)
     r2=(w1*d4+w*a4)/(w1+w)*d*d1^2/h^2
     r3=(w*a4+w2*u4)/(w2+w)*d^2*d1/h^2
     y=r1+r2+r3
     goto WIG4

WIG4: nout=nout+1 : a2(nout)=y
      if nout = nintmax then goto WIG5
     next x
WIG5:
   end sub

'===================================================================
sub m5detrnd (nin,dt,y(1),a,b) static

' detrend by fitting line through data and removing it
'  Y() is the input amplitude array with NIN values sampled at DT increment
'  A is the Y-intercept and B is the slope of the best-fit line

     sumx=0: sumy=0: sumxy=0: sumx2=0: sumy2=0
     for i=1 to nin
      x=(i-1)*dt : y=y(i)
      sumx=sumx+x : sumy=sumy+y
      sumxy=sumxy+x*y : sumx2=sumx2+x*x
      sumy2=sumy2+y*y
     next
      t2=sumx2-sumx*sumx/nin
      t3=sumxy-sumx*sumy/nin
      t4=sumy2-sumy*sumy/nin
      b=t3/t2
      a=((sumx2*sumy-sumx*sumxy)/nin)/t2
     for i=1 to nin
      ycor=a+b*(i-1)*dt
      y(i)=y(i)-ycor
     next
end sub
'===================================================================


