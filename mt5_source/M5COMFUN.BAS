'M5COMFUN
'======================================================================
' complex arithmetic
  
sub comadd (ar,ai,br,bi,cr,ci) static
' adds complex numbers A + B = C
cr=ar+br
ci=ai+bi
 end sub

sub comult (ar,ai,br,bi,cr,ci) static
' multiplies complex numbers A * B = C
cr=ar*br-ai*bi
ci=ar*bi+ai*br
 end sub

sub comexp (ar,ai,br,bi) static
' complex exp,  B = EXP(A)
br=exp(ar)*cos(ai)
bi=exp(ar)*sin(ai)
 end sub

sub cominv (ar,ai,br,bi) static
' complex inverse of A,  B=1/A
  x=ar*ar+ai*ai
  br= ar/x
  bi= -ai/x
 end sub

sub compwr (ar,ai,n,br,bi) static
' complex power B=A^N, N is integer
  if n=1 then br=ar : bi=ai : exit sub
  cr=ar : ci=ai
     for i=2 to n
       call comult (cr,ci,ar,ai,br,bi)
       cr=br : ci=bi
     next
 end sub

