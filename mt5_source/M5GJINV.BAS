sub m5gjinv (ncol,b(2),z(1),determ,ising) static
'---- Subroutine for Gauss-Jordan matrix inversion
'         Solves Bx=Z, output is Binv in B, and x in Z
'  ISING becomes true (-1) for singular matrix

' By McCaffrey and Abers
' Modified by Zwick to include double precision & scaled pivoting - Sept 1991

  REDIM INDX(NCOL,3),VMAX#(NCOL),A#(NCOL,NCOL),ZD#(NCOL)
  ISING=0
  DBLMIN#=1D-12

'---  make B# the full matrix
  FOR L=1 TO NCOL-1 : FOR K=L+1 TO NCOL : B(L,K)=B(K,L) : NEXT K,L

'--- get scaling factor for each row and overall max element
  T2=0
  FOR J=1 TO NCOL
    T1=0
    FOR IT=1 TO NCOL
      TST=ABS(B(J,IT))
      IF TST>T2 THEN
        T1=TST
        T2=TST : I3=J : I4=IT
      ELSEIF TST>T1 THEN
        T1=TST
      END IF
    NEXT IT
    VMAX#(J)=T1
    IF VMAX#(J)<DBLMIN# THEN GOTO GJ6
  NEXT J

'--- duplicate B & Z into double precision A# & ZD#
  FOR I=1 TO NCOL
    ZD#(I)=Z(I)
    FOR K=1 TO NCOL : A#(I,K) = B(I,K) : NEXT K
  NEXT I

'--- find largest element (maximum scaled full pivoting)
  FOR I=1 TO NCOL   
    IF I>1 THEN 'already have max overall element for I=1
      B1#=0
      FOR J=1 TO NCOL
        IF INDX(J,3)<>1 THEN
          t2#=B1#*VMAX#(j)
          FOR IT=1 TO NCOL
            IF INDX(IT,3)>1 THEN goto gj6
            IF INDX(IT,3)<>1 THEN
              IF t2#<=ABS(A#(J,IT)) THEN
                I3=J: I4=IT
                t2#=ABS(A#(J,IT))
              END IF
            END IF  'INDX(IT,3)<>1
          NEXT IT
          B1#=t2#/VMAX#(j)
        END IF  'INDX(J,3)<>1
      NEXT J
    END IF

    INDX(I4,3)=INDX(I4,3)+1
    INDX(I,1)=I3
    INDX(I,2)=I4

'--- interchange rows to put pivot on diagonal
    IF I3<>I4 THEN
      FOR L=1 TO NCOL : SWAP A#(I3,L),A#(I4,L) : NEXT L
      SWAP ZD#(I3),ZD#(I4) : SWAP VMAX#(I3),VMAX#(I4)
    END IF  'I3<>I4

'--- divide pivot row by pivot element
    IF ABS(A#(I4,I4))<DBLMIN# THEN GOTO gj6
    PIVOTINV#=1.D0/A#(I4,I4)
    A#(I4,I4)=1
    FOR L=1 TO NCOL: A#(I4,L)=A#(I4,L)*PIVOTINV#: NEXT L
    ZD#(I4)=ZD#(I4)*PIVOTINV#

'---reduce nonpivot rows
    FOR L1=1 TO NCOL
      IF L1<>I4 THEN
        T#=A#(L1,I4)
        A#(L1,I4)=0
        FOR L=1 TO NCOL : A#(L1,L)=A#(L1,L)-A#(I4,L)*T# : NEXT L
        ZD#(L1)=ZD#(L1)-ZD#(I4)*T#
      END IF  'L1<>I4
    NEXT L1
  NEXT I

  '---- interchange columns back to restore original order
  FOR I=NCOL TO 1 STEP -1
    IF INDX(I,1)<>INDX(I,2) THEN
      FOR IT=1 TO NCOL : SWAP A#(IT,INDX(I,1)),A#(IT,INDX(I,2)) : NEXT IT
    END IF
  NEXT I

  FOR IT=1 TO NCOL
    IF INDX(IT,3)<>1 THEN GOTO GJ6
  NEXT IT

  FOR I=1 TO NCOL
    Z(I)=ZD#(I)
    FOR K=1 TO NCOL : B(I,K)=A#(I,K) : NEXT K
  NEXT I
  EXIT SUB

gj6:  ISING=-1 'singular matrix
END SUB

