'M5SUB1.BAS
'Compile and link separately with MT5IT and MT5MAIN

DEFINT I-N: OPTION BASE 1

' additional common block needed by m5iresp
dim zero(100), pole(100)
COMMON /INSTRESP/ nzero,npole,zero(),pole(),ts,ds,tg,dg,tgain

rem $INCLUDE: 'M5COMMON'
rem $INCLUDE: 'M5FUNCTN'

rem $INCLUDE: 'm5rapars'
rem $INCLUDE: 'm5iresp'

