'-- This first group of COMMON variables needs to be the same in all files
'    that use M5SUB.

'-- six text color combinations: foreground & background for each

 TYPE COLORTYPE
    IFORE AS INTEGER
    IBACK AS INTEGER
 END TYPE
 DIM KTXTCOL(6) AS COLORTYPE,KGRCOL(0 to 3)
 CONST KTRUE = -1, KFALSE = 0, PII = 3.14159265 , PIIHALF = 1.57079633 , SQR2=1.41421356
 CONST KCOLDAT=1, KCOLSYN=2, KCOLBAK=3      'elements of kgrcol

 COMMON KEGA, KSCRX, KSCRY, KGRCOL(), KPENHP
 COMMON KSCREEN, KOKCOL, KTXTCOL() AS COLORTYPE
 COMMON SCRATIO, NCOLCHAR, NROWCHAR, PRGNM$, VERSION$, SYNDIR$, DRIVE$
 COMMON ITIMEROW, ITIMECOL  'position for display of time

'---related to number of seismograms NSTAMAX = 50
DIM ndat(NSTAMAX), ipol(NSTAMAX), INDEX(NSTAMAX)
DIM INDEXES&(NSTAMAX), KPHTYPE(NSTAMAX), REFMAG(NSTAMAX)
'--- number of phase data types, =7
DIM PHMAG(7)
'--- limiters and constants
COMMON nsamax, INVFLAG, nsynmax
'--- data parameters
COMMON ndat(), ipol(), index(), indexes&(), refmag()
'--- I/O parameters
COMMON seistmp$, gfxixt$
'--- dynamic arrays in common
COMMON ES(1),SEISTMP(1)
'--- instrument and phase type parameters
COMMON kphtype(), phmag()
'--- ipmi can be used to allow printout of debug info 
' (currently activated by an environment variable of the same name)
COMMON IPMI

'--- end of variables defined as COMMON in M4SUB
