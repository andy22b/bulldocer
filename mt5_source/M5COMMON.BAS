'  M5COMMON.BAS
'  commons for MT5.BAS
'  revised Mar, 1992
'  increased nesmax limit to handle broadband data,   PZ-11/9/93
'  added COMMON maxnes to track actual maximum of nes PZ-5/10/93
' This include file should be first and should be preceded by:
' DEFINT I-N: OPTION BASE 1

CONST NSTAMAX = 50  'seismogram limit also used in M5COMSCR

'--- screen parameters
rem $INCLUDE: 'M5COMSCR'

'------- array dimensions
'--- fixed sizes, do not adjust, for geometrical spreading
   DIM CGS(16), DGS(16)

CONST NOFFMT = 7 'offset of moment tensor to double couple elements in sorpar()

'---related to number of inversion parameters (=13)
'--  dimensions are NTAUMAX=28, NSORMAX=6
   DIM ATAU(28,6), DAMP(13,6), SORPAR(13,6), ERPAR(4,6), STFDAMP(6)
   DIM TAM(6), ERMOM(6), NTAU(6), ERTAU(28,6), ERT95(6), ERMECH(6,6)
   DIM TAMISO(6),EVCT(3,3,6)
'-- limits and controls for variable range parameters
   DIM VARBEG(4,6), VAREND(4,6), VARINC(4,6), VRP$(6)

'---related to number of layers in structure, NLYRMAX=2
   DIM ZTOP(2), ALPH(2), BETA(2), RHO(2)

'---related to number of seismograms NSTAMAX = 50
   DIM WT(NSTAMAX), RAPAR(NSTAMAX), IPH(NSTAMAX), IFLAG(NSTAMAX), PTOA(NSTAMAX)
   DIM STAVAR(NSTAMAX), REFVAR(NSTAMAX), DEL(NSTAMAX), AZD(NSTAMAX), IPF1(NSTAMAX)
   DIM IPF2(NSTAMAX), IPHES(NSTAMAX), TSTIME(NSTAMAX), THTIME(NSTAMAX), IXFINE(NSTAMAX)
   DIM IYFINE(NSTAMAX), DEXT(NSTAMAX), GCOR(NSTAMAX), NINV(NSTAMAX), ISPOS(NSTAMAX,3)
   DIM XMAG(NSTAMAX),  INST(NSTAMAX), INPTS(NSTAMAX), WT0(NSTAMAX), TSTAR(NSTAMAX)
   DIM DT(NSTAMAX), TAU(3,NSTAMAX), NES(NSTAMAX), NZERO(NSTAMAX), NPOLE(NSTAMAX)
   DIM CALCONS(NSTAMAX), ZEROPOLE(2000), INDEXZP(NSTAMAX), HPFILT(NSTAMAX), NODAT(NSTAMAX)
   DIM SYNRMS(NSTAMAX), DATRMS(NSTAMAX), INSTIN(NSTAMAX), INSFILE$(NSTAMAX)

'--- number of phase data types, =7
   DIM KPHFLAG(7), PHTIME(7), PHSCL(7), PHRLEN(7), TWINDOW(7)
   DIM PHDT(7), PHTSTAR(7), IPHLONG(7), IPWAVE(7)
   DIM DATVAR(7), RESVAR(7)

'--- number of phases, =2
   DIM KPLOT(2), KFOC(2), IFOCX0(2), IFOCY0(2), IFOCRAD0(2), IFOCPOS(2)

'--- plot box dimensions
   DIM BOXLEN0(2), BOXHIGH0(2)

'--****** COMMONS ********
'--- limiters and constants
 COMMON NES(), NESSUM, DT, ITF, XMAG0, DEL0, ICONT, IPCFLAG
 COMMON ISTFLAG, IPARFLG, ISPFLAG, ILPFLAG, IGDSNFLAG
 COMMON NODAT, ICHOICE, IMT5, MTFLAG, LASTIT, MAXNES
'--- screen plotting
 COMMON NPPS, NPIX, IPLOFF, KPDAT, KPSYN
'--- structure parameters
 COMMON NREVCOM, ALPHR, BETAR, RHOR, ALPHF, BETAF, RHOF
 COMMON WTHK, ALPHW, BETAW, RHOW
 COMMON CTHK, ZTOP(), ALPH(), BETA(), RHO(), NLYR, LYRS
 COMMON ALPH2, BETA2, RHO2, VPS
'--- data parameters
 COMMON FHDR$, DF$, ROF$, NSTA, IPF2(), ISPOS(), INST(), HYPOC$
 COMMON STACODE$, IPH(), DEL(), AZD(), WT(), IPF1(), DEXT()
 COMMON RAPAR(), PTOA(), STAVAR(), REFVAR(), GCOR(), NINV(), DT()
 COMMON IFLAG(), IPHES(), XMAG(), WT0()
 COMMON TSTIME(), THTIME(), HPFILT(), NODAT(), IXFINE(), IYFINE()
 COMMON INPTS(), TSTAR(), SYNRMS(), DATRMS(), INSTIN(), INSFILE$()
'--- instrument and phase type parameters
 COMMON NZERO(), NPOLE(), CALCONS(), ZEROPOLE(), INDEXZP()
 COMMON KPHFLAG(), PHTIME(), PHSCL(), PHRLEN(), TWINDOW()
 COMMON PHDT(), PHTSTAR(), NPHTYPE, IPHLONG(), IPWAVE()
 COMMON DATVAR(), RESVAR()
'--- herrin tables and geometrical spreading curve
 COMMON CGS(), DGS()
'--- source parameters
 COMMON TAU(), TAU, TAU0, NORPS, NTAU(), TAM(), SORPAR(), DEPTOT, TAMISO()
 COMMON NSOR, ATAU(), STFDAMP(), ISORM, ERTAU(), VRUPT, AZRUPT, PLRUPT, EVCT()
'--- inversion parameters
 COMMON THETA, DAMP(), STFMIN, DAMPFAC, DAMPMAX, VARMIN
'--- variable range parameters
 COMMON VARBEG(), VAREND(), VARINC(), VRP$()
'--- statistics variables
 COMMON VARCOR, ERPAR(), ERMOM(), KRMSC, ERT95(),VARDAT, VARRES, VARNOWT
 COMMON ERMECH()
'--- plot parameters
 COMMON KTSC, KDLY, LABPOS, NORMAMP, LABCODE
 COMMON SINCR, LSYN, IAMPOS, GS$, KPLOT(), IPLFORM, MTPLFORM
 COMMON TITLE$, KAPPEND, SPCSP, SPCLP
 COMMON KPMECH, KLETTER, KPLABEL, ISPOSIT, YAPPOFF, I2PAGE, GSROT
 COMMON KFOC(), IFOCX0(), IFOCY0(), IFOCRAD0(), IFOCPOS()
 COMMON BOXLEN0(), BOXHIGH0(), IGRAPH, SCRORIG, XMAGP, XOFFHP, YOFFHP, IGSFN
'--- dynamic arrays in common
 COMMON SP()
'--- I/O parameters
 COMMON GF$, GFXI$
'--- Variables set in SETUP required by MT5MAIN module  (added Jan 1991)
 COMMON MKDAT, MKGRN, MKTAU, MKINST, ISET, RMSCALE
'--- end of COMMON block



' DECLARES for M5SUB.BAS
rem $INCLUDE: 'M5FNDECL'

'-------- limiters
'  NSYNMAX =15000  'data plus synthetic points stored in array SEISTMP
  NPARMAX =   40  'parameters for inversion at one time
  NESMAX  = 2896  'points in one seismogram, increased from 500 11/9/93 PZ
  NINVMAX =  500  'points to invert for one seismogram, not used anywhere in code
  NTAUMAX =   25  'source time function elements
  NLYRMAX =    2  'layers in structure
  NPHSMAX =    2  'phases (i.e. P and SH)
  NSORMAX =    6  'number of sources

'-- dimension SS to 2*NESMAX, ODAT to NESMAX AS STATIC ARRAYS
 DIM SS(2*NESMAX),ODAT(NESMAX)

'-- functions handling phase and instrument types

'-- return station name and * if not being used in the inversion
DEF FNSTANM$(I)
  y$=" "
  if not iflag(i) or ninv(i)=0 or wt(i)=0 or nodat(i) then y$="*"
  fnstanm$=ucase$(mid$(STACODE$,1+(i-1)*4,4))+y$
end def

'-- return 4-letter code for phase-data type
DEF FNKPHLAB$(I)
   fnkphlab$="    " : y$="    "
   if (i mod 10)>0 and (i mod 10)<8 then y$=mid$("WWLPWWLSWWSPGDLPGDLSGDSPGDBP",1+((i mod 10)-1)*4,4)
   fnkphlab$=y$ 
'-- if non-standard instrument add prefix of instrument code
'   if int(i/10)=0 then fnkphlab$=y$ : exit def
'   fnkphlab$=mid$("TSCOPZ",1+(int(i/10)-1)*2,2)+mid$(y$,3,2)
end def

'-- single letter code for phase-data type
DEF FNPHMARK$(I)= mid$("wwsddxb",i,1)

'-- flag data types
DEF FNPWAVE(i) = (kphtype(i)=1) OR (kphtype(i)=3) OR (kphtype(i)=4) OR (kphtype(i)=6) OR (kphtype(i)=7)
DEF FNILONG(i) = (kphtype(i)=1) OR (kphtype(i)=2) OR (kphtype(i)=4) OR (kphtype(i)=5) OR (kphtype(i)=7)
DEF FNISHORT(i) = (kphtype(i)=3) OR (kphtype(i)=6)
DEF FNIBB(i) = (kphtype(i)=7)
DEF FNINGDSN(i) = (kphtype(i)=4) OR (kphtype(i)=5) OR (kphtype(i)=6) OR (kphtype(i)=7)
DEF FNIWWSSN(i) = (kphtype(i)=1) OR (kphtype(i)=2) OR (kphtype(i)=3)
