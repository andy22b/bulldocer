' MT5 (MT5.00) Startup module

' teleseismic body waveform inversion program for PC's written by
'   Rob McCaffrey and Geoff Abers
'              SYN3.01 finished April 1988,
'     revision SYN3.02 June 1989
'     revision SYN3.04 June 1990
'     revision SYN4.00 September 1990 (with Peter Zwick)
'     revision MT4.00 (includes inversion for moment tensor) Nov 1991 
'                       by Peter Zwick
'     revisions V1.1 and V1.11 are minor bug fixes
'     revision MT5 V1 May - 1994 by Peter Zwick includes grid search
'        and input file sizes up to ~350 Kb for more broadband data
'   Link with m5sub and specify QB as library during link

DEFINT I-N: OPTION BASE 1

rem $INCLUDE: 'M5COMMON'
rem $INCLUDE: 'M5FUNCTN'


'------ start --------------------------------------------------------------
prgnm$="MT5"
authors$="P. Zwick, R. McCaffrey, & G. Abers"
version$="1.0"

if (not icont) and (not iset) then

  call welcome
  call scprint(13,3,"Setting up...")
  gosub SETUP

elseif ichoice=101 then           '-- E for edit \MT5.SET

  call GETSETUP
  call READSET

else                              '   just reinitialize all parameters

  CLEAR   ' initializes all variables to zero
  prgnm$="MT5"
  gosub SETUP

end if

chain syndir$+"MT5MAIN"

'--------SETUP---------
SETUP:  '--- set up values and defaults

nsynmax=15000
redim seistmp(nsynmax)

'--- read setup file

alphr=6.0   : betar=3.4   : rhor=2.7    '--- receiver structure
alphf=0.335 : betaf=0.001 : rhof=0.001  '--- air
alphw=1.50  : betaw=0.001 : rhow=1.03   '--- water

varmin = 0.      'minimum variance reduction to end iterating
varcor=1.0e8     'correct variance from cm^2 to micron^2 for printing
stfmin=-0.05     'minimum value of STF element
dampfac=100.     'damping increase when STF element amplitude < STFMIN
dampmax=1.0e37   'maximum damping allowed on STF element
theta=0.0        'least squares damper
nrevcom=2        'number of complete reverberations for upper layer 
del0=40.         'reference delta for seismogram distance corrections

imt5=ktrue      'flag that program was started with MT5.EXE
invflag=ktrue    'invert seismogram amplitudes
istflag=kfalse   'T-statistics off
ipcflag=ktrue    'positivity constraint on
mkdat=ktrue      'flag to read in data
mkgrn=ktrue      'flag to calculate Greens functions
mktau=ktrue      'flag to convolve TAU with elementary seismograms
mkinst=ktrue     'flag to calculate instrument responses
iset=ktrue       'flag this setup routine
nodat=ktrue      'no data has been read yet

'--- screen plotting defaults 
npps=4           'number of seismograms plotted per screen
npix=0           'number of screen pixels per second for plotting
iploff=kfalse    'seismograms not offset on screen
kpdat=ktrue      'plot data on screen
kpsyn=ktrue      'plot synthetics on screen

'--- hardcopy plotting defaults
ktsc=ktrue       'plot time and amplitude scales
normamp=kfalse   'normalize seismograms in plots
rmscale=ktrue    'scale seismograms to RMS if NORMAMP on
kappend=kfalse   'appending to .PLT file
kletter=ktrue    'plot letters for station positions in focal sphere
kplabel=ktrue    'write label on plot
sincr=1.5        'increment in mm of dashes in synthetic seismograms
labpos=0         'position of title
isposit=31       'position of time scales and STF
iampos=0         'position of amplitude axes
iplform=1        'plot format number
boxlen0(2)=35    'length of plot box for format 2
boxhigh0(2)=30   'height of plot box for format 2
yappoff=30       'offset between lines of seismograms in format 2
igraph=1         'output in Golden Software PLOTCALL format

for i=1 to 2 
 kfoc(i)=ktrue   'plot focal sphere, i=1 for P, i=2 for SH
 ifocpos(i)=0    'position for focal sphere
 ifocx0(i)=125   'X position for focal sphere
 ifocy0(i)=90    'Y position for focal sphere
 ifocrad0(i)=30  'radius of focal sphere
next

'-- default source parameter values, set for MT5TEST.INV
sorpar(1,1)=fnrad(20.) : sorpar(2,1)=fnrad(30.) : sorpar(3,1)=fnrad(40.)
for j=1 to nsormax
 ntau(j)=1       'number of elements in STF
 atau(1,j)=1     'amplitude of STF element
 tam(j)=1.0e18   'seismic moment for event in Nm
 sorpar(4,j)=20. 'depth of source
 stfdamp(j)=1.   'STF damper
   for i=1 to 4 : damp(i,j)=1.0 : next
   for i=5 to 7 : damp(i,j)=0.0 : next
   for i=8 to 13: damp(i,j)=1.0 : next
next

'--- parameter defaults 
tau0=1.0     'length of STF triangular element in s
nsor=1       'number of sources
isorm=1      'start on source #1
vrupt=0.0    'rupture velocity in km/s
azrupt=0.0   'rupture azimuth in degrees
plrupt=0.0   'rupture plunge in degrees
nsta=0       'number of stations
ztop(1)=0.0  'top of layer in km
alph(1)=6.5  'P velocity
beta(1)=3.7  'S velocity
rho(1)=2.8   'density of layer
deptot=20.   'depth of source in halfspace layer
wthk=0       'thickness of water layer
nlyr=1       'number of layers in structure

'--- setup defaults for data-phase types
'-- 6 types are used now, the 7th is reserved for broadband P data
'     1=WWSSN LP P, 2=WWSSN LP SH, 3=WWSSN SP P, 4=GDSN LP P, 5=GDSN LP SH, 
'     6=GDSN SP P, 7=BB P

nphtype=7
for i=1 to nphtype: phscl(i)=1 : next     'plotting scale factor

'--- time lengths to be plotted
  phtime(1)=60. : phtime(2)=60. : phtime(3)=12. : phtime(4)=60.
  phtime(5)=60. : phtime(6)=12. : phtime(7)=60.

'--- time lengths to be calculated
  phrlen(1)=50. : phrlen(2)=60. : phrlen(3)=15. : phrlen(4)=60.
  phrlen(5)=70. : phrlen(6)=15. : phrlen(7)=50.

'--- default window lengths for phases
  twindow(1)=30. : twindow(2)=30. : twindow(3)=10. : twindow(4)=40.
  twindow(5)=50. : twindow(6)=10. : twindow(7)=30.

'--- reference magnifications
  phmag(1)=3000. : phmag(2)=3000. : phmag(3)=25000. : phmag(4)=6000.
  phmag(5)=6000. : phmag(6)=25000. : phmag(7)=3000.

'--- default dt
  phdt(1)=0.5 : phdt(2)=0.5 : phdt(3)=0.1 : phdt(4)=1.0
  phdt(5)=1.0 : phdt(6)=0.1 : phdt(7)=0.1

'--- default t*
  phtstar(1)=1.0 : phtstar(2)=4.0 : phtstar(3)=0.7 : phtstar(4)=1.0
  phtstar(5)=4.0 : phtstar(6)=0.7 : phtstar(7)=1.0

'--- phase types, P=true, SH=false
  ipwave(1)=ktrue : ipwave(2)=kfalse: ipwave(3)=ktrue: ipwave(4)=ktrue
  ipwave(5)=kfalse: ipwave(6)=ktrue : ipwave(7)=ktrue

'--- phase-data types, long-period=true, short-period=false
  iphlong(1)=ktrue : iphlong(2)=ktrue : iphlong(3)=kfalse: iphlong(4)=ktrue
  iphlong(5)=ktrue : iphlong(6)=kfalse: iphlong(7)=ktrue

'--- setup I/O files
fhdr$="*"
gf$=drive$+"gf.tmp"
gfxi$=drive$+"gfxi.tmp"
gfxixt$=drive$+"gfxixt.tmp"
seistmp$=drive$+"seis.tmp"

'---get geometrical spreading curve
'----- from Langston and Helmberger, GJRAS, 42, p. 124, 1975.
data 8.40E-5,7.95E-5,7.65E-5,7.40E-5,7.15E-5,6.90E-5,6.70E-5,6.35E-5
data 5.95E-5,5.60E-5,5.30E-5,5.00E-5,4.70E-5,4.40E-5,4.10E-5,3.85E-5
'--read constants
for i=1 to 16: read cgs(i): next 
'--read delta indices
data 30,32.5,35,37.5,40,42.5,45,50,55,60,65,70,75,80,85,90
for i=1 to 16: read dgs(i): next

call scpause(13,3)

return

'--------------------------------------------------------------------
rem $INCLUDE: 'm5rdset'

