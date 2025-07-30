* SPICE modeli elementov za laboratorijske vaje pri predmetu AEV 
* in nekaj dodatkov za raziskavo avtomatske sinteze topologij AEV (Ziga)


********
* Model opampa z vgrajenim napetostnim virom na vhodnih sponkah za merjenje vhodne upornosti.
* Napetostni vir bo nastavljen na fiksno napetost, zaporedno z njim pa spremenljiva upornost,
* ki bo imela vlogo stikala.
* 
* LM348 Opamp "MACROMODEL" SUBCIRCUIT
* (REV N/A)      	SUPPLY VOLTAGE: +/-15V
* CONNECTIONS:   	NON-INVERTING INPUT
*                	| INVERTING INPUT
*                	| | POSITIVE POWER SUPPLY
*                	| | | NEGATIVE POWER SUPPLY
*                	| | | | OUTPUT
*                	| | | | |
.SUBCKT LM348T_plus   	1 2 3 4 5
*
*Built-in voltage source
.PARAM rinmeas_dc = 0
.PARAM inres = 1E12
Vinres 222 2 DC=rinmeas_dc
Rinres 1 222 inres

C1   11 12 9.461E-12
C2    6  7 30.00E-12
DC    5 53 DX
DE   54  5 DX
DLP  90 91 DX
DLN  92 90 DX
DP    4  3 DX
EGND 99  0 POLY(2) 3 0 4 0 0 .5 .5
FB    7 99 POLY(5) VB VC VE VLP VLN 0 4.715E6 -5E6 5E6 5E6 -5E6
GA 6  0 11 12 256.2E-6
GCM   0  6 10 99 4.023E-9
IEE  10  4 DC 15.06E-6
HLIM 90  0 VLIM 1K
Q1   11  2 13 QX
Q2   12  1 14 QX
R2    6  9 100.0E3
RC1   3 11 4.420E3
RC2   3 12 4.420E3
RE1  13 10 968
RE2  14 10 968
REE  10 99 13.28E6
RO1   8  5 150
RO2   7 99 150
RP    3  4 51.28E3
VB    9  0 DC 0
VC    3 53 DC 3.600
VE   54  4 DC 3.600
VLIM  7  8 DC 0
VLP  91  0 DC 25
VLN   0 92 DC 25
.MODEL DX D(IS=800.0E-18)
.MODEL QX NPN(IS=800.0E-18 BF=250)
.ENDS
**********










.MODEL T2N2222 NPN (IS=15.2F NF=1 BF=105 VAF=98.5 IKF=.5
+ ISE=8.2P NE=2 BR=4 NR=1 VAR=20 IKR=.225 RE=.373 RB=1.49
+ RC=.149 XTB=1.5 CJE=35.5P CJC=12.2P TF=500P TR=85N)

.model 2n2907 pnp (is=1.1P bf=200 nf=1.2 vaf=50 ikf=0.1 ise=13p
+ ne=1.9 br=6 rc=0.6 cje=23p vje=0.85 mje=1.25 tf=0.5n cjc=19p
+ vjc=0.5 mjc=0.2 tr=34n xtb=1.5)
.end 



.MODEL BC238B NPN IS=1.8E-14 ISE=5.0E-14 NF=.9955 NE=1.46 BF=400 BR=35.5
+ IKF=.14 IKR=.03 ISC=1.72E-13 NC=1.27 NR=1.005 RB=.56 RE=.6 RC=.25 VAF=80
+ VAR=12.5 CJE=13E-12 TF=.64E-9 CJC=4E-12 TR=50.72E-9 VJC=.54 MJC=.33

.MODEL BC308B PNP IS=1.149E-14 ISE=5E-14 NF=.9872 NE=1.4 BF=330 BR=13
+ IKF=.1 IKR=.012 ISC=1.43E-14 NC=1.1 NR=.996 RB=.2 RE=.4 RC=.95 VAR=8.15
+ VAF=84.56 CJE=16E-12 TF=.493E-9 CJC=10.5E-12 TR=73.55E-9 MJC=.415 VJC=.565


.model 2n3510 npn
  + bf=100 br=1.35e-4 xtb=1.5 is=8.35e-14 eg=1.11 cjc=9.63e-12
  + cje=9.47e-12 rb=16.7 rc=1.66 vaf=90 tf=1e-10 tr=1.27e-4
  + cjs=1e-15 vjs=0.8 mjs=0.5 var=100 ise=4.77e-11 isc=1e-16
  + ikf=0.18 ikr=1000 irb=1 rbm=0 vtf=1000

.end


.SUBCKT par3pnp 1 2 3
q_1 1 2 3 BC308B
q_2 1 2 3 BC308B
q_3 1 2 3 BC308B
.ENDS

* 4.7V, 5mA, 300pF, 70Ohm
.model zd4v7 d is=1e-10 rs=70 bv=4.3085 cjo=300p

* 1N4148
.MODEL D1N4148 D 
+ IS = 4.352E-9 
+ N = 1.906 
+ BV = 110 
+ IBV = 0.0001 
+ RS = 0.6458 
+ CJO = 7.048E-13 
+ VJ = 0.869 
+ M = 0.03 
+ FC = 0.5 
+ TT = 3.48E-9 

********
* LM358 OPERATIONAL AMPLIFIER
* "MACROMODEL" SUBCIRCUIT
*
* CONNECTIONS:
*	1 -	NON-INVERTING INPUT
*	2 -	INVERTING INPUT
*	3 -	POSITIVE POWER SUPPLY
*	4 -	NEGATIVE POWER SUPPLY
*	5 -	OUTPUT
*
.SUBCKT LM358 1 2 3 4 5
*
C1   11 12 3.000E-12
C2    6  7 6.000E-12
CEE  10 99 388.6E-15
DC    5 53 DX
DE   54  5 DX
DLP  90 91 DX
DLN  92 90 DX
DP    4  3 DX
EGND 99  0 POLY(2) (3,0) (4,0) 0 .5 .5
FB    7 99 POLY(5) VB VC VE VLP VLN 0 48.23E6 
+ -50E6 50E6 50E6 -50E6
GA    6  0 11 12 41.47E-6
GCM   0  6 10 99 414.7E-12
IEE   3 10 DC 2.543E-6
HLIM 90  0 VLIM 1K
Q1   11  2 13 QX
Q2   12  1 14 QX
R2    6  9 100.0E3
RC1   4 11 24.11E3
RC2   4 12 24.11E3
RE1  13 10 3.043E3
RE2  14 10 3.043E3
REE  10 99 78.65E6
RO1   8  5 50
RO2   7 99 25
RP    3  4 37.62E3
VB    9  0 DC 0
VC    3 53 DC 1.800
VE   54  4 DC .88
VLIM  7  8 DC 0
VLP  91  0 DC 40
VLN   0 92 DC 40
DESD1  1  3 DX
DESD2  2  3 DX
.MODEL DX D(IS=800.0E-18)
.MODEL QX PNP(IS=800.0E-18 BF=32.03)
.ENDS
***********

**********
* LM348 Opamp "MACROMODEL" SUBCIRCUIT
* (REV N/A)      SUPPLY VOLTAGE: +/-15V
* CONNECTIONS:   NON-INVERTING INPUT
*                | INVERTING INPUT
*                | | POSITIVE POWER SUPPLY
*                | | | NEGATIVE POWER SUPPLY
*                | | | | OUTPUT
*                | | | | |
.SUBCKT LM348T   1 2 3 4 5
*
C1   11 12 9.461E-12
C2    6  7 30.00E-12
DC    5 53 DX
DE   54  5 DX
DLP  90 91 DX
DLN  92 90 DX
DP    4  3 DX
EGND 99  0 POLY(2) 3 0 4 0 0 .5 .5
FB    7 99 POLY(5) VB VC VE VLP VLN 0 4.715E6 -5E6 5E6 5E6 -5E6
GA 6  0 11 12 256.2E-6
GCM   0  6 10 99 4.023E-9
IEE  10  4 DC 15.06E-6
HLIM 90  0 VLIM 1K
Q1   11  2 13 QX
Q2   12  1 14 QX
R2    6  9 100.0E3
RC1   3 11 4.420E3
RC2   3 12 4.420E3
RE1  13 10 968
RE2  14 10 968
REE  10 99 13.28E6
RO1   8  5 150
RO2   7 99 150
RP    3  4 51.28E3
VB    9  0 DC 0
VC    3 53 DC 3.600
VE   54  4 DC 3.600
VLIM  7  8 DC 0
VLP  91  0 DC 25
VLN   0 92 DC 25
.MODEL DX D(IS=800.0E-18)
.MODEL QX NPN(IS=800.0E-18 BF=250)
.ENDS
**********

**********
* UA741 Opamp "MACROMODEL" SUBCIRCUIT
* (REV N/A)      SUPPLY VOLTAGE: +/-15V
* CONNECTIONS:   NON-INVERTING INPUT
*                | INVERTING INPUT
*                | | POSITIVE POWER SUPPLY
*                | | | NEGATIVE POWER SUPPLY
*                | | | | OUTPUT
*                | | | | |
.SUBCKT UA741T   1 2 3 4 5
*
C1   11 12 4.664E-12
C2    6  7 20.00E-12
DC    5 53 DX
DE   54  5 DX
DLP  90 91 DX
DLN  92 90 DX
DP    4  3 DX
EGND 99  0 POLY(2) 3 0 4 0 0 .5 .5
FB    7 99 POLY(5) VB VC VE VLP VLN 0 10.61E6 -10E6 10E6 10E6 -10E6
GA 6  0 11 12 137.7E-6
GCM 0  6 10 99 2.574E-9
IEE  10  4 DC 10.16E-6
HLIM 90  0 VLIM 1K
Q1   11  2 13 QX
Q2   12  1 14 QX
R2    6  9 100.0E3
RC1   3 11 7.957E3
RC2   3 12 7.957E3
RE1  13 10 2.740E3
RE2  14 10 2.740E3
REE  10 99 19.69E6
RO1   8  5 150
RO2   7 99 150
RP    3  4 18.11E3
VB    9  0 DC 0
VC 3 53 DC 2.600
VE   54  4 DC 2.600
VLIM  7  8 DC 0
VLP  91  0 DC 25
VLN   0 92 DC 25
.MODEL DX D(IS=800.0E-18)
.MODEL QX NPN(IS=800.0E-18 BF=62.50)
.ENDS
**********

**********
*LM741/A/C/E OPERATIONAL AMPLIFIER MACRO-MODEL
*
* connections:      non-inverting input
*                   |   inverting input
*                   |   |   positive power supply
*                   |   |   |   negative power supply
*                   |   |   |   |   output
*                   |   |   |   |   |
*                   |   |   |   |   |
.SUBCKT LM741N      1   2  99  50  28
*
*Features:
*Improved performance over industry standards
*Plug-in replacement for LM709,LM201,MC1439,748
*Input and output overload protection
*
*INPUT STAGE
*
IOS 2 1 20N
*^Input offset current
R1 1 3 250K
R2 3 2 250K
I1 4 50 100U
R3 5 99 517
R4 6 99 517
Q1 5 2 4 QX
Q2 6 7 4 QX
*Fp2=2.55 MHz
C4 5 6 60.3614P
*
*COMMON MODE EFFECT
*
I2 99 50 1.6MA
*^Quiescent supply current
EOS 7 1 POLY(1) 16 49 1E-3 1
*Input offset voltage.^
R8 99 49 40K
R9 49 50 40K
*
*OUTPUT VOLTAGE LIMITING
V2 99 8 1.63
D1 9 8 DX
D2 10 9 DX
V3 10 50 1.63
*
*SECOND STAGE
*
EH 99 98 99 49 1
G1 98 9 5 6 2.1E-3
*Fp1=5 Hz
R5 98 9 95.493MEG
C3 98 9 333.33P
*
*POLE STAGE
*
*Fp=30 MHz
G3 98 15 9 49 1E-6
R12 98 15 1MEG
C5 98 15 5.3052E-15
*
*COMMON-MODE ZERO STAGE
*
*Fpcm=300 Hz
G4 98 16 3 49 3.1623E-8
L2 98 17 530.5M
R13 17 16 1K
*
*OUTPUT STAGE
*
F6 50 99 POLY(1) V6 450U 1
E1 99 23 99 15 1
R16 24 23 25
D5 26 24 DX
V6 26 22 0.65V
R17 23 25 25
D6 25 27 DX
V7 22 27 0.65V
V5 22 21 0.18V
D4 21 15 DX
V4 20 22 0.18V
D3 15 20 DX
L3 22 28 100P
RL3 22 28 100K
*
*MODELS USED
*
.MODEL DX D(IS=1E-15)
.MODEL QX NPN(BF=625)
*
.ENDS
**********
