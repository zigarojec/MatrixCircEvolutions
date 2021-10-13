 * This is a list of models that can be used for sythesis of resilent circuit topologies.
 * xint (internal) nodes are created for measurement purposes (avoiding global node naming of outer-terminals in Spice!)
 

 * NPN *************************************
.subckt T2N2222_resil_nom c b e
vc c cint dc=0
vb b bint dc=0
ve e eint dc=0
q1 (cint bint eint) T2N2222
.ends
 
.subckt T2N2222_resil_sck C B E
* Failure anti-model
* Short-circuit over collector-emittor
q1 (C B E) T2N2222
r1 (C E)   r=1
.ends
 
.subckt T2N2222_resil_himp C B E
* Failure anti-model
* High impedance triangle (device removed)
rc  (C mid) r=1meg
rb  (B mid) r=1meg
re  (E mid) r=1meg
*q1 (C B E) T2N2222
.ends
 
 
 .MODEL T2N2222 NPN (IS=15.2F NF=1 BF=105 VAF=98.5 IKF=.5
+ ISE=8.2P NE=2 BR=4 NR=1 VAR=20 IKR=.225 RE=.373 RB=1.49
+ RC=.149 XTB=1.5 CJE=35.5P CJC=12.2P TF=500P TR=85N)
* END NPN *********************************

* PNP *************************************
.subckt 2N2907_resil_nom c b e
vc c cint dc=0
vb b bint dc=0
ve e eint dc=0
q1 (cint bint eint) 2N2907
.ends
 
.subckt 2N2907_resil_sck C B E
* Failure anti-model
* Short-circuit over collector-emittor
q1 (C B E) 2N2907
r1 (C E)   r=1
.ends
 
.subckt 2N2907_resil_himp C B E
* Failure anti-model
* High impedance triangle (device removed)
rc  (C mid) r=1meg
rb  (B mid) r=1meg
re  (E mid) r=1meg
*q1 (C B E) 2N2907
.ends
* END PNP *********************************


* PNP 2N2907 (complementary transistor to T2N2222)
.model 2N2907 pnp (is=1.1P bf=200 nf=1.2 vaf=50 ikf=0.1 ise=13p
+ ne=1.9 br=6 rc=0.6 cje=23p vje=0.85 mje=1.25 tf=0.5n cjc=19p
+ vjc=0.5 mjc=0.2 tr=34n xtb=1.5)
*.end 

* DIODES (zd4v7) **************************
* 4.7V, 5mA, 300pF, 70Ohm
* BARE MODEL (ZENER)
.model zd4v7 d is=1e-10 rs=70 bv=4.3085 cjo=300p

* RESILENCE MODEL
.subckt zd4v7_resil_nom p n 
vp p pint dc=0
vn n nint dc=0

d pint nint zd4v7
.ends

.subckt zd4v7_resil_himp p n 
vp p pint dc=0
vn n nint dc=0

*d1 pint nint zd4v7
r1 pint nint r=1meg
.ends


.subckt zd4v7_resil_sck p n 
vp p pint dc=0
vn n nint dc=0

d1 pint nint zd4v7
r1 pint nint r=1m
.ends
* END DIODES (zd4v7) ***********************


* RESISTOR RESILIENCE MODEL ***************
.subckt resistor_nom p n param:r=10k
vp p pint dc=0
vn n nint dc=0

r1 pint nint r={r}
.ends

.subckt resistor_himp p n param:r=10k
vp p pint dc=0
vn n nint dc=0

r1 pint nint r=10meg
.ends

.subckt resistor_sck p n param:r=10k
vp p pint dc=0
vn n nint dc=0

r1 pint nint r=1m
.ends


* END RESISTOR ****************************

* VDC *************************************
.subckt vdc p n param: v=1.69
vdc p n dc={v}
.ends

* END VDC**********************************