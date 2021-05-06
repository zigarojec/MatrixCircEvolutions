 * This is a list of models that can be used for sythesis of resilent circuit topologies.
 
.subckt T2N2222_resil_nom C B E
q1 (C B E) T2N2222
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

