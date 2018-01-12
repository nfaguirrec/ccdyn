De = 5.0
Re = 5.0
alpha = 1.0
c6 = 7.0
c8 = 7.0
R0 = 7.0
w = 3.0

damp(x,R0,w) = 0.5*(1.0-tanh(6.0*(x-R0)/w))

MorseDispDamp( r, De, alpha, Re, R0, w, c6, c8 ) = damp(r,R0,w)*Morse( r, De, alpha, Re ) - (1.0-damp(r,R0,w))*( sgn(c6)*(c6/r)**6 + sgn(c8)*(c8/r)**8 )

De = 0.0
fit [6.0:] MorseDispDamp( x, De, alpha, Re, R0, w, c6, c8 ) "MgHePotential.dat" u 1:($2+$3) via c6,c8
De = 5.0
fit [4.2:6.0] Morse( x, De, alpha, Re ) "MgHePotential.dat" u 1:($2+$3) via De,alpha,Re
fit [4.2:] MorseDispDamp( x, De, alpha, Re, R0, w, c6, c8 ) "MgHePotential.dat" u 1:($2+$3) via De, alpha, Re, c6, c8
fit [4.2:] MorseDispDamp( x, De, alpha, Re, R0, w, c6, c8 ) "MgHePotential.dat" u 1:($2+$3) via De, alpha, Re, R0, w, c6, c8

# fit [4.2:] MorseDisp( x, De, alpha, Re, c6, c8 ) "MgHePotential.dat" u 1:($2+$3) via De,alpha,Re,c6,c8
# 
plot [] [-6:10] \
MorseDispDamp( x, De, alpha, Re, R0, w, c6, c8 ), \
"MgHePotential.dat" u 1:($2+$3) w p pt 7, \
damp(x,R0,w)

pause -1