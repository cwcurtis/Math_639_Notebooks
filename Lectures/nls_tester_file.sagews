︠e55f61df-23de-4c58-88e3-3250f39c67d9︠
def Ht(f,y,ep,s):
    Htf = f.coefficient(y,0)
    for jj in range(1,4):
        Htf += I*s*f.coefficient(y,jj)*y^jj
        Htf += -I*s*f.coefficient(y,-jj)/y^(jj)
    return Htf

def Dx(f,y,k0):
    Dxf = 0
    for jj in range(1,4):
        Dxf += jj*k0*I*f.coefficient(y,jj)*y^jj
        Dxf += -jj*k0*I*f.coefficient(y,-jj)/y^(jj)
    return Dxf

def Dt(f,y,Om):
    Dtf = 0
    for jj in range(1,4):
        Dtf += jj*Om*I*f.coefficient(y,jj)*y^jj
        Dtf += -jj*Om*I*f.coefficient(y,-jj)/y^(jj)
    return Dtf

# y is a place holder for e^{i*k0*x+i*Om*t}
var('ep k0 s sig Om n0 n1 n1X n1XX n1T n2 n2X n2T q0 q1 q1X q1T q2 q2X q2T y a1 a2')

eta = ep^2*n0 + n1*y + n1/y + ep*(n2*y^2 + n2/y^2)
Q = ep^2*q0 + q1*y + q1/y + ep*(q2*y^2 + q2/y^2)

etax = Dx(eta,y,k0)
etaxx = Dx(etax.expand(),y,k0)
Qx = Dx(Q,y,k0)

etat = Dt(eta,y,Om)
etatx = Dx(etat.expand(),y,k0)

Qt = Dt(Q,y,Om)

etaX = ep*(n1X*y + n1X/y) + ep^2*(n2X*y^2 + n2X/y^2)
etaXX = ep^2*(n1XX*y + n1XX/y)
QX = ep*(q1X*y + q1X/y) + ep^2*(q2X*y^2 + q2X/y^2)

etaT = ep*(n1T*y + n1T/y) + ep^2*(n2T*y^2 + n2T/y^2)
QT = ep*(q1T*y + q1T/y) + ep^2*(q2T*y^2 + q2T/y^2)

etatX = Dt(etaX,y,Om)
etatX = etatX.expand()

etatT = Dt(etaT,y,Om)
etatT = etatT.expand()

etaxT = Dx(etaT,y,k0)
etaxT = etaxT.expand()

eta2 = (eta^2).expand()

Heta2 = Ht(eta2,y,ep,s)
Heta2 = Heta2.expand()

dteta2 = Dt(eta2,y,Om)
dteta2 = dteta2.expand()

eta_etat = (eta*etat).expand()
eta_etaT = (eta*etaT).expand()
eta_etaX = (eta*etaX).expand()
etat_etaT = (etat*etaT).expand()
etat_etatX = (etat*etatX).expand()
etaxt_etaT = (etatx*etaT).expand()
etat_etaxT = (etat*etaxT).expand()

eta3 = (n1*y + n1/y)^3
eta3 = eta3.expand()
dxeta3 = Dx(eta3,y,k0)
dxeta3 = dxeta3.expand()

Q2 = (Q^2).expand()

etaQ = (eta*Q).expand()
HetaQ = Ht(etaQ,y,ep,s)
HetaQ = HetaQ.expand()

eta2Q = (eta2*Q).expand()

pretot1 = ep/2*(-etat^2 + Q2) + ep^2*(-etat*etax*Q + 3/2*sig*etax^2*etaxx)
pretot2 = ep*(-etat_etatX - etaxt_etaT/2 -etat_etaxT/2 + Q*QX)
bnl1 = Dx(pretot1.expand(),y,k0) + pretot2

eq1 = Qt + QT - sig*( Dx(Dx(Dx(eta,y,k0),y,k0),y,k0) + 3*Dx(Dx(etaX,y,k0),y,k0) + 3*Dx(etaXX,y,k0) ) + etax + etaX + bnl1

pretot3 = etaX*etat + eta*etatX + etax*etaT + eta*etaxT
pretot4 = eta2*etat/2 + Ht(eta2Q/2,y,ep,s)

eq2 = etat + etaT + Ht(Q,y,ep,s) + ep*( Dx(etaQ - Ht(dteta2/2,y,ep,s),y,k0) + QX*eta + Q*etaX - Ht(pretot3.expand(),y,ep,s) ) - ep^2*Dx(Dx(pretot4.expand(),y,k0),y,k0)

fstfac1 = (eq1.expand().coefficient(ep,1).collect(y).coefficient(y,2)/I).subs({q1:-s*Om*n1}).collect(n1)

print fstfac1.coefficient(n1,0).collect(n2)
print fstfac1.coefficient(n1,2).subs({s^2:1})

print (eq1.expand().coefficient(ep,1).collect(y).coefficient(y,2)/I).subs({q1:-s*Om*n1}).collect(n2)
print eq1.expand().coefficient(ep,1).collect(y).coefficient(y,1)
print eq1.expand().coefficient(ep,1).collect(y).coefficient(y,0)
print
print (eq2.expand().coefficient(ep,1).collect(y).coefficient(y,2)/I).subs({q1:-s*Om*n1})
print eq2.expand().coefficient(ep,1).collect(y).coefficient(y,1)
print eq2.expand().coefficient(ep,1).collect(y).coefficient(y,0)

alpha1 = fstfac1.coefficient(n1,0).collect(n2).coefficient(n2,1)
alpha2 = -fstfac1.coefficient(n1,2).subs({s^2:1})

︡cef0b52e-9244-4e35-8a2d-fe024bcc9d0a︡{"stdout":"(ep, k0, s, sig, Om, n0, n1, n1X, n1XX, n1T, n2, n2X, n2T, q0, q1, q1X, q1T, q2, q2X, q2T, y, a1, a2)\n"}︡{"stdout":"(8*k0^3*sig + 2*k0)*n2 + 2*Om*q2\n"}︡{"stdout":"2*Om^2*k0\n"}︡{"stdout":"Om^2*k0*n1^2*s^2 + Om^2*k0*n1^2 + (8*k0^3*sig + 2*k0)*n2 + 2*Om*q2\n"}︡{"stdout":"3*k0^2*n1X*sig + n1X + q1T\n"}︡{"stdout":"0\n"}︡{"stdout":"\n"}︡{"stdout":"2*Om*n2 + q2*s\n"}︡{"stdout":"n1T\n"}︡{"stdout":"0\n"}︡{"done":true}︡
︠f04424a4-ed0f-4806-92f4-9aca2edbc360s︠

fincoeff1 = (((( eq1.expand().coefficient(ep,2).collect(y).coefficient(y,1).subs({q1:-s*Om*n1})).subs({n2:a1*n1^2})).subs({q2:a2*n1^2})).subs({s^2:1})).subs({s^3:s}).collect(n1^3)

fincoeff2 = (((( eq2.expand().coefficient(ep,2).collect(y).coefficient(y,1).subs({q1:-s*Om*n1})).subs({n2:a1*n1^2})).subs({q2:a2*n1^2})).subs({s^2:1})).subs({s^3:s}).collect(n1^3)

print eq1.expand().coefficient(ep,2).collect(y).coefficient(y,3)
print eq1.expand().coefficient(ep,2).collect(y).coefficient(y,2)
print fincoeff1
print eq1.expand().coefficient(ep,2).collect(y).coefficient(y,0)
print
print eq2.expand().coefficient(ep,2).collect(y).coefficient(y,3)
print eq2.expand().coefficient(ep,2).collect(y).coefficient(y,2)
print fincoeff2
print eq2.expand().coefficient(ep,2).collect(y).coefficient(y,0)


︡86da2520-5186-4dd3-b2a6-e82f48c062a0︡{"stdout":"9/2*I*k0^5*n1^3*sig + 3*I*Om*k0^2*n1^2*q1 + 6*I*Om^2*k0*n1*n2 + 3*I*k0*q1*q2\n"}︡{"stdout":"Om*k0*n1*n1T + Om^2*n1*n1X + 12*k0^2*n2X*sig + q1*q1X + n2X + q2T\n"}︡{"stdout":"-1/2*(3*I*k0^5*sig - 2*I*Om^2*k0^2*s + 4*I*Om^2*a1*k0 + 2*I*Om*a2*k0*s)*n1^3 - 3*I*k0*n1XX*sig\n"}︡{"stdout":"-2*Om^2*n1*n1X + 2*q1*q1X\n"}︡{"stdout":"\n"}︡{"stdout":"9/2*I*Om*k0^2*n1^3 + 9/2*I*k0^2*n1^2*q1*s + 9*I*Om*k0*n1*n2*s + 3*I*k0*n2*q1 + 3*I*k0*n1*q2\n"}︡{"stdout":"2*k0*n1*n1T*s + 2*Om*n1*n1X*s + n1X*q1 + n1*q1X + n2T\n"}︡{"stdout":"(-I*Om*k0^2 + I*a2*k0)*n1^3\n"}︡{"stdout":"2*n1X*q1 + 2*n1*q1X + q0\n"}︡{"done":true}︡
︠b9df322d-6785-47d6-8334-af3c7bdb680a︠
a1 = (s*alpha2)/(s*alpha1-4*Om^2)
a2 = (-2*Om*alpha2)/(s*alpha1-4*Om^2)
alpha3 = (fincoeff1.coefficient(n1,3)/I-s*Om*(fincoeff2.coefficient(n1,3)/I)).expand().collect(Om).subs({s^2:1})

print a1
print a2
print alpha3
︡0dd16f2d-9cb5-4da7-92c9-e4d6a3f93cb8︡{"stdout":"(2*Om*k0*w + (2*Om*k0*s*w + 2*Om^2*k0 + k0*w^2)*s)/(4*Om^2 - (8*k0^3*sig - 2*Om*w + 2*k0)*s)\n"}︡{"stdout":"-((8*k0^3*sig - 2*Om*w + 2*k0)*k0*w + 2*(2*Om*k0*s*w + 2*Om^2*k0 + k0*w^2)*Om)/(4*Om^2 - (8*k0^3*sig - 2*Om*w + 2*k0)*s)\n"}︡{"stdout":"-3/2*k0^5*sig + a1*k0*w^2 + (2*k0^2*s - 2*a1*k0)*Om^2 - a2*k0*w + (a1*k0*s*w - 2*a2*k0*s + 2*k0^2*w)*Om\n"}︡{"done":true}︡
︠fe5c20df-97f8-40dd-b1d3-5b55b39c3b4c︠
︡05b513fe-2cc5-45fe-9767-afdf015d7273︡{"done":true}︡









