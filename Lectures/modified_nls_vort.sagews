︠17039b5a-5b08-4ccc-b3d6-1b7ad138c220s︠
#reset()
%auto
typeset_mode(True)

def Ht(f,y,xi,tau,s,Hv):
    Htf(xi,tau) = function('Htf')(xi,tau)
    Htf(xi,tau) = 0
    f = f.expand().collect(y)
    Htf(xi,tau) = Hv*f.coefficient(y,0).collect(ep)

    for jj in range(1,5):
        Htf(xi,tau) = Htf(xi,tau) + I*s*f.coefficient(y,jj)*y^jj
        Htf(xi,tau) = Htf(xi,tau) - I*s*f.coefficient(y,-jj)/y^jj
    return Htf.expand()

def Dx(f,y,xi,tau,k0):
    Dxf = function('Dxf')(xi,tau)
    Dxf(xi,tau) = 0
    f = f.expand().collect(y)
    for jj in range(1,5):
        Dxf(xi,tau) = Dxf(xi,tau) + jj*k0*I*(f.coefficient(y,jj))*y^jj
        Dxf(xi,tau) = Dxf(xi,tau) - jj*k0*I*(f.coefficient(y,-jj))/y^jj
    return Dxf.expand()

def Dt(f,y,xi,tau,Om):
    Dtf = function('Dtf')(xi,tau)
    Dtf(xi,tau) = 0
    f = f.expand().collect(y)
    for jj in range(1,5):
        Dtf(xi,tau) = Dtf(xi,tau) + jj*Om*I*(f.coefficient(y,jj))*y^jj
        Dtf(xi,tau) = Dtf(xi,tau) - jj*Om*I*(f.coefficient(y,-jj))/y^jj
    return Dtf.expand()
# y is a place holder for e^{i*k0*x+i*Om*t}
# b denotes the Bond number
var('ep epsilon k0 Omega c_g x t y xi tau s H_s w sigma a0 a1 a_d')
n0(xi,tau) = function('n0')(xi,tau)

Hn0(xi,tau) = function('Hn0')(xi,tau)
Hn2(xi,tau) = function('Hn2')(xi,tau) # H(n1*conj(n1))
Hn2v2(xi,tau) = function('Hn2v2')(xi,tau) # H(n2*conj(n2))
Hn02(xi,tau) = function('Hn02')(xi,tau) # H(n0^2)

n1(xi,tau) = function('n1')(xi,tau)
v1(xi,tau) = function('v1')(xi,tau)

n2(xi,tau) = function('n2')(xi,tau)
v2(xi,tau) = function('v2')(xi,tau)

n3(xi,tau) = function('n3')(xi,tau)
v3(xi,tau) = function('v3')(xi,tau)

n4(xi,tau) = function('n4')(xi,tau)
v4(xi,tau) = function('v4')(xi,tau)

n(xi,tau) = ep*n0(xi,tau) + n1(xi,tau)*y + v1(xi,tau)/y + ep*( n2(xi,tau)*y^2 + v2(xi,tau)/y^2 ) + ep^2*( n3(xi,tau)*y^3 + v3(xi,tau)/y^3 ) + ep^3*( n4(xi,tau)*y^4 + v4(xi,tau)/y^4 )
nzav(xi,tau) = n(xi,tau) - ep*n0(xi,tau)
Hn(xi,tau) = ep*Hn0(xi,tau) + Ht(nzav,y,xi,tau,s,H_s)

np2(xi,tau) = (n(xi,tau)*n(xi,tau)).expand().collect(ep).series(ep,3).truncate()
np3(xi,tau) = (n(xi,tau)*np2(xi,tau)).expand().collect(ep).series(ep,2).truncate()
np4(xi,tau) = (n(xi,tau)*np3(xi,tau)).expand().collect(ep).series(ep,1).truncate()

nt(xi,tau) = (Dt(n,y,xi,tau,Omega) + ep*c_g*diff(n(xi,tau),xi) + ep^2*diff(n(xi,tau),tau)).expand().collect(ep).series(ep,4).truncate()
nx(xi,tau) = (Dx(n,y,xi,tau,k0) + ep*diff(n(xi,tau),xi)).expand().collect(ep).series(ep,4).truncate()
nxx(xi,tau) = (Dx(nx,y,xi,tau,k0) + ep*diff(nx(xi,tau),xi)).expand().collect(ep).series(ep,4).truncate()

Htnt(xi,tau) = (Dt(Hn,y,xi,tau,Omega) + ep*c_g*diff(Hn(xi,tau),xi) + ep^2*diff(Hn(xi,tau),tau)).expand().collect(ep).series(ep,4).truncate()

zavnp2(xi,tau) = (np2(xi,tau) - np2(xi,tau).collect(y).coefficient(y,0)).expand()
Htnp2(xi,tau) = ( 2*Hn2(xi,tau) + ep^2*(2*Hn2v2+Hn02) + Ht(zavnp2,y,xi,tau,s,H_s) ).expand().collect(ep).series(ep,3).truncate()

nlsub1(xi,tau) = (n(xi,tau)*Htnt(xi,tau) + w*np2(xi,tau)/2 - (Dt(Htnp2,y,xi,tau,Omega) + ep*c_g*diff(Htnp2(xi,tau),xi) + ep^2*diff(Htnp2(xi,tau),tau))/2).expand().collect(ep).series(ep,3).truncate()
R1(xi,tau) = (Dx(nlsub1,y,xi,tau,k0) + ep*diff(nlsub1(xi,tau),xi)).expand().collect(ep).series(ep,3).truncate()
HtR1(xi,tau) = (Ht(R1,y,xi,tau,s,H_s)).expand().collect(ep).series(ep,3).truncate().subs({s^2:1}).subs({s^3:s}).subs({s^4:1})

HtR1za(xi,tau) = (HtR1(xi,tau) - HtR1(xi,tau).collect(y).coefficient(y,0)).expand()
Htr1zm(xi,tau) = -ep*(2*Omega*s-w)*diff(Hn2(xi,tau),xi) + ep^2*c_g*diff(diff(Hn2(xi,tau),xi),xi) + ep^2*c_g*s/a_d*diff(Hn2(xi,tau),tau)
HtR1(xi,tau) = (Htr1zm(xi,tau) + HtR1za(xi,tau))

nltm1(xi,tau) = (np2(xi,tau)*Htnt(xi,tau)/2 + w*np3(xi,tau)/3).expand().collect(ep).series(ep,2).truncate()
nlsub2(xi,tau) = (Dt(np3,y,xi,tau,Omega)/6 + ep*c_g*diff(np3(xi,tau),xi)/6 + Ht(nltm1,y,xi,tau,s,H_s)).expand().collect(ep).series(ep,2).truncate()
nltm2(xi,tau) = (Dx(nlsub2,y,xi,tau,k0)+ep*diff(nlsub2(xi,tau),xi)).expand().collect(ep).series(ep,2).truncate()
nlsub3(xi,tau) = (HtR1(xi,tau)*n(xi,tau) - nltm2(xi,tau)).expand().collect(ep).series(ep,2).truncate()
R2(xi,tau) = Dx(nlsub3,y,xi,tau,k0) + ep*diff(nlsub3(xi,tau),xi)
HtR2(xi,tau) = (Ht(R2,y,xi,tau,s,H_s)).expand().collect(ep).series(ep,2).truncate().subs({s^2:1}).subs({s^3:s}).subs({s^4:1})

np2R1(xi,tau) = (np2(xi,tau)*HtR1(xi,tau)).expand().collect(ep).series(ep,1).truncate()
nlsub4(xi,tau) = Dt(Ht(np4,y,xi,tau,s,H_s),y,xi,tau,Omega)/24 - np3(xi,tau)*Htnt(xi,tau)/6 - w*np4(xi,tau)/8
nltm3(xi,tau) = Dx(nlsub4,y,xi,tau,k0).expand().collect(ep).series(ep,1).truncate()
nlsub5(xi,tau) = -Ht(np2R1,y,xi,tau,s,H_s)/2 + nltm3(xi,tau)
R3(xi,tau) = (n(xi,tau)*HtR2(xi,tau) + Dx(nlsub5,y,xi,tau,k0)).expand().collect(ep).series(ep,1).truncate()
HtR3(xi,tau) = (Dx(Ht(R3,y,xi,tau,s,H_s),y,xi,tau,k0)).expand().subs({s^2:1}).subs({s^3:s}).subs({s^4:1})

Q(xi,tau) = (Htnt(xi,tau) + ep*HtR1(xi,tau) + ep^2*HtR2(xi,tau) + ep^3*HtR3(xi,tau)).expand()

onQ(xi,tau) = (w*n(xi,tau) + Q(xi,tau)).expand().collect(y)
onQ2(xi,tau) = ((onQ(xi,tau))^2-(nt(xi,tau))^2).expand().subs({s^2:1}).subs({s^3:s}).subs({s^4:1}).collect(ep).series(ep,3).truncate().collect(y)

nxsq(xi,tau) = (nx(xi,tau)*nx(xi,tau)).expand().series(ep,2).truncate()
ntnx(xi,tau) = (nt(xi,tau)*nx(xi,tau)).expand().series(ep,2).truncate()

nlfin(xi,tau) = (ep*onQ2(xi,tau)/2 + ep^2*(3/2*sigma*(nxx(xi,tau)*nxsq(xi,tau)).expand()-(ntnx(xi,tau)*onQ(xi,tau)).expand()) - ep^3*(nxsq(xi,tau)*onQ2(xi,tau)/2).expand()).expand().collect(ep).series(ep,4).truncate().collect(y)

dnlf(xi,tau) = (Dx(nlfin,y,xi,tau,k0) + ep*diff(nlfin(xi,tau),xi)).expand().collect(ep).series(ep,4).truncate()

bly(xi,tau) = (Dt(onQ,y,xi,tau,Omega).expand() + ep*c_g*diff(onQ(xi,tau),xi).expand() + ep^2*diff(onQ(xi,tau),tau).expand() + nx(xi,tau) - sigma*(Dx(nxx,y,xi,tau,k0)+ep*diff(nxx(xi,tau),xi)) + dnlf(xi,tau)).expand().collect(ep).series(ep,4).truncate().subs({s^2:1}).subs({s^3:s}).subs({s^4:1})

phix(xi,tau) = (Q(xi,tau) - ep*nx(xi,tau)*nt(xi,tau) - ep^2*nxsq(xi,tau)*onQ(xi,tau) ).expand().series(ep,3).truncate().collect(y)
phiz(xi,tau) = (nt(xi,tau) + ep*nx(xi,tau)*onQ(xi,tau) - ep^2*(ntnx(xi,tau))^2).expand().series(ep,3).truncate().collect(y)

dxphix(xi,tau) = (Dx(phix,y,xi,tau,k0) + ep*diff(phix(xi,tau),xi)).expand().series(ep,2).truncate()
dxphiz(xi,tau) = (Dx(phiz,y,xi,tau,k0) + ep*diff(phiz(xi,tau),xi)).expand().series(ep,2).truncate()
dx2phix(xi,tau) = Dx(dxphix,y,xi,tau,k0).expand().series(ep,1).truncate()
dx2phiz(xi,tau) = Dx(dxphiz,y,xi,tau,k0).expand().series(ep,1).truncate()

Hdxphix(xi,tau) = (Ht(dxphix,y,xi,tau,s,H_s)).expand().series(ep,2).truncate()
Hdxphiz(xi,tau) = (Ht(dxphiz,y,xi,tau,s,H_s)).expand().series(ep,2).truncate()
nHdxphix(xi,tau) = (n(xi,tau)*Hdxphix(xi,tau)).expand().series(ep,2).truncate()
nHdxphiz(xi,tau) = (n(xi,tau)*Hdxphiz(xi,tau)).expand().series(ep,2).truncate()
dxnHdxphix(xi,tau) = Dx(nHdxphix,y,xi,tau,k0).expand().series(ep,1).truncate()
dxnHdxphiz(xi,tau) = Dx(nHdxphiz,y,xi,tau,k0).expand().series(ep,1).truncate()

At(xi,tau) = (phix(xi,tau) + ep*nHdxphix(xi,tau) + ep^2*n(xi,tau)*( Ht(dxnHdxphix,y,xi,tau,s,H_s) + n(xi,tau)*dx2phix(xi,tau)/2 ) ).expand().series(ep,3).truncate().subs({s^3:s}).subs({s^2:1}).collect(y)
Bt(xi,tau) = (phiz(xi,tau) + ep*nHdxphiz(xi,tau) + ep^2*n(xi,tau)*( Ht(dxnHdxphiz,y,xi,tau,s,H_s) + n(xi,tau)*dx2phiz(xi,tau)/2 ) ).expand().series(ep,3).truncate().subs({s^3:s}).subs({s^2:1}).collect(y)
︡392d054d-e2c2-46c1-a9b3-8aa815f68484︡{"html":"<div align='center'>($\\displaystyle \\mathit{ep}$, $\\displaystyle \\epsilon$, $\\displaystyle k_{0}$, $\\displaystyle \\Omega$, $\\displaystyle c_{g}$, $\\displaystyle x$, $\\displaystyle t$, $\\displaystyle y$, $\\displaystyle \\xi$, $\\displaystyle \\tau$, $\\displaystyle s$, $\\displaystyle H_{s}$, $\\displaystyle w$, $\\displaystyle \\sigma$, $\\displaystyle a_{0}$, $\\displaystyle a_{1}$, $\\displaystyle a_{d}$)</div>"}︡{"done":true}
︠a4686cc7-2ac8-4a76-b64b-2c84a1d9bb8ds︠

bly0(xi,tau) =  bly.coefficient(y,0).collect(ep)
bly1(xi,tau) =  bly.coefficient(y,1).collect(ep)
bly2(xi,tau) =  bly.coefficient(y,2).collect(ep)
bly3(xi,tau) =  bly.coefficient(y,3).collect(ep)
#bly4(xi,tau) =  bly.coefficient(y,4).collect(ep)

fsorder(xi,tau) = bly1(xi,tau).coefficient(ep,0)
snorder(xi,tau) = bly1(xi,tau).coefficient(ep,1)
nls(xi,tau) = bly1(xi,tau).coefficient(ep,2)
nlsho(xi,tau) = bly1(xi,tau).coefficient(ep,3)

show(bly0(xi,tau).coefficient(ep,2))
print
show(bly0(xi,tau).coefficient(ep,3))
print

show(fsorder(xi,tau))
print
show(snorder(xi,tau))
print
show(nls(xi,tau))
print
show(nlsho(xi,tau))
print

show(bly2(xi,tau).coefficient(ep,1))
print
show(bly2(xi,tau).coefficient(ep,2))
print
#show(bly2(xi,tau).coefficient(ep,3))
#print
#show(nlfin(xi,tau).coefficient(y,3).collect(ep))
#print
#show((bly3(xi,tau).coefficient(ep,2)/I).collect(n3))
#print
#show(bly3(xi,tau).coefficient(ep,3))
#print
#print "phix coefficients"
#show(phix(xi,tau).coefficient(y,0).collect(ep))
#print
#show(phix(xi,tau).coefficient(y,1).collect(ep))
#print
#show(phix(xi,tau).coefficient(y,2).collect(ep))
#print
#show(phix(xi,tau).coefficient(y,3).collect(ep))


#print "At coefficients"
#show(At(xi,tau).coefficient(y,0))
#print
#show(At(xi,tau).coefficient(y,1).collect(ep))
#print
#show(At(xi,tau).coefficient(y,2).collect(ep))
#print
#show(At(xi,tau).coefficient(y,3).collect(ep))

#print "Bt coefficients"
#show(Bt(xi,tau).coefficient(y,0))
#print
#show(Bt(xi,tau).coefficient(y,1))
#print
#show(Bt(xi,tau).coefficient(y,2))
#print
#show(Bt(xi,tau).coefficient(y,3))
#print
#show(Bt(xi,tau).coefficient(y,4))
︡7e09b9b1-f64c-40ed-b8f7-56569a89e125︡{"html":"<div align='center'>$\\displaystyle -2 \\, \\Omega s w v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 2 \\, \\Omega s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) + w^{2} v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + w^{2} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) + c_{g} w \\frac{\\partial}{\\partial \\xi}n_{0}\\left(\\xi, \\tau\\right) + \\frac{\\partial}{\\partial \\xi}n_{0}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle i \\, c_{g} s w v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial^{2}}{(\\partial \\xi)^{2}}n_{1}\\left(\\xi, \\tau\\right) - i \\, c_{g} s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial^{2}}{(\\partial \\xi)^{2}}v_{1}\\left(\\xi, \\tau\\right) - 2 \\, \\Omega c_{g} s \\frac{\\partial^{2}}{(\\partial \\xi)^{2}}{\\rm Hn}_{2}\\left(\\xi, \\tau\\right) + c_{g}^{2} \\frac{\\partial^{2}}{(\\partial \\xi)^{2}}{\\rm Hn}_{0}\\left(\\xi, \\tau\\right) + c_{g} w \\frac{\\partial^{2}}{(\\partial \\xi)^{2}}{\\rm Hn}_{2}\\left(\\xi, \\tau\\right) + w \\frac{\\partial}{\\partial \\tau}n_{0}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle i \\, k_{0}^{3} \\sigma n_{1}\\left(\\xi, \\tau\\right) - i \\, \\Omega^{2} s n_{1}\\left(\\xi, \\tau\\right) + i \\, \\Omega w n_{1}\\left(\\xi, \\tau\\right) + i \\, k_{0} n_{1}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle -2 \\, \\Omega c_{g} s \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 3 \\, k_{0}^{2} \\sigma \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + c_{g} w \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{3}{2} i \\, k_{0}^{5} \\sigma n_{1}\\left(\\xi, \\tau\\right)^{2} v_{1}\\left(\\xi, \\tau\\right) + 2 i \\, \\Omega^{2} k_{0}^{2} s n_{1}\\left(\\xi, \\tau\\right)^{2} v_{1}\\left(\\xi, \\tau\\right) - i \\, k_{0}^{2} s w^{2} n_{1}\\left(\\xi, \\tau\\right)^{2} v_{1}\\left(\\xi, \\tau\\right) - 2 i \\, \\Omega k_{0} s w n_{0}\\left(\\xi, \\tau\\right) n_{1}\\left(\\xi, \\tau\\right) - 4 i \\, \\Omega k_{0} s w n_{2}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) + i \\, k_{0} w^{2} n_{0}\\left(\\xi, \\tau\\right) n_{1}\\left(\\xi, \\tau\\right) + 2 i \\, \\Omega^{2} k_{0} n_{2}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) + i \\, k_{0} w^{2} n_{2}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) + i \\, c_{g}^{2} s \\frac{\\partial^{2}}{(\\partial \\xi)^{2}}n_{1}\\left(\\xi, \\tau\\right) - 3 i \\, k_{0} \\sigma \\frac{\\partial^{2}}{(\\partial \\xi)^{2}}n_{1}\\left(\\xi, \\tau\\right) - 2 \\, \\Omega s \\frac{\\partial}{\\partial \\tau}n_{1}\\left(\\xi, \\tau\\right) + w \\frac{\\partial}{\\partial \\tau}n_{1}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle 4 \\, \\Omega c_{g} k_{0}^{2} s n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 9 \\, k_{0}^{4} \\sigma n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - \\frac{3}{2} \\, k_{0}^{4} \\sigma n_{1}\\left(\\xi, \\tau\\right)^{2} \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) + 6 \\, \\Omega^{2} k_{0} s n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 3 \\, k_{0} s w^{2} n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 2 \\, \\Omega^{2} k_{0} s n_{1}\\left(\\xi, \\tau\\right)^{2} \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) - k_{0} s w^{2} n_{1}\\left(\\xi, \\tau\\right)^{2} \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) - 2 i \\, \\Omega c_{g} k_{0} s n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{0}\\left(\\xi, \\tau\\right) - 4 i \\, \\Omega k_{0} s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{2}\\left(\\xi, \\tau\\right) - c_{g} k_{0} s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{0}\\left(\\xi, \\tau\\right) - 2 \\, c_{g} k_{0} s w n_{0}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 2 \\, c_{g} k_{0} s w v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{2}\\left(\\xi, \\tau\\right) + i \\, c_{g} k_{0} w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{0}\\left(\\xi, \\tau\\right) + 4 i \\, \\Omega^{2} k_{0} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{2}\\left(\\xi, \\tau\\right) + i \\, k_{0} w^{2} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{2}\\left(\\xi, \\tau\\right) - 2 \\, \\Omega s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{0}\\left(\\xi, \\tau\\right) - 2 \\, \\Omega s w n_{0}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 2 \\, \\Omega c_{g} k_{0} v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{2}\\left(\\xi, \\tau\\right) - 4 \\, \\Omega s w v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{2}\\left(\\xi, \\tau\\right) - 4 \\, \\Omega s w n_{2}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) + w^{2} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{0}\\left(\\xi, \\tau\\right) + w^{2} n_{0}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 2 \\, \\Omega^{2} v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{2}\\left(\\xi, \\tau\\right) + w^{2} v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{2}\\left(\\xi, \\tau\\right) + 2 \\, \\Omega^{2} n_{2}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) + w^{2} n_{2}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) + 2 i \\, c_{g} s \\frac{\\partial^{2}}{\\partial \\xi\\partial \\tau}n_{1}\\left(\\xi, \\tau\\right) - \\sigma \\frac{\\partial^{3}}{(\\partial \\xi)^{3}}n_{1}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle -4 i \\, \\Omega k_{0} s w n_{1}\\left(\\xi, \\tau\\right)^{2} + 2 i \\, \\Omega^{2} k_{0} n_{1}\\left(\\xi, \\tau\\right)^{2} + i \\, k_{0} w^{2} n_{1}\\left(\\xi, \\tau\\right)^{2} + 8 i \\, k_{0}^{3} \\sigma n_{2}\\left(\\xi, \\tau\\right) - 4 i \\, \\Omega^{2} s n_{2}\\left(\\xi, \\tau\\right) + 2 i \\, \\Omega w n_{2}\\left(\\xi, \\tau\\right) + 2 i \\, k_{0} n_{2}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle -4 \\, c_{g} k_{0} s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 4 \\, \\Omega c_{g} k_{0} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 4 \\, \\Omega s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 2 \\, \\Omega^{2} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + w^{2} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 4 \\, \\Omega c_{g} s \\frac{\\partial}{\\partial \\xi}n_{2}\\left(\\xi, \\tau\\right) + 12 \\, k_{0}^{2} \\sigma \\frac{\\partial}{\\partial \\xi}n_{2}\\left(\\xi, \\tau\\right) + c_{g} w \\frac{\\partial}{\\partial \\xi}n_{2}\\left(\\xi, \\tau\\right) + \\frac{\\partial}{\\partial \\xi}n_{2}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"done":true}︡
︠e1102ca0-e97c-4088-b73b-d0b7bb01e35as︠

︡4440dd0a-8cc4-4208-b5d7-e3e5dcd0dac2︡{"done":true}︡
︠e51e0758-ac9b-4044-8ed0-b6f643d08258s︠
l20 = 2*(k0 + w*Omega -2*s*Omega^2 + 4*sigma*k0^3)
a20 = w^2*k0 + 2*Omega^2*k0 - 4*Omega*s*w*k0
l21 = 1 + w*(c_g) - 4*s*(c_g)*Omega + 12*sigma*k0^2
a21 = (w^2)/2 + Omega^2 - 2*Omega*s*w + 2*(Omega-w*s)*(c_g)*k0
#n20(xi,tau) = -a20/l20*(n1(xi,tau))^2
#n21(xi,tau) =  2*I*(a21 - l21*a20/l20)/l20*n1(xi,tau)*diff(n1(xi,tau),xi)
n20(xi,tau) = a0*(n1(xi,tau))^2
n21(xi,tau) = I*a1*n1(xi,tau)*diff(n1(xi,tau),xi)

nlsex(xi,tau) = nls(xi,tau).subs({n2(xi,tau):(n20(xi,tau)+epsilon*n21(xi,tau))}).expand().subs({s^2:1})
nlsexho(xi,tau) = nlsho(xi,tau).subs({n2(xi,tau):n20(xi,tau)}).subs({diff(n2(xi,tau,xi),xi):2*a0*n1(xi,tau)*diff(n1(xi,tau),xi)}).expand().subs({s^2:1})
nlsfin(xi,tau) = (nlsex(xi,tau)+epsilon*nlsexho(xi,tau)).expand()
nlszero(xi,tau) = nlsfin(xi,tau).collect(epsilon).coefficient(epsilon,0)
nlsone(xi,tau) = nlsfin(xi,tau).collect(epsilon).coefficient(epsilon,1)

n12v1 = nlszero(xi,tau).coefficient((n1(xi,tau))^2*v1(xi,tau))
n1n0 = nlszero(xi,tau).coefficient(n1(xi,tau)*n0(xi,tau))

n1dn0 = nlsone(xi,tau).coefficient(n1(xi,tau)*diff(n0(xi,tau),xi))
dn1n0 = nlsone(xi,tau).coefficient(n0(xi,tau)*diff(n1(xi,tau),xi))
n1dHn0 = nlsone(xi,tau).coefficient(n1(xi,tau)*diff(Hn0(xi,tau),xi))

n1dn2 = nlsone(xi,tau).coefficient(n1(xi,tau)*v1(xi,tau)*diff(n1(xi,tau),xi)).collect(a0)
n12dv1 = nlsone(xi,tau).coefficient((n1(xi,tau))^2*diff(v1(xi,tau),xi)).collect(a0)
n1dHn2 = nlsone(xi,tau).coefficient(n1(xi,tau)*diff(Hn2(xi,tau),xi))

show(nlszero(xi,tau))
print
show(nlsone(xi,tau))
print
print "The ternary interactions."
show(n12v1)
print
show(n1dn2)
print
show(n12dv1)
print
show(n1dHn2)
print
print "The interactions with the mean."
show(n1n0)
print
show(n1dn0)
print
show(dn1n0)
print
show(n1dHn0)
print
︡9b30ec1d-de17-4822-a2e9-2fc2b5215fdd︡{"html":"<div align='center'>$\\displaystyle -\\frac{3}{2} i \\, k_{0}^{5} \\sigma n_{1}\\left(\\xi, \\tau\\right)^{2} v_{1}\\left(\\xi, \\tau\\right) + 2 i \\, \\Omega^{2} k_{0}^{2} s n_{1}\\left(\\xi, \\tau\\right)^{2} v_{1}\\left(\\xi, \\tau\\right) - 4 i \\, \\Omega a_{0} k_{0} s w n_{1}\\left(\\xi, \\tau\\right)^{2} v_{1}\\left(\\xi, \\tau\\right) - i \\, k_{0}^{2} s w^{2} n_{1}\\left(\\xi, \\tau\\right)^{2} v_{1}\\left(\\xi, \\tau\\right) + 2 i \\, \\Omega^{2} a_{0} k_{0} n_{1}\\left(\\xi, \\tau\\right)^{2} v_{1}\\left(\\xi, \\tau\\right) + i \\, a_{0} k_{0} w^{2} n_{1}\\left(\\xi, \\tau\\right)^{2} v_{1}\\left(\\xi, \\tau\\right) - 2 i \\, \\Omega k_{0} s w n_{0}\\left(\\xi, \\tau\\right) n_{1}\\left(\\xi, \\tau\\right) + i \\, k_{0} w^{2} n_{0}\\left(\\xi, \\tau\\right) n_{1}\\left(\\xi, \\tau\\right) + i \\, c_{g}^{2} s \\frac{\\partial^{2}}{(\\partial \\xi)^{2}}n_{1}\\left(\\xi, \\tau\\right) - 3 i \\, k_{0} \\sigma \\frac{\\partial^{2}}{(\\partial \\xi)^{2}}n_{1}\\left(\\xi, \\tau\\right) - 2 \\, \\Omega s \\frac{\\partial}{\\partial \\tau}n_{1}\\left(\\xi, \\tau\\right) + w \\frac{\\partial}{\\partial \\tau}n_{1}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle 4 \\, \\Omega c_{g} k_{0}^{2} s n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 9 \\, k_{0}^{4} \\sigma n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 4 \\, \\Omega a_{1} k_{0} s w n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 4 \\, a_{0} c_{g} k_{0} s w n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - \\frac{3}{2} \\, k_{0}^{4} \\sigma n_{1}\\left(\\xi, \\tau\\right)^{2} \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) - 2 \\, \\Omega^{2} a_{1} k_{0} n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 4 \\, \\Omega a_{0} c_{g} k_{0} n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 6 \\, \\Omega^{2} k_{0} s n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 8 \\, \\Omega a_{0} s w n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - a_{1} k_{0} w^{2} n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) - 3 \\, k_{0} s w^{2} n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 2 \\, \\Omega^{2} k_{0} s n_{1}\\left(\\xi, \\tau\\right)^{2} \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) - 4 \\, \\Omega a_{0} s w n_{1}\\left(\\xi, \\tau\\right)^{2} \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) - k_{0} s w^{2} n_{1}\\left(\\xi, \\tau\\right)^{2} \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) - 2 i \\, \\Omega c_{g} k_{0} s n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{0}\\left(\\xi, \\tau\\right) - 4 i \\, \\Omega k_{0} s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{2}\\left(\\xi, \\tau\\right) - c_{g} k_{0} s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{0}\\left(\\xi, \\tau\\right) - 2 \\, c_{g} k_{0} s w n_{0}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 4 \\, \\Omega^{2} a_{0} n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 2 \\, a_{0} w^{2} n_{1}\\left(\\xi, \\tau\\right) v_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 2 \\, \\Omega^{2} a_{0} n_{1}\\left(\\xi, \\tau\\right)^{2} \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) + a_{0} w^{2} n_{1}\\left(\\xi, \\tau\\right)^{2} \\frac{\\partial}{\\partial \\xi}v_{1}\\left(\\xi, \\tau\\right) + i \\, c_{g} k_{0} w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{0}\\left(\\xi, \\tau\\right) + 4 i \\, \\Omega^{2} k_{0} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{2}\\left(\\xi, \\tau\\right) + i \\, k_{0} w^{2} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}{\\rm Hn}_{2}\\left(\\xi, \\tau\\right) - 2 \\, \\Omega s w n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{0}\\left(\\xi, \\tau\\right) - 2 \\, \\Omega s w n_{0}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + w^{2} n_{1}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{0}\\left(\\xi, \\tau\\right) + w^{2} n_{0}\\left(\\xi, \\tau\\right) \\frac{\\partial}{\\partial \\xi}n_{1}\\left(\\xi, \\tau\\right) + 2 i \\, c_{g} s \\frac{\\partial^{2}}{\\partial \\xi\\partial \\tau}n_{1}\\left(\\xi, \\tau\\right) - \\sigma \\frac{\\partial^{3}}{(\\partial \\xi)^{3}}n_{1}\\left(\\xi, \\tau\\right)$</div>"}︡{"stdout":"\n"}︡{"stdout":"The ternary interactions.\n"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{3}{2} i \\, k_{0}^{5} \\sigma + 2 i \\, \\Omega^{2} k_{0}^{2} s - 4 i \\, \\Omega a_{0} k_{0} s w - i \\, k_{0}^{2} s w^{2} + 2 i \\, \\Omega^{2} a_{0} k_{0} + i \\, a_{0} k_{0} w^{2}$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle 4 \\, \\Omega c_{g} k_{0}^{2} s - 9 \\, k_{0}^{4} \\sigma + 4 \\, \\Omega a_{1} k_{0} s w - 2 \\, \\Omega^{2} a_{1} k_{0} + 6 \\, \\Omega^{2} k_{0} s - a_{1} k_{0} w^{2} - 3 \\, k_{0} s w^{2} - {\\left(4 \\, c_{g} k_{0} s w - 4 \\, \\Omega c_{g} k_{0} + 8 \\, \\Omega s w - 4 \\, \\Omega^{2} - 2 \\, w^{2}\\right)} a_{0}$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{3}{2} \\, k_{0}^{4} \\sigma + 2 \\, \\Omega^{2} k_{0} s - k_{0} s w^{2} - {\\left(4 \\, \\Omega s w - 2 \\, \\Omega^{2} - w^{2}\\right)} a_{0}$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle -4 i \\, \\Omega k_{0} s w + 4 i \\, \\Omega^{2} k_{0} + i \\, k_{0} w^{2}$</div>"}︡{"stdout":"\n"}︡{"stdout":"The interactions with the mean.\n"}︡{"html":"<div align='center'>$\\displaystyle -2 i \\, \\Omega k_{0} s w + i \\, k_{0} w^{2}$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle -c_{g} k_{0} s w - 2 \\, \\Omega s w + w^{2}$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle -2 \\, c_{g} k_{0} s w - 2 \\, \\Omega s w + w^{2}$</div>"}︡{"stdout":"\n"}︡{"html":"<div align='center'>$\\displaystyle -2 i \\, \\Omega c_{g} k_{0} s + i \\, c_{g} k_{0} w$</div>"}︡{"stdout":"\n"}︡{"done":true}︡
︠78deadd1-2358-474d-b301-2de5a5dba1f9︠









