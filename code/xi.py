import numpy as np, sympy as sp
# SM, symmetric phase, 3 generations, one Higgs doublet, all Yukawas + QCD sphalerons in equilibrium.
# unknowns (flavour-universal): muQ, muL, muH ; muu=muQ+muH, mud=muQ-muH, mue=muL-muH
mQ,mL,mH,mcs=sp.symbols('mQ mL mH mcs')
mu_=mQ+mH; md_=mQ-mH; me_=mL-mH
ng=3
# number densities n = (T^2/6) g mu for fermions (g internal dof), bosons (T^2/3) g mu
fer={'Q':(6,mQ,sp.Rational(1,6)),'u':(3,mu_,sp.Rational(2,3)),'d':(3,md_,-sp.Rational(1,3)),
     'L':(2,mL,-sp.Rational(1,2)),'e':(1,me_,-1)}
Y_neutral=sum(ng*g*m*Y for g,m,Y in fer.values())/6 + (2*mH*sp.Rational(1,2))/3   # in units of T^2
# weak sphaleron bias: sum over 12 doublet fermions created = mcs  => ng*(3 mQ + mL) = mcs
eqs=[sp.Eq(Y_neutral,0), sp.Eq(ng*(3*mQ+mL),mcs)]
sol=sp.solve(eqs,[mQ,mL],dict=True)[0]
# H free? need third condition: SU(2)xU(1) neutrality of T3 not needed in symmetric phase; Higgs chem pot fixed by Yukawa+hypercharge -> mH from hypercharge:
# we already used hypercharge neutrality; remaining unknown mH is fixed by requiring zero net B-L-independent? In the symmetric phase mH is determined by hypercharge neutrality;
# so solve jointly for (mQ,mL,mH) using also zero SU(2)-isospin charge (automatic) -> use B-L fixed instead? Here B-L is NOT fixed; add constraint of conserved B-L = 0 initial:
BmL=sp.Rational(1,3)*ng*(6*mQ+3*mu_+3*md_)/6 - ng*(2*mL+me_)/6
eqs2=[sp.Eq(Y_neutral,0), sp.Eq(ng*(3*mQ+mL),mcs), sp.Eq(BmL,0)]
s2=sp.solve(eqs2,[mQ,mL,mH],dict=True)[0]
B=sp.Rational(1,3)*ng*(6*mQ+3*mu_+3*md_)/6
Bv=sp.simplify(B.subs(s2)); print("B/(T^2) per unit mu_CS (B-L=0):",Bv)
# also equilibrium B for fixed B-L: ratio B/(B-L) with mcs=0
Bl=sp.symbols('Bl')
s3=sp.solve([sp.Eq(Y_neutral,0),sp.Eq(ng*(3*mQ+mL),0),sp.Eq(BmL,Bl)],[mQ,mL,mH],dict=True)[0]
print("B/(B-L) in symmetric phase:",sp.simplify(B.subs(s3)/Bl),"(expected 28/79=%.4f)"%(28/79))
g=106.75
Bv=sp.simplify(Bv/mcs)        # n_B/(T^2 mu_CS) = 12/79
# n_B = Bv T^2 mu ; s=(2 pi^2/45) g T^3 -> Y = Bv*mu/T / (2 pi^2 g/45)
xi_exact=float(Bv)/(2*np.pi**2*g/45)
print("n_B = %.4f T^2 mu_CS ;  xi = Y_B/(mu_CS/T) = %.3e (g*=106.75)"%(float(Bv),xi_exact))
print("paper-style xi = %.3e"%(45/(2*np.pi**2*g)*18/79))
