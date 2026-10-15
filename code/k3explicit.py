import sympy as sp, mpmath as mp, json, sys
mp.mp.dps=60
t1,Vv,t3,t2,dup=sp.symbols('t1 Vv t3 t2 dup',positive=True)
al,ga,xi,gs,W0,a3,A3,A,B,C=sp.symbols('al ga xi gs W0 a3 A3 A B C',positive=True)
xh=xi/gs**sp.Rational(3,2)
mu1=8*a3**2*A3**2/(3*al*ga); mu2=4*a3*A3; mu3=3*xh/4
# CBQ eq (3.56): V in M_P=1 units, before overall normalisation g_s e^{K_cs}/(8 pi)
Vhat=mu1*sp.sqrt(t3)/Vv*sp.exp(-2*a3*t3)-mu2*W0*t3*sp.exp(-a3*t3)/Vv**2+mu3*W0**2/Vv**3+dup/Vv**sp.Rational(4,3)+(A/t1**2-B/(Vv*sp.sqrt(t1))+C*t1/Vv**2)*W0**2/Vv**2
# Kahler geometry: x=(t1,t2,t3), Volume = al(sqrt(t1) t2 - ga t3^{3/2}); K=-2 ln(Volume + xh/2)
x=(t1,t2,t3)
Vol=al*(sp.sqrt(t1)*t2-ga*t3**sp.Rational(3,2))
K=-2*sp.log(Vol+xh/2)
Kij=sp.Matrix(3,3,lambda i,j: sp.diff(K,x[i],x[j])/4)
# change of coordinates y=(t1,Vv,t3): t2 = (Vv/al + ga t3^{3/2})/sqrt(t1)
t2y=(Vv/al+ga*t3**sp.Rational(3,2))/sp.sqrt(t1)
y=(t1,Vv,t3)
J=sp.Matrix(3,3,lambda i,j: sp.diff([t1,t2y,t3][i],y[j]))
Kyy=(J.T*Kij.subs(t2,t2y)*J)
pars=(al,ga,xi,gs,W0,a3,A3,A,B,C)
fV=sp.lambdify((t1,Vv,t3,dup)+pars,Vhat,'mpmath')
grads=[sp.lambdify((t1,Vv,t3,dup)+pars,sp.diff(Vhat,v),'mpmath') for v in y]
hess=[[sp.lambdify((t1,Vv,t3,dup)+pars,sp.diff(Vhat,v,w),'mpmath') for w in y] for v in y]
fK=sp.lambdify((t1,Vv,t3)+pars,Kyy,'mpmath')
SETS={
 'LV': dict(al=.5,ga=.385,xi=.409,gs=.1,W0=1,a3=mp.pi,A3=1,A=1e-4,B=10,C=5e-5,V0=2.75e13,t30=10.46),
 'SV1':dict(al=.1291,ga=3.651,xi=.934,gs=.3,W0=100,a3=mp.pi/5,A3=1,A=2e-3,B=.52,C=1.9e-5,V0=1709.55,t30=4.28),
 'SV2':dict(al=.1543,ga=3.055,xi=.755,gs=.3,W0=100,a3=mp.pi/4,A3=1,A=2.9e-3,B=.93,C=4.3e-5,V0=1626.12,t30=3.73)}
M=2.435e18
def solve(name):
    S=SETS[name]; P=tuple(mp.mpf(S[k]) for k in ('al','ga','xi','gs','W0','a3','A3','A','B','C'))
    t10=(4*P[7]*S['V0']/P[8])**(mp.mpf(2)/3)
    def F(l1,lV,l3,d):
        a,b,c=mp.e**l1,mp.e**lV,mp.e**l3
        sc=mp.mpf(1)/ (abs(fV(a,b,c,0,*P))+mp.mpf(10)**-300)
        return [grads[0](a,b,c,d,*P)*a*sc, grads[1](a,b,c,d,*P)*b*sc, grads[2](a,b,c,d,*P)*c*sc, fV(a,b,c,d,*P)*sc]
    # initial dup from V=0 at start
    a,b,c=t10,mp.mpf(S['V0']),mp.mpf(S['t30'])
    d0=-(fV(a,b,c,0,*P))*b**(mp.mpf(4)/3)
    sol=mp.findroot(F,(mp.log(a),mp.log(b),mp.log(c),d0),tol=1e-30,maxsteps=100)
    a,b,c,d=mp.e**sol[0],mp.e**sol[1],mp.e**sol[2],sol[3]
    H=mp.matrix([[hess[i][j](a,b,c,d,*P) for j in range(3)] for i in range(3)])
    G=2*mp.matrix(fK(a,b,c,*P))
    Mm=G**-1*H
    ev=sorted([mp.re(e) for e in mp.eig(Mm)[0]])
    gs_=S['gs']; NV=gs_*mp.e**(3*mp.log(gs_))/(8*mp.pi)    # g_s e^{K_cs}/(8 pi), K_cs=3 ln g_s
    masses=[float(mp.sqrt(e*NV))*M for e in ev]
    # analytic fibre mass from CBQ: m^2 = 4(1+7R/6) C2/V^{10/3}, C2=A W0^2 zeta^{-4/3}, zeta=4A/B
    zeta=4*P[7]/P[8]; C2=P[7]*P[4]**2*zeta**(-mp.mpf(4)/3); R=16*P[7]*P[9]/P[8]**2
    m_an=float(mp.sqrt(4*(1+7*R/6)*C2/b**(mp.mpf(10)/3)*NV))*M
    m32=float(mp.sqrt(NV)*P[4]/b)*M
    return dict(set=name,tau1=float(a),V=float(b),tau3=float(c),NV=float(NV),m_fibre=masses[0],m_volume=masses[1],m_tau3=masses[2],m_fibre_analytic=m_an,m32=m32,R=float(R),eig_raw=[float(e) for e in ev])
if __name__=="__main__":
    out=[]
    for n in ('LV','SV1','SV2'):
        try:
            r=solve(n); out.append(r); print(r)
        except Exception as e: print(n,"fail",e)
    json.dump(out,open('k3_explicit.json','w'))
