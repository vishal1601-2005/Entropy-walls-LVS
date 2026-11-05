"""Higgs-sector couplings of the light volume modulus in the explicit Swiss-cheese vacua.
Kane-Winkler (arXiv:1909.04705) Eqs. (12)-(16): Gamma_hh = m |d_phi mu|^2/(4 pi), Gamma_HH = (|d_phi mH^2|^2+|d_phi B mu|^2)/(4 pi m),
with mu = (m32 Z - F^mbar d_mbar Z)/Zhat etc.  Matter metric Zhat ~ V^(-2/3); Z = z0 * Zhat.  d_phi = v^k d/d tau_k along the light eigenvector.
All quantities in M_Pl = 1 units; ratios to Gamma_aa = m^3/(48 pi)."""
from lvs_decay_widths import *
from mpmath import mpf, diff, log, sqrt, pi, exp
def lnZh(tb,ts): return -mpf(2)/3*log(volf(tb,ts))
def d_lnZ(tb,ts):   # holomorphic derivatives d_m = (1/2) d/dtau_m of ln Zhat
    return [diff(lambda x:lnZh(x,ts),tb)/2, diff(lambda y:lnZh(tb,y),ts)/2]
def dd_lnZ(tb,ts): # d_mbar d_n ln Zhat = (1/4) d^2/dtau_m dtau_n
    H=[[diff(lnZh,(tb,ts),(2,0)),diff(lnZh,(tb,ts),(1,1))],[diff(lnZh,(tb,ts),(1,1)),diff(lnZh,(tb,ts),(0,2))]]
    return [[H[i][j]/4 for j in range(2)] for i in range(2)]
def m32f(tb,ts,th,a):
    Kb,Kmat,Kinv,W,Wi,DW=pieces(tb,ts,th,a); return exp(Kc(tb,ts)/2)*abs(W)
def mu_hat(tb,ts,th,a):          # mu/z0
    F=Fterm(tb,ts,th,a); d=d_lnZ(tb,ts); m32=m32f(tb,ts,th,a)
    return (m32-sum(F[i].real*d[i] for i in range(2)))
def mH2_hat(tb,ts,th,a):         # (m_H^2 - |mu|^2 ) without V0 : m32^2 - F F d d ln Zhat
    F=Fterm(tb,ts,th,a); dd=dd_lnZ(tb,ts); m32=m32f(tb,ts,th,a)
    return m32**2-sum(F[i].real*F[j].real*dd[i][j] for i in range(2) for j in range(2))
def Bmu_hat(tb,ts,th,a):         # B mu / z0
    F=Fterm(tb,ts,th,a); d=d_lnZ(tb,ts); dd=dd_lnZ(tb,ts); m32=m32f(tb,ts,th,a)
    t1=2*m32**2-2*m32*sum(F[i].real*d[i] for i in range(2))
    t2=-sum(F[i].real*F[j].real*(dd[i][j]-d[i]*d[j]) for i in range(2) for j in range(2))
    return t1+t2
def analyze_h(aval):
    a,ts,tb=solve_a(aval); th=pi/a; Vt=volf(tb,ts)
    Ki,Kij,Yi,Y=Kder(tb,ts)
    G=matrix(2,2)
    for i in range(2):
        for j in range(2): G[i,j]=2*Kij[i][j]/4
    f=lambda x,y:Vpot(x,y,th,a); H=matrix(2,2)
    H[0,0]=diff(f,(tb,ts),(2,0)); H[1,1]=diff(f,(tb,ts),(0,2)); H[0,1]=H[1,0]=diff(f,(tb,ts),(1,1))
    M=inverse(G)*H; tr=M[0,0]+M[1,1]; det=M[0,0]*M[1,1]-M[0,1]*M[1,0]; disc=sqrt(tr**2-4*det); l0=(tr-disc)/2
    v=matrix([M[0,1],l0-M[0,0]]); v=v/sqrt((v.T*G*v)[0,0]); m=sqrt(l0)
    dphi=lambda fn: v[0]*diff(lambda x:fn(x,ts,th,a),tb)+v[1]*diff(lambda y:fn(tb,y,th,a),ts)
    m32=m32f(tb,ts,th,a); V0=Vpot(tb,ts,th,a)
    mu=mu_hat(tb,ts,th,a); dmu=dphi(mu_hat); mH2=mH2_hat(tb,ts,th,a)+V0; dmH2=dphi(mH2_hat); dB=dphi(Bmu_hat); Bmu=Bmu_hat(tb,ts,th,a)
    Gaa=m**3/(48*pi)
    Gh=m*dmu**2/(4*pi); GH=(dmH2**2+dB**2)/(4*pi*m)
    return dict(V=Vt,m=m,m32=m32,mu=mu,dmu=dmu,mH2=mH2,Bmu=Bmu,Gh_over_aa=Gh/Gaa,GH_over_aa=GH/Gaa,V0=V0)
if __name__=="__main__":
    print("units M_Pl=1; z0=Z/Zhat; entries are per unit z0 (rates scale as z0^2)")
    print("%9s %10s %10s | %12s %12s %12s | %12s %12s | %10s"%("V","m_Phi/m32","mu/m32","mu*V/m32","dphi mu/m_Phi","mH/m32 *sqrtV","Gam_hh/Gam_aa","Gam_HH/Gam_aa","(mH^2+|mu|^2)>0"))
    for aval in (1.0,1.2,1.3,1.4,1.5):
        r=analyze_h(aval); V=r['V']
        print("%9.2e %10.2e %10.2e | %12.3f %12.3e %12.3f | %12.3e %12.3e | %s"%(float(V),float(r['m']/r['m32']),float(r['mu']/r['m32']),float(r['mu']*V/r['m32']),float(r['dmu']/r['m']),float(sqrt(abs(r['mH2']))/r['m32']*sqrt(V)),float(r['Gh_over_aa']),float(r['GH_over_aa']),"m_H^2=%.2e m32^2"%float(r['mH2']/r['m32']**2)))
