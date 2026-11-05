"""Affleck-Dine (udd, n=6) maximal volume as a function of the operator scale M_* = zeta*M_Pl.
Case I (m_AD > m_mod): Y_B = eta * eps*|phi0|^2*T_AD/(4 m_AD M^2) * D,  D = (T_R/T_i)(1-r)/r (kappa_Y=3/4).
Inputs: explicit Swiss-cheese vacua (Table 2 of the paper) and the K3 fibre-modulus scaling."""
import numpy as np
from scipy.optimize import brentq
M=2.435e18; Yobs=8.6e-11; r=1/6.; g=106.75; c=1.0; lam=1.0; eta=28/79; gs=0.3
V2=np.array([1.0e16,4.3e14,8.9e12,4.1e11,7.9e9,4.2e8,7.1e6,4.8e5,1.8e4,2.5e3])
m2=np.array([2.4e-6,2.9e-4,1.0e-1,1.1e1,4.4e3,3.8e5,2.0e8,1.2e10,2.0e12,4.0e13])
T2=np.array([6.4e-19,8.6e-16,5.6e-12,6.2e-9,5.0e-5,4.1e-2,2.7e2,1.3e5,2.7e8,2.5e10])
lv=np.log(V2[::-1])
def swiss(V): return (np.exp(np.interp(np.log(V),lv,np.log(m2[::-1]))), np.exp(np.interp(np.log(V),lv,np.log(T2[::-1]))))
def k3(V,cf=1.446):
    m=cf*M*V**(-5/3); Gam=m**3/(4*np.pi*M**2); b=np.sqrt(np.pi**2*g/90); return m,np.sqrt(Gam*M/b)
def Ms(V): return gs**0.25*M/np.sqrt(4*np.pi*V)
def YB(V,zeta,eps,geom='swiss',rr=r):
    mm,TR=(swiss(V) if geom=='swiss' else k3(V))
    mAD=0.1*M/V; Mst=zeta*M
    phi2=np.sqrt((c/5)**0.5*mAD*Mst**3/lam)**1   # |phi0|^2 = ((c/5)^{1/2} mAD M*^3/lam)^{1/2}
    phi2=((c/5)**0.5*mAD*Mst**3/lam)**0.5
    TAD=(90/(np.pi**2*g))**0.25*np.sqrt(mAD*M)
    Ti=(90*(1-rr)/(np.pi**2*g))**0.25*np.sqrt(mm*M)
    YAD=eps*phi2*TAD/(4*mAD*M**2)
    D=(TR/Ti)*(1-rr)/rr
    return eta*YAD*D
def Vmax(zeta,eps,geom='swiss'):
    f=lambda lv_: np.log(YB(np.exp(lv_),zeta(np.exp(lv_)) if callable(zeta) else zeta,eps,geom)/Yobs)
    lo,hi=np.log(2.5e3 if geom=='swiss' else 1e3),np.log(1e9)
    if f(lo)<0: return np.nan
    if f(hi)>0: return np.exp(hi)
    return np.exp(brentq(f,lo,hi))
zc=lambda V: (Ms(V)**2/(((1/5)**0.5*0.1*(M/V)*1.0)**0.5*1.0))**(2/3)/M   # |phi0|=M_s  => zeta_c(V)
# |phi0|^2=((c/5)^.5 mAD M*^3)^.5 = Ms^2 -> M*^3 = Ms^4/((c/5)^.5 mAD)
zc=lambda V: (Ms(V)**4/(((1/5)**0.5)*0.1*(M/V)))**(1/3)/M
zs=lambda V: Ms(V)/M
if __name__=="__main__":
    print("check zeta=1: swiss eps=1.2,0.51,0.16 ->",["%.2e"%Vmax(1,e) for e in (1.2,0.51,0.16)],"(paper 1.1e6, 5.9e5, 2.8e5)")
    print("check zeta=1: K3    eps=1.2,0.51,0.16 ->",["%.2e"%Vmax(1,e,'k3') for e in (1.2,0.51,0.16)],"(paper 2.5e5, 1.5e5, 7.3e4)")
    for z in (1,0.3,0.1,0.03,0.01,0.003):
        print("zeta=%g swiss %s  k3 %s"%(z,["%.1e"%Vmax(z,e) for e in (1.2,0.51,0.16)],["%.1e"%Vmax(z,e,'k3') for e in (1.2,0.51,0.16)]))
    for name,fz in (("|phi0|=Ms",zc),("M*=Ms",zs)):
        print(name,"swiss",["%.1e"%Vmax(fz,e) for e in (1.2,0.51,0.16)],"k3",["%.1e"%Vmax(fz,e,'k3') for e in (1.2,0.51,0.16)])
    for V in (1e3,1e4,1e5,5.9e5): print("V=%.1e zeta_c=%.4f zeta_s=%.2e |phi0|/Ms(zeta=1)=%.1f"%(V,zc(V),zs(V),(((1/5)**0.5*0.1*(M/V)*M**3)**0.5)**0.5/Ms(V)))
