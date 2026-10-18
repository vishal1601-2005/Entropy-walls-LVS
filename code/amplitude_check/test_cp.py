import numpy as np, time
from test_kw import *
r2=np.sqrt(2)
def buildN(m1,m2,mt,c,s,k1R,k1L,k2R,k2L,phi=(0,0),rho=r2):
    A=[[0,r2*k1L],[0,r2*k2L]]; B=[[r2*k1R,0],[r2*k2R,0]]
    u=np.array([c,s]); D=np.outer(u,u)/mt**2
    return Model(m1,m2,A,B,D,rho,phi)
if __name__=='__main__':
    m1,m2,mt=1.0,0.5,3.0
    x=m2**2/m1**2
    t0=time.time()
    for label,(c,s,k1R,k1L,k2R,k2L) in {"pure R":(1,0,1.0,0,1.3,0),"R,L mixed (f1 and f2)":(np.cos(.7),np.sin(.7),0.7,-0.4,1.3,-1.3)}.items():
        kt,kr,kd=kw(m1,m2,mt,c,s,k1R,k1L,k2R,k2L)
        M=buildN(m1,m2,mt,c,s,k1R,k1L,k2R,k2L,phi=(0,np.pi/4)); Gt=M.gamma_tsb(); Gd=M.gamma_diff()
        # KW pieces separately
        kd_f2=(m1**5/(768*pi**3*mt**4))*(m1*m2/(16*pi*mt**2))*(k1R**2*k2R**2*c**6+k1L**2*k2L**2*c**2*s**4)*f2(x)
        kd_f1=(m1**5/(768*pi**3*mt**4))*(m1**2/(20*pi*mt**2))*k1R*k1L*k2R*k2L*c**4*s**2*f1(x)
        print(label)
        print("  Gamma_tsb engine/KW (should be 1): %.4f"%(Gt/kt))
        print("  Gamma_diff engine = %.4e | KW f2-term %.4e  f1-term %.4e  sum %.4e | engine/KWsum = %.3f"%(Gd,kd_f2,kd_f1,kd_f2+kd_f1,Gd/(kd_f2+kd_f1)))
        for ph in (0.0,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2):
            Mp=buildN(m1,m2,mt,c,s,k1R,k1L,k2R,k2L,phi=(0,ph)); print("     phi2=%.3f  Gamma_diff=%.4e  (sin2phi=%.3f)"%(ph,Mp.gamma_diff(),np.sin(2*ph)))
    print("time %.1fs"%(time.time()-t0))
