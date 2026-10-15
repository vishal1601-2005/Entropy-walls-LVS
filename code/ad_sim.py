import numpy as np, json
from scipy.integrate import solve_ivp
# units: M=1 (reduced Planck), time in 1/m_AD, field phi = M f ; mu=m/M.  EOM (L=|dphi|^2):
# f'' + 3h f' + (1-c h^2) f + (n-1) lam^2 |f|^(2n-4) f/mu^2 + aA*lam f*^(n-1) e^{-i alpha}/mu = 0
def simulate(n,c=1.0,lam=1.0,aA=1.0,alpha=np.pi/2,theta0=0.0,mu=1e-2,hi=30.,hf=0.02,rd=True):
    def f0(h): return (np.sqrt(c/(n-1))*h*mu/lam)**(1/(n-2))   # |f| tracking minimum (M=1,m=1 units)
    # radiation domination: a ~ t^{1/2}, h=1/(2t)
    ti=1/(2*hi); tf=1/(2*hf)
    def rhs(t,y):
        x,yy,vx,vy=y; h=1/(2*t)
        f=complex(x,yy); fp=complex(vx,vy)
        acc=-3*h*fp-(1-c*h*h)*f-(n-1)*lam**2*abs(f)**(2*n-4)*f/mu**2-aA*lam*np.conj(f)**(n-1)*np.exp(-1j*alpha)/mu
        return [vx,vy,acc.real,acc.imag]
    f_i=f0(hi)*np.exp(1j*theta0)
    y0=[f_i.real,f_i.imag,0.0,0.0]
    ts=np.geomspace(ti,tf,4000)
    sol=solve_ivp(rhs,[ti,tf],y0,t_eval=ts,method='DOP853',rtol=1e-10,atol=1e-14)
    x,yy,vx,vy=sol.y; f=x+1j*yy; fp=vx+1j*vy
    n_dens=1j*(np.conj(fp)*f-np.conj(f)*fp)         # i(phi_dot^* phi - phi^* phi_dot)
    a3=(sol.t)**1.5                                    # a ~ t^{1/2} -> a^3 ~ t^{3/2}
    Q=(a3*n_dens.real)
    # reference: m |phi0|^2 at h=1 (t=1/2) with a^3(t=1/2)
    ref=(0.5)**1.5*(f0(1.0)**2)
    return sol.t,Q/ref
if __name__=="__main__":
    res={}
    for n in (6,4):
        for alpha in (np.pi/2,np.pi/4,np.pi/8):
            t,eps=simulate(n,alpha=alpha)
            late=eps[-400:]
            print("n=%d alpha=%.2f: eps (late, comoving) mean=%.3f  spread=%.1e  max|eps|=%.3f"%(n,alpha,late.mean(),late.std(),np.abs(eps).max()))
            res[f"{n}_{alpha:.2f}"]=float(late.mean())
    # dependence on mu and aA, c for n=6
    for mu in (1e-2,1e-3):
        t,eps=simulate(6,mu=mu); print("n=6 mu=%g eps=%.3f"%(mu,eps[-400:].mean()))
    for aA in (0.3,1.0,3.0):
        t,eps=simulate(6,aA=aA); print("n=6 aA=%g eps=%.3f"%(aA,eps[-400:].mean()))
    for c in (0.3,1.0,3.0):
        t,eps=simulate(6,c=c); print("n=6 c=%g eps=%.3f"%(c,eps[-400:].mean()))
    json.dump(res,open('ad_sim.json','w'))
