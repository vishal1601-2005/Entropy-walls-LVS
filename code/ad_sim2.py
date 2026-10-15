import numpy as np
from scipy.integrate import solve_ivp
# EOM (L=|dphi|^2, units M=1, t in 1/m, mu=m/M), A-terms: low-energy a e^{i alpha} and Hubble-induced aH*h e^{i alphaH}
# f'' + 3h f' + (1-c h^2) f + (n-1) lam^2 |f|^(2n-4) f/mu^2 + lam f*^(n-1) conj(a_tot(h))/mu = 0 ,  a_tot=a e^{i alpha}+aH h e^{i alphaH}
def simulate2(n,c=1.0,lam=1.0,a=1.0,alpha=np.pi/2,aH=1.0,alphaH=0.0,mu=1e-2,hi=30.,hf=0.02):
    def f0(h): return (np.sqrt(c/(n-1))*h*mu/lam)**(1/(n-2))
    ti=1/(2*hi); tf=1/(2*hf)
    def atot(h): return a*np.exp(1j*alpha)+aH*h*np.exp(1j*alphaH)
    def rhs(t,y):
        h=1/(2*t); f=complex(y[0],y[1]); fp=complex(y[2],y[3])
        acc=-3*h*fp-(1-c*h*h)*f-(n-1)*lam**2*abs(f)**(2*n-4)*f/mu**2-lam*np.conj(f)**(n-1)*np.conj(atot(h))/mu
        return [y[2],y[3],acc.real,acc.imag]
    # start on the phase minimising Re(f^n atot(h)) : n*theta + arg(atot) = pi
    th0=(np.pi-np.angle(atot(hi)))/n
    # exact scale-free attractor f = R t^{-1/4} e^{i theta_min} (valid for h>>1), R from the algebraic radial equation
    from scipy.optimize import brentq
    def radial(R):
        k=(5/16-3/8-c/4) if n==6 else None
        # general n: f=R t^{-p}, p=1/(n-2)? use n=6 (p=1/4) and n=4 (p=1/2) separately
        return None
    p=1.0/(n-2)*1.0 if n!=6 else 0.25
    # f=R t^{-p}: f'=-p R t^{-p-1}, f''=p(p+1) R t^{-p-2}; h=1/(2t): 3h f'=-(3p/2) R t^{-p-2}; -c h^2 f=-(c/4) R t^{-p-2}
    # potential/Hubble-A terms scale as t^{-p(2n-3)} and t^{-p(n-1)-1}; scale-free requires p(2n-3)=p+2 -> p=1/(n-2)
    p=1.0/(n-2)
    k=p*(p+1)-1.5*p-c/4.0
    def G(R): return k*R+(n-1)*lam**2*R**(2*n-3)/mu**2-lam*aH*R**(n-1)/(2*mu)
    # find positive root bracket
    Rs=np.geomspace(1e-6,10,4000); vals=np.array([G(r) for r in Rs]); idx=np.where(np.sign(vals[:-1])!=np.sign(vals[1:]))[0]
    R=brentq(G,Rs[idx[-1]],Rs[idx[-1]+1]) if len(idx)>0 else f0(hi)*ti**p
    th0=(np.pi-np.angle(atot(hi)))/n
    fi=R*ti**(-p)*np.exp(1j*th0); vi=-p*fi/ti
    ts=np.geomspace(ti,tf,4000)
    sol=solve_ivp(rhs,[ti,tf],[fi.real,fi.imag,vi.real,vi.imag],t_eval=ts,method='DOP853',rtol=1e-11,atol=1e-15)
    f=sol.y[0]+1j*sol.y[1]; fp=sol.y[2]+1j*sol.y[3]
    Q=(sol.t**1.5)*(1j*(np.conj(fp)*f-np.conj(f)*fp)).real
    ref=(0.5)**1.5*f0(1.0)**2
    return sol.t,Q/ref
if __name__=="__main__":
    for hi in (15.,30.,60.,100.,200.):
        t,e=simulate2(6,hi=hi); print("benchmark a=aH=c=1, delta=pi/2: hi=%g eps=%.4f std=%.1e"%(hi,e[-300:].mean(),e[-300:].std()))
