import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
# Same EOM as ad_sim2 but with a general power-law background H = k/t (k=1/2 radiation, k=2/3 matter/inflaton-oscillation domination).
def simulate3(n,k=2/3,c=1.0,lam=1.0,a=1.0,alpha=np.pi/2,aH=1.0,alphaH=0.0,mu=1e-2,hi=400.,hf=0.02):
    p=1.0/(n-2)
    def f0(h): return (np.sqrt(c/(n-1))*h*mu/lam)**(1/(n-2))
    ti=k/hi; tf=k/hf
    atot=lambda h: a*np.exp(1j*alpha)+aH*h*np.exp(1j*alphaH)
    def rhs(t,y):
        h=k/t; f=complex(y[0],y[1]); fp=complex(y[2],y[3])
        acc=-3*h*fp-(1-c*h*h)*f-(n-1)*lam**2*abs(f)**(2*n-4)*f/mu**2-lam*np.conj(f)**(n-1)*np.conj(atot(h))/mu
        return [y[2],y[3],acc.real,acc.imag]
    kc=p*(p+1)-3*k*p-c*k*k
    G=lambda R: kc*R+(n-1)*lam**2*R**(2*n-3)/mu**2-lam*aH*k*R**(n-1)/mu
    Rs=np.geomspace(1e-6,10,4000); v=np.array([G(r) for r in Rs]); idx=np.where(np.sign(v[:-1])!=np.sign(v[1:]))[0]
    R=brentq(G,Rs[idx[-1]],Rs[idx[-1]+1])
    th0=(np.pi-np.angle(atot(hi)))/n
    fi=R*ti**(-p)*np.exp(1j*th0); vi=-p*fi/ti
    ts=np.geomspace(ti,tf,4000)
    sol=solve_ivp(rhs,[ti,tf],[fi.real,fi.imag,vi.real,vi.imag],t_eval=ts,method='DOP853',rtol=1e-11,atol=1e-15)
    f=sol.y[0]+1j*sol.y[1]; fp=sol.y[2]+1j*sol.y[3]
    Q=(sol.t**(3*k))*(1j*(np.conj(fp)*f-np.conj(f)*fp)).real
    ref=(k)**(3*k)*f0(1.0)**2          # a^3 m |phi0|^2 at H=m  (t=k)
    return sol.t,Q/ref
if __name__=="__main__":
    for n in (6,4):
        for hi in (100.,400.,1600.):
            t,e=simulate3(n,hi=hi); print("MD background n=%d hi=%g eps=%.4f std=%.1e"%(n,hi,e[-300:].mean(),e[-300:].std()),flush=True)
