"""Entropy-wall quantities for the UV constructions of Kane & Winkler (arXiv:1909.04705, Table 1) and sequestered LVS.
Gamma = c m^3/M_Pl^2 (their convention), g_*(T_R) and T_*(T_R) iterated self-consistently.
Output reproduces Table 'Entropy-wall quantities in four moduli-dominated cosmologies' of the paper."""
import numpy as np
M=2.435e18; Yobs=8.6e-11
Tg=np.array([1e-4,0.05,0.2,0.5,1,3,10,100,170]); gg=np.array([10.75,10.75,50,65,75,80,86,96,106.75])
gR=lambda T: np.interp(np.log(T),np.log(Tg),gg)
TRs=np.array([0.018,0.056,1.0,10.0,100.0]); Ts=np.array([156,153,145,139,131.7]); Tst=lambda T: np.interp(np.log(T),np.log(TRs),Ts)
def TR(m,c):
    T=0.02
    for _ in range(60):
        g=gR(T); b=np.sqrt(np.pi**2*g/90); T=np.sqrt(c*m**3/M/b)
    return T,gR(T)
def row(m,c):
    T,g=TR(m,c); Y=0.75*T/m; D=0.16*g/106.75*(T/131.7)**5; eL=28/79*4/15*g/106.75*(T/Tst(T))**4
    return T,Y,D,eL,Yobs/Y
cases={"IIB/KKLT":[(1e6,.13),(8e6,.13)],"heterotic":[(2e6,.21),(6e6,.21)],"G2":[(6e4,.1),(1.6e5,.1),(6e4,36),(1.6e5,36)]}
print("%-10s  T_R[GeV] (min,max)   D(131.7 GeV) (min,max)   eta_L (min,max)   eps_req (min,max)"%"model")
for name,cs in cases.items():
    R_=np.array([row(*c) for c in cs])
    rng=lambda k:(R_[:,k].min(),R_[:,k].max())
    print("%-10s  %.2g--%.2g         %.0e--%.0e         %.0e--%.0e        %.1e--%.1e"%((name,)+rng(0)+rng(2)+rng(3)+rng(4)))
# threshold modulus mass for T_R = 6.1 GeV
for c in (0.1,0.13,0.21,36):
    m=1e6
    for _ in range(80):
        T,g=TR(m,c); m*=(6.1/T)**(2/3)
    print("c_tot=%.2f: m_Phi(T_R=6.1 GeV) = %.2e GeV"%(c,m))
# sequestered LVS window from the explicit vacua (Table 2 of the paper)
V2=np.array([1.0e16,4.3e14,8.9e12,4.1e11,7.9e9,4.2e8,7.1e6,4.8e5,1.8e4,2.5e3]);m2=np.array([2.4e-6,2.9e-4,1.0e-1,1.1e1,4.4e3,3.8e5,2.0e8,1.2e10,2.0e12,4.0e13]);T2=np.array([6.4e-19,8.6e-16,5.6e-12,6.2e-9,5.0e-5,4.1e-2,2.7e2,1.3e5,2.7e8,2.5e10])
lv=np.log(V2[::-1]); f=lambda y,v:np.exp(np.interp(np.log(v),lv,np.log(y[::-1])))
for v in (1.5e7,1.5e8):
    m=f(m2,v);T=f(T2,v);Y=0.75*T/m;g=gR(T);D=0.16*g/106.75*(T/131.7)**5
    print("sequestered LVS V=%.1e: m=%.1e GeV, T_R=%.2g GeV, D=%.1e, eps_req=%.1e"%(v,m,T,D,Yobs/Y))
