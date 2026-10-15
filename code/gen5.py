import numpy as np, json, warnings
warnings.filterwarnings("ignore")
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
exec(open('emde.py').read().split("def run(")[0])
M=2.435e18; Yobs=8.6e-11
def sci(x,p=1):
    e=int(np.floor(np.log10(abs(x)))); return r"%.*f\times10^{%d}"%(p,x/10**e,e)
rows=sorted([r for r in json.load(open('lvs_explicit.json')) if r['V']>1e3],key=lambda r:r['V'])
lV=np.log([r['V'] for r in rows]); lm=np.log([r['m_light'] for r in rows])
m_swiss=lambda V: np.exp(np.interp(np.log(V),lV,lm))
m_k3=lambda c: (lambda V: c*M*V**(-5/3))
def TR_of(mfun,V,cG=1.0):
    mm=mfun(V); TR=0.2
    for _ in range(60): TR=(90/(np.pi**2*g(max(TR,1e-4))))**0.25*np.sqrt(cG*mm**3/(4*np.pi*M))
    return mm,TR
def Vmax_TR(mfun,TRt):
    f=lambda lv: np.log(TR_of(mfun,10**lv)[1]/TRt)
    lo,hi=3.4,17
    if f(hi)>0: return 10**hi
    return 10**brentq(f,lo,hi)
def kappaY(TR,Tstart=3000.):
    G=Gamma_of_TR(TR); rR=np.pi**2/30*g(Tstart)*Tstart**4; H0=rR/(1.2*G*M**2); rphi=3*H0**2*M**2
    def f(N,y):
        a=rphi*np.exp(min(y[0],50)); b=rR*np.exp(min(y[1],50)); H=np.sqrt((a+b)/3)/M
        return [-3-G/H,-4+G*a/(H*b)]
    Nend=(8/3)*np.log(Tstart/TR)+14; N=np.linspace(0,Nend,60001)
    sol=solve_ivp(f,[0,Nend],[0,0],t_eval=N,method='Radau',rtol=1e-9,atol=1e-14)
    rp=rphi*np.exp(sol.y[0]); rr=rR*np.exp(sol.y[1]); T=T_of_rho(rr)
    S=(4/3)*rr/T*np.exp(3*N); i=np.where(rp/rr<1e-7)[0][0]; return rphi/(S[i]*TR)
TRk=np.logspace(-3,2,9); kY=np.array([kappaY(t) for t in TRk]); kYf=lambda TR: np.interp(np.log(TR),np.log(TRk),kY)
W=json.load(open('washout.json')); TRw=np.array([w[0] for w in W]); etaB=np.array([w[1] for w in W]); etaBf=lambda TR: np.interp(np.log(TR),np.log(TRw),etaB)
def Vmax_eps(mfun,eps):
    def f(lv):
        mm,TR0=TR_of(mfun,10**lv); TR=min(max(TR0,TRw[0]),TRw[-1])
        return np.log(eps*etaBf(TR)*kYf(TR)*TR0/mm/Yobs)
    lo,hi=3.4,17
    if f(lo)<0: return np.nan
    if f(hi)>0: return 10**hi
    return 10**brentq(f,lo,hi)
cases=[("Swiss-cheese (explicit vacua)",m_swiss),("K3-fibred, $c_f=0.7$",m_k3(0.695)),("K3-fibred, $c_f=1.45$ (reference)",m_k3(1.446)),("K3-fibred, $c_f=3.0$",m_k3(3.008)),("K3-fibred, $c_f=17.1$ (CBQ LV loop data)",m_k3(17.1))]
lines=[]
for name,mf in cases:
    vals=[Vmax_TR(mf,6.1),Vmax_TR(mf,100.),Vmax_TR(mf,4e-3)]+[Vmax_eps(mf,e) for e in (1,1e-2,1e-4)]
    lines.append(name+" & "+" & ".join("$%s$"%sci(v,1) for v in vals)+r" \\")
    print(name,["%.1e"%v for v in vals])
open('k3_table.tex','w').write("\n".join(lines)+"\n")
# K3 reheating at representative volumes
for V in [1e5,1e6,1e7,1e8]:
    print("V=%.0e"%V,"K3(c=1): m_f=%.2e T_R=%.2e | Swiss: m=%.2e T_R=%.2e"%(m_k3(1)(V),TR_of(m_k3(1),V)[1],m_swiss(V),TR_of(m_swiss,V)[1]))

# figure T_R(V)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Vs=np.logspace(3.5,15,200)
fig,ax=plt.subplots(figsize=(5.6,3.8))
ax.loglog(Vs,[TR_of(m_swiss,V)[1] if V<=1e16 else np.nan for V in Vs],label='Swiss-cheese: volume modulus (explicit)')
lo=[TR_of(m_k3(0.695),V)[1] for V in Vs]; hi=[TR_of(m_k3(3.008),V)[1] for V in Vs]
ax.fill_between(Vs,lo,hi,alpha=.3,color='C1'); ax.loglog(Vs,[TR_of(m_k3(1.446),V)[1] for V in Vs],'C1-',label='K3-fibred: fibre modulus (explicit)')
for y,l in [(4e-3,'BBN'),(6.1,r'pre-decay EW'),(100.,'leptogenesis')]:
    ax.axhline(y,color='gray',ls=':'); ax.text(2e3,y*1.4,l,fontsize=7)
ax.set_xlabel(r'$\mathcal{V}$'); ax.set_ylabel(r'$T_R$ [GeV]'); ax.set_ylim(1e-30,1e12); ax.legend(fontsize=7,loc='lower left'); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig('fig_TRV.pdf'); plt.close()
# table of published sets
Kx=json.load(open('k3_explicit.json'))
def sc(x): 
    e=int(np.floor(np.log10(abs(x)))); return r"%.1f\times10^{%d}"%(x/10**e,e)
tabs="\n".join(r"%s & $%s$ & $%s$ & $%.1f$ & $%s$ & $%s$ & $%s$ & $%s$ \\"%(r['set'],sc(r['V']),sc(r['tau1']),r['tau3'],sc(r['m32']),sc(r['m_fibre']),sc(r['m_volume']),sc(r['m_tau3'])) for r in Kx)
open('k3_sets.tex','w').write(tabs+"\n")
print(tabs)
