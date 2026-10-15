import numpy as np, json, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
exec(open('emde.py').read().split("def run(")[0])
M=2.435e18; Yobs=8.6e-11; aW=0.0338; Hstar=2.4e-14
def sci(x,p=1):
    e=int(np.floor(np.log10(abs(x)))); return r"%.*f\times10^{%d}"%(p,x/10**e,e)
mac=[]
def m_(n,v): mac.append(r"\newcommand{\%s}{%s}"%(n,v))
rows=[r for r in json.load(open('lvs_explicit.json')) if r['V']>1e3]
rows.sort(key=lambda r:r['V'])
lV=np.log(np.array([r['V'] for r in rows])); lm=np.log(np.array([r['m_light'] for r in rows]))
m_of_V=lambda V: np.exp(np.interp(np.log(V),lV,lm))
cm=np.array([r['m_light']*r['V']**1.5/M for r in rows]); print("c_m=m V^1.5/M:",np.round(cm,2))
m_('cmmin',"%.1f"%cm.min()); m_('cmmax',"%.1f"%cm.max())
def TR_of_V(V,cG=1.0):
    mm=m_of_V(V); TR=0.2
    for _ in range(60): TR=(90/(np.pi**2*g(max(TR,1e-4))))**0.25*np.sqrt(cG*mm**3/(4*np.pi*M))
    return mm,TR
# ---- explicit table
tab=[]
for r in sorted(rows,key=lambda r:-r['V']):
    V=r['V']; mm,TR=TR_of_V(V)
    Nsm=(1/aW)/r['ts']; Fs=M*np.sqrt(2*r['Kss'])/(2*np.pi*Nsm)
    Gam=(3/8)*aW**2*r['m_rho']**3/(32*np.pi**3*Fs**2)
    tab.append((V,r['ts'],r['m32'],mm,r['m_rho'],TR,Fs,Gam/Hstar))
tex="\n".join(r"$%s$ & $%.1f$ & $%s$ & $%s$ & $%s$ & $%s$ \\"%(sci(V),ts,sci(m32),sci(mm),sci(mr),sci(TR)) for V,ts,m32,mm,mr,TR,Fs,gh in tab)
open('lvs_explicit_table.tex','w').write(tex+"\n")
tex2="\n".join(r"$%s$ & $%s$ & $%s$ \\"%(sci(V),sci(Fs),sci(gh)) for V,ts,m32,mm,mr,TR,Fs,gh in tab)
open('fs_explicit_table.tex','w').write(tex2+"\n")
m_('GHminExp',sci(min(t[7] for t in tab),1)); m_('FsExpOne',sci(tab[0][6],1))
# ---- V_max bounds from explicit spectrum
def Vmax_for_TR(TRt,cG=1.0):
    lo,hi=np.log10(rows[0]['V']),np.log10(rows[-1]['V'])
    f=lambda lv: np.log(TR_of_V(10**lv,cG)[1]/TRt)
    if f(hi)>0: return 10**hi
    if f(lo)<0: return np.nan
    return 10**brentq(f,lo,hi)
for nm,tr in [('a',6.1),('b',37.7),('BBN',4e-3),('lep',100.),('ten',10.)]:
    m_('VexpTR'+nm,sci(Vmax_for_TR(tr),1))
m_('VexpTRaLo',sci(Vmax_for_TR(6.1,8.0),1)); m_('VexpTRaHi',sci(Vmax_for_TR(6.1,0.3),1))
# ---- decay baryogenesis V_max(eps) explicit
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
def Vmax_eps(eps,cG=1.0):
    def f(lv):
        mm,TR0=TR_of_V(10**lv,cG); TR=min(max(TR0,TRw[0]),TRw[-1])
        return np.log(eps*etaBf(TR)*kYf(TR)*TR0/mm/Yobs)
    lo,hi=np.log10(rows[0]['V']),np.log10(rows[-1]['V'])
    if f(lo)<0: return np.nan
    if f(hi)>0: return 10**hi
    return 10**brentq(f,lo,hi)
eps=np.logspace(-8,0,25)
V1=np.array([Vmax_eps(e) for e in eps]); Vlo=np.array([Vmax_eps(e,8.0) for e in eps]); Vhi=np.array([Vmax_eps(e,0.3) for e in eps])
for e,nm in [(1,'one'),(1e-2,'two'),(1e-4,'four'),(1e-6,'six')]:
    v=Vmax_eps(e); m_('VexpEps'+nm,sci(v,1) if v==v else r"\text{none}")
fig,ax=plt.subplots(figsize=(5.4,3.7))
ax.fill_between(eps,Vlo,Vhi,alpha=.25,color='C0'); ax.loglog(eps,V1,'C0-',label=r'$B$-violating decays: $\mathcal{V}_{\max}(\epsilon)$')
ax.axhline(Vmax_for_TR(6.1),color='C3',ls='--',label=r'pre-decay EW baryogenesis ($\mathcal{V}\lesssim3\times10^7$)')
ax.axhline(1e15,color='k',ls=':',label=r'TeV-scale LVS ($\mathcal{V}\sim10^{15}$)')
ax.set_xlabel(r'CP-violating baryon number per decay, $\epsilon$'); ax.set_ylabel(r'$\mathcal{V}_{\max}$'); ax.set_ylim(1e3,1e17); ax.legend(fontsize=6.5,loc='lower right'); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig('fig_map.pdf'); plt.close()
# closed forms check
g1,gR=106.75,10.75; kap=(np.pi**2*g1/36)/np.sqrt(np.pi**2*gR/90); c=0.118; Tc=159.; aw=0.0338
TReq=np.sqrt(c*kap*Tc**3/(18*aw**5*M)); m_('TReqclosed',"%.1f"%(TReq*1e3))
chk=[]
for TRq,Ts in [(0.018,156.),(1.0,145.),(10.0,139.)]:
    eta=(28/79)*(4/15)*(g(TRq)/g1)*(TRq/Ts)**4; chk.append((TRq,eta))
    print("eta_L closed form TR=%g: %.2e"%(TRq,eta))
m_('etaClosedA',sci(chk[0][1],1)); m_('etaClosedB',sci(chk[1][1],1)); m_('etaClosedC',sci(chk[2][1],1))
open('numbers4.tex','w').write("\n".join(mac)+"\n"); print(open('numbers4.tex').read()); print(tex)
