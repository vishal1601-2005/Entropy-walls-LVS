import numpy as np, json, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
exec(open('emde.py').read().split("def run(")[0])   # g, T_of_rho, M, Gamma_of_TR
M=2.435e18; Yobs=8.6e-11; aW=0.0338; Tgen=131.7; gstar=106.75
xi=(12/79)/(2*np.pi**2*gstar/45)
mac=[]
def m(n,v): mac.append(r"\newcommand{\%s}{%s}"%(n,v))
def sci(x,p=1):
    e=int(np.floor(np.log10(abs(x)))); return r"%.*f\times10^{%d}"%(p,x/10**e,e)
# ---- kappa_Y: Y_dec = kappa_Y T_R/m  (decays per entropy)
def kappaY(TR,Tstart=3000.):
    G=Gamma_of_TR(TR); rR=np.pi**2/30*g(Tstart)*Tstart**4; H0=rR/(1.2*G*M**2); rphi=3*H0**2*M**2
    def f(N,y):
        a=rphi*np.exp(min(y[0],50)); b=rR*np.exp(min(y[1],50)); H=np.sqrt((a+b)/3)/M
        return [-3-G/H,-4+G*a/(H*b)]
    Nend=(8/3)*np.log(Tstart/TR)+14
    N=np.linspace(0,Nend,100001)
    sol=solve_ivp(f,[0,Nend],[0,0],t_eval=N,method='Radau',rtol=1e-10,atol=1e-14)
    rp=rphi*np.exp(sol.y[0]); rr=rR*np.exp(sol.y[1]); T=T_of_rho(rr)
    S=(4/3)*rr/T*np.exp(3*N); i=np.where(rp/rr<1e-7)[0][0]
    return rphi/(S[i]*TR)
TRk=np.logspace(-3,2,11); kY=np.array([kappaY(t) for t in TRk])
print("kappa_Y:",np.round(kY,3))
kYf=lambda TR: np.interp(np.log(TR),np.log(TRk),kY)
# ---- eta_B, eta_L from washout.json
W=json.load(open('washout.json')); TRw=np.array([w[0] for w in W]); etaB=np.array([w[1] for w in W]); etaL=np.abs([w[3] for w in W])
etaBf=lambda TR: np.interp(np.log(TR),np.log(TRw),etaB)
fig,ax=plt.subplots(figsize=(5.4,3.7))
ax.loglog(TRw,etaB,'o-',label=r'$B$-violating decays ($\eta_B$)')
ax.loglog(TRw,etaL,'s-',label=r'$L$-violating decays ($|\eta_L|$)')
ax.axhline(28/79,color='gray',ls=':'); ax.text(0.0035,0.4,r'$28/79$',fontsize=7)
ax.set_xlabel(r'$T_R$ [GeV]'); ax.set_ylabel('conversion efficiency'); ax.set_ylim(1e-23,3); ax.grid(alpha=.3); ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig('fig_eta.pdf'); plt.close()
# fit eta_L around 0.2-10 GeV
mk=(TRw>=1)&(TRw<=10); p=np.polyfit(np.log(TRw[mk]),np.log(etaL[mk]),1)
m('etaLexp',"%.1f"%p[0]); m('etaLone',sci(np.exp(p[1]),1))
m('etaLa',sci(etaL[3],1)); m('etaLb',sci(etaL[6],1)); m('etaLc',sci(etaL[8],1)); m('etaLd',"%.2f"%etaL[10]); m('etaLe',"%.2f"%etaL[11])
m('etaBc',"%.2f"%etaB[10]); m('etaBd',"%.2f"%etaB[11])
# ---- LVS: T_R(V), m(V); V_max(eps) from  eps * etaB * kappaY * TR/m = Yobs
def TRm(V,cG=1.0):
    mm=M*V**-1.5; TR=0.2
    for _ in range(60): TR=(90/(np.pi**2*g(max(TR,1e-4))))**0.25*np.sqrt(cG*mm**3/(4*np.pi*M))
    return mm,TR
def Vmax_eps(eps,cG=1.0):
    f=lambda lv: np.log(eps*etaBf(max(TRm(10**lv,cG)[1],TRw[0]))*kYf(max(TRm(10**lv,cG)[1],TRk[0]))*TRm(10**lv,cG)[1]/TRm(10**lv,cG)[0]/Yobs)
    # decreasing in V
    lo,hi=3,16
    if f(lo)<0: return np.nan
    if f(hi)>0: return 10**hi
    return 10**brentq(f,lo,hi)
eps=np.logspace(-8,0,25)
V1=np.array([Vmax_eps(e) for e in eps]); Vlo=np.array([Vmax_eps(e,1/(4*np.pi)) for e in eps]); Vhi=np.array([Vmax_eps(e,4*np.pi) for e in eps])
fig,ax=plt.subplots(figsize=(5.4,3.7))
ax.fill_between(eps,Vlo,Vhi,alpha=.25,color='C0'); ax.loglog(eps,V1,'C0-',label=r'$B$-violating decays: $\mathcal{V}_{\max}(\epsilon)$')
ax.axhline(4e7,color='C3',ls='--',label=r'pre-decay EW baryogenesis ($\mathcal{V}\lesssim 4\times10^7$)')
ax.axhline(1e15,color='k',ls=':',label=r'TeV-scale LVS ($\mathcal{V}\sim10^{15}$)')
ax.set_xlabel(r'CP-violating baryon number per decay, $\epsilon$'); ax.set_ylabel(r'$\mathcal{V}_{\max}$'); ax.set_ylim(1e3,1e17); ax.legend(fontsize=6.5,loc='lower right'); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig('fig_map.pdf'); plt.close()
for e,nm in [(1,'one'),(1e-2,'two'),(1e-4,'four'),(1e-6,'six')]:
    m('Veps'+nm,sci(Vmax_eps(e),1))
m('kYsmall',"%.2f"%kY[0]); m('kYlarge',"%.2f"%kY[-2])
# ---- axion kinetic-energy bound
fN_max=xi*np.pi*np.sqrt(gstar/15)*Tgen/Yobs
fN_LVS=np.sqrt(2)*aW*M/(4*np.pi)
m('xiv',"%.2f\\times10^{-3}"%(xi*1e3)); m('fNmax',sci(fN_max,1)); m('fNLVS',sci(fN_LVS,1)); m('fNshort',sci(fN_LVS/fN_max,1))
m('Ymaxaxion',sci(xi*np.pi*np.sqrt(gstar/15)*Tgen/fN_LVS,1))
# EMDE version of the axion bound: Y<= xi*sqrt(6)*(H_EMDE/T)*M/(f/N)*D, maximise over TR
Dtab=np.load('Ds.npy'); Df=lambda TR: np.exp(np.interp(np.log(TR),np.log(Dtab[0]),np.log(Dtab[1])))
H_E=lambda T,TR: (np.pi**2*g(T)/36)/np.sqrt(np.pi**2*g(TR)/90)*T**4/(TR**2*M)
best=max(xi*np.sqrt(6)*H_E(Tgen,TR)/Tgen*M*Df(TR) for TR in np.logspace(-3,1.5,200))
rd=xi*np.sqrt(6)*(np.sqrt(np.pi**2*gstar/90)*Tgen**2/M)/Tgen*M
m('axEMDEratio',sci(best/rd,1))
# dark radiation
for TR,nm in [(0.005,'a'),(10.,'b')]:
    fac=(43/7)*(10.75/g(TR))**(1/3)
    for dN,tag in [(0.3,'x')]:
        B=dN/fac/(1+dN/fac); m('Bamax'+nm,"%.2f"%B); m('cGmin'+nm,"%.1f"%((1/12)/B))
open('numbers2.tex','w').write("\n".join(mac)+"\n"); print(open('numbers2.tex').read())
print("V1:",np.round(np.log10(V1),1))
