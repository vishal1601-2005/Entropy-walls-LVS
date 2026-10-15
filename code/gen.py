import numpy as np, json, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
exec(open('emde.py').read().split("out=[]")[0])          # run(), g(), M ...
exec(open('sph.py').read().split("rows=[]")[0].split("M=2.435e18; Tc=159.0; aW=0.0338")[1].join(["M=2.435e18; Tc=159.0; aW=0.0338\n",""]) if False else "")
M=2.435e18; Tc=159.0; aW=0.0338
def Gam_over_T4(T):
    T=np.asarray(T,float); return np.where(T<Tc,np.exp(np.minimum(0.83*T-147.7,50)),18*aW**5)
def H_RD(T): return np.sqrt(np.pi**2*g(T)/90)*T**2/M
def H_EMDE(T,TR): return (np.pi**2*g(T)/36)/np.sqrt(np.pi**2*g(TR)/90)*T**4/(TR**2*M)
Tstar=131.7; c=Gam_over_T4(Tstar)*Tstar/H_RD(Tstar)
def edges(TR,cc=c):
    Tg=np.linspace(100,3000,200000); ok=Gam_over_T4(Tg)*Tg/(cc*H_EMDE(Tg,TR))>=1
    if not ok.any(): return None,None
    i=np.where(ok)[0]; return Tg[i[0]],Tg[i[-1]]
Tgen=131.7
# ---------- dilution grid
TRs=np.logspace(-3,2,26); Ds=[]
for TR in TRs: Ds.append(run(TR)[0])
Ds=np.array(Ds); Cs=Ds*(Tgen/TRs)**5
np.save('Ds.npy',np.vstack([TRs,Ds]))
fig,ax=plt.subplots(1,2,figsize=(9,3.4))
ax[0].loglog(TRs,Ds,'o-',label='numerical')
ax[0].loglog(TRs,(TRs/Tgen)**5,'--',label=r'$(T_R/T_{\mathrm{gen}})^5$')
ax[0].set_xlabel(r'$T_R$ [GeV]'); ax[0].set_ylabel(r'$\mathcal{D}=S_{\mathrm{gen}}/S_{\mathrm{fin}}$'); ax[0].legend(); ax[0].grid(alpha=.3)
ax[1].semilogx(TRs,Cs,'o-'); ax[1].set_xlabel(r'$T_R$ [GeV]'); ax[1].set_ylabel(r'$\mathcal{C}=\mathcal{D}\,(T_{\mathrm{gen}}/T_R)^5$'); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig('fig_dilution.pdf'); plt.close()
Dfun=lambda TR: np.exp(np.interp(np.log(TR),np.log(TRs),np.log(Ds)))
# ---------- sphaleron window
TRg=np.logspace(-2.6,2,60); lo=[];hi=[]
for TR in TRg:
    a,b=edges(TR); lo.append(np.nan if a is None else a); hi.append(np.nan if b is None else b)
lo=np.array(lo);hi=np.array(hi)
lo10=np.array([ (edges(TR,c*10)[0] or np.nan) for TR in TRg]); lod=np.array([ (edges(TR,c/10)[0] or np.nan) for TR in TRg])
fig,ax=plt.subplots(figsize=(5.2,3.6))
ax.fill_between(TRg,lo10,lod,alpha=.25,color='C0',label=r'lower edge, criterion $\times10^{\pm1}$')
ax.semilogx(TRg,lo,'C0-',label=r'$T_*(T_R)$ (lower edge)')
ax.semilogx(TRg,np.minimum(hi,3000),'C3-',label=r'$T_{\mathrm{up}}(T_R)$ (upper edge, capped at 3 TeV)')
ax.axhline(131.7,color='k',ls=':',label=r'RD: $T_*=131.7$ GeV'); ax.axhline(159,color='gray',ls='--',label=r'$T_c=159$ GeV')
ax.set_xlabel(r'$T_R$ [GeV]'); ax.set_ylabel('Temperature [GeV]'); ax.set_yscale('log'); ax.set_ylim(100,4000); ax.legend(fontsize=6.5,loc='center right'); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig('fig_window.pdf'); plt.close()
def TRmin_eq(cc):
    lo_,hi_=-4,-0.5
    for _ in range(40):
        mid=(lo_+hi_)/2
        if edges(10**mid,cc)[0] is not None: hi_=mid
        else: lo_=mid
    return 1e3*10**hi_
TReq=[TRmin_eq(c),TRmin_eq(c/10),TRmin_eq(c*10)]
# ---------- lower bound on T_R vs maximal generated asymmetry
pref=(12/79)/(2*np.pi**2*106.75/45)
Ymax=np.logspace(-9,-2,200)*1.0; Yobs=8.6e-11
from scipy.optimize import brentq
def TRmin(Yg):
    f=lambda l: np.log(Dfun(10**l)*Yg/Yobs)
    return 10**brentq(f,-3,2)
rows=[]
for mu_over_T in [1,1e-1,1e-2,1e-3,1e-4]:
    Yg=pref*mu_over_T; rows.append((mu_over_T,Yg,TRmin(Yg)))
# LVS dictionary
def TR_of_V(V,cG=1.0):
    m=M*V**-1.5
    TR=0.2
    for _ in range(50):
        TR=(90/(np.pi**2*g(max(TR,1e-4))))**0.25*np.sqrt(cG*m**3/(4*np.pi*M))
    return m,TR
lvs=[(V,M/V)+TR_of_V(V) for V in [1e5,1e6,1e7,1e8,1e9,1e10,1e12,1e14,1e15]]
def Vmax(TRtarget,cG=1.0):
    lo_,hi_=4,16
    for _ in range(60):
        mid=(lo_+hi_)/2
        if TR_of_V(10**mid,cG)[1]>TRtarget: lo_=mid
        else: hi_=mid
    return 10**lo_
# figure: T_R^min vs Ymax with secondary V axis
fig,ax=plt.subplots(figsize=(5.2,3.6))
muT=np.logspace(-6,0,100); ax.semilogx(muT,[TRmin(pref*x) for x in muT]); ax.set_yscale('log')
ax.set_xlabel(r'maximal $\mu_{B+L}/T$ at generation'); ax.set_ylabel(r'$T_R^{\min}$ [GeV]'); ax.axhline(4e-3,color='gray',ls=':',label='BBN floor'); ax.legend(); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig('fig_TRmin.pdf'); plt.close()
# slow-source bound
Nalpha=30.
HT_RD=H_RD(Tgen)/Tgen
Y_RD=pref*Nalpha*HT_RD
best=0
for TR in np.logspace(-3,1.7,120):
    best=max(best,pref*Nalpha*(H_EMDE(Tgen,TR)/Tgen)*Dfun(TR))
need=Yobs/Y_RD
# macros
mac=[]
def m(n,v): mac.append(r"\newcommand{\%s}{%s}"%(n,v))
def sci(x,p=1):
    e=int(np.floor(np.log10(abs(x)))); mant=x/10**e
    return r"%.*f\times10^{%d}"%(p,mant,e)
m('ccal',"%.2f"%c)
m('Csmall',"%.3f"%Cs[2]); m('Clarge',"%.2f"%Cs[-4])
m('Hratio',"%.1f"%(np.pi**2*106.75/36/np.sqrt(np.pi**2*10.75/90)/np.sqrt(np.pi**2*106.75/90)))
m('TReqc',"%.1f"%TReq[0]); m('TReqlo',"%.1f"%TReq[1]); m('TReqhi',"%.1f"%TReq[2])
m('YRD',sci(Y_RD)); m('YEMDEbest',sci(best)); m('shortfall',sci(Yobs/best)); m('thetaneed',sci(need/ Nalpha*Nalpha/ (1.0) ,1))
m('HTRD',sci(HT_RD)); m('HRD',sci(H_RD(Tgen)))
for i,(x,Yg,t) in enumerate(rows):
    m('muT%s'%"abcde"[i],"%g"%x); m('TRmin%s'%"abcde"[i],("%.1f"%t) if t>=1 else ("%.0f\\,\\mathrm{MeV}"%(t*1e3)))
for TRq,name in [(0.018,'a'),(0.056,'b'),(1.0,'c'),(10.0,'d')]:
    a,b=edges(TRq); m('Tstar%s'%name,"%.0f"%a); m('Tup%s'%name,("%.0f"%b if b<3000 else r"\gtrsim 3000"))
for k,(V,m32,mV,TR) in enumerate(lvs):
    pass
m('Vmaxfew',sci(Vmax(4.0),1)); m('VmaxBBN',sci(Vmax(4e-3),1)); m('Vmaxten',sci(Vmax(10.0),1)); m('VmaxfewHi',sci(Vmax(4.0,4*np.pi),1)); m('VmaxfewLo',sci(Vmax(4.0,1/(4*np.pi)),1))
m('Vfifteen',sci(lvs[-1][3] if False else lvs[-1][3],1)); m('mVfifteen',sci(lvs[-1][2],1)); m('mthreeh',sci(lvs[-1][1],1))
open('numbers.tex','w').write("\n".join(mac)+"\n")
# LVS table rows
tab="\n".join(r"$10^{%d}$ & $%s$ & $%s$ & $%s$ \\"%(round(np.log10(V)),sci(m32,1),sci(mV,1),sci(TR,1)) for V,m32,mV,TR in lvs)
open('lvstable.tex','w').write(tab+"\n")
print(open('numbers.tex').read()); print(tab)
print("Cs:",np.round(Cs,4))
print("rows",rows)
print("TReq",TReq)
