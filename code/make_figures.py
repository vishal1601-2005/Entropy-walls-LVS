"""Regenerates the figures that are produced from formulas and tables in the paper:
 fig_eff.png    conversion efficiencies (closed forms + quoted numerical values)
 fig_eps.png    required CP-violating baryon number per decay in sequestered LVS (uses Table 2 vacua)
 fig_admap.png  Affleck-Dine plane (volume, operator scale)   [needs ad_scan_Mstar.py]
fig_dilution.png (Fig. 1) comes from the numerical solution of the Boltzmann system and is not regenerated here."""
import numpy as np, matplotlib, os
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'mathtext.fontset':'cm','font.family':'serif'})
out=os.path.join(os.path.dirname(__file__),'..','paper','figures')
# ---- Table 2 (explicit Swiss-cheese vacua)
V=np.array([1.0e16,4.3e14,8.9e12,4.1e11,7.9e9,4.2e8,7.1e6,4.8e5,1.8e4,2.5e3])
m=np.array([2.4e-6,2.9e-4,1.0e-1,1.1e1,4.4e3,3.8e5,2.0e8,1.2e10,2.0e12,4.0e13])
TR=np.array([6.4e-19,8.6e-16,5.6e-12,6.2e-9,5.0e-5,4.1e-2,2.7e2,1.3e5,2.7e8,2.5e10])
lV=np.log(V[::-1]); f=lambda y:(lambda v:np.exp(np.interp(np.log(v),lV,np.log(y[::-1]))))
mf,Tf=f(m),f(TR); Yobs=8.6e-11
# ---- fig_eff
TRs=np.array([0.018,0.056,1.0,10.0,100.0]); Ts=np.array([156,153,145,139,131.7]); Tstar=lambda T:np.interp(np.log(T),np.log(TRs),Ts)
Tg=np.array([1e-4,0.05,0.2,0.5,1,3,10,100,170]); gg=np.array([10.75,10.75,50,65,75,80,86,96,106.75]); gR=lambda T:np.interp(np.log(T),np.log(Tg),gg)
T=np.logspace(np.log10(3e-3),np.log10(30),200); fw=(2/3)*(2*gR(T)/(5*106.75))*(T/Tstar(T))**4
etaL=28/79*fw; etaB=1-51/79*fw
fig,ax=plt.subplots(figsize=(5.6,3.9))
ax.plot(T,etaB,'C0',lw=1.8,label=r'$B$-violating decays ($\eta_B$)'); ax.plot(T,etaL,'C1',lw=1.8,label=r'$L$-violating decays ($|\eta_L|$)')
ax.plot([0.018,1.0,100.0,300.0],[1.6e-18,1.4e-10,0.02,0.33],'s',color='C1',ms=6,mfc='none',label='numerical solution')
ax.plot([100,300],[0.96,0.40],'o',color='C0',ms=6,mfc='none')
ax.plot([30,100,300],[etaB[-1],0.96,0.40],':',color='C0',lw=1.2); ax.plot([30,100,300],[etaL[-1],0.02,0.33],':',color='C1',lw=1.2)
ax.axhline(28/79,color='gray',ls='--',lw=0.9); ax.text(4e-3,0.45,r'$28/79$',color='gray',fontsize=9)
ax.set_xscale('log');ax.set_yscale('log');ax.set_ylim(1e-23,5);ax.set_xlim(2.5e-3,500);ax.set_xlabel(r'$T_R$ [GeV]');ax.set_ylabel('conversion efficiency');ax.grid(alpha=0.3)
ax.legend(fontsize=8,loc='lower right',frameon=True); plt.tight_layout(pad=0.4); plt.savefig(os.path.join(out,'fig_eff.png'),dpi=220); plt.close()
# ---- fig_eps
er=lambda v,k=0.75: Yobs/(k*Tf(v)/mf(v))
vv=np.logspace(np.log10(1.2e7),np.log10(3e8),200)
fig,ax=plt.subplots(figsize=(5.4,3.6))
ax.axvspan(1.5e7,1.5e8,color='C1',alpha=0.18,label=r'$m_{\rm soft}=M_{\rm Pl}/\mathcal{V}^2$: 10 TeV to 0.1 TeV')
ax.fill_between(vv,[er(v,0.9) for v in vv],[er(v,0.7) for v in vv],color='C0',alpha=0.3)
ax.plot(vv,[er(v) for v in vv],'C0',lw=1.8,label=r'required $\epsilon$ ($\kappa_Y=0.75$; band $0.7$--$0.9$)')
ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(1.2e7,3e8);ax.set_ylim(5e-5,1.2e-3);ax.set_xlabel(r'$\mathcal{V}$');ax.set_ylabel(r'$\epsilon$')
ax.grid(alpha=0.3,which='both');ax.legend(fontsize=8,loc='upper left'); plt.tight_layout(pad=0.4); plt.savefig(os.path.join(out,'fig_eps.png'),dpi=220); plt.close()
# ---- fig_admap
from ad_scan_Mstar import YB,zc,zs
Vv=np.logspace(2.7,6.3,300); need=lambda e,g:(Yobs/np.array([YB(v,1.0,e,g) for v in Vv]))**(2/3)
fig,ax=plt.subplots(figsize=(5.6,4.0))
ax.fill_between(Vv,need(1.2,'swiss'),need(0.16,'swiss'),color='C0',alpha=0.2)
ax.plot(Vv,need(0.51,'swiss'),'C0',lw=1.6,label=r'$Y_B=Y_B^{\rm obs}$, Swiss-cheese ($\epsilon=0.51$)'); ax.plot(Vv,need(0.51,'k3'),'C1',lw=1.6,label='K3-fibred')
ax.plot(Vv,[zc(v) for v in Vv],'C3',lw=1.6,label=r'$|\phi_0|=M_s$'); ax.plot(Vv,[zs(v) for v in Vv],'k--',lw=1.2,label=r'$M_*=M_s$')
ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(5e2,2e6);ax.set_ylim(1e-3,1.5);ax.set_xlabel(r'$\mathcal{V}$');ax.set_ylabel(r'$\zeta=M_*/M_{\rm Pl}$');ax.grid(alpha=0.3);ax.legend(fontsize=7.5,loc='lower left')
plt.tight_layout(pad=0.4); plt.savefig(os.path.join(out,'fig_admap.png'),dpi=220); plt.close()
print("figures written to",os.path.abspath(out))
