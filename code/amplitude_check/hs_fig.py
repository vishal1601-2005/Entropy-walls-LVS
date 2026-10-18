import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from hs_mixed import run
plt.rcParams.update({'font.size':10,'mathtext.fontset':'cm','font.family':'serif'})
ab=np.array([0,0.05,0.1,0.2,0.3,0.4,0.5,0.6])
fig,ax=plt.subplots(figsize=(5.6,3.9))
ax.axhspan(1.6e-4,7e-4,color='C1',alpha=0.2,label=r'required ($\mathrm{Br}_{\tilde h\tilde h}=0.4$)')
data={}
for ms,X,st in ((2.0,0.865,'C0'),(3.0,1.298,'C3')):
    e=[abs(run(np.sqrt(0.5*(1-a**2)),a,+1,0.4,ms,X)) for a in ab]; data[ms]=e
    ax.plot(ab,np.maximum(e,1e-9),'o-',color=st,ms=4,lw=1.5,label=r'$m_{\tilde t}=%d\mu$, $X_t=2.5\,m_{\tilde t}$'%int(ms))
    print(ms,["%.2e"%v for v in e])
ax.set_yscale('log');ax.set_ylim(1e-6,3e-3);ax.set_xlabel(r'bino admixture $a_B$');ax.set_ylabel(r'$|\epsilon_{\rm CP}|$')
ax.grid(alpha=0.3,which='both');ax.legend(fontsize=8,loc='lower right');plt.tight_layout(pad=0.4)
plt.savefig('../../paper/figures/fig_higgsino.png',dpi=220)
