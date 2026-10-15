exec(open('ad_lib.py').read())
def sci(x,p=1):
    e=int(np.floor(np.log10(abs(x)))); return r"%.*f\times10^{%d}"%(p,x/10**e,e)
mac=[]
def m_(n,v): mac.append(r"\newcommand{\%s}{%s}"%(n,v))
res={}
for nm,mf in [('S',m_swiss),('K',m_k3)]:
    for lab,n,kind in [('B',6,'B'),('L',4,'L')]:
        for eps,en in [(1.0,'one'),(0.1,'tenth')]:
            v=Vmax(mf,n,kind,eps); res[(nm,lab,en)]=v
            m_('AD'+nm+lab+en, sci(v,1) if v==v else r"\text{none}")
# details at Swiss udd eps=1 Vmax
V0=res[('S','B','one')]; mm=m_swiss(V0); TR=TR_of(mm)
m_('ADmmod',sci(mm,1)); m_('ADmthree',sci(M/V0,1)); m_('ADTR',sci(TR,1)); m_('ADmmodreq',sci(Yobs/(0.2*(28/79)*1.0*M),1))
# required eps at a few volumes (Swiss, udd): eps_req = Yobs/Y_B(eps=1)
for V,nm in [(1e4,'a'),(1e5,'b'),(1e6,'c')]:
    y=Y_B(V,m_swiss,6,'B',1.0)[0]; m_('ADeps'+nm,sci(Yobs/y,1))
# case II: m_AD = x m_mod
def Y2(V,mfun,x,n=6,kind='B',eps=1.0):
    mm=mfun(V); mAD=x*mm; TR=TR_of(mm); p2=phi0sq(mAD,n)
    return eta(TR,kind)*eps*kY*TR*p2/(3*mAD*M**2)
def Vmax2(mfun,x,n=6,kind='B',eps=1.0):
    f=lambda lv: np.log(Y2(10**lv,mfun,x,n,kind,eps)/Yobs)
    lo,hi=3.4,17
    if f(lo)<0: return np.nan
    if f(hi)>0: return 10**hi
    return 10**brentq(f,lo,hi)
t2=[]
for x in (1e-1,1e-2,1e-4,1e-6):
    vs=Vmax2(m_swiss,x); vk=Vmax2(m_k3,x); t2.append((x,vs,vk)); print("x=%g Swiss %.2e K3 %.2e"%(x,vs,vk))
open('ad_case2_table.tex','w').write("\n".join(r"$10^{%d}$ & $%s$ & $%s$ \\"%(round(np.log10(x)),sci(vs,1),sci(vk,1)) for x,vs,vk in t2)+"\n")
# figure
Vs=np.logspace(3.5,9,160)
fig,ax=plt.subplots(figsize=(5.6,3.8))
for nm,mf,col in [('Swiss-cheese',m_swiss,'C0'),('K3-fibred',m_k3,'C1')]:
    ax.loglog(Vs,[Y_B(V,mf,6,'B',1.0)[0] for V in Vs],col+'-',label=nm+r', $udd$ ($\epsilon=1$)')
    ax.loglog(Vs,[Y_B(V,mf,6,'B',0.1)[0] for V in Vs],col+'--')
    ax.loglog(Vs,[Y_B(V,mf,4,'L',1.0)[0] for V in Vs],col+':',label=nm+r', $LH_u$ ($\epsilon=1$)')
ax.axhline(Yobs,color='k',lw=.8); ax.text(4e3,Yobs*1.6,r'$Y_B^{\mathrm{obs}}$',fontsize=7)
ax.set_xlabel(r'$\mathcal{V}$'); ax.set_ylabel(r'$Y_B^{\mathrm{AD}}$'); ax.set_ylim(1e-30,1e-5); ax.legend(fontsize=6.5,loc='lower left'); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig('fig_AD.pdf'); plt.close()
open('numbers6.tex','w').write("\n".join(mac)+"\n"); print(open('numbers6.tex').read())
