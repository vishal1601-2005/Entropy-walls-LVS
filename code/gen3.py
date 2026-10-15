import numpy as np
M=2.435e18; aW=0.0338; g=106.75; Tg=131.7
gam=1/(9*np.sqrt(2)); mac=[]
def sci(x,p=1):
    e=int(np.floor(np.log10(abs(x)))); return r"%.*f\times10^{%d}"%(p,x/10**e,e)
def m(n,v): mac.append(r"\newcommand{\%s}{%s}"%(n,v))
tau_s=1/aW; N=1
Fs=lambda V: M/(2*np.pi*N)*np.sqrt(3*gam/(4*V*np.sqrt(tau_s)))
m('FsCoef',sci(Fs(1.0),1))                     # F_s = coef / sqrt(V)
Fmax=4.163e10
Vthm=(Fs(1.0)/Fmax)**2; m('Vthm',sci(Vthm,1))
H=np.sqrt(np.pi**2*g/90)*Tg**2/M; m('HTstar',sci(H,1))
def Gam_ratio(V):
    m32=M/V; mrho=m32*np.log(V); F=Fs(V)
    G=(3/8)*aW**2*mrho**3/(32*np.pi**3*F**2)
    return mrho,G,G/H
rows=[]
for V in [1e6,1e8,1e10,1e12,1e15]:
    mr,G,r=Gam_ratio(V); rows.append((V,Fs(V),mr,r)); print("V=%.0e F_s=%.2e m_rho=%.1e Gamma=%.1e Gamma/H=%.1e"%(V,Fs(V),mr,G,r))
m('GHmin',sci(min(r for *_,r in rows),1))
tab="\n".join(r"$10^{%d}$ & $%s$ & $%s$ & $%s$ \\"%(round(np.log10(V)),sci(F,1),sci(mr,1),sci(r,1)) for V,F,mr,r in rows)
open('fstable.tex','w').write(tab+"\n")
open('numbers3.tex','w').write("\n".join(mac)+"\n"); print("\n".join(mac))
