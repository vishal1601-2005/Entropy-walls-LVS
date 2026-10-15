import numpy as np, json, warnings
warnings.filterwarnings("ignore")
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open('emde.py').read().split("def run(")[0])
M=2.435e18; Yobs=8.6e-11
rows=sorted([r for r in json.load(open('lvs_explicit.json')) if r['V']>1e3],key=lambda r:r['V'])
lV=np.log([r['V'] for r in rows]); lm=np.log([r['m_light'] for r in rows])
m_swiss=lambda V: np.exp(np.interp(np.log(V),lV,lm))
m_k3=lambda V: 1.446*M*V**(-5/3)
W=json.load(open('washout.json')); TRw=np.array([w[0] for w in W]); etaB=np.array([w[1] for w in W]); etaL=np.abs([w[3] for w in W])
def eta(TR,kind):
    if TR>=TRw[-1]: return 28/79 if kind=='B' else 28/79
    arr=etaB if kind=='B' else etaL
    return float(np.exp(np.interp(np.log(max(TR,TRw[0])),np.log(TRw),np.log(arr))))
# kappa_Y ~ 0.75 (sudden-decay value); numerical 0.7-0.9
kY=0.75
def TR_of(mm,cG=1.0):
    TR=0.2
    for _ in range(60): TR=(90/(np.pi**2*g(max(TR,1e-4))))**0.25*np.sqrt(cG*mm**3/(4*np.pi*M))
    return TR
def phi0sq(mAD,n,c=1.0,lam=1.0):
    return (np.sqrt(c/(n-1))*mAD*M**(n-3)/lam)**(2/(n-2))
def Y_B(V,mfun,n,kind,eps=1.0,ks=0.1,ri=1/6,c=1.0,lam=1.0):
    mm=mfun(V); m32=M/V; mAD=ks*m32; TR=TR_of(mm)
    p2=phi0sq(mAD,n,c,lam); nB=eps*mAD*p2
    if mAD>mm:   # case I: AD in radiation domination before modulus oscillates, then diluted
        gA=106.75; TAD=(90*mAD**2*M**2/(np.pi**2*gA))**0.25; sAD=(2*np.pi**2/45)*gA*TAD**3
        YAD=nB/sAD
        Ti=(90*mm**2*M**2*(1-ri)/(np.pi**2*g(100.)))**0.25 if False else (90*mm**2*M**2*(1-ri)/(np.pi**2*106.75))**0.25
        D=(4/3)*kY*(TR/Ti)*(1-ri)/ri
        Y=YAD*min(D,1.0)
    else:        # case II: AD during modulus domination
        Y=eps*kY*TR*p2/(3*mAD*M**2)
    return eta(TR,kind)*Y, TR, mAD>mm
def Vmax(mfun,n,kind,eps,ks=0.1,lo=3.4,hi=17):
    f=lambda lv: np.log(Y_B(10**lv,mfun,n,kind,eps,ks)[0]/Yobs)
    if f(lo)<0: return np.nan
    if f(hi)>0: return 10**hi
    return 10**brentq(f,lo,hi)
out={}
for nm,mf in [('Swiss',m_swiss),('K3',m_k3)]:
    for lab,n,kind in [('udd (n=6, B)',6,'B'),('LH_u (n=4, L)',4,'L')]:
        for eps in (1.0,0.1):
            v=Vmax(mf,n,kind,eps); out[(nm,lab,eps)]=v
            print(nm,lab,"eps=%g"%eps,"Vmax=%.2e"%v)
# print Y_B at a few volumes for sanity
for V in [1e4,1e5,1e6,1e7,1e8]:
    y,TR,c1=Y_B(V,m_swiss,6,'B'); y4,_,_=Y_B(V,m_swiss,4,'L')
    print("V=%.0e Swiss m_mod=%.2e m_AD=%.2e TR=%.2e caseI=%s  Y_B(udd)=%.2e  Y_B(LHu)=%.2e"%(V,m_swiss(V),0.1*M/V,TR,c1,y,y4))
