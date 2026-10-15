import numpy as np, json
from scipy.optimize import brentq
M=2.435e18; Tc=159.0; aW=0.0338
Tt=np.array([1e-4,0.05,0.2,0.5,1,3,10,100,170,1e4]); gt=np.array([10.75,10.75,50,65,75,80,86,96,106.75,106.75])
g=lambda T: np.interp(np.log(T),np.log(Tt),gt)
def Gam_over_T4(T):
    T=np.asarray(T,float)
    brok=np.exp(0.83*T-147.7)
    sym=18*aW**5
    return np.where(T<Tc,brok,sym)
def H_RD(T): return np.sqrt(np.pi**2*g(T)/90)*T**2/M
def H_EMDE(T,TR): return (np.pi**2*g(T)/36)/np.sqrt(np.pi**2*g(TR)/90)*T**4/(TR**2*M)
# calibrate c so that RD freeze-out is 131.7 GeV (D'Onofrio et al.)
Tstar=131.7
c=Gam_over_T4(Tstar)*Tstar/H_RD(Tstar)   # Gamma/T^3 = c H
print("calibrated c = %.3f"%c)
def ratio(T,TR,c=c): return Gam_over_T4(T)*T/(c*H_EMDE(T,TR))
def edges(TR,c=c):
    Tg=np.linspace(100,3000,300000); r=ratio(Tg,TR,c)
    ok=r>=1
    if not ok.any(): return None,None
    idx=np.where(ok)[0]
    return Tg[idx[0]],Tg[idx[-1]]
rows=[]
for TR in [0.001,0.002,0.003,0.005,0.01,0.018,0.056,0.1,0.3,1,3,10,30,100]:
    lo,hi=edges(TR)
    lo2,hi2=edges(TR,c*10); lo3,hi3=edges(TR,c/10)
    rows.append((TR,lo,hi,lo2,lo3))
    print("TR=%7.3f GeV : T_*=%s  T_up=%s   (c x10: %s ; c/10: %s)"%(TR,None if lo is None else "%.1f"%lo,None if hi is None else "%.0f"%hi,None if lo2 is None else "%.1f"%lo2,None if lo3 is None else "%.1f"%lo3))
# minimal TR for any equilibrium window
f=lambda lTR: (lambda e: (-1 if e[0] is None else 1))(edges(10**lTR))
lo,hi=-3.5,-0.5  # log10 TR range
for _ in range(40):
    mid=(lo+hi)/2
    if f(mid)>0: hi=mid
    else: lo=mid
print("TR,min for sphaleron equilibrium window: %.2f MeV"%(1e3*10**hi))
for cc in [c/10,c*10]:
    lo,hi=-4,-0.5
    for _ in range(40):
        mid=(lo+hi)/2
        e=edges(10**mid,cc)
        if e[0] is not None: hi=mid
        else: lo=mid
    print("  c=%.3f: TR,min=%.2f MeV"%(cc,1e3*10**hi))
json.dump({"c":c,"rows":rows},open('sph.json','w'))
