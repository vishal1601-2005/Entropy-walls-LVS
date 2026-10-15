import numpy as np, json
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
M=2.435e18
# approximate g_*(T) (rho and s identified), log-T interpolation; O(10-20%) uncertainty near QCD
Tt=np.array([1e-4,0.05,0.2,0.5,1,3,10,100,170,1e4]); gt=np.array([10.75,10.75,50,65,75,80,86,96,106.75,106.75])
g=lambda T: np.interp(np.log(T),np.log(Tt),gt)
Tgrid=np.logspace(-12,4,8000); rhogrid=np.pi**2/30*g(Tgrid)*Tgrid**4
T_of_rho=lambda r: np.exp(np.interp(np.log(r),np.log(rhogrid),np.log(Tgrid)))
def Gamma_of_TR(TR): return np.sqrt(np.pi**2*g(TR)/90)*TR**2/M
def run(TR,Tstart=3000.,Tgen=131.7):
    G=Gamma_of_TR(TR)
    rR=np.pi**2/30*g(Tstart)*Tstart**4
    H=rR/(1.2*G*M**2); rphi=3*H**2*M**2
    def f(N,y):
        a=rphi*np.exp(min(y[0],50)); b=rR*np.exp(min(y[1],50))
        H=np.sqrt((a+b)/3)/M
        return [-3-G/H, -4+G*a/(H*b)]
    Nend=(8/3)*np.log(Tstart/TR)+14
    N=np.linspace(0,Nend,200001)
    sol=solve_ivp(f,[0,Nend],[0.,0.],t_eval=N,method='Radau',rtol=1e-10,atol=1e-12)
    u,v=sol.y; rho_phi=rphi*np.exp(u); rho_R=rR*np.exp(v)
    T=T_of_rho(rho_R); S=(4/3)*rho_R/T*np.exp(3*N)
    # first crossing of Tgen from above
    idx=np.where(T<=Tgen)[0][0]
    ifin=np.where(rho_phi/rho_R<1e-7)[0][0]
    Sgen=S[idx]; Sfin=S[ifin]
    # check phi decayed
    
    Hn=np.sqrt((rho_phi+rho_R)/3)/M
    return Sgen/Sfin, Hn[idx], T[idx]
out=[]
for TR in [0.003,0.01,0.018,0.056,0.2,1,3,10,30]:
    D,Hg,Tg=run(TR)
    C=D*(131.7/TR)**5
    Hanalytic=(np.pi**2*g(131.7)/36)/np.sqrt(np.pi**2*g(TR)/90)*131.7**4/(TR**2*M)
    out.append((TR,D,C,Hg,Hanalytic))
    print("TR=%7.3f GeV  D=%.3e  C=D(Tgen/TR)^5=%.3f  H_num(Tgen)=%.3e  H_steady=%.3e"%(TR,D,C,Hg,Hanalytic))
json.dump(out,open('dilution.json','w'))
