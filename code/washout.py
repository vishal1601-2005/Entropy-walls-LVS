import numpy as np, json, warnings
warnings.filterwarnings("ignore")
from scipy.integrate import solve_ivp
M=2.435e18; Tc=159.0; aW=0.0338
Tt=np.array([1e-4,0.05,0.2,0.5,1,3,10,100,170,1e4]); gt=np.array([10.75,10.75,50,65,75,80,86,96,106.75,106.75])
g=lambda T: np.interp(np.log(T),np.log(Tt),gt)
Tgrid=np.logspace(-12,4,8000); rhogrid=np.pi**2/30*g(Tgrid)*Tgrid**4
T_of_rho=lambda r: np.exp(np.interp(np.log(r),np.log(rhogrid),np.log(Tgrid)))
Gam_T4=lambda T: (np.exp(min(0.83*T-147.7,50)) if T<Tc else 18*aW**5)
c_cal=0.118; kw=1/c_cal; EQ=28/79
def run(TR,epsB,epsL,Tstart=3000.):
    G=np.sqrt(np.pi**2*g(TR)/90)*TR**2/M
    rR=np.pi**2/30*g(Tstart)*Tstart**4; H0=rR/(1.2*G*M**2); rphi=3*H0**2*M**2
    def f(N,y):
        u,v,xB,xL=y
        a=rphi*np.exp(min(u,50)); b=rR*np.exp(min(v,50)); H=np.sqrt((a+b)/3)/M
        T=T_of_rho(b); lam=kw*Gam_T4(T)*T   # Gamma_sph/T^3 * kw  [GeV]
        dec=G*np.exp(min(u+3*N,60))/H        # decays per e-fold, units rho_phi0/m
        dxB=epsB*dec-(lam/H)*(xB-EQ*xL)
        dxL=(epsB-epsL)*dec
        return [-3-G/H,-4+G*a/(H*b),dxB,dxL]
    Nend=(8/3)*np.log(Tstart/TR)+14
    sol=solve_ivp(f,[0,Nend],[0,0,0,0],method='Radau',rtol=1e-9,atol=1e-35,dense_output=False)
    return sol.y[2,-1],sol.y[3,-1]
rows=[]
for TR in [0.003,0.005,0.01,0.018,0.056,0.2,1,3,10,30,100,300]:
    xb,xl=run(TR,1,0); xb2,xl2=run(TR,0,1)
    rows.append((TR,xb,xl,xb2,xl2)); print("TR=%7.3f  B-source: eta_B=%.4f (B-L=%.3f) | L-source: eta=%.3e (B-L=%.3f)"%(TR,xb,xl,xb2,xl2))
json.dump(rows,open('washout.json','w'))
