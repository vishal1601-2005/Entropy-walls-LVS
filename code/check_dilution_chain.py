import numpy as np
from scipy.integrate import solve_ivp
exec(open('emde.py').read().split("def run(")[0])      # g(), T_of_rho(), M, Gamma_of_TR
M=2.435e18
Tgrid=np.logspace(-12,16,12000); rhogrid=np.pi**2/30*g(Tgrid)*Tgrid**4
T_of_rho=lambda r: np.exp(np.interp(np.log(r),np.log(rhogrid),np.log(Tgrid)))
def run_chain(m_mod,cG=1.0,ri=1/6):
    """Independent check of D = S_i/S_fin for an asymmetry created before the modulus oscillates.
       Start at H=m_mod with rho_phi=ri*3H^2M^2, rho_R=(1-ri)*3H^2M^2 (a=1); integrate (rho_phi,rho_R,N_B) until phi decays."""
    Gam=cG*m_mod**3/(4*np.pi*M**2)
    H0=m_mod; rtot=3*H0**2*M**2; rphi0=ri*rtot; rR0=(1-ri)*rtot
    def f(N,y):
        a=rphi0*np.exp(min(y[0],60)); b=rR0*np.exp(min(y[1],60)); H=np.sqrt((a+b)/3)/M
        return [-3-Gam/H,-4+Gam*a/(H*b)]
    # time to decay: rho_phi ~ a^-3 and Gamma t~1 ; integrate far enough
    Nend=14.0+ (2/3)*np.log(H0/Gam)*1.2
    N=np.linspace(0,Nend,200001)
    sol=solve_ivp(f,[0,Nend],[0.,0.],t_eval=N,method='Radau',rtol=1e-9,atol=1e-13)
    rp=rphi0*np.exp(sol.y[0]); rr=rR0*np.exp(sol.y[1]); T=T_of_rho(rr)
    S=(4/3)*rr/T*np.exp(3*N)
    Si=S[0]; ifin=np.where(rp/rr<1e-7)[0][0]; Sf=S[ifin]
    Ti=T[0]
    TR=(90/(np.pi**2*g(max(1e-3,0.5))))**0.25*0  # placeholder
    return Si/Sf, Ti, Gam
for m_mod in (4e9,8.5e9,3e10):
    D_num,Ti,Gam=run_chain(m_mod)
    # analytic: D = (4/3) kappaY (T_R/T_i) (1-ri)/ri with T_R from Gamma (g_R at T_R), kappaY=0.75
    gR=106.75; TR=(90/(np.pi**2*gR))**0.25*np.sqrt(Gam*M)
    gR=g(TR); TR=(90/(np.pi**2*gR))**0.25*np.sqrt(Gam*M)
    ri=1/6; Ti_an=(90*m_mod**2*M**2*(1-ri)/(np.pi**2*106.75))**0.25
    D_an=(4/3)*0.75*(TR/Ti_an)*(1-ri)/ri
    print("m_mod=%.1e  T_R=%.2e  D_numeric=%.3e  D_analytic=%.3e  ratio=%.3f"%(m_mod,TR,D_num,D_an,D_num/D_an))
