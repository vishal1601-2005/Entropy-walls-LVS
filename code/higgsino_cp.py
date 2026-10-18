"""CP asymmetry of a singlet fermion (hidden bino / neutral higgsino) decaying to t s b through the RPV operator
U^c D^c D^c, with a gluino in the loop.  Uses the effective-operator rates of Kane & Winkler, arXiv:1909.04705, Eqs. (39)-(42),
generalised from one stop mass eigenstate to the full stop propagator matrix D = (M^2_stop)^-1 in the (t_R, t_L) basis."""
import numpy as np
pi=np.pi
f1=lambda x:(1-x)**5
f2=lambda x:1-8*x+8*x**3-x**4-12*x**2*np.log(x)
f3=lambda x:1+9*x-9*x**2-x**3+6*(x+x**2)*np.log(x)
g3sq=1.25
Ca,Cb,Cc=6.,4.,8.          # singlet decaying fermion, gluino in loop (KW Table 2, last row)
def kw(m1,m2,mt,c,s,k1R,k1L,k2R,k2L,lam,sin2=1.0):
    """Original Kane-Winkler Eqs. (39)-(41) for a single stop eigenstate t1 = c tR + s tL."""
    x=m2**2/m1**2
    Gt=Ca*lam**2*m1**5/(768*pi**3*mt**4)*(k1L**2*s**2*c**2+k1R**2*c**4)
    Gr=Cb*m1**5/(512*pi**3*mt**4)*((8/3)*(m2/m1)*c**2*s**2*k1L*k1R*k2L*k2R*f3(x)
        +(k1R**2*k2R**2*c**4+k1L**2*k2L**2*s**4+(2/3)*c**2*s**2*(k1L**2*k2R**2+k1R**2*k2L**2))*f2(x))
    Gd=Cc*lam**2*m1**5/(768*pi**3*mt**4)*((m1**2/(20*pi*mt**2))*k1R*k1L*k2R*k2L*c**4*s**2*f1(x)
        +(m1*m2/(16*pi*mt**2))*(k1R**2*k2R**2*c**6+k1L**2*k2L**2*c**2*s**4)*f2(x))*sin2
    return Gd/(Gt+Gr),Gd,Gt,Gr
def dmat(mR,mL,mtX):
    M2=np.array([[mR**2,mtX],[mtX,mL**2]]); D=np.linalg.inv(M2)
    return D[0,0],D[0,1],D[1,1]      # D_RR, D_RL, D_LL
def higgsino(m1,m2,mR,mL,mtX,lam,yt=0.9,a2=0.5,f3opt=0,f1opt=0,sin2=1.0):
    """Neutral higgsino mass eigenstate: couples via y_t to (psi_Uc, phi_L) and (psi_Q, phi_R*) -- i.e. cross-multiplet.
    f1opt selects the D-combination of the (unverified) f1 term: 0 -> D_RR D_LR^2, 1 -> D_RR^2 D_LL, 2 -> drop f1 term.
    f3opt selects the D-combination of the f3 term: 0 -> D_RR D_LL, 1 -> D_LR D_RL, 2 -> drop."""
    DRR,DLR,DLL=dmat(mR,mL,mtX); x=m2**2/m1**2
    kH2=a2*yt**2; k2sq=2*g3sq
    Gt=Ca*lam**2*m1**5/(768*pi**3)*kH2*(DLR**2+DRR**2)
    D3=[DRR*DLL,DLR*DLR,0.0][f3opt]
    Gr=Cb*m1**5/(512*pi**3)*kH2*(k2sq*((2*DLR**2)*f2(x)+(2/3)*(DLL**2+DRR**2)*f2(x))-(8/3)*(m2/m1)*k2sq*D3*f3(x))
    P2=kH2*k2sq*DRR*2*DLR**2
    Dc=[DRR*DLR**2,DRR**2*DLL,0.0][f1opt]
    P1=-kH2*k2sq*Dc
    Gd=Cc*lam**2*m1**5/(768*pi**3)*((m1**2/(20*pi))*P1*f1(x)+(m1*m2/(16*pi))*P2*f2(x))*sin2
    return Gd/(Gt+Gr),Gd,Gt,Gr
if __name__=="__main__":
    # --- validation against Kane-Winkler Eq. (47) (hidden bino, gluino loop, purely right-handed stop, lambda >> g3)
    m1,m2,mt=6e3,3e3,1e4; g3=np.sqrt(g3sq); k=np.sqrt(2)*g3
    e,_,_,_=kw(m1,m2,mt,1.0,0.0,2*np.sqrt(2)/3*0.36,-np.sqrt(2)/6*0.36,k,-k,lam=1e3)
    x=m2**2/m1**2; e47=g3sq*m1*m2*f2(x)/(6*pi*mt**2)
    print("KW (47) benchmark: implementation %.3e, closed form %.3e, KW quoted 1.2e-2*f2 = %.3e"%(e,e47,1.2e-2*f2(x)))
    # --- validation against Kane-Winkler Eq. (43) (gluino decaying, bino in loop), lambda >> g1
    g1=0.36; eB,_,_,_=None,None,None,None
    # single-stop limit of the generalised higgsino code must equal KW with kappa_R=kH s/c, kappa_L=kH c/s
    th=0.6; c,s=np.cos(th),np.sin(th); mtt=1e4
    kH=np.sqrt(0.5)*0.9
    ref=kw(m1,m2,mtt,c,s,kH*s/c,kH*c/s,k,-k,lam=2.0)[0]
    # one-stop propagator matrix: D = u u^T / mt^2  -> build via huge second stop mass
    DRR,DLR,DLL=c*c/mtt**2,c*s/mtt**2,s*s/mtt**2
    import types
    def hs_single(opt1,opt3):
        kH2=0.5*0.81;k2sq=2*g3sq;x=m2**2/m1**2
        Gt=Ca*4*m1**5/(768*pi**3)*kH2*(DLR**2+DRR**2)
        D3=[DRR*DLL,DLR**2,0][opt3]
        Gr=Cb*m1**5/(512*pi**3)*kH2*(k2sq*(2*DLR**2*f2(x)+(2/3)*(DLL**2+DRR**2)*f2(x))-(8/3)*(m2/m1)*k2sq*D3*f3(x))
        Dc=[DRR*DLR**2,DRR**2*DLL][opt1]
        Gd=Cc*4*m1**5/(768*pi**3)*((m1**2/(20*pi))*(-kH2*k2sq*Dc)*f1(x)+(m1*m2/(16*pi))*kH2*k2sq*DRR*2*DLR**2*f2(x))
        return Gd/(Gt+Gr)
    print("single-stop limit: KW mapped %.4e | D-matrix version (f1:0/1, f3:0/1) %.4e %.4e %.4e %.4e"%(ref,hs_single(0,0),hs_single(1,0),hs_single(0,1),hs_single(1,1)))
