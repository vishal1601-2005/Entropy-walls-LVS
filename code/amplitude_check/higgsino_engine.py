import numpy as np
from engine import *
g3=np.sqrt(1.25); k2=np.sqrt(2)*g3; yt=0.9; kH=yt/np.sqrt(2)
Ca,Cc=6.0,8.0
def hs(m1,m2,mR,mL,mtX,phi=(0,np.pi/2),ny=8,nc=5):
    """neutral higgsino (A_R,B_L = -kH) with gluino (A_L=-k2, B_R=+k2); stop propagator D=(M^2)^-1 in (R,L) basis; units m1=1"""
    M2=np.array([[mR**2,mtX],[mtX,mL**2]]); D=np.linalg.inv(M2)
    A=[[-kH,0],[0,-k2]]; B=[[0,-kH],[k2,0]]
    M=Model(m1,m2,A,B,D,2.0,phi)
    Gt=M.gamma_tsb(10); Gd=M.gamma_diff(ny,nc)
    return Gd,Gt,Gd/Gt*Cc/Ca,D
if __name__=="__main__":
    import sys
    m1=1.0
    print("units mu=1.  stops degenerate m_stop, m_t X_t = mtX (in mu^2)")
    for (ms,mtX) in ((2.0,0.0),(2.0,0.865),(2.0,3.0),(3.0,0.385)):
        for m2 in (0.3,0.5):
            Gd,Gt,eps,D=hs(m1,m2,ms,ms,mtX)
            x=m2**2
            print("m_stop=%.1f mtX=%.3f m2=%.1f : eps_CP(lam''->inf) = %+.3e   [D_RR=%.4f D_LR=%.5f]"%(ms,mtX,m2,eps,D[0,0],D[0,1]))
