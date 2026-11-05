import numpy as np
from engine import *
g3=np.sqrt(1.25); k2=np.sqrt(2)*g3; yt=0.9; g1=0.36; Ca,Cc=6.0,8.0
def run(aH,aB,sgn,m2,ms,X,phi=(0,np.pi/2),ny=8,nc=5):
    kH=yt*aH/np.sqrt(2)*np.sqrt(2)/np.sqrt(2)   # y_t * a_H with a_H the H_u^0 component of the mass eigenstate
    kH=yt*aH
    M2=np.array([[ms**2,X],[X,ms**2]]); D=np.linalg.inv(M2)
    kBR=np.sqrt(2)*g1*(2/3)*aB; kBL=-np.sqrt(2)*g1*(1/6)*aB*sgn
    A=[[-kH,-kBL if False else kBL],[0,-k2]]; B=[[kBR,-kH],[k2,0]]
    # chi1: A_R=-kH, A_L=kBL ; B_R=kBR, B_L=-kH
    A=[[-kH,kBL],[0,-k2]]; B=[[kBR,-kH],[k2,0]]
    M=Model(1.0,m2,A,B,D,2.0,phi); return M.gamma_diff(ny,nc)/M.gamma_tsb(10)*Cc/Ca
print("neutralino = higgsino (H_u^0 component 1/sqrt2) + bino admixture a_B ; m_stop=2 mu, m_t X_t=0.865 mu^2, m_gluino=0.4 mu")
for aB in (0.0,0.1,0.2,0.4):
    aH=np.sqrt(0.5*(1-aB**2))
    print("  a_B=%.1f :"%aB,["sgn=%+d eps=%+.2e"%(sg,run(aH,aB,sg,0.4,2.0,0.865)) for sg in (+1,-1)])
