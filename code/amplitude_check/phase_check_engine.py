import numpy as np
from engine import *
g3=np.sqrt(1.25); k2=np.sqrt(2)*g3; yt=0.9; g1=0.36
aB=0.4; aH=np.sqrt(0.5*(1-aB**2)); kH=yt*aH
kBR=np.sqrt(2)*g1*(2/3)*aB; kBL=-np.sqrt(2)*g1*(1/6)*aB
A=[[-kH,kBL],[0,-k2]]; B=[[kBR,-kH],[k2,0]]
def run(phi2,th):
    mtX=0.865*np.exp(1j*th); M2=np.array([[4.0,np.conj(mtX)],[mtX,4.0]]); D=np.linalg.inv(M2)
    return Model(1.0,0.4,A,B,D,2.0,(0,phi2)).gamma_diff(8,5)
print("engine: bino-higgsino mixture (a_B=0.4); invariant Phi = phi2 - phi1 + 2 theta_X  (theta_X = phase of D_LR up to pi)")
for phi2,th,lab in ((np.pi/2,0,"Phi=pi/2 (theta=0)"),(0,np.pi/4,"Phi=pi/2 (phi2=0,theta=pi/4)"),(np.pi/4,np.pi/8,"Phi=pi/2 (mixed)"),(np.pi/2,np.pi/4,"Phi=pi"),(np.pi/2,-np.pi/4,"Phi=0")):
    print("  %-34s Gamma_diff = %+.4e"%(lab,run(phi2,th)))
