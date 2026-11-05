import numpy as np
from engine import *
g3=np.sqrt(1.25); k2=np.sqrt(2)*g3; yt=0.9; g1=0.36; Ca,Cb,Cc=6.0,4.0,8.0
aB=0.4; aH=np.sqrt(0.5*(1-aB**2)); kH=yt*aH
kBR=np.sqrt(2)*g1*(2/3)*aB; kBL=-np.sqrt(2)*g1*(1/6)*aB
D=np.linalg.inv([[4.0,0.865],[0.865,4.0]])
A=[[-kH,kBL],[0,-k2]]; B=[[kBR,-kH],[k2,0]]
M=Model(1.0,0.4,A,B,D,1.0,(0,np.pi/2))        # rho=1  (KW's lambda'' corresponds to rho=2*lambda'')
Gd=M.gamma_diff(8,5); Gt=M.gamma_tsb(10); Gr=M.gamma_rpc(14)
print("rho=1: Gdiff=%.3e Gtsb=%.3e Grpc(engine)=%.3e"%(Gd,Gt,Gr))
print(" lambda''   eps_CP(engine Grpc)   eps_CP(Grpc x4, Kane-Winkler normalisation)")
for lam in (0.2,0.5,1.0,2.0,5.0,1e6):
    rho2=(2*lam)**2
    e1=Cc*rho2*Gd/(Ca*rho2*Gt+Cb*Gr); e2=Cc*rho2*Gd/(Ca*rho2*Gt+Cb*4*Gr)
    print(" %7.1f   %+.3e            %+.3e"%(lam,e1,e2))
