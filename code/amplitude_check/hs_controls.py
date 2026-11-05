import numpy as np
from engine import *
g3=np.sqrt(1.25); k2=np.sqrt(2)*g3; yt=0.9; kH=yt/np.sqrt(2); Ca,Cc=6.0,8.0
def run(AR,BL,m2,mR,mL,mtX,phi,ny=8,nc=5):
    M2=np.array([[mR**2,np.conj(mtX)],[mtX,mL**2]],dtype=complex); D=np.linalg.inv(M2)
    A=[[AR,0],[0,-k2]]; B=[[0,BL],[k2,0]]
    M=Model(1.0,m2,A,B,D,2.0,phi); Gt=M.gamma_tsb(10); Gd=M.gamma_diff(ny,nc)
    return Gd/Gt*Cc/Ca
m2=0.4; ms=2.0; X=0.865
print("control 1: break the L<->R symmetry of the higgsino couplings (phi2=pi/2, real stop mixing)")
for lab,(AR,BL) in {"A_R=-kH, B_L=-kH (physical higgsino)":(-kH,-kH),"A_R only":(-kH,0),"B_L only":(0,-kH),"A_R=-kH, B_L=-0.5kH":(-kH,-0.5*kH)}.items():
    print("   %-40s eps = %+.3e"%(lab,run(AR,BL,m2,ms,ms,X,(0,np.pi/2))))
print("control 2: complex stop mixing phase (physical higgsino); scan arg(X_t) and gluino mass phase")
best=0
for th in (0,np.pi/4,np.pi/2,3*np.pi/4,np.pi):
    row=[]
    for ph in (0,np.pi/4,np.pi/2,np.pi):
        e=run(-kH,-kH,m2,ms,ms,X*np.exp(1j*th),(0,ph)); row.append(e); best=max(best,abs(e))
    print("   arg X_t=%.2f :"%th,["%+.1e"%v for v in row])
print("max |eps| over this phase scan: %.2e"%best)
