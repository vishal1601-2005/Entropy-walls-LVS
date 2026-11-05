import numpy as np, time
from engine import *
g3=np.sqrt(1.25); k2=np.sqrt(2)*g3; yt=0.9; g1=0.36; Ca,Cc=6.0,8.0
aB=0.4; aH=np.sqrt(0.5*(1-aB**2)); kH=yt*aH
kBR=np.sqrt(2)*g1*(2/3)*aB; kBL=-np.sqrt(2)*g1*(1/6)*aB
A=[[-kH,kBL],[0,-k2]]; B=[[kBR,-kH],[k2,0]]
def model(phi2,th,mtXabs=0.865):
    mtX=mtXabs*np.exp(1j*th); D=np.linalg.inv(np.array([[4.0,np.conj(mtX)],[mtX,4.0]]))
    return Model(1.0,0.4,A,B,D,2.0,(0,phi2))
M0=model(np.pi/2,0.0); Gt=M0.gamma_tsb(10)
t0=time.time()
ref=M0.gamma_diff(8,5); ref2=M0.gamma_diff(6,4)
print("accuracy: (8,5) %.4e vs (6,4) %.4e"%(ref,ref2))
grid=np.linspace(0,np.pi,6,endpoint=False)
best=0;arg=None;tab=[]
for p in grid:
    row=[]
    for t in grid:
        e=Cc/Ca*model(p,t).gamma_diff(6,4)/Gt; row.append(e)
        if abs(e)>abs(best): best=e;arg=(p,t)
    tab.append(row)
print("eps_CP (a_B=0.4, m_stop=2mu, lambda''->inf) over phi2 (rows) and theta_X (cols), grid of pi/6:")
for p,row in zip(grid,tab): print("  phi2=%.2f :"%p,["%+.1e"%v for v in row])
print("benchmark (phi2=pi/2, theta=0): %+.3e ; maximum over grid: %+.3e at (phi2,theta)=(%.2f,%.2f) ; time %.0fs"%(Cc/Ca*ref/Gt,best,arg[0],arg[1],time.time()-t0))
