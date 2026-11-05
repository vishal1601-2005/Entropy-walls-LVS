import numpy as np
from engine import *
A=[[0,-0.4],[0,-1.3]]; B=[[0.7,0],[1.3,0]]
c,s=np.cos(.7),np.sin(.7); D=np.outer([c,s],[c,s])/9
print("rephasing invariance of Gamma_diff (R+L mixed couplings; depends only on phi2-phi1?)")
for (p1,p2) in ((0.0,0.8),(0.5,1.3),(-0.7,0.1),(1.0,1.8),(0.8,0.0),(0.8,1.6)):
    M=Model(1.0,0.5,A,B,D,2.0,(p1,p2)); print("  phi1=%+.2f phi2=%+.2f  phi2-phi1=%+.2f : Gamma_diff = %+.5e"%(p1,p2,p2-p1,M.gamma_diff(8,5)))
