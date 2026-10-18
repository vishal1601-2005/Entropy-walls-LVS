import numpy as np
from test_cp import buildN
from test_kw import f1,f2,pi
m1,mt=1.0,3.0
c,s=np.cos(.7),np.sin(.7)
print("   m2    x    |  f2-piece: engine/KW   |  f1-piece: engine/KW")
for m2 in (0.3,0.5,0.7):
    x=m2**2/m1**2
    res={}
    for sg in (+1,-1):
        k1R,k1L,k2R,k2L=0.7,-0.4,1.3,sg*1.3
        M=buildN(m1,m2,mt,c,s,k1R,k1L,k2R,k2L,phi=(0,pi/2)); res[sg]=M.gamma_diff()
    e2=(res[1]+res[-1])/2; e1=(res[1]-res[-1])/2      # f2-piece (even in k2L), f1-piece (odd)
    k1R,k1L,k2R=0.7,-0.4,1.3; k2L=1.3
    pref=m1**5/(768*pi**3*mt**4)
    kw2=pref*(m1*m2/(16*pi*mt**2))*(k1R**2*k2R**2*c**6+k1L**2*k2L**2*c**2*s**4)*f2(x)
    kw1=pref*(m1**2/(20*pi*mt**2))*k1R*k1L*k2R*k2L*c**4*s**2*f1(x)
    print("  %.2f %.3f |  %.4f                |  %.4f   (engine f1-piece %.3e, KW %.3e)"%(m2,x,e2/kw2,e1/kw1,e1,kw1))
