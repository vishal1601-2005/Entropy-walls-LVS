import numpy as np, sys
from engine import *
f1=lambda x:(1-x)**5
f2=lambda x:1-8*x+8*x**3-x**4-12*x**2*np.log(x)
f3=lambda x:1+9*x-9*x**2-x**3+6*(x+x**2)*np.log(x)
def kw(m1,m2,mt,c,s,k1R,k1L,k2R,k2L,lam=1.0):
    x=m2**2/m1**2
    Gt=lam**2*m1**5/(768*pi**3*mt**4)*(k1L**2*s**2*c**2+k1R**2*c**4)
    Gr=m1**5/(512*pi**3*mt**4)*((8/3)*(m2/m1)*c**2*s**2*k1L*k1R*k2L*k2R*f3(x)+(k1R**2*k2R**2*c**4+k1L**2*k2L**2*s**4+(2/3)*c**2*s**2*(k1L**2*k2R**2+k1R**2*k2L**2))*f2(x))
    Gd=lam**2*m1**5/(768*pi**3*mt**4)*((m1**2/(20*pi*mt**2))*k1R*k1L*k2R*k2L*c**4*s**2*f1(x)+(m1*m2/(16*pi*mt**2))*(k1R**2*k2R**2*c**6+k1L**2*k2L**2*c**2*s**4)*f2(x))
    return Gt,Gr,Gd
def build(m1,m2,mt,c,s,k1R,k1L,k2R,k2L,phi=(0,0)):
    A=[[0,k1L],[0,k2L]]; B=[[k1R,0],[k2R,0]]
    u=np.array([c,s]); D=np.outer(u,u)/mt**2
    return Model(m1,m2,A,B,D,1.0,phi)
if __name__=="__main__":
    m1,m2,mt=1.0,0.5,3.0
    cases=[(1.0,0.0,0.7,-0.4,1.3,-1.3),(np.cos(.6),np.sin(.6),0.7,-0.4,1.3,-1.3),(np.cos(.9),np.sin(.9),1.0,-1.0,1.0,-1.0)]
    for (c,s,k1R,k1L,k2R,k2L) in cases:
        M=build(m1,m2,mt,c,s,k1R,k1L,k2R,k2L)
        Gt=M.gamma_tsb(); Gr=M.gamma_rpc()
        kt,kr,kd=kw(m1,m2,mt,c,s,k1R,k1L,k2R,k2L)
        print("c=%.3f s=%.3f: Gamma_tsb engine/KW = %.5f   Gamma_rpc engine/KW = %.5f"%(c,s,Gt/kt,Gr/kr))
