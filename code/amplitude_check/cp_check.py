"""CP-invariance check of the relative sign eta between the f and fbar cut contributions.
At zero Majorana-mass phases all couplings are real and the theory is CP invariant, so S_fbar(true) must equal S_f at every kinematic point.
We test whether the RAW engine amplitudes give S_fbar = -S_f for every real coupling set (i.e. whether a universal sign eta=-1 is required)."""
import numpy as np
from engine import *
def Spoint(M,yv=0.6,ct=0.3,nth=4,nph=6):
    m1,m2=M.m1,M.m2
    pt,q,sq,bz=lorentz_q_frame(m1,yv,ct); st=np.sqrt(1-ct**2)
    b=bz(np.array([sq/2,sq/2*st,0,sq/2*ct])); sp=bz(np.array([sq/2,-sq/2*st,0,-sq/2*ct]))
    w1,wt,wb,ws=wf(np.array([m1,0,0,0]),m1),wf(pt,0),wf(b,0),wf(sp,0)
    T0f=M.T0f(w1,wt,wb,ws); T0b=M.T0fbar(w1,wt,wb,ws)
    lmag=(yv*m1**2-m2**2)/(2*sq); Cf=0; Cb=0
    xs,ws_=np.polynomial.legendre.leggauss(nth)
    for c_,wc in zip(xs,ws_):
        for ph in np.arange(nph)*2*np.pi/nph:
            sn=np.sqrt(1-c_**2); l2=lmag*np.array([sn*np.cos(ph),sn*np.sin(ph),c_])
            L2=bz(np.array([lmag,*l2])); L1=bz(np.array([np.sqrt(m2**2+lmag**2),*(-l2)]))
            wx,wl2=wf(L1,m2),wf(L2,0); g=wc/2/nph
            Cf=Cf+g*np.einsum('abxz,xzef->abef',M.Tdec(w1,wt,wx,wl2),M.Tresc_f(wx,wl2,wb,ws))
            Cb=Cb+g*np.einsum('atxz,xtef->azef',M.Tdec(w1,wl2,wx,wf(pt,0)),M.Tresc_fbar(wx,wl2,wb,ws))
    return np.sum(np.conj(T0f)*Cf),np.sum(np.conj(T0b)*Cb)
cases={
 "pure R, c=1":dict(A=[[0,0],[0,0]],B=[[0.7,0],[1.3,0]],D=[[1/9,0],[0,0]]),
 "R+L mixed, real D":dict(A=[[0,-0.4],[0,-1.3]],B=[[0.7,0],[1.3,0]],D=np.outer([np.cos(.7),np.sin(.7)],[np.cos(.7),np.sin(.7)])/9),
 "higgsino A_R=B_L":dict(A=[[-0.6,0],[0,-1.6]],B=[[0,-0.6],[1.6,0]],D=np.linalg.inv([[4,.8],[.8,4]])),
 "bino-higgsino mix":dict(A=[[-0.5,-0.1],[0,-1.6]],B=[[0.2,-0.5],[1.6,0]],D=np.linalg.inv([[4,.8],[.8,4]])),
}
for name,c in cases.items():
    M=Model(1.0,0.4,c['A'],c['B'],c['D'],2.0,(0,0))
    Sf,Sb=Spoint(M)
    print("%-22s S_f=%+.4e%+.1ei   S_fbar(raw)=%+.4e%+.1ei   ratio S_fbar/S_f = %+.4f"%(name,Sf.real,Sf.imag,Sb.real,Sb.imag,(Sb/Sf).real))
