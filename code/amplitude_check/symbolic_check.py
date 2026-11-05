"""Independent analytic check of the relative sign of the right- and left-handed mass-insertion contributions.
Part A: coupling bookkeeping (sympy): phases of the tree*loop coupling products of the R- and L-chains.
Part B: spin sums from the DHM identities (sympy, explicit epsilon tensors; no explicit spinors) for the kinematic factors S_R, S_L.
Part C: the same S_R, S_L from explicit spinors of the engine."""
import sympy as sp, numpy as np
from sympy import symbols, Matrix, Rational, I, conjugate, simplify, Eq
# ---------------- Part A: couplings ----------------
yR,yL,kR2,kL2,rho = symbols('y_R y_L k_R2 k_L2 rho',positive=True)
phi1,phi2,thX = symbols('phi1 phi2 theta_X',real=True)
DRR=symbols('D_RR',positive=True); DLR=symbols('d',positive=True)*sp.exp(I*thX)      # D_LR = |d| e^{i thX};  D_RL = conj
DRL=sp.conjugate(DLR)
# higgsino chi1: A_R=-yR e^{-i phi1/2}, B_L=-yL e^{+i phi1/2};  gluino chi2: A_L=-kL2 e^{-i phi2/2}, B_R=+kR2 e^{+i phi2/2}
A1R=-yR*sp.exp(-I*phi1/2); B1L=-yL*sp.exp(I*phi1/2); A2L=-kL2*sp.exp(-I*phi2/2); B2R=kR2*sp.exp(I*phi2/2)
c=sp.conjugate
# R-chain (engine, crossed Case II, B-structures): tree: rho D_{L R} conj(B1L) ; loop: -D_{R L} conj(B2R) B1L * rho D_{R R} conj(B2R)
treeR=rho*DLR*c(B1L); loopR=-DRL*c(B2R)*B1L*rho*DRR*c(B2R)
CR=sp.simplify(c(treeR)*loopR)
# L-chain (A-structures): tree: rho D_{R R} conj(A1R) ; loop: -D_{L R} conj(A2L) A1R * rho D_{L R} conj(A2L)
treeL=rho*DRR*c(A1R); loopL=-DLR*c(A2L)*A1R*rho*DLR*c(A2L)
CL=sp.simplify(c(treeL)*loopL)
print("Part A")
print("  C_R =",sp.simplify(CR))
print("  C_L =",sp.simplify(CL))
print("  C_L - conj(C_R) =",sp.simplify(CL-c(CR)),"  (zero => the two chains carry complex-conjugate phase factors)")
print("  Im C_R + Im C_L =",sp.simplify(sp.im(sp.expand_complex(CR))+sp.im(sp.expand_complex(CL))))
# ---------------- Part B: spin sums with identities ----------------
eps_up=Matrix([[0,1],[-1,0]]); eps_dn=Matrix([[0,-1],[1,0]])     # eps^{12}=1, eps_{12}=-1 (DHM)
s0=sp.eye(2); sx=Matrix([[0,1],[1,0]]); sy=Matrix([[0,-I],[I,0]]); sz=Matrix([[1,0],[0,-1]])
def psigma(p): return p[0]*s0-p[1]*sx-p[2]*sy-p[3]*sz          # (p.sigma)_{alpha betadot}
m1,m2=symbols('m1 m2',positive=True)
pt=symbols('pt0:4',real=True); l2=symbols('l0:4',real=True)
P=psigma(pt); L=psigma(l2)
# identities: sum x_a y^b = m d_a^b ; sum x^dag_ad y^dag^bd = m d ; sum x_a x^dag_bd = (p sigma)_{a bd} = sum y_a y^dag_bd
def R_chain():
    # S_R = sum [x1dag^{ad} ytdag_{ad}] [ychi^b yt_b] [y1dag_{gd} xtbdag^{gd}] [xchi^d xtb_d]
    tot=0
    for ad in range(2):
      for gd in range(2):
        s1=0   # chi1: sum x1dag^{ad} y1dag_{gd} = eps^{ad bd} sum x1dag_bd y1dag_gd ; sum x_bd y_gd = eps_{gd dd} sum x_bd y^dd = eps_{gd bd} m1
        for bd in range(2):
            s1+=eps_up[ad,bd]*m1*eps_dn[gd,bd]
        for b in range(2):
          for d in range(2):
            st=P[b,ad]                                  # sum_t ytdag_{ad} y_{tb} = (pt sigma)_{b ad}
            sx2=m2*eps_up[d,b]                          # sum_chi2 y^b x^d = m2 eps^{d b}
            # sum_tb x_{tb d} xtbdag^{gd} = eps^{gd ed} (l2 sigma)_{d ed}
            stb=sum(eps_up[gd,ed]*L[d,ed] for ed in range(2))
            tot+=s1*st*sx2*stb
    return sp.simplify(tot)
def L_chain():
    # S_L = sum [y1_a xt^a] [xchi_{ed} xt^{ed}] [x1^b ytb_b] [ychi_{zd} ytb^{zd}]   (tree* ; loop pieces)
    tot=0
    for a in range(2):
      for b in range(2):
        # chi1: sum y1_a x1^b = eps^{b g} sum x_g y_a ; sum x_g y_a = eps_{a d} sum x_g y^d = eps_{a g} m1
        s1=sum(eps_up[b,g]*m1*eps_dn[a,g] for g in range(2))
        for ed in range(2):
          for zd in range(2):
            # top: sum xt^a xtdag^{ed} = eps^{a g} eps^{ed dd} (pt sigma)_{g dd}
            st=sum(eps_up[a,g]*eps_up[ed,dd]*P[g,dd] for g in range(2) for dd in range(2))
            # chi2: sum xchidag_{ed} ychidag_{zd} = m2 eps_{zd ed}
            sx2=m2*eps_dn[zd,ed]
            # tbar: sum ytb_b ytbdag^{zd} = eps^{zd dd} (l2 sigma)_{b dd}
            stb=sum(eps_up[zd,dd]*L[b,dd] for dd in range(2))
            tot+=s1*st*sx2*stb
    return sp.simplify(tot)
SR=R_chain(); SL=L_chain()
dotpl=sp.simplify(pt[0]*l2[0]-pt[1]*l2[1]-pt[2]*l2[2]-pt[3]*l2[3])
print("Part B")
print("  S_R / (2 m1 m2 pt.l2) =",sp.simplify(SR/(2*m1*m2*dotpl)))
print("  S_L / (2 m1 m2 pt.l2) =",sp.simplify(SL/(2*m1*m2*dotpl)))
print("  S_R - S_L =",sp.simplify(SR-SL))
# ---------------- Part C: explicit spinors ----------------
from engine import *
rng=np.random.default_rng(7)
def rnd(m): 
    pv=rng.normal(size=3); return np.array([np.sqrt(pv@pv+m*m),*pv])
p1,ptn,l1,l2n=rnd(0.9),rnd(0),rnd(0.6),rnd(0)
w1,wt,wx,wtb=wf(p1,0.9),wf(ptn,0),wf(l1,0.6),wf(l2n,0)
SRn=0;SLn=0
for a in range(2):
  for b in range(2):
    for c_ in range(2):
      for d in range(2):
        SRn+=np.conj(und(w1['x'][a],wt['y'][b]))*und(wx['y'][c_],wt['y'][b])*dot(w1['yd'][a],wtb['xd'][d])*und(wx['x'][c_],wtb['x'][d])
        SLn+=np.conj(dot(w1['yd'][a],wt['xd'][b]))*dot(wx['xd'][c_],wt['xd'][b])*und(w1['x'][a],wtb['y'][d])*dot(wx['yd'][c_],wtb['yd'][d])
dpl=ptn[0]*l2n[0]-ptn[1:]@l2n[1:]
print("Part C (explicit spinors, random kinematics, m1=0.9, m2=0.6)")
print("  S_R/(2 m1 m2 pt.l2) = %.6f   S_L/(2 m1 m2 pt.l2) = %.6f"%((SRn/(2*0.9*0.6*dpl)).real,(SLn/(2*0.9*0.6*dpl)).real))
