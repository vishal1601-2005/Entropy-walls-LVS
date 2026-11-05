import sympy as sp, numpy as np
from sympy import symbols, I
# f1-type coupling product (single stop, kappa notation): tree R (conj B1R) x loop (top vertex conj(B2R), antitop vertex A1L, resc conj(A2L))
phi1,phi2=symbols('phi1 phi2',real=True); k1R,k1L,k2R,k2L=symbols('k1R k1L k2R k2L',positive=True)
B1R=k1R*sp.exp(I*phi1/2); A1L=k1L*sp.exp(-I*phi1/2); B2R=k2R*sp.exp(I*phi2/2); A2L=k2L*sp.exp(-I*phi2/2)
c=sp.conjugate
Cf1=sp.simplify(c(c(B1R))*c(B2R)*A1L*c(A2L))
print("f1-type coupling product (chirality-mixed chain):",Cf1,"  -> independent of phi1, phi2 (cannot carry sin phi12)")
# R-type mass-insertion product for comparison
CRR=sp.simplify(c(c(B1R))*c(B2R)*B1R*c(B2R)); print("R-R mass-insertion product:",CRR)
