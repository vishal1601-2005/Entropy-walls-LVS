import sympy as sp
# exact Kahler geometry of K3-fibred LVS: V = alpha*sqrt(t1)*t2 - gamma*t3^(3/2), K=-2 ln(V)
t1,t2,t3,al,ga=sp.symbols('t1 t2 t3 alpha gamma',positive=True)
V=al*sp.sqrt(t1)*t2-ga*t3**sp.Rational(3,2)
K=-2*sp.log(V)
xs=(t1,t2,t3)
Kij=sp.Matrix(3,3,lambda i,j: sp.diff(K,xs[i],xs[j])/4)
# leading order in t3^{3/2}/V -> 0
lead=sp.simplify(Kij.subs(ga,0))
print("K_ij (gamma->0):",lead)
# flat direction at fixed V: dt2/t2 = -(1/2) dt1/t1  -> tangent vector
v=sp.Matrix([t1, -t2/2, 0])
G=sp.simplify(2*(v.T*lead*v)[0])   # G = 2 M^2 K_ij v^i v^j  in units M=1, per (d ln t1)^2
print("G_flat per (dln t1)^2 =",G)
# volume direction at fixed t1: dt2/t2 = dV/V
w=sp.Matrix([0,t2,0]); print("G_vol per (dlnV)^2 =",sp.simplify(2*(w.T*lead*w)[0]))
# axion decay constant for fibre axion: K_11 = 1/(4 t1^2)
print("K_11 =",lead[0,0]," K_22=",lead[1,1])
