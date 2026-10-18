from mpmath import mp, mpf, mpc, exp, log, sqrt, matrix, diff, findroot, pi, zeta, cos, sin, conj, inverse, eigsy, lu_solve
import sys
mp.dps=45
gs=mpf('0.3'); chi=-540
xi=-zeta(3)*chi/(2*(2*pi)**3*gs**mpf('1.5'))
gam=1/(9*sqrt(mpf(2))); W0=mpf(1); A=mpf(1)
def Yf(tb,ts): return gam*(tb**mpf('1.5')-ts**mpf('1.5'))+xi/2
def volf(tb,ts): return gam*(tb**mpf('1.5')-ts**mpf('1.5'))
def Kc(tb,ts): return -2*log(Yf(tb,ts))
def Kder(tb,ts):
    Y=Yf(tb,ts)
    Yi=[mpf('1.5')*gam*sqrt(tb), -mpf('1.5')*gam*sqrt(ts)]
    Yij=[[mpf('0.75')*gam/sqrt(tb),0],[0,-mpf('0.75')*gam/sqrt(ts)]]
    Ki=[-2*Yi[i]/Y for i in range(2)]   # dK/dtau_i
    Kij=[[ -2*Yij[i][j]/Y+2*Yi[i]*Yi[j]/Y**2 for j in range(2)] for i in range(2)]
    return Ki,Kij,Yi,Y
def pieces(tb,ts,th,a):
    Ki,Kij,Yi,Y=Kder(tb,ts)
    Kb=[Ki[0]/2,Ki[1]/2]                      # K_i (holomorphic deriv)
    Kmat=matrix(2,2)
    for i in range(2):
        for j in range(2): Kmat[i,j]=Kij[i][j]/4   # K_{i jbar}
    Kinv=inverse(Kmat)
    Ts=mpc(ts,th); eW=exp(-a*Ts)
    W=W0+A*eW; Ws=-a*A*eW
    Wi=[mpc(0),Ws]
    DW=[Wi[i]+Kb[i]*W for i in range(2)]
    return Kb,Kmat,Kinv,W,Wi,DW
def Vpot(tb,ts,th,a):
    Kb,Kmat,Kinv,W,Wi,DW=pieces(tb,ts,th,a)
    s=0
    for i in range(2):
        for j in range(2): s+=Kinv[i,j]*DW[i]*conj(DW[j])
    return (exp(Kc(tb,ts))*(s-3*abs(W)**2)).real
def Fterm(tb,ts,th,a):
    Kb,Kmat,Kinv,W,Wi,DW=pieces(tb,ts,th,a)
    eG2=exp(Kc(tb,ts)/2)*abs(W)*(conj(W)/abs(W))  # e^{G/2} with phase of Wbar? use |W|e^{K/2}
    eG2=exp(Kc(tb,ts)/2)*abs(W)
    Gbar=[conj(DW[j])/conj(W) for j in range(2)]  # G_jbar = K_j + conj(W_j/W)
    F=[sum(Kinv[i,j]*Gbar[j] for j in range(2))*eG2 for i in range(2)]
    return F

def reduced_min(a):
    # reduced 3-term potential: minimise over tau_s then V
    xi_=xi
    def Vred_of_ts(ts):
        al=8*sqrt(ts)/(3*gam)*(a*A)**2*exp(-2*a*ts); be=4*a*A*W0*ts*exp(-a*ts); de=3*xi_*W0**2/4
        d=be**2-3*al*de
        if d<0: return None
        Vv=(be-sqrt(d))/al
        return al/Vv-be/Vv**2+de/Vv**3, Vv
    best=None
    import numpy as np
    for ts in np.linspace(8,40,321):
        r=Vred_of_ts(mpf(ts))
        if r is None: continue
        if best is None or r[0]<best[0]: best=(r[0],mpf(ts),r[1])
    return best
def solve_a(a):
    a=mpf(a)
    b=reduced_min(a)
    ts0=b[1]; V0=b[2]
    tb0=((V0/gam)+ts0**mpf('1.5'))**(mpf(2)/3)
    S=V0**3
    th=pi/a
    def eqs(tb,ts):
        f1=diff(lambda x:Vpot(x,ts,th,a),tb)*S
        f2=diff(lambda x:Vpot(tb,x,th,a),ts)*S
        return [f1,f2]
    tb,ts=findroot(eqs,(tb0,ts0),tol=1e-20,maxsteps=100)
    return a,ts,tb
def solve(Vt): raise NotImplementedError
def analyze(aval):
    a,ts,tb=solve_a(aval); th=pi/a; Vt=volf(tb,ts)
    Ki,Kij,Yi,Y=Kder(tb,ts)
    G=matrix(2,2)
    for i in range(2):
        for j in range(2): G[i,j]=2*Kij[i][j]/4   # 2 K_{ij bar}
    H=matrix(2,2)
    f=lambda x,y:Vpot(x,y,th,a)
    H[0,0]=diff(f,(tb,ts),(2,0)); H[1,1]=diff(f,(tb,ts),(0,2)); H[0,1]=H[1,0]=diff(f,(tb,ts),(1,1))
    # axion theta_s mass: second derivative wrt th
    Hth=diff(lambda t:Vpot(tb,ts,t,a),th,2)
    M=inverse(G)*H
    # eigen of M (2x2)
    tr=M[0,0]+M[1,1]; det=M[0,0]*M[1,1]-M[0,1]*M[1,0]
    disc=sqrt(tr**2-4*det)
    lam=[(tr-disc)/2,(tr+disc)/2]
    vecs=[]
    for l in lam:
        v=matrix([M[0,1], l-M[0,0]])
        n=sqrt((v.T*G*v)[0,0]); vecs.append(v/n)
    v=vecs[0]; vh=vecs[1]
    m=sqrt(lam[0]); mh=sqrt(lam[1])
    Kb,Kmat,Kinv,W,Wi,DW=pieces(tb,ts,th,a)
    m32=exp(Kc(tb,ts)/2)*abs(W)
    F=Fterm(tb,ts,th,a)
    # axion coupling: g = d_phi K_bb / K_bb  (K_bb = Kij[0][0]/4)
    def Kbb(x,y): return Kder(x,y)[1][0][0]/4
    dK=v[0]*diff(Kbb,(tb,ts),(1,0))+v[1]*diff(Kbb,(tb,ts),(0,1))
    g=dK/Kbb(tb,ts)
    Gaa=g**2*m**3/(128*pi)
    # f = T_s
    Ng=12
    c=v[1]
    Ggg=Ng*m**3*c**2/(64*pi*ts**2)
    # gaugino: X=F^s
    dX=v[0]*diff(lambda x:Fterm(x,ts,th,a)[1].real,tb)+v[1]*diff(lambda y:Fterm(tb,y,th,a)[1].real,ts)
    gl=abs(dX)/(2*ts)
    Gll=Ng*gl**2*m/(16*pi)
    Mlam=abs(F[1])/(2*ts)
    # Higgs parametric: Gamma_HH+Gamma_hh = 8 z^2 Gaa-like (z=Z/Zhat, d ln Z = -d ln V^{2/3}) compute via kappa
    dlnV=v[0]*(Yi[0]/(Y-xi/2))*(Y-xi/2)/(Y-xi/2) if False else None
    return dict(V=Vt,a=a,ts=ts,tb=tb,m32=m32,m=m,mh=mh,mtheta=sqrt(Hth/G[1,1]) if False else None,
        vb=v[0],vs=v[1],Gaa=Gaa,Ggg=Ggg,Gll=Gll,Mlam=Mlam,g2=g**2,Fs=F[1],Fb=F[0])

if __name__=="__main__":
    import sys
    for aval in [1.00,1.05,1.10,1.15,1.20]:
        try:
            r=analyze(aval)
        except Exception as e:
            print("a=",aval,"fail",e); continue
        Vt=r['V']; m=r['m']
        print("a=%.2f V=%.3e ts=%.3f tb=%.1f"%(aval,float(Vt),float(r['ts']),float(r['tb'])))
        mm=float(abs(m)); 
        print("  m32=%.3e  m_light=%.3e  c_m=m*V^1.5=%.3f  m_heavy=%.3e"%(float(r['m32']),mm,float(abs(m)*Vt**1.5),float(abs(r['mh']))))
        print("  v=(%.4f,%.3e)  g^2=%.4f (8/3=2.667)"%(float(r['vb'].real),float(r['vs'].real),float(r['g2'].real)))
        Gaa=float(abs(r['Gaa']));Ggg=float(abs(r['Ggg']));Gll=float(abs(r['Gll']))
        print("  Gaa/(m^3/48pi)=%.4f  Ggg/Gaa=%.3e  Gll/Gaa=%.3e  Mlam/m=%.3e  |Fs|/m32=%.3e |Fb|/m32=%.3f"%(Gaa/(mm**3/(48*3.14159265)),Ggg/Gaa,Gll/Gaa,float(r['Mlam'])/mm,float(abs(r['Fs'])/r['m32']),float(abs(r['Fb'])/r['m32'])))
