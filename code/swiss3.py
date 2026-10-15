import sympy as sp, mpmath as mp, sys, time, json
mp.mp.dps=90
tb,ts,xh=sp.symbols('tb ts xh',positive=True)
gam=1/(9*sp.sqrt(2)); xs=(tb,ts)
K=-2*sp.log(gam*(tb**sp.Rational(3,2)-ts**sp.Rational(3,2))+xh/2)
Kij=sp.Matrix(2,2,lambda i,j: sp.diff(K,xs[i],xs[j])/4)
Ki=[sp.diff(K,x)/2 for x in xs]
fK=sp.lambdify((tb,ts,xh),K,'mpmath'); fKi=sp.lambdify((tb,ts,xh),Ki,'mpmath'); fKij=sp.lambdify((tb,ts,xh),Kij,'mpmath')
M=2.435e18; gf=1/(9*mp.sqrt(2))
def V(tbv,tsv,phi,P):
    W0,A,a,xhv=P
    e=mp.e**(-a*tsv); ph=mp.e**(-1j*phi)
    Wc=W0+A*e*ph
    Kd=fKi(tbv,tsv,xhv); Km=mp.matrix(fKij(tbv,tsv,xhv)); Kin=Km**-1
    D=[Kd[0]*Wc, -a*A*e*ph+Kd[1]*Wc]
    s=sum(Kin[i,j]*D[i]*mp.conj(D[j]) for i in range(2) for j in range(2))
    return mp.re(mp.e**fK(tbv,tsv,xhv)*(s-3*abs(Wc)**2))
def derivs(f,x,hs):
    n=len(x); g=[]; H=mp.matrix(n,n); f0=f(*x)
    for i in range(n):
        xp=list(x); xm=list(x); xp[i]+=hs[i]; xm[i]-=hs[i]
        fp,fm=f(*xp),f(*xm); g.append((fp-fm)/(2*hs[i])); H[i,i]=(fp-2*f0+fm)/hs[i]**2
    for i in range(n):
        for j in range(i+1,n):
            def sh(si,sj):
                y=list(x); y[i]+=si*hs[i]; y[j]+=sj*hs[j]; return f(*y)
            H[i,j]=H[j,i]=(sh(1,1)-sh(1,-1)-sh(-1,1)+sh(-1,-1))/(4*hs[i]*hs[j])
    return mp.matrix(g),H
def solve(gs,W0v=1,Av=1,av=2*mp.pi):
    chi=-540; xi=-mp.zeta(3)*chi/(2*(2*mp.pi)**3); xhv=xi/mp.mpf(gs)**mp.mpf(1.5)
    P=(mp.mpf(W0v),mp.mpf(Av),mp.mpf(av),xhv)
    for tsg in [mp.mpf(x) for x in (4,6,8,10,13,17,22,28,35,45)]:
        Vg=(3*P[0]*mp.sqrt(tsg)/(4*P[2]*abs(P[1])))*mp.e**(P[2]*tsg)
        if Vg>mp.mpf('1e26') or Vg<50: continue
        tbg=(Vg/gf)**(mp.mpf(2)/3)+tsg**1.5
        x=[mp.log(tbg),mp.log(tsg)]
        sc=abs(V(tbg,tsg,0,P))
        if sc==0: continue
        f=lambda u,v: V(mp.e**u,mp.e**v,0,P)/sc
        ok=False
        for it in range(60):
            g,H=derivs(f,x,[mp.mpf('1e-12')]*2)
            try: dx=mp.lu_solve(H,-g)
            except: break
            # damp
            nrm=mp.norm(dx)
            if nrm>1: dx=dx/nrm
            x=[x[0]+dx[0],x[1]+dx[1]]
            if mp.norm(dx)<mp.mpf('1e-20'): ok=True; break
        if not ok: continue
        tbv,tsv=mp.e**x[0],mp.e**x[1]
        if tbv<=tsv or tsv<1: continue
        Hs=derivs(lambda a_,b_,c_: V(a_,b_,c_,P),[tbv,tsv,mp.mpf(0)],[tbv*mp.mpf('1e-12'),tsv*mp.mpf('1e-12'),mp.mpf('1e-6')])[1]
        # require positive-definite
        try:
            ev=mp.eigsy(Hs)[0]
        except: continue
        if all(e>0 for e in ev): return tbv,tsv,P,Hs
    return None
def spectrum(res):
    tbv,tsv,P,H=res
    Kv=mp.matrix(fKij(tbv,tsv,P[3])); G=2*Kv
    Hm=mp.matrix([[H[0,0],H[0,1]],[H[1,0],H[1,1]]])
    ev=sorted([mp.re(x) for x in mp.eig(G**-1*Hm)[0]])
    mrho2=P[2]**2*H[2,2]/(2*Kv[1,1])
    Vv=gf*(tbv**1.5-tsv**1.5); m32=mp.e**(fK(tbv,tsv,P[3])/2)*P[0]
    return dict(V=float(Vv),tb=float(tbv),ts=float(tsv),m_light=float(mp.sqrt(ev[0]))*M,m_heavy=float(mp.sqrt(ev[1]))*M,
                m_rho=float(mp.sqrt(mrho2))*M,m32=float(m32)*M,Kss=float(Kv[1,1]),Kbb=float(Kv[0,0]),a=float(P[2]))
if __name__=="__main__":
    t=time.time(); r=solve(float(sys.argv[1]),av=2*mp.pi/float(sys.argv[2])); print("found",r is not None,round(time.time()-t,1))
    if r: print(spectrum(r))

def solve2(gs,N,W0v=1,Av=1):
    chi=-540; xi=-mp.zeta(3)*chi/(2*(2*mp.pi)**3); xhv=xi/mp.mpf(gs)**mp.mpf(1.5)
    a=2*mp.pi/N; P=(mp.mpf(W0v),mp.mpf(Av),a,xhv)
    ts0=(xhv/2)**(mp.mpf(2)/3)
    def tb_of(tsv,tbg):
        f=lambda lx: mp.diff(lambda u: V(mp.e**u,tsv,0,P),lx)*mp.e**(-0)  # d/dlntb
        sc=abs(V(mp.e**tbg,tsv,0,P))
        return mp.findroot(lambda lx: mp.diff(lambda u: V(mp.e**u,tsv,0,P)/sc,lx),tbg,tol=1e-30,maxsteps=80)
    best=None
    for tsg in [ts0*x for x in (1.0,0.8,1.2,0.6,1.5)]:
        Vg=3*P[0]*mp.sqrt(tsg)/(4*a*P[1])*mp.e**(a*tsg)
        tbg=mp.log((Vg/gf)**(mp.mpf(2)/3)+tsg**1.5)
        try:
            def gfun(lt):
                global lasttb
                tsv=mp.e**lt
                lx=tb_of(tsv,gfun.tb); gfun.tb=lx
                return V(mp.e**lx,tsv,0,P)
            gfun.tb=tbg
            sc0=abs(gfun(mp.log(tsg)))
            lt=mp.findroot(lambda u: mp.diff(lambda w: gfun(w)/sc0,u),mp.log(tsg),tol=1e-25,maxsteps=60)
            tsv=mp.e**lt; lx=gfun.tb; tbv=mp.e**lx
            Hs=derivs(lambda a_,b_,c_: V(a_,b_,c_,P),[tbv,tsv,mp.mpf(0)],[tbv*mp.mpf('1e-12'),tsv*mp.mpf('1e-12'),mp.mpf('1e-6')])[1]
            ev=mp.eigsy(Hs)[0]
            if all(e>0 for e in ev): return tbv,tsv,P,Hs
        except Exception as ex:
            print("fail",tsg,ex); continue
    return None

def gold(f,lo,hi,n=70):
    r=(mp.sqrt(5)-1)/2; c=hi-r*(hi-lo); d=lo+r*(hi-lo); fc=f(c); fd=f(d)
    for _ in range(n):
        if fc<fd: hi,d,fd=d,c,fc; c=hi-r*(hi-lo); fc=f(c)
        else: lo,c,fc=c,d,fd; d=lo+r*(hi-lo); fd=f(d)
    return (lo+hi)/2
def solve3(gs,N,W0v=1,Av=1):
    chi=-540; xi=-mp.zeta(3)*chi/(2*(2*mp.pi)**3); xhv=xi/mp.mpf(gs)**mp.mpf(1.5)
    a=2*mp.pi/N; P=(mp.mpf(W0v),mp.mpf(Av),a,xhv)
    ts0=(xhv/2)**(mp.mpf(2)/3)
    Vg=3*mp.sqrt(ts0)/(4*a*P[1])*mp.e**(a*ts0); lb0=mp.log((Vg/gf)**(mp.mpf(2)/3)+ts0**1.5)
    def inner(ts_):
        lb=gold(lambda u: V(mp.e**u,ts_,0,P),lb0-2.5,lb0+2.5,n=60); return lb
    # outer: minimise over ts of V(tb*(ts))
    def outer(ts_): return V(mp.e**inner(ts_),ts_,0,P)
    tsv=gold(outer,ts0*0.6,ts0*1.5,n=45)
    lb=inner(tsv); tbv=mp.e**lb
    Hs=derivs(lambda a_,b_,c_: V(a_,b_,c_,P),[tbv,tsv,mp.mpf(0)],[tbv*mp.mpf('1e-14'),tsv*mp.mpf('1e-14'),mp.mpf('1e-7')])[1]
    ev=mp.eigsy(Hs)[0]
    return (tbv,tsv,P,Hs) if all(e>0 for e in ev) else None
