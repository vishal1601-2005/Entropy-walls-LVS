import time, json, sys
from swiss3 import *
mp.mp.dps=110
def reduced(gs,N,W0v=1,Av=1):
    chi=-540; xi=-mp.zeta(3)*chi/(2*(2*mp.pi)**3); xhv=xi/mp.mpf(gs)**mp.mpf(1.5); a=2*mp.pi/N
    P=(mp.mpf(W0v),mp.mpf(-abs(Av)),a,xhv)   # A<0 gives the negative cross term at phi=0
    A=abs(Av); W0=mp.mpf(W0v)
    al=lambda t:(8/(3*gf))*mp.sqrt(t)*(a*A)**2*mp.e**(-2*a*t)
    be=lambda t:4*a*A*W0*t*mp.e**(-a*t)
    de=3*xhv*W0**2/4
    tmin=(3*xhv/(8*gf))**(mp.mpf(2)/3)
    def vstar(t):
        d=be(t)**2-3*al(t)*de
        if d<0: return None
        return (be(t)-mp.sqrt(d))/al(t)
    def f(t):
        v=vstar(t)
        if v is None: return mp.mpf(10)**50
        return al(t)/v-be(t)/v**2+de/v**3
    lo=tmin*mp.mpf('1.0001'); hi=tmin*4
    t=gold(f,lo,hi,n=80); v=vstar(t)
    return P,t,v,float(f(t)),float(tmin)
def refine(P,tsv,Vv):
    gfm=gf
    tbv=(Vv/gfm+tsv**mp.mpf(1.5))**(mp.mpf(2)/3)
    x=[mp.log(tbv),mp.log(tsv)]
    sc=abs(V(tbv,tsv,0,P)); f=lambda u,v: V(mp.e**u,mp.e**v,0,P)/sc
    for it in range(40):
        g,H=derivs(f,x,[mp.mpf('1e-15')]*2)
        dx=mp.lu_solve(H,-g); n=mp.norm(dx)
        if n>0.5: dx=dx*(0.5/n)
        x=[x[0]+dx[0],x[1]+dx[1]]
        if n<mp.mpf('1e-25'): break
    tbv,tsv=mp.e**x[0],mp.e**x[1]
    Hs=derivs(lambda a_,b_,c_: V(a_,b_,c_,P),[tbv,tsv,mp.mpf(0)],[tbv*mp.mpf('1e-15'),tsv*mp.mpf('1e-15'),mp.mpf('1e-8')])[1]
    ev=mp.eigsy(Hs)[0]
    return tbv,tsv,P,Hs,all(e>0 for e in ev),float(V(tbv,tsv,0,P))
if __name__=="__main__":
    gs=float(sys.argv[1]); N=float(sys.argv[2])
    t0=time.time(); P,t,v,fv,tmin=reduced(gs,N)
    print("reduced: tau_s=%.3f V=%.3e Vmin=%.2e (tau_min=%.2f)"%(float(t),float(v),fv,tmin),round(time.time()-t0,1))
    r=refine(P,t,v); print("refined posdef:",r[4],"Vmin=%.3e"%r[5],round(time.time()-t0,1))
    if r[4]: print(spectrum(r[:4]))
