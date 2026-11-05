import numpy as np
from spin import *
from numpy.polynomial.legendre import leggauss
R_,L_=0,1
def boostz(p,beta):
    g=1/np.sqrt(1-beta**2); return np.array([g*(p[0]+beta*p[3]),p[1],p[2],g*(p[3]+beta*p[0])])
def lorentz_q_frame(m1,y,cth,phi=0.0):
    """chi1 at rest (mass m1). returns p_t, q, and a function boosting q-frame vectors to the chi1 frame"""
    q2=y*m1**2; sq=np.sqrt(q2); Et=(m1**2-q2)/(2*m1); Eq=m1-Et
    pt=np.array([Et,0,0,Et]); q=np.array([Eq,0,0,-Et]); beta=-Et/Eq
    return pt,q,sq,(lambda v:boostz(v,beta))
class Model:
    def __init__(s,m1,m2,A,B,D,rho=1.0,phi=(0.0,0.0)):
        s.m1,s.m2=m1,m2; s.D=np.array(D,dtype=complex); s.rho=rho; s.eta=-1
        s.A=[np.array(A[i],dtype=complex)*np.exp(-1j*phi[i]/2) for i in range(2)]
        s.B=[np.array(B[i],dtype=complex)*np.exp(+1j*phi[i]/2) for i in range(2)]
    # ---------------- tree amplitudes (arrays over spins) ----------------
    def T0f(s,w1,wt,wb,ws):
        out=np.zeros((2,2,2,2),dtype=complex)
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    for d in range(2):
                        top=sum(s.D[al,R_]*(np.conj(s.A[0][al])*dot(w1['yd'][a],wt['xd'][b])+np.conj(s.B[0][al])*und(w1['x'][a],wt['y'][b])) for al in range(2))
                        out[a,b,c,d]=s.rho*top*und(wb['y'][c],ws['y'][d])
        return out
    def T0fbar(s,w1,wtb,wbb,wsb):
        out=np.zeros((2,2,2,2),dtype=complex)
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    for d in range(2):
                        v=sum(s.D[R_,be]*(s.A[0][be]*und(w1['x'][a],wtb['y'][b])+s.B[0][be]*dot(w1['yd'][a],wtb['xd'][b])) for be in range(2))
                        out[a,b,c,d]=np.conj(s.rho)*v*dot(wbb['xd'][c],wsb['xd'][d])
        return out
    def Tdec(s,w1,wt,wx,wtb):
        out=np.zeros((2,2,2,2),dtype=complex)   # [s1,st,sx,stb]
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    for d in range(2):
                        tI=[np.conj(s.A[0][al])*dot(w1['yd'][a],wt['xd'][b])+np.conj(s.B[0][al])*und(w1['x'][a],wt['y'][b]) for al in range(2)]
                        aI=[s.A[1][be]*und(wx['y'][c],wtb['y'][d])+s.B[1][be]*dot(wx['xd'][c],wtb['xd'][d]) for be in range(2)]
                        tII=[np.conj(s.A[1][al])*dot(wx['xd'][c],wt['xd'][b])+np.conj(s.B[1][al])*und(wx['y'][c],wt['y'][b]) for al in range(2)]
                        aII=[s.A[0][be]*und(w1['x'][a],wtb['y'][d])+s.B[0][be]*dot(w1['yd'][a],wtb['xd'][d]) for be in range(2)]
                        v=0
                        for al in range(2):
                            for be in range(2):
                                v+=s.D[al,be]*(tI[al]*aI[be]-tII[al]*aII[be])
                        out[a,b,c,d]=v
        return out
    def Tresc_f(s,wx,wtb,wb,ws):       # chi2 + tbar -> b s ; array [sx,stb,sb,ss]
        out=np.zeros((2,2,2,2),dtype=complex)
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    for d in range(2):
                        v=sum(s.D[al,R_]*(np.conj(s.A[1][al])*dot(wx['yd'][a],wtb['yd'][b])+np.conj(s.B[1][al])*und(wx['x'][a],wtb['x'][b])) for al in range(2))
                        out[a,b,c,d]=s.rho*v*und(wb['y'][c],ws['y'][d])
        return out
    def Tresc_fbar(s,wx,wt,wbb,wsb):   # chi2 + t -> bbar sbar
        out=np.zeros((2,2,2,2),dtype=complex)
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    for d in range(2):
                        v=sum(s.D[R_,be]*(s.A[1][be]*und(wx['x'][a],wt['x'][b])+s.B[1][be]*dot(wx['yd'][a],wt['yd'][b])) for be in range(2))
                        out[a,b,c,d]=np.conj(s.rho)*v*dot(wbb['xd'][c],wsb['xd'][d])
        return out
    # ---------------- rates ----------------
    def gamma_tsb(s,n=10):
        m1=s.m1; xs,ws_=leggauss(n); tot=0
        for yv,wy in zip((xs+1)/2,ws_/2):
            for ct,wc in zip(xs,ws_):
                pt,q,sq,bz=lorentz_q_frame(m1,yv,ct); st=np.sqrt(1-ct**2)
                b=bz(np.array([sq/2,sq/2*st,0,sq/2*ct])); sp=bz(np.array([sq/2,-sq/2*st,0,-sq/2*ct]))
                p1=np.array([m1,0,0,0])
                w1,wt,wb,ws=wf(p1,m1),wf(pt,0),wf(b,0),wf(sp,0)
                M2=0.5*(np.abs(s.T0f(w1,wt,wb,ws))**2).sum()+0.5*(np.abs(s.T0fbar(w1,wt,wb,ws))**2).sum()
                tot+=wy*m1**2*(m1**2*(1-yv))/2*wc*M2
        return tot/((2*pi)**3*32*m1**3)
    def gamma_rpc(s,n=14):
        m1,m2=s.m1,s.m2; M=m1; xs,ws_=leggauss(n); tot=0
        s12lo,s12hi=m2**2,M**2    # t(1) chi2(2) tbar(3): s12=m(t chi2)^2
        for u,wu in zip((xs+1)/2,ws_/2):
            s12=s12lo+(s12hi-s12lo)*u
            E2=(s12-0+m2**2)/(2*np.sqrt(s12)) ; # energy of chi2 in (t,chi2) frame: (s12 - m_t^2 + m2^2)/(2 sqrt s12)
            E3=(M**2-s12)/(2*np.sqrt(s12))
            pa=np.sqrt(max(E2**2-m2**2,0)); pb=E3
            lo=(E2+E3)**2-(pa+pb)**2; hi=(E2+E3)**2-(pa-pb)**2
            for v,wv in zip((xs+1)/2,ws_/2):
                s23=lo+(hi-lo)*v
                E1=(M**2+0-s23)/(2*M); E3_=(M**2+0-s12)/(2*M); E2_=M-E1-E3_
                p1_,p3_,p2_=E1,E3_,np.sqrt(max(E2_**2-m2**2,0))
                c13=(p2_**2-p1_**2-p3_**2)/(2*p1_*p3_)
                if abs(c13)>1: continue
                P1=np.array([E1,0,0,p1_]); P3=np.array([E3_,p3_*np.sqrt(1-c13**2),0,p3_*c13]); P2=np.array([E2_,-P1[1]-P3[1],0,-P1[3]-P3[3]])
                w1=wf(np.array([m1,0,0,0]),m1); wt=wf(P1,0); wx=wf(P2,m2); wtb=wf(P3,0)
                M2=0.5*(np.abs(s.Tdec(w1,wt,wx,wtb))**2).sum()
                tot+=wu*wv*(s12hi-s12lo)*(hi-lo)*M2
        return tot/((2*pi)**3*32*m1**3)
    def gamma_diff(s,ny=10,nc=6,nth=6,nph=8):
        m1,m2=s.m1,s.m2; x=m2**2/m1**2; xs,ws_=leggauss(ny); xc,wc_=leggauss(nc); xt,wt_=leggauss(nth)
        tot=0; ph=np.arange(nph)*2*pi/nph
        for yv,wy in zip(x+(1-x)*(xs+1)/2,(1-x)*ws_/2):
            for ct,wct in zip(xc,wc_):
                pt,q,sq,bz=lorentz_q_frame(m1,yv,ct); stt=np.sqrt(1-ct**2)
                b=bz(np.array([sq/2,sq/2*stt,0,sq/2*ct])); sp=bz(np.array([sq/2,-sq/2*stt,0,-sq/2*ct]))
                p1=np.array([m1,0,0,0]); w1,wt_f,wb,ws_f=wf(p1,m1),wf(pt,0),wf(b,0),wf(sp,0)
                T0f=s.T0f(w1,wt_f,wb,ws_f); T0b=s.T0fbar(w1,wt_f,wb,ws_f)
                lmag=(yv*m1**2-m2**2)/(2*sq); pref=(1-m2**2/(yv*m1**2))/(8*pi)
                Cf=np.zeros((2,2,2,2),dtype=complex); Cb=np.zeros((2,2,2,2),dtype=complex); wsum=0
                for cth_,wth in zip(xt,wt_):
                    for phi in ph:
                        sn=np.sqrt(1-cth_**2); l2=lmag*np.array([sn*np.cos(phi),sn*np.sin(phi),cth_]);
                        L2=bz(np.array([lmag,*l2])); L1=bz(np.array([np.sqrt(m2**2+lmag**2),*(-l2)]))
                        wx,wl2=wf(L1,m2),wf(L2,0)
                        wgt=wth/2/nph    # angular average
                        # f channel: top spectator = external t ; antitop internal
                        Td=s.Tdec(w1,wt_f,wx,wl2); Tr=s.Tresc_f(wx,wl2,wb,ws_f)
                        Cf+=wgt*np.einsum('aTxz,xzef->aTef',Td,Tr)
                        # fbar channel: tbar spectator = external; top internal (momentum L2), chi2 (L1)
                        wtbx=wf(pt,0)
                        Td2=s.Tdec(w1,wl2,wx,wtbx); Tr2=s.Tresc_fbar(wx,wl2,wb,ws_f)
                        Cb+=wgt*np.einsum('atxz,xtef->azef',Td2,Tr2)
                Cf*=pref; Cb*=pref
                dlt=-np.imag(np.sum(np.conj(T0f)*Cf))+s.eta*np.imag(np.sum(np.conj(T0b)*Cb))
                M2=0.5*dlt
                tot+=wy*m1**2*(m1**2*(1-yv))/2*wct*M2
        return tot/((2*pi)**3*32*m1**3)
