import numpy as np
pi=np.pi
I2=np.eye(2,dtype=complex)
sx=np.array([[0,1],[1,0]],dtype=complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.array([[1,0],[0,-1]],dtype=complex)
sig=[I2,sx,sy,sz]; sigb=[I2,-sx,-sy,-sz]
Z=np.zeros((2,2),dtype=complex)
gam=[np.block([[Z,sig[m]],[sigb[m],Z]]) for m in range(4)]
g0=gam[0]
C=1j*gam[2]@gam[0]                      # charge conjugation, chiral rep
def slash(p): return p[0]*gam[0]-p[1]*gam[1]-p[2]*gam[2]-p[3]*gam[3]
def uv(p,m):
    """u(p,s), v(p,s)=C ubar^T for s=0,1 ; returns arrays (2,4)"""
    E=p[0]; ps=p[1]*sx+p[2]*sy+p[3]*sz; N=np.sqrt(2*(E+m))
    A=(E+m)*I2-ps; B=(E+m)*I2+ps
    U=[];V=[]
    for s in range(2):
        xi=I2[s]; u=np.concatenate([A@xi/N,B@xi/N]); U.append(u)
        ub=u.conj()@g0; V.append(-C@ub)          # v = -C ubar^T (DHM sign)
    return np.array(U),np.array(V)
def wf(p,m):
    """DHM 2-component wavefunctions per spin: x_a (lower), ydag^{adot}, y_a, xdag^{adot}"""
    U,V=uv(p,m)
    return dict(x=U[:,:2],yd=U[:,2:],y=V[:,:2],xd=V[:,2:])
def und(a,b): return a[1]*b[0]-a[0]*b[1]           # a^alpha b_alpha, eps^{12}=+1
def dot(a,b): return -a[1]*b[0]+a[0]*b[1]          # a_adot b^adot (upper-index dotted inputs), eps_{12}=-1
if __name__=="__main__":
    rng=np.random.default_rng(1)
    for m in (0.0,0.7):
        pv=rng.normal(size=3); E=np.sqrt(pv@pv+m*m); p=np.array([E,*pv])
        U,V=uv(p,m)
        print("m=",m,"Dirac eq u:",max(np.abs((slash(p)-m*np.eye(4))@u).max() for u in U),
              " v:",max(np.abs((slash(p)+m*np.eye(4))@v).max() for v in V))
        S1=sum(np.outer(u,u.conj()@g0) for u in U)-(slash(p)+m*np.eye(4))
        S2=sum(np.outer(v,v.conj()@g0) for v in V)-(slash(p)-m*np.eye(4))
        print("  spin sums:",np.abs(S1).max(),np.abs(S2).max())
        # DHM consistency: xdag^{adot}= eps^{adot bdot} (x_b)^* ?   eps^{12}=+1 -> xdag^1=(x_2)^*, xdag^2=-(x_1)^*
        w=wf(p,m)
        for s in range(2):
            x=w['x'][s]; xd=w['xd'][s]
            pred=np.array([x[1].conj(),-x[0].conj()])
            print("  s=%d  xdag vs eps*x^*: diff %.2e (or opposite sign: %.2e)"%(s,np.abs(xd-pred).max(),np.abs(xd+pred).max()))
        # Majorana-type sum rules: sum_s x_a y^b = m delta  (y^b = eps^{bc} y_c)
        S=np.zeros((2,2),dtype=complex)
        for s in range(2):
            x=w['x'][s]; y=w['y'][s]; yup=np.array([y[1],-y[0]])
            S+=np.outer(x,yup)
        print("  sum_s x_a y^b =\n",np.round(S,6))
