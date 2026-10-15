import json, time
from k3explicit import *
def solve_generic(S):
    P=tuple(mp.mpf(S[k]) for k in ('al','ga','xi','gs','W0','a3','A3','A','B','C'))
    xh=P[2]/P[3]**mp.mpf(1.5)
    t30=(xh/(2*P[0]*P[1]))**(mp.mpf(2)/3)
    V0=(3*P[0]*P[1]/(4*P[5]*P[6]))*P[4]*mp.sqrt(t30)*mp.e**(P[5]*t30)
    t10=(4*P[7]*V0/P[8])**(mp.mpf(2)/3)
    def F(l1,lV,l3,d):
        a,b,c=mp.e**l1,mp.e**lV,mp.e**l3
        sc=mp.mpf(1)/(abs(fV(a,b,c,0,*P))+mp.mpf(10)**-300)
        return [grads[0](a,b,c,d,*P)*a*sc, grads[1](a,b,c,d,*P)*b*sc, grads[2](a,b,c,d,*P)*c*sc, fV(a,b,c,d,*P)*sc]
    a,b,c=t10,V0,t30
    d0=-(fV(a,b,c,0,*P))*b**(mp.mpf(4)/3)
    sol=mp.findroot(F,(mp.log(a),mp.log(b),mp.log(c),d0),tol=1e-28,maxsteps=120)
    a,b,c,d=mp.e**sol[0],mp.e**sol[1],mp.e**sol[2],sol[3]
    H=mp.matrix([[hess[i][j](a,b,c,d,*P) for j in range(3)] for i in range(3)])
    G=2*mp.matrix(fK(a,b,c,*P)); ev=sorted([mp.re(e) for e in mp.eig(G**-1*H)[0]])
    if min(ev)<=0: return None
    gsv=float(P[3]); NV=gsv*gsv**3/(8*mp.pi)
    return dict(V=float(b),tau1=float(a),tau3=float(c),m_fibre=float(mp.sqrt(ev[0]*NV))*M,m_volume=float(mp.sqrt(ev[1]*NV))*M,m_tau3=float(mp.sqrt(ev[2]*NV))*M,m32=float(mp.sqrt(NV)*P[4]/b)*M)
base=dict(al=.1291,ga=3.651,xi=.934,gs=.3,A3=1,A=2e-3,B=.52,C=1.9e-5)
NV=0.3*0.3**3/(8*mp.pi); W0=float(1/mp.sqrt(NV))
rows=[]
variants={'ref':{}, 'B_x3':{'B':base['B']*3}, 'B_/3':{'B':base['B']/3}, 'A_x3':{'A':base['A']*3,'C':base['C']*9}}
for vn,mod in variants.items():
    for a3 in [0.9,1.1,1.4,1.7,2.0,2.4,2.8,3.3,3.9,4.6,5.4,6.2]:
        S=dict(base); S.update(mod); S.update(W0=W0,a3=a3)
        try:
            r=solve_generic(S)
            if r is None: continue
            r['a3']=a3; r['variant']=vn; r['cf']=r['m_fibre']*r['V']**(5/3)/M; r['cm']=r['m_volume']*r['V']**1.5/M; r['m32ratio']=r['m32']*r['V']/M
            rows.append(r); print(vn,"a3=%.1f V=%.2e tau1=%.2e m_f=%.2e m_V=%.2e m_tau3=%.2e c_f=%.3f c_m=%.2f m32V/M=%.2f"%(a3,r['V'],r['tau1'],r['m_fibre'],r['m_volume'],r['m_tau3'],r['cf'],r['cm'],r['m32ratio']))
        except Exception as e: print(vn,a3,"fail",str(e)[:60])
json.dump(rows,open('k3_scan.json','w'))
