import json, time
from swiss4 import *
rows=[]
for gs,N in [(0.3,2.2),(0.3,2.4),(0.3,2.7),(0.3,3.0),(0.3,3.5),(0.3,4.0),(0.3,5.0),(0.3,6.0),(0.3,8.0),(0.3,10.0),(0.4,12.0),(0.4,16.0)]:
    try:
        P,t,v,fv,tmin=reduced(gs,N)
        r=refine(P,t,v)
        if not r[4]: print("not posdef",gs,N); continue
        sp=spectrum(r[:4]); sp.update(gs=gs,N=N,Vmin=r[5]); rows.append(sp)
        print("gs=%.2f N=%.1f V=%.2e ts=%.1f m_Vol=%.2e m_ts=%.2e m_rho=%.2e m32=%.2e"%(gs,N,sp['V'],sp['ts'],sp['m_light'],sp['m_heavy'],sp['m_rho'],sp['m32']))
    except Exception as e: print("fail",gs,N,e)
json.dump(rows,open('lvs_explicit.json','w'))
