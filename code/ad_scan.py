import numpy as np, json, time, sys
from ad_sim import simulate
t0=time.time()
aS=[0.1,0.3,1.0,2.0,3.0,5.0]; alS=[np.pi/8,np.pi/4,3*np.pi/8,np.pi/2,3*np.pi/4,np.pi]; cS=[0.3,1.0,3.0]
res=[]
best=(0,None)
for c in cS:
    for a in aS:
        for al in alS:
            try:
                t,eps=simulate(6,c=c,aA=a,alpha=al)
                late=eps[-300:]
                stable=late.std()<1e-3*max(1e-9,abs(late.mean()))+1e-4
                e=float(late.mean())
                res.append((c,a,float(al),e,bool(stable)))
                if stable and abs(e)>abs(best[0]): best=(e,(c,a,float(al)))
            except Exception as ex:
                res.append((c,a,float(al),None,False))
    print("c=%g done, elapsed %.0fs, best so far %s"%(c,time.time()-t0,best),flush=True)
json.dump(res,open('ad_scan.json','w'))
ok=[r for r in res if r[4] and r[3] is not None]
print("n stable:",len(ok),"of",len(res)); print("max |eps| =",max(abs(r[3]) for r in ok),"at",best[1]); 
print("percentiles |eps|:",np.percentile([abs(r[3]) for r in ok],[50,90,100]))
