import numpy as np, json, itertools
from ad_sim2 import simulate2
pi=np.pi
out={}
def run(**kw):
    t,e=simulate2(6,hi=400.,**kw); late=e[-300:]; return float(late.mean()),float(late.std())
# slices around benchmark (c=1,a=1,aH=1,delta=pi/2)
sl={}
sl['a']=[(v,run(a=v)) for v in (0.3,1.0,3.0)]
sl['delta']=[(v,run(alpha=v)) for v in (pi/8,pi/4,pi/2,3*pi/4)]
sl['c']=[(v,run(c=v)) for v in (0.3,1.0,3.0)]
sl['aH']=[(v,run(aH=v)) for v in (0.3,1.0,3.0)]
for k,v in sl.items(): print(k,[(round(x,2),round(e[0],3)) for x,e in v])
# full grid
vals=[]
for c,a,aH,d in itertools.product((0.3,1,3),(0.3,1,3),(0.3,1,3),(pi/8,pi/4,3*pi/8,pi/2,3*pi/4,pi)):
    e,s=run(c=c,a=a,aH=aH,alpha=d); vals.append((c,a,aH,d,e,s))
ab=np.array([abs(v[4]) for v in vals]); st=np.array([v[5] for v in vals])
print("grid n=%d  median|eps|=%.3f  90th=%.3f  max=%.3f  (max std %.1e)"%(len(ab),np.median(ab),np.percentile(ab,90),ab.max(),st.max()))
top=sorted(vals,key=lambda v:-abs(v[4]))[:5]
for v in top: print("c=%.1f a=%.1f aH=%.1f delta=%.2f eps=%.3f"%(v[0],v[1],v[2],v[3],v[4]))
json.dump({'slices':{k:[(x,e[0]) for x,e in v] for k,v in sl.items()},'grid':vals,'median':float(np.median(ab)),'p90':float(np.percentile(ab,90)),'max':float(ab.max())},open('ad_scan2.json','w'))
