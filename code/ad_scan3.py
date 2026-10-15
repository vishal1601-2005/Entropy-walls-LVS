import numpy as np, itertools, json
from ad_sim3 import simulate3
pi=np.pi
vals=[]
for c,a,aH,d in itertools.product((0.3,1,3),(0.3,1,3),(0.3,1,3),(pi/8,pi/4,3*pi/8,pi/2,3*pi/4,pi)):
    t,e=simulate3(6,c=c,a=a,aH=aH,alpha=d); vals.append((c,a,aH,d,float(e[-300:].mean()),float(e[-300:].std())))
ab=np.array([abs(v[4]) for v in vals])
print("MD grid n=%d median=%.3f p90=%.3f max=%.3f"%(len(ab),np.median(ab),np.percentile(ab,90),ab.max()))
json.dump({'grid':vals,'median':float(np.median(ab)),'p90':float(np.percentile(ab,90)),'max':float(ab.max())},open('ad_scan3.json','w'))
