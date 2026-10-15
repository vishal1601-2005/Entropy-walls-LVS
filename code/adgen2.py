exec(open('ad_lib.py').read())
def sci(x,p=1):
    e=int(np.floor(np.log10(abs(x)))); return r"%.*f\times10^{%d}"%(p,x/10**e,e)
def Vm(mf,n,kind,eps):
    f=lambda lv: np.log(Y_B(10**lv,mf,n,kind,eps)[0]/Yobs)
    lo,hi=3.4,17
    if f(lo)<0: return np.nan
    if f(hi)>0: return 10**hi
    return 10**brentq(f,lo,hi)
mac=[]
for nm,mf in [('S',m_swiss),('K',m_k3)]:
    for lab,n,kind in [('B',6,'B'),('L',4,'L')]:
        v=Vm(mf,n,kind,0.3); mac.append(r"\newcommand{\AD%s%sthree}{%s}"%(nm,lab,sci(v,1) if v==v else r"\text{none}")); print(nm,lab,v)
# required mass for eps=0.3
mac.append(r"\newcommand{\ADmreqthree}{%s}"%sci(Yobs*M/(0.133*(28/79)*0.3),1))
open('numbers7.tex','w').write("\n".join(mac)+"\n"); print(open('numbers7.tex').read())
