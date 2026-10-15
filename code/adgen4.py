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
        for tag,eps in [('ceil',1.23),('bench',0.51),('med',0.16)]:
            v=Vm(mf,n,kind,eps); mac.append(r"\newcommand{\AD%s%s%s}{%s}"%(nm,lab,tag,sci(v,1) if v==v else r"\text{none}")); print(nm,lab,tag,eps,"%.2e"%v if v==v else v)
mac.append(r"\newcommand{\ADmreqbench}{%s}"%sci(Yobs*M/(0.133*(28/79)*0.51),1))
mac.append(r"\newcommand{\ADmreqceil}{%s}"%sci(Yobs*M/(0.133*(28/79)*1.23),1))
open('numbers9.tex','w').write("\n".join(mac)+"\n")
