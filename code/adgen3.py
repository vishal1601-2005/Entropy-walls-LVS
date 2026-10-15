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
    for tag,eps in [('med',0.21),('ninety',0.77),('ceil',2.1)]:
        v=Vm(mf,6,'B',eps); mac.append(r"\newcommand{\AD%sB%s}{%s}"%(nm,tag,sci(v,1))); print(nm,tag,eps,"%.2e"%v)
mac.append(r"\newcommand{\ADmreqceil}{%s}"%sci(Yobs*M/(0.133*(28/79)*2.1),1))
open('numbers8.tex','w').write("\n".join(mac)+"\n")
