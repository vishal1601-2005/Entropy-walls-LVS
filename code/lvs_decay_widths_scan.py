from lvs_decay_widths import *
for aval in [1.30,1.40,1.50,1.60]:
    try:
        r=analyze(aval)
    except Exception as e:
        print("a=",aval,"fail",e); continue
    Vt=r['V']; m=abs(r['m']); tb=r['tb']; ts=r['ts']
    # d_phi ln V
    vb,vs=r['vb'].real,r['vs'].real
    Vv=volf(tb,ts)
    dlnV=(vb*mpf('1.5')*gam*sqrt(tb)-vs*mpf('1.5')*gam*sqrt(ts))/Vv
    kap=mpf(2)/3*dlnV
    ratio_H=2*kap**2*2/ (mpf(2)/3)*mpf(1)  # placeholder
    # Gamma_HH/Gamma_aa = [z^2 kap^2/(16pi)] / [1/(48pi)] = 3 z^2 kap^2 ; two channels (HH + hh)
    Gaa=abs(r['Gaa'])
    print("a=%.2f V=%.3e ts=%.3f m32=%.3e(M) m=%.3e GeV(m*2.435e18) c_m=%.3f vs=%.3f Mlam/m=%.1f Ggg/Gaa=%.4f kappa=%.4f  (Gamma_HH+hh)/Gaa=%.3f z^2"%(
        a if False else aval,float(Vt),float(ts),float(r['m32']),float(m)*2.435e18,float(m*Vt**1.5),float(vs),float(r['Mlam']/m),float(abs(r['Ggg'])/Gaa),float(kap),float(2*3*kap**2)))
