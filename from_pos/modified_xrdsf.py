#!/usr/bin/env python
######################################################
# Original XRDSF url: https://github.com/rraffiu/XRDSF
# Author. : Amsera
# Date.   : 21st Dec. 2020
######################################################
import time, numpy as np
from sub import load_cons, read_vasp, load_argv, oput, lpdf, dplot, calc_ints, otext



def xrd_sf(initime):
    time_rec = []
    #load global parameters
    prms = load_cons()
    if prms.ffmat == 'vasp':
        s_obj = read_vasp(prms.fname)
    time_rec.append(time.time()-initime)
    #Map K range  and G points
    oput("[Info]: Starting Map K and G...")
    kmap = np.array([x for x in prms.g_cut/s_obj.ampvec],dtype=int)
    gmap = np.array([(np.matmul(np.array([h,k,l],dtype=float),s_obj.recvec),h,k,l) for h in range(-kmap[0],kmap[0]+1) for k in range(-kmap[1],kmap[1]+1) for l in range(-kmap[2],kmap[2]+1)])
    time_rec.append(time.time()-initime)
    #Calc Structure factor
    oput("[Info]: Starting Calc Per atom...")
    ints = np.array([calc_ints(ig,s_obj,prms.dfmsd) for ig in gmap])
    ints = ints[np.argsort(ints[:,3]),:]
    oput("[Info]: Starting Calc Structure Factor...")
    gdat, indx, cnts = np.unique(ints[:,3].round(decimals=2),return_index=True,return_counts=True)
    qtmp = np.linspace(0,prms.highq,prms.reslq*prms.highq)
    ktmp = ints[indx]
    ksum = np.array([np.sum(ints[cind:cnts[iind]+cind,4]) for iind,cind in enumerate(indx)])
    ksum[ksum<=4] = 0.0 ; ksum[0]       = 0.0
    ktmp[:,4]     = ksum
    sfac = np.array([ik[4]*lpdf(qtmp,ik[3],prms.gamma) for ik in ktmp]).T
    sfac = np.array([np.sum(line) for line in sfac])
    time_rec.append(time.time()-initime)
    oput("[Info]: Processing Time : {0[0]:.4f}s {0[1]:.4f}s {0[2]:.4f}s".format(time_rec))
    #plot out S(q)
    if prms.ishow :
        oput("[Info]: Plot Structure Factor...")
        pobj = dplot(qtmp,sfac/prms.maxvl)
        if prms.verbs > 1 :
            [pobj.pinand('[{0[0]:.0f} {0[1]:.0f} {0[2]:.0f}]'.format(ikt[:3]), (ikt[3]-0.1, ikt[4]*1.01),90) for ikt in ktmp[:13,:] if ikt[4] > 20.0 ]
        pobj.show()
    else:
        adat = np.array([qtmp,sfac/prms.maxvl]).T
        oput("[Info]: Output Data in Text mode...")
        otext(adat)

if __name__=="__main__":
    initime = time.time()
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    test = xrd_sf(initime)
