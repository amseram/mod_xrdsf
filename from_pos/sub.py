######################################################
# Original XRDSF url: https://github.com/rraffiu/XRDSF
# Author. : Amsera
# Date.   : 21st Dec. 2020
######################################################
import os, sys
import numpy as np
import pandas as pd
import argparse as arg
import matplotlib
from matplotlib import rc
import matplotlib.pyplot as plt
import xrdsflib as x
_verbs = 1

class load_cons:
    def __init__(self):
        global _verbs
        self.isint = True 
        self.argvs = load_argv()
        self.fname = self.argvs[0]
        self.ffmat = self.argvs[1]
        self.dfmsd = self.argvs[2]
        self.verbs = self.argvs[3] 
        self.ishow = self.argvs[4]
        self.highq = self.argvs[5]
        self.reslq = self.argvs[6]
        self.gamma = self.argvs[7]
        self.g_cut = self.argvs[8]
        self.maxvl = lpdf(0.0,0.0,self.gamma)
        #avoid recursion
        if self.isint :
            _verbs     = self.verbs
            self.isint = False

class read_vasp:
    def __init__(self,fname,isxdat=False):
        self.f_name = fname
        if not os.path.isfile(fname):
            exit("Can't find {}!! PLZ check it!".format(fname))
        self.load_info()
        if isxdat:
            self.load_xdat()
        else:
            self.load_pscr()

    def load_info(self):
        self.inilne = 7 
        self.atmlne = 5
        self.rawdat = open(self.f_name,'r').readlines()
        self.celvec = np.array([line.strip().split() for line in self.rawdat[2:5]],dtype=float)
        tmp = self.rawdat[self.atmlne].strip().split()[0].isalpha() 
        if not tmp: 
            self.atmlne = 4
            self.inilne = 6
        if tmp : self.atmlst = self.rawdat[self.atmlne].strip().split()
        self.atmnum = np.array([int(item) for item in self.rawdat[self.atmlne+1].strip().split()])
        self.totatm = np.sum(self.atmnum)
        self.chk_vers()
        self.atmrnk = [self.atmlst[ia] for ia in range(len(self.atmnum)) for j in range(self.atmnum[ia])]
        self.recvec = 2.*np.pi*np.linalg.inv(self.celvec).T
        self.ampvec = np.linalg.norm(self.recvec,axis=1)
        self.chk_frac()

    def chk_vers(self):
        if self.atmlne < 5:
            oput("[Warning]: VASP4.x POSCAR are using!! Load Atoms' info from 'SYSTEM' line")
            self.atmlst = self.rawdat[0].strip().split()
            if len(self.atmlst) != len(self.atmnum):
                oput("[Warning]: Load Atom's info from 'SYSTEM' line failed! Set all Atoms as 'He'!!")
                self.atmlst = ['He' for item in self.atmnum]
    
    def chk_frac(self):
        tmp = self.rawdat[self.inilne][0].upper()
        if tmp == "S" :
            self.inilne = self.inilne + 1
            tmp = self.rawdat[self.inilne][0].upper()
        if tmp == "C" :
            self.isfrac = False
        else:
            self.isfrac = True

    def load_pscr(self):
        if self.isfrac:
            self.atmloc = np.array([np.matmul(np.array(ilne.strip().split()[:3],dtype=float),self.celvec) for ilne in self.rawdat[self.inilne+1:self.inilne+1+self.totatm]])
        else:
            self.atmloc = np.array([ilne.strip().split()[:3] for ilne in self.rawdat[self.inilne+1:self.inilne+1+self.totatm]],dtype=float)
        if self.inilne > self.atmlne + 2 :
            self.select = [ilne.strip().split()[3:] for ilne in self.rawdat[self.inilne+1:self.inilne+1+self.totatm]]
        self.gen_atom_array()

    def gen_atom_array(self):
        self.atmary = [(self.atmrnk[iatm],self.atmloc[iatm]) for iatm in range(self.totatm)]
        tmd         = load_data()
        self.atmobj = dict([(ia,ff(ia,tmd)) for ia in set(self.atmlst)])

    def load_xdat(self):
        self.sepdat = dict([ (str(1+int((nline-self.inilne)/self.totatm)),self.rawdat[nline+1:nline+self.totatm]) for nline in range(len(self.rawdat)) if "configuration=" in self.rawdat[nline]])
        if self.isfrac :
            self.atmtrj = dict([tuple([str(ia),[self.atmrnk[ia],np.array([np.matmul(np.array(self.sepdat[ikey][ia].strip().split()[:3],dtype=float),self.celvec) for ikey in list(self.sepdat.keys())])]]) for ia in range(len(self.atmrnk))])
        else:
            self.atmtrj = dict([tuple([str(ia),[self.atmrnk[ia],np.array([self.sepdat[ikey][ia].strip().split()[:3] for ikey in list(self.sepdat.keys())])]]) for ia in range(len(self.atmrnk))])

        

def load_argv():
    parser = arg.ArgumentParser(description='Calculate Structure Factor from (VASP inputfile) POSCAR-like stucture files')
    parser.add_argument(dest='fname',metavar='filename', nargs='*')
    parser.add_argument('-f',      action="store",dest="fmat", metavar="format",            default='vasp',help='Define stucture file format: such as vasp')
    parser.add_argument('-m',      action="store",dest='msd',  metavar='msd_values',        default='0.1', help='Define MSD: such as 0.01')
    parser.add_argument('-vb',     action="store",dest='verb', metavar='verbers',           default='1',   help='Define verbs: select in 0 1 2')
    parser.add_argument('-is_plot',action="store",dest='ishow',metavar='plot or not',       default="T",   help='Plot data or Output as Text: T[rue] of F[alse]')
    parser.add_argument('-qmax',   action="store",dest='qmax', metavar='Max Q',             default="7",   help='Define Max Q')
    parser.add_argument('-qres',   action="store",dest='qres', metavar='Q resolution',      default="1000",help='Define Q resolution')
    parser.add_argument('-gam',    action="store",dest='gamm', metavar='Gamma',             default='0.01',help='Define Gamma')
    parser.add_argument('-gcut',   action="store",dest='gcut', metavar='G_cutoff',          default='10',  help='Define G_cutoff')
    args   = parser.parse_args()
    fname  = args.fname[0] if args.fname != [] else "POSCAR"
    if args.fname != fname :
        oput("[Info]: Undefined Structure Filename!!! Use default Filename : POSCAR!!!")
    if not os.path.isfile(fname):
        if not os.path.isfile(fname+".vasp"):
            exit("[Error]: PLZ echeck the Structure Filename!!! {} is not here!!!".format(fname))
        else:
            oput("[Info]: Complete Structure .fname with '.vasp' Filename Extention") 
            fname = fname + ".vasp"
    ishow = True if args.ishow[0].upper() == 'T' else False
    return fname, args.fmat, float(args.msd), int(args.verb), ishow, int(args.qmax), int(args.qres), float(args.gamm), float(args.gcut) 

def exit(hlt):
    sys.exit("[Error]: {}".format(hlt))


class load_data:
    def __init__(self):
        self.prmtab = pd.read_csv('trimed.csv',index_col=0)

    def chksym(self):
        try:
            self.elmtab = self.prmtab.loc[self.cursym]
        except:
            exit("[Error]: Incorrect symbol ({}) or no data available for this element.".format(self.cursym))

class ff:
    def __init__(self,sym,data_base=load_data()):
        data_base.cursym = sym
        data_base.chksym()
        self.cursym = sym
        self.prmtab = data_base.prmtab
        self.elmtab = data_base.elmtab
        self.atmprm = np.array([np.array([self.elmtab.loc['a{}'.format(ind)],self.elmtab.loc['b{}'.format(ind)]],dtype=float) for ind in range(1,5)]).T
        self.atmc   = float(self.elmtab.loc['c'])
        self.func   = np.frompyfunc(self.eaff_func,3,1)
        self.atma   = self.atmprm[0]
        self.atmb   = self.atmprm[1]

    def eaff_func(self,x,y,l):
        return x*np.exp(y*l)

    def calc_eaff(self,k):
        tmpl = -1.*np.power(k/(4.*np.pi), 2.0)
        return self.atmc+np.sum(self.func(self.atma,self.atmb,tmpl))

def gpdf(k,k0=0.0,sig=0.05):
    return np.exp(-np.power((k-k0)/sig,2)/2)

def lpdf(k,k0=0.0,sig=0.01):
    return  sig/(2*np.pi)*(1/((k-k0)**2+(0.5*sig)**2))

def oput(strs):
    global _verbs 
    if _verbs == 0 :
        return
    elif _verbs == 1 :
        if "Warn" not in strs:
            return
    print(strs)

def otext(sf):
    with open('sfact.dat','w') as fobj:
        [print("{0[0]:.6f} {0[1]:.12f}".format(isf),file=fobj) for isf in sf]
        [print("{0[0]:.6f} {0[1]:.12f}".format(isf)) for isf in sf]


class dplot:
    def __init__(self,q,sf):
        matplotlib.use('TkAgg')
        rc('text',usetex=True)
        rc('font',family='serif', size=12)
        fig = plt.figure()
        ax  = fig.add_subplot(1,1,1)
        ax.plot(q,sf,color='k')
        plt.ion()
        plt.xlim(q[0],q[-1])
        plt.ylim(0,1.2*np.max(sf))
        plt.yticks([])
        plt.xlabel(r"q (\AA$^{-1}$)")
        plt.ylabel('Intensity (arb. units)')
        plt.margins(0.05,0.1)
        self.fline = []
        self.sline = []
        self.inilog()
    
    def inilog(self):
        self.infobj = open("xrdsf.log",'w')

    def pinand(self,tag,loc,rot=0):
        plt.annotate(tag, loc,rotation=rot)
        self.fline.append('{0[0]:12.6f} '.format(loc))
        self.sline.append('{0[1]:12.6f} '.format(loc))

    def show(self):
        if len(self.fline) > 0:
            print("<"+"="*150+">")
            print("\t".join(self.fline))
            print("\t".join(self.fline),file=self.infobj)
            print("\t".join(self.sline))
            print("\t".join(self.sline),file=self.infobj)
        plt.ioff()
        plt.show()

def calc_ints(cur_g,s_obj,msd):
    nrm_g = np.linalg.norm(cur_g[0])
    s_atm = np.sum(np.array([s_obj.atmobj[ia[0]].calc_eaff(nrm_g)*np.exp(np.dot(cur_g[0],ia[1])*1j) for ia in s_obj.atmary]))
    dwf_g = np.exp(-msd*np.power(nrm_g,2.)/3.)
    #print(nrm_g,s_atm,dwf_g)
    return [cur_g[1],cur_g[2],cur_g[3],nrm_g,np.power(np.abs(s_atm*dwf_g),2.)]

def f_calc_ints(cur_g,s_obj,msd):
    atmloc = np.array(s_obj.atmloc,order='F')
    curgpt = np.array(cur_g[0],order='F')
    atmpra = np.array([[s_obj.atmobj[ia[0]].atmprm[0]] for ia in s_obj.atmary],order='F') 
    atmprb = np.array([[s_obj.atmobj[ia[0]].atmprm[1]] for ia in s_obj.atmary],order='F') 
    atmprc = np.array([[s_obj.atmobj[ia[0]].atmc] for ia in s_obj.atmary],order='F') 
    tmp    = x.f_calc_ints(msd,atmloc,curgpt,atmpra,atmprb,atmprc)
    return [cur_g[1],cur_g[2],cur_g[3],tmp[0],tmp[1]]

def f_loop_g(knds,sobj,msd):
    mh     = knds[0]
    mk     = knds[1]
    ml     = knds[2]
    ng     = (2*mh+1)*(2*mk+1)*(2*ml+1)
    rvec   = np.array(sobj.recvec,order='F')
    atmpra = np.array([sobj.atmobj[ia[0]].atma for ia in sobj.atmary],order='F') 
    atmprb = np.array([sobj.atmobj[ia[0]].atmb for ia in sobj.atmary],order='F') 
    atmprc = np.array([sobj.atmobj[ia[0]].atmc for ia in sobj.atmary],order='F') 
    atmloc = np.array(sobj.atmloc,order='F')
    tmp    = x.f_loop(msd,mh,mk,ml,ng,rvec,atmloc,atmpra,atmprb,atmprc)
    return tmp

