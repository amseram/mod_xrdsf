#!/usr/bin/env python
import os, sys
import copy as cp 
import numpy as np
from ase import units as u
from ase.atoms import Atoms
from ase.calculators.vasp import Vasp as v
from ase.io.vasp import read_vasp as rdv, write_vasp as wrv 

class sqs:
    def __init__(self, fdat="./sq.dat",qmin=1.4,qmax=10.0):
        self.edat = fdat 
        self.esqs = np.loadtxt(fdat).T 
        self.qrng = self.esqs[0]
        self.sfts = self.esqs[1]
        temp_inds = np.where((self.qrng<qmax)&(self.qrng>qmin)) 
        self.expq = self.qrng[temp_inds].copy()
        self.exps = self.sfts[temp_inds].copy()
        self.expp = np.square(self.exps)
        self.mcsq = self.expq*0.
        self.escl = self.exps/np.average(self.exps)

    def calc(self, pos, kv, c):
        tmprad    = np.array([np.linalg.norm(ip-c) for ip in pos])
        self.mcsq = np.array([np.abs(np.sum(np.exp(-1.j*iq*tmprad)))**2. for iq in self.expq])
        self.mcsq = self.escl*self.mcsq/np.average(self.mcsq)

    def cmp(self, pos, kv, c):
        self.calc(pos, kv, c)
        dvec = np.abs(np.square(self.mcsq)-self.expp)
        return np.sum(dvec[:-300]/self.exps[:-300])

    def save(self, no=0):
        with open("{}.sq".format(no),'w') as fobj:
            [print("{0:.8f}  {1:.8f} ".format(self.expq[ind],isq),file=fobj) for ind,isq in enumerate(self.mcsq)]


'''
Structure Object 
MS: I like POSCAR format
'''
class vasp:
    #Interface to vasp via ase 
    def __init__(self):
        self.poscar = "./POSCAR"
        self.system = "Simple MC"
        self.totstp = 1000000
        self.tmpcnt = 2.
        self.radius = 1.9
        self.init_setup()
        self.curitr = 0
        self.trstat = self.atmobj.copy()
        self.crstat = self.atmobj.copy()
        self.tpstat = self.atmobj.copy()

    def init_setup(self):
        if os.path.isfile(self.poscar):
            self.atmobj = rdv(self.poscar)
        else:
            tmpatmobj   = Atoms('Ga',
                                cell=[(2.43,0.,0.),(0.,2.43,0.),(0.,0.,2.43)],
                                pbc=True)
            self.atmobj = tmpatmobj.repeat((10,10,10))
        self.totatm = self.atmobj.get_global_number_of_atoms()
        self.atmlst = np.array(range(self.totatm))
        self.atmwgt = np.ones(self.totatm)
        self.center = np.matmul([1.,1.,1.],self.atmobj.cell*0.5)
        self.repvec = np.array(2.*np.pi*self.atmobj.cell.reciprocal())

    def storage_pos(self,x=""):
        if x=="":
            fname = "POSCAR.cnf_{}.vasp".format(self.curitr)
        else:
            fname = "POSCAR.{}.vasp".format(x)
        self.crstat.set_positions(self.crstat.get_positions(wrap=True))
        wrv(fname, self.crstat)


'''
Simple MC Engine
'''
from scipy import constants as c
from scipy.signal import convolve
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

class spmc:
    def __init__(self,vobj="",batch=0.0,mvobj=""):
        self.var = 10.0
        if vobj != "" :
            self.nsmp = int(np.ceil(batch*vobj.totatm)) if batch > 0 else 1
            self.atom = vobj
            self.mvls = mvobj if mvobj != "" else vobj.atmlst 
            self.qobj = sqs()
            self.mind = []
            self.aind = []
            self.traj = []
            self.sqtj = []
            self.chit = 0.01
            self.citr = 0 
        else:
            print("Error Call!!!! Check the arguments!!!")
    
    def _adv(self):
        self.citr = self.citr + 1 
        np.random.seed()
        move_ind  = np.random.choice(self.mvls, self.nsmp)
        self.mind.append(move_ind)
        curr_pos  = self.move(move_ind)
        if self.detect_prob():
            self.atom.trstat = cp.deepcopy(self.atom.crstat) #Backup structure
            self.atom.crstat = cp.deepcopy(self.atom.tpstat) #Update structure 
            self.tchi        = self.cchi 
            self.aind.append(move_ind)
            self.traj.append(curr_pos)
            self.atom.storage_pos(self.citr)
            self.qobj.save(self.citr)



    def move(self, aind):
        self.smv(aind)
        while self.isclash(): 
            self.smv(aind)

    def smv(self,aind):
       #[self.atom.tpstat[ia].set_position(np.random.normal(self.atom.crstat.positions[ia], self.var)) for ia in aind]
        tmv       = np.zeros((self.atom.totatm,3))
        tmv[aind] = np.random.RandomState(42).normal(scale=self.var,size=(self.nsmp,3))
        self.atom.tpstat.set_positions(self.atom.crstat.get_positions(wrap=True)+tmv) 

    def isclash(self):
        dmat = self.atom.tpstat.get_all_distances(mic=True, vector=False)
        mind = np.min(dmat)
        if 0.0 < mind <= self.atom.radius:
            return True # Atom clashed
        else :
            return False 
        
    def detect_prob(self):
        self.calc_sq()
        print("Old chi : {0:.12f} --- New chi : {1:.12f}".format(self.tchi, self.cchi))
        if self.cchi < self.tchi : return True
        return np.random.rand() < min([1, np.exp(self.dchi)])

    def calc_sq(self,typ=2):
        if typ > 1 :
            self.cchi = self.qobj.cmp(self.atom.tpstat.get_positions(wrap=True),self.atom.repvec,self.atom.center)
            self.dchi = (self.tchi - self.cchi)/self.atom.tmpcnt

    def runme(self):
        self.tchi = self.qobj.cmp(self.atom.trstat.get_positions(wrap=True), self.atom.repvec,self.atom.center)
        self.cchi = self.qobj.cmp(self.atom.crstat.get_positions(wrap=True), self.atom.repvec,self.atom.center)
        while self.citr < self.atom.totstp and self.cchi > self.chit :
            print("Start {} Loop!!!".format(self.citr))
            self._adv()
            sys.stdout.flush()

if __name__ == "__main__":
       vobj = vasp()
       mrmc = spmc(vobj, 0.0010)
       mrmc.runme()
