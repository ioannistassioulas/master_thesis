
import numpy as np
import os
import shutil
from dmrg.MPS import MPS 
from dmrg.cont import CONT
from dmrg.dmrg import dmrg
from utils import *

class MPO_main():
    
    Id = np.identity(2)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, 0-1j], [0+1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    
    def __init__(self, h, k, J, pol=None, d=2):
        self.h = h
        self.k = k
        self.J = J

        self.pol = pol
        self.d = d

    def Wl(self):
        Wleft = np.zeros((2, 2, 10),dtype='complex')
        Wleft[:, :, 0] = MPO_main.Id
        Wleft[:, :, 1] = MPO_main.X
        Wleft[:, :, 2] = MPO_main.X
        Wleft[:, :, 3] = MPO_main.Y
        Wleft[:, :, 4] = MPO_main.Z
        Wleft[:, :, 9] = self.h * MPO_main.Z

        if self.pol == 'tot':
            Wleft[:,:,9] += 25*MPO_main.X

        return Wleft
    
    def Wr(self):
        Wright = np.zeros((2, 2, 10),dtype='complex')
        Wright[:, :, 0] =  self.h * MPO_main.Z
        Wright[:, :, 1] =  self.J * MPO_main.X
        Wright[:, :, 5] =  self.k * MPO_main.Y
        Wright[:, :, 6] =  self.k * MPO_main.Z
        Wright[:, :, 7] = -self.k * MPO_main.X
        Wright[:, :, 8] = -self.k * MPO_main.X
        Wright[:, :, 9] =  MPO_main.Id

        if self.pol == 'tot':
            Wright[:,:,0] -= 25*MPO_main.X

        return Wright
    
    def mpo(self, p=None):
        MPO = np.zeros((2, 2, 10, 10),dtype='complex')

        # All interactions
        MPO[:,:,0, 0] =  MPO_main.Id
        MPO[:,:,0, 1] =  MPO_main.X
        MPO[:,:,0, 2] =  MPO_main.X
        MPO[:,:,0, 3] =  MPO_main.Y
        MPO[:,:,0, 4] =  MPO_main.Z
        MPO[:,:,0, 9] =  self.h * MPO_main.Z
        

        MPO[:,:,1, 5] =  MPO_main.Id
        MPO[:,:,2, 6] =  MPO_main.Y
        MPO[:,:,3, 7] =  MPO_main.Id
        MPO[:,:,4, 8] =  MPO_main.Y

        MPO[:,:,1, 9] =  self.J * MPO_main.X
        MPO[:,:,5, 9] =  self.k * MPO_main.Y
        MPO[:,:,6, 9] =  self.k * MPO_main.Z
        MPO[:,:,7, 9] = -self.k * MPO_main.X
        MPO[:,:,8, 9] = -self.k * MPO_main.X 
        MPO[:,:,9, 9] =  MPO_main.Id

        return MPO

# all of my own defined classes (including observables and mpos) I leave here
def dmrg_main(par, pol, task_id = '0_0', J = 1):

    # Base output folder for this task
    base_path = task_id
    os.makedirs(base_path, exist_ok=True)

    path_out = os.path.join(base_path, 'OUT')
    os.makedirs(path_out, exist_ok=True)
    path_mps = os.path.join(base_path, 'MPS')
    os.makedirs(path_mps, exist_ok=True)
    path_cont = os.path.join(base_path, 'CONT')
    os.makedirs(path_cont, exist_ok=True)

    

    # Define the system size and bond dimension
    L = 20  # check
    chi = 100

    for h_x, k_x in par:
        # define the par output
        path_par = os.path.join(path_out, f'out_{h_x:.3f}_{k_x:.3f}/')
        os.makedirs(path_par, exist_ok=True)
        os.chdir(base_path)
        # initialise the MPS for the indicated chain length
        mps = MPS(L)

        # define the MPO 
        h = MPO_main(h=h_x,k=k_x, J=J, pol=pol)

        # define the contractions (it needs an mps and a MPO class as imputs)
        cont = CONT(mps=mps,H=h)

        # randomize initial MPS and CONT
        mps.random()
        cont.random()

        # Initialize your dmrg (set low bond dimension to make the system grow faster)
        sys = dmrg(cont=cont,chi=10,cut=1e-12)

        # Grow the system up to the desired dimension
        # En = sys.infinite()

        En = [0.0]
        # open the energy sweep file and write the starting energy
        with open(path_par + 'E_sweep.txt','w') as f:
            f.write(f'{En[-1]} \n')

        # Increase the bond dimension to the desired one 
        sys.chi = chi

        # run the first one and half sweep - do not record!
        for site,dir in mps.first_sweep():
            _,_ = sys.step2sites(site,dir=dir)
            
            # write sweep energy
            # with open(path_par + 'E_sweep.txt','a') as f:
            #     f.write(f'{E} \n')

        # Set up counter and energy check
        counter = 0
        En_temp = np.zeros(2*L-8)
        En_temp[0] = 1

        # Now we can sweep the system  (ideally until convergence)
        while np.abs(En_temp[0] - En_temp[-1]) > 1e-9:
            j = 0
            for site,dir in mps.sweep():
                En_temp[j],S = sys.step2sites(site,dir=dir)

                # write sweep energy
                with open(path_par + 'E_sweep.txt','a') as f:
                    f.write(f'{En_temp[j]} \n')

                j +=1

            # set maximum number of sweeps
            if counter > 5:
                break
            
            # increase system bond dimension for more accuracy
            sys.chi += 50
            counter += 1

        # define observables
        obs = observables(mps)

        # Final sweep to store observables
        for site,dir in mps.right_sweep():
            _,S = sys.step2sites(site,dir=dir,stage='Final')

            # Store local magnetization
            with open(path_par + 'X.txt','a') as fz:
                    fz.write(f'{site} {obs.single_site(site,h.X).real} \n')
            with open(path_par + 'Y.txt', 'a') as fz:
                    fz.write(f'{site} {obs.single_site(site,h.Y).real} \n')
            with open(path_par + 'Z.txt', 'a') as fz:
                    fz.write(f'{site} {obs.single_site(site,h.Z).real} \n')

            # Store entanglement entropy
            with open(path_par + 'S.txt','a') as fz:
                    fz.write(f'{site} {site+1} {S} \n')

            # Store all two point correlations from site
            obs.all_corr(path_par + 'XX.txt', site, string=h.Id , obs1=h.X)
            obs.all_corr(path_par + 'YY.txt', site, string=h.Id , obs1=h.Y)
            obs.all_corr(path_par + 'ZZ.txt', site, string=h.Id , obs1=h.Z)

            obs.all_corr(path_par + 'XZ.txt', site, string=h.Id , obs1=h.X, obs2=h.Z)
            obs.all_corr(path_par + 'XY.txt', site, string=h.Id , obs1=h.X, obs2=h.Y)
            obs.all_corr(path_par + 'YZ.txt', site, string=h.Id , obs1=h.Y, obs2=h.Z)

            obs.all_corr(path_par + 'XYZ.txt',site, string=h.Y, obs1=h.X, obs2=h.Z)

        print(f'Parameter ({h_x:.3f}, {k_x:.3f}) done !!!')

# class MPO_main():
    
    # Id = np.identity(2)
    # X = np.array([[0, 1], [1, 0]])
    # Y = np.array([[0, 0-1j], [0+1j, 0]])
    # Z = np.array([[1, 0], [0, -1]])
    
    # def __init__(self, h, k, J, pol=None, d=2):
    #     self.h = h
    #     self.k = k
    #     self.J = J

    #     self.pol = pol
    #     self.d = d

    # def Wl(self):
    #     Wleft = np.zeros((2, 2, 9),dtype='complex')
    #     Wleft[:, :, 0] = MPO_main.Id
    #     Wleft[:, :, 1] = MPO_main.X
    #     Wleft[:, :, 2] = MPO_main.Y
    #     Wleft[:, :, 3] = MPO_main.Z
    #     Wleft[:, :, 8] = self.h * MPO_main.Z

    #     if self.pol == 'tot':
    #         Wleft[:,:,8] += 25*MPO_main.X

    #     return Wleft
    
    # def Wr(self):
    #     Wright = np.zeros((2, 2, 9),dtype='complex')
    #     Wright[:, :, 0] =  self.h * MPO_main.Z
    #     Wright[:, :, 1] =  self.J * MPO_main.X
    #     Wright[:, :, 4] =  self.k * MPO_main.Y
    #     Wright[:, :, 5] = -self.k * MPO_main.X
    #     Wright[:, :, 6] =  self.k * MPO_main.Z
    #     Wright[:, :, 7] = -self.k * MPO_main.X
    #     Wright[:, :, 8] =  MPO_main.Id

    #     if self.pol == 'tot':
    #         Wright[:,:,0] -= 25*MPO_main.X

    #     return Wright
    
    # def mpo(self, p=None):
    #     MPO = np.zeros((2, 2, 9, 9),dtype='complex')

    #     # All interactions
    #     MPO[:,:,0, 0] =  MPO_main.Id
    #     MPO[:,:,0, 1] =  MPO_main.X
    #     MPO[:,:,0, 2] =  MPO_main.Y
    #     MPO[:,:,0, 3] =  MPO_main.Z
    #     MPO[:,:,0, 8] =  self.h * MPO_main.Z

    #     MPO[:,:,1, 4] =  MPO_main.Id
    #     MPO[:,:,1, 6] =  MPO_main.Y
    #     MPO[:,:,2, 5] =  MPO_main.Id
    #     MPO[:,:,3, 7] =  MPO_main.Y

    #     MPO[:,:,1, 8] =  self.J * MPO_main.X
    #     MPO[:,:,4, 8] =  self.k * MPO_main.Y
    #     MPO[:,:,5, 8] = -self.k * MPO_main.X
    #     MPO[:,:,6, 8] =  self.k * MPO_main.Z
    #     MPO[:,:,7, 8] = -self.k * MPO_main.X 
    #     MPO[:,:,8, 8] =  MPO_main.Id

    #     return MPO

class observables():
    def __init__(self,MPS):
        self.mps = MPS
        self.L = MPS.L 

    def single_site(self,site,obs):
        ten = np.tensordot(self.mps.read(site),self.mps.readS(site),(2,0)) 
        return np.tensordot(np.tensordot(obs,ten,(0,0)),np.conj(ten),((0,1,2),(0,1,2)))

    def all_corr(self,path,site,string,obs1,obs2=None):
        if obs2 is None:
            obs2 = obs1
        ten = np.tensordot(self.mps.read(site),self.mps.readS(site),(2,0)) 
        cont1 = np.tensordot(np.tensordot(obs1,ten,(0,0)),np.conj(ten),((0,1),(0,1)))
        
        
        for i in range(site+1,self.L-1):
            cont2 = np.tensordot(np.tensordot(obs2,self.mps.read(i),(0,0)),np.conj(self.mps.read(i)),((0,2),(0,2)))
            if i > site + 1:
                cont1 = np.tensordot(cont1,np.tensordot(string,self.mps.read(i-1),(0,0)),(0,1))
                cont1 = np.tensordot(cont1,np.conj(self.mps.read(i-1)),((0,1),(1,0))) 
            
            res = np.tensordot(cont1,cont2,((0,1),(0,1)))
        
            with open(path,'a') as f:
                f.write(f'{site} {i} {res}\n')
