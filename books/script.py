# Python script to submit to HPC

# from dmrg import MPS, MPO_TFI, CONT, dmrg, observables
from dmrg.MPS import MPS 
# from dmrg.MPO import MPO_TFI
from dmrg.cont import CONT
from dmrg.dmrg import dmrg
from dmrg.obs import observables

# from utils import *
import os
import shutil
import numpy as np
import argparse
from itertools import product

# Define folders where to store simulation data
path_mps = 'MPS'
path_cont = 'CONT'
path_out = 'OUT/'

class MPO_TFI():
    
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
        Wleft = np.zeros((2, 2, 9))
        Wleft[:, :, 0] = MPO_TFI.Id
        Wleft[:, :, 1] = MPO_TFI.X
        Wleft[:, :, 2] = MPO_TFI.Y
        Wleft[:, :, 3] = MPO_TFI.Z
        Wleft[:, :, 8] = self.h * MPO_TFI.Z

        if self.pol == 'tot':
            Wleft[:,:,8] -= 25*MPO_TFI.X

        return Wleft
    
    def Wr(self):
        Wright = np.zeros((2, 2, 9))
        Wright[:, :, 0] = self.h * MPO_TFI.Z
        Wright[:, :, 5] = MPO_TFI.Z
        Wright[:, :, 6] = MPO_TFI.Y
        Wright[:, :, 7] = MPO_TFI.X
        Wright[:, :, 8] = MPO_TFI.Id

        if self.pol == 'tot':
            Wright[:,:,0] += 25*MPO_TFI.X

        return Wright
    
    def mpo(self):
        MPO = np.zeros((2, 2, 9, 9))

        # All interactions
        MPO[:,:,0, 0] = MPO_TFI.Id
        MPO[:,:,0, 1] =  self.k * MPO_TFI.X
        MPO[:,:,0, 2] = -self.k * MPO_TFI.Z
        MPO[:,:,0, 3] =  self.k * MPO_TFI.X
        MPO[:,:,0, 4] = -self.k * MPO_TFI.Y
        MPO[:,:,0, 7] =  self.J * MPO_TFI.X
        MPO[:,:,0, 8] =  self.h * MPO_TFI.Z

        MPO[:,:,1, 5] = MPO_TFI.Y
        MPO[:,:,2, 7] = MPO_TFI.Y
        MPO[:,:,3, 6] = MPO_TFI.Id
        MPO[:,:,4, 7] = MPO_TFI.Id
        MPO[:,:,5, 8] = MPO_TFI.Z
        MPO[:,:,6, 8] = MPO_TFI.Y
        MPO[:,:,7, 8] = MPO_TFI.X
        MPO[:,:,8, 8] = MPO_TFI.Id

        return MPO


def TFIM_DMRG(h, k, count, pol, J = 1):
    # create folder out if not present
    if os.path.isdir(path_out):
        shutil.rmtree(path_out)
    os.mkdir(path_out)

    # Define the system size and bond dimension
    L = 50
    chi = 200

    # Define parameters of h and k (look at all combinations)
    par = list(product(np.linspace(0, h, count), np.linspace(-k, k, count+10)))


    for h_x, k_x in par:
        # define the par output
        path_par = path_out + f'out_{h_x:.2f}_{k_x:.2f}/' 
        os.mkdir(path_par)

        # initialise the MPS for the indicated chain length
        mps = MPS(L)

        # define the MPO 
        h = MPO_TFI(h=h_x,k=k_x, J=J, pol=pol)
        print(k_x)
        # define the contractions (it needs an mps and a MPO class as imputs)
        cont = CONT(mps=mps,H=h)

        # Initialize your dmrg (set low bond dimension to make the system grow faster)
        sys = dmrg(cont=cont,chi=10,cut=1e-12)

        # Grow the system up to the desired dimension
        En = sys.infinite()

        # open the energy sweep file and write the starting energy
        with open(path_par + 'E_sweep.txt','w') as f:
            f.write(f'{En[-1]} \n')

        # Increase the bond dimension to the desired one 
        sys.chi = chi

        # run the first one and half sweep
        for site,dir in mps.first_sweep():
            E,_ = sys.step2sites(site,dir=dir)
            
            # write sweep energy
            with open(path_par + 'E_sweep.txt','a') as f:
                f.write(f'{E} \n')

        # Set up counter and energy check
        counter = 0
        En_temp = np.zeros(2*L-8)
        En_temp[0] = En[-1]

        # Now we can sweep the system  (ideally until convergence)
        while np.abs(En_temp[0] - En_temp[-1]) > 1e-10:
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
            
            counter += 1

        # define observables
        obs = observables(mps)

        # Final sweep to store observables
        for site,dir in mps.right_sweep():
            _,S = sys.step2sites(site,dir=dir,stage='Final')

            # Store local magnetization
            with open(path_par + 'Z.txt','a') as fz:
                    fz.write(f'{site} {obs.single_site(site,h.Z).real} \n')

            # Store entanglement entropy
            with open(path_par + 'S.txt','a') as fz:
                    fz.write(f'{site} {site+1} {S} \n')

            # Store all two point correlations from site

            obs.all_corr(path_par + 'ZZ.txt',site,obs1=h.Z)


        print(f'Parameter ({h_x:.2f}, {k_x:.2f}) done !!!')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("h", type=float)
    parser.add_argument("k", type=float)
    parser.add_argument("count", type=int)
    parser.add_argument("pol", type=str)
    
    args = parser.parse_args()

    TFIM_DMRG(args.h, args.k, args.count, args.pol)