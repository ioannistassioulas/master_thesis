# Python script to submit to HPC

from dmrg.MPS import MPS 
# from dmrg.MPO import MPO_TFI
# from dmrg.obs import observables
from dmrg.cont import CONT
from dmrg.dmrg import dmrg
from utils import *

import os
import shutil
import numpy as np
import argparse
# from itertools import product

# Define folders where to store simulation data
path_mps = 'MPS'
path_cont = 'CONT'
path_out = 'OUT/'

def dmrg_main_lines(h, k, step, pol, axis, J = 1):
    # create folder out if not present
    if os.path.isdir(path_out):
        shutil.rmtree(path_out)
    os.mkdir(path_out)

    # Define the system size and bond dimension
    L = 21
    chi = 200

    # Define parameters of h and k (look at all combinations)
    par = [(i, k) for i in np.arange(0, h+step, step)] if axis is "h" else [(h, i) for i in np.arange(0, k+step, step)]

    for h_x, k_x in par:
        # define the par output
        path_par = path_out + f'out_{h_x:.2f}_{k_x:.2f}/' 
        os.mkdir(path_par)

        # initialise the MPS for the indicated chain length
        mps = MPS(L)

        # define the MPO 
        h = MPO_main(h=h_x,k=k_x, J=J, pol=pol)

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
            if counter > 8:
                break
            
            counter += 1

        # define observables
        obs = observables(mps)

        # Final sweep to store observables
        for site,dir in mps.right_sweep():
            _,S = sys.step2sites(site,dir=dir,stage='Final')

            # Store local magnetization
            with open(path_par + 'X.txt','a') as fz:
                    fz.write(f'{site} {obs.single_site(site,h.X).real} \n')

            # Store entanglement entropy
            with open(path_par + 'S.txt','a') as fz:
                    fz.write(f'{site} {site+1} {S} \n')

            # Store all two point correlations from site
            obs.all_corr(path_par + 'XX.txt', site, obs1=h.X)
            obs.all_corr(path_par + 'XY.txt', site, obs1=h.X, obs2=h.Y)
            obs.all_corr(path_par + 'YZ.txt', site, obs1=h.Y, obs2=h.Z)
            obs.all_corr(path_par + 'XZ.txt', site, obs1=h.X, obs2=h.Z)

        print(f'Parameter ({h_x:.2f}, {k_x:.2f}) done !!!')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("h", type=float)
    parser.add_argument("k", type=float)
    parser.add_argument("steps", type=float)
    parser.add_argument("pol", type=str)
    parser.add_argument("axis", type=str)
    
    args = parser.parse_args()

    dmrg_main_lines(args.h, args.k, args.steps, args.pol, args.axis)