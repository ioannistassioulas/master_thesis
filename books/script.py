# Python script to submit to HPC

# from dmrg import MPS, MPO_TFI, CONT, dmrg, observables
from dmrg.MPS import MPS 
from dmrg.MPO import MPO_TFI
from dmrg.cont import CONT
from dmrg.dmrg import dmrg
from dmrg.obs import observables

# from utils import *
import os
import shutil
import numpy as np
import argparse

# Define folders where to store simulation data
path_mps = 'MPS'
path_cont = 'CONT'
path_out = 'OUT/'


def TFIM_DMRG(h_min, h_max):
    # create folder out if not present
    if os.path.isdir(path_out):
        shutil.rmtree(path_out)
    os.mkdir(path_out)
    # Define the parameters, system size and bond dimension
    par = np.linspace(h_min,h_max, int(10 * (h_max - h_min) + 1))
    L = 50
    chi = 200

    for h_x in par:
        # define the par output
        path_par = path_out + f'out_{h_x:.2f}/' 
        os.mkdir(path_par)

        # initialise the MPS for the indicated chain length
        mps = MPS(L)

        # define the MPO 
        h = MPO_TFI(J=1,h_x=h_x,pol='tot')

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
        k = 0 
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
            if k > 5:
                break
            
            k += 1

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


        print(f'Parameter {h_x:.2f} done !!!')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("h_min", type=float)
    parser.add_argument("h_max", type=float)
    
    args = parser.parse_args()

    TFIM_DMRG(args.h_min, args.h_max)