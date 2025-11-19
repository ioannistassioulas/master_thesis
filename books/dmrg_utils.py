import numpy as np
import os
import shutil
from dmrg.MPS import MPS 
from dmrg.cont import CONT
from dmrg.dmrg import dmrg
from dmrg.obs import observables
from datetime import datetime
from tqdm import tqdm

class MPO_main():
    Id = np.identity(2)
    X = np.array([[0, 1], [1, 0]])
    Y = np.array([[0, 0-1j], [0+1j, 0]])
    Z = np.array([[1, 0], [0, -1]])
    
    def __init__(self, h, k_left, k_right, J, pol=None, d=2):
        self.h = h
        self.k_left = k_left
        self.k_right = k_right
        self.J = J

        self.pol = pol
        self.d = d

    def Wl(self):
        Wleft = np.zeros((2, 2, 9),dtype='complex')
        Wleft[:, :, 0] = MPO_main.Id
        Wleft[:, :, 1] = MPO_main.X
        Wleft[:, :, 2] = MPO_main.Y
        Wleft[:, :, 3] = MPO_main.Z
        Wleft[:, :, 8] = self.h * MPO_main.Z

        if self.pol == 'tot':
            Wleft[:,:,8] += 25*MPO_main.X

        return Wleft
    
    def Wr(self):
        Wright = np.zeros((2, 2, 9),dtype='complex')
        Wright[:, :, 0] =  self.h * MPO_main.Z
        Wright[:, :, 1] =  self.J * MPO_main.X
        Wright[:, :, 4] =  self.k_right * MPO_main.Y
        Wright[:, :, 5] = -self.k_left * MPO_main.X
        Wright[:, :, 6] =  self.k_left * MPO_main.Z
        Wright[:, :, 7] = -self.k_right * MPO_main.X
        Wright[:, :, 8] =  MPO_main.Id

        if self.pol == 'tot':
            Wright[:,:,0] -= 25*MPO_main.X

        return Wright
    
    def mpo(self, p=None):
        MPO = np.zeros((2, 2, 9, 9),dtype='complex')

        # All interactions
        MPO[:,:,0, 0] =  MPO_main.Id
        MPO[:,:,0, 1] =  MPO_main.X
        MPO[:,:,0, 2] =  MPO_main.Y
        MPO[:,:,0, 3] =  MPO_main.Z
        MPO[:,:,0, 8] =  self.h * MPO_main.Z

        MPO[:,:,1, 4] =  MPO_main.Id
        MPO[:,:,1, 6] =  MPO_main.Y
        MPO[:,:,2, 5] =  MPO_main.Id
        MPO[:,:,3, 7] =  MPO_main.Y

        MPO[:,:,1, 8] =  self.J * MPO_main.X
        MPO[:,:,4, 8] =  self.k_right * MPO_main.Y
        MPO[:,:,5, 8] = -self.k_left * MPO_main.X
        MPO[:,:,6, 8] =  self.k_left * MPO_main.Z
        MPO[:,:,7, 8] = -self.k_right * MPO_main.X 
        MPO[:,:,8, 8] =  MPO_main.Id

        return MPO


def dmrg_main(L, par, pol, task_id = '0_0', J = 1):

    start_time = datetime.now()

    # Base output folder for this task
    base_path = task_id
    # os.makedirs(base_path, exist_ok=True)

    path_out = os.path.join(base_path, 'OUT')
    os.makedirs(path_out, exist_ok=True)
    path_mps = os.path.join(base_path, 'MPS')
    os.makedirs(path_mps, exist_ok=True)
    path_cont = os.path.join(base_path, 'CONT')
    os.makedirs(path_cont, exist_ok=True)

    # Define the system size and bond dimension
    L = L  # check
    chi = 200

    h_x, k_l, k_r = par

    # define the par output
    path_par = os.path.join(path_out, f'out_{h_x:.3f}_{k_l:.3f}_{k_r:.3f}/')
    os.makedirs(path_par, exist_ok=True)

    # initialise the MPS for the indicated chain length
    mps = MPS(L)

    # define the MPO 
    h = MPO_main(h=h_x,k_left=k_l, k_right=k_r, J=J, pol=pol)

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
        _,_,_ = sys.step2sites(site,dir=dir)
        
        # write sweep energy
        # with open(path_par + 'E_sweep.txt','a') as f:
        #     f.write(f'{E} \n')

    # Set up counter and energy check
    En_temp = np.zeros(2*L-8)
    En_temp[0] = 1

    # Now we can sweep the system  (ideally until convergence)
    sweep = 0
    while np.abs(En_temp[0] - En_temp[-1]) > 1e-7:
        j = 0
        for site,dir in mps.sweep():
            En_temp[j],S,_ = sys.step2sites(site,dir=dir)
            # write sweep energy
            with open(path_par + 'E_sweep.txt','a') as f:
                f.write(f'{En_temp[j]} \n')
            j +=1
        # increase system bond dimension for more accuracy
        # sys.chi += 100 - adds too much computation time
        if sweep > 4:
            break
        sweep += 1

    # define observables
    obs = observables(mps)

    # Final sweep to store observables
    for site,dir in mps.right_sweep():
        En,S,_ = sys.step2sites(site,dir=dir,stage='Final')

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

    end_time = datetime.now()
    diff = end_time - start_time
    total_seconds = int(diff.total_seconds())

    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60

    print(f'Parameter ({h_x:.3f}, {k_l:.3f}, {k_r:.3f}) done !!! Time elapsed: {hours:02}:{minutes:02}:{seconds:02}')

def dmrg_line(L, pol, scan_var, values, opp, set, home, J=1):
    ''' Perform a single line of parameter values for DMRG'''

    folder = f"dmrg_results/L{L}_{opp}{set}_{scan_var}{np.min(values)}-{np.max(values)}"
    if os.path.exists(os.path.join(os.getcwd(), folder)):
        return 1
    # determine h and k
    for ind, val in tqdm(enumerate(values)):
        if scan_var == "h":
            h = val
            k_l = k_r = set
        else:
            h = set     # fixed field, change if needed
            k_l = k_r = val

        # create folder for results
        task_id = f"L{L}_{opp}{set}_{scan_var}{np.min(values)}-{np.max(values)}/{h:.2f}_{k_l:.2f}_{k_r:.2f}"
        workdir = f"run_{task_id}"
        os.makedirs(workdir, exist_ok=True)
        os.chdir(workdir)

        parameter = (h,k_l,k_r)
        dmrg_main(L, parameter, pol, task_id='.')

        os.chdir("../..")
    
    # go back home
    os.chdir(home)
