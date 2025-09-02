
import numpy as np
import os
import shutil

# all of my own defined classes (including observables and mpos) I leave here

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
        Wright[:, :, 0] = self.h * MPO_main.Z
        Wright[:, :, 5] = MPO_main.Z
        Wright[:, :, 6] = MPO_main.Y
        Wright[:, :, 7] = MPO_main.X
        Wright[:, :, 8] = MPO_main.Id

        if self.pol == 'tot':
            Wright[:,:,0] += 25*MPO_main.X

        return Wright
    
    def mpo(self, p=None):
        MPO = np.zeros((2, 2, 9, 9),dtype='complex')

        # All interactions
        MPO[:,:,0, 0] = MPO_main.Id
        MPO[:,:,0, 1] =  self.k * MPO_main.X
        MPO[:,:,0, 2] = -self.k * MPO_main.Z
        MPO[:,:,0, 3] =  self.k * MPO_main.X
        MPO[:,:,0, 4] = -self.k * MPO_main.Y
        MPO[:,:,0, 7] =  self.J * MPO_main.X
        MPO[:,:,0, 8] =  self.h * MPO_main.Z

        MPO[:,:,1, 5] = MPO_main.Y
        MPO[:,:,2, 7] = MPO_main.Y
        MPO[:,:,3, 6] = MPO_main.Id
        MPO[:,:,4, 7] = MPO_main.Id
        MPO[:,:,5, 8] = MPO_main.Z
        MPO[:,:,6, 8] = MPO_main.Y
        MPO[:,:,7, 8] = MPO_main.X
        MPO[:,:,8, 8] = MPO_main.Id

        return MPO

class observables():
    def __init__(self,MPS):
        self.mps = MPS
        self.L = MPS.L 

    def single_site(self,site,obs):
        ten = np.tensordot(self.mps.read(site),self.mps.readS(site),(2,0)) 
        np.tensordot(obs,ten,(0,0))
        return np.tensordot(np.tensordot(obs,ten,(0,0)),np.conj(ten),((0,1,2),(0,1,2)))

    def all_corr(self,path,site,obs1,obs2=None):
        if obs2 == None:
            obs2 = obs1
        ten = np.tensordot(self.mps.read(site),self.mps.readS(site),(2,0)) 
        
        cont1 = np.tensordot(np.tensordot(obs1,ten,(0,0)),np.conj(ten),((0,1),(0,1)))
        
        for i in range(site+1,self.L-1):
            cont2 = np.tensordot(np.tensordot(obs2,self.mps.read(i),(0,0)),np.conj(self.mps.read(i)),((0,2),(0,2)))
            if i > site + 1:
                cont1 = np.tensordot(cont1,self.mps.read(i-1),(0,1))
                cont1 = np.tensordot(cont1,np.conj(self.mps.read(i-1)),((0,1),(1,0))) 
            
            res = np.tensordot(cont1,cont2,((0,1),(0,1)))
        
            with open(path,'a') as f:
                f.write(f'{site} {i} {res}\n')