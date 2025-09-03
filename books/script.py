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




if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("h", type=float)
    parser.add_argument("k", type=float)
    parser.add_argument("steps", type=float)
    parser.add_argument("pol", type=str)
    parser.add_argument("axis", type=str)
    
    args = parser.parse_args()

    #
    dmrg_main_lines(args.h, args.k, args.steps, args.pol, args.axis)