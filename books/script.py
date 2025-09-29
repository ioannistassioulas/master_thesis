# Python script to submit to HPC

from dmrg.MPS import MPS 
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

    # for lines
    # parser.add_argument("param", type=float)
    # parser.add_argument("step", type=float)
    # parser.add_argument("axis", type=str)
    
    ## For sanity checks
    parser.add_argument("J", type=float)
    parser.add_argument("h", type=float)
    parser.add_argument("k", type=float)

    parser.add_argument("pol", type=str)

    args = parser.parse_args()

    # for lines
    # parameter = dmrg_lines(args.param, args.step, args.axis) # define whicg parameters the simulation will run
    # dmrg_main(parameter, args.pol)

    ## for sanity checks
    parameter = [(args.h, args.k)]
    dmrg_main(parameter, args.pol, J=args.J)

