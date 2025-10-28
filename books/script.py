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
    k=0.1
    parameter = [(h, k) for h in np.arange(0, 1.1, 0.1)]
    task_id = "0_1"
    dmrg_main(parameter, "tot", task_id=task_id)

