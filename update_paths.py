import os
import glob
import numpy as np
import subprocess as sp

if __name__ == "__main__":
    paths = np.genfromtxt('path_config.dat', dtype = 'str', delimiter = ',')
    for old_path, new_path in paths:
        print old_path
        print new_path
        sed_str = 'sed -n \'s:%s:%s:p\' text_path_config.txt >text_path_config.tx' %(old_path, new_path)
        print sed_str
        sp.call(sed_str, shell = True)
