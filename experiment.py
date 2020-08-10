import numpy as np
from utils import *


data = np.loadtxt('../BASSET-noEmbds_hyperParams.txt',dtype=str,delimiter='\t')
get_params_dict(data)