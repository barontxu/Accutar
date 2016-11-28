import numpy as np
import sys
sys.path.append("..")
import config
from kabsch import kabsched_Q, rmsd, kabsch_rmsd

from init_Q import q_config


def optimize_structure_by_pos(P_origin, Q_origin, k_pos):
	assert k_pos > 1
	assert k_pos < P_origin.shape[0] - 1
	
	P = np.copy(P_origin)
    Q = np.copy(Q_origin)
    k = k_pos

    v1 = (Q[k]-Q[k-1])
    u1 = v1 / np.linalg.norm(v1)

    u2 = (1, 0, 0)

    if np.abs(u1[1]) < np.abs(u1[0]):
    	
