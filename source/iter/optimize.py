import numpy as np
import sys
sys.path.append("..")
import config
from kabsch import kabsched_Q, rmsd, kabsch_rmsd

from init_Q import q_config


def optimize_structure_by_pos(Q_origin, P_origin, k_pos):
	assert k_pos > 1
	assert k_pos < P_origin.shape[0] - 1
	
	nr_P_row = P_origin.shape[0]

	P = np.copy(P_origin)
	Q = np.copy(Q_origin)
	k = k_pos

	v1 = (Q[k]-Q[k-1])
	u1 = v1 / np.linalg.norm(v1)

	u2 = np.array([1, 0, 0])

	if np.abs(u1[1]) < np.abs(u1[0]):
		u2 = np.array([0, 1, 0])
	
	u2 = np.array([u1[1] * u2[2] - u1[2] * u2[1], u1[2] * u2[0] - u1[0] * u2[2], u1[0] * u2[1] - u1[1] * u2[0]])
	u2 = u2 / np.linalg.norm(u2)

	u3 = np.array([u1[1] * u2[2] - u1[2] * u2[1], u1[2] * u2[0] - u1[0] * u2[2], u1[0] * u2[1] - u1[1] * u2[0]])

	U = np.array([u1, u2, u3])

	Qk = Q[k+1: nr_P_row]
	Pk = P[k+1: nr_P_row]

	q_point_k_wrap = np.copy(np.array([Q[k]]))
	Qk = Qk - np.repeat(q_point_k_wrap, [Qk.shape[0]], axis = 0)
	Pk = Pk - np.repeat(q_point_k_wrap, [Qk.shape[0]], axis = 0)

	Qk_U = np.dot(Qk, U.transpose())
	Pk_U = np.dot(Pk, U.transpose())

	alpha = np.dot(Qk_U[:, 1], Pk_U[:, 2]) - np.dot(Qk_U[:, 2], Pk_U[:, 1])
	beta  = np.dot(Qk_U[:, 1], Pk_U[:, 1]) + np.dot(Qk_U[:, 2], Pk_U[:, 2])

	cos_theta = beta / np.sqrt(alpha**2 + beta**2)
	sin_theta = alpha / np.sqrt(alpha**2 + beta**2)

	pre_A = np.array([[1, 0, 		 0],
			 		  [0, cos_theta,  sin_theta],
			 		  [0, -sin_theta, cos_theta]])

	A = np.dot(np.dot(U.transpose(), pre_A), U)

	Qk = np.dot(Qk, A)
	Qk = Qk + np.repeat(q_point_k_wrap, [Qk.shape[0]], axis = 0)
	Q[k+1: nr_P_row] = Qk

	Q = kabsched_Q(Q, P)
	# from IPython import embed; embed()
	return Q


def optimize(Q_origin, P_origin):
	rounds = 1
	P = np.copy(P_origin)
	Q = np.copy(Q_origin)
	print rmsd(Q, P)
	print "--------------------------------------"
	for _ in range(rounds):
		for k in range(2, P.shape[0]-1):
			Q = optimize_structure_by_pos(Q, P, k)
			print rmsd(Q, P)

	print rmsd(Q, P)
	return Q







