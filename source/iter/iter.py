"""

iterate to get better result of Q structure 
approach the best structure by multi_time iterations

"""

import numpy as np
import sys
sys.path.append("..")
import config
from kabsch import kabsched_Q, rmsd


def wrap_config():
	cf = config.config
	wrapped_config = []
	for i in range(3):
		paras = {}
		paras["side_lengths"] = [cf[i%3], cf[(i+1)%3]]
		paras["rad"] 		  = cf[i%3 + 3]

		wrapped_config.append(paras)
	return wrapped_config

q_config = wrap_config()


def construct_Q_from(P):
	'''
	construct a init Q for further refinement
	'''
	lengths = q_config[0]["side_lengths"]
	p0 		= P[0]
	p1 		= P[1]
	q1 		= p0 + (p1 - p0) / np.linalg.norm(p1-p0) * lengths[0]
	
	Q = np.array([P[0], q1])
	Q = kabsched_Q(Q, P[0:2])

	for i in range(2, P.shape[0]):
		tmp_p = P[0:i+1]
		Q = expand_Q_from(Q, tmp_p)
		for _ in range(5):
			Q = expand_Q_from(Q[0:-1], tmp_p)

	print rmsd(Q, P)
	from IPython import embed; embed()

	return Q

def expand_Q_from(Q, P):
	'''
	Q 		 : k * 3
	P 		 : k+1 * 3
	return 
	Q 		 : k+1 * 3
	'''
	nr_points = Q.shape[0]
	lengths = q_config[(nr_points-1)%3]["side_lengths"]
	rad 	= q_config[(nr_points-1)%3]["rad"]

	center_para = lengths[1] * np.cos(rad) / lengths[0]

	center_of_next_Q = Q[-1] * (1 - center_para) + Q[-2] * center_para

	proj_para = np.dot(P[-1]-Q[-1], Q[-1]-Q[-2]) / (lengths[0]*lengths[0])

	proj = Q[-1] + proj_para * (Q[-1]-Q[-2])

	new_q_point = center_of_next_Q + (P[-1]-proj) * lengths[1]*np.sin(rad) / np.linalg.norm(P[-1]-proj)


	Q = np.vstack([Q, new_q_point])

	Q = kabsched_Q(Q, P)
	return Q







if __name__ == "__main__":
	import argparse
	import file_reader as fr

	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-p', action='store', dest='p_filename', type=str, help='')
	args = parser.parse_args()

	p_filename = args.p_filename

	PMatrix = fr.get_matrix_in_file(p_filename)

	QMatrix = construct_Q_from(PMatrix)



