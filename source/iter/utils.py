import numpy as np
import sys
sys.path.append("..")
import config
from kabsch import kabsched_Q, rmsd, kabsch_rmsd

from init_Q import q_config, lengths_config, rads_config

np_norm = np.linalg.norm

def structure_test(QMatrix):
	Q = np.copy(QMatrix)
	nr_Q_points = QMatrix.shape[0]
	for i in range(nr_Q_points-1):
		l = np.linalg.norm(QMatrix[i+1]-QMatrix[i])
		if np.abs(l-lengths_config[i%3]) > 1e-4:
			return False

	for i in range(nr_Q_points-2):
		p1 = QMatrix[i] - QMatrix[i+1]
		p2 = QMatrix[i+2] - QMatrix[i+1]
		cos_theta = np.dot(p1, p2) / (np_norm(p1) * np_norm(p2))
		theta = np.arccos(cos_theta)
		if np.abs(theta - rads_config[i%3]) > 1e-4:
			return False
	return True


def save_structure(filename, Matrix, rmsd_value):
	f = open(filename, 'w')
	f.write('rmsd: '+str(rmsd_value)+'\n')
	f.write(str(Matrix.shape[0])+'\n')
	for i in range(Matrix.shape[0]):
		for each in Matrix[i]:
			f.write(str(each) + ' ')
		f.write('\n')
	f.close()


def load_structure(filename):
	f = open(filename, 'r').read()
	fl = f.splitlines()

	points_count = int(fl[1])

	NpArray = []

	for i in range(2, points_count+2):
		np_point = np.array(map(lambda x: float(x), fl[i].split()))
		NpArray.append(np_point)

	return np.array(NpArray)


