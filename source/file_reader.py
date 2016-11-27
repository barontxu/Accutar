"""

read file in and parse the floats into numpy array

"""

import numpy as np

def get_matrix_in_file(filename):
	f = open(filename, 'r').read()
	fl = f.splitlines()

	points_count = 3 * int(fl[0])

	NpArray = []

	for i in range(1, points_count+1):
		np_point = np.array(map(lambda x: float(x), fl[i].split()))
		NpArray.append(np_point)

	return np.array(NpArray)
	