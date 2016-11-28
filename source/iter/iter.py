"""

iterate to get better result of Q structure 
approach the best structure by multi_time iterations

"""

import numpy as np
import sys
sys.path.append("..")
import config
from kabsch import kabsched_Q, rmsd, kabsch_rmsd
from init_Q import construct_Q_from
from optimize import optimize


if __name__ == "__main__":
	import argparse
	import file_reader as fr

	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-p', action='store', dest='p_filename', type=str, help='')
	args = parser.parse_args()

	p_filename = args.p_filename

	PMatrix = fr.get_matrix_in_file(p_filename)

	QMatrix = construct_Q_from(PMatrix)

	QMatrix = optimize(QMatrix, PMatrix)
	from IPython import embed; embed()


