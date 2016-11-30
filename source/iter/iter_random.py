"""

iterate to get better result of Q structure 
approach the best structure by multi_time iterations

"""

import numpy as np
import sys
sys.path.append("..")
import config
from kabsch import kabsched_Q, rmsd, kabsch_rmsd
from init_Q import construct_Q_from, randomly_construct_Q_from
from optimize import optimize, optimize_accelerate
from utils import structure_test, save_structure, load_structure


if __name__ == "__main__":
	import argparse
	import time
	import file_reader as fr

	'''
	-p pmatrix file
	-c continue optimize model or not
	-m if need continue, the input save model
	-r the rounds QMatrix need to optimize
	-s QMatrix save file

	example:
	python iter.py -p 0.txt  -r 100
	python iter.py -p 0.txt  -r 100 -c -m default_model.txt
	'''
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-p', action='store', dest='p_filename', type=str, help='')
	parser.add_argument('-c', action='store_true', dest='need_continue', default=False, help='')
	parser.add_argument('-m', action='store', dest='model_file', default='', type=str, help='')
	parser.add_argument('-r', action='store', dest='rounds', default=100, type=int, help='')
	parser.add_argument('-s', action='store', dest='save_filename', default='default_model.txt', type=str, help='')
	args = parser.parse_args()

	p_filename 	  = args.p_filename
	need_continue = args.need_continue
	model_file 	  = args.model_file
	rounds		  = args.rounds
	save_filename = args.save_filename

	PMatrix = fr.get_matrix_in_file(p_filename)

	Q_min = 0
	rmsd_min = 1000

	# if not need_continue:
	for _ in range(1):
		QMatrix = randomly_construct_Q_from(PMatrix)
		QMatrix = optimize(QMatrix, PMatrix, rounds)
		rmsd_value = rmsd(QMatrix, PMatrix)
		if rmsd_value < rmsd_min:
			rmsd_min = rmsd_value
			Q_min = QMatrix
			print "min: ", rmsd_min
	# else:
		# QMatrix = load_structure(model_file)
		# QMatrix = optimize(QMatrix, PMatrix, rounds)

	if structure_test(QMatrix):
		print 'valid structure'
		save_structure(save_filename, QMatrix, rmsd(QMatrix, PMatrix))
	else:
		print 'invalid structure !!!!!!'


	from IPython import embed; embed()


