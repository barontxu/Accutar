"""

molecule:
restore XYZ coordinates as a numpy matrix
consturct complicated structure by inheritance
provide matrix

"""

import numpy as np

class Molecule(object):
	"""docstring for Molecule"""
	def __init__(self, xyz_matrix):
		self.points_matrix = xyz_matrix

	@property
	def matrix(self):
		return self.points_matrix


class MoleculeQ(object):
	"""docstring for MoleculeQ"""
	def __init__(self, arg):
		super(MoleculeQ, self).__init__()
		self.arg = arg
		