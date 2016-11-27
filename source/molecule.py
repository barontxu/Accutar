"""

molecule:
restore XYZ coordinates as a numpy matrix
consturct complicated structure by inheritance
provide matrix

"""

import numpy as np

class Molecule(object):
	"""docstring for Molecule"""
	points_matrix = np.asarray([0])

	def __init__(self, xyz_matrix):
		self.points_matrix = xyz_matrix

	@property
	def matrix(self):
		return self.points_matrix


class MoleculeQ(Molecule):
	"""docstring for Molecule"""
	def __init__(self, lengths, rads):
		self.lengths = lengths
		self.rads = rads

	def set_matrix(self, xyz_matrix):
		self.points_matrix = xyz_matrix

	
		
