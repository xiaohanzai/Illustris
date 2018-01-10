import numpy as np
from Illustris.utils.util_general import xyz2cyl, xyz2sph
from Illustris.utils.util_illustris import getData, getVc, findCenter, calcShape

Rb_all = 20.

class StarData(object):
	def __init__(self, x, v, mpart, Rb = Rb_all, **kwargs):
		'''
		Positions, velocities, and masses of the particles must be provided.
		Other quantities like metallicity, luminosity, etc. can be provided using
		  in **kwargs.
		Rb is used for calculating shape.
		'''
		# initialize x, v, particle mass, and other things
		self.x = x
		self.v = v
		self.mpart = mpart
		for name in kwargs.keys():
			setattr(self, name, kwargs[name])

		self.xc = kwargs.get('xc', findCenter(self.x, mpart=self.mpart))
		self.x_ori = self.x - self.xc
		self.vc = kwargs.get('vc', getVc(self.x_ori, self.v, mpart=self.mpart))
		self.v_ori = self.v - self.vc

		# initialize shape
		self.ba, self.ca, angle_l, angle_s, self.Tiv = calcShape(self.x_ori, Rb = Rb)
		self.shape = 'prolate' # currently self.Tiv makes z axis the longest axis

	def addStarAttrs(self, **kwargs):
		for name in kwargs.keys():
			setattr(self, name, kwargs[name])

	def convert2cyl(self, shape):
		'''
		Given galaxy shape, calculate position vectors in principal axes coordinates
		  and cylindrical coordinates.
		'''
		if shape not in ['oblate', 'prolate']:
			raise ValueError("Shape must be oblate or prolate.")
		if self.shape != shape:
			self.shape = shape
			self.Tiv = self.Tiv[[2,1,0],:]
		self.x_prin = np.dot(self.Tiv,(self.x_ori).T).T
		self.v_prin = np.dot(self.Tiv,(self.v_ori).T).T
		self.xcyl, self.vcyl = xyz2cyl(self.x_prin, self.v_prin)

	def convert2sph(self):
		'''
		Calculate position and velocities in spherical coordinates.
		'''
		if not hasattr(self, 'x_prin'):
			print('Galaxy shape not determined. z axis may not be the symmetry axis.')
			x = self.x_ori
			v = self.v_ori
		else:
			x = self.x_prin
			v = self.v_prin
		self.xsph, self.vsph = xyz2sph(x, v)

