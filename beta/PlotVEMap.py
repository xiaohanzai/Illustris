'''
Plot the velocity ellipsoids.
'''
import numpy as np
import numpy.linalg as LA
import os
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
rc('xtick', labelsize=12)
rc('ytick', labelsize=12)
rc('axes', labelsize=15)
from Illustris.beta.StarData_Beta import StarData_Beta, Rb_all
from Illustris.utils.util_illustris import getData
from Illustris.utils.paths import illustris_samplepath, illustris_path
outpath = illustris_samplepath

class PlotVEMap(StarData_Beta):
	def __init__(self, x, v, mpart, shape, Rb = Rb_all, **kwargs):
		super(PlotVEMap, self).__init__(x, v, mpart, shape, Rb, **kwargs)
	
	def plotVE(self, rmin = 2., rmax = Rb_all, N_radius = 10, 
		       thetamin = 10., thetamax = 80., N_theta = 8, 
		       fig = None, ax = None, 
		       massweight = True, Nmin = 100, s = 1.):
		'''
		rmin rmax are in units of kpc.
		s is the scaling, used to scale the size of the velocity ellipsoids.
		Only plot the velocity ellipsoid for bins containing > Nmin particles.
		'''
		# make bins
		theta = np.linspace(thetamin, thetamax, N_theta)*np.pi/180. # in units of rad
		step = (thetamax - thetamin)/N_theta
		radius = np.logspace(np.log10(rmin), np.log10(rmax), N_radius)

		if fig == None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		# make a grid for plotting
		xR = yz = np.arange(0, rmax+1, 0.05)
		xR, yz = np.meshgrid(xR, yz)
		
		# calculate and plot velocity ellipsoids
		for r in radius:
			ax.contour(xR, yz, xR**2 + yz**2 - r**2, [0], colors = 'k')
			dR = dz = min(step*r, 0.8)

			for ang in theta:
				Rc = r*np.cos(ang)
				zc = r*np.sin(ang)

				a11, a12, a13, a22, a23, a33, v1, v2, v3, M, V, Npart = \
				    self.measureV2Tensor(Rc, zc, np.pi, massweight = massweight, 
	                    dR = dR, dz = dz, dphi = 2*np.pi, abs_z = True)
				if Npart < Nmin:
					continue

				a11 -= v1**2
				a12 -= v1*v2
				a13 -= v1*v3
				a22 -= v2**2
				a23 -= v2*v3
				a33 -= v3**2

				# draw the graph
				mat = LA.inv(np.array([[a11,a12,a13],[a12,a22,a23],[a13,a23,a33]]))
				a2 = 0.5*(mat[0,0]+mat[2,2]+np.sqrt((mat[0,0]-mat[2,2])**2+4*mat[0,2]**2))
				b2 = 0.5*(mat[0,0]+mat[2,2]-np.sqrt((mat[0,0]-mat[2,2])**2+4*mat[0,2]**2))
				rotateTheta = np.arccos((mat[0,0]-mat[2,2])/(b2-a2))/2
				if mat[0,2]>0:
					rotateTheta = np.pi-rotateTheta

				b2_scaled = (r*np.pi/2/N_theta*0.4)**2*s
				a2_scaled = a2/b2*b2_scaled

				ax.contour(xR, yz, \
					(xR-Rc)**2*mat[0,0]+(yz-zc)**2*mat[2,2]+2*(xR-Rc)*(yz-zc)*mat[0,2],\
					[a2*b2*b2_scaled/b2], colors = 'k')
				ax.plot([Rc-np.sqrt(a2_scaled)*np.cos(rotateTheta),Rc+np.sqrt(a2_scaled)*np.cos(rotateTheta)], 
					    [zc-np.sqrt(a2_scaled)*np.sin(rotateTheta),zc+np.sqrt(a2_scaled)*np.sin(rotateTheta)], 
					    'r') # plot the major axis
				ax.plot([Rc-np.sqrt(b2_scaled)*np.sin(rotateTheta),Rc+np.sqrt(b2_scaled)*np.sin(rotateTheta)], 
					    [zc+np.sqrt(b2_scaled)*np.cos(rotateTheta),zc-np.sqrt(b2_scaled)*np.cos(rotateTheta)], 
					    'b') # plot the minor axis
				
		ax.set_xlabel(r'$R\ \mathrm{(kpc)}$')
		ax.set_ylabel(r'$z\ \mathrm{(kpc)}$')
		ax.set_aspect(1)

		return fig

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--subhaloID')
	parser.add_argument('--shape')
	parser.add_argument('--snapNum', default=135)
	args = parser.parse_args()

	subhaloID = args.subhaloID
	subhaloNum = int(subhaloID[7:])
	shape = args.shape
	snapNum = int(args.snapNum)

	os.system('mkdir -p {}/{}'.format(outpath,subhaloID))
	path = outpath + subhaloID + '/'

	fig = plt.figure(figsize = (8,8))
	ax = fig.add_subplot(111)
	fig.subplots_adjust(bottom = 0.1)

	data = getData(illustris_path, snapNum, subhaloNum, 4)
	plotdata = PlotVEMap(data['Coordinates'], data['Velocities'], data['Masses'], shape)
	fig = plotdata.plotVE(fig = fig, ax = ax)
	fig.text(0.5, 0.9, 'subhalo' + str(subhaloNum) + ' (%s)' % shape, 
		horizontalalignment = 'center', fontsize = 20)
	fig.savefig(path + 'VE.eps', dpi = 300)
	plt.close(fig)

if __name__ == '__main__':
	main()


# def plotVE(xcyl, vcyl, mpart, Rb, massweight = True, 
# 	       rmin = 2., N_radius = 10, N_theta = 8, Nmin = 100, s = 1.):
# 	'''
# 	rmin is in units of kpc.
# 	s is the scaling, used to scale the size of the velocity ellipsoids.
# 	'''
# 	id_inRb = xcyl[:,0]**2 + xcyl[:,2]**2 <= Rb**2
# 	xcyl, vcyl = xcyl[id_inRb].copy(), vcyl[id_inRb].copy()
# 	mpart = mpart[id_inRb].copy()
# 	weights = mpart.copy()
# 	if not massweight:
# 		weights = np.ones_like(weights)

# 	theta = np.linspace(10, 80, N_theta)*np.pi/180
# 	radius = np.logspace(np.log10(rmin), np.log10(Rb-0.5), N_radius)

# 	xR = yz = np.arange(0, Rb+1, 0.05)
# 	xR, yz = np.meshgrid(xR, yz)
# 	fig = plt.figure(figsize = (8,8))
# 	ax = fig.add_subplot(111)
# 	fig.subplots_adjust(bottom = 0.1)

# 	for r in radius:
# 		ax.contour(xR, yz, xR**2 + yz**2 - r**2, [0], colors = 'k')
# 		dR = dz = min(70/N_theta*r*np.pi/180, 0.8)

# 		for ang in theta:
# 			Rc = r*np.cos(ang)
# 			zc = r*np.sin(ang)

# 			id_inring = (np.abs(xcyl[:,0]-Rc)<=dR/2) & \
# 			            (np.abs(np.abs(xcyl[:,2])-zc)<=dz/2)

# 			if id_inring.sum()>=Nmin:
				
# 				a11, a12, a13, a22, a23, a33, v1, v2, v3 = \
# 				    calcV2Tensor(vcyl[id_inring], weights[id_inring])
# 				a11 -= v1**2
# 				a12 -= v1*v2
# 				a13 -= v1*v3
# 				a22 -= v2**2
# 				a23 -= v2*v3
# 				a33 -= v3**2

# 				# draw the graph
# 				mat = LA.inv(np.array([[a11,a12,a13],[a12,a22,a23],[a13,a23,a33]]))
# 				a2 = 0.5*(mat[0,0]+mat[2,2]+np.sqrt((mat[0,0]-mat[2,2])**2+4*mat[0,2]**2))
# 				b2 = 0.5*(mat[0,0]+mat[2,2]-np.sqrt((mat[0,0]-mat[2,2])**2+4*mat[0,2]**2))
# 				rotateTheta = np.arccos((mat[0,0]-mat[2,2])/(b2-a2))/2
# 				if mat[0,2]>0:
# 					rotateTheta = np.pi-rotateTheta

# 				b2_scaled = (r*np.pi/2/N_theta*0.4)**2*s
# 				a2_scaled = a2/b2*b2_scaled

# 				ax.contour(xR, yz, \
# 					(xR-Rc)**2*mat[0,0]+(yz-zc)**2*mat[2,2]+2*(xR-Rc)*(yz-zc)*mat[0,2],\
# 					[a2*b2*b2_scaled/b2], colors = 'k')
# 				ax.plot([Rc-np.sqrt(a2_scaled)*np.cos(rotateTheta),Rc+np.sqrt(a2_scaled)*np.cos(rotateTheta)], 
# 					    [zc-np.sqrt(a2_scaled)*np.sin(rotateTheta),zc+np.sqrt(a2_scaled)*np.sin(rotateTheta)], 
# 					    'r') # plot the major axis
# 				ax.plot([Rc-np.sqrt(b2_scaled)*np.sin(rotateTheta),Rc+np.sqrt(b2_scaled)*np.sin(rotateTheta)], 
# 					    [zc+np.sqrt(b2_scaled)*np.cos(rotateTheta),zc-np.sqrt(b2_scaled)*np.cos(rotateTheta)], 
# 					    'b') # plot the minor axis
			
# 	ax.set_xlabel(r'$R\ \mathrm{(kpc)}$')
# 	ax.set_ylabel(r'$z\ \mathrm{(kpc)}$')
# 	ax.set_aspect(1)

# 	return fig




