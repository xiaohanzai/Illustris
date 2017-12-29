'''
Plot beta distributions in the equatorial plane and the meridional plane.
'''
import sys
sys.path.append('../utils/')
import numpy as np
import os
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
rc('xtick', labelsize=10)
rc('ytick', labelsize=10)
rc('axes', labelsize=15)
from StarData_Beta import StarData_Beta, Rb_all
from util_illustris import getData
from paths import illustris_samplepath, illustris_path
outpath = illustris_samplepath

colormap = plt.get_cmap('jet')

def _vminvmax(c):
	'''
	c can either be 1D or 2D.
	Determine the proper vmin and vmax by keeping only the middle 96% points.
	'''
	c = c.reshape(-1)
	c = c[~np.isnan(c)]
	c.sort()
	tmp = max(int(0.02*c.size),1)
	vmin = c[tmp]
	vmax = c[c.size-tmp]
	return vmin, vmax

class PlotBetaMap(StarData_Beta):
	def __init__(self, x, v, mpart, shape, Rb = Rb_all, **kwargs):
		super(PlotBetaMap, self).__init__(x, v, mpart, shape, Rb = Rb_all, **kwargs)
		self.calcGlobalAnisotropy()

	def map_equ(self, rmax, massweight = True, 
		        dR = 1., dz = 4., N_phi = 6, s = 30, Nmin = 50, 
		        fig = None, ax = None, cax = None):
		'''
		Plot beta map in the equatorial plane.
		  dR:  step for binning in R
		  dz: thinkness of the equatorial plane
		  N_phi: number of bins in phi, for the innermost ring.
		  s: size of the points.
		'''
		R_bin = np.arange(dR, rmax - dR, dR)

		R = np.zeros(N_phi*(1+len(R_bin))*len(R_bin)//2)
		phi = R.copy()
		dphi = R.copy()

		# make bins
		N = 0
		for j in range(len(R_bin)):
			N += j*N_phi
			N_tmp = N_phi*(j+1)

			step = 2*np.pi/N_tmp
			dphi[N:N+N_tmp] = step
			phi[N:N+N_tmp] = np.linspace(step/2., 2*np.pi-step/2., N_tmp)
			R[N:N+N_tmp] = R_bin[j]		

		# calculate beta
		a11, a12, a13, a22, a23, a33, v1, v2, v3, M, V, Npart = \
			self.measureV2Tensor(R, np.zeros_like(R), phi, massweight = massweight, 
		        dR = dR, dz = dz, dphi = dphi, abs_z = False)
	
		# we still use velocity dispersion to calculate beta
		beta = 1 - (a33 - v3**2)/(a11 - v1**2)
		beta[Npart < Nmin] = np.nan
		
		# make plots
		if fig == None:
			fig = plt.figure()
			ax = fig.add_axes([0.05,0.1,0.85,0.9], projection='polar')
			cax = fig.add_axes([0.87, 0.1, 0.02, 0.9])
		vmin, vmax = _vminvmax(beta)
		norm = mpl.colors.Normalize(vmin, vmax)
		sc = ax.scatter(phi, R, c = beta, s = s, cmap = colormap,
				        norm = norm, linewidths = 0)
		ax.set_xticks([])
		ax.set_ylim(0, rmax)
		cbar = fig.colorbar(sc, cax = cax)

		return fig

	def map_meri(self, rmax, massweight = True, 
		         N_bin = 20, Nmin = 100, 
		         fig = None, ax = None, cax = None):
		'''
		Plot beta map in the meridional plane.
		  N_bin: number of bins in R and z.
		'''
		# make bins
		step = rmax/N_bin
		R_bin = z_bin = np.linspace(step/2., rmax - step/2., N_bin)
		z, R = np.meshgrid(z_bin, R_bin)
		R = R.reshape(-1)
		z = z.reshape(-1)

		# calculate beta
		a11, a12, a13, a22, a23, a33, v1, v2, v3, M, V, Npart = \
			self.measureV2Tensor(R, z, np.zeros_like(R), massweight = massweight, 
		        dR = step, dz = step, dphi = 2*np.pi, abs_z = True)
	
		# we still use velocity dispersion to calculate beta
		beta = 1 - (a33 - v3**2)/(a11 - v1**2)
		beta[Npart < Nmin] = np.nan

		img = np.zeros([N_bin, N_bin])
		for j in range(len(R_bin)):
			for k in range(len(z_bin)):
				img[N_bin-1-k,j] = beta[j*N_bin+k]
		
		# make plots
		if fig == None:
			fig = plt.figure()
			ax = fig.add_axes([0.05,0.1,0.85,0.9])
			cax = fig.add_axes([0.87, 0.1, 0.02, 0.9])
		vmin, vmax = _vminvmax(img)
		im = ax.imshow(img, extent = (0, rmax, 0, rmax), 
			           vmin = vmin, vmax = vmax, 
			           interpolation = 'nearest', cmap = colormap,)
		ax.text(0.6*rmax, 1.02*rmax, r'$\beta = %.2f$' % self.Beta, 
			     verticalalignment = 'bottom', fontsize = 15)
		ax.set_aspect(1)
		ax.set_xlabel(r'$R\ \mathrm{(kpc)}$')
		ax.set_ylabel(r'$z\ \mathrm{(kpc)}$')
		cbar = fig.colorbar(im, cax = cax)

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

	fig = plt.figure(figsize = (8,4))
	ax1 = fig.add_axes([0.02, 0.08, 0.35, 0.7], projection='polar')
	cax1 = fig.add_axes([0.4, 0.12, 0.02, 0.65])
	ax2 = fig.add_axes([0.39, 0.12, 0.65, 0.65])
	cax2 = fig.add_axes([0.9, 0.12, 0.02, 0.65])

	data = getData(illustris_path, snapNum, subhaloNum, 4)
	plotdata = PlotBetaMap(data['Coordinates'], data['Velocities'], data['Masses'], shape)
	fig = plotdata.map_equ(20., fig = fig, ax = ax1, cax = cax1)
	fig = plotdata.map_meri(20., fig = fig, ax = ax2, cax = cax2)
	fig.text(0.5, 0.9, 'subhalo' + str(subhaloNum) + ' (%s)' % shape,
	        horizontalalignment = 'center', fontsize = 20)
	fig.savefig(path + 'beta.eps', dpi=300)
	plt.close(fig)

if __name__ == '__main__':
	main()


# def map_equVSmeri(xcyl, vcyl, mpart, Rb, Beta = -99, massweight = True, 
# 	              dR = 1., dz = 4., N_phi = 6, s = 30, Nmin1 = 50, 
# 	              N_bin = 20, Nmin2 = 100):
# 	'''
# 	LHS: beta map in the equatorial plane.
# 	  dR:  step for binning in R
# 	  dz: thinkness of the equatorial plane
# 	  N_phi: number of bins in phi, for the innermost ring.
# 	  s: size of the points
# 	RHS: beta map in the meridional plane.
# 	  N_bin: number of bins in R and z.
# 	'''
# 	id_inRb = xcyl[:,0]**2 + xcyl[:,2]**2 <= Rb**2
# 	xcyl, vcyl = xcyl[id_inRb].copy(), vcyl[id_inRb].copy()
# 	mpart = mpart[id_inRb].copy()
# 	weights = mpart.copy()
# 	if not massweight:
# 		weights = np.ones_like(weights)

# 	if Beta < -10:
# 		Beta = calcGlobalAnisotropy(xcyl, vcyl, mpart, Rb)[0]

# 	fig = plt.figure(figsize = (8,4))

# 	###### RHS ######
# 	step = Rb/N_bin
# 	R_bin = z_bin = np.linspace(step/2., Rb-step/2., N_bin)
# 	img = np.zeros([N_bin, N_bin])

# 	for j in range(len(R_bin)):
# 		ii1 = np.abs(xcyl[:,0]-R_bin[j])<=step/2.
# 		for k in range(len(z_bin)):
# 			ii2 = np.abs(np.abs(xcyl[ii1,2])-z_bin[k])<=step/2.
# 			if ii2.sum() < Nmin2:
# 				img[N_bin-1-k,j] = np.nan
# 				continue

# 			a11_tmp, a12_tmp, a13_tmp, a22_tmp, a23_tmp, a33_tmp, v1, v2, v3 = \
# 				calcV2Tensor(vcyl[ii1][ii2], weights[ii1][ii2])

# 			img[N_bin-1-k,j] = 1-(a33_tmp - v3**2)/(a11_tmp - v1**2)
	
# 	ax2 = fig.add_axes([0.39, 0.12, 0.65, 0.65])
# 	vmin, vmax = _vminvmax(img)
# 	im = ax2.imshow(img, extent = (0, Rb, 0, Rb), 
# 		            vmin = vmin, vmax = vmax, 
# 		            interpolation = 'nearest', cmap = colormap,)
# 	ax2.text(0.6*Rb, 1.02*Rb, r'$\beta = %.2f$' % Beta, 
# 		     verticalalignment = 'bottom', fontsize = 15)
# 	ax2.set_aspect(1)
# 	ax2.set_xlabel(r'$R\ \mathrm{(kpc)}$')
# 	ax2.set_ylabel(r'$z\ \mathrm{(kpc)}$')
# 	cax2 = fig.add_axes([0.9, 0.12, 0.02, 0.65])
# 	cbar2 = fig.colorbar(im, cax = cax2)
# 	# cax2 = fig.add_axes([0.53, 0.7, 0.35, 0.02])
# 	# cbar2 = fig.colorbar(im, cax = cax2, orientation = 'horizontal')

# 	###### LHS ######
# 	id_equ = np.abs(xcyl[:,2])<dz/2
# 	xcyl, vcyl = xcyl[id_equ], vcyl[id_equ]
# 	weights = weights[id_equ]

# 	R_bin = np.arange(dR, Rb-dR, dR)

# 	a11 = np.zeros(N_phi*(1+len(R_bin))*len(R_bin)//2)
# 	a33 = a11.copy()
# 	xval1 = a11.copy()
# 	yval1 = a11.copy()

# 	for j in range(len(R_bin)):
# 		N_tmp = N_phi*(j+1)
# 		dphi = 2*np.pi/N_phi/R_bin[j]*dR
		
# 		N = j*(j+1)//2*N_phi
# 		xval1[N:N+N_tmp] = \
# 		     np.linspace(np.pi/N_tmp, 2*np.pi-np.pi/N_tmp, N_tmp)
# 		yval1[N:N+N_tmp] = R_bin[j]
		
# 		ii1 = np.abs(xcyl[:,0]-R_bin[j])<=dR/2

# 		for k in range(N_tmp):
# 			ii2 = np.abs(xcyl[ii1,1]-dphi*(k+0.5))<dphi/2
# 			a11_tmp, a12_tmp, a13_tmp, a22_tmp, a23_tmp, a33_tmp, v1, v2, v3 = \
# 				calcV2Tensor(vcyl[ii1][ii2], weights[ii1][ii2])
		
# 			a11[N+k] = (a11_tmp - v1**2)
# 			a33[N+k] = (a33_tmp - v3**2)

# 	beta = 1-a33/a11
	
# 	ax1 = fig.add_axes([0.02, 0.08, 0.35, 0.7], projection='polar')
# 	vmin, vmax = _vminvmax(beta)
# 	norm = mpl.colors.Normalize(vmin, vmax)
# 	# norm = mpl.colors.Normalize(np.min(beta[~np.isnan(beta)]), 
# 		                        # np.max(beta[~np.isnan(beta)]))
# 	sc1 = ax1.scatter(xval1, yval1, c = beta, s = s, cmap = colormap,
# 			        norm = norm, linewidths = 0)
# 	ax1.set_xticks([])
# 	ax1.set_ylim(0, Rb)
# 	# ax1.set_ylabel('$R$ (kpc)')
# 	cax1 = fig.add_axes([0.4, 0.12, 0.02, 0.65])
# 	cbar1 = fig.colorbar(sc1, cax = cax1)

# 	return fig




