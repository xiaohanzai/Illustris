import numpy as np
import numpy.linalg  as LA
from Illustris.data.StarData import StarData, Rb_all

def calcV2Tensor(vcyl, weights):
	if len(vcyl) == 0:
		return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
	a11 = np.average(vcyl[:,0]**2, weights=weights)
	a12 = np.average(vcyl[:,0]*vcyl[:,1], weights=weights)
	a13 = np.average(vcyl[:,0]*vcyl[:,2], weights=weights)
	a22 = np.average(vcyl[:,1]**2, weights=weights)
	a23 = np.average(vcyl[:,1]*vcyl[:,2], weights=weights)
	a33 = np.average(vcyl[:,2]**2, weights=weights)
	v1 = np.average(vcyl[:,0], weights=weights)
	v2 = np.average(vcyl[:,1], weights=weights)
	v3 = np.average(vcyl[:,2], weights=weights)
	return a11, a12, a13, a22, a23, a33, v1, v2, v3

class StarData_Beta(StarData):
	'''
	An extension of the StarData class. Added methods to calculate velocity ellipsoids
	  and anisotropy parameter.
	'''
	def __init__(self, x, v, mpart, shape, Rb = Rb_all, **kwargs):
		super(StarData_Beta, self).__init__(x, v, mpart, Rb = Rb, **kwargs)

		# initialize shape
		self.convert2cyl(shape)

	def measureV2Tensor(self, R, z, phi, massweight = True, 
	                    dR = 0.5, dz = 0.5, dphi = 2*np.pi, abs_z = True):
		'''
		Given (R, phi, z), measure the velocity (dispersion) tensor at that point, using 
		  particles in (R-dR/2, R+dR/2) x (phi-dphi/2, phi+dphi/2) x (z-dz/2, z+dz/2).
		If abs_z, then calculate using abs(z) and abs(xcyl[:,2]).
		Return the measured values of VR2, Vz2, Vphi2, and other things.
		R z phi can be vectors. So can dR, dz, dphi.
		'''
		R = np.atleast_1d(R)
		z = np.atleast_1d(z)
		phi = np.atleast_1d(phi)
		if (len(R) != len(z)) | (len(R) != len(phi)) | (len(z) != len(phi)):
			print('R z phi must be all scalars or vectors of the same length.')
			return
		if np.isscalar(dR):
			dR = dR*np.ones_like(R)
		if np.isscalar(dz):
			dz = dz*np.ones_like(z)
		if np.isscalar(dphi):
			dphi = dphi*np.ones_like(phi)
		dR.clip(0.)
		dz.clip(0.)
		dphi.clip(0.,2*np.pi)
		if (len(R) != len(dR)) | (len(phi) != len(dphi)) | (len(z) != len(dz)):
			print('R z phi must be have the same length as dR dz dphi respectively.')
			return

		xcyl = self.xcyl.copy()
		vcyl = self.vcyl.copy()
		mpart = self.mpart.copy()
		if abs_z:
			z = np.abs(z)
			ii = xcyl[:,2] < 0.
			vcyl[ii,2] *= -1.
			xcyl[ii,2] *= -1.
		Rmax = np.max(R + dR/2.)
		Rmin = np.min(R - dR/2.)
		zmax = np.max(z + dz/2.)
		zmin = np.min(z - dz/2.)
		ii = (xcyl[:,0] >= Rmin) & (xcyl[:,0] <= Rmax) & \
		     (xcyl[:,2] >= zmin) & (xcyl[:,2] <= zmax)
		xcyl = xcyl[ii]
		vcyl = vcyl[ii]
		mpart = mpart[ii]
		
		weights = mpart.copy()
		if not massweight:
			weights = np.ones_like(weights)

		a11 = np.zeros(len(R))
		a12 = a11.copy()
		a13 = a11.copy()
		a22 = a11.copy()
		a23 = a11.copy()
		a33 = a11.copy()

		v1 = a11.copy()
		v2 = a11.copy()
		v3 = a11.copy()

		M = a11.copy()
		V = a11.copy()
		Npart = a11.copy()

		for j in range(len(R)):
			ii = (np.abs(xcyl[:,0] - R[j]) < dR[j]/2.) & \
			     (np.abs(xcyl[:,2] - z[j]) < dz[j]/2.) & \
			     ((np.abs(xcyl[:,1] - phi[j]) < dphi[j]/2.) | \
			      (np.abs(xcyl[:,1] - phi[j] - 2*np.pi) < dphi[j]/2.) | \
			      (np.abs(xcyl[:,1] - phi[j] + 2*np.pi) < dphi[j]/2.))
			Npart[j] = ii.sum()

			a11[j],a12[j],a13[j],a22[j],a23[j],a33[j],v1[j],v2[j],v3[j] = \
			    calcV2Tensor(vcyl[ii,:], weights[ii])

			V[j] = dphi[j]*R[j]*dR[j]*dz[j]
			M[j] = mpart[ii].sum()

			xcyl = xcyl[~ii,:]
			vcyl = vcyl[~ii,:]
			mpart = mpart[~ii]
			weights = weights[~ii]

		if len(R) == 1:
			return a11[0], a12[0], a13[0], a22[0], a23[0], a33[0], \
			       v1[0], v2[0], v3[0], M[0], V[0], Npart[0]

		return a11, a12, a13, a22, a23, a33, v1, v2, v3, M, V, Npart

	def calcGlobalAnisotropy(self, Rb = Rb_all, bintype = 'linear', 
	                         N_Rbin = 10, N_phibin = 1, N_zbin = 10, 
	                         Nmin = 50, massweight = True):
		'''
		Calculate the global anisotropy for a galaxy. See the definitions 
		  in section 4.2 of Cappellari et al. 2007.
		Sample the stars inside Rb.
		bintype = 'linear' or 'log'.
		N_*bin are the number of bins on each axis.
		Calculate Beta and set it as an attribute to this class.
		'''
		if bintype not in ['linear', 'log']:
			print('Bintype must be either linear or log.')
			return
		if N_Rbin*N_phibin*N_zbin < 1e-4:
			print('The number of bins must all be at least 1')
			return

		# make bins
		phi = np.linspace(0., 2*np.pi, N_phibin+1)
		dphi = phi[1:] - phi[0:-1]
		phi = (phi[0:-1] + phi[1:])/2.
		if bintype == 'linear':
			R = np.linspace(0., Rb, N_Rbin+1)
			z = np.linspace(0., Rb, N_zbin+1)
		else:
			R = np.logspace(-1.5, np.log10(Rb), N_Rbin+1)
			z = np.logspace(-1.5, np.log10(Rb), N_zbin+1)
		dR = R[1:] - R[0:-1]
		dz = z[1:] - z[0:-1]
		R = (R[1:] + R[0:-1])/2.
		z = (z[1:] + z[0:-1])/2.

		phi, z = np.meshgrid(phi, z)
		phi = phi.reshape(-1)
		z = z.reshape(-1)
		phi = np.meshgrid(phi, R)[0].reshape(-1)
		z, R = np.meshgrid(z, R)
		z = z.reshape(-1)
		R = R.reshape(-1)

		dphi, dz = np.meshgrid(dphi, dz)
		dphi = dphi.reshape(-1)
		dz = dz.reshape(-1)
		dphi = np.meshgrid(dphi, dR)[0].reshape(-1)
		dz, dR = np.meshgrid(dz, dR)
		dz = dz.reshape(-1)
		dR = dR.reshape(-1)

		# calculate V2 tensor
		a11, a12, a13, a22, a23, a33, v1, v2, v3, M, V, Npart = \
		    self.measureV2Tensor(R, z, phi, massweight, dR, dz, dphi, abs_z = True)
		
		# calculate Beta
		ii = Npart > Nmin
		if ii.sum() == 0:
			print('Number of bins not appropriate or Nmin too big. Exiting...')
			return
		Beta = 1 - np.sum((a33[ii] - v3[ii]**2)*M[ii])/np.sum((a11[ii] - v1[ii]**2)*M[ii])
		self.Beta = Beta

	# def calcLambdaR(self, phi=None, inc=None, resolution = 0.5):
	# 	'''
	# 	Calculate lambdaR. phi and theta are the two Euler angles. If they are
	# 	not given, then use self.x_ori to do the calculations.
	# 	phi and theta should be in units of radians.
	# 	'''
	# 	if phi != None:
	# 		xprime = utils.rotateCoordinates(self.x_prin, phi, inc, 0)
	# 		vprime = utils.rotateCoordinates(self.v_prin, phi, inc, 0)
	# 	else:
	# 		xprime = self.x_ori*1.
	# 		vprime = self.v_ori*1.
		
	# 	ii = xprime[:,0]**2 + xprime[:,1]**2 < self.Rb**2
	# 	xprime = xprime[ii,:]
	# 	vprime = vprime[ii,:]
	# 	mpart = self.mpart[ii]*1.

	# 	xbin = ybin = np.arange(-self.Rb, self.Rb, resolution)
	# 	numerator = denominator = 0.

	# 	xtmp = xprime
	# 	vtmp = vprime
	# 	mtmp = mpart

	# 	for j in range(len(xbin)-1):
	# 		ii1 = (xtmp[:,0]>xbin[j]) & (xtmp[:,0]<=xbin[j+1])
	# 		if np.sum(ii1) < 100:
	# 			continue

	# 		xtmp1 = xtmp[ii1,:]
	# 		vtmp1 = vtmp[ii1,:]
	# 		mtmp1 = mtmp[ii1]
	# 		for k in range(len(ybin)-1):
	# 			ii2 = (xtmp1[:,1]>ybin[k]) & (xtmp1[:,1]<=ybin[k+1])
	# 			# print(ii2.sum())
	# 			if np.sum(ii2) < 50:
	# 				continue

	# 			numerator += np.sum(mtmp1[ii2])*\
	# 			             abs(np.mean(vtmp1[ii2,2]))
	# 			denominator += np.sum(mtmp1[ii2])*\
	# 			    np.sqrt(calcDispersion(vtmp1[ii2,2],vtmp1[ii2,2]))

	# 			xtmp1 = xtmp1[~ii2,:]
	# 			vtmp1 = vtmp1[~ii2,:]
	# 			mtmp1 = mtmp1[~ii2]

	# 		xtmp = xtmp[~ii1,:]
	# 		vtmp = vtmp[~ii1,:]
	# 		mtmp = mtmp[~ii1]

	# 	return numerator/denominator

#########################################################################

# def measureV2Tensor(xcyl, vcyl, mpart, R, z, phi, massweight = True, 
# 	                dR = 0.5, dz = 0.5, dphi = 2*np.pi, abs_z = True):
# 	'''
# 	Given (R, phi, z), measure the velocity (dispersion) tensor at that point, using 
# 	  particles in (R-dR/2, R+dR/2) x (phi-dphi/2, phi+dphi/2) x (z-dz/2, z+dz/2).
# 	If abs_z, then calculate using abs(z) and abs(xcyl[:,2]).
# 	R z phi cannot be vectors.
# 	Return the measured values of sigmaR2, sigmaz2, sigmaphi2, and other things.
# 	'''
# 	weights = mpart.copy()
# 	if not massweight:
# 		weights = np.ones_like(weights)

# 	xcyl = xcyl.copy()
# 	if abs_z:
# 		z = np.abs(z)
# 		xcyl[:,2] = np.abs(xcyl[:,2])

# 	ii = (np.abs(xcyl[:,0] - R) < dR/2.) & \
# 	     (np.abs(xcyl[:,2] - z) < dz/2.) & \
# 	     ((np.abs(xcyl[:,1] - phi) < dphi/2.) | \
# 	      (np.abs(xcyl[:,1] - phi - 2*np.pi) < dphi/2.) | \
# 	      (np.abs(xcyl[:,1] - phi + 2*np.pi) < dphi/2.))
# 	Npart = ii.sum()

# 	return calcV2Tensor(vcyl[ii], weights[ii]), Npart

def measureV2map(xcyl, vcyl, mpart, mode = 'linspace', N_Rbin = 20, N_zbin = 20, 
	             Rb = 30, N_phibin = 5, usebin = False, massweight = True):
	'''
	Output V2 vs R phi z. Mode can be set to linspace or logspace.
	If usebin, then use the equal number binning method.
	'''
	xcyl = xcyl.copy()
	vcyl = vcyl.copy()
	ii = xcyl[:,2] < 0.
	vcyl[ii,2] *= -1.
	xcyl[ii,2] *= -1.
	mpart = mpart.copy()
	weights = mpart.copy()
	if not massweight:
		weights = np.ones_like(weights)

	if not usebin:
		if mode == 'linspace':
			R_bin = np.linspace(0, Rb, N_Rbin+1)
			z_bin = np.linspace(0, Rb, N_zbin+1)
		elif mode == 'logspace':
			R_bin = np.logspace(-1, np.log10(Rb), N_Rbin+1)
			z_bin = np.logspace(-1, np,log10(Rb), N_zbin+1)
		z_bin = np.array([z_bin for i in range(len(R_bin)-1)])
	else:
		R_bin, z_bin = \
		    equalNumBin.bin2D(xcyl[:,0],xcyl[:,2],N_Rbin,N_zbin,0,0,Rb,Rb)

	phi_bin = np.linspace(0,2*np.pi,N_phibin+1)

	a11 = np.zeros(N_Rbin*N_zbin*N_phibin)
	a12 = a11.copy()
	a13 = a11.copy()
	a22 = a11.copy()
	a23 = a11.copy()
	a33 = a11.copy()

	v1 = a11.copy()
	v2 = a11.copy()
	v3 = a11.copy()

	nu = a11.copy()

	NperBin = a11.copy()

	R_bin0 = a11.copy()
	z_bin0 = a11.copy()
	phi_bin0 = a11.copy()

	n = 0
	for j in range(N_Rbin):
		ii1 = (xcyl[:,0]>R_bin[j]) & (xcyl[:,0]<=R_bin[j+1])

		xcyl_tmp = xcyl[ii1,:]
		vcyl_tmp = vcyl[ii1,:]
		m_tmp = mpart[ii1]
		weights_tmp = weights[ii1]

		for k in range(N_zbin):
			ii2 = (xcyl_tmp[:,2]>z_bin[j,k]) & (xcyl_tmp[:,2]<z_bin[j,k+1])

			xcyl_tmp_tmp = xcyl_tmp[ii2,:]
			vcyl_tmp_tmp = vcyl_tmp[ii2,:]
			m_tmp_tmp = m_tmp[ii2]
			weights_tmp_tmp = weights_tmp[ii2]

			for l in range(N_phibin):
				ii3 = (xcyl_tmp_tmp[:,1]>phi_bin[l]) & (xcyl_tmp_tmp[:,1]<phi_bin[l+1])

				a11[n],a12[n],a13[n],a22[n],a23[n],a33[n],v1[n],v2[n],v3[n] = \
				    calcV2Tensor(vcyl_tmp_tmp[ii3,:], weights_tmp_tmp[ii3])

				NperBin[n] = ii3.sum()

				V = np.pi*(R_bin[j+1]**2 - R_bin[j]**2)*(z_bin[j,k+1] - z_bin[j,k])*1e9 / N_phibin
				nu[n] = m_tmp_tmp[ii3].sum()/2./V

				xcyl_tmp_tmp = xcyl_tmp_tmp[~ii3,:]
				vcyl_tmp_tmp = vcyl_tmp_tmp[~ii3,:]
				m_tmp_tmp = m_tmp_tmp[~ii3]
				weights_tmp_tmp = weights_tmp_tmp[~ii3]

				R_bin0[n] = (R_bin[j] + R_bin[j+1]) / 2
				z_bin0[n] = (z_bin[j,k] + z_bin[j,k+1]) / 2
				phi_bin0[n] = (phi_bin[l] + phi_bin[l+1]) / 2

				n += 1

			xcyl_tmp = xcyl_tmp[~ii2,:]
			vcyl_tmp = vcyl_tmp[~ii2,:]
			m_tmp = m_tmp[~ii2]
			weights_tmp = weights_tmp[~ii2]

		xcyl = xcyl[~ii1,:]
		vcyl = vcyl[~ii1,:]
		mpart = mpart[~ii1]
		weights = weights[~ii1]

	return R_bin0, z_bin0, phi_bin0, a11, a12, a13, a22, a23, a33, v1, v2, v3, nu, NperBin

# def calcGlobalAnisotropy(xcyl, vcyl, mpart, Rb, bintype = 'linear', 
# 	                     N_Rbin = 10, N_phibin = 1, N_zbin = 10, 
# 	                     Nmin = 50, massweight = True):
# 	'''
# 	Calculate the global anisotropy for a galaxy. See the definitions 
# 	  in section 4.2 of Cappellari et al. 2007.
# 	The input should be in cylindrical coordinates.
# 	Sample the stars inside Rb.
# 	bintype = 'linear' or 'log'.
# 	N_*bin are the number of bins in each axis.
# 	z should be the symmetry axis.
# 	Return Beta, Gamma, Delta.
# 	'''
# 	if bintype not in ['linear', 'log']:
# 		raise ValueError("bintype must be linear or log.")
# 	if N_Rbin*N_phibin*N_zbin < 1e-4:
# 		raise RuntimeError("The number of bins must all be at least 1")

# 	id_inRb = xcyl[:,0]**2 + xcyl[:,2]**2 <= Rb**2
# 	xcyl = xcyl[id_inRb].copy()
# 	xcyl[:,2] = np.abs(xcyl[:,2])
# 	vcyl = vcyl[id_inRb].copy()
# 	mpart = mpart[id_inRb].copy()

# 	# make bins
# 	phi_bin = np.linspace(0., 2*np.pi, N_phibin+1)
# 	if bintype == 'linear':
# 		R_bin = np.linspace(0., Rb, N_Rbin+1)
# 		z_bin = np.linspace(0., Rb, N_zbin+1)
# 	else:
# 		R_bin = np.logspace(-1.5, np.log10(Rb), N_Rbin+1)
# 		z_bin = np.logspace(-1.5, np.log10(Rb), N_zbin+1)

# 	# weights
# 	weights = mpart.copy()
# 	if not massweight:
# 		weights = np.ones_like(weights)

# 	#Calculate PIRR, PIphiphi, PIzz
# 	PIRR = PIphiphi = PIzz = 0.
# 	for j in range(N_Rbin):
# 		ii1 = (xcyl[:,0]>R_bin[j]) & (xcyl[:,0]<=R_bin[j+1])
# 		if ii1.sum() < Nmin:
# 			continue
		
# 		for k in range(N_zbin):
# 			ii2 = (xcyl[ii1,2]>z_bin[k]) & (xcyl[ii1,2]<=z_bin[k+1])
# 			if ii2.sum() < Nmin:
# 				continue

# 			for l in range(N_phibin):
# 				ii3 = (xcyl[ii1][ii2,1]>phi_bin[l]) & (xcyl[ii1][ii2,1]<=phi_bin[l+1])
# 				if ii3.sum() < Nmin:
# 					continue

# 				a11, a12, a13, a22, a23, a33, v1, v2, v3 = \
# 					calcV2Tensor(vcyl[ii1][ii2][ii3], weights[ii1][ii2][ii3])

# 				M = np.sum(mpart[ii1][ii2][ii3])
# 				PIRR += M*(a11 - v1**2)
# 				PIphiphi += M*(a22 - v2**2)
# 				PIzz += M*(a33 - v3**2)

# 	Beta = 1 - PIzz/PIRR
# 	Gamma = 1 - PIphiphi/PIRR
# 	Delta = (2*Beta-Gamma)/(2-Gamma)

# 	return Beta, Gamma, Delta
