'''
Get the measured V2 tensors at a number of grid points and store them in a .npy file.
'''
import numpy as np
import os
import argparse
from Illustris.beta.StarData_Beta import StarData_Beta, measureV2map
from Illustris.utils.util_illustris import getData
import Illustris.utils.paths as paths
from illustris_python.groupcat import loadSingle

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--subhaloID')
	parser.add_argument('--shape')
	parser.add_argument('--snapNum')
	parser.add_argument('--TNG', default = False, type = bool)
	parser.add_argument('--outpath', default = None)
	args = parser.parse_args()

	subhaloID = args.subhaloID
	subhaloNum = int(subhaloID[7:])
	shape = args.shape
	snapNum = int(args.snapNum)

	outpath = args.outpath
	if outpath is None:
		if args.TNG:
			outpath = paths.TNG_samplepath
		else:
			outpath = paths.illustris_samplepath

	illustris_path = paths.illustris_path
	if args.TNG:
		illustris_path = paths.TNG_path

	os.system('mkdir -p {}/{}'.format(outpath,subhaloID))
	path = outpath + subhaloID + '/'

	data = getData(illustris_path, snapNum, subhaloNum, 4)
	V2data = StarData_Beta(data['Coordinates'], data['Velocities'], data['Masses'], shape, xc=data['SubhaloPos'])
	# rmax = 30
	rmax = data['HalfMassRad']*2.5
	# R_bin, z_bin, phi_bin, a11, a12, a13, a22, a23, a33, v1, v2, v3, M, V, Npart = \
	#     measureV2map(V2data.xcyl, V2data.vcyl, V2data.mpartz)
	# beta = 1 - a33/a11
	# np.save(path + 'V2tensor.npy', 
	#          [R_bin, z_bin, phi_bin, a11, a12, a13, a22, a23, a33, v1, v2, v3, beta, M, V, Npart])
	R_bin, z_bin, phi_bin, a11, a12, a13, a22, a23, a33, v1, v2, v3, M, V, Npart = \
	    measureV2map(V2data.xcyl, V2data.vcyl, V2data.mpart, N_phibin = 1, Rb = rmax)
	beta = 1 - a33/a11
	f = open(path + 'V2tensor.txt', 'w')
	for i in range(len(R_bin)):
		f.write('%.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2e  %.2e  %d \n'
			% (R_bin[i], z_bin[i], a11[i], a12[i], a13[i], a22[i], a23[i], a33[i], 
				v1[i], v2[i], v3[i], beta[i], M[i], V[i], Npart[i]))
	f.close()

if __name__ == '__main__':
	main()



# '''
# Get the measured V2 tensors at a number of grid points and store them in a .npy file.
# '''
# import sys
# sys.path.append('../utils/')
# import numpy as np
# import os
# import argparse
# from StarData_Beta import StarData_Beta
# import equalNumBin
# from util_illustris import getData
# from paths import illustris_samplepath, illustris_path
# outpath = illustris_samplepath

# ###########################################

# def main():
# 	parser = argparse.ArgumentParser()
# 	parser.add_argument('--subhaloID')
# 	parser.add_argument('--shape')
# 	parser.add_argument('--snapNum', default=135)
# 	args = parser.parse_args()

# 	subhaloID = args.subhaloID
# 	subhaloNum = int(subhaloID[7:])
# 	shape = args.shape
# 	snapNum = int(args.snapNum)

# 	os.system('mkdir -p {}/{}'.format(outpath,subhaloID))
# 	path = outpath + subhaloID + '/'

# 	# make bins
# 	N_Rbin = 20
# 	N_zbin = 20
# 	N_phibin = 5
# 	usebin = False
# 	mode = 'linspace'
# 	Rb = 30

# 	# make bins
# 	phi = np.linspace(0., 2*np.pi, N_phibin+1)
# 	dphi = phi[1:] - phi[0:-1]
# 	phi = (phi[0:-1] + phi[1:])/2.
# 	if not usebin:
# 		if mode == 'linspace':
# 			R = np.linspace(0., Rb, N_Rbin+1)
# 			z = np.linspace(0., Rb, N_zbin+1)
# 		elif mode == 'logspace':
# 			R = np.logspace(-1, np.log10(Rb), N_Rbin+1)
# 			z = np.logspace(-1, np,log10(Rb), N_zbin+1)
# 		z = np.array([z for i in range(len(R)-1)])
# 	else:
# 		R, z = \
# 		    equalNumBin.bin2D(xcyl[:,0],xcyl[:,2],N_Rbin,N_zbin,0,0,Rb,Rb)
# 	dR = R[1:] - R[0:-1]
# 	dz = z[1:] - z[0:-1]
# 	R = (R[1:] + R[0:-1])/2.
# 	z = (z[1:] + z[0:-1])/2.

# 	phi, z = np.meshgrid(phi, z)
# 	phi = phi.reshape(-1)
# 	z = z.reshape(-1)
# 	phi = np.meshgrid(phi, R)[0].reshape(-1)
# 	z, R = np.meshgrid(z, R)
# 	z = z.reshape(-1)
# 	R = R.reshape(-1)

# 	dphi, dz = np.meshgrid(dphi, dz)
# 	dphi = dphi.reshape(-1)
# 	dz = dz.reshape(-1)
# 	dphi = np.meshgrid(dphi, dR)[0].reshape(-1)
# 	dz, dR = np.meshgrid(dz, dR)
# 	dz = dz.reshape(-1)
# 	dR = dR.reshape(-1)

# 	# calculate V2 tensors
# 	data = getData(illustris_path, snapNum, subhaloNum, 4)
# 	V2data = StarData_Beta(data['Coordinates'], data['Velocities'], data['Masses'], shape)
# 	a11, a12, a13, a22, a23, a33, v1, v2, v3, M, V, Npart = \
# 	    V2data.measureV2Tensor(R, z, phi, massweight = True, 
# 	                    dR = dR, dz = dz, dphi = dphi, abs_z = True)
# 	nu = M/V
#	beta = 1 - a33/a11

# 	np.save(path + 'V2tensor.npy', 
# 		   [R, z, phi, a11, a12, a13, a22, a23, a33, beta, nu, Npart])

# if __name__ == '__main__':
# 	main()
