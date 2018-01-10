'''
Some useful utilities.
'''
import numpy as np

def xyz2cyl(x, v):
	'''
	Convert from (x,y,z) to (R,phi,z).
	z axis should be the symmetry axis.
	'''
	R = np.sqrt(x[:,0]**2 + x[:,1]**2)
	phi = np.arccos(x[:,0]/R)
	id_phi = np.where(x[:,1]<0)[0]
	phi[id_phi] = 2*np.pi - phi[id_phi]
	
	xcyl = np.zeros((len(x),3))
	xcyl[:,0] = R; xcyl[:,1] = phi; xcyl[:,2] = x[:,2]+0.

	vR = v[:,0]*x[:,0]/xcyl[:,0] + v[:,1]*x[:,1]/xcyl[:,0]
	vphi = -v[:,0]*x[:,1]/xcyl[:,0] + v[:,1]*x[:,0]/xcyl[:,0]
	
	vcyl = np.zeros((len(x),3))
	vcyl[:,0] = vR; vcyl[:,1] = vphi; vcyl[:,2] = v[:,2]+0.

	return xcyl, vcyl

def xyz2sph(x, v):
	'''
	Convert from (x,y,z) to (r,phi,theta).
	'''
	r = np.sqrt(np.sum(x**2, axis=1))
	R = np.sqrt(x[:,0]**2 + x[:,1]**2)
	theta = np.arccos(x[:,2]/r)
	phi = np.arccos(x[:,0]/R)
	id_phi = np.where(x[:,1]<0)[0]
	phi[id_phi] = 2*np.pi - phi[id_phi]

	xsph = np.zeros_like(x)
	xsph[:,0] = r; xsph[:,1] = phi; xsph[:,2] = theta

	vR = v[:,0]*x[:,0]/R + v[:,1]*x[:,1]/R
	vphi = -v[:,0]*x[:,1]/R + v[:,1]*x[:,0]/R
	vr = vR*np.sin(theta) + v[:,2]*np.cos(theta)
	vtheta = -vR*np.cos(theta) + v[:,2]*np.sin(theta)

	vsph = np.zeros_like(v)
	vsph[:,0] = vr; vsph[:,1] = vphi; vsph[:,2] = vtheta

	return xsph, vsph

def rotateCoordinates(x, phi, inc, pa):
	'''
	Rotate the coordinates x. See Monnet et al.1992.
	The angles should all be in units of radians.
	'''
	R1 = np.array([[np.cos(phi),  np.sin(phi), 0],
		           [-np.sin(phi), np.cos(phi), 0],
		           [    0,            0,       1]])
	R2 = np.array([[1,      0,            0      ],
		           [0,  np.cos(inc),  np.sin(inc)],
		           [0,  -np.sin(inc), np.cos(inc)]])
	R3 = np.array([[np.cos(pa),  np.sin(pa), 0],
		           [-np.sin(pa), np.cos(pa), 0],
		           [     0,          0,      1]])

	xprime = np.dot(np.dot(R3, np.dot(R2, R1)), x.T)
	return xprime.T

def generateCircle(xc, yc, r):
	theta = np.linspace(0, 2*np.pi, 100)
	return xc+r*np.cos(theta), yc+r*np.sin(theta)




