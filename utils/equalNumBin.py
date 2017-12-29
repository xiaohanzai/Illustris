'''
Try to make grids that contain equal number of particles inside each cell.
Basicly aimed at binning illustris galaxies. 2D grids of (R,z).
Although this is not entirely doable... We have to make a compromise.
'''
import numpy as np

def bin1D(x1, Nbin1, start1, end1):
	'''
	x1 is the coordinate vector to be binned.
	Nbin1 is the number of bins you want to make.
	start1 and end1 specify the ranges of points to be binned.
	'''
	ii = (x1>start1) & (x1<end1)
	x1 = x1[ii]

	# These two vectors will tell you where the bin boundaries are.
	x1_bin = np.zeros(Nbin1+1)+start1
	x1_bin[-1] = end1

	# # These two vectors are dictionaries. If the ith element of x1 fall into
	# # the nth bin, then x1_idbin[i] = n. n starts from 0.
	# x1_idbin = np.zeros(len(x1)) + Nbin1

	# eps is the tolerance: the deviation of the number of particles in a bin 
	# from perfect equal-number binning. In units of number.
	eps1 = max(len(x1)/Nbin1**2, 1)

	N1 = len(x1)/Nbin1

	########## begin calculations ##########
	# bin x1
	left = start1
	right = end1
	for j in range(Nbin1-1):
		# first guess a boundary at tmp_bd
		tmp_bd = left + (right - left) / (Nbin1-j)

		ii = (x1>left) & (x1<tmp_bd)
		N = ii.sum()

		while abs(N - N1) > eps1:
			if N - N1 > eps1: # shift tmp_bd to the left
				right = tmp_bd
				tmp_bd = (tmp_bd + left) / 2.
			else:            # shift tmp_bd to the right
				tmp_bd = (tmp_bd + right) / 2.

			ii = (x1>left) & (x1<tmp_bd)
			N = ii.sum()

		left = tmp_bd
		right = end1
		x1_bin[1+j] = tmp_bd
		#x1_idbin[ii] = j

	return x1_bin#, x1_idbin

def bin2D(x1, x2, Nbin1, Nbin2, start1, start2, end1, end2):
	x1_bin = bin1D(x1, Nbin1, start1, end1)
	x2_bin = np.zeros((Nbin1, Nbin2+1))

	for j in range(Nbin1):
		ii = (x1>x1_bin[j]) & (x1<x1_bin[j+1])
		tmp2_bin = bin1D(x2[ii], Nbin2, start2, end2)
		x2_bin[j] = tmp2_bin

	return x1_bin, x2_bin

def main():
	a = np.random.rand(10000)
	x1_bin = bin1D(a, 10, 0, 1)

	b = np.random.randn(10000)
	x2_bin = bin1D(b, 10, -2, 2)

	print(x1_bin)
	print(x2_bin)

	ii = (b>x2_bin[0]) & (b<x2_bin[1])
	print(ii.sum())
	ii = (b>x2_bin[4]) & (b<x2_bin[5])
	print(ii.sum())
	ii = (b>x2_bin[-2]) & (b<x2_bin[-1])
	print(ii.sum())

	x1_bin, x2_bin = bin2D(b, b, 10, 10, -2, -2, 2, 2)
	print(x1_bin)
	print(x2_bin)

if __name__ == '__main__':
	main()
