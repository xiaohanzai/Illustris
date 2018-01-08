import numpy as np
import Illustris.utils.util_illustris as ui
import Illustris.utils.paths as paths
import illustris_python as il
import argparse
import os
import pickle
illustris_path = paths.illustris_path
outpath = paths.illustris_samplepath

def shape(R, xpart):
    '''
    Given a series of values of R, calculate the axis ratios and eigenvectors
      for particles inside each R.
    '''
    axisRatios = np.zeros([len(R), 2])
    eigenVectors = np.zeros([len(R), 3, 3])

    for i in range(len(R)):
        ba, ca, tmp, tmp, Tiv = ui.calcShape(xpart, Rb=R[i])
        if np.dot(Tiv[2, :], np.cross(Tiv[0, :], Tiv[1, :])) < 0: # right hand coordinate
            Tiv[1, :] *= -1
        for j in range(3):
            if Tiv[j, 2] < 0:
                Tiv[j, :] *= -1

        axisRatios[i, :] = np.array([ba, ca])
        eigenVectors[i, :, :] = Tiv # eigen vectors Tiv[0,:], Tiv[1,:], Tiv[2,:]

    return axisRatios, eigenVectors


def run(xpart_star, mpart_star, hmr_star, xpart_dark, mpart_dark, hmr_dark, path, 
        Nbin_star = 30, Nbin_dark = 100):
    '''
    Save relavent files under "path".
    '''
    xpart = np.concatenate((xpart_star, xpart_dark))
    mpart = np.concatenate((mpart_star, mpart_dark))
    xc = ui.findCenter(xpart, mpart=mpart)

    Rstar = np.linspace(3.0, 2.5*hmr_star, Nbin_star)
    Rdark = np.linspace(3.0, 2.5*hmr_dark, Nbin_dark)
    axisRatiosStar, eigenVectorsStar = shape(Rstar, xpart_star-xc)
    axisRatiosDark, eigenVectorsDark = shape(Rdark, xpart_dark-xc)
    rst = {'Rstar': Rstar, 'hmr_star': hmr_star, 'axisRatiosStar':
           axisRatiosStar, 'eigenVectorsStar': eigenVectorsStar,
           'Rdark': Rdark, 'hmr_dark': hmr_dark, 'axisRatiosDark':
           axisRatiosDark, 'eigenVectorsDark': eigenVectorsDark}
    with open('{}/shape.dat'.format(path), 'wb') as f:
        pickle.dump(rst, f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subhaloID')
    parser.add_argument('--snapNum')
    args = parser.parse_args()

    subhaloID = args.subhaloID
    subhaloNum = int(subhaloID[7:])
    snapNum = int(args.snapNum)

    foldername = subhaloID
    os.system('mkdir -p {}/{}'.format(outpath, foldername))
    path = outpath + foldername
    os.system('mkdir -p {}/shape'.format(path))

    subhalo = il.groupcat.loadSingle(illustris_path, snapNum, subhaloID=subhaloNum)
    # z = ui.snap2z(snapNum)
    hmr_star = subhalo['SubhaloHalfmassRadType'][4] / ui.h0 # / (1 + z)
    hmr_dark = subhalo['SubhaloHalfmassRadType'][1] / ui.h0 # / (1 + z)

    data_star = ui.getData(illustris_path, snapNum, subhaloNum, 4)
    data_dark = ui.getData(illustris_path, snapNum, subhaloNum, 1)
    
    run(data_star['Coordinates'], data_star['Masses'], hmr_star, 
        data_dark['Coordinates'], data_dark['Masses'], hmr_dark, path+'/shape/')


if __name__ == '__main__':
    main()
