import numpy as np
import Illustris.shape 
import Illustris.utils.util_illustris as ui
import Illustris.utils.paths as paths
import illustris_python as il
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subhaloID')
    parser.add_argument('--snapNum')
    parser.add_argument('--TNG', default = False, type = bool)
    parser.add_argument('--outpath', default = None)
    args = parser.parse_args()

    subhaloID = args.subhaloID
    subhaloNum = int(subhaloID[7:])
    snapNum = int(args.snapNum)

    illustris_path = paths.illustris_path
    if args.TNG:
        illustris_path = paths.TNG_path

    outpath = args.outpath
    if outpath is None:
        if args.TNG:
            outpath = paths.TNG_samplepath
        else:
            outpath = paths.illustris_samplepath

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
    
    Illustris.shape.getShape.run(data_star['Coordinates'] - data_star['SubhaloPos'], data_star['Masses'], hmr_star, 
        data_dark['Coordinates'] - data_dark['SubhaloPos'], data_dark['Masses'], hmr_dark, path+'/shape/')
    Illustris.shape.plotShape.run(path+'/shape/')

if __name__ == '__main__':
	main()
