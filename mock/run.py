import numpy as np
import argparse
import os
import Illustris.mock
import Illustris.utils.util_illustris as ui
import Illustris.utils.paths as paths

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subhaloID')
    parser.add_argument('--snapNum')
    parser.add_argument('--rotation') # oblate or prolate rotation
    parser.add_argument('--phi') # let's input in units of degrees
    parser.add_argument('--inc')
    parser.add_argument('--pa', default = 500)
    parser.add_argument('--TNG', default = False, type = bool)
    parser.add_argument('--outpath', default = None)
    parser.add_argument('-i', action='store_false', dest='iteration',
                      default=True, help='Iteration fit')
    parser.add_argument('-c', action='store_true', dest='circle',
                      default=False, help='use circlular aperture')
    args = parser.parse_args()

    subhaloID = args.subhaloID
    subhaloNum = int(subhaloID[7:])
    snapNum = int(args.snapNum)
    
    rotation = args.rotation
    phi = float(args.phi)
    inc = float(args.inc)
    pa = float(args.pa)

    illustris_path = paths.illustris_path
    if args.TNG:
        illustris_path = paths.TNG_path

    outpath = args.outpath
    if outpath is None:
        if args.TNG:
            outpath = paths.TNG_savepath+'/data/snap%03d/' % snapNum
        else:
            outpath = paths.illustris_savepath+'/data/snap%03d/' % snapNum

    foldername = subhaloID + '_p' + args.phi + '_i' + args.inc
    os.system('mkdir -p {}/{}'.format(outpath, foldername))
    path = outpath + foldername
    os.system('mkdir -p {}/imgs'.format(path))
    os.system('mkdir -p {}/mge'.format(path))
    os.system('mkdir -p {}/ifu'.format(path))

    data = ui.getData(illustris_path, snapNum, subhaloNum, 4)

    Illustris.mock.makeImg.run(data['Coordinates'], data['Velocities'], data['Masses'], 
        data['GFM_StellarPhotometrics'], data['GFM_Metallicity'], 
        rotation, phi, inc, pa, path+'/imgs/')
    Illustris.mock.makeMGE.run(path+'/imgs/img_M.npy', path+'/mge/', args.iteration)
    Illustris.mock.makeIFUbins.run(path+'/imgs/img_ifu.npy', path+'/imgs/img_M.npy', 
        path+'/mge/mge.npy', path+'/ifu/', circle = args.circle)
    Illustris.mock.makeIFUdata.run(path+'/ifu/', path+'/imgs/')
    Illustris.mock.dumpIFUdata.run(path+'/ifu/')
    Illustris.mock.v_plot.run(path+'/ifu/')

if __name__ == '__main__':
    main()

