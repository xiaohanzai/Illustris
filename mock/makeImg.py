'''
Modified from Hongyu Li's script: https://github.com/HongyuLi2016/illustris-tools
'''
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.switch_backend('agg')
from Illustris.data.StarData import StarData, Rb_all
import Illustris.utils.util_illustris as ui
from Illustris.utils.util_general import rotateCoordinates
import warnings
import Illustris.utils.paths as paths
illustris_path = paths.illustris_path
outpath = paths.illustris_savepath

warnings.simplefilter("ignore")

boxsize_img = ui.boxsize_img
scale_img = ui.scale_img
boxsize_ifu = ui.boxsize_ifu
scale_ifu = ui.scale_ifu

def _makeImg(x, L, scale=0.5, boxsize=50.0):
    '''
    Make a mock image, pixle size = 0.5 kpc/pix. Image size = 50.0 kpc * 50.0 kpc.
    x is the coordinate vector of the particles. L = light, or mass, or whatever.
    '''
    ngrid = int(boxsize/scale)
    img = np.zeros([ngrid, ngrid])
    for n in range(len(x)):
        # assume z is LOS
        i = int(np.rint(x[n, 0]/scale) + int(ngrid/2))
        k = int(np.rint(x[n, 1]/scale) + int(ngrid/2))
        if i >= 0 and i < ngrid and k >= 0 and k < ngrid:
            img[k, i] += L[n]
    return img

def _plotImg(img, boxsize):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.imshow(np.log10(img), origin='lower',
              extent=(-boxsize/2., boxsize/2., -boxsize/2., boxsize/2.))
    ax.set_xlabel('kpc')
    ax.set_ylabel('kpc')
    return fig

def makeImg(x, v, mpart, lrpart, rotation, phi, inc, pa = 500, Rb = Rb_all, **kwargs):
    '''
    rotation = oblate or prolate.
    phi inc are be in units of degrees.
    You can also provide pa in kwargs.
    '''
    data = StarData(x, v, mpart, Rb = Rb, **kwargs)
        
    if rotation not in ['oblate', 'prolate']:
        print('rotation must be oblate or prolate!')
        return

    # convert to principal axes coordinates
    data.convert2cyl(rotation)

    # rotate x v
    xpart = rotateCoordinates(data.x_prin, np.radians(phi), np.radians(inc), 0)
    if pa > 360.:
        pa = ui.calcPa(xpart, Rb = Rb)[1]
    xpart = rotateCoordinates(xpart, 0, 0, np.radians(pa))
    vpart = rotateCoordinates(data.v_prin, np.radians(phi), np.radians(inc), np.radians(pa))

    # make images
    img_M = _makeImg(xpart, data.mpart, scale=scale_img, boxsize=boxsize_img)
    img_L = _makeImg(xpart, lrpart, scale=scale_img, boxsize=boxsize_img)
    img_ifu = _makeImg(xpart, np.ones_like(lrpart), scale=scale_ifu, boxsize=boxsize_ifu)

    fig_M = _plotImg(img_M, boxsize_img)
    fig_L = _plotImg(img_L, boxsize_img)
    fig_ifu = _plotImg(img_ifu, boxsize_ifu)

    return img_M, img_L, img_ifu, fig_M, fig_L, fig_ifu, xpart, vpart, data.ba, data.ca, pa

def run(x, v, mpart, lpart, Z, rotation, phi, inc, pa, path, Rb = Rb_all, **kwargs):
    '''
    lpart = 8 band magnitude.
    Z = metallicity.
    phi inc are in units of degrees.
    Save all relavent files to "path".
    '''
    lrpart = 10**((4.58 - lpart[:, 5])/2.5)/1e10  # r band luminosity
    mpart /= 1e10
    img_M, img_L, img_ifu, fig_M, fig_L, fig_ifu, xpart, vpart, ba, ca, pa = \
        makeImg(x, v, mpart, lrpart, rotation, phi, inc, pa, Rb = Rb, **kwargs)

    np.save('{}/img_M.npy'.format(path), img_M)
    np.save('{}/img_L.npy'.format(path), img_L)
    np.save('{}/img_ifu.npy'.format(path), img_ifu)

    fig_M.savefig('{}/img_M.eps'.format(path), dpi=300)
    fig_L.savefig('{}/img_L.eps'.format(path), dpi=300)
    fig_ifu.savefig('{}/img_ifu.eps'.format(path), dpi=300)

    with open('{}/makeImg.log'.format(path), 'w') as ff:
        print('inc: {:.2f}'.format(inc), file = ff)
        print('phi: {:.2f}'.format(phi), file = ff)
        print('pa: {:.2f}'.format(pa), file = ff)
        print('boxsize_img: {:.2f} kpc'.format(boxsize_img), file = ff)
        print('boxsize_ifu: {:.2f} kpc'.format(boxsize_ifu), file = ff)
        print('scale_img: {:.2f} kpc/pixel'.format(scale_img), file = ff)
        print('scale_ifu: {:.2f} kpc/pixel'.format(scale_ifu), file = ff)
        print('Stellar mass: {:.4e} M_solar'.format(mpart.sum()*1e10), file = ff)
        print(('r band luminosity: {:.4e} L_solar'
                      .format(lrpart.sum()*1e10)), file = ff)
        print(('Averaged M*/L: {:.3f} M_solar/L_solar'
                     .format(mpart.sum()/lrpart.sum())), file = ff)
        print('ba={:.2f}  ca={:.2f}'.format(ba, ca), file = ff)

    # x in kpc, v in km/s, mpart in 10^10 M_solar, lrpart in 10^10 L_solar
    data = np.zeros([xpart.shape[0], 9])
    data[:, 0:3] = xpart
    data[:, 3:6] = vpart
    data[:, 6] = mpart
    data[:, 7] = lrpart
    data[:, 8] = Z
    np.save('{}/coordinates_star.npy'.format(path), data)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subhaloID')
    parser.add_argument('--snapNum', default = '135')
    parser.add_argument('--rotation') # oblate or prolate rotation
    parser.add_argument('--phi') # let's input in units of degrees
    parser.add_argument('--inc')
    parser.add_argument('--pa', default = '500')
    args = parser.parse_args()

    subhaloID = args.subhaloID
    subhaloNum = int(subhaloID[7:])
    snapNum = int(args.snapNum)
    
    rotation = args.rotation
    phi = float(args.phi)
    inc = float(args.inc)
    pa = float(args.pa)

    os.system('mkdir -p {}/{}'.format(outpath, subhaloID))
    path = outpath + subhaloID
    os.system('mkdir -p {}/imgs'.format(path))

    data = ui.getData(illustris_path, snapNum, subhaloNum, 4)
    run(data['Coordinates'], data['Velocities'], data['Masses'], 
        data['GFM_StellarPhotometrics'], data['GFM_Metallicity'], 
        rotation, phi, inc, pa, path+'/imgs/')

if __name__ == '__main__':
    main()





# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--subhaloID')
#     parser.add_argument('--rotation') # oblate or prolate rotation
#     parser.add_argument('--phi')
#     parser.add_argument('--inc')
#     parser.add_argument('--pa')
#     args = parser.parse_args()

#     subhaloID = args.subhaloID
#     subhaloNum = int(subhaloID[7:])
#     rotation = args.rotation
#     os.system('mkdir -p {}/{}'.format(outpath, subhaloID))
#     path = outpath + subhaloID
#     os.system('mkdir -p {}/imgs'.format(path))

#     phi = float(args.phi)
#     inc = float(args.inc)
#     pa = float(args.pa)

#     tmp = ui.getData(illustris_path, 135, subhaloNum, 4)
#     xpart_star = tmp['Coordinates']
#     vpart_star = tmp['Velocities']
#     mpart_star = tmp['Masses']
#     Meta = tmp['GFM_Metallicity']
#     lpart = tmp['GFM_StellarPhotometrics']
#     lrpart = 10**((4.58 - lpart[:, 5])/2.5)/1e10  # r band luminosity
#     xcenter_star = ui.findCenter(xpart_star, mpart=mpart_star, percent=20.0)
#     xpart_star -= xcenter_star
#     vsys_star = ui.getVc(xpart_star, vpart_star, mpart=mpart_star, R=10.0)
#     vpart_star -= vsys_star

#     tmp = ui.getData(illustris_path, 135, subhaloNum, 1)
#     xpart_dark = tmp['Coordinates']
#     vpart_dark = tmp['Velocities']
#     xcenter_dark = ui.findCenter(xpart_dark, percent=20.0)
#     xpart_dark -= xcenter_dark
#     vsys_dark = ui.getVc(xpart_dark, vpart_dark, R=10.0)
#     vpart_dark -= vsys_dark

#     # difference between halo and galaxy
#     centerDiff = xcenter_dark - xcenter_star
#     vsysDiff = vsys_dark - vsys_star

#     # get galaxy principal axes
#     ba, ca, a1, a2, Tiv = ui.calcShape(xpart_star, Rb=10.)
#     if rotation == 'oblate':
#         Tiv = Tiv[[2,1,0],:]

#     xpart_axis = np.dot(Tiv, xpart_star.T).T
#     vpart_axis = np.dot(Tiv, vpart_star.T).T

#     xpart = rotateCoordinates(xpart_axis, np.radians(phi), np.radians(inc), np.radians(pa)) # rotate x
#     vpart = rotateCoordinates(vpart_axis, np.radians(phi), np.radians(inc), np.radians(pa)) # rotate v

#     img_M = _makeImg(xpart, mpart_star, scale=scale_img, boxsize=boxsize_img)
#     img_L = _makeImg(xpart, lrpart, scale=scale_img, boxsize=boxsize_img)
#     img_ifu = _makeImg(xpart, lrpart*0.0+1, scale=scale_ifu,
#                        boxsize=boxsize_ifu)

#     # save necessary files
#     os.system('mkdir -p {}/imgs'.format(path))

#     np.save('{}/imgs/img_M.npy'.format(path), img_M)
#     np.save('{}/imgs/img_L.npy'.format(path), img_L)
#     np.save('{}/imgs/img_ifu.npy'.format(path), img_ifu)

#     fig = plt.figure(figsize=(6, 6))
#     ax = fig.add_subplot(1, 1, 1)
#     ax.imshow(np.log10(img_M), origin='lower',
#               extent=(-boxsize_img/2, boxsize_img/2, -boxsize_img/2,
#                       boxsize_img/2))
#     ax.set_xlabel('kpc')
#     ax.set_ylabel('kpc')
#     fig.savefig('{}/imgs/img_M.eps'.format(path), dpi=300)

#     fig = plt.figure(figsize=(6, 6))
#     ax = fig.add_subplot(1, 1, 1)
#     ax.imshow(np.log10(img_L), origin='lower',
#               extent=(-boxsize_img/2, boxsize_img/2, -boxsize_img/2,
#                       boxsize_img/2))
#     ax.set_xlabel('kpc')
#     ax.set_ylabel('kpc')
#     fig.savefig('{}/imgs/img_L.eps'.format(path), dpi=300)

#     fig = plt.figure(figsize=(6, 6))
#     ax = fig.add_subplot(1, 1, 1)
#     ax.imshow(np.log10(img_ifu), origin='lower',
#               extent=(-boxsize_ifu/2, boxsize_ifu/2, -boxsize_ifu/2,
#                       boxsize_ifu/2))
#     ax.set_xlabel('kpc')
#     ax.set_ylabel('kpc')
#     fig.savefig('{}/imgs/img_ifu.eps'.format(path), dpi=300)

#     with open('{}/imgs/makeImg.log'.format(path), 'w') as ff:
#         print('inclination: {:.2f}'.format(inc), file = ff)
#         print('phi: {:.2f}'.format(phi), file = ff)
#         print('pa: {:.2f}'.format(pa), file = ff)
#         print('boxsize_img: {:.2f} kpc'.format(boxsize_img), file = ff)
#         print('boxsize_ifu: {:.2f} kpc'.format(boxsize_ifu), file = ff)
#         print('scale_img: {:.2f} kpc/pixel'.format(scale_img), file = ff)
#         print('scale_ifu: {:.2f} kpc/pixel'.format(scale_ifu), file = ff)
#         print('Stellar mass: {:.4e} M_solar'.format(mpart_star.sum()*1e10), file = ff)
#         print(('r band luminosity: {:.4e} L_solar'
#                       .format(lrpart.sum()*1e10)), file = ff)
#         print(('Averaged M*/L: {:.3f} M_solar/L_solar'
#                      .format(mpart_star.sum()/lrpart.sum())), file = ff)
#         print(('Mass centre offset: {:.2f} {:.2f} {:.2f} kpc'
#                      .format(centerDiff[0], centerDiff[1], centerDiff[2])), file = ff)
#         print(('System velocity offset: {:.2f} {:.2f} {:.2f} km/s'
#                      .format(vsysDiff[0], vsysDiff[1], vsysDiff[2])), file = ff)
#         print('ba={:.2f}  ca={:.2f}'.format(ba, ca), file = ff)
#     # x in kpc, v in km/s, mpart_star in 10^10 M_solar, lrpart in 10^10 L_solar
#     data = np.zeros([xpart.shape[0], 9])
#     data[:, 0:3] = xpart
#     data[:, 3:6] = vpart
#     data[:, 6] = mpart_star
#     data[:, 7] = lrpart
#     data[:, 8] = Meta
#     np.save('{}/imgs/coordinates_star.npy'.format(path), data)

# if __name__ == '__main__':
#     main()









