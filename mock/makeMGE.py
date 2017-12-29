'''
Modified from Hongyu Li's script: https://github.com/HongyuLi2016/illustris-tools
'''
import sys
sys.path.append('../utils/')
import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import argparse
from sectors_photometry import sectors_photometry
from mge_print_contours import mge_print_contours
from find_galaxy import find_galaxy
from mge_fit_sectors import mge_fit_sectors
import util_illustris as ui
import paths
outpath = paths.illustris_savepath

pix2kpc = ui.scale_img   # 1 pixel = 0.5 kpc
kpc2arcsec = ui.kpc2arcsec   # 1 kpc = 1.612 arcsec
pix2arcsec = pix2kpc * kpc2arcsec   # pixel / arcsec

def makeMGE(img, path, iteration = True):
    '''
    Save the relavant files to "path".
    '''
    part_numb = img/1.5e6*1e10
    ii = part_numb > 10
    
    level = img[ii].min()
    level *= 0.1
    old_img = img.copy()
    # img[55:85,115:145]=0
    lhy_f = find_galaxy(img, plot=1, fraction=0.05,  # level=level,
                        quiet=1, path=path)
    eps = lhy_f.eps
    theta = lhy_f.theta
    pa = np.mod(270 - theta, 180)
    xmed = lhy_f.xmed
    ymed = lhy_f.ymed
    xpeak = lhy_f.xpeak
    ypeak = lhy_f.ypeak
    lhy_s = sectors_photometry(img, eps, theta, xpeak, ypeak,
                               minlevel=level, plot=1, path=path)
    radius = lhy_s.radius
    angle = lhy_s.angle
    counts = lhy_s.counts
    qbound = 0.06
    lhy_mge = mge_fit_sectors(radius, angle, counts, eps, ngauss=15,
                              sigmapsf=0, scale=pix2arcsec,
                              qbounds=[qbound, 0.999], linear=False, quiet=True,
                              outer_slope=4, bulge_disk=False,
                              plot=0)
    sol = lhy_mge.sol.T
    absdev = lhy_mge.absdev

    if iteration:
        sol_old = sol
        absdev_old = absdev
        absdev_init = absdev
        while True:
            qbound += 0.05
            lhy_mge = mge_fit_sectors(radius, angle, counts, eps, ngauss=15,
                                      sigmapsf=0, scale=pix2arcsec,
                                      qbounds=[qbound, 0.999], linear=False,
                                      quiet=True, outer_slope=4,
                                      bulge_disk=False, plot=0)
            sol = lhy_mge.sol.T
            absdev = lhy_mge.absdev

            if (absdev/absdev_old > 1.10 or absdev/absdev_init >
                2.0 or qbound > 0.8):
                break
            absdev_old = absdev
            sol_old = sol
        sol = sol_old

    mge_print_contours(old_img, theta, xpeak, ypeak, sol.T,
                       sigmapsf=0, scale=pix2arcsec, magrange=5,
                       path=path)

    print('total mass: %.4e'%(sol[:, 0].sum()*1e10))
    ml = 5.0  # mock stellar mass-to-light ratio
    sol[:, 1] *= pix2arcsec
    sol[:, 0] = sol[:, 0]*1e10 / \
        (2*np.pi*(sol[:, 1]/kpc2arcsec*1e3)**2*sol[:, 2]) / ml

    return sol, pa, eps, xmed, ymed, absdev

def run(img_file, path, iteration = True):
    try:
        img = np.load(img_file, encoding = 'bytes')
    except:
        print('Error - no image file found')
        sys.exit(1)

    sol, pa, eps, xmed, ymed, absdev = makeMGE(img, path, iteration)

    ff = open('%s/mge.dat'%path, 'w+')
    print('Pa: %5.2f'%pa, file = ff)
    print('Eps: %5.3f'%eps, file = ff)
    print('Absdev: %4.3f'%absdev, file = ff)
    print('Xc Yc: %.3f %.3f'%(xmed, ymed), file = ff)
    for i in range(len(sol[:, 0])):
        print('%10.4e %10.2f %10.3f'%(sol[i, 0], sol[i, 1], sol[i, 2]), file = ff)
    ff.close()
    np.save('%s/mge.npy'%path, [sol, pa, eps, xmed, ymed])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subhaloID')
    parser.add_argument('-i', action='store_false', dest='iteration',
                      default=True, help='Iteration fit')
    args = parser.parse_args()

    subhaloID = args.subhaloID
    subhaloNum = int(subhaloID[7:])

    os.system('mkdir -p {}/{}'.format(outpath, subhaloID))
    path = outpath + subhaloID
    os.system('mkdir -p {}/mge'.format(path))

    run(path+'/imgs/img_M.npy', path+'/mge/', args.iteration)

if __name__ == '__main__':
    main()












# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--subhaloID')
#     parser.add_argument('-i', action='store_false', dest='iteration',
#                       default=True, help='Iteration fit')
#     args = parser.parse_args()

#     subhaloID = args.subhaloID
#     subhaloNum = int(subhaloID[7:])

#     os.system('mkdir -p {}/{}'.format(outpath, subhaloID))
#     path = outpath + subhaloID
#     os.system('mkdir -p {}/mge'.format(path))

#     try:
#         img = np.load('{}/imgs/img_M.npy'.format(path))
#     except:
#         print('Error - no image file found')
#         sys.exit(1)
    
#     part_numb = img/1.5e6*1e10
#     ii = part_numb > 10
#     # plt.imshow(np.log10(img))
#     # img[~ii]=0.
#     # plt.imshow(np.log10(img))
#     # plt.show()

#     level = img[ii].min()
#     level *= 0.1
#     old_img = img.copy()
#     # img[55:85,115:145]=0
#     lhy_f = find_galaxy(img, plot=1, fraction=0.05,  # level=level,
#                         quiet=1, path='%s/mge/'%path)
#     eps = lhy_f.eps
#     theta = lhy_f.theta
#     pa = np.mod(270 - theta, 180)
#     xmed = lhy_f.xmed
#     ymed = lhy_f.ymed
#     xpeak = lhy_f.xpeak
#     ypeak = lhy_f.ypeak
#     lhy_s = sectors_photometry(img, eps, theta, xpeak, ypeak,
#                                minlevel=level, plot=1, path='%s/mge/'%path)
#     radius = lhy_s.radius
#     angle = lhy_s.angle
#     counts = lhy_s.counts
#     qbound = 0.06
#     lhy_mge = mge_fit_sectors(radius, angle, counts, eps, ngauss=15,
#                               sigmapsf=0, scale=pix2arcsec,
#                               qbounds=[qbound, 0.999], linear=False, quiet=True,
#                               outer_slope=4, bulge_disk=False,
#                               plot=0)
#     sol = lhy_mge.sol.T
#     absdev = lhy_mge.absdev

#     if args.iteration:
#         sol_old = sol
#         absdev_old = absdev
#         absdev_init = absdev
#         while True:
#             qbound += 0.05
#             lhy_mge = mge_fit_sectors(radius, angle, counts, eps, ngauss=15,
#                                       sigmapsf=0, scale=pix2arcsec,
#                                       qbounds=[qbound, 0.999], linear=False,
#                                       quiet=True, outer_slope=4,
#                                       bulge_disk=False, plot=0)
#             sol = lhy_mge.sol.T
#             absdev = lhy_mge.absdev

#             if (absdev/absdev_old > 1.10 or absdev/absdev_init >
#                 2.0 or qbound > 0.8):
#                 break
#             absdev_old = absdev
#             sol_old = sol
#         sol = sol_old

#     mge_print_contours(old_img, theta, xpeak, ypeak, sol.T,
#                        sigmapsf=0, scale=pix2arcsec, magrange=5,
#                        path='%s/mge/'%path)

#     print('total mass: %.4e'%(sol[:, 0].sum()*1e10))
#     ml = 5.0  # mock stellar mass-to-light ratio
#     sol[:, 1] *= pix2arcsec
#     sol[:, 0] = sol[:, 0]*1e10 / \
#         (2*np.pi*(sol[:, 1]/kpc2arcsec*1e3)**2*sol[:, 2]) / ml

#     ff = open('%s/mge/mge.dat'%path, 'w+')
#     print('Pa: %5.2f'%pa, file = ff)
#     print('Eps: %5.3f'%eps, file = ff)
#     print('Absdev: %4.3f'%absdev, file = ff)
#     print('Xc Yc: %.3f %.3f'%(xmed, ymed), file = ff)
#     for i in range(len(sol[:, 0])):
#         print('%10.4e %10.2f %10.3f'%(sol[i, 0], sol[i, 1], sol[i, 2]), file = ff)
#     ff.close()
#     np.save('%s/mge/mge.npy'%path, [sol, pa, eps, xmed, ymed])


# if __name__ == '__main__':
#     main()
