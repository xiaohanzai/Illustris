'''
Modified from Hongyu Li's script: https://github.com/HongyuLi2016/illustris-tools

read img_ifu.npy (particel number map), mge.npy (mge fitting)
and img_M.npy (mass image) to create voronoi bin for IFU
'''
import sys
sys.path.append('../utils/')
import numpy as np
from scipy.spatial import ConvexHull
import argparse
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from voronoi_2d_binning import voronoi_2d_binning
import util_illustris as ui
import paths
outpath = paths.illustris_savepath

pix2kpc = ui.scale_ifu   # 1 pixel = 1 kpc
kpc2arcsec = ui.kpc2arcsec
pix2arcsec = pix2kpc * kpc2arcsec  # pixel / arcsec

def getRe(sol):
    '''
    Effective radius.
    '''
    L = 2*np.pi*sol[:, 0]*sol[:, 1]**2*sol[:, 2]
    R = np.logspace(np.log10(0.5), np.log10(50), 5000)
    Lu = R.copy()
    for i in range(R.size):
        Lu[i] = np.sum(L*(1-np.exp(-(R[i])**2/(2.*sol[:, 1]**2*sol[:, 2]))))
    tLu = np.sum(L)/2.
    ii = np.where(np.abs(Lu-tLu) == np.min(np.abs(Lu-tLu)))
    return R[ii]

def makeIFUbins(img_ifu, img_M, sol, pa, eps, path, circle = True):
    '''
    sol, pa, eps are MGE fitting results.
    Save relavent files to "path".
    '''
    ba = 1.0-eps
    theta = -((90.-pa)/180.*np.pi)
    Re = getRe(sol)
    print('Re: %.2f'%Re)
    aa = Re/np.sqrt(ba)
    bb = Re*np.sqrt(ba)

    naxis_ifu = img_ifu.shape[0]
    naxis_img = img_M.shape[0]
    # !!! be careful, x,y should be in arcsec, not kpc!
    x = np.linspace(-naxis_ifu*pix2arcsec/2,
                    naxis_ifu*pix2arcsec/2, naxis_ifu)

    y = np.linspace(-naxis_ifu*pix2arcsec/2,
                    naxis_ifu*pix2arcsec/2, naxis_ifu)
    X, Y = np.meshgrid(x, y)

    Xbin = np.cos(theta)*X+np.sin(theta)*Y
    Ybin = -np.sin(theta)*X+np.cos(theta)*Y
    if circle:
        aa = Re
        bb = Re
    R = ((Xbin/aa)**2+(Ybin/bb)**2)**0.5
    ii = (R < 2.5)  # *(~(((Xbin+73.0)**2 + (Ybin+15.0)**2)**0.5<15.0))
    xbin = X[ii].reshape(-1)  # x position for each spaxel
    ybin = Y[ii].reshape(-1)  # y position for each spaxel
    signal = img_ifu[ii].reshape(-1)
    signal = signal.clip(1.0)
    print((signal == 0).sum())
    noise = np.sqrt(signal)
    targetSN = 20.0
    binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale, fig = \
        voronoi_2d_binning(xbin, ybin, signal, noise, targetSN, plot=1, quiet=1)
    fig.savefig('%s/voronoi_bin.eps'%path, dpi=300)

    # Re and 2.5 Re ellipse
    phi = np.linspace(0, 2*np.pi, 100)
    xx_tem = aa * np.cos(phi) / kpc2arcsec
    yy_tem = bb * np.sin(phi) / kpc2arcsec
    xx = np.cos(-theta)*xx_tem+np.sin(-theta)*yy_tem  # xx, yy in kpc
    yy = -np.sin(-theta)*xx_tem+np.cos(-theta)*yy_tem

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.axis('equal')
    ax.set_adjustable('box-forced')
    ax.imshow(np.log10(img_M), origin='lower',
              extent=(-ui.boxsize_img/2, ui.boxsize_img/2,
                      -ui.boxsize_img/2, ui.boxsize_img/2))
    ax.plot(xNode/kpc2arcsec, yNode/kpc2arcsec, 'ok',
            markersize=1.0, alpha=0.3)
    ax.plot(xx, yy, 'c', alpha=0.6, lw=2)
    ax.plot(2.5*xx, 2.5*yy, 'c', alpha=0.6, lw=2)
    ax.set_xlim([-ui.boxsize_img/2, ui.boxsize_img/2])
    ax.set_ylim([-ui.boxsize_img/2, ui.boxsize_img/2])
    plt.xlabel("kpc")
    plt.ylabel("kpc")
    fig.savefig('%s/img.eps'%path, dpi=300)

    # output some necessary files and figures
    # convex hull
    hull = ConvexHull(np.array([xbin, ybin]).T)
    x_hull = xbin[hull.vertices][::-1]
    y_hull = ybin[hull.vertices][::-1]  # must be clockwise!
    R = (x_hull**2 + y_hull**2)**0.5
    Rmax = np.mean(R.max())
    r = R.copy() * 0.0
    for i in range(len(r)):
        position = np.array([[x_hull[i-1], x_hull[i]],
                             [y_hull[i-1], y_hull[i]]])
        vector1 = np.array([position[0, 1], position[1, 1]])
        vector2 = np.array([position[0, 1]-position[0, 0],
                            position[1, 1]-position[1, 0]])
        project_length = abs(np.dot(vector1, vector2) /
                             np.sqrt(np.dot(vector2, vector2)))
        r[i] = (np.dot(vector1, vector1) - project_length**2)**0.5
    Rmin = np.mean(r.min())
    Rect_x_min = np.mean(x_hull.min())
    Rect_x_max = np.mean(x_hull.max())
    Rect_y_min = np.mean(y_hull.min())
    Rect_y_max = np.mean(y_hull.max())
    with open('%s/IFU_hull'%path, 'w') as ff:
        print(('{0:d}  {1:+e}  {2:+e}  {3:+e}  {4:+e}  {5:+e}  {6:+e}'
            .format(len(x_hull)+1, Rmin, Rmax, Rect_x_min,
                    Rect_x_max, Rect_y_min, Rect_y_max)), file = ff)
        for i in range(len(x_hull)):
            print('{0:+e}  {1:+e}'.format(x_hull[i], y_hull[i]), file = ff)
        print('{0:+e}  {1:+e}'.format(x_hull[0], y_hull[0]), file = ff)
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(xNode, yNode, 'o', markersize=2.0)
    ax.plot(xbin, ybin, '.r', markersize=0.5, alpha=0.5)
    for simplex in hull.simplices:
        ax.plot(hull.points[simplex, 0], hull.points[simplex, 1], 'k-')
    circle = Circle(xy=(0.0, 0.0), fc='none', radius=Rmin,
                    ec='c', zorder=1, lw=2.)
    ax.add_artist(circle)
    circle = Circle(xy=(0.0, 0.0), fc='none', radius=Rmax,
                    ec='c', zorder=1, lw=2.)
    ax.add_artist(circle)
    squre = Rectangle(xy=(Rect_x_min, Rect_y_min), fc='none',
                      width=Rect_x_max-Rect_x_min,
                      height=Rect_y_max-Rect_y_min, ec='r', zorder=1, lw=2.)
    ax.add_artist(squre)
    ax.set_aspect(1)
    ax.set_aspect('equal', adjustable='box', anchor='C')
    xlim = [Rect_x_min, Rect_x_max]
    ylim = [Rect_y_min, Rect_y_max]
    lim = [min(xlim[0], ylim[0]), max(xlim[1], ylim[1])]
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('x arcsec')
    ax.set_ylabel('y arcsec')
    fig.savefig('%s/IFU_hull.eps'%path, dpi=300)
    # spaxel bins
    with open('%s/spaxel_bins.dat'%path, 'w') as ff:
        for i in range(len(xbin)):
            print('{:12f} {:12f} {:6d}'.format(xbin[i],ybin[i],binNum[i]), file = ff)
    # voronoi bins
    with open('%s/voronoi_bins.dat'%path, 'w') as ff:
        for i in range(len(xNode)):
            print(('{:6d} {:12f} {:12f} {:12f} {:6d}'
                .format(i, xNode[i], yNode[i], 0.1, 1)), file = ff)

def run(img_ifu_file, img_M_file, mge_file, path, circle = True):
    try:
        img_ifu = np.load(img_ifu_file, encoding = 'bytes')
        img_M = np.load(img_M_file, encoding = 'bytes')
    except:
        print('Error - no image file found')
        sys.exit(1)

    try:
        sol, pa, eps, xmed, ymed = np.load(mge_file, encoding = 'bytes')
    except:
        print('Error - no MGE file found')
        sys.exit(1)

    makeIFUbins(img_ifu, img_M, sol, pa, eps, path, circle = circle)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--subhaloID')
    parser.add_argument('-c', action='store_true', dest='circle',
                      default=False, help='use circlular aperture')
    args = parser.parse_args()

    subhaloID = args.subhaloID
    subhaloNum = int(subhaloID[7:])
    path = outpath + subhaloID
    os.system('mkdir -p %s/ifu/'%path)
    
    run(path+'/imgs/img_ifu.npy', path+'/imgs/img_M.npy', path+'/mge/mge.npy',
        path+'/ifu/', circle = args.circle)

if __name__ == '__main__':
    main()






# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--subhaloID')
#     parser.add_argument('-c', action='store_true', dest='circle',
#                       default=False, help='use circlular aperture')
#     args = parser.parse_args()

#     subhaloID = args.subhaloID
#     subhaloNum = int(subhaloID[7:])
#     path = outpath + subhaloID
#     os.system('mkdir -p %s/ifu/'%path)

#     try:
#         img_ifu = np.load('%s/imgs/img_ifu.npy'%path)
#         img_M = np.load('%s/imgs/img_M.npy'%path)
#     except:
#         print('Error - no image file found')
#         sys.exit(1)

#     try:
#         sol, pa, eps, xmed, ymed = np.load('%s/mge/mge.npy'%path)
#         ba = 1.0-eps
#         theta = -((90.-pa)/180.*np.pi)
#         Re = getRe(sol)
#         print('Re: %.2f'%Re)
#         aa = Re/np.sqrt(ba)
#         bb = Re*np.sqrt(ba)
#     except:
#         print('Error - no MGE file found')
#         sys.exit(1)
    
#     naxis_ifu = img_ifu.shape[0]
#     naxis_img = img_M.shape[0]
#     # !!! be careful, x,y should be in arcsec, not kpc!
#     x = np.linspace(-naxis_ifu*pix2arcsec/2,
#                     naxis_ifu*pix2arcsec/2, naxis_ifu)

#     y = np.linspace(-naxis_ifu*pix2arcsec/2,
#                     naxis_ifu*pix2arcsec/2, naxis_ifu)
#     X, Y = np.meshgrid(x, y)

#     Xbin = np.cos(theta)*X+np.sin(theta)*Y
#     Ybin = -np.sin(theta)*X+np.cos(theta)*Y
#     if args.circle:
#         aa = Re
#         bb = Re
#     R = ((Xbin/aa)**2+(Ybin/bb)**2)**0.5
#     ii = (R < 2.5)  # *(~(((Xbin+73.0)**2 + (Ybin+15.0)**2)**0.5<15.0))
#     xbin = X[ii].reshape(-1)  # x position for each spaxel
#     ybin = Y[ii].reshape(-1)  # y position for each spaxel
#     signal = img_ifu[ii].reshape(-1)
#     signal = signal.clip(1.0)
#     print((signal == 0).sum())
#     noise = np.sqrt(signal)
#     targetSN = 20.0
#     binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale, fig = \
#         voronoi_2d_binning(xbin, ybin, signal, noise, targetSN, plot=1, quiet=1)
#     fig.savefig('%s/ifu/voronoi_bin.eps'%path, dpi=300)

#     # Re and 2.5 Re ellipse
#     phi = np.linspace(0, 2*np.pi, 100)
#     xx_tem = aa * np.cos(phi) / kpc2arcsec
#     yy_tem = bb * np.sin(phi) / kpc2arcsec
#     xx = np.cos(-theta)*xx_tem+np.sin(-theta)*yy_tem  # xx, yy in kpc
#     yy = -np.sin(-theta)*xx_tem+np.cos(-theta)*yy_tem

#     fig = plt.figure(figsize=(6, 6))
#     ax = fig.add_subplot(1, 1, 1)
#     ax.axis('equal')
#     ax.set_adjustable('box-forced')
#     ax.imshow(np.log10(img_M), origin='lower',
#               extent=(-ui.boxsize_img/2, ui.boxsize_img/2,
#                       -ui.boxsize_img/2, ui.boxsize_img/2))
#     ax.plot(xNode/kpc2arcsec, yNode/kpc2arcsec, 'ok',
#             markersize=1.0, alpha=0.3)
#     ax.plot(xx, yy, 'c', alpha=0.6, lw=2)
#     ax.plot(2.5*xx, 2.5*yy, 'c', alpha=0.6, lw=2)
#     ax.set_xlim([-ui.boxsize_img/2, ui.boxsize_img/2])
#     ax.set_ylim([-ui.boxsize_img/2, ui.boxsize_img/2])
#     plt.xlabel("kpc")
#     plt.ylabel("kpc")
#     fig.savefig('%s/ifu/img.eps'%path, dpi=300)

#     # output some necessary files and figures
#     # convex hull
#     hull = ConvexHull(np.array([xbin, ybin]).T)
#     x_hull = xbin[hull.vertices][::-1]
#     y_hull = ybin[hull.vertices][::-1]  # must be clockwise!
#     R = (x_hull**2 + y_hull**2)**0.5
#     Rmax = np.mean(R.max())
#     r = R.copy() * 0.0
#     for i in range(len(r)):
#         position = np.array([[x_hull[i-1], x_hull[i]],
#                              [y_hull[i-1], y_hull[i]]])
#         vector1 = np.array([position[0, 1], position[1, 1]])
#         vector2 = np.array([position[0, 1]-position[0, 0],
#                             position[1, 1]-position[1, 0]])
#         project_length = abs(np.dot(vector1, vector2) /
#                              np.sqrt(np.dot(vector2, vector2)))
#         r[i] = (np.dot(vector1, vector1) - project_length**2)**0.5
#     Rmin = np.mean(r.min())
#     Rect_x_min = np.mean(x_hull.min())
#     Rect_x_max = np.mean(x_hull.max())
#     Rect_y_min = np.mean(y_hull.min())
#     Rect_y_max = np.mean(y_hull.max())
#     with open('%s/ifu/IFU_hull'%path, 'w') as ff:
#         print(('{0:d}  {1:+e}  {2:+e}  {3:+e}  {4:+e}  {5:+e}  {6:+e}'
#             .format(len(x_hull)+1, Rmin, Rmax, Rect_x_min,
#                     Rect_x_max, Rect_y_min, Rect_y_max)), file = ff)
#         for i in range(len(x_hull)):
#             print('{0:+e}  {1:+e}'.format(x_hull[i], y_hull[i]), file = ff)
#         print('{0:+e}  {1:+e}'.format(x_hull[0], y_hull[0]), file = ff)
#     fig = plt.figure(figsize=(6, 6))
#     ax = fig.add_subplot(1, 1, 1)
#     ax.plot(xNode, yNode, 'o', markersize=2.0)
#     ax.plot(xbin, ybin, '.r', markersize=0.5, alpha=0.5)
#     for simplex in hull.simplices:
#         ax.plot(hull.points[simplex, 0], hull.points[simplex, 1], 'k-')
#     circle = Circle(xy=(0.0, 0.0), fc='none', radius=Rmin,
#                     ec='c', zorder=1, lw=2.)
#     ax.add_artist(circle)
#     circle = Circle(xy=(0.0, 0.0), fc='none', radius=Rmax,
#                     ec='c', zorder=1, lw=2.)
#     ax.add_artist(circle)
#     squre = Rectangle(xy=(Rect_x_min, Rect_y_min), fc='none',
#                       width=Rect_x_max-Rect_x_min,
#                       height=Rect_y_max-Rect_y_min, ec='r', zorder=1, lw=2.)
#     ax.add_artist(squre)
#     ax.set_aspect(1)
#     ax.set_aspect('equal', adjustable='box', anchor='C')
#     xlim = [Rect_x_min, Rect_x_max]
#     ylim = [Rect_y_min, Rect_y_max]
#     lim = [min(xlim[0], ylim[0]), max(xlim[1], ylim[1])]
#     ax.set_xlim(lim)
#     ax.set_ylim(lim)
#     ax.set_xlabel('x arcsec')
#     ax.set_ylabel('y arcsec')
#     fig.savefig('%s/ifu/IFU_hull.eps'%path, dpi=300)
#     # spaxel bins
#     with open('%s/ifu/spaxel_bins.dat'%path, 'w') as ff:
#         for i in range(len(xbin)):
#             print('{:12f} {:12f} {:6d}'.format(xbin[i],ybin[i],binNum[i]), file = ff)
#     # voronoi bins
#     with open('%s/ifu/voronoi_bins.dat'%path, 'w') as ff:
#         for i in range(len(xNode)):
#             print(('{:6d} {:12f} {:12f} {:12f} {:6d}'
#                 .format(i, xNode[i], yNode[i], 0.1, 1)), file = ff)

# if __name__ == '__main__':
#     main()


