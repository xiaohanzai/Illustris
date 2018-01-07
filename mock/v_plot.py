import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.io import fits
import argparse
from JAM.utils.velocity_plot import velocity_plot
from scipy import stats
import Illustris.utils.paths as paths
outpath = paths.illustris_savepath

def v_plot(x0, y0, v0, v0_err, vd, vd_err):
    r = np.sqrt(x0**2+y0**2)
    ii = np.where(r < 3.0)
    mv0 = np.mean(v0[ii])
    vel = v0-mv0
    good = (vel**2+vd**2)**0.5 < 800.
    x0 = x0[good]
    y0 = y0[good]
    vel = vel[good]
    v0_err = v0_err[good]
    vd = vd[good]
    vd_err = vd_err[good]

    fig = plt.figure(figsize=(4*1.5, 3.3*1.5))
    fig.subplots_adjust(left=0.07, bottom=0.05, right=0.88,
                        top=0.98, wspace=0.6, hspace=0.01)

    ax1 = fig.add_subplot(2, 2, 1)
    nans = np.isnan(vel)
    vmax = stats.scoreatpercentile(vel[~nans], 98.0)
    norm = colors.Normalize(vmin=-vmax, vmax=vmax)
    velocity_plot(x0, y0, vel, markersize=0.2, norm=norm,
                  ax=ax1, text='$\mathbf{V^{*}}$', equal=True,
                  xreverse=False)
    ax2 = fig.add_subplot(2, 2, 2)
    velocity_plot(x0, y0, v0_err.clip(0., 100.), markersize=0.2,
                  ax=ax2, text='$\mathbf{V^{*}_{err}}$', equal=True,
                  xreverse=False)
    ax3 = fig.add_subplot(2, 2, 3)
    velocity_plot(x0, y0, vd.clip(0, 600.), markersize=0.2,
                  ax=ax3, text='$\mathbf{\sigma^{*}}$',
                  xreverse=False)
    ax4 = fig.add_subplot(2, 2, 4)
    velocity_plot(x0, y0, vd_err.clip(0, 100.), markersize=0.2,
                  ax=ax4, text='$\mathbf{\sigma^{*}_{err}}$',
                  xreverse=False)
    
    return fig

def run(ifu_path):
    if not os.path.isfile('%s/IFU.fits'%ifu_path):
        print('No IFU.fits in %s'%ifu_path)
        sys.exit(1)
    hdulist = fits.open('%s/IFU.fits'%ifu_path)[1]
    tem = hdulist.data
    x0 = tem['xbin']
    y0 = tem['ybin']
    v0 = tem['v0']
    v0_err = tem['v0_err']
    vd = tem['vd']
    vd_err = tem['vd_err']

    fig = v_plot(x0, y0, v0, v0_err, vd, vd_err)
    fig.savefig('%s/IFU.eps'%ifu_path, dpi=300)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--foldername')
    args = parser.parse_args()
    
    foldername = args.foldername
    path = outpath + foldername
    run(path+'/ifu/')

if __name__ == '__main__':
    main()






# def v_plot(path):
#     if not os.path.isfile('%s/ifu/IFU.fits'%path):
#         print('No IFU.fits in %s'%path)
#         sys.exit(1)
#     hdulist = fits.open('%s/ifu/IFU.fits'%path)[1]
#     tem = hdulist.data
#     x0 = tem['xbin']
#     y0 = tem['ybin']
#     v0 = tem['v0']
#     v0_err = tem['v0_err']
#     vd = tem['vd']
#     vd_err = tem['vd_err']
#     r = np.sqrt(x0**2+y0**2)
#     ii = np.where(r < 3.0)
#     mv0 = np.mean(v0[ii])
#     vel = v0-mv0
#     good = (vel**2+vd**2)**0.5 < 800.
#     x0 = x0[good]
#     y0 = y0[good]
#     vel = vel[good]
#     v0_err = v0_err[good]
#     vd = vd[good]
#     vd_err = vd_err[good]

#     fig = plt.figure(figsize=(4*1.5, 3.3*1.5))
#     fig.subplots_adjust(left=0.07, bottom=0.05, right=0.88,
#                         top=0.98, wspace=0.6, hspace=0.01)

#     ax1 = fig.add_subplot(2, 2, 1)
#     nans = np.isnan(vel)
#     vmax = stats.scoreatpercentile(vel[~nans], 98.0)
#     norm = colors.Normalize(vmin=-vmax, vmax=vmax)
#     velocity_plot(x0, y0, vel, markersize=0.2, norm=norm,
#                   ax=ax1, text='$\mathbf{V^{*}}$', equal=True,
#                   xreverse=False)
#     ax2 = fig.add_subplot(2, 2, 2)
#     velocity_plot(x0, y0, v0_err.clip(0., 100.), markersize=0.2,
#                   ax=ax2, text='$\mathbf{V^{*}_{err}}$', equal=True,
#                   xreverse=False)
#     ax3 = fig.add_subplot(2, 2, 3)
#     velocity_plot(x0, y0, vd.clip(0, 600.), markersize=0.2,
#                   ax=ax3, text='$\mathbf{\sigma^{*}}$',
#                   xreverse=False)
#     ax4 = fig.add_subplot(2, 2, 4)
#     velocity_plot(x0, y0, vd_err.clip(0, 100.), markersize=0.2,
#                   ax=ax4, text='$\mathbf{\sigma^{*}_{err}}$',
#                   xreverse=False)
#     fig.savefig('%s/ifu/IFU.eps'%path, dpi=300)

# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--subhaloID')
#     args = parser.parse_args()
    
#     subhaloID = args.subhaloID
#     path = outpath + subhaloID
#     v_plot(path)

# if __name__ == '__main__':
#     main()





