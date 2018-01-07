'''
Modified from Hongyu Li's script: https://github.com/HongyuLi2016/illustris-tools
'''
import numpy as np
from astropy.io import fits
import argparse
import Illustris.utils.paths as paths
outpath = paths.illustris_savepath

def run(ifu_path):
    spax_bin_data = np.loadtxt('{}/spaxel_bins.dat'.format(ifu_path))

    spax_bin_X = spax_bin_data[:, 0]
    spax_bin_Y = spax_bin_data[:, 1]
    spax_bin_id = spax_bin_data[:, 2].astype(int)

    bin_data = np.loadtxt('{}/voronoi_bins.dat'.format(ifu_path))

    bin_id = bin_data[:, 0].astype(int)
    bin_X = bin_data[:, 1]
    bin_Y = bin_data[:, 2]
    bin_area = bin_data[:, 3]
    bin_inuse = bin_data[:, 4].astype(int)

    ifu_data = np.loadtxt('{}/IFU_data'.format(ifu_path))
    bin_index = ifu_data[:, 0].astype(int)
    gh_v0 = ifu_data[:, 1]
    gh_v0_err = ifu_data[:, 2]
    gh_vd = ifu_data[:, 3]
    gh_vd_err = ifu_data[:, 4]
    gh_h3 = ifu_data[:, 5]
    gh_h3_err = ifu_data[:, 6]
    gh_h4 = ifu_data[:, 7]
    gh_h4_err = ifu_data[:, 8]
    v0 = ifu_data[:, 9]
    v0_err = ifu_data[:, 10]
    vd = ifu_data[:, 11]
    vd_err = ifu_data[:, 12]
    metal = ifu_data[:, 13]
    flux = ifu_data[:, 14]

    c1 = fits.Column(name='xbin', format='D', array=bin_X)
    c2 = fits.Column(name='ybin', format='D', array=bin_Y)
    c3 = fits.Column(name='v0', format='D', array=gh_v0)
    c4 = fits.Column(name='v0_err', format='D', array=gh_v0_err)
    c5 = fits.Column(name='vd', format='D', array=gh_vd)
    c6 = fits.Column(name='vd_err', format='D', array=gh_vd_err)
    c7 = fits.Column(name='h3', format='D', array=gh_h3)
    c8 = fits.Column(name='h3_err', format='D', array=gh_h3_err)
    c9 = fits.Column(name='h4', format='D', array=gh_h4)
    c10 = fits.Column(name='h4_err', format='D', array=gh_h4_err)
    c11 = fits.Column(name='metal', format='D', array=metal)
    c12 = fits.Column(name='flux', format='D', array=flux)

    c13 = fits.Column(name='rebin_x', format='D', array=spax_bin_X)
    c14 = fits.Column(name='rebin_y', format='D', array=spax_bin_Y)
    c15 = fits.Column(name='binid', format='D', array=spax_bin_id)

    coldefs1 = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12])
    coldefs2 = fits.ColDefs([c13, c14, c15])
    hdu = fits.PrimaryHDU()
    tbhdu1 = fits.BinTableHDU.from_columns(coldefs1)
    tbhdu2 = fits.BinTableHDU.from_columns(coldefs2)
    hdulist = fits.HDUList([hdu, tbhdu1, tbhdu2])
    hdulist.writeto('{}/IFU_gh.fits'.format(ifu_path), clobber=True)

    c1 = fits.Column(name='xbin', format='D', array=bin_X)
    c2 = fits.Column(name='ybin', format='D', array=bin_Y)
    c3 = fits.Column(name='v0', format='D', array=v0)
    c4 = fits.Column(name='v0_err', format='D', array=v0_err)
    c5 = fits.Column(name='vd', format='D', array=vd)
    c6 = fits.Column(name='vd_err', format='D', array=vd_err)
    c7 = fits.Column(name='metal', format='D', array=metal)
    c8 = fits.Column(name='flux', format='D', array=flux)

    c12 = fits.Column(name='rebin_x', format='D', array=spax_bin_X)
    c13 = fits.Column(name='rebin_y', format='D', array=spax_bin_Y)
    c14 = fits.Column(name='binid', format='D', array=spax_bin_id)

    coldefs1 = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8])
    coldefs2 = fits.ColDefs([c12, c13, c14])
    hdu = fits.PrimaryHDU()
    tbhdu1 = fits.BinTableHDU.from_columns(coldefs1)
    tbhdu2 = fits.BinTableHDU.from_columns(coldefs2)
    hdulist = fits.HDUList([hdu, tbhdu1, tbhdu2])
    hdulist.writeto('{}/IFU.fits'.format(ifu_path), clobber=True)

    # plt.plot(bin_X,bin_Y,'.r')
    # plt.plot(spax_bin_X,spax_bin_Y,'.k',alpha=0.1)
    # plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--foldername')
    args = parser.parse_args()

    foldername = args.foldername
    path = outpath + foldername

    run(path+'/ifu/')

if __name__ == '__main__':
    main()



# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--subhaloID')
#     args = parser.parse_args()

#     subhaloID = args.subhaloID
#     subhaloNum = int(subhaloID[7:])
#     path = outpath + subhaloID

#     spax_bin_data = np.loadtxt('{}/ifu/spaxel_bins.dat'.format(path))

#     spax_bin_X = spax_bin_data[:, 0]
#     spax_bin_Y = spax_bin_data[:, 1]
#     spax_bin_id = spax_bin_data[:, 2].astype(int)

#     bin_data = np.loadtxt('{}/ifu/voronoi_bins.dat'.format(path))

#     bin_id = bin_data[:, 0].astype(int)
#     bin_X = bin_data[:, 1]
#     bin_Y = bin_data[:, 2]
#     bin_area = bin_data[:, 3]
#     bin_inuse = bin_data[:, 4].astype(int)

#     ifu_data = np.loadtxt('{}/ifu/IFU_data'.format(path))
#     bin_index = ifu_data[:, 0].astype(int)
#     gh_v0 = ifu_data[:, 1]
#     gh_v0_err = ifu_data[:, 2]
#     gh_vd = ifu_data[:, 3]
#     gh_vd_err = ifu_data[:, 4]
#     gh_h3 = ifu_data[:, 5]
#     gh_h3_err = ifu_data[:, 6]
#     gh_h4 = ifu_data[:, 7]
#     gh_h4_err = ifu_data[:, 8]
#     v0 = ifu_data[:, 9]
#     v0_err = ifu_data[:, 10]
#     vd = ifu_data[:, 11]
#     vd_err = ifu_data[:, 12]
#     metal = ifu_data[:, 13]
#     flux = ifu_data[:, 14]

#     c1 = fits.Column(name='xbin', format='D', array=bin_X)
#     c2 = fits.Column(name='ybin', format='D', array=bin_Y)
#     c3 = fits.Column(name='v0', format='D', array=gh_v0)
#     c4 = fits.Column(name='v0_err', format='D', array=gh_v0_err)
#     c5 = fits.Column(name='vd', format='D', array=gh_vd)
#     c6 = fits.Column(name='vd_err', format='D', array=gh_vd_err)
#     c7 = fits.Column(name='h3', format='D', array=gh_h3)
#     c8 = fits.Column(name='h3_err', format='D', array=gh_h3_err)
#     c9 = fits.Column(name='h4', format='D', array=gh_h4)
#     c10 = fits.Column(name='h4_err', format='D', array=gh_h4_err)
#     c11 = fits.Column(name='metal', format='D', array=metal)
#     c12 = fits.Column(name='flux', format='D', array=flux)

#     c13 = fits.Column(name='rebin_x', format='D', array=spax_bin_X)
#     c14 = fits.Column(name='rebin_y', format='D', array=spax_bin_Y)
#     c15 = fits.Column(name='binid', format='D', array=spax_bin_id)

#     coldefs1 = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12])
#     coldefs2 = fits.ColDefs([c13, c14, c15])
#     hdu = fits.PrimaryHDU()
#     tbhdu1 = fits.BinTableHDU.from_columns(coldefs1)
#     tbhdu2 = fits.BinTableHDU.from_columns(coldefs2)
#     hdulist = fits.HDUList([hdu, tbhdu1, tbhdu2])
#     hdulist.writeto('{}/ifu/IFU_gh.fits'.format(path), clobber=True)

#     c1 = fits.Column(name='xbin', format='D', array=bin_X)
#     c2 = fits.Column(name='ybin', format='D', array=bin_Y)
#     c3 = fits.Column(name='v0', format='D', array=v0)
#     c4 = fits.Column(name='v0_err', format='D', array=v0_err)
#     c5 = fits.Column(name='vd', format='D', array=vd)
#     c6 = fits.Column(name='vd_err', format='D', array=vd_err)
#     c7 = fits.Column(name='metal', format='D', array=metal)
#     c8 = fits.Column(name='flux', format='D', array=flux)

#     c12 = fits.Column(name='rebin_x', format='D', array=spax_bin_X)
#     c13 = fits.Column(name='rebin_y', format='D', array=spax_bin_Y)
#     c14 = fits.Column(name='binid', format='D', array=spax_bin_id)

#     coldefs1 = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8])
#     coldefs2 = fits.ColDefs([c12, c13, c14])
#     hdu = fits.PrimaryHDU()
#     tbhdu1 = fits.BinTableHDU.from_columns(coldefs1)
#     tbhdu2 = fits.BinTableHDU.from_columns(coldefs2)
#     hdulist = fits.HDUList([hdu, tbhdu1, tbhdu2])
#     hdulist.writeto('{}/ifu/IFU.fits'.format(path), clobber=True)

#     # plt.plot(bin_X,bin_Y,'.r')
#     # plt.plot(spax_bin_X,spax_bin_Y,'.k',alpha=0.1)
#     # plt.show()

# if __name__ == '__main__':
#     main()

