import numpy as np
import matplotlib
matplotlib.use('Agg')
from JAM.utils import util_fig
import matplotlib.pyplot as plt
import argparse
import pickle
from Illustris.utils.util_shape import plotShape_profile, plotShape_dot
from Illustris.utils.util_shape import plotOrientation_profile, plotOrientation_dot
util_fig.label_font.set_size(15)
util_fig.ticks_font.set_size(10)


def plot_all(Rstar, axisRatiosStar, hmr_star, eigenVectorsStar,
             Rdark, axisRatiosDark, hmr_dark, eigenVectorsDark,
             path):
    # plot shape profile for star
    fig = plt.figure(figsize=(4, 3))
    fig.subplots_adjust(left=0.155, bottom=0.15, right=0.97, top=0.97)
    ax = fig.add_subplot(111)
    plotShape_profile(Rstar, axisRatiosStar, axes=ax, color='r')
    ax.axvline(hmr_star, color='c', lw=2.0)
    ax.axvline(2.0*hmr_star, color='c', lw=2.0)
    ax.set_xlim([0.0, 2.5*hmr_star])
    ax.set_xlabel(r'$R$ [kpc]', fontproperties=util_fig.label_font)
    ax.set_ylabel('Axis ratios', fontproperties=util_fig.label_font)
    util_fig.set_labels(ax)
    fig.savefig('{}/shape_star.eps'.format(path), dpi=500)

    # plot shape dot for star
    fig = plt.figure(figsize=(4, 4))
    fig.subplots_adjust(left=0.14, bottom=0.12, right=0.97, top=0.97)
    ax = fig.add_subplot(111)
    plotShape_dot(Rstar, axisRatiosStar, axes=ax, color='r')
    ax.set_xlabel(r'$b/a$', fontproperties=util_fig.label_font)
    ax.set_ylabel(r'$c/a$', fontproperties=util_fig.label_font)
    lim = np.array([0, 1])
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.plot(lim, lim, '--k', lw=1)
    util_fig.set_labels(ax)
    fig.savefig('{}/shape_star_dot.eps'.format(path), dpi=500)

    # plot shape profile for dark + star
    fig = plt.figure(figsize=(4, 3))
    fig.subplots_adjust(left=0.155, bottom=0.15, right=0.97, top=0.97)
    ax = fig.add_subplot(111)
    plotShape_profile(Rstar, axisRatiosStar, axes=ax, color='r')
    plotShape_profile(Rdark, axisRatiosDark, axes=ax, color='k')
    ax.axvline((hmr_dark), color='c', lw=2.0)
    ax.axvline(2.0*hmr_dark, color='c', lw=2.0)
    ax.set_xlim([0.0, 2.5*hmr_dark])
    ax.set_xlabel(r'$R$ [kpc]', fontproperties=util_fig.label_font)
    ax.set_ylabel('Axis ratios', fontproperties=util_fig.label_font)
    util_fig.set_labels(ax)
    fig.savefig('{}/shape_dark.eps'.format(path), dpi=500)

    # plot shape dot for dark + star
    fig = plt.figure(figsize=(4, 4))
    fig.subplots_adjust(left=0.14, bottom=0.12, right=0.97, top=0.97)
    ax = fig.add_subplot(111)
    plotShape_dot(Rstar, axisRatiosStar, axes=ax, color='r')
    plotShape_dot(Rdark, axisRatiosDark, axes=ax, color='k')
    ax.set_xlabel(r'$b/a$', fontproperties=util_fig.label_font)
    ax.set_ylabel(r'$c/a$', fontproperties=util_fig.label_font)
    lim = np.array([0, 1])
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.plot(lim, lim, '--k', lw=1)
    util_fig.set_labels(ax)
    fig.savefig('{}/shape_dark_dot.eps'.format(path), dpi=500)

    # plot orientation profile for star
    fig, axes = plt.subplots(2, 3)
    plotOrientation_profile(Rstar, eigenVectorsStar,
                             axes=axes, color='r')
    axes[1, 1].set_xlabel(r'$R$ [kpc]', fontproperties=util_fig.label_font)
    axes[0, 0].set_ylabel(r'$\theta$', fontproperties=util_fig.label_font)
    axes[1, 0].set_ylabel(r'$\phi$', fontproperties=util_fig.label_font)
    for ax in axes.reshape(-1):
        util_fig.set_labels(ax)
    for ax in axes[0, :]:
        ax.set_ylim([0, 90])
    for ax in axes[1, :]:
        ax.set_ylim([-180, 180])
    fig.savefig('{}/orientation_star.eps'.format(path), dpi=500)

    # plot orientation dot for star
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    plotOrientation_dot(Rstar, eigenVectorsStar,
                         axes=axes, color='r')
    axes[1].set_xlabel(r'$\phi$', fontproperties=util_fig.label_font)
    axes[0].set_ylabel(r'$\theta$', fontproperties=util_fig.label_font)
    for ax in axes[:]:
        ax.set_ylim([0, 90])
        ax.set_xlim([-180, 180])
        util_fig.set_labels(ax)
    fig.subplots_adjust(left=0.04, bottom=0.13, right=0.97, top=0.97)
    fig.savefig('{}/orientation_star_dot.eps'.format(path), dpi=500)

    # plot orientation profile for dark + star
    fig, axes = plt.subplots(2, 3)
    plotOrientation_profile(Rstar, eigenVectorsStar,
                             axes=axes, color='r')
    plotOrientation_profile(Rdark, eigenVectorsDark,
                             axes=axes, color='k')
    axes[1, 1].set_xlabel(r'$R$ [kpc]', fontproperties=util_fig.label_font)
    axes[0, 0].set_ylabel(r'$\theta$', fontproperties=util_fig.label_font)
    axes[1, 0].set_ylabel(r'$\phi$', fontproperties=util_fig.label_font)
    for ax in axes.reshape(-1):
        util_fig.set_labels(ax, xrotate=45)
    for ax in axes[0, :]:
        ax.set_ylim([0, 90])
    for ax in axes[1, :]:
        ax.set_ylim([-180, 180])
    fig.savefig('{}/orientation_dark.eps'.format(path), dpi=500)

    # plot orientation dot for dark + star
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    plotOrientation_dot(Rstar, eigenVectorsStar,
                         axes=axes, color='r')
    plotOrientation_dot(Rdark, eigenVectorsDark,
                         axes=axes, color='k')
    axes[1].set_xlabel(r'$\phi$', fontproperties=util_fig.label_font)
    axes[0].set_ylabel(r'$\theta$', fontproperties=util_fig.label_font)
    for ax in axes[:]:
        ax.set_ylim([0, 90])
        ax.set_xlim([-180, 180])
        util_fig.set_labels(ax)
    fig.subplots_adjust(left=0.04, bottom=0.13, right=0.97, top=0.97)
    fig.savefig('{}/orientation_dark_dot.eps'.format(path), dpi=500)


def run(path):
    '''
    Save all figs to "path".
    '''
    with open('{}/shape.dat'.format(path), 'rb') as f:
        data = pickle.load(f)
    Rstar = data['Rstar']
    hmr_star = data['hmr_star']
    axisRatiosStar = data['axisRatiosStar']
    eigenVectorsStar = data['eigenVectorsStar']
    Rdark = data['Rdark']
    hmr_dark = data['hmr_dark']
    axisRatiosDark = data['axisRatiosDark']
    eigenVectorsDark = data['eigenVectorsDark']

    plot_all(Rstar, axisRatiosStar, hmr_star, eigenVectorsStar,
             Rdark, axisRatiosDark, hmr_dark, eigenVectorsDark,
             path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--foldername')
    args = parser.parse_args()

    foldername = args.foldername
    path = outpath + foldername

    run(path+'/shape/')


if __name__ == '__main__':
    main()
