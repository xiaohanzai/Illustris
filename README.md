# Illustris
Some functions used to handle Illustris and IllustrisTNG data. Especially designed for my prolate JAM project. Some scripts are modified from Hongyu Li's code: 
https://github.com/HongyuLi2016/illustris-tools

This is a package, so please add the corresponding path into $PYTHONPATH.

## Contents:
* `data`: 
  * `StarData.py`: contains a most general class _StarData_.
* `beta`: contains scripts used for measuring the distribution of the anisotropy parameter and velocity ellipsoids.
  * `StarData_Beta.py`: contains class _StarData_Beta_, which used _StarData_ as the base class. It has a method to calculate the velocity dispersion tensor in a given spatial bin, and a method to calculate the global anisotropy. There is also a function _measureV2map_. This function is much faster than _StarData_Beta.measureV2Tensor_ if you want to calculate the dispersion tensor at a grid of linspace or logspace points.
  * `PlotBetaMap.py`: contains class _PlotBetaMap_, which used _StarData_Beta_ as the base class. It has methods to plot beta distributions in the equatorial plane and the meridional plane. Can be run directly as __python PlotBetaMap.py --subhaloID subhalo3 --shape oblate --snapNum 135__.
  * `PlotVEMap.py`: contains class _PlotVEMap_, which used _StarData_Beta_ as the base class. It has a method to plot velocity ellipsoids in the meridional plane. Can be run directly as __python PlotVEMap.py --subhaloID subhalo3 --shape oblate --snapNum 135__.
  * `getV2.py`: this is a script to get the measured dispersion tensors at a grid of linspace points in the meridional plane. Can be run as __python getV2.py --subhaloID subhalo3 --shape oblate --snapNum 135__.
* `mock`: make mock images, MGEs and Voronoi binned IFU data for a subhalo. See Hongyu's code: https://github.com/HongyuLi2016/illustris-tools. You can simply run __python run.py --subhaloID subhalo3 --rotation oblate --phi 0 --inc 90__, or use run.sh to obtain a full output.
  * `makeImg.py`: __python makeImg.py --subhaloID subhalo3 --rotation oblate --phi 0 --inc 90__ will create mock images for a subhalo, after rotating the principal axes cooridinates into axes specified by phi, inc, pa (see Figure 1 of Monnet et al.1992). --rotation = oblate or prolate. "oblate rotation" means z axis is the shortest axis in principal axes coordinates, and "prolate rotation" means z is the longest.  The output files will be in subhalo3/imgs.
  * `makeMGE.py`: __python makeMGE.py --subhaloID subhalo3__ will create MGE using the mock image (img_M.npy) generated by the first script. The output files will be in subhalo3/mge.
  * `makeIFUbins.py`: __python makeIFUbins.py --subhaloID subhalo3__ will create Voronoi bins for the IFU data. The output files will be in subhalo3/ifu.
  * `makeIFUdata.py`: __python makeIFUdata.py --subhaloID subhalo3__ will calculate mean velocity, velocity dispersion and their errors for the particles in each Voronoi bin. The output files will be in subhalo3/ifu
  * `dumpIFUdata.py`: __python dumpIFUdata.py --subhaloID subhalo3__ will dump the IFU data into a fits file IFU.fits in subhalo3/ifu.
  * `v_plot.py`: __python v_plot.py --subhaloID subhalo3__ will read the IFU.fits file and plot the velocity maps. The output figure will be in subhalo3/ifu.
