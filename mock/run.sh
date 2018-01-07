#!/bin/bash
# Modified from Hongyu Li's script: https://github.com/HongyuLi2016/illustris-tools

# make mock image -p phi -i inclination -o oblate rotatorion
python makeImg.py --subhaloID $1 --rotation $2 --phi $3 --inc $4 #--pa $5
foldername="${1}_p${3}_i$4"
# make stellar mge using the mock image created by illustris-make_img.py
python makeMGE.py --foldername $foldername
# Voronoi bin
python makeIFUbins.py --foldername $foldername
# Calculate Vel, Vdisp and their error in each Voronoi bin
python makeIFUdata.py --foldername $foldername
# dump IFU data to a fits file and plot velocity map
python dumpIFUdata.py --foldername $foldername
python v_plot.py --foldername $foldername
