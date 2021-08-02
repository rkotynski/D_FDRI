#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27 Jul 2021

This D-FDRI code in Python enables to generate a measurement matrix with sampling patterns and a reconstruction matrix for singlepixel imaging (SPI)
The description of D-FDRI is included in the following paper. If you find the code useful, please cite this work:

Anna Pastuszczak, Rafał Stojek, Piotr Wróbel, and Rafał Kotyński, "Differential real-time single-pixel imaging with Fourier domain regularization - applications to VIS-IR imaging and polarization imaging,"
Optics Express, https://doi.org/10.1364/OE.433199 (open accesss)


@author: anna.pastuszczak@fuw.edu.pl, rafal.stojek@fuw.edu.pl, piotr.wrobel@fuw.edu.pl, rafal.kotynski@fuw.edu.pl


  GPL LICENSE INFORMATION

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import io
import datetime
from PIL import Image
from dfdri import DFDRI
from copy import copy

m=7
p=1
μ=0.5
ϵ=1e-8
CR=0.03 # for 256x256 images, CR=0.03 requires 8GB RAM, CR=0.1 requires 32GB RAM, etc
dim=(256,256)
fname=f'd_fdri_cr{round(CR*100)}proc_m{m}_p{p}_mu{μ}_eps{ϵ}.mat'


dfdri=DFDRI(m=m, p=p, μ=μ, ϵ=ϵ,CR=CR,dim=dim,verbose=True)
print('\nI. PREPARATION STAGE\n')
if os.path.isfile(fname):
    print(f'Reading precalculated matrices from file {fname}')
    data=io.loadmat(fname)
    M=data['M']
    Pg=data['P']
    print('Done')
else:    
    M_dct,SM=dfdri.dct_sampling_functions(CR=CR) # calculate the continuous low-frequency DCT patterns
    A=dfdri.auxiliary_matrix_a(m=m,p=p) # find an auxilliary matrix A_m^p,  see https://doi.org/10.1364/OE.433199 
    M=dfdri.binary_measurement_matrix(M_dct,A) # binarize the DCT patterns, and combine all the sampling patterns into the binary measurement matrix M
    Pg=dfdri.d_fdri(Mbin=M, p=p,μ=μ,ϵ=ϵ) # calculate the reconstruction matrix (takes a lot of time)
    
    print('Saving the measurement and reconstruction matrices to a mat-file (which may be e.g. read in Matlab)')
    io.savemat(fname, {'M':M,'P':Pg,'p':p,'m':m,'eps':ϵ,'mu':μ},do_compression=True,appendmat=True)
    print('Done')
Rec1=dfdri.reconstruct(Pg,channel=1)
Rec2=dfdri.reconstruct(Pg,channel=2)



print("Current Working Directory " , os.getcwd())

tstimg_path='tst_images/' # path to test images
tst_images=['lena512.bmp','bird512.jpg','fox512.gif','FUWchart512.jpg']

print('\nII. COMPRESSIVE MEASUREMENTS\n')


for testnr in range(len(tst_images)): 
    fname=f'{tstimg_path}{tst_images[testnr]}'
    print(f'1. Preparing scene {fname}\n');
    x = np.array(Image.open(fname).resize(size=dim),dtype=np.single).reshape((-1,1))/255
    print(f'2. Taking the compressive measurement {testnr+1}...\n'); 


    y=M@x# This line of code simulates the compressive measurement
    y+=np.single(np.random.rand()*1000) # an added random bias shows that the reconstruction result does not depend on the detector offset

    print('2. Reconstructing the image...\n'); 
    t0 = datetime.datetime.now()
    x0=Rec1(y) # This line of code reconstructs the image from the compressive measurement y[testnr]
    t1 = datetime.datetime.now()
    dt=(t1-t0).microseconds/1e3
    psnr=dfdri.psnr(x,x0)
    print(f'Reconstruction time: {dt}ms, PSNR={round(psnr,2)}dB\nDone...\n'); 
    
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(8,4))

    ax1.semilogx(y,'.b')
    ax1.semilogx(y[:p+m],'.m')
    ax1.set_xlabel(f'{M.shape[0]} binary patterns')
    ax1.set_title(f'Compressive measurement\nCR={round(100*M.shape[0]/M.shape[1],1)}%')
    ax2.imshow(x0.reshape(dim),cmap='gray')
    ax2.set_title(f'Reconstructed image\nPSNR={round(psnr,1)}dB\n reconstr. time  dt={round(dt,1)}ms')
 