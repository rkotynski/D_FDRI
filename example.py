#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27 Jul 2021, updated Aug 5 2021

This D-FDRI code in Python enables to generate a measurement matrix with sampling patterns and a reconstruction matrix for single pixel imaging (SPI)
The description of D-FDRI is included in the following paper. 
This file shows how to use the dfdri module. It assumes that  sample images are put in a subdirectory tstimg_path.

If you find the code useful, please cite this work:

Citation: [1] A. Pastuszczak, R. Stojek, P. Wróbel, and R. Kotyński, "Differential real-time single-pixel imaging with Fourier domain regularization: applications to VIS-IR imaging and polarization imaging," Opt. Express 29, 26685-26700 (2021).
https://doi.org/10.1364/OE.433199 (open access)
Download: https://github.com/rkotynski/D_FDRI/
Contact (authors): anna.pastuszczak@fuw.edu.pl, rafal.stojek@fuw.edu.pl, piotr.wrobel@fuw.edu.pl, rafal.kotynski@fuw.edu.pl
Acknowledgement: National Science Center (Poland), UMO-2017/27/B/ST7/00885 (RS,PW,RK), UMO-2019/35/D/ST7/03781 (AP).




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


'''
μ,ϵ,p, m : D-FDRI parameters (See sect. 2.2 of ref. [1]).
p=1,2 is the order of finite difference operator D
m is an odd intteger defining the number of pixel areas used to evaluate the zeroth spatial frequency of the image
μ,ϵ describe the properties of generalized inversion 
''' 
m=7
p=1
μ=0.5
ϵ=1e-7
CR=0.03 #Compression ratio. For 256x256 images, CR=0.03 requires 8GB RAM, CR=0.1 requires 32GB RAM, for larger dimensions use >=128GB
dim=(256,256)
fname=f'd_fdri_cr{round(CR*100)}proc_m{m}_p{p}_mu{μ}_eps{ϵ}.mat' # save the matrices to a mat-file


dfdri=DFDRI(m=m, p=p, μ=μ, ϵ=ϵ,CR=CR,dim=dim,verbose=True) # initialize the object and set the deafult parameters
print('\nI. PREPARATION STAGE\n')
if os.path.isfile(fname): # check if the matrices are already precalculated
    print(f'Reading precalculated matrices from file {fname}')
    data=io.loadmat(fname)
    M=data['M']
    Pg=data['P']
    print('Done')
    SM=None
else:    # calculate the measurement and reconstruction matrices 
    M_dct,SM=dfdri.dct_sampling_functions(CR=CR) # calculate the continuous low-frequency DCT patterns
    A=dfdri.auxiliary_matrix_a(m=m,p=p) # find an auxilliary matrix A_m^p,  see sect. 2.1 of ref. [1]
    M=dfdri.binary_measurement_matrix(M_dct,A) # binarize the DCT patterns, and combine all the sampling patterns into the binary measurement matrix M (see sect. 2.1 of [1] and figs. 1-4)
    Pg=dfdri.d_fdri(Mbin=M, p=p,μ=μ,ϵ=ϵ) # calculate the reconstruction matrix using Eq. (6) of Ref. [1] (takes a lot of time)
    print('Saving the measurement and reconstruction matrices to a mat-file (which may be e.g. read in Matlab)')
    io.savemat(fname, {'M':M,'P':Pg,'p':p,'m':m,'eps':ϵ,'mu':μ},do_compression=True,appendmat=True)
    print('Done')
Rec1=dfdri.reconstruct(Pg,channel=1) # create an image  reconstruction function (Eq. (11) of ref. [1] for l=1)
#Rec2=dfdri.reconstruct(Pg,channel=2) # for a complementary measurement (Eq. (11) of ref. [1] for l=2)



if SM is not None: # Plot the positions of selected DCT basis
    fig,ax=plt.subplots(1,1,figsize=(3,3))
    ax.imshow(np.double(SM),origin='upper')
    ax.set_title('Selected DCT patterns')
print("Current Working Directory " , os.getcwd())

tstimg_path='tst_images/' # path to test images
tst_images=['lena512.bmp','bird512.jpg','fox512.gif','FUWchart512.jpg']
print('\nII. COMPRESSIVE MEASUREMENTS\n')

fig,ax=plt.subplots(len(tst_images),3,figsize=(12,5*len(tst_images)))
for testnr in range(len(tst_images)): 
    fname=f'{tstimg_path}{tst_images[testnr]}'
    print(f'1. Preparing scene {fname}\n');
    x = np.array(Image.open(fname).resize(size=dim),dtype=np.single).reshape((-1,1))/255
    print(f'2. Taking the compressive measurement {testnr+1}...\n'); 


    y=M@x# This line of code simulates the compressive measurement
    y+=np.single(np.random.rand()*1000) # an added random bias shows that the reconstruction result does not depend on the detector offset
                                # and that the DAQ may use AC-coupling

    print('2. Reconstructing the image...\n'); 
    t0 = datetime.datetime.now()
    x0=Rec1(y) # This line of code reconstructs the image from the compressive measurement y[testnr] using Eq. (11)
    t1 = datetime.datetime.now()
    dt=(t1-t0).microseconds/1e3
    psnr=dfdri.psnr(x,x0)
    print(f'Reconstruction time: {dt}ms, PSNR={round(psnr,2)}dB\nDone...\n'); 
    
    vmax=max(x0.max(),x.max())
    vmin=min(x0.min(),x.min())
    cmap='bone'#'cubehelix','cividis'#
    ax0=ax[testnr,0]
    ax1=ax[testnr,1]
    ax2=ax[testnr,2]
    ax0.imshow(x.reshape(dim),cmap=cmap,vmin=vmin,vmax=vmax)
    ax0.set_title('Ground truth')
    ax1.semilogx(y,'.b')
    ax1.semilogx(y[:p+m],'.m')
    ax1.set_xlabel(f'{M.shape[0]} binary patterns')    
    ax1.set_title(f'Compr. measurem. CR={round(100*M.shape[0]/M.shape[1],1)}%')
    ax1.set_aspect('auto')
    ax2.imshow(x0.reshape(dim),cmap=cmap,vmin=vmin,vmax=vmax)
    ax2.set_title(f'Reconstructed image\nPSNR={round(psnr,1)}dB\n reconstr. time  dt={round(dt,1)}ms')
fig.savefig(f'reconstr_{round(100*CR,1)}proc.jpg')