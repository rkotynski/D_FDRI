# D_FDRI
D-FDRI (Differential Fourier Domain Regularized Inversion) - an efficient, fast, single-step reconstruction method for single-pixel imaging (DOI: 10.1364/OE.433199)

This is a pre-release. The Python module is in the file dfdri.py, while example.py is an example showing how to use it for generating binary patterns to be used with SPI, as well as how to reconstruct the image from a compressive measurement. The measurement and reconstruction matrices may be easily saved and used with other software, sucha as Matlab, afterwards. Further examples and documentation will be included soon.


-- A typical output from the example.py file is:
Differential Fourier Domain Regularized Inversion (D-FDRI)
This is a python module for compressive differential Single Pixel Imaging (SPI).
It calculates binary sampling patterns (a binary measurement matrix) and reconstructs images from compressive measurements at a very high speed.
For details please see https://doi.org/10.1364/OE.433199 (open access)
If you find this code useful, please cite the paper in your work.


Default parameters:
dim=(256, 256)	-image resolution
p=1				-order of the finite difference operator
m=7				- (p+m) is the number of binary patterns used to measure the 0th spatial freq.
μ=0.5, ϵ=1e-08	-FDRI parameters
CR=0.03			-compression ratio

I. PREPARATION STAGE

Reading precalculated matrices from file d_fdri_cr3proc_m7_p1_mu0.5_eps1e-08.mat
Done
Current Working Directory  /home/rkotynski/Programy/D-FDRI (github)

II. COMPRESSIVE MEASUREMENTS

1. Preparing scene tst_images/lena512.bmp

2. Taking the compressive measurement 1...

2. Reconstructing the image...

Reconstruction time: 31.703ms, PSNR=24.06dB
Done...

1. Preparing scene tst_images/bird512.jpg

2. Taking the compressive measurement 2...

2. Reconstructing the image...

Reconstruction time: 32.41ms, PSNR=23.08dB
Done...

1. Preparing scene tst_images/fox512.gif

2. Taking the compressive measurement 3...

2. Reconstructing the image...

Reconstruction time: 28.325ms, PSNR=24.32dB
Done...

1. Preparing scene tst_images/FUWchart512.jpg

2. Taking the compressive measurement 4...

2. Reconstructing the image...

Reconstruction time: 31.005ms, PSNR=16.6dB
Done...
