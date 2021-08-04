# D-FDRI - an efficient, fast, single-step reconstruction method for single-pixel imaging



**Single pixel imaging (SPI)** is an indirect image aquisition technique which enables to capture an image using a single detector rather than with a usual camera containing a pixel array and sophisticated imaging optics. SPI is especially usefull for multispectral imaging, or imaging in the IR or THz ranges, imaging through scattering media etc. Usually image reconstruction from the indirect measurement is computationally intensive. It is ofnen based on the methods of compressive sensing.

**D-FDRI** is an acronym for **Differential Fourier Domain Regularized Inversion**,  https://10.1364/OE.433199, (For non-differential FDRI see also https://doi.org/10.1364/OE.26.020009).

The  <code>dfdri</code> Python module is in the file <code>dfdri.py</code>, while example.py is an example showing how to use  <code>dfdri</code> for generating binary patterns, as well as how to reconstruct the image from a compressive measurement using <code>dfdri</code>. The measurement and reconstruction matrices may be easily saved to a ".mat" file and later used with other software, such as Matlab or Octave.

Basic features of D-FDRI are:
1. **Works with highly compressive measurements** - for the compression ratio of CR=2% images with 256x256 pixels may be reconstructed with mean PSNR of 23dB. For CR=10%, PSNR=27dB.
1. **Fast reconstruction times**  - image reconstruction is based on a single matrix-vector product, for which optimized routines exist (e.g. openblas, mkl). These packages are often used in other software . In effect, real time image reconstruction is usually possible without GPU acceleration (See e.g. an example with VIS-IR imaging, https://doi.org/10.6084/m9.figshare.14718318). At the same time, preparation of matrices  takes a lot of time and memory but could be done once, and before the actual measuremnt.
3. **Differential treatment of data** - image reconstruction is based on the first or second order differences of subsequent measurements so the recovered images are invariant to a DC component of the measurement data. This allows to use AC-coupling for the detector DAQ. 
4. **Distributed differential measurement of the zeroth spatial frequency of the image** - the zeroth spatial frequency (mean value) is measured using a sequence of binary patterns. The measurement is differential and more accurate than in a classical approach.
5. **Binary nonorthogonal patterns may be used for sampling** - arbitrary sampling patterns may be used in the measurement. In particular the patterns may be binary (many optical modulators such as DMD are binary) and they do not have to be orthogonal. In the presented code the patterns are taken as binarized low spatial frequency basis of the DCT transform.
6. **The DAQ may be operated at a reduced bit-resolution** - by using binary patterns with the average number of pixels in the "on"-state varied randomly within a given range, we keep the detection signal approximatley within an assumed range. This improves the noise robustness of the signal, in particular to discretization noise of the DAQ.

**Citation:** <em>A. Pastuszczak, R. Stojek, P. Wróbel, and R. Kotyński, "Differential real-time single-pixel imaging with Fourier domain regularization: applications to VIS-IR imaging and polarization imaging," Opt. Express 29, 26685-26700 (2021).  
 https://doi.org/10.1364/OE.433199 (open access)  
**Download code:** https://github.com/rkotynski/D_FDRI/
  
-- A typical output from the example.py file is:

Differential Fourier Domain Regularized Inversion (D-FDRI)  
This is a python module for compressive differential Single Pixel Imaging (SPI).  
It calculates binary sampling patterns (a binary measurement matrix) and reconstructs images from compressive measurements at a very high speed.  



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
