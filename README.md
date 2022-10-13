# D-FDRI - an efficient, fast, single-step reconstruction method for single-pixel imaging
A. Pastuszczak, R. Stojek, P. Wróbel, and R. Kotyński, Aug 2021, https://github.com/rkotynski/D_FDRI/


**Single pixel imaging (SPI)** is an indirect image aquisition technique which enables capturing an image using a single detector rather than with a usual camera containing a pixel array and sophisticated imaging optics. SPI is especially usefull for multispectral imaging, or imaging in the IR or THz ranges, imaging through scattering media etc. Usually image reconstruction from the indirect measurement is computationally intensive. It is ofnen based on the methods of compressive sensing (CS).

**D-FDRI** is an acronym for **Differential Fourier Domain Regularized Inversion**,  https://10.1364/OE.433199, (For non-differential FDRI see also https://doi.org/10.1364/OE.26.020009 Note that D-FDRI makes use of FDRI internally but includes a number of modifications concerning both the sampling patterns and the image reconstruction method. D-FDRI and FDRI perform comparably well in a noise-free situation, but D-FDRI has a better noise robustness. It has also other experimental advantages - the reconstructed image is insensitive  to a constant additive bias of the detector and the measured signal has a more uniform distribution so it may be aquired even with a low bit-resolution DAQ using AC-coupling).

The  <code>dfdri</code> Python module is in the file <code>dfdri.py</code>, while example.py is an example showing how to use  <code>dfdri</code> for generating binary patterns, as well as how to reconstruct the image from a compressive measurement using <code>dfdri</code>. The measurement and reconstruction matrices may be easily saved to a ".mat" file and later used with other software, such as Matlab or Octave.

Basic features of D-FDRI are:
1. **Works with highly compressive measurements** - for the compression ratio of CR=2% images with 256x256 pixels may be reconstructed with mean PSNR of 23dB. For CR=10%, PSNR=26.5dB.
1. **Fast reconstruction times**  - image reconstruction is based on a single matrix-vector product, for which optimized routines exist (e.g. openblas, mkl). These packages are often used in other software . In effect, real time image reconstruction is usually possible without GPU acceleration (See e.g. an example with VIS-IR imaging, https://doi.org/10.6084/m9.figshare.14718318). At the same time, preparation of matrices  takes a lot of time and memory but could be done once, and before the actual measuremnt. Image quality is comparable to that obtained with CS but the reconstruction times are typically by two orders of magnitude shorter.
3. **Differential treatment of data** - image reconstruction is based on the first or second order differences of subsequent measurements so the recovered images are invariant to a DC component of the measurement data. This allows to use AC-coupling for the detector DAQ. 
4. **Distributed differential measurement of the zeroth spatial frequency of the image** - the zeroth spatial frequency (mean value) is measured using a sequence of binary patterns. The measurement is differential and more accurate than in a classical approach.
5. **Binary nonorthogonal patterns may be used for sampling** - arbitrary sampling patterns may be used in the measurement. In particular the patterns may be binary (many optical modulators such as DMD are binary) and they do not have to be orthogonal. In the included code the patterns are taken as binarized low spatial frequency basis of the DCT transform but other sets of patterns could be used as well. Image reconstruction makes use of the actual measurements with these patterns. In many other approches several patterns represent a single baisis function, eg. because the basis are real-valued and the actualy displayed patterns are binary.
6. **The DAQ may be operated at a reduced bit-resolution** - by using binary patterns with the average number of pixels in the "on"-state varied randomly within a given range, we keep the detection signal approximatley within an assumed range. This improves the noise robustness of the signal, in particular to discretization noise of the DAQ.
7. **Complementary measurement** - measurements obtained with binary complementary sampling (e.g. using two beams reflected from a single DMD by mirrors in the "on" and "off"-states) may be reconstructed using the same reconstruction matrix. This makes it simple to optical build set-ups with two independent measurement channels working concurrently and eg. measuring images in different spectral ranges or polarizations. 

**Citation:** <em>A. Pastuszczak, R. Stojek, P. Wróbel, and R. Kotyński, "Differential real-time single-pixel imaging with Fourier domain regularization: applications to VIS-IR imaging and polarization imaging," Opt. Express 29, 26685-26700 (2021).  
 https://doi.org/10.1364/OE.433199 (open access)  
**Download:** https://github.com/rkotynski/D_FDRI/  (GPL license)
**See also:** https://github.com/rkotynski/MD-FDRI/

**Contact:** rafal.stojek@fuw.edu.pl, apast@igf.fuw.edu.pl ,  piotr.wrobel@fuw.edu.pl, rafalk@fuw.edu.pl
**Acknowledgement:** National Science Center (Poland), UMO-2017/27/B/ST7/00885 (RS,PW,RK), UMO-2019/35/D/ST7/03781 (AP).  
 
 
**A typical output from the example.py for compression ratios of CR=2% and CR=20%  is:**
![Image reconstruction at CR=20%](reconstr_2.0proc.jpg?raw=true "D-FDRI image reconstruction at CR=2%")
![Image reconstruction at CR=20%](reconstr_20.0proc.jpg?raw=true "D-FDRI image reconstruction at CR=20%")
 
**As another example we also show below the 128x128 pixel cameraman image reconstructions for sampling at CR=8%, 4% and 2% (m=7,p=1,μ=0.7,ϵ=1e-8):**
![128x128 image reconstruction at CR=8%]( reconstr_cameraman128x128_8proc.jpg?raw=true "D-FDRI 128x128 image reconstruction at CR=8%")
![128x128 image reconstruction at CR=4%]( reconstr_cameraman128x128_4proc.jpg?raw=true "D-FDRI 128x128 image reconstruction at CR=4%")
![128x128 image reconstruction at CR=2%]( reconstr_cameraman128x128_2proc.jpg?raw=true "D-FDRI 128x128 image reconstruction at CR=2%")
