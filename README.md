# Single-Pixel Imaging (SPI) at high pixel resolutions (MD-FDRI Python code)
R. Stojek, A. Pastuszczak, P. Wróbel, and R. Kotyński, Oct 2022, https://github.com/rkotynski/MD-FDRI/

**MD-FDRI** stands for **Map-based, Differential, Fourier Domain Regularized Inversion**. MD-FDRI is a framework for Single-Pixel Imaging (SPI) applicable at high resolutions and high compression.
**Single pixel imaging (SPI)** is an indirect image acquisition technique which allows to capture images using a single detector rather than with a usual camera containing a pixel array and sophisticated imaging optics. SPI is especially useful for multi-spectral imaging, or imaging in the IR or THz ranges, imaging through scattering media etc. Usually image reconstruction from the indirect measurement is computationally intensive. It is often based on the methods of compressive sensing (CS).
 MDFDRI code accompanies the Opt. Express 30, 22730, 2022 paper
by R.Stojek, A. Pastuszczak, P. Wróbel and R. Kotyński on single-pixel imaging at high resolutions:

https://doi.org/10.1364/OE.460025

https://doi.org/10.6084/m9.figshare.19863556

The MDFDRI class definition is included in the _mdfdri.py_ file, and 
_example_mdfdri_animation.py_ is an example program which demonstrates the use of the MDFRI class. This example includes a simulation of a sequence of SPI measurements with compressive measurements conducted through a varied aperture. As a result it produces an animation and figures with examples of the sampling functions and image maps which were used.

MDFDRI needs huge (2x9.4GB) image sampling and reconstruction matrices, which we will first try to load from a compressed file in the current directory, secondly to download from a repository, and third to recalculate. We recommend having a fast SDD and 32GB RAM for executing this example program, and 128GB memory and a reasonable swap file for matrix recalculation. When the matrices are loaded to RAM, the SPI image reconstruction at the resolution of 1024x768 and sampling ratio of 0.4% is relatively fast, and depending on the CPU takes between 0.3-1s per frame.

The advantages of MD-FDRI are that it combines several features important for realistic optical SPI at high resolutions. Sampling functions are binary and can be directly displayed on optical modulators such as digital micromirror devices (DMD) at their native resolution. The measurement is inherently differential and the reconstructed images are not affected by a constant bias added to the measurement signal for instance due to a bias signal from background illumination or due to the detector dark current. Sampling is non-adaptive (and adaptive sampling at the DMD operating frequency of >20kHz is in practice problematic). It works at very strong compression (eg. <0.5%) so the image acquisition time may be kept at over 5fps, and still the reconstruction time is below 1s even without GPU acceleration. Extremely strong compression implies poor image quality but the MD-FDRI framework works well after limiting the field of view (with the same nonadaptive sampling sequence kept in each measurement). So we can point the aperture at a point of interest and get an improved image quality at that location. On top of this, the distribution of measurement points within a single measurement assures a better entropy than is obtained with many other methods. This in particular enables using simpler A/D converters (eg. 8bit DACs) to capture the measurement.


Sample 31-level 1024x768 image maps used during the creation of the sampling functions:
![](readme/fig_image_maps.jpg)

Examples of 1024x768 binary sampling patterns used in the compressive measurement:
![](readme/fig_sampling_patterns.jpg)

Animated sequence of SPI measurements conducted at the resolution of 1024x768 with a changing field of view. Reducing the field of view improves the reconstruction quality (reconstruction times are obtained with a i7-11700K CPU.)
![](readme/mdfdri_animation_768_1024.gif)


**Citation:** *R. Stojek, A. Pastuszczak,  P. Wróbel, and R. Kotyński, "Single-pixel imaging at high pixel resolutions," Optics Express, vol. 30(13) , pp. 22730-22745* , 10.1364/OE.460025 (open access)

**Download**: https://github.com/rkotynski/MD-FDRI/ (GPL license) Contact: rafal.stojek@fuw.edu.pl, apast@igf.fuw.edu.pl, piotr.wrobel@fuw.edu.pl,  rafalk@fuw.edu.pl, 

**Acknowledgement**: National Science Center (Poland), UMO-2017/27/B/ST7/00885 (RS,PW,RK), UMO-2019/35/D/ST7/03781 (AP).
