
# MATLAB software for image processing

[![DOI](https://zenodo.org/badge/97616448.svg)](https://zenodo.org/badge/latestdoi/97616448)

The software has been developed at tested on MATLAB 9.2, Release 2017a.

## MATLAB Software implementations of the Shape and Affine Invariant descriptor

For example usage refer to the MATLAB scripts in directory 'Testing\comparision'.

Before running any tests, the test datasets 'Affine regions' and 'OxFrei' should be downloaded as explained in 
[this README.md file](https://github.com/NLeSC/LargeScaleImaging/tree/master/Data)

The funciton 'config.m' contains important common parameters for the rest of  the software. 

### Reference
Ranguelova, E., [Local Shape and Moment Invariant Descriptor for Structured Images](http://eprints.maynoothuniversity.ie/8841/1/IMVIP2017_Proceedings.pdf), Proceedings of the 19th Irish Machine Vision and Image Processing (IMVIP) conference, Maynooth, Ireland, 2017, pp. 245-248

### Also of interest
[Salient detectors-matlab](https://github.com/NLeSC/SalientDetector-matlab)



The data used to test this software are the structured scenes of the test data from 

http://www.robots.ox.ac.uk/~vgg/research/affine/

where you can find all data and their full description.

The original files are in PPM format named 'image1..6.ppm' and the homographies between each image and the original one are given.

Here the same data are converted to PNG format and named after the scene sequence, i.e. 'boat1..6.png' for the 'boat' sequence.

Each sequence is used to test a single transformaiton. The mapping sequence <-> transformaiton is:

graffiti <----> viewpoint
bikes    <----> blur
boat     <----> zoom + rotation
leuven   <----> lighting
The OxFrei data can be obtained at

https://github.com/NLeSC/LargeScaleImaging/tree/master/Data/OxFrei
