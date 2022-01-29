# ReDCM
Dynamic Causal Modelling in R. This package is an R reimplementation of the DCM algorithm found in the SPM toolbox (https://www.fil.ion.ucl.ac.uk/spm/). Current version includes the bilinear deterministic DCM algorithm.<br/>

## Available methods
**Model type:**<br/>
bilinear<br/>
nonlinear (untested)</br>
two-state (untested)</br>
<br/>
**Data fitting**<br/>
deterministic (time-series)<br/>
stochastic (NOT IMPLEMENTED)<br/>
cross-spectra (NOT IMPLEMENTED)<br/>

## Requirements for Linux

R ( >= 3.2.0 )<br/>
CMake ( >= 2.8 )<br/>
GSL (libgsl-dev)

### Required R packages
methods<br/>
R.matlab<br/>
Matrix<br/>
matrixcalc<br/>
exmp<br/>
numDeriv<br/>
pracma<br/>
plyr<br/>
stringr


## Examples

### Running ReDCM on built-in data

    > library('ReDCM')

    > Amx = array( c(1, 1, 0, 1, 1, 1, 0, 1, 1), dim=c(3,3) )
    > Amx
         [,1] [,2] [,3]
    [1,]    1    1    0
    [2,]    1    1    1
    [3,]    0    1    1

    > DCMe = ReDCM_estimate( ReDCM::DCM, A=Amx ) #change priors of matrix A
    > DCMe@Fe
    [1] -3577.417

    > class(DCMe@Ep[[1]])
    [1] "Params"
    attr(,"package")
    [1] "ReDCM"

    > DCMe@Ep@A
               [,1]      [,2]      [,3]
    [1,] 0.50456649 0.3927084  0.000000
    [2,] 0.02413916 0.3689886 -0.305103
    [3,] 0.00000000 0.3392079  0.125141


## Short documentation

### 1. Estimating DCM models
ReDCM can currently estimate models specified in Matlab, using the SPM toolbox (https://www.fil.ion.ucl.ac.uk/spm/) and exported using the **-v6** option, i.e. save('DCM.mat', 'DCM', '-v6').<br/>

**ReDCM_estimate = function(DCM.mat, A=NULL, B=NULL, C=NULL)**<br/>

*Function:*<br/>
Estimates DCM model parameters.<br/>

*Arguments:*<br/>
**DCM.mat** - A .mat file of a DCM structure exported from Matlab, OR a ReDCM structure built with ReDCM_prepare_dcm() function.<br/>
**A** - A numeric array of zeros and ones specifying the endogenous connectivity matrix of DCM. Size of 'A' have to match the size of DCM.n * DCM.n, where DCM.n is the number of regions in the DCM model.<br/>
**B** - A numeric array of zeros and ones specifying the modulatory effects of DCM. Size of 'B' have to match the size of DCM.n * DCM.n * DCM.uN, where DCM.n is the number of regions and DCM.uN is the number of inputs in the DCM model.<br/>
**C** - A numeric array of zeros and ones specifying the direct regional effects of DCM. Size of 'C' have to match the size of DCM.n * DCM.uN, where DCM.n is the number of regions and DCM.uN is the number of inputs in the DCM model.<br/>
<br/>

**ReDCM_prepare_dcm = function(DCM.mat)**<br/>

*Function:*<br/>
Prepares DCM structure for ReDCM_estimate() using a .mat file of a DCM structure exported from Matlab.<br/>

*Arguments:*<br/>
**DCM.mat** - A .mat file of a DCM structure exported from Matlab.<br/>


### 2. Bayesian Model Selection

**ReDCM_run_BMS = function( Fe, method, model.names=NULL )**<br/>

*Function:*<br/>
Performs Bayesian Model Selection by fixed or random effects analysis.<br/>

*Arguments:*<br/>
**Fe** -  An N x M numeric array or matrix containing estimated Bayesian model evidence (in terms of free-energy) for N subjects and M models.<br/>
**method** - A character string selecting BMS method. Use 'ffx' for fixed effects and 'rfx' for random effects analysis.<br/>
**model.names** - An array of M character strings, where M is the number of models. This must match the number of columns in Fe.<br/>
