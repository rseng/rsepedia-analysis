# AstroData

Set of C++ classes to operate on radio astronomical data.

# Installation

Set the `INSTALL_ROOT` environment variable to the location where you want the library and include files installed.
If this package is installed in `$HOME/Code/APERTIF/AstroData` this would be:

```bash
 $ export INSTALL_ROOT=$HOME/Code/APERTIF
```

Then build as follows:

```bash
 $ make install
```

## Dependencies

 * [utils](https://github.com/isazi/utils) - master branch

## Optional dependencies

 * [HDF5](https://support.hdfgroup.org/HDF5/)
 Set the LOFAR environment variable to *true*.

 * [PSRDADA](http://psrdada.sourceforge.net/)
 Set the PSRDADA environment variable to location of the psrdada source, preferably `${INSTALL_ROOT}/psrdada`.

# Included classes

## ReadData.hpp

Data io functions:

 * *readZappedChannels* Zapped channels (excluded from computation)
 * *readIntegrationSteps* Integration steps
 * *readSIGPROC* SIGPROC data
 * *readLOFAR* LOFAR data
 * *readPSRDadaHeader* PSRDADA buffer
 * *readPSRDada* PSRDADA data

## Platform.hpp

Classes and readers for:

 * *paddingConf* Memory padding
 * *readPaddingConf* 
 * *vectorWidthConf* Vector unit width
 * readVectorWidthConf

## Observation.hpp

A class to hold physical observations parameters and search configuration.

## Generator.hpp

Generator for fake data, useful for for testing.

 * *generatePulsar* Generates a periodic single signal, not too relastic.
 * *generateSinglePulse* Generates a single pulse

# License

Licensed under the Apache License, Version 2.0.
