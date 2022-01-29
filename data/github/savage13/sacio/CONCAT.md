libsacio Documentation
======================

[![Build Status](https://travis-ci.com/savage13/sacio.svg?branch=master)](https://travis-ci.com/savage13/sacio)

Overview
--------

The sacio library, libsacio, provides an interface to read and write sac files and manipulate their contents.

SAC Files
---------
SAC (Seismic Analysis Code) files are binary or alphanumeric data files for storing time series data, primarily ground motion recorded by seismometers. Format of the SAC file consists of a metadata header followed by the ground motion stored as equally or unequally spaced in time. The data section may also be the Fourier Transform of an input ground motion or a 2D data set.


Documentation
-------------

[libsacio Documention](https://savage13.github.io/sacio/)

IRIS-SAC Relationship
---------------------

The code and library here was originally written and licensed under the [2-Clause BSD License](https://choosealicense.com/licenses/bsd-2-clause/) and subsequently included for the [IRIS-SAC distribution](http://ds.iris.edu/ds/nodes/dmc/forms/sac/). This distribution also includes other, essentially, public domain knowledge on the sac format.  The code found here is the same included within the IRIS-SAC distribution, but does **not** included the licensed code from LLNL/IRIS; this is primarily the processing, filtering, and data operations. If you require the SAC program and this functionality, make that request to [IRIS](http://ds.iris.edu/ds/nodes/dmc/forms/sac/)


Examples
--------

Examples of using the library can be found within the [documenation](https://savage13.github.io/sacio/). Of particular note are:

   - [sac_read()](https://savage13.github.io/sacio/html/group__sac.html#gab2623928ccd5cac1a4cc94be5dd89273)
   - [sac_write()](https://savage13.github.io/sacio/html/group__sac.html#gafe8cd1cadc546ea9d67a5ca14fd3886d)
   - [sac_get_string()](https://savage13.github.io/sacio/html/group__sac.html#ga58de657b18177e79b5fa001b21c55e32) / [sac_set_string()](https://savage13.github.io/sacio/html/group__sac.html#ga9acd9d129945c8fdf5a21434ce1d3359)
   - [sac_get_int()](https://savage13.github.io/sacio/html/group__sac.html#ga65ddb9d01a8d1ea66bc7ad1024b30534) / [sac_set_int()](https://savage13.github.io/sacio/html/group__sac.html#ga65ddb9d01a8d1ea66bc7ad1024b30534)
   - [sac_get_float()](https://savage13.github.io/sacio/html/group__sac.html#ga94ddfd21929cd9cb6faa508d8b6d1460) / [sac_set_float()](https://savage13.github.io/sacio/html/group__sac.html#gaa5cd583512409156b09d6a4b0ec4b683)

Documentation for each function includes a selection of example usage (these also function as tests).

Within the test directory `t` are a set of example programs that use the library.  The best programs to get started are:

   - [t/iotest.c](t/iotest.c) - Reading of all data files and accessing header and data values
   - [t/cut.c](t/cut.c)    - Reading data files with a imposed time window (cut)
   - [t/cutim.c](cutim.c)  - Windowing data within memory (cut)
   - [t/compat.c](t/compat.c) - Usage of the library with original sacio library function names
   - [t/ver.c](t/ver.c)    - Examples of behavior differences between v6 and v7 SAC header files

### Read, Change and Write a sac file

Compile and run the example below using (assuming you have compiled the library):

    gcc -I. -o example example.c libsacio_bsd.a
    ./example t/imp.sac

The file `t/imp.sac` is an impulse function and distributed with `sacio`.

If you have not compiled the library yet, see the [Compiling](#compiling) section below. 
The library does not need to be installed to compile/run the example code, it only needs 
to know where the header file `sacio.h` is.

If you have compiled and installed the library, the example below can be compiled using:

    gcc -I/usr/local/include/sacio -o example  example.c -L/usr/local/lib -lsacio_bsd
    ./example t/imp.sac

```c

#include <stdio.h>
#include <sacio.h>

int
main(int argc, char *argv[]) {
    // Variable Declaration
    int i;
    int nerr = 0;
    char org[128] = {0};
    char id[64] = {0};
    double delta;
    sac *s = NULL;

    if(argc < 2) {
        printf("Usage: example file.sac\n");
        return -1;
    }

    // Read a SAc file, could be evenly or unevenly spaced
    s = sac_read(argv[1], &nerr);
    if(nerr != 0) {
        printf("Error: %d reading file\n", nerr);
        return -1;
    }

    // Get a floating point value: delta
    sac_get_float(s, SAC_DELTA, &delta);

    // Print out some header values
    printf("delta: %f\n", delta);
    printf("npts:  %d\n", s->h->npts);
    printf("comps: %d\n", sac_comps(s));

    // Print out the first 5 data values
    for(i = 0; i < 5; i++) {
        printf("y[%4d]  %f\n", i, s->y[i]);
    }

    // Set some header character string values
    sac_set_string(s, SAC_NET, "IU");
    sac_set_string(s, SAC_STA, "BORG");
    sac_set_string(s, SAC_LOC, "00");
    sac_set_string(s, SAC_CHA, "BHZ");

    // Set the origin time to 1994/160 00:33:16.123
    sac_set_time(s, timespec64_from_yjhmsf(1994, 160, 0, 33, 16, 123));

    // Format the origin time in absolute time
    sac_fmt(org, sizeof(org), "%TO", s);
    printf("Origin Time: %s [Absolute]\n", org);

    // Output the Station id: NET.STA.LOC.CHA
    sac_fmt(id, sizeof(id), "%Z", s);
    printf("Station ID: %s\n", id);

    // Write out the modified sac file
    sac_write(s, "output.sac", &nerr);

    return 0;
}
```


Downloading and installing
--------------------------

### Downloading 

This library can be downloaded directly by either going to **Code->Download Zip** or by using git as:

    git clone https://github.com/savage13/sacio.git

This library is self contained and does not require any additional dependencies. 
 
### Compiling

Once downloaded, the library can be compiled from within the `sacio` or `sacio-master` directory using:
  
    ./configure
    make 
    
### Installation

Installation to the default location `/usr/local` can be completed using:

    make install

This will install `libsacio_bsd.a` into `/usr/local/lib/libsacio_bsd.a` and 
`sacio.h` and `timespec.h` into the `/usr/local/inclucde/sacio` directory.  Creating 
a symbolic link to this library will assist in existing programs that require
the sacio library, e.g.:

    ln -s /usr/local/lib/libsacio_bsd.a /usr/local/lib/libsacio.a

Passing the `--prefix` option to the configure command allows a different installation location.  

### Testing

Tests for the library can be run if desired using

    make test
  
Community Guidelines
--------------------

### Issues / Bugs

Please report issues, possible bugs, inconsistencies, and improvements to the project in the Issue Tracker or using a Pull Request.

### Contributions

If you would like to contribute to the project please file Pull Requests and/or create issues for discussion at the libsacio project.

### Support

The best place to find support for this library and the IRIS version of SAC (and where you can find me), 
is on the [sac-help mailing list](http://ds.iris.edu/message-center/topic/sac-help/).

Included Libraries
------------------

- 64-bit time https://github.com/evalEmpire/y2038 (MIT License)
- GeographicLib https://geographiclib.sourceforge.io/ version 1.49 (MIT/X11 License)

License
-------

The code here is licensed under the 2-Clause BSD License, except where specified.
---
title: 'sacio: A library for Seismic Analysis Code data files'
tags:
  - C
  - Fortran
  - seismology
  - earth science
  - geophysics
authors:
  - name: Brian Savage
    orcid: 0000-0002-9252-6282
    affiliation: 1
affiliations:
  - name: University of Rhode Island
    index: 1
date: 8 July 2021
bibliography: paper.bib
---

# Summary

Since nearly the inception of digital seismological data, the Seismic Analysis Code (SAC) has been utilized as processing software and a file format [@goldstein1995status; @goldstein2003sac2000; @goldstein2005sac]. Having a well-defined format is essential for the distribution and processing of seismological data.  Easy and quick access to seismological data over the internet or within local networks is important for the development of seismological tools that impact hazards (earthquakes, volcanoes, and tsunami), national security, and geophysical imaging on a variety of length scales. Here we implement an improved and open source version of the SAC format with backwards compatibility.

# Statement of need

The sacio library allows reading and writing of binary data files containing time series data, typically ground motion recorded by seismometers.  This library is used in seismological research including, but not limited to, earthquake location [@michael2013precise] and magnitude determination [@miao2007empirical], seismic tomography [@tape2010seismic], earthquake early warning activities [@colombelli2013application], earthquake source process discrimination [@templeton2018seismic], and hazard probability estimates [@al2008improving; @cramer2014update].

A number of deficiencies exist with the original implementation of the routines to read and write SAC binary files, written in Fortran 77 in the 1980s and directly translated into C in the late 1990s.  First, it was tightly coupled with the processing routines included with the SAC program and difficult to use in external programs.  Second, the routines were closed source and are covered by export restrictions by the US Department of Energy (US-DOE) and Lawrence Livermore National Laboratory (LLNL). Third, the routines and associated metadata in the data files did not meet existing needs of current seismological datasets, including very high samples rates, > 100 Hz, long duration data series (multiple days to months), and dense sampling arrays.

This open source library that reads and writes SAC binary and alphanumeric files fixes these outstanding issues.  The library is fully decoupled and released under an open source license, BSD 2-Clause.  It is straightforward to link this library with existing software projects without need to request the closed source, export-restricted version of SAC.

This library was designed as a drop-in replacement with strong backwards compatibility for the original closed source version to facilitate an easy transition for scientists; it is currently included within the official SAC distribution and used by seismologists globally. Moreover, the library was written to allow simple linking from C (for language interoperability) and Fortran, still used by scientists and seismologists.

Finally, the library adds routines to handle 64-bit metadata, or header values, to the binary file format to allow for very high sampling rates, long duration data series, and high precision station locations. 

# Acknowledgments

We acknowledge the contributions of Arthur Snoke for the detailed discussions on version 7 of the header and the support of the SAC user community.

# References
