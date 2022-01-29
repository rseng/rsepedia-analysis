---
title: 'EDO.js: A comprehensive JavaScript library for interaction with musical set theory in any tuning'
tags:
- JavaScript
- music theory
- music analysis
- music perception
- set theory
- diatonic set theory
- equal divisions of the octave
- TET

authors:
- name: Michael Seltenreich
  orcid: 0000-0002-6061-0674
  affiliation: 1

affiliations:
- name: Department of Music, New York University, New York, New York 10003
  index: 1
  
date: 23 September 2021
bibliography: paper.bib

---

# Summary

EDO.js is a JavaScript library for generating, analyzing, and visualizing pitch-collections within any tuning system of Equal 
Divisions of an Octave (EDO). The current focus is on psychological experiments pertaining to musical 
scale perception and pitch-centric music research more broadly. EDO.js implements major algorithms presented in the
music theory and music cognition literature [e.g., @Balzano:1982; @clampitt; @rahn1991coordination; @vassilakis2001auditory; @lerdahl2001tonal; @carey2007coherence], alongside a wide array of additional tools for analyzing and manipulating midi-data and melodies in the framework of musical set theory [@forte1973structure] and diatonic set theory. EDO.js sets itself apart from existing music analysis libraries (such as, music21 [@music21]), in that while most libraries are centered around engaging with the canon of Western music through the lens of applied (and adopted) music-theoretical concepts (e.g., figured bass, roman numerals, etc.), EDO.js is approaching music from a more fundamental level and is aimed at engaging with music's primitives from a cognitive and psycho-acoustical point of view. As such, EDO.js provides analytical tools that can be applied also to non-Western music, making the library more exploratory in nature. Most critically, EDO.js is, to the best of our knowledge, the most comprehensive library to support non-standard tuning to date.

Two classes are at the core of the package. The EDO class, and a daughter class, the Scale class. The EDO class and the Scale class conceptually diverge only on a single aspect: pitches in the Scale class are regarded as pitch-classes (exist only within an octave) and they have an enharmonic meaning. However, in the EDO class pitches are regarded as specific, and are not confined to a single octave. That is to say, while the Scale class treats C3 and C4 simply as “C”, the EDO class would regard them as separate entities. As such, the Scale class is better suited to engage with scales, their properties, and their structures, while the EDO class is an implementation of useful functions for engaging with more abstract musical structures.

### The EDO Class

The class contains eight sets of functions:

 * Functions used for converting between equivalent representations of the input in various formats (cents, ratios, pitches, intervals, frequencies, etc.)
 * Functions used for quantifying various parameters of a given input (e.g., the number of common tones between two inputs).
 * Boolean assertations on input data (e.g., is a set of pitches a transposition of another?) 
 * Functions used for visualization (e.g., contours, fractal trees, necklaces, fractal necklaces, nested necklaces, etc.). 
 * Functions used for importing and processing of midi files.
 * Functions used for importing and processing MusicXml files.
 * Functions used for exporting data in various formats (e.g., images).
 * Other functions for analysis, and manipulation (e.g., extracting motives, extracting contour, generating harmonic progressions, extracting pitch distribution, generating pseudo-random melodies, and more). 

### The Scale Class

The class contains six sets of functions:

 * Functions for converting between equivalent representations of the scale (e.g., for converting pitches to step sizes). 
 * Functions for quantifying various parameters of a given scale (cardinality, trichords, dissonance, coherence, and others).
 * Boolean assertations about the scale (e.g., is it invertible? Is it a mode-of-limited-transposition?)
 * Functions used for exporting (e.g., Scala files). 
 * Functions used for scale structure visualization (e.g., scalar fractal trees). 
 * Analysis, manipulation, and generation functions (e.g., calculating the coherence quotient, extracting diatonic motives, returning the Rothernberg Propriety, calculating Vassilakis Roughness, retrieving available n-chords, etc.). 

In addition to the sets of functions described above, the Scale class also contains five chainable methods commonly used in set-theory.

 * `Scale.invert()` returns the inversion of the original set.
 * `Scale.mode(n)` returns the nth mode of the original set.
 * `Scale.normal()` returns the set in normal order.
 * `Scale.prime()` returns the set in its prime form.
 * `Scale.complement()` returns the complement of the set in the current EDO.

For instance, `Scale.mode(3).invert().complement()` will return an instance of a Scale that is equivalent to the complement scale of the inversion of the 3rd mode of the original scale. 

# Statement of need

`Edo.js` is written in JavaScript for deployment on a Node.js server or in a web-browser. The full documentation for the package is available [here](https://michaelsel.github.io/edoJS), and the package is available on [GitHub](https://github.com/MichaelSel/edoJS) or in Node Package Manager ([NPM](https://www.npmjs.com/package/edo.js)).

This library is aimed at music theorists, musicologists, and cognitive scientists working on pitch-centric musical research, for the creation of stimuli for experiments, and for the analysis of musical structures within and outside of standard tuning systems.

# References
# edoJS

A Set-Theory based JavaScript library pitch manipulation and analysis in any tuning system that is based on equal divisions of the octave (i.e., EDO, and also known as "TET").

As such, it allows to describe collections of pitches in psycho-acoustical terms (e.g., roughness or dissonance), it implements experimental algorithms proposed in the music theory literature, and allows for standard and novel set-theory manipulations.  

This library is aimed at music theorists, musicologists, and cognitive scientists working on musical research, for the creation of stimuli for experiments, and for the analysis of musical structures.

## Installation and Usage

To install npm (its dependencies)
```
npm i edo.js
```





Similarly, you can do any type of import.

#### Import Library
Client-Side
CDN source example - https://www.jsdelivr.com/package/npm/edo.js?version=1.2.14&path=dist
```xhtml
<script src="https://cdn.jsdelivr.net/npm/edo.js@1.2.14/dist/edo.js"></script>```
```
Server-Side
```Javascript
// NOTE: here no relative path so node will use the edojs installed from NPM
const EDO = require("edo.js").EDO;
```

#### Basic Usage
```javascript
let edo = new EDO(12) //create a new EDO context with 12 divisions.
 
//once the object has been created, you can access its functions.

//invert pitches 
edo.get.inversion([0,2,4,5,7,9,11]) //returns [0, 2,  4, 6, 7, 9, 11] 
 
edo.convert.ratio_to_interval(3/2) //returns [7]
 
edo.count.pitches([0, 3, 3, 2, 4, 3, 4]) //returns [[3,3],[4,2], [2,1], [0,1]] 
// (3 appears 3 times, 4 appears 2 times, etc.)
 
edo.is.subset([2,4],[1,2,3,4,5]) //returns true (the set [2,4] IS a subset of [1,2,3,4,5])
```

## Development

To install development dependencies:
```
npm install --include=dev
```

To build the project

```
npm run build
```

To regenerate the docs

```
npm run docs
```

To test: 

```
npm run test
```

## Author
[Michael Seltenreich](http://www.michaelseltenreich.com) 

## License
[GNU AGPLv3](https://choosealicense.com/licenses/agpl-3.0/)

## Some demos
[DEMOS](https://michaelsel.github.io/edoJS/demos/index.html)

## Full Documentation
[Documentation](https://michaelsel.github.io/edoJS/)
A good place to start: 
[EDO Class](https://michaelsel.github.io/edoJS/EDO.html)

`EDO.js` would love some community contributions. The [project](https://github.com/MichaelSel/edoJS) is hosted at GitHub.


## Issue Tickets ##

Currently, there are no special instructions for creating an issue. For bugs, feature suggestions, and general discussions, feel free to create an issue. 


## Submitting Pull Requests ##

Anyone is welcome to contribute. Please open an issue to propose a feature or report a bug before raising a pull request.
At this time, there is not a well-established style-guide. However, make sure any contribution is well documented, and hopefully includes at least one example along with the expected return. 
