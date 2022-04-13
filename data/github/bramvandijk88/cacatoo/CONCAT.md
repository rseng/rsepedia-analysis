# **CACATOO** [![DOI](https://zenodo.org/badge/317027997.svg)](https://zenodo.org/badge/latestdoi/317027997) [![DOI](https://joss.theoj.org/papers/10.21105/joss.03948/status.svg)](https://doi.org/10.21105/joss.03948)

Cacatoo is a toolbox for modelling and exploration of spatially structured simulations. Because it is written in 100% javascript it requires no installation and works on any machine, making building, sharing, and exploring your model easier than it ever was! With a web-based interface that is ideal for students to learn how to program, and a NodeJS-mode which allows one to rapidly run simulations directly from the command line, it is suited for beginners and advanced programmers alike!

<center>
<img src="https://bramvandijk88.github.io/cacatoo/images/elephant_cacatoo.png" width="400"
     alt="Riding on the shoulders of giants"
/></center>

*Cacatoo was peer reviewed and published in JOSS on February 6th, 2022*


## **Getting started**

You don't need to install anything. You can either immediately start playing with one of the many [JSFiddle examples](https://bramvandijk88.github.io/cacatoo/examples_jsfiddle.html), or dowload this repository and explore dozens of more examples! 
If you want help with Javascript, I recommend [this](https://youtu.be/W6NZfCO5SIk) tutorial from "Programming with Mosh". The first hour is free, and the other 6 hours are definitely worth the money. I also recommend setting up a nice coding environment, as I explain in [this blog post](https://www.bramvandijk.com/blog/2020/11/20/javascript-programming-part-ii-my-setup).


## **How to Cacatoo**
If you prefer to learn by example, I have made dozens of different models for you to start playing (found in the "examples" folder of the repository), ranging from beginner to expert! That said, Cacatoo has been extenstively documented. Tutorials and overviews of its many functions can be found [here](https://bramvandijk88.github.io/cacatoo). Briefly, setting up a Cacatoo model consists of roughly three steps: 1) setup, 2) defining the rules, and 3) setting up a main simulation loop (see image below). I have added comments to some of the example files, denoting which step is to be executed where in the code. If you want to get in debt on what each function specifically does, check out the [JSDocs](https://bramvandijk88.github.io/cacatoo/jsdocs/index.html).

<center>
<img src="images/cacatoo_recipe.png" width="700"
     alt="The basic recipe of a Cacatoo simulation contains three ingredients: 1) setup, 2) defining the rules, and 3) setting up the main simulation loop."
/></center>

<br><br>

# **Motivation**

Complex systems like microbial communities, the human immune system, and the Earth's climate, have many emergent properties which arise from the interactions between individual components. As such, predicting exactly how these systems will behave and respond to various stimuli is difficult. Simulation offers a solution by allowing a modeller to simply put in what they deem important, and observe the outcome. As such, direct visual feedback is an important part of the process of "getting to know" your model. Not only are you more likely to detect programming mistakes, but it also aids exploration of parameter space. Over the last decade, most of my models were implemented by using the CASH C-library, which would allow direct visual feedback by means of the X11 library (R.J. de Boer & A.D. Staritsk). However, as it is based in C, developing models with CASH requires a lot of programming experience, and even an experienced user like myself can sometimes take days to track down a simple bug. Moreover, sharing your model with other users can be a pain in the neck, as installation is slightly different depending on the operating system. The amazing toolbox [Artistoo](https://artistoo.net/) has illustrated how Javascript can be both beginner-friendly and versatile. Hence, I decided that Javascript was the way forward! <br><br>

# **FAQ**

**Q**: My model isn't running, what can I do to find out what's wrong?\
**A**: Open the developer console (CTRL+SHIFT+I in Google Chrome), and see which line of your code causes the problem. If it is an error in Cacatoo, please submit an issue on the Github repository. <br> <br>

**Q**: Grid points / individuals are dissappearing every update!\
**A**: By far the most common reason for this is the synchronous updating of the grid. When synchronously updating the grid,
make sure only the focal cell is modified by the nextState function. If not, the order of updating may cause individuals to appear / dissapear. If you want the nextState to be able to modify neighbouring grid points, only use *asynchronous* updating (see examples). <br><br>

**Q**: My model is too slow, how do I find what is slowing things down?\
**A**: Open the developer console (CTRL+SHIFT+I in Google Chrome), and go to the 'profile' tab. You can run the code for a bit and observe which functions are taking most time. As a general rule, try and avoid build-in functions of javascript (e.g. 'reduce' to get the sum of an array), because despite being easier to use, they are quite a bit slower than a C-style piece of code where you loop over the values yourself. After you understand the model, and you simply want to rapidly run many simulations, I recommend running the code from the command line with nodeJS (also see the "cheater" example).  <br> <br>
**Q**: Why is this toolbox called Cacatoo?\
**A**: Cacatoo is an acronym for CAsh-like Cellular Automaton TOOlbox. However, as Cacatoo developed it became much more than just a tool for cellular automata. The name stuck, however. <br><br>

# **Contributing to Cacatoo**

As the sole developer of Cacatoo, I am eager to get help, suggestions, or use cases from others. For organisational reasons, I encourage everyone to always submit an official Github issue, and be sure to include the following details:

## **Reporting bugs**

When describing a bug, make sure to:
* Describe the unexpected behaviour
* Describe the desired behaviour
* Include a reproducable example

## **Suggesting additions**

When making suggestions, make sure to:
* Describe the missing feature
* Describe or sketch the expected output

## **Pull requests**

When you want to actively contribute to Cacatoo, you can suggest to become a collaborator on Github.
Implemented new branches may be merged with the main branch if the changes are expected to be useful to other. 

* Please give an extensive description of what your branch adds to Cacatoo
* Please run the default unit test (see below) and add the output to the pull request. 
* (optional, but appreciated) Please add a branch-specific unit test for your branch to the unit_test directory

## **Mocha unit testing**

If you want to test your code, or want to issue a pull request, please use Mocha to run the tests provided in the directory unit_test:

* Install mocha (npm install mocha)
* E.g. run the default tests (./node_modules/mocha/bin/mocha unit_test/)

## **Submitting JS fiddle examples**

If you made a nice Cacatoo model which you would like to see on the [JS fiddle examples page](https://bramvandijk88.github.io/cacatoo/examples_jsfiddle.html), be sure to:
* Make sure your code works as intended in JS fiddle (should not require any rewriting, just some copy-pasting)
* Give a title and description of your model 
* Include your name so I can credit you

## **Useful developer commands for contributors**


### **Installing cacatoo locally**

To install cacatoo locally, make sure you have node and npm installed. Then simply run 'npm install' from the main cacatoo directory. Alternatively, you can run 'npm install cacatoo' without cloning all the code first, allowing
you to use cacatoo in your other projects. 

If you have run the above command, a directory "node_modules" should appear which contain the developer dependencies that are explained below. 

### **Bundling the cacatoo package**

The source code in the directory 'src' can be bundled into a single package. This is most convenient when testing your code, and quickly switching between your HTML page and code. The published package is bundled with rollup, which I recommend for these purposes. Below, I illustrate how to use the locally installed developer tools to bundle and test the code. However, you can also do global installs of these tools if you prefer.

To bundle the package using rollups, run:

> ./node_modules/rollup/dist/bin/rollup src/model.js -o dist/cacatoo.js -f cjs  -w

When rollup is installed globally, you can ommit the path and use "rollup src/model.js -o dist/cacatoo.js -f cjs -w" instead. (for a global install of rollup, run npm install -g rollup)

Documentation was compiled with the locally provided version of jsdoc 
> ./node_modules/.bin/jsdoc dist/cacatoo.js -d docs/jsdocs
(alternatively, for a global install, run npm install jsdoc -g)

Unit testing is done with the locally provided version of Mocha 
>./node_modules/mocha/bin/mocha unit_test/
(alternatively, for a global install, run npm install -g mocha)


## **License**
This library is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation. 

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

---
title: 'Cacatoo: building, exploring, and sharing spatially structured models of biological systems'
tags:
  - javascript
  - spatial structure
  - dynamics
  - individual-based models
  - microbial ecology and evolution
authors:
  - name: Bram van Dijk
    orcid: 0000-0002-6330-6934
    affiliation: 1 
affiliations:
 - name: Max Planck Institute for Evolutionary Biology
   index: 1 
date: 11 November 2021
bibliography: paper.bib
---


# Summary

In the field of ecology and evolution, a treasure trove of data has revealed the importance of spatial structure and biogeography. Despite these rich data sets, our conceptual understanding of how spatial structure shapes biodiversity, pathogenicity, and microbial pangenomes is lagging behind. For example, we only have a limited understanding of how interactions at the microscale (molecular machinery, bacteriophages, metabolism) scale up to define eco-evolutionary dynamics of microbial communities [@rainey2020toward] and metaorganisms [@jaspers2019resolving]. To develop our intuition on these systems, I argue we need to embrace multiple levels and structural complexity in our models. `Cacatoo` is a toolbox developed to make it easy to design, explore, and share simulations of such systems. Simulations can be interactively explored from a web browser, allowing the user to change parameters and observe graphs in real-time. `Cacatoo` is designed to be both easy-to-use and extendable, making it suitable for beginners and experts alike. Because it requires no installation and works on practically every computer, it is also ideal for teaching purposes and student projects. In summary, `Cacatoo` provides opportunities for everyone to get involved in spatially structured modelling.

# Statement of need

Complex systems like microbial communities have many emergent properties which arise from the interactions between individual components. As such, predicting exactly how these systems will behave and respond to various stimuli is difficult. Simulation offers a solution by allowing a modeller to simply put in what they deem important, and observe the outcome. As such, direct visual feedback is an important part of exploring the model. Not only is one more likely to detect programming mistakes this way, but it also aids in communicating the models to peers and the general public. Many current modelling frameworks require knowledge of C (e.g. [Cash](https://tbb.bio.uu.nl/rdb/software.html) and [Morpheus](https://academic.oup.com/bioinformatics/article-abstract/30/9/1331/234757) [@starruss2014morpheus]), motivated by its unparalleled speed. However, the learning curve for programming in C is steep, and even an experienced user may take days to track down a simple bug. Moreover, sharing your model with other users can be a pain, as installation is slightly different depending on each operating system. With Javascript tools like d3.js [@zhu2013data] and Artistoo [@wortel2021artistoo] paving the way, `Cacatoo` was developed to overcome these issues, making spatially structured models easy, fast, sharable, and customisable. Other than Artistoo, which first and foremost specialises in tissue simulations, `Cacatoo` provides easy-to-use functions for individual-based modelling  (roulette wheel selection, random walks, neighbourhood retrieval), reaction-diffusion systems (ordinary differential equations), and tools to build your own user interface (displays, graphs, sliders, interaction with mouse and keyboard, *etc.*). The basic recipe of a Cacatoo simulation is simple (\autoref{fig:recipe}), and the extensive documentation provides ample opportunity for students and advanced modellers to get started right away! <br><br>

![The basic recipe of a Cacatoo simulation contains three ingredients: 1) setup, 2) defining the rules, and 3) setting up the main simulation loop.\label{fig:recipe}](../images/cacatoo_recipe.png)

# Use cases

Potential use cases for Cacatoo range from exploring the consequences of [mutations in space](https://bramvandijk88.github.io/cacatoo/example_mutational_jackpot.html) [@fusco2016excess] to setting up a multi-level eco-evolutionary system where [selfish genetic elements co-evolve with their cellular hosts](https://bramvandijk88.github.io/cacatoo/TEs_streamlining/). The latter model is published as part of the special issue "The secret live of mobile genetic elements" [@hall2022introduction],  which revealed that transposable elements can promote genome streamlining [@van2021transposable]. Fullmer is currently in the process of exploring the impact of horizontal gene transfer on the black queen hypothesis [@fullmer2015pan].

# Financial support / Acknowledgements

BD acknowledges support from the Deutsche Forschungsgemeinschaft (DFG) Collaborative Research Center 1182 ‘Origin and Function of Metaorganisms’ (grant no. SFB1182, Project C4 to P.R.).

The author acknowledges Inge Wortel and Johannes Textor for paving the way with their toolbox Artistoo, and for helping me with the early stages of development. Next, I would like to thank Jeroen Meijer for testing and debugging the library.

# References
