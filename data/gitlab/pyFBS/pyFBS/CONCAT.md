
---
title: 'pyFBS: A Python package for Frequency Based Substructuring'
tags:
  - Python
  - Structural dynamics
  - Frequency Based Substructuring
  - System Equivalent Model Mixing
  - Transfer Path Analysis
authors:
  - name: Tomaž Bregar
    affiliation: 1 
  - name: Ahmed El Mahmoudi
    affiliation: 2
  - name: Miha Kodrič
    affiliation: 3
  - name: Domen Ocepek
    affiliation: 3
  - name: Francesco Trainotti
    affiliation: 2
  - name: Miha Pogačar
    affiliation: 3
  - name: Mert Göldeli
    affiliation: 2	
  - name: Gregor Čepon
    affiliation: 3
  - name: Miha Boltežar
    affiliation: 3
  - name: Daniel J. Rixen
    affiliation: 2
affiliations:
 - name: Gorenje d.o.o., Partizanska 12, 3503 Velenje, Slovenia
   index: 1
 - name: Technical University of Munich, Institute of Applied Mechanics, Boltzmannstr.  15, 85748 Garching, Germany
   index: 2
 - name: Faculty of Mechanical Engineering, University of Ljubljana, Aškerčeva 6, 1000 Ljubljana, Slovenia
   index: 3
date: 13 August 2017
bibliography: paper.bib
---

# Summary

In science, engineering and technology complex problems are often decomposed into smaller, simpler subsystems. 
Each subsystem can then be analyzed and evaluated separately. 
This approach can often reduce the complexity of the overall problem and provide invaluable insight into the optimization and troubleshooting of each individual component. 
The subsystems can also be assembled back together and with that the system can be analyzed as a whole.

Dynamic Substructuring (DS) is an engineering concept where dynamic systems are modeled and analyzed in terms of their components or so-called substructures. 
There are several ways of formulating the dynamics of substructures. One of them is with Frequency Response Functions (FRFs), which describe the response as the result of a unit harmonic force. 
The method is well suited for experimental approaches where FRFs are obtained from measurement of components. Such approaches were already investigated in the 1970s [@KLOSTERMAN_1971_PHD]  and 1980s [@MARTINEZ_1984_COMBINEDEXPANALYTICAL; @KLOSTERMAN_1984_SMURF; @JETMUNDSEN_1988_FBS; @URGUEIRA_1989_DYNAMIC]. 
Due to complicated formulations and difficulties in obtaining good measurements, the method was hardly applicable. 
Thanks to better measurement hardware and proper formulation of the problem, Frequency Based Substructuring (FBS) has gained popularity in recent years [@deKlerk2008; @vanderSeijs2016; @RIXEN_2006_GUITAR]. 
With this approach, it is also possible to build hybrid models in which experimentally characterized and numerically modelled parts are combined.

pyFBS is a Python package for Frequency Based Substructuring. The package implements an object-oriented approach for dynamic substructuring. 
Current state-of-the-art methodologies in frequency based substructuring are available in pyFBS. Each method can be used as a standalone or interchangeably with others. 
Also a 3D display is available so a user can simply place and orient associated input/outputs used with each method.
Furthermore, basic and application examples are provided with the package together with real experimental and numerical data [@ahmed]. 


# Statement of need

Evaluating structural dynamics is a necessary step in the development of any complex mechanical system. 
Vibro-acoustic character together with the visual design contributes to the customer's perception of a premium product.
With DS the vibro-acoustic performance of a product can be analysed in terms of its subcomponents. This approach is highly beneficial, as almost all complex products are designed modularly.
pyFBS helps the user to perform experimental modelling in DS. It enables an intuitive way to position sensors and impact locations on the analysed structures.
The position and orientation of sensors/impacts can be then be used in each implemented DS method. 
Furthermore, experimental measurements can be simulated with ease from a numerical model, where the same positional information can be used.   

To the best of authors knowledge there is currently no open source software alternative, which would enable the user to use dynamic substructuring methodologies. 
pyFBS has been designed to be used for scientific research in the field of dynamic substructuring. 
It is currently being used by a number of undergraduate students and postgraduate researchers. 


# Features

pyFBS enables the user to use state-of-the-art dynamic substructuring methodologies in an intuitive manner.
Currently implemented features are listed below. 

## 3D display

Structures and positions of impacts, sensors and channels can be visualized in 3D in \autoref{fig:3D}. 
The 3D display is built on top of PyVista [@sullivan2019pyvista] and enables an intuitive way to display relevant data. 
Sensors and impacts can be interactively positioned on structures and the updated positions can be directly used within the pyFBS.
With this feature the experimental setup can be prepared in advance, to avoid possible mistakes in experimental modelling.
Furthermore, various animations can be performed directly in the 3D display, such as the animation of mode shapes or operational deflection shapes.

![An example of a simple structure depicted in the pyFBS 3D display.\label{fig:3D}](./images/figure.png)

## FRF synthetization

Frequency Response Functions can be synthetized for predefined positions of channels and impacts in a numerical model. 
Currently, mode superposition FRF synthetization is supported, where mass and stiffness matrices are imported from FEM software. 
Damping can be introduced as modal damping for each mode shape. 
Additionally, noise can be introduced to the response so a realistic set of FRFs, representing experimental measurements, can be obtained.

## Virtual Point Transformation (VPT) 

VPT projects measured dynamics on the predefined interface displacement modes (IDMs) [@solvingRDOF]. 
The interface is usually considered to be rigid; therefore, only 6 rigid IDMs are used in the transformation. 
After applying the transformation, a collocated set of FRFs is obtained, which can afterwards directly be used in DS. 
Expanded VPT is also supported, where directly measured rotational response is included in the transformation [@Bregar2020].

## System Equivalent Model Mixing (SEMM)

SEMM enables mixing of two equivalent frequency-based models into a hybrid model [@Klaassen2018; @semm_svd]. 
The models used can either be of numerical or experimental nature. One of the models provides the dynamic properties (overlay model) and the second model provides a set of degrees of freedom. 
A numerical model is commonly used as a parent model and an experimental model is used as an overlay model. 

# Discussion

The development of the pyFBS is an ongoing effort. Currently the package is used as a research tool. 
In the future, more examples on the topic of DS and applications of Transfer Path Analysis (TPA) are going to be introduced in the documentation. 
Furthermore, implementation of Operational Source Identification (OSI) is going to be integrated into the pyFBS in the near future.

# Acknowledgments

The pyFBS package was developed as a part of collaboration between the Laboratory for Dynamics of Machines and Structures (LADISK), Faculty of Mechanical Engineering, University of Ljubljana (UL FME) 
and the Chair of Applied Mechanics (AM), Technical University of Munich (TUM).

# References
