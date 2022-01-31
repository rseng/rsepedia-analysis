cRacklet: a spectral boundary integral method library for interfacial rupture simulation
===========================================================================================

cRacklet is a C++ boundary integral library based on a spectral formulation of the dynamic wave
equations in two semi-infinite linearly elastic solids. Its implementation is specially tailored
for the modeling of dynamic crack/rupture propagation along a planar interface bonding the two
half-spaces. The main benefit of this spectral method is a numerical discretization limited to
the mid-plane, thereby providing a very fine description of the dynamic rupture processes,
unattainable with finite-element or finite-difference schemes.
For more details about the method and its applications, refer to "related publications" section.

## Dependencies

The following dependencies are required to compile cRacklet:

- a C++ compiler
- [cMake](https://cmake.org/)
- [FFTW3](http://www.fftw.org/)
- [GSL](https://www.gnu.org/software/gsl/)

Optional dependencies are:

- a Fortran compiler to generate new kernels
- [python 3+](https://www.python.org/), for python binding.
- [pybind11](https://github.com/pybind/pybind11), for python binding, automatically installed if not found on the system.
- FFTW3 with OpenMP support, for multi-threaded parallel computing
- [pytest](https://docs.pytest.org/en/latest/) (for tests)
- [Doxygen](http://doxygen.nl/) (for documentation)
- [Sphinx](https://www.sphinx-doc.org/en/stable/) (for documentation)

## Compile and build

You can compile cRacklet using cMake:

    mkdir build
    cd build
    ccmake .. or cmake ..
    make

## Tests

You need to activate the test options with cMake (`CRACKLET_TESTS`). Go inside `build/` and then run the test with:

    make test

## Documentation

Online documentation is available at https://cracklet.gitlab.io/cracklet/ .

You can generate the documentation locally using sphinx and doxygen. You first have to activate the option `CRACKLET_DOCUMENTATION` with cmake. Then go inside `build/` and run

    make dev_doc

The documentation will be generated in html format and store inside build/doc/html

## Usage - Examples

Exemples of simulations are available in the folder `examples`.

You can create your own CMake project and link it to cRacklet by adding

    find_package(cRacklet REQUIRED)

to your `CMakeLists`. After doing so, you can create an executable with the command

    add_cracklet_simulation(executable source_file
    NU_TOP nu_t
    NU_BOT nu_b)

with `NU_TOP` and `NU_BOT` optional arguments that, if provided, will copy the corresponding kernels to your simulation folder. Note that the convolution kernels describing the elastodynamic response
of the surrounding solids (as function of Poisson's ratio) shall be pre-computed using the Fortran routine
`inverse_serial.f` provided in the folder `pre-computed-kernels/`. The kernels for $`\nu = 0.33`$ and $`\nu = 0.35`$ are available by default.

Basically, the cRacklet engine is composed of seven major types of object:

### SpectralModel:   
This class is the core of cRacklet librairy and contains the methods processing the different steps
required to solve the elastodynamic response of the two semi-infinite half space.

### InterfaceLaw:
This abstract class contains the law describing interface conditions (e.g. cohesive law, frictional interface).
New interface laws can be simply added by creating classes inheriting from InterfaceLaw framework. 

### ContactLaw:
This abstract class contains law describing contact law in case overlapping is prevented at the interface.
New contact laws can be added by creating classes inheriting from ContactLaw.

### Interfacer:
This object helps the user to create and define initial interface properties.
Interfacer is templated by the chosen type of InterfaceLaw.

### DataRegister:
This class is the static centralized memory register which contains most of the simulation data

### SimulationDriver:
High-level class providing a UI to drive simulations in different standard situations.

### DataDumper: 
This class is used to generate various types of output file during simulation.

## PYTHON INTERFACE

In order to use the python interface build for cRacklet, you need pybind11 and python3.

During configuration, activate `CRACKLET_PYTHON_INTERFACE`. The make command will create a python library in the `build/python/` folder. Please add this path to your python path.

You can activate the python example with the option `CRACKLET_EXAMPLES_PYTHON`

An example of the use of the python interface is provided in `build/examples/python/001-modeI/modeI.py`

## Tutorials

Tutorials with the python interface are also available on Binder with a pre-installed version of cRacklet.

[Supershear transition along heterogeneous interface](https://mybinder.org/v2/gl/cracklet%2Ftutorials/72d3a17295b091c74cd57fc5647e9e91dcfc5e51?filepath=supershear%2Fsupershear.ipynb)

[Frictional crack nucleation and propagation for rate and state friction](https://mybinder.org/v2/gl/cracklet%2Ftutorials/2cbe12707784f451c14ffe401e95ab5d6b415b8b?filepath=rate-and-state%2Frate_and_state.ipynb)

The source for the notebook tutorials is available [here](https://gitlab.com/cracklet/tutorials)

## Contributing

Contributions to cRacklet are welcome! Please follow the guidelines below.

### Report an issue

If you have an account on [gitlab](https://gitlab.com), you can [submit an
issue](https://gitlab.com/cracklet/cracklet/-/issues/new). The full list of issues
is available [here](https://gitlab.com/cracklet/cracklet/-/issues).

### Submit a patch / merge-request

Follow [this
guide](https://docs.gitlab.com/ee/user/project/merge_requests/creating_merge_requests.html#new-merge-request-from-a-fork)
to create a merge request on GitLab. Please target the repository's `develop`
branch.

## Citing

To give proper credit to cRacklet and the researchers who have developed it, please cite cRacklet as:
Roch et al., (2022). cRacklet: a spectral boundary integral method library for interfacial rupture simulation. Journal of Open Source Software, 7(69), 3724, https://doi.org/10.21105/joss.03724

## Authors

- Fabian Barras <fabian.barras@epfl.ch>
- Thibault Roch <thibault.roch@epfl.ch>
- Damien Spielmann
- David Kammer
- Guillaume Anciaux
- Nicolas Richart
- Philippe H Geubelle
- Jean-François Molinari

## License

cRacklet is distributed under the terms of the [GNU General Public License
v3.0](https://www.gnu.org/licenses/gpl.html).

## RELATED PUBLICATIONS

The following publications have been made possible with cRacklet:

- [Barras, F., Kammer, D. S., Geubelle, P. H., Molinari, J.-F. (2014) A study of frictional contact in dynamic fracture along bimaterial interfaces. International Journal of Fracture 189(2), 149–162](https://doi.org/10.1007/s10704-014-9967-z)

- [Barras, F., Geubelle, P. H., Molinari, J.-F. (2017) Interplay between Process Zone and Material Heterogeneities for Dynamic Cracks. Physical Review Letters 119(14)](https://doi.org/10.1103/PhysRevLett.119.144101)

- [Barras, F., Carpaij, R., Geubelle, P. H., & Molinari, J.-F. (2018). Supershear bursts in the propagation of a tensile crack in linear elastic material. Physical Review E, 98(6), 063002](https://doi.org/10.1103/PhysRevE.98.063002)

- [Brener, E. A., Aldam, M., Barras, F., Molinari, J.-F., & Bouchbinder, E. (2018). Unstable Slip Pulses and Earthquake Nucleation as a Nonequilibrium First-Order Phase Transition. Physical Review Letters, 121(23), 234302](https://doi.org/10.1103/PhysRevLett.121.234302)

- [Barras, F., Aldam, M., Roch, T., Brener, E. A., Bouchbinder, E., & Molinari, J.-F. (2019). Emergence of Cracklike Behavior of Frictional Rupture: The Origin of Stress Drops. Physical Review X, 9(4), 041043](https://doi.org/10.1103/PhysRevX.9.041043)

- [Barras, F., Aldam, M., Roch, T., Brener, E. A., Bouchbinder, E., & Molinari, J.-F. (2020). The emergence of crack-like behavior of frictional rupture: Edge singularity and energy balance. Earth and Planetary Science Letters, 531, 115978](https://doi.org/10.1016/j.epsl.2019.115978)

- [Fekak, F., Barras, F., Dubois, A., Spielmann, D., Bonamy, D., Geubelle, P. H., & Molinari J. F. (2020). Crack front waves: A 3D dynamic response to a local perturbation of tensile and shear cracks. Journal of the Mechanics and Physics of Solids, 135, 103806](https://doi.org/10.1016/j.jmps.2019.103806)

- [Rezakhani, R., Barras, F., Brun, M., & Molinari, J.-F. (2020). Finite element modeling of dynamic frictional rupture with rate and state friction. Journal of the Mechanics and Physics of Solids, 141, 103967](https://doi.org/10.1016/j.jmps.2020.103967)

- [Brener, E. A., & Bouchbinder, E. (2021). Unconventional singularities and energy balance in frictional rupture. Nature Communications, 12(1), 2585](https://doi.org/10.1038/s41467-021-22806-9)

- [Lebihain, M., Roch, T., Violay, M., & Molinari, J.-F. (2021). Earthquake Nucleation Along Faults With Heterogeneous Weakening Rate. Geophysical Research Letters, 48, e2021GL094901.](https://doi.org/10.1029/2021GL094901)

- [Roch, T., Brener, E. A., Molinari, J.-F., & Bouchbinder, E. (2022). Velocity-driven frictional sliding: Coarsening and steady-state pulses. Journal of the Mechanics and Physics of Solids, 128, 104607](https://doi.org/10.1016/j.jmps.2021.104607)

                                                                                                  
                        ´.-,,,...´´                                                                     
                       .,,,,,,,-----,,,....´´                                                           
                     ´.,,,,,,,,,,,,,,,,,,,,,,,,,....´´´                                                 
                    .,.....,.,,,,,,,,,,,,,,,,,,,,,,,,,,,,.....´´´                                       
                   ´}}[(l>+:°¹,,..,,,,,,,,,,,-,,,,,,,,,,,,,,,,,,¹,......´´                              
                   ´zhaaaaabttahh})i+~¹-,,.,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....´´´                    
                   ´{bbahaabbbttttttmmmtau}?<+;~¹-,..,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,....´´           
                    }abaaaaaabbbtttttttttbbaaahhhau{(!<+:°¹,,.,,,,,,,,,,,,,,,,,,,,,,,,,,,---¹)          
                    ]aaaahaaaabbbttmmmtttbbahhhhaaaaaaaaaabahz]?i+:¹,,......,,,,,,,,,,,,---°OY          
                    )aaaaaa>¹++abbbtttttttbbaaaaaaabaaaaabbabbttmmYVVYbz[?<;°-,...,,,,,,,-°QBY          
                    !haaa{-....-)bbttttttbbaaaaaaabbbbbbbbbbbbttmtmYYYYVVOOOOQOYaz[?i+:°,~MWWY          
                    ¹;>l)-.......~utbttttbbbaaaaaaaaaabbbbbbbbbbtttmmYYVVVOOOOOQQQMRXXXBXRBWXY          
                         ....´.´´,.<bbttttbbbbabbaabbbbbbbbtttttbttttmmYYVVVOOOOOQZZRXWBBXBQXY          
                        ´;,´.´´´´...¹:+![hbbbbbbbbbbbbbbbtttttttttttttmYYVVVOOOOQQZZRWWBBXBmRt          
                        .++¹....´´´...´   ´.,°;>![ubbbtttttmtmmtttttttttmYYYVVOOOOZZMXWBBWWMMa          
                        .+<>~..´´´´....´           ´.-°;i){bmmmmmmmmmYmmmmYYVYVVOOQQZMXWBXXZZi          
                        ,+i<<;,´´´.´,°,..´                  .,°;>?}hYYmtmmmmYYVVVOOOQZMXBXRZz           
                        ¹<lll<+°.....¹¹´´..                         ´.¹~+!]utYYVVVOOQZZMXRMM]           
                        ¹>ii><<i;,....´.¹~...                               ´.-~+?{bVQQZRXQV[           
                        ,>i>~;+<i>°.´....°,...                                      +QQZZRYV[           
                        .><;:~;+ii<;,...°-..´.´.                                    +VOOQQYV}           
                        ´+>::~;++>ll>°..,.....,,..                                  ;VOOQQYV}           
                         :i+;~;+i<l!!i;,.....,.....´                                ;VOOOOYV{           
                         °i>;;++>lll!ll>~....,,.,,.,.                               ;OOOOOYVz           
                         ,i>>++illl!!+ii<+,...-,-..,.´´                             :OOOOOYYz           
                         .<>i<iiilli+i<ili>~,...¹.,..,,.                            ~OOQOQVVu           
                          +<llillli<;~+illli;,..,...,¹,.,´                          ~OQQQQVYh           
                        .i]ll!il!!!l+:++>i!li>°....,°,.¹,,,´                        °QQQQQVYh           
                       ,)tM?!lii!!l!i++<<l!llli:,....,.,,,°,.                       °QQQQQOYa           
                       mVZBu)(?l<i>ii>+iiilll<ll+¹....,,.,,¹-.´                     ¹QQQZQOVa           
                       RQQWR}ut[(?!i+~:;>lllii>ll<;,..,-.,.,;,,.                    ¹QQQZZOVt           
                       iZZBmWXOz]?ll>;;~il<<>+>!li>¹..,,,,,°,,,.´                   -ZZQZZOVt           
                        mZMR##WMb[)ll>;+:;+><>+ll!li:,..,,¹:,,,,,´                  ,ZZZZZQVt           
                         ;QMRR$#$RV}(i!l<<;>+i>>:>+<ll>°...,..,,,,,.                ,ZZZZZQVm           
                          aZMRRB##XQu())!!l<;:;~~:+;;lli;-....,,,,,,,´              ,MMZZZQYY           
                          ¹MMRWBX##BMb}])?!+;;::~~;;:+ill>:,..,,,.¹°,-´             .MMMMMZmm´          
                           (MMRXBW$#$RO{[))!+;;;;;::;+;+ill+¹,..,,-¹,°,´            .MMMMMZmm´          
                           .tMMMXBBB@#XZu[?ll;;;;:~:;+>><+ll<:,.,,,--¹¹.            .MMMMMZmm´          
                           °MRMMRXXBB##BMt[)l<+;+:~::+<li+>!?!+°,,,,,--.            .RMMMMMtt´          
                           ,R#MMRXRXBBBXXRO{(l!!i++;;>ilil+<!??l;-,,---´            .MRRRRMbt´          
                          ,[BBMRRRZMXW#XB#XZh[l)?!!!!?)?li<i!))?l+°,,¹.             ´MRRRRRbh.          
                         ;QQ;+i[VB#$B$$@BB#BMm{])!??))))))!??))?il>::¹              ´MRRRRRbb.          
                         [Omtmtat@YM#W$BBBB##RQu[()]((((((??!!?llli>°               ´MRRRRXtu~          
                        +ttbbbttO#??RBWRXRXW$#WZb{}[[]]((((!ll!ill+¹                ´MRXXRXYmi          
                      ´!tbbbbbbbb#YVW#WRRMRXBB#$MYu{[((()()?!l!??+.                  MXXXXXRh?          
                     .(btbaaahaaam$XRRbOMMXXBBW$#RQh}}](()????li°´                   MXXXXXRX<          
                    .}bthhahhahhhhattuhumRRXBBWXW#WYz{}[(((((l;.                     MXXXXWRt,          
                   -uhhhhahuhuuhuuhuuuuuzhOXBWWRXX$bzz}}}[(([}bu?>~-.                ZWXXXWV{,          
                  ~hhhuuzzuzzuzhzu{{{{zzzz{mRXWWWWRMRZQQMYuuuhhhhhbhhz(i;¹,´         ZWWWXWmb,          
                 +huu{z{{z}{z{}{}{{{z{{u}uzhhZWB$$$WR$$#Mzhuhuhhuuhuhhhzuuzz(l;~,.   ZWWWWWYa°.´        
               ´lzu}uzz}}{}{}{{[}{[}[}z}z}u{zzVOM$#$BMRZa}{{uzzuzzzzzz{zz{{{{{uzuu{)iRWWWWBVOYBO,       
              .(zzz{}{{z}}}[]{[}[}}}{{}}[}}[{{z!)}zaZYbt}z{}{}{}}{{[}}}}}}{}}{}}}{uzzXBBWWBR$@#W+       
             .}{{}[{}{]][][[{][(]}]][[)]}]{}{[h()([}]}[}}{[}][[}}]}{}{[}[}[}[{{[z{z}{RWWWWBXB@RX;       
            ¹}u{{}}[[][[}](([)[((](([[(](]((()h]!]((}][]](]][[}}}[[]}[[[[[][]][[[[[{zRBBWWBX#@BX¹       
           ~uz{[[[[]]((()(](]](?((](}]([?[[?)!(];(])]][[(](](([][[[([[]](([())[[]([[}ZBBBWWZV#$X.       
          ;zz[](}(??]()(])()]))[][]}[[z][(]}[]<:<[](!!)()(l?]()((()ll!)()]?)(]?((((](QWWWBBObVW(        
        ´+{[[](?(?)???!)))()(()][{}uhuzzz}[[(}}i!(((??)!(?l???(??i¹:+;+;<l))?!()?)??!Y#$$BBYt),         
       .¹,°;++!!?)!l?i?l!l)?(!([}ubVVYYmaau{}zau}](??!!??!??!)??+,,..-~~~:;>l)?)?l??)hR#@@@X)°          
       )<:;°,-!!?ll?l!llll!???]}{tVMWBWXMOVttbthzz]())()?!???l!!;,´´´´´.,°°~:;+il!?i?uOXB##[i           
       OZOth}!<>li<i>iill?ll!(]}bYZRB$$###BWRZQtuuz{[}[]]])?)?!()!l~.    ´.,¹°°~:>il)}mZWBa[´           
       YMMMMZZZVbz]l+<>+?<?l?]]htQMXWBB$$$###$RYmVYtbh}z[(]]????l)l?)l:.´   ´.,¹°°~;<[tZRQ}-            
       ~)aZRMRRMMMMMZYa{!<>>l[uaVZRXWWWB$$$##$XQQXWMZOmbaaz}]]])!lllili[(~´    ´.,¹°:;)YR}+             
           .~i}YXXRXRRRRRMZtb[?(}tMWWWWB$$$##$XMMW###$WXZOtbhu}](!??!i>i<[b)-     ´-?>¹°i}.             
                ´¹+)bMXXRRXXXXXRQtu}{hYMB$$###BXXB#@@@@@#$XMQVauz}[)l?<<ii>iht<.  ,zu{[i+;+;,           
                      ,~lzOXXWXWXWWWBXQthubVW$BWWB#@@@@@@@##$XZVbz[)l!i<i<!>>>}Za(X$MYaz{(lil<°.        
                           .¹>[mX$$$$$$$$BWRVYmYVX#@@@@@@@@##BWMYh{]?!i<>i>>+><>)bh{R$bzbaz[)?!?;,´     
                                ´,;?aMBB$$B$$$BBWZVOOMM$@@@##BRQYu{]i<><<i><>+>+>>>i]bb..;{bz[(?l>+°.   
                                      .°<}VB$$$$$$$$$WXQOOZZRXRYtz])lli<<<>i>><+<il<th~    -lahz])l>).  
                                           ´-;)bR$$$$$$#$$BBRQQObz{)il<ii<>>>>+>>>l]a]       ´:[buzh]~  
                                                 .~izOB#$##$$$$$BWMQYh[?i><<<><>>>?Yt.          ,>uYu,  
                                                      ´¹+]tX######$$$BBXMOa}?i+><<tY°              .´   
                                                            .~luQ$#$$$$$$BBBWMZY{zVl                    
                                                                 .¹+]mX##$$$$$$BBXh´                    
                                                                       ,~!uZ$#$#$B,                     
                                                                            .¹+[mM                      
                                                                                                    
---
title: 'cRacklet: a spectral boundary integral method library for interfacial rupture simulation'
tags:
  - boundary integral
  - dynamic rupture
  - elastodynamics
  - friction
  - c++
  - python
authors:
  - name: Thibault Roch
    orcid: 0000-0002-2495-8841
    affiliation: 1
  - name: Fabian Barras
    orcid: 0000-0003-1109-0200
    affiliation: "1, 2"
  - name: Philippe H Geubelle
    orcid: 0000-0002-4670-5474
    affiliation: 3
  - name: Jean-François Molinari
    orcid: 0000-0002-1728-1844
    affiliation: 1
affiliations:
 - name: Civil Engineering Institute, École Polytechnique Fédérale de Lausanne, Switzerland
   index: 1
 - name: The Njord Centre Department of Physics, Department of Geosciences, University of Oslo, Norway
   index: 2
 - name: Department of Aerospace Engineering of the University of Illinois at Urbana-Champaign, United States of America
   index: 3
date: 18 May 2021
bibliography: paper.bib

---

# Summary

The study of dynamically propagating rupture along interfaces is of prime importance in various fields and system sizes, including tribology ($nm$ to $\mu m$), engineering ($mm$ to $m$) and geophysics ($m$ to $km$) [@vanossi_colloquium_2013; @ben-zion_collective_2008; @armstrong-helouvry_survey_1994]. Numerical simulations of these phenomena are computationally costly and challenging, as they usually require the coupling of two different spatio-temporal scales. A fine spatial discretization is needed to represent accurately the singular fields associated with the rupture edges. Besides, the problems of interest usually involve a larger length scale along which rupture will propagate driven by long-range traveling elastic waves. The physical phenomena at play also occur at different timescales, from the slow process of rupture nucleation to the fast propagation of crack front close the elastic wave speeds. Large and finely discretized spatio-temporal domains are required, which are computationally costly. In addition, the behavior of such interfaces can be highly non-linear thus increasing the problem complexity. The use of boundary integral methods reduces the dimensionality of the problem. This enables to focus the computational efforts on the fracture plane and allows for a detailed description of the interfacial failure processes.

# Statement of need

`cRacklet` is a C++ library with a Python interface [@pybind11] initiated as a collaboration between the Computational Solid Mechanics Laboratory at EPFL and the Department of Aerospace Engineering of the University of Illinois at Urbana-Champaign.  `cRacklet` implements a spectral formulation of the elastodynamics boundary integral relations between the displacements and the corresponding traction stress acting at a planar interface between two homogeneous elastic solids [@geubelle_spectral_1995; @breitenfeld_numerical_1998]. The formulation implemented is the *independent* one, which considers the top and bottom solids separately [@breitenfeld_numerical_1998]. The stresses acting at the interface are related to the history of interfacial displacements via a time convolution evaluated in the Fourier domain. The convolutions are efficiently computed within a shared-memory parallel framework using FFTW3/OpenMP. The prescription of an interfacial behavior allows for solving the continuity of tractions and displacements through the interface. Time integration is achieved using an explicit time-stepping scheme. `cRacklet` is aimed at researchers interested in interfacial dynamics, ranging from nucleation problems to dynamic propagation of rupture fronts. While the spectral boundary integral formulation is a well-established method that has been extensively referenced in the literature [@MDSBI; @UGUCA], we believe that `cRacklet` will be a useful addition to the community by gathering in the same framework various kinds of interfacial problems and constitutive laws, and by offering an easy to handle software thanks to its python interface. `cRacklet` is efficient, accessible (C++ or Python), and suited to study a broad class of [problems](https://gitlab.com/cracklet/cracklet/-/tree/master/examples) (fracture and friction). We wish that `cRacklet` will become a link between model developers and users by providing both adaptability and usability.

# Features

1. `cRacklet` is versatile and can be used to study a broad class of problems focused on the behavior of an interface between two semi-infinite solids. The code is particularly suited to study planar dynamic fracture and friction. The interface can be either between two or three-dimensional solids. It can be loaded in any combination of normal traction, in-plane, and out-of-plane shear solicitations. `cRacklet` handles the simulation of interfaces bonded between dissimilar elastic solids. Any stress or material heterogeneity along the fracture plane can be resolved. Several interfacial behaviors are included in the library, such as:

   - Cohesive fracture law [@dugdale_yielding_1960; @barenblatt_mathematical_1962]: the cohesive strength is a linearly decreasing function of the opening gap. This law can be coupled with a friction law to handle surface interactions in the case of post-failure contact between the solids. Two implementations are available, the classical Coulomb friction law and a regularized one [@prakash_frictional_1998].

   - Rate and state dependant friction laws: the frictional resistance is a function of the slip velocity and the history of the interface (the state variable). Several formulations are implemented, including the original ones by @dieterich_modeling_1979 and @ruina_slip_1983. More novel formulations such as rate and state friction with velocity-strengthening behaviors (i.e., N-shaped) are also available, see @barsinai_2014 for example.

2. `cRacklet` is accessible and adaptable. It provides access through both its C++ and Python API to several options to design the various kind of problems mentioned before. `cRacklet` is adaptable due to its object-oriented implementation: it is simple to implement additional behavior for the interface without having to deal with the technical core of the code that handles the computation of the stresses in the Fourier domain. `cRacklet` can also be loaded as an external library to easily interact with other existing computational software. `cRacklet` also has tutorials available on Binder [@binder] which allows for a quick and easy introduction to its functionality.

3. `cRacklet` is efficient: the Fourier transforms and the convolutions are computed within a shared-memory parallel framework using FFTW3/OpenMP. We illustrate in \autoref{fig:scalability} the scaling capability of `cRacklet` and compare it to Amdahl's law. The scaling study shows that approximately $85\%$ to $90\%$ of the program is parallelized: this includes the computation of the Fourier transform of the displacements, the convolution, and the invert transform of the stresses back to the real domain.

![Time required to solve $1\text{e}5$ time steps with $2^{12}$ discretization points, as a function of the number of threads. The code uses `cRacklet` 1.0.0-pre and FFTW 3.3.8, is compiled using GCC [@GCC] and run on the computational facilities of EPFL, here on a node (2 Intel Broadwell processors running at $2.6\,\text{GHz}$ with 14 cores each and 128 GB RAM) of the computing cluster \textit{Fidis}. The dashed grey lines correspond to Amdahl's law for the theoretical speedup, respectively with $90\%$ (upper bound) and $85\%$ (lower bound) of the program parallelized. \label{fig:scalability}](scalability.png){ width=100% }

# Example

The onset of sliding between two rough surfaces in frictional contact is an illustrative example of a multiscale rupture problem. Macroscopic shearing is resisted by the microcontacts, i.e., by the sparse contacting junctions existing between the asperities of the two surfaces.

The successive panels of \autoref{fig:evolution} illustrate the nucleation and propagation of a frictional rupture at the interface between two solids, from the individual failure of the microcontacts in pannel (b) to the propagation of a macroscopic circular rupture in panel (d). The spatially heterogeneous strength used in this example is a representation of the heterogeneous map of contact between two rough surfaces. In \autoref{fig:evolution} (a), the initial configuration of the system is shown. The areas in white are sticking (i.e. no velocity) and correspond to asperities in contact. Colored areas are sliding (blue is for low slip velocity and red for larger ones). The shear load is increased with time in the following panels. The slip velocities increase and previously sticking parts of the interface start sliding (micro-contacts are broken). The inset of \autoref{fig:evolution} (b) is a zoomed view of the interface where rupture starts at the asperity scale. In \autoref{fig:evolution} (d), frictional cracks have expanded over almost the entire interface.

![Snapshot of the slip velocity at the interface between two elastic solids under shear loading. The initial strength is highly heterogeneous. Loading and time have increased between the snapshots, starting from (a) to (d). White areas correspond to sticking conditions (no velocity) while colored ones are sliding. Low velocities are in blue and large ones in red. The code is compiled using [@Intel]. This simulation involve $2^{24}$ points and was run on one node (with two 16-core Intel E5-2683v4 2.1 GHz and 512 GiB RAM) of the computing cluster \textit{Fram} from the Norwegian e-infrastructure for research and education. \label{fig:evolution}](evolution.png){ width=95% }

# Publications

The following publications have been made possible with `cRacklet`: @barras_study_2014, @barras_interplay_2017, @barras_supershear_2018, @brener_unstable_2018, @barras_emergence_2019, @barras_emergence_2020, @fekak_crack_2020, @rezakhani_finite_2020, @brener_unconventional_2021, @lebihain_earthquake_2021, and @roch_velocity_2022.

# Acknowledgments

T.R, F. B., and J-F. M. acknowledge the financial support from the Swiss National Science Foundation (grants #162569 "Contact mechanics of rough surfaces) and from the Rothschild Caesarea Foundation. F.B. acknowledges  support  of  the  Swiss  National  Science  Foundation through the fellowship No. P2ELP2/188034. F.B. acknowledges the Norwegian e-infrastructure for research and education (UNINETT Sigma2) for computing resources through grant NN9814K.

# References
Citing
======

To give proper credit to cRacklet and the researchers who have developed it, please cite cRacklet as:
Roch et al., (2022). cRacklet: a spectral boundary integral method library for interfacial rupture simulation. Journal of Open Source Software, 7(69), 3724, https://doi.org/10.21105/joss.03724
Interface Law
=============

cRacklet has several interfacial behavior implemented which are described in this section. These law are listed in :cpp:enum:`FractureLawType <FractureLawType>`.

Mixed-mode cohesive law
-----------------------

This cohesive law has its strength decreasing linearly with the opening. The opening considered here is the norm of the opening displacement, thus taking into account both normal and shear componenets. The critical displacement :math:`\delta_c` is the critical displacement required to transition from the peak stress :math:`\tau_c` to the residual value :math:`\tau_r` (By default is 0). The shear and normal components can be prescribed independantly. While the norm of the opening :math:`||\delta||` is lower than :math:`\delta_c`, the strength is given by:  

.. math::
   \tau^{str} = \left( \tau_c - \tau_r \right) \left(1- ||\delta||/\delta_c \right)

If :math:`||\delta||` is larger than :math:`\delta_c`, the strength is given by:

.. math::
   \tau^{str} = \tau_r

Friction law for cohesive law
-----------------------------

In addition to the mixed-mode cohesive law, one can defined a friction law to handle cases when the surfaces are in contact with each others. 

Coulomb friction
^^^^^^^^^^^^^^^^

The shear strength is given by the normal stress :math:`\sigma_{yy}` multiplied by a constant friction coefficient :math:`\mu`

.. math::
   \tau^{str} = \sigma_{yy} \mu
   
Regularized Coulomb friction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The classical formulation of Coulomb friction might result in ill-posedness of the friction problem. A simplified regularization is thus implemented, based on `Prakash (1998) <https://asmedigitalcollection.asme.org/tribology/article/120/1/97/439195/Frictional-Response-of-Sliding-Interfaces>`_ and `Rubin and Ampuero (2007) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2006JB004337>`_ . The contact pressure is regularized as:

.. math::
   \frac{\tilde{\sigma}_{yy}}{dt} = -\frac{1}{t^*}\left(\tilde{\sigma}_{yy} - \sigma{yy} \right)

with :math:`t^*` a regularization parameter. The strength is computed as

.. math::
   \tau^{str} = \tilde{\sigma}_{yy} \mu

Rate and state friction
-----------------------

The rate and state framework involve two functionnals: one for the friction coefficient :math:`f(v,\phi)` and one evolution law for the state variable :math:`\dot\phi = g(v,\phi)`. Several variations of the friction law and the evolution law are available in cRacklet. Their exact formulations are summarized here: 

Pure Weakening Formulation
^^^^^^^^^^^^^^^^^^^^^^^^^^

This formulation is based on the one originaly proposed by `Dieterich (1979) <https://doi.org/10.1029/JB084iB05p02161>`_ and `Ruina (1983) <https://doi.org/10.5194/npg-15-1-2008>`_ :

.. math::
   f(v,\phi) = f_0 + a \log \left(v/v_* \right) + b f_0 \log \left(\phi / \phi_* \right)

Standard Formulation
^^^^^^^^^^^^^^^^^^^^

This formulation is a modified version of the original formulation proposed by `Dieterich (1979) <https://doi.org/10.1029/JB084iB05p02161>`_ and `Ruina (1983) <https://doi.org/10.5194/npg-15-1-2008>`_ , with the addition of :math:`+1` in the logarithm.

.. math::
   f(v,\phi) = f_0 + a \log \left(1+v/v_* \right) + b f_0 \log \left(1 + \phi / \phi_* \right)

Regularized Formulation
^^^^^^^^^^^^^^^^^^^^^^^

This formulation is a generic N-shape friction law, as supported by experimental observations `Bar-Sinai & al. (2014) <https://doi.org/10.1002/2013JB010586>`_ . See `Brener & al. (2018) <https://doi.org/10.1103/PhysRevLett.121.234302>`_ for a discussion of the physical sense of each term and the comparison with more conventional rate and state friction law.. The steady-state friction is strenghtening at low and high velocities, and has a velocity-weakening branch at intermediate velocities.

.. math::
   f(v,\phi) = \left( 1 + b \log \left(1+ \frac{\phi}{\phi_*} \right) \right) \left( \frac{f_0}{  \sqrt{\left(  1+v_0^2 / v^2 \right)} } + a \log \left( 1+\frac{v}{v_*} \right) \right)

Regularized Weakening
^^^^^^^^^^^^^^^^^^^^^

This formulation is derived from the N-shape one, with the omission of the :math:`+1` term in the :math:`log` of :math:`\phi` term, resulting in the absence of the velocity-strengthening branch at high velocites.
   
.. math::
   f(v,\phi) = \left( 1 + b \log \left(\frac{\phi}{\phi_*} \right) \right) \left( \frac{f_0}{  \sqrt{\left(  1+v_0^2 / v^2 \right)} } + a \log \left( 1+\frac{v}{v_*} \right) \right)

Aging Law
^^^^^^^^^

The original aging law proposed by `Dieterich (1979) <https://doi.org/10.1029/JB084iB05p02161>`_.

.. math::
   g(\phi) = 1 - \frac{v \phi}{D}

Regularized Aging Law
^^^^^^^^^^^^^^^^^^^^^

A regularization of the slip law, ensuring that for vanishingly small steady-state velocity the value of the state :math:`\phi` saturates to a finite value :math:`D / v_*` after long times instead of diverging.

.. math::
   g(\phi) = 1 - \frac{v \phi}{D} \sqrt{1 + \left(\frac{v_*}{v}\right)^2}

Slip Law
^^^^^^^^

The slip law proposed by `Ruina (1983) <https://doi.org/10.5194/npg-15-1-2008>`_.

.. math::
   g(\phi) = - \frac{v \phi}{D} \log \left( \frac{v \phi}{D} \right)
   
Custom interface law
--------------------

If you want to use a constitutive law that is not implemented in cRacklet, you will need to code it in C++. Please consider creating a merge request for your contribution in GitLab with appropriate documentation if you think that the constitutive behavior you implemented could be used by the community.


If you want to add an evolution law or a friction law that is part of the rate and state framework (depends on a state variable and the sliding velocity), you can define a simple *Struct* as done in the file *rate_and_state_formulations.hh*.

Adding a rate and state evolution law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For an evolution law, you need to define the evolution law :math:`g(\phi)` as the operator () and the evaluation of :math:`g(\phi)` in steady state as the method getSteadyState. In addition, you will need to create an initialization function for the class :cpp:class:`RateAndStateLaw <RateAndStateLaw>` that allocates a shared pointer to an object of your new "struct" (should be similar to :cpp:func:`initStateEvolution <RateAndStateLaw::initStateEvolution>`).

Adding a rate and state friction law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a friction law, you need to define the friction law :math:`f(v,\phi)` as the operator (), the derivative of the friction coefficient with respect to the velocity :math:`\partial f(v,\phi) / \partial v` as the method getTangent, the derivative of the steady-state friction coefficient with respect to the steady velocity :math:`\partial f(v_{ss}) / \partial v_{ss}` as the method getSteadyTangent, and a method that computes the state value satisfying steady state as getStableState.
In addition, an initialization function for the class :cpp:class:`RateAndStateLaw <RateAndStateLaw>` that allocates a shared pointer to an object of your new "struc" should be defined, similarly to :cpp:func:`initStandardFormulation <RateAndStateLaw::initStandardFormulation>`.

To create the interface itself, you should create a new entry in the enum :cpp:enum:`FractureLawType <FractureLawType>` with your new law, and create a templated version of the method :cpp:func:`createUniformInterface <Interfacer::createUniformInterface>` and :cpp:func:`createHeterogeneousInterface <Interfacer::createHeterogeneousInterface>` for the class :cpp:class:`Interfacer`. 

Creating a custom constitutive law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to define a constitutive law that is neither a cohesive law (strength depends only on opening) nor is defined in the rate and state framework (friction depends on velocity and a state variable), you will need to create a new class. This class should inherit from :cpp:class:`InterfaceLaw` and you need to define the following method:

- :cpp:func:`initInterfaceConditions <InterfaceLaw::initInterfaceConditions>`: computes the interface fields for the initial conditions.
- :cpp:func:`updateInterfaceConditions <InterfaceLaw::updateInterfaceConditions>`: computes the interface fields after the displacement has been updated according to the explicit time stepping and the dynamic stresses computed based on the displacement history.
- :cpp:func:`restart <InterfaceLaw::restart>`: save or load the fields required to restart a simulation. This should only concern additional fields used for the constitutive behavior implemented, as the displacement and velocity fields are handled by the mother class.

To solve for the interface fields, you can take inspiration from the classes :cpp:class:`CohesiveLaw` or :cpp:class:`RateAndStateLaw`.
With the class :cpp:class:`CohesiveLaw`, the system is well defined and the velocity is directly expressed in terms of the loading, the dynamic stresses and the material parameters. The normal and shear velocities are computed independently. The strength of the interface has to be evaluated against the stress to assess if interface opening is occuring.
The :cpp:class:`RateAndStateLaw` however uses an iterative Newton-Raphson procedure to find the value of the velocity and the state variable to satisfy that the interfacial stress is always equal to the strength.

To create the interface itself, you should create a new entry in the enum :cpp:enum:`FractureLawType <FractureLawType>` with your new law, and create a templated version of the method :cpp:func:`createUniformInterface <Interfacer::createUniformInterface>` and :cpp:func:`createHeterogeneousInterface <Interfacer::createHeterogeneousInterface>` for the class :cpp:class:`Interfacer`.
Examples
========

The directory ``examples/`` in cRacklet' root repository contains example scripts
dedicated to various aspects of cRacklet:

Bi-material interface
---------------------

Dynamic debonding at a planar interface btw Aluminium (mtl) and Homalite (poly). This example is a 1D interface loaded in mixed mode (Mode I and Mode II). 

Super-shear transition and constant driving speed
-------------------------------------------------

This folder contains two script, dealing with a 1D interface between two similar bodies. The interface follows a linear coupled cohesive law. The properties of the interface are spatially heterogeneous and consist of successive weak and strong stripes.

andrews_transition.cc
^^^^^^^^^^^^^^^^^^^^^

Study super-shear transition following Andrews(1976) formalism (accelerating mode II shear crack that transition to super-shear velocity via the nucleation of a daughter crack ahead of the primary crack)

cst_speed_fracture.cc
^^^^^^^^^^^^^^^^^^^^^

Dynamic fracture at a constant crack speed in presence of heterogeneities. The loading is updated dynamically such that the rupture propagates at a given velocity (sub-Rayleigh)

3D Interface in the presence of a strong asperity
-------------------------------------------------

3d interface in the presence of one asperity.

Rate and state friction laws
----------------------------

Rupture nucleation on an interface following a generic rate and state friction law. The interface is initially sliding at a steady-state velocity and nucleation occurs via a slight perturbation of the velocity field.
API Reference
=============

Python API
----------

.. autoclass:: py11_cracklet.SpectralModel
   :members:
   :show-inheritance:

   .. automethod:: __init__

.. autoclass:: py11_cracklet.SimulationDriver
   :members:
   :show-inheritance:

   .. automethod:: __init__
      
cracklet Interface Laws
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: py11_cracklet.CohesiveLaw
   :members:
   :show-inheritance:

   .. automethod:: __init__
      
.. autoclass:: py11_cracklet.InterfacerLinearCoupledCohesive
   :members:
   :show-inheritance:
      
.. autoclass:: py11_cracklet.RateAndStateLaw
   :members:
   :show-inheritance:
      
cracklet Interfacer
^^^^^^^^^^^^^^^^^^^
      
.. autoclass:: py11_cracklet.InterfacerRateAndState
   :members:
   :show-inheritance:

.. autoclass:: py11_cracklet.InterfacerRegularizedRateAndState
   :members:
   :show-inheritance:
      
.. autoclass:: py11_cracklet.InterfacerWeakeningRateAndState
   :members:
   :show-inheritance:
      
cracklet Contact Laws
^^^^^^^^^^^^^^^^^^^^^
      
.. autoclass:: py11_cracklet.CoulombLaw
   :members:
   :show-inheritance:

.. autoclass:: py11_cracklet.RegularizedCoulombLaw
   :members:
   :show-inheritance:

cracklet Dumper and Data Management
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: py11_cracklet.DataRegister
   :members:
   :show-inheritance:

.. autoclass:: py11_cracklet.DataDumper
   :members:
   :show-inheritance:
Quickstart
==========

Here is a quick introduction to get started with cRacklet.

Installation from source
------------------------

The following dependencies are required for cRacklet:

 - a **C++ compiler**
 - **CMake**
 - **FFTW3**
 - **GSL** library

If you have Debian/Ubuntu based system you can install them with::

  apt install cmake gcc libfftw3-dev libgsl-dev

If you have a MacOS::

  brew install cmake gcc fftw gsl

f you have a windows system you can use windows sub system Linux and check apt instruction.
  
Optional dependencies are:

 - **pybind11** (for python binding, automatically installed)
 - **python3** (for python binding)
 - a compiler that supports **OpenMP** (for multi-threaded parallel computing)
 - **pytest** (for tests)
 - **Doxygen** and **Sphinx** (for documentation)
 - a **Fortran** compiler to generate new kernels
  
First clone the git repository::

  git clone https://gitlab.com/cracklet/cracklet
  
You can then compile cRacklet using cMake::

  cd cracklet
  mkdir build
  cd build
  ccmake .. or cmake ..
  [ Set the options that you need ]
  make
  
Building the docs
-----------------

To build the documentation locally, activate the documentation option with **CMake** (CRACKLET_DOCUMENTATION). Then from the build directory::
  
  make dev_doc

If you want to generate the documentation for the python interface, you first need to generate the python module. See :ref:`python-interface`.
  
Running the test
----------------

You need to activate the test options with **CMake** (CRACKLET_TESTS). You can then run the test with::

  make test
Authors
=======

- Fabian Barras <fabian.barras@mn.uio.no>
- Thibault Roch <thibault.roch@epfl.ch>
- Damien Spielmann
- David Kammer
- Guillaume Anciaux
- Nicolas Richart
- Philippe H Geubelle
- Jean-François Molinari
Publications
============

The following publications have been made possible with cRacklet:

- `Barras, F., Kammer, D. S., Geubelle, P. H., Molinari, J.-F. (2014) A study of frictional contact in dynamic fracture along bimaterial interfaces. International Journal of Fracture 189(2), 149–162 <https://doi.org/10.1007/s10704-014-9967-z>`_

- `Barras, F., Geubelle, P. H., Molinari, J.-F. (2017) Interplay between Process Zone and Material Heterogeneities for Dynamic Cracks. Physical Review Letters 119(14) <https://doi.org/10.1103/PhysRevLett.119.144101>`_

- `Barras, F., Carpaij, R., Geubelle, P. H., & Molinari, J.-F. (2018). Supershear bursts in the propagation of a tensile crack in linear elastic material. Physical Review E, 98(6), 063002. <https://doi.org/10.1103/PhysRevE.98.063002>`_

- `Brener, E. A., Aldam, M., Barras, F., Molinari, J.-F., & Bouchbinder, E. (2018). Unstable Slip Pulses and Earthquake Nucleation as a Nonequilibrium First-Order Phase Transition. Physical Review Letters, 121(23), 234302. <https://doi.org/10.1103/PhysRevLett.121.234302>`_

- `Barras, F., Aldam, M., Roch, T., Brener, E. A., Bouchbinder, E., & Molinari, J.-F. (2019). Emergence of Cracklike Behavior of Frictional Rupture: The Origin of Stress Drops. Physical Review X, 9(4), 041043. <https://doi.org/10.1103/PhysRevX.9.041043>`_

- `Barras, F., Aldam, M., Roch, T., Brener, E. A., Bouchbinder, E., & Molinari, J.-F. (2020). The emergence of crack-like behavior of frictional rupture: Edge singularity and energy balance. Earth and Planetary Science Letters, 531, 115978. <https://doi.org/10.1016/j.epsl.2019.115978>`_

- `Fekak, F., Barras, F., Dubois, A., Spielmann, D., Bonamy, D., Geubelle, P. H., & Molinari J. F. (2020). Crack front waves: A 3D dynamic response to a local perturbation of tensile and shear cracks. Journal of the Mechanics and Physics of Solids, 135, 103806. <https://doi.org/10.1016/j.jmps.2019.103806>`_

- `Rezakhani, R., Barras, F., Brun, M., & Molinari, J.-F. (2020). Finite element modeling of dynamic frictional rupture with rate and state friction. Journal of the Mechanics and Physics of Solids, 141, 103967. <https://doi.org/10.1016/j.jmps.2020.103967>`_

- `Brener, E. A., & Bouchbinder, E. (2021). Unconventional singularities and energy balance in frictional rupture. Nature Communications, 12(1), 2585. <https://doi.org/10.1038/s41467-021-22806-9>`_

- `Lebihain, M., Roch, T., Violay, M., & Molinari, J.-F. (2021). Earthquake Nucleation Along Fault With Heterogeneous Weakening Rate. Geophysical Research Letters, 48 e2021GL094901. <https://doi.org/10.1029/2021GL094901>`_

- `Roch, T., Brener, E. A., Molinari, J.-F., & Bouchbinder, E. (2022). Velocity-driven frictional sliding: Coarsening and steady-state pulse trains. Journal of the Mechanics and Physics of Solids, 158, 104607. <https://doi.org/10.1016/j.jmps.2021.104607>`_
Tutorials
=========

Tutorials for cRacklet in the form of Jupyter notebook are available in the following git repository: https://gitlab.com/cracklet/tutorials . You can either clone the repository and then follow the instruction in the **README.md** to install the tutorials on your machine, or you can directly do the tutorials directly on `Binder <https://doi.org/10.25080/Majora-4af1f417-011>`_ .

The following tutorials are available:

`Supershear transition along heterogeneous interface <https://mybinder.org/v2/gl/cracklet%2Ftutorials/72d3a17295b091c74cd57fc5647e9e91dcfc5e51?filepath=supershear%2Fsupershear.ipynb>`_ 

`Frictional crack nucleation and propagation for rate and state friction <https://mybinder.org/v2/gl/cracklet%2Ftutorials/2cbe12707784f451c14ffe401e95ab5d6b415b8b?filepath=rate-and-state%2Frate_and_state.ipynb>`_ 
.. cRacklet documentation master file, created by
   sphinx-quickstart on Wed Feb  3 18:37:50 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

cRacklet - A boundary integral library for interfacial rupture simulations
##########################################################################

.. image:: https://gitlab.com/cracklet/cracklet/badges/master/pipeline.svg
   :target: https://gitlab.com/cracklet/cracklet/-/commits/master

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gl/cracklet%2Ftutorials/master?filepath=supershear%2Fsupershear.ipynb

.. image:: https://joss.theoj.org/papers/10.21105/joss.03724/status.svg
   :target: https://doi.org/10.21105/joss.03724
	    
*cRacklet* is a C++ boundary integral library based on a spectral formulation of the dynamic wave equations in two semi-infinite linearly elastic solids `Geubelle and Rice (1995) <http://www.sciencedirect.com/science/article/pii/002250969500043I>`_  `Breitenfeld and Geubelle (1998) <https://link.springer.com/article/10.1023/A:1007535703095>`_ . Its implementation is specially tailored for the modeling of dynamic crack/rupture propagation along a planar interface bonding the two half-spaces. The main benefit of this spectral method is a numerical discretization limited to the mid-plane, thereby providing a very fine description of the dynamic rupture processes, unattainable with finite-element or finite-difference schemes. For more details about the method and its applications, refer to the publications section.

.. toctree::
   :maxdepth: 2
   :caption: Table of Contents:

   ./quickstart
   ./simulation
   ./interface_law
   ./python_interface
   ./examples
   ./tutorials
   ./api_reference
   ./publications
   ./citing
   ./authors

Seeking help - Bug reports
--------------------------

You can ask your questions or report any bugs you find here https://gitlab.com/cracklet/cracklet/-/issues

Code contribution
-----------------

Any contribution to cRacklet, bugfixes or new features, is welcome with a merge request on GitLab.
   
Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
.. _python-interface:

Python interface
================

Compile the library from sources
--------------------------------

In order to use the python interface build for cRacklet, you need pybind11 and python3.

During configuration, activate **CRACKLET_PYTHON_INTERFACE**. Running::

  make

in **build** folder will compile the python library called cracklet in **build/python/** folder. Please add this path to your python path. Otherwise you can do ::

  make install
  
You can activate the python examples with the option **CRACKLET_EXAMPLES_PYTHON**.

Install via pip
---------------

Alternatively, one can install cracklet using pip. The last version of the package can be installed with::

  pip install cracklet --extra-index-url https://gitlab.com/api/v4/projects/25914378/packages/pypi/simple

Use the python interface
------------------------

Once the library is compiled and its path added to your python path, you can import cRacklet with::

  import cracklet

One can run a complete simulation using the python interface. The required steps to run a simulation are identical to the **c++** counterpart. First a :py:class:`SpectralModel <py11_cracklet.SpectralModel>` needs to be created::

  model = cracklet.SpectralModel(nb_elements,nb_time_steps,dom_size,nu,E,cs,tcut,simulation_summary,output_dir)

A :py:class:`SimulationDriver <py11_cracklet.SimulationDriver>` can then be associated to the model::

  sim_driver = cracklet.SimulationDriver(model)
  
The parameters can be registered using :py:func:`RegisterParameterReal <py11_cracklet.DataRegister::registerParameterReal>`, and the :py:class:`Interfacer <py11_cracklet.InterfacerLinearCoupledCohesive>` can then be used to create the interface (here, an example with an homogeneous interface following a linear coupled cohesive law, with a broken area)::

  cracklet.DataRegister.registerParameterReal("critical_normal_opening",crit_n_open)
  cracklet.DataRegister.registerParameterReal("critical_shear_opening",crit_s_open)
  cracklet.DataRegister.registerParameterReal("max_normal_strength",max_n_str)
  cracklet.DataRegister.registerParameterReal("max_shear_strength",max_s_str)
  interfacer = cracklet.InterfacerLinearCoupledCohesive(model)    
  interfacer.createUniformInterface()
  interfacer.createThroughCrack(crack_start,crack_end)

The loading can then be initiated with::

   sim_driver.initConstantLoading(load,psi,phi)
  
The :py:class:`DataDumper <py11_cracklet.DataDumper>` is used to dump the interfacial fields that are available in :py:obj:`DataFields <py11_cracklet.DataFields>` ::
   
  dumper = DataDumper(model)
  dumper.initDumper(output_file,cracklet.DataFields._top_displacements)

One can then solve a step with ::
  
  sim_driver.solveStep()

and dump the results with::

  dumper.dumpAll()
Simulation
==========

Find cRacklet and add a simulation
----------------------------------

You can create your own **cMake** project and link it to **cRacklet** by adding::

  find_package(cRacklet REQUIRED)

to your **CMakeLists.txt**. After doing so you can create an executable with::

  add_cracklet_simulation(executable source
                          NU_TOP nu_t
                          NU_BOT nu_b)

with NU_TOP and NU_BOT optional arguments that, if provided, will copy the corresponding kernels to your simulation folder. nu_t and nu_b are the value of the Poisson modulus associated with respectively the top and the bottom material.

Only specific kernels are already generated with the distributed version of cRacklet. You can generate new kernels using the Fortran routine **invert_serial.f** found in the **pre-computed-kernels/** folder in the root of **cRacklet** repository. You will need to change the inputs used by **invert_serial.f** in **invert_serial.in**.

Simulation structure
--------------------

Create the model
^^^^^^^^^^^^^^^^

In order to run a **cRacklet** simulation, you will need first to create a :cpp:class:`SpectralModel <SpectralModel>` object. This class is the core of cRacklet librairy and contains the methods processing the different steps required to solve the elastodynamic response of the two semi-infinite half space as desribed in [2].

The model for a 2D simulation with similar material on top and bottom is the following::

  #include "spectral_model.hh"

  int main(){
  
  SpectralModel model(nb_elements,nb_time_steps,dom_size,nu,E,cs,tcut,simulation_summary,output_dir);

  // your code ....

  }
  
with ``nb_elements`` the number of discretization points, ``nb_time_steps`` the number of time steps, ``dom_size`` the real length of the interface, ``nu``, ``E``, ``cs`` the characteristics of the elastic bulk (respectively the Poisson ratio, the Young Modulus and the shear wave speed velocity), ``t_cut`` the cut in the kernels, ``simulation_summary`` a string containing a description of your simmulation and ``output_dir`` the folder in which the output will be written.

It is then required to init the model and setting the loading case::

  model.initModel();
  model.setLoadingCase(load,psi,phi);

with ``psi`` and ``phi`` the orientation of the load and ``load`` its absolute value.

Create the interface
^^^^^^^^^^^^^^^^^^^^

The second step is to create and define the initial interface properties, with the help of the :cpp:class:`Interfacer <Interfacer>` object which is templated by one of the possible :cpp:enum:`FractureLawType <FractureLawType>`. The corresponding headers have to be included (here an example for the linear coupled cohesive law)::
  
  #include "interfacer.hh"
  #include "cohesive_law.hh"

  int main(){

  // Model creation...

  Interfacer<_linear_coupled_cohesive> interfacer(model);
  
  }
  
Each interface law has build-in methods method allowing for the creation of specific interface layout, see the API documentation for more details. The simplest case is to create a spatially homogeneous interface. It requires some parameters to be registered in the model first: this can either be done manually with :cpp:func:`registerParameter <DataRegister::registerParameter>` of the :cpp:class:`DataRegister <DataRegister>` class or using the built-in reader :cpp:func:`readInputFile <DataRegister::readInputFile>` that parse an input file and register every parameters. As an example, to create an uniform interface with a cohesive law behavior, one can write::

  DataRegister::registerParameter("critical_normal_opening",crit_n_open);
  DataRegister::registerParameter("critical_shear_opening",crit_s_open);
  DataRegister::registerParameter("max_normal_strength",max_n_str);
  DataRegister::registerParameter("max_shear_strength",max_s_str);
  interfacer.createUniformInterface();

with ``critical_normal_opening``, ``critical_shear_opening``, ``max_normal_strength`` and ``max_shear_strength`` being parameters describing the interface response.

Additionaly, one can add a pre-inserted crack with the following call, setting the resistance to 0 between ``start_crack`` and ``end_crack``  ::

  interfacer.createThroughCrack(start_crack,end_crack);

Then, the loads can be properly applied and the interface fields (stresses, velocities, displacements) can be initiated::

  model.updateLoads();
  model.initInterfaceFields();

Setting up the dumper
^^^^^^^^^^^^^^^^^^^^^

A :cpp:class:`DataDumper <DataDumper>` object needs to be created. One can then select the fields to dump with the various method associted to the DataDumper class. Here we give an example where the dumper is configure to output the displacement field of the top elastic body in the ``output`` file. The fields that can be dumped are found in :cpp:enum:`DataFields <DataFields>` ::

  DataDumper dumper(model);
  dumper.initDumper("output", _top_displacements);

Solving steps
^^^^^^^^^^^^^

Solving the equilibrium of stresses at the interface requires to execute several functions from the object::

  model.updateDisplacements();
  model.fftOnDisplacements();
  model.computeStress();
  model.computeInterfaceFields();
  model.increaseTimeStep();

The output can be generated by calling::

  dumper.dumpAll()

Simulation Driver
^^^^^^^^^^^^^^^^^
  
Alternatively, one can ressort to :cpp:class:`SimulationDriver <SimulationDriver>` object to handle all the necessary operations required to solve a step. The constructor of the :cpp:class:`SimulationDriver <SimulationDriver>` object takes a :cpp:class:`SpectralModel <SpectralModel>` as an argument. The construction of the simulation driver include the model initialisation. The loading can be handle directly from the :cpp:class:`SimulationDriver <SimulationDriver>`. The whole solve step routine can then be reduce to a single call, as shown in the example below for a homogenous constant loading ::

  SimulationDriver sim_driver(model);  
  sim_driver.initConstantLoading(load,psi,phi);  
  sim_driver.solveStep();
