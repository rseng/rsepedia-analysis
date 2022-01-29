---
title: 'COMPAS: A rapid binary population synthesis suite'
tags:
  - Python
  - C++
  - astronomy
  - gravitational waves
  - binary evolution
authors:
  - name: Team COMPAS
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Jeff Riley 
    affiliation: "2,3"
  - name: Poojan Agrawal
    affiliation: "4,3"
  - name: Jim W. Barrett
    affiliation: 5
  - name: Kristian N. K. Boyett
    affiliation: 6
  - name: Floor S. Broekgaarden
    affiliation: 7
  - name: Debatri Chattopadhyay
    affiliation: "8,4,3"
  - name: Sebastian M. Gaebel
    affiliation: 9
  - name: Fabian Gittins
    affiliation: 10
  - name: Ryosuke Hirai
    affiliation: "2,3"
  - name: George Howitt
    affiliation: 11
  - name: Stephen Justham
    affiliation: "12,13,14"
  - name: Lokesh Khandelwal
    affiliation: 12
  - name: Floris Kummer
    affiliation: 12
  - name: Mike Y. M. Lau
    affiliation: "2,3"
  - name: Ilya Mandel
    affiliation: "2,3,5"
  - name: Selma E. de Mink
    affiliation: "14,12,7"
  - name: Coenraad Neijssel
    affiliation: "5,3"
  - name: Tim Riley
    affiliation: "2,3"
  - name: Lieke van Son
    affiliation: "7,12,14"
  - name: Simon Stevenson
    affiliation: "4,3"
  - name: Alejandro Vigna-Gómez
    affiliation: "15,16"
  - name: Serena Vinciguerra
    affiliation: 12
  - name: Tom Wagg
    affiliation: "7,14"
  - name: Reinhold Willcox
    affiliation: "2,3"
affiliations:
 - name: The public COMPAS code is a product of work by the entire COMPAS collaboration over many years; we therefore kindly request that, in recognition of this team effort, the paper is cited as Team COMPAS - J. Riley et al.
   index: 1
 - name: School of Physics and Astronomy, Monash University, Clayton, Victoria 3800, Australia
   index: 2
 - name: OzGrav, Australian Research Council Centre of Excellence for Gravitational Wave Discovery, Australia
   index: 3
 - name: Centre for Astrophysics and Supercomputing, Swinburne University of Technology, Hawthorn, VIC 3122, Australia
   index: 4
 - name: Institute of Gravitational Wave Astronomy and School of Physics and Astronomy, University of Birmingham, Birmingham, B15 2TT
   index: 5
 - name: Department of Physics, University of Oxford, Denys Wilkinson Building, Keble Road, Oxford OX1 3RH, UK
   index: 6
 - name: Center for Astrophysics |Harvard & Smithsonian, 60 Garden St., Cambridge, MA 02138, USA
   index: 7
 - name: School of Physics and Astronomy, Cardiff University, Cardiff, CF24 3AA, United Kingdom
   index: 8
 - name: Max Planck Institute for Gravitational Physics (Albert Einstein Institute), Callinstrasse 38, D-30167 Hannover, Germany
   index: 9
 - name: Mathematical Sciences and STAG Research Centre, University of Southampton, Southampton SO17 1BJ, UK
   index: 10
 - name: School of Physics, University of Melbourne, Parkville, Victoria, 3010, Australia
   index: 11
 - name: Anton Pannekoek Institute of Astronomy and GRAPPA, Science Park 904, University of Amsterdam, 1098XH Amsterdam, The Netherlands
   index: 12
 - name: School of Astronomy & Space Science, University of the Chinese Academy of Sciences, Beijing 100012, China
   index: 13
 - name: Max Planck Institute for Astrophysics, Karl-Schwarzschild-Str. 1, 85748 Garching, Germany
   index: 14
 - name: DARK, Niels Bohr Institute, University of Copenhagen, Jagtvej 128, 2200, Copenhagen, Denmark
   index: 15
 - name: Niels Bohr International Academy, The Niels Bohr Institute, Blegdamsvej 17, 2100 Copenhagen, Denmark
   index: 16
date: xx Month 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/1538-4365/ac416c
aas-journal: Astrophysical Journal Supplements 
---

# Summary

Most massive stars---those with initial masses greater than 8 $M_\odot$---are born with another massive star as a companion [@Sana:2012Sci;@Moe:2017ApJS]. Massive binary stars are responsible for producing many exotic astrophysical phenomena, such as the observed diversity of supernovae, binary pulsars, X-ray binaries and merging compact objects. The latter are now regularly observed by the ground-based gravitational wave observatories Advanced LIGO and Virgo [@abbott2016observation;@GWTC3].  Population models of massive binary evolution make it possible to interpret existing observations and to make predictions for future observing campaigns.  

# Statement of need

Binary population synthesis generates population models of isolated stellar binaries under a set of parametrized assumptions.  These models permit comparisons against observational data sets, such as X-ray binaries of gravitational-wave mergers.   

In particular, rapid binary population synthesis is needed in order to efficiently explore a broad parameter space of uncertain assumptions about the physics of stellar and binary evolution, including supernova remnant masses and natal kicks, mass transfer efficiency and stability, and the outcome of common-envelope events.  

A range of binary population synthesis codes have been developed over the last three decades.  These include the Scenario Machine [@Scenario], IBiS [@IBiS], SeBa [@SeBa],  BSE [@Hurley:2002rf], StarTrack [@Belczynski:2008],  binary$\_$c [@BinaryC],  MOBSE [@2018MNRAS.474.2959G]  and COSMIC [@2019arXiv191100903B].  These codes range from private to semi-public to fully public, and differ in the range of available tools, computational complexity, and speed of execution.

[COMPAS](https://compas.science) is a rapid binary population synthesis suite. It parametrizes complex astrophysical processes with prescriptions calibrated to detailed models.  COMPAS is designed to allow for flexible modifications as evolutionary models improve.  All code is fully public and, including pre-processing and post-processing tools.  COMPAS is computationally efficient, with a focus on the statistical analysis of large populations, particularly but not exclusively in the context of gravitational-wave astronomy.  


# Details

The core engine of COMPAS---responsible for calculating the evolution of single [@Hurley:2000pk] and binary [@Hurley:2002rf] stars---is written in object oriented C++ for speed and flexibility. COMPAS is able to simulate the evolution of a typical binary over 10 Gyr in approximately 10 milliseconds.

A detailed description of the implementation of the COMPAS suite can be found in @COMPAS:2021methodsPaper.

In addition to the core stellar and binary evolution engine, we provide Python scripts for both pre- and post-processing COMPAS outputs. Post-processing can account for integrating populations formed throughout cosmic history [@2019MNRAS.490.3740N] and methods to account for gravitational-wave selection effects [@Barrett:2017fcw]. A set of examples is also provided.

COMPAS is *embarrassingly* parallel and can be trivially run on high performance computers and distributed on cloud computing.

COMPAS was initially designed to focus on studies of merging binaries containing neutron stars and black holes that are being observed through gravitational waves [@Stevenson2017FormationEvolution;@2018MNRAS.481.4009V]. 
In recent years, the scope of systems investigated with COMPAS has expanded to incorporate, e.g., Be X-ray binaries [@Vinciguerra:2020] and luminous red novae [@Howitt:2020] (see @COMPAS:2021methodsPaper or [the COMPAS collaboration website](https://compas.science) for a summary of COMPAS publications to date.)

COMPAS development happens on [Github](https://github.com/TeamCOMPAS/COMPAS). We maintain a [Zenodo community](https://zenodo.org/communities/compas/) where data from many publications using COMPAS is publicly available. 


# Acknowledgements

Multiple authors are supported by the Australian Research Council Centre of Excellence for Gravitational Wave Discovery (OzGrav), through project number CE170100004. Multiple authors were funded in part by the National Science Foundation under Grant No. (NSF grant number 2009131), the Netherlands Organization for Scientific Research (NWO) as part of the Vidi research program BinWaves with project number 639.042.728 and by the European Union’s Horizon 2020 research and innovation program from the European Research Council (ERC, Grant agreement No. 715063).  FSB is supported in part by the Prins Bernard Cultuurfonds studiebeurs. IM is a recipient of an Australian Research Council Future Fellowship (FT190100574).  AVG acknowledges funding support by the Danish National Research Foundation (DNRF132)


# References
[//]: ## (grip -b README.md)

![COMPASlogo](docs/media/COMPASlogo.png)

# Compact Object Mergers: Population Astrophysics & Statistics

[//]: ## (Outline features)
COMPAS is a publicly available rapid binary population synthesis code (https://compas.science/) that is designed so that evolution prescriptions and model parameters are easily 
adjustable.  COMPAS draws properties for a binary star system from a set of initial distributions, and evolves it from zero-age main sequence to the end of its life as two compact 
remnants.  It has been used for inference from observations of gravitational-wave mergers, Galactic neutron stars, X-ray binaries, and luminous red novae.

## Documentation
https://compas.science/docs

## Contact
Please email your queries to compas-user@googlegroups.com. You are also welcome to join the [COMPAS User Google Group](https://groups.google.com/forum/#!members/compas-user) to engage in discussions with COMPAS users and developers.

## Acknowledgements
If you use this code or parts of this code for results presented in a scientific publication, we would greatly appreciate if you send us your paper reference and make your input settings and output data publicly available by uploading it to the [COMPAS Zenodo community](https://zenodo.org/communities/compas/). Please also kindly include citations to our COMPAS methods paper https://ui.adsabs.harvard.edu/abs/2021arXiv210910352T/abstract. As the public COMPAS code is a product of work by the entire COMPAS collaboration over many years, we kindly request that, in recognition of this team effort, the paper is cited as “Team COMPAS: J. Riley et al.”. An example bibtex code is:


@ARTICLE{COMPAS:2021methodsPaper,
       author = {{\noopsort{Team COMPAS}}{Team COMPAS: J. Riley} and  {Riley}, Jeff and {Agrawal}, Poojan and {Barrett}, Jim W. and {Boyett}, Kristan N.~K. and {Broekgaarden}, Floor S. and {Chattopadhyay}, Debatri and {Gaebel}, Sebastian M. and {Gittins}, Fabian and {Hirai}, Ryosuke and {Howitt}, George and {Justham}, Stephen and {Khandelwal}, Lokesh and {Kummer}, Floris and {Lau}, Mike Y.~M. and {Mandel}, Ilya and {de Mink}, Selma E. and {Neijssel}, Coenraad and {Riley}, Tim and {van Son}, Lieke and {Stevenson}, Simon and {Vigna-Gomez}, Alejandro and {Vinciguerra}, Serena and {Wagg}, Tom and {Willcox}, Reinhold},
        title = "{Rapid stellar and binary population synthesis with COMPAS}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Solar and Stellar Astrophysics},
         year = 2021,
        month = sep,
archivePrefix = {arXiv},
       eprint = {2109.10352},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021arXiv210910352T},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

Note that the preferred acknowledgement relies on \noopsort; to make it work, you'll have to include the following line at the start of your bibtex file:
@PREAMBLE{ {\providecommand{\noopsort}[1]{}} }



In addition, we suggest to kindly include the two following papers:

1. Stevenson S., Vigna-Gómez A., Mandel I., Barrett J. W., Neijssel C. J., Perkins D., de Mink S. E., 2017, [Nature Communications, 8, 14906](https://ui.adsabs.harvard.edu/abs/2017NatCo...814906S/abstract)
2. Vigna-Gómez A., Neijssel C. J., Stevenson S., Barrett J. W., Belczynski K., Justham S., de Mink S., M&uuml;ller B., Podsiadlowski Ph., Renzo M., Szécsi D., Mandel I., 2018, [MNRAS, 481, 4009](https://ui.adsabs.harvard.edu/abs/2018MNRAS.481.4009V/abstract)

We also greatly appreciate an acknowledgement of the form: 

>_Simulations in this paper made use of the COMPAS rapid binary population synthesis code (version X.X.X), which is freely available at http://github.com/TeamCOMPAS/COMPAS_.

Furthermore,

  * If you use COMPAS's importance sampling algorithm STROOPWAFEL, please cite 

     Broekgaarden F. S., Justham S., de Mink S. E., Gair J., Mandel I., Stevenson S., Barrett J. W., Vigna-Gómez A., Neijssel C. J., 2019, [MNRAS, 490, 5228](https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.5228B/abstract)

  * If using the COMPAS model of gravitational wave selection effects, please cite

     Barrett J. W., Gaebel S. M., Neijssel C. J., Vigna-Gómez A., Stevenson S., Berry C. P. L., Farr W. M., Mandel I., 2018, [MNRAS, 477, 4685](https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.4685B/abstract)

  * If using COMPAS's integration over cosmic star formation history, please cite 

     Neijssel C. J., Vigna-Gómez A., Stevenson S., Barrett J. W., Gaebel S. M., Broekgaarden F. S., de Mink S. E., Szécsi D., Vinciguerra S., Mandel I., 2019, [MNRAS, 490, 3740](https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3740N/abstract)

  * If using the COMPAS model of (pulsational) pair instability supernova, please cite 

     Stevenson S., Sampson M., Powell J., Vigna-Gómez A., Neijssel C. J., Szécsi D., Mandel I., 2019, [ApJ, 882, 121](https://ui.adsabs.harvard.edu/abs/2019ApJ...882..121S/abstract)
     
  * If evolving pulsar spins and magnetic fields with COMPAS, please cite
  
     Chattopadhyay D., Stevenson S., Hurley J. R., Rossi L. J., Flynn C., 2020,  [MNRAS](https://ui.adsabs.harvard.edu/abs/2020MNRAS.tmp..697C/abstract)

## License
[MIT](https://choosealicense.com/licenses/mit/)



### Highlighted papers that have made use of COMPAS are listed at https://compas.science/science.html ; see https://ui.adsabs.harvard.edu/public-libraries/gzRk1qpbRUy4cP2GydR36Q for a full ADS library


# COMPAS Code of Conduct

Dear colleagues:

I want to take this opportunity to remind all of us (myself included) 
how important it is to maintain the highest standards of respectful and tolerant 
behavior in our professional interactions.  

We are fortunate that Team COMPAS is a diverse group, and we want to benefit 
from that diversity and make sure that everyone in the team feels comfortable and welcome. 

There are numerous codes of conduct describing what is considered appropriate 
professional behaviour;  since many of us are based in Australia and connected to OzGrav, 
let me point to these two in particular:

ASA code of conduct: <https://asa.astronomy.org.au/wp-content/uploads/2020/09/ASA_Code_of_Conduct.pdf>

OzGrav code of conduct: <https://www.ozgrav.org/code-of-conduct.html>

Fundamentally, all of these codes are focused on one thing: 
maintaining a professional environment where everyone is treated with respect. 


If you notice that this is not happening, and are comfortable with doing so, please speak up or let me know offline*.  We (myself included) should aim to be held accountable to each other, respect being challenged for our behavior, and remain willing to get advice on how we as individuals can help create a welcoming community where we can all thrive. Making mistakes is how we learn, both as individuals and as a community. I and others will be very grateful if you let us know about any suggestions and will endeavor to learn from them.  

But there may be situations where this can be very uncomfortable.  
With that in mind, I’ve asked Prof. Elena Maria Rossi from Leiden University in the 
Netherlands to act as the COMPAS Ombudsperson, and Elena kindly agreed to do so.  
In addition to being a brilliant astrophysicist, 
Elena has a wealth of experience in championing equality and inclusivity; 
for example, she is the founder and chair of the LISA consortium Diversity and Inclusion committee.  
If you are experiencing any issues in interactions with COMPAS colleagues and would 
like to discuss them confidentially with someone who is outside of any local structures, 
please contact Elena, her contact details can be found at <https://www.universiteitleiden.nl/en/staffmembers/elena-rossi#tab-1> 

Best wishes,
Ilya 

* see <https://ilyamandel.github.io/> for my contact information
# How to recreate/edit the main sequence lifetime plot from the COMPAS methods paper

This folder contains everything you need to reproduce (or edit!) the main sequence plot from the methods paper. This one is a little simpler than some of the other plots as it comes directly from the fitting formulae in Hurley+00 and hence can be created without running COMPAS.

Therefore all you need to do is open `make_fig_7.ipynb` in Jupyter Lab/notebook and run the whole thing to create the figure (and learn how to change it).# How to recreate/edit the maximum radius plot from the COMPAS methods paper

This folder contains everything you need to reproduce (or edit!) the maximum radius plot from the methods paper. You'll need to do the following:

1. Run `python create_fig_6_grid.py` to create a grid of stars for the plot (you can change this to a custom range of masses or metallicities if you like)
2. Run `python pythonSubmit.py` to run COMPAS for this grid (note this is probably going to take about 10 minutes!)
3. Open `make_fig_6.ipynb` in Jupyter Lab/notebook and run the whole thing to create the figure (and learn how to change it)

For reference, the changes from the *default* pythonSubmit.py are just:

- use detailed output
- use smaller time steps (to get smoother lines for the plot)# How to recreate/edit the HR diagram from the COMPAS methods paper

This folder contains everything you need to reproduce (or edit!) the HRD from the methods paper. You'll need to do the following:

1. Run `python create_fig_5_grid.py` to create a grid of stars for the plot (you can change this to a custom range of masses or metallicities if you like)
2. Run `python pythonSubmit.py` to run COMPAS for this grid
3. Open `make_fig_5.ipynb` in Jupyter Lab/notebook and run the whole thing to create the figure (and learn how to change it)

For reference, the changes from the *default* pythonSubmit.py are just:

- use detailed output
- use smaller time steps (to get smoother lines for the plot)# How to recreate/edit the initial-core-final mass relation plot from the COMPAS methods paper

This folder contains everything you need to reproduce (or edit!) the initial-core-final mass relation plot from the methods paper. You'll need to do the following:

1. Run `python create_fig_8_grids.py` to create a grid of stars for the plot
2. Run three different python submit files: `python pythonSubmitDefaults.py`,`python pythonSubmitRapid.py` and `python pythonSubmitMandelMueller.py` to run COMPAS for these three grids(note this is probably going to take quite a few minutes!)
3. Open `make_fig_8.ipynb` in Jupyter Lab/notebook and run the whole thing to create the figure (and learn how to change it)

For reference, the changes from the *default* pythonSubmit.py are just:

- change the output file name (to distinguish the 3 grids)
- change the logfiles definitions so that only the variables needed are included
- change grid name
- (for some) change the remnant mass prescription# How to recreate/edit the BBH chirp mass distribution from the COMPAS methods paper

This folder contains everything you need to reproduce the BBH chirpmass distribution from the COMPAS methods paper (_arXiv_: https://arxiv.org/abs/2109.10352;  Team COMPAS et al. 2021). 


The data used for this figure is publicly available at: https://zenodo.org/record/5655483

This data set contains the output of 10,000,000 binaries evolved using COMPAS 02.21.00, using adaptive importance sampling (STROOPWAFEL, Broekgaarden et al. 2019), sampling from a metallicity uniform in $\log(Z) \in [10^{-4},0.03]$. More details can be found in `Run_Details.txt`.

### Data reporduction
The data can be reproduced by running version `02.21.00` of COMPAS, 

1. Run `stroopwafel_interface.py`, that reads in the `Fig16_pythonSubmit.py.py` (both contained in this folder).
2. Calculate the rates by running ```FastCosmicIntegration.py```  from COMPAS's post-processing tools, with the following flags altered from their default values:


```:::bash
    python FastCosmicIntegration.py  --mu0 0.035 --muz -0.23 --sigma0 0.39 --sigmaz 0.0 --alpha 0.0 --weight mixture_weight --zstep 0.01 --sens O3 --m1min 10. --aSF 0.01 --bSF 2.77 --cSF 2.9 --dSF 4.7 
```

---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**Label the issue**
Please label the 'severity' and 'urgency' of this issue. You can choose:

`urgency_high`     - This is a very urgent issue and should be resolved as soon as possible

`urgency_moderate` - This is a moderately urgent issue

`urgency_low`      - This issue is not urgent

`severity_major`    - This is a severe bug

`severity_moderate` - This is a moderately severe bug

`severity_minor`    - This is a minor bug with minimal impact

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Versioning (please complete the following information):**
 - OS: [e.g. Ubuntu 18.04]
 - COMPAS [e.g. v02.09.02]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
COMPAS online documentation
===========================

The COMPAS online documentation is available at this URL:

    https://compas.science/docs

which redirects to the underlying readthedocs page at:

    https://compas.readthedocs.io/en/latest/index.html


The source files for the online documentation are in the 'docs' directory at:

    https://github.com/TeamCOMPAS/COMPAS


The online documentation is provided via readthedocs.  To learn more about readthedocs, visit:

    https://readthedocs.org/


readthedocs uses the Sphinx document generator.  To learn more about Sphinx, visit:

https://www.sphinx-doc.org/en/master/index.html


The file 'requirements.txt', in the 'docs' directory, lists the software/python modules required
to build the documentation.  If you plan to build the documentation locally (either so you can 
test your changes before you push them to the repository, or just to view the documentation locally
rather than online), you will need to install these requirements.

The requirements.txt file informs readthedocs what dependencies need to be installed in order for
readthedocs to build the online documentation. Some modules are commented ('#') in requirements.txt
because readthedocs either doesn't require them, or installs them by default - you will need to install
the commented modules to build the documentation locally.

Python is required, so install that if it is not already installed.  Then install 'sphinx', and the
modules listed in requirements.txt (including the ones commented).


Updating the documentation
--------------------------

The documentation source file are ReST (Restructured text) files - similar to markdown.  With readthedocs,
eash .rst file is compiled into a .html file, so viewable with a web browser. The documentation source
files (the .rst files) are in 'docs/online-docs' (only the 'requirements.txt' file is required by
readthedocs to be in the 'docs' directory - if that file is not in the 'docs' directory, or top-level
directory of the repo, readthedocs will not build the online documentation).

The 'online-docs' directory contains files and sub-directories that are structured to match the structured
of the online documentation. Find the .rst file you want to modify, and make the changes in your favourite
text editor. If you need to add new .rst files, follow the existing structure, and make sure you link the
files into an index (toctree) somewhere.


Building the documentation locally
----------------------------------

Once you have Sphinx and the dependencies (from 'requirements.txt') installed, navigate to the 'docs/online-docs'
directory in you local COMPAS repo, and type:

    make clean
    make html

If everything has been installed correctly this will first remove the existing .html files for the documentation
(make clean), and then recreate them (make html).  During the build process you may see some warnings that some
documents (.rst files) are not included in any toctree - that's ok, not all our pages are accessible via a table
of contents (toctree).

If the build completes successfully, there will be .html files in the 'docs/online-docs/\_build/html' directory.
To view the newly build documentation, open 'docs/online-docs/\_build/html/index.html' with your web browser
(e.g. type 'file://path-to-compas/docs/online-docs/\_build/html/index.html' into your web browser address bar, 
where 'path-to-compas' is the path to your local COMPAS repo).

Aside: the '\_build' directory is not required on the remote repo (and will only bloat the repo), so you should
add it to your .gitignore.


Pushing the changes online
--------------------------

Once you are satisfied with your changes, push the updated source files to the COMPAS repo as you would any source
changes.  If things work properly, that's all you need to do: readthedocs has webhooks that will notice the change
and automatically rebuild the online documentation.  The process (noticing the change, rebuilding the documentation,
then deploying the updated web pages) could take up to 15 minutes or so. If something goes wrong and the changes are
not noticed by readthedocs, or you just don't want to wait 15 minutes for your changes to appear on the web, you can
log into the readthedocs project and initiate a rebuild manually.

The COMPAS readthedocs project name is 'compas', and the project page is:

    https://readthedocs.org/projects/compas/

You need to login to the readthedocs project to do anything other than look at the dashboard. Once you are logged in 
you can rebuild the project. To manually start the rebuild, make sure you are on the 'Overview' page (select the
'Overview' button if you are not) and select the 'Build' button. The build may sometimes fail with an "environment"
error - the solution is to wait a few minutes and try again (readthedocs has some concurrency and timing limits).
Logon details can be found on the COMPAS slack workspace (devel\_compas\_documentation channel).

