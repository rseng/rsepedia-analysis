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

.. things that are needed globally ...  JR


.. newline - inserts '<br />' in the html

.. |br| raw:: html

   <br />


.. nbsp - insertt a non-breaking space in the html

.. |_| unicode:: 0xA0 
   :trim:


.. rst doesn't seem to support nesting decorations (e.g. can't do '**text**' to get 'text' in bold inside quote marks
.. this implements bold text and italic text via CSS (along with the corresponding entry in the CSS file) which can then
.. be wrapped in quote marks (or whatever else you want to wrap it in) e.g. some text ':boldtext:`bold_text`'
.. (I wanted to just use role :bold: and :italic:, but I think at least :bold: is defined elsewhere in readthedocs... 
.. so I made :bolditalictext: consistent)

.. role:: boldtext
  :class: boldtext

.. role:: italictext
  :class: italictext


.. rst doesn't support bold italics - this implements it via CSS (along with the corresponding entry in the CSS file)

.. role:: bolditalictext
  :class: bolditalictext

.. COMPAS documentation master file, created by
   sphinx-quickstart on Thu Sep  2 08:36:56 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.





COMPAS
======

**C**\ ompact **O**\ bject **M**\ ergers: **P**\ opulation **A**\ strophysics and **S**\ tatistics
is a publicly available rapid stellar / binary population synthesis code designed so that evolution prescriptions and model parameters are 
easily adjustable. COMPAS draws properties for a binary star system from a set of initial distributions, and evolves it from 
zero-age main sequence to the end of its life as two compact remnants. It has been used for inference from observations of 
gravitational-wave mergers, Galactic neutron stars, X-ray binaries, and luminous red novae.

`... by` `TeamCOMPAS <https://compas.science/index.html>`_



Contents
--------

.. toctree::
   :maxdepth: 1

   ./pages/Getting started/getting-started

   ./pages/User guide/user-guide
   ./pages/Developer guide/developer-guide

    Code repository <https://github.com/TeamCOMPAS/COMPAS>
    Request an enhancement <https://github.com/TeamCOMPAS/COMPAS/issues/new>
    Report a problem <https://github.com/TeamCOMPAS/COMPAS/issues/new>
   
   ./pages/references
   ./pages/contact-us


If you use COMPAS in the preparation of a publication, please :doc:`cite COMPAS <./pages/how-to-cite>`.

Licence
-------

COMPAS is available under the `MIT <https://choosealicense.com/licenses/mit/>`_ licence.


Citing COMPAS
-------------

If you use this code or parts of this code for results presented in a scientific publication, we would greatly appreciate you sending
us your paper reference and making your input settings and output data publicly available by uploading it to the COMPAS Zenodo community. 

Please also cite: 

.. _cite-compas:

    Team COMPAS: J. Riley `et al.` [:cite:year:`compas2021`]. |_| |_| |_| |_| |_| |_| :download:`Bibtex citation <../COMPAS-2021methodsPaper.bib>`

|br|
We would also greatly appreciate an acknowledgement of the form:

"Simulations in this paper made use of the COMPAS rapid binary population synthesis code (version x.y.z), which is freely available at
http://github.com/TeamCOMPAS/COMPAS."

|br|
Furthermore,

If using the COMPAS model of gravitational wave selection effects, please cite :cite:t:`Barrett2018`.

If you use COMPAS's importance sampling algorithm STROOPWAFEL, please cite :cite:t:`Broekgaarden2019`.

If using COMPAS's integration over cosmic star formation history, please cite :cite:t:`Neijssel2019`.

If using the COMPAS model of (pulsational) pair instability supernova, please cite :cite:t:`Stevenson2019`.

If evolving pulsar spins and magnetic fields with COMPAS, please cite :cite:t:`Chattopadhyay2020`.

If you use the COMPAS model of chemically homogeneous evolution, please cite :cite:t:`Riley2021`.

Contact us
==========

Join the `COMPAS User Google Group <https://groups.google.com/forum/#!members/compas-user>`__ to engage in discussions with 
COMPAS users and developers, or email your queries to compas-user@googlegroups.com.

References
==========

.. bibliography::
   :all:

Grid files
==========

A grid file allows users to specify initial values for multiple systems for both Single Star Evolution (SSE) and Binary Star Evolution 
(BSE).  Each line of a grid file is used by COMPAS to set the initial conditions and evolutionary parameters for an individual single 
star (SSE) or binary star (BSE), and each single star or binary star defined by a grid file line is evolved using those values.

Each line of a grid file is a set of program option specifications, with the specifications being exactly as they would appear on the 
command line if running COMPAS from the command line.

For example, a grid file could contain the following two lines:

    --metallicity 0.001 --eccentricity 0.0 --remnant-mass-prescription fryer2012 |br|
    --remnant-mass-prescription mullermandel --metallicity 0.02 --semi-major-axis 45.678

in which case COMPAS would evolve two binaries, with the option values set per the grid file lines.

Grid files can have blank lines and comments. Comments begin with a hash/pound character ('#') - the hash/pound character and text 
beyond it are ignored by COMPAS.

Not all program options can be specified in a grid file. Options that should remain constant for a single execution of COMPAS, such as
options that specify the mode of evolution (e.g. ``--mode``), or the name or path of output files (e.g. ``--output-path``, 
``--logfile-detailed-output`` etc.) can only be specified on the command line.  COMPAS will issue an error message if an option that is
not supported in a grid file is specified on a grid file line.

COMPAS imposes no limit on the number of grid file lines: the size of a grid file is limited only by the filesystem of the host system.


Specifying a Subset of the Grid File to be Processed
----------------------------------------------------

Users can instruct COMPAS to process only a subset of a specified grid file. This is achieved via the program options ``--grid-start-line``,
and ``--grid-lines-to-process``:

    The ``--grid-start-line`` program option takes a single parameter: an integer specifying the zero-based line number of the first line of 
    the grid file to be processed. The default value is 0 - the first line of the grid file.  Specifying a start line beyond the end of the 
    grid file will result in an UNEXPECTED-END-OF-FILE error.

    The ``--grid-lines-to-process`` program option takes a single parameter: an integer specifying the number of lines of the grid file to be 
    processed. The default is to process all lines in the grid file from the start line (which may have been specified by the 
    ``--grid-start-line option``) through to the end of the grid file. Specifying a number of lines to be processed that, when coupled with 
    the start line, would result in attempting to process lines beyond the end of the grid file will result in an INCOMPLETE-GRID error.
    
    Note that blank lines and comments count towards the number of grid lines processed by COMPAS when deciding if the number specified by the
    user has been reached.

Both option ``--grid-start-line`` and ``--grid-lines-to-process`` are ignored if no grid file is specified via the ``--grid`` program option.
Random seed
===========

The ``--random-seed`` option allows users to specify the initial value to be used to seed the pseudo-random number generator. Once set, the
random seed values increments from its initial value for each star, or binary star, evolved. How the random seed increments depends upon the
context.

The ``--random-seed`` option can be specified on either, or both, the command line and a :doc:`grid file <./grid-files>` line. If the option 
is not specified on one or the other, the default value is used (see :doc:`./Program options/program-options-list-defaults`).

In general, if the ``--random-seed`` option is specified, the pseudo-random number generator will be seeded using the specified value for 
the first star, or binary star, evolved, then for each subsequent star or binary star, the seed value will be incremented by one and the 
pseudo-random number generator re-seeded. Seeding the pseudo-random number generator with a known seed for each star, or binary star, 
evolved ensures that the evolution of specific stars, or binary stars, can be reproduced.

Consider a single execution of COMPAS effected with the command::

    ./COMPAS --random-seed 15 --number-of-systems 100 --metallicity 0.015

This would evolve 100 binary stars, each with `metallicity = 0.015`, and other initial attributes set to their defaults. The first of the 
100 binary stars will be evolved using the random seed 15, the second 16, the third 17, and so on - each binary star will evolve using
a unique random seed.

In the example shown above (see Section :doc:`./Program options/program-options-mixing-ranges-sets`), all 104 binary stars would evolve with
unique random seed values, ranging from 0 (the default, since the option was not specified on either the command line or in the grid file), to 103.

In both these examples, the random seed was incremented in the context of the command line. In the first example, the random seed was 
explicitly specified on the command line, and in the second example the random seed defaulted to the command line default.

Consider now a single execution of COMPAS, using the grid file ``mygrid.txt``::

    ./COMPAS --random-seed 12 --grid mygrid.txt

where the contents of the grid file ``mygrid.txt`` are::

    --allow-rlof-at-birth true --metallicity 0.1
    --semi-major-axis 23.4 --random-seed 107
    --random-seed 63 --metallicity 0.12 --eccentricity s[0.1,0.2,0.3,0.4]
    --initial-mass-1 12.3

This would evolve 7 binary stars with random seed values 12, 107, 63, 64, 65, 66, and 18.

The first binary star evolved is the first line of the grid file. This line does not specify the ``--random-seed`` option, so the random seed
defaults to the command line value. The command line did specify a value for the random seed (12), so that value is used. Since the first
line of the grid file is the first binary star evolved, the random seed is not incremented, and the value of 12 is used.

The second binary star evolved is the second line of the grid file. This line does specify the ``--random-seed`` option. Since this is the 
first binary star evolved in the context of the random seed specified on the grid file line, the random seed is not incremented, and the value 
of 107 is used.

The third binary star evolved is the third line of the grid file. This line does specify the ``--random-seed`` option. Since this is the first
binary star evolved in the context of the random seed specified on the grid file line, the random seed is not incremented, and the value of 63
is used.

The fourth, fifth, and sixth binary stars evolved are also from the third line of the grid file - a set of four values for eccentricity was 
specified. Since these are subsequent to the first binary star evolved in the context of the random seed specified on the grid file line, the
random seed is incremented, and the values of 64, 65, and 66 are used.

The seventh binary star evolved is the fourth line of the grid file. This line does not specify the ``--random-seed`` option, so the random 
seed defaults to the command line value. The command line did specify a value for the random seed (12), so that value is used, but since this 
binary star is subsequent to the first binary star evolved in the context of the random seed specified on the command line, the random seed is 
incremented. This is the sixth subsequent binary star evolved in the context of the command line (all stars, or binary stars, evolved in a 
single execution of COMPAS are evolved in the context of the command line), so the random seed is incremented from 12 to 18 (by 1 for each 
binary star evolved), and the value used for this binary star is 18.

Note that in this example, all binary stars were evolved using a unique random seed. This is because the values specified for the random seeds 
via the ``--random-seed`` option were ’well-behaved’. Unfortunately there is no reasonable way to protect the user against specifying duplicate 
random seeds – especially since the random seed increments for each star or binary star. If the user chooses to specify multiple grid file lines 
with the same random seed, or initial random seeds that would collide with other random seed values and cause duplicates as they increment 
through ranges and sets, then there will be duplicate random seeds in the output files. Users should take care when specifying random seeds in 
grid files via the ``--random-seed`` option.
User guide
==========

This section contains the basic user guide for COMPAS.


.. toctree::
   :maxdepth: 1

   ./configuration
   ./Program options/program-options
   ./grid-files
   ./random-seed
   ./Running COMPAS/running-compas
   ./COMPAS output/output
   ./Post-processing/post-processing
   ./Tutorial/example-compas-run
   ./sampling
   ./docker
Configuration
=============

Run-time Configuration
----------------------

COMPAS is configured at run-time via command-line options, referred to as "program options" in this documentation, and, optionally, 
grid files.

Configuring COMPAS via command-line options and grid files gives users the flexibility to specify both the initial conditions and the 
evolutionary parameters that define single and binary stars at birth, and the conditions under which the stars and systems evolve over 
their lifetime.

COMPAS has a rich set of command-line options that can be set by the user, and these, coupled with the flexibility afforded by the 
COMPAS grid file implementation, allow users to configure COMPAS over a broad range of initial conditions and evolutionary parameters.
This provides users with the means to explore a broad range of physics and physical assumptions.


Compile-time Configuration
--------------------------

The values of some physical constants, bounds for some initial conditions, evolutionary parameters, and physical processes, etc., are 
specified in the COMPAS source file `constants.h`.  While it is unlikely that these constants would need to be changed in most ordinary 
COMPAS runs, the possibility exists that users may want to change some of them.  Should that be the case, the user should change the 
value(s) required in ``constants.h`` and rebuild the COMPAS executable. A makefile is provided in the source directory.

See :doc:`../Getting started/building-COMPAS` for details of how to build the COMPAS executable.
Sampling in COMPAS
==================



Here are some basic instructions for efficient sampling of the COMPAS
input parameters, using the python sampling package Stroopwafel.

Note that the intended Stroopwafel functionality for "Adaptive
Importance Sampling" is not yet implemented, but is currently in
development.

Requirements




If you have not already, you will need to install Stroopwafel. If you
have admin rights, Stroopwafel can be installed on your system with
``pip install stroopwafel``.

Instructions




To use Stroopwafel sampling, copy
``preProcessing/stroopwafelInterface.py`` into your working directory.

Settings
~~~~~~~~

NOTE: This sampling method is currently being updated as part of an upgrade
in our method to parse user-defined options. We plan to address this shortly. 
Please bear with us and contact the COMPAS team if an urgent solution is needed.

1. runSubmit


If you are running COMPAS on default settings, skip this section.

If you have many non-default COMPAS arguments, you may want to set
them in the ``compasConfigDefault.yaml``, that is read and executed by the 
``runSubmit.py`` file in the same directory. For now, the file must
be named this way and placed in the same directory as the ``stroopwafelInterface.py``
file.

A configurable runSubmit file can be found in the ``preProcessing/``
directory.

Set your desired options, then set the ``usePythonSubmit`` parameter to ``True``
in the ``stroopwafelInterface.py``.

2. Stroopwafel inputs


The lines below ``usePythonSubmit`` represent stroopwafel inputs.
These are treated as
defaults, but can be overriden by command-line arguments to
stroopwafel.
See ``python3 stroopwafelInterface.py --help``.

``num_systems`` is the total number of binaries you would like to
evolve.
This value overrides the value set in the ``compasConfigDefault.yaml`` file.

``num_cores`` is the number of cores you would like to use to
parallelize your run. More cores means your run will finish sooner, but
may reduce your ability to run other tasks while you wait. On linux
systems, the command ``echo $(nproc)`` will tell you how many (virtual)
CPUs you have available.

``num_per_core`` is the number of systems to run on a core at a given
time. This translates to the number of systems in a single batch file.
This is more relevent for adaptive importance sampling.

``mc_only`` specifies if you would like to do naive MC sampling only.
Currently, this option must be set to True

``run_on_hpc`` specifies if you are running on a High-Performance
Computer (HPC).
If so, see `docs/compasHPC.md <compasHPC.md>`__ for assistance.

``output_folder`` a string specifying the output folder. Relative paths
will be appended onto the current directory path.

``output_filename`` a string specifying the name of the output samples
file.

``debug`` whether to print the COMPAS output/error.

3. Sampling parameters


Sampled parameters will be combined into grid files which COMPAS then
reads in.
Users should choose which parameters they would like to be sampled
over, as well as
the relevent distributions.

See the `COMPAS
Documentation <https://github.com/TeamCOMPAS/COMPAS/blob/Documentation/COMPAS_Documentation.pdf>`__
for details on which sets of
parameters are allowed/required.

See `Stroopwafel
Documentation <https://github.com/lokiysh/stroopwafel>`__ for details on
which distributions are available.

4. Run Stroopwafel


When your satisfied with your settings, simply run with
``python3 stroopwafelInterface.py``. The output will be collected into
batch containers in your output folder.
To postprocess the output, see
`getting\_started.md <getting_started.md>`__
COMPAS and Docker
=================

Docker has been added to COMPAS to reduce time and effort required to
set up the COMPAS deployment environment.

Instead of having to install and configure several libraries and tools
(e.g. python/pip, numpy, g++, boost) which can vary considerably beween
operating systems and existing toolchains, users can instead opt to
install Docker and run COMPAS with a single command.

This also gives users the ability to run COMPAS on cloud solutions like
`AWS EC2 <https://aws.amazon.com/ec2/>`__ or `Google Compute
Engine <https://cloud.google.com/compute>`__ where hundreds of cores can
be provisioned without having to manually configure the environment.

Docker works by creating an isolated and standalone environment known
as a `container <https://www.docker.com/resources/what-container>`__.
Containers can be created or destroyed without affecting the host
machine or other containers\*.

Containers are instances of images. An image is a pre-defined
setup/environment that is instantiated when started as a container
(containers are to images what objects are to classes). More
`here <https://stackoverflow.com/questions/23735149/what-is-the-difference-between-a-docker-image-and-a-container#:~:text=An%20instance%20of%20an%20image,of%20layers%20as%20you%20describe.&text=You%20can%20see%20all%20your,an%20image%20is%20a%20container.>`__
on the relationship between images and container.

Containers are (almost) always run as a Linux environment. A major
benefit of this is the ability to run Linux applications in a Windows or
MacOS environment without having to jump through hoops or have a
diminished experience.

Image definitions can be defined by users (e.g. Dockerfiles); there are
also standard images publicly available on `Docker
Hub <https://hub.docker.com/>`__

All that is required to start using COMPAS with Docker is the "Usage"
section (the "CI/CD" section is also highly recommended).
The other sections are provided for extra info.

\* Containers can still interact with each other and the host machine
through mounted directories/files or exposed ports.



Usage


N.B. This section assumes `Docker <https://www.docker.com/>`__ has
been installed and is running.
For Windows and MacOS users, see
`here <https://www.docker.com/products/docker-desktop>`__.

Installing


The latest compiled version of COMPAS (dev branch) can be retrieved by
running
``docker pull teamcompas/compas``

Other versions can be used by adding a version
`tag <https://docs.docker.com/engine/reference/commandline/tag/>`__.
For example, COMPAS version 2.12.0 would be
``teamcompas/compas:2.12.0``.
To see all available versions, go to the TeamCOMPAS docker hub page
`here <https://hub.docker.com/u/teamcompas>`__.

Running


COMPAS can still be configured via command line arguments passed to the
COMPAS executable or via a ``runSubmit.py`` file.

Running runSubmit.py


To run COMPAS via a ``runSubmit.py`` file, the command is a little
more complex.

::

    docker run                                                  \
    --rm                                                    \
    -it                                                     \
    -v $(pwd)/compas-logs:/app/COMPAS/logs                  \
    -v $(pwd)/runSubmit.py:/app/starts/runSubmit.py   \
    -e COMPAS_EXECUTABLE_PATH=/app/COMPAS/bin/COMPAS        \
    -e COMPAS_LOGS_OUTPUT_DIR_PATH=/app/COMPAS/logs         \
    teamcompas/compas                                       \
    python3 /app/starts/runSubmit.py                     

Breaking down this command:

``docker run``
creates a container

``--rm``
`Clean
up <https://docs.docker.com/engine/reference/run/#clean-up---rm>`__
destroy the container once it finishes running the command

``-it``
short for `-i and
-t <https://docs.docker.com/engine/reference/run/#foreground>`__ -
provides an interactive terminal

``-v <path-on-host>:<path-in-container>``
`Bind mounts <https://docs.docker.com/storage/bind-mounts/>`__
mount ``<path-on-host>`` to ``<path-in-container``>
This time we not only want to get the output from COMPAS on the host
machine, we also want to supply a ``runSubmit.py`` to the container
from the host machine.

NOTE: if you decide to execute using ``runSubmit.py``, you will need 
a ``compasConfigDefault.yaml``  file in the same directory. This file 
can be find in the same directory as the ``runSubmit.py``, and contains
the default COMPAS choices for stellar and binary physics. These choices
can be changed by modifying the options availabe in the ``compasConfigDefault.yaml`` 
file.

``-e VAR_NAME=value``
`Environment
variables <https://docs.docker.com/engine/reference/run/#env-environment-variables>`__
set the environment variable ``VAR_VAME`` to ``value``

``teamcompas/compas``
the image to run

``python3 /app/starts/runSubmit.py``
the command to run when the container starts

Run the COMPAS executable


To run the COMPAS executable directly (i.e. without ``runSubmit.py``)

::

    docker run                                  \
    --rm                                    \
    -it                                     \
    -v $(pwd)/compas-logs:/app/COMPAS/logs  \
    teamcompas/compas                       \
    bin/COMPAS                              \
    --number-of-binaries=5                  \
    --outputPath=/app/COMPAS/logs

Breaking down this command:

``docker run``
creates a container

``--rm``
`Clean
up <https://docs.docker.com/engine/reference/run/#clean-up---rm>`__
destroy the container once it finishes running the command

``-it``
short for `-i and
-t <https://docs.docker.com/engine/reference/run/#foreground>`__ -
provides an interactive terminal

``-v <path-on-host>:<path-in-container>``
`Bind mounts <https://docs.docker.com/storage/bind-mounts/>`__
mount ``<path-on-host>`` to ``<path-in-container>``
In this instance, make it so
``$(pwd)/compas-logs on my machine is the same as``/app/COMPAS/logs\`
inside the container

``teamcompas/compas``
the image to run

``bin/COMPAS``
the command to run when the container starts

``--number-of-binaries``
anything after the given start command is passed to that command, in
this case, the flag to set the number of binaries

``--outputPath /app/COMPAS/logs``
same as above, anthing after the start command is given to that start
command, here it forces logs to go to the directory that is mapped to
the host machine

More info on ``docker run``
`here <https://docs.docker.com/engine/reference/run/>`__

NOTE 1:

Two new environment variables have been added, both of these apply to
``runSubmit.py`` only and are non-breaking changes.

``COMPAS_EXECUTABLE_PATH`` is an addition to the default
``runSubmit.py`` that overrides where ``runSubmit.py`` looks for
the compiled COMPAS.
This override exists purely for ease-of-use from the command line.

``COMPAS_LOGS_OUTPUT_DIR_PATH`` is also an addition to the default
``runSubmit.py`` that overrides where logs are placed.
The override exists because the mounted directory (option ``-v``) is
created before COMPAS runs. COMPAS sees that the directory where it's
supposed to put logs already exists, so it created a different (i.e.
non-mapped) directory to deposit logs in.

NOTE 2:

The ``docker run ...`` examples above both use the ``-it`` options.
If you want to run multiple instances of COMPAS, I would highly
recommend using `detached
mode <https://docs.docker.com/engine/reference/run/#detached--d>`__
(``-d``) instead.
All container output will be hidden.

An example where this would be useful is if you were running 4
instances of COMPAS at once.
You could copy/paste the following into the terminal...

::

    docker run --rm -d -v $(pwd)/compas-logs/run_0:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_01.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    
    docker run --rm -d -v $(pwd)/compas-logs/run_1:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_02.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    
    docker run --rm -d -v $(pwd)/compas-logs/run_2:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_03.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    
    docker run --rm -d -v $(pwd)/compas-logs/run_3:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_04.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py

...which would run 4 separate instances of COMPAS, each with its own
``runSubmit.py`` file and logging directory, and all console output
supressed.

You may want to check the console output to see how far into the run
COMPAS is.
The command for this is ``docker logs <container_id>``.
You can get the container id by running ``docker ps``.



CI/CD


The latest version of COMPAS (dev branch) is available at
``teamcompas/compas``.
This is provided automatically by CI/CD.

Whenever a push to
`TeamCOMPAS/dev <https://github.com/TeamCOMPAS/COMPAS/tree/dev>`__ a
continuous deployment process automatically
`builds <https://docs.docker.com/engine/reference/commandline/build/>`__
a new image and deploys it to DockerHub with a ``tag`` that corresponds
to the value of ``VERSION_STRING`` in ``constants.h``.

At time of writing, `GitHub
Actions <https://github.com/features/actions>`__ is facilitating the
above process. While this is convenient (because it's free and well
supported) it is quite slow. I have plans to create a
`runner <https://help.github.com/en/actions/getting-started-with-github-actions/core-concepts-for-github-actions#runner>`__
locally with a high core count that can be used to compile COMPAS
quickly, but haven't gotten around to it yet.

You can realistically expect the latest COMPAS docker image to be
available 5 - 10 minutes after pushing/merging.

The Github Actions configuration is in
``/.github/workflows/dockerhub-ci.yml``.

Atlassian has a `good
writeup <https://www.atlassian.com/continuous-delivery/principles/continuous-integration-vs-delivery-vs-deployment>`__
about what CI/CD is.



Bonus Info


Dockerfile


The `Dockerfile <https://docs.docker.com/engine/reference/builder/>`__
defines how the docker image is constructed.

Images are created as a combination of layers.
During the build process each layer is cached and only updated on
subsequent builds if that layer would change.

The Dockerfile for COMPAS is made up of 8 layers.

``FROM ubuntu:18.04``
Use `Ubuntu 18.04 <https://hub.docker.com/_/ubuntu>`__ as a base
(provided by Docker Hub)
`https://docs.docker.com/engine/reference/builder/#from <FROM>`__ docs

``WORKDIR /app/COMPAS``
Effectively ``cd /app/COMPAS`` within the container.
`WORKDIR <https://docs.docker.com/engine/reference/builder/#workdir>`__
docs

``RUN apt-get update && apt-get install -y ...``
Install the required dependencies.
``-y`` so there's no prompt to install any of the packages.
``update`` and ``install`` are in the same layer because now if there
are any updates, it will force all of the dependencies to be
re-installed
`RUN <https://docs.docker.com/engine/reference/builder/#run>`__ docs

``RUN pip3 install numpy``
Install numpy.
`RUN <https://docs.docker.com/engine/reference/builder/#run>`__ docs

``COPY src/ src/``
Copy ``./src/`` directory from the local machine to ``./src`` in the
container (remembering that ``WORKDIR`` changes the cwd).
`COPY <https://docs.docker.com/engine/reference/builder/#copy>`__ docs

``RUN mkdir obj bin logs``
Create the directories required by COMPAS.
`RUN <https://docs.docker.com/engine/reference/builder/#run>`__ docs

``ENV COMPAS_ROOT_DIR /app/COMPAS``
Set the required environment variable(s).
`ENV <https://docs.docker.com/engine/reference/builder/#env>`__ docs

``RUN cd src && make -f Makefile.docker -j $(nproc)``
Make COMPAS using a specific makefile (more below) and as many cores
as possible.
`RUN <https://docs.docker.com/engine/reference/builder/#run>`__ docs

Dockerfiles will usually end with a ``CMD`` directive that specifies
what command should run when the container is started.
COMPAS doesn't have a ``CMD`` directive because some users will want
to run the executable directly and some will want to use
``runSubmit.``.
`CMD <https://docs.docker.com/engine/reference/builder/#cmd>`__ docs

Makefile.docker


A separate makefile is required for Docker in this scenario for two
reasons.

#. To separate compiled files from source files
#. To prevent the usage of ``-march=native``

``-march=native`` is a fantastic optimisation for users who compile
and run COMPAS on the same machine, however it causes fatal errors when
running COMPAS on a machine that it was not compiled for.
`Docs <https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html>`__ for
``-march``.

This selects the CPU to generate code for at compilation time by
determining the processor type of the **compiling machine**.

Using -march=native enables all instruction subsets supported by the
local machine (hence the result might not run on different
machines).



Program option properties
=========================

When specifying known properties in a log file record specification record, the property name must be prefixed with 
the property type. Currently there is a single binary property type available for use: PROGRAM_OPTION.

For example, to specify the program option property ``RANDOM_SEED``, use::

    PROGRAM OPTION::RANDOM_SEED


.. _spec-options-props-top:

Following is an alphabetical list of program option properties available for inclusion in log file record specifiers.

**Jump to**
:ref:`A <spec-options-props-A>` :ref:`B <spec-options-props-B>` :ref:`C <spec-options-props-C>` :ref:`D <spec-options-props-D>`
:ref:`E <spec-options-props-E>` :ref:`F <spec-options-props-F>` :ref:`G <spec-options-props-G>` :ref:`H <spec-options-props-H>`
:ref:`I <spec-options-props-I>` :ref:`J <spec-options-props-J>` :ref:`K <spec-options-props-K>` :ref:`L <spec-options-props-L>`
:ref:`M <spec-options-props-M>` :ref:`N <spec-options-props-N>` :ref:`O <spec-options-props-O>` :ref:`P <spec-options-props-P>`
:ref:`Q <spec-options-props-Q>` :ref:`R <spec-options-props-R>` :ref:`S <spec-options-props-S>` :ref:`T <spec-options-props-T>`
:ref:`U <spec-options-props-U>` :ref:`V <spec-options-props-V>` :ref:`W <spec-options-props-W>` :ref:`X <spec-options-props-X>`
:ref:`Y <spec-options-props-Y>` :ref:`Z <spec-options-props-Z>`

.. _spec-options-props-A:

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ADD_OPTIONS_TO_SYSPARMS**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_AddOptionsToSysParms
   * - Description:
     - Value of program option ``--add-options-to-sysparms``
   * - Header String:
     - Add_Options_To_SysParms

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ALLOW_MS_STAR_TO_SURVIVE_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_AllowMainSequenceStarToSurviveCommonEnvelope
   * - Description:
     - Value of program option ``--common-envelope-allow-main-sequence-survive``
   * - Header String:
     - Allow_MS_To_Survive_CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ALLOW_RLOF_AT_BIRTH**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_AllowRLOFAtBirth
   * - Description:
     - Value of program option ``--allow-rlof-at-birth``
   * - Header String:
     - Allow_RLOF@\ Birth

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ALLOW_TOUCHING_AT_BIRTH**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_AllowTouchingAtBirth
   * - Description:
     - Value of program option ``--allow-touching-at-birth``
   * - Header String:
     - Allow_Touching@\ Birth

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ANG_MOM_CONSERVATION_DURING_CIRCULARISATION**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_AngularMomentumConservationDuringCircularisation
   * - Description:
     - Value of program option ``--angular-momentum-conservation-during-circularisation``
   * - Header String:
     - Conserve_AngMom@\ Circ

.. _spec-options-props-B:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BLACK_HOLE_KICKS**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_BlackHoleKicks
   * - Description:
     - Value of program option ``--black-hole-kicks``
   * - Header String:
     - BH_Kicks

.. _spec-options-props-C:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CASE_BB_STABILITY_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_CaseBBStabilityPrescription
   * - Description:
     - Value of program option ``--case-BB-stability-prescription``
   * - Header String:
     - BB_Mass_xFer_Stblty_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CHE_MODE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_CheMode
   * - Description:
     - Value of program option ``--chemically-homogeneous-evolution``
   * - Header String:
     - CHE_Mode

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CIRCULARISE_BINARY_DURING_MT**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - Options::m_CirculariseBinaryDuringMassTransfer
   * - Description:
     - Value of program option ``--circularise-binary-during-mass-transfer``
   * - Header String:
     - Circularise@\ MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_ALPHA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeAlpha
   * - Description:
     - Value of program option ``--common-envelope-alpha``
   * - Header String:
     - CE_Alpha

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_ALPHA_THERMAL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeAlphaThermal
   * - Description:
     - Value of program option ``--common-envelope-alpha-thermal``
   * - Header String:
     - CE_Alpha_Thermal

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_LAMBDA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeLambda
   * - Description:
     - Value of program option ``--common-envelope-lambda``
   * - Header String:
     - CE_Lambda

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_LAMBDA_MULTIPLIER**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeLambdaMultiplier
   * - Description:
     - Value of program option ``--common-envelope-lambda-multiplier``
   * - Header String:
     - CE_Lambda_Multiplier

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_LAMBDA_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_CommonEnvelopeLambdaPrescription
   * - Description:
     - Value of program option ``--common-envelope-lambda-prescription``
   * - Header String:
     - CE_Lambda_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_MASS_ACCRETION_CONSTANT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeMassAccretionConstant
   * - Description:
     - Value of program option ``--common-envelope-mass-accretion-constant``
   * - Header String:
     - CE_Mass_Accr_Constant

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_MASS_ACCRETION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeMassAccretionMax
   * - Description:
     - Value of program option ``--common-envelope-mass-accretion-max``
   * - Header String:
     - CE_Mass_Accr_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_MASS_ACCRETION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeMassAccretionMin
   * - Description:
     - Value of program option ``--common-envelope-mass-accretion-min``
   * - Header String:
     - CE_Mass_Accr_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_MASS_ACCRETION_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_CommonEnvelopeMassAccretionPrescription
   * - Description:
     - Value of program option ``--common-envelope-mass-accretion-prescription``
   * - Header String:
     - CE_Mass_Accr_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_RECOMBINATION_ENERGY_DENSITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeRecombinationEnergyDensity
   * - Description:
     - Value of program option ``--common-envelope-recombination-energy-density``
   * - Header String:
     - CE_Recomb_Enrgy_Dnsty

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_SLOPE_KRUCKOW**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_CommonEnvelopeSlopeKruckow
   * - Description:
     - Value of program option ``--common-envelope-slope-kruckow``
   * - Header String:
     - CE_Slope_Kruckow

.. _spec-options-props-D:

.. _spec-options-props-E:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_Eccentricity
   * - Description:
     - Value of program option ``--eccentricity``
   * - Header String:
     - Eccentricity

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_EccentricityDistribution
   * - Description:
     - Value of program option ``--eccentricity-distribution``
   * - Header String:
     - Eccentricity_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_DISTRIBUTION MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_EccentricityDistributionMax
   * - Description:
     - Value of program option ``--eccentricity-max``
   * - Header String:
     - Eccentricity_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_EccentricityDistributionMin
   * - Description:
     - Value of program option ``--eccentricity-min``
   * - Header String:
     - Eccentricity_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EDDINGTON_ACCRETION_FACTOR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_EddingtonAccretionFactor
   * - Description:
     - Value of program option ``--eddington-accretion-factor``
   * - Header String:
     - Eddington_Accr_Factor

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ENVELOPE_STATE_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_EnvelopeStatePrescription
   * - Description:
     - Value of program option ``--envelope-state-prescription``
   * - Header String:
     - Envelope_State_Prscrptn

.. _spec-options-props-F:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **FRYER_SUPERNOVA_ENGINE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_FryerSupernovaEngine
   * - Description:
     - Value of program option ``--fryer-supernova-engine``
   * - Header String:
     - Fryer_SN_Engine

.. _spec-options-props-G:

.. _spec-options-props-H:

.. _spec-options-props-I:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMass
   * - Description:
     - Value of program option ``--initial-mass``
   * - Header String:
     - Initial_Mass

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMass1
   * - Description:
     - Value of program option ``--initial-mass-1``
   * - Header String:
     - Initial_Mass(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMass2
   * - Description:
     - Value of program option ``--initial-mass-2``
   * - Header String:
     - Initial_Mass(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_FUNCTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_InitialMassFunction
   * - Description:
     - Value of program option ``--initial-mass-function``
   * - Header String:
     - Initial_Mass_Function

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_FUNCTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMassFunctionMax
   * - Description:
     - Value of program option ``--initial-mass-max``
   * - Header String:
     - Initial_Mass_Func_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_FUNCTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMassFunctionMin
   * - Description:
     - Value of program option ``--initial-mass-min``
   * - Header String:
     - Initial_Mass_Func_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_MASS_FUNCTION_POWER**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_InitialMassFunctionPower
   * - Description:
     - Value of program option ``--initial-mass-power``
   * - Header String:
     - Initial_Mass_Func_Power

.. _spec-options-props-J:

.. _spec-options-props-K:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_DIRECTION_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_KickDirectionDistribution
   * - Description:
     - Value of program option ``--kick-direction``
   * - Header String:
     - Kick_Direction_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_DIRECTION_POWER**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickDirectionPower
   * - Description:
     - Value of program option ``--kick-direction-power``
   * - Header String:
     - Kick_Direction_Power

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_SCALING_FACTOR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickScalingFactor
   * - Description:
     - Value of program option ``--kick-scaling-factor``
   * - Header String:
     - Kick_Scaling_Factor

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitude
   * - Description:
     - Value of program option ``--kick-magnitude``
   * - Header String:
     - Kick_Magnitude

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitude1
   * - Description:
     - Value of program option ``--kick-magnitude-1``
   * - Header String:
     - Kick_Magnitude(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitude2
   * - Description:
     - Value of program option ``--kick-magnitude-2``
   * - Header String:
     - Kick_Magnitude(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_KickMagnitudeDistribution
   * - Description:
     - Value of program option ``--kick-magnitude-distribution``
   * - Header String:
     - Kick_Magnitude_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_MAXIMUM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitudeDistributionMaximum
   * - Description:
     - Value of program option ``--kick-magnitude-max``
   * - Header String:
     - Kick_Magnitude_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_kickMagnitudeDistributionSigmaCCSN BH
   * - Description:
     - Value of program option ``--kick-magnitude-sigma-CCSN-BH``
   * - Header String:
     - Sigma_Kick_CCSN_BH

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_kickMagnitudeDistributionSigmaCCSN NS
   * - Description:
     - Value of program option ``--kick-magnitude-sigma-CCSN-NS``
   * - Header String:
     - Sigma_Kick_CCSN_NS

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_kickMagnitudeDistributionSigmaForECSN
   * - Description:
     - Value of program option ``--kick-magnitude-sigma-ECSN``
   * - Header String:
     - Sigma_Kick_ECSN

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_kickMagnitudeDistributionSigmaForUSSN
   * - Description:
     - Value of program option ``--kick-magnitude-sigma-USSN``
   * - Header String:
     - Sigma_Kick_USSN

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MEAN_ANOMALY_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMeanAnomaly1
   * - Description:
     - Value of program option ``--kick-mean-anomaly-1``
   * - Header String:
     - Kick_Mean_Anomaly(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MEAN_ANOMALY_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMeanAnomaly2
   * - Description:
     - Value of program option ``--kick-mean-anomaly-2``
   * - Header String:
     - Kick_Mean_Anomaly(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_RANDOM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitudeRandom
   * - Description:
     - Value of program option ``--kick-magnitude-random``
   * - Header String:
     - Kick_Magnitude_Random

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_RANDOM_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitudeRandom1
   * - Description:
     - Value of program option ``--kick-magnitude-random-1``
   * - Header String:
     - Kick_Magnitude_Random(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE_RANDOM_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickMagnitudeRandom2
   * - Description:
     - Value of program option ``--kick-magnitude-random-2``
   * - Header String:
     - Kick_Magnitude_Random(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_PHI_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickPhi1
   * - Description:
     - Value of program option ``--kick-phi-1``
   * - Header String:
     - Kick_Mean_Phi(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_PHI_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickPhi2
   * - Description:
     - Value of program option ``--kick-phi-2``
   * - Header String:
     - Kick_Mean_Phi(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_THETA_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickTheta1
   * - Description:
     - Value of program option ``--kick-theta-1``
   * - Header String:
     - Kick_Mean_Theta(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_THETA_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_KickTheta2
   * - Description:
     - Value of program option ``--kick-theta-2``
   * - Header String:
     - Kick_Mean_Theta(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LBV_FACTOR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_LuminousBlueVariableFactor
   * - Description:
     - Value of program option ``--luminous-blue-variable-multiplier``
   * - Header String:
     - LBV_Factor

.. _spec-options-props-L:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LBV_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_LuminousBlueVariablePrescription
   * - Description:
     - Value of program option ``--luminous-blue-variable-prescription``
   * - Header String:
     - LBV_Mass_Loss_Prscrptn

.. _spec-options-props-M:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_LOSS_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassLossPrescription
   * - Description:
     - Value of program option ``--mass-loss-prescription``
   * - Header String:
     - Mass_Loss_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_RATIO**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassRatio
   * - Description:
     - Value of program option ``-``-mass-ratio``
   * - Header String:
     - Mass_Ratio

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_RATIO_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassRatioDistribution
   * - Description:
     - Value of program option ``--mass-ratio-distribution``
   * - Header String:
     - Mass_Ratio_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_RATIO_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassRatioDistributionMax
   * - Description:
     - Value of program option ``--mass-ratio-max``
   * - Header String:
     - Mass_Ratio_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_RATIO_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassRatioDistributionMin
   * - Description:
     - Value of program option ``--mass-ratio-min``
   * - Header String:
     - Mass_Ratio_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MAXIMUM_EVOLUTION_TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MaxEvolutionTime
   * - Description:
     - Value of program option ``--maximum-evolution-time``
   * - Header String:
     - Max_Evolution_Time

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MAXIMUM_DONOR_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MaximumMassDonorNandezIvanova
   * - Description:
     - Value of program option ``--maximum-mass-donor-nandez-ivanova``
   * - Header String:
     - Max_Donor_Mass

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MAXIMUM_NEUTRON_STAR_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MaximumNeutronStarMass
   * - Description:
     - Value of program option ``--maximum-neutron-star-mass``
   * - Header String:
     - Max_NS_Mass

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MAXIMUM_TIMESTEPS**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MaxNumberOfTimestepIterations
   * - Description:
     - Value of program option ``--maximum-number-timestep-iterations``
   * - Header String:
     - Max_Timesteps

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MCBUR1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_mCBUR1
   * - Description:
     - Value of program option ``--mcbur1``
   * - Header String:
     - MCBUR1

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_Metallicity
   * - Description:
     - Value of program option ``--metallicity``
   * - Header String:
     - Metallicity

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MetallicityDistribution
   * - Description:
     - Value of program option ``--metallicity-distribution``
   * - Header String:
     - Metallicity_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MetallicityDistributionMax
   * - Description:
     - Value of program option ``--metallicity-max``
   * - Header String:
     - Metallicity_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MetallicityDistributionMin
   * - Description:
     - Value of program option ``--metallicity-min``
   * - Header String:
     - Metallicity_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MINIMUM_MASS_SECONDARY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MinimumMassSecondary
   * - Description:
     - Value of program option ``--minimum-secondary-mass``
   * - Header String:
     - Min_Secondary_Mass

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_ACCRETION_EFFICIENCY_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassTransferAccretionEfficiencyPrescription
   * - Description:
     - Value of program option ``--mass-transfer-accretion-efficiency-prescription``
   * - Header String:
     - MT_Acc_Efficiency_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_ANG_MOM_LOSS_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassTransferAngularMomentumLossPrescription
   * - Description:
     - Value of program option ``--mass-transfer-angular-momentum-loss-prescription``
   * - Header String:
     - MT_AngMom_Loss_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_FRACTION_ACCRETED**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassTransferFractionAccreted
   * - Description:
     - Value of program option ``--mass-transfer-fa``
   * - Header String:
     - MT_Fraction_Accreted

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_JLOSS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassTransferJloss
   * - Description:
     - Value of program option ``--mass-transfer-jloss``
   * - Header String:
     - MT_JLoss

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_THERMAL_LIMIT_C**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MassTransferCParameter
   * - Description:
     - Value of program option ``--mass-transfer-thermal-limit-C``
   * - Header String:
     - MT_Thermal_Limit_C

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_REJUVENATION_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassTransferRejuvenationPrescription
   * - Description:
     - Value of program option ``--mass-transfer-rejuvenation-prescription``
   * - Header String:
     - MT_Rejuvenation_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_THERMALLY_LIMITED_VARIATION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_MassTransferThermallyLimitedVariation
   * - Description:
     - Value of program option ``--mass-transfer-thermal-limit-accretor``
   * - Header String:
     - MT_Thermally_Lmtd_Variation

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MULLER_MANDEL_KICK_MULTIPLIER_BH**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MullerMandelKickBH
   * - Description:
     - Value of program option ``--muller-mandel-kick-multiplier-BH``
   * - Header String:
     - MM_Kick_Multiplier_BH

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MULLER_MANDEL_KICK_MULTIPLIER_NS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_MullerMandelKickNS
   * - Description:
     - Value of program option ``--muller-mandel-kick-multiplier-NS``
   * - Header String:
     - MM_Kick_Multiplier_NS

.. _spec-options-props-N:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NEUTRINO_MASS_LOSS_ASSUMPTION_BH**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_NeutrinoMassLossAssumptionBH
   * - Description:
     - Value of program option ``--neutrino-mass-loss-BH-formation``
   * - Header String:
     - Neutrino_Mass_Loss_Assmptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NEUTRINO_MASS_LOSS_VALUE_BH**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_NeutrinoMassLossValueBH
   * - Description:
     - Value of program option ``--neutrino-mass-loss-BH-formation-value``
   * - Header String:
     - Neutrino_Mass_Loss_Value

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NOTES**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - Options::m_Notes
   * - Description:
     - Value of program option ``--notes``
   * - Header String:
     - as specified by program option ``--Notes-Hdrs``

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NS_EOS**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_NeutronStarEquationOfState
   * - Description:
     - Value of program option ``--neutron-star-equation-of-state``
   * - Header String:
     - NS_EOS

.. _spec-options-props-O:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_PERIOD**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_OrbitalPeriod
   * - Description:
     - Value of program option ``--orbital-period``
   * - Header String:
     - Orbital_Period

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_PERIOD_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_OrbitalPeriodDistribution
   * - Description:
     - Value of program option ``--orbital-period-distribution``
   * - Header String:
     - Orbital_Period_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_PERIOD_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_OrbitalPeriodDistributionMax
   * - Description:
     - Value of program option ``--orbital-period-max``
   * - Header String:
     - Orbital_Period_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_PERIOD_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_OrbitalPeriodDistributionMin
   * - Description:
     - Value of program option ``--orbital-period-min``
   * - Header String:
     - Orbital_Period_Min

.. _spec-options-props-P:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PISN_LOWER_LIMIT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PairInstabilityLowerLimit
   * - Description:
     - Value of program option ``--PISN-lower-limit``
   * - Header String:
     - PISN_Lower_Limit

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PISN_UPPER_LIMIT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PairInstabilityUpperLimit
   * - Description:
     - Value of program option ``--PISN-upper-limit``
   * - Header String:
     - PISN_Upper_Limit

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PPI_LOWER_LIMIT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsationalPairInstabilityLowerLimit
   * - Description:
     - Value of program option ``--PPI-lower-limit``
   * - Header String:
     - PPI_Lower_Limit

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PPI_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_PulsationalPairInstabilityPrescription
   * - Description:
     - Value of program option ``--pulsational-pair-instability-prescription``
   * - Header String:
     - PPI_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PPI_UPPER_LIMIT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsationalPairInstabilityUpperLimit
   * - Description:
     - Value of program option ``--PPI-upper-limit``
   * - Header String:
     - PPI_Upper_Limit

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_PulsarBirthMagneticFieldDistribution
   * - Description:
     - Value of program option ``--pulsar-birth-magnetic-field-distribution``
   * - Header String:
     - Pulsar_Mag_Field_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarBirthMagneticFieldDistributionMax
   * - Description:
     - Value of program option ``--pulsar-birth-magnetic-field-distribution-max``
   * - Header String:
     - Pulsar_Mag_Field_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarBirthMagneticFieldDistributionMin
   * - Description:
     - Value of program option ``--pulsar-birth-magnetic-field-distribution-min``
   * - Header String:
     - Pulsar_Mag_Field_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_PulsarBirthSpinPeriodDistribution
   * - Description:
     - Value of program option ``--pulsar-birth-spin-period-distribution``
   * - Header String:
     - Pulsar_Spin_Period_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_BIRTH_SPIN_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarBirthSpinPeriodDistributionMax
   * - Description:
     - Value of program option ``--pulsar-birth-spin-period-distribution-max``
   * - Header String:
     - Pulsar_Spin_Period_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_BIRTH_SPIN_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarBirthSpinPeriodDistributionMin
   * - Description:
     - Value of program option ``--pulsar-birth-spin-period-distribution-min``
   * - Header String:
     - Pulsar_Spin_Period_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DECAY_MASS_SCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options:m_PulsarMagneticFieldDecayMassscale
   * - Description:
     - Value of program option ``--pulsar-magnetic-field-decay-massscale``
   * - Header String:
     - Pulsar_Mag_Field_Decay_mScale

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD_DECAY_TIME_SCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options:m_PulsarMagneticFieldDecayTimescale
   * - Description:
     - Value of program option ``--pulsar-magnetic-field-decay-timescale``
   * - Header String:
     - Pulsar_Mag_Field_Decay_tScale

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MINIMUM_MAGNETIC_FIELD**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_PulsarLog10MinimumMagneticField
   * - Description:
     - Value of program option ``--pulsar-minimum-magnetic-field``
   * - Header String:
     - Pulsar_Minimum_Mag_Field

.. _spec-options-props-Q:

.. _spec-options-props-R:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RANDOM_SEED**
     -
   * - Data type:
     - UNSIGNED LONG INT
   * - COMPAS variable:
     - Options::m_RandomSeed
   * - Description:
     - Value of program option ``--random-seed``
   * - Header String:
     - SEED(OPTION)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RANDOM_SEED_CMDLINE**
     -
   * - Data type:
     - UNSIGNED LONG INT
   * - COMPAS variable:
     - Options::m_FixedRandomSeed
   * - Description:
     - Value of program option ``--random-seed`` (specified on the commandline)
   * - Header String:
     - SEED(CMDLINE)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **REMNANT_MASS_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_RemnantMassPrescription
   * - Description:
     - Value of program option ``--remnant-mass-prescription``
   * - Header String:
     - Remnant_Mass_Prscrptn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROTATIONAL_VELOCITY_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_RotationalVelocityDistribution
   * - Description:
     - Value of program option ``--rotational-velocity-distribution``
   * - Header String:
     - Rotational_Velocity_Dstrbtn

.. _spec-options-props-S:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_SemiMajorAxis
   * - Description:
     - Value of program option ``--semi-major-axis``
   * - Header String:
     - Semi-Major_Axis

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_DISTRIBUTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_SemiMajorAxisDistribution
   * - Description:
     - Value of program option ``--semi-major-axis-distribution``
   * - Header String:
     - Semi-Major_Axis_Dstrbtn

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_DISTRIBUTION_MAX**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_SemiMajorAxisDistributionMax
   * - Description:
     - Value of program option ``--semi-major-axis-max``
   * - Header String:
     - Semi-Major_Axis_Dstrbtn_Max

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_DISTRIBUTION_MIN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_SemiMajorAxisDistributionMin
   * - Description:
     - Value of program option ``--semi-major-axis-min``
   * - Header String:
     - Semi-Major_Axis_Dstrbtn_Min

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_ZETA_PRESCRIPTION**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - Options::m_StellarZetaPrescription
   * - Description:
     - Value of program option ``--stellar-zeta-prescription``
   * - Header String:
     - Stellar_Zeta_Prscrptn

.. _spec-options-props-T:

.. _spec-options-props-U:

.. _spec-options-props-V:

.. _spec-options-props-W:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **WR_FACTOR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_WolfRayetFactor
   * - Description:
     - Value of program option ``--wolf-rayet-multiplier``
   * - Header String:
     - WR_Factor

.. _spec-options-props-X:

.. _spec-options-props-Y:

.. _spec-options-props-Z:

:ref:`Back to Top <spec-options-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_ADIABATIC_ARBITRARY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_ZetaAdiabaticArbitrary
   * - Description:
     - Value of program option ``--zeta-adiabatic-arbitrary``
   * - Header String:
     - Zeta_Adiabatic_Arbitrary

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_MS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_ZetaMainSequence
   * - Description:
     - Value of program option ``--zeta-main-sequence``
   * - Header String:
     - Zeta_Main_Sequence_Giant

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_RADIATIVE_ENVELOPE_GIANT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - Options::m_ZetaRadiativeEnvelopeGiant
   * - Description:
     - Value of program option ``--zeta-radiative-envelope-giant``
   * - Header String:
     - Zeta_Radiative_Envelope_Giant

SSE supernovae
==============

Default record definition for the SSE Supernovae log file::

    const ANY_PROPERTY_VECTOR SSE_SUPERNOVAE_REC = {
        STAR_PROPERTY::RANDOM_SEED,
        STAR_PROPERTY::DRAWN_KICK_MAGNITUDE,
        STAR_PROPERTY::KICK_MAGNITUDE,
        STAR_PROPERTY::FALLBACK_FRACTION,
        STAR_PROPERTY::MEAN_ANOMALY,				
        STAR_PROPERTY::SN_TYPE,
        STAR_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION,
        STAR_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
        STAR_PROPERTY::MASS,
        STAR_PROPERTY::STELLAR_TYPE,
        STAR_PROPERTY::STELLAR_TYPE_PREV,
        STAR_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
        STAR_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
        STAR_PROPERTY::TIME,
        STAR_PROPERTY::IS_HYDROGEN_POOR
    };

Standard log files
==================

COMPAS defines several standard log files that may be produced depending upon the simulation mode (Single Star Evolution (SSE), 
or Binary Star Evolution (BSE), see the ``--mode`` program option)), and the value of various program options.

The standard log files are:

    .. list-table::
       :widths: 22 78 
       :header-rows: 0
       :class: aligned-text

       * - System Parameters
         - Records summary information for all stars, or binary stars, during evolution.
       * -
         -
       * - Supernovae
         - Records summary information for all stars that experience a SN event during evolution.
       * -
         -
       * - Detailed Output
         - Records detailed information for a star, or binary star, during evolution.
       * -
         - Enable with program option ``--detailed-output``.
       * -
         -
       * - SwitchLog
         - Records detailed information for all stars, or binary stars, at the time of each stellar type switch during evolution.
       * - 
         - Enable with program option ``--switch-log``.
       * -
         -
       * - Double Compact Objects
         - Records summary information for all binary systems that form DCOs during BSE.
       * -
         -
       * - Common Envelopes
         - Records summary information for all binary systems that experience CEEs during BSE.
       * -
         -
       * - Pulsar Evolution
         - Records detailed Pulsar evolution information during BSE.
       * -
         - Enable with program option ``--evolve-pulsars``.
       * -
         -
       * - RLOF
         - Records detailed information RLOF events during BSE.
       * - 
         - Enable with program option ``--rlof-printing``.

BSE pulsar evolution
====================

Default record definition for the BSE Pulsar Evolution log file::

    const ANY_PROPERTY_VECTOR BSE_PULSAR_EVOLUTION_REC = {
        BINARY_PROPERTY::RANDOM_SEED,
        STAR_1_PROPERTY::MASS,
        STAR_2_PROPERTY::MASS,
        STAR_1_PROPERTY::STELLAR_TYPE,
        STAR_2_PROPERTY::STELLAR_TYPE,
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL,
        BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY,
        STAR_1_PROPERTY::PULSAR_MAGNETIC_FIELD,
        STAR_2_PROPERTY::PULSAR_MAGNETIC_FIELD,
        STAR_1_PROPERTY::PULSAR_SPIN_FREQUENCY,
        STAR_2_PROPERTY::PULSAR_SPIN_FREQUENCY,
        STAR_1_PROPERTY::PULSAR_SPIN_DOWN_RATE,
        STAR_2_PROPERTY::PULSAR_SPIN_DOWN_RATE,
        BINARY_PROPERTY::TIME,
        BINARY_PROPERTY::DT
    };

    Standard log file format
========================

COMPAS can produce log files in several formats:

    - Hierarchical Data Format version 5 (``HDF5``)\ [#f1]_
    - Comma Separated Values (``CSV``)
    - Tab Separated Values (``TSV``)
    - Plain text: space separated values (``TXT``)
    
The log file type is set using the ``--logfile-type`` program option.

Standard ``CSV``, ``TSV``, and ``TXT`` log files are human-readable files, and formatted in a similar fashion. Each standard
``CSV``, ``TSV``, and ``TXT`` log file consists of three header records followed by data records. Header records and data records
are delimiter separated fields, and the fields as specified by the log file record specifier.

The header records for all standard ``CSV``, ``TSV``, and ``TXT`` log files are:

    Header record 1: Column Data Type Names |br|
    Header record 2: Column Units (where applicable) |br|
    Header record 3: Column Headings

Column Data Type Names are taken from the set **{ BOOL, INT, FLOAT, STRING }**, where

    .. list-table::
       :widths: 12 88 
       :header-rows: 0
       :class: aligned-text

       * - **BOOL**
         - indicates the data value will be a boolean value.
       * - 
         - Boolean data values will be recorded in the log file in either numerical format (1 or 0, where 1 = TRUE and 0 = FALSE), or string format ("TRUE" or "FALSE"), depending upon the value of the ``--print-bool-as-string`` program option.
       * - **INT**
         - indicates the data value will be an integer number.
       * - **FLOAT**
         - indicates the data value will be a floating-point number.
       * - **STRING**
         - indicates the data value will be a text string.

Column Units is a string indicating the units of the corresponding data values (e.g. "Msol\*AU\ :sup:`2`\ \*yr\ :sup:`-1`\ ",
"Msol", "AU", etc.). The Column Units value may be blank where units are not applicable, or may be one of:

    .. list-table::
       :widths: 12 88 
       :header-rows: 0
       :class: aligned-text

       * - **Count**
         - indicates the data value is the total of a counted entity.
       * - **State**
         - indicates the data value describes a state (e.g. "Unbound" state is "TRUE" or "FALSE").
       * - **Event**
         - the data value describes an event status (e.g. "Simultaneous_RLOF" is "TRUE").

Column Headings are string labels that describe the corresponding data values. The heading strings for stellar properties of
constituent stars of a binary will have appropriate identifiers appended. That is, heading strings for:

    .. list-table::
       :widths: 38 62 
       :header-rows: 0
       :class: aligned-text

       * - STAR_1 PROPERTY::properties
         - will have ":boldtext:`(1)`" appended
       * - STAR_2 PROPERTY::properties
         - will have ":boldtext:`(2)`" appended
       * - SUPERNOVA_PROPERTY::properties
         - will have ":boldtext:`(SN)`" appended: any column with a header with a suffix of ":boldtext:`(SN)`" represents an attribute of the star undergoing a supernova event, either just before the supernova [e.g., `Mass_Total@CO(SN)`] or just after the supernovae [e.g., `Mass(SN)`].
       * - COMPANION_PROPERTY::properties
         - will have ":boldtext:`(CP)`" appended: any column with a header with a suffix of ":boldtext:`(CP)`" represents an attribute of the the companion after the supernova event.

``HDF5`` files are not human-readable. The ``HDF5`` file format supports large, complex, heterogeneous data, enabling the data to be stored
in a structured way in a single file. When the ``HDF5`` format is specified for COMPAS log files, a single ``HDF5`` file is produced for 
non-detailed output log files, containing all non-detailed output log files described above. Detailed output files are created, as for
other logfile types, as individual files (in this case, ``HDF5`` files), in the ’Detailed_Output’ container directory.

Each file described above is created as a `group` within the ``HDF5`` file, with the name of the group set to the name of the file
(e.g. "BSE_System_Parameters"). Each column in the files described above is created as a `dataset` within its corresponding group in the
``HDF5`` file, with the name of the datset set to the column header as described above (e.g. "Mass(1)"). Each dataset in an ``HDF5`` file
is typed, and the dataset data types are set to the column data types as described above. The column units described above are attached to
their corresponding datasets in the ``HDF5`` file as `attributes`.


.. rubric:: Footnotes

.. [#f1] https://www.hdfgroup.org/
Binary properties
=================

When specifying known properties in a log file record specification record, the property name must be prefixed with 
the property type. Currently there is a single binary property type available for use: BINARY_PROPERTY.

For example, to specify the property ``SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE`` for a binary star being evolved in ``BSE``, use::

    BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE


.. _binary-props-top:

Following is an alphabetical list of binary properties available for inclusion in log file record specifiers.

**Jump to**
:ref:`A <binary-props-A>` :ref:`B <binary-props-B>` :ref:`C <binary-props-C>` :ref:`D <binary-props-D>`
:ref:`E <binary-props-E>` :ref:`F <binary-props-F>` :ref:`G <binary-props-G>` :ref:`H <binary-props-H>`
:ref:`I <binary-props-I>` :ref:`J <binary-props-J>` :ref:`K <binary-props-K>` :ref:`L <binary-props-L>`
:ref:`M <binary-props-M>` :ref:`N <binary-props-N>` :ref:`O <binary-props-O>` :ref:`P <binary-props-P>`
:ref:`Q <binary-props-Q>` :ref:`R <binary-props-R>` :ref:`S <binary-props-S>` :ref:`T <binary-props-T>`
:ref:`U <binary-props-U>` :ref:`V <binary-props-V>` :ref:`W <binary-props-W>` :ref:`X <binary-props-X>`
:ref:`Y <binary-props-Y>` :ref:`Z <binary-props-Z>`


Following is the list of binary properties available for inclusion in log file record specifiers.
Binary Properties


.. _binary-props-A:

.. _binary-props-B:

.. _binary-props-C:

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CIRCULARIZATION_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CircularizationTimescale
   * - Description:
     - Tidal circularisation timescale (Myr)
   * - Header String:
     - Tau_Circ

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_AT_LEAST_ONCE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_CEDetails.CEEcount
   * - Description:
     - Flag to indicate if there has been at least one common envelope event.
   * - Header String:
     - CEE    

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **COMMON_ENVELOPE_EVENT_COUNT**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.CEEcount
   * - Description:
     - The number of common envelope events.
   * - Header String:
     - CE_Event_Counter

.. _binary-props-D:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DIMENSIONLESS_KICK_MAGNITUDE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_uK
   * - Description:
     - Dimensionless kick magnitude supplied by user (see option --fix-dimensionless-kick-magnitude).
   * - Header String:
     - Kick_Magnitude(uK)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DOUBLE_CORE_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.doubleCoreCE
   * - Description:
     - Flag to indicate double-core common envelope.
   * - Header String:
     - Double_Core_CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_Dt
   * - Description:
     - Current timestep (Myr).
   * - Header String:
     - dT

.. _binary-props-E:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_Eccentricity
   * - Description:
     - Orbital eccentricity.
   * - Header String:
     - Eccentricity

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_AT_DCO_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_EccentricityAtDCOFormation
   * - Description:
     - Orbital eccentricity at DCO formation.
   * - Header String:
     - Eccentricity@\ DCO

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_INITIAL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_EccentricityInitial
   * - Description:
     - Supplied by user via grid file or sampled from distribution (see ``--eccentricity-distribution`` option) (default).
   * - Header String:
     - Eccentricity@\ ZAMS

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.postCEE.eccentricity
   * - Description:
     - Eccentricity immediately following common envelope event.
   * - Header String:
     - Eccentricity>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_PRE_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_EccentricityPreSN
   * - Description:
     - Eccentricity of the binary immediately prior to supernova event.
   * - Header String:
     - Eccentricity<SN

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRICITY_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.preCEE.eccentricity
   * - Description:
     - Eccentricity at the onset of RLOF leading to the CE.
   * - Header String:
     - Eccentricity<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ERROR**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_Error
   * - Description:
     - Error number (if error condition exists, else 0).
   * - Header String:
     - Error

.. _binary-props-F:

.. _binary-props-G:

.. _binary-props-H:

.. _binary-props-I:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ID**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_ObjectId
   * - Description:
     - Unique object identifier for ``C++`` object – used in debugging to identify objects.
   * - Header String:
     - ID

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IMMEDIATE_RLOF_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.immediateRLOFPostCEE
   * - Description:
     - Flag to indicate if either star overflows its Roche lobe immediately following common envelope event.
   * - Header String:
     - Immediate_RLOF>CE

.. _binary-props-J:

.. _binary-props-K:

.. _binary-props-L:

.. _binary-props-M:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.mass
   * - Description:
     - Mass of the primary star immediately following common envelope event (\ :math:`M_\odot`).
   * - Header String:
     - Mass(1)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.mass
   * - Description:
     - Mass of the primary star immediately prior to common envelope event (\ :math:`M_\odot`).
   * - Header String:
     - Mass(1)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.mass
   * - Description:
     - Mass of the secondary star immediately following common envelope event (\ :math:`M_\odot`).
   * - Header String:
     - Mass(2)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.mass
   * - Description:
     - Mass of the secondary star immediately prior to common envelope event (\ :math:`M_\odot`).
   * - Header String:
     - Mass(2)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_ENV_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_MassEnv1
   * - Description:
     - Envelope mass of the primary star (\ :math:`M_\odot`).
   * - Header String:
     - Mass_Env(1)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_ENV_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_MassEnv2
   * - Description:
     - Envelope mass of the secondary star (\ :math:`M_\odot`).
   * - Header String:
     - Mass_Env(2)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASSES_EQUILIBRATED**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_MassesEquilibrated
   * - Description:
     - Flag to indicate whether chemically homogeneous stars had masses equilibrated and orbit circularised due to Roche lobe overflow during evolution.
   * - Header String:
     - Equilibrated

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASSES_EQUILIBRATED_AT_BIRTH**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_MassesEquilibratedAtBirth
   * - Description:
     - Flag to indicate whether stars had masses equilibrated and orbit circularised at birth due to Roche lobe overflow.
   * - Header String:
     - Equilibrated_At_Birth

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_TRANSFER_TRACKER_HISTORY**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_MassTransferTrackerHistory
   * - Description:
     - Indicator of mass transfer history for the binary. Will be printed as one of:

        .. list-table::
           :widths: 35 5
           :header-rows: 0
           :class: aligned-text
           
           * - NO MASS TRANSFER
             - = 0
           * - STABLE FROM 1 TO 2
             - = 1
           * - STABLE FROM 2 TO 1
             - = 2
           * - CE FROM 1 TO 2
             - = 3
           * - CE FROM 2 TO 1
             - = 4
           * - CE DOUBLE CORE
             - = 5
           * - CE BOTH MS
             - = 6
           * - CE MS WITH CO
             - = 7

   * - Header String:
     - MT_History

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MERGES_IN_HUBBLE_TIME**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_MergesInHubbleTime
   * - Description:
     - Flag to indicate if the binary compact remnants merge within a Hubble time.
   * - Header String:
     - Merges_Hubble_Time

.. _binary-props-N:

.. _binary-props-O:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **OPTIMISTIC_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.optimisticCE
   * - Description:
     - Flag that returns TRUE if we have a Hertzsprung-gap star, and we allow it to survive the CE.
   * - Header String:
     - Optimistic_CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_VELOCITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_OrbitalVelocity
   * - Description:
     - Orbital velocity (\ :math:`km s^{-1}`).
   * - Header String:
     - Orbital_Velocity

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_VELOCITY_PRE_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_OrbitalVelocityPreSN
   * - Description:
     - Orbital velocity immediately prior to supernova event (\ :math:`km s^{-1}`).
   * - Header String:
     - Orbital_Velocity<SN

.. _binary-props-P:

.. _binary-props-Q:

.. _binary-props-R:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.radius
   * - Description:
     - Radius of the primary star immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - Radius(1)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.radius
   * - Description:
     - Radius of the primary star at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - Radius(1)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.radius
   * - Description:
     - Radius of the secondary star immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - Radius(2)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.radius
   * - Description:
     - Radius of the secondary star at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - Radius(2)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RANDOM_SEED**
     -
   * - Data type:
     - UNSIGNED LONG
   * - COMPAS variable:
     - BaseBinaryStar::m_RandomSeed
   * - Description:
     - Seed for random number generator for this binary star. Optionally supplied by user via program option ``--random-seed``; default generated from system time.
   * - Header String:
     - SEED

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→isCE
   * - Description:
     - Flag to indicate if the RLOF leads to a common-envelope event 
   * - Header String:
     - CEE>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_ECCENTRICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→eccentricity
   * - Description:
     - Eccentricity immediately after RLOF.
   * - Header String:
     - Eccentricity>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_EVENT_COUNTER**
     -
   * - Data type:
     - UNSIGNED INT
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→eventCounter
   * - Description:
     - The number of times the binary has overflowed its Roche lobe up to and including this episode
   * - Header String:
     - MT_Event_Counter

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_SEMI_MAJOR_AXIS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→semiMajorAxis
   * - Description:
     - Semi-major Axis(\ :math:`R_\odot`) immediately after RLOF.
   * - Header String:
     - SemiMajorAxis>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→mass1
   * - Description:
     - Mass (\ :math:`M_\odot`) of the primary immediately after RLOF.
   * - Header String:
     - Mass(1)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→mass2
   * - Description:
     - Mass (\ :math:`M_\odot`) of the secondary immediately after RLOF.
   * - Header String:
     - Mass(2)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→radius1
   * - Description:
     - Radius (\ :math:`R_\odot`) of the primary immediately after RLOF.
   * - Header String:
     - Radius(1)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→radius2
   * - Description:
     - Radius (\ :math:`R_\odot`) of the secondary immediately after RLOF.
   * - Header String:
     - Radius(2)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→isRLOF1
   * - Description:
     - Flag to indicate whether the primary is overflowing its Roche Lobe.
   * - Header String:
     - RLOF(1)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→isRLOF2
   * - Description:
     - Flag to indicate whether the secondary is overflowing its Roche Lobe.
   * - Header String:
     - RLOF(2)>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m RLOFDetails.propsPostMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star immediately after RLOF.
   * - Header String:
     - Stellar_Type(1)>MT

`Note that this property has the same header string as RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR1_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m RLOFDetails.propsPostMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star immediately after RLOF.
   * - Header String:
     - Stellar_Type(1)>MT

`Note that this property has the same header string as RLOF_POST_MT_STAR1_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star immediately after RLOF.
   * - Header String:
     - Stellar_Type(2)>MT

`Note that this property has the same header string as RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_STAR2_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_RLOFDetails.propsPostMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star immediately after RLOF.
   * - Header String:
     - Stellar_Type(2)>MT

`Note that this property has the same header string as RLOF_POST_MT_STAR2_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_POST_MT_TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPostMT→time
   * - Description:
     - Time since ZAMS (Myr) immediately after RLOF.
   * - Header String:
     - Time>MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_ECCENTRICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→eccentricity
   * - Description:
     - Eccentricity at the onset of RLOF.
   * - Header String:
     - Eccentricity<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_SEMI_MAJOR_AXIS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→semiMajorAxis
   * - Description:
     - Semi-major Axis (\ :math:`R_\odot`) at the onset of RLOF.
   * - Header String:
     - SemiMajorAxis<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→mass1
   * - Description:
     - Mass (\ :math:`M_\odot`) of the primary at the onset of RLOF.
   * - Header String:
     - Mass(1)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→mass2
   * - Description:
     - Mass (\ :math:`M_\odot`) of the secondary at the onset of RLOF.
   * - Header String:
     - Mass(2)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→radius1
   * - Description:
     - Radius (\ :math:`R_\odot`) of the primary at the onset of RLOF.
   * - Header String:
     - Radius(1)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→radius2
   * - Description:
     - Radius (\ :math:`R_\odot`) of the secondary at the onset of RLOF.
   * - Header String:
     - Radius(2)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→isRLOF1
   * - Description:
     - Flag to indicate whether the primary is overflowing its Roche Lobe.
   * - Header String:
     - RLOF(1)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→isRLOF2
   * - Description:
     - Flag to indicate whether the secondary is overflowing its Roche Lobe.
   * - Header String:
     - RLOF(2)<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star at the onset of RLOF.
   * - Header String:
     - Stellar_Type(1)<MT

`Note that this property has the same header string as RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR1_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_RLOFDetails.propsPreMT→stellarType1
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star at the onset of RLOF.
   * - Header String:
     - Stellar_Type(1)<MT

`Note that this property has the same header string as RLOF_PRE_MT_STAR1_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→stellarType2
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star at the onset of RLOF.
   * - Header String:
     - Stellar_Type(2)<MT

`Note that this property has the same header string as RLOF_PRE_MTvSTAR2_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_STAR2_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_RLOFDetails.propsPreMT→stellarType2
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star at the onset of RLOF.
   * - Header String:
     - Stellar_Type(2)<MT

`Note that this property has the same header string as RLOF_PRE_MT_STAR2_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_PRE_MT_TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.propsPreMT→time
   * - Description:
     - Time since ZAMS (Myr) at the onset of RLOF.
   * - Header String:
     - Time<MT

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_RocheLobeRadius
   * - Description:
     - Roche radius of the primary star (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(1)|a

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_RocheLobeRadius
   * - Description:
     - Roche radius of the secondary star (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(2)|a

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.postCEE.rocheLobe1to2
   * - Description:
     - Roche radius of the primary star immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(1)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.postCEE.rocheLobe2to1
   * - Description:
     - Roche radius of the secondary star immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(2)>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.preCEE.rocheLobe1to2
   * - Description:
     - Roche radius of the primary star at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(1)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.preCEE.rocheLobe2to1
   * - Description:
     - Roche radius of the secondary star at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - RocheLobe(2)<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_TRACKER_1**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_RocheLobeTracker
   * - Description:
     - Ratio of the primary star’s stellar radius to Roche radius (R/RL), evaluated at periapsis.
   * - Header String:
     - Radius|RL

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ROCHE_LOBE_TRACKER_2**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_RocheLobeTracker
   * - Description:
     - Ratio of the secondary star’s stellar radius to Roche radius (R/RL), evaluated at periapsis.
   * - Header String:
     - Radius|RL

.. _binary-props-S:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_AT_DCO_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SemiMajorAxisAtDCOFormation
   * - Description:
     - Semi-major axis at DCO formation (AU).
   * - Header String:
     - SemiMajorAxis@\ DCO

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_INITIAL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SemiMajorAxisInitial
   * - Description:
     - Semi-major axis at ZAMS (AU).
   * - Header String:
     - SemiMajorAxis@\ ZAMS

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.postCEE.semiMajorAxis
   * - Description:
     - Semi-major axis immediately following common envelope event (\ :math:`R_\odot`).
   * - Header String:
     - SemiMajorAxis>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_PRE_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SemiMajorAxisPreSN
   * - Description:
     - Semi-major axis immediately prior to supernova event (AU).
   * - Header String:
     - SemiMajorAxis<SN

`Note that this property has the same header string as SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_SemiMajorAxisPreSN
   * - Description:
     - Semi-major axis immediately prior to supernova event (\ :math:`R_\odot`).
   * - Header String:
     - SemiMajorAxis<SN

`Note that this property has the same header string as SEMI_MAJOR_AXIS_PRE_SUPERNOVA. It is expected that one or the other is printed in any file, but 
not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_CEDetails.preCEE.semiMajorAxis
   * - Description:
     - Semi-major axis at the onset of RLOF leading to the common-envelope episode (\ :math:`R_\odot`).
   * - Header String:
     - SemiMajorAxis<CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SemiMajorAxis
   * - Description:
     - Semi-major axis at ZAMS (AU).
   * - Header String:
     - SemiMajorAxis@\ ZAMS

`Note that this property has the same header string as SEMI_MAJOR_AXIS_RSOL. It is expected that one or the other is printed in any file, but not both. 
If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SEMI_MAJOR_AXIS_RSOL**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_SemiMajorAxis
   * - Description:
     - Semi-major axis at ZAMS (\ :math:`R_\odot`).
   * - Header String:
     - SemiMajorAxis@\ ZAMS

`Note that this property has the same header string as SEMI_MAJOR_AXIS. It is expected that one or the other is printed in any file, but not both. If both 
are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SIMULTANEOUS_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.simultaneousRLOF
   * - Description:
     - Flag to indicate that both stars are undergoing RLOF.
   * - Header String:
     - Simultaneous_RLOF

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STABLE_RLOF_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_RLOFDetails.stableRLOFPostCEE
   * - Description:
     - Flag to indicate stable mass transfer after common envelope event.
   * - Header String:
     - Stable_RLOF>CE

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_MERGER**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_StellarMerger
   * - Description:
     - Flag to indicate the stars merged (were touching) during evolution.
   * - Header String:
     - Merger

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_MERGER_AT_BIRTH**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_StellarMergerAtBirth
   * - Description:
     - Flag to indicate the stars merged (were touching) at birth.
   * - Header String:
     - Merger_At_Birth

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.stellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star immediately following common envelope event.
   * - Header String:
     - Stellar_Type(1)>CE

`Note that this property has the same header string as STELLAR_TYPE_NAME_1_POST_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.stellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the primary star at the onset of RLOF leading to the common-envelope episode.
   * - Header String:
     - Stellar_Type(1)<CE

`Note that this property has the same header string as STELLAR_TYPE_NAME_1_PRE_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.stellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star immediately following common envelope event.
   * - Header String:
     - Stellar_Type(2)>CE

`Note that this property has the same header string as STELLAR_TYPE_NAME_2_POST_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.stellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) of the secondary star at the onset of RLOF leading to the common-envelope episode.
   * - Header String:
     - Stellar_Type(2)<CE

`Note that this property has the same header string as STELLAR_TYPE_NAME_2_PRE_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME_1_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_CEDetails.postCEE.stellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) of the primary star immediately following common envelope event. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header String:
     - Stellar_Type(1)>CE

`Note that this property has the same header string as STELLAR_TYPE_1_POST_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME_1_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_CEDetails.preCEE.stellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) of the primary star at the onset of RLOF leading to the common-envelope episode. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header String:
     - Stellar_Type(1)<CE

`Note that this property has the same header string as STELLAR_TYPE_1_PRE_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, but not 
both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME_2_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_CEDetails.postCEE.stellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) of the secondary star immediately following common envelope event. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header String:
     - Stellar_Type(2)>CE

`Note that this property has the same header string as STELLAR_TYPE_2_POST_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, but not 
both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME_2_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_CEDetails.preCEE.stellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) of the secondary star at the onset of RLOF leading to the common-envelope episode. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header String:
     - Stellar_Type(2)<CE

`Note that this property has the same header string as STELLAR_TYPE_2_PRE_COMMON_ENVELOPE. It is expected that one or the other is printed in any file, but not 
both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SUPERNOVA_STATE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseBinaryStar::m_SupernovaState
   * - Description:
     - Indicates which star(s) went supernova. Will be printed as one of:

        .. list-table::
           :widths: 45 5
           :header-rows: 0
           :class: aligned-text

           * - No supernova
             - = 0
           * - Star 1 is the supernova
             - = 1
           * - Star 2 is the supernova
             - = 2
           * - Both stars are supernovae
             - = 3

   * - Header String:
     - Supernova_State

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SYNCHRONIZATION_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SynchronizationTimescale
   * - Description:
     - Tidal synchronisation timescale (Myr).
   * - Header String:
     - Tau_Sync

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SYSTEMIC_SPEED**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_SystemicVelocity
   * - Description:
     - Post-supernova systemic (centre-of-mass) velocity (\ :math:`km s^{-1}`).
   * - Header String:
     - Systemic_Velocity

.. _binary-props-T:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_Time
   * - Description:
     - Time since ZAMS (Myr).
   * - Header String:
     - Time

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TIME_TO_COALESCENCE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_TimeToCoalescence
   * - Description:
     - Time between formation of double compact object and gravitational-wave driven merger (Myr).
   * - Header String:
     - Coalescence_Time

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TOTAL_ANGULAR_MOMENTUM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_TotalAngularMomentum
   * - Description:
     - Total angular momentum calculated using regular conservation of energy (\ :math:`M_\odot  AU^2 y^{r−1}`).
   * - Header String:
     - Ang_Momentum_Total

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TOTAL_ENERGY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_TotalEnergy
   * - Description:
     - Total energy calculated using regular conservation of energy (\ :math:`M_\odot  AU^2 y^{r−1}`).
   * - Header String:
     - Energy_Total

.. _binary-props-U:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **UNBOUND**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseBinaryStar::m_Unbound
   * - Description:
     - Flag to indicate the binary is unbound (or has become unbound after a supernova event).
   * - Header String:
     - Unbound

.. _binary-props-V:

.. _binary-props-W:

.. _binary-props-X:

.. _binary-props-Y:

.. _binary-props-Z:

:ref:`Back to Top <binary-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_LOBE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_ZetaLobe
   * - Description:
     - The logarithmic derivative of Roche lobe radius with respect to donor mass for :math:`q=\frac{Md}{Ma}` at the onset of the RLOF.
   * - Header String:
     - Zeta_Lobe

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_STAR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseBinaryStar::m_ZetaStar
   * - Description:
     - Mass-radius exponent of the star at the onset of the RLOF. Calculated differently based on the value of program option ``--zeta-prescription``
   * - Header String:
     - Zeta_Star

Standard log file record specifiers
===================================

Each standard log file has an associated log file record specifier that defines what data are to be written to the log files. Each record 
specifier is a list of known properties that are to be written as the log record for the log file associated with the record specifier. 
Default record specifiers for each of the standard log files are shown in :doc:`standard-logfiles-default-record-specifications`. 
The standard log file record specifiers can be defined by the user at run-time (see :doc:`standard-logfiles-record-specification`).

When specifying known properties, the property name must be prefixed with the property type. The current list of valid property types 
available for use is:

    - STAR_PROPERTY
    - STAR_1_PROPERTY
    - STAR_2_PROPERTY
    - SUPERNOVA_PROPERTY
    - COMPANION_PROPERTY
    - BINARY_PROPERTY
    - PROGRAM_OPTION

The stellar property types (all types except BINARY_PROPERTY AND PROGRAM_OPTION) must be paired with properties from the stellar property list, 
the binary property type BINARY_PROPERTY with properties from the binary property list, and the program option type PROGRAM_OPTION with properties 
from the program option property list.
SSE system parameters
=====================

 Default record definition for the SSE System Parameters log file::

    const ANY_PROPERTY_VECTOR SSE_SYSTEM_PARAMETERS_REC = {
        STAR_PROPERTY::RANDOM_SEED,
        STAR_PROPERTY::MZAMS,
        STAR_PROPERTY::RZAMS,
        STAR_PROPERTY::METALLICITY,
        STAR_PROPERTY::OMEGA_ZAMS,
        STAR_PROPERTY::INITIAL_STELLAR_TYPE,
        STAR_PROPERTY::STELLAR_TYPE,
        STAR_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,
        STAR_PROPERTY::MASS,
        STAR_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,
        PROGRAM_OPTION::LBV_FACTOR,
        PROGRAM_OPTION::WR_FACTOR,
        PROGRAM_OPTION::NOTES
    };

BSE detailed output
===================

Default record definition for the BSE Detailed Output log file::

    const ANY_PROPERTY_VECTOR BSE_DETAILED_OUTPUT_REC = {
        BINARY_PROPERTY::RANDOM_SEED,
        BINARY_PROPERTY::DT,
        BINARY_PROPERTY::TIME,
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL,
        BINARY_PROPERTY::ECCENTRICITY,
        STAR_1_PROPERTY::MZAMS,
        STAR_2_PROPERTY::MZAMS,
        STAR_1_PROPERTY::MASS_0,
        STAR_2_PROPERTY::MASS_0,
        STAR_1_PROPERTY::MASS,
        STAR_2_PROPERTY::MASS,
        STAR_1_PROPERTY::ENV_MASS,
        STAR_2_PROPERTY::ENV_MASS,
        STAR_1_PROPERTY::CORE_MASS,
        STAR_2_PROPERTY::CORE_MASS,
        STAR_1_PROPERTY::HE_CORE_MASS,
        STAR_2_PROPERTY::HE_CORE_MASS,
        STAR_1_PROPERTY::CO_CORE_MASS,
        STAR_2_PROPERTY::CO_CORE_MASS,
        STAR_1_PROPERTY::RADIUS,
        STAR_2_PROPERTY::RADIUS,
        BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1,
        BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2,
        BINARY_PROPERTY::ROCHE_LOBE_TRACKER_1,
        BINARY_PROPERTY::ROCHE_LOBE_TRACKER_2,
        STAR_1_PROPERTY::OMEGA,
        STAR_2_PROPERTY::OMEGA,
        STAR_1_PROPERTY::OMEGA_BREAK,
        STAR_2_PROPERTY::OMEGA_BREAK,
        STAR_1_PROPERTY::INITIAL_STELLAR_TYPE,
        STAR_2_PROPERTY::INITIAL_STELLAR_TYPE,
        STAR_1_PROPERTY::STELLAR_TYPE,
        STAR_2_PROPERTY::STELLAR_TYPE,
        STAR_1_PROPERTY::AGE,
        STAR_2_PROPERTY::AGE,
        STAR_1_PROPERTY::LUMINOSITY,
        STAR_2_PROPERTY::LUMINOSITY,
        STAR_1_PROPERTY::TEMPERATURE,
        STAR_2_PROPERTY::TEMPERATURE,
        STAR_1_PROPERTY::ANGULAR_MOMENTUM,
        STAR_2_PROPERTY::ANGULAR_MOMENTUM,
        STAR_1_PROPERTY::DYNAMICAL_TIMESCALE,
        STAR_2_PROPERTY::DYNAMICAL_TIMESCALE,
        STAR_1_PROPERTY::THERMAL_TIMESCALE,
        STAR_2_PROPERTY::THERMAL_TIMESCALE,
        STAR_1_PROPERTY::NUCLEAR_TIMESCALE,
        STAR_2_PROPERTY::NUCLEAR_TIMESCALE,
        STAR_1_PROPERTY::ZETA_SOBERMAN,
        STAR_2_PROPERTY::ZETA_SOBERMAN,
        STAR_1_PROPERTY::ZETA_SOBERMAN_HE,
        STAR_2_PROPERTY::ZETA_SOBERMAN_HE,
        STAR_1_PROPERTY::ZETA_HURLEY,
        STAR_2_PROPERTY::ZETA_HURLEY,
        STAR_1_PROPERTY::ZETA_HURLEY_HE,
        STAR_2_PROPERTY::ZETA_HURLEY_HE,
        STAR_1_PROPERTY::MASS_LOSS_DIFF,
        STAR_2_PROPERTY::MASS_LOSS_DIFF,
        STAR_1_PROPERTY::DOMINANT_MASS_LOSS_RATE,
        STAR_2_PROPERTY::DOMINANT_MASS_LOSS_RATE,
        STAR_1_PROPERTY::MASS_TRANSFER_DIFF,
        STAR_2_PROPERTY::MASS_TRANSFER_DIFF,
        BINARY_PROPERTY::TOTAL_ANGULAR_MOMENTUM,
        BINARY_PROPERTY::TOTAL_ENERGY,
        STAR_1_PROPERTY::METALLICITY,
        STAR_2_PROPERTY::METALLICITY,
        BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY,
        STAR_1_PROPERTY::PULSAR_MAGNETIC_FIELD,
        STAR_2_PROPERTY::PULSAR_MAGNETIC_FIELD,
        STAR_1_PROPERTY::PULSAR_SPIN_FREQUENCY,
        STAR_2_PROPERTY::PULSAR_SPIN_FREQUENCY,
        STAR_1_PROPERTY::PULSAR_SPIN_DOWN_RATE,
        STAR_2_PROPERTY::PULSAR_SPIN_DOWN_RATE,
        STAR_1_PROPERTY::RADIAL_EXPANSION_TIMESCALE,
        STAR_2_PROPERTY::RADIAL_EXPANSION_TIMESCALE
    };

BSE system parameters
=====================

Default record definition for the BSE System Parameters log file::

    const ANY_PROPERTY_VECTOR BSE_SYSTEM_PARAMETERS_REC = {
        BINARY_PROPERTY::RANDOM_SEED,
        STAR_1_PROPERTY::MZAMS,
        STAR_2_PROPERTY::MZAMS,
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL,
        BINARY_PROPERTY::ECCENTRICITY_INITIAL,
        STAR_1_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,
        STAR_2_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER,
        STAR_1_PROPERTY::OMEGA_ZAMS,
        STAR_2_PROPERTY::OMEGA_ZAMS,
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_BH,
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_ECSN,
        PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_FOR_USSN,
        PROGRAM_OPTION::LBV_FACTOR,
        PROGRAM_OPTION::WR_FACTOR,
        PROGRAM_OPTION::COMMON_ENVELOPE_ALPHA,
        STAR_1_PROPERTY::METALLICITY,
        STAR_2_PROPERTY::METALLICITY,
        BINARY_PROPERTY::UNBOUND,
        BINARY_PROPERTY::STELLAR_MERGER,
        BINARY_PROPERTY::STELLAR_MERGER_AT_BIRTH,
        BINARY_PROPERTY::MASSES_EQUILIBRATED_AT_BIRTH,
        STAR_1_PROPERTY::INITIAL_STELLAR_TYPE,
        STAR_1_PROPERTY::STELLAR_TYPE,
        STAR_2_PROPERTY::INITIAL_STELLAR_TYPE,
        STAR_2_PROPERTY::STELLAR_TYPE,
        STAR_1_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,
        STAR_2_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE,
        BINARY_PROPERTY::ERROR,
        PROGRAM_OPTION::NOTES
    };

BSE double compact objects
==========================

Default record definition for the BSE Double Compact Objects log file::

    const ANY_PROPERTY_VECTOR BSE_DOUBLE_COMPACT_OBJECTS_REC = {
        BINARY_PROPERTY::RANDOM_SEED,
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_AT_DCO_FORMATION, 
        BINARY_PROPERTY::ECCENTRICITY_AT_DCO_FORMATION,
        STAR_1_PROPERTY::MASS,
        STAR_1_PROPERTY::STELLAR_TYPE,
        STAR_2_PROPERTY::MASS, 
        STAR_2_PROPERTY::STELLAR_TYPE,
        BINARY_PROPERTY::TIME_TO_COALESCENCE,
        BINARY_PROPERTY::TIME,
        BINARY_PROPERTY::MERGES_IN_HUBBLE_TIME, 
        STAR_1_PROPERTY::RECYCLED_NEUTRON_STAR,  
        STAR_2_PROPERTY::RECYCLED_NEUTRON_STAR
    };

    COMPAS output
=============

Summary and status information during the evolution of stars is written to stdout; how much is written depends upon the value of the 
``--quiet`` program option.

Detailed information is written to log files (described below). All COMPAS output files are created inside a container directory, specified 
by the ``--output-container`` program option.

If detailed output log files are created (see the ``--detailed-output`` program option), they will be created inside a containing directory 
named ``Detailed_Output`` within the COMPAS output container directory.

Also created in the COMPAS container directory is a file named ``Run_Details`` in which COMPAS records some details of the run (COMPAS 
version, start time, command line program option values etc.). Note that the option values recorded in the Run details file are the values
specified on the command line, not the values specified in a ``grid`` file (if used).


.. toctree::
   :maxdepth: 1

   standard-logfiles
   standard-logfiles-record-specifiers
   standard-logfiles-record-specification
   standard-logfiles-format
   standard-logfiles-annotations
   standard-logfiles-example-definitions-file
   standard-logfiles-default-record-specifications
Log file annotations
====================

COMPAS provides functionality that allows users to annotate standard log files.

The original motivation for log file annotation functionality was to enable users to describe the contents of custom grid files,
but annotations can be used for any reason.  With the ability to annotate the log files, users can indicate the origin of various
input data - e.g. the user could indicate what IMF was used to draw initial mass values, or what distribution was used to draw 
mass ratio (q), etc.  

Annotations are written to log files as columns of data (aka datasets in ``HDF5`` files). Annotations are specified via program
options, and so can be specified on the command line as well as in grid files. Because annotations are written to log files
as columns of data, functionality is provided to allow the user to specify header strings for the columns, as well as the
actual data for the columns (the annotations).

Two program options are provided:

``--notes-hdrs``

Allows users to specify header strings for annotation columns.  Can only be specified on the command line: not valid in grid
files, and will be ignored, and elicit a warning, if included.

The ``notes-hdrs`` program option is a vector program option, and provides for the specification of one or more annotation
header strings.  Usage is::

    --notes-hdrs headerstr1 headerstr2 headerstr3 ... headerstrN

Header strings are separated by a space character. There is no limit to the number of header strings specified, and the number
specified defines the maximum number of annotations allowed.


``--notes``

Allows users to specify annotations.  Can be specified on the command line and in a grid file.

The ``notes`` program option is a vector program option, and provides for the specification of one or more annotation strings.
Usage is::

    --notes annotation1 annotation2 annotation3 ... annotationN

Annotation strings are separated by a space character.  The number of annotation strings is limited to the number of annotation
header strings specified (via the ``--notes-hdrs`` program option).  If more annotation strings are specified than header strings,
the excess annotation strings will be ignored (and a warning displayed).

When using this notation, all annotation strings must be provided: there is no mechanism to allow a default annotation using the
fallback method for program options (to the command-line value, then to the COMPAS default) - leaving an annotation string blank
would be ambiguous (as to which annotation string had been left blank), and specifying an empty string ("") as an annotation
string would be ambiguous (as to whether the user wanted the annotation string to default, or just be an empty string).

Neither header strings, nor annotation strings, may begin with the dash character ('-'), because the shell parser will parse them
as option names before passing them through to COMPAS.

COMPAS imposes no limit to the length of an individual annotation header string or annotation string (number of characters), but
there may be practical limits imposed by the underlying system. 


Using vector option shorthand notation
--------------------------------------

Because the notation described above could become awkward, and to allow for default annotations, a shorthand notation for vector
program options is provided (see :doc:`../Program options/program-options-vector-options` for details).

Usage using the shorthand notation is::

    --notes-hdrs [headerstr1,headerstr2,headerstr3,...,headerStrN]

    --notes [annotation1,annotation2,annotation3,...,annotationN]

Because the parameters are bounded by the brackets, and delimited by commas (and so are now positional), users can omit specific
annotations::

    --notes [,,annotation3,,annotation5]

In the example above, annotations 1, 2, 4, and those beyond annotation 5 have been omitted. Annotations 1, 2 & 4 will default - if
they are specified in this manner on a grid line they will default to the correspodning annotation specified on the command line; if
they are specified in this manner on the command line they will default to the COMPAS default annotation (the empty string).  In the
example above, if the number of annotations expected, as defined by the number of annotation headers specified via the ``--notes-hdrs``
program option is more than 5, then anotations beyond annotation 5 (the last annotation actually specified by the user) will default
in the same manner as described above.

Spaces in annotation header strings and annotation strings need to be enclosed in quotes, or the shell parser will parse them as 
separate arguments before passing them through to COMPAS.  If the log file type is specified as TXT, then any spaces in annotation
header strings and annotation strings need to be enclosed in quotes to avoid the shell parser parsing them as separate arguments, but
also need to have enclosing quotes propagated to the logfile, or the spaces will be interpreted as delimiters in the log file - in this
case, the user will need to add enclosing escaped quote characters ('\"') before adding the enclosing quotes.  e.g.::

    --notes-hdrs [headerstr1,"\"header str 2\"",headerstr3,...,headerStrN]

The shorthand notation is expanded to the notation described above (the COMPAS code just fills in the omitted annotation strings with
the required defaults), so the caveat mentioned above that neither header strings, nor annotation strings, may begin with the dash 
character ('-') applies to the shorthand notation.

Shorthand notation is optional: users may choose to use the notation described above rather than shorthand notation, but in that case
all annotation strings must be specified (no omissions, no defaults).


Including annotations in log files
----------------------------------

Annotations can be included in the record specifiers for any of the standard log files in the same way that other program options
are included: they can be included in the default log file record specifiers in ``constants.h`` (compile-time configuration: see
:doc:`../../Developer guide/constants-dot-h`), and/or they can be added to, or removed from, the default log file record specifiers
via the use of a log file definitions file (run-time configuration: see :doc:`./standard-logfiles-record-specification`).


Compile-time configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

The property specifier ``PROGRAM_OPTION::NOTES`` is available to be included in any of the default log file record specifiers in
``constants.h``.  If the property specifier ``PROGRAM_OPTION::NOTES`` is included in a default log file record specifier:

    - a column will be included in the log file for each of the annotation header strings specfied by the user (via the ``notes-hdrs``
      program option), with each column having the respective header string specified by the user
    - each record in the log file will include an annotation string for each annotation column, the contents of which is determined by
      what the user specified using the ``--notes`` program option.

Adding the property specifier ``PROGRAM_OPTION::NOTES`` to the default record specifier in ``constants.h`` for a log file causes
*all* annotation columns (as specified via the ``notes-hdrs`` program option) to be included in the log file: including annotations
in a log file via compile-time configuration is an all-or-nothing proposition.


Run-time configuration
~~~~~~~~~~~~~~~~~~~~~~

COMPAS provides functionality to allow users to change which properties are to be written to the standard log files at run-time:
see :doc:`./standard-logfiles-record-specification`. The property specifier ``PROGRAM_OPTION::NOTES`` can be added to, or removed from,
any of the log file record specifiers by the use of this functionality.

Furthermore, because at run-time the number of annotation columns is known (information not known at compile-time), the 
``PROGRAM_OPTION::NOTES`` property specifier can (optionally) be indexed to allow the specification of a particular annotation column.
Thus, ``PROGRAM_OPTION::NOTES`` (with no index) indicates *all* annotations columns, whereas ``PROGRAM_OPTION::NOTES[2]`` indicates 
annotation column 2 (1-based: the first annotation column is indicated by ``PROGRAM_OPTION::NOTES[1]``). By using the optional index, 
users can add specific annotations columns to, or remove them from, any of the log files.

For example, this log file definitions file entry::

    SSE_SYSPARMS_REC -= { PROGRAM_OPTION::NOTES }

removes all annotations columns from the SSE System Parameters log file.

This entry::

    BSE_PULSARS_REC += { PROGRAM_OPTION::NOTES[1] }

adds annotations column 1 to the BSE Pulsar Evolution log file.

These entries::

    BSE_SYSPARMS_REC -= { PROGRAM_OPTION::NOTES }
    BSE_SYSPARMS_REC += { PROGRAM_OPTION::NOTES[1], PROGRAM_OPTION::NOTES[3] }

    BSE_SNE_REC += { PROGRAM_OPTION::NOTES[2], PROGRAM_OPTION::NOTES[3],  PROGRAM_OPTION::NOTES[7] }

    BSE_RLOF_REC -= { PROGRAM_OPTION::NOTES[5] }

    BSE_CEE_REC += { PROGRAM_OPTION::NOTES }

will:

    - remove all annotations columns from the BSE System Parameters log file
    - add annotations columns 1 and 3 to the BSE System Parameters log file
    - add annotations columns 2, 3, and 7 to the BSE Supernovae log file
    - remove annotations column 5 from the BSE RLOF log file
    - add all annotation columns to the BSE Common envelopes files

Specifying an index value less than 1 or greater than the number of annotation headers specified will result in an error.

See :doc:`./standard-logfiles-record-specification` for more details.

|br|
The property ``PROGRAM_OPTION::NOTES`` is included in the default record specifier in ``constants.h`` for both the SSE System Parameters
log file, and the BSE System Parameters log file.

Note that whichever configuration method is used to include annotations in log files, if no annotation headers are specified via the 
``notes-hdrs`` program option, no annotations will be included in any log file.


Log file space considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Annotations are strings of (effectively) unlimited length, so if several annotations are recorded in a log file for each system evolved,
the storage space required for the log file will grow accordingly. Sometimes annotations will apply to the population rather than to each
individual system, and in those cases it may not be necessary to record the same annotation string for possibly millions of systems in the
log file(s). In such situations, it is likely that the annotation headers and strings will be specified on the command line rather than
for each system in a grid file, and so will be recorded in the ``Run_Details`` file - so removing the ``PROGRAM_OPTION::NOTES`` property
from the log file record specifications will prevent them from being repeated needlessly in the log files, and they can be retrieved as
required from the ``Run_Details`` file.

Standard log file record specification
======================================

The standard log file record specifiers can be changed at run-time by supplying a log file definitions file via the ``--logfile-definitions`` 
program option. This allows users to change what appears in any of the standard COMPAS log files without the need to change the source code
or rebuild the COMPAS executable.

The syntax of the definitions file is fairly simple. The definitions file is expected to contain zero or more log file record specifications, as 
explained below.

Complete lists of available properties selectable for inclusion (or exclusion) from COMPAS log file are availabe:


.. toctree::
   :maxdepth: 1

   standard-logfiles-record-specification-stellar
   standard-logfiles-record-specification-binary
   standard-logfiles-record-specification-options


For the following specification::

     ::=   means `expands to` or `is defined as`
    { x }  means (possible) repetition: x may appear zero or more times
    [ x ]  means x is optional: x may appear, or not
    <name> is a term (expression)
    "abc"  means literal string "abc"
      n    means integer number
      |    means `OR`
      #    indicates the start of a comment


Log file definitions file specification
---------------------------------------

::

    <def_file>   ::= {<rec_spec>}

    <rec_spec>   ::= <rec_name> <op> "{" { [ <props_list> ] } "}" <spec_delim>

    <rec_name>   ::= "SSE_SYSPARMS_REC"    |   # SSE only
                     "SSE_DETAILED_REC"    |   # SSE only
                     "SSE_SNE_REC"         |   # SSE only
                     "SSE_SWITCH_REC"      |   # SSE only
                     "BSE_SYSPARMS_REC"    |   # BSE only
                     "BSE_SWITCH_REC"      |   # BSE only
                     "BSE_DCO_REC"         |   # BSE only
                     "BSE_SNE_REC"         |   # BSE only
                     "BSE_CEE_REC"         |   # BSE only
                     "BSE_PULSARS_REC"     |   # BSE only
                     "BSE_RLOF_REC"        |   # BSE only
                     "BSE_DETAILED_REC"    |   # BSE only
   
    <op>         ::= "=" | "+=" | "-="

    <props_list> ::= <prop_spec> [ <props_delim> <props_list> ]

    <prop_spec>  ::= <prop_type> "::" <prop_name> [ <prop_index> ] <prop_delim>

    <spec_delim> ::= "::" | "EOL"

    <prop_delim> ::= "," | <spec_delim>

    <prop_type>  ::= "STAR_PROPERTY"       |   # SSE only
                     "STAR1_PROPERTY"      |   # BSE only
                     "STAR2_PROPERTY"      |   # BSE only
                     "SUPERNOVA_PROPERTY"  |   # BSE only
                     "COMPANION_PROPERTY"  |   # BSE only
                     "BINARY_PROPERTY"     |   # BSE only
                     "PROGRAM_OPTION"      |   # SSE_or BSE

    <prop_index> ::= "[" n "]"

    <prop_name>  ::= valid property name for specified property type (see definitions in constants.h)


The log file definitions may contain comments. Comments are denoted by the hash/pound character ("#"). The hash/pound character and any 
text following it on the line in which the hash character appears is ignored by the parser. The hash/pound character can appear anywhere 
on a line - if it is the first character then the entire line is a comment and ignored by the parser, or it can follow valid symbols on 
a line, in which case the symbols before the hash/pound character are parsed and interpreted by the parser.

A log file specification record is initially set to its default value (see :doc:`standard-logfiles-default-record-specifications`). 
The log file definitions file informs COMPAS as to the modifications to the default values the user wants. This means that the log file 
definitions log file is not mandatory, and if the log file definitions file is not present, or contains no valid record specifiers, the log 
file record definitions will remain at their default values.

The assignment operator given in a record specification (``<op>`` in the specification above) can be one of "=", "+=", and "-=".  The meanings of 
these are:

    "=" means the record specifier should be assigned the list of properties specified in the braced-list following the "=" operator. The value of the record specifier prior to the assignment is discarded, and the new value set as described.
    
    "+=" means the list of properties specified in the braced-list following the "+=" operator should be appended to the existing value of the record specifier. Note that the new properties are appended to the existing list, so will appear at the end of the list (properties are printed in the order they appear in the list).
    
    "-=" means the list of properties specified in the braced-list following the "-=" operator should be subtracted from the existing value of the record specifier.

The `prop-index` qualifier may only be used for vector program options.


Example definitions file entries
--------------------------------

Some example log file definitions file entries are::

    SSE_SYSPARMS_REC = { STAR_PROPERTY::RANDOM_SEED, STAR_PROPERTY::RADIUS,
                         STAR_PROPERTY::MASS,
                         STAR_PROPERTY::LUMINOSITY
                       }

    BSE_PULSARS_REC += { STAR_1_PROPERTY::LUMINOSITY, STAR_2_PROPERTY::CORE_MASS, 
                         BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL, COMPANION PROPERTY::RADIUS
                       }

    BSE_PULSARS_REC -= { SUPERNOVA_PROPERTY::TEMPERATURE }

    BSE_PULSARS_REC += { PROGRAM_OPTION::KICK_MAGNITUDE_DISTRIBUTION_SIGMA_CCSN_NS,
                         BINARY_PROPERTY::ORBITAL_VELOCITY
                       }

    BSE_SYSPARMS_REC -= { PROGRAM_OPTION::NOTES }
    BSE_SYSPARMS_REC += { PROGRAM_OPTION::NOTES[1], PROGRAM_OPTION::NOTES[3] }
    
A full example log file record specifications file is shown in :doc:`standard-logfiles-example-definitions-file`.

The record specifications in the definitions file are processed individually in the sequence they appear in the file, and are cumulative: 
for record specifications pertaining to the same record name, the output of earlier specifications is input to later specifications.

For each record specification:

- Properties requested to be added to an existing record specification that already exist in that record specification are ignored. Properties will not appear in a record specification twice.
- Properties requested to be subtracted from an existing record specification that do not exist in that record specification are ignored.

Note that neither of those circumstances will cause a parse error for the definitions file – in both cases the user’s intent is satisfied.
BSE RLOF parameters
===================

Default record definition for the BSE RLOF Parameters log file::

    const ANY_PROPERTY_VECTOR BSE_RLOF_PARAMETERS_REC = {
        BINARY_PROPERTY::RLOF_POST_MT_RANDOM_SEED,
        BINARY_PROPERTY::RLOF_POST_MT_STAR1_MASS,
        BINARY_PROPERTY::RLOF_POST_MT_STAR2_MASS,
        BINARY_PROPERTY::RLOF_POST_MT_STAR1_RADIUS,
        BINARY_PROPERTY::RLOF_POST_MT_STAR2_RADIUS,
        BINARY_PROPERTY::RLOF_POST_MT_STAR1_STELLAR_TYPE,
        BINARY_PROPERTY::RLOF_POST_MT_STAR2_STELLAR_TYPE,
        BINARY_PROPERTY::RLOF_POST_MT_SEMI_MAJOR_AXIS,
        BINARY_PROPERTY::RLOF_POST_MT_ECCENTRICITY,
        BINARY_PROPERTY::RLOF_POST_MT_EVENT_COUNTER,
        BINARY_PROPERTY::RLOF_POST_MT_TIME,
        BINARY_PROPERTY::RLOF_POST_MT_STAR1_RLOF,
        BINARY_PROPERTY::RLOF_POST_MT_STAR2_RLOF,
        BINARY_PROPERTY::RLOF_POST_MT_COMMON_ENVELOPE,
        BINARY_PROPERTY::RLOF_PRE_MT_STAR1_MASS,
        BINARY_PROPERTY::RLOF_PRE_MT_STAR2_MASS,
        BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RADIUS,
        BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RADIUS,
        BINARY_PROPERTY::RLOF_PRE_MT_STAR1_STELLAR_TYPE,
        BINARY_PROPERTY::RLOF_PRE_MT_STAR2_STELLAR_TYPE,
        BINARY_PROPERTY::RLOF_PRE_MT_SEMI_MAJOR_AXIS,
        BINARY_PROPERTY::RLOF_PRE_MT_ECCENTRICITY,
        BINARY_PROPERTY::RLOF_PRE_MT_TIME,
        BINARY_PROPERTY::RLOF_PRE_MT_STAR1_RLOF,
        BINARY_PROPERTY::RLOF_PRE_MT_STAR2_RLOF,
        STAR_1_PROPERTY::ZETA_SOBERMAN,
        STAR_1_PROPERTY::ZETA_SOBERMAN_HE,
        STAR_1_PROPERTY::ZETA_HURLEY,
        STAR_1_PROPERTY::ZETA_HURLEY_HE,
        STAR_2_PROPERTY::ZETA_SOBERMAN,
        STAR_2_PROPERTY::ZETA_SOBERMAN_HE,
        STAR_2_PROPERTY::ZETA_HURLEY,
        STAR_2_PROPERTY::ZETA_HURLEY_HE
    };

Example log file definitions file
=================================

Following is an example log file definitions file. COMPAS can be configured to use this file via the ``--logfile-definitions`` program option.

This file (``COMPAS_Output_Definitions.txt``) is also delivered as part of the COMPAS github repository.

::

    # sample standard log file specifications file
    
    # the ’#’ character and anything following it on a single line is considered a comment
    # (so, lines starting with ’#’ are comment lines)
    
    # case is not significant
    # specifications can span several lines
    # specifications for the same log file are cumulative
    # if a log file is not specified in this file, the default specification is used
    
    
    # SSE Parameters
    
    # start with the default SSE Parameters specification and add ENV_MASS
    
    sse_sysparms_rec += { STAR_PROPERTY::ENV_MASS }
    
    # take the updated SSE Parameters specification and add ANGULAR_MOMENTUM
    
    sse_sysparms_rec += { STAR_PROPERTY::ANGULAR_MOMENTUM }
    
    # take the updated SSE Parameters specification and subtract MASS_0 and MDOT
    
    sse_sysparms_rec -= { STAR_PROPERTY::MASS_0, STAR_PROPERTY::MDOT }
    
    
    # BSE System Parameters
    
    bse_sysparms_rec = {                  # set the BSE System Parameters specification to:
        BINARY_PROPERTY::ID,              # ID of the binary
        BINARY_PROPERTY::RANDOM_SEED,     # RANDOM_SEED for the binary
        STAR_1_PROPERTY::MZAMS,           # MZAMS for Star1
        STAR_2_PROPERTY::MZAMS            # MZAMS for Star2
    }
    
    # ADD to the BSE System Parameters specification:
    # SEMI_MAJOR_AXIS_INITIAL for the binary
    # ECCENTRICITY_INITIAL for the binary
    # SUPERNOVA_THETA for Star1 and SUPERNOVA_PHI for Star1
    
    bse_sysparms_rec += {
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_INITIAL,
        BINARY_PROPERTY::ECCENTRICITY_INITIAL,
        STAR_1_PROPERTY::SUPERNOVA_THETA, STAR_1_PROPERTY::SUPERNOVA_PHI
    }
    
    bse_sysparms_rec += {                          # ADD to the BSE System Parameters specification:
        SUPERNOVA_PROPERTY::IS_ECSN,               # IS_ECSN for the supernova star
        SUPERNOVA_PROPERTY::IS_SN,                 # IS_SN for the supernova star
        SUPERNOVA_PROPERTY::IS_USSN,               # IS_USSN for the supernova star
        SUPERNOVA_PROPERTY::EXPERIENCED_PISN,      # EXPERIENCED_PISN for the supernova star
        SUPERNOVA_PROPERTY::EXPERIENCED_PPISN,     # EXPERIENCED_PPISN for the supernova star
        BINARY_PROPERTY::UNBOUND,                  # UNBOUND for the binary
        SUPERNOVA_PROPERTY::MZAMS,                 # MZAMS for the supernova star
        COMPANION_PROPERTY::MZAMS                  # MZAMS for the companion star
    }
    
    # SUBTRACT from the BSE System Parameters specification:
    # RANDOM_SEED for the binary
    # ID for the binary
    # all ANNOTATIONS
    
    bse_sysparms_rec -= {                 # SUBTRACT from the BSE System Parameters specification:
        BINARY_PROPERTY::RANDOM_SEED,     # RANDOM_SEED for the binary
        BINARY_PROPERTY::ID,              # ID for the binary
        PROGRAM_OPTION::NOTES             # all ANNOTATIONS
    }
    
    bse_sysparms_rec += { PROGRAM_OPTION::NOTES[1], PROGRAM_OPTION::NOTES[3] } # ADD ANNOTATIONS 1 & 3 to the BSE System Parameters specification


    # BSE Double Compact Objects
    
    # set the BSE Double Compact Objects specification to MZAMS for Star1, and MZAMS for Star2
    
    BSE_DCO_Rec = { STAR_1_PROPERTY::MZAMS, STAR_2_PROPERTY::MZAMS }
    
    # set the BSE Double Compact Objects specification to empty - nothing will be printed
    # (file will not be created)
    
    BSE_DCO_Rec = {}
    
    
    # BSE Supernovae
    
    BSE_SNE_Rec = {}     # set spec empty - nothing will be printed (file will not be created)
    
    
    # BSE Common Envelopes
    
    BSE_CEE_Rec = {}     # set spec empty - nothing will be printed (file will not be created)
    
    
    # BSE Pulsars
    
    # line ignored (comment). BSE Pulsars specification will be default
    
    # BSE_Pulsars_Rec = { STAR_1_PROPERTY::MASS, STAR_2_PROPERTY::MASS }
    
    
    # BSE Detailed Output
    
    BSE_Detailed_Rec = {} # set spec empty - nothing will be printed (file will not be created)
BSE switchlog
=============

Default record definition for the BSE SwitchLog log file::

    const ANY_PROPERTY_VECTOR BSE_SWITCH_LOG_REC = {
        BINARY_PROPERTY::RANDOM_SEED,
        BINARY_PROPERTY::TIME
    };


The default record specification can be modified at runtime via a logfile record specifications file (program option ``--logfile-definitions``).
See :doc:`standard-logfiles-record-specification` for details.

Note that the BSE SwitchLog file has the following columns automatically appended to each record:

    - The constituent star switching stellar type: 1 = Primary, 2 = Secondary.
    - The stellar type from which the star is switching.
    - The stellar type to which the star is switching.

|br|
**STAR_SWITCHING**

.. list-table::
   :widths: 20 80 
   :header-rows: 0
   :class: aligned-text

   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` ``BaseBinaryStar::m_Star1/m_Star2``
   * - Description:
     - The constituent star switching stellar type, where 1 = Primary, and 2 = Secondary.
   * - Header String:
     - "STAR_SWITCHING"

**SWITCHING_FROM**

.. list-table::
   :widths: 20 80 
   :header-rows: 0
   :class: aligned-text

   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` ``BaseStar::m_StellarType``
   * - Description:
     - The stellar type of the constituent star immediately prior to the switch.
   * - Header String:
     - "SWITCHING_FROM"

**SWITCHING_TO**

.. list-table::
   :widths: 20 80 
   :header-rows: 0
   :class: aligned-text

   * - Data type:
     - INT
   * - COMPAS variable:
     - Not applicable
   * - Description:
     - The stellar type to which the constituent star will switch (i.e. the stellar type immediately following the switch).
   * - Header String:
     - "SWITCHING_TO"

These columns will always be automatically appended to each BSE Switch Log record: they cannot be removed via the log file record 
specifications file.
Stellar properties
==================

When specifying known properties in a log file record specification record, the property name must be prefixed with 
the property type. The current list of valid stellar property types available for use is:

    - STAR_PROPERTY for all stars for ``SSE``
    - STAR_1_PROPERTY for the primary star of a binary for ``BSE``
    - STAR_2_PROPERTY for the secondary star of a binary for ``BSE``
    - SUPERNOVA_PROPERTY for the exploding star in a supernova event for ``BSE``
    - COMPANION_PROPERTY for the companion star in a supernova event for ``BSE``

For example, to specify the property ``TEMPERATURE`` for an individual star being evolved for ``SSE``, use::

    STAR_PROPERTY::TEMPERATURE


.. _stellar-props-top:

Following is an alphabetical list of stellar properties available for inclusion in log file record specifiers.

**Jump to**
:ref:`A <stellar-props-A>` :ref:`B <stellar-props-B>` :ref:`C <stellar-props-C>` :ref:`D <stellar-props-D>`
:ref:`E <stellar-props-E>` :ref:`F <stellar-props-F>` :ref:`G <stellar-props-G>` :ref:`H <stellar-props-H>`
:ref:`I <stellar-props-I>` :ref:`J <stellar-props-J>` :ref:`K <stellar-props-K>` :ref:`L <stellar-props-L>`
:ref:`M <stellar-props-M>` :ref:`N <stellar-props-N>` :ref:`O <stellar-props-O>` :ref:`P <stellar-props-P>`
:ref:`Q <stellar-props-Q>` :ref:`R <stellar-props-R>` :ref:`S <stellar-props-S>` :ref:`T <stellar-props-T>`
:ref:`U <stellar-props-U>` :ref:`V <stellar-props-V>` :ref:`W <stellar-props-W>` :ref:`X <stellar-props-X>`
:ref:`Y <stellar-props-Y>` :ref:`Z <stellar-props-Z>` |_| |_| |_| |_| :ref:`Supernova events/states <supernova-events-states>`

.. _stellar-props-A:

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **AGE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Age
   * - Description:
     - Effective age (changes with mass loss/gain) (Myr)
   * - Header Strings:
     - Age, Age(1), Age(2), Age(SN), Age(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ANGULAR_MOMENTUM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_AngularMomentum
   * - Description:
     - Angular momentum (\ :math:`M_\odot AU^2 yr^{−1}`)
   * - Header Strings:
     - Ang_Momentum, Ang_Momentum(1), Ang_Momentum(2), Ang_Momentum(SN), Ang_Momentum(CP)

.. _stellar-props-B:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.bindingEnergy
   * - Description:
     - Absolute value of the envelope binding energy at the onset of unstable RLOF (erg). Used for calculating post-CE separation.
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     -  Binding_Energy@CE(1), Binding_Energy@CE(2), Binding_Energy@CE(SN), Binding_Energy@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_FIXED**
     -
   * - Data type:
     -  DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.fixed
   * - Description:
     - Absolute value of the envelope binding energy calculated using a fixed lambda parameter (erg). Calculated using lambda = m_Lambdas.fixed.
   * - Header Strings:
     - BE_Fixed, BE_Fixed(1), BE_Fixed(2), BE_Fixed(SN), BE_Fixed(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_KRUCKOW**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.kruckow
   * - Description:
     - Absolute value of the envelope binding energy calculated using the fit by :cite:`Vigna-Gomez2018` to :cite:`Kruckow2016` (erg). Calculated using alpha = OPTIONS→CommonEnvelopeSlopeKruckow().
   * - Header Strings:
     - BE_Kruckow, BE_Kruckow(1), BE_Kruckow(2), BE_Kruckow(SN), BE_Kruckow(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_LOVERIDGE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.loveridge
   * - Description:
     - Absolute value of the envelope binding energy calculated as per :cite:`Loveridge2011` (erg). Calculated using lambda = m_Lambdas.loveridge.
   * - Header Strings:
     - BE_Loveridge, BE_Loveridge(1), BE_Loveridge(2), BE_Loveridge(SN), BE_Loveridge(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_LOVERIDGE_WINDS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.loveridgeWinds
   * - Description:
     - Absolute value of the envelope binding energy calculated as per :cite:`Webbink1984` & :cite:`Loveridge2011` including winds (erg). Calculated using lambda = m_Lambdas.loveridgeWinds.
   * - Header Strings:
     - BE_Loveridge_Winds, BE_Loveridge_Winds(1), BE_Loveridge_Winds(2), BE_Loveridge_Winds(SN), BE_Loveridge_Winds(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_NANJING**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_BindingEnergies.nanjing
   * - Description:
     - Absolute value of the envelope binding energy calculated as per :doc:`Xu & Li (2010) <../../references>` (erg). Calculated using lambda = m_Lambdas.nanjing.
   * - Header Strings:
     - BE_Nanjing, BE_Nanjing(1), BE_Nanjing(2), BE_Nanjing(SN), BE_Nanjing(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.bindingEnergy
   * - Description:
     - Absolute value of the binding energy immediately after CE (erg).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Binding_Energy>CE(1), Binding_Energy>CE(2), Binding_Energy>CE(SN), Binding_Energy>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **BINDING_ENERGY_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.bindingEnergy
   * - Description:
     - Absolute value of the binding energy at the onset of unstable RLOF leading to the CE (erg). 
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Binding_Energy<CE(1), Binding_Energy<CE(2), Binding_Energy<CE(SN), Binding_Energy<CE(CP)

.. _stellar-props-C:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable: BaseStar::m_CHE
   * - Description:
     - Flag to indicate whether the star evolved as a ``CH`` star for its entire MS lifetime.
   * -
     - TRUE indicates star evolved as ``CH`` star for entire MS lifetime.
   * -
     - FALSE indicates star spun down and switched from ``CH`` to a ``MS_gt_07`` star.
   * - Header Strings:
     - CH_on_MS, CH_on_MS(1), CH_on_MS(2), CH_on_MS(SN), CH_on_MS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CO_CORE_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_COCoreMass
   * - Description:
     - Carbon-Oxygen core mass (\ :math:`M\odot`).
   * - Header Strings:
     - Mass_CO_Core, Mass_CO_Core(1), Mass_CO_Core(2), Mass_CO_Core(SN), Mass_CO_Core(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CO_CORE_MASS_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.COCoreMass
   * - Description:
     - Carbon-Oxygen core mass at the onset of unstable RLOF leading to the CE (\ :math:`M\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Mass_CO_Core@CE(1), Mass_CO_Core@CE(2), Mass_CO_Core@CE(SN), Mass_CO_Core@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.COCoreMassAtCOFormation
   * - Description:
     - Carbon-Oxygen core mass immediately prior to a supernova (\ :math:`M\odot`).
   * - Header Strings:
     - Mass CO_Core@\ CO, Mass_CO_Core@CO(1), Mass_CO_Core@CO(2), Mass_CO_Core@CO(SN), Mass_CO_Core@CO(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CORE_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_CoreMass
   * - Description:
     - Core mass (\ :math:`M\odot`).
   * - Header Strings:
     - Mass_Core, Mass_Core(1), Mass_Core(2), Mass_Core(SN), Mass_Core(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CORE_MASS_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.CoreMass
   * - Description:
     - Core mass at the onset of unstable RLOF leading to the CE (\ :math:`M\odot`)
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Mass_Core@CE(1), Mass_Core@CE(2), Mass_Core@CE(SN), Mass_Core@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **CORE_MASS_AT_COMPACT_OBJECT_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.CoreMassAtCOFormation
   * - Description:
     - Core mass immediately prior to a supernova (\ :math:`M\odot`).
   * - Header Strings:
     - Mass_Core@\ CO, Mass_Core@CO(1), Mass_Core@CO(2), Mass_Core@CO(SN), Mass_Core@CO(CP)

.. _stellar-props-D:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DRAWN_KICK_MAGNITUDE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.drawnKickMagnitude
   * - Description:
     - Magnitude of natal kick without accounting for fallback (\ :math:`km s^{−1}`). Supplied by user in grid file or drawn from distribution (default). This value is used to calculate the actual kick magnitude.
   * - Header Strings:
     - Drawn_Kick_Magnitude, Drawn_Kick_Magnitude(1), Drawn_Kick_Magnitude(2), Drawn_Kick_Magnitude(SN), Drawn_Kick_Magnitude(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DOMINANT_MASS_LOSS_RATE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseStar::m_DMLR
   * - Description:
     - Current dominant mass loss rate printed as one of

        .. list-table::
           :widths: 85 15
           :header-rows: 0
           :class: aligned-text

           * - None 
             - = 0
           * - Nieuwenhuijzen and de Jager 
             - = 1
           * - Kudritzki and Reimers 
             - = 2
           * - Vassiliadis and Wood 
             - = 3
           * - Wolf-Rayet-like (Hamann, Koesterke and de Koter) 
             - = 4
           * - Vink 
             - = 5
           * - Luminous Blue Variable 
             - = 6

   * - Header Strings:
     - Dominant_Mass_Loss_Rate, Dominant_Mass_Loss_Rate(1), Dominant_Mass_Loss_Rate(2), Dominant_Mass_Loss_Rate(SN), Dominant_Mass_Loss_Rate(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Dt
   * - Description: 
     - Current timestep (Myr).
   * - Header Strings: 
     - dT, dT(1), dT(2), dT(SN), dT(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DYNAMICAL_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_DynamicalTimescale
   * - Description:
     - Dynamical time (Myr).
   * - Header Strings:
     - Tau_Dynamical, Tau_Dynamical(1), Tau_Dynamical(2), Tau_Dynamical(SN), Tau_Dynamical(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DYNAMICAL_TIMESCALE_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.dynamicalTimescale
   * - Description:
     - Dynamical time immediately following common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Dynamical>CE(1), Tau_Dynamical>CE(2), Tau_Dynamical>CE(SN), Tau_Dynamical>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.dynamicalTimescale
   * - Description:
     - Dynamical timescale immediately prior to common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Dynamical<CE(1), Tau_Dynamical<CE(2), Tau_Dynamical<CE(SN), Tau_Dynamical<CE(CP)

.. _stellar-props-E:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ECCENTRIC_ANOMALY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.eccentricAnomaly
   * - Description:
     - Eccentric anomaly calculated using Kepler’s equation.
   * - Header Strings:
     - Eccentric_Anomaly, Eccentric_Anomaly(1), Eccentric_Anomaly(2), Eccentric_Anomaly(SN), Eccentric_Anomaly(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ENV_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_EnvMass
   * - Description:
     - Envelope mass calculated using :cite:`Hurley2000` (\ :math:`M\odot`).
   * - Header Strings:
     - Mass_Env, Mass_Env(1), Mass_Env(2), Mass_Env(SN), Mass_Env(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ERROR**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseStar::m_Error
   * - Description:
     - Error number (if error condition exists, else 0).
   * - Header Strings:
     - Error, Error(1), Error(2), Error(SN), Error(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_CCSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as a core-collapse supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_CCSN, Experienced_CCSN(1), Experienced_CCSN(2), Experienced_CCSN(SN), Experienced_CCSN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_ECSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as an electron-capture supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_ECSN, Experienced_ECSN(1), Experienced_ECSN(2), Experienced_ECSN(SN), Experienced_ECSN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_PISN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as an pair-instability supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_PISN, Experienced_PISN(1), Experienced_PISN(2), Experienced_PISN(SN), Experienced_PISN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_PPISN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as a pulsational pair-instability supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_PPISN, Experienced_PPISN(1), Experienced_PPISN(2), Experienced_PPISN(SN), Experienced_PPISN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BinaryConstituentStar::m_RLOFDetails.experiencedRLOF
   * - Description:
     - Flag to indicate whether the star has overflowed its Roche Lobe at any time prior to the current timestep.
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Experienced_RLOF(1), Experienced_RLOF(2), Experienced_RLOF(SN), Experienced_RLOF(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_SN_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - The type of supernova event experienced by the star prior to the current timestep. Printed as one of

        .. list-table::
           :widths: 10 5
           :header-rows: 0
           :class: aligned-text

           * - NONE
             - = 0
           * - CCSN
             - = 1
           * - ECSN
             - = 2
           * - PISN
             - = 4
           * - PPISN
             - = 8
           * - USSN
             - = 16

   * -
     - (see :ref:`Supernova events/states <supernova-events-states>` for explanation).

   * - Header Strings:
     - Experienced_SN_Type, Experienced_SN_Type(1), Experienced_SN_Type(2), Experienced_SN_Type(SN), Experienced_SN_Type(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **EXPERIENCED_USSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star exploded as an ultra-stripped supernova at any time prior to the current timestep.
   * - Header Strings:
     - Experienced_USSN, Experienced_USSN(1), Experienced_USSN(2), Experienced_USSN(SN), Experienced_USSN(CP)

.. _stellar-props-F:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **FALLBACK_FRACTION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.fallbackFraction
   * - Description:
     - Fallback fraction during a supernova.
   * - Header Strings:
     - Fallback_Fraction, Fallback_Fraction(1), Fallback_Fraction(2), Fallback_Fraction(SN), Fallback_Fraction(CP)

.. _stellar-props-G:

.. _stellar-props-H:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **HE_CORE_MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_HeCoreMass
   * - Description:
     - Helium core mass (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass_He_Core, Mass_He_Core(1), Mass_He_Core(2), Mass_He_Core(SN), Mass_He_Core(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **HE_CORE_MASS_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.HeCoreMass
   * - Description:
     - Helium core mass at the onset of unstable RLOF leading to the CE (\ :math:`M_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Mass_He_Core@CE(1), Mass_He_Core@CE(2), Mass_He_Core@CE(SN), Mass_He_Core@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.HeCoreMassAtCOFormation
   * - Description:
     - Helium core mass immediately prior to a supernova (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass_He_Core@\ CO, Mass_He_Core@CO(1), Mass_He_Core@CO(2), Mass_He_Core@CO(SN), Mass_He_Core@CO(CP)

.. _stellar-props-I:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ID**
     -
   * - Data type:
     - UNSIGNED LONG INT
   * - COMPAS variable:
     - BaseStar::m_ObjectId
   * - Description:
     - Unique object identifier for ``C++`` object – used in debugging to identify objects.
   * - Header Strings:
     - ID, ID(1), ID(2), ID(SN), ID(CP)

`Note that this property has the same header string as BINARY_PROPERTY::ID & BINARY_PROPERTY::RLOF_CURRENT_ID. It is expected that one or 
the other is printed in any file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseStar::m_StellarType
   * - Description:
     - Stellar type at zero age main-sequence (per :cite:`Hurley2000`).
   * - Header Strings:
     - Stellar_Type@\ ZAMS, Stellar_Type@ZAMS(1), Stellar_Type@ZAMS(2), Stellar_Type@ZAMS(SN), Stellar_Type@ZAMS(CP)

`Note that this property has the same header string as INITIAL_STELLAR_TYPE_NAME. It is expected that one or the other is printed in any 
file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **INITIAL_STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseStar::m_StellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) at zero age main-sequence. e.g. "First Giant Branch", "Core Helium Burning", "Helium White Dwarf", etc.
   * - Header Strings:
     - Stellar_Type@\ ZAMS, Stellar_Type@ZAMS(1), Stellar_Type@ZAMS(2), Stellar_Type@ZAMS(SN), Stellar_Type@ZAMS(CP)

`Note that this property has the same header string as INITIAL_STELLAR_TYPE. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_CCSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently a core-collapse supernova.
   * - Header Strings:
     - CCSN, CCSN(1), CCSN(2), CCSN(SN), CCSN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_ECSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently an electron-capture supernova.
   * - Header Strings:
     - ECSN, ECSN(1), ECSN(2), ECSN(SN), ECSN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_HYDROGEN_POOR**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.isHydrogenPoor
   * - Description:
     - Flag to indicate if the star is hydrogen poor.
   * - Header Strings:
     - Is_Hydrogen_Poor, Is_Hydrogen_Poor(1), Is_Hydrogen_Poor(2), Is_Hydrogen_Poor(SN), Is_Hydrogen_Poor(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_PISN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently a pair-instability supernova.
   * - Header Strings:
     - PISN, PISN(1), PISN(2), PISN(SN), PISN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_PPISN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently a pulsational pair-instability supernova.
   * - Header Strings:
     - PPISN, PPISN(1), PPISN(2), PPISN(SN), PPISN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_RLOF**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BinaryConstituentStar::m_RLOFDetails.isRLOF
   * - Description:
     - Flag to indicate whether the star is currently undergoing Roche Lobe overflow.
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - RLOF(1), RLOF(2), RLOF(SN), RLOF(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **IS_USSN**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - Flag to indicate whether the star is currently an ultra-stripped supernova.
   * - Header Strings:
     - USSN, USSN(1), USSN(2), USSN(SN), USSN(CP)

.. _stellar-props-J:

.. _stellar-props-K:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **KICK_MAGNITUDE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.kickMagnitude
   * - Description:
     - Magnitude of natal kick received during a supernova (\ :math:`km s^{−1}`). Calculated using the drawn kick magnitude.
   * - Header Strings:
     - Applied_Kick_Magnitude, Applied_Kick_Magnitude(1), Applied_Kick_Magnitude(2), Applied_Kick_Magnitude(SN), Applied_Kick_Magnitude(CP)

.. _stellar-props-L:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_AT_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.lambda
   * - Description:
     - Common-envelope lambda parameter calculated at the unstable RLOF leading to the CE.
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Lambda@CE(1), Lambda@CE(2), Lambda@CE(SN), Lambda@CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_DEWI**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.dewi
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Dewi2000` using the fit from Appendix A of :doc:`Claeys et al. (2014) <../../references>`.
   * - Header Strings:
     - Dewi, Dewi(1), Dewi(2), Dewi(SN), Dewi(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_FIXED**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.fixed
   * - Description:
     - Universal common envelope lambda parameter specified by the user (program option ``--common-envelope-lambda``).
   * - Header Strings:
     - Lambda_Fixed, Lambda_Fixed(1), Lambda_Fixed(2), Lambda_Fixed(SN), Lambda_Fixed(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_KRUCKOW**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.kruckow
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Kruckow2016` with the alpha exponent set by program option ``--common-envelope-slope-Kruckow``. Spectrum fit to the region bounded by the upper and lower limits as shown in :cite:`Kruckow2016`, Fig. 1.
   * - Header Strings:
     - Kruckow, Kruckow(1), Kruckow(2), Kruckow(SN), Kruckow(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_KRUCKOW_BOTTOM**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.kruckowBottom
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Kruckow2016` with the alpha exponent set to −1. Spectrum fit to the region bounded by the upper and lower limits as shown in :cite:`Kruckow2016`, Fig. 1.
   * - Header Strings:
     - Kruckow_Bottom, Kruckow_Bottom(1), Kruckow_Bottom(2), Kruckow_Bottom(SN), Kruckow_Bottom(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_KRUCKOW_MIDDLE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.kruckowMiddle
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Kruckow2016` with the alpha exponent set to :math:`-\frac{4}{5}`. Spectrum fit to the region bounded by the upper and lower limits as shown in :cite:`Kruckow2016`, Fig. 1.
   * - Header Strings:
     - Kruckow_Middle, Kruckow_Middle(1), Kruckow_Middle(2), Kruckow_Middle(SN), Kruckow_Middle(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_KRUCKOW_TOP**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.kruckowTop
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Kruckow2016` with the alpha exponent set to :math:`-\frac{2}{3}`. Spectrum fit to the region bounded by the upper and lower limits as shown in :cite:`Kruckow2016`, Fig. 1.
   * - Header Strings:
     - Kruckow_Top, Kruckow_Top(1), Kruckow_Top(2), Kruckow_Top(SN), Kruckow_Top(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_LOVERIDGE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.loveridge
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Webbink1984` & :cite:`Loveridge2011`.
   * - Header Strings:
     - Loveridge, Loveridge(1), Loveridge(2), Loveridge(SN), Loveridge(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_LOVERIDGE_WINDS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.loveridgeWinds
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :cite:`Webbink1984` & :cite:`Loveridge2011` including winds.
   * - Header Strings:
     - Loveridge_Winds, Loveridge_Winds(1), Loveridge_Winds(2), Loveridge_Winds(SN), Loveridge_Winds(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LAMBDA_NANJING**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Lambdas.nanjing
   * - Description:
     - Envelope binding energy parameter lambda calculated as per :doc:`Xu & Li (2010) <../../references>`.
   * - Header Strings:
     - Lambda_Nanjing, Lambda_Nanjing(1), Lambda_Nanjing(2), Lambda_Nanjing(SN), Lambda_Nanjing(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LBV_PHASE_FLAG**
     -
   * - Data type:
     - BOOL
   * - COMPAS variable:
     - BaseStar::m_LBVphaseFlag
   * - Description:
     - Flag to indicate if the star ever entered the luminous blue variable phase.
   * - Header Strings:
     - LBV_Phase_Flag, LBV_Phase_Flag(1), LBV_Phase_Flag(2), LBV_Phase_Flag(SN), LBV_Phase_Flag(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LUMINOSITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Luminosity
   * - Description:
     - Luminosity (\ :math:`L_\odot`).
   * - Header Strings:
     - Luminosity, Luminosity(1), Luminosity(2), Luminosity(SN), Luminosity(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LUMINOSITY_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.luminosity
   * - Description:
     - Luminosity immediately following common envelope event (\ :math:`L_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Luminosity>CE(1), Luminosity>CE(2), Luminosity>CE(SN), Luminosity>CE(CP)


.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **LUMINOSITY_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.luminosity
   * - Description:
     - Luminosity at the onset of unstable RLOF leading to the CE (\ :math:`L_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Luminosity<CE(1), Luminosity<CE(2), Luminosity<CE(SN), Luminosity<CE(CP)

.. _stellar-props-M:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Mass
   * - Description:
     - Mass (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass, Mass(1), Mass(2), Mass(SN), Mass(CP)

`Note that this property has the same header string as RLOF_CURRENT_STAR1_MASS & RLOF_CURRENT_STAR2_MASS. It is expected that one or 
the other is printed in any file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_0**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Mass0
   * - Description:
     - Effective initial mass (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass_0, Mass_0(1), Mass_0(2), Mass_0(SN), Mass_0(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_LOSS_DIFF**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_MassLossDiff
   * - Description:
     - The amount of mass lost due to winds (\ :math:`M_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - dmWinds(1), dmWinds(2), dmWinds(SN), dmWinds(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MASS_TRANSFER_DIFF**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_MassTransferDiff
   * - Description:
     - The amount of mass accreted or donated during a mass transfer episode (\ :math:`M_\odot`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - dmMT(1), dmMT(2), dmMT(SN), dmMT(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MDOT**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Mdot
   * - Description:
     - Mass loss rate (\ :math:`M_\odot yr^{−1}`).
   * - Header Strings:
     - Mdot, Mdot(1), Mdot(2), Mdot(SN), Mdot(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MEAN_ANOMALY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.meanAnomaly
   * - Description:
     - Mean anomaly of supernova kick. Supplied by user in ``grid`` file, default = random number drawn from [0..2π).
   * -
     - See https://en.wikipedia.org/wiki/Mean_anomaly for explanation.
   * - Header Strings:
     - SN Kick Mean Anomaly, SN_Kick_Mean_Anomaly(1), SN_Kick_Mean_Anomaly(2), SN_Kick_Mean_Anomaly(SN), SN_Kick_Mean_Anomaly(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **METALLICITY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Metallicity
   * - Description:
     - ZAMS Metallicity.
   * - Header Strings:
     - Metallicity@\ ZAMS, Metallicity@ZAMS(1), Metallicity@ZAMS(2), Metallicity@ZAMS(SN), Metallicity@ZAMS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MT_DONOR_HIST**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - BaseStar::m_MassTransferDonorHistory
   * - Description:
     - A list of all of the stellar types from which the current star was a Mass Transfer donor. This can be readily converted into the different cases of Mass Transfer, depending on the working definition. The output string is formatted as #-#-#... where each # represents a stellar type, in chronological order. E.g, 2-8 means the star was an MT donor as a ``HG`` (stellar type 2) star, and later as a ``HeHG`` (stellar type 8) star.
   * - Header Strings:
     - MT_Donor_Hist, MT_Donor_Hist(1), MT_Donor_Hist(2), MT_Donor_Hist(SN), MT_Donor_Hist(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **MZAMS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_MZAMS
   * - Description:
     - ZAMS Mass (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass@\ ZAMS, Mass@ZAMS(1), Mass@ZAMS(2), Mass@ZAMS(SN), Mass@ZAMS(CP)

.. _stellar-props-N:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NUCLEAR_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_NuclearTimescale
   * - Description:
     - Nuclear timescale (Myr).
   * - Header Strings:
     - Tau_Nuclear, Tau_Nuclear(1), Tau_Nuclear(2), Tau_Nuclear(SN), Tau_Nuclear(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NUCLEAR_TIMESCALE_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.nuclearTimescale
   * - Description:
     - Nuclear timescale immediately following common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau Nuclear>CE(1), Tau Nuclear>CE(2), Tau Nuclear>CE(SN), Tau Nuclear>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.nuclearTimescale
   * - Description:
     - Nuclear timescale at the onset of unstable RLOF leading to the CE (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Nuclear<CE(1), Tau_Nuclear<CE(2), Tau_Nuclear<CE(SN), Tau_Nuclear<CE(CP)

.. _stellar-props-O:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **OMEGA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Omega
   * - Description:
     - Angular frequency (\ :math:`yr^{−1}`).
   * - Header Strings:
     - Omega, Omega(1), Omega(2), Omega(SN), Omega(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **OMEGA_BREAK**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_OmegaBreak
   * - Description:
     - Break-up angular frequency (\ :math:`yr^{−1}`).
   * - Header Strings:
     - Omega_Break, Omega_Break(1), Omega_Break(2), Omega_Break(SN), Omega_Break(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **OMEGA_ZAMS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_OmegaZAMS
   * - Description:
     - Angular frequency at ZAMS (\ :math:`yr^{−1}`).
   * - Header Strings:
     - Omega@\ ZAMS, Omega@ZAMS(1), Omega@ZAMS(2), Omega@ZAMS(SN), Omega@ZAMS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_ENERGY_POST_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_PostSNeOrbitalEnergy
   * - Description:
     - Absolute value of orbital energy immediately following supernova event (\ :math:`M_\odot AU^2 yr^{−2}`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Orbital_Energy>SN(1), Orbital_Energy>SN(2), Orbital_Energy>SN(SN), Orbital_Energy>SN(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ORBITAL_ENERGY_PRE_SUPERNOVA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_PreSNeOrbitalEnergy
   * - Description:
     - Orbital energy immediately prior to supernova event (\ :math:`M_\odot AU^2 yr^{−2}`).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Orbital_Energy<SN(1), Orbital_Energy<SN(2), Orbital_Energy<SN(SN), Orbital_Energy<SN(CP)

.. _stellar-props-P:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_MAGNETIC_FIELD**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_PulsarDetails.magneticField
   * - Description:
     - Pulsar magnetic field strength (G).
   * - Header Strings:
     - Pulsar_Mag_Field, Pulsar_Mag_Field(1), Pulsar_Mag_Field(2), Pulsar_Mag_Field(SN), Pulsar_Mag_Field(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_SPIN_DOWN_RATE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_PulsarDetails.spinDownRate
   * - Description:
     - Pulsar spin-down rate.
   * - Header Strings:
     - Pulsar_Spin_Down, Pulsar_Spin_Down(1), Pulsar_Spin_Down(2), Pulsar_Spin_Down(SN), Pulsar_Spin_Down(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_SPIN_FREQUENCY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_PulsarDetails.spinFrequency
   * - Description:
     - Pulsar spin angular frequency (\ :math:`rads s^{−1}`).
   * - Header Strings:
     - Pulsar_Spin_Freq, Pulsar_Spin_Freq(1), Pulsar_Spin_Freq(2), Pulsar_Spin_Freq(SN), Pulsar_Spin_Freq(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **PULSAR_SPIN_PERIOD**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_PulsarDetails.spinPeriod
   * - Description:
     - Pulsar spin period (ms).
   * - Header Strings:
     - Pulsar_Spin_Period, Pulsar_Spin_Period(1), Pulsar_Spin_Period(2), Pulsar_Spin_Period(SN), Pulsar_Spin_Period(CP)

.. _stellar-props-Q:

.. _stellar-props-R:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIAL_EXPANSION_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_RadialExpansionTimescale
   * - Description:
     - e-folding time of stellar radius (Myr).
   * - Header Strings:
     - Tau_Radial, Tau_Radial(1), Tau_Radial(2), Tau_Radial(SN), Tau_Radial(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIAL_EXPANSION_TIMESCALE_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.radialExpansionTimescale
   * - Description:
     - e-folding time of stellar radius immediately following common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Radial<CE(1), Tau_Radial<CE(2), Tau_Radial<CE(SN), Tau_Radial<CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.radialExpansionTimescale
   * - Description:
     - e-folding time of stellar radius at the onset of unstable RLOF leading to the CE (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Radial<CE(1), Tau_Radial<CE(2), Tau_Radial<CE(SN), Tau_Radial<CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RADIUS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Radius
   * - Description:
     - Radius (\ :math:`R_\odot`).
   * - Header Strings:
     - Radius, Radius(1), Radius(2), Radius(SN), Radius(CP)

`Note that this property has the same header string as RLOF_CURRENT_STAR1_RADIUS & RLOF_CURRENT_STAR2_RADIUS. It is expected that one or 
the other is printed in any file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RANDOM_SEED**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_RandomSeed
   * - Description:
     - Seed for random number generator for this star.
   * - Header Strings:
     - SEED, SEED(1), SEED(2), SEED(SN), SEED(CP)

`Note that this property has the same header string as BINARY_PROPERTY::RANDOM_SEED & BINARY_PROPERTY::RLOF_CURRENT_RANDOM_SEED. It is 
expected that one or the other is printed in any file, but not both. If both are printed then the file will contain two columns with the 
same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RECYCLED_NEUTRON_STAR**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the object was a recycled neutron star at any time prior to the current timestep (was a neutron star accreting mass).
   * - Header Strings:
     - Recycled_NS, Recycled_NS(1), Recycled_NS(2), Recycled_NS(SN), Recycled_NS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RLOF_ONTO_NS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star transferred mass to a neutron star at any time prior to the current timestep.
   * - Header Strings:
     - RLOF->NS, RLOF->NS(1), RLOF->NS(2), RLOF->NS(SN), RLOF->NS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RUNAWAY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.past
   * - Description:
     - Flag to indicate whether the star was unbound by a supernova event at any time prior to the current timestep. (i.e Unbound after supernova event and not a WD, NS, BH or MR).
   * - Header Strings:
     - Runaway, Runaway(1), Runaway(2), Runaway(SN), Runaway(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **RZAMS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_RZAMS
   * - Description:
     - ZAMS Radius (\ :math:`R_\odot`).
   * - Header Strings:
     - R@\ ZAMS, R@ZAMS(1), R@ZAMS(2), R@ZAMS(SN), R@ZAMS(CP)

.. _stellar-props-S:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SN_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` BaseStar::m_SupernovaDetails.events.current
   * - Description:
     - The type of supernova event currently being experienced by the star. Printed as one of

        .. list-table::
           :widths: 20 10
           :header-rows: 0
           :class: aligned-text
    
           * - NONE
             - = 0
           * - CCSN
             - = 1
           * - ECSN
             - = 2
           * - PISN
             - = 4
           * - PPISN
             - = 8
           * - USSN
             - = 16
   * -
     - (see :ref:`Supernova events/states <supernova-events-states>` for explanation).
   * - Header Strings:
     - SN_Type, SN_Type(1), SN_Type(2), SN_Type(SN), SN_Type(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseStar::m_StellarType
   * - Description:
     - Stellar type (per :cite:`Hurley2000`).
   * - Header Strings:
     - Stellar_Type, Stellar_Type(1), Stellar_Type(2), Stellar_Type(SN), Stellar_Type(CP)

`Note that this property has the same header string as STELLAR_TYPE_NAME. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseStar::m_StellarType
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`). e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header Strings:
     - Stellar_Type, Stellar_Type(1), Stellar_Type(2), Stellar_Type(SN), Stellar_Type(CP)

`Note that this property has the same header string as STELLAR_TYPE. It is expected that one or the other is printed in any file, but 
not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_PREV**
     -
   * - Data type:
     - INT
   * - COMPAS variable:
     - BaseStar::m_StellarTypePrev
   * - Description:
     - Stellar type (per :cite:`Hurley2000`) at previous timestep.
   * - Header Strings:
     - Stellar_Type_Prev, Stellar_Type_Prev(1), Stellar_Type_Prev(2), Stellar_Type_Prev(SN), Stellar_Type_Prev(CP)

`Note that this property has the same header string as STELLAR_TYPE_PREV_NAME. It is expected that one or the other is printed in any 
file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **STELLAR_TYPE_PREV_NAME**
     -
   * - Data type:
     - STRING
   * - COMPAS variable:
     - `derived from` BaseStar::m_StellarTypePrev
   * - Description:
     - Stellar type name (per :cite:`Hurley2000`) at previous timestep. e.g. "First_Giant_Branch", "Core_Helium_Burning", "Helium_White_Dwarf", etc.
   * - Header Strings:
     - Stellar_Type_Prev, Stellar_Type_Prev(1), Stellar_Type_Prev(2), Stellar_Type_Prev(SN), Stellar_Type_Prev(CP)

`Note that this property has the same header string as STELLAR_TYPE_PREV. It is expected that one or the other is printed in any file, 
but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SUPERNOVA_KICK_MAGNITUDE_MAGNITUDE_RANDOM_NUMBER**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.kickMagnitudeRandom
   * - Description:
     - Random number for drawing the supernova kick magnitude (if required). Supplied by user in grid file, default = random number drawn from [0..1).
   * - Header Strings:
     - SN_Kick_Magnitude_Random_Number, SN_Kick_Magnitude_Random_Number(1), SN_Kick_Magnitude_Random_Number(2), SN_Kick_Magnitude_Random_Number(SN), SN_Kick_Magnitude_Random_Number(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SUPERNOVA_PHI**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.phi
   * - Description:
     - Angle between 'x' and 'y', both in the orbital plane of supernovae vector (rad). Supplied by user in grid file, default = random number drawn from [0..2π).
   * - Header Strings:
     - SN_Kick_Phi, SN_Kick_Phi(1), SN_Kick_Phi(2), SN_Kick_Phi(SN), SN_Kick_Phi(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **SUPERNOVA_THETA**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.theta
   * - Description:
     - Angle between the orbital plane and the ’z’ axis of supernovae vector (rad). Supplied by user in grid file, default = drawn from distribution specified by program option ``--kick direction``.
   * - Header Strings:
     - SN_Kick_Theta, SN_Kick_Theta(1), SN_Kick_Theta(2), SN_Kick_Theta(SN), SN_Kick_Theta(CP)

.. _stellar-props-T:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TEMPERATURE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Temperature
   * - Description:
     - Effective temperature (K).
   * - Header Strings:
     - Teff, Teff(1), Teff(2), Teff(SN), Teff(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TEMPERATURE_POST_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.temperature
   * - Description:
     - Effective temperature immediately following common envelope event (K).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Teff>CE(1), Teff>CE(2), Teff>CE(SN), Teff>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TEMPERATURE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.temperature
   * - Description:
     - Effective temperature at the unstable RLOF leading to the CE (K).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Teff<CE(1), Teff<CE(2), Teff<CE(SN), Teff<CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **THERMAL_TIMESCALE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_ThermalTimescale
   * - Description: 
     - Thermal timescale (Myr).
   * - Header Strings:
     - Tau_Thermal, Tau_Thermal(1), Tau_Thermal(2), Tau_Thermal(SN), Tau_Thermal(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **THERMAL_TIMESCALE_POST_COMMON ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.postCEE.thermalTimescale
   * - Description:
     - Thermal timescale immediately following common envelope event (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Thermal>CE(1), Tau_Thermal>CE(2), Tau_Thermal>CE(SN), Tau_Thermal>CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BinaryConstituentStar::m_CEDetails.preCEE.thermalTimescale
   * - Description:
     - Thermal timescale at the onset of the unstable RLOF leading to the CE (Myr).
   * -
     - `Applies only to constituent stars of a binary system (i.e. does not apply to` ``SSE``\ `).`
   * - Header Strings:
     - Tau_Thermal<CE(1), Tau_Thermal<CE(2), Tau_Thermal<CE(SN), Tau_Thermal<CE(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TIME**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Time
   * - Description:
     - Time since ZAMS (Myr).
   * - Header Strings:
     - Time, Time(1), Time(2), Time(SN), Time(CP)

`Note that this property has the same header string as BINARY_PROPERTY::TIME & BINARY_PROPERTY::RLOF_CURRENT_TIME. It is expected that one 
or the other is printed in any file, but not both. If both are printed then the file will contain two columns with the same header string.`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TIMESCALE_MS**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Timescales[tMS]
   * - Description:
     - Main Sequence timescale (Myr).
   * - Header Strings: tMS, tMS(1), tMS(2), tMS(SN), tMS(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.totalMassAtCOFormation
   * - Description:
     - Total mass of the star at the beginning of a supernova event (\ :math:`M_\odot`).
   * - Header Strings:
     - Mass_Total@\ CO, Mass_Total@CO(1), Mass_Total@CO(2), Mass_Total@CO(SN), Mass_Total@CO(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **TRUE_ANOMALY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_SupernovaDetails.trueAnomaly
   * - Description:
     - True anomaly calculated using Kepler’s equation (rad).
   * -
     - See https://en.wikipedia.org/wiki/True anomaly for explanation.
   * - Header Strings:
     - True_Anomaly(psi), True_Anomaly(psi)(1), True_Anomaly(psi)(2), True_Anomaly(psi)(SN), True_Anomaly(psi)(CP)

.. _stellar-props-U:

.. _stellar-props-V:

.. _stellar-props-W:

.. _stellar-props-X:

.. _stellar-props-Y:

.. _stellar-props-Z:

:ref:`Back to Top <stellar-props-top>`

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_HURLEY**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Zetas.hurley
   * - Description:
     - Adiabatic exponent calculated per :cite:`Hurley2000` using core mass.
   * - Header Strings:
     - Zeta_Hurley, Zeta_Hurley(1), Zeta_Hurley(2), Zeta_Hurley(SN), Zeta_Hurley(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_HURLEY_HE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Zetas.hurleyHe
   * - Description:
     - Adiabatic exponent calculated per :cite:`Hurley2000` using He core mass.
   * - Header Strings:
     - Zeta_Hurley_He, Zeta_Hurley_He(1), Zeta_Hurley_He(2), Zeta_Hurley_He(SN), Zeta_Hurley_He(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_SOBERMAN**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Zetas.soberman
   * - Description:
     - Adiabatic exponent calculated per :cite:`Soberman1997` using core mass.
   * - Header Strings:
     - Zeta_Soberman, Zeta_Soberman(1), Zeta_Soberman(2), Zeta_Soberman(SN), Zeta_Soberman(CP)

.. flat-table::
   :widths: 25 75 1 1
   :header-rows: 0
   :class: aligned-text

   * - :cspan:`2` **ZETA_SOBERMAN_HE**
     -
   * - Data type:
     - DOUBLE
   * - COMPAS variable:
     - BaseStar::m_Zetas.sobermanHe
   * - Description:
     - Adiabatic exponent calculated per :cite:`Soberman1997` using He core mass.
   * - Header Strings:
     - Zeta_Soberman_He, Zeta_Soberman_He(1), Zeta_Soberman_He(2), Zeta_Soberman_He(SN), Zeta_Soberman_He(CP)

.. _supernova-events-states:

:ref:`Back to Top <stellar-props-top>`

Supernova events/states
-----------------------

Supernova events/states, both current ("is") and past ("experienced"), are stored within COMPAS as bitmaps. That means different values 
can be ORed or ANDed into the bit map, so that various events or states can be set concurrently.

The values shown below for the ``SN_EVENT`` type are powers of 2 so that they can be used in a bitmap and manipulated with bit-wise logical 
operators. Any of the individual supernova event/state types that make up the ``SN_EVENT`` type can be set independently of any other event/state.

`constants.h` defines an enum class for ``SN_EVENT``, and an associated label map, ``SN_EVENT_LABEL``, to provide labels for the events.  These 
are shown below::

    enum class SN_EVENT: int {
        NONE         = 0,
        CCSN         = 1,
        ECSN         = 2,
        PISN         = 4,
        PPISN        = 8,
        USSN         = 16,
        RUNAWAY      = 32,
        RECYCLED NS  = 64,
        RLOF ONTO NS = 128
    };


    const COMPASUnorderedMap<SN_EVENT, std::string> SN_EVENT_LABEL = {
        { SN EVENT::NONE,         "No Supernova" },
        { SN EVENT::CCSN,         "Core Collapse Supernova" },
        { SN EVENT::ECSN,         "Electron Capture Supernova" },
        { SN EVENT::PISN,         "Pair Instability Supernova" },
        { SN EVENT::PPISN,        "Pulsational Pair Instability Supernova" },
        { SN EVENT::USSN,         "Ultra Stripped Supernova" },
        { SN EVENT::RUNAWAY,      "Runaway Companion" },
        { SN EVENT::RECYCLED NS,  "Recycled Neutron Star" },
        { SN EVENT::RLOF ONTO NS, "Donated Mass to Neutron Star through RLOF" }
    };

A convenience function (shown below) is provided in ``utils.cpp`` to interpret the bit map.

::

    /*
    * Returns a single SN type based on the SN EVENT parameter passed
    *
    * Returns (in priority order):
    *
    * SN EVENT::CCSN iff CCSN bit is set and USSN bit is not set
    * SN EVENT::ECSN iff ECSN bit is set
    * SN EVENT::PISN iff PISN bit is set
    * SN EVENT::PPISN iff PPISN bit is set
    * SN EVENT::USSN iff USSN bit is set
    * SN EVENT::NONE otherwise
    *
    *
    * @param [IN] p SNEvent SN EVENT mask to check for SN event type
    * @return SN EVENT
    /
    SN EVENT SNEventType(const SN EVENT p SNEvent) {
        if ((p SNEvent & (SN EVENT::CCSN | SN EVENT::USSN)) == SN EVENT::CCSN ) return SN EVENT::CCSN;
        if ((p SNEvent & SN EVENT::ECSN )                   == SN EVENT::ECSN ) return SN EVENT::ECSN;
        if ((p SNEvent & SN EVENT::PISN )                   == SN EVENT::PISN ) return SN EVENT::PISN;
        if ((p SNEvent & SN EVENT::PPISN)                   == SN EVENT::PPISN) return SN EVENT::PPISN;
        if ((p SNEvent & SN EVENT::USSN )                   == SN EVENT::USSN ) return SN EVENT::USSN;

        return SN EVENT::NONE;
    }

BSE common envelopes
====================

Default record definition for the BSE Common Envelopes log file::

    const ANY_PROPERTY_VECTOR BSE_COMMON_ENVELOPES_REC = {
        BINARY_PROPERTY::RANDOM_SEED,
        BINARY_PROPERTY::TIME,
        STAR_1_PROPERTY::LAMBDA_AT_COMMON_ENVELOPE,
        STAR_2_PROPERTY::LAMBDA_AT_COMMON_ENVELOPE,
        STAR_1_PROPERTY::BINDING_ENERGY_PRE_COMMON_ENVELOPE,
        STAR_2_PROPERTY::BINDING_ENERGY_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::ECCENTRICITY_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::ECCENTRICITY_POST_COMMON_ENVELOPE,
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_POST_COMMON_ENVELOPE,
        BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::ROCHE_LOBE_RADIUS_1_POST_COMMON_ENVELOPE,
        BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::ROCHE_LOBE_RADIUS_2_POST_COMMON_ENVELOPE,
        BINARY_PROPERTY::MASS_1_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::MASS_1_POST_COMMON_ENVELOPE,
        BINARY_PROPERTY::MASS_ENV_1,
        BINARY_PROPERTY::RADIUS_1_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::RADIUS_1_POST_COMMON_ENVELOPE,
        BINARY_PROPERTY::STELLAR_TYPE_1_PRE_COMMON_ENVELOPE,
        STAR_1_PROPERTY::STELLAR_TYPE,
        STAR_1_PROPERTY::LAMBDA_FIXED,
        STAR_1_PROPERTY::LAMBDA_NANJING,
        STAR_1_PROPERTY::LAMBDA_LOVERIDGE,
        STAR_1_PROPERTY::LAMBDA_LOVERIDGE_WINDS,
        STAR_1_PROPERTY::LAMBDA_KRUCKOW,
        STAR_1_PROPERTY::BINDING_ENERGY_FIXED,
        STAR_1_PROPERTY::BINDING_ENERGY_NANJING,
        STAR_1_PROPERTY::BINDING_ENERGY_LOVERIDGE,
        STAR_1_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS,
        STAR_1_PROPERTY::BINDING_ENERGY_KRUCKOW,
        BINARY_PROPERTY::MASS_2_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::MASS_2_POST_COMMON_ENVELOPE,
        BINARY_PROPERTY::MASS_ENV_2,
        BINARY_PROPERTY::RADIUS_2_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::RADIUS_2_POST_COMMON_ENVELOPE,
        BINARY_PROPERTY::STELLAR_TYPE_2_PRE_COMMON_ENVELOPE,
        STAR_2_PROPERTY::STELLAR_TYPE,
        STAR_2_PROPERTY::LAMBDA_FIXED,
        STAR_2_PROPERTY::LAMBDA_NANJING,
        STAR_2_PROPERTY::LAMBDA_LOVERIDGE,
        STAR_2_PROPERTY::LAMBDA_LOVERIDGE_WINDS,
        STAR_2_PROPERTY::LAMBDA_KRUCKOW,
        STAR_2_PROPERTY::BINDING_ENERGY_FIXED,
        STAR_2_PROPERTY::BINDING_ENERGY_NANJING,
        STAR_2_PROPERTY::BINDING_ENERGY_LOVERIDGE,
        STAR_2_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS,
        STAR_2_PROPERTY::BINDING_ENERGY_KRUCKOW,
        BINARY_PROPERTY::MASS_TRANSFER_TRACKER_HISTORY,
        BINARY_PROPERTY::STELLAR_MERGER,
        BINARY_PROPERTY::OPTIMISTIC_COMMON_ENVELOPE,
        BINARY_PROPERTY::COMMON_ENVELOPE_EVENT_COUNT,
        BINARY_PROPERTY::DOUBLE_CORE_COMMON_ENVELOPE,
        STAR_1_PROPERTY::IS_RLOF,
        STAR_1_PROPERTY::LUMINOSITY_PRE_COMMON_ENVELOPE,
        STAR_1_PROPERTY::TEMPERATURE_PRE_COMMON_ENVELOPE,
        STAR_1_PROPERTY::DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE,
        STAR_1_PROPERTY::THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE,
        STAR_1_PROPERTY::NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE,
        STAR_2_PROPERTY::IS_RLOF,
        STAR_2_PROPERTY::LUMINOSITY_PRE_COMMON_ENVELOPE,
        STAR_2_PROPERTY::TEMPERATURE_PRE_COMMON_ENVELOPE,
        STAR_2_PROPERTY::DYNAMICAL_TIMESCALE_PRE_COMMON_ENVELOPE,
        STAR_2_PROPERTY::THERMAL_TIMESCALE_PRE_COMMON_ENVELOPE,
        STAR_2_PROPERTY::NUCLEAR_TIMESCALE_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::ZETA_STAR,
        BINARY_PROPERTY::ZETA_LOBE,
        BINARY_PROPERTY::SYNCHRONIZATION_TIMESCALE,
        BINARY_PROPERTY::CIRCULARIZATION_TIMESCALE,
        STAR_1_PROPERTY::RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE,
        STAR_2_PROPERTY::RADIAL_EXPANSION_TIMESCALE_PRE_COMMON_ENVELOPE,
        BINARY_PROPERTY::IMMEDIATE_RLOF_POST_COMMON_ENVELOPE,
        BINARY_PROPERTY::SIMULTANEOUS_RLOF
    };

BSE supernovae
==============

Default record definition for the BSE Supernovae log file::

    const ANY_PROPERTY_VECTOR BSE_SUPERNOVAE_REC = {
        BINARY_PROPERTY::RANDOM_SEED,
        SUPERNOVA_PROPERTY::DRAWN_KICK_MAGNITUDE,
        SUPERNOVA_PROPERTY::KICK_MAGNITUDE,
        SUPERNOVA_PROPERTY::FALLBACK_FRACTION,
        BINARY_PROPERTY::ORBITAL_VELOCITY_PRE_SUPERNOVA,
        BINARY_PROPERTY::DIMENSIONLESS_KICK_MAGNITUDE, 
        SUPERNOVA_PROPERTY::MEAN_ANOMALY,
        SUPERNOVA_PROPERTY::SUPERNOVA_THETA,
        SUPERNOVA_PROPERTY::SUPERNOVA_PHI,
        SUPERNOVA_PROPERTY::SN_TYPE,
        BINARY_PROPERTY::ECCENTRICITY_PRE_SUPERNOVA,  
        BINARY_PROPERTY::ECCENTRICITY,
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_PRE_SUPERNOVA_RSOL,
        BINARY_PROPERTY::SEMI_MAJOR_AXIS_RSOL,
        BINARY_PROPERTY::TIME,
        BINARY_PROPERTY::SUPERNOVA_STATE,
        BINARY_PROPERTY::UNBOUND,
        COMPANION_PROPERTY::STELLAR_TYPE,
        SUPERNOVA_PROPERTY::STELLAR_TYPE,
        SUPERNOVA_PROPERTY::STELLAR_TYPE_PREV,
        COMPANION_PROPERTY::MASS,
        SUPERNOVA_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION,
        SUPERNOVA_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
        SUPERNOVA_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
        SUPERNOVA_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION,
        SUPERNOVA_PROPERTY::MASS,
        SUPERNOVA_PROPERTY::EXPERIENCED_RLOF,
        SUPERNOVA_PROPERTY::MASS_TRANSFER_DONOR_HISTORY,
        SUPERNOVA_PROPERTY::SPEED,
        COMPANION_PROPERTY::SPEED,
        BINARY_PROPERTY::SYSTEMIC_SPEED,
        SUPERNOVA_PROPERTY::IS_HYDROGEN_POOR
    };

Default log file record specifications
======================================

Following are the default log file record specifications for each of the standard log files. These specifications
can be overridden by the use of a log file specifications file via the ``logfile-definitions`` program option.

For Single Star Evolution (SSE):

.. toctree::
   :maxdepth: 1

   standard-logfiles-default-record-specifications-SSE-sysparms
   standard-logfiles-default-record-specifications-SSE-detailed
   standard-logfiles-default-record-specifications-SSE-supernovae
   standard-logfiles-default-record-specifications-SSE-switchlog

For Binary Star Evolution (BSE):

.. toctree::
   :maxdepth: 1

   standard-logfiles-default-record-specifications-BSE-sysparms
   standard-logfiles-default-record-specifications-BSE-detailed
   standard-logfiles-default-record-specifications-BSE-supernovae
   standard-logfiles-default-record-specifications-BSE-dco
   standard-logfiles-default-record-specifications-BSE-ce
   standard-logfiles-default-record-specifications-BSE-pulsars
   standard-logfiles-default-record-specifications-BSE-rlof
   standard-logfiles-default-record-specifications-BSE-switchlog
SSE switchlog
=============

Default record definition for the SSE SwitchLog log file::

    const ANY_PROPERTY_VECTOR SSE_SWITCH_LOG_REC = {
        STAR_PROPERTY::RANDOM_SEED,
        STAR_PROPERTY::TIME
    };


The default record specification can be modified at runtime via a logfile record specifications file (program option ``--logfile-definitions``).
See :doc:`standard-logfiles-record-specification` for details.

Note that the SSE SwitchLog file has the following columns automatically appended to each record:

    - The stellar type from which the star is switching.
    - The stellar type to which the star is switching.

|br|
**SWITCHING_FROM**

.. list-table::
   :widths: 20 80 
   :header-rows: 0
   :class: aligned-text

   * - Data type:
     - INT
   * - COMPAS variable:
     - `derived from` ``BaseStar::m_StellarType``
   * - Description:
     - The stellar type of the star immediately prior to the switch.
   * - Header String:
     - "SWITCHING_FROM"

**SWITCHING_TO**

.. list-table::
   :widths: 20 80 
   :header-rows: 0
   :class: aligned-text

   * - Data type:
     - INT
   * - COMPAS variable:
     - Not applicable
   * - Description:
     - The stellar type to which the star will switch (i.e. the stellar type immediately following the switch).
   * - Header String:
     - "SWITCHING_TO"

These columns will always be automatically appended to each SSE Switch Log record: they cannot be removed via the log file record 
specifications file.
SSE detailed output
===================

Default record definition for the SSE Detailed Output log file::

    const ANY_PROPERTY_VECTOR SSE_DETAILED_OUTPUT_REC = {
        STAR_PROPERTY::AGE,
        STAR_PROPERTY::DT,
        STAR_PROPERTY::TIME,
        STAR_PROPERTY::STELLAR_TYPE,
        STAR_PROPERTY::METALLICITY,
        STAR_PROPERTY::MASS_0,
        STAR_PROPERTY::MASS,
        STAR_PROPERTY::RADIUS,
        STAR_PROPERTY::RZAMS,
        STAR_PROPERTY::LUMINOSITY,
        STAR_PROPERTY::TEMPERATURE,
        STAR_PROPERTY::CORE_MASS,
        STAR_PROPERTY::CO_CORE_MASS,
        STAR_PROPERTY::HE_CORE_MASS,
        STAR_PROPERTY::MDOT,
        STAR_PROPERTY::DOMINANT_MASS_LOSS_RATE,
        STAR_PROPERTY::TIMESCALE_MS
    };

Program options
===============

COMPAS provides a rich set of configuration parameters via program options, allowing users to vary many
parameters that affect the evolution of single and binary stars, and the composition of the population of
stars being evolved. Furthermore, COMPAS allows some parameters to be specified as ranges, or sets, of
values via the program options, allowing users to specify a grid of parameter values on the command line.

Combining command-line program options ranges and sets with a :doc:`grid file <../grid-files>` allows users
more flexibility and the ability to specify more complex combinations of parameter values.

Not all program options can be specified as ranges or sets of values. Options for which mixing different
values in a single execution of COMPAS would either not be meaningful, or might cause undesirable results,
such as options that specify the mode of evolution (e.g. ``--mode``), or the name or path of output files
(e.g. ``--output-path``, ``--logfile-detailed-output`` etc.), cannot be specified as a range or set of values.
COMPAS will issue an error message if ranges or sets are specified for options for which they are not supported.

.. toctree::
   :maxdepth: 1

   program-options-ranges
   program-options-sets
   program-options-mixing-ranges-sets
   program-options-vector-options
   program-options-list-defaults

See :doc:`./program-options-list-defaults` for a full list of available program options and their default valaues.

Program option ranges
=====================

A range of values can be specified for any numeric options (i.e. integer (or integer variant), and floating
point (or floating point variant) data types) that are not excluded from range specifications.

Option value ranges are specified by:

    --option-name range-specifier

where `range-specifier` is defined as:

    range-identifier[start,count,increment]

    .. list-table::
       :widths: 19 81 
       :header-rows: 0
       :class: aligned-text

       * - `range-identifier`
         - is one of {’r’, ’range’} (case is not significant)
       * - `start`
         - is the starting value of the range
       * - `count`
         - is the number of values in the range (must be an unsigned long int)
       * - `increment`
         - is the amount by which the value increments for each value in the range

    Note that:

        `range-identifier` is optional for range-specifier. |br|
        `start` and `increment` must be the same data type as option-name. |br|
        `count` must be a positive integer value.

There should be no spaces inside the brackets ([]). Spaces on the command line are interpreted as argument delimiters
by the shell parser before passing the command-line arguments to the COMPAS executable, so if spaces are present inside
the brackets the shell parser breaks the range specification into multiple command-line arguments.

To specify a range of values for the ``--metallicity`` option, a user, if running COMPAS from the command line
and with no grid file, would type any of the following::

    ./COMPAS --metallicity [0.0001,5,0.0013]

    ./COMPAS --metallicity r[0.0001,5,0.0013]
    
    ./COMPAS --metallicity range[0.0001,5,0.0013]

In each of the examples above the user has specified, by the use of the `range-specifier`, that five binary stars
should be evolved, with constituent star metallicities = 0.0001, 0.0014, 0.0027, 0.0040, and 0.0053.

To evolve a grid of binaries with ten different metallicities, starting at 0.0001 and incrementing by 0.0002,
and five different common envelope alpha values, starting at 0.1 and incrementing by 0.2, the user would
type::

    ./COMPAS --metallicity [0.0001,10,0.0013] --common-envelope-alpha [0.1,5,0.2]

and COMPAS would evolve a grid of 50 binaries using the 10 metallicity values and 5 common envelope alpha values.

Note that when a range is, or ranges are, specified on the command line, the ``--number-of-systems`` command-line option is ignored.
This is to avoid multiple systems with identical initial values being evolved.  Ranges and sets can be mixed with grid files, and
in that case ranges and sets specified on the command line will be played out for each grid file line.

   Program option list and default values
======================================

Any program options that are not specified take default values.

- On the command line, program options that are not explicitly specified default to the COMPAS default value for the option (as specified in the COMPAS code - may be sampled from a distribution).

- On a :doc:`grid file <../grid-files>` line, program options that are not explicitly specified default to the value specified for that option on the command line. If the program option was not explicitly specified on the command line, it will default to the COMPAS default value for the option, as described above. That is, the value for any option not specified on a grid file line option falls back to the value specified on the command line, which falls back to the COMPAS default if it was not specified on the command line.


.. _options-props-top:

The full list of program options with brief explanations and their default values is shown below.  We also include a listing of options (this time, by name only) grouped by category.

**Alphabetical listing: jump to**
:ref:`A <options-props-A>` :ref:`B <options-props-B>` :ref:`C <options-props-C>` :ref:`D <options-props-D>`
:ref:`E <options-props-E>` :ref:`F <options-props-F>` :ref:`G <options-props-G>` :ref:`H <options-props-H>`
:ref:`I <options-props-I>` :ref:`J <options-props-J>` :ref:`K <options-props-K>` :ref:`L <options-props-L>`
:ref:`M <options-props-M>` :ref:`N <options-props-N>` :ref:`O <options-props-O>` :ref:`P <options-props-P>`
:ref:`Q <options-props-Q>` :ref:`R <options-props-R>` :ref:`S <options-props-S>` :ref:`T <options-props-T>`
:ref:`U <options-props-U>` :ref:`V <options-props-V>` :ref:`W <options-props-W>` :ref:`X <options-props-X>`
:ref:`Y <options-props-Y>` :ref:`Z <options-props-Z>`

**Category listing: jump to**
:ref:`Initial conditions <options-initial-conditions>`
:ref:`Stellar evolution and winds <options-stellar-evolution>`
:ref:`Mass transfer physics <options-mass-transfer>`
:ref:`Supernovae <options-supernovae>`
:ref:`Administrative <options-admin>`

COMPAS information
------------------

**--help [ -h ]** |br|
Prints COMPAS help.

**--version [ -v ]** |br|
Prints COMPAS version string.


Alphabetical listing
--------------------

.. _options-props-A:

**--add-options-to-sysparms** |br|
Add columns for program options to SSE System Parameters/BSE System Parameters file (mode dependent). |br|
Options: { ALWAYS, GRID, NEVER } |br|
Default = GRID

.. list-table::
   :widths: 11 80 
   :header-rows: 0
   :class: aligned-text

   * - ALWAYS
     - indicates that the program options should be added to the sysparms file
   * - GRID
     - indicates that the program options should be added to the sysparms file `only if`
   * -  
     - a GRID file is specified, or RANGEs or SETs are specified for options
   * - NEVER
     - indicates that the program options should `not` be added to the sysparms file

**--allow-rlof-at-birth** |br|
Allow binaries that have one or both stars in RLOF at birth to evolve as over-contact systems. |br|
Default = TRUE

**--allow-touching-at-birth** |br|
Allow binaries that are touching at birth to be included in the sampling. |br|
Default = FALSE

**--angular-momentum-conservation-during-circularisation** |br|
Conserve angular momentum when binary is circularised when entering a Mass Transfer episode. |br|
Default = FALSE

.. _options-props-B:

:ref:`Back to Top <options-props-top>`

**--black-hole-kicks** |br|
Black hole kicks relative to NS kicks. |br|
Options: { FULL, REDUCED, ZERO, FALLBACK } |br|
Default = FALLBACK

.. _options-props-C:

:ref:`Back to Top <options-props-top>`

**--case-bb-stability-prescription** |br|
Prescription for the stability of case BB/BC mass transfer. |br|
Options: { ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE } |br|
Case BB mass transfer is treated as always stable, always stable only for mass transfer onto neutron stars or black holes, with stability as determined for all other mass transfer, or always unstable, respectively |br|
Default = ALWAYS_STABLE

**--check-photon-tiring-limit** |br|
Check the photon tiring limit is not exceeded during mass loss. |br|
Default = FALSE

**--chemically-homogeneous-evolution** |br|
Chemically Homogeneous Evolution mode. See :cite:`Riley2021` for details of the implementation
of Chemically Homogeneous Evolution in COMPAS |br|
Options: { NONE, OPTIMISTIC, PESSIMISTIC } |br|
Default = PESSIMISTIC

**--circulariseBinaryDuringMassTransfer** |br|
Circularise binary when it enters a Mass Transfer episode. |br|
Default = TRUE

**--common-envelope-allow-immediate-RLOF-post-CE-survive** |br|
Allow binaries that experience Roche lobe overflow immediately at the end of the CE phase to survive. |br|
Default = FALSE

**--common-envelope-allow-main-sequence-survive** |br|
Allow main sequence accretors to survive common envelope evolution if other criteria point to survival. |br|
Default = TRUE

**--common-envelope-allow-radiative-envelope-survive** |br| 
Allow binaries with an evolved component with a radiative envelope to survive the common envelope phase. |br|
Deafult = FALSE

**--common-envelope-alpha** |br|
Common Envelope efficiency alpha. |br|
Default = 1.0

**--common-envelope-alpha-thermal** |br|
Thermal energy contribution to the total envelope binding energy. |br|
Defined such that :math:`λ = \alpha_{th} \times \lambda_{b} + (1.0 - \alpha_{th}) \times \lambda_{g}`. |br|
Default = 1.0

**--common-envelope-lambda** |br|
Common Envelope lambda. |br|
Only used when ``--common-envelope-lambda-prescription = LAMBDA_FIXED``. |br|
Default = 0.1

**--common-envelope-lambda-multiplier** |br|
Multiplicative constant to be applied to the common envelope lambda parameter for any prescription. |br|
Default = 1.0

**--common-envelope-lambda-nanjing-enhanced** |br|
Continuous extrapolation beyond maximum radius range in Nanjing lambda's as implemented in StarTrack. Only used when ``--common-envelope-lambda-prescription = LAMBDA_NANJING``. |br|
Default = FALSE

**--common-envelope-lambda-nanjing-interpolate-in-mass** |br|
Interpolate Nanjing lambda parameters across different mass models. Only used when ``--common-envelope-lambda-prescription = LAMBDA_NANJING``. |br|
Default = FALSE

**--common-envelope-lambda-nanjing-interpolate-in-metallicity** |br|
Interpolate Nanjing lambda parameters across population I and population II metallicity models. Only used when ``--common-envelope-lambda-prescription = LAMBDA_NANJING``. |br|
Default = FALSE

**--common-envelope-lambda-nanjing-use_rejuvenated-mass** |br|
Use rejuvenated or effective ZAMS mass instead of true birth mass when computing Nanjing lambda parameters. Only used when ``--common-envelope-lambda-prescription = LAMBDA_NANJING``. |br|
Default = FALSE

**--common-envelope-lambda-prescription** |br|
CE lambda (envelope binding energy) prescription. |br|
Options: { LAMBDA_FIXED, LAMBDA_LOVERIDGE, LAMBDA_NANJING, LAMBDA_KRUCKOW, LAMBDA_DEWI } |br|
``LAMBDA_FIXED`` is a constant; ``LAMBDA_LOVERIDGE`` is the prescription from Loveridge et al., 2011; ``LAMBDA_NANJING`` is from Xu & Li, 2010; ``LAMBDA_KRUCKOW`` is from Kruckow et al., 2016; and ``LAMBDA_DEWI`` is the fit from Appendix A in Claeys et al. 2014, based on Dewi & Tauris 2000 |br|
Default = LAMBDA_NANJING

**--common-envelope-mass-accretion-constant** |br|
Value of mass accreted by NS/BH during common envelope evolution if assuming all NS/BH accrete same amount of mass. |br|
Used when ``--common-envelope-mass-accretion-prescription = CONSTANT``, ignored otherwise. |br|
Default = 0.0

**--common-envelope-mass-accretion-max** |br|
Maximum amount of mass accreted by NS/BHs during common envelope evolution (:math:`M_\odot`). |br|
Default = 0.1

**--common-envelope-mass-accretion-min** |br|
Minimum amount of mass accreted by NS/BHs during common envelope evolution (:math:`M_\odot`). |br|
Default = 0.04

**--common-envelope-mass-accretion-prescription** |br|
Assumption about whether NS/BHs can accrete mass during common envelope evolution. |br|
``ZERO`` is no accretion; ``CONSTANT`` means a fixed amount of accretion determined by ``--common-envelope-mass-accretion-constant``; ``UNIFORM`` means a uniform random draw between ``--common-envelope-mass-accretion-min`` and ``--common-envelope-mass-accretion-max`` (Oslowski et al., 2011); and ``MACLEOD`` follows the prescription of MacLeod et al., 2015 |br|
Options: { ZERO, CONSTANT, UNIFORM, MACLEOD } |br|
Default = ZERO

**--common-envelope-recombination-energy-density** |br|
Recombination energy density (erg g−1). |br|
Default = :math:`1.5 \times 10^{13}`

**--common-envelope-slope-kruckow** |br|
Slope for the Kruckow lambda (see Kruckow et al. 2016 as implemented by Vigna-Gomez et al. 2018). |br|
Default = −0.833333

**--cool-wind-mass-loss-multiplier** |br|
Multiplicative constant for wind mass loss of cool stars, i.e. those with temperatures below the
VINK_MASS_LOSS_MINIMUM_TEMP (default 12500K). |br|
Only applicable when ``--mass-loss-prescription = VINK``. |br|
Default = 1.0

.. _options-props-D:

:ref:`Back to Top <options-props-top>`

**--debug-classes** |br|
Developer-defined debug classes to enable (vector). |br|
Default = `All debug classes enabled (e.g. no filtering)`

**--debug-level** |br|
Determines which print statements are displayed for debugging. |br|
Default = 0

**--debug-to-file** |br|
Write debug statements to file. |br|
Default = FALSE

**--detailed-output** |br|
Print BSE detailed information to file. |br|
Default = FALSE

.. _options-props-E:

:ref:`Back to Top <options-props-top>`

**--eccentricity [ -e ]** |br|
Initial eccentricity for a binary star when evolving in BSE mode.
Default = 0.0 |br|

**--eccentricity-distribution** |br|
Initial eccentricity distribution. |br|
Options: { ZERO, FLAT, GELLER+2013, THERMAL, DUQUENNOYMAYOR1991, SANA2012 } |br|
``ZERO`` always circular, ``FLAT`` is uniform on [``--eccentricity-min``,``--eccentricity-max``], ``THERMAL`` is p(e) proportional to e, and the other options refer to the distributions of Geller et al. 2013, Duqennoy & Mayor 1991, and Sana et al. 2012. |br|
Default = ZERO

**--eccentricity-max** |br|
Maximum eccentricity to generate. |br|
Default = 1.0

**--eccentricity-min** |br|
Minimum eccentricity to generate. |br|
Default = 0.0

**--eddington-accretion-factor** |br|
Multiplication factor for Eddington accretion for NS & BH (i.e. > 1 is super-eddington and 0 is no accretion). |br|
Default = 1.0

**--enable-warnings** |br|
Display warning messages to stdout. |br|
Default = FALSE

**--envelope-state-prescription** |br|
Prescription for determining whether the envelope of the star is convective or radiative. |br|
Options: { LEGACY, HURLEY, FIXED_TEMPERATURE } |br|
``LEGACY`` refers to the model used in Stevenson et al., 2017; ``HURLEY`` refers to the model of Hurley, Pols, Tout, 2002; and ``FIXED_TEMPERATURE`` assumes that a deep convective envelope developes only when the temperature drops below ``CONVECTIVE_BOUNDARY_TEMPERATURE`` (Klencki et al., 2020) |br|
Default = LEGACY


**--errors-to-file** |br|
Write error messages to file. |br|
Default = FALSE

**--evolve-pulsars** |br|
Evolve pulsar properties of Neutron Stars. |br|
Default = FALSE

**--evolve-unbound-systems** |br|
Continue evolving stars even if the binary is disrupted. |br|
Default = FALSE

.. _options-props-F:

:ref:`Back to Top <options-props-top>`

**--fix-dimensionless-kick-magnitude** |br|
Fix dimensionless kick magnitude to this value. |br|
Default = n/a (not used if option not present)

**--fryer-supernova-engine** |br|
Supernova engine type if using the fallback prescription from :cite:`Fryer2012`. |br|
Options: { DELAYED, RAPID }
Default = DELAYED

.. _options-props-G:

:ref:`Back to Top <options-props-top>`

**--grid** |br|
Grid filename. |br|
Default = ’’ (None)

**--grid-lines-to-process** |br|
The number of grid file lines to be processed. |br|
Default = Process to EOF

**--grid-start-line** |br|
The first line of the grid file to be processed. |br|
Default = 0

.. _options-props-H:

:ref:`Back to Top <options-props-top>`

**--hdf5-buffer-size** |br|
The ``HDF5`` IO buffer size for writing to ``HDF5`` logfiles (number of ``HDF5`` chunks). |br|
Default = 1

**--hdf5-chunk-size** |br|
The ``HDF5`` dataset chunk size to be used when creating ``HDF5`` logfiles (number of logfile entries). |br|
Default = 100000

**--help [ -h ]** |br|
Prints COMPAS help (-h is short form, --help includes more information).

.. _options-props-I:

:ref:`Back to Top <options-props-top>`

**--initial-mass** |br|
Initial mass for a single star when evolving in SSE mode (:math:`M_\odot`). |br|
Default = Sampled from IMF

**--initial-mass-1** |br|
Initial mass for the primary star when evolving in BSE mode (:math:`M_\odot`). |br|
Default = Sampled from IMF

**--initial-mass-2** |br|
Initial mass for the secondary star when evolving in BSE mode (:math:`M_\odot`). |br|
Default = Sampled from IMF

**--initial-mass-function [ -i ]** |br|
Initial mass function. |br|
Options: { SALPETER, POWERLAW, UNIFORM, KROUPA } |br|
``SALPETER`` and ``KROUPA`` use the IMFs of Salpeter 1955 and Kroupa 2001, ``POWERLAW`` samples from a single power law with slope ``--initial-mass-power``, and ``UNIFORM`` samples uniformly between ``--initial-mass-min`` and ``--initial-mass-min`` |br|
Default = KROUPA

**--initial-mass-max** |br|
Maximum mass to generate using given IMF (:math:`M_\odot`). |br|
Default = 150.0

**--initial-mass-min** |br|
Minimum mass to generate using given IMF (:math:`M_\odot`). |br|
Default = 5.0

**--initial-mass-power** |br|
Single power law power to generate primary mass using ``POWERLAW`` IMF. |br|
Default = 0.0

.. _options-props-J:

.. _options-props-K:

:ref:`Back to Top <options-props-top>`

**--kick-direction** |br|
Natal kick direction distribution. |br|
Options: { ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES } |br|
Kick angles are defined relative to the spin axis.  ``INPLANE`` and ``PERPENDICULAR`` are strictly in the equatorial plane or in polar directions, while ``WEDGE`` and ``POLES`` are preferentially but exactly in the equatorial plane or in polar directions with 1 degree scales, respectively;  ``POWERLAW`` quantifies the preference for polar vs planar kicks with the ``--kick-direction-power`` parameter. |br|
Default = ISOTROPIC

**--kick-direction-power** |br|
Power for power law kick direction distribution, where 0.0 = isotropic, +ve = polar, -ve = in plane. |br|
Default = 0.0 (isotropic)

**--kick-magnitude** |br|
Value to be used as the (drawn) kick magnitude for a single star when evolving in SSE mode, should the star
undergo a supernova event (:math:`km s^{−1}`). |br|
If a value for option ``--kick-magnitude-random`` is specified, it will be used in preference to ``--kick-magnitude``. |br|
Default = 0.0

**--kick-magnitude-1** |br|
Value to be used as the (drawn) kick magnitude for the primary star of a binary system when evolving in
BSE mode, should the star undergo a supernova event (:math:`km s^{−1}`). |br|
If a value for option ``--kick-magnitude-random-1`` is specified, it will be used in preference to ``--kick-magnitude-1``. |br|
Default = 0.0

**--kick-magnitude-2** |br|
Value to be used as the (drawn) kick magnitude for the secondary star of a binary system when evolving in
BSE mode, should the star undergo a supernova event (:math:`km s^{−1}`). |br|
If a value for option ``--kick-magnitude-random-2`` is specified, it will be used in preference to ``--kick-magnitude-2``. |br|
Default = 0.0

**--kick-magnitude-distribution** |br|
Natal kick magnitude distribution. |br|
Options: { ZERO, FIXED, FLAT, MAXWELLIAN, BRAYELDRIDGE, MULLER2016, MULLER2016MAXWELLIAN, MULLERMANDEL } |br|
``ZERO`` assigns kick magnitudes of 0, ``FIXED`` always sets the magnitude to a fixed value based on supernova type, ``FLAT`` and ``MAXWELLIAN`` draw kicks from uniform or Maxwellian (e.g., Hobbs et al., 2005) distributions, respectively, ``BRAYELDRIDGE`` and ``MULLERMANDEL`` use momenum-preserving kicks from Bray & Eldrigde 2018 and Mandel & Mueller 2020, respectively, and ``MULLER2016`` and ``MULLER2016MAXWELLIAN`` use kicks from Mueller 2016 as implemented in Vigna-Gomez et al., 2018 (reduced by a factor of sqrt(3) in the latter case). |br|
Default = MAXWELLIAN

**--kick-magnitude-max** |br|
Maximum drawn kick magnitude (:math:`km s^{−1}`). |br|
Must be > 0 if using ``--kick-magnitude-distribution = FLAT``. |br|
Default = −1.0

**--kick-magnitude-random** |br|
CDF value to be used to draw the kick magnitude for a single star when evolving in SSE mode, should the star undergo a supernova event and should the chosen distribution sample from a cumulative distribution function. |br|
Must be a floating-point number in the range :math:`[0.0, 1.0)`. |br|
The specified value for this option will be used in preference to any specified value for ``--kick-magnitude``. |br|
Default = Random number drawn uniformly from :math:`[0.0, 1.0)`

**--kick-magnitude-random-1** |br|
CDF value to be used to draw the kick magnitude for the primary star of a binary system when evolving in BSE mode, should the star undergo a supernova event and should the chosen distribution sample from a cumulative distribution function. |br|
Must be a floating-point number in the range :math:`[0.0, 1.0)`. |br|
The specified value for this option will be used in preference to any specified value for ``--kick-magnitude-1``. |br|
Default = Random number drawn uniformly from :math:`[0.0, 1.0)`

**--kick-magnitude-random-2** |br|
CDF value to be used to draw the kick magnitude for the secondary star of a binary system when evolving in BSE mode, should the star undergo a supernova event and should the chosen distribution sample from a cumulative distribution function. |br|
Must be a floating-point number in the range :math:`[0.0, 1.0)`. |br|
The specified value for this option will be used in preference to any specified value for ``--kick-magnitude-2``. |br|
Default = Random number drawn uniformly from :math:`[0.0, 1.0)`

**--kick-magnitude-sigma-CCSN-BH** |br|
Sigma for chosen kick magnitude distribution for black holes (:math:`km s^{−1}`); ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 265.0

**--kick-magnitude-sigma-CCSN-NS** |br|
Sigma for chosen kick magnitude distribution for neutron stars (:math:`km s^{−1}`); ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 265.0

**--kick-magnitude-sigma-ECSN** |br|
Sigma for chosen kick magnitude distribution for ECSN (:math:`km s^{−1}`); ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 30.0

**--kick-magnitude-sigma-USSN** |br|
Sigma for chosen kick magnitude distribution for USSN (:math:`km s^{−1}`); ignored if not needed for the chosen kick magnitude distribution. |br|
Default = 30.0

**--kick-mean-anomaly-1** |br|
The mean anomaly at the instant of the supernova for the primary star of a binary system when evolving in
BSE mode, should it undergo a supernova event. |br|
Must be a floating-point number in the range :math:`[0.0, 2\pi)`. |br|
Default = Random number drawn uniformly from :math:`[0.0, 2\pi)`

**--kick-mean-anomaly-2** |br|
The mean anomaly at the instant of the supernova for the secondary star of a binary system when evolving
in BSE mode, should it undergo a supernova event. |br|
Must be a floating-point number in the range :math:`[0.0, 2\pi)`. |br|
Default = Random number drawn uniformly from :math:`[0.0, 2\pi)`

**--kick-phi-1** |br|
The angle between ’x’ and ’y’, both in the orbital plane of the supernova vector, for the primary star of a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

**--kick-phi-2** |br|
The angle between ’x’ and ’y’, both in the orbital plane of the supernova vector, for the secondary
star of a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

**--kick-scaling-factor** |br|
Arbitrary factor used to scale kicks. |br|
Default = 1.0

**--kick-theta-1** |br|
The angle between the orbital plane and the ’z’ axis of the supernova vector for the primary star of a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

**--kick-theta-2** |br|
The angle between the orbital plane and the ’z’ axis of the supernova vector for the secondary star of a binary system when evolving in BSE mode, should it undergo a supernova event (radians). |br|
Default = Drawn according to specified ``--kick-direction`` distribution

.. _options-props-L:

:ref:`Back to Top <options-props-top>`

**--log-classes** |br|
Logging classes to be enabled (vector). |br|
Default = `All debug classes enabled (e.g. no filtering)`

**--logfile-common-envelopes** |br|
Filename for Common Envelopes logfile (BSE mode). |br|
Default = ’BSE_Common_Envelopes’

**--logfile-definitions** |br|
Filename for logfile record definitions file. |br|
Default = ’’ (None)

**--logfile-detailed-output** |br|
Filename for the Detailed_Output logfile. |br|
Default = ’SSE_Detailed_Output’ for SSE mode; ’BSE_Detailed_Output’ for BSE mode |br|

**--logfile-double-compact-objects** |br|
Filename for the Double Compact Objects logfile (BSE mode). |br|
Default = ’BSE_Double_Compact_Objects’

**--logfile-name-prefix** |br|
Prefix for logfile names. |br|
Default = ’’ (None)

**--logfile-pulsar-evolution** |br|
Filename for the Pulsar Evolution logfile (BSE mode). |br|
Default = ’BSE_Pulsar_Evolution’

**--logfile-rlof-parameters** |br|
Filename for the RLOF Printing logfile (BSE mode). |br|
Default = ’BSE_RLOF’

**--logfile-supernovae** |br|
Filename for the Supernovae logfile. |br|
Default = ’SSE_Supernovae’ for SSE mode; ’BSE_Supernovae’ for BSE mode |br|

**--logfile-switch-log** |br|
Filename for the Switch Log logfile. |br|
Default = ’SSE_Switch_Log’ for SSE mode; ’BSE_Switch_Log’ for BSE mode |br|

**--logfile-system-parameters** |br|
Filename for the System Parameters logfile (BSE mode). |br|
Default = ’SSE_System_Parameters’ for SSE mode; ’BSE_System_Parameters’ for BSE mode |br|

**--logfile-type** |br|
The type of logfile to be produced by COMPAS. |br|
Default = ’HDF5’

**--log-level** |br|
Determines which print statements are included in the logfile. |br|
Default = 0

**--luminous-blue-variable-multiplier** |br|
Multiplicative constant for LBV mass loss. (Use 10 for Mennekens & Vanbeveren (2014)). |br|
Note that wind mass loss will also be multiplied by the ``--overall-wind-mass-loss-multiplier``. |br|
Default = 1.5

**--luminous-blue-variable-prescription** |br|
Luminous blue variable mass loss prescription. |br|
Options: { NONE, HURLEY, HURLEY_ADD, BELCZYNSKI } |br|
No LBV winds for ``NONE``,  Hurley, Pols, Tout (2000) LBV winds only for ``HURLEY`` LBV stars (or in addition to other winds for ``HURLEY_ADD``, Belzcynski et al. 2010 winds for ``BELCZYNSKI`` |br|
Default = HURLEY_ADD

.. _options-props-M:

:ref:`Back to Top <options-props-top>`

**--mass-loss-prescription** |br|
Mass loss prescription. |br|
Options: { NONE, HURLEY, VINK } |br|
Default = VINK

**--mass-ratio [ -q ]** |br|
Mass ratio :math:`\frac{m2}{m1}` used to determine secondary mass if not specified via ``--initial-mass-2``. |br|
Default = Value is sampled if option not specified.

**--mass-ratio-distribution** |br|
Initial mass ratio distribution for :math:`q = \frac{m2}{m1}`. |br|
Options: { FLAT, DUQUENNOYMAYOR1991, SANA2012 } |br|
``FLAT`` is uniform in the mass ratio between ``--mass-ratio-min`` and ``--mass-ratio-max``, the other prescriptions follow Duquennoy & Mayor 1991 and Sana et al. 2012 |br|
Default = FLAT

**--mass-ratio-max** |br|
Maximum mass ratio :math:`\frac{m2}{m1}` to generate. |br|
Default = 1.0

**--mass-ratio-min** |br|
Minimum mass ratio :math:`\frac{m2}{m1}` to generate. |br|
Default = 0.01

**--mass-transfer** |br|
Enable mass transfer. |br|
Default = TRUE

**--mass-transfer-accretion-efficiency-prescription** |br|
Mass transfer accretion efficiency prescription. |br|
Options: { THERMAL, FIXED, CENTRIFUGAL } |br|
Default = THERMAL

**--mass-transfer-angular-momentum-loss-prescription** |br|
Mass Transfer Angular Momentum Loss prescription. |br|
Options: { JEANS, ISOTROPIC, CIRCUMBINARY, ARBITRARY } |br|
Default = ISOTROPIC

**--mass-transfer-fa** |br|
Mass Transfer fraction accreted. |br|
Used when ``--mass-transfer-accretion-efficiency-prescription = FIXED_FRACTION``. |br|
Default = 0.5

**--mass-transfer-jloss** |br|
Specific angular momentum with which the non-accreted system leaves the system. |br|
Used when ``--mass-transfer-angular-momentum-loss-prescription = ARBITRARY``, ignored otherwise. |br|
Default = 1.0

**--mass-transfer-rejuvenation-prescription** |br|
Mass Transfer Rejuvenation prescription. |br|
Options: { NONE, STARTRACK } |br|
``NONE`` uses the Hurley, Pols, Tout (2000) model, ``STARTRACK`` uses the model from Belczynski et al. 2008 |br|
Default = STARTRACK

**--mass-transfer-thermal-limit-accretor** |br|
Mass Transfer Thermal Accretion limit multiplier. |br|
Options: { CFACTOR, ROCHELOBE } |br|
Default = CFACTOR

**--mass-transfer-thermal-limit-C** |br|
Mass Transfer Thermal rate factor for the accretor. |br|
Default = 10.0

**--maximum-evolution-time** |br|
Maximum time to evolve binaries (Myr). Evolution of the binary will stop if this number is reached. |br|
Default = 13700.0

**--maximum-mass-donor-nandez-ivanova** |br|
Maximum donor mass allowed for the revised common envelope formalism of Nandez & Ivanova (:math:`M_\odot`). |br|
Default = 2.0

**--maximum-neutron-star-mass** |br|
Maximum mass of a neutron star (:math:`M_\odot`). |br|
Default = 2.5

**--maximum-number-timestep-iterations** |br|
Maximum number of timesteps to evolve binary. Evolution of the binary will stop if this number is reached. |br|
Default = 99999

**--mcbur1** |br|
Minimum core mass at base of AGB to avoid fully degnerate CO core formation (:math:`M_\odot`). |br|
e.g. 1.6 in :cite:`Hurley2000` presciption; 1.83 in :cite:`Fryer2012` and :doc:`Belczynski et al. (2008) <../../references>` models. |br|
Default = 1.6

**--metallicity [ -z ]** |br|
Metallicity. |br|
The value specified for metallicity is applied to both stars for BSE mode. |br|
Default = 0.0142

**--metallicity-distribution** |br|
Metallicity distribution. |br|
Options: { ZSOLAR, LOGUNIFORM } |br|
``ZSOLAR`` uses ``ZSOL_ASPLUND`` for all initial metallicities, ``LOGUNIFORM`` draws the metallicity uniformly in the log between ``--metallicity-min`` and ``--metallicity-max`` |br|
Default = ZSOLAR

**--metallicity-max** |br|
Maximum metallicity to generate. |br|
Default = 0.03

**--metallicity-min** |br|
Minimum metallicity to generate. |br|
Default = 0.0001

**--minimum-secondary-mass** |br|
Minimum mass of secondary to generate (:math:`M_\odot`). |br|
Default = 0.1 if ``--initial-mass-2`` specified; value of ``--initial-mass-min`` if ``--initial-mass-2`` not specified.

**--mode** |br|
The mode of evolution. |br|
Options: { SSE, BSE } |br|
Default = BSE

**--muller-mandel-kick-multiplier-BH** |br|
Scaling prefactor for BH kicks when using the `MULLERMANDEL` kick magnitude distribution |br|
Default = 200.0

**--muller-mandel-kick-multiplier-NS** |br|
Scaling prefactor for NS kicks when using the `MULLERMANDEL` kick magnitude distribution |br|
Default = 400.0

.. _options-props-N:

:ref:`Back to Top <options-props-top>`

**--neutrino-mass-loss-BH-formation** |br|
Assumption about neutrino mass loss during BH formation. |br|
Options: { FIXED_FRACTION, FIXED_MASS } |br|
Default = FIXED_MASS

**--neutrino-mass-loss-BH-formation-value** |br|
Amount of mass lost in neutrinos during BH formation (either as fraction or in solar masses, depending on the value of ``--neutrino-mass-loss-bh-formation``). |br|
Default = 0.1

**--neutron-star-equation-of-state** |br|
Neutron star equation of state. |br|
Options: { SSE, ARP3 } |br|
Default = SSE

**--notes** |br|
Annotation strings (vector). |br|
Default = ""

**--notes-hdrs** |br|
Annotations header strings (vector). |br|
Default = `No annotations`

**--number-of-systems [ -n ]** |br|
The number of systems to simulate. |br|
Single stars for SSE mode; binary stars for BSE mode. |br|
This option is ignored if either of the following is true: |br|

    - the user specified a grid file |br|
    - the user specified a range or set for any options - this implies a grid |br|

In both cases the number of objects evolved will be the number specified by the grid. |br|
Default = 10

.. _options-props-O:

:ref:`Back to Top <options-props-top>`

**--orbital-period** |br|
Initial orbital period for a binary star when evolving in BSE mode (days). |br|
Used only if the semi-major axis is not specified via ``--semi-major-axis``. |br|
Default = Value is sampled if option not specified.

**--orbital-period-distribution** |br|
Initial orbital period distribution. |br|
Options: { FLATINLOG } |br|
Default = FLATINLOG

**--orbital-period-max** |br|
Maximum period to generate (days). |br|
Default = 1000.0

**--orbital-period-min** |br|
Minimum period to generate (days). |br|
Default = 1.1

**--output-container [ -c ]** |br|
Container (directory) name for output files. |br|
Default = ’COMPAS_Output’

**--output-path [ -o ]** |br|
Path to which output is saved (i.e. directory in which the output container is created). |br|
Default = Current working directory (CWD)

**--overall-wind-mass-loss-multiplier** |br|
Multiplicative constant for overall wind mass loss. |br|
Note that this multiplication factor is applied after the ``luminous-blue-variable-multiplier``,
the ``wolf-rayet-multiplier`` and the ``cool-wind-mass-loss-multiplier``. |br|
Default = 1.0

.. _options-props-P:

:ref:`Back to Top <options-props-top>`

**--pair-instability-supernovae** |br|
Enable pair instability supernovae (PISN). |br|
Default = TRUE

**--PISN-lower-limit** |br|
Minimum core mass for PISN (:math:`M_\odot`). |br|
Default = 60.0

**--PISN-upper-limit** |br|
Maximum core mass for PISN (:math:`M_\odot`). |br|
Default = 135.0

**--population-data-printing** |br|
Print details of population. |br|
Default = FALSE

**--PPI-lower-limit** |br|
Minimum core mass for PPI (:math:`M_\odot`). |br|
Default = 35.0

**--PPI-upper-limit** |br|
Maximum core mass for PPI (:math:`M_\odot`). |br|
Default = 60.0

**--print-bool-as-string** |br|
Print boolean properties as ’TRUE’ or ’FALSE’. |br|
Default = FALSE

**--pulsar-birth-magnetic-field-distribution** |br|
Pulsar birth magnetic field distribution. |br|
Options: { ZERO, FIXED, FLATINLOG, UNIFORM, LOGNORMAL } |br|
Default = ZERO

**--pulsar-birth-magnetic-field-distribution-max** |br|
Maximum (:math:`log_{10}`) pulsar birth magnetic field. |br|
Default = 13.0

**--pulsar-birth-magnetic-field-distribution-min** |br|
Minimum (:math:`log_{10}`) pulsar birth magnetic field. |br|
Default = 11.0

**--pulsar-birth-spin-period-distribution** |br|
Pulsar birth spin period distribution. |br|
Options: { ZERO, FIXED, UNIFORM, NORMAL } |br|
Default = ZERO

**--pulsar-birth-spin-period-distribution-max** |br|
Maximum pulsar birth spin period (ms). |br|
Default = 100.0

**--pulsar-birth-spin-period-distribution-min** |br|
Minimum pulsar birth spin period (ms). |br|
Default = 10.0

**--pulsar-magnetic-field-decay-massscale** |br|
Mass scale on which magnetic field decays during accretion (:math:`M_\odot`). |br|
Default = 0.025

**--pulsar-magnetic-field-decay-timescale** |br|
Timescale on which magnetic field decays (Myr). |br|
Default = 1000.0

**--pulsar-minimum-magnetic-field** |br|
:math:`log_{10}` of the minimum pulsar magnetic field (Gauss). |br|
Default = 8.0

**--pulsational-pair-instability** |br|
Enable mass loss due to pulsational-pair-instability (PPI). |br|
Default = TRUE

**--pulsational-pair-instability-prescription** |br|
Pulsational pair instability prescription. |br|
Options: { COMPAS, STARTRACK, MARCHANT, FARMER } |br|
``COMPAS``, ``STARTRACK`` and ``MARCHANT`` follow Woosley 2017, Belczynski et al. 2016, and Marchant et al. 2018, all as implemented in Stevenson et al. 2019, ``FARMER`` follows Farmer et al. 2019 |br|
Default = MARCHANT

.. _options-props-Q:

:ref:`Back to Top <options-props-top>`

**--quiet** |br|
Suppress printing to stdout. |br|
Default = FALSE

.. _options-props-R:

:ref:`Back to Top <options-props-top>`

**--random-seed** |br|
Value to use as the seed for the random number generator. |br|
Default = 0

**--remnant-mass-prescription** |br|
Remnant mass prescription. |br|
Options: { HURLEY2000, BELCZYNSKI2002, FRYER2012, MULLER2016, MULLERMANDEL, SCHNEIDER2020, SCHNEIDER2020ALT } |br|
Remnant mass recipes from Hurley, Pols, Tout (2000) for ``HURLEY2000``, Belczynski et al. 2002, Fryer et al. 2012, Mueller 2016, Mandel & Mueller 2020, and Schneider et al. 2020 (with the alternative prescription for effectively single stars from the same paper in the ``SCHNEIDER2020ALT`` case) |br|
Default = FRYER2012

**--revised-energy-formalism-nandez-ivanova** |br|
Enable revised energy formalism of Nandez & Ivanova. |br|
Default = FALSE

**--rlof-printing** |br|
Print RLOF events to logfile. |br|
Default = TRUE

**--rotational-frequency** |br|
Initial rotational frequency of the star for SSE (Hz). |br|
Default = 0.0 (``--rotational-velocity-distribution`` used if ``--rotational-frequency`` not specified)

**--rotational-frequency-1** |br|
Initial rotational frequency of the primary star for BSE (Hz). |br|
Default = 0.0 (``--rotational-velocity-distribution`` used if ``--rotational-frequency-1`` not specified)

**--rotational-frequency-2** |br|
Initial rotational frequency of the secondary star for BSE (Hz). |br|
Default = 0.0 (``--rotational-velocity-distribution`` used if ``--rotational-frequency-2`` not specified)

**--rotational-velocity-distribution** |br|
Initial rotational velocity distribution. |br|
Options: { ZERO, HURLEY, VLTFLAMES } |br|
``ZERO`` sets all initial rotational velocities to 0, while ``HURLEY`` and ``VLTFLAMES`` sample them from the Hurley, Pols, Tout (2000) and Ramirez-Agudelo et al. (2013,2015), respectively |br|
Default = ZERO

.. _options-props-S:

:ref:`Back to Top <options-props-top>`

**--semi-major-axis** |br|
Initial semi-major axis for a binary star when evolving in BSE mode (AU). |br|
Default = 0.1

**--semi-major-axis-distribution [ -a ]** |br|
Initial semi-major axis distribution. |br|
Options: { FLATINLOG, DUQUENNOYMAYOR1991, SANA2012 } |br|
Flat-in-log (Opik 1924), Duquennoy & Mayor (1991) or Sana et al. (2012) distributions |BR|
Default = FLATINLOG

**--semi-major-axis-max** |br|
Maximum semi-major axis to generate (AU). |br|
Default = 1000.0

**--semi-major-axis-min** |br|
Minimum semi-major axis to generate (AU). |br|
Default = 0.01

**--stellar-zeta-prescription** |br|
Prescription for stellar zeta. |br|
Options: { SOBERMAN, HURLEY, ARBITRARY } |br|
Use Soberman, Phinney, and van den Heuvel (1997) or Hurley, Pols, Tout (2002) or the fixed value specified via ``--zeta-adiabatic-arbitrary`` for the stellar radial response to mass loss for convective-envelope giant-like stars |br|
Default = SOBERMAN

**--store-input-files** |br|
Enables copying of any specified grid file and/or logfile-definitios file to the COMPAS output container. |br|
Default = TRUE

**--switch-log** |br|
Enables printing of the Switch Log logfile. |br|
Default = FALSE

.. _options-props-T:

:ref:`Back to Top <options-props-top>`

**--timestep-multiplier** |br|
Multiplicative factor for timestep duration. |br|
Default = 1.0

.. _options-props-U:

:ref:`Back to Top <options-props-top>`

**--use-mass-loss** |br|
Enable mass loss. |br|
Default = TRUE

.. _options-props-V:

:ref:`Back to Top <options-props-top>`

**--version [ -v ]** |br|
Prints COMPAS version string.

.. _options-props-W:

:ref:`Back to Top <options-props-top>`

**--wolf-rayet-multiplier** |br|
Multiplicative constant for Wolf Rayet winds. Note that wind mass loss will also be multiplied by the
``overall-wind-mass-loss-multiplier``. |br|
Default = 1.0

.. _options-props-X:
.. _options-props-Y:
.. _options-props-Z:

:ref:`Back to Top <options-props-top>`

**--zeta-adiabatic-arbitrary** |br|
Value of logarithmic derivative of radius with respect to mass, :math:`\zeta` adiabatic. |br|
Default = :math:`1.0 \times 10^4`

**--zeta-main-sequence** |br|
Value of logarithmic derivative of radius with respect to mass, :math:`\zeta` on the main sequence. |br|
Default = 2.0

**--zeta-radiative-giant-star** |br|
Value of logarithmic derivative of radius with respect to mass, :math:`\zeta` for radiative-envelope giant-like stars
(including Hertzsprung Gap (HG) stars). |br|
Default = 6.5


Category listing
--------------------

Go to :ref:`the top of this page <options-props-top>` for the full alphabetical list of options with explanations and default values

.. _options-initial-conditions:

**Initial conditions**

--initial-mass-function, --initial-mass, --initial-mass-1, --initial-mass-2, --initial-mass-min, --initial-mass-max, --initial-mass-power

--mass-ratio-distribution, --mass-ratio, --mass-ratio-min, --mass-ratio-max, --minimum-secondary-mass

--eccentricity-distribution, --eccentricity, --eccentricity-min, --eccentricity-max

--metallicity-distribution, --metallicity, --metallicity-min, --metallicity-max

--orbital-period-distribution, --orbital-period, --orbital-period-min, --orbital-period-max, --semi-major-axis-distribution, --semi-major-axis, --semi-major-axis-min, --semi-major-axis-max, --allow-rlof-at-birth, --allow-touching-at-birth

--rotational-velocity-distribution, --rotational-frequency, --rotational-frequency-1, --rotational-frequency-2

:ref:`Back to Top <options-props-top>`

.. _options-stellar-evolution:

**Stellar evolution and winds**

--use-mass-loss, --check-photon-tiring-limit, --cool-wind-mass-loss-multiplier, --luminous-blue-variable-prescription, --luminous-blue-variable-multiplier, --mass-loss-prescription, --overall-wind-mass-loss-multiplier, --wolf-rayet-multiplier

--chemically-homogeneous-evolution

:ref:`Back to Top <options-props-top>`

.. _options-mass-transfer:

**Mass transfer physics**

--mass-transfer, --mass-transfer-accretion-efficiency-prescription, --mass-transfer-angular-momentum-loss-prescription, --mass-transfer-fa, --mass-transfer-jloss, --mass-transfer-rejuvenation-prescription, --mass-transfer-thermal-limit-accretor, --mass-transfer-thermal-limit-C, --stellar-zeta-prescription, --zeta-adiabatic-arbitrary, --zeta-main-sequence, --zeta-radiative-giant-star, --case-bb-stability-prescription, --eddington-accretion-factor

--circulariseBinaryDuringMassTransfer, --angular-momentum-conservation-during-circularisation

--envelope-state-prescription, --common-envelope-alpha, --common-envelope-alpha-thermal, --common-envelope-lambda-prescription, --common-envelope-lambda, --common-envelope-slope-kruckow, --common-envelope-lambda-multiplier, --common-envelope-lambda-nanjing-enhanced, --common-envelope-lambda-nanjing-interpolate-in-mass, --common-envelope-lambda-nanjing-interpolate-in-metallicity, --common-envelope-lambda-nanjing-use_rejuvenated-mass, --common-envelope-allow-main-sequence-survive, --common-envelope-allow-radiative-envelope-survive*, --common-envelope-allow-immediate-RLOF-post-CE-survive, --common-envelope-mass-accretion-prescription, --common-envelope-mass-accretion-constant, --common-envelope-mass-accretion-min, --common-envelope-mass-accretion-max, --common-envelope-recombination-energy-density, --maximum-mass-donor-nandez-ivanova, --revised-energy-formalism-nandez-ivanova

:ref:`Back to Top <options-props-top>`

.. _options-supernovae:

**Supernovae**

--remnant-mass-prescription, --fryer-supernova-engine, --maximum-neutron-star-mass, --mcbur1, --neutrino-mass-loss-BH-formation, --neutrino-mass-loss-BH-formation-value, --neutron-star-equation-of-state

--pair-instability-supernovae, --PISN-lower-limit, --PISN-upper-limit, --PPI-lower-limit, --PPI-upper-limit, --pulsational-pair-instability, --pulsational-pair-instability-prescription

--pulsar-birth-magnetic-field-distribution, --pulsar-birth-magnetic-field-distribution-min, --pulsar-birth-magnetic-field-distribution-max, --pulsar-birth-spin-period-distribution, --pulsar-birth-spin-period-distribution-min, --pulsar-birth-spin-period-distribution-max, --pulsar-magnetic-field-decay-massscale, --pulsar-magnetic-field-decay-timescale, --pulsar-minimum-magnetic-field

--kick-magnitude-distribution, --kick-magnitude-sigma-CCSN-BH, --kick-magnitude-sigma-CCSN-NS, --kick-magnitude-sigma-ECSN, --kick-magnitude-sigma-USSN, --black-hole-kicks, --fix-dimensionless-kick-magnitude, --kick-magnitude, --kick-magnitude-1, --kick-magnitude-2, --kick-magnitude-min, --kick-magnitude-max, --kick-magnitude-random, --kick-magnitude-random-1, --kick-magnitude-random-2, --kick-scaling-factor, -muller-mandel-kick-multiplier-BH, --muller-mandel-kick-multiplier-NS

--kick-direction, --kick-direction-power, --kick-mean-anomaly-1, --kick-mean-anomaly-2, --kick-phi-1, --kick-phi-2, --kick-theta-1, --kick-theta-2

:ref:`Back to Top <options-props-top>`

.. _options-admin:

**Administrative**

--mode, --number-of-systems, --evolve-pulsars, --evolve-unbound-systems, --maximum-evolution-time, --maximum-number-timestep-iterations, --random-seed, --timestep-multiplier

--grid, --grid-start-line, --grid-lines-to-process

--add-options-to-sysparms, --debug-classes, --debug-level, --debug-to-file, --detailed-output, --enable-warnings, --errors-to-file, --help, --notes, --notes-hdrs, --population-data-printing, --print-bool-as-string, --quiet, --version

--log-classes, --logfile-definitions, --logfile-name-prefix, --logfile-type, --log-level, --logfile-common-envelopes, --logfile-detailed-output, --logfile-double-compact-objects, --logfile-pulsar-evolution, --logfile-rlof-parameters, --logfile-supernovae, --logfile-switch-log, --logfile-system-parameters, --output-container, --output-path, --rlof-printing, --store-input-files, --switch-log, --hdf5-buffer-size, --hdf5-chunk-size

:ref:`Back to Top <options-props-top>`


Vector program options
======================

Most of the program options available in COMPAS allow users to specify a `single` value for the option (e.g. ``--initial-mass-1 7.5``,
or ``--metallicity 0.0142``, etc.), while others allow users to specify multiple values (e.g. ``--log-classes class1 class2 class3``,
``--notes "this is note 1" note2 a_third_note``, etc.). The latter are `vector` options, because they allow users to specify a vector
of values.

Currently, the only vector program options are:

- --debug-classes
- --log-classes
- --notes-hdrs
- --notes


The notation for vector program options provides for the specification of one or more values. e.g.::

    --log-classes class1 class2 class3 ... classN

Option values are separated by a space character. Unless otherwise specified in the documentation for the program option, there is no
limit to the number of values specified.

When using this notation, all options required must be provided: there is no mechanism to allow a default value using the fallback method
for program options (to the command-line value, then to the COMPAS default) - leaving a value blank would be ambiguous (as to which e.g.
`log-class` had been left blank), and specifying an empty string ("") for a value would be ambiguous (as to whether the user wanted the
option value to default, or just be an empty string).

Option values (in general, but also specifically for vector options) may not begin with the dash character ('-'), because the shell parser
will parse them as option names before passing them through to COMPAS.

COMPAS imposes no limit to the length (number of characters) of an individual option values that are specified as strings, but there may 
be practical limits imposed by the underlying system.


Shorthand notation
------------------

Because the notation described above could become awkward, and to allow for default values for vector options, a shorthand notation for
vector options is provided. Usage using the shorthand notation is::

    --debug-classes [class1,headerclass2str2,class3,...,classN]

    --notes [annotation1,annotation2,annotation3,...,annotationN]

Because the parameters are bounded by the brackets, and delimited by commas (and so are now positional), users can omit specific values::

    --notes [,,annotation3,,annotation5]

In the example above, annotations 1, 2, 4, and those beyond annotation 5 have been omitted. Annotations 1, 2 & 4 will default - if they are
specified in this manner on a grid line they will default to the correspodning annotation specified on the command line; if they are specified
in this manner on the command line they will default to the COMPAS default annotation.

Spaces in option values (in general, but also specifically for vector options) strings need to be enclosed in quotes, or the shell parser will
parse them as separate arguments before passing them through to COMPAS.  If the log file type is specified as TXT, then any spaces in option
values need to be enclosed in quotes to avoid the shell parser parsing them as separate arguments, but also need to have enclosing quotes 
propagated to the logfile, or the spaces will be interpreted as delimiters in the log file.  e.g. - in the following example, enclosing escaped
quote characters ('\"') are added before adding the enclosing quotes::

    --notes-hdrs [headerstr1,"\"header str 2\"",headerstr3,...,headerStrN]


The shorthand notation is expanded to the notation described above (the COMPAS code just fills in the omitted values with the required defaults),
so the caveat mentioned above that option values may not begin with the dash character ('-') applies to the shorthand notation.

Shorthand notation is optional: users may choose to use the notation described above rather than shorthand notation, but in that case all option
values must be specified (no omissions, no defaults).

Mixing ranges and sets
======================

Ranges and sets can be specified together, and there is no limit to the number of ranges or sets that can be
specified on the command line, or in the :doc:`grid file <../grid-files>`.

Running COMPAS with the command::

    ./COMPAS --metallicity r[0.0001,5,0.0013] --common-envelope-alpha s[0.1,0.2,0.6,0.9]

would result in 20 binaries being evolved: 5 for the range of metallicities, times 4 for the set of CE alpha values.


Consider the following grid file, named `gridfile.txt`::

    --metallicity r[0.0001,5,0.0013] --common-envelope-alpha s[0.1,0.2,0.6,0.9]
    --fryer-supernova-engine s[rapid,delayed] --eccentricity r[0.0001,3,0.0003]

Running COMPAS with the command::

    ./COMPAS --grid gridfile.txt

would result in 26 binaries being evolved:

- 20 for the first grid line (5 for the range of metallicities, times 4 for the set of CE alpha values), and |br|
- 6 for the second grid line (2 for the set of Fryer SN engine values, and 3 for the range of eccentricities)


Running COMPAS with the command::

    ./COMPAS --remnant-mass-prescription s[mullermandel,fryer2012,hurley2000,muller2016] --grid gridfile.txt

would result in 104 binaries being evolved: the grid file would be ‘executed’ for each of the four remnant
mass prescriptions specified on the command line.

Multiple ranges and/or sets can be specified on the command line, and on each line of the grid file – so very
large numbers of stars/binaries can be evolved with just a few range/set specifications.
Program option sets
===================

A set of values can be specified for options of any data type that are not excluded from set specifications (see
note above).

Option value sets are specified by:

    --option-name set-specifier

where `set-specifier` is defined as:

    set-identifier[value\ :sub:`1` ,value\ :sub:`2` ,value\ :sub:`3` , ... ,value\ :sub:`n`]

    and

    .. list-table::
       :widths: 24 76 
       :header-rows: 0
       :class: aligned-text

       * - `set-identifier`
         - is one of {’s’, ’set’} (case is not significant)
       * - `value`:sub:`i`
         - is a value for the option

    Note that:

        `set-identifier` is mandatory for `set-specifier`. |br|
        `value`:sub:`i` must be the same data type as `option-name`.
        
There should be no spaces inside the brackets ([]). Spaces on the command line are interpreted as argument delimiters
by the shell parser before passing the command-line arguments to the COMPAS executable, so if spaces are present inside
the brackets the shell parser breaks the set specification into multiple command-line arguments.

Valid values for boolean options are {1|0, TRUE|FALSE, YES|NO, ON|OFF}, and all set values must be of
the same type (i.e. all 1|0, or all YES|NO etc.).

There is no limit to the number of values specified for a set, values can be repeated, and neither order nor
case is significant.

To specify a set of values for the ``--eccentricity-distribution`` option, a user, if running COMPAS from the
command line and with no grid file, would type any of the following::

    ./COMPAS --eccentricity-distribution s[THERMALISED,FIXED,FLAT]

    ./COMPAS --eccentricity-distribution set[THERMALISED,FIXED,FLAT]

In each of the examples above the user has specified, by the use of the `set-specifier`, that three binary stars
should be evolved, using the eccentricity distributions ’THERMALISED’, ’FIXED’, and ’FLAT’.

Note that when a set is, or sets are, specified on the command line, the ``--number-of-systems`` command-line option is ignored.
This is to avoid multiple systems with identical initial values being evolved.  Ranges and sets can be mixed with grid files, and
in that case ranges and sets specified on the command line will be played out for each grid file line.

   Post-processing tools
=====================

The COMPAS suite includes some useful post-processing tools that are located in the `postProcessing` directory.


.. toctree::
   :maxdepth: 1

   HDF5: Basics and COMPAS command line tools <post-processing-hdf5>
   Jupyter Notebooks: Working with HDF5, Data Analysis, and Cosmic Integration <../../../notebooks/Overview.ipynb>
HDF5 basics
===========

Here we provide some basic information regarding the COMPAS ``HDF5`` log files, and how COMPAS interacts with ``HDF5``.
Interested readers can learn more about the ``HDF5`` file format at `The HDF Group <https://www.hdfgroup.org/>`__.


.. _HDF5-chunking:

HDF5 file chunking and IO
-------------------------

Following is a brief description of ``HDF5`` files and ``chunking``, in the COMPAS context.

Data in ``HDF5`` files are arranged in ``groups`` and ``datasets``:

    A COMPAS output file (e.g. `BSE_System_Parameters`, `BSE_RLOF`, etc.) maps to an ``HDF5 group`` where the group name is
    the name of the COMPAS output file.

    A column in a COMPAS output file (e.g. `SEED`, `Mass(1)`, `Radius(2)`, etc.) maps to an ``HDF5 dataset``, where the
    dataset name is the column heading string.

    COMPAS column datatype strings are encoded in the ``HDF5 dataset`` meta-details (``dataset.dtype``).

    COMPAS column units strings are attached to ``HDF5 datasets`` as ``attributes``.

Each dataset in an ``HDF5`` file is broken into `chunks`, where a chunk is defined as a number of dataset entries. In COMPAS,
all datasets are 1-d arrays (columns), so a chunk is defined as a number of values in the 1-d array (or column). Chunking can 
be enabled or not, but if chunking is not enabled a dataset cannot be resized - so if chunking is not enabled the size of the 
dataset must be known at the time of creation, and the entire datset created in one go. That doesn't work for COMPAS - even 
though we know the number of systems being evolved, we don't know the number of entries we'll have in each of the output log 
files (and therefore the ``HDF5`` datasests if we're logging to ``HDF5`` files).  So, we enable chunking.
 
Chunking can improve, or degrade, performance depending upon how it is implemented - mostly related to the chunk size chosen.
 
``Datasets`` are stored inside an ``HDF5`` file as a number of chunks - the chunks are not guaranteed (not even likely) to be
contiguous in the file or on the storage media (HDD, SSD etc.). Chunks are mapped/indexed in the ``HDF5`` file using a B-tree,
and the size of the B-tree, and therefore the traversal time, depends directly upon the number of chunks allocated for a 
dataset - so the access time for a chunk increases as the number of chunks in the dataset increases. So many small chunks will
degrade performance.
 
``Chunks`` are the unit of IO for ``HDF5`` files - all IO to ``HDF5`` is performed on the basis of chunks. This means that 
whenever dataset values are accessed (read or written (i.e. changed)), if the value is not already in memory, the entire chunk
containing the value must be read from, or written to, the storage media - even if the dataset value being accessed is the only
value in the chunk. So few large chunks could cause empty, "wasted", space in the ``HDF5`` files (at the end of datasets) - but
they could also adversely affect performance by causing unecessary IO traffic (although probably not much in the way we access 
data in COMPAS files).
 
``HDF5`` files implement a chunk cache on a per-dataset basis. The default size of the chunk cache is $\small{1MB}$, and its 
maximum size is $\small{32MB}$. The purpose of the chunk cache is to reduce storage media IO - even with SSDs, memory access is 
much faster than storage media access, so the more of the file data that can be kept in memory and maipulated there, the better.  
Assuming the datatype of a particular dataset is ``DOUBLE``, and therefore consumes $\small{8}$ bytes of storage space, at its 
maximum size the chunk cache for that dataset could hold $\small{4,000,000}$ values - so a single chunk with $\small{4,000,000}$ 
values, two chunks with $\small{2,000,000}$ values, four with $\small{1,000,000}$, and so on. Caching a single chunk defeats the 
purpose of the cache, so chunk sizes somewhat less that $\small{4,000,000}$ would be most appropriate if the chunk cache is to be 
utilised. Chunks too big to fit in the cache simply bypass the cache and are read from, or written to, the storage media directly.
 
However, the chunk cache is really only useful for random access of the dataset. Most, if not all, of the access in the COMPAS 
context (including post-creation analyses) is serial - the COMPAS code writes the datasets from top to bottom, and later analyses
(generally) read the datasets the same way. Caching the chunks for serial access just introduces overhead that costs memory (not 
much, to be sure: up to $\small{32MB}$ per open dataset), and degrades performace (albeit it a tiny bit). For that reason we disable
the chunk cache in COMPAS - so all IO to/from an ``HDF5`` file in COMPAS is directly to/from the storage media. (To be clear, 
post-creation analysis software can disable the cache or not when accessing ``HDF5`` files created by COMPAS - disabling the cache
here does not affect how other software accesses the files post-creation).
 
So many small chunks is not so good, and neither is just a few very large chunks. So what's the optimum chunk size? That depends 
upon several things, and probably the most important of those are the final size of the dataset and the access pattern.
 
As mentioned above, we tend to access datasets serially, and generally from top to bottom, so larger chunks would seem appropriate, 
but not so large that we generate ``HDF5`` files with lots of unused space. However, disk space, even SSD space, is cheap, so 
trading space against performance is probably a good trade.
 
Also as mentioned above, we don't know the final size of (most of) the datasets when creating the ``HDF5`` in COMPAS - though we do
know the number of systems being generated, which allows us to determine an upper bound for at least some of the datasets (though 
not for groups such as `BSE_RLOF`).
 
One thing we need to keep in mind is that when we create the ``HDF5`` file we write each dataset of a group in the same 
iteration - this is analagous to writing a single record in (e.g.) a ``CSV`` log file (the ``HDF5 group`` corresponds to the ``CSV``
file, and the ``HDF5 datasets`` in the group correspond to the columns in the ``CSV`` file). So for each iteration - typically each
system evolved (though each timestep for detailed output files) we do as many IOs to the ``HDF5`` file as there are datasets in the
group (columns in the file). We are not bound to reading or writing a single chunk at a time - but we are bound to reading or writing
an integral multiple of whole chunks at a time.
 
We want to reduce the number of storage media accesses when writing (or later reading) the ``HDF5`` files, so larger chunk sizes are 
appropriate, but not so large that we create excessively large ``HDF5`` files that have lots of unused space (bearing in mind the 
trade-off mentioned above), especially when we're evolving just a few systems (rather than millions).
 
To really optimise IO performance for ``HDF5`` files we'd choose chunk sizes that are close to multiples of storage media block sizes,
but that would be too problematic given the number of disparate systems COMPAS could be run on...
 
Based on everything written above, and some tests, we've chosen a default chunk size of $\small{100,000}$ (dataset entries) for all 
datasets (``HDF5_DEFAULT_CHUNK_SIZE`` in ``constants.h``) for the COMPAS C++ code. This clearly trades performance against storage space.
For the (current) default logfile record specifications, per-binary logfile space is about $\small{1K}$ bytes, so in the very worst case
we will waste some space at the end of a COMPAS ``HDF5`` output file, but the performance gain, especially for post-creation analyses, is 
significant. Ensuring the number of systems evolved is an integral multiple of this fixed chunk size will minimise storage space waste.

We have chosen a minimum chunk size of $\small{1,000}$ (``HDF5_MINIMUM_CHUNK_SIZE`` in ``constants.h``) for the COMPAS C++ code. If the
number of systems being evolved is not less than ``HDF5_MINIMUM_CHUNK_SIZE`` the chunk size used will be the value of the ``hdf5-chunk-size``
program option (either ``HDF5_DEFAULT_CHUNK_SIZE`` or a value specified by the user), but if the number of systems being evolved is less
than ``HDF5_MINIMUM_CHUNK_SIZE`` the chunk size used will be ``HDF5_MINIMUM_CHUNK_SIZE``. This is just so we don't waste too much storage 
space when running small tests - and if they are that small performance is probably not going to be much of an issue, so no real trade-off 
against storage space.  


Copying and concatenating HDF5 files with h5copy.py
---------------------------------------------------

The chunk size chosen for the COMPAS C++ code determines the chunk size of the logfiles produced by the COMPAS C++ code.  If those files
are only ever given to programs such as ``h5copy`` as input files, their chunk size only matters in that it affects the read performance 
of the files (the more chunks, and the more smaller chunks, in a dataset of an input file means locating the chunks and reading them takes 
longer).  That may not be a huge problem depending upon how many input files there are and how big they are.  If a COMPAS logfile is used 
as a base file and other files are being appended to it via ``h5copy``, then the chunk size of the base output file will be the chunk size 
used for writing to the file - that could affect performance if it is too small.  We provide command-line options to specify the chunk size 
in both ``h5copy`` and the COMPAS C++ code so that users have some control over chunksize and performance.

Writing to the output HDF5 file is buffered in both ``h5copy`` and the COMPAS C++ code - we buffer a number of chunks for each open dataset 
and write the buffer to the file when the buffer fills (or a partial buffer upon file close if the buffer is not full). This IO buffering is 
not ``HDF5`` or filesystem buffering - this is a an internal implementation of ``h5copy`` and the COMPAS C++ code to improve performance.  
The IO buffer size can be changed via command-line options in both ``h5copy`` and the COMPAS C++ code.
 
Users should bear in mind that the combination of ``HDF5`` chunk size and ``HDF5`` IO buffer size affect performance, storage space, and 
memory usage - so they may need to experiment to find a balance that suits their needs.


A note on string values
-----------------------

COMPAS writes string data to its ``HDF5`` output files as C-type strings.  Python interprets C-type strings in ``HDF5`` files as byte 
arrays - regardless of the specified datatype when written (COMPAS writes the strings as ASCII data (as can be seen with ``h5dump``), but 
Python ignores that).  Note that this affects the values in datasets (and attributes) only, not the dataset names (or group names,
attribute names, etc.).

The only real impact of this is that if the byte array is printed by Python, it will be displayed as (e.g.) "b'abcde'" rather than just 
"abcde".  All operations on the data work as expected - it is just the output that is impacted.  If that's an issue, use ``.decode('utf-8')``
to decode the byte array as a Python string variable.

For example::

    str = h5File[Group][Dataset][0]

Here str is a byte array and ``print(str)`` will display (e.g.) ``b'abcde'``, but::

    str = h5File[Group][Dataset][0].decode('utf-8')

Here str is a Python string and ``print(str)`` will display (e.g.) ``abcde``

Note: ``HDF5`` files not created by COMPAS will not (necessarily) exhibit this behaviour, so for ``HDF5`` files created by the existing 
post-processing Python scripts the use of ``.decode()`` is not only not necessary, it will fail (because the strings in ``HDF5`` files 
created by Python are already Python strings, not byte arrays).

h5copy.py
=========

This program copies ``COMPAS_Output.h5`` ``HDF5`` file(s) [but not ``Detailed_Ouput`` files] to a designated output ``HDF5`` file. 
If the output file is an existing ``HDF5`` file, the user can specify whether the existing content should be erased before copying 
begins, or whether the copied data should be appended to the existing data. If multiple files are given as input files, the 
resultant ``HDF5`` file is the concatenation of the input files.


Some nomenclature
-----------------

Data in ``HDF5`` files are arranged in ``groups`` and ``datasets``:

    A COMPAS output file (e.g. `BSE_System_Parameters`, `BSE_RLOF`, etc.) maps to an ``HDF5 group`` where the group name is
    the name of the COMPAS output file.

    A column in a COMPAS output file (e.g. `SEED`, `Mass(1)`, `Radius(2)`, etc.) maps to an ``HDF5 dataset``, where the
    dataset name is the column heading string.

    COMPAS column datatype strings are encoded in the ``HDF5 dataset`` meta-details (``dataset.dtype``).

    COMPAS column units strings are attached to ``HDF5 datasets`` as ``attributes``.


h5copy usage
------------

::

    h5copy.py [-h] [-b BUFFER_SIZE] [-c CHUNK_SIZE] [-e] [-f FILENAME_FILTER]
                   [-o OUTPUT_FILENAME] [-r [RECURSION_DEPTH]] [-s]
                   [-x EXCLUDE_GROUP [EXCLUDE_GROUP ...]]
                   input [input ...]

    HDF5 file copier.

    positional arguments:
      input
        input directory and/or file name(s)

    optional arguments:
      -h, --help
        show this help message and exit
      -b BUFFER_SIZE, --buffer-size BUFFER_SIZE
        IO buffer size (number of HDF5 chunks, default = 10)
      -c CHUNK_SIZE, --chunk-size CHUNK_SIZE
        HDF5 output file dataset chunk size (default = 100000)
      -e, --erase-output
        erase existing output file before copying input files (default = False)
      -f FILENAME_FILTER, --filter FILENAME_FILTER
        input filename filter (default = *)
      -o OUTPUT_FILENAME, --output OUTPUT_FILENAME
        output file name (default = h5out.h5)
      -r [RECURSION_DEPTH], --recursive [RECURSION_DEPTH]
        recursion depth (default is no recursion)
      -s, --stop-on-error
        stop all copying if an error occurs (default is skip to next file and continue)
      -x EXCLUDE_GROUP [EXCLUDE_GROUP ...], --exclude EXCLUDE_GROUP [EXCLUDE_GROUP ...]
        list of input groups to be excluded (default is all groups will be copied)


    Note: if the -x option is specified, it should be specified at the end of the options 
          list (i.e. the list of input files can't follow the -x option or the list of input 
          files will be subsumed by the list of groups to be excluded)


h5copy functionality overview
-----------------------------

Output file
~~~~~~~~~~~

- If the specified ``HDF5`` output file does not exist, it will be created, and data from the input file(s) will be appended 
  to the new output file.

- If the specified ``HDF5`` output file does exist, and the ``--erase-ouput [-e]`` command-line option is not specified, the 
  existing output file will be preserved and data from the input file(s) will be appended to the existing output file if the 
  output file was created with chunking enabled (see :ref:`HDF5-chunking` - only files that were created with chunking enabled 
  can be extended).

  Attempting to append data to an existing ``HDF5`` file that was not created with chunking enabled will result in the following 
  message being displayed for each dataset::

      Only chunked datasets can be resized

  If that happens, using ``h5copy.py`` to copy the output file to a new file will copy the existing data and create the new file 
  with chunking enabled - the newly created file can then be used as a base file to which other files can be appended.

- If the specified ``HDF5`` output file does exist, and the ``--erase-ouput [-e]`` command-line option is specified, the existing
  output file will be deleted, and new, empty, file created, and data from the input file(s) will be appended to the new output file.


##############
Appending data
##############

   (a) if a group in an input file already exists in the output file, then the group data (the datasets within the group) will only 
       be copied if the number of datasets in the input file group matches the number of datasets in the output file group. If there 
       is a mismatch a warning will be issued and the group will not be copied (but the file copy will continue, just as though the 
       group had been excluded) 
   
   (b) if a dataset in an input file does not exist in the output file, the dataset will be created, otherwise the data copied will 
       be appended to the existing dataset.


Input files
~~~~~~~~~~~

A list of input filenames and or directory names must be supplied. The list can be a single name. Each name in the list is processed 
in order. If a name is the name of a file with the file extension `.h5`, the contents of the file will be copied to the output file, 
otherwise it will be ignored. If a name is the name of a directory and the specified recursion level requires that the directory be 
processed (see command-line option ``--recursive [-r]``), the program will descend into the directory and process all files and 
directories there (directories will be processed depending upon the value of the ``--recursive [-r]`` option), otherwise it will be 
ignored.

The command-line option ``--recursive [-r]`` specifies whether recursion is enabled for directory processing, and if it is, to what 
depth:

    - If the ``--recursive [-r]`` option is not specified, recursion is not enabled and only files in the specified working directory 
      will be candidates for copying.

    - if ``--recursive [-r]`` is specified with no ``depth`` value, recursion is enabled and the depth is not limited - that is, all 
      files in the specified working directory, and all files in all directories below the specified working directory, will be 
      candidates for copying.
    
    - If ``--recursive [-r]`` is specified with a specified ``depth`` value, recursion is enabled and the depth is limited to the
      depth specified - that is, all files in the specified working directory, and all files in all directories `depth` levels below
      the specified working directory, will be candidates for copying.
    

#####################
Input filename filter
#####################

If the ``--filter [-f]`` command-line option is specified, the names of all candidate files will be checked against the specified filter,
and only files whose names match the filter will be copied.  The specified filter is a filename-only filter - the file's path (i.e. its
location) will not be matched to the filter. The specified filter should not include a file extension, but the program adds the extension
`.h5` to the specified filter - only files that have the file extension `.h5` will match the filter.

If ``--filter [-f]`` is not specified, the program uses a default filter value of `*`, then adds the `.h5` file extension - so all candidate
files with the `.h5` extension will be copied.


Excluding HDF5 groups
~~~~~~~~~~~~~~~~~~~~~

If the ``--exclude [-x]`` command-line option is specified, the specified list of groupnames will be excluded from data copied from all input
files.  If ``--exclude [-x]`` is not specified, all groups in all candidate files will be copied.


Erase output
~~~~~~~~~~~~

If the ``--erase-ouput [-e]`` command-line option is specified and the output file (specified or default) exists, it will be erased before 
copying begins.  The ``--erase-ouput [-e]`` command-line option is ignored if the output file does not exist.

If ``--erase-ouput [-e]`` is not specified and the output file (specified or default) exists, the existing content will be preserved and any
data copied to the file will be appended to the existing data.

Processing HDF5 files
=====================

One of the log file formats COMPAS provides is the ``HDF5``\ [#f1]_ file format.  Since ``HDF5`` files are not
human-readable, we provide some tools to help with processing the ``HDF5`` log files produced by COMPAS.


.. toctree::
   :maxdepth: 1

   post-processing-hdf5-info.rst 
   Copying and concatenating: h5copy.py <post-processing-h5copy>
   View summary, header details, and contents: h5view.py <post-processing-h5view>

.. rubric:: Footnotes

.. [#f1] https://www.hdfgroup.org/
h5view.py
=========

This program displays summary, header, and content information for specified COMPAS ``HDF5`` file(s). It's fairly rudimetary - the 
``HDF5`` package provides ``h5dump`` and ``h5ls`` which have far more options than this program - but this program is somewhat 
COMPAS aware.


h5view usage
------------

::

    h5view.py [-h] [-f FILENAME_FILTER] [-r [RECURSION_DEPTH]] [-S] [-H]
                   [-C [CONTENTS]] [-s] [-x EXCLUDE_GROUP [EXCLUDE_GROUP ...]]
                   [-V SEED_LIST [SEED_LIST ...]]
                   input [input ...]

    HDF5 file content viewer.

    positional arguments:
      input
        input directory and/or file name(s)

    optional arguments:
      -h, --help
        show this help message and exit
      -f FILENAME_FILTER, --filter FILENAME_FILTER
        input filename filter (default = *)
      -r [RECURSION_DEPTH], --recursive [RECURSION_DEPTH]
        recursion depth (default is no recursion)
      -S, --summary
        display summary output for HDF5 file (default is not to displat summary)
      -H, --headers
        display file headers for HDF5 file (default is not to display headers)
      -C [CONTENTS], --contents [CONTENTS]
        display file contents for HDF5 file: argument is number of entries (+ve from top, -ve
        from bottom) (default is not to display contents)
      -s, --stop-on-error
        stop all copying if an error occurs (default is skip to next file and continue)
      -x EXCLUDE_GROUP [EXCLUDE_GROUP ...], --exclude EXCLUDE_GROUP [EXCLUDE_GROUP ...]
        list of input groups to be excluded (default is all groups will be copied)
      -V SEED_LIST [SEED_LIST ...], --seeds SEED_LIST [SEED_LIST ...]
        list of seeds to be printed (for content printing) (default is print all seeds)


Example
-------

Typing::

    python3 h5view.py compas-output-file.h5
    
will result in summary output of the ``HDF5`` file `compas-output-file.h5` that looks something like this::

    Summary of HDF5 file /d/compas/h5out.h5
    =======================================

    File size    : 2.1520 GB
    Last modified: 2021-07-26 16:25:12.928401

    COMPAS Filename              Columns   Entries   Unique SEEDs
    --------------------------   -------   -------   ------------
    Run_Details                      346        30
    BSE_Common_Envelopes              73    582485         476514
    BSE_Double_Compact_Objects        12      8725           8725
    BSE_RLOF                          34   2997481         536332
    BSE_Supernovae                    32    103000          87162
    BSE_Switch_Log                    13   7472234         956623
    BSE_System_Parameters             33   1050000        1050000


Other ``h5view.py`` options (listed above) display headers and file contents.



h5view functionality overview
-----------------------------

``h5view.py`` displays summary, header and content information for specified COMPAS ``HDF5`` file(s). If none of the command-line
options ``--summary [-S]``, ``--headers [-H]``, or ``--contents [-C]`` are specified, ``--summary [-s]`` is assumed. If any of 
``--summary [-S]``, ``--headers [-H]``, or ``--contents [-C]`` are specified, then only the option(s) specified are actioned.

Displaying summary information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Summary information displays, for each COMPAS file in the ``HDF5`` file:
   - the name of the COMPAS file
   - the number of columns in the COMPAS file
   - the number of entries in the COMPAS file (actually, the maximum number of entries in any column in the COMPAS file)


Displaying header information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Header information displays, for each COMPAS file in the ``HDF5`` file:
   - the name of each column in the COMPAS file
   - the number of entries in each column of the COMPAS file
   - the data type of each column of the COMPAS file
   - the units associated with each column of the COMPAS file
     (with the exception of the ``Run_Details`` file - there are no units associated with columns in the ``Run_Details`` file)


Displaying contents
~~~~~~~~~~~~~~~~~~~

Contents information displays, for each COMPAS file in the ``HDF5`` file:
   - a header showing the column names in the COMPAS file
   - a row for each entry in the COMPAS file, showing the column values for that row (comma delimited)

   The contents display can be limited in two ways:

      (a) The ``--contents [-C]`` option takes and optional argument: an integer number of rows to display. The argement to 
          ``--contents [-C]`` can be positive or negative: a positive value indicates that the number of rows specified by the 
          argument should be displayed from the start of the file; a negative value indicates that the number of rows specified
          by the (absolute value of the) argument should be displayed from the end of the file.  The +ve and -ve arguments to 
          the ``--contents [-C]`` option are akin the the Unix ``head`` and ``tail`` commands.

      (b) The ``--seeds [-V]`` option allows the user to specify a list of SEED values that should be printed. If the 
          ``--seeds [-V]`` option is specified, only rows containing the seeds specified by the user will be printed - and only 
          if they are in the entries printed if limited by the ``--contents [-C]`` argument  described in (a).

          Note that printing only seeds specified in a list of seeds could be slow - we effectively have to look through the 
          entire dataset looking for the seeds required.

Tutorial: simple COMPAS run
===========================

NOTE: we are currently updating our documentation and will include a ``compasConfigDefault.yaml`` for the demo asap.

..
    This tutorial assumes that you have already built the COMPAS executable as described in :doc:`../../Getting started/building-COMPAS`.

    For this example you will need the python script ``pythonSubmitDemo.py``, which specifies all the program options (physics assumptions, 
    output types) and runs COMPAS in the terminal. Although the primary functionality of COMPAS is to evolve a whole population of binary 
    stars rapidly, for now, let's focus on evolving a single stellar system and examining the detailed output.

    If you haven't yet defined the ``COMPAS_ROOT_DIR`` environment variable, do that now::

        export COMPAS_ROOT_DIR=path-to-compas

    where `path-to-compas` should be replaced with the path to the parent directory of the COMPAS `src` directory. Depending upon your system,
    for the ``export`` command to take effect, it may be necessary to either restart your session or execute the following command::

        source ~/.bashrc

    To start, change to the ``examples/methods_paper_plots/detailed_evolution/`` directory::

      cd $COMPAS_ROOT_DIR/examples/methods_paper_plots/detailed_evolution/

    where you will find the script ``pythonSubmitDemo.py`` for this demo.


.. toctree::
   :maxdepth: 1

   ./example-compas-run-grid
   ./example-compas-run-detailed-output
   
Running COMPAS using a grid file
================================

In population synthesis, the initial stellar population is usually generated by drawing the primary mass, secondary mass, semi-major axis, 
and eccentricity from their respective distributions specified in the program options. However, we illustrate COMPAS's ability to specify 
a grid of initial values for single and binary star evolution using COMPAS's grid functionality.

An example grid file, ``Grid_demo.txt``, is included in the ``detailed_evolution`` directory. Open it with a text editor to view it::

    # Demo BSE Grid file

    --initial-mass-1 35.4 --initial-mass-2 29.3 --metallicity 0.001  --eccentricity 0.000000e+00 --semi-major-axis 1.02

It should be clear that this grid file specifies a binary of zero-age main sequence stars with primary mass 
35.4\ :math:`\small M_\odot`, secondary mass 29.3\ :math:`\small M_\odot`, metallicity 0.001, zero eccentricity, and semi-major axis of 
1.02AU. See :doc:`../grid-files` for detailed information regarding COMPAS's grid functionality for both single and binary stars.

We will execute COMPAS via the ``runSubmit.py`` script, but first we need to edit the companion ``compasConfigDefault.yaml`` script to instruct COMPAS to read the grid file
(via the ``grid`` program option).

Open ``$COMPAS_ROOT_DIR/preProcessing/compasConfigDefault.yaml`` with a text editor, and specify the grid filename::

    grid_filename = 'Grid_demo.txt'
    
Note the quotes around the filename. 

If the filename specified is not fully-qualified, and the shell environment variable ``COMPAS_INPUT_DIR_PATH`` exists and is not empty,
the value of ``COMPAS_INPUT_DIR_PATH`` will be prepended to the specified grid filename. 


To print the detailed evolution of binary properties over time, we need to turn on detailed output, by specifying::

    detailed_output = True

in ``compasConfigDefault.yaml``.

COMPAS can produce logfiles of different types: ``HDF5``, ``CSV``, ``TSV``, and ``TXT``, which can be chosen by editing the line::

    logfile_type = 'HDF5'

in ``compasConfigDefault.yaml``. The default type is ``HDF5`` - we'll leave the default.

NOTE: we are currently updating our documentation and will include a ``compasConfigDefault.yaml`` for the demo asap.

..
    For this turtorial, this has all been done for you in the file ``pythonSubmitDemo.py`` found in the ``examples/methods_paper_plots/detailed_evolution/`` directory.

..
    Now let's run COMPAS!

..
    ::

        $ python3 pythonSubmitDemo.py

        COMPAS v02.18.06
        Compact Object Mergers: Population Astrophysics and Statistics
        by Team COMPAS (http://compas.science/index.html)
        A binary star simulator

        Start generating binaries at Thu Feb 25 14:42:05 2021

        Evolution of current binary stopped: Double compact object
        0: Evolution stopped: (Main_Sequence_>_0.7 -> Black_Hole) + (Main_Sequence_>_0.7 -> Black_Hole)

        Generated 1 of 1 binaries requested

        Simulation completed

        End generating binaries at Thu Feb 25 14:42:05 2021

        Clock time = 0.108338 CPU seconds
        Wall time  = 00:00:00 (hh:mm:ss)

..
    Congratulations! You've just made a binary black hole. And it didn't even take a second.
Examining detailed output
=========================

The COMPAS run from the tutorial creates a new directory ``COMPAS_Output``, inside which you will find the following files/directories 
(here we assume ``logfile_type = 'HDF5'`` in the python submit file):

**Run_Details** |br|
A record of the COMPAS command-line program options specified for this tutorial (these are the values set by ``compasConfigDefault.yaml`` when using 
``runSubmit.py``, or the COMPAS default values if not executing via ``runSubmit.py``).

**COMPAS_Output.h5** |br|
The primary output file, containing ``HDF5`` data groups for the relevant output physics. By default, and for a sufficiently large simulation, 
this will include:

    - BSE_Common_Envelopes
    - BSE_Double_Compact_Objects
    - BSE_RLOF
    - BSE_Supernovae
    - BSE_System_Parameters

**Detailed_Output** |br|
This directory contains the detailed output file, ``BSE_Detailed_Output_0.h5``, which records the detailed time evolution of binary. 
This file, and directory, is only produced if ``detailed_output = True`` in the python submit file.

We examine ``BSE_Detailed_Output_0.h5`` to look at the evolution of the two stars. A default python plotting script has been included to 
visualise the data. Let's run the script::

  python3 detailed_evol_plotter.py

This should produce the plot shown in :ref:`Figure 5 <fig-5>`:

.. _fig-5:

.. figure:: ../../../images/example-plot-compressed.svg
    :width: 1100px
    :height: 625px
    :align: center
    :figclass: align-center
    :alt: Example COMPAS run demo plot

    Figure 5 Example COMPAS run.

COMPAS provides many tools for analysing and post-processing the data - see :doc:`../Post-processing/post-processing` for more details.

Installing the COMPAS Docker image
----------------------------------

COMPAS ``Docker`` images for all releases of COMPAS are hosted on ``dockerHub``, and tagged\ [#f1]_ with the COMPAS version number.

The latest COMPAS ``Docker`` compiled version of COMPAS (dev branch) can be retrieved by executing either of the following commands::

    docker pull teamcompas/compas
    docker pull teamcompas/compas:latest

Other versions can be retrieved by specifying the tag that corresponds to the COMPAS version required. For example, to retrieve the
image for COMPAS version 2.12.0, type::

    docker pull teamcompas/compas:2.12.0


To see all available versions, go to the `TeamCOMPAS docker hub page <https://hub.docker.com/u/teamcompas>`_.


.. rubric:: Footnotes

.. [#f1] https://docs.docker.com/engine/reference/commandline/tag/
Running COMPAS from the command line
====================================

COMPAS is a command-line application.  Interaction with COMPAS is entirely through the terminal and shell - there is no
visual or graphical user interface (GUI).  COMPAS interacts with the user by accepting input via the keyboard, and providing
ouput by writing plain text to the terminal. COMPAS reads input files where necessary: a ``grid`` file (see :doc:`../grid-files`),
and a log file definitions file (see :doc:`../COMPAS output/standard-logfiles-record-specification`), and produces output files
(see :doc:`../COMPAS output/output`), but these are not interactive.

Command-line applications accept interactive input from the user in a number of ways: one of those is via command-line switches and 
arguments, or, more generally, command-line options. This is the method COMPAS uses to interact with the user.  

A few example COMPAS runs, using a small sample of available program option, are shown below. For detailed information regarding program 
option use, and a full list of program options available including their default values see :doc:`../Program options/program-options`.


To run COMPAS from the command line, and have COMPAS output its version string, type either of::

    ./compas -v
    ./compas --version

This should produce output that looks something like::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator


To run COMPAS from the command line, and have COMPAS output a list of the command-line options available, type::

    ./compas -h

Note the single dash before the option name.  This should produce output that looks something like::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Options:
    --PISN-lower-limit
    --PISN-upper-limit
    --PPI-lower-limit
    --PPI-upper-limit
    --add-options-to-sysparms
    --allow-rlof-at-birth
    --allow-touching-at-birth
    --angular-momentum-conservation-during-circularisation
    --black-hole-kicks
    --case-BB-stability-prescription
    --check-photon-tiring-limit
    --chemically-homogeneous-evolution

    ...

The options are listed in `ASCIIbetical` order (digits come before uppercase alpha characters, which come before lowercase alpha characters).


To run COMPAS from the command line, and have COMPAS output a list of the command-line options available, together with a brief description
of  the option, and its default value in the COMPAS code, type::

    ./compas --help

Note two dashes before the option name.  This should produce output that looks something like::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Options:
    --PISN-lower-limit
      Minimum core mass for PISN (default = 60.000000)
    --PISN-upper-limit
      Maximum core mass for PISN (default = 135.000000)
    --PPI-lower-limit
      Minimum core mass for PPI (default = 35.000000)
    --PPI-upper-limit
      Maximum core mass for PPI (default = 60.000000)
    --add-options-to-sysparms
      Add program options columns to BSE/SSE SysParms file (options: [ALWAYS, GRID, NEVER], default = GRID)
    --allow-rlof-at-birth
      Allow binaries that have one or both stars in RLOF at birth to evolve (default = TRUE)
    --allow-touching-at-birth
      Allow binaries that are touching at birth to evolve (default = FALSE)
    --angular-momentum-conservation-during-circularisation
      Conserve angular momentum when binary is circularised when entering a Mass Transfer episode (default = FALSE)
    --black-hole-kicks
      Black hole kicks relative to NS kicks (options: [FULL, REDUCED, ZERO, FALLBACK], default = FALLBACK)
    --case-BB-stability-prescription
      Case BB/BC mass transfer stability prescription (options: [ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE], default = ALWAYS_STABLE)
    --check-photon-tiring-limit
      Check the photon tiring limit hasn't been exceeded by wind mass loss (default = FALSE)
    --chemically-homogeneous-evolution
      Chemically Homogeneous Evolution (options: [NONE, OPTIMISTIC, PESSIMISTIC], default = PESSIMISTIC)

    ...

Again, the options are listed in `ASCIIbetical` order.


To run a default COMPAS run of 10 binary systems with default initial conditions and evolutionary parameters, type::

    ./compas

This should produce an output put similar to::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Start generating binaries at Tue Sep  7 17:34:49 2021

    0: Stars merged: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    1: Stars merged: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    2: Double White Dwarf: (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf) + (Main_Sequence_>_0.7 -> Helium_White_Dwarf)
    3: Stars merged: (Main_Sequence_>_0.7 -> Naked_Helium_Star_MS) + (Main_Sequence_<=_0.7 -> Main_Sequence_<=_0.7)
    4: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    5: Unbound binary: (Main_Sequence_>_0.7 -> Naked_Helium_Star_MS) + (Main_Sequence_>_0.7 -> Neutron_Star)
    6: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    7: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Core_Helium_Burning)
    8: Allowed time exceeded: (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf) + (Main_Sequence_>_0.7 -> Neutron_Star)
    9: Unbound binary: (Main_Sequence_>_0.7 -> Black_Hole) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)

    Generated 10 of 10 binaries requested

    Simulation completed

    End generating binaries at Tue Sep  7 17:34:49 2021

    Clock time = 0.078125 CPU seconds
    Wall time  = 0000:00:00 (hhhh:mm:ss)


To run COMPAS and evolve five binary system with somespecific initial conditions and evolutionary parameters (default values for the remainder), type::

    ./compas --number-of-systems 5 --initial-mass-1 8.5 --initial-mass-2 13.7 --metallicity 0.015 --mass-loss-prescription VINK --common-envelope-alpha 0.8 --common-envelope-lambda 0.2

This should produce an output put similar to::

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Start generating binaries at Tue Sep  7 17:49:40 2021

    0: Unbound binary: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Neutron_Star)
    1: Unbound binary: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Neutron_Star)
    2: Stars merged: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Naked_Helium_Star_MS)
    3: Unbound binary: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Neutron_Star)
    4: Unbound binary: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Neutron_Star)

    Generated 5 of 5 binaries requested

    Simulation completed

    End generating binaries at Tue Sep  7 17:49:40 2021

    Clock time = 0.0625 CPU seconds
    Wall time  = 0000:00:00 (hhhh:mm:ss)

Running COMPAS via Python
=========================

A convenient method of managing the many program options provided by COMPAS is to run COMPAS via Python, using a script to manage and 
specify the values of the program options.

An example Python script is provided in the COMPAS suite on github: ``runSubmit.py``. Additionally, the default COMPAS options are specified on ``compasConfigDefault.yaml``. Users should copy the ``runSubmit.py`` and ``runSubmit.py`` scripts and modify the ``compasConfigDefault.yaml`` copy to match their experimental requirements. Refer to the :doc:`Getting started guide <../../Getting started/getting-started>` for more details.

To run COMPAS via Python using the ``runSubmit.py`` script provided, set the shell environment variable ``COMPAS-ROOT-DIR``
to the parent directory of the directory in which the COMPAS executable resides, then type `python /path-to-runSubmit/runSubmit.py`. 
For example, for Ubuntu Linux, type::

    export COMPAS_ROOT_DIR=/path-to-dir

    python /path-to-runSubmit/runSubmit.py

This should produce an output similar to::

    python_version = 3
    ../src/COMPAS --enable-warnings False --use-mass-loss  --mass-transfer  --detailed-output False --evolve-unbound-systems False --population-data-printing False --rlof-printing  --circularise-binary-during-mass-transfer  --angular-momentum-conservation-during-circularisation False --pair-instability-supernovae  --pulsational-pair-instability  --quiet False --common-envelope-allow-main-sequence-survive  --evolve-pulsars False --debug-to-file False --errors-to-file False --allow-rlof-at-birth  --allow-touching-at-birth False --store-input-files  --switch-log False --check-photon-tiring-limit False --number-of-systems 10 --metallicity 0.0142 --common-envelope-alpha 1.0 --common-envelope-lambda 0.1 --common-envelope-slope-kruckow -0.8333333333333334 --common-envelope-alpha-thermal 1.0 --common-envelope-lambda-multiplier 1.0 --luminous-blue-variable-multiplier 1.5 --overall-wind-mass-loss-multiplier 1.0 --wolf-rayet-multiplier 1.0 --cool-wind-mass-loss-multiplier 1.0 --mass-transfer-fa 0.5 --mass-transfer-jloss 1.0 --maximum-evolution-time 13700.0 --maximum-number-timestep-iterations 99999 --timestep-multiplier 1.0 --initial-mass-min 5.0 --initial-mass-max 150.0 --initial-mass-power 0.0 --semi-major-axis-min 0.01 --semi-major-axis-max 1000.0 --mass-ratio-min 0.01 --mass-ratio-max 1.0 --minimum-secondary-mass 0.1 --eccentricity-min 0.0 --eccentricity-max 1.0 --metallicity-min 0.0001 --metallicity-max 0.03 --pulsar-birth-magnetic-field-distribution-min 11.0 --pulsar-birth-magnetic-field-distribution-max 13.0 --pulsar-birth-spin-period-distribution-min 10.0 --pulsar-birth-spin-period-distribution-max 100.0 --pulsar-magnetic-field-decay-timescale 1000.0 --pulsar-magnetic-field-decay-massscale 0.025 --pulsar-minimum-magnetic-field 8.0 --orbital-period-min 1.1 --orbital-period-max 1000 --kick-magnitude-sigma-CCSN-NS 265.0 --kick-magnitude-sigma-CCSN-BH 265.0 --fix-dimensionless-kick-magnitude -1 --kick-direction-power 0.0 --random-seed 0 --mass-transfer-thermal-limit-C 10.0 --eddington-accretion-factor 1 --pisn-lower-limit 60.0 --pisn-upper-limit 135.0 --ppi-lower-limit 35.0 --ppi-upper-limit 60.0 --maximum-neutron-star-mass 2.5 --kick-magnitude-sigma-ECSN 30.0 --kick-magnitude-sigma-USSN 30.0 --kick-scaling-factor 1.0 --maximum-mass-donor-nandez-ivanova 2.0 --common-envelope-recombination-energy-density 15000000000000.0 --common-envelope-mass-accretion-max 0.1 --common-envelope-mass-accretion-min 0.04 --zeta-main-sequence 2.0 --zeta-radiative-envelope-giant 6.5 --kick-magnitude-max -1.0 --muller-mandel-kick-multiplier-BH 200.0 --muller-mandel-kick-multiplier-NS 400.0 --log-level 0 --debug-level 0 --hdf5-chunk-size 100000 --hdf5-buffer-size 1 --neutrino-mass-loss-BH-formation-value 0.1 --mode BSE --case-BB-stability-prescription ALWAYS_STABLE --chemically-homogeneous-evolution PESSIMISTIC --luminous-blue-variable-prescription HURLEY_ADD --mass-loss-prescription VINK --mass-transfer-angular-momentum-loss-prescription ISOTROPIC --mass-transfer-accretion-efficiency-prescription THERMAL --mass-transfer-rejuvenation-prescription STARTRACK --initial-mass-function KROUPA --semi-major-axis-distribution FLATINLOG --orbital-period-distribution FLATINLOG --mass-ratio-distribution FLAT --eccentricity-distribution ZERO --metallicity-distribution ZSOLAR --rotational-velocity-distribution ZERO --remnant-mass-prescription FRYER2012 --fryer-supernova-engine DELAYED --black-hole-kicks FALLBACK --kick-magnitude-distribution MAXWELLIAN --kick-direction ISOTROPIC --output-path /d/Jeff/User_Files/compas/dev/my_fork/compas/src --common-envelope-lambda-prescription LAMBDA_NANJING --stellar-zeta-prescription SOBERMAN --mass-transfer-thermal-limit-accretor CFACTOR --pulsational-pair-instability-prescription MARCHANT --neutron-star-equation-of-state SSE --pulsar-birth-magnetic-field-distribution ZERO --pulsar-birth-spin-period-distribution ZERO --common-envelope-mass-accretion-prescription ZERO --envelope-state-prescription LEGACY --logfile-type HDF5 --neutrino-mass-loss-BH-formation FIXED_MASS

    COMPAS v02.22.00
    Compact Object Mergers: Population Astrophysics and Statistics
    by Team COMPAS (http://compas.science/index.html)
    A binary star simulator

    Start generating binaries at Tue Sep  7 18:14:40 2021

    0: Stars merged: (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    1: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Thermally_Pulsing_Asymptotic_Giant_Branch)
    2: Unbound binary: (Main_Sequence_>_0.7 -> Black_Hole) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    3: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    4: Unbound binary: (Main_Sequence_>_0.7 -> Black_Hole) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    5: Allowed time exceeded: (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf) + (Main_Sequence_<=_0.7 -> Main_Sequence_<=_0.7)
    6: Unbound binary: (Main_Sequence_>_0.7 -> Neutron_Star) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)
    7: Double White Dwarf: (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf) + (Main_Sequence_>_0.7 -> Carbon-Oxygen_White_Dwarf)
    8: Unbound binary: (Main_Sequence_>_0.7 -> Black_Hole) + (Main_Sequence_<=_0.7 -> Main_Sequence_<=_0.7)
    9: Stars merged: (Main_Sequence_>_0.7 -> Naked_Helium_Star_Giant_Branch) + (Main_Sequence_>_0.7 -> Main_Sequence_>_0.7)

    Generated 10 of 10 binaries requested

    Simulation completed

    End generating binaries at Tue Sep  7 18:14:40 2021

    Clock time = 0.109375 CPU seconds
    Wall time  = 0000:00:00 (hhhh:mm:ss)

Note that Python prints the Python version, the executes the command to run COMPAS.  The command exceuted is echoed to the stdout.  COMPAS
then runs and produces its usual output.

When using Python and a script file (such as `runSubmit.py`) to run COMPAS, care must be taken to specify program option values correctly in the ``compasConfigDefault.yaml`` file.
For example, ranges and sets can be specified for options in the ``compasConfigDefault.yaml`` file, but the range or set parameter must be enclosed in quotes – 
otherwise python tries to parse the construct. For example, to specify a set of metallicity values in the Python script file, use::

    metallicity = 's[0.001,0.002,0.003,0.007,0.01,0.015,0.02]'

If the set parameter is not enclosed in quotes, Python will attempt to parse it, and will fail.
Running COMPAS
==============

There are several ways to run COMPAS: each method suiting different needs:


.. toctree::
   :maxdepth: 1

   running-via-cmdline
   running-via-python
   running-via-docker

Running the COMPAS Docker image
===============================

COMPAS can be configured as usual via command line arguments passed to the COMPAS executable or via a ``runSubmit.py`` file in the 
``Docker`` environment.


via runSubmit.py
-------------------

To run COMPAS via a ``runSubmit.py`` file, type::

    docker run                                                  \
        --rm                                                    \
        -it                                                     \
        -v $(pwd)/compas-input:/app/COMPAS/config               \
        -v $(pwd)/compas-logs:/app/COMPAS/logs                  \
        -v $(pwd)/runSubmit.py:/app/starts/runSubmit.py   \
        -e COMPAS_EXECUTABLE_PATH=/app/COMPAS/bin/COMPAS        \
        -e COMPAS_INPUT_DIR_PATH=/app/COMPAS/config             \
        -e COMPAS_LOGS_OUTPUT_DIR_PATH=/app/COMPAS/logs         \
        teamcompas/compas                                       \
        python3 /app/starts/runSubmit.py                     


NOTE: if you decide to execute using ``runSubmit.py``, you will need 
a ``compasConfigDefault.yaml``  file in the same directory. This file 
can be find in the same directory as the ``runSubmit.py``, and contains
the default COMPAS choices for stellar and binary physics. These choices
can be changed by modifying the options availabe in the ``compasConfigDefault.yaml`` 
file.

Breaking down this command:

**docker run** |br|
creates a container.

**--rm** |br|
destroy the container once it finishes running the command\ [#f1]_.

**-it** |br|
short for [-i and -t] - provides an interactive terminal\ [#f2]_.

**-v <path-on-host>:<path-in-container>** |br|
mount ``<path-on-host>`` to ``<path-in-container>``\ [#f3]_. |br|

This time we not only want to read the COMPAS input files (i.e. grid file and/or logfile-definitions file) on the
host from the container, and get the output from COMPAS in the container onto the host machine, we also want to 
supply a ``runSubmit.py`` to the container from the host machine.

**-e VAR_NAME=value** |br|
set the environment variable ``VAR_VAME`` to `value`\ [#f4]_.

**teamcompas/compas** |br|
the image to run.

**python3 /app/starts/runSubmit.py** |br|
the command to run when the container starts.


via the command line
--------------------

To run the COMPAS executable from the command line (i.e. without ``runSubmit.py``), type::

    docker run                                      \
        --rm                                        \
        -it                                         \
        -v $(pwd)/compas-input:/app/COMPAS/config   \
        -v $(pwd)/compas-logs:/app/COMPAS/logs      \
        teamcompas/compas                           \
        bin/COMPAS                                  \
        --number-of-systems=5                       \
        --output-path=/app/COMPAS/logs


Breaking down this command:

**docker run** |br|
creates a container\ [#f5]_.

**--rm** |br|
destroy the container once it finishes running the command\ [#f1]_.

**-it** |br|
short for [-i and -t] - provides an interactive terminal\ [#f2]_.

**-v <path-on-host>:<path-in-container>** |br|
mount ``<path-on-host>`` to ``<path-in-container>``\ [#f3]_. |br|

In this instance, make it so |br|
   `$(pwd)/compas-input` on my machine is the same as `/app/COMPAS/config` inside the container. |br|
   `$(pwd)/compas-logs` on my machine is the same as `/app/COMPAS/logs` inside the container.

**teamcompas/compas** |br|
the image to run.

**bin/COMPAS** |br|
the command to run when the container starts.

**--number-of-systems=5** |br|
the flag to set the number of binaries.

**--output-path=/app/COMPAS/logs** |br|
forces logs to go to the directory that is mapped to the host machine.



Environment variables
---------------------

Three new environment variables are used in ``runSubmit.py``.  These environment variables are used primarily in the ``Docker``
environment, and are non-breaking changes (i.e. benign to other environments).

``COMPAS_EXECUTABLE_PATH`` specifies where ``runSubmit.py`` looks for the COMPAS executable. This override exists purely for 
ease-of-use from the command line.

`COMPAS_LOGS_OUTPUT_DIR_PATH` specifies where COMPAS output log files are created. The override exists because the mounted directory 
(option `-v`) is created before COMPAS runs. COMPAS sees that the directory where it's supposed to put logs already exists, so it 
creates a different (i.e. non-mapped) directory for the output log files.

`COMPAS_INPUT_DIR_PATH` specifies where input files (such as the ``grid`` file, or ``logfile-definitions`` file are located.


Detached mode
-------------

The ``docker run`` examples above use the ``-it`` option.
To run multiple instances of COMPAS, an alternative is to use detached mode (`-d`)\ [#f6]_. In detached mode, containers are run in 
the background of the current shell - they do not receive input or display output.

Typing::

    docker run --rm -d -v $(pwd)/compas-logs/run_0:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_01.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    docker run --rm -d -v $(pwd)/compas-logs/run_1:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_02.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    docker run --rm -d -v $(pwd)/compas-logs/run_2:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_03.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &
    docker run --rm -d -v $(pwd)/compas-logs/run_3:/app/COMPAS/logs -v $(pwd)/runSubmitMMsolar_04.py:/app/starts/runSubmit.py teamcompas/compas python3 /app/starts/runSubmit.py &

runs 4 separate instances of COMPAS, each with its own ``runSubmit.py`` file and logging directory, and all local console output supressed.

To see the console output of detached containers to check progress, first type::

  docker ps

to get the container id of interest, then type::

    docker logs container_id


.. rubric:: Footnotes

.. [#f1] https://docs.docker.com/engine/reference/run/#clean-up---rm
.. [#f2] https://docs.docker.com/engine/reference/run/#foreground
.. [#f3] https://docs.docker.com/storage/bind-mounts/
.. [#f4] https://docs.docker.com/engine/reference/run/#env-environment-variables
.. [#f5] https://docs.docker.com/engine/reference/run/
.. [#f6] https://docs.docker.com/engine/reference/run/#detached--d

   
Running COMPAS via Docker
=========================

Docker functionality has been added to COMPAS to reduce time and effort required to set up the COMPAS deployment environment.

Instead of having to install and configure several libraries and tools (e.g. ``python``/``pip``, ``numpy``, ``g++``, ``boost``, ``hdf5``) 
which can vary considerably beween operating systems and existing toolchains, users can instead opt to install ``Docker`` and run COMPAS 
with a single command.

This also gives users the ability to run COMPAS on cloud solutions like ``AWS EC2``\ [#f1]_ or ``Google Compute Engine``\ [#f2]_ where hundreds 
of cores can be provisioned without having to manually configure the environment.

``Docker`` works by creating an isolated and standalone environment known as a `container`\ [#f3]_. Containers can be created or destroyed 
without affecting the host machine or other containers (containers can still interact with each other and the host machine through mounted 
directories/files or exposed ports).

Containers are instances of images. An image is a pre-defined setup/environment that is instantiated when started as a container (containers 
are to images what objects are to classes in the OO paradigm)\ [#f4]_. 

Containers are (almost) always run as a ``Linux`` environment. A major benefit of this is the ability to run Linux applications in a ``Windows`` 
or ``macOS`` environment without having to jump through hoops or have a diminished experience.

Image definitions can be defined by users (e.g. ``Dockerfiles``); there are also standard images publicly available on ``dockerHub``\ [#f5]_

This following sections assume ``Docker`` has been installed and is running. For Windows and MacOS users, see 
`Docker Desktop <https://www.docker.com/products/docker-desktop>`_.


.. toctree::
   :maxdepth: 1

   running-via-docker-installation
   running-via-docker-running


.. rubric:: Footnotes

.. [#f1] https://aws.amazon.com/ec2/
.. [#f2] https://cloud.google.com/compute
.. [#f3] https://www.docker.com/resources/what-container
.. [#f4] `https://stackoverflow.com/questions/23735149 <https://stackoverflow.com/questions/23735149/what-is-the-difference-between-a-docker-image-and-a-container#:~:text=An%20instance%20of%20an%20image,of%20layers%20as%20you%20describe.&text=You%20can%20see%20all%20your,an%20image%20is%20a%20container>`_
.. [#f5] https://hub.docker.com/

Floating-point comparisons
==========================

Floating-point comparisons are inherently problematic. Testing floating-point numbers for equality, or even inequality, is fraught with
problems due to the internal representation of floating-point numbers: floatingpoint numbers are stored with a fixed number of binary 
digits, which limits their precision and accuracy. The problems with floating-point comparisons are even more evident if one or both of 
the numbers being compared are the results of (perhaps several) floating-point operations (rather than comparing constants).

To avoid the problems associated with floating-point comparisons it is (almost always) better to do any such comparisons with a tolerance
rather than an absolute comparison. To this end, a floating-point comparison function has been provided, and (almost all of) the 
floating-point comparisons in the code have been changed to use that function. The function uses both an absolute tolerance and a relative 
tolerance, which are both declared in constants.h. Whether the function uses a tolerance or not can be changed by ``#define``-ing or 
``#undef``-ing the ``COMPARE_WITH_TOLERANCE`` flag in ``constants.h`` (so the change is a compile-time change, not run-time).


The compare function is defined in utils.h and is implemented as follows::

    static int Compare(const double p_X, const double p_Y) { |br|
    #ifdef COMPARE WITH TOLERANCE
    
        return (fabs(p X – p Y) <= max( FLOAT_TOLERANCE_ABSOLUTE,
                                        FLOAT_TOLERANCE_RELATIVE * 
                                        max( fabs(p_X), fabs(p Y)))) ? 0 
                                                                     : (p_X < p_Y ? –1 : 1);
    #else

        return (p_X == p_Y) ? 0 : (p_X < p_Y ? –1 : 1);

    #endif



If ``COMPARE_WITH_TOLERANCE`` is defined, ``p_X`` and ``p_Y`` are compared with tolerance values, whereas if ``COMPARE_WITH_TOLERANCE`` is
not defined the comparison is an absolute comparison.

The function returns an integer indicating the result of the comparison:
    .. list-table::
       :widths: 8 92 
       :header-rows: 0
       :class: aligned-text

       * - –1 
         - indicates that ``p_X`` is considered to be less than ``p_Y``
       * - |_| |_| 0
         - indicates ``p_X`` and ``p_Y`` are considered to be equal
       * - +1
         - indicates that ``p_X`` is considered to be greater than ``p_Y``

The comparison is done using both an absolute tolerance and a relative tolerance. The tolerances can be defined to be the same number, or
different numbers. If the relative tolerance is defined as 0.0, the comparison is done using the absolute tolerance only, and if the 
absolute tolerance is defined as 0.0 the comparison is done with the relative tolerance only.

Absolute tolerances are generally more effective when the numbers being compared are small – so using an absolute tolerance of (say) 
0.0000005 is generally effective when comparing single-digit numbers (or so), but is less effective when comparing numbers in the thousands
or millions. For comparisons of larger numbers a relative tolerance is generally more effective (the actual tolerance is wider because the 
relative tolerance is multiplied by the larger absolute value of the numbers being compared).

The tolerances used for the comparison are defined in ``constants.h`` as ``FLOAT_TOLERANCE_ABSOLUTE`` and ``FLOAT_TOLERANCE_RELATIVE``.

There is a little overhead in the comparisons even when the tolerance comparison is disabled, but it shouldn’t be prohibitive.
Constants source file
=====================

``constants.h`` is the COMPAS ``C++`` constants source file.

As well as plain constant values, many distribution and prescription identifiers are declared in ``constants.h``. These are mostly 
declared as enum classes, with each enum class having a corresponding map of labels. The benefit is that the values of a particular
(e.g.) prescription are limited to the values declared in the enum class, rather than any integer value, so the compiler will complain 
if an incorrect value is inadvertently used to reference that prescription.

For example, the Common_Envelope Accretion Prescriptions are declared in ``constants.h`` thus::

    enum class CE_ACCRETION_PRESCRIPTION: int { ZERO, CONSTANT, UNIFORM, MACLEOD };

    const COMPASUnorderedMap<CE_ACCRETION_PRESCRIPTION, std::string> CE_ACCRETION_PRESCRIPTION_LABEL = {
        { CE_ACCRETION_PRESCRIPTION::ZERO, ”ZERO” },
        { CE_ACCRETION_PRESCRIPTION::CONSTANT, ”CONSTANT” },
        { CE_ACCRETION_PRESCRIPTION::UNIFORM, ”UNIFORM” },
        { CE_ACCRETION_PRESCRIPTION::MACLEOD, ”MACLEOD” },
    };

Refer to ``constants.h`` for the definition of ``COMPASUnorderedMap``.

Note that the values allowed for variables of type ``CE_ACCRETION_PRESCRIPTION`` are limited to ZERO, CONSTANT, UNIFORM, and MACLEOD – 
anything else will cause a compilation error.

The unordered map ``CE_ACCRETION_PRESCRIPTION_LABEL`` declares a string label for each ``CE_ACCRETION_PRESCRIPTION``, and is indexed by 
``CE_ACCRETION_PRESCRIPTION``. The strings declared in ``CE_ACCRETION_PRESCRIPTION_LABEL`` are used by the Options service to match user
input to the required ``CE_ACCRETION_PRESCRIPTION``. These strings can also be used if an English description of the value of a variable
is required: instead of just printing an integer value that maps to a ``CE_ACCRETION_PRESCRIPTION``, the string label associated with the
prescription can be printed.


Stellar types are also declared in ``constants.h`` via an enum class and associated label map. This allows stellar types to be referenced
using symbolic names rather than an ordinal number. The stellar types enum class is ``STELLAR_TYPE``, and is declared as::

    enum class STELLAR_TYPE: int {
        MS_LTE_07,
        MS_GT_07,
        HERTZSPRUNG_GAP,
        FIRST_GIANT_BRANCH,
        CORE_HELIUM_BURNING,
        EARLY_ASYMPTOTIC_GIANT_BRANCH,
        THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH,
        NAKED_HELIUM_STAR_MS,
        NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,
        NAKED_HELIUM_STAR_GIANT_BRANCH,
        HELIUM_WHITE_DWARF,
        CARBON_OXYGEN_WHITE_DWARF,
        OXYGEN_NEON_WHITE_DWARF,
        NEUTRON_STAR,
        BLACK_HOLE,
        MASSLESS_REMNANT,
        CHEMICALLY_HOMOGENEOUS,
        STAR,
        BINARY_STAR,
        NONE
    };

Ordinal numbers can still be used to reference the stellar types, and because of the order of definition in the enum class the ordinal numbers
match those given in :cite:`Hurley2000`.

The label map ``STELLAR_TYPE_LABEL`` can be used to print text descriptions of the stellar types, and is declared as::

    const std::unordered map<STELLAR_TYPE, std::string> STELLAR_TYPE_LABEL = {
        { STELLAR TYPE::MS_LTE_07,                                 ”Main Sequence <= 0.7” },
        { STELLAR_TYPE::MS_GT_07,                                  ”Main Sequence > 0.7” },
        { STELLAR_TYPE::HERTZSPRUNG_GAP,                           ”Hertzsprung Gap” },
        { STELLAR_TYPE::FIRST_GIANT_BRANCH,                        ”First Giant Branch” },
        { STELLAR_TYPE::CORE_HELIUM_BURNING,                       ”Core Helium Burning” },
        { STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,             ”Early Asymptotic Giant Branch” },
        { STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH, ”Thermally Pulsing Asymptotic Giant Branch” },
        { STELLAR_TYPE::NAKED_HELIUM_STAR_MS,                      ”Naked Helium Star MS” },
        { STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,         ”Naked Helium Star Hertzsprung Gap” },
        { STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,            ”Naked Helium Star Giant Branch” },
        { STELLAR_TYPE::HELIUM_WHITE_DWARF,                        ”Helium White Dwarf” },
        { STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,                 ”Carbon-Oxygen White Dwarf” },
        { STELLAR_TYPE::OXYGEN NEON WHITE DWARF,                   ”Oxygen-Neon White Dwarf” },
        { STELLAR_TYPE::NEUTRON_STAR,                              ”Neutron Star” },
        { STELLAR_TYPE::BLACK_HOLE,                                ”Black Hole” },
        { STELLAR_TYPE::MASSLESS_REMNANT,                          ”Massless Remnant” },
        { STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS,                    ”Chemically Homogeneous” },
        { STELLAR_TYPE::STAR,                                      ”Star” },
        { STELLAR_TYPE::BINARY_STAR,                               ”Binary Star” },
        { STELLAR_TYPE::NONE,                                      ”Not a Star!” }
    };

To print the ordinal number of a stellar type as an integer (sometimes referred to as the "Hurley type"), use ``static_cast``.  For example::

    std::cout << "CHeB ordinal number = " << static_cast<int>(STELLAR_TYPE::CORE_HELIUM_BURNING) << "\n";

To print the text label of a stellar type::

    std::cout << "CHeB label = " << STELLAR_TYPE_LABEL.at(STELLAR_TYPE::CORE_HELIUM_BURNING) << "\n";
Object identifiers
==================

All objects (instantiations of a class) are assigned unique object identifiers of type ``OBJECT_ID`` (unsigned long int - see 
``constants.h`` for the typedef). The purpose of the unique object id is to aid in object tracking and debugging.

**Note that the object id is not the same as, nor it superseded by, the RANDOM_SEED value assigned to each single or binary star**.

The RANDOM_SEED is used to seed the random number generator, and can be used to uniquely identify a single or binary star. The object
id is more granular the the RANDOM_SEED. For example, each binary star is comprised of multiple objects: the ``BinaryStar`` object, 
which contains two ``BaseBinaryStar`` objects (the object undergoing evolution, and a saved copy); each ``BaseBinaryStar`` object 
contains two ``BinaryConstituentStar`` objects (one for each of the constituent stars), and each ``BinaryConstituentStar`` object 
inherits from the ``Star`` class, which contains two ``BaseStar`` objects (the underlying star and a saved copy). Whereas the RANDOM_SEED
uniquely identifies (for example) a binary star, and so identifies the collection of objects that comprise the binary star, the object ids
uniquely identify the constituent objects of the binary star.  Object tracking at this lower level cannot be achieved using the RANDOM_SEED,
hence the need for object ids when debugging.

As well as unique object ids, all objects are assigned an object type (of type ``OBJECT_TYPE`` – see ``constants.h`` for the enum class
declaring ``OBJECT_TYPE``), and a stellar type where applicable (of type ``STELLAR_TYPE`` – see ``constants.h`` for the enum class declaring
``STELLAR_TYPE``).

Objects should expose the following functions::

    OBJECT_ID ObjectId() const { return m ObjectId; }
    OBJECT_TYPE ObjectType() const { return m ObjectType; }
    STELLAR_TYPE StellarType() const { return m StellarType; }

If any of the functions are not applicable to the object, then they must return "\*::NONE" (all objects should implement ``ObjectId()``
correctly).

Any object that uses the Errors service (i.e. the ``ERRORS`` and ``WARNINGS`` macros) must expose these functions: the functions are 
called by the ``ERRORS`` and ``WARNINGS`` macros (:doc:`./Services/services-error-handling`).Developer guide
===============

TeamCOMPAS welcomes the active involvement of colleagues and others interested in the ongoing development and improvement of
the COMPAS software. We hope this Developer Guide helps anyone interested in contributing to the COMPAS software. We expect
this guide to be a living document and improve along with the improvements made to the software.

This guide covers the C++ COMPAS population synthesis code.  Post-processing code (Python scrpts) are not (yet) covered
by this guide.


.. toctree::
   :maxdepth: 1

   ./SSE/SSE
   ./BSE/BSE
   ./object-ids
   ./Services/services
   ./floating-point-comparisons
   ./constants-dot-h
   ./changelog
   ./Programming style/programming-style
   ./Developer build/developer-building-compas
   
Change log
==========

The COMPAS code includes a change log header file - ``changelog.h``.  The change log records changes to the COMPAS ``C++`` code, as well
as the COMPAS version string.

Following is a fragment of the change log::

    // 02.09.06      JR - Apr 07, 2020 - Defect repair:
    //                                      - corrected calculation in return statement for Rand::Random(const double p_Lower, const double p_Upper) (issue #201)
    //                                      - corrected calculation in return statement for Rand::RandomInt(const double p_Lower, const double p_Upper) (issue #201)
    // 02.09.07      SS - Apr 07, 2020 - Change eccentricity, semi major axis and orbital velocity pre-2nd supernove to just pre-supernova everywhere in the code
    // 02.09.08      SS - Apr 07, 2020 - Update zetaMainSequence=2.0 and zetaHertzsprungGap=6.5 in Options::SetToFiducialValues
    // 02.09.09      JR - Apr 11, 2020 - Defect repair:
    //                                      - restored property names in COMPASUnorderedMap<STAR_PROPERTY, std::string> STAR_PROPERTY_LABEL in constants.h (issue #218) (was causing logfile definitions files to be parsed incorrectly)
    // 02.09.10	     IM - Apr 12, 2020 - Minor enhancement: added Mueller & Mandel 2020 remnant mass and kick prescription, MULLERMANDEL
    //  			                     Defect repair: corrected spelling of output help string for MULLER2016 and MULLER2016MAXWELLIAN
    // 02.10.01	     IM - Apr 14, 2020 - Minor enhancement: 
    //  				                            - moved code so that SSE will also sample SN kicks, following same code branch as BSE 
    // 02.10.02      SS - Apr 16, 2020 - Bug Fix for issue #105 ; core and envelope masses for HeHG and TPAGB stars
    // 02.10.03      JR - Apr 17, 2020 - Defect repair:
    //                                      - added LBV and WR winds to SSE (issue #223)
    // 02.10.04	     IM - Apr 25, 2020 - Minor enhancement: moved Mueller & Mandel prescription constants to constants.h, other cleaning of this option

Recorded in each entry of the  changelog is:

    - the new COMPAS version number
    - initials of the developer making the changes
    - date of the changes
    - a brief description of the changes - including github issue number where appropriate


COMPAS version number
---------------------

Currently the COMPAS version number is set manually whenever changes are made to the code.  A planned enhancement is to have the v ersion number
increment automaticall whenever a github pull request is mmerged.

The version number is formatted as: 

    `major`\.\ `minor`\ .\ `defect`

where each of the components is a 2-digit integer.

The `major` number will increment whenever significant new functionality is added to COMPAS. |br|
The `minor` number will increment whenever minor new functionality is added to COMPAS. |br|
The `defect` number will increment whenever defect repairs are made to COMPAS.

`major` and `minor` (and `significant`) are somewhat subjective terms. Some guidance on what constitures a `major` release vs a `minor` release is
given below. `defect` releases are releases that repair known defects.

The COMPAS version string is recorded towards the end of the change log file, and should be incremented whenver a change is made to the code::

    const std::string VERSION_STRING = "02.22.00";


Major releases
--------------

A major release typically includes substantial changes to the way the application functions, and introduces key improvements in functionality. A major
release might include:

    - refactoring the code base `(as we did when we moved from COMPAS v1 to COMPAS v2)`
    - significantly improved performance
    - significant changes to the user interface `(the new grid file format probably should have been a major release)`
    - significant new features
    - removing deprecated features
    - integration with other applications
    
Major releases typically occur somewhat infrequently.


Minor releases
--------------

Minor releases introduce new features to the application. Minor releases are smaller than major 
releases - they can be regarded as `edits` to the current version of the application. Minor releases are not a total overhaul - they enhance and 
improve existing functionality. A minor release might include:

    - limited new features and functionality
    - small updates to existing features

Minor releases typically occur fairly frequently.
SSE evolution model
===================

The high-level stellar evolution model is shown in :ref:`Figure 2 <fig-2>`.

.. _fig-2:

.. figure:: ../../../images/SSE-flow-chart-compressed.svg
    :width: 530px
    :height: 800px
    :align: center
    :figclass: align-center
    :alt: SSE flow chart

    Figure 2 High-level SSE evolution.

The stellar evolution model is driven by the ``Evolve()`` function in the ``Star`` class, which evolves the star through its entire lifetime
by doing the following::

    DO:
        1. calculate time step

            a) calculate the giant branch parameters (as necessary)
            b) calculate the timescales|
            c) choose time step

        2. save the state of the underlying star object

        3. DO:

               a) evolve a single time step
               b) if too much change

                    i) revert to the saved state
                   ii) reduce the size of the time step

           UNTIL timestep not reduced

        4. resolve any mass loss

           a) update initial mass (mass0)
           b) update age after mass loss
           c) apply mass transfer rejuvenation factor

        5. evolve to the next stellar type if necessary

    WHILE the underlying star object is not one of: { HeWD, COWD, ONeWD, NS, BH, MR }

Evolving the star through a single time step (step 3a above) is driven by the ``UpdateAttributesAndAgeOneTimestep()`` function in the
``BaseStar`` class which does the following::

    1. check if the star should be a massless remnant
    2. check if the star is a supernova
    3. if evolution on the phase should be performed
         a) evolve the star on the phase – update stellar attributes
         b) check if the star should evolve off the current phase to a different stellar type
       else
         c) ready the star for the next time step

Evolving the star on its current phase, and off the current phase and preparing to evolve to a different stellar type, is handled by
two functions in the ``BaseStar`` class: ``EvolveOnPhase()`` and ``ResolveEndOfPhase()``.

The ``EvolveOnPhase()`` function does the following::

     1. Calculate Tau
     2. Calculate CO Core Mass
     3. Calculate Core Mass
     4. Calculate He Core Mass
     5. Calculate Luminosity
     6. Calculate Radius
     7. Calculate Perturbation Mu
     8. Perturb Luminosity and Radius
     9. Calculate Temperature
    10. Resolve possible envelope loss

Each of the calculations in the ``EvolveOnPhase()`` function is performed in the context of the star evolving on its current phase.
Each of the classes implements their own version of the calculations (via member functions) – some may inherit functions from the
inheritance chain, while others might just return the value unchanged if the calculation is not relevant to their stellar type.


The ``ResolveEndOfPhase()`` function does the following::

     1. Resolve possible envelope loss
     2. Calculate Tau
     3. Calculate CO Core Mass
     4. Calculate Core Mass
     5. Calculate He Core Mass
     6. Calculate Luminosity
     7. Calculate Radius
     8. Calculate Perturbation Mu
     9. Perturb Luminosity and Radius
    10. Calculate Temperature
    11. Evolve star to next phase

Each of the calculations in the ``ResolveEndOfPhase()`` function is performed in the context of the star evolving off its current 
phase to the next phase. The remainder of the code (in general terms) supports these main driver functions.
Single star evolution
=====================

The Single Star Evolution (``SSE``) architecture is based on the classification of individual stars, with each stellar classification
being described by a separate ``C++`` class.


.. toctree::
   :maxdepth: 1

   SSE-class-hierarchy
   SSE-evolution-model

SSE class hierarchy
===================

:ref:`Figure 1 <fig-1>` shows the ``SSE`` class diagram, where the arrows indicate inheritance. The COMPAS ``C++`` code is implemented
using multiple inheritance, and all stellar classes also inherit directly from the ``BaseStar`` class (arrows not shown in 
:ref:`Figure 1 <fig-1>` for clarity). Each of the stellar classes encapsulates data structures and algorithms specific to the evolutionary
phase corresponding to the class.

The main class for ``SSE`` is the ``Star`` class. The ``Star`` class is a wrapper that abstracts away the details of the star and the 
evolution. Internally the ``Star`` class maintains a pointer to an object representing the star being evolved, with that object being 
an instance of one of the following classes:

    .. list-table::
       :widths: 25 75 
       :header-rows: 0
       :class: aligned-text

       * - ``MS_lte_07``
         - Main Sequence :math:`\small \leq 0.07 M_\odot`
       * - ``MS_gt_07``
         - Main Sequence :math:`\small \gt 0.07 M_\odot`
       * - ``CH``
         - Chemically Homogeneous
       * - ``HG``
         - Hertzsprung Gap
       * - ``FGB``
         - First Giant Branch
       * - ``CHeB``
         - Core Helium Burning
       * - ``EAGB``
         - Early Asymptotic Giant Branch
       * - ``TPAGB``
         - Thermally Pulsing Asymptotic Giant Branch
       * - ``HeMS``
         - Helium Main Sequence
       * - ``HeHG``
         - Helium Hertzsprung Gap
       * - ``HeGB``
         - Helium Giant Branch
       * - ``HeWD``
         - Helium White Dwarf
       * - ``COWD``
         - Carbon-Oxygen White Dwarf
       * - ``ONeWD``
         - Oxygen-Neon White Dwarf
       * - ``NS``
         - Neutron Star
       * - ``BH``
         - Black Hole
       * - ``MR``
         - Massless Remnant

which track the phases from :cite:`Hurley2000`, with the exception of the ``CH`` class for Chemically Homogeneous stars,
which are not described in :cite:`Hurley2000`.

Several other ``SSE`` classes are defined:

    ``BaseStar`` |br|
    ``MainSequence`` |br|
    ``GiantBranch`` |br|
    ``Remnants`` |br|
    ``WhiteDwarfs`` |br|

These extra classes are included to allow inheritance of common functionality.

The ``BaseStar`` class is the main class for the underlying star object held by the ``Star`` class. The ``BaseStar`` class defines all member
variables, and many member functions that provide common functionality. Similarly, the ``MainSequence`` and ``GiantBranch`` classes provide 
repositories for common functionality for main sequence and giant branch stars respectively, and the the ``Remnants`` and ``WhiteDwarfs`` classes
provide repositories for common functionality for remnant and white dwarf stars respectively.

.. _fig-1:

.. figure:: ../../../images/SSE-class-diagram-compressed.svg
    :width: 650px
    :height: 455px
    :align: center
    :figclass: align-center
    :alt: SSE class diagram

    Figure 1 SSE class & container diagram.

``CH`` (Chemically Homogeneous) class stars inherit from the ``MS_gt_07`` class because (in this implementation) they are just (large) main
sequence stars that have a static radius.

``HG`` (Hertzsprung Gap) class stars inherit from the ``GiantBranch`` class because they share the giant branch parameters described in 
:cite:`Hurley2000`, section 5.2.

Each class has its own set of member functions that calculate various attributes of the star according to the phase the class represents (using
the equations and parameters from :cite:`Hurley2000` where applicable).
Programming style & conventions
===============================

The goal of coding to a suggested style is readability and maintainability – if many developers implement code in COMPAS with their own coding
style, readability and maintainability will be more difficult than if a consistent style is used throughout the code. Strict adherence isn’t 
really necessary, but it will make it easier on all COMPAS developers if the coding style is consistent throughout.

Some elements of programming style:


.. toctree::
   :maxdepth: 1

   programming-style-style-comments
   programming-style-style-braces
   programming-style-style-indentation
   programming-style-style-func-parms
   programming-style-style-performance
   programming-style-style-naming
Polymorphism
============

Polymorphism means having many forms. In OOP, polymorphism occurs when there is a hierarchy of classes and they are related by inheritance.

Following the discussion earlier regarding inheritance, in the OOP paradigm, and ``C++`` specifically, derived classes can override methods 
defined by ancestor classes, allowing a derived class to implement functions specific to its circumstances. This means that a call to a class 
member function will cause a different function to be executed depending on the type of object that invokes the function. Descendent classes 
of a class that has overridden a base class member function inherit the overridden function (but can override it themselves).

COMPAS makes heavy use of inheritance and polymorphism, especially for the implementation of the different stellar types.Comments
========

An old, but good, rule-of-thumb is that any file that contains computer code should be about one-third code, one-third comments, and one-third
white space. Adhering to this rule-of-thumb just makes the code a bit easier on the eye, and provides some description (at least of the intention)
of the implementation.
Programming style
=================

Everyone has their own preferences and style, and the nature of a project such as COMPAS will reflect that. However, there is a need
to suggest some guidelines for programming style, naming conventions etc. Following is a description of some of the elements of 
programming style and naming conventions used to develop COMPAS v2. These may evolve over time.


.. toctree::
   :maxdepth: 1

   programming-style-oo
   programming-style-style
Performance & optimisation
==========================

In general, COMPAS developers should code for performance – within reason. Bear in mind that many functions will be called many, many thousands of 
times (in some cases, millions) in one execution of the program.

- Avoid calculating values inside loops that could be calculated once outside the loop.
- Try to use constants where possible.
- Use multiplication in preference to functions such as pow() and sqrt() (note that pow() is very expensive computationally; sqrt() is expensive, but much less expensive than pow()).
- Don’t optimise to the point that readability and maintainability is compromised. Bear in mind that most compilers are good at optimising, and are very forgiving of less-than-optimally-written code (though they are not miracle workers...).
Object-oriented programming
===========================

COMPAS is written in ``C++``, an object-oriented programming (OOP) language, and OOP concepts and conventions should apply throughout
the code. There are many texts and web pages devoted to understanding ``C++`` and OOP – following is a brief description of the key OOP
concepts:


.. toctree::
   :maxdepth: 1

   programming-style-oo-abstraction
   programming-style-oo-encapsulation
   programming-style-oo-inheritance
   programming-style-oo-polymorphism
Braces
======

The placement of braces in ``C++`` code (actually, any code that uses braces to enclose scope) is a contentious issue, with many developers having 
long-held, often dogmatic preferences. COMPAS (so far) uses the K&R style (”the one true brace style”) - the style used in the original Unix kernel
and Kernighan and Ritchie’s book :doc:`The C Programming Language <../../references>`.

The K&R style puts the opening brace on the same line as the control statement:

::

    while (x == y) {
        call_something();
        var1 = var2
        call_somethingelse();
    }

Note also the space between the keyword while and the opening parenthesis, surrounding the ``==`` operator, and between the closing parenthesis 
and the opening brace. Spaces in those places help with code readability. Surrounding all arithmetic operators with spaces is preferred.
Abstraction
===========

For any entity, product, or service, the goal of abstraction is to handle the complexity of the implementation by hiding details that 
don’t need to be known in order to use, or consume, the entity, product, or service. In the OOP paradigm, hiding details in this way 
enables the consumer to implement more complex logic on top of the provided abstraction without needing to understand the hidden 
implementation details and complexity. (There is no suggestion that consumers shouldn’t understand the implementation details, but they
shouldn’t need to in order to consume the entity, product, or service).

Abstraction in ``C++`` is achieved via the use of objects – an object is an instance of a class, and typically corresponds to a 
real-world object or entity (in COMPAS, usually a star or binary star). An object maintains the state of an object (via class member 
variables), and provides all necessary means of changing the state of the object (by exposing public class member functions (methods)). 
A class may expose public functions to allow consumers to determine the value of class member variables (“getters”), and to set the value
of class member variables (“setters”).
Encapsulation
=============

Encapsulation binds together the data and functions that manipulate the data in an attempt to keep both safe from outside interference and
accidental misuse. An encapsulation paradigm that does not allow calling code to access internal object data and permits access through 
functions only is a strong form of abstraction. ``C++`` allows developers to enforce access restrictions explicitly by defining class member
variables and functions as private, protected, or public. These keywords are used throughout COMPAS to enforce encapsulation.

There are very few circumstances in which a consumer should change the value of a class member variable directly (via the use of a setter 
function) – almost always consumers should present new situational information to an object (via a public member function), and allow the 
object to respond to the new information. For example, in COMPAS, there should be almost no reason for a consumer of a star object to 
directly change (say) the radius of the star – the consumer should inform the star object of new circumstances or events, and allow the star 
object to respond to those events (perhaps changing the value of the radius of the star). Changing a single class member variable directly 
introduces the possibility that related class member variables (e.g. other attributes of stars) will not be changed accordingly. Moreover, 
developers changing the code in the future should, in almost all cases, expect that the state of an object is maintained consistently by the
object, and that there should be no unexpected side-effects caused by calling non class-member functions.

In short, changing the state of an object outside the object is potentially unsafe and should be avoided where possible.
Naming conventions
------------------

COMPAS (so far) uses the following naming conventions:

- All variable names should be in camelCase – don’t use underscore to separate words.
- Function names should be in camelCase, beginning with an uppercase letter. Function names should be descriptive.
- Class member variable names are prefixed with 'm\_', and the character immediately following the prefix should be uppercase (in most cases – sometimes, for well-known names or words that are always written in lowercase, lowercase might be used).
- Local variable names are just camelCase, beginning with a lowercase letter (again, with the caveat that sometimes, for well-known names or words that are always written in uppercase, uppercase might be used).
- Function parameter names are prefixed with 'p\_', and the character immediately following the prefix should be uppercase (again, with the caveat that sometimes, for well-known names or words that are always written in lowercase, lowercase might be used).
Function parameters
===================

In most cases, function parameters should be input only – meaning that the values of function parameters should not be changed by the function. 
Anything that needs to be changed and returned to the caller should be returned as a functional return. There are a few exceptions to this in COMPAS – 
all were done for performance reasons, and are documented in the code.

To avoid unexpected side-effects, developers should expect (in most cases) that any variables they pass to a function will remain unchanged – all 
changes should be returned as a functional return.
Inheritance
===========

Inheritance allows classes to be arranged in a hierarchy that represents is-a-type-of relationships. All nonprivate class member variables and
functions of the parent (base) class are available to the child (derived) class (and, therefore, child classes of the child class). This allows
easy re-use of the same procedures and data definitions, in addition to describing real-world relationships in an intuitive way. ``C++`` allows
multiple inheritance – a class may inherit from multiple parent classes.

Derived classes can define additional class member variables (using the private, protected, and public access restrictions), which will be 
available to any descendent classes (subject to inheritance rules), but will only be available to ancestor classes via the normal access methods
(getters and setters).
Indentation
===========

There is ongoing debate in the programming community as to whether indentation should be achieved using spaces or tabs (strange, but true...). 
The use of spaces is more common. Unfortunately a mix of spaces and tabs doesn’t work well with some editors - we should settle on one method and 
try to stick to it.  COMPAS (so far) has a mix of both – no matter which convention we choose, once there are instances of both in the code it's 
hard to recover.  So while using tabs is probably preferred, use whatever is convenient (pragmatism is your friend...). 

COMPAS (mostly) uses an indentation size of 4 spaces - again we should settle on a size and stick to it.
Binary star evolution
=====================

The Binary Star Evolution (``BSE``) architecture implements is based on the classification of individual stars, with each stellar classification
being described by a separate ``C++`` class.


.. toctree::
   :maxdepth: 1

   BSE-class-hierarchy
   BSE-evolution-model

BSE evolution model
-------------------

The high-level binary evolution model is shown in :ref:`Figure 4 <fig-4>`.

.. _fig-4:

.. figure:: ../../../images/BSE-flow-chart-compressed.svg
    :width: 650px
    :height: 1574px
    :align: center
    :figclass: align-center
    :alt: BSE flow chart

    Figure 4 High-level BSE evolution.


The binary evolution model is driven by the ``Evolve()`` function in the ``BaseBinaryStar`` class, which evolves the star through its entire 
lifetime by doing the following::

    if touching
        STOP = true
    else
        calculate initial time step
        STOP = false
    
    DO WHILE NOT STOP AND NOT max iterations:
    
        evolve a single time step
            evolve each constituent star a single time step (see the :ref:`SSE evolution model`)
        
        if error OR unbound OR touching OR Massless Remnant
            STOP = true
        else
            evaluate the binary
                calculate mass transfer
                calculate winds mass loss
    
                if common envelope
                    resolve common envelope
                else if supernova
                    resolve supernova
                else
                    resolve mass changes
    
                evaluate supernovae
                calculate total energy and angular momentum
                update magnetic field and spin: both constituent stars
    
            if unbound OR touching OR merger
                STOP = true
            else
                if NS+BH
                    resolve coalescence
                    STOP = true
                else
                    if WD+WD OR max time
                        STOP = true
                    else
                        if NOT max iterations
                            calculate new time step 
BSE class hierarchy
-------------------

The ``BSE`` class diagram is shown in :ref:`Figure 3 <fig-3>`.

.. _fig-3:

.. figure:: ../../../images/BSE-class-diagram-compressed.svg
    :width: 600px
    :height: 184px
    :align: center
    :figclass: align-center
    :alt: BSE class diagram

    Figure 3 BSE class & container diagram.

The ``BSE`` inheritance hierarchy is as follows:

    ``BinaryStar``
        ``BaseBinaryStar``
        
            ``Star`` → ``BinaryConstituentStar`` (star1) |br|
            ``Star`` → ``BinaryConstituentStar`` (star2)

The main class for binary star evolution is the ``BinaryStar`` class. The ``BinaryStar`` class is a wrapper that abstracts away the details
of the binary star and the evolution. Internally the ``BinaryStar`` class maintains a pointer to an object representing the binary star being
evolved, with that object being an instance of the ``BaseBinaryStar`` class. The ``BinaryStar`` class maintains a pointer to a second object,
representing the saved state of the binary star being evolved, with that object also being an instance of the ``BaseBinaryStar`` class. The 
second object is a copy of the binary star being evolved at some earlier timestep, facilitating reverting the binary star to a previous state.

The ``BaseBinaryStar`` class is a container class for the objects that represent the component stars of a binary system. An instance of the
``BaseBinaryStar`` class is a binary system being evolved by COMPAS, and contains a ``BinaryConstituentStar`` class object for each of the
component stars (i.e. the primary and secondary stars), as well as data structures and algorithms specific to the evolution of a binary system.
The ``BaseBinaryStar`` class also maintains pointers to the ``BinaryConstituentStar`` class objects considered to be the current donor and 
accretor during a mass transfer event, as well as pointers to the ``BinaryConstituentStar`` class objects considered to be the current 
supernova and companion star, should one of the stars undergo a supernova event.

The ``BinaryConstituentStar`` class inherits from the ``Star`` class, so objects instantiated from the ``BinaryConstituentStar`` class inherit
the characteristics of the ``Star`` class, particularly the stellar evolution model. The ``BinaryConstituentStar`` class defines member variables
and functions that pertain specifically to a constituent star of a binary system but that do not (generally) pertain to single stars that are not
part of a binary system (there are some functions that are defined in the ``BaseStar`` class and its derived classes that deal with binary star 
attributes and behaviour – in some cases the stellar attributes that are required to make these calculations reside in the ``BaseStar`` class so
it is easier and cleaner to define the functions there).

An instance of the ``BinaryConstituentStar`` class is a single component star of a binary system being evolved by COMPAS, and inhertis from the
``Star``, so will evolve over time through various SSE classes shown in :ref:`fig-1`. The ``BinaryConstituentStar`` class defines additional data
structures and algorithms (to the data structures and algorithms provided by the ``SSE`` classes) required to support the evolution of a binary 
system component star.

The ``BaseBinaryStar`` class is the main class for the underlying binary star object held by the ``BinaryStar`` class. The ``BaseBinaryStar`` 
class defines all member variables that pertain specifically to a binary star, and many member functions that provide binary-star specific 
functionality. Internally, the ``BaseBinaryStar`` class maintains pointers to the two BinaryConstituentStar class objects that constitute the 
binary star.

Building COMPAS locally
=======================

A makefile is provided to build COMPAS locally (``Mkaefile``).  The Makefile provided defines a number of variables that can be 
specified on the command line when ``make`` is run, including variables that allow the user to specify the compiler, the `include` 
directory (for source header files), the `lib` directory (for shared libraries) for each external library required by COMPAS, and 
the COMPAS executable name:

    **CPP** |br|
    Specifies the compiler to be used.  The default value is 'g++'.

    **GSLINCDIR** |br|
    Specifies the `include` directory for the GNU scientific library, ``gsl``.  The default value is '/include'.

    **GSLLIBDIR** |br|
    Specifies the `lib` directory for the GNU scientific library, ``gsl``.  The default value is '/lib'.

    **BOOSTINCDIR** |br|
    Specifies the `include` directory for the BOOST library.  The default value is '/include'
    
    **BOOSTLIBDIR** |br|
    Specifies the `lib` directory for the BOOST library.  The default value is '/lib'.
    
    **HDF5INCDIR** |br|
    Specifies the `include` directory for the HDF5 library.  The default value is '/usr/include/hdf5/serial'

    **HDF5LIBDIR** |br|
    Specifies the `lib` directory for the HDF5 library.  The default value is '/usr/lib/x86_64-linux-gnu/hdf5/serial'
    
    **EXE** |br|
    Specifies the name of the COMPAS exectuable to be built.  The default is 'COMPAS'


For example, typing::

    make GCC=c++ EXE=mycompas -j$(nproc)
    
will cause the `c++` compiler to be used to create the executable file 'mycompas', using all available cores.


The makefile provided also defines several entry points:

    **clean** |br|
    Instructs ``make`` to remove all existing object files (.o), and the COMPAS executable.  A subsequent ``make`` is then forced to 
    compile all source files and link the resultant object files (and external libraries) into a new executable.

    **static** |br|
    Specifies that functions in the external libraries referenced by COMPAS should be statically linked - that is, they are copied into
    the COMPAS executable.  The default executable name for the *static* entry point is `COMPAS_STATIC`.

    **fast** |br| 
    Adds `-march=native` and `-O3` to the compiler flags.

        - specifying `-march=native` enables all instruction subsets supported by the compiling machine, thus producing an executable
          file that will perform better than it would if the full native instruction set was not available.  Care should be taken when
          using the executable created from this entry point - some of the native instructions may not be available for use on machines
          of different architectures, so the resultant executable file is not necessarily portable.
        - specifying `-O3` causes the compiler to perform many optimisations to produce an executable that is optimised for performance
          (at the expense of compile time).\ [#f1]_

    **staticfast** |br|
    The functionality of the **static** and **fast** entry points combined.



.. rubric:: Footnotes

.. [#f1] https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html

COMPAS build methods
====================

We provide two methods to build the COMPAS executable: a method to build COMPAS locally, and a method to create a ``Docker`` image.

Note that each method has a separate, and different, makefile: if functionality is added, modified, or removed from either of the
provided makefiles (``Makefile`` locally, ``Makefile.docker`` for the docker image), the corresponding changes should also be made
to the other makefile.


.. toctree::
   :maxdepth: 1

   ./COMPAS-local-build
   ./docker-developer

   Creating a Docker image
=======================


GitHub CI/CD
------------
 
``Docker`` images of all versions of COMPAS (``dev`` branch) are available at the 
`TeamCOMPAS dockerHub page <https://hub.docker.com/u/teamcompas>`__. Building of these images is performed automatically by the 
``GitHub`` CI/CD process.

Whenever a push to `TeamCOMPAS/dev <https://github.com/TeamCOMPAS/COMPAS/tree/dev>`__ occurs, a continuous deployment process 
automatically
builds\ [#f1]_ a new image and deploys it to ``dockerHub`` with a `tag`\ [#f2]_ that corresponds to the value of ``VERSION_STRING`` 
in ``changelog.h`` (see :doc:`../changelog` for detailed information regarding ``changelog.h``).

At time of writing, `GitHub Actions`\ [#f3]_ facilitates the above process. While this is convenient (because it's free and well 
supported) it is somewhat slow - the COMPAS ``Docker`` image is available 5 - 10 minutes after pushing/merging. A future improvement 
may be to create a `runner`\ [#f4]_ locally with a high core count that can be used to compile COMPAS quickly.

The Github Actions configuration is in ``/.github/workflows/dockerhub-ci.yml``.

See the `Atlassian CI/CD <https://www.atlassian.com/continuous-delivery/principles/continuous-integration-vs-delivery-vs-deployment>`__ 
documentation for detailed information regarding the ``GitHub`` CI/CD process.


Dockerfile
----------

The Dockerfile\ [#f5]_ defines how the docker image is constructed.

Images are created as a combination of layers. During the build process, each layer is cached and only updated on subsequent builds 
if that layer would change. 

The Dockerfile for COMPAS is made up of 8 layers:

    **FROM ubuntu:18.04**\ [#f6]_ |br|
    Use `Ubuntu 18.04 <https://hub.docker.com/_/ubuntu>`__ as a base (provided by Docker Hub)
    
    **WORKDIR /app/COMPAS**\ [#f7]_ |br|
    Effectively ``cd /app/COMPAS`` within the container.
    
    **RUN apt-get update && apt-get install -y ...**\ [#f8]_ |br|
    Install the required dependencies.

        - `-y` so there's no prompt to install any of the packages.
        - `update` and `install` are in the same layer because now if there are any updates, it will force all of the dependencies
          to be re-installed

    **RUN pip3 install numpy**\ [#f8]_ |br|
    Install numpy.

    **COPY src/ src**\ [#f9]_ |br|
    Copy ``./src/`` directory from the local machine to ``./src`` in the container (remembering that `WORKDIR` changes the cwd).

    **RUN mkdir obj bin logs**\ [#f8]_ |br|
    Create the directories required by COMPAS.

    **ENV COMPAS_ROOT_DIR /app/COMPAS**\ [#f10]_ |br|
    Set the required environment variable(s).

    **RUN cd src && make -f Makefile.docker -j $(nproc)**\ [#f8]_ |br|
    Change to the ``src`` directory; make COMPAS using a specific makefile (see below) and as many cores as possible.

A Dockerfile usually ends with a ``CMD`` directive that specifies what command should run when the container is started\ [#f11]_. 
The COMPAS Dockerfile doesn't have a ``CMD`` directive because some users will want to run the executable directly and some will 
want to use ``runSubmit.py``.


Makefile.docker
---------------

A separate makefile is required for ``Docker`` in this scenario for two reasons:

    #. To separate compiled files from source files
    #. To prevent the usage of the ``gcc`` ``-march=native`` compiler option

``-march=native`` is a very useful optimisation for users who compile and run COMPAS on the same machine, however it causes fatal 
errors when running COMPAS on a machine for which it was not compiled. The ``-march=native`` compiler option selects the CPU for
which code should be generated by determining the processor type of the `compiling machine`. Using ``-march=native`` enables all
instruction subsets supported by the compiling machine, thus producing an executable file that will perform better than it would
if the full native instruction set was not available, but some of those instructions may not be available for use on machines of
different architectures - hence the possible fatal run-time errors if the code is run on machines with different architectures).
See `gcc x86-Options <https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html>`__ for detailed information regarding ``-march``.

The Docker makefile provided (``Makefile.docker``) functions exactly like the local makefile provided (``Makefile``) in all other
respects.  See :doc:`./COMPAS-local-build` for a detailed description of the local makefile functionality.


.. rubric:: Footnotes

.. [#f1] https://docs.docker.com/engine/reference/commandline/build/
.. [#f2] https://docs.docker.com/engine/reference/commandline/tag/
.. [#f3] https://github.com/features/actions
.. [#f4] https://help.github.com/en/actions/getting-started-with-github-actions/core-concepts-for-github-actions#runner
.. [#f5] https://docs.docker.com/engine/reference/builder/
.. [#f6] https://docs.docker.com/engine/reference/builder/#from
.. [#f7] https://docs.docker.com/engine/reference/builder/#workdir
.. [#f8] https://docs.docker.com/engine/reference/builder/#run
.. [#f9] https://docs.docker.com/engine/reference/builder/#copy
.. [#f10] https://docs.docker.com/engine/reference/builder/#env
.. [#f11] https://docs.docker.com/engine/reference/builder/#cmd

Error handling service
======================

An error handling service is provided encapsulated in a singleton object (an instantiation of the ``Errors`` class).

The Errors service provides global error handling functionality. Following is a brief description of the Errors service
(full documentation coming soon...):

Errors are defined in the error catalog in ``constants.h`` (see ``ERROR_CATALOG``). |br|

Errors defined in the error catalog have a scope and message text. The scope is used to determine if and when an error should
be printed.

The current values for scope are:

    .. list-table::
       :widths: 30 70 
       :header-rows: 0
       :class: aligned-text

       * - **NEVER**
         - the error will not be printed.
       * - **ALWAYS**
         -  the error will always be printed.
       * - **FIRST**
         -  the error will be printed only on the first time it is encountered anywhere in the program in the current execution of COMPAS.
       * - **FIRST_IN_OBJECT_TYPE**
         - the error will be printed only on the first time it is encountered anywhere in objects of the same type (e.g. Binary Star objects) in the current execution of COMPAS.
       * - **FIRST_IN_STELLAR_TYPE**
         - the error will be printed only on the first time it is encountered anywhere in objects of the same stellar type (e.g. HeWD Star objects) in the current execution of COMPAS.
       * - **FIRST_IN_OBJECT_ID**
         - the error will be printed only on the first time it is encountered anywhere in an object instance in the current execution of COMPAS (useful for debugging).
       * - **FIRST_IN_FUNCTION**
         - the error will be printed only on the first time it is encountered anywhere in the same function of an object instance in the current execution of COMPAS (i.e. will print more than once if encountered in the same function name in different objects; useful for debugging).

The Errors service provides methods to print both warnings and errors – essentially the same thing, but warning messages are prepended with 
':boldtext:`WARNING:` ', whereas error messages are prepended with ':boldtext:`ERROR:` '.

Errors and warnings are printed by using the macros defined in ``ErrorsMacros.h``. They are:

Error macros
------------

::
    
    SHOW_ERROR(error_number)                         // prints the error message associated with error
                                                     // number (from the error catalog) prepended by
                                                     // 'ERROR: '

    SHOW_ERROR(error_number, error_string)           // prints the error message associated with error
                                                     // number (from the error catalog) prepended by
                                                     // 'ERROR: ', and appends 'error_string'

    SHOW_ERROR_IF(cond, error_number)                // if 'cond' is TRUE, prints the error message
                                                     // associated with error number (from the error
                                                     // catalog) prepended by 'ERROR: '

    SHOW_ERROR_IF(cond, error_number, error_string)  // if 'cond' is TRUE, prints the error message
                                                     // associated with error number (from the error
                                                     // catalog) prepended by 'ERROR: ', and appends 
                                                     // 'error_string'


Warning macros
--------------

The ``WARNING`` macros function in the same way as the ``ERROR`` macros, with the exception that instead of prepending the
message with ':boldtext:`ERROR:` ', the ``WARNING`` macros prepend the message with ':boldtext:`WARNING:` '.

The ``WARNING`` macros are:

::

    SHOW_WARN(error_number)
    SHOW_WARN(error_number, error_string)
    SHOW_WARN_IF(cond, error_number)
    SHOW_WARN_IF(cond, error_number, error_string)

|br|
Error and warning messages always contain:

- The object id of the calling object.
- The object type of the calling object.
- The stellar type of the calling object (will be ”NONE” if the calling object is not a star-type object).
- The function name of the calling function.

Any object that uses the Errors service (i.e. any of the ``ERROR`` and ``WARNING`` macros) must expose the following functions:

::

    OBJECT_ID ObjectId()** const { return m ObjectId; }
    OBJECT_TYPE ObjectType()** const { return m ObjectType; }
    STELLAR_TYPE StellarType()** const { return m StellarType; }

These functions are called by the ``ERROR`` and ``WARNING`` macros. If any of the functions are not applicable to the object, 
then they must return ':italictext:`type`::NONE' (all objects should implement ObjectId() correctly).

The filename to which error records are written when ``Log::Start()`` parameter ``p_ErrorsToFile`` is TRUE is declared in 
``constants.h`` – see the enum class ``LOGFILE`` and associated descriptor map ``LOGFILE_DESCRIPTOR``. Currently the name is 'Error_Log'.
Services
========

A number of services have been provided to help simplify the code. The code for each service is encapsulated in a singleton object
(an instantiation of the relevant class). The singleton design pattern allows the definition of a class that can only be instantiated
once, and that instance effectively exists as a global object available to all the code without having to be passed around as a 
parameter. Singletons are a little anti-OO, but provided they are used judiciously are not necessarily a bad thing, and can be very 
useful in certain circumstances.

Services provided are:


.. toctree::
   :maxdepth: 1

   ./Program options/services-program-options
   ./Random numbers/services-random-numbers
   ./Logging and debugging/services-logging-debugging
   ./services-error-handlingProgram options service
=======================

A Program Options service is provided, encapsulated in a singleton object (an instantiation of the ``Options`` class).

The ``Options`` class member variables are private, and public getter functions have been created for the program options currently
used in the code.

The Options service can be accessed by referring to the ``Options::Instance()`` object. For example, to retrieve the value of 
the ``--quiet`` program option, call the ``Options::Quiet()`` getter function::

    bool quiet = Options::Instance()→Quiet();

Since that could become unwieldy, there is a convenience macro to access the Options service. The macro just defines "OPTIONS" as
"Options::Instance()", so retrieving the value of the ``--quiet`` program option can be written as::

    bool quiet = OPTIONS→Quiet();

The Options service must be initialised before use. Initialise the Options service by calling the ``Options::Initialise()`` function::

    COMMANDLINE_STATUS programStatus = OPTIONS→Initialise(argc, argv);

(see ``constants.h`` for details of the ``COMMANDLINE_STATUS`` type)

See :doc:`../../../User guide/Program options/program-options-list-defaults` for a full list of available program options and their default valaues.
Logging & debugging service
===========================

A logging and debugging service is provided encapsulated in a singleton object (an instantiation of the ``Log`` class).


.. toctree::
   :maxdepth: 1

   ./Base-level logging/services-base-level-logging
   ./services-extended-logging
   ./Logging macros/services-logging-debugging-macros-logging
   ./services-logging-debugging-macros-debugging


:bolditalictext:`Historical note:`

The logging functionality was first implemented when the Single Star Evolution code was refactored, and the base-level of logging
was sufficient for the needs of the ``SSE`` code. Refactoring the Binary Star Evolution code highlighted the need for expanded logging
functionality. To provide for the logging needs of the ``BSE`` code, new functionality was added almost as a wrapper around the original,
base-level logging functionality. Some of the original base-level logging functionality has almost been rendered redundant by the new
functionality implemented for ``BSE`` code, but it remains (almost) in its entirety because it may still be useful in some circumstances.

When the base-level logging functionality was created, debugging functionality was also provided, as well as a set of macros to make
debugging and the issuing of warning messages easier. A set of logging macros was also provided to make logging easier. The debug
macros are still useful, and their use is encouraged (rather than inserting print statements using ``std::cout`` or ``std::cerr``).

When the ``BSE`` code was refactored, some rudimentary error handling functionality was also provided in the form of the 
:doc:`Errors service <../services-error-handling>` an attempt at making error handling easier. Some of the functionality provided by the
:doc:`Errors service <../services-error-handling>` supersedes the ``DBG_WARN*`` macros provided as part of the Log class, but the ``DBG_WARN*``
macros are still useful in some circumstances (and in fact are still used in various places in the code). The ``LOG*`` macros are somewhat less
useful, but remain in case the original base-level logging functionality (that which underlies the expanded logging functionality) is used in 
the future (as mentioned above, it could still be useful in some circumstances).

The expanded logging functionality introduces Standard Log Files - described in :doc:`services-extended-logging`.
Debugging macros
================

A set of macros similar to the :doc:`logging macros <./Logging macros/services-logging-debugging-macros-logging>` is also
provided for debugging purposes. These macros are analogous to their logging counterparts.

The debugging macros write directly to stdout rather than the log file, but their output can also be written to
the log file if desired (see the ``p_DbgToFile`` parameter of ``Log::Start()``, and the ``--debug-to-file`` program
option).

A major difference between the logging macros and the debugging macros is that the debugging macros can be defined
away. The debugging macro definitions are enclosed in an ``#ifdef`` enclosure, and are only present in the source code
if ``#DEBUG`` is defined. This means that if ``#DEBUG`` is not defined (``#undef``), all debugging statements using
the debugging macros will be removed from the source code by the preprocessor before the source is compiled. Un-defining
``#DEBUG`` not only prevents bloat of unused code in the executable, it improves performance. Many of the functions in
the code are called hundreds of thousands, if not millions, of times as the stellar evolution proceeds. Even if the
debugging classes and debugging level are set so that no debug statement is displayed, just checking the debugging level
every time a function is called increases the run time of the program. The suggested use is to enable the debugging
macros (``#define DEBUG``) while developing new code, and disable them (``#undef DEBUG``) to produce a production version
of the executable.

The debugging macros provided are::

    DBG(...)              // analogous to LOG()
    DBG_ID(...)           // analogous to LOG_ID()
    DBG_IF(cond, ...)     // analogous to LOG_IF()
    DBG_ID_IF(cond, ...)  // analogous to LOG_ID_IF()


Two further debugging macros are provided::

    DBG_WAIT(...)
    DBG_WAIT_IF(cond, ...)

The ``DBG_WAIT`` macros function in the same way as their non-wait counterparts (``DBG(...)`` and ``DBG_IF(cond, ...)``),
with the added functionality that they will pause execution of the program and wait for user input before proceeding.

A set of macros for printing warning message is also provided. These are the ``DBG_WARN`` macros::

    DBG_WARN(...)              // analogous to LOG()
    DBG_WARN_ID(...)           // analogous to LOG_ID()
    DBG_WARN_IF(cond, ...)     // analogous to LOG_IF()
    DBG_WARN_ID_IF(cond, ...)  // analogous to LOG_ID_IF()

The ``DBG_WARN`` macros write to stdout via the ``SAY`` macro, so honour the logging classes and level, and are not written
to the debug or errors files.

Note that the ``id`` parameter of the ``LOG`` macros (to specify the logfileId) is not required for the ``DBG`` macros (the
filename to which debug records are written is declared in ``constants.h`` – see the enum class ``LOGFILE`` and associated 
descriptor map ``LOGFILE_DESCRIPTOR``).
Extended logging
================

The extended logging service supports standard log files for both Single Star Evolution ``(SSE``) and Binary Star 
Evolution (``BSE``).

The standard log files defined are:

For ``SSE``:

- SSE_System_Parameters log file
- SSE_Supernovae log file
- SSE_Detailed_Output log file
- SSE_Switchlog log file

For ``BSE``:

- BSE_System_Parameters log file
- BSE_Double_Compact Objects log file
- BSE_Common_Envelopes log file
- BSE_Supernovae log file
- BSE_Pulsar_Evolution log file
- BSE_RLOF_Parameters log file
- BSE_Detailed_Output log file
- BSE_Switchlog log file

The Logging service maintains information about each of the standard log files, and will handle creating, opening, writing and
closing the files. For each execution of the COMPAS program, one (and only one) of each of the log files listed above that
pertain to the mode of evolution (``--mode`` option, ``SSE`` or ``BSE``) will be created, except for the ``Detailed_Output`` 
log files, in which case there will be one log file created for each system (single star or binary star) evolved.

The Logging service provides the following public member functions specifically for managing standard log files:

For ``SSE`` log files::

    BOOL LogSSESystemParameters(CONST T* CONST p_Star, CONST string p_Rec)
    BOOL LogSSESupernovaDetails(CONST T* CONST p_Star, CONST string p_Rec)
    BOOL LogSSEDetailedOutput(CONST T* CONST p_Star, CONST int p_Id, CONST string p_Rec)
    BOOL LogSSESwitchLog(CONST T* CONST p_Star, CONST string p_Rec)

Each ``SSE`` function is passed a pointer to the single star for which details are to be logged (``p_Star``), and a string to be 
written to the log file (``p_Rec``). If ``p_Rec`` is an empty string, the function constructs the log record from the current 
attributes of the star and the default record specifier for the log file (see property vectors in ``constants.h``,  e.g. 
``SSE_DETAILED_OUTPUT_REC``). ``LogSSEDetailedOutput()`` is also passed an integer identifier (typically the loop index of the
star) that is appended to the log file name (``p_Id``).

For ``BSE`` log files::

    BOOL LogBSESystemParameters(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogDoubleCompactObject(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogCommonEnvelope(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogBSESupernovaDetails(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogBSEPulsarEvolutionParameters(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogRLOFParameters(CONST T* CONST p_Binary, CONST string p_Rec)
    BOOL LogBSEDetailedOutput(CONST T* CONST p_Binary, CONST long int p_Id, CONST string p_Rec)
    BOOL LogBSESwitchLog(CONST T* CONST p_Binary, CONST bool p_PrimarySwitching)

Each ``BSE`` function is passed a pointer to the binary star for which details are to be logged (``p_Binary``), and a string to 
be  written to the log file (``p_Rec``). If ``p_Rec`` is an empty string, the function constructs the log record from the current 
attributes of the binary and the default record specifier for the log file (see property vectors in ``constants.h``,  e.g. 
``BSE_DETAILED_OUTPUT_REC``). ``LogBSEDetailedOutput()`` is also passed an integer identifier (typically the loop index of the
binary) that is appended to the log file name (``p_Id``).

Each of the functions listed above will, if necessary, create and open the appropriate log file. Internally the Log service opens
(creates first if necessary) once at first use, and keeps the files open for the life of the program.

The Log service provides a further two functions to manage standard log files::

    BOOL CloseStandardFile()
    BOOL CloseAllStandardFiles()

``CloseStandardFile()`` flushes and closes a standard log file. The function returns a boolean indicating whether the log file was
closed successfully.

``CloseAllStandardFiles()`` flushes and closes all currently open standard log files. The function returns a boolean indicating
whether all standard log files were closed successfully.

Standard log file names are supplied via program options (e.g. ``--logfile-system-parameters``), with default values declared in
``constants.h``.

The extended logging service always sets the log record ``class`` to the name of the standard log file being written to, and the
log record ``level`` to 0.  See :doc:`./Base-level logging/services-base-level-logging` for details regarding log record ``class``
and ``level``.

LOG() macro
===========

::

    LOG(id, ...)

Writes a log record to the log file specified by ``id``. Usage::

    LOG(id, string)                // writes 'string' to the log file

    LOG(id, level, string)         // writes 'string' to the log file if 'level' <= 'p_LogLevel' in
                                   // Log::Start()

    LOG(id, class, level, string)  // writes 'string to the log file if 'class' is in 'p_LogClasses' in
                                   // Log::Start(), and 'level' <= 'p_LogLevel' in Log::Start()

Default class is ””; default level is 0.

Examples::

    LOG(SSEfileId, ”This is a log record”);
    LOG(OutputFile2Id, ”The value of x is ” << x << ” km”);
    LOG(MyLogfileId, 2, ”Log string”);
    LOG(SSEfileId, ”CHeB”, 4, ”This is a CHeB only log record”); Logging macros
==============

The logging macros provided by the :doc:`logging service <../services-logging-debugging>` are:


.. toctree::
   :maxdepth: 1

   ./services-logging-debugging-macros-log
   ./services-logging-debugging-macros-log_id
   ./services-logging-debugging-macros-log_if
   ./services-logging-debugging-macros-log_id_if

|br|
The logging macros described above are also provided in a verbose variant. The verbose macros function the same way as their non-verbose
counterparts, with the added functionality that the log records written to the log file will also be written to stdout. The verbose 
logging macros are::

    LOGV(id, ...)
    LOGV_ID(id, ...)
    LOGV_IF(id, cond, ...)
    LOGV_ID_IF(id, cond, ...)


|br|
A further four macros are provided that allow writing directly to stdout rather than a log file. These are::

    SAY( ...)
    SAY_ID( ...)
    SAY_IF(cond, ...)
    SAY_ID IF(cond, ...)

The ``SAY`` macros function the same way as their ``LOG`` counterparts, but write directly to stdout instead of a log file. The ``SAY`` 
macros honour the logging classes and level.
LOG_ID() macro
==============

::

    LOG_ID(id, ...)

Writes a log record prepended with calling function name to the log file specified by ``id``. Usage::

    LOG_ID(id)                        // writes the name of calling function to the log file

    LOG_ID(id, string)                // writes 'string' prepended with name of calling function to the
                                      // log file

    LOG_ID(id, level, string)         // writes 'string' prepended with name of calling function to the 
                                      // log file if 'level' <= 'p_LogLevel' in Log::Start()

    LOG_ID(id, class, level, string)  // writes 'string' prepended with name of calling function to the 
                                      // log file if 'class' is in 'p_LogClasses' in Log::Start(), and 
                                      // 'level' <= 'p_LogLevel' in Log::Start().

Default class is ””; default level is 0.

Examples::

    LOG_ID(Outf1Id);
    LOG_ID(Outf2Id, ”This is a log record”);
    LOG_ID(MyLogfileId, ”The value of x is ” << x << ” km”);
    LOG_ID(OutputFile2Id, 2, ”Log string”);
    LOG_ID(CHeBfileId, ”CHeB”, 4, ”This is a CHeB only log record”);
LOG_ID_IF() macro
=================

::

    LOG_ID_IF(id, cond, ...)

Writes a log record prepended with calling function name to the log file specified by ``id`` if the condition given by ``cond``
is met.

Usage: see :doc:`LOG_ID() <services-logging-debugging-macros-log_id>` and :doc:`LOG_IF() <services-logging-debugging-macros-log_if>`.LOG_IF() macro
==============

::

    LOG_IF(id, cond, ...)

Writes a log record to the log file specified by ``id`` if the condition given by ``cond`` is met. Usage::

    LOG_IF(id, cond, string)                // writes 'string' to the log file if 'cond' is TRUE.
    
    LOG_IF(id, cond, level, string)         // writes 'string' to the log file if 'cond' is TRUE and 
                                            // 'level' <= 'p_LogLevel' in Log::Start().
    
    LOG_IF(id, cond, class, level, string)  // writes 'string' to the log file if 'cond' is TRUE, 
                                            // 'class' is in 'p_LogClasses' in Log::Start(), and 
                                            // 'level' <= 'p_LogLevel' in Log::Start().

Default class is ””; default level is 0.

Examples::

    LOG_IF(MyLogfileId, a > 1.0, ”This is a log record”);
    LOG_IF(SSEfileId, (b == c && a > x), ”The value of x is ” << x << ” km”);
    LOG_IF(CHeBfileId, flag, 2, ”Log string”);
    LOG_IF(SSEfileId, (x >= y), ”CHeB”, 4, ”This is a CHeB only log record”);
Log::Open(...)
==============

::

    INT Log::Open(
        STRING  p_LogFileName,     // the name of the log file to be created and opened.
                                   // This should be the filename only – the path, prefix and
                                   // extensions are added by the logging service. If the file 
                                   // already exists, the logging service will append a version 
                                   // number to the name if necessary (see p_Append parameter below).
        BOOL    p_Append,          // flag indicating whether the file should be opened in append mode
                                   // (i.e. existing data is preserved) and new records written to the
                                   // file appended, or whether a new file should be opened (with 
                                   // version number if necessary).
        BOOL    p_TimeStamp,       // flag indicating whether timestamps should be written with each 
                                   // log file record.
        BOOL    p_Label,           // flag indicating whether a label should be written with each log 
                                   // record. This is useful when different types of logging data are 
                                   // being written to the same log file file.
        LOGFILE p_StandardLogfile  // (optional) the standard log file type for this log file.
                                   // If not provided LOGFILE::NONE is used.
    )

Opens a log file. If the ``p_Append`` parameter is TRUE and a file named ``p_LogFileName`` exists, the existing file will be opened and
the existing contents retained, otherwise a new file will be created and opened (not a Standard Log File – see 
:doc:`../services-extended-logging`).

The log file container directory is created at the path specified by the ``p_LogBasePath`` parameter passed to the 
:doc:`./services-base-level-logging-func-start` function. New log files are created in the log file container directory. ``BSE`` detailed 
log files are created in the ``Detailed_Output`` directory, which is created in the log file container directory (if required).

The filename is prefixed by the ``p_LogNamePrefix`` parameter passed to the :doc:`./services-base-level-logging-func-start` function.

If a file with the name as given by the ``p_LogFileName`` parameter already exists, and the ``p_Append`` parameter is false, a version
number will be appended to the filename before the extension (this functionality is largely redundant since the implementation of the 
log file container directory).

The integer log file identifier is returned to the caller. A value of −1 indicates the log file was not opened successfully. Multiple 
log files can be open simultaneously – referenced by the identifier returned.
Base-level logging
==================

The Log class member variables are private, and public functions have been created for logging and debugging functionality
required by the code.

The Log service can be accessed by referring to the Log::Instance() object. For example, to check if the logging service is
enabled, call the ``Log::Enabled()`` function::

    Log::Instance()→Enabled();

Since that could become unwieldy, there is a convenience macro to access the Log service. The macro just defines “LOGGING”
as “Log::Instance()”, so calling the ``Log::Enabled()`` function can be written as::

    LOGGING→Enabled();

The Log service must be initialised and started before logging and debugging functionality can be used. Initialise and start 
logging by calling the ``Log::Start()`` function::

    LOGGING→Start(...)

Refer to the description of the ``Log::Start()`` function below for parameter definitions.

The Log service should be stopped before exiting the program – this ensures all open log files are flushed to disk and closed
properly. Stop logging by calling the ``Log::Stop()`` function::

    LOGGING→Stop(...)

Refer to the description of the ``Log::Stop()`` function below for parameter definitions.


Log & debug record filtering
----------------------------

The Log service provides a set of functions and macros to manage log files, and to write log and debug records to the log files,
``stdout``, and ``stderr``. Base-level logging allows developers to tag log and debug records with a string ``class``, and an 
integer ``level``. The Log service will filter log and debug records by ``class`` and ``level``, and only write those records 
that meet the ``class`` and ``level`` filters specified by the users via the ``--log-level``,  ``--log-class``, ``--debug-level``,
and ``--debug-classes`` program options.


Log service functions:
----------------------

The Log service provides the following public member functions:


.. toctree::
   :maxdepth: 1

   ./services-base-level-logging-func-start
   ./services-base-level-logging-func-stop
   ./services-base-level-logging-func-open
   ./services-base-level-logging-func-close
   ./services-base-level-logging-func-write
   ./services-base-level-logging-func-put
   ./services-base-level-logging-func-debug
   ./services-base-level-logging-func-debugWait
   ./services-base-level-logging-func-error
   ./services-base-level-logging-func-squawk
   ./services-base-level-logging-func-say
   ./services-base-level-logging-func-enabled

Log::Debug(...)
===============

::

  BOOL Log::Debug(
      STRING p_DbgClass,  // the log class of the record to be written. An empty string (””) satisfies
                          // all checks against enabled classes.
      INT    p_DbgLevel,  // the log level of the record to be written. A value of 0 satisfies all 
                          // checks against enabled levels.
      STRING p_DbgStr     // the string to be written to stdout (and optionally to file).
  )

Writes ``p_DbgStr`` to stdout and, if logging is active and so configured (via program option ``--debug-to-file``), writes ``p_DbgStr`` 
to the debug log file.

Returns a boolean indicating whether the record was written successfully. If an error occurred writing to the debug log file, 
the log file will be disabled.
Log::Close(p_LogfileId)
=======================

::

    BOOL Log::Close(
        INT p_LogfileId  // the identifier of the log file to be closed (as returned by Log::Open()).
    )

Closes the log file specified by the ``p_LogfileId`` parameter. If the log file specified by the logFileId parameter is open, it is flushed
to disk and closed.

Returns a boolean indicating whether the file was closed successfully.
Log::Say(...)
=============

::

    VOID Log::Say(
        STRING p_SayClass,  // the log class of the record to be written. An empty string (””) satisfies
                            // all checks against enabled classes.
        INT p_SayLevel,     // the log level of the record to be written. A value of 0 satisfies all 
                            // checks against enabled levels.
        STRING p_SayStr     // the string to be written to stdout.
    )

Writes ``p_SayStr`` to stdout.

This function does not return a value.
Log::Put(...)
=============

::

    BOOL Log::Put(
        INT    p_LogfileId,  // the identifier of the log file to be written.
        STRING p_LogClass,   // the log class of the record to be written. An empty string (””) satisfies
                             // all checks against enabled classes.
        INT    p_LogLevel,   // the log level of the record to be written. A value of 0 satisfies all 
                             // checks against enabled levels.
        STRING p_LogStr      // the string to be written to the log file.
    )

Writes a minimally formatted record to the specified log file. If the Log service is enabled, the specified log file is active, and the log
record ``class`` and ``level`` passed are enabled, the string is written to the file. See :doc:`./services-base-level-logging` for details
regarding log record ``class`` and ``level``.

If labels are enabled for the log file, a label will be prepended to the record. The label text will be the ``p_LogClass`` parameter.

If timestamps are enabled for the log file, a formatted timestamp is prepended to the record. The timestamp format is **yyyymmdd hh:mm:ss**.

Returns a boolean indicating whether the record was written successfully. If an error occurred the log file will be disabled.
Log::Start(...)
===============

::

    VOID Log::Start(**
        STRING      p_LogBasePath,       // the name of the top-level directory in which log files will
                                         // be created.
        STRING      p_LogContainerName,  // the name of the directory to be created at p_LogBasePath to
                                         // hold all log files.
        STRING      p_LogNamePrefix,     // prefix for log file names (can be blank).
        INT         p_LogLevel,          // logging level (see below).
        STRING[]    p_LogClasses,        // enabled logging classes (see below).
        INT         p_DbgLevel,          // debug level (see below).
        STRING[]    p_DbgClasses,        // enabled debug classes (see below).
        BOOL        p_DbgToFile,         // flag indicating whether debug statements should also be 
                                         // written to the log file.
        BOOL        p_ErrorsToFile,      // flag indicating whether error messages should also be 
                                         // written to the log file.
        LOGFILETYPE p_LogfileType        // the log file type (see below).
    )

Initialises and starts the logging and debugging service. Logging parameters are set per the program options specified (using
default values if no options are specified by the user). The log file container directory is created. If a directory with the
name as given by the containerName parameter already exists, a version number will be appended to the directory name. The 
``Run_Details`` file is created within the log file container directory. Any input files specified by the user, such as a grid
file and/or a log file definitions file, are copied to the log file container directory.  Log files to which debug statements
and error messages should be written will be created and opened if required.

The filename to which debug records are written when parameter ``p_ErrorsToFile`` is TRUE is declared in ``constants.h`` – see
enum class ``LOGFILE`` and associated descriptor map ``LOGFILE_DESCRIPTOR``. Currently the name is ”Debug_Log”.

Log file types are defined in the enum class ``LOGFILETYPE`` in ``constants.h``. The log file type and file extension are defined
by the ``p_LogfileType`` parameter:

    .. list-table::
       :widths: 25 75 
       :header-rows: 0
       :class: aligned-text

       * - LOGFILETYPE::TXT
         - will result in a plain text file, delimited by the space character (' '), with a file extension of ”.txt”
       * - LOGFILETYPE::TSV
         - will result in a plain text file, delimited by the tab character ('\\t'), with a file extension of ”.tsv”
       * - LOGFILETYPE::CSV
         - will result in a plain text file, delimited by the comma character (','), with a file extension of ”.csv”
       * - LOGFILETYPE::HDF5
         - will result in an ``HDF5``\ [#f1]_ file, with a file extension of ”.h5”

This function does not return a value.


.. rubric:: Footnotes

.. [#f1] https://www.hdfgroup.org/
Log::Error(p_ErrStr)
====================

::

    BOOL Log::Error(
        STRING p_ErrStr  // the string to be written to stdout (and optionally to file).
    )

Writes ``p_ErrStr`` to stdout and, if logging is active and so configured (via program option ``--errors-to-file``), writes ``p_ErrStr`` 
to the error log file.

Returns a boolean indicating whether the record was written successfully. If an error occurred writing to the error log file, 
the log file will be disabled.
Log::Write(...)
===============

::

    BOOL Log::Write(
        INT    p_LogfileId,  // the identifier of the log file to be written.
        STRING p_LogClass,   // the log class of the record to be written. An empty string (””) satisfies
                             // all checks against enabled classes.
        INT    p_LogLevel,   // the log level of the record to be written. A value of 0 satisfies all 
                             // checks against enabled levels.
        STRING p_LogStr      // the string to be written to the log file.
    )

Writes an unformatted record to the specified log file. If the Log service is enabled, the specified log file is active, and the log 
record ``class`` and ``level`` passed are enabled, the string is written to the file. See :doc:`./services-base-level-logging` for details
regarding log record ``class`` and ``level``.

Returns a boolean indicating whether the record was written successfully. If an error occurred the log file will be disabled.
Log::Squawk(p_SquawkStr)
========================

::

    VOID Log::Squawk(
        STRING p_SquawkStr  // the string to be written to stderr.
    )

Writes ``p_SquawkStr`` to stderr.

This function does not return a value.
Log::DebugWait(...)
===================

::

    BOOL Log::DebugWait(
        STRING p_DbgClass,  // the log class of the record to be written. An empty string (””) satisfies
                            // all checks against enabled classes.
        INT    p_DbgLevel,  // the log level of the record to be written. A value of 0 satisfies all 
                            // checks against enabled levels.
        STRING p_DbgStr     // the string to be written to stdout (and optionally to file).
    )

Writes ``p_DbgStr`` to stdout and, if logging is active and so configured (via program option ``--debug-to-file``), writes ``p_DbgStr`` 
to the debug log file, then waits for user input.

Returns a boolean indicating whether the record was written successfully. If an error occurred writing to the debug log file, 
the log file will be disabled.
Log::Enabled()
==============

::

    BOOL Log::Enabled()

Returns a boolean indicating whether the Log service is enabled – TRUE indicates the Log service is enabled and available;
FALSE indicates the Log service is not enabled and so not available.
Log::Stop(p_ObjectStats)
========================

::

    VOID Log::Stop(
        TUPLE<INT, INT> p_ObjectStats  // number of stars or binaries requested, count created.
    )

Stops the logging and debugging service. All open log files are flushed to disk and closed (including any Standard Log Files
open - see description of Standard Log Files in :doc:`../services-extended-logging`. The ``Run_Details`` file is populated and closed.

This function does not return a value.
Random numbers service
======================

A Random Number service is provided, with the ``gsl`` Random Number Generator encapsulated in a singleton object (an instantiation of
the ``Rand`` class).

The ``Rand`` class member variables are private, and public functions have been created for random number functionality required by 
the code.

The Rand service can be accessed by referring to the ``Rand::Instance()`` object. For example, to generate a uniform random floating 
point number in the range [0, 1), call the ``Rand::Random()`` function::

     double u = Rand::Instance()→Random();

Since that could become unwieldy, there is a convenience macro to access the Rand service. The macro just defines "RAND" as 
"Rand::Instance()", so calling the ``Rand::Random()`` function can be written as::

    double u = RAND→Random();

The Rand service must be initialised before use. Initialise the Rand service by calling the ``Rand::Initialise()`` function::

    RAND→Initialise();

Dynamically allocated memory associated with the ``gsl`` random number generator should be returned to the system by calling the 
``Rand::Free()`` function::

    RAND→Free();

before exiting the program.

The Rand service provides the following public member functions:


.. toctree::
   :maxdepth: 1

   ./services-random-numbers-initialise
   ./services-random-numbers-free
   ./services-random-numbers-seed
   ./services-random-numbers-defaultSeed
   ./services-random-numbers-random1
   ./services-random-numbers-random2
   ./services-random-numbers-randomGaussian
   ./services-random-numbers-randomInt1
   ./services-random-numbers-randomInt2

Rand::RandomInt(p_Upper)
========================

::

    INT Rand::RandomInt(const INT p_Upper)
    
Returns a random integer number uniformly distributed in the range [0, ``p_Upper``), where 0 :math:`\small \leq` ``p_Upper``.
    
Returns 0 if ``p_Upper`` :math:`\small \lt` 0.Rand::RandomGaussian(p_Sigma)
=============================

::

    DOUBLE Rand::RandomGaussian(const DOUBLE p_Sigma)

Returns a Gaussian random variate, with mean 0.0 and standard deviation ``p_Sigma``.Rand:Random()
=============

::

    DOUBLE Rand::Random()

Returns a random floating point number uniformly distributed in the range [0.0, 1.0).Rand::Seed(p_Seed)
==================

::

    UNSIGNED LONG INT Rand::Seed(const UNSIGNED LONG p_Seed)

Sets the seed for the ``gsl`` random number generator to ``p_Seed``. The return value is the seed.Rand::Initialise()
==================

::

    VOID Rand::Initialise()

Initialises the ``gsl`` random number generator. If the environment variable ``GSL_RNG_SEED`` exists, the ``gsl`` random number generator
is seeded with the value of the environment variable, otherwise it is seeded with the current time.Rand::Free()
============

::

    VOID Rand::Free()

Frees any dynamically allocated memory.Rand::Random(p_Lower, p_Upper)
==============================

::

    DOUBLE Rand::Random(const DOUBLE p_Lower, const DOUBLE p_Upper)

Returns a random floating point number uniformly distributed in the range [``p_Lower``, ``p_Upper``), where ``p_Lower`` 
:math:`\small \leq` ``p_Upper``.

(``p_Lower`` and ``p_Upper`` will be swapped if ``p_Lower`` :math:`\small \gt` ``p_Upper`` as passed)Rand::DefaultSeed()
===================

::

    UNSIGNED LONG INT Rand::DefaultSeed()

Returns the ``gsl`` default seed ``gsl_rng_default_seed``.Rand::RandomInt(p_Lower, p_Upper)
=================================

::

    INT Rand::RandomInt(const INT p_Lower, const INT p_Upper)

Returns a random integer number uniformly distributed in the range [``p_Lower``, ``p_Upper``), where  ``p_Lower`` 
:math:`\small \leq` ``p_Upper``.

(``p_Lower`` and ``p_Upper`` will be swapped if ``p_Lower`` :math:`\small \gt` ``p_Upper`` as passed)Linux installation
==================

You will need to install the following packages (and their prerequisites) using your package manager:


                .. list-table::
                   :widths: 25 40 35
                   :header-rows: 1
                   :class: bordered
                
                   * - Package
                     - Ubuntu (apt)
                     - RHEL (yum)
                   * - g++
                     - g++
                     - gcc
                   * - GSL
                     - libgsl-dev
                     - gsl gsl-devel
                   * - Boost
                     - libboost-all-dev
                     - boost-devel
                   * - hdf5
                     - libhdf5-serial-dev
                     - hdf5-devel

For Ubuntu, type::

    sudo apt-get install g++ libboost-all-dev libgsl-dev libhdf5-serial-dev

For Fedora::

    sudo dnf install gcc boost-devel gsl gsl-devel hdf5-devel

For RHEL or CentOS::

    sudo yum install gcc boost-devel gsl gsl-devel hdf5-devel

Getting started
===============

COMPAS is a platform for the exploration of populations of single stars and compact binaries formed through isolated binary evolution. 
The COMPAS population synthesis code is flexible, fast and modular, allowing rapid simulation of binary star evolution. The complete 
COMPAS suite includes the population synthesis code together with a collection of tools for sophisticated statistical treatment of the 
synthesised populations.

To start using COMPAS, get a copy of the code, and install the libraries and tools required to build and run COMPAS:

.. toctree::
   :maxdepth: 1

   ./git-details
   ./COMPAS-dependencies
   ./building-COMPAS
   ./dev-git-workflow

Once you have completed the steps shown above, you're ready to run COMPAS. The :doc:`COMPAS User Guide <../User guide/user-guide>`
explains in detail how to run COMPAS, but to check that COMPAS is installed correctly, and to get a taste of what running COMPAS looks
like, you could try the :doc:`COMPAS Tutorial <../User guide/Tutorial/example-compas-run>`.

COMPAS code repository
======================

The public COMPAS code resides in a ``GitHub``\ [#f1]_ repository.  You will need access to ``GitHub``, and the ``git`` version 
control system installed.

If you do not have ``git`` installed, visit `Install Git <https://www.atlassian.com/git/tutorials/install-git>`__ and follow the instructions there.


Getting your first copy of the COMPAS source code
-------------------------------------------------

While you don't need a ``GitHub`` account to read the COMPAS ``GitHub`` repository, you will need an account to push your code 
changes to the repository. If you plan to contribute to the COMPAS code base, you will need to create a ``GitHub`` account.


To Fork or Not to Fork?
-----------------------

You can clone the COMPAS repository directly, or first create a fork of the repository and clone the fork. 

Forking a ``GitHub`` repository creates a copy of the repository on ``GitHub``, in your account.  The original repository, and the fork 
created, are linked on ``GitHub``.

Cloning a repository creates a copy of the repository on your local system.

See `The difference between forking and cloning a repository <https://github.community/t/the-difference-between-forking-and-cloning-a-repository/10189>`__ 
to learn more about the pros and cons of forking a ``GitHub`` repository.


If you choose to fork the COMPAS repository, use the “fork” button at the top right of the ``GitHub`` repository page.


Creating a clone of the GitHub repository
-----------------------------------------

Whether you forked the COMPAS repository or chose to work directly with the repository, you will need to clone (copy) the repository to 
your local system.

First, change to the directory on your local system within which you wish to store your local copy of the COMPAS source code
(throughout this documentation we use the example directory `~/codes`):

::

    cd ~/codes


If you have configured your ``GitHub`` account for ``SSH``, you can clone with:

::

    git clone git@github.com:USERNAME/REPONAME.git


If you have not yet configured your ``GitHub`` account with ``SSH``, you can clone over ``HTTPS``:

::

    git clone https://github.com/USERNAME/REPONAME.git


(Subsititute your GitHub username for "USERNAME", and the GitHub repository name for "REPONAME"
("COMPAS" if you did not create a fork))

At the completion of the command you will have a copy (clone) of the COMPAS ``GitHub`` repository in the `~/codes` directory on your 
local system.

You can also use the green "Code" dropdown on the ``GitHub`` repository page to clone the repository.



Setting your username and email address
---------------------------------------

Before you can push changes, you must ensure that your name and email address are set:

::

   cd ~/codes
   git config --global user.name "Fred Nurk"
   git config --global user.email "fred.nurk@aplace.adomain"


You should  now be ready to start using COMPAS!


.. rubric:: Footnotes

.. [#f1] https://github.com/
Installing dependencies
=======================

COMPAS requires a ``C++`` compiler, and the ``gsl``, ``boost``, and ``hdf5`` libraries.  ``Python`` is required for the COMPAS post-processing tools.


.. toctree::
   :maxdepth: 1

   ./COMPAS-dependencies-linux
   ./COMPAS-dependencies-macOS
   ./COMPAS-dependencies-python

Git Workflow for COMPAS software developers
===========================================



Contents of this document


`Introduction <#introduction>`__


`Getting Set Up <#getting-set-up>`__


`Day to Day Commands <#day-to-day-commands>`__


`Lifetime of a Project <#lifetime-of-a-project>`__


`COMPAS Git Workflow <#the-compas-git-workflow>`__


`Terminology <#terminology>`__ 




Introduction
============

Git & Github for COMPAS developers


For those who are unfamiliar, git and github are popular tools in the
software development community for sharing and collaborating on software
projects.

Git is a light-weight command line tool for maintaining different
versions of software locally, and sharing those versions to remote
servers.
Github is a website that stores git-managed projects and enables
developers to collaborate centrally on many projects.

It is a bigger topic than we can get into here, but if you are curious
you should `read more
here. <https://www.atlassian.com/git/tutorials/what-is-version-control>`__

Learning git is somewhat similar to learning a new language, and it
can be difficult to fully grasp the vocabulary when starting out (which
makes searching the internet for help significantly more challenging!).
Some of the most fundamental terms are `described
below <#terminology>`__ to assist new users.



Purpose of this document


The purpose of this document is to:

-  Help COMPAS users who are new to git to get set up,
-  Outline a consistent workflow for COMPAS developers in their
   day-to-day use of git, and
-  Provide some of the commands that are required for this workflow

Git is very powerful, so this is only a very small subset of the
available git commands.

This is, in some sense, a living document, meaning we are always open
to suggestions and criticism with the workflow, and seek only to find
the best option for everybody.
Please send any feedback to the `COMPAS User Google
Group <mailto:compas-user@googlegroups.com>`__.

With that said, all developers should commit to learning the agreed upon
workflow, to ensure consistency and protect against conflicts which may
derail development.



Outline of the COMPAS code repository


*Note:* If anything below doesn't make sense, try looking at the end of this document for relevant `Terminology. <#terminology>`__

COMPAS users who are not developers can download the source code from
the Main Repository, found at
`github.com/TeamCOMPAS/COMPAS <github.com/TeamCOMPAS/COMPAS>`__ (details
can be found below).
You will only need the default ``production`` branch and do not need
to worry about what branches are.

For developers, this repository (or 'repo') is considered "pristine",
meaning that any work done here should be in a mature stage.

The repository contains 2 permanent branches, ``production`` and
``dev``.
All other branches are either feature or hotfix branches, whose
purpose is to either introduce some new functionality or fix a bug,
respectively, and then be deleted.

Feature branches on the Main Repository (also called the Main Fork or
simply Main) should be ready to be tested by others.
The Main Fork is not a "sandbox" for new, experimental ideas.
You should `create your own fork <#fork-the-main-repo>`__ off of the
Main Repository if you want to have public-facing experimental work.

This approach to the repository and workflow below are based on the
`Feature Branch
Workflow <https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow>`__
which is in common use in industry.



Getting Set Up
==============

**Step-by-step directions for how to configure your local and remote git
repositories**

*COMPAS Users and Developers*


Set up git and create a Github account


#. If you have not already, go to `github.com <https://github.com/>`__
   and setup an account.

#. Check that you have a `working install of
   git. <https://www.atlassian.com/git/tutorials/install-git>`__

#. It is recommended, though not necessary, that you configure `Github
   with
   ssh <https://help.github.com/en/articles/connecting-to-github-with-ssh>`__
   as well.
   This will allow you to bypass frequent login verification requests.



Clone the COMPAS repository to your personal computer


Cloning a git repository is slightly different to downloading a code
database.
Git includes a history of changes made to a repo, and these changes
are included with a clone.
Cloning also includes pointers to the original (remote) repo, so that
changes to the remote can be easily imported to your local machine.

#. To clone the COMPAS repo, go to the public COMPAS Github page and
   click on the green ``code`` button (see below).
   You can copy the repo address for ssh if you have it configured,
   otherwise click ``Use HTTPS`` and copy that address.

.. raw:: html

    <p align="center">
    <img src="./media/git_clone_button.png" width="600" />
    </p>

#. In a terminal window, change directory into the location where you
   plan to install COMPAS, e.g ``cd ~/git_repos/``

#. Type the appropriate of the following two commands into a terminal
   window (pasting the repo name copied above).

   -  SSH: ``git clone git@github.com:TeamCOMPAS/COMPAS.git``

   -  HTTPS: ``git clone https://github.com/TeamCOMPAS/COMPAS.git``

#. To confirm that it worked, run the following two commands:

::

    cd COMPAS
    git branch

#. If the clone finished without error, you should see as output:

``* production``

At this point, if you do not plan to do any COMPAS development, you're
all set.
See `getting\_started.md <getting_started.md>`__ to see how to compile
and run COMPAS.
If you run into issues or would like to see new features implemented,
you can `contact us here. <mailto:compas-user@googlegroups.com>`__.
You can read on if you are curious, but if you are not invited to be a
collaborator, you will only have read-access to the repository.

*COMPAS Developers Only*


*Note:* This section is very technical.  Take a look at the section below on `Terminology, <#terminology>`__ if you get stuck!



Join as a collaborator


In order to contribute to COMPAS, you will need to be added as a
collaborator.
Non-collaborators have read-only access to all of the branches.

`Contact us here <mailto:compas-dev@googlegroups.com>`__ to inquire
about collaborating, or reach out to one of us directly (see the `COMPAS
homepage <https://compas.science/>`__ for an up-to-date list).



Fork the main repo


As a COMPAS developer, you are highly encouraged to create your own
personal github fork of the Main repo.
This is a second, remote repository (distinct from your local repo),
but is managed by your github account.
It serves as a public-facing 'sandbox' of your current work, where you
can share partially-developed ideas and projects with others who might
be interested in assisting, without interferring with or clogging up the
Main repo.

On Github, go to the TeamCOMPAS/COMPAS repo and click on ``Fork`` in
the upper-right corner.
This will create a copy of the current state of the TeamCOMPAS/COMPAS
repo, including all branches and all commit histories, and place it on
your profile, identified as ``<your-username>/COMPAS``.

Since this is your personal repo, you can be as organized or
scatter-brained as you wish here.
If you work best with 50 branches, obscure names, and code scraps
everywhere, have at it.
You can also give or take away access to any other collaborators who
you might wish to contribute.
Note that for public repositories, your code will still be read-only
for everyone who is not a collaborator.



Connect to your remote fork from your local repo


Once your fork is created, you'll want to connect it to your local
repository.
In the terminal, navigate to your COMPAS git repo and type:

``git remote add <fork-nickname> <remote-fork-url>``

The ``<remote-fork-url>`` can be found on your remote repo under the
same green 'Clone or Download' button as before.
If you have ssh configured, it will be similar to
``git@github.com:reinhold-willcox/COMPAS.git``.
The ``<fork-nickname>`` is your choice, but should be informative, e.g
``reinhold_fork``.



Day to Day commands
===================

Basic commands for navigating local git


Branches allow a developer to experiment with multiple new features
simultaneously on the same code-base.
In git, branches are very lightweight and easy to manage, making them
incredibly useful.

To view, switch, and create branches (akin to ``ls``, ``cd``, and
``mkdir``), use:

::

    git branch
    git checkout <branch-name>
    git checkout -b <new-branch>

*Note:* Many git commands require that you are on the correct branch before executing the command - using these 3 commands regularly before running more complicated commands will save you headaches down the road!

**Important:** A new branch is already created as a copy of current
branch, so you always need to double check that you're on the branch you
want to copy (typically, ``dev``).



Committing changes


What git does best is to record all the small changes and edits that
accumulate as we modify code.
After many small changes, you might have a feature that you decide
isn't actually what you want, and you want to get rid of it.
Or you might have introduced a bug at some point that spans many
files, and you need to remove it without undoing all the work you've
accomplished since then.
Git makes this incredibly easy by storing small edits as "commits".
Commits, like branches, are incredibly versatile and powerful, but can
be conceptually tricky to grasp at first.

Committing is the process of adding a selection of changes to the
history of your branch.
It is effectively saving your work, and should be done often (every
time any small fix has been made).
To perform a commit:

#. Check that you're on the correct branch!

``git branch``

#. Add the relevant files to your "index" (whatever files you've just
   edited)

``git add <file1> <file2> <...>``

#. Then submit the commit with a commit message. The message should be
   clear and concise, to help identify exactly when certain changes were
   made and undo them if necessary.

``git commit -m "really clear message indicating all the changes you made in this commit."``

*Note:* A single commit should capture an entire "fix" of one kind.

*Example:* Say you want to add a function to
a C file and its header, and you also want to update the internal
contents of a completely different function in the same C file, you
should do 2 commits.

#. First, make the edits to the first function and header, then add and
   commit

::

    git add file.C file.h
    git commit -m "created function myFunction to do someStuff and added it to the header file"

#. Then do the same for the second edits

::

    git add file.C
    git commit -m "updated internal contents of thisOtherFunction to allow for specificUseCase"

You can (and should) check the status of the current index regularly
with:

``git status``

The printout is pretty self explanatory: it tells you which files have
been added, and which have been changed that you might consider adding
before committing.

If you accidentally staged a file to your index, you can undo a
``git add`` before you have done a ``git commit`` with:

``git reset <file>``

You can also use ``git commit --amend`` to fix the most recent,
erroneous commit.

::

    git commit -m "previous commit which had the wrong files staged"
    git add <fogotten-file>
    git reset <file-that-does-not-belong>
    git commit --amend

which will open an editor where you can modify the commit message.

The takeaway message is that branches are made up of many commits
strung together, one after another, forming a history of minor edits to
a given branch.
You can view the commit history of a branch with any of:

::

    git log
    git log --pretty=oneline
    git log --pretty=medium --graph
    git log --all --decorate --oneline --graph

(If you have some spare time/ interest, there are actually quite a few
elaborate git log setups online you can look at for inspiration).

Looking through your git log, you may begin to appreciate the value of
clear, concise, commit messages.



Merging branches


Creating a branch for every new idea is great, but at some point
you'll have two somewhat-finalized, distinct features on different
branches that you will want to combine into one.
To do that, you need to merge the branches.
Merging a separate branch onto your current branch adds a 'merge
commit' to the tip of your current branch, and leaves the other branch
as it was.
The two original branches are called parent branches, and the result,
appropriately, the child.
Typically, once you successfully merge, it is desirable to delete the
separate branch to keep things tidy.

::

    git checkout branch1
    git merge branch2
    git branch -d branch2

Merging can be difficult at first because, unless you are very good at
thinking ahead or very lucky, you probably have some overlap in the two
branches that you were working on.
Git has some pretty clever tools to figure out which changes to pull
into the merge commit, but if it is ambiguous (e.g if you've edited the
same part of a file in both parents), you will get a merge conflict.
You will have to manually edit the files to choose how to resolve the
conflict.
You are encouraged to make backup copies of both parent branches until
you are more comfortable.
Git has several `ways to deal with merge
conflicts; <https://www.atlassian.com/git/tutorials/using-branches/merge-conflicts>`__
the best option for you may depend on the particular IDE you are using.



Comparing branches


Often it is useful to see differences between branches and workspaces
without actually making any changes to either.
To accomplish this, we use the ``git diff`` command.
The arguments (or lack thereof) determine which objects are compared.

To see all the recent changes in your working directory:

::

    git diff            # compare working directory to index
    git diff HEAD       # compare working directory to tip of the current branch
    git diff --cached   # compare index to tip of the current branch

To compare two branches ``<b1>`` and ``<b2>`` (or even a single file on
separate branches):

::

    git diff <b1> <b2>                         # compare the tips of two branches
    git diff <b1> <remote-fork>/<b2>           # compare local branch to a remote branch
    git diff <b1>:./file/path <b2>:./file/path # compare the same file on different branches

For even more flexibility and control over branch/file comparisons, you
should checkout ``git difftool`` and its customizations for your
preferred text editor.



Deleting branches


You should become comfortable deleting branches, or else your repos
might pile up with old branches that are no longer active.
Branches are also very easy to manage in git (relative to other
version control systems), so you should practice creating new branches,
making quick edits, committing, and deleting again without worry.
To delete a branch,

#. Navigate to any other branch

``git checkout <unrelated-branch>``

#. Try deleting the branch

``git branch -d <branch-name>``

#. If that throws an error, likely there were some uncommited/unmerged
   changes (work that would be completely lost if the branch gets
   deleted).
   Either commit/merge the branch before deleting, or if you don't want
   to keep the changes, you can force the delete with:

``git branch -D <branch-name>``



Fetch other branches from a remote


If you followed the above workflow, you can verify that the COMPAS
repo is a designated remote fork in your local repo, nicknamed
``origin``.
You can also see any other remote forks that you have linked from your
local repo:

``git remote -v``

should output something like:

::

    origin  git@github.com:TeamCOMPAS/COMPAS.git (fetch)
    origin  git@github.com:TeamCOMPAS/COMPAS.git (push)
    reinhold_fork   git@github.com:reinhold-willcox/COMPAS.git (fetch)
    reinhold_fork   git@github.com:reinhold-willcox/COMPAS.git (push)
    another_fork    git@github.com:another-user/COMPAS.git (fetch)
    another_fork    git@github.com:another-user/COMPAS.git (push)

To see all of the available branches across all your linked forks:

``git branch -a``

should output something similar to

::

    * production
    local_feature_branch
    remotes/another_fork/dev
    remotes/another_fork/production
    remotes/another_fork/runSubmit
    remotes/origin/HEAD -> origin/production
    remotes/origin/dev
    remotes/origin/production
    remotes/origin/release
    remotes/reinhold_fork/dev
    remotes/reinhold_fork/git_workflow
    remotes/reinhold_fork/production

where anything not starting with "remotes/" is a local branch, and the
\* indicates your current branch.

*Note:* The remote branch named ``origin/HEAD`` is a pointer to the ``origin/production`` branch.  HEAD, when used locally, is a pointer to the most recent commit, or "tip", of the current branch.  `Read more. <https://stackoverflow.com/questions/2529971/what-is-the-head-in-git>`__

All of the remote branches are available to be copied locally with:

``git checkout -b <new-local-branch-name> <remote-name>/<remote-branch-name>``

*Example:*

``git checkout -b myPySubmit another_fork/runSubmit``



Configuring remote tracking branches - pushing & pulling


**Important:** This section is crucially important, but it contains
some of the more confusing subtleties of git.
I tried to make these explicit throughout, but as a result this
section is a bit dense (sorry about that).
I highly recommend trying the commands yourself as you read through.

It's often useful, though not required, to point local branches to a
branch on a remote repo, from which it will inherit changes.
For example, when changes occur on the ``dev`` branch of the Main
repo, you will probably want to pull them into your local ``dev`` branch
to keep up to date.

If changes occur on the remote, your local git repo will not
automatically know about it (git does not regularly ping the remote
server with update requests like, e.g, most phone apps).
You can check for remote changes on a fork with:

``git fetch <remote-fork>``

*Warning:* This is a bit subtle - ``git fetch`` only updates
git's "local knowledge" of the remote branches, it does not affect your
local branches.
That makes it very "safe" - you can't overwrite any of your own work
with ``fetch``.
This is not true of ``git pull`` `(see below). <#git-pull>`__

To see which local branches are tracking remote branches, use:

``git branch -vv``

which will have an output that looks similar to:

::

    * compas_hpc_updates eea656f [origin/compas_hpc_updates: behind 14] Removed references to dead files:
    dev                  a110d38 [origin/dev: ahead 2, behind 12] Remove unwanted demo files (#150)
    production           d379be5 [origin/production] Jeff's defect repairs from previous commits that had to be readded (#82)
    new_branch           b6aee96 generic branch to test git branch -vv, don't keep this

#. The first column lists your local branches (the \* indicates your
   current branch).
#. The second column is the unique hash that identifies the commit of
   the tip of that branch (technically, it's only the beginning of the
   hash, but it suffices to identify the commit).
#. If the local branch is tracking a remote branch, this will be
   specified in brackets in the third column as
   ``[<remote_repo>/<remote_tracking_branch>]``.

   -  If there is a colon after the branch name with either "ahead N" or
      "behind M" (or both), this describes whether the tip of the local
      branch has additional commits that the remote does not, and vice
      versa.

#. If there are no brackets, the branch is not tracking anything.



git pull


If you have a branch which is "behind" the remote branch it is tracking
by some number of commits, then yours is out of date and you should
update it with:

::

    git checkout <outdated-branch>
    git pull

The ``git pull`` command defaults to the remote tracking branch of the
current branch (whatever was in the brackets above).
If the current branch is not tracking anything, or if you want to pull
from a different remote branch
(e.g, if ``origin/dev`` was updated and you want your
``<local-feature-branch>`` to pull in those updates), you can set it
explicitly:

::

    git checkout <local-feature-branch>
    git pull <remote-fork> <remote-branch>

*Note:* You should regularly check that your branches are updated. If not, you should pull to avoid larger conflicts later on.



git push


To share your local work with the other collaborators, you need to
"push" your changes to a remote repository.
Similar to ``git pull``, ``git push`` defaults to the designated
remote tracking branch, if it exists.
If not, or if you want to push to a different remote branch, you can
set it manually:

::

    git checkout <local-branch-to-push>
    git push <remote-fork> <remote-branch>

Pushing to your personal remote repository is a way to save all of
your commits (i.e the history of edits) somewhere off your local
computer.
This is good practice because it acts as a backup in the event
something happens to your local machine, and it also allows other
collaborators to see your work
(without having to explicitly send them your work all the time).
This should also be done often, but not necessarily for every commit.
A good rule of thumb is to push any updated branches at the end of the
day.



pull requests


We will briefly introduce here the concept of pull requests. If
working on a remote repo, especially a shared one, it is often desirable
to block direct push access, as this could
potentially lead to bad code being introduced without proper vetting.
The solution is pull requests: the user who wrote the new code will
submit the changes as a pull request,
for another developer to review. If they pass inspection, the reviewer
can then approve the pull request and merge the changes into the remote
repo.

Clarification of the difference between push, pull, and pull requests
can be found in the `Terminology <#terminology>`__ section below.



set remote tracking branch


You can add or update a branch's remote tracking branch (sometimes
called the "upstream" branch) with:

::

    git checkout <branch-to-update>`
    git branch --set-upstream-to=<remote-fork>/<remote-branch-to-track>

*Note:* The syntax may vary slightly depending on your version of git.  ``man git branch`` should be able to shed some light.



Lifetime of a New Feature
=========================

New feature branches


When beginning a new feature, you will typically want to branch off of
the most updated version of the ``dev`` branch.
Ultimately, the feature will be merged back into ``dev`` (or else
abandoned), and this will facilitate the merge later on.

::

    git checkout dev 
    git pull
    git checkout -b <new-project>

The name of your branch should *clearly* describe the feature you plan
to implemented.
This will help you to keep track of where different bits of code live
once the number of branches gets large.



Ongoing feature branches


Commit regularly as you make changes.

::

    git status
    git add <file1> <file2> <...>
    git commit -m "useful message"

When you have made many commits and want to push your work up to the
remote, first check that you have the correct current and target
branches

::

    git branch -vv
    git push

If you are working on a shared remote branch, you should also pull
regularly to keep up with any changes that are made there. A safe way to
check if there are any changes, without risking overwriting your local
work, is to fetch and diff.

::

    git fetch <remote-fork>
    git diff HEAD <remote-fork>/<remote-branch>



Finalized features


When a feature branch is nearing completion (e.g when the code is
nearly ready to be submitted into the Main Repository and tested), you
will want to ensure that it is fully up-to-date with the Main repo.
Then, push your branches up to your personal remote repo before
submitting a Pull Request.

#. Ensure that your branch has the latest updates from ``dev``.

::

    git checkout dev
    git pull
    git checkout <mature-branch>
    git merge dev

#. Push to your personal remote repo

::

    git checkout <mature-branch>
    git push --set-upstream <your-remote-repo>

#. Submit a Pull Request to the Main repo

   -  Login to github and go to your personal remote repo
      ``<your-username>/COMPAS``.

   -  Click ``Pull request`` (If you recently pushed your branch, you
      could also click on ``Compare & pull request``)

   -  Double check that you have selected the correct feature and target
      branches. In almost all cases, the base should be
      ``TeamCOMPAS/COMPAS`` with branch ``dev``, which will probably not
      be the default. Then click ``Create pull request``

   -  Add a comment describing your feature and what changes you made.
      If you have any particular reviewers in mind, or your feature
      solves one of the Git Issues, you should link those here. Then
      click ``Create pull request``, and you're all set!

.. raw:: html

    <p align="center">
    <img src="./media/git_pr_button.png" width="600" />
    </p>

Once you have created the pull request, it is up to the other team
members to review it (see below). They may ask you to fix some parts
before accepting it, so keep an eye on the pull request conversation.



The COMPAS Git Workflow
=======================

The above sections go over many of the available git commands that you
might find useful.
This section delves into how we apply these specifically to the COMPAS
workflow.

Overview


There should always be only 2 branches on the Main Repo:
``production`` and ``dev``.
They are both permanent, and both can only be modified with pull
requests which must be approved by another COMPAS developer.

The ``production`` branch is the current "long term service release",
meaning that it should be well-tested.
Of course, code is never truly bug-free, but this branch is the one
that the public will use, so updates should be extensively tested.

The ``dev`` branch is where new features are joined together in
preparation for the next release.
Pull requests to ``dev`` should be made from feature branches sitting
on other remote repos (e.g the personal repo of the author).
Presumably, these new features have each been tested in isolation and
correctly do what they propose to do.
But ``dev`` is a place to confirm that all the new features combined
together still produce sensible output.



Reviewing Pull Requests


Typically, a new feature branch will be formally reviewed when it is
submitted as a pull request into ``dev``.
Reviewers have a responsibility to check the following:

-  The code compiles without error on the usual assortment of Operating
   Systems.
-  The code runs without error using all default values (``./COMPAS``).
-  The code runs without error on a medium-sized population of binaries.
-  The new feature(s) do what they propose to do.
-  All new features are explicitly mentioned (i.e nothing is fixed
   quietly).
-  Documentation has been updated appropriately.
-  Formatting conforms to the rest of COMPAS.

This does not all have to be done by one reviewer, but there should be a
consensus among all reviewers that all tests have been passed.

A new release is defined by a pull request from ``dev`` to
``production`` and should involve most of the active developers.
The ``dev`` branch should be tested heavily for a variety of potential
bugs, including speed tests, different package and OS versions, and
comparisons of key plots from different papers.



Terminology
===========

-  **Commit**: A single commit records a collection of edits to one or
   more files, with an associated commit message.
   You can make and undo many changes before making a commit, and you
   can similarly revert commits which are later deemed unnecessary.
   As a verb, committing changes means to create a commit of the
   changes and append that commit onto a sequence of previous commits (a
   "branch", see below).

-  **Branch**: A single branch is an ordered sequence of commits.
   A new commit is always appended onto the tip of a branch, and the
   name of the branch is really just a pointer to this most updated
   branch tip.
   When a new branch is created from an old one, they initially still
   point at the same commit, the tip of both (currently identical)
   branches.
   New commits can be applied to one branch or the other, leading to a
   divergent history (which is not a bad thing).
   The imagery of the shared history of commits, followed by the split
   into two separate histories, readily leads to the name "branches".
   A branch will often represent a place to experiment with changes in
   a way that doesn't risk destroying the existing code.
   Major branches will add some new functionality or some new physical
   prescription, while sub-branches may pop-up to quickly test some
   variation to the new functionality.
   These sub-branches might be merged in to the major feature branch,
   destroyed, or possibly continue on their own to be expanded into a
   more major feature (and then merged in later on).
   Whether the branch is merged or scrapped, it should always
   `ultimately be deleted <#deleting-branches>`__
   `[1] <https://rickandmorty.fandom.com/wiki/Mr._Meeseeks>`__ (aside
   from the permanent ``production`` and ``dev`` branches).

.. raw:: html

    <p align="center">
    <a href="https://nvie.com/posts/a-successful-git-branching-model/">
    <img src="./media/git_branches.png" width="600" />
    </a>
    </p>

-  **Repository**: A Repository (or Repo) is a single storage location
   for a given code base.
   A single github user may have many repos for all of their different
   software projects.
   In our case, we have the Main Repository hosted by on Github at
   `TeamCOMPAS/COMPAS. <github.com/TeamCOMPAS/COMPAS>`__
   There are often many repositories for a given development project -
   these can be local or remote repositories (see below), each (usually)
   hosted by one the developers.
   Each repo can contain different branches each with slight
   variations on the code base, and these branches can be readily shared
   between repos, along with their history of commits.
   A Repo can be public (often called Open Source) or private.
   COMPAS is Open Source, but the general public has only read-access.
   Prospective contributors need to be added as a collaborator in
   order to make changes and submit pull requests.

-  **Local/Remote**: Local refers to the repository on your personal
   computer, while Remote refers to any repo that isn't.
   Github repos (whether Main or someone else's) will be remote for
   everyone.
   My local computer is only local to me; from a purely git
   perspective, it would be considered remote to anyone else, though
   this should not come up often because other users should never have
   even remote access to your personal computer.
   The purpose of your personal remote fork is to be a public proxy
   for your local fork, where you can add things you've worked on that
   you wish to share around.

-  **Fork**: A Fork is a full copy of a repo, including all its
   branches, to another location.
   Most of the time, "another location" will mean elsewhere on the
   github servers, since we will be Forking from the Main Repo to our
   Personal Repo when we are setting up.
   In our case, Forks will distinguish different users, or perhaps
   groups of users (e.g Copenhagen/COMPAS).
   All core developers should have a personal fork.
   If you are familiar with the ``git clone`` command, this is
   identical to Forking from a remote server onto your own personal
   computer.

-  **Origin**: Origin is the name commonly used for the primary remote
   repository.
   It is configured by default whenever you clone from a repository,
   so yours will probably point to the Main Repo.
   If you track multiple remote forks, you should give them all
   helpful, distinguishing names (e.g ``jeff_fork``, ``reinhold_fork``,
   etc.)

-  **Working Directory**: The Working Directory is where a user makes
   edits to files.
   It has meaning in git only in reference to the Index and the most
   recent commit (i.e the tip of the current branch).
   Files are editted in the working directory, before being added to
   the Index (or "staged"), and then finally committed to the current
   branch, or HEAD (see below).

-  **Index**: The Index (aka Staging Area) exists only in the
   intermediate step between editing local files and committing those
   files.
   Historically, other Version Control systems only allowed editting
   files, and then committing those files one by one.
   The issue with that is that sometimes a collection of edits of
   different files logically make up one full "commit-worthy-edit".
   The classic example of this is adding a function to a .C file and
   it's header .h file.
   If you need to revert this commit back for any reason, it makes
   sense to remove both of those edits at once - you would virtually
   never need to remove the function from the C file but leave it in the
   header.
   Adding files to the index is the way to collect all of the files
   that were involved in a given series of edits that you want to treat
   as one big Edit.

-  **Tracking**: The word tracking has two meanings, and could refer to
   either tracked remote repositories, or tracked local files in the
   current branch, and they have slightly different meanings.

   - A tracked repository is one which contains a branch which is
     currently being tracked, or "upstream", of a branch in your local
     repository.
     By default, all the branches on a forked repository track the
     branches they were forked from.
     You can modify the upstream branch of a given branch to point at
     any other branch you like, whether local or remote. You can also
     have multiple tracked remote repositories, though any given branch
     can only track at most one other branch at a time.
     This is useful if you want to check out and keep up-to-date with a
     branch that sits on a colleague's fork.
     You can view all tracked repositories with ``git remote -v``
   - A tracked file is one that git "knows about", meaning it was
     included in the last commit.
     You can have other files in the same folders as your git repo
     which are not tracked (if, e.g, you want to have output files from
     COMPAS runs but do not want to share those around).
     If you make modifications to a tracked file but don't commit it,
     git will not let you leave the branch.

-  **Push, Pull, and Pull Request**: These commands form the backbone
   of file-sharing across repositories.
   They all cover the same conceptual idea of "taking a branch and
   copying it over to a different branch on another repo." The
   difference is where you are relative to the target.
   You ``pull`` from a remote into your local, and you ``push`` from
   your local into a remote.
   For many remotes, there are protections in place to keep arbitrary
   users from pushing changes ad hoc.
   ``Pull-requests`` are the polite version of a ``push`` - instead of
   forcing your changes onto a remote, you are asking the manager of the
   remote to review your changes, and hopefully pull them into the
   remote if they approve.

-  **Revert**: A revert is used when it is decided that a particular
   previous commit (or perhaps several) have introduced bugs or are
   otherwise no longer undesired, and we want to remove them from the
   branch.
   A ``git revert`` will attempt to identify the changes made in those
   commits, and create a new commit to undo them.
   This is a fairly advanced git command and can easily become quite
   complicated, so make sure to use this one with caution, make backups
   of your work, and do lots of testing before you try anything.

-  **HEAD**: HEAD is a pointer to a commit, but the specific commit it
   points to moves around regularly.
   In general, it refers to the tip of whichever is the current
   branch.
   When you make a commit to the branch, HEAD updates to the new tip.

-  **Log**: The log of a branch is the history of that branch in terms
   of its commits.
   The log shows when the commits occured, who authored them, and what
   the commit message stated.


macOS installation
==================

It is strongly recommended to update to the latest version of macOS through the App Store. You can find what macOS version you are 
using by clicking on the Apple symbol on the top left of your screen and clicking "About This Mac".

The next step is to install or update Xcode. You can find it directly in the App Store or at `Xcode <https://developer.apple.com/xcode/>`__\ . 
Note: Xcode installation requires around 20 GB of disk space. If you are low on disk space, you may consider installing a ``C++`` 
compiler directly.
 
Once Xcode is installed, open a Terminal, and execute the following command to install the required command line developer tools::
 
    xcode-select --install

Next, you need to install several extra libraries and python modules. Popular ways of installing them are via package managers MacPorts and Homebrew. 
We give instructions for installing ``gsl``, ``boost``, and ``hdf5`` with Homebrew. To install Homebrew, run::

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

If the installation was successful, the following should run without error::

    brew --version

Now install ``gsl``, ``boost``, and ``hdf5`` using Homebrew by running::

    brew install gsl
    brew install boost
    brew install hdf5

Building COMPAS
===============

We first need to define an environment variable for the root directory of COMPAS in your shell start-up file for COMPAS to run properly. For example, 
if you use bash as your shell, open `~/.bashrc` with a text editor and put in the following::

    export COMPAS_ROOT_DIR=~/codes/COMPAS

where `~/codes` should be replaced with the path to the directory where you cloned the COMPAS repository. For this to take effect, either restart your 
bash session or run::

    source ~/.bashrc

If your shell is ``zsh`` (which is the default of macOS 10.15), set the environment variable as above in `~/.zshrc` instead of `~/.bashrc`. If your shell
is ``csh``, set the environment variable in `~/.cshrc` using::

    setenv COMPAS_ROOT_DIR ~/codes/COMPAS
    
Now go to the COMPAS source code directory::

    cd $COMPAS_ROOT_DIR/src

In this directory you will find the file ``Makefile``, which you need to edit to point to your ``gsl``, ``boost``, and ``hdf5`` include files and libraries. 

If you installed the packages with Homebrew, the package files are likely to be found in /usr/local/opt (in directories gsl, boost, and hdf5 respectively),
but if they are not found there you will need to use Homebrew to locate the files::

    $ brew info boost
    boost: stable 1.72.0 (bottled), HEAD
    Collection of portable C++ source libraries
    https://www.boost.org/
    /usr/local/Cellar/boost/1.72.0 (14,466 files, 648.5MB) *
    ...

Copy the path, which in this case is `/usr/local/Cellar/boost/1.72.0`, and add it to the appropriate lines of the Makefile::

    BOOSTINCDIR = /usr/local/Cellar/boost/1.72.0/include
    BOOSTLIBDIR = /usr/local/Cellar/boost/1.72.0/lib
 
To build the COMPAS executable (compile and link) type::

    make -f Makefile

The build process will run much faster if multiple processors/cores are available. To build the COMPAS executable using (e.g.) 4 cores, type::

    make -j 4 -f Makefile

Note that both ``make`` commands shown above will conduct incremental builds: they will only compile source files that have changed. To ensure a clean build
in which all source files are compiled, type::

    make clean
    make -j 4 -f Makefile

The `clean` option instructs ``make`` to remove all existing object files (.o), and the COMPAS executable.  A subsequent ``make`` is then forced to compile
all source files and link the resultant object files (and external libraries) into a new executable.

The executable can be tested with, e.g.,

    ./COMPAS -v

which will display the code version.

See :doc:`../Developer guide/Developer build/COMPAS-local-build` for a detailed description of ``Makefile`` functionality.


:bolditalictext:`A note for Mac users:`

If you are using MacOS and running into linking issues with the boost libraries, try::

    make clean
    make CPP=clang++ -j$(sysctl -n hw.ncpu)

In some Mac installations, the GNU C++ compiler is not installed how we might expect, so trying to compile and link with ``clang++`` might help.

Installing Python
=================

Python and some selected libraries are required for interfacing with the code, and also for post-processing. We recommend using ``python3``. The 
``matplotlib`` and ``numpy`` libraries should also be installed. The libraries ``scipy``, ``astropy``, and ``pandas`` are also used in some other scripts.

First check if you have ``python3`` installed. If you do, typing the following should give you the version number::

    python3 --version

If you do not have ``python3`` installed, install it by following the instructions below for your OS:

- For macOS, We recommend installing ``Python`` and its libraries using MacPorts. You can follow the instructions on `MacPorts Python installation on Mac <https://astrofrog.github.io/macports-python/>`__.
- For Linux, install `python3` using your package manager (e.g. in Ubuntu, run `sudo apt-get install python3`). We recommend installing the required python libraries using the package installer ``pip``. E.g. To install ``numpy``, run `pip install numpy`; for ``h5py``, run `pip install h5py`.

