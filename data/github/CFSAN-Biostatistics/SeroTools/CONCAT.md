---
title: 'SeroTools: a Python package for *Salmonella* serotype data analysis' 
tags:
  - Python
  - Salmonella
  - serotype
  - serovar
  - White-Kauffmann-Le Minor
  - bioinformatics
authors:
  - name: Joseph D. Baugher, Ph.D.
    orcid: 0000-0003-2466-3047
    affiliation: 1
affiliations:
 - name: Center for Food Safety and Applied Nutrition, U.S. Food and Drug Administration
   index: 1
date: 01 July 2020
bibliography: paper.bib

---

# Summary

Subtyping, the ability to differentiate and characterize closely related microorganisms, 
has historically been a critical component of successful outbreak identification and 
traceback efforts employed by public health researchers and regulatory agencies for 
foodborne pathogens. Serological subtyping (or serotyping) has been the standard approach, largely based on antibody binding to surface antigens [@methods-micro]. The identification of specific antigenic factors has facilitated the creation of serotyping schemes, which define each serovar using a specific (generally unique) combination of antigenic factors. Serotyping schemes have been developed to assist in characterization of many microorganisms, including pathogens such as *Salmonella*, *E. coli*, *Shigella* [@clin-micro37], *Streptococcus* [@clin-micro22], and *H. influenzae* [@clin-micro36].

*Salmonella* is a major foodborne pathogen for which serotyping has played a fundamental monitoring role for over 50 years[@CDC2015]. *Salmonella* serotyping is generally based on antibody binding to the O antigen (a surface antigen) and one or more H antigen phases (flagellar antigens) [@clin-micro37; @BAM]. The White-Kauffmann-Le Minor (WKL) *Salmonella* scheme specifies the naming and formatting conventions for *Salmonella* serotyping data and the antigenic factors (and other characteristics) which define each serovar [@Grimont2007]. SeroTools includes the 2007 WKL scheme [@Grimont2007] and updates [@supp47; @supp48; @Bugarel2015]. 

The WKL scheme currently recognizes two species of *Salmonella*,  *S. enterica* and *S. bongori*. *S. enterica* is comprised of six subspecies (subsp.): *enterica* (I), *salamae* (II), *arizonae* (IIIa), *diarizonae* (IIIb), *houtenae* (IV) and *indica* (VI). Note that *S. bongori* is still frequently designated as subsp. V for scheme consistency, although it is no longer considered a subspecies of *S. enterica*. The WKL scheme assigns a unique name (e.g. serovar Enteritidis) to each of the serovars of *S. enterica* subsp. *enterica* (I), while the serovars representing the other subspecies are referred to by their antigenic formulae. The antigenic formula formatting is defined by the WKL scheme and is demonstrated for serovar Agona in Figure 1. The formula contains a subspecies designation and a colon-separated list of antigenic factors for which the following fields are required: O antigen, phase 1 H antigen, and phase 2 H antigen. The field for ‘Other H' antigen includes R phases and  third phases and is present only when populated. An antigenic formula may include additional annotation such as:

1. *Square* brackets to indicate optional factors, (e.g. I 1,4,**[5]**,12:f,g,s:**[1,2]**:**[z27]**,**[z45]**).
2. *Underlining* to indicate O factors present only in the presence of the converting phage, represented here and in SeroTools as optional (with *square* brackets) due to the inability to capture typographical formatting in plain text, (e.g. I **[1]**,9,12:e,h:1,5).
3. *Curly* brackets to indicate mutually exclusive factors, (e.g. I 3,**{10}{[15]}**:k:1,5).
4. *Parentheses* to indicate factors which are weakly agglutinable, (e.g. IIIb **(6)**,14:k:z53).
5. A *dash* to indicate a missing antigen, (e.g. I 1,9,12:g,m:**–**).

These additional annotations are captured in the SeroTools repository and employed for determination of congruence between serovars.

![Standard formatting of the antigenic formula.](antigenic_profile.png)

# Statement of Need

SeroTools addresses multiple critical needs for the efficient analysis of *Salmonella* serotyping data within the public health community. In recent years, significant technological advances have resulted in a wide range of molecular-based subtyping options, including highly sensitive approaches based on whole genome sequencing (WGS). One such approach involves the application of software tools to WGS data for *in silico* serovar prediction [@SeqSero; @SeqSero2; @shigatyper; @SerotypeFinder; @ECTyper; @hicap], including real-time prediction [@Feng2020]. SeqSero (a *Salmonella*-specific tool) and other *in silico* serovar designation tools have been adopted by U.S. public health agencies as an alternative to serological testing and for quality control applications [@Dowdy2017; @Timme2019]. The advent of new methodologies for serovar determination has engendered a need for method-comparison studies, and has sparked a growing collection of recent publications comparing various laboratory-based and *in silico* serovar predictions [@Cooper2020; @Banerji2020; @Diep2019; @Tang2019; @Ibrahim2018; @Yachison2017; @SeqSero2; @SeqSero]. In light of the growing interest in *in silico* serovar prediction and serotyping method-comparison studies, SeroTools provides unique tools which fill multiple gaps in the analysis process. It serves as the only multiformat WKL repository accessible for software development. Currently the WKL scheme is available only as a pdf document [@Grimont2007] and as Python lists in SeqSero [@SeqSero] and SeqSero2 [@SeqSero2]. SeroTools also provides the only existing tools for querying the WKL scheme, comparing serovars for congruence, and predicting the most abundant serovar for clusters of isolates.

# Functionality and Features

The SeroTools Python package provides the following functionality:

1. Repository –
    * SeroTools includes an updated WKL repository in multiple formats, including Python data structures (a pandas DataFrame, dictionaries, and lists) and spreadsheets (Excel and tab-delimited). The repository includes fields representing serovar name, antigenic formula, species, subspecies, O antigen, phase 1 H antigen, phase 2 H antigen, other H antigens, the new O group designation (e.g. O:2), and the old O group designation (e.g. A).
2. Toolkit –
    * **query** - SeroTools provides the ability to easily query the WKL repository with 
                serovar names or antigenic formulas.
    * **compare** - SeroTools provides a convenient method for automated comparison of
                  serovar designations, including increased differentiation for levels
                  of congruence.
    * **cluster** - SeroTools includes methods for robust determination of the most
                  abundant serovar for a cluster of isolates.
3. Additional functionality –
    * SeroTools includes Pythonic data structures and a host of utility functions for analyzing and manipulating large *Salmonella* serovar datasets. Other functionality includes the ability to determine the antigenic factors common to a group of serovars.

SeroTools defines four levels of congruence for use in querying the repository 
and comparing serovars. Note - *optional* factors as referenced below include optional, 
exclusive, and weakly agglutinable factors, as specified in the WKL scheme.

1. **Exact** matches must meet **one** of the following criteria:
    - The serovar designations are the identical string.

        For example:
	
            Corvallis                   Corvallis
            I 8,[20]:z4,z23:[z6]        I 8,[20]:z4,z23:[z6]

    - Every antigenic factor (*required* or *optional*) matches.

        For example:
	    
            Corvallis                   I 8,[20]:z4,z23:[z6]
            I 8,[20]:z4,z23:[z6]        I 8,20:z4,z23:z6
            I 1,3,10,19:f,g,t:1,(2),7   I 1,3,10,19:f,g,t:1,2,7

    - The subspecies designations are identical and neither serovar designation includes any antigenic factors.

        For example:
	
            I ::                        I –:–:–
            II :                        II –:

2. **Congruent** matches must meet **all** of the following criteria:
    - The subspecies field must be present either for both serovars or for neither.
    - All *required* antigenic factors match.
    - Any differences are due to the presence/absence of *optional* factors.

        For example:
	
            I 6,7,14:g,m,s:–            I 6,7,[14],[54]:g,m,[p],s:–
            I 6,7:g,m,s:–               I 6,7,[14],[54]:g,m,[p],s:[1,2,7]
            Amager var. 15+             Amager
            I 3,15:y:1,2:[z45]          I 3,{10}{15}:y:1,2:[z45]
            6,7:k:[z6]                  6,7:k:–	
	
3. **Minimally congruent** matches must meet the following criteria:
    - Every antigen of at least one serovar can be considered a formal subset of the 
corresponding antigen (no direct conflicts). Note - the empty set (–) is a subset of 
every set.

        For example:
	
            I 6,7,14,[54]:g,m,[p],s:–   6,7,[14],[54]:g,m,[p],s:–
            I                           I 6,7,8,[14],[54]:g,m,[p],s:–
            I 7:g:–                     I 6,7:g,m,s:–
            Gallinarum                  Enteritidis

4. **Incongruent** matches must meet the following criteria:
    - Any comparison which is not at least minimally congruent.

        For example:
	
            I                           II
            I 1:                        1 2:
            Javiana                     Saintpaul
            I 7,8:g,m,s:–               I 6,7,[14],[54]:g,m,[p],s:[1,2,7]
            I 4,5:a,b:6,7               I 5:a,b,c:6,7

The 'minimally congruent' designation is unique to SeroTools and is useful for 
distinguishing between two scenarios: serovars which differ due to sample misannotation 
(truly incongruent) and serovars derived from correctly annotated samples with variation 
based solely on missing information. When comparing serovar predictions, minor 
differences may be expected due to method-specific irregularities, for example, reagent 
variation for laboratory-based techniques or sequencing read coverage for *in silico* 
techniques. Our assumption is that these minor method-specific differences are more likely manifested as missing data (e.g. all but one of the correct factors were detected) than direct conflicts.


# Links

Documentation:
https://serotools.readthedocs.io/en/latest/readme.html

Source Code:
https://github.com/CFSAN-Biostatistics/serotools

PyPI Distribution:
https://pypi.python.org/pypi/serotools


# References
============
Contributing
============

.. highlight:: bash

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/CFSAN-Biostatistics/SeroTools/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/CFSAN-Biostatistics/SeroTools/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `SeroTools` for local development.

1. Fork the `SeroTools` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/serotools.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv serotools
    $ cd serotools/
    $ pip install sphinx_rtd_theme    # the documentation uses the ReadTheDocs theme
    $ pip install pytest
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the tests, including testing other Python versions with tox::

    $ flake8 serotools tests
    $ pytest -v
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Update the documentation and review the changes locally with sphinx::

    $ make docs

7. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

8. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.5, 3.6, 3.7, and 3.8.

Tips
----

To run a subset of tests::

    $ pytest -v tests/test_serotools.py
=======
Credits
=======

Author
----------------

* Joseph D. Baugher, Ph.D. <joseph.baugher@fda.hhs.gov>

Development Lead
----------------

* Joseph D. Baugher, Ph.D. <joseph.baugher@fda.hhs.gov>

CFSAN Bioinformatics Team
-------------------------

* Joseph D. Baugher, Ph.D. <joseph.baugher@fda.hhs.gov>

External Contributors
---------------------


.. :changelog:

History
=======

0.2.1 (2020-09-04)
---------------------

* Updated documentation
* Added JOSS manuscript


0.2.0 (2020-02-17)
---------------------

Significant updates in this version - not backwards compatible.

* The underlying data structures have been converted to pandas Series and DataFrames.
* New 'cluster' subcommand functionality provides the most abundant serovar(s) for clusters of isolates. 
* The 'predict' subcommand functionality has been merged into the 'query' subcommand, such that the default query will return any exact, congruent, and minimally congruent matches unless only exact matches are desired.
* The WKL repository is now available as a pandas DataFrame, in addition to dictionaries and lists.


0.1.1 (2019-11-27)
---------------------

* Corrected a variable name in cli.py
* Updated the algorithm for minimally congruent serovars


0.1.0 (2019-11-19)
---------------------

* Initial version.
===============================
SeroTools
===============================


.. Image showing the PyPI version badge - links to PyPI
.. image:: https://img.shields.io/pypi/v/serotools.svg
        :target: https://pypi.python.org/pypi/serotools

.. Image showing the Travis Continuous Integration test status, commented out for now
.. .. image:: https://img.shields.io/travis/CFSAN-Biostatistics/serotools.svg
..        :target: https://travis-ci.org/CFSAN-Biostatistics/serotools

.. Image showing the JOSS paper badge
.. image:: https://joss.theoj.org/papers/10.21105/joss.02556/status.svg
   :target: https://doi.org/10.21105/joss.02556

This package serves as a toolkit and repository for the White-Kauffmann-Le Minor scheme for *Salmonella* serotyping, which is made available in multiple formats, along with methods for querying and comparing serovar names and antigenic formulae, as well as determining the most abundant serovar for a cluster of isolates.

SeroTools was developed by the United States Food and Drug Administration, Center for Food 
Safety and Applied Nutrition.

* Free software
* Documentation: https://serotools.readthedocs.io
* Source Code: https://github.com/CFSAN-Biostatistics/serotools
* PyPI Distribution: https://pypi.python.org/pypi/serotools

Introduction
------------

*Salmonella* bacteria are major foodborne pathogens estimated by the U.S. Centers for Disease Control and Prevention to cause 1.35 million infections annually in the United States [1]_. Serological subtyping (serotyping) of *Salmonella* has historically been a critical component of characterization and successful outbreak identification and traceback efforts employed by public health researchers and regulatory agencies. The White-Kauffmann-Le Minor (WKL) *Salmonella* serotyping scheme specifies the commonly accepted naming and formatting conventions for *Salmonella* serotyping data and the antigenic factors (and other characteristics) which define each serovar. *Salmonella* serotyping data is routinely employed by a broad range of scientific researchers, physicians, public health professionals, food safety experts, etc.

SeroTools addresses multiple critical needs for the efficient analysis of *Salmonella* serotyping data. As technological advances continue to produce a range of high resolution subtyping options, including *in silico* serovar prediction based on whole genome sequencing, new tools are necessary for efficient method-comparison studies and quality control applied to increasingly large numbers of isolates. SeroTools serves as the only multiformat WKL repository accessible for software development and provides the only existing tools for querying the WKL scheme, comparing serovars for congruence, and predicting the most abundant serovar for clusters of isolates.

.. [1] The U.S. Centers for Disease Control and Prevention. <https://www.cdc.gov/salmonella/index.html>.


Features
--------

* Query the White-Kauffmann-Le Minor *Salmonella* serotyping repository

* Compare serovar predictions for state of congruence

* Determine the most abundant serovar for a cluster of isolates


Citing SeroTools
--------------------------------------

If you use SeroTools, please cite the following publication:
   
Baugher, J. D., (2020). SeroTools: a Python package for Salmonella serotype data analysis. Journal of Open Source Software, 5(53), 2556, https://doi.org/10.21105/joss.02556.



License
-------

See the LICENSE file included in the SeroTools distribution.




========
Notes
========

.. highlight:: none

1. Phage conversion factors denoted with underlining in [1]_ are here denoted as optional '[]' with the exception of the exclusive factors (e.g. {15} and {15,34}).
2. Serovar Montevideo is listed twice in [1]_:  O:7 I 6,7,[14]:g,m,[p],s:[1,2,7] and O:54 I {6,7,[14]}{54}:g,m,s:–. The profile from the 'Alphabetical List' p. 137 will be used here - I 6,7,[14],[54]:g,m,[p],s:[1,2,7].
3. As in [1]_, although S. bongori is not a subspecies of S. enterica, symbol 'V' was retained in order to avoid formatting confusion. 

.. [1] Grimont PA, Weill FX. Antigenic Formulae of the Salmonella Serovars. 9th. Paris, France: WHO Collaborating Center for Reference and Research on Salmonella, Institut Pasteur; 2007 <https://www.pasteur.fr/sites/default/files/veng_0.pdf>.
.. include:: ../HISTORY.rst
========
Usage
========

.. highlight:: none

SeroTools provides methods for querying and comparing serovar names and antigenic formulae, 
as well as determining the most abundant serovar for a cluster of isolates.

.. _query-label:

query
-----

Query the White-Kauffmann-Le Minor (WKL) repository by submitting one of more 
serovar names or antigenic formulas in an input file composed of a single query per line:: 

    $ serotools query -i <input_file>
    
or as a command line argument::

    $ serotools query -s 'Paratyphi A'
    
Output::

    Input        Name         Formula             Match
    Paratyphi A  Paratyphi A  I [1],2,12:a:[1,5]  exact

.. _compare-label:

compare
-------

Compare serovar predictions by evaluating multiple states of congruence (exact, congruent,
minimally congruent, incongruent). Serovar names and/or antigenic formulae may be submitted 
in a tab-delimited input file composed of two columns of serovar predictions::  

    $ serotools compare -i <input_file>

or as command line arguments::

    $ serotools compare -1 'Hull' -2 'I 16:b:1,2'

Output::

    Serovar1    Name    Formula     Serovar2    Name    Formula     Result
    Hull        Hull    I 16:b:1,2  I 16:b:1,2  Hull    I 16:b:1,2  exact

.. _cluster-label:

cluster
-------
Determine the most abundant serovar(s) for one or more clusters of isolates. Input data must be 
submitted in the form of a tab-delimited file in which each line consists of a cluster ID and one serovar as follows::

Input File - example.txt::

    cluster1	Dunkwa
    cluster1	Dunkwa
    cluster1	Utah
    cluster2	Hull
    
::

    $ serotools cluster -i example.txt
    
Output::

    ClusterID   ClusterSize Input   Name    Formula      P_Exact  P_Congruent	P_MinCon
    cluster1    2           Dunkwa  Dunkwa  I 6,8:d:1,7  0.6667   0.6667        0.6667
    cluster2    1           Hull    Hull    I 16:b:1,2   1.0      1.0           1.0
    
===========
Repository
===========

.. highlight:: none

SeroTools provides a repository of the White-Kauffmann-Le Minor (WKL) Salmonella serotyping scheme based on these :ref:`references-label` in the following formats:

- Python data structures (`serotools.py <https://github.com/CFSAN-Biostatistics/SeroTools/blob/master/serotools/serotools.py>`__)

  - pandas DataFrame:: 
  
      wklm_df
    
  - Dictionaries::
  
      wklm_name_to_formula
      wklm_formula_to_name
    
  - Lists with common indexing::
  
      wklm_name
      std_wklm_name (standardized for matching)
      wklm_formula
      std_wklm_formula (standardized for matching)
      wklm_sp (species)
      wklm_subsp (subspecies)
      wklm_O
      wklm_P1
      wklm_P2
      wklm_other_H
      wklm_group (O group)
      wklm_old_group (previous O group)
    
- An Excel spreadsheet (`White-Kauffman-LeMinor-Scheme.xlsx <https://github.com/CFSAN-Biostatistics/SeroTools/blob/master/wklm_scheme/White-Kauffman-LeMinor-Scheme.xlsx>`__)

- A tab-delimited text file (`White-Kauffman-LeMinor_scheme.tsv <https://github.com/CFSAN-Biostatistics/SeroTools/blob/master/wklm_scheme/White-Kauffman-LeMinor_scheme.tsv>`__)
==================
Statement of Need
==================

.. highlight:: none

SeroTools addresses multiple critical needs for the efficient analysis of *Salmonella* serotyping data within the public health community. In recent years, significant technological advances have resulted in a wide range of molecular-based subtyping options, including highly sensitive approaches based on whole genome sequencing which are being adopted by public health agencies for quality control and as an alternative to serological testing. In light of the growing interest in *in silico* serovar prediction and serotyping method-comparison studies, SeroTools provides unique tools which fill multiple gaps in the analysis process. It serves as the only multiformat White-Kauffmann-Le Minor (WKL) repository accessible for software development. SeroTools also provides the only existing tools for querying the WKL scheme, comparing serovars for congruence, and predicting the most abundant serovar for clusters of isolates.
===========
Congruence
===========

.. highlight:: none

SeroTools evaluates multiple levels of congruence for comparisons between serovar designations. 


exact
-----

Exact matches must meet one of the following criteria:

- Two serovar designations are the identical string::

    Corvallis                  Corvallis
    I 8,[20]:z4,z23:[z6]       I 8,[20]:z4,z23:[z6]
- Every antigenic factor (**required** or **optional**) matches::

    Corvallis                  I 8,[20]:z4,z23:[z6]
    I 8,[20]:z4,z23:[z6]       I 8,20:z4,z23:z6
    I 1,3,10,19:f,g,t:1,(2),7  I 1,3,10,19:f,g,t:1,2,7
- Neither serovar designation includes any antigenic factors, and the subspecies designations match::

    I ::                       I –:–:–
    II :                       II –:  


congruent
---------

Congruent matches must meet the following criteria:

- The subspecies field must be present for both serovars or neither.

- All **required** antigenic factors match. For example::

    I 6,7,14:g,m,s:–          I 6,7,[14],[54]:g,m,[p],s:–
    I 6,7:g,m,s:–             I 6,7,[14],[54]:g,m,[p],s:[1,2,7]
    Amager var. 15+           Amager
    I 3,15:y:1,2:[z45]        I 3,{10}{15}:y:1,2:[z45]
    6,7:k:[z6]                6,7:k:–                       


minimally congruent
-------------------

Minimally congruent matches must meet the following criteria:

- Every antigen of at least one serovar can be considered a formal subset of the corresponding antigen (no direct conflicts). For example::

    I 6,7,14,[54]:g,m,[p],s:–     6,7,[14],[54]:g,m,[p],s:–
    I                             I 6,7,8,[14],[54]:g,m,[p],s:–
    I 7:g:–                       I 6,7:g,m,s:–
    Gallinarum                    Enteritidis
* Note - the empty set (–) is a subset of every set

The minimally congruent designation is unique to SeroTools and is useful for distinguishing between two scenarios: 

- Serovars which differ due to sample misannotation (incongruent)

- Serovars derived from correctly annotated samples with variation based solely on missing information. When comparing serovar designations, minor differences may be expected due to method-specific irregularities, for example reagent variation for laboratory-based techniques or the presence of nonproductive genomic data when comparing antigenic agglutination to *in silico*-based techniques. Our assumption is that these minor method-specific differences are more likely manifested as missing data (e.g. all but one of the correct factors were detected) than direct conflicts. 


incongruent
-----------
Any comparison which is not minimally congruent. For example::

    I                             II     
    I 1:                          1 2:
    Javiana                       Saintpaul
    I 7,8:g,m,s:–                 I 6,7,[14],[54]:g,m,[p],s:[1,2,7]
    I 4,5:a,b:6,7                 I 5:a,b,c:6,7
.. include:: ../AUTHORS.rst
.. include:: ../README.rst
============
Installation
============

.. highlight:: bash

At the command line::

    $ pip install --user serotools

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv serotools
    $ pip install serotools


Upgrading SeroTools
-----------------------------------------

If you previously installed with pip, you can upgrade to the newest version from the command line::

    $ pip install --user --upgrade serotools


Uninstalling SeroTools
--------------------------------------------

If you installed with pip, you can uninstall from the command line::

    $ pip uninstall serotools
.. SeroTools documentation master file, created by
   sphinx-quickstart on Tue Jul  9 22:26:36 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SeroTools
===============================

Contents:

.. toctree::
   :maxdepth: 2

   readme
   installation
   statement_of_need
   formatting
   usage
   congruence
   repository
   notes
   references
   authors
   history

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. _references-label:

========
References
========

.. highlight:: none

SeroTools includes serovar names and antigenic formulas as specified in the following publications:

1. Grimont PA, Weill FX. Antigenic Formulae of the Salmonella Serovars. 9th. Paris, France: WHO Collaborating Center for Reference and Research on Salmonella, Institut Pasteur; 2007 <https://www.pasteur.fr/sites/default/files/veng_0.pdf>.
2. Guibourdenche M, Roggentin P, Mikoleit M, Fields PI, Bockemühl J, Grimont PA, Weill FX. Supplement 2003-2007 (No. 47) to the White-Kauffmann-Le Minor scheme. Res Microbiol. 2010 Jan-Feb;161(1):26-9 <https://doi.org/10.1016/j.resmic.2009.10.002>.
3. Issenhuth-Jeanjean S, Roggentin P, Mikoleit M, Guibourdenche M, de Pinna E, Nair S, Fields PI, Weill FX. Supplement 2008-2010 (no. 48) to the White-Kauffmann-Le Minor scheme. Res Microbiol. 2014 Sep;165(7):526-30 <https://doi.org/10.1016/j.resmic.2014.07.004>.
4. Bugarel M, den Bakker HC, Nightingale KK, Brichta-Harhay DM, Edrington TS, Loneragan GH. Two Draft Genome Sequences of a New Serovar of Salmonella enterica, Serovar Lubbock. Genome Announc. 2015 Apr 16;3(2) <https://doi.org/10.1128/genomeA.00215-15>. 
    
=================
Input Formatting
=================

.. highlight:: none

 
The input data must follow the nomenclature conventions as specified in [1]_.

- Named serovars (subsp. *enterica*) generally compose a single word or concatenation (e.g. Saintpaul) with no whitespace, with a few exceptions (e.g. Paratyphi B, II Alsterdorf), and should not contain hyphens or additional information (e.g. Choleraesuis is correct, while Cholerae-suis and Cholerae suis are incorrect).

- Named variants specified in [1]_ must adhere to the format below. For example, Westhampton var. 15+ or Westhampton var. 15+,34+:: 

    <SerovarName> var. <factor>+,<factor>+

- The proper convention for an antigenic formula is as follows: a space must separate the subspecies symbol from the antigens; colons must separate the antigens; commas must separate the antigenic factors (e.g. I 1,4,[5],12:e,h:1,5:[R1…])::

    <Subspecies> <O_factor1,O_factor2>:<P1_factor1,P1_factor2>:<P2_factor1,P2_factor2>
    <Subspecies> <O_factor1,O_factor2>:<P1_factor1,P1_factor2>:<P2_factor1,P2_factor2>:<otherH_factor1,otherH_factor2>

- Missing antigens should be specified using '–' (e.g. I 4,12,27:b:– or I 1,9,12:–:–). 

- Optional, exclusive, and weakly agglutinable factors should be designated as follows::

    optional            '[]'
    exclusive           '{}'
    weakly agglutinable '()'
    
- Input may contain the terms Nonmotile, Rough, or Mucoid, which will be converted into the appropriate format::

    Nonmotile  :-:-
    Rough      -:
    Mucoid     -:

- Phage conversion factors denoted by underlining in [1]_ are denoted as optional '[]' in SeroTools, with the exception of exclusive phage conversion factors (e.g. {15} and {15,34}). 


.. [1] Grimont PA, Weill FX. Antigenic Formulae of the Salmonella Serovars. 9th. Paris, France: WHO Collaborating Center for Reference and Research on Salmonella, Institut Pasteur; 2007 <https://www.pasteur.fr/sites/default/files/veng_0.pdf>.
