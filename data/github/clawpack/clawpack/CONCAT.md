# Citing Clawpack

If you use Clawpack in publications, please cite the following:

    @misc{clawpack,
          title={Clawpack software}, 
          author={Clawpack Development Team}, 
          url={http://www.clawpack.org}, 
          note={Version 5.1},
          year={2014}}

If you are using a version of the software other than the one above, please make
sure to cite that instead.  Please also cite at least one of the following 
regarding the algorithms used in Clawpack:


#### Basic algorithms in 1D and 2D

    R. J. LeVeque, 1997. Wave propagation algorithms for multi-dimensional 
    hyperbolic systems. J. Comput. Phys. 131, 327–353.


    @article{rjl:wpalg,
        Author = {R. J. LeVeque},
        Title = {Wave propagation algorithms for multi-dimensional hyperbolic 
                 systems},
        Journal = {J. Comput. Phys.},
        Pages = {327--353},
        Volume = {131},
        Year = {1997}
    }

    R. J. LeVeque. Finite Volume Methods for Hyperbolic Problems. Cambridge 
    University Press, Cambridge, UK, 2002.

    @book{LeVeque-FVMHP,
          Author = {R. J. LeVeque},
          Title = {Finite Volume Methods for Hyperbolic Problems},
          Publisher = {Cambridge University Press},
          Year = {2002},
          Url = {http://www.clawpack.org/book.html}
    }

#### 3D algorithms
    
    J. O. Langseth and R. J. LeVeque. 2000. A wave-propagation method for 
    three-dimensional hyperbolic conservation laws. J. Comput. Phys. 165, 
    126–166.
    
    @article{LangsethLeVeque00,
             Author = {J. O. Langseth and R. J. LeVeque},
             Title = {A wave-propagation method for three-dimensional hyperbolic
                      conservation laws},
             Journal = {J. Comput. Phys.},
             Pages = {126--166},
             Volume = {165},
             Year = {2000}
    }

#### Adaptive Mesh Refinement (AMR)

    M. J. Berger and R. J. LeVeque. 1998. Adaptive Mesh Refinement using Wave-Propagation Algorithms for Hyperbolic Systems. SIAM J. Numer. Anal. 35, 2298–2316.

    @article{BergerLeVeque98,
             Author = {M. J. Berger and R. J. LeVeque},
             Journal = {SIAM J. Numer. Anal.},
             Pages = {2298--2316},
             Title = {Adaptive Mesh Refinement using Wave-Propagation Algorithms 
                      for Hyperbolic Systems},
             Volume = {35},
             Year = {1998}
    }
 
#### F-wave Algorithms

    D. S. Bale, R. J. LeVeque, S. Mitran, and J. A. Rossmanith. A wave-propagation 
    method for conservation laws with spatially varying flux functions, SIAM J. 
    Sci. Comput 24 (2002), 955-978.

    @article{BaleLevMitRoss02,
        Author = {D. Bale and R. J. LeVeque and S. Mitran and J. A. Rossmanith},
        Title = {A wave-propagation method for conservation laws and balance laws
                 with spatially varying flux functions},
        Journal = {SIAM J. Sci. Comput.},
        Pages = {955--978},
        Volume = {24},
        Year = {2002}
    }

#### GeoClaw

    M. J. Berger, D. L. George, R. J. LeVeque and K. M. Mandli, The GeoClaw 
    software for depth-averaged flows with adaptive refinement, Advances in Water 
    Resources 34 (2011), pp. 1195-1206.

    @article{BergerGeorgeLeVequeMandli11,
             Author = {M. J. Berger and D. L. George and R. J. LeVeque and K. T.  
             Mandli},
             Journal = {Adv. Water Res.},
             Pages = {1195-1206},
             Title = {The {GeoClaw} software for depth-averaged flows with 
                      adaptive refinement},
             Volume = {34},
             Year = {2011},
             Url = {\url{www.clawpack.org/links/papers/awr11}}
    }

    R. J. LeVeque, D. L. George, and M. J. Berger, 2011, Tsunami modelling with 
    adaptively refined finite volume methods, Acta Numerica, pp. 211-289.

    @article{mjb-dg-rjl:actanum2011,
             Author = {R.J. LeVeque  and D. L. George and M. J. Berger},
             Title = {Adaptive Mesh Refinement Techniques for Tsunamis and Other
                     Geophysical Flows Over Topography},
             Journal = {Acta Numerica},
             Pages = {211-289},
             Year = {2011}
    }

#### PyClaw

Please change the version number and year to the version you have used.

    @misc{pyclaw,
          title={PyClaw software}, 
          url={http://www.pyclaw.org}, 
          author={Mandli, Kyle T. and Ketcheson, David I. and others}, 
          note={Version 5.1}
          year={2014}
    }

    David I. Ketcheson, Kyle T. Mandli, Aron J. Ahmadia, Amal Alghamdi, Manuel 
    Quezada de Luna, Matteo Parsani, Matthew G. Knepley, and Matthew Emmett, 
    2012, PyClaw: Accessible, Extensible, Scalable Tools for Wave Propagation 
    Problems, SIAM Journal on Scientific Computing, 34(4):C210-C231

    @article{pyclaw-sisc,
             Author = {Ketcheson, David I. and Mandli, Kyle T. and Ahmadia, Aron 
             J. and Alghamdi, Amal and {Quezada de Luna}, Manuel and Parsani, 
             Matteo and Knepley, Matthew G. and Emmett, Matthew},
             Journal = {SIAM Journal on Scientific Computing},
             Month = nov,
             Number = {4},
             Pages = {C210--C231},
             Title = {{PyClaw: Accessible, Extensible, Scalable Tools for Wave    
             Propagation Problems}},
             Volume = {34},
             Year = {2012}
    }

#### SharpClaw (High Order WENO)

    D. I. Ketcheson, Matteo Parsani, and R J LeVeque, 2013, High-order Wave 
    Propagation Algorithms for Hyperbolic Systems, SIAM Journal on Scientific 
    Computing, 35(1):A351-A377 (2013)

    @article{Ketcheson2011,
             Author = {Ketcheson, David I. and Parsani, Matteo and LeVeque,
             Randall J.},
             Journal = {SIAM Journal on Scientific Computing},
             Number = {1},
             Pages = {A351--A377},
             Title = {{High-order Wave Propagation Algorithms for Hyperbolic  
                       Systems}},
             Volume = {35},
             Year = {2013}
    }
# Changes in Clawpack 5.0 (after 4.6.3)

See (http://clawpack.github.io/doc/clawpack5.html) for a changes to Classic, AMRClaw, and GeoClaw.

See (http://clawpack.github.io/doc/pyclaw/) for changes to PyClaw.
# Installing the Python Clawpack tools

The following command is sufficient to get up and running with PyClaw, VisClaw, and the Python interfaces to other packages:

    pip install clawpack

If you want to install clawpack in a manner that allows also using the
Fortran codes (Classic, AMRClaw, and GeoClaw), see
http://www.clawpack.org/installing.html

# Installation Tools

## Git

If you are just interested in working with development repositories, those are available from a checked out repository with the following command:

    python setup.py git-dev

Scroll down to "Working with a Development Version of Clawpack" for more details.  

## Python

We use the [numpy distutils](http://docs.scipy.org/doc/numpy/reference/distutils.html) as our fundamental install tool because we are using extended f2py extension module installers.

We really like the pip installer.  If you don't have pip installed and will be working with the full Clawpack stack, we recommend that you use [Anaconda](http://www.continuum.io/downloads) to quickly get up and running with a user-space Scientific Python stack.  

# Working with a Development Version of Clawpack

## Clawpack Packages as Git Submodules

Although the Clawpack packages have been set up as independent repositories, occasionally changes in one module will need accompanying changes in other modules.  In this situation, it is a good idea for maintainers to make a "meta-commit" to this top-level repository that contains the other repositories as submodules.  This allows maintainers to coordinate changes across multiple Clawpack packages.  An example of how to do this is shown later.

## Creating a read-only development version of Clawpack

```
git clone git://github.com/clawpack/clawpack.git
cd clawpack
python setup.py git-dev
# optionally, for installation of Python components
pip install -e .
```

This downloads all of the clawpack modules as subrepositories checked out at specific commits (as opposed to the tip of a branch).  

## Contributing to Clawpack!

If you'd like to contribute to Clawpack, you do not need to worry about the submodules.  Just checkout a branch in the package you are updating and push your changes as a pull request.  Just be sure to post your pull requests at the same time (and refer to the other pull requests) when submitting coordinated changes.

### Steps for contributors
* make local changes in one or more packages
* test your changes, if they're not testable, please add a test!
* commit and publish (push) your changes
* submit pull requests for each package you've modified

Optionally, you may also create a top-level clawpack branch and submit a pull request showing how your changes work together.  This is especially recommended if you are unable to create working tests otherwise.  Remember that you need to push commits in your submodules to GitHub before you can open a pull request that refers to those commits.

## Creating a development version of Clawpack for modification on GitHub/submitting pull requests

First, head to http://github.com/clawpack and fork the repositories that you will be working with.  You will probably need to fork one or more of the following:

* amrclaw
* clawpack
* clawutil
* geoclaw
* pyclaw
* riemann
* visclaw

Then, for each submodule you'd like to modify, add your own repository as a remote.  Here's an example for publishing changes to geoclaw (replace `username` below with your own username).

```
cd geoclaw
git remote add username git@github.com:username
# you may want to start from upstream
git checkout -b new_feature origin/master
# make some changes
git commit
# push your changes to your GitHub repository
git push username 
# open pull request on GitHub
```

## Working with a Maintenance Version of Clawpack

(This section is a work-in-progress)

# Maintaining Clawpack

The Clawpack maintainers are responsible for reviewing and accepting pull requests, although community support is always welcome in reviewing incoming pull requests.  Once a pull request has been accepted and tested, the maintainers need to decide if they'd like to update the top-level clawpack repository, which should be treated as the stable development branch coordinating the different repositories.  Changes in the individual repositories are not visible to somebody using `pip install` or the `python setup.py git-dev` commands unless they manually fetch/checkout changes or the top-level clawpack master branch has been updated.

## Complete Submodule Workflow Example

```
cd pyclaw
# make local changes
git commit -a -m "Updated pyclaw *cross-changes* with riemann."
# push to your GitHub fork of pyclaw
git push
cd ../riemann
# make local changes
git commit -a -m "Updated riemann *cross-changes* with pyclaw."
# push to your GitHub fork of riemann
git push
cd ..
# now at top-level clawpack directory
git add pyclaw riemann
git commit -m "Cross-update pyclaw and riemann."
# push to your GitHub fork of clawpack
git push

# open PRs on GitHub for all branches that need to be merged

```

More reading about submodules.
* [Submodules in the Git Community Book](http://book.git-scm.com/5_submodules.html)
* [Submodules in the Pro Git book](http://progit.org/book/ch6-6.html)
# Clawpack Community Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting one of the project team,
either [David Ketcheson](david.ketcheson@kaust.edu.sa), [Randy
LeVeque](rjl@uw.edu), or [Kyle Mandli](kyle.mandli@columbia.edu). All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org
