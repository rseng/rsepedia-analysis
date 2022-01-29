# Eir: Simulate Epidemics Using Spatial Models in Python

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03247/status.svg)](https://doi.org/10.21105/joss.03247)

![alt text](https://cdn.shopify.com/s/files/1/1879/3511/files/Eir_grande.jpg?v=1552614189)


Eir, named after the Norse valkyrie with great medical skill, is an API that allows the user to conduct stochastic simulations of epidemics, primarily using spatial models. With this software, one can simulate not only how epidemics relate to the distances between an infectious and susceptible indivdual, but also how the movement on infectious individuals plays a role in the spread of a disease. Eir also offers a lot of variety to the user, containing many more compartmental models that is present in any of the existing packages similar to Eir, including hospitalizations and vaccinations. Eir's usefulness can clearly be seen in modern day, where simulations and models are constantly used to form policy to combat COVID-19.
## Dependencies
Eir depends on ```numpy```, ```pandas```, ```matplotlib```, and ```multipledispatch```.
## Installation

One can install Eir via PyPI by running the following command via the command line:

```pip install Eir ```
The dependencies will automatically be installed as well.
## Notable Features
Eir offers countless different compartmental models, including:
- SIS
- SIR
- SIRS
- SIRD
- SIRV
- SIRSD
- SIRSV
- SIRDV
- SIRSDV
- SEIR
- SEIRS
- SEIRD
- SEIRV
- SEIRSD
- SEIRSV
- SEIRDV
- SEIRSDV
- ICU models. 

Eir also offers these models in different spatial models, some with mobility and some static.

## Examples

If one were to model the ICU hospitalizations using the Hub Model, the code could look as follows:

```python
from Eir import PeriodicICUV

test = PeriodicICUV(S0=999, E0=0, I0=1, R0=0, V0=0, rho=.3, ioda=.3, gamma=.25, mu=0.007, omega=.14, phi = .42, chi=.15, kappa=.05, eta=.02, spread_r=2, sigma_r=.25, move_R=4, sigma_R=.75, side=33, days=31)       
test.run()
test.plot()
```
In the above code segment:
  S0 : int
            The starting number of susceptible individuals in the simulation.
        
        E0: int
            The starting number of exposed individuals in the simulation.
        
        I0: int
            The starting number of infected individuals in the simulation.

        R0: int
            The starting number of recovered individuals in the simulation.
        
        V0: int
            The starting number of vaccinated individuals in the simulation.
        
        rho: float
            The probability of an individual leaving the E compartment.
        
        ioda: float
            The probability that, given an individual is leaving the E compartment, he goes to L compartment. The probability of that person going to I compartment is (1-ioda).
        
        gamma: float
            The probability of a person in I compartment going to the R compartment
        
        mu: float
            The probability of going from I to D, given that the person didn't go from I to R.
        
        phi: float
            The probability of going from L compartment to ICU compartment.
        
        chi: float
            The probability of going from ICU compartment to R compartment.
        
        omega: float
            The probability of going from ICU compartment to D compartment, given that individual didn't go from ICU compartment to R compartment.
        
        kappa: float
            The probability of going from R compartment to S compartment.
        
        eta: float 
            The probability of going from S compartment to V compartment, given that the individual didn't go from S compartment to E compartment. 
        
        timeDelay: float
            The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. 
        
        spread_r: the mean of the normal distribution of spreading radii that is use to generate spreading radii for each individual in the simulation.

        sigma_r: the standard deviation of the normal distribution of spreading radii that is used to generate spreading raddi for each individual in the simulation.

        move_R: the mean of the normal distribution of spreading radii that is use to generate movement radii for each individual's periodic movement in the simulation.

        sigma_R: the standard deviation of the normal distribution of spreading radii that is use to generate movement radii for each individual's periodic movement in the simulation.

        side: the length of the side of the square plane that individuals are confined to during the simulation.

        days: the number of days being simulated. 


To understand the variables and their meaning for different models, the documentation can be found in the docs folder in this repository, or looking at the docstrings in python. Additionally, if more detailed information about transmission chains and state histories was required, the methods from the Simul_Details class would allow the user to get a more in-depth look at the dynamics of the simulation.

## Contributors
The author welcomes and encourages new contributors to help test ``` Eir``` and add new functionality. If one wishes to contact the author, they may do so by emailing mjacob1002@gmail.com. Response times may vary.

## How to Cite

BibTex:
@article{Jacob2021,
  doi = {10.21105/joss.03247},
  url = {https://doi.org/10.21105/joss.03247},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {62},
  pages = {3247},
  author = {Mathew Jacob},
  title = {Eir: A Python Package for Epidemic Simulation},
  journal = {Journal of Open Source Software}
}

Citation String: Jacob, M., (2021). Eir: A Python Package for Epidemic Simulation. Journal of Open Source Software, 6(62), 3247, https://doi.org/10.21105/joss.03247


# Eir: Python Package for Simulating Epidemics

Eir is a package that focuses on allowing the user to simulate epidemics using spatial models. Eir mainly focuses on stochastic simulations, but offers a few deterministic models as well. The four spatial model types are: Strong Infectious, Hub Model, Periodic Mobility Model, and Random Movement model. Each model makes different assumptions about the spread and movements of the people in the simulation. The Strong Infectious and Hub Model are taken from researchers Fujie & Odagaki, so it is a good idea to read their paper at https://pubmed.ncbi.nlm.nih.gov/32288078/. The Random Movement model assumes that at each state, an individual is equally likely to move in any direction for a random distance, pulled from a normal distribution. The Periodic Mobility models attemtpts to capture the routine in people's lives by modelling each individuals movement as movement along a circle.
---
title: 'Eir: A Python Package for Epidemic Simulation'
tags:
  - Python
  - epidemiology
  - super spreaders
  - spatial models
  - stochastic
  - Python
authors:
  - name: Mathew Jacob
    orcid: 0000-0001-5513-432X
    affiliation: "1" # (Multiple affiliations must be quoted)

affiliations:
 - name: RxCovea, New York University, New York, USA
   index: 1

date: 8 April 2021

bibliography: paper.bib

---

# Summary

` Eir ` is a Python package that simulates epidemics with various assumptions and compartmental models. By mainly using spatial models, ` Eir ` allows for better modeling of the spread of an epidemic in smaller communities and localities. There are four different spatial models used: the Strong Infectious Model [@Fujie&Odagaki], Hub Model [@Fujie&Odagaki], the Random Movement Model, and the Periodic Mobility Model. While the Strong Infectious Model and Hub Model are static models, where no movement pattern is simulated, the Random Movement and Periodic Mobility Model simulates the movement pattern of people using different assumptions. The Random Movement Model assumes that an individual is equally likely to move in any direction for a random distance, generated from a user-defined normal distribution. The Periodic Mobility Model assumes that individuals move in a circle to capture the routine nature of people's movements. Therefore, after randomly generating locations, the model calculates a center of a circle and individuals move along that circle. With all of these spatial models, researchers can use the `Simul_Details` objects to interact with the simulation data and get a deeper understanding of the simulation's dynamics, such as transmission chains, state histories, etc. The data from these simulations can be plotted using Matplotlib or be returned in the form of a pandas DataFrame. It is important to note that these models are stochastic and are discrete-time Markov chains. 

` Eir ` leverages many different compartmental models in addtition to different spatial models in order to give researchers flexibility and robustness when simulating the spread of an epidemic. The models include but are not limited to: Susceptible-Infected-Recovered (SIR), Susceptible-Infected-Recovered-Dead (SIRD), Susceptible-Exposed-Infected-Recovered (SEIR), Susceptible-Exposed-Infected-Recovered-Vaccinated (SEIRV), Susceptible-Infected-Recovered-Susceptible-Vaccinated (SIRSV), and ICU models. The ICU models, in particular, are based on the work of Mühlpfordt et al [@ICU]. Additionally, for any model that includes vaccinations, there is a timeDelay variable that allows the user to delay vaccine distribution by a particular amount of days. 

While not the primary focus of the Python package, `Eir` also allows for researchers to use determinstic models to simulate basic compartmental models, including SIR, SEIR, SIRV, and more. 

`Eir` is currently being used for a publication involving simulating hospitalizations and studying vaccinations strategies. It is intended to be used by epidemiologists and researchers who are interested in understanding the dynamics of epidemics when operating under different modelling assumptions, such as periodic movement or random movement, which can be evaluated by using different models in Eir.

# Statement of Need

When modelling the spread of diseases, epidemiologists often use compartmental models, which classify people into states, the most common one being the SIR Model [@SIR]. There are certain rules that govern these compartmental models, such as that an individual can only exist in one compartment at any given time. With these models, epidemiologists are able to understand the dynamics of how an epidemic spreads. However, a lot of these models often use ODEs to approximate parameters for large populations and make assumptions to simplify calculations, which don't capture the complexities of smaller communities, where national averages may not be accurate. Then, researchers often use discrete-time Markov chains [@Allen], which generate stochastic events at certain time steps, such as every day or every hour. Additionally, in these smaller localities, spatial models are often deployed to get more accurate transmission probabilities. The sheer number of spatial models in Eir, as well as Eir's inclusion of models with movement, will make the process of modelling local disease transmission far easier. These spatial models are known as the Hub model [@Fujie&Odagaki], Strong Infectious Model [@Fujie&Odagaki], Random Movement Model, and Periodic Mobility Model. The Random Movement Model assumes that an individual is equally likely to move in every direction. The Periodic Mobility Model assumes that each person moves along a circle with a normally distributed radius. These models offer different structure to the spatial models, which can be used by epidemiologists. 


# Dependencies
NumPy, Matplotlib, multipledispatch, pandas

# Related Packages
There are different software package that enable users to simulate epidemics in different ways. A couple of such packages are:

### Epidemics on Networks

 Epidemics on Networks [@Miller2019], abbreviated EoN, is a Python package that allows users to simulate SIS and SIR epidemics on networks. It handles many different simulation types, as well. More can be found by reading the documentation at https://epidemicsonnetworks.readthedocs.io/en/latest/. Eir allows for a wider variety of models than EoN, which only includes SIS and SIR. Additionally, another difference is that Eir allows for the movement of individuals within the simulation, while the structure in EoN does not. 

### EpiModel

 EpiModel [@EpiModel] is an R package that allows the user to simulate models such as SIS and SIR. The details of this model can be found at https://www.epimodel.org/. EpiModel offers network based transmission, which is a feature that is not included in Eir. However, similar to EoN, there is no clear way to incorporate the movement of individuals within the simulation. Additionally, Eir offers a more detailed look at the simulation through the Simul_Details object, which can include information about individual state history, transmission chains, and more. 

# Acknowledgements
The author would like to acknowledge the contributions of Xueyao Guo to the models. Additionally, the author would like to acknowledge and thank Dr. Ernest Battifarano, Dr. Jantine Broek, and Professor Shirshendu Chatterjee for their mentorship and guidance through the author's research. Finally, the author would like to thank Professor Bud Mishra and the RxCovea group for guiding the author through his research. 

This project received no funding. 

# References

# Overview

These static models are ones in which the randomly generated (x,y) coordinate of each person in the simulation does not change throughout the entire simulation. Therefore, because there is no call to a _move() function, these simulations tend to run a little faster than the movement models. The two static spatial models are derived from the work of Ryo Fujie and Takashi Odagaki, whose paper proposed two different models: the Hub Model and the Strong Infectious Model. These two models attempt to make assumptions about super spreaders and how they propogate the disease differently than a normal spreader.

The default formula for a normal spreader's probability of infecting a susceptible is:

w(r) = 
    w0(1-r/r_0)^alpha, 0 <= r <= r_0,
    0, r > r_0.

In this equation, w0 is the probability that an infectious person infects someone when they are 0 units away, alpha is a constant, r is the distance between the infectious individual and the susceptible, and r_0 is the spreading radius of the infectious individual.

In the following two models, however, the formula for a super spreader changes.
## Hub Model
The Hub Model assumes that super spreaders are those who interact with more people, and therefore have a higher a spreading radius that is a scaled by a constant, allowing them to reach more people. Therefore, if we call k the scaling factor and r_0 the spreading distance for a normal spreader, the equation for a super spreader's probability of infecting a susceptible is:

w(r) = 
    w0(1-r/r_n)^alpha, 0 <= r <= r_n,
    0, r > r_n
, where r_n = k*r_0. 

### Parameters

These parameters of all Hub objects. don't include compartment-specific parameters, which can be found in the Compartments section of the documentation.

    pss: float
        probability someone is considered a super spreader upon generating the simulation.
    
    rstart: float
        the spreading radius of every normal spreader.
    
    side: float
        size of one side of the square plane that individuals are confined to.
    
    days: int
        The number of days that are simulated.
    
    w0: float (optional)
        The probability of infection if an infectious and susceptible individual are in the same location. Default is 1.0.
    
    hubConstant: float (optional)
        The factor k multliplied to the rstart if the person is a super spreader. Default is sqrt(6).
    
    alpha: float optional
        constant used in the P(infection) formula. Default is 2.0.


## Strong Infectious Model
The Strong Infectious Model assumes that super spreaders are intrinsically more infectious, and therefore have a fixed probability of spreading the disease over spreading radius identical to that of a normal spreader. Therefore, the formula for a super spreader's probability of propogating an infectious disease is: 
\begin{align*}
w(r) = \begin{cases}
    w0, 0 <= r <= r_0,
    0, r > r_0.
\end{cases}
\end{align*}.

### Parameters

The parameters of the Strong Infectious Model are almost identical to that of the Hub Model, except for two differences. The first is that there is no scaling factor, k, because the radius of the super spreader hasn't increased. Secondly, the default for w0 is now 0.7 rather than 1.0.# Overview

Unlike the two static spatial models, the Random Movement and Periodic Mobility models change the (x,y) coordinates of each individual at every time step according to the rules for that model. 

The same formula used for normal spreaders in static spatial models is the same for all individuals in the movement models. However, the spreading radii and movement radii for everyone in the simulation is picked from normal disributions generated from user input. However, what they do with these parameters depends on the model.

## Random Movement

In the Random Movement model, at each time step, the new (x,y) coordiante is picked by randomly picking a distance from a normal distribution and then randomly picking an angle from a uniform distribution because we assume that people in this model are equally likely to move in any direction. 

## Periodic Mobility

In the Periodic Mobility Model, at each step, the new (x,y) coordiante is picked by randomly picking an angle to move from a normal distribution, and then using the individual's movement radius R, move along the individual's circular motion, calculated upon generating the individual. This is aimed to capture the routine nature of people's movements. 

## Parameters

move_R: float
    The mean of the distribution of movement radii of a a person in the simulation. Used when genereating the movement radius of each individual in the simulation.

sigma_R: float
    The standard deviation of the distribution of movement radii of a a person in the simulation. Used when genereating the movement radius of each individual in the simulation.

spread_r: float
    The mean of the distribution of spreading radii of a person in simulation. Used when generating the spreading radius of each individaul in the simultion. 

sigma_r: float
    The standard deviation of the distribution of spreading radii of a normal spreader.


side: float
    The length of one side of the square plane that the individuals are confined to.

days: int
    The number of days that the simulations lasts for.

alpha: float, optional
    A constant used in the formula for calculating probability of infection between infectious person and susceptible person. Default is 2.0.

w0: float, optional
    The probability of a susceptible being infected by an infectious individual if they are 0 units apart. Default is 1.0.

There are also special parameters for the Periodic Mobility Model:

k: float
    Determines the mean of the distribution of thetas. The mean of that distribution is 2π/k.

std: float
    The standard deviation of the normal distribution of thetas. # Spatial Models

## Note
Because most of the functions for all the spatial models operate the same way after initializing the object, the examples will done using HubSIS. 

## Examples

### Spatial Model object
After initializing the spatial model, the simulation must be run in order to get the final outputs. For this, the run() function is used.

```python
    from Eir import HubSIS
    # initialize the object
    test = HubSIS(S0=999, I0=1, pss=.2, rstart=3, side=25, days=31, gamma=.3)
    # run the object
    details = test.run()
```

By default, run() returns a Simul_Details object, which allows the user to get deeper information about the simulation, such as transmission chains and state histories. If that isn't wanted, set getDetails=False in the run function.

```python
    from Eir import HubSIS
    test = HubSIS(S0=999, I0=1, pss=.2, rstart=3, side=25, days=31, gamma=.3)
    test.run(getDetails=False)
```

Now that the simulation has been run, there are different things that can be done with the simulation. If the data is wanted in a pandas DataFrame, the toDataFrame() method will return the data in the form of a pandas DataFrame.

```python 

    df = test.toDataFrame()
```

When printing the dataframe, the output will be something like this:

```shell
>>> print(df)
    Days  Susceptible  Infected
0    0.0        999.0       1.0
1    1.0        996.0       4.0
2    2.0        982.0      18.0
3    3.0        949.0      51.0
4    4.0        830.0     170.0
5    5.0        647.0     353.0
6    6.0        463.0     537.0
7    7.0        275.0     725.0
8    8.0        246.0     754.0
9    9.0        224.0     776.0
10  10.0        256.0     744.0
11  11.0        215.0     785.0
12  12.0        213.0     787.0
13  13.0        222.0     778.0
14  14.0        240.0     760.0
15  15.0        224.0     776.0
16  16.0        225.0     775.0
17  17.0        247.0     753.0
18  18.0        261.0     739.0
19  19.0        240.0     760.0
20  20.0        235.0     765.0
21  21.0        223.0     777.0
22  22.0        230.0     770.0
23  23.0        225.0     775.0
24  24.0        225.0     775.0
25  25.0        210.0     790.0
26  26.0        225.0     775.0
27  27.0        214.0     786.0
28  28.0        228.0     772.0
29  29.0        213.0     787.0
30  30.0        226.0     774.0
31  31.0        250.0     750.0

```

### Simul_Details

In the case that getDetails=True, interacting with the Simul_Details object can prove insightful. One useful function is the personHistory function, which will get the state history of the person number on each day. The persons are number from 0 to one less than the population size. For example, if there were 1000 people, the people would be identified as Person 0, Person 1, Person 2... Person 998, Person 999. 

```python

>>> from Eir import HubSIS
>>> test = HubSIS(S0=999, I0=1, pss=.2, rstart=3, side=25, days=31, gamma=.3)
>>> details = test.run()
>>> print(details.personHistory(9))
[(0, 'S'), (4, 'I'), (7, 'S'), (8, 'I'), (9, 'S'), (10, 'I'), (12, 'S'), (13, 'I'), (18, 'S'), (19, 'I'), (21, 'S'), (22, 'I'), (23, 'S'), (24, 'I'), (25, 'S'), (26, 'I'), (27, 'S'), (28, 'I'), (29, 'S'), (30, 'I'), (31, 'S')]
>>> 

```
The personHistory method returns a list of tuples. The 0th element of the tuple is the day, and the second element of the tuple is a string denoting the state they switched to that day. 

There is a default variable movementHistory which can be set to True, which then returns a list containing the state history in the 0th element, and a list of (x,y) coordinates of the particular person in the 1st element. Below is a script that one can run to access that information

```python
>>> from Eir import HubSIS

>>> # initialize object
>>> test = HubSIS(S0=999, I0=1, pss=.2, side=25, days=31, gamma=.3, rstart=3)
>>> # run the simulation
>>> d = test.run()
>>> # ge the list in the 1st element of the value returned by personHistory
>>> a = d.personHistory(9, True)[1]
>>> print(a)
[(13.215897854222156, 18.688142252818587)]
```

One can also view the people that a Person u infected by using the personTransmissionHistory() method. By inputting a number, u, the method will return of tuples that represent the person that u infected in the 0th element of each tuple, as well as the day it occured, the 1st element of each tuple.

```python 

>>>b = d.personTransmissionHistory(9)
>>>print(b)
[(305, 7), (369, 7), (636, 7), (544, 8), (897, 8), (897, 12), (369, 13), (636, 15), (897, 15), (82, 16), (224, 16), (169, 17), (11, 18), (408, 19), (564, 22), (224, 24), (95, 26), (601, 27), (636, 27), (957, 27), (95, 28), (601, 29), (957, 30)]
```
If we look at this list, we can see that perosn 9 infected person 305 on day 7, person 369 on 7, etc.

If one is interested in the transmission history on a particular day of the entire system rather than an individual, the transmissionHistoryOnDay() method will prove useful. 

```python

>>> c = d.transmissionHistoryOnDay(9)
>>> print(c)
[(2, 274), (2, 470), (3, 255), (3, 536), (6, 568), (9, 155), (9, 520), (9, 528), (9, 535), (9, 588), (9, 648), (9, 997), (10, 194), (10, 260), (10, 943), (11, 347), (11, 421), (11, 644), (11, 650), (11, 985), (13, 387), (14, 1), (16, 267), (16, 298), (16, 610), (16, 899), (17, 8), (17, 174), (17, 186), (17, 191), (17, 195), (17, 512), (17, 530), (17, 573), (17, 661), (17, 729), (17, 861), (19, 557), (19, 915), (21, 958), (21, 983), (22, 643), (22, 698), (23, 18), (23, 60), (23, 83), (23, 137), (23, 673), (24, 342), (24, 442), (24, 612), (24, 718), (24, 773), (24, 809), (24, 816), (24, 857), (24, 863), (25, 7), (25, 140), (25, 345), (26, 203), (26, 844), (27, 821), (28, 89), (28, 178), (28, 538), (28, 575), (30, 497), (30, 609), (31, 111), (31, 166), (31, 376), (31, 814), (33, 653), (33, 741), (35, 756), (35, 895), (36, 426), (36, 508), (37, 122), (37, 128), (37, 169), (37, 196), (37, 304), (37, 965), (38, 364), (38, 772), (39, 822), (39, 930), (40, 201), (40, 606), (40, 611), (40, 807), (40, 840), (40, 992), (41, 386), (41, 790), (42, 443), (42, 923), (43, 391), (44, 129), (44, 764), (49, 542), (50, 56), (50, 124), (50, 135), (50, 216), (50, 268), (50, 486), (50, 724), (50, 738), (52, 71), (52, 305), (52, 554), (54, 45), (54, 357), (54, 371), (54, 377), (54, 655), (54, 910), (55, 394), (58, 626), (58, 675), (58, 734), (58, 869), (62, 330), (62, 980), (66, 205), (66, 973), (68, 412), (68, 441), (70, 549), (72, 107), (72, 296), (72, 502), (72, 815), (72, 927), (73, 988), (74, 663), (77, 707), (80, 34), (80, 113), (80, 900), (81, 161), (81, 189), (88, 0), (90, 665), (90, 744), (92, 158), (92, 995), (93, 117), (93, 192), (93, 352), (93, 550), (93, 692), (93, 702), (95, 320), (96, 258), (100, 200), (100, 471), (100, 926), (108, 289), (108, 479), (112, 901), (118, 168), (118, 343), (118, 960), (119, 338), (119, 545), (131, 299), (133, 327), (136, 546), (143, 847), (144, 605), (144, 703), (146, 220), (147, 736), (152, 820), (154, 595), (160, 4), (162, 722), (173, 465), (177, 395), (177, 630), (179, 921), (180, 358), (202, 775), (206, 831), (210, 416), (210, 419), (211, 526), (222, 404), (224, 491), (225, 244), (228, 614), (229, 642), (237, 681), (240, 427), (246, 506), (249, 931), (265, 689), (270, 451), (281, 667), (282, 253), (292, 105), (292, 577), (324, 582), (331, 523), (399, 695), (462, 193), (507, 714), (566, 308), (651, 151)]
```

When interpreting this data, each tuple has 2 persons: the person in element 0 of each tuple is the infectious person, while the person in element 1 is the susceptible person that got infected. Therefore, looking at the first tuple, person 2 infected person 274 on day 9.

If one wants the number of transmissions, the sortedTransmissions() will return the transmissions, sorted in descending order, as well as the persons that got that many transmissions.

```python

>>>e = d.sortedTransmissions()
>>>print(e)
[(193.0, [28]), (180.0, [19]), (178.0, [15]), (176.0, [30]), (158.0, [65]), (154.0, [41]), (152.0, [6, 58]), (136.0, [43]), (129.0, [79]), (116.0, [56]), (103.0, [158]), (100.0, [53, 70]), (95.0, [69]), (94.0, [39]), (86.0, [102, 127]), (81.0, [109]), (77.0, [61]), (76.0, [131]), (67.0, [90]), (66.0, [148]), (57.0, [145]), (55.0, [83]), (53.0, [71]), (50.0, [17]), (49.0, [163]), (48.0, [27]), (47.0, [184, 185]), (46.0, [40, 161]), (45.0, [48]), (44.0, [8, 129, 142]), (43.0, [218]), (42.0, [0, 1]), (41.0, [110, 169]), (40.0, [7, 18]), (39.0, [2]), (38.0, [13]), (37.0, [3, 187]), (36.0, [16]), (35.0, [162]), (34.0, [10, 202]), (33.0, [20]), (32.0, [95]), (31.0, [29, 101]), (30.0, [11, 25, 213]), (29.0, [9, 23, 166]), (28.0, [103, 155, 172]), (27.0, [52]), (26.0, [4]), (25.0, [80, 87]), (24.0, [59, 159]), (23.0, [5, 32, 55, 126]), (22.0, [34, 188, 211, 251, 296]), (21.0, [21, 66, 96, 100, 235]), (20.0, [31, 37, 51, 147, 234]), (19.0, [36, 46, 119, 186, 203]), (18.0, [26, 47, 50, 60, 72, 82]), (17.0, [63, 84, 86, 89, 111, 264, 299]), (16.0, [12, 14, 22, 42, 45, 74, 106, 201, 317, 436]), (15.0, [49, 73, 75, 122, 140, 284, 314, 607]), (14.0, [62, 64, 85, 250, 339]), (13.0, [76, 91, 123, 133, 365]), (12.0, [38, 67, 77, 99, 116, 236, 287, 355, 721, 964]), (11.0, [35, 78, 107, 136, 154, 200, 239, 276, 305, 361]), (10.0, [54, 68, 94, 97, 112, 115, 121, 138, 167, 171, 196, 999]), (9.0, [57, 92, 104, 149, 156, 206, 219, 246, 315, 864]), (8.0, [88, 125, 130, 151, 152, 194, 277, 301, 348, 354, 378, 572]), (7.0, [24, 81, 93, 98, 117, 124, 128, 139, 141, 227, 295, 446]), (6.0, [113, 134, 137, 192, 238, 268, 379, 421, 430, 640, 658, 715, 811]), (5.0, [44, 191, 247, 254, 313, 328, 334, 343, 395, 404, 522, 605, 678, 966]), (4.0, [153, 157, 164, 165, 175, 181, 182, 190, 198, 199, 208, 220, 223, 245, 252, 258, 288, 345, 359, 370, 384, 429, 454, 488, 489, 490, 610, 871]), (3.0, [33, 105, 108, 135, 143, 144, 150, 168, 173, 176, 177, 180, 195, 226, 232, 265, 304, 306, 332, 337, 352, 362, 371, 381, 385, 387, 390, 394, 464, 483, 517, 533, 534, 697, 837, 939, 954]), (2.0, [118, 132, 178, 179, 183, 204, 207, 209, 210, 212, 224, 230, 237, 241, 242, 274, 278, 297, 318, 319, 325, 327, 335, 342, 363, 366, 374, 402, 411, 412, 413, 449, 450, 452, 461, 479, 492, 495, 501, 512, 514, 540, 547, 571, 582, 598, 600, 602, 604, 616, 619, 686, 691, 718, 809, 823, 835, 838, 842, 899, 916, 946]), (1.0, [114, 160, 189, 205, 222, 225, 231, 233, 240, 243, 244, 253, 256, 257, 259, 260, 263, 267, 271, 280, 281, 282, 285, 286, 290, 298, 308, 309, 311, 322, 329, 346, 347, 349, 351, 353, 373, 388, 391, 393, 409, 418, 420, 426, 428, 431, 432, 440, 447, 451, 455, 475, 476, 477, 482, 485, 493, 498, 500, 504, 511, 515, 519, 548, 554, 559, 562, 564, 568, 570, 576, 577, 581, 583, 584, 587, 593, 603, 629, 630, 643, 655, 662, 664, 675, 676, 684, 696, 704, 705, 710, 728, 754, 804, 819, 865, 880, 891, 913, 925, 931, 958, 961, 968]), (0.0, [120, 146, 170, 174, 193, 197, 214, 215, 216, 217, 221, 228, 229, 248, 249, 255, 261, 262, 266, 269, 270, 272, 273, 275, 279, 283, 289, 291, 292, 293, 294, 300, 302, 303, 307, 310, 312, 316, 320, 321, 323, 324, 326, 330, 331, 333, 336, 338, 340, 341, 344, 350, 356, 357, 358, 360, 364, 367, 368, 369, 372, 375, 376, 377, 380, 382, 383, 386, 389, 392, 396, 397, 398, 399, 400, 401, 403, 405, 406, 407, 408, 410, 414, 415, 416, 417, 419, 422, 423, 424, 425, 427, 433, 434, 435, 437, 438, 439, 441, 442, 443, 444, 445, 448, 453, 456, 457, 458, 459, 460, 462, 463, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 478, 480, 481, 484, 486, 487, 491, 494, 496, 497, 499, 502, 503, 505, 506, 507, 508, 509, 510, 513, 516, 518, 520, 521, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 535, 536, 537, 538, 539, 541, 542, 543, 544, 545, 546, 549, 550, 551, 552, 553, 555, 556, 557, 558, 560, 561, 563, 565, 566, 567, 569, 573, 574, 575, 578, 579, 580, 585, 586, 588, 589, 590, 591, 592, 594, 595, 596, 597, 599, 601, 606, 608, 609, 611, 612, 613, 614, 615, 617, 618, 620, 621, 622, 623, 624, 625, 626, 627, 628, 631, 632, 633, 634, 635, 636, 637, 638, 639, 641, 642, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 656, 657, 659, 660, 661, 663, 665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 677, 679, 680, 681, 682, 683, 685, 687, 688, 689, 690, 692, 693, 694, 695, 698, 699, 700, 701, 702, 703, 706, 707, 708, 709, 711, 712, 713, 714, 716, 717, 719, 720, 722, 723, 724, 725, 726, 727, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745, 746, 747, 748, 749, 750, 751, 752, 753, 755, 756, 757, 758, 759, 760, 761, 762, 763, 764, 765, 766, 767, 768, 769, 770, 771, 772, 773, 774, 775, 776, 777, 778, 779, 780, 781, 782, 783, 784, 785, 786, 787, 788, 789, 790, 791, 792, 793, 794, 795, 796, 797, 798, 799, 800, 801, 802, 803, 805, 806, 807, 808, 810, 812, 813, 814, 815, 816, 817, 818, 820, 821, 822, 824, 825, 826, 827, 828, 829, 830, 831, 832, 833, 834, 836, 839, 840, 841, 843, 844, 845, 846, 847, 848, 849, 850, 851, 852, 853, 854, 855, 856, 857, 858, 859, 860, 861, 862, 863, 866, 867, 868, 869, 870, 872, 873, 874, 875, 876, 877, 878, 879, 881, 882, 883, 884, 885, 886, 887, 888, 889, 890, 892, 893, 894, 895, 896, 897, 898, 900, 901, 902, 903, 904, 905, 906, 907, 908, 909, 910, 911, 912, 914, 915, 917, 918, 919, 920, 921, 922, 923, 924, 926, 927, 928, 929, 930, 932, 933, 934, 935, 936, 937, 938, 940, 941, 942, 943, 944, 945, 947, 948, 949, 950, 951, 952, 953, 955, 956, 957, 959, 960, 962, 963, 965, 967, 969, 970, 971, 972, 973, 974, 975, 976, 977, 978, 979, 980, 981, 982, 983, 984, 985, 986, 987, 988, 989, 990, 991, 992, 993, 994, 995, 996, 997, 998])]

```
In this list, there are tuples, where the first element is an number, and the second element is a list of integers. The 0th element of each tuple is the number of transmissions, and the list contains all the people who had that many transmissions. For example, when interpreting the above output, person 28 had 193 transmissions, person 19 had 180 transmissions, etc.# Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

E0: the number of exposed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

V0: the number of vaccinated individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.


rho: the probability that an individuals goes from E compartment to I compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

kappa: the probability that an individuals goes from R compartment to S compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

mu: the probability that an individual goes from I compartment to D compartment, given that the individual didn't go from I to R in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

eta: the probability that an individual goes from S compartment to V compartment, given that the individual didn't go from S to E in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

timeDelay: float
    The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. If not float, will throw ```python NotFloatException```. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

E0: the number of exposed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

V0: the number of vaccinated individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.


rho: the probability that an individuals goes from E compartment to I compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

mu: the probability that an individual goes from I compartment to D compartment, given that the individual didn't go from I to R in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

eta: the probability that an individual goes from S compartment to V compartment, given that the individual didn't go from S to E in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

timeDelay: float
    The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. If not float, will throw ```python NotFloatException```. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

E0: the number of exposed individuals at the start of the simulation. Must be an itn. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

rho: the probability that an individuals goes from E compartment to I compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

kappa: the probability that an individual goes from R compartment to S compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

mu: the probability that an individual goes from I compartment to D compartment, given that the person didn't go from I to R in the same timestep. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

V0: the number of vaccinated individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

eta: the probability that an individual goes from S compartment to V compartment, given that the individual didn't go from S to I in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

mu: the probability that an individual goes from I compartment to D compartment, given that the individual didn't go from I to R in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

timeDelay: float
    The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. If not float, will throw ```python NotFloatException```. # Note
In all compartmental models with an E compartment, it is assumed that those in E cannot propogate the disease.

# Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

E0: the number of exposed individuals at the start of the simulation. Must be an itn. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

rho: the probability that an individuals goes from E compartment to I compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

E0: the number of exposed individuals at the start of the simulation. Must be an itn. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

rho: the probability that an individuals goes from E compartment to I compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

kappa: the probability that an individual goes from R compartment to S compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

E0: the number of exposed individuals at the start of the simulation. Must be an itn. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

V0: the number of vaccinated individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

rho: the probability that an individuals goes from E compartment to I compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

eta: the probability that an individual goes from S compartment to V compartment, given that the person didn't go from S to E. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

timeDelay: float
    The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. If not float, will throw ```python NotFloatException```. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

E0: the number of exposed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

rho: the probability that an individuals goes from E compartment to I compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

mu: the probability that an individual goes from I compartment to D compartment, given that the individual didn't go from I to R in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. # Overview

In the ICU model, there are many compartments. There are the standard compartments that were seen before, namely S, E, I, R and V. However, there is now two new compartments: L and ICU. The L compartment, or "Lag compartment", represents those who, after leaving E compartment, will be destined to eventually go to the ICU, but aren't there yet. It is assumed those in L can still propogate the infection. The ICU compartment represents those who are hospitalized and will either go to the D compartment or R compartment. 

# Parameters

        S0 : int
            The starting number of susceptible individuals in the simulation.
        
        E0: int
            The starting number of exposed individuals in the simulation.
        
        I0: int
            The starting number of infected individuals in the simulation.

        R0: int
            The starting number of recovered individuals in the simulation.
        
        V0: int
            The starting number of vaccinated individuals in the simulation.
        
        rho: float
            The probability of an individual leaving the E compartment.
        
        ioda: float
            The probability that, given an individual is leaving the E compartment, he goes to L compartment. The probability of that person going to I compartment is (1-ioda).
        
        gamma: float
            The probability of a person in I compartment going to the R compartment
        
        mu: float
            The probability of going from I to D, given that the person didn't go from I to R.
        
        phi: float
            The probability of going from L compartment to ICU compartment.
        
        chi: float
            The probability of going from ICU compartment to R compartment.
        
        omega: float
            The probability of going from ICU compartment to D compartment, given that individual didn't go from ICU compartment to R compartment.
        
        kappa: float
            The probability of going from R compartment to S compartment.
        
        eta: float 
            The probability of going from S compartment to V compartment, given that the individual didn't go from S compartment to E compartment. 
        
        timeDelay: float
            The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. 
        # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int.
Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int.
Throws ```python NotIntException``` otherwise.

gamma: the probability that an individuals goes from I compartment to S compartment.
Throws ```python ProbabilityException``` otherwise.# Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

V0: the number of vaccinated individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

eta: the probability that an individual goes from S compartment to V compartment, given that the individual didn't go from S to I in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

timeDelay: float
    The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. If not float, will throw ```python NotFloatException```. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

V0: the number of vaccinated individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

kappa: the probability that an individual goes from R compartment to S compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

eta: the probability that an individual goes from S compartment to V compartment, given that the individual didn't go from S to I in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

mu: the probability that an individual goes from I compartment to D compartment, given that the individual didn't go from I to R in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

timeDelay: float
    The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. If not float, will throw ```python NotFloatException```. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

V0: the number of vaccinated individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

kappa: the probability that an individual goes from R compartment to S compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

eta: the probability that an individual goes from S compartment to V compartment, given that the individual didn't go from S to I in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

timeDelay: float
    The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. If not float, will throw ```python NotFloatException```. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

E0: the number of exposed individuals at the start of the simulation. Must be an itn. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

V0: the number of vaccinated individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

rho: the probability that an individuals goes from E compartment to I compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

kappa: the probability that an individual goes from R compartment to S compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

eta: the probability that an individual goes from S compartment to V compartment, given that the person didn't go from S to E. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

timeDelay: float
    The number of days that vaccine rollout is delayed. If negative or 0, then there is no delay in vaccine rollout. Default value is -1. If not float, will throw ```python NotFloatException```. 

# Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

kappa: the probability that an individual goes from R compartment to S compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

mu: the probability that an individual goes from I compartment to D compartment, given that the individual didn't go from I to R in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. # Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

kappa: the probability that an individual goes from R compartment to S compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 
# Parameters

S0 : the number of susceptible individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

I0: the number of infected individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

R0: the number of removed individuals at the start of the simulation. Must be an int. Throws ```python NotIntException``` otherwise.

gamma: the probability that an individuals goes from I compartment to R compartment. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. 

mu: the probability that an individual goes from I compartment to D compartment, given that the individual didn't go from I to R in that same time step. Must be a value that belongs to [0,1]. Throws ```python ProbabilityException``` or ```python NotFloatException``` otherwise. # Examples

To run the simulation with a certain number of days, you can run code like the following:

```python

from Eir import SEIR

sim = SEIR(S0=9999999, E0=10000, I0=1, R0=0, beta=1.5, rho=.25, gamma=.15)
df = sim.run(31, .1, plot=False)
```

run() takes the number of simulated as first parameter, and the differential, or step, used for the Euler approximation of the ODE. 

If a plot is wanted, simply run the following:

```python

from Eir import SEIR

sim = SEIR(S0=9999999, E0=10000, I0=1, R0=0, beta=1.5, rho=.25, gamma=.15)
df, fig = sim.run(31, .1)

```

This will display a plot of all variables, which can be further customized using the default boolean parameters that represent all of the variables within the model.# Examples

If one wanted to get a plot of a particular SIS model, the run function would be used. run() takes the number of simulated as first parameter, and the differential, or step, used for the Euler approximation of the ODE. A way to run it would be as follows:

```python 

from Eir import SIS

sim = SIS(S0=9999999, I0=1, beta=1.5, gamma=.15)
df = sim.run(31, .1, plot=False)

```

When printing the dataframe, the output will be the following:

```
>>> print(df)
     Days   Susceptible      Infected
0     0.0  9.999999e+06  1.000000e+00
1     0.1  9.999999e+06  1.135000e+00
2     0.2  9.999999e+06  1.288225e+00
3     0.3  9.999999e+06  1.462135e+00
4     0.4  9.999998e+06  1.659524e+00
..    ...           ...           ...
306  30.6  1.000000e+06  9.000000e+06
307  30.7  1.000000e+06  9.000000e+06
308  30.8  1.000000e+06  9.000000e+06
309  30.9  1.000000e+06  9.000000e+06
310  31.0  1.000000e+06  9.000000e+06

[311 rows x 3 columns]
>>> 
```

However, if plot=True in the run method, then a pyplot figure object will be returned with the dataframe in a tuple, with element 0 being the dataframe and element 1 being the figure object. A plot will also be displayed that represent all of the variables within the model.

# Examples

To run the simulation with a certain number of days, you can run code like the following:

```python

from Eir import SIRV

sim = SIRV(S0=9999999, I0=1, R0=0, V0=0, beta=1.5, gamma=.15, eta=.01)
df = sim.run(31, .1, plot=False)
```

run() takes the number of simulated as first parameter, and the differential, or step, used for the Euler approximation of the ODE.

If a plot is wanted, simply run the following:

```python

from Eir import SIRV

sim = SIRV(S0=9999999, I0=1, R0=0, V0=0, beta=1.5, gamma=.15, eta=.01)
df, fig = sim.run(31, .1)

```

This will display a plot of all variables, which can be further customized using the default boolean parameters that represent all of the variables within the model.# Examples

To run the simulation with a certain number of days, you can run code like the following:

```python

from Eir import SIR

sim = SIR(S0=9999999, I0=1, R0=0, beta=1.5, gamma=.15)
df = sim.run(31, .1, plot=False)
```

run() takes the number of simulated as first parameter, and the differential, or step, used for the Euler approximation of the ODE.

If a plot is wanted, simply run the following:

```python

from Eir import SIR

sim = SIR(S0=9999999, I0=1, R0=0, beta=1.5, gamma=.15)
df, fig = sim.run(31, .1)

```

This will display a plot of all variables, which can be further customized using the default boolean parameters that represent all of the variables within the model.# Examples

To run the simulation with a certain number of days, you can run code like the following:

```python

from Eir import SIRS

sim = SIRS(S0=9999999, I0=1, R0=0, beta=1.5, gamma=.15, kappa=.05)
df = sim.run(31, .1, plot=False)
```

run() takes the number of simulated as first parameter, and the differential, or step, used for the Euler approximation of the ODE.

If a plot is wanted, simply run the following:

```python

from Eir import SIRS

sim = SIRS(S0=9999999, I0=1, R0=0, beta=1.5, gamma=.15, kappa=.15)
df, fig = sim.run(31, .1)

```

This will display a plot of all variables, which can be further customized using the default boolean parameters that represent all of the variables within the model.# Examples

To run the simulation with a certain number of days, you can run code like the following:

```python

from Eir import SIRD

sim = SIRD(S0=9999999, I0=1, R0=0, beta=1.5, gamma=.15, omega=.01)
df = sim.run(31, .1, plot=False)
```

run() takes the number of simulated as first parameter, and the differential, or step, used for the Euler approximation of the ODE.

If a plot is wanted, simply run the following:

```python

from Eir import SIRD

sim = SIRD(S0=9999999, I0=1, R0=0, beta=1.5, gamma=.15, omega=.01)
df, fig = sim.run(31, .1)

```

This will display a plot of all variables, which can be further customized using the default boolean parameters that represent all of the variables within the model.In order to run each test, after forking the repository, go to each of the python file for each of the tests and run the following in the shell:

```shell

python3 <NAME_OF_FILE>
```

For example, if I wanted to run the test for the HubSEIRSV model, the following command would be run:

```shell

python3 test_HubSEIR.py

```
Additionally, the shell files in each subdirectory(Deterministic, Hub, Periodic). For example, to run all of the periodic model tests from the tests/Periodic/ directory:

```shell

./test_Periodic.sh
```

In order to test every test for every model type, run the following from the tests/ directory:
```shell
./all_tests.sh
```
# Structure

Because Eir is an epidemiological simulation/modeling package, the outputs are stochastic. Therefore, to control for this, the random number generator in each of the tests are seeded. Therefore, the RNG will produce the same numbers, which should therefore lead to reproducible results. The CSVs contain the outputs that each test file should output with the seed they have. 

