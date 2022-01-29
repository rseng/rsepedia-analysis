
# `caffsim` R package: Simulation of Plasma Caffeine Concentrations by Using Population Pharmacokinetic Model

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.842649.svg)](https://doi.org/10.5281/zenodo.842649)
![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/caffsim)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/caffsim)](https://cran.r-project.org/package=caffsim)

> Simulate plasma caffeine concentrations using population
> pharmacokinetic model described in Lee, Kim, Perera, McLachlan and Bae
> (2015) <doi:10.1007/s00431-015-2581-x> and the package was published
> <doi:10.12793/tcp.2017.25.3.141>.

![](inst/doc/cover.png)

  - Github: <https://github.com/asancpt/caffsim>
  - Package vignettes and references by `pkgdown`:
    <http://asancpt.github.io/caffsim>

## Installation

``` r
install.pacakges("devtools")
devtools::install_github("asancpt/caffsim")

# Simply create single dose dataset
caffsim::caffPkparam(Weight = 20, Dose = 200, N = 20) 

# Simply create multiple dose dataset
caffsim::caffPkparamMulti(Weight = 20, Dose = 200, N = 20, Tau = 12) 
```

## Single dose

### Create a PK dataset for caffeine single dose

``` r
library(caffsim)
MyDataset <- caffPkparam(Weight = 20, Dose = 200, N = 20)
head(MyDataset)
```

<div class="kable-table">

| subjid |      Tmax |      Cmax |       AUC | Half\_life |       CL |        V |         Ka |        Ke |
| -----: | --------: | --------: | --------: | ---------: | -------: | -------: | ---------: | --------: |
|      1 | 0.9858664 | 15.006075 | 115.44184 |   4.594639 | 1.732474 | 11.48643 |  3.2719387 | 0.1508280 |
|      2 | 3.8148250 |  7.699503 | 107.02942 |   6.354802 | 1.868645 | 17.13545 |  0.5169846 | 0.1090514 |
|      3 | 0.8826191 | 12.147005 | 107.60372 |   5.491894 | 1.858672 | 14.72962 |  4.0586205 | 0.1261860 |
|      4 | 0.2567103 | 13.950749 | 139.11091 |   6.730028 | 1.437702 | 13.96215 | 20.7765914 | 0.1029713 |
|      5 | 1.1563495 |  8.302007 |  65.44492 |   4.587332 | 3.056005 | 20.22931 |  2.6177408 | 0.1510682 |
|      6 | 1.2026953 |  9.445008 |  55.68603 |   3.130858 | 3.591565 | 16.22609 |  2.0869246 | 0.2213451 |

</div>

### Create a dataset for concentration-time curve

``` r
MyConcTime <- caffConcTime(Weight = 20, Dose = 200, N = 20)
head(MyConcTime)
```

<div class="kable-table">

| Subject | Time |      Conc |
| ------: | ---: | --------: |
|       1 |  0.0 |  0.000000 |
|       1 |  0.1 |  4.129452 |
|       1 |  0.2 |  7.174502 |
|       1 |  0.3 |  9.410497 |
|       1 |  0.4 | 11.042950 |
|       1 |  0.5 | 12.225252 |

</div>

### Create a concentration-time curve

``` r
caffPlot(MyConcTime)
```

![](assets/figures/MyPlotMyConcTime-1.png)<!-- -->

### Create plots for publication (according to the amount of caffeine)

  - `cowplot` package is required

<!-- end list -->

``` r
#install.packages("cowplot") # if you don't have it
library(cowplot)

MyPlotPub <- lapply(
  c(seq(100, 800, by = 100)), 
  function(x) caffPlotMulti(caffConcTime(20, x, 20)) + 
    theme(legend.position="none") + 
    labs(title = paste0("Single Dose ", x, "mg")))

plot_grid(MyPlotPub[[1]], MyPlotPub[[2]],
          MyPlotPub[[3]], MyPlotPub[[4]],
          MyPlotPub[[5]], MyPlotPub[[6]],
          MyPlotPub[[7]], MyPlotPub[[8]],
          labels=LETTERS[1:8], ncol = 2, nrow = 4)
```

![](assets/figures/MyPlotPub-1.png)<!-- -->

## Multiple dose

### Create a PK dataset for caffeine multiple doses

``` r
MyDatasetMulti <- caffPkparamMulti(Weight = 20, Dose = 200, N = 20, Tau = 12)
head(MyDatasetMulti)
```

<div class="kable-table">

| subjid |     TmaxS |    CmaxS |      AUCS |       AI |     Aavss |     Cavss |   Cmaxss |   Cminss |
| -----: | --------: | -------: | --------: | -------: | --------: | --------: | -------: | -------: |
|      1 | 0.9548818 | 14.35410 |  78.14388 | 1.068908 |  72.79768 |  6.511990 | 19.08368 | 1.230238 |
|      2 | 1.1467783 | 11.52173 |  89.29910 | 1.187191 | 108.04592 |  7.441592 | 16.31939 | 2.573172 |
|      3 | 0.6768031 | 16.69770 | 136.35979 | 1.250128 | 124.04012 | 11.363316 | 22.85721 | 4.573316 |
|      4 | 0.7753311 | 14.44837 | 119.92576 | 1.251847 | 124.46355 |  9.993813 | 20.06162 | 4.036002 |
|      5 | 0.2400863 | 13.36207 |  65.42138 | 1.081913 |  77.33312 |  5.451782 | 15.22267 | 1.152520 |
|      6 | 0.6535510 | 10.73536 |  76.56109 | 1.183885 | 107.17374 |  6.380091 | 14.06610 | 2.184791 |

</div>

### Create a dataset for concentration-time curve

``` r
MyConcTimeMulti <- caffConcTimeMulti(Weight = 20, Dose = 200, N = 20, Tau = 12, Repeat = 10)
head(MyConcTimeMulti)
```

<div class="kable-table">

| Subject | Time |      Conc |
| ------: | ---: | --------: |
|       1 |  0.0 |  0.000000 |
|       1 |  0.1 |  4.915431 |
|       1 |  0.2 |  8.455652 |
|       1 |  0.3 | 10.989413 |
|       1 |  0.4 | 12.786772 |
|       1 |  0.5 | 14.045507 |

</div>

### Create a concentration-time curve

``` r
caffPlotMulti(MyConcTimeMulti)
```

![](assets/figures/MyPlotMultiMyConcTimeMulti-1.png)<!-- -->

### Create plots for publication (according to dosing interval)

  - `cowplot` package is required

<!-- end list -->

``` r
#install.packages("cowplot") # if you don't have it
library(cowplot)

MyPlotMultiPub <- lapply(
  c(seq(4, 32, by = 4)), 
  function(x) caffPlotMulti(caffConcTimeMulti(20, 250, 20, x, 15)) + 
    theme(legend.position="none") + 
    labs(title = paste0("q", x, "hr" )))

plot_grid(MyPlotMultiPub[[1]], MyPlotMultiPub[[2]],
          MyPlotMultiPub[[3]], MyPlotMultiPub[[4]],
          MyPlotMultiPub[[5]], MyPlotMultiPub[[6]],
          MyPlotMultiPub[[7]], MyPlotMultiPub[[8]],
          labels=LETTERS[1:8], ncol = 2, nrow = 4)
```

![](assets/figures/MyPlotMultiPub-1.png)<!-- -->

## Interactive shiny app

``` r
caffShiny()
```
caffsum 0.2.4
===========

* Update `shiny` app

caffsim 0.2.3
===========

* Update manual


caffsim 0.2.2
===========

* Introduce a `shiny` example - `Caffeine Concentration Predictor`
* Add several functions

caffsim 0.2.0
===========

* Initial CRAN release


caffsim 0.1.0
===========

* Initial private beta release!
### Parameters

- CMAX: maximum concentration, Cmax
- CMAXD: CMAX / Dose, Cmax / Dose
- TMAX: time of maximum concentration, Tmax
- TLAG: time until first nonzero concentration, for extravascular administration only
- CLST: last positive concentration observed, Clast
- CLSTP: last positive concentration predicted, Clast_pred
- TLST: time of last positive concentration, Tlast
- LAMZHL: half-life by lambda z, ln(2)/LAMZ
- LAMZ: lambda_z negative of best fit terminal slope
- LAMZLL: earlist time for LAMZ
- LAMZUL: last time for LAMZ
- LAMZNPT: number of points for LAMZ
- CORRXY: correlaton of log(concentration) and time
- R2: R-squared
- R2ADJ: R-squared adjusted
- C0: back extrapolated concentration at time 0, for bolus intravascular administration only
- AUCLST: AUC from 0 to TLST
- AUCALL: AUC using all the given points, including trailing zero concentrations
- AUCIFO: AUC infinity observed
- AUCIFOD: AUCIFO / Dose
- AUCIFP: AUC infinity predicted using CLSTP instead of CLST
- AUCIFPD: AUCIFP / Dose
- AUCPEO: AUC % extrapolation observed
- AUCPEP: AUC % extrapolated for AUCIFP
- AUCPBEO: AUC % back extrapolation observed, for bolus IV administration only
- AUCPBEP: AUC % back extrapolation predicted with AUCIFP, for bolus IV administration only
- AUMCLST: AUMC to the TLST
- AUMCIFO: AUMC infinity observed using CLST
- AUMCIFP: AUMC infinity determiend by CLSTP
- AUMCPEO: AUMC % extrapolated observed
- AUMCPEP: AUMC % extrapolated predicted
- MRTIVLST: mean residence time(MRT) to TLST, for intravascular administration
- MRTIVIFO: mean residence time(MRT) infinity using CLST, for intravascular administration
- MRTIVIFP: mean residence time(MRT) infinity ucinsg CLSTP, for intravascular administration
- MRTEVLST: mean residence time(MRT) to TLST, for extravascular administration
- MRTEVIFO: mean residence time(MRT) infinity using CLST, for extravascular administration
- MRTEVIFP: mean residence time(MRT) infinity ucinsg CLSTP, for extravascular administration
- VZO: volume of distribution determined by LAMZ and AUCIFO, for intravascular administration
- VZP: volume of distribution determined by LAMZ and AUCIFP, for intravascular administration
- VZFO: VZO for extravascular administration, VZO/F, F is bioavailability
- VZFP: VZP for extravascular administration, VZP/F, F is bioavailability
- CLO: clearance using AUCIFO, for intravascular administration
- CLP: clearance using AUCIFP, for intravascular administration
- CLFO: CLO for extravascular administration, CLO/F, F is bioavailability
- CLFP: CLP for extravascular administration, CLP/F, F is bioavailability
- VSSO: volume of distribution at steady state using CLST, for intravascular administration only
- VSSP: volume of distribution at stead state using CLSTP, for intravascular administration only
### Contact

- `Caffeine Concentration Predictor` Shiny developer: Sungpil Han <shan@acp.kr> / <https://shanmdphd.github.io>
- Acknowledgements : The main idea and intellectual resources were mainly provided by Professor Kyun-Seop Bae <k@acr.kr>.
- Copyright: 2016, Sungpil Han
- License: GPL-3
## Caffeine Concentration Predictor

- `Caffeine Concentration Predictor` <https://asan.shinyapps.io/caff>
- You can run it locally by entering `caffsim::shinyCaff()` in a R console.
- `Caffeine Concentration Predictor` is open to everyone. We are happy to take your input. Please fork the repo, modify the codes and submit a pull request. <https://github.com/asancpt/caff>

### Background

The pharmacokinetic parameters from the paper were derived and used in the app as follows:

$$ 
\begin{split}
\begin{bmatrix}
     \eta_1 \newline
     \eta_2 \newline
     \eta_3
\end{bmatrix}
& \sim MVN \bigg(
\begin{bmatrix}
     0 \newline
     0 \newline
     0
\end{bmatrix}
, 
\begin{bmatrix}
     0.1599 & 6.095 \cdot 10^{-2} & 9.650 \cdot 10^{-2} \newline
     6.095 \cdot 10^{-2} & 4.746 \cdot 10^{-2} & 1.359 \cdot 10^{-2} \newline
     9.650 \cdot 10^{-2} & 1.359 \cdot 10^{-2} & 1.004
\end{bmatrix}
\bigg) \newline
\newline
CL\ (mg/L) & = 0.09792 \cdot W \cdot e^{\eta1} \newline
V\ (L) & = 0.7219 \cdot W \cdot e^{\eta2} \newline
k_a\ (1/hr) & = 4.268 \cdot e^{\eta3} \newline
\newline
k\ (1/hr) & = \frac{CL}{V} \newline
t_{1/2}\ (hr) & = \frac{0.693}{k} \newline
t_{max}\ (hr) & = \frac{ln(k_a) - ln(k)}{k_a - k} \newline
C_{max}\ (mg/L) & = \frac{Dose}{V} \cdot \frac{k_a}{k_a - k} \cdot (e^{-k \cdot  t_{max}} - e^{-k_a \cdot t_{max}}) \newline
AUC\ (mg \cdot hr / L)  & = \frac{Dose}{CL} \newline
\newline
C_{av,ss} & = \frac{Dose}{CL \cdot \tau} \newline 
AI & = \frac{1}{1-e^{-k_e \cdot \tau}} \newline
\end{split}
$$
(Abbreviation: $AI$, accumulation index; $AUC$, area under the plasma drug concentration-time curve; $CL$, total clearance of drug from plasma; $C_{av,ss}$, average drug concentration in plasma during a dosing interval at steady state on administering a fixed dose at equal dosing intervals; $C_{max}$, highest drug concentration observed in plasma; $MVN$, multivariate normal distribution; $V$, Volume of distribution (apparent) based on drug concentration in plasma; $W$, body weight (kg); $\eta$, interindividual random variability parameter; $k$, elimination rate constant;  $k_a$, absorption rate constant; $\tau$, dosing interval; $t_{1/2}$, elimination half-life)

### R Packages
- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.
- Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2016). shiny: Web Application Framework for R. R package version 0.14.2. https://CRAN.R-project.org/package=shiny
- JJ Allaire, Jeffrey Horner, Vicent Marti and Natacha Porte (2015). markdown: 'Markdown' Rendering for R. R package version 0.7.7. https://CRAN.R-project.org/package=markdown
- Hadley Wickham and Romain Francois (2016). dplyr: A Grammar of Data Manipulation. R package version 0.5.0. https://CRAN.R-project.org/package=dplyr


### Reference

This work is solely dependent on the paper published in Eur J Pediatr in 2015. The package is published in Translational Clinical Pharmacology in 2017.

- "Prediction of plasma caffeine concentrations in young adolescents following ingestion of caffeinated energy drinks: a Monte Carlo simulation." Eur J Pediatr. 2015 Dec;174(12):1671-8. doi: 10.1007/s00431-015-2581-x <https://www.ncbi.nlm.nih.gov/pubmed/26113286>
-  <doi:10.12793/tcp.2017.25.3.141>. 
"Caffsim: simulation of plasma caffeine concentrations implemented as an R package and Web-applications" Transl Clin Pharmacol. 2017 Sep;25(3):141-146. doi: 10.12793/tcp.2017.25.3.141 <https://doi.org/10.12793/tcp.2017.25.3.141>
- "Clinical pharmacokinetics and pharmacodynamics: concepts and applications, 4th edition" Lippincott Williams & Wilkins. 2011. ISBN 978-0-7817-5009-7

### Caffeine contents

<div align=center><img src=http://graphs.net/wp-content/uploads/2013/01/Caffeine-Content-in-Energy-Drinks.jpg width = 450 /></div>
**Reference range**

- Below 10 mg/L: generally considered safe (Green horizontal line)
- Over 40 mg/L: several fatalities (Blue horizontal line)
- Over 80 mg/L: fatal caffeine poisoning (Red horizontal line)
- Reference: de Wijkerslooth LR et al.(2008), Seifert et al.(2013), Banerjee et al. (2014), Cannon et al. (2001)
### Apps by AMC CPT

#### Vancomycin TDM KR
<img src=http://i.imgur.com/M8R3Me8.png width = 500 />

- <https://asan.shinyapps.io/vtdm>
- 2016-12 

#### Caffeine Concentration Predictor
<img src=http://i.imgur.com/RYJHxNq.png width = 500 />

- <https://asan.shinyapps.io/caff>
- 2016-12

#### Online NonCompart
<img src=http://i.imgur.com/k6VqHp2.png width = 500 />

- <https://asan.shinyapps.io/noncompart>
- 2016-11 
