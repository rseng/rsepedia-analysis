---
title: 'ACEInhibPKPD: An open-source MATLAB app for a pharmacokinetic/pharmacodynamic model of ACE inhibition'
tags:
  - pharmacology
  - biomedical engineering
  - mathematical model
  - MATLAB
authors:
 - name: Ashlee N. Ford Versypt
   orcid: 0000-0001-9059-5703
   affiliation: 1
 - name: Grace K. Harrell
   affiliation: 1
 - name: Alexandra N. McPeak
   affilitation: 1
   
affiliations:
 - name: School of Chemical Engineering, Oklahoma State University
   index: 1
date: 15 June 2017
bibliography: paper.bib
---

# Summary

This software features an MATLAB app install file and the source code for an interactive computer simulation packaged as a MATLAB app that can be used to design the best dosage of pharmaceuticals (ACE inhibitors) for reducing high blood pressure. The app uses concepts from chemical engineering and tools from mathematics, pharmacology, and computational science to describe and solve the dynamics of the chemical reactions in the human body involved in the absorption, metabolism, and excretion of the drug and how the blood pressure-regulating hormone angiotensin II is affected by the drug concentration as a function of time. The scienfic significance is that the model solved in the software fills the gap in pharmacokinetic/pharmacodynamic models for ACE inhibitors in impaired renal function, given the prevalence of using this class of pharmaceuticals for treating chronic kidney disease in addition to high blood pressure. The model used in this software is described in detail in [@FordVersypt2017]. Use of the software in the context of engineering education is described in [@Harrell2017].

![GUI screenshot](thumbnail.png)

# References
# ACEInhibPKPD
PKPD model for ACE inhibition with benazepril or cilazapril; educational drug dosing app for treatment of high blood pressure

Note: SB17CHE.mlappinstall is an updated version that includes 5 drugs: benazepril, cilazapril, lisinopril, delapril, and ramipril.

Code v.1.1 JOSS reviewed
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.885362.svg)](https://doi.org/10.5281/zenodo.885362)

Code v.2.1 Expanded GUI
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1318682.svg)](https://doi.org/10.5281/zenodo.1318682)

Paper:
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00340/status.svg)](https://doi.org/10.21105/joss.00340)

![GUI screenshot](thumbnail.png)

## Overview
This app is an interactive computer simulation that can be used to design the best dosage of two pharmaceuticals for reducing high blood pressure. The app uses concepts from chemical engineering and tools from mathematics, pharmacology, and computational science to describe and solve the the dynamics of the chemical reactions in the human body involved in the absorption, metabolism, and excretion of the drug and how the blood pressure-regulating homone angiotensin II (Ang II) is affected by the drug concentration as a function of time. The model solved in the software is a pharmacokinetic/pharmacodynamic model for ACE inhibitors that reduce the levels of Ang II in the body to lower blood pressure. A novel feature of the model is the incoporation of normal and impaired renal function patient data and scenarios. This is imporant given the prevalence of using this class of pharmaceuticals for treating chronic kidney disease in addition to high blood pressure.

## Target Audience
This software is intended to be used by engineering educators working with college freshmen or high school students and by researchers investigating pharmaceuticals to reduce Ang II with ACE inhibitors in hypertension and chronic kidney disease. 

## PKPD ACE Inhibition Model
### Authors
Ashlee N. Ford Versypt, Grace K. Harrell, and Alexandra N. McPeak, 
School of Chemical Engineering,
Oklahoma State University.
Corresponding author: A. N. Ford Versypt, ashleefv@okstate.edu

v.2.1 co-author Carley V. Eastep, School of Chemical Engineering,
Oklahoma State University.

### Software publication
A. N. Ford Versypt, G. K. Harrell, A. N. McPeak, (2017). ACEInhibPKPD: An open-source MATLAB app for a pharmacokinetic/pharmacodynamic model of ACE inhibition. Journal of Open Source Software, 2(17), 340, https://doi.org/10.21105/joss.00340

### Related publication for model details
A. N. Ford Versypt, G. K. Harrell, A. N. McPeak, A pharmacokinetic/pharmacodynamic model of ACE inhibition of the renin-angiotensin system for normal and impaired renal function, Computers & Chemical Engineering, 104 (2017) 311–322. https://doi.org/10.1016/j.compchemeng.2017.03.027

### Related publication for educational use details
G. K. Harrell, A. N. McPeak, and A. N. Ford Versypt, A pharmacokinetic simulation-based module to introduce mass balances and chemical engineering design concepts to engineering freshmen, Proceedings of the ASEE Annual Conference, Columbus, OH, 2017. https://peer.asee.org/27493

### Data sources for ACE inhibitor drugs
* H. Shionoiri, S. Ueda, K. Minamisawa, Pharmacokinetics and Pharmacodynamics 
   of Benazepril in Hypertensive Patients with Normal and Impaired Renal-
   Function, Journal of Cardiovascular Pharmacology 20 (3) (1992) 348–57.

* H. Shionoiri, E. Gotoh, N. Takagi, Antihypertensive Effects and 
   Pharmacokinetics of Single and Consecutive Doses of Cilazapril in 
   Hypertensive Patients with Normal and Impaired Renal-Function, Journal of 
   Cardiovascular Pharmacology 11 (2) (1988) 242–9.

for v.2.1
* M. Neubeck, D. Fliser, M. Pritsch, K. Weißer, M. Fliser, J.
    Nussberger, E. Ritz, E. Mutschler, Pharmacokinetics and
   pharmacodynamics of lisinopril in advanced renal failure, European
   Journal of Clinical Pharmacology 46 (6) (1994) 537-43.
*  H.Shionoiri, G. Yasuda, A. Ikeda, T. Ohta, E. Miyajima, and
   Y. Kaneko, Pharmacokinetics and depressor effect of delapril
   in patients with essential hypertension, Clinical Pharmacology & 
   Therapeutics 41 (1) (1987) 74–79.
* ER. Debusmann, J. Ocon-Pujadas, W. Lahn, et al, Influence of renal
    function on the pharmacokinetics of ramipril (HOE 498), American Journal Cardiology  59 (1987) D70-D78. 
*    J. Ocon-Pujadas, E. R. Debusmann, F. Jane, W. Lahn, R. Irmisch, J.
    Mora, and H. Grotsch, Pharmacodynamic effects of a single 10-mg dose
    of the angiotensin converting enzyme-inhibitor ramipril in patients
    with impaired renal-function, Journal of Cardiovascular Pharmacology
    13 (3) (1989) S45-S48.

### Installation file

The app install file ACEInhibPKPD.mlappinstall works much like a zip file with the added features of automated extraction in MATLAB and creation of a menu button to run the GUI. Click the menu button to open the GUI. Alternatively, the code may be opened directly in MATLAB. Descriptions of the main files and dependent & supplemental files are described in the following sections.

Note: SB17CHE.mlappinstall is an updated version that includes 5 drugs with normal and/or impaired kidney functions.

### Main files

* GUI_ACE_Inhib_PKPD.m & .fig.
   Runs the model interactively with a graphical user interface. Designed for 
   educational use for students in grades 6-12 and undergraduate students.
 
OR
  
* run_PKPD_without_GUI.m.
   Runs the model manually either in the Editor or by passing arguments in the
   command line. Can run more cases that allowable in the GUI and can specify
   plot_mode to show more plots than available through GUI or to produce and 
   save publication-quality plots.

### Dependent & supplemental files

* PKPD_ACE_Inhibition_AngII.m.
   The central file that uses case-specific parameters and dosing information 
   as input to the PKPD model of ACE inhibitor dose impact on AngII plasma 
   concentration. Output includes the vectors of time and the concentrations of 
   the drug, diacid form, AngI, AngII, and Renin and the percent of inhibition 
   (not in this order). Note: do not run this file directly. Instead let it be called by either 
   call_PKPD_model_scalar.m or call_PKPD_model_vector_param_est.m
   (described subsequently) that correctly define all the input values.

* setup_mat_for_PK_params.m.
   All PK parameters for the drug/renalfunction combinations are edited and 
   written to .mat files for use in other routines.

* setup_mat_for_data_estimated_params.m.
   All data for parameter estimation and the estimated parameters are edited and
   written to .mat files for use in other routines.
   
* run_params.m.
   Sets the input values that GUI_ACE_Inhib_PKPD.m sets interactively. It is
   used as a default input through a .mat file for run_PKPD_without_GUI.m and 
   param_estimation.m.

* call_PKPD_model_vector_param_est.m.
   Sets up parameters once before a loop through tfinal_dosing vector calling
   PKPD_ACE_Inhibition_AngII.m at each loop iteration. Output is the vector of  
   Y values used for parameter estimation comparison to data.
   
* call_PKPD_model_scalar.m.
   Used by GUI_ACE_Inhib_PKPD.m, run_PKPD_without_GUI.m, and param_estimation.m. 
   Sets up parameters before calling PKPD_ACE_Inhibition_AngII.m for a scalar 
   value of tfinal_dosing. Output is all of the output from 
   PKPD_ACE_Inhibition_AngII.m with options for plotting and saving figures
   through plot_mode string.
   
* .mat files. 
   Needed to run the programs to pass data, parameters, and calculated values.
   
* /pub_plots/.
   Subdirectory where plots are stored when plot_mode = 'pub_plots'
   
* export_fig.m
   MATLAB package to nicely export figures.
   Download: https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
   Tips for usage: https://github.com/altmany/export_fig/blob/master/README.md
   
## Example Usage
Download and run the app install file ACEInhibPKPD.mlappinstall. Begin using the GUI by entering dose size and frequency values and select the ACE inhibitor drug and the patient kidney function from the drop down menus. Run the simulation to generate the plots. Use the plots to analyze the behavior of the tested dosing schedule. Determine what dose size and frequency will be the most effective in lowering the blood pressure by lowering the angiotensin II hormone concentration of the individual to a normal level while managing tradeoffs like efficacy, convenience, price, and side effects. Advanced users may refer to the MATLAB run_PKPD_without_GUI.m to bypass the GUI to show additional results or edit the code directly.
