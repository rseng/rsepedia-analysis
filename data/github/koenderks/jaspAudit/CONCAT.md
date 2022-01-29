<img src="inst/icons/readmeBackground.png" width="100%" />

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02733/status.svg)](https://doi.org/10.21105/joss.02733)
[![R_build_status](https://github.com/jasp-stats/jaspAudit/workflows/unit-tests/badge.svg)](https://github.com/jasp-stats/jaspAudit/actions)
[![jfa](https://img.shields.io/cran/v/jfa?color=yellow&label=jfa&logo=r)](https://cran.r-project.org/package=jfa)
[![digitTests](https://img.shields.io/cran/v/digitTests?color=yellow&label=digitTests&logo=r)](https://cran.r-project.org/package=digitTests)

# The Audit Module

JASP for Audit (the Audit module) is an add-on module for JASP that facilitates statistical auditing. Among other things, the module provides functionality for planning, performing, evaluating, and reporting a statistical audit sample. More concretely, it contains analyses for calculating sample sizes, selecting itms according to standard audit sampling techniques, and performing inference about the population misstatement on the basis of a data sample or summary statistics of a sample. The module also features Bayesian equivalents of these analyses that enable the user to easily incorporate prior information into the statistical procedure. In all analyses, the Audit module offers explanatory text that helps the auditor in interpreting, explaining, and reporting the analysis.

- [Manual](#manual)
- [Structure](#structure)
- [Analyses](#analyses)
- [R packages](#r-packages)
- [Translations](#translations)

## Manual

For an introduction to the Audit module, please [download](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf) the manual or view it [online](https://github.com/jasp-stats/jaspAudit/blob/master/man/manual.pdf). See this [link](https://doi.org/10.21105/joss.02733) for the paper about the Audit module.

## Structure

The analyses in the Audit module are structured in JASP in the following way:

```
--- Audit
    -- Workflow
       - Sampling Workflow
       - Bayesian Sampling Workflow
    -- Planning
       - Planning
       - Bayesian Planning
    -- Selection
       - Selection
    -- Evaluation
       - Evaluation
       - Bayesian Evaluation
    -- Digit Analysis
       - Benford's Law
       - Repeated Values Analysis
    -- Other
       - True Value Estimation
```

## Analyses

### (Bayesian) Sampling Workflow

The sampling workflow is a four-stage analysis that guides the user through the process of planning, selecting, annotating, and evaluating a statistical audit sample. To make this process as easy as possible, the workflow automatically selects the appropriate options according to the type of data and audit question at hand. At the end, the sampling workflow produces a downloadable report containing the statistical results and their interpretation.

<p align="center">
  <img src='https://github.com/jasp-stats/jaspAudit/raw/master/inst/help/img/workflow.png' width='500' height='50'>
</p>

### (Bayesian) Planning

The planning analysis allows the user to, given a set of sampling objectives, calculate the minimum sample size for a statistical audit sample. By specifying a (set of) sampling objective(s), a sample size can be calculated which (when the intended sample goed according to plan) allows for a statement about the population misstatement with a prespecified amount of assurance.

<p align="center">
  <img src='https://github.com/jasp-stats/jaspAudit/raw/master/inst/help/img/workflowPlanning.png' width='500' height='50'>
</p>

### Selection

The selection analysis is an interface for the most commonly used sampling methods in auditing. The analysis allows the user to select a specified number of sampling units from the population according to one of three sampling algorithms: fixed interval sampling, cell sampling, and random sampling. Sampling units can be items (rows) or monetary units. The sample can be saved and exported to a file for futher annotation.

<p align="center">
  <img src='https://github.com/jasp-stats/jaspAudit/raw/master/inst/help/img/workflowSelection.png' width='500' height='50'>
</p>

### (Bayesian) Evaluation

The evaluation analysis allows the user to perform inference about the population misstatement on the basis of a data sample or summary statistics of a sample. The analysis contains a variety of statistical methods to evaluate the misstatement.

<p align="center">
  <img src='https://github.com/jasp-stats/jaspAudit/raw/master/inst/help/img/workflowEvaluation.png' width='500' height='50'>
</p>

### Digit Analysis
The two types of digits analyses (Benford's law and Repeated values) provide functionality for detecting irregular digit patterns in numerical data.

## `R` Packages

The functionality of the Audit module heavily mirrors that of two `R` packages:

- For the sampling analyses, the Audit module uses the [`jfa`](https://cran.r-project.org/package=jfa) package [[package website](https://koenderks.github.io/jfa/)].
- For the digit analyses, the Audit module uses the [`digitTests`](https://cran.r-project.org/package=digitTests) package [[package website](https://koenderks.github.io/digitTests/)].

<p align="center">
  <img src='https://github.com/koenderks/jfa/raw/development/man/figures/logo.png' width='149' height='173'>
  <img src='https://github.com/koenderks/digitTests/raw/development/man/figures/logo.png' width='149' height='173'>
</p>

## Translations

| Interface | Results |
| :---: | :---: |
| [![image](https://hosted.weblate.org/widgets/jasp/-/jaspaudit-qml/multi-auto.svg)](https://hosted.weblate.org/engage/jasp/) | [![image](https://hosted.weblate.org/widgets/jasp/-/jaspaudit-r/multi-auto.svg)](https://hosted.weblate.org/engage/jasp/) |
Planning
===

The planning analysis allows the user to calculate a minimum sample size given a set of sampling objectives and summary statistics of the population. Note that when you have access to the raw population data you may want to use the audit workflow, an analysis that guides you through the sampling process.

<img src="%HELP_FOLDER%/img/workflowPlanning.png" />

Please see the manual of the Audit module (download [here](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) for more detailed information about this analysis.

### Input
---

#### Sampling Objectives
- Performance materiality: Also called the upper error limit, the tolerable deviation rate, or the tolerable misstatement, the performance materiality is the upper bound of tolerable misstatement in the population to be tested. By testing against a performance materiality, you are able to plan a sample in order to collect evidence for or against the statement that the population as a whole does not contain misstatements that are considered material (i.e., are greater than the upper bound of tolerable misstatement). You should enable this objective when you want to find out whether the population contains misstatement above or below a certain limit (the performance materiality) using a sample of the population. A lower performance materiality will result in a higher required sample size. Vice versa, a higher performance materiality will result in a lower required sample size.
- Minimum precision: The precision is the the difference between the estimated most likely error and the upper bound on the misstatement. By enabling this sampling objective, you are be able to plan a sample so that the difference between the estimated most likely error and the upper bound on the misstatement is reduced to a minimum percentage. You should enable this objective if you are interested in making an estimate of the population misstatement with a certain accuracy. A lower minimum required precision will result in a higher required sample size. Vice versa, a higher minimum required precision will result in a lower required sample size.

#### Population
- No. units: The total number of units in the population. Note that the units can be items (rows) or monetary units (values) depending on the audit question.

#### Confidence
The confidence level used. The confidence level is the complement of the audit risk: the risk that the user is willing to take to give an incorrect judgment about the population. For example, if you want to have an audit risk of 5%, this equals 95% confidence.

#### Audit Risk Model
- Inherent risk: A category or probability for the inherent risk. Inherent risk is defined as the risk of material misstatement posed by an error or omission in a financial statement due to a factor other than a failure of internal control.
- Control risk: A category or probability for the internal control risk. Control risk is defined as the risk of a material misstatement in the financial statements arising due to absence or failure in the operation of relevant controls of the auditee.

When the auditor has information that indicates a low-risk profile on the population, they can use this information to reduce their required sample size via the Audit Risk Model (ARM) provided that there are no errors in the population. According to the ARM, the audit risk (AR) is a function of the inherent risk (IR), the internal control risk (CR), and the detection risk (DR).

*AR = IR x CR x DR*

The auditor assesses inherent risk and internal control risk generally on a 3-point scale to determine the appropriate detection risk. Using the ARM and zero errors the sample size depends on the risk factor *R* and the performance materiality. The risk factor *R* is a function of the detection risk (Stewart 2012).

*R = -ln(DR)*

The following table presents values of *R* as a function of the detection risk, provided that there are zero errors (Touw and Hoogduin 2012).

| Detection risk (%) | 1 | 4 | 5 | 10 | 14 |
| :---: | :---: | :---: | :---: | :---: | :---: |
| R | 4.6 | 3.2 | 3 | 2.3 | 2 |

The risk factor *R* can be adjusted using the assessments of the inherent risk and the internal control risk. By default, the standard method of setting the probabilities of IR and CR is by following the table below for a detection risk of 5%:

|  | High | Medium | Low | 
| :---: | :---: | :---: |
| R | 3 | 2 | 1 |

These values of *R* are used to set default percentages for IR and CR. The Audit module handles the following default values for IR and CR:

- High: 100%
- Medium: 60%
- Low: 36%

You can manually adjust the value of IR and CR by selecting the Custom option under the corresponding risk assessment, thus adjusting the risk factor *R*.

#### Expected errors in Sample
The expected errors are the tolerable errors that can be found in the sample while still achieving the specified sampling objectives. A sample size is calculated so that, when the number of expected errors is found in the sample, the desired confidence is retained.

*Note:* It is advised to set this value conservatively to minimize the probability of the observed errors exceeding the expected errors, which would imply that insufficient work has been done.

- Relative: Enter your expected errors as a percentage relative to the total size of the selection.
- Absolute: Enter your expected errors as the sum of (proportional) errors.

#### Probability distribution
- Hypergeometric: The hypergeometric distribution assumes a finite population size and is therefore generally used when the population size is small. It is a probability distribution that models the number of errors (*K*) in the population as a function of the population size (*N*), the number of observed found errors (*k*) and the number of correct transactions (*n*).
- Binomial: The binomial distribution assumes an infinite population size and is therefore generally used when the population size is large. It is a probability distribution that models the rate of misstatement (*\u03B8*) as a function of the observed number of errors (*k*) and the number of correct transactions (*n - k*). Because the binomial distribution strictly does not accommodate partial errors, it is generally used when you are not planning a monetary unit sample.
- Poisson: The Poisson distribution assumes an infinite population size and is therefore generally used when the population size is large. It is a probability distribution that models the rate of misstatement (*\u03B8*) as a function of the observed sample size (*n*) and the sum of the proportional errors (*t*). Because the Poisson distribution accommodates partial errors it is generally used when you are planning a monetary unit sample.

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.

#### Plots
- Compare sample sizes: Produces a plot that compares the sample size 1) across probability distributions, and 2) across the number of expected errors in the sample.
- Assumed error distribution: Produces a plot that displays the probability distribution implied by the input options and the calculated sample size.

#### Format Tables
- Numbers: Display table output as numbers.
- Percentages: Display table output as percentages.

#### Iterations
- Increment: The increment alows you to limit the possible sample sizes to a multiple of its value. For example, an increment of 5 allows only sample sizes of 5, 10, 15, 20, 25, etc.
- Maximum: The maximum allows you to limit the sample size with a maximum.

### Output
---

#### Planning Summary
- Performance materiality: When provided, the performance materiality.
- Min. precision: When provided, the minimum precision.
- Expected errors: The number (sum of proportional taints) of expected / tolerable errors in the sample.
- Minimum sample size: The minimum sample size.

#### Plots
- Compare sample sizes: Produces a plot that compares the sample size 1) across probability distributions, and 2) across the number of expected errors in the sample.
- Assumed error distribution: Produces a plot that displays the probability distribution implied by the input options and the calculated sample size.

### References
---
- AICPA (2017). <i>Audit Guide: Audit Sampling</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesian and Classical Audit Sampling. R package version 0.6.2.
- Stewart, T. (2012). <i>Technical notes on the AICPA audit guide Audit Sampling</i>. American Institute of Certified Public Accountants, New York.

### R Packages
---
- jfa
Bayesian Planning
===

The Bayesian planning analysis allows the user to calculate a minimum sample size given a set of sampling objectives and summary statistics of the population. Note that when you have access to the raw population data you may want to use the audit workflow, an analysis that guides you through the sampling process.

<img src="%HELP_FOLDER%/img/workflowPlanning.png" />

Please see the manual of the Audit module (download [here](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) for more detailed information about this analysis.

### Input
---

#### Sampling Objectives
- Performance materiality: Also called the upper error limit, the tolerable deviation rate, or the tolerable misstatement, the performance materiality is the upper bound of tolerable misstatement in the population to be tested. By testing against a performance materiality, you are able to plan a sample in order to collect evidence for or against the statement that the population as a whole does not contain misstatements that are considered material (i.e., are greater than the upper bound of tolerable misstatement). You should enable this objective when you want to find out whether the population contains misstatement above or below a certain limit (the performance materiality) using a sample of the population. A lower performance materiality will result in a higher required sample size. Vice versa, a higher performance materiality will result in a lower required sample size.
- Minimum precision: The precision is the the difference between the estimated most likely error and the upper bound on the misstatement. By enabling this sampling objective, you are be able to plan a sample so that the difference between the estimated most likely error and the upper bound on the misstatement is reduced to a minimum percentage. You should enable this objective if you are interested in making an estimate of the population misstatement with a certain accuracy. A lower minimum required precision will result in a higher required sample size. Vice versa, a higher minimum required precision will result in a lower required sample size.

#### Population
- No. units: The total number of units in the population. Note that the units can be items (rows) or monetary units (values) depending on the audit question.

#### Confidence
The confidence level used. The confidence level is the complement of the audit risk: the risk that the user is willing to take to give an incorrect judgment about the population. For example, if you want to have an audit risk of 5%, this equals 95% confidence.

#### Expected errors in Sample
The expected errors are the tolerable errors that can be found in the sample while still achieving the specified sampling objectives. A sample size is calculated so that, when the number of expected errors is found in the sample, the desired confidence is retained.

*Note:* It is advised to set this value conservatively to minimize the probability of the observed errors exceeding the expected errors, which would imply that insufficient work has been done.

- Relative: Enter your expected errors as a percentage relative to the total size of the selection.
- Absolute: Enter your expected errors as the sum of (proportional) errors.

#### Probability Distribution
- Beta-binomial: The beta-binomial distribution accompanies the hypergeometric likelihood (Dyer & Pierce, 1993). The hypergeometric likelihood assumes a finite population size and is therefore generally used when the population size is small. It is a likelihood that models the number of errors (*K*) in the population as a function of the population size (*N*), the number of observed found errors (*k*) and the number of correct transactions (*n*).
- Beta: The beta distribution accompanies the binomial likelihood. The binomial likelihood assumes an infinite population size and is therefore generally used when the population size is large. It is a likelihood that models the rate of misstatement (*\u03B8*) as a function of the observed number of errors (*k*) and the number of correct transactions (*n - k*). Because the binomial distribution strictly does not accommodate partial errors, it is generally used when you are not planning a monetary unit sample. However, the beta distribution does accommodate partial errors, and may also be used for monetary unit sampling (de Swart, Wille & Majoor, 2013).
- Gamma: The gamma distribution accompanies the Poisson likelihood. The Poisson likelihood assumes an infinite population size and is therefore generally used when the population size is large. It is a likelihood that models the rate of misstatement (*\u03B8*) as a function of the observed sample size (*n*) and the sum of the proportional errors found (*t*). Because the gamma distribution accommodates partial errors it is generally used when you are planning a monetary unit sample (Stewart, 2013).

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.

#### Tables
- Equivalent prior sample: Produces a table that displays the implicit sample on which the prior distribution is based.
- Prior and posterior: Produces a table in which the prior and expected posterior distribution are summarized through several statistics, such as their functional form, their prior and expected posterior probabilities and odds, and the shift between these.

#### Plots
- Prior and posterior: Produces a plot that shows the prior distribution and the posterior distribution after observing the intended sample.
  - Additional info: Produces dots on the materiality.
- Prior predictive: Produces a plot of the predictions of the prior distribution.
- Compare sample sizes: Produces a plot that compares the sample size 1) across probability distributions, and 2) across the number of expected errors in the sample.

#### Prior
- Default: This option does not incorporate any information into the statistical analysis and therefore assumes a negligible and conservative prior distribution.
- Manual: Provide the parameters of the prior distribution.
- Earlier sample: Create a prior distribution on the basis of an earlier sample.
  - Size: Earlier sample size.
  - Errors: Earlier found errors.
- Impartial: Create a prior distribution that is impartial with respect to the tested hypotheses.
- Risk assessments: Translate information from the audit risk model into a prior distribution.
  - Inherent risk: A category or probability for the inherent risk. Inherent risk is defined as the risk of material misstatement posed by an error or omission in a financial statement due to a factor other than a failure of internal control.
  - Control risk: A category or probability for the internal control risk. Control risk is defined as the risk of a material misstatement in the financial statements arising due to absence or failure in the operation of relevant controls of the auditee.

#### Format Tables
- Numbers: Display table output as numbers.
- Percentages: Display table output as percentages.

#### Iterations
- Increment: The increment alows you to limit the possible sample sizes to a multiple of its value. For example, an increment of 5 allows only sample sizes of 5, 10, 15, 20, 25, etc.
- Maximum: The maximum allows you to limit the sample size with a maximum.

### Output
---

#### Planning Summary
- Performance materiality: When provided, the performance materiality.
- Min. precision: When provided, the minimum precision.
- Expected errors: The number (sum of proportional taints) of expected / tolerable errors in the sample.
- Minimum sample size: The minimum sample size.

#### Equivalent prior sample
- Equivalent sample size: The sample size equivalent to the prior information.
- Equivalent errors: The number of errors equivalent to the prior information.

#### Prior and Posterior
- Functional form: The functional form of the distribution.
- Support H-: Total probability in the range of H- under the distribution. Only displayes when testing against a performance materiality.
- Support H+: Total probability in the range of H+ under the distribution. Only displayes when testing against a performance materiality.
- Ratio H- / H+: Odds in favor of H- under the distribution. Only displayes when testing against a performance materiality.
- Mean: Mean of the distribution.
- Median: Median of the distribution.
- Mode: Mode of the distribution.
- Upper bound: x-% percentile of the distribution.
- Precision: Difference between the upper bound and the mode of the distribution.

#### Plots
- Prior and posterior: Produces a plot that shows the prior distribution and the posterior distribution after observing the intended sample.
  - Additional info: Produces dots on the materiality.
- Prior predictive: Produces a plot of the predictions of the prior distribution.
- Compare sample sizes: Produces a plot that compares the sample size 1) across probability distributions, and 2) across the number of expected errors in the sample.

### References
---
- AICPA (2017). <i>Audit Guide: Audit Sampling</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesian and Classical Audit Sampling. R package version 0.6.2.
- Dyer, D., & Pierce, R. L. (1993). On the choice of the prior distribution in hypergeometric sampling. <i>Communications in Statistics-Theory and Methods</i>, 22(8), 2125-2146.
- Stewart, T. R. (2013). A Bayesian audit assurance model with application to the component materiality problem in group audits (Doctoral dissertation).
- de Swart, J., Wille, J., & Majoor, B. (2013). Het 'Push Left'-Principe als Motor van Data Analytics in de Accountantscontrole [The 'Push-Left'-Principle as a Driver of Data Analytics in Financial Audit]. <i>Maandblad voor Accountancy en Bedrijfseconomie</i>, 87, 425-432.

### R Packages
---
- jfa
Benford's Law
===

Benford's law states that the distribution of leading digits in a population naturally follows a certain distribution. In auditing, assessing whether a distribution of digits in the population conforms to Benford's law may provide additional evidence that the transactions in the population might need further investigation.

*Note:* Non-conformity to Benford's law does not necessarily indicate fraud. A Benford's law analysis should therefore only be used to acquire insight into whether a population might need further investigation.

### Input
---

#### Assignment Box
- Variable: In this box the variable is selected whose digits should be tested against the reference distribution. The value zero (0) will be omitted from the data.

#### Reference
- Benford's law: Test the digits against Benford's law.
- Uniform distribution: Test the digits against the uniform distribution.

#### Digits
- First: Checks only the first digit of the items against the specified distribution.
- First and second: Checks the first and second digit of the items against the specified distribution.
- Last: Checks only the last digit of the items against the specified distribution.

#### Bayes Factor
- BF10 : Bayes factor to quantify evidence for the alternative hypothesis relative to the null hypothesis.
- BF01 : Bayes factor to quantify evidence for the null hypothesis relative to the alternative hypothesis.
- Log(BF10) : Natural logarithm of BF10.

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.
  - Confidence: The confidence level used in the explanatory text.

### Output
---

#### Goodness-of-fit table
- n: The total number of observations in the data set.
- X<sup>2</sup>: The value of the Chi-squared test statistic.
- df: Degrees of freedom associated with the Chi-squared test.
- p: The *p* value associated with the Chi-squared test.
- BF: The Bayes factor resulting from a non-informative prior.

#### Frequency Table
- Leading / Last digit: The digit for which the information in the row applies.
- Count: The observed counts of the digits.
- Relative frequency: The observed relative frequency of the digits.
- Benford's Law / Uniform distribution: The expected relative frequency of the digits.

#### Plots
- Observed vs. expected: Produces a plot that shows the observed distribution of digits in the population compared to the expected distribution under the Benford's law or the uniform distribution.

### References
---
- Derks, K (2021). digitTests: Tests for Detecting Irregular Digit Patterns. R package version 0.1.0.

### R Packages
---
- digitTests
Waarde Schatting
===

De schattingsanalyse stelt de gebruiker in staat om de werkelijke waarde van een populatie te schatten op basis van een steekproef.

### Invoer
---

#### Opdrachtbox
- Item-ID: een unieke niet-ontbrekende identifier voor elk item in de populatie. Het rijnummer van de items is voldoende.
- Boekwaarde: De variabele die de boekwaarden van de items in de populatie bevat.
- Auditwaarde: de variabele die de audit (true) waarden bevat, of de binaire classificatie van correct (0) of incorrect (1).

#### Betrouwbaarheid
Het gebruikte betrouwbaarheidsniveau. Het betrouwbaarheidsniveau is het complement van het auditrisico: het risico dat de gebruiker bereid is te nemen om een ​​onjuist oordeel over de populatie te geven. Als u bijvoorbeeld een controlerisico van 5% wilt hebben, staat dit gelijk aan 95% betrouwbaarheid.

#### Populatie
- Aantal items: Het totale aantal items (rijen) in de populatie.
- Aantal eenheden: Het totale aantal eenheden in de populatie. Let op dat de eenheden items (rijen) of monetaire eenheden (waarden) kunnen zijn, afhankelijk van het controlevraagstuk.

#### Methode
- Directe schatter: Deze methode gebruikt alleen de controlewaarden om de afwijking te schatten (Touw en Hoogduin, 2011).
- Verschilschatter: Deze methode gebruikt het verschil tussen de boekwaarden en de controlewaarden om de afwijking in te schatten (Touw en Hoogduin, 2011).
- Ratio schatter: Deze methode gebruikt de verhouding van correctheid tussen de boekwaarden en de controlewaarden om de afwijking te schatten (Touw en Hoogduin, 2011).
- Regressieschatter: Deze methode gebruikt de lineaire relatie tussen de boekwaarden en de controlewaarden om de afwijking te schatten (Touw en Hoogduin, 2011).

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.

#### Tabellen
- Vereiste steekproefomvang: Produceert een tabel die, voor een gegeven onzekerheid, de vereiste steekproefomvang laat zien.

#### Figuren
- Spreidingsplot: Produceert een spreidingsplot die boekwaarden van de selectie vergelijkt met hun controlewaarden. Waarnemingen die fout zijn, zijn rood gekleurd.

### Uitgang
---

#### Schattingstabel
- Schatting W: De puntschatting van de totale fout in de populatie.
- Onzekerheid: De onzekerheid die samenhangt met het betrouwbaarheidsinterval.
- x-% betrouwbaarheidsinterval: het betrouwbaarheidsinterval dat bij de schatting hoort.

#### Vereiste steekproefgrootte
- Estimator: De gebruikte methode.
- Onzekerheid: het verschil tussen de meest waarschijnlijke waarde en de grenzen.
- Vereist n: Vereiste steekproefomvang.
- Extra n: Extra steekproefomvang.

### Referenties
---
- AICPA (2017). <i>Auditgids: controlesteekproeven</i>. American Institute of Certified Public Accountants.
- Touw, P., & Hoogduin, L. (2011). Statistiek voor audit en controlling.

### R-pakketten
---
- Basis RTrue Value Estimation
===

The estimation analysis allows the user to estimate the true value of a population on the basis of a sample.

### Input
---

#### Assignment Box
- Item ID: A unique non-missing identifier for every item in the population. The row number of the items is sufficient.
- Book Value: The variable that contains the book values of the items in the population.
- Audit Value: The variable that contains the audit (true) values, or the binary classification of correct (0) or incorrect (1).

#### Confidence
The confidence level used. The confidence level is the complement of the audit risk: the risk that the user is willing to take to give an incorrect judgment about the population. For example, if you want to have an audit risk of 5%, this equals 95% confidence.

#### Population
- No. items: The total number of items (rows) in the population.
- No. units: The total number of units in the population. Note that the units can be items (rows) or monetary units (values) depending on the audit question.

#### Method
- Direct estimator: This method uses only the audit values to estimate the misstatement (Touw and Hoogduin, 2011).
- Difference estimator: This method uses the difference between the book values and the audit values to estimate the misstatement (Touw and Hoogduin, 2011).
- Ratio estimator: This method uses the ratio of correctness between the book values and the audit values to estimate the misstatement (Touw and Hoogduin, 2011).
- Regression estimator: This method uses the linear relation between the book values and the audit values to estimate the misstatement (Touw and Hoogduin, 2011).

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.

#### Tables
- Required sample size: Produces a table that shows, for a given uncertainty, the required sample size.

#### Plots
- Scatter plot: Produces a scatter plot comparing book values of the selection against their audit values. Observations that are in error are colored in red.

### Output
---

#### Estimation Table
- Estimate W: The point estimate of the total error in the population.
- Uncertainty: The uncertainty associated with the confidence interval.
- x-% confidence interval: The confidence interval associated with the estimate.

#### Required Sample Size
- Estimator: The used method.
- Uncertainty: The difference between the most likely value and the bounds.
- Required n: Required sample size.
- Additional n: Additional sample size.

### References
---
- AICPA (2017). <i>Audit Guide: Audit Sampling</i>. American Institute of Certified Public Accountants.
- Touw, P., & Hoogduin, L. (2011). Statistiek voor audit en controlling.

### R Packages
---
- Base R
Bayesiaanse Steekproef Workflow
===

De taak van een auditor is om een ​​oordeel te vellen over de billijkheid van de gepresenteerde transacties in een populatie. Wanneer de auditor toegang heeft tot de onbewerkte populatiegegevens, kan hij de *auditworkflow* gebruiken om te berekenen hoeveel monsters moeten worden geëvalueerd om een ​​zeker vertrouwen in zijn oordeel te krijgen. De gebruiker kan vervolgens steekproeven nemen van deze items uit de populatie, deze items inspecteren en controleren, en statistische conclusies trekken over de afwijking in de populatie. De workflow voor steekproeven leidt de auditor door het auditproces en maakt onderweg de juiste keuzes van berekeningen.

Zie de handleiding van de Audit module (download [hier](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) voor meer gedetailleerde informatie over deze analyse.

### Werkstroom
---

- Planning: Bereken de minimale steekproefomvang om uw steekproefdoelstellingen met het gespecificeerde vertrouwen te bereiken.
- Selectie: Selecteer de gewenste steekproefeenheden uit de populatie.
- Uitvoering: Annoteer de selectie met uw beoordeling van de eerlijkheid van de geselecteerde items.
- Evaluatie: maak een bevolkingsverklaring op basis van uw geannoteerde selectie.

<img src="%HELP_FOLDER%/img/workflow.png" />

### Invoer - Planning
---

#### Steekproef Doelstellingen
- Uitvoeringsmaterialiteit: ook wel de bovengrens voor de fout, het maximaal toelaatbare afwijkingspercentage of de maximaal toelaatbare afwijking genoemd, de uitvoeringsmaterialiteit is de bovengrens van de toelaatbare afwijking in de te testen populatie. Door te toetsen aan een uitvoeringsmaterialiteit bent u in staat een steekproef te plannen om bewijs te verzamelen voor of tegen de stelling dat de populatie als geheel geen afwijkingen bevat die als materieel worden beschouwd (d.w.z groter zijn dan de uitvoeringsmaterialiteit). U moet deze doelstelling inschakelen als u met een steekproef uit de populatie wilt weten of de populatie een afwijking boven of onder een bepaalde limiet (de prestatiematerialiteit) bevat. Een lagere uitvoeringsmaterialiteit zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere uitvoeringsmaterialiteit resulteren in een lagere steekproefomvang.
- Minimale precisie: de precisie is het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de fout. Door deze steekproefdoelstelling in te schakelen, kunt u een steekproef plannen zodat het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de afwijking tot een minimumpercentage wordt teruggebracht. U moet deze doelstelling inschakelen als u geïnteresseerd bent in het maken van een schatting van de populatieafwijking met een bepaalde nauwkeurigheid. Een lagere minimaal vereiste precisie zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere minimale precisie resulteren in een vereiste steekproefomvang.

#### Betrouwbaarheid
Het gebruikte betrouwbaarheidsniveau. Het betrouwbaarheidsniveau is het complement van het auditrisico: het risico dat de gebruiker bereid is te nemen om een ​​onjuist oordeel over de populatie te geven. Als u bijvoorbeeld een controlerisico van 5% wilt hebben, staat dit gelijk aan 95% betrouwbaarheid.

#### Opdrachtbox
- Item-ID: een unieke niet-ontbrekende identifier voor elk item in de populatie. Het rijnummer van de items is voldoende.
- Boekwaarden: de variabele die de boekwaarden van de items in de populatie bevat.

#### Verwachte fouten in steekproef
De verwachte fouten zijn de toelaatbare fouten die in de steekproef kunnen worden gevonden terwijl de gespecificeerde steekproefdoelstellingen nog steeds worden bereikt. Er wordt een steekproefomvang berekend zodat, wanneer het aantal verwachte fouten in de steekproef wordt gevonden, het gewenste vertrouwen behouden blijft.

*Opmerking:* Het wordt aangeraden om deze waarde conservatief in te stellen om de kans te minimaliseren dat de waargenomen fouten de verwachte fouten overschrijden, wat zou betekenen dat er onvoldoende werk is verricht.

- Relatief: voer uw verwachte fouten in als een percentage ten opzichte van de totale grootte van de selectie.
- Absoluut: Voer uw verwachte fouten in als de som van (proportionele) fouten.

#### Kansverdeling
- Gamma: De gammaverdeling gaat samen met de Poisson-verdeling. De Poisson-verdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een verdeling waarmee het percentage afwijkingen (*\u03B8*) modelleerd als functie van de waargenomen steekproefomvang (*n*) en de som van de gevonden proportionele fouten (*t*). Omdat de gammaverdeling deelfouten mogelijk maakt, wordt deze over het algemeen gebruikt bij het plannen van een steekproef in munteenheden (Stewart, 2013).
- Beta: de bètaverdeling gaat samen met de binominaalverdeling. De binominaalverdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een verdeling waarmee het percentage afwijkingen (*\u03B8*) wordt gemodelleerd als functie van het waargenomen aantal fouten (*k*) en het aantal correcte transacties (*n - k*). Omdat de binominaalverdeling strikt genomen geen rekening houdt met gedeeltelijke fouten, wordt deze over het algemeen gebruikt wanneer u geen steekproef in een munteenheid plant. De betaverdeling is echter geschikt voor gedeeltelijke fouten en kan ook worden gebruikt voor steekproeven op monetaire eenheden (de Swart, Wille & Majoor, 2013).
- Beta-binomiaal: De beta-binomiaalverdeling gaat samen met de hypergeometrische verdeling (Dyer & Pierce, 1993). De hypergeometrische verdeling gaat uit van een eindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang klein is. Het is een verdeling waarmee het aantal fouten (*K*) in de populatie wordt gemodelleerd als functie van de populatieomvang (*N*), het aantal geobserveerde gevonden fouten (*k*) en het aantal correcte transacties (*N*).

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.

#### Tabellen
- Beschrijvende statistieken: Produceert een tabel met beschrijvende statistieken van de boekwaarden in de populatie.
- Impliciete eerdere steekproef: Produceert een tabel die de impliciete steekproef weergeeft waarop de prior verdeling is gebaseerd.
- Prior en posterior: Produceert een tabel waarin de eerdere en verwachte posterieure distributie worden samengevat door middel van verschillende statistieken, zoals hun functionele vorm, hun eerdere en verwachte posterieure kansen en kansen, en de verschuiving daartussen.

#### Figuren
- Prior en posterior: Produceert een plot die de prior verdeling en de posterieure distributie toont na observatie van het beoogde monster.
  - Extra info: Produceert stippen op de materialiteit.
- Prior predictive: Produceert een plot van de voorspellingen van de prior verdeling.
- Vergelijk steekproefomvang: Produceert een plot die de steekproefomvang vergelijkt 1) over kansverdelingen, en 2) over het aantal verwachte fouten in de steekproef.
- Verdeling van boekwaarden: Produceert een histogram van de boekwaarden in de populatie.

#### Prior
- Standaard: deze optie neemt geen informatie op in de statistische analyse en gaat daarom uit van een verwaarloosbare en conservatieve eerdere verdeling.
- Handmatig: geef de parameters van de prior verdeling op.
- Eerdere steekproef: Maak een eerdere verdeling op basis van een eerdere steekproef.
  - Grootte: eerdere steekproefomvang.
  - Fouten: eerder gevonden fouten.
- Onpartijdig: maak een voorafgaande verdeling die onpartijdig is met betrekking tot de geteste hypothesen.
- Risicobeoordelingen: Vertaal informatie uit het auditrisicomodel naar een eerdere verspreiding.
  - Inherent risico: Een categorie of waarschijnlijkheid voor het inherente risico. Inherent risico wordt gedefinieerd als het risico op een afwijking van materieel belang als gevolg van een fout of weglating in een financieel overzicht als gevolg van een andere factor dan een falen van de interne beheersing.
  - Beheersingsrisico: Een categorie of waarschijnlijkheid voor het internecontrolerisico. Het interne beheersingsrisico wordt gedefinieerd als het risico van een afwijking van materieel belang in de financiële overzichten die voortvloeit uit het ontbreken of falen van de relevante interne beheersingsmaatregelen van de gecontroleerde.

#### Kritieke items
- Negatieve boekwaarden: Isoleert negatieve boekwaarden van de populatie.
  - Bewaren: houdt negatieve boekwaarden die moeten worden geïnspecteerd in het monster.
  - Verwijderen: verwijdert negatieve boekwaarden.

#### Tabellen opmaken
- Cijfers: geef tabeluitvoer weer als getallen.
- Percentages: geef de tabeluitvoer weer als percentages.

#### Stapgrootte
Met de stapgrootte kunt u de mogelijke steekproefomvang beperken tot een veelvoud van de waarde. Een stapgrootte van 5 staat bijvoorbeeld alleen steekproefomvang van 5, 10, 15, 20, 25, enz. toe.

#### Veronderstel homogene taints
Als u op dit vakje klikt, kunt u de bekende en onbekende afwijking in de populatie scheiden om efficiënter te werken. Let op dat hiervoor de aanname vereist is dat de taints in de steekproef representatief zijn voor de taints in het onzichtbare deel van de populatie.

### Uitgang - Planning
---

#### Planningsoverzicht
- Uitvoeringsmaterialiteit: indien aanwezig, de uitvoeringsmaterialiteit.
- Min. precisie: indien aanwezig, de minimale precisie.
- Verwachte fouten: Het aantal (som van proportionele taints) verwachte / toelaatbare fouten in de steekproef.
- Minimale steekproefomvang: de minimale steekproefomvang.

#### Beschrijvende statistieken
- Populatiegrootte: aantal items in de populatie.
- Waarde: Totale waarde van de boekwaarden.
- Absolute waarde: Absolute waarde van de boekwaarden.
- Gemiddelde: gemiddelde van de boekwaarden.
- Soa. afwijking: Standaarddeviatie van de boekwaarden.
- Kwartiel: Kwartielen van de boekwaarden.

#### Impliciete voorafgaande steekproef
- Impliciete steekproefomvang: de steekproefomvang die overeenkomt met de voorafgaande informatie.
- Impliciete fouten: het aantal fouten dat gelijk is aan de eerdere informatie.

#### Prior en posterieur
- Functionele vorm: De functionele vorm van de distributie.
- Ondersteuning H-: Totale kans in het bereik van H- onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Ondersteuning H+: Totale kans in het bereik van H+ onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Verhouding H- / H+: Kansen in het voordeel van H- onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Gemiddelde: gemiddelde van de verdeling.
- Mediaan: Mediaan van de verdeling.
- Mode: Mode van de distributie.
- Bovengrens: x-% percentiel van de verdeling.
- Precisie: verschil tussen de bovengrens en de wijze van verdeling.

#### Figuren
- Prior en posterior: Produceert een plot die de prior verdeling en de posterieure distributie toont na observatie van het beoogde monster.
  - Extra info: Produceert stippen op de materialiteit.
- Prior predictive: Produceert een plot van de voorspellingen van de prior verdeling.
- Vergelijk steekproefomvang: Produceert een plot die de steekproefomvang vergelijkt 1) over kansverdelingen, en 2) over het aantal verwachte fouten in de steekproef.
- Verdeling van boekwaarden: Produceert een histogram van de boekwaarden in de populatie.

### Invoer - Selectie
---

#### Opdrachtbox
- Rangschikkingsvariabele: indien opgegeven, wordt de populatie eerst gerangschikt in oplopende volgorde met betrekking tot de waarden van deze variabele.
- Aanvullende variabelen: alle andere variabelen die in de steekproef moeten worden opgenomen.

#### Steekproefgrootte
Het vereiste aantal steekproefeenheden dat uit de populatie moet worden geselecteerd. Houd er rekening mee dat de steekproefeenheden worden bepaald door de optie *eenheden*. Als er geen boekwaarden worden opgegeven, zijn de steekproefeenheden standaard items (rijen). Wanneer boekwaarden worden verstrekt, zijn de ideale steekproefeenheden om te gebruiken monetaire eenheden.

#### Bemonsteringseenheden
- Items: voert selectie uit met behulp van de items in de populatie als steekproefeenheden.
- Monetaire eenheden: voert een selectie uit waarbij de monetaire eenheden in de populatie als steekproefeenheden worden gebruikt. Deze methode heeft de voorkeur als u meer items met een hoge waarde in de steekproef wilt opnemen.

#### Methode
- Steekproef met vast interval: voert selectie uit door de populatie in gelijke intervallen te verdelen en in elk interval een vaste eenheid te selecteren. Elk item met een waarde groter dan het interval wordt altijd in de steekproef opgenomen.
  - Startpunt: Selecteert welke bemonsteringseenheid wordt geselecteerd uit elk interval.
- Celbemonstering: voert selectie uit door de populatie in gelijke intervallen te verdelen en in elk interval een variabele eenheid te selecteren. Elk item met een waarde groter dan tweemaal het interval wordt altijd in de steekproef opgenomen.
  - Seed: Selecteert de seed voor de generator van willekeurige getallen om resultaten te reproduceren.
- Willekeurige steekproef: Voert willekeurige selectie uit waarbij elke steekproefeenheid een gelijke kans heeft om geselecteerd te worden.
  - Seed: Selecteert de seed voor de generator van willekeurige getallen om resultaten te reproduceren.

#### Items willekeurig maken
Randomiseert de items in de populatie voordat de selectie wordt uitgevoerd.

#### Tabellen
- Beschrijvende statistiek: Produceert een tabel met beschrijvende informatie over numerieke variabelen in de selectie. Statistieken die zijn opgenomen zijn het gemiddelde, de mediaan, de standaarddeviatie, de variantie, het minimum, het maximum en het bereik.
- Ruwe steekproef: produceert een tabel met de geselecteerde transacties samen met eventuele aanvullende waarnemingen in het veld met aanvullende variabelen.

### Uitvoer - Selectie
---

#### Selectie Samenvatting
- Aantal eenheden: Het aantal geselecteerde steekproefeenheden uit de populatie.
- Aantal items: Het aantal geselecteerde items uit de populatie.
- Selectiewaarde: De totale waarde van de geselecteerde items. Alleen weergegeven wanneer steekproeven op munteenheid worden gebruikt.
- % van populatieomvang / waarde: Het geselecteerde deel van de totale omvang of waarde van de populatie.

#### Beschrijvende statistieken
- Geldig: aantal geldige gevallen.
- Gemiddelde: rekenkundig gemiddelde van de gegevenspunten.
- Mediaan: Mediaan van de gegevenspunten.
- Soa. deviatie: Standaarddeviatie van de gegevenspunten.
- Variantie: Variantie van de datapunten.
- Bereik: bereik van de datapunten.
- Minimum: Minimum van de datapunten.
- Maximum: Maximum van de datapunten.

#### Ruw monster
- Rij: het rijnummer van het item.
- Geselecteerd: het aantal keren (een eenheid in) dat het item is geselecteerd.

### Invoer - Uitvoering
---

#### Annotatie
- Auditwaarde: Annoteer de items in de selectie met hun audit (true) waarden. Deze aanpak wordt aanbevolen (en automatisch geselecteerd) wanneer de items een geldwaarde hebben.
- Correct/Incorrect: Annoteer de items in de selectie met correct (0) of incorrect (1). Deze aanpak wordt aanbevolen (en automatisch geselecteerd) wanneer uw artikelen geen geldwaarde hebben.

### Invoer - Evaluatie
---

#### Opdrachtbox
- Auditresultaat / waarden: De variabele die de audit (true) waarden bevat, of de binaire classificatie van correct (0) of incorrect (1).

#### Methode
Zie *Kansverdeling*.

#### Gebied onder posterior
- Eenzijdige bovengrens: Geeft een (bovenste) schatting van de afwijking in de populatie.
- Tweezijdig interval: Geeft een (bovenste en onderste) schatting van de afwijking in de populatie.

#### Tabellen
- Correcties op populatie: Produceert een tabel die de vereiste correcties op de populatiewaarde bevat om de steekproefdoelstellingen te bereiken.
- Prior en posterior: Produceert een tabel waarin de eerdere en verwachte posterieure distributie worden samengevat door middel van verschillende statistieken, zoals hun functionele vorm, hun eerdere en verwachte posterieure kansen en kansen, en de verschuiving daartussen.
- Aannamecontroles: Produceert een tabel die de correlatie weergeeft tussen de boekwaarden in de steekproef en hun invloeden.
  - Betrouwbaarheidsinterval: Breedte van het betrouwbaarheidsinterval voor de correlatie.

#### Figuren
- Prior en posterior: Produceert een plot die de prior verdeling en de posterieure distributie toont na observatie van het beoogde monster.
  - Extra info: Produceert stippen op de materialiteit.
- Posterior predictive: Produceert een plot van de voorspellingen van de posterieure verdeling.
- Steekproefdoelstellingen: Produceert een staafdiagram waarin de materialiteit, de bovengrens van de afwijking en de meest waarschijnlijke fout (MLE) worden vergeleken.
- Spreidingsplot: Produceert een spreidingsplot die boekwaarden van de selectie vergelijkt met hun controlewaarden. Waarnemingen die fout zijn, zijn rood gekleurd.
  - Correlatie weergeven: Voegt de correlatie tussen de boekwaarden en de controlewaarden toe aan de plot.
  - Item-ID's weergeven: voegt de item-ID's toe aan de plot.

### Uitvoer - Evaluatie
---

#### Evaluatieoverzicht
- Materialiteit: indien aanwezig, de uitvoeringsmaterialiteit.
- Min. precisie: indien aanwezig, de minimale precisie.
- Steekproefomvang: De steekproefomvang (aantal eenheden).
- Fouten: het aantal foutieve elementen in de selectie.
- Taint: De som van de proportionele fouten. Gecontroleerde items kunnen worden geëvalueerd terwijl de omvang van de afwijking wordt meegenomen door hun taint te berekenen. De taint van een item *i* is het proportionele verschil tussen de boekwaarde van dat item (*y*) en de controlewaarde (true) van het item (*x*). Positieve taint worden geassocieerd met te hoge bedragen, terwijl negatieve taints optreden wanneer items worden onderschat.
<img src="%HELP_FOLDER%/img/taints.png" />
- Meest waarschijnlijke fout: De meest waarschijnlijke fout in de populatie.
- x-% Betrouwbaarheidsgrens: Bovengrens van de afwijking in de populatie.
- Precisie: verschil tussen bovengrens en meest waarschijnlijke fout.
- BF-+: De Bayes-factor voor de test.

#### Prior en posterieur
- Functionele vorm: De functionele vorm van de distributie.
- Ondersteuning H-: Totale kans in het bereik van H- onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Ondersteuning H+: Totale kans in het bereik van H+ onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Verhouding H- / H+: Kansen in het voordeel van H- onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Gemiddelde: gemiddelde van de verdeling.
- Mediaan: Mediaan van de verdeling.
- Mode: Mode van de distributie.
- Bovengrens: x-% percentiel van de verdeling.
- Precisie: verschil tussen de bovengrens en de wijze van verdeling.

#### Correcties voor de bevolking
- Correctie: Het van de populatie af te trekken bedrag of percentage.

#### Aannamecontroles
- n: steekproefomvang.
- Pearsons r: Pearson-correlatiecoëfficiënt.
- x-% bovengrens: bovengrens voor correlatiecoëfficiënt.
- p: p-waarde voor de test.
- BF-0: Bayes-factor voor de test.

#### Figuren
- Prior en posterior: Produceert een plot die de prior verdeling en de posterieure distributie toont na observatie van het beoogde monster.
  - Extra info: Produceert stippen op de materialiteit.
- Posterior predictive: Produceert een plot van de voorspellingen van de posterieure verdeling.
- Steekproefdoelstellingen: Produceert een staafdiagram waarin de materialiteit, de bovengrens van de afwijking en de meest waarschijnlijke fout (MLE) worden vergeleken.
- Spreidingsplot: Produceert een spreidingsplot die boekwaarden van de selectie vergelijkt met hun controlewaarden. Waarnemingen die fout zijn, zijn rood gekleurd.
  - Correlatie weergeven: Voegt de correlatie tussen de boekwaarden en de controlewaarden toe aan de plot.
  - Item-ID's weergeven: voegt de item-ID's toe aan de plot.

### Referenties
---
- AICPA (2017). <i>Auditgids: controlesteekproeven</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesian and Classical Audit Sampling. R package version 0.6.2.
- Dyer, D., & Pierce, R.L. (1993). Over de keuze van de voorafgaande verdeling bij hypergeometrische steekproeven. <i>Communicatie in statistiek-theorie en methoden</i>, 22(8), 2125-2146.
- Stewart, TR (2013). Een Bayesiaans audit assurance-model met toepassing op het component materialiteitsprobleem bij groepsaudits (proefschrift).
- de Swart, J., Wille, J., & Majoor, B. (2013). Het 'Push-Left'-Principe als Motor van Data Analytics in de Accountantcontrole [Het 'Push-Left'-Principe als aanjager van Data Analytics in Financial Audit]. <i>Maandblad voor Accountancy en Bedrijfseconomie</i>, 87, 425-432.

### R-pakketten
---
- jfaSelectie
===

De selectieanalyse stelt de gebruiker in staat om een ​​aantal steekproefeenheden uit een populatie te selecteren met behulp van een combinatie van steekproeftechnieken (recordsteekproef versus geldeenheidsteekproef) en steekproefmethoden (willekeurige steekproef, celsteekproef, vaste intervalsteekproef) die standaard zijn in een audit context.

<img src="%HELP_FOLDER%/img/workflowSelection.png" />

Zie de handleiding van de Audit module (download [hier](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) voor meer gedetailleerde informatie over deze analyse.

### Invoer
---

#### Opdrachtbox
- Item-ID: een unieke niet-ontbrekende identifier voor elk item in de populatie. Het rijnummer van de items is voldoende.
- Boekwaarden: de variabele die de boekwaarden van de items in de populatie bevat.
- Rangschikkingsvariabele: indien opgegeven, wordt de populatie eerst gerangschikt in oplopende volgorde met betrekking tot de waarden van deze variabele.
- Aanvullende variabelen: alle andere variabelen die in de steekproef moeten worden opgenomen.

#### Steekproefgrootte
Het vereiste aantal steekproefeenheden dat uit de populatie moet worden geselecteerd. Houd er rekening mee dat de steekproefeenheden worden bepaald door de optie *eenheden*. Als er geen boekwaarden worden opgegeven, zijn de steekproefeenheden standaard items (rijen). Wanneer boekwaarden worden verstrekt, zijn de ideale steekproefeenheden om te gebruiken monetaire eenheden.

#### Bemonsteringseenheden
- Items: voert selectie uit met behulp van de items in de populatie als steekproefeenheden.
- Monetaire eenheden: voert een selectie uit waarbij de monetaire eenheden in de populatie als steekproefeenheden worden gebruikt. Deze methode heeft de voorkeur als u meer items met een hoge waarde in de steekproef wilt opnemen.

#### Methode
- Steekproef met vast interval: voert selectie uit door de populatie in gelijke intervallen te verdelen en in elk interval een vaste eenheid te selecteren. Elk item met een waarde groter dan het interval wordt altijd in de steekproef opgenomen.
  - Startpunt: Selecteert welke bemonsteringseenheid wordt geselecteerd uit elk interval.
- Celbemonstering: voert selectie uit door de populatie in gelijke intervallen te verdelen en in elk interval een variabele eenheid te selecteren. Elk item met een waarde groter dan tweemaal het interval wordt altijd in de steekproef opgenomen.
  - Seed: Selecteert de seed voor de generator van willekeurige getallen om resultaten te reproduceren.
- Willekeurige steekproef: Voert willekeurige selectie uit waarbij elke steekproefeenheid een gelijke kans heeft om geselecteerd te worden.
  - Seed: Selecteert de seed voor de generator van willekeurige getallen om resultaten te reproduceren.

#### Items willekeurig maken
Randomiseert de items in de populatie voordat de selectie wordt uitgevoerd.

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.

#### Kolomnaam Selectie Resultaat
Wanneer een naam is opgegeven, wordt het resultaat van de selectieanalyse in een nieuwe kolom aan de gegevens toegevoegd. De nieuwe kolom geeft weer hoe vaak elke transactie in de steekproef is opgenomen.

#### Tabellen
- Beschrijvende statistiek: Produceert een tabel met beschrijvende informatie over numerieke variabelen in de selectie. Statistieken die zijn opgenomen zijn het gemiddelde, de mediaan, de standaarddeviatie, de variantie, het minimum, het maximum en het bereik.
- Ruwe steekproef: produceert een tabel met de geselecteerde transacties samen met eventuele aanvullende waarnemingen in het veld met aanvullende variabelen.

### Uitgang
---

#### Selectie Samenvatting
- Aantal eenheden: Het aantal geselecteerde steekproefeenheden uit de populatie.
- Aantal items: Het aantal geselecteerde items uit de populatie.
- Selectiewaarde: De totale waarde van de geselecteerde items. Alleen weergegeven wanneer steekproeven op munteenheid worden gebruikt.
- % van populatieomvang / waarde: Het geselecteerde deel van de totale omvang of waarde van de populatie.

#### Informatie over monetaire intervalselectie
- Items: Het aantal items in de populatie.
- Waarde: De waarde van de items in de populatie.
- Geselecteerde items: het aantal items in de steekproef.
- Geselecteerde eenheden: Het aantal geselecteerde eenheden uit de populatie.
- Selectiewaarde: De waarde van de items in de steekproef.
- % van totale waarde: het geselecteerde aandeel van de totale waarde van de items in vergelijking met de items in de populatie.

#### Beschrijvende statistieken
- Geldig: aantal geldige gevallen.
- Gemiddelde: rekenkundig gemiddelde van de gegevenspunten.
- Mediaan: Mediaan van de gegevenspunten.
- Soa. deviatie: Standaarddeviatie van de gegevenspunten.
- Variantie: Variantie van de datapunten.
- Bereik: bereik van de datapunten.
- Minimum: Minimum van de datapunten.
- Maximum: Maximum van de datapunten.

#### Ruw monster
- Rij: het rijnummer van het item.
- Geselecteerd: het aantal keren (een eenheid in) dat het item is geselecteerd.

### Referenties
---
- AICPA (2017). <i>Auditgids: controlesteekproeven</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesiaanse en klassieke auditsteekproeven. R-pakket versie 0.6.2.

### R-pakketten
---
- jfaEvaluation
===

The evaluation analysis allows the user to perform inference about the total misstatement in the population on the basis of an audit sample.

<img src="%HELP_FOLDER%/img/workflowEvaluation.png" />

Please see the manual of the Audit module (download [here](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) for more detailed information about this analysis.

### Input
---

#### Sampling Objectives
- Performance materiality: Also called the upper error limit, the tolerable deviation rate, or the tolerable misstatement, the performance materiality is the upper bound of tolerable misstatement in the population to be tested. By testing against a performance materiality, you are able to plan a sample in order to collect evidence for or against the statement that the population as a whole does not contain misstatements that are considered material (i.e., are greater than the upper bound of tolerable misstatement). You should enable this objective when you want to find out whether the population contains misstatement above or below a certain limit (the performance materiality) using a sample of the population. A lower performance materiality will result in a higher required sample size. Vice versa, a higher performance materiality will result in a lower required sample size.
- Minimum precision: The precision is the the difference between the estimated most likely error and the upper bound on the misstatement. By enabling this sampling objective, you are be able to plan a sample so that the difference between the estimated most likely error and the upper bound on the misstatement is reduced to a minimum percentage. You should enable this objective if you are interested in making an estimate of the population misstatement with a certain accuracy. A lower minimum required precision will result in a higher required sample size. Vice versa, a higher minimum required precision will result in a lower required sample size.

#### Population
- No. items: The total number of items (rows) in the population.
- No. units: The total number of units in the population. Note that the units can be items (rows) or monetary units (values) depending on the audit question.

#### Confidence
The confidence level used. The confidence level is the complement of the audit risk: the risk that the user is willing to take to give an incorrect judgment about the population. For example, if you want to have an audit risk of 5%, this equals 95% confidence.

#### Assignment Box
- Item ID: A unique non-missing identifier for every item in the population. The row number of the items is sufficient.
- Book Values: The variable that contains the book values of the items in the population.
- Audit result / values: The variable that contains the audit (true) values, or the binary classification of correct (0) or incorrect (1).
- Selection counter: The variable that contains how many times each observation should be evaluated.

#### Data
- Raw: Use raw data.
- Summary statistics: Use summary statistics.

#### Audit Risk Model
- Inherent risk: A category or probability for the inherent risk. Inherent risk is defined as the risk of material misstatement posed by an error or omission in a financial statement due to a factor other than a failure of internal control.
- Control risk: A category or probability for the internal control risk. Control risk is defined as the risk of a material misstatement in the financial statements arising due to absence or failure in the operation of relevant controls of the auditee.

When the auditor has information that indicates a low-risk profile on the population, they can use this information to reduce their required sample size via the Audit Risk Model (ARM) provided that there are no errors in the population. According to the ARM, the audit risk (AR) is a function of the inherent risk (IR), the internal control risk (CR), and the detection risk (DR).

*AR = IR x CR x DR*

The auditor assesses inherent risk and internal control risk generally on a 3-point scale to determine the appropriate detection risk. Using the ARM and zero errors the sample size depends on the risk factor *R* and the performance materiality. The risk factor *R* is a function of the detection risk (Stewart 2012).

*R = -ln(DR)*

The following table presents values of *R* as a function of the detection risk, provided that there are zero errors (Touw and Hoogduin 2012).

| Detection risk (%) | 1 | 4 | 5 | 10 | 14 |
| :---: | :---: | :---: | :---: | :---: | :---: |
| R | 4.6 | 3.2 | 3 | 2.3 | 2 |

The risk factor *R* can be adjusted using the assessments of the inherent risk and the internal control risk. By default, the standard method of setting the probabilities of IR and CR is by following the table below for a detection risk of 5%:

|  | High | Medium | Low | 
| :---: | :---: | :---: |
| R | 3 | 2 | 1 |

These values of *R* are used to set default percentages for IR and CR. The Audit module handles the following default values for IR and CR:

- High: 100%
- Medium: 60%
- Low: 36%

You can manually adjust the value of IR and CR by selecting the Custom option under the corresponding risk assessment, thus adjusting the risk factor *R*.

#### Method
- Poisson: Uses the Poisson likelhood to evaluate the sample.
- Binomial: Uses the binomial likelhood to evaluate the sample.
- Hypergeometric: Uses the hypergeometric likelhood to evaluate the sample.
- Stringer: The Stringer bound to evaluate the sample (Stringer, 1963).
  - LTA adjustment: LTA adjustment for the stringer bound to incorporate understatements (Leslie, Teitlebaum, & Anderson, 1979).
- Mean-per-unit estimator: Uses the mean-per-unit estimator.
- Direct estimator: This method uses only the audit values to estimate the misstatement (Touw and Hoogduin, 2011).
- Difference estimator: This method uses the difference between the book values and the audit values to estimate the misstatement (Touw and Hoogduin, 2011).
- Ratio estimator: This method uses the ratio of correctness between the book values and the audit values to estimate the misstatement (Touw and Hoogduin, 2011).
- Regression estimator: This method uses the linear relation between the book values and the audit values to estimate the misstatement (Touw and Hoogduin, 2011).

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.

#### Tables
- Corrections to population: Produces a table that contains the required corrections to the population value to achieve the sampling objectives.

#### Plots
- Sampling objectives: Produces a bar chart comparing the materiality, maximum misstatement and most likely error (MLE).
- Scatter plot: Produces a scatter plot comparing book values of the selection against their audit values. Observations that are in error are colored in red.

#### Critical Items
- Negative book values: Isolates negative book values from the population.
  - Keep: Keeps negative book values to be inspected in the sample.
  - Remove: Removes negative book values.

#### Format Tables
- Numbers: Display table output as numbers.
- Percentages: Display table output as percentages.
- Monetary values: Display table output as monetary values.

### Output
---

#### Evaluation summary
- Materiality: When provided, the performance materiality.
- Min. precision: When provided, the minimum precision.
- Sample size: The sample size (number of units).
- Errors: The number of erroneous elements in the selection.
- Taint: The sum of the proportional errors. Audited items can be evaluated while incorporating the magnitude of the misstatement by calculating their taints. The taint of an item *i* is the proportional difference between that item's book value (*y*) and the item's audit (true) value (*x*). Positive taints are associated with overstatements, while negative taints occur when items are understated.
<img src="%HELP_FOLDER%/img/taints.png" />
- Most likely error: The most likely error in the population.
- x-% Confidence bound: Upper bound on the misstatement in the population.
- Precision: Difference between upper bound and most likely error.
- p: The p-value for the test.

#### Corrections to Population
- Correction: The amount or percentage to be deducted from the population.

#### Plots
- Sampling objectives: Produces a bar chart comparing the materiality, upper bound on the misstatement and most likely error (MLE).
- Scatter plot: Produces a scatter plot comparing book values of the selection against their audit values. Observations that are in error are colored in red.
  - Display correlation: Adds the correlation between the book values and the audit values to the plot.
  - Display item ID's: Adds the item ID's to the plot.

### References
---
- AICPA (2017). <i>Audit Guide: Audit Sampling</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesian and Classical Audit Sampling. R package version 0.6.2.
- Leslie, D. A., Teitlebaum, A. D., Anderson, R. J. (1979). <i>Dollar-unit Sampling: A Practical Guide for Auditors</i>. Toronto: Copp Clark Pitman.
- Stringer, K. W. (1963) Practical aspects of statistical sampling in auditing. <i>Proceedings of Business and Economic Statistics Section</i>, American Statistical Association.
- Touw, P., & Hoogduin, L. (2011). Statistiek voor audit en controlling.

### R Packages
---
- jfa
Steekproef Workflow
===

De taak van een auditor is om een ​​oordeel te vellen over de billijkheid van de gepresenteerde transacties in een populatie. Wanneer de auditor toegang heeft tot de onbewerkte populatiegegevens, kan hij de *auditworkflow* gebruiken om te berekenen hoeveel monsters moeten worden geëvalueerd om een ​​zeker vertrouwen in zijn oordeel te krijgen. De gebruiker kan vervolgens steekproeven nemen van deze items uit de populatie, deze items inspecteren en controleren, en statistische conclusies trekken over de afwijking in de populatie. De workflow voor steekproeven leidt de auditor door het auditproces en maakt onderweg de juiste keuzes van berekeningen.

Zie de handleiding van de Audit module (download [hier](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) voor meer gedetailleerde informatie over deze analyse.

### Werkstroom
---

- Planning: Bereken de minimale steekproefomvang om uw steekproefdoelstellingen met het gespecificeerde vertrouwen te bereiken.
- Selectie: Selecteer de gewenste steekproefeenheden uit de populatie.
- Uitvoering: Annoteer de selectie met uw beoordeling van de eerlijkheid van de geselecteerde items.
- Evaluatie: maak een bevolkingsverklaring op basis van uw geannoteerde selectie.

<img src="%HELP_FOLDER%/img/workflow.png" />

### Invoer - Planning
---

#### Steekproef Doelstellingen
- Uitvoeringsmaterialiteit: ook wel de bovengrens voor de fout, het maximaal toelaatbare afwijkingspercentage of de maximaal toelaatbare afwijking genoemd, de uitvoeringsmaterialiteit is de bovengrens van de toelaatbare afwijking in de te testen populatie. Door te toetsen aan een uitvoeringsmaterialiteit bent u in staat een steekproef te plannen om bewijs te verzamelen voor of tegen de stelling dat de populatie als geheel geen afwijkingen bevat die als materieel worden beschouwd (d.w.z groter zijn dan de uitvoeringsmaterialiteit). U moet deze doelstelling inschakelen als u met een steekproef uit de populatie wilt weten of de populatie een afwijking boven of onder een bepaalde limiet (de prestatiematerialiteit) bevat. Een lagere uitvoeringsmaterialiteit zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere uitvoeringsmaterialiteit resulteren in een lagere steekproefomvang.
- Minimale precisie: de precisie is het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de fout. Door deze steekproefdoelstelling in te schakelen, kunt u een steekproef plannen zodat het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de afwijking tot een minimumpercentage wordt teruggebracht. U moet deze doelstelling inschakelen als u geïnteresseerd bent in het maken van een schatting van de populatieafwijking met een bepaalde nauwkeurigheid. Een lagere minimaal vereiste precisie zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere minimale precisie resulteren in een vereiste steekproefomvang.

#### Betrouwbaarheid
Het gebruikte betrouwbaarheidsniveau. Het betrouwbaarheidsniveau is het complement van het auditrisico: het risico dat de gebruiker bereid is te nemen om een ​​onjuist oordeel over de populatie te geven. Als u bijvoorbeeld een controlerisico van 5% wilt hebben, staat dit gelijk aan 95% betrouwbaarheid.

#### Opdrachtbox
- Item-ID: een unieke niet-ontbrekende identifier voor elk item in de populatie. Het rijnummer van de items is voldoende.
- Boekwaarden: de variabele die de boekwaarden van de items in de populatie bevat.

#### Auditrisicomodel
- Inherent risico: Een categorie of waarschijnlijkheid voor het inherente risico. Inherent risico wordt gedefinieerd als het risico op een afwijking van materieel belang als gevolg van een fout of weglating in een financieel overzicht als gevolg van een andere factor dan een falen van de interne beheersing.
- Beheersingsrisico: Een categorie of waarschijnlijkheid voor het internecontrolerisico. Het interne beheersingsrisico wordt gedefinieerd als het risico van een afwijking van materieel belang in de financiële overzichten die voortvloeit uit het ontbreken of falen van de relevante interne beheersingsmaatregelen van de gecontroleerde.

Wanneer de auditor informatie heeft die wijst op een laag risicoprofiel van de populatie, kan hij deze informatie gebruiken om zijn vereiste steekproefomvang te verkleinen via het Audit Risk Model (ARM), op voorwaarde dat er geen fouten in de populatie zitten. Volgens de ARM is het auditrisico (AR) een functie van het inherente risico (IR), het internecontrolerisico (CR) en het detectierisico (DR).

*AR = IR x CR x DR*

De auditor beoordeelt het inherente risico en het internecontrolerisico in het algemeen op een 3-puntsschaal om het juiste ontdekkingsrisico te bepalen. Met behulp van de ARM en nul fouten hangt de steekproefomvang af van de risicofactor *R* en de uitvoeringsmaterialiteit. De risicofactor *R* is een functie van het detectierisico (Stewart 2012).

*R = -ln(DR)*

De volgende tabel geeft de waarden van *R* weer als functie van het detectierisico, op voorwaarde dat er geen fouten zijn (Touw en Hoogduin 2012).

| Detectierisico (%) | 1 | 4 | 5 | 10 | 14 |
| :---: | :---: | :---: | :---: | :---: | :---: |
| R | 4.6 | 3.2 | 3 | 2.3 | 2 |

Deze waarden van *R* worden gebruikt om standaardpercentages in te stellen voor IR en CR. De Audit-module verwerkt de volgende standaardwaarden voor IR en CR:

- Hoog: 100%
- Gemiddeld: 60%
- Laag: 36%

U kunt de waarde van IR en CR handmatig aanpassen door de optie Aangepast te selecteren onder de bijbehorende risicobeoordeling, waardoor de risicofactor *R* wordt aangepast.

#### Verwachte fouten in steekproef
De verwachte fouten zijn de toelaatbare fouten die in de steekproef kunnen worden gevonden terwijl de gespecificeerde steekproefdoelstellingen nog steeds worden bereikt. Er wordt een steekproefomvang berekend zodat, wanneer het aantal verwachte fouten in de steekproef wordt gevonden, het gewenste vertrouwen behouden blijft.

*Opmerking:* Het wordt aangeraden om deze waarde conservatief in te stellen om de kans te minimaliseren dat de waargenomen fouten de verwachte fouten overschrijden, wat zou betekenen dat er onvoldoende werk is verricht.

- Relatief: voer uw verwachte fouten in als een percentage ten opzichte van de totale grootte van de selectie.
- Absoluut: Voer uw verwachte fouten in als de som van (proportionele) fouten.

#### Kansverdeling
- Poisson: De Poisson-verdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een kansverdeling die het percentage afwijkingen (*\u03B8*) modelleert als functie van de waargenomen steekproefomvang (*n*) en de som van de proportionele fouten (*t*). Omdat de Poisson-verdeling rekening houdt met gedeeltelijke fouten, wordt deze over het algemeen gebruikt wanneer u een steekproef in een munteenheid plant.
- Binomiaal: De binomiale verdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een kansverdeling die het afwijkingspercentage (*\u03B8*) modelleert als functie van het waargenomen aantal fouten (*k*) en het aantal correcte transacties (*n - k*). Omdat de binominale verdeling strikt geen rekening houdt met gedeeltelijke fouten, wordt deze over het algemeen gebruikt wanneer u geen steekproef in een munteenheid plant.
- Hypergeometrische verdeling: De hypergeometrische verdeling gaat uit van een eindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang klein is. Het is een kansverdeling die het aantal fouten (*K*) in de populatie modelleert als functie van de populatieomvang (*N*), het aantal waargenomen gevonden fouten (*k*) en het aantal correcte transacties ( *N*).

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.

#### Tabellen
- Beschrijvende statistieken: Produceert een tabel met beschrijvende statistieken van de boekwaarden in de populatie.

#### Figuren
- Verdeling van boekwaarden: Produceert een histogram van de boekwaarden in de populatie.
- Vergelijk steekproefomvang: Produceert een plot die de steekproefomvang vergelijkt 1) over kansverdelingen, en 2) over het aantal verwachte fouten in de steekproef.
- Veronderstelde foutenverdeling: Produceert een plot die de kansverdeling weergeeft die wordt geïmpliceerd door de invoeropties en de berekende steekproefomvang.

#### Kritieke items
- Negatieve boekwaarden: Isoleert negatieve boekwaarden van de populatie.
  - Bewaren: houdt negatieve boekwaarden die moeten worden geïnspecteerd in het monster.
  - Verwijderen: verwijdert negatieve boekwaarden.

#### Tabellen opmaken
- Cijfers: geef tabeluitvoer weer als getallen.
- Percentages: geef de tabeluitvoer weer als percentages.

#### Stapgrootte
Met de stapgrootte kunt u de mogelijke steekproefomvang beperken tot een veelvoud van de waarde. Een stapgrootte van 5 staat bijvoorbeeld alleen steekproefomvang van 5, 10, 15, 20, 25, enz. toe.

### Uitgang - Planning
---

#### Planningsoverzicht
- Uitvoeringsmaterialiteit: indien aanwezig, de uitvoeringsmaterialiteit.
- Min. precisie: indien aanwezig, de minimale precisie.
- Verwachte fouten: Het aantal (som van proportionele taints) verwachte / toelaatbare fouten in de steekproef.
- Minimale steekproefomvang: de minimale steekproefomvang.

#### Beschrijvende statistieken
- Populatiegrootte: aantal items in de populatie.
- Waarde: Totale waarde van de boekwaarden.
- Absolute waarde: Absolute waarde van de boekwaarden.
- Gemiddelde: gemiddelde van de boekwaarden.
- Soa. afwijking: Standaarddeviatie van de boekwaarden.
- Kwartiel: Kwartielen van de boekwaarden.

#### Figuren
- Verdeling van boekwaarden: Produceert een histogram van de boekwaarden in de populatie.
- Vergelijk steekproefomvang: Produceert een plot die de steekproefomvang vergelijkt 1) over kansverdelingen, en 2) over het aantal verwachte fouten in de steekproef.
- Veronderstelde foutenverdeling: Produceert een plot die de kansverdeling weergeeft die wordt geïmpliceerd door de invoeropties en de berekende steekproefomvang.

### Invoer - Selectie
---

#### Opdrachtbox
- Rangschikkingsvariabele: indien opgegeven, wordt de populatie eerst gerangschikt in oplopende volgorde met betrekking tot de waarden van deze variabele.
- Aanvullende variabelen: alle andere variabelen die in de steekproef moeten worden opgenomen.

#### Steekproefgrootte
Het vereiste aantal steekproefeenheden dat uit de populatie moet worden geselecteerd. Houd er rekening mee dat de steekproefeenheden worden bepaald door de optie *eenheden*. Als er geen boekwaarden worden opgegeven, zijn de steekproefeenheden standaard items (rijen). Wanneer boekwaarden worden verstrekt, zijn de ideale steekproefeenheden om te gebruiken monetaire eenheden.

#### Bemonsteringseenheden
- Items: voert selectie uit met behulp van de items in de populatie als steekproefeenheden.
- Monetaire eenheden: voert een selectie uit waarbij de monetaire eenheden in de populatie als steekproefeenheden worden gebruikt. Deze methode heeft de voorkeur als u meer items met een hoge waarde in de steekproef wilt opnemen.

#### Methode
- Steekproef met vast interval: voert selectie uit door de populatie in gelijke intervallen te verdelen en in elk interval een vaste eenheid te selecteren. Elk item met een waarde groter dan het interval wordt altijd in de steekproef opgenomen.
  - Startpunt: Selecteert welke bemonsteringseenheid wordt geselecteerd uit elk interval.
- Celbemonstering: voert selectie uit door de populatie in gelijke intervallen te verdelen en in elk interval een variabele eenheid te selecteren. Elk item met een waarde groter dan tweemaal het interval wordt altijd in de steekproef opgenomen.
  - Seed: Selecteert de seed voor de generator van willekeurige getallen om resultaten te reproduceren.
- Willekeurige steekproef: Voert willekeurige selectie uit waarbij elke steekproefeenheid een gelijke kans heeft om geselecteerd te worden.
  - Seed: Selecteert de seed voor de generator van willekeurige getallen om resultaten te reproduceren.

#### Items willekeurig maken
Randomiseert de items in de populatie voordat de selectie wordt uitgevoerd.

#### Tabellen
- Beschrijvende statistiek: Produceert een tabel met beschrijvende informatie over numerieke variabelen in de selectie. Statistieken die zijn opgenomen zijn het gemiddelde, de mediaan, de standaarddeviatie, de variantie, het minimum, het maximum en het bereik.
- Ruwe steekproef: produceert een tabel met de geselecteerde transacties samen met eventuele aanvullende waarnemingen in het veld met aanvullende variabelen.

### Uitvoer - Selectie
---

#### Selectie Samenvatting
- Aantal eenheden: Het aantal geselecteerde steekproefeenheden uit de populatie.
- Aantal items: Het aantal geselecteerde items uit de populatie.
- Selectiewaarde: De totale waarde van de geselecteerde items. Alleen weergegeven wanneer steekproeven op munteenheid worden gebruikt.
- % van populatieomvang / waarde: Het geselecteerde deel van de totale omvang of waarde van de populatie.

#### Informatie over monetaire intervalselectie
- Items: Het aantal items in de populatie.
- Waarde: De waarde van de items in de populatie.
- Geselecteerde items: het aantal items in de steekproef.
- Geselecteerde eenheden: Het aantal geselecteerde eenheden uit de populatie.
- Selectiewaarde: De waarde van de items in de steekproef.
- % van totale waarde: het geselecteerde aandeel van de totale waarde van de items in vergelijking met de items in de populatie.

#### Beschrijvende statistieken
- Geldig: aantal geldige gevallen.
- Gemiddelde: rekenkundig gemiddelde van de gegevenspunten.
- Mediaan: Mediaan van de gegevenspunten.
- Soa. deviatie: Standaarddeviatie van de gegevenspunten.
- Variantie: Variantie van de datapunten.
- Bereik: bereik van de datapunten.
- Minimum: Minimum van de datapunten.
- Maximum: Maximum van de datapunten.

#### Ruw monster
- Rij: het rijnummer van het item.
- Geselecteerd: het aantal keren (een eenheid in) dat het item is geselecteerd.

### Invoer - Uitvoering
---

#### Annotatie
- Auditwaarde: Annoteer de items in de selectie met hun audit (true) waarden. Deze aanpak wordt aanbevolen (en automatisch geselecteerd) wanneer de items een geldwaarde hebben.
- Correct/Incorrect: Annoteer de items in de selectie met correct (0) of incorrect (1). Deze aanpak wordt aanbevolen (en automatisch geselecteerd) wanneer uw artikelen geen geldwaarde hebben.

### Invoer - Evaluatie
---

#### Opdrachtbox
- Auditresultaat / waarden: De variabele die de audit (true) waarden bevat, of de binaire classificatie van correct (0) of incorrect (1).

#### Methode
- Poisson: gebruikt de Poisson-waarschijnlijkheid om het monster te evalueren.
- Binomiaal: gebruikt de binominale waarschijnlijkheid om het monster te evalueren.
- Hypergeometrische: gebruikt de hypergeometrische waarschijnlijkheid om het monster te evalueren.
- Stringer: de Stringer moet het monster evalueren (Stringer, 1963).
  - LTA-aanpassing: LTA-aanpassing voor de stringer die onvermijdelijk understatements bevat (Leslie, Teitlebaum, & Anderson, 1979).
- Gemiddelde-per-eenheid schatter: Gebruikt de gemiddelde-per-eenheid schatter.
- Directe schatter: Deze methode gebruikt alleen de controlewaarden om de afwijking te schatten (Touw en Hoogduin, 2011).
- Verschilschatter: Deze methode gebruikt het verschil tussen de boekwaarden en de controlewaarden om de afwijking in te schatten (Touw en Hoogduin, 2011).
- Ratio schatter: Deze methode gebruikt de verhouding van correctheid tussen de boekwaarden en de controlewaarden om de afwijking te schatten (Touw en Hoogduin, 2011).
- Regressieschatter: Deze methode gebruikt de lineaire relatie tussen de boekwaarden en de controlewaarden om de afwijking te schatten (Touw en Hoogduin, 2011).

#### Tabellen
- Correcties op populatie: Produceert een tabel die de vereiste correcties op de populatiewaarde bevat om de steekproefdoelstellingen te bereiken.

#### Figuren
- Steekproefdoelstellingen: Produceert een staafdiagram waarin de materialiteit, de maximale afwijking en de meest waarschijnlijke fout (MLE) worden vergeleken.
- Spreidingsplot: Produceert een spreidingsplot die boekwaarden van de selectie vergelijkt met hun controlewaarden. Waarnemingen die fout zijn, zijn rood gekleurd.

### Uitvoer - Evaluatie
---

#### Evaluatieoverzicht
- Materialiteit: indien aanwezig, de uitvoeringsmaterialiteit.
- Min. precisie: indien aanwezig, de minimale precisie.
- Steekproefomvang: De steekproefomvang (aantal eenheden).
- Fouten: het aantal foutieve elementen in de selectie.
- Taint: De som van de proportionele fouten. Gecontroleerde items kunnen worden geëvalueerd terwijl de omvang van de afwijking wordt meegenomen door hun taint te berekenen. De taint van een item *i* is het proportionele verschil tussen de boekwaarde van dat item (*y*) en de controlewaarde (true) van het item (*x*). Positieve taint worden geassocieerd met te hoge bedragen, terwijl negatieve taints optreden wanneer items worden onderschat.
<img src="%HELP_FOLDER%/img/taints.png" />
- Meest waarschijnlijke fout: De meest waarschijnlijke fout in de populatie.
- x-% Betrouwbaarheidsgrens: Bovengrens van de afwijking in de populatie.
- Precisie: verschil tussen bovengrens en meest waarschijnlijke fout.
- p: De p-waarde voor de test.

#### Correcties voor de bevolking
- Correctie: Het van de populatie af te trekken bedrag of percentage.

#### Figuren
- Steekproefdoelstellingen: Produceert een staafdiagram waarin de materialiteit, de bovengrens van de afwijking en de meest waarschijnlijke fout (MLE) worden vergeleken.
- Spreidingsplot: Produceert een spreidingsplot die boekwaarden van de selectie vergelijkt met hun controlewaarden. Waarnemingen die fout zijn, zijn rood gekleurd.
  - Correlatie weergeven: Voegt de correlatie tussen de boekwaarden en de controlewaarden toe aan de plot.
  - Item-ID's weergeven: voegt de item-ID's toe aan de plot.

### Referenties
---
- AICPA (2017). <i>Auditgids: controlesteekproeven</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesiaanse en klassieke auditsteekproeven. R-pakket versie 0.6.2.
- Leslie, D. A., Teitlebaum, A.D., Anderson, R.J. (1979). <i>Sampling in dollars: een praktische gids voor auditors</i>. Toronto: Copp Clark Pitman.
- Stringer, K. W. (1963) Praktische aspecten van statistische steekproeven bij auditing. <i>Proceedings of Business and Economic Statistics Section</i>, American Statistical Association.
- Touw, P., & Hoogduin, L. (2011). Statistiek voor audit en controlling.

### R-pakketten
---
- jfaPlanning
===

Met de planningsanalyse kan de gebruiker een minimale steekproefomvang berekenen op basis van een reeks steekproefdoelstellingen en samenvattende statistieken van de populatie. Houd er rekening mee dat wanneer u toegang heeft tot de onbewerkte populatiegegevens, u misschien de auditworkflow wilt gebruiken, een analyse die u door het steekproefproces leidt.

<img src="%HELP_FOLDER%/img/workflowPlanning.png" />

Zie de handleiding van de Audit module (download [hier](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) voor meer gedetailleerde informatie over deze analyse.

### Invoer
---

#### Steekproef Doelstellingen
- Uitvoeringsmaterialiteit: ook wel de bovengrens voor de fout, het maximaal toelaatbare afwijkingspercentage of de maximaal toelaatbare afwijking genoemd, de uitvoeringsmaterialiteit is de bovengrens van de toelaatbare afwijking in de te testen populatie. Door te toetsen aan een uitvoeringsmaterialiteit bent u in staat een steekproef te plannen om bewijs te verzamelen voor of tegen de stelling dat de populatie als geheel geen afwijkingen bevat die als materieel worden beschouwd (d.w.z groter zijn dan de uitvoeringsmaterialiteit). U moet deze doelstelling inschakelen als u met een steekproef uit de populatie wilt weten of de populatie een afwijking boven of onder een bepaalde limiet (de prestatiematerialiteit) bevat. Een lagere uitvoeringsmaterialiteit zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere uitvoeringsmaterialiteit resulteren in een lagere steekproefomvang.
- Minimale precisie: de precisie is het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de fout. Door deze steekproefdoelstelling in te schakelen, kunt u een steekproef plannen zodat het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de afwijking tot een minimumpercentage wordt teruggebracht. U moet deze doelstelling inschakelen als u geïnteresseerd bent in het maken van een schatting van de populatieafwijking met een bepaalde nauwkeurigheid. Een lagere minimaal vereiste precisie zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere minimale precisie resulteren in een vereiste steekproefomvang.

#### Populatie
- Aantal eenheden: Het totale aantal eenheden in de populatie. Let op dat de eenheden items (rijen) of monetaire eenheden (waarden) kunnen zijn, afhankelijk van het controlevraagstuk.

#### Betrouwbaarheid
Het gebruikte betrouwbaarheidsniveau. Het betrouwbaarheidsniveau is het complement van het auditrisico: het risico dat de gebruiker bereid is te nemen om een ​​onjuist oordeel over de populatie te geven. Als u bijvoorbeeld een controlerisico van 5% wilt hebben, staat dit gelijk aan 95% betrouwbaarheid.

#### Auditrisicomodel
- Inherent risico: Een categorie of waarschijnlijkheid voor het inherente risico. Inherent risico wordt gedefinieerd als het risico op een afwijking van materieel belang als gevolg van een fout of weglating in een financieel overzicht als gevolg van een andere factor dan een falen van de interne beheersing.
- Beheersingsrisico: Een categorie of waarschijnlijkheid voor het internecontrolerisico. Het interne beheersingsrisico wordt gedefinieerd als het risico van een afwijking van materieel belang in de financiële overzichten die voortvloeit uit het ontbreken of falen van de relevante interne beheersingsmaatregelen van de gecontroleerde.

Wanneer de auditor informatie heeft die wijst op een laag risicoprofiel van de populatie, kan hij deze informatie gebruiken om zijn vereiste steekproefomvang te verkleinen via het Audit Risk Model (ARM), op voorwaarde dat er geen fouten in de populatie zitten. Volgens de ARM is het auditrisico (AR) een functie van het inherente risico (IR), het internecontrolerisico (CR) en het detectierisico (DR).

*AR = IR x CR x DR*

De auditor beoordeelt het inherente risico en het internecontrolerisico in het algemeen op een 3-puntsschaal om het juiste ontdekkingsrisico te bepalen. Met behulp van de ARM en nul fouten hangt de steekproefomvang af van de risicofactor *R* en de uitvoeringsmaterialiteit. De risicofactor *R* is een functie van het detectierisico (Stewart 2012).

*R = -ln(DR)*

De volgende tabel geeft de waarden van *R* weer als functie van het detectierisico, op voorwaarde dat er geen fouten zijn (Touw en Hoogduin 2012).

| Detectierisico (%) | 1 | 4 | 5 | 10 | 14 |
| :---: | :---: | :---: | :---: | :---: | :---: |
| R | 4.6 | 3.2 | 3 | 2.3 | 2 |

De risicofactor *R* kan worden aangepast met behulp van de beoordelingen van het inherente risico en het internecontrolerisico. Standaard is de standaardmethode voor het instellen van de kansen op IR en CR door de onderstaande tabel te volgen voor een detectierisico van 5%:

| | Hoog | Gemiddeld | Laag |
| :---: | :---: | :---: |
| R | 3 | 2 | 1 |

Deze waarden van *R* worden gebruikt om standaardpercentages in te stellen voor IR en CR. De Audit-module verwerkt de volgende standaardwaarden voor IR en CR:

- Hoog: 100%
- Gemiddeld: 60%
- Laag: 36%

U kunt de waarde van IR en CR handmatig aanpassen door de optie Aangepast te selecteren onder de bijbehorende risicobeoordeling, waardoor de risicofactor *R* wordt aangepast.

#### Verwachte fouten in steekproef
De verwachte fouten zijn de toelaatbare fouten die in de steekproef kunnen worden gevonden terwijl de gespecificeerde steekproefdoelstellingen nog steeds worden bereikt. Er wordt een steekproefomvang berekend zodat, wanneer het aantal verwachte fouten in de steekproef wordt gevonden, het gewenste vertrouwen behouden blijft.

*Opmerking:* Het wordt aangeraden om deze waarde conservatief in te stellen om de kans te minimaliseren dat de waargenomen fouten de verwachte fouten overschrijden, wat zou betekenen dat er onvoldoende werk is verricht.

- Relatief: voer uw verwachte fouten in als een percentage ten opzichte van de totale grootte van de selectie.
- Absoluut: Voer uw verwachte fouten in als de som van (proportionele) fouten.

#### Kansverdeling
- Poisson: De Poisson-verdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een kansverdeling die het percentage afwijkingen (*\u03B8*) modelleert als functie van de waargenomen steekproefomvang (*n*) en de som van de proportionele fouten (*t*). Omdat de Poisson-verdeling rekening houdt met gedeeltelijke fouten, wordt deze over het algemeen gebruikt wanneer u een steekproef in een munteenheid plant.
- Binomiaal: De binomiale verdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een kansverdeling die het afwijkingspercentage (*\u03B8*) modelleert als functie van het waargenomen aantal fouten (*k*) en het aantal correcte transacties (*n - k*). Omdat de binominale verdeling strikt geen rekening houdt met gedeeltelijke fouten, wordt deze over het algemeen gebruikt wanneer u geen steekproef in een munteenheid plant.
- Hypergeometrische verdeling: De hypergeometrische verdeling gaat uit van een eindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang klein is. Het is een kansverdeling die het aantal fouten (*K*) in de populatie modelleert als functie van de populatieomvang (*N*), het aantal waargenomen gevonden fouten (*k*) en het aantal correcte transacties ( *N*).

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.

#### Figuren
- Vergelijk steekproefomvang: Produceert een plot die de steekproefomvang vergelijkt 1) over kansverdelingen, en 2) over het aantal verwachte fouten in de steekproef.
- Veronderstelde foutenverdeling: Produceert een plot die de kansverdeling weergeeft die wordt geïmpliceerd door de invoeropties en de berekende steekproefomvang.

#### Tabellen opmaken
- Cijfers: geef tabeluitvoer weer als getallen.
- Percentages: geef de tabeluitvoer weer als percentages.

#### Stapgrootte
Met de stapgrootte kunt u de mogelijke steekproefomvang beperken tot een veelvoud van de waarde. Een stapgrootte van 5 staat bijvoorbeeld alleen steekproefomvang van 5, 10, 15, 20, 25, enz. toe.

### Uitgang
---

#### Planningsoverzicht
- Uitvoeringsmaterialiteit: indien aanwezig, de uitvoeringsmaterialiteit.
- Min. precisie: indien aanwezig, de minimale precisie.
- Verwachte fouten: Het aantal (som van proportionele taints) verwachte / toelaatbare fouten in de steekproef.
- Minimale steekproefomvang: de minimale steekproefomvang.

#### Figuren
- Vergelijk steekproefomvang: Produceert een plot die de steekproefomvang vergelijkt 1) over kansverdelingen, en 2) over het aantal verwachte fouten in de steekproef.
- Veronderstelde foutenverdeling: Produceert een plot die de kansverdeling weergeeft die wordt geïmpliceerd door de invoeropties en de berekende steekproefomvang.

### Referenties
---
- AICPA (2017). <i>Auditgids: controlesteekproeven</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesiaanse en klassieke auditsteekproeven. R-pakket versie 0.6.2.
-Stewart, T. (2012). <i>Technische opmerkingen over de AICPA-auditgids Auditsteekproeven</i>. American Institute of Certified Public Accountants, New York.

### R-pakketten
---
- jfaRepeated Value Analysis
===

This analysis analyzes the frequency with which values get repeated within a dataset (called “number-bunching”) to statistically identify whether the data were likely tampered with. Unlike Benford’s law this approach examines the entire number at once, not only the first or last digit (Simonsohn, 2019).

To determine whether the data show an excessive amount of bunching, the null hypothesis that the data do not contain an unexpected amount of repeated values is tested. To quantify what is expected, this test requires the assumption that the integer portions of the numbers are not associated with their decimal portions.

### Input
---

#### Assignment Box
- Variable: In this box the variable is selected whose digits should be analyzed for repeated values.

#### Tests
- Average frequency: Compute the average frequency of the data.
- Entropy: Compute the entropy of the data.

#### Shuffle Decimal Digits
- Last: Last decimal digit is shuffled.
- Last two: Last two decimal digits are shuffled.
- All: All decimal digits are shuffled.

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.
  - Confidence: The confidence level used in the explanatory text.

#### Tables
- Assumption checks: This table shows the correlation between the integer portions of the numbers and their decimal counterparts. To meet the required assumptions for this procedure, this correlation must be non-existent. This table also displays the correlation between the samples of the two simulation runs (average frequency and entropy).
- Frequency table: Produces a table containing the count and the percentage for every unique value in the data set.

#### Plots
- Observed vs. expected: Produces a histogram of the expected average frequencies and / or entropy vs. the observed average frequency and / or entropy.
- Histogram: Produces a histogram with a single bin for each observed value.

#### Advanced Options
- Number of samples: The number of samples to use for simulating the p value.
- Seed: Selects the seed for the random number generator in order to reproduce results.

### Output
---

#### Repeated Values Test
- n: The number of observations in the data.
- Frequency: The average frequency with which numbers are repeated in the data. The formula for the average frequency is *AF = &#8721; f&#7522;&#178; / &#8721; f&#7522;* where f&#7522; is the frequency of each unique value *i* in the data set.
- Entropy: The entropy is the average level of information inherent in the variable's outcomes. The entropy is calculated as *S = - &#8721; (p&#7522; &#215; log(p&#7522;))* where p&#7522; is the proportion of observations with each value (so *p&#7522; = f&#7522; / N*).

#### Assumption Checks
- n: Sample size.
- r: Pearson correlation coefficient.
- t: t-value.
- df: Degrees of freedom.
- p: p-value.

#### Frequency Table
- Value: The value in the row.
- Count: The number of times the value is observed.
- Percentage: The percentage of times the value is observed.

#### Plots
- Observed vs. expected: Displays the observed vs. the simulated value(s).
- Histogram: Displays a histogram with a single bin for each observed value.

### References
---
- Derks, K (2021). digitTests: Tests for Detecting Irregular Digit Patterns. R package version 0.1.0.
- Simohnsohn, U. (2019, May 25). Number-Bunching: A New Tool for Forensic Data Analysis. Retrieved from [http://datacolada.org/77](http://datacolada.org/77).

### R Packages
---
- digitTests
Selection
===

The selection analysis allows the user to select a number of sampling units from a population using a combination of sampling techniques (record sampling versus monetary unit sampling) and sampling methods (random sampling, cell sampling, fixed interval sampling) that are standard in an audit context.

<img src="%HELP_FOLDER%/img/workflowSelection.png" />

Please see the manual of the Audit module (download [here](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) for more detailed information about this analysis.

### Input
---

#### Assignment Box
- Item ID: A unique non-missing identifier for every item in the population. The row number of the items is sufficient.
- Book Values: The variable that contains the book values of the items in the population.
- Ranking Variable: When provided, the population is first ranked in ascending order with respect to the values of this variable.
- Additional Variables: Any other variables that should be included in the sample.

#### Sample Size
The required number of sampling units that should be selected from the population. Be aware that the sampling units are determined by the *units* option. By default, when no book values are provided, the sampling units are items (rows). When book values are provided, the ideal sampling units to use are monetary units.

#### Sampling Units
- Items: Performs selection using the items in the population as sampling units.
- Monetary units: Performs selection using the monetary units in the population as sampling units. This method is preferred when you want to include more items with a high value in the sample.

#### Method
- Fixed interval sampling: Performs selection by dividing the population in equal intervals and selecting a fixed unit in each interval. Any item with a value larger than the interval will always be included in the sample.
  - Starting point: Selects which sampling unit is selected from each interval.
- Cell sampling: Performs selection by dividing the population in equal intervals and selecting a variable unit in each interval. Any item with a value larger than twice the interval will always be included in the sample.
  - Seed: Selects the seed for the random number generator in order to reproduce results.
- Random sampling: Performs random selection in which each sampling unit has an equal chance of being selected.
  - Seed: Selects the seed for the random number generator in order to reproduce results.

#### Randomize Item Order
Randomizes the items in the population before selection is performed.

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.

#### Column Name Selection Result
When a name is provided, adds the result from the selection analysis in a new column to the data. The new column reflects how many times each transaction is included in the sample.

#### Tables
- Descriptive statistics: Produces a table containing descriptive information about numerical variables in the selection. Statistics that are included are the mean, the median, the standard deviation, the variance, the minimum, the maximum, and the range.
- Selected items: Produces a table containing the selected transactions along with any additional observations provided in the additional variables field.

### Output
---

#### Selection Summary
- No. units: The number of selected sampling units from the population.
- No. items: The number of selected items from the population.
- Selection value: The total value of the selected items. Only displayed when monetary unit sampling is used.
- % of population size / value: The selected proportion of the total size or value of the population.

#### Information about Monetary Interval Selection
- Items: The number of items in the population.
- Value: The value of the items in the population.
- Selected items: The number of items in the sample.
- Selected units: The number of selected units from the population.
- Selection value: The value of the items in the sample.
- % of total value: The selected proportion of the total value of the items compared to the items in the population.

#### Descriptive Statistics
- Valid: Number of valid cases.
- Mean: Arithmetic mean of the data points.
- Median: Median of the data points.
- Std. deviation: Standard deviation of the data points.
- Variance: Variance of the data points.
- Range: Range of the data points.
- Minimum: Minimum of the data points.
- Maximum: Maximum of the data points.

#### Selected Items
- Row: The row number of the item.
- Selected: The number of times (a unit in) the item is selected.

### References
---
- AICPA (2017). <i>Audit Guide: Audit Sampling</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesian and Classical Audit Sampling. R package version 0.6.2.

### R Packages
---
- jfa
Bayesian Sampling Workflow
===

The task of an auditor is to make a judgment regarding the fairness of the presented transactions in a population. When the auditor has access to the raw population data, they can use the *audit workflow* to calculate how many samples need to be evaluated in order to meet a certain confidence in their judgment. The user can then sample these items from the population, inspect and audit these items, and perform statistical inference about the misstatement in the population. The sampling workflow guides the auditor through the audit process, making the correct choices of calculations along the way.

Please see the manual of the Audit module (download [here](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) for more detailed information about this analysis.

### Workflow
---

- Planning: Calculate the minimum sample size to achieve your sampling objectives with the specified confidence.
- Selection: Select the required sampling units from the population.
- Execution: Annotate the selection with your assessment of the fairness of the selected items.
- Evaluation: Make a population statement based on your annotated selection.

<img src="%HELP_FOLDER%/img/workflow.png" />

### Input - Planning
---

#### Sampling Objectives
- Performance materiality: Also called the upper error limit, the tolerable deviation rate, or the tolerable misstatement, the performance materiality is the upper bound of tolerable misstatement in the population to be tested. By testing against a performance materiality, you are able to plan a sample in order to collect evidence for or against the statement that the population as a whole does not contain misstatements that are considered material (i.e., are greater than the upper bound of tolerable misstatement). You should enable this objective when you want to find out whether the population contains misstatement above or below a certain limit (the performance materiality) using a sample of the population. A lower performance materiality will result in a higher required sample size. Vice versa, a higher performance materiality will result in a lower required sample size.
- Minimum precision: The precision is the the difference between the estimated most likely error and the upper bound on the misstatement. By enabling this sampling objective, you are be able to plan a sample so that the difference between the estimated most likely error and the upper bound on the misstatement is reduced to a minimum percentage. You should enable this objective if you are interested in making an estimate of the population misstatement with a certain accuracy. A lower minimum required precision will result in a higher required sample size. Vice versa, a higher minimum required precision will result in a lower required sample size.

#### Confidence
The confidence level used. The confidence level is the complement of the audit risk: the risk that the user is willing to take to give an incorrect judgment about the population. For example, if you want to have an audit risk of 5%, this equals 95% confidence.

#### Assignment Box
- Item ID: A unique non-missing identifier for every item in the population. The row number of the items is sufficient.
- Book Values: The variable that contains the book values of the items in the population.

#### Expected errors in Sample
The expected errors are the tolerable errors that can be found in the sample while still achieving the specified sampling objectives. A sample size is calculated so that, when the number of expected errors is found in the sample, the desired confidence is retained.

*Note:* It is advised to set this value conservatively to minimize the probability of the observed errors exceeding the expected errors, which would imply that insufficient work has been done.

- Relative: Enter your expected errors as a percentage relative to the total size of the selection.
- Absolute: Enter your expected errors as the sum of (proportional) errors.

#### Probability Distribution
- Beta-binomial: The beta-binomial distribution accompanies the hypergeometric likelihood (Dyer & Pierce, 1993). The hypergeometric likelihood assumes a finite population size and is therefore generally used when the population size is small. It is a likelihood that models the number of errors (*K*) in the population as a function of the population size (*N*), the number of observed found errors (*k*) and the number of correct transactions (*n*).
- Beta: The beta distribution accompanies the binomial likelihood. The binomial likelihood assumes an infinite population size and is therefore generally used when the population size is large. It is a likelihood that models the rate of misstatement (*\u03B8*) as a function of the observed number of errors (*k*) and the number of correct transactions (*n - k*). Because the binomial distribution strictly does not accommodate partial errors, it is generally used when you are not planning a monetary unit sample. However, the beta distribution does accommodate partial errors, and may also be used for monetary unit sampling (de Swart, Wille & Majoor, 2013).
- Gamma: The gamma distribution accompanies the Poisson likelihood. The Poisson likelihood assumes an infinite population size and is therefore generally used when the population size is large. It is a likelihood that models the rate of misstatement (*\u03B8*) as a function of the observed sample size (*n*) and the sum of the proportional errors found (*t*). Because the gamma distribution accommodates partial errors it is generally used when you are planning a monetary unit sample (Stewart, 2013).

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.

#### Tables
- Descriptive statistics: Produces a table with descriptive statistics of the book values in the population.
- Equivalent prior sample: Produces a table that displays the implicit sample on which the prior distribution is based.
- Prior and posterior: Produces a table in which the prior and expected posterior distribution are summarized through several statistics, such as their functional form, their prior and expected posterior probabilities and odds, and the shift between these.

#### Plots
- Prior and posterior: Produces a plot that shows the prior distribution and the posterior distribution after observing the intended sample.
  - Additional info: Produces dots on the materiality.
- Prior predictive: Produces a plot of the predictions of the prior distribution.
- Compare sample sizes: Produces a plot that compares the sample size 1) across probability distributions, and 2) across the number of expected errors in the sample.
- Distribution of book values: Produces a histogram of the book values in the population.

#### Prior
- Default: This option does not incorporate any information into the statistical analysis and therefore assumes a negligible and conservative prior distribution.
- Manual: Provide the parameters of the prior distribution.
- Earlier sample: Create a prior distribution on the basis of an earlier sample.
  - Size: Earlier sample size.
  - Errors: Earlier found errors.
- Impartial: Create a prior distribution that is impartial with respect to the tested hypotheses.
- Risk assessments: Translate information from the audit risk model into a prior distribution.
  - Inherent risk: A category or probability for the inherent risk. Inherent risk is defined as the risk of material misstatement posed by an error or omission in a financial statement due to a factor other than a failure of internal control.
  - Control risk: A category or probability for the internal control risk. Control risk is defined as the risk of a material misstatement in the financial statements arising due to absence or failure in the operation of relevant controls of the auditee.

#### Critical Items
- Negative book values: Isolates negative book values from the population.
  - Keep: Keeps negative book values to be inspected in the sample.
  - Remove: Removes negative book values.

#### Format Tables
- Numbers: Display table output as numbers.
- Percentages: Display table output as percentages.

#### Iterations
- Increment: The increment alows you to limit the possible sample sizes to a multiple of its value. For example, an increment of 5 allows only sample sizes of 5, 10, 15, 20, 25, etc.
- Maximum: The maximum allows you to limit the sample size with a maximum.

#### Assume Homogeneous Taints
Clicking this box will allow you to separate the known and the unknown misstatement in the population to be more efficient. Note that this requires the assumption that the taints in the sample are representative of the taints in the unseen part of the population.

### Ouput - Planning
---

#### Planning Summary
- Performance materiality: When provided, the performance materiality.
- Min. precision: When provided, the minimum precision.
- Expected errors: The number (sum of proportional taints) of expected / tolerable errors in the sample.
- Minimum sample size: The minimum sample size.

#### Descriptive Statistics
- Population size: Number of items in the population.
- Value: Total value of the book values.
- Absolute value: Absolute value of the book values.
- Mean: Mean of the book values.
- Std. deviation: Standard deviation of the book values.
- Quartile: Quartiles of the book values.

#### Equivalent prior sample
- Equivalent sample size: The sample size equivalent to the prior information.
- Equivalent errors: The number of errors equivalent to the prior information.

#### Prior and Posterior
- Functional form: The functional form of the distribution.
- Support H-: Total probability in the range of H- under the distribution. Only displayes when testing against a performance materiality.
- Support H+: Total probability in the range of H+ under the distribution. Only displayes when testing against a performance materiality.
- Ratio H- / H+: Odds in favor of H- under the distribution. Only displayes when testing against a performance materiality.
- Mean: Mean of the distribution.
- Median: Median of the distribution.
- Mode: Mode of the distribution.
- Upper bound: x-% percentile of the distribution.
- Precision: Difference between the upper bound and the mode of the distribution.

#### Plots
- Prior and posterior: Produces a plot that shows the prior distribution and the posterior distribution after observing the intended sample.
  - Additional info: Produces dots on the materiality.
- Prior predictive: Produces a plot of the predictions of the prior distribution.
- Compare sample sizes: Produces a plot that compares the sample size 1) across probability distributions, and 2) across the number of expected errors in the sample.
- Distribution of book values: Produces a histogram of the book values in the population.

### Input - Selection
---

#### Assignment Box
- Ranking Variable: When provided, the population is first ranked in ascending order with respect to the values of this variable.
- Additional Variables: Any other variables that should be included in the sample.

#### Sample Size
The required number of sampling units that should be selected from the population. Be aware that the sampling units are determined by the *units* option. By default, when no book values are provided, the sampling units are items (rows). When book values are provided, the ideal sampling units to use are monetary units.

#### Sampling Units
- Items: Performs selection using the items in the population as sampling units.
- Monetary units: Performs selection using the monetary units in the population as sampling units. This method is preferred when you want to include more items with a high value in the sample.

#### Method
- Fixed interval sampling: Performs selection by dividing the population in equal intervals and selecting a fixed unit in each interval. Any item with a value larger than the interval will always be included in the sample.
  - Starting point: Selects which sampling unit is selected from each interval.
- Cell sampling: Performs selection by dividing the population in equal intervals and selecting a variable unit in each interval. Any item with a value larger than twice the interval will always be included in the sample.
  - Seed: Selects the seed for the random number generator in order to reproduce results.
- Random sampling: Performs random selection in which each sampling unit has an equal chance of being selected.
  - Seed: Selects the seed for the random number generator in order to reproduce results.

#### Randomize Item Order
Randomizes the items in the population before selection is performed.

#### Tables
- Descriptive statistics: Produces a table containing descriptive information about numerical variables in the selection. Statistics that are included are the mean, the median, the standard deviation, the variance, the minimum, the maximum, and the range.
- Selected items: Produces a table containing the selected transactions along with any additional observations provided in the additional variables field.

### Output - Selection
---

#### Selection Summary
- No. units: The number of selected sampling units from the population.
- No. items: The number of selected items from the population.
- Selection value: The total value of the selected items. Only displayed when monetary unit sampling is used.
- % of population size / value: The selected proportion of the total size or value of the population.

#### Information about Monetary Interval Selection
- Items: The number of items in the population.
- Value: The value of the items in the population.
- Selected items: The number of items in the sample.
- Selected units: The number of selected units from the population.
- Selection value: The value of the items in the sample.
- % of total value: The selected proportion of the total value of the items compared to the items in the population.

#### Descriptive Statistics
- Valid: Number of valid cases.
- Mean: Arithmetic mean of the data points.
- Median: Median of the data points.
- Std. deviation: Standard deviation of the data points.
- Variance: Variance of the data points.
- Range: Range of the data points.
- Minimum: Minimum of the data points.
- Maximum: Maximum of the data points.

#### Selected Items
- Row: The row number of the item.
- Selected: The number of times (a unit in) the item is selected.

### Input - Execution
---

#### Annotation
- Audit value: Annotate the items in the selection with their audit (true) values. This approach is recommended (and automatically selected) when the items have a monetary value.
- Correct / Incorrect: Annotate the items in the selection with correct (0) or incorrect (1). This approach is recommended (and automatically selected) when your items do not have a monetary value.

### Input - Evaluation
---

#### Assignment Box
- Audit result / values: The variable that contains the audit (true) values, or the binary classification of correct (0) or incorrect (1).

#### Method
See *Probability Distribution*.

#### Area Under Posterior
- One-sided upper bound: Gives an (upper) estimate of the misstatement in the population.
- Two-sided interval: Gives an (upper and lower) estimate of the misstatement in the population.

#### Tables
- Corrections to population: Produces a table that contains the required corrections to the population value to achieve the sampling objectives.
- Prior and posterior: Produces a table in which the prior and expected posterior distribution are summarized through several statistics, such as their functional form, their prior and expected posterior probabilities and odds, and the shift between these.
- Assumption checks: Produces a table that displays the correlation between the book values in the sample and their taints.
  - Confidence interval: Width of the confidence interval for the correlation.

#### Plots
- Prior and posterior: Produces a plot that shows the prior distribution and the posterior distribution after observing the intended sample.
  - Additional info: Produces dots on the materiality.
- Posterior predictive: Produces a plot of the predictions of the posterior distribution.
- Sampling objectives: Produces a bar chart comparing the materiality, upper bound on the misstatement and most likely error (MLE).
- Scatter plot: Produces a scatter plot comparing book values of the selection against their audit values. Observations that are in error are colored in red.
  - Display correlation: Adds the correlation between the book values and the audit values to the plot.
  - Display item ID's: Adds the item ID's to the plot.

### Output - Evaluation
---

#### Evaluation summary
- Materiality: When provided, the performance materiality.
- Min. precision: When provided, the minimum precision.
- Sample size: The sample size (number of units).
- Errors: The number of erroneous elements in the selection.
- Taint: The sum of the proportional errors. Audited items can be evaluated while incorporating the magnitude of the misstatement by calculating their taints. The taint of an item *i* is the proportional difference between that item's book value (*y*) and the item's audit (true) value (*x*). Positive taints are associated with overstatements, while negative taints occur when items are understated.
<img src="%HELP_FOLDER%/img/taints.png" />
- Most likely error: The most likely error in the population.
- x-% Confidence bound: Upper bound on the misstatement in the population.
- Precision: Difference between upper bound and most likely error.
- BF-+: The Bayes factor for the test.

#### Prior and Posterior
- Functional form: The functional form of the distribution.
- Support H-: Total probability in the range of H- under the distribution. Only displayes when testing against a performance materiality.
- Support H+: Total probability in the range of H+ under the distribution. Only displayes when testing against a performance materiality.
- Ratio H- / H+: Odds in favor of H- under the distribution. Only displayes when testing against a performance materiality.
- Mean: Mean of the distribution.
- Median: Median of the distribution.
- Mode: Mode of the distribution.
- Upper bound: x-% percentile of the distribution.
- Precision: Difference between the upper bound and the mode of the distribution.

#### Corrections to Population
- Correction: The amount or percentage to be deducted from the population.

#### Assumption Checks
- n: Sample size.
- Pearsons r: Pearson correlation coefficient.
- x-% upper bound: Upper bound for correlation coefficient.
- p: p-value for the test.
- BF-0: Bayes factor for the test.

#### Plots
- Prior and posterior: Produces a plot that shows the prior distribution and the posterior distribution after observing the intended sample.
  - Additional info: Produces dots on the materiality.
- Posterior predictive: Produces a plot of the predictions of the posterior distribution.
- Sampling objectives: Produces a bar chart comparing the materiality, upper bound on the misstatement and most likely error (MLE).
- Scatter plot: Produces a scatter plot comparing book values of the selection against their audit values. Observations that are in error are colored in red.
  - Display correlation: Adds the correlation between the book values and the audit values to the plot.
  - Display item ID's: Adds the item ID's to the plot.

### References
---
- AICPA (2017). <i>Audit Guide: Audit Sampling</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesian and Classical Audit Sampling. R package version 0.6.2.
- Dyer, D., & Pierce, R. L. (1993). On the choice of the prior distribution in hypergeometric sampling. <i>Communications in Statistics-Theory and Methods</i>, 22(8), 2125-2146.
- Stewart, T. R. (2013). A Bayesian audit assurance model with application to the component materiality problem in group audits (Doctoral dissertation).
- de Swart, J., Wille, J., & Majoor, B. (2013). Het 'Push Left'-Principe als Motor van Data Analytics in de Accountantscontrole [The 'Push-Left'-Principle as a Driver of Data Analytics in Financial Audit]. <i>Maandblad voor Accountancy en Bedrijfseconomie</i>, 87, 425-432.

### R Packages
---
- jfa
Bayesian Evaluation
===

The Bayesian evaluation analysis allows the user to perform inference about the total misstatement in the population on the basis of an audit sample.

<img src="%HELP_FOLDER%/img/workflowEvaluation.png" />

Please see the manual of the Audit module (download [here](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) for more detailed information about this analysis.

### Input
---

#### Sampling Objectives
- Performance materiality: Also called the upper error limit, the tolerable deviation rate, or the tolerable misstatement, the performance materiality is the upper bound of tolerable misstatement in the population to be tested. By testing against a performance materiality, you are able to plan a sample in order to collect evidence for or against the statement that the population as a whole does not contain misstatements that are considered material (i.e., are greater than the upper bound of tolerable misstatement). You should enable this objective when you want to find out whether the population contains misstatement above or below a certain limit (the performance materiality) using a sample of the population. A lower performance materiality will result in a higher required sample size. Vice versa, a higher performance materiality will result in a lower required sample size.
- Minimum precision: The precision is the the difference between the estimated most likely error and the upper bound on the misstatement. By enabling this sampling objective, you are be able to plan a sample so that the difference between the estimated most likely error and the upper bound on the misstatement is reduced to a minimum percentage. You should enable this objective if you are interested in making an estimate of the population misstatement with a certain accuracy. A lower minimum required precision will result in a higher required sample size. Vice versa, a higher minimum required precision will result in a lower required sample size.

#### Population
- No. items: The total number of items (rows) in the population.
- No. units: The total number of units in the population. Note that the units can be items (rows) or monetary units (values) depending on the audit question.

#### Confidence
The confidence level used. The confidence level is the complement of the audit risk: the risk that the user is willing to take to give an incorrect judgment about the population. For example, if you want to have an audit risk of 5%, this equals 95% confidence.

#### Assignment Box
- Item ID: A unique non-missing identifier for every item in the population. The row number of the items is sufficient.
- Book Values: The variable that contains the book values of the items in the population.
- Audit result / values: The variable that contains the audit (true) values, or the binary classification of correct (0) or incorrect (1).
- Selection counter: The variable that contains how many times each observation should be evaluated.

#### Data
- Raw: Use raw data.
- Summary statistics: Use summary statistics.

#### Probability Distribution
- Gamma: The gamma distribution accompanies the Poisson likelihood. The Poisson likelihood assumes an infinite population size and is therefore generally used when the population size is large. It is a likelihood that models the rate of misstatement (*\u03B8*) as a function of the observed sample size (*n*) and the sum of the proportional errors found (*t*). Because the gamma distribution accommodates partial errors it is generally used when you are planning a monetary unit sample (Stewart, 2013).
- Beta: The beta distribution accompanies the binomial likelihood. The binomial likelihood assumes an infinite population size and is therefore generally used when the population size is large. It is a likelihood that models the rate of misstatement (*\u03B8*) as a function of the observed number of errors (*k*) and the number of correct transactions (*n - k*). Because the binomial distribution strictly does not accommodate partial errors, it is generally used when you are not planning a monetary unit sample. However, the beta distribution does accommodate partial errors, and may also be used for monetary unit sampling (de Swart, Wille & Majoor, 2013).
- Beta-binomial: The beta-binomial distribution accompanies the hypergeometric likelihood (Dyer & Pierce, 1993). The hypergeometric likelihood assumes a finite population size and is therefore generally used when the population size is small. It is a likelihood that models the number of errors (*K*) in the population as a function of the population size (*N*), the number of observed found errors (*k*) and the number of correct transactions (*n*).

#### Area Under Posterior
- One-sided upper bound: Gives an (upper) estimate of the misstatement in the population.
- Two-sided interval: Gives an (upper and lower) estimate of the misstatement in the population.

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.

#### Tables
- Prior and posterior: Produces a table in which the prior and expected posterior distribution are summarized through several statistics, such as their functional form, their prior and expected posterior probabilities and odds, and the shift between these.
- Corrections to population: Produces a table that contains the required corrections to the population value to achieve the sampling objectives.
- Assumption checks: Produces a table that displays the correlation between the book values in the sample and their taints.
  - Confidence interval: Width of the confidence interval for the correlation.

#### Plots
- Prior and posterior: Produces a plot that shows the prior distribution and the posterior distribution after observing the intended sample.
  - Additional info: Produces dots on the materiality.
- Posterior predictive: Produces a plot of the predictions of the posterior distribution.
- Sampling objectives: Produces a bar chart comparing the materiality, maximum misstatement and most likely error (MLE).
- Scatter plot: Produces a scatter plot comparing book values of the selection against their audit values. Observations that are in error are colored in red.
  - Display correlation: Adds the correlation between the book values and the audit values to the plot.
  - Display item ID's: Adds the item ID's to the plot.

#### Prior
- Default: This option does not incorporate any information into the statistical analysis and therefore assumes a negligible and conservative prior distribution.
- Manual: Provide the parameters of the prior distribution.
- Earlier sample: Create a prior distribution on the basis of an earlier sample.
  - Size: Earlier sample size.
  - Errors: Earlier found errors.
- Impartial: Create a prior distribution that is impartial with respect to the tested hypotheses.
- Risk assessments: Translate information from the audit risk model into a prior distribution.
  - Inherent risk: A category or probability for the inherent risk. Inherent risk is defined as the risk of material misstatement posed by an error or omission in a financial statement due to a factor other than a failure of internal control.
  - Control risk: A category or probability for the internal control risk. Control risk is defined as the risk of a material misstatement in the financial statements arising due to absence or failure in the operation of relevant controls of the auditee.

#### Expected errors in Sample
The expected errors are the tolerable errors that can be found in the sample while still achieving the specified sampling objectives. A sample size is calculated so that, when the number of expected errors is found in the sample, the desired confidence is retained.

*Note:* It is advised to set this value conservatively to minimize the probability of the observed errors exceeding the expected errors, which would imply that insufficient work has been done.

- Relative: Enter your expected errors as a percentage relative to the total size of the selection.
- Absolute: Enter your expected errors as the sum of (proportional) errors.

#### Critical Items
- Negative book values: Isolates negative book values from the population.
  - Keep: Keeps negative book values to be inspected in the sample.
  - Remove: Removes negative book values.

#### Format Tables
- Numbers: Display table output as numbers.
- Percentages: Display table output as percentages.
- Monetary values: Display table output as monetary values.

#### Assume Homogeneous Taints
Clicking this box will allow you to separate the known and the unknown misstatement in the population to be more efficient. Note that this requires the assumption that the taints in the sample are representative of the taints in the unseen part of the population.

### Output
---

#### Evaluation summary
- Materiality: When provided, the performance materiality.
- Min. precision: When provided, the minimum precision.
- Sample size: The sample size (number of units).
- Errors: The number of erroneous elements in the selection.
- Taint: The sum of the proportional errors. Audited items can be evaluated while incorporating the magnitude of the misstatement by calculating their taints. The taint of an item *i* is the proportional difference between that item's book value (*y*) and the item's audit (true) value (*x*). Positive taints are associated with overstatements, while negative taints occur when items are understated.
<img src="%HELP_FOLDER%/img/taints.png" />
- Most likely error: The most likely error in the population.
- x-% Confidence bound: Upper bound on the misstatement in the population.
- Precision: Difference between upper bound and most likely error.
- BF-+: The Bayes factor for the test.

#### Prior and Posterior
- Functional form: The functional form of the distribution.
- Support H-: Total probability in the range of H- under the distribution. Only displayes when testing against a performance materiality.
- Support H+: Total probability in the range of H+ under the distribution. Only displayes when testing against a performance materiality.
- Ratio H- / H+: Odds in favor of H- under the distribution. Only displayes when testing against a performance materiality.
- Mean: Mean of the distribution.
- Median: Median of the distribution.
- Mode: Mode of the distribution.
- Upper bound: x-% percentile of the distribution.
- Precision: Difference between the upper bound and the mode of the distribution.

#### Corrections to Population
- Correction: The amount or percentage to be deducted from the population.

#### Assumption Checks
- n: Sample size.
- Pearsons r: Pearson correlation coefficient.
- x-% upper bound: Upper bound for correlation coefficient.
- p: p-value for the test.
- BF-0: Bayes factor for the test.

#### Plots
- Prior and posterior: Produces a plot that shows the prior distribution and the posterior distribution after observing the intended sample.
  - Additional info: Produces dots on the materiality.
- Posterior predictive: Produces a plot of the predictions of the posterior distribution.
- Sampling objectives: Produces a bar chart comparing the materiality, maximum misstatement and most likely error (MLE).
- Scatter plot: Produces a scatter plot comparing book values of the selection against their audit values. Observations that are in error are colored in red.
  - Display correlation: Adds the correlation between the book values and the audit values to the plot.
  - Display item ID's: Adds the item ID's to the plot.

### References
---
- AICPA (2017). <i>Audit Guide: Audit Sampling</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesian and Classical Audit Sampling. R package version 0.6.2.
- Dyer, D., & Pierce, R. L. (1993). On the choice of the prior distribution in hypergeometric sampling. <i>Communications in Statistics-Theory and Methods</i>, 22(8), 2125-2146.
- Stewart, T. R. (2013). A Bayesian audit assurance model with application to the component materiality problem in group audits (Doctoral dissertation).
- de Swart, J., Wille, J., & Majoor, B. (2013). Het 'Push Left'-Principe als Motor van Data Analytics in de Accountantscontrole [The 'Push-Left'-Principle as a Driver of Data Analytics in Financial Audit]. <i>Maandblad voor Accountancy en Bedrijfseconomie</i>, 87, 425-432.

### R Packages
---
- jfa
Bayesiaanse Evaluatie
===

De Bayesiaanse evaluatie-analyse stelt de gebruiker in staat om op basis van een steekproef conclusies te trekken over de totale fout in de populatie.

<img src="%HELP_FOLDER%/img/workflowEvaluation.png" />

Zie de handleiding van de Audit module (download [hier](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) voor meer gedetailleerde informatie over deze analyse.

### Invoer
---

#### Steekproef Doelstellingen
- Uitvoeringsmaterialiteit: ook wel de bovengrens voor de fout, het maximaal toelaatbare afwijkingspercentage of de maximaal toelaatbare afwijking genoemd, de uitvoeringsmaterialiteit is de bovengrens van de toelaatbare afwijking in de te testen populatie. Door te toetsen aan een uitvoeringsmaterialiteit bent u in staat een steekproef te plannen om bewijs te verzamelen voor of tegen de stelling dat de populatie als geheel geen afwijkingen bevat die als materieel worden beschouwd (d.w.z groter zijn dan de uitvoeringsmaterialiteit). U moet deze doelstelling inschakelen als u met een steekproef uit de populatie wilt weten of de populatie een afwijking boven of onder een bepaalde limiet (de prestatiematerialiteit) bevat. Een lagere uitvoeringsmaterialiteit zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere uitvoeringsmaterialiteit resulteren in een lagere steekproefomvang.
- Minimale precisie: de precisie is het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de fout. Door deze steekproefdoelstelling in te schakelen, kunt u een steekproef plannen zodat het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de afwijking tot een minimumpercentage wordt teruggebracht. U moet deze doelstelling inschakelen als u geïnteresseerd bent in het maken van een schatting van de populatieafwijking met een bepaalde nauwkeurigheid. Een lagere minimaal vereiste precisie zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere minimale precisie resulteren in een vereiste steekproefomvang.

#### Populatie
- Aantal items: Het totale aantal items (rijen) in de populatie.
- Aantal eenheden: Het totale aantal eenheden in de populatie. Let op dat de eenheden items (rijen) of monetaire eenheden (waarden) kunnen zijn, afhankelijk van het controlevraagstuk.

#### Betrouwbaarheid
Het gebruikte betrouwbaarheidsniveau. Het betrouwbaarheidsniveau is het complement van het auditrisico: het risico dat de gebruiker bereid is te nemen om een ​​onjuist oordeel over de populatie te geven. Als u bijvoorbeeld een controlerisico van 5% wilt hebben, staat dit gelijk aan 95% betrouwbaarheid.

#### Opdrachtbox
- Item-ID: een unieke niet-ontbrekende identifier voor elk item in de populatie. Het rijnummer van de items is voldoende.
- Boekwaarden: de variabele die de boekwaarden van de items in de populatie bevat.
- Auditresultaat / waarden: De variabele die de audit (true) waarden bevat, of de binaire classificatie van correct (0) of incorrect (1).
- Selectieteller: De variabele die aangeeft hoe vaak elke waarneming moet worden geëvalueerd.

#### Gegevens
- Raw: gebruik onbewerkte gegevens.
- Overzichtsstatistieken: gebruik overzichtsstatistieken.

#### Kansverdeling
- Gamma: De gammaverdeling gaat samen met de Poisson-verdeling. De Poisson-verdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een verdeling waarmee het percentage afwijkingen (*\u03B8*) modelleerd als functie van de waargenomen steekproefomvang (*n*) en de som van de gevonden proportionele fouten (*t*). Omdat de gammaverdeling deelfouten mogelijk maakt, wordt deze over het algemeen gebruikt bij het plannen van een steekproef in munteenheden (Stewart, 2013).
- Beta: de bètaverdeling gaat samen met de binominaalverdeling. De binominaalverdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een verdeling waarmee het percentage afwijkingen (*\u03B8*) wordt gemodelleerd als functie van het waargenomen aantal fouten (*k*) en het aantal correcte transacties (*n - k*). Omdat de binominaalverdeling strikt genomen geen rekening houdt met gedeeltelijke fouten, wordt deze over het algemeen gebruikt wanneer u geen steekproef in een munteenheid plant. De betaverdeling is echter geschikt voor gedeeltelijke fouten en kan ook worden gebruikt voor steekproeven op monetaire eenheden (de Swart, Wille & Majoor, 2013).
- Beta-binomiaal: De beta-binomiaalverdeling gaat samen met de hypergeometrische verdeling (Dyer & Pierce, 1993). De hypergeometrische verdeling gaat uit van een eindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang klein is. Het is een verdeling waarmee het aantal fouten (*K*) in de populatie wordt gemodelleerd als functie van de populatieomvang (*N*), het aantal geobserveerde gevonden fouten (*k*) en het aantal correcte transacties (*N*).

#### Gebied onder Posterior
- Eenzijdige bovengrens: Geeft een (bovenste) schatting van de afwijking in de populatie.
- Tweezijdig interval: Geeft een (bovenste en onderste) schatting van de afwijking in de populatie.

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.

#### Tabellen
- Prior en posterior: Produceert een tabel waarin de eerdere en verwachte posterieure distributie worden samengevat door middel van verschillende statistieken, zoals hun functionele vorm, hun prior en posterior kansen, en de verschuiving daartussen.
- Correcties op populatie: Produceert een tabel die de vereiste correcties op de populatiewaarde bevat om de steekproefdoelstellingen te bereiken.
- Aannamecontroles: Produceert een tabel die de correlatie weergeeft tussen de boekwaarden in de steekproef en hun taints.
  - Betrouwbaarheidsinterval: Breedte van het betrouwbaarheidsinterval voor de correlatie.

#### Figuren
- Prior en posterior: Produceert een plot die de prior verdeling en de posterieure distributie toont na observatie van het beoogde monster.
  - Extra info: Produceert stippen op de materialiteit.
- Posterior predictive: Produceert een plot van de voorspellingen van de posterieure verdeling.
- Steekproefdoelstellingen: Produceert een staafdiagram waarin de materialiteit, de maximale afwijking en de meest waarschijnlijke fout (MLE) worden vergeleken.
- Spreidingsplot: Produceert een spreidingsplot die boekwaarden van de selectie vergelijkt met hun controlewaarden. Waarnemingen die fout zijn, zijn rood gekleurd.
  - Correlatie weergeven: Voegt de correlatie tussen de boekwaarden en de controlewaarden toe aan de plot.
  - Item-ID's weergeven: voegt de item-ID's toe aan de plot.

#### Prior
- Standaard: deze optie neemt geen informatie op in de statistische analyse en gaat daarom uit van een verwaarloosbare en conservatieve eerdere verdeling.
- Handmatig: geef de parameters van de prior verdeling op.
- Eerdere steekproef: Maak een eerdere verdeling op basis van een eerdere steekproef.
  - Grootte: eerdere steekproefomvang.
  - Fouten: eerder gevonden fouten.
- Onpartijdig: maak een voorafgaande verdeling die onpartijdig is met betrekking tot de geteste hypothesen.
- Risicobeoordelingen: Vertaal informatie uit het auditrisicomodel naar een eerdere verspreiding.
  - Inherent risico: Een categorie of waarschijnlijkheid voor het inherente risico. Inherent risico wordt gedefinieerd als het risico op een afwijking van materieel belang als gevolg van een fout of weglating in een financieel overzicht als gevolg van een andere factor dan een falen van de interne beheersing.
  - Beheersingsrisico: Een categorie of waarschijnlijkheid voor het internecontrolerisico. Het interne beheersingsrisico wordt gedefinieerd als het risico van een afwijking van materieel belang in de financiële overzichten die voortvloeit uit het ontbreken of falen van de relevante interne beheersingsmaatregelen van de gecontroleerde.

#### Verwachte fouten in steekproef
De verwachte fouten zijn de toelaatbare fouten die in de steekproef kunnen worden gevonden terwijl de gespecificeerde steekproefdoelstellingen nog steeds worden bereikt. Er wordt een steekproefomvang berekend zodat, wanneer het aantal verwachte fouten in de steekproef wordt gevonden, het gewenste vertrouwen behouden blijft.

*Opmerking:* Het wordt aangeraden om deze waarde conservatief in te stellen om de kans te minimaliseren dat de waargenomen fouten de verwachte fouten overschrijden, wat zou betekenen dat er onvoldoende werk is verricht.

- Relatief: voer uw verwachte fouten in als een percentage ten opzichte van de totale grootte van de selectie.
- Absoluut: Voer uw verwachte fouten in als de som van (proportionele) fouten.

#### Kritieke items
- Negatieve boekwaarden: Isoleert negatieve boekwaarden van de populatie.
  - Bewaren: houdt negatieve boekwaarden die moeten worden geïnspecteerd in het monster.
  - Verwijderen: verwijdert negatieve boekwaarden.

#### Tabellen opmaken
- Cijfers: geef tabeluitvoer weer als getallen.
- Percentages: geef de tabeluitvoer weer als percentages.
- Geëxtrapoleerde bedragen: geef tabeluitvoer weer als geldwaarden.

#### Veronderstel homogene taints
Als u op dit vakje klikt, kunt u de bekende en onbekende afwijking in de populatie scheiden om efficiënter te werken. Let op dat hiervoor de aanname vereist is dat de taints in de steekproef representatief zijn voor de taints in het onzichtbare deel van de populatie.

### Uitgang
---

#### Evaluatieoverzicht
- Materialiteit: indien aanwezig, de uitvoeringsmaterialiteit.
- Min. precisie: indien aanwezig, de minimale precisie.
- Steekproefomvang: De steekproefomvang (aantal eenheden).
- Fouten: het aantal foutieve elementen in de selectie.
- Taint: De som van de proportionele fouten. Gecontroleerde items kunnen worden geëvalueerd terwijl de omvang van de afwijking wordt meegenomen door hun taint te berekenen. De taint van een item *i* is het proportionele verschil tussen de boekwaarde van dat item (*y*) en de controlewaarde (true) van het item (*x*). Positieve taint worden geassocieerd met te hoge bedragen, terwijl negatieve taints optreden wanneer items worden onderschat.
<img src="%HELP_FOLDER%/img/taints.png" />
- Meest waarschijnlijke fout: De meest waarschijnlijke fout in de populatie.
- x-% Betrouwbaarheidsgrens: Bovengrens van de afwijking in de populatie.
- Precisie: verschil tussen bovengrens en meest waarschijnlijke fout.
- BF-+: De Bayes-factor voor de test.

#### Prior en posterieur
- Functionele vorm: De functionele vorm van de distributie.
- Ondersteuning H-: Totale kans in het bereik van H- onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Ondersteuning H+: Totale kans in het bereik van H+ onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Verhouding H- / H+: Kansen in het voordeel van H- onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Gemiddelde: gemiddelde van de verdeling.
- Mediaan: Mediaan van de verdeling.
- Mode: Mode van de distributie.
- Bovengrens: x-% percentiel van de verdeling.
- Precisie: verschil tussen de bovengrens en de wijze van verdeling.

#### Correcties voor de bevolking
- Correctie: Het van de populatie af te trekken bedrag of percentage.

#### Aannamecontroles
- n: steekproefomvang.
- Pearsons r: Pearson-correlatiecoëfficiënt.
- x-% bovengrens: bovengrens voor correlatiecoëfficiënt.
- p: p-waarde voor de test.
- BF-0: Bayes-factor voor de test.

#### Figuren
- Prior en posterior: Produceert een plot die de prior verdeling en de posterieure distributie toont na observatie van het beoogde monster.
  - Extra info: Produceert stippen op de materialiteit.
- Posterior predictive: Produceert een plot van de voorspellingen van de posterieure verdeling.
- Steekproefdoelstellingen: Produceert een staafdiagram waarin de materialiteit, de maximale afwijking en de meest waarschijnlijke fout (MLE) worden vergeleken.
- Spreidingsplot: Produceert een spreidingsplot die boekwaarden van de selectie vergelijkt met hun controlewaarden. Waarnemingen die fout zijn, zijn rood gekleurd.
  - Correlatie weergeven: Voegt de correlatie tussen de boekwaarden en de controlewaarden toe aan de plot.
  - Item-ID's weergeven: voegt de item-ID's toe aan de plot.

### Referenties
---
- AICPA (2017). <i>Auditgids: controlesteekproeven</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesiaanse en klassieke auditsteekproeven. R-pakket versie 0.6.2.
- Dyer, D., & Pierce, R.L. (1993). Over de keuze van de voorafgaande verdeling bij hypergeometrische steekproeven. <i>Communicatie in statistiek-theorie en methoden</i>, 22(8), 2125-2146.
- Stewart, TR (2013). Een Bayesiaans audit assurance-model met toepassing op het component materialiteitsprobleem bij groepsaudits (proefschrift).
- de Swart, J., Wille, J., & Majoor, B. (2013). Het 'Push-Left'-Principe als Motor van Data Analytics in de Accountantcontrole [Het 'Push-Left'-Principe als aanjager van Data Analytics in Financial Audit]. <i>Maandblad voor Accountancy en Bedrijfseconomie</i>, 87, 425-432.

### R-pakketten
---
- jfaEvaluatie
===

De evaluatieanalyse stelt de gebruiker in staat om op basis van een controlesteekproef conclusies te trekken over de totale afwijking in de populatie.

<img src="%HELP_FOLDER%/img/workflowEvaluation.png" />

Zie de handleiding van de Audit module (download [hier](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) voor meer gedetailleerde informatie over deze analyse.

### Invoer
---

#### Steekproef Doelstellingen
- Uitvoeringsmaterialiteit: ook wel de bovengrens voor de fout, het maximaal toelaatbare afwijkingspercentage of de maximaal toelaatbare afwijking genoemd, de uitvoeringsmaterialiteit is de bovengrens van de toelaatbare afwijking in de te testen populatie. Door te toetsen aan een uitvoeringsmaterialiteit bent u in staat een steekproef te plannen om bewijs te verzamelen voor of tegen de stelling dat de populatie als geheel geen afwijkingen bevat die als materieel worden beschouwd (d.w.z groter zijn dan de uitvoeringsmaterialiteit). U moet deze doelstelling inschakelen als u met een steekproef uit de populatie wilt weten of de populatie een afwijking boven of onder een bepaalde limiet (de prestatiematerialiteit) bevat. Een lagere uitvoeringsmaterialiteit zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere uitvoeringsmaterialiteit resulteren in een lagere steekproefomvang.
- Minimale precisie: de precisie is het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de fout. Door deze steekproefdoelstelling in te schakelen, kunt u een steekproef plannen zodat het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de afwijking tot een minimumpercentage wordt teruggebracht. U moet deze doelstelling inschakelen als u geïnteresseerd bent in het maken van een schatting van de populatieafwijking met een bepaalde nauwkeurigheid. Een lagere minimaal vereiste precisie zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere minimale precisie resulteren in een vereiste steekproefomvang.

#### Populatie
- Aantal items: Het totale aantal items (rijen) in de populatie.
- Aantal eenheden: Het totale aantal eenheden in de populatie. Let op dat de eenheden items (rijen) of monetaire eenheden (waarden) kunnen zijn, afhankelijk van het controlevraagstuk.

#### Betrouwbaarheid
Het gebruikte betrouwbaarheidsniveau. Het betrouwbaarheidsniveau is het complement van het auditrisico: het risico dat de gebruiker bereid is te nemen om een ​​onjuist oordeel over de populatie te geven. Als u bijvoorbeeld een controlerisico van 5% wilt hebben, staat dit gelijk aan 95% betrouwbaarheid.

#### Opdrachtbox
- Item-ID: een unieke niet-ontbrekende identifier voor elk item in de populatie. Het rijnummer van de items is voldoende.
- Boekwaarden: de variabele die de boekwaarden van de items in de populatie bevat.
- Auditresultaat / waarden: De variabele die de audit (true) waarden bevat, of de binaire classificatie van correct (0) of incorrect (1).
- Selectieteller: De variabele die aangeeft hoe vaak elke waarneming moet worden geëvalueerd.

#### Gegevens
- Raw: gebruik onbewerkte gegevens.
- Overzichtsstatistieken: gebruik overzichtsstatistieken.

#### Auditrisicomodel
- Inherent risico: Een categorie of waarschijnlijkheid voor het inherente risico. Inherent risico wordt gedefinieerd als het risico op een afwijking van materieel belang als gevolg van een fout of weglating in een financieel overzicht als gevolg van een andere factor dan een falen van de interne beheersing.
- Beheersingsrisico: Een categorie of waarschijnlijkheid voor het internecontrolerisico. Het interne beheersingsrisico wordt gedefinieerd als het risico van een afwijking van materieel belang in de financiële overzichten die voortvloeit uit het ontbreken of falen van de relevante interne beheersingsmaatregelen van de gecontroleerde.

Wanneer de auditor informatie heeft die wijst op een laag risicoprofiel van de populatie, kan hij deze informatie gebruiken om zijn vereiste steekproefomvang te verkleinen via het Audit Risk Model (ARM), op voorwaarde dat er geen fouten in de populatie zitten. Volgens de ARM is het auditrisico (AR) een functie van het inherente risico (IR), het internecontrolerisico (CR) en het detectierisico (DR).

*AR = IR x CR x DR*

De auditor beoordeelt het inherente risico en het internecontrolerisico in het algemeen op een 3-puntsschaal om het juiste ontdekkingsrisico te bepalen. Met behulp van de ARM en nul fouten hangt de steekproefomvang af van de risicofactor *R* en de uitvoeringsmaterialiteit. De risicofactor *R* is een functie van het detectierisico (Stewart 2012).

*R = -ln(DR)*

De volgende tabel geeft de waarden van *R* weer als functie van het detectierisico, op voorwaarde dat er geen fouten zijn (Touw en Hoogduin 2012).

| Detectierisico (%) | 1 | 4 | 5 | 10 | 14 |
| :---: | :---: | :---: | :---: | :---: | :---: |
| R | 4.6 | 3.2 | 3 | 2.3 | 2 |

De risicofactor *R* kan worden aangepast met behulp van de beoordelingen van het inherente risico en het internecontrolerisico. Standaard is de standaardmethode voor het instellen van de kansen op IR en CR door de onderstaande tabel te volgen voor een detectierisico van 5%:

| | Hoog | Gemiddeld | Laag |
| :---: | :---: | :---: |
| R | 3 | 2 | 1 |

Deze waarden van *R* worden gebruikt om standaardpercentages in te stellen voor IR en CR. De Audit-module verwerkt de volgende standaardwaarden voor IR en CR:

- Hoog: 100%
- Gemiddeld: 60%
- Laag: 36%

U kunt de waarde van IR en CR handmatig aanpassen door de optie Aangepast te selecteren onder de bijbehorende risicobeoordeling, waardoor de risicofactor *R* wordt aangepast.

#### Methode
- Poisson: gebruikt de Poisson-waarschijnlijkheid om het monster te evalueren.
- Binomiaal: gebruikt de binominale waarschijnlijkheid om het monster te evalueren.
- Hypergeometrische: gebruikt de hypergeometrische waarschijnlijkheid om het monster te evalueren.
- Stringer: de Stringer moet het monster evalueren (Stringer, 1963).
  - LTA-aanpassing: LTA-aanpassing voor de stringer die onvermijdelijk understatements bevat (Leslie, Teitlebaum, & Anderson, 1979).
- Gemiddelde-per-eenheid schatter: Gebruikt de gemiddelde-per-eenheid schatter.
- Directe schatter: Deze methode gebruikt alleen de controlewaarden om de afwijking te schatten (Touw en Hoogduin, 2011).
- Verschilschatter: Deze methode gebruikt het verschil tussen de boekwaarden en de controlewaarden om de afwijking in te schatten (Touw en Hoogduin, 2011).
- Ratio schatter: Deze methode gebruikt de verhouding van correctheid tussen de boekwaarden en de controlewaarden om de afwijking te schatten (Touw en Hoogduin, 2011).
- Regressieschatter: Deze methode gebruikt de lineaire relatie tussen de boekwaarden en de controlewaarden om de afwijking te schatten (Touw en Hoogduin, 2011).

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.

#### Tabellen
- Correcties op populatie: Produceert een tabel die de vereiste correcties op de populatiewaarde bevat om de steekproefdoelstellingen te bereiken.

#### Figuren
- Steekproefdoelstellingen: Produceert een staafdiagram waarin de materialiteit, de maximale afwijking en de meest waarschijnlijke fout (MLE) worden vergeleken.
- Spreidingsplot: Produceert een spreidingsplot die boekwaarden van de selectie vergelijkt met hun controlewaarden. Waarnemingen die fout zijn, zijn rood gekleurd.

#### Kritieke items
- Negatieve boekwaarden: Isoleert negatieve boekwaarden van de populatie.
  - Bewaren: houdt negatieve boekwaarden die moeten worden geïnspecteerd in het monster.
  - Verwijderen: verwijdert negatieve boekwaarden.

#### Tabellen opmaken
- Cijfers: geef tabeluitvoer weer als getallen.
- Percentages: geef de tabeluitvoer weer als percentages.
- Geëxtrapoleerde bedragen: geef tabeluitvoer weer als geldwaarden.

### Uitgang
---

#### Evaluatieoverzicht
- Materialiteit: indien aanwezig, de uitvoeringsmaterialiteit.
- Min. precisie: indien aanwezig, de minimale precisie.
- Steekproefomvang: De steekproefomvang (aantal eenheden).
- Fouten: het aantal foutieve elementen in de selectie.
- Taint: De som van de proportionele fouten. Gecontroleerde items kunnen worden geëvalueerd terwijl de omvang van de afwijking wordt meegenomen door hun taint te berekenen. De taint van een item *i* is het proportionele verschil tussen de boekwaarde van dat item (*y*) en de controlewaarde (true) van het item (*x*). Positieve taint worden geassocieerd met te hoge bedragen, terwijl negatieve taints optreden wanneer items worden onderschat.
<img src="%HELP_FOLDER%/img/taints.png" />
- Meest waarschijnlijke fout: De meest waarschijnlijke fout in de populatie.
- x-% Betrouwbaarheidsgrens: Bovengrens van de afwijking in de populatie.
- Precisie: verschil tussen bovengrens en meest waarschijnlijke fout.
- p: De p-waarde voor de test.

#### Correcties voor de bevolking
- Correctie: Het van de populatie af te trekken bedrag of percentage.

#### Figuren
- Steekproefdoelstellingen: Produceert een staafdiagram waarin de materialiteit, de bovengrens van de afwijking en de meest waarschijnlijke fout (MLE) worden vergeleken.
- Spreidingsplot: Produceert een spreidingsplot die boekwaarden van de selectie vergelijkt met hun controlewaarden. Waarnemingen die fout zijn, zijn rood gekleurd.
  - Correlatie weergeven: Voegt de correlatie tussen de boekwaarden en de controlewaarden toe aan de plot.
  - Item-ID's weergeven: voegt de item-ID's toe aan de plot.

### Referenties
---
- AICPA (2017). <i>Auditgids: controlesteekproeven</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesiaanse en klassieke auditsteekproeven. R-pakket versie 0.6.2.
- Leslie, D. A., Teitlebaum, A.D., Anderson, R.J. (1979). <i>Sampling in dollars: een praktische gids voor auditors</i>. Toronto: Copp Clark Pitman.
- Stringer, K. W. (1963) Praktische aspecten van statistische steekproeven bij auditing. <i>Proceedings of Business and Economic Statistics Section</i>, American Statistical Association.
- Touw, P., & Hoogduin, L. (2011). Statistiek voor audit en controlling.

### R-pakketten
---
- jfaWet van Benford
===

De wet van Benford stelt dat de verdeling van de eerste cijfers in een populatie van nature een bepaalde verdeling volgt. Bij audits kan het beoordelen of een verdeling van cijfers in de populatie in overeenstemming is met de wet van Benford aanvullend bewijs leveren dat de transacties in de populatie mogelijk nader onderzoek behoeven.

*Opmerking:* Niet-conformiteit met de wet van Benford duidt niet noodzakelijk op fraude. Een analyse van de wet van Benford mag daarom alleen worden gebruikt om inzicht te krijgen of een populatie nader onderzoek nodig heeft.

### Invoer
---

#### Opdrachtbox
- Variabele: In dit vak wordt de variabele geselecteerd waarvan de cijfers moeten worden getoetst aan de referentieverdeling. De waarde nul (0) wordt weggelaten uit de gegevens.

#### Verwijzing
- Wet van Benford: toets de cijfers aan de wet van Benford.
- Uniforme verdeling: toets de cijfers aan de uniforme verdeling.

#### Cijfers
- Eerste: controleert alleen het eerste cijfer van de items ten opzichte van de opgegeven verdeling.
- Eerste en tweede: Controleert het eerste en tweede cijfer van de items tegen de opgegeven verdeling.
- Laatste: controleert alleen het laatste cijfer van de items ten opzichte van de opgegeven distributie.

#### Bayes-factor
- BF10 : Bayes-factor om bewijs te kwantificeren voor de alternatieve hypothese ten opzichte van de nulhypothese.
- BF01 : Bayes-factor om bewijs voor de nulhypothese te kwantificeren ten opzichte van de alternatieve hypothese.
- Log (BF10) : Natuurlijke logaritme van BF10.

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.
  - Betrouwbaarheid: Het betrouwbaarheidsniveau dat in de verklarende tekst wordt gebruikt.

### Uitgang
---

#### Goodness-of-fit tafel
- n: Het totaal aantal waarnemingen in de dataset.
- X<sup>2</sup>: de waarde van de chi-kwadraat-teststatistiek.
- df: vrijheidsgraden geassocieerd met de Chi-kwadraattoets.
- p: de *p*-waarde die is gekoppeld aan de Chi-kwadraattoets.
- BF: De Bayes-factor als gevolg van een niet-informatieve prior.

#### Frequentietabel
- Voorloop / Laatste cijfer: Het cijfer waarvoor de informatie in de rij geldt.
- Telling: De waargenomen tellingen van de cijfers.
- Relatieve frequentie: De waargenomen relatieve frequentie van de cijfers.
- Wet van Benford / Uniforme verdeling: De verwachte relatieve frequentie van de cijfers.

#### Figuren
- Waargenomen vs. verwacht: Produceert een grafiek die de waargenomen verdeling van cijfers in de populatie laat zien in vergelijking met de verwachte verdeling volgens de wet van Benford of de uniforme verdeling.

### Referenties
---
- Derks, K (2021). digitTests: Tests voor het detecteren van onregelmatige cijferpatronen. R-pakket versie 0.1.0.

### R-pakketten
---
- cijfertoetsenSampling Workflow
===

The task of an auditor is to make a judgment regarding the fairness of the presented transactions in a population. When the auditor has access to the raw population data, they can use the *audit workflow* to calculate how many samples need to be evaluated in order to meet a certain confidence in their judgment. The user can then sample these items from the population, inspect and audit these items, and perform statistical inference about the misstatement in the population. The sampling workflow guides the auditor through the audit process, making the correct choices of calculations along the way.

Please see the manual of the Audit module (download [here](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) for more detailed information about this analysis.

### Workflow
---

- Planning: Calculate the minimum sample size to achieve your sampling objectives with the specified confidence.
- Selection: Select the required sampling units from the population.
- Execution: Annotate the selection with your assessment of the fairness of the selected items.
- Evaluation: Make a population statement based on your annotated selection.

<img src="%HELP_FOLDER%/img/workflow.png" />

### Input - Planning
---

#### Sampling Objectives
- Performance materiality: Also called the upper error limit, the tolerable deviation rate, or the tolerable misstatement, the performance materiality is the upper bound of tolerable misstatement in the population to be tested. By testing against a performance materiality, you are able to plan a sample in order to collect evidence for or against the statement that the population as a whole does not contain misstatements that are considered material (i.e., are greater than the upper bound of tolerable misstatement). You should enable this objective when you want to find out whether the population contains misstatement above or below a certain limit (the performance materiality) using a sample of the population. A lower performance materiality will result in a higher required sample size. Vice versa, a higher performance materiality will result in a lower required sample size.
- Minimum precision: The precision is the the difference between the estimated most likely error and the upper bound on the misstatement. By enabling this sampling objective, you are be able to plan a sample so that the difference between the estimated most likely error and the upper bound on the misstatement is reduced to a minimum percentage. You should enable this objective if you are interested in making an estimate of the population misstatement with a certain accuracy. A lower minimum required precision will result in a higher required sample size. Vice versa, a higher minimum required precision will result in a lower required sample size.

#### Confidence
The confidence level used. The confidence level is the complement of the audit risk: the risk that the user is willing to take to give an incorrect judgment about the population. For example, if you want to have an audit risk of 5%, this equals 95% confidence.

#### Assignment Box
- Item ID: A unique non-missing identifier for every item in the population. The row number of the items is sufficient.
- Book Values: The variable that contains the book values of the items in the population.

#### Audit Risk Model
- Inherent risk: A category or probability for the inherent risk. Inherent risk is defined as the risk of material misstatement posed by an error or omission in a financial statement due to a factor other than a failure of internal control.
- Control risk: A category or probability for the internal control risk. Control risk is defined as the risk of a material misstatement in the financial statements arising due to absence or failure in the operation of relevant controls of the auditee.

When the auditor has information that indicates a low-risk profile on the population, they can use this information to reduce their required sample size via the Audit Risk Model (ARM) provided that there are no errors in the population. According to the ARM, the audit risk (AR) is a function of the inherent risk (IR), the internal control risk (CR), and the detection risk (DR).

*AR = IR x CR x DR*

The auditor assesses inherent risk and internal control risk generally on a 3-point scale to determine the appropriate detection risk. Using the ARM and zero errors the sample size depends on the risk factor *R* and the performance materiality. The risk factor *R* is a function of the detection risk (Stewart 2012).

*R = -ln(DR)*

The following table presents values of *R* as a function of the detection risk, provided that there are zero errors (Touw and Hoogduin 2012).

| Detection risk (%) | 1 | 4 | 5 | 10 | 14 |
| :---: | :---: | :---: | :---: | :---: | :---: |
| R | 4.6 | 3.2 | 3 | 2.3 | 2 |

The risk factor *R* can be adjusted using the assessments of the inherent risk and the internal control risk. By default, the standard method of setting the probabilities of IR and CR is by following the table below for a detection risk of 5%:

|  | High | Medium | Low | 
| :---: | :---: | :---: |
| R | 3 | 2 | 1 |

These values of *R* are used to set default percentages for IR and CR. The Audit module handles the following default values for IR and CR:

- High: 100%
- Medium: 60%
- Low: 36%

You can manually adjust the value of IR and CR by selecting the Custom option under the corresponding risk assessment, thus adjusting the risk factor *R*.

#### Expected errors in Sample
The expected errors are the tolerable errors that can be found in the sample while still achieving the specified sampling objectives. A sample size is calculated so that, when the number of expected errors is found in the sample, the desired confidence is retained.

*Note:* It is advised to set this value conservatively to minimize the probability of the observed errors exceeding the expected errors, which would imply that insufficient work has been done.

- Relative: Enter your expected errors as a percentage relative to the total size of the selection.
- Absolute: Enter your expected errors as the sum of (proportional) errors.

#### Probability distribution
- Hypergeometric: The hypergeometric distribution assumes a finite population size and is therefore generally used when the population size is small. It is a probability distribution that models the number of errors (*K*) in the population as a function of the population size (*N*), the number of observed found errors (*k*) and the number of correct transactions (*n*).
- Binomial: The binomial distribution assumes an infinite population size and is therefore generally used when the population size is large. It is a probability distribution that models the rate of misstatement (*\u03B8*) as a function of the observed number of errors (*k*) and the number of correct transactions (*n - k*). Because the binomial distribution strictly does not accommodate partial errors, it is generally used when you are not planning a monetary unit sample.
- Poisson: The Poisson distribution assumes an infinite population size and is therefore generally used when the population size is large. It is a probability distribution that models the rate of misstatement (*\u03B8*) as a function of the observed sample size (*n*) and the sum of the proportional errors (*t*). Because the Poisson distribution accommodates partial errors it is generally used when you are planning a monetary unit sample.

#### Display
- Explanatory Text: When checked, enables explanatory text in the analysis to help interpret the procedure and the statistical results.

#### Tables
- Descriptive statistics: Produces a table with descriptive statistics of the book values in the population.

#### Plots
- Distribution of book values: Produces a histogram of the book values in the population.
- Compare sample sizes: Produces a plot that compares the sample size 1) across probability distributions, and 2) across the number of expected errors in the sample.
- Assumed error distribution: Produces a plot that displays the probability distribution implied by the input options and the calculated sample size.

#### Critical Items
- Negative book values: Isolates negative book values from the population.
  - Keep: Keeps negative book values to be inspected in the sample.
  - Remove: Removes negative book values.

#### Format Tables
- Numbers: Display table output as numbers.
- Percentages: Display table output as percentages.

#### Iterations
- Increment: The increment alows you to limit the possible sample sizes to a multiple of its value. For example, an increment of 5 allows only sample sizes of 5, 10, 15, 20, 25, etc.
- Maximum: The maximum allows you to limit the sample size with a maximum.

### Ouput - Planning
---

#### Planning Summary
- Performance materiality: When provided, the performance materiality.
- Min. precision: When provided, the minimum precision.
- Expected errors: The number (sum of proportional taints) of expected / tolerable errors in the sample.
- Minimum sample size: The minimum sample size.

#### Descriptive Statistics
- Population size: Number of items in the population.
- Value: Total value of the book values.
- Absolute value: Absolute value of the book values.
- Mean: Mean of the book values.
- Std. deviation: Standard deviation of the book values.
- Quartile: Quartiles of the book values.

#### Plots
- Distribution of book values: Produces a histogram of the book values in the population.
- Compare sample sizes: Produces a plot that compares the sample size 1) across probability distributions, and 2) across the number of expected errors in the sample.
- Assumed error distribution: Produces a plot that displays the probability distribution implied by the input options and the calculated sample size.

### Input - Selection
---

#### Assignment Box
- Ranking Variable: When provided, the population is first ranked in ascending order with respect to the values of this variable.
- Additional Variables: Any other variables that should be included in the sample.

#### Sample Size
The required number of sampling units that should be selected from the population. Be aware that the sampling units are determined by the *units* option. By default, when no book values are provided, the sampling units are items (rows). When book values are provided, the ideal sampling units to use are monetary units.

#### Sampling Units
- Items: Performs selection using the items in the population as sampling units.
- Monetary units: Performs selection using the monetary units in the population as sampling units. This method is preferred when you want to include more items with a high value in the sample.

#### Method
- Fixed interval sampling: Performs selection by dividing the population in equal intervals and selecting a fixed unit in each interval. Any item with a value larger than the interval will always be included in the sample.
  - Starting point: Selects which sampling unit is selected from each interval.
- Cell sampling: Performs selection by dividing the population in equal intervals and selecting a variable unit in each interval. Any item with a value larger than twice the interval will always be included in the sample.
  - Seed: Selects the seed for the random number generator in order to reproduce results.
- Random sampling: Performs random selection in which each sampling unit has an equal chance of being selected.
  - Seed: Selects the seed for the random number generator in order to reproduce results.

#### Randomize Item Order
Randomizes the items in the population before selection is performed.

#### Tables
- Descriptive statistics: Produces a table containing descriptive information about numerical variables in the selection. Statistics that are included are the mean, the median, the standard deviation, the variance, the minimum, the maximum, and the range.
- Selected items: Produces a table containing the selected transactions along with any additional observations provided in the additional variables field.

### Output - Selection
---

#### Selection Summary
- No. units: The number of selected sampling units from the population.
- No. items: The number of selected items from the population.
- Selection value: The total value of the selected items. Only displayed when monetary unit sampling is used.
- % of population size / value: The selected proportion of the total size or value of the population.

#### Information about Monetary Interval Selection
- Items: The number of items in the population.
- Value: The value of the items in the population.
- Selected items: The number of items in the sample.
- Selected units: The number of selected units from the population.
- Selection value: The value of the items in the sample.
- % of total value: The selected proportion of the total value of the items compared to the items in the population.

#### Descriptive Statistics
- Valid: Number of valid cases.
- Mean: Arithmetic mean of the data points.
- Median: Median of the data points.
- Std. deviation: Standard deviation of the data points.
- Variance: Variance of the data points.
- Range: Range of the data points.
- Minimum: Minimum of the data points.
- Maximum: Maximum of the data points.

#### Selected Items
- Row: The row number of the item.
- Selected: The number of times (a unit in) the item is selected.

### Input - Execution
---

#### Annotation
- Audit value: Annotate the items in the selection with their audit (true) values. This approach is recommended (and automatically selected) when the items have a monetary value.
- Correct / Incorrect: Annotate the items in the selection with correct (0) or incorrect (1). This approach is recommended (and automatically selected) when your items do not have a monetary value.

### Input - Evaluation
---

#### Assignment Box
- Audit result / values: The variable that contains the audit (true) values, or the binary classification of correct (0) or incorrect (1).

#### Method
- Poisson: Uses the Poisson likelhood to evaluate the sample.
- Binomial: Uses the binomial likelhood to evaluate the sample.
- Hypergeometric: Uses the hypergeometric likelhood to evaluate the sample.
- Stringer: The Stringer bound to evaluate the sample (Stringer, 1963).
  - LTA adjustment: LTA adjustment for the stringer bound to incorporate understatements (Leslie, Teitlebaum, & Anderson, 1979).
- Mean-per-unit estimator: Uses the mean-per-unit estimator.
- Direct estimator: This method uses only the audit values to estimate the misstatement (Touw and Hoogduin, 2011).
- Difference estimator: This method uses the difference between the book values and the audit values to estimate the misstatement (Touw and Hoogduin, 2011).
- Ratio estimator: This method uses the ratio of correctness between the book values and the audit values to estimate the misstatement (Touw and Hoogduin, 2011).
- Regression estimator: This method uses the linear relation between the book values and the audit values to estimate the misstatement (Touw and Hoogduin, 2011).

#### Tables
- Corrections to population: Produces a table that contains the required corrections to the population value to achieve the sampling objectives.

#### Plots
- Sampling objectives: Produces a bar chart comparing the materiality, maximum misstatement and most likely error (MLE).
- Scatter plot: Produces a scatter plot comparing book values of the selection against their audit values. Observations that are in error are colored in red.

### Output - Evaluation
---

#### Evaluation summary
- Materiality: When provided, the performance materiality.
- Min. precision: When provided, the minimum precision.
- Sample size: The sample size (number of units).
- Errors: The number of erroneous elements in the selection.
- Taint: The sum of the proportional errors. Audited items can be evaluated while incorporating the magnitude of the misstatement by calculating their taints. The taint of an item *i* is the proportional difference between that item's book value (*y*) and the item's audit (true) value (*x*). Positive taints are associated with overstatements, while negative taints occur when items are understated.
<img src="%HELP_FOLDER%/img/taints.png" />
- Most likely error: The most likely error in the population.
- x-% Confidence bound: Upper bound on the misstatement in the population.
- Precision: Difference between upper bound and most likely error.
- p: The p-value for the test.

#### Corrections to Population
- Correction: The amount or percentage to be deducted from the population.

#### Plots
- Sampling objectives: Produces a bar chart comparing the materiality, upper bound on the misstatement and most likely error (MLE).
- Scatter plot: Produces a scatter plot comparing book values of the selection against their audit values. Observations that are in error are colored in red.
  - Display correlation: Adds the correlation between the book values and the audit values to the plot.
  - Display item ID's: Adds the item ID's to the plot.

### References
---
- AICPA (2017). <i>Audit Guide: Audit Sampling</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesian and Classical Audit Sampling. R package version 0.6.2.
- Leslie, D. A., Teitlebaum, A. D., Anderson, R. J. (1979). <i>Dollar-unit Sampling: A Practical Guide for Auditors</i>. Toronto: Copp Clark Pitman.
- Stringer, K. W. (1963) Practical aspects of statistical sampling in auditing. <i>Proceedings of Business and Economic Statistics Section</i>, American Statistical Association.
- Touw, P., & Hoogduin, L. (2011). Statistiek voor audit en controlling.

### R Packages
---
- jfa
Bayesiaanse Planning
===

Met de Bayesiaanse planningsanalyse kan de gebruiker een minimale steekproefomvang berekenen op basis van een reeks steekproefdoelstellingen en samenvattende statistieken van de populatie. Houd er rekening mee dat wanneer u toegang heeft tot de onbewerkte populatiegegevens, u misschien de auditworkflow wilt gebruiken, een analyse die u door het steekproefproces leidt.

<img src="%HELP_FOLDER%/img/workflowPlanning.png" />

Zie de handleiding van de Audit module (download [hier](https://github.com/jasp-stats/jaspAudit/raw/master/man/manual.pdf)) voor meer gedetailleerde informatie over deze analyse.

### Invoer
---

#### Steekproef Doelstellingen
- Uitvoeringsmaterialiteit: ook wel de bovengrens voor de fout, het maximaal toelaatbare afwijkingspercentage of de maximaal toelaatbare afwijking genoemd, de uitvoeringsmaterialiteit is de bovengrens van de toelaatbare afwijking in de te testen populatie. Door te toetsen aan een uitvoeringsmaterialiteit bent u in staat een steekproef te plannen om bewijs te verzamelen voor of tegen de stelling dat de populatie als geheel geen afwijkingen bevat die als materieel worden beschouwd (d.w.z groter zijn dan de uitvoeringsmaterialiteit). U moet deze doelstelling inschakelen als u met een steekproef uit de populatie wilt weten of de populatie een afwijking boven of onder een bepaalde limiet (de prestatiematerialiteit) bevat. Een lagere uitvoeringsmaterialiteit zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere uitvoeringsmaterialiteit resulteren in een lagere steekproefomvang.
- Minimale precisie: de precisie is het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de fout. Door deze steekproefdoelstelling in te schakelen, kunt u een steekproef plannen zodat het verschil tussen de geschatte meest waarschijnlijke fout en de bovengrens van de afwijking tot een minimumpercentage wordt teruggebracht. U moet deze doelstelling inschakelen als u geïnteresseerd bent in het maken van een schatting van de populatieafwijking met een bepaalde nauwkeurigheid. Een lagere minimaal vereiste precisie zal resulteren in een hogere steekproefomvang. Omgekeerd zal een hogere minimale precisie resulteren in een vereiste steekproefomvang.

#### Populatie
- Aantal eenheden: Het totale aantal eenheden in de populatie. Let op dat de eenheden items (rijen) of monetaire eenheden (waarden) kunnen zijn, afhankelijk van het controlevraagstuk.

#### Betrouwbaarheid
Het gebruikte betrouwbaarheidsniveau. Het betrouwbaarheidsniveau is het complement van het auditrisico: het risico dat de gebruiker bereid is te nemen om een ​​onjuist oordeel over de populatie te geven. Als u bijvoorbeeld een controlerisico van 5% wilt hebben, staat dit gelijk aan 95% betrouwbaarheid.

#### Verwachte fouten in steekproef
De verwachte fouten zijn de toelaatbare fouten die in de steekproef kunnen worden gevonden terwijl de gespecificeerde steekproefdoelstellingen nog steeds worden bereikt. Er wordt een steekproefomvang berekend zodat, wanneer het aantal verwachte fouten in de steekproef wordt gevonden, het gewenste vertrouwen behouden blijft.

*Opmerking:* Het wordt aangeraden om deze waarde conservatief in te stellen om de kans te minimaliseren dat de waargenomen fouten de verwachte fouten overschrijden, wat zou betekenen dat er onvoldoende werk is verricht.

- Relatief: voer uw verwachte fouten in als een percentage ten opzichte van de totale grootte van de selectie.
- Absoluut: Voer uw verwachte fouten in als de som van (proportionele) fouten.

#### Kansverdeling
- Gamma: De gammaverdeling gaat samen met de Poisson-verdeling. De Poisson-verdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een verdeling waarmee het percentage afwijkingen (*\u03B8*) modelleerd als functie van de waargenomen steekproefomvang (*n*) en de som van de gevonden proportionele fouten (*t*). Omdat de gammaverdeling deelfouten mogelijk maakt, wordt deze over het algemeen gebruikt bij het plannen van een steekproef in munteenheden (Stewart, 2013).
- Beta: de bètaverdeling gaat samen met de binominaalverdeling. De binominaalverdeling gaat uit van een oneindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang groot is. Het is een verdeling waarmee het percentage afwijkingen (*\u03B8*) wordt gemodelleerd als functie van het waargenomen aantal fouten (*k*) en het aantal correcte transacties (*n - k*). Omdat de binominaalverdeling strikt genomen geen rekening houdt met gedeeltelijke fouten, wordt deze over het algemeen gebruikt wanneer u geen steekproef in een munteenheid plant. De betaverdeling is echter geschikt voor gedeeltelijke fouten en kan ook worden gebruikt voor steekproeven op monetaire eenheden (de Swart, Wille & Majoor, 2013).
- Beta-binomiaal: De beta-binomiaalverdeling gaat samen met de hypergeometrische verdeling (Dyer & Pierce, 1993). De hypergeometrische verdeling gaat uit van een eindige populatieomvang en wordt daarom over het algemeen gebruikt wanneer de populatieomvang klein is. Het is een verdeling waarmee het aantal fouten (*K*) in de populatie wordt gemodelleerd als functie van de populatieomvang (*N*), het aantal geobserveerde gevonden fouten (*k*) en het aantal correcte transacties (*N*).

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.

#### Tabellen
- Impliciete eerdere steekproef: Produceert een tabel die de impliciete steekproef weergeeft waarop de prior verdeling is gebaseerd.
- Prior en posterior: Produceert een tabel waarin de eerdere en verwachte posterieure distributie worden samengevat door middel van verschillende statistieken, zoals hun functionele vorm, hun eerdere en verwachte posterieure kansen en kansen, en de verschuiving daartussen.

#### Figuren
- Prior en posterior: Produceert een plot die de prior verdeling en de posterieure distributie toont na observatie van het beoogde monster.
  - Extra info: Produceert stippen op de materialiteit.
- Prior predictive: Produceert een plot van de voorspellingen van de prior verdeling.
- Vergelijk steekproefomvang: Produceert een plot die de steekproefomvang vergelijkt 1) over kansverdelingen, en 2) over het aantal verwachte fouten in de steekproef.

#### Prior
- Standaard: deze optie neemt geen informatie op in de statistische analyse en gaat daarom uit van een verwaarloosbare en conservatieve eerdere verdeling.
- Handmatig: geef de parameters van de prior verdeling op.
- Eerdere steekproef: Maak een eerdere verdeling op basis van een eerdere steekproef.
  - Grootte: eerdere steekproefomvang.
  - Fouten: eerder gevonden fouten.
- Onpartijdig: maak een voorafgaande verdeling die onpartijdig is met betrekking tot de geteste hypothesen.
- Risicobeoordelingen: Vertaal informatie uit het auditrisicomodel naar een eerdere verspreiding.
  - Inherent risico: Een categorie of waarschijnlijkheid voor het inherente risico. Inherent risico wordt gedefinieerd als het risico op een afwijking van materieel belang als gevolg van een fout of weglating in een financieel overzicht als gevolg van een andere factor dan een falen van de interne beheersing.
  - Beheersingsrisico: Een categorie of waarschijnlijkheid voor het internecontrolerisico. Het interne beheersingsrisico wordt gedefinieerd als het risico van een afwijking van materieel belang in de financiële overzichten die voortvloeit uit het ontbreken of falen van de relevante interne beheersingsmaatregelen van de gecontroleerde.

#### Tabellen opmaken
- Cijfers: geef tabeluitvoer weer als getallen.
- Percentages: geef de tabeluitvoer weer als percentages.

#### Stapgrootte
Met de stapgrootte kunt u de mogelijke steekproefomvang beperken tot een veelvoud van de waarde. Een stapgrootte van 5 staat bijvoorbeeld alleen steekproefomvang van 5, 10, 15, 20, 25, enz. toe.

### Uitgang
---

#### Planningsoverzicht
- Uitvoeringsmaterialiteit: indien aanwezig, de uitvoeringsmaterialiteit.
- Min. precisie: indien aanwezig, de minimale precisie.
- Verwachte fouten: Het aantal (som van proportionele taints) verwachte / toelaatbare fouten in de steekproef.
- Minimale steekproefomvang: de minimale steekproefomvang.

#### Impliciete voorafgaande steekproef
- Impliciete steekproefomvang: de steekproefomvang die overeenkomt met de voorafgaande informatie.
- Impliciete fouten: het aantal fouten dat gelijk is aan de eerdere informatie.

#### Prior en posterieur
- Functionele vorm: De functionele vorm van de distributie.
- Ondersteuning H-: Totale kans in het bereik van H- onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Ondersteuning H+: Totale kans in het bereik van H+ onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Verhouding H- / H+: Kansen in het voordeel van H- onder de verdeling. Wordt alleen weergegeven bij toetsing aan een uitvoeringsmaterialiteit.
- Gemiddelde: gemiddelde van de verdeling.
- Mediaan: Mediaan van de verdeling.
- Mode: Mode van de distributie.
- Bovengrens: x-% percentiel van de verdeling.
- Precisie: verschil tussen de bovengrens en de wijze van verdeling.

#### Figuren
- Prior en posterior: Produceert een plot die de prior verdeling en de posterieure distributie toont na observatie van het beoogde monster.
  - Extra info: Produceert stippen op de materialiteit.
- Prior predictive: Produceert een plot van de voorspellingen van de prior verdeling.
- Vergelijk steekproefomvang: Produceert een plot die de steekproefomvang vergelijkt 1) over kansverdelingen, en 2) over het aantal verwachte fouten in de steekproef.

### Referenties
---
- AICPA (2017). <i>Auditgids: controlesteekproeven</i>. American Institute of Certified Public Accountants.
- Derks, K. (2022). jfa: Bayesiaanse en klassieke auditsteekproeven. R-pakket versie 0.6.2.
- Dyer, D., & Pierce, R.L. (1993). Over de keuze van de voorafgaande verdeling bij hypergeometrische steekproeven. <i>Communicatie in statistiek-theorie en methoden</i>, 22(8), 2125-2146.
- Stewart, TR (2013). Een Bayesiaans audit assurance-model met toepassing op het component materialiteitsprobleem bij groepsaudits (proefschrift).
- de Swart, J., Wille, J., & Majoor, B. (2013). Het 'Push-Left'-Principe als Motor van Data Analytics in de Accountantcontrole [Het 'Push-Left'-Principe als aanjager van Data Analytics in Financial Audit]. <i>Maandblad voor Accountancy en Bedrijfseconomie</i>, 87, 425-432.

### R-pakketten
---
- jfaHerhaalde Waardeanalyse
===

Deze analyse analyseert de frequentie waarmee waarden binnen een dataset worden herhaald ('nummerbundeling' genoemd) om statistisch vast te stellen of er waarschijnlijk met de gegevens is geknoeid. In tegenstelling tot de wet van Benford onderzoekt deze benadering het hele getal in één keer, niet alleen het eerste of laatste cijfer (Simonsohn, 2019).

Om te bepalen of de gegevens een overmatige hoeveelheid opeenhoping vertonen, wordt de nulhypothese getest dat de gegevens geen onverwachte hoeveelheid herhaalde waarden bevatten. Om te kwantificeren wat er wordt verwacht, vereist deze test de aanname dat de gehele delen van de getallen niet zijn gekoppeld aan hun decimale delen.

### Invoer
---

#### Opdrachtbox
- Variabele: in dit vak wordt de variabele geselecteerd waarvan de cijfers moeten worden geanalyseerd op herhaalde waarden.

#### Testen
- Gemiddelde frequentie: bereken de gemiddelde frequentie van de gegevens.
- Entropie: bereken de entropie van de gegevens.

#### Shuffle decimale cijfers
- Laatste: het laatste decimaalcijfer wordt geschud.
- Laatste twee: de laatste twee decimale cijfers worden geschud.
- Alles: alle decimale cijfers worden geschud.

#### Weergave
- Verklarende tekst: indien ingeschakeld, wordt verklarende tekst in de analyse ingeschakeld om de procedure en de statistische resultaten te helpen interpreteren.
  - Betrouwbaarheid: Het betrouwbaarheidsniveau dat in de verklarende tekst wordt gebruikt.

#### Tabellen
- Aannamecontroles: deze tabel toont de correlatie tussen de gehele delen van de getallen en hun decimale tegenhangers. Om aan de vereiste aannames voor deze procedure te voldoen, mag deze correlatie niet bestaan. Deze tabel geeft ook de correlatie weer tussen de steekproeven van de twee simulatieruns (gemiddelde frequentie en entropie).
- Frequentietabel: Produceert een tabel met het aantal en het percentage voor elke unieke waarde in de dataset.

#### Figuren
- Waargenomen vs. verwacht: Produceert een histogram van de verwachte gemiddelde frequenties en/of entropie vs. de waargenomen gemiddelde frequentie en/of entropie.
- Histogram: Produceert een histogram met een enkele bak voor elke waargenomen waarde.

#### Geavanceerde mogelijkheden
- Aantal monsters: Het aantal monsters dat moet worden gebruikt voor het simuleren van de p-waarde.
- Seed: Selecteert de seed voor de generator van willekeurige getallen om resultaten te reproduceren.

### Uitgang
---

#### Herhaalde Waarden Test
- n: Het aantal waarnemingen in de gegevens.
- Frequentie: De gemiddelde frequentie waarmee getallen in de gegevens worden herhaald. De formule voor de gemiddelde frequentie is *AF = &#8721; f&#7522;&#178; / &#8721; f&#7522;* waar f&#7522; is de frequentie van elke unieke waarde *i* in de dataset.
- Entropie: De entropie is het gemiddelde informatieniveau dat inherent is aan de uitkomsten van de variabele. De entropie wordt berekend als *S = - &#8721; (p&#7522; &#215; log(p&#7522;))* waarbij p&#7522; is het aandeel waarnemingen met elke waarde (dus *p&#7522; = f&#7522; / N*).

#### Aannamecontroles
- n: steekproefomvang.
- r: Pearson-correlatiecoëfficiënt.
- t: t-waarde.
- df: vrijheidsgraden.
- p: p-waarde.

#### Frequentietabel
- Waarde: de waarde in de rij.
- Telling: het aantal keren dat de waarde wordt waargenomen.
- Percentage: Het percentage keren dat de waarde wordt waargenomen.

#### Figuren
- Waargenomen vs. verwacht: geeft de waargenomen vs. de gesimuleerde waarde(n) weer.
- Histogram: geeft een histogram weer met een enkele bak voor elke waargenomen waarde.

### Referenties
---
- Derks, K (2021). digitTests: Tests voor het detecteren van onregelmatige cijferpatronen. R-pakket versie 0.1.0.
- Simohnsohn, Verenigde Staten (2019, 25 mei). Nummerbundeling: een nieuw hulpmiddel voor forensische gegevensanalyse. Opgehaald van [http://datacolada.org/77](http://datacolada.org/77).

### R-pakketten
---
- cijfertoetsenSampling units
==========================

Sampling from the population requires knowledge of the sampling units; physical representations of the population that need to be audited. Sampling units can be individual <i>transactions</i> or <i>monetary units</i>. For statistical sampling, sampling probabilities are assignned to the sampling units of the population elements. The total collection of all sampling units with an assigned selection probability is named the sampling frame. 

-------

Monetary unit sampling (MUS)
==========================

### Procedure

The sampling frame for monetary unit sampling (MUS) is the total collection of individual monetary units in the population. A sampling unit for MUS is an individual monetary unit within all possible monetary units, like 1$ in a balance of $1000. In the dollar case, the population entry that includes the sampled dollar is selected. For example, if the 10th dollar for a receipt of coffee milk is selected the coffee receipt is audited and its audit (true) value is determined for evaluation. 

### Implications for selection

The implication of a MUS sampling scheme is that transactions are selected with a probability proportional to their value, since larger transactions contain more individual sampling units. It is therefore the preferred selection type for financial populations and is followed up with a Stringer bound as evaluation mechanism in JfA.

### Pitfalls

Be aware that transactions with a book value of $0 also have 0 sampling units. It is therefore recommended that when the population contains transactions that are booked at $0, record sampling is used to sample from the population.Het gebruik van toelichtende tekst
==========================

De toelichtende tekst biedt de gebruiker een uitgebreide uitleg van alle stappen in het inspectieproces. De tekst helpt om de praktische interpretaties van keuzes te begrijpen, en met het verwoorden van statistische bevindingen (in het Nederlands) in communicatie met anderen. Wanneer verklarende tekst aan staat, bevat elk figuur een onderschrift waarin de resultaten van de bijbehorende grafiek worden geïnterpreteerd. 
Fixed interval sampling
==========================

In fixed interval sampling, the population is divided into a set of intervals, the size <i>I</i> of which is computed by dividing the size of the sampling frame by the sample size. A starting point is randomly selected in the first interval and one sampling unit is selected throughout the population at each of the uniform intervals from the starting point. This causes the space between the intervals <i>i</i> to stay the same. In a monetary unit sampling context, the fixed interval method has the property that all transactions that are larger than the interval will always be included in the final sample.

<img src="%HELP_FOLDER%/img/fixedIntervalSampling.png" />Steekproefeenheden
==========================

Enige kennis van steekproefeenheden is vereist om steekproeven uit de populatie te nemen; steekproefeenheden zijn namelijk de fysieke representaties van de populatie die worden geïnspecteerd. Steekproefeenheden kunnen individuele <i>transacties</i> of <i>monetaire eenheden</i> zijn. Bij statistische steekproeven worden selectiekansen toegekend aan de steekproefeenheden van de populatie-elementen. De totale verzameling van alle steekproefeenheden met de toegekende selectiekans heet het steekproefkader. 

-------

Monetaire eenheid sampling (MUS)
==========================

### Procedure

Het steekproefkader voor monetaire eenheid sampling (MUS) is de totale verzameling van alle individuele monetaire eenheden in de populatie. Een steekproefeenheid voor MUS is een individuele monetaire eenheid uit alle mogelijke monetaire eenheden, zoals $1 uit een balans van $1000. In het geval van de dollar, is het populatie aandeel waarin de dollar uit de steekproef zit geselecteerd. Bijvoorbeeld, als de tiende dollar op een bonnetje van koffiemelk wordt geselecteerd, wordt dit koffiemelk bonnetje geïnspecteerd en wordt de geïnspecteerde (ware) waarde bepaald voor evaluatie. 

### Implicaties voor selectie

De bedoeling van MUS is dat transacties worden geselecteerd met een kans die proportioneel staat aan de waarde. Dit aangezien grotere transacties meer individuele steekproefeenheden bevatten. Daarom heeft deze manier van selecteren de voorkeur voor financiële populaties. Hierop volgt een Stringer bound als evaluatiemechanisme in JfA.

### Valkuilen

Let op dat transacties met een inventarisatiewaarde van $0 ook 0 steekproefeenheden bevatten. Daarom wordt aangeraden dat wanneer de populatie transacties bevat van $0, er gebruik wordt gemaakt van record sampling om uit de populatie te selecteren. 
Cell sampling
==========================

Bij cell sampling wordt de populatie in een set van intervallen verdeeld. De grootte <i>I</i> hiervoor wordt berekend door de grootte van het steekproefkader te delen door de steekproefgrootte. In elk interval wordt een steekproefeenheid geselecteerd door aselect een nummer uit ieder interval te nemen. Hierdoor verschilt de ruimte tussen de intervallen <i>i</i>. In het geval van een monetaire steekproefeenheid, worden in de cell sampling methode alle transacties groter dan twee keer het interval meegenomen in de uiteindelijke steekproef.

<img src="%HELP_FOLDER%/img/cellSampling.png" />Vast interval sampling
==========================

Bij vast interval sampling, wordt de populatie verdeeld in een set van intervallen. De grootte <i>I</i> wordt berekend door de grootte van het steekproefkader de delen door de steekproefgrootte. In het eerste interval wordt een willekeurig startpunt geselecteerd, en vervolgens wordt telkens een steekproefeenheid gekozen na een gelijk interval vanaf het startpunt. Hierdoor is de afstand tussen intervallen <i>i</i> altijd hetzelfde. In het geval van een monetaire steekproefeenheid, worden in de vaste interval methode alle transacties groter dan het interval meegenomen in de uiteindelijke steekproef.

<img src="%HELP_FOLDER%/img/fixedIntervalSampling.png" />Incorporating prior information into the statistical anlaysis
==========================

You can incorporate existing information from various sources into the statistical sampling procedure. Possible options are:

- **None:** This option does not incorporate any information into the statistical analysis and therefore assumes a negligible and conservative prior distribution. 

- **Audit Risk Model:** This option can only be selected when testing against a performance materiality. When selected, it incorporates the information from the assessed risks of material misstatement (inherent risk and control risk) from the audit risk model into the statistical analysis. They are mapped to probabilities according to standards, but can also be mapped according to custom preferences.

- High: 100%
- Medium: 60%
- Low: 36%
- Custom

- **Equal prior probabilities:** This option can only be selected when testing against a performance materiality. When selected, it incorporates the information that tolerable misstatement is equally likely to occur a priori as intolerable misstatement. The prior distribution resulting from this method has 50\% of its probability mass below the performance materiality and 50\% of its probability mass above the performance materiality.

- **Custom prior probabilities:** This option can only be selected when testing against a performance materiality. When selected, you can use your own assessment of the a priori probability that the population contains tolerable misstatement. You may also deduct this value as one minus the probability of intolerable misstatement. 

- **Earlier sample:** This option incorporates the information from an earlier sample in the statistical analysis. 

- **Weighted earlier sample:** This option incorporates the information from an earlier sample in the statistical analysis, weighted by a factor *f*, which determines the weight that you assign to the results of the earlier sample.Efficiency technique: Separate known and unknown misstatement
==========================

In estimating the total misstatement in the population the auditor is, in theory, dealing with a part of the population misstatement that is known and a part that is unknown. The known misstatement is divided between two possible sources: the misstatement in the critical transactions and misstatement in the sample transactions. Using the *Separate known and unknown misstatement* technique you choose to extrapolate the uncertainty in your results from the sample over the unseen part of the value of the population alone. 

#### Advantage

This technique will often result in a smaller sample size than is required if the technique is not applied.

#### Cautionary note

This technique requires the assumption that the population taints are homogeneously distributed. You will need to verify this in a post-hoc assumption check (available in the evaluation stage) that tests the correlation between the size of the taints and the size of the transaction. If this correlation does not exist, the use of this efficiency technique is justified. If this correlation is present, you cannot use this technique and will need to evaluate your sample using a different method.Sampling objective: Testing against a performance materiality
==========================

Also called the upper error limit, the tolerable deviation rate, or the tolerable misstatement, the performance materiality is the amount established by the auditor below the normal materiality of the financial reports to decrease the probability that the aggregate of uncorrected and undetectable misstatements exceeds the materiality of financial reports as a whole. 

#### Statistical interpretation
In the statistical analysis, the performance materiality represents the upper bound of tolerable misstatement in the population to be tested. 

#### Effect on procedure
By testing against a performance materiality, you are able to plan a sample in order to collect evidence for or against the statement that the population as a whole does not contain misstatements that are considered material (i.e., are greater than the upper bound of tolerable misstatement).

#### When should you enable this objective?
You should enable this objective when you want to find out whether the population contains misstatements that are greater than a certain limit (the performance materiality) using a sample of the population.

#### Effect on sample size
A lower performance materiality will result in a higher required sample size. Vice versa, a higher performance materiality will result in a lower required sample size.The use of explanatory text
==========================

The explanatory text provides you with an in-depth explanation of each step in the audit process. The text will help you understand what the practical interpretation of your choices is, and helps translate your statistical findings to others (in common English). With explanatory text enabled, each figure will also have its own caption interpreting the result displayed in that graph.Aselecte Steekproef
==========================

De aselecte steekproef is de meest eenvoudige manier van selecteren. Met deze methode heeft iedere steekproefeenheid dezelfde kans om geselecteerd te worden. Een aselecte selectie van de grootte <i>n</i> van de steekproefeenheden wordt gekozen (met de toevalsgenerator beginwaarde optie) en de bijbehorende delen van de populatie worden meegenomen in de steekproef.

<img src="%HELP_FOLDER%/img/randomSampling.png" />Steekproefeenheden
==========================

Enige kennis van steekproefeenheden is vereist om steekproeven uit de populatie te nemen; steekproefeenheden zijn namelijk de fysieke representaties van de populatie die worden geïnspecteerd. Steekproefeenheden kunnen individuele <i>transacties</i> of <i>monetaire eenheden</i> zijn. Bij statistische steekproeven worden selectiekansen toegekend aan de steekproefeenheden van de populatie-elementen. De totale verzameling van alle steekproefeenheden met de toegekende selectiekans heet het steekproefkader. 

-------

Record sampling
==========================

### Procedure

Het steekproefkader voor record sampling is de totale verzameling van alle individuele populatie-eenheden in de populatie. Een steekproefeenheid voor record sampling is een individuele populatie eenheid of transactie, zoals 1 transactie uit een populatie van 1000 transacties. In dit geval wordt het geselecteerde populatie meegenomen in de steekproef. Bijvoorbeeld, als de tiende dollar op een bonnetje van koffiemelk wordt geselecteerd, wordt dit koffiemelk bonnetje geïnspecteerd en wordt de geïnspecteerde (ware) waarde, of volledige juistheid, bepaald voor evaluatie. 

### Implicaties voor selectie 

De bedoeling van record sampling is dat transacties dezelfde kans hebben om te worden geselecteerd, aangezien iedere transactie een individuele steekproefeenheid is. Daarom heeft deze manier van selecteren de voorkeur voor niet-financiële populaties en populaties waarin transacties met een inventarisatiewaarde van $0 voorkomen. Record sampling in niet-monetaire populaties wordt vaak vervolgd met een binaire evaluatie methode, terwijl dit in monetaire populaties vaak een directe/verschil/ratio/regressie schatting volgt.

### Valkuilen

Let op dat record sampling in monetaire populaties een gelijke kans toekent aan elke transactie, wat kan leiden tot een uiteindelijke selectie met maar een klein aandeel van de populatiewaarde. 
Sampling units
==========================

Sampling from the population requires knowledge of the sampling units; physical representations of the population that need to be audited. Sampling units can be individual <i>transactions</i> or <i>monetary units</i>. For statistical sampling, sampling probabilities are assignned to the sampling units of the population elements. The total collection of all sampling units with an assigned selection probability is named the sampling frame.

-------

Record sampling
==========================

### Procedure

The sampling frame for record sampling is the total collection of individual population entries. A sampling unit for record sampling is an individual population entry or transaction, like 1 transaction within a population of 1000 transactions. In this case, the population entry that is selected is included in the sample. For example, if the population entry for a receipt of coffee milk is selected the coffee receipt is audited and its audit (true) value, or complete correctness, is determined for evaluation. 

### Implications for selection

The implication of a record sampling scheme is that transactions are selected with an equal probability, since each transaction is an individual sampling unit. It is therefore the preferred selection type for non-monetary populations and populations where there are booked transactions with a value of $0. Record sampling in a non-monetary population is often followed up with a binary evaluation method, while in monetary populations it is often followed up with an direct/difference/ratio/regression estimator.

### Pitfalls

Be aware that record sampling in monetary populations assigns equal probabilities to each transaction, implying that the sampling procedure may result in a final selection that only covers a small proportion of the population value.Random sampling
==========================

Random sampling is the most straight-forward sampling method. With this method, every sampling unit has the same probability of being selected as every other sampling unit. A random selection of size <i>n</i> of the sampling units is then taken (using the seed option) and the corresponding population entries are included in the sample.

<img src="%HELP_FOLDER%/img/randomSampling.png" />
Sampling objective: Obtaining a required minimum precision
==========================

The precision is a measure of how much certainty there is in the estimate of the misstatement from testing a particular characteristic of a sample at a given level of sampling risk.

#### Statistical interpretation
In the statistical analysis, the precision is represented by the difference between the estimated most likely error and the upper bound on the misstatement. 

#### Effect on procedure
By enabling this sampling objective, you are be able to plan a sample so that the difference between the estimated most likely error and the upper bound on the misstatement is reduced to a minimum percentage. 

#### When should you enable this objective?
You should enable this objective if you are interested in making an estimate of the population misstatement with a certain accuracy.

#### Effect on sample size
A lower minimum required precision will result in a higher required sample size. Vice versa, a higher minimum required precision will result in a lower required sample size.Cell sampling
==========================

In the case of cell sampling, the population is divided into a set of intervals, the size <i>I</i> of which is computed by dividing the size of the sampling frame by the sample size. Within each interval, a sampling unit is selected by randomly drawing a number within each interval. This causes the space between the intervals <i>i</i> to vary. In a monetary unit sampling context, the cell sampling method has the property that all transactions larger than twice the interval will always be included in the final sample.

<img src="%HELP_FOLDER%/img/cellSampling.png" />