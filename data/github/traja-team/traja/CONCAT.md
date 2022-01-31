---
title: 'Traja: A Python toolbox for animal trajectory analysis'
tags:
  - Python
  - animal behavior
  - trajectory
  - multivariate time series
  - neuroscience
authors:
  - name: Justin Shenk
    orcid: 0000-0002-0664-7337
    affiliation: "1, 2"
  - name: Wolf Byttner
    affiliation: 3
    orcid: 0000-0002-9525-9730
  - name: Saranraj Nambusubramaniyan
    affiliation: 1
    orcid: 0000-0002-7314-0261
  - name: Alexander Zoeller
    affiliation: 4
    orcid: 0000-0002-4043-3420
affiliations:
 - name: VisioLab, Berlin, Germany
   index: 1
 - name: Radboud University, Nijmegen, Netherlands
   index: 2
 - name: Rapid Health, London, England, United Kingdom
   index: 3
 - name: Independent researcher
   index: 4
date: 4 June 2021
bibliography: paper.bib
---

# Summary
There are generally four categories of trajectory data: mobility of people, mobility of transportation vehicles, mobility of animals, and mobility of natural phenomena [@zheng-trajectory-2015]. Animal tracking is important for fields as diverse as ethology, optimal foraging theory, and neuroscience. Mouse behavior, for example, is a widely studied in biomedical and brain research in models of neurological disease such as stroke.[^1] 

Several tools exist which allow analyzing mouse locomotion. Tools such as Ethovision [@spink_ethovision_2001] and DeepLabCut [@Mathisetal2018] allow converting video data to pose coordinates, which can further be analyzed by other open source tools. DLCAnalyzer[^2] provides a collection of R scripts for analyzing positional data, in particular visualizing, classifying and plotting movement. B-SOiD [@Hsu770271] allows unsupervised clustering of behaviors, extracted from the pose coordinate outputs of DeepLabCut. SimBA [@sgoldenlab_2021_4521178] provides several classifiers and tools for behavioral analysis in video streams in a Windows-based graphical user interface (GUI) application.

These tools are primarily useful for video data, which is not available for the majority of animal studies. For example, video monitoring of home cage mouse data is impractical today due to housing space constraints. Researchers using Python working with non-visual animal tracking data sources are not able to fully leverage these tools. Thus, a tool that supports modeling in the language of state-of-the-art predictive models  [@amirian_social_2019; @liang_peeking_2019; @chandra_traphic_2019], and which provides animal researchers with a high-level API for multivariate time series feature extraction, modeling and visualization is needed.

Traja is a Python package for statistical analysis and computational modelling of trajectories. Traja extends the familiar pandas [@mckinney-proc-scipy-2010; @reback2020pandas] methods by providing a pandas accessor to the `df.traja` namespace upon import. The API for Traja was designed to provide an object-oriented and user-friendly interface to common methods in analysis and visualization of animal trajectories. Traja also interfaces well with relevant spatial analysis packages in R (e.g., trajr [@mclean_trajr:_2018] and adehabitat [@adehabitat]), Shapely [@shapely], and MovingPandas [@graser_movingpandas_2019] allowing rapid prototyping and comparison of relevant methods in Python. A comprehensive source of documentation is provided on the home page
([http://traja.readthedocs.io](traja.readthedocs.io)).

## Statement of Need
The data used in this project includes animal trajectory data provided by [http://www.tecniplast.it](Tecniplast S.p.A.), manufacturer of laboratory animal equipment based in Varese, Italy, and Radboud University, Nijmegen, Netherlands. Tecniplast provided the mouse locomotion data collected with their Digital Ventilated Cages (DVC). The extracted coordinates of the mice requires further analysis with external tools. Due to lack of access to equipment, mouse home cage data is rather difficult to collect and analyze, thus few studies have been done on home cage data. Furthermore, researchers who are interested in developing novel algorithms must implement from scratch much of the computational and algorithmic infrastructure for analysis and visualization. By packaging a library that is particularly useful for animal locomotion analysis, future researchers can benefit from access to a high-level interface and clearly documented methods for their work.

Other toolkits for animal behavioral analysis either rely on visual data [@Mathisetal2018; @vivek_hari_sridhar_2017_1134016] to estimate the pose of animals or are limited to the R programming language [@mclean_trajr:_2018]. Prototyping analytical approaches and exploratory data analysis is furthered by access to a wide range of methods which existing libraries do not provide. Python is the *de facto* language for machine learning and data science programming, thus a toolkit in Python which provides methods for prototyping multivariate time series data analysis and deep neural network modeling is needed.

## Overview of the Library
Traja targets Python because of its popularity with data scientists. The library leverages the powerful pandas library [@mckinney-proc-scipy-2010], while adding methods specifically for trajectory analysis. When importing Traja, the Traja namespace registers itself within the pandas dataframe namespace via `df.traja`.

The software is structured into three parts. These provide functionality to transform, analyse and visualize trajectories. Full details are available at <https://traja.readthedocs.io/>.  The `trajectory` module provides analytical and preprocessing functionalities. The `models` subpackage provides both traditional and neural network-based tools to determine trajectory properties. The `plotting` module allows visualizing trajectories in various ways. 

Data, e.g., x and y coordinates, are stored as one-dimensional labelled arrays as instances of the pandas native `Series` class. Further, subclassing the pandas `DataFrame` allows providing an API that mirrors the pandas API which is familiar to most data scientists, thus reducing the barrier for entry while providing methods and properties specific to trajectories for rapid prototyping.
Traja depends on Matplotlib [@Hunter:2007] and Seaborn [@Waskom2021] for plotting and NumPy [@harris2020array] for computation.

### Trajectory Data Sources
Trajectory data as time series can be extracted from a wide range of sources, including video processing tools as described above, GPS sensors for large animals or via home cage floor sensors, as described in the section below. The methods presented here are implemented for orthogonal coordinates *(x, y)* primarily to track animal centroids, however with some modification they could be extended to work in 3-dimensions and with body part locations as inputs. Traja is thus positioned at the end of the data scientist's chain of tools with the hope of supporting prototyping novel data processing approaches. A sample dataset of jaguar movement [@morato_jaguar_2018] is provided in the `traja.dataset` subpackage.

## Mouse Locomotion Data
The data samples presented here[^3] are in 2-dimensional location coordinates, reflecting the mouse home cage (25x12.5 cm) dimensions. Analytical methods relevant to 2D rectilinear analysis of highly constrained spatial coordinates are thus primarily considered.

High volume data like animal trajectories has an increased tendency to have missing data due to data collection issues or noise. Filling in the missing data values, referred to as _data imputation_, is achieved with a wide variety of statistical or learning-based methods. As previously observed, data science projects typically require at least _95%_ of the time to be spent on cleaning, pre-processing and managing the data [@bosch_engineering_2021]. Therefore, several methods relevant to preprocessing animal data are demonstrated throughout the following sections.

[^1]: The examples in this paper focus on animal motion, however it is useful for other domains. 

[^2]: <https://github.com/ETHZ-INS/DLCAnalyzer>

[^3]: This dataset has been collected for other studies of our laboratory [@shenk_automated_2020].

## Spatial Trajectory
A *spatial trajectory* is a trace generated by a moving object in geographical space. Trajectories are traditionally modelled as a sequence of spatial points like:

$$T_k = \{P_{k1}, P_{k2},...\}$$

where $P_{ki}(i\geq 1)$ is a point in the trajectory.

Generating spatial trajectory data via a random walk is possible by sampling from a distribution of angles and step sizes [@kareiva_analyzing_1983; @mclean_trajr:_2018]. A correlated random walk (Figure [1](#fig:generated){reference-type="ref" reference="fig:generated"}) is generated with `traja.generate`.

![Generation of a random walk[]{label="fig:generated"}](./images/generate.png){#fig:generated width=80%}

## Spatial Transformations
Transformation of trajectories can be useful for comparing trajectories from various geospatial coordinates, data compression, or simply for visualization purposes.

### Feature Scaling
Feature scaling is common practice for preprocessing data for machine learning [@grus_data_2015] and is essential for even application of methods to attributes. For example, a high dimensional feature vector $\mathbf{x} \in \mathbb{R}^n$ where some attributes are in $(0,100)$ and others are in $(-1,1)$ would lead to biases in the treatment of certain attributes. To limit the dynamic range for multiple data instances simultaneously, scaling is applied to a feature matrix $X = \{\mathbf{x_1}, \mathbf{x_2}, ..., \mathbf{x_N}\} \in \mathbb{R}^{n\times{N}}$, where $n$ is the number of instances.

**Min-Max Scaling** To guarantee that the algorithm applies equally to all attributes, the normalized feature matrix $\hat{X}$ is rescaled into range $(0,1)$ such that

$\hat{X} = \frac{X - X_{min}}{X_{max} - X_{min}}$.

**Standardization** The result of standardization is that the features will be rescaled to have the property of a standard normal distribution with $\mu = 0$ and $\sigma = 1$ where $\mu$ is the mean (average) of the data and $\sigma$ is the standard deviation from the mean. Standard scores (also known as **z**-scores are calculated such that

$z = \frac{x-\mu}{\sigma}$.

**Scaling** Scaling a trajectory is implemented for factor $f$ in `scale` where $f \in R: f \in (-\infty, +\infty)$.

### Rotation
Rotation of a 2D rectilinear trajectory is a coordinate transformation of orthonormal bases x and y at angle $\theta$ (in radians) around the origin defined by

$$\begin{bmatrix} x'\\y' \end{bmatrix} = \begin{bmatrix} cos\theta & i sin\theta\\ sin\theta & cos\theta \end{bmatrix} \begin{bmatrix} x\\y \end{bmatrix} $$

with angle $\theta$ where $\theta \in R : \theta \in [-180,180]$.

### Trip Grid
One strategy for compressing the representation of trajectories is binning the coordinates to produce an image as shown in Figure [2](#fig:tripgridalgo){reference-type="ref" reference="fig:tripgridalgo"}.

![Trip grid image generation from mouse
trajectory.](./images/trip_grid_algo.png){#fig:tripgridalgo width=100%}

Allowing computation on discrete variables rather than continuous ones has several advantages stemming from the ability to store trajectories in a more memory efficient form.[^4] The advantage is that computation is generally faster, more data can fit in memory in the case of complex models, and item noise can be reduced.

[^4]: In this experiment, for example, data can be reduced from single-precision floating point (32 bits) to 8-bit unsigned integer (*uint8*) format.

Creation of an $M * N$ grid allows mapping trajectory $T_k$ onto uniform
grid cells. Generalizing the nomenclature of [@wang_modeling_2017] to rectangular grids, $C_{mn}(1\leq{m}\leq M; 1\leq{n}\leq{N})$ denotes the cell in row $m$ and column $n$ of the grid. Each point $P_{ki}$ is assigned to a cell $C(m,n)$. The result is a two-dimensional image $M*N$ image $I_k$, where the value of pixel $I_k(m,n)(1\leq{m,n}\leq{M})$ indicates the relative number of points assigned to cell $C_{mn}$. Partionining of spatial position into separate grid cells is often followed by generation of hidden Markov models [@jeung_mining_2007] (see below) or visualization of heat maps (Figure [3](#fig:heatmap){reference-type="ref" reference="fig:heatmap"}).

![Visualization of heat map from bins generated with `df.trip_grid`. Note regularly spaced artifacts (bright yellow) in this sample due to a bias in the sensor data interpolation. This type of noise can be minimized by thresholding or using a logarithmic scale, as shown above.[]{label="fig:heatmap"}](./images/tripgrid.png){#fig:heatmap width=50%}

### Smoothing
Smoothing a trajectory can also be achieved with Traja using Savitzky-Golay filtering with `smooth_sg` [@savitzky_smoothing_1964].

## Resampling and Rediscretizing
Trajectories can be resampled by time or rediscretized by an arbitrary step length. This can be useful for aligning trajectories from various data sources and sampling rates or reducing the number of data points to improve computational efficiency. Care must be taken to select a time interval which maintains information on the significant behavior. If the minimal time interval observed is selected for the points, calculations will be computationally intractable for some systems. If too large of an interval is selected, we will fail to capture changes relevant to the target behavior in the data.

Resampling by time is performed with `resample_time` (Figure [4](#fig:sample){reference-type="ref" reference="fig:sample"}). Rediscretizing by step length is performed with `rediscretize`.

![Resampling trajectories by different time scales is performed with `resample_time`.[]{label="fig:sample"}](./images/sample_rate.png){#fig:step width=80%}

For example, the Fortasyn dataset [@shenk_automated_2020] demonstrated in this paper was sampled at 4 Hz and converted to single-precision floating point data. Pandas dataframes store this data in 4 bytes, thus there are approximately 4.15 MB[^5] bytes required to store data for x and y dimensions plus an index reference for a single day. In the case of [@shenk_automated_2020], 24 mice were observed over 35 days. This translates to 3.4 GB ($10^9$) of storage capacity for the uncompressed datasets prior to feature engineering. Thus resampling can be a useful way to reduce the memory footprint for memory constrained processes that have to fit into a standard laptop with 8 GB memory space. A demonstration of how reduction in precision for trajectory data analysis is provided in Figure [4](#fig:step){reference-type="ref" reference="fig:step"}, as applied to a sample from the Fortasyn experiment [@shenk_automated_2020]. Broad effects such as cage crossings, for example, can still be identified while downsampling data to a lower frequency, such as 0.1 Hz, reducing the memory footprint by a factor of 40 (4 Hz/0.1 Hz) and providing significant speedups for processing.

## Movement Analysis
Traja includes traditional as well as advanced methods for trajectory analysis.

### Distance traveled
Distance traveled is a common metric in animal studies - it accounts for the total distance covered by the animal within a given time interval. The distance traveled is typically quantified by summing the square straight-line displacement between discretely sampled trajectories [@rowcliffe_bias_2012; @solla_eliminating_1999]. Alternative distance metrics for the case of animal tracking are discussed in [@noonan_scale-insensitive_2019].

Let $p(t) = [p_x(t), p_y(t)]$ be a $2\times 1$ vector of coordinates on the ground representing the position of the animal at time t. Then, the distance traveled within the time interval $t_1$ and $t_2$ can be computed as a sum of step-wise Euclidean distances

$$p(t_1,t_2) = \Sigma^{t_2}_{t=t_1+1} d(t),$$

where
$$d(t) = \sqrt{(p_x(t) -p_x(t-1))^2 + (p_y(t) - p_y(t-1))^2}  $$

is the Euclidean distance between two positions in adjacent time samples.

[^5]: 4 x 4 Hz x 60 seconds x 60 minutes x 24 hours x 3 features (x, y, and time)

![Velocity histogram from one day of mouse activity.[]{label="fig:velocity-hist"}](./images/velocitylog.png){#fig:velocity-hist width=50%}

### Speed
Speed or velocity is the first derivative of centroids with respect to time. Peak velocity in a home cage environment is perhaps less interesting than a distribution of velocity observations, as in Figure [5](#fig:velocity-hist){reference-type="ref" reference="fig:velocity-hist"}. Additionally, noise can be eliminated from velocity calculations by using a minimal distance moved threshold, as demonstrated in [@shenk_automated_2020]. This allows identifying broad-scale behaviors such as cage crossings.

### Turn Angles
Turn angles are the angle between the movement vectors of two consecutive samples. They can be calculated with `calc_turn_angles`.

### Laterality
Laterality is the preference for left or right turning and a *laterality index* is defined as:
$$LI = \frac{RT}{LT + RT} $$

where RT is the number of right turns observed and LT is the number of left turns observed. Turns are counted within a left turn angle $\in$ ($\theta$, 90) and right turn angle $\in(-\theta,-90)$. A turn is considered to have a minimal step length.

## Periodicity
Periodic behaviors are a consequence of the circadian rhythm as well as observing expression of underlying cognitive traits. Some basic implementations of periodic analysis of mouse cage data are presented.

### Autocorrelation
Autocorrelation is the correlation of a signal with a delayed copy of itself as a function of the decay. Basically, it is similarity of observations as a function of the time lag between them.
An example is shown in Figure [6](#fig:autocorrelation){reference-type="ref" reference="fig:autocorrelation"}.

![Autocorrelation of the y-dimension reveals daily (1440 minutes) periodic behavior[]{label="fig:autocorrelation"}](./images/autocorrelation_E1.png){#fig:autocorrelation width=80%}

### Power Spectrum
Power spectrum of a time series signal can be estimated (Figure [7](#fig:powerspectrum){reference-type="ref" reference="fig:powerspectrum"}). This is useful for analyzing signals, for example, the influence of neuromotor noise on delays in hand movement [@van_galen_effects_1990].

![Power Spectral Density. One day of activity reveals fairly smooth power spectral density.[]{label="fig:powerspectrum"}](./images/spectrum.png){#fig:powerspectrum width=70%}

## Algorithms and Statistical Models

### Machine Learning for Time Series Data
Machine learning methods enable researchers to solve tasks computationally without explicit instructions by detecting patterns or relying on inference. Thus they are particularly relevant for data exploration of high volume datasets such as spatial trajectories and other multivariate time series.

### Principal Component Analysis
Principal Component Analysis projects the data into a linear subspace with a minimum loss of information by multiplying the data by the eigenvectors of the covariance matrix.

![PCA of Fortasyn trajectory data. Daily trajectories (day and night)
were binned into 8x8 grids before applying
PCA.[]{label="fig:pca"}](./images/pca_fortasyn-period.png){#fig:pca
width=80%}

This requires converting the trajectory to a trip grid (see Figure [2(#fig:tripgridalgo){reference-type="ref" reference="fig:tripgridalgo"}]) and performing PCA on the grid in 2D (Figure [8](#fig:pca){reference-type="ref" reference="fig:pca"}) or 3D (Figure [9](#fig:3dpca){reference-type="ref" reference="fig:3dpca"}). Structure in the data is visible if light and dark time periods are compared.

![3D PCA of Fortasyn trajectory data. Daily trajectories (day and night)
were binned into 8x8 grids before applying
PCA.[]{label="fig:3dpca"}](./images/pca_fortasyn-period-3d.png){#fig:3dpca
width=80%}

### Clustering
Clustering of trajectories is an extensive topic with applications in geospatial data, vehicle and pedestrian classification, as well as molecular identification. K-means clustering is an iterative unsupervised learning method that assigns a label to data points based on a distance function [@bishop_pattern_2006] (Figure [10](#fig:kmeans){reference-type="ref" reference="fig:3dpca"}).

![K-means clustering on the results of the PCA shown above reveals a high accuracy
of classification, with a few errors. Cluster labels are generated by
the model.[]{label="fig:kmeans"}](./images/kmeans_pca-fortasyn.png){#fig:kmeans
width=80%}

### Hierarchical Agglomerative Clustering
Clustering spatial trajectories has broad applications for behavioral research, including unsupervised phenotyping [@huang_mapping_2020]. For mice, hierarchical agglomerative clustering can also be used to identify similarities between groups, for example periodic activity and location visit frequency [@clustering_mice].

### Gaussian Processes
Gaussian Processes is a non-parametric method which can be used to model spatial trajectories. This method is not currently implemented in Traja
and is thus outside the scope of the current paper, however the interested reader is directed to the excellent text on Gaussian processes by Rasmussen and Williams [@rasmussen_gaussian_2006] for a complete reference and [@cox_gaussian_2012] for an application to spatial trajectories.

## Other Methods

### Fractal Methods
Fractal (i.e. multiscale) methods are useful for analyzing transitions and clustering in trajectories. For example, search trajectories such as eye movement, hand-eye coordination, and foraging can be analyzed by quantifying the spatial distribution or nesting of temporal point processes using spatial Allen Factor analysis [@kerster_spatial_2016; @huette_drawing_2013]. 

Recurrence plots and derivative recurrence factor analysis can be applied to trajectories to identify multiscale temporal processes to study transition or nonlinear parameters in a system, such as postural fluctuation  [@ross_influence_2016] and synchrony [@shockley] in humans and to movement of animals such as ants [@neves_recurrence_2017] and bees [@ayers]. These methods are not yet implemented in Traja, but are planned for a future release.

### Graph Models
A graph is a pair $G = (V, E)$ comprising a set of vertices and a set of connecting edges. A probabilistic graphical model of a spatial occupancy grid can be used to identify probabilities of state transitions between nodes. A basic example is given with hidden Markov models below.

![Transition matrix. Rows and columns are flattened histogram of a grid
20 cells high and 10 cells wide. Spatially adjacent grid cells are
visible at a spacing of -11, -10, -9, 1, 10, and 11 cells from the
diagonal. The intensity of pixels in the diagonal represents relative
likelihood to stay in the same
position.[]{label="fig:transitionmatrix"}](./images/transition_matrix.png){#fig:transitionmatrix
width=60%}

### Hidden Markov Models
Transition probabilities are most commonly modelled with Hidden Markov Models (HMM) because of their ability to capture spatial and temporal dependencies. A recent introduction to these methods is available provided by [@patterson_statistical_2017]. HMMs have successfully been used to analyze movement of caribou [@franke_analysis_2004], fruit flies [@holzmann_hidden_2006], and tuna [@patterson_migration_2018], among others. Trajectories are typically modelled as bivariate time series consisting of step length and turn angle, regularly spaced in time.

Traja implements the rectangular spatial grid version of HMM with transitions.

The probability of transition from each cell to another cell is stored as a probability within the transition matrix. This can visualized as a heatmap and plotted with `plot_transition_matrix` (Figure [11](#fig:transitionmatrix){reference-type="ref" reference="fig:transitionmatrix"}).

### Convex Hull
The convex hull of a subtrajectory is the set $X$ of points in the Euclidean plane that is the smallest convex set to include $X$. For computational efficiency, a geometric k-simplex can be plotted covering the convex hull by converting to a Shapely object and using Shapely’s `convex_hull` method.

### Recurrent Neural Networks
In recent years, deep learning has transformed the field of machine learning. For example, the current state of the art models for a wide range of tasks, including computer vision, speech to text, and pedestrian trajectory prediction, are achieved with deep neural networks. Neural networks are essentially sequences of matrix operations and elementwise function application based on a collection of computing units known as nodes or neurons. These units perform operations, such as matrix multiplication on input features of a dataset, followed by backpropagation of errors, to identify parameters useful for approximating a function.

![Neural network architectures available in Traja](./images/dnns.jpg){width=100%}

Recurrent Neural Networks (RNNs) are a special type of Neural Networks that use
a state $S(t_{i-1})$ from the previous timestep $t_{i-1}$ alongside X($t_i$) as input. They output a prediction $Y(t_i)$ and a new state $S(t_i)$ at every step. Utilising previous states makes RNNs particularly good at analyzing time series like trajectories, since they can process arbitrarily long inputs. They remember information from previous time steps $X(t_{i-k}), ..., X(t_{i-1})$ when processing the current time step $X(t_i)$.

Trajectory prediction lets researchers forecast the location and trajectory of animals [@wijeyakulasuriya_machine_2020]. Where this technique works well, it is also a sign that the trajectory is highly regular and, fundamentally, follows certain rules and patterns. When tracking an animal live, it would also let researchers predict when it will arrive at a particular location, or where it will go, letting them rig cameras and other equipment ahead of time.

A particularly interesting type of RNN is the Long Short Term Memory (LSTM) architecture. Their layers use stacks of units, each with two hidden variables - one that quickly discards old states and one that slowly does so - to consider relevant information from previous time steps. They can thus look at a trajectory and determine a property of the animal – whether it is sick or injured, say – something that is time-consuming and difficult to do by hand. They can also predict future time steps based on past ones, letting researchers estimate where the animal will go next. LSTMs can also classify trajectories, determining whether a trajectory comes from an animal belonging in a specific category. This lets researchers determine how a controlled or semi-controlled variable (e.g., pregnancy) changes the movement pattern of an animal.

Traja implements neural networks by extending the widely used open source machine learning library PyTorch [@pytorch], primarily developed by Facebook AI Research Group. Traja allows framework-agnostic modeling through data loaders designed for time series. In addition, the Traja package comes with several predefined model architectures which can be configured according to the user’s requirements.

Because RNNs work with time series, the trajectories require special handling. The `traja.dataset.MultiModalDataLoader` efficiently groups subsequent samples and into series and splits these series into training and test data. It represents a Python iterable over the dataset and extends the PyTorch `DataLoader` class, with support for

- random, weighted sampling,
- data scaling,
- data shuffling,
- train/validation/test split.

`MultiModalDataLoader` accepts several important configuration parameters and
allows batched sampling of the data. The two constructor arguments `n_past` and
`n_future` specify the number of samples that the network will be shown and the number that the network will have to guess, respectively. `batch_size` is generally in the dozens and is used to regularise the network.

The RNNs also need to be trained - this is done by the high-level `Trainer` class below. It performs nonlinear optimisation with a Stochastic Gradient Descent-like algorithm. The `Trainer` class by default implements the Huber loss function [@huber_robust_1964], also known as smooth $L_1$ loss, which is a loss function commonly used in robust regression:

$$L_{\delta} (a) = \begin{cases}
 \frac{1}{2}{a^2}                   & \text{for } |a| \le \delta, \\
 \delta (|a| - \frac{1}{2}\delta), & \text{otherwise.}
\end{cases}$$

In comparison to mean-squared error loss, Huber loss is less sensitive to outliers in data: it is quadratic for small values of a, and linear for large values. It extends the PyTorch `SmoothL1Loss` class, where the $d$ parameter is set to 1.[^6] A common optimization algorithm is ADAM and is Traja’s default, but several others are provided as well. Although training with only a CPU is possible, a GPU can provide a $40-100x$ speedup [@Arpteg2018SoftwareEC].

[^6]: [https://pytorch.org/docs/stable/generated/torch.nn.SmoothL1Loss.html](https://pytorch.org/docs/stable/generated/torch.nn.SmoothL1Loss.html)

### Recurrent Autoencoder Networks
Traja can also train autoencoders to either predict the future position of a track or classify the track into a number of categories. Autoencoders embed the time series into a time-invariant latent space, allowing representation of each trajectory or sub-trajectory as a vector. A class of well-separated trajectories would then be restricted to a region of the latent space. The technique is similar to Word2vec [@word2vec], where words are converted to a 100+ dimensional vector. In this approach, forecasting and classification are both preceded by training the data in an autoencoder, which learns an efficient representation of the data for further computation of the target function.

Traja allows training a classifier that works directly on the latent space output - since each class of trajectories converges to a distinct region in the latent space, this technique is often superior to classifying the trajectory itself. Traja trains classifiers for both Autoencoder-style and Variational Autoencoder-style RNNs. When investigating whether animal behavior has changed, or whether two experimental categories of animals behave differently, this unstructured data mining can suggest fruitful avenues for investigation.

# References
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

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

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

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
Contributing to Traja
=====================

Traja is a research library. Functionality must therefore be both
cutting-edge and reliable. Traja is part of a wider project to
increase collaboration in research, through the adoption of
open-source contribution models. It is our hope that traja
remains accessible to researchers and help them do higher-quality
research.

Current status
--------------

Traja is currently undergoing active development with approximately
60 % of features present. Significant interface changes are still
possible, however we avoid these unless absolutely necessary.

The work is currently focussed on reaching version 1.0 with feature
completeness and 95 % test coverage.

The following features are required for feature completeness:

* Latent space visualisers
   * Eigenspace-based
   * Colour-coded to visualise evolution over time
   * Delay coordinate embeddings
* State-space visualisers
* Additional encoder and decoder options in AE and VAE models
   * MLP
   * 1d convolution
* Pituitary gland example dataset
* Regression output visualisers
* VAE GAN models
* Additional VAE latent-space shapes
   * Uniform
   * A shape that works for periodic trajectories (Torus?)
* Delay coordinate embeddings
   * Persistent homology diagrams of the embeddings
* Automatic code formatter
* Tutorials
   * Find time of day based on activity
   * `Recover parameters from Pituitary ODE <https://colab.research.google.com/drive/1Fc8Z6LtjMBk20QPPzJt8uDAxWu5dz_1Y?usp=sharing>`_
   * `Predict stock prices with LSTMs <https://colab.research.google.com/drive/1TAW9Lv0Sm9g8YRgYXaiHIZBqnGtGQn1H?usp=sharing>`_

How to contribute
-----------------

Traja welcomes contributions! To get started, pick up any issue
labeled with `good first issue`! Alternatively you can read some
background material or try a tutorial.

Testing and code quality
------------------------

Since Traja is a library, we strive for sensible tests achieving a
high level of code coverage. Future commits are required to maintain
or improve code quality, test coverage. To aid in this, Travis runs
automated tests on each pull request. Additionally, we run codecov
on PRs.

Background material
-------------------

This is a collection of papers and resources that explain the
main problems we are working on with Traja.

Analysis of mice that have suffered a stroke:

    @article{10.3389/fnins.2020.00518,
      author={Justin Shenk and
              Klara J. Lohkamp and
              Maximilian Wiesmann and
              Amanda J. Kiliaan},
      title={Automated Analysis of Stroke Mouse Trajectory Data With Traja},
      journal={Frontiers in Neuroscience},
      volume={14},
      pages={518},
      year={2020},
      url={https://www.frontiersin.org/article/10.3389/fnins.2020.00518},
      doi={10.3389/fnins.2020.00518},
      issn={1662-453X},
    }


Understanding the parameter space of the pituitary gland ODE (https://www.math.fsu.edu/~bertram/papers/bursting/JCNS_16.pdf):


    @article{10.1007/s10827-016-0600-1,
      author = {Fletcher, Patrick and Bertram, Richard and Tabak, Joel},
      title = {From Global to Local: Exploring the Relationship between Parameters and Behaviors in Models of Electrical Excitability},
      year = {2016},
      publisher = {Springer-Verlag},
      address = {Berlin, Heidelberg},
      volume = {40},
      number = {3},
      issn = {0929-5313},
      url = {https://doi.org/10.1007/s10827-016-0600-1},
      doi = {10.1007/s10827-016-0600-1},
      journal = {J. Comput. Neurosci.},
      month = June,
      pages = {331–345},
    }


Style guide
-----------
TODO
Traja |Python-ver| |Travis| |PyPI| |Conda| |RTD| |Gitter| |Black| |License| |Binder| |Codecov| |DOI| |JOSS|
===========================================================================================================

|Colab|

.. |Python-ver| image:: https://img.shields.io/badge/python-3.6+-blue.svg
    :target: https://www.python.org/downloads/release/python-360/
    :alt: Python 3.6+

.. |Travis| image:: https://travis-ci.org/traja-team/traja.svg?branch=master
    :target: https://travis-ci.org/traja-team/traja

.. |PyPI| image:: https://badge.fury.io/py/traja.svg
    :target: https://badge.fury.io/py/traja

.. |Conda| image:: https://img.shields.io/conda/vn/conda-forge/traja.svg
    :target: https://anaconda.org/conda-forge/traja

.. |Gitter| image:: https://badges.gitter.im/traja-chat/community.svg
    :target: https://gitter.im/traja-chat/community

.. |RTD| image:: https://readthedocs.org/projects/traja/badge/?version=latest
    :target: https://traja.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/ambv/black

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://opensource.org/licenses/MIT
    :alt: License: MIT

.. |Binder| image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/justinshenk/traja/master?filepath=demo.ipynb

.. |Codecov| image:: https://codecov.io/gh/traja-team/traja/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/traja-team/traja

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5069231.svg
   :target: https://doi.org/10.5281/zenodo.5069231

.. |Colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/justinshenk/traja/blob/master/demo.ipynb
   
.. |JOSS| image:: https://joss.theoj.org/papers/0f25dc08671e0ec54714f09597d116cb/status.svg
   :target: https://joss.theoj.org/papers/0f25dc08671e0ec54714f09597d116cb

Traja is a Python library for trajectory analysis. It extends the capability of
pandas DataFrame specific for animal trajectory analysis in 2D, and provides
convenient interfaces to other geometric analysis packages (eg, R and shapely).

Introduction
------------

The traja Python package is a toolkit for the numerical characterization
and analysis of the trajectories of moving animals. Trajectory analysis
is applicable in fields as diverse as optimal foraging theory,
migration, and behavioral mimicry (e.g. for verifying similarities in
locomotion). A trajectory is simply a record of the path followed by a
moving animal. Traja operates on trajectories in the form of a series of
locations (as x, y coordinates) with times. Trajectories may be obtained
by any method which provides this information, including manual
tracking, radio telemetry, GPS tracking, and motion tracking from
videos.

The goal of this package (and this document) is to aid biological
researchers, who may not have extensive experience with Python, to
analyze trajectories without being restricted by a limited knowledge of
Python or programming. However, a basic understanding of Python is
useful.

If you use traja in your publications, please cite the repo 

.. code-block::

    @software{justin_shenk_2019_3237827,
      author       = {Justin Shenk and
                      the Traja development team},
      title        = {justinshenk/traja},
      month        = jun,
      year         = 2019,
      publisher    = {Zenodo},
      version      = {latest},
      doi          = {10.5281/zenodo.3237827},
      url          = {https://doi.org/10.5281/zenodo.3237827}
    }


Installation and setup
----------------------

To install traja with conda, run

``conda install -c conda-forge traja``

or with pip

``pip install traja``.

Import traja into your Python script or via the Python command-line with
``import traja``.

Trajectories with traja
-----------------------

Traja stores trajectories in pandas DataFrames, allowing any pandas
functions to be used.

Load trajectory with x, y and time coordinates:

.. code-block:: python

    import traja

    df = traja.read_file('coords.csv')

Once a DataFrame is loaded, use the ``.traja`` accessor to access the
visualization and analysis methods:

.. code-block:: python

    df.traja.plot(title='Cage trajectory')


Analyze Trajectory
------------------

.. csv-table:: The following functions are available via ``traja.trajectory.[method]``
   :header: "Function", "Description"
   :widths: 30, 80
   
   "``calc_derivatives``", "Calculate derivatives of x, y values "
   "``calc_turn_angles``", "Calculate turn angles with regard to x-axis "
   "``transitions``", "Calculate first-order Markov model for transitions between grid bins"
   "``generate``", "Generate random walk"
   "``resample_time``", "Resample to consistent step_time intervals"
   "``rediscretize_points``", "Rediscretize points to given step length"
   
For up-to-date documentation, see `https://traja.readthedocs.io <https://traja.readthedocs.io>`_.

Random walk
-----------

Generate random walks with

.. code-block:: python

    df = traja.generate(n=1000, step_length=2)
    df.traja.plot()

.. image:: https://raw.githubusercontent.com/justinshenk/traja/master/docs/source/_static/walk_screenshot.png
   :alt: walk\_screenshot.png


Resample time
-------------
``traja.trajectory.resample_time`` allows resampling trajectories by a ``step_time``.


Flow Plotting
-------------

.. code-block:: python

    df = traja.generate()
    traja.plot_surface(df)

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_001.png
   :alt: 3D plot

.. code-block:: python

    traja.plot_quiver(df, bins=32)

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_002.png
   :alt: quiver plot

.. code-block:: python

    traja.plot_contour(df, filled=False, quiver=False, bins=32)

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_003.png
   :alt: contour plot

.. code-block:: python

    traja.plot_contour(df, filled=False, quiver=False, bins=32)

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_004.png
   :alt: contour plot filled

.. code-block:: python

    traja.plot_contour(df, bins=32, contourfplot_kws={'cmap':'coolwarm'})

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_005.png
   :alt: streamplot

Acknowledgements
----------------

traja code implementation and analytical methods (particularly
``rediscretize_points``) are heavily inspired by Jim McLean's R package
`trajr <https://github.com/JimMcL/trajr>`__. Many thanks to Jim for his
feedback.
Clustering and Dimensionality Reduction
=======================================

Clustering Trajectories
-----------------------

Trajectories can be clustered using :func:`traja.plotting.plot_clustermap`.

Colors corresponding to each trajectory can be specified with the ``colors`` argument.

.. autofunction:: traja.plotting.plot_clustermap

PCA
---

Prinicipal component analysis can be used to cluster trajectories based on grid cell occupancy. 
PCA is computed by converting the trajectory to a trip grid (see :meth:`traja.plotting.trip_grid`) followed by PCA (:class:`sklearn.decomposition.PCA`).

.. autofunction:: traja.plotting.plot_pcaContributing to traja
=====================

(Contribution guidelines largely copied from `geopandas <https://geopandas.readthedocs.io/en/latest/contributing.html>`_)

Overview
--------

Contributions to traja are very welcome.  They are likely to
be accepted more quickly if they follow these guidelines.

At this stage of traja development, the priorities are to define a
simple, usable, and stable API and to have clean, maintainable,
readable code. Performance matters, but not at the expense of those
goals.

In general, traja follows the conventions of the pandas project
where applicable.

In particular, when submitting a pull request:

- All existing tests should pass.  Please make sure that the test
  suite passes, both locally and on
  `Travis CI <https://travis-ci.org/justinshenk/traja>`_.  Status on
  Travis will be visible on a pull request.  If you want to enable
  Travis CI on your own fork, please read the pandas guidelines link
  above or the
  `getting started docs <http://about.travis-ci.org/docs/user/getting-started/>`_.

- New functionality should include tests.  Please write reasonable
  tests for your code and make sure that they pass on your pull request.

- Classes, methods, functions, etc. should have docstrings.  The first
  line of a docstring should be a standalone summary.  Parameters and
  return values should be ducumented explicitly.

- traja supports python 3 (3.6+).  Use modern python idioms when possible.

- Follow PEP 8 when possible.

- Imports should be grouped with standard library imports first,
  3rd-party libraries next, and traja imports third.  Within each
  grouping, imports should be alphabetized.  Always use absolute
  imports when possible, and explicit relative imports for local
  imports when necessary in tests.


Seven Steps for Contributing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are seven basic steps to contributing to *traja*:

1) Fork the *traja* git repository
2) Create a development environment
3) Install *traja* dependencies
4) Make a ``development`` build of *traja*
5) Make changes to code and add tests
6) Update the documentation
7) Submit a Pull Request

Each of these 7 steps is detailed below.


1) Forking the *traja* repository using Git
------------------------------------------------

To the new user, working with Git is one of the more daunting aspects of contributing to *traja**.
It can very quickly become overwhelming, but sticking to the guidelines below will help keep the process
straightforward and mostly trouble free.  As always, if you are having difficulties please
feel free to ask for help.

The code is hosted on `GitHub <https://github.com/justinshenk/traja>`_. To
contribute you will need to sign up for a `free GitHub account
<https://github.com/signup/free>`_. We use `Git <http://git-scm.com/>`_ for
version control to allow many people to work together on the project.

Some great resources for learning Git:

* Software Carpentry's `Git Tutorial <http://swcarpentry.github.io/git-novice/>`_
* `Atlassian <https://www.atlassian.com/git/tutorials/what-is-version-control>`_
* the `GitHub help pages <http://help.github.com/>`_.
* Matthew Brett's `Pydagogue <http://matthew-brett.github.com/pydagogue/>`_.

Getting started with Git
~~~~~~~~~~~~~~~~~~~~~~~~~

`GitHub has instructions <http://help.github.com/set-up-git-redirect>`__ for installing git,
setting up your SSH key, and configuring git.  All these steps need to be completed before
you can work seamlessly between your local repository and GitHub.

.. _contributing.forking:

Forking
~~~~~~~~

You will need your own fork to work on the code. Go to the `traja project
page <https://github.com/traja-team/traja>`_ and hit the ``Fork`` button. You will
want to clone your fork to your machine::

    git clone git@github.com:your-user-name/traja.git traja-yourname
    cd traja-yourname
    git remote add upstream git://github.com/traja-team/traja.git

This creates the directory `traja-yourname` and connects your repository to
the upstream (main project) *traja* repository.

The testing suite will run automatically on Travis-CI once your pull request is
submitted.  However, if you wish to run the test suite on a branch prior to
submitting the pull request, then Travis-CI needs to be hooked up to your
GitHub repository.  Instructions for doing so are `here
<http://about.travis-ci.org/docs/user/getting-started/>`__.

Creating a branch
~~~~~~~~~~~~~~~~~~

You want your master branch to reflect only production-ready code, so create a
feature branch for making your changes. For example::

    git branch shiny-new-feature
    git checkout shiny-new-feature

The above can be simplified to::

    git checkout -b shiny-new-feature

This changes your working directory to the shiny-new-feature branch.  Keep any
changes in this branch specific to one bug or feature so it is clear
what the branch brings to *traja*. You can have many shiny-new-features
and switch in between them using the git checkout command.

To update this branch, you need to retrieve the changes from the master branch::

    git fetch upstream
    git rebase upstream/master

This will replay your commits on top of the latest traja git master.  If this
leads to merge conflicts, you must resolve these before submitting your pull
request.  If you have uncommitted changes, you will need to ``stash`` them prior
to updating.  This will effectively store your changes and they can be reapplied
after updating.

.. _contributing.dev_env:

2) Creating a development environment
---------------------------------------
A development environment is a virtual space where you can keep an independent installation of *traja*.
This makes it easy to keep both a stable version of python in one place you use for work, and a development
version (which you may break while playing with code) in another.

An easy way to create a *traja* development environment is as follows:

- Install either `Anaconda <http://docs.continuum.io/anaconda/>`_ or
  `miniconda <http://conda.pydata.org/miniconda.html>`_
- Make sure that you have :ref:`cloned the repository <contributing.forking>`
- ``cd`` to the *traja** source directory

Tell conda to create a new environment, named ``traja_dev``, or any other name you would like
for this environment, by running::

      conda create -n traja_dev

For a python 3 environment::

      conda create -n traja_dev python=3.8

This will create the new environment, and not touch any of your existing environments,
nor any existing python installation.

To work in this environment, Windows users should ``activate`` it as follows::

      activate traja_dev

Mac OSX and Linux users should use::

      source activate traja_dev

You will then see a confirmation message to indicate you are in the new development environment.

To view your environments::

      conda info -e

To return to you home root environment::

      deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

At this point you can easily do a *development* install, as detailed in the next sections.

3) Installing Dependencies
--------------------------

To run *traja* in an development environment, you must first install
*traja*'s dependencies. We suggest doing so using the following commands
(executed after your development environment has been activated)::

    conda install -c conda-forge shapely
    pip install -r requirements/dev.txt

This should install all necessary dependencies.

Next activate pre-commit hooks by running::

    pre-commit install

4) Making a development build
-----------------------------

Once dependencies are in place, make an in-place build by navigating to the git
clone of the *traja* repository and running::

    python setup.py develop


5) Making changes and writing tests
-------------------------------------

*traja* is serious about testing and strongly encourages contributors to embrace
`test-driven development (TDD) <http://en.wikipedia.org/wiki/Test-driven_development>`_.
This development process "relies on the repetition of a very short development cycle:
first the developer writes an (initially failing) automated test case that defines a desired
improvement or new function, then produces the minimum amount of code to pass that test."
So, before actually writing any code, you should write your tests.  Often the test can be
taken from the original GitHub issue.  However, it is always worth considering additional
use cases and writing corresponding tests.

Adding tests is one of the most common requests after code is pushed to *traja*.  Therefore,
it is worth getting in the habit of writing tests ahead of time so this is never an issue.

*traja* uses the `pytest testing system
<http://doc.pytest.org/en/latest/>`_ and the convenient
extensions in `numpy.testing
<http://docs.scipy.org/doc/numpy/reference/routines.testing.html>`_.

Writing tests
~~~~~~~~~~~~~

All tests should go into the ``tests`` directory. This folder contains many
current examples of tests, and we suggest looking to these for inspiration.


Running the test suite
~~~~~~~~~~~~~~~~~~~~~~

The tests can then be run directly inside your Git clone (without having to
install *traja*) by typing::

    pytest

6) Updating the Documentation
-----------------------------

*traja* documentation resides in the `doc` folder. Changes to the docs are
make by modifying the appropriate file in the `source` folder within `doc`.
*traja* docs us reStructuredText syntax, `which is explained here <http://www.sphinx-doc.org/en/stable/rest.html#rst-primer>`_
and the docstrings follow the `Numpy Docstring standard <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.

Once you have made your changes, you can build the docs by navigating to the `doc` folder and typing::

    make html

The resulting html pages will be located in `doc/build/html`.


7) Submitting a Pull Request
------------------------------

Once you've made changes and pushed them to your forked repository, you then
submit a pull request to have them integrated into the *traja* code base.

You can find a pull request (or PR) tutorial in the `GitHub's Help Docs <https://help.github.com/articles/using-pull-requests/>`_.
Plotting Grid Cell Flow
=======================

Trajectories can be discretized into grid cells and the average flow from
each grid cell to its neighbor can be plotted with :func:`~traja.plotting.plot_flow`, eg:

.. code-block:: python

    traja.plot_flow(df, kind='stream')

:func:`~traja.plotting.plot_flow` ``kind`` Arguments
----------------------------------------------------

* `surface` - 3D surface plot extending :meth:`mpl_toolkits.mplot3D.Axes3D.plot_surface``
* `contourf` - Filled contour plot extending :meth:`matplotlib.axes.Axes.contourf`
* `quiver` - Quiver plot extending :meth:`matplotlib.axes.Axes.quiver`
* `stream` - Stream plot extending :meth:`matplotlib.axes.Axes.streamplot`

See the :ref:`gallery<sphx_glr_gallery>` for more examples.

3D Surface Plot
---------------

.. autofunction:: traja.plotting.plot_surface

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_001.png
   :alt: 3D plot

Quiver Plot
-----------

.. autofunction:: traja.plotting.plot_quiver

.. code-block:: python

    traja.plot_quiver(df, bins=32)

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_002.png
   :alt: quiver plot

Contour Plot
------------

.. autofunction:: traja.plotting.plot_contour

.. code-block:: python

    traja.plot_contour(df, filled=False, quiver=False, bins=32)

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_003.png
   :alt: contour plot

Contour Plot (Filled)
---------------------

.. code-block:: python

    traja.plot_contour(df, filled=False, quiver=False, bins=32)

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_004.png
   :alt: contour plot filled

Stream Plot
-----------

.. autofunction:: traja.plotting.plot_stream

.. code-block:: python

    traja.plot_contour(df, bins=32, contourfplot_kws={'cmap':'coolwarm'})

.. image:: https://traja.readthedocs.io/en/latest/_images/sphx_glr_plot_average_direction_005.png
   :alt: streamplot
Reference
=============

Accessor Methods
----------------

The following methods are available via :class:`traja.accessor.TrajaAccessor`:

.. automodule:: traja.accessor
    :members:
    :undoc-members:
    :noindex:

Plotting functions
------------------

The following methods are available via :mod:`traja.plotting`:

.. automethod:: traja.plotting.animate

.. automethod:: traja.plotting.bar_plot

.. automethod:: traja.plotting.color_dark

.. automethod:: traja.plotting.fill_ci

.. automethod:: traja.plotting.find_runs

.. automethod:: traja.plotting.plot

.. automethod:: traja.plotting.plot_3d

.. automethod:: traja.plotting.plot_actogram

.. automethod:: traja.plotting.plot_autocorrelation

.. automethod:: traja.plotting.plot_contour

.. automethod:: traja.plotting.plot_clustermap

.. automethod:: traja.plotting.plot_flow

.. automethod:: traja.plotting.plot_quiver

.. automethod:: traja.plotting.plot_stream

.. automethod:: traja.plotting.plot_surface

.. automethod:: traja.plotting.plot_transition_matrix

.. automethod:: traja.plotting.plot_xy

.. automethod:: traja.plotting.polar_bar

.. automethod:: traja.plotting.plot_prediction

.. automethod:: traja.plotting.sans_serif

.. automethod:: traja.plotting.stylize_axes

.. automethod:: traja.plotting.trip_grid


Analysis
--------

The following methods are available via :mod:`traja.trajectory`:

.. automethod:: traja.trajectory.angles

.. automethod:: traja.trajectory.calc_angle

.. automethod:: traja.trajectory.calc_convex_hull

.. automethod:: traja.trajectory.calc_derivatives

.. automethod:: traja.trajectory.calc_displacement

.. automethod:: traja.trajectory.calc_heading

.. automethod:: traja.trajectory.calc_turn_angle

.. automethod:: traja.trajectory.calc_flow_angles

.. automethod:: traja.trajectory.cartesian_to_polar

.. automethod:: traja.trajectory.coords_to_flow

.. automethod:: traja.trajectory.determine_colinearity

.. automethod:: traja.trajectory.distance_between

.. automethod:: traja.trajectory.distance

.. automethod:: traja.trajectory.euclidean

.. automethod:: traja.trajectory.expected_sq_displacement

.. automethod:: traja.trajectory.fill_in_traj

.. automethod:: traja.trajectory.from_xy

.. automethod:: traja.trajectory.generate

.. automethod:: traja.trajectory.get_derivatives

.. automethod:: traja.trajectory.grid_coordinates

.. automethod:: traja.trajectory.inside

.. automethod:: traja.trajectory.length

.. automethod:: traja.trajectory.polar_to_z

.. automethod:: traja.trajectory.rediscretize_points

.. automethod:: traja.trajectory.resample_time

.. automethod:: traja.trajectory.return_angle_to_point

.. automethod:: traja.trajectory.rotate

.. automethod:: traja.trajectory.smooth_sg

.. automethod:: traja.trajectory.speed_intervals

.. automethod:: traja.trajectory.step_lengths

.. automethod:: traja.trajectory.to_shapely

.. automethod:: traja.trajectory.traj_from_coords

.. automethod:: traja.trajectory.transition_matrix

.. automethod:: traja.trajectory.transitions

io functions
------------

The following methods are available via :mod:`traja.parsers`:

.. automethod:: traja.parsers.read_file

.. automethod:: traja.parsers.from_df


TrajaDataFrame
--------------

A ``TrajaDataFrame`` is a tabular data structure that contains ``x``, ``y``, and ``time`` columns.

All pandas ``DataFrame`` methods are also available, although they may
not operate in a meaningful way on the ``x``, ``y``, and ``time`` columns.

Inheritance diagram:

.. inheritance-diagram:: traja.TrajaDataFrame

TrajaCollection
---------------

A ``TrajaCollection`` holds multiple trajectories for analyzing and comparing trajectories.
It has limited accessibility to lower-level methods.

.. autoclass:: traja.frame.TrajaCollection

.. automethod:: traja.frame.TrajaCollection.apply_all

.. automethod:: traja.frame.TrajaCollection.plot


API Pages
---------

.. currentmodule:: traja
.. autosummary::
  :template: autosummary.rst
  :toctree: reference/

  TrajaDataFrame
  TrajaCollection
  read_file
Periodicity
===========

Several methods for analyzing periodicity are included.

Autocorrelation
---------------

Autocorrelation is plotted using :meth:`pandas.plotting.autocorrelation_plot`.

.. autofunction:: traja.plotting.plot_autocorrelation

Periodogram (Power Spectum)
---------------------------

Convenience wrapper for :meth:`scipy.signal.periodogram`.

.. autofunction:: traja.plotting.plot_periodogramResampling Trajectories
=======================

Rediscretize
------------
Rediscretize the trajectory into consistent step lengths with :meth:`~traja.trajectory.rediscretize` where the `R` parameter is
the new step length.

.. note::

    Based on the appendix in Bovet and Benhamou, (1988) and Jim McLean's
    `trajr <https://github.com/JimMcL/trajr>`_ implementation.


Resample time
-------------
:meth:`~traja.trajectory.resample_time` allows resampling trajectories by a ``step_time``.

.. autofunction:: traja.trajectory.resample_time


For example:

.. ipython:: python :okwarning:

    import traja

    # Generate a random walk
    df = traja.generate(n=1000) # Time is in 0.02-second intervals
    df.head()

.. ipython:: python :okwarning:

    resampled = traja.resample_time(df, "50L") # 50 milliseconds
    resampled.head()

    fig = resampled.traja.plot()

.. image:: https://raw.githubusercontent.com/justinshenk/traja/master/docs/images/resampled.png


Ramer–Douglas–Peucker algorithm
-------------------------------

.. note::

    Graciously yanked from Fabian Hirschmann's PyPI package ``rdp``.

:func:`~traja.contrib.rdp` reduces the number of points in a line using the Ramer–Douglas–Peucker algorithm::

    from traja.contrib import rdp

    # Create dataframe of 1000 x, y coordinates
    df = traja.generate(n=1000)

    # Extract xy coordinates
    xy = df.traja.xy

    # Reduce points with epsilon between 0 and 1:
    xy_ = rdp(xy, epsilon=0.8)


    len(xy_)

    Output:
    317

Plotting, we can now see the many fewer points are needed to cover a similar area.::

    df = traja.from_xy(xy_)
    df.traja.plot()

.. image:: https://raw.githubusercontent.com/justinshenk/traja/master/docs/source/_static/after_rdp.png

Trajectory Collections
======================

TrajaCollection
-------------------

When handling multiple trajectories, Traja allows plotting and analysis simultaneously.

Initialize a :func:`~traja.frame.TrajaCollection` with a dictionary or ``DataFrame`` and ``id_col``.

Initializing with Dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The keys of the dictionary can be used to identify types of objects in the scene, eg, "bus", "car", "person"::

    dfs = {"car0":df0, "car1":df1, "bus0: df2, "person0": df3}


Or, arbitrary numbers can be used to initialize

.. autoclass:: traja.frame.TrajaCollection

.. ipython::

    from traja import TrajaCollection

    dfs = {idx: traja.generate(idx, seed=idx) for idx in range(10,13)}
    trjs = TrajaCollection(dfs)

    print(trjs)

Initializing with a DataFrame
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A dataframe containing an id column can be passed directly to :func:`~traja.frame.TrajaCollection`, as long as the ``id_col`` is specified::

    trjs = TrajaCollection(df, id_col="id")

Grouped Operations
------------------

Operations can be applied to each trajectory with :func:`~traja.frame.TrajaCollection.apply_all`.

.. automethod:: traja.frame.TrajaCollection.apply_all

Plottting Multiple Trajectories
-------------------------------

Plotting multiple trajectories can be achieved with :func:`~traja.frame.TrajaCollection.plot`.

.. automethod:: traja.frame.TrajaCollection.plot

Colors can be specified for ids by supplying ``colors`` with a lookup dictionary:

.. ipython::

    colors = ["10":"red",
              "11":"red",
              "12":"red",
              "13":"orange",
              "14":"orange"]

or with a substring lookup:

    colors = ["car":"red",
              "bus":"orange",
              "12":"red",
              "13":"orange",
              "14":"orange"]


.. image:: https://raw.githubusercontent.com/justinshenk/traja/master/docs/images/collection_plot.png




Pandas Indexing and Resampling
==============================

Traja is built on top of pandas :class:`~pandas.DataFrame`, giving access to low-level pandas indexing functions.

This allows indexing, resampling, etc., just as in pandas::

    from traja import generate, plot
    import pandas as pd

    # Generate random walk
    df = generate(n=1000, fps=30)

    # Select every second row
    df[::2]

    Output:
              x         y      time
    0  0.000000  0.000000  0.000000
    2  2.364589  3.553398  0.066667
    4  0.543251  6.347378  0.133333
    6 -3.307575  5.404562  0.200000
    8 -6.697132  3.819403  0.266667

You can also do resampling to select average coordinate every second, for example::

    # Convert 'time' column to timedelta
    df.time = pd.to_timedelta(df.time, unit='s')
    df = df.set_index('time')

    # Resample with average for every second
    resampled = df.resample('S').mean()
    plot(resampled)

.. image:: https://raw.githubusercontent.com/justinshenk/traja/master/docs/source/_static/resampled.png

.. ipython:: python :okwarning:
   :suppress:

   import matplotlib
   import pandas as pd
   from traja.plotting import trip_grid
   orig = matplotlib.rcParams['figure.figsize']
   matplotlib.rcParams['figure.figsize'] = [orig[0] * 1.5, orig[1]]
   import matplotlib.pyplot as plt
   plt.close('all')


Plotting Paths
==============

Making plots of trajectories is easy using the :meth:`~traja.accessor.TrajaAccessor.plot` method.

See the :ref:`gallery<sphx_glr_gallery>` for more examples.

.. automodule:: traja.plotting
    :members: bar_plot, plot, plot_quiver, plot_contour, plot_surface, plot_stream, plot_flow, plot_actogram, polar_bar

Trip Grid
---------

Trip grid can be plotted for :class:`~traja.frame.TrajaDataFrame`s with :func:`~traja.accessor.TrajaAccessor.trip_grid`:

.. ipython:: python :okwarning:

    import traja
    from traja import trip_grid

    df = traja.TrajaDataFrame({'x':range(10),'y':range(10)})
    @savefig trip_grid.png
    hist, image = trip_grid(df);


If only the histogram is need for further computation, use the `hist_only` option:

.. ipython:: python

    hist, _ = trip_grid(df, hist_only=True)
    print(hist[:5])


Highly dense plots be more easily visualized using the `bins` and `log` argument:

.. ipython:: python :okwarning:

    # Generate random walk
    df = traja.generate(1000)

    @savefig trip_grid_log.png
    trip_grid(df, bins=32, log=True);

The plot can also be normalized into a density function with `normalize`:

.. ipython:: python :okwarning:

    @savefig trip_grid_normalized.png
    hist, _ = trip_grid(df, normalize=True);


Animate
-------

.. autofunction:: traja.plotting.animateGenerate Random Walk
====================

Random walks can be generated using :func:`~traja.trajectory.generate`.


.. ipython:: python :okwarning:

    import traja
    
    # Generate random walk
    df = traja.generate(1000)

.. autofunction:: traja.trajectory.generate

.. image:: https://raw.githubusercontent.com/justinshenk/traja/master/docs/source/_static/walk_screenshot.png
Support for Traja
=================

Bugs
----

Bugs, issues and improvement requests can be logged in `Github Issues <https://github.com/traja-team/traja/issues>`_. 

Community
---------

Community support is provided via `Gitter <https://gitter.im/traja-chat/community>`_. Just ask a question there.
Reading and Writing Files
=========================

Reading trajectory data
-----------------------

traja allows reading files via :func:`traja.parsers.read_file`. For example a CSV file ``trajectory.csv`` with the
following contents::


    x,y
    1,1
    1,2
    1,3

Could be read in like:

.. code-block:: python

    import traja

    df = traja.read_file('trajectory.csv')

``read_file`` returns a `TrajaDataFrame` with access to all pandas and traja methods.

.. automodule:: traja.accessor
   .. automethod::

Any keyword arguments passed to `read_file` will be passed to :meth:`pandas.read_csv`.

Data frames can also be read with pandas :func:`pandas.read_csv` and then converted to TrajaDataFrames
with:

.. code-block:: python

    import traja
    import pandas as pd

    df = pd.read_csv('data.csv')

    # If x and y columns are named different than "x" and "y", rename them, eg:
    df = df.rename(columns={"x_col": "x", "y_col": "y"}) # original column names x_col, y_col
    
    # If the time column doesn't include "time" in the name, similarly rename it to "time"

    trj = traja.TrajaDataFrame(df)



Writing trajectory data
-----------------------

Files can be saved using the built-in pandas :func:`pandas.to_csv`.

.. code-block:: python

    df.to_csv('trajectory.csv')
Turns and Angular Analysis
==========================

Turns
-----

Turns can be calculated using :func:`~traja.trajectory.calc_angle`.

.. autofunction:: traja.trajectory.calc_angle

Heading
-------

Heading can be calculated using :func:`~traja.trajectory.calc_heading`.

.. autofunction:: traja.trajectory.calc_heading

Angles
------

Angles can be calculated using :func:`~traja.trajectory.angles`.

.. autofunction:: traja.trajectory.angles

Smoothing and Analysis
======================

Smoothing
---------

Smoothing can be performed using :func:`~traja.trajectory.smooth_sg`.

.. autofunction:: traja.trajectory.smooth_sg

.. ipython::

    df = traja.generate()
    smoothed = traja.smooth_sg(df, w=101)
    smoothed.traja.plot()

.. image:: https://raw.githubusercontent.com/justinshenk/traja/master/docs/images/smoothed.png


Length
------

Length of trajectory can be calculated using :func:`~traja.trajectory.length`.

.. autofunction:: traja.trajectory.length

Distance
--------

Net displacement of trajectory (start to end) can be calculated using :func:`~traja.trajectory.distance`.

.. autofunction:: traja.trajectory.distance

Displacement
------------

Displacement (distance travelled) can be calculated using :func:`~traja.trajectory.calc_displacement`.

.. autofunction:: traja.trajectory.calc_displacement

Derivatives
-----------

.. autofunction:: traja.trajectory.get_derivatives

Speed Intervals
---------------

.. autofunction:: traja.trajectory.speed_intervals.. traja documentation master file, created by
   sphinx-quickstart on Mon Jan 28 23:36:32 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

traja |version|
===============

Trajectory Analysis in Python

Traja allows analyzing trajectory datasets using a wide range of tools, including pandas and R.
Traja extends the capability of pandas :class:`~pandas.DataFrame` specific for animal or object trajectory
analysis in 2D, and provides convenient interfaces to other geometric analysis packages (eg, shapely).

Description
-----------

The Traja Python package is a toolkit for the numerical characterization and analysis
of the trajectories of moving animals. Trajectory analysis is applicable in fields as
diverse as optimal foraging theory, migration, and behavioural mimicry
(e.g. for verifying similarities in locomotion).
A trajectory is simply a record of the path followed by a moving object.
Traja operates on trajectories in the form of a series of locations (as x, y coordinates) with times.
Trajectories may be obtained by any method which provides this information,
including manual tracking, radio telemetry, GPS tracking, and motion tracking from videos.

The goal of this package (and this document) is to aid biological researchers, who may not have extensive
experience with Python, to analyse trajectories
without being restricted by a limited knowledge of Python or programming.
However, a basic understanding of Python is useful.

If you use Traja in your publications, please cite our paper_ in Journal of Open Source
Software:

.. code-block:: txt

   @article{Shenk2021,
   doi = {10.21105/joss.03202},
   url = {https://doi.org/10.21105/joss.03202},
   year = {2021},
   publisher = {The Open Journal},
   volume = {6},
   number = {63},
   pages = {3202},
   author = {Justin Shenk and Wolf Byttner and Saranraj Nambusubramaniyan and Alexander Zoeller},
   title = {Traja: A Python toolbox for animal trajectory analysis},
   journal = {Journal of Open Source Software}
   }

.. _paper: https://joss.theoj.org/papers/10.21105/joss.03202

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   Installation <install>
   Examples Gallery <gallery/index>

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   Reading and Writing Files <reading>
   Pandas Indexing and Resampling <pandas>
   Generate Random Walk <generate>
   Smoothing and Analysis <calculations>
   Turns <turns>
   Plotting Paths <plots>
   Periodicity <periodicity>
   Plotting Grid Cell Flow <grid_cell>
   Rediscretizing Trajectories <rediscretize>
   Clustering and Dimensionality Reduction <clustering>
   Collections / Scenes <collections>
   Predicting Trajectories <predictions>

.. toctree::
   :maxdepth: 1
   :caption: Reference Guide

   Reference to All Attributes and Methods <reference>
   Bugs and Support <support>

.. toctree::
   :maxdepth: 1
   :caption: Developer

   Contributing to Traja <contributing>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Predicting Trajectories
=======================

Predicting trajectories with traja can be done with a recurrent neural network (RNN). Traja includes
the Long Short Term Memory (LSTM), LSTM Autoencoder (LSTM AE) and LSTM Variational Autoencoder (LSTM VAE)
RNNs. Traja also supports custom RNNs.

To model a trajectory using RNNs, one needs to fit the network to the model. Traja includes the MultiTaskRNNTrainer
that can solve a prediction, classification and regression problem with traja DataFrames.

`Traja` also includes a DataLoader that handles traja dataframes.

Below is an example with a prediction LSTM via :py:class:`traja.models.predictive_models.lstm.LSTM`.

.. code-block:: python

    import traja

    df = traja.dataset.example.jaguar()

.. note::
    LSTMs work better with data between -1 and 1. Therefore the data loader
    scales the data.

.. code-block:: python

    batch_size = 10 # How many sequences to train every step. Constrained by GPU memory.
    num_past = 10 # How many time steps from which to learn the time series
    num_future = 10 # How many time steps to predict
    split_by_id = False # Whether to split data into training, test and validation sets based on
                        # the animal's ID or not. If True, an animal's entire trajectory will only
                        # be used for training, or only for testing and so on.
                        # If your animals are territorial (like Jaguars) and you want to forecast
                        # their trajectories, you want this to be false. If, however, you want to
                        # classify the group membership of an animal, you want this to be true,
                        # so that you can verify that previously unseen animals get assigned to
                        # the correct class.


.. autoclass:: traja.models.predictive_models.lstm.LSTM
    :members:

    dataloaders = traja.dataset.MultiModalDataLoader(df,
                                                batch_size=batch_size,
                                                n_past=num_past,
                                                n_future=num_future,
                                                num_workers=0,
                                                split_by_id=split_by_id)

.. note::

    The width of the hidden layers and depth of the network are the two main ways in which
    one tunes the performance of the network. More complex datasets require wider and deeper
    networks. Below are sensible defaults.

.. code-block:: python

    from traja.models.predictive_models.lstm import LSTM
    input_size = 2 # Number of input dimensions (normally x, y)
    output_size = 2 # Same as input_size when predicting
    num_layers = 2 # Number of LSTM layers. Deeper learns more complex patterns but overfits.
    hidden_size = 32 # Width of layers. Wider learns bigger patterns but overfits. Try 32, 64, 128, 256, 512
    dropout = 0.1 # Ignore some network connections. Improves generalisation.

    model = LSTM(input_size=input_size,
                 hidden_size=hidden_size,
                 num_layers=num_layers,
                 output_size=output_size,
                 dropout=dropout,
                 batch_size=batch_size)

.. code-block:: python

    from traja.models.train import HybridTrainer

    optimizer_type = 'Adam' # Nonlinear optimiser with momentum
    loss_type = 'huber'

    # Trainer
    trainer = HybridTrainer(model=model,
                            optimizer_type=optimizer_type,
                            loss_type=loss_type)
    # Train the model
    trainer.fit(dataloaders, model_save_path='./model.pt', epochs=10, training_mode='forecasting')

After training, you can determine the network's final performance with test data, if you want to pick
the best model, or with validation data, if you want to determine the performance of your model.

The ``dataloaders`` dictionary contains the ``sequential_test_loader`` and ``sequential_validation_loader``,
that preserve the order of the original data. The dictionary also contains the 'test_loader' and
``validation_loader`` data loaders, where the order of the time series is randomised.

.. code-block:: python

    validation_loader = dataloaders['sequential_validation_loader']

    trainer.validate(validation_loader)

Finally, you can display your training results using the built-in plotting libraries.

.. code-block:: python

    from traja.plotting import plot_prediction

    batch_index = 0  # The batch you want to plot
    plot_prediction(model, validation_loader, batch_index)

.. image:: _static/rnn_prediction.png

Parameter searching
-------------------

When optimising neural networks, you often want to change the parameters. When training a forecaster,
you have to reinitialise and retrain your model. However, when training a classifier or regressor, you
can reset these on the fly, since they work directly on the latent space of your model.
VAE models provide utility functions to make this easy.

.. code-block:: python

    from traja.models import MultiModelVAE
    input_size = 2 # Number of input dimensions (normally x, y)
    output_size = 2 # Same as input_size when predicting
    num_layers = 2 # Number of LSTM layers. Deeper learns more complex patterns but overfits.
    hidden_size = 32 # Width of layers. Wider learns bigger patterns but overfits. Try 32, 64, 128, 256, 512
    dropout = 0.1 # Ignore some network connections. Improves generalisation.

    # Classifier parameters
    classifier_hidden_size = 32
    num_classifier_layers = 4
    num_classes = 42

    # Regressor parameters
    regressor_hidden_size = 18
    num_regressor_layers = 1
    num_regressor_parameters = 3

    model = MultiModelVAE(input_size=input_size,
                          hidden_size=hidden_size,
                          num_layers=num_layers,
                          output_size=output_size,
                          dropout=dropout,
                          batch_size=batch_size,
                          num_future=num_future,
                          classifier_hidden_size=classifier_hidden_size,
                          num_classifier_layers=num_classifier_layers,
                          num_classes=num_classes,
                          regressor_hidden_size=regressor_hidden_size,
                          num_regressor_layers=num_regressor_layers,
                          num_regressor_parameters=num_regressor_parameters)

    new_classifier_hidden_size = 64
    new_num_classifier_layers = 2

    model.reset_classifier(classifier_hidden_size=new_classifier_hidden_size,
                           num_classifier_layers=new_num_classifier_layers)

    new_regressor_hidden_size = 64
    new_num_regressor_layers = 2
    model.reset_regressor(regressor_hidden_size=new_regressor_hidden_size,
                          num_regressor_layers=new_num_regressor_layers)traja.contrib package
=====================

Submodules
----------


Module contents
---------------

.. automodule:: traja.contrib
    :members:
    :undoc-members:
    :show-inheritance:
Installation
============

Installing traja
----------------

traja requires Python 3.6+ to be installed. For installing on Windows,
it is recommend to download and install via conda_.

To install via conda::

    conda install -c conda-forge traja

To install via pip::

   pip install traja

To install the latest development version, clone the `GitHub` repository and use the setup script::

   git clone https://github.com/traja-team/traja.git
   cd traja
   pip install .

Dependencies
------------

Installation with pip should also include all dependencies, but a complete list is

- numpy_
- matplotlib_
- scipy_
- pandas_

To install all optional dependencies run::

  pip install 'traja[all]'


.. _GitHub: https://github.com/justinshenk/github

.. _numpy: http://www.numpy.org

.. _pandas: http://pandas.pydata.org

.. _scipy: https://docs.scipy.org/doc/scipy/reference/

.. _shapely: http://toblerity.github.io/shapely

.. _matplotlib: http://matplotlib.org

.. _conda: https://docs.conda.io/en/latest/{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

{% if objtype in ['class', 'method', 'function'] %}

{% endif %}
