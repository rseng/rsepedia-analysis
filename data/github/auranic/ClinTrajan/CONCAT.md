# Methodology and software for quantifying pseudotemporal trajectories in clinical datasets

## Clinical datasets

Large clinical datasets are becoming increasingly available for systematic mining the associations between phenotypic variables,
characterizing a particular disease. Specific properties of the clinical datasets represent challenges for applying the 
standard repertoire of machine learning methods. 

We assume that the real-life clinical data are characterized by the following features:
1) They contain mixed data types (continuous, binary, ordinal, categorical variables, censored data)
2) They typically contain missing values with non-uniform pattern across the data matrix
3) They do not have a uniquely defined labeling (part of the clinical variables can be used to define clinical groups, 
but this can be done in several meaningfull ways)

## Clinical trajectories, pseudotime, and dynamical phenotyping of synchronic data

We assume that any clinical dataset represents a landscape of possible patient phenotypes of various and multivariate 
degree of clinical gravity, which can be accompanied by the details of applied treatment. 
We also assume a possibility of existence of clinical trajectories, i.e. clinically relevant sequences of partly ordered phenotypes 
possibly representing consecutive states of a developing disease phenotype and leading to some final states (i.e., a 
lethal outcome). Each clinical trajectory can be characterized by its proper pseudotime which allows to quantitatively characterize the degree of progression along the trajectory. 
Each clinical variable can be plotted then as a function of pseudotime for a given clinical trajectory.
We also assume that clinical trajectories can be characterized by branching structure, representing
some important bifurcations in the development of a disease. 

Extracting cellular trajectories is a widely used methodology of data analysis in genomics, especially in studying certain highly dynamic phenomena such as differentiation or development. 
Quantifying and visualizing clinical trajectories represents a more challenging data mining problem due to the data specificity.

## Method of elastic principal graphs (ElPiGraph) for extracting and analyzing the clinical trajectories

Here we develop a semi-supervised methodology of clinical data analysis, based on constructing and exploring the properties
of principal trees (PT), which is a non-linear generalization of Principal Component Analysis (PCA). Principal trees are 
constructed using ElPiGraph method, which has been previously exploited in determining branching trajectories in various genoomics 
datasets (in particular, in single cell omics data). 

The methodology takes into account the specificity of clinical data by providing tools for the following steps of clinical data analysis:

1) Univariate and multi-variate quantification of nominal variables, including the original implementation of optimal scaling for ordinal variables
2) Several original methods for missing values imputation based on Singular Value Decomposition
3) Set of state-of-the-art methods for manifold learning
4) Partitioning the data accordingly to the branches of the principal tree (analogue of clustering) and associating the branches to clinical variables.
5) Extracting clinical trajectories using principal tree approach and associating the trajectories to clinical variables.
6) Visualization of clinical variables using principal trees
7) Pseudotime plots of clinical variables along clinical trajectories

The methodology is implemented in Python.

## Installation

For the moment the only way to use the package is to copy the .py files from the 'code' folder and make them available in the Python path

## Tutorial

A simplified tutorial showing the basic steps of ClinTrajan analysis is [available here](https://github.com/auranic/ClinTrajan/blob/master/tutorial/tutorial.md).
This tutorial is just better explained presentation of [this Jupyter notebook](ClinTrajan_tutorial.ipynb).

## Case studies

We demonstrate application of the methodology to two clinical datasets, one of moderate size (1700 patients) and one of relatively large size (100000 patients).

### Complications of myocardial infarction

The database was collected in the Krasnoyarsk Interdistrict Clinical Hospital (Russia) in 1992-1995 years. The original database and its description can be downloaded from https://leicester.figshare.com/articles/Myocardial_infarction_complications_Database/12045261/1. It contains information about 1700 patients and 110 features characterizing the clinical phenotypes and 12 features representing possible complications of the myocardial infarction disease. 

Two Jupyter notebooks provide the exact protocol of the analysis of this database.
In order to use them, download the content of the git and start the notebook from the git folder.

* [QI_Infarctus.ipynb](QI_Infarctus.ipynb) - notebook documenting quantification and imputation of the datatable, which consists of the steps
  1. Removing the columns containing more than 30% of missing values
  2. Removing the rows containing more than 20% of missing values
  3. Determining the complete part of the table
  4. Classifying variables into types (BINARY, ORDINAL, CONTINUOUS). The categorical variables are supposed to be converted using the standard dummy coding.
  5. Univariate variable quantification
  6. Using the quantified complete part of the table, compute SVD of an order corresponding to the intrinsic dimension estimate
  7. Project vectors with missing values into the space of obtained principal components, the imputed values are taken from the projection values.

* [PT_Infarctus.ipynb](PT_Infarctus.ipynb) - notebook documenting the analysis of the imputed table, using the methodology of principal trees. It contains of the following steps:
  1. Pre-processing the dataset by projecting it into the space of the first principal components. Conservative 'elbow rule'-based estimate for the number of top principal components is used in this case.
  2. [Defining the classes of patients](images/definition_of_classes_infarctus.png), by a separate analysis of dependent (complication and lethality cause) variables, using principal trees. Note: when viewed online, the notebook misses the image showing the introduced classification of the patients. This image can be found [here](images/definition_of_classes_infarctus.png).
  3. Constructing the principal tree, and post-processing it (pruning short edges and extending the external branches in order to avoid the border effects).
  4. Computing and visualizing associations of the principal tree branches with the patient classes.
  5. Determining the 'root node' of the tree, most associated to the 'no complication' class.
  6. Compute and visualize all trajectories from the root node to the leaf nodes.
  7. Compute and visualize associations of the trajectories with all variables
  8. Compute and visualize the connection of complication variables with trajectories.
  9. Using principal tree for visualization of various variables
  10. Applying a panel of 12 manifold learning methods to the dataset


### Diabetes readmission data set from UCI Machine Learning Repository

The dataset represents 10 years (1999-2008) of clinical care at 130 US hospitals and integrated delivery networks. It includes over 50 features representing patient and hospital outcomes. The dataset can be downloaded from UCI repository at https://archive.ics.uci.edu/ml/datasets/diabetes+130-us+hospitals+for+years+1999-2008 or from Kaggle at https://www.kaggle.com/brandao/diabetes. The data matrix contains 100000 hospitalization cases with patients suffering from diabetis characterized by 55 attributes.

Three Jupyter notebooks provide the exact protocol of the analysis of this database.
In order to use them, download the content of the git and start the notebook from the git folder.

* [Set_up_Diabetes_Data_for_Modeling.ipynb](Set_up_Diabetes_Data_for_Modeling.ipynb) - notebook providing the encoding of the diabetes readmission dataset
* [QI_Diabetes.ipynb](QI_Diabetes.ipynb) - notebook providing the details of quantification of the diabetes dataset
* [PT_Diabetes.ipynb](PT_Diabetes.ipynb) - notebook providing the details of the analysis of clinical trajectories for the diabetes dataset


## References:

(1) [Albergante L, Mirkes E, Bac J, Chen H, Martin A, Faure L, Barillot E, Pinello L, Gorban A, Zinovyev A. Robust and scalable learning of complex intrinsic dataset geometry via ElPiGraph. 2020. Entropy 22(3):296](https://www.mdpi.com/1099-4300/22/3/296)

(2) [Golovenkin SE, Bac J, Chervov A, Mirkes EM, Orlova YV, Barillot E, Gorban AN, Zinovyev A. Trajectories, bifurcations and pseudotime in large clinical datasets: applications to myocardial infarction and diabetes data. 2020. arXiv:2007.03788](https://arxiv.org/abs/2007.03788)
# ClinTrajan tutorial 

![Trajan emperor, optimus princeps](https://github.com/auranic/ClinTrajan/blob/master/images/trajan.png)

Application of ClinTrajan to a clinical dataset consists in two parts:

1. [Quantification of the data](##quantification-of-the-data)
2. [Application of ElPiGraph](##application-of-elpigraph-to-quantify-branching-pseudotime) and  [downstream analyses](###downstream-analyses-using-elpigraph)

Here we illustrate only the most basic analysis steps using the [dataset of myocardial infarction complications](https://leicester.figshare.com/articles/dataset/Myocardial_infarction_complications_Database/12045261/3).  In order to follow the tutorial, one has to download the [ClinTrajan git](https://github.com/auranic/ClinTrajan) and unpack locally. The easiest way to run the tutorial is to run the code through [this ClinTrajan tutorial Jupyter notebook](../ClinTrajan_tutorial.ipynb). Alternatively, one can copy-paste and run the commands in any convenient Python environment. 

There exist also [complete Jupyter notebooks](https://github.com/auranic/ClinTrajan/), allowing one to reproduce all the analysis steps reported in the [ClinTrajan manuscript](https://arxiv.org/abs/2007.03788).


## Quantification of the data

The quantification functions of ClinTrajan are stored in the module **clintraj_qi**, so first of all we will import it. In addition, we import other modules of ClinTrajan and some other standard functions.

```
import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from clintraj_qi import *
from clintraj_optiscale import *
from clintraj_eltree import *
from clintraj_util import *
from clintraj_ml import *
```

Now, let us import a dataset. It is assumed that the table is tab-delimited, contains only numbers or missing values in any row or column except for the first column (containing observation name) and the first row (containing variable names). Here it is assumed that certain column and rows of the table, containing too many missing values, have been eliminated. If the original dataset contains categorical nominal variables, they must be encoded first, using, for example, dummy encoding:

```
df = pd.read_csv('data/infarction/all_dummies.txt',delimiter='\t')
display(df)
quantify_nans(df)
```

![](https://github.com/auranic/ClinTrajan/blob/master/images/table_import.png)

As one can see, only 34% of rows have complete values for all columns.

We need to know what types different variables have, which can be 'BINARY', 'ORDINAL', 'CONTINUOUS:
```
variable_types, binary, continuous, ordinal = detect_variable_type(df,10,verbose=False)
```
We need to impute missing values in the table, but for this we need some quantification already. Simple univariate quantification can be done via:
```
dfq,replacement_info = quantify_dataframe_univariate(df,variable_types)
with open('temp.txt','w') as fid:
    fid.write(replacement_info)
```
Note that we have written down just in case the quantification parameters such that we could use them after imputation of missing values and restore the data table to the initial variable scales.

Now we can impute the missing values. One of simple ideas is to compute SVD on the complete part of the data matrix and then project the data points with missing variables onto the principal components. The imputed value will be the value of the variable in the projection point.
```
dfq_imputed = SVDcomplete_imputation_method(dfq, variable_types, verbose=True,num_components=-1)
dequant_info = invert_quant_info(load_quantification_info('temp.txt'))
df_imputed = dequantify_table(dfq_imputed,dequant_info)
display(df_imputed)
```

Now, we are ready to quantify the data table. We will do it by applying optimal scaling to the ordinal values. 

```
df = remove_constant_columns_from_dataframe(df_imputed)
variable_names = [str(s) for s in df.columns[1:]]
X = df[df.columns[1:]].to_numpy()
X_original = X
X_before_scaling = X.copy()
X,cik = optimal_scaling(X,variable_types,verbose=True,vmax=0.6)
```
The output looks like this:

![](https://github.com/auranic/ClinTrajan/blob/master/images/optimal_scaling.png)

which means that the sum of squared correlations (Q2 value) between all quantified ordinal variables and between ordinal and numerical variables increased, and in the final correlation table one can see more significant correlations.

Congrats, now we have the data matrix *X* to which we can apply the [ElPiGraph](https://sysbio-curie.github.io/elpigraph/) algorithm! Note that we also kept the original matrix *X_original* and the list with the list of variables *variable_names* and with variable types *variable_types*.

## Application of [ElPiGraph](https://sysbio-curie.github.io/elpigraph/) to quantify branching pseudotime

Before applying ElPiGraph, let us first reduce the dimensionality of the dataset (more than 100!). [Some tests](https://github.com/j-bac/scikit-dimension) showed that the intrinsic dimensionality is much lower (around 12!), so let us apply Principal Component Analysis in order to reduce the dimension to this number. Note that we center all variables to zero mean and scale them to unit variance. Setting *svd_solver* to 'full' helps getting reproducible results (yes, by default, PCA is not reproducible in sckit-learn!)

```
reduced_dimension = 12
X = scipy.stats.zscore(X)
pca = PCA(n_components=X.shape[1],svd_solver='full')
Y = pca.fit_transform(X)
v = pca.components_.T
mean_val = np.mean(X,axis=0)
X = Y[:,0:reduced_dimension]
```
Now we construct the [principal tree](https://www.mdpi.com/1099-4300/22/3/296). You can specify only one parameter for the number of nodes in the tree, but here we explicitly show the values of the other parameters which are close to the default ones.

```
nnodes = 50
tree_elpi = elpigraph.computeElasticPrincipalTree(X,nnodes,drawPCAView=True,
                                                  alpha=0.01,Mu=0.1,Lambda=0.05,
                                                  FinalEnergy='Penalized')
tree_elpi = tree_elpi[0]
# some additional pruning of the graph
prune_the_tree(tree_elpi)
# extend the leafs to reach the extreme data points
tree_extended = ExtendLeaves_modified(X, tree_elpi, Mode = "QuantDists", ControlPar = .5, DoSA = False)
```
This produces a simple plot with a projection of the principal tree on PCA plane. In 12D, the topology of the tree is more complex then it seems from 2D PCA projection!

![](https://github.com/auranic/ClinTrajan/blob/master/images/principal_tree.png)

In particular, in the linear 2D projection, we can see more or less clearly only two branches of the principal tree.

Before moving any further, we will need two partitionings of the data points, by proximity to the node of the graph in the multi-dimensional data space, and by the tree 'segment' (meaning a sequence of nodes in the tree without any branching point).

```
# paritioning the data by tree branches
vec_labels_by_branches = partition_data_by_tree_branches(X,tree_extended)
# paritioning the data by proximity to nodes
partition, dists = elpigraph.src.core.PartitionData(X = X, NodePositions = tree_elpi['NodePositions'], 
                                                    SquaredX = np.sum(X**2,axis=1,keepdims=1),
                                                    MaxBlockSize = 100000000, TrimmingRadius = np.inf
                                                    )
partition_by_node = np.zeros(len(partition))
for i,p in enumerate(partition):
    partition_by_node[i] = p[0]
```

In order to visualize the intrinsic geometry of the principal tree, and use it to project the data points from R<sup>N</sup> to R<sup>2</sup>, we can apply a version of [force-directed layout](https://en.wikipedia.org/wiki/Force-directed_graph_drawing) to the principal tree (remember that a tree is a [planar graph](https://en.wikipedia.org/wiki/Planar_graph)!)

Let us visualize the tree with data points colored by the proximity to the tree segments:

```
fig = plt.figure(figsize=(8, 8))
visualize_eltree_with_data(tree_extended,X,X_original,v,mean_val,'k',variable_names,
                          Color_by_partitioning = True, visualize_partition = vec_labels_by_branches)
plt.show()
```
![](https://github.com/auranic/ClinTrajan/blob/master/images/principal_tree_segments.png)

### Downstream analyses using ElPiGraph

Now let us visualize something more interesting using the principal tree. We will visualize all lethal cases of myocardial infarction complications, and by the thickness of the tree edges, we will visualize the mortality trend along various clinical trajectories. Note that in our table the variable *LET_IS_0* means *'no lethal outcome'*!

```
fig = plt.figure(figsize=(8, 8))
non_lethal_feature = 'LET_IS_0'
visualize_eltree_with_data(tree_extended,X,X_original,v,mean_val,'k',variable_names,
                          Color_by_feature=non_lethal_feature, Feature_Edge_Width=non_lethal_feature,
                           Invert_Edge_Value=True,Min_Edge_Width=10,Max_Edge_Width=50,
                           Visualize_Edge_Width_AsNodeCoordinates=True,cmap='winter')
plt.show()
```
![](https://github.com/auranic/ClinTrajan/blob/master/images/principal_tree_lethality.png )

Ok, let us do some more insightfull visualizations. Not let us highlight all patients with age <65 years having bronchyal asthma in anamnesis:

```
fig = plt.figure(figsize=(8, 8))
inds = np.where((X_original[:,variable_names.index('AGE')]<=65)&(X_original[:,variable_names.index('zab_leg_03')]==1))[0]
colors = ['k' for i in range(len(X))]
for i in inds:
    colors[i] = 'r'
visualize_eltree_with_data(tree_extended,X,X_original,v,mean_val,colors,variable_names,
                          highlight_subset=inds,Big_Point_Size=100,cmap='hot')
plt.show()
```

![](https://github.com/auranic/ClinTrajan/blob/master/images/principal_tree_asthma.png)

Further we want to quantify the pseudotime, but for this we need to define the root node. Here it should be the node corresponding to the least complicated case of myocardial infarction. Let us make yet another visualization in order to find it. For this we will visualize, using pie-charts, the proportion of complications in each node of the tree. The pie-chart will show the fraction of uncomplicated cases by black, and the fraction of complications by red. The size of the chart corresponds to the number of data points it represents (for which this is the closest node of the graph). Note that the lethal outcome is considered as - serious! - complication.

```
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1,1,1)
complication_vars = ['FIBR_PREDS','PREDS_TAH','JELUD_TAH','FIBR_JELUD',
                     'A_V_BLOK','OTEK_LANC','RAZRIV','DRESSLER',
                     'ZSN','REC_IM','P_IM_STEN']
inds_compl = [variable_names.index(a) for a in complication_vars]
lethal = 1-X_original[:,variable_names.index('LET_IS_0')]
has_complication = np.sum(X_original[:,inds_compl],axis=1)>0
inds = np.where((has_complication==0)&(lethal==0))[0]
colors = ['r' for i in range(len(X))]
for i in inds:
    colors[i] = 'k'
visualize_eltree_with_data(tree_extended,X,X_original,v,mean_val,colors,variable_names,
                          highlight_subset=inds,Big_Point_Size=2,Normal_Point_Size=2,showNodeNumbers=True)
add_pie_charts(ax,tree_extended['NodePositions2D'],colors,['r','k'],partition,scale=30)
plt.show()
root_node = 8
print('Root node=',root_node)
```

![](https://github.com/auranic/ClinTrajan/blob/master/images/principal_tree_complications.png)

As a result, we have just identified the node number 8 as the potential root because the proportion of complications in it was the least. Now we know from where to start the clinical trajectories and quantify pseudotime:

```
all_trajectories,all_trajectories_edges = extract_trajectories(tree_extended,root_node)
print(len(all_trajectories),' trajectories found.')
ProjStruct = project_on_tree(X,tree_extended)
PseudoTimeTraj = quantify_pseudotime(all_trajectories,all_trajectories_edges,ProjStruct)
```

Nine clinical trajectories have been finally identified. As in clustering, it is up to us to give some meaning for them though. One of the ways to do it, is to test if a set of clinical variables can be associated to some trajectories via non-linear regression which - in the case of binary variables - can be the logistic regression:

```
vars = ['ritm_ecg_p_01','ritm_ecg_p_02','ritm_ecg_p_04']
for var in vars:
    List_of_Associations = regression_of_variable_with_trajectories(PseudoTimeTraj,var,variable_names,
                                                                    variable_types,X_original,R2_Threshold=0.5,
                                                                    producePlot=True,
                                                                    Continuous_Regression_Type='gpr',
                                                                    verbose=True)
```
Here some of the detected associations are visualized. The actual data values are shown as red points, the probability given by logistic regression is shown as red line, and a simple sliding window smoothing of the data is shown as blue line.

![](https://github.com/auranic/ClinTrajan/blob/master/images/variable_associations.png)

Note that the scale (number of nodes) along each trajectory can be different. In general, the pseudotime quantified along different trajectories can not be directly compared because it can correspond to completely different physical time scale!

The results of regression application for several variables can be shown simultaneously at the same plot:
```
pstt = PseudoTimeTraj[1]
colors = ['r','b','g']
for i,var in enumerate(vars):
    vals = draw_pseudotime_dependence(pstt,var,variable_names,variable_types,X_original,colors[i],
                                               linewidth=3,draw_datapoints=False)
plt.legend()
plt.show()
```

![](https://github.com/auranic/ClinTrajan/blob/master/images/variable_associations_several.png)

Finally, the pseudotime can be used to visualize pseudo-temporal profiles of any quantity that can be derived from the data. For example, we can compute the cumulative function of hazard in order to visualize the lethality risks along different clinical trajectories. We will do this, by using the survival analysis library 'lifelines'.

```
import lifelines
from lifelines import SplineFitter
from lifelines import NelsonAalenFitter
from lifelines import KaplanMeierFitter
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','tab:pink','tab:green']

event_data = np.zeros((len(df),2))
events = 1-np.array(df['LET_IS_0'])
label = 'Death'

for i,pstt in enumerate(PseudoTimeTraj):
    points = pstt['Points']
    times = pstt['Pseudotime']
    for i,p in enumerate(points):
        event_data[p,0] = times[i]
        event_data[p,1] = events[p]

plt.figure(figsize=(8,8))

for i,pstt in enumerate(PseudoTimeTraj):
    TrajName = 'Trajectory:'+str(pstt['Trajectory'][0])+'--'+str(pstt['Trajectory'][-1])
    points = pstt['Points']
    naf = NelsonAalenFitter()
    T = event_data[points,0]
    E = event_data[points,1]
    naf.fit(event_data[points,0], event_observed=event_data[points,1],label=TrajName)  
    naf.plot_hazard(bandwidth=3.0,fontsize=20,linewidth=10,color=colors[i])
```

![](https://github.com/auranic/ClinTrajan/blob/master/images/hazard_pseudotime.png)

We can see that six out of nine trajectories are associated with significantly increasing lethality risk.
