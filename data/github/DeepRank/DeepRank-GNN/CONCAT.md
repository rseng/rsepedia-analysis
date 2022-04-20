# DeepRank-GNN


[![Build Status](https://github.com/DeepRank/DeepRank-GNN/workflows/build/badge.svg)](https://github.com/DeepRank/DeepRank-GNN/actions)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/f3f98b2d1883493ead50e3acaa23f2cc)](https://app.codacy.com/gh/DeepRank/DeepRank-GNN?utm_source=github.com&utm_medium=referral&utm_content=DeepRank/DeepRank-GNN&utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/DeepRank/Deeprank-GNN/badge.svg?branch=master)](https://coveralls.io/github/DeepRank/Deeprank-GNN?branch=master)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5705564.svg)](https://doi.org/10.5281/zenodo.5705564)

![alt-text](./deeprank_gnn.png)

## Installation

Before installing DeepRank-GNN you need to install pytorch_geometric according to your needs. You can find detailled instructions here :
  * pytorch_geometric : https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html

By default the CPU version of pytorch will be installed but you can also customize that installation following the instructions at:
  * pytorch : https://pytorch.org/ 

Once the dependencies installed, you can install the latest release of DeepRank-GNN using the PyPi package manager:

```
pip install DeepRank-GNN
```

Alternatively you can get all the new developments by cloning the repo and installing the code with

```
git clone https://github.com/DeepRank/Deeprank-GNN 
cd DeepRank-GNN
pip install -e ./
```

The documentation can be found here : https://deeprank-gnn.readthedocs.io/ 

## Generate Graphs

All the graphs/line graphs of all the pdb/pssm stored in `data/pdb/` and `data/pssm/` with the `GenGraph.py` script. This will generate the hdf5 file `graph_residue.hdf5` which contains the graph of the different conformations.


```python
from GraphGenMP import GraphHDF5

pdb_path = './data/pdb'
pssm_path = './data/pssm'
ref = './data/ref'

GraphHDF5(pdb_path=pdb_path,ref_path=ref,pssm_path=pssm_path,
	      graph_type='residue',outfile='graph_residue.hdf5')
```

## Graph Interaction Network

Using the graph interaction network is rather simple :


```python
from deeprank_gnn.NeuralNet import NeuralNet
from deeprank_gnn.ginet import GINet

database = './hdf5/1ACB_residue.hdf5'

NN = NeuralNet(database, GINet,
               node_feature=['type', 'polarity', 'bsa',
                             'depth', 'hse', 'ic', 'pssm'],
               edge_feature=['dist'],
               target='irmsd',
               index=range(400),
               batch_size=64,
               percent=[0.8, 0.2])

NN.train(nepoch=250, validate=False)
NN.plot_scatter()
```

## Custom GNN

It is also possible to define new network architecture and to specify the loss and optimizer to be used during the training.

```python


def normalized_cut_2d(edge_index, pos):
    row, col = edge_index
    edge_attr = torch.norm(pos[row] - pos[col], p=2, dim=1)
    return normalized_cut(edge_index, edge_attr, num_nodes=pos.size(0))


class CustomNet(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = SplineConv(d.num_features, 32, dim=2, kernel_size=5)
        self.conv2 = SplineConv(32, 64, dim=2, kernel_size=5)
        self.fc1 = torch.nn.Linear(64, 128)
        self.fc2 = torch.nn.Linear(128, 1)

    def forward(self, data):
        data.x = F.elu(self.conv1(data.x, data.edge_index, data.edge_attr))
        weight = normalized_cut_2d(data.edge_index, data.pos)
        cluster = graclus(data.edge_index, weight)
        data = max_pool(cluster, data)

        data.x = F.elu(self.conv2(data.x, data.edge_index, data.edge_attr))
        weight = normalized_cut_2d(data.edge_index, data.pos)
        cluster = graclus(data.edge_index, weight)
        x, batch = max_pool_x(cluster, data.x, data.batch)

        x = scatter_mean(x, batch, dim=0)
        x = F.elu(self.fc1(x))
        x = F.dropout(x, training=self.training)
        return F.log_softmax(self.fc2(x), dim=1)


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = NeuralNet(database, CustomNet,
               node_feature=['type', 'polarity', 'bsa',
                             'depth', 'hse', 'ic', 'pssm'],
               edge_feature=['dist'],
               target='irmsd',
               index=range(400),
               batch_size=64,
               percent=[0.8, 0.2])
model.optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
model.loss = MSELoss()

model.train(nepoch=50)

```

## h5x support

After installing  `h5xplorer`  (https://github.com/DeepRank/h5xplorer), you can execute the python file `deeprank_gnn/h5x/h5x.py` to explorer the connection graph used by DeepRank-GNN. The context menu (right click on the name of the structure) allows to automatically plot the graphs using `plotly` as shown below.

![alt-text](./h5_deeprank_gnn.png)
# Run DeepRank-GNN pre-trained model on your own data

We herein provide the pretrained model `fold6_treg_yfnat_b128_e20_lr0.001_4.pt` from the DeepRank-GNN paper (doi: 10.1101/2021.12.08.471762).

An example of code to run DeepRank-GNN on new data with the pre-trained model is provided (*test.py*).

## Usage: 
`python test.py `


You can also check Deeprank-GNN documentation: [`use-deeprank-gnn-paper-s-pretrained-model`](https://deeprank-gnn.readthedocs.io/en/latest/tutorial.train_model.html#use-deeprank-gnn-paper-s-pretrained-model)

The data used to train the model are available on [SBGrid](https://data.sbgrid.org/dataset/843/)

## Please cite

M. Réau, N. Renaud, L. C. Xue, A. M. J. J. Bonvin, DeepRank-GNN: A Graph Neural Network Framework to Learn Patterns in Protein-Protein Interfaces
bioRxiv 2021.12.08.471762; doi: https://doi.org/10.1101/2021.12.08.471762

## Details of the code

```
import glob 
import sys 
import time
import datetime 
import numpy as np

from deeprank_gnn.GraphGenMP import GraphHDF5
from deeprank_gnn.NeuralNet import NeuralNet
from deeprank_gnn.ginet import GINet
```
### Graph generation section
```
#path to the docking models in pdb format
pdb_path = '../tests/data/pdb/1ATN/' 
#path to the pssm files
pssm_path = '../tests/data/pssm/1ATN/'

GraphHDF5(pdb_path=pdb_path, pssm_path=pssm_path,
        graph_type='residue', outfile='1ATN_residue.hdf5', nproc=4)
```

### Prediction section
```
pretrained_model = 'fold6_treg_yfnat_b128_e20_lr0.001_4.pt'
gnn = GINet

#path to the graph(s)
database_test = glob.glob('./*.hdf5')

start_time = time.time()
model = NeuralNet(database_test, gnn, pretrained_model = pretrained_model)    
model.test(threshold=None)
end_time = time.time()

print ('Elapsed time: {end_time-start_time}')
```

DeepRank-GNN will generate a HDF5 file with the prediction. We provite a hdf5 to csv converter to easily read it :
`python Deeprank-GNN/deeprank_gnn/tools/hdf5_to_csv.py your.hdf5`

Further information can be found in the DeepRank-GNN online [documentation](https://deeprank-gnn.readthedocs.io/en/latest/index.html)
# Run DeepRank-GNN pre-trained model on your own data

We herein provide the pretrained model `tclass_ybio_interface_b128_e50_lr0.001_26.pth.tar` from the DeepRank-GNN paper (doi: 10.1101/2021.12.08.471762).

An example of code to run DeepRank-GNN on new data with the pre-trained model is provided (*test.py*).

## Usage: 
`python test.py `


The data used to train the model are available on [SBGrid](https://data.sbgrid.org/dataset/843/)

## Please cite

M. Réau, N. Renaud, L. C. Xue, A. M. J. J. Bonvin, DeepRank-GNN: A Graph Neural Network Framework to Learn Patterns in Protein-Protein Interfaces
bioRxiv 2021.12.08.471762; doi: https://doi.org/10.1101/2021.12.08.471762

## Details of the code

```
import glob 
import sys 
import time
import datetime 
import numpy as np

from deeprank_gnn.GraphGenMP import GraphHDF5
from deeprank_gnn.NeuralNet import NeuralNet
from deeprank_gnn.ginet import GINet
```
### Graph generation section
```
#path to the docking models in pdb format
pdb_path = '../DC/pdb/' 
#path to the pssm files
pssm_path = '../DC/pssm/'

#path to the graphs (hdf5) file 
database_test = 'biological_vs_crystal.hdf5'

#generation of the graphs (HDF5 file)
GraphHDF5(pdb_path=pdb, pssm_path=pssm, biopython=False,
              graph_type='residue', outfile=database_test, nproc=8)
```
In a benchmark mode, you can add the target values of your test set to compute the performance metrics
More details are provided in DeepRank-GNN's [online documentation](https://deeprank-gnn.readthedocs.io/en/latest/tutorial.generate_graph.html#add-your-target-values)
```
add_target(graph_path=database_test, target_name='bio_interface',
           target_list='bio_interfaces.txt')

```
### Prediction section
```
pretrained_model = 'tclass_ybio_interface_b128_e50_lr0.001_26.pth.tar'
gnn = GINet

start_time = time.time()
model = NeuralNet(database_test, gnn, pretrained_model = pretrained_model", target='bio_interface')

# if you have no target values (prediction mode), remove the 'target' argument AND THE METRICS COMPUTATION LINES:
# model = NeuralNet(database_test, GINet, pretrained_model = "tclass_ybio_interface_b128_e50_lr0.001_26.pth.tar")
model.test(hdf5, hdf5='prediction_phy_non-phy.hdf5')
end_time = time.time()

print ('Elapsed time: {end_time-start_time}')
```
# Evaluation of the performance 
```
test_metrics = model.get_metrics('test', threshold = 1.0)
print("accuracy:",test_metrics.accuracy)
print("specificity:",test_metrics.specificity)
print("sensitivity:",test_metrics.sensitivity)
print("precision:",test_metrics.precision)
print("FPR:",test_metrics.FPR)
print("FNR:",test_metrics.FNR)
```

DeepRank-GNN will generate a HDF5 file with the prediction. We provite a hdf5 to csv converter to easily read it :
`python Deeprank-GNN/deeprank_gnn/tools/hdf5_to_csv.py your.hdf5`

Further information can be found in the DeepRank-GNN online [documentation](https://deeprank-gnn.readthedocs.io/en/latest/index.html)
############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

#. you have a question;
#. you think you may have found a bug (including unexpected behavior);
#. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

#. use the search functionality `here <https://github.com/NLESC-JCER/pyCHAMP/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/NLESC-JCER/pyCHAMP/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
#. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

#. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
#. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
#. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__ and `here <https://help.github.com/articles/syncing-a-fork/>`__);
#. make sure the existing tests still work by running ``python setup.py test``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the pyCHAMP repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
###############################################################################
Contributor Covenant Code of Conduct
###############################################################################

Our Pledge
**********

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

Our Standards
*************

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

Our Responsibilities
********************

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
*****

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
***********

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at n.renaud@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
***********

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org>`_, version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
Graph and Graph Generator
=====================================

Graphs
--------------------------------------

.. automodule:: deeprank_gnn.Graph
    :members:
    :undoc-members:

Residue Graphs
--------------------------------------

.. automodule:: deeprank_gnn.ResidueGraph
    :members:
    :undoc-members:



Graph Generator
---------------------------------------


.. automodule:: deeprank_gnn.GraphGen
    :members:
    :undoc-members:

Parrallel Graph Generator
---------------------------------------

.. automodule:: deeprank_gnn.GraphGenMP
    :members:
    :undoc-members:
Graph Neural Networks
=====================================

GINet layer
--------------------------------------

Graph Interaction Networks layer

This layer is inspired by Sazan Mahbub et al. "EGAT: Edge Aggregated Graph Attention Networks and Transfer Learning Improve Protein-Protein Interaction Site Prediction", BioRxiv 2020

1) Create edges feature by concatenating node feature

.. math::
    e_{ij} = LeakyReLu (a_{ij} * [W * x_i || W * x_j])
    
2) Apply softmax function, in order to learn to consider or ignore some neighboring nodes

.. math::
    \alpha_{ij}  = softmax(e_{ij})
    
3) Sum over the nodes (no averaging here)

.. math::
    z_i = \sum_j (\alpha_{ij} * Wx_j + b_i)
    
Herein, we add the edge feature to the step 1)

.. math::
    e_{ij} = LeakyReLu (a_{ij} * [W * x_i || W * x_j || We * edge_{attr} ])


.. automodule:: deeprank_gnn.ginet
    :members:
    :undoc-members:


Fout Net layer
---------------------------------------

This layer is described by eq. (1) of "Protein Interface Predition using Graph Convolutional Network", by Alex Fout et al. NIPS 2018

.. math::
    z = x_i * Wc + 1 / Ni Sum_j x_j * Wn + b

.. automodule:: deeprank_gnn.foutnet
    :members:
    :undoc-members:

sGraphAttention (sGAT) layer
---------------------------------------

This is a new layer that is similar to the graph attention network but simpler

.. math::   
    z_i = 1 / Ni Sum_j a_ij * [x_i || x_j] * W + b_i

|| is the concatenation operator: [1,2,3] || [4,5,6] = [1,2,3,4,5,6] Ni is the number of neighbor of node i Sum_j runs over the neighbors of node i :math:`a_ij` is the edge attribute between node i and j

.. automodule:: deeprank_gnn.sGat
    :members:
    :undoc-members:

Neural Network Training Evaluation and Test
=====================================

Neural Network
-------------------------------------
.. automodule:: deeprank_gnn.NeuralNet
    :members:
    :undoc-members:
.. DeepRank-GNN documentation master file, created by
   sphinx-quickstart on Wed May 12 11:56:41 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DeepRank-GNN's documentation!
========================================

.. toctree::
   :maxdepth: 2
   :caption: DeepRank-GNN:
   
   overview
   installation

.. toctree::
   :maxdepth: 2
   :caption: Tutorial:

   tutorial.generate_graph
   tutorial.train_model
   tutorial.advanced
   
.. toctree::
   :maxdepth: 2
   :caption: API

   api.graph_generator
   api.graph_neural_etwork
   api.neural_network

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Overview
=========================

DeepRank-GNN is a framework that converts PPI interfaces into graphs and uses those to learn interaction patterns using Graph Neural Networks.

The framework is designed to be adapted to any PPI related research project.

.. figure:: ../deeprank_gnn.png
    :align: center

Deeprank-GNN works in a two step process:

1) Graph generation 
-------------------------

Use your own protein-protein complexes to generate PPI interface graphs.

Go to :ref:`Creating Graphs<Creating Graphs>` 

2) Model Training
-------------------------

Tune and train your Graph Neural Network and make your own predictions.

Go to :ref:`Training a module<Training a module>` 
 

Motivations
=========================

Protein-protein interactions (PPIs) are essential in all cellular processes of living organisms
including cell growth, structure, communication, protection and death. Acquiring knowledge on PPI is
fundamental to understand normal and altered physiological processes and propose solutions to
restore them. In the past decades, a large number of PPI structures have been solved by experimental
approaches (e.g., X-ray crystallography, nuclear magnetic resonance, cryogenic electron microscopy).
Given the remarkable success of Convolutional Neural Network (CNN) in retrieving patterns in images [1]_,
CNN architectures have been developed to learn interaction patterns in PPI interfaces [2]_, [3]_.

CNNs however come with major limitations: First, they are sensitive to the input PPI
orientation, which may require data augmentation (i.e. multiple rotations of the input data) for the
network to forget about the orientation in the learning process; second, the size of the 3D grid is
unique for all input data, which does not reflect the variety in interface sizes observed in experimental
structures and may be problematic for large interfaces that do not fit inside the predefined grid size.
A solution to this problem is to use instead Graph Neural networks (GNN). 
By definition, graphs are non-structured geometric structures and do not hold orientation information. They are rotational invariant and can easily represent interfaces of varying sizes. 

Building up on our previous tool DeepRank(https://github.com/DeepRank/deeprank) that maps atomic and residue-level features from PPIs to 3D grids and applies 3D CNNs to learn problem-specific interaction patterns, we present here Deeprank-GNN. Deeprank-GNN converts PPI interfaces into graphs and uses those to learn interaction patterns. 

DeepRank-GNN is a framework than can be easily used by the community and adapted to any topic involving 
PPI interactions. The framework allows users to define their own graph neural network, features and target values. 

.. [1] Krizhevsky A, Sutskever I, Hinton GE, ImageNet classification with deep convolutional neural networks. Adv Neural Inf Process Syst 25, 2012

.. [2] Renaud N, Geng C, Georgievska S, Ambrosetti F, Ridder L, Marzella D, Bonvin A, Xue L, DeepRank: A deep learning framework for data mining 3D protein-protein interfaces, bioRxiv, 2021.01.29.425727

.. [3] Wang X, Terashi G, Christoffer CW, Zhu M, Kihara D. Protein docking model evaluation by 3D deep convolutional neural networks. Bioinformatics. 2020 ;36(7):2113-2118.
          
Advanced 
===========================

Design your own GNN layer
---------------------------

Example:

>>> import torch
>>> from torch.nn import Parameter
>>> import torch.nn.functional as F
>>> import torch.nn as nn
>>> 
>>> from torch_scatter import scatter_mean
>>> from torch_scatter import scatter_sum
>>> 
>>> from torch_geometric.utils import remove_self_loops, add_self_loops, softmax
>>> 
>>> # torch_geometric import
>>> from torch_geometric.nn.inits import uniform
>>> from torch_geometric.nn import max_pool_x
>>> 
>>> class GINet_layer(torch.nn.Module):
>>> 
>>>     def __init__(self,
>>>                  in_channels,
>>>                  out_channels,
>>>                  number_edge_features=1,
>>>                  bias=False):
>>> 
>>>         super(EGAT, self).__init__()
>>> 
>>>         self.in_channels = in_channels
>>>         self.out_channels = out_channels
>>> 
>>>         self.fc = nn.Linear(self.in_channels, self.out_channels, bias=bias)
>>>         self.fc_edge_attr = nn.Linear(number_edge_features, 3, bias=bias)
>>>         self.fc_attention = nn.Linear(2 * self.out_channels + 3, 1, bias=bias)
>>>         self.reset_parameters()
>>>         
>>>     def reset_parameters(self):
>>> 
>>>         size = self.in_channels
>>>         uniform(size, self.fc.weight)
>>>         uniform(size, self.fc_attention.weight)
>>>         uniform(size, self.fc_edge_attr.weight)
>>>         
>>>     def forward(self, x, edge_index, edge_attr):
>>> 
>>>         row, col = edge_index
>>>         num_node = len(x)
>>>         edge_attr = edge_attr.unsqueeze(
>>>             -1) if edge_attr.dim() == 1 else edge_attr
>>> 
>>>         xcol = self.fc(x[col])
>>>         xrow = self.fc(x[row])
>>>         
>>>         ed = self.fc_edge_attr(edge_attr)
>>>         # create edge feature by concatenating node features
>>>         alpha = torch.cat([xrow, xcol, ed], dim=1)
>>>         alpha = self.fc_attention(alpha)
>>>         alpha = F.leaky_relu(alpha)
>>>         
>>>         alpha = F.softmax(alpha, dim=1)
>>>         h = alpha * xcol 
>>>         
>>>         out = torch.zeros(
>>>             num_node, self.out_channels).to(alpha.device)
>>>         z = scatter_sum(h, row, dim=0, out=out)
>>> 
>>>         return z
>>>     
>>>     def __repr__(self):
>>>         return '{}({}, {})'.format(self.__class__.__name__,
>>>                                    self.in_channels,
>>>                                    self.out_channels)

Design your  own neural network architecture
---------------------------

The provided example makes use of **internal edges** and **external edges**.

We perform convolutions on the internal and external edges independently by providing the following data to the convolution layers: 

- data_ext.internal_edge_index, data_ext.internal_edge_attr for the **internal** edges 

- data_ext.edge_index, data_ext.edge_attr for the **external** edges

.. note::  
The GNN class **must** be initialized with 3 arguments: 
- input_shape 
- output_shape 
- input_shape_edge 

>>> class GINet(torch.nn.Module):
>>>     def __init__(self, input_shape, output_shape = 1, input_shape_edge = 1):
>>>         super(GINet_layer, self).__init__()
>>>         self.conv1 = GINet_layer(input_shape, 16)
>>>         self.conv2 = GINet_layer(16, 32)
>>> 
>>>         self.conv1_ext = GINet_layer(input_shape, 16, input_shape_edge)
>>>         self.conv2_ext = GINet_layer(16, 32, input_shape_edge)
>>> 
>>>         self.fc1 = nn.Linear(2*32, 128)
>>>         self.fc2 = nn.Linear(128, output_shape)
>>>         self.clustering = 'mcl'
>>>         self.dropout = 0.4
>>> 
>>>     def forward(self, data):
>>>         act = F.relu
>>>         data_ext = data.clone()
>>> 
>>>         # INTER-PROTEIN INTERACTION GRAPH
>>>         # first conv block                                                                                                                                                  
>>>         data.x = act(self.conv1(
>>>             data.x, data.edge_index, data.edge_attr))
>>>         cluster = get_preloaded_cluster(data.cluster0, data.batch)
>>>         data = community_pooling(cluster, data)
>>> 
>>>         # second conv block                                                                                                                                                    
>>>         data.x = act(self.conv2(
>>>             data.x, data.edge_index, data.edge_attr))
>>>         cluster = get_preloaded_cluster(data.cluster1, data.batch)
>>>         x, batch = max_pool_x(cluster, data.x, data.batch)
>>> 
>>>         # INTRA-PROTEIN INTERACTION GRAPH
>>>         # first conv block                                                                                                                                                  
>>>         data_ext.x = act(self.conv1_ext(
>>>             data_ext.x, data_ext.internal_edge_index, data_ext.internal_edge_attr))
>>>         cluster = get_preloaded_cluster(data_ext.cluster0, data_ext.batch)
>>>         data_ext = community_pooling(cluster, data_ext)
>>> 
>>>         # second conv block                                                                                                                                                    
>>>         data_ext.x = act(self.conv2_ext(
>>>             data_ext.x, data_ext.internal_edge_index, data_ext.internal_edge_attr))
>>>         cluster = get_preloaded_cluster(data_ext.cluster1, data_ext.batch)
>>>         x_ext, batch_ext = max_pool_x(cluster, data_ext.x, data_ext.batch)
>>> 
>>>         # FC                                                                                                                                                         
>>>         x = scatter_mean(x, batch, dim=0)
>>>         x_ext = scatter_mean(x_ext, batch_ext, dim=0)
>>> 
>>>         x = torch.cat([x, x_ext], dim=1)
>>>         x = act(self.fc1(x))
>>>         x = F.dropout(x, self.dropout, training=self.training)
>>>         x = self.fc2(x)
>>> 
>>>         return x

Use your GNN architecture in Deeprank-GNN
---------------------------

>>> model = NeuralNet(database, GINet,
>>>                node_feature=node_feature,
>>>                edge_feature=edge_feature,
>>>                target=target,
>>>                task=task,
>>>                lr=lr,
>>>                batch_size=batch_size,
>>>                shuffle=shuffle,
>>>                percent=[0.8, 0.2])
Installation
=========================

Via Python Package
-----------------------------

The latest release of DeepRank-GNN can be installed using the pypi package manager with :

``pip install deeprank-gnn``

If you are planning to only use DeepRank-GNN this is the quickest way to obtain the code


Via GitHub
-------------

For user who would like to contribute, the code is hosted on GitHub_ (https://github.com/DeepRank/DeepRank-GNN)

.. _GitHub: https://github.com/DeepRank/DeepRank-GNN

To install the code

 * clone the repository ``git clone https://github.com/DeepRank/DeepRank-GNN.git``
 * go there ``cd DeepRank-GNN``
 * install the module ``pip install -e ./``

You can then test the installation :

 * ``cd test``
 * ``pytest``

.. note::
  Ensure that at least PyTorch 1.5.0 is installed:
  
  ``python -c "import torch; print(torch.__version__)``
  
  In case PyTorch Geometric installation fails, refer to the installation guide:  https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html 



.. _Creating Graphs:

Creating Graphs
=====================================

Deeprank-GNN automatically generates a residue-level graphs of **protein-protein interfaces** in which nodes correspond to a single residue, and 2 types of edges are defined:
  
- **External edges** connect 2 residues (nodes) of chain A and B if they have at least 1 pairwise atomic distance **< 8.5 A** (Used for to define neighbors)
  
- **Internal edges** connect 2 residues (nodes) within a chain if they have at least 1 pairwise atomic distance **< 3 A** (Used to cluster nodes)


.. warning::
  The graph generation requires an ensemble of PDB files containing two chains: chain **A** and chain **B**. 
  
  You can provide PSSM matrices to compute evolutionary conservation node features. Some pre-calculated PSSM matrices can be downloaded from http://3dcons.cnb.csic.es/.
  A ``3dcons_to_deeprank_pssm.py`` converter can be found in the ``tool`` folder to convert the 3dcons PSSM format into the Deeprank-GNN PSSM format. **Make sure the sequence numbering matches the PDB residues numbering.**
  
  
 
By default, the following features are assigned to each node of the graph :
  
- **pos**: xyz coordinates

- **chain**: chain ID

- **charge**: residue charge

- **polarity**: apolar/polar/neg_charged/pos_charged (one hot encoded)

- **bsa**: buried surface are

- **type**: residue type (one hot encoded)

The following features are computed if PSSM data is provided :

- **pssm**: pssm score for each residues

- **cons**: pssm score of the residue

- **ic**: information content of the PSSM (~Shannon entropy)

The following features are optional, and computed only if Biopython is used (see next example) :

- **depth**: average atom depth of the atoms in a residue (distance to the surface)

- **hse**: half sphere exposure

Generate your graphs 
-------------------------------------

Note that the pssm information is used to compute the **pssm**, **cons** and **ic** node features and is optional.

In this example, all features are computed

>>> from deeprank_gnn.GraphGenMP import GraphHDF5
>>>
>>> pdb_path = './data/pdb/1ATN/'
>>> pssm_path = './data/pssm/1ATN/'
>>>
>>> GraphHDF5(pdb_path=pdb_path, pssm_path=pssm_path, biopython=True,
>>>          graph_type='residue', outfile='1ATN_residue.hdf5', nproc=4)

In this example, the biopython features (hse and depth) are ignored

>>> from deeprank_gnn.GraphGenMP import GraphHDF5
>>>
>>> pdb_path = './data/pdb/1ATN/' # path to the docking model in PDB format
>>> pssm_path = './data/pssm/1ATN/' # path to the pssm files
>>>
>>> GraphHDF5(pdb_path=pdb_path, pssm_path=pssm_path, 
>>>          graph_type='residue', outfile='1ATN_residue.hdf5', nproc=4)

In this example, the biopython features (hse and depth) and the PSSM information are ignored

>>> from deeprank_gnn.GraphGenMP import GraphHDF5
>>>
>>> pdb_path = './data/pdb/1ATN/'
>>>
>>> GraphHDF5(pdb_path=pdb_path, 
>>>          graph_type='residue', outfile='1ATN_residue.hdf5', nproc=4)

Add your target values
-------------------------------------

Use the CustomizeGraph class to add target values to the graphs. 

If you are benchmarking docking models, go to the **next section**.

>>> from deeprank_gnn.GraphGenMP import GraphHDF5
>>> from deeprank_gnn.tools.CustomizeGraph import add_target
>>> 
>>> pdb_path = './data/pdb/1ATN/'
>>> pssm_path = './data/pssm/1ATN/'
>>>
>>> GraphHDF5(pdb_path=pdb_path, pssm_path=pssm_path,
>>>          graph_type='residue', outfile='1ATN_residue.hdf5', nproc=4)
>>>
>>> add_target(graph_path='.', target_name='new_target',
>>>            target_list='list_of_target_values.txt')

.. note::
  The list of target values should respect the following format:
  
  ``model_name_1 0``
  
  ``model_name_2 1``
  
  ``model_name_3 0``
  
  ``model_name_4 0``
  
  if your use other separators (eg. ``,``, ``;``, ``tab``) use the ``sep`` argument:
  
  >>> add_target(graph_path=graph_path, target_name='new_target', 
  >>>            target_list='list_of_target_values.txt', sep=',')
  
  
Docking benchmark mode 
-------------------------------------

In a docking benchmark mode, you can provide the path to the reference structures in the graph generation step. Knowing the reference structure, the following target values will be automatically computed, based on CAPRI quality criteria [1]_,  and assigned to the graphs : 

- **irmsd**: interface RMSD (RMSD between the superimposed interface residues)

- **lrmsd**: ligand RMSD (RMSD between chains B given that chains A are superimposed)

- **fnat**: fraction of native contacts

- **dockQ**: see Basu et al., "DockQ: A Quality Measure for Protein-Protein Docking Models", PLOS ONE, 2016

- **bin_class**: binary classification (0: ``irmsd >= 4 A``, 1: ``RMSD < 4A``)

- **capri_classes**: 1: ``RMSD < 1A``, 2: ``RMSD < 2A``, 3: ``RMSD < 4A``, 4: ``RMSD < 6A``, 0: ``RMSD >= 6A``

>>> from deeprank_gnn.GraphGenMP import GraphHDF5
>>>
>>> pdb_path = './data/pdb/1ATN/'
>>> pssm_path = './data/pssm/1ATN/'
>>> ref = './data/ref/1ATN/'
>>>
>>> GraphHDF5(pdb_path=pdb_path, ref_path=ref, pssm_path=pssm_path,
>>>          graph_type='residue', outfile='1ATN_residue.hdf5', nproc=4)

.. note::  
  The different input files must respect the following nomenclature:
  
   - PDB files: ``1ATN_xxx.pdb`` (xxx may be replaced by anything)
   - PSSM files: ``1ATN.A.pdb.pssm 1ATN.B.pdb.pssm`` or ``1ATN.A.pssm 1ATN.B.pssm``
   - Reference PDB files: ``1ATN.pdb``
   


.. [1] 
  Lensink MF, Méndez R, Wodak SJ, Docking and scoring protein complexes: CAPRI 3rd Edition. Proteins. 2007
.. _Training a module:

Training a module
=============================================


1. Select node and edge features
---------------------------------------------

1.1. Edge feature:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **dist**: distance between nodes

>>> edge_feature=['dist']


.. note::  
  **External edges** connect 2 residues of chain A and B if they have at least 1 pairwise atomic distance **< 8.5 A** (Used for to define neighbors)
  
  **Internal edges** connect 2 residues within a chain if they have at least 1 pairwise atomic distance **< 3 A** (Used to cluster nodes)



1.2. Node features:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **pos**: xyz coordinates

- **chain**: chain ID

- **charge**: residue charge

- **polarity**: apolar/polar/neg_charged/pos_charged (one hot encoded)

- **bsa**: buried surface are

- **pssm**: pssm score for each residues

- **cons**: pssm score of the residue

- **ic**: information content of the PSSM (~Shannon entropy)

- **type**: residue type (one hot encoded)

- **depth** (opt. in graph. gen. step): average atom depth of the atoms in a residue (distance to the surface)

- **hse** (opt. in graph. gen. step): half sphere exposure



>>> node_feature=['type', 'polarity', 'bsa',
>>>               'ic', 'pssm']


2. Select the target (benchmarking mode)
---------------------------------------------

When using Deeprank-GNN in a benchmarking mode, you must specify your target (often referred to as Y).
The target values are pre-calculated during the Graph generation step **if a reference structure is provided**.

**Pre-calculated targets:** 

- **irmsd**: interface RMSD (RMSD between the superimposed interface residues)

- **lrmsd**: ligand RMSD (RMSD between chains B given that chains A are superimposed)

- **fnat**: fraction of native contacts

- **dockQ**: see Basu et al., "DockQ: A Quality Measure for Protein-Protein Docking Models", PLOS ONE, 2016

- **bin_class**: binary classification (0: ``irmsd >= 4 A``, 1: ``RMSD < 4A``)

- **capri_classes**: 1: ``RMSD < 1A``, 2: ``RMSD < 2A``, 3: ``RMSD < 4A``, 4: ``RMSD < 6A``, 0: ``RMSD >= 6A``

.. note::  
 In classification mode (i.e. task="class") you must provide the list of target classes to the NeuralNet (e.g. classes=[1,2,3,4])


>>> target='irmsd'

3. Select hyperparameters
---------------------------------------------

- Regression ('reg') of classification ('class') mode

>>> task='reg' 

- Batch size

>>> batch_size=64

- Shuffle the training dataset

>>> shuffle=True

- Learning rate:

>>> lr=0.001

4. Load the network
---------------------------------------------

This step requires pre-calculated graphs in hdf5 format. 

The user may :

- option 1: input a unique dataset and chose to automatically split it into a training set and an evaluation set

- option 2: input distinct training/evaluation/test sets

4.1. Option 1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>>> from deeprank_gnn.NeuralNet import NeuralNet
>>> from deeprank_gnn.ginet import GINet
>>>
>>> database = './1ATN_residue.hdf5'
>>>
>>> model = NeuralNet(database, GINet,
>>>                node_feature=node_feature,
>>>                edge_feature=edge_feature,
>>>                target=target,
>>>                task=task, 
>>>                lr=lr,
>>>                batch_size=batch_size,
>>>                shuffle=shuffle,
>>>                percent=[0.8, 0.2])
>>>

.. note::  
 The *percent* argument is required to split the input dataset into a training set and a test set. Using ``percent=[0.8, 0.2]``, 80% of the input dataset will constitute the  training set, 20% will constitute the evaluation set. 

4.2. Option 2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>>> from deeprank_gnn.NeuralNet import NeuralNet
>>> from deeprank_gnn.ginet import GINet
>>> import glob 
>>>
>>> # load train dataset
>>> database_train = glob.glob('./hdf5/train*.hdf5')
>>> # load validation dataset
>>> database_eval = glob.glob('./hdf5/eval*.hdf5')
>>> # load test dataset
>>> database_test = glob.glob('./hdf5/test*.hdf5')
>>> 
>>> model = NeuralNet(database_train, GINet,
>>>                node_feature=node_feature,
>>>                edge_feature=edge_attr,
>>>                target=target,
>>>                task=task, 
>>>                lr=lr,
>>>                batch_size=batch_size,
>>>                shuffle=shuffle,
>>>                database_eval = database_eval)

5. Train the model 
---------------------------------------------

- example 1:

train the network, perform 50 epochs

>>> model.train(nepoch=50, validate=False)

- example 2:

train the model, evaluate the model at each epoch, save the best model (i.e. the model with the lowest loss), and write all predictions to ``output.hdf5``

>>> model.train(nepoch=50, validate=True, save_model='best', hdf5='output.hdf5')

.. warning::
 The ``last`` model is saved by default.
 
 When setting ``save_model='best'``, a model that is associated with a lower loss than those generated in the previous epochs will be saved. By default, the epoch number is included in the output name not to write over intermediate models.

6. Analysis
---------------------------------------------

6.1. Plot the loss evolution over the epochs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>>> model.plot_loss(name='plot_loss')

6.2 Analyse the performance in benchmarking conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following analysis only apply if a reference structure was provided during the graph generation step.

6.2.1. **Plot accuracy evolution**

>>> model.plot_acc(name='plot_accuracy')

6.2.2. **Plot hitrate**

A threshold value is required to binarise the target value

>>> model.plot_hit_rate(data='eval', threshold=4.0, mode='percentage', name='hitrate_eval')

6.2.3. **Get various metrics**

The following metrics can be easily computed: 

**Classification metrics:**

- **sensitivity**: Sensitivity, hit rate, recall, or true positive rate

- **specificity**: Specificity or true negative rate

- **precision**: Precision or positive predictive value

- **NPV**: Negative predictive value

- **FPR**: Fall out or false positive rate

- **FNR**: False negative rate

- **FDR**: False discovery rate

- **accuracy**: Accuracy

- **auc()**: AUC

- **hitrate()**: Hit rate

**Regression metrics:**

- **explained_variance**: Explained variance regression score function

- **max_error**: Max_error metric calculates the maximum residual error

- **mean_abolute_error**: Mean absolute error regression loss

- **mean_squared_error**: Mean squared error regression loss

- **root_mean_squared_error**: Root mean squared error regression loss

- **mean_squared_log_error**: Mean squared logarithmic error regression loss

- **median_squared_log_error**: Median absolute error regression loss

- **r2_score**: R^2 (coefficient of determination) regression score function

.. note::  
  All classification metrics can be calculated on continuous targets as soon as a threshold is provided to binarise the data.

>>> train_metrics = model.get_metrics('train', threshold = 4.0)
>>> print('training set - accuracy:', train_metrics.accuracy)
>>> print('training set - sensitivity:', train_metrics.sensitivity)
>>> 
>>> eval_metrics = model.get_metrics('eval', threshold = 4.0)
>>> print('evaluation set - accuracy:', eval_metrics.accuracy)
>>> print('evaluation set - sensitivity:', eval_metrics.sensitivity)

7. Save the model/network
---------------------------------------------

>>> model.save_model("model_backup.pth.tar")

8. Test the model on an external dataset
---------------------------------------------

8.1. On a loaded model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>>> model.test(database_test, threshold=4.0)

8.2. On a pre-trained model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>>> from deeprank_gnn.NeuralNet import NeuralNet
>>> from deeprank_gnn.ginet import GINet
>>>  
>>> database_test = './1ATN_residue.hdf5'
>>>  
>>> model = NeuralNet(database_test, GINet, pretrained_model = "model_backup.pth.tar")
>>> model.test(database_test)
>>>  
>>> test_metrics = model.get_metrics('test', threshold = 4.0)
>>> print(test_metrics.accuracy)

In short 
=============================================

>>> from deeprank_gnn.NeuralNet import NeuralNet
>>> from deeprank_gnn.ginet import GINet
>>>
>>> database = './1ATN_residue.hdf5'
>>>
>>> edge_feature=['dist']
>>> node_feature=['type', 'polarity', 'bsa',
>>>               'depth', 'hse', 'ic', 'pssm']
>>> target='irmsd'
>>> task='reg' 
>>> batch_size=64
>>> shuffle=True
>>> lr=0.001
>>>
>>> model = NeuralNet(database, GINet,
>>>                node_feature=node_feature,
>>>                edge_feature=edge_feature,
>>>                target=target,
>>>                index=None,
>>>                task=task, 
>>>                lr=lr,
>>>                batch_size=batch_size,
>>>                shuffle=shuffle,
>>>                percent=[0.8, 0.2])
>>>
>>> model.train(nepoch=50, validate=True, save_model='best', hdf5='output.hdf5')
>>> model.plot_loss(name='plot_loss')
>>> 
>>> train_metrics = model.get_metrics('train', threshold = 4.0)
>>> print('training set - accuracy:', train_metrics.accuracy)
>>> print('training set - sensitivity:', train_metrics.sensitivity)
>>> 
>>> eval_metrics = model.get_metrics('eval', threshold = 4.0)
>>> print('evaluation set - accuracy:', eval_metrics.accuracy)
>>> print('evaluation set - sensitivity:', eval_metrics.sensitivity)
>>> 
>>> model.save_model("model_backup.pth.tar")
>>> #model.test(database_test, threshold=4.0)

Using default settings 

>>> from deeprank_gnn.NeuralNet import NeuralNet
>>> from deeprank_gnn.ginet import GINet
>>>
>>> database = glob.glob('./hdf5/*_train.hdf5')
>>> dataset_test = glob.glob('./hdf5/*_test.hdf5')
>>>
>>> target='irmsd'
>>>
>>> model = NeuralNet(database, GINet,
>>>                target=target,
>>>                percent=[0.8, 0.2])
>>>
>>> model.train(nepoch=50, validate=True, save_model='best', hdf5='output.hdf5')
>>> model.plot_loss(name='plot_loss')
>>> 
>>> train_metrics = model.get_metrics('train', threshold = 4.0)
>>> print('training set - accuracy:', train_metrics.accuracy)
>>> print('training set - sensitivity:', train_metrics.sensitivity)
>>> 
>>> eval_metrics = model.get_metrics('eval', threshold = 4.0)
>>> print('evaluation set - accuracy:', eval_metrics.accuracy)
>>> print('evaluation set - sensitivity:', eval_metrics.sensitivity)
>>> 
>>> model.save_model("model_backup.pth.tar")
>>> model.test(database_test, threshold=4.0)

Use DeepRank-GNN paper's pretrained model
=============================================

**See**: M. Réau, N. Renaud, L. C. Xue, A. M. J. J. Bonvin, "DeepRank-GNN: A Graph Neural Network Framework to Learn Patterns in Protein-Protein Interfaces", bioRxiv 2021.12.08.471762; doi: https://doi.org/10.1101/2021.12.08.471762

You can get the `pre-trained model <https://github.com/DeepRank/Deeprank-GNN/tree/master/paper_pretrained_models/>`_ from DeepRank-GNN `github repository <https://github.com/DeepRank/Deeprank-GNN/>`_

>>> import glob 
>>> import sys 
>>> import time
>>> import datetime 
>>> import numpy as np
>>> 
>>> from deeprank_gnn.GraphGenMP import GraphHDF5
>>> from deeprank_gnn.NeuralNet import NeuralNet
>>> from deeprank_gnn.ginet import GINet
>>>
>>> ### Graph generation section
>>> pdb_path = '../tests/data/pdb/1ATN/'
>>> pssm_path = '../tests/data/pssm/1ATN/'
>>>
>>> GraphHDF5(pdb_path=pdb_path, pssm_path=pssm_path,
>>>         graph_type='residue', outfile='1ATN_residue.hdf5', nproc=4)
>>> 
>>> ### Prediction section
>>> gnn = GINet
>>> pretrained_model = 'fold6_treg_yfnat_b128_e20_lr0.001_4.pt'
>>> database_test = glob.glob('1ATN_residue.hdf5')
>>> 
>>> start_time = time.time()
>>> model = NeuralNet(database_test, gnn, pretrained_model = pretrained_model)    
>>> model.test(threshold=None)
>>> end_time = time.time()
>>> print ('Elapsed time: {end_time-start_time}')
>>> 
>>> ### The output is automatically stored in **test_data.hdf5** 
         
.. note::  
 For storage convenience, all predictions are stored in a HDF5 file. A converter from HDF5 to csv is provided in the **tools** directory
