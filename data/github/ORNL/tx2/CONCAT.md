# TX2

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPI version](https://badge.fury.io/py/tx2.svg)](https://badge.fury.io/py/tx2)
[![JOSS status](https://joss.theoj.org/papers/b7c161917e5a31af052a597bf98f0e94/status.svg)](https://joss.theoj.org/papers/b7c161917e5a31af052a597bf98f0e94)

Welcome to TX2! This library is intended to aid in the explorability and explainability of
transformer classification networks, or transformer language models with sequence classification
heads. The basic function of this library is to take a trained transformer and
test/train dataset and produce an ipywidget dashboard as seen in the screenshot below,
which can be displayed in a jupyter notebook or jupyter lab.

![screenshot]( https://raw.githubusercontent.com/ORNL/tx2/master/docs/source/screenshot.png)

NOTE: Currently this library's implementation is partially torch-dependent, and so will
not work with tensorflow/keras models - we hope to address this limitation in the future!

## Installation

You can install this package from pypi:

```bash
pip install tx2
```

NOTE: depending on the environment, it may be better to install some of the dependencies separately before
pip installing tx2, e.g. in conda:
```bash
conda install pytorch-gpu pandas scikit-learn matplotlib ipywidgets "numpy<=1.20" -c conda-forge
```

If you do not have access to a GPU on your machine, install the regular pytorch
package:
```bash
conda install pytorch pandas scikit-learn matplotlib ipywidgets "numpy<=1.20"
```

Note that `pytorch-gpu` can only be found in the `conda-forge` channel.

## Examples

Example jupyter notebooks demonstrating and testing the usage of this library can be
found in the [examples
folder](https://github.com/ORNL/tx2/tree/master/examples).

Note that these notebooks can take a while to run the first time, especially
if a GPU is not in use.

Packages you'll need to install for the notebooks to work (in addition to the
conda installs above):

```bash
pip install tqdm transformers~=4.1.1
```

Running through each full notebook will produce the ipywidget dashboard near the
end.

The tests in this repository do not depend on transformers, so raw library
functionality can be tested by running `pytest` in the project root.

## Documentation

The documentation can be viewed at [https://ornl.github.io/tx2/](https://ornl.github.io/tx2/).

The documentation can also be built from scratch with sphinx as needed.

Install all required dependencies: 
```bash
pip install -r requirements.txt
```

Build documentation:

```bash
cd docs
make html
```

The `docs/build/html` folder will now contain an `index.html`

Two notebooks demonstrating the dashboard and how to use TX2 are included
in the `examples` folder, highlighting the default and custom approaches
as discussed in the Basic Usage page of the documentation.

## Citation

To cite usage of TX2 in a publication, the DOI for this code is [https://doi.org/10.21105/joss.03652](https://doi.org/10.21105/joss.03652)

bibtex:
```
@article{Martindale2021,
  doi = {10.21105/joss.03652},
  url = {https://doi.org/10.21105/joss.03652},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3652},
  author = {Nathan Martindale and Scott L. Stewart},
  title = {TX$^2$: Transformer eXplainability and eXploration},
  journal = {Journal of Open Source Software}
}
```
# Contributing to TX2

Help in improving TX2 is welcome! 

If you find a bug or think of an enhancement/improvement you would like to see,
feel free to fill out an appropriate
[issue](https://github.com/ORNL/tx2/issues/new/choose).

If you have a question, double check that it's not covered in our
[documentation](https://ornl.github.io/tx2).

For questions not answered by the docs or anything else that might not fit into one
of the issue templates, you can start a discussion in the [dicussions
tab](https://github.com/ORNL/tx2/discussions).

You are also welcome to contact the developers directly by emailing us at
tx2-help@ornl.gov.

## Submitting a PR

If you have added a useful feature or fixed a bug, open a new pull request with
the changes.  When submitting a pull request, please describe what the pull 
request is addressing and briefly list any significant changes made. If it's in
regards to a specific issue, please include the issue number. Please check and
follow the formatting conventions below!

## Code Formatting

This project uses the [black code formatter](https://github.com/psf/black).

Any public functions and clases should be clearly documented with 
[sphinx-style docstrings](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html).
Local documentation can be generated with

```bash
cd docs
make html
```
---
title: 'TX$^2$: Transformer eXplainability and eXploration'
tags:
  - Python
  - explainability
  - natural language processing
  - deep networks
  - transformers
authors:
  - name: Nathan Martindale^[co-first author] ^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-5036-5433
    affiliation: 1
  - name: Scott L. Stewart^[co-first author] 
    orcid: 0000-0003-4320-5818
    affiliation: 1
affiliations:
 - name: Oak Ridge National Laboratory
   index: 1
date: 21 December 2021
bibliography: paper.bib
---

# Summary

The Transformer eXplainability and eXploration [@tx2], or TX$^2$ software package, is a library designed for artificial intelligence researchers to better understand the performance of transformer models [@vaswani2017attention] used for sequence classification. The tool is capable of integrating with a trained transformer model and a dataset split into training and testing populations to produce an ipywidget [@ipywidgets] dashboard with a number of visualizations to understand model performance with an emphasis on explainability and interpretability. The TX$^2$ package is primarily intended to integrate into a workflow centered around Jupyter Notebooks [@jupyternotebook], and currently assumes the use of PyTorch [@pytorch] and Hugging Face transformers library [@hf-transformers]. The dashboard includes visualization and data exploration features to aid researchers, including an interactive UMAP embedding graph [@mcinnes2018umap] to understand classification clusters, a word salience map that can be updated as researchers alter textual entries in near real time, a set of tools to understand word frequency and importance based on the clusters in the UMAP embedding graph, and a set of traditional confusion matrix analysis tools. 

# Statement of Need

Transformers, although particularly effective on a wide variety of natural language processing tasks, have the same challenge of many deep network approaches in that it is difficult to glean insight into why certain classification decisions are made [@aken2020visbert]. Various works have explored the value of analyzing the attention layers in order to provide explainability in the output of a transformer network [@vig2019multiscale]. However, analyzing attention alone can be insufficient when attempting to gain broader insight into why a transformer is performing a certain way with a specific dataset [@jain2019attention]. TX$^2$ aims to address this challenge by providing a model developer with a number of tools to explore why a certain transformer performs in a certain way for a specific dataset. This tool can help a developer determine, among other things, whether or not a specific transformer has gained a generalized understanding of the semantic meaning behind textual entries in a specific dataset. It can also help with studying the impact of language distribution shifts over time on transformer sequence classification performance.

Existing tools, such as Google PAIR's Language Interpretability Tool [@tenney2020lit], also provide a platform to use multiple visualizations to study transformer model performance. TX$^2$ differs from these tools with its emphasis on cluster analysis and easier customization of both the model interaction and dashboard itself within a Jupyter Notebook. The close integration with Jupyter Notebook is advantageous for those researchers who already rely heavily on the tools within the Jupyter ecosystem. Like the Language Interpretability Tool, TX$^2$ offers a projection map with all of the data points; however it goes further in breaking down the visual clusters and providing separate visualizations for understanding the language per cluster. Additionally, the TX$^2$ design promotes easy modification or customization depending on the researcher's needs, as researchers can completely change the presentation order of plots within the ipywidget and even add additional visualizations if desired. 

# Features

The primary visualization for the widget is a UMAP embedding graph that projects the multidimensional sequence embedding space into 2D clusters. This plots multiple controls that can be used to understand how the sequence classifier is working, including the ability to show or hide training data, highlight certain keywords, and focus on misclassifications. Below the UMAP plot, the dashboard includes a set of tools for exploring textual data including a word salience map that shows information on specific train or test data entries. The salience map serves as a proxy for word importance and is computed by recalculating the soft classifications of a particular entry in the corpus multiple times with each word individually removed. The background coloring in the map indicates the degree of impact word removal has on the classification result, with a darker background highlight corresponding to greater importance. The dashboard also includes a text entry box that is prepopulated with the text from the entry shown in the salience map. The user can use this text box to explore the impact of word addition or removal by modifying the entry. The change is reflected both in the salience map plot as well as with a change in the data point in the UMAP embedding graph.  

The dashboard also includes a set of visual clustering analysis tools. Any clustering algorithm from sklearn's [@scikit-learn] clustering module can be used to assign clusters to the data once it is projected into the UMAP embedding. The dashboard displays cluster labels, along with inter-cluster word frequency, and each words’ importance on the classification result. The salience scores for each word are calculated in aggregate for each cluster, again by iterating with the classifier while individual words are removed. There are also some sampling buttons that allow for a data example to be randomly pulled from a specific cluster so that it can be examined by the entry-specific salience map tool. Finally, it is also possible to output traditional confusion matrices as well as various evaluation scores (e.g., f1-score, accuracy, precision) as part of the dashboard.

# Integration

TX$^2$ includes two main classes: a wrapper class and a dashboard class. The wrapper class wraps around the transformer/classification model and acts as an interface between the dashboard and the transformer. The wrapper is in charge of computing and caching all the necessary data for the dashboard visualizations. The dashboard class is responsible for setting up and rendering the widget layout and handling dashboard interactivity. The flow of interactions between the TX$^2$ library and a Jupyter Notebook can be seen in \autoref{fig:interactivity}.

![Flow of interactions between a Jupyter Notebook and the tx^2 library.\label{fig:interactivity}](interaction_flow.png)

The wrapper communicates with the transformer through a set of four functions as seen in \autoref{fig:integration}. These functions include an embedding function that returns a single sequence of embeddings for each input text, a classification function that returns the predicted output class for each input text, a soft classification function that returns some output value for each class for each input text, and an encoding function that converts the text into model inputs. 

![Example of integrating a transformer with the tx^2 wrapper.\label{fig:integration}](example_interaction.png)

The default implementation for TX$^2$ assumes a huggingface pretrained model. If this use case fits the purposes of the user, they can use the default implementations for these functions. Otherwise, the user will need to redefine the functions to handle their use case while ensuring that the new functions return the necessary data and the correct format.

# Audience

The target audience for the TX$^2$ tool are machine learning or artificial intelligence researchers focused on natural language processing with transformers, and who are comfortable operating within the Jupyter ecosystem for demonstration or exploration. This open-source software is licensed under a BSD-3 clause license, is registered on [DOE Code](https://doi.org/10.11578/dc.20210129.1), and is available on [GitHub](https://github.com/ORNL/tx2). The package is also pip installable with ``pip install tx2`` with Sphinx [@sphinx] built [documentation](https://ornl.github.io/tx2/index.html). Finally, linting for this project is performed using black [@black-linter].

# Acknowledgements

The authors would like to acknowledge the US Department of Energy, National Nuclear Security Administration's Office of Defense Nuclear Nonproliferation Research and Development (NA-22) for supporting this work. 

This manuscript has been authored by UT-Battelle, LLC, under contract DE-AC05-00OR22725 with the US Department of Energy (DOE). The US government retains and the publisher, by accepting the article for publication, acknowledges that the US government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this manuscript, or allow others to do so, for US government purposes. DOE will provide public access to these results of federally sponsored research in accordance with the DOE Public Access Plan (http://energy.gov/downloads/doe-public-access-plan).

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
Visualization
=============


.. automodule:: tx2.visualization
    :autosummary:
    :autosummary-sections: Functions
    :members:
Dashboard Interface
###################

This page contains brief descriptions of what the different buttons and checkboxes in the interface do.


UMAP Sidebar
============

The graph controls and status messages for the currently selected points are in the sidebar of the UMAP plot. 
The top two status indicators respectively show the busy/ready status for backend computation and plotting work.
The model classification section shows the currently predicted label as well as the target (the colors match those used
in the plot itself.)

For the checkboxes under graph controls:

* Show training data - lightly includes the training points in the plot (larger and more transparent points), for helping determine if there's much difference between the training and testing distributions.
* Visual cluster nums - display the numerical labels for the computed clusters of the 2d projections. Theses labels are consistent with the plot titles in the cluster words bar graphs and wordclouds in the tabs below.
* Focus misclassifications - slightly blur out the points the transformer predicted correctly to make it easier to see which ones it got wrong.

The sample misclassified button will randomly select a point that the transformer got wrong.

The keyword search box allows you to type in a word and the plot will highlight in red the instances where that word appears.
The "sample from highlighted" button will randomly select a point from the highlighted instances.

Other controls
==============

The "selected datapoint index" below the UMAP plot is a dropdown menu containing the index of each point in the testing 
set. Changing this selection will update the UMAP plot to highlight the current point, and the entry text and word salience
map will update.

In the Cluster words and Word clouds tabs, the "sampling buttons" will randomly select an instance from the visual cluster
of the corresponding label in square brackets on the button.
Calc
====


.. automodule:: tx2.calc
    :autosummary:
    :autosummary-sections: Functions
    :members:Basic Usage
###########

TX2 consists of two classes: :class:`tx2.wrapper.Wrapper`
and :class:`tx2.dashboard.Dashboard`.


The wrapper class wraps around the transformer/classification
model and acts as an interface between the dashboard and the transformer.
The wrapper is in charge of computing and caching all the necessary
data for the dashboard visualizations.

The dashboard class handles setting up and rendering the widget
layout and handling dashboard interactivity.

Note that this dashboard is primarily for exploring how a transformer responds to a test
set of data, and the larger this test set, the slower the dashboard may respond and the
longer the wrapper's pre-computation steps will take.

The flow of interactions between this library and a jupyter notebook is shown below:

.. image:: interaction_flow.png

All communication between TX2 and the transformer is done entirely through a set of
four interaction functions, discussed further in the sections below.

Wrapper Setup
=============

There are two different general approaches for setting up the transformer
wrapper, depending on the level of customization needed to suit your
model. The wrapper relies on four different functions for computation:

* An **embedding function** - returns a single sequence embedding for each input text.
* A **classification function** - returns the predicted output class for each input text.
* A **soft classification function** - returns some output value for each class for each input text (essentially a non-argmaxed classification output.)
* An **encoding function** - converts text into inputs the model is expecting.

In all cases, the wrapper is instantiated, and then the wrapper's :code:`prepare()` function
must be called. This runs through all necessary data computations that the
dashboard relies on.

An example diagram of a transformer model that provides the expected data for each of these functions is shown here:

.. image:: example_interaction.png

Default Approach
----------------

In the default approach, defaults for the four functions are already handled internally, and
rely on directly passing the necessary model pieces to the :code:`wrapper`
constructor. There are three pieces the constructor expects for this
to work correctly:

1. A huggingface tokenizer (the default **encoding function** will call :code:`encode_plus` on this tokenizer)
2. A calleable huggingface language model (the default **embedding function** will take the final layer outputs of this for the first token, expected to be a :code:`[CLS]` token. Importantly, this means by default it expects a BERT transformer. Any other type will require using the custom approach below)
3. A calleable classifier that returns an output value for each class (this is directly used for the default **soft classification function**, and the default **classification function** argmaxes the output.)


An example model that would work in this approach is shown below, as in the first example notebook:

.. code-block:: python

    import torch
    from transformers import AutoModel

    class BERTClass(torch.nn.Module):
        def __init__(self):
            super(BERTClass, self).__init__()
            self.language_model = AutoModel.from_pretrained("bert-base-cased")
            self.classification_head = torch.nn.Linear(768, 20)

        def forward(self, ids, mask):
            output_1 = self.language_model(ids, mask)
            output = self.classification_head(output_1[0][:, 0, :])
            return output


To instantiate the wrapper, we pass in the data and necessary model pieces, and then call
:code:`prepare()` to run the necesary computations and cache the results.

.. code-block:: python

    from transformers import AutoTokenizer

    from tx2.wrapper import Wrapper

    # initialize
    model = BERTClass()
    tokenizer = AutoTokenizer.from_pretrained("bert-base-cased")
    train_df, test_df, encodings = # load dataframes and encodings dictionary

    # train model

    # create wrapper
    wrapper = Wrapper(
        train_texts=train_df.text,
        train_labels=train_df.target,
        test_texts=test_df.text[:2000]
        test_labels=test_df.target[:2000]
        encodings=encodings,
        classifier=model,
        language_model=model.language_model,
        tokenizer=tokenizer)
    wrapper.prepare()

Note that in the example above, we expect the dataframes to have a "text" column that contains the
input text, and a "target" column that contains the integer target class. :code:`encodings` is a
dictionary that contains class labels/names as keys, with each value as the integer representation for it,
e.g. for the 20 newsgroups dataset:

.. code-block::

    {
        'alt.atheism': 0,
        'comp.graphics': 1,
        'comp.os.ms-windows.misc': 2,
        'comp.sys.ibm.pc.hardware': 3,
        'comp.sys.mac.hardware': 4,
        'comp.windows.x': 5,
        'misc.forsale': 6,
        'rec.autos': 7,
        'rec.motorcycles': 8,
        'rec.sport.baseball': 9,
        'rec.sport.hockey': 10,
        'sci.crypt': 11,
        'sci.electronics': 12,
        'sci.med': 13,
        'sci.space': 14,
        'soc.religion.christian': 15,
        'talk.politics.guns': 16,
        'talk.politics.mideast': 17,
        'talk.politics.misc': 18,
        'talk.religion.misc': 19
    }


Custom Approach
---------------

If a different type of transformer or different way of constructing your model makes
any of the default functions infeasible or incorrect, it is possible to manually specify
any of the four functions the wrapper relies on. This can be done by defining the function
and then assigning it to the corresponding wrapper attributes:

* :attr:`tx2.wrapper.Wrapper.embedding_function`
* :attr:`tx2.wrapper.Wrapper.classification_function`
* :attr:`tx2.wrapper.Wrapper.soft_classification_function`
* :attr:`tx2.wrapper.Wrapper.encode_function`

As an example, one could change the embedding mechanism to average the output token embeddings rather than
expecting a :code:`[CLS]` token.

.. code-block:: python

    import numpy as np

    transformer = # load/train language model

    def average_embedding(inputs):
         return np.mean(transformer(inputs['input_id'], inputs['attention_mask'])[0])

    wrapper = Wrapper(...)
    wrapper.embedding_function = average_embedding
    wrapper.prepare()

Note that while the wrapper's :code:`embed()`, :code:`classify()`, and :code:`soft_clasify()`
all take an array of texts as input, the corresponding backend wrapper attributes are functions
that expect *encoded inputs*, as returned from :attr:`tx2.wrapper.Wrapper.encode_function`.
By default, if you do not specify a custom :code:`encode_function`, the wrapper runs :code:`encode_plus`
on the tokenizer specified in the constructor with the :attr:`tx2.wrapper.Wrapper.encoder_options` passed in.
The results are returned in a dictionary with :code:`"input_ids"` and :code:`"attention_mask"` as keys.

Depending on what custom functions you define determines which model pieces you do or do not need to pass to the
constructor:

* If you define a :code:`encode_function`, you do not need to pass anything to :code:`tokenizer`.
* If you define a :code:`classification_function` **and** :code:`soft_classification_function`, you do not need to pass anything to :code:`classifier`.
* If you define a :code:`embedding_function`, you do not need to pass anything to :code:`language_model`.

Input Data Flow
---------------

To help understand how custom functions fit in, below is an example of how data is converted and passed through
the wrapper when the wrapper's :code:`classify()` is called.

.. image:: wrapper_data_flow.png

1. The :func:`tx2.wrapper.Wrapper.classify` function is called with an array of texts.
2. The input texts are placed into a pytorch dataset child class and dataloader.
3. For each input text the dataset calls the  :attr:`tx2.wrapper.Wrapper.encode_function`.
4. For each batched set in the dataloader (containing the outputs from 2), the batch array of encoded inputs are passed into :attr:`tx2.wrapper.Wrapper.classification_function`.
5. Output predictions are aggregated and sent back up/returned from the :code:`classify()` call.

Dashboard Setup
===============

The dashboard class is relatively straight forward - initialize it with the prepared transformer wrapper and any
settings for which sections to display, make any desired widget alterations, and then call :code:`render()`
or manually pull the components and directly display them with :code:`IPython.display.display()`. (For more details see the
:ref:`Dashboard Widgets`.)


.. code-block:: python

    from tx2.wrapper import Wrapper
    from tx2.dashboard import Dashboard

    # load and train transformer and data

    wrapper = Wrapper(...)
    wrapper.prepare()

    dash = Dashboard(wrapper)
    dash.render()


The dashboard constructor contains six booleans which control what sections get displayed when you call :code:`render()`:

.. code-block:: python

    class Dashboard:
        def __init__(
            self,
            transformer_wrapper: wrapper.Wrapper,
            show_umap=True,
            show_salience=True,
            show_word_count=True,
            show_cluster_salience=True,
            show_cluster_sample_btns=True,
            show_wordclouds=False,
        ):

The :code:`show_wordclouds` option is :code:`False` by default as the cluster-based :code:`show_word_count` and
:code:`show_cluster_salience` tend to convey more useful and representative information than the wordclouds.

Tips
----

Note that for the plots to display correctly, you need to run the :code:`%matplotlib agg` or :code:`%matplotlib inline` magic.

For the matplotlib plots themselves to remain interactive (with zoom/pan controls), you can instead use
:code:`%matplotlib notebook`. To remove the headers from each figure, you can run an HTML magic block to magic
them away:

.. code-block::

    %%html
    <style>
    div.ui-dialog-titlebar {display: none;}
    </style>

Sometimes with :code:`%matplotlib inline`, various graphs will duplicate every time they're re-rendered, which can
be fixed by calling :code:`plt.ioff()` or using :code:`%matplotlib agg` instead.
.. _Dashboard Widgets:

Dashboard Widgets
=================

This page lists and shows the different widgets that compose
the full dashboard. While :code:`Dashboard.render()` returns
the entire dashboard by default, all components can be customized
or manually accessed in order to construct the layout together
yourself.

Initialization
--------------

Note that all widgets are initialized when a new Dashboard
instance is constructed. Any changes made to the class widgets
before calling :code:`render()` will persist in the returned
layout:

.. code-block:: python

    dash = Dashboard(wrapper)
    dash.lbl_projection_graph.value = "<h3>Entry Embeddings</h3>"
    dash.render() # altered 'entry embeddings' label will be in final layout.

The sections below list all of the different components that are set up in the dashboard class.
Any may be individually altered as desired. If you need an entirely different layout for the
dashboard, you can manually reference the components or groups and put them into your own HBox/VBox layouts.
You would manually render the results for this, rather than calling :code:`render()`:

.. code-block:: python

    dash = Dashboard(wrapper)
    my_layout = HBox([
        dash.projection_layout,
        dash.manual_text_entry_and_salience_layout,
        dash.drop_text_picker,
        dash.sampling_group
    ])
    display(my_layout)

Widgets
-------

The below images show the individual widgets in each portion of the dashboard.
All names are attributes that can be directly referenced from the wrapper instance.

.. image:: screen1_components.png

.. image:: screen2_components.png


Sections
--------

Below shows different collections of components that are wrapped in HBox/VBox
layouts - these can alo be directly referenced from the wrapper instance.

.. image:: screen1_sections.png

.. image:: screen2_sections.png
Dashboard
=========
.. automodule:: tx2.dashboard


.. autoclass:: tx2.dashboard.Dashboard
    :autosummary:
    :members:

.. TX2 documentation master file, created by
   sphinx-quickstart on Wed Nov 25 09:02:08 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. TODO -  Visualization. Utils. README

TX\ :sup:`2` Documentation
===============================

Welcome to TX\ :sup:`2`! This library is intended to aid in the explorability and explainability of
transformer classification networks, or transformer language models with sequence classification
heads. The basic function of this library is to take a trained transformer and
test/train dataset and produce an ipywidget dashboard as seen in the screenshot below,
which can be displayed in a jupyter notebook or jupyter lab.

.. image:: screenshot.png

NOTE: Currently this library's implementation is partially torch-dependent, and so will
not work with tensorflow/keras models - we hope to address this limitation in the future!

.. toctree::
   :maxdepth: 2
   :caption: Usage

   basic_usage.rst
   dashboard_interface.rst
   dashboard_widgets.rst

.. toctree::
   :maxdepth: 2
   :caption: API

   wrapper.rst
   dashboard.rst
   calc.rst
   visualization.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Wrapper
=======
.. automodule:: tx2.wrapper

..
    Notable functions
    -----------------
    .. autofunction:: tx2.wrapper.Wrapper.prepare
        :noindex:
    .. autofunction:: tx2.wrapper.Wrapper.classify
        :noindex:
    .. autofunction:: tx2.wrapper.Wrapper.soft_classify
        :noindex:
    .. autofunction:: tx2.wrapper.Wrapper.embed
        :noindex:

    Pre-computed data
    -----------------


    Full API
    --------

.. autoclass:: tx2.wrapper.Wrapper
    :autosummary:
    :members: