[![Python Test](https://github.com/Planet-AI-GmbH/tfaip/actions/workflows/python-test.yml/badge.svg)](https://github.com/Planet-AI-GmbH/tfaip/actions/workflows/python-test.yml)
[![Python Test](https://github.com/Planet-AI-GmbH/tfaip/actions/workflows/python-publish.yml/badge.svg)](https://github.com/Planet-AI-GmbH/tfaip/actions/workflows/python-publish.yml)

# _tfaip_ - A Generic and Powerful Research Framework for Deep Learning based on Tensorflow

*tfaip* is a Python-based research framework for developing, organizing, and deploying Deep Learning models powered by [Tensorflow](https://www.tensorflow.org/).
It enables to implement both simple and complex scenarios that are structured and highly configurable by parameters that can directly be modified by the command line (read the [docs](https://tfaip.readthedocs.io)).
For example, the [tutorial.full](examples/tutorial/full)-scenario for learning MNIST allows to modify the graph during training but also other hyper-parameters such as the optimizer:
```bash
export PYTHONPATH=$PWD  # set the PYTHONPATH so that the examples dir is found
# Change the graph
tfaip-train examples.tutorial.full --model.graph MLP --model.graph.nodes 200 100 50 --model.graph.activation relu
tfaip-train examples.tutorial.full --model.graph MLP --model.graph.nodes 200 100 50 --model.graph.activation tanh
tfaip-train examples.tutorial.full --model.graph CNN --model.graph.filters 40 20 --model.graph.dense 100
# Change the optimizer
tfaip-train examples.tutorial.full --trainer.optimizer RMSprop --trainer.optimizer.beta1 0.01 --trainer.optimizer.clip_global_norm 1
# ...
```

A trained model can then easily be integrated in a workflow to predict provided `data`:
```python
predictor = TutorialScenario.create_predictor("PATH_TO_TRAINED_MODEL", PredictorParams())
for sample in predictor.predict(data):
    print(sample.outputs)
```

In practice, _tfaip_ follows the rules of object orientation, i.e., the code for a scenario (e.g., image-classification (MNIST), text recognition, NLP, etc.) is organized by implementing classes.
By default, each [`Scenario`](https://tfaip.readthedocs.io/en/latest/doc.scenario.html) must implement [`Model`](https://tfaip.readthedocs.io/en/latest/doc.model.html), and [`Data`](https://tfaip.readthedocs.io/en/latest/doc.data.html).
See [here](examples/tutorial/full) for the complete code to run the upper example for MNIST and see [here](examples/tutorial/min) for the minimal setup.


## Setup

To setup _tfaip_ create a virtual Python (at least 3.7) environment and install the `tfaip` pip package: `pip install tfaip`:
```bash
virtualenv -p python3 venv
source venv/bin/activate
pip install tfaip
pip install tfaip[devel]  # to install additional development/test requirements
```
Have a look at the [wiki](https://tfaip.readthedocs.io/en/latest/doc.installation.html) for further setup instructions.

## Run the Tutorial

After the setup succeeded, launch a training of the tutorial which is an implementation of the common MNIST scenario:
```bash
export PYTHONPATH=$PWD  # set the PYTHONPATH so that the examples dir is found
tfaip-train examples.tutorial.full
# If you have a GPU, select it by specifying its ID
tfaip-train examples.tutorial.full --device.gpus 0
```

## Next Steps

Start reading the [Minimum Tutorial](examples/tutorial/min), optionally have a look at the [Full Tutorial](examples/tutorial/full) to see more features.
The [docs](https://tfaip.readthedocs.io/en/latest) provides a full description of `tfaip`.

To set up a _new custom scenario_, copy the [general template](examples/template/general) and implement the abstract methods.
Consider renaming the classes!
Launch the training by providing the path or package-name of the new scenario which _must_ be located in the `PYTHONPATH`!

## Features of _tfaip_

_tfaip_ provides different features which allow designing generic scenarios with maximum flexibility and high performance.

### Code design

* _Fully Object-Oriented_: Implement classes and abstract functions or overwrite any function to extend, adapt, or modify its default functionality.
* _Typing support_: _tfaip_ is fully typed with simplifies working with an IDE (e.g., use PyCharm!).
* Using pythons `dataclasses` module to set up parameters which are automatically converted to parameters of the command line by our [`paiargparse`](https://github.com/Planet-AI-GmbH/paiargparse) package.

### Data-Pipeline
Every scenario requires the setup of a data-pipeline to read and transform data.
*tfaip* offers to easily implement and modify even complex pipelines by defining multiple `DataProcessors` which usually implement a small operation to map an input sample to an output sample.
E.g., one `DataProcessor` loads the data (`input=filename`, `output=image`), another one applies normalization rules, again another one applies data augmentation, etc.
The **great advantage** of this setup is that the data processors run in Python and can automatically be parallelized by *tfaip* for speed up by setting `run_parallel=True`.

### Deep-Learning-Features

Since _tfaip_ is based on Tensorflow the full API are available for designing models, graphs, and even data pipelines.
Furthermore, *tfaip* supports additional common techniques for improving the performance of a Deep-Learning model out of the box:

* Warm-starting (i.e., loading a pretrained model)
* EMA-weights
* Early-Stopping
* Weight-Decay
* various optimizers and learning-rate schedules

## Contributing

We highly encourage users to contribute own scenarios and improvements of _tfaip_.
Please read the [contribution guidelines](https://tfaip.readthedocs.io/en/latest/doc.development.html).

## Benchmarks

All timings were obtained on a Intel Core i7, 10th Gen CPU.

### MNIST

The following Table compares the MNIST Tutorial of Keras to the [Minimum Tutorial](examples/tutorial/min).
The keras code was adopted to use the same network architecture and hyperparemter settings (batch size of 16, 10 epochs of training).

Code | Time Per Epoch | Train Acc | Val Acc | Best Val Acc
:---- | --------------: | ---------: | -------: | ------------: 
Keras |  16 s | 99.65% | 98.24% | 98.60% 
_tfaip_ | 18 s |  99.76% | 98.66% | 98.66% 

_tfaip_ and Keras result in comparable accuracies, as to be expected since the actual code for training the graph is fundamentally identical.
_tfaip_ is however a bit slower due some overhead in the input pipeline and additional functionality (e.g., benchmarks, or automatic tracking of the best model).
This overhead is negligible for almost any real-world scenario because due to a clearly larger network architecture, the computation times for inference and backpropagation become the bottleneck. 

### Data Pipeline

Integrating pure-python operations (e.g., numpy) into a `tf.data.Dataset `to apply high-level preprocessing is slow by default since [tf.data.Dataset.map](https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map) in cooperation with [tf.py_function](https://www.tensorflow.org/api_docs/python/tf/py_function) does not run in parallel and is therefore blocked by Python's GIL.
_tfaip_ curcumvents this issue by providing an (optional) parallelizable input pipeline.
The following table shows the time in seconds for two different tasks:

* PYTHON: applying some pure python functions on the data
* NUMPY: applying several numpy operations on the data


|         Mode        |     Task     |     Threads 1      |     Threads 2      |     Threads 4      |     Threads 6      |
|:---------------------|:--------------|--------------------:|--------------------:|--------------------:|--------------------:|
| tf.py_function |    PYTHON    | 23.47| 22.78 | 24.38  | 25.76  |
|     _tfaip_    |    PYTHON    | 26.68| 14.48 |  8.11  | 8.13  |
| tf.py_function |    NUMPY     | 104.10 | 82.78  | 76.33  | 77.56  |
|     _tfaip_    |    NUMPY     | 97.07  | 56.93  | 43.78 | 42.73  |

The PYTHON task clearly shows that `tf.data.Dataset.map` is not able to utilize multiple threads.
The speed-up in the NUMPY tasks occurs possibly due to paralization in the numpy API to C.

---
title: '_tfaip_ - a Generic and Powerful Research Framework for Deep Learning based on Tensorflow'
tags:
    - Python
    - Deep Learning
    - Tensorflow
    - Keras
    - Research
    - High-Level Framework
    - Generic
authors:
    - name: Christoph Wick
      orcid: 0000-0003-3958-6240
      affiliation: "1"
    - name: Benjamin Kühn
      affiliation: "1"
    - name: Gundram Leifert
      affiliation: "1"
    - name: Konrad Sperfeld
      affiliation: "2"
    - name: Tobias Strauß
      affiliation: "1"
    - name: Jochen Zöllner
      orcid: 0000-0002-3889-6629
      affiliation: "1,2"
    - name: Tobias Grüning
      orcid: 0000-0003-0031-4942
      affiliation: "1"
affiliations:
    - name: Planet AI GmbH, Warnowufer 60, 18059 Rostock, Germany
      index: 1
    - name: Institute of Mathematics, University of Rostock, 18051 Rostock, Germany
      index: 2
date: 29 March 2021
bibliography: paper.bib
---

# Summary

_tfaip_ is a Python-based research framework for developing, structuring, and deploying Deep Learning projects powered by Tensorflow [@tensorflow2015-whitepaper] and is intended for scientists of universities or organizations who research, develop, and optionally deploy Deep Learning models.
_tfaip_ enables both simple and complex implementation scenarios, such as image classification, object detection, text recognition, natural language processing, or speech recognition.
Each scenario is highly configurable by parameters that can directly be modified by the command line or the API.

# Statement of Need

The implementation of a scenario during research and development typically comprises several tasks, for example setting up the graph (e.g., the network architecture), the training (e.g., the optimizer or learning rate schedule), and the data pipeline (e.g., the data sources).
In most current frameworks such as Tensorflow or Keras (see next Section), the implementation is realized by one or several (unstructured) files in Python which can become difficult to maintain and read for bigger scenarios.
Furthermore, comparable frameworks also only provide basic functionality which is why recurrent obstacles during research and development are typically redundant for each scenario.
_tfaip_ resolves the following different aspects in an elegant and robust way.

During research, a scenario is usually configured via a command line interface (CLI) which usually must be implemented by the user itself.
_tfaip_ already provides a powerful CLI which is dynamically created upon runtime by parsing a nested dataclass hierarchy.
To add a new parameter or even a new set of sub parameters, a user simply has to add a new field to the respective dataclass.
Compared to other approached where all possible parameters are simultaneously available in the CLI, the dynamic approach of _tfaip_ only shows the available arguments.
This prevents users from making mistakes by setting parameters without effect for the current configuration (e.g., setting the factor for an Adam optimizer even though RMSprop was selected).
The default CLI of _tfaip_ provides commands to adapt various hyper-parameters such as the learning rate and its schedule, the optimizer, logging, debugging, profiling, or early stopping.
Furthermore, each component of the scenario can itself be fully customized which allows, for example, to dynamically configure the network architecture, e.g., by inserting layers or changing their parameters.
This feature helps researchers to set up various experiments for example to optimize hyper-parameters or test novel ideas.

In comparison to other frameworks such as Tensorflow, _tfaip_ requires users to implement their scenarios in an object-oriented way and encourages them to annotate their code with types.
This is particularly advantageous for larger, more complex scenarios where the codebase grows accordingly. Using _tfaip_ leads to a clean, structured, modularized, and readable code limiting poor coding styles and facilitating maintenance.
In practice, each scenario is created by implementing predefined interfaces (e.g., loss-function or the graph construction).

During research and development, a tedious step is data preparation which often comprises the conversion of data into the format required by the framework.
The Tensorflow-backed _tfaip_ allows integrating Python code in the data pipeline which is however not run (truly) in parallel by multiple processes and results quite often in a bottleneck.
To speed-up Tensorflow, a user has to transform Python into Tensorflow operations which is laborious, and sometimes even impossible, and complicates debugging.
_tfaip_ tackles this issue by providing a sophisticated pipeline setup based on "data processors" which apply simple transformation operations in pure Python code and are automatically executed in parallel.

Another important step which is simplified by _tfaip_ is the deployment of a scenario.
Other frameworks such as plain Tensorflow or Keras allow to easily load a trained model for prediction which does however not include data processing.
The prediction API of _tfaip_ instead automatically applies pre-processing, infer the trained model, and optionally transform the output by a post-processing pipeline in one step.
The information about the pipeline-construction is embedded within the model which enables to store and load models with a different data pipeline even for the same scenario.
This is handy if, for example, certain pre-processing steps are not required for one specific model or other inputs are expected.

Finally, _tfaip_ will automatically log the training process using the Tensorboard and provides utility scripts to resume a crashed or stopped training, or to set up an array of training configurations via an Excel sheet.

# State of the Field

Efficient research in the area of Deep Learning requires the integration of highly sophisticated Open-Source frameworks such as Tensorflow [@tensorflow2015-whitepaper], PyTorch [@pytorch_2019], Caffe [@jia2014caffe], CNTK [@cntk2016], or Trax [@trax2021].
These frameworks provide efficient tools to freely design Deep Learning scenario of any size and complexity.
However, as the number of lines of code grows, a project has to be structured into meaningful components to be maintainable.
These components are almost identical for each Deep Learning scenario: there are modules for the graph, the model, the data, the training, the validation, and the prediction (i.e., the application of a trained model).
Furthermore, support for dynamic parameter adaption for instance via the command line is desirable for efficient research.
Therefore, to obtain a clean code base, it is highly desirable to only implement abstract templates that already set up the interaction among the modules by providing basic functionality that is required in any use-case.
_tfaip_ which is an extension to Tensorflow solves this and thus helps developers to efficiently handle and maintain small but also large-scale projects in research environments.

Recently, several AutoML approaches have emerged, e.g., by Google Cloud or Microsoft Azure.
AutoML targets developers with limited machine learning expertise and enables them to train their own models by automating processes like network construction, feature engineering, or hyperparameter tuning. 
In contrast, _tfaip_ targets researchers with expertise in deep learning who actually design new network architectures, data processing pipelines, and setup training, but with only limited experience in or capacity for software engineering.
_tfaip_ helps to structure and maintain the code bases, and hereby also solves some recurrent problems that will likely occur during development (see next Section).

# _tfaip_ Functionality

In the following, we highlight the main functionality of _tfaip_.

A basic concept of _tfaip_ is to split parameters and their actual object whereby the parameters are used to instantiate the corresponding object.
This allows to build a hierarchical parameter tree where each node can be replaced with other parameters.
Each parameter and also the replacements can be defined via the command line which enables to dynamically adapt, for example, even complete graphs of a model.
In a research environment this simplifies hyper-parameter search and the setup of experiments.

Class-templates define the logical structure of a real-world scenario.
Each scenario requires the definition of a model and a data class whereby basic functionality such as training, exporting, or loading of a model using the data is already provided.

To set up the data pipeline, a data generator and a list of data processors have to be implemented that define how raw data flows.
The prepared data is then automatically fed into the neural network. Since each data processor is written in pure Python, it is simple to set up and debug the data processing pipeline.
To speed up the data sample generation the pipeline can automatically be run in parallel.

The model comprises information about the loss, metrics, and the graph of the scenario.
Its template hereby requires the user to implement methods, the superordinate modules are then connected automatically.

_tfaip_ tracks the training and validation process by utilizing the Tensorboard provided by Tensorflow.
The Tensorboard can be extended by custom data, for example by plotting Precision-Recall (PR) curves or rendering images.

The prediction module loads a trained scenario so that it can easily be applied on new data during deployment.
Since a model stores its pre- and post-processing pipeline no additional data handling has to be performed.

Because each scenario follows a predefined setup, shared research code is clearer and consequently can be easier reviewed, extended, or applied.
For example, this modularity simplifies the process if several users are working together on the same scenario.

An important feature of _tfaip_ is the consistent use of Python's typing module including type checks.
This leads to clean understandable code and fewer errors during development.
Furthermore, this enables IDEs such as PyCharm [@pycharm2021] to perform autocompletion.

In some rare cases, the highly generic API of _tfaip_ might not be sufficient.
To tackle this, each scenario can optionally customize any functionality by implementing the base classes, for example the trainer or the data pipeline.

# _tfaip_ Documentation and Tutorials

To help new users to become familiar with _tfaip_, a comprehensive documentation, several tutorials, and example scenarios with real-world use-cases are available.
First, _tfaip_ provides two tutorials that solve the classification of MNIST:
the _minimal_ scenario shows the minimal implementation that is required to implement the common MNIST-tutorial in _tfaip_,
the _full_ scenario implements the same use-case and highlights different advanced features of _tfaip_.

Some more examples are provided by transferring official Tensorflow tutorials in the _tfaip_ framework.
These scenarios show the power of the framework as complex scenarios are logically split into meaningful parts and components.
The examples comprise Automatic Text Recognition (ATR) of single text line images, Image Classification, and Fine Tuning of a BERT.

Finally, templates for two basic scenarios allow setting up a new scenario by copying basic code and modifying it afterwards.
All required files and classes are already set up, solely the abstract methods have to be implemented and classes should be renamed.

# Usage of _tfaip_ in Research

Diverse active research projects are already based on _tfaip_.
@zoellner2021 integrated _tfaip_ to solve Natural Language Processing (NLP) problems.
Since its 2.0 release, the open-source ATR engine Calamari by @wick_calamari_2020 is based on _tfaip_.
Our research, for example a recent publication on ATR using Transformers, uses _tfaip_ [@wick_bidirectional_2021].


# Acknowledgments

The authors would like to thank the open-source community, especially the developers and maintainers of Python, Tensorflow, and Numpy, since these packages empower _tfaip_.

This work was partially funded by the European Social Fund (ESF) and the Ministry of Education, Science and Culture of Mecklenburg-Western Pomerania (Germany) within the project NEISS under grant no ESF/14-BM-A55-0006/19.

# References
# TFAIP Examples

This directory comprises some tutorials of few use-cases based on _tfaip_.
The use-cases provide only basic implementations which must be clearly extended for an actual real-world application.
The examples orientate at the [tutorials of Tensorflow](https://www.tensorflow.org/tutorials/).

## Requirements

The examples have additional requirements that can be installed by running:

```shell
pip install -re examples/requirements.txt
```

## List of Scenarios
* [MNIST:](tutorial) MNIST scenario designed as tutorial to show the [basic](tutorial/min) or [extended](tutorial/full) features of _tfaip_.
* [Template:](template) Templates to set up new scenarios
* [Image Classification:](imageclassification) Image classification of flowers of five classes (corresponding Tensorflow tutorial is available [here](https://www.tensorflow.org/tutorials/images/classification))
* [ATR:](atr) Line-based Automatic Text Recognition
* [Text/Fine-tuning BERT:](text/finetuningbert) Fine-tuning of BERT on the task that decides whether two sentences are semantically equivalent.
# Minimal Tutorial

Welcome to the minimal tutorial which shows by the example of MNIST how _tfaip_ is structured and what a custom scenario must implement.

## Scenario

The central class is the `TutorialScenario` which is (optionally) parametrized by `TutorialScenarioParams`. Here no additional parameters are required for the overall scenario.
The `TutorialScenario` glues together the `Data`-Pipeline, the training data origin, and the `Model` with is `Graph`.

## Data

The [`TutorialData`](data.py) defines the input and target shapes of the `Scenario`.
In the case of MNIST, the input is an image labelled `img` with a shape of `28x28` and a `dtype` of `uint8`.
The targets (`gt`) are a simple scalar (use `shape=[1]`) with a `dtype` of `uint8`.

```python
class TutorialData(DataBase[DataBaseParams]):
    def _input_layer_specs(self):
        return {'img': tf.TensorSpec(shape=(28, 28), dtype='uint8')}

    def _target_layer_specs(self):
        return {'gt': tf.TensorSpec(shape=[1], dtype='uint8')}
```

## Training Data-Generation

The next step is to load the training data which is done in the `TutorialDataGenerator` (note that the corresponding parameter class will actually instantiate the Generator).
This first loads the `keras.dataset` and selects the desired partition depending on the pipeline mode.
A general `DataGenerator` must overwrite `__len__` to return the number of samples in the dataset and the `generate` to provide the samples.
Here, `generate()` just returns the loaded (and converted) samples.

```python
class TutorialDataGenerator(DataGenerator[TutorialDataGeneratorParams]):
    def __init__(self, mode: PipelineMode, params: 'TutorialDataGeneratorParams'):
        super().__init__(mode, params)
        dataset = getattr(keras.datasets, params.dataset)
        train, test = dataset.load_data()
        data = train if mode == PipelineMode.TRAINING else test
        self.data = to_samples(data)

    def __len__(self):
        return len(self.data)

    def generate(self) -> Iterable[Sample]:
        return self.data
```

Afterwards, parameters that define the `TrainerPipelineParamsBase` must be defined.
This parameter set defines how to create the training data and optionally the validation data.
Here, only one (`train_val`) `TutorialDataGeneratorParams`-set is used for training and validation.
Split this into two, if there are different parameters for training and validation.
```python
@pai_dataclass
@dataclass
class TutorialTrainerGeneratorParams(TrainerPipelineParamsBase[TutorialDataGeneratorParams, TutorialDataGeneratorParams]):
    train_val: TutorialDataGeneratorParams = field(default_factory=TutorialDataGeneratorParams, metadata=pai_meta(mode='flat'))

    def train_gen(self) -> TutorialDataGeneratorParams:
        return self.train_val

    def val_gen(self) -> Optional[TutorialDataGeneratorParams]:
        return self.train_val
```

Now, the training pipeline is set-up and the model and its graph can be defined.

## Model and Graph

The `TutorialModel` and its `TutorialGraph` are described by a corresponding `TutorialModelParams` class:
```python
@pai_dataclass
@dataclass
class TutorialModelParams(ModelBaseParams):
    n_classes: int = field(default=10, metadata=pai_meta(
        help="The number of classes (depends on the selected dataset)"))
```

Here, only the number of classes must be specified (which could be derived from data).

The `TutorialGraph` inherits `GraphBase` and uses the `TutorialModelParams`-structure.
A `GraphBase` is a `keras.Layer` so, as recommended by Tensorflow, first create the layers in the `__init__` function,
then connect them in the `call` method.
Here, a CNN with two conv and pool layers and two dense layers (FF) is created.
The `input` of the layer are of the shape that are defined in [`Data`](#data) ([above](#data)), therefore a dict with one entry of `'img'`.
The `call` function must also return a `dict` of tensors, its keys are later used to access the outputs (e.g. in the loss or the metric).

```python
class TutorialGraph(GraphBase[TutorialModelParams]):
    def __init__(self, params: TutorialModelParams, name='conv', **kwargs):
        super(TutorialGraph, self).__init__(params, name=name, **kwargs)
        self.conv1 = Conv2D(kernel_size=(2, 2), filters=16, strides=(1, 1), padding='same', name='conv1')
        self.pool1 = MaxPool2D(pool_size=(2, 2), strides=(2, 2), name='pool1')
        self.conv2 = Conv2D(kernel_size=(2, 2), filters=32, strides=(1, 1), padding='same', name='conv2')
        self.pool2 = MaxPool2D(pool_size=(2, 2), strides=(2, 2), name='pool2')
        self.flatten = keras.layers.Flatten()
        self.ff = FF(out_dimension=128, name='f_ff', activation='relu')
        self.logits = FF(out_dimension=self._params.n_classes, activation=None, name='classify')

    def call(self, inputs, **kwargs):
        rescaled_img = K.expand_dims(K.cast(inputs['img'], dtype='float32') / 255, -1)
        conv_out = self.pool2(self.conv2(self.pool1(self.conv1(rescaled_img))))
        logits = self.logits(self.ff(self.flatten(conv_out)))
        pred = K.softmax(logits, axis=-1)
        cls = K.argmax(pred, axis=-1)
        return {'pred': pred, 'logits': logits, 'class': cls}
```

The corresponding `TutorialModel` must define how to create the graph, here by just constructing a `TutorialGraph`.
Next, define the `_loss`: the `inputs_targets` are a joined dict of both the inputs and targets coming from the [`TutorialData`](#data).
The `outputs` are the output-dict of the previously defined `TutorialGraph`.
The loss is again an "output" of the graph and must therefore be wrapped in a `keras.Layer`, here, a `sparse_categorical_crossentropy` is wrapped within a `Lambda`-layer.
The metric is computed by selecting the `gt` (targets) and `class` (outputs) and passing them to a `keras.metrics.Accuracy()`.
There are different, also more complex methods to define metrics, see the [full tutorial](../full) for an example.
Both the loss and metric must return a dict. Its keys will be used for displaying information in the log or the TensorBoard.

Finally, optionally define how to determine the best model (`_best_logging_settings_`), here by selecting the model with the `"max"` `"acc"`, and
how to print useful human-readable information during validation.
The `_print_evaluate`-function receives a single sample (unbatched) and displays (here) the target, prediction and if the prediction was correct.

```python
class TutorialModel(ModelBase[TutorialModelParams]):
    def create_graph(self, params: TutorialModelParams) -> 'GraphBase':
        return TutorialGraph(params)

    def _loss(self, inputs_targets, outputs) -> Dict[str, AnyTensor]:
        return {'loss': tf.keras.layers.Lambda(
            lambda x: tf.keras.metrics.sparse_categorical_crossentropy(*x, from_logits=True), name='loss')(
            (inputs_targets['gt'], outputs['logits']))}

    def _metric(self):
        return {'acc': MetricDefinition("gt", "class", keras.metrics.Accuracy())}

    def _best_logging_settings(self):
        return "max", "acc"

    def _print_evaluate(self, inputs, outputs: Dict[str, AnyNumpy], targets: Dict[str, AnyNumpy], data, print_fn=print):
        correct = outputs['class'] == targets['gt']
        print_fn(f"PRED/GT: {outputs['class']}{'==' if correct else '!='}{targets['gt']} (p = {outputs['pred'][outputs['class']]})")
```

## Launch the Training
Training can be started by calling
```bash
tfaip-train tutorial.min
```

## Further reading

After having finished this Tutorial, have a look at the [full tutorial](../full) or the [docs](https://tfaip.readthedocs.io/).
# Full Tutorial

This tutorial shows how to use and implement optional by handy features of *tfaip*.
Go to the [minimal tutorial](../min) to see a scenario that only implements the required classes and functions.

This tutorial sets up training on the MNIST data which can then be used to predict digits of image files.

## Overview

The following features are covered by this tutorial

* setting up of a DataPipeline using DataProcessors, see [here](data).
* setting up of different data generation for [training](data/training_data_generation.py) and [prediction](data/prediction_data_generation.py), see also [here](data).
* selection and configuration of [different dynamic graphs](#dynamic-graphs)
* Writing image files to the tensorboard, see [here](model.py)
* setting up a Predictor that can vote the predictions of multiple individual models, see [here](predictor.py)
* Evaluator


## Dynamic Graphs

Dynamic graphs allow to change and setup layers with parameters that can be set in the command line.

* First, setup a static [Graph](graphs/tutorialgraph.py) which will handle the creation of the dynamic layers.
  For MNIST, this graph also adds the final output as additional layer since it is obligatory.
  Furthermore, the data is normalized and reshaped.
* Next, create a [base class and base params](graphs/backend.py) which is derived from `keras.layers.Layer`
  and must layer be implemented by each variant.
  Add an abstract method to the parameters to define how to create the layer.
  Here (`cls()`), only the class type is returned while assuming that the first and only argument of the `__init__` is the parameter.
  Optionally define a generic `TypeVar` for the parameters that can be used to define the actual parameter type in the actual implemented layer.
* Now, implement the base class and base parameters.
  In the tutorial, a [CNN](graphs/cnn.py) and [MLP](graphs/mlp.py) setup is provided.
* Finally, add a parameter to select the layers to the base params, here called `graph` in the `ModelParams`.
  Optionally, set the `choices` flag of `pai_meta` to provide the list of available parameters that can be selected.
  The static [Graph](graphs/tutorialgraph.py) calls the abstract `cls()` method to retrieve the actual implementation
  and instantiates it.
# Fine Tuning BERT Tutorial

The Fine Tuning BERT Tutorial maps the corresponding [Tensorflow Tutorial](https://www.tensorflow.org/official_models/fine_tuning_bert) to _tfaip_.
Only the training of the BERT is covered in this Tutorial.
Application, i.e. prediction is missing.
Note that instead of the large BERT, this Tutorial used Albert from [Hugging Face](https://huggingface.co/albert-base-v2) as pretrained model.
The task is to decide whether to sentences are semantically equivalent., e.g. the two sentences

> The identical rovers will act as robotic geologists , searching for evidence of past water .

and

> The rovers act as robotic geologists , moving on six wheels .

are not equivalent. 


## Run
To run the training of this scenario execute (in the cloned dir)
```bash
export PYTHONPATH=$PWD  # required so that the scenario is detected

# Training
tfaip-train examples.text.finetuningbert --trainer.output_dir ftbert_model

# TensorBoard can be used to view the training process
tensorboard --logdir .
# Open http://localhost:6006
tfaip-monitor/pythonProject
```

## References
* Official [Tensorflow Tutorial](https://www.tensorflow.org/official_models/fine_tuning_bert)
* A model zoo of diverse pretrained networks is available at [Hugging Face](https://huggingface.co)
# Template

This scenario can be used as blueprint for setting up a new `ListFileScenario` by implementing all recommended functions.

## Test

The `test.py` file comprises a few standard tests that can be run on a newly implemented scenario.
Especially `test_data` is useful to debug the input pipeline.
# Template

This scenario can be used as blueprint for setting up a new scenario by implementing all recommended functions.

## Test

The `test.py` file comprises a few standard tests that can be run on a newly implemented scenario.
Especially `test_data` is useful to debug the input pipeline.
# Image Classification Scenario

The Image Classification Tutorial maps the corresponding [Tensorflow Tutorial](https://www.tensorflow.org/tutorials/images/classification) to _tfaip_.
Improvements such as data augmentation and dropout are omitted in this tutorial and can be integrated as an exercise.

## Run
To run the training of this scenario execute (in the cloned dir)
```bash
export PYTHONPATH=$PWD  # required so that the scenario is detected

# Training
tfaip-train examples.imageclassification --trainer.output_dir ic_model
tfaip-train examples.imageclassification --trainer.output_dir ic_model --device.gpus 0  # to run training on the first GPU, if available
tfaip-train examples.imageclassification --model.conv_filters 30 50 60 --model.dense 200 200 --trainer.output_dir ic_model --device.gpus 0  # try a different (larger) model

# TensorBoard can be used to view the training process
tensorboard --logdir .
# Open http://localhost:6006

# Prediction
tfaip-predict --export_dir ic_model/best --data.image_files examples/imageclassification/examples/592px-Red_sunflower.jpg 
# Possible Output >>> This image most likely belongs to sunflowers with a 0.83 percent confidence.


```

## References
* Official [Tensorflow Tutorial](https://www.tensorflow.org/tutorials/images/classification)
* PLANET AI GmbH offers an [Intelligent Image Analysis](https://planet-ai.de/applications/image-analysis/) which is able to localize and identify many visual categories in images and videos.
  This skill is often utilized by our engines to preprocess images or videos to segment certain objects, but can also be used as a stand-alone solution.
# ATR Scenario

The ATR-Scenario is an example showing how to implement a line-based ATR-engine.
It provides a CNN/LSTM-network architecture which is trained with the CTC-algorithm.
This tutorial shows only the fundamentals and does not include required algorithms for document analysis in general
which are part of a real OCR/Document-Analysis engine.

## Run
To run the training of this scenario execute (in the cloned dir)
```bash
export PYTHONPATH=$PWD  # required so that the scenario is detected

# Training
tfaip-train examples.atr --trainer.output_dir atr_model
tfaip-train examples.atr --trainer.output_dir atr_model --device.gpus 0  # to run training on the first GPU, if available

# Validation (of the best model)
tfaip-lav --export_dir atr_model/best --data.image_files examples/atr/workingdir/uw3_50lines/test/*.png

# Prediction
tfaip-predict --export_dir atr_model/best --data.image_files examples/atr/workingdir/uw3_50lines/test/*.png
```

Note, the prediction will only print the raw output of the network.

## Data
The [working dir](workingdir) provides some example lines of the UW3 dataset which are loaded by default

## References
* PLANET AI GmbH offers an intelligent [Document Analysis Suite (IDA)](https://planet-ai.de/applications/document-analysis/) which is able to read and even understand a broad spectrum of documents from ancient hand-written documents to modern machine-generated ones.# XLSX Experimenter

This tool is designed for scheduling/logging of experiments by using an xlsx file.
Pass an xlsx file to parse `--xlsx`, the desired `--gpus` for scheduling, 
optionally the python executable of the target virtual env via `--python`.
Set `--update` to write the metrics and losses of the lastly calculated epoch 
of each scheduled experiment to a new sheet in the xlsx file.
Launches in the "Working"-Dir.

**Notes:**
* This tool currently requires python >=3.7

Examples:
* Run on gpus 1 and 3, `tfaip-experimenter --xlsx demo.xlsx --gpus 1 3`
* Run on cpu 0 `tfaip-experimenter --xlsx demo.xlsx --cpus 0`

## XLSX file

First row:
* USER (optional): Arbitrary user parameters
* RESULT (optional): Arbitrary result of the experiment (must be manually inserted)
* SKIP (required): Skip this experiment if the xlsx file is parsed
* SCRIPT (required): The script to run, typically `tfaip-train` (expected to be in the scripts dir). 
  If your experiments where stopped for whatever reason, `SKIP` all finished trainings and replace 
  `tfaip-train` with `tfaip-resume-training` for all already started but not finished trainings.   
* ID (required): ID of the run (best practice: simply use a running integer)
* CLEANUP (required): Cleanup the output directory (not implemented!)
* PARAMS (requied): Tag to mark the begin of parameters passed to the "SCRIPT"

Second row (for PARAMS), defines the paramter groups:
* if `default`: no group
* if `NAME` (e.g. `optimizer_params`): pass the following parameters to the `--optimizer_params` flag
* if `EMPTY`: use the previous group definition

Third row (PARAM labels, setup):
* specify the parameter name, e.g. `epochs` (without dashes)
* optionally pass additional flags (e.g. `post_distortions:split:empty`):
  * `:split` split the value of this parameter into individual args
  * `:empty` allow this parameter to be empty when added to the list (note, by default, empty `NULL` values are omitted)

Remarks: 
* Don't use formulas in the xlsx file. 
