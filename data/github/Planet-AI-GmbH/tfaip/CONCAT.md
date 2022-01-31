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
Model
=====

The model glues together several parts that define the setup of the neural network, e.g. the :ref:`loss<doc.model:Loss>` or the :ref:`metrics<doc.model:Metric>`.

The implementation of the model requires to override the base class :ref:`ModelBase<tfaip.model:ModelBase>`
and its parameters :ref:`ModelBaseParams<tfaip.model:ModelBaseParams>` (see the following example of the tutorial):

.. code-block:: python

    @pai_dataclass
    @dataclass
    class TutorialModelParams(ModelBaseParams):
        n_classes: int = field(default=10, metadata=pai_meta(
            help="The number of classes (depends on the selected dataset)"
        ))

        @staticmethod
        def cls():
            return TutorialModel

        def graph_cls(self):
            from examples.tutorial.min.graphs import TutorialGraph
            return TutorialGraph

    class TutorialModel(ModelBase[TutorialModelParams]):
        pass

Parameter Overrides
-------------------

The implementation of the ``ModelBaseParams`` require to override ``cls()`` and ``graph_cls`` to return the class type of the actual model and :ref:`graph<doc.graph:Graph>`.

Loss
----

The loss function defines the optimization target of the model.
There are two ways to define a loss: loss using a ``keras.losses.Loss``, or a loss using a Tensor as output.
Multiple losses can be :ref:`weighted<doc.model:loss weight>`.
The output-values of each loss (and the weighted loss) will be displayed in the console and in the :ref:`Tensorboard<doc.model:tensorboard>`.

Overwrite ``_loss`` and return a dictionary of losses where the key is the (display) label of the metric and the value is the Tensor-valued loss.

To use a ``keras.losses.Loss``, instantiate the loss in the ``__init__``-function and call it in ``_loss``.
Alternatively, return any scalar-valued Tensor.

.. code-block:: python

    def __init(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.scc_loss = keras.losses.SparseCategoricalCrossentropy(from_logits=True, name='keras_loss')

    def _loss(self, inputs, targets, outputs) -> Dict[str, AnyTensor]:
         return {
            'keras_loss': self.scc_loss(targets['gt'], outputs['logits']),  # either call a keras.Loss
            'raw_loss': tf.keras.losses.sparse_categorical_crossentropy(targets['gt'], outputs['logits'], from_logits=True),  # or add a raw loss
        }

Loss Weight
~~~~~~~~~~~

If multiple losses are defined, the ``_loss_weights`` function can be implemented to return weights for the losses.
Here both upper losses are weighted with a factor of 0.5.
If not implemented, each loss is weighted by a factor of 1.

.. code-block:: python

    def _loss_weights(self) -> Optional[Dict[str, float]]:
        return {'keras_loss': 0.5, 'extended_loss': 0.5}

Metric
------

Similar to the loss, a model defines its metrics.
The output-values of each metric will be displayed in the console and in the :ref:`Tensorboard<doc.model:Tensorboard>`.
All metrics are computed on both the training and validation data, except the :ref:`pure Python<doc.model:pure-python metric>` one which is solely computed on the validation set.

Overwrite ``_metric`` and return a list of called ``keras.metric.Metric``.
The ``name`` of the metric is used for display.

.. code-block:: python

    def __init(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.acc_metric = keras.metrics.Accuracy(name='acc')

    def _metric(self, inputs, targets, outputs):
        return [self.acc_metric(targets['gt'], outputs['class'])]

Custom metrics must implement ``keras.metrics.Metric``.
It is also possible to compute the actual value of the metric as Tensor beforehand and wrap it with a ``keras.metrics.Mean``.

Pure-Python Metric
~~~~~~~~~~~~~~~~~~

Pure python metrics are not defined with the Model but instead in the :ref:`Evaluator<doc.scenario:evaluator>`.
They provide a maximum of flexibility since they are computed during :ref:`load and validate<doc.training:extended>` in pure Python.

Logging the best model
----------------------

During :ref:`training<doc.training:training>` the best model will be tracked and automatically exported as "best".
The best model is determined by a models ``_best_logging_settings`` which is by default the minimum loss since every model provides this information.
If you want to track the best model for example by a metric, overwrite this function.
For instance, if a model defines a :ref:`metric<doc.model:Metric>` ``"acc"``, use

.. code-block:: python

    def _best_logging_settings(self):
        return "max", "acc"

The first return value is either ``"max"`` or ``"min"`` while the second argument is the name of a metric or loss.


Output during validation
------------------------

During validation the first few examples are passed to a ``Model``'s ``_print_evaluate`` function which can be used to display the current state of training in a human-readable form.
For MNIST-training this could be the target class and the prediction probabilities, e.g.:

.. code-block:: python

    def _print_evaluate(self, sample: Sample, data, print_fn=print):
        outputs, targets = sample.outputs, sample.targets
        correct = outputs['class'] == targets['gt']
        print_fn(f"PRED/GT: {outputs['class']}{'==' if correct else '!='}{targets['gt']} (p = {outputs['pred'][outputs['class']]})")

Note that a sample is already un-batched.
This function can also access to the ``data``-class if a mapping (e.g. a codec) must be applied.

Tensorboard
-----------

During training, the output of the loss and metrics on the training and validation sets is automatically to the Tensorboard.
The data is stored in the ``output_dir`` defined during [training](07_training.md).

In some cases, additional :ref:`arbitrary data<doc.model:arbitrary data>` such as images, or raw data e.g. such as :ref:`PR-curves<doc.model:pr-curves>` shall be written to the Tensorboard.

Arbitrary Data
~~~~~~~~~~~~~~

To add arbitrary additional data to the Tensorboard ensure that the layer adding the data inherits ``TFAIPLayerBase`` which provides a method ``add_tensorboard`` which must be called with a ``TensorboardWriter`` and the ``value``.

The following examples shows how to write the output of a conv-layer as image to the Tensorboard.
The ``TensorboardWriter`` will receive the raw numpy data and call the provided ``func`` (here ``handle``) to process the raw data and write it to the tensorboard.

.. code-block:: python

    def handle(name: str, value: np.ndarray, step: int):
        # Create the image data as numpy array
        b, w, h, c = value.shape
        ax_dims = int(np.ceil(np.sqrt(c)))
        out_conv_v = np.zeros([b, w * ax_dims, h * ax_dims, 1])
        for i in range(c):
            x = i % ax_dims
            y = i // ax_dims
            out_conv_v[:, x * w:(x + 1) * w, y * h:(y + 1) * h, 0] = value[:, :, :, i]

        # Write the image (use 'name_for_tb' and step)
        tf.summary.image(name, out_conv_v, step=step)

    class Layers(TFAIPLayerBase):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.conv_layer = Conv2D(40)
            self.conv_mat_tb_writer = TensorboardWriter(func=handle, dtype='float32', name='conv_mat')

        def call(self, inputs, **kwargs):
            conv_out = self.conv_layer(inputs)
            self.add_tensorboard(self.conv_mat_tb_writer, conv_out)
            return conv_out

PR-curves
~~~~~~~~~

If a metric (e.g. the PR-curve) returns binary data (already serialized Tensorboard data) it will be automatically written to the Tensorboard.

Exporting additional graphs
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    def _export_graphs(self,
                       inputs: Dict[str, tf.Tensor],
                       outputs: Dict[str, tf.Tensor],
                       targets: Dict[str, tf.Tensor],
                       ) -> Dict[str, tf.keras.Model]:
        # Override this function
        del targets  # not required in the default implementation
        return {"default": tf.keras.Model(inputs=inputs, outputs=outputs)}

This function defines the graphs to export.
By default, this is the graph defined by all inputs and all outputs.
Override this function to export a different or additional graphs, e.g., if you want to only export the encoder in an encoder/decoder setup.
Return a Dict with ``label`` and ``keras.models.Model`` to export.

Root-Graph-Construction
-----------------------

The root graph can be overwritten to have full flexibility when creating a graph.
In most cases this is optional.

.. code-block:: python

    @staticmethod
    def root_graph_cls() Type['RootGraph']:
        from tfaip.model.graphbase import RootGraph
        return RootGraph
Prediction
==========

After a :ref:`Scenario<doc.scenario:scenario>` was :ref:`trained<doc.training:training>` and `evaluated (LAV)<doc.evaluation:load and validate (LAV)`, a ``Predictor`` is used to apply the resulting model on new data.
Similar to LAV, a prediction can be performed with the command line (``tfaip-predict``) or programmatically.

Each resulting ``Sample`` of the predictor passes the :ref:`post_proc-Pipeline<doc.data:data>`.

Command-Line
------------

.. code-block:: shell

    tfaip-predict --export_dir PATH_TO_EXPORTED_MODEL --data.DATA_GENERATOR_PARAMS --predict.PREDICTOR_PARAMS --pipeline.PIPELINE_PARAMS


``tfaip-predict`` requires at least the path to the :ref:`exported<doc.training:Export of Checkpoints and Saved Models>` model (``--export_dir``).
The ``--data`` flag is used to modify the input data, the ``--prediction`` flag sets up the predictor.
Another optional flag is ``--pipeline`` which allow to specify the ``DataPipelineParams``, e.g., the number of threads or batch size.

Programmatically
----------------

Create a ``Predictor`` by calling ``Scenario.create_predictor(model: str, params: PredictorParams)``.
(Overriding ``Scenario.predictor_cls()`` can be used to customize a ``Predictor``).
The resulting object can be used to ``predict`` given ``DataGeneratorParams``, raw input data ``predict_raw`` which is basically the output of a ``DataGenerator``, to predict a ``DataPipeline`` (``predict_pipeline``), or to predict a ``tf.data.Dataset`` (``predict_database``).
The output of each function is a generator of ``Samples``.

A predicted ``Sample`` holds ``inputs`` and ``outputs``.
Optionally, if available in the data and if the ``predictor_params.include_targets`` flag is set, also the ``targets``.

For example:

.. code-block::python

    predictor = MyScenario.create_predictor("PATH_TO_SAVED_MODEL", PredictorParams())

    # Predict on raw data
    for sample in predictor.predict_raw([Sample(inputs=np.zeros([28, 28]))]):
        print(sample)

    # Predict a data generator
    data = MyDataGenerator()
    for sample in predictor.predict(data):
        print(sample)
tfaip.trainer.scheduler
============================

Weight decay
------------
.. automodule:: tfaip.trainer.scheduler.schedule_weightdecay


Learning Rate
-------------
.. automodule:: tfaip.trainer.scheduler.learningrate_params
.. automodule:: tfaip.trainer.scheduler.learningrate

Constant
~~~~~~~~
.. automodule:: tfaip.trainer.scheduler.constant_params
.. automodule:: tfaip.trainer.scheduler.constant

Cosine Decay
~~~~~~~~~~~~
.. automodule:: tfaip.trainer.scheduler.cosine_decay_params
.. automodule:: tfaip.trainer.scheduler.cosine_decay

Exponential Decay
~~~~~~~~~~~~~~~~~
.. automodule:: tfaip.trainer.scheduler.exponential_decay_params
.. automodule:: tfaip.trainer.scheduler.exponential_decay
tfaip.trainer.callbacks
============================

Early Stopping
--------------

.. automodule:: tfaip.trainer.callbacks.earlystopping.params
.. automodule:: tfaip.trainer.callbacks.earlystopping.callback

BenchmarkCallback
-----------------
.. automodule:: tfaip.trainer.callbacks.benchmark_callback

EMACallback
-----------
.. automodule:: tfaip.trainer.callbacks.ema_callback

ExtractLogsCallback
-------------------
.. automodule:: tfaip.trainer.callbacks.extract_logs

LAVCallback
-----------
.. automodule:: tfaip.trainer.callbacks.lav_callback

LoggerCallback
--------------
.. automodule:: tfaip.trainer.callbacks.logger_callback

ProgbarCallback
---------------
.. automodule:: tfaip.trainer.callbacks.progbar

TensorBoardCallback
-------------------
.. automodule:: tfaip.trainer.callbacks.tensor_board_callback

TensorBoardDataHandlerCallback
------------------------------
.. automodule:: tfaip.trainer.callbacks.tensor_board_data_handler

TensorflowFixCallback
---------------------
.. automodule:: tfaip.trainer.callbacks.tensorflow_fix

TrainParamsLoggerCallback
-------------------------
.. automodule:: tfaip.trainer.callbacks.train_params_logger
Debugging
=========

Sooner or later, there will be a point during development where debugging is essential.
Since |tfaip| uses Tensorflow 2 which allows for eager execution, debugging is drastically simplified as each computation of the graph can be manually traced by a debugger.
Unfortunately, there are some rare circumstances that lead to code that runs in eager but not in graph mode.
Usually, the reason is that operations are used that are only allowed in eager-mode.

This file shows how to efficiently debug the :ref:`data-pipeline<doc.debugging:data-pipeline>`, the :ref:`model<doc.debugging:model>`, its :ref:`graph<doc.debugging:graph>`, :ref:`loss<doc.debugging:loss>`, and :ref:`metrics<doc.debugging:metric>`, and how to :ref:`profile<doc.debugging:profiling>` using the Tensorboard which helps to detect bottlenecks.


Data-Pipeline
-------------

Debugging of the data-pipeline is usually the first step since this helps to verify the data integrity.
While ``tf.data.Datasets`` can not be debugged easily, the |tfaip| pipeline based on :ref:`DataProcessors<doc.data:data>` can easily be debugged.

Do the following:

* Set up a data-:ref:`test<doc.installation:tests>` for your scenario and :ref:`Data<doc.data:data>` class
* Disable all multiprocessing (setting ``run_parallel=False`` in the ``pre_proc`` and ``post_proc`` parameters of the ``DataParams``).
* Create a ``DataPipeline``.
* Call ``with pipeline.generate_input_samples() as samples:`` which will return a Generator of samples which are the (un-batched) input of the ``tf.data.Dataset``.
* Optionally call ``with pipeline.generate_input_batches() as batches`` to access the outputs of the ``tf.data.Dataset.as_numpy_iterator()``. Note that this makes debugging of the pipeline impossible this ``tf.data.Dataset`` is accessed.
  Use this only if you want to verify the batched and padded outputs of the dataset not to debug the data-pipeline itself.

Here is an example for the Tutorial:

.. code-block:: python

    class TestTutorialData(unittest.TestCase):
        def test_data_loading(self):
            trainer_params = TutorialScenario.default_trainer_params()
            data = TutorialData(trainer_params.scenario.data)
            with trainer_params.gen.train_data(data).generate_input_samples(auto_repeat) as samples:
                for sample in samples:
                    print(sample)  # un-batched, but can be debugged

            # or
            with trainer_params.gen.train_data(data).generate_input_batches(auto_repeat) as batches:
                for batch in batches:
                    print(batch)  # batched and prepared (inputs, targets) tuple, that can not be debugged. Use prints.

Note that ``generate_input_samples()`` will run infinitely for the ``train_data`` which is why ``auto_repeat=False`` is set to only generate an epoch of data.

Model
-----

To allow for debugging of the model, enable the eager mode (pass ``--trainer.force_eager True`` during :ref:`training<doc.training:training>`, or ``--lav.run_eagerly True`` during :ref:`LAV<doc.evaluation:load and validate (lav)>`)).
Now, the full computations of the graph can be followed.

Graph
~~~~~

During training, additionally pass ``--scenario.debug_graph_construction``.
This will once evaluate the (prediction) graph and compute the :ref:`loss<doc.debugging:loss>` and :ref:`metrics<doc.debugging:metric>` on real data.
It is recommended to use this flag if any error occurs in the graph during construction.

Loss
~~~~

Losses can be fully debugged in eager mode.

Metric
~~~~~~

Metrics of the model can be fully debugged in eager mode.
Also metrics defined in the :ref:`Evaluator<doc.scenario:evaluator>` can always be debugged since they run in pure Python.

Profiling
---------

Profiling is useful to detect bottlenecks in a scenario that slow down training.
Pass the ``--trainer.profile True`` flag to write the full profile of the training (graph mode required) to the Tensorboard.
Also have a look at the `official documentation <https://www.tensorflow.org/tensorboard/tensorboard_profiling_keras#use_the_tensorflow_profiler_to_profile_model_training_performance>`_.

Optimizing the input pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In many cases, the input pipeline to too slow to generate samples for the model. However, there are several parameters for tweaking:

* First, enable parallel processing of the pipeline by setting ``run_parallel`` to ``True``.
* Increase the number of threads for the pipeline ``--train.num_processes 16``.
* Change the default behaviour for prefetching ``--train.prefech 128``.
* Verify that the size of a sample is as small as possible. Python required to pickle the data for parallelization which can drastically slow down the queue-speed.
  We observed crucial problems if the input data size is in the order of more than 50 MB. Consider changing the data type (e.g. ``uint8`` instead of ``int32``)


Optimizing the model
~~~~~~~~~~~~~~~~~~~~

The standard way to increase the throughput of a model is to increase its batch size if the memory of a GPU is not exceeded: ``--train.batch_size 32``.
Installation
============

|tfaip| is available on `pypi <https://pypi.org/project/tfaip>`_ and can thus be installed with `pip`:

.. code-block:: shell

    pip install tfaip


Prerequisites
-------------

You need to install the following on your system:

* Python 3.7 or later, including
* the python development packages (on Ubuntu ``apt install python3.7 python3.7-dev``, if not available in your distribution, use the `deadsnakes repo <https://launchpad.net/~deadsnakes/+archive/ubuntu/ppa>`_ for ubuntu based distros)
* a gcc compiler and all other build essentials (on Ubuntu ``apt install build-essential``)
* (optional) cuda/cudnn libs for GPU support, see `tensorflow <https://www.tensorflow.org/install/source#tested_build_configurations>`_ for the versions which are required/compatible.

Setup in a Virtual Environment
------------------------------

Setup your venv and install:

``<python>`` must be replaced with ``python3.7`` or later (check ``<python> --version`` before):

.. code-block:: shell

    virtualenv -p <python> PATH_TO_VENV
    source PATH_TO_VENV/bin/activate
    pip install -U pip  # recommended to get the latest version of pip
    pip install tfaip


Possible bugs
~~~~~~~~~~~~~

pycocotools
"""""""""""

Currently as of `numpy` 1.20, pycocotools are falsely compiled, run

.. code-block:: shell

    pip uninstall -y pycocotools
    pip install --no-cache-dir pycocotools

to fix your venv by reinstalling and recompiling the pycocotools.

Tests
-----

|tfaip| provides a set of tests which are listed in the `test` directory using the `unittest` framework.

Run the tests
~~~~~~~~~~~~~

Call ``python -m unittest`` or ``pytest`` to run the tests from the command line.

In PyCharm:  Right-click on ``test`` select ``Run 'Unittests in test'``, or "Shift-Shift", write "run" and choose ``Run 'Unittests in test'``


CI/CD
~~~~~
All tests will automatically run for CI/CD.


Development Setup
-----------------

For development support, clone the code, install the requirements in a fresh virtual environment, and link the |tfaip| source to the virtualenv:

.. code-block:: shell

    git clone https://github.com/Planet-AI-GmbH/tfaip.git  # or git clone git@github.com:Planet-AI-GmbH/tfaip.git
    cd tfaip
    virtualenv -p python3 venv
    source venv/bin/activate
    pip install -U pip  # recommended to get the latest version of pip
    pip install -r requirements.txt
    pip install -r devel_requirements.txt  # Requirements only required for testing and development
    python setup.py develop
Development and Contributions
=============================

We highly encourage users of |tfaip| to contribute scenarios but also welcome improvements and bug fixes of the core classes and the documentation!

Create `github issues <https://github.com/Planet-AI-GmbH/tfaip/issues>`_ to report bugs or if you seek for support.
Open a pull request if you want to contribute.

Pull Request Checklist
----------------------

Before sending pull requests, make sure you do the following:

* Read the present contributing guidelines.
* Read the `Code of Conduct <https://github.com/tensorflow/tensorflow/blob/master/CODE_OF_CONDUCT.md>`_
* Check if your changes are consistent with the guidelines.
* Check if your changes are consistent with the Coding Style.
* Run the unit tests.

Guidelines
----------

* Include unit tests when you contribute new features, as they help to a) prove that your code works correctly, and b) guard against future breaking changes to lower the maintenance cost.
* Bug fixes also generally require unit tests, because the presence of bugs usually indicates insufficient test coverage.
* Keep API compatibility in mind when you change code in |tfaip|
* As every PR requires several CPU time of CI testing, we discourage submitting PRs to fix one typo, one warning, etc. We recommend fixing the same issue at the file level at least (e.g.: fix all typos in a file, fix all compiler warning in a file, etc.)
* Tests should follow the testing best practices guide.


Licence
-------

Include the license at the top of new files:

.. code-block:: python

    # Copyright 2021 The tfaip authors. All Rights Reserved.
    #
    # This file is part of tfaip.
    #
    # tfaip is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by the
    # Free Software Foundation, either version 3 of the License, or (at your
    # option) any later version.
    #
    # tfaip is distributed in the hope that it will be useful, but
    # WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    # or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    # more details.
    #
    # You should have received a copy of the GNU General Public License along with
    # tfaip. If not, see http://www.gnu.org/licenses/.
    # ==============================================================================

Coding Style
------------

|tfaip| follows `Black Python style guide <https://black.readthedocs.io>`_ and `PEP 8 <https://pep8.org/>`_, in short:
* Naming of packages: lowercase, no separators, e.g., `ctctransformer`.
* Naming of files: snake case, e.g., `data_params`
* Naming of classes: camelcase, e.g., `DataParamsAtr`
* Naming of variables: snake case

To manually run ``black`` on all files call:

.. code-block:: shell

    black .

To automatically apply the black code style on all files before committing (pre-commit hook) run (once):

.. code-block:: shell

    pre-commit install

To apply all ``pre-commit`` hooks once call:

.. code-block:: shell

    pre-commit run --all-files

For integrating black in your IDE (e.g., pycharm) see `here <https://black.readthedocs.io/en/stable/integrations/editors.html>`_.

Unittests
---------

All unittests can be run by calling

.. code-block:: shell

    pytest

To run tests in parallel install ``pytest-xdist`` via ``pip install pytest-xdist`` and call

.. code-block:: shell

    pytest -n 8

to run 8 tests in parallel. Note this requires some considerable amount of RAM.

Contributing Scenarios
----------------------

We welcome users to share their scenarios that are implemented with |tfaip|.
This helps other users to become familiar with |tfaip| but also to see if there is an already existing solution for a certain problem.
Simply create an issue providing a link to your scenario and a brief explanation.
Similar to the |tfaip| examples, we expect your scenario to ship unittests and a documentation.
tfaip.trainer
==================

.. toctree::
   :maxdepth: 2
   :caption: Modules

    .callbacks.* <tfaip.trainer.callbacks>
    .optimizer.* <tfaip.trainer.optimizer>
    .scheduler.* <tfaip.trainer.scheduler>
    .warmstart.* <tfaip.trainer.warmstart>

Trainer
-------
.. automodule:: tfaip.trainer.trainer


TrainerParams
-------------
.. automodule:: tfaip.trainer.params
Parameter Classes
=================

Parameters are used to control every module in |tfaip| including hyper-paramters of training, the selection of graph, modifying the data pipeline, etc.
Each module is split into a class providing the actual implementation, and a parameter dataclass that only holds its parameters, e.g. ``Trainer`` and ``TrainerParams``.
Each ``Params``-class of `tfaip` is wrapped by ``dataclass`` and ``pai_dataclass`` which allows to easily add new parameters that support typing and will automatically be accessible from the command line.
Command line support is provided by :`paiargparse <https://github.com/Planet-AI-GmbH/paiargparse>`_ which automatically converts a dataclass hierarchy into command line arguments.

In the following, all primary parameter classes are listed

* :ref:`ScenarioBaseParams<doc.scenario:Scenario>`
  * :ref:`DataBaseParams<doc.data:data>`
  * :ref:`ModelParams<doc.model:model>`
  * :ref:`EvaluatorParams<doc.scenario:Evaluator>`
* :ref:`TrainerParams<doc.training:training>`
    * :ref:`DeviceParams<doc.device_config:Device configuration>`
    * :ref:`OptimizerParams<doc.training:optimizer>`
    * :ref:`LearningRateParams<doc.training:Learning rate>`
    * :ref:`WarmstartParams<doc.training:warm-start>`
* :ref:`LAVParams<doc.evaluation:Load and validate (LAV)>`
* :ref:`PredictorParams<doc.prediction:prediction>`


Parameters
----------

All parameters that shall be available during training or loading of a model (e.g., :ref:`resume training<doc.training:resume training>`, or :ref:`LAV<doc.evaluation:load and validate (LAV)>`) must be added to a params class within the hierarchy (see above).

Custom Parameter Classes
~~~~~~~~~~~~~~~~~~~~~~~~

You are allowed to add additional parameter classes which must be wrapped by ``@pai_dataclass`` and ``@dataclass``.
Make sure to use a descriptive name for the name and the field storing that dataclass, e.g.:

.. code-block::python

    @pai_dataclass
    @dataclass
    class BackendParams:
    # This is the additional parameter set that could define a backend (e.g. a convolutional neural net)
    conv_filters: List[int] = field(default_factory=lambda: [8, 16, 32])

    @pai_dataclass
    @dataclass
    class ModelParams(ModelBaseParams):
    # This is the default parameter set for a model
    # other params ...
    backend: BackendParams = field(default_factory=BackendParams)

Subclassing Parameter Classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Subclassing a dataclass will inherit all fields of the parent.
|tfaip| uses this to extend parameters for almost every module.
Note, that this is the reason why every field must have at least a dummy default value.
Required parameters to be set from the command line must be specified by setting the metadata: ``field(..., metadata=pai_meta(required=True))``.

Hidden Parameters
~~~~~~~~~~~~~~~~~

To hide a parameter from the command line, set its field metadata mode to ignore: ``field(..., metadata=pai_meta(mode="ignore")``.
These parameters are considered to be *state* parameters that are set or changed during training (e.g., ``TrainerParams.export_best_state``)
or are set by the scenario to exchange information (e.g., the ``DataParams`` could define ``num_classes`` which are passed to ``ModelParams.num_classes`` during the initialization).
To share parameters implement the ``__post_init__`` function of a parent dataclass if possible, else override ``ScenarioBase.create_model``, ``ScenarioBase.create_data``, etc. to specify a parameter directly before the actual class is instantiated.

Command Line
------------

The parameter hierarchy is parsed and flattened to allow to set the parameters from the command line. The following types are supported, see `here <https://github.com/Planet-AI-GmbH/paiargparse#supported-types>`_ for a full list and `examples <https://github.com/Planet-AI-GmbH/paiargparse#examples>`_:
* Primitive types: ``str``, ``int``, ``float``, ``bool``
* Enums: ``IntEnum``, ``StrEnum``
* Lists: ``List[str]``, ``List[int]``, ``List[float]``, ``List[Enum]``
* Other dataclasses defined with ``@pai_dataclass`` and ``@dataclass``, also in `Dicts` and `Lists`

Naming convention:
* Dataclasses in ``snake`` mode, the *default* (``dc_meta(mode="snake")``) are added as snake mode, e.g. ``--train.batch_size``
* Dataclasses in ``flat`` mode (``pai_meta(mode="flat")``) are added as root parameter, e.g. ``--model`` or ``--trainer``

`Meta data <https://github.com/Planet-AI-GmbH/paiargparse#meta-data>`_ can be specified to for example add a help string or change the mode.

Example:

.. code-block::python

    from dataclasses import dataclass, field
    from typing import List

    @pai_dataclass
    @dataclass
    class ExampleParams:
        example: TYPE = field(default=DEFAULT_VALUE, metadata=pai_meta(help="HELP STR"))

        # For example
        epochs: int = field(default=1000, metadata=pai_meta(help="Number of epochs to train"))
        # Or with factory
        gpus: List[int] = field(default_factory=list, metadata=pai_meta(help="GPUs to use"))


Static Parameters
-----------------

Static parameters are parameters that must be know to create a certain class of ``ModelBase``, ``DataBase``, ``GraphBase``, ``PredictorBase``, ``LAV``, ``Evaluator``, or ``RootGraph``.
Hereto add the parameter as additional argument to the respective ``__init__`` function and override the respective parameter function in ``ScenarioBase``.

See the `how to <https://github.com/Planet-AI-GmbH/tfaip/tree/master/examples/howtos/staticmodelparameters>`_ for a usage for ``ModelBase``, ``DataBase``, and ``GraphBase``.

Use this if you want to:

* pass parameters from ``DataBase`` (instantiated) to the ``ModelBase` or ``Evaluator``, e.g., the size of a loaded codec, or the tokenizer.
  See usage in the `ATR example <https://github.com/Planet-AI-GmbH/tfaip/tree/master/examples/atr/scenario.py>`_.
tfaip.model
================

ModelBase
---------
.. automodule:: tfaip.model.modelbase

ModelBaseParams
---------------
.. automodule:: tfaip.model.modelbaseparams
tfaip.data
===============

DataBase
--------
.. automodule:: tfaip.data.data

DataBaseParams
--------------
.. automodule:: tfaip.data.databaseparams
Data
====


The data class is responsible for loading (e.g., images or texts), and processing (e.g., data augmentation) of samples.
A sample is a single instance of data that is processed in the network.

A preprocessing pipeline handles the data that are fed into the :ref:`Model<doc.model:Model>` by applying a list of ``DataProcessors`` and eventually setting up a ``tf.data.dataset``.
Furthermore, a postprocessing pipeline can be added to process the data (e.g. pruning) that is produced by the network.

Figure 1 shows an overview of the pre-processing pipeline provided by |tfaip|.
The actual pipeline (lower half, purple) is created by a list of ``DataProcessorParams`` (upper half, orange).
The input of the pipeline  (``DataGeneratorParams``, golden) is can depend on the steps training, LAV, or prediction.
Usually, this is ``ListsFileDataGeneratorParams`` which will simply generate filenames parsed from text file (see `Scenario<doc.scenario:ListFileScenario>` for setting up a ``ListFileScenario``).

.. figure:: resources/datapipeline_workflow.png
    :alt: Data-Pipeline Workflow

    Figure 1: The data processing pipeline


Implementation
--------------

A data class has to extend ``from tfaip.data.data import DataBase`` and implement the following

.. code-block:: python

    @pai_dataclass
    @dataclass
    class DataParams(DataBaseParams):
        @staticmethod
        def cls():
            return Data

    class Data(DataBase[DataParams]):
        @abstractmethod
        def _input_layer_specs(self) -> Dict[str, tf.TensorSpec]:
            raise NotImplementedError

        @abstractmethod
        def _target_layer_specs(self) -> Dict[str, tf.TensorSpec]:
            raise NotImplementedError

        @classmethod
        def default_params(cls):
            # can optionally be overwritten
            p = super().default_params()
            p.pre_proc.processors = [
                # Optionally pass a list of DataProcessorParams to define the input pipeline (see orange boxes in Figure 1)
            ]
            return p

Extend DataBaseParams to return the actual implementation in the ``cls`` method. It is recommended to do this in a separate python file (e.g. called ``params.py``).

Required Overrides
~~~~~~~~~~~~~~~~~~

The only two required methods to overwrite return dicts mapping from node names to ``tf.TensorSpec``, e.g.,

.. code-block:: python

    def _input_layer_specs(self) -> Dict[str, tf.TensorSpec]:
        # Shape and type of the input data for the graph
        return {'img': tf.TensorSpec(shape=[28, 28], dtype='float32')}

    def _target_layer_specs(self) -> Dict[str, tf.TensorSpec]:
        # Shape and type of the input data for the graph
        return {'tgt': tf.TensorSpec(shape=[1], dtype='int32')}

for a fixed 28x28 input image, and a class label as target.
The batch dimension must be omitted as batching is performed automatically.

**Important**: No scalars are allowed. Allowed shapes ``[1]`` for a scalar, ``[N]`` for vector, ``[H,W,C]`` for a tensor, ...

Optional Overrides
~~~~~~~~~~~~~~~~~~

By default, all arrays are automatically padded with zeros, overwrite ``_padding_values`` to provide different values.

.. code-block:: python

    def _padding_values(self) -> Dict[str, AnyNumpy]:
       return {'img': 255}

Pre/Post-Processing Pipelines
-----------------------------

``DataProcessors`` enable to modify data by consuming and producing a list of ``Sample``.
There are two types of data processors:

* A ``MappingDataProcessor`` is similar to the default mapping function of Python and produces one output sample for each input sample (a ``1`` to ``1`` mapping).
* A ``GeneratingDataProcessor`` is similar to a generator in Python: ``N`` samples are consumed but ``M`` samples are produces, whereby ``M<N`` and ``M>N`` are equally legitimated.

``DataProcessors`` are defined and constructed by their corresponding ``DataProcessorParams``.

DataProcessor Definition
~~~~~~~~~~~~~~~~~~~~~~~~

The following shows a simple ``DataProcessor`` that will load an image and its GT from the drive based on filenames.
Optionally, add parameters (e.g., the color mode in this case) to the processor params.

.. code-block:: python

    @pai_dataclass
    @dataclass
    class LoadSampleProcessorParams(DataProcessorParams):
        def cls(self):
            return LoadSampleProcessor

    class LoadSampleProcessor(MappingDataProcessor):
        def apply(self, sample):
            return (
                sample.
                    new_inputs(load_image(sample.inputs)).   # Change the inputs
                    new_targets(load_gt(sample.targets))     # Change the targets
            )

Set Up
~~~~~~

Override the ``default_params`` of the ``Data``-class to set up the default preprocessing pipeline:

.. code-block:: python

    @classmethod
    def default_params(cls) -> ListFileDataParams:
        # setup default processing pipeline
        params = super().default_params()
        params.pre_proc = SequentialDataProcessorPipeline(
            run_parallel=True,  # Set to True if this pipeline shall run in parallel
            processors=[
                LoadSampleProcessorParams(),  # Load the sample and its GT
                AugmentProcessorParams(modes={PipelineMode.Training}),  # Apply data augmentation, but only during training
                FinalizeProcessorParams(),  # Finalize the sample, i.e., bring it in the correct form matchin the input and target layer specs
            ],
        )
        return params

Modes of DataProcessor
""""""""""""""""""""""

As shown in the previous code, the ``DataProcessorParams`` provide a field when to apply this processor.
Here, the ``AugmentProcessorParams`` (i.e., data augmentation) shall only be applied on the training pipeline.

Resources
---------

Quite often your model requires resources for training but also for the later application.
A typical resource is a file that needs to be served when exporting the model.
The dataclass automatically handles the export of ``Resources`` via its ``ResourceManager`` by copying the resources to the export dirs and automatically adapting the search path of the resource.
Define a ``Resource`` in your parameters:

.. code-block:: python

    @pai_dataclass
    @dataclass
    class DataParams(DataBaseParams):
        charmap: Resource = field(default=None,
                                  metadata={**pai_meta(help="File specifying the character map used", required=True),
                                            **config(encoder=Resource.encode, decoder=Resource.decode)}
                                  )

In this example, the character map will automatically be copied to the ``resources`` dir in the exported model.

Development
-----------

This section provides additional information about the actual implementation of the data pipeline in |tfaip|.
Read this if you are interested in understanding or extending |tfaip|.

The following image provides an overview of all relevant classes

.. figure:: resources/datapipeline_overview.png
    :alt: Data-Pipeline Overview

    Figure 2: Overview of all classes withing the data module

Have a look at the code documentation for a description of the individual classes, in the following is only a small overview:

* Red: ``DataBaseParams`` and ``DataBase`` define the overall structure of a data pipeline by connecting ``DataProcessorParams``.
  By calling ``get_or_create_pipeline``, or ``create_pipeline`` a ``DataPipeline`` will be prepared.
* Yellow: The ``DataGenerator`` provide the data for the preprocessing pipeline. The actual implementation depends on the scenario and mode (e.g. different data sources for training, lav and prediction)
* Orange: The creation and definition of the different ``DataProcessor`` types. A user should only override ``MappingDataProcessor`` and ``GeneratingDataProcessor``. Multiple ``MappingDataProcessors`` will be joined to a ``SequenceProcessor`` for faster execution (send as a complete block to a spawned process upon parallelization).
* Green: Setting up of the actual Pipeline using parameters. By default, a ``SequenceProcessorPipelineParams`` should suffice, if however `GeneratingDataProcessors` play a role it might be sensible to provide a custom grouping which can be done via ``CompoundProcessorPipelineParams``
  Note, the ``SequenceProcessorPipelineParams`` will be converted automatically to a ``CompoundProcessorPipelineParams`` (see Figure 1).
* Blue: These classes are use for the actual data processing of a set of ``MappingDataProcessors`` in a ``MappingSampleProcessorPipeline`` or one ``GeneratingDataProcessor`` in a ``GeneratingSampleProcessorPipeline``.
  Each class has a corresponding parallel version (see lower half of Figure 1).
  Construction of the actual processors (calling ``DataProcessorParams.create``) is performed within these classes to ensure that only the parameters are passed to a spawned process not the actual class (which might not be serializable via pickle).
* Yellow (lower right): The ``TFDatasetGenerator`` can optionally be overwritten in a DataPipeline to change the creation of the ``tf.data.Dataset`` or to inject additional data mappings performed in Tensorflow.
tfaip.predict
==================

PredictorBase
-------------
.. automodule:: tfaip.predict.predictorbase

PredictorParams
---------------
.. automodule:: tfaip.predict.params

Predictor
---------
.. automodule:: tfaip.predict.predictor

MultiModelPredictor
-------------------
.. automodule:: tfaip.predict.multimodelpredictor
Training
========

This section deals about how to :ref:`start <doc.training:launch training>` or :ref:`resume<doc.training:resume training>` via the command line and how to modify the :ref:`training hyper-parameters<doc.training:parameters>`.

Command-Line
------------

Use the command line to launch or resume a training.

Launch Training
~~~~~~~~~~~~~~~

Training of a scenario is performed by the ``train.py`` script which is launched by ``tfaip-train``:

.. code-block:: shell

    tfaip-train {SCENARIO_MODULE} PARAMS

The ``SCENARIO_MODULE`` must be located in the ``PYTHONPATH``.
Either give the full name to scenario class, the scenario file, or the parent module:

.. code-block:: shell

    tfaip-train tfaip.scenario.tutorial.full.scenario:TutorialScenario ...
    tfaip-train tfaip.scenario.tutorial.full.scenario ...
    tfaip-train tfaip.scenario.tutorial.full ...

The first option is required if there are multiple Scenarios located in a single python file.

Resume Training
~~~~~~~~~~~~~~~

Training can be resumed using ``python resume_training.py`` or ``tfaip-resume-training``.
Specify a checkpoint and that's it, e.g.:

.. code-block:: shell

    tfaip-train tfaip.scenario.tutorial.full --trainer.output_dir model_output
    # CANCELLED
    tfaip-resume-training model_output

Resume training and adapt parameters (e.g., extend the ``epochs``): This is not supported directly, yet, however you can manipulate the ``trainer_params`` in the ``output_dir`` and manually adapt the settings.

Debugging Training
~~~~~~~~~~~~~~~~~~

To use the debugger call ``tfaip/scripts/train.py`` instead of ``tfaip-train`` in a run-configuration (e.g., in PyCharm).
See :ref:`here<doc.debugging:Debugging>` for further information.

Validation During Training
~~~~~~~~~~~~~~~~~~~~~~~~~~

Validation during training is used to compute and print the validation performance and to compute the best performing model (see :ref:`early-stopping<doc.training:early-stopping>`).
Validation is optional and will be performed on the provided validation list, e.g. ``--val.lists`` for a :ref:`ListFileScenario<doc.scenario:ListFileScenario>`.

Default
"""""""

By default, validation will be performed by keras and thereto computes the provided loss and metrics of the respective :ref:`Model<doc.model:Model>`.
The default validation is run every epoch which is defined by the ``TrainerParams.val_every_n`` parameter which defaults to 1.

Note that the default validation will be performed on the training graph.
Therefore, if the scenario requires a different prediction graph (e.g. for a Encoder-Decoder-Model) use the extended validation.

Extended
""""""""

It is possible to run a clean validation via [LAV](10_evaluation.md) during training, i.e., loading the current prediction graph and applying it on the validation list.
To enable set the ``TrainerParams.lav_every_n`` parameter to specify on which epochs to run (e.g., ``--trainer.lav_every_n=1``. First and last epochs are always validated if LAV is enabled).
By default, LAV will use the validation generator, but this can be overwritten in the respective ``TrainingPipelineGeneratorParams``.
Note, a :ref:`ListFileScenario<doc.scenario:listfilescenario>` already provides an additional parameter ``lav.lists`` which defaults to ``val.lists``.

LAV will then evaluate all metrics (including the ones of the :ref:`Evaluator<doc.scenario:evaluator>`) and print them (also to the :ref:`Tensorboard<doc.training:tensorboard>`).

Parameters
----------

The most important parameter during training is the ``output_dir`` which defines where to store the log, the exported models, and checkpoints.
Set with ``--trainer.output_dir MY_OUTPUT_DIR``.
Other parameters are introduced in the following.

Logging
~~~~~~~

``tfaip`` uses the ``logging`` module of python.

Set up the log level using the ``TFAIP_LOG_LEVEL`` environment variable, e.g. ``TFAIP_LOG_LEVEL=debug``.
By default, the log level is set to ``info``.

Logging events written to ``logging`` are printed and also written to the ``train.log``.

Learning rate
~~~~~~~~~~~~~

The learning rate can be adapted using the ``trainer.learning_rate`` field which must be set to a ``LearningRateParams`` structure.
The `LearningRateParams` always provide a ``lr``-field to modify the overall learning rate which defaults to ``0.001``.

Example: ``--learning_rate.lr 0.001``.

To change the schedule, set the learning rate field directly: ``--learning_rate ExponentialDecay``.
The parameters of the schedule can be set similarly to above, e.g., ``--learning_rate.decay 0.95``.

Optimizer
~~~~~~~~~

The optimizer of the trainer can be changed and adapted via the ``trainer.optimizer`` field which is a ``OptimizerParams`` structure.
|tfaip| supports different optimizers by default: ``Adam``, ``Adamax``, ``AdaBelief``, ``RMSprop``, ``SGD``. Each one comes with is custom parameters.

Example: ``--optimizer Adamax``.

To adap the parameters of the optimizer call, e.g., ``--optimizer.epsilon 1e-7``

Gradient-Clipping
"""""""""""""""""

Each optimizer supports gradient-clipping based on the `Tensorflow-Optimizer <https://www.tensorflow.org/api_docs/python/tf/keras/optimizers/Optimizer>`_: ``clip_value``, ``clip_norm``, ``clip_global_norm``.

Example: ``--optimizer.clip_global_norm 5``.

Weight-Decay
""""""""""""

A global weight decay (applied to all weights) is provided by the ``Adam`` and ``SGD`` optimizer with their additional ``weight_decay`` field which defaults to ``0.0``.
Alternatively, you can implement weight-decay directly when define layers, as recommended by `Tensorflow <https://www.tensorflow.org/tutorials/keras/overfit_and_underfit#combined_l2_dropout>`_.

Example: ``--optimizer.weight_decay 0.00002``.

EMA-Weights
~~~~~~~~~~~

|tfaip| supports to compute an exponential moving average (EMA) on the weights which is enabled by ``trainer.calc_ema``.
In many cases this leads to improved results on the drawback that more GPU-memory is required during training.
Note that the exported models are always saved along with the EMAs, while checkpoints comprise both the EMAs and the last weights, i.e. the current state of the trainer.
Furthermore, EMA weights are used for validation.

The parameter expects the ema rate, if ``-1``, the rate is computed automatically.

Example: ``--trainer.calc_ema 0.99``

No train scope
~~~~~~~~~~~~~~

Use ``trainer.no_train_scope`` to pass a regex which defines which layers to exclude from training.
Note, if a parent layer is matched, all children will also be not trained.

Example: ``--trainer.no_train_scope '.*conv.*'``.

Warm-Start
~~~~~~~~~~

Warm-starting a model before training with predefined weights is supported.
See ``WarmstartParams`` for all options.

Example: ``--warmstart.model PATH_TO_WARMSTART_MODEL``

To customise warm starting, pass custom parameters ``--warmstart MyWarmstarter`` that need to inherit ``WarmstartParams``.


Devices
~~~~~~~

See :ref:`DeviceConfig <doc.device_config:Device Configuration>` which is set at ``trainer.device``

Example: ``--device.gpus 0 1``.

Early-Stopping
~~~~~~~~~~~~~~

Setting up early stopping via the ``EarlyStoppingParams`` in ``trainer.early_stopping`` allows the trainer to automatically determine on different constraints when to stop.
Monitoring is based on the best model determined by the :ref:`settings of the model<doc.model:logging the best model>`.

To enable early stopping set ``n_to_go`` to the number of epochs after which to stop training if no improvements was achieved.
To modify not to test every epoch, increase the ``frequency`` (defaults to ``1``, i.e., test every epoch).

Furthermore, you can specify a ``lower_threshold`` and an ``upper_threshold``.
Depending on the ``mode`` during setup, the monitored value must be at least one threshold for early stopping to start, or if the other one is reached, training is stopped immediately.
For example if monitoring the ``max`` ``accurary`` and ``upper_threshold=0.99`` and ``lower_threshold=0.5``, training will not stop until an ``accuracy`` of at least ``0.5`` is reached, then early stopping could kick in.
However, if ``accuracy`` exceeds or equals to ``0.99`` training is stopped immediately. Therefore, setting ``upper_threshold=1`` is sensible if monitoring accuracies.


Export of Checkpoints and Saved Models
--------------------------------------

During training several models/weights are logged:

* ``checkpoint``: stores the complete state (weights and optimizer) used to resume the training.
* ``best``: saved model, that stores the best model on the validation set
* ``export``: saved model, the final state of the model (last model)

There are two formats:

* checkpoint: only the weights are stored. Can only be used to continue the training or for :ref:`warm-start<doc.training:warm-start>`
* saved model (serving): the model and weights are stored. Can be used in LAV and from other apis (e.g., java). Can also be used for :ref:`warm-start<doc.training:warm-start>`

Customizing the Exported Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the exported graph uses all available inputs and outputs tensors, i.e. the inputs and outputs of the :ref:`created graph<doc.graph:graph>`.
To modify this, or to add additional graphs for export, override the `_export_graphs<doc.model:exporting additional graphs>` function of the `ModelBase`.

Tensorboard
-----------

Training can be monitored by the `TensorBoard <https://www.tensorflow.org/tensorboard/>`_.
Hereto, |tfaip| automatically stores the metrics and losses on the train, validation, and lav datasets to the ``output_dir`` which can be displayed by the TensorBoard.
Launch the tensorboard via ``tensorboard --log_dir PATH_TO_THE_MODELS --bind_all``.

Additional output such as images or PR-curves can be setup in the :ref:`model <doc.model:tensorboard>`.

Benchmarking
------------


Profiling Training
~~~~~~~~~~~~~~~~~~

The full training can be profiled using the `Tensorboard`:

* Install requirement ``tensorboard_plugin_profile``
* Set ``--trainer.profile True``


Benchmarking the Input Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Quite often, the bottleneck is not the model but the input data pipeline that is not able to produce enough samples per second.

Replace ``tfaip-train`` with ``tfaip-benchmark-input-pipeline`` and run to profile the number of samples per second.
By default the benchmark will run infinitely, so terminate to process to stop the benchmark.
tfaip.trainer.optimizer
============================


Optimizers
----------
.. automodule:: tfaip.trainer.optimizer.optimizers

Gradient Accumulation Optimizer
-------------------------------
.. automodule:: tfaip.trainer.optimizer.gradient_accumulation_optimizer

Weights Moving Average
----------------------
.. automodule:: tfaip.trainer.optimizer.weights_moving_average
Data Pipeline
=============

The data pipeline can be used as stand-alone without the requirement to setup a full ``Data`` class.
This can be useful to simply apply a data pipeline comprising a list of processors.

.. code-block:: python

    from tfaip.imports import DataBaseParams

    # Obtain or create data params from a scenario
    # these params also include some
    data_params = scenario.data_params

    # Here run the pre-proc pipeline, but this can be an arbitrary pipeline
    pipeline_params = data_params.pre_proc

    # Create some samples (empty here), can also be a DataGeneratorParams instance
    samples = [Sample()]

    # create the pipeline
    pipeline = pipeline_params.create(DataPipelineParams(num_processes=8, limit=100), data_params)

    # apply the pipeline on the samples and retrieve the output
    for i, d in enumerate(pipeline.apply(samples)):
        print(i, d)
Load and Validate (LAV)
=======================

This file describes how to evaluate a trained model performed by Load-And-Validate (LAV).
LAV loads an :ref:`exported model<doc.training:Export of Checkpoints and Saved Models>`, applies it on provided data, and finally compute the :ref:`metrics<doc.model:metric>` and :ref:`losses<doc.model:loss>` defined in the :ref:`model<doc.model:model>` but also the output of an optionally provided :ref:`Evaluator<doc.scenario:evaluator>`.

Calling LAV via the command line
--------------------------------

To load and validate a model use a [saved model](http://gitea.planet-ai.de/pai/tf2_aip/wiki/Model-exporting) (e.g. export or best) and the `tfaip-lav` script, usage:

.. code-block::shell

    tfaip-lav --export_dir PATH_TO_SAVED_MODEL

Parameters
~~~~~~~~~~

LAVParams
"""""""""

The ``PipelineParams`` of the ``LAVParams`` are accessed directly via `--lav`, e.g.

.. code-block::shell

    --lav.batch_size 5
    --lav.num_processes 32

The ``DeviceParams`` can be set by, e.g.:

.. code-block::shell

    --lav.device.gpus 0

DataGeneratorParams
"""""""""""""""""""

``tfaip-lav`` allows to adapt the data generator parameters which is useful to change the evaluation data (by default the validation generator when training is used).
For a :ref:`ListFileScenario<doc.scenario:listfilescenario>`, the evaluation list can be changed by:

.. code-block::shell

    --data.lists OTHER_VAL_LIST


Other Parameters
""""""""""""""""

Specify

* ``--run-eagerly`` to run lav in eager-mode (useful for debugging)
* ``--dump`` to dump the targets and predictions to a pickle file (implement a custom `LAV` and `LAV.extract_dump_data` to modify the dump)


Implement a Custom LAV
----------------------

The default ``LAV`` module can be extended to change the default behaviour.
Do not forget the register it in the :ref:`Scenario<doc.scenario:Scenario>` at ``ScenarioBase.lav_cls()`` or ``ScenarioBase.create_lav()``.
tfaip.evaluator
====================

EvaluatorBase
-------------
.. automodule:: tfaip.evaluator.evaluator

EvaluatorParams
---------------
.. automodule:: tfaip.evaluator.params
tfaip
=====

This documentation contains all necessary information to effectively work with |tfaip|.
It depicts the structural concepts of |tfaip| and provides a step by step tour throughout all (necessary) components.

|tfaip| offers two **tutorials**:
If being new to |tfaip|, have a look at the `minimal tutorial <https://github.com/Planet-AI-GmbH/tfaip/tree/master/examples/tutorial/min>`_ to get an impression how |tfaip| is structured.
To quickly set up a new scenario copy either a `tutorial <https://github.com/Planet-AI-GmbH/tfaip/tree/master/examples/tutorial>`_ or a `Template <https://github.com/Planet-AI-GmbH/tfaip/tree/master/examples/template>`_.
Templates show the most common functions that must be implemented (search for ``raise NotImplemented`` and comments written in squared brackets) to set up a custom scenario.

.. toctree::
   :maxdepth: 1
   :caption: Documentation

    Installation <doc.installation>
    Parameters <doc.parameters>
    Scenario <doc.scenario>
    Data <doc.data>
    Data Pipeline <doc.data.pipeline>
    Model <doc.model>
    Graph <doc.graph>
    Training <doc.training>
    Evaluation <doc.evaluation>
    Prediction <doc.prediction>
    Device Configuration <doc.device_config>
    Debugging <doc.debugging>
    Development <doc.development>

.. toctree::
   :maxdepth: 2
   :caption: Tutorials, Template, Examples

    Minimal Tutorial <https://github.com/Planet-AI-GmbH/tfaip/tree/master/examples/tutorial/min>
    Full Tutorial <https://github.com/Planet-AI-GmbH/tfaip/tree/master/examples/tutorial/full>
    General Template <https://github.com/Planet-AI-GmbH/tfaip/tree/master/examples/template/general>
    ListFile Template <https://github.com/Planet-AI-GmbH/tfaip/tree/master/examples/template/listfile>
    Further Examples <https://github.com/Planet-AI-GmbH/tfaip_example_scenarios>

.. toctree::
   :maxdepth: 2
   :caption: Modules

    tfaip.data.* <tfaip.data>
    tfaip.evaluator.* <tfaip.evaluator>
    tfaip.lav.* <tfaip.lav>
    tfaip.model.* <tfaip.model>
    tfaip.scenario.* <tfaip.scenario>
    tfaip.predict.* <tfaip.predict>
    tfaip.trainer.* <tfaip.trainer>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Graph
=====

Each |tfaip| graph must inherit ``GenericGraphBase``, or its subclass ``GraphBase`` which already encapsulates some optional methods for the graph construction.

GraphBase
---------

An example Multi-Layer-Perceptron (MLP) graph for MNIST is shown in the following code example (excerpt of a Tutorial):

.. code-block:: python

    class TutorialGraph(GraphBase[TutorialModelParams]):
        def __init__(self, params: TutorialModelParams, name='conv', **kwargs):
            super(TutorialGraph, self).__init__(params, name=name, **kwargs)
            # Create all layers
            self.flatten = Flatten()
            self.ff = Dense(128, name='f_ff', activation='relu')
            self.logits = Dense(self._params.n_classes, activation=None, name='classify')

        def build_graph(self, inputs, training=None):
            # Connect all layers and return a dict of the outputs
            rescaled_img = K.cast(inputs['img'], dtype='float32') / 255
            logits = self.logits(self.ff(self.flatten(rescaled_img)))
            pred = K.softmax(logits, axis=-1)
            cls = K.argmax(pred, axis=-1)
            out = {'pred': pred, 'logits': logits, 'class': cls}
            return out


Layers are instantiated in the ``__init__`` function and applied in the ``build_graph`` function which is the sole abstract method of ``GraphBase``


GenericGraphBase
----------------

The ``GenericGraphBase`` class provides more flexibility when creating a graph, including different graphs for training (``build_train_graph``) and prediction (``build_prediction_graph``).
This is for example required to implement a sequence-to-sequence model with a decoder that depends on the model (teacher-forcing during training, decoding during prediction).
Furthermode, the method ``pre_proc_targets`` can be overwritten to apply some preprocessing on the dataset targets that are then fed into the metrics.
tfaip.scenario
===================

ScenarioBase
------------
.. automodule:: tfaip.scenario.scenariobase

ScenarioBaseParams
------------------
.. automodule:: tfaip.scenario.scenariobaseparams

.listfile
---------

ListFileScenario
~~~~~~~~~~~~~~~~
.. automodule:: tfaip.scenario.listfile.listfilescenario

Parameters
~~~~~~~~~~
.. automodule:: tfaip.scenario.listfile.params

Data Generator
~~~~~~~~~~~~~~
.. automodule:: tfaip.scenario.listfile.datagenerator
tfaip.lav
====================

LAV
---
.. automodule:: tfaip.lav.lav

MultiLAV
--------
.. automodule:: tfaip.lav.multilav

LAVParams
---------
.. automodule:: tfaip.lav.params

Callbacks
---------
.. automodule:: tfaip.lav.callbacks.lav_callback
.. automodule:: tfaip.lav.callbacks.dump_results
Device Configuration
====================

The definition which devices are visible and used by |tfaip| and thus Tensorflow is defined via the ``DeviceConfigParams``.
They are used during :ref:`training<doc.training:devices>`, :ref:`validation<doc.evaluation:LAVParams>`, and prediction.
The structure also allows to modify the ``DistributionStrategy`` when training on several GPUs.

Selecting GPUs
--------------

By default, |tfaip| will not use any GPU.
Modify the used GPUs by setting the ``DeviceConfigParams.gpus`` flag which expects a list of GPU indices, e.g.

.. code-block::

    --trainer.device.gpus 0 2

to use the GPUs with index 0 and 2.

Alternatively, if the environment variable ``CUDA_VISIBLE_DEVICES`` is set, |tfaip| will use these GPUs for training.

Multi-GPU setup
---------------

Set the ``DeviceConfigParams.dist_strategy`` flag to ``mirror`` or ``central_storage`` to define how to use multiple devices.
Scenario
========

A scenario is the container for implementing a custom task that shall be solved by |tfaip|.
The scenario glues together all components that must and can be implemented, see Figure 1.

.. figure:: resources/scenario.png
    :alt: Scenario Overview

    Figure 1: Scenario Overview

The ``ScenarioBase`` (green) creates the :ref:`Model<doc.model:model>` (blue) and :ref:`Data<doc.data:data>` (purple).
:ref:`Training<doc.training:Training>` (red), :ref:`Load And Validate<doc.evaluation:Load and Validate (LAV)>` (LAV, orange), and :ref:`Prediction<doc.prediction:Prediction>` (yellow) access the scenario but are also instantiated by the ``ScenarioBase`` to allow to implement custom overrides.


Scenario Directory
------------------

Each scenario should be set up in a directory that comprises the following files:
* ``scenario.py`` which defines the basic scenario information
* ``model.py`` which defines the :ref:`model<doc.model:model>` (neural net)
* ``data.py`` which defines the input :ref:`data<doc.data:data>` pipeline of the model
* ``graphs.py`` which defines the :ref:`graph(s)<doc.graph:graph>` of the model
* ``params.py`` (optional but recommended) which stores the parameters of data, model, and the scenario. (Note, the data params, if implemented should always be implemented in a different file than data itself else there will occurr warnings in the data pipeline).

If ``tfaip`` was set up by cloning the repository, the scenario class should be located at ``tfaip/scenario/XXX/scenario.py`` where ``XXX`` is the scenario name, e.g. ``atr``.
If :ref:`installed<doc.installation:installation>` via ``pypi`` (i.e., ``pip install tfaip``), an arbitrary directory can be used for the location of a scenario.

Implementing a Scenario
-----------------------

To implement a ``ScenarioBase``, first :ref:`setup a directory and the files<doc.scenario:scenario directory>`. In ``scenario.py`` (or ``params.py``) implement ``MyScenarioParams``:

.. code-block:: python

    @pai_dataclass
    @dataclass
    class MyScenarioParams(ScenarioBaseParams[MyDataParams, MyModelParams]):
        pass

The ``MyDataParams`` and ``MyModelParams`` implement ``DataBaseParams`` and ``ModelBaseParams`` to define the parameters for the data and the model.

Next, implement the actual scenario:

.. code-block:: python

    class MyScenario(ScenarioBase[MyScenarioParams, MyTrainerPipelineParams]):
        pass

The ``MyTrainerPipelineParams`` define how the input data source for training and extend either ``TrainerPipelineParamsBase`` or ``TrainerPipelineParams``.
The derived ``ListFileScenario`` already implements the ``TrainerPipelineParams`` by assuming a list file as input (see :ref:`here<doc.scenario:listfilescenario>`).

Development
-----------

The ``Scenario`` defines several ``Generics`` that are used for instantiation of the actual classes of ``TModel``, ``TData``, ``TScenarioParams``, and the ``TTrainerPipelineParams``.
The ``ListFileScenario`` replaces ``TTrainerPipelineParams`` by ``ListFileTrainerPipelineParams``.


Additional Modules
------------------

In the following, additional methods/functionality of a scenario that can optionally be implemented is listed.

Evaluator
~~~~~~~~~

Quite often, defining metrics is difficult in pure Tensorflow-Operations while it is trivial using python and numpy.
Furthermore, some metrics should also first be computed after [post-processing](04_data.md)
For this purpose, |tfaip| offers the ``Evaluator`` which is similar to a ``keras.Metric`` however with the advantage that anything can be computed with most flexibility.
An ``Evaluator`` can optionally be parametrized by ``EvaluatorParams``.

Similar to a ``keras.Metric`` the ``Evaluator`` requires to overwrite two functions, namely ``update_state`` and ``result``.
``update_state`` receives a post-processed (un-batched) Sample and should update an internal state.
Finally, ``result`` shall yield a dictionary of the metrics.

The ``Evaluator`` follows the ``context``-design of Python: A metric is ``__enter__``-ed before the validation, and ``__exit__``-ed after receiving the result.
Use this mechanism to clear the internal state.

The ``Evaluator`` is attached to a ``Scenario`` using the ``evaluator_cls``-method.


Example
"""""""
The full tutorial provides an example:

.. code-block:: python

    class MNISTEvaluator(Evaluator):
        def __init__(self, params):
            super(MNISTEvaluator, self).__init__(params)
            self.true_count = 0
            self.total_count = 0

        def __enter__(self):
            self.true_count = 0
            self.total_count = 0

        def update_state(self, sample: Sample):
            self.total_count += 1
            self.true_count += np.sum(sample.targets['gt'] == sample.outputs['class'])

        def result(self) -> Dict[str, AnyNumpy]:
            return {'eval_acc': self.true_count / self.total_count}

Add this in the ``Scenario``:

.. code-block:: python

    @classmethod
    def evaluator_cls(cls):
        return MNISTEvaluator

Tensorboard
"""""""""""

During training, the computed metrics by ``result`` will be written to the Tensorboard.
This also allows computing custom data (e.g., images or PR-curves) within the ``Evaluator``.
The :ref:`model<doc.model:tensorboard>` defines how to write arbitrary data to the Tensorboard.

ListFileScenario
----------------
The ``ListFileScenario`` is an abstract ``ScenarioBase`` that already provides some additional functionality if using list files as the input source.
A list file is a simple text file where each line is the path to a sample, e.g. an image:

.. code-block::

    path/to/image_001.png
    path/to/image_002.png
    path/to/image_003.png
    ...

The following shows how to extend a ``ListFileScenario``:
Assume the new scenario has the model ``Model`` and corresponding params ``ModelParams``, ``Data`` and corresponding ``DataParams``, and works with list files.
The new scenario, here called ``Scenario`` requires to set up its params the corresponding implementation.
Note, that both classes are empty since in most cases no extra functionality is required.


.. code-block:: python

    @pai_dataclass
    @dataclass
    class ScenarioParams(ScenarioBaseParams[DataParams, ModelParams]):
        pass

    class Scenario(ListFileScenario(ScenarioParams)):
        pass
tfaip.trainer.warmstart
============================

WarmStartParams
---------------
.. automodule:: tfaip.trainer.warmstart.warmstart_params

WarmStarter
-----------
.. automodule:: tfaip.trainer.warmstart.warmstarter
