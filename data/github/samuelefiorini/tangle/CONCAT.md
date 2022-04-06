# TANGLE
**TANGLE** is a timespan-guided neural attention mechanism that can be used to integrate sequences of Medicare services and the time span between them to provide interpretable patients representations, which can be further used for some target prediction task.

Predictive capabilities of **TANGLE** are demonstrated here with an application on the 10\% publicly available sample of deidentified, individual level, linked Medicare Benefits Schedule (MBS) and Pharmaceutical Benefits Scheme (PBS) electronic databases of Australia. Starting from sequences of MBS-items and timespans, **TANGLE** can predict which diabetic patient, currently on metformin only, is likely to be prescribed with a different type of diabetes controlling drug in the near future.

## Installation
**TANGLE** is developed using [keras](https://keras.io/) and currently requires Python 2.7.

To use **TANGLE** for your sequence classification task, you simply need to install the core `tangle` package and its dependencies.

If you have cloned the repository, run the following command from the root of the repository:

`$ cd tangle`

`$ python setup.py install`

If you do not wish to clone the repository, you can install using:

`$ pip install git+https://github.com/samuelefiorini/tangle.git`

## More info
- **Note**: the MBS-PBS 10% is no longer publicly available (see [here](http://www.pbs.gov.au/info/news/2016/08/public-release-of-linkable-10-percent-mbs-and-pbs-data)) and it is not part of this repository.
- For more info see: [pbs.gov.au](http://www.pbs.gov.au/info/news/2016/08/public-release-of-linkable-10-percent-mbs-and-pbs-data).
# Scripts

The scripts in this folder are meant to extract relevant information from the raw MBS-PBS 10% dataset and to fit/cross-validate linear and recurrent models.
To run these scripts you are supposed to have to have organized the MBS-PBS 10% dataset in a folder (*e.g.:* `../../../data`).

Python 2.7 is required.

## Usage and details
### - `get_population_of_interest.py`

Usage example:

`$ python get_population_of_interest.py --root <ROOT> --output <PATH-TO-OUTPUT>`

or, equivalently:

`$ python get_population_of_interest.py -r <ROOT> -o <PATH-TO-OUTPUT>`

This script finds continuously and consistently concessional subjects in the
PBS files (2008 - 2014) that were prescribed to diabetes control drugs.
In order to be considered continuously and consistently concessional, a subject
must:

    1) continuously use the concessional cards, i.e.: they use it for at least
       50% of the observation years,
    2) consistently satisfy the condition at point 1), i.e.: for at least 50%
       of the PBS benefit items each year.

This script will then produce a list of continuously and consistently
concessional subjects that use diabetes control drugs.

### - `assign_labels.py`

Usage example:

`$ python assign_labels.py --root <ROOT> --skip_input_check --source <PATH-TO-SOURCE> --output <PATH-TO-OUTPUT>`

or, equivalently:

`$ python labels_assignment.py -r <ROOT> -sic -s <PATH-TO-SOURCE> -o <PATH-TO-OUTPUT>`

This script aims at finding the patient identifiers of the positive and negative
classes, where:
+ positive class (`y = 1`): subjects that continuously and consistently use
  their concessional card in the observation years to buy diabetes drugs,
- negative class (`y = 0`): subjects that continuously and consistently use
  their concessional card but were never prescribed to diabetes control drugs.

 This script also extracts two other labels:

 `MET_ONLY` - *i.e.*: patients that are using metformin ONLY

 `MET_AFTER` - *i.e.*: patients that after a first metformin prescription started to use another diabetes controlling drug.

### - `extract_sequences.py`

Usage example:

`$ python extract_sequences.py -root <ROOT> --skip_input_check --exclude_pregnancy --source <PATH-TO-SOURCE> --output <PATH-TO-OUTPUT>`

or, equivalently:

`$ python extract_sequences.py -r <ROOT> -sic -ep -s <PATH-TO-SOURCE> --output <PATH-TO-OUTPUT>`

This script extracts the raw sequences from the MBS files. An example of
sequence is `1256 0 56489 12 ...` where odd entries are MBS items
and even entries are days between each visit.

### - `matching_step1.py`

Usage example:

`$ python scripts/matching_step1.py --source <PATH-TO-SOURCE> --output <PATH-TO-OUTPUT>`

or, equivalently:

`$ python scripts/matching_step1.py -s <PATH-TO-SOURCE> -o <PATH-TO-OUTPUT>`

This script prepares the data for the actual matching procedure (see `matching_step2.R`).

### - `matching_step2.R`

Usage example:

`$ Rscript scripts/matching_step2.R`

Run the matching algorithm by CEM package and generate a `matched_CEM_table.csv` file.

### - `single_train.py`

Usage example:

`$ python scripts/single_train.py --labels <PATH-TO matched_CEM_table.csv> --data <PATH-TO raw_data_.pkl> --embedding <PATH-TO embedding.100d.csv> --model tangle --output <PATH-TO-OUTPUT>`

or, equivalently:

`$ python scripts/single_train.py -l <PATH-TO matched_CEM_table.csv> -d <PATH-TO raw_data_.pkl> -e <PATH-TO embedding.100d.csv> -m tangle -o <PATH-TO-OUTPUT>`

Fit **TANGLE** on a random training/validation/test split of the matched dataset.

### - `cross_validate.py`

Usage example:

`$ python scripts/cross_validate.py --n_splits N --labels <PATH-TO matched_CEM_table.csv> --data <PATH-TO raw_data_.pkl> --embedding <PATH-TO embedding.100d.csv> --model tangle --output <PATH-TO-OUTPUT>`

or, equivalently:

`$ python scripts/cross_validate.py -n N -l <PATH-TO matched_CEM_table.csv> -d <PATH-TO raw_data_.pkl> -e <PATH-TO embedding.100d.csv> -m tangle -o <PATH-TO-OUTPUT>`

Evaluate the average predictive performance of **TANGLE** on N random training/validation/test splits.

### - `cross_validate_linear_model.py`

Usage example:

`$ python scripts/cross_validate_linear_model.py --n_splits N --ngram G --labels <PATH-TO matched_CEM_table.csv> --data <PATH-TO raw_data_.pkl> --output <PATH-TO-OUTPUT>`

or, equivalently:

`$ python scripts/cross_validate_linear_model.py -n N -ng G -l <PATH-TO matched_CEM_table.csv> -d <PATH-TO raw_data_.pkl> -o <PATH-TO-OUTPUT>`

Evaluate the average predictive performance of linear models (x-BOW + Logit) on N random training/validation/test splits.

Embeddings generated by averaging `glove6B` embeddings of the MBS items textual description.

See the `embedding.ipynb` notebook.
# Notebooks

The notebooks in this folder are mainly used to test some of the utility functions of this project.

In order to run these notebooks:
1. you are supposed to have organized the MBS-PBS 10% dataset in the folder `../../../data`;
2. you should run `jupyter notebook` or `jupyter lab` on your `<ROOT>` folder
