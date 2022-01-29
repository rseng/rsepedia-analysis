# Research Software Encyclopedia Analysis

We ultimately want to be able to group the different software into categories.
E.g., if I'm an astronomy researcher I want to quickly see what other libraries
are being developed. To do this, we can use [word2vec](https://radimrehurek.com/gensim/models/word2vec.html) on text from the repos
to derive a set of embeddings that represent each code repository.
This means that each repository is going to be represented as a vector of numbers.
We then want to take the label data (where we have topics on the repo) and we will
train a model to derive the probability of each embedding (without labels) of having
the labels. We can try a "separate label" model, or a multiple label model.
On a higher level, once we have the embeddings the topics can be used to train a probabilistic model.

## Analysis Steps

1. For each repository in data, likely we will want to clone and combine all of the markdown or rst files found into a single document.
2. Clean up text, etc., and use word2vec to make a vector.
3. Cluster the vectors to minimally look at similarity (add topics known to see patterns).
4. Then try training the probabilistic models.


## Usage

### 0. Install

Install dependencies (ideally in a virtual environment)

```bash
$ python -m venv env
$ source env/bin/activate
$ pip install -r requirements.txt
$ python -m nltk.downloader popular
```

This will install the research software encyclopedia and Gensim, and nltk
and data we need.

### 1. Data from Repositories

This next step will derive text from each repo, and a master lookup of topics
to use later. If you haven't yet, clone the software repository (or use your own)
somewhere else. E.g.,:

```bash
$ cd ../
$ git clone https://github.com/rseng/software
$ cd rsepedia-analysis
```

Then run the analysis script, targeting the correct rse.ini settings file
for the [software](https://github.com/rseng/software) respository we
just cloned.

```bash
$ mkdir -p data/
$ python 1.download.py --settings-file ../software/rse.ini -o ./data
```

### 2. Preprocess text

Then you'll want to prepare the gensim word2vec corpora. Each subfolder in
data is a unique identifier for a repository, and after this we will generate
a space separated `corpus.txt` in each subfolder.

```bash
$ python 2.preprocess.py ./data
```

**TODO** need to add languages for color!

### 3. Model and Vectors

And then train the model!

```bash
$ python 3.vectors.py ./data
```

This will generate vectors along with embeddings and the distance matrix for
those embeddings that drive the visualization in `index.html`. The plot shows
the different repository embeddings, colored by language.

1. fix jumping of UI 
2. Color by language

### 4. Probabilistic Model

The labels are too distinct I think to be useful, so instead I'm going to try:

1. Generating an embedding for each word across the model
2. For each vector (document) find the K nearest neighbors (KNN) to derive a set of words
3. Associate the words and add to the plot!

**under development**
