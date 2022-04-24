#!/usr/bin/env python

# Given a set of corpus.txt, generate vectors based on Wikipedia topics

__author__ = "Vanessa Sochat"
__copyright__ = "Copyright 2022, Vanessa Sochat"
__license__ = "MPL 2.0"


from rse.utils.file import recursive_find, write_file, read_file

from gensim.models.doc2vec import Doc2Vec, TaggedDocument
from gensim.corpora.wikicorpus import WikiCorpus

from scipy.spatial.distance import pdist, squareform
from sklearn import datasets, preprocessing, manifold
import multiprocessing
from pprint import pprint
import pandas
import tempfile
import shutil
import argparse
import re
import sys
import io
import os



def load_vectors(fname):
    fin = io.open(fname, 'r', encoding='utf-8', newline='\n', errors='ignore')
    n, d = map(int, fin.readline().split())
    data = {}
    for line in fin:
        tokens = line.rstrip().split(' ')
        data[tokens[0]] = map(float, tokens[1:])
    return data


def get_parser():
    parser = argparse.ArgumentParser(
        description="Research Software Encyclopedia Vectorizer",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "data_dir",
        help="Directory with UID subfolders",
        default=os.path.join(os.getcwd(), "data"),
    )
    return parser


class SoftwareCorpus:
    def __init__(self, datadir):
        self.datadir = datadir

    def iter_text(self):
        for filename in recursive_find(self.datadir, "corpus.txt"):
            uid = os.path.dirname(filename.replace(self.datadir, "").strip(os.sep))
            yield uid, read_file(filename, readlines=False).split()

    def __iter__(self):
        for uid, words in self.iter_text():
            print("Training with %s" % uid)
            yield TaggedDocument(words, [uid])


class WikipediaCorpus(object):
    def __init__(self, wiki):
        self.wiki = wiki
        self.wiki.metadata = True

    def __iter__(self):
        for content, (page_id, title) in self.wiki.get_texts():
            yield TaggedDocument([c for c in content], [title])


def main():

    parser = get_parser()
    args, extra = parser.parse_known_args()

    # Make sure output directory exists
    datadir = os.path.abspath(args.data_dir)
    if not os.path.exists(datadir):
        sys.exit("%s does not exist!" % datadir)

    wiki_file = "wiki.en.vec"
    if not os.path.exists(wiki_file):
        sys.exit('%s does not exist, see instructions in README to download and extract.' % wiki_file)
    vectors = load_vectors(wiki_file)

    import IPython
    IPython.embed()
    sys.exit()
    pre_model = 'model.pre.doc2vec'
    pickled_docs = 'documents.pkl'
    if os.path.exists(pre_model) and os.path.exists(pickled_docs):
        print('LOAD MODELS')
        import IPython
        IPython.embed()
        sys.exit()

    else:
        corpus = SoftwareCorpus(datadir)

        wikifile = "enwiki-latest-pages-articles.xml.bz2"
        if not os.path.exists(wikifile):
            sys.exit("You need to download %s first, see README in repository!" % wikifile)
        wiki = WikiCorpus(wikifile)
        documents = WikipediaCorpus(wiki)
        with open(pickle_docs, 'wb') as f:
            pickle.dump(documents, f) 

        # Preprocessing (this also takes a long time)!
        pre = Doc2Vec(min_count=0)
        pre.scan_vocab(documents)
        pre.save(pre_model)
    
    # Calculate the optimal min_count parameter to calculate the vocabulary size
    # NOTE pre didn't have scale_vocab method so we will use 19/20
    cores = multiprocessing.cpu_count()
    models = [
        # PV-DBOW
        Doc2Vec(
            dm=0, dbow_words=1, vector_size=200, window=8, min_count=19, epochs=10, workers=cores
        ),
        # PV-DM w/average
        Doc2Vec(
            dm=1, dm_mean=1, vector_size=200, window=8, min_count=19, epochs=10, workers=cores
        ),
    ]

    import IPython
    IPython.embed()
    sys.exit(0)
    models[0].build_vocab(documents)
    print(str(models[0]))
    models[1].reset_from(models[0])
    print(str(models[1]))
    Doc2Vec(dbow + w, d200, hs, w8, mc19, t8)

    # Train doc2vec on wikipedia
    for model in models:
        model.train(documents)

    # 40 epochs means we do it 40 times
    model = Doc2Vec(corpus, vector_size=50, min_count=100, workers=4, epochs=40)

    # Save the model if we need again
    model.save("model.doc2vec")

    # Create a vector for each document
    # UIDS as id for each row, vectors across columns
    df = pandas.DataFrame(columns=range(50))

    print("Generating vector matrix for documents...")
    for uid, words in corpus.iter_text():
        df.loc[uid] = model.infer_vector(words)

    # Save dataframe to file
    df.to_csv("vectors.csv")

    # Create a distance matrix
    distance = pandas.DataFrame(
        squareform(pdist(df)), index=list(df.index), columns=list(df.index)
    )
    distance.to_csv("software-distances.csv")

    # Make the tsne!
    fit = manifold.TSNE(n_components=2)
    embedding = fit.fit_transform(distance)
    emb = pandas.DataFrame(embedding, index=distance.index, columns=["x", "y"])
    emb.index.name = "name"
    emb.to_csv(os.path.join("docs", "software-embeddings.csv"))


if __name__ == "__main__":
    main()
