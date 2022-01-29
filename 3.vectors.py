#!/usr/bin/env python

# Given a set of corpus.txt, read in each and iteratively train the model.

__author__ = "Vanessa Sochat"
__copyright__ = "Copyright 2022, Vanessa Sochat"
__license__ = "MPL 2.0"


from rse.utils.file import recursive_find, write_file, read_file
from gensim.models.doc2vec import Doc2Vec, TaggedDocument
from scipy.spatial.distance import pdist, squareform
from sklearn import datasets, preprocessing, manifold
import pandas
import tempfile
import shutil
import argparse
import re
import sys
import os


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


def main():

    parser = get_parser()
    args, extra = parser.parse_known_args()

    # Make sure output directory exists
    datadir = os.path.abspath(args.data_dir)
    if not os.path.exists(datadir):
        sys.exit("%s does not exist!" % datadir)

    corpus = SoftwareCorpus(datadir)

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
