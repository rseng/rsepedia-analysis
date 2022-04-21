#!/usr/bin/env python

# Given a set of corpus.txt, read in each and iteratively train the model.

__author__ = "Vanessa Sochat"
__copyright__ = "Copyright 2022, Vanessa Sochat"
__license__ = "MPL 2.0"


from rse.utils.file import recursive_find, write_file, read_file
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix, lil_matrix
from scipy.spatial.distance import euclidean
from sklearn import datasets, preprocessing, manifold
from river import cluster, stream, feature_extraction
import numpy as np
import collections
import pickle
import pandas
import tempfile
import shutil
import argparse
import itertools
import re
import json
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
        
        # We can use this to eliminate repos that have < unique words
        self.counter = feature_extraction.BagOfWords()

    def __iter__(self):
        removed = 0
        for filename in recursive_find(self.datadir, "corpus.txt"):
            uid = os.path.dirname(filename.replace(self.datadir, "").strip(os.sep))
            text = read_file(filename, readlines=False).strip()
            text = [x for x in text.split(" ") if "_" not in x]

            counts = self.counter.process_text(" ".join(text))
            # Don't include repos that have <200 TOTAL words or <100 unique words
            if len(text) < 200 or len(counts) < 100:
                removed += 1
                continue

            print("%s: %s terms" % (uid, len(text)))
            yield uid, " ".join(text), len(counts)
        print("%s repos skipped <100 terms" % removed)


class LimitedBagOfWords(feature_extraction.BagOfWords):
    """
    A custom class to only include words above some threshold
    """
    min_occur = 5

    def transform_one(self, x):
        counts = collections.Counter(self.process_text(x))

        counts = {k:v for k,v in counts.items() if v> self.min_occur}   
        counts = [[k]*v for k,v in counts.items()]
        counts = list(itertools.chain(*counts))
        return collections.Counter(self.process_text(" ".join(counts)))


def main():

    parser = get_parser()
    args, extra = parser.parse_known_args()

    # Make sure output directory exists
    datadir = os.path.abspath(args.data_dir)
    if not os.path.exists(datadir):
        sys.exit("%s does not exist!" % datadir)

    corpus = SoftwareCorpus(datadir)

    # Previous predictions:
    assigned = {}
    #if os.path.exists("predictions.json"):
    #    with open("predictions.json", "r") as fd:
    #        assigned = json.loads(fd.read())

    # Load a previously existing model, if we have it!
    if os.path.exists("model.pkl"):
        with open("model.pkl", "rb") as fd:  
            model = pickle.load(fd)

    #else:
    model = feature_extraction.TFIDF() | cluster.DBSTREAM(
                clustering_threshold=1,
                fading_factor=0.01,
                cleanup_interval=4,
                intersection_factor=0.3,
                minimum_weight=1,
            )

    # Do we have new data?
    new_data = False

    print("--- TRAINING ---")
    for uid, text, _ in corpus:

        # Only learn from samples we haven't seen before
        if uid not in assigned:
            new_data = True
            # Filter out underscores (variables)
            model = model.learn_one(text)

    # Now get a cluster assignment for each repo
    if new_data:
        print("\n--- GENERATING PREDICTIONS ---")
        for uid, text, unique_words in corpus:
            # Filter out underscores (variables)
            assigned[uid] = {"prediction": model.predict_one(text), "unique_words": unique_words}
        print("%s cluster centers." % len(model.steps["DBSTREAM"].centers))

        with open("predictions.json", "w") as fd:
            fd.write(json.dumps(assigned, indent=4))

        with open("model.pkl", "wb") as fd:
            pickle.dump(model, fd)

    # Save lookup of software repos for each cluster
    meta = {}
    for repo, datum in assigned.items():
        cluster_id = datum['prediction']
        if cluster_id not in meta:
            meta[cluster_id] = {"repos": []}
        meta[cluster_id]['repos'].append(repo + ":" + str(datum['unique_words']))

    # Add sizes for easy access
    for cluster_id, repos in meta.items():
        repos['size'] = len(repos['repos'])

    if new_data:
        # Save metadata to file
        with open(os.path.join("docs", "meta.json"), "w") as fd:
            fd.write(json.dumps(meta, indent=4))

    # Generate centroids
    print("\n--- CENTROIDS ---")
    df = pandas.DataFrame(model.steps["DBSTREAM"].centers)
    df = df.fillna(0)

    # Save dataframe to file
    if new_data:
        df.to_csv("centroids.csv")
    
    # Create a sparse distance matrix
    matrix = lil_matrix((df.shape[1], df.shape[1]), dtype=np.int8)

    # Populate manually...
    print("\n--- DISTANCE ---")
    for idx1 in list(df.columns):
        print("parsing %s of %s" %(idx1, df.shape[1]), end="\r")
        for idx2 in list(df.columns):
            if idx2 < idx1:
                matrix[idx1, idx2] = matrix[idx2, idx1] = euclidean(df.loc[:, idx1], df.loc[:, idx2])

    # Make the tsne!
    print("\n--- DIMENSIONALITY REDUCTION ---")
    fit = manifold.TSNE(n_components=2)
    embedding = fit.fit_transform(matrix)
    emb = pandas.DataFrame(embedding, index=df.columns, columns=["x", "y"])
    emb.index.name = "name"
    emb.to_csv(os.path.join("docs", "software-embeddings.csv"))


if __name__ == "__main__":
    main()
