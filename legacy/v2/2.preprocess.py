#!/usr/bin/env python

# Given a data directory with UID folders that have CONCAT.md (markdown files found
# in the repository concatenated together) preprocess into a format we can
# stream into word2vec (or doc2ve)

__author__ = "Vanessa Sochat"
__copyright__ = "Copyright 2022, Vanessa Sochat"
__license__ = "MPL 2.0"


from rse.utils.file import recursive_find, write_file, read_file
from nltk.corpus import stopwords
from nltk.stem import PorterStemmer
from nltk.tokenize import word_tokenize
import tempfile
import shutil
import argparse
import re
import sys
import os


# Derive stop words and stemmer once
stop_words = set(stopwords.words("english"))
stemmer = PorterStemmer()


def get_parser():
    parser = argparse.ArgumentParser(
        description="Research Software Encyclopedia Preprocessor",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "data_dir",
        help="Directory with UID subfolders",
        default=os.path.join(os.getcwd(), "data"),
    )
    return parser


def process_text(text):
    """
    Process text, including:

    1. Lowercase
    2. Remove numbers and punctuation
    3. Strip whitespace
    4. Tokenize and stop word removal
    5. Stemming
    """
    # Make lowercase
    text = text.lower()

    # Remove numbers and punctuation
    text = re.sub(r"\d+", "", text)
    text = re.sub(r"[^\w\s]", "", text)

    # Strip whitespace
    text = text.strip()

    # tokenize and stop word removal
    tokens = [x for x in word_tokenize(text) if not x in stop_words]

    # Split words with underscore into two words
    words = []
    for t in tokens:
        if "_" in t:
            words += [x.strip() for x in t.split("_")]

        # Don't add single letters
        elif len(t) == 1:
            continue
        else:
            words.append(t)

    # Stemming
    words = [stemmer.stem(t) for t in tokens]
    return " ".join(words)


def main():

    parser = get_parser()
    args, extra = parser.parse_known_args()

    # Make sure output directory exists
    datadir = os.path.abspath(args.data_dir)
    if not os.path.exists(datadir):
        sys.exit("%s does not exist!" % datadir)

    for concat in recursive_find(datadir, "CONCAT.md"):
        outfile = os.path.join(os.path.dirname(concat), "corpus.txt")
        # Don't generate existing again!
        if os.path.exists(outfile):
            continue
        uid = os.path.dirname(concat.replace(datadir, "").strip(os.sep))
        text = read_file(concat, readlines=False)
        text = process_text(text)
        # Output file should be corpus.txt in same folder
        print("Writing %s" % outfile)
        write_file(outfile, text)


if __name__ == "__main__":
    main()
