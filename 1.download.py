#!/usr/bin/env python

# This is a basic parsing script that will generate a tree organized by
# repository, and save a concatenated version of all markdown files found
# in the repository under a nested output directory (defaults to data)
# We can then use that for further processing and analysis.

__author__ = "Vanessa Sochat"
__copyright__ = "Copyright 2022, Vanessa Sochat"
__license__ = "MPL 2.0"

from rse.main import Encyclopedia
from rse.utils.command import Command
from rse.utils.file import recursive_find, read_file, mkdir_p, write_file, write_json
import tempfile
import shutil
import argparse
import re
import sys
import os


def clone(url, dest):
    dest = os.path.join(dest, os.path.basename(url))
    cmd = Command("git clone %s %s" % (url, dest))
    cmd.execute()
    if cmd.returncode != 0:
        print("Issue cloning %s" % url)
        return
    return dest


def get_parser():
    parser = argparse.ArgumentParser(
        description="Research Software Encyclopedia Analyzer",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--settings-file",
        dest="settings_file",
        help="custom path to settings file.",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        help="Output directory for data.",
        default=os.path.join(os.getcwd(), "data"),
    )
    return parser


def main():

    parser = get_parser()
    args, extra = parser.parse_known_args()

    # Make sure output directory exists
    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        sys.exit("%s does not exist!" % args.outdir)

    # Create a base temporary folder to work from
    tempdir = tempfile.mkdtemp()

    pedia = Encyclopedia(args.settings_file)
    repos = list(pedia.list())
    total = len(repos)

    # Keep a master lookup of topics
    topics = {}

    # And languages (to color by later)
    languages = {}

    for i, reponame in enumerate(repos):
        repo = pedia.get(reponame[0])
        topics[repo.uid] = repo.data['data'].get('topics', [])

        # TODO we need to add language!
        languages[repo.uid] = repo.data['data'].get('langauge', 'unknown')

        datadir = os.path.join(outdir, repo.uid)
        destfile = os.path.join(datadir, "CONCAT.md")

        # Don't parse twice!
        if os.path.exists(destfile):
            continue

        dest = clone(repo.url, tempdir)
        if not dest:
            continue

        # Concat markdown for the repository
        text = ""
        for md in recursive_find(dest, "*.md"):
            if re.search("LICENSE", md, re.IGNORECASE):
                continue
            try:
                text += "".join(read_file(md))
            except:
                print("Issue parsing file %s" % md)
                continue

        # Don't include empty repos
        if not text:
            continue
        print("Adding data for %s" % repo.uid)
        mkdir_p(datadir)
        write_file(destfile, text)
        shutil.rmtree(dest)

    shutil.rmtree(tempdir)

    # Save topics to file (in docs so we don't overwhelm github pages)
    write_json(topics, os.path.join('docs','topics.json'))

if __name__ == "__main__":
    main()
