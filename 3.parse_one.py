#!/usr/bin/env python

# This is a basic parsing script that will generate output files
# organized / named by repo, along with a "super credit" summary of all repos.

__author__ = "Vanessa Sochat"
__copyright__ = "Copyright 2022, Vanessa Sochat"
__license__ = "MPL 2.0"

from rse.main import Encyclopedia
from rse.utils.command import Command
from rse.utils.file import (
    recursive_find,
    write_json,
    read_json,
)

import citelang.utils as utils
import citelang.main.parser as parser
import citelang.main.packages as packages

import tempfile
import shutil
import argparse
import re
import sys
import os


def clone(url, dest):
    dest = os.path.join(dest, os.path.basename(url))
    cmd = Command("git clone --depth 1 %s %s" % (url, dest))
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
        "repo",
        help="repository UID to parse.",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        help="Output directory for data.",
        default=os.path.join(os.getcwd(), "_repos"),
    )
    return parser


repo_header = """---
title: %s
layout: repo
tipue_search_active: true
exclude_from_search: true
---
"""


def main():

    p = get_parser()
    args, extra = p.parse_known_args()

    # Make sure output directory exists
    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        sys.exit("%s does not exist!" % args.outdir)

    if not args.repo:
        sys.exit("You must provide a repository name in the RSEpedia to update.")
    # Create a base temporary folder to work from
    tempdir = tempfile.mkdtemp()
    pedia = Encyclopedia(args.settings_file)

    repo = pedia.get(args.repo)
    if not repo:
        sys.exit("%s is not known to the RSEPedia")

    destdir = os.path.join(outdir, repo.uid)
    destfile = os.path.join(destdir, "data.json")
    dest = clone(repo.url, tempdir)

    found = None
    for filename in recursive_find(dest, pattern="*"):
        basename = os.path.basename(filename)
        if basename in packages.filesystem_manager_names:
            found = filename
            break

    if not found:
        sys.exit("We could not find a requirements file for %s" % repo.uid)

    reqfile = os.path.join(destdir, os.path.basename(found))
    os.makedirs(os.path.dirname(reqfile))

    # copy found file into folder
    shutil.copyfile(found, reqfile)

    try:                
        cli = parser.RequirementsParser(filename=found, min_credit=0.001)
        result = cli.gen(name=repo.uid, min_credit=0.001)
    except:
        sys.exit("Issue generating result for %s" % repo.uid)

    # If we don't have data (parsing failure) don't include
    if basename == "setup.py" and "pypi" not in cli.data:
        shutil.rmtree(destdir)
        sys.exit("setup.py could not be parsed.")

    utils.write_json(cli.data, destfile)
    destfile = os.path.join(destdir, "README.md")
    utils.write_file(
        repo_header % repo.uid + result.render(start_end_blocks=False), destfile
    )

    # use the root to derive a badge (this is how done internally)
    destfile = os.path.join(destdir, "badge.png")
    os.system(
        "citelang badge %s %s --min-credit 0.001 --outfile %s --force"
        % (repo.uid, found, destfile)
    )

    # We couldn't find a requirements file to parse
    shutil.rmtree(dest)

    # Clean up dest directory if empty
    if os.path.exists(destdir) and len(os.listdir(destdir)) == 0:
        shutil.rmtree(destdir)
    shutil.rmtree(tempdir)
    
if __name__ == "__main__":
    main()
