#!/usr/bin/env python

# This is a basic parsing script that will generate output files
# organized / named by repo, along with a "super credit" summary of all repos.

__author__ = "Vanessa Sochat"
__copyright__ = "Copyright 2022, Vanessa Sochat"
__license__ = "MPL 2.0"

from rse.main import Encyclopedia
from rse.utils.command import Command
from rse.utils.file import recursive_find

import citelang.utils as utils
import citelang.main.parser as parser
import citelang.main.packages as packages

import tempfile
import shutil
import argparse
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
        "--language-file",
        dest="language_files",
        action="append",
        help="Language requirements file basename to re-parse.",
        choices=[
            "DESCRIPTION",
            "requirements.txt",
            "setup.py",
            "go.mod",
            "package.json",
        ],
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

    # Create a base temporary folder to work from
    tempdir = tempfile.mkdtemp()

    # Make sure output directory exists
    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        sys.exit("%s does not exist!" % args.outdir)

    # Create a base temporary folder to work from
    tempdir = tempfile.mkdtemp()

    pedia = Encyclopedia(args.settings_file)
    repos = list(pedia.list())
    for i, reponame in enumerate(repos):

        repo = pedia.get(reponame[0])
        destdir = os.path.join(outdir, repo.uid)
        if not os.path.exists(destdir):
            continue

        destfile = os.path.join(destdir, "data.json")
        found = False
        for lang_file in args.language_files:
            reqfile = os.path.join(destdir, lang_file)
            if os.path.exists(reqfile):
                found = True
                break

        if not found:
            continue

        # Re-clone, because we need reqs files in context of their code!
        dest = clone(repo.url, tempdir)
        if not dest:
            continue

        # By here we've already filtered to requirements files
        found = None
        for filename in recursive_find(dest, pattern="*"):
            basename = os.path.basename(filename)
            if basename in packages.filesystem_manager_names:
                found = filename
                break

        if not found:
            continue

        reqfile = os.path.join(destdir, os.path.basename(found))

        # copy updated found file into folder
        shutil.copyfile(found, reqfile)

        # Files that are weirdly binaries are going to error here
        try:
            cli = parser.RequirementsParser(filename=found, min_credit=0.001)
            result = cli.gen(name=repo.uid, min_credit=0.001)
        except UnicodeDecodeError:
            print("Issue parsing %s because of unicode, skipping." % repo.uid)
            continue

        continue
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
        shutil.rmtree(dest)

    # Clean up
    shutil.rmtree(tempdir)


if __name__ == "__main__":
    main()
