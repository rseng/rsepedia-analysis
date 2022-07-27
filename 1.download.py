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

    # Create a base temporary folder to work from
    tempdir = tempfile.mkdtemp()

    pedia = Encyclopedia(args.settings_file)
    repos = list(pedia.list())

    # Keep track of number we cannot parse
    missing_requirements = set()
    if os.path.exists("missing-requirements.json"):
        missing_requirements = set(utils.read_json("missing-requirements.json"))

    # Load metadata
    # Keep a master lookup of topics and metadata (a lookup)
    meta = {"topics": {}, "language": {}, "url": {}, "description": {}}
    meta_json = os.path.join("docs", "meta.json")

    if os.path.exists(meta_json):
        meta = read_json(meta_json)

    for i, reponame in enumerate(repos):

        # Don't try parsing again if we know cannot parse!
        if reponame[0] in missing_requirements:
            continue

        repo = pedia.get(reponame[0])
        meta["topics"][repo.uid] = repo.data["data"].get("topics", [])
        meta["language"][repo.uid] = repo.data["data"].get("language", "unknown")
        meta["url"][repo.uid] = repo.url
        meta["description"][repo.uid] = repo.description

        # Don't parse twice!
        destdir = os.path.join(outdir, repo.uid)
        if not os.path.exists(destdir):
            os.makedirs(destdir)

        destfile = os.path.join(destdir, "data.json")
        reqfile = None
        for filename in recursive_find(destdir, pattern="*"):
            basename = os.path.basename(filename)
            if basename in packages.filesystem_manager_names:
                reqfile = filename
                break

        # Skip if we already cloned, etc.
        if os.path.exists(destfile) and reqfile:
            continue

        print(f"{reponame}: {i} of {len(repos)}")
        if not repo.url:
            continue
        dest = clone(repo.url, tempdir)
        if not dest:
            continue

        found = None
        for filename in recursive_find(dest, pattern="*"):
            basename = os.path.basename(filename)
            if basename in packages.filesystem_manager_names:
                found = filename
                break

        # gitlab/octopus-code/octopus
        # we shouldn't be looking up numpy twice
        # Add to global and repo-specific parser
        if found:

            reqfile = os.path.join(destdir, os.path.basename(found))

            # copy found file into folder
            shutil.copyfile(found, reqfile)

            if os.path.exists(destfile):
                continue
            try:
                cli = parser.RequirementsParser(filename=found, min_credit=0.001)
                result = cli.gen(name=repo.uid, min_credit=0.001)
            except:
                continue

            # If we don't have data (parsing failure) don't include
            if basename == "setup.py" and "pypi" not in cli.data:
                shutil.rmtree(destdir)
                missing_requirements.add(repo.uid)
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

        # We couldn't find a requirements file to parse
        else:
            shutil.rmtree(destdir)
            missing_requirements.add(repo.uid)

        shutil.rmtree(dest)

        # Clean up dest directory if empty
        if os.path.exists(destdir) and len(os.listdir(destdir)) == 0:
            missing_requirements.add(repo.uid)
            shutil.rmtree(destdir)

        # Always write missing to file
        utils.write_json(list(missing_requirements), "missing-requirements.json")

        # Save topics, etc. to file (in docs so we don't overwhelm github pages)
        write_json(meta, meta_json)

    shutil.rmtree(tempdir)
    write_json(meta, meta_json)
    utils.write_json(list(missing_requirements), "missing-requirements.json")


if __name__ == "__main__":
    main()
