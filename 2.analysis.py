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
    write_file,
    write_json,
    read_json,
)

import citelang.main.parser as parser
import citelang.main.packages as packages

import copy
import argparse
import re
import sys
import os
import tempfile
import shutil

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


header = """---
title: %s
description: "Most used dependencies across the RSEPedia database"
layout: table
permalink: /analysis/%s/
tipue_search_active: true
exclude_from_search: true
---
"""


repo_header = """---
title: %s
layout: repo
tipue_search_active: true
exclude_from_search: true
---
"""


def count_repos(data):
    count = 0
    for _, repos in data.items():
        count += len(repos)
    return count


def remove_emojis(data):
    """
    https://stackoverflow.com/questions/33404752/removing-emojis-from-a-string-in-python
    """
    emoj = re.compile(
        "["
        "\U0001F600-\U0001F64F"  # emoticons
        "\U0001F300-\U0001F5FF"  # symbols & pictographs
        "\U0001F680-\U0001F6FF"  # transport & map symbols
        "\U0001F1E0-\U0001F1FF"  # flags (iOS)
        "\U00002500-\U00002BEF"  # chinese char
        "\U00002702-\U000027B0"
        "\U00002702-\U000027B0"
        "\U000024C2-\U0001F251"
        "\U0001f926-\U0001f937"
        "\U00010000-\U0010ffff"
        "\u2640-\u2642"
        "\u2600-\u2B55"
        "\u200d"
        "\u23cf"
        "\u23e9"
        "\u231a"
        "\ufe0f"  # dingbats
        "\u3030"
        "]+",
        re.UNICODE,
    )
    return re.sub(emoj, "", data)


def main():

    p = get_parser()
    args, extra = p.parse_known_args()

    # Make sure output directory exists
    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        sys.exit("%s does not exist!" % args.outdir)

    pedia = Encyclopedia(args.settings_file)
    repos = list(pedia.list())

    # This is a global client to track ALL credit scores
    global_cli = parser.RequirementsParser()

    # Load metadata
    # Keep a master lookup of topics and metadata (a lookup)
    meta_json = os.path.join("docs", "meta.json")
    meta = read_json(meta_json)
    language_json = os.path.join("docs", "language-counts.json")

    count = 0
    for i, reponame in enumerate(repos):
        repo = pedia.get(reponame[0])

        destdir = os.path.join(outdir, repo.uid)
        destfile = os.path.join(destdir, "data.json")
        if not os.path.exists(destfile):
            continue

        print(f"{reponame}: {i} of {len(repos)}")
        count += 1

        # Find the metadata file in the repo
        found = None
        for filename in recursive_find(destdir, pattern="*"):
            basename = os.path.basename(filename)
            if basename in packages.filesystem_manager_names:
                found = filename
                break

        # This should not happen
        if not found:
            tempdir = tempfile.mkdtemp()
            dest = clone(repo.url, tempdir)
            for filename in recursive_find(tempdir, pattern="*"):
                basename = os.path.basename(filename)
                if basename in packages.filesystem_manager_names:
                    found = filename
                    break

            if not found:
                print('WARNING %s not found!' % repo.uid)
                continue

            # copy found file into folder
            reqfile = os.path.join(destdir, os.path.basename(found))
            shutil.copyfile(found, reqfile)
            shutil.rmtree(dest)

        if not found:
            print('WARNING %s not found!' % repo.uid)
            continue

        try:
            global_cli.gen(repo.uid, filename=found, min_credit=0.001)
        except:
            continue

    write_json(meta, meta_json)

    # Write language counts
    counts = {}
    for repo, language in meta["language"].items():
        if language not in counts:
            counts[language] = 0
        counts[language] += 1
    write_json(counts, language_json)

    # Write results, changing round to include most
    global_cli.round_by = 100
    content = (
        header
        % (
            "RSEPedia Top Dependencies",
            "dependencies",
        )
        + global_cli.render(start_end_blocks=False)
    )

    # Replace for jekyll site
    write_file(
        os.path.join("docs", "all-repos.md"), global_cli.render(start_end_blocks=False)
    )
    write_file(os.path.join("pages", "dependencies.md"), content)

    # Get stats for different languages
    # 'setup.py', 'package.json', 'npm', 'DESCRIPTION', 'cran', 'pypi', 'go.mod', 'go', 'requirements.txt'])
    python_deps = len(
        set(global_cli.data["pypi"])
        .union(global_cli.data["requirements.txt"])
        .union(global_cli.data["pypi"])
    )
    r_deps = len(set(global_cli.data["DESCRIPTION"]).union(global_cli.data["cran"]))
    js_deps = len(set(global_cli.data["npm"]).union(global_cli.data["package.json"]))
    go_deps = len(set(global_cli.data["go"]).union(global_cli.data["go.mod"]))

    stats = {
        "python_deps": python_deps,
        "r_deps": r_deps,
        "js_deps": js_deps,
        "go_deps": go_deps,
        "total_repos": len(repos),
        "total_parsed": count,
    }
    write_json(stats, os.path.join("_data", "stats.json"))
    write_json(counts, os.path.join("_data", "language_counts.json"))

    # Prepare scoped tables to languages
    custom_cli = copy.deepcopy(global_cli)

    # Keep count of repos / deps files for each
    repos_counts = {}

    # Python
    custom_cli.data = custom_cli.prepare_custom_table(
        ["setup.py", "requirements.txt", "pypi"]
    )
    content = (
        header
        % (
            "RSEPedia Top Python Dependencies",
            "python",
        )
        + custom_cli.render(start_end_blocks=False)
    )
    write_file(os.path.join("pages", "python.md"), content)
    repos_counts["Python"] = count_repos(custom_cli.data)

    # R
    custom_cli.data = custom_cli.prepare_custom_table(["cran", "DESCRIPTION"])
    content = header % ("RSEPedia Top R Dependencies", "R") + custom_cli.render(
        start_end_blocks=False
    )
    write_file(os.path.join("pages", "r.md"), content)
    repos_counts["R"] = count_repos(custom_cli.data)

    # Javascript
    custom_cli.data = custom_cli.prepare_custom_table(["package.json", "npm"])
    content = header % ("RSEPedia Top Js Dependencies", "js") + custom_cli.render(
        start_end_blocks=False
    )
    write_file(os.path.join("pages", "js.md"), content)
    repos_counts["Js"] = count_repos(custom_cli.data)

    # Go
    custom_cli.data = custom_cli.prepare_custom_table(["go.mod", "go"])
    content = header % ("RSEPedia Top Go Dependencies", "go") + custom_cli.render(
         start_end_blocks=False
     )
    write_file(os.path.join("pages", "go.md"), content)
    repos_counts["Go"] = count_repos(custom_cli.data)

    # Emoji are problematic for jekyll data
    for repo, desc in meta["description"].items():
        if not desc:
            continue
        meta["description"][repo] = remove_emojis(desc)
    write_json(meta, os.path.join("_data", "repos.json"))
    write_json(repos_counts, os.path.join("_data", "repos_counts_languages.json"))


if __name__ == "__main__":
    main()
