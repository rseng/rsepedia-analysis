#!/usr/bin/env python

# This is a basic parsing script that will generate output files
# organized / named by repo, along with a "super credit" summary of all repos.

__author__ = "Vanessa Sochat"
__copyright__ = "Copyright 2022, Vanessa Sochat"
__license__ = "MPL 2.0"

from rse.main import Encyclopedia
from rse.utils.command import Command
from rse.utils.file import (
    write_file,
    write_json,
    read_json,
)

import citelang.main.parser as parser

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
---
"""


repo_header = """---
title: %s
layout: repo
tipue_search_active: true
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

    # Collect list of data files
    data_files = []
    count = 0
    for i, reponame in enumerate(repos):
        repo = pedia.get(reponame[0])

        destdir = os.path.join(outdir, repo.uid)
        destfile = os.path.join(destdir, "data.json")
        if not os.path.exists(destfile):
            continue
        print(f"{reponame}: {i} of {len(repos)}")
        data_files.append(destfile)
        count += 1

    # This is how to render the custom (loaded) data
    roots = global_cli.load_datafiles(data_files)

    # Write results, changing round to include most
    global_cli.round_by = 100
    content = header % (
        "RSEPedia Top Dependencies",
        "dependencies",
    ) + global_cli.render(start_end_blocks=False, data=roots)
    write_json(meta, meta_json)

    # Write language counts
    counts = {}
    for repo, language in meta["language"].items():
        if language not in counts:
            counts[language] = 0
        counts[language] += 1
    write_json(counts, language_json)

    # Replace for jekyll site
    write_file(
        global_cli.render(start_end_blocks=False, data=roots),
        os.path.join("docs", "all-repos.md"),
    )
    write_file(content, os.path.join("pages", "dependencies.md"))

    # Get stats for different languages
    # 'setup.py', 'package.json', 'npm', 'DESCRIPTION', 'cran', 'pypi', 'go.mod', 'go', 'requirements.txt'])
    python_deps = len(
        set(roots["pypi"]).union(roots["requirements.txt"]).union(roots["pypi"])
    )
    cpp_deps = len(set(roots["spack"]))
    r_deps = len(set(roots["DESCRIPTION"]).union(roots["cran"]))
    js_deps = len(set(roots["npm"]).union(roots["package.json"]))
    go_deps = len(set(roots["go"]).union(roots["go.mod"]))

    stats = {
        "python_deps": python_deps,
        "cpp_deps": cpp_deps,
        "r_deps": r_deps,
        "js_deps": js_deps,
        "go_deps": go_deps,
        "total_repos": len(repos),
        "total_parsed": count,
    }
    write_json(stats, os.path.join("_data", "stats.json"))
    write_json(counts, os.path.join("_data", "language_counts.json"))

    # Prepare scoped tables to languages
    # Keep count of repos / deps files for each
    repos_counts = {}

    # Python
    data = global_cli.load_datafiles(
        data_files, includes=["setup.py", "requirements.txt", "pypi"]
    )
    content = header % (
        "RSEPedia Top Python Dependencies",
        "python",
    ) + global_cli.render(start_end_blocks=False, data=data)
    write_file(content, os.path.join("pages", "python.md"))
    repos_counts["Python"] = count_repos(data)

    # R
    data = global_cli.load_datafiles(data_files, includes=["cran", "DESCRIPTION"])
    content = header % ("RSEPedia Top R Dependencies", "R") + global_cli.render(
        start_end_blocks=False, data=data
    )
    write_file(content, os.path.join("pages", "r.md"))
    repos_counts["R"] = count_repos(data)

    # Cpp
    data = global_cli.load_datafiles(data_files, includes=["spack"])
    content = header % (
        "RSEPedia Top Spack (C++) Dependencies",
        "cpp",
    ) + global_cli.render(start_end_blocks=False, data=data)
    write_file(content, os.path.join("pages", "cpp.md"))
    repos_counts["Cpp"] = count_repos(data)

    # Javascript
    data = global_cli.load_datafiles(data_files, includes=["package.json", "npm"])
    content = header % ("RSEPedia Top Js Dependencies", "js") + global_cli.render(
        start_end_blocks=False, data=data
    )
    write_file(content, os.path.join("pages", "js.md"))
    repos_counts["Js"] = count_repos(data)

    # Go
    data = global_cli.load_datafiles(data_files, includes=["go.mod", "go"])
    content = header % ("RSEPedia Top Go Dependencies", "go") + global_cli.render(
        start_end_blocks=False, data=data
    )
    write_file(content, os.path.join("pages", "go.md"))
    repos_counts["Go"] = count_repos(data)

    # Emoji are problematic for jekyll data
    for repo, desc in meta["description"].items():
        if not desc:
            continue
        meta["description"][repo] = remove_emojis(desc)
    write_json(meta, os.path.join("_data", "repos.json"))
    write_json(repos_counts, os.path.join("_data", "repos_counts_languages.json"))


if __name__ == "__main__":
    main()
