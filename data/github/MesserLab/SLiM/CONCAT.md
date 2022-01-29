<p align="center">
	<img alt="Screenshot of SLiMgui running on OS X." height="45%" width="45%" src="https://messerlab.files.wordpress.com/2021/12/slimgui_screenshot.jpg"/>
</p>



<p align="center">
	SLiM: Selection on Linked Mutations
</p>
<p align="justify">
	SLiM is an evolutionary simulation framework that combines a powerful engine for population genetic simulations with the capability of modeling arbitrarily complex evolutionary scenarios. Simulations are configured via the integrated Eidos scripting language that allows interactive control over practically every aspect of the simulated evolutionary scenarios. The underlying individual-based simulation engine is highly optimized to enable modeling of entire chromosomes in large populations. We also provide a graphical user interface on macOS and Linux for easy simulation set-up, interactive runtime control, and dynamical visualization of simulation output.
</p>

GitHub Actions | Travis CI | Fedora Copr
---|---|---
![SLiM on GitHub Actions:](https://github.com/MesserLab/SLiM/workflows/tests/badge.svg) |![SLiM on Travis-CI:](https://travis-ci.com/MesserLab/SLiM.svg?branch=master) | [![Copr build status](https://copr.fedorainfracloud.org/coprs/bacarson/SLiM-Selection_on_Linked_Mutations/package/SLiM/status_image/last_build.png)](https://copr.fedorainfracloud.org/coprs/bacarson/SLiM-Selection_on_Linked_Mutations/package/SLiM/)

:construction: This GitHub repository hosts the <em>upstream, development head version</em> of SLiM and SLiMgui.

:warning: <strong>End users should generally not use these sources; they may contain serious bugs, or may not even compile</strong>.

:heavy_check_mark: The <strong><em>release</em></strong> version of SLiM and SLiMgui is available at [http://messerlab.org/slim/](http://messerlab.org/slim/).


License
----------

Copyright (c) 2016-2021 Philipp Messer.  All rights reserved.

SLiM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

SLiM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with SLiM.  If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).

Development & Feedback
-----------------------------------

SLiM is under active development, and our goal is to make it as broadly useful as possible.  If you have feedback or feature requests, or if you are interested in contributing to SLiM, please contact Philipp Messer at [messer@cornell.edu](mailto:messer@cornell.edu). Please note that Philipp is also looking for graduate students and postdocs.

Installation
------------
<em>Looking for Binary Packages / Installers?</em>
OS X | Red Hat Enterprise, CentOS, and Fedora | Debian and Ubuntu | Arch | Windows 10 (WSL & GWSL)
---|---|---|---|---
[http://messerlab.org/slim/](http://messerlab.org/slim/) | [Copr Repository](https://copr.fedorainfracloud.org/coprs/bacarson/SLiM-Selection_on_Linked_Mutations/) | [SLiM-Extras Repository](https://github.com/MesserLab/SLiM-Extras/blob/master/installation/DebianUbuntuInstall.sh) | [Arch User Repository](https://aur.archlinux.org/packages/slim-simulator/) | [SLiM-Extras Repository](https://github.com/MesserLab/SLiM-Extras/blob/master/installation/Windows10Installation.md)


Compilation of SLiM from source
----------------------------------

See chapter two of the SLiM manual for more information about building and installing, including instructions on building SLiMgui (the graphical modeling environment for SLiM) on various platforms.  The manual and other SLiM resources can be found at [http://messerlab.org/slim/](http://messerlab.org/slim/).
---
name: Issue
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

Hi!  Please don't ask questions about how to use SLiM on GitHub; instead, please go to the slim-discuss mailing list at https://groups.google.com/g/slim-discuss.  If you have a bug report or feature request, on the other hand, you have found the right place.  Please provide any information that will allow us to reproduce the bug, including a complete SLiM script (if appropriate) that is self-contained, with all defined constants set up in script rather than needing to be passed on the command line.  If the script requires any external files to run, please attach those too, or provide a link to where they can be downloaded.

Incidentally, if you are a new SLiM user and have not yet done the SLiM Workshop, it is highly recommended to get through the initial learning curve.  It is available for free online at http://benhaller.com/workshops/workshops.html.
# Notes on implementation

## Individual tracking

Both individuals and nodes (= genomes) have flags. We use individual flags for
internal bookkeeping in SLiM; while msprime uses the SAMPLE flag of nodes.

1. When an individual is created in the simulation, its genomes (nodes) are automatically
    added to the node table with the `NODE_IS_SAMPLE` flag, because the SAMPLE flag means
    that we have whole-genome information about the node, which at this point we do.
    Node times are recorded as the number of generations before the start of the
    simulation (i.e., -1 * sim.generation). However, the individuals themselves are *not*
    automatically added to the individuals table, in part because much of the information
    in that table (e.g., location, age) may change, so the user should have control over
    when it is recorded. Individuals are only added to the table during output (when all
    the currently alive individuals are added) or if we explictly choose to call the
    RememberIndividuals() function.

2. RememberIndividuals is called with a list of Individual instances. We add each of
    those to to the individual table (but must first check if they are already present).
    If `permanent=T` (default) we remember them permanently, whereas if `permanent=F` we
    simply flag them to be *retained* while their nodes exist in the simulation. When we
    ask an individual to be permanently remembered it has the `INDIVIDUAL_REMEMBERED`
    flag set in the individuals table, and its genomes' node IDs are added to the
    `remembered_genomes_` vector. In constrast, an individual which is retained has the
    `INDIVIDUAL_RETAINED` flag set (which allows round-tripping on input/output).
    However, unless it has also been permanently remembered, a retained individual is
    subject to removal via simplification.

3. Simplify, if it occurs, first constructs the list of nodes to maintain as samples
    by starting with `remembered_genomes_`, and then adding the genomes of the current
    generation.  All permanently remembered individuals in the individual table should be
    maintained, because they are all pointed to by nodes in this list; other individuals
    which have been retained, but not permanently remembered, may get removed by the
    normal simplification algorithm.  We also update the `msp_node_id_` property of the
    currently alive SLiM genomes. After simplify, we reset the `NODE_IS_SAMPLE` flag in
    the node table, leaving it set only on those nodes listed in `remembered_genomes_`.
    Note that we call `simplify` with the `keep_input_roots` option to retain
    lineages back to the first generation, and the `filter_individuals` option to remove
    individuals that no longer need to be retained.

4. On output, we (optionally) simplify, copy the tables, and with these tables:
    add the current generation to the individual table, marked with the `INDIVIDUAL_ALIVE` flag;
    and reorder the individual tables so that the the currently alive ones come first,
    and in the order they currently exist in the population.
    Also, node times have the current generation added to them, so they are in units
    of number of generations before the end of the simulation.
    Note that if simplify is *not* performed before output, all the nodes since the last
    simplify will be marked with the `NODE_IS_SAMPLE` flag, even those not alive. Many of
    these may not have an associated individual. This is deliberate, as we have full
    genomic information for them, even if they are not ALIVE.

5. On input from tables, we use the `INDIVIDUAL_ALIVE` flag to decide who is currently
    alive and therefore must be instantiated into a SLiM individual, and use the
    `INDIVIDUAL_REMEMBERED` flag to decide which genomes to add to `remembered_genomes_`.
    From the individual table, we remove all individuals apart from those that are marked
    as either `INDIVIDUAL_REMEMBERED` or `INDIVIDUAL_RETAINED` (thus we remove from the
    table individuals that were only exported in the first place because they were ALIVE
    at the time of export). `msp_node_id_` is set based on the nodes that point to that
    individual in the tables. Also, node times have the current generation subtracted
    from them, so they are in units of number of generations before the start of the simulation.

Internal state:

- `tables.individuals`:

    * who's in the table: REMEMBERED and (if they have an existing node) RETAINED
    * flags: everyone either REMEMBERED or RETAINED or both (never ALIVE)

- `tables.nodes`:

    * NODE_IS_SAMPLE flag: **unused** by SLiM; but set if they are in
      REMEMBERED individuals, or if they're from a generation since the
      last simplify (or the current generation, if that just happened)

    * individual column: if their individual is REMEMBERED or RETAINED

        - created by AddIndividualsToTable, called by RememberIndividuals

    * time: (-1) * birth generation.

- `remembered_genomes_`: the node IDs of genomes whose individuals are REMEMBERED

    * initialized at the start of each new, not-from-anywhere-else subpopulation
    * added to by RememberIndividuals
    * updated by Simplify

- `msp_node_id_` properties of currently alive SLiM individual `genome1_` and `genome2_`:
    always contains an ID of the node table.

    * created when a new individual is created, by `node_table_add_row()`
    * updated by Simplify

State at output:

- `tables.individuals`:

    * who's in the table: ALIVE and/or REMEMBERED and/or (if node exists) RETAINED
    * flags: everyone is ALIVE and/or REMEMBERED and/or RETAINED
    * reordered by ReorderIndividualTable so that order matches currently alive ones

- `tables.nodes`:

    * NODE_IS_SAMPLE flag: if they are REMEMBERED and/or ALIVE, or all nodes since the
       last simplification if output has been done with simplify=F
    * individual column: if they are REMEMBERED and/or ALIVE
    * time: (current generation - birth generation)


## Mutation recording

### Alleles

We record as the `derived_state` for each mutation the
*concatenation* of all mutations at that site that the individual has.
This is necessary because stacking rules can change dynamically,
and makes sense, because this is what the individual actually passes on to offspring.

### Sites and mutation parents

Whenever a new mutation is encountered, we do the following:

1. Add a new `site` to the site table at this position.
2. Add a new `mutation` to the mutation table at the newly created site.

This is lazy and wrong, because:

a. There might have already been sites in the site table with the same position,
b. and/or a mutation (at the same position) that this mutation should record as it's `parent`.

But, it's all OK because here's what we do:

1. Add rows to the mutation and site tables as described above.
2. Periodically, `sort_tables`, `deduplicate_sites`,  and `simplify_tables`, then return to (1.), except that
3. Sometimes, to output the tables, we `sort_tables`, `compute_parents`,
    (optionally `simplify_tables`), and dump these out to a file.

*Note:* as things are going along we do *not* have to `compute_parents`;
this only has to happen before the final (output) step.
It might also be possible to skip the `deduplicate_sites` as we go along,
if `simplify_tables` is OK with more than one site at the same position.

And, these operations have the following properties:

1. Mutations appear in the mutation table ordered by time of appearance,
    so that a mutation will always appear after the one that it replaced (i.e., it's parent).
2. Furthermore, if mutation B appears after mutation A, but at the same site, 
    then mutation B's site will appear after mutatation A's site in the site table.
3. `sort_tables` sorts sites by position, and then by ID, so that the relative ordering of sites
    at the same position is maintained, thus preserving property (2).
4. `sort_tables` sorts mutations by site, and then by ID, thus preserving property (1);
    if the mutations are at separate sites (but the same position),
    this fact is thanks to property (2).
5. `simplify_tables` also preserves ordering of any rows in the site and mutation tables
    that do not get discarded.
5. `deduplicate_sites` goes through and collapses all sites at the same position to only one site,
    maintaining order otherwise.
6.  `compute_parents` fills in the `mutation.parent` information by using property (1).


## Metadata schemas

tskit provides methods for structured metadata decoding using [JSON schemas](https://json-schema.org/understanding-json-schema/index.html),
to document what the metadata means.
We don't make use of these, but write them to the tree sequence for tskit use.
There's both top-level metadata (ie for the whole tree sequence)
and metadata for every row in every table.
(But, we only use some of these.)

For practical purposes, the metadata schema can be any equivalent JSON.
However, when checking for table equality,
the code checks whether the underlying string representations are equal.
So, we want the metadata schema as written by SLiM to match that written by pyslim (using tskit).
The way that tskit writes out a metadata schema to a tree sequence
is by (a) defining the schema using a dict
and (b) creating a string using `json.dumps(schema_dict, sort_keys=True, indent=4)`.
Happily, the resulting text matches the output of nlohmann::json dump(4),
so - it seems - we can merrily write out JSON ourselves.
Nonetheless, we only actually do JSON parsing and writing in SLiM's code for top-level metadata:
for all the schemas (including the top-level metadata schema)
we just write out the string representation, as output by tskit, saved in `slim_globals.h`.

In the future we may want to *keep* whatever top-level metadata there is
in the tree sequence already (and the associated keys in the schema).
We've not done that yet because making things exactly match seems like a pain,
and no-one else is using the top-level metadata yet.

### Top-level metadata:

Here's an example of the top-level metadata:
```
{
 "SLiM" : {
     "model_type" : "WF",
     "generation" : 123,
     "file_version" : "0.5",
     "spatial_dimensionality" : "xy",
     "spatial_periodicity" : "x",
     "separate_sexes" : true,
     "nucleotide_based" : false
 }
}
```

However, we're currently only using `model_type`, `generation`, and `file_version`.
# Test suite for tree sequence output from SLiM

## Running the tests
Just do `python3 -m pytest` from within this directory to run the tests on the results of SLiM scripts listed in the `testRecipes/` directory. Alternatively, to run just one of the recipes do e.g.
`python3 -m pytest -k test_000_sexual_nonwf`

## Recaching in GitHub Actions to get a new tskit/msprime version
GitHub Actions caches its install of `tskit`, `msprime`, and other software.  When a new version of such software is released, a recache needs to be forced or these tests will likely fail in CI.  This cannot presently be gone in GitHub's UI; see [this GitHub issue](https://github.com/actions/cache/issues/2).  So to trigger a recache, you need to increment the cache version number.  It is found in `.github/workflows/tests.yml` in the line:

`key: ${{ runner.os }}-${{ matrix.python}}-conda-v10-${{ hashFiles('treerec/tests/conda-requirements.txt') }}-${{ hashFiles('treerec/tests/pip-requirements.txt') }}`

For example, change `v10` here to `v11`.  You can also consider putting an explicit version number in `treerec/tests/conda-requirements.txt`; for example, we now have `tskit >= 0.4.0` there, although presumably it would get the most recent shipping version anyway; I guess this expresses the version requirement semantically, for human readers, at least.  Changing either the version number or the `conda-requirements.txt` should apparently trigger a recache.  See [this issue](https://github.com/MesserLab/SLiM/issues/232) for the context that led to this comment.

## Adding new tests


To add a new test:

1. Put the recipe in the `test_recipes/` directory, with a filename like `test_XXX.slim`.
2. Add a config line in `recipe_specs.py` which determines which tests to run.
3. Add the additional necessary stuff to the SLiM recipe:

    * To check that haplotypes agree between SLiM and the tree sequence, just put

      ```
      source("init.slim");
      ```

      in `initialize()` (this defines functions); and then call

      ```
      outputMutationResult();
      ```

      at the end (well, whenever you want, really;).

    * To mark individuals in the initial generation with particular mutation types so we
      can check *something* even with mutation recording turned off, instead of sourcing
      "init.slim", source "init_marked_mutations.slim", i.e. place

      ```
      source("init_marked_mutations.slim");
	  ```

      in `initialize()` and then call

      ```
      initializeMarks(n_marks);
      ```

      in `1` (after adding the subpop), and then, as before,

      ```
      outputMutationResult();
      ```
      at the end. Note that this will only work properly if there is no new mutation.

    * Make sure you call `initializeTreeSeq();` after the `source(...)` line

    * Add `chooseAncestralSamples(5)` to some generations along the way
      to add some individuals as "ancestral samples". This is done by "remembering" them.

    * If you want to check that individuals in the simulation
      are present in generated tree sequence, call

      ```
      outputIndividuals()
      ```

      By default this simply saves a list of the individuals from the simulation at the point
      when it is called, and compares it to the list in the tree sequence. However, you can
      add individuals to the saved list by calling `addIndividuals(individuals)` at any point
      in the simulation.


To temporarily turn off a test just comment out the appropriate line in `recipe_specs.py`
