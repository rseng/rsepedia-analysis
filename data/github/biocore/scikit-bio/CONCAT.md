# scikit-bio changelog

## Version 0.5.6-dev

### Features

* Added support for `np.float32` with `DissimilarityMatrix` objects.
* Added support for method and number_of_dimensions to permdisp reducing the runtime by 100x at 4000 samples, [issue #1769](https://github.com/biocore/scikit-bio/pull/1769).
* OrdinationResults object is now accepted as input for permdisp.

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

* Avoid an implicit data copy on construction of `DissimilarityMatrix` objects.
* Avoid validation on copy of `DissimilarityMatrix` and `DistanceMatrix` objects, see [PR #1747](https://github.com/biocore/scikit-bio/pull/1747)
* Use an optimized version of symmetry check in DistanceMatrix, see [PR #1747](https://github.com/biocore/scikit-bio/pull/1747)
* Avoid performing filtering when ids are identical, see [PR #1752](https://github.com/biocore/scikit-bio/pull/1752)
* center_distance_matrix has been re-implemented in cython for both speed and memory use. Indirectly speeds up pcoa [PR #1749](https://github.com/biocore/scikit-bio/pull/1749)
* Use a memory-optimized version of permute in DistanceMatrix, see [PR #1756](https://github.com/biocore/scikit-bio/pull/1756).
* Refactor pearson and spearman skbio.stats.distance.mantel implementations to drastically improve memory locality. Also cache intermediate results that are invariant across permutations, see [PR #1756](https://github.com/biocore/scikit-bio/pull/1756).
* Refactor permanova to remove intermediate buffers and cythonize the internals, see [PR #1768](https://github.com/biocore/scikit-bio/pull/1768).

### Bug fixes

* Fix windows and 32bit incompatibility in `unweighted_unifrac`.

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous

* Specify build dependencies in pyproject.toml. This allows the package to be installed without having to first manually install numpy.

* Update hdmedians package to a version which doesn't require an initial manual numpy install.

* Now buildable on non-x86 platforms due to use of the [SIMD Everywhere](https://github.com/simd-everywhere/simde) library.

## Version 0.5.6

### Features

* Added option to return a capture group compiled regex pattern to any class inheriting ``GrammaredSequence`` through the ``to_regex`` method. ([#1431](https://github.com/biocore/scikit-bio/issues/1431))

* Added `Dissimilarity.within` and `.between` to obtain the respective distances and express them as a `DataFrame`. ([#1662](https://github.com/biocore/scikit-bio/pull/1662))

* Added Kendall Tau as possible correlation method in the `skbio.stats.distance.mantel` function ([#1675](https://github.com/biocore/scikit-bio/issues/1675)).

* Added support for IUPAC amino acid codes U (selenocysteine), O (pyrrolysine), and J (leucine or isoleucine). ([#1576](https://github.com/biocore/scikit-bio/issues/1576)

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

* Changed `skbio.tree.TreeNode.support` from a method to a property.
* Added `assign_supports` method to `skbio.tree.TreeNode` to extract branch support values from node labels.
* Modified the way a node's label is printed: `support:name` if both exist, or `support` or `name` if either exists.

### Performance enhancements

### Bug fixes
* Require `Sphinx <= 3.0`. Newer Sphinx versions caused build errors. [#1719](https://github.com/biocore/scikit-bio/pull/1719)

* * `skbio.stats.ordination` tests have been relaxed. ([#1713](https://github.com/biocore/scikit-bio/issues/1713))

* Fixes build errors for newer versions of NumPy, Pandas, and SciPy.

* Corrected a criticial bug in `skbio.alignment.StripedSmithWaterman`/`skbio.alignment.local_pairwise_align_ssw` which would cause the formatting of the aligned sequences to misplace gap characters by the number of gap characters present in the opposing aligned sequence up to that point. This was caused by a faulty implementation of CIGAR string parsing, see [#1679](https://github.com/biocore/scikit-bio/pull/1679) for full details.

* Fixes build errors for newer versions of NumPy, Pandas, and SciPy.

* Corrected a criticial bug in `skbio.alignment.StripedSmithWaterman`/`skbio.alignment.local_pairwise_align_ssw` which would cause the formatting of the aligned sequences to misplace gap characters by the number of gap characters present in the opposing aligned sequence up to that point. This was caused by a faulty implementation of CIGAR string parsing, see [#1679](https://github.com/biocore/scikit-bio/pull/1679) for full details.

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous

* `skbio.diversity.beta_diversity` now accepts a pandas DataFrame as input.

* Avoid pandas 1.0.0 import warning ([#1688](https://github.com/biocore/scikit-bio/issues/1688))

* Added support for Python 3.8 and dropped support for Python 3.5.

* This version now depends on `scipy >= 1.3` and `pandas >= 1.0`.

## Version 0.5.5 (2018-12-10)

### Features

* `skbio.stats.composition` now has methods to compute additive log-ratio transformation and inverse additive log-ratio transformation (`alr`, `alr_inv`) as well as a method to build a basis from a sequential binary partition (`sbp_basis`).

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

### Bug fixes

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous
* Python 3.6 and 3.7 compatibility is now supported

* A pytest runner is shipped with every installation ([#1633](https://github.com/biocore/scikit-bio/pull/1633))

* The nosetest framework has been replaced in favor of pytest ([#1624](https://github.com/biocore/scikit-bio/pull/1624))

* The numpy docs are deprecated in favor of [Napoleon](http://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html) ([#1629](https://github.com/biocore/scikit-bio/pull/1629))

* This version is now compatible with NumPy >= 1.9.2 and Pandas >= 0.23. ([#1627](https://github.com/biocore/scikit-bio/pull/1627))

## Version 0.5.4 (2018-08-23)

### Features

* Added `FSVD`, an alternative fast heuristic method to perform Principal Coordinates Analysis, to `skbio.stats.ordination.pcoa`.

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

* Added optimized utility methods `f_matrix_inplace` and `e_matrix_inplace` which perform `f_matrix` and `e_matrix` computations in-place and are used by the new `center_distance_matrix` method in `skbio.stats.ordination`.

### Bug fixes

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous

## Version 0.5.3 (2018-08-07)

### Features
* Added `unpack` and `unpack_by_func` methods to `skbio.tree.TreeNode` to unpack one or multiple internal nodes. The `unpack` operation removes an internal node and regrafts its children to its parent while retaining the overall length. ([#1572](https://github.com/biocore/scikit-bio/pull/1572))
* Added `support` to `skbio.tree.TreeNode` to return the support value of a node.
* Added `permdisp` to `skbio.stats.distance` to test for the homogeniety of groups. ([#1228](https://github.com/biocore/scikit-bio/issues/1228)).

* Added `pcoa_biplot` to `skbio.stats.ordination` to project descriptors into a PCoA plot.

* Fixed pandas to 0.22.0 due to this: https://github.com/pandas-dev/pandas/issues/20527

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

### Bug fixes

* Relaxing type checking in diversity calculations.  ([#1583](https://github.com/biocore/scikit-bio/issues/1583)).

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous


## Version 0.5.2 (2018-04-18)

### Features
* Added ``skbio.io.format.embl`` for reading and writing EMBL files for ``DNA``, ``RNA`` and ``Sequence`` classes.

* Removing ValueError check in `skbio.stats._subsample.subsample_counts` when `replace=True` and `n` is greater than the number of items in counts.  [#1527](https://github.com/biocore/scikit-bio/pull/1527)

* Added ``skbio.io.format.gff3`` for reading and writing GFF3 files for ``DNA``, ``Sequence``, and ``IntervalMetadata`` classes. ([#1450](https://github.com/biocore/scikit-bio/pull/1450))

* `skbio.metadata.IntervalMetadata` constructor has a new keyword argument, `copy_from`, for creating an `IntervalMetadata` object from an existing `IntervalMetadata` object with specified `upper_bound`.

* `skbio.metadata.IntervalMetadata` constructor allows `None` as a valid value for `upper_bound`. An `upper_bound` of `None` means that the `IntervalMetadata` object has no upper bound.

* `skbio.metadata.IntervalMetadata.drop` has a new boolean parameter `negate` to indicate whether to drop or keep the specified `Interval` objects.

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements
* `skbio.tree.nj` wall-clock runtime was decreased by 99% for a 500x500 distance matrix and 93% for a 100x100 distance matrix. ([#1512](https://github.com/biocore/scikit-bio/pull/1512), [#1513](https://github.com/biocore/scikit-bio/pull/1513))

### Bug fixes
* The `include_self` parameter was not being honored in `skbio.TreeNode.tips`. The scope of this bug was that if `TreeNode.tips` was called on a tip, it would always result in an empty `list` when unrolled.

* In `skbio.stats.ordination.ca`, `proportion_explained` was missing in the returned `OrdinationResults` object. ([#1345](https://github.com/biocore/scikit-bio/issues/1345))

* `skbio.diversity.beta_diversity` now handles qualitative metrics as expected such that `beta_diversity('jaccard', mat) == beta_diversity('jaccard', mat > 0)`. Please see [#1549](https://github.com/biocore/scikit-bio/issues/1549) for further detail.

* `skbio.stats.ordination.rda` The occasional column mismatch in output `biplot_scores` is fixed ([#1519](https://github.com/biocore/scikit-bio/issues/1519)).

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous
* scikit-bio now depends on pandas >= 0.19.2, and is compatible with newer pandas versions (e.g. 0.20.3) that were previously incompatible.
* scikit-bio now depends on `numpy >= 1.9.2, < 1.14.0` for compatibility with Python 3.4, 3.5, and 3.6 and the available numpy conda packages in `defaults` and `conda-forge` channels.
* added support for running tests from `setup.py`. Both `python setup.py nosetests` and `python setup.py test` are now supported, however `python setup.py test` will only run a subset of the full test suite. ([#1341](https://github.com/biocore/scikit-bio/issues/1341))

## Version 0.5.1 (2016-11-12)

### Features
* Added `IntervalMetadata` and `Interval` classes in `skbio.metadata` to store, query, and manipulate information of a sub-region of a sequence. ([#1414](https://github.com/biocore/scikit-bio/issues/1414))
* `Sequence` and its child classes (including `GrammaredSequence`, `RNA`, `DNA`, `Protein`) now accept `IntervalMetadata` in their constructor API. Some of their relevant methods are also updated accordingly. ([#1430](https://github.com/biocore/scikit-bio/pull/1430))
* GenBank parser now reads and writes `Sequence` or its subclass objects with `IntervalMetadata`. ([#1440](https://github.com/biocore/scikit-bio/pull/1440))
* `DissimilarityMatrix` now has a new constructor method called `from_iterable`. ([#1343](https://github.com/biocore/scikit-bio/issues/1343)).
* `DissimilarityMatrix` now allows non-hollow matrices. ([#1343](https://github.com/biocore/scikit-bio/issues/1343)).
* `DistanceMatrix.from_iterable` now accepts a `validate=True` parameter. ([#1343](https://github.com/biocore/scikit-bio/issues/1343)).
* ``DistanceMatrix`` now has a new method called ``to_series`` to create a ``pandas.Series`` from a ``DistanceMatrix`` ([#1397](https://github.com/biocore/scikit-bio/issues/1397)).
* Added parallel beta diversity calculation support via `skbio.diversity.block_beta_diversity`. The issue and idea is discussed in ([#1181](https://github.com/biocore/scikit-bio/issues/1181), while the actual code changes are in [#1352](https://github.com/biocore/scikit-bio/pull/1352)).


### Backward-incompatible changes [stable]
* The constructor API for `Sequence` and its child classes (including `GrammaredSequence`, `RNA`, `DNA`, `Protein`) are changed from `(sequence, metadata=None, positional_metadata=None, lowercase=False)` to `(sequence, metadata=None, positional_metadata=None, interval_metadata=None, lowercase=False)`

  The changes are made to allow these classes to adopt `IntervalMetadata` object for interval features on the sequence. The `interval_metadata` parameter is added imediately after `positional_metadata` instead of appended to the end, because it is more natural and logical and, more importantly, because it is unlikely in practice to break user code. A user's code would break only if they had supplied `metadata`, `postional_metadata`, and `lowercase` parameters positionally. In the unlikely event that this happens, users will get an error telling them a bool isn't a valid `IntervalMetadata` type, so it won't silently produce buggy behavior.

### Backward-incompatible changes [experimental]
* Modifying basis handling in `skbio.stats.composition.ilr_inv` prior to checking for orthogonality.  Now the basis is strictly assumed to be in the Aitchison simplex.
* `DistanceMatrix.from_iterable` default behavior is now to validate matrix by computing all pairwise distances. Pass `validate=False` to get the previous behavior (no validation, but faster execution).([#1343](https://github.com/biocore/scikit-bio/issues/1343)).
* GenBank I/O now parses sequence features into the attribute of `interval_metadata` instead of `positiona_metadata`. And the key of `FEATURES` is removed from `metadata` attribute.

### Performance enhancements
* `TreeNode.shear` was rewritten for approximately a 25% performance increase. ([#1399](https://github.com/biocore/scikit-bio/pull/1399))
* The `IntervalMetadata` allows dramatic decrease in memory usage in reading GenBank files of feature rich sequences. ([#1159](https://github.com/biocore/scikit-bio/issues/1159))

### Bug fixes

* `skbio.tree.TreeNode.prune` and implicitly `skbio.tree.TreeNode.shear` were not handling a situation in which a parent was validly removed during pruning operations as may happen if the resulting subtree does not include the root. Previously, an `AttributeError` would raise as `parent` would be `None` in this situation.
* numpy linking was fixed for installation under El Capitan.
* A bug was introduced in #1398 into `TreeNode.prune` and fixed in #1416 in which, under the special case of a single descendent existing from the root, the resulting children parent references were not updated. The cause of the bug was a call made to `self.children.extend` as opposed to `self.extend` where the former is a `list.extend` without knowledge of the tree, while the latter is `TreeNode.extend` which is able to adjust references to `self.parent`.

### Miscellaneous

* Removed deprecated functions from `skbio.util`: `is_casava_v180_or_later`, `remove_files`, and `create_dir`.
* Removed deprecated `skbio.Sequence.copy` method.

## Version 0.5.0 (2016-06-14)

**IMPORTANT**: scikit-bio is no longer compatible with Python 2. scikit-bio is compatible with Python 3.4 and later.

### Features
* Added more descriptive error message to `skbio.io.registry` when attempting to read without specifying `into` and when there is no generator reader. ([#1326](https://github.com/biocore/scikit-bio/issues/1326))
* Added support for reference tags to `skbio.io.format.stockholm` reader and writer. ([#1348](https://github.com/biocore/scikit-bio/issues/1348))
* Expanded error message in `skbio.io.format.stockholm` reader when `constructor` is not passed, in order to provide better explanation to user. ([#1327](https://github.com/biocore/scikit-bio/issues/1327))
* Added `skbio.sequence.distance.kmer_distance` for computing the kmer distance between two sequences. ([#913](https://github.com/biocore/scikit-bio/issues/913))
* Added `skbio.sequence.Sequence.replace` for assigning a character to positions in a `Sequence`. ([#1222](https://github.com/biocore/scikit-bio/issues/1222))
* Added support for `pandas.RangeIndex`, lowering the memory footprint of default integer index objects. `Sequence.positional_metadata` and `TabularMSA.positional_metadata` now use `pd.RangeIndex` as the positional metadata index. `TabularMSA` now uses `pd.RangeIndex` as the default index. Usage of `pd.RangeIndex` over the previous `pd.Int64Index` [should be transparent](http://pandas.pydata.org/pandas-docs/version/0.18.0/whatsnew.html#range-index), so these changes should be non-breaking to users. scikit-bio now depends on pandas >= 0.18.0 ([#1308](https://github.com/biocore/scikit-bio/issues/1308))
* Added `reset_index=False` parameter to `TabularMSA.append` and `TabularMSA.extend` for resetting the MSA's index to the default index after appending/extending.
* Added support for partial pairwise calculations via `skbio.diversity.partial_beta_diversity`. ([#1221](https://github.com/biocore/scikit-bio/issues/1221), [#1337](https://github.com/biocore/scikit-bio/pull/1337)). This function is immediately deprecated as its return type will change in the future and should be used with caution in its present form (see the function's documentation for details).
* `TemporaryFile` and `NamedTemporaryFile` are now supported IO sources for `skbio.io` and related functionality.  ([#1291](https://github.com/biocore/scikit-bio/issues/1291))
* Added `tree_node_class=TreeNode` parameter to `skbio.tree.majority_rule` to support returning consensus trees of type `TreeNode` (the default) or a type that has the same interface as `TreeNode` (e.g. `TreeNode` subclasses) ([#1193](https://github.com/biocore/scikit-bio/pull/1193))
* `TreeNode.from_linkage_matrix` and `TreeNode.from_taxonomy` now support constructing `TreeNode` subclasses. `TreeNode.bifurcate` now supports `TreeNode` subclasses ([#1193](https://github.com/biocore/scikit-bio/pull/1193))
* The `ignore_metadata` keyword has been added to `TabularMSA.iter_positions` to improve performance when metadata is not necessary.
* Pairwise aligners in `skbio.alignment` now propagate per-sequence `metadata` objects (this does not include `positional_metadata`).

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]
* `TabularMSA.append` and `TabularMSA.extend` now require one of `minter`, `index`, or `reset_index` to be provided when incorporating new sequences into an MSA. Previous behavior was to auto-increment the index labels if `minter` and `index` weren't provided and the MSA had a default integer index, otherwise error. Use `reset_index=True` to obtain the previous behavior in a more explicit way.
* `skbio.stats.composition.ancom` now returns two `pd.DataFrame` objects, where it previously returned one. The first contains the ANCOM test results, as before, and the second contains percentile abundances of each feature in each group. The specific percentiles that are computed and returned is controlled by the new `percentiles` parameter to `skbio.stats.composition.ancom`. In the future, this second `pd.DataFrame` will not be returned by this function, but will be available through the [contingency table API](https://github.com/biocore/scikit-bio/issues/848). ([#1293](https://github.com/biocore/scikit-bio/issues/1293))
* `skbio.stats.composition.ancom` now performs multiple comparisons correction by default. The previous behavior of not performing multiple comparisons correction can be achieved by passing ``multiple_comparisons_correction=None``.
* The ``reject`` column in the first ``pd.DataFrame`` returned from `skbio.stats.composition.ancom` has been renamed ``Reject null hypothesis`` for clarity. ([#1375](https://github.com/biocore/scikit-bio/issues/1375))

### Bug fixes
* Fixed row and column names to `biplot_scores` in the `OrdinationResults` object from `skbio.stats.ordination`. This fix affect the `cca` and `rda` methods. ([#1322](https://github.com/biocore/scikit-bio/issues/1322))
* Fixed bug when using `skbio.io.format.stockholm` reader on file with multi-line tree with no id. Previously this raised an `AttributeError`, now it correctly handles this type of tree. ([#1334](https://github.com/biocore/scikit-bio/issues/1334))
* Fixed bug when reading Stockholm files with GF or GS features split over multiple lines. Previously, the feature text was simply concatenated because it was assumed to have trailing whitespace. There are examples of Stockholm files with and without trailing whitespace for multi-line features, so the `skbio.io.format.stockholm` reader now adds a single space when concatenating feature text without trailing whitespace to avoid joining words together. Multi-line trees stored as GF metadata are concatenated as they appear in the file; a space is not added when concatenating. ([#1328](https://github.com/biocore/scikit-bio/issues/1328))
* Fixed bug when using `Sequence.iter_kmers` on empty `Sequence` object. Previously this raised a `ValueError`, now it returns
an empty generator.
* Fixed minor bug where adding sequences to an empty `TabularMSA` with MSA-wide `positional_metadata` would result in a `TabularMSA` object in an inconsistent state. This could happen using `TabularMSA.append` or `TabularMSA.extend`. This bug only affects a `TabularMSA` object *without* sequences that has MSA-wide `positional_metadata` (for example, `TabularMSA([], positional_metadata={'column': []})`).
* `TreeNode.distance` now handles the situation in which `self` or `other` are ancestors. Previosly, a node further up the tree was used resulting in inflated distances. ([#807](https://github.com/biocore/scikit-bio/issues/807))
* `TreeNode.prune` can now handle a root with a single descendent. Previously, the root was ignored from possibly having a single descendent. ([#1247](https://github.com/biocore/scikit-bio/issues/1247))
* Providing the `format` keyword to `skbio.io.read` when creating a generator with an empty file will now return an empty generator instead of raising `StopIteration`. ([#1313](https://github.com/biocore/scikit-bio/issues/1313))
* `OrdinationResults` is now importable from `skbio` and `skbio.stats.ordination` and correctly linked from the documentation ([#1205](https://github.com/biocore/scikit-bio/issues/1205))
* Fixed performance bug in pairwise aligners resulting in 100x worse performance than in 0.2.4.

### Deprecated functionality [stable]
* Deprecated use of the term "non-degenerate", in favor of "definite". `GrammaredSequence.nondegenerate_chars`, `GrammaredSequence.nondegenerates`, and `GrammaredSequence.has_nondegenerates` have been renamed to `GrammaredSequence.definite_chars`, `GrammaredSequence.definites`, and `GrammaredSequence.has_definites`, respectively. The old names will be removed in scikit-bio 0.5.2. Relevant affected public classes include `GrammaredSequence`, `DNA`, `RNA`, and `Protein`.

### Deprecated functionality [experimental]
* Deprecated function `skbio.util.create_dir`. This function will be removed in scikit-bio 0.5.1. Please use the Python standard library
functionality described [here](https://docs.python.org/2/library/os.html#os.makedirs). ([#833](https://github.com/biocore/scikit-bio/issues/833))
* Deprecated function `skbio.util.remove_files`. This function will be removed in scikit-bio 0.5.1. Please use the Python standard library
functionality described [here](https://docs.python.org/2/library/os.html#os.remove). ([#833](https://github.com/biocore/scikit-bio/issues/833))
* Deprecated function `skbio.util.is_casava_v180_or_later`. This function will be removed in 0.5.1. Functionality moved to FASTQ sniffer.
([#833](https://github.com/biocore/scikit-bio/issues/833))

### Miscellaneous
* When installing scikit-bio via `pip`, numpy must now be installed first ([#1296](https://github.com/biocore/scikit-bio/issues/1296))

## Version 0.4.2 (2016-02-17)

Minor maintenance release. **This is the last Python 2.7 compatible release. Future scikit-bio releases will only support Python 3.**

### Features
* Added `skbio.tree.TreeNode.bifurcate` for converting multifurcating trees into bifurcating trees. ([#896](https://github.com/biocore/scikit-bio/issues/896))
* Added `skbio.io.format.stockholm` for reading Stockholm files into a `TabularMSA` and writing from a `TabularMSA`. ([#967](https://github.com/biocore/scikit-bio/issues/967))
* scikit-bio `Sequence` objects have better compatibility with numpy. For example, calling `np.asarray(sequence)` now converts the sequence to a numpy array of characters (the same as calling `sequence.values`).
* Added `skbio.sequence.distance` subpackage for computing distances between scikit-bio `Sequence` objects ([#913](https://github.com/biocore/scikit-bio/issues/913))
* Added ``skbio.sequence.GrammaredSequence``, which can be inherited from to create grammared sequences with custom alphabets (e.g., for use with TabularMSA) ([#1175](https://github.com/biocore/scikit-bio/issues/1175))
* Added ``skbio.util.classproperty`` decorator

### Backward-incompatible changes [stable]
* When sniffing or reading a file (`skbio.io.sniff`, `skbio.io.read`, or the object-oriented `.read()` interface), passing `newline` as a keyword argument to `skbio.io.open` now raises a `TypeError`. This backward-incompatible change to a stable API is necessary because it fixes a bug (more details in bug fix section below).
* When reading a FASTQ or QSEQ file and passing `variant='solexa'`, `ValueError` is now raised instead of `NotImplementedError`. This backward-incompatible change to a stable API is necessary to avoid creating a spin-locked process due to [a bug in Python](https://bugs.python.org/issue25786). See [#1256](https://github.com/biocore/scikit-bio/issues/1256) for details. This change is temporary and will be reverted to `NotImplementedError` when the bug is fixed in Python.

### Backward-incompatible changes [experimental]
* `skbio.io.format.genbank`: When reading GenBank files, the date field of the LOCUS line is no longer parsed into a `datetime.datetime` object and is left as a string. When writing GenBank files, the locus date metadata is expected to be a string instead of a `datetime.datetime` object ([#1153](https://github.com/biocore/scikit-bio/issues/1153))
* `Sequence.distance` now converts the input sequence (`other`) to its type before passing both sequences to `metric`. Previous behavior was to always convert to `Sequence`.

### Bug fixes
* Fixed bug when using `Sequence.distance` or `DistanceMatrix.from_iterable` to compute distances between `Sequence` objects with differing `metadata`/`positional_metadata` and passing `metric=scipy.spatial.distance.hamming` ([#1254](https://github.com/biocore/scikit-bio/issues/1254))
* Fixed performance bug when computing Hamming distances between `Sequence` objects in `DistanceMatrix.from_iterable` ([#1250](https://github.com/biocore/scikit-bio/issues/1250))
* Changed `skbio.stats.composition.multiplicative_replacement` to raise an error whenever a large value of `delta` is chosen ([#1241](https://github.com/biocore/scikit-bio/issues/1241))
* When sniffing or reading a file (`skbio.io.sniff`, `skbio.io.read`, or the object-oriented `.read()` interface), passing `newline` as a keyword argument to `skbio.io.open` now raises a `TypeError`. The file format's `newline` character will be used when opening the file. Previous behavior allowed overriding the format's `newline` character but this could cause issues with readers that assume newline characters are those defined by the file format (which is an entirely reasonable assumption). This bug is very unlikely to have surfaced in practice as the default `newline` behavior is *universal newlines mode*.
* DNA, RNA, and Protein are no longer inheritable because they assume an IUPAC alphabet.
* `DistanceMatrix` constructor provides more informative error message when data contains NaNs ([#1276](https://github.com/biocore/scikit-bio/issues/1276))

### Miscellaneous
* Warnings raised by scikit-bio now share a common subclass ``skbio.util.SkbioWarning``.

## Version 0.4.1 (2015-12-09)

### Features
* The ``TabularMSA`` object was added to represent and operate on tabular multiple sequence alignments. This satisfies [RFC 1](https://github.com/biocore/scikit-bio-rfcs/blob/master/active/001-tabular-msa.md). See the ``TabularMSA`` docs for full details.
* Added phylogenetic diversity metrics, including weighted UniFrac, unweighted UniFrac, and Faith's Phylogenetic Diversity. These are accessible as ``skbio.diversity.beta.unweighted_unifrac``, ``skbio.diversity.beta.weighted_unifrac``, and ``skbio.diversity.alpha.faith_pd``, respectively.
* Addition of the function ``skbio.diversity.alpha_diversity`` to support applying an alpha diversity metric to multiple samples in one call.
* Addition of the functions ``skbio.diversity.get_alpha_diversity_metrics`` and ``skbio.diversity.get_beta_diversity_metrics`` to support discovery of the alpha and beta diversity metrics implemented in scikit-bio.
* Added `skbio.stats.composition.ancom` function, a test for OTU differential abundance across sample categories. ([#1054](https://github.com/biocore/scikit-bio/issues/1054))
* Added `skbio.io.format.blast7` for reading BLAST+ output format 7 or BLAST output format 9 files into a `pd.DataFrame`. ([#1110](https://github.com/biocore/scikit-bio/issues/1110))
* Added `skbio.DissimilarityMatrix.to_data_frame` method for creating a ``pandas.DataFrame`` from a `DissimilarityMatrix` or `DistanceMatrix`. ([#757](https://github.com/biocore/scikit-bio/issues/757))
* Added support for one-dimensional vector of dissimilarities in `skbio.stats.distance.DissimilarityMatrix`
constructor. ([#6240](https://github.com/biocore/scikit-bio/issues/624))
* Added `skbio.io.format.blast6` for reading BLAST+ output format 6 or BLAST output format 8 files into a `pd.DataFrame`. ([#1110](https://github.com/biocore/scikit-bio/issues/1110))
* Added `inner`, `ilr`, `ilr_inv` and `clr_inv`, ``skbio.stats.composition``, which enables linear transformations on compositions ([#892](https://github.com/biocore/scikit-bio/issues/892)
* Added ``skbio.diversity.alpha.pielou_e`` function as an evenness metric of alpha diversity. ([#1068](https://github.com/biocore/scikit-bio/issues/1068))
* Added `to_regex` method to `skbio.sequence._iupac_sequence` ABC - it returns a regex object that matches all non-degenerate versions of the sequence.
* Added ``skbio.util.assert_ordination_results_equal`` function for comparing ``OrdinationResults`` objects in unit tests.
* Added ``skbio.io.format.genbank`` for reading and writing GenBank/GenPept for ``DNA``, ``RNA``, ``Protein`` and ``Sequence`` classes.
* Added ``skbio.util.RepresentationWarning`` for warning about substitutions, assumptions, or particular alterations that were made for the successful completion of a process.
* ``TreeNode.tip_tip_distances`` now supports nodes without an associated length. In this case, a length of 0.0 is assumed and an ``skbio.util.RepresentationWarning`` is raised. Previous behavior was to raise a ``NoLengthError``. ([#791](https://github.com/biocore/scikit-bio/issues/791))
* ``DistanceMatrix`` now has a new constructor method called `from_iterable`.
* ``Sequence`` now accepts ``lowercase`` keyword like ``DNA`` and others. Updated ``fasta``, ``fastq``, and ``qseq`` readers/writers for ``Sequence`` to reflect this.
* The ``lowercase`` method has been moved up to ``Sequence`` meaning all sequence objects now have a ``lowercase`` method.
* Added ``reverse_transcribe`` class method to ``RNA``.
* Added `Sequence.observed_chars` property for obtaining the set of observed characters in a sequence. ([#1075](https://github.com/biocore/scikit-bio/issues/1075))
* Added `Sequence.frequencies` method for computing character frequencies in a sequence. ([#1074](https://github.com/biocore/scikit-bio/issues/1074))
* Added experimental class-method ``Sequence.concat`` which will produce a new sequence from an iterable of existing sequences. Parameters control how positional metadata is propagated during a concatenation.
* ``TreeNode.to_array`` now supports replacing ``nan`` branch lengths in the resulting branch length vector with the value provided as ``nan_length_value``.
* ``skbio.io.format.phylip`` now supports sniffing and reading strict, sequential PHYLIP-formatted files into ``skbio.Alignment`` objects. ([#1006](https://github.com/biocore/scikit-bio/issues/1006))
* Added `default_gap_char` class property to ``DNA``, ``RNA``, and ``Protein`` for representing gap characters in a new sequence.

### Backward-incompatible changes [stable]
* `Sequence.kmer_frequencies` now returns a `dict`. Previous behavior was to return a `collections.Counter` if `relative=False` was passed, and a `collections.defaultdict` if `relative=True` was passed. In the case of a missing key, the `Counter` would return 0 and the `defaultdict` would return 0.0. Because the return type is now always a `dict`, attempting to access a missing key will raise a `KeyError`. This change *may* break backwards-compatibility depending on how the `Counter`/`defaultdict` is being used. We hope that in most cases this change will not break backwards-compatibility because both `Counter` and `defaultdict` are `dict` subclasses.

   If the previous behavior is desired, convert the `dict` into a `Counter`/`defaultdict`:

    ```python
    import collections
    from skbio import Sequence
    seq = Sequence('ACCGAGTTTAACCGAATA')

    # Counter
    freqs_dict = seq.kmer_frequencies(k=8)
    freqs_counter = collections.Counter(freqs_dict)

    # defaultdict
    freqs_dict = seq.kmer_frequencies(k=8, relative=True)
    freqs_default_dict = collections.defaultdict(float, freqs_dict)
    ```

   **Rationale:** We believe it is safer to return `dict` instead of `Counter`/`defaultdict` as this may prevent error-prone usage of the return value. Previous behavior allowed accessing missing kmers, returning 0 or 0.0 depending on the `relative` parameter. This is convenient in many cases but also potentially misleading. For example, consider the following code:

    ```python
    from skbio import Sequence
    seq = Sequence('ACCGAGTTTAACCGAATA')
    freqs = seq.kmer_frequencies(k=8)
    freqs['ACCGA']
    ```

    Previous behavior would return 0 because the kmer `'ACCGA'` is not present in the `Counter`. In one respect this is the correct answer because we asked for kmers of length 8; `'ACCGA'` is a different length so it is not included in the results. However, we believe it is safer to avoid this implicit behavior in case the user assumes there are no `'ACCGA'` kmers in the sequence (which there are!). A `KeyError` in this case is more explicit and forces the user to consider their query. Returning a `dict` will also be consistent with `Sequence.frequencies`.

### Backward-incompatible changes [experimental]
* Replaced ``PCoA``, ``CCA``, ``CA`` and ``RDA`` in ``skbio.stats.ordination`` with equivalent functions ``pcoa``, ``cca``, ``ca`` and ``rda``. These functions now take ``pd.DataFrame`` objects.
* Change ``OrdinationResults`` to have its attributes based on ``pd.DataFrame`` and ``pd.Series`` objects, instead of pairs of identifiers and values. The changes are as follows:
    - ``species`` and ``species_ids`` have been replaced by a ``pd.DataFrame`` named ``features``.
    - ``site`` and ``site_ids`` have been replaced by a ``pd.DataFrame`` named ``samples``.
    - ``eigvals`` is now a ``pd.Series`` object.
    - ``proportion_explained`` is now a ``pd.Series`` object.
    - ``biplot`` is now a ``pd.DataFrame`` object named ``biplot_scores``.
    - ``site_constraints`` is now a ``pd.DataFrame`` object named ``sample_constraints``.
* ``short_method_name`` and ``long_method_name`` are now required arguments of the ``OrdinationResults`` object.
* Removed `skbio.diversity.alpha.equitability`. Please use `skbio.diversity.alpha.pielou_e`, which is more accurately named and better documented. Note that `equitability` by default used logarithm base 2 while `pielou_e` uses logarithm base `e` as described in Heip 1974.
* ``skbio.diversity.beta.pw_distances`` is now called ``skbio.diversity.beta_diversity``. This function no longer defines a default metric, and ``metric`` is now the first argument to this function. This function can also now take a pairwise distances function as ``pairwise_func``.
* Deprecated function ``skbio.diversity.beta.pw_distances_from_table`` has been removed from scikit-bio as scheduled. Code that used this should be adapted to use ``skbio.diversity.beta_diversity``.
* ``TreeNode.index_tree`` now returns a 2-D numpy array as its second return value (the child node index) instead of a 1-D numpy array.
* Deprecated functions `skbio.draw.boxplots` and `skbio.draw.grouped_distributions` have been removed from scikit-bio as scheduled. These functions generated plots that were not specific to bioinformatics. These types of plots can be generated with seaborn or another general-purpose plotting package.
* Deprecated function `skbio.stats.power.bootstrap_power_curve` has been removed from scikit-bio as scheduled. Use `skbio.stats.power.subsample_power` or `skbio.stats.power.subsample_paired_power` followed by `skbio.stats.power.confidence_bound`.
* Deprecated function `skbio.stats.spatial.procrustes` has been removed from scikit-bio as scheduled in favor of `scipy.spatial.procrustes`.
* Deprecated class `skbio.tree.CompressedTrie` and function `skbio.tree.fasta_to_pairlist` have been removed from scikit-bio as scheduled in favor of existing general-purpose Python trie packages.
* Deprecated function `skbio.util.flatten` has been removed from scikit-bio as scheduled in favor of solutions available in the Python standard library (see [here](http://stackoverflow.com/a/952952/3639023) and [here](http://stackoverflow.com/a/406199/3639023) for examples).
* Pairwise alignment functions in `skbio.alignment` now return a tuple containing the `TabularMSA` alignment, alignment score, and start/end positions. The returned `TabularMSA`'s `index` is always the default integer index; sequence IDs are no longer propagated to the MSA. Additionally, the pairwise alignment functions now accept the following input types to align:
    - `local_pairwise_align_nucleotide`: `DNA` or `RNA`
    - `local_pairwise_align_protein`: `Protein`
    - `local_pairwise_align`: `IUPACSequence`
    - `global_pairwise_align_nucleotide`: `DNA`, `RNA`, or `TabularMSA[DNA|RNA]`
    - `global_pairwise_align_protein`: `Protein` or `TabularMSA[Protein]`
    - `global_pairwise_align`: `IUPACSequence` or `TabularMSA`
    - `local_pairwise_align_ssw`: `DNA`, `RNA`, or `Protein`. Additionally, this function now overrides the `protein` kwarg based on input type. `constructor` parameter was removed because the function now determines the return type based on input type.
* Removed `skbio.alignment.SequenceCollection` in favor of using a list or other standard library containers to store scikit-bio sequence objects (most `SequenceCollection` operations were simple list comprehensions). Use `DistanceMatrix.from_iterable` instead of `SequenceCollection.distances` (pass `key="id"` to exactly match original behavior).
* Removed `skbio.alignment.Alignment` in favor of `skbio.alignment.TabularMSA`.
* Removed `skbio.alignment.SequenceCollectionError` and `skbio.alignment.AlignmentError` exceptions as their corresponding classes no longer exist.

### Bug Fixes

* ``Sequence`` objects now handle slicing of empty positional metadata correctly. Any metadata that is empty will no longer be propagated by the internal ``_to`` constructor. ([#1133](https://github.com/biocore/scikit-bio/issues/1133))
* ``DissimilarityMatrix.plot()`` no longer leaves a white border around the
  heatmap it plots (PR #1070).
* TreeNode.root_at_midpoint`` no longer fails when a node with two equal length child branches exists in the tree. ([#1077](https://github.com/biocore/scikit-bio/issues/1077))
* ``TreeNode._set_max_distance``, as called through ``TreeNode.get_max_distance`` or ``TreeNode.root_at_midpoint`` would store distance information as ``list``s in the attribute ``MaxDistTips`` on each node in the tree, however, these distances were only valid for the node in which the call to ``_set_max_distance`` was made. The values contained in ``MaxDistTips`` are now correct across the tree following a call to ``get_max_distance``. The scope of impact of this bug is limited to users that were interacting directly with ``MaxDistTips`` on descendant nodes; this bug does not impact any known method within scikit-bio. ([#1223](https://github.com/biocore/scikit-bio/issues/1223))
* Added missing `nose` dependency to setup.py's `install_requires`. ([#1214](https://github.com/biocore/scikit-bio/issues/1214))
* Fixed issue that resulted in legends of ``OrdinationResult`` plots sometimes being truncated. ([#1210](https://github.com/biocore/scikit-bio/issues/1210))

### Deprecated functionality [stable]
* `skbio.Sequence.copy` has been deprecated in favor of `copy.copy(seq)` and `copy.deepcopy(seq)`.

### Miscellaneous
* Doctests are now written in Python 3.
* ``make test`` now validates MANIFEST.in using [check-manifest](https://github.com/mgedmin/check-manifest). ([#461](https://github.com/biocore/scikit-bio/issues/461))
* Many new alpha diversity equations added to ``skbio.diversity.alpha`` documentation. ([#321](https://github.com/biocore/scikit-bio/issues/321))
* Order of ``lowercase`` and ``validate`` keywords swapped in ``DNA``, ``RNA``, and ``Protein``.

## Version 0.4.0 (2015-07-08)

Initial beta release. In addition to the changes detailed below, the following
subpackages have been mostly or entirely rewritten and most of their APIs are
substantially different (and improved!):

* `skbio.sequence`
* `skbio.io`

The APIs of these subpackages are now stable, and all others are experimental. See the [API stability docs](https://github.com/biocore/scikit-bio/tree/0.4.0/doc/source/user/api_stability.rst) for more details, including what we mean by *stable* and *experimental* in this context. We recognize that this is a lot of backward-incompatible changes. To avoid these types of changes being a surprise to our users, our public APIs are now decorated to make it clear to developers when an API can be relied upon (stable) and when it may be subject to change (experimental).

### Features
* Added `skbio.stats.composition` for analyzing data made up of proportions
* Added new ``skbio.stats.evolve`` subpackage for evolutionary statistics. Currently contains a single function, ``hommola_cospeciation``, which implements a permutation-based test of correlation between two distance matrices.
* Added support for ``skbio.io.util.open_file`` and ``skbio.io.util.open_files`` to pull files from HTTP and HTTPS URLs. This behavior propagates to the I/O registry.
* FASTA/QUAL (``skbio.io.format.fasta``) and FASTQ (``skbio.io.format.fastq``) readers now allow blank or whitespace-only lines at the beginning of the file, between records, or at the end of the file. A blank or whitespace-only line in any other location will continue to raise an error [#781](https://github.com/biocore/scikit-bio/issues/781).
* scikit-bio now ignores leading and trailing whitespace characters on each line while reading FASTA/QUAL and FASTQ files.
* Added `ratio` parameter to `skbio.stats.power.subsample_power`. This allows the user to calculate power on groups for uneven size (For example, draw twice as many samples from Group B than Group A). If `ratio` is not set, group sizes will remain equal across all groups.
* Power calculations (`skbio.stats.power.subsample_power` and `skbio.stats.power.subsample_paired_power`) can use test functions that return multiple p values, like some multivariate linear regression models. Previously, the power calculations required the test to return a single p value.
* Added ``skbio.util.assert_data_frame_almost_equal`` function for comparing ``pd.DataFrame`` objects in unit tests.

### Performance enhancements
* The speed of quality score decoding has been significantly improved (~2x) when reading `fastq` files.
* The speed of `NucleotideSequence.reverse_complement` has been improved (~6x).

### Bug fixes
* Changed `Sequence.distance` to raise an error any time two sequences are passed of different lengths regardless of the `distance_fn` being passed. [(#514)](https://github.com/biocore/scikit-bio/issues/514)
* Fixed issue with ``TreeNode.extend`` where if given the children of another ``TreeNode`` object (``tree.children``), both trees would be left in an incorrect and unpredictable state. ([#889](https://github.com/biocore/scikit-bio/issues/889))
* Changed the way power was calculated in `subsample_paired_power` to move the subsample selection before the test is performed. This increases the number of Monte Carlo simulations performed during power estimation, and improves the accuracy of the returned estimate. Previous power estimates from `subsample_paired_power` should be disregarded and re-calculated. ([#910](https://github.com/biocore/scikit-bio/issues/910))
* Fixed issue where `randdm` was attempting to create asymmetric distance matrices.This was causing an error to be raised by the `DistanceMatrix` constructor inside of the `randdm` function, so that `randdm` would fail when attempting to create large distance matrices. ([#943](https://github.com/biocore/scikit-bio/issues/943))

### Deprecated functionality
* Deprecated `skbio.util.flatten`. This function will be removed in scikit-bio 0.3.1. Please use standard python library functionality
described here [Making a flat list out of lists of lists](http://stackoverflow.com/a/952952/3639023), [Flattening a shallow list](http://stackoverflow.com/a/406199/3639023) ([#833](https://github.com/biocore/scikit-bio/issues/833))
* Deprecated `skbio.stats.power.bootstrap_power_curve` will be removed in scikit-bio 0.4.1. It is deprecated in favor of using ``subsample_power`` or ``sample_paired_power`` to calculate a power matrix, and then the use of ``confidence_bounds`` to calculate the average and confidence intervals.

### Backward-incompatible changes
* Removed the following deprecated functionality:
    - `skbio.parse` subpackage, including `SequenceIterator`, `FastaIterator`, `FastqIterator`, `load`, `parse_fasta`, `parse_fastq`, `parse_qual`, `write_clustal`, `parse_clustal`, and `FastqParseError`; please use `skbio.io` instead.
    - `skbio.format` subpackage, including `fasta_from_sequence`, `fasta_from_alignment`, and `format_fastq_record`; please use `skbio.io` instead.
    - `skbio.alignment.SequenceCollection.int_map`; please use `SequenceCollection.update_ids` instead.
    - `skbio.alignment.SequenceCollection` methods `to_fasta` and `toFasta`; please use `SequenceCollection.write` instead.
    - `constructor` parameter in `skbio.alignment.Alignment.majority_consensus`; please convert returned biological sequence object manually as desired (e.g., `str(seq)`).
    - `skbio.alignment.Alignment.to_phylip`; please use `Alignment.write` instead.
    - `skbio.sequence.BiologicalSequence.to_fasta`; please use `BiologicalSequence.write` instead.
    - `skbio.tree.TreeNode` methods `from_newick`, `from_file`, and `to_newick`; please use `TreeNode.read` and `TreeNode.write` instead.
    - `skbio.stats.distance.DissimilarityMatrix` methods `from_file` and `to_file`; please use `DissimilarityMatrix.read` and `DissimilarityMatrix.write` instead.
    - `skbio.stats.ordination.OrdinationResults` methods `from_file` and `to_file`; please use `OrdinationResults.read` and `OrdinationResults.write` instead.
    - `skbio.stats.p_value_to_str`; there is no replacement.
    - `skbio.stats.subsample`; please use `skbio.stats.subsample_counts` instead.
    - `skbio.stats.distance.ANOSIM`; please use `skbio.stats.distance.anosim` instead.
    - `skbio.stats.distance.PERMANOVA`; please use `skbio.stats.distance.permanova` instead.
    - `skbio.stats.distance.CategoricalStatsResults`; there is no replacement, please use `skbio.stats.distance.anosim` or `skbio.stats.distance.permanova`, which will return a `pandas.Series` object.
* `skbio.alignment.Alignment.majority_consensus` now returns `BiologicalSequence('')` if the alignment is empty. Previously, `''` was returned.
* `min_observations` was removed from `skbio.stats.power.subsample_power` and `skbio.stats.power.subsample_paired_power`. The minimum number of samples for subsampling depends on the data set and statistical tests. Having a default parameter to set unnecessary limitations on the technique.

### Miscellaneous
* Changed testing procedures
    - Developers should now use `make test`
    - Users can use `python -m skbio.test`
    - Added `skbio.util._testing.TestRunner` (available through `skbio.util.TestRunner`). Used to provide a `test` method for each module init file. This class represents a unified testing path which wraps all `skbio` testing functionality.
    - Autodetect Python version and disable doctests for Python 3.
* `numpy` is no longer required to be installed before installing scikit-bio!
* Upgraded checklist.py to check source files non-conforming to [new header style](http://scikit-bio.org/docs/latest/development/new_module.html). ([#855](https://github.com/biocore/scikit-bio/issues/855))
* Updated to use `natsort` >= 4.0.0.
* The method of subsampling was changed for ``skbio.stats.power.subsample_paired_power``. Rather than drawing a paired sample for the run and then subsampling for each count, the subsample is now drawn for each sample and each run. In test data, this did not significantly alter the power results.
* checklist.py now enforces `__future__` imports in .py files.

## Version 0.2.3 (2015-02-13)

### Features
* Modified ``skbio.stats.distance.pwmantel`` to accept a list of filepaths. This is useful as it allows for a smaller amount of memory consumption as it only loads two matrices at a time as opposed to requiring that all distance matrices are loaded into memory.
* Added ``skbio.util.find_duplicates`` for finding duplicate elements in an iterable.

### Bug fixes
* Fixed floating point precision bugs in ``Alignment.position_frequencies``, ``Alignment.position_entropies``, ``Alignment.omit_gap_positions``, ``Alignment.omit_gap_sequences``, ``BiologicalSequence.k_word_frequencies``, and ``SequenceCollection.k_word_frequencies`` ([#801](https://github.com/biocore/scikit-bio/issues/801)).

### Backward-incompatible changes
* Removed ``feature_types`` attribute from ``BiologicalSequence`` and all subclasses ([#797](https://github.com/biocore/scikit-bio/pull/797)).
* Removed ``find_features`` method from ``BiologicalSequence`` and ``ProteinSequence`` ([#797](https://github.com/biocore/scikit-bio/pull/797)).
* ``BiologicalSequence.k_word_frequencies`` now returns a ``collections.defaultdict`` of type ``float`` instead of type ``int``. This only affects the "default" case, when a key isn't present in the dictionary. Previous behavior would return ``0`` as an ``int``, while the new behavior is to return ``0.0`` as a ``float``. This change also affects the ``defaultdict``s that are returned by ``SequenceCollection.k_word_frequencies``.

### Miscellaneous
* ``DissimilarityMatrix`` and ``DistanceMatrix`` now report duplicate IDs in the ``DissimilarityMatrixError`` message that can be raised during validation.

## Version 0.2.2 (2014-12-04)

### Features
* Added ``plot`` method to ``skbio.stats.distance.DissimilarityMatrix`` for creating basic heatmaps of a dissimilarity/distance matrix (see [#684](https://github.com/biocore/scikit-bio/issues/684)). Also added  ``_repr_png_`` and ``_repr_svg_`` methods for automatic display in the IPython Notebook, with ``png`` and ``svg`` properties for direct access.
* Added `__str__` method to `skbio.stats.ordination.OrdinationResults`.
* Added ``skbio.stats.distance.anosim`` and ``skbio.stats.distance.permanova`` functions, which replace the ``skbio.stats.distance.ANOSIM`` and ``skbio.stats.distance.PERMANOVA`` classes. These new functions provide simpler procedural interfaces to running these statistical methods. They also provide more convenient access to results by returning a ``pandas.Series`` instead of a ``CategoricalStatsResults`` object. These functions have more extensive documentation than their previous versions. If significance tests are suppressed, p-values are returned as ``np.nan`` instead of ``None`` for consistency with other statistical methods in scikit-bio. [#754](https://github.com/biocore/scikit-bio/issues/754)
* Added `skbio.stats.power` for performing empirical power analysis. The module uses existing datasets and iteratively draws samples to estimate the number of samples needed to see a significant difference for a given critical value.
* Added `skbio.stats.isubsample` for subsampling from an unknown number of values. This method supports subsampling from multiple partitions and does not require that all items be stored in memory, requiring approximately `O(N*M)`` space where `N` is the number of partitions and `M` is the maximum subsample size.
* Added ``skbio.stats.subsample_counts``, which replaces ``skbio.stats.subsample``. See deprecation section below for more details ([#770](https://github.com/biocore/scikit-bio/issues/770)).

### Bug fixes
* Fixed issue where SSW wouldn't compile on i686 architectures ([#409](https://github.com/biocore/scikit-bio/issues/409)).

### Deprecated functionality
* Deprecated ``skbio.stats.p_value_to_str``. This function will be removed in scikit-bio 0.3.0. Permutation-based p-values in scikit-bio are calculated as ``(num_extreme + 1) / (num_permutations + 1)``, so it is impossible to obtain a p-value of zero. This function historically existed for correcting the number of digits displayed when obtaining a p-value of zero. Since this is no longer possible, this functionality will be removed.
* Deprecated ``skbio.stats.distance.ANOSIM`` and ``skbio.stats.distance.PERMANOVA`` in favor of ``skbio.stats.distance.anosim`` and ``skbio.stats.distance.permanova``, respectively.
* Deprecated ``skbio.stats.distance.CategoricalStatsResults`` in favor of using ``pandas.Series`` to store statistical method results. ``anosim`` and ``permanova`` return ``pandas.Series`` instead of ``CategoricalStatsResults``.
* Deprecated ``skbio.stats.subsample`` in favor of ``skbio.stats.subsample_counts``, which provides an identical interface; only the function name has changed. ``skbio.stats.subsample`` will be removed in scikit-bio 0.3.0.

### Backward-incompatible changes
* Deprecation warnings are now raised using ``DeprecationWarning`` instead of ``UserWarning`` ([#774](https://github.com/biocore/scikit-bio/issues/774)).

### Miscellaneous
* The ``pandas.DataFrame`` returned by ``skbio.stats.distance.pwmantel`` now stores p-values as floats and does not convert them to strings with a specific number of digits. p-values that were previously stored as "N/A" are now stored as ``np.nan`` for consistency with other statistical methods in scikit-bio. See note in "Deprecated functionality" above regarding ``p_value_to_str`` for details.
* scikit-bio now supports versions of IPython < 2.0.0 ([#767](https://github.com/biocore/scikit-bio/issues/767)).

## Version 0.2.1 (2014-10-27)

This is an alpha release of scikit-bio. At this stage, major backwards-incompatible API changes can and will happen. Unified I/O with the scikit-bio I/O registry was the focus of this release.

### Features
* Added ``strict`` and ``lookup`` optional parameters to ``skbio.stats.distance.mantel`` for handling reordering and matching of IDs when provided ``DistanceMatrix`` instances as input (these parameters were previously only available in ``skbio.stats.distance.pwmantel``).
* ``skbio.stats.distance.pwmantel`` now accepts an iterable of ``array_like`` objects. Previously, only ``DistanceMatrix`` instances were allowed.
* Added ``plot`` method to ``skbio.stats.ordination.OrdinationResults`` for creating basic 3-D matplotlib scatterplots of ordination results, optionally colored by metadata in a ``pandas.DataFrame`` (see [#518](https://github.com/biocore/scikit-bio/issues/518)). Also added  ``_repr_png_`` and ``_repr_svg_`` methods for automatic display in the IPython Notebook, with ``png`` and ``svg`` properties for direct access.
* Added ``skbio.stats.ordination.assert_ordination_results_equal`` for comparing ``OrdinationResults`` objects for equality in unit tests.
* ``BiologicalSequence`` (and its subclasses) now optionally store Phred quality scores. A biological sequence's quality scores are stored as a 1-D ``numpy.ndarray`` of nonnegative integers that is the same length as the biological sequence. Quality scores can be provided upon object instantiation via the keyword argument ``quality``, and can be retrieved via the ``BiologicalSequence.quality`` property. ``BiologicalSequence.has_quality`` is also provided for determining whether a biological sequence has quality scores or not. See [#616](https://github.com/biocore/scikit-bio/issues/616) for more details.
* Added ``BiologicalSequence.sequence`` property for retrieving the underlying string representing the sequence characters. This was previously (and still is) accessible via ``BiologicalSequence.__str__``. It is provided via a property for convenience and explicitness.
* Added ``BiologicalSequence.equals`` for full control over equality testing of biological sequences. By default, biological sequences must have the same type, underlying sequence of characters, identifier, description, and quality scores to compare equal. These properties can be ignored via the keyword argument ``ignore``. The behavior of ``BiologicalSequence.__eq__``/``__ne__`` remains unchanged (only type and underlying sequence of characters are compared).
* Added ``BiologicalSequence.copy`` for creating a copy of a biological sequence, optionally with one or more attributes updated.
* ``BiologicalSequence.__getitem__`` now supports specifying a sequence of indices to take from the biological sequence.
* Methods to read and write taxonomies are now available under ``skbio.tree.TreeNode.from_taxonomy`` and ``skbio.tree.TreeNode.to_taxonomy`` respectively.
* Added ``SequenceCollection.update_ids``, which provides a flexible way of updating sequence IDs on a ``SequenceCollection`` or ``Alignment`` (note that a new object is returned, since instances of these classes are immutable). Deprecated ``SequenceCollection.int_map`` in favor of this new method; it will be removed in scikit-bio 0.3.0.
* Added ``skbio.util.cardinal_to_ordinal`` for converting a cardinal number to ordinal string (e.g., useful for error messages).
* New I/O Registry: supports multiple file formats, automatic file format detection when reading, unified procedural ``skbio.io.read`` and ``skbio.io.write`` in addition to OOP interfaces (``read/write`` methods) on the below objects. See ``skbio.io`` for more details.
    - Added "clustal" format support:
        * Has sniffer
        * Readers: ``Alignment``
        * Writers: ``Alignment``
    - Added "lsmat" format support:
        * Has sniffer
        * Readers: ``DissimilarityMatrix``, ``DistanceMatrix``
        * Writers: ``DissimilarityMatrix``, ``DistanceMatrix``
    - Added "ordination" format support:
        * Has sniffer
        * Readers: ``OrdinationResults``
        * Writers: ``OrdinationResults``
    - Added "newick" format support:
        * Has sniffer
        * Readers: ``TreeNode``
        * Writers: ``TreeNode``
    - Added "phylip" format support:
        * No sniffer
        * Readers: None
        * Writers: ``Alignment``
    - Added "qseq" format support:
        * Has sniffer
        * Readers: generator of ``BiologicalSequence`` or its subclasses, ``SequenceCollection``, ``BiologicalSequence``, ``NucleotideSequence``, ``DNASequence``, ``RNASequence``, ``ProteinSequence``
        * Writers: None
    - Added "fasta"/QUAL format support:
        * Has sniffer
        * Readers: generator of ``BiologicalSequence`` or its subclasses, ``SequenceCollection``, ``Alignment``, ``BiologicalSequence``, ``NucleotideSequence``, ``DNASequence``, ``RNASequence``, ``ProteinSequence``
        * Writers: same as readers
    - Added "fastq" format support:
        * Has sniffer
        * Readers: generator of ``BiologicalSequence`` or its subclasses, ``SequenceCollection``, ``Alignment``, ``BiologicalSequence``, ``NucleotideSequence``, ``DNASequence``, ``RNASequence``, ``ProteinSequence``
        * Writers: same as readers

### Bug fixes

* Removed ``constructor`` parameter from ``Alignment.k_word_frequencies``, ``BiologicalSequence.k_words``, ``BiologicalSequence.k_word_counts``, and ``BiologicalSequence.k_word_frequencies`` as it had no effect (it was never hooked up in the underlying code). ``BiologicalSequence.k_words`` now returns a generator of ``BiologicalSequence`` objects instead of strings.
* Modified the ``Alignment`` constructor to verify that all sequences have the same length, if not, raise an ``AlignmentError`` exception.  Updated the method ``Alignment.subalignment`` to calculate the indices only once now that identical sequence length is guaranteed.

### Deprecated functionality
* Deprecated ``constructor`` parameter in ``Alignment.majority_consensus`` in favor of having users call ``str`` on the returned ``BiologicalSequence``. This parameter will be removed in scikit-bio 0.3.0.

* Existing I/O functionality deprecated in favor of I/O registry, old functionality will be removed in scikit-bio 0.3.0. All functionality can be found at ``skbio.io.read``, ``skbio.io.write``, and the methods listed below:
    * Deprecated the following "clustal" readers/writers:
        - ``write_clustal`` -> ``Alignment.write``
        - ``parse_clustal`` -> ``Alignment.read``

    * Deprecated the following distance matrix format ("lsmat") readers/writers:
        - ``DissimilarityMatrix.from_file`` -> ``DissimilarityMatrix.read``
        - ``DissimilarityMatrix.to_file`` -> ``DissimilarityMatrix.write``
        - ``DistanceMatrix.from_file`` -> ``DistanceMatrix.read``
        - ``DistanceMatrix.to_file`` -> ``DistanceMatrix.write``

    * Deprecated the following ordination format ("ordination") readers/writers:
        - ``OrdinationResults.from_file`` -> ``OrdinationResults.read``
        - ``OrdinationResults.to_file`` -> ``OrdinationResults.write``

    * Deprecated the following "newick" readers/writers:
        - ``TreeNode.from_file`` -> ``TreeNode.read``
        - ``TreeNode.from_newick`` -> ``TreeNode.read``
        - ``TreeNode.to_newick`` -> ``TreeNode.write``

    * Deprecated the following "phylip" writers:
        - ``Alignment.to_phylip`` -> ``Alignment.write``

    * Deprecated the following "fasta"/QUAL readers/writers:
        - ``SequenceCollection.from_fasta_records`` -> ``SequenceCollection.read``
        - ``SequenceCollection.to_fasta`` -> ``SequenceCollection.write``
        - ``fasta_from_sequences`` -> ``skbio.io.write(obj, into=<file>, format='fasta')``
        - ``fasta_from_alignment`` -> ``Alignment.write``
        - ``parse_fasta`` -> ``skbio.io.read(<fasta>, format='fasta')``
        - ``parse_qual`` -> ``skbio.io.read(<fasta>, format='fasta', qual=<file>)``
        - ``BiologicalSequence.to_fasta`` -> ``BiologicalSequence.write``

    * Deprecated the following "fastq" readers/writers:
        - ``parse_fastq`` -> ``skbio.io.read(<fastq>, format='fastq')``
        - ``format_fastq_record`` -> ``skbio.io.write(<fastq>, format='fastq')``

### Backward-incompatible changes

* ``skbio.stats.distance.mantel`` now returns a 3-element tuple containing correlation coefficient, p-value, and the number of matching rows/cols in the distance matrices (``n``). The return value was previously a 2-element tuple containing only the correlation coefficient and p-value.
* ``skbio.stats.distance.mantel`` reorders input ``DistanceMatrix`` instances based on matching IDs (see optional parameters ``strict`` and ``lookup`` for controlling this behavior). In the past, ``DistanceMatrix`` instances were treated the same as ``array_like`` input and no reordering took place, regardless of ID (mis)matches. ``array_like`` input behavior remains the same.
* If mismatched types are provided to ``skbio.stats.distance.mantel`` (e.g., a ``DistanceMatrix`` and ``array_like``), a ``TypeError`` will be raised.

### Miscellaneous

* Added git timestamp checking to checklist.py, ensuring that when changes are made to Cython (.pyx) files, their corresponding generated C files are also updated.
* Fixed performance bug when instantiating ``BiologicalSequence`` objects. The previous runtime scaled linearly with sequence length; it is now constant time when the sequence is already a string. See [#623](https://github.com/biocore/scikit-bio/issues/623) for details.
* IPython and six are now required dependencies.

## Version 0.2.0 (2014-08-07)

This is an initial alpha release of scikit-bio. At this stage, major backwards-incompatible API changes can and will happen. Many backwards-incompatible API changes were made since the previous release.

### Features

* Added ability to compute distances between sequences in a ``SequenceCollection`` object ([#509](https://github.com/biocore/scikit-bio/issues/509)), and expanded ``Alignment.distance`` to allow the user to pass a function for computing distances (the default distance metric is still ``scipy.spatial.distance.hamming``) ([#194](https://github.com/biocore/scikit-bio/issues/194)).
* Added functionality to not penalize terminal gaps in global alignment. This functionality results in more biologically relevant global alignments (see [#537](https://github.com/biocore/scikit-bio/issues/537) for discussion of the issue) and is now the default behavior for global alignment.
* The python global aligners (``global_pairwise_align``, ``global_pairwise_align_nucleotide``, and ``global_pairwise_align_protein``) now support aligning pairs of sequences, pairs of alignments, and a sequence and an alignment (see [#550](https://github.com/biocore/scikit-bio/issues/550)). This functionality supports progressive multiple sequence alignment, among other things such as adding a sequence to an existing alignment.
* Added ``StockholmAlignment.to_file`` for writing Stockholm-formatted files.
* Added ``strict=True`` optional parameter to ``DissimilarityMatrix.filter``.
* Added ``TreeNode.find_all`` for finding all tree nodes that match a given name.


### Bug fixes

* Fixed bug that resulted in a ``ValueError`` from ``local_align_pairwise_nucleotide`` (see [#504](https://github.com/biocore/scikit-bio/issues/504)) under many circumstances. This would not generate incorrect results, but would cause the code to fail.

### Backward-incompatible changes

* Removed ``skbio.math``, leaving ``stats`` and ``diversity`` to become top level packages. For example, instead of ``from skbio.math.stats.ordination import PCoA`` you would now import ``from skbio.stats.ordination import PCoA``.
* The module ``skbio.math.gradient`` as well as the contents of ``skbio.math.subsample`` and ``skbio.math.stats.misc`` are now found in ``skbio.stats``. As an example, to import subsample: ``from skbio.stats import subsample``; to import everything from gradient: ``from skbio.stats.gradient import *``.
* The contents of ``skbio.math.stats.ordination.utils`` are now in ``skbio.stats.ordination``.
* Removed ``skbio.app`` subpackage (i.e., the *application controller framework*) as this code has been ported to the standalone [burrito](https://github.com/biocore/burrito) Python package. This code was not specific to bioinformatics and is useful for wrapping command-line applications in general.
* Removed ``skbio.core``, leaving ``alignment``, ``genetic_code``, ``sequence``, ``tree``, and ``workflow`` to become top level packages. For example, instead of ``from skbio.core.sequence import DNA`` you would now import ``from skbio.sequence import DNA``.
* Removed ``skbio.util.exception`` and ``skbio.util.warning`` (see [#577](https://github.com/biocore/scikit-bio/issues/577) for the reasoning behind this change). The exceptions/warnings were moved to the following locations:
 - ``FileFormatError``, ``RecordError``, ``FieldError``, and ``EfficiencyWarning`` have been moved to ``skbio.util``
 - ``BiologicalSequenceError`` has been moved to ``skbio.sequence``
 - ``SequenceCollectionError`` and ``StockholmParseError`` have been moved to ``skbio.alignment``
 - ``DissimilarityMatrixError``, ``DistanceMatrixError``, ``DissimilarityMatrixFormatError``, and ``MissingIDError`` have been moved to ``skbio.stats.distance``
 - ``TreeError``, ``NoLengthError``, ``DuplicateNodeError``, ``MissingNodeError``, and ``NoParentError`` have been moved to ``skbio.tree``
 - ``FastqParseError`` has been moved to ``skbio.parse.sequences``
 - ``GeneticCodeError``, ``GeneticCodeInitError``, and ``InvalidCodonError`` have been moved to ``skbio.genetic_code``
* The contents of ``skbio.genetic_code`` formerly ``skbio.core.genetic_code`` are now in ``skbio.sequence``. The ``GeneticCodes`` dictionary is now a function ``genetic_code``. The functionality is the same, except that because this is now a function rather than a dict, retrieving a genetic code is done using a function call rather than a lookup (so, for example, ``GeneticCodes[2]`` becomes ``genetic_code(2)``.
* Many submodules have been made private with the intention of simplifying imports for users. See [#562](https://github.com/biocore/scikit-bio/issues/562) for discussion of this change. The following list contains the previous module name and where imports from that module should now come from.
 - ``skbio.alignment.ssw`` to ``skbio.alignment``
 - ``skbio.alignment.alignment`` to ``skbio.alignment``
 - ``skbio.alignment.pairwise`` to ``skbio.alignment``
 - ``skbio.diversity.alpha.base`` to ``skbio.diversity.alpha``
 - ``skbio.diversity.alpha.gini`` to ``skbio.diversity.alpha``
 - ``skbio.diversity.alpha.lladser`` to ``skbio.diversity.alpha``
 - ``skbio.diversity.beta.base`` to ``skbio.diversity.beta``
 - ``skbio.draw.distributions`` to ``skbio.draw``
 - ``skbio.stats.distance.anosim`` to ``skbio.stats.distance``
 - ``skbio.stats.distance.base`` to ``skbio.stats.distance``
 - ``skbio.stats.distance.permanova`` to ``skbio.stats.distance``
 - ``skbio.distance`` to ``skbio.stats.distance``
 - ``skbio.stats.ordination.base`` to ``skbio.stats.ordination``
 - ``skbio.stats.ordination.canonical_correspondence_analysis`` to ``skbio.stats.ordination``
 - ``skbio.stats.ordination.correspondence_analysis`` to ``skbio.stats.ordination``
 - ``skbio.stats.ordination.principal_coordinate_analysis`` to ``skbio.stats.ordination``
 - ``skbio.stats.ordination.redundancy_analysis`` to ``skbio.stats.ordination``
 - ``skbio.tree.tree`` to ``skbio.tree``
 - ``skbio.tree.trie`` to ``skbio.tree``
 - ``skbio.util.misc`` to ``skbio.util``
 - ``skbio.util.testing`` to ``skbio.util``
 - ``skbio.util.exception`` to ``skbio.util``
 - ``skbio.util.warning`` to ``skbio.util``
* Moved ``skbio.distance`` contents into ``skbio.stats.distance``.

### Miscellaneous

* Relaxed requirement in ``BiologicalSequence.distance`` that sequences being compared are of equal length. This is relevant for Hamming distance, so the check is still performed in that case, but other distance metrics may not have that requirement. See [#504](https://github.com/biocore/scikit-bio/issues/507)).
* Renamed ``powertrip.py`` repo-checking script to ``checklist.py`` for clarity.
* ``checklist.py`` now ensures that all unit tests import from a minimally deep API. For example, it will produce an error if ``skbio.core.distance.DistanceMatrix`` is used over ``skbio.DistanceMatrix``.
* Extra dimension is no longer calculated in ``skbio.stats.spatial.procrustes``.
* Expanded documentation in various subpackages.
* Added new scikit-bio logo. Thanks [Alina Prassas](http://cargocollective.com/alinaprassas)!

## Version 0.1.4 (2014-06-25)

This is a pre-alpha release. At this stage, major backwards-incompatible API changes can and will happen.

### Features

* Added Python implementations of Smith-Waterman and Needleman-Wunsch alignment as ``skbio.core.alignment.pairwise.local_pairwise_align`` and ``skbio.core.alignment.pairwise.global_pairwise_align``. These are much slower than native C implementations (e.g., ``skbio.core.alignment.local_pairwise_align_ssw``) and as a result raise an ``EfficencyWarning`` when called, but are included as they serve as useful educational examples as they’re simple to experiment with.
* Added ``skbio.core.diversity.beta.pw_distances`` and ``skbio.core.diversity.beta.pw_distances_from_table``. These provide convenient access to the ``scipy.spatial.distance.pdist`` *beta diversity* metrics from within scikit-bio. The ``skbio.core.diversity.beta.pw_distances_from_table`` function will only be available temporarily, until the ``biom.table.Table`` object is merged into scikit-bio (see [#489](https://github.com/biocore/scikit-bio/issues/489)), at which point ``skbio.core.diversity.beta.pw_distances`` will be updated to use that.
* Added ``skbio.core.alignment.StockholmAlignment``, which provides support for parsing [Stockholm-formatted alignment files](http://sonnhammer.sbc.su.se/Stockholm.html) and working with those alignments in the context RNA secondary structural information.
* Added ``skbio.core.tree.majority_rule`` function for computing consensus trees from a list of trees.

### Backward-incompatible changes

* Function ``skbio.core.alignment.align_striped_smith_waterman`` renamed to ``local_pairwise_align_ssw`` and now returns an ``Alignment`` object instead of an ``AlignmentStructure``
* The following keyword-arguments for ``StripedSmithWaterman`` and ``local_pairwise_align_ssw`` have been renamed:
    * ``gap_open`` -> ``gap_open_penalty``
    * ``gap_extend`` -> ``gap_extend_penalty``
    * ``match`` -> ``match_score``
    * ``mismatch`` -> ``mismatch_score``
* Removed ``skbio.util.sort`` module in favor of [natsort](https://pypi.python.org/pypi/natsort) package.

### Miscellaneous

* Added powertrip.py script to perform basic sanity-checking of the repo based on recurring issues that weren't being caught until release time; added to Travis build.
* Added RELEASE.md with release instructions.
* Added intersphinx mappings to docs so that "See Also" references to numpy, scipy, matplotlib, and pandas are hyperlinks.
* The following classes are no longer ``namedtuple`` subclasses (see [#359](https://github.com/biocore/scikit-bio/issues/359) for the rationale):
    * ``skbio.math.stats.ordination.OrdinationResults``
    * ``skbio.math.gradient.GroupResults``
    * ``skbio.math.gradient.CategoryResults``
    * ``skbio.math.gradient.GradientANOVAResults``
* Added coding guidelines draft.
* Added new alpha diversity formulas to the ``skbio.math.diversity.alpha`` documentation.

## Version 0.1.3 (2014-06-12)

This is a pre-alpha release. At this stage, major backwards-incompatible API changes can and will happen.

### Features

* Added ``enforce_qual_range`` parameter to ``parse_fastq`` (on by default, maintaining backward compatibility). This allows disabling of the quality score range-checking.
* Added ``skbio.core.tree.nj``, which applies neighbor-joining for phylogenetic reconstruction.
* Added ``bioenv``, ``mantel``, and ``pwmantel`` distance-based statistics to ``skbio.math.stats.distance`` subpackage.
* Added ``skbio.math.stats.misc`` module for miscellaneous stats utility functions.
* IDs are now optional when constructing a ``DissimilarityMatrix`` or ``DistanceMatrix`` (monotonically-increasing integers cast as strings are automatically used).
* Added ``DistanceMatrix.permute`` method for randomly permuting rows and columns of a distance matrix.
* Added the following methods to ``DissimilarityMatrix``: ``filter``, ``index``, and ``__contains__`` for ID-based filtering, index lookup, and membership testing, respectively.
* Added ``ignore_comment`` parameter to ``parse_fasta`` (off by default, maintaining backward compatibility). This handles stripping the comment field from the header line (i.e., all characters beginning with the first space) before returning the label.
* Added imports of ``BiologicalSequence``, ``NucleotideSequence``, ``DNA``, ``DNASequence``, ``RNA``, ``RNASequence``, ``Protein``, ``ProteinSequence``, ``DistanceMatrix``, ``align_striped_smith_waterman``, `` SequenceCollection``, ``Alignment``, ``TreeNode``, ``nj``, ``parse_fasta``, ``parse_fastq``, ``parse_qual``, ``FastaIterator``, ``FastqIterator``, ``SequenceIterator`` in ``skbio/__init__.py`` for convenient importing. For example, it's now possible to ``from skbio import Alignment``, rather than ``from skbio.core.alignment import Alignment``.

### Bug fixes

* Fixed a couple of unit tests that could fail stochastically.
* Added missing ``__init__.py`` files to a couple of test directories so that these tests won't be skipped.
* ``parse_fastq`` now raises an error on dangling records.
* Fixed several warnings that were raised while running the test suite with Python 3.4.

### Backward-incompatible changes

* Functionality imported from ``skbio.core.ssw`` must now be imported from ``skbio.core.alignment`` instead.

### Miscellaneous

* Code is now flake8-compliant; added flake8 checking to Travis build.
* Various additions and improvements to documentation (API, installation instructions, developer instructions, etc.).
* ``__future__`` imports are now standardized across the codebase.
* New website front page and styling changes throughout. Moved docs site to its own versioned subdirectories.
* Reorganized alignment data structures and algorithms (e.g., SSW code, ``Alignment`` class, etc.) into an ``skbio.core.alignment`` subpackage.

## Version 0.1.1 (2014-05-16)

Fixes to setup.py. This is a pre-alpha release. At this stage, major backwards-incompatible API changes can and will happen.

## Version 0.1.0 (2014-05-15)

Initial pre-alpha release. At this stage, major backwards-incompatible API changes can and will happen.
# Reviewing pull requests

This document provides a high-level, general, and **incomplete** checklist of things that pull request reviewers should be aware of when performing code review. These are guidelines which generally make sense to follow, but they are not intended to be rigid. The checklist mainly consists of things that are specific to the scikit-bio project and that generally apply to incoming pull requests. The checklist is incomplete because it is not possible to describe all things to verify during code review (that depends on what is being reviewed). This document also doesn't attempt to describe *how* to perform code review (there are many online resources for that).

Reviewers are encouraged to keep this document up-to-date as the project evolves and to add anything that's missing.

The checklist is not in any particular order.

## Licensing and attribution

Verify that code being included from external sources has a compatible license and is properly attributed:

- Include the code's license in the top-level `licenses` directory.
- Include a comment with the external code giving attribution and noting the license in the `licenses` directory.
- Any other requirements set by the code's license and/or author.

## Changelog

This is one of the most important points to remember as users will review the changelog to identify changes relevant to them. This is one of the easiest parts to forget in a pull request.

- Note all public (i.e. user-facing) changes in `CHANGELOG.md` under the latest development version section at the top of the file. This includes things like bug fixes, API additions/changes/removal, performance enhancements, etc. The changelog has several subsections for organizing these changes.
- If a corresponding issue exists, it should be linked to from the changelog.
- Use public imports (`skbio.sequence.Sequence` instead of `skbio.sequence._sequence.Sequence`) when documenting import paths in the changelog.
- Internal changes not visible/applicable to users (e.g. refactoring, private methods, etc.) are better suited for commit messages than the changelog.

## Public vs. private API

Be aware of what type of API changes are being made. scikit-bio uses the following conventions to distinguish public vs. private API:

**Public API:** API with an import path that doesn't contain leading underscores in any submodule/subpackage/object names. Users and scikit-bio devs can use public API. Example: `skbio.sequence.Sequence`, `skbio.stats.subsample`

**Package-private API:** API with an import path containing a leading underscore in at least one submodule/subpackage. Users should not import package-private API. Package-private API can be imported anywhere *within* the scikit-bio package. Example: `skbio.util._misc.chunk_str`, `skbio.util._testing.ReallyEqualMixin`

**Private API:** API with object name starting with an underscore. Users should not import private API. Private API should only be imported and used in the module where it is defined. It should not be imported in other parts of the scikit-bio package. Example: `skbio.util._testing._normalize_signs`, `skbio.stats.composition._gram_schmidt_basis`

- Prefer defining private APIs and only promote things to the public API that users need access to.
- Private/package-private API does not need to be decorated with `@experimental`/`@stable`/`@deprecated`, only public API.
- Within scikit-bio, prefer *using* public APIs defined in the package even if private APIs offer the same functionality.

## API stability

Refer to scikit-bio's [API stability docs](http://scikit-bio.org/docs/latest/user/api_stability.html) for the current API lifecycle.

- Prefer making new APIs experimental. This allows for changes to the API without needing to deprecate it first.
- Try to avoid changing stable API if at all possible. If a change is necessary, deprecate the API and remove after 2+ minor releases.
- Under extreme circumstances, a stable API can be changed without deprecation warning to users. This has happened in the past to fix bugs that required an API change. If this happens, discuss with other devs before committing to the change and document this change extensively in the changelog.
- Always document stable/experimental API changes and any deprecated APIs in the changelog.

## Integration/consistency with existing API

When reviewing API changes/additions, look at them in the context of the rest of the codebase and its existing APIs. For example, are there existing parameter names that could be reused for consistency/predictability? Does the new API (or changed API) compose with other relevant APIs? Example: `skbio.stats.distance.anosim` uses a `permutations` parameter, so a new nonparametric function would want to use this name over something like `n` or `num_permutations`.

## Commit messages and merging pull requests

When merging pull requests, use GitHub's "Squash and merge" option to merge a single commit. See [this commit message](https://github.com/biocore/scikit-bio/commit/f3d736aabd717971332781b98d8fde861f354dc3) for an example.

- Rewrite commit message to describe all changes being merged. This usually involves deleting individual commit messages that GitHub includes in the text box.
- Include "fixes #n" text if there's an associated issue to be closed.
- Use [numpy-style commit tags](https://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html#writing-the-commit-message) (ENH, BUG, PERF, etc.).
- Include contributors' and reviewers' GitHub usernames in commit message (attribution will be lost on squash).

## Test changes locally

**This step is extremely important.** Pull down the PR changes locally and try out the API as a user would. Try to break it, make sure the docs are complete, etc. Build the docs locally and verify that they render correctly (this is a common pitfall).

## Docs

- Verify the docs follow the instructions in the scikit-bio [documentation guide](https://github.com/biocore/scikit-bio/blob/master/doc/README.md).
- Verify the docs follow [this page](http://scikit-bio.org/docs/latest/development/new_module.html) when adding a new module or subpackage to scikit-bio.
- Public API should have docstrings conforming to [numpydoc standards](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt). Manual and careful verification of the numpydoc docstrings is currently necessary; they are easy to get wrong and building the docs won't always flag issues. Building the docs and inspecting the rendered output can help with this process.
- Package-private and private APIs do not need to be extensively documented; numpydoc docstrings are not required. Document these APIs as appropriate to help other devs understand the code (code comments are usually better for this anyways).

## CI

- Make sure Travis-CI is passing before merging a pull request.
- Make sure coverage doesn’t drop. Strive to have 100% coverage for new code being merged.

## Unit testing

- Make sure the tests are as complete as possible.
- Check that border cases are tested (e.g. zeros, '', [], None, etc.).
- Check that the base case is tested (`n`), along with the inductive step (`n + 1`).
- Verify that tests cover more than one input data set.
- Make each test case simple, ideally only testing a single thing (follow [Arrange Act Assert](http://wiki.c2.com/?ArrangeActAssert)).
Contributing to scikit-bio
==========================

[scikit-bio](http://scikit-bio.org) is an open source software package and we welcome community contributions. You can find the scikit-bio source code on GitHub [here](https://github.com/biocore/scikit-bio).

This document covers what you should do to get started with contributing to scikit-bio. You should read the entire document before contributing code to scikit-bio. This will save time for both you and the scikit-bio developers.

Types of contributions
----------------------

We're interested in many different types of contributions, including feature additions, bug fixes, continuous integration improvements, and documentation/website updates, additions, and fixes.

When considering contributing to scikit-bio, you should begin by posting an issue to the [scikit-bio issue tracker](https://github.com/biocore/scikit-bio/issues). The information that you include in that post will differ based on the type of contribution. Your contribution will also need to be fully tested where applicable (discussed further below).

* For feature additions, please describe why the functionality that you are proposing to add is relevant. For it to be relevant, it should be demonstrably useful to scikit-bio users and it should also fit within the biology/bioinformatics domain. This typically means that a new analytic method is implemented (you should describe why it's useful, ideally including a link to a paper that uses this method), or an existing method is enhanced (e.g., improved performance). We will request benchmark results comparing your method to the pre-existing methods (which would also be required for publication of your method) so pointing to a paper or other document containing benchmark results, or including benchmark results in your issue, will speed up the process. Before contributing a new feature, it's also a good idea to check whether the functionality exists in other Python packages, or if the feature would fit better in another Python package. For example, low-level statistical methods/tests may fit better in a project that is focused on statistics (e.g., [SciPy](http://scipy.org/) or [statsmodels](http://statsmodels.sourceforge.net/)).

* For bug fixes, please provide a detailed description of the bug so other developers can reproduce it. We take bugs in scikit-bio very seriously. Bugs can be related to errors in code, documentation, or tests. Errors in documentation or tests are usually updated in the next scheduled release of scikit-bio. Errors in code that could result in incorrect results or inability to access certain functionality may result in a bug fix release of scikit-bio that is released ahead of schedule.

 You should include the following information in your bug report:

 1. The exact command(s) necessary to reproduce the bug.
 2. A link to all necessary input files for reproducing the bug. These files should only be as large as necessary to create the bug. For example, if you have an input file with 10,000 FASTA-formatted sequences but the error only arises due to one of the sequences, create a new FASTA file with only that sequence, run the command that was giving you problems, and verify that you still get an error. Then post that command and link to the trimmed FASTA file. This is *extremely* useful to other developers and it is likely that if you don't provide this information you'll get a response asking for it. Often this process helps you to better understand the bug as well.

* For documentation additions, you should first post an issue describing what you propose to add, where you'd like to add it in the documentation, and a description of why you think it's an important addition. For documentation improvements and fixes, you should post an issue describing what is currently wrong or missing and how you propose to address it. For more information about building and contributing to scikit-bio's documentation, see our [documentation guide](doc/README.md).

When you post your issue, the scikit-bio developers will respond to let you know if we agree with the addition or change. It's very important that you go through this step to avoid wasting time working on a feature that we are not interested in including in scikit-bio. **This initial discussion with the developers is important because scikit-bio is rapidly changing, including complete re-writes of some of the core objects. If you don't get in touch first you could easily waste time by working on an object or interface that is deprecated.**

Getting started
---------------

### "quick fixes"

Some of our issues are labeled as ``quick fix``. Working on [these issues](https://github.com/biocore/scikit-bio/issues?q=is%3Aopen+is%3Aissue+label%3A%22quick+fix%22) is a good way to get started with contributing to scikit-bio. These are usually small bugs or documentation errors that will only require one or a few lines of code to fix. Getting started by working on one of these issues will allow you to familiarize yourself with our development process before committing to a large amount of work (e.g., adding a new feature to scikit-bio). Please post a comment on the issue if you're interested in working on one of these "quick fixes".

### Joining development

Once you are more comfortable with our development process, you can check out the [``on deck`` label](https://github.com/biocore/scikit-bio/labels/on%20deck) on our issue tracker. These issues represent what our current focus is in the project. As such, they are probably the best place to start if you are looking to join the conversation and contribute code.

Code review
-----------

When you submit code to scikit-bio, it will be reviewed by one or more scikit-bio developers. These reviews are intended to confirm a few points:

* Your code provides relevant changes or additions to scikit-bio ([Types of contributions](#types-of-contributions)).
* Your code adheres to our coding guidelines ([Coding guidelines](#coding-guidelines)).
* Your code is sufficiently well-tested ([Testing guidelines](#testing-guidelines)).
* Your code is sufficiently well-documented ([Documentation guidelines](#documentation-guidelines)).

This process is designed to ensure the quality of scikit-bio and can be a very useful experience for new developers.

Particularly for big changes, if you'd like feedback on your code in the form of a code review as you work, you should request help in the issue that you created and one of the scikit-bio developers will work with you to perform regular code reviews. This can greatly reduce development time (and frustration) so we highly recommend that new developers take advantage of this rather than submitting a pull request with a massive amount of code. That can lead to frustration when the developer thinks they are done but the reviewer requests large amounts of changes, and it also makes it harder to review.

Submitting code to scikit-bio
-----------------------------

scikit-bio is hosted on [GitHub](http://www.github.com), and we use GitHub's [Pull Request](https://help.github.com/articles/using-pull-requests) mechanism for reviewing and accepting submissions. You should work through the following steps to submit code to scikit-bio.

1. Begin by [creating an issue](https://github.com/biocore/scikit-bio/issues) describing your proposed change (see [Types of contributions](#types-of-contributions) for details).

2. [Fork](https://help.github.com/articles/fork-a-repo) the scikit-bio repository on the GitHub website.

3. Clone your forked repository to the system where you'll be developing with ``git clone``. ``cd`` into the ``scikit-bio`` directory that was created by ``git clone``.

4. Ensure that you have the latest version of all files. This is especially important if you cloned a long time ago, but you'll need to do this before submitting changes regardless. You should do this by adding scikit-bio as a remote repository and then pulling from that repository. You'll only need to run the ``git remote`` command the first time you do this:

 ```
 git remote add upstream https://github.com/biocore/scikit-bio.git
 git checkout master
 git pull upstream master
 ```

5. Install scikit-bio for development. See [Setting up a development environment](#setting-up-a-development-environment).

6. Create a new topic branch that you will make your changes in with ``git checkout -b``:

 ```
 git checkout -b my-topic-branch
 ```

 What you name your topic branch is up to you, though we recommend including the issue number in the topic branch, since there is usually already an issue associated with the changes being made in the pull request. For example, if you were addressing issue number 42, you might name your topic branch ``issue-42``.

7. Run ``make test`` to confirm that the tests pass before you make any changes.

8. Make your changes, add them (with ``git add``), and commit them (with ``git commit``). Don't forget to update associated tests and documentation as necessary. Write descriptive commit messages to accompany each commit. We recommend following [NumPy's commit message guidelines](http://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html#writing-the-commit-message), including the usage of commit tags (i.e., starting commit messages with acronyms such ``ENH``, ``BUG``, etc.).

9. Please mention your changes in [CHANGELOG.md](CHANGELOG.md). This file informs scikit-bio *users* of changes made in each release, so be sure to describe your changes with this audience in mind. It is especially important to note API additions and changes, particularly if they are backward-incompatible, as well as bug fixes. Be sure to make your updates under the section designated for the latest development version of scikit-bio (this will be at the top of the file). Describe your changes in detail under the most appropriate section heading(s). For example, if your pull request fixes a bug, describe the bug fix under the "Bug fixes" section of [CHANGELOG.md](CHANGELOG.md). Please also include a link to the issue(s) addressed by your changes. See [CHANGELOG.md](CHANGELOG.md) for examples of how we recommend formatting these descriptions.

10. When you're ready to submit your code, ensure that you have the latest version of all files in case some changed while you were working on your edits. You can do this by merging master into your topic branch:

 ```
 git checkout master
 git pull upstream master
 git checkout my-topic-branch
 git merge master
 ```

11. Run ``make test`` to ensure that your changes did not cause anything expected to break.

12. Once the tests pass, you should push your changes to your forked repository on GitHub using:

 ```
 git push origin my-topic-branch
 ```

13. Issue a [pull request](https://help.github.com/articles/using-pull-requests) on the GitHub website to request that we merge your branch's changes into scikit-bio's master branch. Be sure to include a description of your changes in the pull request, as well as any other information that will help the scikit-bio developers involved in reviewing your code. Please include ``fixes #<issue-number>`` in your pull request description or in one of your commit messages so that the corresponding issue will be closed when the pull request is merged (see [here](https://help.github.com/articles/closing-issues-via-commit-messages/) for more details). One of the scikit-bio developers will review your code at this stage. If we request changes (which is very common), *don't issue a new pull request*. You should make changes on your topic branch, and commit and push them to GitHub. Your pull request will update automatically.

Setting up a development environment
------------------------------------

**Note:** scikit-bio must be developed in a Python 3.4 or later environment.

The recommended way to set up a development environment for contributing to scikit-bio is using [Anaconda](https://store.continuum.io/cshop/anaconda/) by Continuum Analytics, with its associated command line utility `conda`. The primary benefit of `conda` over `pip` is that on some operating systems (ie Linux), `pip` installs packages from source. This can take a very long time to install Numpy, scipy, matplotlib, etc. `conda` installs these packages using pre-built binaries, so the installation is much faster. Another benefit of `conda` is that it provides both package and environment management, which removes the necessity of using `virtualenv` separately. Not all packages are available using `conda`, therefore our strategy is to install as many packages as possible using `conda`, then install any remaining packages using `pip`.

1. Install Anaconda

 See [Continuum's site](https://store.continuum.io/cshop/anaconda/) for instructions. [Miniconda](http://conda.pydata.org/docs/install/quick.html) provides a fast way to get conda up and running.

2. Create a new conda environment
 ```
 conda create -n env_name python=3.4 pip
 ```

 Note that `env_name` can be any name desired, for example

 ```
 conda create -n skbio python=3.4 pip
 ```

3. Activate the environment

 This may be slightly different depending on the operating system. Refer to the Continuum site to find instructions for your OS.
 ```
 source activate env_name
 ```

4. Navigate to the scikit-bio directory
 See [the section on submitting code](#submitting-code-to-scikit-bio).
 ```
 cd /path/to/scikit-bio
 ```

5. Install `conda` requirements
 ```
 conda install --file ci/conda_requirements.txt
 ```

6. Install `pip` requirements
 ```
 pip install -r ci/pip_requirements.txt
 ```

7. Install scikit-bio
 ```
 pip install --no-deps -e .
 ```

8. Test the installation
 ```
 make test
 ```

Coding guidelines
-----------------

We adhere to the [PEP 8](http://www.python.org/dev/peps/pep-0008/) Python style guidelines. Please see scikit-bio's [coding guidelines](http://scikit-bio.org/docs/latest/development/coding_guidelines.html) for more details. Before submitting code to scikit-bio, you should read this document carefully and apply the guidelines in your code.

Testing guidelines
------------------

All code that is added to scikit-bio must be unit tested, and the unit test code must be submitted in the same pull request as the library code that you are submitting. We will only merge code that is unit tested and that passes the [continuous integration build](https://github.com/biocore/scikit-bio/blob/master/.travis.yml). This build includes, but is not limited to, the following checks:

- Full unit test suite and doctests execute without errors in supported versions of Python 3.
- C code can be correctly compiled.
- Cython code is correctly generated.
- All tests import functionality from the appropriate minimally deep API.
- Documentation can be built.
- Current code coverage is maintained or improved.
- Code passes ``flake8`` checks.

Running ``make test`` locally during development will include a subset of the full checks performed by Travis-CI.

The scikit-bio coding guidelines describe our [expectations for unit tests](http://scikit-bio.org/docs/latest/development/coding_guidelines.html). You should review the unit test section before working on your test code.

Tests can be executed by running ``make test`` from the base directory of the project or from within a Python or IPython session:

``` python
>>> from skbio.test import pytestrunner
>>> pytestrunner()
# full test suite is executed
```

Documentation guidelines
------------------------

We strive to keep scikit-bio well-documented, particularly its public-facing API. See our [documentation guide](doc/README.md) for more details.

Getting help with git
---------------------

If you're new to ``git``, you'll probably find [gitref.org](http://gitref.org/) helpful.
# Releasing a new version of scikit-bio

## Introduction

This guide explains how to release a new version of scikit-bio. To illustrate examples of commands you might run, let's assume that the current version is x.y.y-dev and we want to release version x.y.z.

**Note:** The following commands assume you are in the top-level directory of the scikit-bio repository unless otherwise noted. They also assume that you have [Miniconda3](http://conda.pydata.org/miniconda.html) installed. It is important that you use Miniconda3, and not Miniconda2, because the root `conda` environment needs Python 3 for the build steps below.

## Prep the release

1. Ensure the Travis build is passing against master.

2. Update the version strings (x.y.y-dev) to the new version (x.y.z). This will include `__version__` defined in ``skbio/__init__.py``, as well as any `@experimental/@stable/@deprecated` [API stability decorators](http://scikit-bio.org/docs/latest/user/api_stability.html) with `as_of='x.y.y-dev'`. ``grep`` for the current version string to find all occurrences:

        grep -r 'x\.y\.y-dev' .

3. Remove any deprecated functionality that was scheduled for removal on or before this release. When removing deprecated functionality, make sure the functionality has been in a deprecated state for the appropriate number of releases described in the [API stability docs](http://scikit-bio.org/docs/latest/user/api_stability.html). If there is functionality that shouldn't be removed yet, bump the `until` version to a future version. To find all deprecated functionality, search for `@deprecated` decorators:

        grep -r '@deprecated' .

    Note any deprecated functionality that was removed in the `### Miscellaneous` section of `CHANGELOG.md`.

4. Update ``CHANGELOG.md`` to include descriptions of all **user-facing** changes that made it into this release. Be sure to update the heading to include the new version (x.y.z) and the date of the release. Use the existing structure in the file as a template/guide.

5. Submit a pull request and merge when Travis-CI tests are passing.

## Build website docs

You will need to **fully install** the latest master branch of scikit-bio (including built extensions) and build the docs from this version. **Make sure the version of scikit-bio that is imported by ``import skbio`` is the correct one before building the docs.**

1. Build the documentation locally:

        make -C doc clean html

2. Switch to the ``gh-pages`` branch of the repository.

3. Remove ``docs/latest``:

        git rm -rf docs/latest

4. Copy over the built documentation to ``docs/x.y.z`` and ``docs/latest``:

        cp -r doc/build/html docs/latest
        cp -r doc/build/html docs/x.y.z

5. Add a new list item to ``index.html`` to link to ``docs/x.y.z/index.html``.

6. Port content from ``README.md`` to ``index.html`` if there are any changes that need to be included on the front page.

7. Test out your changes by opening the site locally in a browser. Be sure to check the error console for any errors.

8. Submit a pull request with the website updates and merge. **Note:** Once merged, the live website is updated, so be sure to poke through the live site to make sure things are rendered correctly with the right version strings.

## Tag the release

From the [scikit-bio GitHub page](https://github.com/biocore/scikit-bio), click on the releases tab and draft a new release. Use the version number for the tag name (x.y.z) and create the tag against master. Fill in a release title that is consistent with previous release titles and add a summary of the release (linking to ``CHANGELOG.md`` is a good idea). This release summary will be the primary information that we point users to when we announce the release.

Once the release is created on GitHub, it's a good idea to test out the release tarball before publishing to PyPI:

1. Create a new `conda` environment for testing (fill in a name for `<environment>`):

        conda create -n <environment> python=3.5 numpy
        source activate <environment>

2. Install the release tarball from GitHub and run the tests:

        pip install https://github.com/biocore/scikit-bio/archive/x.y.z.tar.gz
        python -m skbio.test

## Publish the release

Assuming the GitHub release tarball correctly installs and passes its tests, you're ready to create the source distribution (``sdist``) that will be published to PyPI. It is important to test the source distribution because it is created in an entirely different way than the release tarball on GitHub. Thus, there is the danger of having two different release tarballs: the one created on GitHub and the one uploaded to PyPI.

1. Download the release tarball from GitHub, extract it, and ``cd`` into the top-level directory.

2. Build a source distribution:

        python setup.py sdist

3. Create and activate a new `conda` environment, and test the `sdist`:

        pip install dist/scikit-bio-x.y.z.tar.gz
        cd  # cd somewhere outside the extracted scikit-bio directory
        python -m skbio.test

4. If everything goes well, it is finally time to push the release to PyPI:

        python setup.py sdist upload

    You must have the proper login credentials to add a release to PyPI. Currently [@gregcaporaso](https://github.com/gregcaporaso) has these, but they can be shared with other release managers.

5. Once the release is available on PyPI, do a final round of testing. Create a new `conda` environment and run:

        pip install scikit-bio
        cd  # cd somewhere outside the extracted scikit-bio directory
        python -m skbio.test

    If this succeeds, the PyPI release appears to be a success. Make sure the installed version is the correct one.

6. Next, we'll prepare and post the release to [anaconda.org](http://www.anaconda.org).

    You'll need to have ``conda-build`` and ``anaconda-client`` installed to perform these steps. Both can be conda-installed. First, log into anaconda with your anaconda username using the following command. You should have access to push to the ``biocore`` anaconda account through your account (if you don't, get in touch with [@gregcaporaso](https://github.com/gregcaporaso) who is the owner of that account).

        anaconda login

    Due to its C extensions, releasing scikit-bio packages for different platforms will require you to perform the following steps on each of those platforms. For example, an ``osx-64`` package will need to be built on OS X, and a ``linux-64`` package will need to be built on 64-bit Linux. These steps will be the same on all platforms, so you should repeat them for every platform you want to release for.

        conda skeleton pypi scikit-bio
        conda build scikit-bio --python 3.4
        conda build scikit-bio --python 3.5

    **Note:** When building 64-bit Linux packages, it is recommended that you use conda-forge's `linux-anvil` Docker image. This ensures a consistent Linux build environment that has an old enough version of `libc` to be compatible on most Linux systems. To start up a `linux-anvil` Docker container:

        docker run -i -t condaforge/linux-anvil
        # Now you should be in the linux-anvil environment
        sed -i '/conda-forge/d' ~/.condarc
        # Run the build commands from above

    At this stage you have built Python 3.4 and 3.5 packages. The absolute path to the packages will be provided as output from each ``conda build`` commands. You should now create conda environments for each, and run the tests as described above. You can install these local packages as follows:

        conda install --use-local scikit-bio

    If the tests pass, you're ready to upload.

        anaconda upload -u biocore <package-filepath>

    ``<package-filepath>`` should be replaced with the path to the package that was was created above. Repeat this for each package you created (here, the Python 3.4 and 3.5 packages).

    After uploading, you should create new environments for every package you uploaded, install scikit-bio from each package, and re-run the tests. You can install the packages you uploaded as follows:

        conda install -c https://conda.anaconda.org/biocore scikit-bio

## Post-release cleanup

1. Submit and merge a pull request to update the version strings from x.y.z to x.y.z-dev (`skbio.__version__` should be the only thing needing an update). Update ``CHANGELOG.md`` to include a new section for x.y.z-dev (there won't be any changes to note here yet).

2. Close the release milestone on the GitHub issue tracker if there was one.

3. Send an email to the skbio developers list and anyone else who might be interested (e.g., lab mailing lists). You might include links to the GitHub release page.

4. Tweet about the release from `@scikit-bio`, including a link to the GitHub release page (for example, https://github.com/biocore/scikit-bio/releases/tag/x.y.z). Post a similar message to [scikit-bio's Gitter](https://gitter.im/biocore/scikit-bio).

5. :beers:
Please complete the following checklist:

* [ ] I have read the guidelines in [CONTRIBUTING.md](https://github.com/biocore/scikit-bio/blob/master/CONTRIBUTING.md).

* [ ] I have documented all public-facing changes in [CHANGELOG.md](https://github.com/biocore/scikit-bio/blob/master/CHANGELOG.md).

* [ ] **This pull request includes code, documentation, or other content derived from external source(s).** If this is the case, ensure the external source's license is compatible with scikit-bio's [license](https://github.com/biocore/scikit-bio/blob/master/COPYING.txt). Include the license in the `licenses` directory and add a comment in the code giving proper attribution. Ensure any other requirements set forth by the license and/or author are satisfied. **It is your responsibility to disclose code, documentation, or other content derived from external source(s).** If you have questions about whether something can be included in the project or how to give proper attribution, include those questions in your pull request and a reviewer will assist you.

* [ ] **This pull request does not include code, documentation, or other content derived from external source(s).**

**Note:** [REVIEWING.md](https://github.com/biocore/scikit-bio/blob/master/REVIEWING.md) may also be helpful to see some of the things code reviewers will be verifying when reviewing your pull request.
Files in this directory are the QIIME 1.9.1 "tiny test" files. For detail on how these were created, see skbio/diversity/beta/tests/data/qiime-191-tt/README.md.
Files in this directory are the QIIME 1.9.1 "tiny test" files. These data were developed by @gregcaporaso, who gave permission to reproduce them in scikit-bio.

If you have a [QIIME 1.9.1 base installation](http://install.qiime.org), the raw input files in this directory can be obtained by running:

```bash
python -c "from qiime.test import write_test_data; write_test_data('.')"
biom convert -i biom --to-tsv -o otu-table.tsv
```

After converting to tsv, the following OTUs are removed because they are not present in the tree (they're not 16S sequences, so can't be aligned with PyNAST): ``None1``, ``None10``, ``None6``, and ``None2``. The ``not16S.1`` sample is also removed because, after removing those OTUs, it has a total count of 0. This boundary case is tested directly in the ``unifrac_*`` and ``faith_pd`` tests.

Then, in the python interpreter, we midpoint root the tree (since this is a QIIME 1.9.1 installation, this step is performed with scikit-bio 0.2.3):

```python
from skbio import TreeNode
t = TreeNode.read('./tree')
t = t.root_at_midpoint()
t.write('tree', format='newick')
```

The output files (alpha diversity values and beta diversity distance matrices) can then be obtained by running:

```bash
alpha_diversity.py -i biom -t tree -m PD_whole_tree -o pd.txt
beta_diversity.py -m weighted_unifrac,unweighted_unifrac,weighted_normalized_unifrac -i biom -t tree -o o
```
Caporaso lab meeting presentation
=================================

A presentation of scikit-bio, its goals, and a few live demos.

This was presented at the Caporaso lab meeting on 05/13/2014. It was tested
against pre-0.1.0 scikit-bio so may not work with newer versions of scikit-bio,
due to inevitable API changes that will happen. Tested using IPython Notebook
2.0.0 in slideshow mode.
scikit-bio documentation
========================

This guide contains instructions for building the scikit-bio documentation, as
well as guidelines for contributing to the documentation.

**Note:** If you're only interested in viewing the scikit-bio documentation,
visit [scikit-bio.org](http://scikit-bio.org).

Building the documentation
--------------------------

To build the documentation, you'll need a scikit-bio development environment
set up. See [CONTRIBUTING.md](../CONTRIBUTING.md) for instructions. In
addition, you will also need to install Sphinx and the theme for the
documentation, you can do that with:

    pip install 'Sphinx<=3.0' sphinx-bootstrap-theme

**Important:** The documentation will be built for whatever version of
scikit-bio is *currently installed* on your system (i.e., the version imported
by ```import skbio```). This may not match the code located in this repository.
You will need to either install this version of scikit-bio somewhere (e.g., in
a virtualenv) or point your ```PYTHONPATH``` environment variable to this code,
*before* building the documentation.

To build the documentation, assuming you are at the top-level scikit-bio
directory:

    make -C doc clean html

The built HTML documentation will be at ```doc/build/html/index.html```.

Contributing to the documentation
---------------------------------

If you would like to contribute to the documentation, whether by adding
something entirely new or by modifying existing documentation, please first
review our [scikit-bio contribution guide](../CONTRIBUTING.md).

Before submitting your changes, ensure that the documentation builds without
errors or warnings.

### Documentation guidelines

Most of scikit-bio's API documentation is automatically generated from
[docstrings](http://legacy.python.org/dev/peps/pep-0257/#what-is-a-docstring).
The advantage to this approach is that users can access the documentation in an
interactive Python session or from our website as HTML. Other output formats
are also possible, such as PDF.

scikit-bio docstrings follow the [numpydoc conventions](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt).
This ensures that the docstrings are easily readable both from the interpreter
and HTML, PDF, etc. Please read the numpydoc guidelines before continuing.

### Documenting a module in scikit-bio

In addition to following the numpydoc conventions for docstrings, we have a few
more conventions that will ensure your documentation is correctly built and
linked within our website, and that it maintains consistency with the rest of
the scikit-bio docs.

The easiest way to get started with documenting your code is to look at the
docstrings in existing scikit-bio modules. A couple of modules to start with
are ```skbio.sequence``` and ```skbio.stats.distance```. Go ahead and look
through those now. We've structured our docs in a similar way to
[SciPy's documentation](http://docs.scipy.org/doc/scipy/reference/), so that
may be another good place to look for examples.

We'll take a top-down approach by discussing how to document a new module that
you'd like to add to scikit-bio (let's call it ```skbio/example.py```).

#### Module docstring

The first thing you'll need to add is a docstring for the module. The docstring
must start at the first line of the file. It should start with a title for the
module:

    """
    Documentation examples (:mod:`skbio.example`)
    =============================================

It is important to include the ```:mod:``` Sphinx directive in the title, as
this title will be included in the table of contents. Also make sure that the
title underline is the same length as the title.

We also need to include another Sphinx directive below this:

    .. currentmodule:: skbio.example

This directive tells Sphinx that other classes, functions, etc. that we will
reference are located in the ```skbio.example``` module.

Next, include a more detailed description of the module. For example:

    This module consists of several example classes and functions to illustrate
    the scikit-bio documentation system.

Following that, list any classes, functions, and exceptions that you'd like
documentation generated for. Note that you do *not* need to include every
single class, function, or exception that is defined in the module. Also, you
do not need to list class methods, as those will be automatically included in
the generated class documentation. Only include objects that should be exposed
as part of the public API.

For example:

    Classes
    -------

    .. autosummary::
       :toctree: generated/

       ExampleClass1
       ExampleClass2

    Functions
    ---------

    .. autosummary::
       :toctree: generated/

       example_function1
       example_function2

    Exceptions
    ----------

    .. autosummary::
       :toctree: generated/

       ExampleError

The ```autosummary``` directives are important as they generate RST files in
the ```generated/``` directory for each object. A single-line summary and link
to each object is inserted into the page for you.

After listing public module members, we encourage a usage example section
showing how to use some of the module's functionality. Examples should be
written in [doctest](http://docs.python.org/3/library/doctest.html) format so
that they can be automatically tested (e.g., using ```make test``` or
```python -m skbio.test```).

    Examples
    --------

    Run the ``example_function1`` function:

    >>> from skbio.example import example_function1
    >>> example_function1("hello", "world")
    hello world!

You can also embed the plots that an example generates into the built
documentation with the ```.. plot::``` directive. For example:

    .. plot::

       >>> import pandas as pd
       >>> df = pd.DataFrame({'col1': [1, 2, 3, 4], 'col2': [10, 11, 12, 13]})
       >>> fig = df.boxplot()

This will include the plot, a link to the source code used to generate the
plot, and links to different image formats (e.g., PNG and PDF) so that users
can easily download the plot.

You're now ready to document the members of your module.

#### Documenting module members

When documenting the members of a module (e.g., classes, methods, attributes,
functions, and exceptions), follow the numpydoc conventions. In addition to
these conventions, there are a few things to keep in mind:

- When documenting a class, only public methods and attributes are included in
  the built documentation. If a method or attribute starts with an
  underscore, it is assumed to be private.

- When documenting a class, include the ```Parameters``` section in the class
  docstring, instead of in the ```__init__``` docstring. While numpydoc
  technically supports either form, ```__init__``` is not included in the list
  of methods by default and thus should have its documentation included in the
  class docstring.

#### Including the module in the docs

Until now, we've only been editing docstrings, which are attached to Python
code. The final step is to hook up this new module's docstrings to the
documentation build system:

1. Make sure you're within the ```scikit-bio/doc``` directory.
2. Create a new file with the same name as your module under the ```source```
   directory. Do not include ```skbio``` as part of the name, and use
   ```.rst``` as the suffix. For example, ```source/example.rst```.
3. Add the following line to ```source/example.rst``` to have your module's
   docstring pulled into the document:

    ```
    .. automodule:: skbio.example
    ```

4. Add the following line to ```source/index.rst``` to add the new page to the
   top-level table of contents:

    ```
    example
    ```

That's it! You can now try building the documentation, which should include the
documentation for your new module!

### Documenting a subpackage in scikit-bio

The process of documenting a subpackage is very similar to documenting a module
in scikit-bio. The only difference is that the module docstring goes in the
subpackage's ```__init__.py```.

### Troubleshooting

If things aren't working correctly, try running ```make clean``` and then
rebuild the docs. If things still aren't working, try building the docs
*without* your changes, and see if there are any Sphinx errors or warnings.
Make note of these, and then see what new errors or warnings are generated when
you add your changes again.

.. image:: http://scikit-bio.org/assets/logo.svg
   :target: http://scikit-bio.org
   :alt: scikit-bio logo

|Build Status| |Coverage Status| |ASV Benchmarks| |Gitter Badge| |Depsy Badge| |Anaconda Build Platforms| |Anaconda Build Version| |License| |Downloads| |Install|

scikit-bio is an open-source, BSD-licensed Python 3 package providing data structures, algorithms and educational resources for bioinformatics.

To view scikit-bio's documentation, visit `scikit-bio.org
<http://scikit-bio.org>`__.

**Note:** scikit-bio is no longer compatible with Python 2. scikit-bio is compatible with Python 3.6 and later.

scikit-bio is currently in beta. We are very actively developing it, and **backward-incompatible interface changes can and will arise**. To avoid these types of changes being a surprise to our users, our public APIs are decorated to make it clear to users when an API can be relied upon (stable) and when it may be subject to change (experimental). See the `API stability docs <https://github.com/biocore/scikit-bio/blob/master/doc/source/user/api_stability.rst>`_ for more details, including what we mean by *stable* and *experimental* in this context.

Installing
----------

The recommended way to install scikit-bio is via the ``conda`` package manager available in `Anaconda <http://continuum.io/downloads>`_ or `miniconda <http://conda.pydata.org/miniconda.html>`_.

To install the latest release of scikit-bio::

    conda install -c conda-forge scikit-bio

Alternatively, you can install scikit-bio using ``pip``::

    pip install scikit-bio

You can verify your installation by running the scikit-bio unit tests::

    python -m skbio.test

For users of Debian, ``skbio`` is in the Debian software distribution and may
be installed using::

    sudo apt-get install python3-skbio python-skbio-doc


Getting help
------------

To get help with scikit-bio, you should use the `skbio <http://stackoverflow.com/questions/tagged/skbio>`_ tag on StackOverflow (SO). Before posting a question, check out SO's guide on how to `ask a question <http://stackoverflow.com/questions/how-to-ask>`_. The scikit-bio developers regularly monitor the ``skbio`` SO tag.

Projects using scikit-bio
-------------------------

Some of the projects that we know of that are using scikit-bio are:

- `QIIME <http://qiime.org/>`__
- `Emperor <http://biocore.github.io/emperor/>`__
- `An Introduction to Applied
  Bioinformatics <http://readIAB.org>`__
- `tax2tree <https://github.com/biocore/tax2tree>`__
- `Qiita <http://qiita.microbio.me>`__
- `ghost-tree <https://github.com/JTFouquier/ghost-tree>`__
- `Platypus-Conquistador <https://github.com/biocore/Platypus-Conquistador>`__

If you're using scikit-bio in your own projects, feel free to issue a pull request to add them to this list.

scikit-bio development
----------------------

If you're interested in getting involved in scikit-bio development, see `CONTRIBUTING.md <https://github.com/biocore/scikit-bio/blob/master/CONTRIBUTING.md>`__.

See the list of `scikit-bio's contributors
<https://github.com/biocore/scikit-bio/graphs/contributors>`__.

Licensing
---------

scikit-bio is available under the new BSD license. See
`COPYING.txt <https://github.com/biocore/scikit-bio/blob/master/COPYING.txt>`__ for scikit-bio's license, and the
`licenses directory <https://github.com/biocore/scikit-bio/tree/master/licenses>`_ for the licenses of third-party software that is
(either partially or entirely) distributed with scikit-bio.

The pre-history of scikit-bio
-----------------------------

scikit-bio began from code derived from `PyCogent
<http://www.pycogent.org>`__ and `QIIME <http://www.qiime.org>`__, and
the contributors and/or copyright holders have agreed to make the code
they wrote for PyCogent and/or QIIME available under the BSD
license. The contributors to PyCogent and/or QIIME modules that have
been ported to scikit-bio are: Rob Knight (`@rob-knight
<https://github.com/rob-knight>`__), Gavin Huttley (`@gavin-huttley
<https://github.com/gavin-huttley>`__), Daniel McDonald (`@wasade
<https://github.com/wasade>`__), Micah Hamady, Antonio Gonzalez
(`@antgonza <https://github.com/antgonza>`__), Sandra Smit, Greg
Caporaso (`@gregcaporaso <https://github.com/gregcaporaso>`__), Jai
Ram Rideout (`@jairideout <https://github.com/jairideout>`__),
Cathy Lozupone (`@clozupone <https://github.com/clozupone>`__), Mike Robeson
(`@mikerobeson <https://github.com/mikerobeson>`__), Marcin Cieslik,
Peter Maxwell, Jeremy Widmann, Zongzhi Liu, Michael Dwan, Logan Knecht
(`@loganknecht <https://github.com/loganknecht>`__), Andrew Cochran,
Jose Carlos Clemente (`@cleme <https://github.com/cleme>`__), Damien
Coy, Levi McCracken, Andrew Butterfield, Will Van Treuren (`@wdwvt1
<https://github.com/wdwvt1>`__), Justin Kuczynski (`@justin212k
<https://github.com/justin212k>`__), Jose Antonio Navas Molina
(`@josenavas <https://github.com/josenavas>`__), Matthew Wakefield
(`@genomematt <https://github.com/genomematt>`__) and Jens Reeder
(`@jensreeder <https://github.com/jensreeder>`__).

Logo
----

scikit-bio's logo was created by `Alina Prassas <http://cargocollective.com/alinaprassas>`_.

.. |Build Status| image:: https://travis-ci.org/biocore/scikit-bio.svg?branch=master
   :target: https://travis-ci.org/biocore/scikit-bio
.. |Coverage Status| image:: https://coveralls.io/repos/biocore/scikit-bio/badge.png
   :target: https://coveralls.io/r/biocore/scikit-bio
.. |ASV Benchmarks| image:: http://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat
   :target: https://s3-us-west-2.amazonaws.com/scikit-bio.org/benchmarks/master/index.html
.. |Gitter Badge| image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/biocore/scikit-bio
   :target: https://gitter.im/biocore/scikit-bio?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
.. |Depsy Badge| image:: http://depsy.org/api/package/pypi/scikit-bio/badge.svg
   :target: http://depsy.org/package/python/scikit-bio
.. |Anaconda Build Platforms| image:: https://anaconda.org/conda-forge/scikit-bio/badges/platforms.svg
   :target: https://anaconda.org/conda-forge/scikit-bio
.. |Anaconda Build Version| image:: https://anaconda.org/conda-forge/scikit-bio/badges/version.svg
   :target: https://anaconda.org/conda-forge/scikit-bio
.. |License| image:: https://anaconda.org/conda-forge/scikit-bio/badges/license.svg
   :target: https://anaconda.org/conda-forge/scikit-bio
.. |Downloads| image:: https://anaconda.org/conda-forge/scikit-bio/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/scikit-bio
.. |Install| image:: https://anaconda.org/conda-forge/scikit-bio/badges/installer/conda.svg
   :target: https://conda.anaconda.org/conda-forge
.. automodule:: skbio.workflow
.. automodule:: skbio.util
.. automodule:: skbio.tree
.. automodule:: skbio.metadata
.. automodule:: skbio.diversity
.. automodule:: skbio.alignment
scikit-bio
==========

scikit-bio (canonically pronounced *sigh-kit-buy-oh*) is a library for working
with biological data in Python 3. scikit-bio is open source, BSD-licensed
software that is currently under active development.

API Reference
-------------

.. toctree::
   :maxdepth: 1

   io
   sequence
   alignment
   tree
   workflow
   diversity
   stats
   metadata
   util

User Documentation
------------------

The user documentation contains high-level information for users of scikit-bio.

.. toctree::
   :maxdepth: 1

   user/api_stability

Developer Documentation
-----------------------

The developer documentation contains information for how to contribute
to scikit-bio.

.. toctree::
   :maxdepth: 1

   development/coding_guidelines
   development/new_module
.. automodule:: skbio.io
.. automodule:: skbio.stats
.. automodule:: skbio.sequence
.. automodule:: {{ fullname }}

:orphan:

{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. automethod:: {{ objname }}

:orphan:

{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoattribute:: {{ objname }}

{# This template was modified from autosummaries default format #}
{{ fullname | escape | underline}}

{# We need a list of the built-ins that we implemented, not the default ones #}
{% set built_in_methods = []  %}
{% for item in all_methods %}
   {% if (item not in ['__class__',
                        '__delattr__',
                        '__getattribute__',
                        '__init__',
                        '__dir__',
                        '__format__',
                        '__new__',
                        '__reduce__',
                        '__reduce_ex__',
                        '__repr__',
                        '__setattr__',
                        '__sizeof__',
                        '__subclasshook__',
                        '__init_subclass__',
                        '__class_getitem__'] and item.startswith('__')) %}
      {{ built_in_methods.append(item) or '' }}
   {% endif %}
{% endfor %}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}

   {% if built_in_methods %}
   .. rubric:: Built-ins

   .. autosummary::
      :toctree:
   {% for item in built_in_methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree:
   {% for item in methods %}
      {% if item != '__init__' %}
      ~{{ name }}.{{ item }}
      {% endif %}
   {%- endfor %}
   {% endif %}
Coding guidelines
=================

As project size increases, consistency of the code base and documentation becomes more important. We therefore provide guidelines for code and documentation that is contributed to scikit-bio. Our goal is to create a consistent code base where:

* it is easy to find relevant functionality (and to determine when functionality that you're looking for doesn't exist),
* you can trust that the code that you're working with is sufficiently tested, and
* names and interfaces are intuitive.

**As scikit-bio is in beta, our coding guidelines are presented here as a working draft. These guidelines are requirements for all code submitted to scikit-bio, but at this stage the guidelines themselves are malleable. If you disagree with something, or have a suggestion for something new to include, you should** `create an issue`_ **to initiate a discussion.**

.. _`create an issue`: https://github.com/biocore/scikit-bio/issues

What are the naming conventions? and How should I format my code?
-----------------------------------------------------------------

We adhere to the `PEP 8`_ python coding guidelines for code and documentation standards. Before submitting any code to scikit-bio, you should read these carefully and apply the guidelines in your code.

.. _`PEP 8`: http://legacy.python.org/dev/peps/pep-0008/


What should I call my variables?
--------------------------------

- *Choose the name that people will most likely guess.* Make it descriptive, but not too long: ``curr_record`` is better than ``c``, or ``curr``, or ``current_genbank_record_from_database``.

- *Good names are hard to find.* Don't be afraid to change names except when they are part of interfaces that other people are also using. It may take some time working with the code to come up with reasonable names for everything: if you have unit tests, it's easy to change them, especially with global search and replace.

- *Use singular names for individual things, plural names for collections.* For example, you'd expect ``self.name`` to hold something like a single string, but ``self.names`` to hold something that you could loop through like a list or dictionary. Sometimes the decision can be tricky: is ``self.index`` an integer holding a positon, or a dictionary holding records keyed by name for easy lookup? If you find yourself wondering these things, the name should probably be changed to avoid the problem: try ``self.position`` or ``self.look_up``.

- *Don't make the type part of the name.* You might want to change the implementation later. Use ``Records`` rather than ``RecordDict`` or ``RecordList``, etc. Don't use Hungarian Notation either (i.e. where you prefix the name with the type).

- *Make the name as precise as possible.* If the variable is the path of the input file, call it ``input_fp``, not ``input`` or ``file`` (which you shouldn't use anyway, since they're keywords), and not ``infile`` (because that looks like it should be a file object, not just its name).

- *Use* ``result`` *to store the value that will be returned from a method or function.* Use ``data`` for input in cases where the function or method acts on arbitrary data (e.g. sequence data, or a list of numbers, etc.) unless a more descriptive name is appropriate.

- *One-letter variable names should only occur in math functions or as loop iterators with limited scope.* Limited scope covers things like ``for k in keys: print k``, where ``k`` survives only a line or two. Loop iterators should refer to the variable that they're looping through: ``for k in keys, i in items``, or ``for key in keys, item in items``. If the loop is long or there are several 1-letter variables active in the same scope, rename them.

- *Limit your use of abbreviations.* A few well-known abbreviations are OK, but you don't want to come back to your code in 6 months and have to figure out what ``sptxck2`` is. It's worth it to spend the extra time typing ``species_taxon_check_2``, but that's still a horrible name: what's check number 1? Far better to go with something like ``taxon_is_species_rank`` that needs no explanation, especially if the variable is only used once or twice.

Acceptable abbreviations
^^^^^^^^^^^^^^^^^^^^^^^^

The following list of abbreviations can be considered well-known and used with impunity within mixed name variables, but some should not be used by themselves as they would conflict with common functions, python built-in's, or raise an exception. Do not use the following by themselves as variable names: ``dir``,  ``exp`` (a common ``math`` module function), ``in``, ``max``, and ``min``. They can, however, be used as part of a name, eg ``matrix_exp``.

+--------------------+--------------+
|        Full        |  Abbreviated |
+====================+==============+
|          alignment |          aln |
+--------------------+--------------+
|           archaeal |         arch |
+--------------------+--------------+
|          auxiliary |          aux |
+--------------------+--------------+
|          bacterial |         bact |
+--------------------+--------------+
|           citation |         cite |
+--------------------+--------------+
|            current |         curr |
+--------------------+--------------+
|           database |           db |
+--------------------+--------------+
|         dictionary |         dict |
+--------------------+--------------+
|          directory |          dir |
+--------------------+--------------+
|    distance matrix |           dm |
+--------------------+--------------+
|        end of file |          eof |
+--------------------+--------------+
|         eukaryotic |          euk |
+--------------------+--------------+
|          filepath  |           fp |
+--------------------+--------------+
|          frequency |         freq |
+--------------------+--------------+
|           expected |          exp |
+--------------------+--------------+
|              index |          idx |
+--------------------+--------------+
|              input |           in |
+--------------------+--------------+
|            maximum |          max |
+--------------------+--------------+
|            minimum |          min |
+--------------------+--------------+
|      mitochondrial |           mt |
+--------------------+--------------+
|             number |          num |
+--------------------+--------------+
|        observation |          obs |
+--------------------+--------------+
|           observed |          obs |
+--------------------+--------------+
|           original |         orig |
+--------------------+--------------+
|             output |          out |
+--------------------+--------------+
|          parameter |        param |
+--------------------+--------------+
|          phylogeny |        phylo |
+--------------------+--------------+
|           previous |         prev |
+--------------------+--------------+
|        probability |         prob |
+--------------------+--------------+
|            protein |         prot |
+--------------------+--------------+
|             record |          rec |
+--------------------+--------------+
|          reference |          ref |
+--------------------+--------------+
|           sequence |          seq |
+--------------------+--------------+
| standard deviation |        stdev |
+--------------------+--------------+
|         statistics |        stats |
+--------------------+--------------+
|             string |          str |
+--------------------+--------------+
|          structure |       struct |
+--------------------+--------------+
|          temporary |         temp |
+--------------------+--------------+
|               taxa |          tax |
+--------------------+--------------+
|              taxon |          tax |
+--------------------+--------------+
|          taxonomic |          tax |
+--------------------+--------------+
|           taxonomy |          tax |
+--------------------+--------------+
|           variance |          var |
+--------------------+--------------+

How do I organize my modules (source files)?
--------------------------------------------

- *Have a docstring with a description of the module's functions*. If the description is long, the first line should be a short summary that makes sense on its own, separated from the rest by a newline.

- *All code, including import statements, should follow the docstring.* Otherwise, the docstring will not be recognized by the interpreter, and you will not have access to it in interactive sessions (i.e. through ``obj.__doc__``) or when generating documentation with automated tools.

- *Import built-in modules first, followed by third-party modules, followed by any changes to the path and your own modules.* Especially, additions to the path and names of your modules are likely to change rapidly: keeping them in one place makes them easier to find.

- *Don't use* ``from module import *``, *instead use* ``from module import Name, Name2, Name3...`` *or possibly* ``import module``. This makes it *much* easier to see name collisions and to replace implementations.

- If you are importing `NumPy`_, `Matplotlib`_, or another package that encourages a standard style for their import statements use them as needed for example:

::

    import numpy as np
    import numpy.testing as npt
    import pandas as pd

    from matplotlib import pyplot as plt

.. _`NumPy`: http://www.numpy.org/
.. _`Matplotlib`: http://matplotlib.org/

Example of module structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The structure of your module should be similar to the example below. scikit-bio uses the `NumPy doc`_ standard for documentation. Our `doc/README.md`_ explains how to write your docstrings using the `NumPy doc`_ standards for scikit-bio:

.. _`NumPy doc`: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
.. _`doc/README.md`: https://github.com/biocore/scikit-bio/blob/master/doc/README.md

.. code-block:: python

    r"""
    Numbers (:mod:`skbio.numbers`)
    ==============================

    .. currentmodule:: skbio.numbers

    Numbers holds a sequence of numbers, and defines several statistical
    operations (mean, stdev, etc.) FrequencyDistribution holds a mapping from
    items (not necessarily numbers) to counts, and defines operations such as
    Shannon entropy and frequency normalization.


    Classes
    -------

    .. autosummary::
       :toctree: generated/

       Numbers

    """

    # ----------------------------------------------------------------------------
    # Copyright (c) 2013--, scikit-bio development team.
    #
    # Distributed under the terms of the Modified BSD License.
    #
    # The full license is in the file COPYING.txt, distributed with this software.
    # ----------------------------------------------------------------------------

    from random import choice, random

    import numpy as np
    from utils import indices


    class Numbers(list):
        pass


    class FrequencyDistribution(dict):
        pass


How should I write comments?
----------------------------

- *Always update the comments when the code changes.* Incorrect comments are far worse than no comments, since they are actively misleading.

- *Comments should say more than the code itself.* Examine your comments carefully: they may indicate that you'd be better off rewriting your code (especially if *renaming your variables* would allow you to get rid of the comment.) In particular, don't scatter magic numbers and other constants that have to be explained through your code. It's far better to use variables whose names are self-documenting, especially if you use the same constant more than once. Also, think about making constants into class or instance data, since it's all too common for 'constants' to need to change or to be needed in several methods.

    +-------+------------------------------------------------------------+
    | Wrong |       ``win_size -= 20        # decrement win_size by 20`` |
    +-------+------------------------------------------------------------+
    |    OK | ``win_size -= 20        # leave space for the scroll bar`` |
    +-------+------------------------------------------------------------+
    | Right |                             ``self._scroll_bar_size = 20`` |
    +-------+------------------------------------------------------------+
    |       |                      ``win_size -= self._scroll_bar_size`` |
    +-------+------------------------------------------------------------+


- *Use comments starting with #, not strings, inside blocks of code.*
- *Start each method, class and function with a docstring using triple double quotes (""").* Make sure the docstring follows the `NumPy doc`_ standard.

- *Always update the docstring when the code changes.* Like outdated comments, outdated docstrings can waste a lot of time. "Correct examples are priceless, but incorrect examples are worse than worthless." `Jim Fulton`_.

.. _`Jim Fulton`: http://www.python.org/pycon/dc2004/papers/4/PyCon2004DocTestUnit.pdf

How should I test my code ?
---------------------------

There are several different approaches for testing code in python: ``nose``, ``unittest`` and ``numpy.testing``. Their purpose is the same, to check that execution of code given some input produces a specified output. The cases to which the approaches lend themselves are different.

Whatever approach is employed, the general principle is every line of code should be tested. It is critical that your code be fully tested before you draw conclusions from results it produces. For scientific work, bugs don't just mean unhappy users who you'll never actually meet: **they may mean retracted publications**.

Tests are an opportunity to invent the interface(s) you want. Write the test for a method before you write the method: often, this helps you figure out what you would want to call it and what parameters it should take. It's OK to write the tests a few methods at a time, and to change them as your ideas about the interface change. However, you shouldn't change them once you've told other people what the interface is. In the spirit of this, your tests should also import the functionality that they test from the shortest alias possible. This way any change to the API will cause your tests to break, and rightly so!

Never treat prototypes as production code. It's fine to write prototype code without tests to try things out, but when you've figured out the algorithm and interfaces you must rewrite it *with tests* to consider it finished. Often, this helps you decide what interfaces and functionality you actually need and what you can get rid of.

"Code a little test a little". For production code, write a couple of tests, then a couple of methods, then a couple more tests, then a couple more methods, then maybe change some of the names or generalize some of the functionality. If you have a huge amount of code where all you have to do is write the tests', you're probably closer to 30% done than 90%. Testing vastly reduces the time spent debugging, since whatever went wrong has to be in the code you wrote since the last test suite. And remember to use python's interactive interpreter for quick checks of syntax and ideas.

Run the test suite when you change `anything`. Even if a change seems trivial, it will only take a couple of seconds to run the tests and then you'll be sure. This can eliminate long and frustrating debugging sessions where the change turned out to have been made long ago, but didn't seem significant at the time. **Note that tests are executed using Travis CI**, see `this document's section`_ for further discussion.

.. _`this document's section`: https://github.com/biocore/scikit-bio/blob/master/CONTRIBUTING.md#testing-guidelines

Some pointers
^^^^^^^^^^^^^

- *Use the* ``unittest`` *or the* ``nose`` *framework with tests in a separate file for each module.* Name the test file ``test_module_name.py`` and include it inside the tests folder of the module. Keeping the tests separate from the code reduces the temptation to change the tests when the code doesn't work, and makes it easy to verify that a completely new implementation presents the same interface (behaves the same) as the old.

- *Always include an* ``__init__.py`` *file in your tests directory*. This is required for the module to be included when the package is built and installed via ``setup.py``.

- *Always import from a minimally deep API target*. That means you would use ``from skbio import DistanceMatrix`` instead of ``from skbio.stats.distance import DistanceMatrix``. This allows us prevent most cases of accidental regression in our API.

- *Use* ``numpy.testing`` *if you are doing anything with floating point numbers, arrays or permutations* (use ``numpy.testing.assert_almost_equal``). Do *not* try to compare floating point numbers using ``assertEqual`` if you value your sanity.

- *Test the interface of each class in your code by defining at least one* ``TestCase`` *with the name* ``ClassNameTests``. This should contain tests for everything in the public interface.

- *If the class is complicated, you may want to define additional tests with names* ``ClassNameTests_test_type``. These might subclass ``ClassNameTests`` in order to share ``setUp`` methods, etc.

- *Tests of private methods should be in a separate* ``TestCase`` *called* ``ClassNameTests_private``. Private methods may change if you change the implementation. It is not required that test cases for private methods pass when you change things (that's why they're private, after all), though it is often useful to have these tests for debugging.

- *Test `all` the methods in your class.* You should assume that any method you haven't tested has bugs. The convention for naming tests is ``test_method_name``. Any leading and trailing underscores on the method name can be ignored for the purposes of the test; however, *all tests must start with the literal substring* ``test`` *for* ``unittest`` and ``nose`` *to find them.* If the method is particularly complex, or has several discretely different cases you need to check, use ``test_method_name_suffix``, e.g. ``test_init_empty``, ``test_init_single``, ``test_init_wrong_type``, etc. for testing ``__init__``.

- *Docstrings for testing methods should be considered optional*, instead the description of what the method does should be included in the name itself, therefore the name should be descriptive enough such that when running the tests in verbose mode you can immediately see the file and test method that's failing.

.. code-block:: none

    $ python -c "import skbio; skbio.test(verbose=True)"
    skbio.maths.diversity.alpha.tests.test_ace.test_ace ... ok
    test_berger_parker_d (skbio.maths.diversity.alpha.tests.test_base.BaseTests) ... ok

    ----------------------------------------------------------------------
    Ran 2 tests in 0.1234s

    OK

- *Module-level functions should be tested in their own* ``TestCase``\ *, called* ``modulenameTests``. Even if these functions are simple, it's important to check that they work as advertised.

- *It is much more important to test several small cases that you can check by hand than a single large case that requires a calculator.* Don't trust spreadsheets for numerical calculations -- use R instead!

- *Make sure you test all the edge cases: what happens when the input is None, or '', or 0, or negative?* What happens at values that cause a conditional to go one way or the other? Does incorrect input raise the right exceptions? Can your code accept subclasses or superclasses of the types it expects? What happens with very large input?

- *To test permutations, check that the original and shuffled version are different, but that the sorted original and sorted shuffled version are the same.* Make sure that you get *different* permutations on repeated runs and when starting from different points.

- *To test random choices, figure out how many of each choice you expect in a large sample (say, 1000 or a million) using the binomial distribution or its normal approximation.* Run the test several times and check that you're within, say, 3 standard deviations of the mean.

- All tests that depend on a random value should be seeded, for example if using NumPy, `numpy.random.seed(0)` should be used, in any other case the appropriate API should be used to create consistent outputs between runs. It is preferable that you do this for each test case instead of doing it in the `setUp` function/method (if any exists).

- Stochastic failures should occur less than 1/10,1000 times, otherwise you risk adding a significant amount of time to the total running time of the test suite.

Example of a ``nose`` test module structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # ----------------------------------------------------------------------------
    # Copyright (c) 2013--, scikit-bio development team.
    #
    # Distributed under the terms of the Modified BSD License.
    #
    # The full license is in the file COPYING.txt, distributed with this software.
    # ----------------------------------------------------------------------------

    import numpy as np
    from nose.tools import assert_almost_equal, assert_raises

    from skbio.math.diversity.alpha.ace import ace


    def test_ace():
        assert_almost_equal(ace(np.array([2, 0])), 1.0)
        assert_almost_equal(ace(np.array([12, 0, 9])), 2.0)
        assert_almost_equal(ace(np.array([12, 2, 8])), 3.0)
        assert_almost_equal(ace(np.array([12, 2, 1])), 4.0)
        assert_almost_equal(ace(np.array([12, 1, 2, 1])), 7.0)
        assert_almost_equal(ace(np.array([12, 3, 2, 1])), 4.6)
        assert_almost_equal(ace(np.array([12, 3, 6, 1, 10])), 5.62749672)

        # Just returns the number of OTUs when all are abundant.
        assert_almost_equal(ace(np.array([12, 12, 13, 14])), 4.0)

        # Border case: only singletons and 10-tons, no abundant OTUs.
        assert_almost_equal(ace([0, 1, 1, 0, 0, 10, 10, 1, 0, 0]), 9.35681818182)


    def test_ace_only_rare_singletons():
        with assert_raises(ValueError):
            ace([0, 0, 43, 0, 1, 0, 1, 42, 1, 43])


    if __name__ == '__main__':
        import nose
        nose.runmodule()

Git pointers
------------

Commit messages are a useful way to document the changes being made to a project, it additionally documents who is making these changes and when are these changes being made, all of which are relevant when tracing back problems.

Authoring a commit message
^^^^^^^^^^^^^^^^^^^^^^^^^^

The most important metadata in a commit message is (arguably) the author's name and the author's e-mail. GitHub uses this information to attribute your contributions to a project, see for example the `scikit-bio list of contributors`_.

.. _`scikit-bio list of contributors`: https://github.com/biocore/scikit-bio/graphs/contributors

Follow `this guide`_ to set up your system and **make sure the e-mail you use in this step is the same e-mail associated to your GitHub account**.

.. _`this guide`: http://git-scm.com/book/en/Getting-Started-First-Time-Git-Setup

After doing this you should see your name and e-mail when you run the following commands:

.. code-block:: none

    $ git config --global user.name
    Yoshiki Vázquez Baeza
    $ git config --global user.email
    yoshiki89@gmail.com

Writing a commit message
^^^^^^^^^^^^^^^^^^^^^^^^

In general the writing of a commit message should adhere to `NumPy's guidelines`_ which if followed correctly will help you structure your changes better i. e. bug fixes will be in a commit followed by a commit updating the test suite and with one last commit that update the documentation as needed.

GitHub provides a set of handy features that will link together a commit message to a ticket in the issue tracker, this is specially helpful because you can `close an issue automatically`_ when the change is merged into the main repository, this reduces the amount of work that has to be done making sure outdated issues are not open.

.. _`NumPy's guidelines`: http://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html#writing-the-commit-message
.. _`close an issue automatically`: https://help.github.com/articles/closing-issues-via-commit-messages
Adding a new module to skbio
############################

Each module needs an `__init__.py` file and a `tests` folder that also
contains an `__init__.py` file. For a module, a simple one may look
like this::

  r"""
  A module (:mod:`skbio.module`)
  ==============================

  .. currentmodule:: skbio.module

  Documentation for this module.
  """

  # ----------------------------------------------------------------------------
  # Copyright (c) 2013--, scikit-bio development team.
  #
  # Distributed under the terms of the Modified BSD License.
  #
  # The full license is in the file COPYING.txt, distributed with this software.
  # ----------------------------------------------------------------------------

  from skbio.util import TestRunner
  test = TestRunner(__file__).test

Usually, some functionality from the module will be made accessible by
importing it in `__init__.py`. It's convenient to use explicit
relative imports (`from .implementation import compute`), so that
functionality can be neatly separated in different files but the user
doesn't face a deeply nested package: `from skbio.module import
compute` instead of `from skbio.module.implementation import compute`.

Inside the tests folder, a simpler `__init__.py` works fine (it is
necessary so that all tests can be run after installation)::

  # ----------------------------------------------------------------------------
  # Copyright (c) 2013--, scikit-bio development team.
  #
  # Distributed under the terms of the Modified BSD License.
  #
  # The full license is in the file COPYING.txt, distributed with this software.
  # ----------------------------------------------------------------------------

Finally, remember to also follow the `documentation guidelines
<https://github.com/biocore/scikit-bio/blob/master/doc/README.md#documenting-a-module-in-scikit-bio>`_.
API Stability
=============

All public functionality in scikit-bio has a defined stability state.
These states inform users and developers to what extent they can rely on
different APIs in the package.

You can find out the stability state of public functionality by looking at its
docstring, which is formatted based on
`numpydoc <https://github.com/numpy/numpydoc>`_. This information will either
be in the *Extended Summary* section of the docstring, or in the case of
deprecation, this information will appear as a note following the *Short
Summary*.

The following diagram illustrates the API lifecycle in scikit-bio:

.. image:: assets/api-lifecycle.png
   :align: center

Definitions of the stability states and the information associated with each
follow.

Stable
------
Functionality defined as stable is part of scikit-bio's backward-
compatible API. Users can be confident that the API will not change without
first passing through the deprecated state, typically for at least two
release cycles. We make every effort to maintain the API of this code.

The docstrings of stable functionality will indicate the first scikit-bio
version where the functionality was considered stable.

Experimental
------------
Functionality defined as experimental is being considered for addition to
scikit-bio's stable API. Users are encouraged to use this code, but to be
aware that its API may change or be removed. Experimental functionality
will typically pass through the deprecated state before it is removed, but
in rare cases it may be removed directly (for example, if a serious
methodological flaw is discovered that makes the functionality
scientifically invalid).

The docstrings of experimental functionality will indicate the first
scikit-bio version where the functionality was considered experimental.

We aim to move functionality through the experimental phase quickly (for
example, two releases before moving to stable), but we don't make specific
promises about when experimental functionality will become stable. This
aligns with our philosophy that we don't make promises about experimental
APIs, only about stable APIs.

Deprecated
----------
Functionality defined as deprecated is targeted for removal from
scikit-bio. Users should transition away from using it.

The docstrings of deprecated functionality will indicate the first version
of scikit-bio where the functionality was deprecated, the version of
scikit-bio when the functionality will be removed, and the reason for
deprecation of the code (for example, because a function was determined to
be scientifically invalid, or because the API was adapted, and users should
be using a different version of the function).

Using deprecated functionality will raise a ``DeprecationWarning``. Since
Python 2.7, these types of warnings are **silenced by default**. When
developing a tool that uses scikit-bio, we recommend enabling the display of
deprecation warnings to be informed of upcoming API changes. For details on how
to display deprecation warnings, see `Python's deprecation warning docs
<https://docs.python.org/3/whatsnew/2.7.html#changes-to-the-handling-of-deprecation-warnings>`_.
