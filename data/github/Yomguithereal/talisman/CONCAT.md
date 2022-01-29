# Changelog

## 1.1.4

* Fixing `phonetics/nysiis`.

## 1.1.3

* Initial JOSS release.

## 1.1.2

* Fixing `phonetics/french/sonnex` (@tuzepoito).
* Fixing `keyers/name-sig`.

## 1.1.1

* Dropping `.babelrc` from npm package file to avoid bundling issues.

## 1.1.0

* Adding the `tokenizers/words/gersam` namespace.

## 1.0.0

* Dropping some namespaces & consolidating a stable release.

## 0.21.0

* Exposing Punkt Trainer's configuration options.

## 0.20.0

* Adding the `keyers/name-sig` namespace.
* Adding the `metrics/distance/lig` namespace.
* Adding the `phonetics/sound-d` namespace.
* Fixing `keyers/omission` & `keyers/skeleton` and better unit tests.

## 0.19.1

* Improving `tokenizers/sentences/naive`.
* Dropping `generatorics` in favor of [obliterator](https://github.com/Yomguithereal/obliterator).

## 0.19.0

* Adding the `clustering/record-linkage/leader` namespace.
* Adding the `keyers/name-power-set` namespace.
* Adding the `keyers/normalize` namespace.
* Adding the `tokenizers/fingerprint/name` namespace.
* Adding the `tokenizers/skipgrams` namespace.
* Fixing & optimizing `clustering/k-means`.
* Dropping the `helpers/random` namespace in favor of [pandemonium](https://github.com/Yomguithereal/pandemonium).

## 0.18.0

* Adding the `metrics/distance/guth` namespace.
* Fixing a bug related to Levenshtein distance prefix trimming.
* Fixing a bug related to `clustering/k-means`.

## 0.17.0

* Fixing `metrics/distance/jaro-winkler`.
* Improving `metrics/distance/overlap` performance.
* Dropping the `structures` namespace in favor of [mnemonist](https://github.com/Yomguithereal/mnemonist).

## 0.16.0

* Changing the way the fingerprint API.
* Providing index of item in some `clustering/record-linkage` callbacks.
* Adding `merge` option to `clustering/record-linkage/key-collision`.
* Adding the `keyers/fingerprint` namespace back.
* Moving `phonetics/omission` & `phonetics/skeleton` back to the `keyers` namespace.
* Improving `metrics/distance/levenshtein` performance.

## 0.15.0

* Adding the `hash/crc32` namespace.
* Adding the `hash/minhash` namespace.
* Adding the `helpers/random#createRandomIndex` function.
* Adding the `helpers/random#createChoice` function.
* Adding the `helpers/random#createDangerousButPerformantSample` function.
* Adding the `helpers/random#createSuffleInPlace` function.
* Adding the `clustering/record-linkage/blocking` namespace.
* Adding the `clustering/record-linkage/canopy` namespace.
* Adding the `distance/metrics/bag` namespace.
* Adding the `distance/metrics/lcs` namespace.
* Adding the `distance/metrics/length` namespace.
* Adding the `distance/metrics/minhash` namespace.
* Adding the `distance/metrics/mlipns` namespace.
* Adding the `distance/metrics/prefix` namespace.
* Adding the `distance/metrics/ratcliff-obershelp` namespace.
* Adding the `distance/metrics/sift4` namespace.
* Adding the `distance/metrics/smith-waterman` namespace.
* Adding the `distance/metrics/suffix` namespace.
* Adding the `phonetics/french/fonem` namespace.
* Adding the `regexp` namespace.
* Adding the `stemmers/french/carry` namespace.
* Adding the `stemmers/french/eda` namespace.
* Adding the `tokenizers/fingerprint` namespace.
* Adding the `asymmetric` option to `clustering/naive`.
* Adding the `minClusterSize` option to clusterers.
* Adding limited version of `metrics/distance/damerau-levenshtein`.
* Adding limited version of `metrics/distance/levenshtein`.
* Adding bitwise version of `metrics/distance/hamming`.
* Adding normalized version of `metrics/distance/hamming`.
* Moving `stats/ngrams` to `tokenizers/ngrams`.
* Moving `keyers/omission` to `phonetics/omission`.
* Moving `keyers/skeleton` to `phonetics/skeleton`.
* Moving similarity clusterers to the `clustering/record-linkage` namespace.
* Dropping the `stats/tfidf` namespace.
* Dropping the `keyers` namespace.

## 0.14.0

* Dropping the `phonetics/spanish/fonetico` namespace (should use [phonogram](https://github.com/Yomguithereal/phonogram) now).
* Improving `VPTree` performance by building the tree iteratively.
* Found a way to ease CommonJS by getting rid of pesky `.default`.

## 0.13.0

* Adding the `clustering/key-collision` namespace.
* Adding the `clustering/naive` namespace.
* Adding the `clustering/sorted-neighborhood` namespace.
* Adding the `clustering/vp-tree` namespace.
* Adding the `metrics/distance/identity` namespace.
* Adding the `metrics/distance/monge-elkan` namespace.
* Reversing `structures/bk-tree#search` arguments.

## 0.12.0

* Adding the `helpers/random` namespace.
* Adding the `inflectors/spanish/noun` namespace.
* Adding the `keyers/fingerprint` namespace.
* Adding the `keyers/ngram-fingerprint` namespace.
* Adding the `keyers/omission` namespace.
* Adding the `keyers/skeleton` namespace.
* Adding the `tag/averaged-perceptron` namespace.
* Adding the `parsers/brown` namespace.
* Adding the `parsers/conll` namespace.
* Adding the `phonetics/onca` namespace.
* Adding the `stemmers/spanish/unine` namespace.
* Adding the `structures/bk-tree` namespace.
* Adding the `structures/symspell` namespace.
* Adding the `structures/vp-tree` namespace.
* Adding the `sampler` options to `clustering/k-means`.
* Adding the `stats/descriptive#.quantile` function.
* Adding the `stats/descriptive#.median` function.
* Fixing a bug with `clustering/k-means` where k would be superior to the number of vectors.
* Fixing a bug with `clustering/k-means` `initialCentroids` options.
* Fixing a bug with `clustering/k-means` where a vector could end up in several clusters.
* Dropping the internal `regex/classes` namespace.
* Dropping the `hasher` option of the `ngrams` functions.
* Dropping the `set-functions` dependency.

## 0.11.0

* Improving `clustering/k-means` API.
* Improving `tokenizers/syllable/sonoripy` hierarchy definition.

## 0.10.0

* Adding the `metrics/distance/eudex` namespace.
* Adding the `phonetics/eudex` namespace.
* Adding the `tokenizers/syllables/sonoripy` namespace.
* Adding some import shortcuts for naives tokenizers.
* Improving the `tokenizers/syllables/legalipy` API.
* Improving the `tokenizers/sentences/naive` API.
* Fixing `tokenizers/syllables/legalipy` to correctly handle capitalized words.

## 0.9.0

* Adding the `phonetics/alpha-sis` namespace.
* Adding the `phonetics/fuzzy-soundex` namespace.
* Adding the `phonetics/phonex` namespace.
* Adding the `stemmers/uea-lite` namespace.
* Adding the `stats/inferential#sampleCovariance` function.
* Adding the `stats/inferential#sampleCorrelation` function.
* Moving the `metrics` namespace to `metrics/distance`.

## 0.8.0

* Adding the `stemmers/s-stemmer` namespace.
* Adding the `tokenizers/hyphenation/liang` namespace.
* Adding the `tokenizers/lines/naive` namespace.
* Adding the `tokenizers/paragraphs/naive` namespace.
* Adding the `tokenizers/syllables/legalipy` namespace.
* Adding the `tokenizers/tweets/casual` namespace.
* Adding the `stats/frequencies#updateFrequencies` function.
* Adding polymorphism to the `stats/frequencies#relative` function.

## 0.7.0

* Adding the `features/extraction/vectorizer` namespace.
* Adding the `phonetics/lein` namespace.
* Adding the `phonetics/roger-root` namespace.
* Adding the `phonetics/statcan` namespace.

## 0.6.0

* Adding the `stemmers/french/unine` namespace.

## 0.5.0

* Adding the `phonetics/german/phonem` namespace.
* Adding the `phonetics/spanish/fonetico` namespace.
* Adding the `stemmers/german/caumanns` namespace.

## 0.4.0

* Adding the `phonetics/french/phonetic` namespace.
* Adding the `phonetics/french/phonex` namespace.
* Adding the `phonetics/french/sonnex` namespace.
* Adding the `phonetics/french/soundex` namespace.
* Adding the `phonetics/french/soundex2` namespace.
* Adding the `stemmers/french/porter` namespace.
* Adding the `tokenizers/sentences/punkt` namespace.
* Moving the Cologne phonetic algorithm to the `phonetics/german/cologne` namespace.

## 0.3.0

* Adding the `classification/naive-bayes` namespace.
* Adding the `classification/perceptron` namespace.
* Adding the `metrics/damerau-levenshtein` namespace.
* Adding the `stats/descriptive` namespace.
* Adding the `stats/inferential` namespace.
* Improving `metrics/levenshtein` performance.
* Fixing a bug with `stemmers/porter`.

## 0.2.0

* Adding `phonetics/daitch-mokotoff`.
* Adding `metrics/canberra`.
* Adding `metrics/chebyshev`.
* Adding `metrics/minkowski`.
* Adding the refined version of the Soundex.
* Changing `phonetics/doubleMetaphone` to `phonetics/double-metaphone`.

## 0.1.0

* Adding `stemmers/latin/schinke`.

## 0.0.1

* Initial release.
# Notes

## Roadmap

* Minhash (CRC32?) & Fuzzyhash (ssdeep) & simhash
* LSH Binning
* Method to get the number of expected calculations
* Distances
* Phonetics
* MVP Tree
* Higher order VP Tree
* Suffix Tree clustering
* NN-Descent
* NNCTPH
* Fast online K-NN
* Bitap
* KNN clustering
* LSH & MinHash & Rabin-Karp
* Inverted Index (Complex & Simple) / Array to store doc id + weight + number of positions + positions for memory efficiency (Integer Array)
* Write about the rationale behind the naive clustering composition methods

* Create keyers => with phonetics & fingerprint keyers

## Clustering

* Abstract a class for similarity clusterers to enable asynchronous work, possibility to abort computation etc.
* It should be possible to make some optimization to the naive clusterer (whose worst case would perform the same amount of computation) by comparing new elements to only one item of an already existing cluster.
* Method 3 should be possible to do without computing a graph but by holding a hashmap of items.
* Method to return the similarity graph and to get a cluster index rather. (canopy then blocking for instance).
* Method to return the time elapsed to compute.
* Clusterer should hold the number of comparisons made
* Clusterer should have chunk, async & emit progress events
* SNM clustering is quite efficient when using ngram fingerprinting & a really small window.
* The similarity graph must be undirected if you want the clusters to have a full diameter instead of radius somehow.
* Possible to create modes for the naive clusterer `normal`, `minLengthFirst`, `maxLengthFirst`, or even `full`?
* Check Java library about knng cluster extraction.

## UI

* Pre-processing.
* Inverted-Index of unique values to cut computations.
* Should shuffle the values before applying clustering.
* Difference between merge & harmonize.
* Cluster expansion through inverted index sorted by occurrences.
* Suggest methods based on size of the dataset.

## Recipes

* Ngram blocking or SNM.
* Double Metaphone key collision.
* Overlap coefficient on names.
* Minhash + ngrams for documents.
[![Build Status](https://travis-ci.org/Yomguithereal/talisman.svg)](https://travis-ci.org/Yomguithereal/talisman) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02405/status.svg)](https://doi.org/10.21105/joss.02405)

# Talisman

[Full documentation](https://yomguithereal.github.io/talisman/)

Talisman is a JavaScript library collecting algorithms, functions and various building blocks for [fuzzy matching](https://en.wikipedia.org/wiki/Approximate_string_matching), [information retrieval](https://en.wikipedia.org/wiki/Information_retrieval) and [natural language processing](https://en.wikipedia.org/wiki/Natural_language_processing).

## Installation

You can install **Talisman** through npm:

```bash
npm install talisman
```

## Documentation

The library's full documentation can be found [here](https://yomguithereal.github.io/talisman/).

## Bibliography

An extensive bibliography of the methods & functions implemented by the library can be found [here](./BIBLIOGRAPHY.md).

## Goals

* :package: **Modular**: the library is completely modular. This means that if you only need to compute a `levenshtein` distance, you will only load the relevant code.
* :bulb: **Straightforward & simple**: just want to compute a Jaccard index? No need to instantiate a class and use two methods to pass options and then finally succeed in getting the index. Just apply the `jaccard` function and get going.
* :dango: **Consistent API**: the library's API is fully consistent and one should not struggle to understand how to apply two different distance metrics.
* :postal_horn: **Functional**: except for cases where classes might be useful (clustering notably), *Talisman* only uses functions, consumes raw data and order functions' arguments to make partial application & currying etc. as easy as possible.
* :zap: **Performant**: the library should be as performant as possible for a high-level programming language library.
* :globe_with_meridians: **Cross-platform**: the library is cross-platform and can be used both with Node.js and in the browser.

## How to cite

**Talisman** has been published as a [paper](https://joss.theoj.org/papers/10.21105/joss.02405) on the [Journal Of Open Source Software (JOSS)](https://joss.theoj.org/).

## Contribution

Contributions are of course welcome :)

Be sure to lint & pass the unit tests before submitting your pull request.

```bash
# Cloning the repo
git clone git@github.com:Yomguithereal/talisman.git
cd talisman

# Installing the deps
npm install

# Running the tests
npm test

# Linting the code
npm run lint
```

## License

This project is available as open source under the terms of the [MIT License](./LICENSE.txt).
# Talisman Bibliography

[BibTex file](https://raw.githubusercontent.com/Yomguithereal/talisman/master/paper/algorithms.bib)

> Bartolini, I., Ciaccia, P., & Patella, M. (2002). String matching with metric trees using an approximate distance. International Symposium on String Processing and Information Retrieval, 271–283.

> Beider, A. (2008). Beider-Morse phonetic matching: An alternative to Soundex with fewer false hits. Avotaynu: The International Review of Jewish Genealogy, 24(2), 12.

> Black, P. E. (2004). Ratcliff/Obershelp pattern recognition. Dictionary of Algorithms and Data Structures, 17.

> Bouchard, G., Brard, P., & Lavoie, Y. (1981). FONEM: Un code de transcription phonétique pour la reconstitution automatique des familles saguenayennes. Population (French Edition), 1085–1103.

> Broder, A. Z. (1998). On the resemblance and containment of documents. Proceedings. Compression and Complexity of SEQUENCES 1997 (Cat. No.97TB100171), 21–29. https://doi.org/10.1109/SEQUEN.1997.666900

> Caumanns, J. (1999). A fast and simple stemming algorithm for German words.

> Damerau, F. J. (1964). A technique for computer detection and correction of spelling errors. Communications of the ACM, 7(3), 171–176.

> Dice, L. R. (1945). Measures of the amount of ecologic association between species. Ecology, 26(3), 297–302.

> Dong, W., Moses, C., & Li, K. (2011). Efficient k-nearest neighbor graph construction for generic similarity measures. Proceedings of the 20th International Conference on World Wide Web - WWW ’11, 577. https://doi.org/10.1145/1963405.1963487

> Guth, G. J. (1976). Surname spellings and computerized record linkage. Historical Methods Newsletter, 10(1), 10–19.

> Hamming, R. W. (1950). Error detecting and error correcting codes. The Bell System Technical Journal, 29(2), 147–160.

> Harman, D. (1991). How effective is suffixing? Journal of the American Society for Information Science, 42(1), 7–15.

> Holmes, D., & McCabe, C. M . (2002). Improving precision and recall for soundex retrieval. Proceedings. International Conference on Information Technology: Coding and Computing, 22–26.

> Hood, D. (2002). Caverphone: Phonetic matching algorithm. Technical Paper CTP060902, University of Otago, New Zealand.

> Hood, D. (2004). Caverphone revisited. Technical Paper CTP150804, 10.

> Jaccard, P. (1912). The distribution of the flora in the alpine zone. 1. New Phytologist, 11(2), 37–50. https://doi.org/10.1111/j.1469-8137.1912.tb05611.x

> Jaro, M. A. (1989). Advances in record-linkage methodology as applied to matching the 1985 census of Tampa, Florida. Journal of the American Statistical Association, 84(406), 414–420.

> Jaro, M. A. (1995). Probabilistic linkage of large public health data files. Statistics in Medicine, 14(5–7), 491–498.

> Jenkins, M.-C., & Smith, D. (2005). Conservative stemming for search and indexing. Proc. ACM SIGIR.

> Khan, S. I., & Hoque, A. S. Md. L. (2016). Similarity analysis of patients’ data: Bangladesh perspective. 2016 International Conference on Medical Engineering, Health Informatics and Technology (MediTec), 1–5. https://doi.org/10.1109/MEDITEC.2016.7835390

> Kiss, T., & Strunk, J. (2006). Unsupervised multilingual sentence boundary detection. Computational Linguistics, 32(4), 485–525.

> Lait, A. J., & Randell, B. (1996). An assessment of name matching algorithms. Technical Report Series-University of Newcastle Upon Tyne Computing Science.

> Levenshtein, V. I. (1966). Binary codes capable of correcting deletions, insertions, and reversals. Soviet Physics Doklady, 10(8), 707–710.

> Liang, F. M. (1983). Word Hy-phen-a-tion by Com-put-er. Calif. Univ. Stanford. Comput. Sci. Dept.

> Lovins, J. B. (1968). Development of a stemming algorithm. Mech. Translat. &Comp. Linguistics, 11(1–2), 22–31.

> Lynch, B. T., & Arends, W. L. (1977). Selection of a surname coding procedure for the SRS record linkage system. Washington, DC: US Department of Agriculture, Sample Survey Research Branch, Research Division.

> Monge, A. E., Elkan, C., & others. (1996). The Field Matching Problem: Algorithms and Applications. Kdd, 2, 267–270.

> Moore, G. B. (1977). Accessing individual records from personal data files using non-unique identifiers (Vol. 13). US Department of Commerce, National Bureau of Standards.

> Nakache, D. (2007). Extraction automatique des diagnostics à partir des comptes rendus médicaux textuels. Laboratoire CEDRIC-Équipe ISID. Paris, Conservatoire National Des Arts et Métiers, 219.

> Odell, M. K. (1956). The profit in records management. Systems (New York), 20, 20.

> Paice, C. D. (1990). Another stemmer. ACM Sigir Forum, 24(3), 56–61.

> Paternostre, M., Francq, P., Lamoral, J., Wartel, D., & Saerens, M. (2002). Carry, un algorithme de désuffixation pour le français. Rapport Technique Du Projet Galilei.

> Philips, L. (1990). Hanging on the metaphone. Computer Language, 7(12), 39–43.

> Philips, L. (2000). The double metaphone search algorithm. C/C++ Users Journal, 18(6), 38–43.

> Pollock, J. J., & Zamora, A. (1984). Automatic spelling correction in scientific and scholarly text. Communications of the ACM, 27(4), 358–368. https://doi.org/10.1145/358027.358048

> Postel, H. J. (1969). Die Kölner Phonetik. Ein Verfahren zur Identifizierung von Personennamen auf der Grundlage der Gestaltanalyse. IBM-Nachrichten, 19, 925–931.

> Ratcliff, J. W., & Metzener, D. E. (1988). Pattern-matching-the gestalt approach. Dr Dobbs Journal, 13(7), 46.

> Savoy, J. (1993). Stemming of French words based on grammatical categories. Journal of the American Society for Information Science, 44(1), 1–9.

> Savoy, J. (1999). A stemming procedure and stopword list for general French corpora. Journal of the American Society for Information Science, 50(10), 944–952.

> Schinke, R., Greengrass, M., Robertson, A. M., & Willett, P. (1996). A stemming algorithm for Latin text databases. Journal of Documentation.

> Shannaq, B. A., & Alexandrov, V. V. (2010). Using Product Similarity for Adding BusinessValue and Returning Customers. Global Journal of Computer Science and Technology.

> Smith, T. F., Waterman, M. S., & others. (1981). Identification of common molecular subsequences. Journal of Molecular Biology, 147(1), 195–197.

> Snae, C., & Diaz, B. (2002). An interface for mining genealogical nominal data using the concept of linkage and a hybrid name matching algorithm. Journal of 3D-Forum Society, 16(1), 142–147.

> Sørensen, T. (1948). A method of establishing groups of equal amplitude in plant sociology based on similarity of species content and its application to analyses of the vegetation on Danish commons.

> Tversky, A. (1977). Features of similarity. Psychological Review, 84(4), 327.

> Van Rijsbergen, C. J., Robertson, S. E., & Porter, M. F. (1980). New models in probabilistic information retrieval. British Library Research and Development Department London.

> Varol, C., & Bayrak, C. (2012). Hybrid matching algorithm for personal names. Journal of Data and Information Quality (JDIQ), 3(4), 1–18.

> Wilde, G., & Meyer, C. (1999). Doppelgänger gesucht: Ein Programm für kontextsensitive phonetische Textumwandlung. c.

> Winkler, W. E. (1990). String Comparator Metrics and Enhanced Decision Rules in the Fellegi-Sunter Model of Record Linkage.

---
title: 'Talisman: a JavaScript archive of fuzzy matching, information retrieval and record linkage building blocks'
tags:
  - javascript
  - fuzzy matching
  - natural language processing
  - phonetic algorithms
  - stemmers
  - inflectors
  - deduplication
  - record linkage
  - entity resolution
  - similarity metrics
  - information retrieval
  - search engines
  - tokenizers
authors:
  - name: Guillaume Plique
    orcid: 0000-0003-4916-8472
    affiliation: 1
affiliations:
 - name: médialab, SciencesPo Paris
   index: 1
date: 11 June 2020
bibliography: paper.bib
---

# Summary

Information retrieval [@manning2008introduction; @baeza1999modern] and record linkage [@fellegi1969theory; @christen_data_2012; @herzog_data_2007] have always relied on crafty and heuristical routines aimed at implementing what is often called *fuzzy matching*. Indeed, even if fuzzy logic feels natural to humans, one needs to find various strategies to coerce computers into acknowledging that strings, for instance, are not always strictly delimited. But if some of those techniques, such as the Soundex phonetic algorithm [@odell1956profit] invented at the beginning of the 20th century, are still well known and used, a lot of them were unfortunately lost to time.

As such, the **Talisman** JavaScript library aims at being an archive of a wide variety of techniques that have been used throughout computer sciences' history to perform fuzzy comparisons between words, names, sentences etc. Thus, even if **Talisman** obviously provides state-of-the-art functions that are still being used in an industrial context, it also aims at being a safe harbor for less known or clunkier techniques, for historical and archival purposes.

The library therefore compiles a large array of implementations of the following techniques:

* **keyers**: functions used to normalize strings in order to drop or simplify artifacts that could impair comparisons.
* **similarity metrics**: functions used to compute a similarity or distance between two strings, such as the Levenshtein distance [@levenshtein1966binary] or the Jaccard similarity [@jaccard1912distribution], etc.
* **phonetic algorithms**: functions aiming at producing a fuzzy phonetical representation of the given strings to enable comparisons such as the *Kölner phonetik* [@postel1969kolner] or the *Metaphone* [@philips1990hanging], etc.
* **stemmers**: functions reducing given strings to a *stem* to ease comparisons of a word's various inflections such as the Porter stemmer [@van1980new], etc.
* **tokenizers**: functions used to cut strings into relevant pieces such as words, sentences etc.

Those building blocks can then be used to perform and improve the following tasks:

* Building more relevant search engines through fuzzy matching and indexing
* Clustering string by similarity
* Record linkage, entity resolution etc.
* Natural language processing

Finally, this library can also be used to perform some benchmarks of those building blocks, wrt. precision, recall etc. of the fuzzy matches, which is seldom done in the literature because of how hard it can be to find comprehensive archives aggregating many phonetic algorithms, stemmers etc.

# Related works

* [abydos](https://github.com/chrislit/abydos): a python library implementing similar utilities.
* [java-string-similarity](https://github.com/tdebatty/java-string-similarity), [stringdistance](https://github.com/vickumar1981/stringdistance): Java libraries implementing string distance/similarity functions.
* [OpenRefine](https://openrefine.org/): a fully-fledged application designed to apply similar methods to typical data cleaning tasks.
* [clj-fuzzy](https://github.com/Yomguithereal/clj-fuzzy): a Clojure library which stands as an earlier version of **Talisman**

# References
