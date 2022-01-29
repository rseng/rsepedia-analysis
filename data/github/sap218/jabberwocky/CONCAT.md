---
title: 'Jabberwocky: an ontology-aware toolkit for manipulating text'
tags:
  - Python
  - Ontologies
  - Text
authors:
  - name: Samantha C Pendleton
    orcid: 0000-0002-6169-0135
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Georgios V Gkoutos
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Institute of Cancer and Genomic Sciences, University of Birmingham, UK
   index: 1
 - name: University Hospitals Birmingham NHS Foundation Trust, UK
   index: 2
date: 25 February 2020
bibliography: paper.bib

---

# Summary

Unstructured textual data is underused, as extracting the key textual elements is complicated by a lack of structured terms, e.g., collecting the sentences from a corpus that are discussing a particular topic. To extract valuable text about a topic from a corpus, a user will need to gather a set of related terms. For example, when analysing clinical documents we can extract sentences by using specific clinicial terms. However this can miss additional valuable sentences where synonyms are used instead (e.g., physician notes that use shorthand). By considering terms and their synonyms we can extract more sentences from a corpus, making more data available for analysis. One way to do this and represent our knowledge of terms associated with a domain is to create an ontology. Ontologies allow us to formalise our knowledge of a domain in a condensed manner by using controlled terms, called classes [@Hoehndorf2015-qr]. Classes can be annotated with metadata, including synonyms. Ontologies can include relationships between terms, and annotations such as cross-references to other ontologies [@Hoehndorf2015-qr].

Clearly, ontologies are valuable for the analysis of textual data. Unfortunately, despite the existence of many well-established ontologies, such as the "Human Phenotype Ontology" [@Robinson2008-jh] and the "Disease Ontology" [@Schriml2012-qp], there remains a lack of tools that can take advantage of ontologies, especially for general text manipulation. Existing tools for annotating text, such as “spaCy” [@Honnibal2017-dn], “tagtog” [@Cejuela2014-lv], and “Stanford CoreNLP” [@Manning2014-rt] cannot interrogate text with an ontology directly, and require ontologies to be pre-processed into other formats (leaving the time-consuming task of extracting labels and tags from an ontology into a suitable intermediate format as an exercise for the end-user). These are specialist tools, returning all text in the document with every word tagged, as “noun”, “verb”, and other customised tags. There exists a niche for users who want to leverage an ontology to retrieve textual data from a corpus without having to perform any pre-processing, or parse away unwanted tags.

We introduce Jabberwocky, a Python-based [@Van_Rossum1995-ia], open-source toolkit (accessible via https://github.com/sap218/jabberwocky) that allows users to query text in an ontology-aware fashion, and to modify those ontologies based on their findings. For example, with Jabberwocky’s ``catch`` command, a user provides textual data, their chosen ontology, and a set of classes from the ontology to use as search terms. Jabberwocky cleans the input text, collects the annotated synonyms for the user-specified target classes (using “Beautiful Soup” to read the ontology’s XML structure [@Richardson2007-ba]), and then returns the key elements (e.g., lines from a corpus) which match one of the target terms, or a synonym from the ontology. The ``catch`` command will help users retrieve more matches for their chosen terms from the corpus, without users having to explicitly define all the possible synonyms or alternative spellings beforehand.

Jabberwocky also helps ontology developers to iteratively improve their ontology. The ``bite`` command allows a user to provide textual data and rank the important terms by using the term frequency–inverse document frequency (tf-idf) method from “scikit-learn” [@Pedregosa2011-st], which calculates an importance metric for a term based on the frequency of its occurrence and the document size. Providing an ontology will exclude terms already described in the ontology, meaning the result of ``bite`` will be a CSV of candidate terms to potentially be added to the ontology, exported by “pandas” [@McKinney2010-xf]. Once an expert has reviewed the terms and associated them to a class in the ontology, Jabberwocky’s third command, ``arise``, will annotate the classes in the ontology, adding the newly identified synonyms. Iteratively performing multiple rounds of ``bite`` and ``arise`` can help the development and maintenance of ontologies. A user could use the ``catch`` command to confirm the modified ontology now captures more of the corpus.

Jabberwocky’s test repository (see Jabberwocky repo for further instructions), shows examples of each command separately. The ‘process’ directory shows an example that combines all three commands to demonstrate an example workflow. With 24 blog posts, the first use of ``catch`` returned 11 posts with the provided keywords. The example uses ``bite`` to review the CSV of ranked terms and curated new synonyms, simply by adding the corresponding class label from the ontology. It then uses ``arise`` to add the identified synonyms into the ontology. With the second round of ``catch`` the number of posts returned for the same keywords increased to 16. This is a basic and straightforward example, but powerful. With Jabberwocky, users can efficiently search their text and gain more instances, providing new insight.

Jabberwocky leverages the strength of ontologies and text for a wide range of tasks. It will be useful to users who want to manipulate textual data using controlled vocabulary from ontologies.

# Acknowledgements

Project was funded by the Medical Research Council (MRC) (MR/S502431/1) & supported by Health Data Research (HDR) UK (HDRUK/CFC/01).

# References

# Jabberwocky

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02168/status.svg)](https://doi.org/10.21105/joss.02168) [![DOI](https://zenodo.org/badge/227571502.svg)](https://zenodo.org/badge/latestdoi/227571502) 

**see [Jabberwocky site](https://sap218.github.io/jabberwocky/) for in-depth explanation and working scenarios (including test files)**

Jabberwocky is a toolkit for **ontologies**. Since we all know ontologies are "*nonsense*". Not enough tools existsing utilise the power of ontologies. Don't hesitate to create an [`issue`](https://github.com/sap218/jabberwocky/issues) or [`pull request`](https://github.com/sap218/jabberwocky/pulls) (see [**guidelines**](https://github.com/sap218/jabberwocky/blob/master/CONTRIBUTING.md) first).

#### Version

See `setup.py` in your local copy for version number | or `Releases`:
* **v1.0.0.0** [29/06/2020]
* **v2.0.0.0** [10/05/2021]
     - includes `spacy PhraseMatcher`
     - own synonym tags
     - plot output for tf-idf

##### Install
```
$ git clone https://github.com/sap218/jabberwocky
$ cd jabberwocky
$ python3 setup.py install --user
```
**note**: if you are using a virtual environment you can avoid `--user`

##### Prerequisites
```
$ pip3 install click BeautifulSoup4 scikit-learn pandas lxml pytest spacy matplotlib
```
or **after installing**, use the `requirements.txt` file:
```
$ pip3 install -r requirements.txt
```

#### Elements

command | description
------- | -----------
`bandersnatch` | extract synonyms from an RDF/XML syntax `OWL` ontology
`catch` | extract elements / sentences of text using key words
`bite`  | run statistical tf-idf for important words from text
`arise` | adding / updating new synonyms to an ontology

#### Ontology formats
`jabberwocky` works with the `OWL` ontology format: `RDF/XML` - for example, well-known biomedical ontologies such as `doid.owl`, `hpo.owl`, and `uberon.owl` will all work, plus your own created.

#### Examples
for examples of Jabberwocky's commands in use, please see the **[site](https://sap218.github.io/jabberwocky/SCENARIO.html)**.

**OR** to run the automated tests (in the cloned directory):
```
$ git submodule init
$ git submodule update
$ tox
```

---

## bandersnatch
`bandersnatch` curates synonyms for a list of key terms / or words of interest from an ontology of your choice, you provide a list of ontology synonym tags. **note**: it is recommended your list of keywords are exactly the classes from your chosen ontology (all in lowercase).
```
$ jab-bandersnatch -o hpo.owl -s ontology_synonym_tags.txt -k words_of_interest.txt
```

## catch
`catch` essentially "catches" key elements / sentences from textual data using a `.json` of key terms and their synonyms, you can use the outcome from `bandersnatch`. A user will also provide a `.txt` or `.json` of the text data. **note**: if a `.json` of text data is provided, you need specify the parameter for the field that contains the textual data to process.
```
$ jab-catch -k label_with_synonyms.json -t facebook_posts.json -p user-comment -i inner-user-comment-reply
```

## bite
`bite` runs a tf-idf statistical analysis: searching for important terms in a text corpus. a user can use a list of key terms to remove from the text in order to avoid being in the statistical model - meaning other terms may be ranked higher. **note**: again with `catch`, if you provide a `.json` of text data, you need specify the field that contains the textual data to process. Using `-g True` means you'll get a bar plot of the (default) 30-top terms.
```
$ jab-bite -k label_with_synonyms.json -t twitter_posts.txt -g True
```

## arise
`arise` inserts synonyms in an ontology: **you** define these synonyms (e.g. "exact", "broad", "related", or "narrow") - these new synonyms may be based on the tf-idf statistical analysis from `bite`.
```
$ jab-arise -o pocketmonsters.owl -f tfidf_new_synonyms.tsv
```

---

## Thanks! :dragon:

the poem "Jabberwocky" written by Lewis Carrol is described as a "nonsense" poem.

**Contributors** - thank you!
- [@majensen](https://github.com/majensen) for setting up automated testing w/ `pytest` - [see pull request #13 for more details](https://github.com/sap218/jabberwocky/pull/13)

**Citing**
```
@article{Pendleton2020,
  doi = {10.21105/joss.02168},
  url = {https://doi.org/10.21105/joss.02168},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {51},
  pages = {2168},
  author = {Samantha C. Pendleton and Georgios V. Gkoutos},
  title = {Jabberwocky: an ontology-aware toolkit for manipulating text},
  journal = {Journal of Open Source Software}
}
```

---

## ONE LAST THING...

You can combine these commands together to form a process of steps of ontology synonym development and text analysis - see the [SCENARIO](https://sap218.github.io/jabberwocky/SCENARIO.html) for a working example of this process.

![jabberwocky cycle](/images/cycle.jpg)
# Contributing Guildelines / Issues for jabberwocky :dragon_face:


## Contributing Code
* know how to fix a bug or want a new feature, consider either creating an issue or opening a pull request
* if i think it is suited i will merge to master branch
* frequent contributors will be added to a contributors list for thanks and acknowledgement
* **note**: additional code should be commented and w/ your username to acknowledge contribution, e.g. <br> 
```
new_list = []
for word in old_list:      # cycles through old list
  if len(word) > 10:       # if word greater than 10 characters
    new_list.append(word)  # append to new list
    print(word)            # printing word as users may want to see as a reference - @yourusername
``` 


## Issues
* don't hesitate to create an issue
* issues can be regarding anything: error reporting, feature request, or questions for help
* add an issue with a clear description & title
* try to know what error you are getting and make sure the input files are correct
* i will label the issue accordingly (see [`labels`](https://github.com/sap218/jabberwocky/labels) - or below (`bug`, `documentation`, `duplicate`, `help`, `request`, `wontfix`))
* please accept that sometimes information may not be fully understood so there could be follow-up questions
* if you have an issue with the `README` please do say so! i encourage help to make it better, perhaps you have a better ideas than me

#### bug
* error reporting if there is a clear issue
* if you know how to fix it yourself, consider doing it and creating a pull request
* e.g. "my list of ontology classes and synonyms is empty and I'm using HPO"

#### documentation
* any issues which i believe can be fixed by better documentation, i will label with `documentation`

#### duplicate
* if you have a question which you believe could have been asked before, ask anyway! I'll label as `duplicate`
* i will try my best to ensure a link to the original issue will be provided

#### help
* if you have a general question about how it works/etc. - check out [`jabberwocky-tests`](https://github.com/sap218/jabberwocky-tests)
* if the tests repository doesn't help, please ask away - no question is stupid
* e.g. "can i use my own created ontology?" (FYI: yes)

#### request
* have an idea as a new feature? please tell me! 
* if you know how to make it yourself, consider doing it and creating a pull request

#### wontfix
* sometimes people can confuse tools as something they are not - if there is an `bug` or `request` which i believe is not suitably designed for jabberwocky, i will label accordingly
* after labelling as `wontfix` i will comment why, giving you a few days (perhaps to give a rebuttal) before i close the issue from lack of response/activity
* please accept this label, if i think a different tool is doing what you want i will redirect you

# scenario

### >> go back to [main page](https://sap218.github.io/jabberwocky/)


### >> go to [`jabberwocky/test_files`](https://github.com/sap218/jabberwocky/tree/master/test_files) for data of the following examples

## Aim

You have extracted textual data: blog posts from a social media platform. These social media posts include varios individuals discussing a topic, which you are researching. In this scenario the users are talking about [*pokemon*](https://simple.wikipedia.org/wiki/Pok%C3%A9mon).

Some example posts from the text data:

> I think only gen 6 pokemon are on this path, try route 2 - wanderer wendy

> No thanks, I'm, trying to catch a flying type in the mountains with the clear air - trainer penelope

Your aim is to extract particular posts which individuals use specific terms, e.g. "dragon".

This is where **ontologies are useful**. Ontologies are a controlled set of vocabulary with terms logically related to the other, e.g. in anatomy our hand is a part of the arm. The purpose of Jabberwocky is looking at these terms and their synonyms, as in the example above that the arm has a synonym "upper limb". Current tools which exist for NLP don't include an ontology manipulation aspect, which Jabberwocky overcomes.

## Bandersnatch

[test_files](https://github.com/sap218/jabberwocky/tree/master/test_files/bandersnatch) for `bandersnatch`

You have access to `pocketmonsters.owl` an ontology with some concepts of pokemon, e.g. pokemon types. Below is a snippet of the ontology, looking at the class (label) "generation one", which has the synonym "gen one".
```
<owl:Class rdf:about="pocketmonsters#PM_00008">
	<rdfs:subClassOf rdf:resource="pocketmonsters#PM_00001"/>
        <oboInOWL:hasExactSynonym>generation 1</oboInOWL:hasExactSynonym>
        <oboInOWL:hasRelatedSynonym>gen 1</oboInOWL:hasRelatedSynonym>
        <oboInOWL:hasRelatedSynonym>gen one</oboInOWL:hasRelatedSynonym>
        <rdfs:label xml:lang="en">generation one</rdfs:label>
</owl:Class>
```
Notice the relevant synonym tags, to use `bandersnatch` you need to provide a newline-separated list of tags, `ontology_synonym_tags.txt`.
```
oboInOWL:hasExactSynonym
oboInOWL:hasRelatedSynonym
```
Finally, you have some terms of interest, `words_of_interest_for_ontology.txt`. Notice these terms are exactly the same labels from the ontology.
```
generation one
dragon
route
water
small
large
generation six
```
Using the command `bandersnatch`:
```
$ jab-bandersnatch -o pocketmonsters.owl -s ontology_synonym_tags.txt -k words_of_interest_for_ontology.txt
```
The output `output_ontology_label_synonyms.json` includes your terms of interest and their synonyms based on the synonym tags, some synonym tags are in different formats so it is important to investigate. If you don't have an ontology of interest, FOLLOWING THE SAME STYLE YOU SHOULD MAKE YOUR OWN `.json`. Below is an example of the output (the empty lists mean no synonyms for this term).
```
{
    "small": [],
    "large": [],
    "route": [],
    "generation one": [
        "generation 1",
        "gen 1",
        "gen one"
    ],
    "generation six": [
        "generation 6",
        "gen 6",
        "gen six"
    ],
    "water": [],
    "dragon": []
}
```

## Catch

[test_files](https://github.com/sap218/jabberwocky/tree/master/test_files/catch) for `catch`

Next, using the `catch` command, you provide the previous `bandersnatch` output: `output_ontology_label_synonyms.json` **OR** your own created `.json`, in the example you can see `own_created_word_w_synonyms.json` which includes different synonyms not in the ontology, e.g. "small" has the synonym "tiny". **NOTE**: for the remaining scenario, we will be using `output_ontology_label_synonyms.json`.

You have the social media posts, in the `test_files/catch` directory I provide two formats, `social_media_posts.txt` and `social_media_posts.json` - below is an example of the newline-separated unformatted `.txt`:
```
Any small pokemon nearby? I need to catch a Metapod!
I think only gen 6 pokemon are on this path - try route 2.
```
The `.json` example is below - when using this formatted version, you will need to provide the parameter for user comment / text and possible inner-comments / replies. Notice below `post` is used for a user's comment.
```
{
"thread_one":[
	{"name": "bug catcher joe", "post": "Any small pokemon nearby? I need to catch a Metapod!"},
	{"name": "wanderer wendy", "post": "I think only gen 6 pokemon are on this path - try route 2."}
...
```
Below is an example of running `catch`:
```
$ jab-catch -k output_ontology_label_synonyms.json -t social_media_posts.json -p post -i reply
```
The outputs include `output_terms_match_raw.txt` which include a `.txt` file of newline-separated posts which included one of the terms of interest - this output has **9 posts** (remember this number)
```
what route is best for small normal pokemon my skitty needs a friend
any small pokemon nearby i need to catch a metapod
```
Additionally, the posts formatted as with their terms of interest, `output_terms_match.json`: 
```
{
    "small": [
        "any small pokemon nearby i need to catch a metapod",
...
```

## Bite

[test_files](https://github.com/sap218/jabberwocky/tree/master/test_files/bite) for `bite`

But what if there are synonyms which weren't present in the ontology? Or synonyms which you didn't previously consider? The statistical tf-idf method scores each word in a corpus based on the frequency in the document - essentially picking out the "important terms".

Following the same social media input files as `catch`, you provide a `social_media_posts.txt` or `social_media_posts.json` (with parameters) - in addition to `output_ontology_label_synonyms.json` which removes from the posts in order to increase the scores of the tf-idf statistical analysis, however you don't need to include this file and rather you could investigate how these terms present in the results.

Below is an example of running `bite`:
```
$ jab-bite -k output_ontology_label_synonyms.json -t social_media_posts.txt -g True
```
The output includes `tfidf_results.tsv`, which is a tab-separated file of terms and their tf-idf score, as seen below:
```
words	count
route	2.27975328153784
path	1.8866845905817748
pokemon	1.8215522309449206
...
```
You can also request a plot of term scores as a bar plot and set the term limit (`-l`) - the plot is saved as: `tfidf_plot.pdf`.

## Arise

[test_files](https://github.com/sap218/jabberwocky/tree/master/test_files/arise) for `arise`

Looking at the results from the `tfidf_results.tsv`, you curate synonyms of interest and create: `tfidf_new_synonyms.tsv`. This `.tsv` includes the synonym you wish to inject / update in the ontology, with the corresponding class label, and the type of synonym.
```
synonym	class	type
path	route	oboInOWL:hasRelatedSynonym
evolve	generation	oboInOWL:hasBroadSynonym
```
Below is an example of running `arise` whilst providing the original `pocketmonsters.owl`:
```
$ jab-arise -o pocketmonsters.owl -f tfidf_new_synonyms.tsv
```
The output is an updated ontology `updated-ontology.owl`. See below the broad synonym "evolve" addition to "generation".
```
<owl:Class rdf:about="pocketmonsters#PM_00001">
	<oboInOWL:hasBroadSynonym>evolve</oboInOWL:hasBroadSynonym>
	<rdfs:label xml:lang="en">generation</rdfs:label>
</owl:Class>
```

## **Lets go back to the beginning**

## Bandersnatch ROUND 2

Using `bandersnatch` with `updated-ontology.owl`, the `output_ontology_label_synonyms.json` output will include the new synonyms, for example "route" now has the synonym "path". For this example, I renamed it to `roundtwo_output_ontology_label_synonyms.json`.
```
$ jab-bandersnatch -o updated-ontology.owl -s ontology_synonym_tags.txt -k words_of_interest_for_ontology.txt
```

## Catch ROUND 2

Finally, using `catch`, with the updated `roundtwo_output_ontology_label_synonyms.json` - the output posts should be increased due to more synonyms, from **9 posts to 13**.
```
$ jab-catch -k roundtwo_output_ontology_label_synonyms.json -t social_media_posts.txt
```
### >> see [SCENARIO](SCENARIO.md) for a working tutorial

Ontologies are useful for their condense and structured format of a domain of knowledge. Specifically their organised terms and corresponding synonyms. Many NLP tools don't utilize ontologies, Jabberwocky uses ontologies for synonym curation. Here we provide a full-depth explanation, informative scenarios, and working examples for the Jabberwocky toolkit - for installation instructions, see the [Jabberwocky](https://github.com/sap218/jabberwocky) repository. 

## About the Commands

Below is an-indepth explanation of the commands which you can use with Jabberwocky.

## bandersnatch
`bandersnatch` curates synonyms for a list of key terms / or words of interest from an ontology of your choice, you provide a list of ontology synonym tags. **note**: it is recommended your list of keywords are exactly the classes from your chosen ontology (all in lowercase).

#### Usage
```
$ jab-bandersnatch --help
Usage: jab-bandersnatch [OPTIONS]

Options:
  -o, --ontology TEXT     file of ontology.  [required]
  -s, --synonymtags TEXT  list of XML tags for synonym curation.  [required]
  -k, --keywords TEXT     list of class labels you want to use to search.
                          [required]
  --help                  Show this message and exit.
```
#### Running
```
$ jab-bandersnatch --ontology pocketmonsters.owl --synonymtags ontology_synonym_tags.txt --keywords words_of_interest.txt
$ jab-bandersnatch -o pocketmonsters.owl -s ontology_synonym_tags.txt -k words_of_interest.txt
```

###### Output
* a `.json` file: `output_ontology_label_synonyms.json` of the classes and synonyms for your reference - this can be used for the `catch` command

---

## catch
`catch` essentially "catches" key elements / sentences from textual data using a `.json` of key terms and their synonyms, if not using the `output_ontology_label_synonyms.json` from `bandersnatch`, then you can provide your own. The main element of `catch` is the textfile, which can be `.txt` or `.json` - if a `.json` is provided you need specify the parameter for the field that contains the textual data to process.

#### Usage
```
$ jab-catch --help
Usage: jab-catch [OPTIONS]

Options:
  -k, --keywords TEXT        list of terms and synonyms you want for grep, can
                             be from the ontology output.  [required]
  -t, --textfile TEXT        JSON or TXT file of text you want annotate.
                             [required]
  -p, --parameter TEXT       parameter of the the JSON text data.
  -i, --innerparameter TEXT  inner parameter of the the JSON text data if
                             expecting replies.
  --help                     Show this message and exit.
```
#### Running
```
$ jab-catch --keywords output_ontology_label_synonyms.json --textfile example_textfile.txt
$ jab-catch -k own_labels_synonyms.json -t example_tweets.json -p tweet-comment -i tweet-reply
```

###### Output
* a `.txt` file: `output_terms_match_raw.txt` which includes all elements / sentences from the text file which includes a term of interest
* a `.json` file: `output_terms_match.json` which includes the posts for each word of interest

---

## bite
`bite` runs a tf-idf statistical analysis: searching for important terms in a text corpus. a user can use a list of key terms to remove from the text in order to avoid being in the statistical model - meaning other terms may be ranked higher. **note**: with the `.json` input you need specify the field inside the JSON that contains the textual data to process (same as `catch`).

#### Usage
```
$ jab-bite --help
Usage: jab-bite [OPTIONS]

Options:
  -k, --keywords TEXT        list of terms and synonyms you want to remove
                             from tf-idf analysis.
  -t, --textfile TEXT        JSON or TXT file of text you want annotate.
                             [required]
  -p, --parameter TEXT       parameter of the the JSON text data.
  -i, --innerparameter TEXT  inner parameter of the the JSON text data if
                             expecting replies.
  -g, --graph TEXT           make True if you want a plot of top 30 terms.
  -l, --limit TEXT           change if want a different plot limit.
  --help                     Show this message and exit.

```
#### Running
```
$ jab-bite --textfile facebook_posts.txt 
$ jab-bite -k output_ontology_label_synonyms.json -t example_tweets.json -p tweet-comment -i tweet-reply
$ jab-bite -k own_labels_synonyms.json -t facebook_posts.txt -g True
```

###### Output
* a `.tsv` file: `tfidf_results.tsv` of all terms and their tf-idf score
* a `.pdf` file: `tfidf_plot.pdf` the plot output which is requested if a user makes `--graph True` and presents the (default) 30-top scoring terms

---

## arise
`arise` inserts synonyms in an ontology: **you** define these synonyms (e.g. "exact", "broad", "related", or "narrow") - these new synonyms may be based on the tf-idf statistical analysis from `bite`.

#### Usage
```
$ jab-arise --help
Usage: jab-arise [OPTIONS]

Options:
  -o, --ontology TEXT  file of ontology.  [required]
  -f, --tfidf TEXT     TSV file of the synonyms you want to add, can be based
                       from the tf-idf results.  [required]
  --help               Show this message and exit.
```
#### Running
```
$ jab-arise --ontology pocketmonsters.owl --tfidf tfidf_new_synonyms.tsv
$ jab-arise -o pocketmonsters.owl -f tfidf_new_synonyms.tsv
```

###### Output
* a `.owl` file: `updated_ontology.owl`


