# ConTEXT-Explorer
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Build Status](https://travis-ci.com/alicia-ziying-yang/conTEXT-explorer.svg?branch=main)](https://travis-ci.com/alicia-ziying-yang/conTEXT-explorer)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5762643.svg)](https://doi.org/10.5281/zenodo.5762643)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03347/status.svg)](https://doi.org/10.21105/joss.03347)

**ConTEXT Explorer** is an open Web-based system for exploring and visualizing concepts (combinations of co-occurring words and phrases) over time in the text documents. **ConTEXT Explorer** is designed to lower the barriers to applying information retrieval and machine learning for text analysis, including:
- preprocessing text with sentencizer and tokenizer in a **Spacy pipline**;
- building **Gensim** word2vec model for discovering similar terms, which can be used to expand queries;
- indexing the cleaned text, and creating a search engine using **Whoosh**, which allows to rank sentences using the Okapi BM25F function;
- visualizing results across time in interactive plots using **Plotly**.

It is designed to be user-friendly, enabling researchers to make sense of their data without technical knowledge. Users may:

- upload (and save) a text corpus, and customize search fields;
- add terms to the query using input from the word2vec model, sentence ranking, or manually;
- check term frequencies across time;
- group terms with "ALL" or "ANY" operator, and compound the groups to form more complex queries;
- view results across time for each query (using raw counts or proportion of relevant documents);
- save and reload results for further analysis; 
- download a subset of a corpus filtered by user-defined terms.

More details can be found in the user manual below.

## How to install
### Get the app
Clone this repo to your local environment:

    git clone https://github.com/alicia-ziying-yang/conTEXT-explorer.git

## Set up environment
ConTEXT Explorer is developed using Plotly Dash in **Python**. We are using `Python 3.7.5` and all required packages listed in `requirement.txt`. To help you install this application correctly, we provide a conda environment file `ce-env.yml` for you to set up a virtual environment. Simply enter the folder:

    cd conTEXT-explorer
    
and run:

    conda env create -f ce-env.yml
    
To activate this environment, use:

    conda activate ce-env

### Install the application
Then, ConTEXT Explorer can be easily installed by:

    pip install . 
    
### Run the app
- If you want to run ConTEXT Explorer on your **local computer**, comment the code for ubuntu server, and uncomment the last line in `app.py`:

       # app.run_server(debug=False, host="0.0.0.0") # ubuntu server    
       app.run_server(debug=False, port="8010") # local test           

  To start the application, use:

      start-ce
        
  or

      python app.py
    
  The IP address with app access will be displayed in the output.
      
- If you want to run ConTEXT Explorer on an **ubuntu server**, use:

      nohup python app.py &
  


## How to use
A [sample corpus](https://github.com/alicia-ziying-yang/conTEXT-explorer/blob/main/doc/sample_data.csv) with a saved analysis is preset in this app. Feel free to explore the app features using this example. Please check more details in the manual below.

[Click here to view the paged PDF version](https://github.com/alicia-ziying-yang/conTEXT-explorer/blob/main/doc/conTEXT_explorer_ui_manual.pdf)

![alt text](https://github.com/alicia-ziying-yang/conTEXT-explorer/blob/6b1e79e2068eb284a132493815f35e57b3fec409/doc/conTEXT_explorer_ui_manual.png?raw=true)

## Contact and Contribution
This application is designed and developed by Ziying (Alicia) Yang, Gosia Mikolajczak, and Andrew Turpin from the [University of Melbourne](https://www.unimelb.edu.au/) in Australia.

If you encounter any errors while using the app, have suggestions for improvement, or want to contribute to this project by adding new functions or features, please [submit an issue here](https://github.com/alicia-ziying-yang/conTEXT-explorer/issues/new) and pull requests.
---
title: "ConTEXT Explorer: a web-based text analysis tool for exploring and visualizing concepts across time"
tags:
  - Python
  - Dash
  - Data Analysis
  - Data Visulization
authors:
  - name: Ziying Yang
    orcid: 0000-0001-7705-3280
    affiliation: 1
  - name: Gosia Mikolajczak
    orcid: 0000-0002-7386-4155
    affiliation: 1
  - name: Andrew Turpin
    affiliation: 1
affiliations:
  - name: University of Melbourne
    index: 1

date: 19 April 2021
bibliography: paper.bib
---

# Summary

**ConTEXT Explorer** is an open Web-based system that assists in exploring the context of concepts (combinations of co-occurring words and phrases) over time in text documents. It provides a user-friendly interface for the analysis of user-provided text data and integrates functionalities of the Whoosh search engine, Spacy, Gensim, and Plotly Python libraries. By providing suggestions for query expansion and producing interactive plots, `ConTEXT Explorer` facilitates exploratory data analysis, which can serve as the basis for subsequent text classification.

# Statement of need

With the availability of digital sources of data and associated tools, automated text analysis is becoming increasingly popular in the humanities and social sciences. While for very large corpora, unsupervised text mining methods like topic modelling [@lda] can provide some useful summaries of data, many social science and humanities applications require analysis of data in context. That is, simple "bags of words" automatically mined and presented in isolation from the original text are often not meaningful for complex questions involving human behaviour and society. Inevitably, human interpretation is required to make sense of such patterns. For corpora with more than several hundred documents, there is a need for computational tools that can assist researchers in exploring the context in which "bags of words" (we will call them <em>concepts</em> from now on) occurs.

Similarly, there is a need for tools that assist in the manual construction of concepts from text corpora. In particular, manual intervention to judge the semantic intent of words (e.g., word sense disambiguation) is usually needed to filter keywords to add to concepts that might be generated by automatic methods such as query expansion [@QE] or comparison of word embeddings [@word2vec]. For example, if a researcher is interested in finding articles about same sex marriage, they might start with "same_sex marriage"[^1] as a concept. Automated methods processing a corpus of news articles might suggest related words like 'matrimony', 'union', 'erosion', and 'puzzlement'[^2]. Depending on the research question and the context of these words, some might be relevant to the concept and should be included, while others are not. It requires complex human judgement to make the distinction. `ConTEXT Explorer` is a tool to assist the construction of such concepts in context.

[^1]: We use underscore to join multiple words into a single phrase.
[^2]: This is an example where we apply `ConTEXT Explorer` in the Australian Research Council Discovery Project (DP180101711) "Understanding Political Debate and Policy Decisions Using Big Data".

Most existing computational methods underlying automated text processing require at least a working knowledge of relevant methods and programming languages (such as R or Python). `ConTEXT Explorer` is designed to lower these barriers to entry, particularly for humanities and social science researchers, by allowing an application of information retrieval and machine learning methods to text analysis without programming knowledge.

# Comparison with other tools

Current text analysis tools require either previous knowledge of programming (e.g., R, Python), or are commercial products (e.g., **RapidMiner** [@rapidminer], **Google Cloud Natural Language API** [@googlenlp]). One exception that we are aware of is **Voyant Tools** [@voyant], which is an open-source web-based text analysis tool built in Java. It allows the users to explore their data using some basic text analysis techniques such as word frequencies (at the document level), word cloud, and word context (words appearing around a chosen term). `ConTEXT Explorer` provides several functionalities that give a user a deeper understanding of the text, which are currently not available in Voyant Tools, such as concept suggestions, sentence ranking and concept grouping.
It includes models allowing discovery of similar terms, and a search engine allowing retrieval of sentences relevant to concept terms, which can be used for concept expansion, and visualization of concepts over time. One key feature is that a concept can be either conjunction or disjunction of bags of words.

Compared to commercial text analysis systems such as **RapidMiner** [@rapidminer], which include some complex analysis techniques, `ConTEXT Explorer` is open-source (free) and easy to install. It enables researchers to browse text interactively for concepts (bags of words) in their corpus before mining the text in machine learning-driven systems.

`ConTEXT Explorer` is designed to help users interested in defining concepts, and exploring their trends over time. This could be particularly helpful as an input for some popular text analysis systems such as **MonkeyLearn** [@monkeylearn], which enable text classification, tagging, and training AI machine learning models but require prior knowledge of the data.

`ConTEXT Explorer` is engineered to allow for the integration of other Python packages into the analysis process. It can be easily combined with other Python APIs (such as MonkeyLearn), once the concept groups are defined.

# Key features

`ConTEXT Explorer` is developed using **Dash** [@plotly] in Python, and integrates the following packages.

- **Spacy** pipeline [@spacy] - for pre-processing the text corpora uploaded by users.
- **Whoosh** [@whoosh] - for building a search engine, which allows ranking of sentences relevant to the given concepts, and word frequency analysis at the sentence and document level.
- **Gensim** [@gensim] - for training a word2vec [@word2vec] model for the uploaded corpus, which allows the user to find words related to other words for expanding concepts.
- **Plotly** [@plotly] - for visualizing results in interactive graphs, which can be customized and saved as PNG files.

`ConTEXT Explorer` has been tested locally under macOS and as a server running under Ubuntu and continuous integration tests are performed using Travis CI.

## Build a corpus

Users are asked to format their text documents as a CSV file (with each document saved in a separate row), before uploading this file into `ConTEXT Explorer`. At a minimum, users are asked to provide document text and publication year. Users can also upload columns with additional document information (such as document author, title, and so on).

`ConTEXT Explorer` processes the submitted file in the following steps.

1. Sentenize and tokenize English text using Spacy [@spacy]. This allows ranking of documents and speeds up document search.
2. Index the documents, and build a search engine for the corpus using Whoosh [@whoosh] and the Okapi BM25F [@BM25F] ranking function.
3. Remove stop words, lemmatize remaining words, and generate a word2vec [@word2vec] model for the corpus using Gensim [@gensim].

For each corpus, users can create a new analysis, or load a pre-saved analysis to the dashboard.

## Dashboard

![The starting dashboard.](https://paper-attachments.dropbox.com/s_BF58715651395C8B59D508B9A7AFBDF87128C0D6732F3C5CB80FFC81F0067860_1618206868822_overview.png)

As shown in Figure 1, the dashboard interface has two panes. On the left-hand side, users can:

- select the year range of documents to be displayed in search results;
- add or delete query terms (single words or phrases) to create a concept;
- save the current query as a new analysis; and
- download the subset of the corpus filtered by the query terms.

**Overview**. The overview tab summarizes the corpus information such as the total number of documents, year range, document length, most frequent words in the corpus, and most frequent values for selected metadata.

![The sentences tab of the dashboard, with some query terms shown in the left pane.](https://paper-attachments.dropbox.com/s_BF58715651395C8B59D508B9A7AFBDF87128C0D6732F3C5CB80FFC81F0067860_1624346302562_Screen+Shot+2021-06-22+at+5.18.01+pm.png)

**Sentences**. This tab shows the ranking of relevant sentences based on query terms defined in the left pane. Sentences are ranked by the Okapi BM25F ranking function, and the computed similarity score for each sentence is shown in the "SCORE" column. The table can be sorted and filtered by column values. Users can click on each sentence to see its full content in a pop-up window, which also allows checking of the frequency of individual terms and adding them to the query.

![The grouping tab of the dashboard, showing the term frequency of the added terms across time (top), and some examples of query groups (bottom).](https://paper-attachments.dropbox.com/s_BF58715651395C8B59D508B9A7AFBDF87128C0D6732F3C5CB80FFC81F0067860_1618281338713_grouping.png)

**Grouping**. The top part of this tab shows the number of sentences containing each query term within the user-defined year range. In the bottom part, users can group the query terms using "Any" or "All" operators. Groups can be further combined into more complex groups.

!['Graphs' tab, showing the aggregated graph for all groups (top) and individual graph for each group (bottom).](https://paper-attachments.dropbox.com/s_BF58715651395C8B59D508B9A7AFBDF87128C0D6732F3C5CB80FFC81F0067860_1618282796027_graphs.png)

**Graphs**. Based on the query groups generated in the previous tab, this page displays aggregated and individual plots, which allow comparing groups (top) and individual terms within each group (bottom). Users can choose the number of relevant documents, the number of sentences, or the proportion of documents as the y-axis of the graphs. All graphs are plotted by Plotly [@plotly] which allows users to interact with every trace in the graphs.

## Save and reload an analysis

As mentioned in the section above, users are able to save the details of their analysis (including added terms and generated groups) and reload it to view all of the ranking, groups and graphs from the index page.

# Acknowledgements

The development of `ConTEXT Explorer` has been supported by the Australian Research Council Discovery Project (DP180101711) “Understanding Political Debate and Policy Decisions Using Big Data” awarded to the third author.

# References
