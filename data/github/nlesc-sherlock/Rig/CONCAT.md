# Rig

Every data analysis project starts with exploring and cleaning data. For
tabular data, there are some tools available that facilitate data pre-processing, such
as [OpenRefine](http://openrefine.org/) and [Trifacta Wrangler](https://www.trifacta.com/products/wrangler/).
However, these tools one big disadvantage; they don't scale to really big data
sets. Also, Trifacta Wrangler is not open source.
For full text data, no tools for data cleanup exist.

The goal of Rig is to provide a data cleaning tool that is open source, supports both tabular data and full text pre-processing, and
scales to big data sets. The frontend is a web interface that supports loading
the data and creating workflows. Based on user input, it generates scripts that
are send to the backend for execution.
The backend relies on [spark](http://spark.apache.org/) to ensure scalability.

Currently, we ware working on the use case of vocabulary cleanup for full text data. When working with full text data,
pre-processing often consists of tokenizing texts (i.e., splitting them into words),
lematizing or stemming the tokens, and performing some kind of filtering to remove
typo's and other unwanted tokens. Rig should support
customizing these tasks, and provide visualizations that facilitate vocabulary
cleanup.

![Rig screenshot](/app/images/rig-screenshot.png?raw=true)

## Install

### Front end

Install nodejs. See instructions at http://nodejs.org/. Make sure to install nodejs, nodejs-legacy and npm.

Install bower and grunt-cli globally:
```
sudo npm install -g bower grunt-cli
```

Fetch git repository
```
git clone https://github.com/nlesc-sherlock/Rig.git
```

Setup with bower
```
cd Rig
npm install
bower install
```
If you already have a installed the bower packages before, but need to update them for a new version of the code, run
```
bower update
```

### Webserver

Install requirements
```
pip install -r backend/requirements.txt
```

The webserver loads the data (this takes a few seconds). This means that the webserver needs to know where the data is.

1. Download and extract the [data](https://nlesc.sharepoint.com/sites/sherlock/_layouts/15/Group.aspx?GroupId=6aad52c4-7dfc-4076-9772-4f9c9180bde2&AppId=Files&id=%2Fsites%2Fsherlock%2FShared%20Documents%2Fdatasets%2Fenron-rig)
2. Copy `.rigrcEXAMPLE` to `.rigrc`
3. In this file, update `path=/path/to/data/dir` to the directory that  contains the data.

## Run

### Front end

Start development server & open browser
```
grunt serve
```
Changes made to code will automatically reload web page.

### Webserver

Run the webserver:
```
python backend/server.py
```

## Test

### Front end

#### Run unit tests

```
grunt test
```
Generates test report and coverage inside `test/reports` folder.

#### Run end-to-end tests with local browser (chrome)

Tests in Chrome can be run with
```
grunt e2e-local
```

#### Run end-to-end tests on [sauce labs](https://saucelabs.com/)

To connnect to Sauce Labs use sauce connect program. [Here](https://docs.saucelabs.com/reference/sauce-connect/) you can find the details on how to install and run it.

Before tests can be run the sauce labs credentials must be setup

```
export SAUCE_USERNAME=<your sauce labs username>
export SAUCE_ACCESS_KEY=<your sauce labs access key>
```

Tests in Chrome, Firefox on Windows, Linux and OSX can be run with
```
grunt e2e-sauce
```

Travis-ci also runs end-to-end tests on sauce labs.

Note! Running `grunt e2e-sauce` will undo all changes in `app/` folder.
# IPython Notebooks with code blocks
To create a block, add a cell to a pyspark iPython notebook that does the task that it needs to do. Then add the metadata of that cell and copy from the metadata_template.json.
The fields of the metadata are the following:
- `type` : if the type is `block`, the cell is regarded as a block and put in the generated notebook.
- `input`: the variable name of the input spark DataFrame (optional)
- `output`: the variable name of the output spark DataFrame
- `arguments`: list of variable names that are arguments to your block
- `name`: (unique) name of the block


To test if your block works correctly, you can add cells to your Notebook that are not of type `block`. These are ignored in the Notebook generation.
# NLTK Loader Example

All commands should be run in the `nltk-loader-example` directory.

* Build: `docker build -t rig-docker .`
* Run: `docker run -d -p 8888:8888 -v /home/jvdzwaan/code/Rig/nltk-loader-example:/home/jovyan/work -v /home/jvdzwaan/data/enron-seq:/data --name rig rig-docker start-notebook.sh` (Please provide proper directory mappings)
  * Open the notebook server at http://localhost:8888/
