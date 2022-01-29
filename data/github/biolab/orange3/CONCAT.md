Orange development
==================

The source code of [Orange] is versioned in [Git] and hosted on [GitHub]. 
If you want to contribute to this open-source project you will have to use git. However, for minor experimentation with the source code you can also get by without. 

[Orange]: https://orange.biolab.si/
[Git]: https://git-scm.com/
[GitHub]: https://github.com/biolab/orange

Prerequisites
-------------

[Orange] is written mostly in Python, therefore you'll need [Python 3] version 3.6 or newer.

You will also need a C/C++ compiler. On Windows, you can get one by installing [Visual Studio].
A slightly more "minimalistic" option is to install only its [Build Tools].

[Python 3]: https://www.python.org
[Visual Studio]: https://visualstudio.microsoft.com/vs/
[Build Tools]: https://visualstudio.microsoft.com/visual-cpp-build-tools/

Source code
-----------

Get the source code by cloning the git repository

    git clone https://github.com/biolab/orange3.git

or, alternatively, download and unpack the [ZIP archive] of the source code from [GitHub].

[ZIP archive]: https://github.com/biolab/orange3/archive/master.zip

Building
--------

Consider using virtual environments to avoid package conflicts. 

Install the required Python packages

    pip install -r requirements.txt
    
and run the setup script with a development option, which will link to the source code instead of creating a new package in Python's site-packages.

    python setup.py develop
    
Verify the installation by importing the Orange package from Python and loading an example Iris dataset.

    >>> import Orange
    >>> print(Orange.data.Table("iris")[0])
    [5.1, 3.5, 1.4, 0.2 | Iris-setosa]

Using the graphic user interface requires some additional packages.

    pip install -r requirements-gui.txt

To start Orange GUI from the command line, run:

    python3 -m Orange.canvas

Contributing
------------

If you've made improvements that you want to contribute, you'll need your own fork of the [GitHub] repository. After committing and pushing changes to your fork, you can create a pull request. We will review your contribution and hopefully merge it after any potential corrections. 

You can view the list of open [pull requests] and known [issues] on GitHub.

[pull requests]: https://github.com/biolab/orange3/pulls
[issues]: https://github.com/biolab/orange3/issues
# Orange Roadmap

In the next few years, we plan to expand Orange in the following directions:

- __Scalability__. Analyzing very large data in Orange is a challenge. Orange is an interactive data analytic tool with a flexible, user-defined workflow, which makes optimization for large data more complex than in the case of tools, which are either non-interactive or have a fixed data processing flow. We plan to address the following issues:
    - __Memory efficiency__. Orange can currently address data of up to around 500 MB.
    - __Speed-up of complex operations__. Hour-long processing is not acceptable for interactive data analysis.
    - __Interruptibility of processes and progress updates__. Operations that take too long are not always interruptible and may require shutting down the application; those same operations do not provide updates on progress reports.

- __Smoothness of the learning curve for new users__. This is not about fixing a problem but about building on Orange’s existing strength.
    - __Provide a tutorial for first-time users__ upon starting Orange. An excellent example is a guide in Open Street Map editor.
    - __Simplification of the interface__. Base Orange currently has under 200 components, and we keep this number low by design. Still, beginners should see fewer and be provided with assistance in further exploring.
    - __Gamification__. We are already exploring ways to implement self-paced tutorials within Orange and design ways to incorporate guidance for workflow construction, problem-driven exercises, and progress badges.
    - __Video content__. We will update and refresh existing [YouTube videos](http://youtube.com/orangedatamining) and record new ones.
    - __Lecture notes__. We have lecture notes for almost a hundred hands-on workshops carried out around the world. They were typeset in Pages, which worked nicely in the beginning, but now the management of about 300 pages of material became cumbersome. We are converting the lecture notes to LaTeX, enabling fast assembly of notes for a custom-designed hands-on workshop. The current solution needs improvements in design, and text needs updates and additional proofreading.

- __Improved functionality__ in several application fields, with particular priority being
    - __Spectroscopy__: we are making continuous progress with an add-on, thanks to support from [Synchrotron SOLEIL](https://www.synchrotron-soleil.fr/en).
    - __Text mining__: our goal is to design components for document characterization, summarization, and explanation of point-based visualizations. This work is in progress, thanks to the support by the [Slovene Ministry of Public Administration](https://www.gov.si/en/state-authorities/ministries/ministry-of-public-administration/).
    - __Network analysis__: speed-up and beautify existing visualizations and extend the toolbox with new network analysis procedures. This is being done in collaboration with the [International Laboratory for Applied Network Research, Higher Schools of Economics, Moscow](https://anr.hse.ru/en/).
    - __Time series analysis__: replace HighCharts with Qt-based visualizations, enable simultaneous analysis of a group of signals, improve modeling procedures, implement model evaluation techniques, and incorporate deep learning and auto-encoding for time series characterization, feature engineering, and prediction.

We need to address __funding__. The development of Orange relies on a core team sponsored by research and development grants. We will continue searching for new funding opportunities and seeking international collaborations to enrich Orange with new functionalities and improved interfaces. We seek funds both in developing Orange as a tool to democratize data science, assist in training machine learning, and support problem-solving in science and industry.
Change Log
==========

[next] - TBA
------------


[3.31.1] - 2022-01-07
--------------------
##### Bugfixes
* Group by: compute mode when all values in group nan ([#5763](../../pull/5763))
* Support numpy 1.22 ([#5760](../../pull/5760))
* Unpickling pre-3.28.0 Transformation ([#5759](../../pull/5759))


[3.31.0] - 2021-12-17
--------------------
##### Enhancements
* oweditdomain: Indicate variables in error state ([#5732](../../pull/5732))
* Scatterplot: Use opacity for contrast ([#5684](../../pull/5684))
* Feature Constructor: Evaluate categorical variables to strings ([#5637](../../pull/5637))
* New widget: Group By  ([#5541](../../pull/5541))
* Table lock: tests run with tables that are read-only by default ([#5381](../../pull/5381))
* config: sort example workflows ([#5600](../../pull/5600))

##### Bugfixes
* Paint Data: Fix ClearTool's issued commands ([#5718](../../pull/5718))
* pandas_compat: fix table_from_frames for "normal" dataframe ([#5652](../../pull/5652))
* pandas_compat: do not parse column of numbers (object dtype) to datetime ([#5681](../../pull/5681))
* HeatMap: Color gradient center value edit ([#5647](../../pull/5647))
* Distance Matrix: Fix crash on numeric meta vars as labels ([#5664](../../pull/5664))
* Fix running OWS with widgets with WebView in Orange.canvas.run ([#5657](../../pull/5657))
* main: Fix `--clear-widget-settings` parameter declaration ([#5619](../../pull/5619))


[3.30.2] - 2021-10-27
--------------------
##### Bugfixes
* Projections: Fix color density for continuous color palettes ([#5665](../../pull/5665))
* Fixes for scikit-learn 1.0 ([#5608](../../pull/5608))
* table_from_grames: fix indices parsing ([#5620](../../pull/5620))
* Fix overflow in bin calculations for time variables.  ([#5667](../../pull/5667))
* Variable: fix timezone when parsing time variable ([#5617](../../pull/5617))
* Require widget-base 4.15.1 and canvas-core 0.1.23 to fix some bugs/crashes


[3.30.1] - 2021-09-24
--------------------
##### Bugfixes
* OWTable: fix select whole rows regression ([#5605](../../pull/5605))


[3.30.0] - 2021-09-22
--------------------
##### Enhancements
* OWPythonScript: Better text editor ([#5208](../../pull/5208))
* PCA: Output variance of components ([#5513](../../pull/5513))
* Curve Fit: New widget ([#5481](../../pull/5481))
* Hierarchical Clustering: Annotate variables with clusters ([#5514](../../pull/5514))
* Create nodes on canvas drag/drop ([#5031](../../pull/5031))

##### Bugfixes
* setup.py: do not overwrite conda's PyQt5 ([#5593](../../pull/5593))
* Use explicit ordered multiple inputs ([#4860](../../pull/4860))
* Prevent crash when saving in unsupported format ([#5560](../../pull/5560))
* owrocanalysis: Fix test for non empty points array ([#5571](../../pull/5571))
* pandas_compat: fix conversion of datetime series ([#5547](../../pull/5547))
* Fix deepcopy and pickle for classes derived from `np.ndarray` ([#5536](../../pull/5536))
* owheatmap: Fix assertion error when restoring selection ([#5517](../../pull/5517))
* Pivot: Handle empty data, metas only ([#5527](../../pull/5527))
* table_to_frame - handle numeric columns with dtype=object ([#5474](../../pull/5474))
* listfilter: Bypass QListView.dropEvent ([#5477](../../pull/5477))


[3.29.3] - 2021-06-09
--------------------
##### Bugfixes
* Create Class: fix incorrect value assignment


[3.29.2] - 2021-06-08
--------------------
##### Bugfixes
* Bump orange-canvas-core minimum version requirement ([#5472](../../pull/5472))
* owpca: fix component selection when dragging selection line ([#5469](../../pull/5469))
* Save File when workflow basedir is an empty string ([#5459](../../pull/5459))


[3.29.1] - 2021-05-31
--------------------


[3.29.0] - 2021-05-28
--------------------
##### Enhancements
* Scatter plot: Bring discrete attributes functionality back ([#5440](../../pull/5440))
* Caching data delegate ([#5296](../../pull/5296))
* textimport: Mark encoding errors in the preview ([#5438](../../pull/5438))
* DBSCAN: Optional normalization ([#5428](../../pull/5428))
* Automated and better summaries ([#5308](../../pull/5308))
* Transpose: Offload work onto separate thread, remove redundant instance ([#5314](../../pull/5314))
* Domain transformations in batches for less memory use ([#5218](../../pull/5218))
* Feature Statistics: Add median ([#5325](../../pull/5325))
* New widget: Aggregate Columns ([#5256](../../pull/5256))

##### Bugfixes
* Outlier detection: keep instance ids, make thread safe ([#5427](../../pull/5427))
* UrlReader: Support urls with special characters ([#5412](../../pull/5412))
* Speed-up slow table_to_frame ([#5413](../../pull/5413))
* Line Plot: Single instance input fix ([#5408](../../pull/5408))
* Pivot: Assign dataset name to output tables ([#5404](../../pull/5404))
* setup.py: Exclude benchmark directory from the install ([#5392](../../pull/5392))
* Errors when converting negative timestamps on Windows ([#5388](../../pull/5388))
* Nomogram: Retain original compute_value ([#5382](../../pull/5382))
* Radviz VizRank: Implement on_selection_changed ([#5338](../../pull/5338))


[3.28.0] - 2021-03-05
--------------------
##### Enhancements
* Bar Plot: Improve "Group by" visualization ([#5301](../../pull/5301))
* Violin Plot: New widget ([#5252](../../pull/5252))
* Test and Score: Copy selected rows to clipboard ([#5203](../../pull/5203))
* Projections: Allow transparent subset ([#5141](../../pull/5141))
* Gradient Boosting: New widget ([#5160](../../pull/5160))
* Impute: Allow setting a default value for all numeric and time variables ([#5102](../../pull/5102))
* Distribution: Show equal bar widths on unique-valued bins ([#5139](../../pull/5139))
* Implement proper Lift curve; keep Cumulative gains as an option ([#5075](../../pull/5075))
* Create Instance: New widget ([#5033](../../pull/5033))
* Add add_column and concatenate methods to Table ([#5251](../../pull/5251))

##### Bugfixes
* Calibration model: Work with numpy data ([#5159](../../pull/5159))
* Rank: Switch to manual selection on deselect ([#5271](../../pull/5271))
* Create Class: multiple patterns for a class value ([#5283](../../pull/5283))
* Test and Score: Fix stratification warnings ([#5281](../../pull/5281))
* Predictions: Fix crash when clicking on empty left area ([#5222](../../pull/5222))
* Distribution: vectorize variance, speeds up normalization ([#5230](../../pull/5230))
* owimpute: Make `default_numeric` locale independant ([#5209](../../pull/5209))
* Pivot: Display time variable in time format ([#5212](../../pull/5212))
* OWScatterPlotBase: Ignore 'Other' when showing color regions ([#5214](../../pull/5214))
* Fix performance regression in scatterplot ([#5206](../../pull/5206))
* Pivot: Fix table for categorical variables ([#5193](../../pull/5193))
* Distance Matrix: Fix freeze with large selections ([#5176](../../pull/5176))
* Xls reader: Error as missing value ([#5192](../../pull/5192))
* owdataset: Do not capture self in closure ([#5198](../../pull/5198))
* ROC shows all points, including the last ([#5138](../../pull/5138))
* Pivot: Output date for group by table ([#5202](../../pull/5202))
* Enable classification tests ([#5168](../../pull/5168))
* Data Table: Fix freeze with large selections ([#5164](../../pull/5164))
* OWImpute: Preserve default method setting ([#5181](../../pull/5181))
* AdaBoost: Set maximum number of estimators to 10000 ([#5165](../../pull/5165))
* Feature Statistics: Error in time variable display ([#5152](../../pull/5152))
* owpythonscript: Use signal id as is ([#5147](../../pull/5147))
* SqlTable use empty string instead of None for StringVariable ([#5120](../../pull/5120))
* Pivot: Handle big dataset ([#5104](../../pull/5104))
* SQL: Fix the issue with database collation setting when retrieving column values ([#5089](../../pull/5089))
* impute: Remove class vars from input data for ReplaceUnknownsModel ([#5083](../../pull/5083))
* Edit Domain: Preserve renames on categories merge ([#5072](../../pull/5072))
* CSV File Import: sort discrete values naturally ([#5041](../../pull/5041))


[3.27.1] - 2020-10-23
--------------------
##### Bugfixes
* customizableplot.available_font_families: Fix for non-existent default family ([#5037](../../pull/5037))
* Raise canvas-core version to fix some problems with Qt 5.9 ([#5045](../../pull/5045))


[3.27.0] - 2020-10-09
--------------------
##### Enhancements
* Table: Re-add info box about data properties ([#5011](../../pull/5011))
* Rank widget computation in a separate thread ([#4908](../../pull/4908))
* Neighbors: improve exclusion of references, checkbox to (un)limit output data ([#4997](../../pull/4997))
* Bar Plot: New widget ([#4923](../../pull/4923))
* Replace listViews with searchable listViews ([#4924](../../pull/4924))
* Add an option to set intercept in linear regression to 0 ([#4958](../../pull/4958))
* CSV File Import: Add support for explicit workflow relative paths ([#4872](../../pull/4872))
* OWColor: Saving and loading color schemata ([#4977](../../pull/4977))
* Distributions: Add sorting by category size ([#4959](../../pull/4959))
* Proxy support for embeddings ([#4953](../../pull/4953))
* Discretize: Manual cut points entry ([#4929](../../pull/4929))
* Predictions: Allow selecting a subset of rows ([#4871](../../pull/4871))
* Projection plots: Customize labels ([#4828](../../pull/4828))
* get_unique_names: Handle more independent names ([#4866](../../pull/4866))
* Edit Domain: Add option to unlink variable from source variable ([#4863](../../pull/4863))
* Louvain Clustering: Add cosine similarity ([#4864](../../pull/4864))

##### Bugfixes
* Fix line plot's send_report ([#5018](../../pull/5018))
* Scatter Plot: fix unzoom with density plot ([#5004](../../pull/5004))
* ownomogram: Fix wrapped C++ obj error ([#5005](../../pull/5005))
* Fix slicing in from_table ([#4963](../../pull/4963))
* Edit Domain: Multiple item rename/merge ([#4949](../../pull/4949))
* ProjectionWidgetMixinTests: set shape attribute only when discrete var available ([#4946](../../pull/4946))
* Fix variables equality and hashes  ([#4957](../../pull/4957))
* Fix wrong assert in heatmap ([#4955](../../pull/4955))
* Edit Domain (and perhaps other widgets) could cause missing data later in the workflow ([#4922](../../pull/4922))
* OWScatterPlotBase: Reset view before labels update ([#4907](../../pull/4907))
* owcsvimport: Fix a type error in _open for zip archive ([#4921](../../pull/4921))
* Line Plot: Reset axis ticks on data change ([#4873](../../pull/4873))
* MDS: Move lines when points are jittered ([#4920](../../pull/4920))
* owcsvimport: Handle decimal and thousand separator ([#4915](../../pull/4915))
* Select Rows: fix fail when time variable in metas ([#4912](../../pull/4912))
* Re-add TupleList to DiscreteVariable ([#4879](../../pull/4879))
* concurrent: Use a workaround for QObject wrapper deletion ([#4635](../../pull/4635))
* normalize: Adjust number_of_decimals after scaling ([#4779](../../pull/4779))


[3.26.0] - 2020-06-12
--------------------
##### Enhancements
* main: Log to main window output view ([#4842](../../pull/4842))
* Feature statistics report ([#4812](../../pull/4812))
* Distributions: change Histogram Data output ([#4832](../../pull/4832))
* owtable: output sorted data ([#4644](../../pull/4644))
* Add an option to Concatenate to merge columns with different formulae ([#4831](../../pull/4831))
* CSV Import: guess data types ([#4838](../../pull/4838))
* HeatMap: Allow setting the center when using diverging palettes ([#4809](../../pull/4809))
* Heatmap: Split columns, Column annotations ([#4703](../../pull/4703))
* Sort values naturally when reading files ([#4793](../../pull/4793))
* Color widget: Add reset button ([#4718](../../pull/4718))
* Gradient selection/parameters widget ([#4596](../../pull/4596))
* Select Rows: Allow partial context matches ([#4740](../../pull/4740))
* Edit Domain: Add an option to change the output table name ([#4722](../../pull/4722))
* ApplyDomain: data info displayed in the status bar ([#4611](../../pull/4611))
* BoxPlot: data info displayed in the status bar ([#4626](../../pull/4626))
* LinePlot: data info displayed in the status bar ([#4633](../../pull/4633))
* MosaicDisplay: data info displayed in the status bar ([#4634](../../pull/4634))
* CreateClass: data info displayed in the status bar ([#4625](../../pull/4625))

##### Bugfixes
* Variable: Fix cases when equal variables had different hashes ([#4843](../../pull/4843))
* OWBoxPlot: Fix wrong labels position and ordering for values with no items ([#4829](../../pull/4829))
* Select Rows: Fix saving meta variables in context ([#4830](../../pull/4830))
* Select Columns: Fix attributes sorting ([#4827](../../pull/4827))
* Fix and update Softmax regression learner ([#4767](../../pull/4767))
* PCA: fix subsets with the "Data" output ([#4811](../../pull/4811))
* OWContinuize: Fix treatment of continuous features. ([#4806](../../pull/4806))
* Select Rows: Fix incorrectly stored values in settings ([#4798](../../pull/4798))
* Fix colors for discrete variables with >256 values ([#4803](../../pull/4803))
* Unique domain checks ([#4760](../../pull/4760))
* owheatmap: Use 'is' instead of 'eq' for column id comparison ([#4788](../../pull/4788))
* BoxPlot: Fix invalid data range ([#4769](../../pull/4769))
* graphicstextlist: Fix size/spacing adjustment for single item ([#4777](../../pull/4777))
* Feature Statistics: Fix wrong or even crashing selections ([#4741](../../pull/4741))
* UrlReader: shorten TempFile extension ([#4747](../../pull/4747))
* Embedder: catch machine id setting type error ([#4675](../../pull/4675))
* relief: Fix contingency (de)allocation ([#4745](../../pull/4745))
* Test and Score: Improve data errors ([#4738](../../pull/4738))
* PythagoreanTree: Fix crushing when number of classes decreases ([#4743](../../pull/4743))
* Fix report in Predictions ([#4709](../../pull/4709))
* owheatmap: Handle all N/A column for color annotation ([#4742](../../pull/4742))
* Distributions widget's legend: Remove the square from sigma in normal and Rayleigh ([#4739](../../pull/4739))
* Several fixes in learners/models ([#4655](../../pull/4655))
* Heatmap: Split by missing values ([#4686](../../pull/4686))
* Owpaintdata, owpivot: ensure unique domain ([#4578](../../pull/4578))
* Pythagorantrees/forests: change context handler ([#4656](../../pull/4656))
* Color: Fix renaming of variables ([#4669](../../pull/4669))
* Heatmap: Sticky footer ([#4610](../../pull/4610))
* SOM: fix colors for numeric variables ([#4660](../../pull/4660))
* Fixes for deprecations in 3.26, and for changed behaviour of file dialog ([#4643](../../pull/4643))


[3.25.1] - 2020-05-22
--------------------
##### Bugfixes
* Fix compatibility with scikit-learn 0.23 ([#4768](../../pull/4768))


[3.25.0] - 2020-04-10
--------------------
##### Enhancements
* Searchable combo boxes in all evaluate widgets ([#4564](../../pull/4564))
* Searchable combo boxes in all visualize widgets ([#4563](../../pull/4563))
* Searchable combo boxes in all unsupervised widgets ([#4565](../../pull/4565))
* Projections keep colors after merging values ([#4577](../../pull/4577))
* Distance Map: Add a color map legend ([#4593](../../pull/4593))
* Orange.misc.environ config ([#4576](../../pull/4576))
* Searchable combo boxes in all data widgets ([#4562](../../pull/4562))
* Fix printing values with too few decimals ([#4575](../../pull/4575))
* Scatter Plot: Replace combo box with search combo box ([#4447](../../pull/4447))
* Save widgets: Store paths relative to workflow directory ([#4532](../../pull/4532))
* Edit Domain: Option to merge less frequent values ([#4477](../../pull/4477))
* Load Model: Use paths relative to workflow file ([#4534](../../pull/4534))
* Ignore missing values in density plots of projections ([#4525](../../pull/4525))
* Impose a sensible z-order to points in projections  ([#4504](../../pull/4504))
* Use Github actions as a new CI system. ([#4482](../../pull/4482))
* Testing with Tox ([#4481](../../pull/4481))
* Add row side color annotations ([#4443](../../pull/4443))
* Outliers: Offload work onto separate thread ([#4412](../../pull/4412))
* OWScatterPlot: axis displays time specific labels for time variable ([#4434](../../pull/4434))
* Predictions: Update splitter on resize ([#4433](../../pull/4433))
* Import openTSNE lazily for faster loading of Orange ([#4424](../../pull/4424))
* Silhouette Plot: Always output scores ([#4423](../../pull/4423))
* Heatmap: Tighter layout ([#4390](../../pull/4390))
* Reorganize continuous palettes ([#4305](../../pull/4305))
* Outliers: Save model into compute_value ([#4372](../../pull/4372))
* Allow concurrent transformation of tables into new domains ([#4363](../../pull/4363))
* Test & Score: Add comparison of models ([#4261](../../pull/4261))
* Outliers: Widget upgrade ([#4338](../../pull/4338))
* Concatenate: data info displayed in the status bar ([#4617](../../pull/4617))
* Distributions: data info displayed in the status bar ([#4627](../../pull/4627))
* MergeData: data info displayed in the status bar ([#4592](../../pull/4592))
* Neighbors: data info displayed in the status bar ([#4612](../../pull/4612))
* SelectByDataIndex: data info displayed in the status bar ([#4595](../../pull/4595))
* OWOutliers: Data info displayed in the status bar ([#4547](../../pull/4547))
* OWDataProjectionWidget: Upgrade status bar data info ([#4544](../../pull/4544))
* Heatmap: Restore ability to cluster larger datasets ([#4290](../../pull/4290))
* OWImpute: Data info displayed in the status bar ([#4499](../../pull/4499))
* OWFile: Data info displayed in the status bar ([#4506](../../pull/4506))
* OWCSVImport: Data info displayed in the status bar ([#4509](../../pull/4509))
* OWDatasets: Data info displayed in the status bar ([#4512](../../pull/4512))
* OWDataInfo: Data info displayed in the status bar ([#4513](../../pull/4513))
* OWSave: Data info displayed in the status bar ([#4505](../../pull/4505))
* OWPurgeDomain: Data info displayed in the status bar ([#4502](../../pull/4502))
* OWColor: Data info displayed in the status bar ([#4501](../../pull/4501))
* OWRandomize: Data info displayed in the status bar ([#4498](../../pull/4498))
* OWPivotTable: Data info displayed in the status bar ([#4472](../../pull/4472))
* OWContinuize: Data info displayed in the status bar ([#4494](../../pull/4494))
* OWFeatureConstructor: Data info displayed in the status bar ([#4496](../../pull/4496))
* OWSelectRows: Data info displayed in the status bar ([#4471](../../pull/4471))
* OWDiscretize: Data info displayed in the status bar ([#4495](../../pull/4495))
* OWDataSampler: Data info displayed in the status bar ([#4492](../../pull/4492))
* OWRank: Data info displayed in the status bar ([#4473](../../pull/4473))
* OWContinuize: Provide the same options as in Preprocess/Normalize ([#4466](../../pull/4466))
* OWCorrelations: Data info displayed in the status bar ([#4455](../../pull/4455))
* OWEditDomain: Data info displayed in the status bar ([#4456](../../pull/4456))
* OWSelectColumns: Data info displayed in the status bar ([#4454](../../pull/4454))
* OWTranspose: Data info displayed in the status bar ([#4413](../../pull/4413))
* OWFeatureStatistics: data info displayed in the status bar ([#4409](../../pull/4409))
* OWPreporcess: Data info displayed in the status bar ([#4414](../../pull/4414))

##### Bugfixes
* Edit Domain: Fix merge values when missing data ([#4636](../../pull/4636))
* Table: Send correct output when switching between tabs ([#4619](../../pull/4619))
* Give created QGraphicsScenes a parent ([#4352](../../pull/4352))
* Fix dimensionality of probabilities from values ([#4629](../../pull/4629))
* PyTableModel: Allow wrapping empty lists ([#4631](../../pull/4631))
* Edit Domain: Improve Text/Categorical to Time conversion ([#4601](../../pull/4601))
* colorpalettes: fix BinnedContinuousPalette color assignments ([#4609](../../pull/4609))
* Classification models output correct shapes ([#4602](../../pull/4602))
* owtestandscore: Add cancelled tasks to dispose queue ([#4615](../../pull/4615))
* Nomogram: Fix crash on Python 3.8 ([#4591](../../pull/4591))
* K-means slowness ([#4541](../../pull/4541))
* Use new access token in cleanup workflow ([#4590](../../pull/4590))
* Detect duplicate names of variables in projections. ([#4550](../../pull/4550))
* Paint Data: Send correct output after clicking Reset to Input ([#4551](../../pull/4551))
* TimeVariable.parse: Do not modify _ISO_FORMATS ([#4539](../../pull/4539))
* heatmap: Ensure minimim size for color annotations ([#4519](../../pull/4519))
* OWEditDomain: Clear editor when data is disconnected ([#4484](../../pull/4484))
* ContinuousPalettesModel: Disable 'category' items via `flags` ([#4538](../../pull/4538))
* utils/image: Return early when image is empty ([#4520](../../pull/4520))
* graphicstextlist: Use integer font metrics again ([#4524](../../pull/4524))
* Concatenate: Fix wrong merging of categorical features ([#4425](../../pull/4425))
* Ensure unique var names in file ([#4431](../../pull/4431))
* Rank: Fix error with Random forest ([#4457](../../pull/4457))
* owhierclustering: Update the scene's sceneRect ([#4459](../../pull/4459))
* Venn Diagram is slow for big datasets ([#4400](../../pull/4400))
* Fix missing values after purging unused values ([#4432](../../pull/4432))
* File: Construct unique column names. ([#4420](../../pull/4420))
* format_summary_details: Replace 'features' with 'variables' ([#4419](../../pull/4419))
* Feature Constructor: Catch exceptions ([#4401](../../pull/4401))
* Continuize: Disable normalizing sparse data ([#4379](../../pull/4379))
* Python script serialization state ([#4345](../../pull/4345))
* DataProjectionWidget: Update combos on new data ([#4405](../../pull/4405))
* DataProjectionWidget: attribute Selected ([#4393](../../pull/4393))
* Fix slow clear/delete in 'Heat Map' 'Hier. Clustering', 'Distance Map' ([#4365](../../pull/4365))
* Explicitly define the protocol version for pickling ([#4388](../../pull/4388))
* Table.from_table: fix caching with reused ids  ([#4370](../../pull/4370))
* FeatureStatistics: Convert selected rows to list ([#4375](../../pull/4375))
* Error message on tSNE with one variable ([#4364](../../pull/4364))
* Normalize: Integer variable representation ([#4350](../../pull/4350))
* DendrogramWidget: Remove event filters before removing items ([#4361](../../pull/4361))
* Round bhattacharayya ([#4340](../../pull/4340))


[3.24.1] - 2020-01-17
--------------------
##### Enhancements
* OWPreprocess: data info displayed in status bar ([#4333](../../pull/4333))
* OWDiscretize: data info displayed in status bar ([#4331](../../pull/4331))
* OWContinuize: data info displayed in status bar ([#4327](../../pull/4327))
* OWTranspose: data info displayed in status bar ([#4295](../../pull/4295))
* Silhouette plot: Accept distance matrix on input ([#4313](../../pull/4313))
* Preprocess: Add filtering by missing values ([#4266](../../pull/4266))
* Allow add-ons to register file format for the Save widget ([#4302](../../pull/4302))
* Box Plot: Add box for missing group values ([#4292](../../pull/4292))
* clustering/hierarchical: Use optimal\_leaf\_ordering from scipy ([#4288](../../pull/4288))
* Distributions: Add option to hide bars ([#4301](../../pull/4301))

##### Bugfixes
* ExcelReader: Speedup ([#4339](../../pull/4339))
* Nomogram: Adjust scale considering label width ([#4329](../../pull/4329))
* owhierarchicalclustering: Prescale dendrogram geometry ([#4322](../../pull/4322))
* table\_to\_frame: metas lost on conversion ([#4259](../../pull/4259))
* TestOWRank: Setting type should not be np.ndarray ([#4315](../../pull/4315))
* Distances: Fix restoring the cosine distance ([#4311](../../pull/4311))
* Fix saving workflows that contain a Rank widget ([#4289](../../pull/4289))
* utils/textimport: Remove 'exclusive' kwarg from QActionGroup call ([#4298](../../pull/4298))


[3.24.0] - 2019-12-20
--------------------
##### Enhancements
* Remove Variable.make ([#3925](../../pull/3925))
* OWTreeViewer: Bold predicted values in tree nodes ([#4269](../../pull/4269))
* Edit Domain: Reinterpret column type transforms ([#4262](../../pull/4262))
* OWTree: Add 'Classification Tree' keyword ([#4283](../../pull/4283))
* ConcurrentWidgetMixin: Cancel task on input change ([#4219](../../pull/4219))
* PCA: Add a signal with original data + pca ([#4255](../../pull/4255))
* Distances: Offload work to a separate thread ([#4046](../../pull/4046))
* Venn Diagram: Add relations over columns, simplify over rows ([#4006](../../pull/4006))
* Merge Data: Implement context settings ([#4248](../../pull/4248))
* Heatmap: Option to center color palette at 0 ([#4218](../../pull/4218))
* owheatmap: Add Split By combo box ([#4234](../../pull/4234))
* Heatmap: Allow labeling by any variable ([#4209](../../pull/4209))
* Test & Score Widget: Cancellability on input change. ([#4079](../../pull/4079))
* Boxplot no longer stretches bars when this is uninformative ([#4176](../../pull/4176))
* Box plot: Add 'order by importance' checkbox to groups ([#4055](../../pull/4055))
* Add pyproject.toml ([#4179](../../pull/4179))
* Add remove sparse features preprocessor ([#4093](../../pull/4093))
* Output to OWCorrespondence ([#4180](../../pull/4180))
* macOS: Installer python version ([#4130](../../pull/4130))
* OWBoxPlot: Show missing values ([#4135](../../pull/4135))
* Neighbors: Data info displayed in status bar ([#4157](../../pull/4157))
* Preprocess: Tests update
* Bhatthacharayya distance ([#4111](../../pull/4111))
* MergeData: Don't remove duplicate columns with different data ([#4100](../../pull/4100))
* Nice binning of time variables (Distributions, SOM) ([#4123](../../pull/4123))
* OWBaseSql: Base widget for connecting to a database ([#4083](../../pull/4083))
* k-Means: Add normalization checkbox ([#4099](../../pull/4099))
* Datasets: Remove control area ([#4071](../../pull/4071))
* Correlations: Include continuous class and meta variables ([#4067](../../pull/4067))

##### Bugfixes
* SQL: Save selected backend to settings ([#4270](../../pull/4270))
* ExcelReader: Migrate to openpyxl ([#4279](../../pull/4279))
* owcsvimport: Fix last/recent item serialization ([#4272](../../pull/4272))
* owselectcolumns: Fix move up/down type error ([#4271](../../pull/4271))
* Table: Keep pending selection if data is None ([#4281](../../pull/4281))
* Select Rows: Fix crash on changed variable type ([#4254](../../pull/4254))
* Louvain Clustering: Update graph output for compatibility with the new network add-on ([#4258](../../pull/4258))
* Merge data: Migrate settings ([#4263](../../pull/4263))
* Feature Statistics: Do not crash on empty domain ([#4245](../../pull/4245))
* File widget: fix name change ([#4235](../../pull/4235))
* Various fixes of box plot ([#4231](../../pull/4231))
* Fix guessing strategy for date and time variables ([#4226](../../pull/4226))
* owtable: Preserve tab order of updated inputs ([#4225](../../pull/4225))
* Feature Constructor: Compatibility with Python 3.8 ([#4222](../../pull/4222))
* File: Fix domain edit on no changes ([#4232](../../pull/4232))
* Normalizer: Retain attributes of attributes ([#4217](../../pull/4217))
* Hierarchical Clustering: Fix Annotations selection ([#4214](../../pull/4214))
* MDS: Place lines under points and labels ([#4213](../../pull/4213))
* Fix crash in Preprocess' Select Relevant Features when there are no features ([#4207](../../pull/4207))
* Louvain Clustering: fix setting to restore correctly ([#4187](../../pull/4187))
* SOM: Fix crash when clicking on empty canvas ([#4177](../../pull/4177))
* build-conda-installer.sh: Do not use python version from Miniconda ([#4142](../../pull/4142))
* Fix core dump in CI in SOM on sparse data ([#4174](../../pull/4174))
* ROC and LiftCurve: Store context settings ([#4138](../../pull/4138))
* Preprocess: Tests update
* Normalize: Fix crash with nan column ([#4125](../../pull/4125))
* Warning for discrete variable with >100 values in OWFile ([#4120](../../pull/4120))
* owcsvimport: Make the progress update a direct connection ([#4109](../../pull/4109))
* Scatterplot: Disable vizrank when features on input ([#4102](../../pull/4102))
* Sieve: Disable vizrank when features on input ([#4101](../../pull/4101))
* Datetime conversion to string ([#4098](../../pull/4098))
* Self-Organizing Map: Fix restoring width ([#4097](../../pull/4097))
* K-means: Save Silhouette Scores selection ([#4082](../../pull/4082))
* OWTestLearners: Vertically align to center ([#4095](../../pull/4095))
* Model's data_to_model_domain supports sparse matrix ([#4081](../../pull/4081))
* Color widget failed after removing Variable.make ([#4041](../../pull/4041))
* Deprecated processEvents argument in Test&Score ([#4047](../../pull/4047))
* Merge data: Rename variables with duplicated names ([#4076](../../pull/4076))
* k-Means: Impute missing values ([#4068](../../pull/4068))
* Predictions: Handle discrete target with no values ([#4066](../../pull/4066))
* Fix storing and retrieving selection, and unconditional auto commit ([#3957](../../pull/3957))
* Scatterplot: Enable/disable vizrank button ([#4016](../../pull/4016))
* Correlations: Add progress bar, retain responsiveness ([#4011](../../pull/4011))


[3.23.1] - 2019-10-03
--------------------
##### Bugfixes
* Addons install with conda by default


[3.23.0] - 2019-09-05
--------------------
##### Enhancements
* Pull YAML feed of notifications on startup, refactor notifications ([#3933](../../pull/3933))
* Widget for Self-Organizing Maps ([#3928](../../pull/3928))
* DBSCAN widget ([#3917](../../pull/3917))
* BoxPlot: Write the anova/t-test statistic onto the plot. ([#3945](../../pull/3945))
* Feature Constructor: Easier categorical features, allow creation of date/time, easier use of string data ([#3936](../../pull/3936))
* Merge data allows matching by multiple pairs of columns ([#3919](../../pull/3919))
* Sticky graphics header/footer views ([#3930](../../pull/3930))
* Shiny renewed widget Distributions ([#3896](../../pull/3896))
* Calibration plot (add performance curves) and a new Calibrated Learner widget ([#3881](../../pull/3881))
* Added Specificity as a new score in Test&Score ([#3907](../../pull/3907))
* Separate canvas and base widget ([#3772](../../pull/3772))

##### Bugfixes
* Conda Installer: Restore compatibility with latest anaconda python ([#4004](../../pull/4004))
* Scatter plot: Hidden variables fix ([#3985](../../pull/3985))
* Boxplot fix ([#3983](../../pull/3983))
* Heat map: Cannot cluster a single instance ([#3980](../../pull/3980))
* Test and Score: Sort numerically, not alphabetically ([#3951](../../pull/3951))
* OWProjectionWidgetBase: Update when domain is changed ([#3946](../../pull/3946))
* Change normalization to Scaling in SVM ([#3898](../../pull/3898))
* Restore usage tracking ([#3921](../../pull/3921))
* main: Fix widget settings directory path to clear ([#3932](../../pull/3932))
* Backcompatibility stubs ([#3926](../../pull/3926))
* OWNeighbours fix manual apply for some options ([#3911](../../pull/3911))
* Documentation links ([#3897](../../pull/3897))
* Update output on new input, even if auto commit is disabled ([#3844](../../pull/3844))


[3.22.0] - 2019-06-26
--------------------
##### Enhancements
* Unified clustering API ([#3814](../../pull/3814))
* t-SNE: Load openTSNE lazily" ([#3894](../../pull/3894))
* Replace popups with non-intrusive notifications ([#3855](../../pull/3855))
* CSV File Import widget ([#3876](../../pull/3876))
* t-SNE: Load openTSNE lazily ([#3883](../../pull/3883))
* Faster drawing in scatterplot ([#3871](../../pull/3871))
* Mosaic: Wrap legend ([#3866](../../pull/3866))
* Add conditions for all variables, or all numeric or textual variables ([#3836](../../pull/3836))
* Shared namespaces for PythonScript widgets ([#3840](../../pull/3840))
* Pivot: New widget ([#3823](../../pull/3823))
* Reset settings button ([#3795](../../pull/3795))
* WebviewWidget: expose JavaScript timeout limit ([#3811](../../pull/3811))

##### Bugfixes
* OWLinearProjection: limit axis ([#3885](../../pull/3885))
* Development readme ([#3889](../../pull/3889))
* OWRadviz: limit number of vars in RadvizVizRank ([#3886](../../pull/3886))
* OWCreateClass: Reuse variables created with same settings ([#3868](../../pull/3868))
* DistMatrix should return numpy datatypes ([#3865](../../pull/3865))
* OWRadviz: legible axis labels ([#3809](../../pull/3809))
* util: Fix bincount for object arrays ([#3831](../../pull/3831))
* condainstall: Fix env for running conda ([#3843](../../pull/3843))
* DBSCAN: Fix predicted labels ([#3833](../../pull/3833))
* Rare quickmenu crash fix ([#3832](../../pull/3832))
* OWWidget: Fix an error on mouse press when widget has no basic layout ([#3821](../../pull/3821))
* VerticalItemDelegate: Do not cut long vertical labels ([#3803](../../pull/3803))
* Minor improvements to pythagorean trees ([#3777](../../pull/3777))
* owmds: Fix test for error display/activate ([#3813](../../pull/3813))


[3.21.0] - 2019-05-20
---------------------
##### Enhancements
* Error Animations ([#3788](../../pull/3788))
* OWTSNE: Offload computation to separate thread ([#3604](../../pull/3604))
* OWDistributions: add cumulative distribution ([#3766](../../pull/3766))
* Edit Domain: Merge categorical values ([#3725](../../pull/3725))
* Transform: Values of primitive variables as feature names ([#3721](../../pull/3721))
* k-Means: Output centroid labels ([#3695](../../pull/3695))
* Support sparse Jaccard ([#3657](../../pull/3657))
* Offload work to a separate thread ([#3627](../../pull/3627))
* Correlations: Enhancements and fixes ([#3660](../../pull/3660))
* DomainEditor: Indicate changed variables with bold face font ([#3576](../../pull/3576))
* OWPythonScript: dropping and pasting of python scripts ([#3611](../../pull/3611))
* Improve Save widget's gui ([#3545](../../pull/3545))
* Table.from_numpy: Replace infs with nans. ([#3624](../../pull/3624))
* FDR: Calculate FDR using numpy ([#3625](../../pull/3625))
* Correlations: fixes and enhancements ([#3591](../../pull/3591))
* Preprocess: implement Select Relevant Feature's percentile ([#3588](../../pull/3588))
* Added keyboard shortcuts for Align & Freeze/Unfreeze  ([#3601](../../pull/3601))
* PCA: Remove SVD & add normalization for sparse ([#3581](../../pull/3581))
* OwLouvain: Add normalize data checkbox to PCA preprocessing ([#3573](../../pull/3573))
* Use %g (including sci notation) if number of decimals is not set ([#3574](../../pull/3574))
* PCA: Preserve f32s & reduce memory footprint when computing means ([#3582](../../pull/3582))

##### Bugfixes
* OWLinePlot: legible bottom axis labels ([#3768](../../pull/3768))
* OWPythagorasTree: Enable node selection from forests with categorical variables. ([#3775](../../pull/3775))
* canvas/help: Fix a NameError in exception handler ([#3759](../../pull/3759))
* ProjectionWidgetTestMixin: Fix test_plot_once ([#3738](../../pull/3738))
* LinearProjectionVizRank: Add a necessary check ([#3732](../../pull/3732))
* stats: Fix statistics for primitive variables ([#3722](../../pull/3722))
* SQL Table: Restore selected table from settings ([#3703](../../pull/3703))
* OWAnchorProjectionWidget: Retain valid_data when reloading dataset ([#3718](../../pull/3718))
* Compatibility with Logitech's Smart Move ([#3702](../../pull/3702))
* OWTable: Don't set selection when there is no data on input ([#3693](../../pull/3693))
* Save Data: Remove extra file extensions ([#3700](../../pull/3700))
* owheatmap: Fix group label size policy ([#3688](../../pull/3688))
* Test & Score: hide warnings for hidden scores. ([#3676](../../pull/3676))
* VizRankDialog: Use extended thread pool to prevent segfaults ([#3669](../../pull/3669))
* Send usage statistics in a thread at startup ([#3632](../../pull/3632))
* Transform: Replace 'Preprocess' input with 'Template Data' input ([#3673](../../pull/3673))
* OWTable: Include attributes from meta attributes in header row labels ([#3633](../../pull/3633))
* OWSieve: Fix operators in tooltips ([#3602](../../pull/3602))
* MDS: Handle subset data ([#3620](../../pull/3620))
* migrate_context: widgets crash when migrating context without version ([#3603](../../pull/3603))
* Louvain clustering fails when Table or ndarray on input. ([#3618](../../pull/3618))
* Orange: add more informatino on missing compiled libraries ([#3614](../../pull/3614))
* LinearProjection: Disable LDA for less than three classes ([#3615](../../pull/3615))
* Orange restart dialogs: Improve wording ([#3605](../../pull/3605))
* OWDistances: Use only binary features for Jaccard distance ([#3569](../../pull/3569))
* ScatterplotGraph: Use mapToView instead of mapRectFromParent ([#3571](../../pull/3571))
* MDS: Fix crashes when feeding column distances ([#3583](../../pull/3583))
* Naive Bayes: Ignore existing classes in Laplacian smoothing ([#3575](../../pull/3575))


[3.20.1] - 2019-02-12
--------------------
##### Enhancements
* t-SNE: Add Normalize data checkbox ([#3570](../../pull/3570))
* Louvain show number of clusters ([#3572](../../pull/3572))

##### Bugfixes
* t-SNE speed-ups ([#3592](../../pull/3592))
* setup.py: Specify python-louvain version constraint ([#3587](../../pull/3587))


[3.20.0] - 2019-02-01
--------------------
##### Enhancements
* Naive Bayes: Implement predict, fix predict_storage ([#3540](../../pull/3540))
* OWTransform: add option to keep original data #3526 ([#3549](../../pull/3549))
* Implement better randomized PCA ([#3532](../../pull/3532))
* Scatterplot: Draw separate regression lines for colors; add orthonormal regression ([#3518](../../pull/3518))
* Edit Domain: Add support for ordered categorical variables ([#3535](../../pull/3535))
* Warn when data subset is unrelated to data ([#3507](../../pull/3507))
* Label subset ([#3506](../../pull/3506))
* OWLouvain: Ensure deterministic clustering ([#3492](../../pull/3492))
* t-SNE: Updates 2. ([#3475](../../pull/3475))
* Save Data: Support saving to Excel ([#3453](../../pull/3453))
* OWLinePlot: Move from prototypes to core ([#3440](../../pull/3440))

##### Bugfixes
* Warning about a window title without a placeholder ([#3554](../../pull/3554))
* PaintData: Fix control area width ([#3560](../../pull/3560))
* OWWidget: Remove wheelEvent reimplementation ([#3557](../../pull/3557))
* Data sets: Remove class from 'bupa' ([#3556](../../pull/3556))
* Neighbours: Show error when data and reference have different domain. ([#3547](../../pull/3547))
* Feature statistics fixes ([#3480](../../pull/3480))
* Removed header types and flags from .csv and .tab ([#3427](../../pull/3427))
* Python script widget: prevent data loss  ([#3529](../../pull/3529))
* setup.cfg: Change the default for with-htmlhelp config option ([#3536](../../pull/3536))
* Don't add miniconda to path and register as system python ([#3525](../../pull/3525))
* OWDataProjectionWidget: check validity, fix sparse data reloading ([#3485](../../pull/3485))
* OWMergeData removes table meta attributes ([#3474](../../pull/3474))
* [MNT] Remove segfault in tests for building documentation ([#3491](../../pull/3491))
* OWFeatureStatistics: Fix scipy.stats.mode crash on sparse data ([#3488](../../pull/3488))
* OWDataProjectionWidget: Fix applying selection ([#3466](../../pull/3466))


[3.19.0] - 2018-12-11
--------------------
##### Enhancements
* Remove discrete attributes from scatter plot's axes ([#3434](../../pull/3434))
* OWScatterPlotBase: Animate dot resize ([#3436](../../pull/3436))
* Introduce stacking ([#3291](../../pull/3291))
* OWWidget: Input/output summary  ([#2556](../../pull/2556))
* File: Provide percent missing values in Info box ([#3305](../../pull/3305))
* OWHierarchicalClustering: Use selection indices for selection restore ([#3282](../../pull/3282))
* Data Info display data set name ([#3187](../../pull/3187))
* tSNE: Output preprocessor ([#3407](../../pull/3407))
* Pythagorean Tree: children order ([#3393](../../pull/3393))

##### Bugfixes
* RemoveNaNColumns: Move to preprocess ([#3464](../../pull/3464))
* Scatterplot Vizrank: Don't use discrete variables ([#3463](../../pull/3463))
* canvas/widgetsscheme: Remove use of QObject.destroyed for logging ([#3447](../../pull/3447))
* OWFeatureStatistics: Don't attempt to sort when no data on input ([#3449](../../pull/3449))
* Rank: Fix crash on dataset with missing values ([#3458](../../pull/3458))
* Radviz: Enable projection for less than two selected variables ([#3444](../../pull/3444))
* t-SNE: Generate temporary projection data ([#3454](../../pull/3454))
* Scatter Plot: Always setup plot ([#3450](../../pull/3450))
* Mosaic: Always reset discrete_data ([#3445](../../pull/3445))
* Save Data: Reset writer upon changing the extension ([#3437](../../pull/3437))
* Scatter Plot: Replot when input Features ([#3439](../../pull/3439))
* build-conda-installer: Update the included Miniconda installer ([#3429](../../pull/3429))
* OWDataProjectionWidget: Consider tables with nan-s equal ([#3435](../../pull/3435))
* Projections: Retain embedding if non-relevant variables change ([#3428](../../pull/3428))
* PCA: Rename components to PC1, PC2, PC3, ... ([#3423](../../pull/3423))


[3.18.0] - 2018-11-13
--------------------
##### Enhancements
* tSNE: Move from single-cell to core ([#3379](../../pull/3379))
* Transform: Add new widget ([#3346](../../pull/3346))
* Correlations: Move from prototypes to core ([#3362](../../pull/3362))
* Install widget help files ([#3345](../../pull/3345))
* Feature Statistics: Move from prototypes to core ([#3303](../../pull/3303))
* Replace scikit-learn tSNE with faster implementation ([#3192](../../pull/3192))
* Select Columns: Enable filtering of used features ([#3363](../../pull/3363))

##### Bugfixes
* setup.py: Remove trailing slash from directory names in data_files ([#3394](../../pull/3394))
* condainstall.bat: Add conda.bat and activate.bat scripts again ([#3389](../../pull/3389))
* LearnerScorer: fix for preprocessed data ([#3381](../../pull/3381))
* Feature Statistics: Update outputs on new data ([#3382](../../pull/3382))
* python-framework.sh: Fix 'Current' symlink creation ([#3373](../../pull/3373))
* Fix wrong indices in tooltips in projection when some data was invalid ([#3357](../../pull/3357))
* Scatterplot's VizRank no longer crashes in presence of nonprimitive metas ([#3347](../../pull/3347))
* Predictions: Fix failure after failed predictor ([#3337](../../pull/3337))
* Louvain Clustering: Do not invalidate output on PCA slider change with apply disabled ([#3339](../../pull/3339))
* Use minimal keyring implementation for tests ([#3359](../../pull/3359))
* OWFreeViz: Fix optimization for data with missing values ([#3358](../../pull/3358))


[3.17.0] - 2018-10-26
--------------------
##### Enhancements
* OWSelectAttributes: Use input features ([#3299](../../pull/3299))

##### Bugfixes
* OWDataTable: reset selections on domain change ([#3327](../../pull/3327))
* owlouvainclustering: Fix race conditions ([#3322](../../pull/3322))
* Save data widget crash on no data ([#3311](../../pull/3311))
* OWWidget: Preserve widget geometry between hide/show events ([#3304](../../pull/3304))
* Fix OWWidget destruction ([#3296](../../pull/3296))
* OWWidget: Fix size hint propagation ([#3253](../../pull/3253))


[3.16.0] - 2018-09-14
--------------------
##### Enhancements
* ROC analysis: show thresholds ([#3172](../../pull/3172))
* Edit Domain: Record transformations ([#3231](../../pull/3231))
* Data Table: Enable deselection ([#3189](../../pull/3189))
* Empty helper pane message ([#3210](../../pull/3210))
* Matplotlib output for Scatter plot ([#3134](../../pull/3134))
* Scatterplot: indicate overlap of points. ([#3177](../../pull/3177))
* Selection of format and compression in save data widget ([#3147](../../pull/3147))
* OWBoxPlot: add option to sort discrete distributions by size ([#3156](../../pull/3156))
* Table: speed-up computation of basic stats of given columns. ([#3166](../../pull/3166))
* Canvas: 'Window Groups' continued ([#3085](../../pull/3085))
* Combo Box Search Filter ([#3014](../../pull/3014))
* Widget Insertion ([#3179](../../pull/3179))

##### Bugfixes
* Documentation fetching with redirects ([#3248](../../pull/3248))
* DiscreteVariable reconstruction ([#3242](../../pull/3242))
* io: Handle mismatched number of header/data values ([#3237](../../pull/3237))
* OWNeuralNetwork model pickling ([#3230](../../pull/3230))
* Variable: Prevent hashing of Values of DiscreteVariable. ([#3217](../../pull/3217))


[3.15.0] - 2018-08-06
--------------------
##### Enhancements
* Silhouette Plot: Add cosine distance ([#3176](../../pull/3176))
* Add pandas_compat.table_to_frame(tab) ([#3180](../../pull/3180))
* OWSelectByDataIndex: New widget (move from Prototypes) ([#3181](../../pull/3181))
* Make filters available in Orange.data namespace. ([#3170](../../pull/3170))
* Move Louvain clustering from prototypes to core ([#3111](../../pull/3111))
* OWWidget: Collapse/expand the widget on control area toggle ([#3146](../../pull/3146))
* Rank: SklScorer should use the faster SklImputer. ([#3164](../../pull/3164))
* RecentFiles: Check for missing file in workflow dir ([#3064](../../pull/3064))
* Smart widget suggestions ([#3112](../../pull/3112))
* Match Keywords in Widget Search ([#3117](../../pull/3117))
* io: Speedup write_data ([#3115](../../pull/3115))
* OWEditDomain: Enable reordering of discrete variables ([#3119](../../pull/3119))

##### Bugfixes
* oweditdomain: Fix an IndexError when all rows are deselected ([#3183](../../pull/3183))
* OWFreeViz: fix class density size ([#3158](../../pull/3158))
* OWBoxPlot: Fix empty continuous contingency check ([#3165](../../pull/3165))
* OWSql: enforce data download for non PostgreSQL databases ([#3178](../../pull/3178))
* owlouvainclustering: Make the task completion handler a single slot ([#3182](../../pull/3182))
* OWReport: disable save and print on empty report ([#3175](../../pull/3175))
* RemoveConstant: remove constant NaN features. ([#3163](../../pull/3163))
* utils/concurrent: Switch default thread pool ([#3138](../../pull/3138))
* OWBoxPlot: Fix quartiles computation ([#3159](../../pull/3159))
* OWBoxPlot: Plot axis labels and quartiles label layout ([#3162](../../pull/3162))
* [RFC] OWHeatMap: remove empty clusters from visualization ([#3155](../../pull/3155))
* report: Fix the number of hidden rows. ([#3150](../../pull/3150))
* [RFC] KMeans upgrade sparse support ([#3140](../../pull/3140))
* WebView fixes ([#3148](../../pull/3148))
* ci/appveyor: Update test dependencies ([#3139](../../pull/3139))
* Replace use of obsolete QStyle.standardPixmap ([#3127](../../pull/3127))
* BoxPlot: Hide groups with no instances ([#3122](../../pull/3122))


[3.14.0] - 2018-07-04
--------------------
##### Enhancements
* MergeData: add autocommit button ([#3091](../../pull/3091))
* Canvas: Window Groups ([#3066](../../pull/3066))
* Save data with compression ([#3047](../../pull/3047))
* Neural network widget that works in a separate thread ([#2958](../../pull/2958))
* Display Widgets on Top option ([#3038](../../pull/3038))
* Implement multi window editing ([#2820](../../pull/2820))
* Data Info widget displays data attributes ([#3022](../../pull/3022))
* Icon redesign: k-means, clustering, distances ([#3018](../../pull/3018))

##### Bugfixes
* postgres: Fix wrong discrete values ([#3109](../../pull/3109))
* OWRank: Select Attributes fixes and improvements ([#3084](../../pull/3084))
* EditDomain: Editing TimeVariable broke formatting ([#3101](../../pull/3101))
* OWMosaic: Don't offer String meta attributes ([#3072](../../pull/3072))
* owkmeans: fix initialization choice ([#3090](../../pull/3090))
* Workaround for segfaults with Nvidia on Linux ([#3100](../../pull/3100))
* Canvas: Fix 'Widgest on top' ([#3068](../../pull/3068))
* Re-cythonize with Cython 0.28 for Python 3.7 compatibility ([#3067](../../pull/3067))
* BSD compatibility patch ([#3061](../../pull/3061))
* OWScatterOWScatterPlotGraph: Match group colors with marker colors ([#3053](../../pull/3053))
* listfilter: Fix filter line edit completion ([#2896](../../pull/2896))
* VizRank: Fix race condition in `toggle` ([#3042](../../pull/3042))
* Heat Map: Allow labeling by TimeVariable ([#3026](../../pull/3026))
* Select Columns: Drag/drop ([#3032](../../pull/3032))
* gui: Suppress mouse button release on the combo box popup ([#3025](../../pull/3025))
* tests: Fix time tracking in process_events ([#3041](../../pull/3041))
* test_owmosaic: Cleanup/join threads on test tear down ([#3040](../../pull/3040))
* owselectcolumns: Fix performance on filtering with selection ([#3030](../../pull/3030))
* test: Fix tests for 'Datasets' widget ([#3033](../../pull/3033))
* Sorting add-ons in the alphabetical order ([#3013](../../pull/3013))
* owscatterplot: Use correct data for regression line ([#3024](../../pull/3024))
* Add-ons dialog: Restore state ([#3017](../../pull/3017))
* Feature constructor does not restore features when loading from saved workflow ([#2996](../../pull/2996))
* boxplot labels overlap ([#3011](../../pull/3011))
* owdiscretize: Fix quadratic complexitiy in the n variables ([#3016](../../pull/3016))


[3.13.0] - 2018-04-17
--------------------
##### Enhancements
* canvas/add-ons: Add extra packages via a name input dialog ([#3006](../../pull/3006))
* Variable lists (with QListView) optimizations ([#2994](../../pull/2994))

##### Bugfixes
* Add-ons working again (PyPI JSON interface + local official index) ([#3005](../../pull/3005))
* Fix variable type guessing ([#2998](../../pull/2998))
* Addon dialog crashes when site-packages directory does not exist ([#2988](../../pull/2988))
* Fix reading double quoted text fields ([#2989](../../pull/2989))

[3.12.0] - 2018-04-06
--------------------
##### Enhancements
* owselectrows: Change defaults for 'Purging' ([#2969](../../pull/2969))
* statistics: Speed up countnans for sparse matrices ([#2965](../../pull/2965))

##### Bugfixes
* Sieve Diagram: Fix spacing of axis labels ([#2971](../../pull/2971))
* Fix data reading speed ([#2923](../../pull/2923))
* KMeans clear results on k change, do not recluster ([#2915](../../pull/2915))
* gui.ControlledList: Take a weakref to listBox ([#2962](../../pull/2962))
* WidgetManager: Schedule delayed deletion for managed OWWidgets ([#2963](../../pull/2963))
* domaineditor: Give the VarTableModel a parent ([#2961](../../pull/2961))
* scatterplot: limit number of points displayed in tooltip ([#2980](../../pull/2980))
* Speed-up prediction by circumventing Instance-specific prediction. ([#2959](../../pull/2959))
* Vizrank: Properly shutdown/wait when parent deleted ([#2960](../../pull/2960))
* Test & Score: Make the scores view table non-editable ([#2947](../../pull/2947))


[3.11.0] - 2018-03-07
--------------------
##### Enhancements
* Distances: Optimize PearsonR/SpearmanR ([#2852](../../pull/2852))
* Data Table: Optimize performance ([#2905](../../pull/2905))

##### Bugfixes
* Save Image to SVG fixed on Qt5 ([#2930](../../pull/2930))
* Test & Score: Resort scores when data changes ([#2931](../../pull/2931))
* distribution.py: Fix computation on multiclass data ([#2903](../../pull/2903))
* contingency.pyx: Fix out of bound write ([#2924](../../pull/2924))
* Test and Score: Fix averaging over classes for binary scores ([#2918](../../pull/2918))
* sgd: Change deprecated n_iter to max_iter ([#2920](../../pull/2920))
* heatmap: Do not crash on all zero column ([#2916](../../pull/2916))


[3.10.0] - 2018-02-19
--------------------
##### Enhancements
* Select Rows: Add annotated data output ([#2908](../../pull/2908))
* canvas: Open dropped ows files ([#2885](../../pull/2885))
* Settings for HTTP and HTTPS proxies in the canvas ([#2906](../../pull/2906))
* Add-ons: Option to list only trusted add-ons ([#2899](../../pull/2899))

##### Bugfixes
* SPG Legend: Fix vertical spacing ([#2914](../../pull/2914))


[3.9.1] - 2018-02-02
--------------------
##### Enhancements
* Add parameters and similarity measures to tSNE ([#2510](../../pull/2510))
* Canvas: Add zoom ([#2841](../../pull/2841))

##### Bugfixes
* OWWidget: Store quicktip displayed state in non versioned settings dir ([#2875](../../pull/2875))
* Impute: Fix state serialization/restore ([#2830](../../pull/2830))
* Feature Constructor: Make FeatureFunc picklable ([#2873](../../pull/2873))
* Projection widgets: transform data properly ([#2871](../../pull/2871))


[3.9.0] - 2018-01-19
--------------------
##### Enhancements
* Linear Discriminant Analysis: scripting part ([#2823](../../pull/2823))
* owdistances: Add 'Normalize' check box ([#2851](../../pull/2851))
* Variable: Simplify the is_{discrete,continuous,...} implementation ([#2874](../../pull/2874))
* manifold: Use arpack for decomposition in `torgerson` ([#2825](../../pull/2825))
* Radviz: new widget ([#2480](../../pull/2480))
* Canvas: Improve preview rendering ([#2784](../../pull/2784))
* Linear Projection (LDA, PCA) ([#2445](../../pull/2445))
* Scatter Plot Graph: max discrete values colors and shape ([#2804](../../pull/2804))
* Scatter Plot Graph: legend opacity ([#2819](../../pull/2819))

##### Bugfixes
* Add labels to combo when data comes from distance matrix ([#2866](../../pull/2866))
* utils/concurrent: Handle an empty futures list ([#2834](../../pull/2834))
* OWWidget: Move 'splitter' to private members ([#2847](../../pull/2847))
* Radviz: enable VizRank numeric color (class) ([#2853](../../pull/2853))
* Bincount: Fix crash on array with all nans ([#2831](../../pull/2831))
* mds: Fix incorrect assert ([#2844](../../pull/2844))
* VizRank (Linear Projection, Radviz): spin disabled/enabled ([#2846](../../pull/2846))
* Windows installers: Python lookup ([#2827](../../pull/2827))
* Canvas: Palette propagation ([#2760](../../pull/2760))
* mssql: Catch errors due to incorrect connection params ([#2838](../../pull/2838))
* canvas/addons: Fix progress dialog showing up when not necessary ([#2833](../../pull/2833))
* ContextHandler: Merge local context into globals before serialization ([#2837](../../pull/2837))
* Hierarchical Clustering: Fix size constraints ([#2796](../../pull/2796))
* canvas/annotationitem: Use separate item for shadow base ([#2821](../../pull/2821))
* Scatter Plot Graph: vars instead of indices & remove dead code ([#2815](../../pull/2815))
* Table: classes (Y) must be float ([#2822](../../pull/2822))


[3.8.0] - 2017-12-01
--------------------
##### Enhancements
* New signals: Trees, Forest ([#2801](../../pull/2801))
* Scatter Plot: Improve tooltips ([#2703](../../pull/2703))
* Allow custom (generic) names in Transpose Widget  ([#2737](../../pull/2737))
* Scatter Plot VizRank: some fixes and regard to color ([#2787](../../pull/2787))
* Improved Sparsity Handling ([#2341](../../pull/2341))
* Error Reporting: send report even when recursion error is raised
* test_owdatasets: Test files at different (dir tree) depths
* [FIX] Rank: should not fail on data with no attributes
* Domain: Add copy method ([#2734](../../pull/2734))
* Domain Model: order without separators ([#2697](../../pull/2697))

##### Bugfixes
* Test & Learn: do not crash on a data with class only nans ([#2751](../../pull/2751))
* FreeViz: 2 issues when no data ([#2780](../../pull/2780))
* [ENH] Scatter Plot VizRank: some fixes and regard to color ([#2787](../../pull/2787))
* Scatter Plot Graph: crash on metas column with all 0 values ([#2775](../../pull/2775))
* Scatter Plot: subset data ([#2773](../../pull/2773))
* Scatter Plot: VizRank disabled when no class vars ([#2757](../../pull/2757))
* Error Reporting: send report even when recursion error is raised
* Select Rows: None on output when no data ([#2726](../../pull/2726))
* test_owdatasets: Test files at different (dir tree) depths
* [FIX] Rank: should not fail on data with no attributes
* Predictions: space added before bracket in meta name. ([#2742](../../pull/2742))
* Fix AbstractSortTableModel.mapFromSourceRows for empty list or array ([#2730](../../pull/2730))
* Correspondence Analysis: do not crash when no categorical ([#2723](../../pull/2723))
* ScatterPlotGraph: fix zoom CTRL + and CTRL - ([#2716](../../pull/2716))
* errorreporting: Remove use of pip internal api ([#2724](../../pull/2724))


[3.7.1] - 2017-11-17
--------------------
##### Enhancements
* MDS: Support showing textual features as labels ([#2699](../../pull/2699))

##### Bugfixes
* canvas/canvasmain: Fix 'Examples' action handling in the welcome dialog ([#2779](../../pull/2779))
* Nomogram on PyQt4 ([#2763](../../pull/2763))
* Broken installation of Installation of wheels ([#2765](../../pull/2765))
* Add-on installation crashes (when conda not in PATH) ([#2725](../../pull/2725))


[3.7.0] - 2017-10-27
--------------------
##### Enhancements
* Data Sets: Add filter ([#2695](../../pull/2695))
* Add-on installation with Conda ([#2561](../../pull/2561))
* Add Groups column to Selected Data in Scatter plot output ([#2678](../../pull/2678))
* DomainModel: Don't Show Hidden Variables by Default ([#2690](../../pull/2690))
* FreeViz: new widget ([#2512](../../pull/2512))
* FreeViz script ([#2563](../../pull/2563))
* Boxplot: Allow hiding labels ([#2654](../../pull/2654))
* owmds: Support selection/output of multiple groups ([#2666](../../pull/2666))
* Widget status bar buttons ([#2514](../../pull/2514))
* owfile: allow multiple readers with same extension ([#2644](../../pull/2644))

##### Bugfixes
* Tree Viewer: reset view to top left ([#2705](../../pull/2705))
* ScatterPlot Crashes on Data With Infinity Values ([#2709](../../pull/2709))
* Scatter Plot: regression line: show r instead of k ([#2701](../../pull/2701))
* settings: Do not clear schema_only settings on close_context ([#2691](../../pull/2691))
* Statistics.unique: Fix Sparse Return Order For Negative Numbers ([#2572](../../pull/2572))
* Statistics.countnans/bincount: Fix NaN Counting, Consider Implicit Zeros ([#2698](../../pull/2698))
* MDS: No optimization when subset data ([#2675](../../pull/2675))
* Outliers widget no longer checks classes and doesn't crash on singular covariances matrices ([#2677](../../pull/2677))
* OWRank: Fix autocommit ([#2685](../../pull/2685))
* OWScatterPlot: Change output Feature to AttributeList ([#2689](../../pull/2689))
* OWSql does not save selected table/query ([#2659](../../pull/2659))
* Scatter Plot: Scatter Plot automatically sends selection ([#2649](../../pull/2649))
* Silhouette plot rendering ([#2656](../../pull/2656))
* Variable.make returns proxies ([#2667](../../pull/2667))
* owhierarchicalclustering: Fix performance on deselection ([#2670](../../pull/2670))
* Report Table: Make Table Headers Bold ([#2668](../../pull/2668))
* MDS: primitive metas, init_attr_values ([#2661](../../pull/2661))
* MDS: Primitive metas ([#2648](../../pull/2648))
* MDS: similar pairs and combos not cleared ([#2643](../../pull/2643))
* Scatter Plot: remove dead and commented code, tests ([#2627](../../pull/2627))


[3.6.0] - 2017-09-29
--------------------
##### Enhancements
* PythonScript: Multiple inputs ([#2506](../../pull/2506))
* DomainEditor: Add horizontal header ([#2579](../../pull/2579))
* Feature Constructor: Support additional functions () ([#2611](../../pull/2611))
* Miniconda installer: Install conda executable ([#2616](../../pull/2616))
* Datasets: New widget ([#2557](../../pull/2557))
* Neural Network widget ([#2553](../../pull/2553))

##### Bugfixes
* settings: Store settings version in the serialized defaults ([#2631](../../pull/2631))
* canvas/stackedwidget: Check if the new geometry is the same as the old ([#2636](../../pull/2636))
* OWRank: sort NaNs last; fix sort indicator ([#2618](../../pull/2618))
* Schema-only settings in components ([#2613](../../pull/2613))
* OWBaseLearner: Save learner name in workflow ([#2608](../../pull/2608))
* Saving of multiple selections in ScatterPlot ([#2598](../../pull/2598))
* OWBoxPlot: Faster selection ([#2595](../../pull/2595))
* preprocess.randomization: Do not use the same seed for X, Y, and meta ([#2603](../../pull/2603))
* Slow Rank ([#2494](../../pull/2494))
* setup: Increase required setuptools version ([#2602](../../pull/2602))
* Disable pyqtgraph's exit cleanup handler ([#2597](../../pull/2597))
* ScatterPlotGraph: fix labelling when there are missing data ([#2590](../../pull/2590))
* canvas: Fix link runtime state modeling ([#2591](../../pull/2591))
* Tree: Reintroduce preprocessors. ([#2566](../../pull/2566))
* canvas/preview: Fix workflow preview rendering ([#2586](../../pull/2586))
* Fix saving reports on Python 3.6 ([#2584](../../pull/2584))
* Fix failing report tests ([#2574](../../pull/2574))
* widgets/tests: Compatibility with Python 3.5.{0,1} ([#2575](../../pull/2575))


[3.5.0] - 2017-09-04
--------------------
##### Enhancements
* Proper calculation of distances ([#2454](../../pull/2454))
* OWFeatureConstructor: Add new functions from numpy ([#2410](../../pull/2410))
* Widget status bar ([#2464](../../pull/2464))
* Impute widget: Parallel execution in the background ([#2395](../../pull/2395))

##### Bugfixes
* Mosaic Display: subset data ([#2528](../../pull/2528))
* MDS: Optimize similar pairs graphics construction ([#2536](../../pull/2536))
* Error Reporting: read attached schema file as utf8 ([#2416](../../pull/2416))
* Another color palette when too many colors needed ([#2522](../../pull/2522))
* Widget: splitter sizes ([#2524](../../pull/2524))
* Silhouette Plot: another memory error ([#2521](../../pull/2521))
* Fix asynchronous widget tests ([#2520](../../pull/2520))
* Mosaic Vizrank: compute_attr_order is called every step ([#2484](../../pull/2484))
* widgets/model: Restore 'explicit' hint flag for 'Coefficients' output ([#2509](../../pull/2509))


[3.4.5] - 2017-07-27
--------------------
##### Enhancements
* OWMDS, OWLinearProjection: Save selection in workflow ([#2301](../../pull/2301))
* SQL: Save user name and password via credentials manager ([#2403](../../pull/2403))
* Canvas Annotations: Text markup editing ([#2422](../../pull/2422))
* New windows installer scripts ([#2338](../../pull/2338))

##### Bugfixes
* Tree: Fix min_samples_leaf check ([#2507](../../pull/2507))
* Tree: Support classification on sparse data ([#2430](../../pull/2430))
* Trees: Support regression on sparse data ([#2497](../../pull/2497))
* Trees: Fix predictions on sparse data ([#2496](../../pull/2496))
* Change Variable Icons: Discrete -> Categorical, Continuous -> Numeric ([#2477](../../pull/2477))
* Distributions: Show probabilities upon selection ([#2428](../../pull/2428))
* Manifold Learning: Handling out of memory error ([#2441](../../pull/2441))
* CN2 Rule Induction: Handling out of memory error ([#2397](../../pull/2397))
* Hierarchical Clustering: Explicit geometry transform ([#2465](../../pull/2465))
* Scatter Plot Graph: Legend symbols color ([#2487](../../pull/2487))
* Table: Fix printing data with sparse Y ([#2457](../../pull/2457))
* Ensure visible nodes after opening a workflow. ([#2490](../../pull/2490))
* Select Rows: Removing Unused Values for Discrete Variables in Sparse Data ([#2452](../../pull/2452))
* simple_tree.c: Fix mingw compiler compatibility ([#2479](../../pull/2479))
* Add-ons: Fix Installation of Official Add-ons Through Drag & Drop ([#2481](../../pull/2481))
* Mosaic: Clear when data is disconnected ([#2462](../../pull/2462))
* Create Class: Class name cannot be empty ([#2440](../../pull/2440))
* WidgetSignalsMixin: Fix input/output ordering for 'newstyle' signals ([#2469](../../pull/2469))
* Table: Update `ids` in `Table.__del__` ([#2470](../../pull/2470))
* Preprocess: Fix RemoveNaNClasses / Use existing HasClass ([#2450](../../pull/2450))
* SQL: Fixes for Custom SQL option ([#2456](../../pull/2456))
* OWColor: Fix propagating changes to the output ([#2379](../../pull/2379))
* Distances: Prevent inf numbers ([#2380](../../pull/2380))
* Test and Score: Show default columns ([#2437](../../pull/2437))
* Silhouette Plot: Now setting axis range properly ([#2377](../../pull/2377))
* Logistic Regression: Impute ([#2392](../../pull/2392))
* schemeedit: Clear edit focus before removing items ([#2427](../../pull/2427))
* Disable menu and mouse zoom in all evaluate's plotting widgets. ([#2429](../../pull/2429))
* canvas: Fix proposed connection scoring for dynamic signals ([#2431](../../pull/2431))
* ROC Analysis: Color support for more than 9 evaluation learners ([#2394](../../pull/2394))
* Scatter Plot: Two minor errors ([#2381](../../pull/2381))
* Feature Constructor: No fail when no values ([#2417](../../pull/2417))


[3.4.4] - 2017-06-16
--------------------
##### Enhancements
* SimpleTreeLearner: Release GIL & thread safety ([#2398](../../pull/2398))
* Improve support for HiDPI displays ([#2325](../../pull/2325))
* Add a tutorial section on responsive GUI ([#2318](../../pull/2318))
* Check if updates are available upon startup ([#2273](../../pull/2273))

##### Bugfixes
* Vizrank: interact with gui from main thread only ([#2389](../../pull/2389))
* Some preprocessors couldn not be pickled ([#2409](../../pull/2409))
* MDS: Support distances without domain information ([#2335](../../pull/2335))
* Paint Data: Fix crash on empty data ([#2399](../../pull/2399))
* Distributions: do not crash on empty data ([#2383](../../pull/2383))
* Update checker: LooseVersion does not handle str parts ([#2401](../../pull/2401))
* owpreproces: Stable order of continuizers ([#2400](../../pull/2400))
* owmanifoldlearning: Remove `n_jobs=-2` parameter ([#2371](../../pull/2371))
* Scatter Plot: features and no data ([#2384](../../pull/2384))
* tests: Fix test errors when running with numpy 1.13.0 ([#2396](../../pull/2396))
* OWColor: Use DiscreteVariable values for matching contexts ([#2376](../../pull/2376))
* Outliers: handling memory error ([#2374](../../pull/2374))
* score.FCBF: do not segfault on continuous variables w/ <0 values ([#2355](../../pull/2355))
* Rank widget supports Scorer inputs ([#2350](../../pull/2350))
* Silhouette Plot: handling memory error ([#2336](../../pull/2336))
* Distances: handling errors due to too large arrays ([#2315](../../pull/2315))
* Confusion Matrix: do not append extra column if empty ([#2386](../../pull/2386))


[3.4.3] - 2017-06-03
--------------------
##### Enhancements
* Venn diagram: Support sparse data ([#2334](../../pull/2334))
* PCA: Support sparse data ([#2313](../../pull/2313))
* Impute: Support sparse data ([#2357](../../pull/2357))
* Merge: Support sparse data ([#2305](../../pull/2305))
* Scatter Plot: Support sparse data ([#2152](../../pull/2152))
* Manifold: Support t-SNE on sparse data ([#2281](../../pull/2281))
* Mosaic: Selectable color variable ([#2133](../../pull/2133))
* Test & Score: Allow choosing columns ([#2257](../../pull/2257))
* Preprocess: Add all available methods to feature selection ([#2205](../../pull/2205))
* Scatter Plot: Support string metas labels ([#2360](../../pull/2360))

##### Bugfixes
* Fix and improve Precision, Recall, F1 ([#2369](../../pull/2369))
* Paint Data: Stores data in list and not np.array ([#2314](../../pull/2314))
* Paint Data: Save and load labels ([#2259](../../pull/2259))
* File: No domain or empty domain -> no data ([#2337](../../pull/2337))
* File: Support sparse data in Domain Editor ([#2245](../../pull/2245))
* File: Raise and handle Exc. when file bad pickle ([#2232](../../pull/2232))
* Test & Score: Fix migration of old settings ([#2254](../../pull/2254))
* Test & Score: Show correct error ([#2263](../../pull/2263))
* Test & Score: Instantly recognize new input ([#2247](../../pull/2247))
* Test & Score: Handling memory errors ([#2316](../../pull/2316))
* Tree Viewer: Check if there is selected class value ([#2224](../../pull/2224))
* CredentialManager: Handling password credentials error ([#2354](../../pull/2354))
* RowInstance: Fix sparse check ([#2362](../../pull/2362))
* Cross Validation: Cast fold number to string ([#2348](../../pull/2348))
* Silhouette Plot: Elide hover labels if labels are long ([#2278](../../pull/2278))
* Select Rows, Table: Filtering string values ([#2176](../../pull/2176))
* Report: Handle PermissionError when trying to save ([#2225](../../pull/2225))
* Continuize: Prevent crashing - column with equal and NaN values ([#2144](../../pull/2144))
* Add-ons: Handling ValueError due to connection problems ([#2239](../../pull/2239))
* Correspondence: Prevent crashing when cont attr has one value ([#2149](../../pull/2149))
* WebEngineView: Insert JS if loading already started ([#2230](../../pull/2230))
* Manifold Learning: handling numpy LinAlgError ([#2228](../../pull/2228))
* MDS: Fix widget update scheduling ([#2211](../../pull/2211))
* Settings: Handle permission errors when saving settings ([#2147](../../pull/2147))
* Map: Minor fixes ([#2356](../../pull/2356))


[3.4.2] - 2017-04-19
--------------------
##### Enhancements
* Nomogram: Support for sparse data ([#2197](../../pull/2197))
* Add PDF format to image exporters ([#2210](../../pull/2210))
* Reimplement macOS application (.app) build scripts ([#2217](../../pull/2217))
* Canvas: Use 'windowFilePath' to display display current filename instead of the title ([#2206](../../pull/2206))
* OWTestLearners: Cross validation by feature ([#2145](../../pull/2145))
* Pythagorean tree: Make border scale invariant ([#2141](../../pull/2141))

##### Bugfixes
* Scatterplot crashes when loaded from old workflow ([#2241](../../pull/2241))
* Error Report: URL changed ([#2220](../../pull/2220))
* Scatter Plot: update class density ([#2238](../../pull/2238))
* KMeans: should not crash when there is less data rows than k ([#2172](../../pull/2172))
* Edit Domain: Prevent duplicate variable names ([#2146](../../pull/2146))
* Scatter Plot: left margin (axis y) is now adapting ([#2200](../../pull/2200))
* Predictions widget: handle similar but different domains ([#2129](../../pull/2129))
* OWNomogram: Do not paint scene until the widget is not open ([#2202](../../pull/2202))
* Test & Score: crashing prevented when learner disconnects ([#2194](../../pull/2194))
* Widget Logistic Regression: can handle unused values ([#2116](../../pull/2116))
* stats: Open Corpus in OWDataTable after transposing it ([#2177](../../pull/2177))
* Rank: fixes creating Table with infinite numbers ([#2168](../../pull/2168))
* Add-ons: Problems with datetime parsing ([#2196](../../pull/2196))
* OWPredictions: Allow classification when data has no target column ([#2183](../../pull/2183))
* OWDataSampler: Fix typo boostrap to bootstrap ([#2195](../../pull/2195))
* All widgets are set to auto* when they are used for the first time ([#2136](../../pull/2136))
* Preprocess: enums containing function names changed ([#2151](../../pull/2151))
* Fitter: Fix used_vals and params not being set ([#2138](../../pull/2138))
* VizRankDialog: Stop computation when parent widget is deleted ([#2118](../../pull/2118))
* Distributions Report: Visualizations are now fitted ([#2130](../../pull/2130))
* Fitter: Change params uses default if None ([#2127](../../pull/2127))
* Fix invalid settings reuse in File widget  ([#2137](../../pull/2137))
* Scatter Plot - Prevent crash due to missing data ([#2122](../../pull/2122))
* Sieve Diagram: Using datasets with meta data ([#2098](../../pull/2098))


[3.4.1] - 2017-03-16
--------------------
##### Enhancements
* Scatterplot: Implement grouping of selections ([#2070](../../pull/2070))

##### Bugfixes
* Discover widgets when some dependencies are missing ([#2103](../../pull/2103))
* Select Rows: "is defined" fails ([#2087](../../pull/2087))
* report comments and OWFile reporting filename ([#1956](../../pull/1956))
* owcorrespondence: Handle variables with one value ([#2066](../../pull/2066))
* OWTreeViewer: Fix trees being displayed differently for same tree object ([#2067](../../pull/2067))
* Fitter: Properly delegate preprocessors ([#2093](../../pull/2093))


[3.4.0] - 2017-03-06
--------------------
##### Enhancements
* OWSGD: Output coefficients ([#1981](../../pull/1981))
* OWNomogram: Add a new widget ([#1936](../../pull/1936))
* OWRandomize: Add a new widget ([#1863](../../pull/1863))
* Map widget ([#1735](../../pull/1735))
* Table.transpose: Use heuristic to guess data type of attributes of attributes ([#1844](../../pull/1844))
* Create Class widget ([#1766](../../pull/1766))

##### Bugfixes
* Heatmap: Fix crash on data with empty columns ([#2057](../../pull/2057))
* ScatterPlot: Fix crash when coloring by column of unknowns ([#2061](../../pull/2061))
* owpreprocess: Handle columns with only NaN values ([#2064](../../pull/2064))
* File: Disallow changing string columns to datetime ([#2050](../../pull/2050))
* OWKMeans: Auto-commit fix and silhuette optimization ([#2073](../../pull/2073))
* OWDistributions: Fix binning of meta attributes ([#2068](../../pull/2068))
* SelectRows: Fix loading of conditions ([#2065](../../pull/2065))
* OWRandomize: New icon ([#2069](../../pull/2069))
* ZeroDivisionError owmosaic.py ([#2046](../../pull/2046))
* OWMosaic:  Fix crash for empty column ([#2006](../../pull/2006))
* Fitter: Fix infinite recursion in __getattr__ ([#1977](../../pull/1977))
* OWTreeGraph: Update node text when selecting target class ([#2045](../../pull/2045))
* Prevent PickleError (owfile.py) ([#2039](../../pull/2039))
* Fix Chi2 computation for variables with values with no instances ([#2031](../../pull/2031))
* OWDistanceMatrix: Remove quotes with string labels ([#2034](../../pull/2034))
* owheatmap: Prevent sliders to set Low >= High ([#2025](../../pull/2025))
* WebviewWidget: WebEngine Don't Grab Focus on setHtml ([#1983](../../pull/1983))
* OWFile: Show error msg when file doesn't exists ([#2024](../../pull/2024))
* Preprocess Widget: Continuize type error ([#1978](../../pull/1978))
* data/io.py Metadata file not saved anymore when it is empty ([#2002](../../pull/2002))
* Import from AnyQt instead from PyQt4 ([#2004](../../pull/2004))
* OWNomogram: Adjust scene rect ([#1982](../../pull/1982))
* owconcatenate: Fix domain intersection (remove duplicates) ([#1967](../../pull/1967))
* preprocess: Reset number_of_decimals after scaling ([#1914](../../pull/1914))
* Treeviewer sklearn tree compatibility ([#1870](../../pull/1870))
* OWSVR: Update learner when SVR type changes ([#1878](../../pull/1878))
* Tree widget binarization ([#1837](../../pull/1837))


[3.3.12] - 2017-02-14
--------------------
##### Bugfixes
* Highcharts: Fix freezing on Qt5 ([#2015](../../pull/2015))
* Handle KeyError Mosaic Display (owmosaic.py) ([#2014](../../pull/2014))
* Loading iris on C locale ([#1998](../../pull/1998))
* Handle KeyError Sieve Diagram widget (owsieve) when one row ([#2007](../../pull/2007))
* Test Learners: Fix AUC for selected single target class ([#1996](../../pull/1996))
* OWDataSampler: Fix 'Fixed proportion of data' option ([#1995](../../pull/1995))


[3.3.11] - 2017-02-03
--------------------
##### Enhancements
* Widget testing utilities ([#1939](../../pull/1939))

##### Bugfixes
* KMeans: Fix crashes when underlying algorithm fails ([#1974](../../pull/1974))
* owpaintdata: Adjust color model to input dataset ([#1988](../../pull/1988))
* scatterplot: Fix density image ([#1990](../../pull/1990))
* owpaintdata: Fix an error when the input dataset contains NaN ([#1972](../../pull/1972))
* Table: Ensure correct dtype in `_compute_distributions` ([#1969](../../pull/1969))
* Evaluation Results input validation ([#1954](../../pull/1954))
* owimpute: Fix editing of individual imputers ([#1966](../../pull/1966))
* gui: Trigger callback in SpinBoxWFocusOut only if value changed ([#1979](../../pull/1979))
* Python 3.6 compatibility ([#1963](../../pull/1963))
* File: Fix crash when last_path is None ([#1961](../../pull/1961))
* Paint Data: in-place output modification ([#1959](../../pull/1959))
* DataSampler: Fix crash when stratifying unbalanced datasets ([#1952](../../pull/1952))
* Table.__repr__: Fix for sparse data with < 5 instances ([#1951](../../pull/1951))
* Catch errors during learning in learner widgets ([#1949](../../pull/1949))
* OWMosaic: Fix crash for attribute with no values ([#1941](../../pull/1941))
* Impute: Fix crash when model-based imputation fails ([#1937](../../pull/1937))
* OWSieve: Fix crash for attribute with no values ([#1934](../../pull/1934))
* Tree: Fix crash when two attributes equal number of values ([#1928](../../pull/1928))
* Store changed variables in File widget ([#1805](../../pull/1805))


[3.3.10] - 2017-01-18
--------------------
##### Enhancements
* Input/output signal replacement declarations ([#1810](../../pull/1810))

##### Bugfixes
* MDS Widget: Handle NaN values for plot point styling ([#1931](../../pull/1931))
* OWPCA: Fix crash for dataset with no rows or no attributes ([#1915](../../pull/1915))
* OWMosaic: Discretize metas as well ([#1912](../../pull/1912))
* owfeaturecontructor: Fix an IndexError accessing exception's args ([#1905](../../pull/1905))
* owrank: Remove `super()` call from `migrate_settings` ([#1902](../../pull/1902))
* OWBoxPlot: Fix ordering of boxes ([#1900](../../pull/1900))
* canvas/readwrite: Fix byte literal serialization ([#1898](../../pull/1898))
* owpca: Handle the case of 0 total variance in the PCA solution ([#1897](../../pull/1897))
* Copy data attributes for annotated data set ([#1895](../../pull/1895))
* colorpalette: Fix AttributeError ([#1889](../../pull/1889))
* OWDistributions: Reset combos when data is removed ([#1887](../../pull/1887))
* Concatenate bugfix ([#1886](../../pull/1886))
* OWPredictions: Fix crash when opening report ([#1884](../../pull/1884))
* owsilhouetteplot: Fix TypeError when cluster column is an object array ([#1876](../../pull/1876))
* OWSave: Safer Check if Writer Support Sparse ([#1864](../../pull/1864))
* OWImageViewer: Fix selection with missing values ([#1861](../../pull/1861))
* owselectcolumns: Fix auto commit on any change ([#1859](../../pull/1859))
* Table.transpose: Keep metas array two dimensional when no attributes in domain ([#1855](../../pull/1855))
* Select Rows filter enum ([#1854](../../pull/1854))
* Scatter plot: don't crash on report without data ([#1840](../../pull/1840))
* Crash on ctrl-c/cmd-c in widgets without graphs ([#1827](../../pull/1827))
* Fix crash in listview if labels are changed before calling __setitem__ ([#1825](../../pull/1825))
* Scatterplot: Allow labelling by string attributes ([#1812](../../pull/1812))
* Fix copy to clipboard in "Data Table" widget ([#1808](../../pull/1808))
* TreeGraph: Compatibility with old schemas ([#1804](../../pull/1804))


[3.3.9] - 2016-12-02
--------------------
##### Enhancements
* OWTranspose: Add a new widget ([#1738](../../pull/1738))
* Add appveyor configuration ([#1693](../../pull/1693))
* Vizrank indicators and filters ([#1746](../../pull/1746))
* OWManifoldLearning: MDS - enable PCA initialization ([#1702](../../pull/1702))
* Add VizRank to Mosaic ([#1699](../../pull/1699))
* Setting migration ([#1724](../../pull/1724))
* widget: Allow subclasses to disable the default message bar widget ([#1543](../../pull/1543))
* Manifold Learning ([#1624](../../pull/1624))
* SQL Server support in SQL widget ([#1674](../../pull/1674))
* Visualize widgets: Output Annotated data and Fixups ([#1677](../../pull/1677))
* Add support for PyQt5 ([#1171](../../pull/1171))
* Simple benchmarking suite. ([#1510](../../pull/1510))
* Canvas: Always show the link dialog if the user holds Shift ([#1673](../../pull/1673))
* Scatterplot, HeatMap, TreeGraph, ConfusionMatrix and Unsupervised widgets: Output Flagged Data  ([#1655](../../pull/1655))
* CN2RuleViewer: Output sample of training data in absence of separate data ([#1667](../../pull/1667))
* Metadata for data files ([#1603](../../pull/1603))

##### Bugfixes
* owrank: Add migrate_settings ([#1797](../../pull/1797))
* owconfusionmatix: Add migrate_settings ([#1796](../../pull/1796))
* Improve ada boost widget ([#1787](../../pull/1787))
* OWBoxPlot: Fixups ([#1783](../../pull/1783))
* Filter: Fix FilterContinuous eq operator ([#1784](../../pull/1784))
* OWDistanceMatrix: attribute in context ([#1761](../../pull/1761))
* Hierarchical clustering: Make annotation a context setting ([#1748](../../pull/1748))
* Fix varius deprecation (and other) warnings ([#1774](../../pull/1774))
* Fix transformation for non primitive variables ([#1770](../../pull/1770))
* TimeVariable: don't crash Data Table when reloading and Visualize ... ([#1760](../../pull/1760))
* OWDistances: Mahalanobis wrong dimensions notification. ([#1762](../../pull/1762))
* Switch Sieve to DomainModel, which also fixes VizRank crash on meta attributes ([#1642](../../pull/1642))
* Confusion matrix: Map annotated data through row_indices, add probabi… ([#1720](../../pull/1720))
* OWLoadClassifier: Show message on unpickling error ([#1752](../../pull/1752))
* Silhouette Plot: Fixes ([#1747](../../pull/1747))
* canvas/toolgrid: Remove (unused) mouse press event tracking ([#1740](../../pull/1740))
* Box plot: Handle situation when quantiles can't be computed ([#1742](../../pull/1742))
* owfile: Hide apply button after resetting editor_model ([#1711](../../pull/1711))
* oweditdomain: Initialize `var` attribute ([#1731](../../pull/1731))
* FeatureConstructor: Fix crash when new variable is created without data ([#1733](../../pull/1733))
* owpythonscript: Fix QFileDialog.get{Save,Open}FileName usage ([#1726](../../pull/1726))
* DendrogramWidget: Prevent a zero division error ([#1725](../../pull/1725))
* owfile: Skip add_origin if no filename ([#1717](../../pull/1717))
* OWFile: Do not load large files automatically ([#1703](../../pull/1703))
* Do not show messages when data is removed ([#1706](../../pull/1706))
* Confusion Matrix: Show error on regression results ([#1709](../../pull/1709))
* Fix tests ([#1698](../../pull/1698))
* Scatter Plot: Fix a error when restoring from pre DomainModel workflows ([#1672](../../pull/1672))
* Tree Scorers: Change 'int64_t' to 'intp_t' for platform independence ([#1687](../../pull/1687))
* OWTable: Sort Continuous metas as floats; not strings ([#1678](../../pull/1678))
* Error Reporting: Temporary last open/save directory ([#1676](../../pull/1676))
* TableModel: Don't crash on empty sparse data ([#1675](../../pull/1675))
* Statistics.util.stats: Fix negative #nans for sparse ([#1659](../../pull/1659))
* MDS Widget: Fix zero length line, gray square bug ([#1670](../../pull/1670))
* Fix an error when using latest pyqtgraph develop snapshot ([#1662](../../pull/1662))
* OWHeatMap: Resend 'Selected Data' when settings change ([#1664](../../pull/1664))
* Fix pythagoras tree tooltip for regression trees ([#1660](../../pull/1660))
* OWConfusionMatrix: Output None when no data is selected ([#1653](../../pull/1653))
* OWBoxPlot: Reset widget's appearance when data is removed ([#1654](../../pull/1654))


[3.3.8] - 2016-10-11
--------------------
##### Enhancements
* CredentialManager: Store passwords in System Keyring Services ([#1641](../../pull/1641))
* Extend widget creation policy ([#1611](../../pull/1611))
* File widget improvements ([#1607](../../pull/1607))
* Remote reporting of unexpected errors ([#1558](../../pull/1558))
* OWRank: Widget improvements ([#1560](../../pull/1560))
* canvas: Indicate runtime state on links ([#1554](../../pull/1554))
* Rule induction (CN2) ([#1397](../../pull/1397))
* Upgrade OWSvm unittests ([#1499](../../pull/1499))
* Enable Ward clustering in Hierarchical clustering widget  ([#1515](../../pull/1515))
* PCA transformation speedup ([#1539](../../pull/1539))

##### Bugfixes
* owsql: Fix bug when using connection before established ([#1638](../../pull/1638))
* Scatterplot: Reintroduce sliders for size and opacity ([#1622](../../pull/1622))
* Reporting tabular in Data Table and Rank widgets ([#1573](../../pull/1573))
* BoxPlot crashes on variables with no known values (Fixes #1568) ([#1647](../../pull/1647))
* Canvas: Replace illegal file-name characters with _ when saving workf… ([#1644](../../pull/1644))
* OWScatterPlot: Fix progress bar percentages running over 100% ([#1645](../../pull/1645))
* OWFile: Report errors for incorrect file formats instead of crashing ([#1635](../../pull/1635))
* OWFeatureConstructor: Fix domain check for only meta data sets ([#1632](../../pull/1632))
* gui.lineEdit: Restore changed state tracking ([#1630](../../pull/1630))
* errorreporting: Fix an KeyError for a missing 'Widget Module' entry ([#1625](../../pull/1625))
* win-installer: Build scikit-learn in windows installer build script ([#1623](../../pull/1623))
* owlearnerwidget: Fix output initialization ([#1562](../../pull/1562))
* gui: Add a push button class adapted for variable width text ([#1614](../../pull/1614))
* :  Add `Explicit` flag to supplementary Table learner output ([#1617](../../pull/1617))
* Silhouette: Auto-commit on changing checkbox state ([#1606](../../pull/1606))
* Linear regression: Fix Elastic net; Fix Auto-apply buttons ([#1601](../../pull/1601))
* ROC Analysis - Fix roc averaging ([#1595](../../pull/1595))
* OWBaseLearner: Do not re-fit if name has changed ([#1580](../../pull/1580))
* Context attributes with metas in Sieve and Mosaic ([#1545](../../pull/1545))
* Variable: Fix Variable.copy for StringVariable and TimeVariable ([#1589](../../pull/1589))
* Stats: Fix counting of missing values for non-numeric data ([#1585](../../pull/1585))
* Load Classifier widget sends classifier on init ([#1584](../../pull/1584))
* Context settings ([#1577](../../pull/1577))
* Fixed svg function to return svg chart together with container div for highcharts ([#1541](../../pull/1541))
* Fix compatibility with Color widget ([#1552](../../pull/1552))
* Gini impurity: formula and docstring fixed. ([#1495](../../pull/1495))
* owimageviewer: Open local images directly ([#1550](../../pull/1550))
* Fix loading of datasets with paths in variable attributes ([#1549](../../pull/1549))
* Confusion matrix: fix selected_learner setting ([#1523](../../pull/1523))
* canvas/addons: Remove wrong/unnecessary proxy mapping ([#1533](../../pull/1533))
* Scatterplot: Score Plots crash if multiple attributes have the same score ([#1535](../../pull/1535))
* ScatterPlot: Score Plots window title changed to title case ([#1525](../../pull/1525))
* Predictions: column size hint ([#1514](../../pull/1514))


[3.3.7] - 2016-08-05
--------------------
##### Enhancements
* ImageViewer: Add a 'Preview' like window ([#1402](../../pull/1402))
* Pythagorean tree and Pythagorean forest widgets ([#1441](../../pull/1441))
* New workflow examples for the Welcome screen ([#1438](../../pull/1438))
* Test widgets on travis ([#1417](../../pull/1417))
* Save painted data to schema ([#1452](../../pull/1452))
* Welcome screen: New icons for welcome screen ([#1436](../../pull/1436))
* SqlTable: Automatically recognize date/time fields ([#1424](../../pull/1424))
* Proxy support for add-on installation ([#1379](../../pull/1379))
* Automatically create required SQL extensions ([#1395](../../pull/1395))
* Ranking for Sieve, refactoring of Sieve, Mosaic and VizRank ([#1382](../../pull/1382))
* Rank adopted for sparse data ([#1399](../../pull/1399))
* PCA: Add lines and labels showing the explained variance ([#1383](../../pull/1383))
* Implement copying graph to clipboard using Ctrl-C (Cmd-C) ([#1386](../../pull/1386))
* Parallelized cross validation & other evaluation methods ([#1004](../../pull/1004))
* Image viewer thumbnail size ([#1381](../../pull/1381))

##### Bugfixes
* Table names set by readers ([#1481](../../pull/1481))
* OWRandomForest: Fix, refactor and widget tests ([#1477](../../pull/1477))
* KNN: Fix crash when Mahanalobis metric is used ([#1475](../../pull/1475))
* Fix AdaBoost widgets and add some tests ([#1474](../../pull/1474))
* Table: Fix ensure_copy for sparse matrices ([#1456](../../pull/1456))
* statistics.utils: Fix stats for sparse when last column missing ([#1432](../../pull/1432))
* MDS and Distances widges fix ([#1435](../../pull/1435))
* OWBaseLearner: Learner name is changed on output when user changes it and auto apply selected  ([#1453](../../pull/1453))
* Stop advertising support for weights in LogisticRegression. ([#1448](../../pull/1448))
* OWScatterPlot: Fix information message reference ([#1440](../../pull/1440))
* Fix Tree preprocessor order. ([#1447](../../pull/1447))
* SqlTable: Cast to text for unknown and string types ([#1430](../../pull/1430))
* NaiveBayes: Handle degenerate cases ([#1442](../../pull/1442))
* OWBoxPlot: Show corresponding label when ploting discrete variable ([#1400](../../pull/1400))
* Lin and Log Regression: Prevent double commit ([#1401](../../pull/1401))
* KMeans: Silhouette score format precision fixed to integer ([#1434](../../pull/1434))
* Select Rows: skip undefined TimeVariable filters ([#1429](../../pull/1429))
* OWTestLearners: Fix reporting results table ([#1421](../../pull/1421))
* Scatterplot: Sends none in no instance selected ([#1428](../../pull/1428))
* PCA: Fix the variance spin. ([#1396](../../pull/1396))
* overlay: Auto disconnect when the overlay widget is deleted ([#1412](../../pull/1412))
* PaintData: Send None instead of empty table when the plot is empty ([#1425](../../pull/1425))
* Rand Forest Class: Min sample replaces max leaf nodes ([#1403](../../pull/1403))
* SelectRows: Index attrs only by visible in set_data ([#1398](../../pull/1398))
* File: Stores filenames to image attributes ([#1393](../../pull/1393))
* Fix a logging error on windows ([#1390](../../pull/1390))
* OWLearnerWidget: Don't crash when training data contains no features  ([#1389](../../pull/1389))
* TimeVariable: fix repr rounding and repr for nan ([#1387](../../pull/1387))


[3.3.6] - 2016-06-24
--------------------
##### Enhancements
* Automatically discretize continuous variables in Sieve (#1372)
* Nicer reporting of tabular data (e.g. Confusion Matrix) (#1309)
* Match only at the beginning of word in quickMenu search (#1363)
* Univar regression coefficients (#1186)
* Add Mahalanobis distance (#1355)
* Move Data Table higher in the 'suggested widgets' list (#1346)
* Support Distances on Sparse Data (#1345)
* Add auto apply to Test & Score (#1344)
* Support unix timestamps in TimeVariable (#1335)
* Skip Hidden Attributes in SelectRows (#1324)
* Preprocessor widget: add random feature selection (#1112)
* Sort list of available add-ons (#1305)

##### Bugfixes
* Fix and improve handling of input data in PaintData (#1342)
* Mosaic: Output original data, not discretized (#1371)
* Fix report for OWLinearProjection (#1360, #1361)
* Fix auto-apply in SelectColumns (#1353)
* Fix a RuntimeError in ImageViewer when clearing the scene (#1357)
* Use https to query pypi API (#1343)
* OWBaseLearner: Add name attribute to learner (#1340)
* Pass data attributes after preprocess. (#1333)
* Better support for sparse data (#1306)
* gui.auto_commit: Fix crash when caller gives  argument (#1325)
* Fix image viewer runtime error (#1322)
* Better selection indicators in canvas (#1308)
* Open correct help file for add-on widgets (#1307)


[3.3.5] - 2016-06-01
--------------------
##### Bugfixes
* Revert hack that caused missing icons in osx build
* Fix installation in environments without numpy installed (#1291)
* Allow running of library test when PyQt is not available (#1289)


[3.3.4] - 2016-05-27
--------------------
##### Enhancements
* Install add-on by dragging zip/tgz/wheel onto the addons dialog (#1269)
* Added missing reports (#1270, #1271, #1272, #1273, #1279)
* Add auto-apply checkboxes to learner dialogs (#1263)
* Sort numeric values as numbers in Table (#1255)
* Open dragged files in OWFile (#1176)
* Support context cloning in FeatureConstructor and SelectRows (#1196)

##### Bugfixes
* Depend on scikit-learn>=0.17 (#1277)
* Fix installation problem on windows (#1278)
* Fix crash in K-means when silhouette cannot be computed (#1247)
* Fix crash in Distributions on empty data (#1246)
* Reset outputs in MergeData (#1240)
* Compute distances between constructed Instances (#1242)
* Fix links in Changelog (#1244)


[3.3.3] - 2016-05-03
--------------------
##### Bugfixes
* Revert installation of desktop launcher on Linux (#1218)
* Fix a crash when learner is connected to file (#1220)


[3.3.2] - 2016-04-22
--------------------
##### Enhancements
* New preprocessors ReliefF and FCBF (#1133)
* New feature scorers ANOVA, Chi2 and Univariate Linear Regression (#1125)
* Mosaic plot can display numeric attributes (#1165)
* Sheet selection for excel files in File widget (#1164)
* Check code quality with pylint on Travis (#1121)
* Improve PyQt5 forward compatibility (#1029)
* Include default datasets in File widget (#1174)
* Install desktop launcher on linux (#1205)

##### Bugfixes
* Fix a bug in nested transformation of instance (#1192)
* Fix vizrank's crash when pair had no valid data (#1189)
* Save Graph doesn't save axes (#1134)
* Open included tutorials with correct dataset (#1169)
* Disable bsp index on the main canvas scene (#1168, #1173)
* Fix FeatureConstructor crash with Python 3.5 (#1157)
* Include feature names in Rank Widget report (#1022)
* Decrease memory consumption in PCA (#1150)
* Fix dragging of treshold in PCA widget (#1051)
* Save TimeVariable in ISO 8601 format (#1145)
* Allow use of feature metadata in MDS (#1130)
* OWSelectColumns: fix drag and drop for features with Unicode names (#1144)
* Warn when setting values are not present on instance (#1139)
* Fix printing of Table with a TimeVariable (#1141)
* Fix Test & Score widget report (#1138)


[3.3.1] - 2016-03-24
--------------------
##### Enhancements
* Rank widget outputs scores
* SGD Regression widget: Fixed layout and added reporting

##### Bugfixes
* Windows installer: update pip on target system if required

[3.3] - 2016-03-18
------------------
##### Enhancements
*  Changed layout of File widget
*  Distance matrix widget
*  Meta attributes in Sieve and Mosaic
*  New type of variable: Time variable
*  Report for Distance Transformation widget
*  Report for Linear Regression widget
*  Report for Univariate Regression (fixes #1080)
*  score.FCBF: a Fast Correlation-Based Filter for feature selection
*  Sieve enhancements
*  Silhouette Plot widget
*  Venn Diagram: Add option to output unique/all instances.
*  Widgets for saving and loading distances

##### Bugfixes
*  breast-cancer.tab: change type of tumor-size column (fixes #1065)
*  Color: Do not resize columns in continuous table to contents (fixes #1055)
*  Exporting graphs to dot format
*  OWDistributions: Do not remove constant attribute and do not draw if there is no data
*  OWRank: Give name to Scores output table
*  OWRank no longer crashes when additional learners are available
*  ReliefF: Support for missing target values
*  Report: Fix crash on reporting tables
*  RF without pruning by default

##### Documentation
* Update build/install/contributing READMEs
* Update documentation in widget.rst

[3.2] - 2016-02-12
------------------
* Finalized Orange 3.2, switched to stable(r) release cycles

[0.1] - 1996-10-10
----------------
* Initial version based on Python 1.5.2 and Qt 2.3


[next]: https://github.com/biolab/orange3/compare/3.31.1...HEAD
[3.31.1]: https://github.com/biolab/orange3/compare/3.31.0...3.31.1
[3.31.0]: https://github.com/biolab/orange3/compare/3.30.2...3.31.0
[3.30.2]: https://github.com/biolab/orange3/compare/3.30.1...3.30.2
[3.30.1]: https://github.com/biolab/orange3/compare/3.30.0...3.30.1
[3.30.0]: https://github.com/biolab/orange3/compare/3.29.3...3.30.0
[3.29.3]: https://github.com/biolab/orange3/compare/3.29.2...3.29.3
[3.29.2]: https://github.com/biolab/orange3/compare/3.29.1...3.29.2
[3.29.1]: https://github.com/biolab/orange3/compare/3.29.0...3.29.1
[3.29.0]: https://github.com/biolab/orange3/compare/3.28.0...3.29.0
[3.28.0]: https://github.com/biolab/orange3/compare/3.27.1...3.28.0
[3.27.1]: https://github.com/biolab/orange3/compare/3.27.0...3.27.1
[3.27.0]: https://github.com/biolab/orange3/compare/3.26.0...3.27.0
[3.26.0]: https://github.com/biolab/orange3/compare/3.25.1...3.26.0
[3.25.1]: https://github.com/biolab/orange3/compare/3.25.0...3.25.1
[3.25.0]: https://github.com/biolab/orange3/compare/3.24.1...3.25.0
[3.24.1]: https://github.com/biolab/orange3/compare/3.24.0...3.24.1
[3.24.0]: https://github.com/biolab/orange3/compare/3.23.1...3.24.0
[3.23.1]: https://github.com/biolab/orange3/compare/3.23.0...3.23.1
[3.23.0]: https://github.com/biolab/orange3/compare/3.22.0...3.23.0
[3.22.0]: https://github.com/biolab/orange3/compare/3.21.0...3.22.0
[3.21.0]: https://github.com/biolab/orange3/compare/3.20.1...3.21.0
[3.20.1]: https://github.com/biolab/orange3/compare/3.20.0...3.20.1
[3.20.0]: https://github.com/biolab/orange3/compare/3.19.0...3.20.0
[3.19.0]: https://github.com/biolab/orange3/compare/3.18.0...3.19.0
[3.18.0]: https://github.com/biolab/orange3/compare/3.17.0...3.18.0
[3.17.0]: https://github.com/biolab/orange3/compare/3.16.0...3.17.0
[3.16.0]: https://github.com/biolab/orange3/compare/3.15.0...3.16.0
[3.15.0]: https://github.com/biolab/orange3/compare/3.14.0...3.15.0
[3.14.0]: https://github.com/biolab/orange3/compare/3.13.0...3.14.0
[3.13.0]: https://github.com/biolab/orange3/compare/3.12.0...3.13.0
[3.12.0]: https://github.com/biolab/orange3/compare/3.11.0...3.12.0
[3.11.0]: https://github.com/biolab/orange3/compare/3.10.0...3.11.0
[3.10.0]: https://github.com/biolab/orange3/compare/3.9.1...3.10.0
[3.9.1]: https://github.com/biolab/orange3/compare/3.9.0...3.9.1
[3.9.0]: https://github.com/biolab/orange3/compare/3.8.0...3.9.0
[3.8.0]: https://github.com/biolab/orange3/compare/3.7.1...3.8.0
[3.7.1]: https://github.com/biolab/orange3/compare/3.7.0...3.7.1
[3.7.0]: https://github.com/biolab/orange3/compare/3.6.0...3.7.0
[3.6.0]: https://github.com/biolab/orange3/compare/3.5.0...3.6.0
[3.5.0]: https://github.com/biolab/orange3/compare/3.4.5...3.5
[3.4.5]: https://github.com/biolab/orange3/compare/3.4.4...3.4.5
[3.4.4]: https://github.com/biolab/orange3/compare/3.4.3...3.4.4
[3.4.3]: https://github.com/biolab/orange3/compare/3.4.2...3.4.3
[3.4.2]: https://github.com/biolab/orange3/compare/3.4.1...3.4.2
[3.4.1]: https://github.com/biolab/orange3/compare/3.4.0...3.4.1
[3.4.0]: https://github.com/biolab/orange3/compare/3.3.12...3.4.0
[3.3.12]: https://github.com/biolab/orange3/compare/3.3.11...3.3.12
[3.3.11]: https://github.com/biolab/orange3/compare/3.3.10...3.3.11
[3.3.10]: https://github.com/biolab/orange3/compare/3.3.9...3.3.10
[3.3.9]: https://github.com/biolab/orange3/compare/3.3.8...3.3.9
[3.3.8]: https://github.com/biolab/orange3/compare/3.3.7...3.3.8
[3.3.7]: https://github.com/biolab/orange3/compare/3.3.6...3.3.7
[3.3.6]: https://github.com/biolab/orange3/compare/3.3.5...3.3.6
[3.3.5]: https://github.com/biolab/orange3/compare/3.3.4...3.3.5
[3.3.4]: https://github.com/biolab/orange3/compare/3.3.3...3.3.4
[3.3.3]: https://github.com/biolab/orange3/compare/3.3.2...3.3.3
[3.3.2]: https://github.com/biolab/orange3/compare/3.3.1...3.3.2
[3.3.1]: https://github.com/biolab/orange3/compare/3.3...3.3.1
[3.3]: https://github.com/biolab/orange3/compare/3.2...3.3
[3.2]: https://github.com/biolab/orange3/compare/3.1...3.2
[0.1]: https://web.archive.org/web/20040904090723/http://www.ailab.si/orange/
<p align="center">
    <a href="https://orange.biolab.si/download">
    <img src="https://raw.githubusercontent.com/irgolic/orange3/README-shields/distribute/orange-title.png" alt="Orange Data Mining" height="200">
    </a>
</p>
<p align="center">
    <a href="https://orange.biolab.si/download" alt="Latest release">
        <img src="https://img.shields.io/github/v/release/biolab/orange3?label=download" />
    </a>
    <a href="https://orange3.readthedocs.io/en/latest/?badge=latest" alt="Documentation">
        <img src="https://readthedocs.org/projects/orange3/badge/?version=latest">
    </a>
    <a href="https://discord.gg/FWrfeXV" alt="Discord">
        <img src="https://img.shields.io/discord/633376992607076354?logo=discord&color=7389D8&logoColor=white&label=Discord">                                                                                                                                                                                                                                                  </a>
</p>

# Orange Data Mining
[Orange] is a data mining and visualization toolbox for novice and expert alike. To explore data with Orange, one requires __no programming or in-depth mathematical knowledge__. We believe that workflow-based data science tools democratize data science by hiding complex underlying mechanics and exposing intuitive concepts. Anyone who owns data, or is motivated to peek into data, should have the means to do so.

<p align="center">
    <a href="https://orange.biolab.si/download">
    <img src="https://raw.githubusercontent.com/irgolic/orange3/README-shields/distribute/orange-example-tall.png" alt="Example Workflow">
    </a>
</p>

[Orange]: https://orange.biolab.si/


## Installing

### Easy installation

For easy installation, [Download](https://orange.biolab.si/download) the latest released Orange version from our website. To install an add-on, head to `Options -> Add-ons...` in the menu bar.

### Installing with Conda

First, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for your OS. 

Then, create a new conda environment, and install orange3:

```Shell
# Add conda-forge to your channels for access to the latest release
conda config --add channels conda-forge

# Perhaps enforce strict conda-forge priority
conda config --set channel_priority strict

# Create and activate an environment for Orange
conda create python=3 --yes --name orange3
conda activate orange3

# Install Orange
conda install orange3
```

For installation of an add-on, use:
```Shell
conda install orange3-<addon name>
```
[See specific add-on repositories for details.](https://github.com/biolab/)


### Installing with pip

We recommend using our [standalone installer](https://orange.biolab.si/download) or conda, but Orange is also installable with pip. You will need a C/C++ compiler (on Windows we suggest using Microsoft Visual Studio Build Tools).


### Installing with winget (Windows only)

To install Orange with [winget](https://docs.microsoft.com/en-us/windows/package-manager/winget/), run:

```Shell
winget install --id  UniversityofLjubljana.Orange 
```

## Running

Ensure you've activated the correct virtual environment. If following the above conda instructions:

```Shell
conda activate orange3
``` 

Run `orange-canvas` or `python3 -m Orange.canvas`. Add `--help` for a list of program options.

Starting up for the first time may take a while.


## Developing

[![GitHub Actions](https://img.shields.io/endpoint.svg?url=https%3A%2F%2Factions-badge.atrox.dev%2Fbiolab%2Forange3%2Fbadge&label=build)](https://actions-badge.atrox.dev/biolab/orange3/goto) [![codecov](https://img.shields.io/codecov/c/github/biolab/orange3)](https://codecov.io/gh/biolab/orange3) [![Contributor count](https://img.shields.io/github/contributors-anon/biolab/orange3)](https://github.com/biolab/orange3/graphs/contributors) [![Latest GitHub commit](https://img.shields.io/github/last-commit/biolab/orange3)](https://github.com/biolab/orange3/commits/master)

Want to write a widget? [Use the Orange3 example add-on template.](https://github.com/biolab/orange3-example-addon)

Want to get involved? Join us on [Discord](https://discord.gg/FWrfeXV), introduce yourself in #general! 

Take a look at our [contributing guide](https://github.com/irgolic/orange3/blob/README-shields/CONTRIBUTING.md) and [style guidelines](https://github.com/biolab/orange-widget-base/wiki/Widget-UI).

Check out our widget development [docs](https://orange-widget-base.readthedocs.io/en/latest/?badge=latest) for a comprehensive guide on writing Orange widgets.

### The Orange ecosystem

The development of core Orange is primarily split into three repositories:

[biolab/orange-canvas-core](https://www.github.com/biolab/orange-canvas-core) implements the canvas,  
[biolab/orange-widget-base](https://www.github.com/biolab/orange-widget-base) is a handy widget GUI library,  
[biolab/orange3](https://www.github.com/biolab/orange3) brings it all together and implements the base data mining toolbox.	

Additionally, add-ons implement additional widgets for more specific use cases. [Anyone can write an add-on.](https://github.com/biolab/orange3-example-addon) Some of our first-party add-ons:

- [biolab/orange3-text](https://www.github.com/biolab/orange3-text)
- [biolab/orange3-bioinformatics](https://www.github.com/biolab/orange3-bioinformatics)
- [biolab/orange3-timeseries](https://www.github.com/biolab/orange3-timeseries)    
- [biolab/orange3-single-cell](https://www.github.com/biolab/orange3-single-cell)    
- [biolab/orange3-imageanalytics](https://www.github.com/biolab/orange3-imageanalytics)    
- [biolab/orange3-educational](https://www.github.com/biolab/orange3-educational)    
- [biolab/orange3-geo](https://www.github.com/biolab/orange3-geo)    
- [biolab/orange3-associate](https://www.github.com/biolab/orange3-associate)    
- [biolab/orange3-network](https://www.github.com/biolab/orange3-network)
- [biolab/orange3-explain](https://www.github.com/biolab/orange3-explain)

### Setting up for core Orange development

First, fork the repository by pressing the fork button in the top-right corner of this page.

Set your GitHub username,

```Shell
export MY_GITHUB_USERNAME=replaceme
```

create a conda environment, clone your fork, and install it:

```Shell
conda create python=3 --yes --name orange3
conda activate orange3

git clone ssh://git@github.com/$MY_GITHUB_USERNAME/orange3

pip install -e orange3
```

Now you're ready to work with git. See GitHub's guides on [pull requests](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests), [forks](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests/working-with-forks) if you're unfamiliar. If you're having trouble, get in touch on [Discord](https://discord.gg/FWrfeXV).

#### Running

Run Orange with `python -m Orange.canvas` (after activating the conda environment).

`python -m Orange.canvas -l 2 --no-splash --no-welcome` will skip the splash screen and welcome window, and output more debug info. Use `-l 4` for more.

Add `--clear-widget-settings` to clear the widget settings before start.

To explore the dark side of the Orange, try `--style=fusion:breeze-dark`

Argument `--help` lists all available options.

To run tests, use `unittest Orange.tests Orange.widgets.tests`


### Setting up for development of all components

Should you wish to contribute Orange's base components (the widget base and the canvas), you must also clone these two repositories from Github instead of installing them as dependencies of Orange3.

First, fork all the repositories to which you want to contribute. 

Set your GitHub username,

```Shell
export MY_GITHUB_USERNAME=replaceme
```

create a conda environment, clone your forks, and install them:

```Shell
conda create python=3 --yes --name orange3
conda activate orange3

git clone ssh://git@github.com/$MY_GITHUB_USERNAME/orange-widget-base
pip install -e orange-widget-base

git clone ssh://git@github.com/$MY_GITHUB_USERNAME/orange-canvas-core
pip install -e orange-canvas-core

git clone ssh://git@github.com/$MY_GITHUB_USERNAME/orange3
pip install -e orange3

# Repeat for any add-on repositories
```

It's crucial to install `orange-base-widget` and `orange-canvas-core` before `orange3` to ensure that `orange3` will use your local versions.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at orange@biolab.si. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
Contributing
============

Thanks for taking the time to contribute to Orange!

This document outlines our guidelines and standards of
contributing to Orange. If anything is unclear, feel free
to join our [Discord server] for a chat.

Please submit contributions in accordance with the flow explained in the
[GitHub Guides].

[GitHub Guides]: https://guides.github.com/
[Discord server]: https://discord.gg/FWrfeXV


Installing for development
--------------------------
Install Orange as suggested in [README]. Then:

    git clone https://github.com/biolab/orange3.git
    cd orange3
    python setup.py develop

[README]: https://github.com/biolab/orange3/blob/master/README.md


Reporting bugs
--------------
When reporting bugs, please fill out the [issue template] to the best of your ability.
At least, include a method to reproduce the bug (if consistently
reproducible) and a screenshot (if applicable).

[issue template]: https://github.com/biolab/orange3/issues/new?assignees=&labels=bug+report&template=bug_report.md&title=

Coding style
------------

Roughly conform to [PEP-8] style guide for Python code. Whenever PEP-8 is
undefined, adhere to [Google Python Style Guide].

In addition, we add the following guidelines:

* Only ever `import *` to make objects available in another namespace,
  preferably in *\_\_init\_\_.py*. Everywhere else use explicit object
  imports.
* Use [Napoleon]-compatible (e.g. NumPy style) docstrings, preferably with
  [tests].
* When instantiating Qt widgets, pass static property values as
  [keyword args to the constructor] instead of calling separate property
  setters later. For example, do:

      view = QListView(alternatingRowColors=True,
                       selectionMode=QAbstractItemView.ExtendedSelection)

  instead of:

      view = QListView()
      view.setAlternatingRowColors(True)
      view.setSelectionMode(QAbstractItemView.ExtendedSelection)

* Each Orange widget module, or better still, each Python module (within
  reason) should have a `__name__ == '__main__'`-fenced code block that
  shows/tests the gist of that module in a user-friendly way.
* Core library objects should represent (`__repr__`) themselves in accordance
  with the following statement from [Python data model documentation]:

  > If at all possible, \[the string returned by `__repr__`\] should look like
  > a valid Python expression that could be used to recreate an object with
  > the same value (given an appropriate environment).

  To that end, use [`Orange.util.Reprable`] when possible.

Please ensure your commits pass code quality assurance by executing:

    pip install -r requirements-dev.txt
    python setup.py lint

[PEP-8]: https://www.python.org/dev/peps/pep-0008/
[Google Python Style Guide]: https://google.github.io/styleguide/pyguide.html
[Napoleon]: http://www.sphinx-doc.org/en/stable/ext/napoleon.html
[keyword args to the constructor]: http://pyqt.sourceforge.net/Docs/PyQt5/qt_properties.html
[Python data model documentation]: https://docs.python.org/3/reference/datamodel.html#object.__repr__
[`Orange.util.Reprable`]: https://github.com/biolab/orange3/search?q="class+Reprable"&type=Code


Human Interface Guidelines
--------------------------
For UI design, conform to the [OS X Human Interface Guidelines].
In a nutshell, use title case for titles, push buttons, menu titles
and menu options. Elsewhere, use sentence case. Use title case for
combo box options where the item is imperative (e.g. Initialize with Method)
and sentence case otherwise.

[OS X Human Interface Guidelines]: https://developer.apple.com/library/mac/documentation/UserExperience/Conceptual/OSXHIGuidelines/TerminologyWording.html


Testing
-------
[tests]: #tests
If you contribute new code, write [unit tests] for it in _Orange/tests_ or
_Orange/widgets/*/tests_, as appropriate. Ensure the tests pass by running:

    python setup.py test

Additionally, check that the tests for widgets pass:

    python -m unittest -v Orange.tests \
                          Orange.widgets.tests

If testing on GNU/Linux, perhaps install _xvfb_ package and prefix the above
command with `xvfb-run `.

Prefer [doctests] for public APIs. Note, we unit-test doctests with
`NORMALIZE_WHITESPACE` and `ELLIPSIS` options enabled, so you can use them
implicitly.

[unit tests]: https://en.wikipedia.org/wiki/Unit_testing
[doctests]: https://en.wikipedia.org/wiki/Doctest


Environment variables
---------------------
Set these environment variables for value-added behavior:

* `ORANGE_DEBUG=1` - general developing and debugging. Influences stuff like
  DOM Inspector in QWebView right-click menu, etc.
* `ORANGE_DEPRECATIONS_ERROR=1` - whether warnings of type
  `OrangeDeprecationWarning` should be raised as exceptions.


Commit messages
---------------
Make a separate commit for each logical change you introduce. We prefer
short commit messages with descriptive titles. For a general format see
[Commit Guidelines]. E.g.:

> io: Fix reader for XYZ file format
>
> The reader didn't work correctly in such-and-such case.

The commit title (first line) should concisely explain _WHAT_ is the change.
If the reasons for the change aren't reasonably obvious, also explain the
_WHY_ and _HOW_ in the commit body.

The commit title should start with a tag which concisely conveys what
Python package, module, or class the introduced change pertains to.

**ProTip**: Examine project's [commit history] to see examples of commit
messages most probably acceptable to that project.

[Commit Guidelines]: http://git-scm.com/book/ch5-2.html#Commit-Guidelines
[commit history]: https://github.com/biolab/orange3/commits/master


Pull requests
-------------
Implement new features in separate topic branches:

    git checkout master
    git checkout -b my-new-feature   # spin a branch off of current branch

When you are asked to make changes to your pull request, and you add the
commits that implement those changes, squash commits that fit together.

E.g., if your pull request looks like this:

    d43ef09   Some feature I made
    b803d26   reverts part of previous commit
    77d5ad3   Some other bugfix
    9e30343   Another new feature
    1d5b3bc   fix typo (in previous commit)

interactively rebase the commits onto the master branch:

    git rebase --interactive master

and mark `fixup` or `squash` the commits that are just minor patches on
previous commits (interactive rebase also allows you to reword and reorder
commits). The resulting example pull request should look clean:

    b432f18   some_module: Some feature I made
    85d5a0a   other.module: Some other bugfix
    439e303   OWSomeWidget: Another new feature

Read more [about squashing commits].

[about squashing commits]: https://www.google.com/search?q=git+squash+commits


Documentation
-------------
Documentation in located in doc folder. It is split into three parts:
data-mining-library (scripting api), development (development guides),
and visual-programming (widget help files). You can build it with:

    cd doc/<part>
    make html
    # Now open build/html/index.html to see it
<!--
If something's not right with Orange, fill out the bug report template at:
https://github.com/biolab/orange3/issues/new?assignees=&labels=bug+report&template=bug_report.md&title=

If you have an idea for a new feature, fill out the feature request template at:
https://github.com/biolab/orange3/issues/new?assignees=&labels=&template=feature_request.md&title=
-->
##### Issue
<!-- E.g. Fixes #1, Closes #2, Resolves #3, etc. -->
<!-- Or a short description, if the issue does not exist. -->


##### Description of changes


##### Includes
- [X] Code changes
- [ ] Tests
- [ ] Documentation
---
name: Bug report
about: Report a bug to help us improve
title: ''
labels: bug report
assignees: ''

---

<!-- 
Thanks for taking the time to report a bug!
If you're raising an issue about an add-on (i.e., installed via Options > Add-ons), raise an issue in the relevant add-on's issue tracker instead. See: https://github.com/biolab?q=orange3
To fix the bug, we need to be able to reproduce it. Please answer the following questions to the best of your ability.
-->


**What's wrong?**
<!-- Be specific, clear, and concise. Include screenshots if relevant. -->
<!-- If you're getting an error message, copy it, and enclose it with three backticks (```). -->





**How can we reproduce the problem?** 
<!-- Upload a zip with the .ows file and data. -->
<!-- Describe the steps (open this widget, click there, then add this...) -->





**What's your environment?**
<!-- To find your Orange version, see "Help → About → Version" or `Orange.version.full_version` in code -->
- Operating system:
- Orange version:
- How you installed Orange:
---
name: Feature request
about: Suggest an idea for Orange
title: ''
labels: ''
assignees: ''

---

<!-- 
Thanks for taking the time to submit a feature request!
For the best chance at our team considering your request, please answer the following questions to the best of your ability.
-->

**What's your use case?**
<!-- In other words, what's your pain point? -->
<!-- Is your request related to a problem, or perhaps a frustration? -->
<!-- Tell us the story that led you to write this request. -->





**What's your proposed solution?**
<!-- Be specific, clear, and concise. -->





**Are there any alternative solutions?**





# Computation of distances in Orange 3

This document describes and justifies how Orange 3 computes distances between data rows or columns from the data that can include discrete (nominal) and numeric features with missing values.

The aim of normalization is to bring all numeric features onto the same scale and on the same scale as discrete features. The meaning of *the same scale* is rather arbitrary. We gauge the normalization so that missing values have the same effect for all features and their types.

For missing values, we compute the expected difference given the probability distribution of the feature that is estimated from the data.

Two nominal values are treated as same or different, that is, the difference between them is 0 and 1.

Difference between values of two distinct nominal features does not make sense, so Orange reports an error when the user tries to compute column-wise distance in data with some non-numeric features.

## Euclidean distance

#### Normalization of numeric features

Orange 2 used to normalize by subtracting the minimum and dividing by the span (the difference between the maximum and the minimum) to bring the range of differences into interval $[0, 1]$. This however did not work since the data on which the distances were computed could include more extreme values than the training data.

Normalization in Orange 3 is based on mean and variance due to other desired effects described below. A value $x$ is normalized as

$$ x' = \frac{x - \mu}{\sqrt{2\sigma^2}},$$

where $\mu$ and $\sigma^2$ are the mean and the variance of that feature (e.g. estimated accross the column).

Normalized values thus have a mean of $\mu'=0$ and a variance of $\sigma'^2 = 1/2$.

#### Missing values of numeric features

If one value (denoted by $v$) is known and one missing, the expected difference along this dimension is

$$\int_{-\infty}^{\infty}(v - x)^2p(x)dx = \\
v^2\int_{-\infty}^{\infty}p(x)dx- 2v\int_{-\infty}^{\infty}xp(x) + \int_{-\infty}^{\infty}x^2p(x) = \\
v^2 - 2v\mu + (\sigma^2 + \mu^2) = \\
(v - \mu)^2 + \sigma^2.$$

If both values are unknown and we compute the difference between rows so that both values come from the same distirbutions, we have

$$\int_{-\infty}^{\infty}\int_{-\infty}^{\infty}(x - y)^2p(x)p(y)dxdy = \\
 \int_{-\infty}^{\infty}\int_{-\infty}^{\infty}x^2p(x)p(y)dxdy
 + \int_{-\infty}^{\infty}\int_{-\infty}^{\infty}y^2p(x)p(y)dxdy
 - 2\int_{-\infty}^{\infty}\int_{-\infty}^{\infty}xyp(x)p(y)dxdy = \\
 (\sigma^2 + \mu^2) + (\sigma^2 + \mu^2) - 2\int_{-\infty}^{\infty}xp(x)dx\int_{-\infty}^{\infty}yp(y)dxdy = \\
 (\sigma^2 + \mu^2) + (\sigma^2 + \mu^2) - 2\mu\mu = \\
 2\sigma^2.$$

When computing the difference between columns, the derivation is similar except that the two distributions are not the same. For one missing value we get

$$(v - \mu_x)^2 + \sigma_x^2$$

where $\mu_x$ and $\sigma_x$ correspond to the distribution of the unknown value. For two missing values, we get

$$\sigma_x^2 + \sigma_y^2 + \mu_x^2 + \mu_y^2 - 2\mu_x\mu_y = \\
\sigma_x^2 + \sigma_y^2 + (\mu_x - \mu_y)^ 2.$$

For normalized data, the difference between $v$ and unknown value is $(v' - \mu)^2 + \sigma^2 = v'^2 + 1/2$. The difference between two missing values is $2\sigma^2 = 1$ (or $\sigma_x^2 + \sigma_y^2 + (\mu_x - \mu_y)^2 = 1$).

#### Missing values of discrete features

The difference between a known and a missing value is

$$\sum_x \mbox{I}_{v\ne x}^2p(x) = 1 - p(x).$$

The difference between two missing values is

$$\sum_x\sum_y \mbox{I}_{y\ne x}^2p(x)p(y) = 1 - \sum_x p(x)^2.$$

This is the Gini index. Also, if the number of values goes to infinity and the distribution towards the uniform, the difference goes towards 1, which brings it, in some sense, to the same scale as continuous features.

This case assumes that $x$ and $y$ come from the same distribution. The case when these are missing values of two distinct discrete features is not covered since Orange does not support such distances (see the introduction).

## Manhattan distance

The derivation here may be more difficult, so we do it by analogy with Euclidean distances and hope for the best.

We use the median ($m$) and median of absolute distances to the median (mad) ($a$) instead of the mean and the variance.

Normalization of numeric features: $x' = (x - m)\,/\,(2a)$

#### Missing values of numeric features

**Between a known and a missing**: $|v - m| + a.$

**Between two unknowns**: $2a$ (same features), $a_x + a_y$ (different features).

**For normalized data**: $|v'| + 1/2$ (one unknown), $1$ (both unknown).

#### Missing values of discrete features

Same as for Euclidean because $I_{v\ne x}^2 = I_{v\ne x}$.


## Cosine distance

The cosine similarity is normalized, so we do not normalize the data.

Cosine distance treats discrete features differently from the Euclidean and Manhattan distance. The latter needs to consider only whether two values of a discrete feature is different or not, while the cosine distance normalizes by dividing by vector lengths. For this, we need the notion of absolute magnitude of a (single) discrete value -- as compared to some "base value".

For this reason, cosine distance treats discrete attributes as boolean, that is, all non-zero values are treated as 1. This may be incorrect in some scenarios, especially those in which cosine distance is inappropriate anyway. How to (and whether to) use cosine distances on discrete data is up to the user.

#### Distances between rows

For a continuous variable $x$, product of $v$ and a missing value is computed as

$$\int_{-\infty}^{\infty}vxp(x)dx = v\mu_x$$

The product of two unknown values is

$$\int_{-\infty}^{\infty}xp(x)yp(y)dxdy=\mu_x^2$$

For discrete values, we compute the probabilities $p(x=0)$ and $p_x(x\ne 0)$. The product of known value $v$ and a missing value is 0 if v=0 and $p(x \ne 1)$ otherwise. The product of two missing values is $p(x\ne 1)^2$.

When computing the absolute value of a row, a missing value of continuous variable theoretically contributes

$$\int_{-infty}^{\infty}x^2p(x)dx = \mu_x^2 + \sigma_x^2$$

However, since we essentially impute the mean in the dot product, the actual contribution of the missing value is $\mu_x^2$. We therefore use this value, which also simplifies the computation which is reduced to simple imputation of means.

A missing value of discrete variable contributes

$$1\cdot 1\;p(x\ne 0) = p(x\ne 0)$$


#### Distances between columns

The product of a known value $v$ and a missing value of $x$ is $v\mu_x$. The contribution of a missing value to the absolute value of the column is $\mu_x^2+\sigma_x^2$. All derivations are same as above.

## Jaccard distance

Jaccard index, whose computation is described below, is a measure of similarity. The distance is computed by subtracting the similarity from 1.

Let $p(A_i)$ be the probability (computed from the training data) that a random data instance belongs to set $A_i$, i.e., have a non-zero value for feature $A_i$.

### Similarity between rows (instances)

Let $M$ and $N$ be two data instances. $M$ and $N$ can belong to $A_i$ (or not). The Jaccard similarity between $M$ and $N$ is the number of the common sets to which $M$ and $N$ belong, divided by the number of sets with either $M$ or $N$ (or both).

Let $\mbox{I}_M$ be 1 if $M\in A_i$ and 0 otherwise. Similarly, $\mbox{I}_{M'}$ will indicate that $M\not\in A_i$ and $\mbox{I}_{M?}$ will indicate that it is unknown whether $M$ belongs to $A_i$ or not. We will also use conjuctions and disjunctions by adding more indices to $\mbox{I}$; e.g. $\mbox{I}_{M\wedge N'}$ indicates that $M$ belongs to $A_i$ and $N$ does not. $A_i$ is omitted for clarity as it can be deduced from the context.

$$Jaccard(M, N) = \frac{\sum_i\mbox{I}_{M\wedge N}}{\sum_i\mbox{I}_{M\vee N}}$$

Consider that $\mbox{I}_{M\wedge N} = \mbox{I}_M \mbox{I}_N$ and $\mbox{I}_{M\vee N} = \max(\mbox{I}_M, \mbox{I}_N)$. If the data whether $M\in A_i$ or $N\in A_i$ is missing, we replace indicator function with the probability. In the denominator we add a few terms to $\mbox{I}_{M\wedge N}$

$$\mbox{I}_{M\wedge N}
  + p(A_i)\mbox{I}_{M\wedge N?}
  + p(A_i)\mbox{I}_{M?\wedge N}
  + p(A_i)^2\mbox{I}_{M?\wedge N?},$$

and in the nominator we get

$$\mbox{I}_{M\vee N} + p(A_i)\mbox{I}_{M'\wedge N?} + p(A_i)\mbox{I}_{M?\wedge N'} + \left(1 - (1 - p(A_i)^2\right)\mbox{I}_{I?\wedge J?}$$

Note that the denominator counts cases $\mbox{I}_{M'\wedge N?}$ and not $\mbox{I}_{N?}$, since those for which $M\in A_i$ are already covered in $\mbox{I}_{M\vee N}$. The last term refers to the probability that at least one (that is, not none) of the two instances is in $A_i$.

### Similarity between columns

$\mbox{I}_{i}$ will now denote that a data instance $M$ belongs to $A_i$, $\mbox{I}_{i'}$ will denote it does not, and $\mbox{I}_{i?}$ will denote that it is unknown whether $M$ belongs to $A_i$.

Without considering missing data, Jaccard index between two columns is

$$Jaccard(A_i, A_j) = \frac{\sum_M\mbox{I}_{i\wedge j}}{\sum_M\mbox{I}_{i\vee j}}$$

By the same reasoning as above, the denominator becomes

$$\mbox{I}_{i\wedge j}
  + p(A_j)\mbox{I}_{i\wedge j?}
  + p(A_i)\mbox{I}_{i?\wedge j}
  + p(A_i)p(A_j)\mbox{I}_{i?\wedge j?},$$

and the nominator is

$$\mbox{I}_{i\vee j} + p(A_j)\mbox{I}_{i'\wedge j?} + p(A_i)\mbox{I}_{i?\wedge j'} + \left(1 - [1 - p(A_i)][1 - p(A_j)]\right)\mbox{I}_{I?\wedge J?}.$$

The sums runs over instances, $M$, so the actual implementation can work by counting the cases and multipying by probabilities at the end. Let $N_c$ represent the number of cases that match condition $c$, i.e. $N_c = \sum_M\mbox{I}_c$. Then

$$Jaccard(A_i, A_j) = \frac{
      N_{i\wedge j}
    + p(A_j)N_{i\wedge j?}
    + p(A_i)N_{i?\wedge j}
    + p(A_i) p(A_j)N_{i?\wedge j?}
    }{
      N_{i\vee j}
    + p(A_j)N_{i'\wedge j?}
    + p(A_i)N_{i?\wedge j'}
    + \left(1 - [1 - p(A_i)][1 - p(A_j]\right)N_{i?\wedge j?}
    }$$
This directory contains IPython Notebooks (http://ipython.org/notebook.html)
with a selection of Orange tutorials.

The current list includes:
* [Custom learners in Orange](https://github.com/biolab/orange3/blob/master/tutorials/learners.ipynb)
# Learners as Scorers

Certain learners can be used as feature scorers in Orange. Here's a quick example with [Random Forest](../widgets/model/randomforest.md).

We are using the *iris* data for the example. Connect [File](../widgets/data/file.md) with [Rank](../widgets/data/rank.md). Then connect **Random Forest** to Rank. Random Forest will be used as a Scorer in this case. Rank will use Random Forest's feature importance to rank the attributes.

![](scoring-with-RF.png)

Passing additional scorers works for both, classification and regression:

- [Logistic Regression](../widgets/model/logisticregression.md) (classification) / [Linear Regression](../widgets/model/linearregression.md) (regression)
- [Stochastic Gradient Descent](../widgets/model/stochasticgradient.md)
- [Gradient Boosting](../widgets/model/gradientboosting.md)
- Random Forest
# Loading your Data

Orange comes with its [own data format](https://docs.biolab.si/3/data-mining-library/tutorial/data.html#data-input), but can also handle native Excel, comma- or tab-delimited data files. The input data set is usually a table, with data instances (samples) in rows and data attributes in columns. Attributes can be of different *types* (numeric, categorical, datetime, and text) and have assigned *roles* (input features, meta attributes, and class). Data attribute type and role can be provided in the data table header. They can also be changed in the [File](../widgets/data/file.md) widget, while data role can also be modified with [Select Columns](../widgets/data/selectcolumns.md) widget.

### In a Nutshell

- Orange can import any comma- or tab-delimited data file, or Excel's native files or Google Sheets document. Use [File](../widgets/data/file.md) widget to load the data and, if needed, define the class and meta attributes.
- Types and roles can be set in the File widget.
- Attribute names in the column header can be preceded with a label followed by a hash. Use c for class and m for meta attribute, i to ignore a column, w for weights column, and C, D, T, S for continuous, discrete, time, and string attribute types. Examples: C\#mph, mS\#name, i\#dummy.
- An alternative to the hash notation is Orange's native format with three header rows: the first with attribute names, the second specifying the type (**continuous**, **discrete**, **time**, or **string**), and the third proving information on the attribute role (**class**, **meta**, **weight** or **ignore**).

## Data from Excel

Here is an example dataset ([sample.xlsx](http://file.biolab.si/datasets/sample.xlsx)) as entered in Excel:

![](spreadsheet1.png)

The file contains a header row, eight data instances (rows) and seven data attributes (columns). Empty cells in the table denote missing data entries. Rows represent genes; their function (class) is provided in the first column and their name in the second. The remaining columns store measurements that characterize each gene. With this data, we could, say, develop a classifier that would predict gene function from its characteristic measurements.

Let us start with a simple workflow that reads the data and displays it in a table:

![](file-data-table-workflow.png)

To load the data, open the File widget (double click on the icon of the widget), click on the file browser icon ("...") and locate the downloaded file (called [sample.xlsx](http://file.biolab.si/datasets/sample.xlsx)) on your disk:

![](File.png)

### File Widget: Setting the Attribute Type and Role

The **File** widget sends the data to the **Data Table**. Double click the **Data Table** to see its contents:

![](table-widget.png)

Orange correctly assumed that a column with gene names is meta information, which is displayed in the **Data Table** in columns shaded with light-brown. It has not guessed that *function*, the first non-meta column in our data file, is a class column. To correct this in Orange, we can adjust attribute role in the column display of File widget (below). Double-click the *feature* label in the *function* row and select *target* instead. This will set *function* attribute as our target (class) variable.

![](File-set-feature-kind.png)

You can also change attribute type from nominal to numeric, from string to datetime, and so on. Naturally, data values have to suit the specified attribute type. Datetime accepts only values in [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601) format, e.g. 2016-01-01 16:16:01. Orange would also assume the attribute is numeric if it has several different values, else it would be considered nominal. All other types are considered strings and are as such automatically categorized as meta attributes.

Change of attribute roles and types should be confirmed by clicking the **Apply** button.

### Select Columns: Setting the Attribute Role

Another way to set the data role is to feed the data to the [Select Columns](../widgets/data/selectcolumns.md) widget:

![](select-columns-schema.png)

Opening [Select Columns](../widgets/data/selectcolumns.md) reveals Orange's classification of attributes. We would like all of our continuous attributes to be data features, gene function to be our target variable and gene names considered as meta attributes. We can obtain this by dragging the attribute names around the boxes in **Select Columns**:

![](select-columns-start.png)

To correctly reassign attribute types, drag attribute named *function* to a **Class** box, and attribute named *gene* to a **Meta Attribute** box. The [Select Columns](../widgets/data/selectcolumns.md) widget should now look like this:

![](select-columns-reassigned.png)

Change of attribute types in *Select Columns* widget should be confirmed by clicking the **Apply** button. The data from this widget is fed into [Data Table](../widgets/data/datatable.md) that now renders the data just the way we intended:

![](data-table-with-class1.png)

We could also define the domain for this dataset in a different way. Say, we could make the dataset ready for regression, and use *heat 0* as a continuous class variable, keep gene function and name as meta variables, and remove *heat 10* and *heat 20* from the dataset:

![](select-columns-regression.png)

By setting the attributes as above, the rendering of the data in the
Data Table widget gives the following output:

![](data-table-regression1.png)

## Header with Attribute Type Information

Consider again the [sample.xlsx](http://file.biolab.si/datasets/sample.xlsx) dataset. This time we will augment the names of the attributes with prefixes that define attribute type (continuous, discrete, time, string) and role (class or meta attribute). Prefixes are separated from the attribute name with a hash sign ("\#"). Prefixes for attribute roles are:

- c: class attribute
- m: meta attribute
- i: ignore the attribute
- w: instance weights

and for the type:

- C: Continuous
- D: Discrete
- T: Time
- S: String

This is how the header with augmented attribute names looks like in Excel ([sample-head.xlsx](http://file.biolab.si/datasets/sample-head.xlsx)):

![](spreadsheet-simple-head1.png)

We can again use a **File** widget to load this dataset and then render it in the **Data Table**:

![](select-cols-simplified-header.png)

Notice that the attributes we have ignored (label "i" in the attribute name) are not present in the dataset.

## Three-Row Header Format

Orange's legacy native data format is a tab-delimited text file with three header rows. The first row lists the attribute names, the second row defines their type (continuous, discrete, time and string, or abbreviated c, d, t, and s), and the third row an optional role (class, meta, weight, or ignore). Here is an example:

![](excel-with-tab1.png)

Data from Google Sheets
-----------------------

Orange can read data from Google Sheets, as long as it conforms to the data presentation rules we have presented above. In Google Sheets, copy the shareable link (Share button, then Get shareable link) and paste it in the *Data File / URL* box of the File widget. For a taste, here's one such link you can use: [http://bit.ly/1J12Tdp](http://bit.ly/1J12Tdp), and the way we have entered it in the **File** widget:

![](File-Google-Sheet.png)

## Data from LibreOffice

If you are using LibreOffice, simply save your files in Excel (.xlsx) format (available from the drop-down menu under *Save As Type*).

![](saving-tab-delimited-files.png)

## Datetime Format

To avoid ambiguity, Orange supports date and/or time formatted in one of the [ISO 8601](https://en.wikipedia.org/wiki/ISO_8601) formats. For example, the following values are all valid:

    2016
    2016-12-27
    2016-12-27 14:20:51
    16:20
# Building Workflows

The core principle of Orange is visual programming, which means each analytical step in contained within a widget. Widgets are placed on the canvas and connected into an analytical workflow, which is executed from left to right. Orange never passes data backwards.

## Simple workflow

Let us start with a simple workflow. We will load the data with the File widget, say the famous *Iris* data set. Right-click on the canvas. A menu will appear. Start typing "File", then press Enter to confirm the selection. [File](../widgets/data/file.md) widget will be placed on the canvas.

![](file.gif)

**File** widget has an "ear" on its right side – this is the output of the widget. Click on the "ear" and drag a connection out of it. Upon releasing the connection, a menu will appear. Start typing the name of the widget to connect with the File widget, say Data Table. Select the widget and press enter. The widget is added to the canvas.

![](file-datatable.gif)

This is a simple workflow. The File widget loads the data and sends it to the output. Data Table receives the data and displays it in a table. Please note that Data Table is a viewer and passes onwards only the selection. The data is always available at the source - in the File widget.

![](DataTable-wrong.png)

## Workflows with subsets

Visualizations in Orange are interactive, which means the user can select data instances from the plot and pass them downstream. Let us look at two examples with subsets.

### Selecting subsets

Place **File** widget on the canvas. Then connect [Scatter Plot](../widgets/visualize/scatterplot.md) to it. Click and drag a rectangle around a subset of points. Connect [Data Table](../widgets/data/datatable.md) to Scatter Plot. Data Table will show selected points.

![](subset-selection.gif)

### Highlighting workflows

Place **File** widget on the canvas. Then connect **Scatter Plot** to it and a **Data Table**. Connect Data Table to Scatter Plot. Select a subset of points from the Data Table. Scatter Plot will highlight selected points.

![](subset-highlight.gif)

## Workflows with models

Predictive models are evaluated in [Test and Score](../widgets/evaluate/testandscore.md) widget, while predictions on new data are done in [Predictions](../widgets/evaluate/predictions.md). Test and Score accepts several inputs: data (data set for evaluating models), learners (algorithms to use for training the model), and an optional preprocessor (for normalization or feature selection).

![](prediction-workflow.png)

For prediction, the training data is first passed to the model. Once the model is trained, it is passed to **Predictions**. The Predictions widget also needs data to predict on, which are passed as a second input.

![](prediction-workflow2.png)
Test and Score
==============

Tests learning algorithms on data.

**Inputs**

- Data: input dataset
- Test Data: separate data for testing
- Learner: learning algorithm(s)

**Outputs**

- Evaluation Results: results of testing classification algorithms

The widget tests learning algorithms. Different sampling schemes are available, including using separate test data. The widget does two things. First, it shows a table with different classifier performance measures, such as [classification accuracy](https://en.wikipedia.org/wiki/Accuracy_and_precision) and [area under the curve](https://en.wikipedia.org/wiki/Receiver_operating_characteristic#Area_under_the_curve). Second, it outputs evaluation results, which can be used by other widgets for analyzing the performance of classifiers, such as [ROC Analysis](../evaluate/rocanalysis.md) or [Confusion Matrix](../evaluate/confusionmatrix.md).

The *Learner* signal has an uncommon property: it can be connected to more than one widget to test multiple learners with the same procedures.

![](images/TestAndScore-stamped.png)

1. The widget supports various sampling methods.
   - [Cross-validation](https://en.wikipedia.org/wiki/Cross-validation_\(statistics\)) splits the data into a given number of folds (usually 5 or 10). The algorithm is tested by holding out examples from one fold at a time; the model is induced from other folds and examples from the held out fold are classified. This is repeated for all the folds.
   - **Cross validation by feature** performs cross-validation but folds are defined by the selected categorical feature from meta-features.   
   - **Random sampling** randomly splits the data into the training and testing set in the given proportion (e.g. 70:30); the whole procedure is repeated for a specified number of times.
   - **Leave-one-out** is similar, but it holds out one instance at a time, inducing the model from all others and then classifying the held out instances. This method is obviously very stable, reliable... and very slow.
   - **Test on train data** uses the whole dataset for training and then for testing. This method practically always gives wrong results.
   - **Test on test data**: the above methods use the data from *Data* signal only. To input another dataset with testing examples (for instance from another file or some data selected in another widget), we select *Separate Test Data* signal in the communication channel and select Test on test data.
2. For classification, *Target class* can be selected at the bottom of the widget. When *Target class* is (Average over classes), methods return scores that are weighted averages over all classes. For example, in case of the classifier with 3 classes, scores are computed for class 1 as a target class, class 2 as a target class, and class 3 as a target class. Those scores are averaged with weights based on the class size to retrieve the final score.
3. The widget will compute a number of performance statistics. A few are shown by default. To see others, right-click on the header and select the desired statistic.
   - Classification
   ![](images/TestAndScore-Classification.png)
        - [Area under ROC](http://gim.unmc.edu/dxtests/roc3.htm) is the area under the receiver-operating curve.
        - [Classification accuracy](https://en.wikipedia.org/wiki/Accuracy_and_precision) is the proportion of correctly classified examples.
        - [F-1](https://en.wikipedia.org/wiki/F1_score) is a weighted harmonic mean of precision and recall (see below).
        - [Precision](https://en.wikipedia.org/wiki/Precision_and_recall) is the proportion of true positives among instances classified as positive, e.g. the proportion of *Iris virginica* correctly identified as Iris virginica.
        - [Recall](https://en.wikipedia.org/wiki/Precision_and_recall) is the proportion of true positives among all positive instances in the data, e.g. the number of sick among all diagnosed as sick.
        - [Specificity](https://en.wikipedia.org/wiki/Sensitivity_and_specificity) is the proportion of true negatives among all negative instances, e.g. the number of non-sick among all diagnosed as non-sick.
        - [LogLoss](https://en.wikipedia.org/wiki/Cross_entropy) or cross-entropy loss takes into account the uncertainty of your prediction based on how much it varies from the actual label. 
        - Train time - cumulative time in seconds used for training models.
        - Test time - cumulative time in seconds used for testing models.
   - Regression
   ![](images/TestAndScore-Regression.png)
      - [MSE](https://en.wikipedia.org/wiki/Mean_squared_error) measures the average of the squares of the errors or deviations (the difference between the estimator and what is estimated).
      - [RMSE](https://en.wikipedia.org/wiki/Root_mean_square) is the square root of the arithmetic mean of the squares of a set of numbers (a measure of imperfection of the fit of the estimator to the data)
      - [MAE](<https://en.wikipedia.org/wiki/Mean_absolute_error>) is used to measure how close forecasts or predictions are to eventual outcomes.
      - [R2](<https://en.wikipedia.org/wiki/Coefficient_of_determination>) is interpreted as the proportion of the variance in the dependent variable that is predictable from the independent variable.
      - [CVRMSE](https://en.wikipedia.org/wiki/Root-mean-square_deviation) is RMSE normalized by the mean value of actual values.
      - Train time - cumulative time in seconds used for training models.
      - Test time - cumulative time in seconds used for testing models.
4. Choose the score for pairwise comparison of models and the region of practical equivalence (ROPE), in which differences are considered negligible.
5. Pairwise comparison of models using the selected score (available only for cross-validation). The number in the table gives the probability that the model corresponding to the row has a higher score than the model corresponding to the column. What the higher score means depends on the metric: a higher score can either mean a model is better (for example, CA or AUC) or the opposite (for example, RMSE). If negligible difference is enabled, the smaller number below shows the probability that the difference between the pair is negligible. The test is based on the [Bayesian interpretation of the t-test](https://link.springer.com/article/10.1007/s10994-015-5486-z) ([shorter introduction](https://baycomp.readthedocs.io/en/latest/introduction.html)).
6. Get help and produce a report.

Preprocessing for predictive modeling
--------------------------------------

When building predictive models, one has to be careful about how to preprocess the data. There are two possible ways to do it in Orange, each slightly different:

1. Connect [Preprocess](../data/preprocess.md) to the learner. This will override the default preprocessing pipeline for the learner and apply only custom preprocessing pipeline (default preprocessing steps are described in each learner's documentation). The procedure might lead to errors within the learner.

   ![](../data/images/Preprocess-Models1.png)

2. Connect **Preprocess** to Test and Score. This will apply the preprocessors to each batch within cross-validation. Then the learner's preprocessors will be applied to the preprocessed subset.

   ![](../data/images/Preprocess-Models2.png)

Finally, there's a wrong way to do it. Connecting **Preprocess** directly to the original data and outputting preprocessed data set will likely overfit the model. Don't do it.

   ![](../data/images/Preprocess-Models3.png)

Example
-------

In a typical use of the widget, we give it a dataset and a few learning algorithms and we observe their performance in the table inside the **Test & Score** widget and in the [ROC](../evaluate/rocanalysis.md). The data is often preprocessed before testing; in this case we did some manual feature selection ([Select Columns](../data/selectcolumns.md) widget) on *Titanic* dataset, where we want to know only the sex and status of the survived and omit the age.

In the bottom table, we have a pairwise comparison of models. We selected that comparison is based on the _area under ROC curve_ statistic. The number in the table gives the probability that the model corresponding to the row is better than the model corresponding to the column. We can, for example, see that probability for the tree to be better than SVM is almost one, and the probability that tree is better than Naive Bayes is 0.001. Smaller numbers in the table are probabilities that the difference between the pair is negligible based on the negligible threshold 0.1.

![](images/TestAndScore-Example.png)

Another example of using this widget is presented in the documentation for the [Confusion Matrix](../evaluate/confusionmatrix.md) widget.
Confusion Matrix
================

Shows proportions between the predicted and actual class.

**Inputs**

- Evaluation results: results of testing classification algorithms

**Outputs**

- Selected Data: data subset selected from confusion matrix
- Data: data with the additional information on whether a data instance was selected

The [Confusion Matrix](https://en.wikipedia.org/wiki/Confusion_matrix) gives the number/proportion of instances between the predicted and actual class. The selection of the elements in the matrix feeds the corresponding instances into the output signal. This way, one can observe which specific instances were misclassified and how.

The widget usually gets the evaluation results from [Test & Score](../evaluate/testandscore.md); an example of the schema is shown below.

![](images/ConfusionMatrix-stamped.png)

1. When evaluation results contain data on multiple learning algorithms, we have to choose one in the *Learners* box.
   The snapshot shows the confusion matrix for [Tree](../model/tree.md) and [Naive Bayesian](../model/naivebayes.md) models trained and tested on the *iris* data. The right-hand side of the widget contains the matrix for the naive Bayesian model (since this model is selected on the left). Each row corresponds to a correct class, while columns represent the predicted classes. For instance, four instances of *Iris-versicolor* were misclassified as *Iris-virginica*. The rightmost column gives the number of instances from each class (there are 50 irises of each of the three classes) and the bottom row gives the number of instances classified into each class (e.g., 48 instances were classified into virginica).
2. In *Show*, we select what data we would like to see in the matrix.
   - **Number of instances** shows correctly and incorrectly classified instances numerically.
   - **Proportions of predicted** shows how many instances classified as, say, *Iris-versicolor* are in which true class; in the table we can read the 0% of them are actually setosae, 88.5% of those classified as versicolor are versicolors, and 7.7% are virginicae.
   - **Proportions of actual** shows the opposite relation: of all true versicolors, 92% were classified as versicolors and 8% as virginicae.
   ![](images/ConfusionMatrix-propTrue.png)
3. In *Select*, you can choose the desired output.
   - **Correct** sends all correctly classified instances to the output by selecting the diagonal of the matrix.
   - **Misclassified** selects the misclassified instances.
   - **None** annuls the selection.
   As mentioned before, one can also select individual cells of the table to select specific kinds of misclassified instances (e.g. the versicolors classified as virginicae).
4. When sending selected instances, the widget can add new attributes, such as predicted classes or their probabilities, if the corresponding options *Predictions* and/or *Probabilities* are checked.
5. The widget outputs every change if *Send Automatically* is ticked. If not, the user will need to click *Send Selected* to commit the changes.
6. Produce a report.

Example
-------

The following workflow demonstrates what this widget can be used for.

![](images/ConfusionMatrix-Schema.png)

[Test & Score](../evaluate/testandscore.md) gets the data from [File](../data/file.md) and two learning algorithms from [Naive Bayes](../model/naivebayes.md) and [Tree](../model/tree.md). It performs cross-validation or some other train-and-test procedures to get class predictions by both algorithms for all (or some) data instances. The test results are fed into the **Confusion Matrix**, where we can observe how many instances were misclassified and in which way.

In the output, we used [Data Table](../data/datatable.md) to show the instances we selected in the confusion matrix. If we, for instance, click *Misclassified*, the table will contain all instances which were misclassified by the selected method.

The [Scatter Plot](../visualize/scatterplot.md) gets two sets of data. From the [File](../data/file.md) widget it gets the complete data, while the confusion matrix sends only the selected data, misclassifications for instance. The scatter plot will show all the data, with bold symbols representing the selected data.

![](images/ConfusionMatrix-Example.png)
Calibration Plot
================

Shows the match between classifiers' probability predictions and actual class probabilities.

**Inputs**

- Evaluation Results: results of testing classification algorithms

The [Calibration Plot](https://en.wikipedia.org/wiki/Calibration_curve)plots class probabilities against those predicted by the classifier(s).

![](images/CalibrationPlot-stamped.png)

1. Select the desired target class from the drop down menu.
2. Choose which classifiers to plot. The diagonal represents optimal behavior; the closer the classifier's curve gets, the more accurate its prediction probabilities are. Thus we would use this widget to see whether a classifier is overly optimistic (gives predominantly positive results) or pessimistic (gives predominantly negative results).
3. If *Show rug* is enabled, ticks are displayed at the bottom and the top of the graph, which represent negative and positive examples respectively. Their position corresponds to the classifier's probability prediction and the color shows the classifier. At the bottom of the graph, the points to the left are those which are (correctly) assigned a low probability of the target class, and those to the right are incorrectly assigned high probabilities. At the top of the graph, the instances to the right are correctly assigned high probabilities and vice versa.
4. Press *Save Image* if you want to save the created image to your computer in a .svg or .png format.
5. Produce a report.

Example
-------

At the moment, the only widget which gives the right type of signal needed by the **Calibration Plot** is [Test & Score](../evaluate/testandscore.md). The Calibration Plot will hence always follow Test & Score and, since it has no outputs, no other widgets follow it.

Here is a typical example, where we compare three classifiers (namely [Naive Bayes](../model/naivebayes.md), [Tree](../model/tree.md) and [Constant](../model/constant.md)) and input them into [Test & Score](../evaluate/testandscore.md). We used the *Titanic* dataset. Test & Score then displays evaluation results for each classifier. Then we draw **Calibration Plot** and [ROC Analysis](../evaluate/rocanalysis.md) widgets from Test & Score to further analyze the performance of classifiers. **Calibration Plot** enables you to see prediction accuracy of class probabilities in a plot.

![](images/CalibrationPlot-example.png)
Lift Curve
==========

Measures the performance of a chosen classifier against a random classifier.

**Inputs**

- Evaluation Results: results of testing classification algorithms

The **Lift curve** shows the curves for analysing the proportion of true positive data instances in relation to the classifier's threshold or the number of instances that we classify as positive.

Cumulative gains chart shows the proportion of true positive instances (for example, the number of clients who accept the offer) as a function of the number of positive instances (the number of clients contacted), assuming the the instances are ordered according to the model's probability of being positive (e.g. ranking of clients).

![](images/LiftCurve-cumulative-gain.png)

Lift curve shows the ratio between the proportion of true positive instances in the selection and the proportion of customers contacted. See [a tutorial for more details](https://medium.com/analytics-vidhya/understanding-lift-curve-b674d21e426).

![](images/LiftCurve-stamped.png)

1. Choose the desired *Target class*. The default is chosen alphabetically.
2. Choose whether to observe lift curve or cumulative gains.
3. If test results contain more than one classifier, the user can choose which curves she or he wants to see plotted. Click on a classifier to select or deselect the curve.
4. *Show lift convex hull* plots a convex hull over lift curves for all classifiers (yellow curve). The curve shows the optimal classifier (or combination thereof) for each desired lift or cumulative gain.
5. Press *Save Image* to save the created image in a .svg or .png format.
6. Produce a report.
7. A plot with **lift** or **cumulative gain** vs. **positive rate**. The dashed line represents the behavior of a random classifier.


Example
-------

The widgets that provide the right type of the signal needed by the **Lift Curve** (evaluation data) are [Test & Score](../evaluate/testandscore.md) and [Predictions](../evaluate/predictions.md).

In the example below, we observe the lift curve and cumulative gain for the bank marketing data, where the classification goal is to predict whether the client will accept a term deposit offer based on his age, job, education, marital status and similar data. The data set is available in the Datasets widget. We run the learning algorithms in the Test and Score widget and send the results to Lift Curve to see their performance against a random model. Of the two algorithms tested, logistic regression outperforms the naive Bayesian classifier. The curve tells us that by picking the first 20 % of clients as ranked by the model, we are going to hit four times more positive instances than by selecting a random sample with 20 % of clients.
Predictions
===========

Shows models' predictions on the data.

**Inputs**

- Data: input dataset
- Predictors: predictors to be used on the data

**Outputs**

- Predictions: data with added predictions
- Evaluation Results: results of testing classification algorithms

The widget receives a dataset and one or more predictors (predictive models, not learning algorithms - see the example below). It outputs the data and the predictions.

![](images/Predictions-stamped.png)

1. Information on the input, namely the number of instances to predict, the number of predictors and the task (classification or regression). If you have sorted the data table by attribute and you wish to see the original view, press *Restore Original Order*.
2. You can select the options for classification. If *Predicted class* is ticked, the view provides information on predicted class. If *Predicted probabilities for* is ticked, the view provides information on probabilities predicted by the classifier(s). You can also select the predicted class displayed in the view. The option *Draw distribution bars* provides a visualization of probabilities.
3. By ticking the *Show full dataset*, you can view the entire data table (otherwise only class variable will be shown).
4. Select the desired output.
5. Predictions.

The widget show the probabilities and final decisions of [predictive models](https://en.wikipedia.org/wiki/Predictive_modelling). The output of the widget is another dataset, where predictions are appended as new meta attributes. You can select which features you wish to output (original data, predictions, probabilities). The result can be observed in a [Data Table](../data/datatable.md). If the predicted data includes true class values, the result of prediction can also be observed in a [Confusion Matrix](../evaluate/confusionmatrix.md).

Examples
--------

In the first example, we will use *Attrition - Train* data from the [Datasets](../data/datasets.md) widget. This is a data on attrition of employees. In other words, we wish to know whether a certain employee will resign from the job or not. We will construct a predictive model with the [Tree](../model/tree.md) widget and observe probabilities in **Predictions**.

For predictions we need both the training data, which we have loaded in the first **Datasets** widget and the data to predict, which we will load in another [Datasets](../data/datasets.md) widget. We will use *Attrition - Predict* data this time. Connect the second data set to **Predictions**. Now we can see predictions for the three data instances from the second data set.

The [Tree](../model/tree.md) model predicts none of the employees will leave the company. You can try other model and see if predictions change. Or test the predictive scores first in the [Test & Score](../evaluate/testandscore.md) widget.

![](images/Predictions-Example1.png)

In the second example, we will see how to properly use [Preprocess](../data/preprocess.md) with **Predictions** or [Test & Score](../evaluate/testandscore.md).

This time we are using the *heart disease.tab* data from the [File](../data/file.md) widget. You can access the data through the dropdown menu. This is a dataset with 303 patients that came to the doctor suffering from a chest pain. After the tests were done, some patients were found to have diameter narrowing and others did not (this is our class variable).

The heart disease data have some missing values and we wish to account for that. First, we will split the data set into train and test data with [Data Sampler](../data/datasampler.md).

Then we will send the *Data Sample* into [Preprocess](../data/preprocess.md). We will use *Impute Missing Values*, but you can try any combination of preprocessors on your data. We will send preprocessed data to [Logistic Regression](../model/logisticregression.md) and the constructed model to **Predictions**.

Finally, **Predictions** also needs the data to predict on. We will use the output of [Data Sampler](../data/datasampler.md) for prediction, but this time not the *Data Sample*, but the *Remaining Data*, this is the data that wasn't used for training the model.

Notice how we send the remaining data directly to **Predictions** without applying any preprocessing. This is because Orange handles preprocessing on new data internally to prevent any errors in the model construction. The exact same preprocessor that was used on the training data will be used for predictions. The same process applies to [Test & Score](../evaluate/testandscore.md).

![](images/Predictions-Example2.png)
ROC Analysis
============

Plots a true positive rate against a false positive rate of a test.

**Inputs**

- Evaluation Results: results of testing classification algorithms

The widget shows ROC curves for the tested models and the corresponding convex hull. It serves as a mean of comparison between classification models. The curve plots a false positive rate on an x-axis (1-specificity; probability that target=1 when true value=0) against a true positive rate on a y-axis (sensitivity; probability that target=1 when true value=1). The closer the curve follows the left-hand border and then the top border of the ROC space, the more accurate the classifier. Given the costs of false positives and false negatives, the widget can also determine the optimal classifier and threshold.

![](images/ROCAnalysis-basic-stamped.png)

1. Choose the desired *Target Class*. The default class is chosen alphabetically.
2. If test results contain more than one classifier, the user can choose which curves she or he wants to see plotted. Click on a classifier to select or deselect it.
3. When the data comes from multiple iterations of training and testing, such as k-fold cross validation, the results can be (and usually are) averaged.
   ![](images/ROC-Comparison.png)
   The averaging options are:
     - **Merge predictions from folds** (top left), which treats all the test data as if they came from a single iteration
     - **Mean TP rate** (top right) averages the curves vertically, showing the corresponding confidence intervals
     - **Mean TP and FP at threshold** (bottom left) traverses over threshold, averages the positions of curves and shows horizontal and vertical confidence intervals
     - **Show individual curves** (bottom right) does not average but prints all the curves instead
4. Option *Show convex ROC curves* refers to convex curves over each individual classifier (the thin lines positioned over curves). *Show ROC convex hull* plots a convex hull combining all classifiers (the gray area below the curves). Plotting both types of convex curves makes sense since selecting a threshold in a concave part of the curve cannot yield optimal results, disregarding the cost matrix. Besides, it is possible to reach any point on the convex curve by combining the classifiers represented by the points on the border of the concave region.
  ![](images/ROCAnalysis-AUC.png)
The diagonal dotted line represents the behavior of a random classifier. The full diagonal line represents iso-performance. A black "*A*" symbol at the bottom of the graph proportionally readjusts the graph.
5. The final box is dedicated to the analysis of the curve. The user can specify the cost of false positives (FP) and false negatives (FN), and the prior target class probability.

   - *Default threshold (0.5) point* shows the point on the ROC curve achieved by the classifier if it predicts the target class if its probability equals or exceeds 0.5.
   - *Show performance line* shows iso-performance in the ROC space so that all the points on the line give the same profit/loss. The line further to the upper left is better than the one down and right. The direction of the line depends upon costs and probabilities. This gives a recipe for depicting the optimal threshold for the given costs: this is the point where the tangent with the given inclination touches the curve and it is marked in the plot. If we push the iso-performance higher or more to the left, the points on the iso-performance line cannot be reached by the learner. Going down or to the right, decreases the performance.
   - The widget allows setting the costs from 1 to 1000. Units are not important, as are not the magnitudes. What matters is the relation between the two costs, so setting them to 100 and 200 will give the same result as 400 and 800.
   Defaults: both costs equal (500), Prior target class probability 50%(from the data).
    ![](images/ROCAnalysis-Plain.png)
    False positive cost: 830, False negative cost 650, Prior target
    class probability 73%.
    ![](images/ROCAnalysis.png)
6. Press *Save Image* if you want to save the created image to your
    computer in a .svg or .png format.
7. Produce a report.

Example
-------

At the moment, the only widget which gives the right type of signal needed by the **ROC Analysis** is [Test & Score](../evaluate/testandscore.md). Below, we compare two classifiers, namely [Tree](../model/tree.md) and [Naive Bayes](../model/naivebayes.md), in **Test\&Score** and then compare their performance in **ROC Analysis**, [Life Curve](../evaluate/liftcurve.md) and [Calibration Plot](../evaluate/calibrationplot.md).

![](images/ROCAnalysis-example.png)
CN2 Rule Viewer
===============

CN2 Rule Viewer

**Inputs**

- Data: dataset to filter
- CN2 Rule Classifier: CN2 Rule Classifier, including a list of induced rules

**Outputs**

- Filtered Data: data instances covered by all selected rules

A widget that displays [CN2 classification](https://en.wikipedia.org/wiki/CN2_algorithm) rules. If data is also connected, upon rule selection, one can analyze which instances abide to the conditions.

![](images/CN2RuleViewer-stamped.png)

1. Original order of induced rules can be restored.
2. When rules are many and complex, the view can appear packed. For this reason, *compact view* was implemented, which allows a flat presentation and a cleaner inspection of rules.
3. Click *Report* to bring up a detailed description of the rule induction algorithm and its parameters, the data domain, and induced rules.

Additionally, upon selection, rules can be copied to clipboard by pressing the default system shortcut (ctrl+C, cmd+C).

Examples
--------

In the schema below, the most common use of the widget is presented. First, the data is read and a CN2 rule classifier is trained. We are using *titanic* dataset for the rule construction. The rules are then viewed using the [Rule Viewer](../visualize/cn2ruleviewer.md). To explore different CN2 algorithms and understand how adjusting parameters influences the learning process, **Rule Viewer** should be kept open and in sight, while setting the CN2 learning algorithm (the presentation will be updated promptly).

![](images/CN2-Viewer-Example1.png)

Selecting a rule outputs filtered data instances. These can be viewed in a [Data Table](../data/datatable.md).
Radviz
======

Radviz vizualization with explorative data analysis and intelligent data
visualization enhancements.

**Inputs**

- Data: input dataset
- Data Subset: subset of instances

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected
- Components: Radviz vectors

Radviz (Hoffman et al. 1997) is a non-linear multi-dimensional visualization technique that can display data defined by three or more variables in a 2-dimensional projection. The visualized variables are presented as anchor points equally spaced around the perimeter of a unit circle. Data instances are shown as points inside the circle, with their positions determined by a metaphor from physics: each point is held in place with springs that are attached at the other end to the variable anchors. The stiffness of each spring is proportional to the value of the corresponding variable and the point ends up at the position where the spring forces are in equilibrium. Prior to visualization, variable values are scaled to lie between 0 and 1. Data instances that are close to a set of variable anchors have higher values for these variables than for the others.

The snapshot shown below shows a Radviz widget with a visualization of the dataset from functional genomics (Brown et al. 2000). In this particular visualization the data instances are colored according to the corresponding class, and the visualization space is colored according to the computed class probability. Notice that the particular visualization very nicely separates data instances of different class, making the visualization interesting and potentially informative.

![](images/Radviz-Brown.png)

Just like all point-based visualizations, this widget includes tools for intelligent data visualization (VizRank, see Leban et al. 2006) and an interface for explorative data analysis - selection of data points in visualization. Just like the [Scatter Plot](../visualize/scatterplot.md) widget, it can be used to find a set of variables that would result in an interesting visualization. The Radviz graph above is according to this definition an example of a very good visualization, while the one below - where we show an VizRank's interface (*Suggest features* button)
with a list of 3-attribute visualizations and their scores - is not.

![](images/Radviz-Brown-2.png)

References
----------

Hoffman, P. E. et al. (1997) DNA visual and analytic data mining. In the Proceedings of the IEEE Visualization. Phoenix, AZ, pp. 437-441.

Brown, M. P., W. N. Grundy et al. (2000). "Knowledge-based analysis of microarray gene expression data by using support vector machines." Proc Natl Acad Sci U S A 97(1): 262-7.

Leban, G., B. Zupan et al. (2006). "VizRank: Data Visualization Guided by Machine Learning." Data Mining and Knowledge Discovery 13(2): 119-136.

Mramor, M., G. Leban, J. Demsar, and B. Zupan. Visualization-based cancer microarray data classification analysis. Bioinformatics 23(16): 2147-2154, 2007.
Mosaic Display
==============

Display data in a mosaic plot.

**Inputs**

- Data: input dataset
- Data subset: subset of instances

**Outputs**

- Selected data: instances selected from the plot

The **Mosaic plot** is a graphical representation of a two-way frequency table or a contingency table. It is used for visualizing data from two or more qualitative variables and was introduced in 1981 by Hartigan and Kleiner and expanded and refined by Friendly in 1994. It provides the user with the means to more efficiently recognize relationships between different variables. If you wish to read up on the history of Mosaic Display, additional reading is available [here](http://www.datavis.ca/papers/moshist.pdf).

![](images/Mosaic-Display-stamped.png)

1. Select the variables you wish to see plotted.
2. Select interior coloring. You can color the interior according to class or you can use the *Pearson residual*, which is the difference between observed and fitted values, divided by an estimate of the standard deviation of the observed value. If *Compare to total* is clicked, a comparison is made to all instances.
3. *Save image* saves the created image to your computer in a .svg or .png format.
4. Produce a report.

Example
-------

We loaded the *titanic* dataset and connected it to the **Mosaic Display** widget. We decided to focus on two variables, namely status, sex and survival. We colored the interiors according to Pearson residuals in order to demonstrate the difference between observed and fitted values.

![](images/Mosaic-Display-Example.png)

We can see that the survival rates for men and women clearly deviate from the fitted value.
Silhouette Plot
===============

A graphical representation of consistency within clusters of data.

**Inputs**

- Data: input dataset

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected

The **Silhouette Plot** widget offers a graphical representation of consistency within clusters of data and provides the user with the means to visually assess cluster quality. The silhouette score is a measure of how similar an object is to its own cluster in comparison to other clusters and is crucial in the creation of a silhouette plot. The silhouette score close to 1 indicates that the data instance is close to the center of the cluster and instances possessing the silhouette scores close to 0 are on the border between two clusters.

![](images/SilhouettePlot-stamped.png)

1. Choose the distance metric. You can choose between:
   - [Euclidean](https://en.wikipedia.org/wiki/Euclidean_distance) ("straight line" distance between two points)
   - [Manhattan](https://en.wiktionary.org/wiki/Manhattan_distance) (the sum of absolute differences for all attributes)
   - [Cosine](https://en.wiktionary.org/wiki/Cosine_similarity) (1 - cosine of the angle between two vectors)
2. Select the cluster label. You can decide whether to group the instances by cluster or not.
3. Display options:
   - *Choose bar width*.
   - *Annotations*: annotate the silhouette plot.
4. *Save Image* saves the created silhouette plot to your computer in a *.png* or *.svg* format.
5. Produce a report.
6. Output:
   - *Add silhouette scores* (good clusters have higher silhouette scores)
   - By clicking *Commit*, changes are communicated to the output of the widget. Alternatively, tick the box on the left and changes will be communicated automatically.
7. The created silhouette plot.

Example
-------

In the snapshot below, we have decided to use the **Silhouette Plot** on the *iris* dataset. We selected data instances with low silhouette scores and passed them on as a subset to the [Scatter Plot](../visualize/scatterplot.md) widget. This visualization only confirms the accuracy of the **Silhouette Plot** widget, as you can clearly see that the subset lies in the border between two clusters.

![](images/SilhouettePlot-Example.png)

If you are interested in other uses of the **Silhouette Plot** widget, feel free to explore our [blog post](http://blog.biolab.si/2016/03/23/all-i-see-is-silhouette/).
Violin Plot
===========

Visualize the distribution of feature values in a violin plot.

**Inputs**

- Data: input dataset

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected

The **Violin Plot** widget plays a similar role as a [Box Plot](boxplot.md). It shows the distribution of quantitative data across several levels of a categorical variable such that those distributions can be compared. Unlike the Box Plot, in which all of the plot components correspond to actual data points, the Violin Plot features a kernel density estimation of the underlying distribution.

![](images/ViolinPlot-stamped.png)


1. Select the variable you want to plot. Tick *Order by relevance to subgroups* to order variables by Chi2 or ANOVA over the selected subgroup.
2. Choose *Subgroups* to see [violin plots](https://en.wikipedia.org/wiki/Violin_plot) displayed by a discrete subgroup. Tick *Order by relevance to variable* to order subgroups by Chi2 or ANOVA over the selected variable.
3.  *Box plot*: Tick to show the underlying box plot.
    ![](images/ViolinPlot-boxplot.png)

    *Strip plot*: Tick to show the underlying data represented by points.
    
    *Rug plot*: Tick to show the underlying data represented by lines.
    
    *Order subgroups*: Tick to order violins by *median* (ascending).
    
    *Orientation*: Determine violin orientation.
4.  *Kernel*: Select the kernel used to estimate the density. Possible kernels are: *Normal*, *Epanechnikov* and *Linear*.

    *Scale*: Select the method used to scale the width of each violin. If *area* is selected, each violin will have the same area. If *count* is selected, the width of the violins will be scaled by the number of observations in that bin. If *width* is selected, each violin will have the same width.

Examples
--------

The **Violin Plot** widget is most commonly used immediately after the [File](../data/file.md) widget to observe the statistical properties of a dataset. In the first example, we have used *heart-disease* data to inspect our variables.

![](images/ViolinPlot-example1.png)

The **Violin Plot** could also be used for *outlier detection*. In the next example we eliminate the outliers by selecting only instances that fall inside the [Q1 − 1.5  and Q3 + 1.5 IQR](https://en.wikipedia.org/wiki/Interquartile_range).

![](images/ViolinPlot-example2.png)
Sieve Diagram
=============

Plots a sieve diagram for a pair of attributes.

**Inputs**

- Data: input dataset

A **Sieve Diagram** is a graphical method for visualizing frequencies in a two-way contingency table and comparing them to [expected frequencies](http://cnx.org/contents/d396c4ad-2fd7-47cd-be84-152b44880feb@2/What-is-an-expected-frequency) under assumption of independence. It was proposed by Riedwyl and Schüpbach in a technical report in 1983 and later called a parquet diagram (Riedwyl and Schüpbach 1994). In this display, the area of each rectangle is proportional to the expected frequency, while the observed frequency is shown by the number of squares in each rectangle. The difference between observed and expected frequency (proportional to the standard Pearson residual) appears as the density of shading, using color to indicate whether the deviation from independence is positive (blue) or negative (red).

![](images/SieveDiagram-stamped.png)

1. Select the attributes you want to display in the sieve plot.
2. Score combinations enables you to fin the best possible combination of attributes.
3. *Save Image* saves the created image to your computer in a .svg or .png format.
4. Produce a report.

The snapshot below shows a sieve diagram for the *Titanic* dataset and has the attributes *sex* and *survived* (the latter is a class attribute in this dataset). The plot shows that the two variables are highly associated, as there are substantial differences between observed and expected frequencies in all of the four quadrants. For example, and as highlighted in the balloon, the chance for surviving the accident was much higher for female passengers than expected (0.06 vs. 0.15).

![](images/SieveDiagram-Titanic.png)

Pairs of attributes with interesting associations have a strong shading, such as the diagram shown in the above snapshot. For contrast, a sieve diagram of the least interesting pair (age vs. survival) is shown below.

![](images/SieveDiagram-Titanic-age-survived.png)

Example
-------

Below, we see a simple schema using the *Titanic* dataset, where we use the
[Rank](../data/rank.md) widget to select the best attributes (the ones with the highest information gain, gain ratio or Gini index) and feed them into the **Sieve Diagram**. This displays the sieve plot for the two best attributes, which in our case are sex and status. We see that the survival rate on the Titanic was very high for women of the first class and very low for female crew members.

![](images/SieveDiagram-Example2.PNG)

The **Sieve Diagram** also features the *Score Combinations* option, which makes the ranking of attributes even easier.

![](images/SieveDiagram-Example1.PNG)

References
----------

Riedwyl, H., and Schüpbach, M. (1994). Parquet diagram to plot contingency tables. In Softstat '93: Advances in Statistical Software, F. Faulbaum (Ed.). New York: Gustav Fischer, 293-299.
Tree Viewer
===========

A visualization of classification and regression trees.

**Inputs**

- Tree: decision tree

**Outputs**

- Selected Data: instances selected from the tree node
- Data: data with an additional column showing whether a point is selected

This is a versatile widget with 2-D visualization of [classification and regression trees](https://en.wikipedia.org/wiki/Decision_tree_learning). The user can select a node, instructing the widget to output the data associated with the node, thus enabling explorative data analysis.

![](images/TreeViewer-stamped.png)

1. Information on the input.
2. Display options:
   - Zoom in and zoom out
   - Select the tree width. The nodes display information bubbles when hovering over them.
   - Select the depth of your tree.
   - Select edge width. The edges between the nodes in the tree graph are drawn based on the selected edge width.
      - All the edges will be of equal width if *Fixed* is chosen.
      - When *Relative to root* is selected, the width of the edge will
         correspond to the proportion of instances in the corresponding
         node with respect to all the instances in the training data. Under
         this selection, the edge will get thinner and thinner when
         traversing toward the bottom of the tree.
      - *Relative to parent* makes the edge width correspond to the proportion
         of instances in the nodes with respect to the instances in their
         parent node.
   - Define the target class, which you can change based on classes in the data.
3. Press *Save image* to save the created tree graph to your computer as a *.svg* or *.png* file.
4. Produce a report.

Examples
--------

Below, is a simple classification schema, where we have read the data, constructed the decision tree and viewed it in our **Tree Viewer**. If both the viewer and [Tree](../model/tree.md) are open, any re-run of the tree induction algorithm will immediately affect the visualization. You can thus use this combination to explore how the parameters of the induction algorithm influence the structure of the resulting tree.

![](images/TreeViewer-classification.png)

Clicking on any node will output the related data instances. This is explored in the schema below that shows the subset in the data table and in the [Scatter Plot](../visualize/scatterplot.md). Make sure that the tree data is passed as a data subset; this can be done by connecting the **Scatter Plot** to the [File](../data/file.md) widget first, and connecting it to the **Tree Viewer** widget next. Selected data will be displayed as bold dots.

**Tree Viewer** can also export labeled data. Connect [Data Table](../data/datatable.md) to **Tree Viewer** and set the link between widgets to *Data* instead of *Selected Data*. This will send the entire data to **Data Table** with an additional meta column labeling selected data instances (*Yes* for selected and *No* for the remaining).

![](images/TreeViewer-selection.png)

Finally, **Tree Viewer** can be used also for visualizing regression trees. Connect [Random Forest](../model/randomforest.md) to [File](../data/file.md) widget using *housing.tab* dataset. Then connect [Pythagorean Forest](../visualize/pythagoreanforest.md) to **Random Forest**. In **Pythagorean Forest** select a regression tree you wish to further analyze and pass it to the **Tree Viewer**. The widget will display the constructed tree. For visualizing larger trees, especially for regression, [Pythagorean Tree](../visualize/pythagoreantree.md) could be a better option.

![](images/TreeViewer-regression.png)
Nomogram
========

Nomograms for visualization of Naive Bayes and Logistic Regression classifiers.

**Inputs**

- Classifier: trained classifier
- Data: input dataset

**Outputs**

- Features: selected variables, 10 by default

The **Nomogram** enables some classifier's (more precisely Naive Bayes classifier and Logistic Regression classifier) visual representation. It offers an insight into the structure of the training data and effects of the attributes on the class probabilities. Besides visualization of the classifier, the widget offers interactive support for prediction of class probabilities. A snapshot below shows the nomogram of the Titanic dataset, that models the probability for a passenger not to survive the disaster of the Titanic.

When there are too many attributes in the plotted dataset, only best ranked ones can be selected for display. It is possible to choose from 'No sorting', 'Name', 'Absolute importance', 'Positive influence' and 'Negative influence' for Naive Bayes representation and from 'No sorting', 'Name' and 'Absolute importance' for Logistic Regression representation.

The probability for the chosen target class is computed by '1-vs-all' principle, which should be taken in consideration when dealing with multiclass data (alternating probabilities do not sum to 1). To avoid this inconvenience, you can choose to normalize probabilities.

![](images/Nomogram-NaiveBayes.png)

1. Select the target class you want to model the probability for. Select, whether you want to normalize the probabilities or not.
2. By default Scale is set to Log odds ration. For easier understanding and interpretation option *Point scale* can be used. The unit is obtained by re-scaling the log odds so that the maximal absolute log odds ratio in the nomogram represents 100 points.
3. Display all attributes or only the best ranked ones. Sort them and set the projection type.

Continuous attributes can be plotted in 2D (only for Logistic Regression).

![logreg](images/Nomogram-LogisticRegression.png)

Examples
--------

The **Nomogram** widget should be used immediately after trained classifier widget (e.g. [Naive Bayes](../model/naivebayes.md) or [Logistics Regression](../model/logisticregression.md)). It can also be passed a data instance using any widget that enables selection (e.g. [Data Table](../data/datatable.md)) as shown in the workflow below.

![](images/Nomogram-Example.png)

Referring to the Titanic dataset once again, 1490 (68%) passengers on Titanic out of 2201 died. To make a prediction, the contribution of each attribute is measured as a point score and the individual point scores are summed to determine the probability. When the value of the attribute is unknown, its contribution is 0 points. Therefore, not knowing anything about the passenger, the total point score is 0 and the corresponding probability equals the unconditional prior. The nomogram in the example shows the case when we know that the passenger is a male adult from the first class. The points sum to -0.36, with a corresponding probability of not surviving of about 53%.

#### Features output

The second example shows how to use the Features output. Let us use *heart_disease* data for this exercise and load it in the File widget. Now connect File to [Naive Bayes](../model/naivebayes.md) (or [Logistic Regression](../model/logisticregression.md)) and add Nomogram to Naive Bayes. Finally, connect File to [Select Columns](../data/selectcolumns.md).

Select Columns selects a subset of variables, while Nomogram shows the top scoring variables for the trained classifier. To filter the data by the variables selected in the Nomogram, connect Nomogram to Select Columns as shown below. Nomogram will pass a list of selected variables to Select Columns, which will retain only the variables from the list. For this to work, you have to press *Use input features* in Select Columns (or tick it to always apply it).

We have selected the top 5 variables in Nomogram and used Select Columns to retain only those variables.

![](images/Nomogram-Features.png)
Scatter Plot
============

Scatter plot visualization with exploratory analysis and intelligent data visualization enhancements.

**Inputs**

- Data: input dataset
- Data Subset: subset of instances
- Features: list of attributes

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected

The **Scatter Plot** widget provides a 2-dimensional scatter plot visualization. The data is displayed as a collection of points, each having the value of the x-axis attribute determining the position on the horizontal axis and the value of the y-axis attribute determining the position on the vertical axis. Various properties of the graph, like color, size and shape of the points, axis titles, maximum point size and jittering can be adjusted on the left side of the widget. A snapshot below shows the scatter plot of the *Iris* dataset with the coloring matching of the class attribute.

![](images/Scatterplot-Iris-stamped.png)

1. Select the x and y attribute. Optimize your projection with **Find Informative Projections**. This feature scores attribute pairs by average classification accuracy and returns the top scoring pair with a simultaneous visualization update.
2. *Attributes*: Set the color of the displayed points (you will get colors for categorical values and blue-green-yellow points for numeric). Set label, shape and size to differentiate between points. *Label only selected points* allows you to select individual data instances and label only those.
3. Set symbol size and opacity for all data points. Set [jittering](https://en.wikipedia.org/wiki/Jitter) to prevent the dots overlapping. Jittering will randomly scatter point only around categorical values. If *Jitter numeric values* is checked, points are also scattered around their actual numeric values.
   - *Show color regions* colors the graph by class (see the screenshot below).
   - *Show legend* displays a legend on the right. Click and drag the legend to move it.
   - *Show gridlines* displays the grid behind the plot.
   - *Show all data on mouse hover* enables information bubbles if the cursor is placed on a dot.
   - *Show regression line* draws the regression line for pair of numeric attributes. If a categorical variable is selected for coloring the plot, individual regression lines for each class value will be displayed. The reported r value corresponds to the `rvalue` from [linear least-squares regression](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html), which is equal to the Pearson's correlation coefficient.
   - *Treat variables as independent* fits regression line to a group of points (minimize distance from points), rather than fitting y as a function of x (minimize vertical distances).
4. *Select, zoom, pan and zoom to fit* are the options for exploring the graph. The manual selection of data instances works as an angular/square selection tool. Double click to move the projection. Scroll in or out for zoom.
5. If *Send automatically* is ticked, changes are communicated automatically. Alternatively, press *Send*.

Here is an example of the **Scatter Plot** widget if the *Show color regions* and *Show regression line* boxes are ticked.

![](images/Scatterplot-ClassDensity.png)

Intelligent Data Visualization
------------------------------

If a dataset has many attributes, it is impossible to manually scan through all the pairs to find interesting or useful scatter plots. Orange implements intelligent data visualization with the **Find Informative Projections** option in the widget.

If a categorical variable is selected in the Color section, the [score](http://eprints.fri.uni-lj.si/210/) is computed as follows. For each data instance, the method finds 10 nearest neighbors in the projected 2D space, that is, on the combination of attribute pairs. It then checks how many of them have the same color. The total score of the projection is then the average number of same-colored neighbors.

Computation for numeric colors is similar, except that the [coefficient of determination](https://en.wikipedia.org/wiki/Coefficient_of_determination) is used for measuring the local homogeneity of the projection.

To use this method, go to the *Find Informative Projections* option in the widget, open the subwindow and press *Start Evaluation*. The feature will return a list of attribute pairs by average classification accuracy score.

Below, there is an example demonstrating the utility of ranking. The first scatter plot projection was set as the default sepal width to sepal length plot (we used the Iris dataset for simplicity). Upon running *Find Informative Projections* optimization, the scatter plot converted to a much better projection of petal width to petal length plot.

![](images/ScatterPlotExample-Ranking.png)

Selection
---------

Selection can be used to manually defined subgroups in the data. Use Shift modifier when selecting data instances to put them into a new group. Shift + Ctrl (or Shift + Cmd on macOs) appends instances to the last group.

Signal data outputs a data table with an additional column that contains group indices.

![](images/ScatterPlot-selection.png)

Exploratory Data Analysis
-------------------------

The **Scatter Plot**, as the rest of Orange widgets, supports zooming-in and out of part of the plot and a manual selection of data instances. These functions are available in the lower left corner of the widget.

The default tool is *Select*, which selects data instances within the chosen rectangular area. *Pan* enables you to move the scatter plot around the pane. With *Zoom* you can zoom in and out of the pane with a mouse scroll, while *Reset zoom* resets the visualization to its optimal size. An example of a simple schema, where we selected data instances from a rectangular region and sent them to the [Data Table](../data/datatable.md) widget, is shown below. Notice that the scatter plot doesn't show all 52 data instances, because some data instances overlap (they have the same values for both attributes used).

![](images/ScatterPlotExample-Explorative.png)

Example
-------

The **Scatter Plot** can be combined with any widget that outputs a list of selected data instances. In the example below, we combine [Tree](../model/tree.md) and **Scatter Plot** to display instances taken from a chosen decision tree node (clicking on any node of the tree will send a set of selected data instances to the scatter plot and mark selected instances with filled symbols).

![](images/ScatterPlotExample-Classification.png)

References
----------

Gregor Leban and Blaz Zupan and Gaj Vidmar and Ivan Bratko (2006) VizRank: Data Visualization Guided by Machine Learning. Data Mining and Knowledge Discovery, 13 (2). pp. 119-136. Available [here](http://eprints.fri.uni-lj.si/210/).
FreeViz
=======

Displays FreeViz projection.

**Inputs**

- Data: input dataset
- Data Subset: subset of instances

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected
- Components: FreeViz vectors

**FreeViz** uses a paradigm borrowed from particle physics: points in the same class attract each other, those from different class repel each other, and the resulting forces are exerted on the anchors of the attributes, that is, on unit vectors of each of the dimensional axis. The points cannot move (are projected in the projection space), but the attribute anchors can, so the optimization process is a hill-climbing optimization where at the end the anchors are placed such that forces are in equilibrium. The button Optimize is used to invoke the optimization process. The result of the optimization may depend on the initial placement of the anchors, which can be set in a circle, arbitrary or even manually. The later also works at any stage of optimization, and we recommend to play with this option in order to understand how a change of one anchor affects the positions of the data points. In any linear projection, projections of unit vector that are very short compared to the others indicate that their associated attribute is not very informative for particular classification task. Those vectors, that is, their corresponding anchors, may be hidden from the visualization using Radius slider in Show anchors box.

![](images/freeviz-zoo-stamped.png)

1. Two initial positions of anchors are possible: random and circular. Optimization moves anchors in an optimal position.
2. Set the color of the displayed points (you will get colors for discrete values and grey-scale points for continuous). Set label, shape and size to differentiate between points. Set symbol size and opacity for all data points.
3. Anchors inside a circle are hidden. Circle radius can be be changed using a slider.
4. Adjust plot properties:
   - Set [jittering](https://en.wikipedia.org/wiki/Jitter) to prevent the dots from overlapping (especially for discrete attributes).
   - *Show legend* displays a legend on the right. Click and drag the legend to move it.
   - *Show class density* colors the graph by class (see the screenshot below).
   - *Label only selected points* allows you to select individual data instances and label them.
5. *Select, zoom, pan and zoom to fit* are the options for exploring the graph. The manual selection of data instances works as an angular/square selection tool. Double click to move the projection. Scroll in or out for zoom.
6. If *Send automatically* is ticked, changes are communicated automatically. Alternatively, press *Send*.
7. *Save Image* saves the created image to your computer in a .svg or .png format.
8. Produce a report.

Manually move anchors
---------------------

![](images/freeviz-moveanchor.png)

One can manually move anchors. Use a mouse pointer and hover above the end of an anchor. Click the left button and then you can move selected anchor where ever you want.

Selection
---------

Selection can be used to manually defined subgroups in the data. Use Shift modifier when selecting data instances to put them into a new group. Shift + Ctrl (or Shift + Cmd on macOs) appends instances to the last group.

Signal data outputs a data table with an additional column that contains group indices.

![](images/FreeViz-selection.png)

Explorative Data Analysis
-------------------------

The **FreeViz**, as the rest of Orange widgets, supports zooming-in and out of part of the plot and a manual selection of data instances. These functions are available in the lower left corner of the widget. The default tool is *Select*, which selects data instances within the chosen rectangular area. *Pan* enables you to move the plot around the pane. With *Zoom* you can zoom in and out of the pane with a mouse scroll, while *Reset zoom* resets the visualization to its optimal size. An example of a simple schema, where we selected data instances from a rectangular region and sent them to the [Data Table](../data/datatable.md) widget, is shown below.

![](images/FreeViz-Example-Explorative.png)
Box Plot
========

Shows distribution of attribute values.

**Inputs**

- Data: input dataset

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected

The **Box Plot** widget shows the distributions of attribute values. It is a good practice to check any new data with this widget to quickly discover any anomalies, such as duplicated values (e.g., gray and grey), outliers, and alike. Bars can be selected - for example, values for categorical data or the quantile range for numeric data.

![](images/BoxPlot-Continuous.png)

1. Select the variable you want to plot. Tick *Order by relevance to subgroups* to order variables by Chi2 or ANOVA over the selected subgroup.
2. Choose *Subgroups* to see [box plots](https://en.wikipedia.org/wiki/Box_plot) displayed by a discrete subgroup. Tick *Order by relevance to variable* to order subgroups by Chi2 or ANOVA over the selected variable.
3. When instances are grouped by a subgroup, you can change the display mode. Annotated boxes will display the end values, the mean and the median, while comparing medians and compare means will, naturally, compare the selected value between subgroups.
![continuous](images/BoxPlot-Continuous-small.png)
4. The mean (the dark blue vertical line). The thin blue line represents the [standard deviation](http://mathworld.wolfram.com/StandardDeviation.html).
5. Values of the first (25%) and the third (75%) quantile. The blue highlighted area represents the values between the first and the third quartile.
6. The median (yellow vertical line).

For discrete attributes, the bars represent the number of instances with each particular attribute value. The plot shows the number of different animal types in the *Zoo* dataset: there are 41 mammals, 13 fish, 20 birds, and so on.

Display shows:
- *Stretch bars*: Shows relative values (proportions) of data instances. The unticked box shows absolute values.
- *Show box labels*: Display discrete values above each bar.
- *Sort by subgroup frequencies*: Sort subgroups by their descending frequency.

![](images/BoxPlot-Discrete.png)

Examples
--------

The **Box Plot** widget is most commonly used immediately after the [File](../data/file.md) widget to observe the statistical properties of a dataset. In the first example, we have used *heart-disease* data to inspect our variables.

![](images/BoxPlot-Example1.png)

**Box Plot** is also useful for finding the properties of a specific dataset, for instance, a set of instances manually defined in another widget (e.g. [Scatter Plot](../visualize/scatterplot.md) or instances belonging to some cluster or a classification tree node. Let us now use *zoo* data and create a typical clustering workflow with [Distances](../unsupervised/distances.md) and [Hierarchical Clustering](../unsupervised/hierarchicalclustering.md).

Now define the threshold for cluster selection (click on the ruler at the top). Connect **Box Plot** to **Hierarchical Clustering**, tick *Order by relevance*, and select *Cluster* as a subgroup. This will order attributes by how well they define the selected subgroup, in our case, a cluster. It seems like our clusters indeed correspond very well with the animal type!

![](images/BoxPlot-Example2.png)
Line Plot
=========

Visualization of data profiles (e.g., time series).

**Inputs**

- Data: input dataset
- Data Subset: subset of instances

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected

[Line plot](https://en.wikipedia.org/wiki/Line_chart) a type of plot which displays the data as a series of points, connected by straight line segments. It only works for numerical data, while categorical can be used for grouping of the data points.

![](images/LinePlot-stamped.png)

1. Information on the input data.
2. Select what you wish to display:
   - Lines show individual data instances in a plot.
   - Range shows the range of data points between 10th and 90th percentile.
   - Mean adds the line for mean value. If group by is selected, means will be displayed per each group value.
   - Error bars show the standard deviation of each attribute.
3. Select a categorical attribute to use for grouping of data instances. Use None to show ungrouped data.
4. *Select, zoom, pan and zoom to fit* are the options for exploring the graph. The manual selection of data instances works as a line selection, meaning the data under the selected line plots will be sent on the output. Scroll in or out for zoom. When hovering over an individual axis, scrolling will zoom only by the hovered-on axis (vertical or horizontal zoom).
5. If *Send Automatically* is ticked, changes are communicated automatically. Alternatively, click *Send*.

Example
-------

**Line Plot** is a standard visualization widget, which displays data profiles, normally of ordered numerical data. In this simple example, we will display the *iris* data in a line plot, grouped by the iris attribute. The plot shows how petal length nicely separates between class values.

If we observe this in a [Scatter Plot](../visualize/scatterplot.md), we can confirm this is indeed so. Petal length is an interesting attribute for separation of classes, especially when enhanced with petal width, which is also nicely separated in the line plot.

![](images/LinePlot-Example.png)
Distributions
=============

Displays value distributions for a single attribute.

**Inputs**

- Data: input dataset

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether an instance is selected
- Histogram Data: bins and instance counts from the histogram

The **Distributions** widget displays the [value distribution](https://en.wikipedia.org/wiki/Frequency_distribution) of discrete or continuous attributes. If the data contains a class variable, distributions may be conditioned on the class.

The graph shows how many times (e.g., in how many instances) each attribute value appears in the data. If the data contains a class variable, class distributions for each of the attribute values will be displayed (like in the snapshot below). To create this graph, we used the *Zoo* dataset.

![](images/Distributions-Discrete.png)

1. A list of variables for display. *Sort categories by frequency* orders displayed values by frequency.
2. Set *Bin width* with the slider. Precision scale is set to sensible intervals. *Fitted distribution* fits selected distribution to the plot. Options are [Normal](https://en.wikipedia.org/wiki/Normal_distribution), [Beta](https://en.wikipedia.org/wiki/Beta_distribution), [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution), [Rayleigh](https://en.wikipedia.org/wiki/Rayleigh_distribution), [Pareto](https://en.wikipedia.org/wiki/Pareto_distribution), [Exponential](https://en.wikipedia.org/wiki/Exponential_distribution), [Kernel density](https://en.wikipedia.org/wiki/Kernel_density_estimation).
3. Columns:

- *Split by* displays value distributions for instances of a certain class.
- *Stack columns* displays one column per bin, colored by proportions of class values.
- *Show probabilities* shows probabilities of class values at selected variable.
- *Show cumulative distribution* cumulatively stacks frequencies.

4. If *Apply Automatically* is ticked, changes are communicated automatically. Alternatively, click *Apply*.

For continuous attributes, the attribute values are also displayed as a histogram. It is possible to fit various distributions to the data, for example, a Gaussian kernel density estimation. *Hide bars* hides histogram bars and shows only distribution (old behavior of Distributions).

For this example, we used the *Iris* dataset.

![](images/Distributions-Continuous.png)

In class-less domains, the bars are displayed in blue. We used the *Housing* dataset.

![](images/Distributions-NoClass.png)Pythagorean Tree
================

Pythagorean tree visualization for classification or regression trees.

**Inputs**

- Tree: tree model
- Selected Data: instances selected from the tree

**Pythagorean Trees** are plane fractals that can be used to depict general tree hierarchies as presented in an article by [Fabian Beck and co-authors](http://publications.fbeck.com/ivapp14-pythagoras.pdf). In our case, they are used for visualizing and exploring tree models, such as [Tree](../model/tree.md).

![](images/Pythagorean-Tree1-stamped.png)

1. Information on the input tree model.
2. Visualization parameters:
    - *Depth*: set the depth of displayed trees.
    - *Target class* (for classification trees): the intensity of the color for nodes of the tree will correspond to the probability of the target class. If *None* is selected, the color of the node will denote the most probable class.
    - *Node color* (for regression trees): node colors can correspond to mean or standard deviation of class value of the training data instances in the node.
    - *Size*: define a method to compute the size of the square representing the node. *Normal* will keep node sizes correspond to the size of training data subset in the node. *Square root* and *Logarithmic* are the respective transformations of the node size.
    - *Log scale factor* is only enabled when *logarithmic* transformation is selected. You can set the log factor between 1 and 10.
3. Plot properties:
    - *Enable tooltips*: display node information upon hovering.
    - *Show legend*: shows color legend for the plot.
4. Reporting:
    - *Save Image*: save the visualization to a SVG or PNG file.
    - *Report*: add visualization to the report.

Pythagorean Tree can visualize both classification and regression trees. Below is an example for regression tree. The only difference between the two is that regression tree doesn't enable coloring by class, but can color by class mean or standard deviation.

![](images/Pythagorean-Tree1-continuous.png)

Example
-------

The workflow from the screenshot below demonstrates the difference between [Tree Viewer](../visualize/treeviewer.md) and Pythagorean Tree. They can both visualize [Tree](../model/tree.md), but Pythagorean visualization takes less space and is more compact, even for a small [Iris flower](https://en.wikipedia.org/wiki/Iris_flower_data_set) dataset. For both visualization widgets, we have hidden the control area on the left by clicking on the splitter between control and visualization area.

![](images/Pythagorean-Tree-comparison.png)

Pythagorean Tree is interactive: click on any of the nodes (squares) to select training data instances that were associated with that node. The following workflow explores these feature.

![](images/Pythagorean-Tree-scatterplot-workflow.png)

The selected data instances are shown as a subset in the [Scatter Plot](../visualize/scatterplot.md), sent to the [Data Table](../data/datatable.md) and examined in the [Box Plot](../visualize/boxplot.md). We have used brown-selected dataset in this example. The tree and scatter plot are shown below; the selected node in the tree has a black outline.

![](images/Pythagorean-Tree-scatterplot.png)

References
----------

Beck, F., Burch, M., Munz, T., Di Silvestro, L. and Weiskopf, D. (2014). [Generalized Pythagoras Trees for Visualizing Hierarchies](http://publications.fbeck.com/ivapp14-pythagoras.pdf). In IVAPP '14 Proceedings of the 5th International Conference on Information Visualization Theory and Applications, 17-28.
Bar Plot
========

Visualizes comparisons among discrete categories.

**Inputs**

- Data: input dataset
- Data Subset: subset of instances

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected

The **Bar Plot** widget visualizes numeric variables and compares them by a categorical variable. The widget is useful for observing outliers, distributions within groups, and comparing categories.

![](images/Bar-Plot-stamped.png)

1. Parameters of the plot. Values are the numeric variable to plot. Group by is the variable for grouping the data. Annotations are categorical labels below the plot. Color is the categorical variable whose values are used for coloring the bars.
2. *Select, zoom, pan and zoom to fit* are the options for exploring the graph. The manual selection of data instances works as an angular/square selection tool. Double click to move the projection. Scroll in or out for zoom.
3. If *Send automatically* is ticked, changes are communicated automatically. Alternatively, press *Send*.
4. Access help, save image, produce a report, or adjust visual settings. On the right, the information on input and output are shown.

Example
-------

The **Bar Plot** widget is most commonly used immediately after the [File](../data/file.md) widget to compare categorical values. In this example, we have used *heart-disease* data to inspect our variables.

![](images/Bar-Plot-Example.png)

First, we have observed cholesterol values of patient from our data set. We grouped them by diameter narrowing, which defines patients with a heart disease (1) and those without (0). We use the same variable for coloring the bars.

Then, we selected patients over 60 years of age with [Select Rows](../data/selectrows.md). We sent the subset to **Bar Plot** to highlight these patients in the widget. The big outlier with a high cholesterol level is apparently over 60 years old.
Pythagorean Forest
==================

Pythagorean forest for visualizing random forests.

**Inputs**

- Random Forest: tree models from random forest

**Outputs**

- Tree: selected tree model

**Pythagorean Forest** shows all learned decision tree models from [Random Forest](../model/randomforest.md) widget. It displays them as Pythagorean trees, each visualization pertaining to one randomly constructed tree. In the visualization, you can select a tree and display it in [Pythagorean Tree](../visualize/pythagoreantree.md) widget. The best tree is the one with the shortest and most strongly colored branches. This means few attributes split the branches well.

Widget displays both classification and regression results. Classification requires discrete target variable in the dataset, while regression requires a continuous target variable. Still, they both should be fed a [Tree](../model/tree.md) on the input.

![](images/Pythagorean-Forest-stamped.png)

1. Information on the input random forest model.
2. Display parameters:
    - *Depth*: set the depth to which the trees are grown.
    - *Target class*: set the target class for coloring the trees. If *None* is selected, the tree will be white. If the input is a classification tree, you can color the nodes by their respective class. If the input is a regression tree, the options are *Class mean*, which will color tree nodes by the class mean value and *Standard deviation*, which will color them by the standard deviation value of the node.
    - *Size*: set the size of the nodes. *Normal* will keep the nodes the size of the subset in the node. *Square root* and *Logarithmic* are the respective transformations of the node size.
    - *Zoom*: allows you to see the size of the tree visualizations.
3. *Save Image*: save the visualization to your computer as a *.svg* or *.png* file. *Report*: produce a report.

Example
-------

**Pythagorean Forest** is great for visualizing several built trees at once. In the example below, we've used *housing* dataset and plotted all 10 trees we've grown with [Random Forest](../model/randomforest.md). When changing the parameters in Random Forest, visualization in Pythagorean Forest will change as well.

Then we've selected a tree in the visualization and inspected it further with [Pythagorean Tree](../visualize/pythagoreantree.md) widget.

![](images/Pythagorean-Forest-Example.png)

References
----------

Beck, F., Burch, M., Munz, T., Di Silvestro, L. and Weiskopf, D. (2014). Generalized Pythagoras Trees for Visualizing Hierarchies. In IVAPP '14 Proceedings of the 5th International Conference on Information Visualization Theory and Applications, 17-28.
Venn Diagram
============

Plots a [Venn diagram](http://en.wikipedia.org/wiki/Venn_diagram) for two or more data subsets.

**Inputs**

- Data: input dataset

**Outputs**

- Selected Data: instances selected from the plot
- Data: entire data with a column indicating whether an instance was selected or not

The **Venn Diagram** widget displays logical relations between datasets by showing the number of common data instances (rows) or the number of shared features (columns). Selecting a part of the visualization outputs the corresponding instances or features.

![](images/venn-workflow.png)

![](images/VennDiagram-stamped.png)

1. Select whether to count common features or instances.
2. Select whether to include duplicates or to output only unique rows; applicable only when matching instances by values of variables.

Rows can be matched
- by their identity, e.g. rows from different data sets match if they came from the same row in a file,
- by equality, if all tables contain the same variables,
- or by values of a string variable that appears in all tables.

Examples
--------

The easiest way to use the **Venn Diagram** is to select data subsets and find matching instances in the visualization. We use the *breast-cancer* dataset to select two subsets with [Select Rows](../data/selectrows.md) widget - the first subset is that of breast cancer patients aged between 40 and 49 and the second is that of patients with a tumor size between 20 and 29. The **Venn Diagram** helps us find instances that correspond to both criteria, which can be found in the intersection of the two circles.

![](images/VennDiagram-Example1.png)

The **Venn Diagram** widget can be also used for exploring different prediction models. In the following example, we analysed 3 prediction methods, namely [Naive Bayes](../model/naivebayes.md), [SVM](../model/svm.md) and [Random Forest](../model/randomforest.md), according to their misclassified instances.

By selecting misclassifications in the three [Confusion Matrix](../evaluate/confusionmatrix.md) widgets and sending them to Venn diagram, we can see all the misclassification instances visualized per method used. Then we open **Venn Diagram** and select, for example, the misclassified instances that were identified by all three methods. This is represented as an intersection of all three circles. Click on the intersection to see this two instances marked in the [Scatter Plot](../visualize/scatterplot.md) widget. Try selecting different diagram sections to see how the scatter plot visualization changes.

![](images/VennDiagram-Example2.png)
Linear Projection
=================

A linear projection method with explorative data analysis.

**Inputs**

- Data: input dataset
- Data Subset: subset of instances
- Projection: custom projection vectors

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected
- Components: projection vectors

This widget displays [linear projections](https://en.wikipedia.org/wiki/Projection_(linear_algebra)) of class-labeled data. It supports various types of projections such as circular, [linear discriminant analysis](https://en.wikipedia.org/wiki/Linear_discriminant_analysis), and [principal component analysis](https://en.wikipedia.org/wiki/Principal_component_analysis).

Consider, for a start, a projection of the *Iris* dataset shown below. Notice that it is the sepal width and sepal length that already separate *Iris setosa* from the other two, while the petal length is the attribute best separating *Iris versicolor* from *Iris virginica*.

![](images/LinearProjection-stamped.png)

1. Axes in the projection that are displayed and other available axes. Optimize your projection by using **Suggest Features**. This feature scores attributes and returns the top scoring attributes with a simultaneous visualization update. Feature scoring computes the classification accuracy (for classification) or MSE (regression) of k-nearest neighbors classifier on the projected, two-dimensional data. The score reflects how well the classes in the projection are separated.
2. Choose the type of projection:
   - Circular Placement
   - [Linear Discriminant Analysis](https://en.wikipedia.org/wiki/Linear_discriminant_analysis)
   - [Principal Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis)
3. Set the color of the displayed points. Set shape, size, and label to differentiate between points.
   *Label only selected points* labels only selected data instances.
4. Adjust plot properties:
   - *Symbol size*: set the size of the points.
   - *Opacity*: set the transparency of the points.
   - *Jittering*: Randomly disperse points with [jittering](https://en.wikipedia.org/wiki/Jitter) to prevent them from overlapping.
   - *Hide radius*: Axes inside the radius are hidden. Drag the slider to change the radius.
5. Additional plot properties:
   - *Show color regions* colors the graph by class.
   - *Show legend* displays a legend on the right. Click and drag the legend to move it.
6. *Select, zoom, pan* and *zoom to fit* are the options for exploring the graph. Manual selection of data instances works as an angular/square selection tool. Double click to move the projection. Scroll in or out for zoom.
7. If *Send automatically* is ticked, changes are communicated automatically. Alternatively, press *Send*.

Example
-------

The **Linear Projection** widget works just like other visualization widgets. Below, we connected it to the [File](../data/file.md) widget to see the set projected on a 2-D plane. Then we selected the data for further analysis and connected it to the [Data Table](../data/datatable.md) widget to see the details of the selected subset.

![](images/LinearProjection-Example.png)

References
----------

Koren Y., Carmel L. (2003). Visualization of labeled data using linear transformations. In Proceedings of IEEE Information Visualization 2003, (InfoVis'03). Available [here](http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=3DDF0DB68D8AB9949820A19B0344C1F3?doi=10.1.1.13.8657&rep=rep1&type=pdf).

Boulesteix A.-L., Strimmer K. (2006). Partial least squares: a versatile tool for the analysis of high-dimensional genomic data. Briefings in Bioinformatics, 8(1), 32-44. Abstract [here](http://bib.oxfordjournals.org/content/8/1/32.abstract).

Leban G., Zupan B., Vidmar G., Bratko I. (2006). VizRank: Data Visualization Guided by Machine Learning. Data Mining and Knowledge Discovery, 13, 119-136. Available [here](http://eprints.fri.uni-lj.si/210/2/1._G._Leban%2C_B._Zupan%2C_G._Vidmar%2C_I._Bratko%2C_Data_Mining_and_Knowledge_Discovery_13%2C_119-36_(2006)..pdf).
Heat Map
========

Plots a heat map for a pair of attributes.

**Inputs**

- Data: input dataset

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected

[Heat map](https://en.wikipedia.org/wiki/Heat_map) is a graphical method for visualizing attribute values in a two-way matrix. It only works on datasets containing numeric variables. The values are represented by color according to the selected color pallette. By combining class variable and attributes on x and y axes, we see where the attribute values are the strongest and where the weakest, thus enabling us to find typical features for each class.

The widget enables row selection with click and drag. One can zoom in with Ctrl++ (Cmd++) and zoom out with Ctrl+- (Cmd+-). Ctrl+0 (Cmd+0) resets zoom to the extended version, while Ctrl+9 (Cmd+9) reset it to the default.

![](images/HeatMap.png)

1. The color pallette. Choose from linear, diverging, color-blind friendly, or other pallettes. **Low** and **High** are thresholds for the color palette (low for attributes with low values and high for attributes with high values). Selecting one of diverging palettes, which have two extreme colors and a neutral (black or white) color at the midpoint, enables an option to set a meaningful mid-point value (default is 0).
2. Merge rows. If there are too many rows in the visualization, one can merge them with k-means algorithm into N selected clusters (default 50).
3. Cluster columns and rows:
   - **None** (lists attributes and rows as found in the dataset)
   - **Clustering** (clusters data by similarity with hierarchical clustering on Euclidean distances and with average linkage)
   - **Clustering with ordered leaves** (same as clustering, but it additionally maximizes the sum of similarities of adjacent elements)
4. Split rows or columns by a categorical variable. If the data contains a class variable, rows will be automatically split by class.
5. Set what is displayed in the plot in **Annotation & Legend**.
   - If *Show legend* is ticked, a color chart will be displayed above the map.
   - If *Stripes with averages* is ticked, a new line with attribute averages will be displayed on the left.
   **Row Annotations** adds annotations to each instance on the right. Color colors the instances with the corresponding value of the selected categorical variable.
   **Column Annotations** adds annotation to each variable at the selected position (default is Top). Color colors the columns with the corresponding value of the selected column annotation.
6. If *Keep aspect ratio* is ticked, each value will be displayed with a square (proportionate to the map).
7. If *Send Automatically* is ticked, changes are communicated automatically. Alternatively, click *Send*.

### Advanced visualization

Heat map enables some neat plot enhancements. Such options are clustering of rows and/or columns for better data organization, row and column annotations, and splitting the data by categorical variables.

Row and column clustering is performed independently. Row clustering is computed from Euclidean distances, while column clustering uses Pearson correlation coefficients. Hierarchical clustering is based on the Ward linkage method. Clustering with optimal leaf ordering reorders left and right branches in the dendrogram to minimize the sum of distances between adjacent leaves (Bar-Joseph et al. 2001).



![](images/HeatMap-advanced.png)

Examples
--------

### Gene expresssions

The **Heat Map** below displays attribute values for the *brown-selected* data set (Brown et al. 2000). Heat maps are particularly appropriate for showing gene expressions and the brown-selected data set contains yeast gene expressions at different conditions.

Heat map shows low expressions in blue and high expressions in yellow and white. For better organization, we added *Clustering (opt. ordering)* to the columns, which puts columns with similar profiles closer together. In this way we can see the conditions that result in low expressions for ribosomal genes in the lower right corner.

Additionally, the plot is enhanced with row color on the right, showing which class the rows belong to.

![](images/HeatMap-Example1.png)

### Sentiment Analysis

Heat maps are great for visualizing any kind of comparable numeric variables, for example sentiment in a collection of documents. We will take *book-excerpts* corpus from the **Corpus** widget and pass it to the **Sentiment Analysis** widget, which computes sentiment scores for each document. The output of sentiment analysis are four columns, positive, negative, and neutral sentiment score, and a compound score that aggregates the previous scores into a single number. Positive compound values (white) represent positive documents, while negative (blue) represent negative documents.

We used row clustering to place similar rows closer together, resulting in clear negative and positive groups. Now we can select negative children's books and explore which are they.

![](images/HeatMap-Example2.png)

References
----------

Bar-Joseph, Z., Gifford, D.K., Jaakkola, T.S. (2001) Fast optimal leaf ordering for hierarchical clustering, Bioinformatics, 17, 22-29.

Brown, M.P., Grundy, W.N., Lin, D., Cristianini, N., Sugnet, C., Furey, T.S., Ares, M., Haussler, D. (2000) Knowledge-based analysis of microarray gene expression data by using support vector machines, Proceedings of the National Academy of Sciences, 1, 262-267.
Distance Map
============

Visualizes distances between items.

**Inputs**

- Distances: distance matrix

**Outputs**

- Data: instances selected from the matrix
- Features: attributes selected from the matrix

The **Distance Map** visualizes distances between objects. The visualization is the same as if we printed out a table of numbers, except that the numbers are replaced by colored spots.

Distances are most often those between instances ("*rows*" in the [Distances](../unsupervised/distances.md) widget) or attributes ("*columns*" in Distances widget). The only suitable input for **Distance Map** is the [Distances](../unsupervised/distances.md) widget. For the output, the user can select a region of the map and the widget will output the corresponding instances or attributes. Also note that the **Distances** widget ignores discrete values and calculates distances only for continuous data, thus it can only display distance map for discrete data if you [Continuize](../data/continuize.md) them first.

The snapshot shows distances between columns in the *heart disease* data, where smaller distances are represented with light and larger with dark orange. The matrix is symmetric and the diagonal is a light shade of orange - no attribute is different from itself. Symmetricity is always assumed, while the diagonal may also be non-zero.

![](images/DistanceMap-stamped.png)

1. *Element sorting* arranges elements in the map by
   - None (lists instances as found in the dataset)
   - **Clustering** (clusters data by similarity)
   - **Clustering with ordered leaves** (maximizes the sum of similarities of adjacent elements)
2. *Colors*
   - **Colors** (select the color palette for your distance map)
   - **Low** and **High** are thresholds for the color palette (low for instances or attributes with low distances and high for instances or attributes with high distances).
3. Select *Annotations*.
4. If *Send Selected Automatically* is on, the data subset is communicated automatically, otherwise you need to press *Send Selected*.
5. Press *Save Image* if you want to save the created image to your computer.
6. Produce a report.

Normally, a color palette is used to visualize the entire range of distances appearing in the matrix. This can be changed by setting the low and high threshold. In this way we ignore the differences in distances outside this interval and visualize the interesting part of the distribution.

Below, we visualized the most correlated attributes (distances by columns) in the *heart disease* dataset by setting the color threshold for high distances to the minimum. We get a predominantly black square, where attributes with the lowest distance scores are represented by a lighter shade of the selected color schema (in our case: orange). Beside the diagonal line, we see that in our example *ST by exercise* and *major vessels colored* are the two attributes closest together.

![](images/DistanceMap-Highlighted.png)

The user can select a region in the map with the usual click-and-drag of the cursor. When a part of the map is selected, the widget outputs all items from the selected cells.

Examples
--------

The first workflow shows a very standard use of the **Distance Map** widget. We select 70% of the original *Iris* data as our sample and view the distances between rows in **Distance Map**.

![](images/DistanceMap-Example1.png)

In the second example, we use the *heart disease* data again and select a subset of women only from the [Scatter Plot](../visualize/scatterplot.md). Then, we visualize distances between columns in the **Distance Map**. Since the subset also contains some discrete data, the [Distances](../unsupervised/distances.md) widget warns us it will ignore the discrete features, thus we will see only continuous instances/attributes in the map.

![](images/DistanceMap-Example.png)
PCA
===

PCA linear transformation of input data.

**Inputs**

- Data: input dataset

**Outputs**

- Transformed Data: PCA transformed data
- Components: [Eigenvectors](https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors).

[Principal Component Analysis](https://en.wikipedia.org/wiki/Principal_component_analysis) (PCA) computes the PCA linear transformation of the input data. It outputs either a transformed dataset with weights of individual instances or weights of principal components.

![](images/PCA-stamped.png)

1. Select how many principal components you wish in your output. It is best to choose as few as possible with variance covered as high as possible. You can also set how much variance you wish to cover with your principal components.
2. You can normalize data to adjust the values to common scale. If checked, columns are divided by their standard deviations.
3. When *Apply Automatically* is ticked, the widget will automatically communicate all changes. Alternatively, click *Apply*.
4. Press *Save Image* if you want to save the created image to your computer.
5. Produce a report.
6. Principal components graph, where the red (lower) line is the variance covered per component and the green (upper) line is cumulative variance covered by components.

The number of components of the transformation can be selected either in the *Components Selection* input box or by dragging the vertical cutoff line in the graph.

Preprocessing
-------------

The widget preprocesses the input data in the following order:

- continuizes categorical variables (with one-hot-encoding)
- imputes missing values with mean values
- if *Normalize variables* is checked, it divides columns by their standard deviation.

Examples
--------

**PCA** can be used to simplify visualizations of large datasets. Below, we used the *Iris* dataset to show how we can improve the visualization of the dataset with PCA. The transformed data in the [Scatter Plot](../visualize/scatterplot.md) show a much clearer distinction between classes than the default settings.

![](images/PCAExample.png)

The widget provides two outputs: transformed data and principal components. Transformed data are weights for individual instances in the new coordinate system, while components are the system descriptors (weights for principal components). When fed into the [Data Table](../data/datatable.md), we can see both outputs in numerical form. We used two data tables in order to provide a more clean visualization of the workflow, but you can also choose to edit the links in such a way that you display the data in just one data table. You only need to create two links and connect the *Transformed data* and *Components* inputs to the *Data* output.

![](images/PCAExample2.png)
Hierarchical Clustering
=======================

Groups items using a hierarchical clustering algorithm.

**Inputs**

- Distances: distance matrix

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether an instance is selected

The widget computes [hierarchical clustering](https://en.wikipedia.org/wiki/Hierarchical_clustering) of arbitrary types of objects from a matrix of distances and shows a corresponding [dendrogram](https://en.wikipedia.org/wiki/Dendrogram).

![](images/HierarchicalClustering-stamped.png)

1. The widget supports four ways of measuring distances between clusters:
   - **Single linkage** computes the distance between the closest elements of the two clusters
   - **Average linkage** computes the average distance between elements of the two clusters
   - **Weighted linkage** uses the [WPGMA](http://research.amnh.org/~siddall/methods/day1.html) method
   - **Complete linkage** computes the distance between the clusters' most distant elements
2. Labels of nodes in the dendrogram can be chosen in the **Annotation** box.
3. Huge dendrograms can be pruned in the *Pruning* box by selecting the maximum depth of the dendrogram. This only affects the display, not the actual clustering.
4. The widget offers three different selection methods:
   - **Manual** (Clicking inside the dendrogram will select a cluster. Multiple clusters can be selected by holding Ctrl/Cmd. Each selected cluster is shown in a different color and is treated as a separate cluster in the output.)
   - **Height ratio** (Clicking on the bottom or top ruler of the dendrogram places a cutoff line in the graph. Items to the right of the line are selected.)
   - **Top N** (Selects the number of top nodes.)
5. Use *Zoom* and scroll to zoom in or out.
6. If the items being clustered are instances, they can be added a cluster index (*Append cluster IDs*). The ID can appear as an ordinary **Attribute**, **Class attribute** or a **Meta attribute**. In the second case, if the data already has a class attribute, the original class is placed among meta attributes.
7. The data can be automatically output on any change (*Auto send is on*) or, if the box isn't ticked, by pushing *Send Data*.
8. Clicking this button produces an image that can be saved.
9. Produce a report.

Examples
--------

The workflow below shows the output of **Hierarchical Clustering** for the *Iris* dataset in [Data Table](../data/datatable.md) widget. We see that if we choose *Append cluster IDs* in hierarchical clustering, we can see an additional column in the **Data Table** named *Cluster*. This is a way to check how hierarchical clustering clustered individual instances.

![](images/HierarchicalClustering-Example.png)

In the second example, we loaded the *Iris* dataset again, but this time we added the [Scatter Plot](../visualize/scatterplot.md), showing all the instances from the [File](../data/file.md) widget, while at the same time receiving the selected instances signal from **Hierarchical Clustering**. This way we can observe the position of the selected cluster(s) in the projection.

![](images/HierarchicalClustering-Example2.png)
Self-Organizing Map
===================

Computation of a self-organizing map.

**Inputs**

- Data: input dataset

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected

A [self-organizing map (SOM)](https://en.wikipedia.org/wiki/Self-organizing_map) is a type of artificial neural network (ANN) that is trained using unsupervised learning to produce a two-dimensional, discretized representation of the data. It is a method to do dimensionality reduction. Self-organizing maps use a neighborhood function to preserve the topological properties of the input space.

The points in the grid represent data instances. By default, the size of the point corresponds to the number of instances represented by the point. The points are colored by majority class (if available), while the intensity of interior color shows the proportion of majority class. To see the class distribution, select *Show pie charts* option.

Just like other visualization widgets, **Self-Organizing Maps** also supports interactive selection of groups. Use Shift key to select a new group and Ctr+Shift to add to the existing group.

![](images/Self-Organizing_Map-stamped.png)

1. SOM properties:
   - Set the grid type. Options are hexagonal or square grid.
   - If *Set dimensions automatically* is checked, the size of the plot will be set automatically. Alternatively, set the size manually.
   - Set the initialization type for the SOM projection. Options are PCA initialization, random initialization and replicable random (random_seed = 0).
   - Once the parameters are set, press *Start* to re-run the optimization.
2. Set the color of the instances in the plot. The widget colors by class by default (if available).
   - *Show pie charts* turns points into pie-charts that show the distributions of the values used for coloring.
   - *Size by number of instances* scales the points according to the number of instances represented by the point.

Example
-------

Self-organizing maps are low-dimensional projections of the input data. We will use the *brown-selected* data and display the data instance in a 2-D projection. Seems like the three gene types are well-separated. We can select a subset from the grid and display it in a Data Table.

![](images/Self-Organizing_Map_Example.png)
Distance Matrix
===============

Visualizes distance measures in a distance matrix.

**Inputs**

- Distances: distance matrix

**Outputs**

- Distances: distance matrix
- Table: distance measures in a distance matrix

The **Distance Matrix** widget creates a distance matrix, which is a two-dimensional array containing the distances, taken pairwise, between the elements of a set. The number of elements in the dataset defines the size of the matrix. Data matrices are essential for hierarchical clustering and they are extremely useful in bioinformatics as well, where they are used to represent protein structures in a coordinate-independent manner.

![](images/DistanceMatrix-stamped.png)

1. Elements in the dataset and the distances between them.
2. Label the table. The options are: *none*, *enumeration*, *according to variables*.
3. Produce a report.
4. Click *Send* to communicate changes to other widgets. Alternatively, tick the box in front of the *Send* button and changes will be communicated automatically (*Send Automatically*).

The only two suitable inputs for **Distance Matrix** are the [Distances](../unsupervised/distances.md) widget and the [Distance Transformation](../unsupervised/distancetransformation.md) widget. The output of the widget is a data table containing the distance matrix. The user can decide how to label the table and the distance matrix (or instances in the distance matrix) can then be visualized or displayed in a separate data table.

Example
-------

The example below displays a very standard use of the **Distance Matrix** widget. We compute the distances between rows in the sample from the *Iris* dataset and output them in the **Distance Matrix**. It comes as no surprise that Iris Virginica and Iris Setosa are the furthest apart.

![](images/DistanceMatrix-Example.png)
k-Means
=======

Groups items using the k-Means clustering algorithm.

**Inputs**

- Data: input dataset

**Outputs**

- Data: dataset with cluster label as a meta attribute
- Centroids: table with initial centroid coordinates

The widget applies the [k-Means clustering](https://en.wikipedia.org/wiki/K-means_clustering) algorithm to the data and outputs a new dataset in which the cluster label is added as a meta attribute. Silhouette scores of clustering results for various k are also shown in the widget. When using the silhouette score option, the higher the silhouette score, the better the clustering.

![](images/kMeans-stamped.png)

1. Select the number of clusters.
   - **Fixed**: algorithm clusters data to a specified number of clusters.
   - **From X to Y**: widget shows clustering scores for the selected cluster range using the [Silhouette](https://en.wikipedia.org/wiki/Silhouette_\(clustering\)) score (contrasts average distance to elements in the same cluster with the average distance to elements in other clusters).
2. **Preprocessing**: If the option is selected, columns are normalized (mean centered to 0 and standard deviation scaled to 1).
3. Initialization method (the way the algorithm begins clustering):
   - [k-Means++](https://en.wikipedia.org/wiki/K-means%2B%2B) (first center is selected randomly, subsequent are chosen from the remaining points with probability proportioned to squared distance from the closest center)
   - **Random initialization** (clusters are assigned randomly at first and then updated with further iterations)

    **Re-runs** (how many times the algorithm is run from random initial positions; the result with the lowest within-cluster sum of squares will be used) and **Maximum iterations** (the maximum number of iterations within each algorithm run) can be set manually.

Preprocessing
-------------

k-Means uses default preprocessing if necessary. It executes it in the following order:

- continuizes categorical variables (with one feature per value)
- imputes missing values with mean values

To override default preprocessing, preprocess the data beforehand with [Preprocess](../data/preprocess.md) widget.

Examples
--------

First, we load the *Iris* dataset, run k-Means with three clusters, and show it in the [Scatter Plot](../visualize/scatterplot.md). To interactively explore the clusters, we can use [Select Rows](../data/selectrows.md) to select the cluster of interest (say, C1) and plot it in the scatter plot using interactive data analysis. That means if we pass a subset to the scatter plot, the subset will be exposed in the plot.

Try the same procedure for 2 or 4 clusters or explore different clusters in the plot (C2, C3).

![](images/kMeans-Example1.png)

But as we used silhouette score to estimate our cluster quality, we can plot the clusters in the [Silhouette Plot](../visualize/silhouetteplot.md) to observe inliers and outliers. Place Silhouette Plot in place of Select Rows.

Silhouette Plot shows silhouette scores for individual data instances. High, positive scores represent instances that are highly representative of the clusters, while negative scores represent instances that are outliers (don't fit well with the cluster). Select negative scores from the green cluster C3 and plot them in a scatter plot as a subset.

It seems like these are mostly iris versicolors, which are bordering the iris virginica region. Note that the green color of the cluster C3 doesn't coincide with the green color of the iris labels - these are two different things.

![](images/kMeans-Example2.png)
Distance Transformation
=======================

Transforms distances in a dataset.

**Inputs**

- Distances: distance matrix

**Outputs**

- Distances: transformed distance matrix

The **Distances Transformation** widget is used for the normalization and inversion of distance matrices. The normalization of data is necessary to bring all the variables into proportion with one another.

![](images/DistanceTransformation-stamped.png)

1. Choose the type of [Normalization](https://en.wikipedia.org/wiki/Normalization_\(statistics\)):
   - **No normalization**
   - **To interval [0, 1]**
   - **To interval [-1, 1]**
   - [Sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function): 1/(1+exp(-X))
2. Choose the type of Inversion:
   - **No inversion**
   - **-X**
   - **1 - X**
   - **max(X) - X**
   - **1/X**
3. Produce a report.
4. After changing the settings, you need to click *Apply* to commit changes to other widgets. Alternatively, tick *Apply automatically*.

Example
-------

In the snapshot below, you can see how transformation affects the distance matrix. We loaded the *Iris* dataset and calculated the distances between rows with the help of the [Distances](../unsupervised/distances.md) widget. In order to demonstrate how **Distance Transformation** affects the [Distance Matrix](../unsupervised/distancematrix.md), we created the workflow below and compared the transformed distance matrix with the "original" one.

![](images/DistanceTransformation-Example.png)
Distance File
=============

Loads an existing distance file.

**Outputs**

- Distance File: distance matrix

![](images/DistanceFile-stamped.png)

1. Choose from a list of previously saved distance files.
2. Browse for saved distance files.
3. Reload the selected distance file.
4. Information about the distance file (number of points,
    labelled/unlabelled).
5. Browse documentation datasets.
6. Produce a report.

Example
-------

When you want to use a custom-set distance file that you've saved before, open the **Distance File** widget and select the desired file with the *Browse* icon. This widget loads the existing distance file. In the snapshot below, we loaded the transformed *Iris* distance matrix from the [Save Distance Matrix](../unsupervised/savedistancematrix.md) example. We displayed the transformed data matrix in the [Distance Map](../unsupervised/distancemap.md) widget. We also decided to display a distance map of the original *Iris* dataset for comparison.

![](images/DistanceFile-Example.png)
t-SNE
=====

Two-dimensional data projection with t-SNE.

**Inputs**

- Data: input dataset
- Data Subset: subset of instances

**Outputs**

- Selected Data: instances selected from the plot
- Data: data with an additional column showing whether a point is selected

The **t-SNE** widget plots the data with a t-distributed stochastic neighbor embedding method. [t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding) is a dimensionality reduction technique, similar to MDS, where points are mapped to 2-D space by their probability distribution.

![](images/tSNE-stamped.png)

1. [Parameters](https://opentsne.readthedocs.io/en/latest/parameters.html) for plot optimization:
   - measure of [perplexity](http://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html). Roughly speaking, it can be interpreted as the number of nearest neighbors to distances will be preserved from each point. Using smaller values can reveal small, local clusters, while using large values tends to reveal the broader, global relationships between data points.
   - *Preserve global structure*: this option will combine two different perplexity values (50 and 500) to try preserve both the local and global structure.
   - *Exaggeration*: this parameter increases the attractive forces between points, and can directly be used to control the compactness of clusters. Increasing exaggeration may also better highlight the global structure of the data. t-SNE with exaggeration set to 4 is roughly equal to UMAP.
   - *PCA components*: in Orange, we always run t-SNE on the principal components of the input data. This parameter controls the number of principal components to use when calculating distances between data points.
   - *Normalize data*: We can apply standardization before running PCA. Standardization normalizes each column by subtracting the column mean and dividing by the standard deviation.
   - Press Start to (re-)run the optimization.
2. Set the color of the displayed points. Set shape, size and label to differentiate between points. If *Label only selection and subset* is ticked, only selected and/or highlighted points will be labelled.
3. Set symbol size and opacity for all data points. Set jittering to randomly disperse data points.
4. *Show color regions* colors the graph by class, while *Show legend* displays a legend on the right. Click and drag the legend to move it.
5. *Select, zoom, pan and zoom to fit* are the options for exploring the graph. The manual selection of data instances works as an angular/square selection tool. Double click to move the projection. Scroll in or out for zoom.
6. If *Send selected automatically* is ticked, changes are communicated automatically. Alternatively, press *Send Selected*.

Preprocessing
-------------

t-SNE uses default preprocessing if necessary. It executes it in the following order:

- continuizes categorical variables (with one feature per value)
- imputes missing values with mean values

To override default preprocessing, preprocess the data beforehand with [Preprocess](../data/preprocess.md) widget.

Examples
--------

The first example is a simple t-SNE plot of *brown-selected* data set. Load *brown-selected* with the [File](../data/file.md) widget. Then connect **t-SNE** to it. The widget will show a 2D map of yeast samples, where samples with similar gene expression profiles will be close together. Select the region, where the gene function is mixed and inspect it in a [Data Table](../data/datatable.md).

![](images/tSNE-Example1.png)

For the second example, use [Single Cell Datasets](https://orangedatamining.com/widget-catalog/single-cell/single_cell_datasets/) widget from the Single Cell add-on to load *Bone marrow mononuclear cells with AML (sample)* data. Then pass it through **k-Means** and select 2 clusters from Silhouette Scores. Ok, it looks like there might be two distinct clusters here.

But can we find subpopulations in these cells? Select a few marker genes with the [Marker Genes](https://orangedatamining.com/widget-catalog/bioinformatics/marker_genes/) widget, for example natural killer cells (NK cells). Pass the marker genes and k-Means results to [Score Cells](https://orangedatamining.com/widget-catalog/single-cell/score_cells/) widget. Finally, add **t-SNE** to visualize the results.

In **t-SNE**, use *Cluster* attribute to color the points and *Score* attribute to set their size. We see that killer cells are nicely clustered together and that t-SNE indeed found subpopulations.

![](images/tSNE-Example2.png)
MDS
===

Multidimensional scaling (MDS) projects items onto a plane fitted to given distances between points.

**Inputs**

- Data: input dataset
- Distances: distance matrix
- Data Subset: subset of instances

**Outputs**

- Selected Data: instances selected from the plot
- Data: dataset with MDS coordinates

[Multidimensional scaling](https://en.wikipedia.org/wiki/Multidimensional_scaling) is a technique which finds a low-dimensional (in our case a two-dimensional) projection of points, where it tries to fit distances between points as well as possible. The perfect fit is typically impossible to obtain since the data is high-dimensional or the distances are not [Euclidean](https://en.wikipedia.org/wiki/Euclidean_distance).

In the input, the widget needs either a dataset or a matrix of distances. When visualizing distances between rows, you can also adjust the color of the points, change their shape, mark them, and output them upon selection.

The algorithm iteratively moves the points around in a kind of a simulation of a physical model: if two points are too close to each other (or too far away), there is a force pushing them apart (or together). The change of the point’s position at each time interval corresponds to the sum of forces acting on it.

![](images/MDS-zoo-stamped.png)

1. The widget redraws the projection during optimization. Optimization is run automatically in the beginning and later by pushing *Start*.
   - **Max iterations**: The optimization stops either when the projection changes only minimally at the last iteration or when a maximum number of iterations has been reached.
   - **Initialization**: PCA (Torgerson) positions the initial points along principal coordinate axes. *Random* sets the initial points to a random position and then readjusts them.
   - **Refresh**: Set how often you want to refresh the visualization. It can be at *Every iteration*, *Every 5/10/25/50 steps* or never (*None*). Setting a lower refresh interval makes the animation more visually appealing, but can be slow if the number of points is high.
2. Defines how the points are visualized. These options are available only when visualizing distances between rows (selected in the [Distances](../unsupervised/distances.md) widget).
   - **Color**: Color of points by attribute (gray for continuous, colored for discrete).
   - **Shape**: Shape of points by attribute (only for discrete).
   - **Size**: Set the size of points (*Same size* or select an attribute) or let the size depend on the value of the continuous attribute the point represents (Stress).
   - **Label**: Discrete attributes can serve as a label.
   - **Symbol size**: Adjust the size of the dots.
   - **Symbol opacity**: Adjust the transparency level of the dots.
   - **Show similar pairs**: Adjust the strength of network lines.
   - **Jitter**: Set [jittering](https://en.wikipedia.org/wiki/Jitter) to prevent the dots from overlapping.
3. Adjust the graph with *Zoom/Select*. The arrow enables you to select data instances. The magnifying glass enables zooming, which can be also done by scrolling in and out. The hand allows you to move the graph around. The rectangle readjusts the graph proportionally.
4. Select the desired output:
      - **Original features only** (input dataset)
      - **Coordinates only** (MDS coordinates)
      - **Coordinates as features** (input dataset + MDS coordinates as regular attributes)
      - **Coordinates as meta attributes** (input dataset + MDS coordinates as meta attributes)
5. Sending the instances can be automatic if *Send selected automatically* is ticked. Alternatively, click *Send selected*.
6. **Save Image** allows you to save the created image either as .svg or .png file to your device.
7. Produce a report.

The MDS graph performs many of the functions of the Visualizations widget. It is in many respects similar to the [Scatter Plot](../visualize/scatterplot.md) widget, so we recommend reading that widget's description as well.

Preprocessing
-------------

When given *Distances* on the input, preprocessing is not applied. When given *Data*, MDS uses default preprocessing if necessary. Preprocessing is executed in the following order:

- continuizes categorical variables (with one feature per value)
- imputes missing values with mean values

To override default preprocessing, preprocess the data beforehand with [Preprocess](../data/preprocess.md) widget.

# Example

The above graphs were drawn using the following simple schema. We used the *iris.tab* dataset. Using the [Distances](../unsupervised/distances.md) widget we input the distance matrix into the **MDS** widget, where we see the *Iris* data displayed in a 2-dimensional plane. We can see the appended coordinates in the [Data Table](../data/datatable.md) widget.

![](images/MDS-Example.png)

# References

Wickelmaier, F. (2003). An Introduction to MDS. Sound Quality Research
Unit, Aalborg University. Available
[here](https://homepages.uni-tuebingen.de/florian.wickelmaier/pubs/Wickelmaier2003SQRU.pdf).
DBSCAN
======

Groups items using the DBSCAN clustering algorithm.

**Inputs**

- Data: input dataset

**Outputs**

- Data: dataset with cluster label as a meta attribute

The widget applies the [DBSCAN clustering](https://en.wikipedia.org/wiki/DBSCAN) algorithm to the data and outputs a new dataset with cluster labels as a meta attribute. The widget also shows the sorted graph with distances to k-th nearest neighbors. With k values set to **Core point neighbors** as suggested in the [methods article](https://www.aaai.org/Papers/KDD/1996/KDD96-037.pdf). This gives the user the idea of an ideal selection for **Neighborhood distance** setting. As suggested by the authors, this parameter should be set to the first value in the first "valley" in the graph.

![](images/DBSCAN.png)

1. **Parameters**:
   - *Core point neighbors*: The number of neighbors for a point to be considered as a core point.
   - *Neighborhood distance*: The maximum distance between two samples for one to be considered as in the neighborhood of the other.
2. Distance metric used in grouping the items (Euclidean, Manhattan, or Cosine). If *Normalize features* is selected, the data will be standardized column-wise (centered to mean and scaled to standard deviation of 1).
3. If *Apply Automatically* is ticked, the widget will commit changes
automatically. Alternatively, click *Apply*.

The graph shows the distance to the k-th nearest neighbor. *k* is
set by the **Core point neighbor** option. With moving the black slider
left and right you can select the right **Neighborhood distance**.

Example
-------

In the following example, we connected the [File](../data/file.md) widget with the Iris dataset to the DBSCAN widget. In the DBSCAN widget, we set **Core points neighbors** parameter to 5. And select the **Neighborhood distance** to the value in the first "valley" in the graph. We show clusters in the [Scatter Plot](../visualize/scatterplot.md) widget.

![](images/DBSCAN-Example.png)
Distances
=========

Computes distances between rows/columns in a dataset.

**Inputs**

- Data: input dataset

**Outputs**

- Distances: distance matrix

The **Distances** widget computes distances between rows or columns in a dataset. By default, the data will be normalized to ensure equal treatment of individual features. Normalization is always done column-wise.

Sparse data can only be used with Euclidean, Manhattan and Cosine metric.

The resulting distance matrix can be fed further to [Hierarchical Clustering](hierarchicalclustering.md) for uncovering groups in the data, to [Distance Map](distancemap.md) or [Distance Matrix](distancematrix.md) for visualizing the distances (Distance Matrix can be quite slow for larger data sets), to [MDS](mds.md) for mapping the data instances using the distance matrix and finally, saved with [Save Distance Matrix](savedistancematrix.md). Distance file can be loaded with [Distance File](distancefile.md).

Distances work well with Orange add-ons, too. The distance matrix can be fed to Network from Distances (Network add-on) to convert the matrix into a graph and to Duplicate Detection (Text add-on) to find duplicate documents in the corpus.

![](images/Distances-stamped.png)

1. Choose whether to measure distances between rows or columns.
2. Choose the *Distance Metric*:
   - [Euclidean](https://en.wikipedia.org/wiki/Euclidean_distance) ("straight line", distance between two points)
   - [Manhattan](https://en.wiktionary.org/wiki/Manhattan_distance) (the sum of absolute differences for all attributes)
   - [Cosine](https://en.wikipedia.org/wiki/Cosine_similarity) (the cosine of the angle between two vectors of an inner product space). Orange computes the cosine distance, which is 1-similarity.
   - [Jaccard](https://en.wikipedia.org/wiki/Jaccard_index) (the size of the intersection divided by the size of the union of the sample sets)
   - [Spearman](https://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient)(linear correlation between the rank of the values, remapped as a distance in a [0, 1] interval)
   - [Spearman absolute](https://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient)(linear correlation between the rank of the absolute values, remapped as a distance in a [0, 1] interval)
   - [Pearson](https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient) (linear correlation between the values, remapped as a distance in a [0, 1] interval)
   - [Pearson absolute](https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient) (linear correlation between the absolute values, remapped as a distance in a [0, 1] interval)
   - [Hamming](https://en.wikipedia.org/wiki/Hamming_distance) (the number of features at which the corresponding values are different)
   - [Bhattacharyya distance](https://en.wikipedia.org/wiki/Bhattacharyya_distance) (Similarity between two probability distributions, not a real distance as it doesn't obey triangle inequality.)

   Normalize the features. Normalization is always done column-wise. Values are zero centered and scaled.
   In case of missing values, the widget automatically imputes the average value of the row or the column.
   The widget works for both numeric and categorical data. In case of categorical data, the distance is 0 if the two values are the same ('green' and 'green') and 1 if they are not ('green' and 'blue').
3. Tick *Apply Automatically* to automatically commit changes to other widgets. Alternatively, press '*Apply*'.

Examples
--------

The first example shows a typical use of the **Distances** widget. We are using the *iris.tab* data from the [File](../data/file.md) widget. We compute distances between data instances (rows) and pass the result to the [Hierarchical Clustering](hierarchicalclustering.md). This is a simple workflow to find groups of data instances.

![](images/Distances-Example1-rows.png)

Alternatively, we can compute distance between columns and find how similar our features are.

![](images/Distances-Example1-columns.png)

The second example shows how to visualize the resulting distance matrix. A nice way to observe data similarity is in a [Distance Map](distancemap.md) or in [MDS](mds.md).

![](images/Distances-Example2.png)
Manifold Learning
=================

Nonlinear dimensionality reduction.

**Inputs**

- Data: input dataset

**Outputs**

- Transformed Data: dataset with reduced coordinates

[Manifold Learning](https://en.wikipedia.org/wiki/Nonlinear_dimensionality_reduction) is a technique which finds a non-linear manifold within the higher-dimensional space. The widget then outputs new coordinates which correspond to a two-dimensional space. Such data can be later visualized with [Scatter Plot](../visualize/scatterplot.md) or other visualization widgets.

![](images/manifold-learning-stamped.png)

1. Method for manifold learning:
   - [t-SNE](http://scikit-learn.org/stable/modules/manifold.html#t-distributed-stochastic-neighbor-embedding-t-sne)
   - [MDS](http://scikit-learn.org/stable/modules/manifold.html#multi-dimensional-scaling-mds), see also [MDS widget](../unsupervised/mds.md)
   - [Isomap](http://scikit-learn.org/stable/modules/manifold.html#isomap)
   - [Locally Linear Embedding](http://scikit-learn.org/stable/modules/manifold.html#locally-linear-embedding)
   - [Spectral Embedding](http://scikit-learn.org/stable/modules/manifold.html#spectral-embedding)
2. Set parameters for the method:
   - t-SNE (distance measures):
     - *Euclidean* distance
     - *Manhattan*
     - *Chebyshev*
     - *Jaccard*
     - *Mahalanobis*
     - *Cosine*
   - MDS (iterations and initialization):
     - *max iterations*: maximum number of optimization interactions
     - *initialization*: method for initialization of the algorithm (PCA or random)
   - Isomap:
     - number of *neighbors*
   - Locally Linear Embedding:
     - *method*:
       - standard
       - modified
       - [hessian eigenmap](http://scikit-learn.org/stable/modules/manifold.html#hessian-eigenmapping)
       - local
     - number of *neighbors*
     - *max iterations*
   - Spectral Embedding:
     - *affinity*:
       - nearest neighbors
       - RFB kernel
3. Output: the number of reduced features (components).
4. If *Apply automatically* is ticked, changes will be propagated automatically. Alternatively, click *Apply*.
5. Produce a report.

**Manifold Learning** widget produces different embeddings for high-dimensional data.

![](images/collage-manifold.png)

From left to right, top to bottom: t-SNE, MDS, Isomap, Locally Linear Embedding and Spectral Embedding.

Preprocessing
-------------

All projections use default preprocessing if necessary. It is executed in the following order:

- continuization of categorical variables (with one feature per value)
- imputation of missing values with mean values

To override default preprocessing, preprocess the data beforehand with [Preprocess](../data/preprocess.md) widget.

Example
-------

*Manifold Learning* widget transforms high-dimensional data into a lower dimensional approximation. This makes it great for visualizing datasets with many features. We used *voting.tab* to map 16-dimensional data onto a 2D graph. Then we used [Scatter Plot](../visualize/scatterplot.md) to plot the embeddings.

![](images/manifold-learning-example.png)
Correspondence Analysis
=======================

Correspondence analysis for categorical multivariate data.

**Inputs**

- Data: input dataset

**Outputs**

- Coordinates: coordinates of all components

[Correspondence Analysis](https://en.wikipedia.org/wiki/Correspondence_analysis) (CA) computes the CA linear transformation of the input data. While it is similar to PCA, CA computes linear transformation on discrete rather than on continuous data.

![](images/CorrespondenceAnalysis-stamped.png)

1. Select the variables you want to see plotted.
2. Select the component for each axis.
3. [Inertia](https://en.wikipedia.org/wiki/Sylvester%27s_law_of_inertia) values (percentage of independence from transformation, i.e. variables are in the same dimension).
4. Produce a report.

Example
-------

Below, is a simple comparison between the **Correspondence Analysis** and [Scatter Plot](../visualize/scatterplot.md) widgets on the *Titanic* dataset. While the [Scatter Plot](../visualize/scatterplot.md) shows fairly well which class and sex had a good survival rate and which one didn't, **Correspondence Analysis** can plot several variables in a 2-D graph, thus making it easy to see the relations between variable values. It is clear from the graph that "no", "male" and "crew" are related to each other. The same goes for "yes", "female" and "first".

![](images/CorrespondenceAnalysis-Example.png)
Louvain Clustering
==================

Groups items using the Louvain clustering algorithm.

**Inputs**

- Data: input dataset

**Outputs**

- Data: dataset with cluster label as a meta attribute
- Graph (with the Network addon): the weighted k-nearest neighbor graph

The widget first converts the input data into a k-nearest neighbor graph. To preserve the notions of distance, the Jaccard index for the number of shared neighbors is used to weight the edges. Finally, a [modularity optimization](https://en.wikipedia.org/wiki/Louvain_Modularity) community detection algorithm is applied to the graph to retrieve clusters of highly interconnected nodes. The widget outputs a new dataset in which the cluster label is used as a meta attribute.

![](images/LouvainClustering.png)

1. Information on the number of clusters found.
2. **Preprocessing**:
   - *Normalize data*: Center to mean and scale to standard deviation of 1.
   - *Apply PCA preprocessing*: PCA processing is typically applied to the original data to remove noise (see [PCA](PCA.md) widget).
   - *PCA Components*: number of principal components used.
3. **Graph parameters**:
   - *Distance metric*: The distance metric is used for finding specified number of nearest neighbors (Euclidean, Manhattan, Cosine).
   - *k neighbors*: The number of nearest neighbors to use to form the KNN graph.
   - *Resolution* is a parameter for the Louvain community detection algorithm that affects the size of the recovered clusters. Smaller resolutions recover smaller clusters and therefore a larger number of them, while, conversely, larger values recover clusters containing more data points.
4. When *Apply Automatically* is ticked, the widget will automatically communicate all changes. Alternatively, click *Apply*.

Preprocessing
-------------

Louvain Clustering uses default preprocessing if necessary. It executes it in the following order:

- continuizes categorical variables (with one feature per value)
- imputes missing values with mean values

To override default preprocessing, preprocess the data beforehand with [Preprocess](../data/preprocess.md) widget.

Example
-------

*Louvain Clustering* converts the dataset into a graph, where it finds highly interconnected nodes. In the example below, we used the iris data set from the [File](../data/file.md) widget, then passed it to **Louvain Clustering**, which found 4 clusters. We plotted the data with [Scatter Plot](../visualize/scatterplot.md), where we colored the data points according to clusters labels.

![](images/LouvainClustering-Example.png)

We can visualize the graph itself using the **Network Explorer** from the Network addon.

References
----------

Blondel, Vincent D., et al. "[Fast unfolding of communities in large networks.](https://arxiv.org/abs/0803.0476)" Journal of statistical mechanics: theory and experiment 2008.10 (2008): P10008.

Lambiotte, Renaud, J-C. Delvenne, and Mauricio Barahona. "Laplacian dynamics and multiscale modular structure in networks." arXiv preprint, [arXiv:0812](https://arxiv.org/abs/0812.1770).1770 (2008).
Save Distance Matrix
====================

Saves a distance matrix.

If the file is saved to the same directory as the workflow or in the subtree of that directory, the widget remembers the relative path. Otherwise it will store an absolute path, but disable auto save for security reasons.

**Inputs**

- Distances: distance matrix

![](images/SaveDistanceMatrix-stamped.png)

1. By clicking *Save*, you choose from previously saved distance matrices. Alternatively, tick the box on the left side of the *Save* button and changes will be communicated automatically.
2. By clicking *Save as*, you save the distance matrix to your computer, you only need to enter the name of the file and click *Save*. The distance matrix will be saved as type *.dst*.

Example
-------

In the snapshot below, we used the [Distance Transformation](../unsupervised/distancetransformation.md) widget to transform the distances in the *Iris* dataset. We then chose to save the transformed version to our computer, so we could use it later on. We decided to output all data instances. You can choose to output just a minor subset of the data matrix. Pairs are marked automatically. If you wish to know what happened to our changed file, see [Distance File](../unsupervised/distancefile.md).

![](images/SaveDistanceMatrix-Example.png)
Merge Data
==========

Merges two datasets, based on values of selected attributes.

**Inputs**

- Data: input dataset
- Extra Data: additional dataset

**Outputs**

- Data: dataset with features added from extra data

The **Merge Data** widget is used to horizontally merge two datasets, based on the values of selected attributes (columns). In the input, two datasets are required, data and extra data. Rows from the two data sets are matched by the values of pairs of attributes, chosen by the user. The widget produces one output. It corresponds to the instances from the input data to which attributes (columns) from input extra data are appended.

If the selected attribute pair does not contain unique values (in other words, the attributes have duplicate values), the widget will give a warning. Instead, one can match by more than one attribute. Click on the plus icon to add the attribute to merge on. The final result has to be a unique combination for each individual row.

![](images/Merge-Data-stamped.png)

1. Information on main data.
2. Information on data to append.
3. Merging type:
   - **Append columns from Extra Data** outputs all rows from the Data, augmented by the columns in the Extra Data. Rows without matches are retained, even where the data in the extra columns are missing.
   - **Find matching pairs of rows** outputs rows from the Data, augmented by the columns in the Extra Data. Rows without matches are removed from the output.
   - **Concatenate tables** treats both data sources symmetrically. The output is similar to the first option, except that non-matched values from Extra Data are appended at the end.
4. List of attributes from Data input.
5. List of attributes from Extra Data input.
6. Produce a report.

Merging Types
-------------

#####Append Columns from Extra Data (left join)

Columns from the Extra Data are added to the Data. Instances with no matching rows will have missing values added.

For example, the first table may contain city names and the second would be a list of cities and their coordinates. Columns with coordinates would then be appended to the data with city names. Where city names cannot be matched, missing values will appear.

In our example, the first Data input contained 6 cities, but the Extra Data did not provide Lat and Lon values for Bratislava, so the fields will be empty.

![](images/MergeData_Append.png)

#####Find matching pairs of rows (inner join)

Only those rows that are matched will be present on the output, with the Extra Data columns appended. Rows without matches are removed.

In our example, Bratislava from the Data input did not have Lat and Lon values, while Belgrade from the Extra Data could not be found in the City column we were merging on. Hence both instances are remove - only the intersection of instances is sent to the output.

![](images/MergeData_Intersection.png)

#####Concatenate tables (outer join)

The rows from both the Data and the Extra Data will be present on the output. Where rows cannot be matched, missing values will appear.

In our example, both Bratislava and Belgrade are now present. Bratislava will have missing Lat and Lon values, while Belgrade will have a missing Population value.

![](images/MergeData_Concatenate.png)

#####Row index

Data will be merged in the same order as they appear in the table. Row number 1 from the Data input will be joined with row number 1 from the Extra Data input. Row numbers are assigned by Orange based on the original order of the data instances.

#####Instance ID

This is a more complex option. Sometimes, data in transformed in the analysis and the domain is no longer the same. Nevertheless, the original row indices are still present in the background (Orange remembers them). In this case one can merge on instance ID. For example if you transformed the data with PCA, visualized it in the Scatter Plot, selected some data instances and now you wish to see the original information of the selected subset. Connect the output of Scatter Plot to Merge Data, add the original data set as Extra Data and merge by Instance ID.

![](images/MergeData-InstanceID.png)

#####Merge by two or more attributes

Sometimes our data instances are unique with respect to a combination of columns, not a single column. To merge by more than a single column, add the *Row matching* condition by pressing plus next to the matching condition. To remove it, press the x.

In the below example, we are merging by *student* column and *class* column.

![](images/MergeData-multiple.png)

Say we have two data sets with student names and the class they're in. The first data set has students' grades and the second on the elective course they have chosen. Unfortunately, there are two Jacks in our data, one from class A and the other from class B. Same for Jane.

To distinguish between the two, we can match rows on both, the student's name and her class.

![](images/MergeData-multiple2.png)

Examples
--------

Merging two datasets results in appending new attributes to the original file, based on a selected common attribute. In the example below, we wanted to merge the **zoo.tab** file containing only factual data with [zoo-with-images.tab](http://file.biolab.si/datasets/zoo-with-images.tab) containing images. Both files share a common string attribute *names*. Now, we create a workflow connecting the two files. The *zoo.tab* data is connected to **Data** input of the **Merge Data** widget, and the *zoo-with-images.tab* data to the **Extra Data** input. Outputs of the **Merge Data** widget is then connected to the [Data Table](../data/datatable.md) widget. In the latter, the **Merged Data** channels are shown, where image attributes are added to the original data.

![](images/MergeData-Example.png)

The case where we want to include all instances in the output, even those where no match by attribute *names* was found, is shown in the following workflow.

![](images/MergeData-Example2.png)

The third type of merging is shown in the next workflow. The output consists of both inputs, with unknown values assigned where no match was found.

![](images/MergeData-Example3.png)
CSV File Import
===============

Import a data table from a CSV formatted file.

**Outputs**

- Data: dataset from the .csv file
- Data Frame: pandas DataFrame object

The **CSV File Import** widget reads comma-separated files and sends the dataset to its output channel. File separators can be commas, semicolons, spaces, tabs or manually-defined delimiters. The history of most recently opened files is maintained in the widget.

*Data Frame* output can be used in the [Python Script](../data/pythonscript.md) widget by connecting it to the `in_object` input (e.g. `df = in_object`). Then it can be used a regular DataFrame.

### Import Options

The import window where the user sets the import parameters. Can be re-opened by pressing *Import Options* in the widget.

Right click on the column name to set the column type. Right click on the row index (on the left) to mark a row as a header, skipped or a normal data row.

![](images/CSVFileImport-ImportOptions-stamped.png)

1. File encoding. Default is UTF-8. See Encoding subchapter for details.
2. Import settings:
   - *Cell delimiter*:
      - Tab
      - Comma
      - Semicolon
      - Space
      - Other (set the delimiter in the field to the right)
   - *Quote character*: either " or '. Defines what is considered a text.
   - *Number separators*:
      - Grouping: delimiters for thousands, e.g. 1,000
      - Decimal: delimiters for decimals, e.g. 1.234
3. Column type: select the column in the preview and set its type. Column type can be set also by right-clicking on the selected column.
   - *Auto*: Orange will automatically try to determine column type. (default)
   - *Numeric*: for continuous data types, e.g. (1.23, 1.32, 1.42, 1.32)
   - *Categorical*: for discrete data types, e.g. (brown, green, blue)
   - *Text*: for string data types, e.g. (John, Olivia, Mike, Jane)
   - *Datetime*: for time variables, e.g. (1970-01-01)
   - *Ignore*: do not output the column.
4. Pressing *Reset* will return the settings to the previously set state (saved by pressing OK in the Import Options dialogue). *Restore Defaults* will set the settings to their default values. *Cancel* aborts the import, while *OK* imports the data and saves the settings.

### Widget

The widget once the data is successfully imported.

![](images/CSVFileImport-widget-stamped.png)

1. The folder icon opens the dialogue for import the local .csv file. It can be used to either load the first file or change the existing file (load new data). The *File* dropdown stores paths to previously loaded data sets.
2. Information on the imported data set. Reports on the number of instances (rows), variables (features or columns) and meta variables (special columns).
3. *Import Options* re-opens the import dialogue where the user can set delimiters, encodings, text fields and so on. *Cancel* aborts data import. *Reload* imports the file once again, adding to the data any changes made in the original file.

### Encoding

The dialogue for settings custom encodings list in the Import Options - Encoding dropdown. Select *Customize Encodings List...* to change which encodings appear in the list. To save the changes, simply close the dialogue. Closing and reopening Orange (even with Reset widget settings) will not re-set the list. To do this, press *Restore Defaults*. To have all the available encodings in the list, press *Select all*.

![](images/CSVFileImport-encodings.png)

Example
-------

**CSV File Import** works almost exactly like the [File](../data/file.md) widget, with the added options for importing different types of .csv files. In this workflow, the widget read the data from the file and sends it to the [Data Table](../data/datatable.md) for inspection.

![](images/CSVFileImport-Example.png)

File
====

Reads attribute-value data from an input file.

**Outputs**

- Data: dataset from the file

The **File** widget [reads the input data file](../../loading-your-data/index.md) (data table with data instances) and sends the dataset to its output channel. The history of most recently opened files is maintained in the widget. The widget also includes a directory with sample datasets that come pre-installed with Orange.

The widget reads data from Excel (**.xlsx**), simple tab-delimited (**.txt**), comma-separated files (**.csv**) or URLs. For other formats see Other Formats section below.

![](images/File-stamped.png)

1. Browse through previously opened data files, or load any of the sample ones.  
2. Browse for a data file.
3. Reloads currently selected data file.
4. Insert data from URL addresses, including data from Google Sheets.
5. Information on the loaded dataset: dataset size, number and types of data features.
6. Additional information on the features in the dataset. Features can be edited by double-clicking on them. The user can change the attribute names, select the type of variable per each attribute (*Continuous*, *Nominal*, *String*, *Datetime*), and choose how to further define the attributes (as *Features*, *Targets* or *Meta*). The user can also decide to ignore an attribute.
7. Browse documentation datasets.
8. Produce a report.

Example
-------

Most Orange workflows would probably start with the **File** widget. In the schema below, the widget is used to read the data that is sent to both the [Data Table](../data/datatable.md) and the [Box Plot](../visualize/boxplot.md) widget.

![](images/File-Workflow.png)

### Loading your data

- Orange can import any comma, .xlsx or tab-delimited data file or URL. Use the **File** widget and then, if needed, select class and meta attributes.
- To specify the domain and the type of the attribute, attribute names can be preceded with a label followed by a hash. Use c for class and m for meta attribute, i to ignore a column, and C, D, S for continuous, discrete and string attribute types. Examples: C#mpg, mS#name, i#dummy.
- Orange's native format is a tab-delimited text file with three header rows. The first row contains attribute names, the second the type (*continuous*, *discrete* or *string*), and the third the optional element (*class*, *meta* or *time*).

![](images/spreadsheet-simple-head1.png)

Read more on loading your data [here](../../loading-your-data/index.md).

### Other Formats

Supported formats and the widgets to load them:

- distance matrix: [Distance File](../unsupervised/distancefile.md)
- predictive model: [Load Model](../model/loadmodel.md)
- network: Network File from Network add-on
- images: Import Images from Image Analytics add-on
- text/corpus: Corpus or Import Documents from Text add-on
- single cell data: Load Data from Single Cell add-on
- several spectroscopy files: Multifile from Spectroscopy add-on
Melt
=========

Transform [wide data to narrow](https://en.wikipedia.org/wiki/Wide_and_narrow_data).

**Inputs**

- Data: wide data table

**Outputs**

- Data: narrow data table

The **Melt** widget receives a dataset in the more common wide format and outputs a table of (row_id, variable, value) triplets.


![](images/Melt-Default-stamped.png)

1. Select the variable used as id. The widget offers only columns without duplicated values. Alternatively, row number can be used as id.
2. Select whether to include non-numeric variables, and whether to exclude zero values.
3. Set the names of the columns with name of the variable ("item") and the corresponding value.

Example
-------

In the following workflow we play with the Zoo data set, in which we convert all variables to numeric by treating them as ordinal. All variables except the number of legs boolean (e.g. the animal lays or does not lay eggs), so a value of 1 will correspond to an animal having a particular feature. In data table we select all rows (Ctrl-A or Cmd-A) and deselect the duplicate description of the frog in order to avoid duplicate values in the "name" column.

We pass it to Melt, where we designate the name as the row id, and discard zero values. The resulting table has multiple rows for each animal: one for each of animals features.

An interesting immediate use for this is to pass this data to Distributions and see what are the most and the least common features of animals.

![](images/Melt-Workflow.png)

In the next example we show how shuffling class values influences model performance on the same dataset as above.

![](images/Melt-Distribution.png)
Correlations
============

Compute all pairwise attribute correlations.

**Inputs**

- Data: input dataset

**Outputs**

- Data: input dataset
- Features: selected pair of features
- Correlations: data table with correlation scores

**Correlations** computes Pearson or Spearman correlation scores for all pairs of features in a dataset. These methods can only detect monotonic relationship.

![](images/Correlations-stamped.png)

1. Correlation measure:
   - Pairwise [Pearson](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) correlation.
   - Pairwise [Spearman](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient) correlation.
2. Filter for finding attribute pairs.
3. A list of attribute pairs with correlation coefficient. Press *Finished* to stop computation for large datasets.
4. Access widget help and produce report.

Example
-------

Correlations can be computed only for numeric (continuous) features, so we will use *housing* as an example data set. Load it in the [File](file.md) widget and connect it to **Correlations**. Positively correlated feature pairs will be at the top of the list and negatively correlated will be at the bottom.

![](images/Correlations-links.png)

Go to the most negatively correlated pair, DIS-NOX. Now connect [Scatter Plot](../visualize/scatterplot.md) to **Correlations** and set two outputs, Data to Data and Features to Features. Observe how the feature pair is immediately set in the scatter plot. Looks like the two features are indeed negatively correlated.

![](images/Correlations-Example.png)
Data Info
=========

Displays information on a selected dataset.

**Inputs**

- Data: input dataset

A simple widget that presents information on dataset size, features,
targets, meta attributes, and location.

![](images/data-info-stamped.png)

1. Information on dataset size
2. Information on discrete and continuous features
3. Information on targets
4. Information on meta attributes
5. Information on where the data is stored
6. Produce a report.

Example
-------

Below, we compare the basic statistics of two **Data Info** widgets - one with information on the entire dataset and the other with information on the (manually) selected subset from the [Scatter Plot](../visualize/scatterplot.md) widget. We used the *Iris* dataset.

![](images/DataInfo-Example.png)
Feature Constructor
===================

Add new features to your dataset.

**Inputs**

- Data: input dataset

**Outputs**

- Data: dataset with additional features

The **Feature Constructor** allows you to manually add features (columns) into your dataset. The new feature can be a computation of an existing one or a combination of several (addition, subtraction, etc.). You can choose what type of feature it will be (discrete, continuous or string) and what its parameters are (name, value, expression). For continuous variables you only have to construct an expression in Python.

![](images/feature-constructor1-stamped.png)

1. List of constructed variables
2. Add or remove variables
3. New feature name
4. Expression in Python
5. Select a feature
6. Select a function
7. Produce a report
8. Press *Send* to communicate changes

For discrete variables, however, there's a bit more work. First add or remove the values you want for the new feature. Then select the base value and the expression. In the example below, we have constructed an expression with 'if lower than' and defined three conditions; the program ascribes 0 (which we renamed to lower) if the original value is lower than 6, 1 (mid) if it is lower than 7 and 2 (higher) for all the other values. Notice that we use an underscore for the feature name (e.g. petal\_length).

![](images/feature-constructor2-stamped.png)

1. List of variable definitions
2. Add or remove variables
3. New feature name
4. Expression in Python
5. Select a feature
6. Select a function
7. Assign values
8. Produce a report
9. Press *Send* to communicate changes

Example
-------

With the **Feature Constructor** you can easily adjust or combine existing features into new ones. Below, we added one new discrete feature to the *Titanic* dataset. We created a new attribute called *Financial status* and set the values to be *rich* if the person belongs to the first class (status = first) and *not rich* for everybody else. We can see the new dataset with [Data Table](../data/datatable.md) widget.

![](images/FeatureConstructor-Example.png)

Hints
-----

If you are unfamiliar with Python math language, here's a quick introduction.

- +, - to add, subtract
- \* to multiply
- / to divide
- % to divide and return the remainder
- \*\* for exponent (for square root square by 0.5)
- // for floor division
- <, >, <=, >= less than, greater than, less or equal, greater or equal
- == for equal
- != for not equal

As in the example: (*value*) if (*feature name*) < (*value*), else (*value*) if (*feature name*) < (*value*), else (*value*)

[Use value 1 if feature is less than specified value, else use value 2 if feature is less than specified value 2, else use value 3.]

See more [here](http://www.tutorialspoint.com/python/python_basic_operators.htm).
Apply Domain
============

Given dataset and template transforms the dataset.

**Inputs**

- Data: input dataset
- Template Data: template for transforming the dataset

**Outputs**

- Transformed Data: transformed dataset

**Apply Domain** maps new data into a transformed space. For example, if we transform some data with PCA and wish to observe new data in the same space, we can use Apply Domain to map the new data into the PCA space created from the original data.

![](images/ApplyDomain.png)

The widget receives a dataset and a template dataset used to transform the dataset.

Side note
--------

Domain transformation works by using information from the template data. For example, for PCA, Components are not enough. Transformation requires information on the center of each column, variance (if the data is normalized), and if and how the data was preprocessed (continuized, imputed, etc.).

Example
-------

We will use iris data from the [File](../data/file.md) widget for this example. To create two separate data sets, we will use [Select Rows](../data/selectrows.md) and set the condition to *iris is one of iris-setosa, iris-versicolor*. This will output a data set with a 100 rows, half of them belonging to iris-setosa class and the other half to iris-versicolor.

We will transform the data with [PCA](../unsupervised/PCA.md) and select the first two components, which explain 96% of variance. Now, we would like to apply the same preprocessing on the 'new' data, that is the remaining 50 iris virginicas. Send the unused data from **Select Rows** to **Apply Domain**. Make sure to use the *Unmatched Data* output from **Select Rows** widget. Then add the *Transformed data* output from **PCA**.

**Apply Domain** will apply the preprocessor to the new data and output it. To add the new data to the old data, use [Concatenate](../data/concatenate.md). Use *Transformed Data* output from **PCA** as *Primary Data* and *Transformed Data* from **Apply Domain** as *Additional Data*.

Observe the results in a [Data Table](../data/datatable.md) or in a [Scatter Plot](../visualize/scatterplot.md) to see the new data in relation to the old one.

![](images/ApplyDomain-Example.png)
Data Table
==========

Displays attribute-value data in a spreadsheet.

**Inputs**

- Data: input dataset

**Outputs**

- Selected Data: instances selected from the table

The **Data Table** widget receives one or more datasets in its input and presents them as a spreadsheet. Data instances may be sorted by attribute values. The widget also supports manual selection of data instances.

![](images/DataTable-stamped.png)

1. The name of the dataset (usually the input data file). Data
   instances are in rows and their attribute values in columns. In this
   example, the dataset is sorted by the attribute "sepal length".
2. Info on current dataset size and number and types of attributes
3. Values of continuous attributes can be visualized with bars; colors
   can be attributed to different classes.
4. Data instances (rows) can be selected and sent to the widget's output
   channel.
5. Use the *Restore Original Order* button to reorder data instances after
   attribute-based sorting.
6. Produce a report.
7. While auto-send is on, all changes will be automatically communicated
   to other widgets. Otherwise, press *Send Selected Rows*.

Example
-------

We used two [File](../data/file.md) widgets to read the *Iris* and *Glass* dataset (provided in Orange distribution), and send them to the **Data Table** widget.

![](images/DataTable-Schema.png)

Selected data instances in the first **Data Table** are passed to the second **Data Table**. Notice that we can select which dataset to view (iris or glass). Changing from one dataset to another alters the communicated selection of data instances if *Commit on any change* is selected.

![](images/DataTable-Example.png)
Edit Domain
===========

Rename features and their values.

**Inputs**

- Data: input dataset

**Outputs**

- Data: dataset with edited domain

This widget can be used to edit/change a dataset's domain - rename features, rename or merge values of categorical features, add a categorical value, and assign labels.

![](images/EditDomain-stamped.png)

1. All features (including meta attributes) from the input dataset are listed in the *Variables* list. Selecting one feature displays an editor on the right.
2. Editing options:
   - Change the name of the feature.
   - Change the type of the feature. For example, convert a string variable to categorical.
   - *Unlink variable from its source variable*. This option removes existing computation for a variable (say for Cluster how clustering was computed), making it 'plain'. This enables merging variables with same names in [Merge Data](../data/mergedata.md).
   - Change the value names for discrete features in the *Values* list box. Double-click to edit the name.
   - Add, remove or edit additional feature annotations in the *Labels* box. Add a new label with the + button and add the *Key* and *Value* for the new entry. Key will be displayed in the top left corner of the [Data Table](../data/datatable.md), while values will appear below the specified column. Remove an existing label with the - button.
3. Reorder or merge values of categorical features. To reorder the values (for example, to display them in [Distributions](../visualize/distributions.md), use the up and down keys at the bottom of the box. To add or remove a value, use + and - buttons. Select two or more variables and click = to merge them into a single value. Use the M button to merge variables on condition.
4. Rename the output table. Useful for displaying table names in [Venn Diagram](../visualize/venndiagram.md).
5. To revert the changes made to the selected feature, press the *Reset Selected* button while the feature is selected in the *Variables* list. Pressing *Reset All* will remove all the changes to the domain. Press *Apply* to send the new domain to the output.

**Merging options**

![](images/EditDomain-merge.png)

- *Group selected values*: selected cateogorical values become a single variable.
- *Group values with less than N occurrences*: values which appear less than N times in the data, will be grouped into a single value.
- *Group values with less than % occurrences*: values which appear less then X % of the time in the data, will be grouped into a single value.
- *Group all except N most frequent values*: all values but the N most frequent will be grouped into a single variable.
- *New value name*: the name of the grouped value.

Example
-------

Below, we demonstrate how to simply edit an existing domain. We selected the *heart_disease.tab* dataset and edited the *gender* attribute. Where in the original we had the values *female* and *male*, we changed it into *F* for female and *M* for male. Then we used the down key to switch the order of the variables. Finally, we added a label to mark that the attribute is binary. We can observe the edited data in the [Data Table](../data/datatable.md) widget.

![](images/EditDomain-Example.png)
Purge Domain
============

Removes unused attribute values and useless attributes, sorts the remaining values.

**Inputs**

- Data: input dataset

**Outputs**

- Data: filtered dataset

Definitions of nominal attributes sometimes contain values which don’t appear in the data. Even if this does not happen in the original data, filtering the data, selecting exemplary subsets and alike can remove all examples for which the attribute has some particular value. Such values clutter data presentation, especially various visualizations, and should be removed.

After purging an attribute, it may become single-valued or, in extreme case, have no values at all (if the value of this attribute was undefined for all examples). In such cases, the attribute can be removed.

A different issue is the order of attribute values: if the data is read from a file in a format in which values are not declared in advance, they are sorted “in order of appearance”. Sometimes we would prefer to have them sorted alphabetically.

![](images/PurgeDomain-stamped.png)

1. Purge attributes.
2. Purge classes.
3. Purge meta attributes.
4. Information on the filtering process.
5. Produce a report.
6. If *Apply automatically* is ticked, the widget will output data at
   each change of widget settings.

Such purification is done by the widget **Purge Domain**. Ordinary attributes and class attributes are treated separately. For each, we can decide if we want the values sorted or not. Next, we may allow the widget to remove attributes with less than two values or remove the class attribute if there are less than two classes. Finally, we can instruct the widget to check which values of attributes actually appear in the data and remove the unused values. The widget cannot remove values if it is not allowed to remove the attributes, since having attributes without values makes no sense.

The new, reduced attributes get the prefix “R”, which distinguishes them from the original ones. The values of new attributes can be computed from the old ones, but not the other way around. This means that if you construct a classifier from the new attributes, you can use it to classify the examples described by the original attributes. But not the opposite: constructing a classifier from the old attributes and using it on examples described by the reduced ones won’t work. Fortunately, the latter is seldom the case. In a typical setup, one would explore the data, visualize it, filter it, purify it… and then test the final model on the original data.

Example
-------

The **Purge Domain** widget would typically appear after data filtering, for instance when selecting a subset of visualized examples.

In the above schema, we play with the *adult.tab* dataset: we visualize it and select a portion of the data, which contains only four out of the five original classes. To get rid of the empty class, we put the data through **Purge Domain** before going on to the [Box Plot](../visualize/boxplot.md) widget. The latter shows only the four classes which are in the **Purge Data** output. To see the effect of data purification, uncheck *Remove unused class variable values* and observe the effect this has on [Box Plot](../visualize/boxplot.md).

![](images/PurgeDomain-example.png)
Select by Data Index
====================

Match instances by index from data subset.

**Inputs**

- Data: reference data set
- Data Subset: subset to match

**Outputs**

- Matching data: subset from reference data set that matches indices from subset data
- Unmatched data: subset from reference data set that does not match indices from subset data
- Annotated data: reference data set with an additional column defining matches

**Select by Data Index** enables matching the data by indices. Each row in a data set has an index and given a subset, this widget can match these indices to indices from the reference data. Most often it is used to retrieve the original data from the transformed data (say, from PCA space).

![](images/Select-by-Data-Index-stamped.png)

1. Information on the reference data set. This data is used as index reference.
2. Information on the data subset. The indices of this data set are used to find matching data in the reference data set. Matching data are on the output by default.

Example
-------

A typical use of **Select by Data Index** is to retrieve the original data after a transformation. We will load *iris.tab* data in the [File](../data/file.md) widget. Then we will transform this data with [PCA](../unsupervised/PCA.md). We can project the transformed data in a [Scatter Plot](../visualize/scatterplot.md), where we can only see PCA components and not the original features.

Now we will select an interesting subset (we could also select the entire data set). If we observe it in a [Data Table](../data/datatable.md), we can see that the data is transformed. If we would like to see this data with the original features, we will have to retrieve them with **Select by Data Index**.

Connect the original data and the subset from [Scatter Plot](../visualize/scatterplot.md) to **Select by Data Index**. The widget will match the indices of the subset with the indices of the reference (original) data and output the matching reference data. A final inspection in another [Data Table](../data/datatable.md) confirms the data on the output is from the original data space.

![](images/Select-by-Data-Index-Example1.png)
Concatenate
===========

Concatenates data from multiple sources.

**Inputs**

- Primary Data: data set that defines the attribute set
- Additional Data: additional data set

**Outputs**

- Data: concatenated data

The widget concatenates multiple sets of instances (data sets). The merge is “vertical”, in a sense that two sets of 10 and 5 instances yield a new set of 15 instances.

![](images/Concatenate-stamped.png)

1. Set the attribute merging method.
2. Add the identification of source data sets to the output data set.
3. Produce a report.
4. If *Apply automatically* is ticked, changes are communicated automatically. Otherwise, click *Apply*.

If one of the tables is connected to the widget as the primary table, the resulting table will contain its own attributes. If there is no primary table, the attributes can be either a union of all attributes that appear in the tables specified as *Additional Tables*, or their intersection, that is, a list of attributes common to all the connected tables.

Example
-------

As shown below, the widget can be used for merging data from two separate files. Let's say we have two data sets with the same attributes, one containing instances from the first experiment and the other instances from the second experiment and we wish to join the two data tables together. We use the **Concatenate** widget to merge the data sets by attributes (appending new rows under existing attributes).

Below, we used a modified *Zoo* data set. In the [first](http://file.biolab.si/datasets/zoo-first.tab) [File](../data/file.md) widget, we loaded only the animals beginning with the letters A and B and in the [second](http://file.biolab.si/datasets/zoo-second.tab) one only the animals beginning with the letter C. Upon concatenation, we observe the new data in the [Data Table](../data/datatable.md) widget, where we see the complete table with animals from A to C.

![](images/Concatenate-Example.png)
Unique
======

Remove duplicated data instances.

**Inputs**

- Data: data table

**Outputs**

- Data: data table without duplicates

The widget removes duplicated data instances. The user can choose a subset of observed variables, so two instances are considered as duplicates although they may differ in values of other, ignored variables.

![](images/Unique-stamped.png)

1. Select the variables that are considered in comparing data instances.
2. Data instance that is kept. The options are to use the first, last, middle or random instance, or to keep none, that is, to remove duplicated instances altogether.

Example
-------

Data set *Zoo* contains two frogs. This workflow keeps only one by removing instances with the same names.

![](images/Unique-Example.png)
SQL Table
=========

Reads data from an SQL database.

**Outputs**

- Data: dataset from the database

The **SQL** widget accesses data stored in an SQL database. It can connect to PostgreSQL (requires [psycopg2](http://initd.org/psycopg/) module) or [SQL Server](https://www.microsoft.com/en-us/sql-server/) (requires [pymssql](http://pymssql.org/en/stable/) module).

To handle large databases, Orange attempts to execute a part of the computation in the database itself without downloading the data. This only works with PostgreSQL database and requires quantile and tsm_system_time [extensions](https://github.com/biolab/orange3/wiki/Installation-of-SQL-extensions) installed on server. If these extensions are not installed, the data will be downloaded locally.

![](images/SQLTable-stamped.png)

1. Database type (can be either PostgreSQL or MSSQL).
2. Host name.
3. Database name.
4. Username.
5. Password.
6. Press the blue button to connect to the database. Then select the table in the dropdown.
7. *Auto-discover categorical variables* will cast INT and CHAR columns with less than 20 distinct values as categorical variables (finding all distinct values can be slow on large tables). When not selected, INT will be treated as numeric and CHAR as text. *Download to local memory* downloads the selected table to your local machine.

##Installation Instructions

###PostgreSQL

Install the backend.

    pip install psycopg2

Alternatively, you can follow [these instructions](https://blog.biolab.si/2018/02/16/how-to-enable-sql-widget-in-orange/) for installing the backend.

If the installation of `psycopg2` fails, follow to instructions in the error message you get (it explains how to solve the error) or install an already compiled version of `psycopg2-binary` package:

    pip install psycopg2-binary

Note: `psycopg2-binary` comes with own versions of a few C libraries, among which libpq and libssl, which will be used regardless of other libraries available on the client: upgrading the system libraries will not upgrade the libraries used by psycopg2. Please build psycopg2 from source if you want to maintain binary upgradeability.

[Install the extensions](https://github.com/biolab/orange3/wiki/Installation-of-SQL-extensions). [optional]

###MSSQL

Install the backend.

    pip install pymssql

If you are encountering issues, follow [these instructions](https://github.com/biolab/orange3/wiki/Installation-of-SQL-extensions#mssql).

##Example

Here is a simple example on how to use the **SQL Table** widget. Place the widget on the canvas, enter your database credentials and connect to your database. Then select the table you wish to analyse.

Connect **SQL Table** to [Data Table](../data/datatable.md) widget to inspect the output. If the table is populated, your data has transferred correctly. Now, you can use the **SQL Table** widget in the same way as the [File](../data/file.md) widget.

![](images/SQLTable-Example.png)
Outliers
========

Outlier detection widget.

**Inputs**

- Data: input dataset

**Outputs**

- Outliers: instances scored as outliers
- Inliers: instances not scored as outliers
- Data: input dataset appended *Outlier* variable

The **Outliers** widget applies one of the four methods for outlier detection. All methods apply classification to the dataset. *One-class SVM with non-linear kernels (RBF)* performs well with non-Gaussian distributions, while *Covariance estimator* works only for data with Gaussian distribution. One efficient way to perform outlier detection on moderately high dimensional datasets is to use the *Local Outlier Factor* algorithm. The algorithm computes a score reflecting the degree of abnormality of the observations. It measures the local density deviation of a given data point with respect to its neighbors. Another efficient way of performing outlier detection in high-dimensional datasets is to use random forests (*Isolation Forest*).

![](images/Outliers-stamped.png)

1. Method for outlier detection:
   - [One Class SVM](http://scikit-learn.org/stable/modules/generated/sklearn.svm.OneClassSVM.html)
   - [Covariance Estimator](http://scikit-learn.org/stable/modules/generated/sklearn.covariance.EllipticEnvelope.html)
   - [Local Outlier Factor](http://scikit-learn.org/stable/modules/generated/sklearn.neighbors.LocalOutlierFactor.html)
   - [Isolation Forest](http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.IsolationForest.html)
2. Set parameters for the method:
   - **One class SVM with non-linear kernel (RBF)**: classifies data as similar or different from the core class:
      - *Nu* is a parameter for the upper bound on the fraction of training errors and a lower bound of the fraction of support vectors
      - *Kernel coefficient* is a gamma parameter, which specifies how much influence a single data instance has
    - **Covariance estimator**: fits ellipsis to central points with Mahalanobis distance metric:
      - *Contamination* is the proportion of outliers in the dataset
      - *Support fraction* specifies the proportion of points included in the estimate
   - **Local Outlier Factor**: obtains local density from the k-nearest neighbors:
      - *Contamination* is the proportion of outliers in the dataset
      - *Neighbors* represents number of neighbors
      - *Metric* is the distance measure
   - **Isolation Forest**: isolates observations by randomly selecting a feature and then randomly selecting a split value between the maximum and minimum values of the selected feature:
     - *Contamination* is the proportion of outliers in the dataset
     - *Replicabe training* fixes random seed
3. If *Apply automatically* is ticked, changes will be propagated automatically. Alternatively, click *Apply*.
4. Produce a report.
5. Number of instances on the input, followed by number of instances scored as inliers.


Example
-------

Below is an example of how to use this widget. We used subset (*versicolor* and *virginica* instances) of the *Iris* dataset to detect the outliers. We chose the *Local Outlier Factor* method, with *Euclidean* distance. Then we observed the annotated instances in the [Scatter Plot](../visualize/scatterplot.md) widget. In the next step we used the *setosa* instances to demonstrate novelty detection using [Apply Domain](../data/applydomain.md) widget. After concatenating both outputs we examined the outliers in the *Scatter Plot (1)*.

![](images/Outliers-Example.png)
Create Class
============

Create class attribute from a string attribute.

**Inputs**

- Data: input dataset

**Outputs**

- Data: dataset with a new class variable

**Create Class** creates a new class attribute from an existing discrete or string attribute. The widget matches the string value of the selected attribute and constructs a new user-defined value for matching instances.

![](images/CreateClass-stamped.png)

1. The attribute the new class is constructed from.
2. Matching:
   - Name: the name of the new class value
   - Substring: regex-defined substring that will match the values from the above-defined attribute
   - Instances: the number of instances matching the substring
   - Press '+' to add a new class value
3. Name of the new class column.
4. Match only at the beginning will begin matching from the beginning of the string. Case sensitive will match by case, too.
5. Produce a report.
6. Press *Apply* to commit the results.

Example
-------

Here is a simple example with the *auto-mpg* dataset. Pass the data to **Create Class**. Select *car_name* as a column to create the new class from. Here, we wish to create new values that match the car brand. First, we type *ford* as the new value for the matching strings. Then we define the substring that will match the data instances. This means that all instances containing *ford* in their *car_name*, will now have a value *ford* in the new class column. Next, we define the same for *honda* and *fiat*. The widget will tell us how many instance are yet unmatched (remaining instances). We will name them *other*, but you can continue creating new values by adding a condition with '+'.

We named our new class column *car_brand* and we matched at the beginning of the string.

![](images/CreateClass-example.png)

Finally, we can observe the new column in a [Data Table](../data/datatable.md) or use the value as color in the [Scatter Plot](../visualize/scatterplot.md).
Feature Statistics
==================

Show basic statistics for data features.

**Inputs**

- Data: input data

**Outputs**
    
- Reduced data: table containing only selected features
- Statistics: table containing statistics of the selected features

The **Feature Statistics** widget provides a quick way to inspect and find interesting features in a given data set.

![](images/feature_statistics-stamped.png)

The Feature Statistics widget on the *heart-disease* data set. The feature *exerc ind ang* was manually changed to a meta variable for illustration purposes.

1. Info on the current data set size and number and types of features
2. The histograms on the right can be colored by any feature. If the selected feature is categorical, a discrete color palette is used (as shown in the example). If the selected feature is numerical, a continuous color palette is used. The table on the right contains statistics about each feature in the data set. The features can be sorted by each statistic, which we now describe.
3. The feature type - can be one of categorical, numeric, time and string.
4. The name of the feature.
5. A histogram of feature values. If the feature is numeric, we appropriately discretize the values into bins. If the feature is categorical, each value is assigned its own bar in the histogram.
6. The central tendency of the feature values. For categorical features, this is the [mode](https://en.wikipedia.org/wiki/Mode_(statistics)). For numeric features, this is [mean](https://en.wikipedia.org/wiki/Mean) value.
7. The dispersion of the feature values. For categorical features, this is the [entropy](https://en.wikipedia.org/wiki/Entropy_(information_theory)) of the value distribution. For numeric features, this is the [coefficient of variation](https://en.wikipedia.org/wiki/Coefficient_of_variation).
8. The minimum value. This is computed for numerical and ordinal categorical features.
9. The maximum value. This is computed for numerical and ordinal categorical features.
10. The number of missing values in the data.

Notice also that some rows are colored differently. White rows indicate regular features, gray rows indicate class variables and the lighter gray indicates meta variables.

Example
-------

The Feature Statistics widget is most often used after the [File](../data/file.md) widget to inspect and find potentially interesting features in the given data set. In the following examples, we use the *heart-disease* data set.

![](images/feature_statistics_workflow.png)

Once we have found a subset of potentially interesting features, or we have found features that we would like to exclude, we can simply select the features we want to keep. The widget outputs a new data set with only these features.

![](images/feature_statistics_example1.png)

Alternatively, if we want to store feature statistics, we can use the *Statistics* output and manipulate those values as needed. In this example, we simply select all the features and display the statistics in a table.

![](images/feature_statistics_example2.png)
Select Columns
==============

Manual selection of data attributes and composition of data domain.

**Inputs**

- Data: input dataset

**Outputs**

- Data: dataset with columns as set in the widget

The **Select Columns** widget is used to manually compose your [data domain](https://en.wikipedia.org/wiki/Data_domain). The user can decide which attributes will be used and how. Orange distinguishes between ordinary attributes, (optional) class attributes and meta attributes. For instance, for building a classification model, the domain would be composed of a set of attributes and a discrete class attribute. Meta attributes are not used in modeling, but several widgets can use them as instance labels.

Orange attributes have a type and are either discrete, continuous or a character string. The attribute type is marked with a symbol appearing before the name of the attribute (D, C, S, respectively).

![](images/SelectColumns-stamped.png)

1. Left-out data attributes that will not be in the output data file
2. Data attributes in the new data file
3. Target variable. If none, the new dataset will be without a target variable.
4. Meta attributes of the new data file. These attributes are included in the dataset but are, for most methods, not considered in the analysis.
5. Produce a report.
6. Reset the domain composition to that of the input data file.
7. Tick if you wish to auto-apply changes of the data domain. 
8. Apply changes of the data domain and send the new data file to the output channel of the widget.

Examples
--------

In the workflow below, the *Iris* data from the [File](../data/file.md) widget is fed into the **Select Columns** widget, where we select to output only two attributes (namely petal width and petal length). We view both the original dataset and the dataset with selected columns in the [Data Table](../data/datatable.md) widget.

![](images/SelectColumns-Example1.png)

For a more complex use of the widget, we composed a workflow to redefine the classification problem in the *heart-disease* dataset. Originally, the task was to predict if the patient has a coronary artery diameter narrowing. We changed the problem to that of gender classification, based on age, chest pain and cholesterol level, and informatively kept the diameter narrowing as a meta attribute.

![](images/SelectColumns-Example2.png)
Pivot Table
===========

Reshape data table based on column values.

**Inputs**

- Data: input data set

**Outputs**

- Pivot Table: contingency matrix as shown in the widget
- Filtered Data: subset selected from the plot
- Grouped Data: aggregates over groups defined by row values

**Pivot Table** summarizes the data of a more extensive table into a table of statistics. The statistics can include sums, averages, counts, etc. The widget also allows selecting a subset from the table and grouping by row values, which have to be a discrete variable. Data with only numeric variables cannot be displayed in the table.

![](images/Pivot-stamped.png)

1. Discrete or numeric variable used for row values. Numeric variables are considered as integers.
2. Discrete variable used for column values. Variable values will appear as columns in the table.
3. Values used for aggregation. Aggregated values will appear as cells in the table.
4. Aggregation methods:
   - For any variable type:
      - *Count*: number of instances with the given row and column value.
      - *Count defined*: number of instances where the aggregation value is defined.
   - For numeric variables:
      - *Sum*: sum of values.
      - *Mean*: average of values.
      - *Mode*: most frequent value of the subset.
      - *Min*: smallest value.
      - *Max*: highest value.
      - *Median*: middle value.
      - *Var*: variance of the subset.
   - For discrete variables:
      - *Majority*: most frequent value of the subset.
5. Tick the box on the left to automatically output any changes. Alternatively, press *Apply* .

Discrete variables
------------------

![](images/Pivot-discrete.png)

Example of a pivot table with only discrete variables selected. We are using *heart-disease* data set for this example. Rows correspond to values of *diameter narrowing* variable. Our columns are values of *gender*, namely female and male. We are using *thal* as values in our cells.

We have selected *Count* and *Majority* as aggregation methods. In the pivot table, we can see the number of instances that do not have diameter narrowing and are female. There are 72 such patients. Concurrently, there are 92 male patients that don't have diameter narrowing. Thal values don't have any effect here, we are just counting occurrences in the data.

The second row shows majority. This means most female patients that don't have diameter narrowing have normal thal results. Conversely, female patients that have diameter narrowing most often have reversable defect.

Numeric variables
-----------------

![](images/Pivot-continuous.png)

Example of a pivot table with numeric variables. We are using *heart-disease* data set for this example. Rows correspond to values of *diameter narrowing* variable. Our columns are values of *gender*, namely female and male. We are using *rest SBP* as values in our cells.

We have selected *Count*, *Sum* and *Median* as aggregation methods. Under *Count*, we see there are 72 female patients that don't have diameter narrowing, same as before for discrete values. What is different are the sum and median aggregations. We see that the sum of resting systolic blood pressure for female patients that don't have diameter narrowing is 9269 and the median value is 130.

Example
-------

We are using *Forest Fires* for this example. The data is loaded in the [Datasets](../data/datasets.md) widget and passed to **Pivot Table**. *Forest Fires* datasets reports forest fires by the month and day they happened. We can aggregate all occurrences of forest fires by selecting *Count* as aggregation method and using *month* as row and *day* as column values. Since we are using *Count*, *Values* variable will have no effect.

We can plot the counts in [Line Plot](../visualize/lineplot.md). But first, let us organize our data a bit. With [Edit Domain](../data/editdomain.md), we will reorder rows values so that months will appear in the correct order, namely from January to December. To do the same for columns, we will use [Select Columns](../data/selectcolumns.md) and reorder day to go from Monday to Sunday.

Finally, our data is ready. Let us pass it to **Line Plot**. We can see that forest fires are most common in August and September, while their frequency is higher during the weekend than during weekdays.

![](images/Pivot-example.png)
Save Data
=========

Saves data to a file.

**Inputs**

- Data: input dataset

The **Save Data** widget considers a dataset provided in the input channel and saves it to a data file with a specified name. It can save the data as:

- a tab-delimited file (.tab)
- comma-separated file (.csv)
- pickle (.pkl), used for storing preprocessing of [Corpus](https://orange.biolab.si/widget-catalog/text-mining/corpus-widget/) objects
- Excel spreadsheets (.xlsx)
- spectra ASCII (.dat)
- hyperspectral map ASCII (.xyz)
- compressed formats (.tab.gz, .csv.gz, .pkl.gz)

The widget does not save the data every time it receives a new signal in the input as this would constantly (and, mostly, inadvertently) overwrite the file. Instead, the data is saved only after a new file name is set or the user pushes the *Save* button.

If the file is saved to the same directory as the workflow or in the subtree of that directory, the widget remembers the relative path. Otherwise, it will store an absolute path but disable auto save for security reasons.

![](images/SaveData.png)

- *Add type annotations to header*: Include Orange's three-row header in the output file.
- *Autosave when receiving new data*: Always save new data. Be careful! This will overwrite existing data on your system.
- *Save* by overwriting the existing file.
- *Save as* to create a new file.

Example
-------

In the workflow below, we used the *Zoo* dataset. We loaded the data into the [Scatter Plot](../visualize/scatterplot.md) widget, with which we selected a subset of data instances and pushed them to the **Save Data** widget to store them in a file.

![](images/Save-Workflow.png)
Aggregate Columns
=================

Compute a sum, max, min ... of selected columns.

**Inputs**

- Data: input dataset

**Outputs**

- Data: extended dataset

**Aggregate Columns** outputs an aggregation of selected columns, for example a sum, min, max, etc.

![](images/AggregateColumns.png)

1. Selected attributes.
2. Operator for aggregation:
   - sum
   - product
   - min
   - max
   - mean
   - variance
   - median
3. Set the name of the computed attribute.
4. If *Apply automatically* is ticked, changes will be communicated automatically. Alternatively, click *Apply*.

Example
-------

We will use iris data from the [File](../data/file.md) widget for this example and connect it to **Aggregate Columns**.

Say we wish to compute a sum of *sepal_length* and *sepal_width* attributes. We select the two attributes from the list.

![](images/AggregateColumns-Example.png)
Preprocess
==========

Preprocesses data with selected methods.

**Inputs**

- Data: input dataset

**Outputs**

- Preprocessor: preprocessing method
- Preprocessed Data: data preprocessed with selected methods

Preprocessing is crucial for achieving better-quality analysis results. The **Preprocess** widget offers several preprocessing methods that can be combined in a single preprocessing pipeline. Some methods are available as separate widgets, which offer advanced techniques and greater parameter tuning.

![](images/preprocess-stamped.png)

1. List of preprocessors. Double click the preprocessors you wish to use and shuffle their order by dragging them up or down. You can also add preprocessors by dragging them from the left menu to the right.
2. Preprocessing pipeline.
3. When the box is ticked (*Send Automatically*), the widget will communicate changes automatically. Alternatively, click *Send*.

Preprocessors
-------------

![](images/Preprocess1.png)

1. List of preprocessors.
2. Discretization of continuous values:
   - [Entropy-MDL discretization](http://sci2s.ugr.es/keel/pdf/algorithm/congreso/fayyad1993.pdf) by Fayyad and Irani that uses [expected information](http://kevinmeurer.com/a-simple-guide-to-entropy-based-discretization/) to determine bins.
   - *Equal frequency discretization* splits by frequency (same number of instances in each bin.
   - *Equal width discretization* creates bins of equal width (span of each bin is the same).
   - *Remove numeric features* altogether.
3. Continuization of discrete values:
   - *Most frequent as base* treats the most frequent discrete value as 0 and others as 1. The discrete attributes with more than 2 values, the most frequent will be considered as a base and contrasted with remaining values in corresponding columns.
   - *One feature per value* creates columns for each value, place 1 where an instance has that value and 0 where it doesn't. Essentially [One Hot Encoding](http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.OneHotEncoder.html).
   - *Remove non-binary features* retains only categorical features that have values of either 0 or 1 and transforms them into continuous.
   - *Remove categorical features* removes categorical features altogether.
   - *Treat as ordinal* takes discrete values and treats them as numbers. If discrete values are categories, each category will be assigned a number as they appear in the data.
   - *Divide by number of values* is similar to treat as ordinal, but the final values will be divided by the total number of values and hence the range of the new continuous variable will be [0, 1].
4. Impute missing values:
   - *Average/Most frequent* replaces missing values (NaN) with the average (for continuous) or most frequent (for discrete) value.
   - *Replace with random value* replaces missing values with random ones within the range of each variable.
   - *Remove rows with missing values*.
5. Select relevant features:
   - Similar to [Rank](../data/rank.md), this preprocessor outputs only the most informative features. Score can be determined by information gain, [gain ratio](https://en.wikipedia.org/wiki/Information_gain_ratio), [gini index](https://en.wikipedia.org/wiki/Gini_coefficient), [ReliefF](https://en.wikipedia.org/wiki/Relief_(feature_selection)), [fast correlation based filter](https://www.aaai.org/Papers/ICML/2003/ICML03-111.pdf), [ANOVA](https://en.wikipedia.org/wiki/One-way_analysis_of_variance), [Chi2](https://en.wikipedia.org/wiki/Chi-squared_distribution), [RReliefF](http://lkm.fri.uni-lj.si/rmarko/papers/robnik03-mlj.pdf), and [Univariate Linear Regression](http://scikit-learn.org/stable/modules/feature_selection.html#feature-selection-using-selectfrommodel).
   - *Strategy* refers to how many variables should be on the output. *Fixed* returns a fixed number of top scored variables, while *Percentile* return the selected top percent of the features.
6. *Select random features* outputs either a fixed number of features from the original data or a percentage. This is mainly used for advanced testing and educational purposes.

![](images/Preprocess2.png)

1. Normalize adjusts values to a common scale. Center values by mean or median or omit centering altogether. Similar for scaling, one can scale by SD (standard deviation), by span or not at all.
2. Randomize instances. Randomize classes shuffles class values and destroys connection between instances and class. Similarly, one can randomize features or meta data. If replicable shuffling is on, randomization results can be shared and repeated with a saved workflow. This is mainly used for advanced testing and educational purposes.
3. *Remove sparse features* retains features that have more than a number/percentage of non-zero/missing values. The rest are discarded.
4. Principal component analysis outputs results of a PCA transformation. Similar to the [PCA](../unsupervised/PCA.md) widget.
5. [CUR matrix decomposition](https://en.wikipedia.org/wiki/CUR_matrix_approximation) is a dimensionality reduction method, similar to SVD.

Preprocessing for predictive modeling
--------------------------------------

When building predictive models, one has to be careful about how to do preprocessing. There are two possible ways to do it in Orange, each slightly different:

1. Connect **Preprocess** to the learner. This will override the default preprocessing pipeline for the learner and apply only custom preprocessing pipeline (default preprocessing steps are described in each learner's documentation).

   ![](images/Preprocess-Models1.png)

2. Connect **Preprocess** to Test and Score. This will apply the preprocessors to each batch within cross-validation. Then the learner's preprocessors will be applied to the preprocessed subset.

   ![](images/Preprocess-Models2.png)

Finally, there's a wrong way to do it. Connecting **Preprocess** directly to the original data and outputting preprocessed data set will likely overfit the model. Don't do it.

   ![](images/Preprocess-Models3.png)

Examples
--------

In the first example, we have used the *heart_disease.tab* dataset available in the dropdown menu of the [File](../data/file.md) widget. then we used **Preprocess** to impute missing values and normalize features. We can observe the changes in the [Data Table](../data/datatable.md) and compare it to the non-processed data.

![](images/Preprocess-Example1.png)

In the second example, we show how to use **Preprocess** for predictive modeling.

This time we are using the *heart_disease.tab* data from the [File](../data/file.md) widget. You can access the data in the dropdown menu. This is a dataset with 303 patients that came to the doctor suffering from a chest pain. After the tests were done, some patients were found to have diameter narrowing and others did not (this is our class variable).

Some values are missing in our data set, so we would like to impute missing values before evaluating the model. We do this by passing a preprocessor directly to [Test and Score](../evaluate/testandscore.md). In **Preprocess**, we set the correct preprocessing pipeline (in our example only a single preprocessor with *Impute missing values*), then connect it to the Preprocessor input of Test and Score.

We also pass the data and the learner (in this case, a [Logistic Regression](../model/logisticregression.md)). This is the correct way to pass a preprocessor to cross-validation as each fold will independently get preprocessed in the training phase. This is particularly important for feature selection.

![](images/Preprocess-Example2.png)
Group by
========

Groups data by selected variables and aggregate columns with selected aggregations.

**Inputs**

- Data: input data table

**Outputs**

- Data: aggregated data

Group By widget first identifies groups based on selected variables in the **Group by** list. Groups are defined by all distinct combinations of values in selected variables.

In the second step, the widget computes aggregations defined in the table on the right side of the widget for each group.


![](images/Group-by-stamped.png)

1. Select variables that define groups
2. View variables and their aggregations. To change aggregation for one or more variables, select them in the table.
3. Change aggregations for variables selected in the view above.
4. When the *Send automatically* box is ticked, all changes will be automatically communicated to other widgets.
5. Get documentation, observe a number of items on input or output

Examples
--------

We first load **heart_disease** dataset in the **File** widget. In the **Group By** widget, we set variables that define groups -- **diameter narrowing** and **gender**. Each group includes items (rows) that belong to one combination of both variables. 

In the table on the right-hand side of the widget, we set that we want to compute **mean** and  **median** for values of **rest SBP** variable in each group, **median** for values of **cholesterol** variable, and **mean** for **major vessels colored**.

In the **Data Table** widget, we can see that both females and males have lower average values for **rest SBP** when **diameter narrowing** is 0. The difference is greater for females. The median of **rest SBP** is different only for females, while for males is the same.

You can also observe differences between median **cholesterol** level and mean value of **major vessel colored** between groups.


![](images/Group-by-example.png)
Transpose
=========

Transposes a data table.

**Inputs**

- Data: input dataset

**Outputs**

- Data: transposed dataset

**Transpose** widget transposes data table.

![](images/transpose-stamped.png)

Example
-------

This is a simple workflow showing how to use **Transpose**. Connect the widget to [File](../data/file.md) widget. The output of **Transpose** is a transposed data table with rows as columns and columns as rows. You can observe the result in a [Data Table](../data/datatable.md).

![](images/transpose-example.png)
Color
=====

Set color legend for variables.

**Inputs**

- Data: input data set

**Outputs**

- Data: data set with a new color legend

The **Color** widget sets the color legend for visualizations.

![](images/Color-stamped.png)

1. A list of discrete variables. Set the color of each variable by double-clicking on it. The widget also enables renaming variables by clicking on their names.
2. A list of continuous variables. Click on the color strip to choose a different palette. To use the same palette for all variables, change it for one variable and click *Copy to all* that appears on the right. The widget also enables renaming variables by clicking on their names.
3. Produce a report.
4. Apply changes. If *Apply automatically* is ticked, changes will be communicated automatically. Alternatively, just click *Apply*.

![](images/Color-Continuous_unindexed.png)

Palettes for numeric variables are grouped and tagged by their properties.

- Diverging palettes have two colors on its ends and a central color (white or black) in the middle. Such palettes are particularly useful when the the values can be positive or negative, as some widgets (for instance the Heat map) will put the 0 at the middle point in the palette.

- Linear palettes are constructed so that human perception of the color change is linear with the change of the value.

- Color blind palettes cover different types of color blindness, and can also be linear or diverging.

- In isoluminant palettes, all colors have equal brightness.

- Rainbow palettes are particularly nice in widgets that bin numeric values in visualizations.

Example
-------

We chose to work with the *heart_disease* data set. We opened the color palette and selected two new colors for diameter narrowing variable. Then we opened the [Scatter Plot](../visualize/scatterplot.md) widget and viewed the changes made to the scatter plot.

![](images/Color-Example-Discrete.png)

To see the effect of color palettes for numeric variables, we color the points in the scatter plot by cholesterol and change the palette for this attribute in the Color widget.

![](images/Color-Example-Continuous.png)
Continuize
==========

Turns discrete variables (attributes) into numeric ("continuous") dummy variables.

**Inputs**

- Data: input data set

**Outputs**

- Data: transformed data set

The **Continuize** widget receives a data set in the input and outputs the same data set in which the discrete variables (including binary variables) are replaced with continuous ones.

![](images/Continuize-stamped.png)

1. Define the treatment of non-binary categorical variables.

    Examples in this section will assume that we have a discrete attribute status with the values low, middle and high, listed in that order. Options for their transformation are:

   - **First value as base**: a N-valued categorical variable will be transformed into N-1 numeric variables, each serving as an indicator for one of the original values except for the base value. The base value is the first value in the list. By default, the values are ordered alphabetically; their order can be changed in [Edit Domain](../data/editdomain).

       In the above case, the three-valued variable *status* is transformed into two numeric variables, *status=middle* with values 0 or 1 indicating whether the original variable had value *middle* on a particular example, and similarly, *status=high*.

   - **Most frequent value as base**: similar to the above, except that the most frequent value is used as a base. So, if the most frequent value in the above example is *middle*, then *middle* is considered as the base and the two newly constructed variables are *status=low* and *status=high*.

   - **One attribute per value**: this option constructs one numeric variable per each value of the original variable. In the above case, we would get variables *status=low*, *status=middle* and *status=high*.

   - **Ignore multinomial attributes**: removes non-binary categorical variables from the data.

   - **Treat as ordinal**: converts the variable into a single numeric variable enumerating the original values. In the above case, the new variable would have the value of 0 for *low*, 1 for *middle* and 2 for *high*. Again note that the order of values can be set in  [Edit Domain](../data/editdomain).

   - **Divide by number of values**: same as above, except that values are normalized into range 0-1. In our example, the values of the new variable would be 0, 0.5 and 1.

2. Define the treatment of continuous attributes. Besised the option to *Leave them as they are*, we can *Normalize by span*, which will subtract the lowest value found in the data and divide by the span, so all values will fit into [0, 1]. Option *Normalize by standard deviation* subtracts the average and divides by the standard deviation.

3. Define the treatment of class attributes (outcomes, targets). Besides leaving it as it is, the available options mirror those for multinomial attributes, except for those that would split the outcome into multiple outcome variables.

4. This option defines the ranges of new variables. In the above text, we supposed the range *from 0 to 1*.

5. Produce a report.

6. If *Apply automatically* is ticked, changes are committed automatically. Otherwise, you have to press *Apply* after each change.

Examples
--------

First, let's see what is the output of the **Continuize** widget. We feed the original data (the *Heart disease* data set) into the [Data Table](../data/datatable) and see how they look like. Then we continuize the discrete values and observe them in another [Data Table](../data/datatable).

![](images/Continuize-Example1.png)

In the second example, we show a typical use of this widget - in order to properly plot the linear projection of the data, discrete attributes need to be converted to continuous ones and that is why we put the data through the **Continuize** widget before drawing it. The attribute "*chest pain*" originally had four values and was transformed into three continuous attributes; similar happened to gender, which was transformed into a single attribute "*gender=female*".

![](images/Continuize-Example2.png)
Python Script
=============

Extends functionalities through Python scripting.

**Inputs**

- Data (Orange.data.Table): input dataset bound to ``in_data`` variable
- Learner (Orange.classification.Learner): input learner bound to ``in_learner`` variable
- Classifier (Orange.classification.Learner): input classifier bound to ``in_classifier`` variable
- Object: input Python object bound to ``in_object`` variable

**Outputs**

- Data (Orange.data.Table): dataset retrieved from ``out_data`` variable
- Learner (Orange.classification.Learner): learner retrieved from ``out_learner`` variable
- Classifier (Orange.classification.Learner): classifier retrieved from ``out_classifier`` variable
- Object: Python object retrieved from ``out_object`` variable

**Python Script** widget can be used to run a python script in the input, when a suitable functionality is not implemented in an existing widget. The script has ``in_data``, ``in_distance``, ``in_learner``, ``in_classifier`` and ``in_object`` variables (from input signals) in its local namespace. If a signal is not connected or it did not yet receive any data, those variables contain ``None``.

After the script is executed variables from the script’s local namespace are extracted and used as outputs of the widget. The widget can be further connected to other widgets for visualizing the output.

For instance the following script would simply pass on all signals it receives:

    out_data = in_data
    out_distance = in_distance
    out_learner = in_learner
    out_classifier = in_classifier
    out_object = in_object

Note: You should not modify the input objects in place.

![](images/PythonScript-stamped.png)

1. Info box contains names of basic operators for Orange Python script.
2. The *Library* control can be used to manage multiple scripts. Pressing "+" will add a new entry and open it in the *Python script* editor. When the script is modified, its entry in the *Library* will change to indicate it has unsaved changes. Pressing *Update* will save the script (keyboard shortcut "Ctrl+S"). A script can be removed by selecting it and pressing the "-" button.
3. Pressing *Execute* in the *Run* box executes the script (keyboard shortcut "Ctrl+R"). Any script output (from ``print``) is captured and displayed in the *Console* below the script.
4. The *Python script* editor on the left can be used to edit a script (it supports some rudimentary syntax highlighting).
5. Console displays the output of the script.

Examples
--------

Python Script widget is intended to extend functionalities for advanced users. Classes from Orange library are described in the [documentation](https://docs.biolab.si/3/data-mining-library/#reference). To find further information about orange Table class see [Table](https://docs.biolab.si/3/data-mining-library/reference/data.table.html), [Domain](https://docs.biolab.si/3/data-mining-library/reference/data.domain.html), and [Variable](https://docs.biolab.si/3/data-mining-library/reference/data.variable.html) documentation.

One can, for example, do batch filtering by attributes. We used zoo.tab for the example and we filtered out all the attributes that have more than 5 discrete values. This in our case removed only 'leg' attribute, but imagine an example where one would have many such attributes.

    from Orange.data import Domain, Table
    domain = Domain([attr for attr in in_data.domain.attributes
                     if attr.is_continuous or len(attr.values) <= 5],
                    in_data.domain.class_vars)
    out_data = Table(domain, in_data)

![](images/PythonScript-filtering.png)

The second example shows how to round all the values in a few lines of code. This time we used wine.tab and rounded all the values to whole numbers.

    import numpy as np
    out_data = in_data.copy()
    #copy, otherwise input data will be overwritten
    np.round(out_data.X, 0, out_data.X)

![](images/PythonScript-round.png)

The third example introduces some Gaussian noise to the data. Again we make a copy of the input data, then walk through all the values with a double for loop and add random noise.

    import random
    from Orange.data import Domain, Table
    new_data = in_data.copy()
    for inst in new_data:
      for f in inst.domain.attributes:
        inst[f] += random.gauss(0, 0.02)
    out_data = new_data

![](images/PythonScript-gauss.png)

The final example uses Orange3-Text add-on. **Python Script** is very useful for custom preprocessing in text mining, extracting new features from strings, or utilizing advanced *nltk* or *gensim* functions. Below, we simply tokenized our input data from *deerwester.tab* by splitting them by whitespace.

    print('Running Preprocessing ...')
    tokens = [doc.split(' ') for doc in in_data.documents]
    print('Tokens:', tokens)
    out_object = in_data
    out_object.store_tokens(tokens)

You can add a lot of other preprocessing steps to further adjust the output. The output of **Python Script** can be used with any widget that accepts the type of output your script produces. In this case, connection is green, which signalizes the right type of input for Word Cloud widget.

![](images/PythonScript-Example3.png)
Select Rows
===========

Selects data instances based on conditions over data features.

**Inputs**

- Data: input dataset

**Outputs**

- Matching Data: instances that match the conditions
- Non-Matching Data: instances that do not match the conditions
- Data: data with an additional column showing whether a instance is selected

This widget selects a subset from an input dataset, based on user-defined conditions. Instances that match the selection rule are placed in the output *Matching Data* channel.

Criteria for data selection are presented as a collection of conjunct terms (i.e. selected items are those matching all the terms in '*Conditions*').

Condition terms are defined through selecting an attribute, selecting an operator from a list of operators, and, if needed, defining the value to be used in the condition term. Operators are different for discrete, continuous and string attributes.

![](images/SelectRows-stamped.png)

1. Conditions you want to apply, their operators and related values
2. Add a new condition to the list of conditions.
3. Add all the possible variables at once.
4. Remove all the listed variables at once.
5. Information on the input dataset and information on instances that match the condition(s)
6. Purge the output data.
7. When the *Send automatically* box is ticked, all changes will be automatically communicated to other widgets.
8. Produce a report.

Any change in the composition of the condition will update the information pane (*Data Out*).

If *Send automatically* is selected, then the output is updated on any change in the composition of the condition or any of its terms.

Example
-------

In the workflow below, we used the *Zoo* data from the [File](../data/file.md) widget and fed it into the **Select Rows** widget. In the widget, we chose to output only two animal types, namely fish and reptiles. We can inspect both the original dataset and the dataset with selected rows in the [Data Table](../data/datatable.md) widget.

![](images/SelectRows-Example.png)

In the next example, we used the data from the *Titanic* dataset and similarly fed it into the [Box Plot](../visualize/boxplot.md) widget. We first observed the entire dataset based on survival. Then we selected only first class passengers in the **Select Rows** widget and fed it again into the [Box Plot](../visualize/boxplot.md). There we could see all the first class passengers listed by their survival rate and grouped by gender.

![](images/SelectRows-Workflow.png)
Create Instance
===============

Interactively creates an instance from a sample dataset.

**Inputs**

- Data: input dataset
- Reference: refrence dataset

**Outputs**

- Data: input dataset appended the created instance

The **Create Instance** widget creates a new instance, based on the input data. The widget displays all variables of the input dataset in a table of two columns. The column *Variable* represents the variable's name, meanwhile the column *Value* enables setting the variable's value. Each value is initially set to median value of the variable. The values can be manually set to *Median*, *Mean*, *Random* or *Input* by clicking the corresponding button. For easier searching through the variables, the table has filter attached. When clicking upon one of the mentioned buttons, only filtered variables are considered. One can also set the value by right-clicking a row and selecting an option in a context menu.

![](images/CreateInstance-stamped.png)

1. Filter table by variable name.
2. The column represents a variable's name and type. The table can be sorted by clicking the columns header. 
3. Provides controls for value editing.
4. Set filtered variables' values to:
   - *Median*: median value of variable in the input dataset
   - *Mean*: mean value of variable in the input dataset
   - *Random*: random value in a range of variable in the input dataset
   - *Input*: median value of variable in the reference dataset
5. If *Append this instance to input data* is ticked, the created instance is appended to the input dataset. Otherwise, a single instance appears on the output. To distinguish between created and original data, *Source ID* variable is added.
5. If *Apply automatically* is ticked, changes are committed automatically. Otherwise, you have to press *Apply* after each change.
6. Produce a report.
7. Information on input and reference dataset.
8. Information on output dataset.

Example
-------

The **Create Instance** is usually used to examine a model performance on some arbitrary data. The basic usage is shown in the following workflow, where a (*Housing*) dataset is used to fit a [Linear Regression](../model/linearregression.md) model, which is than used to [predict](../evaluate/predictions.md) a target value for data, created by the *Create Instance* widget. Inserting a [Rank](../data/rank.md) widget between [File](../data/file.md) and *Create Instance* enables outputting (and therefore making predictions on) the most important features. 
A [Select Column](../data/selectcolumns.md) widget is inserted to omit the actual target value.

![](images/CreateInstance-example.png)

The next example shows how to check whether the created instance is some kind of outlier. The creates instance is feed to [PCA](../unsupervised/PCA.md) whose first and second componens are then examined in a [Scatter Plot](../visualize/scatterplot.md). The created instance is colored red in the plot and it could be considered as an outlier if it appears far from the original data (blue).

![](images/CreateInstance-example2.png)

Discretize
==========

Discretizes continuous attributes from an input dataset.

**Inputs**

- Data: input dataset

**Outputs**

- Data: dataset with discretized values

The **Discretize** widget [discretizes](https://en.wikipedia.org/wiki/Discretization) continuous attributes with a selected method.

![](images/Discretize-All-stamped.png)

1. The basic version of the widget is rather simple. It allows choosing between three different discretizations.
   - [Entropy-MDL](http://ijcai.org/Past%20Proceedings/IJCAI-93-VOL2/PDF/022.pdf), invented by Fayyad and Irani is a top-down discretization, which recursively splits the attribute at a cut maximizing information gain, until the gain is lower than the minimal description length of the cut. This discretization can result in an arbitrary number of intervals, including a single interval, in which case the attribute is discarded as useless (removed).
   - [Equal-frequency](http://www.saedsayad.com/unsupervised_binning.htm) splits the attribute into a given number of intervals, so that they each contain approximately the same number of instances.
   - [Equal-width](https://en.wikipedia.org/wiki/Data_binning) evenly splits the range between the smallest and the largest observed value. The *Number of intervals* can be set manually.
   - The widget can also be set to leave the attributes continuous or to remove them.
2. To treat attributes individually, go to **Individual Attribute Settings**. They show a specific discretization of each attribute and allow changes. First, the top left list shows the cut-off points for each attribute. In the snapshot, we used the entropy-MDL discretization, which determines the optimal number of intervals automatically; we can see it discretized the age into seven intervals with cut-offs at 21.50, 23.50, 27.50, 35.50, 43.50, 54.50 and 61.50, respectively, while the capital-gain got split into many intervals with several cut-offs. The final weight (fnlwgt), for instance, was left with a single interval and thus removed.
On the right, we can select a specific discretization method for each attribute. Attribute *“fnlwgt”* would be removed by the MDL-based discretization, so to prevent its removal, we select the attribute and choose, for instance, **Equal-frequency discretization**. We could also choose to leave the attribute continuous.
3. Produce a report.
4. Tick *Apply automatically* for the widget to automatically commit changes. Alternatively, press *Apply*.

Example
-------

In the schema below, we show the *Iris* dataset with continuous attributes
(as in the original data file) and with discretized attributes.

![](images/Discretize-Example.png)
Paint Data
==========

Paints data on a 2D plane. You can place individual data points or use a brush to paint larger datasets.

**Outputs**

- Data: dataset as painted in the plot

The widget supports the creation of a new dataset by visually placing data points on a two-dimension plane. Data points can be placed on the plane individually (*Put*) or in a larger number by brushing (*Brush*). Data points can belong to classes if the data is intended to be used in supervised learning.

![](images/PaintData-stamped.png)

1. Name the axes and select a class to paint data instances. You can add or remove classes. Use only one class to create classless, unsupervised datasets.
2. Drawing tools. Paint data points with *Brush* (multiple data instances) or *Put* (individual data instance). Select data points with *Select* and remove them with the Delete/Backspace key. Reposition data points with [Jitter](https://en.wikipedia.org/wiki/Jitter) (spread) and *Magnet* (focus). Use *Zoom* and scroll to zoom in or out. Below, set the radius and intensity for Brush, Put, Jitter and Magnet tools.
3. Reset to Input Data.
4. *Save Image* saves the image to your computer in a .svg or .png format.
5. Produce a report.
6. Tick the box on the left to automatically commit changes to other widgets. Alternatively, press *Send* to apply them.

Example
-------

In the example below, we have painted a dataset with 4 classes. Such dataset is great for demonstrating k-means and hierarchical clustering methods. In the screenshot, we see that [k-Means](../unsupervised/kmeans.md), overall, recognizes clusters better than [Hierarchical Clustering](../unsupervised/hierarchicalclustering.md). It returns a score rank, where the best score (the one with the highest value) means the most likely number of clusters. Hierarchical clustering, however, doesn’t group the right classes together. This is a great tool for learning and exploring statistical concepts.

![](images/PaintData-Example.png)
Datasets
========

Load a dataset from an online repository.

**Outputs**

- Data: output dataset

**Datasets** widget retrieves selected dataset from the server and sends it to the output. File is downloaded to the local memory and thus instantly available even without the internet connection. Each dataset is provided with a description and information on the data size, number of instances, number of variables, target and tags.

![](images/Datasets-stamped.png)

1. Information on the number of datasets available and the number of them downloaded to the local memory.
2. Content of available datasets. Each dataset is described with the size, number of instances and variables, type of the target variable and tags.
3. Formal description of the selected dataset.
4. If *Send Data Automatically* is ticked, selected dataset is communicated automatically. Alternatively, press *Send Data*.

Example
-------

Orange workflows can start with **Datasets** widget instead of **File** widget. In the example below, the widget retrieves a dataset from an online repository (Kickstarter data), which is subsequently sent to both the [Data Table](../data/datatable) and the [Distributions](../visualize/distributions).

![](images/Datasets-Workflow.png)
Randomize
=========

Shuffles classes, attributes and/or metas of an input dataset.

**Inputs**

- Data: input dataset

**Outputs**

- Data: randomized dataset

The **Randomize** widget receives a dataset in the input and outputs the same dataset in which the classes, attributes or/and metas are shuffled.

![](images/Randomize-Default.png)

1. Select group of columns of the dataset you want to shuffle.
2. Select proportion of the dataset you want to shuffle.
3. Produce replicable output.
4. If *Apply automatically* is ticked, changes are committed automatically. Otherwise, you have to press *Apply* after each change.
5. Produce a report.

Example
-------

The **Randomize** widget is usually placed right after (e.g. [File](../data/file.md) widget. The basic usage is shown in the following workflow, where values of class variable of Iris dataset are randomly shuffled.

![](images/Randomize-Example1.png)

In the next example we show how shuffling class values influences model performance on the same dataset as above.

![](images/Randomize-Example2.png)
Data Sampler
============

Selects a subset of data instances from an input dataset.

**Inputs**

- Data: input dataset

**Outputs**

- Data Sample: sampled data instances
- Remaining Data: out-of-sample data

The **Data Sampler** widget implements several data sampling methods. It outputs a sampled and a complementary dataset (with instances from the input set that are not included in the sampled dataset). The output is processed after the input dataset is provided and *Sample Data* is pressed.

![](images/DataSampler-stamped.png)

1. Information on the input and output dataset.
2. The desired sampling method:
   - **Fixed proportion of data** returns a selected percentage of the entire data (e.g. 70% of all the data)
   - **Fixed sample size** returns a selected number of data instances with a chance to set *Sample with replacement*, which always samples from the entire dataset (does not subtract instances already in the subset). With replacement, you can generate more instances than available in the input dataset.
   - [Cross Validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics)) partitions data instances into the specified number of complementary subsets. Following a typical validation schema, all subsets except the one selected by the user are output as Data Sample, and the selected subset goes to Remaining Data. (Note: In older versions, the outputs were swapped. If the widget is loaded from an older workflow, it switches to compatibility mode.)
   - [Bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)) infers the sample from the population statistic.
3. *Replicable sampling* maintains sampling patterns that can be carried
   across users, while *stratify sample* mimics the composition of the
   input dataset.
4. Press *Sample Data* to output the data sample.

If all data instances are selected (by setting the proportion to 100 % or setting the fixed sample size to the entire data size), output instances are still shuffled.

Examples
--------

First, let's see how the **Data Sampler** works. We will use the *iris* data from the [File](../data/file.md) widget. We see there are 150 instances in the data. We sampled the data with the **Data Sampler** widget and we chose to go with a fixed sample size of 5 instances for simplicity. We can observe the sampled data in the [Data Table](../data/datatable.md) widget (Data Table (in-sample)). The second [Data Table](../data/datatable.md) (Data Table (out-of-sample)) shows the remaining 145 instances that weren't in the sample. To output the out-of-sample data, double-click the connection between the widgets and rewire the output to *Remaining Data --> Data*.

![](images/DataSampler-Example1.png)

Now, we will use the **Data Sampler** to split the data into training and testing part. We are using the *iris* data, which we loaded with the [File](../data/file.md) widget. In **Data Sampler**, we split the data with *Fixed proportion of data*, keeping 70% of data instances in the sample.

Then we connected two outputs to the [Test & Score](../evaluate/testandscore.md) widget, *Data Sample --> Data* and *Remaining Data --> Test Data*. Finally, we added [Logistic Regression](../model/logisticregression.md) as the learner. This runs logistic regression on the Data input and evaluates the results on the Test Data.

![](images/DataSampler-Example2.png)

Over/Undersampling
------------------

**Data Sampler** can also be used to oversample a minority class or undersample majority class in the data. Let us show an example for oversampling. First, separate the minority class using a [Select Rows](../data/selectrows.md) widget. We are using the *iris* data from the [File](../data/file.md) widget. The data set has 150 data instances, 50 of each class. Let us oversample, say, *iris-setosa*.

In **Select Rows**, set the condition to *iris is iris-setosa*. This will output 50 instances of the *iris-setosa* class. Now, connect *Matching Data* into the **Data Sampler**, select *Fixed sample size*, set it to, say, 100 and select *Sample with replacement*. Upon pressing *Sample Data*, the widget will output 100 instances of *iris-setosa* class, some of which will be duplicated (because we used *Sample with replacement*).

Finally, use [Concatenate](../data/concatenate) to join the oversampled instances and the *Unmatched Data* output of the **Select Rows** widget. This outputs a data set with 200 instances. We can observe the final results in the [Distributions](../visualize/distributions).

![](images/DataSampler-Example-OverUnderSampling.png)
Rank
====

Ranking of attributes in classification or regression datasets.

**Inputs**

- Data: input dataset
- Scorer: models for feature scoring

**Outputs**

- Reduced Data: dataset with selected attributes
- Scores: data table with feature scores
- Features: list of attributes

The **Rank** widget scores variables according to their correlation with discrete or numeric target variable, based on applicable internal scorers (like information gain, chi-square and linear regression) and any connected external models that supports scoring, such as linear regression, logistic regression, random forest, SGD, etc. The widget can also handle unsupervised data, but only by external scorers, such as PCA.

![](images/Rank-stamped.png)

1. Select scoring methods. See the options for classification, regression and unsupervised data in the **Scoring methods** section.
2. Select attributes to output. *None* won't output any attributes, while *All* will output all of them. With manual selection, select the attributes from the table on the right. *Best ranked* will output n best ranked attributes.
   If *Send Automatically* is ticked, the widget automatically communicates changes to other widgets.
3. Status bar. Produce a report by clicking on the file icon. Observe input and output of the widget. On the right, warnings and errors are shown.

Scoring methods (classification)
--------------------------------

1. Information Gain: the expected amount of information (reduction of entropy)
2. [Gain Ratio](https://en.wikipedia.org/wiki/Information_gain_ratio): a ratio of the information gain and the attribute's intrinsic information, which reduces the bias towards multivalued features that occurs in information gain
3. [Gini](https://en.wikipedia.org/wiki/Gini_coefficient): the inequality among values of a frequency distribution
4. [ANOVA](https://en.wikipedia.org/wiki/One-way_analysis_of_variance): the difference between average values of the feature in different classes
5. [Chi2](https://en.wikipedia.org/wiki/Chi-squared_distribution): dependence between the feature and the class as measured by the chi-square statistic
6. [ReliefF](https://en.wikipedia.org/wiki/Relief_(feature_selection)): the ability of an attribute to distinguish between classes on similar data instances
7. [FCBF (Fast Correlation Based Filter)](https://www.aaai.org/Papers/ICML/2003/ICML03-111.pdf): entropy-based measure, which also identifies redundancy due to pairwise correlations between features

Additionally, you can connect certain learners that enable scoring the features according to how important they are in models that the learners build (e.g. [Logistic Regression](../model/logisticregression.md), [Random Forest](../model/randomforest.md), [SGD](../model/stochasticgradient.md)). Please note that the data is normalized before ranking.

Scoring methods (regression)
----------------------------

1. [Univariate Regression](https://en.wikipedia.org/wiki/Simple_linear_regression): linear regression for a single variable
2. [RReliefF](http://www.clopinet.com/isabelle/Projects/reading/robnik97-icml.pdf): relative distance between the predicted (class) values of the two instances.

Additionally, you can connect regression learners (e.g. [Linear Regression](../model/linearregression.md), [Random Forest](../model/randomforest.md), [SGD](../model/stochasticgradient.md)). Please note that the data is normalized before ranking.

Scoring method (unsupervised)
-----------------------------

Currently, only [PCA](../unsupervised/PCA.md) is supported for unsupervised data. Connect PCA to Rank to obtain the scores. The scores correspond to the correlation of a variable with the individual principal component.

Scoring with learners
---------------------

Rank can also use certain learners for feature scoring. See [Learners as Scorers](../../learners-as-scorers/index.md) for an example.

Example: Attribute Ranking and Selection
----------------------------------------

Below, we have used the **Rank** widget immediately after the [File](../data/file.md) widget to reduce the set of data attributes and include only the most informative ones:

![](images/Rank-Select-Schema.png)

Notice how the widget outputs a dataset that includes only the best-scored attributes:

![](images/Rank-Select-Widgets.png)

Example: Feature Subset Selection for Machine Learning
------------------------------------------------------

What follows is a bit more complicated example. In the workflow below, we first split the data into a training set and a test set. In the upper branch, the training data passes through the **Rank** widget to select the most informative attributes, while in the lower branch there is no feature selection. Both feature selected and original datasets are passed to their own [Test & Score](../evaluate/testandscore.md) widgets, which develop a *Naive Bayes* classifier and score it on a test set.

![](images/Rank-and-Test.png)

For datasets with many features, a naive Bayesian classifier feature selection, as shown above, would often yield a better predictive accuracy.
Neighbors
=========

Compute nearest neighbors in data according to reference.

**Inputs**

- Data: An input data set.
- Reference: A reference data for neighbor computation. 

**Outputs**

- Neighbors: A data table of nearest neighbors according to reference.

The **Neighbors** widget computes nearest neighbors for a given reference and for a given distance measure. The reference can be either one instance or more instances. In the case with one reference widget outputs closest `n` instances from data where `n` is set by the **Number of neighbors** option in the widget. When reference contains more instances widget computes the combined distance for each data instance as a minimum of distances to each reference. Widget outputs `n` data instances with lowest combined distance.

![](images/neighbours-stamped.png)

1. Distance measure for computing neighbors. Supported measures are: Euclidean, Manhattan, Mahalanobis, Cosine, Jaccard, Spearman, absolute Spearman, Pearson, absolute Pearson. 
2. Number of neighbors on the output.
3. If *Exclude rows (equal to) references* is ticked, data instances that are highly similar to the reference (distance < 1e-5), will be excluded.
4. Click *Apply* to commit the changes. To communicate changes automatically tick *Apply Automatically*.
5. Status bar with access to widget help and information on the input and output data.

Examples
--------

In the first example, we used *iris* data and passed it to **Neighbors** and to [Data Table](../data/datatable.md). In **Data Table**, we selected an instance of iris, that will serve as our reference, meaning we wish to retrieve 10 closest examples to the select data instance. We connect **Data Table** to **Neighbors** as well.

We can observe the results of neighbor computation in **Data Table (1)**, where we can see 10 closest images to our selected iris flower.

![](images/neighbours-example1.png)

Now change the selection **Data Table** to multiple examples. As a result, we get instances with closest combined distances to the references. The method computes the combined distance as a minimum of distances to each reference.

![](images/neighbours-example-multiple.png)

Another example requires the installation of Image Analytics add-on. We loaded 15 paintings from famous painters with **Import Images** widget and passed them to **Image Embedding**, where we selected *Painters* embedder.

Then the procedure is the same as above. We passed embedded images to **Image Viewer** and selected a painting from Monet to serve as our reference image. We passed the image to **Neighbors**, where we set the distance measure to *cosine*, ticked off *Exclude reference* and set the neighbors to 2. This allows us to find the actual closest neighbor to a reference painting and observe them side by side in **Image Viewer (1)**.

![](images/neighbours-example2.png)
Impute
======

Replaces unknown values in the data.

**Inputs**

- Data: input dataset
- Learner: learning algorithm for imputation

**Outputs**

- Data: dataset with imputed values

Some Orange's algorithms and visualizations cannot handle unknown values in the data. This widget does what statisticians call imputation: it substitutes missing values by values either computed from the data or set by the user. The default imputation is (1-NN).

![](images/impute-stamped.png)

1. In the top-most box, *Default method*, the user can specify a general imputation technique for all attributes.
   - **Don't Impute** does nothing with the missing values.
   - **Average/Most-frequent** uses the average value (for continuous attributes) or the most common value (for discrete attributes).
   - **As a distinct value** creates new values to substitute the missing ones.
   - **Model-based imputer** constructs a model for predicting the missing value, based on values of other attributes; a separate model is constructed for each attribute. The default model is 1-NN learner, which takes the value from the most similar example (this is sometimes referred to as hot deck imputation). This algorithm can be substituted by one that the user connects to the input signal Learner for Imputation. Note, however, that if there are discrete and continuous attributes in the data, the algorithm needs to be capable of handling them both; at the moment only 1-NN learner can do that. (In the future, when Orange has more regressors, the Impute widget may have separate input signals for discrete and continuous models.)
   - **Random values** computes the distributions of values for each attribute and then imputes by picking random values from them.
   - **Remove examples with missing values** removes the example containing missing values. This check also applies to the class attribute if *Impute class values* is checked.

2. It is possible to specify individual treatment for each attribute, which overrides the default treatment set. One can also specify a manually defined value used for imputation. In the screenshot, we decided not to impute the values of "*normalized-losses*" and "*make*", the missing values of "*aspiration*" will be replaced by random values, while the missing values of "*body-style*" and "*drive-wheels*" are replaced by "*hatchback*" and "*fwd*",respectively. If the values of "*length*", "*width*" or "*height*" are missing, the example is discarded. Values of all other attributes use the default method set above (model-based imputer, in our case).
3. The imputation methods for individual attributes are the same as default methods.
4. *Restore All to Default* resets the individual attribute treatments to default.
5. Produce a report.
6. All changes are committed immediately if *Apply automatically* is checked. Otherwise, *Apply* needs to be ticked to apply any new settings.

Example
-------

To demonstrate how the **Impute** widget works, we played around with the *Iris* dataset and deleted some of the data. We used the **Impute** widget and selected the *Model-based imputer* to impute the missing values. In another [Data Table](../data/datatable.md), we see how the question marks turned into distinct values ("Iris-setosa, "Iris-versicolor").

![](images/Impute-Example.png)
AdaBoost
========

An ensemble meta-algorithm that combines weak learners and adapts to the 'hardness' of each training sample.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)
- Learner: learning algorithm

**Outputs**

- Learner: AdaBoost learning algorithm
- Model: trained model

The [AdaBoost](https://en.wikipedia.org/wiki/AdaBoost) (short for "Adaptive boosting") widget is a machine-learning algorithm, formulated by [Yoav Freund and Robert Schapire](https://cseweb.ucsd.edu/~yfreund/papers/IntroToBoosting.pdf). It can be used with other learning algorithms to boost their performance. It does so by tweaking the weak learners.

**AdaBoost** works for both classification and regression.

![](images/AdaBoost-stamped.png)

1. The learner can be given a name under which it will appear in other widgets. The default name is "AdaBoost".
2. Set the parameters. The base estimator is a tree and you can set:
   - *Number of estimators*
   - *Learning rate*: it determines to what extent the newly acquired information will override the old information (0 = the agent will not learn anything, 1 = the agent considers only the most recent information)
   - *Fixed seed for random generator*: set a fixed seed to enable reproducing the results.
3. Boosting method.
   - *Classification algorithm* (if classification on input): SAMME (updates base estimator's weights with classification results) or SAMME.R (updates base estimator's weight with probability estimates).
   - *Regression loss function* (if regression on input): Linear (), Square (), Exponential ().
4. Produce a report.
5. Click *Apply* after changing the settings. That will put the new learner in the output and, if the training examples are given, construct a new model and output it as well. To communicate changes automatically tick *Apply Automatically*.

Preprocessing
-------------

AdaBoost uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- continuizes categorical variables (with one-hot-encoding)
- removes empty columns
- imputes missing values with mean values

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Examples
--------

For classification, we loaded the *iris* dataset. We used *AdaBoost*, [Tree](../model/tree.md) and [Logistic Regression](../model/logisticregression.md) and evaluated the models' performance in [Test & Score](../evaluate/testandscore.md).

![](images/AdaBoost-classification.png)

For regression, we loaded the *housing* dataset, sent the data instances to two different models (**AdaBoost** and [Tree](../model/tree.md)) and output them to the [Predictions](../evaluate/predictions.md) widget.

![](images/AdaBoost-regression.png)
Random Forest
=============

Predict using an ensemble of decision trees.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: random forest learning algorithm
- Model: trained model

[Random forest](https://en.wikipedia.org/wiki/Random_forest) is an ensemble learning method used for classification, regression and other tasks. It was first proposed by Tin Kam Ho and further developed by Leo Breiman (Breiman, 2001) and Adele Cutler.

**Random Forest** builds a set of decision trees. Each tree is developed from a bootstrap sample from the training data. When developing individual trees, an arbitrary subset of attributes is drawn (hence the term "Random"), from which the best attribute for the split is selected. The final model is based on the majority vote from individually developed trees in the forest.

**Random Forest** works for both classification and regression tasks.

![](images/RandomForest.png)

1. Specify the name of the model. The default name is "Random Forest".
2. Basic properties:
   - *Number of trees*: Specify how many decision trees will be included in the forest.
   - *Number of trees considered at each split*: Specify how many attributes will be arbitrarily drawn for consideration at each node. If the latter is not specified (option *Number of attributes...* left unchecked), this number is equal to the square root of the number of attributes in the data.
   - *Replicable training*: Fix the seed for tree generation, which enables replicability of the results.
   - *Balance class distribution*: [Weigh classes](https://scikit-learn.org/stable/modules/generated/sklearn.utils.class_weight.compute_class_weight.html?highlight=sklearn%20utils%20class_weight) inversely proportional to their frequencies.
3. Growth control:
   - *Limit depth of individual trees*: Original Breiman's proposal is to grow the trees without any pre-pruning, but since pre-pruning often works quite well and is faster, the user can set the depth to which the trees will be grown.
   - *Do not split subsets smaller than*: Select the smallest subset that can be split.
4. Click *Apply* to communicate the changes to other widgets. Alternatively, tick the box on the left side of the *Apply* button and changes will be communicated automatically.

Preprocessing
-------------

Random Forest uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- continuizes categorical variables (with one-hot-encoding)
- removes empty columns
- imputes missing values with mean values

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Feature Scoring
---------------

Random Forest can be used with Rank for feature scoring. See [Learners as Scorers](../../learners-as-scorers/index.md) for an example.

Examples
--------

For classification tasks, we use *iris* dataset. Connect it to [Predictions](../evaluate/predictions.md). Then, connect [File](../data/file.md) to **Random Forest** and [Tree](../model/tree.md) and connect them further to [Predictions](../evaluate/predictions.md). Finally, observe the predictions for the two models.

![](images/RandomForest-classification.png)

For regressions tasks, we will use *housing* data. Here, we will compare different models, namely **Random Forest**, [Linear Regression](../model/linearregression.md) and [Constant](../model/constant.md), in the [Test & Score](../evaluate/testandscore.md) widget.

![](images/RandomForest-regression.png)

References
----------

Breiman, L. (2001). Random Forests. In Machine Learning, 45(1), 5-32. Available [here](https://www.stat.berkeley.edu/~breiman/randomforest2001.pdf).
Stacking
========

Stack multiple models.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)
- Learners: learning algorithm
- Aggregate: model aggregation method

**Outputs**

- Learner: aggregated (stacked) learning algorithm
- Model: trained model

**Stacking** is an ensemble method that computes a meta model from several base models. The **Stacking** widget has the **Aggregate** input, which provides a method for aggregating the input models. If no aggregation input is given the default methods are used. Those are **Logistic Regression** for classification and **Ridge Regression** for regression problems.

![](images/Stacking-stamped.png)

1. The meta learner can be given a name under which it will appear in other widgets. The default name is “Stack”.
2. Click *Apply* to commit the aggregated model. That will put the new learner in the output and, if the training examples are given, construct a new model and output it as well. To communicate changes automatically tick *Apply Automatically*.
3. Access help and produce a report.

Example
-------

We will use [Paint Data](../data/paintdata.md) to demonstrate how the widget is used. We painted a complex dataset with 4 class labels and sent it to [Test & Score](../evaluate/testandscore.md). We also provided three [kNN](../model/knn.md) learners, each with a different parameters (number of neighbors is 5, 10 or 15). Evaluation results are good, but can we do better?

Let's use **Stacking**. **Stacking** requires several learners on the input and an aggregation method. In our case, this is [Logistic Regression](../model/logisticregression.md). A constructed meta learner is then sent to **Test & Score**. Results have improved, even if only marginally. **Stacking** normally works well on complex data sets.

![](images/Stacking-Example.png)
kNN
===

Predict according to the nearest training instances.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: kNN learning algorithm
- Model: trained model

The **kNN** widget uses the [kNN algorithm](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm) that searches for k closest training examples in feature space and uses their average as prediction.

![](images/kNN-stamped.png)

1. A name under which it will appear in other widgets. The default name is "kNN".
2. Set the number of nearest neighbors, the distance parameter (metric) and weights as model criteria.
    - Metric can be:
        - [Euclidean](https://en.wikipedia.org/wiki/Euclidean_distance) ("straight line", distance between two points)
        - [Manhattan](https://en.wikipedia.org/wiki/Taxicab_geometry) (sum of absolute differences of all attributes)
        - [Maximal](https://en.wikipedia.org/wiki/Chebyshev_distance) (greatest of absolute differences between attributes)
        - [Mahalanobis](https://en.wikipedia.org/wiki/Mahalanobis_distance) (distance between point and distribution).
    - The *Weights* you can use are:
        - **Uniform**: all points in each neighborhood are weighted equally.
        - **Distance**: closer neighbors of a query point have a greater influence than the neighbors further away.
3. Produce a report.
4. When you change one or more settings, you need to click *Apply*, which will put a new learner on the output and, if the training examples are given, construct a new model and output it as well. Changes can also be applied automatically by clicking the box on the left side of the *Apply* button.

Preprocessing
-------------

kNN uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- continuizes categorical variables (with one-hot-encoding)
- removes empty columns
- imputes missing values with mean values
- normalizes the data by centering to mean and scaling to standard deviation of 1

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Examples
--------

The first example is a classification task on *iris* dataset. We compare the results of [k-Nearest neighbors](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm) with the default model [Constant](../model/constant.md), which always predicts the majority class.

![](images/Constant-classification.png)

The second example is a regression task. This workflow shows how to use the *Learner* output. For the purpose of this example, we used the *housing* dataset. We input the **kNN** prediction model into [Predictions](../evaluate/predictions.md) and observe the predicted values.

![](images/kNN-regression.png)
Tree
====

A tree algorithm with forward pruning.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: decision tree learning algorithm
- Model: trained model

**Tree** is a simple algorithm that splits the data into nodes by class purity (information gain for categorical and MSE for numeric target variable). It is a precursor to [Random Forest](../model/randomforest.md). Tree in Orange is designed in-house and can handle both categorical and numeric datasets.

It can also be used for both classification and regression tasks.

![](images/Tree-stamped.png)

1. The learner can be given a name under which it will appear in other widgets. The default name is "Tree".
2. Tree parameters:
   - **Induce binary tree**: build a binary tree (split into two child nodes)
   - **Min. number of instances in leaves**: if checked, the algorithm will never construct a split which would put less than the specified number of training examples into any of the branches.
   - **Do not split subsets smaller than**: forbids the algorithm to split the nodes with less than the given number of instances.
   - **Limit the maximal tree depth**: limits the depth of the classification tree to the specified number of node levels.
3. **Stop when majority reaches [%]**: stop splitting the nodes after a specified majority threshold is reached
4. Produce a report. After changing the settings, you need to click *Apply*, which will put the new learner on the output and, if the training examples are given, construct a new classifier and output it as well. Alternatively, tick the box on the left and changes will be communicated automatically.

Preprocessing
-------------

Tree does not use any preprocessing.

Examples
--------

There are two typical uses for this widget. First, you may want to induce a model and check what it looks like in [Tree Viewer](../visualize/treeviewer.md).

![](images/Tree-classification-visualize.png)

The second schema trains a model and evaluates its performance against [Logistic Regression](../model/logisticregression.md).

![](images/Tree-classification-model.png)

We used the *iris* dataset in both examples. However, **Tree** works for regression tasks as well. Use *housing* dataset and pass it to **Tree**. The selected tree node from [Tree Viewer](../visualize/treeviewer.md) is presented in the [Scatter Plot](../visualize/scatterplot.md) and we can see that the selected examples exhibit the same features.

![](images/Tree-regression-subset.png)
Logistic Regression
===================

The logistic regression classification algorithm with LASSO (L1) or ridge (L2) regularization.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: logistic regression learning algorithm
- Model: trained model
- Coefficients: logistic regression coefficients

**Logistic Regression** learns a [Logistic Regression](https://en.wikipedia.org/wiki/Logistic_regression) model from the data. It only works for classification tasks.

![](images/LogisticRegression-stamped.png)

1. A name under which the learner appears in other widgets. The default name is "Logistic Regression".
2. [Regularization](https://en.wikipedia.org/wiki/Regularization_(mathematics)) type (either [L1](https://en.wikipedia.org/wiki/Least_squares#Lasso_method) or [L2](https://en.wikipedia.org/wiki/Tikhonov_regularization)). Set the cost strength (default is C=1).
3. Press *Apply* to commit changes. If *Apply Automatically* is ticked, changes will be communicated automatically.

Preprocessing
-------------

Logistic Regression uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- continuizes categorical variables (with one-hot-encoding)
- removes empty columns
- imputes missing values with mean values

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Feature Scoring
---------------

Logistic Regression can be used with Rank for feature scoring. See [Learners as Scorers](../../learners-as-scorers/index.md) for an example.

Example
-------

The widget is used just as any other widget for inducing a classifier. This is an example demonstrating prediction results with logistic regression on the *hayes-roth* dataset. We first load *hayes-roth_learn* in the [File](../data/file.md) widget and pass the data to **Logistic Regression**. Then we pass the trained model to [Predictions](../evaluate/predictions.md).

Now we want to predict class value on a new dataset. We load *hayes-roth_test* in the second **File** widget and connect it to **Predictions**. We can now observe class values predicted with **Logistic Regression** directly in **Predictions**.

![](images/LogisticRegression-classification.png)
Calibrated Learner
==================

Wraps another learner with probability calibration and decision threshold optimization.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)
- Base Learner: learner to calibrate

**Outputs**

- Learner: calibrated learning algorithm
- Model: trained model using the calibrated learner

This learner produces a model that calibrates the distribution of class probabilities and optimizes decision threshold. The widget works only for binary classification tasks.

![](images/Calibrated-Learner-stamped.png)

1. The name under which it will appear in other widgets. Default name is composed of the learner, calibration and optimization parameters.
2. Probability calibration:

   - [Sigmoid calibration](http://citeseer.ist.psu.edu/viewdoc/summary?doi=10.1.1.41.1639)
   - [Isotonic calibration](https://scikit-learn.org/stable/auto_examples/plot_isotonic_regression.html)
   - No calibration

3. Decision threshold optimization:

   - Optimize classification accuracy
   - Optimize F1 score
   - No threshold optimization

4. Press *Apply* to commit changes. If *Apply Automatically* is ticked, changes are committed automatically.

Example
-------

A simple example with **Calibrated Learner**. We are using the *titanic* data set as the widget requires binary class values (in this case they are 'survived' and 'not survived').

We will use [Logistic Regression](logisticregression.md) as the base learner which will we calibrate with the default settings, that is with sigmoid optimization of distribution values and by optimizing the CA.

Comparing the results with the uncalibrated **Logistic Regression** model we see that the calibrated model performs better.

![](images/Calibrated-Learner-Example.png)
Curve Fit
=========

Fit a function to data.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: curve fit learning algorithm
- Model: trained model
- Coefficients: fitted coefficients

The **Curve Fit** widget fits an arbitrary function to the input data. It only works for regression tasks.
The widget uses [scipy.curve_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html) to find the optimal values of the parameters.

The widget works only on regression tasks and only numerical features can be used for fitting.

![](images/CurveFit-stamped.png)

1. The learner/predictor name.
2. Introduce model parameters.
3. Input an expression in Python. The expression should consist of at least one fitting parameter.
4. Select a feature to include into the expression. Only numerical features are available.
5. Select a parameter. Only the introduced parameters are available.
6. Select a function.
7. Press *Apply* to commit changes. If *Apply Automatically* is ticked, changes are committed automatically.
8. Show help, produce a report, input/output info.

Preprocessing
-------------

Curve fit uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- removes empty columns
- imputes missing values with mean values

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Example
-------

Below, is a simple workflow with *housing* dataset. Due to example simplicity we used only a single feature. Unlike the other modelling widgets, the Curve Fit needs data on the input. We trained **Curve Fit** and [Linear Regression](../model/linearregression.md) and evaluated their performance in [Test & Score](../evaluate/testandscore.md).

![](images/CurveFit-example.png)
Stochastic Gradient Descent
===========================

Minimize an objective function using a stochastic approximation of gradient descent.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: stochastic gradient descent learning algorithm
- Model: trained model

The **Stochastic Gradient Descent** widget uses [stochastic gradient descent](https://en.wikipedia.org/wiki/Stochastic_gradient_descent) that minimizes a chosen loss function with a linear function. The algorithm approximates a true gradient by considering one sample at a time, and simultaneously updates the model based on the gradient of the loss function. For regression, it returns predictors as minimizers of the sum, i.e. M-estimators, and is especially useful for large-scale and sparse datasets.

![](images/StochasticGradientDescent-stamped.png)

1. Specify the name of the model. The default name is "SGD".
2. Algorithm parameters:
    - Classification loss function:
        - [Hinge](https://en.wikipedia.org/wiki/Hinge_loss) (linear SVM)
        - [Logistic Regression](http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html#sklearn.linear_model.LogisticRegression) (logistic regression SGD)
        - [Modified Huber](https://en.wikipedia.org/wiki/Huber_loss) (smooth loss that brings tolerance to outliers as well as probability estimates)
        - *Squared Hinge* (quadratically penalized hinge)
        - [Perceptron](http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Perceptron.html#sklearn.linear_model.Perceptron) (linear loss used by the perceptron algorithm)
        - [Squared Loss](https://en.wikipedia.org/wiki/Mean_squared_error#Regression) (fitted to ordinary least-squares)
        - [Huber](https://en.wikipedia.org/wiki/Huber_loss) (switches to linear loss beyond ε)
        - [Epsilon insensitive](http://kernelsvm.tripod.com/) (ignores errors within ε, linear beyond it)
        - *Squared epsilon insensitive* (loss is squared beyond ε-region).
    - Regression loss function:
        - [Squared Loss](https://en.wikipedia.org/wiki/Mean_squared_error#Regression) (fitted to ordinary least-squares)
        - [Huber](https://en.wikipedia.org/wiki/Huber_loss) (switches to linear loss beyond ε)
        - [Epsilon insensitive](http://kernelsvm.tripod.com/) (ignores errors within ε, linear beyond it)
        - *Squared epsilon insensitive* (loss is squared beyond ε-region).
3. Regularization norms to prevent overfitting:
   - None.
   - [Lasso (L1)](https://en.wikipedia.org/wiki/Taxicab_geometry) (L1 leading to sparse solutions)
   - [Ridge (L2)](https://en.wikipedia.org/wiki/Norm_(mathematics)#p-norm) (L2, standard regularizer)
   - [Elastic net](https://en.wikipedia.org/wiki/Elastic_net_regularization) (mixing both penalty norms).

   Regularization strength defines how much regularization will be applied (the less we regularize, the more we allow the model to fit the data) and the mixing parameter what the ratio between L1 and L2 loss will be (if set to 0 then the loss is L2, if set to 1 then it is L1).
4. Learning parameters.
   - Learning rate:
      - *Constant*: learning rate stays the same through all epochs (passes)
      - [Optimal](http://leon.bottou.org/projects/sgd): a heuristic proposed by Leon Bottou
      - [Inverse scaling](http://users.ics.aalto.fi/jhollmen/dippa/node22.html): earning rate is inversely related to the number of iterations
   - Initial learning rate.
   - Inverse scaling exponent: learning rate decay.
   - Number of iterations: the number of passes through the training data.
   - If *Shuffle data after each iteration* is on, the order of data instances is mixed after each pass.
   - If *Fixed seed for random shuffling* is on, the algorithm will use a fixed random seed and enable replicating the results.
5. Produce a report.
6. Press *Apply* to commit changes. Alternatively, tick the box on the left side of the *Apply* button and changes will be communicated automatically.

Preprocessing
-------------

SGD uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- continuizes categorical variables (with one-hot-encoding)
- removes empty columns
- imputes missing values with mean values
- normalizes the data by centering to mean and scaling to standard deviation of 1

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Feature Scoring
---------------

Stochastic Gradient Descent can be used with Rank for feature scoring. See [Learners as Scorers](../../learners-as-scorers/index.md) for an example.

Examples
--------

For the classification task, we will use *iris* dataset and test two models on it. We connected [Stochastic Gradient Descent](../model/stochasticgradient.md) and [Tree](../model/tree.md) to [Test & Score](../evaluate/testandscore.md). We also connected [File](../data/file.md) to **Test & Score** and observed model performance in the widget.

![](images/StochasticGradientDescent-classification.png)

For the regression task, we will compare three different models to see which predict what kind of results. For the purpose of this example, the *housing* dataset is used. We connect the [File](../data/file.md) widget to **Stochastic Gradient Descent**, [Linear Regression](../model/linearregression.md) and [kNN](../model/knn.md) widget and all four to the [Predictions](../evaluate/predictions.md) widget.

![](images/StochasticGradientDescent-regression.png)
Constant
========

Predict the most frequent class or mean value from the training set.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: majority/mean learning algorithm
- Model: trained model

This learner produces a model that always predicts the[majority](https://en.wikipedia.org/wiki/Predictive_modelling#Majority_classifier) for classification tasks and [mean value](https://en.wikipedia.org/wiki/Mean) for regression tasks.

For classification, when predicting the class value with [Predictions](../evaluate/predictions.md), the widget will return relative frequencies of the classes in the training set. When there are two or more majority classes, the classifier chooses the predicted class randomly, but always returns the same class for a particular example.

For regression, it *learns* the mean of the class variable and returns a predictor with the same mean value.

The widget is typically used as a baseline for other models.

![](images/Constant-stamped.png)

This widget provides the user with two options:

1. The name under which it will appear in other widgets. Default name is "Constant".
2. Produce a report.

If you change the widget's name, you need to click *Apply*. Alternatively, tick the box on the left side and changes will be communicated automatically.

Preprocessing
-------------

Constant does not use any preprocessing.

Examples
--------

In a typical classification example, we would use this widget to compare the scores of other learning algorithms (such as kNN) with the default scores. Use *iris* dataset and connect it to [Test & Score](../evaluate/testandscore.md). Then connect **Constant** and [kNN](../model/knn.md) to [Test & Score](../evaluate/testandscore.md) and observe how well [kNN](../model/knn.md) performs against a constant baseline.

![](images/Constant-classification.png)

For regression, we use **Constant** to construct a predictor in [Predictions](../evaluate/predictions.md). We used the *housing* dataset. In **Predictions**, you can see that *Mean Learner* returns one (mean) value for all instances.

![](images/Constant-regression.png)
Neural Network
==============

A multi-layer perceptron (MLP) algorithm with backpropagation.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: multi-layer perceptron learning algorithm
- Model: trained model

The **Neural Network** widget uses sklearn's [Multi-layer Perceptron algorithm](http://scikit-learn.org/stable/modules/neural_networks_supervised.html) that can learn non-linear models as well as linear.

![](images/NeuralNetwork-stamped.png)

1. A name under which it will appear in other widgets. The default name is "Neural Network".
2. Set model parameters:
   - Neurons per hidden layer: defined as the ith element represents the number of neurons in the ith hidden layer. E.g. a neural network with 3 layers can be defined as 2, 3, 2.
   - Activation function for the hidden layer:
      - Identity: no-op activation, useful to implement linear bottleneck
      - Logistic: the logistic sigmoid function
      - tanh: the hyperbolic tan function
      - ReLu: the rectified linear unit function
   - Solver for weight optimization:
      - L-BFGS-B: an optimizer in the family of quasi-Newton methods
      - SGD: stochastic gradient descent
      - Adam: stochastic gradient-based optimizer
   - Alpha: L2 penalty (regularization term) parameter
   - Max iterations: maximum number of iterations

   Other parameters are set to [sklearn's defaults](http://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPClassifier.html).
3. Produce a report.
4. When the box is ticked (*Apply Automatically*), the widget will communicate changes automatically. Alternatively, click *Apply*.

Preprocessing
-------------

Neural Network uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- continuizes categorical variables (with one-hot-encoding)
- removes empty columns
- imputes missing values with mean values
- normalizes the data by centering to mean and scaling to standard deviation of 1

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Examples
--------

The first example is a classification task on *iris* dataset. We compare the results of **Neural Network** with the [Logistic Regression](../model/logisticregression.md).

![](images/NN-Example-Test.png)

The second example is a prediction task, still using the *iris* data. This workflow shows how to use the *Learner* output. We input the **Neural Network** prediction model into [Predictions](../evaluate/predictions.md) and observe the predicted values.

![](images/NN-Example-Predict.png)
CN2 Rule Induction
==================

Induce rules from data using CN2 algorithm.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: CN2 learning algorithm
- CN2 Rule Classifier: trained model

The CN2 algorithm is a classification technique designed for the efficient induction of simple, comprehensible rules of form "if *cond* then predict *class*", even in domains where noise may be present.

**CN2 Rule Induction** works only for classification.

![](images/CN2-stamped.png)

1. Name under which the learner appears in other widgets. The default name is *CN2 Rule Induction*.
2. *Rule ordering*:
   - **Ordered**: induce ordered rules (decision list). Rule conditions are found and the majority class is assigned in the rule head.
   - **Unordered**: induce unordered rules (rule set). Learn rules for each class individually, in regard to the original learning data.
3. *Covering algorithm*:
   - **Exclusive**: after covering a learning instance, remove it from further consideration.
   - **Weighted**: after covering a learning instance, decrease its weight (multiplication by *gamma*) and in-turn decrease its impact on further iterations of the algorithm.
4. *Rule search*:
   - **Evaluation measure**: select a heuristic to evaluate found hypotheses:
     - [Entropy](https://en.wikipedia.org/wiki/Entropy_(information_theory)) (measure of unpredictability of content)
     - [Laplace Accuracy](https://en.wikipedia.org/wiki/Laplace%27s_method)
     - Weighted Relative Accuracy
   - **Beam width**; remember the best rule found thus far and monitor a fixed number of alternatives (the beam).
5. *Rule filtering*:
   - **Minimum rule coverage**: found rules must cover at least the minimum required number of covered examples. Unordered rules must cover this many target class examples.
   - **Maximum rule length**: found rules may combine at most the maximum allowed number of selectors (conditions).
   - **Default alpha**: significance testing to prune out most specialised (less frequently applicable) rules in regard to the initial distribution of classes.
   - **Parent alpha**: significance testing to prune out most specialised (less frequently applicable) rules in regard to the parent class distribution.
6. Tick 'Apply Automatically' to auto-communicate changes to other widgets and to immediately train the classifier if learning data is connected. Alternatively, press ‘Apply‘ after configuration.

Preprocessing
-------------

CN2 Rule Induction uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes empty columns
- removes instances with unknown target values
- imputes missing values with mean values

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Examples
--------

For the example below, we have used *zoo* dataset and passed it to **CN2 Rule Induction**. We can review and interpret the built model with [CN2 Rule Viewer](../visualize/cn2ruleviewer.md) widget.

![](images/CN2-visualize.png)

The second workflow tests evaluates **CN2 Rule Induction** and [Tree](../model/tree.md) in [Test & Score](../evaluate/testandscore.md).

![](images/CN2-classification.png)

References
----------

1. Fürnkranz, Johannes. "Separate-and-Conquer Rule Learning", Artificial Intelligence Review 13, 3-54, 1999.
2. Clark, Peter and Tim Niblett. "The CN2 Induction Algorithm", Machine Learning Journal, 3 (4), 261-283, 1989.
3. Clark, Peter and Robin Boswell. "Rule Induction with CN2: Some Recent Improvements", Machine Learning - Proceedings of the 5th European Conference (EWSL-91),151-163, 1991.
4. Lavrač, Nada et al. "Subgroup Discovery with CN2-SD",Journal of Machine Learning Research 5, 153-188, 2004
Save Model
==========

Save a trained model to an output file.

If the file is saved to the same directory as the workflow or in the subtree of that directory, the widget remembers the relative path. Otherwise it will store an absolute path, but disable auto save for security reasons.

**Inputs**

- Model: trained model

![](images/SaveModel-stamped.png)

1. Choose from previously saved models.
2. Save the created model with the *Browse* icon. Click on the icon and enter the name of the file. The model will be saved to a pickled file.
![](images/SaveModel-save.png)
3. Save the model.

Example
-------

When you want to save a custom-set model, feed the data to the model (e.g. [Logistic Regression](../model/logisticregression.md)) and connect it to **Save Model**. Name the model; load it later into workflows with [Load Model](../model/loadmodel.md). Datasets used with **Load Model** have to contain compatible attributes.

![](images/SaveModel-example.png)
SVM
===

Support Vector Machines map inputs to higher-dimensional feature spaces.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: linear regression learning algorithm
- Model: trained model
- Support Vectors: instances used as support vectors

[Support vector machine](https://en.wikipedia.org/wiki/Support_vector_machine) (SVM) is a machine learning technique that separates the attribute space with a hyperplane, thus maximizing the margin between the instances of different classes or class values. The technique often yields supreme predictive performance results. Orange embeds a popular implementation of SVM from the [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) package. This widget is its graphical user interface.

For regression tasks, **SVM** performs linear regression in a high dimension feature space using an ε-insensitive loss. Its estimation accuracy depends on a good setting of C, ε and kernel parameters. The widget outputs class predictions based on a [SVM Regression](https://en.wikipedia.org/wiki/Support_vector_machine#Regression).

The widget works for both classification and regression tasks.

![](images/SVM-stamped.png)

1. The learner can be given a name under which it will appear in other widgets. The default name is "SVM".
2. SVM type with test error settings. *SVM* and *ν-SVM* are based on different minimization of the error function. On the right side, you can set test error bounds:
   - [SVM](http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVR.html):
      - [Cost](http://www.quora.com/What-are-C-and-gamma-with-regards-to-a-support-vector-machine): penalty term for loss and applies for classification and regression tasks.
      - ε: a parameter to the epsilon-SVR model, applies to regression tasks. Defines the distance from true values within which no penalty is associated with predicted values.
   - [ν-SVM](http://scikit-learn.org/stable/modules/generated/sklearn.svm.NuSVR.html#sklearn.svm.NuSVR):
      - [Cost](http://www.quora.com/What-are-C-and-gamma-with-regards-to-a-support-vector-machine): penalty term for loss and applies only to regression tasks
      - ν: a parameter to the ν-SVR model, applies to classification and regression tasks. An upper bound on the fraction of training errors and a lower bound of the fraction of support vectors.
3. Kernel is a function that transforms attribute space to a new feature space to fit the maximum-margin hyperplane, thus allowing the algorithm to create the model with [Linear](https://en.wikipedia.org/wiki/Linear_model), [Polynomial](https://en.wikipedia.org/wiki/Polynomial_kernel), [RBF](https://en.wikipedia.org/wiki/Radial_basis_function_kernel) and [Sigmoid](http://crsouza.com/2010/03/kernel-functions-for-machine-learning-applications/#sigmoid) kernels. Functions that specify the kernel are presented upon selecting them, and the constants involved are:
   - **g** for the gamma constant in kernel function (the recommended value is 1/k, where k is the number of the attributes, but since there may be no training set given to the widget the default is 0 and the user has to set this option manually),
   - **c** for the constant c0 in the kernel function (default 0), and
   - **d** for the degree of the kernel (default 3).
4. Set permitted deviation from the expected value in *Numerical Tolerance*. Tick the box next to *Iteration Limit* to set the maximum number of iterations permitted.
5. Produce a report.
6. Click *Apply* to commit changes. If you tick the box on the left side of the *Apply* button, changes will be communicated automatically.

Preprocessing
-------------

SVM uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- continuizes categorical variables (with one-hot-encoding)
- removes empty columns
- imputes missing values with mean values

For classification, SVM also normalizes dense and scales sparse data.

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Examples
--------

In the first (regression) example, we have used *housing* dataset and split the data into two data subsets (*Data Sample* and *Remaining Data*) with [Data Sampler](../data/datasampler.md). The sample was sent to SVM which produced a *Model*, which was then used in [Predictions](../evaluate/predictions.md) to predict the values in *Remaining Data*. A similar schema can be used if the data is already in two separate files; in this case, two [File](../data/file.md) widgets would be used instead of the [File](../data/file.md) - [Data Sampler](../data/datasampler.md) combination.

![](images/SVM-Predictions.png)

The second example shows how to use **SVM** in combination with [Scatter Plot](../visualize/scatterplot.md). The following workflow trains a SVM model on *iris* data and outputs support vectors, which are those data instances that were used as support vectors in the learning phase. We can observe which are these data instances in a scatter plot visualization. Note that for the workflow to work correctly, you must set the links between widgets as demonstrated in the screenshot below.

![](images/SVM-support-vectors.png)

References
----------

[Introduction to SVM on StatSoft](http://www.statsoft.com/Textbook/Support-Vector-Machines).
Load Model
==========

Load a model from an input file.

**Outputs**

- Model: trained model

![](images/LoadModel-stamped.png)

1. Choose from a list of previously used models.
2. Browse for saved models.
3. Reload the selected model.

Example
-------

When you want to use a custom-set model that you've saved before, open the **Load Model** widget and select the desired file with the *Browse* icon. This widget loads the existing model into [Predictions](../evaluate/predictions.md) widget. Datasets used with **Load Model** have to contain compatible attributes!

![](images/LoadModel-example.png)
Naive Bayes
===========

A fast and simple probabilistic classifier based on Bayes' theorem with the assumption of feature independence.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: naive bayes learning algorithm
- Model: trained model

**Naive Bayes** learns a [Naive Bayesian](https://en.wikipedia.org/wiki/Naive_Bayes_classifier) model from the data. It only works for classification tasks.

![](images/NaiveBayes-stamped.png)

This widget has two options: the name under which it will appear in other widgets and producing a report. The default name is *Naive Bayes*. When you change it, you need to press *Apply*.

Preprocessing
-------------

Naive Bayes uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes empty columns
- discretizes numeric values to 4 bins with equal frequency

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Examples
--------

Here, we present two uses of this widget. First, we compare the results of the
**Naive Bayes** with another model, the [Random Forest](../model/randomforest.md). We connect *iris* data from [File](../data/file.md) to [Test & Score](../evaluate/testandscore.md). We also connect **Naive Bayes** and [Random Forest](../model/randomforest.md) to **Test & Score** and observe their prediction scores.

![](images/NaiveBayes-classification.png)

The second schema shows the quality of predictions made with **Naive Bayes**. We feed the [Test & Score](../evaluate/testandscore.md) widget a Naive Bayes learner and then send the data to the [Confusion Matrix](../evaluate/confusionmatrix.md). We also connect [Scatter Plot](../visualize/scatterplot.md) with **File**. Then we select the misclassified instances in the **Confusion Matrix** and show feed them to [Scatter Plot](../visualize/scatterplot.md). The bold dots in the scatterplot are the misclassified instances from **Naive Bayes**.

![](images/NaiveBayes-visualize.png)
Gradient Boosting
=================

Predict using gradient boosting on decision trees.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: gradient boosting learning algorithm
- Model: trained model

[Gradient Boosting](https://en.wikipedia.org/wiki/Gradient_boosting) is a machine learning technique for regression and classification problems, which produces a prediction model in the form of an ensemble of weak prediction models, typically decision trees.

![](images/GradientBoosting-stamped.png)

1. Specify the name of the model. The default name is "Gradient Boosting".
2. Select a gradient boosting method:
   - [Gradient Boosting (scikit-learn)](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.GradientBoostingClassifier.html)
   - [Extreme Gradient Boosting (xgboost)](https://xgboost.readthedocs.io/en/latest/index.html)
   - [Extreme Gradient Boosting Random Forest (xgboost)](https://xgboost.readthedocs.io/en/latest/index.html)
   - [Gradient Boosting (catboost)](https://catboost.ai/docs/concepts/python-quickstart.html)
3. Basic properties:
   - *Number of trees*: Specify how many gradient boosted trees will be included. A large number usually results in better performance.
   - *Learning rate*: Specify the boosting learning rate. Learning rate shrinks the contribution of each tree.
   - *Replicable training*: Fix the random seed, which enables replicability of the results.
   - *Regularization*: Specify the L2 regularization term. Available only for *xgboost* and *catboost* methods.
4. Growth control:
   - *Limit depth of individual trees*: Specify the maximum depth of the individual tree.
   - *Do not split subsets smaller than*: Specify the smallest subset that can be split. Available only for *scikit-learn* methods.
5. Subsampling:
   - *Fraction of training instances*: Specify the percentage of the training instances for fitting the individual tree. Available for *scikit-learn* and *xgboost* methods.
   - *Fraction of features for each tree*: Specify the percentage of features to use when constructing each tree. Available for *xgboost* and *catboost* methods.
   - *Fraction of features for each level*: Specify the percentage of features to use for each level. Available only for *xgboost* methods.
   - *Fraction of features for each split*: Specify the percentage of features to use for each split. Available only for *xgboost* methods.
6. Click *Apply* to communicate the changes to other widgets. Alternatively, tick the box on the left side of the *Apply* button and changes will be communicated automatically.

Preprocessing
-------------

Gradient Boosting uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- continuizes categorical variables (with one-hot-encoding)
- removes empty columns
- imputes missing values with mean values

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Feature Scoring
---------------

Gradient Boosting can be used with Rank for feature scoring. See [Learners as Scorers](../../learners-as-scorers/index.md) for an example.

Example
-------

For a classification tasks, we use the *heart disease* data. Here, we compare all available methods in the [Test & Score](../evaluate/testandscore.md) widget.

![](images/GradientBoosting-example.png)
Linear Regression
=================

A linear regression algorithm with optional L1 (LASSO), L2 (ridge) or L1L2 (elastic net) regularization.

**Inputs**

- Data: input dataset
- Preprocessor: preprocessing method(s)

**Outputs**

- Learner: linear regression learning algorithm
- Model: trained model
- Coefficients: linear regression coefficients

The **Linear Regression** widget constructs a learner/predictor that learns a [linear function](https://en.wikipedia.org/wiki/Linear_regression) from its input data. The model can identify the relationship between a predictor xi and the response variable y. Additionally, [Lasso](https://en.wikipedia.org/wiki/Least_squares#Lasso_method) and [Ridge](https://en.wikipedia.org/wiki/Least_squares#Lasso_method) regularization parameters can be specified. Lasso regression minimizes a penalized version of the least squares loss function with L1-norm penalty and Ridge regularization with L2-norm penalty.

Linear regression works only on regression tasks.

![](images/LinearRegression-stamped.png)

1. The learner/predictor name
2. Choose a model to train:
   - no regularization
   - a [Ridge](https://en.wikipedia.org/wiki/Least_squares#Lasso_method) regularization (L2-norm penalty)
   - a [Lasso](https://en.wikipedia.org/wiki/Least_squares#Lasso_method) bound (L1-norm penalty)
   - an [Elastic net](https://en.wikipedia.org/wiki/Elastic_net_regularization) regularization
3. Produce a report.
4. Press *Apply* to commit changes. If *Apply Automatically* is ticked, changes are committed automatically.

Preprocessing
-------------

Linear Regression uses default preprocessing when no other preprocessors are given. It executes them in the following order:

- removes instances with unknown target values
- continuizes categorical variables (with one-hot-encoding)
- removes empty columns
- imputes missing values with mean values

To remove default preprocessing, connect an empty [Preprocess](../data/preprocess.md) widget to the learner.

Feature Scoring
---------------

Linear Regression can be used with Rank for feature scoring. See [Learners as Scorers](../../learners-as-scorers/index.md) for an example.

Example
-------

Below, is a simple workflow with *housing* dataset. We trained **Linear Regression** and [Random Forest](../model/randomforest.md) and evaluated their performance in [Test & Score](../evaluate/testandscore.md).

![](images/LinearRegression-regression.png)
# Report

It is possible to compile a report in Orange. We can save the report in .html, .pdf or .report format. Reports allow us to trace back analytical steps as it saves the workflow at which each report segment was created.

Each widget has a report button in the status bar at the bottom. Pressing on the the File icon adds a new section to the report.

![](report-button.png)

Report can be examined with View - Show report.

## Simple example

We built a simple workflow with File and Scatter Plot, adding a section to the report at each step. Widgets report parameters, visualizations, and other settings. Each section includes a comment for extra explanation.

![](report.png)

To remove a report section, hover on the section in the list on the left. A Trash and an Orange icon will appear. The trash icon removes the section from the report list. Orange icon loads the workflow as it was at the time of creating the section. This is very handy if a colleague wishes to inspect the results. This option is available only if the report is saved in .report format.
# Exporting Models

Predictive models can be saved and re-used. Models are saved in Python [pickle](https://docs.python.org/3/library/pickle.html) format.

![](load-save-model.png)

## Save model

Models first require data for training. They output a trained model, which can be saved with [Save Model](../widgets/model/savemodel.md) widget in the pickle format.

## Load model

Models can be reused in different Orange workflows. [Load Model](../widgets/model/loadmodel.md) loads a trained model, which can be used in [Predictions](../widgets/evaluate/predictions.md) and elsewhere.

## Load in Python

Models can also be imported directly into Python and used in a script.

```python
import pickle

with open('model.pkcls', 'rb') as model:
    lr = pickle.loads(model)

lr
>> LogisticRegressionClassifier(skl_model=LogisticRegression(C=1,
                                class_weight=None, dual=False, 
                                fit_intercept=True, intercept_scaling=1.0, 
                                l1_ratio=None, max_iter=10000, 
                                multi_class='auto', n_jobs=1, penalty='l2', 
                                random_state=0, solver='lbfgs', tol=0.0001, 
                                verbose=0, warm_start=False))
```
# Exporting Visualizations

Visualizations are an essential part of data science, and analytical reports are incomplete without them. Orange provides a couple of options for saving and modifying visualizations.

At the bottom of each widget, there is a status bar. Visualization widgets have a Save icon (second from the left) and a Palette icon (fourth from the left). Save icon saves the plot to the computer. Palette icon opens a dialogue for modifying visualizations.

![](statusbar-viz.png)

## Saving a plot

Visualizations in Orange can be saved in several formats, namely .png, .svg, .pdf, .pdf from matplotlib and as a matplotlib Python code. A common option is saving in .svg (scalable vector graphic), which you can edit with a vector graphics software such as [Inkscape](https://inkscape.org/). Ctrl+C (cmd+C) will copy a .png plot, which you can import with ctrl+V (cmd+V) into Word, PowerPoint, or other software tools.

![](plot-format.png)

[Matplotlib](https://matplotlib.org/) Python code is ideal for detailed editing and a high customization level. Below is an example of the Python code. It is possible to adjust the colors, size of the symbols, markers, etc.

```python
import matplotlib.pyplot as plt
from numpy import array

plt.clf()

# data
x = array([1.4, 1.4, 1.3, 1.5, 1.4])
y = array([0.2, 0.7, 0.9, 0.2, 0.1])
# style
sizes = 13.5
edgecolors = ['#3a9ed0ff', '#c53a27ff']
edgecolors_index = array([0, 0, 1, 1, 1], dtype='int')
facecolors = ['#46befa80', '#ed462f80']
facecolors_index = array([0, 0, 1, 1, 1], dtype='int')
linewidths = 1.5
plt.scatter(x=x, y=y, s=sizes**2/4, marker='o',
            facecolors=array(facecolors)[facecolors_index], 
            edgecolors=array(edgecolors)[edgecolors_index],
            linewidths=linewidths)
plt.xlabel('petal length')
plt.ylabel('petal width')

plt.show()
```

## Modifying a plot

It is possible to modify certain parameters of a plot without digging into the code. Click on the Palette icon to open visual settings. One can change various attributes of the plot, such as fonts, font sizes, titles and so on.

![](plot-options.png)
