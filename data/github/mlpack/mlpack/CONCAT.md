### mlpack ?.?.?
###### ????-??-??
  * Migrate from boost tests to Catch2 framework (#2523), (#2584).

  * Bump minimum armadillo version from 8.400 to 9.800 (#3043), (#3048).

  * Adding a copy constructor in the Convolution layer (#3067).

  * Replace `boost::spirit` parser by a local efficient implementation (#2942).

  * Disable correctly the autodownloader + fix tests stability (#3076).

  * Replace `boost::any` with `core::v2::any` or `std::any` if available (#3006).

  * Remove old non used Boost headers (#3005).

  * Replace `boost::enable_if` with `std::enable_if` (#2998).

  * Replace `boost::is_same` with `std::is_same` (#2993).

  * Remove invalid option for emsmallen and STB (#2960).

  * Check for armadillo dependencies before downloading armadillo (#2954).

  * Disable the usage of autodownloader by default (#2953).

  * Install dependencies downloaded with the autodownloader (#2952).
 
  * Download older Boost if the compiler is old (#2940). 

  * Add support for embedded systems (#2531).

  * Build mlpack executable statically if the library is statically linked (#2931).

  * Fix cover tree loop bug on embedded arm systems (#2869).

  * Fix a LAPACK bug in `FindArmadillo.cmake` (#2929).

  * Add an autodownloader to get mlpack dependencies (#2927).

  * Remove Coverage files and configurations from CMakeLists (#2866).

  * Added `Multi Label Soft Margin Loss` loss function for neural networks
   (#2345).

  * Added Decision Tree Regressor (#2905). It can be used using the class
    `mlpack::tree::DecisionTreeRegressor`. It is accessible only though C++.

  * Added dict-style inspection of mlpack models in python bindings (#2868).

  * Added Extra Trees Algorithm (#2883). Currently, it can be used using the
    class `mlpack::tree::ExtraTrees`, but only through C++.

  * Add Flatten T Swish activation function (`flatten-t-swish.hpp`)

  * Added warm start feature to Random Forest (#2881); this feature is
    accessible from mlpack's bindings to different languages.

  * Added Pixel Shuffle layer (#2563).

  * Add "check_input_matrices" option to python bindings that checks
    for NaN and inf values in all the input matrices (#2787).

  * Add Adjusted R squared functionality to R2Score::Evaluate (#2624).

  * Disabled all the bindings by default in CMake (#2782).

  * Added an implementation to Stratify Data (#2671).

  * Add `BUILD_DOCS` CMake option to control whether Doxygen documentation is
    built (default ON) (#2730).

  * Add Triplet Margin Loss function (#2762).

  * Add finalizers to Julia binding model types to fix memory handling (#2756).

  * HMM: add functions to calculate likelihood for data stream with/without
    pre-calculated emission probability (#2142).

  * Replace Boost serialization library with Cereal (#2458).

  * Add `PYTHON_INSTALL_PREFIX` CMake option to specify installation root for
    Python bindings (#2797).

  * Removed `boost::visitor` from model classes for `knn`, `kfn`, `cf`,
    `range_search`, `krann`, and `kde` bindings (#2803).

  * Add k-means++ initialization strategy (#2813).

  * `NegativeLogLikelihood<>` now expects classes in the range `0` to
    `numClasses - 1` (#2534).

  * Add `Lambda1()`, `Lambda2()`, `UseCholesky()`, and `Tolerance()` members to
    `LARS` so parameters for training can be modified (#2861).

  * Remove unused `ElemType` template parameter from `DecisionTree` and
    `RandomForest` (#2874).

  * Fix Python binding build when the CMake variable `USE_OPENMP` is set to
    `OFF` (#2884).

  * The `mlpack_test` target is no longer built as part of `make all`.  Use
    `make mlpack_test` to build the tests.

  * Fixes to `HoeffdingTree`: ensure that training still works when empty
    constructor is used (#2964).

  * Fix Julia model serialization bug (#2970).

  * Fix `LoadCSV()` to use pre-populated `DatasetInfo` objects (#2980).

  * Add `probabilities` option to softmax regression binding, to get class
    probabilities for test points (#3001).

  * Fix thread safety issues in mlpack bindings to other languages (#2995).

  * Fix double-free of model pointers in R bindings (#3034).

  * Fix Julia, Python, R, and Go handling of categorical data for
    `decision_tree()` and `hoeffding_tree()` (#2971).

  * Depend on `pkgbuild` for R bindings (#3081).
 
  * Replaced Numpy deprecated code in Python bindings (#3126).

### mlpack 3.4.2
###### 2020-10-26
  * Added Mean Absolute Percentage Error.

  * Added Softmin activation function as layer in ann/layer.

  * Fix spurious ARMA_64BIT_WORD compilation warnings on 32-bit systems (#2665).

### mlpack 3.4.1
###### 2020-09-07
  * Fix incorrect parsing of required matrix/model parameters for command-line
    bindings (#2600).

  * Add manual type specification support to `data::Load()` and `data::Save()`
    (#2084, #2135, #2602).

  * Remove use of internal Armadillo functionality (#2596, #2601, #2602).

### mlpack 3.4.0
###### 2020-09-01
  * Issue warnings when metrics produce NaNs in KFoldCV (#2595).

  * Added bindings for _R_ during Google Summer of Code (#2556).

  * Added common striptype function for all bindings (#2556).

  * Refactored common utility function of bindings to bindings/util (#2556).

  * Renamed InformationGain to HoeffdingInformationGain in
    methods/hoeffding_trees/information_gain.hpp (#2556).

  * Added macro for changing stream of printing and warnings/errors (#2556).

  * Added Spatial Dropout layer (#2564).

  * Force CMake to show error when it didn't find Python/modules (#2568).

  * Refactor `ProgramInfo()` to separate out all the different
    information (#2558).

  * Add bindings for one-hot encoding (#2325).

  * Added Soft Actor-Critic to RL methods (#2487).

  * Added Categorical DQN to q_networks (#2454).

  * Added N-step DQN to q_networks (#2461).

  * Add Silhoutte Score metric and Pairwise Distances (#2406).

  * Add Go bindings for some missed models (#2460).

  * Replace boost program_options dependency with CLI11 (#2459).

  * Additional functionality for the ARFF loader (#2486); use case sensitive
    categories (#2516).

  * Add `bayesian_linear_regression` binding for the command-line, Python,
    Julia, and Go.  Also called "Bayesian Ridge", this is equivalent to a
    version of linear regression where the regularization parameter is
    automatically tuned (#2030).

  * Fix defeatist search for spill tree traversals (#2566, #1269).

  * Fix incremental training of logistic regression models (#2560).

  * Change default configuration of `BUILD_PYTHON_BINDINGS` to `OFF` (#2575).

### mlpack 3.3.2
###### 2020-06-18
  * Added Noisy DQN to q_networks (#2446).

  * Add Go bindings (#1884).

  * Added Dueling DQN to q_networks, Noisy linear layer to ann/layer
    and Empty loss to ann/loss_functions (#2414).

  * Storing and adding accessor method for action in q_learning (#2413).

  * Added accessor methods for ANN layers (#2321).

  * Addition of `Elliot` activation function (#2268).

  * Add adaptive max pooling and adaptive mean pooling layers (#2195).

  * Add parameter to avoid shuffling of data in preprocess_split (#2293).

  * Add `MatType` parameter to `LSHSearch`, allowing sparse matrices to be used
    for search (#2395).

  * Documentation fixes to resolve Doxygen warnings and issues (#2400).

  * Add Load and Save of Sparse Matrix (#2344).

  * Add Intersection over Union (IoU) metric for bounding boxes (#2402).

  * Add Non Maximal Supression (NMS) metric for bounding boxes (#2410).

  * Fix `no_intercept` and probability computation for linear SVM bindings
    (#2419).

  * Fix incorrect neighbors for `k > 1` searches in `approx_kfn` binding, for
    the `QDAFN` algorithm (#2448).

  * Fix serialization of kernels with state for FastMKS (#2452).

  * Add `RBF` layer in ann module to make `RBFN` architecture (#2261).

### mlpack 3.3.1
###### 2020-04-29
  * Minor Julia and Python documentation fixes (#2373).

  * Updated terminal state and fixed bugs for Pendulum environment (#2354,
    #2369).

  * Added `EliSH` activation function (#2323).

  * Add L1 Loss function (#2203).

  * Pass CMAKE_CXX_FLAGS (compilation options) correctly to Python build
    (#2367).

  * Expose ensmallen Callbacks for sparseautoencoder (#2198).

  * Bugfix for LARS class causing invalid read (#2374).

  * Add serialization support from Julia; use `mlpack.serialize()` and
    `mlpack.deserialize()` to save and load from `IOBuffer`s.

### mlpack 3.3.0
###### 2020-04-07
  * Added `Normal Distribution` to `ann/dists` (#2382).

  * Templated return type of `Forward function` of loss functions (#2339).

  * Added `R2 Score` regression metric (#2323).

  * Added `poisson negative log likelihood` loss function (#2196).

  * Added `huber` loss function (#2199).

  * Added `mean squared logarithmic error` loss function for neural networks
    (#2210).

  * Added `mean bias loss function` for neural networks (#2210).

  * The DecisionStump class has been marked deprecated; use the `DecisionTree`
    class with `NoRecursion=true` or use `ID3DecisionStump` instead (#2099).

  * Added `probabilities_file` parameter to get the probabilities matrix of
    AdaBoost classifier (#2050).

  * Fix STB header search paths (#2104).

  * Add `DISABLE_DOWNLOADS` CMake configuration option (#2104).

  * Add padding layer in TransposedConvolutionLayer (#2082).

  * Fix pkgconfig generation on non-Linux systems (#2101).

  * Use log-space to represent HMM initial state and transition probabilities
    (#2081).

  * Add functions to access parameters of `Convolution` and `AtrousConvolution`
    layers (#1985).

  * Add Compute Error function in lars regression and changing Train function to
    return computed error (#2139).

  * Add Julia bindings (#1949).  Build settings can be controlled with the
    `BUILD_JULIA_BINDINGS=(ON/OFF)` and `JULIA_EXECUTABLE=/path/to/julia` CMake
    parameters.

  * CMake fix for finding STB include directory (#2145).

  * Add bindings for loading and saving images (#2019); `mlpack_image_converter`
    from the command-line, `mlpack.image_converter()` from Python.

  * Add normalization support for CF binding (#2136).

  * Add Mish activation function (#2158).

  * Update `init_rules` in AMF to allow users to merge two initialization
    rules (#2151).

  * Add GELU activation function (#2183).

  * Better error handling of eigendecompositions and Cholesky decompositions
    (#2088, #1840).

  * Add LiSHT activation function (#2182).

  * Add Valid and Same Padding for Transposed Convolution layer (#2163).

  * Add CELU activation function (#2191)

  * Add Log-Hyperbolic-Cosine Loss function (#2207).

  * Change neural network types to avoid unnecessary use of rvalue references
    (#2259).

  * Bump minimum Boost version to 1.58 (#2305).

  * Refactor STB support so `HAS_STB` macro is not needed when compiling against
    mlpack (#2312).

  * Add Hard Shrink Activation Function (#2186).

  * Add Soft Shrink Activation Function (#2174).

  * Add Hinge Embedding Loss Function (#2229).

  * Add Cosine Embedding Loss Function (#2209).

  * Add Margin Ranking Loss Function (#2264).

  * Bugfix for incorrect parameter vector sizes in logistic regression and
    softmax regression (#2359).

### mlpack 3.2.2
###### 2019-11-26
  * Add `valid` and `same` padding option in `Convolution` and `Atrous
    Convolution` layer (#1988).

  * Add Model() to the FFN class to access individual layers (#2043).

  * Update documentation for pip and conda installation packages (#2044).

  * Add bindings for linear SVM (#1935); `mlpack_linear_svm` from the
    command-line, `linear_svm()` from Python.

  * Add support to return the layer name as `std::string` (#1987).

  * Speed and memory improvements for the Transposed Convolution layer (#1493).

  * Fix Windows Python build configuration (#1885).

  * Validate md5 of STB library after download (#2087).

  * Add `__version__` to `__init__.py` (#2092).

  * Correctly handle RNN sequences that are shorter than the value of rho (#2102).

### mlpack 3.2.1
###### 2019-10-01
  * Enforce CMake version check for ensmallen (#2032).

  * Fix CMake check for Armadillo version (#2029).

  * Better handling of when STB is not installed (#2033).

  * Fix Naive Bayes classifier computations in high dimensions (#2022).

### mlpack 3.2.0
###### 2019-09-25
  * Fix some potential infinity errors in Naive Bayes Classifier (#2022).

  * Fix occasionally-failing RADICAL test (#1924).

  * Fix gcc 9 OpenMP compilation issue (#1970).

  * Added support for loading and saving of images (#1903).

  * Add Multiple Pole Balancing Environment (#1901, #1951).

  * Added functionality for scaling of data (#1876); see the command-line
    binding `mlpack_preprocess_scale` or Python binding `preprocess_scale()`.

  * Add new parameter `maximum_depth` to decision tree and random forest
    bindings (#1916).

  * Fix prediction output of softmax regression when test set accuracy is
    calculated (#1922).

  * Pendulum environment now checks for termination. All RL environments now
    have an option to terminate after a set number of time steps (no limit
    by default) (#1941).

  * Add support for probabilistic KDE (kernel density estimation) error bounds
    when using the Gaussian kernel (#1934).

  * Fix negative distances for cover tree computation (#1979).

  * Fix cover tree building when all pairwise distances are 0 (#1986).

  * Improve KDE pruning by reclaiming not used error tolerance (#1954, #1984).

  * Optimizations for sparse matrix accesses in z-score normalization for CF
    (#1989).

  * Add `kmeans_max_iterations` option to GMM training binding `gmm_train_main`.

  * Bump minimum Armadillo version to 8.400.0 due to ensmallen dependency
    requirement (#2015).

### mlpack 3.1.1
###### 2019-05-26
  * Fix random forest bug for numerical-only data (#1887).

  * Significant speedups for random forest (#1887).

  * Random forest now has `minimum_gain_split` and `subspace_dim` parameters
    (#1887).

  * Decision tree parameter `print_training_error` deprecated in favor of
    `print_training_accuracy`.

  * `output` option changed to `predictions` for adaboost and perceptron
    binding. Old options are now deprecated and will be preserved until mlpack
    4.0.0 (#1882).

  * Concatenated ReLU layer (#1843).

  * Accelerate NormalizeLabels function using hashing instead of linear search
    (see `src/mlpack/core/data/normalize_labels_impl.hpp`) (#1780).

  * Add `ConfusionMatrix()` function for checking performance of classifiers
    (#1798).

  * Install ensmallen headers when it is downloaded during build (#1900).

### mlpack 3.1.0
###### 2019-04-25
  * Add DiagonalGaussianDistribution and DiagonalGMM classes to speed up the
    diagonal covariance computation and deprecate DiagonalConstraint (#1666).

  * Add kernel density estimation (KDE) implementation with bindings to other
    languages (#1301).

  * Where relevant, all models with a `Train()` method now return a `double`
    value representing the goodness of fit (i.e. final objective value, error,
    etc.) (#1678).

  * Add implementation for linear support vector machine (see
    `src/mlpack/methods/linear_svm`).

  * Change DBSCAN to use PointSelectionPolicy and add OrderedPointSelection (#1625).

  * Residual block support (#1594).

  * Bidirectional RNN (#1626).

  * Dice loss layer (#1674, #1714) and hard sigmoid layer (#1776).

  * `output` option changed to `predictions` and `output_probabilities` to
    `probabilities` for Naive Bayes binding (`mlpack_nbc`/`nbc()`).  Old options
    are now deprecated and will be preserved until mlpack 4.0.0 (#1616).

  * Add support for Diagonal GMMs to HMM code (#1658, #1666).  This can provide
    large speedup when a diagonal GMM is acceptable as an emission probability
    distribution.

  * Python binding improvements: check parameter type (#1717), avoid copying
    Pandas dataframes (#1711), handle Pandas Series objects (#1700).

### mlpack 3.0.4
###### 2018-11-13
  * Bump minimum CMake version to 3.3.2.

  * CMake fixes for Ninja generator by Marc Espie.

### mlpack 3.0.3
###### 2018-07-27
  * Fix Visual Studio compilation issue (#1443).

  * Allow running local_coordinate_coding binding with no initial_dictionary
    parameter when input_model is not specified (#1457).

  * Make use of OpenMP optional via the CMake 'USE_OPENMP' configuration
    variable (#1474).

  * Accelerate FNN training by 20-30% by avoiding redundant calculations
    (#1467).

  * Fix math::RandomSeed() usage in tests (#1462, #1440).

  * Generate better Python setup.py with documentation (#1460).

### mlpack 3.0.2
###### 2018-06-08
  * Documentation generation fixes for Python bindings (#1421).

  * Fix build error for man pages if command-line bindings are not being built
    (#1424).

  * Add 'shuffle' parameter and Shuffle() method to KFoldCV (#1412).  This will
    shuffle the data when the object is constructed, or when Shuffle() is
    called.

  * Added neural network layers: AtrousConvolution (#1390), Embedding (#1401),
    and LayerNorm (layer normalization) (#1389).

  * Add Pendulum environment for reinforcement learning (#1388) and update
    Mountain Car environment (#1394).

### mlpack 3.0.1
###### 2018-05-10
  * Fix intermittently failing tests (#1387).

  * Add big-batch SGD (BBSGD) optimizer in
    src/mlpack/core/optimizers/bigbatch_sgd/ (#1131).

  * Fix simple compiler warnings (#1380, #1373).

  * Simplify NeighborSearch constructor and Train() overloads (#1378).

  * Add warning for OpenMP setting differences (#1358/#1382).  When mlpack is
    compiled with OpenMP but another application is not (or vice versa), a
    compilation warning will now be issued.

  * Restructured loss functions in src/mlpack/methods/ann/ (#1365).

  * Add environments for reinforcement learning tests (#1368, #1370, #1329).

  * Allow single outputs for multiple timestep inputs for recurrent neural
    networks (#1348).

  * Add He and LeCun normal initializations for neural networks (#1342).
    Neural networks: add He and LeCun normal initializations (#1342), add FReLU
    and SELU activation functions (#1346, #1341), add alpha-dropout (#1349).

### mlpack 3.0.0
###### 2018-03-30
  * Speed and memory improvements for DBSCAN.  --single_mode can now be used for
    situations where previously RAM usage was too high.

  * Bump minimum required version of Armadillo to 6.500.0.

  * Add automatically generated Python bindings.  These have the same interface
    as the command-line programs.

  * Add deep learning infrastructure in src/mlpack/methods/ann/.

  * Add reinforcement learning infrastructure in
    src/mlpack/methods/reinforcement_learning/.

  * Add optimizers: AdaGrad, CMAES, CNE, FrankeWolfe, GradientDescent,
    GridSearch, IQN, Katyusha, LineSearch, ParallelSGD, SARAH, SCD, SGDR,
    SMORMS3, SPALeRA, SVRG.

  * Add hyperparameter tuning infrastructure and cross-validation infrastructure
    in src/mlpack/core/cv/ and src/mlpack/core/hpt/.

  * Fix bug in mean shift.

  * Add random forests (see src/mlpack/methods/random_forest).

  * Numerous other bugfixes and testing improvements.

  * Add randomized Krylov SVD and Block Krylov SVD.

### mlpack 2.2.5
###### 2017-08-25
  * Compilation fix for some systems (#1082).

  * Fix PARAM_INT_OUT() (#1100).

### mlpack 2.2.4
###### 2017-07-18
  * Speed and memory improvements for DBSCAN. --single_mode can now be used for
    situations where previously RAM usage was too high.

  * Fix bug in CF causing incorrect recommendations.

### mlpack 2.2.3
###### 2017-05-24
  * Bug fix for --predictions_file in mlpack_decision_tree program.

### mlpack 2.2.2
###### 2017-05-04
  * Install backwards-compatibility mlpack_allknn and mlpack_allkfn programs;
    note they are deprecated and will be removed in mlpack 3.0.0 (#992).

  * Fix RStarTree bug that surfaced on OS X only (#964).

  * Small fixes for MiniBatchSGD and SGD and tests.

### mlpack 2.2.1
###### 2017-04-13
  * Compilation fix for mlpack_nca and mlpack_test on older Armadillo versions
    (#984).

### mlpack 2.2.0
###### 2017-03-21
  * Bugfix for mlpack_knn program (#816).

  * Add decision tree implementation in methods/decision_tree/.  This is very
    similar to a C4.5 tree learner.

  * Add DBSCAN implementation in methods/dbscan/.

  * Add support for multidimensional discrete distributions (#810, #830).

  * Better output for Log::Debug/Log::Info/Log::Warn/Log::Fatal for Armadillo
    objects (#895, #928).

  * Refactor categorical CSV loading with boost::spirit for faster loading
    (#681).

### mlpack 2.1.1
###### 2016-12-22
  * HMMs now use random initialization; this should fix some convergence issues
    (#828).

  * HMMs now initialize emissions according to the distribution of observations
    (#833).

  * Minor fix for formatted output (#814).

  * Fix DecisionStump to properly work with any input type.

### mlpack 2.1.0
###### 2016-10-31
  * Fixed CoverTree to properly handle single-point datasets.

  * Fixed a bug in CosineTree (and thus QUIC-SVD) that caused split failures for
    some datasets (#717).

  * Added mlpack_preprocess_describe program, which can be used to print
    statistics on a given dataset (#742).

  * Fix prioritized recursion for k-furthest-neighbor search (mlpack_kfn and the
    KFN class), leading to orders-of-magnitude speedups in some cases.

  * Bump minimum required version of Armadillo to 4.200.0.

  * Added simple Gradient Descent optimizer, found in
    src/mlpack/core/optimizers/gradient_descent/ (#792).

  * Added approximate furthest neighbor search algorithms QDAFN and
    DrusillaSelect in src/mlpack/methods/approx_kfn/, with command-line program
    mlpack_approx_kfn.

### mlpack 2.0.3
###### 2016-07-21
  * Added multiprobe LSH (#691).  The parameter 'T' to LSHSearch::Search() can
    now be used to control the number of extra bins that are probed, as can the
    -T (--num_probes) option to mlpack_lsh.

  * Added the Hilbert R tree to src/mlpack/core/tree/rectangle_tree/ (#664).  It
    can be used as the typedef HilbertRTree, and it is now an option in the
    mlpack_knn, mlpack_kfn, mlpack_range_search, and mlpack_krann command-line
    programs.

  * Added the mlpack_preprocess_split and mlpack_preprocess_binarize programs,
    which can be used for preprocessing code (#650, #666).

  * Added OpenMP support to LSHSearch and mlpack_lsh (#700).

### mlpack 2.0.2
###### 2016-06-20
  * Added the function LSHSearch::Projections(), which returns an arma::cube
    with each projection table in a slice (#663).  Instead of Projection(i), you
    should now use Projections().slice(i).

  * A new constructor has been added to LSHSearch that creates objects using
    projection tables provided in an arma::cube (#663).

  * Handle zero-variance dimensions in DET (#515).

  * Add MiniBatchSGD optimizer (src/mlpack/core/optimizers/minibatch_sgd/) and
    allow its use in mlpack_logistic_regression and mlpack_nca programs.

  * Add better backtrace support from Grzegorz Krajewski for Log::Fatal messages
    when compiled with debugging and profiling symbols.  This requires libbfd
    and libdl to be present during compilation.

  * CosineTree test fix from Mikhail Lozhnikov (#358).

  * Fixed HMM initial state estimation (#600).

  * Changed versioning macros __MLPACK_VERSION_MAJOR, __MLPACK_VERSION_MINOR,
    and __MLPACK_VERSION_PATCH to MLPACK_VERSION_MAJOR, MLPACK_VERSION_MINOR,
    and MLPACK_VERSION_PATCH.  The old names will remain in place until
    mlpack 3.0.0.

  * Renamed mlpack_allknn, mlpack_allkfn, and mlpack_allkrann to mlpack_knn,
    mlpack_kfn, and mlpack_krann.  The mlpack_allknn, mlpack_allkfn, and
    mlpack_allkrann programs will remain as copies until mlpack 3.0.0.

  * Add --random_initialization option to mlpack_hmm_train, for use when no
    labels are provided.

  * Add --kill_empty_clusters option to mlpack_kmeans and KillEmptyClusters
    policy for the KMeans class (#595, #596).

### mlpack 2.0.1
###### 2016-02-04
  * Fix CMake to properly detect when MKL is being used with Armadillo.

  * Minor parameter handling fixes to mlpack_logistic_regression (#504, #505).

  * Properly install arma_config.hpp.

  * Memory handling fixes for Hoeffding tree code.

  * Add functions that allow changing training-time parameters to HoeffdingTree
    class.

  * Fix infinite loop in sparse coding test.

  * Documentation spelling fixes (#501).

  * Properly handle covariances for Gaussians with large condition number
    (#496), preventing GMMs from filling with NaNs during training (and also
    HMMs that use GMMs).

  * CMake fixes for finding LAPACK and BLAS as Armadillo dependencies when ATLAS
    is used.

  * CMake fix for projects using mlpack's CMake configuration from elsewhere
    (#512).

### mlpack 2.0.0
###### 2015-12-24
  * Removed overclustering support from k-means because it is not well-tested,
    may be buggy, and is (I think) unused.  If this was support you were using,
    open a bug or get in touch with us; it would not be hard for us to
    reimplement it.

  * Refactored KMeans to allow different types of Lloyd iterations.

  * Added implementations of k-means: Elkan's algorithm, Hamerly's algorithm,
    Pelleg-Moore's algorithm, and the DTNN (dual-tree nearest neighbor)
    algorithm.

  * Significant acceleration of LRSDP via the use of accu(a % b) instead of
    trace(a * b).

  * Added MatrixCompletion class (matrix_completion), which performs nuclear
    norm minimization to fill unknown values of an input matrix.

  * No more dependence on Boost.Random; now we use C++11 STL random support.

  * Add softmax regression, contributed by Siddharth Agrawal and QiaoAn Chen.

  * Changed NeighborSearch, RangeSearch, FastMKS, LSH, and RASearch API; these
    classes now take the query sets in the Search() method, instead of in the
    constructor.

  * Use OpenMP, if available.  For now OpenMP support is only available in the
    DET training code.

  * Add support for predicting new test point values to LARS and the
    command-line 'lars' program.

  * Add serialization support for Perceptron and LogisticRegression.

  * Refactor SoftmaxRegression to predict into an arma::Row<size_t> object, and
    add a softmax_regression program.

  * Refactor LSH to allow loading and saving of models.

  * ToString() is removed entirely (#487).

  * Add --input_model_file and --output_model_file options to appropriate
    machine learning algorithms.

  * Rename all executables to start with an "mlpack" prefix (#229).

  * Add HoeffdingTree and mlpack_hoeffding_tree, an implementation of the
    streaming decision tree methodology from Domingos and Hulten in 2000.

### mlpack 1.0.12
###### 2015-01-07
  * Switch to 3-clause BSD license (from LGPL).

### mlpack 1.0.11
###### 2014-12-11
  * Proper handling of dimension calculation in PCA.

  * Load parameter vectors properly for LinearRegression models.

  * Linker fixes for AugLagrangian specializations under Visual Studio.

  * Add support for observation weights to LinearRegression.

  * MahalanobisDistance<> now takes the root of the distance by default and
    therefore satisfies the triangle inequality (TakeRoot now defaults to true).

  * Better handling of optional Armadillo HDF5 dependency.

  * Fixes for numerous intermittent test failures.

  * math::RandomSeed() now sets the random seed for recent (>=3.930) Armadillo
    versions.

  * Handle Newton method convergence better for
    SparseCoding::OptimizeDictionary() and make maximum iterations a parameter.

  * Known bug: CosineTree construction may fail in some cases on i386 systems
    (#358).

### mlpack 1.0.10
###### 2014-08-29
  * Bugfix for NeighborSearch regression which caused very slow allknn/allkfn.
    Speeds are now restored to approximately 1.0.8 speeds, with significant
    improvement for the cover tree (#347).

  * Detect dependencies correctly when ARMA_USE_WRAPPER is not being defined
    (i.e., libarmadillo.so does not exist).

  * Bugfix for compilation under Visual Studio (#348).

### mlpack 1.0.9
###### 2014-07-28
  * GMM initialization is now safer and provides a working GMM when constructed
    with only the dimensionality and number of Gaussians (#301).

  * Check for division by 0 in Forward-Backward Algorithm in HMMs (#301).

  * Fix MaxVarianceNewCluster (used when re-initializing clusters for k-means)
    (#301).

  * Fixed implementation of Viterbi algorithm in HMM::Predict() (#303).

  * Significant speedups for dual-tree algorithms using the cover tree (#235,
    #314) including a faster implementation of FastMKS.

  * Fix for LRSDP optimizer so that it compiles and can be used (#312).

  * CF (collaborative filtering) now expects users and items to be zero-indexed,
    not one-indexed (#311).

  * CF::GetRecommendations() API change: now requires the number of
    recommendations as the first parameter.  The number of users in the local
    neighborhood should be specified with CF::NumUsersForSimilarity().

  * Removed incorrect PeriodicHRectBound (#58).

  * Refactor LRSDP into LRSDP class and standalone function to be optimized
    (#305).

  * Fix for centering in kernel PCA (#337).

  * Added simulated annealing (SA) optimizer, contributed by Zhihao Lou.

  * HMMs now support initial state probabilities; these can be set in the
    constructor, trained, or set manually with HMM::Initial() (#302).

  * Added Nyström method for kernel matrix approximation by Marcus Edel.

  * Kernel PCA now supports using Nyström method for approximation.

  * Ball trees now work with dual-tree algorithms, via the BallBound<> bound
    structure (#307); fixed by Yash Vadalia.

  * The NMF class is now AMF<>, and supports far more types of factorizations,
    by Sumedh Ghaisas.

  * A QUIC-SVD implementation has returned, written by Siddharth Agrawal and
    based on older code from Mudit Gupta.

  * Added perceptron and decision stump by Udit Saxena (these are weak learners
    for an eventual AdaBoost class).

  * Sparse autoencoder added by Siddharth Agrawal.

### mlpack 1.0.8
###### 2014-01-06
  * Memory leak in NeighborSearch index-mapping code fixed (#298).

  * GMMs can be trained using the existing model as a starting point by
    specifying an additional boolean parameter to GMM::Estimate() (#296).

  * Logistic regression implementation added in methods/logistic_regression (see
    also #293).

  * L-BFGS optimizer now returns its function via Function().

  * Version information is now obtainable via mlpack::util::GetVersion() or the
    __MLPACK_VERSION_MAJOR, __MLPACK_VERSION_MINOR, and  __MLPACK_VERSION_PATCH
    macros (#297).

  * Fix typos in allkfn and allkrann output.

### mlpack 1.0.7
###### 2013-10-04
  * Cover tree support for range search (range_search), rank-approximate nearest
    neighbors (allkrann), minimum spanning tree calculation (emst), and FastMKS
    (fastmks).

  * Dual-tree FastMKS implementation added and tested.

  * Added collaborative filtering package (cf) that can provide recommendations
    when given users and items.

  * Fix for correctness of Kernel PCA (kernel_pca) (#270).

  * Speedups for PCA and Kernel PCA (#198).

  * Fix for correctness of Neighborhood Components Analysis (NCA) (#279).

  * Minor speedups for dual-tree algorithms.

  * Fix for Naive Bayes Classifier (nbc) (#269).

  * Added a ridge regression option to LinearRegression (linear_regression)
    (#286).

  * Gaussian Mixture Models (gmm::GMM<>) now support arbitrary covariance matrix
    constraints (#283).

  * MVU (mvu) removed because it is known to not work (#183).

  * Minor updates and fixes for kernels (in mlpack::kernel).

### mlpack 1.0.6
###### 2013-06-13
  * Minor bugfix so that FastMKS gets built.

### mlpack 1.0.5
###### 2013-05-01
  * Speedups of cover tree traversers (#235).

  * Addition of rank-approximate nearest neighbors (RANN), found in
    src/mlpack/methods/rann/.

  * Addition of fast exact max-kernel search (FastMKS), found in
    src/mlpack/methods/fastmks/.

  * Fix for EM covariance estimation; this should improve GMM training time.

  * More parameters for GMM estimation.

  * Force GMM and GaussianDistribution covariance matrices to be positive
    definite, so that training converges much more often.

  * Add parameter for the tolerance of the Baum-Welch algorithm for HMM
    training.

  * Fix for compilation with clang compiler.

  * Fix for k-furthest-neighbor-search.

### mlpack 1.0.4
###### 2013-02-08
  * Force minimum Armadillo version to 2.4.2.

  * Better output of class types to streams; a class with a ToString() method
    implemented can be sent to a stream with operator<<.

  * Change return type of GMM::Estimate() to double (#257).

  * Style fixes for k-means and RADICAL.

  * Handle size_t support correctly with Armadillo 3.6.2 (#258).

  * Add locality-sensitive hashing (LSH), found in src/mlpack/methods/lsh/.

  * Better tests for SGD (stochastic gradient descent) and NCA (neighborhood
    components analysis).

### mlpack 1.0.3
###### 2012-09-16

  * Remove internal sparse matrix support because Armadillo 3.4.0 now includes
    it.  When using Armadillo versions older than 3.4.0, sparse matrix support
    is not available.

  * NCA (neighborhood components analysis) now support an arbitrary optimizer
    (#245), including stochastic gradient descent (#249).

### mlpack 1.0.2
###### 2012-08-15
  * Added density estimation trees, found in src/mlpack/methods/det/.

  * Added non-negative matrix factorization, found in src/mlpack/methods/nmf/.

  * Added experimental cover tree implementation, found in
    src/mlpack/core/tree/cover_tree/ (#157).

  * Better reporting of boost::program_options errors (#225).

  * Fix for timers on Windows (#212, #211).

  * Fix for allknn and allkfn output (#204).

  * Sparse coding dictionary initialization is now a template parameter (#220).

### mlpack 1.0.1
###### 2012-03-03
  * Added kernel principal components analysis (kernel PCA), found in
    src/mlpack/methods/kernel_pca/ (#74).

  * Fix for Lovasz-Theta AugLagrangian tests (#182).

  * Fixes for allknn output (#185, #186).

  * Added range search executable (#192).

  * Adapted citations in documentation to BibTeX; no citations in -h output
    (#195).

  * Stop use of 'const char*' and prefer 'std::string' (#176).

  * Support seeds for random numbers (#177).

### mlpack 1.0.0
###### 2011-12-17
  * Initial release.  See any resolved tickets numbered less than #196 or
    execute this query:
    http://www.mlpack.org/trac/query?status=closed&milestone=mlpack+1.0.0
<h2 align="center">
  <a href="http://mlpack.org"><img
src="https://cdn.rawgit.com/mlpack/mlpack.org/e7d36ed8/mlpack-black.svg" style="background-color:rgba(0,0,0,0);" height=230 alt="mlpack: a fast, flexible machine learning library"></a>
  <br>a fast, flexible machine learning library<br>
</h2>

<h5 align="center">
  <a href="https://mlpack.org">Home</a> |
  <a href="https://www.mlpack.org/docs.html">Documentation</a> |
  <a href="https://www.mlpack.org/doc/mlpack-git/doxygen/index.html">Doxygen</a> |
  <a href="https://www.mlpack.org/community.html">Community</a> |
  <a href="https://www.mlpack.org/questions.html">Help</a> |
  <a href="https://webchat.freenode.net/?channels=mlpack">IRC Chat</a>
</h5>

<p align="center">
  <a href="https://dev.azure.com/mlpack/mlpack/_build?definitionId=1"><img alt="Azure DevOps builds (job)" src="https://img.shields.io/azure-devops/build/mlpack/84320e87-76e3-4b6e-8b6e-3adaf6b36eed/1/master?job=Linux&label=Linux%20Build&style=flat-square"></a>
  <a href="https://opensource.org/licenses/BSD-3-Clause"><img src="https://img.shields.io/badge/License-BSD%203--Clause-blue.svg?style=flat-square" alt="License"></a>
  <a href="http://numfocus.org/donate-to-mlpack"><img src="https://img.shields.io/badge/sponsored%20by-NumFOCUS-orange.svg?style=flat-square&colorA=E1523D&colorB=007D8A" alt="NumFOCUS"></a>
</p>

<p align="center">
  <em>
    Download:
    <a href="https://www.mlpack.org/files/mlpack-3.4.2.tar.gz">current stable version (3.4.2)</a>
  </em>
</p>

**mlpack** is an intuitive, fast, and flexible C++ machine learning library with
bindings to other languages.  It is meant to be a machine learning analog to
LAPACK, and aims to implement a wide array of machine learning methods and
functions as a "swiss army knife" for machine learning researchers.  In addition
to its powerful C++ interface, mlpack also provides command-line programs,
Python bindings, Julia bindings, Go bindings and R bindings.

[//]: # (numfocus-fiscal-sponsor-attribution)

mlpack uses an [open governance model](./GOVERNANCE.md) and is fiscally
sponsored by [NumFOCUS](https://numfocus.org/).  Consider making a
[tax-deductible donation](https://numfocus.org/donate-to-mlpack) to help the
project pay for developer time, professional services, travel, workshops, and a
variety of other needs.

<div align="center">
  <a href="https://numfocus.org/">
    <img height="60px"
         src="https://raw.githubusercontent.com/numfocus/templates/master/images/numfocus-logo.png"
         align="center">
  </a>
</div>
<br>

### 0. Contents

  1. [Introduction](#1-introduction)
  2. [Citation details](#2-citation-details)
  3. [Dependencies](#3-dependencies)
  4. [Building mlpack from source](#4-building-mlpack-from-source)
  5. [Running mlpack programs](#5-running-mlpack-programs)
  6. [Using mlpack from Python](#6-using-mlpack-from-python)
  7. [Further documentation](#7-further-documentation)
  8. [Bug reporting](#8-bug-reporting)

###  1. Introduction

The mlpack website can be found at https://www.mlpack.org and it contains
numerous tutorials and extensive documentation.  This README serves as a guide
for what mlpack is, how to install it, how to run it, and where to find more
documentation. The website should be consulted for further information:

  - [mlpack homepage](https://www.mlpack.org/)
  - [mlpack documentation](https://www.mlpack.org/docs.html)
  - [Tutorials](https://www.mlpack.org/doc/mlpack-git/doxygen/tutorials.html)
  - [Development Site (Github)](https://www.github.com/mlpack/mlpack/)
  - [API documentation (Doxygen)](https://www.mlpack.org/doc/mlpack-git/doxygen/index.html)

### 2. Citation details

If you use mlpack in your research or software, please cite mlpack using the
citation below (given in BibTeX format):

    @article{mlpack2018,
        title     = {mlpack 3: a fast, flexible machine learning library},
        author    = {Curtin, Ryan R. and Edel, Marcus and Lozhnikov, Mikhail and
                     Mentekidis, Yannis and Ghaisas, Sumedh and Zhang,
                     Shangtong},
        journal   = {Journal of Open Source Software},
        volume    = {3},
        issue     = {26},
        pages     = {726},
        year      = {2018},
        doi       = {10.21105/joss.00726},
        url       = {https://doi.org/10.21105/joss.00726}
    }

Citations are beneficial for the growth and improvement of mlpack.

### 3. Dependencies

mlpack has the following dependencies:

      Armadillo      >= 9.800
      Boost (math_c99, spirit) >= 1.58.0
      CMake          >= 3.6
      ensmallen      >= 2.10.0
      cereal         >= 1.1.2

All of those should be available in your distribution's package manager.  If
not, you will have to compile each of them by hand.  See the documentation for
each of those packages for more information.

If you would like to use or build the mlpack Python bindings, make sure that the
following Python packages are installed:

      setuptools
      cython >= 0.24
      numpy
      pandas >= 0.15.0

If you would like to build the Julia bindings, make sure that Julia >= 1.3.0 is
installed.

If you would like to build the Go bindings, make sure that Go >= 1.11.0 is
installed with this package:

     Gonum

If you would like to build the R bindings, make sure that R >= 4.0 is
installed with these R packages.

     Rcpp >= 0.12.12
     RcppArmadillo >= 0.8.400.0
     RcppEnsmallen >= 0.2.10.0
     BH >= 1.58
     roxygen2

If the STB library headers are available, image loading support will be
compiled.

If you are compiling Armadillo by hand, ensure that LAPACK and BLAS are enabled.

### 4. Building mlpack from source

This document discusses how to build mlpack from source. These build directions
will work for any Linux-like shell environment (for example Ubuntu, macOS,
FreeBSD etc). However, mlpack is in the repositories of many Linux distributions
and so it may be easier to use the package manager for your system.  For example,
on Ubuntu, you can install the mlpack library and command-line executables (e.g.
mlpack_pca, mlpack_kmeans etc.) with the following command:

    $ sudo apt-get install libmlpack-dev mlpack-bin

On Fedora or Red Hat (EPEL):

    $ sudo dnf install mlpack-devel mlpack-bin

*Note*: Older Ubuntu versions may not have the most recent version of mlpack
available---for instance, at the time of this writing, Ubuntu 16.04 only has
mlpack 3.4.2 available.  Options include upgrading your Ubuntu version, finding
a PPA or other non-official sources, or installing with a manual build.

*Note*: If you are using RHEL7/CentOS 7, gcc 4.8 is too old to compile mlpack.
One option is to use `devtoolset-8`; see
[here](https://www.softwarecollections.org/en/scls/rhscl/devtoolset-8/) for more
information.

There are some useful pages to consult in addition to this section:

  - [Building mlpack From Source](https://www.mlpack.org/doc/mlpack-git/doxygen/build.html)
  - [Building mlpack From Source on Windows](https://www.mlpack.org/doc/mlpack-git/doxygen/build_windows.html)

mlpack uses CMake as a build system and allows several flexible build
configuration options. You can consult any of the CMake tutorials for
further documentation, but this tutorial should be enough to get mlpack built
and installed.

First, unpack the mlpack source and change into the unpacked directory.  Here we
use mlpack-x.y.z where x.y.z is the version.

    $ tar -xzf mlpack-x.y.z.tar.gz
    $ cd mlpack-x.y.z

Then, make a build directory.  The directory can have any name, but 'build' is
sufficient.

    $ mkdir build
    $ cd build

The next step is to run CMake to configure the project.  Running CMake is the
equivalent to running `./configure` with autotools. If you run CMake with no
options, it will configure the project to build with no debugging symbols and
no profiling information:

    $ cmake ../

Options can be specified to compile with debugging information and profiling information:

    $ cmake -D DEBUG=ON -D PROFILE=ON ../

Options are specified with the -D flag.  The allowed options include:

    DEBUG=(ON/OFF): compile with debugging symbols
    PROFILE=(ON/OFF): compile with profiling symbols
    ARMA_EXTRA_DEBUG=(ON/OFF): compile with extra Armadillo debugging symbols
    BOOST_ROOT=(/path/to/boost/): path to root of boost installation
    ARMADILLO_INCLUDE_DIR=(/path/to/armadillo/include/): path to Armadillo headers
    ARMADILLO_LIBRARY=(/path/to/armadillo/libarmadillo.so): Armadillo library
    BUILD_CLI_EXECUTABLES=(ON/OFF): whether or not to build command-line programs
    BUILD_PYTHON_BINDINGS=(ON/OFF): whether or not to build Python bindings
    PYTHON_EXECUTABLE=(/path/to/python_version): Path to specific Python executable
    PYTHON_INSTALL_PREFIX=(/path/to/python/): Path to root of Python installation
    BUILD_JULIA_BINDINGS=(ON/OFF): whether or not to build Julia bindings
    JULIA_EXECUTABLE=(/path/to/julia): Path to specific Julia executable
    BUILD_GO_BINDINGS=(ON/OFF): whether or not to build Go bindings
    GO_EXECUTABLE=(/path/to/go): Path to specific Go executable
    BUILD_GO_SHLIB=(ON/OFF): whether or not to build shared libraries required by Go bindings
    BUILD_R_BINDINGS=(ON/OFF): whether or not to build R bindings
    R_EXECUTABLE=(/path/to/R): Path to specific R executable
    BUILD_TESTS=(ON/OFF): whether or not to build tests
    BUILD_SHARED_LIBS=(ON/OFF): compile shared libraries and executables as
        opposed to static libraries
    DISABLE_DOWNLOADS=(ON/OFF): whether to disable all downloads during build
    ENSMALLEN_INCLUDE_DIR=(/path/to/ensmallen/include): path to include directory
       for ensmallen
    STB_IMAGE_INCLUDE_DIR=(/path/to/stb/include): path to include directory for
       STB image library
    USE_OPENMP=(ON/OFF): whether or not to use OpenMP if available
    BUILD_DOCS=(ON/OFF): build Doxygen documentation, if Doxygen is available
       (default ON)

For example, to build mlpack library and CLI bindings statically the following
command can be used:

    $ cmake  -D BUILD_SHARED_LIBS=OFF ../

Other tools can also be used to configure CMake, but those are not documented
here.  See [this section of the build guide](https://www.mlpack.org/doc/mlpack-git/doxygen/build.html#build_config)
for more details, including a full list of options, and their default values.

By default, command-line programs will be built, and if the Python dependencies
(Cython, setuptools, numpy, pandas) are available, then Python bindings will
also be built.  OpenMP will be used for parallelization when possible by
default.

Once CMake is configured, building the library is as simple as typing 'make'.
This will build all library components and bindings.

    $ make

If you do not want to build everything in the library, individual components
of the build can be specified:

    $ make mlpack_pca mlpack_knn mlpack_kfn

If you want to build the tests, just make the `mlpack_test` target, and use
`ctest` to run the tests:

    $ make mlpack_test
    $ ctest .

If the build fails and you cannot figure out why, register an account on Github
and submit an issue. The mlpack developers will quickly help you figure it out:

[mlpack on Github](https://www.github.com/mlpack/mlpack/)

Alternately, mlpack help can be found in IRC at `#mlpack` on chat.freenode.net.

If you wish to install mlpack to `/usr/local/include/mlpack/`, `/usr/local/lib/`,
and `/usr/local/bin/`, make sure you have root privileges (or write permissions
to those three directories), and simply type

    $ make install

You can now run the executables by name; you can link against mlpack with
    `-lmlpack`
and the mlpack headers are found in
    `/usr/local/include/mlpack/`
and if Python bindings were built, you can access them with the `mlpack`
package in Python.

If running the programs (i.e. `$ mlpack_knn -h`) gives an error of the form

    error while loading shared libraries: libmlpack.so.2: cannot open shared object file: No such file or directory

then be sure that the runtime linker is searching the directory where
`libmlpack.so` was installed (probably `/usr/local/lib/` unless you set it
manually).  One way to do this, on Linux, is to ensure that the
`LD_LIBRARY_PATH` environment variable has the directory that contains
`libmlpack.so`.  Using bash, this can be set easily:

    export LD_LIBRARY_PATH="/usr/local/lib/:$LD_LIBRARY_PATH"

(or whatever directory `libmlpack.so` is installed in.)

### 5. Running mlpack programs

After building mlpack, the executables will reside in `build/bin/`.  You can call
them from there, or you can install the library and (depending on system
settings) they should be added to your PATH and you can call them directly.  The
documentation below assumes the executables are in your PATH.

Consider the 'mlpack_knn' program, which finds the k nearest neighbors in a
reference dataset of all the points in a query set.  That is, we have a query
and a reference dataset. For each point in the query dataset, we wish to know
the k points in the reference dataset which are closest to the given query
point.

Alternately, if the query and reference datasets are the same, the problem can
be stated more simply: for each point in the dataset, we wish to know the k
nearest points to that point.

Each mlpack program has extensive help documentation which details what the
method does, what each of the parameters is, and how to use them:

```shell
$ mlpack_knn --help
```

Running `mlpack_knn` on one dataset (that is, the query and reference
datasets are the same) and finding the 5 nearest neighbors is very simple:

```shell
$ mlpack_knn -r dataset.csv -n neighbors_out.csv -d distances_out.csv -k 5 -v
```

The `-v (--verbose)` flag is optional; it gives informational output.  It is not
unique to `mlpack_knn` but is available in all mlpack programs.  Verbose
output also gives timing output at the end of the program, which can be very
useful.

### 6. Using mlpack from Python

If mlpack is installed to the system, then the mlpack Python bindings should be
automatically in your PYTHONPATH, and importing mlpack functionality into Python
should be very simple:

```python
>>> from mlpack import knn
```

Accessing help is easy:

```python
>>> help(knn)
```

The API is similar to the command-line programs.  So, running `knn()`
(k-nearest-neighbor search) on the numpy matrix `dataset` and finding the 5
nearest neighbors is very simple:

```python
>>> output = knn(reference=dataset, k=5, verbose=True)
```

This will store the output neighbors in `output['neighbors']` and the output
distances in `output['distances']`.  Other mlpack bindings function similarly,
and the input/output parameters exactly match those of the command-line
programs.

### 7. Further documentation

The documentation given here is only a fraction of the available documentation
for mlpack.  If doxygen is installed, you can type `make doc` to build the
documentation locally.  Alternately, up-to-date documentation is available for
older versions of mlpack:

  - [mlpack homepage](https://www.mlpack.org/)
  - [mlpack documentation](https://www.mlpack.org/docs.html)
  - [Tutorials](https://www.mlpack.org/doc/mlpack-git/doxygen/tutorials.html)
  - [Development Site (Github)](https://www.github.com/mlpack/mlpack/)
  - [API documentation (Doxygen)](https://www.mlpack.org/doc/mlpack-git/doxygen/index.html)

To learn about the development goals of mlpack in the short- and medium-term
future, see the [vision document](https://www.mlpack.org/papers/vision.pdf).

### 8. Bug reporting

   (see also [mlpack help](https://www.mlpack.org/questions.html))

If you find a bug in mlpack or have any problems, numerous routes are available
for help.

Github is used for bug tracking, and can be found at
https://github.com/mlpack/mlpack/.
It is easy to register an account and file a bug there, and the mlpack
development team will try to quickly resolve your issue.

In addition, mailing lists are available.  The mlpack discussion list is
available at

  [mlpack discussion list](http://lists.mlpack.org/mailman/listinfo/mlpack)

and the git commit list is available at

  [commit list](http://lists.mlpack.org/mailman/listinfo/mlpack-git)

Lastly, the IRC channel `#mlpack` on Freenode can be used to get help.
# mlpack governance structure

Revised Oct. 21st, 2019.

## Introduction

mlpack has grown much since its initial inception as a small project out of a
university research lab.  Now that there are over 150 contributors, it's
important that we have a clearly defined process for making our decisions and
organizing ourselves.

This document aims to clarify the governance of mlpack.  This is a living
document: it may change over time.  The process for making these changes is
detailed in the "Governance Changes" section.

## Code of Conduct

mlpack aims to be an open and welcoming environment, and as such, we have a
code of conduct that helps foster this environment.  See
[here](https://github.com/mlpack/mlpack/blob/master/CODE_OF_CONDUCT.md) for
more information.

## Teams & Roles

To keep overhead minimal, mlpack's teams and roles are simple: there is only
the [Committers](https://github.com/orgs/mlpack/teams/contributors) team, and
the [NumFOCUS leadership team](TODO:link).

Members of the Committers team have commit access to all mlpack repositories
and help guide the development directions and goals of mlpack.  Committers
should be familiar with the [contribution
process](https://github.com/mlpack/mlpack/blob/master/CONTRIBUTING.md) and
follow it when merging code and reviewing pull requests; this is important for
the continued stability and quality of mlpack's codebase.  Responsibilities and
activities of Committers team members can include:

 * Welcoming new members to the community: helping support users and point
   potential contributors in the correct direction.

 * Reviewing pull requests and approving them when they are ready.

 * Merging pull requests after they have been approved by others for merge.

 * Communicating and coordinating with contributors to help get code merged and
   improve the software.

 * Helping map out mlpack's development directions and processes.

 * Maintaining mlpack infrastructure (build systems, continuous integration,
   etc.).

Membership on the Committers team does not expire.  Contributors who have
repeatedly shown that their code quality is high, demonstrated adherence to the
code of conduct, and shown that they have a strong interest in the project can
be added to the Committers team using the organizational decision process in
the next section.

The NumFOCUS leadership team is a subset of the Committers team whose
additional responsibilities are to coordinate with NumFOCUS and maintain this
governance document.  Membership in the NumFOCUS leadership team is limited to
five people, and does not confer any special voting power or decision rights.

## Voting and Organizational Decisions

Historically, mlpack organizational decisions have not been controversial and
this has allowed efficient decision making.  Therefore, a vote on a proposal is
not required unless there is any explicit disagreement or concern with the
proposal.  The topics of a proposal might be:

 * Adding/removing a new member to/from the Committers team.

 * Participating in a program such as Google Summer of Code or Outreachy.

 * A change to some part of the mlpack infrastructure or contribution process.

 * Refactoring or change of an important public part of the API.

 * Use of funds for a particular project.

That list is not inclusive.  Introducing a proposal or idea can be done
informally in a public place, such as the mlpack mailing list or on Github as
an issue.  It's a good idea (but not mandatory) to make the proposal discussion
fully public so that people who are not on the Committers team can also comment
and provide opinions---after all, this is a community-led project so we should
be sure to include the *entire* community whenever possible.

If there is any disagreement or concern with the proposal, the person who
introduced the proposal should work to try and find a resolution or compromise
if possible.  If that is not possible, then the proposal can be brought to a
vote.

For a proposal to pass, a simple majority vote suffices.  Each Committer has
one equal vote, and they may choose to abstain from voting if they do prefer.
Since some Committers may be inactive or busy, it is not required for every
Committer to participate in every vote; instead, someone who has a proposal
should make a good-faith effort to post the proposal in a public location so
that interested and active Committers can respond.  Voting for any proposal
should be open for at least five days to allow sufficient time.

If a proposal passes despite votes against it, it is generally a good idea for
the Committer who introduced the proposal to spend some time considering and
understanding the arguments that were presented against the proposal, or if
appropriate, for the Committer to try and find an acceptable compromise or
alternate strategy that addresses the given feedback.

## Governance Changes

The NumFOCUS leadership team is responsible for this governance document, and
thus any changes to this document, NumFOCUS membership, or the NumFOCUS
leadership team must be approved by that team, also by a simple majority vote.
Because every member of the NumFOCUS leadership team should be an active
Committer, efforts should be made to collect votes from all five members.
Voting for any proposal should be open for at least five days to allow
sufficient time. If necessary, every active member of the NumFOCUS leadership
team can ask for an extension of the voting deadline.
# mlpack Code of Conduct

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and
expression, level of experience, education, socio-economic status, nationality,
personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission

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

This Code of Conduct applies within all project spaces, and it also applies when
an individual is representing the project or its community in public spaces.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at conduct@mlpack.org. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Reporting

If you believe someone is violating the code of conduct we ask that you report
it by emailing conduct@mlpack.org. All reports will be kept confidential. In
some cases we may determine that a public statement will need to be made. If
that's the case, the identities of all victims and reporters will remain
confidential unless those individuals instruct us otherwise.

If you are unsure whether the incident is a violation, or whether the space
where it happened is covered by this Code of Conduct, we encourage you to still
report it. We would much rather have a few extra reports where we decide to take
no action, rather than miss a report of an actual violation. We do not look
negatively on you if we find the incident is not a violation. And knowing about
incidents that are not violations, or happen outside our spaces, can also help
us to improve the Code of Conduct or the processes surrounding it.

In your report please include:

* Your contact info (so we can get in touch with you if we need to follow up)
* Names (real, nicknames, or pseudonyms) of any individuals involved. If there
  were other witnesses besides you, please try to include them as well.
* When and where the incident occurred. Please be as specific as possible.
* Your account of what occurred. If there is a publicly available record
  (e.g. a mailing list archive or a public IRC logger) please include a link.
* Any extra context you believe existed for the incident.
* If you believe this incident is ongoing.
* Any other information you believe we should have.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.4, available at
https://www.contributor-covenant.org/version/1/4/code-of-conduct.html, and
includes some aspects of the Drupal Code of Conduct.
# Contributing to mlpack

mlpack is a community-led project; that means that anyone is welcome to
contribute to mlpack and join the community!  If you would like to make
improvements to the library, add new features that are useful to you and others,
or have found a bug that you know how to fix, please submit a pull request!

If you would like to learn more about how to get started contributing, see the
[Community](http://www.mlpack.org/community.html) page, and if you are
interested in participating in Google Summer of Code, see
[mlpack and Google Summer of Code](http://www.mlpack.org/gsoc.html).

## Pull request process

Once a pull request is submitted, it must be approved by at least one member of
mlpack's Contributors team, to ensure that (if applicable):

 * the design meshes with the rest of mlpack
 * the style matches the
   [Style Guide](http://github.com/mlpack/mlpack/wiki/DesignGuidelines)
 * any new functionality is tested and working

The pull request can be merged as soon as it receives two approvals; 24 hours
after the first approval, mlpack-bot will provide a second approval.  This is to
leave time for anyone to comment on the PR before it is merged.

Members of the Contributors team are encouraged to review pull requests that
have already been reviewed, and pull request contributors are encouraged to seek
multiple reviews.  Reviews from anyone not on the Contributors team are always
appreciated and encouraged!

## Reviewing Pull Requests

All mlpack contributors who choose to review and provide feedback on pull
requests have a responsibility to both the project and the individual making
the contribution. 

Reviews and feedback should be
[helpful, insightful, and geared towards improving the contribution](
  https://www.youtube.com/watch?v=NNXk_WJzyMI).
If there are reasons why you feel the PR should not be merged, explain
what those are. Be open to having your mind changed. Be open to
working with the contributor to make the pull request better.

Please don't leave dismissive or disrespectful reviews!  It's not helpful for
anyone.

When reviewing a pull request, the primary goals are:

- For the codebase/project to improve
- For the person submitting the request to succeed

Even if a pull request does not get merged, the submitters should come away
from the experience feeling like their effort was not wasted or unappreciated.
Every pull request from a new contributor is an opportunity to grow the community. 

When changes are necessary, request them, do not demand them, and do not assume
that the contributor already knows how to do that. Be there to lend a helping
hand in case of need.

Since there can sometimes be a lot more pull requests being opened than
reviewed, we highly encourage everyone to review each others pull request
keeping in mind all the above mentioned points.

Let's welcome new contributors with ❤️.

## Pull Request Waiting Time 

mlpack is a community-driven project, so everyone only works on it in their
free time; this means it may take some time for them to review pull requests.
While gentle reminders are welcome, please be patient and avoid constantly
messaging contributors or tagging them on pull requests.

Typically small PRs will be reviewed within a handful of days; larger PRs might
take a few weeks for an initial review, and it may be a little bit longer in 
times of high activity.
# mlpack Tests

This directory contains code and data used to test all the algorithms and functions implemented in mlpack.

## Test Directories Structure

- *_test.cpp - methods tests
- main_tests/*_test.cpp - binding tests
- data - data needed to run the tests

## Add tests 

We have a rich test suite, consisting of almost 2000 tests (and still counting). It is suggested to add tests when:

- Adding new functionality.
- Fixing regressions and bugs.

## Building Tests

To build the test suite you can simply run `make mlpack_test`.

## To run Tests

We use `Catch2` to write our tests. To run all tests, you can simply run:

`./bin/mlpack_test`

To run all tests in a particular file you can run:

`./bin/mlpack_test "[testname]"`

where `testname` is the name of the test suite. For example to run all
collaborative filtering tests implemented in `cv_test.cpp` you can run:

`./bin/mlpack_test "[CVTest]"`

Now similarly you can run all the binding related tests using:

`./bin/mlpack_test "[BindingTests]"`

To run a single test, you can explicitly provide the name of the test, for example, to run
`BinaryClassificationMetricsTest` implemented in `cv_test.cpp` you can run the following:

`./bin/mlpack_test BinaryClassificationMetricsTest`

Catch2 provides many other features like filter, checkout the [Catch2 reference section](https://github.com/catchorg/Catch2/blob/devel/docs/Readme.md#top) - for more details.
The files in this directory are taken from MNMLSTC Core 1.1.0 in order to 
backport features from C++ 17 standard library:

 * C++17 STL algorithms such as std::any and std::basic_string_view. 
 * Dependencies files that are used to implement these features.

These files are licensed under the Apache 2.0 License, available in LICENSE.txt
in this directory.

If you want a copy of mlpack without a dependence on the Apache License or 
without the backported version then you will need to

 * Remove this entire directory.
 * Remove the line "std_backport" from src/mlpack/core/CMakeLists.txt.
 * Use the C++17 standard by modifying the mlpack/CMakeLists.txt.
---
name: Documentation issue
about: Use this template to report an issue you've found with the documentation.
title: ''
labels: 't: bug report, c: documentation, s: unanswered'
assignees: ''

---

<!--

Welcome!  Unfortunately not all documentation is perfect, and if you're opening
a documentation issue we are interested in fixing it.  Please fill out the
template below so that we can solve the problem more quickly; or, alternately,
open a PR with a fix, if you like.

-->

#### Problem location

<!-- Link to incorrect website or location of source file with bad
documentation. -->

#### Description of problem

<!-- Tell us what is wrong with the documentation so we can fix it. -->
---
name: Bug report
about: Use this template for reporting a bug that you have found in mlpack.
title: ''
labels: 't: bug report, s: unanswered'
assignees: ''

---

<!--

Welcome!  Please fill out the template below; that makes it easier for us to 
quickly figure out what the issue is and solve it.  Thanks!

-->

#### Issue description

<!-- Describe your issue here. -->

#### Your environment

 * version of mlpack:
 * operating system:
 * compiler:
 * version of dependencies (Boost/Armadillo):
 * any other environment information you think is relevant:

#### Steps to reproduce

<!-- Tell us how to reproduce the issue; please provide a working example if
possible! -->

#### Expected behavior

<!-- Tell us what should happen. -->

#### Actual behavior

<!-- Tell us what happened instead. -->
---
name: "🚀 Feature request"
about: Use this issue for suggesting new features
---

<!--

Welcome!  If you have a feature you would like us to implement you can suggest it with this template.
-->

### What is the desired addition or change?

<!--
  Provide a clear and concise description of what the feature is.
  For example, "I'm always frustrated when..."
-->

### What is the motivation for this feature?


### If applicable, describe how this feature would be implemented.

<!--
  You might want to link related issues or research papers here
-->

### Additional information?

<!--
  Is there anything else you can add about the proposal?
  You might want to link to related issues here, if you haven't already.
-->
---
name: Question
about: Use this template for other problems, requests, or questions.
title: ''
labels: ''
assignees: ''

---

<!--

Welcome!  If you have a question you'd like to ask, you can do it here or on the
mlpack mailing list; see also https://mlpack.org/help.html.

If you're looking for how to get involved and contribute, there's no need to
open an issue---you can see https://www.mlpack.org/community.html instead.

-->
---
title: 'mlpack 3: a fast, flexible machine learning library'
tags:
- machine learning
- deep learning
- c++
- optimization
- template metaprogramming

authors:
- name: Ryan R. Curtin
  orcid: 0000-0002-9903-8214
  affiliation: 1

- name: Marcus Edel
  orcid: 0000-0001-5445-7303
  affiliation: 2

- name: Mikhail Lozhnikov
  orcid: 0000-0002-8727-0091
  affiliation: 3

- name: Yannis Mentekidis
  orcid: 0000-0003-3860-9885
  affiliation: 5

- name: Sumedh Ghaisas
  orcid: 0000-0003-3753-9029
  affiliation: 5

- name: Shangtong Zhang
  orcid: 0000-0003-4255-1364
  affiliation: 4

affiliations:
- name: Center for Advanced Machine Learning, Symantec Corporation
  index: 1
- name: Institute of Computer Science, Free University of Berlin
  index: 2
- name: Moscow State University, Faculty of Mechanics and Mathematics
  index: 3
- name: University of Alberta
  index: 4
- name: None
  index: 5

date: 5 April 2018
bibliography: paper.bib
---

# Summary

In the past several years, the field of machine learning has seen an explosion
of interest and excitement, with hundreds or thousands of algorithms developed
for different tasks every year.  But a primary problem faced by the field is the
ability to scale to larger and larger data---since it is known that training on
larger datasets typically produces better results [@halevy2009unreasonable].
Therefore, the development of new algorithms for the continued growth of the
field depends largely on the existence of good tooling and libraries that enable
researchers and practitioners to quickly prototype and develop solutions
[@sonnenburg2007need].  Simultaneously, useful libraries must also be efficient
and well-implemented.  This has motivated our development of mlpack.

mlpack is a flexible and fast machine learning library written in C++ that has
bindings that allow use from the command-line and from Python, with support for
other languages in active development.  mlpack has been developed actively for
over 10 years [@mlpack2011, @mlpack2013], with over 100 contributors from
around the world, and is a frequent mentoring organization in the Google Summer
of Code program (\url{https://summerofcode.withgoogle.com}).  If used in C++,
the library allows flexibility with no speed penalty through policy-based design
and template metaprogramming [@alexandrescu2001modern]; but bindings are
available to other languages, which allow easy use of the fast mlpack codebase.

For fast linear algebra, mlpack is built on the Armadillo C++ matrix library
[@sanderson2016armadillo], which in turn can use an optimized BLAS
implementation such as OpenBLAS [@xianyi2018openblas] or even NVBLAS
[@nvblas] which would allow mlpack algorithms to be run on the GPU.  In
order to provide fast code, template metaprogramming is used throughout the
library to reduce runtime overhead by performing any possible computations and
optimizations at compile time.  An automatic benchmarking system is developed
and used to test the efficiency of mlpack's algorithms [@edel2014automatic].

mlpack contains a number of standard machine learning algorithms, such as
logistic regression, random forests, and k-means clustering, and also contains
cutting-edge techniques such as a compile-time optimized deep learning and
reinforcement learning framework, dual-tree algorithms for nearest neighbor
search and other tasks [@curtin2013tree], a generic optimization framework with
numerous optimizers [@curtin2017generic], a generic hyper-parameter tuner, and
other recently published machine learning algorithms.

For a more comprehensive introduction to mlpack, see the website at
\url{http://www.mlpack.org/} or a recent paper detailing the design and
structure of mlpack [@curtin2017designing].

# References

## Tutorials

Tutorials for mlpack can be found [here : mlpack tutorials](https://www.mlpack.org/doc/mlpack-git/doxygen/tutorials.html).


### General mlpack tutorials

These tutorials introduce the basic concepts of working with mlpack, aimed at developers who want to use and contribute to mlpack but are not sure where to start.

* [Building mlpack from source](https://www.mlpack.org/doc/mlpack-git/doxygen/build.html)
* [File Formats in mlpack](https://www.mlpack.org/doc/mlpack-git/doxygen/formatdoc.html)
* [Matrices in mlpack](https://www.mlpack.org/doc/mlpack-git/doxygen/matrices.html)
* [mlpack input and output](https://www.mlpack.org/doc/mlpack-git/doxygen/iodoc.html)
* [mlpack timers](https://www.mlpack.org/doc/mlpack-git/doxygen/timer.html)
* [Simple sample mlpack programs](https://www.mlpack.org/doc/mlpack-git/doxygen/sample.html)


### Method-specific tutorials

These tutorials introduce the various methods mlpack offers, aimed at users who want to get started quickly. These tutorials start with simple examples and progress to complex, extensible uses.

* [NeighborSearch tutorial (mlpack_knn / mlpack_kfn)](https://www.mlpack.org/doc/mlpack-git/doxygen/nstutorial.html)
* [LinearRegression tutorial (mlpack_linear_regression)](https://www.mlpack.org/doc/mlpack-git/doxygen/lrtutorial.html)
* [RangeSearch tutorial (mlpack_range_search)](https://www.mlpack.org/doc/mlpack-git/doxygen/rstutorial.html)
* [Density Estimation Trees tutorial (mlpack_det)](https://www.mlpack.org/doc/mlpack-git/doxygen/dettutorial.html)
* [K-Means tutorial (mlpack_kmeans)](https://www.mlpack.org/doc/mlpack-git/doxygen/kmtutorial.html)
* [FastMKS tutorial (mlpack_fastmks)](https://www.mlpack.org/doc/mlpack-git/doxygen/fmkstutorial.html)
* [Euclidean Minimum Spanning Trees tutorial (mlpack_emst)](https://www.mlpack.org/doc/mlpack-git/doxygen/emst_tutorial.html)
* [Alternating Matrix Factorization Tutorial](https://www.mlpack.org/doc/mlpack-git/doxygen/amftutorial.html)
* [Collaborative Filtering Tutorial](https://www.mlpack.org/doc/mlpack-git/doxygen/cftutorial.html)


### Policy Class Documentation

mlpack uses templates to achieve its genericity and flexibility. Some of the template types used by mlpack are common across multiple machine learning algorithms. The links below provide documentation for some of these common types.

* [The MetricType policy in mlpack](https://www.mlpack.org/doc/mlpack-git/doxygen/metrics.html)
* [The KernelType policy in mlpack](https://www.mlpack.org/doc/mlpack-git/doxygen/kernels.html)
* [The TreeType policy in mlpack](https://www.mlpack.org/doc/mlpack-git/doxygen/trees.html)
