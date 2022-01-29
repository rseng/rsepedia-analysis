# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2022-01-27

### Added

### Changed

### Deprecated

### Removed

### Fixed

### CI/CD

### Other


## [4.2.1] - 2022-01-27

### Fixed
* Implemented the optional method, minCost, of the OptimizationProblem and
  IntegerCostOptimizationProblem classes in all of the various TSP implementations
  to return a simple, extremely loose, lower bound of 0, to enable using the
  InverseCostFitnessFunction class that require a finite lower bound. The default
  implementation of that method otherwise returns negative infinity.


## [4.2.0] - 2022-01-24

### Added
* Class RandomTSPMatrix that enables generating random instances of the Traveling Salesperson
  Problem (TSP) as well as the Asymmetric TSP (ATSP) by generating a random distance matrix. This
  new class includes the option to generate a random distance matrix that enforces the triangle
  inequality. The library had previously introduced a TSP class in a prior release that bases the
  instances on points in 2D space with a user defined metric. This new RandomTSPMatrix class
  compliments that class by additionally providing the ability to specify arbitrary distance matrices
  for a TSP instance.

### Changed
* Refactored existing Traveling Salesperson Problem classes (all non-breaking changes)
  to ease the addition of other variations.

### CI/CD
* Integrated Snyk code scanning on pushes/PRs.


## [4.1.0] - 2022-01-13

### Added
* Transformation between Permutation problem and BitVector problem, enabling using operators
  for BitVectors when solving Permutation problems. 

### CI/CD
* Revised CI/CD workflows to deploy API website changes automatically upon each new release.


## [4.0.0] - 2022-01-05

### CONTAINS BREAKING CHANGES
This release contains breaking changes. See the Removed list below for details. The most significant
breaking changes relate to functionality initially introduced in v3.1.0, released on 12/21/2021, and 
thus likely affects very few users. Other breaking changes relate to removal of previously deprecated
methods.

### Added
* Added GenerationalMutationOnlyEvolutionaryAlgorithm, moving the mutation-only EA
  functionality into this new class from the existing GenerationalEvolutionaryAlgorithm
  class.
* Added GenerationalNANDOperatorsEvolutionaryAlgorithm, moving the EA
  functionality related to mutually-exclusive genetic operators (i.e., generation variation where
  each child is result of either mutation, or crossover, or identical copy of a parent, but never
  result of both mutation and crossover) into this new class from the existing 
  GenerationalEvolutionaryAlgorithm class.
* Enhanced all of the following Evolutionary Algorithm and Genetic Algorithm classes with the option to use elitism:
  * GenerationalEvolutionaryAlgorithm
  * GenerationalMutationOnlyEvolutionaryAlgorithm
  * GenerationalNANDOperatorsEvolutionaryAlgorithm
  * GeneticAlgorithm
  * MutationOnlyGeneticAlgorithm
* Crossover operators for IntegerVector and BoundedIntegerVector classes, including:
  * Single-point crossover
  * Two-point crossover
  * K-point crossover
  * Uniform crossover
* Crossover operators for RealVector and BoundedRealVector classes, including:
  * Single-point crossover
  * Two-point crossover
  * K-point crossover
  * Uniform crossover
* Enhancements to IntegerVector and BoundedIntegerVector, including:
  * IntegerVector.exchange method which exchanges a subsequence between two IntegerVectors.
  * BoundedIntegerVector.sameBounds method which checks if two bounded integer vectors are
    subject to the same min and max values.
* Enhancements to RealVector and BoundedRealVector, including:
  * RealVector.exchange method which exchanges a subsequence between two RealVectors.
  * BoundedRealVector.sameBounds method which checks if two bounded real vectors are
    subject to the same min and max values.
* Added constructors to TSP (traveling salesperson problem) classes to enable specifying
  problem instance from arrays of x and y coordinates.

### Changed
* Refactored evolutionary algorithm classes to improve maintainability, and ease
  integration of planned future functionality.
* Refactored ProgressTracker.update(SolutionCostPair) to simplify logic.

### Removed
* Removed all mutation-only EA functionality from the GenerationalEvolutionaryAlgorithm,
  moving that functionality into the new GenerationalMutationOnlyEvolutionaryAlgorithm. This
  was done for maintainability. It is a breaking-change, although it should affect minimal
  users as the GenerationalEvolutionaryAlgorithm was just introduced in v3.1.0. For those using it
  simply use the new GenerationalMutationOnlyEvolutionaryAlgorithm class where you will
  find all the same constructors and methods necessary for mutation-only EAs.
* Removed all mutually-exclusive genetic operator functionality from the GenerationalEvolutionaryAlgorithm,
  moving that functionality into the new GenerationalNANDOperatorsEvolutionaryAlgorithm. This
  was done for maintainability. It is a breaking-change, although it should affect minimal
  users as the GenerationalEvolutionaryAlgorithm was just introduced in v3.1.0. For those using it
  simply use the new GenerationalNANDOperatorsEvolutionaryAlgorithm class where you will
  find that functionality.
* All methods/constructors that were previously deprecated in earlier releases have 
  now been removed, including:
  * The following deprecated methods of the ProgressTracker class have been removed:
    * `update(int, T)` in favor of using `update(int, T, boolean)`.
    * `update(double, T)` in favor of using `update(double, T, boolean)`.
    * `setFoundBest()` in favor of using either `update(int, T, boolean)` 
      or `update(double, T, boolean)`.
  * The following deprecated constructors of the SolutionCostPair class have been removed:
    * `SolutionCostPair(T, int)` in favor of using `SolutionCostPair(T, int, boolean)`
    * `SolutionCostPair(T, double)` in favor of using `SolutionCostPair(T, double, boolean)`


## [3.1.0] - 2021-12-21

### Added
* Enhancements to BitVector class, including:
  * A new method to exchange a sequence of bits between two BitVectors.
  * A new method to exchange a selection of bits, specified with a bit mask, between two BitVectors.
  * A new constructor, BitVector(length, p), that generates random bit masks given probability p of a 1-bit.
* Generational evolutionary algorithms, including the following features and functionality:
  * A generic base class, GenerationalEvolutionaryAlgorithms, with type parameter to specify type of structure we
    are evolving, and which supports the following EA variations:
    * Standard generational model, parents replaced by children, crossover and mutation.
    * Generational with mutation-only (no crossover).
    * Generational, parents replaced by children, with mutually exclusive crossover and mutation (i.e.,
      no child is result of both crossover and mutation).
  * Subclasses of GenerationalEvolutionaryAlgorithm for Genetic Algorithm specific cases for convenience
    when optimizing BitVectors, including:
    * GeneticAlgorithm: a standard generational GA with bit flip mutation, and choice of crossover operator
      and selection operator.
    * MutationOnlyGeneticAlgorithm: a generational GA with only bit flip mutation, and choice of selection operator.
    * SimpleGeneticAlgorithm: a generational GA with bit flip mutation, single-point crossover, and fitness proportional
      selection.
  * Crossover features (mirrors the existing features of mutation operators):
    * CrossoverOperator interface for defining custom crossover operators.
    * Support for hybrid crossover operators (e.g., picking randomly from set of crossover operators).
  * Crossover operators for BitVectors, including:
    * Single-point crossover
    * Two-point crossover
    * K-point crossover
    * Uniform crossover
  * The following selection operators:
    * FitnessProportionalSelection
    * StochasticUniversalSampling
    * TournamentSelection
    * TruncationSelection
    * RandomSelection
    * BiasedFitnessProportionalSelection, which is fitness proportional selection but which
      enables transforming the fitness values by a bias function.
    * BiasedStochasticUniversalSampling, which is stochastic universal sampling but which
      enables transforming the fitness values by a bias function.
    * ShiftedFitnessProportionalSelection, which is a variation of FitnessProportionalSelection
      that uses transformed fitness values that are shifted so that the minimum is equal to 1,
      enabling safe use with negative fitness values.
    * ShiftedStochasticUniversalSampling, which is a variation of StochasticUniversalSampling
      that uses transformed fitness values that are shifted so that the minimum is equal to 1,
      enabling safe use with negative fitness values.
    * BiasedShiftedFitnessProportionalSelection, which is a variation of BiasedFitnessProportionalSelection
      that uses transformed fitness values that are shifted so that the minimum is equal to 1,
      enabling safe use with negative fitness values.
    * BiasedShiftedStochasticUniversalSampling, which is a variation of BiasedStochasticUniversalSampling
      that uses transformed fitness values that are shifted so that the minimum is equal to 1,
      enabling safe use with negative fitness values.
  * Classes to transform cost values to fitness values:
    * InverseCostFitnessFunction: Transforms optimization cost functions to fitnesses with a
      division.
    * NegativeCostFitnessFunction: Transforms cost to fitness with fitness(s)=-cost(s).
    * NegativeIntegerCostFitnessFunction: Transforms cost to fitness with fitness(s)=-cost(s),
      same as above but for integer cost problems.
* Constructive heuristics for the TSP for use by stochastic sampling algorithms, including:
  * Nearest city heuristic
  * Nearest city pair heuristic, which prefers the first city of the nearest 2-city combination
* New methods added to ProgressTracker class:
  * `update(int, T, boolean)`
  * `update(double, T, boolean)`
  * `update(SolutionCostPair<T>)`
* New constructors in SolutionCostPair class:
  * `SolutionCostPair(T, int, boolean)`
  * `SolutionCostPair(T, double, boolean)`

### Deprecated
* The following methods of the ProgressTracker class have been deprecated:
  * `update(int, T)` in favor of using `update(int, T, boolean)`.
  * `update(double, T)` in favor of using `update(double, T, boolean)`.
  * `setFoundBest()` in favor of using either `update(int, T, boolean)` 
    or `update(double, T, boolean)`.
* The following constructors of the SolutionCostPair class have been deprecated:
  * `SolutionCostPair(T, int)` in favor of using `SolutionCostPair(T, int, boolean)`
  * `SolutionCostPair(T, double)` in favor of using `SolutionCostPair(T, double, boolean)`

### CI/CD
* Updated CI/CD workflow to comment on PRs with the coverage and branches coverage percentages.


## [3.0.0] - 2021-10-25

### Added
* Implementation of the Traveling Salesperson Problem (class TSP), with the following features:
  * Generates random instances with cities distributed uniformly at random within a square.
  * Defaults to Euclidean distance, but also supports specifying a function for edge distance.
  * Two variations with a precomputed matrix of edge costs: floating-point costs, and integer costs, 
    where the integer cost version by default rounds each edge cost to nearest int, but which can be 
    customized. These use quadratic memory.
  * Two variations where edge costs are computed as needed: floating-point costs, and integer costs, 
    where the integer cost version by default rounds each edge cost to nearest int, but which can be 
    customized. These use linear memory.

### Changed
* Beginning with release 3.0.0, the minimum supported Java version is now Java 11+.
* The randomization utilities for generating Gaussian-distributed random numbers, previously 
  contained in the package org.cicirello.math.rand has been moved to a new 
  library [&rho;&mu;](https://rho-mu.cicirello.org/), which is a transitive dependency (via 
  our dependency on [JPT](https://github.com/cicirello/JavaPermutationTools)).
* Refactored GaussianMutation to improve maintainability by eliminating dependence upon a 
  specific algorithm for generating Gaussian distributed random numbers.
* The library now uses Java modules, providing the module `org.cicirello.chips_n_salsa`, which
  includes the package `org.cicirello.search` and all of its very many subpackages.
  * The required dependent module `org.cicirello.jpt` is declared with `requires transitive`
    because the `Permutation` class from one of the packages of that module is a parameter
	and/or returned by a variety of methods such as operators that manipulate permutations.
	As a consequence, projects that include Chips-n-Salsa as a dependency will also include JPT
	and its dependencies. User's dependency manager should handle this.
* Changed the default annealing schedule in the SimulatedAnnealing class to the Self-Tuning Lam
  adaptive schedule of Cicirello (2021), which is implemented in class SelfTuningLam.

### Other
* Reorganized directory hierarchy to the Maven default (we had been using a custom directory
  hierarchy for quite some time from before we switched to Maven).
  

## [2.13.0] - 2021-09-16

### Added
* The Self-Tuning Lam annealing schedule (class SelfTuningLam), as described by the
  paper: Vincent A. Cicirello. 2021. Self-Tuning Lam Annealing: Learning Hyperparameters 
  While Problem Solving. Applied Sciences, 11(21), Article 9828 (November 2021). https://doi.org/10.3390/app11219828.
* Implementation of Forrester et al (2008) continuous function optimization problem.
* Implementation of Gramacy and Lee (2012) continuous function optimization problem.

### Other
* Added CITATION.cff file to configure GitHub's citation link.
* Edited Zenodo metadata.


## [2.12.1] - 2021-7-30

### Fixed
* Corrected acceptance rate calculations in AcceptanceTracker class for
  cases where some runs terminated early upon finding solution with cost
  equal to theoretical minimum cost for problem.


## [2.12.0] - 2021-7-28

### Added
* Implementations of various Royal Road problems (optimization problems on bit vectors).
  * RoyalRoad class which implements Mitchell et al's original two variations (with
    and without stepping stones) of the Royal Road problem, but generalized to BitVectors
    of any length as well as other block sizes (Mitchell et al's original problem was on
    bit vectors of length 64 with 8-bit blocks).
  * HollandRoyalRoad class implements Holland's Royal Road function.

### Changed
* Completely redesigned the website for the library, including
  redesigning with AMP to accelerate mobile browsing.
* Upgraded JPT dependency to JPT 2.6.1.


## [2.11.1] - 2021-5-13

### Changed
* Code improvements based on report from initial run of MuseDev's 
  static code analysis tool (non-functional changes to improve code).

### Fixed
* Fixed null pointer bug in stochastic sampling implementations for the case when
  the given heuristic doesn't use an incremental evaluation.


## [2.11.0] - 2021-5-11

### Added
* New mutation operators for permutations:
    * UniformScrambleMutation: Like ScrambleMutation, but scrambled elements
      are not necessarily contiguous.
    * UndoableUniformScrambleMutation: Exactly like UniformScrambleMutation (see above),
      but with support for the undo method of the UndoableMutationOperator interface.
    * TwoChangeMutation: Mutation operator on permutations that is the equivalent of
      the classic two change operator, assuming that the permutation represents a cyclic
      sequence of edges.
    * ThreeOptMutation: Mutation operator on permutations that is the equivalent of
      the classic 3-Opt operator, which is a combination of two-changes and three-changes,
      assuming that the permutation represents a cyclic sequence of edges.
    * CycleMutation: Generates a random permutation cycle.

### Changed
* Updated dependency to JPT, v2.6.0.
* Minor optimization to UndoableScrambleMutation, utilizing the JPT updates.

### CI/CD
* Started using CodeQL code scanning on all push/pull-request events.
* Started using Muse.dev for static analysis.
* Upgraded coverage reporting to JaCoCo 0.8.7.


## [2.10.0] - 2021-3-27

### Added
* Enhanced support for parallel metaheuristics, including:
    * ParallelMetaheuristic class, which enables running multiple copies of a
      metaheuristic, or multiple metaheuristics, in parallel using multiple threads.
    * ParallelReoptimizableMetaheuristic, which enables running multiple copies of a
      metaheuristic that supports the reoptimize method, or multiple such metaheuristics, 
      in parallel using multiple threads.
* Additional artificial search landscapes, including:
    * OneMaxAckley class implements the original version of the One Max
      problem as described by Ackley (1985), whereas the existing OneMax 
      class in the library implements a variation. Ackley defined the problem
      as maximize 10 * number of one bits.
    * Plateaus class implements the Plateaus problem, originally described by Ackley (1987),
      and is an optimization problem over the space of bit strings, characterized by large
      flat areas.
    * Mix class implements the Mix problem, originally described by Ackley (1987),
      and is an optimization problem over the space of bit strings, that mixes the
      characteristics of OneMax, TwoMax, Trap, Porcupine, and Plateau problems.
* Added the following predicate methods to BitVector: allOnes, allZeros, anyOnes, anyZeros.
* Enhancements to BitVector.BitIterator, including:
    * BitVector.BitIterator.skip() method enabling skipping over bits.
    * BitVector.BitIterator.numRemaining() method to get number of remaining bits.
    * BitVector.BitIterator.nextLargeBitBlock() for getting more than 32 bits at a time.
* RotationMutation, a mutation operator for Permutations.

### Changed
* Enhanced support for parallel metaheuristics, including:
    * Refactored ParallelMultistarter class to reduce redundancy among its constructors.
    * Refactored ParallelMultistarter class to utilize the new ParallelMetaheuristic as
      a base class.
    * Refactored ParallelReoptimizableMultistarter class to utilize the 
      new ParallelReoptimizableMetaheuristic as its super class.

### Fixed
* Bug in definition of Ackley's Trap function (had been missing a 
  floor due to an error in a secondary source of description of problem).
  

## [2.9.0] - 2021-3-20

### Added
* Added an implementation of the TwoMax problem.
* Added an implementation of a variation of the TwoMax problem, but
  with two global optima, unlike the original version which has one 
  global and one sub-optimal local optima.
* Added an implementation of Ackley's Trap function, an artificial
  search landscape with one global optima, and one sub-optimal local
  optima, where most of the search space is within the attraction basin 
  of the local optima.
* Added an implementation of Ackley's Porcupine function, an artificial
  search landscape with one global optima, and an exponential number of
  local optima, a very rugged landscape.


## [2.8.0] - 2021-3-5

### Added
* Added AcceptanceTracker, which provides a mechanism for extracting 
  information about the behavior of an annealing schedule across a
  set of runs of simulated annealing.
* Added CostFunctionScaler, which enables scaling the cost function of
  an optimization problem by a positive constant, such as if you want to
  explore the effects of the range of a cost function on a stochastic local
  search algorithm by scaling a cost function with a known range.
* Added IntegerCostFunctionScaler, which enables scaling the cost function of
  an optimization problem by a positive constant, such as if you want to
  explore the effects of the range of a cost function on a stochastic local
  search algorithm by scaling a cost function with a known range. This class
  is just like CostFunctionScaler but for problems with integer costs.


## [2.7.0] - 2021-2-25

### Added
* Added HybridConstructiveHeuristic, which provides the ability
  to use multiple heuristics with a stochastic sampling search,
  where a heuristic is chosen from a set of constructive heuristics
  at the start of each iteration of the stochastic sampler and used
  for all decisions made during that iteration. The class supports the
  following strategies for selecting the next heuristic:
    * Heuristic uniformly at random from among the available heuristics.
    * Heuristic chosen using a round robin strategy that systematically cycles
      over the heuristics.
    * Heuristic chosen via a weighted random selection process.
* Added versions of a few constructive heuristics for 
  scheduling problems that precompute heuristic values 
  upon construction of the heuristic. Many of the scheduling 
  heuristic implementations already do this, where it is 
  possible to do so, provided the memory and time requirements 
  for precomputing heuristic values is linear in the number 
  of jobs. Some heuristics that consider setup times, however, 
  need quadratic memory if the heuristic is precomputed prior 
  to running the search. The existing implementations of these 
  compute the heuristic as it is needed during the search, 
  which means for a stochastic sampler with many iterations 
  that the same heuristic values will be computed repeatedly. 
  We've added versions of these that precompute the heuristic 
  values so that if your problem is small enough to afford the 
  quadratic memory, and if you will be running enough iterations 
  of the stochastic sampler to afford the quadratic time of 
  this precomputation step, then you can choose to use a 
  precomputed version.
    * A scheduling heuristic, SmallestSetupPrecompute, which is 
      a version of SmallestSetup, but which precomputes a table 
      of heuristic values to avoid recomputing the same heuristic 
      values repeatedly across multiple iterations of a stochastic 
      sampling search.
    * A scheduling heuristic, ShortestProcessingPlusSetupTimePrecompute, which is 
      a version of ShortestProcessingPlusSetupTime, but which precomputes a table 
      of heuristic values to avoid recomputing the same heuristic 
      values repeatedly across multiple iterations of a stochastic 
      sampling search.
    * A scheduling heuristic, WeightedShortestProcessingPlusSetupTimePrecompute, which is 
      a version of WeightedShortestProcessingPlusSetupTime, but which precomputes a table 
      of heuristic values to avoid recomputing the same heuristic 
      values repeatedly across multiple iterations of a stochastic 
      sampling search.

### Changed
* Dependency updated to JPT v2.4.0.


## [2.6.0] - 2021-1-25

### Added
* A factory method, ConstantRestartSchedule.createRestartSchedules, for creating multiple restart schedules.
* Factory methods, VariableAnnealingLength.createRestartSchedules, for creating multiple restart schedules.
* Factory methods, LubyRestarts.createRestartSchedules, for creating multiple restart schedules.

### Changed
* Refactored TimedParallelMultistarter and TimedParallelReoptimizableMultistarter to eliminate redundancy in common with these two classes.
* Refactored ParallelMultistarter and ParallelReoptimizableMultistarter to eliminate redundancy in common with these two classes.
* Eliminated code redundancy between ParallelReoptimizableMultistarter and TimedParallelReoptimizableMultistarter.

### Removed
* Javadoc's zipped versions of the api doc search index (and added to gitignore). The search feature on the API website will still be functional as it uses the unzipped versions as a fallback. The web server should compress these anyway. From Java 15 onward, they eliminated the generation of the zipped versions all together. We're removing them here since the api website is hosted from GitHub Pages, and git detects these as changed every time javadoc runs even if the index has not actually changed, and binary files are not efficiently stored in git.

### Fixed
* Minor documentation edits to clarify details of annealing schedules (ModifiedLam, OriginalModifiedLam, and SimulatedAnnealing classes).


## [2.5.2] - 2020-11-11

### Fixed
* Various minor bug fixes:
    * Fixed bugs in factory methods of the uniform mutation operator, the Cauchy mutation operator, and the Gaussian mutation operator for real-valued representations related to checking for invalid parameters.
    * Fixed bug in toArray method of uniform mutation operator, the Cauchy mutation operator, and the Gaussian mutation operator in case when null passed to method.
    * Fixed bug in IntegerVectorInitializer.split().


## [2.5.1] - 2020-11-5

### Changed
* Refactored ParallelVariableAnnealingLength to eliminate redundant code.
* Refactored constructors of window limited permutation mutation operators to eliminate redundant code.
* Refactored BlockMoveIterator, WindowLimitedBlockMoveIterator, and BlockInterchangeIterator to remove unnecessary conditions.
* Refactored internal methods of IntegerVectorInitializer to remove unnecessary checks.

### Fixed
* Various minor bug fixes:
    * Fixed bug in WeightedHybridMutation's and WeightedHybridUndoableMutation's constructors in the check for non-positive weights.
    * Fixed bug in window limited permutation mutation operators in handling of limit greater than or equal to permutation length.
    * Fixed missing parameter bounds check in LubyRestarts constructor.
    * Fixed bugs in factory methods of the uniform mutation operator for integer-valued representations related to checking for invalid parameters.
    * Fixed bugs in equals methods for bounded integer-valued representations in case when vectors contain same elements, but constrained by different bounds.
    * Fixed bug in IntegerVectorInitializer.split().


## [2.5.0] - 2020-11-2

### Added
* EXPET heuristic, a constructive heuristic for common duedate scheduling problems.

### Changed
* Minor optimizations in scheduling problem classes (heuristics and cost functions).

### Fixed
* Fixed bug in parallel searches involving handling exceptions thrown by component searches.
* Added missing parameter checks in constructors of parallel metaheuristic implementations.
* Added missing parameter bounds check in LINET scheduling heuristic implementation.
* Added test cases to improve coverage.

### CI/CD
* Code coverage reporting, via JaCoCo, added to build process.


## [2.4.1] - 2020-10-15

### Fixed
* Bug in release workflow: attaching jars to GitHub Release.


## [2.4.0] - 2020-10-15

### Added
* VBSS support for exponential bias functions.

### Changed
* Build process migrated from Ant to Maven.
* Improved formatting of description field in pom.xml.

### Removed
* The example programs have been removed from the Chips-n-Salsa repository, and now reside in their own dedicated repository: https://github.com/cicirello/chips-n-salsa-examples
* Distribution jar files removed from repository.  The distribution jars are available via GitHub Releases, GitHub Packages, and Maven Central.


## [2.3.0] - 2020-10-13

### Added
* Problem instance generator for Common Duedate Scheduling, based on the method used to generate the instances from the OR-Lib.
* Parser for the file format specified in OR-Lib for the Common Duedate Scheduling benchmarks.
* WLPT heuristic, a constructive heuristic for common duedate scheduling problems.
* LINET heuristic, a constructive heuristic for common duedate scheduling problems.



## [2.2.0] - 2020-09-24

### Added
* ModifiedLamOriginal: This class is a direct implementation of the Modified Lam annealing schedule as originally described by Swartz and Boyan. Introduced to the library to enable comparing with the optimized version already included in the library in the class ModifiedLam.
* Instance generator for static scheduling with weights and duedates.

### Changed
* Upgraded dependency to JPT 2.2.0.
* Modifications to API documentation website to improve browsing on mobile web browsers.


## [2.1.0] - 2020-09-18

### Added
* Parser for benchmark scheduling instance data files for single machine scheduling problems with weights, duedates, and sequence-dependent setup times.

### Changed
* Revised README to include instructions for importing from Maven Central
* API documentation website (https://chips-n-salsa.cicirello.org) updated to html5


## [2.0.0] - 2020-09-16

First release available through Maven Central. New major release, 2.0.0, due 
to incompatible interface changes within the org.cicirello.search.ss package. Code 
that doesn't depend upon any of the stochastic sampling search algorithm 
implementations should be unaffected by these changes, and in those cases should 
be safe to simply upgrade to this new version without need to change anything in 
dependent code.

### Added
* An interface Partial to enable generalizing the various stochastic sampling search algorithms from optimizing the space of permutations to more general types.
* A method was added to ConstructiveHeuristic for creating empty Partials of the appropriate length needed by the heuristic.
* Class HeuristicSolutionGenerator, which is a generalization of existing class HeuristicPermutationGenerator, for constructing solutions to optimization problems via constructive heuristics (HeuristicPermutationGenerator now extends this new class).
* PartialIntegerVector class to support using stochastic sampling search algorithms to generate vectors of integers.
* BoundedIntegerVector and BoundedRealVector classes added to the package org.cicirello.search.representations

### Changed
* ConstructiveHeuristic and IncrementalEvaluation now have a type parameter.
* ValueBiasedStochasticSampling, HeuristicBiasedStochasticSampling, and AcceptanceBandSampling also now have a type parameter.
* All of the scheduling heuristics were modified based on interface changes.
* Renamed ConstructiveHeuristic.completePermutationLength to ConstructiveHeuristic.completeLength to be more general.
* Minor edits to the example program SchedulingWithVBSS related to other changes (added missing type parameter to construction of the ValueBiasedStochasticSampling instance to avoid compiler warning).
* Workflows and pom.xml to enable publishing to both Maven Central as well as Github Packages


## [1.4.0] - 2020-08-18
### Added
* .zenodo.json to provide project metadata to Zenodo during releases
* Enabled DependaBot to keep dependencies up to date.

### Changed
* Refactored SimulatedAnnealing to simplify implementation.
* Minor modifications to example programs that involve the SimulatedAnnealing class to use simplified constructors.
* Refactored HBSS, VBSS, AcceptanceBandSampling, and IterativeSampling to use a common abstract base class to remove code redundancy.
* Refactored HeuristicPermutationGenerator to simplify implementation.
* Updated dependency versions.

## [1.3.0] - 2020-08-06
### Added
* ATCS scheduling heuristic (as described by the heuristic's authors).
* DynamicATCS scheduling heuristic (dynamically updates heuristic's parameters during schedule construction).
* An example program for an industrial scheduling problem, demonstrating using a stochastic sampling search algorithm, VBSS, as well as a hybrid of VBSS and hill climbing.
* CHANGELOG.md
* Automated sitemap.xml generation for the API website.

### Changed
* Improved representation of setup times for scheduling problems with sequence-dependent setups.
* Refactored ATC and Montagne scheduling heuristics, and scheduling heuristic abstract base class for performance and style.
* Refactored hill climbing classes, introducing an abstract base class, to eliminate redundant code.

### Removed
* Unnecessary maven settings.xml.

### Fixed
* Fixed classpath bug in Makefile (Windows use of semicolon vs Linux, etc use of colon as delimiter). Now adjusts based on OS.

## [1.2.0]
### Added
* Several constructive heuristics for scheduling problems for use with stochastic sampling search algorithms.
* Heuristics include: EDD, MST, SPT, WSPT, Smallest Setup, Smallest Two-Job Setup, Montagne, Weighted COVERT, Weighted Critical Ratio, ATC.
* Heuristics also include variations of several of the above that adjust for setup times for scheduling problems with sequence-dependent setups.
* Contributing guidelines.

### Changed
* Refactored scheduling heuristic classes to improve readability.
* Some code optimizations within the scheduling heuristic classes.
* Improvements to README.

## [1.1.0]
### Added
* Constructive heuristics for single machine scheduling problems.
* Prebuilt jars now included with release.

## [1.0.2]
### Changed
* Minor non-functional changes related to code style.
* Edits to JOSS paper submission.

## [1.0.1]
### Other
* Added a workflow to the repo to automatically publish to the GitHub Package Registry on release. This release is simply to confirm that it works and to publish the package to the GPR.

## [1.0.0]
### First public release
Chips-n-Salsa is a Java library of customizable, hybridizable, iterative, parallel, stochastic, and self-adaptive local search algorithms. The library includes implementations of several stochastic local search algorithms, including simulated annealing, hill climbers, as well as constructive search algorithms such as stochastic sampling. The library most extensively supports simulated annealing. It includes several classes for representing solutions to a variety of optimization problems. For example, the library includes a BitVector class that implements vectors of bits, as well as classes for representing solutions to problems where we are searching for an optimal vector of integers or reals. For each of the built-in representations, the library provides the most common mutation operators for generating random neighbors of candidate solutions. Additionally, the library provides extensive support for permutation optimization problems, including implementations of many different mutation operators for permutations, and utilizing the efficiently implemented Permutation class of the JavaPermutationTools (JPT) library.

Chips-n-Salsa is customizable, making extensive use of Java's generic types, enabling using the library to optimize other types of representations beyond what is provided in the library. It is hybridizable, providing support for integrating multiple forms of local search (e.g., using a hill climber on a solution generated by simulated annealing), creating hybrid mutation operators (e.g., local search using multiple mutation operators), as well as support for running more than one type of search for the same problem concurrently using multiple threads as a form of algorithm portfolio. Chips-n-Salsa is iterative, with support for multistart metaheuristics, including implementations of several restart schedules for varying the run lengths across the restarts. It also supports parallel execution of multiple instances of the same, or different, stochastic local search algorithms for an instance of a problem to accelerate the search process. The library supports self-adaptive search in a variety of ways, such as including implementations of adaptive annealing schedules for simulated annealing, such as the Modified Lam schedule, implementations of the simpler annealing schedules but which self-tune the initial temperature and other parameters, and restart schedules that adapt to run length.
# Chips-n-Salsa - A Java library of customizable, hybridizable, iterative, parallel, stochastic, and self-adaptive local search algorithms

[<img alt="Chips-n-Salsa - A Java library of customizable, hybridizable, iterative, parallel, stochastic, and self-adaptive local search algorithms" 
     src="../gh-pages/images/chips-n-salsa.png" width="640">](#chips-n-salsa---a-java-library-of-customizable-hybridizable-iterative-parallel-stochastic-and-self-adaptive-local-search-algorithms)

Copyright (C) 2002-2022 [Vincent A. Cicirello](https://www.cicirello.org/).

Website: https://chips-n-salsa.cicirello.org/

API documentation: https://chips-n-salsa.cicirello.org/api/

| __Publications About the Library__ | [![DOI](https://joss.theoj.org/papers/10.21105/joss.02448/status.svg)](https://doi.org/10.21105/joss.02448) |
| :--- | :--- |
| __Packages and Releases__ | [![Maven Central](https://img.shields.io/maven-central/v/org.cicirello/chips-n-salsa.svg?label=Maven%20Central&logo=apachemaven)](https://search.maven.org/artifact/org.cicirello/chips-n-salsa) [![GitHub release (latest by date)](https://img.shields.io/github/v/release/cicirello/Chips-n-Salsa?logo=GitHub)](https://github.com/cicirello/Chips-n-Salsa/releases) |
| __Build Status__ | [![build](https://github.com/cicirello/Chips-n-Salsa/workflows/build/badge.svg)](https://github.com/cicirello/Chips-n-Salsa/actions/workflows/build.yml) [![docs](https://github.com/cicirello/Chips-n-Salsa/workflows/docs/badge.svg)](https://chips-n-salsa.cicirello.org/api/) [![CodeQL](https://github.com/cicirello/Chips-n-Salsa/actions/workflows/codeql-analysis.yml/badge.svg)](https://github.com/cicirello/Chips-n-Salsa/actions/workflows/codeql-analysis.yml) |
| __JaCoCo Test Coverage__ | [![coverage](../badges/jacoco.svg)](https://github.com/cicirello/Chips-n-Salsa/actions/workflows/build.yml) [![branches coverage](../badges/branches.svg)](https://github.com/cicirello/Chips-n-Salsa/actions/workflows/build.yml) |
| __Security__ | [![Snyk security score](https://snyk-widget.herokuapp.com/badge/mvn/org.cicirello/chips-n-salsa/badge.svg)](https://snyk.io/vuln/maven%3Aorg.cicirello%3Achips-n-salsa) [![Snyk Known Vulnerabilities](https://snyk.io/test/github/cicirello/Chips-n-Salsa/badge.svg)](https://snyk.io/test/github/cicirello/Chips-n-Salsa) |
| __DOI__ | [![DOI](https://zenodo.org/badge/273074441.svg)](https://zenodo.org/badge/latestdoi/273074441) |
| __License__ | [![GitHub](https://img.shields.io/github/license/cicirello/Chips-n-Salsa)](https://github.com/cicirello/Chips-n-Salsa/blob/master/LICENSE) |

## How to Cite

If you use this library in your research, please cite the following paper:

> Cicirello, V. A., (2020). Chips-n-Salsa: A Java Library of Customizable, Hybridizable, Iterative, Parallel, Stochastic, and Self-Adaptive Local Search Algorithms. *Journal of Open Source Software*, 5(52), 2448, https://doi.org/10.21105/joss.02448 .

## Overview

Chips-n-Salsa is a Java library of customizable, hybridizable, iterative, parallel, stochastic, and self-adaptive local search algorithms. The library includes implementations of several stochastic local search algorithms, including simulated annealing, hill climbers, as well as constructive search algorithms such as stochastic sampling. Chips-n-Salsa now also includes genetic algorithms as well as evolutionary algorithms more generally. The library very extensively supports simulated annealing. It includes several classes for representing solutions to a variety of optimization problems. For example, the library includes a BitVector class that implements vectors of bits, as well as classes for representing solutions to problems where we are searching for an optimal vector of integers or reals. For each of the built-in representations, the library provides the most common mutation operators for generating random neighbors of candidate solutions, as well as common crossover operators for use with evolutionary algorithms. Additionally, the library provides extensive support for permutation optimization problems, including implementations of many different mutation operators for permutations, and utilizing the efficiently implemented Permutation class of the [JavaPermutationTools (JPT)](https://jpt.cicirello.org/) library.

Chips-n-Salsa is customizable, making extensive use of Java's generic types, enabling using the library to optimize other types of representations beyond what is provided in the library. It is hybridizable, providing support for integrating multiple forms of local search (e.g., using a hill climber on a solution generated by simulated annealing), creating hybrid mutation operators (e.g., local search using multiple mutation operators), as well as support for running more than one type of search for the same problem concurrently using multiple threads as a form of algorithm portfolio. Chips-n-Salsa is iterative, with support for multistart metaheuristics, including implementations of several restart schedules for varying the run lengths across the restarts. It also supports parallel execution of multiple instances of the same, or different, stochastic local search algorithms for an instance of a problem to accelerate the search process. The library supports self-adaptive search in a variety of ways, such as including implementations of adaptive annealing schedules for simulated annealing, such as the Modified Lam schedule, implementations of the simpler annealing schedules but which self-tune the initial temperature and other parameters, and restart schedules that adapt to run length.

## Java 11+

Beginning with version 3.0.0, the library supports Java 11+. Our development process utilizes OpenJDK 11, and
all jar files released to Maven Central, GitHub Packages, and GitHub Releases are built for a Java 11 target.

## Versioning Scheme

Chips-n-Salsa uses [Semantic Versioning](https://semver.org/) with 
version numbers of the form: MAJOR.MINOR.PATCH, where differences 
in MAJOR correspond to incompatible API changes, differences in MINOR 
correspond to introduction of backwards compatible new functionality, 
and PATCH corresponds to backwards compatible bug fixes.

## Building the Library (with Maven)

The Chips-n-Salsa library is built using Maven. The root of the
repository contains a Maven `pom.xml`.  To build the library, 
execute `mvn package` at the root of the repository, which
will compile all classes, run all tests, run javadoc, and generate 
jar files of the library, the sources, and the javadocs.  In addition
to the jar of the library itself, this will also generate a
jar of the library that includes all dependencies.  The file names
make this distinction explicit.  All build outputs will then
be found in the directory `target`.

To include generation of a code coverage report during the build,
execute `mvn package -Pcoverage` at the root of the repository to 
enable a Maven profile that executes JaCoCo during the test phase.

## Example Programs

There are several example programs available in a separate repository:
[cicirello/chips-n-salsa-examples](https://github.com/cicirello/chips-n-salsa-examples).
The examples repository contains example usage of several of the 
classes of the library. Each of the examples contains detailed
comments within the source code explaining the example. Running the 
examples without reading the source comments is not advised.   

## Java Modules

This library provides a Java module, `org.cicirello.chips_n_salsa`. If your Java 11+
project utilizes modules, then add the following to your `module-info.java`:

```Java
module your.module.name.here {
	requires org.cicirello.chips_n_salsa;
}
```

The `org.cicirello.chips_n_salsa` module, in turn, transitively requires
the `org.cicirello.jpt` module because the functionality for optimizing
permutational problems uses the Permutation class as both parameters and
return types of methods, so if you need the functionality of that dependency
in your application, the Chips-n-Salsa dependencies, including the module
requirements will handle this for you.

## Importing the Library from Maven Central

Add this to the dependencies section of your pom.xml, replacing 
the version number with the version that you want to use.

```XML
<dependency>
  <groupId>org.cicirello</groupId>
  <artifactId>chips-n-salsa</artifactId>
  <version>3.0.0</version>
</dependency>
```

## Importing the Library from GitHub Packages

If you'd prefer to import from GitHub Packages, rather than Maven Central,
then: (1) add the dependency as indicated in previous section above, and (2) add 
the following to the repositories section of your pom.xml:

```XML
<repository>
  <id>github</id>
  <name>GitHub cicirello Apache Maven Packages</name>
  <url>https://maven.pkg.github.com/cicirello/Chips-n-Salsa</url>
  <releases><enabled>true</enabled></releases>
  <snapshots><enabled>true</enabled></snapshots>
</repository>
```

## Downloading Jar Files

If you don't use a dependency manager that supports importing from Maven Central,
or if you simply prefer to download manually, prebuilt jars are also attached to 
each [GitHub Release](https://github.com/cicirello/Chips-n-Salsa/releases).

In addition to the regular jar of the library, we also regularly publish 
a `jar-with-dependencies`. The `jar-with-dependencies` is built for a Java 11 target,
but does not contain any module declarations (unlike the regular jar file). Therefore, the
`jar-with-dependencies` should only be used on the classpath.

## License

The Chips-n-Salsa library is licensed under the [GNU General Public License 3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Contribute

If you would like to contribute to Chips-n-Salsa in any way, such 
as reporting bugs, suggesting new functionality, or code contributions 
such as bug fixes or implementations of new functionality, then start 
by reading 
the [contribution guidelines](https://github.com/cicirello/.github/blob/main/CONTRIBUTING.md).
This project has adopted 
the [Contributor Covenant Code of Conduct](https://github.com/cicirello/.github/blob/main/CODE_OF_CONDUCT.md).
---
title: 'Chips-n-Salsa: A Java Library of Customizable, Hybridizable, Iterative, Parallel, Stochastic, and Self-Adaptive Local Search Algorithms'
tags:
  - stochastic local search
  - parallel metaheuristics
  - self-adaptive
  - simulated annealing
  - discrete optimization
  - combinatorial optimization
  - Java
authors:
  - name: Vincent A. Cicirello
    orcid: 0000-0003-1072-8559
    affiliation: 1
affiliations:
  - name: Computer Science, School of Business, Stockton University, Galloway, NJ 08205
    index: 1
date: 23 July 2020
bibliography: paper.bib
---

# Summary

Discrete optimization problems are often of practical real-world importance as well as computationally intractable. For example, the traveling salesperson, bin packing, and longest common subsequence problems are NP-Hard, as is resource constrained scheduling, and many single-machine scheduling problems [@Garey1979]. Polynomial time algorithms for such problems are unlikely to exist, and the best known algorithms that guarantee optimal solutions have a worst-case exponential runtime. It is thus common to use stochastic local search and other metaheuristics [@Gonzalez2018]. Stochastic local search algorithms begin at a random search state, and apply a sequence of neighbor transitions to nearby search states. This includes perturbative [@Hoos2018] algorithms like simulated annealing [@Delahaye2019] and hill climbers [@Hoos2018], where each search state is a complete candidate feasible solution, and a mutation operator makes a small random modification to move to another local candidate solution; and also includes constructive [@Hoos2018] algorithms like stochastic samplers [@Grasas2017; @ReyesRubiano2020; @Cicirello2005; @Bresina1996; @Langley1992], where each search state is a partial solution that is iteratively transformed into a complete solution. Stochastic local search algorithms do not guarantee optimal solutions. However, they often find near-optimal solutions in much less time than systematic search. They also offer an anytime property [@Jesus2020; @Zilberstein1996], where solution quality improves with runtime.

Chips-n-Salsa is a Java library of customizable, hybridizable, iterative, parallel, stochastic, and self-adaptive local search algorithms. Its focus is discrete optimization, but also supports continuous optimization. The library provides a variety of solution representations, including BitVector and IntegerVector classes, indexable vectors of bits and integers, respectively, along with corresponding mutation operators. The library utilizes our significant prior research on permutation optimization problems, providing an extensive set of mutation operators for permutations [@Cicirello2016; @Cicirello2013], including window-limited mutation [@Cicirello2014]. It uses the JavaPermutationTools (JPT) [@Cicirello2018b] library for efficiently representing solutions to such problems. For continuous problems, Chips-n-Salsa provides a RealVector class, Gaussian mutation [@Petrowski2017], Cauchy mutation [@Petrowski2017], and uniform mutation. The library includes optimization problems useful for benchmarking metaheuristic implementations, such as the well-known OneMax problem [@Doerr2019], BoundMax (generalization of OneMax to integers), the Permutation in a Haystack problem [@Cicirello2016], and polynomial root finding.

The repository (https://github.com/cicirello/Chips-n-Salsa) contains the library source code, and programs with examples of key functionality. API and other documentation is hosted on the web (https://chips-n-salsa.cicirello.org/). The library can be integrated into projects using maven via GitHub Packages.

# Chips-n-Salsa Features and Functionality:

* __Stochastic Local Search Algorithms__: The library supports simulated annealing [@Delahaye2019], with all of the common annealing schedules, as well as advanced annealing schedules. It includes common forms of hill climbing [@Hoos2018], and several stochastic sampling [@Grasas2017] algorithms, such as iterative sampling [@Langley1992], Heuristic Biased Stochastic Sampling (HBSS) [@Bresina1996], Value Biased Stochastic Sampling (VBSS) [@Cicirello2005; @Cicirello2003], and stochastic sampling with acceptance bands [@Oddi1997; @Gomes1998]. We optimized random number generation based on our prior research to minimize the runtime impact of one of the more costly operations of a metaheuristic [@Cicirello2018a].
* __Customizable__: Chips-n-Salsa uses generic types enabling using the library to optimize other representations beyond what is provided in the library. It also enables easily integrating new mutation operators, and customizing all stages of the local search (e.g., choice of annealing schedule for simulated annealing).
* __Hybridizable__: Chips-n-Salsa supports integrating multiple forms of local search (e.g., hybrids of hill climbing with simulated annealing), creating hybrid mutation operators (e.g., combining multiple mutation operators), and running more than one type of search for the same problem concurrently using multiple threads as a form of algorithm portfolio [@Tong2019; @Gomes2001].
* __Iterative__: Chips-n-Salsa supports multistart metaheuristics, including several restart schedules [@Luby1993; @Cicirello2017] for varying the run lengths across the restarts. 
* __Parallel__: Chips-n-Salsa enables parallel execution of multiple instances of the same, or different, stochastic local search algorithms for a problem instance to accelerate search by exploiting multicore architectures [@Cicirello2017]. 
* __Self-Adaptive__: Chips-n-Salsa includes adaptive annealing schedules [@Hubin2019; @Stefankovic2009] for simulated annealing, such as Modified Lam [@Boyan1998; @Cicirello2007], self-tuning variations of the simpler annealing schedules, and adaptive restart schedules [@Cicirello2017].

# Statement of Need

The target audience of Chips-n-Salsa includes those conducting computational research in diverse domains, wherever applications of NP-Hard optimization problems is important. Existing local search open source implementations are often intimately tied to a specific problem, and usually a specific form of local search. For example, a GitHub search finds countless solvers for problems like Boolean Satisfiability or the Traveling Salesperson. Such problem-dependent implementations are not easily adapted to new problems. There are notable exceptions. In particular, ECJ [@Scott2019] and Jenetics [@Jenetics] are both mature and well-maintained Java libraries for genetic algorithms and related forms of evolutionary computation. They are problem-independent, support multiple representations, and include multithreaded capabilities. Their focus, however, is on population-based evolutionary search; whereas Chips-n-Salsa is focused on single-solution algorithms. Another especially noteworthy project is emili [@Pagnozzi2019], which is a C++ framework supporting hybrid, single-solution, stochastic local search algorithms. Chips-n-Salsa is problem-independent, representation-independent, and supports hybrid local search, as well as parallel execution to easily exploit multicore architectures to speed up problem solving. Furthermore, the self-adaptive and self-tuning features streamline development, eliminating the need for the developer to tune control parameters.

# References
