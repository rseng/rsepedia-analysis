# mctreesearch4j

![Static](https://img.shields.io/static/v1?label=docs&message=latest&color=blue&style=flat-square)
![kotlin-version](https://kotlin-version.aws.icerock.dev/kotlin-version?group=ca.aqtech&name=mctreesearch4j)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/ca.aqtech/mctreesearch4j/badge.svg?style=plastic)](https://maven-badges.herokuapp.com/maven-central/ca.aqtech/mctreesearch4j)
[![DOI](https://zenodo.org/badge/311016146.svg)](https://zenodo.org/badge/latestdoi/311016146)

An robust Implementation of the Monte Carlo Tree Search (MCTS) algorithm for the JVM.

## Objective

Flexible implementations of Monte Carlo Tree Search (MCTS), combined with domain specific knowledge and hybridization with other search algorithms, can be a very powerful for the solution of problems in complex planning. We introduce *mctreesearch4j*, a standard MCTS implementation written as a standard JVM library following key design principles of object oriented programming. We define key class abstractions allowing the MCTS library to flexibly adapt to any well defined Markov Decision Process or turn-based adversarial game. Furthermore, our library is designed to be modular and extensible, utilizing class inheritance and generic typing to standardize custom algorithm definitions. We demonstrate that the design of the MCTS implementation provides ease of adaptation for unique heuristics and customization across varying Markov Decision Process (MDP) domains. In addition, the implementation is reasonably performant and accurate for standard MDP’s. In addition, via the implementation of *mctreesearch4j*, the nuances of different types of MCTS algorithms are discussed.

## Structure

The library contains the implementation of the MCTS algorithm for single player and adversarial games. The project is compiled using `gradle`.

```
lib
└──
    ...
    └── src.main.kotlin
        └── ca.aqtech.mctreesearch4j
            ├── MDP.kt
            ├── Node.kt
            ├── Solver.kt
            ├── StateNode.kt
            ├── StatefulSolver.kt
            ├── ActionNode.kt
            └── StatelessSolver.kt

```

The MCTS lib can be implemented by creating an interface implementing a Markov Decision Process (MDP) class, to communicate with the MCTS solver. The MCTS solver comes in two flavours `StatefulSolver.kt` and `StatelessSolver.kt`. The `StatefulSolver.kt` records both the action history and state, whereas the `StatelessSolver.kt` only records actions taken, and infers the state from actions taken.

## Defining an MDP

`MDP.kt` defines the interface of the MDP. To run the MCTS solver, you must extend this `MDP.kt` class in any JVM language, and provide your specific definitions.

```kotlin
abstract class MDP<StateType, ActionType> {
    abstract fun transition(state: StateType, action: ActionType) : StateType
    /* Class method that returns a legal state given current state (StateType) and action taken (ActionType) */

    abstract fun reward(previousState: StateType?, action: ActionType?, state: StateType) : Double
    /* Class method to return a reward (Double) given state transitions parameters */

    abstract fun initialState() : StateType
    /* Return the initial state of the MCTS (StateType) */

    abstract fun isTerminal(state: StateType) : Boolean
    /* Return boolean indicating if the state is terminal. */

    abstract fun actions(state: StateType) : Iterable<ActionType>
    /* Return an Iterable of legal actions given a current state. */
}
```

## Examples

In the `app/` folder of this repository, you may see various examples of game engines interfacing with the MDP controller.

### Running Gradle examples via Gradle

It is possible to run the examples directly via Gradle. Simply select from the list of MDP / Games to run. Here is the list of arguments `args`:

- `ReversiGame`: Initiates a game of Reversi between the user and AI, complete with user interface.
- `ReversiAdversarialSim`: Generates a simulation of competitive Reversi between 2 AI's.
- `PushYourLuckSim`: Initiates an MDP (Push Your Luck) and produces a policy for a single state.
- `GridWorld`: Initiates an MDP and produces a policy for a state in the game of GridWorld.
- `Game2048`: Initiates an MDP and produces a policy for a state in the game of the Game 2048.

To run the tasks via Gradle simply run,

```bash
gradle run --args="{arg1} {arg2}"
```

Where `arg1`, `arg2` for example are the strings representing the task to run, seperated by spaces. For example:

```bash
gradle run --args="ReversiGame GridWorld"
```

Or simply to run one example,

```bash
gradle run --args="ReversiGame"
```

Furthermore, there is an example of *mctreesearch4j* integrated with the game of [Connect-4 in an example written in Scala](https://github.com/larkz/connect4-scala).


### Compiling the Core Library

To perform a local compilation of the core library located in the `lib/` folder, the user must omit `uploadArchives` step defined in the `gradle.build` file, which signs the .jar artifcat via GPG key and uploads to Maven Central Repository. The GPG key is tied to the Maven release of this repository, and is not intended for local compilation. To compile a local .jar without signing it for release, simply run,

```bash
gradle clean build -x uploadArchives
```

This will effectively exclude the upload to Maven step which is automated in this codebase, and allow the user to compile a local .jar for testing and experimentation. For a complete release to the Maven repository as a signed artifact, a PR to the main branch must be submitted, reviewed, and approved.


## Maven Central

The JVM artifact is available directly via the [Maven Central repository - mctreesearch4j](https://search.maven.org/artifact/ca.aqtech/mctreesearch4j). In principle, it
has full compatibility with any JVM language.

### Kotlin Integration (groovy)
```groovy
dependencies {
    implementation "ca.aqtech.mctreesearch4j:0.0.4"
}
```

### Scala Integration (SBT)

```sbt
libraryDependencies ++= Seq(
  "ca.aqtech" % "mctreesearch4j" % "0.0.4"
)
```

### Cite mctreesearch4j

Currently, to cite our framework, please cite our preprint:

```
@article{mctreesearch4j,
  author = {Liu, Larkin and Luo, Jun Tao},
  title = {An Extensible and Modular Design and Implementation of Monte Carlo Tree Search for the JVM},
  journal = {arXiv},
  year = {2021},
  url = {https://arxiv.org/abs/2108.10061}
 }
```
# Contributing to mctreesearch4j
We love your input! We want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

## We Develop with GitHub
We use GitHub to host code, to track issues and feature requests, as well as accept pull requests.

## All Code Changes Happen Through Pull Requests
Pull requests are the best way to propose changes to the codebase. We actively welcome your pull requests:

1. Fork the repo and create your branch from `main` or `dev`.
2. If you've added code that should be tested, add tests.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Make sure your code lints.
6. Issue that pull request to `dev` branch.
7. Code will be put into a trial period before a Maven release is ready to be released, and merged with `main` branch.

## Report bugs using GitHub's [issues](https://github.com/aqtech-ca/mctreesearch4j/issues)
We use GitHub issues to track public bugs. Report a bug by [opening a new issue](); it's that easy!

## Write bug reports with detail, background, and sample code
[This is an example](http://stackoverflow.com/q/12488905/180626) of a good bug report. Here's [another example from Craig Hockenberry](http://www.openradar.me/11905408), an great app developer.

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give sample code if you can!
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

## License
By contributing, you agree that your contributions will be licensed under the MIT License. In short, when you submit code changes, your submissions are understood to be under the same [MIT License](http://choosealicense.com/licenses/mit/) that covers the project. Feel free to contact the maintainers if that's a concern.

## Use a Consistent Coding Style

We use standard [Kotlin coding conventions](https://kotlinlang.org/docs/coding-conventions.html).


## Continuing Support

For additional support or questions related to the software you may reach out to `larkin@aqtech.ca`.

## References
This document was adapted from the open-source contribution guidelines from [briandk](https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62).

---
title: 'mctreesearch4j: A Monte Carlo Tree Search Implementation for the JVM'
tags:
  - Java Virtual Machine
  - Monte Carlo Tree Search
  - Markov Decision Process
  - Software Design
authors:
  - name: Larkin Liu
    orcid: 0000-0002-9375-1035
    affiliation: 1
  - name: Jun Tao Luo
    orcid: 0000-0002-2681-8922
    affiliation: 2
affiliations:
 - name: University of Toronto
   index: 1
 - name: Carnegie Mellon University
   index: 2
date: 30 August 2021
bibliography: paper.bib
---

# Summary

We introduce *mctreesearch4j*, a Monte Carlo Tree Search (MCTS) implementation written as a standard JVM (Java Virtual Machine) library following key object oriented programming design principles. This modular and extensibile implementation of MCTS provides a powerful tool to enable the discovery of approximate solutions to complex planning problems via rapid experimentation. *mctreesearch4j* utilizes class inheritance and generic types to standardize custom algorithm definitions. In addition, key class abstractions are designed for the library to flexibly adapt to any well-defined Markov Decision Process (MDP) or turn-based adversarial games. Furthermore, *mctreesearch4j* is capable of customization across a variety of MDP domains, consequently enabling the adoption of MCTS heuristics and customization into the core library with ease.

# Statement of Need

Open source implementations of Monte Carlo Tree Search exist, but have not gained widespread adoption. Implementations of MCTS have been written for Java Virtual Machine (JVM), C/C++, and Python, yet many do not provide easy access to heuristics, nor do they implement extensible and modular state action space abstractions [@Cowling:2012] [@Lucas:2018]. *mctreesearch4j* aims to build on the performance advantages of a compiled language families, such as JVM, while simultaneously being extensible and modular. Therefore, many improvements and research experiments could be performed by extending and modifying the capabilities of the base software. *mctreesearch4j* is designed to enable researchers to experiment with various MCTS strategies while standardizing the functionality of the MCTS solver and ensuring reliability, where the experiments are reproducible and the solver is compatible with common JVM runtimes.

In the domain of probabilistic planning, extensible, reproducible, and modular software packages do exist in the solution space of Reinforcement Learning. We take inspiration from software libraries such as Acme [@acme:2020] and Tensorforce [@tensorforce:2018], notably from its design concept. Such libraries function as a domain independent software library to generate solutions for a variety of MDP’s. The MDP’s are defined separate from the core library, and communicate via a common interface. Such a design concept has also been incorporated to *mctreesearch4j* within the solution space of MCTS for probabilistic and deterministic planning.

# Monte Carlo Tree Search

Monte Carlo Tree Search primarily makes use of a deterministic selection of actions and stochastic simulations to estimate the reward function of well defined Markov Decision Process (MDP). MCTS is a tree search adaptation of the UCB1 Multi-Armed Bandit Strategy [@Auer:2002]. In more detail, the MAB does not change state with each action taken, and is primarily focused on optimizing exploration - exploring new actions - versus exploitation - obtaining rewards from known actions. MCTS, is the extension, where we optimize for exploration vs. exploitation for dynamic states contingent on any actions taken in the action space. We present a more detailed discussion in [@Liu:2021].

The MCTS algorithm is distinctly divided into 4-phases, *Selection*, *Expansion*, *Simulation*, and *Backpropagation*, which are illustrated in Fig. \ref{fig:mcts-diagram}. In *Selection*, a policy deterministically selects which action to play, based on previously expanded states, guided by some exploration constant, *C*. In the *Expansion* phase, states that are unexplored, represented by a leaf node, are added to the search tree. Subsequently, in the *Simulation* phase, a simulation is stochastically played out. Finally, *Backpropagation* propagates the final reward of either a terminal state, or a node at a specified depth limit, back to the root node. This 4-phase process is repeated until a maximum number of iterations or a convergence criteria is established.

![Monte Carlo Tree Search Key Mechanisms. [@Browne:2012] \label{fig:mcts-diagram}](mcts-diagram-v2.png?raw=true "Title")

# Design Principles

*mctreesearch4j* is designed follow three key design principles.

- **Adaptability**: Adaptability is defined as the ability for MDP domain to be easily integrated into the *mctreesearch4j* framework using provided class abstractions. Our implementation seeks to simplify the adoption of MCTS solutions for a variety of domains.
- **JVM Compatibility**: We design a library that is fully compatible with the Java Virtual Machine (JVM), and consequently functional with any JVM languages, ie. Java, Scala, Kotlin etc.
- **Extensibility** We design to achieve a high degree of extensibility and modularity within the framework. Extensibility is defined as the ability for key mechanisms to be reused, redefined, and enhanced, without sacrificing interoperability.

*mctreesearch4j* serves as a framework to standardize the development of Monte Carlo Tree Search development in a research setting, while simultaneously serving as a production-ready software for deployment of MCTS for the JVM ecosystem. *mctreesearch4j* is specifically designed to be extensible, standardized, and modular. This standardization is key to enabling scientists and engineers to better reproduce results of AI research, and furthermore reduce the iteration speed of research into novel algorithms in AI.

# Domain Abstraction

The main abstraction that is used to define an MDP problem is the abstract class defined in *MDP Abstraction*, using generic type variables. Each of the methods correspond to specific behaviour of a discrete MDP. In *mctreesearch4j* we use generic types ``StateType`` and ``ActionType`` to represent the MDP states and actions respectively. This abstract class has five members that must be implemented. These abstract class methods define the functionality of an MDP. The MDP abstraction will be used by MCTS solvers to compute the optimal policy. The MDP interface can be written in any JVM language, we use Kotlin and Scala for this paper, with the Scala implementation from [@Liu:2021-connect4].

### MDP Abstraction
```kotlin
abstract class MDP<StateType, ActionType> {
    abstract fun transition(StateType, ActionType) : StateType
    abstract fun reward(StateType, ActionType?, StateType) : Double
    abstract fun initialState() : StateType
    abstract fun isTerminal(StateType) : Boolean
    abstract fun actions(StateType) : Collection<ActionType>
```

## Solver Design

*mctreesearch4j* provides a default implementation known as ``GenericSolver``, and an alternate ``StatefulSolver``. The abstract ``Solver`` serves as the base class for both versions, and defines a set of functionalities that all solver implementations must provide as well as a set of extensible utilities. Similar to the MDP abstraction, the solver uses a set of type parameters to provide strongly typed methods that unify and simplify the solver implementation. The programmer is free to select the data type or data structure that best defines how states and actions are represented in their MDP domain.

The difference between solvers lies in their respective memory utilization of abstract node types to track states during MCTS iterations. The default ``GenericSolver`` provides a leaner implementation, where actions are tracked and no explicit states are stored permanently. The states tracked with ``GenericSolver`` are dynamic and the MDP transitions must be rerun when traversing the search tree during selection in order to infer the state. The ``StatefulSolver`` keeps an explicit record of the visited states, and additional information, such as terminality and availability of downstream actions. The extra overhead of storing the state explicitly in the MCTS node, allows the algorithm to optimize its search using information from previously visited states. This is particularly useful for deterministic games, where a re-computation of the MDP transition is not necessary to determine the state of the agent after taking a specific action. This results in different implementations of the *Selection* step, while maintaining identical implementations of *Expansion*, *Simulate* and *Backpropagation*.

![Software Architecture.](software_design_mcts.png?raw=true "Title")


# Customization

*mctreesearch4j* maintains a high degree of modularity. For example, the key mechanisms highlighted in Fig. \ref{fig:mcts-diagram} can be rearranged, allowing scientists to experiment with different versions of MCTS. Furthermore, the internal functionality of each of the key mechanisms can be redefined. And complex meta heuristics and search algorithms can be implemented to add specific advantages to the solver for specific MDP domains. All of these features are important to give researchers the flexibility to design bespoke MCTS algorithms for their respective research.

Though the default MCTS implementation works well in many scenarios, there are situations where knowledge about specific problem domains can be applied to improve the MCTS performance. Improvements to MCTS, such as heuristics driven simulation, exploit domain knowledge to improve solver performance. We demonstrate that a Reversi AI that uses heuristics derived from [@Guenther:2004] is able to outperform the basic MCTS implementation contained in *mctreesearch4j*. These heuristics are programmed via extensibility points in the *mctreesearch4j* solver implementation, where the key mechanisms can be altered or augmented. In our Heuristic Implementation Example, we introduce the ``heuristicWeight`` array, a 2D array storing domain specific ratings of every position on a Reversi board representing the desirability of that position on the board. The negative numbers represent a likely loss and positive numbers representing a likely win, again as represented in Fig. \ref{fig:reversi-heu}.

This value is taken into consideration when traversing down the simulation tree. The ``heuristicWeight`` array adjusts the propensity to explore any position for both agents based on the heurisitc's belief about the desirability of the position. To alter the MCTS simulation phase we override the ``simulate()`` method and create a new definition for it. The application of this ``heuristicWeight`` only requires minor alterations to the ``simulate()`` method, as illustrated in *Heuristic Implementation Example*.

### Heuristic Implementation Example

```kotlin
override fun simulate(node: NodeType): Double {
    /*... Original Simulation code ...*/
    while(true) {
        val validActions = mdp.actions(currentState)
        var bestActionScore = Int.MIN_VALUE // Begin heuristic
            var bestActions = mutableListOf<Point>()
            for (action in validActions) {
                val score = heuristicWeight[action.x][action.y]
                if (score > bestActionScore) {
                    bestActionScore = score
                    bestActions = mutableListOf(action)
                }
                if (score == bestActionScore) {bestActions.add(action)}
            }
        val randomAction = bestActions.random() // End heuristic
        val newState = mdp.transition(currentState, randomAction)
        /*... Original Simulation code ...*/
    }
}
```

![Reversi Heuristic. \label{fig:reversi-heu} ](reversi-heu.png?raw=true "Title")

# Results

When the MCTS solver is accurately selecting the optimal solutions, it will continue to compel the agent to explore in the optimal subset of the state space, and reduce its exploration in the non-optimal subset of the state space. We provide a quick example in the MDP Domain of GridWorld, detailed in [@Liu:2021]. The cumulative number of visits corresponding to the optimal policy is proportionally increasing with respect to the number of MCTS iterations. Whereas for non-optimal solutions, the cumulative visits are significantly lower because the solver will avoid visiting the non-optimal state subspace.

![Convergence of visits.](gw_visits.png?raw=true "Title")

# Real World Applications

*mctreesearch4j* is designed to run on the JVM, and therefore is useful in a research setting where JVM is intended to be used alongside either the implementation of the game or the core AI engine. Typically, for state-of-the-art AI, MCTS is used in combination with other Deep Learning algorithms, where MCTS drives the exploration of hypothetical games to be played out in order for the agent to make better decisions.

MCTS is typically combined with Deep Learning strategies to achieve state-of-the-art performance for complex games whose solutions are typically intractable via standard tree-search methods. For example in AlphaGo [silver:2016], MCTS is combined with two types of neural networks, namely value networks and policy networks, to achieve state-of-the-art performance in the game of Go, never achieved before. Similarly, in the area of control theory, MCTS acting in combination with Value and Policy Networks to produce state and action space distributions, have shown impressive performance in probabilistic planning for self-driving cars, as demonstrated in Tesla’s AI day event [tesla:2021]. The utility of MCTS illustrates the acceptance of MCTS into the standard repetoire of both probablisitc and deterministic planning. And as the need for MCTS grows among industry and research, the impetus to create MCTS software libraries that are standardized, reliable, and can integrate well into general software ecosystems, such as JVM, is increasing.

One key industry application of *mctreesearch4j* is the fact that it is a JVM software, which allows it to integrate into both legacy and modern JVM software. One key application can be in mobile gaming, which has a large developer base in Android, largely written in Kotlin within the JVM. *mctreesearch4j*’s lightweight design with guaranteed reliability and JVM compatibility, tested extensively on both Java 8 and Java 11, make it a prime candidate for robust AI in turn-based games on mobile devices.

# Future Contributions

Future work on *mctreesearch4j* could potentially involve further optimization, particularly of the memory consumption of the algorithm. Additional speed improvements can be improved via multi-threading and parallelization, which are fully supported on the JVM.

The fact that *mctreesearch4j* is based in the JVM can however also serve as a downside, since Java based programming languages are not typically used for research problems within AI. That, combined with JVM language’s added complexity when compared to scripting languages such as Python, can result in *mctreesearch4j* experiencing a slow adoption in a pure research setting. Nevertheless, the key design principles that *mctreesearch4j* was designed on, could easily be ported into another scripting based language, more commonly found in AI research, namely Python. *mctreesearch4j* is designed to serve as a foundational template for research in MCTS for the JVM, and it can also serve as a reliable importable library for JVM applications which implement MCTS.

# Conclusion

In closing, *mctreesearch4j* presents a framework which enables programmers to adapt an MCTS solver to a variety of MDP domains. This is important because software application was a main focus of *mctreesearch4j*. Furthermore, *mctreesearch4j* is fully compatible with JVM, and this design decision was made due to the excellent support of class structure and generic variable typing in Kotlin, and other JVM languages, as well as support for mobile applications. Yet most importantly, *mctreesearch4j* is modular and extensible, the key mechanism of MCTS are broken down, and the programmer is able to inherit class members, redefine and/or re-implement certain sections of the algorithm while maintaining a high degree of MCTS standardization.

# Acknowledgements

This research was conducted without external funding, resulting in no conflicts of interests. We would like to acknowledge the Maven Central Repository, for providing an important service to make our JVM artifacts accessible to all JVM projects.

*mctreesearch4j* is available here: [mvnrepository.com/artifact/ca.aqtech/mctreesearch4j](https://mvnrepository.com/artifact/ca.aqtech/mctreesearch4j)


# References