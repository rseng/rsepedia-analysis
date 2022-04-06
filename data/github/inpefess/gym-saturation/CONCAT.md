---
title: 'gym-saturation: an OpenAI Gym environment for saturation provers'
tags:
  - Python
  - OpenAI Gym
  - automated theorem prover
  - saturation prover
  - reinforcement learning
authors:
  - name: Boris Shminke
    orcid: 0000-0002-1291-9896
    affiliation: 1
affiliations:
 - name: Laboratoire J.A. Dieudonné, CNRS and Université Côte d'Azur, France
   index: 1
date: 1 October 2021
bibliography: paper.bib
# Summary
---

`gym-saturation` is an OpenAI Gym [@DBLP:journals/corr/BrockmanCPSSTZ16] environment for reinforcement learning (RL) agents capable of proving theorems. Currently, only theorems written in a formal language of the Thousands of Problems for Theorem Provers (TPTP) library [@Sut17] in clausal normal form (CNF) are supported. `gym-saturation` implements the 'given clause' algorithm (similar to the one used in Vampire [@DBLP:conf/cav/KovacsV13] and E Prover [@DBLP:conf/cade/0001CV19]). Being written in Python, `gym-saturation` was inspired by PyRes [@DBLP:conf/cade/0001P20]. In contrast to the monolithic architecture of a typical Automated Theorem Prover (ATP), `gym-saturation` gives different agents opportunities to select clauses themselves and train from their experience. Combined with a particular agent, `gym-saturation` can work as an ATP. Even with a non trained agent based on heuristics, `gym-saturation` can find refutations for 688 (of 8257) CNF problems from TPTP v7.5.0.

# Statement of need

Current applications of RL to saturation-based ATPs like Enigma [@DBLP:conf/cade/JakubuvCOP0U20] or Deepire [@DBLP:conf/cade/000121a] are similar in that the environment and the agent are not separate pieces of software but parts of larger systems that are hard to disentangle. The same is true for non saturation-based RL-friendly provers too (e.g. lazyCoP, @DBLP:conf/tableaux/RawsonR21). This monolithic approach hinders free experimentation with novel machine learning (ML) models and RL algorithms and creates unnecessary complications for ML and RL experts willing to contribute to the field. In contrast, for interactive theorem provers, projects like HOList [@DBLP:conf/icml/BansalLRSW19] or GamePad [@DBLP:conf/iclr/HuangDSS19] separate the concepts of environment and agent. Such modular architecture may lead to the development of easily comparable agents based on diverse approaches (see, e.g. @DBLP:conf/aaai/PaliwalLRBS20 or @DBLP:journals/corr/abs-1905-10501). `gym-saturation` is an attempt to implement a modular environment-agent architecture of an RL-based ATP. In addition, some RL empowered saturation ATPs are not accompanied with their source code [@9669114], while `gym-saturation` is open-source software.

# Usage example

Suppose we want to prove an extremely simple theorem with a very basic agent. We can do that in the following way:

```python
# first we create and reset a OpenAI Gym environment
from importlib.resources import files
import gym

env = gym.make(
    "gym_saturation:saturation-v0",
    # we will try to find a proof shorter than 10 steps
    step_limit=10,
    # for a classical syllogism about Socrates
    problem_list=[
        files("gym_saturation").joinpath(
            "resources/TPTP-mock/Problems/TST/TST003-1.p"
        )
    ],
)
env.reset()
# we can render the environment (that will become the beginning of the proof)
print("starting hypotheses:")
print(env.render("human"))
# our 'age' agent will always select clauses for inference
# in the order they appeared in current proof attempt
action = 0
done = False
while not done:
    observation, reward, done, info = env.step(action)
    action += 1
# SaturationEnv has an additional method
# for extracting only clauses which became parts of the proof
# (some steps were unnecessary to find the proof)
print("refutation proof:")
print(env.tstp_proof)
print(f"number of attempted steps: {action}")
```

The output of this script includes a refutation proof found:

```
starting hypotheses:
cnf(p_imp_q, hypothesis, ~man(X0) | mortal(X0)).
cnf(p, hypothesis, man(socrates)).
cnf(q, hypothesis, ~mortal(socrates)).
refutation proof:
cnf(_0, hypothesis, mortal(socrates), inference(resolution, [], [p_imp_q, p])).
cnf(_2, hypothesis, $false, inference(resolution, [], [q, _0])).
number of attempted steps: 6
```

# Architecture

`gym-saturation` includes several sub-packages:

* parsing (happens during `env.reset()` in example code snippet)
* logic operations (happen during `env.step(action)` in the example)
* AI Gym environment implementation
* agent testing (a bit more elaborated version of the `while` loop from the examle)

`gym-saturation` relies on a deduction system of four rules which is known to be refutationally complete [@doi:10.1137/0204036]:

\begin{align*}
{\frac{C_1\vee A_1,C_2\vee\neg A_2}{\sigma\left(C_1\vee C_2\right)}},\sigma=mgu\left(A_1,A_2\right)\quad\text{(resolution)}
\end{align*}
\begin{align*}
{\frac{C_1\vee s\approx t,C_2\vee L\left[r\right]}{\sigma\left(L\left[t\right]\vee C_1\vee C_2\right)}},\sigma=mgu\left(s,r\right)\quad\text{(paramodulation)}
\end{align*}
\begin{align*}
{\frac{C\vee A_1\vee A_2}{\sigma\left(C\vee L_1\right)}},\sigma=mgu\left(A_1,A_2\right)\quad\text{(factoring)}
\end{align*}
\begin{align*}
\frac{C\vee s\not\approx t}{\sigma\left(C\right)},\sigma=mgu\left(s,t\right)\quad\text{(reflexivity resolution)}
\end{align*}

where $C,C_1,C_2$ are clauses, $A_1,A_2$ are atomic formulae, $L$ is a literal, $r,s,t$ are terms, and $\sigma$ is a substitution (most general unifier). $L\left[t\right]$ is a result of substituting the term $t$ in $L\left[r\right]$ for the term $r$ at only one chosen position.

For parsing, we use the LARK parser [@LARK]. We represent the clauses as Python classes forming tree-like structures. `gym-saturation` also includes a JSON serializer/deserializer for those trees. For example, a TPTP clause

```
cnf(a2,hypothesis,
    ( ~ q(a) | f(X) = X )).
``` 
becomes

```python
Clause(
	literals=[
		Literal(
			negated=True,
			atom=Predicate(
				name="q", arguments=[Function(name="a", arguments=[])]
			),
		),
		Literal(
			negated=False,
			atom=Predicate(
				name="=",
				arguments=[
					Function(name="f", arguments=[Variable(name="X")]),
					Variable(name="X"),
				],
			),
		),
	],
	label="a2",
)
```

This grammar serves as the glue for `gym-saturation` sub-packages, which are, in principle, independent of each other. After switching to another parser or another deduction system, the agent testing script won't break, and RL developers won't need to modify their agents for compatibility (for them, the environment will have the same standard OpenAI Gym API).

![A diagram showing interactions between four main subpackages of `gym-saturation`: 1) parsing; 2) logic operations (including the given clause algorithm); 3) OpenAI Gym Env implementation; 4) the agent testing script.\label{fig:architecture}](architecture.png)

Agent testing is a simple episode pipeline (see \autoref{fig:architecture}). It is supposed to be run in parallel (e.g. using GNU Parallel, @tange_2021_5233953) for a testing subset of problems. See the following table for the testing results of two popular heuristic-based agents on TPTP v7.5.0 (trained RL agents should strive to be more successful than those primitive baselines):

| | __size agent__ | __age agent__ | __size&age agent__ |
|-|-|-|-|
| __proof found__ | 509 | 206 | 688 |
| __step limit__ | 1385 | 35 | 223 |
| __out of memory__ | 148 | 149 | 148 |
| __5 min time out__ | 6215 | 7867 | 7198 |
| __total__ | 8257 | 8257 | 8257 |

`size agent` is an agent which always selects the shortest clause.
     
`age agent` is an agent which always selects the clause which arrived first to the set of unprocessed clauses ('the oldest one').
     
`size&age agent` is an agent which selects the shortest clause five times in a row and then one time --- the oldest one.
     
'Step limit' means an agent didn't find proof after 1000 steps (the longest proof found consists of 287 steps). This can work as a 'soft timeout'.

# Mentions

At the moment of writing this paper, `gym-saturation` was used by its author during their PhD studies for creating experimental RL-based ATPs.

# Acknowledgements

This work has been supported by the French government, through the 3IA Côte d'Azur Investments in the Future project managed by the National Research Agency (ANR) with the reference number ANR-19-P3IA-0002. This work was performed using HPC resources from GENCI-IDRIS (Grant 2021-AD011013125).

# References
..
  Copyright 2021-2022 Boris Shminke

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

|PyPI version| |Anaconda| |CircleCI| |Documentation Status| |codecov| |Binder| |JOSS|

gym-saturation
==============

``gym-saturation`` is an `OpenAI Gym <https://gym.openai.com/>`__
environment for reinforcement learning (RL) agents capable of proving
theorems. Currently, only theorems written in `TPTP
library <http://tptp.org>`__ formal language in clausal normal form
(CNF) are supported. ``gym-saturation`` implements the ‘given clause’
algorithm (similar to one used in
`Vampire <https://github.com/vprover/vampire>`__ and `E
Prover <https://github.com/eprover/eprover>`__). Being written in
Python, ``gym-saturation`` was inspired by
`PyRes <https://github.com/eprover/PyRes>`__. In contrast to monolithic
architecture of a typical Automated Theorem Prover (ATP),
``gym-saturation`` gives different agents opportunities to select
clauses themselves and train from their experience. Combined with a
particular agent, ``gym-saturation`` can work as an ATP.

``gym-saturation`` can be interesting for RL practitioners willing to
apply their experience to theorem proving without coding all the
logic-related stuff themselves. It also can be useful for automated
deduction researchers who want to create an RL-empowered ATP.

How to Install
==============

The best way to install this package is to use ``pip``:

.. code:: sh

   pip install gym-saturation

Another option is to use ``conda``:

.. code:: sh

   conda install -c conda-forge gym-saturation
   
One can also run it in a Docker container:

.. code:: sh

   docker build -t gym-saturation https://github.com/inpefess/gym-saturation.git
   docker run -it --rm -p 8888:8888 gym-saturation jupyter-lab --ip=0.0.0.0 --port=8888 --no-browser

How to use
==========

See `the
notebook <https://github.com/inpefess/gym-saturation/blob/master/examples/example.ipynb>`__
or run it in
`Binder <https://mybinder.org/v2/gh/inpefess/gym-saturation/HEAD?labpath=example.ipynb>`__
for more information.

How to Contribute
=================

`Pull requests <https://github.com/inpefess/gym-saturation/pulls>`__ are
welcome. To start:

.. code:: sh

   git clone https://github.com/inpefess/gym-saturation
   cd gym-saturation
   # activate python virtual environment with Python 3.7+
   pip install -U pip
   pip install -U setuptools wheel poetry
   poetry install
   # recommended but not necessary
   pre-commit install

All the tests in this package are
`doctests <https://docs.python.org/3/library/doctest.html>`__. One can
run them with the following command:

.. code:: sh

   pytest --doctest-modules gym-saturation

To check the code quality before creating a pull request, one might run
the script ``local-build.sh``. It locally does nearly the same as the CI
pipeline after the PR is created.

Reporting issues or problems with the software
==============================================

Questions and bug reports are welcome on `the
tracker <https://github.com/inpefess/gym-saturation/issues>`__.

More documentation
==================

More documentation can be found
`here <https://gym-saturation.readthedocs.io/en/latest>`__.

.. |PyPI version| image:: https://badge.fury.io/py/gym-saturation.svg
   :target: https://badge.fury.io/py/gym-saturation
.. |CircleCI| image:: https://circleci.com/gh/inpefess/gym-saturation.svg?style=svg
   :target: https://circleci.com/gh/inpefess/gym-saturation
.. |Documentation Status| image:: https://readthedocs.org/projects/gym-saturation/badge/?version=latest
   :target: https://gym-saturation.readthedocs.io/en/latest/?badge=latest
.. |codecov| image:: https://codecov.io/gh/inpefess/gym-saturation/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/inpefess/gym-saturation
.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/inpefess/gym-saturation/HEAD?labpath=example.ipynb
.. |JOSS| image:: https://joss.theoj.org/papers/c4f36ec7331a0dde54d8c3808fbff9c3/status.svg
   :target: https://joss.theoj.org/papers/c4f36ec7331a0dde54d8c3808fbff9c3
.. |Anaconda| image:: https://anaconda.org/conda-forge/gym-saturation/badges/version.svg
   :target: https://anaconda.org/conda-forge/gym-saturation
..
  Copyright 2021-2022 Boris Shminke

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

#################
Testing an Agent
#################

Suppose you already have a trained agent implemented. Then you can use `an agent testing script`_ in parallel like that::

  find $TPTP_HOME/Problems/*/*-*.p | parallel --bar --jobs 80% --timeout 30000% python agent_testing.py --problem_file {} --output_folder TPTP_CNF --step_limit 20

Or you can use `Slurm Workload Manager <https://slurm.schedmd.com/>`__. See `an example <https://github.com/inpefess/gym-saturation/tree/master/slurm-jobs>`__ from the project's repo.
  
You can write your own agent testing script based on ``agent_testing.py`` by calling ``episode`` function with your agent as an argument.

After processing problems, you can get a report about its performance::

  import os
  from gym_saturation.agent_testing import agent_testing_report
  from glob import glob
  import sys

  sys.setrecursionlimit(10000)
  problem_list = sorted(glob(os.path.join(
      os.environ["TPTP_HOME"], "Problems", "*", "*-*.p")
  ))
  report = agent_testing_report(problem_list, "TPTP_CNF_20")

Applying different policies ``gym-saturation`` leads to the following results on all CNF problems from TPTP-v7.5.0:

.. list-table:: Numbers of problems
   :header-rows: 1

   * - 
     - size agent
     - age agent
     - size&age agent
   * - **proof found**
     - 509
     - 206
     - 688
   * - **step limit**
     - 1385
     - 35
     - 223
   * - **out of memory**
     - 148
     - 149
     - 148
   * - **5 min time out**
     - 6215
     - 7867
     - 7198
   * - **total**
     - 8257
     - 8257
     - 8257

:ref:`Size agent<size_agent>` is an agent which always selects the shortest clause.
     
:ref:`Age agent<age_agent>` is an agent which always selects the clause which arrived first to the set of unprocessed clauses ('the oldest one').
     
:ref:`Size&age agent<size_age_agent>` is an agent which selects the shortest clause five times in a row and then one time --- the oldest one.
     
'Step limit' means an agent didn't find proof after 1000 steps (the longest proof found consists of 287 steps). This can work as a 'soft timeout'.

.. _an agent testing script: https://github.com/inpefess/gym-saturation/tree/master/gym_saturation/agent_testing.py
..
  Copyright 2021-2022 Boris Shminke

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

.. include:: ../../README.rst
    
.. toctree::
   :maxdepth: 2
   :caption: Contents:
	     
   what-is-going-on
   testing-an-agent
   package-documentation

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
  
..
  Copyright 2021-2022 Boris Shminke

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

#################  
What is going on
#################

One can write theorems in a machine-readable form. This package uses the `CNF`_ sublanguage of `TPTP`_. Before using the environment, you will need to download a recent TPTP archive (ca 600MB).

A statement of a theorem becomes a list of clauses. In a given clause algorithm, one divides the clauses in processed and not processed yet. Then at each step, one selects a not processed yet clause as a given clause. If it's empty (we arrived at a contradiction, i.e. found a refutation proof), the algorithm stops with success. If not, one applies all possible deduction rules to the given clause and all processed clauses. Then we add deduction results to the unprocessed set, and the given clause goes into the processed. The algorithm iterates if we didn't run out of time and unprocessed clauses.

The deduction rules are the following (this deductive system is known to be refutation complete):

* :ref:`resolution <resolution>`
* :ref:`factoring <factoring>`
* :ref:`paramodulation <paramodulation>`
* :ref:`reflexivity resolution <reflexivity_resolution>`

For the choice of a given clause, one usually employs a clever combination of heuristics. Of course, we can reformulate the same process as a reinforcement learning task.

What is a State
****************

(More or less resembles `ProofState class of PyRes`_)

The environment's state is a list of logical clauses. Each clause is a list of literals and also has several :ref:`properties <Clause>`.

Literal is a predicate, negated or not. A predicate can have arguments, which can be functions or variables. Functions can have arguments, which in turn can be functions or variables.

Grammar is encoded in Python objects in a self-explanatory way. Each grammar object is a dictionary with an obligatory key ``class`` (:ref:`Clause <Clause>`, :ref:`Literal <Literal>`, :ref:`Predicate <Predicate>`, :ref:`Function <Function>`, :ref:`Variable <Variable>`), and other keys representing this object's properties (such as being negated or having a list of arguments).

What is an Observation
***********************

An observation visible by an agent is a Python dictionary having two keys: `action_mask` and `real_obs`. Action mask is a `numpy` array of zeros and ones of some fixed length. A user can change a default value (100000) for this length by passing a `max_clauses` argument to the environment constructor. If at some step there are more than `max_clauses` clauses in the state, the environment returns ``done == True``. For any index in `action_mask`, if there is no clause with such an index in the state, the mask value is zero. It's also zero if the clause is marked as processed. For the indices of the clauses available to become a so-called 'given clause', the mask equals one.

`real_obs` is the state (a list of clauses). Since in OpenAI Gym observations have to live in some pre-defined space, there is a OpenAI compatible :ref:`space class<clause_space>` for a list of clauses.

What is an Action
******************

Action is an index of a clause from the state. Valid actions are only indices of not processed clauses.

What is a Reward
*****************

``1.0`` if the proof is found (a clause with an empty list of literals is selected as an action).

``0.0`` otherwise

Important notice
*****************

Usually, saturation provers use a timeout in seconds since they work in real-time mode. Here, we live in a discrete time, so we limit a prover by the number of saturation algorithm steps taken, not wall-clock time.

.. _CNF: https://en.wikipedia.org/wiki/Clausal_normal_form
.. _TPTP: http://www.tptp.org/
.. _ProofState class of PyRes: https://github.com/eprover/PyRes/blob/master/saturation.py
.. _resolution: https://en.wikipedia.org/wiki/Resolution_(logic)#Resolution_in_first_order_logic
.. _factoring: https://en.wikipedia.org/wiki/Resolution_(logic)#Factoring
.. _paramodulation: https://en.wikipedia.org/wiki/Resolution_(logic)#Paramodulation
..
  Copyright 2021-2022 Boris Shminke

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

######################
Package Documentation
######################

Environment and Agent Testing
******************************
.. automodule:: gym_saturation.envs.saturation_env
   :members:
.. automodule:: gym_saturation.agent_testing
   :members:
.. automodule:: gym_saturation.clause_space
   :members:

Common modules
***************

.. automodule:: gym_saturation.grammar
   :members:

Logical Operations
*******************

.. automodule:: gym_saturation.logic_ops.factoring
   :members:
.. automodule:: gym_saturation.logic_ops.paramodulation
   :members:
.. automodule:: gym_saturation.logic_ops.reflexivity_resolution
   :members:
.. automodule:: gym_saturation.logic_ops.resolution
   :members:
.. automodule:: gym_saturation.logic_ops.substitution
   :members:
.. automodule:: gym_saturation.logic_ops.unification
   :members:
.. automodule:: gym_saturation.logic_ops.utils
   :members:

Parsing
********

.. automodule:: gym_saturation.parsing.cnf_parser
   :members:
.. automodule:: gym_saturation.parsing.tptp_parser
   :members:
