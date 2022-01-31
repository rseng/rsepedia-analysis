# Bluesky — An Experiment Specification & Orchestration Engine

[![Build Status](https://img.shields.io/github/workflow/status/bluesky/bluesky/Unit%20Tests)](https://github.com/bluesky/bluesky/actions?query=workflow%3A%22Unit+Tests%22+branch%3Amaster)
[![PyPI](https://img.shields.io/pypi/v/bluesky)](https://pypi.org/project/bluesky/)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/bluesky)](https://anaconda.org/conda-forge/bluesky)

The Bluesky Python Package is an experiment specification and orchestration engine. 
- Specify the logic of an experiment in a high-level, hardware-abstracted way.
- First-class support for adaptive feedback between analysis and acquisition.
- Data is emitted in a streaming fashion in standard Python data structures.
- Pause/resume, robust error handling, and rich metadata capture are built in.

[**Bluesky Documentation**](http://blueskyproject.io/bluesky).

The Bluesky Project enables experimental science at the lab-bench or facility scale. It is a collection of Python libraries that are co-developed but independently useful and may be adopted *a la carte*.

[**Bluesky Project Documentation**](http://blueskyproject.io).

<!--- Provide a general summary of the issue in the Title above -->

## Expected Behavior
<!--- If you're describing a bug, tell us what should happen -->
<!--- If you're suggesting a change/improvement, tell us how it should work -->

## Current Behavior
<!--- If describing a bug, tell us what happens instead of the expected behavior -->
<!--- If suggesting a change/improvement, explain the difference from current behavior -->

## Possible Solution
<!--- Not obligatory, but suggest a fix/reason for the bug, -->
<!--- or ideas how to implement the addition or change -->

## Steps to Reproduce (for bugs)
<!--- Provide a link to a live example, or an unambiguous set of steps to -->
<!--- reproduce this bug. Include code to reproduce, if relevant -->
1.
2.
3.

## Context
<!--- How has this issue affected you? What are you trying to accomplish? -->
<!--- Providing context helps us come up with a solution that is most useful in the real world -->

## Your Environment
<!--- Include as many relevant details about the environment you experienced the bug in -->
# Contributing

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Submit a ticket for your issue, assuming one does not already exist.
  * Clearly describe the issue including steps to reproduce when it is a bug.
  * Make sure you fill in the earliest version that you know has the issue.
* Fork the repository on GitHub


## Making Changes

* Create a topic branch from where you want to base your work.
  * This is usually the master branch.
  * Only target release branches if you are certain your fix must be on that
    branch.
  * To quickly create a topic branch based on master; `git checkout -b
    fix/master/my_contribution master`. Please avoid working directly on the
    `master` branch.
* Make commits of logical units.
* Check for unnecessary whitespace with `git diff --check` before committing.
* Make sure your commit messages are in the proper format (see below)
* Make sure you have added the necessary tests for your changes.
* Run _all_ the tests to assure nothing else was accidentally broken.

### Writing the commit message

Commit messages should be clear and follow a few basic rules. Example:

```
ENH: add functionality X to bluesky.<submodule>.

The first line of the commit message starts with a capitalized acronym
(options listed below) indicating what type of commit this is.  Then a blank
line, then more text if needed.  Lines shouldn't be longer than 72
characters.  If the commit is related to a ticket, indicate that with
"See #3456", "See ticket 3456", "Closes #3456" or similar.
```

Describing the motivation for a change, the nature of a bug for bug fixes 
or some details on what an enhancement does are also good to include in a 
commit message. Messages should be understandable without looking at the code 
changes. 

Standard acronyms to start the commit message with are:
```
API: an (incompatible) API change
BLD: change related to building numpy
BUG: bug fix
CI : continuous integration
DEP: deprecate something, or remove a deprecated object
DEV: development tool or utility
DOC: documentation
ENH: enhancement
MNT: maintenance commit (refactoring, typos, etc.)
REV: revert an earlier commit
STY: style fix (whitespace, PEP8)
TST: addition or modification of tests
REL: related to releases
```
## The Pull Request

* Now push to your fork
* Submit a [pull request](https://help.github.com/articles/using-pull-requests) to this branch. This is a start to the conversation.

At this point you're waiting on us. We like to at least comment on pull requests within three business days 
(and, typically, one business day). We may suggest some changes or improvements or alternatives.

Hints to make the integration of your changes easy (and happen faster):
- Keep your pull requests small
- Don't forget your unit tests
- All algorithms need documentation, don't forget the .rst file
- Don't take changes requests to change your code personally
<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Motivation and Context
<!--- Why is this change required? What problem does it solve? -->
<!--- If it fixes an open issue, please link to the issue here. -->

## How Has This Been Tested?
<!--- Please describe in detail how you tested your changes. -->
<!--- Include details of your testing environment, and the tests you ran to -->
<!--- see how your change affects other areas of the code, etc. -->

<!--
## Screenshots (if appropriate):
-->
============
Contributing
============

Getting Started
===============

* Make sure you have a GitHub account
* Submit a ticket for your issue, assuming one does not already exist.

   * Clearly describe the issue including steps to reproduce when it is a bug.
   * Make sure you fill in the earliest version that you know has the issue.
* Fork the repository on GitHub


Making Changes
==============

* Create a topic branch from where you want to base your work.

   * This is usually the master branch.
   * Only target release branches if you are certain your fix must be on that
     branch.
   * To quickly create a topic branch based on master; `git checkout -b
     fix/master/my_contribution master`. Please avoid working directly on the
     `master` branch.
* Make commits of logical units.
* Check for unnecessary whitespace with `git diff --check` before committing.
* Make sure your commit messages are in the proper format (see below)
* Make sure you have added the necessary tests for your changes.
* Run _all_ the tests to assure nothing else was accidentally broken.

Writing the commit message
--------------------------

Commit messages should be clear and follow a few basic rules. Example:

.. code-block::

   ENH: add functionality X to bluesky.<submodule>.

   The first line of the commit message starts with a capitalized acronym
   (options listed below) indicating what type of commit this is.  Then a blank
   line, then more text if needed.  Lines shouldn't be longer than 72
   characters.  If the commit is related to a ticket, indicate that with
   "See #3456", "See ticket 3456", "Closes #3456" or similar.

Describing the motivation for a change, the nature of a bug for bug fixes
or some details on what an enhancement does are also good to include in a
commit message. Messages should be understandable without looking at the code
changes.

Standard acronyms to start the commit message with are:

* API: an (incompatible) API change
* BLD: change related to building the project
* BUG: bug fix
* CI : continuous integration
* DEP: deprecate something, or remove a deprecated object
* DEV: development tool or utility
* DOC: documentation
* ENH: enhancement
* MNT: maintenance commit (refactoring, typos, etc.)
* REV: revert an earlier commit
* STY: style fix (whitespace, PEP8)
* TST: addition or modification of tests
* REL: related to releases

The Pull Request
----------------

* Now push to your fork
* Submit a `pull request <https://help.github.com/articles/using-pull-requests>`_

At this point you're waiting on us. We like to at least comment on pull requests within three business days (and, typically, one business day).
We may suggest some changes or improvements or alternatives.

Hints to make the integration of your changes easy (and happen faster):
* Keep your pull requests small
* Don't forget your unit tests
* All features need documentation, don't forget the .rst file
* Don't take changes requests to change your code personally

Releasing
=========

To release Bluesky:

* Make a pull request adding release notes to the API changes document:
  :file:`docs/source/api_changes.rst`.
  It might be helpful to `review the changes since last release <https://docs.github.com/en/github/administering-a-repository/releasing-projects-on-github/comparing-releases>`_.
* Once the changelog is merged, mint a new release.

   * tag: :code:`v<VERSION>`
   * title: :code:`bluesky v<VERSION>`
   * description: copy content from api_changes
* Bluesky will be released to PyPI via our `github action <https://github.com/bluesky/bluesky/blob/master/.github/workflows/python-publish.yml>`_.
* The conda-forge bot will make a PR to update the `bluesky-feedstock <https://github.com/conda-forge/bluesky-feedstock>`_.
Interruptions
*************

The RunEngine can be safely interrupted and resumed. All plans get this
feature "for free."

.. _pausing_interactively:

Pausing Interactively
=====================

.. note::

    Looking for a quick refresher on pausing, resuming, or aborting
    interactively? Skip to the :ref:`interactive_pause_summary`.

While the RunEngine is executing a plan, it captures SIGINT (Ctrl+C).

Pause Now: Ctrl+C twice
-----------------------

.. code-block:: python

    In [14]: RE(scan([det], motor, 1, 10, 10))
    Transient Scan ID: 2     Time: 2018/02/12 12:43:12
    Persistent Unique Scan ID: '33a16823-e214-4952-abdd-032a78b8478f'
    New stream: 'primary'
    +-----------+------------+------------+------------+
    |   seq_num |       time |      motor |        det |
    +-----------+------------+------------+------------+
    |         1 | 12:43:13.3 |      1.000 |      0.607 |
    |         2 | 12:43:14.3 |      2.000 |      0.135 |
    |         3 | 12:43:15.3 |      3.000 |      0.011 |
    ^C
    A 'deferred pause' has been requested. The RunEngine will pause at the next checkpoint. To pause immediately, hit Ctrl+C again in the next 10 seconds.
    Deferred pause acknowledged. Continuing to checkpoint.
    ^C
    Pausing...
    ---------------------------------------------------------------------------
    RunEngineInterrupted                      Traceback (most recent call last)
    <ipython-input-14-826ee9dfb918> in <module>()
    ----> 1 RE(scan([det], motor, 1, 10, 10))

    ~/Documents/Repos/bluesky/bluesky/run_engine.py in __call__(self, *args, **metadata_kw)
        670
        671             if self._interrupted:
    --> 672                 raise RunEngineInterrupted(self.pause_msg) from None
        673
        674         return tuple(self._run_start_uids)

    RunEngineInterrupted:
    Your RunEngine is entering a paused state. These are your options for changing
    the state of the RunEngine:

    RE.resume()    Resume the plan.
    RE.abort()     Perform cleanup, then kill plan. Mark exit_stats='aborted'.
    RE.stop()      Perform cleanup, then kill plan. Mark exit_status='success'.
    RE.halt()      Emergency Stop: Do not perform cleanup --- just stop.

Before returning the prompt to the user, the RunEngine ensures that all motors
that it has touched are stopped. It also performs any device-specific cleanup
defined in the device's (optional) ``pause()`` method.

If execution is later resumed, the RunEngine will "rewind" through the plan to
the most recent :ref:`checkpoint <checkpoints>`, the last safe place to restart.

Pause Soon: Ctrl+C once
-----------------------

Pause at the next :ref:`checkpoint <checkpoints>`: typically, the next step in
a step scan. We call this "deferred pause." It avoids having to repeat any work
when the plan is resumed.

Notice that this time when Ctrl+C (^C) is hit, the current step (4) is allowed
to complete before execution is paused.

.. code-block:: python

    In [12]: RE(scan([det], motor, 1, 10, 10))
    Transient Scan ID: 1     Time: 2018/02/12 12:40:36
    Persistent Unique Scan ID: 'c5db9bb4-fb7f-49f4-948b-72fb716d1f67'
    New stream: 'primary'
    +-----------+------------+------------+------------+
    |   seq_num |       time |      motor |        det |
    +-----------+------------+------------+------------+
    |         1 | 12:40:37.6 |      1.000 |      0.607 |
    |         2 | 12:40:38.7 |      2.000 |      0.135 |
    |         3 | 12:40:39.7 |      3.000 |      0.011 |
    ^CA 'deferred pause' has been requested. The RunEngine will pause at the next checkpoint. To pause immediately, hit Ctrl+C again in the next 10 seconds.
    Deferred pause acknowledged. Continuing to checkpoint.
    |         4 | 12:40:40.7 |      4.000 |      0.000 |
    Pausing...
    ---------------------------------------------------------------------------
    RunEngineInterrupted                      Traceback (most recent call last)
    <ipython-input-12-826ee9dfb918> in <module>()
    ----> 1 RE(scan([det], motor, 1, 10, 10))

    ~/Documents/Repos/bluesky/bluesky/run_engine.py in __call__(self, *args, **metadata_kw)
        670
        671             if self._interrupted:
    --> 672                 raise RunEngineInterrupted(self.pause_msg) from None
        673
        674         return tuple(self._run_start_uids)

    RunEngineInterrupted:
    Your RunEngine is entering a paused state. These are your options for changing
    the state of the RunEngine:

    RE.resume()    Resume the plan.
    RE.abort()     Perform cleanup, then kill plan. Mark exit_stats='aborted'.
    RE.stop()      Perform cleanup, then kill plan. Mark exit_status='success'.
    RE.halt()      Emergency Stop: Do not perform cleanup --- just stop.

What to do after pausing
------------------------

After being paused, the RunEngine holds on to information that it might need in
order to resume later. It "knows" that it is in a paused state, and you can
check that at any time:

.. code-block:: python

    In [2]: RE.state
    Out[2]: 'paused'


During the pause, we can do anything: check readings, move motors, etc. It will
not allow you to execute a new plan until the current one is either resumed or
terminated. Your options are:

Resume
^^^^^^

.. code-block:: python

    In [3]: RE.resume()
    |         4 | 07:21:29.5 |     -5.714 |      0.000 |
    |         5 | 07:21:29.5 |     -4.286 |      0.000 |
    |         6 | 07:21:29.6 |     -2.857 |      0.017 |
    |         7 | 07:21:29.7 |     -1.429 |      0.360 |
    (etc.)

Depending on the plan, it may "rewind" to safely continue on and ensure all
data is collected correctly.

Abort
^^^^^

Allow the plan to perform any final cleanup. For example, some plans move
motors back to their starting positions. Mark the data as having been aborted,
so that this fact can be noted (if desired) in later analysis. All of the data
collected up this point will be saved regardless.

From a paused state:

.. code-block:: python

    In [3]: RE.abort()
    Aborting...
    Out[3]: ['8ef9388c-75d3-498c-a800-3b0bd24b88ed']

Stop
^^^^

``RE.stop()`` is functionally identical to ``RE.abort()``. The only
difference is that aborted runs are marked with ``exit_status: 'abort'``
instead of ``exit_status: 'success'``. This may be a useful distinction
during analysis.

Halt
^^^^

Aborting or stopping allows the plan to perform cleanup. We already mentioned
the example of a plan moving motors back to their starting positions at the
end.

In some situations, you may wish to prevent the plan from doing *anything*
--- you want to halt immediately, skipping cleanup. For this, use
``RE.halt()``.

.. _interactive_pause_summary:

Summary
-------

Interactively Interrupt Execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======================= ===========
Command                 Outcome
======================= ===========
Ctrl+C                  Pause soon.
Ctrl+C twice            Pause now.
======================= ===========
    
From a paused state
^^^^^^^^^^^^^^^^^^^

============== ===========
Command        Outcome
============== ===========
RE.resume()    Safely resume plan.
RE.abort()     Perform cleanup. Mark as aborted.
RE.stop()      Perform cleanup. Mark as success.
RE.halt()      Do not perform cleanup --- just stop.
RE.state       Check if 'paused' or 'idle'.
============== ===========

.. _suspenders:

Automated Suspension
====================

It can also be useful to interrupt execution automatically in response to some
condition (e.g., shutter closed, beam dumped, temperature exceeded some limit).
We use the word *suspension* to mean an unplanned pause initialized by some
agent running the background. The agent (a "suspender") monitors some condition
and, if it detects a problem, it suspends execution. When it detects that
conditions have returned to normal, it gives the RunEngine permission to resume
after some interval. This can operate unattended.

.. ipython::
    :verbatim:

    In [1]: RE(scan([det], motor, -10, 10, 15), LiveTable([motor, det]))
    +------------+-------------------+----------------+----------------+
    |   seq_num  |             time  |         motor  |           det  |
    +------------+-------------------+----------------+----------------+
    |         1  |  16:46:08.953815  |          0.03  |        290.00  |
    Suspending....To get prompt hit Ctrl-C to pause the scan
    |         2  |  16:46:20.868445  |          0.09  |        279.00  |
    |         3  |  16:46:29.077690  |          0.16  |        284.00  |
    |         4  |  16:46:33.540643  |          0.23  |        278.00  |
    +------------+-------------------+----------------+----------------+

A *suspended* plan does not return the prompt to the user. Like a paused plan,
it stops executing new instructions and rewinds to the most recent checkpoint.
But unlike a paused plan, it resumes execution automatically when conditions
return to normal.

To take manual control of a suspended plan, pause it by hitting Ctrl+C twice.
You will be given the prompt. When conditions are good again, you may manually
resume using ``RE.resume()``.

.. _installing_suspenders:

Installing Suspenders
---------------------

Bluesky includes several "suspenders" that work with ophyd Signals to monitor
conditions and suspend execution. It's also possible to write suspenders
from scratch to monitor anything at all.

We'll start with an example.

Example: Suspend a plan if the beam current dips low
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This defines a suspender and installs it on the RunEngine. With this, plans
will be automatically suspended when the ``beam_current`` signal goes below 2
and resumed once it exceeds 3.

.. code-block:: python

    from ophyd import EpicsSignal
    from bluesky.suspenders import SuspendFloor

    beam_current = EpicsSignal('...PV string...')
    sus = SuspendFloor(beam_current, 2, resume_thresh=3)
    RE.install_suspender(sus)

In the following example, the beam current dipped below 2 in the middle of
taking the second data point. It later recovered.

.. ipython::
    :verbatim:

    In [6]: RE(my_scan)
    +------------+-------------------+----------------+----------------+
    |   seq_num  |             time  |         theta  |    sclr_chan4  |
    +------------+-------------------+----------------+----------------+
    |         1  |  16:46:08.953815  |          0.03  |        290.00  |
    Suspending....To get prompt hit Ctrl-C to pause the scan
    |         2  |  16:46:20.868445  |          0.09  |        279.00  |
    |         3  |  16:46:29.077690  |          0.16  |        284.00  |
    |         4  |  16:46:33.540643  |          0.23  |        278.00  |
    +------------+-------------------+----------------+----------------+

Notice that the plan was suspended and then resumed. When it resumed, it went
back to the last checkpoint and re-took the second data point cleanly.

See the API documentation (follow the links in the table below) for other
suspender types and options, including a waiting period and cleanup
procedures to run pre-suspend and pre-resume.

Built-in Suspenders
-------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   bluesky.suspenders.SuspendBoolHigh
   bluesky.suspenders.SuspendBoolLow
   bluesky.suspenders.SuspendFloor
   bluesky.suspenders.SuspendCeil
   bluesky.suspenders.SuspendWhenOutsideBand
   bluesky.suspenders.SuspendWhenChanged

.. _checkpoints:

Checkpoints
===========

Plans are specified as a sequence of :ref:`messages <msg>`, granular
instructions like 'read' and 'set'. The messages can optionally include one
or more 'checkpoint' messages, indicating a place where it is safe to resume
after an interruption. For example, checkpoints are placed before each step of a
:func:`bluesky.plans.scan`.

Some experiments are not resumable: for example, the sample may be melting or
aging. Incorporating :func:`bluesky.plan_stubs.clear_checkpoint` in a plan
makes it un-resuming. If a pause or suspension are requested, the plan will
abort instead.

.. note::

    Some details about checkpoints and when they are allowed:

    It is not legal to create a checkpoint in the middle of a data point
    (between 'create' and 'save'). Checkpoints are implicitly created after
    actions that it is not safe to replay: staging a device, adding a
    monitor, or adding a subscription.

.. _planned_pauses:

Planned Pauses
==============

Pausing is typically done :ref:`interactively <pausing_interactively>` (Ctrl+C)
but it can also be incorporated into a plan. The plan can pause the RunEngine,
requiring the user to type ``RE.resume()`` to continue or ``RE.stop()``
(or similar) to clean up and stop.

.. code-block:: python

    import bluesky.plan_stubs as bps

    def pausing_plan():
        while True:
            yield from some_plan(...)
            print("Type RE.resume() to go again or RE.stop() to stop.")
            yield from bps.checkpoint()  # marking where to resume from
            yield from bps.pause()

Associated RunEngine Interface
==============================

State
-----

The RunEngine has a state machine defining its phases of operation and the
allowed transitions between them. As illustrated above, it can be inspected via
the ``state`` property.

The states are:

* ``'idle'``: RunEngine is waiting for instructions.
* ``'running'``: RunEngine is executing instructions.
* ``'paused'``: RunEngine is waiting for user input.

Suspender-related Methods
-------------------------

.. automethod:: bluesky.run_engine.RunEngine.install_suspender
    :noindex:

.. automethod:: bluesky.run_engine.RunEngine.remove_suspender
    :noindex:

.. automethod:: bluesky.run_engine.RunEngine.clear_suspenders
    :noindex:

The RunEngine also has a ``suspenders`` property, a collection of the
currently-installed suspenders.

Request Methods
---------------

This method is called when Ctrl+C is pressed or when a 'pause' Message is
processed. It can also be called by user-defined agents. See the next example.

.. automethod:: bluesky.run_engine.RunEngine.request_pause
    :noindex:

This method is used by the ``PVSuspend*`` classes above. It can also be called
by user-defined agents.

.. automethod:: bluesky.run_engine.RunEngine.request_suspend
    :noindex:


Example: Requesting a pause from the asyncio event loop
-------------------------------------------------------

Since the user does not have control of the prompt, calls to
``RE.request_pause`` must be planned in advance. Here is a example that pauses
the plan after 5 seconds.

.. code-block:: python

    from bluesky.plan_stubs import null

    def loop_forever():
        "a silly plan"
        while True:
            yield from null()

    import asyncio
    loop = asyncio.get_event_loop()
    # Request a pause 5 seconds from now.
    loop.call_later(5, RE.request_pause)

    # Execute the plan.
    RE(loop_forever())

    # Five seconds after ``call_later`` was run, the plan is paused.
    # Observe that the RunEngine is in a 'paused' state.
    RE.state

Above, we passed ``True`` to ``RE.request_pause`` to request a deferred pause.

Experimental: Record Interruptions
==================================

In the analysis stage, it can be useful to know if and when a run was
interrupted.  This experimental feature creates a special event stream
recording the time and nature of any interruptions.

.. warning::

    This is an experimental feature. It is tested but not yet widely used. It
    might be changed or removed in the future.

Activate this feature by setting

.. code-block:: python

    RE.record_interruptions = True

In this mode, the RunEngine emits a special event descriptor after opening a
new run. This name field in the descriptor is 'interruptions'. It has a single
data key:

.. code-block:: python

    {'interruptions': {'dtype': 'string',
                       'shape': None,
                       'source': 'RunEngine'}}

Each time the RunEngine is paused, suspended, or resumed during the run, an
Event document for that descriptor is created. The data payload
``event['data']['interruptions']`` is ``'pause'``, ``'suspend'``, or
``'resume'``. The associated time notes when the interruptions/resume was
processed.

To see this in action, try this example:

.. code-block:: python

    from bluesky.plans import count
    from bluesky.preprocessors import pchain
    from bluesky.plan_stubs import pause
    from ophyd.sim import det

    RE.record_interruptions = True

    RE(pchain(count([det]), pause(), count([det])), print)
    # ... RunEngine pauses
    RE.resume()

In the text that ``print`` dumps to the screen, look for the special
'interruptions' event descriptor and associated events.
===============================================
Translating Direct PyEpics Code to Bluesky Code
===============================================

.. warning:

    This section is still a work in progress.

How?
====

===========================   ======================================
interactive (blocking)        re-write for BlueSky plan()
===========================   ======================================
some.device.put("config")     yield from mv(some.device, "config")
motor.move(52)                yield from mv(motor, 52)
motor.velocity.put(5)         yield from mv(motor.velocity, 5)
===========================   ======================================


Why?
====
.. currentmodule:: bluesky.plans

Documents
=========

A primary design goal of bluesky is to enable better research by recording
rich metadata alongside measured data for use in later analysis. Documents are
how we do this.

A *document* is our term for a Python dictionary with a schema --- that is,
organized in a
`formally specified <https://github.com/NSLS-II/event-model>`_ way --- created
by the RunEngine during plan execution.  All of the metadata and data generated
by executing the plan is organized into documents.

A :doc:`later section <callbacks>` describes how outside functions can
"subscribe" to a stream of these documents, visualizing, processing, or saving
them. This section provides an outline of documents themselves, aiming to give
a sense of the structure and familiarity with useful components.

.. _run_overview:

Overview of a "Run"
-------------------

Each document belongs to a *run* --- loosely speaking, a dataset. Executing any
of the :ref:`built-in pre-assembled plans <preassembled_plans>`, like
:func:`scan` and :func:`count`, creates one run.

.. note::

    Fundamentally, the scope of a run is intentionally vague and flexible. One
    plan might generate many runs or one long run. It just depends on how you
    want to organize your data, both at collection time and analysis time.

    The tutorial's :ref:`tutorial_capture_data` section explores this.

The documents in each run are:

- A **Run Start document**, containg all of the metadata known at the start of
  the run. Highlights:

    - time --- the start time
    - plan_name --- e.g., ``'scan'`` or ``'count'``
    - uid --- unique ID that identifies this run
    - scan_id --- human-friendly integer scan ID (not necessarily unique)
    - any other :doc:`metadata captured at execution time <metadata>` from the
      plan or the user

- **Event documents**, containing the actual measurements. These are your data.

    - time --- a timestamp for this group of readings
    - seq_num --- sequence number, counting up from 1
    - data --- a dictionary of readings like
      ``{'temperature': 5.0, 'position': 3.0}``
    - timestamps --- a dictionary of individual timestamps for each reading,
      from the hardware

- **Event Descriptor documents** provide a schema for the data in the Event
  documents. They list all of the keys in the Event's data and give useful
  information about them, such as units and precision. They also contain
  information about the configuration of the hardware.

- A **Run Stop document**, containing metadata known only at the end of the
  run. Highlights:

    - time --- the time when the run was completed
    - exit_status --- "success", "abort", or "fail"

Every document has a ``time`` (its creation time) and a separate ``uid`` to
identify it. The Event documents also have a ``descriptor`` field linking them
to the Event Descriptor with their metadata. And the Event Descriptor and
Run Stop documents have a ``run_start`` field linking them to their Run
Start. Thus, all the documents in a run are linked back to the Run Start.

Documents in Detail
-------------------

Run Start
+++++++++

Again, a 'start' document marks the beginning of the run. It comprises
everything we know before we start taking data, including all metadata provided
by the user and the plan. (More on this in the :doc:`next section <metadata>`.)

All built-in plans provide some useful metadata like the names of the
detector(s) and motor(s) used. (User-defined plans may also do this; see
:ref:`this section <tutorial_plan_metadata>` of the tutorial.)

The command:

.. code-block:: python

    from bluesky.plans import scan
    from ophyd.sim import det, motor  # simulated detector, motor

    # Scan 'motor' from -3 to 3 in 10 steps, taking readings from 'det'.
    RE(scan([det], motor, -3, 3, 16), purpose='calibration',
       sample='kryptonite')

generates a 'start' document like this:

.. code-block:: python

    # 'start' document
    {'purpose': 'calibration',
     'sample': 'kryptonite',
     'detectors': ['det'],
     'motors': ['motor'],
     'plan_name': 'scan',
     'plan_type': 'generator',
     'plan_args': {'detectors': '[det]',
                   'motor': 'Mover(...)',
                   'num': '16',
                   'start': '-3',
                   'stop': '3'},
     'scan_id': 282,
     'time': 1442521005.6099606,
     'uid': '<randomly-generated unique ID>',
    }

.. note::

    Time is given in UNIX time (seconds since 1970). Software for looking at
    the data would, of course, translate that into a more human-readable form.

Event
+++++

An 'event' records one or more measurements with an associated time.

.. code-block:: python

    # 'event' document
    {'data':
        {'temperature': 5.0,
          'x_setpoint': 3.0,
          'x_readback': 3.05},
     'timestamps':
        {'temperature': 1442521007.9258342,
         'x_setpoint': 1442521007.5029348,
         'x_readback': 1442521007.5029348},
     'time': 1442521007.3438923,
     'seq_num': 1
     'uid': '<randomly-generated unique ID>',
     'descriptor': '<reference to a descriptor document>'}

From a data analysis perspective, these readings were simultaneous, but in
actuality the occurred at separate times.  The separate times of the individual
readings are not thrown away (they are recorded in 'timestamps') but the
overall event 'time' is often more useful.

Run Stop
++++++++

A 'stop' document marks the end of the run. It contains metadata that is not
known until the run completes.

The most commonly useful fields here are 'time' and 'exit_status'.

.. code-block:: python

    # 'stop' document
    {'exit_status': 'success',  # or 'fail' or 'abort'
     'reason': '',  # The RunEngine can provide reason for failure here.
     'time': 1442521012.1021606,
     'uid': '<randomly-generated unique ID>',
     'start': '<reference to the start document>',
     'num_events': {'primary': 16}
    }

Event Descriptor
++++++++++++++++

As stated above, a 'descriptor' document provides a schema for the data in the
Event documents. It provides useful information about each key in the data and
about the configuration of the hardware. The layout of a descriptor is detailed
and takes some time to cover, so we defer it to a
:doc:`later section <event_descriptors>`.
.. currentmodule:: bluesky.plans

Asynchronous Acquisition
========================

This section encompasses "fly scans," "monitoring," and in general handling
data acquisition that is occurring at different rates.

.. note::

    If you are here because you just want to "move two motors at once" or
    something in that category, you're in luck: you don't need anything as
    complex as what we present in this section. Read about multidimensional
    plans in the section on :doc:`plans`.

In short, "flying" is for acquisition at high rates and "monitoring" is for
acquisition at an irregular or slow rate. Monitoring does not guarantee that
all readings will be captured; i.e. monitoring is lossy. It is susceptible to
network glitches. But flying, by contract, is not lossy if correctly
implementated.

**Flying** means: "Let the hardware take control, cache data externally, and
then transfer all the data to the RunEngine at the end." This is essential when
the data acquisition rates are faster than the RunEngine or Python can go.

.. note::

    As a point of reference, the RunEngine processes message at a rate of
    about 35k/s (not including any time added by whatever the message *does*).


    .. code-block:: python

        In [3]: %timeit RE(Msg('null') for j in range(1000))
        10 loops, best of 3: 26.8 ms per loop

**Monitoring** a means acquiring readings whenever a new reading is available,
at a device's natural update rate. For example, we might monitor background
condition (e.g., beam current) on the side while executing the primary logic of
a plan. The documents are generated in real time --- not all at the end, like
flying --- so if the update rate is too high, monitoring can slow down the
execution of the plan. As mentioned above, monitoring is also lossy: if network
traffic is high, some readings may be missed.

Flying
------

In bluesky's view, there are three steps to "flying" a device during a scan.

1. **Kickoff**: Begin accumulating data. A 'kickoff' command completes once
   acquisition has successfully started.
2. **Complete**: This step tells the device, "I am ready whenever you are
   ready." If the device is just collecting until it is told to stop, it will
   report that it is ready immediately. If the device is executing some
   predetermined trajectory, it will finish before reporting ready.
3. **Collect**: Finally, the data accumulated by the device is transferred to
   the RunEngine and processed like any other data.

To "fly" one or more "flyable" devices during a plan, bluesky provides a
`preprocessor <preprocessors>`. It is available as a wrapper,
:func:`fly_during_wrapper`

.. code-block:: python

    from ophyd.sim import det, flyer1, flyer2  # simulated hardware
    from bluesky.plans import count
    from bluesky.preprocessors import fly_during_wrapper

    RE(fly_during_wrapper(count([det], num=5), [flyer1, flyer2]))

and as a decorator, :func:`fly_during_decorator`.

.. code-block:: python

    from ophyd.sim import det, flyer1, flyer2  # simulated hardware
    from bluesky.plans import count
    from bluesky.preprocessors import fly_during_wrapper

    # Define a new plan for future use.
    fly_and_count = fly_during_decorator([flyer1, flyer2])(count)

    RE(fly_and_count([det]))

Alternatively, if you are using :ref:`supplemental_data`, simply
append to or extend its list of flyers to kick off during every run:

.. code-block:: python

    from ophyd.sim import flyer1, flyer2

    # Assume sd is an instance of the SupplementalData set up as
    # descripted in the documentation linked above.
    sd.flyers.extend([flyer1, flyer2])

They will be included with all plans until removed.

.. _async_monitoring:

Monitoring
----------

To monitor some device during a plan, bluesky provides a
`preprocessor <preprocessors>`. It
is available as a wrapper, :func:`monitor_during_wrapper`

.. code-block:: python

    from ophyd.sim import det, det1
    from bluesky.plans import count
    from bluesky.preprocessors import monitor_during_wrapper

    # Record any updates from det1 while 'counting' det 5 times.
    RE(monitor_during_wrapper(count([det], num=5), [det1]))

and as a decorator, :func:`monitor_during_decorator`.

.. code-block:: python

    from ophyd.sim import det, det1
    from bluesky.plans import count
    from bluesky.preprocessors import monitor_during_wrapper

    # Define a new plan for future use.
    monitor_and_count = monitor_during_decorator([det1])(count)

    RE(monitor_and_count([det]))

Alternatively, if you are using :ref:`supplemental_data`, simply
append to or extend its list of signals to monitor:

.. code-block:: python

    from ophyd.sim import det1

    # Assume sd is an instance of the SupplementalData set up as
    # descripted in the documentation linked above.
    sd.monitors.append(det1)

They will be included with all plans until removed.
Live Visualization and Processing
*********************************

.. ipython:: python
   :suppress:
   :okwarning:

   from bluesky import RunEngine
   RE = RunEngine({})

.. _callbacks:

Overview of Callbacks
---------------------

As the RunEngine executes a plan, it organizes metadata and data into
*Documents,* Python dictionaries organized in a
:doc:`specified but flexible <documents>` way. Each time a new Document is
created, the RunEngine passes it to a list of functions. These functions can do
anything: store the data to disk, print a line of text to the screen, add a
point to a plot, or even transfer the data to a cluster for immediate
processing. These functions are called "callbacks."

We "subscribe" callbacks to the live stream of Documents coming from the
RunEngine. You can think of a callback as a self-addressed stamped envelope: it
tells the RunEngine, "When you create a Document, send it to this function for
processing."

Callback functions are run in a blocking fashion: data acquisition cannot
continue until they return. For light tasks like simple plotting or critical
tasks like sending the data to a long-term storage medium, this behavior is
desirable. It is easy to debug and it guarantees that critical errors will be
noticed immediately. But heavy computational tasks --- anything that takes more
than about 0.2 seconds to finish --- should be executed in a separate process
or server so that they do not hold up data acquisition. Bluesky provides nice
tooling for this use case --- see :ref:`zmq_callback`.

Simplest Working Example
------------------------

This example passes every Document to the ``print`` function, printing
each Document as it is generated during data collection.

.. code-block:: python

    from bluesky.plans import count
    from ophyd.sim import det

    RE(count([det]), print)

The ``print`` function is a blunt instrument; it dumps too much information to
the screen.  See :ref:`LiveTable <livetable>` below for a more refined option.

Ways to Invoke Callbacks
------------------------

Interactively
+++++++++++++

As in the simple example above, pass a second argument to the RunEngine.
For some callback function ``cb``, the usage is:

.. code-block:: python

    RE(plan(), cb))

A working example:

.. code-block:: python

    from ophyd.sim import det, motor
    from bluesky.plans import scan
    from bluesky.callbacks import LiveTable
    dets = [det]
    RE(scan(dets, motor, 1, 5, 5), LiveTable(dets))

A *list* of callbacks --- ``[cb1, cb2]`` --- is also accepted; see
:ref:`filtering`, below, for additional options.

Persistently
++++++++++++

The RunEngine keeps a list of callbacks to apply to *every* plan it executes.
For example, the callback that saves the data to a database is typically
invoked this way. For some callback function ``cb``, the usage is:

.. code-block:: python

    RE.subscribe(cb)

This step is usually performed in a startup file (i.e., IPython profile).

.. automethod:: bluesky.run_engine.RunEngine.subscribe
    :noindex:

.. automethod:: bluesky.run_engine.RunEngine.unsubscribe
    :noindex:

.. _subs_decorator:

Through a plan
++++++++++++++

Use the ``subs_decorator`` :ref:`plan preprocessor <preprocessors>` to attach
callbacks to a plan so that they are subscribed every time it is run.

In this example, we define a new plan, ``plan2``, that adds some callback
``cb`` to some existing plan, ``plan1``.

.. code-block:: python

    from bluesky.preprocessors import subs_decorator

    @subs_decorator(cb)
    def plan2():
        yield from plan1()

or, equivalently,

.. code-block:: python

    plan2 = subs_decorator(cb)(plan1)

For example, to define a variant of ``scan`` that includes a table by default:

.. code-block:: python

    from bluesky.plans import scan
    from bluesky.preprocessors import subs_decorator

    def my_scan(detectors, motor, start, stop, num, *, per_step=None, md=None):
        "This plan takes the same arguments as `scan`."

        table = LiveTable([motor] + list(detectors))

        @subs_decorator(table)
        def inner():
            yield from scan(detectors, motor, start, stop, num,
                            per_step=per_step, md=md)

        yield from inner()

Callbacks for Visualization & Fitting
-------------------------------------

.. _livetable:

LiveTable
+++++++++

As each data point is collected (i.e., as each Event Document is generated) a
row is added to the table. Demo:

.. ipython:: python

    from bluesky.plans import scan
    from ophyd.sim import motor, det
    from bluesky.callbacks import LiveTable

    RE(scan([det], motor, 1, 5, 5), LiveTable([motor, det]))

Pass an empty list of columns to show simply 'time' and 'seq_num' (sequence
number).

.. code-block:: python

    LiveTable([])

In the demo above, we passed in a list of *device(s)*, like so:

.. code-block:: python

    LiveTable([motor])

Internally, ``LiveTable`` obtains the name(s) of the field(s) produced by
reading ``motor``. You can do this yourself too:

.. ipython:: python

    list(motor.describe().keys())

In the general case, a device can produce tens or even hundreds of separate
readings, and it can be useful to spell out specific fields rather than a whole
device.

.. code-block:: python

    # the field 'motor', in quotes, not the device, motor
    LiveTable(['motor'])

In fact, almost all other callbacks (including :ref:`LivePlot`) *require* a
specific field. They will not accept a device because it may have more than one
field.

.. autoclass:: bluesky.callbacks.LiveTable

.. warning

   The data in the table is formatted according to its type and the
   precision reported by the control system.  If you are seeing too
   many or too few decimal places, this should be adjusted at the
   controls system level.  In EPICS, this is typically the ``.PREC``
   field on the record.

.. _kickers:

Aside: Making plots update live
+++++++++++++++++++++++++++++++

.. note::

    If you are a user working with a pre-configured setup, you can probably
    skip this. Come back if your plots are not appearing / updating.

    This configuration is typically performed in an IPython profile startup
    script so that is happens automatically at startup time.

To make plots live-update while the RunEngine is executing a plan, you have run
this command once. In an IPython terminal, the command is:

.. code-block:: python

    %matplotlib qt
    from bluesky.utils import install_qt_kicker
    install_qt_kicker()

If you are using a Jupyter notebook, the command is:

.. code-block:: python

    %matplotlib notebook
    from bluesky.utils import install_nb_kicker
    install_nb_kicker()

Why? The RunEngine and matplotlib (technically, matplotlib's Qt backend) both
use an event loop. The RunEngine takes control of the event loop while it is
executing a plan. The kicker function periodically "kicks" the Qt event loop so
that the plots can re-draw while the RunEngine is running.

The ``%matplotlib ...`` command is standard setup, having nothing to do with
bluesky in particular. See
`the relevant section of the IPython documentation <https://ipython.readthedocs.io/en/stable/interactive/magics.html?highlight=matplotlib#magic-matplotlib>`_
for details.

.. autofunction:: bluesky.utils.install_kicker
.. autofunction:: bluesky.utils.install_qt_kicker
.. autofunction:: bluesky.utils.install_nb_kicker

.. _liveplot:

LivePlot (for scalar data)
++++++++++++++++++++++++++

Plot scalars. Example:

.. code-block:: python

    from bluesky.plans import scan
    from ophyd.sim import det, motor
    from bluesky.callbacks.mpl_plotting import LivePlot

    RE(scan([det], motor, -5, 5, 30), LivePlot('det', 'motor'))

.. plot::

    from bluesky import RunEngine
    from bluesky.plans import scan
    from ophyd.sim import det, motor
    from bluesky.callbacks.mpl_plotting import LivePlot
    RE = RunEngine({})
    RE(scan([det], motor, -5, 5, 30), LivePlot('det', 'motor'))

To customize style, pass in any
`matplotlib line style keyword argument <http://matplotlib.org/api/lines_api.html#module-matplotlib.lines>`_.
(``LivePlot`` will pass it through to ``Axes.plot``.) Example:

.. code-block:: python

    RE(scan([det], motor, -5, 5, 30),
       LivePlot('det', 'motor', marker='x', markersize=10, color='red'))

.. plot::

    from bluesky import RunEngine
    from bluesky.plans import scan
    from ophyd.sim import det, motor
    from bluesky.callbacks.mpl_plotting import LivePlot
    RE = RunEngine({})
    RE(scan([det], motor, -5, 5, 30),
       LivePlot('det', 'motor', marker='x', markersize=10, color='red'))

.. autoclass:: bluesky.callbacks.mpl_plotting.LivePlot

Live Image
++++++++++

.. autoclass:: bluesky.callbacks.broker.LiveImage

.. _liveraster:

LiveGrid (gridded heat map)
+++++++++++++++++++++++++++

Plot a scalar value as a function of two variables on a regular grid. Example:

.. code-block:: python

    from bluesky.plans import grid_scan
    from ophyd.sim import det4, motor1, motor2
    from bluesky.callbacks.mpl_plotting import LiveGrid

    RE(grid_scan([det4], motor1, -3, 3, 6, motor2, -5, 5, 10, False),
       LiveGrid((6, 10), 'det4'))

.. plot::

    from bluesky import RunEngine
    from bluesky.plans import grid_scan
    from ophyd.sim import det4, motor1, motor2
    from bluesky.callbacks.mpl_plotting import LiveGrid
    motor1.delay = 0
    motor2.delay = 0
    RE = RunEngine({})
    RE(grid_scan([det4], motor1, -3, 3, 6, motor2, -5, 5, 10, False),
       LiveGrid((6, 10), 'det4'))

.. autoclass:: bluesky.callbacks.mpl_plotting.LiveGrid

LiveScatter (scattered heat map)
++++++++++++++++++++++++++++++++

Plot a scalar value as a function of two variables. Unlike
:class:`bluesky.callbacks.mpl_plotting.LiveGrid`, this does not assume a regular grid.
Example:

.. code-block:: python

    from bluesky.plans import grid_scan
    from ophyd.sim import det5, jittery_motor1, jittery_motor2
    from bluesky.callbacks.mpl_plotting import LiveScatter

    # The 'jittery' example motors won't go exactly where they are told to go.

    RE(grid_scan([det5],
                          jittery_motor1, -3, 3, 6,
                          jittery_motor2, -5, 5, 10, False),
       LiveScatter('jittery_motor1', 'jittery_motor2', 'det5',
                xlim=(-3, 3), ylim=(-5, 5)))

.. plot::

    from bluesky import RunEngine
    from bluesky.plans import grid_scan
    from ophyd.sim import det5, jittery_motor1, jittery_motor2
    from bluesky.callbacks.mpl_plotting import LiveScatter
    RE = RunEngine({})
    RE(grid_scan([det5],
                          jittery_motor1, -3, 3, 6,
                          jittery_motor2, -5, 5, 10, False),
       LiveScatter('jittery_motor1', 'jittery_motor2', 'det5',
                xlim=(-3, 3), ylim=(-5, 5)))

.. autoclass:: bluesky.callbacks.mpl_plotting.LiveScatter

LiveFit
+++++++

Perform a nonlinear least squared best fit to the data with a user-defined
model function. The function can depend on any number of independent variables.
We integrate with the package
`lmfit <https://lmfit.github.io/lmfit-py/model.html>`_, which provides a nice
interface for NLS minimization.

In this example, we fit a Gaussian to detector readings as a function of motor
position. First, define a Gaussian function, create an ``lmfit.Model`` from it,
and provide initial guesses for the parameters.

.. code-block:: python

    import numpy as np
    import lmfit

    def gaussian(x, A, sigma, x0):
        return A*np.exp(-(x - x0)**2/(2 * sigma**2))

    model = lmfit.Model(gaussian)
    init_guess = {'A': 2,
                  'sigma': lmfit.Parameter('sigma', 3, min=0),
                  'x0': -0.2}

The guesses can be given as plain numbers or as ``lmfit.Parameter`` objects, as
in the case of 'sigma' above, to specify constraints.

To integrate with the bluesky we need to provide:

* the field with the dependent variable (in this example, ``'noisy_det'``)
* a mapping between the name(s) of independent variable(s) in
  the function (``'x'``) to the corresponding field(s) in the data
  (``'motor'``)
* any initial guesses expected by the model (defined above)

.. code-block:: python

    from bluesky.plans import scan
    from ophyd.sim import motor, noisy_det
    from bluesky.callbacks import LiveFit

    lf = LiveFit(model, 'noisy_det', {'x': 'motor'}, init_guess)

    RE(scan([noisy_det], motor, -1, 1, 100), lf)
    # best-fit values for 'A', 'sigma' and 'x0' are in lf.result.values

The fit results are accessible in the ``result`` attribute of the callback.
For example, the center of the Gaussian is ``lf.result.values['x0']``. This
could be used in a next step, like so:

.. code-block:: python

    x0 = lf.result.values['x0']
    RE(scan([noisy_det], x0 - 1, x0 + 1, 100))

Refer the
`lmfit documentation <https://lmfit.github.io/lmfit-py/model.html#the-modelresult-class>`_
for more about ``result``.

This example uses a model with two independent variables, x and y.

.. code-block:: python

    from ophyd.sim import motor1, motor2, det4

    def gaussian(x, y, A, sigma, x0, y0):
        return A*np.exp(-((x - x0)**2 + (y - y0)**2)/(2 * sigma**2))

    # Specify the names of the independent variables to Model.
    model = lmfit.Model(gaussian, ['x', 'y'])

    init_guess = {'A': 2,
                  'sigma': lmfit.Parameter('sigma', 3, min=0),
                  'x0': -0.2,
                  'y0': 0.3}

    lf = LiveFit(model, 'det4', {'x': 'motor1', 'y': 'motor2'}, init_guess)

    # Scan a 2D mesh.
    RE(grid_scan([det4], motor1, -1, 1, 20, motor2, -1, 1, 20, False),
       lf)

By default, the fit is recomputed every time a new data point is available. See
the API documentation below for other options. Fitting does not commence until
the number of accumulated data points is equal to the number of free parameters
in the model.

.. autoclass:: bluesky.callbacks.LiveFit

LiveFitPlot
+++++++++++

This is a variation on ``LivePlot`` that plots the best fit curve from
``LiveFit``. It applies to 1D model functions only.

Repeating the example from ``LiveFit`` above, adding a plot:

.. code-block:: python

    # same as above...

    import numpy as np
    import lmfit
    from bluesky.plans import scan
    from ophyd.sim import motor, noisy_det
    from bluesky.callbacks import LiveFit

    def gaussian(x, A, sigma, x0):
        return A*np.exp(-(x - x0)**2/(2 * sigma**2))

    model = lmfit.Model(gaussian)
    init_guess = {'A': 2,
                  'sigma': lmfit.Parameter('sigma', 3, min=0),
                  'x0': -0.2}

    lf = LiveFit(model, 'noisy_det', {'x': 'motor'}, init_guess)

    # now add the plot...

    from bluesky.callbacks.mpl_plotting import LiveFitPlot
    lpf = LiveFitPlot(lf, color='r')

    RE(scan([noisy_det], motor, -1, 1, 100), lfp)

    # Notice that we did'nt need to subscribe lf directly, just lfp.
    # But, as before, the results are in lf.result.

.. plot::

    import numpy as np
    import lmfit
    from bluesky.plans import scan
    from ophyd.sim import motor, noisy_det
    from bluesky.callbacks import LiveFit
    from bluesky.callbacks.mpl_plotting import LiveFitPlot
    from bluesky import RunEngine

    RE = RunEngine({})

    def gaussian(x, A, sigma, x0):
        return A*np.exp(-(x - x0)**2/(2 * sigma**2))

    model = lmfit.Model(gaussian)
    init_guess = {'A': 2,
                  'sigma': lmfit.Parameter('sigma', 3, min=0),
                  'x0': -0.2}

    lf = LiveFit(model, 'noisy_det', {'x': 'motor'}, init_guess)
    lfp = LiveFitPlot(lf, color='r')

    RE(scan([noisy_det], motor, -1, 1, 100), lfp)

We can use the standard ``LivePlot`` to show the data on the same axes.
Notice that they can styled independently.

.. code-block:: python

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()  # explitly create figure, axes to use below
    lfp = LiveFitPlot(lf, ax=ax, color='r')
    lp = LivePlot('noisy_det', 'motor', ax=ax, marker='o', linestyle='none')

    RE(scan([noisy_det], motor, -1, 1, 100), [lp, lfp])

.. plot::

    import numpy as np
    import lmfit
    from bluesky.plans import scan
    from ophyd.sim import motor, noisy_det
    from bluesky.callbacks import LiveFit
    from bluesky.callbacks.mpl_plotting import LivePlot, LiveFitPlot
    from bluesky import RunEngine

    RE = RunEngine({})

    def gaussian(x, A, sigma, x0):
        return A*np.exp(-(x - x0)**2/(2 * sigma**2))

    model = lmfit.Model(gaussian)
    init_guess = {'A': 2,
                  'sigma': lmfit.Parameter('sigma', 3, min=0),
                  'x0': -0.2}

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    lf = LiveFit(model, 'noisy_det', {'x': 'motor'}, init_guess)
    lfp = LiveFitPlot(lf, ax=ax, color='r')
    lp = LivePlot('noisy_det', 'motor', ax=ax, marker='o', linestyle='none')

    RE(scan([noisy_det], motor, -1, 1, 50), [lfp, lp])
    plt.draw()

.. autoclass:: bluesky.callbacks.mpl_plotting.LiveFitPlot

PeakStats
++++++++++

Compute statistics of peak-like data. Example:

.. code-block:: python

    from bluesky.callbacks.fitting import PeakStats
    from ophyd.sim import motor, det
    from bluesky.plans import scan

    ps = PeakStats('motor', 'det')
    RE(scan([det], motor, -5, 5, 10), ps)

Now attributes of ``ps``, documented below, contain various peak statistics.
There is also a convenience function for plotting:

.. code-block:: python

    from bluesky.callbacks.mpl_plotting import plot_peak_stats

    plot_peak_stats(ps)

.. plot::

    from bluesky import RunEngine
    from bluesky.callbacks.fitting import PeakStats
    from bluesky.callbacks.mpl_plotting import plot_peak_stats
    from ophyd.sim import motor, det
    from bluesky.plans import scan

    RE = RunEngine({})
    ps = PeakStats('motor', 'det')
    RE(scan([det], motor, -5, 5, 10), ps)
    plot_peak_stats(ps)

.. autoclass:: bluesky.callbacks.fitting.PeakStats
.. autofunction:: bluesky.callbacks.mpl_plotting.plot_peak_stats

.. _best_effort_callback:

Best-Effort Callback
--------------------

.. warning::

    This is a new, experimental feature. It will likely be changed in future
    releases in a way that is not backward-compatible.

This is meant to be permanently subscribed to the RunEngine like so:

.. code-block:: python

    # one-time configuration
    from bluesky.callbacks.best_effort import BestEffortCallback
    bec = BestEffortCallback()
    RE.subscribe(bec)

It provides best-effort plots and visualization for *any* plan. It uses the
'hints' key provided by the plan, if present. (See the source code of the
plans in :mod:`bluesky.plans` for examples.)

.. ipython:: python
    :suppress:

    from bluesky.callbacks.best_effort import BestEffortCallback
    bec = BestEffortCallback()
    RE.subscribe(bec)

.. ipython:: python

    from ophyd.sim import det1, det2
    from bluesky.plans import scan

    dets = [det1, det2]

    RE(scan(dets, motor, 1, 5, 5))  # automatically prints table, shows plot

.. plot::

    from bluesky import RunEngine
    from bluesky.plans import scan
    from ophyd.sim import det, motor
    from bluesky.callbacks.best_effort import BestEffortCallback
    RE = RunEngine({})
    bec = BestEffortCallback()
    RE.subscribe(bec)
    RE(scan([det], motor, 1, 5, 5))

Use these methods to toggle on or off parts of the functionality.

.. currentmodule:: bluesky.callbacks.best_effort

.. autosummary::
    :toctree: generated

    BestEffortCallback
    BestEffortCallback.enable_heading
    BestEffortCallback.disable_heading
    BestEffortCallback.enable_table
    BestEffortCallback.disable_table
    BestEffortCallback.enable_baseline
    BestEffortCallback.disable_baseline
    BestEffortCallback.enable_plots
    BestEffortCallback.disable_plots

Blacklist plotting certain streams using the ``bec.noplot_streams`` attribute,
which is a list of stream names.  The blacklist is set to ``['baseline']`` by
default.

The attribute ``bec.overplot`` can be used to control whether line plots for
subsequent runs are plotted on the same axes. It is ``True`` by default.
Overplotting only occurs if the names of the axes are the same from one plot
to the next.

Peak Stats
++++++++++

For each plot, simple peak-fitting is performed in the background. Of
course, it may or may not be applicable depending on your data, and it is
not shown by default. To view fitting annotations in a plot, click the
plot area and press Shift+P. (Lowercase p is a shortcut for
"panning" the plot.)

To access the peak-fit statistics programmatically, use ``bec.peaks``.

.. _hints:

Hints
+++++

The best-effort callback aims to print and plot useful information without
being overwhelmingly comprehensive. Its usefulness is improved and tuned by the
``hints`` attribute on devices (if available) and ``hints`` metadata injected
by plans (if available). If either or both of these are not available, the
best-effort callback still makes a best effort to display something useful.

The contents of hints *do not at all affect what data is saved*. The content
only affect what is displayed automatically by the best-effort callback and
other tools that opt to look at the hints. Additional callbacks may still be
set up for live or *post-facto* visualization or processing that do more
specific things without relying on hints.

The ``hints`` attribute or property on devices is a dictionary with the key
``'fields'`` mapped to a list of fields.

On movable devices such as motors or temperature controllers, these fields are
expected to comprise the independent axes of the device. A motor that reads
the fields ``['x', 'x_setpoint']`` might provide the hint ``{'fields': ['x']}``
to indicate that it has one independent axis and that the field ``x`` is the best
representation of its value.

A readable device might report many fields like
``['chan1', 'chan2', 'chan3', 'chan4', 'chan5']`` but perhaps only a couple are
usually interesting. A useful hint might narrow them down to
``{'fields': ['chan1', 'chan2']}`` so that a "best-effort" plot does not
display an overwhelming amount of information.

The hints provided by the devices are read by the RunEngine and collated in the
:doc:`Event Descriptor documents <event_descriptors>`.

The plans generally know which devices are being used as dependent and
independent variables (i.e., which are being "scanned" over), and they may
provide this information via a ``'hints'`` metadata key that they inject into
the start document along with the rest of their metadata. Examples:

.. code-block:: python

    # The pattern is
    # {'dimensions': [(fields, stream_name), (fields, stream_name), ...]}

    # a scan over time
    {'dimensions': [(('time',), 'primary')]}

    # a one-dimensional scan
    {'dimensions': [(motor.hints['fields'], 'primary')]}

    # a two-dimensional scan
    {'dimensions': [(x_motor.hints['fields'], 'primary'),
                    (y_motor.hints['fields'], 'primary')]}

    # an N-dimensional scan
    {'dimensions': [(motor.hints['fields'], 'primary') for motor in motors]}

It's possible to adjust hints interactively, but they are generally intended to
be set in a startup file. Err on the side of displaying more information than
you need to see, and you will rarely need to adjust them.

Plans may also hint that their data is sampled on a regular rectangular grid
via the hint ``{'gridding': 'rectilinear'}``. This is useful, for example, for
decided whether to visualize 2D data with LiveGrid or with LiveScatter.

.. _export:

Callback for Export
-------------------

Exporting Image Data as TIFF Files
++++++++++++++++++++++++++++++++++

First, compose a filename template. The template can include metadata or event
data from the scan.

.. code-block:: python

    # a template that includes the scan ID and sequence number in each filename
    template = "output_dir/{start[scan_id]}_{event[seq_num]}.tiff"

    # a template that sorts files into directories based user and scan ID
    template = "output_dir/{start[user]}/{start[scan_id]}/{event[seq_num]}.tiff"

    # a more complex template includes actual measurements in the filenames
    template = ("output_dir/{start[scan_id]}_{start[sample_name]}_"
                "{event[data][temperature]}_{event[seq_num]}.tiff")

Above, we are using a Python language feature called format strings. Notice
that inside the curly brackets we don't use quotes around the key names; it's
``{event[seq_num]}`` not ``{event['seq_num']}``.

If each image data point is actually a stack of 2D image planes, the template
must also include ``{i}``, which will count through the image planes in the
stack.

.. note::

    Most metadata comes from the "start" document, hence ``start.scan_id``
    above.  Review the :doc:`documents` section for details.

Create a callback that exports TIFFs using your template.

.. code-block:: python

    from bluesky.callbacks.broker import LiveTiffExporter

    exporter = LiveTiffExporter('image', template)

Finally, to export all the images from a run when it finishes running, wrap the
exporter in ``post_run`` and subscribe.

.. code-block:: python

    from bluesky.callbacks.broker import post_run

    RE.subscribe(post_run(exporter))

It also possible to write TIFFs live, hence the name ``LiveTiffExporter``, but
there is an important disadvantage to doing this subscription in the same
process: progress of the experiment may be intermittently slowed while data is
written to disk. In some circumstances, this affect on the timing of the
experiment may not be acceptable.

.. code-block:: python

    RE.subscribe(exporter)

There are more configuration options available, as given in detail below. It is
recommended to use these expensive callbacks in a separate process.

.. autoclass:: bluesky.callbacks.broker.LiveTiffExporter

Export All Data and Metadata in an HDF5 File
++++++++++++++++++++++++++++++++++++++++++++

A Stop Document is emitted at the end of every run. Subscribe to it, using it
as a cue to load the dataset via the DataBroker and export an HDF5 file
using `suitcase <https://nsls-ii.github.io/suitcase>`_.


Working example:

.. code-block:: python

    from databroker import DataBroker as db
    import suitcase

    def suitcase_as_callback(name, doc):
        if name != 'stop':
            return
        run_start_uid = doc['run_start']
        header = db[run_start_uid]
        filename = '{}.h5'.format(run_start_uid)
        suitcase.export(header, filename)

    RE.subscribe(suitcase_as_callback, 'stop')

Export Metadata to the Olog
+++++++++++++++++++++++++++

The `Olog <http://olog.github.io/2.2.7-SNAPSHOT/>`_ ("operational log") is an
electronic logbook. We can use a callback to automatically generate log entries
at the beginning of a run. The Python interface to Olog is not straightforward,
so there is some boilerplate:

.. code-block:: python

    from functools import partial
    from pyOlog import SimpleOlogClient
    from bluesky.callbacks.olog import logbook_cb_factory

    # Set up the logbook. This configures bluesky's summaries of
    # data acquisition (scan type, ID, etc.).

    LOGBOOKS = ['Data Acquisition']  # list of logbook names to publish to
    simple_olog_client = SimpleOlogClient()
    generic_logbook_func = simple_olog_client.log
    configured_logbook_func = partial(generic_logbook_func, logbooks=LOGBOOKS)

    cb = logbook_cb_factory(configured_logbook_func)
    RE.subscribe(cb, 'start')

The module ``bluesky.callbacks.olog`` includes some templates that format the
data from the 'start' document into a readable log entry. You can also write
customize templates and pass them to ``logbook_cb_factory``.

You may specify a custom template. Here is a very simple example; see the
`source code <https://github.com/NSLS-II/bluesky/blob/master/bluesky/callbacks/olog.py>`_
for a more complex example (the default template).

.. code-block:: python

    CUSTOM_TEMPLATE = """
    My Log Entry

    {{ start.plan_name }}
    Detectors: {{ start.detectors }}
    """

    # Do same boilerplate above to set up configured_logbook_func. Then:
    cb = logbook_cb_factory(configured_logbook_func,
                            desc_template=CUSTOM_TEMPLATE)

You may also specify a variety of different templates that are suitable for
different kinds of plans. The callback will use the ``'plan_name'`` field to
determine which template to use.

.. code-block:: python

    # a template for a 'count' plan (which has no motors)
    COUNT_TEMPLATE = """
    Plan Name: {{ start.plan_name }}
    Detectors: {{ start.detectors }}
    """

    # a template for any plan with motors
    SCAN_TEMPLATE = """
    Plan Name: {{ start.plan_name }}
    Detectors: {{ start.detectors }}
    Motor(s): {{ start.motors }}
    """

    templates = {'count': COUNT_TEMPLATE,
                 'scan': SCAN_TEMPLATE,
                 'rel_scan': SCAN_TEMPLATE}

    # Do same boilerplate above to set up configured_logbook_func. Then:
    cb = logbook_cb_factory(configured_logbook_func,
                            desc_dispatch=templates)

.. autofunction:: bluesky.callbacks.olog.logbook_cb_factory

Verify Data Has Been Saved
--------------------------

The following verifies that all Documents and external files from a run have
been saved to disk and are accessible from the DataBroker.  It prints a message
indicating success or failure.

Note: If the data collection machine is not able to access the machine where
some external data is being saved, it will indicate failure. This can be a
false alarm.

.. code-block:: python

    from bluesky.callbacks.broker import post_run, verify_files_saved

    RE.subscribe(post_run(verify_files_saved))

.. _debugging_callbacks:

Ignoring Callback Exceptions
----------------------------

If an exception is raised while processing a callback, the error can interrupt
data collection. Sometimes, this is good: if, for example, the callback that is
saving your data encounters an error, you want to know immediately rather than
continuing to *think* you are collecting data when in fact it is being lost.
But in many situations, such as visualization or first-pass data processing, it
is usually better for data collection to proceed even if a callback fails.
These decorators may be used to wrap callbacks so that any errors they
encounter are converted to log messages.

.. autofunction:: bluesky.callbacks.core.make_callback_safe

.. autofunction:: bluesky.callbacks.core.make_class_safe

It is also possible to configure the RunEngine to ignore *all* callback
exceptions globally, but this feature is not recommended.

.. code-block:: python

    RE.ignore_callback_exceptions = False

.. versionchanged:: 0.6.4

   In bluesky version 0.6.4 (September 2016) the default value was changed from
   ``True`` to ``False``.

.. _filtering:

Filtering by Document Type
--------------------------

There are four "subscriptions" that a callback to receive documents from:

* 'start'
* 'stop'
* 'event'
* 'descriptor'

Additionally, there is an 'all' subscription.

The command:

.. code-block:: python

    RE(plan(), cb)

is a shorthand that is normalized to ``{'all': [cb]}``. To receive only certain
documents, specify the document routing explicitly. Examples:

.. code-block:: python

    RE(plan(), {'start': [cb]}
    RE(plan(), {'all': [cb1, cb2], 'start': [cb3]})

The ``subs_decorator``, presented above, accepts the same variety of inputs.

Writing Custom Callbacks
------------------------

Any function that accepts a Python dictionary as its argument can be used as
a callback. Refer to simple examples above to get started.

Two Simple Custom Callbacks
+++++++++++++++++++++++++++

These simple examples illustrate the concept and the usage.

First, we define a function that takes two arguments

#. the name of the Document type ('start', 'stop', 'event', or 'descriptor')
#. the Document itself, a dictionary

This is the *callback*.

.. ipython:: python

    def print_data(name, doc):
        print("Measured: %s" % doc['data'])

Then, we tell the RunEngine to call this function on each Event Document.
We are setting up a *subscription*.

.. ipython:: python

    from ophyd.sim import det
    from bluesky.plans import count

    RE(count([det]), {'event': print_data})

Each time the RunEngine generates a new Event Document (i.e., data point)
``print_data`` is called.

There are five kinds of subscriptions matching the four kinds of Documents plus
an 'all' subscription that receives all Documents.

* 'start'
* 'descriptor'
* 'event'
* 'stop'
* 'all'

We can use the 'stop' subscription to trigger automatic end-of-run activities.
For example:

.. code-block:: python

    def celebrate(name, doc):
        # Do nothing with the input; just use it as a signal that run is over.
        print("The run is finished!")

Let's use both ``print_data`` and ``celebrate`` at once.

.. code-block:: python

    RE(plan(), {'event': print_data, 'stop': celebrate})

Using multiple document types
+++++++++++++++++++++++++++++

Some tasks use only one Document type, but we often need to use more than one.
For example, LiveTable uses 'start' kick off the creation of a fresh table,
it uses 'event' to see the data, and it uses 'stop' to draw the bottom border.

A convenient pattern for this kind of subscription is a class with a method
for each Document type.

.. code-block:: python

    from bluesky.callbacks import CallbackBase

    class MyCallback(CallbackBase):
        def start(self, doc):
            print("I got a new 'start' Document")
            # Do something
        def descriptor(self, doc):
            print("I got a new 'descriptor' Document")
            # Do something
        def event(self, doc):
            print("I got a new 'event' Document")
            # Do something
        def stop(self, doc):
            print("I got a new 'stop' Document")
            # Do something

The base class, ``CallbackBase``, takes care of dispatching each Document to
the corresponding method. If your application does not need all four, you may
simple omit methods that aren't required.

.. _zmq_callback:

Subscriptions in Separate Processes or Host with 0MQ
----------------------------------------------------

Because subscriptions are processed during a scan, it's possible that they can
slow down data collection. We mitigate this by making the subscriptions run in
a separate process.

In the main process, where the RunEngine is executing the plan, a ``Publisher``
is created. It subscribes to the RunEngine. It serializes the documents it
receives and it sends them over a socket to a 0MQ proxy which rebroadcasts the
documents to any number of other processes or machines on the network.

These other processes or machines set up a ``RemoteDispatcher`` which connects
to the proxy receives the documents, and then runs callbacks just as they would
be run if they were in the local ``RunEngine`` process.

Multiple Publishers (each with its own RunEngine) can send documents to the
same proxy. RemoteDispatchers can filter the document stream based a byte
prefix.

Minimal Example
+++++++++++++++

Start a 0MQ proxy using the CLI packaged with bluesky. It requires two ports as
arguments.

.. code-block:: bash

    bluesky-0MQ-proxy 5577 5578

Alternatively, you can start the proxy using a Python API:

.. code-block:: python

    from bluesky.callbacks.zmq import Proxy
    proxy = Proxy(5577, 5578)
    proxy.start()

Start a callback that will receive documents from the proxy and, in this
simple example, just print them.

.. code-block:: python

    from bluesky.callbacks.zmq import RemoteDispatcher
    d = RemoteDispatcher('localhost:5578')
    d.subscribe(print)

    # when done subscribing things and ready to use:
    d.start()  # runs event loop forever

As `described above <kickers>`_, if you want to use any live-updating plots,
you will need to install a "kicker". It needs to be installed on the same
event loop used by the RemoteDispatcher, like so, and it must be done before
calling ``d.start()``.

.. code-block:: python

    from bluesky.utils import install_qt_kicker
    install_qt_kicker(loop=d.loop)

In a Jupyter notebook, replace ``install_qt_kicker`` with
``install_nb_kicker``.

On the machine/process where you want to collect data, hook up a subscription
to publish documents to the proxy.

.. code-block:: python

    # Create a RunEngine instance (or, of course, use your existing one).
    from bluesky import RunEngine, Msg
    RE = RunEngine({})

    from bluesky.callbacks.zmq import Publisher
    publisher = Publisher('localhost:5577')
    RE.subscribe(publisher)

Finally, execute a plan with the RunEngine. As a result, the callback in the
RemoteDispatcher should print the documents generated by this plan.

Publisher / RemoteDispatcher API
++++++++++++++++++++++++++++++++

.. autoclass:: bluesky.callbacks.zmq.Proxy
.. autoclass:: bluesky.callbacks.zmq.Publisher
.. autoclass:: bluesky.callbacks.zmq.RemoteDispatcher


Secondary Event Stream
----------------------
For certain applications, it may desirable to interpret event documents as
they are created instead of waiting for them to reach offline storage. In order
to keep this information completely quarantined from the raw data, the
:class:`.LiveDispatcher` presents a completely unique stream that can be
subscribed to using the same syntax as the RunEngine.

In the majority of applications of :class:`.LiveDispatcher`, it is expected
that subclasses are created to implement online analysis. This secondary event
stream can be displayed and saved offline using the same callbacks that you
would use to display the raw data.

Below is an example using the `streamz
<https://streamz.readthedocs.io/en/latest>`_ library to average a number of
events together. The callback can be configured by looking at the start
document metadata, or at initialization time. Events are then received and
stored by the ``streamz`` network and a new averaged event is emitted when the
correct number of events are in the cache. The important thing to note here is
that the analysis only handles creating new ``data`` keys, but the descriptors,
sequence numbering and event ids are all handled by the base `LiveDispatcher`
class.

.. code-block:: python

    class AverageStream(LiveDispatcher):
        """Stream that averages data points together"""
        def __init__(self, n=None):
            self.n = n
            self.in_node = None
            self.out_node = None
            self.averager = None
            super().__init__()

        def start(self, doc):
            """
            Create the stream after seeing the start document

            The callback looks for the 'average' key in the start document to
            configure itself.
            """
            # Grab the average key
            self.n = doc.get('average', self.n)
            # Define our nodes
            if not self.in_node:
                self.in_node = streamz.Source(stream_name='Input')

            self.averager = self.in_node.partition(self.n)

            def average_events(cache):
                average_evt = dict()
                desc_id = cache[0]['descriptor']
                # Check that all of our events came from the same configuration
                if not all([desc_id == evt['descriptor'] for evt in cache]):
                    raise Exception('The events in this bundle are from '
                                    'different configurations!')
                # Use the last descriptor to avoid strings and objects
                data_keys = self.raw_descriptors[desc_id]['data_keys']
                for key, info in data_keys.items():
                    # Information from non-number fields is dropped
                    if info['dtype'] in ('number', 'array'):
                        # Average together
                        average_evt[key] = np.mean([evt['data'][key]
                                                    for evt in cache], axis=0)
                return {'data': average_evt, 'descriptor': desc_id}

            self.out_node = self.averager.map(average_events)
            self.out_node.sink(self.process_event)
            super().start(doc)

        def event(self, doc):
            """Send an Event through the stream"""
            self.in_node.emit(doc)

        def stop(self, doc):
            """Delete the stream when run stops"""
            self.in_node = None
            self.out_node = None
            self.averager = None
            super().stop(doc)


LiveDispatcher API
++++++++++++++++++
.. autoclass:: bluesky.callbacks.stream.LiveDispatcher
   :members:
How Bluesky Interfaces with Hardware
====================================

.. _hardware_interface:

Overview
--------

Bluesky interacts with hardware through a high-level abstraction, leaving the
low-level details of communication as a separate concern. In bluesky's view,
*all* devices are in a sense "detectors," in that they can be read. A subset
of these devices are "positioners" that can also be set (i.e., written to or
moved).

In short, each device is represented by a Python object that has attributes and
methods with certain established names. We have taken pains to make this
interface as slim as possible, while still being general enough to address
every kind of hardware we have encountered.

Specification
-------------

.. _status_obj_api:

Status object
+++++++++++++

The interface of a "status" object, which the ``RunEngine`` uses to
asynchronously monitor the compeletion of having triggered or set a device.

.. autoclass:: bluesky.protocols.Status
   :members:
   :undoc-members:

If ``success`` is ``False`` when the Status is marked done, this is taken
to mean, "We have given up." For example, "The motor is stuck and will
never get where it is going." A ``FailedStatus`` exception will be raised
inside the RunEngine.

Additionally, ``Status`` objects may (optionally) add a watch function that
conforms to the following definition

   .. method:: watch(func)

        Subscribe to notifications about progress. Useful for progress bars.

        **Parameters**

        func : callable
            Expected to accept the keyword arguments:

                * ``name``
                * ``current``
                * ``initial``
                * ``target``
                * ``unit``
                * ``precision``
                * ``fraction``
                * ``time_elapsed``
                * ``time_remaining``

            Any given call to ``func`` may only include a subset of these
            parameters, depending on what the status object knows about its own
            progress.

            The ``fraction`` argument accepts a single float representing fraction
            remaining.
            A fraction of zero indicates completion.
            A fraction of one indicates progress has not started.

Readable Device
+++++++++++++++

The interface of a readable device:

.. autoclass:: bluesky.protocols.Readable
    :members:
    :undoc-members:

Movable (or "Settable")  Device
+++++++++++++++++++++++++++++++

The interface of a movable device extends the interface of a readable device
with the following additional methods and attributes.


.. autoclass:: bluesky.protocols.Movable
    :members:
    :show-inheritance:

    .. attribute:: position

        A heuristic that describes the current position of a device as a
        single scalar, as opposed to the potentially multi-valued description
        provided by ``read()``.

        Optional: bluesky itself does not use the position attribute, but other
        parts of the ecosystem might.
        Developers are encouraged to implement this attribute where possible.


"Flyer" Interface
+++++++++++++++++

*For context on what we mean by "flyer", refer to the section on :doc:`async`.*

The interface of a "flyable" device is separate from the interface of a readable
or settable device, though there is some overlap.


.. autoclass:: bluesky.protocols.Flyable
    :members:
    :undoc-members:


Optional Interfaces
-------------------

These are additional interfaces for providing *optional* behavior to ``Readable``, ``Movable``,
and ``Flyable`` devices.

The methods described here are either hooks for various plans/RunEngine messages which are
ignored if not present or required by only a subset of RunEngine messages.
In the latter case, the RunEngine may error if it tries to use a device which does not define
the required method.

.. autoclass:: bluesky.protocols.Stageable
    :members:

.. autoclass:: bluesky.protocols.Subscribable
    :members:

.. autoclass:: bluesky.protocols.Pausable
    :members:

.. autoclass:: bluesky.protocols.Stoppable
    :members:

.. autoclass:: bluesky.protocols.Checkable
    :members:

.. autoclass:: bluesky.protocols.Hinted
    :members:


Implementations
---------------

Real Hardware
+++++++++++++

The `ophyd
<https://nsls-ii.github.io/ophyd>`_ package implements this interface for
a wide variety of hardware, communicating using
`EPICS <http://www.aps.anl.gov/epics/>`_ via the Python bindings
`pyepics <http://cars9.uchicago.edu/software/python/pyepics3/>`_.Other control
systems (Tango, LabView, etc.) could be integrated with bluesky in the future
by implementing this same interface.

Simulated Hardware
++++++++++++++++++

A toy "test" implementation the interface is included in the
:mod:`ophyd.sim` module. These implementations act as simulated hardware,
and we use them extensively in examples, demos, and the test suite. They can
also be useful for exercising analysis workflows before running a real
experiment. API documentation is below.
.. _msg:

Message Protocol
================

*Note: This is a technical document not optimized for user readability.*

Overview
--------

A *plan* is a sequence of atomic operations describing a data acquisition
procedure. Each operation is represented by a ``bluesky.Msg`` ("message")
object. A plan may be implemented as a simple list of messages:

.. code-block:: python

    from bluesky import Msg

    # (Behold, the most boring data acquisition ever conducted!)
    plan = [Msg('open_run'), Msg('close_run')]

or as a generator the yields messages one at time:

.. code-block:: python

    def plan():
        yield Msg('open_run')
        yield Msg('close_run')

The above examples are equivalent. For more sophisticated uses, the second one
is more powerful, as it can incorporate loops, conditionals, adaptive logic ---
generally any Python code.

But, crucially, the plan code itself must not communicate with hardware.
(You should never put ``epics.caput(...)`` in a plan!) Rather, each operation
is represented by a ``Msg`` object that *describes* what should be done. This
makes it safe to introspect the plan for error-checking, simulation, and
visualization purposes --- without touching real hardware. For example, we
could print each message in the plan like so:

.. code-block:: python

    plan = [Msg('open_run'), Msg('close_run')]

    # a very, very simple 'plan simulator'
    for msg in plan:
        print(msg)

A ``Msg`` has five members, accessible as attributes:

- command
- obj
- args
- kwargs
- run

where ``command`` must be one of a controlled list of commands, ``obj`` is the
object (i.e. Device) to apply the command to, if applicable, ``args`` and
``kwargs`` are arguments to the command and ``run`` is a user-defined run key.
The run key is used by Run Engine to associate each message with one of the open runs,
manage the state of each open run, and route run data to a separate set of callbacks
(see documentation on Multi-Run Plans).

To execute the plan, the :doc:`RunEngine <run_engine>` consumes it, one message at a time.

.. code-block:: python

    def very_simple_run_engine(plan):
        for msg in plan:
            # Process the msg.

The ``RunEngine`` has a registry which is used to dispatch the ``Msg`` objects
based on the value of the ``Msg.command``. For example, if the RunEngine
receives the message ``Msg('set', motor, 5)``, the RunEngine will:

1. Identify that the command for this message is ``'set'``.
2. Look up ``'set'`` in its command registry and find that it is mapped to
   ``RunEngine._set``.
3. Pass ``Msg('set', motor, 5)`` to its ``_set`` method.
4. Inside ``_set``, call ``motor.set(5)``. (This is where the actual
   communication with hardware occurs.)
5. Update some internal caches that will be useful later. For example, it will
   keep track of that fact that ``motor`` may be in motion so that it can stop
   it safely if an error occurs. This illustrates another important reason that
   plans must always yield messages to interact with hardware and absolutely
   never communicate with hardware directly. Calling ``epics.caput`` inside a
   plan prevents the RunEngine from knowing about it and thus circumvents
   its facilities for putting devices in a safe state in the event of an
   unexpected exit or error.

A standard set of commands are registered by default.  By convention, a ``Msg``
with the command ``'name'`` is mapped to a coroutine method on the RunEngine
named ``_name``, as in ``'set'`` -> ``RunEngine._set`` in the example above.
Users can register their own coroutines to add custom commands, though this is
very rarely necessary.

Some commands do not involve communication with hardware. For example,
``Msg('sleep', None, 5)`` causes the RunEngine to sleep for 5 seconds. ``None``
is a placeholder for the "object" (Device) which is not applicable for a
``'sleep'`` command. Just as plans should never communicate with hardware
directly, they should also never employ long blocking calls like
``time.sleep()``. Instead, the ``'sleep'`` command, mapped to
``RunEngine._sleep``, integrates with the RunEngine's event loop to sleep in a
non-blocking way that allows for the RunEngine to stay responsive in the
meantime --- watching for user interruptions and possibility collecting data
asynchronously in the background.

Other commands are used to control metadata and I/O. For example,
``Msg('open_run')`` and ``Msg('close_run')`` delineate the scope of one run.
Any keyword arguments passed to the ``'open_run'`` message are interpreted as
metadata, encoded into the RunStart document.

The following is a comprehensive overview of the built-in commands.

.. _commands:

Commands
--------

.. warning::

    This section of the documentation is incomplete.

These are the 'built in' commands, some of which are deeply tied to the
state of the `RunEngine` instance.

create
++++++

This command tells the run engine that it should start to collect the results
of ``read`` to create an event.  If this is called twice without a ``save`` or
``drop`` between them it is an exception (as you can not have more than one
open event going at a time).

This relies very heavily on the internal state of the run engine and should not
be overridden by the user.

This call returns `None` back to the co-routine.

This ignores all parts of the `Msg` except the command.

save
++++

This is the pair to ``create`` which bundles and causes ``Event`` documents to
be emitted.  This must be called after a ``create`` or a the scan will die and
raise `IllegalMessageSequence`.

This relies very heavily on the internal state of the run engine and should not
be messed with.

This call returns `None` back to the co-routine.

This ignores all parts of the `Msg` except the command.

read
++++

This causes `read` to be called on the ``obj`` in the message ::

  msg.obj.read(*msg.args, **msg.kwargs)

Anything that is read between a ``create`` and ``save`` will be bundled into
a single event.

This relies very heavily on the internal state of the run engine and should not
be messed with.

Returns the dictionary returned by `read` to the co-routine.

The ``args`` and ``kwargs`` parts of the message are passed to the `read`
method.


null
++++

This is a null message and is ignored by the run engine.  This exists to make
the algebra work.

Returns `None` to the co-routine.

Ignores all values in the `Msg` except the command.

set
+++

Tells a ``Mover`` object to move.  Currently this mimics the epics-like logic
of immediate motion.

stage and unstage
+++++++++++++++++
Instruct the RunEngine to stage/unstage the object. This calls
``obj.stage()``/``obj.unstage``.

Expected message objects are::

    Msg('stage', object)
    Msg('unstage', object)

which results in these calls::

    staged_devices = object.stage()
    unstaged_devices = object.unstage()

where ``staged_devices``/``unstaged_devices`` are a list of the
``ophyd.Device`` (s) that were (un)staged, not status objects.

One may wonder why the return is a list of Devices as opposed to Status
objects, such as in ``set`` and similar ``Msg`` s.
This was debated for awhile. Operations performed during staging are supposed
to involve twiddling configuration, and should happen fast. Staging should not
involve lengthy set calls.

Why a list of the objects staged? Staging a Device causes that Device's
component Devices (if any) to also be staged. All of these children are added
to a list, along with [self], and returned by Device.stage(), so that the plan
can keep track of what has been staged, like so::

    devices_staged = yield Msg('stage', device)

Why would the plan want to know that? It needs to avoid accidentally trying to
stage something twice, such as a staging a parent and then trying to also stage
its child. It's important to avoid that because staging something redundantly
raises an error.


trigger
+++++++

This will call the ``obj.trigger`` method and cache the returned status object
and caches the returned status object.


sleep
+++++

Sleep the event loop.

wait
++++

Block progress until every object that was triggered or set the keyword
argument `group=<GROUP>` is done.

Expected message object is:

Msg('wait', group=<GROUP>)

where ``<GROUP>`` is any hashable key.

wait_for
++++++++
Instruct the ``RunEngine`` to wait for this ``asyncio.Future`` object to be
done. This allows for external arbitrary control of the ``RunEngine``.
Ex ::

    from asyncio.futures import Future
    future = Future()
    future.done() # will give false
    RE(Msg('wait_for', [lambda : future ,]))
    # this sets the future to done
    future.set_result(3)
    future.done() # will give True


input
+++++
Process an input. Allows for user input during a run.

Examples::

    Msg('input', None)
    Msg('input', None, prompt='>')  # customize prompt


checkpoint
++++++++++

Instruct the RunEngine to create a checkpoint so that we can rewind to this
point if necessary.

clear_checkpoint
++++++++++++++++
Clear a set checkpoint.

rewindable
++++++++++

pause
+++++

Request the run engine to pause

Expected message object is::

    Msg('pause', defer=False, name=None, callback=None)


kickoff
+++++++

Start a flyscan object.

collect
+++++++

Collect data cached by a flyer and emit descriptor and event documents.
This calls the ``obj.collect()`` method.

complete
++++++++

Tell a flyer, 'stop collecting, whenever you are ready'.

This calls the method ``obj.complete()`` of the given object. The flyer returns
a status object. Some flyers respond to this command by stopping collection and
returning a finished status object immediately. Other flyers finish their given
course and finish whenever they finish, irrespective of when this command is
issued.


configure
+++++++++

Configure an object.

Expected message object is::

    Msg('configure', object, *args, **kwargs)

which results in this call::

    object.configure(*args, **kwargs)


subscribe
+++++++++
Add a subscription after the run has started.

This, like subscriptions passed to __call__, will be removed at the
end by the RunEngine.

Expected message object is:

    Msg('subscribe', None, callback_function, document_name)

where `document_name` is one of:

    {'start', 'descriptor', 'event', 'stop', 'all'}

and `callback_function` is expected to have a signature of:

    ``f(name, document)``

    where name is one of the ``document_name`` options and ``document``
    is one of the document dictionaries in the event model.

See the docstring of bluesky.run_engine.Dispatcher.subscribe() for more
information.

unsubscribe
+++++++++++

Remove a subscription during a call -- useful for a multi-run call
where subscriptions are wanted for some runs but not others.

Expected message object is::

    Msg('unsubscribe', None, TOKEN)
    Msg('unsubscribe', token=TOKEN)

where ``TOKEN`` is the return value from ``RunEngine._subscribe()``

open_run
++++++++
Instruct the RunEngine to start a new "run"

Expected message object is::

    Msg('open_run', None, **kwargs)

where ``**kwargs`` are any additional metadata that should go into the RunStart
document

close_run
+++++++++

Instruct the RunEngine to write the RunStop document

Expected message object is::

    Msg('close_run', None, exit_status=None, reason=None)

if *exit_stats* and *reason* are not provided, use the values
stashed on the RE.


drop
++++

Drop a bundle of readings without emitting a completed Event document.

This is a command that abandons previous ``create`` and ``read`` commands
without emitting an event. This can be used to drop known bad events
(e.g. no beam) and keep the event document stream clean. It is safe to start
another ``create``, ``read``, ``save`` sequence after a ``drop``.

This must be called after a ``create`` or a the scan will die and raise
`IllegalMessageSequence`.

This call returns `None` back to the co-routine.

This ignores all parts of the `Msg` except the command.


monitor
+++++++
Monitor a signal. Emit event documents asynchronously.

A descriptor document is emitted immediately. Then, a closure is
defined that emits Event documents associated with that descriptor
from a separate thread. This process is not related to the main
bundling process (create/read/save).

Expected message object is::

    Msg('monitor', obj, **kwargs)
    Msg('monitor', obj, name='event-stream-name', **kwargs)

where kwargs are passed through to ``obj.subscribe()``


unmonitor
+++++++++

Stop monitoring; i.e., remove the callback emitting event documents.

Expected message object is::

    Msg('unmonitor', obj)


stop
++++

Stop a device.

Expected message object is::

    Msg('stop', obj)

This amounts to calling ``obj.stop()``.


Registering Custom Commands
---------------------------

The RunEngine can be taught any new commands. They can be registered using the
following methods.

.. automethod:: bluesky.run_engine.RunEngine.register_command
    :noindex:

.. automethod:: bluesky.run_engine.RunEngine.unregister_command
    :noindex:

.. autoattribute:: bluesky.run_engine.RunEngine.commands
    :noindex:

.. automethod:: bluesky.run_engine.RunEngine.print_command_registry
    :noindex:
Progress Bar
************

Bluesky provides a progress bar add-on. For example, two motors moving
simulateously make a display like this:

.. code-block:: none

    mtr1  9%|███▊                                       | 0.09/1.0 [00:00<00:01,  1.21s/deg]
    mtr2100%|████████████████████████████████████████████| 1.0/1.0 [00:01<00:00,  1.12s/deg]

This display includes:

* the name of the device (motor, temperature controller, etc.)
* the distance (or degrees, etc.) traveled so far
* the total distance to be covered
* the time elapsed
* the estimated time remaining
* the rate (determined empirically)

The progress bar relies on the device to report its progress. If a device does
not provide comprehensive information, a simpler progress bar will be shown,
listing the names of devices being waited on and reporting which have
completed.

.. code-block:: none

    mtr1 [No progress bar available.]
    mtr2 [Complete.]

Any time the RunEngine waits on hardware the progress bar is notified. This
includes, for example, waiting for a motor to move or waiting for a detector to
trigger. (In bluesky jargon, the progress bar is notified any time the
RunEngine processes a 'wait' command).

The progress bar is not set up by default. It must be attached to a RunEngine.
This need only be done once (say, in a startup file).

.. code-block:: python

    from bluesky.utils import ProgressBarManager
    
    RE.waiting_hook = ProgressBarManager()

Some motions are very quick and not worth displaying a progress bar for. By
default, a progress bar is only drawn after 0.2 seconds. If an action completes
before then, the progress bar is never shown. To choose a shorter or longer
delay---say 5 seconds---use the parameter ``ProgressBarManager(delay_draw=5)``.

For more technical detail about communication between the device, the
RunEngine, and the ProgressBarManager, read about the ``watch`` method in the
:ref:`status_obj_api` and ``waiting_hook`` in the :doc:`run_engine_api`.

The implementation of the progress bar itself makes use of
`tqdm <https://github.com/tqdm/tqdm/>`_, a lovely Python package for making
a progress bar out of any iterable.
=================
 Release History
=================


v1.8.2 (2021-12-20)
===================

Fixed
-----

* Changed from using ``SafeConfigParser`` to ``ConfigParser`` in
  ``versioneer.py`` (fix to support Python 3.11).

Enhancements
------------

* Added public ``call_returns_result`` property.
* Implemented human-readable printable representation for ``PeakStats``.

Documentation
-------------

* Updated ``RunEngine`` docstring with ``call_returns_result`` property.
* Fixed a small mistake in ``bp.scan`` docstring.
* Added documentation about intended behavior of fraction in the ``watch``
  method of the status object.


v1.8.1 (2021-10-11)
===================

Fixed
-----

* More fixes for Python 3.10 to propagate the ``loop`` kwarg correctly.

Enhancements
------------

* Added optional calculation of the derivative and its statistics (``min``,
  ``max``, ``fwhm``, etc.) to ``PeakStats`` and ``BestEffortCallback``.

Added
-----

* Read-only property ``RunEngine.deferred_pause_requested`` which may be useful
  for `bluesky-queueserver <https://github.com/bluesky/bluesky-queueserver>`_.

Documentation
-------------

* Unpin ``sphinx_rtd_theme``.


v1.8.0 (2021-09-15)
===================

Fixed
-----

* Updated the tests to use databroker.temp instead of sqlite databroker.
* ``xfail`` test that uses removed API.
* Fix ``list_grid_scan`` metadata for ``plan_pattern_args``.
* Fix descriptors for flyers that do not implement ``read_configuration``.

Enhancements
------------

* Do not pass the ``loop`` kwarg to ``RunEngine`` and ``RunBundler`` if we do
  not have to.
* ``RunEngine``'s ``__call__`` now may return plan value, as toggled by new
  ``call_returns_result`` flag.  Default behavior has not changed, but may
  change in a future release.

Added
-----

* Enabled support of Python 3.9 and added it to the test matrix.

Documentation
-------------

* Update TOC links to blueskyproject.io.
* Added release instructions.
* Filled out ``README.md`` and added the ``description`` and
  ``long_description`` fields to ``setup.py``.


v1.7.0 (2021-07-14)
===================

Fixed
-----

* Fixed missing log output for CLI ZMQ proxy.
* Depreciated argument `logfile` of
  :func:`bluesky.commandline.zmq_proxy.start_dispatcher`.
* Better behavior when zmq RemoteDispatcher receives malformed messages.

Enhancements
------------

* Reorganized utils into subpackage, no API changes.
* Added :class:`bluesky.utils.jupyter.NotebookProgressBar`.
* :class:`bluesky.utils.PersistentDict` now inherits from
  :class:`collections.abc.MutableMapping`.
* New module :mod:`bluesky.protocols` designed for type checking devices.
  See PEP 544.


v1.6.7 (2020-11-04)
===================

Fixed
-----

* Tweak layout of plots produced by the Best-Effort Callback when showing
  many LiveGrids.
* The :func:`bluesky.simulators.check_limits` simulator now calls
  ``obj.check_value()`` instead of looking at ``obj.limits``.
* When a document is emitted from a RunEngine, a log message is always issued.
  Previously, Resource and Datum documents were missed.
* Various docstrings were fixed to match the actual function signatures.
* The utility :func:`bluesky.utils.is_movable` for checking with an object
  satifies the expected interfaced for a "movable" object now correctly treats
  the ``stop`` method and ``position`` attribute as optional.
* Documentation about the expected interface for "movable" objects was
  incomplete and has been revised to match reality.

v1.6.6 (2020-08-26)
===================

Fixed
-----

* :class:`bluesky.utils.PersistentDict` has new methods
  :meth:`bluesky.utils.PersistentDict.reload` and
  :meth:`bluesky.utils.PersistentDict.flush` to syncing from and to disk. It
  flushes at garbage collection or system exit, which ensures that any values
  that have been mutated are updated on disk.

v1.6.5 (2020-08-06)
===================

Fixed
-----

* LiveGrid and LiveScatter failed to update

Enhancements
------------

* Expand the class of objects considered "moveable" to include those with expected
  attributes defined as instance attributes

v1.6.4 (2020-07-08)
===================

Fixed
-----

* Allow ``:`` to be used in keynames and still format LiveTable.
* Address use of ``loop`` argument deprecated in Python 3.8.
* Ensure that ``bluesky.utils`` is importable from a background thread. (Do
  not create an instance of `~bluesky.utils.DefaultDuringTask` at import time.)

v1.6.3 (2020-06-25)
===================

Fixed
-----

* Incorrect implementation of :func:`~bluesky.bundlers.RunBundler.collect` has been corrected.

v1.6.2 (2020-06-05)
===================

Fixed
-----

* Missing implementation details of :func:`~bluesky.bundlers.RunBundler.collect` have been added.

v1.6.1 (2020-05-08)
===================

Added
-----

* The plans :func:`~bluesky.plans.grid_scan` and
  :func:`~bluesky.plans.rel_grid_scan` accept a new ``snake_axes`` parameter,
  now matching what :func:`~bluesky.plans.list_grid_scan` and
  :func:`~bluesky.plans.rel_list_grid_scan` do. This can be used to control
  which axes follow a back-and-forth "snake-like" trajectory.

  .. code:: python

     # Default - snaking is disabled
     grid_scan([hw.det], hw.motor, 1, 2, 5, hw.motor1, 7, 2, 10, hw.motor2, 3, 5, 4)

     # Snaking is explicitely disabled
     grid_scan([hw.det], hw.motor, 1, 2, 5, hw.motor1, 7, 2, 10, hw.motor2, 3, 5, 4, snake_axes=False)

     # Snaking can also be disabled by providing empty list of motors
     grid_scan([hw.det], hw.motor, 1, 2, 5, hw.motor1, 7, 2, 10, hw.motor2, 3, 5, 4, snake_axes=[])

     # Snaking is enabled for all motors except the slowest hw.motor
     grid_scan([hw.det], hw.motor, 1, 2, 5, hw.motor1, 7, 2, 10, hw.motor2, 3, 5, 4, snake_axes=True)

     # Snaking is enabled only for hw.motor1
     grid_scan([hw.det], hw.motor, 1, 2, 5, hw.motor1, 7, 2, 10, hw.motor2, 3, 5, 4, snake_axes=[hw.motor1])

     # Snaking is enabled only for hw.motor1 and hw.motor2
     grid_scan([hw.det], hw.motor, 1, 2, 5, hw.motor1, 7, 2, 10, hw.motor2, 3, 5, 4, snake_axes=[hw.motor1, hw.motor2])

  The old (harder to read) way of specifying "snake" parameters, interleaved
  with the other parameters, is still supported for backward-compatibility.

  .. code:: python

     grid_scan([hw.det], hw.motor, 1, 2, 5, hw.motor1, 7, 2, 10, True, hw.motor2, 3, 5, 4, False)

  The two styles---interleaved parameters vs. the new ``snake_axes``
  parameter---cannot be mixed. Mixing them will cause a ``ValueError`` to be
  raised.

Fixed
-----

* Fixed a regression in v1.6.0 which accidentally broke some usages of the
  ``per_step`` parameter in scans.
* The plan :func:`bluesky.plans.fly` returned ``None`` by mistake. It now
  returns the Run Start uid, as do all the other plans that module.

v1.6.0 (2020-03-16)
===================

The most important change in this release is a complete reworking of how
bluesky interacts with the asyncio event loop. This resolves a long-running
issue of bluesky being incompatible with ``tornado >4``, which often tripped up
users in the context of using bluesky from Jupyter notebooks.

There are several other new features and fixes, including new plans and more
helpful error messages, enumerated further below.

Event loop re-factor
--------------------

Previously, the :class:`~bluesky.run_engine.RunEngine` had been repeatedly starting and
stopping the asyncio event loop in :meth:`~bluesky.run_engine.RunEngine.__call__`,
:meth:`~bluesky.run_engine.RunEngine.request_pause`, :meth:`~bluesky.run_engine.RunEngine.stop`, in
:meth:`~bluesky.run_engine.RunEngine.abort`, :meth:`~bluesky.run_engine.RunEngine.halt`, and
:meth:`~bluesky.run_engine.RunEngine.resume`.  This worked, but is bad practice.  It
complicates attempts to integrate with the event loop with other tools.
Further, because as of tornado 5, tornado reports its self as an asyncio event
loop so attempts to start another asyncio event loop inside of a task fails
which means bluesky will not run in a jupyter notebook.  To fix this we now
continuously run the event loop on a background thread and the
:class:`~bluesky.run_engine.RunEngine` object manages the interaction with creating tasks
on that event loop.  To first order, users should not notice this change,
however details of how manage both blocking the user prompt and how we
suspend processing messages from a plan are radically different.
One consequence of running the event loop on a background thread is
that the code in plans and the callbacks is executed in that thread as well.
This means that plans and callbacks must now be threadsafe.

API Changes
~~~~~~~~~~~

``install_qt_kicker`` deprecated
++++++++++++++++++++++++++++++++

Previously, we were running the asyncio event loop on the main thread
and blocked until it returned.  This meant that if you were using
Matplotlib and Qt for plots they would effectively be "frozen" because
the Qt event loop was not being given a chance to run.  We worked
around this by installing a 'kicker' task onto the asyncio event loop
that would periodically spin the Qt event loop to keep the figures
responsive (both to addition of new data from callbacks and from user
interaction).

Now that we are running the event loop on a background thread this no
longer works because the Qt event loop must be run on the main thread.
Instead we use *during_task* to block the main thread by running the
Qt event loop directly.


``during_task`` kwarg to ``RunEngine.__init__``
+++++++++++++++++++++++++++++++++++++++++++++++

We need to block the main thread in :meth:`~bluesky.run_engine.RunEngine.__call__` (and
:meth:`~bluesky.run_engine.RunEngine.resume`) until the user supplied plan is complete.
Previously, we would do this by calling ``self.loop.run_forever()`` to
start the asyncio event loop.  We would then stop the event loop an
the bottom of ``RunEngine._run`` and in
:meth:`~bluesky.run_engine.RunEngine.request_pause` to un-block the main thread and return
control to the user terminal.  Now we must find an alternative way to achieve
this effect.

There is a a :class:`threading.Event` on the :class:`~bluesky.run_engine.RunEngine` that
will be set when the task for ``RunEngine._run`` in completed,
however we can not simple wait on that event as that would again cause the Qt
windows to freeze.  We also do not want to bake a Matplotlib / Qt dependency
directly into the :class:`~bluesky.run_engine.RunEngine` so we added a hook, set at init
time, for an object expected to implement the method ``block(event)``.
While the RunEngine executes a plan, it is passed the :class:`threading.Event`
and is responsible for blocking until the Event is set.  This function can do
other things (such as run the Qt event loop) during that time.  The required
signature is ::

  def block(ev: Threading.Event) -> None:
      "Returns when ev is set"


The default hook will handle the case of the Matplotilb Qt backend and
the case of Matplotlib not being imported.


``'wait_for'`` Msg now expects callables rather than futures
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Messages are stashed and re-run when plans are interrupted which would
result in re-using the coroutines passed through.  This has always
been broken, but due to the way were stopping the event loop to pause
the scan it was passing tests.

Instead of directly passing the values passed into :func:`asyncio.wait`, we
now expect that the iterable passed in is callables with the signature::

  def fut_fac() -> awaitable:
      'This must work multiple times'

The persistent dict used by ``RE.md`` must be thread-safe
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

By default, ``RE.md`` is an ordinary dictionary, but any dict-like object may
be used. It is often convenient for the contents of that dictionary to persist
between sessions. To achieve this, we formerly recommended using
``~historydict.HistoryDict``. Unfortunately,
``~historydict.HistoryDict`` is not threadsafe and is not compatible with
bluesky's new concurrency model. We now recommend using
:class:`~bluesky.utils.PersistentDict`. See :ref:`md_persistence` for
instructions on how to migrate existing metadata.

Callbacks must be thread-safe
+++++++++++++++++++++++++++++

Because callbacks now run on the background thread they must be
thread-safe.  The place where this is most likely to come up is in the
context of plotting which generally creates a GUI window.  Almost all
GUI frameworks insist that they only be interacted with only on the
main thread.  In the case of Qt we provide
:class:`~bluesky.callbacks.mpl_plotting.QtAwareCallback` to manage
moving Qt work back to the main thread (via a Qt ``Signal``).


Plans must be thread-safe
+++++++++++++++++++++++++

Because the plans now execute on the background thread they must be
thread-safe if the touch any external state.  Similarly the callbacks,
we expect that the most likely place for this to fail is with
plotting.  In most cases this can be addressed by using a thread-safe
version of the callback.


Features
--------

* Added support for :doc:`multi_run_plans`.
* Added better logging and convenience functions for managing it more easily.
  See :doc:`debugging`.
* Generalized :func:`~bluesky.plans.list_scan` to work on any number of motors,
  not just one. In v1.2.0, :func:`~bluesky.plans.scan` was generalized in the
  same way.
* Added :func:`~bluesky.plans.list_grid_scan`.
* Added :func:`~bluesky.plan_stubs.rd`.
* Added :class:`~bluesky.suspenders.SuspendWhenChanged`.
* Added :func:`~bluesky.callbacks.core.make_callback_safe` and
  :func:`~bluesky.callbacks.core.make_class_safe`.
* Added a ``per_shot`` parameter to :func:`bluesky.plans.count`, analogous to
  the ``per_step`` parameter supported by plans that do scans.
* Accept ``**kwargs`` to :func:`~bluesky.plan_stubs.mv` and
  :func:`~bluesky.plan_stubs.mvr`. Pass them through to all motors involved in
  the move. Notably, this allows plans to pass a ``timeout`` parameter through
  the ``obj.set()``.
* Added a new built-in RunEngine command, ``RE_class``, which sends the type of
  the ``RunEngine`` into the generator. This allows the plan to know if it is
  being consumed by the usual ``RunEngine``, a subclass, or some
  non-responsive consumer like ``list``.
* Raise a more helpful error message if the ``num`` parameter given to
  :func:`~bluesky.plans.scan` is not a whole number, as can happen if ``num`` is
  mistaken to mean "step size".
* Report the version of bluesky and (if available) ophyd in the metadata.
* Add a more helpful error message if the value returned from some call to
  ``obj.read()`` returns ``None`` instead of the expected dict.
* If the user tries to start a :class:`~bluesky.callbacks.zmq.RemoteDispatcher`
  after it has been stopped, raise a more helpful error message.

Bug Fixes
---------

* The ``state`` attribute of the ``RunEngine`` is now a read-only property, as
  it should have always been.
* In the Best-Effort Callback, do not assume that the RunStart document
  includes ``'scan_id'``, which is an optional key.
* The commandline utility ``bluesky-0MQ-proxy`` now works on Windows.
* The IPython integrations have been updated for compatibility with IPython 7.
* Added support for "adaptive fly scans" by enabling the ``'collect'`` message
  to (optionally) return the Events it emitted.
* Fixed bug in tqdm-based progress bar where tqdm could be handed a value it
  considered invalid.

Other API Changes
-----------------

* Removed attribute ``nnls`` from
  :class:`bluesky.callbacks.best_effort.PeakResults`. It has always been
  ``None`` (never implemented) and only served to cause confusion.

v1.5.7 (2020-05-01)
===================

Bug Fixes
---------

This release fixes a bug that resulted in no configuration data related
to fly scans being added to descriptors.


v1.5.6 (2020-03-11)
===================

Added support for Python 3.8 and the following for forward-compatibility with
1.6.0.

* :class:`bluesky.utils.PersistentDict`
* :class:`bluesky.callbacks.mpl_plotting.QtAwareCallback`

See
`the 1.5.6 GH milestone <https://github.com/bluesky/bluesky/milestone/19?closed=1>`_
for the complete list of changes.

v1.5.5 (2019-08-16)
===================

Support fix ``bluesky.utils.register_transform`` with IPython >= 7


v1.5.4 (2019-08-09)
===================

Support Maplotlib 3.1 and above. (Do not use deprecated and removed aspect
adjustable values.)

v1.5.3 (2019-05-27)
===================

This release removes the dependency on an old version of the ``jsonschema``
library and requires the latest version of the ``event-model`` library.


v1.5.2 (2019-03-11)
===================

This release fixes compatibility with matplotlib 2.x; at least some matplotlib
2.x releases are not compatible with the matplotlib plotting callbacks in
bluesky v1.5.1. This release of bluesky is compatible with all 2.x and 3.x
releases.

v1.5.1 (2019-03-08)
===================

This release contains bug fixes and documentation updates.

Features
--------

* Use the ISO8601 delimiters for date in RE scans.

Bug Fixes
---------

* Pin jsonschema <3 due to its deprecations.
* Stop using deprecated API in Matplotlib.


v1.5.0 (2019-01-03)
===================

This release includes many documentation fixes and handful of new features,
especially around improved logging.

Features
--------

* Logging has been increased and improved.
* A default handler is added to the ``'bluesky'`` logger at import time. A new
  convenience function, :func:`~bluesky.set_handler`, addresses common cases
  such as directing the log output to a file.
* The ``bluesky-0MQ-proxy`` script now supports a ``-v, --verbose`` option,
  which logs every start and stop document received and a ``-vvv`` ("very
  verbose") option, which logs every document of every type.
* The prefix on messages sent by :class:`bluesky.callbacks.zmq.Publisher` can
  be set to arbitrary bytes. (In previous versions, the prefix was hardcoded to
  an encoded combination of the hostname, process ID, and the Python object ID
  of a RunEngine instance.)
* The RunEngine includes a human-readable, not-necessarily-unique ``scan_id``
  key in the RunStart document. The source of the ``scan_id`` is now pluggable
  via a new parameter, ``scan_id_source``. See :doc:`run_engine_api` for
  details.
* The convenience function, :func:`bluesky.utils.ts_msg_hook` accepts new
  parameter ``file`` for directing the output to a file instead of the standard
  out.
* It is possible to use those callbacks that do not require matplotlib without
  importing it.

Bug Fixes
---------

* Fixed BestEffortCallback's handling of integer data in plots.
* Fixed invalid escape sequence that produced a warning in Python 3.6.

Breaking Changes
----------------

* The signature of :class:`bluesky.callbacks.zmq.RemoteDispatcher` has been
  changed in a non-backward-compatible way. The parameters for filtering
  messages by ``hostname``, ``pid``, and ``run_engine_id`` have been replaced
  by one new parameter, ``prefix``.
* The default value of ``RunEngine.verbose`` is now ``True``, meaning that the
  ``RunEngine.log`` is *not* disabled by default.

Deprecations
------------

* The :class:`bluesky.callbacks.zmq.Publisher` accepts an optional RunEngine
  instance, which the Publisher subscribes to automatically. This parameter has
  been deprecated; users are now encouraged to subscribe the publisher to the
  RunEngine manually, in the normal way (``RE.subscribe(publisher)``). The
  parameter may be removed in a future release of bluesky.

v1.4.1 (2018-09-24)
===================

This release fixes a single regression introduced in v1.4.0. We recommend all
users upgrade.

Bug Fixes
---------

* Fix a critical typo that made
  :class:`~bluesky.callbacks.mpl_plotting.LiveGrid` unusable.

Note that the 1.4.x series is not compatible with newer versions of matplotlib;
it needs a version lower than 3.1.0 due to an API change in matplotlib. The
1.5.x series is compatible with matplotlib versions before and after the
change.

v1.4.0 (2018-09-05)
===================

Features
--------

* Added ability to control 'sense' of
  :class:`~bluesky.callbacks.mpl_plotting.LiveGrid` (ex "positive goes
  down and to the right") to match the coordinates in the hutch.
* Learned how to specify the serializer / deserializer for the zmq
  publisher / client.
* Promoted the inner function from :func:`~bluesky.plan_stubs.one_nd_step`
  to a top-level plan :func:`bluesky.plan_stubs.move_per_step`.
* Added flag to :func:`~bluesky.plans.ramp_plan` to control if a
  data point is taken before the ramp starts.

Bug Fixes
---------

* Ensure order stability in :func:`~bluesky.magics.get_labeled_devices`
  on all supported versions of Python.
* Fixed typos, dev requirements, and build details.


v1.3.3 (2018-06-06)
===================

Bug Fixes
---------

* Fixed show-shopping RunEngine bug in flyer asset collection. (The impact of
  this bug is expected to be low, as there *are* no flyers with asset
  collection yet and the bug was discovered while writing the first one.)
* Fixed packaging issue where certain important files (notably
  ``requirements.txt``) were not included in the source tarball.
* Made BestEffortCallback swallow errors related to matplotlib's "tight layout"
  if the occur --- better to show a messy plot than error out.

v1.3.2 (2018-05-24)
===================

Bug Fixes
---------

* Revised behavior of magics that integrate with ophyd's experimental
  "labels" feature. The most important difference is that the ``%wa`` magic now
  traverses the children of labeled devices to find any sub-devices that are
  positioners.

v1.3.1 (2018-05-19)
===================

Bug Fixes
---------

* Fixed race condition where monitored signals could emit an Event document
  before the corresponding Event Descriptor document.
* Addressed incompatibilities with upcoming release of Python, 3.7.

v1.3.0 (2018-05-15)
===================

Features
--------

* When used with ophyd v1.2.0 or later, emit Resource and Datum documents
  through the RunEngine. Previously, ophyd would insert these documents
  directly into a database. This left other consumers with only partial
  information (for example, missing file paths to externally-stored data) and
  no guarantees around synchronization. Now, ophyd need not interact with a
  database directly. All information flows through the RunEngine and out to any
  subscribed consumers in a deterministic order.
* New Msg commands, ``install_suspender`` and ``remove_suspender``, allow plans
  to temporarily add and remove Suspenders.
* The RunEngine's signal handling (i.e. Ctrl+C capturing) is now configurable.
  The RunEngine accepts a list of ``context_managers`` that it will enter and
  exit before and after running. By default, it has one context manager that
  handles Ctrl+C. To disable Ctrl+C handling, pass in an empty list instead.
  This can also be used to inject other custom behavior.
* Add new plans: :func:`~bluesky.plans.x2x_scan`,
  :func:`~bluesky.plans.spiral_square_plan`, and
  :func:`~bluesky.plans.rel_spiral_square_plan`.
* Add convenience methods for reviewing the available commands,
  :meth:`~bluesky.run_engine.RunEngine.commands` and
  :meth:`~bluesky.run_engine.RunEngine.print_command_registry`.
* Add a ``crossings`` attribute to ``PeakStats``.

Bug Fixes
---------

* When resuming after a suspender, call ``resume()`` on all devices (if
  present).
* Fixed BEC LiveGrid plot for a motor with one step.
* A codepath in ``LiveFit`` that should have produced a warning produced an
  error instead.

Breaking Changes
----------------

* User-defined callbacks subscribed to the RunEngine ``'all'`` stream must
  accept documents with names ``'resource'``, ``'datum'`` and ``'bulk_datum'``.
  It does not necessarily have to heed their contents, but it must not fall
  over if it receives one.

Deprecations
------------

* The IPython "magics", always marked as experimental, have been reworked.
  Instead of relying on the singleton lists, ``BlueskyMagics.positioners`` and
  ``BlueskyMagics.detectors``, the magics now scrape the user namespace for
  objects that implement the ``_ophyd_labels_`` interface. See :doc:`magics`
  for the new usage. The magics will revert to their old behavior if the
  singleton lists are non-empty, but they will produce a warning. The old
  behavior will be removed in a future release.

v1.2.0 (2018-02-20)
===================

Features
--------

* Refreshed documentation with a new :doc:`tutorial` section.
* Extend :func:`.scan` and :func:`.rel_scan` to
  handle multiple motors, rendering :func:`.inner_product_scan` and
  :func:`relative_inner_product_scan` redundant.
* A new plan stub, :func:`~bluesky.plan_stubs.repeat`, repeats another plan N
  times with optional interleaved delays --- a kind of customizable version of
  :func:`~bluesky.plans.count`.
* Better validation of user-defined ``per_step`` functions and more informative
  error messages to match.

Bug Fixes
---------

* Fix axes orientation in :class:`.LiveRaster`.
* Make :class:`.BestEffortCallback` display multi-motor scans properly.
* Fix bug in :func:`.ts_msg_hook` where it conflated month and minute. Also,
  include sub-second precision.
* Avoid situation where plan without hints caused the
  :class:`.BestEffortCallback` to error instead of do its best to guess useful
  behavior.
* Skip un-filled externally-stored data in :class:`.LiveTable`. This fixes a
  bug where it is expecting array data but gets UUID (``datum_id``) and errors
  out.

Deprecations
------------

* The :func:`~bluesky.plan_stubs.caching_repeater` plan has been deprecated
  because it is incompatible with some preprocessors. It will be removed in
  a future release of bluesky. It was not documented in any previous releases
  and rarely if ever used, so the impact of this removal is expected to be low.

v1.1.0 (2017-12-19)
===================

This release fixes small bugs in v1.0.0 and introduces one new feature. The
API changes or deprecations are not expected to affect many users.

Features
--------

* Add a new command to the :class:`~bluesky.run_engine.RunEngine`, ``'drop'``,
  which jettisons the currently active event bundle without saving. This is
  useful for workflows that generate many readings that can immediately be
  categorized as not useful by the plan and summarily discarded.
* Add :func:`~bluesky.utils.install_kicker`, which dispatches automatically to
  :func:`~bluesky.utils.install_qt_kicker` or
  :func:`~bluesky.utils.install_nb_kicker` depending on the current matplotlib
  backend.

Bug Fixes
---------

* Fix the hint for :func:`~bluesky.plans.inner_product_scan`, which previously
  used a default hint that was incorrect.

Breaking Changes and Deprecations
---------------------------------

* In :func:`~bluesky.plans.tune_centroid`, change the meaning of the
  ``step_factor`` parameter to be the factor to reduce the range of each
  successive iteration. Enforce bounds on the motion, and determine the
  centroid from each pass separately.
* The :class:`~bluesky.preprocessors.SupplementalData` preprocessor inserts its
  instructions in a more logical order: first baseline readings, then
  monitors, then flyers. Previously, the order was reversed.
* The suspender :class:`~bluesky.suspenders.SuspendInBand` has been renamed to
  :class:`~bluesky.suspenders.SuspendWhenOutsideBand` to make its meaning more
  clear. Its behavior has not changed: it suspends when a value exits a given
  range. The original, confusing name now issues a warning.
* The suspender :class:`~bluesky.suspenders.SuspendOutBand`, which
  counter-intuitively suspends *when a value enters a given range*, has been
  deprecated. (If some application is found for this unusual scenario, the user
  can always implement a custom suspender to handle it.)

v1.0.0 (2017-11-07)
===================

This tag marks an important release for bluesky, signifying the conclusion of
the early development phase. From this point on, we intend that this project
will be co-developed between multiple facilities. The 1.x series is planned to
be a long-term-support release.

Bug Fixes
---------

* :func:`~bluesky.plan_stubs.mv` and :func:`~bluesky.plan_stubs.mvr` now works
  on pseudopositioners.
* :func:`~bluesky.preprocessors.reset_positions_wrapper` now works on
  pseudopositioners.
* Plans given an empty detectors list, such as ``count([])``, no longer break
  the :class:`~bluesky.callbacks.best_effort.BestEffortCallback`.

v0.11.0 (2017-11-01)
====================

This is the last release before 1.0.0. It contains major restructurings and
general clean-up.

Breaking Changes and Deprecations
---------------------------------

* The :mod:`bluesky.plans` module has been split into

    * :mod:`bluesky.plans` --- plans that create a run, such as :func:`count`
      and :func:`scan`
    * :mod:`bluesky.preprocessors` --- plans that take in other plans and
      motify them, such as :func:`baseline_wrapper`
    * :mod:`bluesky.plan_stubs` --- small plans meant as convenient building
      blocks for creating custom plans, such as :func:`trigger_and_read`
    * :mod:`bluesky.object_plans` and :mod:`bluesky.cntx`, containing
      legacy APIs to plans that were deprecated in a previous release and
      will be removed in a future release.

* The RunEngine raises a ``RunEngineInterrupted`` exception when interrupted
  (e.g. paused). The optional argument ``raise_if_interrupted`` has been
  removed.
* The module :mod:`bluesky.callbacks.scientific` has been removed.
* ``PeakStats`` has been moved to :mod:`bluesky.callbacks.fitting`, and
  :func:`plot_peak_stats` has been moved to `bluesky.callbacks.mpl_plotting`.
* The synthetic 'hardware' objects in ``bluesky.examples`` have been relocated
  to ophyd (bluesky's sister package) and aggressively refactored to be more
  closely aligned with the behavior of real hardware. The ``Reader`` and
  ``Mover`` classes have been removed in favor of ``SynSignal``,
  ``SynPeriodicSignal``, ``SynAxis``, ``SynSignalWithRegistry``.

Features
--------

* Add :func:`stub_wrapper` and :func:`stub_decorator` that strips
  open_run/close_run and stage/unstage messages out of a plan, so that it can
  be reused as part of a larger plan that manages the scope of a run manually.
* Add :func:`tune_centroid` plan that iteratively finds the centroid of a
  single peak.
* Allow devices with couple axes to be used in N-dimensional scan plans.
* Add :func:`contingency_wrapper` and :func:`contingency_decorator` for
  richer cleanup specification.
* The number of events in each event stream is recorded in the RunStop document
  under the key 'num_events'.
* Make the message shown when the RunEngine is paused configurable via the
  attribute ``RunEngine.pause_msg``.

Bug Fixes
---------

* Fix ordering of dimensions in :func:`grid_scan` hints.
* Show Figures created internally.
* Support a negative direction for adaptive scans.
* Validate that all descriptors with a given (event stream) name have
  consistent data keys.
* Correctly mark ``exit_status`` field in RunStop metadata based on which
  termination method was called (abort, stop, halt).
* ``LiveFitPlot`` handles updates more carefully.

Internal Changes
----------------

* The :mod:`bluesky.callbacks` package has been split up into more modules.
  Shim imports maintain backward compatibility, except where noted in the
  section on API Changes above.
* Matplotlib is now an optional dependency. If it is not importable,
  plotting-related callbacks will not be available.
* An internal change to the RunEngine supports ophyd's new Status object API
  for adding callbacks.

v0.10.3 (2017-09-12)
====================

Bug Fixes
---------

* Fix critical :func:`baseline_wrapper` bug.
* Make :func:`plan_mutator` more flexible. (See docstring.)

v0.10.2 (2017-09-11)
====================

This is a small release with bug fixes and UI improvements.

Bug Fixes
---------

* Fix bug wherein BestEffortCallback tried to plot strings as floats. The
  intended behavior is to skip them and warn.

Features
--------

* Include a more informative header in BestEffortCallback.
* Include an 'Offset' column in %wa output.

v0.10.1 (2017-09-11)
====================

This release is equivalent to v0.10.2. The number was skipped due to packaging
problems.

v0.10.0 (2017-09-06)
====================

Highlights
----------

* Automatic best-effort visualization and peak-fitting is available for all
  plans, including user-defined ones.
* The "SPEC-like" API has been fully removed, and its most useful features have
  been applied to the library in a self-consistent way. See the next section
  for detailed instructions on migrating.
* Improved tooling for streaming documents over a network for live processing
  and visualization in a different process or on a different machine.

Breaking Changes
----------------

* The modules implementing what was loosely dubbed a "SPEC-like" interface
  (``bluesky.spec_api`` and ``bluesky.global_state``) have been entirely
  removed. This approach was insufficently similar to SPEC to satisfy SPEC
  users and confusingly inconsistent with the rest of bluesky.

  The new approach retains the good things about that interface and makes them
  available for use with *all* plans consistently, including user defined ones.
  Users who have been fully utilitzing these "SPEC-like" plans will notice four
  differences.

  1. No ``gs.DETS``. Just use your own variable for detectors. Instead of:

     .. code-block:: python

         # OLD ALTERNATIVE, NO LONGER SUPPORTED

         from bluesky.global_state import gs
         from bluesky.spec_api import ct

         gs.DETS = # a list of some detectors
         RE(ct())

     do:

     .. code-block:: python

        from bluesky.plans import count

        dets = # a list of some detectors
        RE(count(dets))

     Notice that you can use multiple lists to enable easy task switching.
     Instead of continually updating one global list like this:

     .. code-block:: python

         # OLD ALTERNATIVE, NO LONGER SUPPORTED

         gs.DETS = # some list of detectors
         RE(ct())

         gs.DETS.remove(some_detector)
         gs.DETS.append(some_other_detector)
         RE(ct())

     you can define as many lists as you want and call them whatever you want.

     .. code-block:: python

        d1 = # a list of some detectors
        d2 = # a list of different detectors
        RE(count(d1))
        RE(count(d2))

  2. Automatic baseline readings, concurrent monitoring, and "flying"
     can be set up uniformly for all plans.

     Formerly, a list of devices to read at the beginning and the end of each
     run ("baseline" readings), a list of signals to concurrent monitor, and
     a list of "flyers" to run concurrently were configured like so:

     .. code-block:: python

        # OLD ALTERNATIVE, NO LONGER SUPPORTED

        from bluesky.spec_api import ct

        gs.BASELINE_DEVICES = # a list of devices to read at start and end
        gs.MONTIORS = # a list of signals to monitor concurrently
        gs.FLYERS = # a list of "flyable" devices

        gs.DETS = # a list of detectors

        RE(ct())  # monitoring, flying, and baseline readings are added

     And formerly, those settings only affected the behavior of the "SPEC-like"
     plans, such as ``ct`` and ``ascan``. They were ignored by their
     counterparts ``count`` and ``scan``, as well as user-defined plans. This
     was not desirable!

     This scheme has been replaced by the
     :ref:`supplemental data <supplemental_data>`, which can be
     used to globally modify *all* plans, including user-defined ones.

     .. code-block:: python

        from bluesky.plans import count

        # one-time configuration
        from bluesky import SupplementalData
        sd = SupplementalData()
        RE.preprocessors.append(sd)

        # interactive use
        sd.monitors = # a list of signals to monitor concurrently
        sd.flyers = # a list of "flyable" devices
        sd.baseline = # a list of devices to read at start and end

        dets = # a list of detectors
        RE(count(dets))  # monitoring, flying, and baseline readings are added

  3. Automatic live visualization and peak analysis can be set up uniformly for
     all plans.

     Formerly, the "SPEC-like" plans such as ``ct`` and ``ascan`` automatically
     set up a suitable table and a plot, while their "standard" vanilla
     counterparts, :func:`bluesky.plans.count` and :func:`bluesky.plans.scan`
     required explicit, detailed instructions to do so. Now, a best-effort
     table and plot can be made for *all* plans, including user-defined ones,
     by invoking this simple configuration:

     .. code-block:: python

        from bluesky.plans import count

        # one-time configuration
        from bluesky.callbacks.best_effort import BestEffortCallback
        bec = BestEffortCallback()
        RE.subscribe(bec)

        # interactive use
        dets = # a list of detectors
        RE(count(dets), num=5))  # automatically prints table, shows plot

     Use ``bec.disable()`` and ``bec.enable()`` to temporarily toggle the
     output off and on.

  4. Peak anallysis, now computed automatically by the BestEffortCallback
     above, can be viewed with a keyboard shortcut. The peak statistics,
     formerly encapsulated in ``gs.PS``, are now organized differently.

     For each plot, simple peak-fitting is performed in the background. Of
     course, it may or may not be applicable depending on your data, and it is
     not shown by default. To view fitting annotations in a plot, click the
     plot area and press Shift+P. (Lowercase p is a shortcut for
     "panning" the plot.)

     To access the peak-fit statistics programmatically, use ``bec.peaks``. For
     convenience, you may alias this like:

     .. code-block:: python

        peaks = bec.peaks

     Inside ``peaks``, access various statistics like:

     .. code-block:: python

        peaks.com
        peaks.cen
        peaks.max
        peaks.min

     Each of these is a dictionary with an entry for each field that was fit.
     For example, the 'center of mass' peak statistics for a field named
     ``'ccd_stats1_total'`` would be accessed like
     ``peaks.com['ccd_stats1_total']``.
* The functions and classes in the module ``bluesky.callbacks.broker`` require
  a instance of ``Broker`` to be passed in as an argument. They used to default
  to the 'singleton' instance via ``from databroker import db``, which is now a
  deprecated usage in databroker.
* The plan preprocessors ``configure_count_time_wrapper`` and
  ``configure_count_time_decorator`` were moved to ``bluesky.plans`` from
  ``bluesky.spec_api``, reverting a change made in v0.9.0.
* The 0MQ pubsub integration classes ``Publisher`` and ``RemoteDispatcher``
  have been overhauled. They have been moved from
  :mod:`bluesky.callbacks.zmqpub` and :mod:`bluesky.callbacks.zmqsub` to
  :mod:`bluesky.callbacks.zmq` and their signatures have been changed to match
  similar utilities in the pydata ecosystem. See the Enhancements section for
  more information.
* The module ``bluesky.qt_kicker`` has been removed. Its former contents are
  avaiable in ``bluesky.utils``. The module was originally deprecated in April
  2016, and it has been issuing warnings about this change since.
* The plan ``bluesky.plans.input`` has been renamed
  ``bluesky.plans.input_plan`` to avoid shadowing a builtin if the module is
  bulk-imported. The plan was previously undocumented and rarely used, so the
  impact of this change on users is expected to be small.

Deprecations
------------

* The module :mod:`bluesky.plan_tools` has been renamed
  :mod:`bluesky.simualtors`.  In the new module,
  :func:`bluesky.plan_tools.print_summary`` has been renamed
  :func:`bluesky.simulators.summarize_plan`.
  The old names are supported in this release, with a warning, but will be
  removed in a future release.
* The Object-Orientated plans (``Count``, ``Scan``, etc.) have been deprecated
  and will be removed in a future release. Their documentation has been
  removed.
* The plan context managers (``run_context``, ``stage_context``, etc.) have
  been deprecated and will be removed in a future release. They were never
  documented or widely used.
* The method :meth:`bluesky.Dispatcher.subscribe` (which is encapsulated into
  :class:`bluesky.run_engine.RunEngine` and inherited by
  :class:`bluesky.callbacks.zmq.RemoteDispatcher`) has a new signature. The
  former signature was ``subscribe(name, func)``. The new signature is
  ``subscribe(func, name='all')``. Because the meaning of the arguments is
  unambigious (they must be a callable and a string, respectively) the old
  order will be supported indefeinitely, with a warning.

Features
--------

* A :doc:`progress bar <progress-bar>` add-on is available.
* As addressed above:
    * The new :ref:`supplemental data <supplemental_data>` feature make it
      easy to set up "baseline" readings and asynchronous acquisition in a way
      that applies automatically to all plans.
    * The new :ref:`best-effort callback <best_effort_callback>` sets up
      automatic terminal output and plots for all plans, including user-defined
      ones.
* ``LivePlot`` now accepts ``x='time'``. It can set t=0 to the UNIX epoch or to
  the start of the run. It also accepts ``x='seq_num'``---a synonym for
  ``x=None``, which remains the default.
* A new simulator :func:`bluesky.simulators.check_limits` verifies that a plan
  will not try to move a movable device outside of its limits.
* A new plan, :func:`bluesky.plan.mvr`, has been added as a relative counterpart
  to :func:`bluesky.plan.mv`.
* The 0MQ pubsub integration classes :class:`bluesky.callbacks.zmq.Publisher``
  and :class:`bluesky.callbacks.zmq.RemoteDispatcher` have been simplified.
  A new class :class:`bluesky.callbacks.zmq.Proxy` and command-line utility
  ``bluesky-0MQ-proxy`` has been added to streamline configuration.
* Metadata recorded by many built-in plans now includes a new item,
  ``'hints'``, which is used by the best-effort callback to produce useful
  visualizations. Additionally, the built-in examples devices have
  :ref:`a new hints attribute <hints>`. Its content may change or expand in
  future releases as this new feature is explored.
* Some :doc:`IPython magics <magics>` mimicing the SPEC API have been added.
  These are experimental and may be altered or removed in the future.

Bug Fixes
---------

* Using the "fake sleep" feature of simulated Movers (motors) caused them to
  break.
* The ``requirements.txt`` failed to declare that bluesky requires matplotlib.

v0.9.0 (2017-05-08)
===================

Breaking Changes
----------------

* Moved ``configure_count_time_wrapper`` and
  ``configure_count_time_decorator`` to ``bluesky.spec_api`` from
  ``bluesky.plans``.
* The metadata reported by step scans that used to be labeled ``num_steps``
  is now renamed ``num_points``, generally considered a less ambiguous name.
  Separately, ``num_interals`` (which one might mistakenly assume is what was
  meant by ``num_steps``) is also stored.


v0.8.0 (2017-01-03)
===================

Features
--------

* If some plan or callback has hung the RunEngine and blocked its normal
  ability to respond to Ctrl+C by pausing, it is not possible to trigger a
  "halt" (emergency stop) by hammering Ctrl+C more than ten times.

Bug Fixes
---------

* Fix bug where failed or canceled movements could cause future executions of
  the RunEngine to error.
* Fix bug in ``plan_mutator`` so that it properly handles return values. One
  effect of this fix is that ``baseline_wrapper`` properly passed run uids
  through.
* Fix bug in ``LiveFit`` that broke multivariate fits.
* Minor fixes to example detectors.

Breaking Changes
----------------

* A ``KeyboardInterrupt`` exception captured during a run used to cause the
  RunEngine to pause. Now it halts instead.

v0.7.0 (2016-11-01)
===================

Features
--------

* Nonlinear least-squares minimization callback ``LiveFit`` with
  ``LiveFitPlot``
* Added ``RunEngine.clear_suspenders()`` convenience method.
* New ``RunEngine.preprocessors`` list that modifies all plans passed to the
  RunEngine.
* Added ``RunEngine.state_hook`` to monitor state changes, akin to ``msg_hook``.
* Added ``pause_for_debug`` options to ``finalize_wrapper`` which allows pauses
  the RunEngine before performing any cleanup, making it easier to debug.
* Added many more examples and make it easier to create simulated devices that
  generate interesting simulated data. They have an interface closer to the
  real devices implemented in ophyd.
* Added ``mv``, a convenient plan for moving multiple devices in parallel.
* Added optional ``RunEngine.max_depth`` to raise an error if the RunEngine
  thinks it is being called from inside a function.

Bug Fixes
---------

* The 'monitor' functionality was completely broken, packing configuration
  into the wrong structure and starting seq_num from 0 instead of 1, which is
  the (regrettable) standard we have settled on.
* The RunEngine coroutines no longer mutate the messages they receive.
* Fix race condition in ``post_run`` callback.
* Fix bugs in several callbacks that caused them not to work on saved documents
  from the databroker. Also, make them call ``super()`` to play better with
  multiple inheritance in user code.


Breaking Changes
----------------

* The flag ``RunEngine.ignore_callback_exceptions`` now defaults to False.
* The plan ``complete``, related to fly scans, previously had ``wait=True`` by
  default, although its documentation indicated that ``False`` was the default.
  The code has been changed to match the documentation. Any calls to
  ``complete`` that are expected to be blocking should be updated with the
  keyword ``wait=True``.
* Completely change the API of ``Reader`` and ``Mover``, the classes for
  definding simulated devices.
* The bluesky interface now expects the ``stop`` method on a device to accept
  an optional ``success`` argument.
* The optional, undocumented ``fig`` argument to ``LivePlot`` has been
  deprecated and will be removed in a future release.  An ``ax`` argument has
  been added. Additionally, the axes used by ``LiveGrid`` and ``LiveScatter`` is
  configurable through a new, optional ``ax`` argument.
* The "shortcut" where mashing Ctrl+C three times quickly ran ``RE.abort()``
  has been removed.
* Change the default stream name for monitors to ``<signal_name>_monitor`` from
  ``signal_name>-monitor`` (underscore vs. dash). The impact of this change is
  minimal because, as noted above, the monitor functionality was completely
  broken in previous releases.

v0.6.4 (2016-09-07)
===================

Features
--------

* Much-expanded and overhauled documentation.
* Add ``aspect`` argument to ``LiveGrid``.
* Add ``install_nb_kicker`` to get live-updating matplotlib figures in the
  notebook while the RunEngine is running.
* Simulated hardware devices ``Reader`` and ``Mover`` can be easily customized
  to mock a wider range of behaviors, for testing and demos.
* Integrate the SPEC API with mew global state attribute ``gs.MONITORS``.
* Callbacks that use the databroker accept an optional ``Broker`` instance
  as an argument.

Bug Fixes
---------

* Minor fix in the tilt computation for spiral scans.
* Expost 'tilt' option through SPEC-like API
* The "infinite count" (``ct`` with ``num=None``) should spawn a LivePlot.
* ``finalize_decorator`` accepts a callable (e.g., generator function)
  and does not accept an iterable (e.g., generator instance)
* Restore ``gs.FLYERS`` integration to the SPEC API (accidentally removed).

Breaking Changes
----------------

* The API for the simulated hardware example devices ``Reader`` and ``Mover``
  has been changed to make them more general.
* Remove ``register_mds`` metadatastore integration.

v0.6.3 (2016-08-16)
===================

Features
--------

* Change how "subscription factories" are handled, making them configurable
  through global state.
* Make PeakStats configurable through global state.
* Add an experimental utility for passing documents over a network and
  processing them on a separate process or host, using 0MQ.
* Add ``monitor_during_wrapper`` and corresponding decorator.
* Add ``stage_wrapper`` and corresponding decorator.
* Built-in plans return the run uid that they generated.
* Add a new ``ramp_plan`` for taking data while polling the status of a
  movement.

Bug Fixes
---------

* Boost performance by removing unneeded "sleep" step in message processing.
* Fix bug related to rewinding in preparation for resuming.

Breaking Changes
----------------

* Remove the ``planify`` decorator and the plan context managers. These were
  experimental and ultimately proved problematic because they make it difficult
  to pass through return values cleanly.
* Remove "lossy" subscriptions feature, rendered unnecessary by the utility for
  processing documents in separate processes (see Enhancements, above).

v0.6.2 (2016-07-26)
===================

Bug Fixes
---------

* Make ``make_decorator`` return proper decorators. The original implementation
  returned functions that could not actually be used as decorators.

v0.6.1 (2016-07-25)
===================

This release contained only a minor UX fix involving more informative error
reporting related to Area Detector plugin port configuration.

v0.6.0 (2016-07-25)
===================

Features
--------

* Address the situation where plan "rewinding" after a pause or suspension
  interacted badly with some devices. There are now three ways to temporarily
  turn off rewinding: a Msg with a new 'rewindable' command; a special
  attribute on the device that the ``trigger_and_read`` plan looks for;
  and a special exception that devices can raise when their ``pause`` method
  is called. All three of these features should be considered experimental.
  They will likely be consolidated in the future once their usage is tested
  in the wild.
* Add new plan wrappers and decorators: ``inject_md_wrapper``, ``run_wrapper``,
  ``rewindable_wrapper``.

Bug Fixes
---------

* Fix bug where RunEngine was put in the "running" state, encountered an
  error before starting the ``_run`` coroutine, and thus never switch back to
  "idle."
* Ensure that plans are closed correctly and that, if they fail to close
  themselves, a warning is printed.
* Allow plan to run its cleanup messages (``finalize``) when the RunEngine is
  stopped or aborted.
* When an exception is raised, give each plan in the plan stack an opportunity
  to handle it. If it is handled, carry on.
* The SPEC-style ``tw`` was not passing its parameters through to the
  underlying ``tweak`` plan.
* Silenced un-needed suspenders warnings
* Fix bug in separating devices

Internal Changes
----------------

* Reduce unneeded usage of ``bluesky.plans.single_gen``.
* Don't emit create/save messages with no reads in between.
* Re-work exception handling in main run engine event loop.

v0.5.3 (2016-06-06)
===================

Breaking Changes
----------------

* ``LiveTable`` only displays data from one event stream.
* Remove used global state attribute ``gs.COUNT_TIME``.

Bug Fixes
---------

* Fix "infinite count", ``ct(num=None)``.
* Allow the same data keys to be present in different event streams. But, as
  before, a given data key can only appear once per event.
* Make SPEC-style plan ``ct`` implement baseline readings, referring to
  ``gs.BASELINE_DETS``.
* Upon resuming after a deferred pause, clear the deferred pause request.
* Make ``bluesky.utils.register_transform`` character configurable.

v0.5.2 (2016-05-25)
===================

Features
--------

* Plans were reimplemented as simple Python generators instead of custom Python
  classes. The old "object-oriented" plans are maintained for
  back-compatibility. See plans documentation to review new capabilities.

Breaking Changes
----------------

* SPEC-style plans are now proper generators, not bound to the RunEngine.

v0.5.0 (2016-05-11)
===================

Breaking Changes
----------------

* Move ``bluesky.scientific_callbacks`` to ``bluesky.callbacks.scientific``
  and ``bluesky.broker_callbacks`` to ``bluesky.callbacks.broker``.
* Remove ``bluesky.register_mds`` whose usage can be replaced by:
  ``import metadatastore.commands; RE.subscribe_lossless('all', metadatastore.commands.insert)``
* In all occurrences, the argument ``block_group`` has been renamed ``group``
  for consistency. This affects the 'trigger' and 'set' messages.
* The (not widely used) ``Center`` plan has been removed. It may be
  distributed separately in the future.
* Calling a "SPEC-like" plan now returns a generator that must be passed
  to the RunEngine; it does not execute the plan with the global RunEngine in
  gs.RE. There is a convenience wrapper available to restore the old behavior
  as desired. But since that usage renders the plans un-composable, it is
  discouraged.
* The 'time' argument of the SPEC-like plans is a keyword-only argument.
* The following special-case SPEC-like scans have been removed

    * hscan
    * kscan
    * lscan
    * tscan
    * dtscan
    * hklscan
    * hklmesh

  They can be defined in configuration files as desired, and in that location
  they will be easier to customize.
* The ``describe`` method on flyers, which returns an iterable of dicts of
  data keys for one or more descriptors documents, has been renamed to
  ``describe_collect`` to avoid confusion with ``describe`` on other devices,
  which returns a dict of data keys for one descriptor document.
* An obscure feature in ``RunEngine.request_pause`` has been removed, which
  involved removing the optional arguments ``callback`` and ``name``.

v0.4.3 (2016-03-03)
===================

Bug Fixes
---------

* Address serious performance problem in ``LiveTable``.

v0.4.2 (2016-03-02)
===================

Breaking Changes
----------------

* Stage the ultimate parent ("root") when a device is staging its child, making
  it impossible to leave a device in a partially-staged state.

v0.4.1 (2016-02-29)
===================

Features
--------

* Give every event stream a ``name``, using ``'primary'`` by default.
* Record a mapping of device/signal names to ordered data keys in the
  EventDescriptor.
* Let ``LiveRaster`` account for "snaked" trajectories.

Bug Fixes
---------

* ``PeakStats.com`` is a scalar, not a single-element array.
* Restore Python 3.4 compatibility.

v0.4.0 (2016-02-23)
===================

(TO DO)

v0.3.2 (2015-10-28)
===================

(TO DO)

v0.3.1 (2015-10-15)
===================

(TO DO)

v0.3.0 (2015-10-14)
===================

Breaking Changes
----------------

* Removed ``RunEngine.persistent_fields``; all fields in ``RE.md`` persist
  between runs by default.
* No metadata fields are "reserved"; any can be overwritten by the user.
* No metadata fields are absolutely required. The metadata validation function
  is user-customizable. The default validation function behaves the same
  as previous versions of bluesky, but it is no longer manditory.
* The signature of ``RunEngine`` has changed. The ``logbook`` argument is now
  keyword-only, and there is a new keyword-only argument, ``md_validator``.
  See docstring for details.
* The ``configure`` method on readable objects now takes a single optional
  argument, a dictionary that the object can use to configure itself however
  it sees fit. The ``configure`` method always has a new return value, a tuple
  of dicts describing its old and new states:
  ``old, new = obj.configure(state)``
* Removed method ``increment_scan_id``
* `callbacks.broker.post_run` API and docstring brought into agreement.
  The API is change to expect a callable with signature
  ``foo(doc_name, doc)`` rather than

    - a callable which takes a document (as documented)
    - an object with ``start``, ``descriptor``, ``event`` and ``stop``
      methods (as implemented).

  If classes derived from ``CallbackBase`` are being used this will not
  not have any effect on user code.

v0.2.3 (2015-09-29)
===================

(TO DO)

v0.2.2 (2015-09-24)
===================

(TO DO)

v0.2.1 (2015-09-24)
===================

(TO DO)

v0.2.0 (2015-09-22)
===================

(TO DO)

v0.1.0 (2015-06-25)
===================

Initial release

Utility classes and functions
=============================

.. automodule:: bluesky.utils


Msg
---
.. autosummary::
   :nosignatures:
   :toctree: generated

   Msg


Persistent metadata
-------------------

To maintain a peristent set of meta-data between Python sessions
we include a dictionary duck-type based on `zict.Func`.

.. autosummary::
   :nosignatures:
   :toctree: generated

   PersistentDict
   PersistentDict.directory



Internal exceptions
-------------------

We define a number of `Exception` sub-classes for internal signaling.

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunEngineControlException
   RequestAbort
   RequestStop
   RunEngineInterrupted
   NoReplayAllowed
   IllegalMessageSequence
   FailedPause
   FailedStatus
   InvalidCommand
   PlanHalt
   RampFail


Progress bars
-------------

These are used by the RunEngine to display progress bars and
are the clients of the :obj:`~ophyd.status.MoveStatus.watch` API



.. autosummary::
   :nosignatures:
   :toctree: generated

   ProgressBar
   ProgressBar.update
   ProgressBar.draw
   ProgressBar.clear

   ProgressBarManager


During tasks
------------

These objects encapsulate what the  RunEngine should do on its thread while
waiting for the plan to complete in the background thread

.. autosummary::
   :nosignatures:
   :toctree: generated

   DuringTask
   DuringTask.block

   DefaultDuringTask
The RunEngine run loop
======================

*Note: This is a technical document not optimized for user readability.*

In this document, we start with a simplified version of the bluesky RunEngine.
We add more complexity step by step, with commentary.

The heart of bluesky is the ``RunEngine._run`` co-routine which dispatches the
``Msg`` in the plan to functions that actually carry out the requested task.
The core operation is obscured by the layers of exception handling, state
management, and clean up the RunEngine is responsible for. (Some of this may
be refactored in the near future). This document is only going to discuss the
run loop, not Document generation or hardware clean up.

Minimal RunEngine
-----------------

A minimal (run-able) RunEngine is

.. code:: python

    from time import sleep
    import datetime
    now = datetime.datetime.now
    from bluesky import Msg

    function_map = {'print':
                    lambda msg: print('-- {!s:10.10s} : {: <25.25s} --'.format(now().time(), msg.obj)),
                    'sleep':
                    lambda msg: sleep(msg.args[0])}


    def RE_v0(plan):
        for msg in plan:
            func = function_map[msg.command]
            func(msg)

    welcome_plan = [Msg('print', 'hello'), Msg('sleep', None, 1), Msg('print', 'world!')]

    RE_v0(welcome_plan)

which captures one of the key abstractions of bluesky: A plan is
just an iterable of messages. This abstraction means that the to plan an
experiment you only need to generate a stream of ``Msg`` objects and the
RunEngine will take care of actually executing the code.

Adaptive Plans
--------------

Simply having a stream of commands is not quite enough, you may want to
have the code generating the stream of messages be aware of the return
value of a previous ``Msg`` to decide what to do next. This sort of
thing is supported in python using
`generators <https://docs.python.org/3.5/reference/expressions.html#generator-iterator-methods>`__
which 'suspend' their execution at a ``yield`` statement. When you
iterate over a generator, it runs until the next ``yield``
statement, suspends, and yields the value to the code which is iterating
over it.

Switching to generators requires we change our minimal RE to

.. code:: python

    from bluesky.utils import ensure_generator



    def RE_v1(plan):
        plan = ensure_generator(plan)
        last_result = None

        while True:
            try:
                msg = plan.send(last_result)
            except StopIteration:
                # generators indicate they are done by raising
                # StopIteration
                break
            func = function_map[msg.command]
            last_result = func(msg)


which still works with the ``welcome_plan``

.. code:: python

    RE_v1([Msg('print', 'hello'), Msg('sleep', None, 1), Msg('print', 'world!')])

but we can also do more sophisticated things like

.. code:: python

    function_map['sum'] = lambda msg: sum(msg.args)

    def adding_plan(a, b):
        yield Msg('print', '{} + {} = ??'.format(a, b))
        ret = yield Msg('sum', None, a, b)
        yield Msg('print', '{} + {} = {}'.format(a, b, ret))
        yield Msg('print', 'thanks for adding')

Which gives

.. code:: python

    RE_v1(adding_plan(1, 2))
    RE_v1(adding_plan(5, 2))

This is obviously overkill for simple addition, but enables this like an
adaptive dscan that changes the step size based on the local slope.

Exception Handling
------------------

In addition to ``generator.send`` (which inserts a value into the
generator) you can also use ``generator.throw`` which raises an
exception at the point where the generator is paused. If the generator
handles the exception (via a ``try...except`` block) then generator
runs until the next ``yield`` and ``throw`` returns the yielded
value. If the generator does not handle the exception (or raises a
different exception) then it is (re)raised by ``throw``.

We want to be able to capture any exceptions raised by the ``RE``
and pass those back to the plan.

.. code:: python


    def RE_v2(plan):
        plan = ensure_generator(plan)
        last_result = None
        _exception = None
        while True:
            try:
                if _exception is not None:
                    msg = plan.throw(_exception)
                    _exception = None
                else:
                    msg = plan.send(last_result)

            except StopIteration:
                break
            try:
                func = function_map[msg.command]
                last_result = func(msg)
            except Exception as e:
                _exception = e


We can now write plans that handle exception from the RE, in this case
reporting that the addition failed due to a ``TypeError``

.. code:: python

    def safe_adding_plan(a, b):
        yield Msg('print', '{} + {} = ??'.format(a, b))
        try:
            ret = yield Msg('sum', None, a, b)
        except TypeError:
            yield Msg('print', 'can not add {} + {}!'.format(a, b))
        else:
            yield Msg('print', '{} + {} = {}'.format(a, b, ret))
        finally:
            yield Msg('print', 'thanks for adding')

Compare the behavior of between ``adding_plan`` and ``addingplan`` in cases
where they succeed

.. code:: python

    RE_v2(safe_adding_plan(1, 2))
    RE_v2(adding_plan(1, 2))

and fail

.. code:: python

    RE_v2(safe_adding_plan('a', 2))
    RE_v2(adding_plan('a', 2))

Again, this is overkill for these simple cases, but this mechanism
allows us to write delta scans that always return the motors to their
original position, shut shutters, etc even if the plan fails or is
canceled.

Turn into a callable class
--------------------------

We are going to want to have access to the internal state of the
``_run`` loop very soon. An way to do this, while maintaining
the API we have above is to write a callable class instead of a
function.

.. code:: python

    class RunEngine_v3:
        def _sleep(self, msg):
            sleep(msg.args[0])

        def _print(self, msg):
            print('-- {!s:10.10s} : {: <25.25s} --'.format(now().time(), msg.obj)),

        def _sum(self, msg):
            return sum(msg.args)

        def __init__(self):
            self._command_registry = {
                'print': self._print,
                'sum': self._sum,
                'sleep': self._sleep}

        def __call__(self, plan):
            self._run(plan)

        def _run(self, plan):
            plan = ensure_generator(plan)
            last_result = None
            _exception = None
            while True:
                try:
                    if _exception is not None:
                        msg = plan.throw(_exception)
                        _exception = None
                    else:
                        msg = plan.send(last_result)

                except StopIteration:
                    break
                try:
                    func = self._command_registry[msg.command]
                    last_result = func(msg)
                except Exception as e:
                    _exception = e


    RE_v3 = RunEngine_v3()

In doing this we also pulled the function the commands dispatched to
into the class. While these methods are almost trivial, we will soon
have methods that alter the internal state of the ``RunEngine``.

``asyncio`` integration
-----------------------

So far all of these RE implementations have been synchronous functions,
that is they run straight through the plan. However, at a beamline we
need to be able to support asynchronous functionality and gracefully
interrupt the plan.

To enable this we are using ``asyncio`` from the python standard library
(new in 3.4) to provide the outer event loop. At this point we are
integrating together two event loops: the RE loop which is processing
the plan and the ``asyncio`` event loop which is managing multiple
frames of execution. The event loop may switch between execution frames
when a coroutine is suspended by a ``yield from`` statement. Thus we
change the methods we dispatch to and the main ``_run`` method to
co-routines by adding the ``@asyncio.coroutine`` decorator and calling
the dispatched functions via ``yield from`` rather than with a direct
function call.

We also added a ``msg_hook`` attribute to the ``RunEngine``
which is a super handy debugging tool to see exactly what messages are
being processed by the RunEngine. It can be set to any callable which
takes a single ``Msg`` as input (ex ``print``)

.. code:: python

    import asyncio


    class RunEngine_v4:
        def __init__(self, *, loop=None):
            # map messages to coro
            self._command_registry = {
                'print': self._print,
                'sum': self._sum,
                'sleep': self._sleep}

            # debugging hook
            self.msg_hook = None


            # bind RE to a specific loop
            if loop is None:
                loop = asyncio.get_event_loop()
            self.loop = loop

            # The RunEngine keeps track of a *lot* of state.
            # All flags and caches are defined here with a comment. Good luck.
            self._task = None  # asyncio.Task associated with call to self._run

        def __call__(self, plan):
            self._task = self.loop.create_task(self._run(plan))
            self.loop.run_until_complete(self._task)

            if self._task.done() and not self._task.cancelled():
                exc = self._task.exception()
                if exc is not None:
                    raise exc

        @asyncio.coroutine
        def _run(self, plan):
            plan = ensure_generator(plan)
            last_result = None
            _exception = None
            while True:
                try:
                    yield from asyncio.sleep(0.0001, loop=self.loop)
                    if _exception is not None:
                        msg = plan.throw(_exception)
                        _exception = None
                    else:
                        msg = plan.send(last_result)

                except StopIteration:
                    break

                if self.msg_hook:
                    self.msg_hook(msg)

                try:
                    func = self._command_registry[msg.command]
                    last_result = yield from func(msg)
                except Exception as e:
                    _exception = e

        @asyncio.coroutine
        def _sleep(self, msg):
            yield from asyncio.sleep(msg.args[0])

        @asyncio.coroutine
        def _print(self, msg):
            print('-- {!s:10.10s} : {: <25.25s} --'.format(now().time(), msg.obj)),

        @asyncio.coroutine
        def _sum(self, msg):
            return sum(msg.args)



    RE_v4 = RunEngine_v4()

Pausing, Resuming, and Rewinding
--------------------------------

Adding the ability to pause/resume/rewind a scan adds a fair amount of
complexity as now the ``RunEngine`` must keep track of a stack of plans
rather than a single plan, cache ``Msg`` as they go by and expose enough
API to control the behavior.

.. code:: python

    from collections import deque
    import asyncio

    import datetime
    import functools
    from bluesky.utils import (AsyncInput, FailedPause, InvalidCommand, Msg,
                               ensure_generator)
    from bluesky.run_engine import RunEngineStateMachine, PropertyMachine
    from super_state_machine.errors import TransitionError


    class RunEngine_v5:
        state = PropertyMachine(RunEngineStateMachine)
        _UNCACHEABLE_COMMANDS = ['pause', ]

        def __init__(self, *, loop=None):
            # map messages to coro
            self._command_registry = {
                'print': self._print,
                'sum': self._sum,
                # coros on real RE
                'sleep': self._sleep,
                'checkpoint': self._checkpoint,
                'clear_checkpoint': self._clear_checkpoint,
                'rewindable': self._rewindable,
                'pause': self._pause,
                'input': self._input,
                'null': self._null, }

            # debugging hook
            self.msg_hook = None

            # bind RE to a specific loop
            if loop is None:
                loop = asyncio.get_event_loop()
            self.loop = loop

            # The RunEngine keeps track of a *lot* of state.
            # All flags and caches are defined here with a comment. Good luck.
            self._task = None  # asyncio.Task associated with call to self._run

            self._deferred_pause_requested = False  # pause at next 'checkpoint'
            self._msg_cache = deque()  # history of processed msgs for rewinding
            self._rewindable_flag = True  # if the RE is allowed to replay msgs
            self._plan = None  # the scan plan instance from __call__
            self._plan_stack = deque()  # stack of generators to work off of
            self._response_stack = deque([None])  # resps to send into the plans
            self._interrupted = False  # True if paused, aborted, or failed

        def __call__(self, plan):
            # First thing's first: if we are in the wrong state, raise.
            if not self.state.is_idle:
                raise RuntimeError("The RunEngine is in a %s state" % self.state)

            self._clear_call_cache()

            self._plan = plan
            gen = ensure_generator(plan)

            self._plan_stack.append(gen)
            self._response_stack.append(None)

            self._task = self.loop.create_task(self._run())
            self.loop.run_forever()

            if self._task.done() and not self._task.cancelled():
                exc = self._task.exception()
                if exc is not None:
                    raise exc

        def _clear_call_cache(self):
            self._deferred_pause_requested = False
            self._plan_stack = deque()
            self._msg_cache = deque()
            self._response_stack = deque([None])
            self._exception = None
            self._task = None
            self._plan = None
            self._interrupted = False

        @property
        def rewindable(self):
            return self._rewindable_flag

        @rewindable.setter
        def rewindable(self, v):
            cur_state = self._rewindable_flag
            self._rewindable_flag = bool(v)
            if self.resumable and self._rewindable_flag != cur_state:
                self._reset_checkpoint_state()

        @property
        def resumable(self):
            "i.e., can the plan in progress by rewound"
            return self._msg_cache is not None

        @asyncio.coroutine
        def _run(self):
            pending_cancel_exception = None
            try:
                self.state = 'running'
                while True:
                    try:
                        yield from asyncio.sleep(0.0001, loop=self.loop)
                        # The case where we have a stashed exception
                        if self._exception is not None:
                            # throw the exception at the current plan
                            try:
                                msg = self._plan_stack[-1].throw(
                                    self._exception)
                            except Exception as e:
                                # The current plan did not handle it,
                                # maybe the next plan (if any) would like
                                # to try
                                self._plan_stack.pop()
                                if len(self._plan_stack):
                                    self._exception = e
                                    continue
                                # no plans left and still an unhandled exception
                                # re-raise to exit the infinite loop
                                else:
                                    raise
                            # clear the stashed exception, the top plan
                            # handled it.
                            else:
                                self._exception = None
                        # The normal case of clean operation
                        else:
                            resp = self._response_stack.pop()
                            try:
                                msg = self._plan_stack[-1].send(resp)
                            # We have exhausted the top generator
                            except StopIteration:
                                # pop the dead generator go back to the top
                                self._plan_stack.pop()
                                if len(self._plan_stack):
                                    continue
                                # or reraise to get out of the infinite loop
                                else:
                                    raise
                            # Any other exception that comes out of the plan
                            except Exception as e:
                                # pop the dead plan, stash the exception and
                                # go to the top of the loop
                                self._plan_stack.pop()
                                if len(self._plan_stack):
                                    self._exception = e
                                    continue
                                # or reraise to get out of the infinite loop
                                else:
                                    raise

                        if self.msg_hook:
                            self.msg_hook(msg)

                        # if this message can be cached for rewinding, cache it
                        if (self._msg_cache is not None and
                                self._rewindable_flag and
                                msg.command not in self._UNCACHEABLE_COMMANDS):
                            # We have a checkpoint.
                            self._msg_cache.append(msg)

                        # try to look up the coroutine to execute the command
                        try:
                            coro = self._command_registry[msg.command]
                        # replace KeyError with a local sub-class and go
                        # to top of the loop
                        except KeyError:
                            # TODO make this smarter
                            self._exception = InvalidCommand(msg.command)
                            continue

                        # try to finally run the command the user asked for
                        try:
                            # this is one of two places that 'async'
                            # exceptions (coming in via throw) can be
                            # raised
                            response = yield from coro(msg)
                        # special case `CancelledError` and let the outer
                        # exception block deal with it.
                        except asyncio.CancelledError:
                            raise
                        # any other exception, stash it and go to the top of loop
                        except Exception as e:
                            self._exception = e
                            continue
                        # normal use, if it runs cleanly, stash the response and
                        # go to the top of the loop
                        else:
                            self._response_stack.append(response)
                            continue

                    except KeyboardInterrupt:
                        # This only happens if some external code captures SIGINT
                        # -- overriding the RunEngine -- and then raises instead
                        # of (properly) calling the RunEngine's handler.
                        # See https://github.com/NSLS-II/bluesky/pull/242
                        print("An unknown external library has improperly raised "
                              "KeyboardInterrupt. Intercepting and triggering "
                              "a hard pause instead.")
                        self.loop.call_soon(self.request_pause, False)
                        print(PAUSE_MSG)
                    except asyncio.CancelledError as e:
                        # if we are handling this twice, raise and leave the plans
                        # alone
                        if self._exception is e:
                            raise e
                        # the case where FailedPause, RequestAbort or a coro
                        # raised error is not already stashed in _exception
                        if self._exception is None:
                            self._exception = e
                        pending_cancel_exception = e
            except StopIteration:
                pass
            finally:
                self.loop.stop()
                self.state = 'idle'
            # if the task was cancelled
            if pending_cancel_exception is not None:
                raise pending_cancel_exception
        @asyncio.coroutine
        def _sleep(self, msg):
            yield from asyncio.sleep(msg.args[0])

        @asyncio.coroutine
        def _print(self, msg):
            now = datetime.datetime.now
            print('-- {!s:10.10s} : {: <25.25s} --'.format(now().time(), msg.obj))

        @asyncio.coroutine
        def _sum(self, msg):
            return sum(msg.args)

        @asyncio.coroutine
        def _input(self, msg):
            """
            Process a 'input' Msg. Expected Msg:

                Msg('input', None)
                Msg('input', None, prompt='>')  # customize prompt
            """
            prompt = msg.kwargs.get('prompt', '')
            async_input = AsyncInput(self.loop)
            async_input = functools.partial(async_input, end='', flush=True)
            return (yield from async_input(prompt))

        @asyncio.coroutine
        def _pause(self, msg):
            """Request the run engine to pause

            Expected message object is:

                Msg('pause', defer=False, name=None, callback=None)

            See RunEngine.request_pause() docstring for explanation of the three
            keyword arguments in the `Msg` signature
            """
            self.request_pause(*msg.args, **msg.kwargs)

        def request_pause(self, defer=False):
            """
            Command the Run Engine to pause.

            This function is called by 'pause' Messages. It can also be called
            by other threads. It cannot be called on the main thread during a run,
            but it is called by SIGINT (i.e., Ctrl+C).

            If there current run has no checkpoint (via the 'clear_checkpoint'
            message), this will cause the run to abort.

            Parameters
            ----------
            defer : bool, optional
                If False, pause immediately before processing any new messages.
                If True, pause at the next checkpoint.
                False by default.
            """
            if defer:
                self._deferred_pause_requested = True
                print("Deferred pause acknowledged. Continuing to checkpoint.")
                return

            # We are pausing. Cancel any deferred pause previously requested.
            self._deferred_pause_requested = False
            self._interrupted = True
            print("Pausing...")
            self.state = 'paused'
            if not self.resumable:
                # cannot resume, so we cannot pause.  Abort the scan
                print("No checkpoint; cannot pause.")
                print("Aborting: running cleanup and marking "
                      "exit_status as 'abort'...")
                self._exception = FailedPause()
                self._task.cancel()
                for task in self._failed_status_tasks:
                    task.cancel()
                return
            # stop accepting new tasks in the event loop (existing tasks will
            # still be processed)
            self.loop.stop()

        def resume(self):
            """Resume a paused plan from the last checkpoint.

            Returns
            -------
            uids : list
                list of Header uids (a.k.a RunStart uids) of run(s)
            """
            # The state machine does not capture the whole picture.
            if not self.state.is_paused:
                raise TransitionError("The RunEngine is the {0} state. "
                                      "You can only resume for the paused state."
                                      "".format(self.state))

            self._interrupted = False
            new_plan = self._rewind()
            self._plan_stack.append(new_plan)
            self._response_stack.append(None)

            self._resume_event_loop()
            return []

        def _rewind(self):
            '''Clean up in preparation for resuming from a pause or suspension.

            Returns
            -------
            new_plan : generator
                 A new plan made from the messages in the message cache

            '''
            new_plan = ensure_generator(list(self._msg_cache))
            self._msg_cache = deque()
            # This is needed to 'cancel' an open bundling (e.g. create) if
            # the pause happens after a 'checkpoint', after a 'create', but before
            # the paired 'save'.
            return new_plan

        def _resume_event_loop(self):
            # may be called by 'resume' or 'abort'
            self.state = 'running'
            self._last_sigint_time = None
            self._num_sigints_processed = 0

            if self._task.done():
                return
            self.loop.run_forever()
            if self._task.done() and not self._task.cancelled():
                exc = self._task.exception()
                if exc is not None:
                    raise exc

        @asyncio.coroutine
        def _checkpoint(self, msg):
            """Instruct the RunEngine to create a checkpoint so that we can rewind
            to this point if necessary

            Expected message object is:

                Msg('checkpoint')
            """
            yield from self._reset_checkpoint_state_coro()

            if self._deferred_pause_requested:
                # We are at a checkpoint; we are done deferring the pause.
                # Give the _check_for_signals coroutine time to look for
                # additional SIGINTs that would trigger an abort.
                yield from asyncio.sleep(0.5, loop=self.loop)
                self.request_pause(defer=False)

        def _reset_checkpoint_state(self):
            if self._msg_cache is None:
                return

            self._msg_cache = deque()

        _reset_checkpoint_state_coro = asyncio.coroutine(_reset_checkpoint_state)

        @asyncio.coroutine
        def _clear_checkpoint(self, msg):
            """Clear a set checkpoint

            Expected message object is:

                Msg('clear_checkpoint')
            """
            # clear message cache
            self._msg_cache = None
            # clear stashed
            self._teed_sequence_counters.clear()

        @asyncio.coroutine
        def _rewindable(self, msg):
            '''Set rewindable state of RunEngine

            Expected message object is:

                Msg('rewindable', None, bool or None)
            '''

            rw_flag, = msg.args
            if rw_flag is not None:
                self.rewindable = rw_flag

            return self.rewindable

        @asyncio.coroutine
        def _null(self, msg):
            """
            A no-op message, mainly for debugging and testing.
            """
            pass


    RE_v5 = RunEngine_v5()
    RE_v5.msg_hook = print


    def pausing_plan():
        yield Msg('null')
        yield Msg('null')
        yield Msg('pause')
        yield Msg('null')

Stop, Abort, Halt
-----------------

Suspending
----------

Object/hardware clean up
------------------------

Document creation and emission
------------------------------

SIGINT interception
-------------------
.. currentmodule:: bluesky.plans

====================
 Recording Metadata
====================

Capturing useful metadata is the main objective of bluesky. The more
information you can provide about what you are doing and why you are
doing it, the more useful bluesky and downstream data search and
analysis tools can be.

The term "metadata" can be a controversial term, one scientist's
"data" is another's "metadata" and classification is context- dependent.
The same exact information can be "data" in one
experiment, but "metadata" in a different experiment done on the exact
same hardware.
The `Document Model
<https://blueskyproject.io/event-model/data-model.html>`_ provides a framework
for deciding _where_ to record a particular piece of information.

There are some things that we know *a priori* before doing an experiment;
where are we? who is the user? what sample are we looking at? what did
the user just ask us to do?  These are all things that we can, in
principle, know independent of the control system.  These are the
prime candidates for inclusion in the `Start Document
<https://blueskyproject.io/event-model/data-model.html#run-start-document>`_.
Downstream DataBroker provides tools to do rich searches on this data.
The more information you can include the better.

There is some information that we need that is nominally independent of
any particular device but we need to consult the controls system
about.  For example the location of important, but un-scanned motors
or the configuration of beam attenuators.  If the values *should* be fixed over
the course of the experiment then this it is a good candidate for
being a "baseline device" either via the `Supplemental pre-processor
<https://blueskyproject.io/bluesky/tutorial.html#baseline-readings-and-other-supplemental-data>`_
or explicitly in custom plans.  This will put the readings in a separate stream
(which is a peer to the "primary" data).  In principle, these values *could* be
read from the control system once and put into the Start document along with
the *a priori* information, however that has several draw backs:

1. There is only ever 1 reading of the values so if they do drift during
   data acquisition, you will never know.
2. We cannot automatically capture information about the device like
   we do for data in Events.  This includes things like the datatype,
   units, and shape of the value and any configuration information about the
   hardware it is being read from.

A third class of information that can be called "metadata" is
configuration information of pieces of hardware.  These are things
like the velocity of a motor or the integration time of a detector.
These readings are embedded in the `Descriptor
<https://blueskyproject.io/event-model/data-model.html#event-descriptor>`_
and are extracted from the hardware via the `read_configuration
<https://blueskyproject.io/bluesky/hardware.html#ReadableDevice.read_configuration>`_
method of the hardware.  We expect that these values will not change over
the course of the experiment so only read them once.

Information that does not fall into one of these categories, because
you expect it to change during the experiment,
should be treated as "data", either as an explicit part of the
experimental plan or via :ref:`async_monitoring`.


Adding to the Start Document
============================

When the RunEngine mints a Start document it includes structured data.  That
information can be injected in via several mechanisms:

1. entered interactively by the user at execution time
2. provided in the code of the *plan*
3. automatically inferred
4. entered by user once and stashed for reuse on all future plans

If there is a conflict between these sources, the higher entry in this
list wins.  The "closer" to a user the information originated the
higher priority it has.


1. Interactively, for One Use
-----------------------------

Suppose we are executing some custom plan called ``plan``.

.. code-block:: python

    RE(plan())

If we give arbitrary extra keyword arguments to ``RE``, they will be
interpreted as metadata.

.. code-block:: python

    RE(plan(), sample_id='A', purpose='calibration', operator='Dan')

The :ref:`run(s) <run_overview>` --- i.e., datasets --- generated by ``plan()``
will include the custom metadata:

.. code-block:: python

    ...
    'sample_id': 'A',
    'purpose': 'calibration'.
    'operator': 'Dan',
    ...

If ``plan`` generates more that one run, all the runs will get this metadata.
For example, this plan generates three different runs.

.. code-block:: python

    from bluesky.plans import count, scan
    from ophyd.sim det1, det2, motor  # simulated detectors, motor

    def plan():
        yield from count([det])
        yield from scan([det], motor, 1, 5, 5)
        yield from count([det])

If executed as above:

.. code-block:: python

    RE(plan(), sample_id='A', purpose='calibration', operator='Dan')

each run will get a copy of the sample_id, purpose and operator metadata.

2. Through a plan
-----------------

Revisiting the previous example:

.. code-block:: python

    def plan():
        yield from count([det])
        yield from scan([det], motor, 1, 5, 5)
        yield from count([det])

we can pass different metadata for each run. Every
:ref:`built-in pre-assembled plan <preassembled_plans>` accepts a parameter
``md``, which you can use to inject metadata that applies only to that plan.

.. code-block:: python

    def plan():
        yield from count([det], md={'purpose': 'calibration'})  # one
        yield from scan([det], motor, 1, 5, 5, md={'purpose': 'good data'})  # two
        yield from count([det], md={'purpose': 'sanity check'})  # three

The metadata passed into ``RE`` is combined with the metadata passed in to each
plan. Thus, calling

.. code-block:: python

    RE(plan(), sample_id='A', operator='Dan')

generates these three sets of metadata:

.. code-block:: python

    # one
    ...
    'sample_id': 'A',
    'purpose': 'calibration'.
    'operator': 'Dan',
    ...

    # two
    ...
    'sample_id': 'A',
    'purpose': 'good data'.
    'operator': 'Dan',
    ...

    # three
    ...
    'sample_id': 'A',
    'purpose': 'sanity check'.
    'operator': 'Dan',
    ...

If there is a conflict, ``RE`` keywords takes precedence. So

.. code-block:: python

    RE(plan(), purpose='test')

would override the individual 'purpose' metadata from the plan, marking all
three as purpose=test.

If you define your own plans, it is best practice have them take a keyword only
argument ``md=None``.  This allows the hard-coded meta-data to be over-ridden
later:

.. code-block:: python

    def plan(*, md=None):
        md = md or {}  # handle the default case
        # putting unpacking **md at the end means it "wins"
        # and if the user calls
        #    yield from plan(md={'purpose': bob})
        # it will over-ride these values
        yield from count([det], md={'purpose': 'calibration', **md})
        yield from scan([det], motor, 1, 5, 5, md={'purpose': 'good data', **md})
        yield from count([det], md={'purpose': 'sanity check', **md})

This is consistent with all of the :ref:`preassembled_plans`.

For more on injecting metadata via plans, refer to
:ref:`this section <tutorial_plan_metadata>` of the tutorial.

.. note::

    All of the built-in plans provide certain metadata automatically. Custom
    plans are not *required* to provide any of this, but it is a nice pattern
    to follow.

    * plan_name --- e.g., ``'scan'``
    * detectors --- a list of the names of the detectors
    * motors --- a list of the names of the motors
    * plan_args --- dict of keyword arguments passed to the plan
    * plan_pattern -- function used to create the trajectory
    * plan_pattern_module --- Python module where ``plan_pattern`` is defined
    * plan_pattern_args --- dict of keyword arguments passed to
      ``plan_pattern`` to create the trajectory

    The ``plan_name`` and ``plan_args`` together should provide sufficient
    information to recreate the plan. The ``detectors`` and ``motors`` are
    convenient keys to search on later.

    The ``plan_pattern*`` entries provide lower-level, more explicit
    information about the *trajectory* ("pattern") generated by the plan,
    separate from the specific detectors and motors involved. For complex
    trajectories like spirals, this is especially useful. As a simple example,
    here is the pattern-related metadata for :func:`scan`.

    .. code-block:: python

        ...
        'plan_pattern': 'linspace',
        'plan_pattern_module': 'numpy',
        'plan_pattern_args': dict(start=start, stop=stop, num=num)
        ...

    Thus, one can re-create the "pattern" (trajectory) like so:

    .. code-block:: python

        numpy.linspace(**dict(start=start, stop=stop, num=num))

3. Automatically
----------------

For each run, the RunEngine automatically records:

* 'time' --- In this context, the start time. (Other times are also recorded.)
* 'uid' --- a globally unique ID for this run
* 'plan_name' --- the function or class name of ``plan`` (e.g., 'count')
* 'plan_type'--- e.g., the Python type of ``plan`` (e.g., 'generator')

The last two can be overridden by any of the methods above. The first two
cannot be overridden by the user.

.. note::

    If some custom plan does not specify a 'plan_name' and 'plan_type', the
    RunEngine infers them as follows:

    .. code-block:: python

        plan_name = type(plan).__name__
        plan_type = getattr(plan, '__name__', '')

    These may be more or less informative depending on what ``plan`` is. They
    are just heuristics to provide *some* information by default if the plan
    itself and the user do not provide it.

4. Interactively, for Repeated Use
----------------------------------

Each time a plan is executed, the current contents of ``RE.md`` are copied into
the metadata for all runs generated by the plan.  To enter metadata once to
reuse on all plans, add it to ``RE.md``.

.. code-block:: python

    RE.md['proposal_id'] = 123456
    RE.md['project'] = 'flying cars'
    RE.md['dimensions'] = (5, 3, 10)

View its current contents,

.. code-block:: python

    RE.md

delete a key you want to stop using,

.. code-block:: python

    del RE.md['project']   # delete a key

or use any of the standard methods that apply to
`dictionaries in Python <https://docs.python.org/3/library/stdtypes.html#typesmapping>`_.

.. warning::


   In general we recommend against putting device readings in the Start
   document. (The Start document is for who/what/why/when, things you
   know before you start communicating with hardware.) It is *especially*
   critical that you do not put device readings in the ``RE.md`` dictionary.
   The value will remain until you change it and not track the state of the
   hardware.  This will result in recording out-of-date, incorrect data!

   This can be particularly dangerous if ``RE.md`` is backed by a
   persistent data store (see next section) because out-of-date readings will
   last across sessions.


The ``scan_id``, an integer that the RunEngine automatically increments at the
beginning of each scan, is stored in ``RE.md['scan_id']``.

.. warning::

    Clearing all keys, like so:

    .. code-block:: python

        RE.md.clear()  # clear *all* keys

    will reset the ``scan_id``. The next time a plan is executed, the
    RunEngine will start with a ``scan_id`` of 1 and set

    .. code-block:: python

        RE.md['scan_id'] = 1

    Some readers may prefer to reset the scan ID to 1 at the beginning of a new
    experiment; others way wish to maintain a single unbroken sequence of scan
    IDs forever.

    From a technical standpoint, it is fine to have duplicate scan IDs. All
    runs also have randomly-generated 'uid' ("unique ID") which is globally
    unique forever.

.. _md_persistence:

Persistence Between Sessions
----------------------------

We provide a way to save the contents of the metadata stash ``RE.md`` between
sessions (e.g., exiting and re-opening IPython).

In general, the ``RE.md`` attribute may be anything that supports the
dictionary interface. The simplest is just a plain Python dictionary.

.. code-block:: python

    RE.md = {}

To persist metadata between sessions, bluesky recommends
:class:`bluesky.utils.PersistentDict` --- a Python dictionary synced with a
directory of files on disk. Any changes made to ``RE.md`` are synced to the
file, so the contents of ``RE.md`` can persist between sessions.

.. code-block:: python

    from bluesky.utils import PersistentDict
    RE.md = PersistentDict('some/path/here')

Bluesky does not provide a strong recommendation on that path; that a detail
left to the local deployment.

Bluesky formerly recommended using :class:`~historydict.HistoryDict` --- a
Python dictionary backed by a sqlite database file. This approach proved
problematic with the threading introduced in bluesky v1.6.0, so it is no longer
recommended. If you have been following that recommendation, you should migrate
your metadata from `~historydict.HistoryDict` to
:class:`~bluesky.utils.PersistentDict`. First, update your configuration to
make ``RE.md`` a :class:`~bluesky.utils.PersistentDict` as shown above. Then,
migrate like so:

.. code-block:: python

   from bluesky.utils import get_history
   old_md = get_history()
   RE.md.update(old_md)

The :class:`~bluesky.utils.PersistentDict` object has been back-ported to
bluesky v1.5.6 as well. It is not available in 1.4.x or older, so once you move
to the new system, you must run bluesky v1.5.6 or higher.

.. warning::

    The ``RE.md`` object can also be set when the RunEngine is instantiated:

    .. code-block:: python

        # This:
        RE = RunEngine(...)

        # is equivalent to this:
        RE = RunEngine({})
        RE.md = ...

    As we stated
    :ref:`at the start of the tutorial <tutorial_run_engine_setup>`, if you are
    using bluesky at a user facility or with shared configuration, your
    ``RE`` may already be configured, and defining a new ``RE`` as above can
    result in data loss! If you aren't sure, it's safer to use ``RE.md = ...``.


Allowed Data Types
------------------

Custom metadata keywords can be mapped to:

* strings --- e.g., ``task='calibration'``
* numbers --- e.g., ``attempt=5``
* lists or tuples --- e.g., ``dimensions=[1, 3]``
* (nested) dictionaries --- e.g., ``dimensions={'width': 1, 'height': 3}``


Required Fields
---------------

The fields:

* **uid**
* **time**

are reserved by the document model and cannot be set by the user.

In current versions of bluesky, **no fields are universally required by bluesky
itself**. It is possible specify your own required fields in local
configuration. See :ref:`md_validator`. (At NSLS-II, there are facility-wide
requirements coming soon.)


Special Fields
--------------

Arbitrary custom fields are allowed --- you can invent any names that are
useful to you.

But certain fields are given special significance by bluesky's document model,
and are either disallowed are required to be a certain type.

The fields:

* **owner**
* **group**
* **project**

are optional but, to facilitate searchability, if they are not blank they must
be strings. A non-string, like ``owner=5`` will produce an error that will
interrupt scan execution immediately after it starts.

Similarly, the keyword **sample** has special significance. It must be either a
string or a dictionary.

The **scan_id** field is expected to be an integer, and it is automatically
incremented between runs. If a scan_id is not provided by the user or stashed
in the persistent metadata from the previous run, it defaults to 1.


.. _md_validator:

Validation
----------

Additional, customized metadata validation can be added to the RunEngine.
For example, to ensure that a run will not be executed unless the parameter
'sample_number' is specified, define a function that accepts a dictionary
argument and raises if 'sample_number' is not found.

.. code-block:: python

    def ensure_sample_number(md):
        if 'sample_number' not in md:
            raise ValueError("You forgot the sample number.")

Apply this function by setting

.. code-block:: python

    RE.md_validator = ensure_sample_number

The function will be executed immediately before each new run in opened.
Appendix
========

This section covers Python language features that may be new to some readers.
They are used by bluesky but not *unique* Wherever possible, we to bluesky.

.. _yield_from_primer:

A Primer on ``yield`` and ``yield from``
----------------------------------------

This is a very brief primer on the Python syntax ``yield`` and ``yield from``,
a feature of the core language that we will use extensively.

A Python *function* returns once:

.. ipython:: python

    def f():
        return 1

    f()

A Python *generator* is like a function with multiple exit points. Calling a
generator produces an *iterator* that yields one value at a time. After
each ``yield`` statement, its execution is suspended.

.. ipython:: python

    def f():
        yield 1
        yield 2

We can exhaust the generator (i.e., get all its values) by calling ``list()``.

.. ipython:: python

    list(f())

We can get one value at a time by calling ``next()``

.. ipython:: python

    it = f()
    next(it)
    next(it)

or by looping through the values.

.. ipython:: python

    for val in f():
        print(val)

To examine what is happening when, we can add prints.

.. ipython:: python

    def verbose_f():
        print("before 1")
        yield 1
        print("before 2")
        yield 2

.. ipython:: python

    it = verbose_f()
    next(it)
    next(it)

Notice that execution is suspended after the first yield statement. The
second ``print`` is not run until we resume execution by requesting a second
value. This is a useful feature of generators: they can express "lazy"
execution.

Generators can delegate to other generators using ``yield from``. This is
syntax we commonly use to combine plans.

.. ipython:: python

    def double_f():
        yield from f()
        yield from f()

The above is equivalent to:

.. ipython:: python

    def double_f():
        for val in f():
            yield val
        for val in f():
            yield val

The ``yield from`` syntax is just more succinct.

.. ipython:: python

    list(double_f())
.. currentmodule:: bluesky.plans

=====
Plans
=====

A *plan* is bluesky's concept of an experimental procedure. A plan may be any
iterable object (list, tuple, custom iterable class, ...) but most commonly it
is implemented as a Python generator. For a more technical discussion we refer
you :doc:`msg`.

A variety of pre-assembled plans are provided. Like sandwiches on a deli menu,
you can use our pre-assembled plans or assemble your own from the same
ingredients, catalogued under the heading :ref:`stub_plans` below.

.. note::

    In the examples that follow, we will assume that you have a RunEngine
    instance named ``RE``. This may have already been configured for you if you
    are a user at a facility that runs bluesky. See
    :ref:`this section of the tutorial <tutorial_run_engine_setup>` to sort out
    if you already have a RunEngine and to quickly make one if needed.

.. _preassembled_plans:

Pre-assembled Plans
===================

Below this summary table, we break the down the plans by category and show
examples with figures.

Summary
-------

Notice that the names in the left column are links to detailed API
documentation.

.. autosummary::
   :toctree: generated
   :nosignatures:

   count
   scan
   rel_scan
   list_scan
   rel_list_scan
   list_grid_scan
   rel_list_grid_scan
   log_scan
   rel_log_scan
   grid_scan
   rel_grid_scan
   scan_nd
   spiral
   spiral_fermat
   spiral_square
   rel_spiral
   rel_spiral_fermat
   rel_spiral_square
   adaptive_scan
   rel_adaptive_scan
   tune_centroid
   tweak
   ramp_plan
   fly


Time series ("count")
---------------------

Examples:

.. code-block:: python

    from ophyd.sim import det
    from bluesky.plans import count

    # a single reading of the detector 'det'
    RE(count([det]))

    # five consecutive readings
    RE(count([det], num=5))

    # five sequential readings separated by a 1-second delay
    RE(count([det], num=5, delay=1))

    # a variable delay
    RE(count([det], num=5, delay=[1, 2, 3, 4]))

    # Take readings forever, until interrupted (e.g., with Ctrl+C)
    RE(count([det], num=None))

.. code-block:: python

    # We'll use the 'noisy_det' example detector for a more interesting plot.
    from ophyd.sim import noisy_det

    RE(count([noisy_det], num=5))


.. plot::

    from bluesky import RunEngine
    from bluesky.plans import count
    from ophyd.sim import noisy_det
    from bluesky.callbacks.best_effort import BestEffortCallback
    bec = BestEffortCallback()
    RE = RunEngine({})
    RE.subscribe(bec)
    RE(count([noisy_det], num=5))

.. note::

   Why doesn't :func:`count` have an ``exposure_time`` parameter?

   Modern CCD detectors typically parametrize exposure time with *multiple*
   parameters (acquire time, acquire period, num exposures, ...) as do scalers
   (preset time, auto count time). There is no one "exposure time" that can be
   applied to all detectors.

   Additionally, when using multiple detectors as in ``count([det1, det2]))``,
   the user would need to provide a separate exposure time for each detector in
   the general case, which would grow wordy.

   One option is to set the time-related parameter(s) as a separate step.

   For interactive use:

   .. code-block:: python

      # Just an example. Your detector might have different names or numbers of
      # exposure-related parameters---which is the point.
      det.exposure_time.set(3)
      det.acquire_period.set(3.5)

   From a plan:

   .. code-block:: python

      # Just an example. Your detector might have different names or numbers of
      # exposure-related parameters---which is the point.
       yield from bluesky.plan_stubs.mv(
           det.exposure_time, 3,
           det.acquire_period, 3.5)

   Another is to write a custom plan that wraps :func:`count` and sets the
   exposure time. This plan can encode the details that bluesky in general
   can't know.

   .. code-block:: python

      def count_with_time(detectors, num, delay, exposure_time, *, md=None):
          # Assume all detectors have one exposure time component called
          # 'exposure_time' that fully specifies its exposure.
          for detector in detectors:
              yield from bluesky.plan_stubs.mv(detector.exposure_time, exposure_time)
          yield from bluesky.plans.count(detectors, num, delay, md=md)

.. autosummary::
   :toctree: generated
   :nosignatures:

   count

Scans over one dimension
------------------------

The "dimension" might be a physical motor position, a temperature, or a
pseudo-axis. It's all the same to the plans. Examples:

.. code-block:: python

    from ophyd.sim import det, motor
    from bluesky.plans import scan, rel_scan, list_scan

    # scan a motor from 1 to 5, taking 5 equally-spaced readings of 'det'
    RE(scan([det], motor, 1, 5, 5))

    # scan a motor from 1 to 5 *relative to its current position*
    RE(rel_scan([det], motor, 1, 5, 5))

    # scan a motor through a list of user-specified positions
    RE(list_scan([det], motor, [1, 1, 2, 3, 5, 8]))

.. code-block:: python

    RE(scan([det], motor, 1, 5, 5))

.. plot::

    from bluesky import RunEngine
    from bluesky.plans import scan
    from ophyd.sim import det, motor
    RE = RunEngine({})
    from bluesky.callbacks.best_effort import BestEffortCallback
    bec = BestEffortCallback()
    RE.subscribe(bec)
    RE(scan([det], motor, 1, 5, 5))

.. note::

   Why don't scans have a ``delay`` parameter?

   You may have noticed that :func:`count` has a ``delay`` parameter but none
   of the scans do. This is intentional.

   The common reason for wanting a delay in a scan is to allow a motor to
   settle or a temperature controller to reach equilibrium. It is better to
   configure this on the respective devices, so that scans will always add the
   appropriate delay for the particular device being scanned.

   .. code-block:: python

      motor.settle_time = 1
      temperature_controller.settle_time = 10

   For many cases, this is more convenient and more robust than typing a delay
   parameter in every invocation of the scan. You only have to set it once, and
   it applies thereafter.

   This is why bluesky leaves ``delay`` out of the scans, to guide users toward
   an approach that will likely be a better fit than the one that might occur
   to them first. For situations where a ``delay`` parameter really is the
   right tool for the job, it is of course always possible to add a ``delay``
   parameter yourself by writing a custom plan. Here is one approach, using a
   :ref:`per_step hook <per_step_hook>`.

   .. code-block:: python

      import bluesky.plans
      import bluesky.plan_stubs

      def scan_with_delay(*args, delay=0, **kwargs):
          "Accepts all the normal 'scan' parameters, plus an optional delay."

          def one_nd_step_with_delay(detectors, step, pos_cache):
              "This is a copy of bluesky.plan_stubs.one_nd_step with a sleep added."
              motors = step.keys()
              yield from bluesky.plan_stubs.move_per_step(step, pos_cache)
              yield from bluesky.plan_stubs.sleep(delay)
              yield from bluesky.plan_stubs.trigger_and_read(list(detectors) + list(motors))

          kwargs.setdefault('per_step', one_nd_step_with_delay)
          yield from bluesky.plans.scan(*args, **kwargs)

.. autosummary::
   :toctree: generated
   :nosignatures:

   scan
   rel_scan
   list_scan
   rel_list_scan
   log_scan
   rel_log_scan

.. _multi-dimensional_scans:

Multi-dimensional scans
-----------------------

See :ref:`tutorial_multiple_motors` in the tutorial for an introduction to the
common cases of moving multiple motors in coordination (i.e. moving X and Y
along a diagonal) or in a grid. The key examples are reproduced here. Again,
see the section linked for further explanation.

.. code-block:: python

    from ophyd.sim import det, motor1, motor2, motor3
    from bluesky.plans import scan, grid_scan, list_scan, list_grid_scan

    RE(scan(dets,
            motor1, -1.5, 1.5,  # scan motor1 from -1.5 to 1.5
            motor2, -0.1, 0.1,  # ...while scanning motor2 from -0.1 to 0.1
            11))  # ...both in 11 steps

    # Scan motor1 and motor2 jointly through a 5-point trajectory.
    RE(list_scan(dets, motor1, [1, 1, 3, 5, 8], motor2, [25, 16, 9, 4, 1]))

    # Scan a 3 x 5 x 2 grid.
    RE(grid_scan([det],
                 motor1, -1.5, 1.5, 3,  # no snake parameter for first motor
                 motor2, -0.1, 0.1, 5, False))
                 motor3, -200, 200, 5, False))

    # Scan a grid with abitrary spacings given as specific positions.
    RE(list_grid_scan([det],
                      motor1, [1, 1, 2, 3, 5],
                      motor2, [25, 16, 9]))

All of these plans are built on a more general-purpose plan,
:func:`~bluesky.plan.scan_nd`, which we can use for more specialized cases.

Some jargon: we speak of :func:`~bluesky.plans.scan`-like joint movement as an
"inner product" of trajectories and :func:`~bluesky.plans.grid_scan`-like
movement as an "outer product" of trajectories. The general case, moving some
motors together in an "inner product" against another motor (or motors) in an
"outer product," can be addressed using a ``cycler``.  Notice what happens when
we add or multiply ``cycler`` objects.

.. ipython:: python

    from cycler import cycler
    from ophyd.sim import motor1, motor2, motor3

    traj1 = cycler(motor1, [1, 2, 3])
    traj2 = cycler(motor2, [10, 20, 30])
    list(traj1)  # a trajectory for motor1
    list(traj1 + traj2)  # an "inner product" trajectory
    list(traj1 * traj2)  # an "outer product" trajectory

We have reproduced inner product and outer product. The real power comes in
when we combine them, like so. Here, motor1 and motor2 together in a mesh
against motor3.

.. ipython:: python

    traj3 = cycler(motor3, [100, 200, 300])
    list((traj1 + traj2) * traj3)

For more on cycler, we refer you to the
`cycler documentation <http://matplotlib.org/cycler/>`_. To build a plan
incorporating these trajectories, use our general N-dimensional scan plan,
:func:`scan_nd`.

.. code-block:: python

    RE(scan_nd([det], (traj1 + traj2) * traj3))

.. autosummary::
   :toctree: generated
   :nosignatures:

   scan
   rel_scan
   grid_scan
   rel_grid_scan
   list_scan
   rel_list_scan
   list_grid_scan
   rel_list_grid_scan
   scan_nd

Spiral trajectories
-------------------

We provide two-dimensional scans that trace out spiral trajectories.

A simple spiral:

.. plot::
   :include-source:

    from bluesky.simulators import plot_raster_path
    from ophyd.sim import motor1, motor2, det
    from bluesky.plans import spiral

    plan = spiral([det], motor1, motor2, x_start=0.0, y_start=0.0, x_range=1.,
                  y_range=1.0, dr=0.1, nth=10)
    plot_raster_path(plan, 'motor1', 'motor2', probe_size=.01)


A fermat spiral:

.. plot::
   :include-source:

    from bluesky.simulators import plot_raster_path
    from ophyd.sim import motor1, motor2, det
    from bluesky.plans import spiral_fermat

    plan = spiral_fermat([det], motor1, motor2, x_start=0.0, y_start=0.0,
                         x_range=2.0, y_range=2.0, dr=0.1, factor=2.0, tilt=0.0)
    plot_raster_path(plan, 'motor1', 'motor2', probe_size=.01, lw=0.1)


A square spiral:

.. plot::
   :include-source:

    from bluesky.simulators import plot_raster_path
    from ophyd.sim import motor1, motor2, det
    from bluesky.plans import spiral_square

    plan = spiral_square([det], motor1, motor2, x_center=0.0, y_center=0.0,
                         x_range=1.0, y_range=1.0, x_num=11, y_num=11)
    plot_raster_path(plan, 'motor1', 'motor2', probe_size=.01)


.. autosummary::
   :toctree: generated
   :nosignatures:

   spiral
   spiral_fermat
   spiral_square
   rel_spiral
   rel_spiral_fermat
   rel_spiral_square

Adaptive scans
--------------

These are one-dimension scans with an adaptive step size tuned to move quickly
over flat regions can concentrate readings in areas of high variation by
computing the local slope aiming for a target delta y between consecutive
points.

This is a basic example of the power of adaptive plan logic.

.. code-block:: python

    from bluesky.plans import adaptive_scan
    from ophyd.sim import motor, det

    RE(adaptive_scan([det], 'det', motor,
                     start=-15,
                     stop=10,
                     min_step=0.01,
                     max_step=5,
                     target_delta=.05,
                     backstep=True))

.. plot::

    from bluesky import RunEngine
    from bluesky.plans import adaptive_scan
    from bluesky.callbacks.best_effort import BestEffortCallback
    bec = BestEffortCallback()
    from ophyd.sim import motor, det

    RE = RunEngine({})
    RE.subscribe(bec)

    RE(adaptive_scan([det], 'det', motor,
                     start=-15.5,
                     stop=10,
                     min_step=0.01,
                     max_step=5,
                     target_delta=.05,
                     backstep=True))

From left to right, the scan lengthens its stride through the flat region. At
first, it steps past the peak. The large jump causes it to double back and then
sample more densely through the peak. As the peak flattens, it lengthens its
stride again.

.. autosummary::
   :toctree: generated
   :nosignatures:

   adaptive_scan
   rel_adaptive_scan

Misc.
-----

.. autosummary::
   :toctree: generated
   :nosignatures:

   tweak
   fly

.. _stub_plans:

Stub Plans
==========
.. currentmodule:: bluesky.plan_stubs

These are the aforementioned "ingredients" for remixing, the pieces from which
the pre-assembled plans above were made. See :ref:`tutorial_custom_plans` in
the tutorial for a practical introduction to these components.

Plans for interacting with hardware:

.. autosummary::
   :nosignatures:
   :toctree: generated

    abs_set
    rel_set
    mv
    mvr
    trigger
    read
    rd
    stage
    unstage
    configure
    stop

Plans for asynchronous acquisition:

.. autosummary::
   :nosignatures:
   :toctree: generated

    monitor
    unmonitor
    kickoff
    complete
    collect

Plans that control the RunEngine:

.. autosummary::
   :nosignatures:
   :toctree: generated

    open_run
    close_run
    create
    save
    drop
    pause
    deferred_pause
    checkpoint
    clear_checkpoint
    sleep
    input_plan
    subscribe
    unsubscribe
    install_suspender
    remove_suspender
    wait
    wait_for
    null

Combinations of the above that are often convenient:

.. autosummary::
   :toctree: generated

    trigger_and_read
    one_1d_step
    one_nd_step
    one_shot
    move_per_step

Special utilities:

.. autosummary::
   :toctree: generated

   repeat
   repeater
   caching_repeater
   broadcast_msg

.. _preprocessors:

Plan Preprocessors
==================
.. currentmodule:: bluesky.preprocessors

.. _supplemental_data:

Supplemental Data
-----------------

Plan preprocessors modify a plans contents on the fly. One common use of a
preprocessor is to take "baseline" readings of a group of devices at the
beginning and end of each run. It is convenient to apply this to *all* plans
executed by a RunEngine using the :class:`SupplementalData`.

.. autoclass:: SupplementalData
    :members:

We have installed a "preprocessor" on the RunEngine. A preprocessor modifies
plans, supplementing or altering their instructions in some way. From now on,
every time we type ``RE(some_plan())``, the RunEngine will silently change
``some_plan()`` to ``sd(some_plan())``, where ``sd`` may insert some extra
instructions. Envision the instructions flow from ``some_plan`` to ``sd`` and
finally to ``RE``. The ``sd`` preprocessors has the opportunity to inspect
he
instructions as they go by and modify them as it sees fit before they get
processed by the RunEngine.

Preprocessor Wrappers and Decorators
------------------------------------

Preprocessors can make arbirary modifcations to a plan, and can get quite
devious. For example, the :func:`relative_set_wrapper` rewrites all positions
to be relative to the initial position.

.. code-block:: python

    def rel_scan(detectors, motor, start, stop, num):
        absolute = scan(detectors, motor, start, stop, num)
        relative = relative_set_wrapper(absolute, [motor])
        yield from relative

This is a subtle but remarkably powerful feature.

Wrappers like :func:`relative_set_wrapper` operate on a generator *instance*,
like ``scan(...)``. There are corresponding decorator functions like
``relative_set_decorator`` that operate on a generator
*function* itself, like :func:`scan`.

.. code-block:: python

    # Using a decorator to modify a generator function
    def rel_scan(detectors, motor, start, stop, num):

        @relative_set_decorator([motor])  # unfamiliar syntax? -- see box below
        def inner_relative_scan():
            yield from scan(detectors, motor, start, stop, num)

        yield from inner_relative_scan()

Incidentally, the name ``inner_relative_scan`` is just an internal variable,
so why did we choose such a verbose name? Why not just name it ``f``? That
would work, of course, but using a descriptive name can make debugging easier.
When navigating gnarly, deeply nested tracebacks, it helps if internal variables
have clear names.

.. note::

    The decorator syntax --- the ``@`` --- is a succinct way of passing a
    function to another function.

    This:

    .. code-block:: python

        @g
        def f(...):
            pass

        f(...)

    is equivalent to

    .. code-block:: python

        g(f)(...)

Built-in Preprocessors
----------------------
.. currentmodule:: bluesky.preprocessors

Each of the following functions named ``<something>_wrapper`` operates on
a generator instance. The corresponding functions named
``<something_decorator>`` operate on a generator function.

.. autosummary::
   :nosignatures:
   :toctree: generated

    baseline_decorator
    baseline_wrapper
    contingency_wrapper
    finalize_decorator
    finalize_wrapper
    fly_during_decorator
    fly_during_wrapper
    inject_md_decorator
    inject_md_wrapper
    lazily_stage_decorator
    lazily_stage_wrapper
    monitor_during_decorator
    monitor_during_wrapper
    relative_set_decorator
    relative_set_wrapper
    reset_positions_decorator
    reset_positions_wrapper
    run_decorator
    run_wrapper
    stage_decorator
    stage_wrapper
    subs_decorator
    subs_wrapper
    suspend_decorator
    suspend_wrapper

Custom Preprocessors
--------------------

The preprocessors are implemented using :func:`msg_mutator` (for altering
messages in place) and :func:`plan_mutator` (for inserting
messages into the plan or removing messages).

It's easiest to learn this by example, studying the implementations of the built-in
processors (catalogued above) in the
`the source of the plans module <https://github.com/NSLS-II/bluesky/blob/master/bluesky/plans.py>`_.

.. _per_step_hook:

Customize Step Scans with ``per_step``
======================================

The one-dimensional and multi-dimensional plans are composed (1) setup,
(2) a loop over a plan to perform at each position, (3) cleanup.

We provide a hook for customizing step (2). This enables you to write a
variation of an existing plan without starting from scratch.

For one-dimensional plans, the default inner loop is:

.. code-block:: python

    from bluesky.plan_stubs import checkpoint, abs_set, trigger_and_read

    def one_1d_step(detectors, motor, step):
        """
        Inner loop of a 1D step scan

        This is the default function for ``per_step`` param in 1D plans.
        """
        yield from checkpoint()
        yield from abs_set(motor, step, wait=True)
        return (yield from trigger_and_read(list(detectors) + [motor]))

Some user-defined function, ``custom_step``, with the same signature can be
used in its place:

.. code-block:: python

    scan([det], motor, 1, 5, 5, per_step=custom_step)

For convenience, this could be wrapped into the definition of a new plan:

.. code-block:: python

    def custom_scan(detectors, motor, start, stop, step, *, md=None):
        yield from scan([det], motor, start, stop, step, md=md
                        per_step=custom_step)

For multi-dimensional plans, the default inner loop is:

.. code-block:: python

    from bluesky.utils import short_uid
    from bluesky.plan_stubs import checkpoint, abs_set, wait, trigger_and_read

    def one_nd_step(detectors, step, pos_cache):
        """
        Inner loop of an N-dimensional step scan

        This is the default function for ``per_step`` param in ND plans.

        Parameters
        ----------
        detectors : iterable
            devices to read
        step : dict
            mapping motors to positions in this step
        pos_cache : dict
            mapping motors to their last-set positions
        """
        def move():
            yield from checkpoint()
            grp = short_uid('set')
            for motor, pos in step.items():
                if pos == pos_cache[motor]:
                    # This step does not move this motor.
                    continue
                yield from abs_set(motor, pos, group=grp)
                pos_cache[motor] = pos
            yield from wait(group=grp)

        motors = step.keys()
        yield from move()
        yield from trigger_and_read(list(detectors) + list(motors))

Likewise, a custom function with the same signature may be passed into the
``per_step`` argument of any of the multi-dimensional plans.

Asynchronous Plans: "Fly Scans" and "Monitoring"
================================================

See the section on :doc:`async` for some context on these terms and, near the
end of the section, some example plans.

.. _plan_utils:

Plan Utilities
==============

These are useful utilities for defining custom plans and plan preprocessors.

.. autosummary::
   :toctree: generated
   :nosignatures:

    pchain
    msg_mutator
    plan_mutator
    single_gen
    make_decorator
.. currentmodule:: bluesky.run_engine

RunEngine API Documentation
===========================

The ``RunEngine``
-----------------

.. autosummary::
   :nosignatures:
   :toctree: generated

    RunEngine

The main user entry point tho the RunEngine is ``RE(my_plan(...))``

.. autosummary::
   :nosignatures:
   :toctree: generated

    RunEngine.__call__

The RunEngine maintains a callback registry of functions that receive any
:doc:`documents` generated by plan execution. These methods add and remove
functions from that registry.


.. autosummary::
   :nosignatures:
   :toctree: generated

   RunEngine.subscribe
   RunEngine.unsubscribe

When the RunEngine is in a paused state, it can be resumed or stopped in
:ref:`various ways <interactive_pause_summary>` using these methods:


.. autosummary::
   :nosignatures:
   :toctree: generated

   RunEngine.resume
   RunEngine.abort
   RunEngine.stop
   RunEngine.halt

The RunEngine can suspend and resume plan execution in response to external
changes. See :ref:`suspenders`.


.. autosummary::
   :nosignatures:
   :toctree: generated

   RunEngine.install_suspender
   RunEngine.remove_suspender
   RunEngine.clear_suspenders

These methods are used internally to pause or suspend the RunEngine.
Typically the user accomplishes this with Ctrl+C or by installing
suspenders, respectively. For special applicaitons, they can be called
directly.


.. autosummary::
   :nosignatures:
   :toctree: generated

   RunEngine.request_pause
   RunEngine.request_suspend

These methods may be used to register custom commands to supplement or
replace the built-in commands recognized by the RunEngine.


.. autosummary::
   :nosignatures:
   :toctree: generated

   RunEngine.register_command
   RunEngine.unregister_command

These methods may be used to list the commands available, or print
a summary as to what they do.

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunEngine.print_command_registry
   RunEngine.commands

The ``RunEngineResult``
-----------------------

A :class:`RunEngineResult` will be returned
if the RunEngine was created with ``call_returns_result=True``. The default
is currently ``False``, so this behavior is opt-in for now, but will change
to opt-out in the future.

The :class:`RunEngineResult` class encapsulates useful information about
the plan that was run, including the ultimate plan return value, run uids,
exit status, and exceptions.

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunEngineResult
   


The ``Dispatcher``
------------------

A RunEngine encapsulates a :class:`Dispatcher` for emitting any
:doc:`documents` generated by plan execution. The methods
:meth:`RunEngine.subscribe` and :meth:`RunEngine.unsubscribe`, documented
above, are aliases to the corresponding methods on the RunEngine's Dispatcher.

.. autosummary::
   :nosignatures:
   :toctree: generated

   Dispatcher
   Dispatcher.subscribe
   Dispatcher.unsubscribe
   Dispatcher.unsubscribe_all
   Dispatcher.process


The ``RunBundler``
------------------


The RunEngine also creates `RunBundler`\s instances to encapsulate the
logic and book keeping for generating events (and allow multiple runs
to be open at once).  In general you should not be directly working with
a `RunBundler`.

.. currentmodule:: bluesky.bundlers

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunBundler

The co-routines for opening and closing a run.

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunBundler.open_run
   RunBundler.close_run

The co-routines for opening / filling / closing / dropping an Event

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunBundler.create
   RunBundler.read
   RunBundler.drop
   RunBundler.save

The co-routines for managing flyers

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunBundler.kickoff
   RunBundler.complete
   RunBundler.collect
   RunBundler.backstop_collect

The co-routines for changing a device configuration

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunBundler.configure

The co-routines for managing monitors

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunBundler.monitor
   RunBundler.suspend_monitors
   RunBundler.restore_monitors
   RunBundler.clear_monitors
   RunBundler.unmonitor

The co-routines for checkpoint management

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunBundler.clear_checkpoint
   RunBundler.reset_checkpoint_state
   RunBundler.reset_checkpoint_state_coro
   RunBundler.rewind

The co-routines to record interruptions

.. autosummary::
   :nosignatures:
   :toctree: generated

   RunBundler.record_interruption
Multi-Run Plans
===============

Introduction
------------

This section is a brief tutorial on multi-run plans (introduced in Bluesky v1.6.0).
A traditional single-run plan contains a set of instructions for performing only one run,
which is assigned a scan ID and a UID. When a multi-run plan is executed by the Run Engine, multiple
runs can be performed as part of a single plan. Data from each run can be independently
displayed and saved to the database via Databroker. Prior versions of Bluesky supported
only sequential execution of multiple runs within a plan: building larger plans by creating
a sequence of smaller plans and preassembled plans shipped with Bluesky is a standard
practice. In Bluesky v1.6.0 a number of features were introduced to allow plans
with nested runs. Two runs are considered nested if one 'outer' run is interrupted, another
'inner' run is executed, and then the first run is resumed and completed. The number of levels
of nesting is not limited by Bluesky. Interruptions can be initiated by the plan itself
(simply by opening another run before closing currently executed run) or externally (e.g.
by triggering a suspender and causing execution of pre- or post-plan). This tutorial includes
a brief explanation of the new Bluesky features for supporting multi-run plans and several
examples that demonstrate the implementation of plans that contain sequential, nested and recursive
runs.

Definition of a 'Run'
---------------------

From the point of view of Bluesky, a run is a sequence of instructions (messages) for controlling
the instrumental equipment that starts with `open_run` and ends with `close_run` message.
We may also apply the term 'run' to a block of code which generates such a sequence of messages.
Data from each run is bundled together via an assigned distinct Scan ID and UID. The set of documents
is also generated for each run, including mandatory 'start' and 'stop' documents. The documents
can be processed by callbacks (such as BestEffortCallback) and saved to the database via Databroker.

In the plan, the run may be defined by explicitely enclosing the code in `bps.open_run()` and
`bps.close_run()` stubs:

.. code-block:: python

    # Using 'bps.open_run()' and 'bps.close_run()' stubs to define a run

    import bluesky.plan_stubs as bps
    from bluesky import RunEngine

    RE = RunEngine({})

    def sample_plan():
        ...
        yield from bps.open_run(md={})  # 'md' - metadata to be added to the 'start' document
        ...
        < code that controls execution of the scan >
        ...
        yield from bps.close_run()

    RE(sample_plan())

or using `@bpp.run_decorator`, which inserts `open_run` and `close_run` control messages
before and after the sequnce generated by the enclosed code:

.. code-block:: python

    # Using 'bpp.run_decorator' to define a run

    import bluesky.preprocessors as bpp
    from bluesky import RunEngine

    RE = RunEngine({})

    @bpp.run_decorator(md={})  # 'md' - metadata to be added to the 'start' document
    def sample_plan():
        ...
        < code that controls execution of the scan >
        ...

    RE(sample_plan())

The rules for basic Bluesky plans require that the currently running scan is closed before
the next scan is opened, therefore the following code works:

.. code-block:: python

    # This code works, since the first run is closed before the second one is opened

    import bluesky.plan_stubs as bps
    from bluesky import RunEngine

    RE = RunEngine({})

    def sample_plan():
        yield from bps.open_run(md={})
        < code that controls execution of the scan >
        yield from bps.close_run()  # Closing the first run (scan)
        yield from bps.open_run(md={})  # Opening the second run (scan)
        < code that controls execution of the scan >
        yield from bps.close_run()

    RE(sample_plan())

but the following code fails:

.. code-block:: python

    # This code fails, since the second run is opened before the first run is closed

    import bluesky.plan_stubs as bps
    from bluesky import RunEngine

    RE = RunEngine({})

    def sample_plan():
        yield from bps.open_run(md={})  # Opening the first run
        < code that controls execution of the scan >
        yield from bps.open_run(md={})  # Opening the second run before the first one is closed
        < code that controls execution of the scan >
        yield from bps.close_run()
        yield from bps.close_run()

    RE(sample_plan())


Note, that the preassembled plans, such as `bluesky.plans.count` or `bluesky.plans.list_scan`,
are complete single-run plans, enclosed in `open_run` and `close_run` messages, therefore
the following code fails as well:

.. code-block:: python

    # This code fails while attempting to start a preassembled plan from an open run

    import bluesky.plan_stubs as bps
    from bluesky.plans import count
    from bluesky import RunEngine

    RE = RunEngine({})

    def sample_plan():
        yield from bps.open_run(md={})  # Starting the first run
        < code that controls execution of the scan >
        yield from bpp.count(<some arguments>)  # Attempting to run a preassembled plan from an open run
        yield from bps.close_run()

    RE(sample_plan())

An example of the situation when a preassembled plan is called from another open run is
when a preassembled plan is included in a suspender pre- or post-plan. When the suspender is
triggered, the current run is interrupted (not closed) and the pre- or post-plan attempts to open
another run (the mechanism is the same as in the case of nested runs, see below). As a result,
Run Engine fails for the same reason as in the two previous code examples. The new multi-run plan
Bluesky features allow to implement nested plans, as well as include full-featured scans
in pre- and post-plans.

Bluesky Features for Support of Multi-run Plans
-----------------------------------------------

In order to handle simultaneously open runs within a plan, Run Engine is looking at the run key attribute
of each control message to decide which scan is currently being executed. The default value for the run key
is `None`, but it could be manually set in the plan for any block of code which define the run. A run key
value may be of any type, but it is **strongly** recommended that manually assigned run keys are
human-readable informative strings.

The new 'inner' run can be opened from within the 'outer' run only if the run keys of the 'inner' and
'outer' scans are different. Otherwise the plan exectuion fails.

The run key is used by Run Engine

* to maintain the state of each run independently from other open runs;

* to include run metadata, such as scan ID and UID, into the emitted documents. (Metadata is then used
  to route the documents to the appropriate callbacks. If documents are saved using Databroker, the metadata
  allows to associate documents with runs and retrieve run data from the database.)

Run key is assigned to a block of code using `bpp.set_run_key_wrapper` or `@bpp.set_run_key_decorator`:

.. code-block:: python

    import bluesky.preprocessors as bpp
    from bluesky import RunEngine

    # Using decorator
    @bpp.set_run_key_decorator("run_key_example_1")
    @bpp.run_decorator(md={})
    def sample_plan():
        ...
        < code that controls execution of the run >
        ...

    RE(sample_plan())

    from bluesky.plans import scan
    from ophyd.sim import hw
    det, motor = hw().det, hw().motor

    # Using wrapper
    s = scan([det], motor, -1, 1, 10)
    s_wrapped = bpp.set_run_key_wrapper(s, "run_key_example_2")
    RE(s_wrapped)

The implementation of `@bpp.set_run_key_decorator` and `bpp.set_run_key_wrapper` is
replacing the default value `None` of the attribute `run` in each message generated within
the enclosed block with the user-defined run key.

The `@bpp.set_run_key_decorator` and `bpp.set_run_key_wrapper` are primarily intended
to be applied to a function that contains a run implementation, but may be also used
with any block of plan code. For example, one may write a plan that simultaneously
opens multiple runs and executes them in parallel by generating groups of messages
with run ids of the open scans. This is currently not recommended and should be attempted
only at the developer's own risk.

Plans with Sequential Runs
---------------------------

Sequential calling of multiple runs is supported by older versions of Bluesky. There is no need
to use multi-run plan features if runs are not overlapping (the next run is opened only after
the previous run is closed), but run keys still can be assigned to all or some runs if needed.

In the following example, two preassembled plans are called in sequence. Run Engine is subscribed to
a single instance of BestEffortCallback, which is set up to display data specific for each run
when the run opened.

.. literalinclude:: examples/multi_run_plans_sequential.py

.. ipython:: python
    :suppress:

    %run -m multi_run_plans_sequential

.. ipython:: python

    RE(plan_sequential_runs(10))

Plans with Nested Runs
----------------------

The following example illustrates the use of `@bpp.set_run_key_decorator` to implement two nested runs:
the 'outer' run interrupts measurements, calls the 'inner' run and then completes the measurements.
The 'outer' and 'inner' runs are assigned different run ids ('run_1' and 'run_2'). Note that
the `@bpp.set_run_key_decorator` for the 'outer' run does not overwrite the run id of the 'inner' scan,
despite the fact that it is generated inside the enclosed code, since the decorator is designed to replace
the run id attribute of the message only if it has the default value of `None`, i.e. the run id of
a message can be replaced by the decorator only the first time it is processed by the decorator.

If multiple runs are to be opened simultaneously, each run needs to be subscribed to its own instance
of callback. Standard RunEngine subscription mechanism does not provide this capability. Instead,
subscription should be performed via `RunRouter`. The code in the following example demonstrates how
to use `BestEffortCallback` to monitor data from multiple nested runs.

.. literalinclude:: examples/multi_run_plans_nested.py

The output of the plan contains data from two runs with each run assigned its own ID and UID. The tables
for the runs are printed by two separate instances of `BestEffortCallback`. The data from two tables
is printed in the order of acquisition: the table for the 'inner' run is printed in the gap of
the table for the 'outer' run.

.. ipython:: python
    :suppress:

    %run -m multi_run_plans_nested

.. ipython:: python

    RE(sim_plan_outer(10))

The wrapper `bpp.set_run_key_wrapper` can be used instead of the decorator. For example
the run `sim_plan_inner` from the previous example can be rewritten as follows:

.. code-block:: python

    def sim_plan_inner(npts):
        def f():
            for j in range(npts):
                yield from bps.mov(hw.motor1, j * 0.1 + 1, hw.motor2, j * 0.2 - 2)
                yield from bps.trigger_and_read([hw.motor1, hw.motor2, hw.det2])
        f = bpp.run_wrapper(f(), md={})
        return bpp.set_run_key_wrapper(f, "run_2")

Subscription to callbacks via RunRouter provides flexibility to subscribe each run
to its own set of callbacks. In the following example `run_key` is added to the start
document metadata and used to distinguish between two runs in the function factory that
performs callback subscriptions.

.. literalinclude:: examples/multi_run_plans_select_cb.py

.. ipython:: python
    :suppress:

    %run -m multi_run_plans_select_cb

.. ipython:: python

    RE(sim_plan_outer(10))

In some cases it may be necessary to implement a run that could be interrupted
and a new instance of the same run started. For example, the suspender pre- or post-plan
may contain a run, which takes substantial time to execute. Such run may be interrupted
if the suspender is repeatedly triggered. This will cause another instance of the pre-
or post-plan to be started while the first one is still in the open state. This process
is similar to recursive calling of the run (run which includes instructions to call
itself). Recursive calls are possible if unique run key is assigned to a run each
time it is started.

The following example illustrates dynamic generation of run keys. The plan may have no practical purpose
besides demonstration of the principle. The plan is calling itself recursively multiple times until
the global counter `n_calls` reaches the maximum value of `n_calls_max`. The unique run key is generated
before at each call.

.. literalinclude:: examples/multi_run_plans_recursive.py

.. ipython:: python
    :suppress:

    %run -m multi_run_plans_recursive

.. ipython:: python

    RE(sim_plan_recursive(4))

The identical result can be achieved by using `bpp.set_run_key_wrapper()`:

.. code-block:: python

    # Call counter and the maximum number calls
    n_calls, n_calls_max = 0, 3

    def sim_plan_recursive(npts):
        global n_calls, n_calls_max

        n_calls += 1  # Increment counter
        if n_calls <= n_calls_max:
            # Generate unique key for each run. The key generation algorithm
            #   must only guarantee that execution of the runs that are assigned
            #   the same key will never overlap in time.
            run_key = f"run_key_{n_calls}"

            @bpp.run_decorator(md={})
            def plan(npts):

                for j in range(int(npts/2)):
                    yield from bps.mov(hw.motor1, j * 0.2)
                    yield from bps.trigger_and_read([hw.motor1, hw.det1])

                # Different parameter values may be passed to the recursively called plans
                yield from sim_plan_recursive(npts + 2)

                for j in range(int(npts/2), npts):
                    yield from bps.mov(hw.motor1, j * 0.2)
                    yield from bps.trigger_and_read([hw.motor1, hw.det1])

            yield from bpp.set_run_key_wrapper(plan(npts), run_key)
Comparison with SPEC
====================

`SPEC <https://www.certif.com/content/spec/>`_ is a popular software package
for instrument control and data acquisition. Many users in the synchrotron
community, from which bluesky originated, know SPEC and ask what differentiates
and motivates bluesky. Answering this question in an informed and unbiased way
is difficult, and we welcome any corrections.

There are many good features of SPEC that have been incorporated into
bluesky, including:

* Simple commands for common experiments, which can be used as building blocks
  for more complex procedures
* Easy hardware configuration
* Interruption handling (Ctrl+C)
* Integration with EPICS (and potentially other instrument control systems)
* "Pseudomotors" presenting virtual axes
* Integration with reciprocal space transformation code

Bluesky has also addressed certain limitations of SPEC. In fairness to SPEC, we
have the benefit of learning from its decades of use, and we are standing on
the shoulders of the modern open-source community.

* Bluesky is free and open-source in all aspects. Macros from the SPEC user
  community are open-source, but the core SPEC C source code is closed and not
  free.
* Bluesky provides more control over the console output and plotting.
* SPEC was designed before large area detectors existed. Ingesting area
  detector data is possible, but *ad hoc*. In bluesky, area detectors and other
  higher-dimensional inputs are integrated naturally.
* SPEC writes to a custom text-based format (a "SPEC file"). Bluesky can
  write---in real time or *post facto*---to any format.
* SPEC has a simulation mode. Bluesky allows users to incorporate much richer
  simulation capabilities (about which more below) but, as of version 0.9.0,
  provides less than SPEC out of the box.

Using Python, a general-purpose programming language, gives several immediate
advantages:

* Python checks syntax automatically.
* Python provides tools for interactive debugging.
* There are many more resources for learning Python.
* The language is more flexible.
* It's easy to integrate with the scientific Python ecosystem.

Bluesky tries to go further than SPEC in some regards:

* Complex custom procedures are easier to express.
* Automated "suspension" (pausing and resuming) is consistent and easier to
  manage.
* The prevailing model in SPEC is to collect data as a step scan. Other types
  of scans---such as fly scans or asynchronous monitoring---can be
  done, but they are *ad hoc*. Bluesky supports several modalities of data
  acquisition with equal ease.
* Bluesky can acquire multiple asynchronous, uncoordinated streams of data and
  represent them in a simple :doc:`event-based data model <documents>`.
* It is easy to build tools that inspect/simulate a procedure before it is run
  to check for safety, estimate time to completion, or visualize its behavior.
* Bluesky is a library that works well interactively but can also be used
  programmatically in scripts or other libraries.
* Users can add arbitrary metadata with rich semantics, including large arrays
  (such as masks) or nested mappings.
* Bluesky is a holistic solution for data acquisition and management. Users can
  push live streaming data directly into their data processing and analysis
  pipelines and/or export it into a file.

On the other hand, one major advantage of SPEC over bluesky is its maturity.
SPEC is battle-hardened from decades of use at many facilities, and it has a
large user community. Bluesky is a young project.

A Remark About Syntax
---------------------

SPEC users immediately notice that simple bluesky commands are more verbose
than their counterparts in SPEC. This is a trade-off we have made in choosing a
more expressive, general-purpose language over a single-purpose command line
interface. That "easy integration with scientific libraries" comes at the cost
of some parentheses and commas. Some of the difference is also due to the
richer abstractions required to capture the complexity of modern hardware.  The
simplest commands are made less terse, but more interesting commands are made
much easier to express. Of course, users can save time by using tab-completion
and by accessing previous commands with the up arrow key.
Event Descriptors
=================

In the section on :doc:`documents`, we gave an overview of the four kinds of
document. We presented an example Run Start, Event, and Run Stop, but we
deferred detailed discussion of the Event Descriptor.

Recall our example 'event' document.

.. code-block:: python

    # 'event' document (same as above, shown again for reference)
    {'data':
        {'temperature': 5.0,
          'x_setpoint': 3.0,
          'x_readback': 3.05},
     'timestamps':
        {'temperature': 1442521007.9258342,
         'x_setpoint': 1442521007.5029348,
         'x_readback': 1442521007.5029348},
     'time': 1442521007.3438923,
     'seq_num': 1,
     'uid': '<randomly-generated unique ID>', 
     'descriptor': '<reference to a descriptor document>'}

Typically, an experiment generates multiple event documents with the same data
keys. For example, there might be ten sequential readings, generating ten event
documents like the one above --- with different readings and timestamps but
identical data keys. All these events refer back to a 'descriptor' with
metadata about the data keys and the configuration of the devices involved.

.. note:: 

    We got the term "data keys" from ``event['data'].keys()``. Again, in our
    example, the data keys are ``['temperature', 'x_setpoint', 'x_readback']``

Data Keys
---------

First, the descriptor provides metadata about each data key.

* dtype --- 'number', 'string', 'array', or 'object' (dict)
* shape --- ``None`` or a list of dimensions like ``[5, 5]`` for a 5x5 array
* source --- a description of the hardware that uniquely identifies it, such as
  an EPICS Process Variable
* (optional) external --- a string specifying where external data, such as a
  large image array, is stored

Arbitrary additional fields are allowed, such as precision or units.
The RunEngine obtains this information from each device it sees by calling
``device.describe()``.

.. code-block:: python

    # excerpt of a 'descriptor' document
    {'data_keys':
        {'temperature':
            {'dtype': 'number',
             'source': '<descriptive string>',
             'shape': [],
             'units': 'K',
             'precision': 3},
         'x_setpoint':
            {'dtype': 'number',
             'source': '<descriptive string>',
             'shape': [],
             'units': 'mm',
             'precision': 2},
         'x_readback':
            {'dtype': 'number',
             'source': '<descriptive string>',
             'shape': [],
             'units': 'mm',
             'precision': 2}},
     ...}

Object Keys
-----------

The ``object_keys`` provide an association between each device and its data keys.

This is needed because a given device can produce multiple data keys. For
example, suppose the ``x_readback`` and ``x_setpoint`` data keys in our example
came from the same device, a motor named ``'x'``.

.. code-block:: python

    # excerpt of a 'descriptor' document
    {'object_keys':
        {'x': ['x_setpoint', 'x_readback'],
         'temp_ctrl': ['temperature']},
     ...}

Specifically, it maps ``device.name`` to ``list(device.describe())``.

Configuration
-------------

Complex devices often have many parameters that do not need to be read anew
with every data point. They are "configuration," by which we mean they don't
typically change in the middle of a run. A detector's exposure time is usually
(but not always) in this category.

Devices delineate between the two by providing two different methods that the
RunEngine can call: ``device.read()`` returns normals readings that are *not*
considered configuration; ``device.read_configuration()`` returns the readings
that are considered configuration.

The first time during a run that the RunEngine is told to read a device, it
reads the device's configuration also. The return value of
``device.describe_configuration()`` is recorded in
``configuration[device.name]['data_keys']``. The return value of
``device.read_configuration()`` is collated into
``configuration[device.name]['data']`` and
``configuration[device.name]['timestamps']``.

In this example, ``x`` has one configuration data key, and ``temp_ctrl``
happens to provide no configuration information.

.. code-block:: python

    # excerpt of a 'descriptor' document
    {'configuration':
        {'x':
           {'data': {'offset': 0.1},
            'timestamps': {'offset': 1442521007.534918},
            'data_keys':
               {'offset':
                   {'dtype': 'number',
                    'source': '<descriptive string>',
                    'shape': [],
                    'units': 'mm',
                    'precision': 2}}},
         'temp_ctrl':
            {'data': {},
             'timestamps': {}
             'data_keys': {}}}
     ...}

Hints
-----

This is an experimental feature. Devices can provide information via a
``hints`` attribute that is stored here. See :ref:`hints`.

.. code-block:: python

    # excerpt of a 'descriptor' document
     {'hints':
        {'x' {'fields': ['x_readback']},
         'temp_ctrl': {'fields': ['temperature']}}
      ...}


Complete Sample
---------------

Taken together, our example 'descriptor' document looks like this.

.. code-block:: python

    # complete 'descriptor' document
    {'data_keys':
        {'temperature':
            {'dtype': 'number',
             'source': '<descriptive string>',
             'shape': [],
             'units': 'K',
             'precision': 3},
         'x_setpoint':
            {'dtype': 'number',
             'source': '<descriptive string>',
             'shape': [],
             'units': 'mm',
             'precision': 2}},
         'x_readback':
            {'dtype': 'number',
             'source': '<descriptive string>',
             'shape': [],
             'units': 'mm',
             'precision': 2}},

     'object_keys':
        {'x': ['x_setpoint', 'x_readback'],
         'temp_ctrl': ['temperature']},

     'configuration':
         {'x':
            {'data': {'offset': 0.1},
             'timestamps': {'offset': 1442521007.534918},
             'data_keys':
                {'offset':
                    {'dtype': 'number',
                     'source': '<descriptive string>',
                     'shape': [],
                     'units': 'mm',
                     'precision': 2}
          'temp_ctrl':
            {'data': {},
             'timestamps': {}
             'data_keys': {}}}
         }

     'hints':
        {'x' {'fields': ['x_readback']},
         'temp_ctrl': {'fields': ['temperature']}}

     'time': 1442521007.3438923,
     'uid': '<randomly-generated unique ID>', 
     'run_start': '<reference to the start document>'}
*******************************
IPython 'Magics' [Experimental]
*******************************

.. warning::

    This section covers an experimental feature of bluesky. It may be altered
    or removed in the future.

What Are 'Magics'?
------------------

IPython is an interactive Python interpreter designed for and by scientists. It includes a feature called "magics" --- convenience commands that aren't part of the Python language itself. For example, ``%history`` is a magic:

.. ipython:: python

    a = 1
    b = 2
    %history

The IPython documentation documents the
`complete list of built-in magics <https://ipython.readthedocs.io/en/stable/interactive/magics.html>`_
and, further,
`how to define custom magics <https://ipython.readthedocs.io/en/stable/config/custommagics.html>`_.

Bluesky Magics
--------------

Bundled with bluesky are some IPython magics. They are intended for maintenance
tasks or casual sanity checks.  **Intentionally, none of the magics save data;
for that you should use the RunEngine and plans.**

To use the magics, first register them with IPython:

.. ipython:: python

    from bluesky.magics import BlueskyMagics
    get_ipython().register_magics(BlueskyMagics)

For this example we'll use some simulated hardware.

.. ipython:: python

    from ophyd.sim import motor1, motor2

Moving a Motor
~~~~~~~~~~~~~~

Suppose you want to move a motor interactively. You can use the ``%mov`` magic:

.. ipython:: python

    %mov motor1 42

Where ``motor1`` refers to the actual ophyd object itself.
This is equivanent to:

.. code-block:: python

    from bluesky.plan_stubs import mv

    RE(mv(motor1, 42))

but less to type. There is also a ``%movr`` magic for "relative move". They can
move multiple devices in parallel like so:

.. ipython:: python

    %mov motor1 -3 motor2 3


Note: Magics has changed from version v1.3.0 onwards. The previous method will
be described in the next section.

Taking a reading using ``%ct`` (Post v1.3.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before we may make use of the power of magics for counting, we must "label"
this hardware. To add a label, we must give hardware a ``labels={'mylabel'}``
keyword argument. For example, here we initialize five simulated signals: two
motors, a shutter motor, an area detector and a point detector:

.. ipython:: python

    import numpy as np
    from ophyd.sim import SynAxis, SynSignal
    motor1 = SynAxis(name='motor1', labels={'motors', 'scan_motors'})
    motor2 = SynAxis(name='motor2', labels={'motors', 'scan_motors'})
    shutter_motor = SynAxis(name='shutter_motor', labels={'motors', 'shutter_motors'})
    # create a fake area detector that returns a 2x2 array
    area_detector = SynSignal(func=lambda: np.random.random((2, 2)),
                              name='adet1', labels={'detectors', 'area_detectors'})
    point_detector = SynSignal(func=lambda: np.random.random((1,)),
                               name='pointdet1', labels={'detectors', 'point_detectors'})

Now we have detectors and motors, with proper labels.

Now suppose you want to take a quick reading of some devices and print the
results to the screen without saving them or doing any fancy processing. Use
the ``%ct`` magic:

.. ipython:: python

    %ct area_detectors

Where the names after count are a list of whitespace separated labels. In this
case, only ``area_detector`` will be counted.

Running ``%ct`` without arguments looks for the ``detectors`` label by default:

.. ipython:: python

    %ct

In this case, we count both on the area detector and the point detector.


Aside on the automagic feature in IPython
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If IPython’s ‘automagic’ feature is enabled, IPython will even let you drop the
``%`` as long as the meaning is unambiguous:

.. ipython:: python

    ct
    ct = 3  # Now ct is a variable so automagic will not work...
    ct
    # ... but the magic still works.
    %ct

For what it’s worth, we recommend disabling 'automagic'. The ``%`` is useful
for flagging what follows as magical, non-Python code.

Listing available motors using ``%wa`` (Post v1.3.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Finally, the ``%wa`` magic displays the a list of labeled devices.

.. ipython:: python

    %wa scan_motors

will display all motors used for a scan.
If blank, will print all labeled devices.

.. ipython:: python

    %wa

Note: It is possible to give a device more than one label. Thus it is possible
to have the same device in more than one list when calling ``%wa``. It is up to
the user to decide whether they want overlapping labels or not.

    
Comparison with SPEC
~~~~~~~~~~~~~~~~~~~~

The names of these magics, and the order of the parameters they take, are meant
to feel familiar to users of :doc:`SPEC <comparison-with-spec>`.

Again, they must be registered with IPython before they can be used:

.. code-block:: python

    from bluesky.magics import BlueskyMagics
    get_ipython().register_magics(BlueskyMagics)


Taking a reading using ``%ct`` (Pre v1.3.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Previously, you could set a default list of detectors and them use ``%ct``
without any parameters. This behaviour is deprecated. Do not use this:

.. ipython:: python
    :okwarning:

    BlueskyMagics.detectors = [area_detector, point_detector]
    %ct

This is no longer supported.

Listing available motors using ``%wa`` (Pre v1.3.0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Previously, it was possible to supply a list of motors. This feature is also
deprecated. Do not use this:

.. ipython:: python
    :okwarning:

    BlueskyMagics.positioners = [motor1, motor2]
    %wa

======================================================================= ==============================
Magic                                                                   Plan Invoked
======================================================================= ==============================
``%mov``                                                                :func:`~bluesky.plan_stubs.mv`
``%movr``                                                               :func:`~bluesky.plan_stubs.mvr`
``%ct``                                                                 :func:`~bluesky.plans.count`
``%wa``                                                                 ("where all") Survey positioners*
======================================================================= ==============================
********
Tutorial
********

Before You Begin
================

.. note::

    NSLS-II deploys a free, public "sandbox" for trying the software in the
    browser using Jupyter notebooks. There will be no need to install any
    software, and you can skip the rest of this section.  Go to
    `https://try.nsls2.bnl.gov <https://try.nsls2.bnl.gov>`_.

* You will need Python 3.5 or newer. From a shell ("Terminal" on OSX,
  "Command Prompt" on Windows), check your current Python version.

  .. code-block:: bash

    python3 --version

  If that version is less than 3.5, you must update it.

  We recommend install bluesky into a "virtual environment" so this
  installation will not interfere with any existing Python software:

  .. code-block:: bash

    python3 -m venv ~/bluesky-tutorial
    source ~/bluesky-tutorial/bin/activate

  Alternatively, if you are a
  `conda <https://conda.io/docs/user-guide/install/download.html>`_ user,
  you can create a conda environment:

  .. code-block:: bash

    conda create -n bluesky-tutorial "python>=3.5"
    conda activate bluesky-tutorial

* Install the latest versions of bluesky and ophyd. Also install the databroker
  unless you plan to skip the sections about accessing saved data. If you want
  to follow along with the visualization examples, install matplotlib and
  PyQt5. Finally, install IPython (a Python interpreter designed by scientists
  for scientists).

  .. code-block:: bash

     python3 -m pip install --upgrade bluesky ophyd databroker matplotlib pyqt5 ipython pyepics

  Alternatively, if you are a conda user and you prefer conda packages, you can
  use:

  .. code-block:: bash

     conda install -c nsls2forge bluesky ophyd databroker matplotlib pyqt=5 ipython

* Start IPython:

  .. code-block:: python

     ipython --matplotlib=qt5

  The flag ``--matplotlib=qt5`` is necessary for live-updating plots to work.

  Or, if you wish you use bluesky from a Jupyter notebook, install a kernel like
  so:

  .. code-block:: python

     ipython kernel install --user --name=bluesky-tutorial --display-name "Python (bluesky)"

  You may start Jupyter from any environment where it is already installed, or
  install it in this environment alongside bluesky and run it from there:

  .. code-block:: python

     pip install notebook
     jupyter notebook

If you get lost or confused...
==============================

...then we want to know! We have a friendly
`chat channel <https://gitter.im/NSLS-II/DAMA>`_, or you can
`file a bug <https://github.com/NSLS-II/Bug-Reports/issues>`_ to let us know
where our documentation could be made more clear.

.. _tutorial_run_engine_setup:

The RunEngine
=============

Bluesky encodes an experimental procedure as a *plan*, a sequence of
atomic instructions. The *RunEngine* is an interpreter for plans. It lets
us focus on the logic of our experimental procedure while it handles important
technical details consistently: it communicates with hardware, monitors for
interruptions, organizes metadata and data, coordinates I/O, and ensures that
the hardware is left in a safe state at exit time.

This separation of the executor (the RunEngine) from the instruction set (the
plan) pays off in several ways, as we will see in the examples that follow.

.. note::

    If you are a visiting user at a facility that runs bluesky, you can skip
    this section and go straight to :ref:`common_experiments`. A RunEngine will
    have already been configured for you. **If you ignore this and define your
    own, you may be overriding pre-configured defaults, which can result in
    data loss.**

    To check, type ``RE``. If a RunEngine has already been configured, you
    should get something like:

    .. ipython::
        :verbatim:

        In [1]: RE
        Out[1]: <bluesky.run_engine.RunEngine at 0x10fd1d978>

    and you should skip the rest of this section. But if this gives you a
    ``NameError``, you'll need to finish this section.

Create a RunEngine:

.. code-block:: python

    from bluesky import RunEngine

    RE = RunEngine({})

.. ipython:: python
    :suppress:

    # for use in later demos
    from bluesky import RunEngine
    RE = RunEngine({})


This RunEngine is ready to use --- but if you care about visualizing or saving
your data, there is more to do first....

During data acquisition, the RunEngine dispatches a live stream of metadata and
data to one or more consumers ("callbacks") for in-line data processing and
visualization and long-term storage. Example consumers include a live-updating
plot, a curve-fitting algorithm, a database, a message queue, or a file in your
preferred format. See :doc:`callbacks` for more detail.

Prepare Live Visualization
--------------------------

To start, let's use the all-purpose
:class:`~bluesky.callback.best_effort.BestEffortCallback`.

.. code-block:: python

    from bluesky.callbacks.best_effort import BestEffortCallback
    bec = BestEffortCallback()

    # Send all metadata/data captured to the BestEffortCallback.
    RE.subscribe(bec)

    # Make plots update live while scans run.
    from bluesky.utils import install_kicker
    install_kicker()

.. ipython:: python
    :suppress:

    # for use in later demos
    from bluesky.callbacks.best_effort import BestEffortCallback
    bec = BestEffortCallback()
    RE.subscribe(bec)

The :class:`~bluesky.callback.best_effort.BestEffortCallback` will receive the
metadata/data in real time and produce plots and text, doing its best to
provide live feedback that strikes the right balance between "comprehensive"
and "overwhelming."

For more tailored feedback, customized to a particular experiment, you may
configure custom callbacks. Start by reading up on :doc:`documents`, the
structure into which bluesky organized metadata and data captured during an
experiment. But for this tutorial and for many real experiments, the
:class:`~bluesky.callback.best_effort.BestEffortCallback` will suffice.

Prepare Data Storage
--------------------

.. _databroker_setup:

The `databroker <https://nsls-ii.github.io>`_, a library developed in tandem
with bluesky, is an interface to searchable storage for metadata and data
generated by bluesky. For this tutorial, we will spin up a databroker backed by
temporary files.

.. code-block:: python

    from databroker import Broker
    db = Broker.named('temp')

    # Insert all metadata/data captured into db.
    RE.subscribe(db.insert)

.. ipython:: python
    :suppress:

    # for use in later demos
    from databroker import Broker
    db = Broker.named('temp')
    RE.subscribe(db.insert)

.. warning::

    **This example makes a temporary database. Do not use it for important
    data.** The data will become difficult to access once Python exits or the
    variable ``db`` is deleted. Running ``Broker.named('temp')`` a second time
    creates a fresh, separate temporary database.

Add a Progress Bar
------------------

Optionally, you can configure a progress bar.

.. code-block:: python

    from bluesky.utils import ProgressBarManager
    RE.waiting_hook = ProgressBarManager()

See :doc:`progress-bar` for more details and configuration.

Let's take some data!

.. _common_experiments:

Common Experiments ("Plans")
============================

Read Some Detectors
-------------------

Begin with a very simple experiment: trigger and read some detectors. Bluesky
calls this "counting", a term of art inherited from the spectroscopy
community.

For this tutorial, we will not assume that you have access to real detectors or
motors. In the examples that follow, we will use simulated hardware from
`ophyd <https://nsls-ii.github.io/ophyd>`_, a library developed in tandem with
bluesky. In a :ref:`later section <tutorial_device>` we will see what it looks
like to configure *real* hardware with ophyd.

.. code-block:: python

    from ophyd.sim import det1, det2  # two simulated detectors

Using the RunEngine, ``RE``, "count" the detectors:

.. code-block:: python

    from bluesky.plans import count
    dets = [det1, det2]   # a list of any number of detectors

    RE(count(dets))

Demo:

.. ipython:: python
    :suppress:

    from bluesky.plans import count
    from ophyd.sim import det1, det2
    dets = [det1, det2]

.. ipython:: python

    RE(count(dets))

A key feature of bluesky is that these detectors could be simple photodiodes or
complex CCDs. All of those details are captured in the implementation of the
Device. From the point of view of bluesky, detectors are just Python objects
with certain methods.

See :func:`~bluesky.plans.count` for more options. You can also view this
documentation in IPython by typing ``count?``.

Try the following variations:

.. code-block:: python

    # five consecutive readings
    RE(count(dets, num=5))

    # five sequential readings separated by a 1-second delay
    RE(count(dets, num=5, delay=1))

    # a variable delay
    RE(count(dets, num=5, delay=[1, 2, 3, 4]))

The :func:`~bluesky.plans.count` function (more precisely, Python *generator
function*) is an example of a *plan*, a sequence of instructions encoding an
experimental procedure. We'll get a better sense for why this design is useful
as we continue. Briefly, it empowers us to:

* Introspect the instructions before we execute them, checking for accuracy,
  safety, estimated duration, etc.
* Interrupt and "rewind" the instructions to a safe point to resume from,
  both interactively and automatically (e.g. in the middle of the night).
* Reuse a generic set of instructions on different hardware.
* Modify the instructions programmatically, such as inserting a set of
  baseline readings to be taken automatically before every experiment.

.. warning::

    Notice that entering a plan by itself doesn't do anything:

    .. ipython:: python
        :suppress:

        from bluesky.plans import count
        from ophyd.sim import det
        dets = [det]

    .. ipython:: python

        count(dets, num=3)

    If we mean to *execute* the plan, we must use the RunEngine:

    .. ipython:: python

        RE(count(dets, num=3))

Scan
----

Use :func:`~bluesky.plans.scan` to scan ``motor`` from ``-1`` to ``1`` in ten
equally-spaced steps, wait for it to arrive at each step, and then trigger and
read some detector, ``det``.

.. code-block:: python

    from ophyd.sim import det, motor
    from bluesky.plans import scan
    dets = [det]   # just one in this case, but it could be more than one

    RE(scan(dets, motor, -1, 1, 10))

.. ipython:: python
    :suppress:

    from bluesky.plans import scan
    from ophyd.sim import det, motor
    dets = [det]

.. ipython:: python

    RE(scan(dets, motor, -1, 1, 10))

.. plot::

    from bluesky.plans import scan
    from ophyd.sim import det, motor
    dets = [det]
    RE(scan(dets, motor, -1, 1, 10))

Again, a key feature of bluesky is that ``motor`` may be any "movable" device,
including a temperature controller, a sample changer, or some pseudo-axis. From
the point of view of bluesky and the RunEngine, all of these are just objects
in Python with certain methods.

In addition the producing a table and plot, the
:class:`~bluesky.callback.best_effort.BestEffortCallback` computes basic peak
statistics. Click on the plot area and press Shift+P ("peaks") to visualize
them over the data. The numbers (center of mass, max, etc.) are available in a
dictionary stashed as ``bec.peaks``. This is updated at the end of each run.
Of course, if peak statistics are not applicable, you may just ignore this
feature.

Use :func:`~bluesky.plans.rel_scan` to scan from ``-1`` to ``1`` *relative to
the current position*.

.. code-block:: python

    from bluesky.plans import rel_scan

    RE(rel_scan(dets, motor, -1, 1, 10))

Use :func:`~bluesky.plans.list_scan` to scan points with some arbitrary
spacing.

.. code-block:: python

    from bluesky.plans import list_scan

    points = [1, 1, 2, 3, 5, 8, 13]

    RE(list_scan(dets, motor, points))

For a complete list of scan variations and other plans, see :doc:`plans`.

.. _tutorial_multiple_motors:

Scan Multiple Motors Together
-----------------------------

There are two different things we might mean by the phrase "scan multiple
motors 'together'". In this case we mean that we move N motors along a line in
M steps, such as moving X and Y motors along a diagonal. In the other case, we
move N motors through an (M_1 x M_2 x ... x M_N) grid; that is addressed in the
next section.

SPEC users may recognize this case as analogous to an "a2scan" or "d2scan", but
with an arbitrary number of dimensions, not just two.

We'll use the same plans that we used in the previous section. (If you already
imported them, there is no need to do so again.)

.. code-block:: python

    from bluesky.plans import scan, rel_scan

We'll use two new motors and a new detector that is coupled to them via
a simulation. It simulates a 2D Gaussian peak centered at ``(0, 0)``.
Again, we emphasize that these "motors" could be anything that can be "set"
(temperature controller, pseudo-axis, sample changer).

.. code-block:: python

    from ophyd.sim import det4, motor1, motor2
    dets = [det4]   # just one in this case, but it could be more than one

The plans :func:`~bluesky.plans.scan` and  :func:`~bluesky.plans.rel_scan`
accept multiple motors.

.. code-block:: python

    RE(scan(dets,
            motor1, -1.5, 1.5,  # scan motor1 from -1.5 to 1.5
            motor2, -0.1, 0.1,  # ...while scanning motor2 from -0.1 to 0.1
            11))  # ...both in 11 steps

The line breaks are intended to make the command easier to visually parse. They
are not technically meaningful; you may take them or leave them.

Demo:

.. ipython:: python
    :suppress:

    from bluesky.plans import scan
    from ophyd.sim import det4, motor1, motor2
    dets = [det4]

.. ipython:: python

    RE(scan(dets,
            motor1, -1.5, 1.5,  # scan motor1 from -1.5 to 1.5
            motor2, -0.1, 0.1,  # ...while scanning motor2 from -0.1 to 0.1
            11))  # ...both in 11 steps

.. plot::

    from bluesky.plans import scan
    from ophyd.sim import det4, motor1, motor2
    dets = [det4]
    RE(scan(dets,
            motor1, -1.5, 1.5,  # scan motor1 from -1.5 to 1.5
            motor2, -0.1, 0.1,  # ...while scanning motor2 from -0.1 to 0.1
            11))  # ...both in 11 steps

This works for any number of motors, not just two. Try importing ``motor3``
from ``ophyd.sim`` and running a 3-motor scan.

To move motors along arbitrary trajectories instead of equally-spaced points,
use :func:`~bluesky.plans.list_scan` and :func:`~bluesky.plans.rel_list_scan`.

.. code-block:: python

    from bluesky.plans import list_scan

    # Scan motor1 and motor2 jointly through a 5-point trajectory.
    RE(list_scan(dets, motor1, [1, 1, 3, 5, 8], motor2, [25, 16, 9, 4, 1]))

Demo:

.. ipython:: python
   :suppress:

   from bluesky.plans import list_scan

.. ipython:: python

    RE(list_scan(dets,
                 motor1, [1, 1, 3, 5, 8],
                 motor2, [25, 16, 9, 4, 1]))

.. plot::

    from bluesky.plans import list_scan
    from ophyd.sim import det4, motor1, motor2
    dets = [det4]
    RE(list_scan(dets,
                 motor1, [1, 1, 3, 5, 8],
                 motor2, [25, 16, 9, 4, 1]))

Scan Multiple Motors in a Grid
------------------------------

In this case scan N motors through an N-dimensional rectangular grid. We'll use
the same simulated hardware as in the previous section:

.. code-block:: python

    from ophyd.sim import det4, motor1, motor2
    dets = [det4]   # just one in this case, but it could be more than one

We'll use a new plan, named :func:`~bluesky.plans.grid_scan`.

.. code-block:: python

    from bluesky.plans import grid_scan

Let's start with a 3x5x5 grid.

.. code-block:: python

    RE(grid_scan(dets,
                 motor1, -1.5, 1.5, 3,  # scan motor1 from -1.5 to 1.5 in 3 steps
                 motor2, -0.1, 0.1, 5,  # scan motor2 from -0.1 to 0.1 in 5 steps
                 motor3, 10, -10, 5))  # scan motor3 from 10 to -10 in 5 steps

The order of the motors controls how the grid is traversed. The "slowest" axis
comes first. Numpy users will appreciate that this is consistent with numpy's
convention for indexing multidimensional arrays.

The optional parameter ``snake_axes`` can be used to control which motors'
trajectories "snake" back and forth. A snake-like path is usually more
efficient, but it is not suitable for certain hardware, so it is disabled by
default. To enable snaking for specific axes, give a list like
``snake_axes=[motor2]``.  Since the first (slowest) axis is only traversed
once, it is not eligible to be included in ``snake_axes``. As a convenience,
you may use ``snake_axes=True`` to enable snaking for all except that first
axis.

.. plot::

    from bluesky.simulators import plot_raster_path
    from ophyd.sim import motor1, motor2, det
    from bluesky.plans import grid_scan
    import matplotlib.pyplot as plt

    snaked = grid_scan([det], motor1, -5, 5, 10, motor2, -7, 7, 15, snake_axes=True)
    not_snaked = grid_scan([det], motor1, -5, 5, 10, motor2, -7, 7, 15)

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    plot_raster_path(snaked, 'motor1', 'motor2', probe_size=.3, ax=ax1)
    plot_raster_path(not_snaked, 'motor1', 'motor2', probe_size=.3, ax=ax2)
    ax1.set_title('True')
    ax2.set_title('False')
    ax1.set_xlim(-6, 6)
    ax2.set_xlim(-6, 6)

Demo:

.. ipython:: python
    :suppress:

    from bluesky.plans import grid_scan
    from ophyd.sim import motor1, motor2, det4
    dets = [det4]

.. ipython:: python

    RE(grid_scan(dets,
                 motor1, -1.5, 1.5, 3,  # scan motor1 from -1.5 to 1.5 in 3 steps
                 motor2, -0.1, 0.1, 5))  # scan motor2 from -0.1 to 0.1 in 5 steps

.. plot::

    from bluesky.plans import grid_scan
    from ophyd.sim import motor1, motor2, det4
    dets = [det4]
    RE(grid_scan(dets,
                 motor1, -1.5, 1.5, 3,  # scan motor1 from -1.5 to 1.5 in 3 steps
                 motor2, -0.1, 0.1, 5))  # scan motor2 from -0.1 to 0.1 in 5 steps

To move motors along arbitrary trajectories instead of equally-spaced points,
use :func:`~bluesky.plans.list_grid_scan` and
:func:`~bluesky.plans.rel_list_grid_scan`.

.. code-block:: python

    from bluesky.plans import list_grid_scan

    RE(list_grid_scan(dets,
                      motor1, [1, 1, 2, 3, 5],
                      motor2, [25, 16, 9]))

Demo:

.. ipython:: python
   :suppress:

   from bluesky.plans import list_grid_scan

.. ipython:: python

    RE(list_grid_scan(dets,
                      motor1, [1, 1, 2, 3, 5],
                      motor2, [25, 16, 9]))

.. plot::

    from bluesky.plans import list_grid_scan
    from ophyd.sim import det4, motor1, motor2
    dets = [det4]
    RE(list_grid_scan(dets,
                      motor1, [1, 1, 2, 3, 5],
                      motor2, [25, 16, 9]))

See :ref:`multi-dimensional_scans` to handle more specialized cases, including
combinations of :func:`~bluesky.plans.scan`-like and
:func:`~bluesky.plans.grid_scan`-like movement.

More generally, the :doc:`plans` documentation includes more exotic
trajectories, such as spirals, and plans with adaptive logic, such as
efficient peak-finders.

Aside: Access Saved Data
========================

At this point it is natural to wonder, "How do I access my saved data?"
From the point of view of *bluesky*, that's really not bluesky's concern, but
it's a reasonable question, so we'll address a typical scenario.

.. note::

    This section presumes that you are using the databroker. (We configured
    one in :ref:`an earlier section of this tutorial <databroker_setup>`.)
    You don't have to use the databroker to use bluesky; it's just
    one convenient way to capture the metadata and data generated by the
    RunEngine.

Very briefly, you can access saved data by referring to a dataset (a "run") by
its unique ID, which is returned by the RunEngine at collection time.

.. ipython:: python

    from bluesky.plans import count
    from ophyd.sim import det
    uid, = RE(count([det], num=3))
    header = db[uid]

Alternatively, perhaps more conveniently, you can access it by recency:

.. ipython:: python

    header = db[-1]  # meaning '1 run ago', i.e. the most recent run

.. note::

    We assumed above that the plan generated one "run" (dataset), which is
    typical for simple plans like :func:`~bluesky.plans.count`. In the
    *general* case, a plan can generate multiple runs, returning multiple uids,
    which in turn causes ``db`` to return a list of headers, not just one.

    .. code-block:: python

        uids = RE(some_plan(...))
        headers = db[uids]  # list of Headers

Most of the useful metadata is in this dictionary:

.. ipython:: python

    header.start

And the ("primary") stream of data is accessible like so:

.. ipython:: python

    header.table()  # return a table (a pandas.DataFrame)

From here we refer to the
`databroker tutorial <https://nsls-ii.github.io/databroker/tutorial.html>`_.

.. _tutorial_simple_customization:

Simple Customization
====================

Save Some Typing with 'Partial'
-------------------------------

Suppose we nearly always use the same detector(s) and we tire of typing out
``count([det])``. We can write a custom variant of :func:`~bluesky.plans.count`
using a built-in function provided by Python itself, :func:`functools.partial`.

.. code-block:: python

    from functools import partial
    from bluesky.plans import count
    from ophyd.sim import det

    my_count = partial(count, [det])
    RE(my_count())  # equivalent to RE(count([det]))

    # Additional arguments to my_count() are passed through to count().
    RE(my_count(num=3, delay=1))

Plans in Series
---------------

A custom plan can dispatch out to other plans using the Python syntax
``yield from``. (See :ref:`appendix <yield_from_primer>` if you want to know
why.) Examples:

.. code-block:: python

    from bluesky.plans import scan

    def coarse_and_fine(detectors, motor, start, stop):
        "Scan from 'start' to 'stop' in 10 steps and then again in 100 steps."
        yield from scan(detectors, motor, start, stop, 10)
        yield from scan(detectors, motor, start, stop, 100)

    RE(coarse_and_fine(dets, motor, -1, 1))

All of the plans introduced thus far, which we imported from
:mod:`bluesky.plans`, generate data sets ("runs"). Plans in the
:mod:`bluesky.plan_stubs` module do smaller operations. They can be used alone
or combined to build custom plans.

The :func:`~bluesky.plan_stubs.mv` plan moves one or more devices and waits for
them all to arrive.

.. code-block:: python

    from bluesky.plan_stubs import mv
    from ophyd.sim import motor1, motor2

    # Move motor1 to 1 and motor2 to 10, simultaneously. Wait for both to arrive.
    RE(mv(motor1, 1, motor2, 10))

We can combine :func:`~bluesky.plan_stubs.mv` and :func:`~bluesky.plans.count`
into one plan like so:

.. code-block:: python

    def move_then_count():
        "Move motor1 and motor2 into position; then count det."
        yield from mv(motor1, 1, motor2, 10)
        yield from count(dets)

    RE(move_then_count())

It's very important to remember the ``yield from``. The following plan does
nothing at all! (The plans inside it will be *defined* but never executed.)

.. code-block:: python

    # WRONG EXAMPLE!

    def oops():
        "Forgot 'yield from'!"
        mv(motor1, 1, motor2, 10)
        count(dets)

Much richer customization is possible, but we'll leave that for a
:ref:`a later section of this tutorial <tutorial_custom_plans>`. See also the
complete list of :ref:`plan stubs <stub_plans>`.

.. warning::

    **Never put ``RE(...)`` inside a loop or a function. You should always call
    it directly --- typed by the user at the terminal --- and only once.**

    You might be tempted to write a script like this:

    .. code-block:: python

        from bluesky.plans import scan
        from ophyd.sim import motor, det

        # Don't do this!
        for j in [1, 2, 3]:
            print(j, 'steps')
            RE(scan([det], motor, 5, 10, j)))

    Or a function like this:

    .. code-block:: python

        # Don't do this!
        def bad_function():
            for j in [1, 2, 3]:
                print(j, 'steps')
                RE(scan([det], motor, 5, 10, j)))

    But, instead, you should do this:

    .. code-block:: python

        from bluesky.plans import scan
        from ophyd.sim import motor, det

        def good_plan():
            for j in [1, 2, 3]:
                print(j, 'steps')
                yield from scan([det], motor, 5, 10, j)

        RE(my_plan())

    If you try to hide ``RE`` inside a function, someone later might
    use that function inside another function, and now we're entering and
    exiting the RunEngine multiple times from a single prompt. This can lead
    to unexpected behavior, especially around handling interruptions and
    errors.

    To indulge a musical metaphor, the plan is the sheet music, the hardware is
    the orchestra, and the RunEngine is the conductor. There should be only
    one conductor and she needs to run whole show, start to finish.

"Baseline" Readings (and other Supplemental Data)
=================================================

In addition to the detector(s) and motor(s) of primary interest during an
experiment, it is commonly useful to take a snapshot ("baseline reading") of
other hardware. This information is typically used to check consistency over
time. ("Is the temperature of the sample mount roughly the same as it was last
week?") Ideally, we'd like to *automatically* capture readings from these
devices during all future experiments without any extra thought or typing per
experiment. Bluesky provides a specific solution for this.

Configure
---------

.. note::

    If you are visiting user at a facility that runs bluesky, you may not need
    to do this configuration, and you can skip the next subsection just below
    --- :ref:`choose_baseline_devices`.

    You can type ``sd`` to check. If you get something like:

    .. ipython::
        :verbatim:

        In [1]: sd
        Out[1]: SupplementalData(baseline=[], monitors=[], flyers=[])

    you should skip this configuration.

Before we begin, we have to do a little more RunEngine configuration, like what
we did in the :ref:`tutorial_run_engine_setup` section with ``RE.subscribe``.

.. code-block:: python

    from bluesky.preprocessors import SupplementalData

    sd = SupplementalData()
    RE.preprocessors.append(sd)

.. ipython:: python
    :suppress:

    from bluesky.preprocessors import SupplementalData
    sd = SupplementalData()
    RE.preprocessors.append(sd)

.. _choose_baseline_devices:

Choose "Baseline" Devices
-------------------------

We'll choose the detectors/motors that we want to be read automatically at the
beginning and end of each dataset ("run"). If you are using a shared
configuration, this also might already have been done, so you should check the
content of ``sd.baseline`` before altering it.

.. ipython:: python

    sd.baseline  # currently empty

Suppose that we want to take baseline readings from three detectors and two
motors. We'll import a handful of simulated devices for this purpose, put them
into a list, and assign ``sd.baseline``.

.. ipython:: python

    from ophyd.sim import det1, det2, det3, motor1, motor2
    sd.baseline = [det1, det2, det3, motor1, motor2]

Notice that we can put a mixture of detectors and motors in this list. It
doesn't matter to bluesky that some are movable and some are not because it's
just going to be *reading* them, and both detectors and motors can be read.

Use
---

Now we can just do a scan with the detector and motor of primary interest. The
RunEngine will automatically take baseline readings before and after each run.
Demo:

.. ipython:: python

    from ophyd.sim import det, motor
    from bluesky.plans import scan
    RE(scan([det], motor, -1, 1, 5))

We can clear or update the list of baseline detectors at any time.

.. ipython:: python

    sd.baseline = []

As an aside, this is one place where the design of bluesky really pays off. By
separating the executor (the RunEngine) from the instruction sets (the plans)
it's easy to apply global configuration without updating every plan
individually.

Access Baseline Data
--------------------

If you access the data from our baseline scan, you might think that the
baseline data is missing!

.. ipython:: python

    header = db[-1]
    header.table()

Looking again at the output when we executed this scan, notice these lines:

.. code-block:: none

    New stream: 'baseline'
    ...
    New stream: 'primary'

By default, ``header.table()`` gives us the "primary" data stream:

.. ipython:: python

    header.table('primary')  # same result as header.table()

We can access other streams by name.

.. ipython:: python

    header.table('baseline')

A list of the stream names in a given run is available as
``header.stream_names``. From here we refer to the
`databroker tutorial <https://nsls-ii.github.io/databroker/tutorial.html>`_.

Other Supplemental Data
-----------------------

Above, we used ``sd.baseline``. There is also ``sd.monitors`` for signals to
monitor asynchronously during a run and ``sd.flyers`` for devices to "fly-scan"
during a run. See :ref:`supplemental_data` for details.

.. _tutorial_pause_resume_suspend:

Pause, Resume, Suspend
======================

Interactive Pause & Resume
--------------------------

Sometimes it is convenient to pause data collection, check on some things, and
then either resume from where you left off or quit. The RunEngine makes it
possible to do this cleanly and safely on *any* plan, including user-defined
plans, with minimal effort by the user. Of course, experiments on systems
that evolve with time can't be arbitrarily paused and resumed. It's up to the
user to know that and use this feature only when applicable.

Take this example, a step scan over ten points.

.. code-block:: python

    from ophyd.sim import det, motor
    from bluesky.plans import scan

    motor.delay = 1  # simulate slow motor movement
    RE(scan([det], motor, 1, 10, 10))

Demo:

.. ipython::
    :verbatim:

    In [1]: RE(scan([det], motor, 1, 10, 10))
    Transient Scan ID: 1     Time: 2018/02/12 12:40:36
    Persistent Unique Scan ID: 'c5db9bb4-fb7f-49f4-948b-72fb716d1f67'
    New stream: 'primary'
    +-----------+------------+------------+------------+
    |   seq_num |       time |      motor |        det |
    +-----------+------------+------------+------------+
    |         1 | 12:40:37.6 |      1.000 |      0.607 |
    |         2 | 12:40:38.7 |      2.000 |      0.135 |
    |         3 | 12:40:39.7 |      3.000 |      0.011 |

At this point we decide to hit **Ctrl+C** (SIGINT). The RunEngine will catch
this signal and react like so. We will examine this output piece by piece.

.. code-block:: none

    ^C
    A 'deferred pause' has been requested.The RunEngine will pause at the next
    checkpoint. To pause immediately, hit Ctrl+C again in the next 10 seconds.
    Deferred pause acknowledged. Continuing to checkpoint.
    <...a few seconds later...>
    |         4 | 12:40:40.7 |      4.000 |      0.000 |
    Pausing...

    ---------------------------------------------------------------------------
    RunEngineInterrupted                      Traceback (most recent call last)
    <ipython-input-14-826ee9dfb918> in <module>()
    ----> 1 RE(scan([det], motor, 1, 10, 10))
    <...snipped details...>

    RunEngineInterrupted:
    Your RunEngine is entering a paused state. These are your options for changing
    the state of the RunEngine:
    RE.resume()    Resume the plan.
    RE.abort()     Perform cleanup, then kill plan. Mark exit_stats='aborted'.
    RE.stop()      Perform cleanup, then kill plan. Mark exit_status='success'.
    RE.halt()      Emergency Stop: Do not perform cleanup --- just stop.

When it pauses, the RunEngine immediately tells all Devices that it has touched
so far to "stop". (Devices define what that means to them in their ``stop()``
method.) This is not a replacement for proper equipment protection; it is just
a convenience.

Now, at our leisure, we may:

* pause to think
* investigate the state of our hardware, such as the detector's exposure time
* turn on more verbose logging  (see :doc:`debugging`)
* decide whether to stop here or resume

Suppose we decide to resume. The RunEngine will pick up from the last
"checkpoint". Typically, this means beginning of each step in a scan, but
plans may specify checkpoints anywhere they like.

.. ipython::
    :verbatim:

    In [13]: RE.resume()
    |         5 | 12:40:50.1 |      5.000 |      0.000 |
    |         6 | 12:40:51.1 |      6.000 |      0.000 |
    |         7 | 12:40:52.1 |      7.000 |      0.000 |
    |         8 | 12:40:53.1 |      8.000 |      0.000 |
    |         9 | 12:40:54.1 |      9.000 |      0.000 |
    |        10 | 12:40:55.1 |     10.000 |      0.000 |
    +-----------+------------+------------+------------+
    generator scan ['c5db9bb4'] (scan num: 1)

The scan has completed successfully.

If you go back and read the output from when we hit Ctrl+C, you will notice
that the RunEngine didn't pause immediately: it finished the current step of
the scan first. Quoting an excerpt from the demo above:

.. code-block:: none

    ^C
    A 'deferred pause' has been requested.The RunEngine will pause at the next
    checkpoint. To pause immediately, hit Ctrl+C again in the next 10 seconds.
    Deferred pause acknowledged. Continuing to checkpoint.
    <...a few seconds later...>
    |         4 | 12:40:40.7 |      4.000 |      0.000 |
    Pausing...

Observe that hitting Ctrl+C *twice* pauses immediately, without waiting to
finish the current step.

.. code-block:: none

    In [2]: RE(scan([det], motor, 1, 10, 10))
    Transient Scan ID: 2     Time: 2018/02/15 12:31:14
    Persistent Unique Scan ID: 'b342448f-6a64-4f26-91a6-37f559cb5537'
    New stream: 'primary'
    +-----------+------------+------------+------------+
    |   seq_num |       time |      motor |        det |
    +-----------+------------+------------+------------+
    |         1 | 12:31:15.8 |      1.000 |      0.607 |
    |         2 | 12:31:16.8 |      2.000 |      0.135 |
    |         3 | 12:31:17.8 |      3.000 |      0.011 |
    ^C^C
    Pausing...

When resumed, the RunEngine will *rewind* to the last checkpoint (the beginning
of the fourth step in the scan) and repeat instructions as needed.

Quoting again from the demo, notice that ``RE.resume()`` was only one of our
options. If we decide not to continue we can quit in three different ways:

.. code-block:: none

    Your RunEngine is entering a paused state. These are your options for changing
    the state of the RunEngine:
    RE.resume()    Resume the plan.
    RE.abort()     Perform cleanup, then kill plan. Mark exit_stats='aborted'.
    RE.stop()      Perform cleanup, then kill plan. Mark exit_status='success'.
    RE.halt()      Emergency Stop: Do not perform cleanup --- just stop.

"Aborting" and "stopping" are almost the same thing: they just record different
metadata about why the experiment was ended. Both signal to the plan that it
should end early, but they still let it specify more instructions so that it
can "clean up." For example, a :func:`~bluesky.plans.rel_scan` moves the motor
back to its starting position before quitting.

In rare cases, if we are worried that the plan's cleanup procedure might be
dangerous, we can "halt". Halting circumvents the cleanup instructions.

Try executing ``RE(scan([det], motor, 1, 10, 10))``, pausing, and exiting in
these various ways. Observe that the RunEngine won't let you run a new plan
until you have resolved the paused plan using one of these methods.

Automated Suspend & Resume
--------------------------

The RunEngine can be configured in advance to *automatically* pause and resume
in response to external signals. To distinguish automatic pause/resume from
interactive, user-initiated pause and resume, we call this behavior
"suspending."

For details, see :ref:`suspenders`.

.. _tutorial_metadata:

Metadata
========

If users pass extra keyword arguments to ``RE``, they are interpreted as
metadata

.. code-block:: python

    RE(count([det]), user='Dan', mood='skeptical')
    RE(count([det]), user='Dan', mood='optimistic')

and they can be used for searching later:

.. code-block:: python

    headers = db(user='Dan')
    headers = db(mood='skeptical')

Metadata can also be added *persistently* (i.e. applied to all future runs
until removed) by editing the dictionary ``RE.md``.

.. code-block:: python

    RE.md
    RE.md['user'] = 'Dan'

No need to specify ``user`` every time now....

.. code-block:: python

    RE(count([det]))  # automatically includes user='Dan'

The key can be temporarily overridden:

.. code-block:: python

    RE(count([det]), user='Tom')  # overrides the setting in RE.md, just once

or deleted:

.. code-block:: python

    del RE.md['user']

In addition to any user-provided metadata, the RunEngine, the devices, and the
plan capture some metadata automatically. For more see, :doc:`metadata`.

Simulate and Introspect Plans
=============================

We have referred to a *plan* as a "sequence of instructions encoding an
experimental procedure." But what's inside a plan really? Bluesky calls each
atomic instruction inside a plan a *message*.  Handling the messages directly
is only necessary when debugging or doing unusually deep customization, but
it's helpful to see them at least once before moving on to more practical
tools.

Try printing out every message in a couple simple plans:

.. code-block:: python

    from bluesky.plans import count
    from ophyd.sim import det

    for msg in count([]):
        print(msg)

    for msg in count([det]):
        print(msg)

See the :doc:`msg` section for more.

Bluesky includes some tools for producing more useful, human-readable summaries
to answer the question, "What will this plan do?"

.. ipython:: python

    from bluesky.simulators import summarize_plan
    from bluesky.plans import count, rel_scan
    from ophyd.sim import det, motor
    # Count a detector 3 times.
    summarize_plan(count([det], 3))
    # A 3-step scan.
    summarize_plan(rel_scan([det], motor, -1, 1, 3))

For more possibilities, see :doc:`simulation`.

.. _tutorial_device:

Devices
=======

Theory
------

The notion of a "Device" serves two goals:

* Provide a **standard interface** to all hardware for the sake of generality
  and code reuse.
* **Logically group** individual signals into composite "Devices" that can be
  read together, as a unit, and configured in a coordinated way. Provide a
  human-readable name to this group, with an eye toward later data analysis.

In bluesky's view of the world, there are only three different kinds of devices
used in data acquisition.

* Some devices can be **read**. This includes simple points detectors that
  produce a single number and large CCD detectors that produce big arrays.
* Some devices can be both **read and set**. Setting a motor physically moves
  it to a new position. Setting a temperature controller impels it to gradually
  change its temperature. Setting the exposure time on some detector promptly
  updates its configuration.
* Some devices produce data at a rate too high to be read out in real time, and
  instead **buffer their data externally** in separate hardware or software
  until it can be read out.

Bluesky interacts with all devices via a :doc:`specified interface <hardware>`.
Each device is represented by a Python object with certain methods and
attributes (with names like ``read`` and ``set``). Some of these methods are
asynchronous, such as ``set``, which allows for the concurrent movement of
multiple devices.

Implementation
--------------

`Ophyd <https://nsls-ii.github.io/ophyd>`_, a Python library that was
developed in tandem with bluesky, implements this interface for devices that
speak `EPICS <http://www.aps.anl.gov/epics/>`_. But bluesky is not tied to
ophyd or EPICS specifically: any Python object may be used, so long as it
provides the specified methods and attributes that bluesky expects.
For example, an experimental implementation of the bluesky interface for LabView has been written.
See :ref:`Hardware Interface Packages <hardware_interface_packages>` for more examples.
And the simulated hardware that we have been using in this
tutorial is all based on pure-Python constructs unconnected from hardware or
any specific hardware control protocol.

To get a flavor for what it looks like to configure hardware in ophyd,
connecting to an EPICS motor looks like this:

.. code-block:: python

    from ophyd import EpicsMotor

    nano_top_x = EpicsMotor('XF:31ID-ES{Dif:Nano-Ax:TopX}Mtr', name='nano_top_x')

We have provided both the machine-readable address of the motor on the network,
``'XF:31ID-ES{Dif:Nano-Ax:TopX}Mtr'`` (in EPICS jargon, the "PV" for
"Process Variable"), and a human-readable name, ``'nano_top_x'``, which will be
used to label the data generated by this motor. When it comes time to analyze
the data, we will be grateful to be dealing with the human-readable label.

The ``EpicsMotor`` device is a logical grouping of many signals. The most
important are the readback (actual position) and setpoint (target position).
All of the signals are summarized thus. The details here aren't important at
this stage: the take-away message is, "There is a lot of stuff to keep track of
about a motor, and a Device helpfully groups that stuff for us."

.. code-block:: none

    In [3]: nano_top_x.summary()
    data keys (* hints)
    -------------------
    *nano_top_x
    nano_top_x_user_setpoint

    read attrs
    ----------
    user_readback        EpicsSignalRO       ('nano_top_x')
    user_setpoint        EpicsSignal         ('nano_top_x_user_setpoint')

    config keys
    -----------
    nano_top_x_acceleration
    nano_top_x_motor_egu
    nano_top_x_user_offset
    nano_top_x_user_offset_dir
    nano_top_x_velocity

    configuration attrs
    ----------
    motor_egu            EpicsSignal         ('nano_top_x_motor_egu')
    velocity             EpicsSignal         ('nano_top_x_velocity')
    acceleration         EpicsSignal         ('nano_top_x_acceleration')
    user_offset          EpicsSignal         ('nano_top_x_user_offset')
    user_offset_dir      EpicsSignal         ('nano_top_x_user_offset_dir')

    Unused attrs
    ------------
    offset_freeze_switch EpicsSignal         ('nano_top_x_offset_freeze_switch')
    set_use_switch       EpicsSignal         ('nano_top_x_set_use_switch')
    motor_is_moving      EpicsSignalRO       ('nano_top_x_motor_is_moving')
    motor_done_move      EpicsSignalRO       ('nano_top_x_motor_done_move')
    high_limit_switch    EpicsSignal         ('nano_top_x_high_limit_switch')
    low_limit_switch     EpicsSignal         ('nano_top_x_low_limit_switch')
    direction_of_travel  EpicsSignal         ('nano_top_x_direction_of_travel')
    motor_stop           EpicsSignal         ('nano_top_x_motor_stop')
    home_forward         EpicsSignal         ('nano_top_x_home_forward')
    home_reverse         EpicsSignal         ('nano_top_x_home_reverse')


.. _tutorial_custom_plans:

Write Custom Plans
==================

As mentioned in the :ref:`tutorial_simple_customization` section above, the
"pre-assembled" plans with :func:`~bluesky.plans.count` and
:func:`~bluesky.plans.scan` are built from smaller "plan stubs". We can
mix and match the "stubs" and/or "pre-assembled" plans to build custom plans.

There are many of plan stubs, so it's convenient to import the whole module and
work with that.

.. code-block:: python

    import bluesky.plan_stubs as bps

Move in Parallel
----------------

Before writing a custom plan to coordinate the motion of multiple devices,
consider whether your use case could be addressed with one of the built-in
:ref:`multi-dimensional_scans`.

We previously introduced the :func:`~bluesky.plan_stubs.mv` plan that moves one
or more devices and waits for them all to arrive. There is also
:func:`~bluesky.plans.mvr` for moving *relative* to the current position.

.. code-block:: python

    from ophyd.sim import motor1, motor2

    # Move motor1 to 1 and motor2 10 units in the positive direction relative
    # to their current positions. Wait for both to arrive.
    RE(bps.mvr(motor1, 1, motor2, 10))

Some scenarios require more low-level control over when the waiting occurs.
For these, we employ :func:`~bluesky.plan_stubs.wait` and
:func:`~bluesky.plan_stubs.abs_set` ("absolute set") or
:func:`~bluesky.plan_stubs.rel_set` ("relative set").

Here is a scenario that does require a custom solution: we want to set several
motors in motion at once, including multiple fast motors and one slow motor. We
want to wait for the fast motors to arrive, print a message, then wait for the
slow motor to arrive, and print a second message.

.. code-block:: python

    def staggered_wait(fast_motors, slow_motor):
        # Start all the motors, fast and slow, moving at once.
        # Put all the fast_motors in one group...
        for motor in fast_motors:
            yield from bps.abs_set(motor, 5, group='A')
        # ...but put the slow motor is separate group.
        yield from bps.abs_set(slow_motor, 5, group='B')

        # Wait for all the fast motors.
        print('Waiting on the fast motors.')
        yield from bps.wait('A')
        print('Fast motors are in place. Just waiting on the slow one now.')

        # Then wait for the slow motor.
        yield from bps.wait('B')
        print('Slow motor is in place.')

Sleeping (Timed Delays)
-----------------------

.. note::

    If you need to wait for your motor to finish moving, temperature to finish
    equilibrating, or shutter to finish opening, inserting delays into plans
    isn't the best way to do that. It should be the *Device's* business to
    report accurately when it is done, including any extra padding for settling
    or equilibration. On some devices, such as ``EpicsMotor``, this can be
    configured like ``motor.settle_time = 3``.

For timed delays, bluesky has a special plan, which allows the RunEngine to
continue its business during the sleep.

.. code-block:: python

    def sleepy_plan(motor, positions):
        "Step a motor through a list of positions with 1-second delays between steps.")
        for position in positions:
            yield from bps.mv(motor, position)
            yield from bps.sleep(1)

**You should always use this plan, *never* Python's built-in function
:func:`time.sleep`.** Why?
The RunEngine uses an event loop to concurrently manage many tasks. It assumes
that none of those tasks blocks for very long. (A good figure for "very long"
is 0.2 seconds.) Therefore, you should never incorporate long blocking function
calls in your plan, such as ``time.sleep(1)``.

.. _tutorial_capture_data:

Capture Data
------------

.. ipython:: python
    :suppress:

    # Define a examples that we will use interactively below.
    import bluesky.plan_stubs as bps
    def one_run_one_event(detectors):
        yield from bps.open_run()
        yield from bps.trigger_and_read(detectors)
        yield from bps.close_run()
    def one_run_multi_events(detectors, num):
        yield from bps.open_run()
        for i in range(num):
            yield from bps.trigger_and_read(detectors)
        yield from bps.close_run()
    def multi_runs_multi_events(detectors, num, num_runs):
        for i in range(num_runs):
            yield from one_run_multi_events(detectors, num)

Any plan that generates data must include instructions for grouping readings
into *Events* (i.e. rows in a table) and grouping those Events into *Runs*
(datasets that are given a "scan ID"). This is best explained by example.

.. code-block:: python

    import bluesky.plan_stubs as bps

    def one_run_one_event(detectors):
        # Declare the beginning of a new run.
        yield from bps.open_run()

        # Trigger each detector and wait for triggering to complete.
        # Then read the detectors and bundle these readings into an Event
        # (i.e. one row in a table.)
        yield from bps.trigger_and_read(detectors)

        # Declare the end of the run.
        yield from bps.close_run()

Execute the plan like so:

.. ipython:: python

    RE(one_run_one_event([det1, det2]))

We observe:

* one table (one Run)
* one row (one Event)
* two columns (a column for each detector)

Here's the same plan again, with :func:`~bluesky.plan_stubs.trigger_and_read`
moved inside a for loop.

.. code-block:: python

    def one_run_multi_events(detectors, num):
        yield from bps.open_run()

        for i in range(num):
            yield from bps.trigger_and_read(detectors)

        yield from bps.close_run()

Execute the plan like so:

.. ipython:: python

    RE(one_run_multi_events([det1, det2], 3))

We observe:

* one table (one Run)
* three rows (three Events)
* two columns (a column for each detector)

Finally, add another loop re-using ``one_run_multi_events`` inside that loop.

.. code-block:: python

    def multi_runs_multi_events(detectors, num, num_runs):
        for i in range(num_runs):
            yield from one_run_multi_events(detectors, num)

.. ipython:: python

    RE(multi_runs_multi_events([det1, det2], num=3, num_runs=2))

We observe:

* two tables (two Runs)
* three rows (three Events)
* two columns (a column for each detector)

We also notice that the return value output from the RunEngine is a tuple with
two unique IDs, one per Run generated by this plan.

In order to focus on the scope of an Event and a Run, we have left out an
important detail, addressed in the next section, which may be necessary to
incorporate before trying these plans on real devices.

Stage and Unstage
-----------------

Complex devices often require some preliminary setup before they can be used
for data collection, moving them from a resting state into a state where they
are ready to acquire data. Bluesky accommodates this in a general way by
allowing every Device to implement an optional ``stage()`` method, with a
corresponding ``unstage()`` method. Plans should stage every device that they
touch exactly once and unstage every device at the end. If a Device does not
have a ``stage()`` method the RunEngine will just skip over it.

Revising our simplest example above, ``one_run_one_event``,

.. code-block:: python

    import bluesky.plan_stubs as bps

    def one_run_one_event(detectors):
        yield from bps.open_run()
        yield from bps.trigger_and_read(detectors)
        yield from bps.close_run()

we incorporate staging like so:

.. code-block:: python

    def one_run_one_event(detectors):

        # 'Stage' every device.
        for det in detectors:
            yield from bps.stage(det)

        yield from bps.open_run()
        yield from bps.trigger_and_read(detectors)
        yield from bps.close_run()

        # 'Unstage' every device.
        for det in detectors:
            yield from bps.unstage(det)

This is starting to get verbose. At this point, we might want to accept some
additional complexity in exchange for brevity --- and some assurance that we
don't forget to use these plans in matching pairs. To that end, this plan is
equivalent:

.. code-block:: python

    import bluesky.preprocessors as bpp

    def one_run_one_event(detectors):

        @bpp.stage_decorator(detectors)
        def inner():
            yield from bps.open_run()
            yield from bps.trigger_and_read(detectors)
            yield from bps.close_run()

        return (yield from inner())

The :func:`~bluesky.preprocessors.stage_decorator` is a *plan preprocessor*, a
plan which consumes another plan and modifies its instructions. In this case,
it adds inserts 'stage' and 'unstage' messages, supplanting
:func:`~bluesky.plan_stubs.stage` and :func:`~bluesky.plan_stubs.unstage`. We
can trim the verbosity down yet more by employing
:func:`~bluesky.preprocessors.run_decorator`, supplanting
:func:`~bluesky.plan_stubs.open_run` and :func:`~bluesky.plan_stubs.close_run`.
The result:

.. code-block:: python

    import bluesky.preprocessors as bpp

    def one_run_one_event(detectors):

        @bpp.stage_decorator(detectors)
        @bpp.run_decorator()
        def inner():
            yield from bps.trigger_and_read(detectors)

        return (yield from inner())

Incidentally, recall that we have already encountered a preprocessor in this
tutorial, in the section on baseline readings.
:class:`~bluesky.preprocessors.SupplementalData` is a preprocessor.

.. _tutorial_plan_metadata:

Add Metadata
------------

To make it easier to search for data generated by the plan and to inspect what
was done afterward, we should include some metadata. We create a dictionary and
pass it to :func:`~bluesky.preprocessors.run_decorator` (or, in the more
verbose formulation, to :func:`~bluesky.plan_stubs.open_run`). The RunEngine
will combine this metadata with any information provided by the user, as shown
in the :ref:`the earlier section on metadata <tutorial_metadata>`.

.. code-block:: python

    def one_run_one_event(detectors):

        md = {
            # Human-friendly names of detector Devices (useful for searching)
            'detectors': [det.name for det in detectors],

            # The Python 'repr's each argument to the plan
            'plan_args': {'detectors': list(map(repr, detectors))},

            # The name of this plan
            'plan_name': 'one_run_one_event',
        }

        @bpp.stage_decorator(detectors)
        @bpp.run_decorator(md)
        def inner():
            yield from bps.trigger_and_read(detectors)

        return (yield from inner())

.. warning::

    The values in the metadata dictionary must be strings, numbers,
    lists/arrays, or dictionaries only. Metadata cannot contain arbitrary
    Python types because downstream consumers (like databases) do not know what
    to do with those and will error.

To be polite, we should allow the user to override this metadata. All of
bluesky's "pre-assembled" plans (:func:`~bluesky.plans.count`,
:func:`~bluesky.plans.scan`, etc.) provide an optional ``md`` argument for this
purpose, implemented like so:

.. code-block:: python

    def one_run_one_event(detectors, md=None):

        _md = {
            'detectors': [det.name for det in detectors],
            'plan_args': {'detectors': list(map(repr, detectors))},
            'plan_name': 'one_run_one_event',
        }

        # If a key exists in md, it overwrites the default in _md.
        _md.update(md or {})

        @bpp.stage_decorator(detectors)
        @bpp.run_decorator(_md)
        def inner():
            yield from bps.trigger_and_read(detectors)

        return (yield from inner())

Add "Hints" in Metadata
-----------------------

The metadata dictionary may optionally include a key named ``'hints'``. This
key has special significance to the
:class:`~bluesky.callback.best_effort.BestEffortCallback` and potentially
other downstream consumers, which use it to try to infer useful ways to
present the data. Currently, it solves two specific problems.

1. Narrow the potentially large set of readings to a manageable number of most
   important ones that fit into a table.
2. Identify the dimensionality of the data (1D scan? 2D grid? N-D grid?) and
   the dependent and independent parameters, for visualization and peak-fitting
   purposes.

It's up to each device to address (1). The plan has no role in that.
Each device has an optional ``hints`` attribute with a value like
``{'fields': [...]}`` to answer the question, "Of all the readings you
produce, what are the names of the most important ones?"

We need the plan to help us with (2). Only the plan can sort out which devices
are being employed as "independent" axes and which are being measured as
dependent variables. This isn't clear just from looking at the Devices alone
because any given movable device can be used as an axis or as a "detector"
depending on the context --- ``count([motor])`` is a perfectly valid thing to
do!

The schema of the plan's hint metadata is:

.. code-block:: python

    {'dimensions': [([<FIELD>, ...], <STREAM_NAME>),
                    ([<FIELD>, ...], <STREAM_NAME>),
                    ...
                   ]}

Examples:

.. code-block:: python

    # a 1-D scan over x
    {'dimensions': [(['x'], 'primary')]}

    # a 2-D grid_scan over x and y
    {'dimensions': [(['x'], 'primary'),
                    (['y'], 'primary')]}

    # a scan moving x and y together along a diagonal
    {'dimensions': [(['x', 'y'], 'primary')]}

    # a 1-D scan over temperature, represented in C and K units
    {'dimensions': [(['C', 'K'], 'primary')]}

    # a 1-D scan over energy, as measured in energy and diffractometer position
    {'dimensions': [(['E', 'dcm'], 'primary')]}

    # special case: a sequence of readings where the independent axis is just time
    {'dimensions': [(['time'], 'primary')]}

Each entry in the outer list represents one independent dimension. A dimension
might be represented by multiple fields, either from different devices moved in
a coordinated fashion by the plan (``['x', 'y']``), presented as fully redundant
information from one device (``['C', 'K']``), or coupled information from two
sub-devices (``['E', 'dcm']``).

The second element in each entry is the stream name: ``'primary'`` in every
example above.  This should correspond to the ``name`` passed into
:func:`~bluesky.plan_stubs.trigger_and_read` or
:func:`~bluesky.plan_stubs.create` inside the plan. The default name is
``primary``.

Putting it all together, the plan asks the device(s) being used as independent
axes for their important field(s) and builds a list of dimensions like so:

.. code-block:: python

   dimensions = [(motor.hints['fields'], 'primary')]

We must account for the fact that ``hints`` is optional. A given Device
might not have a ``hints`` attribute at all and, even if it does, the
hints might not contain the ``'fields'`` key that we are interested in. This
pattern silently omits the dimensions hint if the necessary information is not
provided by the Device:

.. code-block:: python

    def scan(..., md=None):
        _md = {...}
        _md.update(md or {})

        try:
            dimensions = [(motor.hints['fields'], 'primary')]
        except (AttributeError, KeyError):
            pass
        else:
            _md['hints'].setdefault('dimensions', dimensions)

        ...

Finally, by using ``setdefault``, we have allowed user to override these hints
if they know better by passing in ``scan(..., md={'hints': ...})``.

.. _tutorial_adaptive:

Adaptive Logic in a Plan
------------------------

Two-way communication is possible between the generator and the RunEngine.
For example, the :func:`~trigger_and_read` plan responds with its readings. We
can use it to make an on-the-fly decision about whether to continue or stop.

.. code-block:: python

    import bluesky.preprocessors as bpp
    import bluesky.plan_stubs as bps
    from ophyd.sim import det, motor
    def conditional_break(threshold):
        """Set, trigger, read until the detector reads intensity < threshold"""

        @bpp.stage_decorator([det, motor])
        @bpp.run_decorator()
        def inner():
            i = 0
            while True:
                yield from bps.mv(motor, i)
                readings = yield from bps.trigger_and_read([det])
                if readings['det']['value'] < threshold:
                    break
                i += 1
        return (yield from inner())

.. ipython:: python
    :suppress:

    import bluesky.preprocessors as bpp
    import bluesky.plan_stubs as bps
    from bluesky import Msg
    from ophyd.sim import det, motor
    def conditional_break(threshold):
        def inner():
            i = 0
            while True:
                yield from bps.mv(motor, i)
                readings = yield from bps.trigger_and_read([det])
                if readings['det']['value'] < threshold:
                    break
                i += 1
        # Decorators do not work in IPython sphinx directive!
        # Using wrapper instead...
        return (yield from bpp.stage_wrapper(bpp.run_wrapper(inner()), [det, motor]))

Demo:

.. ipython:: python

    RE(conditional_break(0.2))

The important line in this example is

.. code-block:: python

    reading = yield from bps.trigger_and_read([det])

The action proceeds like this:

1. The plan yields a 'read' message to the RunEngine.
2. The RunEngine reads the detector.
3. The RunEngine sends that reading *back to the plan*, and that response is
   assigned to the variable ``reading``.

The response, ``reading``, is formatted like:

.. code-block:: python

     {<name>: {'value': <value>, 'timestamp': <timestamp>}, ...}

For a detailed technical description of the messages and their responses,
see :ref:`msg`.

.. _tutorial_exception_handling:

Plan "Cleanup" (Exception Handling)
-----------------------------------

If an exception is raised, the RunEngine gives the plan the opportunity to
catch the exception and either handle it or merely yield some "clean up"
messages before re-raising the exception and killing plan execution. (Recall
this from :ref:`tutorial_pause_resume_suspend` above.)

This is the general idea:

.. code-block:: python

    # This example is illustrative, but it is not completely correct.
    # Use `finalize_wrapper` instead (or read its source code).

    def plan_with_cleanup():
        def main_plan():
            # do stuff...

        def cleanup_plan():
            # do other stuff...

        try:
            yield from main_plan()
        finally:
            # Do this even if an Exception is raised.
            yield from cleanup_plan()

The exception in question may originate from the plan itself or from the
RunEngine when it attempts to execute a given command.

The :func:`~bluesky.preprocessors.finalize_wrapper` preprocessor provides a
succinct and fully correct way of applying this general pattern.

.. code-block:: python

    import bluesky.preprocessors as bpp

    def plan_with_cleanup():
        yield from bpp.finalize_wrapper(main_plan(), cleanup_plan())

Further Reading
---------------

* :ref:`per_step_hook`
* Specifying checkpoints (TODO)
* Monitoring (TODO)
* Fly Scanning (TODO)
* :ref:`Pausing from a plan <planned_pauses>`
* :func:`~bluesky.plans.input_plan` (TODO)
* Going deeper than :func:`~bluesky.plan_stubs.trigger_and_read` (TODO)
Bluesky Data Collection Framework
=================================

Bluesky is a library for experiment control and collection of scientific data
and metadata. It emphasizes the following virtues:

* **Live, Streaming Data:** Available for inline visualization and processing.
* **Rich Metadata:** Captured and organized to facilitate reproducibility and
  searchability.
* **Experiment Generality:** Seamlessly reuse a procedure on completely
  different hardware.
* **Interruption Recovery:** Experiments are "rewindable," recovering cleanly
  from interruptions.
* **Automated Suspend/Resume:** Experiments can be run unattended,
  automatically suspending and resuming if needed.
* **Pluggable I/O:** Export data (live) into any desired format or database.
* **Customizability:** Integrate custom experimental procedures and commands,
  and get the I/O and interruption features for free.
* **Integration with Scientific Python:** Interface naturally with numpy and
  Python scientific stack.

How to Use This Documentation
-----------------------------

Start with the :doc:`tutorial`. It's a good place to start for everyone, and it
gives a good overview of the project in a narrative style. Read as far as you
need to solve your problem, and come back again if your needs change. Each
section of the tutorial adds a piece of complexity in exchange for deeper
customization.

The remaining sections document bluesky's behavior in a less narrative style,
providing clear API documentation intermixed with some examples and explanation
of design and intent.

Index
-----

.. toctree::
   :caption: User Documentation
   :maxdepth: 1

   tutorial
   plans
   documents
   metadata
   callbacks
   state-machine
   simulation
   progress-bar
   event_descriptors
   async
   multi_run_plans
   debugging
   run_engine_api
   utils
   magics
   from-pyepics-to-bluesky
   comparison-with-spec
   hardware-interfaces
   appendix

.. toctree::
   :caption: Developer Documentation
   :maxdepth: 1

   hardware
   msg
   run_engine
   api_changes
   contributing

.. toctree::
   :hidden:
   :caption: Data Collection

   bluesky <https://blueskyproject.io/bluesky>
   ophyd <https://blueskyproject.io/ophyd>

.. toctree::
   :hidden:
   :caption: Data Access and Management

   databroker <https://blueskyproject.io/databroker>
   amostra <https://nsls-ii.github.io/amostra>
   datamuxer <https://nsls-ii.github.io/datamuxer>
   suitcase <https://blueskyproject.io/suitcase>

.. toctree::
   :hidden:
   :caption: GitHub Links

   NSLS-II Repositories <https://github.com/NSLS-II/>
   Bug Reports <https://github.com/NSLS-II/Bug-Reports/issues>
*********************
Debugging and Logging
*********************

.. versionchanged:: 1.6.0

   Bluesky's use of Python's logging framework has been completely reworked to
   follow Python's documented best practices for libraries.

Bluesky uses Python's logging framework, which enables sophisticated log
management. For common simple cases, including viewing logs in the terminal or
writing them to a file, the next section illustrates streamlined,
copy/paste-able examples. Users who are familiar with that framework or who
need to route logs to multiple destinations may wish to skip ahead to
:ref:`logger_api`.

Useful Snippets
===============

Log warnings
------------

This is the recommended standard setup.

.. code-block:: python

   from bluesky import config_bluesky_logging
   config_bluesky_logging()

It will display ``'bluesky'`` log records of ``WARNING`` level or higher in the
terminal (standard out) with a format tailored to bluesky.

Maximum verbosity
-----------------

If the RunEngine is "hanging," running slowly, or repeatedly encountering an
error, it is useful to know exactly where in the plan the problem is occurring.
To follow the RunEngine's progress through the plan, crank up the verbosity of
the logging.

This will display each message from the plan just before the RunEngine
processes it, giving a clear indication of when plan execution is stuck.

.. code-block:: python

   from bluesky import config_bluesky_logging
   config_bluesky_logging(level='DEBUG')

Log to a file
-------------

This will direct all log messages to a file instead of the terminal (standard
out).

.. code-block:: python

    from bluesky import config_bluesky_logging
    config_bluesky_logging(file='/tmp/bluesky.log', level='DEBUG')

.. important::

   We strongly recommend setting levels on *handlers* not on *loggers*.
   In previous versions of bluesky, we recommended adjusting the level on the
   *logger*, as in ``RE.log.setLevel('DEBUG')``. We now recommended
   that you *avoid* setting levels on loggers because it would affect all
   handlers downstream, potentially inhibiting some other part of the program
   from collecting the records it wants to collect.

.. _logger_api:

Bluesky's Logging-Related API
=============================

Logger Names
------------

Here are the primary loggers used by bluesky.

* ``'bluesky'`` --- the logger to which all bluesky log records propagate
* ``'bluesky.emit_document'`` --- A log record is emitted whenever a Document
  is emitted. The log record does not contain the full content of the
  Document.
* ``'bluesky.RE'`` --- Records from a RunEngine. INFO-level notes state
  changes. DEBUG-level notes when each message from a plan is about to be
  processed and when a status object has completed.
* ``'bluesky.RE.msg`` --- A log record is emitted when each
  :class:`~bluesky.utils.Msg` is about to be processed.
* ``'bluesky.RE.state`` --- A log record is emitted when the RunEngine's state
  changes.

There are also some module-level loggers for specific features.

Formatter
---------

.. autoclass:: bluesky.log.LogFormatter

Global Handler
---------------

Following Python's recommendation, bluesky does not install any handlers at
import time, but it provides a function to set up a basic useful configuration
in one line, similar to Python's :py:func:`logging.basicConfig` but with some
additional options---and scoped to the ``'bluesky'`` logger with bluesky's
:class:`bluesky.log.LogFormatter`. It streamlines common use cases without
interfering with more sophisticated use cases.

We recommend that facilities using bluesky leave this function for users and
configure any standardized, facility-managed logging handlers separately, as
described in the next section.

.. autofunction:: bluesky.log.config_bluesky_logging
.. autofunction:: bluesky.log.get_handler

Advanced Example
================

The flow of log event information in loggers and handlers is illustrated in the
following diagram:

.. image:: https://docs.python.org/3/_images/logging_flow.png

For further reference, see the Python 3 logging howto:
https://docs.python.org/3/howto/logging.html#logging-flow

As an illustrative example, we will set up two handlers using the Python
logging framework directly, ignoring bluesky's convenience function.

Suppose we set up a handler aimed at a file:

.. code-block:: python

    import logging
    file_handler = logging.FileHandler('bluesky.log')

And another aimed at `Logstash <https://www.elastic.co/products/logstash>`_:

.. code-block:: python

    import logstash  # requires python-logstash package
    logstash_handler = logstash.TCPLogstashHandler(<host>, <port>, version=1)

We can attach the handlers to the bluesky logger, to which all log records
created by bluesky propagate:

.. code-block:: python

    logger = logging.getLogger('bluesky')
    logger.addHandler(logstash_handler)
    logger.addHandler(file_filter)

We can set the verbosity of each handler. Suppose want maximum verbosity in the
file but only medium verbosity in logstash.

.. code-block:: python

    logstash_handler.setLevel('INFO')
    file_handler.setLevel('DEBUG')

Finally, ensure that "effective level" of ``logger`` is at least as verbose as
the most verbose handler---in this case, ``'DEBUG'``. By default, at import,
its level is not set

.. ipython:: python
   :verbatim:

    logging.getLevelName(logger.level)
    'NOTSET'

and so it inherits the level of Python's default
"handler of last resort," :py:obj:`logging.lastResort`, which is ``'WARNING'``.

.. ipython:: python
   :verbatim:

    logging.getLevelName(logger.getEffectiveLevel())
    'WARNING'

In this case we should set it to ``'DEBUG'``, to match the most verbose level
of the handler we have added.

.. code-block:: python

   logger.setLevel('DEBUG')

This makes DEBUG-level records *available* to all handlers. Our logstash
handler, set to ``'INFO'``, will filter out DEBUG-level records.

To globally disable the generation of any log records at or below a certain
verbosity, which may be helpful for optimizing performance, Python provides
:py:func:`logging.disable`.
.. _hardware_interface_packages:

Hardware Interface Packages
===========================

The Bluesky library does not provide direct support for communicating with real hardware.
Instead, we define a high-level abstraction: the :ref:`Bluesky Hardware Interface <hardware_interface>`.
This allows different experimentalists to use different hardware control systems.

The following packages provide support for real hardware communication from Bluesky:

=============  ================================================================================
Ophyd_         EPICS_ integration for Bluesky. Reference implementation for hardware interface.
Instrbuilder_  Lightweight package with a focus on SCPI_.
Ophyd-Tango_   Tango_ integration for Bluesky. Incomplete and experimental early work.
pycertifspec_  Communication with SPEC_ instruments.
yaqc-bluesky_  yaq_ integration for Bluesky.
=============  ================================================================================

Importantly, you may mix hardware interfaces from multiple different packages within the same RunEngine.
Please note that the above packages are developed and maintained separately from Bluesky itself.

Are you maintaining a Python Package which provides hardware communication functionality?
The :ref:`Bluesky Hardware Interface <hardware_interface>` is a simple set of attributes and methods that can easily be added to your existing classes.
Please consider supporting our interface specification to unlock the full capabilities of the Bluesky ecosystem for your supported hardware.
Let us know if you add Bluesky support so we can add you to the above list.

.. _Instrbuilder: https://lucask07.github.io/instrbuilder/build/html/
.. _SCPI: https://en.wikipedia.org/wiki/Standard_Commands_for_Programmable_Instruments
.. _Ophyd: https://blueskyproject.io/ophyd/
.. _EPICS: https://epics-controls.org/
.. _Ophyd-Tango: https://github.com/bluesky/ophyd-tango
.. _pycertifspec: https://github.com/SEBv15/pycertifspec
.. _SPEC: https://www.certif.com/content/spec/
.. _Tango: https://www.tango-controls.org/
.. _yaqc-bluesky: https://github.com/bluesky/yaqc-bluesky
.. _yaq: https://yaq.fyi/
.. currentmodule:: bluesky.simulators

Simulation and Error Checking
=============================

Bluesky provides three different approaches for simulating a plan without
actually executing it:

1. Introspect a plan by passing it to a "simulator" instead of a RunEngine.
2. Execute a plan with the real RunEngine, but use simulated hardware objects.
3. Redefine the RunEngine commands to change their meanings.

Approaches (1) and (2) are the most straightforward and most common.

Introspection
-------------

Recall that plans yield messages that *describe* what should be done; they
do not communicate with hardware directly. Therefore it's easy to use (or
write) a simple function that iterates through the plan and summarizes or
analyzes its actions.

.. autosummary::
   :toctree: generated
   :nosignatures:

   summarize_plan
   plot_raster_path
   check_limits

Summarize
^^^^^^^^^

The simulator :func:`summarize_plan` print a summary of what a plan would do if
executed by the RunEngine.

.. ipython:: python

    from bluesky.simulators import summarize_plan
    from ophyd.sim import det, motor
    from bluesky.plans import scan
    summarize_plan(scan([det], motor, 1, 3 ,3))

To see the unabridged contents of a plan, simply use the builtin Python
function :func:`list`. Note that it is not possible to summarize plans that
have adaptive logic because their contents are determined dynamically during
plan executation.

.. ipython:: python

    list(scan([det], motor, 1, 3 ,3))

Check Limits
^^^^^^^^^^^^

.. ipython:: python
    :suppress:

    motor.limits = (-1000, 1000)

Suppose that this motor is configured with limits on its range of motion at +/-
1000. The :func:`check_limits` simulator can verify whether or not a plan will
violate these limits, saving you from discovering this part way through a long
experiment.

.. ipython:: python
    :okexcept:

    from bluesky.simulators import check_limits

    check_limits(scan([det], motor, 1, 3 ,3))  # no problem here
    check_limits(scan([det], motor, 1, -3000, 3000))  # should raise an error

Simulated Hardware
------------------

.. warning::

    This feature has recently been changed, and it has yet to be documented.
    
Customizing RunEngine Methods
-----------------------------

The RunEngine allows you to customize the meaning of commands (like 'set' and
'read'). One could use this feature to create a dummy RunEngine that, instead
of actually reading and writing to hardware, merely reports what it *would*
have done.

.. automethod:: bluesky.run_engine.RunEngine.register_command
   :noindex:

.. automethod:: bluesky.run_engine.RunEngine.unregister_command
   :noindex:
