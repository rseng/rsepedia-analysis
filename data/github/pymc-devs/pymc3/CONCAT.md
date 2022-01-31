# Architecture
This document describes the high-level architecture of PyMC.

# Bird's Eye View
[comment]: <> (https://drive.google.com/file/d/1lfEzokkNUJr_JIeSDQfha5a57pokz0qI)
![Architecture](docs/Architecture.png)
Lets you define probabilistic graphs or models that can be easily used to compute log probabilities for posterior
inference or to draw random samples for prior and posterior prediction.

PyMC includes a few inference techniques, in particular:
* Markov chain Monte Carlo
* Variational Inference
* Sequential Monte Carlo

It also contains numerous others pieces of functionality such as GraphviZ model visualization tools
as well as various mathematical helper functions.

The most central pieces functionality of PyMC are shown visually below, as well as their
relation to other major packages. Not all modules are shown, either because
they are smaller or self explanatory in scope, or they're pending
deprecation

## Functionality not in PyMC
It is easier to start with functionality that is not present in PyMC but
rather deferred to outside libraries. If seeking to understand any
of the topics below refer to that specific library

### Aesara
* Gradient computation
* Random number generation
* Low level tensor operation definition
* Low level operation graphs

### ArviZ
* Plotting e.g. Trace plots, rank plots, posterior plots
* MCMC sampling diagnostics e.g. Rhat, Effective Sample Size.
* Model comparison, particularly efficient leave-one-out cross-validation approximation
* Inference Data structure


# Modules
The codebase of PyMC is split among single Python file modules at the root
level, as well as directories with Python code for logical groups of functionality.
Admittedly the split between single `.py` module or directory is not defined by a strict
criteria but tends to occur when single `.py` files would be "too big".
We will with the modules needed implement "simple MCMC" model shown below
before detailing the remaining modules, such as Variational Inference, Ordinary Differential Equations,
or Sequential Monte Carlo.

```python
with pm.Model() as model:
  theta = pm.Beta("theta", alpha=1, beta=2)
  p = pm.Beta("n", p=theta, n=2, observed=[1,2])
  inf_data = pm.sample()


```

## {mod}`pymc.model`
Contains primitives related model definition and methods used for evaluation of the model.
In no particular order they are

* `ContextMeta`: The context manager that enables the `with pm.Model() as model` syntax
* {class}`~pymc.Factor`: Defines the methods for the various logprobs for models
* `ValueGrad` which handles the value and gradient and is the main connection point to Aesara
* `Deterministic` and `Potential`: Definitions for two pieces of functionality useful in some model definitions

## distributions/
Contains multiple submodules that define distributions,  as well as logic that aids in distributions usage.
Important modules to note are

* `distribution.py`: This contains parent class for all PyMC distributions.
  Notably the `distribution.distribution` class contains the `observed` argument which in PyMC differentiates
  a random variable distribution from a likelihood distribution.

* `logprob.py`: This contains the log probability logic for the distributions themselves.
  The log probability calculation is deferred to Aesara

* `dist_math.py`: Various convenience operators for distributions.
  This includes mathematical operators such as `logpower` or `all_true`methods.
  It also contains a suite of lognormal methods and transformation methods

## /sampling.py
Interface to posterior, prior predictive, and posterior sampling as well as various methods to identify and initialize
stepper methods. Also contains logic to check for "all continuous" variables and initialize NUTS

## step_methods/
Contains various step methods for various sampling algorithms, such as MCMC, and SMC. `step_methods.hmc` includes
the Hamiltonian Monte Carlo sampling methods as well as helper functions such as the integrators used for those methods

## tests/
All tests for testing functionality of codebase. All modules prefixed with `test_` are tests themselves, whereas all
other modules contain various supporting code such as fixtures, configurations, etc
# Main Governance Document

## The Project

The PyMC Project (The Project) is an open source software project
affiliated with the 501c3 NumFOCUS Foundation. The goal of The Project is to
develop open source software and deploy open and public websites and services
for reproducible, exploratory and interactive computing.
The main focus of The Project is in scientific and statistical computing.
The Software developed
by The Project is released under OSI approved open source licenses,
developed openly and hosted in public GitHub repositories under the
[pymc-devs GitHub organization](https://github.com/pymc-devs). Examples of
Project Software include the PyMC library and its documentation, etc.
The Services run by The Project consist of public websites and web-services
that are hosted at [http://docs.pymc.io](https://docs.pymc.io)

The Project is developed by a team of distributed developers, called
Contributors. Contributors are individuals who have contributed code,
documentation, designs or other work to one or more Project repositories,
or who have done significant work to empower the Community,
participating on [Discourse](https://discourse.pymc.io),
organizing [PyMCon](https://pymcon.com) or helped on other platforms and events.
Anyone can be a Contributor.
Contributors can be affiliated with any legal entity or none.
The foundation of Project participation is openness and transparency.

There have been over 250 Contributors to the Project, their contributions are listed in the
logs of the PyMC GitHub repositories as well as those of associated projects and venues.

The Project Community consists of all Contributors and Users of the Project.
Contributors work on behalf of and are responsible to the larger Project
Community and we strive to keep the barrier between Contributors and Users as
low as possible.

The Project is formally affiliated with the 501c3
[NumFOCUS Foundation](http://numfocus.org), which serves as its fiscal
sponsor, may hold project trademarks and other intellectual property, helps
manage project donations and acts as a parent legal entity. NumFOCUS is the
only legal entity that has a formal relationship with the project (see
Institutional Partners section below).

## Governance

This section outlines the governance and leadership model of The Project.

The foundations of Project governance are:

- Openness & Transparency
- Active Contribution
- Institutional Neutrality

Traditionally, Project leadership was provided by a BDFL (Chris Fonnesbeck) and
subset of Contributors, called Core Developers, whose active and consistent
contributions have been recognized by their receiving “commit rights” to the
Project GitHub repositories. In general all Project decisions were made through
consensus among the Core Developers with input from the Community. The BDFL
could, but rarely chose to, override the Core Developers and make a final
decision on a matter.

While this approach has served us well, as the Project grows and faces more
legal and financial decisions and interacts with other institutions, we see a
need for a more formal governance and organization model.
We view this governance model as the formalization of what we are already doing,
rather than a change in direction.

## Community and Team Architecture
The PyMC community is organized in an onion-like fashion.
The tiers relevant to the project governance are listed below sorted by
increasing responsibility. Due to the onion-like structure, members of a group are
also members of all the groups listed above:

* Contributors
* Recurring Contributors
* Core Contributors
* Steering Council
* BDFL

Recurring Contributors comprise what we understand as the PyMC Team.
The Team will generally act as a single unit, except for some specific
questions where dedicated teams will prevail.
The PyMC project currently has Developer, Documentation and Discourse teams.
Team members can be part of one, some or none of these dedicated teams.
The diagram below should help illustrate this idea.

<img src="docs/community_diagram.png" alt="community diagram" width="600" height="400">

Anyone working with The Project has the responsibility to personally uphold
the Code of Conduct. Recurrent Contributors have the additional responsibility
of _enforcing_ the Code of Conduct to maintain a safe community.

## Recurring Contributors
Recurring Contributors are those individuals who contribute recurrently to the
project and can provide valuable insight on the project.
They are therefore actively consulted and can participate in the same communication
channels as Core Contributors. However, unlike Core Contributors,
Recurrent Contributors don't have voting, managing or writing rights.

In practice, this translates in participating from private team discussions
(i.e. in Slack or live meetings) but not being able to vote Steering Council
members nor having commit rights on GitHub.

The Recurrent Contributor position will often be an intermediate step for people
in becoming Core Contributors once their contributions are frequent enough
and during a sustained period of time.
But it is also an important role by itself for people who want to be part of
the project on a more advisory-like role, as they for example might not have
the time availability or don't want the responsibilities that come
with being a Core Contributor.

### Recurring Contributor membership
Recurring Contributors can nominate any Contributor to participate in the
Project private communication channels (i.e. Slack public channel)
and become a Recurring Contributor.
For the nomination to go forward, it has to be ratified by the Steering Council.
For a nomination to be rejected, clear reasoning behind the decision must be
shared with the rest of the team. People whose nomination has been rejected can
be nominated at any time again in the future, three months after the previous
nomination at the earliest. The nomination process is explained below
in more detail in a section of its own.

Interns and contractors are added to the team as Recurrent Contributors.
We consider the selection/hiring process to replace the nomination process.

#### Current Recurring Contributors
Contributors who are also part of a dedicated team or are institutional
contributors will have so indicated after their name.
Dedicated teams only cover a small part of the work needed to
get the project going, tasks like fundraising, outreach and marketing,
or organizing events for example don't (yet) have a dedicated team.
Contributors don't need to be part of any dedicated team.

* Abhipsha Das (docs)
* Benjamin Vincent (docs - PyMC Labs)
* Jon Sedar
* Kaustubh Chaudhari (dev)
* Larry Dong (dev)
* Lorenzo Toniazzi (docs)
* Martin Ingram (community)
* Olga Khan (docs)
* Peadar Coyle
* Raul Maldonado (docs)

## Core Contributors
Core Contributors are those individuals entrusted with the development and
well being of the Project due to their frequency of quality contributions over
a sustained period of time.

They are the main governing and decision body
of the Project and are therefore given voting and managing rights to the Project
services (i.e. commit rights on GitHub or moderation rights on Discourse).

Team memberships for Core Contributors refer to their Core Contributor
role, as experienced and entrusted members of the team they are considered
Recurrent Contributors for the rest of the teams and given permissions accordingly

The exact permissions of all Core Contributors may therefore not be the same
and depend on their team memberships. Even if they have commit rights,
Core Contributors should still have their pull requests reviewed by at least
one other Core Contributor before merging.
In rare circumstances, a Core Contributor may bypass this review
if in their best judgement it is appropriate to do so,
but such expedited PR merges must be justifiable and
ultimately subject to review post hoc, if necessary.
If overstepping, Core Contributors can also be subject to a vote
of no confidence (see below) and see their permissions revoked.

### Core Contributor membership
To become a Core Contributor, one must already be a Recurring Contributor.
Core Contributors can nominate any Recurring Contributor to become a
Core Contributor. For the nomination to go forward, it has to be
ratified by the Steering Council.
For a nomination to be rejected, clear reasoning behind the decision must be
shared with the rest of the team. People whose nomination has been rejected can
be nominated at any time again in the future, three months after the previous
nomination at the earliest. The nomination process is explained below
in more detail in a section of its own.

### Current Core Contributors
Contributors who are also part of a dedicated team or are institutional
contributors will have so indicated after their name.

Dedicated teams only cover a small part of the work needed to
get the project going, tasks like fundraising, outreach and marketing,
or organizing events for example don't (yet) have a dedicated team.
Contributors don't need to be part of any dedicated team.

* Adrian Seyboldt (dev - PyMC Labs)
* Alex Andorra (dev - PyMC Labs)
* Austin Rochford
* Bill Engels (dev)
* Brandon T. Willard (dev)
* Chris Fonnesbeck (dev, docs)
* Christian Luhmann (discourse)
* Colin Carroll (dev)
* Eelke Spaak (dev)
* Eric Ma (dev - PyMC Labs)
* George Ho (dev)
* Junpeng Lao (dev, discourse)
* Luciano Paz (dev - PyMC Labs)
* Martina Cantaro (docs)
* Maxim Kochurov (dev - PyMC Labs)
* Meenal Jhajharia (docs)
* Michael Osthege (dev)
* Oriol Abril-Pla (docs, discourse)
* Osvaldo Martin (dev, docs)
* Ravin Kumar (dev, discourse, docs)
* Ricardo Vieira (dev, discourse)
* Robert P. Goldman (dev)
* Sayam Kumar (dev, docs)
* Thomas Wiecki (dev, discourse - PyMC Labs)

## Steering Council

The Project will have a Steering Council that consists of Project Contributors
who have produced contributions that are substantial in quality and quantity,
and sustained over at least one year. The overall role of the Council is to
ensure, through working with the BDFL and taking input from the Community, the
long-term well-being of the project, both technically and as a community.

The Steering Council will have between 4 and 7 members with at least one member
per dedicated team.
No more than 2 Council Members can report to one person or company
(including Institutional Partners) through employment or
contracting work (including the reportee, i.e. the reportee + 1 is the max).


During the everyday project activities, council members participate in all
discussions, code review and other project activities as peers with all other
Contributors and the Community. In these everyday activities, Council Members
do not have any special power or privilege through their membership on the
Council. However, it is expected that because of the quality and quantity of
their contributions and their expert knowledge of the Project Software and
Services that Council Members will provide useful guidance, both technical and
in terms of project direction, to potentially less experienced contributors.

The Steering Council and its Members play a special role in certain situations.
In particular, the Council may:

- Make decisions about the overall scope, vision and direction of the
  project.
- Make decisions about strategic collaborations with other organizations or
  individuals.
- Make decisions about specific technical issues, features, bugs and pull
  requests.
- Make decisions about the Services that are run by The Project and manage
  those Services for the benefit of the Project and Community.
- Make decisions when regular community discussion doesn’t produce consensus
  on an issue in a reasonable time frame.

### Current Steering Council

The current Steering Council membership comprises:

- Chris Fonnesbeck (dev, docs)
- Junpeng Lao (dev, discourse)
- Oriol Abril-Pla (docs, discourse)
- Ravin Kumar (dev, discourse, docs)
- Thomas Wiecki (dev, discourse - PyMC Labs)

Note that as explained in the [community architecture section](#community-and-team-architecture)
and as indicated again in the description of the Steering Council above,
Council members are first and foremost Core Contributors, and have no special
power or privilege in everyday activities.
To emphasize that, Steering Council members are listed both in the current core
contributors section and in this section even if redundant.

### Council membership

To become eligible for being a Steering Council Member an individual must be a
Core Contributor who has produced contributions that are substantial in
quality and quantity, and sustained over at least one year.

Similarly to when nominating new team members, when considering potential
Council Members one should look at candidates with a
comprehensive view of their contributions. This will include but is not limited
to code, code review, infrastructure work, forum and chat participation,
community help/building, education and outreach, design work, etc. We are
deliberately not setting arbitrary quantitative metrics (like “100 commits in
this repo”) to avoid encouraging behavior that plays to the metrics rather than
the project’s overall well-being. We want to encourage a diverse array of
backgrounds, viewpoints and talents in our team, which is why we explicitly do
not define code as the sole metric on which council membership will be
evaluated. See the section on election process for more details.

Council membership is assigned for a two year period, with no limit on how many
periods can be served.

Council members can renounce at any time and are
encouraged to do so if they foresee they won't be able to attend their
responsibilities for an extended interval of time.

If a Council member becomes inactive in the project for a period of six months,
they will be considered for removal from the Council. Before removal, inactive
Member will be approached by the BDFL to see if they plan on returning to
active participation. If not they will be removed immediately, as they
are effectively renouncing to their position.
If they plan on returning to active participation soon, they will be
given a grace period of six months. If they don’t return to active participation
within that time period they will be removed without
further grace period. All former Council members can be considered for
membership again at any time in the future.
Retired Council members will be listed on the project website, acknowledging
the period during which they were active in the Council.

The Council reserves the right to eject current Members, other than the BDFL,
if they are deemed to be actively harmful to the project’s well-being, and
attempts at communication and conflict resolution have failed. See
the section on votes of no-confidence for details on the process.

### Private communications of the Council

Unless specifically required, all Council discussions and activities will be
public and done in collaboration and discussion with the Project Team
and also the Community when possible.
The Council will have a private mailing list that will be used
sparingly and only when a specific matter requires privacy. When private
communications and decisions are needed, the Council will do its best to
summarize those to the Team after eliding personal/private/sensitive
information that should not be posted to the public internet.

### Subcommittees

The Council can create subcommittees that provide leadership and guidance for
specific aspects of the project. Like the Council as a whole, subcommittees
should conduct their business in an open and public manner unless privacy is
specifically called for.

Even if the BDFL does not sit on a specific subcommittee, they still retain
override authority on the subcommittee's decisions. However, it is expected that
they will appoint a delegate to oversee the subcommittee's decisions, and
explicit intervention from the BDFL will only be sought if the committee
disagrees with the delegate's decision and no resolution is possible within the
subcommittee. This is a different situation from a BDFL delegate for a specific
decision, or a recusal situation, in which the BDFL gives up their authority
to someone else in full.

### NumFOCUS Subcommittee

The Council will maintain one narrowly focused subcommittee to manage its
interactions with NumFOCUS.

- The NumFOCUS Subcommittee is comprised of 5 persons who manage project
  funding that comes through NumFOCUS. It is expected that these funds will
  be spent in a manner that is consistent with the non-profit mission of
  NumFOCUS and the direction of the Project as determined by the full
  Council.
- This Subcommittee shall NOT make decisions about the direction, scope,
  technical or financial direction of the Project.

#### NumFOCUS subcommittee membership
This Subcommittee will have 5 members. With at least
2 members being on the Steering Council. No more
than 2 Subcommitee Members can report to one person or company through
employment or contracting work (including the reportee, i.e.
the reportee + 1 is the max).
This avoids effective majorities resting on one person.

Any Core Contributor is eligible for the NumFOCUS subcommittee.

#### Current NumFOCUS Subcommitee
The current NumFOCUS Subcommittee consists of:

- Peadar Coyle
- Chris Fonnesbeck
- John Salvatier
- Jon Sedar
- Thomas Wiecki

## BDFL

The Project will have a BDFL (Benevolent Dictator for Life), who is currently
Chris Fonnesbeck. As Dictator, the BDFL has the authority to make all final
decisions for The Project. As Benevolent, the BDFL, in practice chooses to
defer that authority to the consensus of the community discussion channels and
the Steering Council. It is expected, and in the past has been the
case, that the BDFL will only rarely assert their final authority. Because
rarely used, we refer to BDFL’s final authority as a “special” or “overriding”
vote. When it does occur, the BDFL override typically happens in situations
where there is a deadlock in the Steering Council or if the Steering Council
asks the BDFL to make a decision on a specific matter. To ensure the
benevolence of the BDFL, The Project encourages others to fork the project if
they disagree with the overall direction the BDFL is taking. The BDFL is chair
of the Steering Council (see below) and may delegate their authority on a
particular decision or set of decisions to any other Council member at their
discretion.

The BDFL can appoint their successor, but it is expected that the Steering
Council would be consulted on this decision. If the BDFL is unable to appoint a
successor, the Steering Council will make a suggestion or suggestions to the
Main NumFOCUS Board. While the Steering Council and Main NumFOCUS Board will
work together closely on the BDFL selection process, the Main NUMFOCUS Board
will make the final decision.


## Conflict of interest

It is expected that the BDFL, Council Members and Contributors will be
employed at a wide range of companies, universities and non-profit organizations.
Because of this, it is possible that Members will have conflict of interests.
Such conflict of interests include, but are not limited to:

- Financial interests, such as investments, employment or contracting work,
  outside of The Project that may influence their work on The Project.
- Access to proprietary information of their employer that could potentially
  leak into their work with The Project.

All members of the Council, BDFL included, shall disclose to the rest of the
Council any conflict of interest they may have. Members with a conflict of
interest in a particular issue may participate in Council discussions on that
issue, but must recuse themselves from voting on the issue. If the BDFL has
recused themselves for a particular decision, they will appoint a substitute
BDFL for that decision.

## Voting processes
### Nomination process
> Used when adding members to the team as recurrent or core contributors.

A nomination process is triggered automatically whenever a team member
requests so on one of the team's communication channels
(public Slack channels at the day of writing, preferably `#general`).
Nomination should be explicit regarding which roles and teams are
suggested, but the council makes the final decision on
dedicated team membership.
Again, note that team members don't need to be part of any
dedicated team to be recurrent nor core contributors.

After this happens, the Steering Council will reach out to the candidate
to see if they accept the nomination. If the nomination is accepted
it will be considered by the Steering Council.
At their earliest convenience and no later than two weeks, the Steering
Council will vote on the nominee using the process below on
Steering Council decisions.

Voting will be private to the Steering Council only with results published
on the nomination request.
In the case of a rejection, results must include the reasons behind
the decision (i.e. the time since starting to contribute is deemed
too short for now).

### Steering Council decisions
By and large we expect the decisions in PyMC to be made _ad hoc_
and require little formal coordination and with the community at large.
However, for controversial proposals and new team members the council can
intervene to make the final decision in a group vote.

#### Call for a vote
Core Contributors can call for a vote to resolve a target issue
they feel has been stale for too long and for which
informal consensus appears unlikely.
For a vote to be called, the target issue or discussion post (i.e. on Discourse)
must be at least 1 month old.

To do so, they have to open a proposal issue ticket labeled "Council Vote".
The proposal issue should contain a link to the target issue and
a proposal on how to resolve it.
Proposals should include a statement making clear what it means to
"agree" or to "disagree".

Before voting starts, at least 3 days will be left for Core Contributors
to raise doubts about the proposal's phrasing, no extra discussion will
take place in the proposal issue.
Proposal issues should be locked from creation to prevent attracting
discussion from people not familiar with the decision process.

A vote is also called automatically whenever someone is nominated to
be added to the team.

The Steering Council can also call a vote on their own in order
to eject a Core contributor.

Upon ejecting a core contributor the council must publish an issue ticket,
or public document detailing the
* Violations
* Evidence if available
* Remediation plan (if necessary)
* Signatures majority of council members to validate correctness and accuracy

#### Voting

* Each Council Member will vote either "Yes", "No", or "Abstain".
* It is recommended that all Council Members expose their reasons when voting.
  "No" votes, however, must list the reasons for disagreement.
  Any "No" vote with no reason listed will be considered a "Abstain" vote.
* An absence of vote is considered as "Abstain".
* Voting will remain open for at least 3 days.
* For the proposal to pass, at least 60% of the council must vote "Yes", and no more than 20% can vote "No".

For decisions about the project the Council will perform it directly
on the proposal issue. For decisions about people,
such as electing or ejecting Team Members, the Council will vote privately.
However the decision will be posted publicly in an issue ticket.

### Vote of no confidence
In exceptional circumstances, Council Members as well as Core Contributors
may remove a sitting council member via a vote of no confidence.
Core contributors can also call for a vote to remove the entire council
-- in which case, Council Members do not vote.
A no-confidence vote is triggered when a Core Contributor calls for one
publicly on an appropriate project communication channel,
and two other core team members second the proposal.
The initial call for a no-confidence vote must specify which type is intended.

The vote lasts for two weeks, and the people taking part in it vary:
* If this is a single-member vote
  all Core contributors (including council members) vote,
  and the vote is deemed successful if at least two thirds of voters
  express a lack of confidence.

  Such votes can target removing from the Council
  (while continuing to be a Core Contributor) or a Core Contributor removal,
  which should be clear from the initial vote call. Council members
  can be called for a Core Contributor removal.
* If this is a whole-council vote, then it was necessarily called by
  Core contributors (since Council members can’t remove the whole Council)
  and only Core contributors vote.
  The vote is deemed successful if at least two thirds of voters
  express a lack of confidence.

After voting:
* If a single-member vote on a council member succeeds, then that member is
  removed from the council and the resulting vacancy can be handled in the usual way.
* If a single-member vote on a core contributor succeeds, their permissions are
  revoked and would have to wait six months to be eligible for core contributor
  nomination again.
* If a whole-council vote succeeds, the council is dissolved
  and a new council election is triggered immediately.
  In such case, members of the dissolved council continue to be Core Contributors.

### Election process
> Used when choosing the steering council and it's subcommittees

The log of past election processes is kept on [Discourse](https://discourse.pymc.io/tag/elections)

#### Nominations
* Nominations are taken over a single Discourse topic (with at least the `elections` tag)
  over the course of 2 weeks.
* Only Core Contributors may nominate candidates for Council membership
* Self Nominations are allowed
* There are no limits to the number of people that can be nominated by a single Core Contributor
* Once someone has been nominated, extra nominations are ignored. All candidates are treated equally
  in the election independently of potential differences (i.e. number of likes) in their respective nominations.
* At the conclusion of the 2 weeks, the list of nominations is posted on the ticket and this ticket is closed.

#### Voting

* Voting occurs over a period of at least 1 week, at the conclusion of the nominations.
  Voting is blind and mediated by either an application or a third party like NumFOCUS.
  Each voter can vote zero or more times, once per each candidate.
  As this is not about ranking but about capabilities,
  voters vote on a yes/abstain/no basis per candidate --
  “Should this person be on the Steering Council?
  Would I trust this person to lead PyMC?”.
  The absence of vote is considered abstain.
* Candidates are evaluated independently,
  each candidate having 60% or more of yes votes and less or
  equal than 20% of no votes is chosen.
  If the number of chosen candidates matches the number or range set for the
  council/subcommittee being chosen and all extra constrains are met,
  all candidates are confirmed and the election process stops here.
* In the event that not enough/too many candidates were confirmed or
  the membership constraints were not met,
  candidates are ranked by interpreting yes=+1, abstain=0 and no=-1.
  * If too many candidates were confirmed,
    the {max_range} number of candidates with higher rank are elected.
  * If not enough candidates were chosen,
    the {min_range} number of candidates with higher rank are elected.
  * If membership constraints were not met, candidates are selected
    progressively by rank if they meet the membership requirements.
    If for example out of 7 candidates for the NumFOCUS subcommittee
    only three are on the Steering Council and they were ranked 5th-7th,
    in order to meet the membership constraints, the person ranked 4th
    would not be elected, as their election would prevent membership
    requirements from being met.
* In the event of a tie there will be a runoff election for the tied candidates.
  To avoid further ties and discriminate more among the tied candidates,
  this vote will be held by Majority Judgment:
  for each candidate, voters judge their suitability for office as either
  "Excellent", "Very Good", "Good", "Acceptable", "Poor", or "Reject".
  Multiple candidates may be given the same grade by a voter.
  The candidate with the highest median grade is the winner.
* If more than one candidate has the same highest median-grade,
  the Majority Judgment winner is discovered by removing (one-by-one)
  any grades equal in value to the shared median grade from
  each tied candidate's total.
  This is repeated until only one of the previously tied candidates
  is currently found to have the highest median-grade.
* If ties are still present after this second round, the winner will be chosen at random. First we make a alphanumerically sorted list of the names in the tie. Then we will draw one prior predictive sample from a `pm.Categorical` distribution over the elements in the list to determine the winner.
* At the conclusion of voting, all the results will be posted. And at least 24 hours will be left to challenge the election result in case there were suspicions of irregularities or the process had not been correctly carried out.

## Leaving the project
Any contributor can also voluntarily leave the project by notifying the community through a public means or by notifying the entire council. When doing so, they can add themselves to the alumni section below if desired.

People who leave the project voluntarily can rejoin at any time.

## Team Organization
As stated previously, The Team will generally act as a single unit,
except for some specific questions where dedicated teams will prevail.
These dedicated teams have no difference in how they are governed.
Decisions should be reached by consensus within the team with the Steering
Council and the BDFL acting if necessary.

The dedicated teams are work units with two main objectives: better
distributing the work related to The Project, and to better showcase all the different tasks
involved in The Project to attract more diverse Contributors.

The PyMC project currently has Developer, Documentation and Discourse teams.
Team members can be part of one, some or none of these dedicated teams.

Team members are expected to participate and join these dedicated teams
organically. That is, the Steering Council will take part actively
in team assignments if they are part of a role change with the respective
nomination.

### Developer Team
The focus of the developer team is the probabilistic programming library
and flagship of The Project, [PyMC](https://github.com/pymc-devs/pymc).

For current members of the developer team, refer to the recurrent and
core contributor membership sections.

### Documentation Team
The focus of the documentation team is ensuring the PyMC library
is well documented, building and maintaining the infrastructure needed
for that aim, and making sure there are resources to learn
Bayesian statistics with PyMC.
It is not the goal nor responsibility of the Documentation team to
write all the documentation for the PyMC library.

For current members of the documentation team, refer to the recurrent and
core contributor membership sections.

### Discourse team
The focus of the Discourse team is managing and energizing the PyMC Discourse.

For current members of the discourse team, refer to the recurrent and
core contributor membership sections.

### "No-team" tasks
All tasks related to the project that are not specifically listed in the
description of a dedicated team are the responsibility of the PyMC team
as a whole. At the time of writing, this includes but is not limited to:
fundraising, issue triaging, running PyMC related events like PyMCon or
sprints, outreach, or presence on social networks.

### Team structure in practice
This section describes how members of the PyMC team are given
permissions to the relevant services based on their roles
and dedicated team affiliations.

#### GitHub
Two of the teams are currently structured about GitHub-centric tasks, so the
permissions on GitHub repositories is mapped to team membership and role
within the team. The team defines to which repositories the permissions
are given, the role defines the type of permissions given:

Role:
- Recurring Contributors are given triage permissions
- Core Contributors are given write permissions

Team:
* Development team members are given permissions to the [pymc](https://github.com/pymc-devs/pymc)
  and [pymc-experimental](https://github.com/pymc-devs/pymc-experimental) repository.
* Documentation team members are given permissions to [pymc-examples](https://github.com/pymc-devs/pymc-examples)
  and [resources](https://github.com/pymc-devs/resources) repositories.

In addition, Council members are given admin rights to all repositories within
the [pymc-devs](https://github.com/pymc-devs) organization.

#### Discourse
Similarly to the above section, Discourse permissions are also mapped to the discourse team
and the two contributor roles.

Role:
- Recurring Contributors are added to the [PyMC_team](https://discourse.pymc.io/g/PyMC_devs)
  group and are given the "leader" trust level.
- Core Contributors are given [moderator permissions](https://discourse.pymc.io/g/moderators)
  if possible.

#### Accounts and services ownership and administration
The PyMC Project also has accounts and hosts services on several platforms
such as GitHub, Discourse, Twitter, ReadTheDocs or Medium.

If possible, all Council Members and relevant Core Contributors should have
admin rights on those platforms.
If this were not possible, admin rights should be distributed between
Council Members and relevant Core Contributors and establish a rotation
of the admin rights every 1-2 years.

#### Permission examples
This section lists some imaginary contributors with their teams and roles to
provide examples on how to assign permissions:

<details><summary>See permission examples</summary>

- Arnau, recurrent contributor, discourse team
  * Added to the Discourse PyMC_team and given "leader" trust level
  * Added to all private communication channels
  * No permissions on any GitHub

- Berta, recurrent contributor, dev and doc teams
  * No permissions on Discourse
  * Added to all private communication channels
  * Triage permissions on pymc, pymc-experimental, pymc-examples and resources repositories
    of the pymc-devs organization

- Carme, core contributor, doc team
  * Added to the Discourse PyMC_team and given "leader" trust level
  * Added to all private communication channels
  * Write permissions on pymc-examples and resources repositories, triage permissions
    to pymc and pymc-experimental repositories
  * Admin access to ReadTheDocs accounts

- Dolors, core contributor, dev and discourse teams
  * Added to the Discourse PyMC_team and given "leader" trust level,
    given moderator permissions on Discourse (that might be rotating and not always active)
  * Added to all private communication channels
  * Write permissions on pymc and pymc-experimental repositories, triage permissions
    to pymc-examples and resources repositories

- Eudald, core contributor, no dedicated team membership
  * Added to the Discourse PyMC_team and given "leader" trust level
  * Added to all private communication channels
  * Triage permissions on all repositories
  * Access to pymc_devs twitter account as they are the main manager

</details>

## Institutional Partners and Funding

The PyMC Core Contributors (together with the BDFL and Steering Council)
are the primary leadership for the project. No
outside institution, individual or legal entity has the ability to own,
control, usurp or influence the project other than by participating in the
Project as Contributors and Council Members. However, because institutions are
the primary funding mechanism for the project, it is important to formally
acknowledge institutional participation in the project. These are Institutional
Partners.

An Institutional Contributor is any individual Core Contributor who
contributes to the project as part of their official duties at an Institutional
Partner. Likewise, an Institutional Council Member is any Project Steering
Council Member who contributes to the project as part of their official duties
at an Institutional Partner.

With these definitions, an Institutional Partner is any recognized legal entity
in the United States or elsewhere that employs at least one Institutional
Contributor or Institutional Council Member. Institutional Partners can be
for-profit or non-profit entities.

Institutions become eligible to become an Institutional Partner by
employing individuals who actively contribute to The Project as part
of their official duties. To state this another way, the only way for
an Institutional Partner to influence the project is by actively
contributing to the open development of the project, on equal terms
with any other member of the community of Contributors and Council
Members. Merely using PyMC Software or Services in an
institutional context does not allow an entity to become an
Institutional Partner. Financial gifts do not enable an entity to
become an Institutional Partner (see Sponsors below for financial gift recognition).
Once an institution becomes eligible
for Institutional Partnership, the Steering Council must nominate and
approve the Partnership.

If an existing Institutional Partner no longer has a contributing employee,
they will be given a one-year grace period for other employees to begin
contributing.

An Institutional Partner is free to pursue funding for their work on The
Project through any legal means. This could involve a non-profit organization
raising money from private foundations and donors or a for-profit company
building proprietary products and services that leverage Project Software and
Services. Funding acquired by Institutional Partners to work on The Project is
called Institutional Funding. However, no funding obtained by an Institutional
Partner can override The Project BDFL and Steering Council. If a Partner has
funding to do PyMC work and the Council decides to not pursue that
work as a project, the Partner is free to pursue it on their own. However in
this situation, that part of the Partner’s work will not be under the
PyMC banner and cannot use the Project trademarks in a way that
suggests a formal relationship.

To acknowledge institutional contributions, there are two level of Institutional
Partners, with associated benefits:

**Tier 1** = an institution with at least one Institutional Council Member

- Acknowledged on the PyMC websites, in talks and T-shirts.
- Ability to acknowledge their own funding sources on the PyMC
  websites, in talks and T-shirts.
- Unlimited participation in the annual Institutional Partners Workshop, held
  during the (planned) annual PyMC Project Retreat. This allows the
  Institutional Partner to invite as many of their own employees and funding
  sources and collaborators as they want, even if they are not project
  Contributors or Council Members.
- Ability to influence the project through the participation of their Council
  Member.
- Council Members are invited to the bi-annual PyMC Developer Meeting.

**Tier 2** = an institution with at least one Institutional Contributor

- Same benefits as Tier 1 level Partners, but:
- Only Institutional Contributors are invited to the Institutional Partners
  Workshop and bi-annual PyMC Developer Meeting

The PyMC project currently recognizes PyMC Labs as a Tier 1 Institutional Partner,
with Thomas Wiecki and Adrian Seyboldt as their institutional contributors
and council members.

## Sponsors
Sponsors are organizations that provide significant funding to the PyMC project
directly. Interested sponsors are encouraged to reach out
to the Steering Council to arrange the sponsorship and recognition.

The PyMC project reserves the right to not approve a sponsorship if
the goals or culture of the prospective sponsor are deemed incompatible
with the goals of the project. In such case, like with any negative vote
from the Steering Council, reasoning behind the decision will be provided.

Sponsors will be recognized by placing their logo on the PyMC website but will have
no extra benefits related to The Project. Note that PyMCon sponsors may have
extra benefits but those will be related to the conference, not the Project.

## Team Alumni

* Agustina Arroyuelo (GSoC 2018)
* Anand Patil
* David Huard
* Demetri Pananos (GSoC 2019)
* John Salvatier
* Joseph Willard (GSoC 2019)
* Juan Martín Loyola (GSoC 2019)
* Rasul Karimov (GSoC 2020)
* Sharan Yalburgi (GSoC 2018)
* Taku Yoshioka
* Tirth Patel (GSoC 2020)
# Release Notes

<!--
⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠
Do not create individual sections for the beta releases.
Instead update the vNext section until 4.0.0 is out.
⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠⚠
-->

## PyMC vNext (4.0.0b1 → 4.0.0b2 → 4.0.0b3 → 4.0.0)
⚠ The changes below are the delta between the upcoming releases `v3.11.5` →...→ `v4.0.0`.

### No-yet working features
We plan to get these working again, but at this point their inner workings have not been refactored.
- Timeseries distributions (see [#4642](https://github.com/pymc-devs/pymc/issues/4642))
- Mixture distributions (see [#4781](https://github.com/pymc-devs/pymc/issues/4781))
- Cholesky distributions (see WIP PR [#4784](https://github.com/pymc-devs/pymc/pull/4784))
- Variational inference submodule (see WIP PR [#4582](https://github.com/pymc-devs/pymc/pull/4582))
- Elliptical slice sampling (see [#5137](https://github.com/pymc-devs/pymc/issues/5137))
- `BaseStochasticGradient` (see [#5138](https://github.com/pymc-devs/pymc/issues/5138))
- `pm.sample_posterior_predictive_w` (see [#4807](https://github.com/pymc-devs/pymc/issues/4807))
- Partially observed Multivariate distributions (see [#5260](https://github.com/pymc-devs/pymc/issues/5260))

Also check out the [milestones](https://github.com/pymc-devs/pymc/milestones) for a potentially more complete list.

### Unexpected breaking changes (action needed)
+ New API is not available in `v3.11.5`.
+ Old API does not work in `v4.0.0`.

All of the above apply to:

- ⚠ The library is now named, installed and imported as "pymc". For example: `pip install pymc`.
- ⚠ Theano-PyMC has been replaced with Aesara, so all external references to `theano`, `tt`, and `pymc3.theanof` need to be replaced with `aesara`, `at`, and `pymc.aesaraf` (see [4471](https://github.com/pymc-devs/pymc/pull/4471)).
- `pm.Distribution(...).logp(x)` is now `pm.logp(pm.Distribution(...), x)`.
- `pm.Distribution(...).logcdf(x)` is now `pm.logcdf(pm.Distribution(...), x)`.
- `pm.Distribution(...).random(size=x)` is now `pm.draw(pm.Distribution(...), draws=x)`.
- `pm.draw_values(...)` and `pm.generate_samples(...)` were removed.
- `pm.fast_sample_posterior_predictive` was removed.
- `pm.sample_prior_predictive`, `pm.sample_posterior_predictive` and `pm.sample_posterior_predictive_w` now return an `InferenceData` object by default, instead of a dictionary (see [#5073](https://github.com/pymc-devs/pymc/pull/5073)).
- `pm.sample_prior_predictive` no longer returns transformed variable values by default. Pass them by name in `var_names` if you want to obtain these draws (see [4769](https://github.com/pymc-devs/pymc/pull/4769)).
- `pm.sample(trace=...)` no longer accepts `MultiTrace` or `len(.) > 0` traces ([see 5019#](https://github.com/pymc-devs/pymc/pull/5019)).
- `logpt`, `logpt_sum`, `logp_elemwiset` and `nojac` variations were removed. Use `Model.logpt(jacobian=True/False, sum=True/False)` instead.
- `dlogp_nojact` and `d2logp_nojact` were removed. Use `Model.dlogpt` and `d2logpt` with `jacobian=False` instead.
- `logp`, `dlogp`, and `d2logp` and `nojac` variations were removed. Use `Model.compile_logp`, `compile_dlgop` and `compile_d2logp` with `jacobian` keyword instead.
- `model.makefn` is now called `Model.compile_fn`, and `model.fn` was removed.
- Methods starting with `fast_*`, such as `Model.fast_logp`, were removed. Same applies to `PointFunc` classes
- The GLM submodule was removed, please use [Bambi](https://bambinos.github.io/bambi/) instead.
- `pm.Bound` interface no longer accepts a callable class as argument, instead it requires an instantiated distribution (created via the `.dist()` API) to be passed as an argument. In addition, Bound no longer returns a class instance but works as a normal PyMC distribution. Finally, it is no longer possible to do predictive random sampling from Bounded variables. Please, consult the new documentation for details on how to use Bounded variables (see [4815](https://github.com/pymc-devs/pymc/pull/4815)).
- `Model(model=...)` kwarg was removed
- `Model(theano_config=...)` kwarg was removed
- `Model.size` property was removed (use `Model.ndim` instead).
- `dims` and `coords` handling:
    - `Model.RV_dims` and `Model.coords` are now read-only properties. To modify the `coords` dictionary use `Model.add_coord`.
    - `dims` or coordinate values that are `None` will be auto-completed (see [#4625](https://github.com/pymc-devs/pymc/pull/4625)).
    - Coordinate values passed to `Model.add_coord` are always converted to tuples (see [#5061](https://github.com/pymc-devs/pymc/pull/5061)).
- `Model.update_start_values(...)` was removed. Initial values can be set in the `Model.initial_values` dictionary directly.
- Test values can no longer be set through `pm.Distribution(testval=...)` and must be assigned manually.
- `Transform.forward` and `Transform.backward` signatures changed.
- `pm.DensityDist` no longer accepts the `logp` as its first position argument. It is now an optional keyword argument. If you pass a callable as the first positional argument, a `TypeError` will be raised (see [5026](https://github.com/pymc-devs/pymc/pull/5026)).
- `pm.DensityDist` now accepts distribution parameters as positional arguments. Passing them as a dictionary in the `observed` keyword argument is no longer supported and will raise an error (see [5026](https://github.com/pymc-devs/pymc/pull/5026)).
- The signature of the `logp` and `random` functions that can be passed into a `pm.DensityDist` has been changed (see [5026](https://github.com/pymc-devs/pymc/pull/5026)).
- Changes to the Gaussian process (`gp`) submodule:
  - The `gp.prior(..., shape=...)` kwarg was renamed to `size`.
  - Multiple methods including `gp.prior` now require explicit kwargs.
- Changes to the BART implementation:
  - A BART variable can be combined with other random variables. The `inv_link` argument has been removed (see [4914](https://github.com/pymc-devs/pymc3/pull/4914)).
  - Moved BART to its own module (see [5058](https://github.com/pymc-devs/pymc3/pull/5058)).
- Changes to the Gaussian Process (GP) submodule (see [5055](https://github.com/pymc-devs/pymc/pull/5055)):
  - For all implementations, `gp.Latent`, `gp.Marginal` etc., `cov_func` and `mean_func` are required kwargs.
  - In Windows test conda environment the `mkl` version is fixed to verison 2020.4, and `mkl-service` is fixed to `2.3.0`.  This was required for `gp.MarginalKron` to function properly.
  - `gp.MvStudentT` uses rotated samples from `StudentT` directly now, instead of sampling from `pm.Chi2` and then from `pm.Normal`.
  - The "jitter" parameter, or the diagonal noise term added to Gram matrices such that the Cholesky is numerically stable, is now exposed to the user instead of hard-coded.  See the function `gp.util.stabilize`.
  - The `is_observed` arguement for `gp.Marginal*` implementations has been deprecated.
  - In the gp.utils file, the `kmeans_inducing_points` function now passes through `kmeans_kwargs` to scipy's k-means function.
  - The function `replace_with_values` function has been added to `gp.utils`.
  - `MarginalSparse` has been renamed `MarginalApprox`.
- ...

### Expected breaks
+ New API was already available in `v3`.
+ Old API had deprecation warnings since at least `3.11.0` (2021-01).
+ Old API stops working in `v4` (preferably with informative errors).

All of the above apply to:

- `pm.sample(return_inferencedata=True)` is now the default (see [#4744](https://github.com/pymc-devs/pymc/pull/4744)).
- ArviZ `plots` and `stats` *wrappers* were removed. The functions are now just available by their original names (see [#4549](https://github.com/pymc-devs/pymc/pull/4471) and `3.11.2` release notes).
- `pm.sample_posterior_predictive(vars=...)` kwarg was removed in favor of `var_names` (see [#4343](https://github.com/pymc-devs/pymc/pull/4343)).
- `ElemwiseCategorical` step method was removed (see [#4701](https://github.com/pymc-devs/pymc/pull/4701))
- `LKJCholeskyCov` `compute_corr` keyword argument is now set to `True` by default (see[#5382](https://github.com/pymc-devs/pymc/pull/5382))

### Ongoing deprecations
- Old API still works in `v4` and has a deprecation warning.
- Prefereably the new API should be available in `v3` already.

This includes API changes we did not warn about since at least `3.11.0` (2021-01).

- Setting initial values through `pm.Distribution(testval=...)` is now `pm.Distribution(initval=...)`.


### New features
- The length of `dims` in the model is now tracked symbolically through `Model.dim_lengths` (see [#4625](https://github.com/pymc-devs/pymc/pull/4625)).
- The `CAR` distribution has been added to allow for use of conditional autoregressions which often are used in spatial and network models.
- The dimensionality of model variables can now be parametrized through either of `shape`, `dims` or `size` (see [#4696](https://github.com/pymc-devs/pymc/pull/4696)):
  - With `shape` the length of dimensions must be given numerically or as scalar Aesara `Variables`. Numeric entries in `shape` restrict the model variable to the exact length and re-sizing is no longer possible.
  - `dims` keeps model variables re-sizeable (for example through `pm.Data`) and leads to well defined coordinates in `InferenceData` objects.
  - The `size` kwarg behaves like it does in Aesara/NumPy. For univariate RVs it is the same as `shape`, but for multivariate RVs it depends on how the RV implements broadcasting to dimensionality greater than `RVOp.ndim_supp`.
  - An `Ellipsis` (`...`) in the last position of `shape` or `dims` can be used as short-hand notation for implied dimensions.
- Added a `logcdf` implementation for the Kumaraswamy distribution (see [#4706](https://github.com/pymc-devs/pymc/pull/4706)).
- The `OrderedMultinomial` distribution has been added for use on ordinal data which are _aggregated_ by trial, like multinomial observations, whereas `OrderedLogistic` only accepts ordinal data in a _disaggregated_ format, like categorical
  observations (see [#4773](https://github.com/pymc-devs/pymc/pull/4773)).
- The `Polya-Gamma` distribution has been added (see [#4531](https://github.com/pymc-devs/pymc/pull/4531)). To make use of this distribution, the [`polyagamma>=1.3.1`](https://pypi.org/project/polyagamma/) library must be installed and available in the user's environment.
- A small change to the mass matrix tuning methods jitter+adapt_diag (the default) and adapt_diag improves performance early on during tuning for some models. [#5004](https://github.com/pymc-devs/pymc/pull/5004)
- New experimental mass matrix tuning method jitter+adapt_diag_grad. [#5004](https://github.com/pymc-devs/pymc/pull/5004)
- `pm.DensityDist` can now accept an optional `logcdf` keyword argument to pass in a function to compute the cummulative density function of the distribution (see [5026](https://github.com/pymc-devs/pymc/pull/5026)).
- `pm.DensityDist` can now accept an optional `get_moment` keyword argument to pass in a function to compute the moment of the distribution (see [5026](https://github.com/pymc-devs/pymc/pull/5026)).
- New features for BART:
  - Added partial dependence plots and individual conditional expectation plots [5091](https://github.com/pymc-devs/pymc3/pull/5091).
  - Modify how particle weights are computed. This improves accuracy of the modeled function (see [5177](https://github.com/pymc-devs/pymc3/pull/5177)).
  - Improve sampling, increase default number of particles [5229](https://github.com/pymc-devs/pymc3/pull/5229).
- The new `pm.find_constrained_prior` function can be used to find optimized prior parameters of a distribution under some
  constraints (e.g lower and upper bound). See [#5231](https://github.com/pymc-devs/pymc/pull/5231).
- New features for `pm.Data` containers:
  - With `pm.Data(..., mutable=True/False)`, or by using `pm.MutableData` vs. `pm.ConstantData` one can now create `TensorConstant` data variables. They can be more performant and compatible in situations where a variable doesn't need to be changed via `pm.set_data()`. See [#5295](https://github.com/pymc-devs/pymc/pull/5295).
  - New named dimensions can be introduced to the model via `pm.Data(..., dims=...)`. For mutable data variables (see above) the lengths of these dimensions are symbolic, so they can be re-sized via `pm.set_data()`.
  - `pm.Data` now passes additional kwargs to `aesara.shared`/`at.as_tensor`. [#5098](https://github.com/pymc-devs/pymc/pull/5098).
- Univariate censored distributions are now available via `pm.Censored`. [#5169](https://github.com/pymc-devs/pymc/pull/5169)
- Nested models now inherit the parent model's coordinates. [#5344](https://github.com/pymc-devs/pymc/pull/5344)
- `softmax` and `log_softmax` functions added to `math` module (see [#5279](https://github.com/pymc-devs/pymc/pull/5279)).
- ...

## Documentation
- Switched to the [pydata-sphinx-theme](https://pydata-sphinx-theme.readthedocs.io/en/latest/)
- Updated our documentation tooling to use [MyST](https://myst-parser.readthedocs.io/en/latest/), [MyST-NB](https://myst-nb.readthedocs.io/en/latest/), sphinx-design, notfound.extension,
  sphinx-copybutton and sphinx-remove-toctrees.
- Separated the builds of the example notebooks and of the versioned docs.
- Restructured the documentation to facilitate learning paths
- Updated API docs to document objects at the path users should use to import them

### Internal changes
- ⚠ PyMC now requires Scipy version `>= 1.4.1` (see [4857](https://github.com/pymc-devs/pymc/pull/4857)).
- Removed float128 dtype support (see [#4514](https://github.com/pymc-devs/pymc/pull/4514)).
- Logp method of `Uniform` and `DiscreteUniform` no longer depends on `pymc.distributions.dist_math.bound` for proper evaluation (see [#4541](https://github.com/pymc-devs/pymc/pull/4541)).
- We now include `cloudpickle` as a required dependency, and no longer depend on `dill` (see [#4858](https://github.com/pymc-devs/pymc/pull/4858)).
- The `incomplete_beta` function in `pymc.distributions.dist_math` was replaced by `aesara.tensor.betainc` (see [4857](https://github.com/pymc-devs/pymc/pull/4857)).
- `math.log1mexp` and `math.log1mexp_numpy` will expect negative inputs in the future. A `FutureWarning` is now raised unless `negative_input=True` is set (see [#4860](https://github.com/pymc-devs/pymc/pull/4860)).
- Changed name of `Lognormal` distribution to `LogNormal` to harmonize CamelCase usage for distribution names.
- Attempt to iterate over MultiTrace will raise NotImplementedError.
- Removed silent normalisation of `p` parameters in Categorical and Multinomial distributions (see [#5370](https://github.com/pymc-devs/pymc/pull/5370)).
- ...


## PyMC 3.11.2 (14 March 2021)

### New Features
+ `pm.math.cartesian` can now handle inputs that are themselves >1D (see [#4482](https://github.com/pymc-devs/pymc/pull/4482)).
+ Statistics and plotting functions that were removed in `3.11.0` were brought back, albeit with deprecation warnings if an old naming scheme is used (see [#4536](https://github.com/pymc-devs/pymc/pull/4536)). In order to future proof your code, rename these function calls:
  + `pm.traceplot` → `pm.plot_trace`
  + `pm.compareplot` → `pm.plot_compare` (here you might need to rename some columns in the input according to the [`arviz.plot_compare` documentation](https://arviz-devs.github.io/arviz/api/generated/arviz.plot_compare.html))
  + `pm.autocorrplot` → `pm.plot_autocorr`
  + `pm.forestplot` → `pm.plot_forest`
  + `pm.kdeplot` → `pm.plot_kde`
  + `pm.energyplot` → `pm.plot_energy`
  + `pm.densityplot` → `pm.plot_density`
  + `pm.pairplot` → `pm.plot_pair`

### Maintenance
- ⚠ Our memoization mechanism wasn't robust against hash collisions ([#4506](https://github.com/pymc-devs/pymc/issues/4506)), sometimes resulting in incorrect values in, for example, posterior predictives. The `pymc.memoize` module was removed and replaced with `cachetools`.  The `hashable` function and `WithMemoization` class were moved to `pymc.util` (see [#4525](https://github.com/pymc-devs/pymc/pull/4525)).
- `pm.make_shared_replacements` now retains broadcasting information which fixes issues with Metropolis samplers (see [#4492](https://github.com/pymc-devs/pymc/pull/4492)).

**Release manager** for 3.11.2: Michael Osthege ([@michaelosthege](https://github.com/michaelosthege))

## PyMC 3.11.1 (12 February 2021)

### New Features
- Automatic imputations now also work with `ndarray` data, not just `pd.Series` or `pd.DataFrame` (see[#4439](https://github.com/pymc-devs/pymc/pull/4439)).
- `pymc.sampling_jax.sample_numpyro_nuts` now returns samples from transformed random variables, rather than from the unconstrained representation (see [#4427](https://github.com/pymc-devs/pymc/pull/4427)).

### Maintenance
- We upgraded to `Theano-PyMC v1.1.2` which [includes bugfixes](https://github.com/pymc-devs/aesara/compare/rel-1.1.0...rel-1.1.2) for...
  - ⚠ a problem with `tt.switch` that affected the behavior of several distributions, including at least the following special cases (see [#4448](https://github.com/pymc-devs/pymc/pull/4448))
    1.  `Bernoulli` when all the observed values were the same (e.g., `[0, 0, 0, 0, 0]`).
    2.  `TruncatedNormal` when `sigma` was constant and `mu` was being automatically broadcasted to match the shape of observations.
  - Warning floods and compiledir locking (see [#4444](https://github.com/pymc-devs/pymc/pull/4444))
- `math.log1mexp_numpy` no longer raises RuntimeWarning when given very small inputs. These were commonly observed during NUTS sampling (see [#4428](https://github.com/pymc-devs/pymc/pull/4428)).
- `ScalarSharedVariable` can now be used as an input to other RVs directly (see [#4445](https://github.com/pymc-devs/pymc/pull/4445)).
- `pm.sample` and `pm.find_MAP` no longer change the `start` argument (see [#4458](https://github.com/pymc-devs/pymc/pull/4458)).
- Fixed `Dirichlet.logp` method to work with unit batch or event shapes (see [#4454](https://github.com/pymc-devs/pymc/pull/4454)).
- Bugfix in logp and logcdf methods of `Triangular` distribution (see [#4470](https://github.com/pymc-devs/pymc/pull/4470)).

**Release manager** for 3.11.1: Michael Osthege ([@michaelosthege](https://github.com/michaelosthege))

## PyMC3 3.11.0 (21 January 2021)

This release breaks some APIs w.r.t. `3.10.0`. It also brings some dreadfully awaited fixes, so be sure to go through the (breaking) changes below.

### Breaking Changes
- ⚠ Many plotting and diagnostic functions that were just aliasing ArviZ functions were removed (see [4397](https://github.com/pymc-devs/pymc/pull/4397/files#)). This includes `pm.summary`, `pm.traceplot`, `pm.ess` and many more!
- ⚠ We now depend on `Theano-PyMC` version `1.1.0` exactly (see [#4405](https://github.com/pymc-devs/pymc/pull/4405)). Major refactorings were done in `Theano-PyMC` 1.1.0. If you implement custom `Op`s or interact with Theano in any way yourself, make sure to read the [Theano-PyMC 1.1.0 release notes](https://github.com/pymc-devs/Theano-PyMC/releases/tag/rel-1.1.0).
- ⚠ Python 3.6 support was dropped (by no longer testing) and Python 3.9 was added (see [#4332](https://github.com/pymc-devs/pymc/pull/4332)).
- ⚠ Changed shape behavior: __No longer collapse length 1 vector shape into scalars.__ (see [#4206](https://github.com/pymc-devs/pymc/issue/4206) and [#4214](https://github.com/pymc-devs/pymc/pull/4214))
  - __Applies to random variables and also the `.random(size=...)` kwarg!__
  - To create scalar variables you must now use `shape=None` or `shape=()`.
  - __`shape=(1,)` and `shape=1` now become vectors.__ Previously they were collapsed into scalars
  - 0-length dimensions are now ruled illegal for random variables and raise a `ValueError`.
- In `sample_prior_predictive` the `vars` kwarg was removed in favor of `var_names` (see [#4327](https://github.com/pymc-devs/pymc/pull/4327)).
- Removed `theanof.set_theano_config` because it illegally changed Theano's internal state (see [#4329](https://github.com/pymc-devs/pymc/pull/4329)).


### New Features
- Option to set `check_bounds=False` when instantiating `pymc.Model()`. This turns off bounds checks that ensure that input parameters of distributions are valid. For correctly specified models, this is unneccessary as all parameters get automatically transformed so that all values are valid. Turning this off should lead to faster sampling (see [#4377](https://github.com/pymc-devs/pymc/pull/4377)).
- `OrderedProbit` distribution added (see [#4232](https://github.com/pymc-devs/pymc/pull/4232)).
- `plot_posterior_predictive_glm` now works with `arviz.InferenceData` as well (see [#4234](https://github.com/pymc-devs/pymc/pull/4234))
- Add `logcdf` method to all univariate discrete distributions (see [#4387](https://github.com/pymc-devs/pymc/pull/4387)).
- Add `random` method to `MvGaussianRandomWalk` (see [#4388](https://github.com/pymc-devs/pymc/pull/4388))
- `AsymmetricLaplace` distribution added (see [#4392](https://github.com/pymc-devs/pymc/pull/4392)).
- `DirichletMultinomial` distribution added (see [#4373](https://github.com/pymc-devs/pymc/pull/4373)).
- Added a new `predict` method to `BART` to compute out of sample predictions (see [#4310](https://github.com/pymc-devs/pymc/pull/4310)).

### Maintenance
- Fixed bug whereby partial traces returns after keyboard interrupt during parallel sampling had fewer draws than would've been available [#4318](https://github.com/pymc-devs/pymc/pull/4318)
- Make `sample_shape` same across all contexts in `draw_values` (see [#4305](https://github.com/pymc-devs/pymc/pull/4305)).
- The notebook gallery has been moved to https://github.com/pymc-devs/pymc-examples (see [#4348](https://github.com/pymc-devs/pymc/pull/4348)).
- `math.logsumexp` now matches `scipy.special.logsumexp` when arrays contain infinite values (see [#4360](https://github.com/pymc-devs/pymc/pull/4360)).
- Fixed mathematical formulation in `MvStudentT` random method. (see [#4359](https://github.com/pymc-devs/pymc/pull/4359))
- Fix issue in `logp` method of `HyperGeometric`. It now returns `-inf` for invalid parameters (see [4367](https://github.com/pymc-devs/pymc/pull/4367))
- Fixed `MatrixNormal` random method to work with parameters as random variables. (see [#4368](https://github.com/pymc-devs/pymc/pull/4368))
- Update the `logcdf` method of several continuous distributions to return -inf for invalid parameters and values, and raise an informative error when multiple values cannot be evaluated in a single call. (see [4393](https://github.com/pymc-devs/pymc/pull/4393) and [#4421](https://github.com/pymc-devs/pymc/pull/4421))
- Improve numerical stability in `logp` and `logcdf` methods of `ExGaussian` (see [#4407](https://github.com/pymc-devs/pymc/pull/4407))
- Issue UserWarning when doing prior or posterior predictive sampling with models containing Potential factors (see [#4419](https://github.com/pymc-devs/pymc/pull/4419))
- Dirichlet distribution's `random` method is now optimized and gives outputs in correct shape (see [#4416](https://github.com/pymc-devs/pymc/pull/4407))
- Attempting to sample a named model with SMC will now raise a `NotImplementedError`. (see [#4365](https://github.com/pymc-devs/pymc/pull/4365))

**Release manager** for 3.11.0: Eelke Spaak ([@Spaak](https://github.com/Spaak))

## PyMC3 3.10.0 (7 December 2020)

This is a major release with many exciting new features. The biggest change is that we now rely on our own fork of [Theano-PyMC](https://github.com/pymc-devs/Theano-PyMC). This is in line with our [big announcement about our commitment to PyMC3 and Theano](https://pymc-devs.medium.com/the-future-of-pymc-or-theano-is-dead-long-live-theano-d8005f8a0e9b).

When upgrading, make sure that `Theano-PyMC` and not `Theano` are installed (the imports remain unchanged, however). If not, you can uninstall `Theano`:
```
conda remove theano
```

And to install:
```
conda install -c conda-forge theano-pymc
```

Or, if you are using pip (not recommended):
```
pip uninstall theano
```
And to install:
```
pip install theano-pymc
```

This new version of `Theano-PyMC` comes with an experimental JAX backend which, when combined with the new and experimental JAX samplers in PyMC3, can greatly speed up sampling in your model. As this is still very new, please do not use it in production yet but do test it out and let us know if anything breaks and what results you are seeing, especially speed-wise.

### New features
- New experimental JAX samplers in `pymc.sample_jax` (see [notebook](https://docs.pymc.io/notebooks/GLM-hierarchical-jax.html) and [#4247](https://github.com/pymc-devs/pymc/pull/4247)). Requires JAX and either TFP or numpyro.
- Add MLDA, a new stepper for multilevel sampling. MLDA can be used when a hierarchy of approximate posteriors of varying accuracy is available, offering improved sampling efficiency especially in high-dimensional problems and/or where gradients are not available (see [#3926](https://github.com/pymc-devs/pymc/pull/3926))
- Add Bayesian Additive Regression Trees (BARTs) [#4183](https://github.com/pymc-devs/pymc/pull/4183))
- Added `pymc.gp.cov.Circular` kernel for Gaussian Processes on circular domains, e.g. the unit circle (see [#4082](https://github.com/pymc-devs/pymc/pull/4082)).
- Added a new `MixtureSameFamily` distribution to handle mixtures of arbitrary dimensions in vectorized form for improved speed (see [#4185](https://github.com/pymc-devs/pymc/issues/4185)).
- `sample_posterior_predictive_w` can now feed on `xarray.Dataset` - e.g. from `InferenceData.posterior`. (see [#4042](https://github.com/pymc-devs/pymc/pull/4042))
- Change SMC metropolis kernel to independent metropolis kernel [#4115](https://github.com/pymc-devs/pymc/pull/4115))
- Add alternative parametrization to NegativeBinomial distribution in terms of n and p (see [#4126](https://github.com/pymc-devs/pymc/issues/4126))
- Added semantically meaningful `str` representations to PyMC3 objects for console, notebook, and GraphViz use (see [#4076](https://github.com/pymc-devs/pymc/pull/4076), [#4065](https://github.com/pymc-devs/pymc/pull/4065), [#4159](https://github.com/pymc-devs/pymc/pull/4159), [#4217](https://github.com/pymc-devs/pymc/pull/4217), [#4243](https://github.com/pymc-devs/pymc/pull/4243), and [#4260](https://github.com/pymc-devs/pymc/pull/4260)).
- Add Discrete HyperGeometric Distribution (see [#4249](https://github.com/pymc-devs/pymc/pull/#4249))

### Maintenance
- Switch the dependency of Theano to our own fork, [Theano-PyMC](https://github.com/pymc-devs/Theano-PyMC).
- Removed non-NDArray (Text, SQLite, HDF5) backends and associated tests.
- Use dill to serialize user defined logp functions in `DensityDist`. The previous serialization code fails if it is used in notebooks on Windows and Mac. `dill` is now a required dependency. (see [#3844](https://github.com/pymc-devs/pymc/issues/3844)).
- Fixed numerical instability in ExGaussian's logp by preventing `logpow` from returning `-inf` (see [#4050](https://github.com/pymc-devs/pymc/pull/4050)).
- Numerically improved stickbreaking transformation - e.g. for the `Dirichlet` distribution. [#4129](https://github.com/pymc-devs/pymc/pull/4129)
- Enabled the `Multinomial` distribution to handle batch sizes that have more than 2 dimensions. [#4169](https://github.com/pymc-devs/pymc/pull/4169)
- Test model logp before starting any MCMC chains (see [#4211](https://github.com/pymc-devs/pymc/pull/4211))
- Fix bug in `model.check_test_point` that caused the `test_point` argument to be ignored. (see [PR #4211](https://github.com/pymc-devs/pymc/pull/4211#issuecomment-727142721))
- Refactored MvNormal.random method with better handling of sample, batch and event shapes. [#4207](https://github.com/pymc-devs/pymc/pull/4207)
- The `InverseGamma` distribution now implements a `logcdf`. [#3944](https://github.com/pymc-devs/pymc/pull/3944)
- Make starting jitter methods for nuts sampling more robust by resampling values that lead to non-finite probabilities. A new optional argument `jitter-max-retries` can be passed to `pm.sample()` and `pm.init_nuts()` to control the maximum number of retries per chain. [4298](https://github.com/pymc-devs/pymc/pull/4298)

### Documentation
- Added a new notebook demonstrating how to incorporate sampling from a conjugate Dirichlet-multinomial posterior density in conjunction with other step methods (see [#4199](https://github.com/pymc-devs/pymc/pull/4199)).
- Mentioned the way to do any random walk with `theano.tensor.cumsum()` in `GaussianRandomWalk` docstrings (see [#4048](https://github.com/pymc-devs/pymc/pull/4048)).

**Release manager** for 3.10.0: Eelke Spaak ([@Spaak](https://github.com/Spaak))

## PyMC3 3.9.3 (11 August 2020)

### New features
- Introduce optional arguments to `pm.sample`: `mp_ctx` to control how the processes for parallel sampling are started, and `pickle_backend` to specify which library is used to pickle models in parallel sampling when the multiprocessing context is not of type `fork` (see [#3991](https://github.com/pymc-devs/pymc/pull/3991)).
- Add sampler stats `process_time_diff`, `perf_counter_diff` and `perf_counter_start`, that record wall and CPU times for each NUTS and HMC sample (see [ #3986](https://github.com/pymc-devs/pymc/pull/3986)).
- Extend `keep_size` argument handling for `sample_posterior_predictive` and `fast_sample_posterior_predictive`, to work on ArviZ `InferenceData` and xarray `Dataset` input values (see [PR #4006](https://github.com/pymc-devs/pymc/pull/4006) and issue [#4004](https://github.com/pymc-devs/pymc/issues/4004)).
- SMC-ABC: add the Wasserstein and energy distance functions. Refactor API, the distance, sum_stats and epsilon arguments are now passed `pm.Simulator` instead of `pm.sample_smc`. Add random method to `pm.Simulator`. Add option to save the simulated data. Improved LaTeX representation [#3996](https://github.com/pymc-devs/pymc/pull/3996).
- SMC-ABC: Allow use of potentials by adding them to the prior term. [#4016](https://github.com/pymc-devs/pymc/pull/4016).

### Maintenance
- Fix an error on Windows and Mac where error message from unpickling models did not show up in the notebook, or where sampling froze when a worker process crashed (see [#3991](https://github.com/pymc-devs/pymc/pull/3991)).
- Require Theano >= 1.0.5 (see [#4032](https://github.com/pymc-devs/pymc/pull/4032)).

### Documentation
- Notebook on [multilevel modeling](https://docs.pymc.io/notebooks/multilevel_modeling.html) has been rewritten to showcase ArviZ and xarray usage for inference result analysis (see [#3963](https://github.com/pymc-devs/pymc/pull/3963)).

_NB: The `docs/*` folder is still removed from the tarball due to an upload size limit on PyPi._

**Release manager** for 3.9.3: Kyle Beauchamp ([@kyleabeauchamp](https://github.com/kyleabeauchamp))

## PyMC3 3.9.2 (24 June 2020)

### Maintenance
- Warning added in GP module when `input_dim` is lower than the number of columns in `X` to compute the covariance function (see [#3974](https://github.com/pymc-devs/pymc/pull/3974)).
- Pass the `tune` argument from `sample` when using `advi+adapt_diag_grad` (see issue [#3965](https://github.com/pymc-devs/pymc/issues/3965), fixed by [#3979](https://github.com/pymc-devs/pymc/pull/3979)).
- Add simple test case for new coords and dims feature in `pm.Model` (see [#3977](https://github.com/pymc-devs/pymc/pull/3977)).
- Require ArviZ >= 0.9.0 (see [#3977](https://github.com/pymc-devs/pymc/pull/3977)).
- Fixed issue [#3962](https://github.com/pymc-devs/pymc/issues/3962) by making a change in the `_random()` method of `GaussianRandomWalk` class (see PR [#3985](https://github.com/pymc-devs/pymc/pull/3985)). Further testing revealed a new issue which is being tracked by [#4010](https://github.com/pymc-devs/pymc/issues/4010).

_NB: The `docs/*` folder is still removed from the tarball due to an upload size limit on PyPi._

**Release manager** for 3.9.2: Alex Andorra ([@AlexAndorra](https://github.com/AlexAndorra))

## PyMC3 3.9.1 (16 June 2020)
The `v3.9.0` upload to PyPI didn't include a tarball, which is fixed in this release.
Though we had to temporarily remove the `docs/*` folder from the tarball due to a size limit.

**Release manager** for 3.9.1: Michael Osthege ([@michaelosthege](https://github.com/michaelosthege))

## PyMC3 3.9.0 (16 June 2020)

### New features
- Use [fastprogress](https://github.com/fastai/fastprogress) instead of tqdm [#3693](https://github.com/pymc-devs/pymc/pull/3693).
- `DEMetropolis` can now tune both `lambda` and `scaling` parameters, but by default neither of them are tuned. See [#3743](https://github.com/pymc-devs/pymc/pull/3743) for more info.
- `DEMetropolisZ`, an improved variant of `DEMetropolis` brings better parallelization and higher efficiency with fewer chains with a slower initial convergence. This implementation is experimental. See [#3784](https://github.com/pymc-devs/pymc/pull/3784) for more info.
- Notebooks that give insight into `DEMetropolis`, `DEMetropolisZ` and the `DifferentialEquation` interface are now located in the [Tutorials/Deep Dive](https://docs.pymc.io/nb_tutorials/index.html) section.
- Add `fast_sample_posterior_predictive`, a vectorized alternative to `sample_posterior_predictive`.  This alternative is substantially faster for large models.
- GP covariance functions can now be exponentiated by a scalar. See PR [#3852](https://github.com/pymc-devs/pymc/pull/3852)
- `sample_posterior_predictive` can now feed on `xarray.Dataset` - e.g. from `InferenceData.posterior`. (see [#3846](https://github.com/pymc-devs/pymc/pull/3846))
- `SamplerReport` (`MultiTrace.report`) now has properties `n_tune`, `n_draws`, `t_sampling` for increased convenience (see [#3827](https://github.com/pymc-devs/pymc/pull/3827))
- `pm.sample(..., return_inferencedata=True)` can now directly return the trace as `arviz.InferenceData` (see [#3911](https://github.com/pymc-devs/pymc/pull/3911))
- `pm.sample` now has support for adapting dense mass matrix using `QuadPotentialFullAdapt` (see [#3596](https://github.com/pymc-devs/pymc/pull/3596), [#3705](https://github.com/pymc-devs/pymc/pull/3705), [#3858](https://github.com/pymc-devs/pymc/pull/3858), and [#3893](https://github.com/pymc-devs/pymc/pull/3893)). Use `init="adapt_full"` or `init="jitter+adapt_full"` to use.
- `Moyal` distribution added (see [#3870](https://github.com/pymc-devs/pymc/pull/3870)).
- `pm.LKJCholeskyCov` now automatically computes and returns the unpacked Cholesky decomposition, the correlations and the standard deviations of the covariance matrix (see [#3881](https://github.com/pymc-devs/pymc/pull/3881)).
- `pm.Data` container can now be used for index variables, i.e with integer data and not only floats (issue [#3813](https://github.com/pymc-devs/pymc/issues/3813), fixed by [#3925](https://github.com/pymc-devs/pymc/pull/3925)).
- `pm.Data` container can now be used as input for other random variables (issue [#3842](https://github.com/pymc-devs/pymc/issues/3842), fixed by [#3925](https://github.com/pymc-devs/pymc/pull/3925)).
- Allow users to specify coordinates and dimension names instead of numerical shapes when specifying a model. This makes interoperability with ArviZ easier. ([see #3551](https://github.com/pymc-devs/pymc/pull/3551))
- Plots and Stats API sections now link to ArviZ documentation [#3927](https://github.com/pymc-devs/pymc/pull/3927)
- Add `SamplerReport` with properties `n_draws`, `t_sampling` and `n_tune` to SMC. `n_tune` is always 0 [#3931](https://github.com/pymc-devs/pymc/issues/3931).
- SMC-ABC: add option to define summary statistics, allow to sample from more complex models, remove redundant distances [#3940](https://github.com/pymc-devs/pymc/issues/3940)

### Maintenance
- Tuning results no longer leak into sequentially sampled `Metropolis` chains (see #3733 and #3796).
- We'll deprecate the `Text` and `SQLite` backends and the `save_trace`/`load_trace` functions, since this is now done with ArviZ. (see [#3902](https://github.com/pymc-devs/pymc/pull/3902))
- ArviZ `v0.8.3` is now the minimum required version
- In named models, `pm.Data` objects now get model-relative names (see [#3843](https://github.com/pymc-devs/pymc/pull/3843)).
- `pm.sample` now takes 1000 draws and 1000 tuning samples by default, instead of 500 previously (see [#3855](https://github.com/pymc-devs/pymc/pull/3855)).
- Moved argument division out of `NegativeBinomial` `random` method. Fixes [#3864](https://github.com/pymc-devs/pymc/issues/3864) in the style of [#3509](https://github.com/pymc-devs/pymc/pull/3509).
- The Dirichlet distribution now raises a ValueError when it's initialized with <= 0 values (see [#3853](https://github.com/pymc-devs/pymc/pull/3853)).
- Dtype bugfix in `MvNormal` and `MvStudentT` (see [3836](https://github.com/pymc-devs/pymc/pull/3836)).
- End of sampling report now uses `arviz.InferenceData` internally and avoids storing
  pointwise log likelihood (see [#3883](https://github.com/pymc-devs/pymc/pull/3883)).
- The multiprocessing start method on MacOS is now set to "forkserver", to avoid crashes (see issue [#3849](https://github.com/pymc-devs/pymc/issues/3849), solved by [#3919](https://github.com/pymc-devs/pymc/pull/3919)).
- The AR1 logp now uses the precision of the whole AR1 process instead of just the innovation precision (see issue [#3892](https://github.com/pymc-devs/pymc/issues/3892), fixed by [#3899](https://github.com/pymc-devs/pymc/pull/3899)).
- Forced the `Beta` distribution's `random` method to generate samples that are in the open interval $(0, 1)$, i.e. no value can be equal to zero or equal to one (issue [#3898](https://github.com/pymc-devs/pymc/issues/3898) fixed by [#3924](https://github.com/pymc-devs/pymc/pull/3924)).
- Fixed an issue that happened on Windows, that was introduced by the clipped beta distribution rvs function ([#3924](https://github.com/pymc-devs/pymc/pull/3924)). Windows does not support the `float128` dtype, but we had assumed that it had to be available. The solution was to only support `float128` on Linux and Darwin systems (see issue [#3929](https://github.com/pymc-devs/pymc/issues/3849) fixed by [#3930](https://github.com/pymc-devs/pymc/pull/3930)).

### Deprecations
- Remove `sample_ppc` and `sample_ppc_w` that were deprecated in 3.6.
- Deprecated `sd` has been replaced by `sigma` (already in version 3.7) in continuous, mixed and timeseries distributions and now raises `DeprecationWarning` when `sd` is used. (see [#3837](https://github.com/pymc-devs/pymc/pull/3837) and [#3688](https://github.com/pymc-devs/pymc/issues/3688)).
- We'll deprecate the `Text` and `SQLite` backends and the `save_trace`/`load_trace` functions, since this is now done with ArviZ. (see [#3902](https://github.com/pymc-devs/pymc/pull/3902))
- Dropped some deprecated kwargs and functions (see [#3906](https://github.com/pymc-devs/pymc/pull/3906))
- Dropped the outdated 'nuts' initialization method for `pm.sample` (see [#3863](https://github.com/pymc-devs/pymc/pull/3863)).

**Release manager** for 3.9.0: Michael Osthege ([@michaelosthege](https://github.com/michaelosthege))

## PyMC3 3.8 (November 29 2019)

### New features
- Implemented robust u turn check in NUTS (similar to stan-dev/stan#2800). See PR [#3605]
- Add capabilities to do inference on parameters in a differential equation with `DifferentialEquation`. See [#3590](https://github.com/pymc-devs/pymc/pull/3590) and [#3634](https://github.com/pymc-devs/pymc/pull/3634).
- Distinguish between `Data` and `Deterministic` variables when graphing models with graphviz. PR [#3491](https://github.com/pymc-devs/pymc/pull/3491).
- Sequential Monte Carlo - Approximate Bayesian Computation step method is now available. The implementation is in an experimental stage and will be further improved.
- Added `Matern12` covariance function for Gaussian processes. This is the Matern kernel with nu=1/2.
- Progressbar reports number of divergences in real time, when available [#3547](https://github.com/pymc-devs/pymc/pull/3547).
- Sampling from variational approximation now allows for alternative trace backends [#3550].
- Infix `@` operator now works with random variables and deterministics [#3619](https://github.com/pymc-devs/pymc/pull/3619).
- [ArviZ](https://arviz-devs.github.io/arviz/) is now a requirement, and handles plotting, diagnostics, and statistical checks.
- Can use GaussianRandomWalk in sample_prior_predictive and sample_prior_predictive [#3682](https://github.com/pymc-devs/pymc/pull/3682)
- Now 11 years of S&P returns in data set[#3682](https://github.com/pymc-devs/pymc/pull/3682)

### Maintenance
- Moved math operations out of `Rice`, `TruncatedNormal`, `Triangular` and `ZeroInflatedNegativeBinomial` `random` methods. Math operations on values returned by `draw_values` might not broadcast well, and all the `size` aware broadcasting is left to `generate_samples`. Fixes [#3481](https://github.com/pymc-devs/pymc/issues/3481) and [#3508](https://github.com/pymc-devs/pymc/issues/3508)
- Parallelization of population steppers (`DEMetropolis`) is now set via the `cores` argument. ([#3559](https://github.com/pymc-devs/pymc/pull/3559))
- Fixed a bug in `Categorical.logp`. In the case of multidimensional `p`'s, the indexing was done wrong leading to incorrectly shaped tensors that consumed `O(n**2)` memory instead of `O(n)`. This fixes issue [#3535](https://github.com/pymc-devs/pymc/issues/3535)
- Fixed a defect in `OrderedLogistic.__init__` that unnecessarily increased the dimensionality of the underlying `p`. Related to issue issue [#3535](https://github.com/pymc-devs/pymc/issues/3535) but was not the true cause of it.
- SMC: stabilize covariance matrix [3573](https://github.com/pymc-devs/pymc/pull/3573)
- SMC: is no longer a step method of `pm.sample` now it should be called using `pm.sample_smc` [3579](https://github.com/pymc-devs/pymc/pull/3579)
- SMC: improve computation of the proposal scaling factor [3594](https://github.com/pymc-devs/pymc/pull/3594) and [3625](https://github.com/pymc-devs/pymc/pull/3625)
- SMC: reduce number of logp evaluations [3600](https://github.com/pymc-devs/pymc/pull/3600)
- SMC: remove `scaling` and `tune_scaling` arguments as is a better idea to always allow SMC to automatically compute the scaling factor [3625](https://github.com/pymc-devs/pymc/pull/3625)
- Now uses `multiprocessong` rather than `psutil` to count CPUs, which results in reliable core counts on Chromebooks.
- `sample_posterior_predictive` now preallocates the memory required for its output to improve memory usage. Addresses problems raised in this [discourse thread](https://discourse.pymc.io/t/memory-error-with-posterior-predictive-sample/2891/4).
- Fixed a bug in `Categorical.logp`. In the case of multidimensional `p`'s, the indexing was done wrong leading to incorrectly shaped tensors that consumed `O(n**2)` memory instead of `O(n)`. This fixes issue [#3535](https://github.com/pymc-devs/pymc/issues/3535)
- Fixed a defect in `OrderedLogistic.__init__` that unnecessarily increased the dimensionality of the underlying `p`. Related to issue issue [#3535](https://github.com/pymc-devs/pymc/issues/3535) but was not the true cause of it.
- Wrapped `DensityDist.rand` with `generate_samples` to make it aware of the distribution's shape. Added control flow attributes to still be able to behave as in earlier versions, and to control how to interpret the `size` parameter in the `random` callable signature. Fixes [3553](https://github.com/pymc-devs/pymc/issues/3553)
- Added `theano.gof.graph.Constant` to type checks done in `_draw_value` (fixes issue [3595](https://github.com/pymc-devs/pymc/issues/3595))
- `HalfNormal` did not used to work properly in `draw_values`, `sample_prior_predictive`, or `sample_posterior_predictive` (fixes issue [3686](https://github.com/pymc-devs/pymc/pull/3686))
- Random variable transforms were inadvertently left out of the API documentation. Added them. (See PR [3690](https://github.com/pymc-devs/pymc/pull/3690)).
- Refactored `pymc.model.get_named_nodes_and_relations` to use the ancestors and descendents, in a way that is consistent with `theano`'s naming convention.
- Changed the way in which `pymc.model.get_named_nodes_and_relations` computes nodes without ancestors to make it robust to changes in var_name orderings (issue [#3643](https://github.com/pymc-devs/pymc/issues/3643))

## PyMC3 3.7 (May 29 2019)

### New features

- Add data container class (`Data`) that wraps the theano SharedVariable class and let the model be aware of its inputs and outputs.
- Add function `set_data` to update variables defined as `Data`.
- `Mixture` now supports mixtures of multidimensional probability distributions, not just lists of 1D distributions.
- `GLM.from_formula` and `LinearComponent.from_formula` can extract variables from the calling scope. Customizable via the new `eval_env` argument. Fixing [#3382](https://github.com/pymc-devs/pymc/issues/3382).
- Added the `distributions.shape_utils` module with functions used to help broadcast samples drawn from distributions using the `size` keyword argument.
- Used `numpy.vectorize` in `distributions.distribution._compile_theano_function`. This enables `sample_prior_predictive` and `sample_posterior_predictive` to ask for tuples of samples instead of just integers. This fixes issue [#3422](https://github.com/pymc-devs/pymc/issues/3422).

### Maintenance

- All occurances of `sd` as a parameter name have been renamed to `sigma`. `sd` will continue to function for backwards compatibility.
- `HamiltonianMC` was ignoring certain arguments like `target_accept`, and not using the custom step size jitter function with expectation 1.
- Made `BrokenPipeError` for parallel sampling more verbose on Windows.
- Added the `broadcast_distribution_samples` function that helps broadcasting arrays of drawn samples, taking into account the requested `size` and the inferred distribution shape. This sometimes is needed by distributions that call several `rvs` separately within their `random` method, such as the `ZeroInflatedPoisson` (fixes issue [#3310](https://github.com/pymc-devs/pymc/issues/3310)).
- The `Wald`, `Kumaraswamy`, `LogNormal`, `Pareto`, `Cauchy`, `HalfCauchy`, `Weibull` and `ExGaussian` distributions `random` method used a hidden `_random` function that was written with scalars in mind. This could potentially lead to artificial correlations between random draws. Added shape guards and broadcasting of the distribution samples to prevent this (Similar to issue [#3310](https://github.com/pymc-devs/pymc/issues/3310)).
- Added a fix to allow the imputation of single missing values of observed data, which previously would fail (fixes issue [#3122](https://github.com/pymc-devs/pymc/issues/3122)).
- The `draw_values` function was too permissive with what could be grabbed from inside `point`, which lead to an error when sampling posterior predictives of variables that depended on shared variables that had changed their shape after `pm.sample()` had been called (fix issue [#3346](https://github.com/pymc-devs/pymc/issues/3346)).
- `draw_values` now adds the theano graph descendants of `TensorConstant` or `SharedVariables` to the named relationship nodes stack, only if these descendants are `ObservedRV` or `MultiObservedRV` instances (fixes issue [#3354](https://github.com/pymc-devs/pymc/issues/3354)).
- Fixed bug in broadcast_distrution_samples, which did not handle correctly cases in which some samples did not have the size tuple prepended.
- Changed `MvNormal.random`'s usage of `tensordot` for Cholesky encoded covariances. This lead to wrong axis broadcasting and seemed to be the cause for issue [#3343](https://github.com/pymc-devs/pymc/issues/3343).
- Fixed defect in `Mixture.random` when multidimensional mixtures were involved. The mixture component was not preserved across all the elements of the dimensions of the mixture. This meant that the correlations across elements within a given draw of the mixture were partly broken.
- Restructured `Mixture.random` to allow better use of vectorized calls to `comp_dists.random`.
- Added tests for mixtures of multidimensional distributions to the test suite.
- Fixed incorrect usage of `broadcast_distribution_samples` in `DiscreteWeibull`.
- `Mixture`'s default dtype is now determined by `theano.config.floatX`.
- `dist_math.random_choice` now handles nd-arrays of category probabilities, and also handles sizes that are not `None`. Also removed unused `k` kwarg from `dist_math.random_choice`.
- Changed `Categorical.mode` to preserve all the dimensions of `p` except the last one, which encodes each category's probability.
- Changed initialization of `Categorical.p`. `p` is now normalized to sum to `1` inside `logp` and `random`, but not during initialization. This could hide negative values supplied to `p` as mentioned in [#2082](https://github.com/pymc-devs/pymc/issues/2082).
- `Categorical` now accepts elements of `p` equal to `0`. `logp` will return `-inf` if there are `values` that index to the zero probability categories.
- Add `sigma`, `tau`, and `sd` to signature of `NormalMixture`.
- Set default lower and upper values of -inf and inf for pm.distributions.continuous.TruncatedNormal. This avoids errors caused by their previous values of None (fixes issue [#3248](https://github.com/pymc-devs/pymc/issues/3248)).
- Converted all calls to `pm.distributions.bound._ContinuousBounded` and `pm.distributions.bound._DiscreteBounded` to use only and all positional arguments (fixes issue [#3399](https://github.com/pymc-devs/pymc/issues/3399)).
- Restructured `distributions.distribution.generate_samples` to use the `shape_utils` module. This solves issues [#3421](https://github.com/pymc-devs/pymc/issues/3421) and [#3147](https://github.com/pymc-devs/pymc/issues/3147) by using the `size` aware broadcating functions in `shape_utils`.
- Fixed the `Multinomial.random` and `Multinomial.random_` methods to make them compatible with the new `generate_samples` function. In the process, a bug of the `Multinomial.random_` shape handling was discovered and fixed.
- Fixed a defect found in `Bound.random` where the `point` dictionary was passed to `generate_samples` as an `arg` instead of in `not_broadcast_kwargs`.
- Fixed a defect found in `Bound.random_` where `total_size` could end up as a `float64` instead of being an integer if given `size=tuple()`.
- Fixed an issue in `model_graph` that caused construction of the graph of the model for rendering to hang: replaced a search over the powerset of the nodes with a breadth-first search over the nodes. Fix for [#3458](https://github.com/pymc-devs/pymc/issues/3458).
- Removed variable annotations from `model_graph` but left type hints (Fix for [#3465](https://github.com/pymc-devs/pymc/issues/3465)). This means that we support `python>=3.5.4`.
- Default `target_accept`for `HamiltonianMC` is now 0.65, as suggested in Beskos et. al. 2010 and Neal 2001.
- Fixed bug in `draw_values` that lead to intermittent errors in python3.5. This happened with some deterministic nodes that were drawn but not added to `givens`.

### Deprecations

- `nuts_kwargs` and `step_kwargs` have been deprecated in favor of using the standard `kwargs` to pass optional step method arguments.
- `SGFS` and `CSG` have been removed (Fix for [#3353](https://github.com/pymc-devs/pymc/issues/3353)). They have been moved to [pymc-experimental](https://github.com/pymc-devs/pymc-experimental).
- References to `live_plot` and corresponding notebooks have been removed.
- Function `approx_hessian` was removed, due to `numdifftools` becoming incompatible with current `scipy`. The function was already optional, only available to a user who installed `numdifftools` separately, and not hit on any common codepaths. [#3485](https://github.com/pymc-devs/pymc/pull/3485).
- Deprecated `vars` parameter of `sample_posterior_predictive` in favor of `varnames`.
-  References to `live_plot` and corresponding notebooks have been removed.
- Deprecated `vars` parameters of `sample_posterior_predictive` and `sample_prior_predictive` in favor of `var_names`.  At least for the latter, this is more accurate, since the `vars` parameter actually took names.

### Contributors sorted by number of commits
    45  Luciano Paz
    38  Thomas Wiecki
    23  Colin Carroll
    19  Junpeng Lao
    15  Chris Fonnesbeck
    13  Juan Martín Loyola
    13  Ravin Kumar
     8  Robert P. Goldman
     5  Tim Blazina
     4  chang111
     4  adamboche
     3  Eric Ma
     3  Osvaldo Martin
     3  Sanmitra Ghosh
     3  Saurav Shekhar
     3  chartl
     3  fredcallaway
     3  Demetri
     2  Daisuke Kondo
     2  David Brochart
     2  George Ho
     2  Vaibhav Sinha
     1  rpgoldman
     1  Adel Tomilova
     1  Adriaan van der Graaf
     1  Bas Nijholt
     1  Benjamin Wild
     1  Brigitta Sipocz
     1  Daniel Emaasit
     1  Hari
     1  Jeroen
     1  Joseph Willard
     1  Juan Martin Loyola
     1  Katrin Leinweber
     1  Lisa Martin
     1  M. Domenzain
     1  Matt Pitkin
     1  Peadar Coyle
     1  Rupal Sharma
     1  Tom Gilliss
     1  changjiangeng
     1  michaelosthege
     1  monsta
     1  579397

## PyMC3 3.6 (Dec 21 2018)

This will be the last release to support Python 2.

### New features

- Track the model log-likelihood as a sampler stat for NUTS and HMC samplers
  (accessible as `trace.get_sampler_stats('model_logp')`) (#3134)
- Add Incomplete Beta function `incomplete_beta(a, b, value)`
- Add log CDF functions to continuous distributions: `Beta`, `Cauchy`, `ExGaussian`, `Exponential`, `Flat`, `Gumbel`, `HalfCauchy`, `HalfFlat`, `HalfNormal`, `Laplace`, `Logistic`, `LogNormal`, `Normal`, `Pareto`, `StudentT`, `Triangular`, `Uniform`, `Wald`, `Weibull`.
- Behavior of `sample_posterior_predictive` is now to produce posterior predictive samples, in order, from all values of the `trace`. Previously, by default it would produce 1 chain worth of samples, using a random selection from the `trace` (#3212)
- Show diagnostics for initial energy errors in HMC and NUTS.
- PR #3273 has added the `distributions.distribution._DrawValuesContext` context
  manager. This is used to store the values already drawn in nested `random`
  and `draw_values` calls, enabling `draw_values` to draw samples from the
  joint probability distribution of RVs and not the marginals. Custom
  distributions that must call `draw_values` several times in their `random`
  method, or that invoke many calls to other distribution's `random` methods
  (e.g. mixtures) must do all of these calls under the same `_DrawValuesContext`
  context manager instance. If they do not, the conditional relations between
  the distribution's parameters could be broken, and `random` could return
  values drawn from an incorrect distribution.
- `Rice` distribution is now defined with either the noncentrality parameter or the shape parameter (#3287).

### Maintenance

- Big rewrite of documentation (#3275)
- Fixed Triangular distribution `c` attribute handling in `random` and updated sample codes for consistency (#3225)
- Refactor SMC and properly compute marginal likelihood (#3124)
- Removed use of deprecated `ymin` keyword in matplotlib's `Axes.set_ylim` (#3279)
- Fix for #3210. Now `distribution.draw_values(params)`, will draw the `params` values from their joint probability distribution and not from combinations of their marginals (Refer to PR #3273).
- Removed dependence on pandas-datareader for retrieving Yahoo Finance data in examples (#3262)
- Rewrote `Multinomial._random` method to better handle shape broadcasting (#3271)
- Fixed `Rice` distribution, which inconsistently mixed two parametrizations (#3286).
- `Rice` distribution now accepts multiple parameters and observations and is usable with NUTS (#3289).
- `sample_posterior_predictive` no longer calls `draw_values` to initialize the shape of the ppc trace. This called could lead to `ValueError`'s when sampling the ppc from a model with `Flat` or `HalfFlat` prior distributions (Fix issue #3294).
- Added explicit conversion to `floatX` and `int32` for the continuous and discrete probability distribution parameters (addresses issue #3223).


### Deprecations

- Renamed `sample_ppc()` and `sample_ppc_w()` to `sample_posterior_predictive()` and `sample_posterior_predictive_w()`, respectively.

## PyMC3 3.5 (July 21 2018)

### New features

- Add documentation section on survival analysis and censored data models
- Add `check_test_point` method to `pm.Model`
- Add `Ordered` Transformation and `OrderedLogistic` distribution
- Add `Chain` transformation
- Improve error message `Mass matrix contains zeros on the diagonal. Some derivatives might always be zero` during tuning of `pm.sample`
- Improve error message `NaN occurred in optimization.` during ADVI
- Save and load traces without `pickle` using `pm.save_trace` and `pm.load_trace`
- Add `Kumaraswamy` distribution
- Add `TruncatedNormal` distribution
- Rewrite parallel sampling of multiple chains on py3. This resolves long standing issues when transferring large traces to the main process, avoids pickling issues on UNIX, and allows us to show a progress bar for all chains. If parallel sampling is interrupted, we now return partial results.
- Add `sample_prior_predictive` which allows for efficient sampling from the unconditioned model.
- SMC: remove experimental warning, allow sampling using `sample`, reduce autocorrelation from final trace.
- Add `model_to_graphviz` (which uses the optional dependency `graphviz`) to plot a directed graph of a PyMC3 model using plate notation.
- Add beta-ELBO variational inference as in beta-VAE model (Christopher P. Burgess et al. NIPS, 2017)
- Add `__dir__` to `SingleGroupApproximation` to improve autocompletion in interactive environments

### Fixes

- Fixed grammar in divergence warning, previously `There were 1 divergences ...` could be raised.
- Fixed `KeyError` raised when only subset of variables are specified to be recorded in the trace.
- Removed unused `repeat=None` arguments from all `random()` methods in distributions.
- Deprecated the `sigma` argument in `MarginalSparse.marginal_likelihood` in favor of `noise`
- Fixed unexpected behavior in `random`. Now the `random` functionality is more robust and will work better for `sample_prior` when that is implemented.
- Fixed `scale_cost_to_minibatch` behaviour, previously this was not working and always `False`

## PyMC3 3.4.1 (April 18 2018)

### New features

- Add `logit_p` keyword to `pm.Bernoulli`, so that users can specify the logit of the success probability. This is faster and more stable than using `p=tt.nnet.sigmoid(logit_p)`.
- Add `random` keyword to `pm.DensityDist` thus enabling users to pass custom random method which in turn makes sampling from a `DensityDist` possible.
- Effective sample size computation is updated. The estimation uses Geyer's initial positive sequence, which no longer truncates the autocorrelation series inaccurately. `pm.diagnostics.effective_n` now can reports N_eff>N.
- Added `KroneckerNormal` distribution and a corresponding `MarginalKron` Gaussian Process implementation for efficient inference, along with lower-level functions such as `cartesian` and `kronecker` products.
- Added `Coregion` covariance function.
- Add new 'pairplot' function, for plotting scatter or hexbin matrices of sampled parameters. Optionally it can plot divergences.
- Plots of discrete distributions in the docstrings
- Add logitnormal distribution
- Densityplot: add support for discrete variables
- Fix the Binomial likelihood in `.glm.families.Binomial`, with the flexibility of specifying the `n`.
- Add `offset` kwarg to `.glm`.
- Changed the `compare` function to accept a dictionary of model-trace pairs instead of two separate lists of models and traces.
- add test and support for creating multivariate mixture and mixture of mixtures
- `distribution.draw_values`, now is also able to draw values from conditionally dependent RVs, such as autotransformed RVs (Refer to PR #2902).

### Fixes

- `VonMises` does not overflow for large values of kappa. i0 and i1 have been removed and we now use log_i0 to compute the logp.
- The bandwidth for KDE plots is computed using a modified version of Scott's rule. The new version uses entropy instead of standard deviation. This works better for multimodal distributions. Functions using KDE plots has a new argument `bw` controlling the bandwidth.
- fix PyMC3 variable is not replaced if provided in more_replacements (#2890)
- Fix for issue #2900. For many situations, named node-inputs do not have a `random` method, while some intermediate node may have it. This meant that if the named node-input at the leaf of the graph did not have a fixed value, `theano` would try to compile it and fail to find inputs, raising a `theano.gof.fg.MissingInputError`. This was fixed by going through the theano variable's owner inputs graph, trying to get intermediate named-nodes values if the leafs had failed.
- In `distribution.draw_values`, some named nodes could be `theano.tensor.TensorConstant`s or `theano.tensor.sharedvar.SharedVariable`s. Nevertheless, in `distribution._draw_value`, these would be passed to `distribution._compile_theano_function` as if they were `theano.tensor.TensorVariable`s. This could lead to the following exceptions `TypeError: ('Constants not allowed in param list', ...)` or `TypeError: Cannot use a shared variable (...)`. The fix was to not add `theano.tensor.TensorConstant` or `theano.tensor.sharedvar.SharedVariable` named nodes into the `givens` dict that could be used in `distribution._compile_theano_function`.
- Exponential support changed to include zero values.

### Deprecations

- DIC and BPIC calculations have been removed
- df_summary have been removed, use summary instead
- `njobs` and `nchains` kwarg are deprecated in favor of `cores` and `chains` for `sample`
- `lag` kwarg in `pm.stats.autocorr` and `pm.stats.autocov` is deprecated.


## PyMC3 3.3 (January 9, 2018)

### New features

- Improve NUTS initialization `advi+adapt_diag_grad` and add `jitter+adapt_diag_grad` (#2643)
- Added `MatrixNormal` class for representing vectors of multivariate normal variables
- Implemented `HalfStudentT` distribution
- New benchmark suite added (see http://pandas.pydata.org/speed/pymc/)
- Generalized random seed types
- Update loo, new improved algorithm (#2730)
- New CSG (Constant Stochastic Gradient) approximate posterior sampling algorithm (#2544)
- Michael Osthege added support for population-samplers and implemented differential evolution metropolis (`DEMetropolis`).  For models with correlated dimensions that can not use gradient-based samplers, the `DEMetropolis` sampler can give higher effective sampling rates. (also see [PR#2735](https://github.com/pymc-devs/pymc/pull/2735))
- Forestplot supports multiple traces (#2736)
- Add new plot, densityplot (#2741)
- DIC and BPIC calculations have been deprecated
- Refactor HMC and implemented new warning system (#2677, #2808)

### Fixes

- Fixed `compareplot` to use `loo` output.
- Improved `posteriorplot` to scale fonts
- `sample_ppc_w` now broadcasts
- `df_summary` function renamed to `summary`
- Add test for `model.logp_array` and `model.bijection` (#2724)
- Fixed `sample_ppc` and `sample_ppc_w` to iterate all chains(#2633, #2748)
- Add Bayesian R2 score (for GLMs) `stats.r2_score` (#2696) and test (#2729).
- SMC works with transformed variables (#2755)
- Speedup OPVI (#2759)
- Multiple minor fixes and improvements in the docs (#2775, #2786, #2787, #2789, #2790, #2794, #2799, #2809)

### Deprecations

- Old (`minibatch-`)`advi` is removed (#2781)


## PyMC3 3.2 (October 10, 2017)

### New features

This version includes two major contributions from our Google Summer of Code 2017 students:

* Maxim Kochurov extended and refactored the variational inference module. This primarily adds two important classes, representing operator variational inference (`OPVI`) objects and `Approximation` objects. These make it easier to extend existing `variational` classes, and to derive inference from `variational` optimizations, respectively. The `variational` module now also includes normalizing flows (`NFVI`).
* Bill Engels added an extensive new Gaussian processes (`gp`) module. Standard GPs can be specified using either `Latent` or `Marginal` classes, depending on the nature of the underlying function. A Student-T process `TP` has been added. In order to accomodate larger datasets, approximate marginal Gaussian processes (`MarginalSparse`) have been added.

Documentation has been improved as the result of the project's monthly "docathons".

An experimental stochastic gradient Fisher scoring (`SGFS`) sampling step method has been added.

The API for `find_MAP` was enhanced.

SMC now estimates the marginal likelihood.

Added `Logistic` and `HalfFlat` distributions to set of continuous distributions.

Bayesian fraction of missing information (`bfmi`) function added to `stats`.

Enhancements to `compareplot` added.

QuadPotential adaptation has been implemented.

Script added to build and deploy documentation.

MAP estimates now available for transformed and non-transformed variables.

The `Constant` variable class has been deprecated, and will be removed in 3.3.

DIC and BPIC calculations have been sped up.

Arrays are now accepted as arguments for the `Bound` class.

`random` method was added to the `Wishart` and `LKJCorr` distributions.

Progress bars have been added to LOO and WAIC calculations.

All example notebooks updated to reflect changes in API since 3.1.

Parts of the test suite have been refactored.

### Fixes

Fixed sampler stats error in NUTS for non-RAM backends

Matplotlib is  no longer a hard dependency, making it easier to use in settings where installing Matplotlib is problematic. PyMC3 will only complain if plotting is attempted.

Several bugs in the Gaussian process covariance were fixed.

All chains are now used to calculate WAIC and LOO.

AR(1) log-likelihood function has been fixed.

Slice sampler fixed to sample from 1D conditionals.

Several docstring fixes.

### Contributors

The following people contributed to this release (ordered by number of commits):

Maxim Kochurov <maxim.v.kochurov@gmail.com>
Bill Engels <w.j.engels@gmail.com>
Chris Fonnesbeck <chris.fonnesbeck@vanderbilt.edu>
Junpeng Lao <junpeng.lao@unifr.ch>
Adrian Seyboldt <adrian.seyboldt@gmail.com>
AustinRochford <arochford@monetate.com>
Osvaldo Martin <aloctavodia@gmail.com>
Colin Carroll <colcarroll@gmail.com>
Hannes Vasyura-Bathke <hannes.bathke@gmx.net>
Thomas Wiecki <thomas.wiecki@gmail.com>
michaelosthege <thecakedev@hotmail.com>
Marco De Nadai <me@marcodena.it>
Kyle Beauchamp <kyleabeauchamp@gmail.com>
Massimo <mcavallaro@users.noreply.github.com>
ctm22396 <ctm22396@gmail.com>
Max Horn <maexlich@gmail.com>
Hennadii Madan <madanh2014@gmail.com>
Hassan Naseri <h.nasseri@gmail.com>
Peadar Coyle <peadarcoyle@googlemail.com>
Saurav R. Tuladhar <saurav@fastmail.com>
Shashank Shekhar <shashank.f1@gmail.com>
Eric Ma <ericmjl@users.noreply.github.com>
Ed Herbst <ed.herbst@gmail.com>
tsdlovell <dlovell@twosigma.com>
zaxtax <zaxtax@users.noreply.github.com>
Dan Nichol <daniel.nichol@univ.ox.ac.uk>
Benjamin Yetton <bdyetton@gmail.com>
jackhansom <jack.hansom@outlook.com>
Jack Tsai <jacksctsai@gmail.com>
Andrés Asensio Ramos <aasensioramos@gmail.com>


## PyMC3 3.1 (June 23, 2017)

### New features

* New user forum at http://discourse.pymc.io

* [Gaussian Process submodule](http://pymc-devs.github.io/pymc/notebooks/GP-introduction.html)

* Much improved variational inference support:

  - [Add Operator Variational Inference (experimental).](http://pymc-devs.github.io/pymc/notebooks/bayesian_neural_network_opvi-advi.html)

  - [Add Stein-Variational Gradient Descent as well as Amortized SVGD (experimental).](https://github.com/pymc-devs/pymc/pull/2183)

  - [Add pm.Minibatch() to easily specify mini-batches.](http://pymc-devs.github.io/pymc/notebooks/bayesian_neural_network_opvi-advi.html#Minibatch-ADVI)

  - Added various optimizers including ADAM.

  - Stopping criterion implemented via callbacks.

* sample() defaults changed: tuning is enabled for the first 500 samples which are then discarded from the trace as burn-in.

* MvNormal supports Cholesky Decomposition now for increased speed and numerical stability.

* Many optimizations and speed-ups.

* NUTS implementation now matches current Stan implementation.

* Add higher-order integrators for HMC.

* [Add sampler statistics.](http://pymc-devs.github.io/pymc/notebooks/sampler-stats.html)

* [Add live-trace to see samples in real-time.](http://pymc-devs.github.io/pymc/notebooks/live_sample_plots.html)

* ADVI stopping criterion implemented.

* Improved support for theano's floatX setting to enable GPU computations (work in progress).

* MvNormal supports Cholesky Decomposition now for increased speed and numerical stability.

* [Add Elliptical Slice Sampler.](http://pymc-devs.github.io/pymc/notebooks/GP-slice-sampling.html)

* Added support for multidimensional minibatches

* [Sampled posteriors can now be turned into priors for Bayesian updating with a new interpolated distribution.](https://github.com/pymc-devs/pymc/pull/2163)

* Added `Approximation` class and the ability to convert a sampled trace into an approximation via its `Empirical` subclass.

* `Model` can now be inherited from and act as a base class for user specified models (see pymc.models.linear).

* Add MvGaussianRandomWalk and MvStudentTRandomWalk distributions.

* GLM models do not need a left-hand variable anymore.

* Refactored HMC and NUTS for better readability.

* Add support for Python 3.6.

### Fixes

* Bound now works for discrete distributions as well.

* Random sampling now returns the correct shape even for higher dimensional RVs.

* Use theano Psi and GammaLn functions to enable GPU support for them.


## PyMC3 3.0 (January 9, 2017)

We are proud and excited to release the first stable version of PyMC3, the product of more than [5 years](https://github.com/pymc-devs/pymc/commit/85c7e06b6771c0d99cbc09cb68885cda8f7785cb) of ongoing development and contributions from over 80 individuals. PyMC3 is a Python module for Bayesian modeling which focuses on modern Bayesian computational methods, primarily gradient-based (Hamiltonian) MCMC sampling and variational inference. Models are specified in Python, which allows for great flexibility. The main technological difference in PyMC3 relative to previous versions is the reliance on Theano for the computational backend, rather than on Fortran extensions.

### New features

Since the beta release last year, the following improvements have been implemented:

* Added `variational` submodule, which features the automatic differentiation variational inference (ADVI) fitting method. Also supports mini-batch ADVI for large data sets. Much of this work was due to the efforts of Taku Yoshioka, and important guidance was provided by the Stan team (specifically Alp Kucukelbir and Daniel Lee).

* Added model checking utility functions, including leave-one-out (LOO) cross-validation, BPIC, WAIC, and DIC.

* Implemented posterior predictive sampling (`sample_ppc`).

* Implemented auto-assignment of step methods by `sample` function.

* Enhanced IPython Notebook examples, featuring more complete narratives accompanying code.

* Extensive debugging of NUTS sampler.

* Updated documentation to reflect changes in code since beta.

* Refactored test suite for better efficiency.

* Added von Mises, zero-inflated negative binomial, and Lewandowski, Kurowicka and Joe (LKJ) distributions.

* Adopted `joblib` for managing parallel computation of chains.

* Added contributor guidelines, contributor code of conduct and governance document.

### Deprecations

* Argument order of tau and sd was switched for distributions of the normal family:
- `Normal()`
- `LogNormal()`
- `HalfNormal()`

Old: `Normal(name, mu, tau)`
New: `Normal(name, mu, sd)` (supplying keyword arguments is unaffected).

* `MvNormal` calling signature changed:
Old: `MvNormal(name, mu, tau)`
New: `MvNormal(name, mu, cov)` (supplying keyword arguments is unaffected).

We on the PyMC3 core team would like to thank everyone for contributing and now feel that this is ready for the big time. We look forward to hearing about all the cool stuff you use PyMC3 for, and look forward to continued development on the package.

### Contributors

The following authors contributed to this release:

Chris Fonnesbeck <chris.fonnesbeck@vanderbilt.edu>
John Salvatier <jsalvatier@gmail.com>
Thomas Wiecki <thomas.wiecki@gmail.com>
Colin Carroll <colcarroll@gmail.com>
Maxim Kochurov <maxim.v.kochurov@gmail.com>
Taku Yoshioka <taku.yoshioka.4096@gmail.com>
Peadar Coyle (springcoil) <peadarcoyle@googlemail.com>
Austin Rochford <arochford@monetate.com>
Osvaldo Martin <aloctavodia@gmail.com>
Shashank Shekhar <shashank.f1@gmail.com>

In addition, the following community members contributed to this release:

A Kuz <for.akuz@gmail.com>
A. Flaxman <abie@alum.mit.edu>
Abraham Flaxman <abie@alum.mit.edu>
Alexey Goldin <alexey.goldin@gmail.com>
Anand Patil <anand.prabhakar.patil@gmail.com>
Andrea Zonca <code@andreazonca.com>
Andreas Klostermann <andreasklostermann@googlemail.com>
Andres Asensio Ramos
Andrew Clegg <andrew.clegg@pearson.com>
Anjum48
Benjamin Edwards <bedwards@cs.unm.edu>
Boris Avdeev <borisaqua@gmail.com>
Brian Naughton <briannaughton@gmail.com>
Byron Smith
Chad Heyne <chadheyne@gmail.com>
Corey Farwell <coreyf@rwell.org>
David Huard <david.huard@gmail.com>
David Stück <dstuck@users.noreply.github.com>
DeliciousHair <mshepit@gmail.com>
Dustin Tran
Eigenblutwurst <Hannes.Bathke@gmx.net>
Gideon Wulfsohn <gideon.wulfsohn@gmail.com>
Gil Raphaelli <g@raphaelli.com>
Gogs <gogitservice@gmail.com>
Ilan Man
Imri Sofer <imrisofer@gmail.com>
Jake Biesinger <jake.biesinger@gmail.com>
James Webber <jamestwebber@gmail.com>
John McDonnell <john.v.mcdonnell@gmail.com>
Jon Sedar <jon.sedar@applied.ai>
Jordi Diaz
Jordi Warmenhoven <jordi.warmenhoven@gmail.com>
Karlson Pfannschmidt <kiudee@mail.uni-paderborn.de>
Kyle Bishop <citizenphnix@gmail.com>
Kyle Meyer <kyle@kyleam.com>
Lin Xiao
Mack Sweeney <mackenzie.sweeney@gmail.com>
Matthew Emmett <memmett@unc.edu>
Michael Gallaspy <gallaspy.michael@gmail.com>
Nick <nalourie@example.com>
Osvaldo Martin <aloctavodia@gmail.com>
Patricio Benavente <patbenavente@gmail.com>
Raymond Roberts
Rodrigo Benenson <rodrigo.benenson@gmail.com>
Sergei Lebedev <superbobry@gmail.com>
Skipper Seabold <chris.fonnesbeck@vanderbilt.edu>
Thomas Kluyver <takowl@gmail.com>
Tobias Knuth <mail@tobiasknuth.de>
Volodymyr Kazantsev
Wes McKinney <wesmckinn@gmail.com>
Zach Ploskey <zploskey@gmail.com>
akuz <for.akuz@gmail.com>
brandon willard <brandonwillard@gmail.com>
dstuck <dstuck88@gmail.com>
ingmarschuster <ingmar.schuster.linguistics@gmail.com>
jan-matthis <mail@jan-matthis.de>
jason <JasonTam22@gmailcom>
kiudee <quietdeath@gmail.com>
maahnman <github@mm.maahn.de>
macgyver <neil.rabinowitz@merton.ox.ac.uk>
mwibrow <mwibrow@gmail.com>
olafSmits <o.smits@gmail.com>
paul sorenson <paul@metrak.com>
redst4r <redst4r@web.de>
santon <steven.anton@idanalytics.com>
sgenoud <stevegenoud+github@gmail.com>
stonebig <stonebig>
Tal Yarkoni <tyarkoni@gmail.com>
x2apps <x2apps@yahoo.com>
zenourn <daniel@zeno.co.nz>

## PyMC3 3.0b (June 16th, 2015)

Probabilistic programming allows for flexible specification of Bayesian statistical models in code. PyMC3 is a new, open-source probabilistic programmer framework with an intuitive, readable and concise, yet powerful, syntax that is close to the natural notation statisticians use to describe models. It features next-generation fitting techniques, such as the No U-Turn Sampler, that allow fitting complex models with thousands of parameters without specialized knowledge of fitting algorithms.

PyMC3 has recently seen rapid development. With the addition of two new major features: automatic transforms and missing value imputation, PyMC3 has become ready for wider use. PyMC3 is now refined enough that adding features is easy, so we don't expect adding features in the future will require drastic changes. It has also become user friendly enough for a broader audience. Automatic transformations mean NUTS and find_MAP work with less effort, and friendly error messages mean its easy to diagnose problems with your model.

Thus, Thomas, Chris and I are pleased to announce that PyMC3 is now in Beta.

### Highlights
* Transforms now automatically applied to constrained distributions
* Transforms now specified with a `transform=` argument on Distributions. `model.TransformedVar` is gone.
* Transparent missing value imputation support added with MaskedArrays or pandas.DataFrame NaNs.
* Bad default values now ignored
* Profile theano functions using `model.profile(model.logpt)`

### Contributors since 3.0a
* A. Flaxman <abie@alum.mit.edu>
* Andrea Zonca <code@andreazonca.com>
* Andreas Klostermann <andreasklostermann@googlemail.com>
* Andrew Clegg <andrew.clegg@pearson.com>
* AustinRochford <arochford@monetate.com>
* Benjamin Edwards <bedwards@cs.unm.edu>
* Brian Naughton <briannaughton@gmail.com>
* Chad Heyne <chadheyne@gmail.com>
* Chris Fonnesbeck <chris.fonnesbeck@vanderbilt.edu>
* Corey Farwell <coreyf@rwell.org>
* John Salvatier <jsalvatier@gmail.com>
* Karlson Pfannschmidt <quietdeath@gmail.com>
* Kyle Bishop <citizenphnix@gmail.com>
* Kyle Meyer <kyle@kyleam.com>
* Mack Sweeney <mackenzie.sweeney@gmail.com>
* Osvaldo Martin <aloctavodia@gmail.com>
* Raymond Roberts <rayvroberts@gmail.com>
* Rodrigo Benenson <rodrigo.benenson@gmail.com>
* Thomas Wiecki <thomas.wiecki@gmail.com>
* Zach Ploskey <zploskey@gmail.com>
* maahnman <github@mm.maahn.de>
* paul sorenson <paul@metrak.com>
* zenourn <daniel@zeno.co.nz>
# PyMC Code of Conduct

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
reported by contacting PyMC developer Christopher Fonnesbeck via email
(chris.fonnesbeck@vanderbilt.edu) or phone (615-955-0380). Alternatively, you
may also contact NumFOCUS Executive Director Leah Silen (512-222-5449), as PyMC
is a member of NumFOCUS and subscribes to their code of conduct as a
precondition for continued membership. All complaints will be reviewed and
investigated and will result in a response that is deemed necessary and
appropriate to the circumstances. The project team is obligated to maintain
confidentiality with regard to the reporter of an incident. Further details of
specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org
# Guidelines for Contributing

As a scientific community-driven software project, PyMC welcomes contributions from interested individuals or groups. These guidelines are provided to give potential contributors information to make their contribution compliant with the conventions of the PyMC project, and maximize the probability of such contributions to be merged as quickly and efficiently as possible.

There are 4 main ways of contributing to the PyMC project (in descending order of difficulty or scope):

* Adding new or improved functionality to the existing codebase
* Fixing outstanding issues (bugs) with the existing codebase. They range from low-level software bugs to higher-level design problems.
* Contributing or improving the documentation (`docs`) or examples (`pymc/examples`)
* Submitting issues related to bugs or desired enhancements

# Opening issues

We appreciate being notified of problems with the existing PyMC code. We prefer that issues be filed the on [Github Issue Tracker](https://github.com/pymc-devs/pymc/issues), rather than on social media or by direct email to the developers.

Please verify that your issue is not being currently addressed by other issues or pull requests by using the GitHub search tool to look for key words in the project issue tracker.

Filter on the ["beginner friendly"](https://github.com/pymc-devs/pymc/issues?q=is%3Aopen+is%3Aissue+label%3A%22beginner+friendly%22) label for issues which are good for new contributors.

# Contributing code via pull requests

While issue reporting is valuable, we strongly encourage users who are inclined to do so to submit patches for new or existing issues via pull requests. This is particularly the case for simple fixes, such as typos or tweaks to documentation, which do not require a heavy investment of time and attention.

Contributors are also encouraged to contribute new code to enhance PyMC's functionality, also via pull requests. Please consult the [PyMC documentation](https://pymc-devs.github.io/pymc/) to ensure that any new contribution does not strongly overlap with existing functionality.

The preferred workflow for contributing to PyMC is to fork the [GitHub repository](https://github.com/pymc-devs/pymc/), clone it to your local machine, and develop on a feature branch.

## Steps:

1. Fork the [project repository](https://github.com/pymc-devs/pymc/) by clicking on the 'Fork' button near the top right of the main repository page. This creates a copy of the code under your GitHub user account.

2. Clone your fork of the PyMC repo from your GitHub account to your local disk, and add the base repository as a remote:

   ```bash
   $ git clone git@github.com:<your GitHub handle>/pymc.git
   $ cd pymc
   $ git remote add upstream git@github.com:pymc-devs/pymc.git
   ```

3. Create a ``feature`` branch to hold your development changes:

   ```bash
   $ git checkout -b my-feature
   ```

   Always use a ``feature`` branch. It's good practice to never routinely work on the ``main`` branch of any repository.

4. Project requirements are in ``requirements.txt``, and libraries used for development are in ``requirements-dev.txt``. The easiest (and recommended) way to set up a development environment is via [miniconda](https://docs.conda.io/en/latest/miniconda.html):

   ```bash
   $ conda env create -f conda-envs/environment-dev-py37.yml  # or py38 or py39
   $ conda activate pymc-dev-py37
   $ pip install -e .
   ```

   _Alternatively_ you may (probably in a [virtual environment](https://docs.python-guide.org/dev/virtualenvs/)) run:

   ```bash
   $ pip install -e .
   $ pip install -r requirements-dev.txt
   ```

   Yet another alternative is to create a docker environment for development. See: [Developing in Docker](#Developing-in-Docker).

5. Develop the feature on your feature branch. Add changed files using ``git add`` and then ``git commit`` files:

   ```bash
   $ git add modified_files
   $ git commit
   ```

   to record your changes locally.
   After committing, it is a good idea to sync with the base repository in case there have been any changes:
   ```bash
   $ git fetch upstream
   $ git rebase upstream/main
   ```

   Then push the changes to your GitHub account with:

   ```bash
   $ git push -u origin my-feature
   ```

6. Go to the GitHub web page of your fork of the PyMC repo. Click the 'Pull request' button to send your changes to the project's maintainers for review. This will send an email to the committers.

## Pull request checklist

We recommended that your contribution complies with the following guidelines before you submit a pull request:

*  If your pull request addresses an issue, please use the pull request title to describe the issue and mention the issue number in the pull request description. This will make sure a link back to the original issue is created.

*  All public methods must have informative docstrings with sample usage when appropriate.

*  Please prefix the title of incomplete contributions with `[WIP]` (to indicate a work in progress). WIPs may be useful to (1) indicate you are working on something to avoid duplicated work, (2) request broad review of functionality or API, or (3) seek collaborators.

*  All other tests pass when everything is rebuilt from scratch.  See
[Developing in Docker](#Developing-in-Docker) for information on running the test suite locally.

*  When adding additional functionality, provide at least one example script or Jupyter Notebook in the ``pymc/examples/`` folder. Have a look at other examples for reference. Examples should demonstrate why the new functionality is useful in practice and, if possible, compare it to other methods available in PyMC.

* Documentation and high-coverage tests are necessary for enhancements to be accepted.

* Run any of the pre-existing examples in ``docs/source/notebooks`` that contain analyses that would be affected by your changes to ensure that nothing breaks. This is a useful opportunity to not only check your work for bugs that might not be revealed by unit test, but also to show how your contribution improves PyMC for end users.

You can also check for common programming errors with the following
tools:

* Code with good test **coverage** (at least 80%), check with:

  ```bash
  $ pip install pytest pytest-cov coverage
  $ pytest --cov=pymc pymc/tests/<name of test>.py
  ```

* No `pre-commit` errors: see the [Python code style](https://github.com/pymc-devs/pymc/wiki/Python-Code-Style) and [Jupyter Notebook style](https://github.com/pymc-devs/pymc/wiki/PyMC-Jupyter-Notebook-Style-Guide) page from our Wiki on how to install and run it.

## Developing in Docker

We have provided a Dockerfile which helps for isolating build problems, and local development.
Install [Docker](https://www.docker.com/) for your operating system, clone this repo, then
run `./scripts/start_container.sh`. This should start a local docker container called `pymc`,
as well as a [`jupyter`](http://jupyter.org/) notebook server running on port 8888. The
notebook should be opened in your browser automatically (you can disable this by passing
`--no-browser`). The repo will be running the code from your local copy of `pymc`,
so it is good for development.

You may also use it to run the test suite, with

```bash
$  docker exec -it pymc  bash # logon to the container
$  cd ~/pymc/tests
$  . ./../../scripts/test.sh # takes a while!
```

This should be quite close to how the tests run on TravisCI.

If the container was started without opening the browser, you
need the notebook instances token to work with the notebook. This token can be
accessed with

```
docker exec -it pymc jupyter notebook list
```

## Style guide

We have configured a pre-commit hook that checks for `black`-compliant code style.
We encourage you to configure the pre-commit hook as described in the [PyMC Python Code Style Wiki Page](https://docs.pymc.io/en/latest/contributing/python_style.html), because it will automatically enforce the code style on your commits.

Similarly, consult the [PyMC's Jupyter Notebook Style](https://docs.pymc.io/en/latest/contributing/jupyter_style.html) guides for notebooks.

For documentation strings, we *prefer* [numpy style](https://numpydoc.readthedocs.io/en/latest/format.html) to comply with the style that predominates in our upstream dependencies.

#### This guide was derived from the [scikit-learn guide to contributing](https://github.com/scikit-learn/scikit-learn/blob/master/CONTRIBUTING.md)

**Thank your for opening a PR!**

Depending on what your PR does, here are a few things you might want to address in the description:
+ [ ] what are the (breaking) changes that this PR makes?
+ [ ] important background, or details about the implementation
+ [ ] are the changes—especially new features—covered by tests and docstrings?
+ [ ] [linting/style checks have been run](https://github.com/pymc-devs/pymc3/wiki/PyMC3-Python-Code-Style)
+ [ ] [consider adding/updating relevant example notebooks](https://github.com/pymc-devs/pymc-examples)
+ [ ] right before it's ready to merge, mention the PR in the RELEASE-NOTES.md
---
name: 'Bug Report'
about: Inform about bugs in the PyMC/PyMC3 library

---

## Description of your problem

**Please provide a minimal, self-contained, and reproducible example.**
```python
[Your code here]
```

**Please provide the full traceback.**

<details><summary>Complete error traceback</summary>

```python
[The complete error output here]
```

</details>

**Please provide any additional information below.**


## Versions and main components

* PyMC/PyMC3 Version:
* Aesara/Theano Version:
* Python Version:
* Operating system:
* How did you install PyMC/PyMC3: (conda/pip)
---
sd_hide_title: true
---
(index)=
# PyMC Documentation

:::{card}
:img-top: https://raw.githubusercontent.com/pymc-devs/pymc/main/docs/logos/PyMC.jpg
:margin: 0 2 auto auto
:width: 50%
:text-align: center
:shadow: none

+++
Probabilistic Programming in Python
:::

::::{div} sd-d-flex-row sd-align-major-center
:::{button-ref} learning
:color: primary
:ref-type: ref
:class: sd-fs-2 sd-px-5

Get started!
:::
::::


:::{card} Announcements: library name change and launching PyMC 4.0!
:width: 75%
:margin: auto

We have two major announcements that we're excited to share. First of all, a new name for our library: the PyMC3 library you know and love is now called PyMC. PyMC3 still exists, as a specific major release between PyMC2 and PyMC 4.0. Read more about the renaming and how to solve related issues you might experience from this update [here]().

This ties into our second announcement, which is that we are hereby launching the newest version of PyMC: PyMC 4.0! Read more about this new release [here]().
:::

---

# Main Features & Benefits

::::::{grid} 1 1 2 2
:gutter: 1

:::::{grid-item}

::::{grid} 1
:gutter: 1

:::{grid-item}

**Friendly modelling API**

PyMC allows you to write down models using an intuitive syntax to describe a data generating process.
:::
:::{grid-item}

**Cutting edge algorithms and model building blocks**

Fit your model using gradient-based MCMC algorithms like NUTS, using ADVI for fast approximate inference — including minibatch-ADVI for scaling to large datasets — or using Gaussian processes to build Bayesian nonparametric models.
:::
::::
:::::

:::::{grid-item}

```{code-block} python
---
---
    import numpy as np
    import pymc as pm

    X = np.random.normal(size=100)
    y = np.random.normal(X) * 1.2

    with pm.Model() as linear_model:
        weights = pm.Normal("weights", mu=0, sigma=1)
        noise = pm.Gamma("noise", alpha=2, beta=1)
        y_observed = pm.Normal(
            "y_observed",
            mu=X @ weights,
            sigma=noise,
            observed=y,
        )

        prior = pm.sample_prior_predictive()
        posterior = pm.sample()
        posterior_pred = pm.sample_posterior_predictive(posterior)
```
:::::

::::::

# Support

PyMC is a non-profit project under NumFOCUS umbrella. If you value PyMC and want to support its development, consider donating to the project.

::::{div} sd-d-flex-row sd-align-major-center
:::{button-link} https://numfocus.org/donate-to-pymc3
:color: secondary
:class: sd-fs-2 sd-px-5


Donate
:::
::::

## Our sponsors

::::{grid} 2 4 4 4

:::{grid-item}
:::

:::{grid-item-card}
:img-background: _static/sponsors/numfocus.png
:link: https://numfocus.org/
:shadow: none
:::

:::{grid-item-card}
:img-background: _static/sponsors/pymc-labs.png
:link: https://www.pymc-labs.io/
:shadow: none
:::

:::{grid-item}
:::

::::

# Testimonials

:::::{card-carousel} 2

::::{include} about/featured_testimonials.md
::::

:::::

Find more testimonials {ref}`here <testimonials>`

# Citing PyMC

Use this to cite the library:

Salvatier J., Wiecki T.V., Fonnesbeck C. (2016) Probabilistic programming in Python using PyMC3. PeerJ Computer Science 2:e55 [DOI: 10.7717/peerj-cs.55.](https://doi.org/10.7717/peerj-cs.55)

For detailed advise on citing PyMC go to {ref}`citing_pymc`.
See [Google Scholar](https://scholar.google.de/scholar?oi=bibs&hl=en&authuser=1&cites=6936955228135731011) for a continuously updated list of papers citing PyMC3.

:::{toctree}
:maxdepth: 1
:hidden:

installation
learning
Examples <https://docs.pymc.io/projects/examples/en/latest/>
api
community
contributing/index
about/index
:::
(installation)=
# Installation

To install PyMC, select the operating system where you want to perform the installation.

## Linux
[Linux installation guide](https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(Linux))

## MacOS
[MacOS installation guide](https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(MacOS))

## Windows
[Windows installation guide](https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(Windows))

# Migrate to PyMC version 4

TODO: paste migration guide here
(learning)=
# Learning

## Getting started

Start here to get acquainted with the core concepts of Bayesian analysis and PyMC. The following resources only assume a very basic knowledge of code and statistics.

### {octicon}`book;1em;sd-text-info` Introductory books

#### Bayesian Methods for Hackers

By Cameron Davidson-Pilon

The "hacker" in the title  means learn-as-you-code. This hands-on introduction teaches intuitive definitions of the Bayesian approach to statistics, worklflow and decision-making by applying them using PyMC.

[Github repo](https://github.com/CamDavidsonPilon/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers)

[Project homepage](http://camdavidsonpilon.github.io/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers/)

#### Bayesian Analysis with Python

By Osvaldo Martin

A great introductory book written by a maintainer of PyMC. It provides a hands-on introduction to the main concepts of Bayesian statistics using synthetic and real data sets. Mastering the concepts in this book is a great foundation to pursue more advanced knowledge.

[Book website](https://www.packtpub.com/big-data-and-business-intelligence/bayesian-analysis-python-second-edition)

[Code and errata in PyMC](https://github.com/aloctavodia/BAP)

### {octicon}`mortar-board;1em;sd-text-info` Tutorial notebooks

#### Getting started
The {ref}`pymc_overview` notebook in our documentation shows the PyMC 4.0 code in action

#### General Linear Models: Linear regression

The {ref}`GLM_linear` notebook provides a gentle introduction to Bayesian linear regression and how it differs from the frequentist approach, and showcases how to implement it using PyMC.

#### Comparing models

The {ref}`model_comparison` notebook demonstrates the use of model comparison criteria in PyMC.

#### Validating a model using prior and posterior predictive checks

The {ref}`posterior_predictive` notebooks explains what prior and posterior predictive checks are and how to implement them in PyMC to validate your model.

### {octicon}`list-unordered;1em;sd-text-info` Glossary

PyMC's own {doc}`glossary` defines many core terms and provides useful references.

---

## Using PyMC

(the path in pymc-examples repo excluding the `examples` folder, or using manual anchors if we decide to create them).

::::{grid} 2
:gutter: 4

:::{grid-item-card} TODO


TODO
:::

::::

---
## Diving deeper

Links to intermediate notebooks. {doc}`More about step 3... <learn/step3>`


:::{toctree}
:hidden:

learn/examples
learn/books
learn/videos_and_podcasts
glossary
:::
(glossary)=
# Glossary

A glossary of common terms used throughout the PyMC documentation and examples.

:::::{glossary}
:sorted:

Functional Programming
  Functional programming is a programming style that prefers the use of basic functions with explicit and distinct inputs and outputs.
  This contrasts with functions or methods that depend on variables that are not explicitly passed as an input (such as accessing `self.variable` inside a method) or that alter the inputs or other state variables in-place, instead of returning new distinct variables as outputs.

Dispatching
  Choosing which function or method implementation to use based on the type of the input variables (usually just the first variable). For some examples, see Python's documentation for the [singledispatch](https://docs.python.org/3/library/functools.html#functools.singledispatch) decorator.

[Dispersion](https://en.wikipedia.org/wiki/Statistical_dispersion)
  In statistics, dispersion (also called variability, scatter, or spread) is the extent to which a distribution is stretched or squeezed

[Overdispersion](https://en.wikipedia.org/wiki/Overdispersion)
  In statistics, overdispersion is the presence of greater {term}`variability <dispersion>` in a data set than would be expected based on a given statistical model.

Underdispersion
  In statistics, underdispersion is the presence of lower {term}`variability <dispersion>` in a data set than would be expected based on a given statistical model.

Bayesian Workflow
  The Bayesian workflow involves all the steps needed for model building. This includes {term}`Bayesian inference` but also other tasks such as i) diagnoses of the quality of the inference, ii) model criticism, including evaluations of both model assumptions and model predictions, iii) comparison of models, not
just for the purpose of model selection or model averaging but more importantly to better understand these models and iv) Preparation of the results for a particular audience. These non-inferencial tasks require both numerical and visual summaries to help practitioners analyse their models. And they are sometimes collectively known as [Exploratory Analysis of Bayesian Models](https://joss.theoj.org/papers/10.21105/joss.01143).
  - For a compact overview, see Bayesian statistics and modelling by van de Schoot, R., Depaoli, S., King, R. et al in Nat Rev Methods - Primers 1, 1 (2021).
  - For an in-depth overview, see Bayesian Workflow by Andrew Gelman, Aki Vehtari, Daniel Simpson, Charles C. Margossian, Bob Carpenter, Yuling Yao, Lauren Kennedy, Jonah Gabry, Paul-Christian Bürkner, Martin Modrák
  - For an exercise-based material, see Think Bayes 2e: Bayesian Statistics Made Simple by Allen B. Downey
  - For an upcoming textbook that uses PyMC, Tensorflow Probability, and ArviZ libraries, see Bayesian Modeling and Computation by Osvaldo A. Martin, Ravin Kumar and Junpeng Lao

Bayesian inference
  Once we have defined the statistical model, Bayesian inference processes the data and model to produce a {term}`posterior` distribution. That is a joint distribution of all parameters in the model. This distribution is used to represent plausibility, and is the logical consequence of the model and data.

Bayesian model
  A Bayesian model is a composite of variables and distributional definitions for these variables. Bayesian models have two defining characteristics: i) Unknown quantities are described using probability distributions and ii) Bayes' theorem is used to update the values of the parameters conditioned on the data

Prior
  Bayesian statistics allow us, in principle, to include all information we have about the structure of the problem into a model. We can do this via assuming prior distributions of the model’s parameters. Priors represent the plausibility of the value of the parameters before accounting for the data. Priors multiplied by {term}`likelihood` produce the {term}`posterior`.

  Priors’ informativeness can fall anywhere on the complete uncertainty to relative certainty continuum. An informative prior might encode known restrictions on the possible range of values of that parameter.

  To understand the implications of a prior and likelihood we can simulate predictions from the model, before seeing any data. This can be done by taking samples from the prior predictive distribution.

  - For an in-depth guide to priors, consider Statistical Rethinking 2nd Edition By Richard McElreath, especially chapters 2.3

Likelihood
  There are many perspectives on likelihood, but conceptually we can think about it as the probability of the data, given the parameters. Or in other words, as the relative number of ways the data could have been produced.

  - For an in-depth unfolding of the concept, refer to Statistical Rethinking 2nd Edition By Richard McElreath, particularly chapter 2.
  - For the problem-based material, see Think Bayes 2e: Bayesian Statistics Made Simple by Allen B. Downey
  - For univariate, continuous scenarios, see the calibr8 paper: Bayesian calibration, process modeling and uncertainty quantification in biotechnology by Laura Marie Helleckes,  Michael Osthege, Wolfgang Wiechert, Eric von Lieres, Marco Oldiges

Posterior
  The outcome of Bayesian inference is a posterior distribution, which describes the relative plausibilities of every possible combination of parameter values, given the observed data. We can think of the posterior as the updated {term}`priors` after the model has seen the data.

  When the posterior is obtained using numerical methods we generally need to first diagnose the quality of the computed approximation. This is necessary as, for example, methods like {term}`MCMC` has only asymptotic guarantees. In a Bayesian setting predictions can be simulated by sampling from the posterior predictive distribution. When such predictions are used to check the internal consistency of the models by comparing it with the observed data used for inference, the process is known as the posterior predictive checks.

  Once you are satisfied with the model, posterior distribution can be summarized and interpreted. Common questions for the posterior include: intervals of defined boundaries, intervals of defined probability mass, and point estimates. When the posterior is very similar to the prior, the available data does not contain much information about a parameter of interest.

  - For more on generating and interpreting the posterior samples, see Statistical Rethinking 2nd Edition By Richard McElreath, chapter 3.

Generalized Linear Model
GLM
  In a Generalized Linear Model (GLM), we assume the response variable $y_i$ to follow an
  exponential family distribution with mean $\mu_i$, which is assumed to be some (often nonlinear)
  function of $x_i^T\beta$. They're considered linear because the covariates affect the distribution
  of $Y_i$ only through the linear combination $x_i^T\beta$. Some examples of Generalized Linear
  Models are: Linear Regression, ANOVA, Logistic Regression and Poisson Regression

  :::{note} Do not confuse these with general linear models
  :::

[Probability Mass Function](https://en.wikipedia.org/wiki/Probability_mass_function)
[PMF](https://en.wikipedia.org/wiki/Probability_mass_function)
  A function that gives the probability that a discrete random variable is exactly equal to some value.

[Maximum a Posteriori](https://en.wikipedia.org/wiki/Maximum_a_posteriori_estimation)
[MAP](https://en.wikipedia.org/wiki/Maximum_a_posteriori_estimation)
  It is a point-estimate of an unknown quantity, that equals the mode of the posterior distribution.

  If the prior distribution is a flat distribution, the MAP method is numerically equivalent to the {term}`Maximum Likelihood Estimate` (MLE).
  When the prior is not flat the MAP estimation can be seen as a regularized version of the MLE.

  - For a concise comparison between {term}`MLE` and {term}`MAP`, consider Deep Learning by Ian Goodfellow, chapter 5.6.1 or [Machine Learning: a Probabilistic Perspective](https://probml.github.io/pml-book/book1.html) by Kevin Murphy.

[No-U-Turn Sampler](https://arxiv.org/abs/1111.4246)
[NUTS](https://arxiv.org/abs/1111.4246)
  An extension of {term}`Hamiltonian Monte Carlo` that algorithmically sets likely candidate points that spans a wide swath of the target distribution, stopping automatically when it starts to double back and retrace its steps.

[Hamiltonian Monte Carlo](https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo)
[HMC](https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo)
  A {term}`Markov Chain Monte Carlo` method for obtaining a sequence of random samples which converge to being distributed according to a target probability distribution.

[Credibility](https://en.wikipedia.org/wiki/Credibility_theory)
  A form of statistical inference used to forecast an uncertain future event

[Ordinary Differential Equation](https://en.wikipedia.org/wiki/Ordinary_differential_equation)
[ODE](https://en.wikipedia.org/wiki/Ordinary_differential_equation)
  A type of differential equation containing one or more functions of one independent variable and the derivatives of those functions

Hierarchical Ordinary Differential Equation
  Individual, group, or other level types calculations of {term}`Ordinary Differential Equation`'s.

[Generalized Poisson Distribution](https://doi.org/10.2307/1267389)
  A generalization of the {term}`Poisson distribution`, with two parameters X1, and X2, is obtained as a limiting form of the generalized negative binomial distribution. The variance of the distribution is greater than, equal to or smaller than the mean according as X2 is positive, zero or negative. For formula and more detail, visit the link in the title.

[Bayes' theorem](https://en.wikipedia.org/wiki/Bayes%27_theorem)
  Describes the probability of an event, based on prior knowledge of conditions that might be related to the event. For example, if the risk of developing health problems is known to increase with age, Bayes' theorem allows the risk to an individual of a known age to be assessed more accurately (by conditioning it on their age) than simply assuming that the individual is typical of the population as a whole.
  Formula:

  $$
  P(A|B) = \frac{P(B|A) P(A)}{P(B)}
  $$

  Where $A$ and $B$ are events and $P(B) \neq 0$


[Markov Chain](https://en.wikipedia.org/wiki/Markov_chain)
  A Markov chain or Markov process is a stochastic model describing a sequence of possible events in which the probability of each event depends only on the state attained in the previous event.

[Markov Chain Monte Carlo](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo)
[MCMC](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo)
  Markov chain Monte Carlo (MCMC) methods comprise a class of algorithms for sampling from a probability distribution. By constructing a {term}`Markov Chain` that has the desired distribution as its equilibrium distribution, one can obtain a sample of the desired distribution by recording states from the chain.  Various algorithms exist for constructing chains, including the Metropolis–Hastings algorithm.

:::::
---
orphan: true
---

# Page not found

**Sorry, we could not find this page**

We are working on the next major release for PyMC,
which will come with faster sampling, more flexible
model building, multiple computational backends...

::::{grid} 3
:::{grid-item}
![old banner](https://raw.githubusercontent.com/pymc-devs/pymc/v3/docs/logos/svg/PyMC3_banner.svg)
:::
:::{grid-item}
:class: text-center
{fas}`arrow-alt-circle-right;fa-5x`
:::
:::{grid-item}
![new banner](https://raw.githubusercontent.com/pymc-devs/pymc/main/docs/logos/svg/PyMC_banner.svg)
:::
::::

...and much better documentation too! However, to
do so we have moved some files around and we have
modified the base url in order to support multi version docs
(you'll see the version switcher at the bottom right of the page).

Back to {{ '[documentation homepage](https://docs.pymc.io/en/{}/)'.format(version_slug) }}
(community)=
# Community

## Discourse

The [PyMC Discourse Forum](https://discourse.pymc.io/) is a great place to ask general questions about Bayesian statistics, or more specific ones about PyMC usage.

:::{tip}
Discourse is the best place to interact with the PyMC community and its team
:::

## Conferences

The PyMC community organizes [PyMCon](https://pymc-devs.github.io/pymcon/), whose last edition was in 2020.
Follow Discourse or Twitter to be updated about the next edition!

PyMC talks have been given at a number of conferences, including [PyCon](https://us.pycon.org/),
[PyData](https://pydata.org/events/), and [ODSC](https://odsc.com/) events.

## Twitter

You can also follow us on Twitter at [@pymc_devs](https://twitter.com/pymc_devs)

## Governance

Learn how the PyMC project is governed and how to become a member of PyMC [here](https://github.com/pymc-devs/pymc/blob/main/GOVERNANCE.md).

## Code of Conduct

In order to foster a welcoming environment, participation in our project and community is guided by this [code of conduct](https://github.com/pymc-devs/pymc/blob/main/CODE_OF_CONDUCT.md).
(videos_and_podcasts)=
# Videos and Podcasts

:::{card} PyMC Developers Youtube channel

[See all videos here](https://www.youtube.com/c/PyMCDevelopers/videos)
:::

:::{card} PyMC talks

Actively curated [YouTube playlist](https://www.youtube.com/playlist?list=PL1Ma_1DBbE82OVW8Fz_6Ts1oOeyOAiovy) of PyMC talks
:::

:::{card} Learning Bayesian Statistics podcast

[See all videos here](https://www.youtube.com/channel/UCAwVseuhVrpJFfik_cMHrhQ/videos)
:::
(books)=
# Books

TODO add books
(examples)=
# Notebooks

TODO add "core" notebooks:

- https://docs.pymc.io/en/stable/pymc-examples/examples/generalized_linear_models/GLM-linear.html
  - https://github.com/pymc-devs/pymc-examples/blob/main/examples/generalized_linear_models/GLM-linear.ipynb

- https://docs.pymc.io/en/stable/pymc-examples/examples/case_studies/multilevel_modeling.html
  - https://github.com/pymc-devs/pymc-examples/blob/main/examples/case_studies/multilevel_modeling.ipynb

- https://docs.pymc.io/en/stable/pymc-examples/examples/diagnostics_and_criticism/model_comparison.html
  - https://github.com/pymc-devs/pymc-examples/blob/main/examples/diagnostics_and_criticism/model_comparison.ipynb
  TODO: Add reference to ArviZ docs since it was taken from there

- https://docs.pymc.io/en/stable/Probability_Distributions.html
  - https://github.com/pymc-devs/pymc/blob/1048b69ac1410f6dec64eb1f26778d89b28a62c4/docs/source/Probability_Distributions.rst

- https://docs.pymc.io/en/stable/PyMC3_and_Theano.html
  update needed: from "pymc and theano" to "pymc and aesara" tutorial (not just auto-renaming from theano to aesara but showing correct use). Very related to migration guide as well, though the migration guide already covers it.
  - https://github.com/pymc-devs/pymc/blob/37ba9a3e3a19b738f48cb30007f4d70c33bdd0f6/docs/source/PyMC_and_Aesara.rst

- https://docs.pymc.io/en/stable/pymc-examples/examples/diagnostics_and_criticism/posterior_predictive.html
  - https://github.com/pymc-devs/pymc-examples/blob/2002ebd815a199be89b011039906b197bca42361/examples/diagnostics_and_criticism/posterior_predictive.ipynb

- https://docs.pymc.io/en/stable/Gaussian_Processes.html
  - https://github.com/pymc-devs/pymc/blob/37ba9a3e3a19b738f48cb30007f4d70c33bdd0f6/docs/source/Gaussian_Processes.rst

TODO link to the pymc-examples
TODO add categories to all notebooks

:::{toctree}
:hidden:
:maxdepth: 1

examples/pymc_overview
examples/GLM_linear
examples/model_comparison
examples/posterior_predictive
:::
# Contributing

PyMC is an open source, collective effort.
There are many ways in which you can help make it better.
And all of them are welcome!

## Contribute as an individual
PyMC is a joint effort of many people, each contributing to the areas they like
and have some expertise in, coordinating to try and cover all tasks.

Coding and documentation are the most common types of contributions, but
there are many more things that you can do to help PyMC which are just as
important. Moreover, both code and docs require submitting PRs via GitHub
to some of the repositories under the [pymc-devs](https://github.com/pymc-devs) organization, and
while we have a {ref}`pr_tutorial` guide available, GitHub might not be
everyone's cup of tea. If that is your case, don't worry, you will be
more than welcome if you want to help.

:::{tip}
Contact us on [Discourse](https://discourse.pymc.io/) if you want to contribute to the project but are not sure where you can contribute or how to start
:::

Below there are some examples of non code nor doc contributions that could serve as an inspiration.
If you have other ideas let us know on [Discourse](https://discourse.pymc.io/) to see if we can make it happen too.

* Report a bug or make a suggestion for improvement by [opening an issue in Github](https://github.com/pymc-devs/pymc/issues/new/choose)
* Answer questions on [Discourse](https://discourse.pymc.io/)
* Teach about PyMC and advertise best practices by writing blogs or giving talks
* Help plan PyMCon
* Help with outreach and marketing. This could include for example reaching out to potential sponsor
  companies, to people who could use PyMC in their work or making sure that academics who use PyMC
  cite it correctly in their work
* Help with our fundraising efforts

### Code related contributions
Join the discussion or submit a solution for an open issue. [See open issues](https://github.com/pymc-devs/pymc/issues)


### Documentation related contributions

See all open issues in documentation [here](https://github.com/pymc-devs/pymc/issues?q=is%3Aissue+is%3Aopen+label%3A%22docs%22+)

New to the open source space? Find a good beginner friendly documentation issue you can help with [here](https://github.com/pymc-devs/pymc/issues?q=is%3Aissue+is%3Aopen+label%3A%22beginner+friendly%22+label%3A%22docs%22) or [here](https://github.com/pymc-devs/pymc-examples/issues?q=is%3Aopen+label%3Adocs+label%3A%22good+first+issue%22).

You will need to express your interest in the issue by commenting on it using your GitHub account, and you will need to use some beginner level Git to complete the Pull Request process, but don't worry! There are plenty of good [Github tutorials](https://guides.github.com/activities/hello-world/) out there, and PyMC reviewers are happy to help.


## Contribute as an institution

Institutions can contribute in the following ways:

- By becoming [Institutional Partners](https://github.com/pymc-devs/pymc/blob/main/GOVERNANCE.md#institutional-partners-and-funding)
- By becoming [Sponsors](https://github.com/pymc-devs/pymc/blob/main/GOVERNANCE.md#sponsors)
- By subscribing to {ref}`Tidelift <pymc_for_enterprise>`

Contact PyMC at pymc.devs@gmail.com for more information.


TODO: Add https://github.com/pymc-devs/pymc/blob/main/CONTRIBUTING.md

:::{toctree}
:hidden:
:maxdepth: 1
:caption: Tutorials

pr_tutorial
:::

:::{toctree}
:hidden:
:maxdepth: 1
:caption: How-to guides

developer_guide_implementing_distribution
build_docs
running_the_test_suite
:::

:::{toctree}
:hidden:
:maxdepth: 1
:caption: Reference content

python_style
jupyter_style
pr_checklist
release_checklist
:::


:::{toctree}
:hidden:
:maxdepth: 1
:caption: In depth explanations

developer_guide
:::
(python_style)=
# Python style guide

## Pre commit checks

Some code-quality checks are performed during continuous integration. The easiest way to check that they pass locally,
before submitting your pull request, is by using [pre-commit](https://pre-commit.com/).

Steps to get set up are (run these within your virtual environment):

1. install:

    ```bash
    pip install pre-commit
    ```

2. enable:

    ```bash
    pre-commit install
    ```

Now, whenever you stage some file, when you run `git commit -m "<some descriptive message>"`, `pre-commit` will run
the checks defined in `.pre-commit-config.yaml` and will block your commit if any of them fail. If any hook fails, you
should fix it (if necessary), run `git add <files>` again, and then re-run `git commit -m "<some descriptive message>"`.

You can skip `pre-commit` using `--no-verify`, e.g.

```bash
git commit -m "wip lol" --no-verify
```

To skip one particular hook, you can set the `SKIP` environment variable. E.g. (on Linux):

```bash
SKIP=pyupgrade git commit -m "<descriptive message>"
```

You can manually run all `pre-commit` hooks on all files with

```bash
pre-commit run --all-files
```

or, if you just want to manually run them on a subset of files,

```bash
pre-commit run --files <file_1> <file_2> ... <file_n>
```
(pr_checklist)=
# Pull request checklist

We strongly recommended that all contribution comply with the following guidelines before being merged:

*  If your pull request addresses an issue, please use the pull request title to describe the issue and mention the issue number in the pull request description. This will make sure a link back to the original issue is created.

   :::{caution}
   Adding the related issue in the PR title generates no link and is therefore
   not useful as nobody knows issue numbers. Please mention all related
   issues in the PR but do so only in the PR description.
   :::

*  All public methods must have informative docstrings with sample usage when appropriate.

*  Please prefix the title of incomplete contributions with `[WIP]` (to indicate a work in progress). WIPs may be useful to (1) indicate you are working on something to avoid duplicated work, (2) request broad review of functionality or API, or (3) seek collaborators.

*  All other tests pass when everything is rebuilt from scratch. See {ref}`running_the_test_suite`

*  When adding additional functionality, consider adding also one example notebook at [pymc-examples](https://github.com/pymc-devs/pymc-examples). Open a [proposal issue](https://github.com/pymc-devs/pymc-examples/issues/new/choose) in the example repo to discuss the specific scope of the notebook.

* Documentation and high-coverage tests are necessary for enhancements to be accepted.

* Run any of the pre-existing examples in [pymc-examples](https://github.com/pymc-devs/pymc-examples) that contain analyses that would be affected by your changes to ensure that nothing breaks. This is a useful opportunity to not only check your work for bugs that might not be revealed by unit test, but also to show how your contribution improves PyMC for end users.

You can also check for common programming errors with the following
tools:

* Code with good test **coverage** (at least 80%), check with:

  ```bash
  $ pip install pytest pytest-cov coverage
  $ pytest --cov=pymc pymc/tests/<name of test>.py
  ```

* No `pre-commit` errors: see the {ref}`python_style` and {ref}`jupyter_style` page on how to install and run it.
(jupyter_style)=
# Jupyter Style Guide

These guidelines should be followed by notebooks in the documentation.
All notebooks in pymc-examples must follow this to the letter, the style
is more permissive for the ones on pymc where not everything is available.

The documentation websites are generated by Sphinx, which uses
{doc}`myst:index` and {doc}`myst-nb:index`
to parse the notebooks.

## General guidelines

* Don't use abbreviations or acronyms whenever you can use complete words. For example, write "random variables" instead of "RVs".

* Explain the reasoning behind each step.

* Attribute quoted text or code, and link to relevant references.

* Keep notebooks short: 20/30 cells for content aimed at beginners or intermediate users, longer notebooks are fine at the advanced level.

### MyST guidelines
Using MyST allows taking advantage of all sphinx features from markdown cells in the notebooks.
All markdown should be valid MyST (note that MyST is a superset of recommonmark).
This guide does not teach nor cover MyST extensively, only gives some opinionated guidelines.

* **Never** use url links to refer to other notebooks, PyMC documentation or other python
  libraries documentations. Use [sphinx cross-references](https://docs.readthedocs.io/en/stable/guides/cross-referencing-with-sphinx.html)
  instead.

  :::{caution}
  Using urls links breaks self referencing in versioned docs! And at the same time they are
  less robust than sphinx cross-references.
  :::

  * When linking to other notebooks, always use a `ref` type cross-reference pointing
    to the target in the {ref}`jupyter_style_first_cell`.

* If the output (or even code and output) of a cell is not necessary to follow the
  notebook or it is very long and can break the flow of reading, consider hiding
  it with a {doc}`toggle button <myst-nb:use/hiding>`

* Consider using {ref}`myst:syntax/figures` to add captions to images used in the notebook.

* Use the glossary whenever possible. If you use a term that is defined in the Glossary, link to it the first time that term appears in a significant manner. Use [this syntax](https://jupyterbook.org/content/content-blocks.html?highlight=glossary#glossaries) to add a term reference. [Link to glossary source](https://github.com/pymc-devs/pymc/blob/main/docs/source/glossary.md) where new terms should be added.

### Variable names

* Above all, stay consistent with variable names within the notebook. Notebooks using multiple names for the same variable will not be merged.

* Use meaningful variable names wherever possible. Our users come from different backgrounds and not everyone is familiar with the same naming conventions.

* Sometimes it makes sense to use Greek letters to refer to variables, for example when writing equations, as this makes them easier to read. In that case, use LaTeX to insert the Greek letter like this `$\theta$` instead of using Unicode like `θ`.

* If you need to use Greek letter variable names inside the code, please spell them out instead of using unicode. For example, `theta`, not `θ`.

* When using non meaningful names such as single letters, add bullet points with a 1-2 sentence description of each variable below the equation where they are first introduced.

Choosing variable names can sometimes be difficult, tedious or annoying.
In case it helps, the dropdown below has some suggestions so you can focus on writing the actual content

:::::::{dropdown} Variable name suggestions
:icon: light-bulb

**Models and sampling results**
* Use `idata` for sampling results, always containing a variable of type InferenceData.
* Store inferecedata groups as variables to ease writing and reading of code operating on sampling results.
  Use underscore separated 3-5 word abbreviations or the group name. Some examples of `abbrebiation`/`group_name`:
  `post`/`posterior`, `const`/`constant_data`, `post_pred`/`posterior_predictive` or `obs_data`/`observed_data`
* For stats and diagnostics, use the ArviZ function name as variable name: `ess = az.ess(...)`, `loo = az.loo(...)`
* If there are multiple models in a notebook, assign a prefix to each model,
  and use it throughout to identify which variables map to each model.
  Taking the famous eight school as example, with a `centered` and `non_centered` model
  to compare parametrizations, use `centered_model` (pm.Model object), `centered_idata`, `centered_post`, `centered_ess`... and `non_centered_model`, `non_centered_idata`...

**Dimension and random variable names**
* Use singular dimension names, following ArviZ `chain` and `draw`.
  For example `cluster`, `axis`, `component`, `forest`, `time`...
* If you can't think of a meaningful name for the dimension representing the number of observations such as time, fall back to `obs_id`.
* For matrix dimensions, as xarray doesn't allow repeated dimension names, add a `_bis` suffix. i.e. `param, param_bis`.
* For the dimension resulting from stacking `chain` and `draw` use `sample`, that is `.stack(sample=("chain", "draw"))`.
* We often need to encode a categorical variable as integers. add `_idx` to the name of the variable it's encoding.
  i.e. from `floor` and `county` to `floor_idx` and `county_idx`.
* To avoid clashes and overwriting variables when using `pm.Data`, use the following pattern:

  ```
  x = np.array(...)
  with pm.Model():
      x_ = pm.Data("x", x)
      ...
  ```

  This avoids overwriting the original `x` while having `idata.constant_data["x"]`,
  and within the model `x_` is still available to play the role of `x`.
  Otherwise, always try to use the same variable name as the string name given to the PyMC random variable.

**Plotting**
* Matplotlib figures and axes. Use:
  * `fig` for matplotlib figures
  * `ax` for a single matplotib axes object
  * `axs` for arrays of matplotlib axes objects

  When manually working with multiple matplotlib axes, use local `ax` variables:

  ::::{tab-set}
  :::{tab-item} Local `ax` variables
  ```{code-block} python
  :emphasize-lines: 3, 7
  fig, axs = pyplot.subplots()

  ax = axs[0, 1]
  ax.plot(...)
  ax.set(...)

  ax = axs[1, 2]
  ax.scatter(...)
  ```
  :::
  :::{tab-item} Instead of subsetting every time
  ```
  fig, axs = pyplot.subplots()

  axs[0, 1].plot(...)
  axs[0, 1].set(...)

  axs[1. 2].scatter(...)
  ```
  :::
  ::::

  This makes editing the code if restructuring the subplots easier, only one change per subplot
  is needed instead of one change per matplotlib function call.

* It is often useful to make a numpy linspace into an {class}`~xarray.DataArray`
  for xarray to handle aligning and broadcasing automatically and ease computation.
  * If a dimension name is needed, use `x_plot`
  * If a variable name is needed for the original array and DataArray to coexist, add `_da` suffix

  Thus, ending up with code like:

  ```
  x = xr.DataArray(np.linspace(0, 10, 100), dims=["x_plot"])
  # or
  x = np.linspace(0, 10, 100)
  x_da = xr.DataArray(x)
  ```

**Looping**
* When using enumerate, take the first letter of the variable as the count:

  ```
  for p, person in enumerate(persons)
  ```

* When looping, if you need to store a variable after subsetting with the loop index,
  append the index variable used for looping to the original variable name:

  ```{code-block} python
  :emphasize-lines: 4, 6
  variable = np.array(...)
  x = np.array(...)
  for i in range(N):
      variable_i = variable[i]
      for j in range(K):
          x_j = x[j]
          ...
  ```

:::::::

(jupyter_style_first_cell)=
## First cell
The first cell of all example notebooks should have a MyST target, a level 1 markdown title (that is a title with a single `#`) followed by the post directive.
The syntax is as follows:

```markdown
(notebook_id)=
# Notebook Title

:::{post} Aug 31, 2021
:tags: tag1, tag2, tags can have spaces, tag4
:category: level
:author: Alice Abat, Bob Barceló
:::
```

The date should correspond to the latest update/execution date, at least roughly (it's not a problem if the date is a few days off due to the review process before merging the PR). This will allow users to see which notebooks have been updated lately and will help the PyMC team make sure no notebook is left outdated for too long.

The {ref}`MyST target <myst:syntax/targets>`
is important to ease referencing and linking notebooks between each other.

Tags can be anything, but we ask you to try to use [existing tags](https://github.com/pymc-devs/pymc/wiki/Categories-and-Tags-for-PyMC-Examples)
to avoid the tag list from getting too long.

Each notebook should have a single category indicating the level of the notebook.
Choose a category from [existing categories](https://github.com/pymc-devs/pymc/wiki/Categories-and-Tags-for-PyMC-Examples#categories).

Authors should list people who authored, adapted or updated the notebook. See {ref}`jupyter_authors`
for more details.

## Extra dependencies
If the notebook uses libraries that are not PyMC dependencies, these extra dependencies should
be indicated together with some advise on how to install them.
This ensures readers know what they'll need to install beforehand and can for example
decide between running it locally or on binder.

To make things easier for notebook writers and maintainers, pymc-examples contains
a template for this that warns about the extra dependencies and provides specific
installation instructions inside a dropdown.

Thus, notebooks with extra dependencies should:

1.  list the extra dependencies as notebook metadata using the `myst_substitutions` category
    and then either the `extra_dependencies` or the `pip_dependencies` and `conda_dependencies`.
    In addition, there is also an `extra_install_notes` to include custom text inside the dropdown.

    * notebook metadata can be edited from the menubar `Edit` -> `Edit notebook metadata`
      in the dropdown

      This will open a window with json formatted text that might look a bit like:

      ::::{tab-set}
      :::{tab-item} No myst_substitutions

      ```json
      {
        "kernelspec": {
          "name": "python3",
          "display_name": "Python 3 (ipykernel)",
          "language": "python"
        },
        "language_info": {
          "name": "python",
          "version": "3.9.7",
          "mimetype": "text/x-python",
          "codemirror_mode": {
            "name": "ipython",
            "version": 3
          },
          "pygments_lexer": "ipython3",
          "nbconvert_exporter": "python",
          "file_extension": ".py"
        }
      }
      ```
      :::

      :::{tab-item} extra_dependencies key

      ```{code-block} json
      :emphasize-lines: 19-21
      {
        "kernelspec": {
          "name": "python3",
          "display_name": "Python 3 (ipykernel)",
          "language": "python"
        },
        "language_info": {
          "name": "python",
          "version": "3.9.7",
          "mimetype": "text/x-python",
          "codemirror_mode": {
            "name": "ipython",
            "version": 3
          },
          "pygments_lexer": "ipython3",
          "nbconvert_exporter": "python",
          "file_extension": ".py"
        },
        "myst_substitutions": {
          "extra_dependencies": "bambi seaborn"
        }
      }
      ```
      :::

      :::{tab-item} pip and conda specific keys
      ```{code-block} json
      :emphasize-lines: 19-22
      {
        "kernelspec": {
          "name": "python3",
          "display_name": "Python 3 (ipykernel)",
          "language": "python"
        },
        "language_info": {
          "name": "python",
          "version": "3.9.7",
          "mimetype": "text/x-python",
          "codemirror_mode": {
            "name": "ipython",
            "version": 3
          },
          "pygments_lexer": "ipython3",
          "nbconvert_exporter": "python",
          "file_extension": ".py"
        },
        "myst_substitutions": {
          "pip_dependencies": "graphviz",
          "conda_dependencies": "python-graphviz",
        }
      }
      ```

      The pip and conda spcific keys overwrite the `extra_installs` one, so it doesn't make
      sense to use `extra_installs` if using them. Either both pip and conda substitutions
      are defined or none of them is.
      :::
      ::::

1.  include the warning and installation advise template with the following markdown right before
    the extra dependencies are imported:

    ```markdown
    :::{include} ../extra_installs.md
    :::
    ```

## Code preamble

In a cell just below the cell where you imported matplotlib and/or ArviZ (usually the first one),
set the ArviZ style to darkgrid (this has to be in another cell than the matplotlib import because of the way matplotlib sets its defaults):

```python
RANDOM_SEED = 8927
rng = np.random.default_rng(RANDOM_SEED)
az.style.use("arviz-darkgrid")
```

A good practice _when generating synthetic data_ is also to set a random seed as above, to improve reproducibility. Also, please check convergence (e.g. `assert all(r_hat < 1.03)`) because we sometime re-run notebooks automatically without carefully checking each one.

## Reading from file

Use a `try... except` clause to load the data and use `pm.get_data` in the except path. This will ensure that users who have cloned pymc-examples repo will read their local copy of the data while also downloading the data from github for those who don't have a local copy. Here is one example:

```python
try:
    df_all = pd.read_csv(os.path.join("..", "data", "file.csv"), ...)
except FileNotFoundError:
    df_all = pd.read_csv(pm.get_data("file.csv"), ...)
```

## pre-commit and code formatting
We run some code-quality checks on our notebooks during Continuous Integration. The easiest way to make sure your notebook(s) pass the CI checks is using [pre-commit](https://github.com/pre-commit/pre-commit). You can install it with

```bash
pip install -U pre-commit
```

and then enable it with

```bash
pre-commit install
```

Then, the code-quality checks will run automatically whenever you commit any changes. To run the code-quality checks manually, you can do, e.g.:

```bash
pre-commit run --files notebook1.ipynb notebook2.ipynb
```

replacing `notebook1.ipynb` and `notebook2.ipynb` with any notebook you've modified.

NB: sometimes, [Black will be frustrating](https://stackoverflow.com/questions/58584413/black-formatter-ignore-specific-multi-line-code/58584557#58584557) (well, who isn't?). In these cases, you can disable its magic for specific lines of code: just write `#fmt: on/off` to disable/re-enable it, like this:

```python
# fmt: off
np.array(
    [
        [1, 0, 0, 0],
        [0, -1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, -1],
    ]
)
# fmt: on
```

(jupyter_authors)=
## Authorship and attribution
After the notebook content finishes, there should be an `## Authors` section with bullet points
to provide attribution to the people who contributed to the the general pattern should be:

```markdown
* <verb> by <author> on <date> ([repo#PR](https://link-to.pr))
```

where `<verb>` must be one listed below, `<author>` should be the name (multiple people allowed)
which can be formatted as hyperlink to personal site or GitHub profile of the person,
and `<date>` should preferably be month and year.

authored
: for notebooks created specifically for pymc-examples

adapted
: for notebooks adapted from other sources such as books or blogposts.
  It will therefore follow a different structure than the example above
  in order to include a link or reference to the original source:

  ```markdown
  Adapted from Alice's [blogpost](blog.alice.com) by Bob and Carol on ...
  ```

re-executed
: for notebooks re-executed with a newer PyMC version without significant changes to the code.
  It can also mention the PyMC version used to run the notebook.

updated
: for notebooks that have not only been re-executed but have also had significant updates to
  their content (either code, explanations or both).

some examples:

```markdown
* Authored by Chris Fonnesbeck in May, 2017 ([pymc#2124](https://github.com/pymc-devs/pymc/pull/2124))
* Updated by Colin Carroll in June, 2018 ([pymc#3049](https://github.com/pymc-devs/pymc/pull/3049))
* Updated by Alex Andorra in January, 2020 ([pymc#3765](https://github.com/pymc-devs/pymc/pull/3765))
* Updated by Oriol Abril in June, 2020 ([pymc#3963](https://github.com/pymc-devs/pymc/pull/3963))
* Updated by Farhan Reynaldo in November 2021 ([pymc-examples#246](https://github.com/pymc-devs/pymc-examples/pull/246))
```

and

```markdown
* Adapted from chapter 5 of Bayesian Data Analysis 3rd Edition {cite:p}`gelman2013bayesian`
  by Demetri Pananos and Junpeng Lao on June, 2018 ([pymc#3054](https://github.com/pymc-devs/pymc/pull/3054))
* Reexecuted by Ravin Kumar with PyMC 3.6 on March, 2019 ([pymc#3397](https://github.com/pymc-devs/pymc/pull/3397))
* Reexecuted by Alex Andorra and Michael Osthege with PyMC 3.9 on June, 2020 ([pymc#3955](https://github.com/pymc-devs/pymc/pull/3955))
* Updated by Raúl Maldonado 2021 ([pymc-examples#24](https://github.com/pymc-devs/pymc-examples/pull/24), [pymc-examples#45](https://github.com/pymc-devs/pymc-examples/pull/45) and [pymc-examples#147](https://github.com/pymc-devs/pymc-examples/pull/147))
```

## References
References should be added to the [`references.bib`](https://github.com/pymc-devs/pymc-examples/blob/main/examples/references.bib) file in bibtex format, and cited with [sphinxcontrib-bibtex](https://sphinxcontrib-bibtex.readthedocs.io/en/latest/) within the notebook text wherever they are relevant.

The references in the `.bib` file should have as id something along the lines `authorlastnameYEARkeyword` or `libraryYEARkeyword` for documentation pages, and they should be alphabetically sorted by this id in order to ease finding references within the file and preventing adding duplicate ones.

References can be cited twice within a single notebook. Two common reference formats are:

```
{cite:p}`bibtex_id`  # shows the reference author and year between parenthesis
{cite:t}`bibtex_id`  # textual cite, shows author and year without parenthesis
```

which can be added inline, within the text itself. At the end of the notebook, add the bibliography with the following markdown

```markdown
## References

:::{bibliography}
:filter: docname in docnames
:::
```

or alternatively, if you wanted to add extra references that have not been cited within the text, use:

```markdown
## References

:::{bibliography}
:filter: docname in docnames

extra_bibtex_id_1
extra_bibtex_id_2
:::
```

## Watermark
Once you're finished with your NB, add a very last cell with [the watermark package](https://github.com/rasbt/watermark). This will automatically print the versions of Python and the packages you used to run the NB -- reproducibility rocks! Here is some example code. Note that the `-p` argument may not be necessary (or it may need to have different libraries as input), but all the other arguments must be present.

```python
%load_ext watermark
%watermark -n -u -v -iv -w -p aesara,aeppl,xarray
```

This second to last code cell should be preceded by a markdown cell with the `## Watermark` title only so it appears in the table of contents.

`watermark` should be in your virtual environment if you installed our `requirements-dev.txt`.
Otherwise, just run `pip install watermark`.
The `p` flag is optional but should be added if Aesara or xarray are not imported explicitly.
This will also be checked by `pre-commit` (because we all forget to do things sometimes 😳).

## Epilogue
The last cell in the notebooks should be a markdown cell with exactly the following content:

```
:::{include} ../page_footer.md
:::
```

The only exception being notebooks that are not on the usual place and therefore need to
update the path to page footer for the include to work.

---

You're all set now 🎉 You can push your changes, open a pull request, and, once it's merged, rest with the feeling of a job well done 👏
Thanks a lot for your contribution to open-source, we really appreciate it!
# PyMC Release workflow
+ Track all relevant issues and PRs via a **version-specific [milestone](https://github.com/pymc-devs/pymc/milestones)**
+ Make sure that there are no major known bugs that should not be released
+ Make a PR to **bump the version number** in `__init__.py` and edit the `RELEASE-NOTES.md`.
  + :::{important}
    Please don't name it after the release itself, and remember to push to your own fork like an ordinary citizen.
    :::
  + Create a new "vNext" section at the top
  + Edit the header with the release version and date
  + Add a line to credit the release manager like in previous releases
+ After merging the PR, check that the CI pipelines on master are all ✔
+ Create a Release with the Tag as ´v1.2.3´ and a human-readable title like the ones on previous releases

After the last step, the [GitHub Action "release-pipeline"](https://github.com/pymc-devs/pymc/blob/master/.github/workflows/release.yml) triggers and automatically builds and publishes the new version to PyPI.

## Troubleshooting
+ If for some reason, the release must be "unpublished", this is possible by manually deleting it on PyPI and GitHub. HOWEVER, PyPI will not accept another release with the same version number!
+ The `release-pipeline` has a `test-install-job`, which can fail if the PyPI index did not update fast enough.

## Post-release steps
+ Head over to [Zenodo](https://zenodo.org/record/4603970) and copy the version specific DOI-bade into the [release notes](https://github.com/pymc-devs/pymc/releases)
+ Rename and close the release milestone and open a new "vNext" milestone
+ Monitor the update the [conda-forge/pymc-feedstock](https://github.com/conda-forge/pymc-feedstock) repository for new PRs. The bots should automatically pick up the new version and open a PR to update it. Manual intervention may be required though (see the repos PR history for examples).
+ Re-run notebooks with the new release (see https://github.com/pymc-devs/pymc-examples)
+ Make sure the new version appears at the website and that [`docs.pymc.io/en/stable`](https://docs.pymc.io/en/stable) points to it.
# Running the test suite

TODO: explain steps to run test suite. This is a how to guide, so assume readers know how to install
things and so on, at most mention what dependencies are needed.
# Build documentation locally

To build the docs, run these commands at pymc repository root:

```bash
$ pip install -r requirements-dev.txt  # Make sure the dev requirements are installed
$ make clean  # clean built docs from previous runs and intermediate outputs
$ make html   # Build docs
$ python -m http.server --directory ../_build/  # Render docs
```

Check the printed url where docs are being served and open it.

The `make clean` step is not always necessary, if you are working on a specific page
for example, you can rebuild the docs without the clean step and everything should
work fine. If you are restructuring the content or editing toctrees, then you'll need
to execute `make clean`.

A good approach is to skip the `make clean`, which makes
the `make html` blazing fast and see how everything looks.
If something looks strange, run `make clean` and `make html` one after the other
to see if it fixes the issue before checking anything else.
# Implementing a Distribution

This guide provides an overview on how to implement a distribution for version 4 of PyMC.
It is designed for developers who wish to add a new distribution to the library.
Users will not be aware of all this complexity and should instead make use of helper methods such as `~pymc.distributions.DensityDist`.

PyMC {class}`~pymc.distributions.Distribution` builds on top of Aesara's {class}`~aesara.tensor.random.op.RandomVariable`, and implements `logp`, `logcdf` and `get_moment` methods as well as other initialization and validation helpers.
Most notably `shape/dims` kwargs, alternative parametrizations, and default `transforms`.

Here is a summary check-list of the steps needed to implement a new distribution.
Each section will be expanded below:

1. Creating a new `RandomVariable` `Op`
1. Implementing the corresponding `Distribution` class
1. Adding tests for the new `RandomVariable`
1. Adding tests for `logp` / `logcdf` and `get_moment` methods
1. Documenting the new `Distribution`.

This guide does not attempt to explain the rationale behind the `Distributions` current implementation, and details are provided only insofar as they help to implement new "standard" distributions.

## 1. Creating a new `RandomVariable` `Op`

{class}`~aesara.tensor.random.op.RandomVariable` are responsible for implementing the random sampling methods, which in version 3 of PyMC used to be one of the standard `Distribution` methods, alongside `logp` and `logcdf`.
The `RandomVariable` is also responsible for parameter broadcasting and shape inference.

Before creating a new `RandomVariable` make sure that it is not offered in the [Numpy library](https://numpy.org/doc/stable/reference/random/generator.html#distributions).
If it is, it should be added to the [Aesara library](https://github.com/aesara-devs/aesara) first and then imported into the PyMC library.

In addition, it might not always be necessary to implement a new `RandomVariable`.
For example if the new `Distribution` is just a special parametrization of an existing `Distribution`.
This is the case of the `OrderedLogistic` and `OrderedProbit`, which are just special parametrizations of the `Categorical` distribution.

The following snippet illustrates how to create a new `RandomVariable`:

```python

from aesara.tensor.var import TensorVariable
from aesara.tensor.random.op import RandomVariable
from typing import List, Tuple

# Create your own `RandomVariable`...
class BlahRV(RandomVariable):
    name: str = "blah"

    # Provide the minimum number of (output) dimensions for this RV
    # (e.g. `0` for a scalar, `1` for a vector, etc.)
    ndim_supp: int = 0

    # Provide the number of (input) dimensions for each parameter of the RV
    # (e.g. if there's only one vector parameter, `[1]`; for two parameters,
    # one a matrix and the other a scalar, `[2, 0]`; etc.)
    ndims_params: List[int] = [0, 0]

    # The NumPy/Aesara dtype for this RV (e.g. `"int32"`, `"int64"`).
    # The standard in the library is `"int64"` for discrete variables
    # and `"floatX"` for continuous variables
    dtype: str = "floatX"

    # A pretty text and LaTeX representation for the RV
    _print_name: Tuple[str, str] = ("blah", "\\operatorname{blah}")

    # If you want to add a custom signature and default values for the
    # parameters, do it like this. Otherwise this can be left out.
    def __call__(self, loc=0.0, scale=1.0, size=None, **kwargs) -> TensorVariable:
        return super().__call__(loc, scale, size=size, **kwargs)

    # This is the Python code that produces samples.  Its signature will always
    # start with a NumPy `RandomState` object, then the distribution
    # parameters, and, finally, the size.
    #
    # This is effectively the PyMC v4.x replacement for `Distribution.random`.
    @classmethod
    def rng_fn(
        cls,
        rng: np.random.RandomState,
        loc: np.ndarray,
        scale: np.ndarray,
        size: Tuple[int, ...],
    ) -> np.ndarray:
        return scipy.stats.blah.rvs(loc, scale, random_state=rng, size=size)

# Create the actual `RandomVariable` `Op`...
blah = BlahRV()

```

Some important things to keep in mind:

1. Everything inside the `rng_fn` method is pure Python code (as are the inputs) and should not make use of other `Aesara` symbolic ops. The random method should make use of the `rng` which is a Numpy    {class}`~numpy.random.RandomState`, so that samples are reproducible.
1. The `size` argument (together with the inputs shape) are the only way for the user to specify non-default `RandomVariable` dimensions. The `rng_fn` will have to take this into consideration for correct output. `size` is the specification used by `Numpy` and `Scipy` and works like PyMC `shape` for univariate distributions, but is different for multivariate distributions. Unfortunately there is no general reference documenting how `size` ought to work for multivariate distributions. This [discussion](https://github.com/numpy/numpy/issues/17669) may be helpful to get more context.
1. `Aesara` tries to infer the output shape of the `RandomVariable` (given a user-specified size) by introspection of the `ndim_supp` and `ndim_params` attributes. However, the default method may not work for more complex distributions. In that case, custom `_shape_from_params` (and less probably, `_infer_shape`) should also be implemented in the new `RandomVariable` class. One simple example is seen in the {class}`~pymc.distributions.multivariate.DirichletMultinomialRV` where it was necessary to specify the `rep_param_idx` so that the `default_shape_from_params` helper method could do its job. In more complex cases, it may not be possible to make use of the default helper, but those have not been found yet!
1. It's okay to use the `rng_fn` `classmethods` of other Aesara and PyMC `RandomVariables` inside the new `rng_fn`. For example if you are implementing a negative HalfNormal `RandomVariable`, your `rng_fn` can simply return `- halfnormal.rng_fn(rng, scale, size)`.

*Note: In addition to `size`, the PyMC API also provides `shape` and `dims` as alternatives to define a distribution dimensionality, but this is taken care of by {class}`~pymc.distributions.Distribution`, and should not require any extra changes.*

For a quick test that your new `RandomVariable` `Op` is working, you can call the `Op` with the necessary parameters and then call `eval()` on the returned object:

```python

# blah = aesara.tensor.random.uniform in this example
blah([0, 0], [1, 2], size=(10, 2)).eval()

# array([[0.83674527, 0.76593773],
#    [0.00958496, 1.85742402],
#    [0.74001876, 0.6515534 ],
#    [0.95134629, 1.23564938],
#    [0.41460156, 0.33241175],
#    [0.66707807, 1.62134924],
#    [0.20748312, 0.45307477],
#    [0.65506507, 0.47713784],
#    [0.61284429, 0.49720329],
#    [0.69325978, 0.96272673]])

```

## 2. Inheriting from a PyMC base `Distribution` class

After implementing the new `RandomVariable` `Op`, it's time to make use of it in a new PyMC {class}`pymc.distributions.Distribution`.
PyMC 4.x works in a very {term}`functional <Functional Programming>` way, and the `distribution` classes are there mostly to facilitate porting the `PyMC3` v3.x code to the new `PyMC` v4.x version, add PyMC API features and keep related methods organized together.
In practice, they take care of:

1. Linking ({term}`Dispatching`) a rv_op class with the corresponding `get_moment`, `logp` and `logcdf` methods.
1. Defining a standard transformation (for continuous distributions) that converts a bounded variable domain (e.g., positive line) to an unbounded domain (i.e., the real line), which many samplers prefer.
1. Validating the parametrization of a distribution and converting non-symbolic inputs (i.e., numeric literals or numpy arrays) to symbolic variables.
1. Converting multiple alternative parametrizations to the standard parametrization that the `RandomVariable` is defined in terms of.

Here is how the example continues:

```python

from pymc.aesaraf import floatX, intX
from pymc.distributions.continuous import PositiveContinuous
from pymc.distributions.dist_math import check_parameters


# Subclassing `PositiveContinuous` will dispatch a default `log` transformation
class Blah(PositiveContinuous):
    # This will be used by the metaclass `DistributionMeta` to dispatch the
    # class `logp` and `logcdf` methods to the `blah` `op`
    rv_op = blah

    # dist() is responsible for returning an instance of the rv_op. We pass
    # the standard parametrizations to super().dist
    @classmethod
    def dist(cls, param1, param2=None, alt_param2=None, **kwargs):
        param1 = at.as_tensor_variable(intX(param1))
        if param2 is not None and alt_param2 is not None:
            raise ValueError('Only one of param2 and alt_param2 is allowed')
        if alt_param2 is not None:
            param2 = 1 / alt_param2
        param2 = at.as_tensor_variable(floatX(param2))

        # The first value-only argument should be a list of the parameters that
        # the rv_op needs in order to be instantiated
        return super().dist([param1, param2], **kwargs)

    # get_moment returns a symbolic expression for the stable moment from which to start sampling
    # the variable, given the implicit `rv`, `size` and `param1` ... `paramN`
    def get_moment(rv, size, param1, param2):
        moment, _ = at.broadcast_arrays(param1, param2)
        if not rv_size_is_none(size):
            moment = at.full(size, moment)
        return moment

    # Logp returns a symbolic expression for the logp evaluation of the variable
    # given the `value` of the variable and the parameters `param1` ... `paramN`
    def logp(value, param1, param2):
        logp_expression = value * (param1 + at.log(param2))

        # A switch is often used to enforce the distribution support domain
        bounded_logp_expression = at.switch(
            at.gt(value >= 0),
            logp_expression,
            -np.inf,
        )

        # We use `check_parameters` for parameter validation. After the default expression,
        # multiple comma-separated symbolic conditions can be added. Whenever
        # a bound is invalidated, the returned expression raises an error with the message
        # defined in the optional `msg` keyword argument.
        return check_parameters(
            logp_expression,
            param2 >= 0,
            msg="param2 >= 0",
    )

    # logcdf works the same way as logp. For bounded variables, it is expected to return
    # `-inf` for values below the domain start and `0` for values above the domain end.
    def logcdf(value, param1, param2):
        ...

```

Some notes:

1. A distribution should at the very least inherit from {class}`~pymc.distributions.Discrete` or {class}`~pymc.distributions.Continuous`. For the latter, more specific subclasses exist: `PositiveContinuous`, `UnitContinuous`, `BoundedContinuous`, `CircularContinuous`, which specify default transformations for the variables. If you need to specify a one-time custom transform you can also override the `__new__` method, as is done for the {class}`~pymc.distributions.multivariate.Dirichlet`.
1. If a distribution does not have a corresponding `random` implementation, a `RandomVariable` should still be created that raises a `NotImplementedError`. This is the case for the {class}`~pymc.distributions.continuous.Flat`. In this case it will be necessary to provide a standard `initval` by
   overriding `__new__`.
1. As mentioned above, `PyMC` v4.x works in a very {term}`functional <Functional Programming>` way, and all the information that is needed in the `logp` and `logcdf` methods is expected to be "carried" via the `RandomVariable` inputs. You may pass numerical arguments that are not strictly needed for the `rng_fn` method but are used in the `logp` and `logcdf` methods. Just keep in mind whether this affects the correct shape inference behavior of the `RandomVariable`. If specialized non-numeric information is needed you might need to define your custom`_logp` and `_logcdf` {term}`Dispatching` functions, but this should be done as a last resort.
1. The `logcdf` method is not a requirement, but it's a nice plus!
1. Currently only one moment is supported in the `get_moment` method, and probably the "higher-order" one is the most useful (that is `mean` > `median` > `mode`)... You might need to truncate the moment if you are dealing with a discrete distribution.
1. When creating the `get_moment` method, we have to be careful with `size != None` and broadcast properly when some parameters that are not used in the moment may nevertheless inform about the shape of the distribution. E.g. `pm.Normal.dist(mu=0, sigma=np.arange(1, 6))` returns a moment of `[mu, mu, mu, mu, mu]`.

For a quick check that things are working you can try the following:

```python

import pymc as pm
from pymc.distributions.distribution import get_moment

# pm.blah = pm.Normal in this example
blah = pm.blah.dist(mu = 0, sigma = 1)

# Test that the returned blah_op is still working fine
blah.eval()
# array(-1.01397228)

# Test the get_moment method
get_moment(blah).eval()
# array(0.)

# Test the logp method
pm.logp(blah, [-0.5, 1.5]).eval()
# array([-1.04393853, -2.04393853])

# Test the logcdf method
pm.logcdf(blah, [-0.5, 1.5]).eval()
# array([-1.17591177, -0.06914345])
```

## 3. Adding tests for the new `RandomVariable`

Tests for new `RandomVariables` are mostly located in `pymc/tests/test_distributions_random.py`.
Most tests can be accommodated by the default `BaseTestDistribution` class, which provides default tests for checking:
1. Expected inputs are passed to the `rv_op` by the `dist` `classmethod`, via `check_pymc_params_match_rv_op`
1. Expected (exact) draws are being returned, via `check_pymc_draws_match_reference`
1. Shape variable inference is correct, via `check_rv_size`

```python

class TestBlah(BaseTestDistribution):

    pymc_dist = pm.Blah
    # Parameters with which to test the blah pymc Distribution
    pymc_dist_params = {"param1": 0.25, "param2": 2.0}
    # Parameters that are expected to have passed as inputs to the RandomVariable op
    expected_rv_op_params = {"param1": 0.25, "param2": 2.0}
    # If the new `RandomVariable` is simply calling a `numpy`/`scipy` method,
    # we can make use of `seeded_[scipy|numpy]_distribution_builder` which
    # will prepare a seeded reference distribution for us.
    reference_dist_params = {"mu": 0.25, "loc": 2.0}
    reference_dist = seeded_scipy_distribution_builder("blah")
    tests_to_run = [
        "check_pymc_params_match_rv_op",
        "check_pymc_draws_match_reference",
        "check_rv_size",
    ]
```

Additional tests should be added for each optional parametrization of the distribution.
In this case it's enough to include the test `check_pymc_params_match_rv_op` since only this differs.

Make sure the tested alternative parameter value would lead to a different value for the associated default parameter.
For instance, if it's just the inverse, testing with `1.0` is not very informative, since the conversion would return `1.0` as well, and we can't be (as) sure that is working correctly.

```python

class TestBlahAltParam2(BaseTestDistribution):

    pymc_dist = pm.Blah
    # param2 is equivalent to 1 / alt_param2
    pymc_dist_params = {"param1": 0.25, "alt_param2": 4.0}
    expected_rv_op_params = {"param1": 0.25, "param2": 2.0}
    tests_to_run = ["check_pymc_params_match_rv_op"]

```

Custom tests can also be added to the class as is done for the {class}`~pymc.tests.test_random.TestFlat`.

### Note on `check_rv_size` test:

Custom input sizes (and expected output shapes) can be defined for the `check_rv_size` test, by adding the optional class attributes `sizes_to_check` and `sizes_expected`:

```python
sizes_to_check = [None, (1), (2, 3)]
sizes_expected = [(3,), (1, 3), (2, 3, 3)]
tests_to_run = ["check_rv_size"]
```

This is usually needed for Multivariate distributions.
You can see an example in {class}`~pymc.test.test_random.TestDirichlet`.

### Notes on `check_pymcs_draws_match_reference` test

The `check_pymcs_draws_match_reference` is a very simple test for the equality of draws from the `RandomVariable` and the exact same python function, given the same inputs and random seed.
A small number (`size=15`) is checked. This is not supposed to be a test for the correctness of the random generator.
The latter kind of test (if warranted) can be performed with the aid of `pymc_random` and `pymc_random_discrete` methods in the same test file, which will perform an expensive statistical comparison between the RandomVariable `rng_fn` and a reference Python function.
This kind of test only makes sense if there is a good independent generator reference (i.e., not just the same composition of numpy / scipy python calls that is done inside `rng_fn`).

Finally, when your `rng_fn` is doing something more than just calling a `numpy` or `scipy` method, you will need to setup an equivalent seeded function with which to compare for the exact draws (instead of relying on `seeded_[scipy|numpy]_distribution_builder`).
You can find an example in {class}`~pymc.tests.test_distributions_random.TestWeibull`, whose `rng_fn` returns `beta * np.random.weibull(alpha, size=size)`.


## 4. Adding tests for the `logp` / `logcdf` methods

Tests for the `logp` and `logcdf` methods are contained in `pymc/tests/test_distributions.py`, and most make use of the `TestMatchesScipy` class, which provides `check_logp`, `check_logcdf`, and
`check_selfconsistency_discrete_logcdf` standard methods.
These will suffice for most distributions.

```python

from pymc.tests.helpers import select_by_precision

R = Domain([-np.inf, -2.1, -1, -0.01, 0.0, 0.01, 1, 2.1, np.inf])
Rplus = Domain([0, 0.01, 0.1, 0.9, 0.99, 1, 1.5, 2, 100, np.inf])

...

def test_blah(self):

  self.check_logp(
      pymc_dist=pm.Blah,
      # Domain of the distribution values
      domain=R,
      # Domains of the distribution parameters
      paramdomains={"mu": R, "sigma": Rplus},
      # Reference scipy (or other) logp function
      scipy_logp = lambda value, mu, sigma: sp.norm.logpdf(value, mu, sigma),
      # Number of decimal points expected to match between the pymc and reference functions
      decimal=select_by_precision(float64=6, float32=3),
      # Maximum number of combinations of domain * paramdomains to test
      n_samples=100,
  )

  self.check_logcdf(
      pymc_dist=pm.Blah,
      domain=R,
      paramdomains={"mu": R, "sigma": Rplus},
      scipy_logcdf=lambda value, mu, sigma: sp.norm.logcdf(value, mu, sigma),
      decimal=select_by_precision(float64=6, float32=1),
      n_samples=-1,
  )

```

These methods will perform a grid evaluation on the combinations of domain and paramdomains values, and check that the pymc methods and the reference functions match.
There are a couple of details worth keeping in mind:

1. By default the first and last values (edges) of the `Domain` are not compared (they are used for other things). If it is important to test the edge of the `Domain`, the edge values can be repeated. This is done by the `Bool`: `Bool = Domain([0, 0, 1, 1], "int64")`
1. There are some default domains (such as `R` and `Rplus`) that you can use for testing your new distribution, but it's also perfectly fine to create your own domains inside the test function if there is a good reason for it (e.g., when the default values lead too many extreme unlikely combinations that are not very informative about the correctness of the implementation).
1. By default, a random subset of 100 `param` x `paramdomain` combinations is tested, in order to keep the test runtime under control. When testing your shiny new distribution, you can temporarily set `n_samples=-1` to force all combinations to be tested. This is important to avoid the your `PR` leading to surprising failures in future runs whenever some bad combinations of parameters are randomly tested.
1. On Github all tests are run twice, under the `aesara.config.floatX` flags of `"float64"` and `"float32"`. However, the reference Python functions will run in a pure "float64" environment, which means the reference and the PyMC results can diverge quite a lot (e.g., underflowing to `-np.inf` for extreme parameters). You should therefore make sure you test locally in both regimes. A quick and dirty way of doing this is to temporariliy add `aesara.config.floatX = "float32"` at the very top of file, immediately after `import aesara`. Remember to set `n_samples=-1` as well to test all combinations. The test output will show what exact parameter values lead to a failure. If you are confident that your implementation is correct, you may opt to tweak the decimal precision with `select_by_precision`, or adjust the tested `Domain` values. In extreme cases, you can mark the test with a conditional `xfail` (if only one of the sub-methods is failing, they should be separated, so that the `xfail` is as narrow as possible):

```python

def test_blah_logp(self):
    ...


@pytest.mark.xfail(
   condition=(aesara.config.floatX == "float32"),
   reason="Fails on float32 due to numerical issues",
)
def test_blah_logcdf(self):
    ...


```

## 5. Adding tests for the `get_moment` method

Tests for the `get_moment` method are contained in `pymc/tests/test_distributions_moments.py`, and make use of the function `assert_moment_is_expected`
which checks if:
1. Moments return the `expected` values
1. Moments have the expected size and shape

```python

import pytest
from pymc.distributions import Blah

@pytest.mark.parametrize(
    "param1, param2, size, expected",
    [
        (0, 1, None, 0),
        (0, np.ones(5), None, np.zeros(5)),
        (np.arange(5), 1, None, np.arange(5)),
        (np.arange(5), np.arange(1, 6), (2, 5), np.full((2, 5), np.arange(5))),
    ],
)
def test_blah_moment(param1, param2, size, expected):
    with Model() as model:
        Blah("x", param1=param1, param2=param2, size=size)
    assert_moment_is_expected(model, expected)

```

Here are some details worth keeping in mind:

1. In the case where you have to manually broadcast the parameters with each other it's important to add test conditions that would fail if you were not to do that. A straightforward way to do this is to make the used parameter a scalar, the unused one(s) a vector (one at a time) and size `None`.
1. In other words, make sure to test different combinations of size and broadcasting to cover these cases.

## 6. Documenting the new `Distribution`

New distributions should have a rich docstring, following the same format as that of previously implemented distributions.
It generally looks something like this:

```python
r"""Univariate blah distribution.

 The pdf of this distribution is

 .. math::

    f(x \mid \param1, \param2) = \exp{x * (param1 + \log{param2})}

 .. plot::

     import matplotlib.pyplot as plt
     import numpy as np
     import scipy.stats as st
     import arviz as az
     x = np.linspace(-5, 5, 1000)
     params1 = [0., 0., 0., -2.]
     params2 = [0.4, 1., 2., 0.4]
     for param1, param2 in zip(params1, params2):
         pdf = st.blah.pdf(x, param1, param2)
         plt.plot(x, pdf, label=r'$\param1$ = {}, $\param2$ = {}'.format(param1, param2))
     plt.xlabel('x', fontsize=12)
     plt.ylabel('f(x)', fontsize=12)
     plt.legend(loc=1)
     plt.show()

 ========  ==========================================
 Support   :math:`x \in [0, \infty)`
 ========  ==========================================

 Blah distribution can be parameterized either in terms of param2 or
 alt_param2. The link between the two parametrizations is
 given by

 .. math::

    \param2 = \dfrac{1}{\alt_param2}


 Parameters
 ----------
 param1: float
     Interpretation of param1.
 param2: float
     Interpretation of param2 (param2 > 0).
 alt_param2: float
     Interpretation of alt_param2 (alt_param2 > 0) (alternative to param2).

 Examples
 --------
 .. code-block:: python

     with pm.Model():
         x = pm.Blah('x', param1=0, param2=10)
 """
```

The new distribution should be referenced in the respective API page in the `docs` module (e.g., `pymc/docs/api/distributions.continuous.rst`).
If appropriate, a new notebook example should be added to [pymc-examples](https://github.com/pymc-devs/pymc-examples/) illustrating how this distribution can be used and how it relates (and/or differs) from other distributions that users are more likely to be familiar with.
# Pull request step-by-step
The preferred workflow for contributing to PyMC is to fork the [GitHub repository](https://github.com/pymc-devs/pymc/), clone it to your local machine, and develop on a feature branch.

## Steps

1. Fork the [project repository](https://github.com/pymc-devs/pymc/) by clicking on the 'Fork' button near the top right of the main repository page. This creates a copy of the code under your GitHub user account.

2. Clone your fork of the PyMC repo from your GitHub account to your local disk, and add the base repository as a remote:

   ```bash
   $ git clone git@github.com:<your GitHub handle>/pymc.git
   $ cd pymc
   $ git remote add upstream git@github.com:pymc-devs/pymc.git
   ```

3. Create a ``feature`` branch to hold your development changes:

   ```bash
   $ git checkout -b my-feature
   ```


   :::{attention}
   Always use a ``feature`` branch. It's good practice to never routinely work on the ``main`` branch of any repository.
   :::

4. Project requirements are in ``requirements.txt``, and libraries used for development are in ``requirements-dev.txt``. The easiest (and recommended) way to set up a development environment is via [miniconda](https://docs.conda.io/en/latest/miniconda.html):

   ```bash
   $ conda env create -f conda-envs/environment-dev-py37.yml  # or py38 or py39
   $ conda activate pymc-dev-py37
   $ pip install -e .
   ```

   _Alternatively_ you may (probably in a [virtual environment](https://docs.python-guide.org/dev/virtualenvs/)) run:

   ```bash
   $ pip install -e .
   $ pip install -r requirements-dev.txt
   ```

   Yet another alternative is to create a docker environment for development. See: [Developing in Docker](#Developing-in-Docker).

5. Develop the feature on your feature branch. Add changed files using ``git add`` and then ``git commit`` files:

   ```bash
   $ git add modified_files
   $ git commit
   ```

   to record your changes locally.
   After committing, it is a good idea to sync with the base repository in case there have been any changes:
   ```bash
   $ git fetch upstream
   $ git rebase upstream/main
   ```

   Then push the changes to your GitHub account with:

   ```bash
   $ git push -u origin my-feature
   ```

6. Go to the GitHub web page of your fork of the PyMC repo. Click the 'Pull request' button to send your changes to the project's maintainers for review. This will send an email to the committers.

   :::{tip}
   Now that your PR is ready, read the {ref}`pr_checklist` to make sure it follows best practices.
   :::
(citing_pymc)=
# Citing PyMC

TODO: extend advise, copy info from readme, explain citing paper plus specific zenodo version...

If you use PyMC in your reseach please cite: Salvatier J., Wiecki T.V., Fonnesbeck C. (2016) Probabilistic programming in Python using PyMC. PeerJ Computer Science 2:e55 [DOI: 10.7717/peerj-cs.55](https://doi.org/10.7717/peerj-cs.55).

The BibTeX entry is:

```bibtex
@article{pymc,
  title={Probabilistic programming in Python using PyMC3},
  author={Salvatier, John and Wiecki, Thomas V and Fonnesbeck, Christopher},
  journal={PeerJ Computer Science},
  volume={2},
  pages={e55},
  year={2016},
  publisher={PeerJ Inc.}
}
```
(index)=
# About PyMC


## Purpose

PyMC is a probabilistic programming package for Python that allows users to fit Bayesian models using a variety of numerical methods, most notably Markov chain Monte Carlo (MCMC) and variational inference (VI). Its flexibility and extensibility make it applicable to a large suite of problems. Along with core model specification and fitting functionality, PyMC includes functionality for summarizing output and for model diagnostics.

## Features

PyMC strives to make Bayesian modeling as simple and painless as possible,  allowing users to focus on their scientific problem, rather than on the methods used to solve it. Here is a partial list of its features:

* Modern methods for fitting Bayesian models, including MCMC and VI.

* Includes a large suite of well-documented statistical distributions.

* Uses Aesara as the computational backend, allowing for fast expression evaluation, automatic gradient calculation, and GPU computing.

* Built-in support for Gaussian process modeling.

* Model summarization and plotting.

* Model checking and convergence detection.

* Extensible: easily incorporates custom step methods and unusual probability
  distributions.

* Bayesian models can be embedded in larger programs, and results can be analyzed
  with the full power of Python.

:::{toctree}
:hidden:

citing_pymc
history_and_versions
pymc_for_enterprise
testimonials
:::
(pymc_for_enterprise)=

# PyMC for enterprise

`PyMC is now available as part of the Tidelift Subscription!`

Tidelift is working with PyMC and the maintainers of thousands of other open source
projects to deliver commercial support and maintenance for the open source dependencies
you use to build your applications. Save time, reduce risk, and improve code health,
while contributing financially to PyMC -- making it even more robust, reliable and,
let's face it, amazing!


TODO fix code

.. raw:: html

    <style>.centered {text-align: center;}</style>
    <p><div class="centered">
    <a href="https://tidelift.com/subscription/pkg/pypi-pymc?utm_source=undefined&utm_medium=referral&utm_campaign=enterprise">
      <button class="ui large orange button" color="orange">Learn more</button>
    </a>
    <a href="https://tidelift.com/subscription/request-a-demo?utm_source=undefined&utm_medium=referral&utm_campaign=enterprise">
      <button class="ui large orange button">Request a demo</button>
    </a>
    </div></p>

## Enterprise-ready open source software — managed for you

The Tidelift Subscription is a managed open source subscription for application
dependencies covering millions of open source projects across JavaScript, Python, Java,
PHP, Ruby, .NET, and more. And now, your favorite probabilistic programming language is included in the Tidelift subscription!

Your subscription includes:

* **Security updates**: Tidelift’s security response team coordinates patches for new breaking security vulnerabilities and alerts immediately through a private channel, so your software supply chain is always secure.

* **Licensing verification and indemnification**: Tidelift verifies license information to enable easy policy enforcement and adds intellectual property indemnification to cover creators and users in case something goes wrong. You always have a 100% up-to-date bill of materials for your dependencies to share with your legal team, customers, or partners.

* **Maintenance and code improvement**: Tidelift ensures the software you rely on keeps working as long as you need it to work. Your managed dependencies are actively maintained and Tidelift recruits additional maintainers where required.

* **Package selection and version guidance**: Tidelift helps you choose the best open source packages from the start—and then guides you through updates to stay on the best releases as new issues arise.

* **Roadmap input**: Take a seat at the table with the creators behind the software you use. PyMC developers and other Tidelift’s participating maintainers earn more income as our software is used by more subscribers, so we’re interested in knowing what you need.

* **Tooling and cloud integration**: Tidelift works with GitHub, GitLab, BitBucket, and more. It supports every cloud platform (and other deployment targets, too).

The end result? All of the capabilities you expect from commercial-grade software, for the full breadth of open source you use. That means less time grappling with esoteric open source trivia, and more time building your own applications — and your business.


TODO fix code
.. raw:: html

    <style>.centered {text-align: center;}</style>
    <p><div class="centered">
    <a href="https://tidelift.com/subscription/pkg/pypi-pymc3?utm_source=undefined&utm_medium=referral&utm_campaign=enterprise">
      <button class="ui large orange button" color="orange">Learn more</button>
    </a>
    <a href="https://tidelift.com/subscription/request-a-demo?utm_source=undefined&utm_medium=referral&utm_campaign=enterprise">
      <button class="ui large orange button">Request a demo</button>
    </a>
    </div></p>
(history_and_versions)=
# History

PyMC began development in 2003, as an effort to generalize the process of
building Metropolis-Hastings samplers, with an aim to making Markov chain Monte
Carlo (MCMC) more accessible to applied scientists.
The choice to develop PyMC as a python module, rather than a standalone
application, allowed the use of MCMC methods in a larger modeling framework. By
2005, PyMC was reliable enough for version 1.0 to be released to the public. A
small group of regular users, most associated with the University of Georgia,
provided much of the feedback necessary for the refinement of PyMC to a usable
state.

In 2006, David Huard and Anand Patil joined Chris Fonnesbeck on the development
team for PyMC 2.0. This iteration of the software strives for more flexibility,
better performance and a better end-user experience than any previous version
of PyMC. PyMC 2.2 was released in April 2012. It contained numerous bugfixes and
optimizations, as well as a few new features, including improved output
plotting, csv table output, improved imputation syntax, and posterior
predictive check plots. PyMC 2.3 was released on October 31, 2013. It included
Python 3 compatibility, improved summary plots, and some important bug fixes.

In 2011, John Salvatier began thinking about implementing gradient-based MCMC samplers, and developed the ``mcex`` package to experiment with his ideas. The following year, John was invited by the team to re-engineer PyMC to accomodate Hamiltonian Monte Carlo sampling. This led to the adoption of Theano as the computational back end, and marked the beginning of PyMC's development. The first alpha version of PyMC was released in June 2015. Over the following 2 years, the core development team grew to 12 members, and the first release, PyMC 3.0, was launched in January 2017.  In 2020 the PyMC developers forked Theano and in 2021 renamed the forked project to Aesara.


# What's new in version 4

:bdg-warning:`TODO`
Add text

# Version 3

The third major version of PyMC benefitted from being re-written from scratch. Substantial improvements in the user interface and performance resulted from this. While PyMC2 relied on Fortran extensions (via f2py) for most of the computational heavy-lifting, PyMC leverages Aesara, a fork of the Theano library from the Montréal Institute for Learning Algorithms (MILA), for array-based expression evaluation, to perform its computation. What this provided, above all else, is fast automatic differentiation, which is at the heart of the gradient-based sampling and optimization methods providing inference for probabilistic programming.

Major changes from previous versions:

* New flexible object model and syntax (not backward-compatible with PyMC2).

* Gradient-based MCMC methods, including Hamiltonian Monte Carlo (HMC), the No U-turn Sampler (NUTS), and Stein Variational Gradient Descent.

* Variational inference methods, including automatic differentiation variational inference (ADVI) and operator variational inference (OPVI).

* An interface for easy formula-based specification of generalized linear models (GLM).

* Elliptical slice sampling.

* Specialized distributions for representing time series.

* A library of Jupyter notebooks that provide case studies and fully developed usage examples.

* Much more!

While the addition of Aesara added a level of complexity to the development of PyMC, fundamentally altering how the underlying computation is performed, the dev team worked hard to maintain the elegant simplicity of the original PyMC model specification syntax.
:::{card}
:class-body: small
:link: https://www.spacex.com/

**SpaceX**
^^^
> At SpaceX PyMC helped us estimate supplier delivery uncertainty and quantify which suppliers were consistent in sending us rocket parts and which suppliers we needed to partner with to understand how we could reduce variability

_Ravin Kumar_
:::

:::{card}
:class-body: small
:link: https://www.salesforce.com

**Salesforce**
^^^
> PyMC is my primary tool for statistical modeling at Salesforce. I use it to combine disparate sources of information and pretty much anywhere that quantifying uncertainty is important. We've also been experimenting with gaussian processes to model time series data for forecasting.

_Eddie Landesberg. Manager, Data Scientist_
:::

:::{card}
:class-body: small
**[Novartis Institutes for Biomedical Research](https://www.novartis.com/our-science/novartis-institutes-biomedical-research)**
^^^
> At the Novartis Institutes for Biomedical Research, we use PyMC for a wide variety of scientific and business use cases. The API is incredible, making it easy to express probabilistic models of our scientific experimental data and business processes, such as models of electrophysiology and molecular dose response.

_Eric J. Ma_
:::

:::{card}
:class-body: small
**[Intercom](https://www.intercom.com)**
^^^
> At Intercom, we've adopted PyMC as part of our A/B testing framework. The API made it easy to integrate into our existing experimentation framework and the methodology has made communication of experiment results much more straightforward for non technical stakeholders.

_Louis Ryan_
:::

:::{card}
:class-body: small
**[Channel 4](http://www.channel4.co.uk)**
^^^
> Used in research code at Channel 4 for developing internal forecasting tools.

_Peader Coyle_
:::
(testimonials)=
# Testimonials

TODO: move testimonials from wiki to here
.. image:: https://cdn.rawgit.com/pymc-devs/pymc/main/docs/logos/svg/PyMC_banner.svg
    :height: 100px
    :alt: PyMC logo
    :align: center

|Build Status| |Coverage| |NumFOCUS_badge| |Binder| |Dockerhub| |DOIzenodo|

PyMC (formerly PyMC3) is a Python package for Bayesian statistical modeling
focusing on advanced Markov chain Monte Carlo (MCMC) and variational inference (VI)
algorithms. Its flexibility and extensibility make it applicable to a
large suite of problems.

Check out the `PyMC overview <https://docs.pymc.io/en/latest/learn/examples/pymc_overview.html>`__,  or
`interact with live examples <https://mybinder.org/v2/gh/pymc-devs/pymc/main?filepath=%2Fdocs%2Fsource%2Fnotebooks>`__
using Binder!
For questions on PyMC, head on over to our `PyMC Discourse <https://discourse.pymc.io/>`__ forum.

Features
========

-  Intuitive model specification syntax, for example, ``x ~ N(0,1)``
   translates to ``x = Normal('x',0,1)``
-  **Powerful sampling algorithms**, such as the `No U-Turn
   Sampler <http://www.jmlr.org/papers/v15/hoffman14a.html>`__, allow complex models
   with thousands of parameters with little specialized knowledge of
   fitting algorithms.
-  **Variational inference**: `ADVI <http://www.jmlr.org/papers/v18/16-107.html>`__
   for fast approximate posterior estimation as well as mini-batch ADVI
   for large data sets.
-  Relies on `Aesara <https://aesara.readthedocs.io/en/latest/>`__ which provides:
    *  Computation optimization and dynamic C or JAX compilation
    *  NumPy broadcasting and advanced indexing
    *  Linear algebra operators
    *  Simple extensibility
-  Transparent support for missing value imputation

Getting started
===============

If you already know about Bayesian statistics:
----------------------------------------------

-  `API quickstart guide <https://docs.pymc.io/en/stable/pymc-examples/examples/pymc3_howto/api_quickstart.html>`__
-  The `PyMC tutorial <https://docs.pymc.io/en/latest/learn/examples/pymc_overview.html>`__
-  `PyMC examples <https://docs.pymc.io/projects/examples/en/latest/>`__ and the `API reference <https://docs.pymc.io/en/stable/api.html>`__

Learn Bayesian statistics with a book together with PyMC
--------------------------------------------------------

-  `Probabilistic Programming and Bayesian Methods for Hackers <https://github.com/CamDavidsonPilon/Probabilistic-Programming-and-Bayesian-Methods-for-Hackers>`__: Fantastic book with many applied code examples.
-  `PyMC port of the book "Doing Bayesian Data Analysis" by John Kruschke <https://github.com/aloctavodia/Doing_bayesian_data_analysis>`__ as well as the `second edition <https://github.com/JWarmenhoven/DBDA-python>`__: Principled introduction to Bayesian data analysis.
-  `PyMC port of the book "Statistical Rethinking A Bayesian Course with Examples in R and Stan" by Richard McElreath <https://github.com/pymc-devs/resources/tree/master/Rethinking>`__
-  `PyMC port of the book "Bayesian Cognitive Modeling" by Michael Lee and EJ Wagenmakers <https://github.com/pymc-devs/resources/tree/master/BCM>`__: Focused on using Bayesian statistics in cognitive modeling.
-  `Bayesian Analysis with Python  <https://www.packtpub.com/big-data-and-business-intelligence/bayesian-analysis-python-second-edition>`__ (second edition) by Osvaldo Martin: Great introductory book. (`code <https://github.com/aloctavodia/BAP>`__ and errata).

Audio & Video
-------------

- Here is a `YouTube playlist <https://www.youtube.com/playlist?list=PL1Ma_1DBbE82OVW8Fz_6Ts1oOeyOAiovy>`__ gathering several talks on PyMC.
- You can also find all the talks given at **PyMCon 2020** `here <https://discourse.pymc.io/c/pymcon/2020talks/15>`__.
- The `"Learning Bayesian Statistics" podcast <https://www.learnbayesstats.com/>`__ helps you discover and stay up-to-date with the vast Bayesian community. Bonus: it's hosted by Alex Andorra, one of the PyMC core devs!

Installation
============

To install PyMC on your system, follow the instructions on the appropriate installation guide:

-  `Installing PyMC on MacOS <https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(MacOS)>`__
-  `Installing PyMC on Linux <https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(Linux)>`__
-  `Installing PyMC on Windows <https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(Windows)>`__


Citing PyMC
===========
Please choose from the following:

- |DOIpaper| *Probabilistic programming in Python using PyMC3*, Salvatier J., Wiecki T.V., Fonnesbeck C. (2016)
- |DOIzenodo| A DOI for all versions.
- DOIs for specific versions are shown on Zenodo and under `Releases <https://github.com/pymc-devs/pymc/releases>`_

.. |DOIpaper| image:: https://img.shields.io/badge/DOI-10.7717%2Fpeerj--cs.55-blue
     :target: https://doi.org/10.7717/peerj-cs.55
.. |DOIzenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4603970.svg
   :target: https://doi.org/10.5281/zenodo.4603970

Contact
=======

We are using `discourse.pymc.io <https://discourse.pymc.io/>`__ as our main communication channel. You can also follow us on `Twitter @pymc_devs <https://twitter.com/pymc_devs>`__ for updates and other announcements.

To ask a question regarding modeling or usage of PyMC we encourage posting to our Discourse forum under the `“Questions” Category <https://discourse.pymc.io/c/questions>`__. You can also suggest feature in the `“Development” Category <https://discourse.pymc.io/c/development>`__.

To report an issue with PyMC please use the `issue tracker <https://github.com/pymc-devs/pymc/issues>`__.

Finally, if you need to get in touch for non-technical information about the project, `send us an e-mail <pymc.devs@gmail.com>`__.

License
=======

`Apache License, Version
2.0 <https://github.com/pymc-devs/pymc/blob/main/LICENSE>`__


Software using PyMC
===================

General purpose
---------------

- `Bambi <https://github.com/bambinos/bambi>`__: BAyesian Model-Building Interface (BAMBI) in Python.
- `SunODE <https://github.com/aseyboldt/sunode>`__: Fast ODE solver, much faster than the one that comes with PyMC.
- `pymc-learn <https://github.com/pymc-learn/pymc-learn>`__: Custom PyMC models built on top of pymc3_models/scikit-learn API
- `fenics-pymc3 <https://github.com/IvanYashchuk/fenics-pymc3>`__: Differentiable interface to FEniCS, a library for solving partial differential equations.

Domain specific
---------------

- `Exoplanet <https://github.com/dfm/exoplanet>`__: a toolkit for modeling of transit and/or radial velocity observations of exoplanets and other astronomical time series.
- `NiPyMC <https://github.com/PsychoinformaticsLab/nipymc>`__: Bayesian mixed-effects modeling of fMRI data in Python.
- `beat <https://github.com/hvasbath/beat>`__: Bayesian Earthquake Analysis Tool.
- `cell2location <https://github.com/BayraktarLab/cell2location>`__: Comprehensive mapping of tissue cell architecture via integrated single cell and spatial transcriptomics.

Please contact us if your software is not listed here.

Papers citing PyMC
==================

See `Google Scholar <https://scholar.google.de/scholar?oi=bibs&hl=en&authuser=1&cites=6936955228135731011>`__ for a continuously updated list.

Contributors
============

See the `GitHub contributor
page <https://github.com/pymc-devs/pymc/graphs/contributors>`__. Also read our `Code of Conduct <https://github.com/pymc-devs/pymc/blob/main/CODE_OF_CONDUCT.md>`__ guidelines for a better contributing experience.

Support
=======

PyMC is a non-profit project under NumFOCUS umbrella. If you want to support PyMC financially, you can donate `here <https://numfocus.salsalabs.org/donate-to-pymc3/index.html>`__.

PyMC for enterprise
===================
`PyMC is now available as part of the Tidelift Subscription!`

Tidelift is working with PyMC and the maintainers of thousands of other open source
projects to deliver commercial support and maintenance for the open source dependencies
you use to build your applications. Save time, reduce risk, and improve code health,
while contributing financially to PyMC -- making it even more robust, reliable and,
let's face it, amazing!

|tidelift_learn| |tidelift_demo|

You can also get professional consulting support from `PyMC Labs <https://www.pymc-labs.io>`__.

Sponsors
========

|NumFOCUS|

|PyMCLabs|

.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/pymc-devs/pymc/main?filepath=%2Fdocs%2Fsource%2Fnotebooks
.. |Build Status| image:: https://github.com/pymc-devs/pymc/workflows/pytest/badge.svg
   :target: https://github.com/pymc-devs/pymc/actions
.. |Coverage| image:: https://codecov.io/gh/pymc-devs/pymc/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/pymc-devs/pymc
.. |Dockerhub| image:: https://img.shields.io/docker/automated/pymc/pymc.svg
   :target: https://hub.docker.com/r/pymc/pymc
.. |NumFOCUS| image:: https://www.numfocus.org/wp-content/uploads/2017/03/1457562110.png
   :target: http://www.numfocus.org/
.. |NumFOCUS_badge| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :target: http://www.numfocus.org/
.. |PyMCLabs| image:: https://raw.githubusercontent.com/pymc-devs/pymc/main/docs/logos/sponsors/pymc-labs.png
   :target: https://pymc-labs.io
.. |tidelift| image:: https://img.shields.io/badge/-lifted!-2dd160.svg?colorA=58595b&style=flat&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABEAAAAOCAYAAADJ7fe0AAAAAXNSR0IArs4c6QAAAAlwSFlzAAAWJQAAFiUBSVIk8AAAAVlpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IlhNUCBDb3JlIDUuNC4wIj4KICAgPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4KICAgICAgPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIKICAgICAgICAgICAgeG1sbnM6dGlmZj0iaHR0cDovL25zLmFkb2JlLmNvbS90aWZmLzEuMC8iPgogICAgICAgICA8dGlmZjpPcmllbnRhdGlvbj4xPC90aWZmOk9yaWVudGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KTMInWQAAAVhJREFUKBV1kj0vBFEUhmd2sdZHh2IlGhKFQuOviEYiNlFodCqtUqPxA%2FwCjUTnDygkGoVERFQaZFlE9nreO%2BdM5u5wkifvuee892Pu3CyEcA0DeIc%2B9IwftJsR6Cko3uCjguZdjuBZhhwmYDjGrOC96WED41UtsgEdGEAPlmAfpuAbFF%2BFZLfoMfRBGzThDtLgePPwBIpdddGzOArhPHUXowbNptE2www6a%2Fm96Y3pHN7oQ1s%2B13pxt1ENaKzBFWyWzaJ%2BRO0C9Jny6VPSoKjLVbMDC5bn5OPuJF%2BBSe95PVEMuugY5AegS9fCh7BedP45hRnj8TC34QQUe9bTZyh2KgvFk2vc8GIlXyTfsvqr6bPpNgv52ynnlomZJNpB70Xhl%2Bf6Sa02p1bApEfnETwxVa%2Faj%2BW%2FFtHltmxS%2FO3krvpTtTnVgu%2F6gvHRFvG78Ef3kOe5PimJXycY74blT5R%2BAAAAAElFTkSuQmCC
   :target: https://tidelift.com/subscription/pkg/pypi-pymc3?utm_source=pypi-pymc3&utm_medium=referral&utm_campaign=enterprise
.. |tidelift_learn| image:: https://img.shields.io/badge/-learn%20more-2dd160.svg?color=orange&labelColor=58595b&style=flat&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABEAAAAOCAYAAADJ7fe0AAAAAXNSR0IArs4c6QAAAAlwSFlzAAAWJQAAFiUBSVIk8AAAAVlpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IlhNUCBDb3JlIDUuNC4wIj4KICAgPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4KICAgICAgPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIKICAgICAgICAgICAgeG1sbnM6dGlmZj0iaHR0cDovL25zLmFkb2JlLmNvbS90aWZmLzEuMC8iPgogICAgICAgICA8dGlmZjpPcmllbnRhdGlvbj4xPC90aWZmOk9yaWVudGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KTMInWQAAAVhJREFUKBV1kj0vBFEUhmd2sdZHh2IlGhKFQuOviEYiNlFodCqtUqPxA%2FwCjUTnDygkGoVERFQaZFlE9nreO%2BdM5u5wkifvuee892Pu3CyEcA0DeIc%2B9IwftJsR6Cko3uCjguZdjuBZhhwmYDjGrOC96WED41UtsgEdGEAPlmAfpuAbFF%2BFZLfoMfRBGzThDtLgePPwBIpdddGzOArhPHUXowbNptE2www6a%2Fm96Y3pHN7oQ1s%2B13pxt1ENaKzBFWyWzaJ%2BRO0C9Jny6VPSoKjLVbMDC5bn5OPuJF%2BBSe95PVEMuugY5AegS9fCh7BedP45hRnj8TC34QQUe9bTZyh2KgvFk2vc8GIlXyTfsvqr6bPpNgv52ynnlomZJNpB70Xhl%2Bf6Sa02p1bApEfnETwxVa%2Faj%2BW%2FFtHltmxS%2FO3krvpTtTnVgu%2F6gvHRFvG78Ef3kOe5PimJXycY74blT5R%2BAAAAAElFTkSuQmCC
   :target: https://tidelift.com/subscription/pkg/pypi-pymc3?utm_source=pypi-pymc3&utm_medium=referral&utm_campaign=enterprise
.. |tidelift_demo| image:: https://img.shields.io/badge/-request%20a%20demo-2dd160.svg?color=orange&labelColor=58595b&style=flat&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABEAAAAOCAYAAADJ7fe0AAAAAXNSR0IArs4c6QAAAAlwSFlzAAAWJQAAFiUBSVIk8AAAAVlpVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IlhNUCBDb3JlIDUuNC4wIj4KICAgPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4KICAgICAgPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIKICAgICAgICAgICAgeG1sbnM6dGlmZj0iaHR0cDovL25zLmFkb2JlLmNvbS90aWZmLzEuMC8iPgogICAgICAgICA8dGlmZjpPcmllbnRhdGlvbj4xPC90aWZmOk9yaWVudGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KTMInWQAAAVhJREFUKBV1kj0vBFEUhmd2sdZHh2IlGhKFQuOviEYiNlFodCqtUqPxA%2FwCjUTnDygkGoVERFQaZFlE9nreO%2BdM5u5wkifvuee892Pu3CyEcA0DeIc%2B9IwftJsR6Cko3uCjguZdjuBZhhwmYDjGrOC96WED41UtsgEdGEAPlmAfpuAbFF%2BFZLfoMfRBGzThDtLgePPwBIpdddGzOArhPHUXowbNptE2www6a%2Fm96Y3pHN7oQ1s%2B13pxt1ENaKzBFWyWzaJ%2BRO0C9Jny6VPSoKjLVbMDC5bn5OPuJF%2BBSe95PVEMuugY5AegS9fCh7BedP45hRnj8TC34QQUe9bTZyh2KgvFk2vc8GIlXyTfsvqr6bPpNgv52ynnlomZJNpB70Xhl%2Bf6Sa02p1bApEfnETwxVa%2Faj%2BW%2FFtHltmxS%2FO3krvpTtTnVgu%2F6gvHRFvG78Ef3kOe5PimJXycY74blT5R%2BAAAAAElFTkSuQmCC
   :target: https://tidelift.com/subscription/request-a-demo?utm_source=pypi-pymc3&utm_medium=referral&utm_campaign=enterprise
:orphan:

..
    _href from docs/source/index.rst

===============
PyMC and Aesara
===============

What is Aesara
==============

Aesara is a package that allows us to define functions involving array
operations and linear algebra. When we define a PyMC model, we implicitly
build up an Aesara function from the space of our parameters to
their posterior probability density up to a constant factor. We then use
symbolic manipulations of this function to also get access to its gradient.

For a thorough introduction to Aesara see the
`aesara docs <https://aesara.readthedocs.io/en/latest/>`_,
but for the most part you don't need detailed knowledge about it as long
as you are not trying to define new distributions or other extensions
of PyMC. But let's look at a simple example to get a rough
idea about how it works. Say, we'd like to define the (completely
arbitrarily chosen) function

.. math::

  f\colon \mathbb{R} \times \mathbb{R}^n \times \mathbb{N}^n \to \mathbb{R}\\
  (a, x, y) \mapsto \sum_{i=0}^{n} \exp(ax_i^3 + y_i^2).


First, we need to define symbolic variables for our inputs (this
is similar to eg SymPy's `Symbol`)::

    import aesara
    import aesara.tensor as at
    # We don't specify the dtype of our input variables, so it
    # defaults to using float64 without any special config.
    a = at.scalar('a')
    x = at.vector('x')
    # `at.ivector` creates a symbolic vector of integers.
    y = at.ivector('y')

Next, we use those variables to build up a symbolic representation
of the output of our function. Note that no computation is actually
being done at this point. We only record what operations we need to
do to compute the output::

    inner = a * x**3 + y**2
    out = at.exp(inner).sum()

.. note::

   In this example we use `at.exp` to create a symbolic representation
   of the exponential of `inner`. Somewhat surprisingly, it
   would also have worked if we used `np.exp`. This is because numpy
   gives objects it operates on a chance to define the results of
   operations themselves. Aesara variables do this for a large number
   of operations. We usually still prefer the Aesara
   functions instead of the numpy versions, as that makes it clear that
   we are working with symbolic input instead of plain arrays.

Now we can tell Aesara to build a function that does this computation.
With a typical configuration, Aesara generates C code, compiles it,
and creates a python function which wraps the C function::

    func = aesara.function([a, x, y], [out])

We can call this function with actual arrays as many times as we want::

    a_val = 1.2
    x_vals = np.random.randn(10)
    y_vals = np.random.randn(10)

    out = func(a_val, x_vals, y_vals)

For the most part the symbolic Aesara variables can be operated on
like NumPy arrays. Most NumPy functions are available in `aesara.tensor`
(which is typically imported as `at`). A lot of linear algebra operations
can be found in `at.nlinalg` and `at.slinalg` (the NumPy and SciPy
operations respectively). Some support for sparse matrices is available
in `aesara.sparse`. For a detailed overview of available operations,
see `the aesara api docs <https://aesara.readthedocs.io/en/latest/library/tensor/index.html>`_.

A notable exception where Aesara variables do *not* behave like
NumPy arrays are operations involving conditional execution.

Code like this won't work as expected::

    a = at.vector('a')
    if (a > 0).all():
        b = at.sqrt(a)
    else:
        b = -a

`(a > 0).all()` isn't actually a boolean as it would be in NumPy, but
still a symbolic variable. Python will convert this object to a boolean
and according to the rules for this conversion, things that aren't empty
containers or zero are converted to `True`. So the code is equivalent
to this::

    a = at.vector('a')
    b = at.sqrt(a)

To get the desired behaviour, we can use `at.switch`::

    a = at.vector('a')
    b = at.switch((a > 0).all(), at.sqrt(a), -a)

Indexing also works similarly to NumPy::

    a = at.vector('a')
    # Access the 10th element. This will fail when a function build
    # from this expression is executed with an array that is too short.
    b = a[10]

    # Extract a subvector
    b = a[[1, 2, 10]]

Changing elements of an array is possible using `at.set_subtensor`::

    a = at.vector('a')
    b = at.set_subtensor(a[:10], 1)

    # is roughly equivalent to this (although aesara avoids
    # the copy if `a` isn't used anymore)
    a = np.random.randn(10)
    b = a.copy()
    b[:10] = 1

How PyMC uses Aesara
====================

Now that we have a basic understanding of Aesara we can look at what
happens if we define a PyMC model. Let's look at a simple example::

    true_mu = 0.1
    data = true_mu + np.random.randn(50)

    with pm.Model() as model:
        mu = pm.Normal('mu', mu=0, sigma=1)
        y = pm.Normal('y', mu=mu, sigma=1, observed=data)

In this model we define two variables: `mu` and `y`. The first is
a free variable that we want to infer, the second is an observed
variable. To sample from the posterior we need to build the function

.. math::

   \log P(μ|y) + C = \log P(y|μ) + \log P(μ) =: \text{logp}(μ)\\

where with the normal likelihood :math:`N(x|μ,σ^2)`

.. math::

    \text{logp}\colon \mathbb{R} \to \mathbb{R}\\
    μ \mapsto \log N(μ|0, 1) + \log N(y|μ, 1),

To build that function we need to keep track of two things: The parameter
space (the *free variables*) and the logp function. For each free variable
we generate an Aesara variable. And for each variable (observed or otherwise)
we add a term to the global logp. In the background something similar to
this is happening::

    # For illustration only, those functions don't actually exist
    # in exactly this way!
    model = pm.Model()

    mu = at.scalar('mu')
    model.add_free_variable(mu)
    model.add_logp_term(pm.Normal.dist(0, 1).logp(mu))

    model.add_logp_term(pm.Normal.dist(mu, 1).logp(data))

So calling `pm.Normal()` modifies the model: It changes the logp function
of the model. If the `observed` keyword isn't set it also creates a new
free variable. In contrast, `pm.Normal.dist()` doesn't care about the model,
it just creates an object that represents the normal distribution. Calling
`logp` on this object creates an Aesara variable for the logp probability
or log probability density of the distribution, but again without changing
the model in any way.

Continuous variables with support only on a subset of the real numbers
are treated a bit differently. We create a transformed variable
that has support on the reals and then modify this variable. For
example::

    with pm.Model() as model:
        mu = pm.Normal('mu', 0, 1)
        sd = pm.HalfNormal('sd', 1)
        y = pm.Normal('y', mu=mu, sigma=sd, observed=data)

is roughly equivalent to this::

    # For illustration only, not real code!
    model = pm.Model()
    mu = at.scalar('mu')
    model.add_free_variable(mu)
    model.add_logp_term(pm.Normal.dist(0, 1).logp(mu))

    sd_log__ = at.scalar('sd_log__')
    model.add_free_variable(sd_log__)
    model.add_logp_term(corrected_logp_half_normal(sd_log__))

    sd = at.exp(sd_log__)
    model.add_deterministic_variable(sd)

    model.add_logp_term(pm.Normal.dist(mu, sd).logp(data))

The return values of the variable constructors are subclasses
of Aesara variables, so when we define a variable we can use any
Aesara operation on them::

    design_matrix = np.array([[...]])
    with pm.Model() as model:
        # beta is a at.dvector
        beta = pm.Normal('beta', 0, 1, shape=len(design_matrix))
        predict = at.dot(design_matrix, beta)
        sd = pm.HalfCauchy('sd', beta=2.5)
        pm.Normal('y', mu=predict, sigma=sd, observed=data)
:orphan:
(Advanced_usage_of_Aesara_in_PyMC)=
..
    _referenced in docs/source/notebooks/table_of_contents_tutorials.js

=================================
Advanced usage of Aesara in PyMC
=================================

Using shared variables
======================

Shared variables allow us to use values in Aesara functions that are
not considered an input to the function, but can still be changed
later. They are very similar to global variables in may ways::

    a = at.scalar('a')
    # Create a new shared variable with initial value of 0.1
    b = aesara.shared(0.1)
    func = aesara.function([a], a * b)
    assert func(2.) == 0.2

    b.set_value(10.)
    assert func(2.) == 20.

Shared variables can also contain arrays, and are allowed to change
their shape as long as the number of dimensions stays the same.

We can use shared variables in PyMC to fit the same model to several
datasets without the need to recreate the model each time (which can
be time consuming if the number of datasets is large)::

    # We generate 10 datasets
    true_mu = [np.random.randn() for _ in range(10)]
    observed_data = [mu + np.random.randn(20) for mu in true_mu]

    data = aesara.shared(observed_data[0])
    with pm.Model() as model:
        mu = pm.Normal('mu', 0, 10)
        pm.Normal('y', mu=mu, sigma=1, observed=data)

    # Generate one trace for each dataset
    idatas = []
    for data_vals in observed_data:
        # Switch out the observed dataset
        data.set_value(data_vals)
        with model:
            idatas.append(pm.sample())

We can also sometimes use shared variables to work around limitations
in the current PyMC api. A common task in Machine Learning is to predict
values for unseen data, and one way to achieve this is to use a shared
variable for our observations::

    x = np.random.randn(100)
    y = x > 0

    x_shared = aesara.shared(x)

    with pm.Model() as model:
      coeff = pm.Normal('x', mu=0, sigma=1)
      logistic = pm.math.sigmoid(coeff * x_shared)
      pm.Bernoulli('obs', p=logistic, observed=y)

      # fit the model
      idata = pm.sample()

      # Switch out the observations and use `sample_posterior_predictive` to predict
      x_shared.set_value([-1, 0, 1.])
      post_pred = pm.sample_posterior_predictive(trace, samples=500)

However, due to the way we handle shapes at the moment, it is
not possible to change the shape of a shared variable if that would
also change the shape of one of the variables.


Writing custom Aesara Ops
=========================

While Aesara includes a wide range of operations, there are cases where
it makes sense to write your own. But before doing this it is a good
idea to think hard if it is actually necessary. Especially if you want
to use algorithms that need gradient information — this includes NUTS and
all variational methods, and you probably *should* want to use those —
this is often quite a bit of work and also requires some math and
debugging skills for the gradients.

Good reasons for defining a custom `Op` might be the following:

- You require an operation that is not available in Aesara and can't
  be build up out of existing Aesara operations. This could for example
  include models where you need to solve differential equations or
  integrals, or find a root or minimum of a function that depends
  on your parameters.
- You want to connect your PyMC model to some existing external code.
- After carefully considering different parametrizations and a lot
  of profiling your model is still too slow, but you know of a faster
  way to compute the gradient than what Aesara is doing. This faster
  way might be anything from clever maths to using more hardware.
  There is nothing stopping anyone from using a cluster via MPI in
  a custom node, if a part of the gradient computation is slow enough
  and sufficiently parallelizable to make the cost worth it.
  We would definitely like to hear about any such examples.

Aesara has extensive `documentation, <https://aesara.readthedocs.io/en/latest/extending/index.html>`_
about how to write new Ops.


Finding the root of a function
------------------------------

We'll use finding the root of a function as a simple example.
Let's say we want to define a model where a parameter is defined
implicitly as the root of a function, that depends on another
parameter:

.. math::

   \theta \sim N^+(0, 1)\\
   \text{$\mu\in \mathbb{R}^+$ such that $R(\mu, \theta)
         = \mu + \mu e^{\theta \mu} - 1= 0$}\\
   y \sim N(\mu, 0.1^2)

First, we observe that this problem is well-defined, because
:math:`R(\cdot, \theta)` is monotone and has the image :math:`(-1, \infty)`
for :math:`\mu, \theta \in \mathbb{R}^+`. To avoid overflows in
:math:`\exp(\mu \theta)` for large
values of :math:`\mu\theta` we instead find the root of

.. math::

    R'(\mu, \theta)
        = \log(R(\mu, \theta) + 1)
        = \log(\mu) + \log(1 + e^{\theta\mu}).

Also, we have

.. math::

    \frac{\partial}{\partial\mu}R'(\mu, \theta)
        = \theta\, \text{logit}^{-1}(\theta\mu) + \mu^{-1}.

We can now use `scipy.optimize.newton` to find the root::

    from scipy import optimize, special
    import numpy as np

    def func(mu, theta):
        thetamu = theta * mu
        value = np.log(mu) + np.logaddexp(0, thetamu)
        return value

    def jac(mu, theta):
        thetamu = theta * mu
        jac = theta * special.expit(thetamu) + 1 / mu
        return jac

    def mu_from_theta(theta):
        return optimize.newton(func, 1, fprime=jac, args=(theta,))

We could wrap `mu_from_theta` with `aesara.compile.ops.as_op` and use gradient-free
methods like Metropolis, but to get NUTS and ADVI working, we also
need to define the derivative of `mu_from_theta`. We can find this
derivative using the implicit function theorem, or equivalently we
take the derivative with respect of :math:`\theta` for both sides of
:math:`R(\mu(\theta), \theta) = 0` and solve for :math:`\frac{d\mu}{d\theta}`.
This isn't hard to do by hand, but for the fun of it, let's do it using
sympy::

    import sympy

    mu = sympy.Function('mu')
    theta = sympy.Symbol('theta')
    R = mu(theta) + mu(theta) * sympy.exp(theta * mu(theta)) - 1
    solution = sympy.solve(R.diff(theta), mu(theta).diff(theta))[0]

We get

.. math::

    \frac{d}{d\theta}\mu(\theta)
        = - \frac{\mu(\theta)^2}{1 + \theta\mu(\theta) + e^{-\theta\mu(\theta)}}

Now, we use this to define a Aesara `Op`, that also computes the gradient::

    import aesara
    import aesara.tensor as at
    import aesara.tests.unittest_tools
    from aesara.graph.op import Op

    class MuFromTheta(Op):
        itypes = [at.dscalar]
        otypes = [at.dscalar]

        def perform(self, node, inputs, outputs):
            theta, = inputs
            mu = mu_from_theta(theta)
            outputs[0][0] = np.array(mu)

        def grad(self, inputs, g):
            theta, = inputs
            mu = self(theta)
            thetamu = theta * mu
            return [- g[0] * mu ** 2 / (1 + thetamu + at.exp(-thetamu))]

If you value your sanity, always check that the gradient is ok::

    aesara.gradient.verify_grad(MuFromTheta(), [np.array(0.2)])
    aesara.gradient.verify_grad(MuFromTheta(), [np.array(1e-5)])
    aesara.gradient.verify_grad(MuFromTheta(), [np.array(1e5)])

We can now define our model using this new `Op`::

    import pymc as pm

    at_mu_from_theta = MuFromTheta()

    with pm.Model() as model:
        theta = pm.HalfNormal('theta', sigma=1)
        mu = pm.Deterministic('mu', at_mu_from_theta(theta))
        pm.Normal('y', mu=mu, sigma=0.1, observed=[0.2, 0.21, 0.3])

        idata = pm.sample()
.. _api:

*************
API Reference
*************

.. toctree::
   :maxdepth: 2

   api/distributions
   api/gp
   api/model
   api/ode
   api/samplers
   api/smc
   api/step_methods
   api/inference

--------------
API extensions
--------------

Plots, stats and diagnostics
----------------------------
Plots, stats and diagnostics are delegated to the
:doc:`ArviZ <arviz:index>`.
library, a general purpose library for
"exploratory analysis of Bayesian models".

* Functions from the `arviz.plots` module are available through ``pymc.<function>`` or ``pymc.plots.<function>``,
but for their API documentation please refer to the :ref:`ArviZ documentation <arviz:plot_api>`.

* Functions from the `arviz.stats` module are available through ``pymc.<function>`` or ``pymc.stats.<function>``,
but for their API documentation please refer to the :ref:`ArviZ documentation <arviz:stats_api>`.

ArviZ is a dependency of PyMC and so, in addition to the locations described above,
importing ArviZ and using ``arviz.<function>`` will also work without any extra installation.

Generalized Linear Models (GLMs)
--------------------------------

Generalized Linear Models are delegated to the
`Bambi <https://bambinos.github.io/bambi>`_.
library, a high-level Bayesian model-building
interface built on top of PyMC.

Bambi is not a dependency of PyMC and should be installed in addition to PyMC
to use it to generate PyMC models via formula syntax.
:orphan:

..
    _href from docs/source/index.rst

.. _prob_dists:

*********************************
Probability Distributions in PyMC
*********************************

The most fundamental step in building Bayesian models is the specification of a full probability model for the problem at hand. This primarily involves assigning parametric statistical distributions to unknown quantities in the model, in addition to appropriate functional forms for likelihoods to represent the information from the data. To this end, PyMC includes a comprehensive set of pre-defined statistical distributions that can be used as model building blocks.

For example, if we wish to define a particular variable as having a normal prior, we can specify that using an instance of the ``Normal`` class.

::

    with pm.Model():

        x = pm.Normal('x', mu=0, sigma=1)

A variable requires at least a ``name`` argument, and zero or more model parameters, depending on the distribution. Parameter names vary by distribution, using conventional names wherever possible. The example above defines a scalar variable. To make a vector-valued variable, a ``shape`` argument should be provided; for example, a 3x3 matrix of beta random variables could be defined with:

::

    with pm.Model():

        p = pm.Beta('p', 1, 1, shape=(3, 3))

Probability distributions are all subclasses of ``Distribution``, which in turn has two major subclasses: ``Discrete`` and ``Continuous``. In terms of data types, a ``Continuous`` random variable is given whichever floating point type is defined by ``aesara.config.floatX``, while ``Discrete`` variables are given ``int16`` types when ``aesara.config.floatX`` is ``float32``, and ``int64`` otherwise.

All distributions in ``pm.distributions`` will have two important methods: ``random()`` and ``logp()`` with the following signatures:

::

    class SomeDistribution(Continuous):

        def random(self, point=None, size=None):
            ...
            return random_samples

        def logp(self, value):
            ...
            return total_log_prob

PyMC expects the ``logp()`` method to return a log-probability evaluated at the passed ``value`` argument. This method is used internally by all of the inference methods to calculate the model log-probability that is used for fitting models. The ``random()`` method is used to simulate values from the variable, and is used internally for posterior predictive checks.


Custom distributions
====================

Despite the fact that PyMC ships with a large set of the most common probability distributions, some problems may require the use of functional forms that are less common, and not available in ``pm.distributions``. One example of this is in survival analysis, where time-to-event data is modeled using probability densities that are designed to accommodate censored data.

An exponential survival function, where :math:`c=0` denotes failure (or non-survival), is defined by:

.. math::

    f(c, t) = \left\{ \begin{array}{l} \exp(-\lambda t), \text{if c=1} \\
               \lambda \exp(-\lambda t), \text{if c=0}  \end{array} \right.

Such a function can be implemented as a PyMC distribution by writing a function that specifies the log-probability, then passing that function as a keyword argument to the ``DensityDist`` function, which creates an instance of a PyMC distribution with the custom function as its log-probability.

For the exponential survival function, this is:

::

    def logp(value, t, lam):
        return (value * log(lam) - lam * t).sum()

    exp_surv = pm.DensityDist('exp_surv', t, lam, logp=logp, observed=failure)

Similarly, if a random number generator is required, a function returning random numbers corresponding to the probability distribution can be passed as the ``random`` argument.


Using PyMC distributions without a Model
========================================

Distribution objects, as we have defined them so far, are only usable inside of a ``Model`` context. If they are created outside of the model context manager, it raises an error.

::

    y = Binomial('y', n=10, p=0.5)


::

    TypeError: No context on context stack

This is because the distribution classes are designed to integrate themselves automatically inside of a PyMC model. When a model cannot be found, it fails. However, each ``Distribution`` has a ``dist`` class method that returns a stripped-down distribution object that can be used outside of a PyMC model.

For example, a standalone binomial distribution can be created by:

::

    y = pm.Binomial.dist(n=10, p=0.5)

This allows for probabilities to be calculated and random numbers to be drawn.

::

    >>> y.logp(4).eval()
    array(-1.5843639373779297, dtype=float32)

    >>> y.random(size=3)
    array([5, 4, 3])


Auto-transformation
===================

To aid efficient MCMC sampling, any continuous variables that are constrained to a sub-interval of the real line are automatically transformed so that their support is unconstrained. This frees sampling algorithms from having to deal with boundary constraints.

For example, the gamma distribution is positive-valued. If we define one for a model:

::

    with pm.Model() as model:
        g = pm.Gamma('g', 1, 1)

We notice a modified variable inside the model ``value_vars`` attribute.  These variables represent the values of each random variable in the model's log-likelihood.

::

    >>> model.value_vars
    [g_log__]

As the name suggests, the variable ``g`` has been log-transformed, and this is the space over which posterior sampling takes place.

The value of the transformed variable is simply back-transformed when a sample is drawn in order to recover the original variable.

By default, auto-transformed variables are ignored when summarizing and plotting model output.
:orphan:

..
    _href from docs/source/index.rst

******************
Gaussian Processes
******************

GP Basics
=========

Sometimes an unknown parameter or variable in a model is not a scalar value or
a fixed-length vector, but a *function*.  A Gaussian process (GP) can be used
as a prior probability distribution whose support is over the space of
continuous functions.  A GP prior on the function :math:`f(x)` is usually written,

.. math::

  f(x) \sim \mathcal{GP}(m(x), \, k(x, x')) \,.

The function values are modeled as a draw from a multivariate normal
distribution that is parameterized by the mean function, :math:`m(x)`, and the
covariance function, :math:`k(x, x')`.  Gaussian processes are a convenient
choice as priors over functions due to the marginalization and conditioning
properties of the multivariate normal distribution.  Usually, the marginal
distribution over :math:`f(x)` is evaluated during the inference step.  The
conditional distribution is then used for predicting the function values
:math:`f(x_*)` at new points, :math:`x_*`.

The joint distribution of :math:`f(x)` and :math:`f(x_*)` is multivariate
normal,

.. math::

  \begin{bmatrix} f(x) \\ f(x_*) \\ \end{bmatrix} \sim
  \text{N}\left(
    \begin{bmatrix} m(x)  \\ m(x_*)    \\ \end{bmatrix} \,,
    \begin{bmatrix} k(x,x')    & k(x_*, x)    \\
                    k(x_*, x) &  k(x_*, x_*')  \\ \end{bmatrix}
          \right) \,.

Starting from the joint distribution, one obtains the marginal distribution
of :math:`f(x)`, as :math:`\text{N}(m(x),\, k(x, x'))`.  The conditional
distribution is

.. math::

  f(x_*) \mid f(x) \sim \text{N}\left( k(x_*, x) k(x, x)^{-1} [f(x) - m(x)] + m(x_*) ,\,
    k(x_*, x_*) - k(x, x_*) k(x, x)^{-1} k(x, x_*) \right) \,.

.. note::

  For more information on GPs, check out the book `Gaussian Processes for
  Machine Learning <http://www.gaussianprocess.org/gpml/>`_ by Rasmussen &
  Williams, or `this introduction <https://www.ics.uci.edu/~welling/teaching/KernelsICS273B/gpB.pdf>`_
  by D. Mackay.

PyMC is a great environment for working with fully Bayesian Gaussian Process
models.  GPs in PyMC have a clear syntax and are highly composable, and many
predefined covariance functions (or kernels), mean functions, and several GP
implementations are included.  GPs are treated as distributions that can be
used within larger or hierarchical models, not just as standalone regression
models.

Mean and covariance functions
=============================

Those who have used the GPy or GPflow Python packages will find the syntax for
construction mean and covariance functions somewhat familiar.  When first
instantiated, the mean and covariance functions are parameterized, but not
given their inputs yet.  The covariance functions must additionally be provided
with the number of dimensions of the input matrix, and a list that indexes
which of those dimensions they are to operate on.  The reason for this design
is so that covariance functions can be constructed that are combinations of
other covariance functions.

For example, to construct an exponentiated quadratic covariance function that
operates on the second and third column of a three column matrix representing
three predictor variables::

    ls = [2, 5] # the lengthscales
    cov_func = pm.gp.cov.ExpQuad(input_dim=3, ls=ls, active_dims=[1, 2])

Here the :code:`ls`, or lengthscale, parameter is two dimensional, allowing the second
and third dimension to have a different lengthscale.  The reason we have to
specify :code:`input_dim`, the total number of columns of :code:`X`, and
:code:`active_dims`, which of those columns or dimensions the covariance
function will act on, is because :code:`cov_func` hasn't actually seen the
input data yet.  The :code:`active_dims` argument is optional, and defaults to
all columns of the matrix of inputs.

Covariance functions in PyMC closely follow the algebraic rules for kernels,
which allows users to combine covariance functions into new ones, for example:

- The sum of two covariance functions is also a covariance function::


    cov_func = pm.gp.cov.ExpQuad(...) + pm.gp.cov.ExpQuad(...)

- The product of two covariance functions is also a covariance function::


    cov_func = pm.gp.cov.ExpQuad(...) * pm.gp.cov.Periodic(...)

- The product (or sum) of a covariance function with a scalar is a
  covariance function::


    cov_func = eta**2 * pm.gp.cov.Matern32(...)



After the covariance function is defined, it is now a function that is
evaluated by calling :code:`cov_func(x, x)` (or :code:`mean_func(x)`).  Since
PyMC is built on top of Aesara, it is relatively easy to define and experiment
with non-standard covariance and mean functons.  For more information check out
the tutorial on covariance functions.


GP Implementations
==================

PyMC includes several GP implementations, including marginal and latent
variable models and also some fast approximations.  Their usage all follows a
similar pattern:  First, a GP is instantiated with a mean function and a
covariance function.  Then, GP objects can be added together, allowing for
function characteristics to be carefully modeled and separated.  Finally, one
of `prior`, `marginal_likelihood` or `conditional` methods is called on the GP
object to actually construct the PyMC random variable that represents the
function prior.

Using :code:`gp.Latent` for the example, the syntax to first specify the GP
is::

    gp = pm.gp.Latent(mean_func, cov_func)

The first argument is the mean function and the second is the covariance
function.  We've made the GP object, but we haven't made clear which function
it is to be a prior for, what the inputs are, or what parameters it will be
conditioned on.

.. note::

  The :code:`gp.Marginal` class and similar don't have a :code:`prior` method.
  Instead they have a :code:`marginal_likelihood` method that is used similarly,
  but has additional required arguments, such as the observed data, noise,
  or other, depending on the implementation.  See the notebooks for examples.
  The :code:`conditional` method works similarly.

Calling the `prior` method will create a PyMC random variable that represents
the latent function :math:`f(x) = \mathbf{f}`::

	f = gp.prior("f", X)

:code:`f` is a random variable that can be used within a PyMC model like any
other type of random variable.  The first argument is the name of the random
variable representing the function we are placing the prior over.
The second argument is the inputs to the function that the prior is over,
:code:`X`.  The inputs are usually known and present in the data, but they can
also be PyMC random variables.  If the inputs are an Aesara tensor or a
PyMC random variable, the :code:`shape` needs to be given.

Usually at this point, inference is performed on the model.  The
:code:`conditional` method creates the conditional, or predictive,
distribution over the latent function at arbitrary :math:`x_*` input points,
:math:`f(x_*)`.  To construct the conditional distribution we write::

	f_star = gp.conditional("f_star", X_star)

Additive GPs
============

The GP implementation in PyMC is constructed so that it is easy to define
additive GPs and sample from individual GP components.  We can write::

    gp1 = pm.gp.Marginal(mean_func1, cov_func1)
    gp2 = pm.gp.Marginal(mean_func2, cov_func2)
    gp3 = gp1 + gp2

The GP objects have to have the same type, :code:`gp.Marginal` cannot
be added to :code:`gp.Latent`.

Consider two independent GP distributed functions, :math:`f_1(x) \sim
\mathcal{GP}\left(m_1(x),\, k_1(x, x')\right)` and :math:`f_2(x) \sim
\mathcal{GP}\left( m_2(x),\, k_2(x, x')\right)`.  The joint distribution of
:math:`f_1,\, f_1^*,\, f_2,\, f_2^*,\, f_1 + f_2 and f_1^* + f_2^*` is

.. math::

  \begin{bmatrix} f_1 \\ f_1^* \\ f_2 \\ f_2^*
               \\ f_1 + f_2    \\ f_1^* + f_2^* \end{bmatrix} \sim
  \text{N}\left(
    \begin{bmatrix} m_1 \\ m_1^* \\ m_2 \\ m_2^* \\
                    m_1 + m_2    \\ m_1^* + m_2^*   \\ \end{bmatrix} \,,\,
    \begin{bmatrix}
      K_1       &  K_1^*     &   0       &    0      & K_1        & K_1^*              \\
      K_1^{*^T} &  K_1^{**}  &   0       &    0      & K_1^*      & K_1^{**}           \\
      0         &  0         & K_2       & K_2^*     & K_2        & K_2^{*}            \\
      0         &  0         & K_2^{*^T} & K_2^{**}  & K_2^{*}    & K_2^{**}           \\
      K_1       &  K_1^{*}   & K_2       & K_2^{*}   & K_1 + K_2  & K_1^{*} + K_2^{*}  \\
      K_1^{*^T} & K_1^{**} & K_2^{*^T} & K_2^{**} & K_1^{*^T}+K_2^{*^T} & K_1^{**}+K_2^{**}
    \end{bmatrix}
  \right) \,.

Using the joint distribution to obtain the conditional distribution of :math:`f_1^*`
with the contribution due to :math:`f_1 + f_2` factored out, we get

.. math::
  f_1^* \mid f_1 + f_2 \sim \text{N}\left(
    m_1^* + K_1^{*^T}(K_1 + K_2)^{-1}\left[f_1 + f_2 - m_1 - m_2\right] \,,\,
    K_1^{**} - K_1^{*^T}(K_1 + K_2)^{-1}K_1^* \right) \,.


These equations show how to break down GP models into individual components to see how each
contributes to the data.  For more information, check out `David Duvenaud's PhD
thesis <https://www.cs.toronto.edu/~duvenaud/thesis.pdf>`_.

The GP objects in PyMC keeps track of these marginals automatically.  The
following code sketch shows how to define the conditional distribution of
:math:`f_2^*`.  We use `gp.Marginal` in the example, but the same works for
other implementations.  The first block fits the GP prior.  We denote
:math:`f_1 + f_2` as just :math:`f` for brevity::

    with pm.Model() as model:
        gp1 = pm.gp.Marginal(mean_func1, cov_func1)
        gp2 = pm.gp.Marginal(mean_func2, cov_func2)

        # gp represents f1 + f2.
        gp = gp1 + gp2

        f = gp.marginal_likelihood("f", X, y, noise)

        idata = pm.sample(1000)


To construct the conditional distribution of :code:`gp1` or :code:`gp2`, we
also need to include the additional arguments, :code:`X`, :code:`y`, and
:code:`noise`::

    with model:
        # conditional distributions of f1 and f2
        f1_star = gp1.conditional("f1_star", X_star,
                                  given={"X": X, "y": y, "noise": noise, "gp": gp})
        f2_star = gp2.conditional("f2_star", X_star,
                                  given={"X": X, "y": y, "noise": noise, "gp": gp})

        # conditional of f1 + f2, `given` not required
        f_star = gp.conditional("f_star", X_star)

This second block produces the conditional distributions.  Notice that extra
arguments are required for conditionals of :math:`f1` and :math:`f2`, but not
:math:`f`.  This is because those arguments are cached when
:code:`.marginal_likelihood` is called on :code:`gp`.

.. note::
  When constructing conditionals, the additional arguments :code:`X`, :code:`y`,
  :code:`noise` and :code:`gp` must be provided as a dict called `given`!

Since the marginal likelihoood method of :code:`gp1` or :code:`gp2` weren't called,
their conditionals need to be provided with the required inputs.  In the same
fashion as the prior, :code:`f_star`, :code:`f1_star` and :code:`f2_star` are random
variables that can now be used like any other random variable in PyMC.

Check the notebooks for detailed demonstrations of the usage of GP functionality
in PyMC.
{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   {% if methods %}

   .. rubric:: Methods

   .. autosummary::
      :toctree: classmethods

   {% for item in methods %}
      {{ objname }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
TODO: incorporate the useful bits of this page into the learning section

**************
Usage Overview
**************

For a detailed overview of building models in PyMC, please read the appropriate sections in the rest of the documentation. For a flavor of what PyMC models look like, here is a quick example.

First, let's import PyMC and `ArviZ <https://arviz-devs.github.io/arviz/>`__ (which handles plotting and diagnostics):

::

    import arviz as az
    import numpy as np
    import pymc as pm

Models are defined using a context manager (``with`` statement). The model is specified declaratively inside the context manager, instantiating model variables and transforming them as necessary. Here is an example of a model for a bioassay experiment:

::

    # Set style
    az.style.use("arviz-darkgrid")

    # Data
    n = np.ones(4)*5
    y = np.array([0, 1, 3, 5])
    dose = np.array([-.86,-.3,-.05,.73])

    with pm.Model() as bioassay_model:

        # Prior distributions for latent variables
        alpha = pm.Normal('alpha', 0, sigma=10)
        beta = pm.Normal('beta', 0, sigma=1)

        # Linear combination of parameters
        theta = pm.invlogit(alpha + beta * dose)

        # Model likelihood
        deaths = pm.Binomial('deaths', n=n, p=theta, observed=y)

Save this file, then from a python shell (or another file in the same directory), call:

::

    with bioassay_model:

        # Draw samples
        idata = pm.sample(1000, tune=2000, cores=2)
        # Plot two parameters
        az.plot_forest(idata, var_names=['alpha', 'beta'], r_hat=True)

This example will generate 1000 posterior samples on each of two cores using the NUTS algorithm, preceded by 2000 tuning samples (these are good default numbers for most models).

::

    Auto-assigning NUTS sampler...
    Initializing NUTS using jitter+adapt_diag...
    Multiprocess sampling (2 chains in 2 jobs)
    NUTS: [beta, alpha]
    |██████████████████████████████████████| 100.00% [6000/6000 00:04<00:00 Sampling 2 chains, 0 divergences]

The sample is returned as arrays inside a ``MultiTrace`` object, which is then passed to the plotting function. The resulting graph shows a forest plot of the random variables in the model, along with a convergence diagnostic (R-hat) that indicates our model has converged.

.. image:: ./images/forestplot.png
   :width: 1000px

See also
========

* `Tutorials <nb_tutorials/index.html>`__
* `Examples <nb_examples/index.html>`__


.. |NumFOCUS| image:: https://numfocus.org/wp-content/uploads/2017/07/NumFocus_LRG.png
   :target: http://www.numfocus.org/
   :height: 120px
.. |PyMCLabs| image:: https://raw.githubusercontent.com/pymc-devs/pymc/main/docs/pymc-labs-logo.png
   :target: https://pymc-labs.io
   :height: 120px
:orphan:

..
    _referenced in html_theme_options docs/source/conf.py

====================
PyMC Developer Guide
====================

`PyMC <https://docs.pymc.io/>`__ is a Python package for Bayesian
statistical modeling built on top of
`Aesara <https://aesara.readthedocs.io/en/latest/index.html>`__. This
document aims to explain the design and implementation of probabilistic
programming in PyMC, with comparisons to other PPL like TensorFlow Probability (TFP)
and Pyro in mind. A user-facing API
introduction can be found in the `API
quickstart <https://docs.pymc.io/notebooks/api_quickstart.html>`__. A more accessible, user facing deep introduction can be found in
`Peadar Coyle's probabilistic programming primer <https://github.com/springcoil/probabilisticprogrammingprimer>`__

Distribution
------------

A high-level introduction of ``Distribution`` in PyMC can be found in
the `documentation <https://docs.pymc.io/Probability_Distributions.html>`__. The source
code of the probability distributions is nested under
`pymc/distributions <https://github.com/pymc-devs/pymc/blob/main/pymc/distributions/>`__,
with the ``Distribution`` class defined in `distribution.py
<https://github.com/pymc-devs/pymc/blob/main/pymc/distributions/distribution.py#L23-L44>`__.
A few important points to highlight in the Distribution Class:

.. code:: python

    class Distribution:
        """Statistical distribution"""
        def __new__(cls, name, *args, **kwargs):
            ...
            try:
                model = Model.get_context()
            except TypeError:
                raise TypeError(...

            if isinstance(name, string_types):
                ...
                dist = cls.dist(*args, **kwargs)
                return model.Var(name, dist, ...)
            ...

In a way, the snippet above represents the unique features of pymc's
``Distribution`` class:

- Distribution objects are only usable inside of a ``Model`` context. If they are created outside of the model context manager, it raises an error.

- A ``Distribution`` requires at least a name argument, and other parameters that defines the Distribution.

- When a ``Distribution`` is initialized inside of a Model context, two things happen:

  1. a stateless distribution is initialized ``dist = {DISTRIBUTION_cls}.dist(*args, **kwargs)``;
  2. a random variable following the said distribution is added to the model ``model.Var(name, dist, ...)``

Thus, users who are building models using ``with pm.Model() ...`` should
be aware that they are never directly exposed to static and stateless
distributions, but rather random variables that follow some density
functions. Instead, to access a stateless distribution, you need to call
``pm.SomeDistribution.dist(...)`` or ``RV.dist`` *after* you initialized
``RV`` in a model context (see
https://docs.pymc.io/Probability\_Distributions.html#using-pymc-distributions-without-a-model).

With this distinction in mind, we can take a closer look at the
stateless distribution part of pymc (see distribution api in `doc
<https://docs.pymc.io/api/distributions.html>`__), which divided into:

- Continuous

- Discrete

- Multivariate

- Mixture

- Timeseries

Quote from the doc:

    All distributions in ``pm.distributions`` will have two important
    methods: ``random()`` and ``logp()`` with the following signatures:

.. code:: python

    class SomeDistribution(Continuous):
        def __init__(...):
            ...

        def random(self, point=None, size=None):
            ...
            return random_samples

        def logp(self, value):
            ...
            return total_log_prob

PyMC expects the ``logp()`` method to return a log-probability
evaluated at the passed value argument. This method is used internally
by all of the inference methods to calculate the model log-probability,
which is then used for fitting models. The ``random()`` method is
used to simulate values from the variable, and is used internally for
posterior predictive checks.

In the PyMC ``Distribution`` class, the ``logp()`` method is the most
elementary. As long as you have a well-behaved density function, we can
use it in the model to build the model log-likelihood function. Random
number generation is great to have, but sometimes there might not be
efficient random number generator for some densities. Since a function
is all you need, you can wrap almost any Aesara function into a
distribution using ``pm.DensityDist``
https://docs.pymc.io/Probability\_Distributions.html#custom-distributions

Thus, distributions that are defined in the ``distributions`` submodule
(e.g. look at ``pm.Normal`` in ``pymc.distributions.continuous``), each
describes a *family* of probabilistic distribution (no different from
distribution in other PPL library). Once it is initialised within a
model context, it contains properties that are related to the random
variable (*e.g.* mean/expectation). Note that if the parameters are
constants, these properties could be the same as the distribution
properties.

Reflection
~~~~~~~~~~

How tensor/value semantics for probability distributions is enabled in pymc:

In PyMC, we treat ``x = Normal('x', 0, 1)`` as defining a random
variable (intercepted and collected under a model context, more on that
below), and x.dist() as the associated density/mass function
(distribution in the mathematical sense). It is not perfect, and now
after a few years learning Bayesian statistics I also realized these
subtleties (i.e., the distinction between *random variable* and
*distribution*).

But when I was learning probabilistic modelling as a
beginner, I did find this approach to be the easiest and most
straightforward. In a perfect world, we should have
:math:`x \sim \text{Normal}(0, 1)` which defines a random variable that
follows a Gaussian distribution, and
:math:`\chi = \text{Normal}(0, 1), x \sim \chi` which define a `probability
density function <https://en.wikipedia.org/wiki/Probability_density_function>`__ that takes input :math:`x`

.. math::
    X:=f(x) = \frac{1}{\sigma \sqrt{2 \pi}} \exp^{- 0.5 (\frac{x - \mu}{\sigma})^2}\vert_{\mu = 0, \sigma=1} = \frac{1}{\sqrt{2 \pi}} \exp^{- 0.5 x^2}

Within a model context, RVs are essentially Aesara tensors (more on that
below). This is different than TFP and pyro, where you need to be more
explicit about the conversion. For example:

**PyMC**

.. code:: python

    with pm.Model() as model:
        z = pm.Normal('z', mu=0., sigma=5.)             # ==> aesara.tensor.var.TensorVariable
        x = pm.Normal('x', mu=z, sigma=1., observed=5.) # ==> aesara.tensor.var.TensorVariable
    x.logp({'z': 2.5})                                  # ==> -4.0439386
    model.logp({'z': 2.5})                              # ==> -6.6973152

**TFP**

.. code:: python

    import tensorflow.compat.v1 as tf
    from tensorflow_probability import distributions as tfd

    with tf.Session() as sess:
        z_dist = tfd.Normal(loc=0., scale=5.)            # ==> <class 'tfp.python.distributions.normal.Normal'>
        z = z_dist.sample()                              # ==> <class 'tensorflow.python.framework.ops.Tensor'>
        x = tfd.Normal(loc=z, scale=1.).log_prob(5.)     # ==> <class 'tensorflow.python.framework.ops.Tensor'>
        model_logp = z_dist.log_prob(z) + x
        print(sess.run(x, feed_dict={z: 2.5}))           # ==> -4.0439386
        print(sess.run(model_logp, feed_dict={z: 2.5}))  # ==> -6.6973152

**pyro**

.. code:: python

    z_dist = dist.Normal(loc=0., scale=5.)           # ==> <class 'pyro.distributions.torch.Normal'>
    z = pyro.sample("z", z_dist)                     # ==> <class 'torch.Tensor'>
    # reset/specify value of z
    z.data = torch.tensor(2.5)
    x = dist.Normal(loc=z, scale=1.).log_prob(5.)    # ==> <class 'torch.Tensor'>
    model_logp = z_dist.log_prob(z) + x
    x                                                # ==> -4.0439386
    model_logp                                       # ==> -6.6973152


``logp`` method, very different behind the curtain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``logp`` method is straightforward - it is an Aesara function within each
distribution. It has the following signature:

.. code:: python

    def logp(self, value):
        # GET PARAMETERS
        param1, param2, ... = self.params1, self.params2, ...
        # EVALUATE LOG-LIKELIHOOD FUNCTION, all inputs are (or array that could be convert to) Aesara tensor
        total_log_prob = f(param1, param2, ..., value)
        return total_log_prob

In the ``logp`` method, parameters and values are either Aesara tensors,
or could be converted to tensors. It is rather convenient as the
evaluation of logp is represented as a tensor (``RV.logpt``), and when
we linked different ``logp`` together (e.g., summing all ``RVs.logpt``
to get the model totall logp) the dependence is taken care of by Aesara
when the graph is built and compiled. Again, since the compiled function
depends on the nodes that already in the graph, whenever you want to generate
a new function that takes new input tensors you either need to regenerate the graph
with the appropriate dependencies, or replace the node by editing the existing graph.
In PyMC we use the second approach by using ``aesara.clone_replace()`` when it is needed.

As explained above, distribution in a ``pm.Model()`` context
automatically turn into a tensor with distribution property (pymc
random variable). To get the logp of a free\_RV is just evaluating the
``logp()`` `on
itself <https://github.com/pymc-devs/pymc/blob/6d07591962a6c135640a3c31903eba66b34e71d8/pymc/model.py#L1212-L1213>`__:

.. code:: python

        # self is a aesara.tensor with a distribution attached
        self.logp_sum_unscaledt = distribution.logp_sum(self)
        self.logp_nojac_unscaledt = distribution.logp_nojac(self)

Or for an observed RV. it evaluate the logp on the data:

.. code:: python

        self.logp_sum_unscaledt = distribution.logp_sum(data)
        self.logp_nojac_unscaledt = distribution.logp_nojac(data)

Model context and Random Variable
---------------------------------

I like to think that the ``with pm.Model() ...`` is a key syntax feature
and *the* signature of PyMC model language, and in general a great
out-of-the-box thinking/usage of the context manager in Python (with
`some
critics <https://twitter.com/_szhang/status/890793373740617729>`__, of
course).

Essentially `what a context manager
does <https://www.python.org/dev/peps/pep-0343/>`__ is:

.. code:: python

    with EXPR as VAR:
        USERCODE

which roughly translates into this:

.. code:: python

    VAR = EXPR
    VAR.__enter__()
    try:
        USERCODE
    finally:
        VAR.__exit__()

or conceptually:

.. code:: python

    with EXPR as VAR:
        # DO SOMETHING
        USERCODE
        # DO SOME ADDITIONAL THINGS

So what happened within the ``with pm.Model() as model: ...`` block,
besides the initial set up ``model = pm.Model()``? Starting from the
most elementary:

Random Variable
~~~~~~~~~~~~~~~

From the above session, we know that when we call eg
``pm.Normal('x', ...)`` within a Model context, it returns a random
variable. Thus, we have two equivalent ways of adding random variable to
a model:


.. code:: python

    with pm.Model() as m:
        x = pm.Normal('x', mu=0., sigma=1.)


.. parsed-literal::

    print(type(x))                              # ==> <class 'aesara.tensor.var.TensorVariable'>
    print(m.free_RVs)                           # ==> [x]
    print(logpt(x, 5.0))                        # ==> Elemwise{switch,no_inplace}.0
    print(logpt(x, 5.).eval({}))                # ==> -13.418938533204672
    print(m.logp({'x': 5.}))                    # ==> -13.418938533204672


In general, if a variable has observations (``observed`` parameter), the RV is
an observed RV, otherwise if it has a ``transformed`` (``transform`` parameter)
attribute, it is a transformed RV otherwise, it will be the most elementary
form: a free RV.  Note that this means that random variables with observations
cannot be transformed.

..
   Below, I will take a deeper look into transformed RV. A normal user
   might not necessarily come in contact with the concept, since a
   transformed RV and ``TransformedDistribution`` are intentionally not
   user facing.

   Because in PyMC there is no bijector class like in TFP or pyro, we only
   have a partial implementation called ``Transform``, which implements
   Jacobian correction for forward mapping only (there is no Jacobian
   correction for inverse mapping). The use cases we considered are limited
   to the set of distributions that are bounded, and the transformation
   maps the bounded set to the real line - see
   `doc
   <https://docs.pymc.io/notebooks/api_quickstart.html#Automatic-transforms-of-bounded-RVs>`__.
   However, other transformations are possible.
   In general, PyMC does not provide explicit functionality to transform
   one distribution to another. Instead, a dedicated distribution is
   usually created in order to optimise performance. But getting a
   ``TransformedDistribution`` is also possible (see also in
   `doc <https://docs.pymc.io/notebooks/api_quickstart.html#Transformed-distributions-and-changes-of-variables>`__):

   .. code:: python


       lognorm = Exp().apply(pm.Normal.dist(0., 1.))
       lognorm


   .. parsed-literal::

       <pymc.distributions.transforms.TransformedDistribution at 0x7f1536749b00>



   Now, back to ``model.RV(...)`` - things returned from ``model.RV(...)``
   are Aesara tensor variables, and it is clear from looking at
   ``TransformedRV``:

   .. code:: python

       class TransformedRV(TensorVariable):
           ...

   as for ``FreeRV`` and ``ObservedRV``, they are ``TensorVariable``\s with
   ``Factor`` as mixin:

   .. code:: python

       class FreeRV(Factor, TensorVariable):
           ...

   ``Factor`` basically `enable and assign the
   logp <https://github.com/pymc-devs/pymc/blob/6d07591962a6c135640a3c31903eba66b34e71d8/pymc/model.py#L195-L276>`__
   (representated as a tensor also) property to an Aesara tensor (thus
   making it a random variable). For a ``TransformedRV``, it transforms the
   distribution into a ``TransformedDistribution``, and then ``model.Var`` is
   called again to added the RV associated with the
   ``TransformedDistribution`` as a ``FreeRV``:

   .. code:: python

           ...
           self.transformed = model.Var(
                       transformed_name, transform.apply(distribution), total_size=total_size)

   note: after ``transform.apply(distribution)`` its ``.transform``
   porperty is set to ``None``, thus making sure that the above call will
   only add one ``FreeRV``. In another word, you *cannot* do chain
   transformation by nested applying multiple transforms to a Distribution
   (however, you can use `Chain
   transformation <https://docs.pymc.io/notebooks/api_quickstart.html?highlight=chain%20transformation>`__).

   .. code:: python

       z = pm.LogNormal.dist(mu=0., sigma=1., transform=tr.Log)
       z.transform           # ==> pymc.distributions.transforms.Log


   .. code:: python

       z2 = Exp().apply(z)
       z2.transform is None  # ==> True



Additional things that ``pm.Model`` does
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a way, ``pm.Model`` is a tape machine that records what is being
added to the model, it keeps track the random variables (observed or
unobserved) and potential term (additional tensor that to be added to
the model logp), and also deterministic transformation (as bookkeeping):
named\_vars, free\_RVs, observed\_RVs, deterministics, potentials,
missing\_values. The model context then computes some simple model
properties, builds a bijection mapping that transforms between
dictionary and numpy/Aesara ndarray, thus allowing the ``logp``/``dlogp`` functions
to have two equivalent versions: one takes a ``dict`` as input and the other
takes an ``ndarray`` as input. More importantly, a ``pm.Model()`` contains methods
to compile Aesara functions that take Random Variables (that are also
initialised within the same model) as input, for example:

.. code:: python

    with pm.Model() as m:
        z = pm.Normal('z', 0., 10., shape=10)
        x = pm.Normal('x', z, 1., shape=10)

    print(m.initial_point)
    print(m.dict_to_array(m.initial_point))  # ==> m.bijection.map(m.initial_point)
    print(m.bijection.rmap(np.arange(20)))


.. parsed-literal::

    {'z': array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]), 'x': array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])}
    [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
    {'z': array([10., 11., 12., 13., 14., 15., 16., 17., 18., 19.]), 'x': array([0., 1., 2., 3., 4., 5., 6., 7., 8., 9.])}


.. code:: python

    list(filter(lambda x: "logp" in x, dir(pm.Model)))


.. parsed-literal::

    ['d2logp',
     'd2logp_nojac',
     'datalogpt',
     'dlogp',
     'dlogp_array',
     'dlogp_nojac',
     'fastd2logp',
     'fastd2logp_nojac',
     'fastdlogp',
     'fastdlogp_nojac',
     'fastlogp',
     'fastlogp_nojac',
     'logp',
     'logp_array',
     'logp_dlogp_function',
     'logp_elemwise',
     'logp_nojac',
     'logp_nojact',
     'logpt',
     'varlogpt']



Logp and dlogp
--------------

The model collects all the random variables (everything in
``model.free_RVs`` and ``model.observed_RVs``) and potential term, and
sum them together to get the model logp:

.. code:: python

    @property
    def logpt(self):
        """Aesara scalar of log-probability of the model"""
        with self:
            factors = [var.logpt for var in self.basic_RVs] + self.potentials
            logp = at.sum([at.sum(factor) for factor in factors])
            ...
            return logp

which returns an Aesara tensor that its value depends on the free
parameters in the model (i.e., its parent nodes from the Aesara
graph).You can evaluate or compile into a python callable (that you can
pass numpy as input args). Note that the logp tensor depends on its
input in the Aesara graph, thus you cannot pass new tensor to generate a
logp function. For similar reason, in PyMC we do graph copying a lot
using aesara.clone_replace to replace the inputs to a tensor.

.. code:: python

    with pm.Model() as m:
        z = pm.Normal('z', 0., 10., shape=10)
        x = pm.Normal('x', z, 1., shape=10)
        y = pm.Normal('y', x.sum(), 1., observed=2.5)

    print(m.basic_RVs)    # ==> [z, x, y]
    print(m.free_RVs)     # ==> [z, x]


.. code:: python

    type(m.logpt)         # ==> aesara.tensor.var.TensorVariable


.. code:: python

    m.logpt.eval({x: np.random.randn(*x.tag.test_value.shape) for x in m.free_RVs})

output:

.. parsed-literal::

    array(-51.25369126)



PyMC then compiles a logp function with gradient that takes
``model.free_RVs`` as input and ``model.logpt`` as output. It could be a
subset of tensors in ``model.free_RVs`` if we want a conditional
logp/dlogp function:

.. code:: python

    def logp_dlogp_function(self, grad_vars=None, **kwargs):
        if grad_vars is None:
            grad_vars = list(typefilter(self.free_RVs, continuous_types))
        else:
            ...
        varnames = [var.name for var in grad_vars]  # In a simple case with only continous RVs,
                                                    # this is all the free_RVs
        extra_vars = [var for var in self.free_RVs if var.name not in varnames]
        return ValueGradFunction(self.logpt, grad_vars, extra_vars, **kwargs)

``ValueGradFunction`` is a callable class which isolates part of the
Aesara graph to compile additional Aesara functions. PyMC relies on
``aesara.clone_replace`` to copy the ``model.logpt`` and replace its input. It
does not edit or rewrite the graph directly.

The important parts of the above function is highlighted and commented.
On a high level, it allows us to build conditional logp function and its
gradient easily. Here is a taste of how it works in action:

.. code:: python

    inputlist = [np.random.randn(*x.tag.test_value.shape) for x in m.free_RVs]

    func = m.logp_dlogp_function()
    func.set_extra_values({})
    input_dict = {x.name: y for x, y in zip(m.free_RVs, inputlist)}
    print(input_dict)
    input_array = func.dict_to_array(input_dict)
    print(input_array)
    print(" ===== ")
    func(input_array)


.. parsed-literal::

    {'z': array([-0.7202002 ,  0.58712205, -1.44120196, -0.53153001, -0.36028732,
           -1.49098414, -0.80046792, -0.26351819,  1.91841949,  1.60004128]), 'x': array([ 0.01490006,  0.60958275, -0.06955203, -0.42430833, -1.43392303,
            1.13713493,  0.31650495, -0.62582879,  0.75642811,  0.50114527])}
    [-0.7202002   0.58712205 -1.44120196 -0.53153001 -0.36028732 -1.49098414
     -0.80046792 -0.26351819  1.91841949  1.60004128  0.01490006  0.60958275
     -0.06955203 -0.42430833 -1.43392303  1.13713493  0.31650495 -0.62582879
      0.75642811  0.50114527]
     =====
    (array(-51.0769075),
     array([ 0.74230226,  0.01658948,  1.38606194,  0.11253699, -1.07003284,
             2.64302891,  1.12497754, -0.35967542, -1.18117557, -1.11489642,
             0.98281586,  1.69545542,  0.34626619,  1.61069443,  2.79155183,
            -0.91020295,  0.60094326,  2.08022672,  2.8799075 ,  2.81681213]))



.. code:: python

    irv = 1
    print("Condition Logp: take %s as input and conditioned on the rest."%(m.free_RVs[irv].name))
    func_conditional = m.logp_dlogp_function(grad_vars=[m.free_RVs[irv]])
    func_conditional.set_extra_values(input_dict)
    input_array2 = func_conditional.dict_to_array(input_dict)
    print(input_array2)
    print(" ===== ")
    func_conditional(input_array2)


.. parsed-literal::

    Condition Logp: take x as input and conditioned on the rest.
    [ 0.01490006  0.60958275 -0.06955203 -0.42430833 -1.43392303  1.13713493
      0.31650495 -0.62582879  0.75642811  0.50114527]
     =====
    (array(-51.0769075),
     array([ 0.98281586,  1.69545542,  0.34626619,  1.61069443,  2.79155183,
            -0.91020295,  0.60094326,  2.08022672,  2.8799075 ,  2.81681213]))



So why is this necessary? One can imagine that we just compile one logp
function, and do bookkeeping ourselves. For example, we can build the
logp function in Aesara directly:

.. code:: python

    import aesara
    func = aesara.function(m.free_RVs, m.logpt)
    func(*inputlist)


.. parsed-literal::

    array(-51.0769075)



.. code:: python

    logpt_grad = aesara.grad(m.logpt, m.free_RVs)
    func_d = aesara.function(m.free_RVs, logpt_grad)
    func_d(*inputlist)


.. parsed-literal::

    [array([ 0.74230226,  0.01658948,  1.38606194,  0.11253699, -1.07003284,
             2.64302891,  1.12497754, -0.35967542, -1.18117557, -1.11489642]),
     array([ 0.98281586,  1.69545542,  0.34626619,  1.61069443,  2.79155183,
            -0.91020295,  0.60094326,  2.08022672,  2.8799075 ,  2.81681213])]



Similarly, build a conditional logp:

.. code:: python

    shared = aesara.shared(inputlist[1])
    func2 = aesara.function([m.free_RVs[0]], m.logpt, givens=[(m.free_RVs[1], shared)])
    print(func2(inputlist[0]))

    logpt_grad2 = aesara.grad(m.logpt, m.free_RVs[0])
    func_d2 = aesara.function([m.free_RVs[0]], logpt_grad2, givens=[(m.free_RVs[1], shared)])
    print(func_d2(inputlist[0]))


.. parsed-literal::

    -51.07690750130328
    [ 0.74230226  0.01658948  1.38606194  0.11253699 -1.07003284  2.64302891
      1.12497754 -0.35967542 -1.18117557 -1.11489642]


The above also gives the same logp and gradient as the output from
``model.logp_dlogp_function``. But the difficulty is to compile
everything into a single function:

.. code:: python

    func_logp_and_grad = aesara.function(m.free_RVs, [m.logpt, logpt_grad])  # ==> ERROR


We want to have a function that return the evaluation and its gradient
re each input: ``value, grad = f(x)``, but the naive implementation does
not work. We can of course wrap 2 functions - one for logp one for dlogp
- and output a list. But that would mean we need to call 2 functions. In
addition, when we write code using python logic to do bookkeeping when
we build our conditional logp. Using ``aesara.clone_replace``, we always have
the input to the Aesara function being a 1d vector (instead of a list of
RV that each can have very different shape), thus it is very easy to do
matrix operation like rotation etc.

Notes
~~~~~

| The current setup is quite powerful, as the Aesara compiled function
  is fairly fast to compile and to call. Also, when we are repeatedly
  calling a conditional logp function, external RV only need to reset
  once. However, there are still significant overheads when we are
  passing values between Aesara graph and numpy. That is the reason we
  often see no advantage in using GPU, because the data is copying
  between GPU and CPU at each function call - and for a small model, the
  result is a slower inference under GPU than CPU.
| Also, ``aesara.clone_replace`` is too convenient (pymc internal joke is that
  it is like a drug - very addictive). If all the operation happens in
  the graph (including the conditioning and setting value), I see no
  need to isolate part of the graph (via graph copying or graph
  rewriting) for building model and running inference.
| Moreover, if we are limiting to the problem that we can solved most
  confidently - model with all continous unknown parameters that could
  be sampled with dynamic HMC, there is even less need to think about
  graph cloning/rewriting.

Inference
---------

MCMC
~~~~

The ability for model instance to generate conditional logp and dlogp
function enable one of the unique feature of PyMC - `CompoundStep
method <https://docs.pymc.io/notebooks/sampling_compound_step.html>`__.
On a conceptual level it is a Metropolis-within-Gibbs sampler. User can
`specify different sampler of different
RVs <https://docs.pymc.io/notebooks/sampling_compound_step.html?highlight=compoundstep#Specify-compound-steps>`__.
Alternatively, it is implemented as yet another interceptor: the
``pm.sample(...)`` call will try to `assign the best step methods to
different
free\_RVs <https://github.com/pymc-devs/pymc/blob/6d07591962a6c135640a3c31903eba66b34e71d8/pymc/sampling.py#L86-L152>`__
(e.g., NUTS if all free\_RVs are continous). Then, (conditional) logp
function(s) are compiled, and the sampler called each sampler within the
list of CompoundStep in a for-loop for one sample circle.

For each sampler, it implements a ``step.step`` method to perform MH
updates. Each time a dictionary (``point`` in ``PyMC`` land, same
structure as ``model.initial_point``) is passed as input and output a new
dictionary with the free\_RVs being sampled now has a new value (if
accepted, see
`here <https://github.com/pymc-devs/pymc/blob/6d07591962a6c135640a3c31903eba66b34e71d8/pymc/step_methods/compound.py#L27>`__
and
`here <https://github.com/pymc-devs/pymc/blob/main/pymc/step_methods/compound.py>`__).
There are some example in the `CompoundStep
doc <https://docs.pymc.io/notebooks/sampling_compound_step.html#Specify-compound-steps>`__.

Transition kernel
^^^^^^^^^^^^^^^^^

The base class for most MCMC sampler (except SMC) is in
`ArrayStep <https://github.com/pymc-devs/pymc/blob/main/pymc/step_methods/arraystep.py>`__.
You can see that the ``step.step()`` is mapping the ``point`` into an
array, and call ``self.astep()``, which is an array in, array out
function. A pymc model compile a conditional logp/dlogp function that
replace the input RVs with a shared 1D tensor (flatten and stack view of
the original RVs). And the transition kernel (i.e., ``.astep()``) takes
array as input and output an array. See for example in the `MH
sampler <https://github.com/pymc-devs/pymc/blob/6d07591962a6c135640a3c31903eba66b34e71d8/pymc/step_methods/metropolis.py#L139-L173>`__.

This is of course very different compare to the transition kernel in eg
TFP, which is a tenor in tensor out function. Moreover, transition
kernels in TFP do not flatten the tensors, see eg docstring of
`tensorflow\_probability/python/mcmc/random\_walk\_metropolis.py <https://github.com/tensorflow/probability/blob/master/tensorflow_probability/python/mcmc/random_walk_metropolis.py>`__:

.. code::

          new_state_fn: Python callable which takes a list of state parts and a
            seed; returns a same-type `list` of `Tensor`s, each being a perturbation
            of the input state parts. The perturbation distribution is assumed to be
            a symmetric distribution centered at the input state part.
            Default value: `None` which is mapped to
              `tfp.mcmc.random_walk_normal_fn()`.


Dynamic HMC
^^^^^^^^^^^

We love NUTS, or to be more precise Dynamic HMC with complex stopping
rules. This part is actually all done outside of Aesara, for NUTS, it
includes: the leapfrog, dual averaging, tunning of mass matrix and step
size, the tree building, sampler related statistics like divergence and
energy checking. We actually have an Aesara version of HMC, but it has never
been used, and has been removed from the main repository. It can still be
found in the `git history
<https://github.com/pymc-devs/pymc/pull/3734/commits/0fdae8207fd14f66635f3673ef267b2b8817aa68>`__,
though.

Variational Inference (VI)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The design of the VI module takes a different approach than
MCMC - it has a functional design, and everything is done within Aesara
(i.e., Optimization and building the variational objective). The base
class of variational inference is
`pymc.variational.Inference <https://github.com/pymc-devs/pymc/blob/main/pymc/variational/inference.py>`__,
where it builds the objective function by calling:

.. code:: python

        ...
        self.objective = op(approx, **kwargs)(tf)
        ...

Where:

.. code::

        op     : Operator class
        approx : Approximation class or instance
        tf     : TestFunction instance
        kwargs : kwargs passed to :class:`Operator`

The design is inspired by the great work `Operator Variational
Inference <https://arxiv.org/abs/1610.09033>`__. ``Inference`` object is
a very high level of VI implementation. It uses primitives: Operator,
Approximation, and Test functions to combine them into single objective
function. Currently we do not care too much about the test function, it
is usually not required (and not implemented). The other primitives are
defined as base classes in `this
file <https://github.com/pymc-devs/pymc/blob/main/pymc/variational/opvi.py>`__.
We use inheritance to easily implement a broad class of VI methods
leaving a lot of flexibility for further extensions.

For example, consider ADVI. We know that in the high-level, we are
approximating the posterior in the latent space with a diagonal
Multivariate Gaussian. In another word, we are approximating each elements in
``model.free_RVs`` with a Gaussian. Below is what happen in the set up:

.. code:: python

    def __init__(self, *args, **kwargs):
        super(ADVI, self).__init__(MeanField(*args, **kwargs))
    # ==> In the super class KLqp
        super(KLqp, self).__init__(KL, MeanField(*args, **kwargs), None, beta=beta)
    # ==> In the super class Inference
        ...
        self.objective = KL(MeanField(*args, **kwargs))(None)
        ...

where ``KL`` is Operator based on Kullback Leibler Divergence (it does
not need any test function).

.. code:: python

        ...
        def apply(self, f):
            return -self.datalogp_norm + self.beta * (self.logq_norm - self.varlogp_norm)

Since the logp and logq are from the approximation, let's dive in
further on it (there is another abstraction here - ``Group`` - that
allows you to combine approximation into new approximation, but we will
skip this for now and only consider ``SingleGroupApproximation`` like
``MeanField``): The definition of ``datalogp_norm``, ``logq_norm``,
``varlogp_norm`` are in
`variational/opvi <https://github.com/pymc-devs/pymc/blob/main/pymc/variational/opvi.py>`__,
strip away the normalizing term, ``datalogp`` and ``varlogp`` are
expectation of the variational free\_RVs and data logp - we clone the
datalogp and varlogp from the model, replace its input with Aesara
tensor that `samples from the variational
posterior <https://github.com/pymc-devs/pymc/blob/6d07591962a6c135640a3c31903eba66b34e71d8/pymc/variational/opvi.py#L1098-L1111>`__.
For ADVI, these samples are from `a
Gaussian <https://github.com/pymc-devs/pymc/blob/6d07591962a6c135640a3c31903eba66b34e71d8/pymc/variational/approximations.py#L84-L89>`__.
Note that the samples from the posterior approximations are usually 1
dimension more, so that we can compute the expectation and get the
gradient of the expectation (by computing the `expectation of the
gradient! <http://blog.shakirm.com/2015/10/machine-learning-trick-of-the-day-4-reparameterisation-tricks/>`__).
As for the ``logq`` since it is a Gaussian `it is pretty
straightforward to evaluate <https://github.com/pymc-devs/pymc/blob/6d07591962a6c135640a3c31903eba66b34e71d8/pymc/variational/approximations.py#L91-L97>`__.

Some challenges and insights from implementing VI.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Graph based approach was helpful, but Aesara had no direct access to
   previously created nodes in the computational graph. you can find a
   lot of ``@node_property`` usages in implementation. This is done to
   cache nodes. TensorFlow has graph utils for that that could
   potentially help in doing this. On the other hand graph management in
   Tensorflow seemed to more tricky than expected. The high level reason
   is that graph is an add only container

-  There were few fixed bugs not obvoius in the first place. Aesara has
   a tool to manipulate the graph (``aesara.clone_replace``) and this tool
   requires extremely careful treatment when doing a lot of graph
   replacements at different level.

-  We coined a term ``aesara.clone_replace`` curse. We got extremely dependent
   on this feature. Internal usages are uncountable:

   -  we use this to `vectorize the
      model <https://github.com/pymc-devs/pymc/blob/main/pymc/model.py#L972>`__
      for both MCMC and VI to speed up computations
   -  we use this to `create sampling
      graph <https://github.com/pymc-devs/pymc/blob/main/pymc/variational/opvi.py#L1483>`__
      for VI. This is the case you want posterior predictive as a part
      of computational graph.

As this is the core of the VI process, we were trying to replicate this pattern
in TF. However, when ``aesara.clone_replace`` is called, Aesara creates a new part of the graph that can
be collected by garbage collector, but TF's graph is add only. So we
should solve the problem of replacing input in a different way.

Forward sampling
----------------

As explained above, in distribution we have method to walk the model
dependence graph and generate forward random sample in scipy/numpy. This
allows us to do prior predictive samples using
``pymc.sampling.sample_prior_predictive`` see `code <https://github.com/pymc-devs/pymc/blob/6d07591962a6c135640a3c31903eba66b34e71d8/pymc/sampling.py#L1303-L1345>`__.
It is a fairly fast batch operation, but we have quite a lot of bugs and
edge case especially in high dimensions. The biggest pain point is the
automatic broadcasting. As in the batch random generation, we want to
generate (n\_sample, ) + RV.shape random samples. In some cases, where
we broadcast RV1 and RV2 to create a RV3 that has one more batch shape,
we get error (even worse, wrong answer with silent error).

The good news is, we are fixing these errors with the amazing works from [lucianopaz](https://github.com/lucianopaz) and
others. The challenge and some summary of the solution could be found in Luciano's [blog post](https://lucianopaz.github.io/2019/08/19/pymc-shape-handling/)

.. code:: python

    with pm.Model() as m:
        mu = pm.Normal('mu', 0., 1., shape=(5, 1))
        sd = pm.HalfNormal('sd', 5., shape=(1, 10))
        pm.Normal('x', mu=mu, sigma=sd, observed=np.random.randn(2, 5, 10))
        trace = pm.sample_prior_predictive(100)

    trace['x'].shape # ==> should be (100, 2, 5, 10)

.. code:: python

    pm.Normal.dist(mu=np.zeros(2), sigma=1).random(size=(10, 4))

There are also other error related random sample generation (e.g.,
`Mixture is currently
broken <https://github.com/pymc-devs/pymc/issues/3270>`__).

Extending PyMC
--------------

-  Custom Inference method
    -  `Inferencing Linear Mixed Model with EM.ipynb <https://github.com/junpenglao/Planet_Sakaar_Data_Science/blob/master/Ports/Inferencing%20Linear%20Mixed%20Model%20with%20EM.ipynb>`__
    -  `Laplace approximation in  pymc.ipynb <https://github.com/junpenglao/Planet_Sakaar_Data_Science/blob/master/Ports/Laplace%20approximation%20in%20pymc.ipynb>`__
-  Connecting it to other library within a model
    -  `Using “black box” likelihood function by creating a custom Aesara Op <https://docs.pymc.io/notebooks/blackbox_external_likelihood.html>`__
    -  Using emcee
-  Using other library for inference
    -  Connecting to Julia for solving ODE (with gradient for solution that can be used in NUTS)

What we got wrong
-----------------

Shape
~~~~~

One of the pain point we often face is the issue of shape. The approach
in TFP and pyro is currently much more rigorous. Adrian’s PR
(https://github.com/pymc-devs/pymc/pull/2833) might fix this problem,
but likely it is a huge effort of refactoring. I implemented quite a lot
of patches for mixture distribution, but still they are not done very
naturally.

Random methods in numpy
~~~~~~~~~~~~~~~~~~~~~~~

There is a lot of complex logic for sampling from random variables, and
because it is all in Python, we can't transform a sampling graph
further. Unfortunately, Aesara does not have code to sample from various
distributions and we didn't want to write that our own.

Samplers are in Python
~~~~~~~~~~~~~~~~~~~~~~

While having the samplers be written in Python allows for a lot of
flexibility and intuitive for experiment (writing e.g. NUTS in Aesara is
also very difficult), it comes at a performance penalty and makes
sampling on the GPU very inefficient because memory needs to be copied
for every logp evaluation.
Gaussian Processes
------------------

.. automodule:: pymc.gp

.. toctree::

   gp/implementations
   gp/mean
   gp/cov
***********
shape_utils
***********

This submodule contains various functions that apply numpy's broadcasting rules to shape tuples, and also to samples drawn from probability distributions.

The main challenge when broadcasting samples drawn from a generative model, is that each random variate has a core shape. When we draw many i.i.d samples from a given RV, for example if we ask for `size_tuple` i.i.d draws, the result usually is a `size_tuple + RV_core_shape`. In the generative model's hierarchy, the downstream RVs that are conditionally dependent on our above sampled values, will get an array with a shape that is incosistent with the core shape they expect to see for their parameters. This is a problem sometimes because it prevents regular broadcasting in complex hierachical models, and thus make prior and posterior predictive sampling difficult.

This module introduces functions that are made aware of the requested `size_tuple` of i.i.d samples, and does the broadcasting on the core shapes, transparently ignoring or moving the i.i.d `size_tuple` prepended axes around.

.. currentmodule:: pymc.distributions.shape_utils

.. autosummary::

   to_tuple
   shapes_broadcasting
   broadcast_dist_samples_shape
   get_broadcastable_dist_samples
   broadcast_distribution_samples
   broadcast_dist_samples_to

.. automodule:: pymc.distributions.shape_utils
   :members:
*************
Distributions
*************

.. toctree::

   distributions/continuous
   distributions/discrete
   distributions/multivariate
   distributions/mixture
   distributions/timeseries
   distributions/censored
   distributions/simulator
   distributions/transforms
   distributions/logprob
   distributions/utilities
*************
Sequential Monte Carlo
*************

.. toctree::

   smc/sample_smc
   smc/smc
*********
Inference
*********

Sampling
--------

.. currentmodule:: pymc.sampling

.. automodule:: pymc.sampling
   :members:


Step-methods
^^^^^^^^^^^^

.. currentmodule:: pymc.sampling

.. autofunction:: pymc.sampling.assign_step_methods


NUTS
~~~~

.. currentmodule:: pymc.step_methods.hmc.nuts

.. automodule:: pymc.step_methods.hmc.nuts
   :members:

Metropolis
~~~~~~~~~~

.. currentmodule:: pymc.step_methods.metropolis
odule:: pymc.step_methods.metropolis
   :members:
.. autom

Slice
~~~~~

.. currentmodule:: pymc.step_methods.slicer

.. automodule:: pymc.step_methods.slicer
   :members:

Hamiltonian Monte Carlo
~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: pymc.step_methods.hmc.hmc

.. autoclass:: pymc.step_methods.hmc.hmc.HamiltonianMC
   :members:

MultiTrace
^^^^^^^^^^

.. currentmodule:: pymc.backends.base
.. autoclass:: pymc.backends.base.MultiTrace
    :members:

.. autoclass:: pymc.backends.base.BaseTrace

Variational Inference
---------------------

OPVI
^^^^

.. currentmodule:: pymc.variational.opvi

.. automodule:: pymc.variational.opvi
   :members:

VI Inference API
^^^^^^^^^^^^^^^^

.. currentmodule:: pymc.variational.inference

.. automodule:: pymc.variational.inference
   :members:

Approximations
^^^^^^^^^^^^^^

.. currentmodule:: pymc.variational.approximations

.. automodule:: pymc.variational.approximations
   :members:

Operators
^^^^^^^^^

.. currentmodule:: pymc.variational.operators

.. automodule:: pymc.variational.operators
   :members:


Sequential Monte Carlo
----------------------

.. currentmodule:: pymc.smc.sample_smc

.. automodule:: pymc.smc.sample_smc
   :members:

Kernels
^^^^^^^

.. currentmodule:: pymc.smc.smc

.. automodule:: pymc.smc.smc
   :members:
====
Samplers
====

This submodule contains functions for MCMC sampling.


.. currentmodule:: pymc.sampling

.. autosummary::

.. automodule:: pymc.sampling
   :members:
**************************************
Ordinary differential equations (ODEs)
**************************************

This submodule contains tools used to perform inference on ordinary differential equations.

.. currentmodule:: pymc.ode

.. autosummary::

.. automodule:: pymc.ode
   :members: DifferentialEquation
Model
-----

.. currentmodule:: pymc.model
.. automodule:: pymc.model
   :members:
====
Math
====

This submodule contains various mathematical functions. Most of them are imported directly
from aesara.tensor (see there for more details). Doing any kind of math with PyMC random
variables, or defining custom likelihoods or priors requires you to use these Aesara
expressions rather than NumPy or Python code.

.. currentmodule:: pymc.math

.. autosummary::
   dot
   constant
   flatten
   zeros_like
   ones_like
   stack
   concatenate
   sum
   prod
   lt
   gt
   le
   ge
   eq
   neq
   switch
   clip
   where
   and_
   or_
   abs_
   exp
   log
   cos
   sin
   tan
   cosh
   sinh
   tanh
   sqr
   sqrt
   erf
   erfinv
   dot
   maximum
   minimum
   sgn
   ceil
   floor
   det
   matrix_inverse
   extract_diag
   matrix_dot
   trace
   sigmoid
   logsumexp
   invlogit
   logit

.. automodule:: pymc.math
   :members:
**********
Timeseries
**********

.. currentmodule:: pymc
.. autosummary::
   :toctree: generated

    AR1
    AR
    GaussianRandomWalk
    GARCH11
    EulerMaruyama
    MvGaussianRandomWalk
    MvStudentTRandomWalk
*******
Mixture
*******

.. currentmodule:: pymc
.. autosummary::
   :toctree: generated

   Mixture
   NormalMixture
   MixtureSameFamily
********
Censored
********

.. currentmodule:: pymc
.. autosummary::
   :toctree: generated

   Censored
************
Multivariate
************

.. currentmodule:: pymc
.. autosummary::
   :toctree: generated

   MvNormal
   MvStudentT
   Dirichlet
   Multinomial
   DirichletMultinomial
   OrderedMultinomial
   Wishart
   WishartBartlett
   LKJCorr
   LKJCholeskyCov
   MatrixNormal
   KroneckerNormal
   CAR
   StickBreakingWeights
**********************
Distribution utilities
**********************

.. currentmodule:: pymc
.. autosummary::
   :toctree: generated/

    Distribution
    SymbolicDistribution
    Discrete
    Continuous
    NoDistribution
    DensityDist
********
Discrete
********

.. currentmodule:: pymc
.. autosummary::
   :toctree: generated

   Binomial
   BetaBinomial
   Bernoulli
   DiscreteWeibull
   Poisson
   NegativeBinomial
   Constant
   ZeroInflatedPoisson
   ZeroInflatedBinomial
   ZeroInflatedNegativeBinomial
   DiscreteUniform
   Geometric
   HyperGeometric
   Categorical
   OrderedLogistic
   OrderedProbit
***************
Transformations
***************

.. currentmodule:: pymc.transforms

Transform Instances
~~~~~~~~~~~~~~~~~~~

Transform instances are the entities that should be used in the
``transform`` parameter to a random variable constructor.

.. autosummary::
   :toctree: generated

    simplex
    logodds
    interval
    log_exp_m1
    ordered
    log
    sum_to_1
    circular

Transform Composition Classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated

    Chain
    CholeskyCovPacked

Specific Transform Classes
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated

    LogExpM1
    Ordered
    SumTo1
**********
Simulator
**********

.. currentmodule:: pymc
.. autosummary::
   :toctree: generated

   SimulatorRV
   Simulator
   KullbackLiebler
**********
Continuous
**********

.. currentmodule:: pymc
.. autosummary::
   :toctree: generated/

   Uniform
   Flat
   HalfFlat
   Normal
   TruncatedNormal
   Beta
   Kumaraswamy
   Exponential
   Laplace
   StudentT
   Cauchy
   HalfCauchy
   Gamma
   Weibull
   HalfStudentT
   LogNormal
   ChiSquared
   HalfNormal
   Wald
   Pareto
   InverseGamma
   ExGaussian
   VonMises
   SkewNormal
   Triangular
   Gumbel
   Logistic
   LogitNormal
   Interpolated
   Rice
   Moyal
   AsymmetricLaplace
   PolyaGamma
****************************
Sequential Monte Carlo (SMC)
****************************

.. currentmodule:: pymc.smc.smc
.. autosummary::
   SMC_KERNEL
   IMH
   MH

.. automodule:: pymc.smc.smc
   :members:
***********************************
Sample Sequential Monte Carlo (SMC)
***********************************

.. currentmodule:: pymc.smc.sample_smc
.. autosummary::
   sample_smc

.. autofunction:: pymc.smc.sample_smc
**************
Mean Functions
**************

.. automodule:: pymc.gp.mean
.. autosummary::
   :toctree: generated

   Zero
   Constant
   Linear
***************
Implementations
***************

.. currentmodule:: pymc.gp
.. autosummary::
   :toctree: generated

   Latent
   LatentKron
   Marginal
   MarginalKron
   MarginalApprox
   TP
********************
Covariance Functions
********************

.. automodule:: pymc.gp.cov
.. autosummary::
   :toctree: generated

    Constant
    WhiteNoise
    ExpQuad
    RatQuad
    Exponential
    Matern52
    Matern32
    Linear
    Polynomial
    Cosine
    Periodic
    WarpedInput
    Gibbs
    Coregion
    ScaledCov
    Kron
