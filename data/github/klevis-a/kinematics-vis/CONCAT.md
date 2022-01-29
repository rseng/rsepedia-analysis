Before proceeding with the test cases below, follow the instructions in the [README](README.md) to install the software and download the sample dataset.

### Test Cases

#### 1. Open a motion trial

**Steps**

1. Click the folder button (top of upper left hand quadrant) to open the trial selection dialog.
2. Select subject `O45_001_F_47_R` and `Scapular Plane Abduction`.
3. Click Analyze.

**Expected Outcome**

The UI should briefly flash `Loading...` then the trial analysis should be loaded.

* The upper left pane should have the `Humerus ISB: yx'y''` view.
* The upper right pane should have the `Scapula ISB: yx'z''` view.
* The lower left pane should have the `Preview` view.
* The lower right panel should contain plots and the default plot selected should be `Humerus ISB: yx'y''`.

#### 2. Orbital Controls

In all non-plot views the following controls should work:

1. Left-button dragging rotates the view.
2. Right-button dragging pans the view.
3. Middle-button scrolling zooms in and out.

#### 3. Motion capture frame preview (slider bar)

**Steps**

Drag the button of the slider bar (middle top of UI) to the right.

**Expected Outcome**

* The frame number indicator to the right of the slider bar should increase.
* All non-plot views should be updated. The slightly transparent scapula and humerus should be elevating.

#### 4. Analyzing a particular motion capture frame (slider bar)

**Steps**

Press the downward facing arrow to the left of the slider bar (once a motion capture frame has been selected via the slider bar).

**Expected Outcome**

* All views that start with Humerus or Scapula should be updated. The scapula view does not have a visual indication that updating completed successfully (subsequent test cases can validate that the update did complete successfully, however). The current longitude and latitude lines (indicated by yellow) on the Humerus views will update.

#### 5. Motion capture frame preview (plots)

**Steps**

Hover the mouse button over any plot curve. Trace the curve with the mouse pointer towards the right of the screen.

**Expected Outcome**

* The frame number indicator to the right of the slider bar should increase.
* All non-plot views should be updated. The slightly transparent scapula and humerus should be elevating.

#### 6. Analyzing a particular motion capture frame (plots)

**Steps**

Once a particular motion capture frame of interest has been determined (by hovering over the plot curves), click the left mouse button.

**Expected Outcome**

* All views that start with Humerus or Scapula should be updated. The scapula view does not have a visual indication that updating completed successfully (subsequent test cases can validate that the update did complete successfully, however). The current longitude and latitude lines (indicated by yellow) on the Humerus views will update.

#### 7. Switching the displayed view in a quadrant

**Steps**

Apart from the lower right quadrant, the views for all other quadrants can be switched via their drop-down selector. For each of the other 3 quadrants, make sure that each of the following views can be selected:

1. `Humerus ISB: yx'y''`
2. `Humerus Phadke: xz'y''`
3. `Humerus Swing Twist`
4. `Humerus Simultaneous`
5. `Preview`
6. `Scapula ISB: yx'z''`

**Expected Outcome**

All of the quadrants apart form the lower right one should update once the view that it displays is updated via the drop-down box.

#### 7. Correct number of steps for each view

**Steps**

This test case assures that the correct number of steps are associated with each view. The number of steps can be seen on the lower left corner of each view.

**Expected Outcome**

The following number of steps should be associated with views:

`Humerus ISB: yx'y''`: 3

`Humerus Phadke xz'y''`: 3

`Humerus Swing Twist`: 2

`Humerus Simultaneous`: 1

`Scapula ISB yxz''`: 3

`Preview`: None

#### 8. Euler angle step animations

**Steps**

Select any view that starts with Humerus or Scapula. Click each Step Number (lower left corner of view) to update the view. Click the Play button to animate a particular step. Click the last step for the particular view and animate until the end.

**Expected Outcome**

* When a step number is clicked the view is updated. The animation bar (to the right of the steps) is in its starting position after a step is clicked.
* When the play button is clicked the step animates and the animation bar slides to the right over a period of roughly 2 seconds.
* In the final step, once the animation has finished, the moving humerus/scapula will completely overlap the transparent humerus/scapula (which indicates the orientation of the humerus/scapula for that particular motion capture frame). This is how you know that the Humerus/Scapula view updated correctly when the downward facing arrow was clicked.

#### 9. Cycling through available plots

**Steps**

Use the drop-down selector of the lower right quadrant to cycle through available plots.

**Expected Outcome**

The following plots are available:

* Humerus Plane of Elevation
* Humerus Angle of Elevation
* Humerus Axial Orientation
* Humerus Axial Rotation
* Humerus ISB: yx'y''
* Humerus Phadke: xz'y''
* Humerus Swing Twist
* Scapula ISB: yx'z''

Selecting a plot updates the lower right quadrant and the title of the plot matches the selected value.

#### 10. Double-clicking to magnify a view

**Steps**

Double-click on any non-plot view. Once the view is magnified, double-click again to return to the standard view layout.

**Expected Outcome**

When a view is in a quadrant, double-clicking on it magnifies the view so it takes up the entire browser window. When the view is magnified, double-clicking restores the standard 4 pane layout.

#### 11. Getting Help

**Steps**

Click the `?` to the right of the frame selection bar (top middle of web app). Click the `X` in the upper right corner of the displayed help to close it.

**Expected Outcome**

* A brief description of how to interact with the software is displayed when the `?` button is clicked.
* The help closes when the `X` button is clicked.

#### 12. Additional options exclusive to Humerus views

**Steps**

Click the `Open Controls` button on the lower right corner of the web app. With a Humerus view displayed on any quadrant change the following options:

* Humerus Base
* Visualize Angles
* Spherical Area
* Show Sphere

Click `Close Controls`.

**Expected Outcome**

The Humerus view(s) will update as each option is changed. The controls should minimize when `Close Controls` is clicked.

#### 13. Additional options affecting both Humerus and Scapula views

**Steps**

Click the `Open Controls` button on the lower right corner of the web app. In both the Humerus and Scapula views select the final step for the view and animate to the end of the step. Then, change the following options:

* Show Triads/Arcs
* Prior Step Bones
* Show Body Planes

**Expected Outcome**

The Humerus and Scapula views will update as each option is changed.

#### 14. Humerus views are linked

**Steps**

In one quadrant select a Humerus view (say `Humerus ISB: yx'y''`). In another quadrant select another Humerus view (say `Humerus Phadke xz'y''`).

Now, rotate, pan, and zoom in/out of one view.

**Expected Outcome**

The other Humerus view will mimic the rotation, panning, and zooming of the controlled Humerus view.## Kinematics-Vis

This JavaScript application enables biomechanics researchers to visualize and analyze shoulder joint kinematics. It is built on top of [three.js](https://threejs.org/), a JavaScript 3D library. Although presently this application is specialized for analyzing and visualizing the shoulder joint, it should be easy to extend its functionality to other joints.

Kinematics-Vis has been peer-reviewed and accepted into the [Journal of Open Source Software](https://joss.theoj.org/). Checkout the associated paper: [![DOI](https://joss.theoj.org/papers/10.21105/joss.03490/status.svg)](https://doi.org/10.21105/joss.03490)

Also, checkout the [live code demo](https://shouldervis.chpc.utah.edu/kinevis/main.html) currently hosted at the [University of Utah Center for High Performance Computing](https://www.chpc.utah.edu/).

### Installation

This repository depends on [Yarn](https://github.com/yarnpkg/yarn) as a package manager. Please [install Yarn](https://yarnpkg.com/en/docs/install) before proceeding.

##### Clone repository
```
git clone https://github.com/klevis-a/kinematics-vis.git
cd kinematics-vis
```

##### Install dependencies and build

```
yarn install
yarn build
```

##### Download sample dataset

```
yarn fetch_data
```

##### Start `webpack` development server

```
yarn webpack serve
```

##### Access web application

[http://localhost:9000/main.html](http://localhost:9000/main.html)

### Usage
Instructions for interacting with the UI are provided within the web application. Once you access the web app click the question mark that appears at the top of the upper right quadrant. A simple way to provide an input dataset for the web app is to utilize the sample dataset (`yarn fetch_data`). Once the sample dataset has been downloaded, click the folder icon (top of the upper left quadrant), and select a trial to analyze.

![Help](help_pointer.png)

### Analyzing your own datasets
To analyze your own data see [INPUT_FILES.md](INPUT_FILES.md) for creating file formats compatible with this web app. To specify your own dataset directory edit the `DATA_DIR` variable within `webpack.config.js`.

### Manual Testing of the User Interface
See [MANUAL_TESTING.md](MANUAL_TESTING.md) for a list of test cases that cover the basic functionality of the user interface.

### Contributing

See the [CONTRIBUTING](CONTRIBUTING.md) document for details on contributing to the project by reporting a bug, submitting a fix, or proposing new features.# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
Klevis Aliaj (<klevis.a@gmail.com>).
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# Contributing to Kinematics-Vis
We want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

## We Develop with Github
We use github to host code, to track issues and feature requests, as well as accept pull requests.

## We Use [Github Flow](https://guides.github.com/introduction/flow/index.html), So All Code Changes Happen Through Pull Requests
Pull requests are the best way to propose changes to the codebase. We actively welcome your pull requests.

## Any contributions you make will be under the MIT Software License
In short, when you submit code changes, your submissions are understood to be under the same [MIT License](http://choosealicense.com/licenses/mit/) that covers the project.

## Report bugs using Github's [issues](https://github.com/klevis-a/kinematics-vis/issues)
We use GitHub issues to track public bugs. Report a bug by [opening a new issue](https://github.com/klevis-a/kinematics-vis/issues/new?assignees=&labels=&template=bug_report.md&title=); it's that easy!

## Write bug reports with detail, background, and sample code (if possible)

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give sample code if you can.
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

## License
By contributing, you agree that your contributions will be licensed under its MIT License.

## References
This document was adapted from the following open-source [contribution guidelines](https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62).### Format for Input Files
To visualize a kinematic trajectory the following files are needed:
1. Humerus ASCII STL 3D model file
2. Scapula ASCII STL 3D model file
3. A CSV (see format below) file specifying the position of the humeral head center (HHC), lateral epicondyle (LE), and medial epicondyle (ME) in millimeters, in the humerus STL model file coordinate system. The specified landmark order must be utilized.
4. A CSV (see format below) file specifying the position of the glenoid center (GC), inferior angle (IA), trigonum spinae (TS), posterior lateral acromion (PLA), and acromion (AC) in millimeters, in the scapula STL model file coordinate system. The specified landmark order must be utilized.
5. A CSV file (see format below) specifying the pose trajectory of the torso, scapula, and humerus.

The provided sample dataset (downloadData.sh and downloadData.bat) contains many examples of the aforementioned files. Furthermore, if you do not have patient-specific bone models you may use any of the STL 3D models (and associated landmarks files) provided in the sample dataset.

##### Humerus landmarks file format

| Landmark | X       | Y       | Z       |
| -------- | ------- | ------- | ------- |
| HHC      | -108.61 | -124.11 | -112.29 |
| LE       | -157.13 | -86.69  | -384.3  |
| ME       | -97.39  | -91.46  | -394.46 |

##### Scapula landmarks file format

| Landmark | X       | Y       | Z       |
| -------- | ------- | ------- | ------- |
| GC       | -87.01  | -112.83 | -116.54 |
| IA       | -21.19  | -42.55  | -211.43 |
| TS       | -8.66   | -57.25  | -106.92 |
| PLA      | -119.27 | -99.63  | -90.29  |
| AC       | -98.14  | -128.33 | -71.56  |

##### Pose trajectory file specification

The pose trajectory file must contain columns for the position and orientation (scalar-last quaternion) for the torso, scapula, and humerus - in that order.

* torso_pos_x
* torso_pos_y
* torso_pos_z
* torso_quat_x
* torso_quat_y
* torso_quat_z
* torso_quat_w
* scapula_pos_x
* scapula_pos_y
* scapula_pos_z
* ...
* humerus_quat_x
* humerus_quat_y
* humerus_quat_z
* humerus_quat_w

Please note that the torso, scapula, and humerus must all be expressed in the same underlying coordinate system (e.g. biplane fluoroscopy or Vicon).

Each row of the pose trajectory file represents a time sample of the motion trajectory. The software assumes that time samples are equally spaced in time - i.e., constant sampling frequency. See details below about how this sampling frequency is defined.

### Specification of Input Files to Analyze

There are two primary ways to specify the input files to the web application: 1) via a database put together via a JSON file, 2) via manual file uploads. This section explores these two methods in more detail. The provided sample dataset (downloadData.sh or downloadData.bat) provides an example of a database put together via a JSON file and can be utilized as a reference for the rest of this section.

##### Specifying Input Files via a JSON file

This method is only possible if you have the ability to upload files to the server where the web app is running or if you are running the web app locally. Under the root application directory, create a `csv` directory and under the `csv` directory create a `db_summary.json` file. Within the `db_summary.json` file, the key for each top-level entry (object) specifies a subject identifier. Within each subject, `config` and `activities` entries must exist. Furthermore, one subject must be set as the default one that is selected when the subject selection dialog is rendered. This is done by setting the `default` entry for this subject to `1`.

The `config` entry specifies files that pertain to a particular subject and do not vary between the activities that the subject performed. Specifically the following files must be specified in relationship to the previously created `csv` directory:

* `humerus_landmarks_file` - the location of the humerus landmarks file in relationship to the `csv` directory.
* `scapula_landmarks_file` -  the location of the scapula landmarks file in relationship to the `csv` directory.
* `humerus_stl_file` - the location of the humerus ASCII STL 3D model file in relationship to the `csv` directory.
* `scapula_stl_file` - the location of the scapula ASCII STL 3D model file in relationship to the `csv` directory.

The `activities` entry contains pose trajectory files for one or more activities that the subject performed. The order that these activities are listed in this file will be preserved in the UI. Special characters for activity names may be specified via their corresponding HTML codes. For each activity, the following entries must be present:

* `trajectory` - the location of the pose trajectory file for this activity in relationship to the `csv` directory.
* `freq` - the capture frequency (in Hz) for this activity.

##### Specifying Input Files via Manual Upload

In the UI, click the folder icon (top of upper left quadrant) to open the database selection dialog. Click **Manual File Selection**. The files discussed in the **Specifying Input Files via a JSON file** section and the capture frequency can be specified here manually.Fixes #

## Proposed Changes

  -
  -
  -
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
