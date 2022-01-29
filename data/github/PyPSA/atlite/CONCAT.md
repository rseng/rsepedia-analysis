<!---
SPDX-FileCopyrightText: 2021 The Atlite Authors

SPDX-License-Identifier: CC0-1.0
--->

Atlite's contributor guidelines can be found in the official [documentation](https://atlite.readthedocs.io/en/master/contributing.html).<!--
SPDX-FileCopyrightText: 2021 The Atlite Authors

SPDX-License-Identifier: CC0-1.0
-->

Closes # (if applicable).

## Change proposed in this Pull Request

<!--- Provide a general, short summary of your changes in the title above -->

## Description
<!--- Describe your changes in detail -->

## Motivation and Context
<!--- Why is this change required? What problem does it solve? -->
<!--- If it fixes an open issue, please link to the issue here. -->

## How Has This Been Tested?
<!--- Please describe in detail how you tested your changes. -->
<!--- Include details of your testing environment, and the tests you ran to -->
<!--- see how your change affects other areas of the code, etc. -->

## Type of change
<!--- What types of changes does your code introduce? Put an `x` in all the boxes that apply: -->
- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to change)

## Checklist
<!--- Go over all the following points, and put an `x` in all the boxes that apply. -->
<!--- If you're unsure about any of these, don't hesitate to ask. We're here to help! -->
- [ ] I tested my contribution locally and it seems to work fine.
- [ ] I locally ran `pytest` inside the repository and no unexpected problems came up.
- [ ] I have adjusted the docstrings in the code appropriately.
- [ ] I have documented the effects of my code changes in the documentation `doc/`.
- [ ] I have added newly introduced dependencies to `environment.yaml` file.
- [ ] I have added a note to release notes `doc/release_notes.rst`.
- [ ] I have used `pre-commit run --all` to lint/format/check my contribution
---
name: Bug report
about: Create a report if something doesn't work quite right.
title: ''
labels: bug
assignees: ''
---
<!---
SPDX-FileCopyrightText: 2021 The Atlite Authors

SPDX-License-Identifier: CC0-1.0
--->


<!-- Provide a general summary of the issue -->

## Description
<!-- Provide a more detailed introduction to the issue itself, and why you consider it to be a bug -->
<!-- If you can, add a minimal example which reproduces the bug -->

## Expected Behavior
<!-- Tell us what should happen -->

## Actual Behavior
<!-- Tell us what goes wrong and happens instead -->

## Error Message
<!-- Paste any terminal output and error message you encounter here to help illustrate the problem -->

## Your Environment
<!-- Include relevant details about the environment you experienced the bug in -->
* The `atlite` version used:
* How you installed `atlite` (`conda`, `pip` or `github`):
* Operating System:
* My environment:
    <details>
      <summary>(output of `conda list`)</summary>
      ```
        <!-- output of `conda list` -->
      ```
    </details>


---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---
<!---
SPDX-FileCopyrightText: 2021 The Atlite Authors

SPDX-License-Identifier: CC0-1.0
--->

<!-- Provide a general summary of the feature you would like to see -->

## Detailed Description
<!-- Provide a detailed description of the change or addition you are proposing -->

## Context
<!-- Why is this change important to you? -->
<!-- How would you use it? -->

## Possible Implementation
<!-- Not obligatory, but suggest an idea for implementing addition or change -->
