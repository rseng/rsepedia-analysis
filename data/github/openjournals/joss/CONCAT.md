## JOSS backlog

> tl;dr This is an attempt to capture a prioritized backlog for JOSS, noting the rough level of difficulty for each item. The intent is that after feedback from @openjournals/joss-editors and @openjournals/dev we'll go ahead and open a bunch of issues to track items here, and manage them in a GitHub Project board.

## Author improvements

[Author improvements across JOSS, Whedon, Whedon-API](https://github.com/search?q=org%3Aopenjournals+label%3Aenhancement%3Aauthor&type=Issues)

### High-priority issues

- Whedon should return a better error when it can't git clone: https://github.com/openjournals/whedon/issues/60
- Whedon should check paper structure: https://github.com/openjournals/whedon-api/issues/54
- Help authors pick from a pre-defined list of tags: https://github.com/openjournals/joss/issues/677
- Author dashboard for a paper to show its status: https://github.com/openjournals/joss/issues/514

## Reviewer improvements

[Reviewer improvements across JOSS, Whedon, Whedon-API](https://github.com/search?q=org%3Aopenjournals+label%3Aenhancement%3Areviewer&type=Issues)

### High-priority issues

- Managing reviewer workload: https://github.com/openjournals/joss/issues/436
- Automate reviewer suggestions: https://github.com/openjournals/whedon-api/issues/2
- Better acknowledgement of reviewers: https://github.com/openjournals/joss/issues/624
- Show reviewer names on papers (rather than GitHub handles): https://github.com/openjournals/joss/issues/667

## Editor improvements

[Editor improvements across JOSS, Whedon, Whedon-API](https://github.com/search?q=org%3Aopenjournals+label%3Aenhancement%3Aeditor&type=Issues)

### High-priority issues

- Whedon should help editors identify reviewers: https://github.com/openjournals/whedon-api/issues/2
- Whedon should be able to update the repository address: https://github.com/openjournals/whedon-api/issues/75
- Whedon should be able to add (or remove) a reviewer mid-review: https://github.com/openjournals/whedon-api/issues/65
- Whedon should not allow arbitrary string for reviewer names: https://github.com/openjournals/whedon-api/issues/30

## General improvements

[Reader improvements across JOSS, Whedon, Whedon-API](https://github.com/search?q=org%3Aopenjournals+label%3Aenhancement%3Areader&type=Issues)  
[General improvements across JOSS, Whedon, Whedon-API](https://github.com/search?q=org%3Aopenjournals+label%3Aenhancement%3Ameta&type=Issues)

### High-priority issues

- JOSS should produce JATS XML for indexers (e.g. PubMed Central) https://github.com/openjournals/whedon/issues/36
- JOSS website should show publication statistics: https://github.com/openjournals/joss/issues/532
- Improved docs on our editorial management policies: https://github.com/openjournals/joss/issues/492
# The Journal of Open Source Software

[![Build Status](https://github.com/openjournals/joss/actions/workflows/tests.yml/badge.svg)](https://github.com/openjournals/joss/actions/workflows/tests.yml)
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](http://numfocus.org)
[![Donate to JOSS](https://img.shields.io/badge/Donate-to%20JOSS-brightgreen.svg)](https://numfocus.org/donate-to-joss)

The [Journal of Open Source Software](https://joss.theoj.org) (JOSS) is a developer friendly journal for research software packages.

### What exactly do you mean by 'journal'

The Journal of Open Source Software (JOSS) is an academic journal with a formal peer review process that is designed to _improve the quality of the software submitted_. Upon acceptance into JOSS, a CrossRef DOI is minted and we list your paper on the JOSS website.

### Don't we have enough journals already?

Perhaps, and in a perfect world we'd rather papers about software weren't necessary but we recognize that for most researchers, papers and not software are the currency of academic research and that citations are required for a good career.

We built this journal because we believe that after you've done the hard work of writing great software, it shouldn't take weeks and months to write a paper<sup>1</sup> about your work.

### You said developer friendly, what do you mean?

We have a simple submission workflow and extensive documentation to help you prepare your submission. If your software is already well documented then paper preparation should take no more than an hour.

<sup>1</sup> After all, this is just advertising.

## The site

The JOSS submission tool is hosted at https://joss.theoj.org

## JOSS Reviews

If you're looking for the JOSS reviews repository head over here: https://github.com/openjournals/joss-reviews/issues

## Code of Conduct

In order to have a more open and welcoming community, JOSS adheres to a code of conduct adapted from the [Contributor Covenant](http://contributor-covenant.org) code of conduct.

Please adhere to this code of conduct in any interactions you have in the JOSS community. It is strictly enforced on all official JOSS repositories, the JOSS website, and resources. If you encounter someone violating these terms, please let the Editor-in-Chief ([@arfon](https://github.com/arfon)) or someone on the [editorial board](https://joss.theoj.org/about#editorial_board) know and we will address it as soon as possible.

## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Added some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create new Pull Request

## ⚙️ Development

[LiveReload](https://github.com/guard/guard-livereload) enables the browser to automatically refresh on change during development.

1. Download the [LiveReload Chrome plugin](https://chrome.google.com/webstore/detail/livereload/jnihajbhpnppcggbcgedagnkighmdlei/)
2. Run `bundle exec guard`
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

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

The JOSS Editors are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

JOSS Editors have the right and responsibility to remove, edit, or
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
reported by contacting the JOSS editor-in-chief at <arfon.smith@gmail.com>. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# JOSS Conflict of Interest Policy

The definition of a conflict of Interest in peer review is a circumstance that makes you "unable to make an impartial scientific judgment or evaluation." ([PNAS Conflict of Interest Policy](http://www.pnas.org/site/authors/coi.xhtml)). JOSS is concerned with avoiding any actual conflicts of interest, and being sufficiently transparent that we avoid the appearance of conflicts of interest as well.

As a reviewer, COIs are your present or previous association with any authors of a submission: recent (past four years) collaborators in funded research or work that is published; and lifetime for the family members, business partners, and thesis student/advisor or mentor. In addition, your recent (past year) association with the same organization of a submitter is a COI, for example, being employed at the same institution.

If you have a conflict of interest with a submission, you should disclose the specific reason to the submissions' editor. This may lead to you not being able to review the submission, but some conflicts may be recorded and then waived, and if you think you are able to make an impartial assessment of the work, you should request that the conflict be waived. For example, if you and a submitter were two of 2000 authors of a high energy physics paper but did not actually collaborate. Or if you and a submitter worked together 6 years ago, but due to delays in the publishing industry, a paper from that collaboration with both of you as authors was published 2 year ago. Or if you and a submitter are both employed by the same very large organization but in different units without any knowledge of each other.

Declaring actual, perceived, and potential conflicts of interest is required under professional ethics. If in doubt: ask the editors.
### 3.4.1 (January 24, 2016)

This includes various SVG viewport refinements.

Refines:

- `thumbs-down`
- `logo-github`

### 3.4.0 (January 22, 2016)

Adds:

- `verified`
- `smiley`

Removes:

- `color-mode`

Refines:

- `primitive-dot`
- `horizontal-rule`
- `triangle-down`
- `triangle-up`
- `triangle-left`
- `triangle-right`
- `globe`
- `flame`
- `comment-discussion`

### 3.3.0 (November 12, 2015)

Adds:

- `logo-gist`

Resizes all our SVG to be 16x16 instead of 1024x1024

### 3.2.0 (November 6, 2015)

Adds:

- `bold`
- `text-size`
- `italic`
- `tasklist`

It also normalizes some styling in:

- `list-unordered`
- `list-ordered`
- `quote`
- `mention`
- `bookmark`
- `threebars`

Removes

- `screen-normal`
- `screen-full`


### 3.1.0 (August 13, 2015)

Adds

- `shield`

This thickens stroke widths slightly on the following icons:

- `circle-slash`
- `clock`
- `cloud-upload`
- `cloud-download`
- `dashboard`
- `info`
- `issue-closed`
- `issue`
- `issue-reopened`
- `history`
- `question`
- `search`

Fills `comment-discussion`

Thickens `x` to match `checkmark`

### 3.0.1 (August 10, 2015)

Some files were missing in `3.0.0`

### 3.0.0 (August 10, 2015)

Removes

- `microscope`
- `beer`
- `split`
- `puzzle`
- `steps`
- `podium`
- `timer`
- all `alignment` icons
- all `move` icons
- all `playback` icons
- all `jump` icons

Adds

- `beaker`
- `bell`
- `desktop-download`
- `watch`

Line-weight changes, sizing normalization, and new drawings

- `circle-slash`
- `lock`
- `cloud-upload`
- `cloud-download`
- `plus`
- `✕`
- `broadcast`
- `lock`
- all `repo` icons
- organization
- person
- all `chevrons` & `triangles`
- all `diff` icons
- `clippy`
- all `issue` and circular icons
- `rss`
- `ruby`
- `cancel`
- `settings`
- `mirror`
- `external-link`
- `history`
- `gear`
- `settings`
- `info`
- `history`
- `package`
- `gist-secret`
- `rocket`
- `law`
- `telescope`
- `search`
- `tag`
- `normal-screen`
- `iphone`
- `no-new-line`
- `desktop`
- all `git` icons
- `circuit-board`
- `heart`
- `home`
- `briefcase`
- `wiki`
- `bookmark`
- `briefcase`
- `calendar`
- `color-mode`
- `comment`
- `discussions`
- `credit-card`
- `dashboard`
- `camera`
- `video`
- `bug`
- `desktop`
- `ellipses`
- `eye`
- all `files` & `folders`
- `fold`
- `unfold`
- `gift`
- `graph`
- `hubot`
- `inbox`
- `jersey`
- `keyboard`
- `light-bulb`
- `link`
- `location`
- `mail`
- `mail-read`
- `marker`
- `plug`
- `mute`
- `pencil`
- `push-pin`
- `fullscreen`
- `unfullscreen`
- `server`
- `sign-in`
- `sign-out`
- `tag`
- `terminal`
- `thumbs-up`
- `thumbs-down`
- `trash`
- `unmute`
- `versions`
- `gist`
- `key`
- `megaphone`
- `checklist`

## 2.4.1 (June 2, 2015)

- Add the scss file I forgot to include

## 2.4.0 (June 2, 2015)

- Add `octicons.scss`
- Revert path changes to `sprockets-octicons.scss`, as they broke octicons in sprockets.

## 2.3.0 (May 28, 2015)

- Add a path variable to `sprockets-octicons.scss` to be consistent with octicons.less`

## 2.2.3 (May 21, 2015)

- Use SPDX license identifiers in package.json

## 2.2.2 (April 1, 2015)

Fixes file icons for

- `file-binary`
- `file-code`
- `file-media`
- `file-pdf`
- `file-symlink-file`
- `file-text`
- `file-zip`

## 2.2.1 (March 30, 2015)

- Fix vector artifact and smooth curves in `mark-github`

## 2.2.0 (Feb 18, 2015)

- Add two new icons: `thumbsup` and `thumbsdown`

## 2.0.1 (June 16, 2014)

- Add mention of github.com/logos to the license

## 2.0.0 (June 16, 2014)

- Hello world
# Octicons!

This is the [Bower][bower] package for [GitHub Octicons][octicons].

## Add Octicons to your project

1. Create a new file called *bower.json* (if you don't have one already).

2. Add a new line for the Octicon dependency, pointing to the correct repository:

  ``` json
  {
    "name": "my_great_project",
    "dependencies": {
      "octicons": "*"
    }
  }
  ```

3. Run `bower install`. The Octicons styles will be downloaded to *bower_components/octicons*.

4. Link to the `octicons.css` stylesheet in the `<head>` of your `<html>` page:

  ``` html
  <link rel="stylesheet" href="bower_components/octicons/octicons/octicons.css">
  ```

4. Simply use an icon in your HTML page:

  ``` html
  <span class="octicon octicon-microscope"></span>
  ```

### Rails' asset pipeline

Octicons includes a stylesheet specifically for [Rails 4/Sprockets][sprockets].

1. Create a new file called *vendor/assets/bower.json* (if you don't have one already).

2. Add a new line for the Octicon dependency, pointing to the correct repository:

  ``` json
  {
    "name": "my_great_project",
    "dependencies": {
      "octicons": "*"
    }
  }
  ```

3. `cd` into `vendor/assets` and run `bower install`. The Octicons styles will be downloaded to *vendor/assets/bower_components/octicons*.

4. Open your config/application.rb, and add this line inside your Application:

  ``` ruby
  config.assets.precompile += %w(*.svg *.eot *.woff *.ttf)
  ```

5. In your application stylesheet, require `sprockets-octicons`:

  ``` css
  /*
  = require sprockets-octicons
  */
  ```

6. Simply use an icon in your HTML page:

  ``` html
  <span class="octicon octicon-flame"></span>
  ```

7. If you want a view helper, add something like this to *app/helpers/application_helper.rb*:

  ``` ruby
  def octicon(code)
    content_tag :span, '', class: "octicon octicon-#{code.to_s.dasherize}"
  end
  ```

## Installing locally

It's easy to install octicons locally if you have [Homebrew](http://brew.sh/) installed. Simply run the following commands:

```
brew install caskroom/cask/brew-cask
brew tap "caskroom/fonts"
brew cask install "font-octicons"
```

## Best practices

- Octicons look best in sizes that are multiples of 16px. You can update the size using the `font-size` CSS property. For example:

  ``` css
  .octicon {
    font-size: 32px;
  }
  ```

- Octicons are not monospaced. This lets them work well next to type, but it means they won’t stack nicely by default. If you intend to stack octicons, such as in navigation, you will want to add some CSS to make them the same width, and centered. For example:

  ``` css
  .navigation .octicon {
    width: 16px;
    text-align: center;
  }
  ```

### Resources

- [octicons.github.com](http://octicons.github.com/) - the Octicons website
- Read why [icon fonts are awesome](http://css-tricks.com/examples/IconFont/)
- How to compose your [HTML for icon font usage](http://css-tricks.com/html-for-icon-font-usage/)
- [sketch-octicons](https://github.com/JuanitoFatas/sketch-octicons) - Octicons icons as Sketch Symbols

## Why can't I see the characters in Font Book??

Give this a try, you should be all set:

![](http://cl.ly/image/2r1B1F2l3Q0D/content#png)

## FAQ

Check out [issues with the FAQ label](https://github.com/github/octicons/issues?q=is%3Aclosed+is%3Aissue+label%3AFAQ).

## Versions

Octicons operates similarly to [Semver](http://semver.org/) with the following version convention:

- Major: Breaking changes — removed icons, markup changes, unicode switches, css renames, icon redesigns
- Minor: Non-breaking changes — new icons, new aliases, minor icon changes
- Patch: Unnoticeable tweaks — slight visual changes, package updates


[octicons]: http://octicons.github.com
[bower]: http://bower.io/
[sprockets]: http://guides.rubyonrails.org/asset_pipeline.html
The contents of */octicons* */svg* are generated by an automated process. Changes to these files may be accepted, but may also be overwritten.

Octicons is GitHub's icon font. At this time, new icons will only add icons when they are needed for GitHub products.
If you intend to install Octicons locally, install `octicons-local.ttf`. It should appear as “github-octicons” in your font list. It is specially designed not to conflict with GitHub's web fonts.
# Maintaining

Steps for updating and releasing changes to Primer and it's site.

## Versioning

Primer follows the semantic versioning approach:

- Bug fixes and docs updates are patch releases, so `1.0.x`.
- New additions are minor updates, so `1.x.x`.
- Deleting or rewriting anything are major updates, so `x.x.x`.

## Changelogs and milestones

Changelogs are handled with dedicated tracking issues ([see example](https://github.com/primer/primer/issues/108)). When starting work on a new release:

1. Open a new milestone.
2. Open a new tracking issue and immediately lock it. (No comments are needed, ship lists are just for us.)
3. As you close issues and merge pull requests, add a link to those threads to the tracking issue.

When the release and milestone are about ready to ship, move on the the releasing flow.

## Releasing

Have a new version to release? Hell yeah, let's do it.

1. Bump the version numbers in `_config.yml` for our docs and `package.json` for dependency management.
2. Run `$ grunt` to generate the latest compiled CSS and Parker stats.
3. Recompile Jekyll for the latest docs changes.
4. Punt any remaining open issues and PRs on the milestone to the next milestone, then close that milestone.
5. Head to <https://github.com/primer/primer/releases/> and create a new release. Title it `vX.X.X` and post the changelog to the body.
6. Run `$ grunt publish` to push the latest docs and CSS changes to <http://primercss.io>.
7. Rejoice!
# Primer

Primer is the CSS toolkit that powers GitHub's front-end design. It's purposefully limited to common components to provide our developers with the most flexibility, and to keep GitHub uniquely *GitHubby*. It's built with SCSS and available via Bower, so it's easy to include all or part of it within your own project.

[**Read the Primer documentation**](http://primercss.io) to learn more.

_**Heads up!** We love open source, but Primer is unlikely to add new features that are not used in GitHub.com. It's first and foremost our CSS toolkit. We really love to share though, so hopefully that means we're still friends <3._

## Contents

- [Install](#install)
- [Usage](#usage)
- [Documentation](#documentation)
  - [Dependencies](#dependencies)
  - [Running locally](#running-locally)
  - [Publishing](#publishing)
  - [Primer stats](#primer-stats)
- [Updating](#updating)
- [Contributing](#contributing)
- [Versioning](#versioning)
- [License](#license)

## Install

### Manually

Download the [latest release](https://github.com/primer/primer/releases/latest) and copy the SCSS files over to your own project. Once your files are in place, jump to the [usage guidelines](#usage) for including Primer into your own CSS.

### Bower

```
$ bower install primer-css --save
```

### Things to know

**Hey, GitHubbers!** For GitHub.com, you'll need to  `cd` into `vendor/assets` and run `bower install` there. Be sure to commit and push all the changes, including the `bower.json` and everything under `bower_components`.

## Usage

Once included, simply `@import` either the master SCSS file, or the individual files as you need them.

```scss
// Example: All of Primer
@import "primer-css/scss/primer";

// Example: Individual files
@import "primer-css/scss/variables";
@import "primer-css/scss/mixins";
@import "primer-css/scss/base";
```

## Documentation

Primer's documentation is built with Jekyll and published to `http://primercss.io` via the `gh-pages` branch.

### Dependencies

You'll need the following installed:

- Latest Jekyll (minimum v2.2.0): `$ gem install jekyll`
- Latest Rouge: `$ gem install rouge`
- Latest Sass: `$ gem install sass`
- Latest Grunt CLI: `$ npm install -g grunt-cli`
- [Node.js and npm](http://nodejs.org/download/)

Chances are you have all this already if you work on `github/github` or similar projects. If you have all those set up, now you can install the dependencies:

```bash
$ npm install
$ bower install
```

### Running locally

From the Terminal, start a local Jekyll server:

```bash
$ jekyll serve
```

Open a second Terminal tab to automatically recompile the Sass files, run autoprefixer, and update our [Primer stats file](#primer-stats):

```bash
$ grunt watch
```

Alternatively, you can manually run `grunt` and `jekyll serve` when needed.

### Publishing

Use the included Grunt task to generate and publish Primer's docs to the `gh-pages` branch.

```bash
$ grunt publish
```

This takes the `_site` directory, generates it's own Git repository there, and publishes the contents to the `gh-pages` branch here on GitHub. Changes are reflected in the hosted docs within a minute or so.

### Primer stats

When compiling or watching the Sass files, Primer will automatically generate a `.primer-stats.md` file. This is tracked in the Git repository to provide us historical and contextual information on the changes we introduce. For example, we'll know when the number of selectors or declarations rises sharply within a single change.

## Updating

Within `bower.json`, update to a new release by changing the version number that follows the `#` in the dependency URL.

```json
{
  "name": "myapp",
  "dependencies": {
    "primer-css": "x.x.x"
  }
}
```

To pull down the updated package, `cd` into `vendor/assets`, and run `bower install`.

```
$ cd vendor/assets
$ bower install
```

Check in `bower.json` and all changes under `vendor/assets/bower_components`.

## Development

Development of Primer happens in our primary branch, `master`. For stable versions, see the [releases page](https://github.com/primer/primer/releases). `master` will always be up to date with the latest changes, including those which have yet to be released.

## Contributing

By contributing to Primer, you agree to the terms presented in [this license agreement](https://cla.github.com/). *More information will be provided here soon.*

When contributing changes to Primer, be sure to do the following steps when opening a pull request:

1. Bump the version number in `bower.json` (it's purely placebo right now, but it's good habit) and `package.json`.
2. Run `grunt` and commit the changes. This compiles the SCSS to CSS so we can do basic analysis on the number of selectors, file size, etc.

In addition, please read through our [contributing guidelines](https://github.com/primer/primer/blob/master/CONTRIBUTING.md). Included are directions for opening issues, coding standards, and notes on development.

All HTML and CSS should conform to the [style guidelines](http://primercss.io/guidelines).

Editor preferences are available in the [editor config](https://github.com/primer/primer/blob/master/.editorconfig) for easy use in common text editors. Read more and download plugins at <http://editorconfig.org>.

## Versioning

For transparency into our release cycle and in striving to maintain backward compatibility, Primer is maintained under [the Semantic Versioning guidelines](http://semver.org/). Sometimes we screw up, but we'll adhere to those rules whenever possible.

## License

Created by and copyright GitHub, Inc. Released under the [MIT license](LICENSE.md).
## Contributing

[fork]: https://github.com/github/primer/fork
[pr]: https://github.com/github/primer/compare
[style]: http://primercss.io/guidelines/

Hi there! We're thrilled that you'd like to contribute to this project. Your help is essential for keeping it great.

After you open your first pull request, you will be asked to accept [this license agreement](https://cla.github.com/). Let us know in the PR if you have any hesitation or concerns.

## Using the issue tracker

The [issue tracker](https://github.com/primer/primer/issues) is the preferred channel for [bug reports](#bug-reports), [features requests](#feature-requests) and [submitting pull requests](#pull-requests), but please respect the following restrictions:

* Please **do not** use the issue tracker for personal support requests.
* Please **do not** derail or troll issues. Keep the discussion on topic and respect the opinions of others.
* Please **do not** open issues or pull requests regarding the code in [`Normalize`](https://github.com/necolas/normalize.css) (open them in their respective repositories).
  
## Bug reports

A bug is a _demonstrable problem_ that is caused by the code in the repository. Good bug reports are extremely helpful, so thanks!

Guidelines for bug reports:

0. **Validate and lint your code** &mdash; [validate your HTML](http://html5.validator.nu) to ensure your problem isn't caused by a simple error in your own code.

1. **Use the GitHub issue search** &mdash; check if the issue has already been reported.

2. **Check if the issue has been fixed** &mdash; try to reproduce it using the latest `master` or development branch in the repository.

3. **Isolate the problem** &mdash; ideally create a [reduced test case](https://css-tricks.com/reduced-test-cases/) and a live example. [This JS Bin](http://jsbin.com/lefey/1/edit?html,output) is a helpful template.

A good bug report shouldn't leave others needing to chase you up for more information. Please try to be as detailed as possible in your report. What is your environment? What steps will reproduce the issue? What browser(s) and OS experience the problem? Do other browsers show the bug differently? What would you expect to be the outcome? All these details will help people to fix any potential bugs.

Example:

> Short and descriptive example bug report title
>
> A summary of the issue and the browser/OS environment in which it occurs. If
> suitable, include the steps required to reproduce the bug.
>
> 1. This is the first step
> 2. This is the second step
> 3. Further steps, etc.
>
> `<url>` - a link to the reduced test case
>
> Any other information you want to share that is relevant to the issue being reported. This might include the lines of code that you have identified as causing the bug, and potential solutions (and your opinions on their merits).

## Feature requests

Feature requests are welcome. But take a moment to find out whether your idea fits with the scope and aims of the project. It's up to *you* to make a strong case to convince the project's developers of the merits of this feature. Please provide as much detail and context as possible.

## Pull requests

Good pull requests—patches, improvements, new features—are a fantastic help. They should remain focused in scope and avoid containing unrelated commits.

**Please ask first** before embarking on any significant pull request (e.g. implementing features, refactoring code, porting to a different language), otherwise you risk spending a lot of time working on something that the project's developers might not want to merge into the project.

Adhering to the following process is the best way to get your work included in the project:

1. Fork and clone the repository
2. Configure and install the dependencies: `bower install`
3. Create a new branch: `git checkout -b my-branch-name`
4. Make your change, add tests, and make sure the tests still pass
5. Push to your fork and [submit a pull request](https://help.github.com/articles/creating-a-pull-request/)
6. Pat your self on the back and wait for your pull request to be reviewed and merged.

Here are a few things you can do that will increase the likelihood of your pull request being accepted:

- Follow the [style guide][style].
- Keep your change as focused as possible. If there are multiple changes you would like to make that are not dependent upon each other, consider submitting them as separate pull requests.
- Write a [good commit message](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html).

## Resources

- [Contributing to Open Source on GitHub](https://guides.github.com/activities/contributing-to-open-source/)
- [Using Pull Requests](https://help.github.com/articles/using-pull-requests/)
- [GitHub Help](https://help.github.com)
# [primer-css]( http://primercss.io )

**Version:** `2.3.5`

> Primer is the CSS toolkit that powers GitHub's front-end design. It's purposefully limited to common components to provide our developers with the most flexibility, and to keep GitHub uniquely *GitHubby*. It's built with SCSS and available via Bower, so it's easy to include all or part of it within your own project.

* * *

## Parker Report

### css/primer.css

- **Total Stylesheets:** 1
- **Total Stylesheet Size:** 27736
- **Total Media Queries:** 1
- **Total Rules:** 363
- **Selectors Per Rule:** 1.4986225895316805
- **Total Selectors:** 544
- **Identifiers Per Selector:** 2.2316176470588234
- **Specificity Per Selector:** 16.762867647058822
- **Top Selector Specificity:** 42
- **Top Selector Specificity Selector:** dl.form.warn dd.warning:after
- **Total Id Selectors:** 0
- **Total Identifiers:** 1214
- **Total Declarations:** 915
- **Total Unique Colors:** 81
- **Total Important Keywords:** 1
Submitting a paper to JOSS
==========================

If you've already developed a fully featured research code, released it under an [OSI-approved license](https://opensource.org/licenses), and written good documentation and tests, then we expect that it should take perhaps an hour or two to prepare and submit your paper to JOSS.
But please read these instructions carefully for a streamlined submission.

## Submission requirements

- The software must be open source as per the [OSI definition](https://opensource.org/osd).
- The software must have an **obvious** research application.
- You must be a major contributor to the software you are submitting, and have a GitHub account to participate in the review process.
- Your paper must not focus on new research results accomplished with the software.
- Your paper (`paper.md` and BibTeX files, plus any figures) must be hosted in a Git-based repository together with your software (although they may be in a short-lived branch which is never merged with the default).

In addition, the software associated with your submission must:

- Be stored in a repository that can be cloned without registration.
- Be stored in a repository that is browsable online without registration.
- Have an issue tracker that is readable without registration.
- Permit individuals to create issues/file tickets against your repository.

### What we mean by research software

JOSS publishes articles about research software. This definition includes software that: solves complex modeling problems in a scientific context (physics, mathematics, biology, medicine, social science, neuroscience, engineering); supports the functioning of research instruments or the execution of research experiments; extracts knowledge from large data sets; offers a mathematical library, or similar. While useful for many areas of research, pre-trained machine learning models and notebooks are not in-scope for JOSS. 

### Substantial scholarly effort

JOSS publishes articles about software that represent substantial scholarly effort on the part of the authors. Your software should be a significant contribution to the available open source software that either enables some new research challenges to be addressed or makes addressing research challenges significantly better (e.g., faster, easier, simpler).

As a rule of thumb, JOSS' minimum allowable contribution should represent **not less than** three months of work for an individual. Some factors that may be considered by editors and reviewers when judging effort include:

- Age of software (is this a well-established software project) / length of commit history.
- Number of commits.
- Number of authors.
- Total lines of code (LOC). Submissions under 1000 LOC will usually be flagged, those under 300 LOC will be desk rejected.
- Whether the software has already been cited in academic papers.
- Whether the software is sufficiently useful that it is _likely to be cited_ by your peer group.

In addition, JOSS requires that software should be feature-complete (i.e., no half-baked solutions), packaged appropriately according to common community standards for the programming language being used (e.g., [Python](https://packaging.python.org), [R](https://r-pkgs.org/index.html)), and designed for maintainable extension (not one-off modifications of existing tools). "Minor utility" packages, including "thin" API clients, and single-function packages are not acceptable.

### Co-publication of science, methods, and software

Sometimes authors prepare a JOSS publication alongside a contribution describing a science application, details of algorithm development, and/or methods assessment. In this circumstance, JOSS considers submissions for which the implementation of the software itself reflects a substantial scientific effort. This may be represented by the design of the software, the implementation of the algorithms, creation of tutorials, or any other aspect of the software. We ask that authors indicate whether related publications (published, in review, or nearing submission) exist as part of submitting to JOSS.

#### Other venues for reviewing and publishing software packages

Authors wishing to publish software deemed out of scope for JOSS have a few options available to them:

- Follow [GitHub's guide](https://guides.github.com/activities/citable-code/) on how to create a permanent archive and DOI for your software. This DOI can then be used by others to cite your work.
- Enquire whether your software might be considered by communities such as [rOpenSci](https://ropensci.org) and [pyOpenSci](https://pyopensci.org).

### Should I write my own software or contribute to an existing package?

While we are happy to review submissions in standalone repositories, we also review submissions that are significant contributions made to existing packages. It is often better to have an integrated library or package of methods than a large number of single-method packages.

### Questions? Open an issue to ask

Authors wishing to make a pre-submission enquiry should [open an issue](https://github.com/openjournals/joss/issues/new?title=Pre-submission%20enquiry) on the JOSS repository.

## Typical paper submission flow

Before you submit, you should:

- Make your software available in an open repository (GitHub, Bitbucket, etc.) and include an [OSI approved open source license](https://opensource.org/licenses).
- Make sure that the software complies with the [JOSS review criteria](review_criteria.html). In particular, your software should be full-featured, well-documented, and contain procedures (such as automated tests) for checking correctness.
- Write a short paper in Markdown format using `paper.md` as file name, including a title, summary, author names, affiliations, and key references. See our [example paper](#example-paper-and-bibliography) to follow the correct format.
- (Optional) create a metadata file describing your software and include it in your repository. We provide [a script](https://gist.github.com/arfon/478b2ed49e11f984d6fb) that automates the generation of this metadata.

## What should my paper contain?

```eval_rst
.. important:: Begin your paper with a summary of the high-level functionality of your software for a non-specialist reader. Avoid jargon in this section.
```

JOSS welcomes submissions from broadly diverse research areas. For this reason, we require that authors include in the paper some sentences that explain the software functionality and domain of use to a non-specialist reader. We also require that authors explain the research applications of the software. The paper should be between 250-1000 words. Authors submitting papers significantly longer than 1000 words may be asked to reduce the length of their paper.

Your paper should include:

- A list of the authors of the software and their affiliations, using the correct format (see the example below).
- A summary describing the high-level functionality and purpose of the software for a diverse, *non-specialist audience*.
- A *Statement of Need* section that clearly illustrates the research purpose of the software.
- A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.
- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.
- Acknowledgement of any financial support.

As this short list shows, JOSS papers are only expected to contain a limited set of metadata (see example below), a Statement of Need, Summary, Acknowledgements, and References sections. You can look at an [example accepted paper](http://bit.ly/2x22gxT). Given this format, a "full length" paper is not permitted, and software documentation such as API (Application Programming Interface) functionality should not be in the paper and instead should be outlined in the software documentation.

```eval_rst
.. important:: Your paper will be reviewed by two or more reviewers in a public GitHub issue. Take a look at the `review checklist <review_checklist.html>`_ and  `review criteria <review_criteria.html>`_ to better understand how your submission will be reviewed.
```

## Example paper and bibliography

This example `paper.md` is adapted from _Gala: A Python package for galactic dynamics_ by Adrian M. Price-Whelan [http://doi.org/10.21105/joss.00388](http://doi.org/10.21105/joss.00388):

```
---
title: 'Gala: A Python package for galactic dynamics'
tags:
  - Python
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Adrian M. Price-Whelan^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0000-0000-0000
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID^[co-first author] # note this makes a footnote saying 'co-first author'
    affiliation: 2
  - name: Author with no affiliation^[corresponding author]
    affiliation: 3
affiliations:
 - name: Lyman Spitzer, Jr. Fellow, Princeton University
   index: 1
 - name: Institution Name
   index: 2
 - name: Independent Researcher
   index: 3
date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References

```

Example `paper.bib` file:

```
@article{Pearson:2017,
  	url = {http://adsabs.harvard.edu/abs/2017arXiv170304627P},
  	Archiveprefix = {arXiv},
  	Author = {{Pearson}, S. and {Price-Whelan}, A.~M. and {Johnston}, K.~V.},
  	Eprint = {1703.04627},
  	Journal = {ArXiv e-prints},
  	Keywords = {Astrophysics - Astrophysics of Galaxies},
  	Month = mar,
  	Title = {{Gaps in Globular Cluster Streams: Pal 5 and the Galactic Bar}},
  	Year = 2017
}

@book{Binney:2008,
  	url = {http://adsabs.harvard.edu/abs/2008gady.book.....B},
  	Author = {{Binney}, J. and {Tremaine}, S.},
  	Booktitle = {Galactic Dynamics: Second Edition, by James Binney and Scott Tremaine.~ISBN 978-0-691-13026-2 (HB).~Published by Princeton University Press, Princeton, NJ USA, 2008.},
  	Publisher = {Princeton University Press},
  	Title = {{Galactic Dynamics: Second Edition}},
  	Year = 2008
}

@article{gaia,
    author = {{Gaia Collaboration}},
    title = "{The Gaia mission}",
    journal = {Astronomy and Astrophysics},
    archivePrefix = "arXiv",
    eprint = {1609.04153},
    primaryClass = "astro-ph.IM",
    keywords = {space vehicles: instruments, Galaxy: structure, astrometry, parallaxes, proper motions, telescopes},
    year = 2016,
    month = nov,
    volume = 595,
    doi = {10.1051/0004-6361/201629272},
    url = {http://adsabs.harvard.edu/abs/2016A%26A...595A...1G},
}

@article{astropy,
    author = {{Astropy Collaboration}},
    title = "{Astropy: A community Python package for astronomy}",
    journal = {Astronomy and Astrophysics},
    archivePrefix = "arXiv",
    eprint = {1307.6212},
    primaryClass = "astro-ph.IM",
    keywords = {methods: data analysis, methods: miscellaneous, virtual observatory tools},
    year = 2013,
    month = oct,
    volume = 558,
    doi = {10.1051/0004-6361/201322068},
    url = {http://adsabs.harvard.edu/abs/2013A%26A...558A..33A}
}

@misc{fidgit,
  author = {A. M. Smith and K. Thaney and M. Hahnel},
  title = {Fidgit: An ungodly union of GitHub and Figshare},
  year = {2020},
  publisher = {GitHub},
  journal = {GitHub repository},
  url = {https://github.com/arfon/fidgit}
}
```

Note that the paper ends with a References heading, and the references are built automatically from the content in the `.bib` file. You should enter in-text citations in the paper body following correct [Markdown citation syntax](https://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html#citation_syntax).  Also note that the references include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.

## Checking that your paper compiles

JOSS uses Pandoc to compile papers from their Markdown form into a PDF. There are a few different ways you can test that your paper is going to compile properly for JOSS:

### JOSS paper preview service

Visit [https://whedon.theoj.org](https://whedon.theoj.org) and enter your repository address (and custom branch if you're using one). Note that your repository must be world-readable (i.e., it cannot require a login to access).

<img width="1348" alt="Screen Shot 2020-11-23 at 12 08 58 PM" src="https://user-images.githubusercontent.com/4483/99960475-b4f7be00-2d84-11eb-83bd-7784e9e23913.png">

### GitHub Action

If you're using GitHub for your repository, you can use the [Open Journals GitHub Action](https://github.com/marketplace/actions/open-journals-pdf-generator) to automatically compile your paper each time you update your repository.

The PDF is available via the Actions tab in your project and click on the latest workflow run. The zip archive file (including the `paper.pdf`) is listed in the run's Artifacts section.

### Docker

If you have Docker installed on your local machine, you can use the same Docker Image to compile a draft of your paper locally. In the example below, the `paper.md` file is in the `paper` directory. Upon successful execution of the command, the `paper.pdf` file will be created in the same location as the `paper.md` file:

```text
docker run --rm \
    --volume $PWD/paper:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/paperdraft
```

## Submitting your paper

Submission is as simple as:

- Filling in the [short submission form](http://joss.theoj.org/papers/new)
- Waiting for the managing editor to start a pre-review issue over in the JOSS reviews repository: https://github.com/openjournals/joss-reviews

## No submission fees

There are no fees for submitting or publishing in JOSS. You can read more about our [cost and sustainability model](http://joss.theoj.org/about#costs).

## Preprint Policy

Authors are welcome to submit their papers to a preprint server ([arXiv](https://arxiv.org/), [bioRxiv](https://www.biorxiv.org/), [SocArXiv](https://socopen.org/), [PsyArXiv](https://psyarxiv.com/) etc.) at any point before, during, or after the submission and review process.

Submission to a preprint server is _not_ considered a previous publication.

## Authorship

Purely financial (such as being named on an award) and organizational (such as general supervision of a research group) contributions are not considered sufficient for co-authorship of JOSS submissions, but active project direction and other forms of non-code contributions are. The authors themselves assume responsibility for deciding who should be credited with co-authorship, and co-authors must always agree to be listed. In addition, co-authors agree to be accountable for all aspects of the work, and to notify JOSS if any retraction or correction of mistakes are needed after publication.

## Submissions using proprietary languages/development environments

We strongly prefer software that doesn't rely upon proprietary (paid for) development environments/programming languages. However, provided _your submission meets our requirements_ (including having a valid open source license) then we will consider your submission for review. Should your submission be accepted for review, we may ask you, the submitting author, to help us find reviewers who already have the required development environment installed.

## The review process

After submission:

- An Associate Editor-in-Chief will carry out an initial check of your submission, and proceed to assign a handling editor.
- The handling editor will assign two or more JOSS reviewers, and the review will be carried out in the [JOSS reviews repository](https://github.com/openjournals/joss-reviews).
- Authors will respond to reviewer-raised issues (if any are raised) on the submission repository's issue tracker. Reviewer and editor contributions, like any other contributions, should be acknowledged in the repository.
- Upon successful completion of the review, authors will make a tagged release of the software, and deposit a copy of the repository with a data-archiving service such as [Zenodo](https://zenodo.org/) or [figshare](https://figshare.com/), get a DOI for the archive, and update the review issue thread with the version number and DOI.
- After we assign a DOI for your accepted JOSS paper, its metadata is deposited with CrossRef and listed on the JOSS website.
- The review issue will be closed, and an automatic tweet from [@JOSS_TheOJ](https://twitter.com/JOSS_TheOJ) will announce it!

If you want to learn more details about the review process, take a look at the [reviewer guidelines](reviewer_guidelines.html).

## Confidential requests

Please write admin@theoj.org with confidential matters such as retraction requests, report of misconduct, and retroactive author name changes.

In case of a name change, the DOI will be unchanged and the paper will be updated without publishing a correction notice or notifying co-authors.

JOSS will also update CrossRef metadata.
Review checklist
===============

JOSS reviews are checklist-driven. That is, there is a checklist for each JOSS reviewer to work through when completing their review. A JOSS review is generally considered incomplete until the reviewer has checked off all of their checkboxes.

Below is an example of the review checklist for the [Yellowbrick JOSS submission](https://github.com/openjournals/joss-reviews/issues/1075).

```eval_rst
.. important:: Note this section of our documentation only describes the JOSS review checklist. Authors and reviewers should consult the `review criteria <review_criteria.html>`_ to better understand how these checklist items should be interpreted.
```
### Conflict of interest

- I confirm that I have read the [JOSS conflict of interest policy](reviewer_guidelines.html#joss-conflict-of-interest-policy) and that: I have no COIs with reviewing this work or that any perceived COIs have been waived by JOSS for the purpose of this review.

### Code of Conduct

- I confirm that I read and will adhere to the [JOSS code of conduct](https://joss.theoj.org/about#code_of_conduct).

### General checks

- **Repository:** Is the source code for this software available at the <a target="_blank" href="https://github.com/DistrictDataLabs/yellowbrick">repository url</a>?
- **License:** Does the repository contain a plain-text LICENSE file with the contents of an [OSI approved](https://opensource.org/licenses/alphabetical) software license?
- **Contribution and authorship:** Has the submitting author made major contributions to the software? Does the full list of paper authors seem appropriate and complete?

### Functionality

- **Installation:** Does installation proceed as outlined in the documentation?
- **Functionality:** Have the functional claims of the software been confirmed?
- **Performance:** If there are any performance claims of the software, have they been confirmed? (If there are no claims, please check off this item.)

### Documentation

- **A statement of need:** Do the authors clearly state what problems the software is designed to solve and who the target audience is?
- **Installation instructions:** Is there a clearly-stated list of dependencies? Ideally these should be handled with an automated package management solution.
- **Example usage:** Do the authors include examples of how to use the software (ideally to solve real-world analysis problems).
- **Functionality documentation:** Is the core functionality of the software documented to a satisfactory level (e.g., API method documentation)?
- **Automated tests:** Are there automated tests or manual steps described so that the functionality of the software can be verified?
- **Community guidelines:** Are there clear guidelines for third parties wishing to 1) Contribute to the software 2) Report issues or problems with the software 3) Seek support

### Software paper

- **Summary:** Has a clear description of the high-level functionality and purpose of the software for a diverse, non-specialist audience been provided?
- **A statement of need:** Does the paper have a section titled 'Statement of Need' that clearly states what problems the software is designed to solve and who the target audience is?
- **State of the field:** Do the authors describe how this software compares to other commonly-used packages?
- **Quality of writing:** Is the paper well written (i.e., it does not require editing for structure, language, or writing quality)?
- **References:** Is the list of references complete, and is everything cited appropriately that should be cited (e.g., papers, datasets, software)? Do references in the text use the proper [citation syntax]( https://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html#citation_syntax)?
# Installing the JOSS application

Any Open Journal (JOSS, JOSE, etc.) can be considered in three parts:

1. The website
2. The Whedon gem
3. A thin API wrapper around the Whedon gem to interface with GitHub

For JOSS, these correspond to the:

1. [JOSS](https://github.com/openjournals/joss),
2. [Whedon](https://github.com/openjournals/whedon), and
3. [Whedon-API](https://github.com/openjournals/whedon-api)

code bases.

## Setting up a local development environment

All Open Journals are coded in Ruby,
with the website and Whedon-API developed as
[Ruby on Rails](https://rubyonrails.org/inst) applications.

If you'd like to develop these locally,
you'll need a working Ruby on Rails development environment.
For more information, please see
[this official guide](https://guides.rubyonrails.org/getting_started.html#creating-a-new-rails-project-installing-rails).

## Deploying your JOSS application

To deploy JOSS, you'll need a [Heroku account](https://signup.heroku.com/).
We also recommend you gather the following information
before deploying your application to Heroku:

1. A [public ORCID API](https://members.orcid.org/api/about-public-api) key
1. A GitHub [Personal Access Token](https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/creating-a-personal-access-token) for the automated account that users will interact with (e.g., `@Whedon`, `@RoboNeuro`). In order to be able to send invitations to reviewers and collaborators, the automated GitHub account must be an admin of the organization the reviews take place at. And the Personal Access Token should include the `admin:org` scope.
1. An email address registered on a domain you control (i.e., not `gmail` or a related service)

```eval_rst
.. warning::
    Do not put these secrets directly into your code base!
    It is important that these keys are not under version control.

    There are different ways to make sure your application has access to these keys,
    depending on whether your code is being developed locally or on Heroku.
    Locally, you can store these locally in a .env file.
    The .gitignore in JOSS is already set to ignore this file type.

    On Heroku, they will be config variables that you can set either with the Heroku CLI or directly on your application's dashboard.
    See `this guide from Heroku <https://devcenter.heroku.com/articles/config-vars#managing-config-vars>`_ for more information.
```

Assuming you [have forked the JOSS GitHub repository](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo)
to your account,
you can [configure Heroku as a git remote](https://devcenter.heroku.com/articles/git#prerequisites-install-git-and-the-heroku-cli) for your code.
This makes it easy to keep your Heroku deployment in sync with your local development copy.

On the JOSS Heroku deployment, you'll need to provision several [add-ons](https://elements.heroku.com/addons).
Specifically, you'll need:

1. [Elasticsearch add-on](https://elements.heroku.com/addons/bonsai)
1. [PostgreSQL add-on](https://elements.heroku.com/addons/heroku-postgresql)
1. [Scheduler add-on](https://devcenter.heroku.com/articles/scheduler)

For the scheduler add-on, you'll need to designate which tasks it should run and when.
These can be found in the `lib/tasks` folder, and involve things such as sending out weekly reminder emails to editors.
Each task should be scheduled as a separate job; for example, `rake send_weekly_emails`.

You can also optionally configure the following add-ons (or simply set their secret keys in your config variables):

1. [SendGrid add-on](https://elements.heroku.com/addons/sendgrid) for sending emails
1. [Honeybadger add-on](https://elements.heroku.com/addons/honeybadger) for error reporting

Once you've pushed your application to Heroku and provisioned the appropriate add-ons,
you're ready to update your config with the appropriate secrets.
For a list of the expected secret key names, see the `app.json` file.

```eval_rst
.. warning::
    One "gotcha" when provisioning the Bonsai add-on is that it may only set the BONSAI_URL variable.
    Make sure that there is also an ELASTICSEARCH_URL which is set to the same address.
```

We will not cover Portico, as this requires that your application is a part of the `openjournals` organization.
If you do not already have access to these keys, you can simply ignore them for now.

```eval_rst
.. note::
    One secret key we have not covered thus far is WHEDON_SECRET.
    This is because it is not one that you obtain from a provide,
    but a secret key that you set yourself.
    We recommend using something like a random SHA1 string.

    It is important to remember this key,
    as you will need it when deploying your Whedon-API application.
```

After pushing your application to Heroku, provisioning the appropriate add-ons,
and confirming that your config variables are set correctly,
you should make sure that your username is registered as an admin on the application.

You can do this on a local Rails console, by logging in and setting the boolean field 'admin'
on your user ID to True.
If you'd prefer to do this on the Heroku deployment, make sure you've logged into the application.
Then you can directly modify this attribute in the deployments Postgres database using SQL.
For more information on accessing your application's Postgres database,
see [the official docs](https://devcenter.heroku.com/articles/heroku-postgresql#pg-psql).

## Making modifications to launch your own site

Some times you may not want to launch an exact copy of JOSS, but a modified version.
This can be especially useful if you're planning to spin up your own platform based on the
Open Journals framework.
[NeuroLibre](https://neurolibre.herokuapp.com) is one such example use-case.

### Modifying your site configuration

In this case, there are several important variables to be aware of and modify.
Most of these are accessible in the `config` folder.

First, there are three files which provide settings for your Rails application in different development contexts:

1. `settings-test.yml`
1. `settings-development.yml`
1. `settings-production.yml`

These each contain site-specific variables that should be modified if you are building off of the Open Journals framework.

Next, you'll need to modify the `repository.yml` file.
This file lists the GitHub repository where you expect papers to be published,
as well as the editor team ID.
For your GitHub organization, make sure you have created and populated a team called `editors`.
Then, you can check its ID number as detailed in [this guide](https://fabian-kostadinov.github.io/2015/01/16/how-to-find-a-github-team-id).
In `config` you should also modify the `orcid.yml` file to list your site as the production site.

Finally, you'll need to set up a [GitHub webhook](https://docs.github.com/en/free-pro-team@latest/developers/webhooks-and-events/about-webhooks) for reviews repository.
This should be a repository that you have write access to,
where you expect most of the reviewer-author interaction to occur.
For JOSS, this corresponds to the `openjournals/joss-reviews` GitHub repository.

In this GitHub repository's settings,
you can add a new webhook with the following configuration:

- Set the `Payload` URL to the `/dispatch` hook for your Heroku application URL.
  For example, https://neurolibre.herokuapp.com/dispatch
- Set the `Content type` to `application/json`
- Set the secret to a high-entropy, random string as detailed in the [GitHub docs](https://docs.github.com/en/free-pro-team@latest/developers/webhooks-and-events/securing-your-webhooks#setting-your-secret-token)
- Set the webhook to deliver `all events`

#### Updating your application database

Optionally, you can edit `seeds.rb`, a file in the `db` folder.
"DB" is short for "database," and this file _seeds_ the database with information about your editorial team.
You can edit the file `seeds.rb` to remove any individuals who are not editors in your organization.
This can be especially useful if you will be creating the database multiple times during development;
for example, if you add in testing information that you'd later like to remove.

You can reinitialize the database from your Heroku CLI using the following commands:

```bash
heroku pg:reset DATABASE_URL
heroku run rails db:schema:load
heroku run rails db:seed
heroku run rails searchkick:reindex:all
```

### Modifying your site contents

You can modify site content by updating files in the `app` and `docs` folders.
For example, in `app/views/notifications` you can change the text for any emails that will be sent by your application.

Note that files which end in `.html.erb` are treated as HTML files, and typical HTML formatting applies.
You can set the HTML styling by modifying the Sass files for your application,
located in `app/assets/stylesheets`.

There are currently a few hard-coded variables in the application which you will also need to update.
Note that these are mostly under `lib/tasks`.
For example, in `stats.rake`, the reviewer sheet ID is hard-coded on line 37.
You should update this to point to your own spreadsheet where you maintain a list of eligible reviewers.

In the same folder, `utils.rake` is currently hard-coded to alternate assignments of editor-in-chief based on weeks.
You should modify this to either set a single editor-in-chief,
or design your own scheme of alternating between members of your editorial board.

## Deploying your Whedon-API Application

Whedon-API can also be deployed on Heroku.
Note that &mdash; for full functionality &mdash; Whedon-API must be deployed on [Hobby dynos](https://devcenter.heroku.com/articles/dyno-types), rather than free dynos.
Hobby dynos allow the Whedon-API application to run continuously, without falling asleep after 30 minutes of inactivity;
this means that Whedon-API can respond to activity on GitHub at any time.
Whedon-API specifically requires two Hobby dynos: one for the `web` service and one for the `worker` service.

On the Whedon-API Heroku deployment, you'll need to provision several [add-ons](https://elements.heroku.com/addons).
Specifically, you'll need:

1. [Cloudinary add-on](https://elements.heroku.com/addons/cloudinary)
1. [Scheduler add-on](https://devcenter.heroku.com/articles/scheduler)
1. [Redis To Go add-on](https://elements.heroku.com/addons/redistogo)

For the scheduler add-on, you'll need to designate which tasks it should run and when.
The only task that needs to be scheduled is the `restart.sh` script,
which should be set to execute every hour.

```eval_rst
.. warn:
    Cloudinary `does not allow free accounts to serve PDFs <https://cloudinary.com/blog/uploading_managing_and_delivering_pdfs#delivering_pdf_files>`_ by default.
    This will prevent your application from offering a paper preview service, as in https://whedon.theoj.org
    To have this restriction lifted, you will need to `contact Cloudinary customer support <https://support.cloudinary.com/hc/en-us/requests/new>`_ directly.
```

As before, once you've pushed your application to Heroku and provisioned the appropriate add-ons,
you're ready to update your config with the appropriate secrets.
For a list of the expected secret key names, see the `app.json` file.
Many of these will be re-used from deploying your JOSS application.

Specifically, the `GH_TOKEN` should be the same personal access token as before.
The `JOSS_API_KEY` should match the `WHEDON_SECRET` key that you created in your JOSS deployment.

You'll also need to provide a `HEROKU_APP_NAME`, `HEROKU_CLI_TOKEN`, and `HEROKU_CLI_USER` that the `restart.sh` script can use when executing.
You can find these directly from the heroku-cli as detailed in [their documentation](https://devcenter.heroku.com/articles/authentication).

## Modifying your Whedon-API deployment

Some times you may not want to launch an exact copy of the Whedon-API, but a modified version.
This can be especially useful if you're planning to spin up your own platform based on the
Open Journals framework.
[RoboNeuro](https://github.com/roboneuro) is one such example use-case.

### Modifying your Whedon-API configuration

Similar to the JOSS deployment described above,
the Whedon-API configuration is controlled through a series of YAML files included in the `config/` folder.
Each of these files provide relevant configuration for a different development context.
Specifically, two files are defined:

1. `settings-test.yml`
1. `settings-production.yml`

which can be used to define testing and production environment variables, respectively.
Much of the information [previously defined for your JOSS deployment](#modifying-your-site-configuration) will carry over,
including the editor team ID.

Finally, you'll need to set up a [GitHub webhook](https://docs.github.com/en/free-pro-team@latest/developers/webhooks-and-events/about-webhooks) for your reviews repository.
As a reminder, this should be a repository that you have write access to.
For JOSS, this corresponds to the `openjournals/joss-reviews` GitHub repository.
**This is in addition to the webhook you previously created for the JOSS deployment,
although it points to the same repository.**

In this GitHub repository's settings,
you can add a new webhook with the following configuration:

- Set the `Payload` URL to the `/dispatch` hook for your Heroku application URL.
  For example, https://roboneuro.herokuapp.com/dispatch
- Set the `Content type` to `application/json`
- Set the secret to a high-entropy, random string as detailed in the [GitHub docs](https://docs.github.com/en/free-pro-team@latest/developers/webhooks-and-events/securing-your-webhooks#setting-your-secret-token)
- Set the webhook to deliver `all events`
Interacting with Whedon
========================

Whedon or `@whedon` on GitHub, is our editorial bot that interacts with authors, reviewers, and editors on JOSS reviews.

`@whedon` can do a bunch of different things. If you want to ask `@whedon` what it can do, simply type the following in a JOSS `review` or `pre-review` issue:

```text
# List all of Whedon's capabilities
@whedon commands

# Assign a GitHub user as the sole reviewer of this submission
@whedon assign @username as reviewer

# Add a GitHub user to the reviewers of this submission
@whedon add @username as reviewer

# Re-invite a reviewer (if they can't update checklists)
@whedon re-invite @username as reviewer

# Remove a GitHub user from the reviewers of this submission
@whedon remove @username as reviewer

# List of editor GitHub usernames
@whedon list editors

# List of reviewers together with programming language preferences and domain expertise
@whedon list reviewers

# Change editorial assignment
@whedon assign @username as editor

# Set the software archive DOI at the top of the issue e.g.
@whedon set 10.0000/zenodo.00000 as archive

# Set the software version at the top of the issue e.g.
@whedon set v1.0.1 as version

# Open the review issue
@whedon start review

EDITORIAL TASKS

# All commands can be run on a non-default branch, to do this pass a custom 
# branch name by following the command with `from branch custom-branch-name`.
# For example:

# Compile the paper
@whedon generate pdf

# Compile the paper from alternative branch
@whedon generate pdf from branch custom-branch-name

# Remind an author or reviewer to return to a review after a
# certain period of time (supported units days and weeks)
@whedon remind @reviewer in 2 weeks

# Ask Whedon to do a  dry run of accepting the paper and depositing with Crossref
@whedon recommend-accept

# Ask Whedon to check the references for missing DOIs
@whedon check references

# Ask Whedon to check repository statistics for the submitted software, for license, and
# for Statement of Need section in paper
@whedon check repository

EiC TASKS

# Flag submission for editoral review, due to size or question about being research software
@whedon query scope

# Invite an editor to edit a submission (sending them an email)
@whedon invite @editor as editor

# Reject a paper
@whedon reject

# Withdraw a paper
@whedon withdraw

# Ask Whedon to actually accept the paper and deposit with Crossref
# (supports custom branches too)
@whedon accept deposit=true
```

## Author commands

A subset of the Whedon commands are available to authors (and reviewers):

### Compiling papers

When a `pre-review` or `review` issue is opened, `@whedon` will try to compile the JOSS paper by looking for a `paper.md` file in the repository specified when the paper was submitted.

If it can't find the `paper.md` file it will say as much in the review issue. If it can't compile the paper (i.e. there's some kind of Pandoc error), it will try and report that error back in the thread too.

```eval_rst
.. note:: If you want to see what command ``@whedon`` is running when compiling the JOSS paper, take a look at the code `here <https://github.com/openjournals/whedon/blob/195e6d124d0fbd5346b87659e695325df9a18334/lib/whedon/processor.rb#L109-L132>`_.
```

Anyone can ask `@whedon` to compile the paper again (e.g. after a change has been made). To do this simply comment on the review thread as follows:

```text
@whedon generate pdf
```

#### Compiling papers from a non-default branch

By default, Whedon will look for papers in the default git branch. If you want to compile a paper from a non-default branch, this can be done as follows:

```text
@whedon generate pdf from branch custom-branch-name
```

### Finding reviewers

Sometimes submitting authors suggest people the think might be appropriate to review their submission. If you want the link to the current list of JOSS reviewers, type the following in the review thread:

```text
@whedon list reviewers
```

## Editorial commands

Most of `@whedon`'s functionality can only be used by the journal editors.

### Assigning an editor

Editors can either assign themselves or other editors as the editor of a submission as follows:

```text
@whedon assign @editorname as editor
```

### Inviting an editor

Whedon can be used by EiCs to send email invites to registered editors as follows:

```text
@whedon invite @editorname as editor
```

This will send an automated email to the editor with a link to the GitHub `pre-review` issue.

### Adding and removing reviewers

Reviewers should be assigned by using the following commands:

```text
# Assign a GitHub user as the sole reviewer of this submission
@whedon assign @username as reviewer

# Add a GitHub user to the reviewers of this submission
@whedon add @username as reviewer

# Remove a GitHub user from the reviewers of this submission
@whedon remove @username as reviewer
```

```eval_rst
.. note:: The ``assign`` command clobbers all reviewer assignments. If you want to add an additional reviewer use the ``add`` command.
```

### Starting the review

Once the reviewer(s) and editor have been assigned in the `pre-review` issue, the editor starts the review with:

```text
@whedon start review
```

```eval_rst
.. important:: If a reviewer recants their commitment or is unresponsive, editors can remove them with the command ``@whedon remove @username as reviewer``. You can also add new reviewers in the ``REVIEW`` issue, but in this case, you need to manually add a review checklist for them by editing the issue body.
```

### Reminding reviewers and authors

Whedon can reminders authors and reviewers after a specified amount of time to return to the review issue. Reminders can only be set by editors, and only for REVIEW issues. For example:

```text
# Remind the reviewer in two weeks to return to the review
@whedon remind @reviewer in two weeks
```

```text
# Remind the reviewer in five days to return to the review
@whedon remind @reviewer in five days
```

```text
# Remind the author in two weeks to return to the review
@whedon remind @author in two weeks
```

```eval_rst
.. note:: Most units of times are understood by Whedon e.g. `hour/hours/day/days/week/weeks`.
```

```eval_rst
.. important:: For reviewers, the reminder will only be triggered if the reviewer's review is outstanding (i.e. outstanding checkboxes).
```

### Setting the software archive

When a submission is accepted, we ask that the authors create an archive (on [Zenodo](https://zenodo.org/), [fig**share**](https://figshare.com/), or other) and post the archive DOI in the `REVIEW` issue. The editor should add the `accepted` label on the issue and ask `@whedon` to add the archive to the issue as follows:

```text
@whedon set 10.0000/zenodo.00000 as archive
```

### Changing the software version

Sometimes the version of the software changes as a consequence of the review process. To update the version of the software do the following:

```text
@whedon set v1.0.1 as version
```

## Accepting a paper (dry run)

Whedon can accept a paper from the review issue. This includes generating the final paper PDF, Crossref metedata, and depositing this metadata with the Crossref API.

JOSS topic editors can ask for the final proofs to be created by Whedon with the following command:

```text
@whedon recommend-accept
```

On issuing this command, Whedon will also check the references of the paper for any missing DOIs. This command can be triggered separately:

### Check references

```text
@whedon check references
```

```eval_rst
.. note:: Whedon can verify that DOIs resolve, but cannot verify that the DOI associated with a paper is actually correct. In addition, DOI suggestions from Whedon are just that - i.e. they may not be correct.
```

## Accepting a paper (for real)

If everything looks good with the draft proofs from the `@whedon accept` command, JOSS editors-in-chief can take the additional step of actually accepting the JOSS paper with the following command:

```text
@whedon accept deposit=true
```

```eval_rst
.. note:: This command is only available to the JOSS editor-in-chief, or associate editors-in-chief.
```
Editorial Guide
===============

The Journal of Open Source Software (JOSS) conducts all peer review and editorial processes in the open, on the GitHub issue tracker.

JOSS editors manage the review workflow with the help of our bot, `@whedon`. The bot is summoned with commands typed directly on the GitHub review issues. For a list of commands, type: `@whedon commands`.

```eval_rst
.. note:: To learn more about ``@whedon``'s functionalities, take a look at our `dedicated guide <whedon.html>`_.
```

## Pre-review

Once a submission comes in, it will be in the queue for a quick check by the Editor-in-chief (EiC). From there, it moves to a `PRE-REVIEW` issue, where the EiC will assign a handling editor, and the author can suggest reviewers. Initial direction to the authors for improving the paper can already happen here, especially if the paper lacks some requested sections.

```eval_rst
.. important:: If the paper is out-of-scope for JOSS, editors assess this and notify the author in the ``PRE-REVIEW`` issue.
```

Editors can flag submissions of questionable scope using the command `@whedon query scope`.

The EiC assigns an editor (or a volunteering editor self-assigns) with the command `@whedon assign @username as editor` in a comment.

```eval_rst
.. note:: If a paper is submitted without a recommended editor, it will show up in the weekly digest email under the category ‘Papers currently without an editor.’ Please review this weekly email and volunteer to edit papers that look to be in your domain. If you choose to be an editor in the issue thread type the command ``@whedon assign @yourhandle as editor`` or simply ``@whedon assign me as editor``
```

### How papers are assigned to editors

By default, unless an editor volunteers, the Associated Editor-in-chief (AEiC) on duty will attempt to assign an incoming paper to the most suitable handling editor. While AEiCs will make every effort to match a submission with the most appropriate editor, there are a number of situations where an AEiC may assign a paper to an editor that doesn't fit entirely within the editor's research domains:

- If there's no obvious fit to _any_ of the JOSS editors
- If the most suitable editor is already handling a large number of papers
- If the chosen editor has a lighter editorial load than other editors

In most cases, an AEiC will ask one or more editors to edit a submission (e.g. `@editor1, @editor 2 - would one of you be willing to edit this submission for JOSS`). If the editor doesn't respond within ~3 working days, the AEiC may assign the paper to the editor regardless.

Editors may also be invited to edit over email when an AEiC runs the  command `@whedon invite @editor1 as editor`.

### Finding reviewers

At this point, the handling editor's job is to identify reviewers who have sufficient expertise in the field of software and in the field of the submission. JOSS papers have to have a minimum of two reviewers per submission, except for papers that have previously been peer-reviewed via rOpenSci. In some cases, the editor also might want to formally add themself as one of the reviewers. If the editor feels particularly unsure of the submission, a third (or fourth) reviewer can be recruited.

To recruit reviewers, the handling editor can mention them in the `PRE-REVIEW` issue with their GitHub handle, ping them on Twitter, or email them. After expressing initial interest, candidate reviewers may need a longer explanation via email. See sample reviewer invitation email, below.

**Reviewer Considerations**

- It is rare that all reviewers have the expertise to cover all aspects of a submission (e.g., knows the language really well and knows the scientific discipline well). As such, a good practice is to try and make sure that between the two or three reviewers, all aspects of the submission are covered.
- Selection and assignment of reviewers should adhere to the [JOSS COI policy](https://joss.theoj.org/about#ethics).

**Potential ways to find reviewers**

Finding reviewers can be challenging, especially if a submission is outside of your immediate area of expertise. Some strategies you can use to identify potential candidates:

- Search the [reviewer spreadsheet](https://bit.ly/joss-reviewers) of volunteer reviewers.
  - When using this spreadsheet, pay attention to the number of reviews this individual is already doing to avoid overloading them.
  - It can be helpful to use the "Data > Filter Views" capability to temporarily filter the table view to include only people with language or domain expertise matching the paper.
- Ask the author(s): You are free to ask the submitting author to suggest possible reviewers by using the [reviewer spreadsheet](https://bit.ly/joss-reviewers) and also people from their professional network. In this situation, the editor still needs to verify that their suggestions are appropriate.
- Use your professional network: You're welcome to invite people you know of who might be able to give a good review.
- Search Google and GitHub for related work, and write to the authors of that related work.
  - You might like to try [this tool](https://github.com/dfm/joss-reviewer) from @dfm.
- Ask on social networks: Sometimes asking on Twitter for reviewers can identify good candidates.
- Check the work being referenced in the submission:
  - Authors of software that is being built on might be interested in reviewing the submission.
  - Users of the the software that is being submission be interested in reviewing the submission
- Avoid asking JOSS editors to review: If at all possible, avoid asking JOSS editors to review as they are generally very busy editing their own papers.

Once a reviewer accepts, the handling editor runs the command `@whedon assign @username as reviewer` in the `PRE-REVIEW` issue. Add more reviewers with the command `@whedon add @username as reviewer`.
Under the uncommon circumstance that a review must be started before all reviewers have been identified (e.g., if finding a second reviewer is taking a long time and the first reviewer wants to get started), an editor may elect to start the review and add the remaining reviewers later. To accomplish this, the editor will need to hand-edit the review checklist to create space for the reviewers added after the review issues is created.

```eval_rst
.. note:: The ``assign`` command clobbers all reviewer assignments. If you want to add an additional reviewer use the ``add`` command.
```

### Starting the review

Next, run the command `@whedon start review`. If you haven't assigned an editor and reviewer, this command will fail and `@whedon` will tell you this. This will open the `REVIEW` issue, with prepared review checklists for each reviewer, and instructions. The editor should close the `PRE-REVIEW` issue, at this point, and move the conversation to the separate `REVIEW` issue.

## Review

The `REVIEW` issue contains some instructions, and reviewer checklists. The reviewer(s) should check off items of the checklist one-by-one, until done. In the meantime, reviewers can engage the authors freely in a conversation aimed at improving the paper.

```eval_rst
.. note:: When a new review is started, Whedon invites the reviewers to the GitHub repository. The reviewers must accept the invitation to the GitHub repository in order for Whedon to be able to assign the review issue to them, and reviewers are unable to check off checklist items until they have accepted the repository invitation.
```

If a reviewer recants their commitment or is unresponsive, editors can remove them with the command `@whedon remove @username as reviewer`. You can also add new reviewers in the `REVIEW` issue, but in this case, you need to manually add a review checklist for them by editing the issue body.

Comments in the `REVIEW` issue should be kept brief, as much as possible, with more lengthy suggestions or requests posted as separate issues, directly in the submission repository. A link-back to those issues in the `REVIEW` is helpful.

When the reviewers are satisfied with the improvements, we ask that they confirm their recommendation to accept the submission.

### Adding a new reviewer once the review has started

Sometimes you'll need to add a new reviewer once the main review (i.e. post pre-review) is already underway. In this situation you should do the following:

- In the review thread, do `@whedon add @newreviewer as reviewer`.
- Manually edit the first message in the thread to add a checklist for @newreviewer.

In the future, this will be more automated but for now, there's some manual work required.

## After reviewers recommend acceptance

When a submission is ready to be accepted, we ask that the authors issue a new tagged release of the software (if changed), and archive it (on [Zenodo](https://zenodo.org/), [fig**share**](https://figshare.com/), or other). The authors then post the version number and archive DOI in the `REVIEW` issue. The handling editor executes the pre-publication steps, and pings the EiCs for final processing.

Pre-publication steps:
- Get a new proof with the `@whedon generate pdf` command.
- Download the proof, check all references have DOIs, follow the links and check the references.
  - Whedon can help check references with the command `@whedon check references`
- Proof-read the paper and ask authors to fix any remaining typos, badly formed citations, awkward wording, etc..
- Ask the author to make a tagged release and archive, and report the version number and archive DOI in the review thread.
- Check the archive deposit has the correct metadata (title and author list), and request the author edit it if it doesn’t match the paper.
- Run `@whedon set <doi> as archive`.
- Run `@whedon set <v1.x.x> as version` if the version was updated.
- Run `@whedon recommend-accept` to generate the final proofs, which has Whedon notify the `@openjournals/joss-eics` team that the paper is ready for final processing.

At this point, the EiC/AEiC will take over to make final checks and publish the paper.

It’s also a good idea to ask the authors to check the proof. We’ve had a few papers request a post-publication change of author list, for example—this requires a manual download/compile/deposit cycle and should be a rare event.

## Handling of papers published together with AAS publishing

JOSS is collaborating with [AAS publishing](https://journals.aas.org/) to offer software review for some of the papers submitted to their journals. A detailed overview of the motivations/background is available in the [announcement blog post](https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing), here we document the additional editorial steps that are necessary for JOSS to follow:

**Before/during review**

- If the paper is a joint publication, make sure you apply the [`AAS`](https://github.com/openjournals/joss-reviews/issues?utf8=%E2%9C%93&q=is%3Aissue+label%3AAAS+) label to both the `pre-review` and the `review` issues.
- Before moving the JOSS paper from `pre-review` to `review`, ensure that you (the JOSS editor) make the reviewers aware that JOSS will be receiving a small financial donation from AAS publishing for this review (e.g. [like this](https://github.com/openjournals/joss-reviews/issues/1852#issuecomment-553203738)).

**After the paper has been accepted by JOSS**

- Once the JOSS review is complete, ask the author for the status of their AAS publication, specifically if they have the AAS paper DOI yet.
- Once this is available, ask the author to add this information to their `paper.md` YAML header as documented in the [submission guidelines](submitting.html#what-should-my-paper-contain).

```
# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
```

- Pause the review (by applying the `paused` label) to await notification that the AAS paper is published.

**Once the AAS paper is published**

- Ask the EiC on rotation to publish the paper as normal (by tagging `@openjournals/joss-eics`).

## Processing of rOpenSci-reviewed and accepted submissions

If a paper has already been reviewed and accepted by rOpenSci, the streamlined JOSS review process is:

- Assign yourself as editor and reviewer
- Add a comment in the pre-review issue pointing to the rOpenSci review
- Add the rOpenSci label to the pre-review issue
- Start the review issue
- Add a comment in the review issue pointing to the rOpenSci review
- Compile the paper and check it looks ok
- Tick off all the review checkboxes
- Go to to the source code repo and grab the Zenodo DOI
- Accept and publish the paper

## Rejecting a paper

If you believe a submission should be rejected, for example, because it is out of scope for JOSS, then you should:

- Ask Whedon to flag the submission as potentially out of scope with the command `@whedon query scope`. This command adds the `query-scope` label to the issue.
- Mention to the author your reasons for flagging the submission as possibly out of scope, and give them an opportunity to defend their submission.
- The EiC on rotation will make a final determination of whether a submission is in scope, taking into account the feedback of other editors.

### Voting on papers flagged as potentially out of scope

Once per week, an email is sent to all JOSS editors with a summary of the papers that are currently flagged as potentially out of scope. Editors are asked to review these submissions and vote on the JOSS website if they have an opinion about a submission.

## Sample messages for authors and reviewers

### Sample email to potential reviewers

```
Dear Dr. Jekyll,

I found you following links from the page of The Super Project and/or on Twitter. This
message is to ask if you can help us out with a submission to JOSS (The Journal of Open
Source Software, https://joss.theoj.org), where I’m an editor.

JOSS publishes articles about open source research software. The submission I'd like you
to review is titled: "great software name here"

and the submission repository is at: https://github.com/< … >

JOSS is a free, open-source, community driven and developer-friendly online journal
(no publisher is seeking to raise revenue from the volunteer labor of researchers!).

The review process at JOSS is unique: it takes place in a GitHub issue, is open,
and author-reviewer-editor conversations are encouraged.

JOSS reviews involve downloading and installing the software, and inspecting the repository
and submitted paper for key elements. See https://joss.readthedocs.io/en/latest/review_criteria.html

Editors and reviewers post comments on the Review issue, and authors respond to the comments
and improve their submission until acceptance (or withdrawal, if they feel unable to
satisfy the review).

Would you be able to review this submission for JOSS? If not, can you recommend
someone from your team to help out?

Kind regards,

JOSS Editor.
```

### Query scope of submission

```
:wave: thanks for your submission to JOSS. From a quick inspection of this submission it's not entirely obvious that it meets our [submission criteria](https://joss.readthedocs.io/en/latest/submitting.html#submission-requirements). In particular, this item:

> - Your software should have an obvious research application

Could you confirm here that there _is_ a research application for this software (and explain what that application is)? The section [_'what should my paper contain'_](https://joss.readthedocs.io/en/latest/submitting.html#what-should-my-paper-contain) has some guidance for the sort of content we're looking to be present in the `paper.md`.

Many thanks!
```

### GitHub invite to potential reviewers

```
:wave: @reviewer1 & @reviewer2, would any of you be willing to review this submission for JOSS? We carry out our checklist-driven reviews here in GitHub issues and follow these guidelines: https://joss.readthedocs.io/en/latest/review_criteria.html
```

### Message to reviewers at the start of a review

```
👋🏼 @authorname @reviewer1 @reviewer2 this is the review thread for the paper. All of our communications will happen here from now on.

Both reviewers have checklists at the top of this thread with the JOSS requirements. As you go over the submission, please check any items that you feel have been satisfied. There are also links to the JOSS reviewer guidelines.

The JOSS review is different from most other journals. Our goal is to work with the authors to help them meet our criteria instead of merely passing judgment on the submission. As such, the reviewers are encouraged to submit issues and pull requests on the software repository. When doing so, please mention `openjournals/joss-reviews#REVIEW_NUMBER` so that a link is created to this thread (and I can keep an eye on what is happening). Please also feel free to comment and ask questions on this thread. In my experience, it is better to post comments/questions/suggestions as you come across them instead of waiting until you've reviewed the entire package.

We aim for reviews to be completed within about 2-4 weeks. Please let me know if any of you require some more time. We can also use Whedon (our bot) to set automatic reminders if you know you'll be away for a known period of time.

Please feel free to ping me (@editorname) if you have any questions/concerns.
```

### Message to authors at the end of a review

```
At this point could you:
- [ ] Make a tagged release of your software, and list the version tag of the archived version here.
- [ ] Archive the reviewed software in Zenodo or a similar service (e.g., figshare, an institutional repository)
- [ ] Check the archival deposit (e.g., in Zenodo) has the correct metadata. This includes the title (should match the paper title) and author list (make sure the list is correct and people who only made a small fix are not on it). You may also add the authors' ORCID.
- [ ] Please list the DOI of the archived version here.

I can then move forward with accepting the submission.
```

###

## Overview of editorial process

**Step 1: An author submits a paper.**

The author can choose to select an preferred editor based on the information available in our biographies. This can be changed later.

**Step 2: If you are selected as an editor you get @-mentioned in the pre-review issue.**

This doesn’t mean that you’re the editor, just that you’ve been suggested by the author.

**Step 3: Once you are the editor, find the link to the code repository in the `pre-review` issue**

**Step 4: The editor looks at the software submitted and checks to see if:**

- There’s a general description of the software
- The software is within scope as research software
- It has an OSI-approved license

**Step 5: The editor responds to the author saying that things look in line (or not) and will search for reviewer**

**Step 6: The editor finds >= 2 reviewers**

- Use the list of reviewers: type the command `@whedon list reviewers` or look at list of reviewers in a Google [spreadsheet](https://docs.google.com/spreadsheets/d/1PAPRJ63yq9aPC1COLjaQp8mHmEq3rZUzwUYxTulyu78/edit?usp=sharing)
- If people are in the review list, the editor can @-mention them on the issue to see if they will review: e.g. `@person1 @person2 can you review this submission for JOSS?`
- Or solicit reviewers outside the list. Send an email to people describing what JOSS is and asking if they would be interested in reviewing.
- If you ask the author to suggest potential reviewers, please be sure to tell the author not to @-tag their suggestions.

**Step 7: Editor tells Whedon to assign the reviewer to the paper**

- Use `@whedon assign @reviewer as reviewer`
- To add a second reviewer use `@whedon add @reviwer2 as reviewer`

```eval_rst
.. note:: The ``assign`` command clobbers all reviewer assignments. If you want to add an additional reviewer use the ``add`` command.
```

**Step 8: Create the actual review issue**

- Use `@whedon start review`
- An issue is created with the review checklist, one per reviewer, e.g. https://github.com/openjournals/joss-reviews/issues/717

**Step 9: Close the pre-review issue**

**Step 10: The actual JOSS review**

- The reviewer reviews the paper and has a conversation with the author. The editor lurks on this conversation and comes in if needed for questions (or CoC issues).
- The reviewer potentially asks for changes and the author makes changes. Everyone agrees it’s ready.

**Step 11: The editor pings the EiC team to get the paper published**

- Make a final check of the paper with `@whedon generate pdf` and that all references have DOIs (where appropriate) with `@whedon check references`.
- If everything looks good, ask the author to make a new release (if possible) of the software being reviewed and deposit a new archive the software with Zenodo/figshare. Update the review thread with this archive DOI: `@whedon set 10.5281/zenodo.xxxxxx` as archive.
- Finally, use `@whedon recommend-accept` on the review thread to ping the `@openjournals/joss-eics` team letting them know the paper is ready to be accepted.

**Step 12: Celebrate publication! Tweet! Thank reviewers! Say thank you on issue.**

## Visualization of editorial flow

![Editorial flow](images/JOSS-flowchart.png)

## Expectations on JOSS editors

### Responding to editorial assignments

As documented above, usually, papers will be assigned to you by one of the AEiCs. We ask that editors do their best to respond in a timely fashion (~ 3 working days) to invites to edit a new submission.

### Continued attention to assigned submissions

As an editor, part of your role is to ensure that submissions you're responsible for are progressing smoothly through the editorial process. This means that once or twice per week we ask that you check your GitHub notifications and/or your editorial dashboard (e.g. `http://joss.theoj.org/dashboard/youreditorname`) for updates to the papers you are handling.

**If reviews go stale**

Sometimes reviews go quiet, either because a reviewer has failed to complete their review or an author has been slow to respond to a reviewer's feedback. **As the editor, we need you to prompt the author/or reviewer(s) to revisit the submission if there has been no response within 7-10 days unless there's a clear statement in the review thread that says an action is coming at a slightly later time, perhaps because a reviewer committed to a review by a certain date, or an author is making changes and says they will be done by a certain date.**

[Whedon has functionality](https://joss.readthedocs.io/en/latest/whedon.html#reminding-reviewers-and-authors) to remind an author or review to return to a review at a certain point in the future. For example:

```
@whedon remind @reviewer in five days
```

## Out of office

Sometimes we need time away from our editing duties at JOSS. The [joss-reviews](https://github.com/openjournals/joss-reviews) repository has the [OoO bot](https://github.com/swinton/probot-ooo) installed which means you can mark yourself as out of the office (and unable to respond to reviews) for a period of time e.g.:

Mark yourself as OoO in one of the reviews you're editing in the [joss-reviews](https://github.com/openjournals/joss-reviews) repository like this:

```
/ooo January 18 until February 2
```

Ooo bot will then respond to any mentions in the [joss-reviews](https://github.com/openjournals/joss-reviews) repository to let people know you're away.

**Note, if you're planning on being out of the office for more than two weeks, please let the JOSS editorial team know.**

## Editorial buddy

New editors are assigned an editorial 'buddy' from the existing editorial team. The buddy is there to help the new editor onboard successfully and to provide a dedicated resource for any questions they might have but don't feel comfortable posting to the editor mailing list.

Buddy assignments don't have a fixed term but generally require a committment for 1-2 months.

Some things you might need to do as a buddy for a new editor:

- Respond to questions via email or on GitHub review issues.
- Check in with the new editor every couple of weeks if there hasn't been any other communication.
- (Optionally) keep an eye on the new editor's submissions.

## Managing notifications

Being on the JOSS editorial team means that there can be a _lot_ of notifications from GitHub if you don't take some proactive steps to minimize noise from the reviews repository.

### Things you should do when joining the editorial team

**Unsubscribe from the reviews repository on GitHub**

When you're added to the editorial team on GitHub, you will almost certainly find yourself subscribed (watching) to the [`joss-reviews`](https://github.com/openjournals/joss-reviews) repository. The first thing you should do is set yourself to 'not watching':

![Repository notifications settings](https://cloud.githubusercontent.com/assets/4483/20250593/64d7ce48-a9de-11e6-9225-d3dfb3e48e68.png)

Please note, that by not watching the reviews repository, you will still receive notifications for issues (reviews) where you are `@mentioned`.

Sometimes another editor might mention you in a review issue (for example to ask you a question). If you've responded and no-longer want to receive messages for that review, you can manually unsubscribe by clicking the button in the right-hand column on the review issue page.

**Curate your GitHub notifications experience**

GitHub has extensive documentation on [managing notifications](https://help.github.com/en/articles/managing-your-notifications) which explains when and why different notifications are sent from a repository.

**Set up email filters**

Email filters can be very useful for managing incoming email notifications, here are some recommended resources:

- A GitHub blog post describing how to set up [email filters](https://github.blog/2017-07-18-managing-large-numbers-of-github-notifications/).

If you use Gmail:

- https://gist.github.com/ldez/bd6e6401ad0855e6c0de6da19a8c50b5
- https://github.com/jessfraz/gmailfilters
- https://hackernoon.com/how-to-never-miss-a-github-mention-fdd5a0f9ab6d

**Use a dedicated tool**

For papers that you are already assigned to edit, the dedicated JOSS dashboard aggregates notifications associated with each paper. The dashboard is available at: `https://joss.theoj.org/dashboard/<yourgithubusername>`

Another tool you might want to try out is [Octobox](https://octobox.io/).
Reviewing for JOSS
=======================

Firstly, thank you so much for agreeing to review for the Journal of Open Source Software (JOSS), we're delighted to have your help. This document is designed to outline our editorial guidelines and help you understand our requirements for accepting a submission into the JOSS. Our review process is based on a tried-and-tested approach of the [rOpenSci collaboration](http://ropensci.org/blog/2016/03/28/software-review).

## Guiding principles

We like to think of JOSS as a 'developer friendly' journal. That is, if the submitting authors have followed best practices (have documentation, tests, continuous integration, and a license) then their review should be rapid.

For those submissions that don't quite meet the bar, please try to give clear feedback on how authors could improve their submission. A key goal of JOSS is to raise the quality of research software generally and you (the experienced reviewer) are well placed to give this feedback.

A JOSS review involves checking submissions against a checklist of essential software features and details in the submitted paper. This should be objective, not subjective; it should be based on the materials in the submission as perceived without distortion by personal feelings, prejudices, or interpretations.

We encourage reviewers to file issues against the submitted repository's issue tracker. **When you have completed your review, please leave a comment in the review issue saying so.**

You can include in your review links to any new issues that you the reviewer believe to be impeding the acceptance of the repository. (Similarly, if the submitted repository is a GitHub repository, mentioning the review issue URL in the submitted repository's issue tracker will create a mention in the review issue's history.)

## JOSS Conflict of Interest Policy

The definition of a conflict of Interest in peer review is a circumstance that makes you "unable to make an impartial scientific judgment or evaluation." ([PNAS Conflict of Interest Policy](http://www.pnas.org/site/authors/coi.xhtml)). JOSS is concerned with avoiding any actual conflicts of interest, and being sufficiently transparent that we avoid the appearance of conflicts of interest as well.

As a reviewer (or editor), COIs are your present or previous association with any authors of a submission: recent (past four years) collaborators in funded research or work that is published; and lifetime for the family members, business partners, and thesis student/advisor or mentor. In addition, your recent (past year) association with the same organization of a submitter is a COI, for example, being employed at the same institution.

If you have a conflict of interest with a submission, you should disclose the specific reason to the submissions' editor (or to an EiC if you are an editor, or to the other EiCs if you are an EiC). This may lead to you not being able to review or edit the submission, but some conflicts may be recorded and then waived, and if you think you are able to make an impartial assessment of the work, you should request that the conflict be waived. For example, if you and a submitter were two of 2000 authors of a high energy physics paper but did not actually collaborate. Or if you and a submitter worked together 6 years ago, but due to delays in the publishing industry, a paper from that collaboration with both of you as authors was published 2 year ago. Or if you and a submitter are both employed by the same very large organization but in different units without any knowledge of each other.

Declaring actual, perceived, and potential conflicts of interest is required under professional ethics. If in doubt: ask the editors.
Review criteria
===============

## The JOSS paper

As outlined in the [submission guidelines provided to authors](submitting.html#what-should-my-paper-contain), the JOSS paper (the compiled PDF associated with this submission) should only include:

- A list of the authors of the software and their affiliations.
- A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.
- A clear statement of need that illustrates the purpose of the software.
- A description of how this software compares to other commonly-used packages in this research area.
- Mentions (if applicable) of any ongoing research projects using the software or recent scholarly publications enabled by it.
- A list of key references including a link to the software archive.

```eval_rst
.. important:: Note the paper *should not* include software documentation such as API (Application Programming Interface) functionality, as this should be outlined in the software documentation.
```

## Review items

```eval_rst
.. important:: Note, if you've not yet been involved in a JOSS review, you can see an example JOSS review checklist `here <review_checklist.html>`_.
```

### Software license

There should be an [OSI approved](https://opensource.org/licenses/alphabetical) license included in the repository. Common licenses such as those listed on [choosealicense.com](https://choosealicense.com) are preferred. Note there should be an actual license file present in the repository not just a reference to the license.

> **Acceptable:** A plain-text LICENSE or COPYING file with the contents of an OSI approved license<br />            
> **Not acceptable:** A phrase such as 'MIT license' in a README file

### Substantial scholarly effort

Reviewers should verify that the software represents substantial scholarly effort. As a rule of thumb, JOSS' minimum allowable contribution should represent **not less than** three months of work for an individual. Signals of effort may include: 

- Age of software (is this a well-established software project) / length of commit history.
- Number of commits.
- Number of authors.
- Lines of code (LOC): These statistics are usually reported by Whedon in the `pre-review` issue thread.
- Whether the software has already been cited in academic papers.
- Whether the software is sufficiently useful that it is _likely to be cited_ by other researchers working in this domain.

These guidelines are not meant to be strictly prescriptive. Recently released software may not have been around long enough to gather citations in academic literature. While some authors contribute openly and accrue a long and rich commit history before submitting, others may upload their software to GitHub shortly before submitting their JOSS paper.  Reviewers should rely on their expert understanding of their domain to judge whether the software is of broad interest (_likely to be cited by other researchers_) or more narrowly focused around the needs of an individual researcher or lab.

```eval_rst
.. note:: The decision on scholarly effort is ultimately one made by JOSS editors. Reviewers are asked to flag submissions of questionable scope during the review process so that the editor can bring this to the attention of the JOSS editorial team.
```

### Documentation

There should be sufficient documentation for you, the reviewer to understand the core functionality of the software under review. A high-level overview of this documentation should be included in a `README` file (or equivalent). There should be:

#### A statement of need

The authors should clearly state what problems the software is designed to solve and who the target audience is.

#### Installation instructions

Software dependencies should be clearly documented and their installation handled by an automated proceedure. Where possible, software installation should be managed with a package manager. However, this criterion depends upon the maturity and availability of software packaging and distribution in the programming language being used. For example, Python packages should be `pip install`able, and have adopted [packaging conventions](https://packaging.python.org), while Fortran submissions with a Makefile may be sufficient.

> **Good:** The software is simple to install, and follows established distribution and dependency management approaches for the language being used<br />
> **OK:** A list of dependencies to install, together with some kind of script to handle their installation (e.g., a Makefile<br />
> **Bad (not acceptable):** Dependencies are unclear, and/or installation process lacks automation

#### Example usage

The authors should include examples of how to use the software (ideally to solve real-world analysis problems).

#### API documentation

Reviewers should check that the software API is documented to a suitable level.

> **Good:** All functions/methods are documented including example inputs and outputs<br />
> **OK:** Core API functionality is documented<br />
> **Bad (not acceptable):** API is undocumented

```eval_rst
.. note:: The decision on API documentation is left largely to the discretion of the reviewer and their experience of evaluating the software.
```

#### Community guidelines

There should be clear guidelines for third-parties wishing to:

- Contribute to the software
- Report issues or problems with the software
- Seek support

### Functionality

Reviewers are expected to install the software they are reviewing and to verify the core functionality of the software.

### Tests

Authors are strongly encouraged to include an automated test suite covering the core functionality of their software.

> **Good:** An automated test suite hooked up to continuous integration (GitHub Actions, Circle CI, or similar)<br />
> **OK:** Documented manual steps that can be followed to objectively check the expected functionality of the software (e.g., a sample input file to assert behavior)<br />
> **Bad (not acceptable):** No way for you, the reviewer, to objectively assess whether the software works

## Other considerations

### Authorship

As part of the review process, you are asked to check whether the submitting author has made a 'substantial contribution' to the submitted software (as determined by the commit history) and to check that 'the full list of paper authors seems appropriate and complete?'

As discussed in the [submission guidelines for authors](submitting.html#authorship), authorship is a complex topic with different practices in different communities.  Ultimately, the authors themselves are responsible for deciding which contributions are sufficient for co-authorship, although JOSS policy is that purely financial contributions are not considered sufficient. Your job as a reviewer is to check that the list of authors appears reasonable, and if it's not obviously complete/correct, to raise this as a question during the review.

### An important note about 'novel' software and citations of relevant work

Submissions that implement solutions already solved in other software packages are accepted into JOSS provided that they meet the criteria listed above and cite prior similar work. Reviewers should point out relevant published work which is not yet cited.

### What happens if the software I'm reviewing doesn't meet the JOSS criteria?

We ask that reviewers grade submissions in one of three categories: 1) Accept 2) Minor Revisions 3) Major Revisions. Unlike some journals we do not reject outright submissions requiring major revisions - we're more than happy to give the author as long as they need to make these modifications/improvements.

### What about submissions that rely upon proprietary languages/development environments?

As outlined in our author guidelines, submissions that rely upon a proprietary/closed source language or development environment are acceptable provided that they meet the other submission requirements and that you, the reviewer, are able to install the software & verify the functionality of the submission as required by our reviewer guidelines.

If an open source or free variant of the programming language exists, feel free to encourage the submitting author to consider making their software compatible with the open source/free variant.
Journal of Open Source Software
===============================

The `Journal of Open Source Software
<http://joss.theoj.org/>`_ (JOSS) is a developer friendly journal for research software packages.

JOSS is an academic journal (ISSN 2475-9066) with a formal peer-review process that is designed to *improve the quality of the software submitted*. Upon acceptance into JOSS, we mint a CrossRef DOI for your paper and we list it on the `JOSS website <http://joss.theoj.org/>`_.

About this site
---------------

This site contains documentation for authors interested in submitting to JOSS, reviewers who have generously volunteered their time to review submissions, and editors who manage the JOSS editorial process.

If you're interested in learning more about JOSS, you might want to read:

- `Our announcement blog post <http://www.arfon.org/announcing-the-journal-of-open-source-software>`_ describing some of the motivations for starting a new journal
- `The paper in Computing in Science and Engineering <https://doi.org/10.1109/MCSE.2018.03221930>`_ introducing JOSS
- `The paper in PeerJ CS <https://doi.org/10.7717/peerj-cs.147>`_ describing the first year of JOSS
- The `about page <http://joss.theoj.org/about>`_ on the main JOSS site

Submitting a paper to JOSS
-----------------------

If you'd like to submit a paper to JOSS, please read the author submission guidelines in the :doc:`submitting` section.

Sponsors and affiliates
-----------------------

.. image:: images/sponsors.png
  :alt: OSI & NumFOCUS

JOSS is a proud affiliate of the `Open Source Initiative <https://opensource.org/>`_. As such, we are committed to public support for open source software and the role OSI plays therein. In addition, Open Journals (the parent entity behind JOSS) is a `NumFOCUS-sponsored project <https://www.numfocus.org/project/open-journals>`_.

.. toctree::
  :caption: Author and Reviewer Guides
  :maxdepth: 5

  submitting
  reviewer_guidelines
  review_criteria
  review_checklist

.. toctree::
  :caption: Editor Guides
  :maxdepth: 2

  editing

.. toctree::
  :caption: The Whedon Editorial Bot
  :maxdepth: 2

  whedon

.. toctree::
  :caption: Developer Guides
  :maxdepth: 2

  installing
