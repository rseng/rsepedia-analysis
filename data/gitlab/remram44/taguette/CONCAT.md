Changelog
=========

1.2.0 (2022-01-22)
------------------

Bugfixes:
* Add some database indexes, improving performance on larger installations
* Reload page when the server updates, avoiding clients on older versions if they never close the tab
* Fix error importing SQLite3 file on Windows (and other permission errors on temporary files)
* Show a helpful message if the imported codebook has rows that are too short, instead of 500 error
* Fix a 403 error when importing a project if no other form has been shown before (single-user mode, no project created)

Enhancements:
* Focus the name (path) field when creating a tag
* Add a keyboard shortcut to create a highlight: 'n' key (localizable)
* Enable portable install in Windows installer (no admin access)

1.1.1 (2021-10-17)
------------------

Bugfixes:
* Don't cache documents forever in desktop mode (only server mode)

1.1.0 (2021-10-12)
------------------

Bugfixes:
* Avoid error "document name too long" if no name was provided and the filename is used
* Don't deselect all tags from the highlight modal if it is open while creating a new tag
* Avoid `TimeoutError` due to the SQLAlchemy connection pool overflowing (mostly affects busy self-hosted non-SQLite deployments)

Enhancements:
* Add support for right-to-left documents
* Add importing tags from a codebook CSV file

1.0.1 (2021-07-22)
------------------

This is a pip-only fix: binary installers, Docker, Poetry install are NOT affected.

Bugfixes:
* Fix incompatibility with SQLAlchemy 1.4

1.0.0 (2021-07-17)
------------------

Bugfixes:
* Fix highlights disappearing when one of their tags is deleted
* Fix highlight modal showing "no tags" multiple times after removing the last tag repeatedly
* Add missing title in exported highlights HTML

Enhancements:
* Allow disabling SQLite3 file import from the configuration file
* Add option to use Redis to coordinate live collaboration if running Taguette in multiple processes

0.11 (2021-07-06)
-----------------

Bugfixes:
* Fix timezone issues on some SQL backends (PostgreSQL, not SQLite). Could cause recorded times to be wrong, and in some cases password reset links could be used multiple times
* Fix the activity detection used to pause/resume polling
* Fix default-config command outputing logs about umask, which could be piped into the new config file (for example when using Docker with `-t`)
* Fix migrations on SQLite failing because of foreign key constraints when recreating tables
* Fix RTF export
* Fix retrying Calibre conversion with --enable-heuristics again
* Fix Calibre leaving behind temporary files if it times out

Enhancements:
* Show more details when document export or import fails
* Show the number of highlights per tag in exported codebooks
* Performance improvement of getting a document or tag with many highlights
* Add terms of service, which you can optionally set in your own instance
* Made document/highlights/codebook exporting functions available under `taguette.export`, for use in scripts and notebooks
* Improve button placement in the left panel (put "add a document" and "create a tag" at the top)
* Add an HTML page for 404 errors
* Improve error message when Calibre is not found
* Add support for MariaDB (MySQL doesn't work because of error 1093, see https://stackoverflow.com/q/4429319)
* Add option to import or export a project as a SQLite3 database
* Allow a user to remove themselves from a project, without having to ask an admin

0.10.1 (2021-02-22)
-------------------

Bugfixes:
* Limit the number of page buttons
* SQL query performance improvements
* Force long document names to wrap (like long tag names)

0.10 (2021-02-17)
-----------------

Bugfixes:
* Check that tags exist in the project before adding them to highlights

Enhancements:
* SEO: Add page titles, robots.txt
* Remind user of their login in the password reset email
* Pause polling for changes when the window has been inactive for 10 minutes
* Improve the error message shown when document import fails
* Add pagination to the highlights view, which greatly improve performances for heavily used tags
* Update to Tornado 6.1, remove workaround for older versions

0.9.2 (2020-08-26)
------------------

Bugfixes:
* Fix document conversion on Windows

Enhancements:
* Show spinner while loading document or tag view
* Show the tag names next to each highlight in document export

0.9.1 (2020-08-24)
------------------

Bugfixes:
* Fix printed link not being flushed, causing it to be inaccessible e.g. running on Docker
* Fix incompatibility with Python 3.8
* Close DB connection during long polls, avoiding overflow of connection pool when using PostgreSQL
* Don't allow a password reset link to be used more than once
* Fix 'new highlight' button on Microsoft Edge
* If Calibre fails, run it again without heuristics
* Fix highlights being off in document export in non-ASCII paragraphs

Enhancements:
* Set umask to 077 by default (and add corresponding command-line options)
* Restore Python 3.5 compatibility
* Add a timeout to document conversion
* Add Calibre output limits to the config
* Fix some important messages being routed through log instead of direct to terminal

0.9 (2019-11-23)
----------------

Bugfixes:
* Fix showing highlights for a tag with a slash in it
* Fix new tags showing NaN highlight count
* Fix export button being on top of text, if the first line is long
* Fix tag names containing '#' not working properly
* Fix "new highlight" button showing at the top of the page when scrolling on Safari
* Sort highlights in the exported highlights view as well
* Fix merging tags that have common highlights (previous 500 error)

Features:
* Export codebook to Excel (.xlsx)
* Show error details in alert messages (in English)
* Highlight the currently selected tag(s) in the left panel
* Don't show 403 if user can't change collaborators: show explanation, hide button
* Add export of highlights to CSV and Excel (.xlsx) formats
* Support PBKDF2 passwords (Taguette no longer requires 'bcrypt')

0.8 (2019-06-15)
----------------

Bugfixes:
* Don't show 500 error on invalid email reset token
* Explicitly close DB connections, which might help with some warnings
* Don't show 'merge' button in modal when creating a new tag
* Fix getting logged out in single-user mode with `--debug`
* Don't scroll to the top of the document when clicking on a disabled link
* Fix taguette --database=filename not working when filename does not contain directories

Features:
* Add limits on converted file size
* Don't have Calibre export image files from PDF, since we don't read them
* Add a scrollbar to modals, since they can grow big in projects with many tags
* Use the file name as document name if left blank
* Show cookie warning before setting any (optional in configuration)
* Add the REFI-QDA Codebook (.qdc) export format
* Improve the collaborator management modal
* Show the number of highlights which each tag in the "highlights" panel

0.7 (2019-05-15)
----------------

Taguette can now be translated! You can help bring Taguette to your language on [Transifex](http://transifex.com/remram44/taguette/).

Bugfixes:
* Fix exporting highlights for non-ASCII tags
* Fix account page not accepting empty optional fields
* Fix document description being validated as its name
* Fix importing documents with completely non-ASCII filenames

Features:
* Merge tags
* Added internationalization
* French translation
* German translation
* Spanish translation
* Show tag names when hovering a highlight

0.6 (2019-04-13)
----------------

Bugfixes:
* Make 'display' headings responsive
* Fix exported highlights being called "path"
* Fix possible weird characters in exported documents on Windows (depending on locale)

Features:
* Convert logins to lower-case (login and collaborator forms will convert too, so it should only affect display). Users with non-lowercase logins will be logged out on update
* Moved the `SECRET_KEY` to the config, no longer writing to `~/.cache`
* Let you know when you have been logged out or removed from a project while working

0.5 (2019-03-23)
----------------

Bugfixes:
* Improve reading of OPF output from Calibre, which might fix compatibility with some combinations of Calibre versions and input formats
* Long tag names no longer stick out of the left pane
* Sort tags in highlight modal, documents in left pane, highlights and their tags in the highlights view

Features:
* Use a configuration file in server mode, rather than a growing list of command-line options
* Expose metrics to Prometheus
* Send errors to Sentry
* Add 'delete project' button
* Add account management page, to update email/password
* Add password recovery feature (if you have an email set)
* "New highlight" button shows up next to selected text rather than mouse, making it work with touch screens (mobile) and screen readers (hopefully)
* Convert old .DOC files (Word 97) using WV if available
* Add collaborator management modal, to add more members to a project
* Changed default port number from `8000` to `7465`
* Add spinning icon while requests are in progress, to prevent multiple submission of forms (document add takes ~10s for example)

0.4.4 (2018-11-29)
------------------

Bugfixes:
* Fix error creating a highlight when a paragraph is selected to the end
* Correctly handle Calibre sometimes writing a `.xhtml` file instead of `.html`
* Also show highlights with no tags when selecting "See all highlights"

0.4.3 (2018-11-17)
------------------

Bugfixes:
* Fix JavaScript error on a brand new project (no recorded Command)
* Fix Commands being sent to the wrong project

Features:
* Add `--xheaders` option for the hosted setup, showing correct IPs in the log
* Show which document in the list is the current one

0.4.2 (2018-11-15)
------------------

Bugfixes:
* Don't show highlights from a different document
* Handle non-ascii text better
* Fix real-time updates pausing if two changes happen in the same second
* Add messages boxes to signal when something goes wrong
* Sanitize name of uploaded files

0.4.1 (2018-11-12)
------------------

Bugfixes:
* Log errors from async handlers to the console instead of hiding them
* Work around a problem computing highlight positions in documents when unicode is present (won't crash anymore, but positions might still be off, fix to come)
* Fix not being able to create tags with names that collide with tags in other projects, and error creating a project if another project still uses default tags
* Fix exporting a document that has 0 highlights
* Fix document export missing highlights
* Fix navbar expand button (shown on smaller screen sizes) not working

0.4 (2018-11-11)
----------------

Bugfixes:
* Make sure to hide auth token from URL bar
* Fix tag description not showing up
* Don't allow two tags to have the same name

Features:
* Add confirmation dialogs before deleting tags or documents
* Create tag from the highlight window
* New theme matching website (thanks to Vicky Steeves)
* Add button to create a tag from the highlight modal
* Show messages from taguette.fr, such as new version available
* Add export options, allowing you to get your highlights or highlighted documents as HTML, DOCX, or PDF
* Add HTML and PDF options for codebook export as well
* Add an option to show (or export) all highlights, rather than only a specific tag

0.3 (2018-10-29)
----------------

Bugfixes:
* Fix having to reload for changes to appear when working on project other than 1.
* Fix tags not being sorted by name

Features:
* Add 'backlight' mode, fading non-tagged text
* Add modal dialog to edit and delete documents
* Add migration system, to automatically upgrade the database to new schema version when required

0.2 (2018-10-21)
----------------

Bugfixes:
* Accept list and numbered lists, as generated from Markdown documents
* Fix tag modal not able to add tags after a tag has been edited

Features:
* Add single-user mode, the default. Multi-user mode now needs `--multiuser`
* Add login and registration pages
* Add codebook export to CSV and DOCX files (contains list of tags with their descriptions)

0.1 (2018-10-21)
----------------

First version, proof of concept. Not very useful, but showcases the app, and can be installed by alpha testers.

* Can create projects
* Can import documents into the database as HTML
* Uses Calibre to convert supported documents into HTML
* Can highlight parts of documents, and assign tags
* Real-time notifications and collaboration
* "Acceptable" UI with bootstrap
# Introduction

Taguette is a web-based Python application. It is meant to run both stand-alone on an end-user machine, where it must run automatically without having to configure databases or accounts; and as a service over the internet.

To learn more about Taguette, see the [Taguette Website](https://www.taguette.org/) and the [README](https://gitlab.com/remram44/taguette/blob/master/README.rst) in this repository. You can find general contributing instructions in [CONTRIBUTING.md](https://gitlab.com/remram44/taguette/blob/master/CONTRIBUTING.md).

# Architecture

![Architecture diagram](architecture.svg)

## The server

A server, written in Python, handles most of the business logic. It is responsible for importing documents, adding tags and highlights, etc. This is the process you are starting when you run `taguette`.

We use Python 3.5+, the [Tornado web framework](https://www.tornadoweb.org/), [SQLAlchemy](https://www.sqlalchemy.org/), and [Jinja](https://jinja.palletsprojects.com/) to render HTML pages from templates (Tornado's built-in template language is quite limited).

On startup, the server asks the system to launch a web browser pointing to itself. The browser downloads "the frontend", HTML pages generated by the server (from templates), and receives requests sent from JavaScript code to the API to perform actions in a project. You can find the templates under `taguette/templates/` and the CSS and JavaScript files under `taguette/static/`.

## The database

The data is stored in a SQL database. Multiple database systems are supported (see [SQLAlchemy's documentation](https://docs.sqlalchemy.org/en/latest/core/engines.html)). SQLite3 is the default when running in single-user mode, PostgreSQL is used for [app.taguette.org](https://app.taguette.org/). MariaDB is tested and known to work. MySQL does NOT work (particularly the query used to merge tags).

![Database diagram](erd.svg)

[Alembic](https://alembic.sqlalchemy.org/) is used to migrate the database, which is to say apply changes to the data or the schema. If you upgrade your version of Taguette, it will read the version of the database from the table `alembic_version`, and apply migrations to take it to the latest version. This is done automatically in single-user mode (there is a warning in the terminal, and a backup of the SQLite3 file is made first), but in server mode, it will refuse to start and instruct you to use the `taguette migrate` command. You can find the migrations in `taguette/migrations/versions/`.

Highlights are represented by their position (`start_offset`, `end_offset`) in the document's text (skipping over HTML tags), counted in bytes in the UTF-8 text (see `taguette/extract.py`). The content of each highlight is stored in the `highlights` table, so we don't have to extract that part of the text from the `documents` table on each access.

## The frontend

The frontend consists of a few different views, rendered from Jinja templates. The most important one is the project view from which a user does all their work. The JavaScript is kept in a single file, `static/js/taguette.js`.

Every user action results in a call to the server over a JSON API. The server then syncs up every browser (including the acting user's) by sending an event over the long-polling endpoint; see [live editing](#live-editing) below.

# Live editing

To provide real-time collaboration between users, a web browser needs to be alerted of other users' changes immediately and update the view. To this end, each web browser maintains a long-lived connection to the server. Instead of replying to the web request immediately, the server keeps that request open until something happens and replies then. The frontend updates using that response and starts a new request. This is called [long-polling](https://en.wikipedia.org/wiki/Long_polling).

The events sent over this endpoint are the `commands` stored in the database. When a browser asks for new changes, it provides the ID of the last command it knows of. The server first checks the database, as there may already be changes to be sent (for example if the browser has been disconnected, the computer was asleep, or the browser tab was inactive) and otherwise waits in the long-polling fashion.

In practice, we also send those events to the browser that initiated those changes, so we don't have to have two separate mechanisms to update the view.

# Localization

Translation files are found in `po/`. The `.pot` files are generated from the code using [Babel](http://babel.pocoo.org/) (see `scripts/update_pot.sh`) and uploaded to [Transifex](https://www.transifex.com/remram44/taguette/) where users can help translate strings to their language. Those files are compiled to `.mo` files which are loaded by Tornado (see `scripts/update_translations.sh`).

By default, a user has no language set and Taguette will use the language preferences selected in the web browser's settings. Additionally, a user can select a language in their account settings, which will be used instead.

The translations are split between two catalogs `main` and `javascript`. The first one is for strings found throughout the server, for example in templates. The second one is for use in the frontend's JavaScript code; JavaScript downloads this catalog from the `/trans.js` endpoint, which gives a dictionary mapping English strings to their counterpart in the user's language.

# Conversion

When a document is added to a project, it first needs to be converted into a format suitable for viewing and highlighting.

The first step is to turn it into HTML. This is done by using [Calibre](https://calibre-ebook.com/)'s `ebook-convert` command. Calibre is run as a subprocess and given a 2-minute timeout (by default).

For old Microsoft Word 97 `.doc` files, `wvHtml` from the [wvWare package](http://wvware.sourceforge.net/)) is used instead.

The second step (unless the document was already HTML) is to sanitize this HTML to remove content that would not be safe or convenient to serve to users, such as scripts, media, etc. This is done using the [bleach](https://github.com/mozilla/bleach) library.

Finally, the document is inserted into the database.

In `taguette/convert.py` you will also find conversion functions going the other way, from HTML to other formats. Those are used when exporting data from Taguette.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at hi@taguette.org. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Contents
* [Notes](#notes)
* [Contributing](#contributing)
* [Resolving Merge Conflicts](#resolving-merge-conflicts)
* [Best Practices for Contributing](#best-practices-for-contributing)
* [Attribution](#attribution)

# Notes

Any contributions received are assumed to be covered by the [BSD 3-Clause license](https://gitlab.com/remram44/taguette/blob/master/LICENSE.txt). We might ask you to sign a Contributor License Agreement before accepting a large contribution. To learn more about Taguette, see the [Taguette Website](https://www.taguette.org/) and the [README](https://gitlab.com/remram44/taguette/blob/master/README.rst) in this repository. You can find an introduction to the codebase in [ARCHITECTURE.md](https://gitlab.com/remram44/taguette/blob/master/ARCHITECTURE.md).

# Contributing

Please follow the [Contributor Covenant](CODE_OF_CONDUCT.md) in all your interactions with the project. If you would like to contribute to this project by modifying/adding to the code, please read the [Best Practices for Contributing](#best-practices-for-contributing) below and feel free to use [GitLab's Web IDE](https://docs.gitlab.com/ee/user/project/web_ide/) for small changes (e.g. fixing a typo on the website) or for contributions to the code base the following workflow:

1. Fork the project.
2. Clone your fork to your computer.
    * From the command line: `git clone https://gitlab.com/<USERNAME>/taguette.git`
3. Change into your new project folder.
    * From the command line: `cd taguette`
4. [optional]  Add the upstream repository to your list of remotes.
    * From the command line: `git remote add upstream https://gitlab.com/remram44/taguette.git`
5. Create a branch for your new feature.
    * From the command line: `git checkout -b my-feature-branch-name`
6. Make your changes.
    * Avoid making changes to more files than necessary for your feature (i.e. refrain from combining your "real" merge request with incidental bug fixes). This will simplify the merging process and make your changes clearer.
7. Commit your changes. From the command line:
    * `git add <FILE-NAMES>`
    * `git commit -m "A descriptive commit message"`
8. While you were working some other changes might have gone in and break your stuff or vice versa. This can be a *merge conflict* but also conflicting behavior or code. Before you test, merge with master.
    * `git fetch upstream`
    * `git merge upstream/master`
9. Test. Run the program and do something related to your feature/fix.
10. Push the branch, uploading it to GitLab.
    * `git push origin my-feature-branch-name`
11. Make a "merge request" from your branch here on GitLab.

# Resolving Merge Conflicts

Depending on the order that merge requests get processed, your MR may result in a conflict and become un-mergable.  To correct this, do the following from the command line:

Switch to your branch: `git checkout my-feature-branch-name`
Pull in the latest upstream changes: `git pull upstream master`
Find out what files have a conflict: `git status`

Edit the conflicting file(s) and look for a block that looks like this:
```
<<<<<<< HEAD
my awesome change
=======
some other person's less awesome change
>>>>>>> some-branch
```

Replace all five (or more) lines with the correct version (yours, theirs, or
a combination of the two).  ONLY the correct content should remain (none of
that `<<<<< HEAD` stuff.)

Then re-commit and re-push the file.

```
git add the-changed-file.cs
git commit -m "Resolved conflict between this and MR #123"
git push origin my-feature-branch-name
```

The merge request should automatically update to reflect your changes.

## Best Practices for Contributing

* Before you start coding, open an issue so that the community can discuss your change to ensure it is in line with the goals of the project and not being worked on by someone else. This allows for discussion and fine-tuning of your feature and results in a more succint and focused addition.
    * If you are fixing a small glitch or bug, you may make a MR without opening an issue.
    * If you are adding a large feature, create an issue so that we may give you feedback and agree on what makes the most sense for the project before making your change and submitting a MR (this will make sure you don't have to do major changes down the line).

* Merge requests are eventually merged into the codebase. Please ensure they are:
    * Well tested by the author. It is the author's job to ensure their code works as expected.
    * Free of unnecessary log calls. Logging is important for debugging, but when a MR is made, log calls should only be present when there is an actual error or to record some important piece of information or progress.

* If your code is untested, log-heavy, or incomplete, prefix your MR with "[WIP]", so others know it is still being tested and shouldn't be considered for merging yet. This way we can still give you feedback or help finalize the feature even if it's not ready for prime time.

That's it! Following these guidelines will ensure that your additions are approved quickly and integrated into the project. Thanks for your contribution!

## Running tests

This projects includes automated tests. Should you send us your changes in the form of a merge request on GitLab, the tests will be run automatically by GitLab CI to check that your version still works. You can also run the tests locally in a terminal:

```
python tests.py
```

Some tests control a web browser. You will need to install Chrome/Chromium or Firefox and get the corresponding webdriver ([geckodriver](https://github.com/mozilla/geckodriver) for Firefox and [chromedriver](https://chromedriver.chromium.org/downloads) for Chrome/Chromium). You can then enable those tests by running:

```
TAGUETTE_TEST_WEBDRIVER=firefox python tests.py
# or
TAGUETTE_TEST_WEBDRIVER=chromium python tests.py
```

## Debugging Taguette with PyCharm

1. Install the poetry plugin for pycharm: https://plugins.jetbrains.com/plugin/14307-poetry
2. In PyCharm's menu: Run/Edit Configurations
3. Add a new configuration by clicking the top-left plus icon or by pressing alt+insert on your keyboard.
4. Select `Python` in the pop-up list (other options that show are `Shell Script`, `Python tests`, etc.)
5. In _Module Name_ type `taguette.main`
6. In _Environment_, _Python interpreter_, choose "Poetry (taguette)"

You can now just click on _Run/Debug_ (if you only have one run configuration) or click on _Run/Debug..._ and choose _taguette.main_ to start the debugging.

# Attribution

This CONTRIBUTING.md was adapted from [ProjectPorcupine](https://github.com/TeamPorcupine/ProjectPorcupine)'s [CONTRIBUTING.md](https://github.com/TeamPorcupine/ProjectPorcupine/blob/master/CONTRIBUTING.md)

# Contact info

You are welcome to chat with us in our [Element room](https://app.element.io/#/room/#taguette:matrix.org) or contact the maintainers directly at [hi@taguette.org](mailto:hi@taguette.org).
Translating Taguette
====================

Taguette is being used globally. Therefore, it needs to be translatable into multiple languages.

What should be translated:

* All the text that is visible on web pages during normal usage (templates, messages, JavaScript messages)
* The messages printed to the terminal and the command-line help

What should NOT be translated:

* Logging and exception messages (those are really destined to developers)
* API messages and errors

There are two catalogs: `main` is used by the server, and `javascript` is loaded client-side as JSON to be used by the JavaScript code.

Run the script `scripts/update_translations.sh` before and after updating PO files. It will update the POT catalogs with the changes from the code and also generate the MO files from the translated PO files.

You will have to update your PO files from the POT files manually.

Standard mapping
================

This table intends to help standardize the terms used within a language translation, and between different languages.

| English   | French      |
| :-------- | :---------- |
| codebook  | *codebook*  |
| tag       | tag         |
| highlight | *marque*    |
| account   | profil      |
| settings  | préférences |
| email     | e-mail      |
Taguette
========

A spin on the phrase "tag it!", `Taguette <https://www.taguette.org/>`__ is a free and open source qualitative research tool that allows users to:

+ Import PDFs, Word Docs (``.docx``), Text files (``.txt``), HTML, EPUB, MOBI, Open Documents (``.odt``), and Rich Text Files (``.rtf``).
+ Highlight words, sentences, or paragraphs and tag them with the codes *you* create.
+ (not yet) Group imported documents together (e.g. as 'Interview' or 'Lit Review').
+ Export tagged documents, highlights for a specific tag, a list of tags with description and colors, and whole projects.

`Check out our website to learn more about how to install and get started. <https://www.taguette.org/>`__

Motivation and goal
-------------------

Qualitative methods generate rich, detailed research materials that leave individuals’ perspectives intact  as well as provide multiple contexts for understanding phenomenon under study. Qualitative methods are used by a wide range of fields, such as anthropology, education, nursing, psychology, sociology, and marketing. Qualitative data has a similarly wide range: observations, interviews, documents, audiovisual materials, and more.

However - the software options for qualitative researchers are either **far too expensive**, don't allow for the seminal method of highlighting and tagging materials, *or actually perform quantitative analysis*, just on text.

**It's not right or fair that qualitative researchers without massive research funds cannot afford the basic software to do their research.**

So, to bolster a fair and equitable entry into qualitative methods, we've made Taguette!

Installation
------------

You can find complete installation instructions on `our website <https://www.taguette.org/install.html>`__, including installers for Windows and MacOS.

Development setup from the repository
-------------------------------------

You can install from a local clone of this repository, which will allow you to easily change the sources to suit your needs:

1. Clone this git repository from the terminal: ``git clone https://gitlab.com/remram44/taguette.git``
2. Navigate on the command line to the repository you've just cloned locally, using the ``cd`` command. To get help using ``cd``, use `this tutorial <https://swcarpentry.github.io/shell-novice/02-filedir/index.html>`__.
3. Taguette uses `Poetry <https://python-poetry.org/>`__ for its packaging and dependency management. You will need to `install Poetry <https://python-poetry.org/docs/#installation>`__.
4. Install Taguette and its dependencies by running ``poetry install``. Poetry will create a virtual environment for you by default, activate it by running ``poetry shell``.
5. Build translation files using ``scripts/update_translations.sh``.
6. You can start taguette in development mode using ``taguette --debug`` (or ``taguette --debug server <config_file>``). This will start Tornado in debug mode, which means in particular that it will auto-restart every time you make changes.
7. Navigate to `localhost:7465 <http://localhost:7465/>`__ to use Taguette!

License
-------

* Copyright (C) 2018, Rémi Rampin and Taguette contributors

Licensed under a **BSD 3-clause "New" or "Revised" License**. See the `LICENSE <LICENSE.txt>`__ for details.
