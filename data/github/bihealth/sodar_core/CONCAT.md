---
title: 'SODAR Core: a Django-based framework for scientific data management and analysis web apps'
tags:
 - scientific data management
 - django
 - framework
 - python
authors:
 - name: Mikko Nieminen
   orcid: 0000-0002-4180-8810
   affiliation: 1
 - name: Oliver Stolpe
   orcid: 0000-0002-8686-2941
   affiliation: 1
 - name: Franziska Schumann
   orcid: 0000-0002-1332-4894
   affiliation: 1
 - name: Manuel Holtgrewe
   orcid: 0000-0002-3051-1763
   affiliation: 1
 - name: Dieter Beule
   orcid: 0000-0002-3284-0632
   affiliation: 1
affiliations:
 - name: "Berlin Institute of Health, Anna-Louisa-Karsch-Straße 2, 10178 Berlin"
   index: 1
date: 26 June 2019
bibliography: paper.bib
---

# Introduction

Modern life science is generating large data sets at a unprecedented speed.
A major source of data are so-called omics (e.g., genomics, metabolomics, or proteomics) experiments.
Consequently, management and analysis of scientific data has become a major challenge.
Further, the heterogeneity of projects makes "one size fits all" data analysis systems infeasible and calls for specialized data analysis platforms.

The authors are actively developing applications for the FAIR (findable, accessible, interoperable, and reuseable, cf. @Wilkinson:2016 data management of omics data and their analysis.
In order to prevent duplication of work, we have extracted the commonly useful components into *SODAR Core*, a Python framework on top of Django.
It is used in our actively developed applications and also proved useful in internal web app prototypes.

Examples for using SODAR Core in the development of scientific data management web apps is Digestiflow [@Holtgrewe:2019] and the Filesfolders module shipping with SODAR Core.
An example for using SODAR Core in the development of scientific data analysis web apps is VarFish [@Holtgrewe:2020].

# Features

SODAR Core provides its features in the form of Django apps.
This section lists the main applications at the time of writing and summarises their features.

## Projects and project-based access control

All SODAR Core features are grouped around projects that can be arranged in a tree structure.
Users can be assigned to projects and have different permission sets in each projects.
This feature is provided in the `projectroles` app.
Further, SODAR Core allows for easy integration of more than one external LDAP authentication source.

SODAR Core currently comes with a fixed set of roles that we found a appropriate for a number of projects.
Making the role set configurable is a part of our future development roadmap.

## Common Projects in Separated Applications

As mentioned above, we found that the complexity of omics datasets requires specialized data analysis applications.
There are several advantages of separating and decoupling these specialized apps from the main data management system.
These include the ability for having different rates of development and stability, e.g., having a separate core data management system and experimental and less well-tested data analysis platforms.
We also found the ability to reduce dependency of development teams from each other very useful.

On the other hand, project and authorization generally should be shared between the data management and data analysis projects.
To prevent users from having to manually replicate projects and authorization configuration SODAR Core provides a *remote site* feature.
An application instance can either be a *source* or a *target* site.
Target sites can be connected to a source site and obtain project and authorization information from the source site.

## Auditing and Timelines

The `timeline` app provides an API for storing events in a project-specific timeline.
This allows for implementing audit trails and a history or changelog.
We have found this very useful to track changes done by users and keep a provenance trail.

## Asynchronous Background Jobs

In many data analysis use cases, there will be a number of data processing tasks that take a longer time, e.g., a complex database query, copying large volumes of data, or longer numeric computation.
In the context of web applications, it is important to implement them as asynchronous jobs as implemented in the `bgjobs` app.
This allows the web app to render HTML documents to the client's browser independently of completing the data analysis task.

## Caching Infrastructure

In many cases, data management and analysis applications will use data stored in external services that might take a longer time to query.
For this use case, the `sodarcache` app provides features for caching query results.

## Miscellaneous Functionality

Further, there are a number of apps in SODAR core that help with various recurrent tasks.
The `userprofile` provide user profiles and together with `projectroles` it allows storing user and/or project-specific settings.
`siteinfo` allows to create a site-wide statistics view while `adminalerts` allows to display and schedule site-wide notifications, e.g., for announcing maintenances.

In addition, we have found the management of files (in the range of several MB) useful for many applications such as sharing documents or storing attachments to messages.
The app `filesfolders` provides the feature that allows users to upload files and manage them in a folder tree structure.

# Summary and Conclusion

SODAR Core addresses many recurring requirements for scientific data management and analysis projects.
The authors are actively using it in their projects.
We believe that other people will find the components useful in their projects and are thus publishing them as free open source software.

# License and Availability

SODAR Core is distributed under the MIT license and available from Github at https://github.com/bihealth/sodar_core.
Each release is also stored to Zenodo.
The current version 0.8.4 is available with the DOI 10.5281/zenodo.4269347.
An example Django site is contained in the SODAR Core Git repository and up-to date documentation is available at https://sodar-core.readthedocs.org.

# References
.. highlight:: shell

Contributing
^^^^^^^^^^^^

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
======================

Report Bugs
-----------

Report bugs at https://github.com/bihealth/sodar-core/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
--------

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
------------------

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
-------------------

SODAR Core could always use more documentation, whether as part of the
official SODAR Core docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
---------------

The best way to send feedback is to file an issue at https://github.com/bihealth/sodar-core/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
============

Ready to contribute? Here's how to set up SODAR Core for local development.

1. Fork the ``sodar-core`` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/sodar-core.git

3. Install your local copy into a virtualenv. Assuming you have
   virtualenvwrapper installed, this is how you set up your fork for local
   development::

    $ mkvirtualenv sodar-core
    $ cd sodar-core/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature dev

   Make sure you base your changes on the ``dev`` branch, which is the current
   active development branch. The ``main`` branch is intended for merging
   stable releases only. Now you can make your changes locally.

5. When you're done making changes, make sure to apply proper formatting using
   Black and the settings specified in the accompanying ``black.sh`` script.
   Next, check that your changes pass flake8 and the tests. It is recommended to
   use the ``Makefile`` to ensure the correct Django configuration is used. For
   testing other Python versions use tox::

    $ ./black.sh
    $ flake8 .
    $ make test
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
=======================

Before you submit a pull request, check that it meets these guidelines:

1. Make sure your pull request is up to date with the ``dev`` branch.
2. The pull request should include tests.
3. Black and flake8 should have been executed without errors using settings
   provided in the repo.
4. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in ``CHANGELOG.rst``.
5. The pull request should work for Python version 3.7, 3.8 and 3.9. Check
   https://github.com/bihealth/sodar-core/actions
   and make sure that the tests pass for supported Python versions.
   From v1.0 onwards SODAR Core no longer supports Python 3.6.

Deploying
=========

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in
``CHANGELOG.rst``). Then run::

$ git tag vX.Y.Z
$ git push
$ git push --tags
$ python setup.py sdist
$ twine upload --repository-url https://test.pypi.org/legacy/ dist/*.tar.gz
$ twine upload dist/*.tar.gz

SODAR Core Changelog
^^^^^^^^^^^^^^^^^^^^

Changelog for the **SODAR Core** Django app package. Loosely follows the
`Keep a Changelog <http://keepachangelog.com/en/1.0.0/>`_ guidelines.

Note that the issue IDs here refer to ones in the private CUBI GitLab.


v0.10.7 (2021-12-14)
====================

Added
-----

- **Adminalerts**
    - UI documentation (#865)
- **Siteinfo**
    - UI documentation (#865)

Changed
-------

- **General**
    - Upgrade minimum Django version to v3.2.10 (#869)
    - Upgrade to python-ldap v3.4.0 (#871)
- **Projectroles**
    - HTTP 403 raised instead of 400 if project type disallowed by API view (#872)
    - Update role list media rules (#863)
    - Add line break for custom email footer (#864)

Fixed
-----

- **Projectroles**
    - ``ManagementCommandLogger`` crash by unset ``LOGGING_LEVEL`` (#862)
    - ``highlight_search_term()`` crash on invalid term input (#867)
    - Search bar allowing invalid input (#868)
    - Wrong project type displayed in project type restriction API response (#872)


v0.10.6 (2021-11-19)
====================

Added
-----

- **General**
    - ``LOGGING_LEVEL`` setting in example configs (#822)
    - ``ProfilingMiddleware`` for cProfile profiling in debug more (#839)
    - ``PROJECTROLES_ENABLE_PROFILING`` setting for profiling (#839)
- **Projectroles**
    - ``cleanup_method`` arg for ``assert_response()`` (#823)
    - Timeline object and data helpers in site and backend plugins (#832)
    - ``ManagementCommandLogger`` helper (#844)
    - ``get_email_user()`` helper (#845)
    - Project type restriction in API views with ``project_type`` attribute (#850)
    - ``Project.has_public_children`` field (#851)
    - Email sending for additional user emails (#861)
    - ``user_email_additional`` app setting (#861)
    - ``email.get_user_addr()`` helper (#861)

Changed
-------

- **General**
    - Upgrade to Chromedriver v96 (#818, #847, #852)
    - Use ``LOGGING_LEVEL`` in example ``set_logging()`` (#822)
    - Upgrade minimum Django version to v3.2.9 (#835, #848)
    - Improve management command output and logging (#844)
    - Optimize project list queries (#851)
- **Filesfolders**
    - Refactor ``checkAll()`` helper (#816)
    - Restrict project type in API views (#850)
- **Projectroles**
    - Upgrade DataTables includes on search results page (#841, #856)
    - Improve email subject prefix formatting (#829)
    - Update user representations in emails (#845)
- **Timeline**
    - Refactor ``TimelineAPI``

Fixed
-----

- **General**
    - Github Actions CI failure by old package version (#821)
    - Codacy code quality badge in README (#815)
- **Appalerts**
    - Random crashes in ``TestTitlebarBadge.test_alert_dismiss_all`` (#811)
- **Projectroles**
    - ``sodar-overflow-container`` failing with certain tables (#830)
    - Sort icons not displayed on search results page (#841)
    - App alert badge content wrapping (#846)
    - Nested categories with public children not displayed correctly for anon users (#853, #855)
    - Public and remote icons breaking project title bar layout (#859)
- **Timeline**
    - Crash from invalid plugin name in ``get_event_description()`` (#831)
    - Redundant database queries in ``get_event_description()`` (#834)
    - Site and backend plugins not supported in ``get_event_description()`` (#832)

Removed
-------

- **Projectroles**
    - ``get_star()`` template tag (#851)
    - ``Project.has_public_children()`` method: use ``has_public_children`` instead (#851)


v0.10.5 (2021-09-20)
====================

Added
-----

- **Appalerts**
    - Display project badge in alert (#790, #801)
    - Dismiss all link in title bar badge (#802)
- **Projectroles**
    - ``exact`` kwarg for ``assert_element_count()`` in UI tests (#798)
    - Custom email header and footer (#789)
    - ``PROJECTROLES_EMAIL_HEADER`` and ``PROJECTROLES_EMAIL_FOOTER`` settings (#789)
    - ``get_all_defs()`` helper in ``AppSettingAPI`` (#808)

Changed
-------

- **General**
    - Unify app settings label notation (#793)
    - Upgrade minimum Django version to v3.2.7 (#800)
- **Appalerts**
    - Improve alert list layout (#790)
- **Projectroles**
    - Improve login button locating in ``login_and_redirect_with_ui()`` (#796)
    - Hide skipped app settings from target remote sync view (#785)
    - Improve app settings layout in target remote sync view (#804)
    - Minor remote sync refactoring (#721, #785, #807)
    - Refactor ``_get_projectroles_settings()`` into ``get_projectroles_defs()`` (#803)

Fixed
-----

- **Appalerts**
    - Redundant HTML anchor in Dismiss All button (#788)
- **Projectroles**
    - Sidebar notch position (#787)
    - ``sodar-overflow-container`` misalignment (#791)
    - App settings recreated if value is identical (#785)
    - Line separators in ``remoteproject_sync.html`` (#805)
    - App settings remote sync only supporting projectroles (#806, #809)
    - Plugin name incorrectly displayed in target remote sync view (#810)
    - Active link check for projectroles URLs ignoring app name (#814)

Removed
-------

- **Projectroles**
    - ``get_plugin_name_by_id()`` template tag (#812)


v0.10.4 (2021-08-19)
====================

Added
-----

- **General**
    - ``LOGGING_APPS`` and ``LOGGING_FILE_PATH`` settings in example site (#762)
    - Siteinfo app to logged apps in base config (#767)
- **Appalerts**
    - "Dismiss All" button in alert list (#770, #781)
    - Update list view with reload link on added alerts (#780)
- **Siteinfo**
    - ``ENABLED_BACKEND_PLUGINS`` in ``CORE_SETTINGS`` (#766)

Changed
-------

- **General**
    - Upgrade to Chromedriver v92 (#772)
    - Upgrade minimum Django version to v3.2.6 (#773)
- **Appalerts**
    - Display no alerts element after clearing list (#779)
- **Projectroles**
    - Refactor view test setup (#769)
- **Siteinfo**
    - UI improvements for empty and unset values

Fixed
-----

- **General**
    - SAML attribute map example in config (#760)
    - Docs layout broken by ``docutils>=0.17`` (#763)
    - Logging level not correctly set for all loggers (#771)
- **Projectroles**
    - HTTP 403 raised instead of 404 in API and UI views if object not found (#774)
    - Incorrect message on ownership transfer email notifications (#778)
    - Project update view loading slowed down by large number of child categories (#765)
- **Siteinfo**
    - Plugin settings not read if ``get_statistics()`` raises exception (#767)
    - List layout broken by empty string values (#768)


v0.10.3 (2021-07-01)
====================

Changed
-------

- **General**
    - Upgrade minimum Django version to v3.2.5 (#744)
    - Upgrade Python dependencies (#744)
- **Userprofile**
    - Hide user update button for non-local users (#748)

Fixed
-----

- **Projectroles**
    - False errors from app settings sync if app not installed on target site (#757)
- **Timeline**
    - Uncaught exceptions in ``get_event_description()`` (#749)
- **Tokens**
    - Expiry date incorrectly displayed in token list (#747)
    - Missing query set ordering in token list (#754)

Removed
-------

- **Tokens**
    - Unused ``admin`` and ``models`` modules


v0.10.2 (2021-06-03)
====================

Changed
-------

- **General**
    - Upgrade to Chromedriver v90 (#731)
    - Rename example site adminalerts URL include (#730)
    - Update documentation screenshots (#734)
    - Reorganize static files in documentation (#734)
    - Rename example ``django-db-file-storage`` URL pattern (#732)
    - Upgrade minimum Django version to v3.2.4 (#727)
    - Upgrade Python dependencies (#727)
    - Reformat with Black v21.5b2
- **Projectroles**
    - Display anonymous icon in titlebar dropdown if not logged in (#726)

Fixed
-----

- **General**
    - Figure aspect ratios in documentation (#735)
- **Projectroles**
    - Unhandled exceptions and missing data in project list extra columns (#733)
    - Project star icon alignment (#736)
    - Project list layout broken by ``FILESFOLDERS_SHOW_LIST_COLUMNS`` setting (#737)
    - Public guest access role not displayed in project list (#739)
- **Timeline**
    - Crash in ``add_event()`` if called with ``AnonymousUser`` (#740)


v0.10.1 (2021-05-06)
====================

Added
-----

- **General**
    - Installation via PyPI (#293)
- **Appalerts**
    - Update alerts in JQuery without page reloading (#701, #723)
    - ``APPALERTS_STATUS_INTERVAL`` setting (#701)

Changed
-------

- **General**
    - Upgrade minimum Django version to v3.2.1 (#696)
    - Upgrade django-debug-toolbar to v3.2.1 (#706)
- **Appalerts**
    - Tweak alert layout (#716)
- **Projectroles**
    - Enforce 3 character minimum limit for terms in multi-term search (#715)
    - Improve remote sync stability

Fixed
-----

- **General**
    - Add ``build/`` to ``.gitignore`` (#707)
    - Invalid operating system qualifier in ``setup.py`` (#708)
- **Projectroles**
    - Uncaught exceptions in app plugin ``search()`` (#713)
    - Broken project icon on search results page (#714)
    - Empty search terms not sanitized (#715)
    - Hardcoded optional ``PROJECTROLES_DISABLE_CATEGORIES`` setting in forms (#719)
    - Remote sync objects referred by database ID instead of ``sodar_uuid`` (#720)
    - Uncaught exceptions in app settings remote sync (#720)
    - Assumed ``sodar_uuid`` match for target app settings in remote sync (#722)


v0.10.0 (2021-04-28)
====================

Added
-----

- **Adminalerts**
    - ``get_statistics()`` implementation
- **Appalerts**
    - Add site app and backend for app alerts (#642)
- **Projectroles**
    - ``geticons`` management command for retrieving Iconify icons (#54)
    - ``spin`` class in ``projectroles.css`` for spinning icon support (#54)
    - Optional public guest access for projects (#574)
    - ``public_guest_access`` and ``set_public()`` in ``Project`` model (#574)
    - Enable allowing anonymous access to site (#574)
    - ``PROJECTROLES_ALLOW_ANONYMOUS`` site setting (#574)
    - ``is_allowed_anonymous`` predicate in ``rules`` (#574)
    - ``site_app_processor`` in ``context_processors`` (#574)
    - ``get_statistics()`` in ``SiteAppPluginPoint``
    - ``info_settings`` in app plugins (#671)
    - ``plugin_type`` argument in ``get_app_plugin()`` (#309)
    - ``handle_project_update()`` in ``ProjectAppPlugin`` (#387, #675)
    - App alerts for project and role updates (#642, #692)
- **Siteinfo**
    - Display selected Django settings in UI (#671)
- **Timeline**
    - Permission tests (#144)
    - Site app plugin for site-wide events (#668)
- **Tokens**
    - Permission tests

Changed
-------

- **General**
    - Upgrade project to Django v3.2 (#194, #695)
    - Upgrade Python dependencies (#194, #678, #685)
    - Rename GitHub repo to ``sodar-core`` (#699)
    - Rename ``master`` branch to ``main``
    - Use Iconify for icons (#54)
    - Use Material Design Icons as default icon set (#54)
    - Bump minimum Python version requirement to v3.7 (#121)
    - Upgraded versioneer (#656)
    - Update views, mixins and tags for anonymous user access (#574)
    - Upgrade recommended development OS version to Ubuntu v20.04 (#640)
    - Do not send redundant emails to users initiating updates (#693)
    - Get all app settings from environment
- **Projectroles**
    - Set parent owner as initial owner in project form (#667)
    - Always show Django admin warning (#677)
    - Modify signature of ``get_history_dropdown()`` template tag (#668)
    - Add default ``superuser`` value to ``LiveUserMixin._make_user()``
    - Include ``select2`` CSS locally (#457)
    - Refactor ``cleanappsettings`` (#673)
- **Siteinfo**
    - Tabbed layout in site info view
- **Timeline**
    - Make ``project`` and ``user`` fields in ``ProjectEvent`` optional (#119, #668)
    - Modify signatures of ``get_object_url()`` and ``get_object_link()`` helpers (#668)
    - Allow custom ``INIT`` status data (#700)
- **Tokens**
    - Refactor view tests

Fixed
-----

- **General**
    - All app settings not properly frozen in test config (#688)
- **Adminalerts**
    - Pagedown widget breaking CSS layout in Firefox (#659)
- **Bgjobs**
    - Plugin queries in template tag module root (#653)
- **Projectroles**
    - Description line spacing in project header (#632)
    - Pagedown widget breaking CSS layout in Firefox (#659)
    - Crash by missing optional ``PROJECTROLES_DELEGATE_LIMIT`` setting (#676)
    - ``cleanappsettings`` deleting defined app settings (#673)
- **Timeline**
    - Double status added when calling ``add_event()`` with ``INIT`` type (#700)

Removed
-------

- **General**
    - Font Awesome support without Iconify (#54)
- **Projectroles**
    - ``get_site_app()`` template tag (#574)
    - Deprecated search functionality with a single ``search_term`` (#618)
    - Deprecated ``get_full_title()`` method from ``Project`` model (#620)


v0.9.1 (2021-03-05)
===================

Added
-----

- **Projectroles**
    - Inline head include from environment variables in base template (#639)
    - ``req_kwargs`` argument in ``SODARAPIPermissionTestMixin.assert_response_api()`` (#662)
    - Display inherited owner note in remote project sync UI (#643)
    - ``is_inherited_owner()`` template tag

Changed
-------

- **General**
    - Improve Codacy support in GitHub Actions
    - Upgrade to Chromedriver v89 (#657)
- **Projectroles**
    - Duplicate ``sodar_uuid`` in ``SODARNestedListSerializer`` (#633)
    - Rename and refactor ``LocalUserForm`` and ``user_form.html`` (#651)

Fixed
-----

- **Filesfolders**
    - File list breadcrumb icon alignment (#660)
    - Cancel link in batch edit view (#647)
    - Batch move folders not displayed in UI (#648)
    - Batch moving objects to project root failing (#661)
- **Projectroles**
    - Login redirect URLs with query strings not properly handled by ``assert_response()`` (#635)
    - Remote project icons in project list not displayed (#664)
    - Version ``0.8.4`` missing from ``CORE_API_ALLOWED_VERSIONS``
- **Userprofile**
    - User update link and template not working as expected (#650)

Removed
-------

- **Userprofile**
    - Unused template ``user_update.html`` (#651)


v0.9.0 (2021-02-03)
===================

Added
-----

- **General**
    - SAML SSO authentication support (#588)
    - REST API example ``HelloExampleProjectAPIView`` in ``example_project_app`` (#518)
- **Projectroles**
    - Projectroles app settings (#532)
    - Remote sync for projectroles app setting (#533, #586)
    - IP address based access restriction for projects (#531)
    - ``is_delegate()`` and ``is_owner_or_delegate()`` helpers for ``Project`` model
    - Remote sync for non-owner category members (#502)
    - ``setting_delete()`` function to ``AppSettingAPI`` (#538)
    - ``cleanappsettings`` management command (#374)
    - ``exclude_inherited`` argument in ``Project.get_delegates()`` (#595)
    - Value options for app settings of type ``STRING`` and ``INTEGER`` (#592)
    - Display placeholders for app setting form fields (#584)
    - Support for local user invites (#548, #613, #615, #621)
    - Local user account creation and updating (#547)
    - ``batchupdateroles`` management command (#15, #602)
    - Project invite REST API views (#15, #598)
    - Advanced search with multiple terms (#609)
    - Search result pagination control (#610)
    - REST API endpoint for retrieving current user info (#626)

Changed
-------

- **General**
    - Replace development helper scripts with ``Makefile`` (#135)
    - Upgrade to Bootstrap v4.5.3 and jQuery v3.5.1 (#563)
    - Upgrade to Chromedriver v87
    - Upgrade general Python requirements (#576)
    - Migrate GitHub CI from Travis to GitHub actions (#577)
    - Refactor example ``PROJECT_USER`` scope app settings (#599)
    - Set logging level in test configurations to ``CRITICAL`` (#604)
- **Filesfolders**
    - Update ``search()`` and ``find()`` for multiple search terms (#609)
- **Projectroles**
    - Allow updating local app settings on a ``TARGET`` site (#545)
    - Refactor project list filtering (#566)
    - Move project list javascript to ``project_list.js`` (#566)
    - Rename owner role transfer URL pattern and timeline event (#590)
    - Add ``sodar_url`` override to ``modify_assignment()``
    - Rename ``ProjectSearchResultsView`` and its template (#609)
    - Implement ``get_full_title()`` as ``Project.full_title`` field (#93)
    - Clarify invite accepting procedure in invite email (#627)
    - Redirect to home view when reusing accepted invite link (#628)
- **Userprofile**
    - Cosmetic updates for user detail template (#600)

Fixed
-----

- **Projectroles**
    - Invite redirect not working in Add Member view (#589)
    - Wrong role label displayed for category owner/delegate in member list (#593)
    - Django settings access in ``forms`` and ``serializers``
    - Delegate limit check broken by existing delegate roles of inherited owners (#595)
    - Crash in project invite if multiple users exist with the same email (#614)
    - Project delegate able to revoke invite for another delegate (#617)
    - Column alignment in invite list (#606)
    - ``get_not_found_alert()`` fails if called with app plugin search type (#624)
- **Taskflowbackend**
    - Django settings access in ``api`` (#605)
    - ``sodar_url`` override not working if ``request`` object is present (#605)

Removed
-------

- **General**
    - Travis CI setup in ``.travis.yml`` (#577)
- **Projectroles**
    - Template ``_project_filter_item.html`` (#566)
    - Template tag ``get_project_list()`` (#566)
    - Deprecate old implementation of ``ProjectAppPluginPoint.search()`` (#609, #618)
    - Deprecate ``Project.get_full_title()`` (#93)


v0.8.4 (2020-11-12)
===================

Changed
-------

- **General**
    - Documentation updates for JOSS submission


v0.8.3 (2020-09-28)
===================

Added
-----

- **General**
    - Missing migration for the ``SODARUser`` model (#581)

Changed
-------

- **General**
    - Upgrade to Chromedriver v85 (#569)
- **Projectroles**
    - Improve project list header legend (#571)
    - Make ``sync_source_data()`` atomic
    - Prevent creation of local projects under remote categories (#583)
- **Siteinfo**
    - Refactor app plugin statistics retrieval (#573)

Fixed
-----

- **General**
    - Invalid statement in ``setup_database.sh`` (#580)
- **Projectroles**
    - Missing exception handling for ``sync_source_data()`` calls (#582)
    - Crash from conflicting local category structure (#582)
- **Siteinfo**
    - Crash from exceptions raised by app plugin ``get_statistics()`` (#572)
- **Timeline**
    - CSS for ``sodar-tl-link-detail`` links (#578)

Removed
-------

- **General**
    - Unused ``Pillow`` dependency (#575)


v0.8.2 (2020-07-22)
===================

Added
-----

- **Bgjobs**
    - Enable site-wide background jobs (#544)
    - Site app plugin for site-wide background jobs (#544)
- **Projectroles**
    - ``sodar-header-button`` CSS class (#550)
    - Logging for ``AppSettingAPI`` (#559)

Changed
-------

- **Projectroles**
    - Upgrade to Chromedriver v83 (#543)
    - Rename ``is_app_link_visible()`` template tag into ``is_app_visible()`` (#546)
    - Refactor project list to reduce queries and template tag use (#551, #567)

Fixed
-----

- **Projectroles**
    - Transferring project ownership to inherited owner not allowed (#534)
    - Uniqueness constraint in ``AppSetting`` incompatible with ``PROJECT_USER`` scope settings (#542)
    - Inherited owner email address not displayed in project member list (#541)
    - App visibility check broken in ``project_detail.html`` (#546)
    - Invite accept for a category invoking Taskflow and causing a crash (#552)
    - Project form ``parent`` forced to wrong value if user lacks role in parent category (#558)
    - Invalid ``app_name`` not handled in ``AppSettingAPI.get_default_setting()`` (#560)
    - Empty JSON and false boolean app settings not set in project form (#557)
    - Minor Javascript errors thrown by ``projectroles.js`` (#536)
    - Long lines breaking email preview layout (#564)


v0.8.1 (2020-04-24)
===================

Added
-----

- **Projectroles**
    - CSS class ``sodar-pr-project-list-custom`` for custom project list items (#525)

Fixed
-----

- **Projectroles**
    - CSS padding issue with ``sodar-list-btn`` and Chrome (#529, sodar#844)
    - Crash from missing optional setting ``PROJECTROLES_DISABLE_CATEGORIES`` (#524)
    - Remote project editing not prevented in REST API views (#523)

Removed
-------

- **Projectroles**
    - Deprecated ``SODARAPIObjectInProjectPermissions`` base class (#527)


v0.8.0 (2020-04-08)
===================

Added
-----

- **General**
    - "For the Impatient" section in docs
- **Filesfolders**
    - API views for file, folder and hyperlink management (#443)
- **Projectroles**
    - Import new REST API view base classes from SODAR (#48, #461)
    - Import base serializers from SODAR (#462)
    - API views for project and role management (#48, #450)
    - ``projectroles.tests.test_views_api.TestAPIViewsBase`` for API view testing (#48)
    - ``SODARAPIPermissionTestMixin`` for API view permission tests
    - New helper methods in ``SODARAPIViewTestMixin``
    - Provide live server URL for Taskflow in ``TestTaskflowBase.request_data`` (#479)
    - ``TestTaskflowAPIBase`` for testing API views with SODAR Taskflow (#488)
    - Permission tests using Knox tokens (#476)
    - Base Ajax view classes in ``projectroles.views_ajax`` (#465)
    - Allow assigning roles for categories (#463)
    - Allow displaying project apps in categories with ``category_enable`` (#447)
    - Allow category delegates and owners to create sub-categories and projects (#464)
    - ``get_role_display_name()`` helper in ``projectroles_common_tags`` (#505)
    - ``get_owners()``, ``is_owner()`` and ``get_all_roles()`` helpers for ``Project`` (#464)
    - Allow using legacy UI test login method with ``PROJECTROLES_TEST_UI_LEGACY_LOGIN`` (#509)
    - Allow moving categories and projects under different categories (#512)
    - ``SODARForm`` and ``SODARModelForm`` base classes for forms
    - Enable retrieving flat recursive list of children objects in ``Project.get_children()``
    - Support for ``data`` in permission test ``assert_response()`` method (#155)
- **Taskflowbackend**
    - ``get_inherited_roles()`` helper (#464)
- **Timeline**
    - ``get_models()`` helper
- **Tokens**
    - Add app from varfish-web (#452)

Changed
-------

- **General**
    - Upgrade minimum Django version to v1.11.29 (#520)
    - Upgrade JQuery to v3.4.1 (#519)
    - Upgrade Bootstrap to v4.4.1 (#460)
    - General upgrade for Python package requirements (#124, #459)
    - Reorganize view classes and URL patterns (#480)
    - Refactor Ajax views (#465, #475)
    - Update ``CONTRIBUTING.rst``
    - Use ``SODARForm`` and ``SODARModelForm`` base classes in forms
- **Projectroles**
    - Suppress peer site removal logging if nothing was removed (#478)
    - Refactor ``SODARCoreAPIBaseView`` into ``SODARCoreAPIBaseMixin`` (#461)
    - Allow providing single user to ``assert_response()`` in permission tests (#474)
    - Move ``SODARAPIViewTestMixin`` into ``test_views_api`` and rename (#471)
    - Move ``KnoxAuthMixin`` functionality into ``SODARAPIViewTestMixin``
    - ``get_accept_header()`` in API tests returns header as dict
    - Refactor base permission test classes (#490)
    - Move ``utils.set_user_group()`` to ``SODARUser.set_group()`` (#483)
    - Call ``set_group()`` in ``SODARUser.save()`` (#483)
    - Replace ``projectroles_tags.is_app_hidden()`` with ``is_app_link_visible()``
    - Inherit owner permissions from parent categories (#464)
    - Refactor project roles template (#505)
    - Disable owner updating in project update form (#508)
    - Allow updating project parent via SODAR Taskflow (#512)
- **Taskflowbackend**
    - Refactor ``synctaskflow`` management command and add logging
- **Timeline**
    - Display app for categories (#447)

Fixed
-----

- **General**
    - Duplicate ``contributing.rst`` redirection file in docs (#481)
    - ``.tox`` not ignored in ``black.sh``
    - Coverage checks in Travis-CI (#507)
- **Projectroles**
    - Swapping owner and delegate roles not allowed if at delegate limit (#477)
    - Remote sync for owner role failing with specific user order in data (#439)
    - Redundant updating of ``Project.submit_status`` during project creation
    - Make ``test_widget_user_options()`` more reliable (#253)
    - Missing permission check by role type in ``RoleAssignmentDeleteView.post()`` (#492)
    - Unordered queryset warnings from the ``User`` model (#494)
    - Incorrect user iteration in ``test_user_autocomplete_ajax()`` (#469)
    - Redundant input validation preventing search with valid characters (#472)
    - Local users disabled in local development configuration (#500)
    - Member link not visible in responsive project dropdown (#466)
    - CSS issues with Bootstrap 4.4.1 in search pagination (#372, #460)
    - Raise ``ImproperlyConfigured`` for missing parameters in ``ProjectAccessMixin`` (#516)
- **Timeline**
    - CSS issues with Bootstrap 4.4.1 (#460)

Removed
-------

- **Projectroles**
    - ``SODARAPIBaseView`` base class, replaced by API view mixins (#461)
    - ``KnoxAuthMixin`` from view tests
    - ``get_selectable_users()`` from ``forms``
    - Redundant render/redirect helpers from ``TestPermissionBase``: use ``assert_response()`` instead (#484)
    - ``APIPermissionMixin`` for API views: use base API/Ajax view classes instead (#467)
    - ``is_app_hidden()`` from ``projectroles_tags``


v0.7.2 (2020-01-31)
===================

Added
-----

- **Projectroles**
    - ``custom_order`` argument in ``get_active_plugins()`` (#431)
    - Enable ordering custom project list columns in project app plugin (#427)
    - ``SODARCoreAPIBaseView`` base API view class for internal SODAR Core apps (#442)
    - API version enforcing in ``RemoteProjectsSyncView`` and ``syncremote.py`` (#444)
    - Allow extra keyword arguments in ``get_backend_api()`` (#397)
    - Example usage of ``get_backend_api()`` extra kwargs in ``example_backend_app`` (#397)
    - ``SODARUserChoiceField`` and ``get_user_widget()`` for user selection in forms (#455)
    - Setting ``reply-to`` headers for role change and invite emails (#446)
    - No reply note and related ``PROJECTROLES_EMAIL_SENDER_REPLY`` setting (#446)
    - Display hidden project app settings to superusers (#424)
- **Sodarcache**
    - Allow limiting ``deletecache`` to a specific project (#448)

Changed
-------

- **General**
    - Upgrade minimum Django version to 1.11.27
    - Base ``RemoteProjectGetAPIView`` on ``SODARCoreAPIBaseView`` (#442)
    - Upgrade to Chromedriver v80 (#510)
- **Bgjobs**
    - Make ``specialize_job()`` more robust (#456)
- **Projectroles**
    - Accept null value for ``AppSetting.value_json`` (#426)
    - Use ``PluginContextMixin`` in ``ProjectContextMixin`` (#430)
    - Move ``get_accept_header()`` to ``SODARAPIViewMixin`` (#445)
    - Allow exceptions to be raised by ``get_backend_plugin()`` (#451)
    - Improve tour help CSS (#438)
    - Field order in ``RoleAssignmentOwnerTransferView`` (#441)
    - Redesign user autocomplete handling in forms (#455)
    - Rename ``SODARUserAutocompleteWidget`` and ``SODARUserRedirectWidget`` (#455)
    - Disable ownership transfer link if owner is the only project user (#454)

Fixed
-----

- **Projectroles**
    - Potential crash in ``_project_header.html`` with ownerless kiosk mode category (#422)
    - Form crash when saving a JSON app setting with ``user_modifiable=False`` (#426)
    - Inconsistent plugin ordering in custom project list columns (#428)
    - Project app plugins included twice in ``HomeView`` (#432)
    - ``ProjectPermissionMixin`` query set override with ``get_project_filter_key()``
    - Search disabled with unchanged input value on search page load (#436)
    - Subprojects queried for non-categories in ``project_detail.html`` (#434)
    - Current owner selectable in ownership transfer form (#440)
- **Taskflowbackend**
    - Potential crash in ``TaskflowAPI`` initialization

Removed
-------

- **Projectroles**
    - Unused backend plugins queried for context data in ``HomeView`` (#433)
    - Unneeded ``UserAutocompleteExcludeMembersAPIView`` (#455)


v0.7.1 (2019-12-18)
===================

Added
-----

- **General**
    - Include CHANGELOG in documentation (#379)
- **Projectroles**
    - ``widget_attrs`` parameter for project and user settings (#404)
    - Remote project member management link for target projects (#382)
    - Current user in ``get_project_list_value()`` arguments (#413)
    - Display category owner in page header (#414)
    - Configuring UI test settings via Django settings or ``TestUIBase`` vars (#417)
    - Initial support for deploying site in kiosk mode (#406)
    - Optional disabling of default CDN Javascript and CSS includes (#418)
    - Defining custom global JS/CSS includes in Django settings (#418)

Changed
-------

- **General**
    - Change "Breaking Changes" doc into "Major Changes" (#201)
    - Refactor and rename ownership transfer classes and template
    - Use RTD theme in documentation (#384)
    - Upgrade to Chromedriver v79
- **Adminalerts**
    - Rename ``INACTIVE`` alert state in UI (#396)
    - Rename URL name and pattern for activation API view (#378)
    - Improve alert detail page layout (#385)
- **Projectroles**
    - Improve unsupported browser warning (#405)
    - Move project list description into tooltip (#388)
- **Siteinfo**
    - Improve page title and heading (#402)
- **Sodarcache**
    - Clarify management command logging (#403)
- **Timeline**
    - Improve extra data status tab legend (#380)

Fixed
-----

- **General**
    - PPA used for Python 3.6 installs no longer available (#416)
- **Filesfolders**
    - Invalid HTML in project list extra columns
- **Projectroles**
    - Dismissing login error alert in ``login.html`` not working (#377)
    - Current owner queries incorrectly filtered in ``RoleAssignmentOwnerTransferView`` (#393)
    - Hardcoded project type display name in sent emails (#398)
    - Silent failing of invalid app setting type in plugin definition (#390)
    - Exception raised by hidden sidebar in sidebar height calculation (#407)
    - Crash in ``get_default_setting()`` if default JSON value was not set (#389)
    - Owner widget hidden in category update view (#394)
    - Project list extra column header alignment not set (#412)
    - ``get_project_list_value()`` template tag displaying "None" on null value (#411)


v0.7.0 (2019-10-09)
===================

Added
-----

- **General**
    - Development env file example ``env.example`` (#297)
    - Postgres database development setup script (#302)
    - ``ENABLE_DEBUG_TOOLBAR`` setting for local development (#349)
    - ``local_target2.py`` config for peer remote site development (#200)
- **Adminalerts**
    - Activate/suspend button in alert list (#42)
- **Bgjobs**
    - Pagination for background job list (#335)
    - ``BGJOBS_PAGINATION`` Django setting (#335)
- **Projectroles**
    - ``get_backend_include()`` common template tag (#261)
    - ``css_url`` member variable in ``BackendPluginPoint`` (#261)
    - Example of on-demand Javascript/CSS inclusion in example apps (#261)
    - Remote project link display toggle for target sites (#276)
    - Project UUID clipboard copying button (#290)
    - Support for app settings in site apps (#308)
    - Initial implemenetation for common clipboard copying visualization (#333)
    - Send email for owner role assignment (#325)
    - Common pagination include template ``_pagination.html`` (#334)
    - Synchronization and display of ``PEER`` sites in remote site management (#200)
    - Link for copying remote site secret token in remote site list (#332)
    - Project ownership transfer from member list (#287)
    - UI notification for disabled member management on target sites (#301)
    - Management command ``addremotesite`` for adding remote sites (#314)
    - JSON support for app settings (#268)
    - ``get_setting_def()`` in app settings API
    - Timeline logging of app settings in project creation (#359)
    - "Project and user" scope for app settings (#266)
    - ``REVOKED`` status for remote projects with revoked access (#327)
    - ``Project.is_revoked()`` helper (#327)
    - Disabling access for non-owner/delegate for revoked projects in ``ProjectPermissionMixin`` (#350)
- **Timeline**
    - Display event extra data as JSON (#6)
- **Userprofile**
    - User setting for project UUID clipboard copying (#290, #308)

Changed
-------

- **General**
    - Upgrade Chromedriver to version 77.0.3865.40
    - Use ``CurrentUserFormMixin`` instead of repeated code (#12)
    - Run tests in parallel where applicable
    - Upgrade minimum Django version to 1.11.25 (#346)
    - General upgrade for Python package requirements (#282)
- **Adminalerts**
    - Use common pagination template
- **Projectroles**
    - Improve user name placeholder in ``login.html`` (#294)
    - Backend app Javascript and CSS included on-demand instead of for all templates (#261)
    - Make sidebar hiding dynamic by content height (#316)
    - Replace ``login_and_redirect()`` in UI tests with a faster cookie based function (#323)
    - Refactor remote project display on details page (#196)
    - Refactor AppSettingAPI (#268)
    - Enable calling ``AppSettingAPI.get_setting_defs()`` with app name instead of plugin object
    - Use ``ProjectPermissionMixin`` on project detail page (#350)
- **Timeline**
    - Use common pagination template (#336)

Fixed
-----

- **Projectroles**
    - Output of template tag ``get_project_link()``
    - Redundant inheritance in ``CurrentUserFormMixin`` (#12)
    - Trailing slashes not parsed correctly in remote project URLs (#319)
    - Crash in ``get_project_column_count()`` with no active project app plugins (#320)
    - UI test helper ``build_selenium_url()`` refactored to work with Chrome v77 (#337)
    - Disallow empty values in ``RemoteSite.name``
    - Remote sync of parent category roles could fail with multiple subprojects
    - ``RemoteProject`` modifications not saved during sync update
    - Timeline events not created in remote project sync (#370)
    - DAL select modifying HTML body width (#365)
    - Footer overflow breaking layout (#367, #375)
- **Timeline**
    - Crash from exception raised by ``get_object_link()`` in a plugin (#328)

Removed
-------

- **Projectroles**
    - Duplicate database indexes from ``RoleAssignment`` (#285)
    - Deprecated ``get_setting()`` tag from ``projectroles_common_tags`` (#283)
    - Project owner change from project updating form (#287)
    - ``ProjectSettingMixin`` from ``projectoles.tests.test_views`` (#357)


v0.6.2 (2019-06-21)
===================

Added
-----

- **General**
    - Badges for Readthedocs documentation and Zenodo DOI (#274)
- **Bgjobs**
    - ``BackgroundJobFactory`` for tests from Varfish-web
- **Projectroles**
    - Unit test to assure owner user creation during project update when using SODAR Taskflow (sodar_taskflow#49)
    - Common template tag ``get_app_setting()`` (#281)
    - Hiding app settings from forms with ``user_modifiable`` (#267)
    - ``AppSetting.value_json`` field (#268)
- **Sodarcache**
    - Logging in ``delete_cache()`` (#279)
- **Userprofile**
    - Support for ``AppSetting.user_modifiable`` (#267)

Changed
-------

- **General**
    - Upgrade minimum Django version to 1.11.21 (#278)
- **Projectroles**
    - ``get_setting()`` template tag renamed into ``get_django_setting()`` (#281)
    - Implement project app descriptions on details page with ``get_info_link()`` (#277)

Fixed
-----

- **General**
    - Documentation sections for Readthedocs


v0.6.1 (2019-06-05)
===================

Added
-----

- **Filesfolders**
    - Example project list columns (#265)
    - Setting ``FILESFOLDERS_SHOW_LIST_COLUMNS`` to manage example project list columns (#265)
- **Projectroles**
    - Optional project list columns for project apps (#265)
- **Sodarcache**
    - ``delete_cache()`` API function (#257)

Changed
-------

- **Projectroles**
    - Refactor ``RemoteProject.get_project()`` (#262)
    - Use ``get_info_link()`` in remote site list (#264)
    - Define ``SYSTEM_USER_GROUP`` in ``SODAR_CONSTANTS`` (#251)
    - Make pagedown textarea element resizeable and increase minimum height (#273)
- **Sodarcache**
    - Handle and log raised exceptions in ``synccache`` management command (#272)
- **Userprofile**
    - Disable user settings link if no settings are available (#260)

Fixed
-----

- **General**
    - Chrome and Chromedriver version mismatch in Travis-CI config (#254)
- **Projectroles**
    - Remove redundant ``get_project_list()`` call from ``project_detail.html``

Removed
-------

- **Projectroles**
    - Unused project statistics in the home view (#269)
    - App settings deprecation protection (#245)
- **Sodarcache**
    - Unused ``TaskflowCacheUpdateAPIView`` (#205)


v0.6.0 (2019-05-10)
===================

Added
-----

- **Filesfolders**
    - Provide app statistics for siteinfo (#18)
- **Projectroles**
    - User settings for settings linked to users instead of projects (#16)
    - ``user_settings`` field in project plugins (#16)
    - Optional ``label`` key for settings
    - Optional "wait for element" args in UI test helpers to ease Javascript testing (#230)
    - ``get_info_link()`` template tag (#239)
    - ``get_setting_defs()`` API function for retrieving project and user setting definitions (#225)
    - ``get_all_defaults()`` API function for retrieving all default setting values (#225)
    - Human readable labels for app settings (#9)
- **Siteinfo**
    - Add app for site info and statistics (#18)
- **Sodarcache**
    - Optional ``--project`` argument for the ``synccache`` command (#232)
- **Timeline**
    - Provide app statistics for siteinfo (#18)
- **Userprofiles**
    - View and form for displaying and updating user settings (#16)

Changed
-------

- **General**
    - Upgrade to ChromeDriver v74 (#221)
- **Bgjobs**
    - Job order to match downstream Varfish
- **Filesfolders**
    - Update app settings (#246)
- **Projectroles**
    - Rename ``project_settings`` module to ``app_settings`` (#225)
    - App settings API updated to support project and user settings (#225)
    - Write an empty dict for ``app_settings`` by default

Fixed
-----

- **Bgjobs**
    - Date formatting in templates (#220)
- **Sodarcache**
    - Crash from ``__repr__()`` if project not set (#223)
    - Broken backend plugin icon (#250)

Removed
-------

- **Timeline**
    - Unused and deprecated project settings (#246)


v0.5.1 (2019-04-16)
===================

Added
-----

- **General**
    - Bgjobs/Celery updates from Kiosc (#175)
    - Default error templates in ``projectroles/error/*.html`` (#210)
- **Projectroles**
    - Optional ``user`` argument in ``ProjectAppPlugin.update_cache()`` (#203)
    - Migration for missing ``RemoteProject`` foreign keys (#197)
- **Sodarcache**
    - API logging (#207)
    - Indexing of identifying fields (#218)

Changed
-------

- **General**
    - Extend ``projectroles/base.html`` for all site app templates, update docs (#217)
    - Use projectroles error templates on the example site (#210)
- **Sodarcache**
    - Make ``user`` field optional in models and API (#204)
    - Rename app configuration into ``SodarcacheConfig`` to follow naming conventions (#202)
    - Rename ``updatecache`` management command to ``synccache`` (#208)

Fixed
-----

- **General**
    - Add missing curl dependency in ``install_os_dependencies.sh`` (#211)
    - Django debug toolbar not displayed when using local configuration (#213)
- **Projectroles**
    - Nested app names not properly returned by ``utils.get_app_names()`` (#206)
    - Forced width set for all Bootstrap modals in ``projectroles.css`` (#209)
    - Long category paths breaking remote project list (#84)
    - Incorrect table rows displayed during project list initialization (#212)
    - Field ``project`` not set for source site ``RemoteProject`` objects (#197)
    - Crash from ``project_base.html`` in site app if not overriding title block (#216)

Removed
-------

- **General**
    - Django debug toolbar workarounds from ``project.css`` and ``project.scss`` (#215)
- **Projectroles**
    - ``PROJECTROLES_ADMIN_OWNER`` deprecation protection: use ``PROJECTROLES_DEFAULT_ADMIN`` (#190)


v0.5.0 (2019-04-03)
===================

Added
-----

- **Projectroles**
    - Warning when using an unsupported browser (#176)
    - Setting ``PROJECTROLES_BROWSER_WARNING`` for unsupported browser warning (#176)
    - Javascript-safe toggle for ``get_setting()`` template tag
    - ID attributes in site containers (#173)
    - Setting ``PROJECTROLES_ALLOW_LOCAL_USERS`` for showing and syncing non-LDAP users (#193)
    - Allow synchronizing existing local target users for remote projects (#192)
    - Allow selecting local users if in local user mode (#192)
    - ``RemoteSite.get_url()`` helper
    - Simple display of links to project on external sites in details page (#182)
- **Sodarcache**
    - Create app (#169)

Changed
-------

- **General**
    - Upgrade to Bootstrap 4.3.1 and Popper 1.14.7 (#181)
- **Projectroles**
    - Improve remote project sync logging (#184, #185)
    - Rename ``PROJECTROLES_ADMIN_OWNER`` into ``PROJECTROLES_DEFAULT_ADMIN`` (#187)
    - Update login template and ``get_login_info()`` to support local user mode (#192)

Fixed
-----

- **Projectroles**
    - Crash in ``get_assignment()`` if called with AnonymousUser (#174)
    - Line breaks in templates breaking ``badge-group`` elements (#180)
    - User autocomplete for users with no group (#199)

Removed
-------

- **General**
    - Deprecated Bootstrap 4 workaround from ``project.js`` (#178)


v0.4.5 (2019-03-06)
===================

Added
-----

- **Projectroles**
    - User autocomplete widgets (#51)
    - Logging in ``syncgroups`` and ``syncremote`` management commands
    - ``PROJECTROLES_DELEGATE_LIMIT`` setting (#21)

Changed
-------

- **General**
    - Upgrade minimum Django version to 1.11.20 (#152)
    - Use user autocomplete in forms in place of standard widget (#51)
- **Filesfolders**
    - Hide parent folder widgets in item creation forms (#159)
- **Projectroles**
    - Enable allowing multiple delegates per project (#21)

Fixed
-----

- **Filesfolders**
    - File upload wiget error not displayed without Bootstrap 4 workarounds (#164)
- **Projectroles**
    - Potential crash in ``syncremote`` if run as Celery job (#160)

Removed
-------

- **General**
    - Old Bootstrap 4 workarounds for django-crispy-forms (#157)


v0.4.4 (2019-02-19)
===================

Changed
-------

- **Projectroles**
    - Modify ``modifyCellOverflow()`` to work with non-table containers (#149)
    - Non-Pagedown form textarea height no longer adjusted automatically (#151)

Fixed
-----

- **Projectroles**
    - Crash in remote project sync caused by typo in ``remoteproject_sync.html`` (#148)
    - Textarea element CSS override breaking layout in third party components (#151)


v0.4.3 (2019-01-31)
===================

Added
-----

- **General**
    - Codacy badge in ``README.rst`` (#140)
- **Projectroles**
    - Category and project display name configuration via ``SODAR_CONSTANTS`` (#141)
    - ``get_display_name()`` utils function and template tag to retrieve ``DISPLAY_NAMES`` (#141)
    - Django admin link warning if taskflowbackend is enabled

Changed
-------

- **General**
    - Use ``get_display_name()`` to display category/project type (#141)
- **Projectroles**
    - Hide immutable fields in forms (#142)
    - Rename Django admin link in user dropdown

Fixed
-----

- **Projectroles**
    - View access control for categories (#143)

Removed
-------

- **General**
    - Redundant ``rules.is_superuser`` predicates from rules (#138)
- **Projectroles**
    - ``get_project_type()`` template tag (use ``get_display_name()`` instead)
    - Unused template ``_roleassignment_import.html``
    - ``PROJECT_TYPE_CHOICES`` from ``SODAR_CONSTANTS``
    - ``force_select_value()`` helper no longer used in forms (#142)


v0.4.2 (2019-01-25)
===================

Added
-----

- **General**
    - Flake8 and Codacy coverage in Travis-CI (#122)
    - Flake8 in GitLab-CI (#127)
- **Projectroles**
    - Automatically pass CSRF token to unsafe Ajax HTTP methods (#116)
    - Queryset filtering in ``ProjectPermissionMixin`` from digestiflow-web (#134)
    - Check for ``get_project_filter_key()`` from digestiflow-web (#134)

Changed
-------

- **General**
    - Upgrade minimum Django version to 1.11.18 (#120)
    - Upgrade Python dependencies (#123)
    - Update .coveragerc
    - Upgrade to Bootstrap 4.2.1 (#23)
    - Upgrade to JQuery 3.3.1 (#23)
    - General code cleanup
    - Code formatting with Black (#133)
- **Filesfolders**
    - Refactor ``BatchEditView`` and ``FileForm.clean()`` (#128)
- **Projectroles**
    - Use ``alert-dismissable`` to dismiss alerts (#13, #130)
    - Update DataTables dependency in ``search.html`` template
    - Refactor ``ProjectModifyMixin`` and ``RemoteProjectAPI`` (#128)
    - Disable ``USE_I18N`` in example site settings (#117)
    - Refactor ``ProjectAccessMixin._get_project()`` into ``get_project()`` (#134)
    - Rename ``BaseAPIView`` into ``SODARAPIBaseView``
- **Timeline**
    - Refactor ``get_event_description()`` (#30, #128)

Fixed
-----

- **General**
    - Django docs references (#131)
- **Projectroles**
    - ``sodar-list-dropdown`` layout broke down with Bootstrap 4.2.1 (#23)
    - ``TASKFLOW_TEST_MODE`` not checked for allowing SODAR Taskflow tests (#126)
    - Typo in ``update_remote`` timeline event description (#129)
    - Textarea height modification (#125)
    - Text wrapping in ``sodar-list-btn`` and ``sodar-list-dropdown`` with Bootstrap 4.2.1 (#132)
- **Taskflowbackend**
    - ``TASKFLOW_TEST_MODE`` not checked for allowing ``cleanup()`` (#126)
    - ``FlowSubmitException`` raised instead of ``CleanupException`` in ``cleanup()``

Removed
-------

- **General**
    - Legacy Python2 ``super()`` calls (#118)
- **Projectroles**
    - Custom alert dismissal script (#13)
- **Example Site App**
    - Example file ``test.py``


v0.4.1 (2019-01-11)
===================

Added
-----

- **General**
    - Travis-CI configuration (#90)
- **Adminalerts**
    - Option to display alert to unauthenticated users with ``require_auth`` (#105)
- **Projectroles**
    - ``TaskflowAPIAuthentication`` for handling Taskflow API auth (#47)
    - Handle ``GET`` requests for Taskflow API views (#47)
    - API version settings ``SODAR_API_ALLOWED_VERSIONS`` and ``SODAR_API_MEDIA_TYPE`` (#111)
    - Site app support in ``change_plugin_status()``
    - ``get_sodar_constants()`` helper (#112)
- **Taskflowbackend**
    - API logging

Changed
-------

- **General**
    - Upgrade minimum Python version requirement to 3.6 (#102)
    - Update and cleanup Gitlab-CI setup (#85)
    - Update Chrome Driver for UI tests
    - Cleanup Chrome setup
    - Enable site message display in login view (#105)
    - Cleanup and refactoring for public GitHub release (#90)
    - Drop support for Ubuntu Jessie and Trusty
    - Update installation utility scripts (#90)
- **Filesfolders**
    - Move inline javascript into ``filesfolders.js``
- **Projectroles**
    - Refactor ``BaseTaskflowAPIView`` (#47)
    - Rename Taskflow specific API views (#104)
    - Unify template tag names in ``projectroles_tags``
    - Change default SODAR API media type into ``application/vnd.bihealth.sodar-core+json`` (#111)
    - Allow importing ``SODAR_CONSTANTS`` into settings for modification (#112)
    - Move ``SODAR_CONSTANTS`` to ``constants.py`` (#112)
- **Timeline**
    - Rename Taskflow specific API views (#104)

Fixed
-----

- **Filesfolders**
    - Overwrite check for zip archive upload if unarchiving was unset (#113)
- **Projectroles**
    - Potential Django crash from auth failure in Taskflow API views
    - Timeline description for updating a remote project
    - Project update with Taskflow failure if description not set (#110)
- **Timeline**
    - ``TaskflowEventStatusSetAPIView`` skipping ``sodar_token`` check (#109)

Removed
-------

- **Filesfolders**
    - Unused dropup app buttons mode in templates (#108)
- **Projectroles**
    - Unused arguments in ``email`` API
    - Unused static file ``shepherd-theme-default.css``
    - Disabled role importing functionality (#61, pending #17)
    - Unused dropup app buttons mode in templates (#108)
- **Timeline**
    - ``ProjectEventStatus.get_timestamp()`` helper


v0.4.0 (2018-12-19)
===================

Added
-----

- **General**
    - ``SODAR_API_DEFAULT_HOST`` setting for server host for API View URLs (sodar#396)
- **Bgjobs**
    - Add app from varfish-web (#95)
- **Filesfolders**
    - Add app from sodar v0.4.0 (#86)
- **Projectroles**
    - Setting ``PROJECTROLES_ENABLE_SEARCH`` (#70)
    - Re-enable "home" link in project breadcrumb (#80)
    - ``get_extra_data_link()`` in ProjectAppPluginPoint for timeline extra data (#6)
    - Allow overriding project class in ProjectAccessMixin
    - Optional disabling of categories and nesting with ``PROJECTROLES_DISABLE_CATEGORIES`` (#87)
    - Optional hiding of apps from project menus using ``PROJECTROLES_HIDE_APP_LINKS`` (#92)
    - Secure SODAR Taskflow API views with ``TASKFLOW_SODAR_SECRET`` (#46)
- **Taskflowbackend**
    - ``test_mode`` flag configured with ``TASKFLOW_TEST_MODE`` in settings (#67)
    - Submit ``sodar_secret`` for securing Taskflow API views (#46)
- **Timeline**
    - Display of extra data using ``{extra-NAME}`` (see documentation) (#6)

Changed
-------

- **General**
    - Improve list button and dropdown styles (#72)
    - Move pagedown CSS overrrides into ``projectroles.css``
    - Reduce default textarea height (#96)
- **Projectroles**
    - Make sidebar resizeable in CSS (#71)
    - Disable search if ``PROJECTROLES_ENABLE_SEARCH`` is set False (#70)
    - Allow appending custom items in project breadcrumb with ``nav_sub_project_extend`` block (#78)
    - Allow replacing project breadcrumb with ``nav_sub_project`` block (#79)
    - Disable remote site access if ``PROJECTROLES_DISABLE_CATEGORIES`` is set (#87), pending #76
    - Disable access to invite views for remote projects (#89)
    - Set "project guest" as the default role for new members (#94)
    - Make noncritical settings variables optional (#14)

Fixed
-----

- **General**
    - Potential inheritance issues in test classes (#74)
    - LDAP dependency script execution (#75)
- **Projectroles**
    - Long words in app names breaking sidebar (#71)
    - Member modification buttons visible for superuser in remote projects (#73)
    - Breadcrumb project detail link display issue in ``base.html`` (#77)
    - "None" string displayed for empty project description (#91)
    - Crash in search from empty project description


v0.3.0 (2018-10-26)
===================

Added
-----

- **General**
    - Test config and script for SODAR Taskflow testing
- **Adminalerts**
    - Add app based on SODAR v0.3.3 (#27)
    - ``TASKFLOW_TARGETS`` setting
- **Projectroles**
    - ``RemoteSite`` and ``RemoteProject`` models (#3)
    - ``RemoteSiteAppPlugin`` site plugin (#3)
    - ``PROJECTROLES_SITE_MODE`` and ``PROJECTROLES_TARGET_CREATE`` settings (#3)
    - Remote site and project management site app (#3)
    - Remote project API (#3)
    - Generic SODAR API base classes
    - ``SodarUserMixin`` for SODAR user helpers in tests
    - Optional ``readme`` and ``sodar_uuid`` args for ``_make_project()`` in tests
    - ``syncremote`` management command for calling ``RemoteProjectAPI.sync_source_data()``
    - ``get_project_by_uuid()`` and ``get_user_by_username()`` template tags
    - ``get_remote_icon()`` template tag (#3)
    - Predicates in rules for handling remote projects (#3)
    - ``ProjectModifyPermissionMixin`` for access control for remote projects (#3)
    - ``is_remote()`` and ``get_source_site()`` helpers in the ``Project`` model (#3)
    - Include template ``_titlebar_nav.html`` for additional title bar links
- **Taskflowbackend**
    - Add app based on SODAR v0.3.3 (#38)
- **Timeline**
    - ``RemoteSite`` model in ``api.get_event_description()`` (#3)

Changed
-------

- **General**
    - Update documentation for v0.3 changes, projectroles usage and fixes to v0.2 docs (#26)
- **Adminalerts**
    - Make ``ADMINALERTS_PAGINATION`` setting optional
- **Projectroles**
    - Allow ``LoggedInPermissionMixin`` to work without a permission object for superusers
    - Enable short/full title selection and remote project icon in ``get_project_link()`` template tag
    - Refactor rules
    - Disable Taskflow API views if Taskflow backend is not enabled (#37)
    - DataTables CSS and JS includes loaded in the search template (#45)
- **Timeline**
    - Minor refactoring of ``api.get_event_description()`` (#30)

Fixed
-----

- **General**
    - Pillow dependency typo in ``requirements/base.txt`` (#33)
    - Login page crash if ``AUTH_LDAP*_DOMAIN_PRINTABLE`` not found (#43)
- **Projectroles**
    - Sidebar create project visible for site apps if URL name was "create" (#36)
    - Enabling LDAP without a secondary backend caused a crash (#39)

Removed
-------

- **General**
    - iRODS specific CSS classes from ``projectroles.css``
    - App content width limit in ``projectroles.css``
    - Domain-specific Login JQuery
    - DataTables CSS and JS includes from base template (#45)


v0.2.1 (2018-09-20)
===================

Changed
-------

- **General**
    - Change ``omics_uuid`` field in all apps' models to ``sodar_uuid`` (sodar#166)
- **Projectroles**
    - Rename abstract ``OmicsUser`` model into ``SODARUser`` (sodar#166)
    - Rename ``OMICS_CONSTANTS`` into ``SODAR_CONSTANTS`` (sodar#166)
    - Rename the ``omics_constant()`` template tag into ``sodar_constant()`` (sodar#166)
    - Rename ``omics_url`` in sodar_taskflow tests to ``sodar_url`` (see sodar_taskflow#36)
    - Rename ``shepherd-theme-omics.css`` to ``shepherd-theme-sodar.css`` (sodar#166)


v0.2.0 (2018-09-19)
===================

Added
-----

- **General**
    - ``example_backend_app`` for a minimal backend app example
    - Backend app usage example in ``example_project_app``
- **Timeline**
    - Add timeline app based on SODAR v0.3.2 (#2)
    - App documentation

Changed
-------

- **General**
    - Update integration documentation (#1)
    - Restructure documentation files and filenames for clarity
- **Timeline**
    - Update CSS classes and overrides
    - Rename list views to ``list_project`` and ``list_objects``
    - Rename list template to ``timeline.html``
    - Refactor ``api.get_event_description()``
    - Make ``TIMELINE_PAGINATION`` optional
    - Improve exception messages in ``api.add_event()``

Fixed
-----

- **Timeline**
    - User model access in ``timeline.api``
    - Misaligned back button (#4)
    - Deprecated CSS in main list
- **Projectroles**
    - Third party apps not correctly recognized in ``get_app_names()``


v0.1.0 (2018-09-12)
===================

Added
-----

- **General**
    - Create app package for Projectroles and other reusable apps based on SODAR release v0.3.1
    - ``example_project_app`` to aid testing and work as a minimal example
    - ``example_site_app`` for demonstrating site apps
    - ``SITE_TITLE`` and ``SITE_INSTANCE_TITLE`` settings
    - ``SITE_PACKAGE`` setting for explicitly declaring site path for code
    - Documentation for integration and development
    - Separate LDAP config in ``install_ldap_dependencies.sh`` and ``requirements/ldap.txt``

- **Projectroles**
    - ``static_file_exists()`` and ``template_exists()`` helpers in common template tags
    - Abstract ``OmicsUser`` model
    - ``get_full_name()`` in abstract OmicsUser model
    - ``auth_backends.py`` file for LDAP backends (sodar#132)
    - Versioneer versioning
    - ``core_version()`` in common template tags
    - Check for footer content in ``include/_footer.html``
    - Example of the site base template in ``projectroles/base_site.html``
    - Example of project footer in ``projectroles/_footer.html``

- **Userprofile**
    - Add site app ``userprofile`` with user details
    - Display user UUID in user profile

Changed
-------

- **Projectroles**
    - Move custom modal into ``projectroles/_modal.html``
    - Check for user.name in user dropdown
    - Move content block structure and sidebar inside ``projectroles/base.html``
    - Move site title bar into optional include template ``projectroles/_site_titlebar.html``
    - Move search form into optional include template ``projectroles/_site_titlebar_search.html``
    - Make title bar dropdown inclueable as ``_site_titlebar_dropdown.html``
    - Title bar CSS and layout tweaks
    - Move ``search.js`` under projectroles
    - Move projectroles specific javascript into ``projectroles.js``
    - Move ``site_version()`` into common template tags
    - Move title bar admin and site app links to user dropdown (sodar#342)
    - Move project specific CSS into optionally includable ``projectroles.css``
    - Refactor and cleanup CSS
    - Move ``set_user_group()`` into ``projectroles.utils``
    - Move ``syncgroups`` management command into projectroles
    - Copy improved multi LDAP backend setup from flowcelltool (sodar#132)
    - Move LDAP authentication backends into projectroles (sodar#132)
    - Move ``login.html`` into projectroles
    - Display ``SITE_INSTANCE_TITLE`` in email instead of a hardcoded string
    - Display the first contact in ``settings.ADMINS`` in email footer
    - Use ``get_full_name()`` in email sending
    - Get site version using ``SITE_PACKAGE``
    - Get LDAP domain names to login template from settings
    - Rename custom CSS classes and HTML IDs from ``omics-*`` into ``sodar-*`` (sodar#166)
    - Move Shepherd theme CSS files into projectroles

Fixed
-----

- **Projectroles**
    - Tests referring to the ``filesfolders`` app not included in this project
    - ``TestHomeView.test_render()`` assumed extra SODAR system user was present (see sodar#367)
    - Tour link setup placing

- **Userprofile**
    - Missing user name if ``name`` field not filled in ``user_detail.html``

Removed
-------

- **Projectroles**
    - Deprecated Javascript variables ``popupWaitHtml`` and ``popupNoFilesHtml``
    - Unused template ``irods_info.html``
.. _code_of_conduct:

Code of Conduct
^^^^^^^^^^^^^^^

Everyone interacting in the SODAR Core project's codebases, issue trackers, chat
rooms, and mailing lists is expected to follow the
`PyPA Code of Conduct <https://www.pypa.io/en/latest/code-of-conduct/>`__.

Credits
=======

* Mikko Nieminen <mikko.nieminen@bih-charite.de>
* Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>
* Oliver Stolpe <oliver.stolpe@bih-charite.de>
* Franziska Schumann <franziska.schumann@fu-berlin.de>
* Hendrik Bomhardt <hendrik.bomhardt@student.hpi.de>
* Tim Garrels <tim.garrels@student.hpi.uni-potsdam.de>
* Raunak Agarwal
SODAR Core
^^^^^^^^^^

.. image:: https://badge.fury.io/py/django-sodar-core.svg
    :target: https://badge.fury.io/py/django-sodar-core

.. image:: https://github.com/bihealth/sodar-core/workflows/build/badge.svg
    :target: https://github.com/bihealth/sodar-core/actions?query=workflow%3ABuild

.. image:: https://app.codacy.com/project/badge/Grade/77c0057a041f4e0c9a0bfc79e9023e04
    :target: https://www.codacy.com/app/bihealth/sodar-core/dashboard

.. image:: https://app.codacy.com/project/badge/Coverage/77c0057a041f4e0c9a0bfc79e9023e04
    :target: https://www.codacy.com/gh/bihealth/sodar-core/dashboard

.. image:: https://img.shields.io/badge/License-MIT-green.svg
    :target: https://opensource.org/licenses/MIT

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/ambv/black

.. image:: https://readthedocs.org/projects/sodar-core/badge/?version=latest
    :target: https://sodar-core.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4269346.svg
    :target: https://doi.org/10.5281/zenodo.4269346

SODAR Core is a framework for Django web application development.

It was conceived to facilitate the creation of scientific data management and
analysis web applications, but can be useful in other contexts as well.
In that it is similar to the CMS or ecommerce frameworks that you can find
`Awesome Django List <https://github.com/wsvincent/awesome-django#content-management-systems>`__ but you will find the components/libraries provided in SODAR Core are more generic and in this reflecting the broader range of applications that we target.


Examples / See It In Action
===========================

SODAR Core is a framework for developing Django web applications. The following
data management and analysis web applications are based on SODAR Core and have
been made available as open source:

- **VarFish** is a web-based tool for the analysis of variants.
  It showcases how to build a complex data warehousing and data analysis web
  appliction using SODAR Core.
  More details are described in the `NAR Web Server Issue publication (doi:10.1093/nar/gkaa241) <https://doi.org/10.1093/nar/gkaa241>`__.
  The source code can be found on `github.com/bihealth/varfish-server <https://github.com/bihealth/varfish-server>`__.
  A demo is available at `varfish-demo.bihealth.org <https://varfish-demo.bihealth.org/login/>`__.
- **Digestiflow** is a web-based data system for the management and
  demultiplexing of Illumina Flow Cells. It further implements various tools for
  sanity checking Illumina sample sheets and quality control (e.g., comparing
  barcode adapter sequence and actual sequence present in the sequencer output).
  You can find out more in our publication in `Bioinformatics (doi:10.1093/bioinformatics/btz850) <https://doi.org/10.1093/bioinformatics/btz850>`__.
  The source code can be found on `github.com/bihealth/digestiflow-server <https://github.com/bihealth/digestiflow-server>`__.
- **Kiosc** is a web application that allows to build scheduler Docker
  containers for "data science" apps and dashboards.
  You can find the source code on `github.com/bihealth/kiosc <https://github.com/bihealth/kiosc>`__.


Quickstart
==========

SODAR Core can only be used from within Django projects. The easiest way to
start out is following the `For the Impatient
<https://sodar-core.readthedocs.io/en/latest/for_the_impatient.html>`__
section in our documentation.


Introduction
============

The SODAR Core repository containes reusable and non-domain-specific apps making
up the base of the SODAR system. These apps can be used for any Django
application which wants to make use of the following features:

- Project-based user access control
- Dynamic app content management
- Advanced project activity logging
- Small file uploading and browsing
- Managing server-side background jobs
- Caching and aggregation of data from external services
- Tracking site information and statistics
- API token management

This repository provides the following installable Django apps:

- **projectroles**: Base app for project access management and
  dynamic app content management. All other apps require the integration of
  projectroles.
- **adminalerts**: Site app for displaying site-wide messages to all users.
- **appalerts**: Site app and backend for raising alerts to users from apps.
- **bgjobs**: Project app for managing background jobs.
- **filesfolders**: Storage and management of small files.
- **siteinfo**: Site app for displaying site information and statistics for
  administrators.
- **sodarcache**: Generic caching and aggregation of data referring to external
  services.
- **taskflowbackend**: Backend app providing an API for the optional
  ``sodar_taskflow`` transaction service.
- **timeline**: Project app for logging and viewing project-related activity.
- **tokens**: Token management for API access.
- **userprofile**: Site app for viewing user profiles.

Also included are resources and examples for developing SODAR compatible apps.


Installation
============

The ``django-sodar-core`` package can be installed into your Django project
from PyPI as follows.

.. code-block:: console

    pip install django-sodar-core==0.10.7

For installing a development version you can point your dependency to a specific
commit ID in GitHub. Note that these versions may not be stable.

.. code-block:: console

    pip install -e git+https://github.com/bihealth/sodar-core.git@xxxxxxxxxxxxx#egg=django-sodar-core

Please note that This package installs a collection of Django apps to be used in
a Django web site project. See
`SODAR Core documentation <https://sodar-core.readthedocs.io/en/latest/?badge=latest>`_
for detailed documentation on use, integration and development.


Repository
==========

The branches of the SODAR Core repository are organized as follows:

- ``main``: The latest stable release.
- ``dev``: The current development branch. Most recent features go here, but
  stability is not guaranteed.

When including SODAR Core in your Django project, it is highly recommended to
freeze your dependency to a specific release tag or commit ID. The project is
under development and breaking changes may be introduced frequently.

Contributions should always be based on the ``dev`` branch.
.. _development:


Development
^^^^^^^^^^^

This document presents instructions and guidelines for developing apps
compatible with the SODAR Core framework, as well as development of the SODAR
Core package itself.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   General <dev_general>
   Project Apps <dev_project_app>
   Site Apps <dev_site_app>
   Backend Apps <dev_backend_app>
   SODAR Core <dev_sodar_core>
.. include:: ../../CONTRIBUTING.rst
.. _app_timeline_install:


Timeline Installation
^^^^^^^^^^^^^^^^^^^^^

This document provides instructions and guidelines for installing the
``timeline`` app to be used with your SODAR Core enabled Django site.

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.


Django Settings
===============

The timeline app is available for your Django site after installing
``django-sodar-core``. Add the app into ``THIRD_PARTY_APPS`` as
follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'timeline.apps.TimelineConfig',
    ]

You also need to add the timeline backend plugin in enabled backend plugins.

.. code-block:: python

    ENABLED_BACKEND_PLUGINS = [
        # ...
        'timeline_backend',
    ]


Optional Settings
=================

To alter default timeline app settings, insert the following **optional**
variables with values of your choosing:

.. code-block:: python

    # Timeline app settings
    TIMELINE_PAGINATION = 15    # Number of events to be shown on one page (int)


URL Configuration
=================

In the Django URL configuration file, add the following line under
``urlpatterns`` to include timeline URLs in your site.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'^timeline/', include('timeline.urls')),
    ]


Migrate Database and Register Plugin
====================================

To migrate the Django database and register the timeline app/backend plugins,
run the following management command:

.. code-block:: console

    $ ./manage.py migrate

In addition to the database migration operation, you should see the following
output:

.. code-block:: console

    Registering Plugin for timeline.plugins.ProjectAppPlugin
    Registering Plugin for timeline.plugins.BackendPlugin
.. _app_timeline_api_django:


Timeline Django API Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This document contains Django API documentation for the ``timeline`` app.
Included are functionalities and classes intended to be used by other
applications.


Backend API
===========

The ``TimelineAPI`` class contains the Timeline backend API. It should be
initialized using the ``Projectroles.plugins.get_backend_api()`` function.

.. autoclass:: timeline.api.TimelineAPI
    :members:

Models
======

.. automodule:: timeline.models
    :members:
.. _app_bgjobs_install:


Bgjobs Installation
^^^^^^^^^^^^^^^^^^^

This document provides instructions and guidelines for installing the
``bgjobs`` app to be used with your SODAR Core enabled Django site.

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.


Django Settings
===============

The bgjobs app is available for your Django site after installing
``django-sodar-core``. Add the app into ``THIRD_PARTY_APPS`` as
follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'bgjobs.apps.BgjobsConfig',
    ]


URL Configuration
=================

In the Django URL configuration file, add the following line under
``urlpatterns`` to include bgjobs URLs in your site.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'^bgjobs/', include('bgjobs.urls')),
    ]


Migrate Database and Register Plugin
====================================

To migrate the Django database and register the bgjobs app and job type plugins,
run the following management command:

.. code-block:: console

    $ ./manage.py migrate

In addition to the database migration operation, you should see the following
output:

.. code-block:: console

    Registering Plugin for bgjobs.plugins.ProjectAppPlugin
    Registering Plugin for bgjobs.plugins.BackgroundJobsPluginPoint


Celery Setup
============

**TODO**
.. _dev_general:


General Development Topics
^^^^^^^^^^^^^^^^^^^^^^^^^^


Guidelines
==========

- Best practices from `Two Scoops <https://www.twoscoopspress.com/>`_
  should be followed where applicable
- To maintain consistency, app packages should be named without delimiting
  characters, e.g. ``projectroles`` and ``userprofile``
- It is recommended to add a *"Projectroles dependency"* comment when directly
  importing e.g. mixins or tags from the ``projectroles`` app
- Hard-coded imports from apps *other than* ``projectroles`` should be avoided
    - Use the plugin structure instead
    - See the ``example_backend_app`` for an example
- Using Bootstrap 4 classes together with SODAR specific overrides and
  extensions provided in ``projectroles.js`` is recommended


Common Helpers
==============

Via the projectroles app, SODAR Core provides optional templates for aiding in
maintaining common functionality and layout. Those are defined here.

App Setting API
---------------

For accessing and modifying app settings for project or site apps, you should
use the ``AppSettingAPI``. Below is an example of invoking the API. For the full
API docs, see :ref:`app_projectroles_api_django`.

.. code-block:: python

    from projectroles.app_settings import AppSettingAPI
    app_settings = AppSettingAPI()
    app_settings.get_app_setting('app_name', 'setting_name', project_object)  # Etc..

Form Base Classes
-----------------

Although not required, it is recommended to use common SODAR Core base classes
with built-in helpers for your Django forms. ``SODARForm`` and
``SODARModelForm`` extend Django's ``Form`` and ``ModelForm`` respectively.
These base classes can be imported from ``projectroles.forms``. Currently they
add logging to ``add_error()`` calls, which helps administrators track form
issues encountered by users. Further improvements are to be added in the future.

Pagination Template
-------------------

A common template for adding navigation for list pagination can be found in
``projectroles/_pagination.html``. This can be included to any Django
``ListView`` template which provides the ``paginate_by`` definition, enabling
pagination. If a smaller layout is desired, the ``pg_small`` argument can be
used. An example can be seen below:

.. code-block:: django

    {% include 'projectroles/_pagination.html' with pg_small=True %}


Management Command Logger
-------------------------

When developing management commands for your apps, you may want to log certain
events while also ensuring relevant output is provided to the administrator
issuing a command. For this SODAR Core provides the
``ManagementCommandLogger`` class. It can be called like the standard Python
logger with shortcut commands such as ``info()``, ``debug()`` etc. If you need
to access the actual Python logger being used, you can access it via
``ManagementCommandLogger.logger``. Example of logger usage can be seen below.

.. code-block:: python

    from projectroles.management.logging import ManagementCommandLogger
    logger = ManagementCommandLogger(__name__)
    logger.info('Testing')

.. note::

    The use of this logger class assumes your site sets up logging simlarly to
    the example site and the SODAR Django Site template, including the use of a
    ``LOGGING_LEVEL`` Django settings variable.


Using Icons
===========

To use icons in your apps, use the ``iconify`` class along with the collection
and icon name into the ``data-icon`` attribute. See
`Iconify <https://docs.iconify.design/implementations/css.html>`_ and
`django-iconify <https://edugit.org/AlekSIS/libs/django-iconify/-/blob/master/README.rst>`_
documentation for further information.

Example:

.. code-block:: HTML

    <i class="iconify" data-icon="mdi:home"></i>

Also make sure to modify the ``icon`` attribute of your app plugins to include
the full ``collection:name`` syntax for Iconify icons.

In certain client side Javascript implementations in which icons are loaded or
replaced dynamically, you may have to refer to these URLs as a direct ``img``
element:

.. code-block:: HTML

    <img src="/icons/mdi/home.svg" />

For modifiers such as color and size when using ``img`` tags,
`see here <https://docs.iconify.design/implementations/css.html>`_.


Testing
=======

SODAR Core provides a range of ready made testing classes and mixins for
different aspects of SODAR app testing, from user permissions to UI testing.
See ``projectroles.tests`` for different base classes.

Test Settings
-------------

SODAR Core provides settings for configuring your UI tests, if using the base
UI test classes found in ``projectroles.tests.test_ui``. Default values for
these settings can be found in ``config/settings/test.py``. The settins are as
follows:

- ``PROJECTROLES_TEST_UI_CHROME_OPTIONS``: Options for Chrome through Selenium.
  Can be used to e.g. enable/disable headless testing mode.
- ``PROJECTROLES_TEST_UI_WINDOW_SIZE``: Custom browser window size.
- ``PROJECTROLES_TEST_UI_WAIT_TIME``: Maximum wait time for UI test operations
- ``PROJECTROLES_TEST_UI_LEGACY_LOGIN``: If set ``True``, use the legacy UI
  login and redirect function for testing with different users. This can be used
  if e.g. issues with cookie-based logins are encountered.

Base Test Classes and Helpers
-----------------------------

For base classes and mixins with useful helpers, see the ``projectroles.tests``
modules. The test cases also provide useful examples on how to set up your own
tests.

.. note::

    For REST API testing, SODAR Core uses separate base test classes for the
    internal SODAR Core API, and the API views implemented in the actual site
    built on SODAR Core. For the API views in your site, make sure to test them
    using e.g. ``TestAPIViewsBase`` and **not** ``TestCoreAPIViewsBase``.


Debugging
=========

Debugging helpers and tips are detailed in this section.

Profiling Middleware
--------------------

SODAR Core provides a cProfile using profiler for tracing back sources of page
loading slowdowns. To enable the profiler middleware, include
``projectroles.middleware.ProfilerMiddleware`` in ``MIDDLEWARE`` under your site
configuration. It is recommended to use a settings variable for this similar to
the example site configuration, where ``PROJECTROLES_ENABLE_PROFILING`` controls
this.

Once enabled, adding the ``?prof`` query string attribute to and URL displays
the profiling information.
.. _app_bgjobs:


Bgjobs App
^^^^^^^^^^

The ``bgjobs`` app allows for the management of project-specific and
asynchronous server-side background jobs.

**TODO:** Docs to be filled out

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   Installation <app_bgjobs_install>
   Usage <app_bgjobs_usage>
.. _app_filesfolders_usage:


Filesfolders Usage
^^^^^^^^^^^^^^^^^^

Usage instructions for the ``filesfolders`` app are detailed in this document.


Filesfolders UI
===============

You can browse and manage files in the app's main view according to your
permissions for each project. The *"File Operations"* menu is used to upload new
files as well as add new folders or links. The menu also contains batch moving
and deletion operations, for which items can be checked using the right hand
side checkboxes.

Updating/deleting operations for single items can be accessed in the dropdown
menus for each item. In the item create/update form, you can also *tag* items
with a choice of icons and stylings to represent the item status.

When uploading a .zip archive, you may choose the *"Extract files from archive"*
option to automatically extract archive files and folders into the filesfolders
app. Note that overwriting of files is not currently allowed.

.. figure:: _static/app_filesfolders/sodar_filesfolders.png
    :align: center
    :figwidth: 100%

    Filesfolders main view


App Settings
============

In the project create/update form, set the boolean setting
``filesfolders.allow_public_links`` true to allow providing public links to
files, for people who can access the site but do not necessarily have a user
account or project rights. Note that public link access still has to be granted
for each file through its create/update form... _app_filesfolders:


Filesfolders App
^^^^^^^^^^^^^^^^

The ``filesfolders`` app enables uploading small files into the Django database
and organizing them in folders. It also permits creating hyperlinks, providing
public links to files and automated unpacking of ZIP archives.

The app is displayed as *"Small Files"* on the SODAR site.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   Installation <app_filesfolders_install>
   Usage <app_filesfolders_usage>
   REST API Documentation <app_filesfolders_api_rest>.. _app_projectroles_usage:


Projectroles Usage
^^^^^^^^^^^^^^^^^^

This document provides instructions for using the ``projectroles`` app which has
been integrated into your Django site.

.. hint::

    Detailed instructions for many pages can be found in an interactive tour by
    clicking the "Help" link in the right side of the top navigation bar.

Before reading this document, be sure to see :ref:`app_projectroles_basics` for
basic concepts regarding the use of this app.


Logging In
==========

Unless anonymous access is specifically set on the site, apart from specific
public or token-enabled views, user login is **mandatory** for using a SODAR
Core based Django site.

One can either log in using a local Django user or, if LDAP/AD is enabled, their
LDAP/AD credentials from a supported site. In the latter case, the user domain
must be appended to the user name in form of ``user@DOMAIN``. Single sign-on
with SAML can also be made available.

.. figure:: _static/app_projectroles/sodar_login.png
    :align: center
    :scale: 75%

    SODAR login form


User Interface
==============

Basics
------

Upon loggin into a SODAR Core based Django site using default templates and CSS,
the general view of your site is split into the following elements:

- **Top navigation bar**: Contains the site logo and title, search element, link
  to advanced search, help link and the user dropdown menu.
- **User dropown menu**: Contains links to user management, admin site and
  site-wide apps the user has access to.
- **Project sidebar**: Shortcuts to project apps and project management pages
- **Project navigation**: Project structure breadcrumb (disabled for site apps)
- **Content**: Actual app content goes in this element
- **Footer**: Optional footer with e.g. site info and version

.. figure:: _static/app_projectroles/sodar_home.png
    :align: center
    :scale: 50%

    Home view

.. figure:: _static/app_projectroles/sodar_user_dropdown.png
    :align: center
    :scale: 75%

    User dropdown

Home View
---------

As content within a SODAR Core based site is by default sorted into projects,
the home view displays a tree view of categories and projects to choose from.
You can filter the list with a search term or restrict display to your starred
projects.

Project Detail View
-------------------

The project detail page dynamically imports elements from installed project
apps, usually showing e.g. overview of latest additions to app data, statistics
and/or shortcuts to app functionalities. Here you can also access project apps
from the project sidebar. For project apps, the sidebar link leads to the app
entry point view as defined in the app plugin.

For each page in a project app which extends the default projectroles template
layout, the **project title bar** is displayed on the top of the page. This
contains the project title and description and a link to "star" the project into
your favourites. Below this, the **project app title bar** with possible
app-specific controls is usually displayed.

.. figure:: _static/app_projectroles/sodar_project_detail.png
    :align: center
    :scale: 50%

    Project detail view


Category and Project Management
===============================

In SODAR based sites, data is split into **categories** and **projects**.
Categories may be freely nested and are used as containers of projects. They
may contain a description and readme, but project apps are disabled for
categories unless explicitly enabled. Projects can not be nested within each
other.

Creating a Top Level Category
-----------------------------

Currently, only users with a superuser status can create a top level category.
This can be done by navigating to the *home view* and clicking the
**Create Category** link. To create a category, a name and owner must be
supplied, along with optional description and/or a readme document. All of these
may be modified later.

.. note::

    Currently, only users already previously logged into the system can be added
    as the owner of a category or project. The ability to invite users not yet
    on the site as owners will be added later.

.. hint::

    When setting up a new site, think about what kind of category and project
    structure makes sense for your team and organization. Moving projects and
    categories under different categories is possible, but is not recommended
    and can currently only be done via the admin view or directly in the Django
    shell.

.. figure:: _static/app_projectroles/sodar_category_create.png
    :align: center
    :scale: 50%

    Category/project creation form

Creating Projects
-----------------

Once navigating into a category, a user with sufficient access will see the
**Create Project or Category** link in the side bar. This opens up a form for
adding a project or a nested category under the current category. The form is
identical to top level category creation, except that you can also choose
between creating a project or a category.

Users with the role of *project contributor* or higher in a category are allowed
to create a project within that category.


Updating Projects
-----------------

An existing project or category can be updated from the
**Update Project/Category** link in the side bar. Again, a similar form as
before will be presented to the user. The owner can not be changed here, but
must be modified in the *Members* view instead. It is possible to move the
current category or project under another category by altering the ``parent``
field. The user who does the updating must have a sufficient role in the
target category or superuser status.

.. figure:: _static/app_projectroles/sodar_project_update.png
    :align: center
    :scale: 50%

    Category/project updating form

.. note::

    For remote project synchronized from another SODAR Core based site, you can
    only edit local application settings in this view.

Public Guest Access
-------------------

From the project create/update view, setting *Public Guest Access* will give
everyone in the system read-only guest access in the project without
specifically granting it.

.. warning::

    If the site has been configured to allow in anonymous users, this will also
    make the project available to anyone who can access the site in your
    network! Please use this feature carefully.

App Settings
------------

Project and site apps may define :term:`app settings<App Settings>`, which can
be either be set with the scope of *project*, *user* or *user within a project*.

Widgets for project specific settings will show up in the project creation and
updating form and can only be modified by users with sufficient project access.
User specific settings will be displayed in the
:ref:`Userpforile app <app_userprofile>`.

By defining the attribute ``user_modifiable=False``, project or user app
settings will not be shown in the respective project/user update views. This is
used e.g. in cases where a project app provides its own UI or updates some
"hidden" setting due to user actions. Superusers will still see these hidden
settings in the Update Project view.

Settings with the scope of user within a project do not currently have a
separate UI of their own. Instead, project apps can produce their own user
specific UIs for this functionality if manual user selection is needed.

.. note::

    Currently, project specific app settings are also enabled for categories but
    do not actually do anything. The behaviour regarding this (remove settings /
    inherit by nested projects / etc) is TBD.

The projectroles app provides the following built-in app settings with the
project scope:

- ``ip_restrict``: Restict project access by an allowed IP list if enabled.
- ``ip_allowlist``: List of allowed IP addresses for project access.

To clean up settings which have been stored in the database but have since
been removed from the plugin app settings definitions, run the following
management command:

.. code-block::

    $ ./manage.py cleanappsettings


Member Management
=================

Project member roles can be viewed and modified through the **Members**
link on the sidebar. Modification requires a sufficient role in the project or
category (owner or delegate) or superuser status.

.. figure:: _static/app_projectroles/sodar_role_list.png
    :align: center
    :scale: 50%

    Project member list view

.. note::

    Owners of categories automatically inherit owner rights to projects
    placed under those categories, starting in SODAR Core v0.8.0. Adding
    separate roles for those users in the inherited projects is not allowed.

.. note::

    At this time, category memberships are not automatically propagated to
    projects created under the category. An inheritance functionality may be
    implemented at a later date.

Adding Members
--------------

There are two ways to add new members to a project or a category:

- **Add Member** is used to add member roles to system users.
- **Invite Member** is used to send email invites to users not yet registered
  in the system.

Addition or modification of users sends an email notification to the user in
question if email sending is enabled on your Django server. The emails can be
previewed in corresponding forms.

.. hint::

    As of SODAR Core v0.4.5, it is also possible to create an invite in the "add
    member" form. Inviting is enabled when inputting an email address not found
    among the system users.

Modifying Members
-----------------

Changing or removing user roles can be done from links next to each role on the
member list. Category or project ownership can be transferred to another user
who currently has a role in the project by using the dropdown next to the owner
role.

Invites
-------

Invites are accepted by the responding user clicking on a link supplied in their
invite email and either logging in to the site with their LDAP/AD credentials or
creating a local user. The latter is only allowed if local users are enabled in
the site's Django settings and the user email domain is not associated with
configured LDAP domains. Invites expire after a certain time and can be reissued
or revoked on the **Project Invites** page.

Batch Member Modifications
--------------------------

Batch member updates can be done either by using REST API views with appropriate
project permissions, or by a site admin using the ``batchupdateroles``
management command. The latter supports multiple projects in one batch. It is
also able to send invites to users who have not yet signed up on the site.


Remote Projects
===============

It is possible to sync project metadata and member roles between multiple SODAR
Core based Django sites. Remote sites and access can be managed in the
**Remote Site Access** site app, found in the user dropdown menu in the top
navigation bar.

Alternatively, remote sites can be created using the following management
command:

.. code-block:: console

    $ ./manage.py addremotesite

In the current implementation, your django site must either be in **source** or
**target** mode. A source site can define one or multiple target sites where
project data can be provided. A target site can define exactly one source site,
from which project data can be retrieved from.

.. note::

    These are arbitrary restrictions which may be relaxed in the future, if use
    cases warrant it.

To enable remote project data reading, you must first set up either a target
or a source site depending on the role of your own SODAR site.

.. figure:: _static/app_projectroles/sodar_remote_sites.png
    :align: center
    :scale: 50%

    Remote site list in source mode

As Source Site
--------------

Navigate to the **Remote Site Access** site app and click on the
*Add Target Site* link. You will be provided with a form for specifying the
remote site. A secret string is generated automatically and you need to provide
this to the administrator of the target site in question for accessing your
site.

Here you also have the option to hide the remote project link from your users.
Users viewing the project on the source site then won't see a link to the target
site. Owners and Superusers will still see the link (greyed out). This is most
commonly used for internal test sites which only needs to be used by admins.

Once created, you can access the list of projects on your site in regards to the
created target site. For each project, you may select an access level, of which
three are currently implemented:

- **No access**: No access on the remote site (default)
- **Read roles**: This allows for the target site to read project metadata *and*
  user roles in order to synchronize project access remotely.
- **Revoked access**: Previously available access which has been revoked. The
  project will still remain in the target site, but only superusers, the project
  owner or the project delegate(s) can acesss it.

.. note::

    The *read roles* access level also provides metadata of the categories above
    the selected project so that the project structure can be maintained.

.. note::

    Only LDAP/AD user roles and local administrator *owner* roles are provided
    to the target site. Other local user roles are ignored.

.. note::

    Access levels for purely checking the existence of the project and only
    reading project metadata (title, description..) without member roles are
    implemented in the data model and backend, but currently disabled in the UI.

Once desired access to specific projects has been granted and confirmed, the
target site will sync the data by sending a request to the source site.

.. figure:: _static/app_projectroles/sodar_remote_projects.png
    :align: center
    :scale: 50%

    Remote project list in source mode

As Target Site
--------------

The source site should be set up as above using the *Set Source Site* link,
using the provided secret string as the access token.

After creating the source site, remote project metadata and member roles (for
which access has been granted) can be accessed using the *Synchronize* link.
Additionaly if the remote Source site is synchronized with multiple Target Sites,
information about those other Target sites will be synchronized as well an displayed
as *Peer Sites*.

Alternatively, the following management command can be used:

.. code-block:: console

    $ ./manage.py syncremote

.. note::

    Creating local projects under a category synchronized from a remote source
    site is **not** allowed from v0.8.3 onwards. For having local projects on a
    target site, you should create and use a local root category.

.. note::

    If a local user is the owner of a synchronized project on the source site,
    the user defined in the ``PROJECTROLES_DEFAULT_ADMIN`` will be given the
    owner role. Hence you **must** have this setting defined if you are
    implementing a SODAR site in target mode.


Search
======

The basic search form is displayed in the top navigation bar if enabled. It
takes one string as a search parameter, followed by optional keyword argument.
At this time, the keyword of ``type`` has been implemented, used to limit the
search to a certain data type as specified in app plugins.

Left to the basic search form is a link to the *Advanced Search* page, where you
can currently search for items using multiple search terms combined with the OR
operator.

Search results are split into results from different apps. For example, entering
``test`` will return all objects from all apps containing this string.
Alternatively, entering ``test type:project`` will provide results from any app
configured to produce results of type *project*. By default, this will result in
the projectroles app listing projects which contain the search string in their
name and/or description.

.. note::

    Additional features such as full-text search and more keywords/operators
    will be defined in the future.


REST API
========

Several SODAR Core functionalities are also available via a HTTP REST API
starting in version 0.8. See :ref:`app_projectroles_api_rest` for instructions
on REST API usage.
.. _user_stories:

User Stories
^^^^^^^^^^^^

This section explains how SODAR Core based web applications are built on a very
high level. We assume that you have read the :ref:`for_the_impatient` section
and are basing your web application on the SODAR Core example site as described
there. Also, we assume that you have intermediate experience with Django and
Python programming.

Please also note that the term *web app* refers to the overall
:term:`Django site<Django Site>`, that is what the user sees and what is
commonly referred to as a "web application" or "dynamic website".

The term :term:`Django app<Django App>` refers to the technical term within
Django development, that is a "Django app Python package".


Flow Cell Data Management
=========================

On a very high level here is how SODAR Core was used to built Digestiflow (cf.
:ref:`for_the_impatient_see_it_in_action`).
You can find the source code in the `Github project <https://github.com/bihealth/digestiflow-server>`__.

The aim is to manage the meta data for sequencers and flow cells.

SODAR Core App Configuration
----------------------------

- Configure your Django site to use the ``projectroles`` SODAR Core Django app.
  Each SODAR Core "project" corresponds to one site and this would allow to have
  multiple groups manage their sequencing meta data in the same web app
  instance.
- Configure the ``timeline`` SODAR Core Django app to provide an audit trail of
  changes to data.
- Configure the ``filesfolders`` SODAR Core Django app to manage small file
  uploads in various places.

Custom Apps
-----------

The following description uses domain-specific language from the high-throughput
sequencing domain. While this makes the section harder to understand to the
layperson, explaining the different terms is out of scope in this manual.

- Write the ``sequencers`` app for management of sequencing machines.
- Write the ``barcodes`` app for management of the barcode adapter sets.
- Write the ``flowcells`` app for managing flow cells and libraries.
  This app depends on ``sequencers`` and ``barcodes`` for these apps' Django
  models.


Variant Analysis
================

On a very high level here is how SODAR Core was used to built VarFish (cf.
:ref:`for_the_impatient_see_it_in_action`).
You can find the source code in the `Github project <https://github.com/bihealth/varfish-server>`__.

The aim is to provide an data analysis web application.

SODAR Core App Configuration
----------------------------

- Configure your Django site to use the ``projectroles`` SODAR Core Django app.
  Each SODAR Core "project" corresponds to one site and this would allow to have
  multiple groups manage their sequencing meta data in the same web app
  instance.
- Configure the ``timeline`` SODAR Core Django app to provide an audit trail of
  changes to data.
- Configure the ``bgjobs`` SODAR Core Django app to manage asynchronous jobs,
  such as long-running queries.

Custom Apps
-----------

The following description uses domain-specific language from the medical
genetics domain. While this makes the section harder to understand to the
layperson, explaining the different terms is out of scope in this manual.

- Write various Django apps for importing background data such as population
  frequencies from the gnomAD project, genes, etc.
- Write the ``variants`` app to implement the actual variant filtration.
  This is the core part of Varfish and provides the different Django models,
  views, and templates to actually perform the variant filtration.
- Write the ``importer`` app to implement efficient bulk import of variants from
  annotated TSV files.
.. _app_filesfolders_install:


Filesfolders Installation
^^^^^^^^^^^^^^^^^^^^^^^^^

This document provides instructions and guidelines for installing the
``filesfolders`` app to be used with your SODAR Core enabled Django site.

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.


Django Settings
===============

The filesfolders app is available for your Django site after installing
``django-sodar-core``. Add the app, along with the prerequisite
``django_db_storage`` app into ``THIRD_PARTY_APPS`` as follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'filesfolders.apps.FilesfoldersConfig',
        'db_file_storage',
    ]

Next set the ``db_file_storage`` app as the default storage app for your site:

.. code-block:: python

    DEFAULT_FILE_STORAGE = 'db_file_storage.storage.DatabaseFileStorage'


Fill out filesfolders app settings to fit your site. The settings variables are
explained below:

* ``FILESFOLDERS_MAX_UPLOAD_SIZE``: Max size for an uploaded file in bytes (int)
* ``FILESFOLDERS_MAX_ARCHIVE_SIZE``: Max size for an archive file to be unpacked
  in bytes (int)
* ``FILESFOLDERS_SERVE_AS_ATTACHMENT``: If true, always serve downloaded files
  as attachment instead of opening them in browser (bool)
* ``FILESFOLDERS_LINK_BAD_REQUEST_MSG``: Message to be displayed for a bad
  public link request (string)

Example of default values:

.. code-block:: python

    # Filesfolders app settings
    FILESFOLDERS_MAX_UPLOAD_SIZE = env.int(
        'FILESFOLDERS_MAX_UPLOAD_SIZE', 10485760)
    FILESFOLDERS_MAX_ARCHIVE_SIZE = env.int(
        'FILESFOLDERS_MAX_ARCHIVE_SIZE', 52428800)
    FILESFOLDERS_SERVE_AS_ATTACHMENT = False
    FILESFOLDERS_LINK_BAD_REQUEST_MSG = 'Invalid request'


URL Configuration
=================

In the Django URL configuration file, add the following lines under
``urlpatterns`` to include filesfolders URLs in your site. The latter line is
required by ``db_file_storage`` and should be obfuscated from actual users.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'^files/', include('filesfolders.urls')),
        url(r'^OBFUSCATED_STRING_HERE/', include('db_file_storage.urls')),
    ]


Migrate Database and Register Plugin
====================================

To migrate the Django database and register the filesfolders app plugin, run the
following management command:

.. code-block:: console

    $ ./manage.py migrate

In addition to the database migration operation, you should see the following
output:

.. code-block:: console

    Registering Plugin for filesfolders.plugins.ProjectAppPlugin
.. _app_sodarcache:


Sodarcache App
^^^^^^^^^^^^^^

The ``sodarcache`` app provides a generic data caching functionality for
a SODAR Core based site. This can be used to e.g. locally cache and aggregate
data referring to external sources in order to speed up commonly repeated
queries to databases other than the local Django PostgreSQL.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   Installation <app_sodarcache_install>
   Usage <app_sodarcache_usage>
   Django API <app_sodarcache_api_django>
.. _app_siteinfo:


Siteinfo App
^^^^^^^^^^^^

The ``siteinfo`` site app enables system administrators and developers to view
site details and statistics gathered from project and backend apps.


Basics
======

The app renders a site which displays information and statistics regarding the
site and installed SODAR apps. Providing app statistics for siteinfo done via
implementing the ``get_statistics()`` function in your app plugins. Currently,
access to the app is limited to site administrators.


Installation
============

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.

Django Settings
---------------

The siteinfo app is available for your Django site after installing
``django-sodar-core``. Add the app into ``THIRD_PARTY_APPS`` as follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'siteinfo.apps.SiteinfoConfig',
    ]


URL Configuration
-----------------

In the Django URL configuration file, add the following line under
``urlpatterns`` to include siteinfo URLs in your site.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'^siteinfo/', include('siteinfo.urls')),
    ]

Migrate Database and Register Plugin
------------------------------------

To migrate the Django database and register the siteinfo site app plugin,
run the following management command:

.. code-block:: console

    $ ./manage.py migrate

In addition to the database migration operation, you should see the following
output:

.. code-block:: console

    Registering Plugin for siteinfo.plugins.SiteAppPlugin


Usage
=====

When logged in as a superuser, you can find the "Site Info" link in your user
dropdown menu in the top right corner of the site. This application is only
available for users with superuser status.

The UI is presented under three tabs:

General
    General site information along with statistics provided by the
    ``get_statistics()`` methods implemented in app plugins.
Apps
    List of installed and enabled project, site and backend app plugins.
Settings
    Django settings for the site. Contains settings from apps as specified in
    the ``info_settings`` member variable in the app plugin.


.. figure:: _static/app_siteinfo/site_info.png
    :align: center
    :scale: 60%

    Siteinfo application with the General tab selected


Providing App Statistics
------------------------

In your project app or backend plugin, implement the ``get_statistics()``
method. It should return a dictionary containing, for each statistics item, a
program friendly key and certain member fields:

``label``
    Human readable label for the statistics item.
``value``
    The value to be rendered.
``url``
    The url to link to from the value for additional information (optional).
``description``
    Additional information (optional).

Example:

.. code-block:: python

    def get_statistics(self):
        return {
            'stat_id': {
                'label': 'Some statistic',
                'value': 9000,
                'url': reverse('home'),
                'description': 'More information here',
            }
        }

Providing Site Settings
-----------------------

The site settings to be presented in the :guilabel:`Settings` tab should be
provided as a list in the ``info_settings`` variable of the app plugin.

Example:

.. code-block:: python

    info_settings = [
        'FILESFOLDERS_LINK_BAD_REQUEST_MSG',
        'FILESFOLDERS_MAX_ARCHIVE_SIZE',
        'FILESFOLDERS_MAX_UPLOAD_SIZE',
        'FILESFOLDERS_SERVE_AS_ATTACHMENT',
        'FILESFOLDERS_SHOW_LIST_COLUMNS',
    ]

.. warning::

    For information security, we recommend against including settings containing
    secret values such as passwords to be displayed in the siteinfo app.
.. _major_changes:


Major Changes
^^^^^^^^^^^^^

This document details highlighted updates and breaking changes in SODAR Core
releases. It is recommended to review these notes whenever upgrading from an
older SODAR Core version. For a complete list of changes in current and previous
releases, see the :ref:`full changelog<changelog>`.


v0.10.7 (2021-12-14)
********************

Release Highlights
==================

- Search bug fixes
- REST API project type restriction fixes
- General bug fixes and minor updates
- Upgraded dependencies

Breaking Changes
================

System Prerequisites
--------------------

The following minimum versions have been bumped:

- ``django>=3.2.10, <3.3``
- ``python-ldap==3.4.0``

API View Invalid Project Type Response
--------------------------------------

If ``project_type`` is set in a REST API view and that view is called with an
disallowed value, the view will return HTTP 403 instead of 400. The cause for
this response is included in the ``detail`` field.


v0.10.6 (2021-11-19)
********************

Release Highlights
==================

- Add additional emails for users
- Add project type restriction for API views
- Add profiling middleware
- Improve management command output
- Improve user representation in email
- Optimize project list queries
- Timeline app bug fixes
- Search results layout fixes
- General bug fixes and minor updates
- Upgraded dependencies

Breaking Changes
================

System Prerequisites
--------------------

The minimum Django version has been bumped to v3.2.9.

Search Results DataTables Upgrade
---------------------------------

DataTables includes on the search results page have been upgraded to version
``bs4/dt-1.11.3/b-2.0.1``. You are advised to review the search results layout
for your own apps to ensure everything looks correct.

Project.has_public_children() Removed
-------------------------------------

The ``Project`` model ``has_public_children()`` helper has been removed. In its
place, you should use the ``Project.has_public_children`` field.


v0.10.5 (2021-09-20)
********************

Release Highlights
==================

- Display project badge in app alerts
- Custom email header and footer
- Fix remote sync of non-projectroles app settings
- Multiple app settings remote sync bug fixes
- General bug fixes and minor updates
- Upgraded dependencies

Breaking Changes
================

System Prerequisites
--------------------

The minimum Django version has been bumped to v3.2.7.

Template Tag Removed
--------------------

The ``get_plugin_name_by_id()`` template tag has been removed from
``projectroles_common_tags``. There should be no reason to query app plugins by
database ID. Please use e.g. the utilities found in ``projectroles.plugins``
instead.


v0.10.4 (2021-08-19)
********************

Release Highlights
==================

- Appalerts list view UI improvements
- Siteinfo app and UI improvements
- Fix API and UI views to return 404 status code if object is not found
- General bug fixes and minor updates
- Upgraded dependencies

Breaking Changes
================

System Prerequisites
--------------------

The minimum Django version has been bumped to v3.2.6.

Base UI and API View 404 Responses
----------------------------------

Base UI and API views have been fixed to correctly return HTTP 404 to authorized
users for resources that are not found. This may affect some test cases which
still operate under the assumption of the views returning 403 instead.


v0.10.3 (2021-07-01)
********************

Release Highlights
==================

- General bug fixes and minor updates
- Upgraded dependencies

Breaking Changes
================

System Prerequisites
--------------------

The minimum Django version has been bumped to v3.2.5.

The following third party Python package requirements have been upgraded:

- ``sphinx-rtd-theme>=0.5.2, <0.6`` (base)
- ``black==21.6b0`` (test)


v0.10.2 (2021-06-03)
********************

Release Highlights
==================

- Project list bug fixes
- General bug fixes and minor updates
- Upgraded dependencies
- Minor changes

Breaking Changes
================

System Prerequisites
--------------------

The minimum Django version has been bumped to v3.2.4.

Third party Python package requirements have been upgraded. See the
``requirements`` directory for up-to-date package versions.


v0.10.1 (2021-05-06)
********************

Release Highlights
==================

- Add JQuery status updating for app alerts
- Make project available in PyPI
- Critical bug fixes for remote sync
- Bug fixes and minor updates

Breaking Changes
================

System Prerequisites
--------------------

The minimum versions of dependencies have been bumped as follows:

- Django: v3.2.1
- Django-debug-toolbar: v3.2.1

Base Template Updated
---------------------

If you are overriding the ``base_site.html`` with your own template and intend
to use the ``appalerts`` app, please add the following snippet into the
``javascript`` block in ``{SITE}/templates/base.html``:

.. code-block:: django

    {% block javascript %}
      {# ... #}
      <!-- App alerts Javascript -->
      {% include 'projectroles/_appalerts_include.html' %}
    {% endblock javascript %}

Remote Sync Bug in v0.9
-----------------------

A bug in remote project sync was recently discovered in SODAR Core v0.9.x and
v0.10.0. The bug has been fixed in this release, but the complete fix requires
for both the ``SOURCE`` and ``TARGET`` sites to be upgraded to v0.10. If you
need to use a site based on SODAR Core v0.9 as a remote sync target, please
upgrade your site to
`this hotfix branch <https://github.com/bihealth/sodar-core/tree/0.9.1/fix-settings-sync>`_.
Note that it is recommended to upgrade all your sites to v0.10 as soon as
possible.


v0.10.0 (2021-04-28)
********************

Release Highlights
==================

- Project upgraded to Django v3.2
- Minimum Python version requirement upgraded to 3.7
- Site icons access via Iconify
- Material Design Icons used as default icon set
- Appalerts app for app-generated user alerts
- Site-wide timeline events
- Timeline events without user
- Allow public guest access to projects for authenticated and anonymous users
- Display Django settings in Site Info app

Breaking Changes
================

System Prerequisites
--------------------

Python version requirements have been upgraded as follows:

- The **minimum** Python version is 3.7
- The **recommended** Python version is 3.8
- CI tests are run on Python 3.7, 3.8 and 3.9
- Support for Python 3.6 has been dropped.

It is recommended to always use the most recent minor version of a Python
release.

Third party Python package requirements have been upgraded. See the
``requirements`` directory for up-to-date package versions.

**Ubuntu 20.04 Focal** is now the recommended OS version for development.

Django v3.2 Upgrade
-------------------

This release updates SODAR Core from Django v1.11 to v3.2. This is a breaking
change which causes many updates and also requires updating several
dependencies.

To upgrade, please update the Django requirement along with your site's other
Python requirements to match ones in ``requirements/*.txt``. See
`Django deprecation documentation <https://docs.djangoproject.com/en/dev/internals/deprecation/>`_
for details about what has been deprecated in Django and which changes are
mandatory.

Common known issues:

- Minimum version of PostgreSQL has been raised to v9.5.
- ``ForeignKey`` fields in models must explicitly declare an ``on_delete``
  argument.
- ``is_authenticated()`` and ``is_anonymous()`` in the user model no longer
  work: use ``is_authenticated`` and ``is_anonymous`` instead.
- Replace imports from ``django.core.urlresolvers`` with ``django.urls``.
- Replace ``django.contrib.postgres.fields.JSONField`` with
  ``django.db.models.JSONField``.
- Add ``DEFAULT_AUTO_FIELD = 'django.db.models.AutoField'`` in
  ``config/settings/base.py`` to get rid of database migration warnings.
- Replace ``{% load staticfiles %}`` with ``{% load static %}``.

In the future, the goal is to keep SODAR Core at the latest stable major version
of Django, except for potential cases in which a critical third party package
has not yet been updated to support a new release.

New Context Processors Required
-------------------------------

The following new context processors are required if you intend to include any
site apps to your projects, or make use of site-wide app alerts, respectively.
To make use of these features, please add the following processors in
``base.py`` under ``TEMPLATES``:

.. code-block:: python

    TEMPLATES = [
        {
            'OPTIONS': {
                'context_processors': {
                    # ...
                    'projectroles.context_processors.site_app_processor',
                    'projectroles.context_processors.app_alerts_processor',
                }
            }
        }

REST API Updates
----------------

The following changes have been made to REST API views:

- ``public_guest_access`` parameter added to project API views.

Site Icons Updated
------------------

Instead of directly including Font Awesome, site icons are now accessed as SVG
using `Iconify <https://iconify.design/>`_. The default icon set has been
changed from Font Awesome to `Material Design Icons <https://materialdesignicons.com>`_.
It is however possible to use other icon sets supported by Iconify for your own
SODAR Core apps.

To make your icons work with SODAR Core v0.10+, you will need to take the
following steps.

First, make sure ``django-iconify`` is installed. Add
``dj_iconify.apps.DjIconifyConfig`` to your Django site settings under
``THIRD_PARTY_APPS`` and ``dj_iconify.urls`` to your site URLs in
``config/urls.py``. See SODAR Core or SODAR Django Site settings for an example.

You will also need to set ``ICONIFY_JSON_ROOT`` in the base Django settings.

.. code-block:: django

    ICONIFY_JSON_ROOT = os.path.join(STATIC_ROOT, 'iconify')

If you are overriding the ``base_site.html`` template, add the following lines
to your base template:

.. code-block:: django

    <script type="text/javascript" src="{% url 'config.js' %}"></script>
    <script type="text/javascript" src="{% static 'projectroles/js/iconify.min.js' %}"></script>

Next, you must download the `Iconify JSON collection files <https://github.com/iconify/collections-json/>`_
required for hosting the icons on your Django server. It is recommended to use
the ``geticons`` management command for this. By default, this downloads the
required ``collections.json`` file along with the ``mdi.json`` file for the MDI
icon collection.

.. code-block:: console

    $ ./manage.py geticons

If you wish to also use other collections than MDI, add them as a list using
the ``-c`` argument. The following example downloads the additional ``carbon``
and ``clarity`` icon sets.

.. code-block:: console

    $ ./manage.py geticons -c carbon clarity

Make sure you run ``collectstatic`` after retrieving the collections for
development.

Before committing your code, it is recommended to update your ``.gitignore``
file with the following lines:

.. code-block::

    */static/iconify/*.json
    */static/iconify/json/*.json

To make the icons in your apps work with this change, you must change the icon
syntax in your Django templates. Use ``iconify`` as the base class of the icon
element. Enter the collection and icon name into the ``data-icon`` attribute.

Example:

.. code-block:: HTML

    <i class="iconify" data-icon="mdi:home"></i>

Also make sure to modify the ``icon`` attribute of your app plugins to include
the full ``collection:name`` syntax for Iconify icons.

You may have to specify icon sizing manually in certain elements. In that
case, use the ``data-height`` and/or ``data-width`` attributes. For spinning
icons, add the ``spin`` class provided in ``projectroles.css``.

Once you have updated all your icons, you can remove the Font Awesome CSS
include from your base template if you are not directly importing it from
``base_site.html``.

In certain client side Javascript implementations in which icons are loaded or
replaced dynamically, you may have to refer to these URLs as a direct ``img``
element:

.. code-block:: HTML

    <img src="/icons/mdi/home.svg" />

For modifiers such as color and size when using ``img`` tags,
`see here <https://docs.iconify.design/implementations/css.html>`_.

Deprecated Features Removed
---------------------------

The following previously deprecated features have been removed in this release:

- ``Project.get_full_title()`` has been removed. Use ``Project.full_title``
  instead.
- Old style search with a single ``search_term`` argument has been removed. Make
  sure your search implementation expects and uses a ``search_terms`` list
  instead.

Timeline API Changes
--------------------

The signatures for ``get_object_url()`` and ``get_object_link()`` helpers have
changed. They now expect the object itself as first argument, followed by an
optional ``Project`` object. The same also applies for
``get_history_dropdown()`` in projectroles common template tags.

Public Guest Access Support
---------------------------

This version adds public guest access support for projects. By setting
``PROJECTROLES_ALLOW_ANONYMOUS`` true, this can be extended to anonymous users.
For your views to properly support anonymous access, please use the override of
``LoginRequiredMixin`` provided in ``projectroles.views`` instead of the
original mixin supplied in Django.

GitHub Repository Updates
-------------------------

The GitHub repository for the project has been renamed from ``sodar_core`` to
``sodar-core``. Otherwise the URL remains the same:
`<https://github.com/bihealth/sodar-core/>`_

GitHub should redirect from the old name indefinitely. However, just to be sure
it is recommend to update your site's dependencies.

Additionally, the former ``master`` branch has been renamed to ``main``.


v0.9.1 (2021-03-05)
*******************

Release Highlights
==================

- Add inline head include from environment variables in base template
- Duplicate object UUIDs in REST API view nested lists

Breaking Changes
================

Base Template Updated
---------------------

The base site template in ``projectroles/base_site.html`` has been updated. If
you have copied the template to your own site's base template to extend it,
please make sure to copy the latest changes to maintain full compatibility. See
diff between templates or search for lines containing ``inline_head_include``.

Duplicate UUIDs in Nested REST API Lists
----------------------------------------

Nested object lists in SODAR Core REST API views are grouped into dictionaries
using each object's ``sodar_uuid``. From this version onwards, the UUID fields
are duplicated within each object as well. While this isn't a breaking change in
itself, if you use ``SODARNestedListSerializer`` it may cause some of your test
cases to fail unless altered.


v0.9.0 (2021-02-03)
*******************

Release Highlights
==================

- Last major update based on Django v1.11
- Enable modifying local app settings in project update form on target sites
- Add projectroles app settings
- Add remote sync for global projectroles app settings
- Add IP address based access restriction for projects
- Add SSO support via SAML
- Add support for local user invites and local user account creation
- Add batch invites and role updates via management command
- Add REST API views for project invite management
- Add advanced search with multiple terms
- Add REST API view for current user info retrieval

Breaking Changes
================

Development Helper Scripts
--------------------------

Development helper scripts (``.sh``) have been replaced by a ``Makefile``.
Get an overview of the available commands via ``make usage``.

System Prerequisites
--------------------

Third party Python package requirements have been upgraded. See the
``requirements`` directory for up-to-date package versions.

The following third party JS/CSS requirements have been updated:

- JQuery v3.5.1
- Bootstrap v4.5.3

.. note::

    This is the last major update of SODAR Core based on and supporting Django
    v1.11, which is now out of long term support. From v0.10 onwards, SODAR Core
    based sites must be implemented on Django v3.x+.

ProjectAppPlugin Search Updates
-------------------------------

The expected signature for ``ProjectAppPluginPoint.search()`` has changed.
Instead of the ``search_term`` string argument, ``search_terms`` is expected.
This argument is a list of strings expected to be combined with ``OR``
operators.

See the ``filesfolders`` app for an example of the new implementation.

In SODAR Core v0.9, the old deprecated implementation still works, but searching
for multiple terms in the "Advanced Search" view will only return results for
the first search term given. This deprecation protection will be removed in the
next major version. Please update the ``search()`` methods in your project app
plugins if you have implemented them.

Project Full Title Field
------------------------

The full title of a project, including the entire category path, can now be
accessed via the ``Project.full_title``. This enables you to use the field
directly in your Django queries and ordering. The value of the field is
auto-populated on ``Project.save()`` and in a database migration accompanied in
this release.

As a result, the ``Project.get_full_title()`` has been deprecated and will be
removed in the next major SODAR Core release. Please refactor your usage of that
helper into referring to ``Project.full_title`` directly.


v0.8.4 (2020-11-12)
*******************

Release Highlights
==================

This release updates documentation for JOSS submission.

Breaking Changes
================

N/A


v0.8.3 (2020-09-28)
*******************

Release Highlights
==================

- Fix issues in remote project synchronization
- Fix crashes in ``siteinfo`` app from exceptions raised by plugins

Breaking Changes
================

Remote Project Sync and Local Categories
----------------------------------------

When working on a ``TARGET`` site, creating local projects under categories
synchronized from a ``SOURCE`` site is no longer allowed. This is done to avoid
synchronization clashes. If you want to enable local projects on your site in
addition to remote ones, you will need to create a local root category for them.

API Changes
-----------

``ProjectCreateAPIView`` now returns status ``403`` if called on a target site
with disabled local projects, instead of ``400`` as before.


v0.8.2 (2020-07-22)
*******************

Release Highlights
==================

- Enable site-wide background jobs
- Critical bug fixes for project member management
- Minor fixes and updates

Breaking Changes
================

N/A


v0.8.1 (2020-04-24)
*******************

Release Highlights
==================

- Fix checking for remote project status in projectroles REST API views
- Miscellaneous bug fixes

Breaking Changes
================

SODARAPIObjectInProjectPermissions Removed
------------------------------------------

The deprecated ``SODARAPIObjectInProjectPermissions`` base class has been
removed from ``projectroles.views_api``. Please base your REST API views to one
of the remaining base classes instead.


v0.8.0 (2020-04-08)
*******************

Release Highlights
==================

- Add API views for the ``projectroles`` and ``filesfolders`` apps
- Add new base view classes and mixins for API/Ajax views
- Import the ``tokens`` API token management app from VarFish
- Allow assigning roles other than owner for categories
- Allow category delegates and owners to create sub-categories and projects
- Allow moving categories and projects under different categories
- Inherit owner permissions from parent categories
- Allow displaying project apps in categories with ``category_enable``
- Reorganization of views in apps

Breaking Changes
================

Owner Permissions Inherited from Categories
-------------------------------------------

Starting in this version of SODAR Core, category owner permissions are
automatically inherited by projects below those categories, as well as possible
subcategories. If this does not fit your use case, it is recommend to reorganize
your project structure and/or give category access to admin users who have
access to all projects anyway.

Projectroles Views Reorganized
------------------------------

Views, base views related mixins for the ``projectroles`` app have been
reorganized in this version. Please review your projectroles imports.

The revised structure is as follows:

- UI views and related mixins **remain** in ``projectroles.views``
- Ajax API view classes were **moved** into ``projectroles.views_ajax``
- REST API view classes **moved** into ``projectroles.views_api``
- Taskflow API view classes **moved** into ``projectroles.views_taskflow``

The same applies to classes and mxins in view tests. See
``projectroles.tests.test_views*`` to update imports in your tests.

Renamed Projectroles View Classes
---------------------------------

In addition to reorganizing classes into different views, certain view classes
intended to be usable by other apps have been renamed. They are listed below.

- ``UserAutocompleteAPIView`` -> ``UserAutocompleteAjaxView``
- ``UserAutocompleteRedirectAPIView`` -> ``UserAutocompleteRedirectAjaxView``

API View Class Changes
----------------------

``SODARAPIBaseView`` and ``APIPermissionMixin`` have been removed. Please use
appropriate classes and mixins found in ``projectroles.views_api`` and
``projectroles.views_ajax`` instead.

Base Test Class and Mixin Changes
---------------------------------

Base test classes and helper mixins in ``projectroles`` have been changed as
detailed below.

- ``SODARAPIViewMixin`` has been moved into ``projectroles.test_views_api`` and
  renamed into ``SODARAPIViewTestMixin``.
- ``KnoxAuthMixin`` has been combined into ``SODARAPIViewTestMixin``.
- ``get_accept_header()`` returns the header as dict instead of a string.
- ``assert_render200_ok()`` and ``assert_redirect()`` have been removed from
  ``TestPermissionBase``. Please use ``assert_response()`` instead.

In addition to the aforementioned changes, certain minor setup details such as
default user rights and may have changed. If you experience unexpected failures
in your tests, please review the SODAR Core base test classes and helper
methods, refactoring your tests where required.

User Group Updating
-------------------

The ``set_user_group()`` helper has been moved from ``projectroles.utils`` into
the ``SODARUser`` model. It is called automatically on ``SODARUser.save()``, so
manual calling of the method is not required for most cases.

System Prerequisites
--------------------

The following third party JS/CSS requirements have been updated:

- JQuery v3.4.1
- Bootstrap v4.4.1
- Popper.js v1.16.0

The minimum supported versions have been upgraded for a number of Python
packages in this release. It is highly recommended to also upgrade these for
your SODAR Core based site. See the ``requirements`` directory for up-to date
dependencies.

The minimum version requirement for Django has been bumped to 1.11.29.

Default Templates Modified
--------------------------

The default template ``base_site.html`` has been modified in this version. If
you override it with your own altered version, please review the difference and
update your templates as appropriate.

SODAR Taskflow v0.4.0 Required
------------------------------

If using SODAR Taskflow, this release requires release v0.4.0 or higher due to
required support for the ``role_update_irods_batch`` flow.

Known Issues
============

- Category roles beyond owner are not synchronized to target sites in remote
  project sync. This was omitted to maintain compatibility in existing APIs in
  this release. The feature is intended to be implemented in SODAR Core v0.9.
- Project/user app settings cannot be set or updated in the project REST API. A
  separate API for this will be developed. Currently the only way to modify
  app settings is via the GUI.


v0.7.2 (2020-01-31)
*******************

Release Highlights
==================

- Enforce API versions in remote project sync
- Separate base API views for SODAR Core API and external SODAR site APIs
- Redesign user autocomplete field
- Set issuing user email to ``reply-to`` header for role and invite emails
- Display hidden project app settings to superusers in project update form
- Allow providing custom keyword arguments for backend plugin ``get_api()``
  through ``get_backend_api()``
- Enable sorting custom project list columns in plugin definition
- Bug fixes for project list columns

Breaking Changes
================

User Autocomplete Field Redesigned
----------------------------------

User autocomplete field for forms with its related widget(s) have been
redesigned with breaking API changes. Please review the :ref:`dev_project_app`
documentation and modify your implementation accordingly.

Remote Project Sync API Version Enforcing
-----------------------------------------

The remote project sync view initiated from a ``TARGET`` site now sends the
version number, making the ``SOURCE`` site enforce allowed API versions in its
request. Hence, when a major breaking change is made on the source site and
version requirements updated, requests from the target site will no longer work
without upgrading to the latest SODAR Core version.

Exceptions Raised by get_backend_api()
--------------------------------------

The ``get_backend_api()`` method for retrieving backend plugin API objects
no longer suppresses potential exceptions raised by API object initialization.
If it is possible for your API object to raise an exception on initialization,
you will need to handle it when calling this method.

System Prerequisites
--------------------

The minimum version requirement for Django has been bumped to 1.11.27.

KnoxAuthMixin in Tests
----------------------

Default API configuration for methods in ``KnoxAuthMixin`` are now set to
internal SODAR Core API values. If you use the mixin in the tests of your site,
please update the arguments in your method calls accordingly. You can also now
supply the `media_type` argument for relevant functions. The
``get_accept_header()`` method has been moved to a separate
``SODARAPIViewMixin`` helper mixin.


v0.7.1 (2019-12-18)
*******************

Release Highlights
==================

- Project list layout and extra column handling improved
- Allow customizing widgets in app settings
- Enable managing global JS/CSS includes in Django settings
- Initial support for deploying site in kiosk mode
- Critical bug fixes for category and project owner management

Breaking Changes
================

Default Templates Modified
--------------------------

The default templates ``base_site.html`` and ``login.html`` have been modified
in this version. If you override them with your own altered versions, please
review the difference and update your templates as appropriate.

User Added to get_project_list_value()
--------------------------------------

The signature of the ``get_project_list_value()`` method implemented by project
app plugins to return data for extra project list columns has changed. The
``user`` argument which provides the current user has been added. If using this
feature, please make sure to update your implementation(s) of the method.

See :ref:`app_projectroles_api_django` to review the API changes.


v0.7.0 (2019-10-09)
*******************

Release Highlights
==================

- Sync peer project information for remote target sites
- Enable revoking access to remote projects
- Allow defining app settings in site apps
- "User in project" scope added into app settings
- Support JSON in app settings
- Project owner management moved to project member views

Breaking Changes
================

System Prerequisites
--------------------

The minimum supported versions have been upgraded for a number of Python
packages in this release. It is highly recommended to also upgrade these for
your SODAR Core based site. See the ``requirements`` directory for up-to date
dependencies.

Backend Javascript Include
--------------------------

The code in ``base.html`` which was including javascript from backend apps to
all templates in projectsroles was removed. Instead, Javascript and CSS
associated to a backend plugin should now be included in app templates as
needed. This is done using the newly introduced ``get_backend_include()``
template tag in ``projectroles_common_tags``.

Deprecated get_setting() Tag Removed
------------------------------------

The deprecated ``get_setting()`` template tag has been removed from
``projectroles_common_tags``. Please use ``get_django_setting()`` in your
templates instead.

ProjectSettingMixin Removed
---------------------------

In ``projectroles.tests.test_views``, the deprecated ``ProjectSettingMixin``
was removed. If you need to populate app settings in your tests, use the
``AppSettingAPI`` instead.

AppSettingAPI get_setting_defs() Signature Changed
--------------------------------------------------

The ``get_settings_defs()`` function in the app settings API now accepts either
a project app plugin or simply the name of the plugin as string. Due to this
change, the signature of the API function including argument order has changed.
Please see the :ref:`API documentation<app_projectroles_api_django>` for more
details and update your function calls accordingly.

Default Footer Styling Changed
------------------------------

The styling of the page footer and the default ``_footer.html`` have changed.
You no longer need an extra ``<div>`` element for the footer content, unless
you need to do styling overrides yourself.


v0.6.2 (2019-06-21)
*******************

Release Highlights
==================

- Allow hiding app settings from UI forms
- Add template tag for retrieving app settings

Breaking Changes
================

System Prerequisites
--------------------

The minimum version requirement for Django has been bumped to 1.11.21.

Template Tag for Django Settings Access Renamed
-----------------------------------------------

The ``get_setting()`` template tag in ``projectroles_common_tags`` has been
renamed into ``get_django_setting()``. In this version the old tag still works,
but this deprecation protection will be removed in the next release. Please
update any references to this tag in your templates.


v0.6.1 (2019-06-05)
*******************

Release Highlights
==================

- Add custom project list columns definable in ProjectAppPlugin
- Add example project list column implementation in the filesfolders app

Breaking Changes
================

App Settings Deprecation Protection Removed
-------------------------------------------

The deprecation protection set up in the previous release has been removed.
Project app plugins are now expected to declare ``app_settings`` in the format
introduced in v0.6.0.


v0.6.0 (2019-05-10)
*******************

Release Highlights
==================

- Add user specific settings
- Refactor project settings into project/user specific app settings
- Add siteinfo app

Breaking Changes
================

App Settings (Formerly Project Settings)
----------------------------------------

The former Project Settings module has been completely overhauled in this
version and requries changes to your app plugins.

The ``projectroles.project_settings`` module has been renamed into
``projectroles.app_settings``. Please update your dependencies accordingly.

Settings must now be defined in ``app_settings``. The format is identical to
the previous ``project_settings`` dictionary, except that a ``scope`` field is
expected for each settings. Currently valid values are "PROJECT" and "USER". It
is recommended to use the related constants from ``SODAR_CONSTANTS``
instead of hard coded strings.

Example of settings:

.. code-block:: python

    #: Project and user settings
    app_settings = {
        'project_bool_setting': {
            'scope': 'PROJECT',
            'type': 'BOOLEAN',
            'default': False,
            'description': 'Example project setting',
        },
        'user_str_setting': {
            'scope': 'USER',
            'type': 'STRING',
            'label': 'String example',
            'default': '',
            'description': 'Example user setting',
        },
    }

.. warning::

    Deprecation protection is place in this version for retrieving settings from
    ``project_settings`` if it has not been changed into ``app_settings`` in
    your project apps. This protection **will be removed** in the next SODAR
    Core release.


v0.5.1 (2019-04-16)
*******************

Release Highlights
==================

- Sodarcache refactoring and improvements for API, models, management and app
  config
- New default error templates

Breaking Changes
================

Site App Templates
------------------

Templates for **site apps** should extend ``projectroles/base.html``. In earlier
versions the documentation erroneously stated ``projectroles/project_base.html``
as the base template to use. Extending that document does work in this version
as long as you override the given template blocks. However, it is not
recommended and may break in the future.

Sodarcache App Changes
----------------------

The following potentially breaking changes have been made to the sodarcache app.

App configuration naming has been changed to
``sodarcache.apps.SodarcacheConfig``. Please update ``config/settings/base.py``
accordingly.

The field ``user`` has been made optional in models and the API.

An optional ``user`` argument has been added to
``ProjectAppPlugin.update_cache()``. Correspondingly, the similar argument in
``ProjectCacheAPI.set_cache_item()`` has been made optional. Please update your
plugin implementations and function calls accordingly.

The ``updatecache`` management command has been renamed to ``synccache``.

Helper get_app_names() Fixed
-----------------------------

The ``projectroles.utils.get_app_names()`` function will now return nested app
names properly instead of omitting everything beyond the topmost module.

Default Admin Setting Deprecation Removed
-----------------------------------------

The ``PROJECTROLES_ADMIN_OWNER`` setting no longer works. Use
``PROJECTROLES_DEFAULT_ADMIN`` instead.


v0.5.0 (2019-04-03)
*******************

Release Highlights
==================

- New sodarcache app for caching and aggregating data from external services
- Local user mode for site UI and remote sync
- Improved display and logging of remote project sync
- Upgrade to Bootstrap 4.3.1

Breaking Changes
================

Default Admin Setting Renamed
-----------------------------

The setting ``PROJECTROLES_ADMIN_OWNER`` has been renamed into
``PROJECTROLES_DEFAULT_ADMIN`` to better reflect its uses. Please rename this
settings variable on your site configuration to prevent issues.

.. note::

    In this release, the old settings value is still accepted in remote project
    management to avoid sudden crashes. This deprecation will be removed in the
    next release.

Bootstrap 4.3.1 Upgrade
-----------------------

The Bootstrap and Popper dependencies have been updated to the latest versions.
Please test your site to make sure this does not result in compatibility issues.
The known issue of HTML content not showing in popovers has already been fixed
in ``projectroles.js``.

Default Templates Modified
--------------------------

The default templates ``base_site.html`` and ``login.html`` have been modified
in this version. If you override them with your own altered versions, please
review the difference and update your templates as appropriate.


v0.4.5 (2019-03-06)
*******************

Release Highlights
==================

- Add user autocomplete in forms
- Allow multiple delegates per project

Breaking Changes
================

System Prerequisites
--------------------

The minimum version requirement for Django has been bumped to 1.11.20.

User Autocomplete Widget Support
--------------------------------

Due to the use of autocomplete widgets for users, the following apps must be
added into ``THIRD_PARTY_APPS`` in ``config/settings/base.py``, regardless of
whether you intend to use them in your own apps:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'dal',
        'dal_select2',
    ]

Project.get_delegate() Helper Renamed
-------------------------------------

As the limit for delegates per project is now arbitrary, the
``Project.get_delegate()`` helper function has been replaced by
``Project.get_delegates()``. The new function returns a ``QuerySet``.

Bootstrap 4 Crispy Forms Overrides Removed
------------------------------------------

Deprecated site-wide Bootstrap 4 theme overrides for ``django-crispy-forms``
were removed from the example site and are no longer supported. These
workarounds were located in ``{SITE_NAME}/templates/bootstrap4/``. Unless
specifically required forms on your site, it is recommended to remove the files
from your project.

.. note::

    If you choose to keep the files or similar workarounds in your site, you
    are responsible of maintaining them and ensuring SODAR compatibility. Such
    site-wide template overrides are outside of the scope for SODAR Core
    components. Leaving the existing files in without maintenance may cause
    undesireable effects in the future.

Database File Upload Widget
---------------------------

Within SODAR Core apps, the only known issue caused by removal of the
aforementioned Bootstrap 4 form overrides in the file upload widget of the
``django-db-file-upload`` package. If you are using the file upload package in
your own SODAR apps and have removed the site-wide Crispy overrides, you can fix
this particular widget by adding the following snippet into your form template.
Make sure to replace ``{FIELD_NAME}`` with the name of your form field.

.. code-block:: django

    {% block css %}
      {{ block.super }}
      {# Workaround for django-db-file-storage Bootstrap4 issue (#164) #}
      <style type="text/css">
        div#div_id_{FIELD_NAME} div p.invalid-feedback {
        display: block;
      }
      </style>
    {% endblock css %}

Alternatively, you can create a common override in your project-wide CSS file.


v0.4.4 (2019-02-19)
*******************

Release Highlights
==================

N/A (maintenance/bugfix release)

Breaking Changes
================

Textarea Height in Forms
------------------------

Due to this feature breaking the layout of certain third party components,
textarea height in forms is no longer adjusted automatically. An exception to
this are Pagedown-specific markdown fields.

To adjust the height of a textarea field in your forms, the easiest way is to
modify the widget of the related field in the ``__init__()`` function of your
form as follows:

.. code-block:: python

    self.fields['field_name'].widget.attrs['rows'] = 4


v0.4.3 (2019-01-31)
*******************

Release Highlights
==================

- Add display name configuration for projects and categories
- Hide immutable fields in projectroles forms

Breaking Changes
================

SODAR Constants
---------------

``PROJECT_TYPE_CHOICES`` has been removed from ``SODAR_CONSTANTS``, as it can
vary depending on implemented ``DISPLAY_NAMES``. If needed, the currently
applicable form structure can be imported from ``projectroles.forms``.


v0.4.2 (2019-01-25)
*******************

Release Highlights
==================

N/A (maintenance/bugfix release)

Breaking Changes
================

System Prerequisites
--------------------

The following minimum version requirements have been upgraded in this release:

- Django 1.11.18+
- Bootstrap 4.2.1
- JQuery 3.3.1
- Numerous required Python packages (see ``requirements/*.txt``)

Please go through your site requirements and update dependencies accordingly.
For project stability, it is still recommended to use exact version numbers for
Python requirements in your SODAR Core based site.

If you are overriding the ``projectroles/base_site.html`` in your site, make
sure to update Javascript and CSS includes accordingly.

.. note::

    Even though the recommended Python version from Django 1.11.17+ is 3.7, we
    only support Python 3.6 for this release. The reason is that some
    dependencies still exhibit problems with the most recent Python release at
    the time of writing.

ProjectAccessMixin
------------------

The ``_get_project()`` function in ``ProjectAccessMixin`` has been renamed into
``get_project()``. Arguments for the function are now optional and may be
removed in a subsequent release: ``self.request`` and ``self.kwargs`` of the
view class will be used if the arguments are not present.

Base API View
-------------

The base SODAR API view has been renamed from ``BaseAPIView`` into
``SODARAPIBaseView``.

Taskflow Backend API
--------------------

The ``cleanup()`` function in ``TaskflowAPI`` now correctly raises a
``CleanupException`` if SODAR Taskflow encounters an error upon calling its
cleanup operation. This change should not affect normally running your site, as
the function in question should only be called during Taskflow testing.


v0.4.1 (2019-01-11)
*******************

Release Highlights
==================

- Configuration updates for API and Projectroles
- Travis-CI setup

Breaking Changes
================

System Prerequisites
--------------------

Changes in system requirements:

- **Ubuntu 16.04 Xenial** is the target OS version.
- **Python 3.6 or newer required**: 3.5 and older releases no longer supported.
- **PostgreSQL 9.6** is the recommended minimum version for the database.

Site Messages in Login Template
-------------------------------

If your site overrides the default login template in
``projectroles/login.html``, make sure your overridden version contains an
include for ``projectroles/_messages.html``. Following the SODAR Core template
conventions, it should be placed as the first element under the
``container-fluid`` div in the ``content`` block. Otherwise, site app messages
not requiring user authorization will not be visible on the login page. Example:

.. code-block:: django

  {% block content %}
    <div class="container-fluid">
      {# Django messages / site app messages #}
      {% include 'projectroles/_messages.html' %}
      {# ... #}
    </div>
  {% endblock content %}


v0.4.0 (2018-12-19)
*******************

Release Highlights
==================

- Add filesfolders app from SODAR v0.4.0
- Add bgjobs app from Varfish-Web
- Secure SODAR Taskflow API views
- Separate test server configuration for SODAR Taskflow
- Extra data variable rendering for timeline
- Additional site settings

Breaking Changes
================

List Button Classes in Templates
--------------------------------

Custom small button and dropdown classes for including buttons within tables and
lists have been modified. The naming has also been unified. The following
classes should now be used:

- Button group: ``sodar-list-btn-group`` (formerly ``sodar-edit-button-group``)
- Button: ``sodar-list-btn``
- Dropdown: ``sodar-list-dropdown`` (formerly ``sodar-edit-dropdown``)

See projectroles templates for examples.

.. warning::

    The standard bootstrap class ``btn-sm`` should **not** be used with these
    custom classes!

SODAR Taskflow v0.3.1 Required
------------------------------

If using SODAR Taskflow, this release requires release v0.3.1 or higher due to
mandatory support of the ``TASKFLOW_SODAR_SECRET`` setting.

Taskflow Secret String
----------------------

If you are using the ``taskflow`` backend app, you **must** set the value of
``TASKFLOW_SODAR_SECRET`` in your Django settings. Note that this must match the
similarly named setting in your SODAR Taskflow instance!


v0.3.0 (2018-10-26)
*******************

Release Highlights
==================

- Add remote project metadata and member synchronization between multiple SODAR
  sites
- Add adminalerts app
- Add taskflowbackend app

Breaking Changes
================

Remote Site Setup
-----------------

For specifying the role of your site in remote project metadata synchronization,
you will need to add two new settings to your Django site configuration:

The ``PROJECTROLES_SITE_MODE`` setting sets the role of your site in remote
project sync and it is **mandatory**. Accepted values are ``SOURCE`` and
``TARGET``. For deployment, it is recommended to fetch this setting from
environment variables.

If your site is set in ``TARGET`` mode, the boolean setting
``PROJECTROLES_TARGET_CREATE`` must also be included to control whether
creation of local projects is allowed. If your site is in ``SOURCE`` mode, this
setting can be included but will have no effect.

Furthermore, if your site is in ``TARGET`` mode you must include the
``PROJECTROLES_ADMIN_OWNER`` setting, which must point to an existing local
superuser account on your site.

Example for a ``SOURCE`` site:

.. code-block:: python

    # Projectroles app settings
    PROJECTROLES_SITE_MODE = env.str('PROJECTROLES_SITE_MODE', 'SOURCE')

Example for a ``TARGET`` site:

.. code-block:: python

    # Projectroles app settings
    PROJECTROLES_SITE_MODE = env.str('PROJECTROLES_SITE_MODE', 'TARGET')
    PROJECTROLES_TARGET_CREATE = env.bool('PROJECTROLES_TARGET_CREATE', True)
    PROJECTROLES_ADMIN_OWNER = env.str('PROJECTROLES_ADMIN_OWNER', 'admin')

General API Settings
--------------------

Add the following lines to your configuration to enable the general API
settings:

.. code-block:: python

    SODAR_API_DEFAULT_VERSION = '0.1'
    SODAR_API_MEDIA_TYPE = 'application/vnd.bihealth.sodar+json'

DataTables Includes
-------------------

Includes for the DataTables Javascript library are no longer included in
templates by default. If you want to use DataTables, include the required CSS
and Javascript in relevant templates. See the ``projectroles/search.html``
template for an example.
.. _app_projectroles_settings:


Projectroles Django Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This document describes the :term:`Django settings<Django Settings>` for the
``projectroles`` app, which also control the configuration of other apps in a
SODAR Core based site.

These settings are usually found in ``config/settings/*.py``, with
``config/settings/base.py`` being the default configuration other files may
override or extend.

If your site is based on ``sodar-django-site``, mandatory settings are already
set to their default values. In that case, you only need to modify or customize
them where applicable.

If you are integrating django-sodar-core with an existing Django site or
building your site from scratch without the recommended template, make sure to
add all mandatory settings into your project.

For values retrieved from environment variables, make sure to configure your
env accordingly. For development and testing, it is highly recommended to set
``DJANGO_READ_DOT_ENV_FILE=1`` in your system's environment variables and
place the env variables into a ``.env`` file in the root directory of your
Django site repository. See ``env.example`` for an example of such a file.


Site Package and Paths
======================

The site package and path configuration should be found at the beginning of the
default configuration file. Substitute {SITE_NAME} with the name of your site
package.

.. code-block:: python

    import environ
    SITE_PACKAGE = '{SITE_NAME}'
    ROOT_DIR = environ.Path(__file__) - 3
    APPS_DIR = ROOT_DIR.path(SITE_PACKAGE)


Apps
====

Apps installed from django-sodar-core are placed in ``THIRD_PARTY_APPS``. The
following apps need to be included in the list in order for SODAR Core to work:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'crispy_forms',
        'rules.apps.AutodiscoverRulesConfig',
        'djangoplugins',
        'pagedown',
        'markupfield',
        'rest_framework',
        'knox',
        'projectroles.apps.ProjectrolesConfig',
        'dal',
        'dal_select2',
        'dj_iconify.apps.DjIconifyConfig',
    ]


Database
========

Under ``DATABASES``, the setting below is recommended:

.. code-block:: python

    DATABASES['default']['ATOMIC_REQUESTS'] = False

.. note::

    If this conflicts with your existing set up, you can modify the code in your
    other apps to use e.g. ``@transaction.atomic``.

.. note::

    This setting mostly is used for the ``sodar_taskflow`` transactions
    supported by projectroles but not commonly used, so having this setting as
    True *may* cause no issues. However, it is not officially supported at this
    time.


Templates
=========

Under ``TEMPLATES['OPTIONS']['context_processors']``, add the required
projectroles processors:

.. code-block:: python

    'projectroles.context_processors.urls_processor',
    'projectroles.context_processors.site_app_processor',


Email
=====

Under ``EMAIL_CONFIGURATION`` or ``EMAIL``, configure email settings:

.. code-block:: python

    EMAIL_SENDER = env('EMAIL_SENDER', default='noreply@example.com')
    EMAIL_SUBJECT_PREFIX = env('EMAIL_SUBJECT_PREFIX', default='')


Authentication
==============

``AUTHENTICATION_BACKENDS`` should contain the following backend classes:

.. code-block:: python

    AUTHENTICATION_BACKENDS = [
        'rules.permissions.ObjectPermissionBackend',
        'django.contrib.auth.backends.ModelBackend',
    ]

.. note::

    The default setup by cookiecutter-django adds the ``allauth`` package. This
    can be left out of the project if not needed, as it mostly provides adapters
    for e.g. social media account logins. If removing allauth, you can also
    remove unused settings variables starting with ``ACCOUNT_*``.

The following settings remain in your auth configuration:

.. code-block:: python

    AUTH_USER_MODEL = 'users.User'
    LOGIN_REDIRECT_URL = 'home'
    LOGIN_URL = 'login'


Icons
=====

The ``ICONIFY_JSON_ROOT`` setting must point to the appropriate path within
your static files directory in order to make icons work on your SODAR Core based
site.

.. code-block:: python

    ICONIFY_JSON_ROOT = os.path.join(STATIC_ROOT, 'iconify')


Django REST Framework
=====================

To enable ``djangorestframework`` API views and ``knox`` authentication, these
values should be added under ``DEFAULT_AUTHENTICATION_CLASSES``:

.. code-block:: python

    REST_FRAMEWORK = {
        'DEFAULT_AUTHENTICATION_CLASSES': (
            'rest_framework.authentication.BasicAuthentication',
            'rest_framework.authentication.SessionAuthentication',
            'knox.auth.TokenAuthentication',
        ),
    }


General Site Settings
=====================

For display in projectroles based templates, customize related variables to
describe your site. ``SITE_INSTANCE_TITLE`` may be used to e.g. differentiate
between site versions used for deployment or staging, use in different
organizations, etc.

.. code-block:: python

    SITE_TITLE = 'Name of Your Project'
    SITE_SUBTITLE = env.str('SITE_SUBTITLE', 'Beta')
    SITE_INSTANCE_TITLE = env.str('SITE_INSTANCE_TITLE', 'Deployment Instance Name')


Projectroles Settings
=====================

**Mandatory** projectroles app settings are explained below:

* ``PROJECTROLES_SITE_MODE``: Site mode for remote project metadata
  synchronization, either ``SOURCE`` (allow others to read local projects) or
  ``TARGET`` (read projects from another site)
* ``PROJECTROLES_TARGET_CREATE``: Whether or not local projects can be created
  if site is in ``TARGET`` mode. If your site is in ``SOURCE`` mode, this
  setting has no effect.
* ``PROJECTROLES_INVITE_EXPIRY_DAYS``: Days until project email invites expire
  (int)
* ``PROJECTROLES_SEND_EMAIL``: Enable/disable email sending (bool)
* ``PROJECTROLES_EMAIL_SENDER_REPLY``: Whether replies are expected to the
  sender address (bool). If set ``False`` and nothing is set in the ``reply-to``
  header, a "do not reply" note is added to the email body.
* ``PROJECTROLES_ENABLE_SEARCH``: Whether you want to enable SODAR search on
  your site (boolean)
* ``PROJECTROLES_DEFAULT_ADMIN``: User name of the default superuser account
  used in e.g. replacing an unavailable user or performing backend admin
  commands (string)

Example:

.. code-block:: python

    # Projectroles app settings
    PROJECTROLES_SITE_MODE = env.str('PROJECTROLES_SITE_MODE', 'TARGET')
    PROJECTROLES_TARGET_CREATE = env.bool('PROJECTROLES_TARGET_CREATE', True)
    PROJECTROLES_INVITE_EXPIRY_DAYS = env.int('PROJECTROLES_INVITE_EXPIRY_DAYS', 14)
    PROJECTROLES_SEND_EMAIL = env.bool('PROJECTROLES_SEND_EMAIL', False)
    PROJECTROLES_EMAIL_SENDER_REPLY = env.bool('PROJECTROLES_EMAIL_SENDER_REPLY', False)
    PROJECTROLES_ENABLE_SEARCH = True
    PROJECTROLES_DEFAULT_ADMIN = env.str('PROJECTROLES_DEFAULT_ADMIN', 'admin')


Optional Projectroles Settings
==============================

The following projectroles settings are **optional**:

* ``PROJECTROLES_EMAIL_HEADER``: Custom email header (string)
* ``PROJECTROLES_EMAIL_FOOTER``: Custom email footer (string)
* ``PROJECTROLES_SECRET_LENGTH``: Character length of secret token used in
  projectroles (int)
* ``PROJECTROLES_SEARCH_PAGINATION``: Amount of search results per each app to
  display on one page (int)
* ``PROJECTROLES_HELP_HIGHLIGHT_DAYS``: Days for highlighting tour help for new
  users (int)
* ``PROJECTROLES_DISABLE_CATEGORIES``: If set True, disable categories and only
  allow a list of projects on the root level (boolean) (see note)
* ``PROJECTROLES_HIDE_APP_LINKS``: Apps hidden from the project sidebar and
  dropdown menus for non-superusers. The app views and URLs are still
  accessible. The names should correspond to the ``name`` property in each
  project app's plugin (list)
* ``PROJECTROLES_DELEGATE_LIMIT``: The number of delegate roles allowed per
  project. The amount is limited to 1 per project if not set, unlimited if set
  to 0. Will be ignored for remote projects synchronized from a source site
  (int)
* ``PROJECTROLES_BROWSER_WARNING``: If true, display a warning to users using
  Internet Explorer (bool)
* ``PROJECTROLES_ALLOW_LOCAL_USERS``: If true, roles for local non-LDAP users
  can be synchronized from a source during remote project sync if they exist on
  the target site. Similarly, local users will be selectable in member dropdowns
  when selecting users (bool)
* ``PROJECTROLES_KIOSK_MODE``: If true, allow accessing certain project views
  *without* user authentication in order to e.g. demonstrate features in a
  kiosk-style deployment. Also hides and/or disables views not intended to be
  used in this mode (bool)
* ``PROJECTROLES_ALLOW_ANONYMOUS``: If true, allow anonymous users to access the
  site and all projects where ``public_guest_access`` is set true (bool)

Example:

.. code-block:: python

    # Projectroles app settings
    # ...
    PROJECTROLES_EMAIL_HEADER = 'This email has been sent by X from Y'
    PROJECTROLES_EMAIL_FOOTER = 'For assistance contact admin@example.com'
    PROJECTROLES_SECRET_LENGTH = 32
    PROJECTROLES_SEARCH_PAGINATION = 5
    PROJECTROLES_HELP_HIGHLIGHT_DAYS = 7
    PROJECTROLES_DISABLE_CATEGORIES = True
    PROJECTROLES_HIDE_APP_LINKS = ['filesfolders']
    PROJECTROLES_DELEGATE_LIMIT = 1
    PROJECTROLES_BROWSER_WARNING = True
    PROJECTROLES_ALLOW_LOCAL_USERS = True
    PROJECTROLES_KIOSK_MODE = False

.. warning::

    Regarding ``PROJECTROLES_DISABLE_CATEGORIES``: In the current SODAR core
    version remote site access and remote project synchronization are disabled
    if this option is used! Use only if a simple project list is specifically
    required in your site.

.. warning::

    Regarding ``PROJECTROLES_ALLOW_LOCAL_USERS``: Please note that this will
    allow synchronizing project roles to local non-LDAP users based on their
    **user name**. You should personally ensure that the users in question are
    authorized for these roles. Furthermore, only roles for **existing** local
    users will be synchronized. New local users will have to be added manually
    through the Django admin or shell on the target site.

.. warning::

    The ``PROJECTROLES_KIOSK_MODE`` setting is under development and considered
    experimental. More implementation, testing and documentation is forthcoming.


Backend App Settings
====================

The ``ENABLED_BACKEND_PLUGINS`` settings lists backend plugins implemented using
``BackendPluginPoint`` which are enabled in the configuration. For more
information see :ref:`dev_backend_app`.

.. code-block:: python

    ENABLED_BACKEND_PLUGINS = env.list('ENABLED_BACKEND_PLUGINS', None, [])


API View Settings (Optional)
============================

If you want to build an API to your site using SODAR Core functionality, it is
recommended to base your API views on ``projectroles.views.SODARAPIBaseView``.
Using this base class also allows you to define your API media type, version
number and allowed versions via Django settings.

The recommended API setup uses accept header versioning. The
``SODAR_API_MEDIA_TYPE`` setting should be changed to your organization and API
identification if API views are introduced. The ``SODAR_API_DEFAULT_HOST``
setting should post to the externally visible host of your server and be
configured in your environment settings.

These settings are **optional**. Default values will be used if they are unset.

Example:

.. code-block:: python

    SODAR_API_DEFAULT_VERSION = '0.1'
    SODAR_API_ACCEPTED_VERSIONS = [SODAR_API_DEFAULT_VERSION]
    SODAR_API_MEDIA_TYPE = 'application/your.application+json'  # Change this
    SODAR_API_DEFAULT_HOST = SODAR_API_DEFAULT_HOST = env.url('SODAR_API_DEFAULT_HOST', 'http://0.0.0.0:8000')


LDAP/AD Configuration (Optional)
================================

If you want to utilize LDAP/AD user logins as configured by projectroles, you
can add the following configuration. Make sure to also add the related env
variables to your configuration.

This part of the setup is **optional**.

.. note::

    In order to support LDAP, make sure you have installed the dependencies from
    ``utility/install_ldap_dependencies.sh`` and ``requirements/ldap.txt``! For
    more information see :ref:`dev_sodar_core`.

.. note::

    If only using one LDAP/AD server, you can leave the "secondary LDAP server"
    values unset.

.. code-block:: python

    ENABLE_LDAP = env.bool('ENABLE_LDAP', False)
    ENABLE_LDAP_SECONDARY = env.bool('ENABLE_LDAP_SECONDARY', False)

    if ENABLE_LDAP:
        import itertools
        import ldap
        from django_auth_ldap.config import LDAPSearch

        # Default values
        LDAP_DEFAULT_CONN_OPTIONS = {ldap.OPT_REFERRALS: 0}
        LDAP_DEFAULT_FILTERSTR = '(sAMAccountName=%(user)s)'
        LDAP_DEFAULT_ATTR_MAP = {
            'first_name': 'givenName',
            'last_name': 'sn',
            'email': 'mail',
        }

        # Primary LDAP server
        AUTH_LDAP_SERVER_URI = env.str('AUTH_LDAP_SERVER_URI', None)
        AUTH_LDAP_BIND_DN = env.str('AUTH_LDAP_BIND_DN', None)
        AUTH_LDAP_BIND_PASSWORD = env.str('AUTH_LDAP_BIND_PASSWORD', None)
        AUTH_LDAP_CONNECTION_OPTIONS = LDAP_DEFAULT_CONN_OPTIONS

        AUTH_LDAP_USER_SEARCH = LDAPSearch(
            env.str('AUTH_LDAP_USER_SEARCH_BASE', None),
            ldap.SCOPE_SUBTREE,
            LDAP_DEFAULT_FILTERSTR,
        )
        AUTH_LDAP_USER_ATTR_MAP = LDAP_DEFAULT_ATTR_MAP
        AUTH_LDAP_USERNAME_DOMAIN = env.str('AUTH_LDAP_USERNAME_DOMAIN', None)
        AUTH_LDAP_DOMAIN_PRINTABLE = env.str(
            'AUTH_LDAP_DOMAIN_PRINTABLE', AUTH_LDAP_USERNAME_DOMAIN
        )

        AUTHENTICATION_BACKENDS = tuple(
            itertools.chain(
                ('projectroles.auth_backends.PrimaryLDAPBackend',),
                AUTHENTICATION_BACKENDS,
            )
        )

        # Secondary LDAP server (optional)
        if ENABLE_LDAP_SECONDARY:
            AUTH_LDAP2_SERVER_URI = env.str('AUTH_LDAP2_SERVER_URI', None)
            AUTH_LDAP2_BIND_DN = env.str('AUTH_LDAP2_BIND_DN', None)
            AUTH_LDAP2_BIND_PASSWORD = env.str('AUTH_LDAP2_BIND_PASSWORD', None)
            AUTH_LDAP2_CONNECTION_OPTIONS = LDAP_DEFAULT_CONN_OPTIONS

            AUTH_LDAP2_USER_SEARCH = LDAPSearch(
                env.str('AUTH_LDAP2_USER_SEARCH_BASE', None),
                ldap.SCOPE_SUBTREE,
                LDAP_DEFAULT_FILTERSTR,
            )
            AUTH_LDAP2_USER_ATTR_MAP = LDAP_DEFAULT_ATTR_MAP
            AUTH_LDAP2_USERNAME_DOMAIN = env.str('AUTH_LDAP2_USERNAME_DOMAIN')
            AUTH_LDAP2_DOMAIN_PRINTABLE = env.str(
                'AUTH_LDAP2_DOMAIN_PRINTABLE', AUTH_LDAP2_USERNAME_DOMAIN
            )

            AUTHENTICATION_BACKENDS = tuple(
                itertools.chain(
                    ('projectroles.auth_backends.SecondaryLDAPBackend',),
                    AUTHENTICATION_BACKENDS,
                )
            )


SAML SSO Configuration (Optional)
=================================

Optional Single Sign-On (SSO) authorization via SAML is also available. To
enable this feature, set ``ENABLE_SAML=1`` in your environment. Configuring SAML
for SSO requires proper configuration of the Keycloak SSO server and the SAML
client library.

Keycloak
--------

Create a new client in Keycloak and configure it as follows. Please note that
**Client ID** can be chosen however you like, but it must match the setting
in the client.

.. figure:: _static/saml/keycloak_client_config.png

To generate the ``metadata.xml`` file required for the client, go to the
**Realm Settings** page and in the **General** tab, click
``SAML 2.0 Identity Provider Metadata`` to download the xml data. Save it
somewhere on the client, the preferred name is ``metadata.xml``.

.. figure:: _static/saml/keycloak_metadata_download.png

For the signing of the request send to the Keycloak server you will require a
certificate and key provided by the Keycloak server and incorporated into the
configuration of the client. Switch to the ``SAML Keys``. Make sure to select
``PKCS12`` as **Archive Format**.

.. figure:: _static/saml/keycloak_saml_key_download1.png
.. figure:: _static/saml/keycloak_saml_key_download2.png

Convert the archive on the commandline with the follow command and store them in
some place on your client.

.. code::

    openssl pkcs12 -in keystore.p12 -password "pass:<PASSWORD>" -nodes | openssl x509 -out cert.pem
    openssl pkcs12 -in keystore.p12 -password "pass:<PASSWORD>" -nodes -nocerts | openssl rsa -out key.pem

SODAR Core
----------

Make sure that your ``config/settings/base.py`` contains the following
configuration:

.. code-block:: python

    ENABLE_SAML = env.bool('ENABLE_SAML', False)
    SAML2_AUTH = {
        # Required setting
        # Pysaml2 Saml client settings
        # See: https://pysaml2.readthedocs.io/en/latest/howto/config.html
        'SAML_CLIENT_SETTINGS': {
            # Optional entity ID string to be passed in the 'Issuer' element of
            # authn request, if required by the IDP.
            'entityid': env.str('SAML_CLIENT_ENTITY_ID', 'SODARcore'),
            'entitybaseurl': env.str(
                'SAML_CLIENT_ENTITY_URL', 'https://localhost:8000'
            ),
            # The auto(dynamic) metadata configuration URL of SAML2
            'metadata': {
                'local': [
                    env.str('SAML_CLIENT_METADATA_FILE', 'metadata.xml'),
                ],
            },
            'service': {
                'sp': {
                    'idp': env.str(
                        'SAML_CLIENT_IPD',
                        'https://sso.hpc.bihealth.org/auth/realms/cubi',
                    ),
                    # Keycloak expects client signature
                    'authn_requests_signed': 'true',
                    # Enforce POST binding which is required by keycloak
                    'binding': 'urn:oasis:names:tc:SAML:2.0:bindings:HTTP-POST',
                },
            },
            'key_file': env.str('SAML_CLIENT_KEY_FILE', 'key.pem'),
            'cert_file': env.str('SAML_CLIENT_CERT_FILE', 'cert.pem'),
            'xmlsec_binary': env.str('SAML_CLIENT_XMLSEC1', '/usr/bin/xmlsec1'),
            'encryption_keypairs': [
                {
                    'key_file': env.str('SAML_CLIENT_KEY_FILE', 'key.pem'),
                    'cert_file': env.str('SAML_CLIENT_CERT_FILE', 'cert.pem'),
                }
            ],
        },
        # Custom target redirect URL after the user get logged in.
        # Defaults to /admin if not set. This setting will be overwritten if you
        # have parameter ?next= specificed in the login URL.
        'DEFAULT_NEXT_URL': '/',
        # # Optional settings below
        # 'NEW_USER_PROFILE': {
        #     'USER_GROUPS': [],  # The default group name when a new user logs in
        #     'ACTIVE_STATUS': True,  # The default active status for new users
        #     'STAFF_STATUS': True,  # The staff status for new users
        #     'SUPERUSER_STATUS': False,  # The superuser status for new users
        # },
        # 'ATTRIBUTES_MAP': env.dict(
        #     'SAML_ATTRIBUTES_MAP',
        #     default={
        #         Change values to corresponding SAML2 userprofile attributes.
        #         'email': 'Email',
        #         'username': 'UserName',
        #         'first_name': 'FirstName',
        #         'last_name': 'LastName',
        #     }
        # ),
        # 'TRIGGER': {
        #     'FIND_USER': 'path.to.your.find.user.hook.method',
        #     'NEW_USER': 'path.to.your.new.user.hook.method',
        #     'CREATE_USER': 'path.to.your.create.user.hook.method',
        #     'BEFORE_LOGIN': 'path.to.your.login.hook.method',
        # },
        # Custom URL to validate incoming SAML requests against
        # 'ASSERTION_URL': 'https://your.url.here',
    }

Add the following settings to your environment variables:

.. code-block::

    ENABLE_SAML=1
    SAML_CLIENT_ENTITY_ID=<Entity ID configured in Keycloak>
    SAML_CLIENT_ENTITY_URL=<Client URL, e.g. https://sodar-core.bihealth.org>
    SAML_CLIENT_METADATA_FILE=<e.g. metadata.xml>
    SAML_CLIENT_IPO=<SSO server URL, e.g. https://sso.hpc.bihealth.org/auth/realms/cubi>
    SAML_CLIENT_KEY_FILE=<e.g. key.pem>
    SAML_CLIENT_CERT_FILE=<e.g. cert.pem>
    SAML_CLIENT_XMLSEC1=<e.g. /usr/bin/xmlsec1>


Global JS/CSS Include Modifications (Optional)
==============================================

It is possible to supplement (or replace, see below) global Javascript and CSS
includes of your SODAR Core site without altering the base template. You can
place a list of custom includes into the list variables
``PROJECTROLES_CUSTOM_JS_INCLUDES`` and ``PROJECTROLES_CUSTOM_CSS_INCLUDES``.
These can either be local static file paths or web URLs to e.g. CDN served
files.

If using the default CDN imports for JQuery, Bootstrap4 etc. are not an optimal
solution in your use case due to e.g. network issues, you can disable these
includes by setting ``PROJECTROLES_DISABLE_CDN_INCLUDES`` to ``True``.

.. warning::

    If disabling the default CDN includes, you **must** provide replacements for
    **all** disabled files in your custom includes. Otherwise your SODAR Core
    based site will not function correctly!

Example:

.. code-block:: python

    PROJECTROLES_DISABLE_CDN_INCLUDES = True
    PROJECTROLES_CUSTOM_JS_INCLUDES = [
        STATIC_ROOT + '/your/path/jquery-3.3.1.min.js',
        STATIC_ROOT + '/your/path/popper.min.js',
        'https://some-cdn.com/bootstrap.min.js',
        # ...
    ]
    PROJECTROLES_CUSTOM_CSS_INCLUDES = [
        STATIC_ROOT + '/your/path/bootstrap.min.css',
        # ...
    ]

It is also possible to define inline HTML in an environment variable and include
it in the ``head`` tag of the base template. To use this feature, add HTML
script as the value of the variable ``PROJECTROLES_INLINE_HEAD_INCLUDE``.

Example:

.. code-block::

    PROJECTROLES_INLINE_HEAD_INCLUDE="<meta name=\"keywords\" content=\"SODAR Core\">"

.. warning::

    Make sure you are inputting valid HTML or you risk breaking the HTML on
    **all** pages of your SODAR Core based site!


Modifying SODAR_CONSTANTS (Optional)
====================================

String identifiers used globally in SODAR project management are defined in the
``SODAR_CONSTANTS`` dictionary. It can be imported into your app code with the
import:

.. code-block:: python

    from projectroles.models import SODAR_CONSTANTS

If you need to update or extend the constants for use your site, you can import
the default dictionary into your Django settings and modify it as necessary with
the following snippet:

.. code-block:: python

    from projectroles.constants import get_sodar_constants
    SODAR_CONSTANTS = get_sodar_constants(default=True)
    # Your changes here..

.. warning::

    Modifying existing default constants may result in unwanted issues,
    especially on a site which already contains created projects. Proceed with
    caution!


Logging (Optional)
==================

It is recommended to add "projectroles" under ``LOGGING['loggers']``. For
production, ``ERROR`` debug level is recommended.

The example site and SODAR Django Site template provide the ``LOGGING_APPS`` and
``LOGGING_FILE_PATH`` helpers for easily adding SODAR Core apps to logging and
providing a system path for optional log file writing.
.. _app_projectroles_integration:


Projectroles Integration
^^^^^^^^^^^^^^^^^^^^^^^^

This document provides instructions and guidelines for integrating projectroles
and other SODAR Core apps into your Django site.


Installation on a New Site
**************************

If you want to set up a new Django site for integrating projectroles, see the
recommended options in this section.


SODAR Django Site Template (Recommended)
========================================

When setting up a new :term:`SODAR Core based site<SODAR Core Based Site>`, it
is strongly recommended to use
`sodar-django-site <https://github.com/bihealth/sodar-django-site>`_ as the
template. The repository contains a minimal :term:`Django site<Django Site>`
pre-configured with projectroles and other
:term:`SODAR Core apps<SODAR Core App>`. The ``main`` branch of this project
always integrates the latest stable release of SODAR Core and projectroles.

To set up your site with this template, clone the repository and follow the
installation instructions in the README.rst file.

To modify default SODAR Core and projectroles settings, see the
:ref:`app_projectroles_settings` document.

Once you have your site set up, you can look into
:ref:`customization tips <app_projectroles_custom>` and start
:ref:`developing your SODAR Core compatible apps <development>`.


Cookiecutter-Django
===================

If the SODAR Django site template does not suit your needs, it is also possible
to set up your site using `cookiecutter-django <https://github.com/pydanny/cookiecutter-django/>`_.
In this case, follow the instructions in the following section as if you were
integrating SODAR Core to an existing Django site.

.. note::

    The project was created using an old version of the cookiecutter script and
    evolved from there. This means the site created by the version currently may
    differ in several ways from how SODAR Core is set up. This method is
    recommended only for experienced Django developers.

.. note::

    For any other issues regarding the cookiecutter-django setup, see the
    cookiecutter-django documentation.


Installation on an Existing Site
********************************

Instructions for setting up projectroles and SODAR Core on an existing Django
site or a fresh site generated with cookiecutter-django are detailed in this
chapter.

.. warning::

    In order to successfully set up projectroles, you are expected to **follow
    all the instructions here in the order they are presented**. Please note
    that leaving out steps may result in a non-working Django site! Attempting
    to run the site before following all of the steps may (and probably will)
    result in errors.

.. warning::

    The rest of this section was originally written for the
    `1.11 release of cookiecutter-django <https://github.com/pydanny/cookiecutter-django/releases/tag/1.11.10>`_.
    Some details such as directory structures and settings variables may differ.

First, add the ``django-plugins`` and ``django-sodar-core`` package requirements
into your ``requirements/base.txt`` file. Make sure you are pointing to the
desired release tag or commit ID.

.. code-block:: console

    -e git+https://github.com/mikkonie/django-plugins.git@42e86e7904e5c09f1da32173862b26843eda5dd8#egg=django-plugins
    django-sodar-core==0.10.7

Install the requirements for development:

.. code-block:: console

    $ pip install -r requirements/local.txt

If any version conflicts arise between django-sodar-core and your existing site,
you will have to resolve them before continuing.

.. hint::

    You can always refer to either the ``sodar-django-site`` repository or
    ``example_site`` in the SODAR Core repository for a working example of a
    Cookiecutter-based Django site integrating SODAR Core. However, note that
    some aspects of the site configuration may vary depending on the
    cookiecutter-django version used on your site.


Django Settings
===============

Next you need to modify your default :term:`Django settings<Django Settings>`
file, usually located in ``config/settings/base.py``. For sites created with an
older cookiecutter-django version the file name may also be ``common.py``.
Naturally, you should make sure no settings in other configuration files
conflict with ones set here.

For values retrieved from environment variables, make sure to configure your
env accordingly. For development and testing, using ``READ_DOT_ENV_FILE`` is
recommended.

Required and optional Django settings are described in the
:ref:`app_projectroles_settings` document.


User Configuration
==================

In order for SODAR Core apps to work on your Django site, you need to extend the
default user model.

Extending the User Model
------------------------

In a cookiecutter-django based project, an extended user model should already
exist in ``{SITE_NAME}/users/models.py``. The abstract model provided by the
projectroles app provides the same model with critical additions, most notably
the ``sodar_uuid`` field used as an unique identifier for SODAR objects
including users.

If you have not added any of your own modifications to the model, you can simply
**replace** the existing model extension with the following code:

.. code-block:: python

    from projectroles.models import SODARUser

    class User(SODARUser):
        pass

If you need to include your own extra fields or functions (or have existing ones
already), you can add them in this model.

After updating the user model, create and run database migrations.

.. code-block:: console

    $ ./manage.py makemigrations
    $ ./manage.py migrate

.. note::

    You probably will need to edit the default unit tests under
    ``{SITE_NAME}/users/tests/`` for them to work after making these changes.
    See ``example_site.users.tests`` in this repository for an example.

Populating UUIDs for Existing Users
-----------------------------------

When integrating projectroles into an existing site with existing users, the
``sodar_uuid`` field needs to be populated. See
`instructions in Django documentation <https://docs.djangoproject.com/en/3.1/howto/writing-migrations/#migrations-that-add-unique-fields>`_
on how to create the required migrations.

Synchronizing User Groups for Existing Users
--------------------------------------------

To set up user groups for existing users, run the ``syncgroups`` management
command.

.. code-block:: console

    $ ./manage.py syncgroups

User Profile Site App
---------------------

The ``userprofile`` site app is installed with SODAR Core. It adds a user
profile page in the user dropdown. Use of the app is not mandatory but
recommended, unless you are already using some other user profile app. See
the :ref:`userprofile app documentation <app_userprofile>` for instructions.

Add Login Template
------------------

You should add a login template to ``{SITE_NAME}/templates/users/login.html``. If
you're OK with using the projectroles login template, the file can consist of
the following line:

.. code-block:: django

    {% extends 'projectroles/login.html' %}

If you intend to use projectroles templates for user management, you can delete
other existing files within the directory.


URL Configuration
=================

In the Django URL configuration file, usually found in ``config/urls.py``, add
the following lines under ``urlpatterns`` to include projectroles URLs in your
site.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'api/auth/', include('knox.urls')),
        url(r'^project/', include('projectroles.urls')),
    ]

If you intend to use projectroles views and templates as the basis of your site
layout and navigation (which is recommended), also make sure to set the site's
home view accordingly:

.. code-block:: python

    from projectroles.views import HomeView

    urlpatterns = [
        # ...
        url(r'^$', HomeView.as_view(), name='home'),
    ]

Finally, make sure your login and logout links are correctly linked. You can
remove any default allauth URLs if you're not using it.

.. code-block:: python

    from django.contrib.auth import views as auth_views

    urlpatterns = [
        # ...
        url(r'^login/$', auth_views.LoginView.as_view(
            template_name='users/login.html'), name='login'),
        url(r'^logout/$', auth_views.logout_then_login, name='logout'),
    ]


Base Template for Your Django Site
==================================

In order to make use of Projectroles views and templates, you should set the
base template of your site accordingly in ``{SITE_NAME}/templates/base.html``.

For a supported example, see ``projectroles/base_site.html``. It is strongly
recommended to use this as the base template for your site, either by extending
it or copying the content into ``{SITE_NAME}/templates/base.html`` and modifying
it to suit your needs.

If you do not need to make any modifications, the most simple way is to replace
the content of the ``{SITE_NAME}/templates/base.html`` file with the following
line:

.. code-block:: django

    {% extends 'projectroles/base_site.html' %}

.. note::

    CSS and Javascript includes in ``site_base.html`` are **mandatory** for
    Projectroles-based views and functionalities.

.. note::

    The container structure defined in the example base.html, along with
    including the ``{STATIC}/projectroles/css/projectroles.css`` are
    **mandatory** for Projectroles-based views to work without modifications.


Site Error Templates
====================

The projectroles app contains default error templates to use on your site.
These are located in the ``projectroles/error/`` template directory. You can
use them by entering ``{% extends 'projectroles/error/*.html %}`` in the
corresponding files found in the ``{SITE_NAME}/templates/`` directory. You have
the options of extending or replacing content on the templates, or simply
implementing your own.


Site Icons
==========

SODAR Core uses `Iconify <https://iconify.design/>`_ to include icons on the
site. SODAR Core apps use the `Material Design Icons <https://materialdesignicons.com>`_
icon collection, however other collections can be added to your site.

To enable the icons on your site, run the following management commands:

.. code-block:: console

    $ ./manage.py geticons
    $ ./manage.py collectstatic

If you want to use additional icon collections, you can add them using the
``-c`` argument as displayed in the following example:

.. code-block:: console

    $ ./manage.py geticons -c carbon clarity

You can view the supported icon collections
`here <https://github.com/iconify/collections-json/tree/master/json>`_.

The Iconify JSON files are rather large and potentially frequently updated, so
it is recommended to ignore them in your Git setup and instead retrieve them
dynamically for CI and deployment. Before committing your code, it is
recommended to update your ``.gitignore`` file with the following lines:

.. code-block::

    */static/iconify/*.json
    */static/iconify/json/*.json

All Done!
=========

After following all the instructions above, you should have a working SODAR Core
based Django site with support for projectroles features and SODAR Core apps. To
test the site locally execute the supplied make command:

.. code-block:: console

    $ make serve

Or, run the standard Django ``runserver`` command:

.. code-block:: console

    $ ./manage.py runserver

You can now browse your site locally at ``http://127.0.0.1:8000``. You are
expected to log in to view the site. Use e.g. the superuser account you created
when setting up your cookiecutter-django site.

You can now continue on to create apps or modify your existing apps to be
compatible with the SODAR Core framework. See the
:ref:`development section <development>` for app development guides. Also see the
:ref:`customization documentation <app_projectroles_custom>` for tips for
modifying the default appearance of SODAR Core.
.. _app_filesfolders_api_rest:


Filesfolders REST API Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This document contains the HTTP REST API documentation for the ``filesfolders``
app. The provided API enpoints allow file operations through HTTP API calls in
addition to the GUI.

For general information on REST API usage in SODAR Core, see
:ref:`app_projectroles_api_rest`.

.. currentmodule:: filesfolders.views_api

.. autoclass:: FolderListCreateAPIView

.. autoclass:: FolderRetrieveUpdateDestroyAPIView

.. autoclass:: FileListCreateAPIView

.. autoclass:: FileRetrieveUpdateDestroyAPIView

.. autoclass:: FileServeAPIView

.. autoclass:: HyperLinkListCreateAPIView

.. autoclass:: HyperLinkRetrieveUpdateDestroyAPIView.. _glossary:


Glossary
^^^^^^^^

.. glossary::

    App Plugin
        Mechanism for defining common properties and operations for dynamically
        including content and functionality from apps in SODAR Core views.

    App Settings
        Project or user specific settings defined in SODAR Core app plugins.
        Different from e.g. Django settings used to configure the web site.

    Backend App
        SODAR Core application which is used to provide additional functionality
        to other SODAR Core apps. Does not have its own GUI entry point. Common
        use cases include APIs to external services or other apps.

    Backend API
        Django API provided by a backend app, to be dynamically imported and
        used by other SODAR Core based Django apps.

    Django API
        Application programming interface offered by an app, to be used by other
        apps within the Django site.

    Django App
        Application built for the Django web framework, including (but not
        limited to) SODAR Core based apps.

    Django Settings
        Django settings used to configure the website. SODAR Core
        apps also use Django settings for configuring framework and app
        behaviour.

    Django Site
        Web site built on the Django framework, including (but not limited to)
        any website based on SODAR Core.

    Peer Site
        A SODAR Core based web site which mirrors one or more projects also
        mirrored on the currrently active site. This allows linking to remote
        projects on other sites where the user would have access.

    Project App
        SODAR Core application with the scope of providing data and
        functionality related to a specific project. Uses project-based access
        control.

    SODAR
        System for Omics Data Access and Retrieval. An omics research data
        management system which is the origin of the reusable SODAR Core
        framework.

    SODAR Core
        Core framework and reusable apps originally built for the SODAR project.

    SODAR Core App
        Django application with additional SODAR Core features. This includes
        one or more app plugin definitions to enable dynamic inclusion of the
        app into the SODAR Core framework, as well as project access control for
        project apps.

    SODAR Core Based Site
        Django-based web site using SODAR Core apps as its framework.

    Site App
        SODAR Core application with does not limit its scope to a single
        project. Common use cases include user account management and
        administrative tools.

    Source Site
        SODAR Core based web site which mirrors project metadata and access
        control to "target" sites.

    Target Site
        SODAR Core based web site which mirrors project metadata and access
        control from a "source" site.
.. _dev_sodar_core:


SODAR Core Development
^^^^^^^^^^^^^^^^^^^^^^

This document details instructions and guidelines for development of the SODAR
Core package.


Repository Contents
===================

In addition to the apps which will be installed by the package, the following
directories are included in the repository for development use and as examples:

- **config**: Example Django site configuration
- **docs**: Usage and development documentation
- **example_backend_app**: Example SODAR Core compatible backend app
- **example_project_app**: Example SODAR Core compatible project app
- **example_site**: Example SODAR Core based Django site for development
- **example_site_app**: Example SODAR Core compatible site-wide app
- **requirements**: Requirements for SODAR Core and development
- **utility**: Setup scripts for development


Installation
============

Instructions on how to install a local development version of SODAR Core are
detailed here. Ubuntu 20.04 LTS (Focal) is the supported OS at this time.
Later Ubuntu versions and Centos 7 have also been proven to to work, but some
system dependencies may vary for different OS versions or distributions.

Installation and development should be possible on most recent versions of
Linux, Mac and Windows, but this may require extra work and your mileage may
vary.

If you need to set up the accompanying example site in Docker, please see online
for up-to-date Docker setup tutorials for Django related to your operating
system of choice.

System Installation
-------------------

First you need to install OS dependencies, PostgreSQL >=9.6 and Python >=3.7.

.. code-block:: console

    $ sudo utility/install_os_dependencies.sh
    $ sudo utility/install_python.sh
    $ sudo utility/install_postgres.sh

Database Setup
--------------

Next you need to setup the database and postgres user. You'll be prompted to
enter a database name, a username and a password.

.. code-block:: console

    $ sudo utility/setup_database.sh

You have to set the database URL and credentials for Django in the environment
variable ``DATABASE_URL``. For development it is recommended to place
environment variables in file ``.env`` located in your project root. To enable
loading the file in Django, set ``DJANGO_READ_DOT_ENV_FILE=1`` in your
environment variables when running SODAR or any of its management commands.
See ``config/settings/base.py`` for more information and the ``env.example``
file for an example environment file.

Example of the database URL variable as set within an ``.env`` file:

.. code-block:: console

    DATABASE_URL=postgres://sodar_core:sodar_core@127.0.0.1/sodar_core

Project Setup
-------------

Clone the repository, setup and activate the virtual environment. Once in
the environment, install Python requirements for the project:

.. code-block:: console

    $ git clone https://github.com/bihealth/sodar-core.git
    $ cd sodar-core
    $ python3.x -m venv .venv
    $ source .venv/bin/activate
    $ utility/install_python_dependencies.sh

LDAP Setup (Optional)
---------------------

If you will be using LDAP/AD auth on your site, make sure to also run:

.. code-block:: console

    $ sudo utility/install_ldap_dependencies.sh
    $ pip install -r requirements/ldap.txt

Final Setup
-----------

Initialize the database (this will also synchronize django-plugins):

.. code-block:: console

    $ ./manage.py migrate

Create a Django superuser for the example_site:

.. code-block:: console

    $ ./manage.py createsuperuser

Now you should be able to run the server:

.. code-block:: console

    $ make serve


App Development
===============

Guidelines for developing **internal** SODAR Core apps (ones included when
installing the django-sodar-core package) are detailed in this section.

REST API Views
--------------

For internal SODAR Core apps, you need to use core counterparts to the mixins
than provided for SODAR Core using sites. The counterparts use different media
type and versioning from views to be implemented on external sites. This is to
prevent version number clashes requiring changes in external APIs. The classes
can be found in ``projectroles.views_api`` and are as follows:

- ``CoreAPIVersioning``
- ``CoreAPIRenderer``
- ``CoreAPIBaseMixin``
- ``CoreAPIBaseProjectMixin``
- ``CoreAPIGenericProjectMixin``

For detailed API descriptions, see docstrings in the ``view_api`` module. The
media type and versioning for these views are **hardcoded** and should not be
changed, except version information upon a new release of SODAR Core.


Projectroles App Development
============================

This section details issues regarding updates to the ``projectroles`` app.

.. warning::

    As all other apps in SODAR Core as well as sites implementing SODAR Core
    are based on projectroles, changes to this app need to be implemented and
    tested with extra care. Also make sure to provide detailed documentation for
    all breaking changes.

Projectroles App Settings
-------------------------

Projectroles defines its own app settings in ``projectroles/app_settings.py``.
These are not expected to be altered by SODAR Core based sites. These settings
add the ``local`` attribute, which allows/disallows editing the value on a
``TARGET`` site.

To alter projectroles app settings when developing the app, update the
``PROJECTROLES_APP_SETTINGS`` dictionary as follows:

.. code-block:: python

     'example_setting': {
         'scope': 'PROJECT',  # PROJECT/USER
         'type': 'STRING',  # STRING/INTEGER/BOOLEAN
         'default': 'example',
         'options': ['example', 'example2'],  # Optional, only for settings of type STRING or INTEGER
         'label': 'Project setting',  # Optional, defaults to name/key
         'placeholder': 'Enter example setting here',  # Optional
         'description': 'Example project setting',  # Optional
         'user_modifiable': True,  # Optional, show/hide in forms
         'local': False,  # Allow editing in target site forms if True
     }


Testing
=======

To run unit tests, you have to install the headless Chrome driver (if not yet
present on your system), followed by the Python test requirements:

.. code-block:: console

    $ sudo utility/install_chrome.sh
    $ pip install -r requirements/test.txt

Now you can run all tests with the following make command:

.. code-block:: console

    $ make test

If you want to only run a certain subset of tests, use e.g.:

.. code-block:: console

    $ make test arg=projectroles.tests.test_views

For running tests with SODAR Taskflow (not currently publicly available), you
can use the supplied make command:

.. code-block:: console

    $ make test_taskflow


Remote Site Development
=======================

For developing remote site features, you will want to run two or more SODAR Core
example sites concurrently: one ``SOURCE`` site and one or more ``TARGET``
sites.

For running a single ``TARGET`` site in addition to the main site, the fastest
way to get started is the following:

First, set up a second database called ``sodar_core_target`` using
``utility/setup_database.sh``.

Next, migrate the new database and create a superuser using
``make manage_target``. It is recommended to use a different admin user name
than on your ``SOURCE`` site, to help debugging.

.. code-block:: console

    $ make manage_target arg=migrate
    $ make manage_target arg=createsuperuser

Launch your site with ``make serve_target``. By default, you can access the site
at Port ``8001`` on localhost. The port can be altered by providing the
``target_port`` parameter, e.g. ``make serve_target target_port=8002``.
Management commands to the target site can be issued with the ``make manage_target``
make command.

Due to how cookies are set by Django, you currently may have to relogin when
switching to a different site on your browser. As a workaround you can launch
one of the sites in a private/incognito window or use different browsers.

If you need to create multiple target sites for testing ``PEER`` synchronization
features, make sure that you have a separate SODAR Core database for each site
and launch each site on a different port on localhost. The configuration
``local_target2.py`` is included for developing with multiple ``TARGET`` sites... _app_timeline_usage:


Timeline Usage
^^^^^^^^^^^^^^

Usage instructions for the ``timeline`` app are detailed in this document.


Timeline UI
===========

You can browse events by navigating to a project or a category and selecting the
"Timeline" app from the project sidebar. Alternatively, a view for site-wide
events (not linked to a project) can be found in the user menu on the top right
corner of the site.

By clicking on the time stamp for each event, you can see details of the event
execution (in case of e.g. asynchronous events).

By clicking on the clock icon next to an object link in the event description,
you can view the event history of that object. The link itself will take you
to the relevant view for the object on your Django site.

Certain events have a file icon in their description. If clicked, a popup
showing a collection of extra data for the given event will appear.
The popup will display extra data of the event itself and of the different
states the event went through, if there are any.

Admin users are able to see certain *"classified"* level events hidden from
normal users.

.. figure:: _static/app_timeline/sodar_timeline.png
    :align: center
    :scale: 50%

    Timeline event list view


Backend API for Event Logging
=============================

The Django backend API for logging events is located in ``timeline.api``. For
the full documentation, see `here <app_timeline_api_django>`_.

Invoking the API
----------------

The API is accessed through a backend plugin. This means you can write calls to
the API without any hard-coded imports and your code should work even if the
timeline app has not been installed on the site.

The most common use case is to save events within the Class-Based Views of your
Django site, but technically this can be done by any part of the code in your
Django apps.

Initialize the API using ``projectroles.plugins.get_backend_api()`` as follows:

.. code-block:: python

    from projectroles.plugins import get_backend_api
    timeline = get_backend_api('timeline_backend')

    if timeline:    # Only proceed if the backend was successfully initialized
        pass        # Save your events here..

Adding an Event
---------------

Once you can access the timeline backend, add the event with
``timeline.add_event()``. A minimal example is as follows:

.. code-block:: python

    tl_event = timeline.add_event(
        project=project,            # Project object (None for site-wide events)
        app_name=APP_NAME,          # Name of the app plugin for the current app
        user=request.user,          # The user triggering the activity (can be None)
        event_name='some_event',    # You can define these yourself, not unique
        description='Description')  # Human readable description

Linking an Object
-----------------

Say you want to link a Django model object to the event for tracking its
history? In this example, let's say it's a SODAR Core compatible ``User`` model
object ``user_obj``.

.. note::

    The given object **must** contain an ``sodar_uuid`` field with an
    auto-generated UUID. For more information, see the
    :ref:`project app development document <dev_project_app>`.

Create the event as in the previous section, but add a label ``target_user`` in
the description. The name of the label is arbitrary:

.. code-block:: python

    tl_event = timeline.add_event(
        project=project,
        app_name=APP_NAME,
        user=request.user,
        event_name='some_event',
        description='Do something to {target_user}')

All you have to do is add an object reference to the created event:

.. code-block:: python

    obj_ref = tl_event.add_object(
        obj=user_obj,
        label='target_user',
        name=user_obj.username)

The ``name`` field specifies which name the object will be referred to when
displaying the event description to a user.

Defining Object References
--------------------------

The example before is all fine and good for a ``User`` object, but what about
your own custom Django model?

When encountering an unknown object model from your app, timeline will call the
``get_object_link()`` function in the ``ProjectAppPlugin`` defined for your app.
Make sure to implement it for all the relevant models in your app.

Displaying Object Links
-----------------------

In order to display object links with timeline history link included, you can
use the ``timeline.api.get_object_link()`` function in your app's template tags.

Defining Status States
----------------------

.. note::

    If your Django apps only deal with normal synchronous requests, you don't
    need to pay attention to this functionality right now.

By default, ``timeline.add_event()`` treats events as synchronous and
automatically saves them with the status of ``OK``. However, in case of e.g.
asynchronous requests, you can alter this by setting the ``status_type`` and
(optionally) ``status_desc`` types upon creation.

.. code-block:: python

    tl_event = timeline.add_event(
        project=project,
        app_name=APP_NAME,
        user=request.user,
        event_name='some_event',
        description='Description',
        status_type='SUBMIT',
        status_desc='Just submitted this')

After that, you can add new status states for the event using the object
returned by ``timeline.add_event()``:

.. code-block:: python

    tl_event.set_status('OK', 'Submission was successful!')

Currently supported status types are listed below, some only applicable to async
events:

- ``OK``: All OK, event successfully performed
- ``INFO``: Used for events which do not change anything, e.g. viewing something
  within an app
- ``INIT``: Initializing the event in progress
- ``SUBMIT``: Event submitted asynchronously
- ``FAILED``: Asynchronous event submission failed
- ``CANCEL``: Event cancelled

Extra Data
----------

Extra data can be added in the JSON format for both events and their status
states with the ``extra_data`` and ``status_extra_data`` parameters.

Speciying a label ``{extra-NAME}`` in the event description will lead to a
callback to ``get_extra_data_link()`` in the app plugin. To support this you
need to make sure to implement the ``get_extra_data_link()`` function in your
plugin.

Classified Events
-----------------

To mark an event *"classified"*, that is, restricting its visibility to project
owners and admins, set the ``classified`` argument to true when invoking
``timeline.add_event()``.

.. note::

    Multiple levels of classification may be introduced to the timeline event
    model in the future.
.. _for_the_impatient:

For the Impatient
^^^^^^^^^^^^^^^^^

This section will give you the essential steps to setup a new SODAR Core based
project. We will link to the parts of the manual where they were taken from such
that you can read more in depth there.

.. _for_the_impatient_see_it_in_action:


See It In Action
================

We have developed the following data management and analysis web applications
using SODAR Core. Although there only is a public demo available for VarFish at
this time, the source code of the applications demonstrate how to use SODAR Core
in complex web applications.

  `VarFish <https://github.com/bihealth/varfish-server>`__
    is a web-based tool for the analysis of variants.
    It showcases how to build a complex data warehousing and data analysis web
    appliction using SODAR Core.
    More details are described in the `NAR Web Server Issue publication (doi:10.1093/nar/gkaa241) <https://doi.org/10.1093/nar/gkaa241>`__.
    The source code can be found on `github.com/bihealth/varfish-server <https://github.com/bihealth/varfish-server>`__.
    A demo is available at `varfish-demo.bihealth.org <https://varfish-demo.bihealth.org/login/>`__.

  `DigestiFlow <https://github.com/bihealth/digestiflow-server>`__
    is a web-based data system for the management and demultiplexing of Illumina
    Flow Cells. It further implements various tools for sanity checking Illumina
    sample sheets and quality control (e.g., comparing barcode adapter sequence
    and actual sequence present in the sequencer output).
    You can find out more in our publication in `Bioinformatics (doi:10.1093/bioinformatics/btz850) <https://doi.org/10.1093/bioinformatics/btz850>`__.
    The source code can be found on `github.com/bihealth/digestiflow-server <https://github.com/bihealth/digestiflow-server>`__.
    There currently is no public demo instance yet.

  `Kiosc <https://github.com/bihealth/kiosc>`__
    is a web application that allows to build scheduler Docker containers for
    "data science" apps and dashboards. There currently is no public demo
    instance yet.


Prerequisites
=============

Linux / Mac
    We use bash syntax on a Unix system and assume that you can adjust this to
    your system when needed.

PostgreSQL
    Please install version 9.6 or above.
    We assume that you have access to the ``postgres`` user or some other
    administrative user.

Development Essentials
    We assume that you have ``git``, Python 3.7 or above, and other essential
    tools installed.
    If you are using a mainstream Unix-like distribution (Mac qualifies) then
    you should be good to go.

.. info:

    In the case that you get an error as follows when installing the
    dependencies, make sure that you have the development libraries of postgres
    installed. E.g., on Debian-based systems install ``postgresql-dev``, for Red
    Hat and CentOS install ``postgresql-devel``.

    ::

        Error: pg_config executable not found.


Isolate Python Environment
==========================

If you use `virtualenv <https://docs.python-guide.org/dev/virtualenvs/>`_,
please create a new virtual environment for the project and activate it.
Otherwise, follow the previous link and do this now or you can follow along us
using `conda <https://www.google.com/search?client=ubuntu&channel=fs&q=conda&ie=utf-8&oe=utf-8>`_.

The following creates a new Miniconda installation on 64 bit Linux or Mac.
The `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ website has
URLs to more.

.. code-block:: bash

    # Linux
    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
    $ source ~/miniconda3/bin/activate
    $ conda install -y python=3.8

    # Mac
    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    $ bash Miniconda3-latest-MacOSX-x86_64.sh -b -p ~/miniconda3
    $ source ~/miniconda3/bin/activate
    $ conda install -y python=3.8

For activating the conda installation, use ``source ~/miniconda3/bin/activate``.


Install SODAR Core
==================

We simply use ``pip`` for this:

.. code-block:: bash

    # you must have your miniconda3 install sourced, skip if done already
    $ source ~/miniconda3/bin/activate
    $ pip install django-sodar-core==0.10.7


Download Example Site
=====================

We maintain a Git repository with a django project using the latest SODAR Core
version here on GitHub: `sodar-django-site <https://github.com/bihealth/sodar-django-site>`_.
We will use this to get hit the ground running.
See :ref:`app_projectroles_integration` on other ways to get started with SODAR
Core.

.. code-block:: bash

    $ git clone https://github.com/bihealth/sodar-django-site.git  sodar-django-site
    $ cd sodar-django-site

From here on, we assume that you are located (a) within the
``sodar-django-site`` directory and (b) have done
``source ~/miniconda3/bin/activate`` such that ``which python`` shows
``~/miniconda3/bin/python``.

To complete this step install the development requirements.

.. code-block:: bash

    $ pip install -r requirements/local.txt


Configure Environment
=====================

The next step is to perform some configuration.
SODAR Core is built on the `12 factor app <https://12factor.net/>`_ principles.
Configuration is done using environment variables.
For development, they are read from the ``.env`` file in your
``sodar-django-site`` checkout.
We are shipping an example setting file that you should copy and then edit.

.. code-block:: bash

    $ cp env.example .env
    # now edit .env

To start out, it will be sufficient to make sure you can connect to the database.
The default value for this is shown below.

.. code-block:: bash

    DATABASE_URL="postgres://sodar-django-site:sodar-django-site@127.0.0.1/sodar-django-site"

To keep it simple, you can use the following commands to create the correct
database, user, and set the password.

.. code-block:: bash

    $ sudo -u postgres createuser -ds sodar-django-site -W
    [sudo] password for USER: <enter your password>
    Password: <enter: sodar-django-site>
    $ sudo -u postgres createdb --owner=sodar-django-site sodar-django-site

Now, we have to make sure that the environment file is read:

.. code-block:: bash

    $ sed -ie "s/^READ_DOT_ENV_FILE.*/READ_DOT_ENV_FILE = env.bool('DJANGO_READ_DOT_ENV_FILE', default=True)/" config/settings/base.py


Database Initialization
=======================

For the final steps, you will initialize the database...

.. code-block:: bash

    $ python manage.py migrate

... and create a new admin user ...

.. code-block:: bash

    $ python manage.py createsuperuser
    Username: root
    Email address: root@example.com
    Password:
    Password (again):
    Superuser created successfully.


The First Login
===============

Now, start the server, and you can then visit http://127.0.0.1:8000/login/?next=/

.. code-block:: bash

    $ make serve
    python manage.py runserver --settings=config.settings.local
    Watching for file changes with StatReloader
    Performing system checks...

    System check identified no issues (0 silenced).
    July 01, 2021 - 13:06:06
    Django version 3.2.5, using settings 'config.settings.local'
    Starting development server at http://127.0.0.1:8000/
    Quit the server with CONTROL-C.

You should see the following:

.. figure:: _static/for_the_impatient/login_screen.png

Now, login with your superuser and you should see the following:

.. figure:: _static/for_the_impatient/project_list.png

Click the little user icon on the top right to access the django admin (here is
where you can create more users, for example) but also the preconfigured
:term:`Site Apps <Site App>` :ref:`adminalerts <app_adminalerts>`,
:ref:`siteinfo <app_siteinfo>`, :ref:`userprofile <app_userprofile>` and
configuration for remote sites.
The plus button on the left allows to create new categories and projects.

Now might also be a good time to read up more on the
:ref:`projectroles <app_projectroles>` app as this is the fundamental app for
most further development.


The First Project
=================

You cannot create projects on the root level but you have to create a new
category first (collections of projects). Use the "create category" button on
the left to create a "example category" first, then create an "example project"
within. The project details view should look as follows.

.. figure:: _static/for_the_impatient/project_details.png

Note that the site search already works, so typing "example" into the text field
on the top and clicking "search" will return your example project. The project
overview shows the "project home card" for the installed project apps
filefolders, timeline, and bgjobs. Usually, the latest five entries are shown
here.

.. note::

    The filesfolders app is an example of the **data management** application of
    SODAR Core based apps. You can easily imagine a more advanced module/app
    that not only allows tagging of files but more structuring data and meta
    data more strongly.

    SODAR Core was extracted as a re-useable library of component from
    our (not yet released) SODAR project that allows the management of
    structured meta data of experiments and files with data from such
    experiments.

Go ahead and try out the filesfolders app by clicking the "small files" icon on
the left. After creating folders and uploading a few files, you will see a trace
of actions in the timeline app:

.. figure:: _static/for_the_impatient/timeline.png

.. note::

    By default, ``sodar-django-site`` will store the files in the PostgreSQL
    database but you can easily configure it to use other storage backends,
    e.g., the S3 protocol, with the
    `django-storage <https://django-storages.readthedocs.io/en/latest/>`_
    package but this goes beyond this documentation.

Also, you will be able to find your uploaded file by name through the search
box. Note that any app that you write can easily provide all the integrations
with the SODAR Core framework (your apps are no different than the built-in
ones). Just have a look how we did it in the apps shipping with SODAR Core.


Summary
=======

Here is a quick summary on how SODAR Core interacts with the built-in and user
apps:

- At the lower most level all content is managed in projects which themselves
  can be assigned into categories.
- Project apps can provide new content types that can be put into projects.
  For example, the filesfolders app allows you to store files, folders, and
  assign meta data to them. As another example, the timelines app stores events
  that occured in a project, and other apps such as the filesfolders app can
  register their events with it.
- SODAR Core defines several plugin extension points that your apps can
  implement and make their content findable, for example.
- Site apps allow to provide features independent of a project.
  For example, the userprofile app allows to access user settings and the
  adminalerts app allows to post global notifications.


Going on From Here
==================

- You can now start exploring your ``sodar-django-site`` and play around with
  it.
- You can read the :ref`user_stories` section to learn how SODAR Core based
  applications are built.
- Continue reading :ref:`getting_started` for a more comprehensive documentation
  and walk-through of SODAR Core and its apps.
- Have a look at the web apps developed by us that are using SODAR Core as shown
  in the :ref:`for_the_impatient_see_it_in_action` section.
.. _app_projectroles_custom:

Projectroles Customization
^^^^^^^^^^^^^^^^^^^^^^^^^^

Here you can find some customization instructions and tips for projectroles and
SODAR Core.


CSS Overrides
=============

If some of the CSS definitions in ``{STATIC}/projectroles/css/projectroles.css``
do not suit your purposes, it is possible to override them in your own includes.
It is still recommended to include the *"Flexbox page setup"* section as
provided.

In this chapter are examples of overrides you can place e.g. in ``project.css``
to change certain defaults.

.. hint::

    While not explicitly mentioned, some parameters may require the
    ``!important`` argument to take effect on your site.

.. warning::

    In the future we may instead offer a full Bootstrap 4 theme, which may
    deprecate current overriding/extending CSS classes.

Static Element Coloring
-----------------------

If you wish to recolor the background of the static elements on the page
(title bar, side bar and project navigation breadcrumb), add the following
CSS overrides.

.. code-block:: css

    .sodar-base-navbar, .sodar-pr-sidebar, .sodar-pr-sidebar-nav {
      background-color: #ff00ff;
    }

    .sodar-pr-navbar {
      background-color: #00ff00;
    }

Sidebar Width
-------------

If the sidebar is not wide enough for your liking or e.g. a name of an app
overflowing, the sidebar can be resized with the following override:

.. code-block:: css

    .sodar-pr-sidebar {
        width: 120px;
    }


Title Bar
=========

You can implement your own title bar by replacing the default base.html include
of ``projectroles/_site_titlebar.html`` with your own HTML or include.

When doing this, it is possible to include elements from the default title bar
separately:

- Search form: ``projectroles/_site_titlebar_search.html``
- Site app and user operation dropdown:
  ``projectroles/_site_titlebar_dropdown.html``

See the templates themselves for further instructions.


Additional Title Bar Links
==========================

If you want to add additional links *not* related to apps in the title bar, you
can implement in the template file
``{SITE_NAME}/templates/include/_titlebar_nav.html``. This can be done for e.g.
documentation links or linking to external sites. Example:

.. code-block:: django

    {# Example extra link #}
    <li class="nav-item">
      <a href="#" class="nav-link" id="site-extra-link-x" target="_blank">
        <i class="iconify" data-icon="mdi:help-circle"></i> Extra Link
      </a>
    </li>


Site Logo
=========

An optional site logo can be placed into ``{STATIC}/images/logo_navbar.png`` to
be displayed in the default Projectroles title bar.


Custom Icon Collections
=======================

To use icons other than the default `Material Design Icons <https://materialdesignicons.com>`_
collection, download the corresponding collection JSON file from the
`Iconify JSON repository <https://github.com/iconify/collections-json/tree/master/json>`_
into `{SITE_NAME}/static/iconify/json`. After that run `collectstatic` and the
icons will be available using the collection identifier.


Project Breadcrumb
==================

To add custom content in the end of the default project breadcrumb, use
``{% block nav_sub_project_extend %}`` in your app template.

The entire breadcrumb element can be overridden by declaring
``{% block nav_sub_project %}`` block in your app template.


Footer
======

Footer content can be specified in the optional template file
``{SITE_NAME}/templates/include/_footer.html``.


Project and Category Display Names
==================================

If the *project* and *category* labels don't match your use case, it is possible
to change the labels displayed to the user by editing ``SODAR_CONSTANTS`` in
your Django site settings file. Example:

.. code-block:: python

    SODAR_CONSTANTS = get_sodar_constants(default=True)
    SODAR_CONSTANTS['DISPLAY_NAMES']['CATEGORY'] = {
        'default': 'not-a-category',
        'plural': 'non-categories',
    }
    SODAR_CONSTANTS['DISPLAY_NAMES']['PROJECT'] = {
        'default': 'not-a-project',
        'plural': 'non-projects',
    }

See more about overriding ``SODAR_CONSTANTS``
:ref:`here <app_projectroles_settings>`.

To print out these values in your views or templates, call the
``get_display_name()`` function, which is available both as a template tag
through ``projectroles_common_tags.py`` and a general utility function in
``utils.py``. Capitalization and pluralization are handled by the function
according to arguments.
See the :ref:`Django API documentation <app_projectroles_api_django>` for
details.

.. note::

    These changes will **not** affect role names or IDs and descriptions of
    Timeline events.

.. _app_taskflow:


Taskflow Backend
^^^^^^^^^^^^^^^^

The ``taskflowbackend`` backend app is an optional add-on used if your site
setup contains the separate **SODAR Taskflow** data transaction service.

If you have not set up a SODAR Taskflow service for any purpose, this backend
is not needed and can be ignored.


Basics
======

The ``taskflowbackend`` backend app is used in the main SODAR site to
communicate with an external SODAR Taskflow service to manage large-scale data
transactions. It has no views or database models and only provides an API for
other apps to use.

.. note::

    At the time of writing, SODAR Taskflow is in development and has not been
    made public.


Installation
============

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.

Django Settings
---------------

The taskflowbackend app is available for your Django site after installing
``django-sodar-core``. Add the app into ``THIRD_PARTY_APPS`` as follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'taskflowbackend.apps.TaskflowbackendConfig',
    ]

Next add the backend to the list of enabled backend plugins:

.. code-block:: python

    ENABLED_BACKEND_PLUGINS = env.list('ENABLED_BACKEND_PLUGINS', None, [
        # ...
        'taskflow',
    ])

The following app settings **must** be included in order to use the backend.
Note that the values for ``TASKFLOW_TARGETS`` depend on your SODAR Taskflow
configuration.

.. code-block:: python

    # Taskflow backend settings
    TASKFLOW_BACKEND_HOST = env.str('TASKFLOW_BACKEND_HOST', 'http://0.0.0.0')
    TASKFLOW_BACKEND_PORT = env.int('TASKFLOW_BACKEND_PORT', 5005)
    TASKFLOW_SODAR_SECRET = env.str('TASKFLOW_SODAR_SECRET', 'CHANGE ME!')
    TASKFLOW_TARGETS = [
        'sodar',
        # ..
    ]

Register Plugin
---------------

To register the taskflowbackend plugin, run the following management command:

.. code-block:: console

    $ ./manage.py syncplugins

You should see the following output:

.. code-block:: console

    Registering Plugin for taskflowbackend.plugins.BackendPlugin


Usage
=====

Once enabled, Retrieve the backend API class with the following in your Django
app python code:

.. code-block:: python

    from projectroles.plugins import get_backend_api
    taskflow = get_backend_api('taskflow')

See the docstrings of the API for more details.

To initiate sync of existing data with your SODAR Taskflow service, you can use
the following management command:

.. code-block:: console

    ./manage.py synctaskflow


Django API Documentation
========================

The ``TaskflowAPI`` class contains the SODAR Taskflow backend API. It should be
initialized using the ``Projectroles.plugins.get_backend_api()`` function.

.. autoclass:: taskflowbackend.api.TaskflowAPI
    :members:
.. _dev_site_app:


Site App Development
^^^^^^^^^^^^^^^^^^^^

This document details instructions and guidelines for developing **site apps**
to be used with the SODAR Core framework.

It is recommended to read :ref:`dev_project_app` before this document.


Site App Basics
===============

Site apps are basically normal Django apps *not* hooked to SODAR projects.
However, they provide a few nice features to be used in a SODAR-enabled Django
site:

- Rules for accessing app data (similar to project apps but without the need for
  being associated with a project)
- Dynamic inclusion into the site and default templates via plugins
- The ability to show site-wide messages to users


Prerequisites
=============

See :ref:`dev_project_app`.


Models
======

No specific model implementation is required. However, it is strongly to refer
to objects using ``sodar_uuid`` fields instead of the database private key.


Rules File
==========

Generate a ``rules.py`` file similar to a project app. However, you should not
use project predicates in this one. Example:

.. code-block:: python

    import rules
    # Allow viewing data
    rules.add_perm('{APP_NAME}.view_data', rules.is_authenticated)

If you allow anonymous users on your site and want to enable anonymous access to
your site app, use the ``is_allowed_anonymous`` predicate:

.. code-block:: python

    from projectroles import rules as pr_rules
    rules.add_perm(
        '{APP_NAME}.view_data',
        rules.is_authenticated | pr_rules.is_allowed_anonymous
    )


SiteAppPlugin
=============

Create a file ``plugins.py`` in your app's directory. In the file, declare a
``SiteAppPlugin`` class implementing
``projectroles.plugins.SiteAppPluginPoint``. Within the class, implement
member variables and functions as instructed in comments and docstrings.

.. code-block:: python

    from projectroles.plugins import SiteAppPluginPoint
    from .urls import urlpatterns

    class SiteAppPlugin(SiteAppPluginPoint):
        """Plugin for registering a site-wide app"""
        name = 'example_site_app'
        title = 'Example Site App'
        urls = urlpatterns
        # ...

The following variables and functions are **mandatory**:

- ``name``: App name (ideally should correspond to the app package name)
- ``title``: Printable app title
- ``urls``: Urlpatterns (usually imported from the app's ``urls.py`` file)
- ``icon``: Iconify collection and icon name (e.g. ``mdi:home``)
- ``entry_point_url_id``: View ID for the app entry point
- ``description``: Verbose description of app
- ``app_permission``: Basic permission for viewing app data in project (see
  above)

Implementing the following is **optional**:

- ``app_settings``: Implement if project or user specific settings for the app
  are needed. See the plugin point definition for an example.
- ``info_settings``: List of names for app-specific Django settings to be
  displayed for administrators in the siteinfo app.
- ``get_messages()``: Implement if your site app needs to display site-wide
  messages for users.
- ``get_statistics()``: Return statistics for the siteinfo app. See details in
  :ref:`the siteinfo documentation <app_siteinfo>`.
- ``get_object_link()``: Return object link for a Timeline event.
- ``get_extra_data_link()``: Return extra data link for a Timeline event.


Views
=====

In your views, you can still use projectroles mixins which are *not* related to
projects. Especially ``LoggedInPermissionMixin`` is useful to ensure users not
allowed to access a view are properly redirected. Example:

.. code-block:: python

    from django.views.generic import TemplateView
    from projectroles.views import LoggedInPermissionMixin

    class ExampleView(LoggedInPermissionMixin, TemplateView):
        """Site app example view"""
        permission_required = 'example_site_app.view_data'
        template_name = 'example_site_app/example.html'

.. note::

    The entry point URL is not expected to have any URL kwargs in the current
    implementation. If you intend to use a view which makes use of URL kwargs,
    you may need to modify it into also accepting a request without any
    parameters (e.g. displaying default content for the view).


Templates
=========

It is recommended for you to extend ``projectroles/base.html`` and put your
actual app content within the ``projectroles`` block. Example:

.. code-block:: django

    {# Projectroles dependency #}
    {% extends 'projectroles/base.html' %}
    {% load projectroles_common_tags %}

    {% block title %}
      Example Site App Page Title
    {% endblock title %}

    {% block projectroles %}

      <div class="container sodar-subtitle-container">
        <h2>
          <i class="iconify" data-icon="mdi:rocket-launch-outline"></i>
          Example Site App
        </h2>
      </div>

      <div class="container-fluid sodar-page-container">
        <div class="alert alert-info">
          This is an example and the entry point for <code>example_site_app</code>.
        </div>
      </div>

    {% endblock projectroles %}


Site App Messages
=================

The site app provides a way to display certain messages to users. For this, you
need to implement ``get_messages()`` in the ``SiteAppPlugin`` class.

If you need to control e.g. which user should see the message or removal of a
message after showing, you need to implement relevant logic in the function.

Example:

.. code-block:: python

    def get_messages(self, user=None):
        """
        Return a list of messages to be shown to users.
        :param user: User object (optional)
        :return: List of dicts or and empty list if no messages
        """
        return [{
            'content': 'Message content in here, can contain html',
            'color': 'info',        # Corresponds to bg-* in Bootstrap
            'dismissable': True     # False for non-dismissable
            'require_auth': True    # Only view for authorized users
        }]
.. _dev_backend_app:


Backend App Development
^^^^^^^^^^^^^^^^^^^^^^^

This document details instructions and guidelines for developing
**backend apps** to be used with the SODAR Core framework.

It is recommended to read :ref:`dev_project_app` before this document.


Backend App Basics
==================

Backend apps are intended as apps used by other apps via their plugin, without
requiring hard-coded imports. These may provide their own views for e.g. Ajax
API functionality, but mostly they're intended to be internal (hence the name).


Prerequisites
=============

See :ref:`dev_project_app`.


Models
======

No specific model implementation is required. However, it is strongly to refer
to objects using ``sodar_uuid`` fields instead of the database private key.


BackendAppPlugin
================

The plugin is detected and retrieved using a ``BackendAppPlugin``.

Declaring the Plugin
--------------------

Create a file ``plugins.py`` in your app's directory. In the file, declare a
``BackendAppPlugin`` class implementing
``projectroles.plugins.BackendPluginPoint``. Within the class, implement
member variables and functions as instructed in comments and docstrings.

.. code-block:: python

    from projectroles.plugins import BackendPluginPoint
    from .urls import urlpatterns

    class BackendAppPlugin(BackendPluginPoint):
        """Plugin for registering a backend app"""
        name = 'example_backend_app'
        title = 'Example Backend App'
        urls = urlpatterns
        # ...

The following variables and functions are **mandatory**:

- ``name``: App name (ideally should correspond to the app package name)
- ``title``: Printable app title
- ``icon``: Iconify collection and icon name (e.g. ``mdi:home``)
- ``description``: Verbose description of app
- ``get_api()``: Function for retrieving the API class for the backend, to be
  implemented

Implementing the following is **optional**:


- ``javascript_url``: Path to on demand includeable Javascript file
- ``css_url``: Path to on demand includeable CSS file
- ``info_settings``: List of names for app-specific Django settings to be
  displayed for administrators in the siteinfo app.
- ``get_statistics()``: Return statistics for the siteinfo app. See details in
  :ref:`the siteinfo documentation <app_siteinfo>`.
- ``get_object_link()``: Return object link for a Timeline event.
- ``get_extra_data_link()``: Return extra data link for a Timeline event.

.. hint::

    If you want to implement a backend API which is closely tied to a project
    app, there's no requirement to declare your backend as a separate Django
    app. You can just include the ``BackendAppPlugin`` in your app's
    ``plugins.py`` along with your ``ProjectAppPlugin``. See the
    :ref:`timeline app <app_timeline>` for an example of this.

Using the Plugin
----------------

To retrieve the API for the plugin, use the
function ``projectroles.plugins.get_backend_api()`` as follows:

.. code-block:: python

    from projectroles.plugins import get_backend_api
    example_api = get_backend_api('example_backend_app')

    if example_api:     # Make sure the API is there, and only after that..
        pass            # ..do stuff with the API

Including Backend Javascript/CSS
--------------------------------

If you want Javascript or CSS files to be associated with your plugin you can
set the ``javascript_url`` or ``css_url`` variables to specify the path to your
file. Note that these should correspond to ``STATIC`` paths under your app
directory.

.. code-block:: python

    class BackendPlugin(BackendPluginPoint):

        name = 'example_backend_app'
        title = 'Example Backend App'
        javascript_url = 'example_backend_app/js/example.js'
        css_url = 'example_backend_app/css/example.css'

The ``get_backend_include`` template-tag will return a ``<script>`` or
``<link>`` html tag with your specific file path, to be used in a template of
your app making use of the backend plugin:

.. code-block:: django

    {% load projectroles_common_tags %}
    {% get_backend_include 'example_backend_app' 'js' as javascript_include_tag %}
    {{ javascript_include_tag|safe }}

    {% get_backend_include 'example_backend_app' 'css' as css_include_tag %}
    {{ css_include_tag|safe }}

This will result in the following HTML:

.. code-block:: html

    <script type="text/javascript" src="/static/example.js"></script>
    <link rel="stylesheet" type="text/css" href="/static/example.css"/>

Be sure to use the backend plugin's name (not the title) as the key and filter
the result by ``safe``, so the tag won't be auto-escaped.
.. _app_projectroles_api_django:


Projectroles Django API Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This document contains the Django API documentation for the ``projectroles``
app. Included are functionalities and classes intended to be used by other
applications when building a SODAR Core based Django site.


Plugins
=======

SODAR plugin point definitions and helper functions for plugin retrieval are
detailed in this section.

.. automodule:: projectroles.plugins
    :members:


Models
======

Projectroles models are used by other apps for project access and metadata
management as well as linking objects to projects.

.. automodule:: projectroles.models
    :members:


App Settings
============

Projectroles provides an API for getting or setting project and user
specific settings.

.. autoclass:: projectroles.app_settings.AppSettingAPI
    :members:


Common Template Tags
====================

These tags can be included in templates with
``{% load 'projectroles_common_tags' %}``.

.. automodule:: projectroles.templatetags.projectroles_common_tags
    :members:


Utilities
=========

General utility functions are stored in ``utils.py``.

.. automodule:: projectroles.utils
    :members:


Base REST API View Classes
==========================

Base view classes and mixins for building REST APIs can be found in
``projectroles.views_api``.

.. currentmodule:: projectroles.views_api

Permissions / Versioning / Rendering
------------------------------------

.. autoclass:: SODARAPIProjectPermission
    :members:
    :show-inheritance:

.. autoclass:: SODARAPIVersioning
    :members:
    :show-inheritance:

.. autoclass:: SODARAPIRenderer
    :members:
    :show-inheritance:

Base API View Mixins
--------------------

.. autoclass:: SODARAPIBaseMixin
    :members:

.. autoclass:: SODARAPIBaseProjectMixin
    :members:
    :show-inheritance:

.. autoclass:: APIProjectContextMixin
    :members:
    :show-inheritance:

.. autoclass:: SODARAPIGenericProjectMixin
    :members:
    :show-inheritance:

.. autoclass:: ProjectQuerysetMixin
    :members:


Base Ajax API View Classes
==========================

Base view classes and mixins for building Ajax API views can be found in
``projectroles.views_ajax``.

.. currentmodule:: projectroles.views_ajax

.. autoclass:: SODARBaseAjaxView
    :members:
    :show-inheritance:

.. autoclass:: SODARBasePermissionAjaxView
    :members:
    :show-inheritance:

.. autoclass:: SODARBaseProjectAjaxView
    :members:
    :show-inheritance:


Base Serializers
================

Base serializers for SODAR Core compatible models are available in
``projectroles.serializers``.

.. currentmodule:: projectroles.serializers

.. autoclass:: SODARModelSerializer
    :members:
    :show-inheritance:

.. autoclass:: SODARProjectModelSerializer
    :members:
    :show-inheritance:

.. autoclass:: SODARNestedListSerializer
    :members:
    :show-inheritance:

.. autoclass:: SODARUserSerializer
    :members:
    :show-inheritance:
.. _app_projectroles_basics:


Projectroles Basics
^^^^^^^^^^^^^^^^^^^

The basic concepts and functionalities of the ``projectroles`` app are detailed
in this document.


Projects
========

The projectroles app groups data into **projects**. Here, a **project** is a
data container object that other objects can be linked to (typically through a
1:n foreign key relationship). A **category** is a sub-type of a project which
is allowed to contain other categories and projects but no other data type.

Using categories and projects, data can be organized in a tree structure of
category and project "containers". Users can be granted access to projects using
roles as described in the next section.


User Roles in Projects
======================

A **role** is a data type that has a string identifier in it score (e.g.,
"project guest"). Roles are assigned to individual users in the context of
individual projects (a n:m relation from user to project with the string
identifier). For example, user "paul" might be assigned the "project guest" role
in one project and another (or no) role in a second project. Users can only have
one role in a given project at any given time. New types of roles can be defined
by extending the default model's database table of the projectroles app.

The default setup of role types used in SODAR sites:

- **Project Owner**
    - Full read/write access to project data and roles
    - Can create sub-projects under owned categories
    - One per project
    - Must be specified upon project creation
- **Project Delegate**
    - Full read/write access to project data
    - Can modify roles except for owner and delegate
    - One per project (by default, the limit can be increased in site settings)
    - Assigned by owner
- **Project Contributor**
    - Can read and write project data
    - Can modify and delete own data
- **Project Guest**
    - Read only access to project data

.. note::

    Django **superuser** status overrides project role access.

The projectroles app provides the following features for managing user roles in
projects:

- Adding/modifying/removing site users as project members
- Inviting people not yet using the site by email
- Automated emailing of users regarding role changes
- Mirroring user roles to/from an external projectroles-enabled site

.. note::

    Currently, only superusers can assign owner roles for top-level categories.


Remote Project Sync
===================

SODAR Core allows optionally reading and synchronizing project metadata between
multiple SODAR-based Django sites. A superuser is able to set desired levels of
remote access for specific sites on a per-project basis.

A SODAR site can have one of three modes: **source**, **target** or **peer**
mode.

A SODAR site can be set by the user in either **source** or **target** mode.

- **Source site** is one expecting to (potentially) serve project metadata to
  an arbitrary number of other SODAR sites.
- **Target site** can be linked with exactly one source site, from which it
  can retrieve project metadata. Creation of local projects can be enabled or
  disabled according to local configuration.
- **Peer** mode is used only if two or more Target sites link to the same Source
  site. If synchronizing a project which has multiple accessing Target sites,
  metadata about those other Target sites is included and stored in Peer mode
  site objects.

Among the data which can be synchronized:

- General project information such as title, description and readme
- Project category structure
- User roles in projects
- User accounts for LDAP/AD users (required for the previous step)
- Information of other Target Sites linking a common project


Rule System
===========

Projectroles uses the `django-rules <https://github.com/dfunckt/django-rules>`_
package to manage permissions for accessing data, apps and functionalities
within projects based on the user role. Predicates for project roles are
provided by the projectroles app and can be used and extended for developing
rules for your other project-specific Django apps.


App Plugins
===========

Projectroles provides a plugin framework to enable integrating apps and
content dynamically to a projectroles-enabled Django site. Types of apps and
corresponding app plugins currently included:

- **Project apps**: Apps related to specific projects, making use of project
  access control and providing data and content within the project's scope
- **Site apps**: Site-wide Django apps which are not project-specific
- **Backend apps**: Backend apps without a GUI entry point, imported and used
  dynamically by other SODAR-based apps for e.g. connectivity to external
  resources.

App plugins are not limited to one per Django app. A single Django app in SODAR
Core may contain one or more of the aforementioned plugin types, depending on
the required functionality.

Existing apps can be modified to conform to the plugin structure by implementing
certain variables, functions, views and templates within the app. For more
details, see the app development documents.


Other Features
==============

Other features in the projectroles app:

- **App settings**: Setting values for project or user specific variables,
  which can be defined in project and site app plugins
- **Project starring**: Ability for users to star projects as their favourites
- **Project search**: Functionality for searching data within projects using
  functions implemented in project app plugins
- **Tour help**: Inline help for pages
- **Project readme**: README document for each project with Markdown support
- **Custom user model**: Additions to the standard Django user model
- **Multi-domain LDAP/AD support**: Support for LDAP/AD users from multiple
  domains
- **SODAR Taskflow and Timeline integration**: Included but disabled unless
  backend apps for Taskflow and Timeline are integrated in the Django site


Templates and Styles
====================

Projectoles provides views and templates for all GUI-related functionalities
described above. The templates utilize the plugin framework to provide content
under projects dynamically. The project also provides default CSS stylings, base
templates and a base layout which can be used or adapted as needed. See the
usage and app development documentation for more details.
.. _getting_started:


Getting Started
^^^^^^^^^^^^^^^

Installation and basic concepts of the SODAR Core framework and its apps are
detailed in this document.


Installation
============

The ``django-sodar-core`` package can be installed into your Django project
from PyPI as follows. It is strongly recommended to specify a version tag, as
the package is under active development and breaking changes are expected.

.. code-block:: console

    pip install django-sodar-core==0.10.7

Please note that the django-sodar-core package only installs
:term:`Django apps<Django App>`, which you need to include in a
:term:`Django web site<Django Site>` project. For instructions for integrating
SODAR Core into an existing Django site or setting up a new site,
see the :ref:`projectroles app documentation <app_projectroles>`.


SODAR Core Apps
===============

The following Django apps will be installed when installing the
``django-sodar-core`` package:

- **projectroles**: Base app for project access management and dynamic app
  content management. All other apps require the integration of projectroles.
- **adminalerts**: Site app for displaying site-wide messages to all users.
- **appalerts**: Site app and backend for raising alerts to users from apps.
- **bgjobs**: Project app for managing background jobs.
- **siteinfo**: Site app for displaying site information and statistics for
  administrators.
- **sodarcache**: Generic caching and aggregation of data referring to external
  services.
- **taskflowbackend**: Backend app providing an API for the optional
  ``sodar_taskflow`` transaction service.
- **timeline**: Project app for logging and viewing project-related activity.
- **tokens**: Token management for API access.
- **userprofile**: Site app for viewing user profiles.


Requirements
============

Major requirements for integrating projectroles and other SODAR Core apps into
your Django site are listed below. For a complete requirement list, see the
``requirements`` and ``utility`` directories in the repository.

- Ubuntu (20.04 Focal recommended and supported) / CentOS 7
- System library requirements (see the ``utility`` directory and/or your own
  Django project)
- Python >=3.7 (**NOTE:** Python 3.6 no longer supported in SODAR Core v0.10+)
- Django 3.2
- PostgreSQL >=9.6 and psycopg2-binary
- Bootstrap 4.x
- JQuery 3.3.x
- Shepherd and Tether
- Clipboard.js
- DataTables

For more details on installation and requirements for local development, see
:ref:`dev_sodar_core`.


Next Steps
==========

To proceed with using the SODAR Core framework in your Django site, you must
first install and integrate the ``projectroles`` app. See the
:ref:`projectroles app documentation <app_projectroles>` for instructions.

Once projectroles has been integrated into your site, you may proceed to
install other apps as needed.
.. _app_sodarcache_install:


Sodarcache Installation
^^^^^^^^^^^^^^^^^^^^^^^

This document provides instructions and guidelines for installing the
``sodarcache`` app to be used with your SODAR Core enabled Django site.

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.


Django Settings
===============

The sodarcache app is available for your Django site after installing
``django-sodar-core``. Add the app into ``THIRD_PARTY_APPS`` as follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'sodarcache.apps.SodarCacheConfig',
    ]

You also need to add the sodarcache backend plugin in enabled backend
plugins.

.. code-block:: python

    ENABLED_BACKEND_PLUGINS = [
        # ...
        'sodar_cache',
    ]


URL Configuration
=================

In the Django URL configuration file, add the following lines under
``urlpatterns`` to include sodarcache URLs in your site.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'^cache/', include('sodarcache.urls')),
    ]


Migrate Database and Register Plugin
====================================

To migrate the Django database and register the sodarcache app plugin, run the
following management command:

.. code-block:: console

    $ ./manage.py migrate

In addition to the database migration operation, you should see the following
output:

.. code-block:: console

    Registering Plugin for sodarcache.plugins.BackendPlugin
.. _app_tokens:


Tokens App
^^^^^^^^^^

The ``tokens`` site app enables users to issue and manage access tokens for REST
API views used on your SODAR Core based Django site.


Basics
======

Users can use this app to create and delete access tokens. These can be
set to expire or work until deleted.


Installation
============

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.

Django Settings
---------------

The siteinfo app is available for your Django site after installing
``django-sodar-core``. Add the app into ``THIRD_PARTY_APPS`` as follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'tokens.apps.TokensConfig',
    ]

URL Configuration
-----------------

In the Django URL configuration file, add the following line under
``urlpatterns`` to include siteinfo URLs in your site.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'^tokens/', include('tokens.urls')),
    ]

Register Plugin
---------------

To register the siteinfo site app plugin, run the following management command:

.. code-block:: console

    $ ./manage.py syncplugins

You should see the following output:

.. code-block:: console

    Registering Plugin for tokens.plugins.SiteAppPlugin


Usage
=====

When logged in to SODAR, you can find the "API Tokens" link in your user
dropdown menu in the top right corner of the site.

Select "Create Token" from the "Token Operations" dropdown to create a new
token. You will only see the token once, so make sure to copy it to clipboard at
this point.

Deleting existing tokens can be done from the token list.
.. _app_appalerts:


Appalerts App
^^^^^^^^^^^^^

The ``appalerts`` site app enables SODAR Core apps to generate alerts for users.
The user will receive a notify about these alerts anywhere on the site even if
they are tied to a specific project.


Basics
======

The app consists of two plugins: the site app plugin for alert displaying and
the backend plugin to create and access alerts from your SODAR Core based apps.


Installation
============

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.

Django Settings
---------------

Add the app into ``THIRD_PARTY_APPS`` in ``config/settings/base.py`` as follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'adminalerts.apps.AdminalertsConfig',
    ]

Next, make sure the context processor for app alerts is included in the
``TEMPLATES`` dictionary:

.. code-block:: python

    TEMPLATES = [
        {
            'OPTIONS': {
                'context_processors': {
                    # ...
                    'projectroles.context_processors.app_alerts_processor',
                }
            }
        }

To set the interval for JQuery status updates in the title bar, use
``APPALERTS_STATUS_INTERVAL`` and provide the value in seconds. Example:

.. code-block:: python

    APPALERTS_STATUS_INTERVAL = env.int('APPALERTS_STATUS_INTERVAL', 5)

URL Configuration
-----------------

In the Django URL configuration file, add the following line under
``urlpatterns`` to include adminalerts URLs in your site.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'^alerts/app/', include('appalerts.urls')),
    ]

Migrate Database and Register Plugin
------------------------------------

To migrate the Django database and register the appalerts site app plugin,
run the following management command:

.. code-block:: console

    $ ./manage.py migrate

In addition to the database migration operation, you should see the following
output:

.. code-block:: console

    Registering Plugin for appalerts.plugins.SiteAppPlugin
    Registering Plugin for appalerts.plugins.BackendPlugin

Base Template Include
---------------------

If your site is overriding the base site template, the following snippet should
be added in the ``javascript`` block of ``{SITE}/templates/base.html`` to enable
appalerts JQuery:

.. code-block:: django

    {% block javascript %}
      {# ... #}
      <!-- App alerts Javascript -->
      {% include 'projectroles/_appalerts_include.html' %}
    {% endblock javascript %}


Usage
=====

When logged in as an authenticated user, you can find the "App Alerts" option in
your user dropdown menu. This displays a list of active alerts you have in the
system. You have the possibility to dismiss alerts or follow a related link. The
latter action will also dismiss the alert.

When you are anywhere on the site, a notification about existing events will
appear on the top bar of the site.

.. figure:: _static/app_appalerts/app_alerts.png
    :align: center
    :scale: 50%

    App alert list and title bar notification


Backend API
===========

For creation and management of alerts, it is recommended to use the backend API
to retrieve and use the plugin, without the need for hard-coded includes. The
``add_alert()`` helper is also provided to simplify alert creation. See the
accompanying API documentation for details.

.. note::

    The logic for alert life cycle management is mostly left to the app issuing
    alerts. Alerts with an accompanying project or plugin will get deleted with
    the accompanying object. Please try to make sure you will not e.g. provide
    a related URL to the event which may no longer be valid when the user
    accesses the alert.


Backend Django API Documentation
================================

This is the backend API, retrievable with
``get_backend_api('appalerts_backend')``.

.. currentmodule:: appalerts.api

.. autoclass:: AppAlertAPI
    :members:
    :show-inheritance:
.. _app_userprofile:


Userprofile App
^^^^^^^^^^^^^^^

The ``userprofile`` app is a site app which provides a user profile view for
projectroles compatible Django users and management of user specific settings.


Installation
============

It is **strongly recommended** to install the userprofile app into your site
when using projectroles, unless you require a specific user profile providing
app of your own.

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.

Django Settings
---------------

The userprofile app is available for your Django site after installing
``django-sodar-core``. Add the app into ``THIRD_PARTY_APPS`` as follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'userprofile.apps.UserprofileConfig',
    ]

URL Configuration
-----------------

In the Django URL configuration file, add the following line under
``urlpatterns`` to include userprofile URLs in your site.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'^user/', include('userprofile.urls')),
    ]

Register Plugin
---------------

To register the app plugin, run the following management command:

.. code-block:: console

    $ ./manage.py syncplugins

You should see the following output:

.. code-block:: console

    Registering Plugin for userprofile.plugins.ProjectAppPlugin


Usage
=====

After successful installation, the link for "User Profile" should be available
in the user dropdown menu in the top-right corner of the website UI after you
have logged in.


User Settings
=============

User settings are configured in the ``app_settings`` dictionary in your project
app plugins.

User settings defined in the ``projectroles`` app, available for all SODAR Core
using sites:

Display Project UUID Copying Link
    If set true, display a link in the project title bar for copying the project
    UUID into the clipboard.
Additional Email
    In addition to the default user email, also send email notifications to
    these addresses.

In the development setup, the SODAR Core example site apps also provide
additional settings for demonstrating settings features.
.. _app_sodarcache_usage:


Sodar Cache Usage
^^^^^^^^^^^^^^^^^

Usage instructions for the ``sodarcache`` app are detailed in this document.


Backend API for Data Caching
============================

The Django backend API for caching data is located in ``sodarcache.api``. For
the full documentation, see `here <app_sodarcache_api_django>`_.

Invoking the API
----------------

The API is accessed through a backend plugin. This means you can write calls to
the API without any hard-coded imports and your code should work even if the
sodarcache app has not been installed on the site.

Initialize the API using ``projectroles.plugins.get_backend_api()`` as follows:

.. code-block:: python

    from projectroles.plugins import get_backend_api
    projectcache = get_backend_api('sodar_cache')

    if projectcache:    # Only proceed if the backend was successfully initialized
        pass

Setting and getting Cache Items
-------------------------------

Once you can access the sodarcache backend, you should set up the
``update_cache()`` function in the ``ProjectAppPlugin`` of the app with which
you want to cache or aggregate data. The update process can be limited by two
parameters: cached item name and project. If neither are specified, the function
should update cached data for all known items within all projects.

.. code-block:: python

        def update_cache(self, name=None, project=None):
        """
        Update cached data for this app, limitable to item ID and/or project.

        :param project: Project object to limit update to (optional)
        :param name: Item name to limit update to (string, optional)
        """
        # TODO: Implement this in your app plugin
        return None

Updating a specific cache item within the ``update_cache()`` function (or
elsewhere) should be done using ``sodarcache.api.set_cache_item()``. A minimal
example is as follows:

.. code-block:: python

    cache_item = projectcache.set_cache_item(
        project=project,            # Project object
        app_name=APP_NAME,          # Name of the current app
        user=request.user,          # The user triggering the cache update
        name='some_item',           # Cached item ID
        data_type='json',           # Data type ("json" currently supported)
        data={'key': 'val'},        # The actual data that should be cached
        )

.. note::

    The item ID in the ``name`` argument is not unique, but it is expected to
    be unique together with the ``project`` and ``app_name`` arguments.

Retrieve items with ``sodarcache.get_cache_item()`` or just check the
time the item was last updated with ``sodarcache.get_update_time()`` like
this:

.. code-block:: python

    projectcache.get_cache_item(
        app_name='yourapp',
        name='some_item',
        project=project,
        data_type='json'
    ) # Returns a JsonCacheItem

    projectcache.get_update_time(
        app_name='yourapp',
        name='some_item',
        project=project
    )

It is also possible to retrieve a Queryset with all cached items for a specific
project with ``sodarcache.get_project_cache()``

.. code-block:: python

    projectcache.get_project_cache(
        project=project,        # Project object
        data_type='json'        # must be 'json' for JsonCacheItem
        )

Using the Management commands
-----------------------------
To create or update the data cache for all apps and projects, you can use a
management command.

.. code-block:: console

    $ ./manage.py synccache

To limit the sync to a specific project, you can provide the ``-p`` or
``--project`` argument with the project UUID.

.. code-block:: console

    $ ./manage.py synccache -p e9701604-4ccc-426c-a67c-864c15aff6e2

Similarly, there is a command to delete all cached data:

.. code-block:: console

    $ ./manage.py deletecache



.. _app_sodarcache_api_django:


Sodarcache Backend API Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This document contains Django API documentation for the backend plugin in the
``sodarcache`` app. Included are functionalities and classes intended to be used
by other applications.


Backend API
===========

The ``SodarCacheAPI`` class contains the Sodar Cache backend API. It should be
initialized with ``Projectroles.plugins.get_backend_api('sodar_cache')``.

.. autoclass:: sodarcache.api.SodarCacheAPI
    :members:


Models
======

.. automodule:: sodarcache.models
    :members:
.. _overview:


SODAR Core Overview and Example Use Case
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This document presents an overview of the SODAR Core package along with an
example use case.

SODAR Core is a relatively complex system and we have created a :ref:`glossary`
to help you with keeping track of the terminology.


SODAR Core Overview
===================

The SODAR Core package provides a suite of :term:`Django apps<Django App>` to be
installed on a :term:`Django-based web site<Django Site>`. The main app in the
package, ``projectroles``, provides core project access, content management
framework and default UI templates for other apps on the site. Those apps must
implement or use specific parts of the projectroles app to enable desired SODAR
Core functionality.

Apps in a SODAR Core based site are separated into :term:`project<Project App>`,
:term:`site<Site App>` and :term:`backend<Backend App>` apps, depending on their
scope and purpose. The SODAR Core package includes optional general purpose apps
of each type, which the user may enable on their site if needed. These apps all
depend on projectroles. More on the general purpose apps can be found in the
:ref:`getting_started` document.

To build their own web based system with SODAR Core, the user will develop
required functionality and UIs as one or more Django apps, using and extending
functionalities offered by the projectroles app and optional backend apps. This
allows integration of the app into the project access management, standardized
layout and other features such as advanced logging. Furthermore, the
projectroles app will call certain functions implemented in the user's apps to
dynamically include app and project content in Django views. In addition to
developing new Django apps, existing apps can be easily modified to gain access
to SODAR Core features.

.. figure:: _static/overview/sodar_core_site.png
    :align: center
    :scale: 50%

    Structure of a SODAR Core based web site


Example Use Case
================

In a typical scenario for SODAR Core use, a research organization wants to
develop a user friendly system for accessing, browsing and/or manipulating
research data. The data may belong to several different projects, with different
research groups or scientists working on it. This data may be confidential in
nature, so access control is required, preferably using the organization's
existing LDAP/AD servers. The organization wants to provide a web-based GUI as
well as programmatic API views.

Site Setup
----------

First the developer will
:ref:`set up a SODAR Core based Django site<app_projectroles_integration>`.
This can be created from a SODAR Core site template or by integrating the
``django-sodar-core`` package on an existing site. The Django site must be
configured according to organization requirements, e.g. setting up user access
via the organization LDAP/AD server.

The initial SODAR Core integration into a Django site adds e.g. a general
layout, organization of apps and content into projects as well as project-based
access control.

App Development
---------------

Next the developer needs to :ref:`develop the Django app(s)<development>` which
contain the actual program logic and user interfaces required for the use case.
The content or number of these apps are not restricted by SODAR Core. Anything
which would go to a typical Django app is OK, as long as certain building blocks
for SODAR Core functionality are used.

Optional applications bundled with SODAR Core can also be enabled or disabled at
this point. for example, if involved projects require sharing and storage of
reports and presentations, the developer may select to enable the
:ref:`filesfolders<app_filesfolders>` app. Likewise, if detailed logging or
audit trails are needed, the developer can enable the
:ref:`timeline<app_timeline>` app.

User and Project Setup
----------------------

Once the site is deployed, the developer should
:ref:`create initial project categories<app_projectroles_usage>` and provide
access to those for high level personnel such as project owners. The owners can
then go on and create relevant projects, grant access to researchers and set up
initial data in the applications.

Using the Site
--------------

The researchers will log in to the site on their web browser, in most cases
using the standard LDAP credentials provided by their organization. They will
see the projects they have been granted access to and can use whichever
applications have been enabled or developed for the site, according to their
assigned user rights. SODAR Core provides common navigation, overview and search
views for all enabled apps, including the one(s) developed by the organization.
The same user access management features are shared for all apps, along with
possible REST APIs developed by the organization.


Next Steps
==========

See the :ref:`getting_started` document for instructions on installing SODAR
Core and how to proceed with developing your own SODAR Core based site.
.. _manual_main:

Welcome to the SODAR Core documentation!
========================================

SODAR Core is a framework for Django web application development.

It was conceived to facilitate the creation of scientific data management and
analysis web applications, but can be useful in other contexts as well.
In that it is similar to the CMS or ecommerce frameworks that you can find
`Awesome Django List <https://github.com/wsvincent/awesome-django#content-management-systems>`__ but you will find the components/libraries provided in SODAR Core are more generic and in this reflecting the broader range of applications that we target.

How to Read This Manual
-----------------------

There are two ways:

Front to Back
  If you have the time and patience, reading the whole manual will teach
  you everything.

Jump Around (recommended)
  Start with :ref:`for_the_impatient` and/or :ref:`user_stories`, skim over the
  summary of each app, and explore what interests you most.

What SODAR Core Is and What It is Not
-------------------------------------

SODAR Core

- is Django-based and you will need some knowledge about Django programming in
  Python for it to be useful,
- provides you with libraries for developing your own applications.

SODAR Core

- is NOT a ready-made web application,
- is NOT for entry-level Python programmers (you will need intermediate Django
  knowledge; you probably do not want to base your first Django web application
  on SODAR Core).

What's Inside SODAR Core
------------------------

The full list of apps are shown in the table of contents (on the left if you
are reading the HTML version of this documentation). Here are some highlights:

- Project-based user access control
- Dynamic app content management
- Advanced project activity logging
- Small file uploading and browsing
- Managing server-side background jobs
- Caching and aggregation of data from external services
- Tracking site information and statistics

What's Inside This Documentation
--------------------------------

Overview & Getting Started
  This part aims at getting you an birds-eye view of SODAR Core and its usage.

SODAR Core Apps
  This part documents each Django app that ships with SODAR. As
  a reminder, in Django development, *apps* are re-useable modules with code
  for supporting a certain use case.

Project Info
  This part of the documentation provides meta information about the project
  and the full changelog.

What's Not Inside This Documentation
------------------------------------

You should know the following before this documentation is useful to you:

Python Programming
  There's tons of documentation on the internet but the `official Python
  documentation <https://docs.python.org/3/>`_ is a good starting point as
  any.

Django Development
  For learning about Django, head over to the `excellent documentation of the
  Django Project <https://docs.djangoproject.com/en/3.2/>`_.

HTML / Javascript / CSS / Bootstrap 4
  Together with Django, SODAR Core provides a framework to plug in your own
  HTML and related front-end code. We assume that you have web development
  experience and in particular know your way around Bootstrap 4.

  We're using the Bootstrap 4 CSS framework and you can learn about it in many
  places including `the official documentation
  <https://getbootstrap.com/docs/4.3/getting-started/introduction/>`_

.. note::

   You can find the official version of this documentation at
   `readthedocs.io <https://sodar-core.readthedocs.io/en/latest/>`_.
   If you view these files on GitHub, beware that their renderer does not
   render the ReStructuredText files correctly and content may be missing.

.. toctree::
    :caption: Overview & Getting started
    :name: overview_getting_started
    :hidden:
    :titlesonly:

    Overview <overview>
    getting_started
    for_the_impatient
    user_stories

.. toctree::
    :maxdepth: 2
    :caption: SODAR Core Apps
    :name: sodar_core_apps
    :hidden:
    :titlesonly:

    app_projectroles
    app_adminalerts
    app_appalerts
    app_bgjobs
    app_filesfolders
    app_userprofile
    app_siteinfo
    app_sodarcache
    app_taskflow
    app_timeline
    app_tokens

.. toctree::
    :maxdepth: 2
    :caption: Project Info
    :name: project_info
    :hidden:
    :titlesonly:

    contributing
    code_of_conduct
    glossary
    development
    major_changes
    Full Changelog <changelog>

Indices and tables
==================

* :ref:`modindex`
* :ref:`genindex`
* :ref:`search`
.. _dev_project_app:


Project App Development
^^^^^^^^^^^^^^^^^^^^^^^

This document details instructions and guidelines for developing
**project apps** to be used with the SODAR Core framework. This also applies for
modifying existing Django apps into project apps.

.. hint::

   The package ``example_project_app`` in the projectroles repository provides
   a concrete minimal example of a working project app.


Project App Basics
==================

**Characteristics** of a project app:

- Provides a functionality related to a project
- Is dynamically included in project views by projectroles using plugins
- Uses the project-based role and access control provided by projectroles
- Is included in projectroles search (optionally)
- Provides a dynamically included element (e.g. content overview) for the
  project details page
- Appears in the project menu sidebar in the default projectroles templates

**Requirements** for setting up a project app:

- Implement project relations and SODAR UUIDs in the app's Django models
- Use provided mixins, keyword arguments and conventions in views
- Extend projectroles base templates in your templates
- Implement specific templates for dynamic inclusion by Projectroles
- Implement ``plugins.py`` with definitions and function implementations
- Implement ``rules.py`` with access rules

Fulfilling these requirements is detailed further in this document.


Prerequisites
=============

This documentation assumes you have a Django project with the ``projectroles``
app set up (see the
:ref:`projectroles integration document <app_projectroles_integration>`).
The instructions can be applied either to modify a previously existing app, or
to set up a fresh app generated in the standard way with
``./manage.py startapp``.

It is also assumed that apps are more or less created according to best
practices defined by `Two Scoops <https://www.twoscoopspress.com/>`_, with the
use of `Class-Based Views <https://docs.djangoproject.com/en/3.2/topics/class-based-views//>`_
being a requirement.


Models
======

In order to hook up your Django models into projects, there are two
requirements: implementing a **project foreign key** and a **UUID field**.

Project Foreign Key
-------------------

Add a ``ForeignKey`` field for the ``projectroles.models.Project`` model,
either called ``project`` or accessible with a ``get_project()`` function
implemented in your model.

If the project foreign key for your is **not** ``project``, make sure to define
a ``get_project_filter_key()`` function. It should return the name of the field
to use as key for filtering your model by project.

.. note::

    If your app contains a complex model structure with e.g. nested models using
    foreign keys, it's not necessary to add this to all your models, just the
    topmost one(s) used e.g. in URL kwargs.

Model UUID Field
----------------

To provide a unique identifier for objects in the SODAR context, add a
``UUIDField`` with the name of ``sodar_uuid`` into your model.

.. note::

    SODAR Core links to objects in URLs, links and forms using UUIDs instead
    of database private keys. This is strongly recommended for all Django models
    in apps using the SODAR Core framework.

.. note::

    When updating an existing Django model with an existing database, the
    ``sodar_uuid`` field needs to be populated. See
    `instructions in Django documentation <https://docs.djangoproject.com/en/3.2/howto/writing-migrations/#migrations-that-add-unique-fields>`_
    on how to create the required migrations.

Model Example
-------------

Below is an example of a projectroles-compatible Django model:

.. code-block:: python

    import uuid
    from django.db import models
    from projectroles.models import Project

    class SomeModel(models.Model):
        some_field = models.CharField(
            help_text='Your own field'
        )
        project = models.ForeignKey(
            Project,
            related_name='some_objects',
            help_text='Project to which this object belongs',
        )
        sodar_uuid = models.UUIDField(
            default=uuid.uuid4,
            unique=True,
            help_text='SomeModel SODAR UUID',
        )

.. note::

    The ``related_name`` field is optional, but recommended as it provides an
    easy way to lookup objects of a certain type related to a project. For
    example the ``project`` foreign key in a model called ``Document`` could
    feature e.g. ``related_name='documents'``.


Rules File
==========

Create a file ``rules.py`` in your app's directory. You should declare at least
one basic rule for enabling a user to view the app data for the project. This
can be named e.g. ``{APP_NAME}.view_data``. Predicates for the rules can be
found in projectroles and they can be extended within your app if needed.

.. code-block:: python

    import rules
    from projectroles import rules as pr_rules

    rules.add_perm(
        'example_project_app.view_data',
        pr_rules.is_project_owner
        | pr_rules.is_project_delegate
        | pr_rules.is_project_contributor
        | pr_rules.is_project_guest,
    )

.. hint::

    The ``rules.is_superuser`` predicate is often redundant, as permission
    checks are skipped for Django superusers. However, it can be handy if you
    e.g. want to define a rule allowing only superuser access for now, with the
    potential for adding other predicates later.


ProjectAppPlugin
================

Create a file ``plugins.py`` in your app's directory. In the file, declare a
``ProjectAppPlugin`` class implementing
``projectroles.plugins.ProjectAppPluginPoint``. Within the class, implement
member variables and functions as instructed in comments and docstrings.

.. code-block:: python

    from projectroles.plugins import ProjectAppPluginPoint
    from .urls import urlpatterns

    class ProjectAppPlugin(ProjectAppPluginPoint):
        """Plugin for registering app with Projectroles"""
        name = 'yourprojectapp'
        title = 'Your Project App'
        urls = urlpatterns
        # ...

The following variables and functions are **mandatory**:

- ``name``: App name (**NOTE:** should correspond to the app package name or
  some functionality may not work as expected)
- ``title``: Printable app title
- ``urls``: Urlpatterns (usually imported from the app's ``urls.py`` file)
- ``icon``: Iconify collection and icon name (e.g. ``mdi:home``)
- ``entry_point_url_id``: View ID for the app entry point (**NOTE:** The view
  **must** take the project ``sodar_uuid`` as a kwarg named ``project``)
- ``description``: Verbose description of app
- ``app_permission``: Basic permission for viewing app data in project (see
  above)
- ``search_enable``: Boolean for enabling/disabling app search
- ``details_template``: Path to template to be included in the project details
  page, usually called ``{APP_NAME}/_details_card.html``
- ``details_title``: Title string to be displayed in the project details page
  for the app details template
- ``plugin_ordering``: Number to define the ordering of the app on the project
  menu sidebar and the details page

Implementing the following is **optional**:

- ``app_settings``: Implement if project, user or project_user (Settings
  specific to a project and user) specific settings for the app are needed. See
  the plugin point definition for an example.
- ``search_types``: Implement if searching the data of the app is enabled
- ``search_template``: Implement if searching the data of the app is enabled
- ``project_list_columns``: Optional custom columns do be shown in the project
  list. See the plugin point definition for an example.
- ``category_enable``: Whether the app should also be made available for
  categories. Defaults to ``False`` and should only be overridden when required.
  For an example of a project app enabled in categories, see
  :ref:`Timeline <app_timeline>`.
- ``info_settings``: List of names for app-specific Django settings to be
  displayed for administrators in the siteinfo app.
- ``get_taskflow_sync_data()``: Applicable only if working with
  ``sodar_taskflow`` and iRODS.
- ``get_object_link()``: Return object link for a Timeline event.
- ``get_extra_data_link()``: Return extra data link for a Timeline event.
- ``search()``: Function called when searching for data related to the app if
  search is enabled.
- ``get_statistics()``: Return statistics for the siteinfo app. See details in
  :ref:`the siteinfo documentation <app_siteinfo>`.
- ``get_project_list_value()``: A function which **must** be implemented if
  ``project_list_columns`` are defined, to retrieve a column cell value for a
  specific project.
- ``handle_project_update()``: A function for enabling carrying out specific
  tasks within your app when the project is updated in projectroles. This is a
  work-in-progress functionality to be expanded later.

Once you have implemented the ``rules.py`` and ``plugins.py`` files and added
the app and its URL patterns to the Django site configuration, you can create
the project app plugin in the Django databse with the following command:

.. code-block:: console

    $ ./manage.py syncplugins

You should see the following output to ensure the plugin was successfully
registered:

.. code-block:: console

    Registering Plugin for {APP_NAME}.plugins.ProjectAppPlugin

For info on how to implement the specific required views/templates, see the end
of this document.


Views
=====

Certain guidelines must be followed in developing Django web UI views for them
to be successfully used with projectroles.

URL Keyword Arguments
---------------------

In order to link a view to project and check user permissions using mixins,
the URL keyword arguments **must** include an argument which matches *one of
the following conditions*:

- Contains a kwarg ``project`` which corresponds to the ``sodar_uuid``
  member value of a ``projectroles.models.Project`` object
- Contains a kwarg corresponding to the ``sodar_uuid`` of another Django
  model, which must contain a member field ``project`` which is a foreign key
  for a ``Projectroles.models.Project`` object. The kwarg **must** be named
  after the Django model of the referred object (in lowercase).
- Same as above, but the Django model provides a
  ``get_project()`` function which returns (you guessed it) a
  ``Projectroles.models.Project`` object.

Examples:

.. code-block:: python

   urlpatterns = [
       # Direct reference to the Project model
       url(
           regex=r'^(?P<project>[0-9a-f-]+)$',
           view=views.ProjectDetailView.as_view(),
           name='detail',
       ),
       # RoleAssignment model has a "project" member which is also OK
       url(
           regex=r'^members/update/(?P<roleassignment>[0-9a-f-]+)$',
           view=views.RoleAssignmentUpdateView.as_view(),
           name='role_update',
       ),
   ]

Mixins
------

The ``projectroles.views`` module provides several useful mixins for augmenting
your view classes to add projectroles functionality. These can be found in the
``projectroles.views`` module.

The most commonly used mixins:

- ``LoginRequiredMixin``: Override of the standard Django mixin
  which may also allow anonymous guests if so configured in SODAR Core. If you
  plan on supporting anonymous users on your site, you **must** use this
  mixing instead of the original one in Django.
- ``LoggedInPermissionMixin``: Ensure correct redirection of users on no
  permissions.
- ``ProjectPermissionMixin``: Provides a ``Project`` object for permission
  checking based on URL kwargs.
- ``ProjectContextMixin``: Provides a ``Project`` object into the view context
  based on URL kwargs.

See ``example_project_app.views.ExampleView`` for an example.


Templates
=========

Template Structure
------------------

It is strongly recommended to extend ``projectroles/project_base.html`` in your
project app templates. Just start your template with the following line:

.. code-block:: django

    {% extends 'projectroles/project_base.html' %}

The following **template blocks** are available for overriding or extending when
applicable:

- ``title``: Page title
- ``css``: Custom CSS (extend with ``{{ block.super }}``)
- ``projectroles_extend``: Your app content goes here!
- ``javascript``: Custom Javascript (extend with ``{{ block.super }}``)
- ``head_extend``: Optional block if you need to include additional content
  inside the HTML ``<head>`` element

Within the ``projectroles_extend`` block, it is recommended to use the
following ``div`` classes, both extending the Bootstrap 4 ``container-fluid``
class:

- ``sodar-subtitle-container``: Container for the page title
- ``sodar-content-container``: Container for the actual content of your app

If you do not want to include the project title header to your project
templates, you can replace the ``projectroles_extend`` block with a
``projectroles`` block.

.. warning::

    When customizing your templates, make sure you are not accidentally nesting
    built-in blocks within eachother by e.g. placing the ``css`` block *inside*
    the ``projectroles`` or ``projectroles_extend`` block. Doing so may cause
    the page to render incorrectly or includes to fail.

Rules
-----

To control user access within a template, just do it as follows:

.. code-block:: django

    {% load rules %}
    {% has_perm 'app.do_something' request.user project as can_do_something %}

This checks if the current user from the HTTP request has permission for
``app.do_something`` in the current project retrieved from the page context.

Template Tags
-------------

General purpose template tags are available in
``projectroles/templatetags/projectroles_common_tags.py``. Include them to your
template as follows:

.. code-block:: django

    {% load projectroles_common_tags %}

Example
-------

Minimal example for a project app template:

.. code-block:: django

    {% extends 'projectroles/project_base.html' %}

    {% load projectroles_common_tags %}
    {% load rules %}

    {% block title %}
      Page Title
    {% endblock title %}

    {% block head_extend %}
      {# OPTIONAL: extra content under <head> goes here #}
    {% endblock head_extend %}

    {% block css %}
      {{ block.super }}
      {# OPTIONAL: Extend or override CSS here #}
    {% endblock css %}

    {% block projectroles_extend %}

      {# Page subtitle #}
      <div class="container-fluid sodar-subtitle-container">
        <h3>
          <i class="iconify" data-icon="mdi:rocket-launch"></i>
          App and/or Page Title
        </h3>
      </div>

      {# App content #}
      <div class="container-fluid sodar-page-container">
        <p>Your app content goes here!</p>
      </div>

    {% endblock projectroles_extend %}

    {% block javascript %}
      {{ block.super }}
      {# OPTIONAL: include additional Javascript here #}
    {% endblock javascript %}

See ``example_project_app/example.html`` for a working and fully commented
example of a minimal template.

.. hint::

    If you include some controls on your ``sodar-subtitle-container`` class and
    want it to remain sticky on top of the page while scrolling, use ``row``
    instead of ``container-fluid`` and add the ``bg-white sticky-top`` classes
    to the element.


General Guidelines for Views and Templates
==========================================

General guidelines and hints for developing views and templates are discussed
in this section.

Referring to Project Type
-------------------------

As of SODAR Core v0.4.3, it is possible to customize the display name for the
project type from the default "project" or "category". For more information, see
:ref:`app_projectroles_custom`.

It is thus recommended that instead of hard coding "project" or "category" in
your views or templates, use the ``get_display_name()`` function to refer to
project type.

In templates, this can be achieved with a custom template tag. Example:

.. code-block:: django

    {% load projectroles_common_tags %}
    {% get_display_name project.type title=True plural=False %}

In views and other Python code, the similar function can be accessed through
``utils.py``:

.. code-block:: python

    from projectroles.utils import get_display_name
    display_name = get_display_name(project.type, plural=False)

.. hint::

    If not dealing with a ``Project`` object, you can provide the
    ``PROJECT_TYPE_*`` constant from ``SODAR_CONSTANTS``. In templates, it's
    most straightforward to use "CATEGORY" and "PROJECT".


Forms
=====

This section contains guidelines for implementing forms.

SODAR User Selection Field
--------------------------

Projectroles offers a custom field, widget and accompanying Ajax API views
for autocomplete-enabled selection of SODAR users in Django forms. The field
will handle providing appropriate choices according to the view context and user
permissions, also allowing for customization.

The recommended way to use the built-in user form field is by using the
``SODARUserChoiceField`` class found in ``projectroles.forms``. The field
extends Django's ``ModelChoiceField`` and takes most of the same keyword
arguments in its init function, with the exception of ``queryset``,
``to_field_name``, ``limit_choices_to`` and ``widget`` which will be overridden.

The init function also takes new arguments which are specified below:

- ``scope``: Scope of users to include (string)
    * ``all``: All users on the site
    * ``project``: Limit search to users in given project
    * ``project_exclude`` Exclude existing users of given project
- ``project``: Project object or project UUID string (optional)
- ``exclude``: List of User objects or User UUIDs to exclude (optional)
- ``forward``: Parameters to forward to autocomplete view (optional)
- ``url``: Autocomplete ajax class override (optional)
- ``widget_class``: Widget class override (optional)

Below is an example of the classes usage. Note that you can also define the
field as a form class member, but the ``project`` or ``exclude`` values are
not definable at that point. The following example assumes you are setting up
your project app form with an extra ``project`` argument.

.. code-block:: python

    from projectroles.forms import SODARUserChoiceField

    class YourForm(forms.ModelForm):
        class Meta:
            # ...
        def __init__(self, project, *args, **kwargs):
            # ...
            self.fields['user'] = SODARUserChoiceField(
                label='User',
                help_text='Select user for your thing here',
                required=True,
                scope='project',
                project=project,
                exclude=[unwanted_user]
            )

For more examples of usage of this field and its widget, see
``projectroles.forms``. If the field class does not suit your needs, you can also
retrieve the related widget to your own field with
``projectroles.forms.get_user_widget()``.

The following ``django-autocomplete-light`` and ``select2`` CSS and Javascript
links have to be added to the HTML template that includes the form with your
user selection field:

.. code-block:: django

    {% block javascript %}
      {{ block.super }}
      <!-- DAL for autocomplete widgets -->
      <script type="text/javascript" src="{% static 'autocomplete_light/jquery.init.js' %}"></script>
      <script type="text/javascript" src="{% static 'autocomplete_light/autocomplete.init.js' %}"></script>
      <script type="text/javascript" src="{% static 'autocomplete_light/vendor/select2/dist/js/select2.full.js' %}"></script>
      <script type="text/javascript" src="{% static 'autocomplete_light/select2.js' %}"></script>
    {% endblock javascript %}

    {% block css %}
      {{ block.super }}
      <!-- Select2 theme -->
      <link href="{% static 'autocomplete_light/vendor/select2/dist/css/select2.min.css' %}" rel="stylesheet" />
    {% endblock css %}

If using a customized widget with its own Javascript, include the corresponding
JS file instead of ``autocomplete_light/select2.js``. See the
``django-autocomplete-light`` documentation for more information on how to
customize your autocomplete-widget.

Markdown
--------

For fields supporting markdown, it is recommended to use the
``SODARPagedownWidget`` found in ``projectroles.models``.


Specific Views and Templates
============================

A few specific views/templates are expected to be implemented.

App Entry Point
---------------

As described in the Plugins chapter, an app entry point view is to be defined
in the ``ProjectAppPlugin``. This is **mandatory**.

The view **must** take a ``project`` URL kwarg which corresponds to a
``Project.sodar_uuid``.

For an example, see ``example_project_app.views.ExampleView`` and the associated
template.

Project Details Element
-----------------------

A sub-template to be included in the project details page (the project's "front
page" provided by projectroles, where e.g. overview of app content is shown).

Traditionally these files are called ``_details_card.html``, but you can name
them as you wish and point to the related template in the ``details_template``
variable of your plugin.

It is expected to have the content in a ``card-body`` container:

.. code-block:: django

   <div class="card-body">
     {# Content goes here #}
   </div>


Project Search Function and Template
====================================

If you want to implement search in your project app, you need to implement the
``search()`` function in your plugin as well as a template for displaying the
results.

.. hint::

   Implementing search *can* be complex. If you have access to the main SODAR
   repository, apps in that project might prove useful examples.

The search() Function
---------------------

See the signature of ``search()`` in
``projectroles.plugins.ProjectAppPluginPoint``. The arguments are as follows:

- ``search_terms``
    - One or more terms to be searched for (list of strings). Expected to be
      combined with OR operators in your search logic.
    - Multiple search terms or phrases containing whitespaces can be provided
      via the Advanced Search view.
- ``user``
    - User object for user initiating search
- ``search_type``
    - The type of object to search for (string, optional)
    - Used to restrict search to specific types of objects
    - You can specify supported types in the plugin's ``search_types`` list.
    - Examples: ``file``, ``sample``..
- ``keywords``
    - Special search keywords, e.g. "exact"
    - **NOTE:** Currently not implemented

.. note::

   Within this function, you are expected to verify appropriate access of the
   seaching user yourself!

.. warning::

    The old expected signature of providing a single ``search_term`` argument
    has been deprecated in v0.9 and will be removed in the next major release!

The return data is a dictionary, which is split by groups in case your app can
return multiple different lists for data. This is useful where e.g. the same
type of HTML list isn't suitable for all returnable types. If only returning one
type of data, you can just use e.g. ``all`` as your only category. Example of
the result:

.. code-block:: python

   return {
       'all': {                     # 1-N categories to be included
           'title': 'List title',   # Title of the result list to be displayed
           'search_types': [],      # Object types included in this category
           'items': []              # The actual objects returned
           }
       }

Search Template
---------------

Projectroles will provide your template context the ``search_results`` object,
which corresponds to the result dict of the aforementioned function. There are
also includes for formatting the results list, which you are encouraged to use.

Example of a simple results template, in case of a single ``all`` category:

.. code-block:: django

   {% if search_results.all.items|length > 0 %}

     {# Include standard search list header here #}
     {% include 'projectroles/_search_header.html' with search_title=search_results.all.title result_count=search_results.all.items|length %}

     {# Set up a table with your results #}
     <table class="table table-striped sodar-card-table sodar-search-table" id="sodar-ff-search-table">
       <thead>
         <tr>
           <th>Name</th>
           <th>Some Other Field</th>
         </tr>
      </thead>
      <tbody>
        {% for item in search_results.all.items %}
          <tr>
            <td>
              <a href="#link_to_somewhere_in your_app">{{ item.name }}</a>
            </td>
            <td>
              {{ item.some_other_field }}
            </td>
          </tr>
        {% endfor %}
      </tbody>
    </table>

    {# Include standard search list footer here #}
    {% include 'projectroles/_search_footer.html' %}

  {% endif %}


Tour Help
=========

SODAR Core uses `Shepherd <https://shipshapecode.github.io/shepherd/docs/welcome/>`_
to present an optional interactive tour for a rendered page. To enable the tour
in your template, set it up inside the ``javascript`` template block. Within an
inline javascript strucure, set the ``tourEnabled`` variable to ``true`` and add
steps according to the `Shepherd documentation <https://shipshapecode.github.io/shepherd>`_.

Example:

.. code-block:: django

    {% block javascript %}
      {{ block.super }}

      {# Tour content #}
      <script type="text/javascript">
        tourEnabled = true;

        /* Normal step */
        tour.addStep('id_of_step', {
            title: 'Step Title',
            text: 'Description of the step',
            attachTo: '#some-element top',
            advanceOn: '.docs-link click',
            showCancelLink: true
        });

        /* Conditional step */
        if ($('.potentially-existing-element').length) {
            tour.addStep('id_of_another_step', {
                title: 'Another Title',
                text: 'Another description here',
                attachTo: '.potentially-existing-element right',
                advanceOn: '.docs-link click',
                showCancelLink: true
            });
        }

      </script>
    {% endblock javascript %}


.. warning::

    Make sure you call ``{{ block.super }}`` at the start of the declared
    ``javascript`` block or you will overwrite the site's default Javascript
    setup!


API Views
=========

API view usage in project apps is detailed in this section.

Rest API Views
--------------

To set up REST API views for project apps, it is recommended to use the base
SODAR API view classes and mixins found in ``projectroles.views_api``. These
set up the recommended authentication methods, versioning through accept headers
and project-based permission checks.

By default, the REST API views built on SODAR Core base classes support two
methods of authentication: Knox tokens and Django session auth. These can of
course be modified by overriding/extending the base classes.

For versioning we strongly recommend using accept header versioning, which is
what is supported by the SODAR Core base classes. For this, supply your custom
media type and version data using the corresponding ``SODAR_API_*`` settings.
For details on these, see :ref:`app_projectroles_settings`.

The base classes provide permission checks via SODAR Core project objects
similar to UI view mixins.

Base REST API classes without a project context can also be used in site apps.

API documentation for each available base class and mixin for REST API views can
be found in :ref:`app_projectroles_api_django`.

An example "hello world" REST API view for SODAR apps is available in
``example_project_app.views.HelloExampleProjectAPIView``.

.. note::

    Internal SODAR Core REST API views, specifically ones used in apps provided
    by the django-sodar-core package, use different media type and versioning
    from views to be implemented on your site. This is to prevent version number
    clashes and not require changes from your API when SODAR Core is updated.

    For implementing your own API views, make sure to use the ``SODARAPI*``
    base classes, **not** the ``CoreAPI`` classes. Similarly, in testing make
    sure to use the base class helpers of the site API instead of the core API.

Ajax API Views
--------------

To set up Ajax API views for the UI, it is recommended to use the base Ajax
view classes found in ``projectroles.views_ajax``. These views only support
Django session authentication by default, so Knox token authentication will not
work. Versioning is omitted. Base views without project permission checks can
also be used in site apps.

API documentation for the base classes Ajax API views can be found in
:ref:`app_projectroles_api_django`.

Example:

.. code-block:: python

    from projectroles.views_api import SODARBaseProjectAjaxView

    class ExampleAjaxAPIView(SODARBaseProjectAjaxView):

    permission_required = 'projectroles.view_project'

    def get(self, request):
        # ...


Serializers
-----------

Base serializers for SODAR Core based API views are available in
``projectroles.serializers``. They provide ``Project`` context where needed, as
well as setting default fields such as ``sodar_uuid`` which should be always
used in place of ``pk``.

API documentation for the base serializers can be found in
:ref:`app_projectroles_api_django`.


Removing a Project App
======================

Removing a project app from your Django site can be slightly more complicated
than removing a normal non-SODAR-supporting Django application. Following the
procedure detailed here you are able to cleanly remove a project app which has
been in use on your site.

The instructions apply to project apps you have created yourself as well as
project apps included in the django-sodar-core package, with the exception of
``projectroles`` which can not be removed from a SODAR based site.

.. warning::

    Make sure to perform these steps **in the order they are presented here**.
    Otherwise you may risk serious problems with your site functionality or your
    database!

.. note::

    Just in case, it is recommended to make a backup of your Django database
    before proceeding.

First you should delete all Timeline references to objects in your app. This is
not done automatically as, by design, the references are kept even after the
original objects are deleted. Go to the Django shell via management command
using ``shell`` or ``shell_plus`` and enter the following. Replace ``app_name``
with the name of your application as specified in its ``ProjectAppPlugin``.

.. code-block:: python

    from timeline.models import ProjectEvent
    ProjectEvent.objects.filter(app='app_name').delete()

Next you should delete existing database objects defined by the models in your
app. This is also most easily done via the Django shell. Example:

.. code-block:: python

    from yourapp.models import YourModel
    YourModel.objects.all().delete()

After the objects have been deleted, reset the database migrations of your
application.

.. code-block:: console

    $ ./manage.py migrate yourapp zero

Once this has been executed successfully, you should delete the plugin object
for your application. Returning to the Django shell, type the following:

.. code-block:: python

    from djangoplugins.models import Plugin
    Plugin.objects.get(name='app_name').delete()

Finally, you should remove the references to the removed app in the Django
configuration.

App dependency in ``config/settings/base.py``:

.. code-block:: python

    LOCAL_APPS = [
    # The app you are removing
    'yourapp.apps.YourAppConfig',
    # ...
    ]

App URL patterns in ``config/urls.py``:

.. code-block:: python

    urlpatterns = [
        # Your app's URLs
        url(r'^yourapp/', include('yourapp.urls')),
        # ...
    ]

Once you have performed the aforementioned database operations and deployed a
version of your Django site with the application dependency and URL patterns
removed, the project app should be cleanly removed from your site.


TODO
====

- Naming conventions
- Examples of recurring template styles (e.g. forms)
.. _app_projectroles:


Projectroles App
^^^^^^^^^^^^^^^^

The ``projectroles`` app is the base app for building a
:term:`SODAR Core based Django site<SODAR Core Based Site>`. It provides a
ramework for project access management, dynamic content retrieval, models and
tools for SODAR-compatible apps plus a default template and CSS layout.

Other Django apps which intend to use aforementioned functionalities depend on
projectroles. While inclusion of other SODAR Core apps can be optional, having
projectroles installed is **mandatory** for working with the SODAR Core project
and app structure.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

    Basics <app_projectroles_basics>
    Integration <app_projectroles_integration>
    Settings <app_projectroles_settings>
    Usage <app_projectroles_usage>
    Customization <app_projectroles_custom>
    REST API Documentation <app_projectroles_api_rest>
    Django API Documentation <app_projectroles_api_django>
.. _app_adminalerts:


Adminalerts App
^^^^^^^^^^^^^^^

The ``adminalerts`` site app enables system administrators to display site-wide
messages to all users with an expiration date.


Basics
======

The app displays un-dismissable small alerts on the top of page content to all
users. They can be used to e.g. warn users of upcoming downtime or highlight
recently deployed changes.

Upon creation, an expiration date is set for each alert. Alerts can also be
freely enabled, disabled or deleted by superuser on the app UI. Additional
information regarding an alert can be provided with Markdown syntax and viewed
on a separate details page.


Installation
============

.. warning::

    To install this app you **must** have the ``django-sodar-core`` package
    installed and the ``projectroles`` app integrated into your Django site.
    See the :ref:`projectroles integration document <app_projectroles_integration>`
    for instructions.

Django Settings
---------------

The adminalerts app is available for your Django site after installing
``django-sodar-core``. Add the app into ``THIRD_PARTY_APPS`` as follows:

.. code-block:: python

    THIRD_PARTY_APPS = [
        # ...
        'adminalerts.apps.AdminalertsConfig',
    ]

Optional Settings
-----------------

To alter default adminalerts app settings, insert the following **optional**
variables with values of your choosing:

.. code-block:: python

    # Adminalerts app settings
    ADMINALERTS_PAGINATION = 15    # Number of alerts to be shown on one page (int)

URL Configuration
-----------------

In the Django URL configuration file, add the following line under
``urlpatterns`` to include adminalerts URLs in your site.

.. code-block:: python

    urlpatterns = [
        # ...
        url(r'^alerts/admin/', include('adminalerts.urls')),
    ]

Migrate Database and Register Plugin
------------------------------------

To migrate the Django database and register the adminalerts site app plugin,
run the following management command:

.. code-block:: console

    $ ./manage.py migrate

In addition to the database migration operation, you should see the following
output:

.. code-block:: console

    Registering Plugin for admimnalert.plugins.SiteAppPlugin


Usage
=====

When logged in as a superuser, you can find the :guilabel:`Admin Alerts` option
in your user dropdown menu in the top right corner of the site. This application
is only available for users with superuser status.

The initial view displays a list of alerts along with their active/inactive
status. For each alert there is a dropdown for updating or deleting the alert.
The :guilabel:`Status` column active/inactive badge doubles as a button for
toggling the alert status. Inactive alerts will not be shown to users.

.. figure:: _static/app_adminalerts/admin_alerts.png
    :align: center
    :scale: 60%

    Admin alert list

To create a new alert, click the :guilabel:`Create Alert` button. This presents
a simple form for creating a new alert. The following fields can be edited:

Message
    The message displayed to users on any page of the site.
Expiry Date
    The date when this alert will automatically expire and no longer be
    displayed.
Active
    Flag for the current state of the alert. Alerts can be activated or
    inactivated at any time, although activating an expired alert will not cause
    it to be displayed to users.
Require Auth
    If set true, this alert will only be shown to users logged in to the site.
    If false, it will also appear in the login screen as well as for anonymous
    users if allowed on the site.
Description
    A longer description, which can be accessed through the :guilabel:`Details`
    link in the alert element. Markdown syntax is supported.

.. figure:: _static/app_adminalerts/alert_form.png
    :align: center
    :scale: 60%

    Admin alert creation form

An example of an active alert along with the alert detail page can be seen in
the following screenshot.

.. figure:: _static/app_adminalerts/alert_example.png
    :align: center
    :scale: 55%

    Admin alert example and details.. include:: ../../CODE_OF_CONDUCT.rst

.. _app_bgjobs_usage:


Bgjobs Usage
^^^^^^^^^^^^

Usage instructions for the ``bgjobs`` app are detailed in this document.

**TODO**.. _app_projectroles_api_rest:


Projectroles REST API Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This document contains the HTTP REST API documentation for the ``projectroles``
app. The provided API enpoints allow project and role operations through
HTTP API calls in addition to the GUI.


API Usage
=========

Usage of the REST API is detailed in this section. These instructions also apply
to REST APIs in any other application within SODAR Core and are recommended
as guidelines for API development in your SODAR Core based Django site.

Authentication
--------------

The API supports authentication through Knox authentication tokens as well as
logging in using your SODAR username and password. Tokens are the recommended
method for security purposes.

For token access, first retrieve your token using the :ref:`app_tokens`. Add
the token in the ``Authorization`` header of your HTTP request as follows:

.. code-block:: console

    Authorization: token 90c2483172515bc8f6d52fd608e5031db3fcdc06d5a83b24bec1688f39b72bcd

Versioning
----------

The SODAR Core REST API uses accept header versioning. While specifying the
desired API version in your HTTP requests is optional, it is
**strongly recommended**. This ensures you will get the appropriate return data
and avoid running into unexpected incompatibility issues.

To enable versioning, add the ``Accept`` header to your request with the
following media type and version syntax. Replace the version number with your
expected version.

.. code-block:: console

    Accept: application/vnd.bihealth.sodar-core+json; version=0.9.1

.. note::

    The media type and version for internal SODAR Core apps are by design
    intended to be different to applications implemented in your Django site.
    Only use the aforementioned values when calling REST API views in
    projectroles or other applications installed from the django-sodar-core
    package.

Model Access and Permissions
----------------------------

Objects in SODAR Core API views are accessed through their ``sodar_uuid`` field.
This is strongly recommended for views implemented in your Django site as well,
as using a field such as ``pk`` may reveal internal database details to users as
well as be incompatible if e.g. mirroring roles between multiple SODAR Core
sites.

In the remainder of this document and other REST API documentation, *"UUID"*
refers to the ``sodar_uuid`` field of each model unless otherwise noted.

For permissions the API uses the same rules which are in effect in the SODAR
Core GUI. That means you need to have appropriate project access for each
operation.

Project Type Restriction
------------------------

IF you want to explicitly restrict access for your API view to a specific
project type, you can set the ``project_type`` attribute of your class to either
``PROJECT_TYPE_PROJECT`` or ``PROJECT_TYPE_CATEGORY`` as defined in
``SODAR_CONSTANTS``. A request to the view using the wrong project type will
result in a ``403 Not Authorized`` response, with the reason displayed in the
``detail`` view.

This works with any API view using ``SODARAPIProjectPermission`` as its
permission class, which includes ``SODARAPIBaseProjectMixin`` and
``SODARAPIGenericProjectMixin``. An example is shown below.

.. code-block:: python

    from rest_framework.generics import RetrieveAPIView
    from projectroles.models import SODAR_CONSTANTS
    from projectroles.views_api import CoreAPIGenericProjectMixin

    class YourAPIView(SODARAPIGenericProjectMixin, RetrieveAPIView):
        # ...
        project_type = SODAR_CONSTANTS['PROJECT_TYPE_PROJECT']

Return Data
-----------

The return data for each request will be a JSON document unless otherwise
specified.

If return data is not specified in the documentation of an API view, it will
return the appropriate HTTP status code along with an optional ``detail`` JSON
field upon a successfully processed request.

For creation views, the ``sodar_uuid`` of the created object is returned
along with other object fields.


API Views
=========

.. currentmodule:: projectroles.views_api

.. autoclass:: ProjectListAPIView

.. autoclass:: ProjectRetrieveAPIView

.. autoclass:: ProjectCreateAPIView

.. autoclass:: ProjectUpdateAPIView

.. autoclass:: RoleAssignmentCreateAPIView

.. autoclass:: RoleAssignmentUpdateAPIView

.. autoclass:: RoleAssignmentDestroyAPIView

.. autoclass:: RoleAssignmentOwnerTransferAPIView

.. autoclass:: ProjectInviteListAPIView

.. autoclass:: ProjectInviteCreateAPIView

.. autoclass:: ProjectInviteRevokeAPIView

.. autoclass:: ProjectInviteResendAPIView

.. autoclass:: UserListAPIView

.. autoclass:: CurrentUserRetrieveAPIView
.. _app_timeline:


Timeline App
^^^^^^^^^^^^

The ``timeline`` app enables the developer of a SODAR Core based site to log
project related or site-wide events and link objects (both existing and deleted)
to those events.

Unlike the standard Django object history accessible in the admin
site, these events are not restricted to creation/modification of objects in the
Django database, but can concern any user-triggered activity.

The events can also have multiple temporal status states in case of e.g. events
requiring async requests.

The app provides front-end views to list timeline events for projects,
categories, objects and site-wide actions. Also included is a backend API for
saving desired activity as timeline events. For details on how to use these, see
the :ref:`timeline usage documentation <app_timeline_usage>`.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   Installation <app_timeline_install>
   Usage <app_timeline_usage>
   Django API Documentation <app_timeline_api_django>
.. _changelog:

.. include:: ../../CHANGELOG.rst
