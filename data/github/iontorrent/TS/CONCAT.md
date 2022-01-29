# TMAP - torrent mapping alignment program

##  General Notes 

TMAP is a fast and accurate alignment software for short and long nucleotide sequences produced by next-generation sequencing technologies.

* The latest TMAP is unsupported.  To use a supported version, please see the TMAP version associated with a Torrent Suite release below.

*  Get the latest source code: 
    <pre lang="bash"><code>git clone git://github.com/iontorrent/TMAP.git
    cd TMAP
    git submodule init
    git submodule update
    </code></pre>
* To download a specific version, please get the latest source code, and then checkout the specific version using the tag listed below:
    <pre lang="bash"><code>git checkout -b tag "tag name below"</code></pre>
    For example: <pre lang="bash"><code>git checkout -b 3.0.1 tmap.3.0.1
    git submodule update
    </code></pre>
    Below is the list of tags associated with a specific Torrent Suite release:<table>
    <tr><td>Torrent Suite 3.0:</td><td>tmap.3.0.1</td></tr>
    <tr><td>Torrent Suite 2.2</td><td>tmap.0.3.7</td></tr>
    <tr><td>Torrent Suite 2.0.1</td><td>tmap.0.2.3</td></tr>
    <tr><td>Torrent Suite 2.0</td><td>tmap.0.2.3</td></tr>
    <tr><td>Torrent Suite 1.5.1</td><td>tmap.0.1.3</td></tr>
    <tr><td>Torrent Suite 1.5</td><td>tmap.0.1.3</td></tr>
    <tr><td>Torrent Suite 1.4.1</td><td>tmap.0.0.28</td></tr>
    <tr><td>Torrent Suite 1.4</td><td>tmap.0.0.25</td></tr>
    <tr><td>Torrent Suite 1.3</td><td>tmap.0.0.19</td></tr>
    <tr><td>Torrent Suite 1.2</td><td>tmap.0.0.9</td></tr>
    </tr>
    </table>
*  See the latest manual: http://github.com/iontorrent/TMAP/blob/master/doc/tmap-book.pdf


##  Pre-requisites
1. Compiler (required):
  The compiler and system must support SSE2 instructions.  

##  To Install

1. Compile TMAP:
  <pre lang="bash"><code>sh autogen.sh && ./configure && make</code></pre>
2. Install
  <pre lang="bash"><code>make install</code></pre>

##  Optional Installs

### TCMalloc (optional)
  TMAP will run approximately 15% faster using the tcmalloc memory allocation
  implementation.  To use tcmalloc, install the Google performance tools:
  http://code.google.com/p/google-perftools
  
  If you have previously compiled TMAP, execute the following command:
  <pre lang="bash"><code>make distclean && sh autogen.sh && ./configure && make clean && make</code></pre>
  After installation, execute the following command:
  <pre lang="bash"><code>sh autogen.sh && ./configure && make clean && make</code></pre>
  The performance improve should occur when using multiple-threads.

##  Developer Notes

There are a number of areas for potential improvement within TMAP for those
that are interested; they will be mentioned here.  A great way to find places
where the code can be improved is to use Google's performance tools:
  http://code.google.com/p/google-perftools
This includes a heap checker, heap profiler, and cpu profiler.  Examining 
performance on large genomes (hg19) is recommended.

### Smith Waterman extensions
  Currently, each hit is examined with Smith Waterman (score only), which
   re-considers the portion of the read that matched during seeding.  We need
   only re-examine the portion of the read that is not matched during seeding.
   This could be tracked during seeding for the Smith Waterman step, though 
   the merging of hits from each algorithm could be complicated by this step.
   Nonetheless, this would improve the run time of the program, especially for
   high-quality data and/or longer reads (>200bp).

### Smith Waterman vectorization
  The vectorized (SSE2) Smith Waterman implemented supports an combination of
    start and end soft-clipping.  To support any type of soft-clipping, some 
    performance trade-offs needed to be made.  In particular, 16-bit integers
	are stored in the 128-bit integers, giving only 8 bytes/values per 128-bit 
    integer.  This could be improved to 16 bytes/values per 128-bit integer by
    using 8-bit integers.  This would require better overflow handling.  Also,
    handling negative infinity for the Smith Waterman initialization would be
    difficult.  Nonetheless, this could significantly improve the performance of
    the most expensive portion of the program.

### Best two-stage mapping
  There is no current recommendation for the best settings for two-stage 
    mapping, which could significantly decrease the running time.  A good 
	project would be to optimize these settings.

### Mapping quality calibration
  The mapping quality is sufficiently calibrated, but can always be improved,
    especially for longer reads.  This is a major area for improvement.

### Better support for paired ends/mate pairs
  There is minimal support for paired ends/mate pairs, which relies on knowing
    a prior the parameters for the insert size distribution.  The insert size 
	could be trained on a subset of the given input data.

### Speeding up lookups in the FM-index/BWT.
  Further implementation improvements or parameter tuning could be made to make
    the lookups in the FM-index/BWT faster.  This includes the occurrence 
	interval, the suffix array interval, and the k-mer occurence hash.  Caching
	these results may also make sense when examining the same sub-strings across
	multiple algorithms.  Speed improvements have already been made to BWA and 
	could be relevant here:
	  http://github.com/RoelKluin/bwa

### Dynamic split read mapping
  It is important to detect Structural Variation (SV), as well as finding splice 
    junctions for RNA-seq.  Support for returning more than one alignment, where
	these alignments do not significantly overlap in terms of which bases they
	consume in the query, could be included.  For example, a 400bp read could span
	a SV breakpiont, with the first 100bp on one side of the breakpoint and the 
	second 300bp on the other.  Currently (with full soft-clipping options turned 
    on), we may produce two alignments for the two parts of the query. Nonetheless,
	the "choice" algorithm will choose the one with the best alignment score 
    (typically the 300bp one), and so only one alignment will be present in the SAM
	file.  A better strategy would be to search for pairs (triples, etc.) of 
	alignments that do not significantly overlap in the query (i.e. consume the same
	query bases).  This would directly find SVs as well as other types of variant
	requiring split read mapping.

### Representative repetitive hits
  If a seeding algorithm finds a large occurence interval that it will save, it 
    could save one of the occurrences (random) as a representative hit for the 
	repetitive interval.  This representative hit could be aligned with Smith-Waterman
	and its alignment score could be compared to the other hits.  If its score is
	better, than the read could be flagged as repetitive and "unmapped".  The 
	algorithm would need to be careful that the repetitive hit is not contained 
	within the returned non-repeititve hits, as to cause many reads to be unmapped.
# Galaxy integration for TMAP

TMAP can be run under galaxy.  This feature is provided for users familiar with 
integrating other tools within Galaxy.  There is no warranty or support offered
for this feature.

## Installation

### Copy the Wrapper Scripts
Place the following files into the "galaxy/tools/sr_mapping/" directory: 
* tmap_wrapper.py
* tmap_wrapper.xml

### Copy and Update the TMAP Index Locations 
If you have pre-built TMAP indexes (recommended), rename the tmap_index.loc.sample
file "tmap_index.loc" and place it int the "galaxy/tools-data/" directory.

### Let Galaxy Know about the TMAP Tool
Modify the "" file to include the TMAP wrapper by adding the following line
to "galaxy/tool_conf.xml" in the "NGS: Mapping" section (search for "bwa" for example):
* <pre lang="xml"><code><tool file="sr_mapping/tmap_wrapper.xml" /></code></pre>

## Running TMAP within Galaxy
The TMAP page includes descriptions of the supported options.  For global, flowspace, 
pairing, and algorithm options, enter them in the given text boxes exactly the same 
as you would on the command line.  If you encounter problems, please review the source
code or contact the galaxy help list.
# Guidelines for plugin development alongside tools

## Rules
- Changes to metric definition require a major revision change for your plugin (major.minor.bug)
- At a minimum, define such major changes for posterity in a README.md file in your plugin

## Updates
- Update to use PluginMixin (then call self.init_plugin())
- Remove block_reshape (and analogs) definitions and leverage block_reshape.py
- Leverage lanes.py for iteration through lanes and other handy multilane tools
- Remove dependencies on average.py, as it is deprecated.# Guidelines for plugin development alongside tools

## Rules
- Changes to metric definition require a major revision change for your plugin (major.minor.bug)
- At a minimum, define such major changes for posterity in a README.md file in your plugin

## Updates
- Update to use PluginMixin (then call self.init_plugin())
- Remove block_reshape (and analogs) definitions and leverage block_reshape.py
- Leverage lanes.py for iteration through lanes and other handy multilane tools
- Remove dependencies on average.py, as it is deprecated.# Info
- This readme file is meant to contain update information to images and metrics that go along with the plugin.
- **At a minimum, major version changes must be made when metric definitions change, and they should be described in this file!**
- Further detail can be handy to retain.  Previously, some of this recent information was included in the docstrings for plugin classes.

---

# Major Revisions

## **6.0.0** [Feb 2019]
- After realization that 560 chips and others did not adhere to the new Valkyrie rules and confirmation that the 'DynamicRange' field is non-deterministic and saves the DR at save time into this field.....so not useful.
- Therefore, reverted to DynamicRangeAfterBF.

## **5.0.X** [Jan 2019]
- Thought it was time to switch back to DynamicRange as the true metric.  Did not turn out to be the right decision.

## **4.0.0** [Sept 2019]
- The primary change here redefined our pixel offset metric definitions based on how we have been pulling in the wrong DR for some time now.  The new value is correctly pulled from the explog file (DynamicRangeAfterBF).
- Previous versions of the plugin still have this information (so we can retroactively correct on ChipDB) so that they can be corrected by the factor of ( DynamicRangeAfterBF / DynamicRange )
- Also created a host of new metrics (usually block-reshape based ones) with a new prefix of 'true_' which are now appropriately calculated due to exclusion of pinned pixels from calculations.  (This is particularly important for studies of the square defect on GX5/GX7 chips)

## **3.0.0** [June 2018]
- Added (radial) edge analysis capabilities for calibration metrics.

# Other Revisions of Note

## 3.3.0

- This was the first revision built to work on either RUO TS or the Valkyrie TS.# Guidelines for plugin development alongside tools

## Rules
- Changes to metric definition require a major revision change for your plugin (major.minor.bug)
- At a minimum, define such major changes for posterity in a README.md file in your plugin

## Updates
- Update to use PluginMixin (then call self.init_plugin())
- Remove block_reshape (and analogs) definitions and leverage block_reshape.py
- Leverage lanes.py for iteration through lanes and other handy multilane tools
- Remove dependencies on average.py, as it is deprecated.# Guidelines for plugin development alongside tools

## Rules
- Changes to metric definition require a major revision change for your plugin (major.minor.bug)
- At a minimum, define such major changes for posterity in a README.md file in your plugin

## Updates
- Update to use PluginMixin (then call self.init_plugin())
- Remove block_reshape (and analogs) definitions and leverage block_reshape.py
- Leverage lanes.py for iteration through lanes and other handy multilane tools
- Remove dependencies on average.py, as it is deprecated.# Guidelines for plugin development alongside tools

## Rules
- Changes to metric definition require a major revision change for your plugin (major.minor.bug)
- At a minimum, define such major changes for posterity in a README.md file in your plugin

## Updates
- Update to use PluginMixin (then call self.init_plugin())
- Remove block_reshape (and analogs) definitions and leverage block_reshape.py
- Leverage lanes.py for iteration through lanes and other handy multilane tools
- Remove dependencies on average.py, as it is deprecated.Bootstrap - Lightbox
==================
This is a plugin for Twitter Bootstrap that adds a lightbox that is based off the modal dialog

Project page and demo: http://jbutz.github.com/bootstrap-lightbox/
FieldSupport
===========

**Do not modify plugins inside the rndplugins directory!**
Plugins should be modified upsteam and then pulled into this repo with the `fab update_plugins` command. 
Plugins can check for the existence of the `TSP_LIMIT_OUTPUT` environment variable if needed. This variable is only set
by the FieldSupport runtime. Plugins may need to run in a limited output mode inside this plugin to keep
the size of the resulting archive small (<10MB).

Updating Plugins
-----------------
Command: `fab update_plugins`

You will need the python fabric library installed. See [http://www.fabfile.org/](fabfile.org)

# Torrent Browser 

## Internalization Support - i18n

### Key Areas
* Dangjo Templates
* Client-Side JS and static files
* Django Server-Side 
    * Views, Models, APIS
    * tbd. 


### Configuration

#### Django Settings - iondb.settings

    ...
    
    LANGUAGES = (
    # English
    # ('en', u'English'), # instead of 'en'
    # English United Stated
    ('en-us', u'English US'), # instead of 'en_US'
    # Russian
    ('ru', u'русский'), # instead of 'ru'
    # Simplified Chinese
    ('zh-cn', u'简体中1'), # instead of 'zh-CN'
    #Traditional Chinese
    # ('zh-tw', u'繁體中文'), # instead of 'zh-TW'
    )
    ...
    # http://django-docs.jarvis.itw/ref/settings.html?highlight=locale_paths#std:setting-LOCALE_PATHS
    LOCALE_PATHS = (
        '/opt/ion/iondb/locale',
        '/var/local/translations/locale' # placeholder for additional translations, e.g. from oem provided .deb
    )
    # If you set this to False, Django will make some optimizations so as not
    # to load the internationalization machinery.
    USE_I18N = True
    ...
    MIDDLEWARE_CLASSES = (
        'iondb.rundb.middleware.ChangeRequestMethodMiddleware',
        'django.middleware.common.CommonMiddleware',
        'iondb.rundb.middleware.DeleteSessionOnLogoutMiddleware',
        'django.contrib.sessions.middleware.SessionMiddleware',
        ...
        # 'django.middleware.locale.LocaleMiddleware',  #uncomment to allow end-user locale to dictate the language used
        'iondb.bin.startup_housekeeping.StartupHousekeeping'
    )

tbd..# [Uni-Form Markup](http://sprawsm.com/uni-form/) : Validation documentation


## Initialize the jQuery plugin

The following code will initialize the jQuery Validation plugin with the default options.

    <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4/jquery.min.js"></script>
    <script type="text/javascript" src="../js/uni-form-validation.jquery.js" charset="utf-8"></script>
    <script>
      $(function(){
        $('form.uniForm').uniform();
      });
    </script>

You may use a global object to hold site wide validation settings. To do this, you should
copy the jQuery.fn.uniform.defaults = {} object from the bottom of the validation javascript
file into a new file that you use throughout your site. You may then edit options there
globally, and will make the Uni-Form library easy to update in the future.

You may also initialize Uni-Form Validation with custom settings by passing a settings object
as a parameter when you call uniform(). 

    <script>
      $(function(){
        $('form.uniForm').uniform({
            prevent_submit : true,
            valid_class    : 'okGo'
        });
      });
    </script>

## Uni-Form Settings

* prevent_submit (false)
  Set this to true to prevent the form from submitting if there are outstanding
  errors in the form
* prevent_submit_callback (false)
  Supply a function here and it will be called instead of the internal handler.
  This function can return true to allow the form to proceed with the commit
* ask_on_leave (false)
  Set this to true to have the browser prompt if the visitor has made changes to
  the form, and then initialized a page unload without submitting the form
* on_leave_callback (false)
  Provide a function and it will be called instead of the internal method
* valid_class ('valid')
  CSS class name used for div.holder_class elements that have passed validation
* invalid_class ('invalid')
  CSS class name used for div.holder_class elements that have failed validation
* error_class ('error')
  Please note that both of these are applied by the validation script.
  You may wish to set them separately at the server perhaps.
* focused_class ('focused')
  CSS class name applied to the .holder_class of the current element
* holder_class ('ctrlHolder')
  CSS class name that you have used as the control holder class
* field_selector ('input, textarea, select')
  List of html elements that will be treated with Uni-Form highlighting and 
  validation (if enabled)
* default_value_color ("#AFAFAF")
  HEX color used to display the default data in the background of empty text inputs
  
## Validators

* required
* validateMinLength
* validateMin
* validateMaxLength
* validateMax
* validateSameAs
* validateEmail
* validateUrl
* validateNumber
* validateInteger
* validateAlpha
* validateAlphaNum
* validatePhrase
* validatePhone
* validateDate
* validateCallback

Validators what require a parameter, such as validateMinLength, take that parameter
as a class name following the validator in the format of _val-{value}_. 

## validateCallback

The validateCallback is a special validator. See the demo/callback.html file for example
use. It allows you to define a custom callback to an input without having to add a new
validator type to the library.

# [Uni-Form Markup](http://sprawsm.com/uni-form/) : Making forms as simple as 1,2,3

## Announcements:

* __Please note that the jQuery plugins no longer automatically initialize.__
  You must init them yourself with the code found in the section below 
  titled "How to use?"
  

## Copyright (c) 2010, Dragan Babic
   
   Permission is hereby granted, free of charge, to any person
   obtaining a copy of this software and associated documentation
   files (the "Software"), to deal in the Software without
   restriction, including without limitation the rights to use,
   copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following
   conditions:
   
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
   OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
   HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
   WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
   OTHER DEALINGS IN THE SOFTWARE.


## About Uni–Form 

Uni-Form is a framework that standardizes form markup and styles it with CSS 
giving you two most widely used layout options to choose from. Anyone can get nice 
looking, well structured, highly customizable, accessible and usable forms. To put 
it simply: it makes a developer's life a lot easier. 

* [Uni-Form Homepage](http://sprawsm.com/uni-form/)
* [Support at Get Satisfaction](http://getsatisfaction.com/uni-form)
* [GitHub repository]()

## How to Use? 

First thing you need to do is to link up the necessary files: 

1.  Link to the main CSS file
    
        <link href="path/to/file/uni-form.css" media="all" rel="stylesheet"/>
    
1.  Link to the Uni–Form style CSS file
    
        <link href="path/to/file/default.uni-form.css" media="all" rel="stylesheet"/>
    
1.  Optionally you'll want to link up jQuery and Uni–Form jQuery files if you'd 
    like Uni–Form to highlight the form rows on focus (it's a usability aid): 
      
        <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.4/jquery.min.js"></script>
        <script type="text/javascript" src="path/to/file/uni-form.jquery.js"></script>
    
1.  You may also want to try out the version of the Uni–Form jQuery plugin that
    supports client side validation, in case replace the regular plugin this this:
    
        <script type="text/javascript" src="path/to/file/uni-form-validation.jquery.js"></script>

1. Please note that this plugin no longer automatically initialize the Uni–Form plugin.
   You must do this yourself, by adding this snippet after you have included
   both jQuery and the plugin you have chosen:
   
       <script type="text/javascript">
        $(function(){
          $('form.uniForm').uniform();
        });
       </script>


Now that you're all set up, all you need to do is add form fields that are formatted
with Uni–Form markup so the CSS and JavaScript will get the “hooks” they need. These
chunks of code are called “units” and all available units can be found within the 
file called fauxform.html that is included in this package. 

Feel free to extend Uni–Form with units of your own and share. 


## Styles 

As of v1.4 Uni–Form supports styles. These are separate CSS files that contain the
presentation aspect of your form (considering that uni-form.css) contains the 
layout. Style CSS files should be used to control how your form looks, spacing… 

Sharing styles is encouraged, and by default Uni–Form is shipped with three: 

 * Default
 * Blue 
 * Dark 
    
Consider these a starting point for making your own. 

## Options and Layout Control 

Uni–Form by default has two form layouts: default and inline. This is controlled 
by adding (or removing) a CSS class .inlineLabels to the fieldset element. 

There is another option in regards to the layout and it concerns what is referred 
to as "multifields". These are fields that contain multiple inputs per unit and 
are usually used for checkboxes and radio buttons. Each layout supports an 
alternate multifield layout. This is achieved by adding (or removing) a CSS class
.alternate to the ul element. 


## Events

Triggering an error event on the form fields will apply the error
class to the controller and overwrite the supplied description of that
controller with the error text, an example would be:

    $(selector).trigger('error',['an error occured']);

Subsequent calls to success on the form field will remove the error
and replace the error text with the originally supplied description,
an example:

    $(selector).trigger('success');

----------------------------------------------------------------------------------

## Form Validation

Uni–Form can be used with the included uni-form-validation.js file for client
side validation. This is accomplished by using class names on the form elements
to trigger validation rules on blur(). It must be noted that these validation rules
should be used to supplement a server side solution.

Required element, cannot be empty:

    <input type="text" class="textInput required" />

Integer with value greater than or equal to 8:

    <input type="text" class="textInput validateInteger validateMin val-8" />

### Available validators:

* required
* validateMinLength
* validateMin
* validateMaxLength
* validateMax
* validateSameAs
* validateEmail
* validateUrl
* validateNumber
* validateInteger
* validateAlpha
* validateAlphaNum
* validatePhrase
* validatePhone
* validateDate
* validateCallback

Validators what require a parameter, such as validateMinLength, take that parameter
as a class name following the validator in the format of _val-{value}_. 



## Give respect and get it back.
[Twitter Bootstrap](http://twitter.github.com/bootstrap) [![Build Status](https://secure.travis-ci.org/twitter/bootstrap.png)](http://travis-ci.org/twitter/bootstrap)
=================

Bootstrap is a sleek, intuitive, and powerful front-end framework for faster and easier web development, created and maintained by [Mark Otto](http://twitter.com/mdo) and [Jacob Thornton](http://twitter.com/fat).

To get started, checkout http://getbootstrap.com!



Quick start
-----------

Clone the repo, `git clone git://github.com/twitter/bootstrap.git`, [download the latest release](https://github.com/twitter/bootstrap/zipball/master), or install with twitter's [Bower](http://twitter.github.com/bower): `bower install bootstrap`.



Versioning
----------

For transparency and insight into our release cycle, and for striving to maintain backward compatibility, Bootstrap will be maintained under the Semantic Versioning guidelines as much as possible.

Releases will be numbered with the following format:

`<major>.<minor>.<patch>`

And constructed with the following guidelines:

* Breaking backward compatibility bumps the major (and resets the minor and patch)
* New additions without breaking backward compatibility bumps the minor (and resets the patch)
* Bug fixes and misc changes bumps the patch

For more information on SemVer, please visit http://semver.org/.



Bug tracker
-----------

Have a bug? Please create an issue here on GitHub that conforms with [necolas's guidelines](https://github.com/necolas/issue-guidelines).

https://github.com/twitter/bootstrap/issues



Twitter account
---------------

Keep up to date on announcements and more by following Bootstrap on Twitter, [@TwBootstrap](http://twitter.com/TwBootstrap).



Blog
----

Read more detailed announcements, discussions, and more on [The Official Twitter Bootstrap Blog](http://blog.getbootstrap.com).



Mailing list
------------

Have a question? Ask on our mailing list!

twitter-bootstrap@googlegroups.com

http://groups.google.com/group/twitter-bootstrap



IRC
---

Server: irc.freenode.net

Channel: ##twitter-bootstrap (the double ## is not a typo)



Developers
----------

We have included a makefile with convenience methods for working with the Bootstrap library.

+ **dependencies**
Our makefile depends on you having recess, connect, uglify.js, and jshint installed. To install, just run the following command in npm:

```
$ npm install recess connect uglify-js jshint -g
```

+ **build** - `make`
Runs the recess compiler to rebuild the `/less` files and compiles the docs pages. Requires recess and uglify-js. <a href="http://twitter.github.com/bootstrap/extend.html#compiling">Read more in our docs &raquo;</a>

+ **test** - `make test`
Runs jshint and qunit tests headlessly in [phantomjs](http://code.google.com/p/phantomjs/) (used for ci). Depends on having phantomjs installed.

+ **watch** - `make watch`
This is a convenience method for watching just Less files and automatically building them whenever you save. Requires the Watchr gem.



Contributing
------------

Please submit all pull requests against *-wip branches. If your unit test contains javascript patches or features, you must include relevant unit tests. Thanks!



Authors
-------

**Mark Otto**

+ http://twitter.com/mdo
+ http://github.com/markdotto

**Jacob Thornton**

+ http://twitter.com/fat
+ http://github.com/fat



Copyright and license
---------------------

Copyright 2012 Twitter, Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this work except in compliance with the License.
You may obtain a copy of the License in the LICENSE file, or at:

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
# Contributing to Bootstrap

Looking to contribute something to Bootstrap? **Here's how you can help.**



## Reporting issues

We only accept issues that are bug reports or feature requests. Bugs must be isolated and reproducible problems that we can fix within the Bootstrap core. Please read the following guidelines before opening any issue.

1. **Search for existing issues.** We get a lot of duplicate issues, and you'd help us out a lot by first checking if someone else has reported the same issue. Moreover, the issue may have already been resolved with a fix available.
2. **Create an isolated and reproducible test case.** Be sure the problem exists in Bootstrap's code with a [reduced test cases](http://css-tricks.com/reduced-test-cases/) that should be included in each bug report.
3. **Include a live example.** Make use of jsFiddle or jsBin to share your isolated test cases.
4. **Share as much information as possible.** Include operating system and version, browser and version, version of Bootstrap, customized or vanilla build, etc. where appropriate. Also include steps to reproduce the bug.



## Key branches

- `master` is the latest, deployed version.
- `gh-pages` is the hosted docs (not to be used for pull requests).
- `*-wip` is the official work in progress branch for the next release.



## Notes on the repo

As of v2.0.0, Bootstrap's documentation is powered by Mustache templates and built via `make` before each commit and release. This was done to enable internationalization (translation) in a future release by uploading our strings to the [Twitter Translation Center](http://translate.twttr.com/). Any edits to the docs should be first done in the Mustache files and then recompiled into the HTML.



## Pull requests

- Try to submit pull requests against the latest `*-wip` branch for easier merging
- Any changes to the docs must be made to the Mustache templates, not just the compiled HTML pages
- CSS changes must be done in .less files first, never just the compiled files
- If modifying the .less files, always recompile and commit the compiled files bootstrap.css and bootstrap.min.css
- Try not to pollute your pull request with unintended changes--keep them simple and small
- Try to share which browsers your code has been tested in before submitting a pull request



## Coding standards: HTML

- Two spaces for indentation, never tabs
- Double quotes only, never single quotes
- Always use proper indentation
- Use tags and elements appropriate for an HTML5 doctype (e.g., self-closing tags)



## Coding standards: CSS

- Adhere to the [Recess CSS property order](http://markdotto.com/2011/11/29/css-property-order/)
- Multiple-line approach (one property and value per line)
- Always a space after a property's colon (.e.g, `display: block;` and not `display:block;`)
- End all lines with a semi-colon
- For multiple, comma-separated selectors, place each selector on it's own line
- Attribute selectors, like `input[type="text"]` should always wrap the attribute's value in double quotes, for consistency and safety (see this [blog post on unquoted attribute values](http://mathiasbynens.be/notes/unquoted-attribute-values) that can lead to XSS attacks).



## Coding standards: JS

- No semicolons
- Comma first
- 2 spaces (no tabs)
- strict mode
- "Attractive"



## License

By contributing your code, you agree to license your contribution under the terms of the APLv2: https://github.com/twitter/bootstrap/blob/master/LICENSE
## 2.0 BOOTSTRAP JS PHILOSOPHY
These are the high-level design rules which guide the development of Bootstrap's plugin apis.

---

### DATA-ATTRIBUTE API

We believe you should be able to use all plugins provided by Bootstrap purely through the markup API without writing a single line of javascript.

We acknowledge that this isn't always the most performant and sometimes it may be desirable to turn this functionality off altogether. Therefore, as of 2.0 we provide the ability to disable the data attribute API by unbinding all events on the body namespaced with `'data-api'`. This looks like this:

    $('body').off('.data-api')

To target a specific plugin, just include the plugins name as a namespace along with the data-api namespace like this:

    $('body').off('.alert.data-api')

---

### PROGRAMATIC API

We also believe you should be able to use all plugins provided by Bootstrap purely through the JS API.

All public APIs should be single, chainable methods, and return the collection acted upon.

    $(".btn.danger").button("toggle").addClass("fat")

All methods should accept an optional options object, a string which targets a particular method, or null which initiates the default behavior:

    $("#myModal").modal() // initialized with defaults
    $("#myModal").modal({ keyboard: false }) // initialized with now keyboard
    $("#myModal").modal('show') // initializes and invokes show immediately afterqwe2

---

### OPTIONS

Options should be sparse and add universal value. We should pick the right defaults.

All plugins should have a default object which can be modified to effect all instance's default options. The defaults object should be available via `$.fn.plugin.defaults`.

    $.fn.modal.defaults = { … }

An options definition should take the following form:

    *noun*: *adjective* - describes or modifies a quality of an instance

examples:

    backdrop: true
    keyboard: false
    placement: 'top'

---

### EVENTS

All events should have an infinitive and past participle form. The infinitive is fired just before an action takes place, the past participle on completion of the action.

    show | shown
    hide | hidden

---

### CONSTRUCTORS

Each plugin should expose it's raw constructor on a `Constructor` property -- accessed in the following way:


    $.fn.popover.Constructor

---

### DATA ACCESSOR

Each plugin stores a copy of the invoked class on an object. This class instance can be accessed directly through jQuery's data API like this:

    $('[rel=popover]').data('popover') instanceof $.fn.popover.Constructor

---

### DATA ATTRIBUTES

Data attributes should take the following form:

- data-{{verb}}={{plugin}} - defines main interaction
- data-target || href^=# - defined on "control" element (if element controls an element other than self)
- data-{{noun}} - defines class instance options

examples:

    // control other targets
    data-toggle="modal" data-target="#foo"
    data-toggle="collapse" data-target="#foo" data-parent="#bar"

    // defined on element they control
    data-spy="scroll"

    data-dismiss="modal"
    data-dismiss="alert"

    data-toggle="dropdown"

    data-toggle="button"
    data-toggle="buttons-checkbox"
    data-toggle="buttons-radio"Bootstrap Modal v2.0
=============

See live demo [here](http://jschr.github.com/bootstrap-modal/).

Extends Bootstrap's native modals to provide additional functionality. Introduces a **ModalManager** class that operates behind the scenes to handle multiple modals by listening on their events. 

A single ModalManager is created by default on body and can be accessed through the jQuery plugin interface.

    $('body').modalmanager('loading');

Bootstrap-Modal can be used as a replacement for Bootstrap's Modal class or as a patch to the library.

Overview
-----------

+ Backwards compatible
+ Responsive
+ Stackable
+ Full width
+ Load content via AJAX
+ Disable background scrolling

Installation 
-----------
+ Include `css/bootstrap-modal.css` after the main bootstrap css files.
+ Include `js/bootstrap-modalmanager.js` and `js/bootstrap-modal.js` after the main bootstrap js files.

	<link href="css/bootstrap.css" rel="stylesheet" />
	<link href="css/bootstrap-responsive.css" rel="stylesheet" />
 	<link href="css/bootstrap-modal.css" rel="stylesheet" />

 	<script src="js/bootstrap.js"></script>
 	<script src="js/bootstrap-modalmanager.js"></script>
 	<script src="js/bootstrap-modal.js"></script>

Options
-----------

In addition to the standard bootstrap options, you now have access to the following options

**Modal**

+ **width**
Set the inital width of the modal.

+ **height**
Set the inital height of the modal.

+ **maxHeight**
Set the max-height of the modal-body.

+ **loading**
Toggle the loading state.

+ **spinner**
Provide a custom image or animation for the loading spinner.

+ **consumeTab**
Used to enable tabindexing for modals with `data-tabindex`. This is set to true by default.

+ **focusOn**
The element or selector to set the focus to once the modal is shown.

+ **attentionAnimation**
Set the animation used by the `attention` method. Any animation in [animate.css](http://daneden.me/animate/) is supported but only the *shake* animation is included by default.

+ **modalOverflow**
Set this property to true for modals with highly dynamic content. This will force the modal to behave as if it is larger than the viewport.

+ **manager**
Set the modal's manager. By default this is set to the `GlobalModalManager` and will most likely not need to be overridden.

**ModalManager**

+ **loading**
Toggle the loading state.

+ **backdropLimit**
Limit the amount of backdrops that will appear on the page at the same time.

+ **spinner**
Provide a custom image or animation for the loading spinner.

Disable Background Scrolling
-----------

If you want to prevent the background page from scrolling (see [demo](http://jschr.github.com/bootstrap-modal/) for example) you must wrap the page contents in a `<div class="page-container">`. For example:

	<body>
		<div class="page-container">
			<div class="navbar navbar-fixed-top">...</div>
			<div class="container">...</div>
		</div>
	</body>

The reason for doing this instead of just simply setting `overflow: hidden` when a modal is open is to avoid having the page shift as a result of the scrollbar appearing/disappearing. This also allows the document to be scrollable when there is a tall modal but only to the height of the modal, not the entire page.

Constrain Modal to Window Size
-----------
	
You can bind the the height of the modal body to the window with something like this:
	
    $.fn.modal.defaults.maxHeight = function(){
        // subtract the height of the modal header and footer
        return $(window).height() - 165; 
    }
	
**Note:** This will be overwritten by the responsiveness and is only set when the modal is displayed, not when the window is resized.
	
Tab Index for Modal Forms
-----------
You can use `data-tabindex` instead of the default `tabindex` to specify the tabindex within a modal.

    <input type="text" data-tabindex="1" />
    <input type="text" data-tabindex="2" />

See the stackable example on the [demo](http://jschr.github.com/bootstrap-modal/) page for an example.


	



# Plupload

Plupload is a cross-browser multi-runtime file uploading API. Basically, a set of tools that will help you to 
build a reliable and visually appealing file uploader in minutes.

Historically, Plupload comes from a dark and hostile age of no HTML5, hence all the alternative fallbacks, 
like Flash, Silverlight and Java (still in development). It is meant to provide an API, that 
will work anywhere and in any case, in one way or another. While having very solid fallbacks, Plupload 
is built with the future of HTML5 in mind.

### Table of Contents
* [Backstory](https://github.com/moxiecode/plupload/blob/master/readme.md#backstory)
* [Structure](https://github.com/moxiecode/plupload/blob/master/readme.md#structure)
  * [File API and XHR L2 pollyfills](https://github.com/moxiecode/moxie/blob/master/README.md)
  * [Plupload API](https://github.com/moxiecode/plupload/wiki/API)
  * [UI Widget](https://github.com/moxiecode/plupload/wiki/UI.Plupload)
  * [Queue Widget](https://github.com/moxiecode/plupload/wiki/pluploadQueue)
* [Demos](https://github.com/jayarjo/plupload-demos/blob/master/README.md)
* [Building Instructions](https://github.com/moxiecode/plupload/blob/master/readme.md#build)
* [Getting Started](https://github.com/moxiecode/plupload/wiki/Getting-Started)
  * [Options](https://github.com/moxiecode/plupload/wiki/Options)
  * [Events](https://github.com/moxiecode/plupload/wiki/Uploader#wiki-events)
  * [Methods](https://github.com/moxiecode/plupload/wiki/Uploader#wiki-methods)
  * [Plupload in Your Language](https://github.com/moxiecode/plupload/wiki/Plupload-in-Your-Language)
  * [File Filters](https://github.com/moxiecode/plupload/wiki/File-Filters) 
  * [Image Resizing on Client-Side](https://github.com/moxiecode/plupload/wiki/Image-Resizing-on-Client-Side) 
  * [Chunking](https://github.com/moxiecode/plupload/wiki/Chunking) 
  * [Upload to Amazon S3](https://github.com/moxiecode/plupload/wiki/Upload-to-Amazon-S3) 
* [FAQ](https://github.com/moxiecode/plupload/wiki/Frequently-Asked-Questions)
* [Support](https://github.com/moxiecode/plupload/blob/master/readme.md##support)
  * [Create a Fiddle](https://github.com/moxiecode/plupload/wiki/Create-a-Fiddle)
* [Contributing](https://github.com/moxiecode/plupload/blob/master/readme.md#contribute)
* [License](https://github.com/moxiecode/plupload/blob/master/readme.md#license)
* [Contact Us](http://www.moxiecode.com/contact.php)

<a name="backstory" />
### Backstory

Plupload started in a time when uploading a file in a responsive and customizable manner was a real pain. 
Internally, browsers only had the `input[type="file"]` element. It was ugly and clunky at the same time. 
One couldn't even change it's visuals, without hiding it and coding another one on top of it from scratch. 
And then there was no progress indication for the upload process... Sounds pretty crazy today.

It was very logical for developers to look for alternatives and writing their own implementations, using 
Flash and Java, in order to somehow extend limited browser capabilities. And so did we, in our search for 
a reliable and flexible file uploader for 
our [TinyMCE](http://www.tinymce.com/index.php)'s
[MCImageManager](http://www.tinymce.com/enterprise/mcimagemanager.php). 

Quickly enough though, Plupload grew big.  It easily split into a standalone project. 
With major *version 2.0* it underwent another huge reconstruction, basically 
[from the ground up](http://blog.moxiecode.com/2012/11/28/first-public-beta-plupload-2/), 
as all the low-level runtime logic has been extracted into separate [File API](http://www.w3.org/TR/FileAPI/) 
and [XHR L2](http://www.w3.org/TR/XMLHttpRequest/) pollyfills (currently known under combined name of [mOxie](https://github.com/moxiecode/moxie)), 
giving Plupload a chance to evolve further.

<a name="structure" />
### Structure

Currently, Plupload may be considered as consisting of three parts: low-level pollyfills, 
Plupload API and Widgets (UI and Queue). Initially, Widgets were meant only to serve as examples 
of the API, but quickly formed into fully-functional API implementations that now come bundled with 
the Plupload API. This has been a source for multiple misconceptions about the API as Widgets were 
easily mistaken for the Plupload itself. They are only implementations, such as any of you can 
build by yourself out of the API.

* [Low-level pollyfills (mOxie)](https://github.com/moxiecode/moxie) - have their own [code base](https://github.com/moxiecode/moxie) and [documentation](https://github.com/moxiecode/moxie/wiki) on GitHub.
* [Plupload API](https://github.com/moxiecode/plupload/wiki/API)
* [UI Widget](https://github.com/moxiecode/plupload/wiki/UI.Plupload)
* [Queue Widget](https://github.com/moxiecode/plupload/wiki/pluploadQueue)

<a name="build" />
### Building instructions

Plupload depends on File API and XHR2 L2 pollyfills that currently have their 
[own repository](https://github.com/moxiecode/moxie) on GitHub. However, in most cases you shouldn't 
care as we bundle the latest build of mOxie, including full and minified JavaScript source and 
pre-compiled `SWF` and `XAP` components, with [every release](https://github.com/moxiecode/plupload/releases). You can find everything you may need under `js/` folder.

There are cases where you might need a custom build, for example free of unnecessary runtimes, half the 
original size, etc. The difficult part of this task comes from mOxie and its set of additional runtimes 
that require special tools on your workstation in order to compile. 
Consider [build instructions for mOxie](https://github.com/moxiecode/moxie#build-instructions) - 
everything applies to Plupload as well.

First of all, if you want to build custom Plupload packages you will require [Node.js](http://nodejs.org/), 
as this is our build environment of choice. Node.js binaries (as well as Source)
[are available](http://nodejs.org/download/) for all major operating systems.

Plupload includes _mOxie_ as a submodule, it also depends on some other repositories for building up it's dev
environment - to avoid necessity of downloading them one by one, we recommended you to simply clone Plupload 
with [git](http://git-scm.com/) recursively (you will require git installed on your system for this operation 
to succeed):

```
git clone --recursive https://github.com/moxiecode/plupload.git
```

And finalize the preparation stage with: `npm install` - this will install all additional modules, including those
required by dev and test environments. In case you would rather keep it minimal, add a `--production` flag.

*Note:* Currently, for an unknown reason, locally installed Node.js modules on Windows, may not be automatically 
added to the system PATH. So, if `jake` commands below are not recognized you will need to add them manually:

```
set PATH=%PATH%;%CD%\node_modules\.bin\
``` 

<a name="support" />
### Support

We are actively standing behind the Plupload and now that we are done with major rewrites and refactoring,
the only real goal that we have ahead is making it as reliable and bulletproof as possible. We are open to 
all the suggestions and feature requests. We ask you to file bug reports if you encounter any. We may not 
react to them instantly, but we constantly bear them in my mind as we extend the code base.

In addition to dedicated support for those who dare to buy our OEM licenses, we got 
[discussion boards](http://www.plupload.com/punbb/index.php), which is like an enormous FAQ, 
covering every possible application case. Of course, you are welcome to file a bug report or feature request, 
here on [GitHub](https://github.com/moxiecode/plupload/issues).

Sometimes it is easier to notice the problem when bug report is accompained by the actual code. Consider providing 
[a Plupload fiddle](https://github.com/moxiecode/plupload/wiki/Create-a-Fiddle) for the troublesome code.

<a name="contribute" />
### Contributing

We are open to suggestions and code revisions, however there are some rules and limitations that you might 
want to consider first.

* Code that you contribute will automatically be licensed under the LGPL, but will not be limited to LGPL.
* Although all contributors will get the credit for their work, copyright notices will be changed to [Moxiecode Systems AB](http://www.moxiecode.com/).
* Third party code will be reviewed, tested and possibly modified before being released.

These basic rules help us earn a living and ensure that code remains Open Source and compatible with LGPL license. All contributions will be added to the changelog and appear in every release and on the site. 

An easy place to start is to [translate Plupload to your language](https://github.com/moxiecode/plupload/wiki/Plupload-in-Your-Language#contribute).

You can read more about how to contribute at: [http://www.plupload.com/contributing](http://www.plupload.com/contributing)

<a name="license" />
### License

Copyright 2013, [Moxiecode Systems AB](http://www.moxiecode.com/)  
Released under [GPLv2 License](https://github.com/moxiecode/plupload/blob/master/license.txt).

We also provide [commercial license](http://www.plupload.com/commercial.php).
[Twitter Bootstrap](http://twitter.github.com/bootstrap) [![Build Status](https://secure.travis-ci.org/twitter/bootstrap.png)](http://travis-ci.org/twitter/bootstrap)
=================

Bootstrap is a sleek, intuitive, and powerful front-end framework for faster and easier web development, created and maintained by [Mark Otto](http://twitter.com/mdo) and [Jacob Thornton](http://twitter.com/fat) at Twitter.

To get started, checkout http://getbootstrap.com!



Quick start
-----------

Clone the repo, `git clone git://github.com/twitter/bootstrap.git`, or [download the latest release](https://github.com/twitter/bootstrap/zipball/master).



Versioning
----------

For transparency and insight into our release cycle, and for striving to maintain backward compatibility, Bootstrap will be maintained under the Semantic Versioning guidelines as much as possible.

Releases will be numbered with the following format:

`<major>.<minor>.<patch>`

And constructed with the following guidelines:

* Breaking backward compatibility bumps the major (and resets the minor and patch)
* New additions without breaking backward compatibility bumps the minor (and resets the patch)
* Bug fixes and misc changes bumps the patch

For more information on SemVer, please visit http://semver.org/.



Bug tracker
-----------

Have a bug? Please create an issue here on GitHub that conforms with [necolas's guidelines](https://github.com/necolas/issue-guidelines).

https://github.com/twitter/bootstrap/issues



Twitter account
---------------

Keep up to date on announcements and more by following Bootstrap on Twitter, [@TwBootstrap](http://twitter.com/TwBootstrap).



Blog
----

Read more detailed announcements, discussions, and more on [The Official Twitter Bootstrap Blog](http://blog.getbootstrap.com).



Mailing list
------------

Have a question? Ask on our mailing list!

twitter-bootstrap@googlegroups.com

http://groups.google.com/group/twitter-bootstrap



IRC
---

Server: irc.freenode.net

Channel: ##twitter-bootstrap (the double ## is not a typo)



Developers
----------

We have included a makefile with convenience methods for working with the Bootstrap library.

+ **dependencies**
Our makefile depends on you having recess, connect, uglify.js, and jshint installed. To install, just run the following command in npm:

```
$ npm install recess connect uglify-js jshint -g
```

+ **build** - `make`
Runs the recess compiler to rebuild the `/less` files and compiles the docs pages. Requires recess and uglify-js. <a href="http://twitter.github.com/bootstrap/less.html#compiling">Read more in our docs &raquo;</a>

+ **test** - `make test`
Runs jshint and qunit tests headlessly in [phantomjs](http://code.google.com/p/phantomjs/) (used for ci). Depends on having phantomjs installed.

+ **watch** - `make watch`
This is a convenience method for watching just Less files and automatically building them whenever you save. Requires the Watchr gem.



Contributing
------------

Please submit all pull requests against *-wip branches. If your unit test contains javascript patches or features, you must include relevant unit tests. Thanks!



Authors
-------

**Mark Otto**

+ http://twitter.com/mdo
+ http://github.com/markdotto

**Jacob Thornton**

+ http://twitter.com/fat
+ http://github.com/fat



Copyright and license
---------------------

Copyright 2012 Twitter, Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this work except in compliance with the License.
You may obtain a copy of the License in the LICENSE file, or at:

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
## 2.0 BOOTSTRAP JS PHILOSOPHY
These are the high-level design rules which guide the development of Bootstrap's plugin apis.

---

### DATA-ATTRIBUTE API

We believe you should be able to use all plugins provided by Bootstrap purely through the markup API without writing a single line of javascript.

We acknowledge that this isn't always the most performant and sometimes it may be desirable to turn this functionality off altogether. Therefore, as of 2.0 we provide the ability to disable the data attribute API by unbinding all events on the body namespaced with `'data-api'`. This looks like this:

    $('body').off('.data-api')

To target a specific plugin, just include the plugins name as a namespace along with the data-api namespace like this:

    $('body').off('.alert.data-api')

---

### PROGRAMATIC API

We also believe you should be able to use all plugins provided by Bootstrap purely through the JS API.

All public APIs should be single, chainable methods, and return the collection acted upon.

    $(".btn.danger").button("toggle").addClass("fat")

All methods should accept an optional options object, a string which targets a particular method, or null which initiates the default behavior:

    $("#myModal").modal() // initialized with defaults
    $("#myModal").modal({ keyboard: false }) // initialized with now keyboard
    $("#myModal").modal('show') // initializes and invokes show immediately afterqwe2

---

### OPTIONS

Options should be sparse and add universal value. We should pick the right defaults.

All plugins should have a default object which can be modified to effect all instance's default options. The defaults object should be available via `$.fn.plugin.defaults`.

    $.fn.modal.defaults = { … }

An options definition should take the following form:

    *noun*: *adjective* - describes or modifies a quality of an instance

examples:

    backdrop: true
    keyboard: false
    placement: 'top'

---

### EVENTS

All events should have an infinitive and past participle form. The infinitive is fired just before an action takes place, the past participle on completion of the action.

    show | shown
    hide | hidden

---

### CONSTRUCTORS

Each plugin should expose it's raw constructor on a `Constructor` property -- accessed in the following way:


    $.fn.popover.Constructor

---

### DATA ACCESSOR

Each plugin stores a copy of the invoked class on an object. This class instance can be accessed directly through jQuery's data API like this:

    $('[rel=popover]').data('popover') instanceof $.fn.popover.Constructor

---

### DATA ATTRIBUTES

Data attributes should take the following form:

- data-{{verb}}={{plugin}} - defines main interaction
- data-target || href^=# - defined on "control" element (if element controls an element other than self)
- data-{{noun}} - defines class instance options

examples:

    // control other targets
    data-toggle="modal" data-target="#foo"
    data-toggle="collapse" data-target="#foo" data-parent="#bar"

    // defined on element they control
    data-spy="scroll"

    data-dismiss="modal"
    data-dismiss="alert"

    data-toggle="dropdown"

    data-toggle="button"
    data-toggle="buttons-checkbox"
    data-toggle="buttons-radio"